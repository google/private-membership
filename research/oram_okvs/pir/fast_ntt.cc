// Copyright 2026 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "research/oram_okvs/pir/fast_ntt.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "shell_encryption/dft_transformations.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE \
  "research/oram_okvs/pir/fast_ntt.cc"
#include "hwy/foreach_target.h"  // IWYU pragma: keep
#include "hwy/highway.h"

HWY_BEFORE_NAMESPACE();
namespace oram_okvs::pir {
namespace HWY_NAMESPACE {

#if HWY_TARGET == HWY_SCALAR

template <typename ModularInt>
void FastIterativeCooleyTukeyHwy(
    std::vector<ModularInt>& coeffs, int log_len,
    const rlwe::NttParameters<ModularInt>& ntt_params,
    const typename ModularInt::Params& mod_params) {
  rlwe::IterativeCooleyTukey<ModularInt>(coeffs, log_len, ntt_params,
                                         mod_params);
}

template <typename ModularInt>
void FastIterativeGentlemanSandeHwy(
    std::vector<ModularInt>& coeffs, int log_len,
    const rlwe::NttParameters<ModularInt>& ntt_params,
    const typename ModularInt::Params& mod_params) {
  rlwe::IterativeGentlemanSande<ModularInt>(coeffs, log_len, ntt_params,
                                            mod_params)
      .IgnoreError();
  ModularInt n_inv = ntt_params.n_inv_ptr.value();
  for (auto& coeff : coeffs) {
    coeff.MulInPlace(n_inv, &mod_params);
  }
}

#else

namespace hn = hwy::HWY_NAMESPACE;

// Reduces x in range [0, 2 * q) to [0, q).
template <class D, class V>
HWY_INLINE V ReduceOnceHwy(D d, V x, V q) {
  auto mask = hn::Ge(x, q);
  return hn::Sub(x, hn::IfThenElseZero(mask, q));
}

// Multiplies vector v by Barrett constant (psi, psi_barrett), returning t in
// [0, 2q).
template <class D, class V>
HWY_INLINE V MulBarrettConstantHwy(D d, V v, V psi, V psi_barrett,
                                   V neg_modulus) {
  V q = hn::MulHigh(psi_barrett, v);
  V psi_v = hn::Mul(psi, v);
  return hn::Add(psi_v, hn::Mul(q, neg_modulus));
}

// Fast Cooley-Tukey Forward NTT with 4x loop unrolling and Harvey reduction
template <typename ModularInt>
void FastIterativeCooleyTukeyHwy(
    std::vector<ModularInt>& coeffs, int log_len,
    const rlwe::NttParameters<ModularInt>& ntt_params,
    const typename ModularInt::Params& mod_params) {
  using Int = typename ModularInt::Int;
  const hn::ScalableTag<Int> d;
  const int num_lanes = hn::Lanes(d);
  const int log_num_lanes = static_cast<int>(std::ceil(std::log2(num_lanes)));

  using Vec = hn::Vec<decltype(d)>;

  const Vec modulus_vec = hn::Set(d, mod_params.modulus);
  const Vec neg_modulus = hn::Set(d, -mod_params.modulus);
  const Vec twice_modulus = hn::Set(d, 2 * mod_params.modulus);

  int index_psi = 1;
  int n = 1 << log_len;
  Int* base_ptr = reinterpret_cast<Int*>(coeffs.data());

  // Vectorized layers (layers log_len - 1 down to log_num_lanes)
  int i = log_len - 1;
  for (; i >= log_num_lanes; --i) {
    const int half_m = 1 << i;
    const int m = half_m << 1;

    for (int k = 0; k < n; k += m) {
      const auto& [psi_scalar, psi_barrett_scalar] =
          ntt_params.psis_bitrev_constant[index_psi++];
      Vec psi = hn::Set(d, psi_scalar);
      Vec psi_barrett = hn::Set(d, psi_barrett_scalar);

      // 4x unrolled vector loop when half_m >= 4 * num_lanes
      int j = 0;
      for (; j + 4 * num_lanes <= half_m; j += 4 * num_lanes) {
        Int* p_u0 = base_ptr + k + j;
        Int* p_u1 = p_u0 + num_lanes;
        Int* p_u2 = p_u1 + num_lanes;
        Int* p_u3 = p_u2 + num_lanes;

        Int* p_v0 = p_u0 + half_m;
        Int* p_v1 = p_u1 + half_m;
        Int* p_v2 = p_u2 + half_m;
        Int* p_v3 = p_u3 + half_m;

        Vec v0 = hn::LoadU(d, p_v0);
        Vec v1 = hn::LoadU(d, p_v1);
        Vec v2 = hn::LoadU(d, p_v2);
        Vec v3 = hn::LoadU(d, p_v3);

        Vec t0 = MulBarrettConstantHwy(d, v0, psi, psi_barrett, neg_modulus);
        Vec t1 = MulBarrettConstantHwy(d, v1, psi, psi_barrett, neg_modulus);
        Vec t2 = MulBarrettConstantHwy(d, v2, psi, psi_barrett, neg_modulus);
        Vec t3 = MulBarrettConstantHwy(d, v3, psi, psi_barrett, neg_modulus);

        Vec u0 = hn::LoadU(d, p_u0);
        Vec u1 = hn::LoadU(d, p_u1);
        Vec u2 = hn::LoadU(d, p_u2);
        Vec u3 = hn::LoadU(d, p_u3);

        u0 = ReduceOnceHwy(d, u0, twice_modulus);
        u1 = ReduceOnceHwy(d, u1, twice_modulus);
        u2 = ReduceOnceHwy(d, u2, twice_modulus);
        u3 = ReduceOnceHwy(d, u3, twice_modulus);

        Vec new_v0 = hn::Add(u0, hn::Sub(twice_modulus, t0));
        Vec new_u0 = hn::Add(u0, t0);
        Vec new_v1 = hn::Add(u1, hn::Sub(twice_modulus, t1));
        Vec new_u1 = hn::Add(u1, t1);
        Vec new_v2 = hn::Add(u2, hn::Sub(twice_modulus, t2));
        Vec new_u2 = hn::Add(u2, t2);
        Vec new_v3 = hn::Add(u3, hn::Sub(twice_modulus, t3));
        Vec new_u3 = hn::Add(u3, t3);

        hn::StoreU(new_u0, d, p_u0);
        hn::StoreU(new_u1, d, p_u1);
        hn::StoreU(new_u2, d, p_u2);
        hn::StoreU(new_u3, d, p_u3);

        hn::StoreU(new_v0, d, p_v0);
        hn::StoreU(new_v1, d, p_v1);
        hn::StoreU(new_v2, d, p_v2);
        hn::StoreU(new_v3, d, p_v3);
      }

      // Remainder vector loop
      for (; j < half_m; j += num_lanes) {
        Int* p_u = base_ptr + k + j;
        Int* p_v = p_u + half_m;

        Vec v = hn::LoadU(d, p_v);
        Vec t = MulBarrettConstantHwy(d, v, psi, psi_barrett, neg_modulus);

        Vec u = hn::LoadU(d, p_u);
        u = ReduceOnceHwy(d, u, twice_modulus);

        Vec new_v = hn::Add(u, hn::Sub(twice_modulus, t));
        Vec new_u = hn::Add(u, t);

        hn::StoreU(new_u, d, p_u);
        hn::StoreU(new_v, d, p_v);
      }
    }
  }

  // Remaining inner layers (i < log_num_lanes)
  bool should_lazy_add = false;
  for (; i >= 0; --i) {
    const int half_m = 1 << i;
    const int m = half_m << 1;

    for (int k = 0; k < n; k += m) {
      const auto& [psi_scalar, psi_barrett_scalar] =
          ntt_params.psis_bitrev_constant[index_psi++];
      for (int j = 0; j < half_m; ++j) {
        ModularInt t = coeffs[k + j + half_m].MulConstant(
            psi_scalar, psi_barrett_scalar, &mod_params);
        ModularInt u = coeffs[k + j];
        if (i != 0 && should_lazy_add) {
          coeffs[k + j].LazyAddInPlace(t, &mod_params);
          coeffs[k + j + half_m] = std::move(u.LazySubInPlace(t, &mod_params));
        } else {
          coeffs[k + j].AddInPlace(t, &mod_params);
          coeffs[k + j + half_m] = std::move(u.SubInPlace(t, &mod_params));
        }
      }
    }
    should_lazy_add = !should_lazy_add;
  }

  // Final reduction pass: ensure all coefficients are in [0, modulus)
  for (int j = 0; j < n; j += num_lanes) {
    Int* p = base_ptr + j;
    Vec x = hn::LoadU(d, p);
    x = ReduceOnceHwy(d, x, twice_modulus);
    x = ReduceOnceHwy(d, x, modulus_vec);
    hn::StoreU(x, d, p);
  }
}

// Fast Gentleman-Sande Inverse NTT
template <typename ModularInt>
void FastIterativeGentlemanSandeHwy(
    std::vector<ModularInt>& coeffs, int log_len,
    const rlwe::NttParameters<ModularInt>& ntt_params,
    const typename ModularInt::Params& mod_params) {
  using Int = typename ModularInt::Int;
  const hn::ScalableTag<Int> d;
  const int num_lanes = hn::Lanes(d);
  const int log_num_lanes = static_cast<int>(std::ceil(std::log2(num_lanes)));

  using Vec = hn::Vec<decltype(d)>;

  const Vec modulus_vec = hn::Set(d, mod_params.modulus);
  const Vec neg_modulus = hn::Set(d, -mod_params.modulus);
  const Vec twice_modulus = hn::Set(d, 2 * mod_params.modulus);

  int index_psi_inv = 0;
  int n = 1 << log_len;
  Int* base_ptr = reinterpret_cast<Int*>(coeffs.data());

  // Lower layers (i = 0 to log_num_lanes - 1)
  int i = 0;
  for (; i < log_num_lanes && i < log_len; ++i) {
    const int half_m = 1 << i;
    const int m = half_m << 1;

    for (int k = 0; k < n; k += m) {
      const auto& [psi_inv_scalar, psi_inv_barrett_scalar] =
          ntt_params.psis_inv_bitrev_constant[index_psi_inv++];
      for (int j = 0; j < half_m; ++j) {
        ModularInt t = coeffs[k + j + half_m];
        ModularInt u = coeffs[k + j];
        coeffs[k + j].AddInPlace(t, &mod_params);
        coeffs[k + j + half_m] = std::move(
            u.LazySubInPlace(t, &mod_params)
                .MulConstantInPlace(psi_inv_scalar, psi_inv_barrett_scalar,
                                    &mod_params));
      }
    }
  }

  // Vectorized layers (i = log_num_lanes to log_len - 1)
  for (; i < log_len; ++i) {
    const int half_m = 1 << i;
    const int m = half_m << 1;

    for (int k = 0; k < n; k += m) {
      const auto& [psi_inv_scalar, psi_inv_barrett_scalar] =
          ntt_params.psis_inv_bitrev_constant[index_psi_inv++];
      Vec psi_inv = hn::Set(d, psi_inv_scalar);
      Vec psi_inv_barrett = hn::Set(d, psi_inv_barrett_scalar);

      // 4x unrolled vector loop
      int j = 0;
      for (; j + 4 * num_lanes <= half_m; j += 4 * num_lanes) {
        Int* p_u0 = base_ptr + k + j;
        Int* p_u1 = p_u0 + num_lanes;
        Int* p_u2 = p_u1 + num_lanes;
        Int* p_u3 = p_u2 + num_lanes;

        Int* p_v0 = p_u0 + half_m;
        Int* p_v1 = p_u1 + half_m;
        Int* p_v2 = p_u2 + half_m;
        Int* p_v3 = p_u3 + half_m;

        Vec u0 = hn::LoadU(d, p_u0);
        Vec u1 = hn::LoadU(d, p_u1);
        Vec u2 = hn::LoadU(d, p_u2);
        Vec u3 = hn::LoadU(d, p_u3);

        Vec v0 = hn::LoadU(d, p_v0);
        Vec v1 = hn::LoadU(d, p_v1);
        Vec v2 = hn::LoadU(d, p_v2);
        Vec v3 = hn::LoadU(d, p_v3);

        Vec sum0 = ReduceOnceHwy(d, hn::Add(u0, v0), twice_modulus);
        Vec sum1 = ReduceOnceHwy(d, hn::Add(u1, v1), twice_modulus);
        Vec sum2 = ReduceOnceHwy(d, hn::Add(u2, v2), twice_modulus);
        Vec sum3 = ReduceOnceHwy(d, hn::Add(u3, v3), twice_modulus);

        Vec diff0 = hn::Sub(hn::Add(u0, twice_modulus), v0);
        Vec diff1 = hn::Sub(hn::Add(u1, twice_modulus), v1);
        Vec diff2 = hn::Sub(hn::Add(u2, twice_modulus), v2);
        Vec diff3 = hn::Sub(hn::Add(u3, twice_modulus), v3);

        Vec prod0 = MulBarrettConstantHwy(d, diff0, psi_inv, psi_inv_barrett,
                                          neg_modulus);
        Vec prod1 = MulBarrettConstantHwy(d, diff1, psi_inv, psi_inv_barrett,
                                          neg_modulus);
        Vec prod2 = MulBarrettConstantHwy(d, diff2, psi_inv, psi_inv_barrett,
                                          neg_modulus);
        Vec prod3 = MulBarrettConstantHwy(d, diff3, psi_inv, psi_inv_barrett,
                                          neg_modulus);

        hn::StoreU(sum0, d, p_u0);
        hn::StoreU(sum1, d, p_u1);
        hn::StoreU(sum2, d, p_u2);
        hn::StoreU(sum3, d, p_u3);

        hn::StoreU(prod0, d, p_v0);
        hn::StoreU(prod1, d, p_v1);
        hn::StoreU(prod2, d, p_v2);
        hn::StoreU(prod3, d, p_v3);
      }

      for (; j < half_m; j += num_lanes) {
        Int* p_u = base_ptr + k + j;
        Int* p_v = p_u + half_m;

        Vec u = hn::LoadU(d, p_u);
        Vec v = hn::LoadU(d, p_v);

        Vec sum = ReduceOnceHwy(d, hn::Add(u, v), twice_modulus);
        Vec diff = hn::Sub(hn::Add(u, twice_modulus), v);
        Vec prod = MulBarrettConstantHwy(d, diff, psi_inv, psi_inv_barrett,
                                         neg_modulus);

        hn::StoreU(sum, d, p_u);
        hn::StoreU(prod, d, p_v);
      }
    }
  }

  // Multiply by n_inv and reduce to [0, modulus)
  ModularInt n_inv = ntt_params.n_inv_ptr.value();
  auto [n_inv_scalar, n_inv_barrett_scalar] = n_inv.GetConstant(&mod_params);
  Vec n_inv_vec = hn::Set(d, n_inv_scalar);
  Vec n_inv_barrett = hn::Set(d, n_inv_barrett_scalar);

  for (int j = 0; j < n; j += num_lanes) {
    Int* p = base_ptr + j;
    Vec x = hn::LoadU(d, p);
    x = MulBarrettConstantHwy(d, x, n_inv_vec, n_inv_barrett, neg_modulus);
    x = ReduceOnceHwy(d, x, modulus_vec);
    hn::StoreU(x, d, p);
  }
}

#endif

}  // namespace HWY_NAMESPACE
}  // namespace oram_okvs::pir
HWY_AFTER_NAMESPACE();

#if HWY_ONCE || HWY_IDE

namespace oram_okvs::pir {

template <typename ModularInt>
absl::Status FastForwardNtt(std::vector<ModularInt>& coeffs,
                            const rlwe::NttParameters<ModularInt>& ntt_params,
                            const typename ModularInt::Params& mod_params) {
  int log_len = static_cast<int>(std::log2(coeffs.size()));
  int expected_len = 1 << log_len;
  if (coeffs.size() != expected_len) {
    return absl::InvalidArgumentError(
        absl::StrCat("`coeffs` must contain ", expected_len, " elements."));
  }

  HWY_EXPORT_AND_DYNAMIC_DISPATCH_T(FastIterativeCooleyTukeyHwy<ModularInt>)
  (coeffs, log_len, ntt_params, mod_params);
  return absl::OkStatus();
}

template <typename ModularInt>
absl::Status FastInverseNtt(std::vector<ModularInt>& coeffs,
                            const rlwe::NttParameters<ModularInt>& ntt_params,
                            const typename ModularInt::Params& mod_params) {
  int log_len = static_cast<int>(std::log2(coeffs.size()));
  int expected_len = 1 << log_len;
  if (coeffs.size() != expected_len) {
    return absl::InvalidArgumentError(
        absl::StrCat("`coeffs` must contain ", expected_len, " elements."));
  }

  HWY_EXPORT_AND_DYNAMIC_DISPATCH_T(FastIterativeGentlemanSandeHwy<ModularInt>)
  (coeffs, log_len, ntt_params, mod_params);
  return absl::OkStatus();
}

template <typename ModularInt>
absl::Status FastConvertToNttForm(
    rlwe::RnsPolynomial<ModularInt>& poly,
    absl::Span<const rlwe::PrimeModulus<ModularInt>* const> moduli) {
  if (poly.IsNttForm()) {
    return absl::InvalidArgumentError("Polynomial already in NTT form");
  }
  int num_moduli = poly.NumModuli();
  if (moduli.size() != num_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
  }
  std::vector<std::vector<ModularInt>> coeff_vectors = poly.Coeffs();
  for (int i = 0; i < num_moduli; ++i) {
    RLWE_RETURN_IF_ERROR(FastForwardNtt(
        coeff_vectors[i], *moduli[i]->NttParams(), *moduli[i]->ModParams()));
  }
  auto res_or = rlwe::RnsPolynomial<ModularInt>::Create(
      std::move(coeff_vectors), /*is_ntt=*/true);
  if (!res_or.ok()) return res_or.status();
  poly = std::move(*res_or);
  return absl::OkStatus();
}

template <typename ModularInt>
absl::Status FastConvertToCoeffForm(
    rlwe::RnsPolynomial<ModularInt>& poly,
    absl::Span<const rlwe::PrimeModulus<ModularInt>* const> moduli) {
  if (!poly.IsNttForm()) {
    return absl::InvalidArgumentError("Polynomial already in Coefficient form");
  }
  int num_moduli = poly.NumModuli();
  if (moduli.size() != num_moduli) {
    return absl::InvalidArgumentError(
        absl::StrCat("`moduli` must contain ", num_moduli, " RNS moduli."));
  }
  std::vector<std::vector<ModularInt>> coeff_vectors = poly.Coeffs();
  for (int i = 0; i < num_moduli; ++i) {
    RLWE_RETURN_IF_ERROR(FastInverseNtt(
        coeff_vectors[i], *moduli[i]->NttParams(), *moduli[i]->ModParams()));
  }
  auto res_or = rlwe::RnsPolynomial<ModularInt>::Create(
      std::move(coeff_vectors), /*is_ntt=*/false);
  if (!res_or.ok()) return res_or.status();
  poly = std::move(*res_or);
  return absl::OkStatus();
}

// Explicit template instantiations for supported ModularInt types
template absl::Status FastForwardNtt(
    std::vector<rlwe::MontgomeryInt<uint64_t>>& coeffs,
    const rlwe::NttParameters<rlwe::MontgomeryInt<uint64_t>>& ntt_params,
    const rlwe::MontgomeryInt<uint64_t>::Params& mod_params);

template absl::Status FastInverseNtt(
    std::vector<rlwe::MontgomeryInt<uint64_t>>& coeffs,
    const rlwe::NttParameters<rlwe::MontgomeryInt<uint64_t>>& ntt_params,
    const rlwe::MontgomeryInt<uint64_t>::Params& mod_params);

template absl::Status FastConvertToNttForm(
    rlwe::RnsPolynomial<rlwe::MontgomeryInt<uint64_t>>& poly,
    absl::Span<const rlwe::PrimeModulus<rlwe::MontgomeryInt<uint64_t>>* const>
        moduli);

template absl::Status FastConvertToCoeffForm(
    rlwe::RnsPolynomial<rlwe::MontgomeryInt<uint64_t>>& poly,
    absl::Span<const rlwe::PrimeModulus<rlwe::MontgomeryInt<uint64_t>>* const>
        moduli);

template absl::Status FastForwardNtt(
    std::vector<rlwe::MontgomeryInt<uint32_t>>& coeffs,
    const rlwe::NttParameters<rlwe::MontgomeryInt<uint32_t>>& ntt_params,
    const rlwe::MontgomeryInt<uint32_t>::Params& mod_params);

template absl::Status FastInverseNtt(
    std::vector<rlwe::MontgomeryInt<uint32_t>>& coeffs,
    const rlwe::NttParameters<rlwe::MontgomeryInt<uint32_t>>& ntt_params,
    const rlwe::MontgomeryInt<uint32_t>::Params& mod_params);

template absl::Status FastConvertToNttForm(
    rlwe::RnsPolynomial<rlwe::MontgomeryInt<uint32_t>>& poly,
    absl::Span<const rlwe::PrimeModulus<rlwe::MontgomeryInt<uint32_t>>* const>
        moduli);

template absl::Status FastConvertToCoeffForm(
    rlwe::RnsPolynomial<rlwe::MontgomeryInt<uint32_t>>& poly,
    absl::Span<const rlwe::PrimeModulus<rlwe::MontgomeryInt<uint32_t>>* const>
        moduli);

}  // namespace oram_okvs::pir

#endif
