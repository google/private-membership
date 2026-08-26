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

#include "research/oram_okvs/pir/onion_pir_client.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/memory/memory.h"
#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "research/oram_okvs/pir/fast_ntt.h"
#include "research/oram_okvs/pir/onion_pir_params.h"
#include "research/oram_okvs/pir/onion_pir_types.h"
#include "research/oram_okvs/pir/pir.pb.h"
#include "shell_encryption/constants.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/prng/hkdf_prng.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/rns/error_distribution.h"
#include "shell_encryption/rns/rns_bfv_ciphertext.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_gadget.h"
#include "shell_encryption/rns/rns_galois_key.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/rns_secret_key.h"
#include "shell_encryption/status_macros.h"

namespace oram_okvs::pir {

OnionPirClient::OnionPirClient(
    OnionPirParams params, std::unique_ptr<const RnsContext> rns_context,
    std::unique_ptr<const RnsErrorParams> error_params,
    std::unique_ptr<const rlwe::RnsGadget<ModularInt>> conversion_gadget,
    std::unique_ptr<const rlwe::RnsGadget<ModularInt>> expansion_gadget,
    std::unique_ptr<const rlwe::RnsGadget<ModularInt>> trace_gadget,
    std::unique_ptr<const rlwe::CoefficientEncoder<ModularInt>> encoder,
    rlwe::RnsRlweSecretKey<ModularInt> secret_key,
    std::vector<rlwe::RnsGaloisKey<ModularInt>> galois_keys,
    RgswCiphertext rgsw_secret_key, std::vector<int64_t> sk_coeffs)
    : params_(std::move(params)),
      rns_context_(std::move(rns_context)),
      error_params_(std::move(error_params)),
      conversion_gadget_(std::move(conversion_gadget)),
      expansion_gadget_(std::move(expansion_gadget)),
      trace_gadget_(std::move(trace_gadget)),
      encoder_(std::move(encoder)),
      secret_key_(std::move(secret_key)),
      galois_keys_(std::move(galois_keys)),
      rgsw_secret_key_(std::move(rgsw_secret_key)),
      sk_coeffs_(std::move(sk_coeffs)) {}

absl::StatusOr<std::unique_ptr<OnionPirClient>> OnionPirClient::Create(
    OnionPirParams params) {
  if (params.poly_degree == 0 ||
      (params.poly_degree & (params.poly_degree - 1)) != 0) {
    return absl::InvalidArgumentError(
        "Polynomial degree must be a power of 2.");
  }
  size_t log_n = 0;
  while ((static_cast<size_t>(1) << log_n) < params.poly_degree) {
    ++log_n;
  }

  RLWE_ASSIGN_OR_RETURN(
      auto rns_context_val,
      RnsContext::CreateForBfv(log_n, params.ciphertext_moduli, /*ps=*/{},
                               params.plaintext_modulus()));
  auto rns_context =
      std::make_unique<const RnsContext>(std::move(rns_context_val));
  auto main_moduli = rns_context->MainPrimeModuli();

  RLWE_ASSIGN_OR_RETURN(
      auto error_params_val,
      RnsErrorParams::Create(log_n, main_moduli, /*aux_moduli=*/{},
                             params.log_t, sqrt(params.variance)));
  auto error_params =
      std::make_unique<const RnsErrorParams>(std::move(error_params_val));

  int level = main_moduli.size() - 1;
  RLWE_ASSIGN_OR_RETURN(auto q_hats,
                        rns_context->MainPrimeModulusComplements(level));
  RLWE_ASSIGN_OR_RETURN(auto q_hat_invs,
                        rns_context->MainPrimeModulusCrtFactors(level));

  std::vector<size_t> conv_log_bs = params.conversion_gadget_log_bs;
  if (conv_log_bs.empty()) {
    conv_log_bs.assign(main_moduli.size(), params.conversion_log_base());
  }
  RLWE_ASSIGN_OR_RETURN(
      auto conv_gadget_val,
      rlwe::RnsGadget<ModularInt>::Create(log_n, conv_log_bs, q_hats,
                                          q_hat_invs, main_moduli));
  auto conversion_gadget = std::make_unique<const rlwe::RnsGadget<ModularInt>>(
      std::move(conv_gadget_val));

  std::vector<size_t> exp_log_bs = params.expansion_gadget_log_bs;
  if (exp_log_bs.empty()) {
    exp_log_bs.assign(main_moduli.size(), params.expansion_log_base());
  }
  RLWE_ASSIGN_OR_RETURN(
      auto exp_gadget_val,
      rlwe::RnsGadget<ModularInt>::Create(log_n, exp_log_bs, q_hats, q_hat_invs,
                                          main_moduli));
  auto expansion_gadget = std::make_unique<const rlwe::RnsGadget<ModularInt>>(
      std::move(exp_gadget_val));

  std::vector<size_t> trace_log_bs = params.trace_gadget_log_bs;
  if (trace_log_bs.empty()) {
    trace_log_bs.assign(main_moduli.size(), params.trace_log_base());
  }
  RLWE_ASSIGN_OR_RETURN(
      auto trace_gadget_val,
      rlwe::RnsGadget<ModularInt>::Create(log_n, trace_log_bs, q_hats,
                                          q_hat_invs, main_moduli));
  auto trace_gadget = std::make_unique<const rlwe::RnsGadget<ModularInt>>(
      std::move(trace_gadget_val));

  RLWE_ASSIGN_OR_RETURN(
      auto encoder_val,
      rlwe::CoefficientEncoder<ModularInt>::Create(rns_context.get()));
  auto encoder = std::make_unique<const rlwe::CoefficientEncoder<ModularInt>>(
      std::move(encoder_val));

  RLWE_ASSIGN_OR_RETURN(std::string prng_seed, rlwe::HkdfPrng::GenerateSeed());
  RLWE_ASSIGN_OR_RETURN(auto prng, rlwe::HkdfPrng::Create(prng_seed));

  RLWE_ASSIGN_OR_RETURN(auto secret_key,
                        rlwe::RnsRlweSecretKey<ModularInt>::Sample(
                            log_n, params.variance, main_moduli, prng.get()));

  // Generate Galois keys for all log2(N) levels (Trace + Expansion).
  // First trace_steps keys use trace_gadget; remaining keys use
  // expansion_gadget.
  size_t total_levels = log_n;
  size_t trace_steps = params.query_trace_steps();
  std::vector<rlwe::RnsGaloisKey<ModularInt>> galois_keys;
  galois_keys.reserve(total_levels);

  for (size_t h = 0; h < total_levels; ++h) {
    size_t k = static_cast<size_t>(1) << h;
    rlwe::Uint64 sub_power = (params.poly_degree / k) + 1;
    const auto* g_ptr =
        (h < trace_steps) ? trace_gadget.get() : expansion_gadget.get();
    RLWE_ASSIGN_OR_RETURN(
        auto gkey, rlwe::RnsGaloisKey<ModularInt>::CreateForBfv(
                       secret_key, static_cast<int>(sub_power), params.variance,
                       g_ptr, rlwe::PRNG_TYPE_HKDF));
    galois_keys.push_back(std::move(gkey));
  }

  // Generate RGSW encryption of secret key s: A = RGSW(s).
  // A has 2 * ell rows of BFV ciphertexts.
  size_t ell = conversion_gadget->Dimension();
  std::vector<RnsBfvCiphertext> rgsw_rows;
  rgsw_rows.reserve(2 * ell);

  // Top ell rows: BFV(0) + (s * g_k, 0)
  for (size_t k = 0; k < ell; ++k) {
    RnsPolynomial s_gk = secret_key.Key();
    RLWE_RETURN_IF_ERROR(
        s_gk.MulInPlace(conversion_gadget->Component(k), main_moduli));

    RLWE_ASSIGN_OR_RETURN(RnsPolynomial a, RnsPolynomial::SampleUniform(
                                               log_n, prng.get(), main_moduli));

    RLWE_ASSIGN_OR_RETURN(RnsPolynomial e,
                          rlwe::SampleError<ModularInt>(
                              log_n, params.variance, main_moduli, prng.get()));

    // c0 = -a * s + e + s * g_k
    RLWE_ASSIGN_OR_RETURN(RnsPolynomial prod_as,
                          a.Mul(secret_key.Key(), main_moduli));
    RLWE_ASSIGN_OR_RETURN(RnsPolynomial neg_prod, prod_as.Negate(main_moduli));
    RLWE_RETURN_IF_ERROR(neg_prod.AddInPlace(e, main_moduli));
    RLWE_RETURN_IF_ERROR(neg_prod.AddInPlace(s_gk, main_moduli));

    std::vector<RnsPolynomial> comps;
    comps.push_back(std::move(neg_prod));
    comps.push_back(std::move(a));

    rgsw_rows.emplace_back(std::move(comps), main_moduli, /*power_of_s=*/1,
                           /*error=*/0.0, error_params.get(),
                           rns_context.get());
  }

  // Bottom ell rows: BFV(0) + (0, s * g_k)
  for (size_t k = 0; k < ell; ++k) {
    RnsPolynomial s_gk = secret_key.Key();
    RLWE_RETURN_IF_ERROR(
        s_gk.MulInPlace(conversion_gadget->Component(k), main_moduli));

    RLWE_ASSIGN_OR_RETURN(RnsPolynomial a, RnsPolynomial::SampleUniform(
                                               log_n, prng.get(), main_moduli));

    RLWE_ASSIGN_OR_RETURN(RnsPolynomial e,
                          rlwe::SampleError<ModularInt>(
                              log_n, params.variance, main_moduli, prng.get()));

    // c0 = -a * s + e
    RLWE_ASSIGN_OR_RETURN(RnsPolynomial prod_as,
                          a.Mul(secret_key.Key(), main_moduli));
    RLWE_ASSIGN_OR_RETURN(RnsPolynomial neg_prod, prod_as.Negate(main_moduli));
    RLWE_RETURN_IF_ERROR(neg_prod.AddInPlace(e, main_moduli));

    // c1 = a + s * g_k
    RLWE_RETURN_IF_ERROR(a.AddInPlace(s_gk, main_moduli));

    std::vector<RnsPolynomial> comps;
    comps.push_back(std::move(neg_prod));
    comps.push_back(std::move(a));

    rgsw_rows.emplace_back(std::move(comps), main_moduli, /*power_of_s=*/1,
                           /*error=*/0.0, error_params.get(),
                           rns_context.get());
  }

  RgswCiphertext rgsw_secret_key(std::move(rgsw_rows));

  auto sk_poly = secret_key.Key();
  if (sk_poly.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(sk_poly.ConvertToCoeffForm(main_moduli));
  }
  const auto* mod_params_q0 = main_moduli[0]->ModParams();
  uint64_t q0 = main_moduli[0]->Modulus();
  uint64_t half_q0 = q0 / 2;
  std::vector<int64_t> sk_coeffs(params.poly_degree);
  for (size_t i = 0; i < params.poly_degree; ++i) {
    uint64_t val = sk_poly.Coeffs()[0][i].ExportInt(mod_params_q0);
    if (val > half_q0) {
      sk_coeffs[i] = static_cast<int64_t>(val) - static_cast<int64_t>(q0);
    } else {
      sk_coeffs[i] = static_cast<int64_t>(val);
    }
  }

  return absl::WrapUnique(new OnionPirClient(
      std::move(params), std::move(rns_context), std::move(error_params),
      std::move(conversion_gadget), std::move(expansion_gadget),
      std::move(trace_gadget), std::move(encoder), std::move(secret_key),
      std::move(galois_keys), std::move(rgsw_secret_key),
      std::move(sk_coeffs)));
}

absl::StatusOr<SerializedEvaluationKeys>
OnionPirClient::GetSerializedEvaluationKeys() const {
  SerializedEvaluationKeys proto;
  *proto.mutable_params() = params_.ToProto();
  for (const auto& gkey : galois_keys_) {
    RLWE_ASSIGN_OR_RETURN(*proto.add_galois_keys(), gkey.Serialize());
  }
  RLWE_ASSIGN_OR_RETURN(*proto.mutable_rgsw_secret_key(),
                        rgsw_secret_key_.Serialize());
  return proto;
}

absl::StatusOr<std::pair<RnsPolynomial, std::string>>
OnionPirClient::CreatePackedQuerySeedComponents(size_t index) const {
  if (index >= params_.num_entries) {
    return absl::InvalidArgumentError(absl::StrCat(
        "Target index ", index, " >= num_entries ", params_.num_entries));
  }

  auto moduli = rns_context_->MainPrimeModuli();

  size_t n = params_.poly_degree;
  size_t n0 = params_.n0;
  size_t d = params_.num_dimensions();
  size_t ell = conversion_gadget_->Dimension();
  size_t query_w = params_.query_w();
  size_t stride = n / query_w;
  size_t tree_height = params_.query_expansion_tree_height();
  rlwe::Uint64 t = params_.plaintext_modulus();

  // Coordinate decomposition:
  size_t m0 = index % n0;
  size_t high_dim_coord = index / n0;

  // Algorithm 1:
  // 1. Create zero plaintext pt across all RNS moduli using CoefficientEncoder.
  // Pre-scale message by N^-1 mod t so that Trace (factor N/query_w) *
  // Expansion (factor query_w) = factor N cancels out!
  std::vector<rlwe::Uint64> messages(params_.poly_degree, 0);
  rlwe::Uint64 inv_n_t = ModularInverse(n, t);
  size_t rev_m0 = BitReverse(static_cast<uint32_t>(m0), tree_height);
  messages[rev_m0 * stride] = inv_n_t;

  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial pt_poly,
      encoder_->EncodeBfv(messages, moduli, /*is_scaled=*/true));
  if (!pt_poly.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(FastConvertToNttForm(pt_poly, moduli));
  }

  // 2. Sample PRG seed for uniform 'a' component
  RLWE_ASSIGN_OR_RETURN(std::string prng_seed, rlwe::HkdfPrng::GenerateSeed());
  RLWE_ASSIGN_OR_RETURN(auto prng, rlwe::HkdfPrng::Create(prng_seed));

  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial a,
      RnsPolynomial::SampleUniform(rns_context_->LogN(), prng.get(), moduli));

  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial e,
      rlwe::SampleError<ModularInt>(rns_context_->LogN(), params_.variance,
                                    moduli, prng.get()));

  RLWE_ASSIGN_OR_RETURN(RnsPolynomial prod_as,
                        a.Mul(secret_key_.Key(), moduli));
  RLWE_ASSIGN_OR_RETURN(RnsPolynomial c0, prod_as.Negate(moduli));
  RLWE_RETURN_IF_ERROR(c0.AddInPlace(e, moduli));
  RLWE_RETURN_IF_ERROR(c0.AddInPlace(pt_poly, moduli));

  // Convert c0 to coefficient form to add RGSW gadget parts for higher
  // dimensions
  RLWE_RETURN_IF_ERROR(FastConvertToCoeffForm(c0, moduli));
  std::vector<std::vector<ModularInt>> c0_all_coeffs = c0.Coeffs();

  // 5. For each higher dimension i = 1 to d - 1 and each modulus limb:
  for (size_t i = 1; i < d; ++i) {
    size_t bit_shift = d - 1 - i;
    size_t m_i = (high_dim_coord >> bit_shift) & 1;
    if (m_i == 1) {
      for (size_t k = 0; k < ell; ++k) {
        const auto& g_k_rns = conversion_gadget_->Component(k);
        size_t idx_in_pack = n0 + ell * (i - 1) + k;
        size_t p = BitReverse(static_cast<uint32_t>(idx_in_pack), tree_height);
        size_t coeff_idx = p * stride;

        for (size_t m_idx = 0; m_idx < moduli.size(); ++m_idx) {
          const auto* prime_mod = moduli[m_idx];
          const auto* m_mod_params = prime_mod->ModParams();
          rlwe::Uint64 q_i = prime_mod->Modulus();
          rlwe::Uint64 inv_n_qi = ModularInverse(n, q_i);

          rlwe::Uint64 g_k_val =
              g_k_rns.Component(m_idx).ExportInt(m_mod_params);
          rlwe::Uint64 g_prime_k = absl::Uint128Low64(
              (static_cast<absl::uint128>(g_k_val) * inv_n_qi) % q_i);

          auto& c0_coeffs = c0_all_coeffs[m_idx];
          rlwe::Uint64 curr = c0_coeffs[coeff_idx].ExportInt(m_mod_params);
          rlwe::Uint64 updated = (curr + g_prime_k) % q_i;
          RLWE_ASSIGN_OR_RETURN(c0_coeffs[coeff_idx],
                                ModularInt::ImportInt(updated, m_mod_params));
        }
      }
    }
  }

  // Convert c0 back to NTT form
  RLWE_ASSIGN_OR_RETURN(
      c0, RnsPolynomial::Create(std::move(c0_all_coeffs), /*is_ntt=*/false));
  RLWE_RETURN_IF_ERROR(FastConvertToNttForm(c0, moduli));

  return std::make_pair(std::move(c0), std::move(prng_seed));
}

absl::StatusOr<RnsBfvCiphertext> OnionPirClient::CreatePackedQueryCiphertext(
    size_t index) const {
  RLWE_ASSIGN_OR_RETURN(auto components,
                        CreatePackedQuerySeedComponents(index));
  auto c0 = std::move(components.first);
  const auto& prng_seed = components.second;

  auto moduli = rns_context_->MainPrimeModuli();
  RLWE_ASSIGN_OR_RETURN(auto prng, rlwe::HkdfPrng::Create(prng_seed));
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial a,
      RnsPolynomial::SampleUniform(rns_context_->LogN(), prng.get(), moduli));

  std::vector<RnsPolynomial> comps;
  comps.reserve(2);
  comps.push_back(std::move(c0));
  comps.push_back(std::move(a));

  return RnsBfvCiphertext(std::move(comps), moduli, /*power_of_s=*/1,
                          /*error=*/0.0, error_params_.get(),
                          rns_context_.get());
}

absl::StatusOr<SerializedPirQuery> OnionPirClient::CreateQuery(
    size_t index) const {
  RLWE_ASSIGN_OR_RETURN(auto components,
                        CreatePackedQuerySeedComponents(index));
  auto c0 = std::move(components.first);
  auto moduli = rns_context_->MainPrimeModuli();
  if (c0.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(FastConvertToCoeffForm(c0, moduli));
  }

  size_t n = params_.poly_degree;
  size_t query_w = params_.query_w();
  size_t stride = n / query_w;

  SerializedPirQuery query;
  query.set_prng_seed(std::move(components.second));
  query.set_num_coeffs(query_w);

  std::string packed_bytes;
  for (size_t m_idx = 0; m_idx < moduli.size(); ++m_idx) {
    const auto* prime_mod = moduli[m_idx];
    const auto* m_mod_params = prime_mod->ModParams();
    uint64_t q_i = prime_mod->Modulus();
    size_t num_bits = NumBitsForModulus(q_i);
    const auto& coeffs = c0.Coeffs()[m_idx];

    std::vector<uint64_t> limb_vals(query_w);
    for (size_t l = 0; l < query_w; ++l) {
      size_t coeff_idx = l * stride;
      limb_vals[l] = coeffs[coeff_idx].ExportInt(m_mod_params);
    }
    std::string limb_packed = BitPack(limb_vals, num_bits);
    packed_bytes.append(limb_packed);
  }
  query.set_compressed_coeffs(std::move(packed_bytes));
  return query;
}

absl::StatusOr<std::vector<uint8_t>> OnionPirClient::DecryptResponseCiphertexts(
    absl::Span<const RnsBfvCiphertext> response_cts,
    size_t entry_size_bytes) const {
  std::vector<uint8_t> result;
  result.reserve(entry_size_bytes);

  for (const auto& ct : response_cts) {
    rlwe::RnsRlweSecretKey<ModularInt> sk = secret_key_;
    while (sk.Level() > ct.Level()) {
      RLWE_RETURN_IF_ERROR(sk.ModReduce());
    }

    RLWE_ASSIGN_OR_RETURN(auto decrypted_messages,
                          sk.DecryptBfv(ct, encoder_.get()));

    for (size_t j = 0; j < decrypted_messages.size(); ++j) {
      if (result.size() >= entry_size_bytes) {
        break;
      }
      result.push_back(static_cast<uint8_t>(decrypted_messages[j] & 0xFF));
    }
  }

  if (result.size() < entry_size_bytes) {
    result.resize(entry_size_bytes, 0);
  }
  return result;
}

absl::StatusOr<std::vector<uint8_t>> OnionPirClient::DecodeResponse(
    const SerializedPirResponse& response, size_t entry_size_bytes) const {
  if (response.scaled_ct().empty()) {
    return absl::InvalidArgumentError(
        "Response contains no scaled ciphertexts.");
  }

  std::vector<uint8_t> result;
  result.reserve(entry_size_bytes);
  size_t n = params_.poly_degree;
  uint64_t t = params_.plaintext_modulus();

  for (const auto& scaled_ct : response.scaled_ct()) {
    uint64_t p = scaled_ct.modulus();
    int64_t p_val = static_cast<int64_t>(p);
    size_t num_bits = NumBitsForModulus(p);

    RLWE_ASSIGN_OR_RETURN(std::vector<uint64_t> c0,
                          BitUnpack(scaled_ct.c0(), n, num_bits));
    RLWE_ASSIGN_OR_RETURN(std::vector<uint64_t> c1,
                          BitUnpack(scaled_ct.c1(), n, num_bits));

    for (size_t k = 0; k < n; ++k) {
      int64_t sum = static_cast<int64_t>(c0[k]);
      for (size_t j = 0; j <= k; ++j) {
        sum += static_cast<int64_t>(c1[j]) * sk_coeffs_[k - j];
      }
      for (size_t j = k + 1; j < n; ++j) {
        sum -= static_cast<int64_t>(c1[j]) * sk_coeffs_[n + k - j];
      }
      int64_t v = (sum % p_val + p_val) % p_val;
      uint64_t m =
          static_cast<uint64_t>(std::lround((static_cast<double>(v) * t) /
                                            static_cast<double>(p_val))) %
          t;
      if (result.size() < entry_size_bytes) {
        result.push_back(static_cast<uint8_t>(m & 0xFF));
      }
    }
  }

  if (result.size() < entry_size_bytes) {
    result.resize(entry_size_bytes, 0);
  }
  return result;
}

}  // namespace oram_okvs::pir
