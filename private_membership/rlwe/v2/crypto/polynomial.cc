#include "private_membership/rlwe/v2/crypto/polynomial.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "private_membership/rlwe/v2/proto/rlwe.pb.h"
#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "shell_encryption/status_macros.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "shell_encryption/constants.h"
#include "shell_encryption/dft_transformations.h"
#include "shell_encryption/int256.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/rns/crt_interpolation.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/status_macros.h"

namespace private_membership {
namespace rlwe {
namespace v2 {

namespace {

// Check if a number is a power of 2.
bool IsPowerOfTwo(size_t n) { return (n != 0) && ((n & (n - 1)) == 0); }

}  // namespace

absl::StatusOr<Context> Context::Create(int log_n) {
  Context params_set;

  // Use three NTT friendly primes.
  std::vector<uint64_t> primes = {::rlwe::kModulus59, ::rlwe::kModulus62,
                                  17592186028033ULL};

  for (uint64_t prime : primes) {
    RLWE_ASSIGN_OR_RETURN(
        auto mod_params,
        ::rlwe::MontgomeryInt<uint64_t>::Params::Create(prime));
    auto mod_params_ptr =
        std::make_unique<const ::rlwe::MontgomeryInt<uint64_t>::Params>(
            std::move(*mod_params));
    RLWE_ASSIGN_OR_RETURN(
        auto ntt_param,
        ::rlwe::InitializeNttParameters<::rlwe::MontgomeryInt<uint64_t>>(
            log_n, mod_params_ptr.get()));

    params_set.ntt_params.push_back(std::move(ntt_param));

    auto modulus = std::make_unique<
        ::rlwe::PrimeModulus<::rlwe::MontgomeryInt<uint64_t>>>();
    modulus->mod_params = std::move(mod_params_ptr);
    // modulus->ntt_params is not used.

    params_set.moduli.push_back(std::move(modulus));
  }

  int num_primes = params_set.moduli.size();
  std::vector<const ::rlwe::PrimeModulus<::rlwe::MontgomeryInt<uint64_t>>*>
      moduli_ptrs(num_primes);
  for (int idx = 0; idx < num_primes; ++idx) {
    moduli_ptrs[idx] = params_set.moduli[idx].get();
  }

  params_set.modulus_hats =
      ::rlwe::RnsModulusComplements<::rlwe::MontgomeryInt<uint64_t>,
                                    ::rlwe::uint256>(moduli_ptrs);

  for (int idx = 0; idx < num_primes; ++idx) {
    auto mod_params_i = params_set.moduli[idx]->ModParams();
    uint64_t hat_mod_i = static_cast<uint64_t>(
        params_set.modulus_hats[idx] % params_set.moduli[idx]->Modulus());
    RLWE_ASSIGN_OR_RETURN(
        auto hat_mod_i_ntt,
        ::rlwe::MontgomeryInt<uint64_t>::ImportInt(hat_mod_i, mod_params_i));
    params_set.modulus_hat_invs.push_back(
        hat_mod_i_ntt.MultiplicativeInverse(mod_params_i));
  }

  ::rlwe::uint256 q = 1;
  for (int idx = 0; idx < num_primes; ++idx) {
    q *= ::rlwe::uint256(params_set.moduli[idx]->Modulus());
  }
  params_set.q = q;
  params_set.q_half = q >> 1;

  return std::move(params_set);
}


absl::StatusOr<NttPolynomial> NttPolynomial::CreateZero(int len,
                                                        const Context& ctx) {
  if (len <= 0 || !IsPowerOfTwo(len)) {
    return absl::InvalidArgumentError(
        "Polynomial length must be a non-zero power of two.");
  }
  std::vector<std::vector<::rlwe::MontgomeryInt<uint64_t>>> ntt_polys;
  ntt_polys.reserve(ctx.moduli.size());
  for (int i = 0; i < ctx.moduli.size(); ++i) {
    auto mod_params_i = ctx.moduli[i]->ModParams();
    RLWE_ASSIGN_OR_RETURN(
        auto zero, ::rlwe::MontgomeryInt<uint64_t>::ImportInt(0, mod_params_i));
    ntt_polys.push_back(
        std::vector<::rlwe::MontgomeryInt<uint64_t>>(len, zero));
  }
  return NttPolynomial(std::move(ntt_polys));
}

absl::StatusOr<NttPolynomial> NttPolynomial::Add(const NttPolynomial& that,
                                                 const Context& ctx) const {
  if (ntt_polys_.size() != that.ntt_polys_.size()) {
    return absl::InvalidArgumentError(
        "Polynomials must have same number of CRT components.");
  }
  std::vector<std::vector<::rlwe::MontgomeryInt<uint64_t>>> result_ntt_polys;
  result_ntt_polys.reserve(ntt_polys_.size());
  for (int i = 0; i < ntt_polys_.size(); ++i) {
    if (ntt_polys_[i].size() != that.ntt_polys_[i].size()) {
      return absl::InvalidArgumentError("Polynomial lengths must match.");
    }
    auto mod_params_i = ctx.moduli[i]->ModParams();
    std::vector<::rlwe::MontgomeryInt<uint64_t>> sum_ntt;
    sum_ntt.reserve(ntt_polys_[i].size());
    for (int j = 0; j < ntt_polys_[i].size(); ++j) {
      sum_ntt.push_back(
          ntt_polys_[i][j].Add(that.ntt_polys_[i][j], mod_params_i));
    }
    result_ntt_polys.push_back(std::move(sum_ntt));
  }
  return NttPolynomial(std::move(result_ntt_polys));
}

absl::StatusOr<NttPolynomial> NttPolynomial::Mult(const NttPolynomial& that,
                                                  const Context& ctx) const {
  if (ntt_polys_.size() != that.ntt_polys_.size()) {
    return absl::InvalidArgumentError(
        "Polynomials must have same number of CRT components.");
  }
  std::vector<std::vector<::rlwe::MontgomeryInt<uint64_t>>> result_ntt_polys;
  result_ntt_polys.reserve(ntt_polys_.size());
  for (int i = 0; i < ntt_polys_.size(); ++i) {
    if (ntt_polys_[i].size() != that.ntt_polys_[i].size()) {
      return absl::InvalidArgumentError("Polynomial lengths must match.");
    }
    auto mod_params_i = ctx.moduli[i]->ModParams();
    std::vector<::rlwe::MontgomeryInt<uint64_t>> prod_ntt = ntt_polys_[i];
    RLWE_RETURN_IF_ERROR(::rlwe::MontgomeryInt<uint64_t>::BatchMulInPlace(
        &prod_ntt, that.ntt_polys_[i], mod_params_i));
    result_ntt_polys.push_back(std::move(prod_ntt));
  }
  return NttPolynomial(std::move(result_ntt_polys));
}

template <typename CoeffType>
absl::StatusOr<Polynomial<CoeffType>> Polynomial<CoeffType>::Create(
    std::vector<CoeffType> coeffs) {
  if (coeffs.empty() || !IsPowerOfTwo(coeffs.size())) {
    return absl::InvalidArgumentError(
        "Polynomial length must be a non-zero power of two.");
  }
  return Polynomial(std::move(coeffs));
}

template <typename CoeffType>
absl::StatusOr<Polynomial<CoeffType>> Polynomial<CoeffType>::CreateZero(
    int len) {
  if (len <= 0 || !IsPowerOfTwo(len)) {
    return absl::InvalidArgumentError(
        "Polynomial length must be a non-zero power of two.");
  }
  return Polynomial(std::vector<CoeffType>(len, 0));
}

template <typename CoeffType>
absl::StatusOr<Polynomial<CoeffType>> Polynomial<CoeffType>::CreateFromProto(
    const proto::Polynomial& proto, int len, int log_modulus) {
  if (len <= 0 || !IsPowerOfTwo(len)) {
    return absl::InvalidArgumentError("Invalid length.");
  }
  if (log_modulus <= 0 || log_modulus > 8 * sizeof(CoeffType)) {
    return absl::InvalidArgumentError("Invalid log_modulus.");
  }
  if (proto.encoded_coeffs().empty()) {
    return absl::InvalidArgumentError("Empty serialized polynomial.");
  }

  absl::string_view bytes = proto.encoded_coeffs();
  size_t expected_bits = static_cast<size_t>(len) * log_modulus;
  size_t expected_bytes = (expected_bits + 7) / 8;
  if (bytes.size() != expected_bytes) {
    return absl::InvalidArgumentError("Unexpected serialized size.");
  }

  std::vector<CoeffType> coeffs;
  coeffs.reserve(len);

  size_t byte_idx = 0;
  int bit_offset = 0;

  for (int j = 0; j < len; ++j) {
    CoeffType coeff = 0;
    for (int i = 0; i < log_modulus; ++i) {
      if (byte_idx >= bytes.size()) {
        return absl::InvalidArgumentError("Serialized polynomial too small.");
      }
      uint8_t current_byte = bytes[byte_idx];
      CoeffType bit = (current_byte >> bit_offset) & 1;
      coeff |= (bit << i);
      bit_offset++;
      if (bit_offset == 8) {
        bit_offset = 0;
        byte_idx++;
      }
    }
    coeffs.push_back(coeff);
  }
  return Polynomial::Create(std::move(coeffs));
}

template <typename CoeffType>
absl::StatusOr<Polynomial<CoeffType>> Polynomial<CoeffType>::Add(
    const Polynomial& that) const {
  if (Len() != that.Len()) {
    return absl::InvalidArgumentError("Polynomial lengths must match.");
  }
  std::vector<CoeffType> result(Len());
  for (int i = 0; i < Len(); ++i) {
    result[i] = coeffs_[i] + that.coeffs_[i];
  }
  return Polynomial::Create(std::move(result));
}

template <typename CoeffType>
absl::Status Polynomial<CoeffType>::AddInPlace(const Polynomial& that) {
  if (Len() != that.Len()) {
    return absl::InvalidArgumentError("Polynomial lengths must match.");
  }
  for (int i = 0; i < Len(); ++i) {
    coeffs_[i] += that.coeffs_[i];
  }
  return absl::OkStatus();
}

template <typename CoeffType>
absl::StatusOr<Polynomial<CoeffType>> Polynomial<CoeffType>::Negate() const {
  std::vector<CoeffType> result(Len());
  for (int i = 0; i < Len(); ++i) {
    result[i] = -coeffs_[i];
  }
  return Polynomial::Create(std::move(result));
}

template <typename CoeffType>
absl::StatusOr<Polynomial<CoeffType>> Polynomial<CoeffType>::Rescale(
    int log_current_modulus, int log_new_modulus) const {
  std::vector<CoeffType> result(Len());

  absl::uint128 current_modulus = static_cast<absl::uint128>(CoeffType{1})
                                  << log_current_modulus;
  absl::uint128 new_modulus = static_cast<absl::uint128>(CoeffType{1})
                              << log_new_modulus;
  absl::uint128 half_current_modulus = current_modulus / 2;
  for (int i = 0; i < Len(); ++i) {
    result[i] = static_cast<CoeffType>(
        (static_cast<absl::uint128>(coeffs_[i]) * new_modulus +
         half_current_modulus) /
        current_modulus);
  }
  return Polynomial::Create(std::move(result));
}

template <typename CoeffType>
absl::StatusOr<Polynomial<CoeffType>> Polynomial<CoeffType>::RightShift(
    int amount) const {
  if (amount < 0) {
    return absl::InvalidArgumentError("Shift amount cannot be negative.");
  }
  if (amount >= 8 * sizeof(CoeffType)) {
    std::vector<CoeffType> result(Len(), 0);
    return Polynomial::Create(std::move(result));
  }
  std::vector<CoeffType> result(Len());
  for (int i = 0; i < Len(); ++i) {
    result[i] = coeffs_[i] >> amount;
  }
  return Polynomial::Create(std::move(result));
}

template <typename CoeffType>
absl::StatusOr<Polynomial<CoeffType>> Polynomial<CoeffType>::LowBits(
    int amount) const {
  if (amount < 0) {
    return absl::InvalidArgumentError("Amount cannot be negative.");
  }
  if (amount >= 8 * sizeof(CoeffType)) {
    return Polynomial::Create(coeffs_);
  }
  std::vector<CoeffType> result(Len());
  CoeffType mask = (CoeffType{1} << amount) - 1;
  for (int i = 0; i < Len(); ++i) {
    result[i] = coeffs_[i] & mask;
  }
  return Polynomial::Create(std::move(result));
}

template <typename CoeffType>
absl::StatusOr<std::vector<Polynomial<CoeffType>>>
Polynomial<CoeffType>::GadgetInv(int log_modulus, int log_digit,
                                 int num_digits) const {
  if (log_digit <= 0 || num_digits <= 0) {
    return absl::InvalidArgumentError(
        "log_digit and num_digits must be positive.");
  }

  if (log_modulus == 0) {
    log_modulus = 8 * sizeof(CoeffType);
  }

  int to_shift = log_modulus - (num_digits * log_digit);
  Polynomial<CoeffType> current;
  if (to_shift > 0) {
    std::vector<CoeffType> rounded_coeffs(Len());
    CoeffType round_val = CoeffType{1} << (to_shift - 1);
    for (int i = 0; i < Len(); ++i) {
      rounded_coeffs[i] = static_cast<CoeffType>(
          (static_cast<absl::uint128>(coeffs_[i]) + round_val) >> to_shift);
    }
    RLWE_ASSIGN_OR_RETURN(current, Polynomial<CoeffType>::Create(rounded_coeffs));
  } else {
    current = *this;
  }
  std::vector<Polynomial<CoeffType>> result;
  result.reserve(num_digits);
  for (int i = 0; i < num_digits; ++i) {
    RLWE_ASSIGN_OR_RETURN(auto lower, current.LowBits(log_digit));
    result.push_back(std::move(lower));
    RLWE_ASSIGN_OR_RETURN(current, current.RightShift(log_digit));
  }
  return result;
}

template <typename CoeffType>
absl::StatusOr<std::vector<Polynomial<CoeffType>>>
Polynomial<CoeffType>::SignedGadgetInv(int log_modulus, int log_digit,
                                       int num_digits) const {
  if (log_digit <= 0 || num_digits <= 0) {
    return absl::InvalidArgumentError(
        "log_digit and num_digits must be positive.");
  }

  if (log_modulus == 0) {
    log_modulus = 8 * sizeof(CoeffType);
  }

  using SignedCoeffType = std::make_signed_t<CoeffType>;

  int to_shift = log_modulus - (num_digits * log_digit);
  std::vector<CoeffType> current_coeffs = coeffs_;
  if (to_shift > 0) {
    CoeffType round_val = CoeffType{1} << (to_shift - 1);
    for (int i = 0; i < Len(); ++i) {
      current_coeffs[i] = static_cast<CoeffType>(
          (static_cast<absl::uint128>(current_coeffs[i]) + round_val) >>
          to_shift);
    }
  }

  std::vector<Polynomial<CoeffType>> result;
  result.reserve(num_digits);

  CoeffType mask = (CoeffType{1} << log_digit) - 1;
  CoeffType half_mask = CoeffType{1} << (log_digit - 1);

  for (int i = 0; i < num_digits; ++i) {
    std::vector<CoeffType> digit_coeffs(Len());
    for (int j = 0; j < Len(); ++j) {
      CoeffType lower = current_coeffs[j] & mask;
      if (lower >= half_mask) {
        digit_coeffs[j] = static_cast<CoeffType>(
            static_cast<SignedCoeffType>(lower) -
            static_cast<SignedCoeffType>(CoeffType{1} << log_digit));
        current_coeffs[j] = (current_coeffs[j] >> log_digit) + 1;
      } else {
        digit_coeffs[j] = lower;
        current_coeffs[j] = (current_coeffs[j] >> log_digit);
      }
    }
    RLWE_ASSIGN_OR_RETURN(auto lower_poly,
                     Polynomial<CoeffType>::Create(std::move(digit_coeffs)));
    result.push_back(std::move(lower_poly));
  }
  return result;
}

template <typename CoeffType>
absl::StatusOr<Polynomial<CoeffType>> Polynomial<CoeffType>::Automorph(
    int power) const {
  if (power < 0 || power % 2 == 0 || power >= 2 * Len()) {
    return absl::InvalidArgumentError(
        "Automorphism power must be a non-negative odd integer less than 2*N.");
  }
  std::vector<CoeffType> result(Len());
  // Automorphism maps X to X^power
  // p(X) = sum a_i X^i
  // p(X^power) = sum a_i X^{i * power}
  // Since X^N = -1, X^{i * power} = (-1)^k X^j
  // where i * power = k * N + j, with 0 <= j < N.
  for (int i = 0; i < Len(); ++i) {
    int index = (i * power) % (2 * Len());
    int sign = 1;
    if (index >= Len()) {
      sign = -1;
      index -= Len();
    }
    if (sign == 1) {
      result[index] += coeffs_[i];
    } else {
      result[index] -= coeffs_[i];
    }
  }
  return Polynomial::Create(std::move(result));
}

template <typename CoeffType>
absl::StatusOr<Polynomial<CoeffType>> Polynomial<CoeffType>::Mult(
    const Polynomial& that, const Context& ctx) const {
  if (Len() != that.Len()) {
    return absl::InvalidArgumentError("Polynomial lengths must match.");
  }

  RLWE_ASSIGN_OR_RETURN(NttPolynomial this_ntt, ToNtt(ctx));
  RLWE_ASSIGN_OR_RETURN(NttPolynomial that_ntt, that.ToNtt(ctx));
  RLWE_ASSIGN_OR_RETURN(NttPolynomial result_ntt, this_ntt.Mult(that_ntt, ctx));
  return FromNtt(result_ntt, ctx);
}

template <typename CoeffType>
absl::StatusOr<NttPolynomial> Polynomial<CoeffType>::ToNtt(
    const Context& ctx) const {
  std::vector<std::vector<::rlwe::MontgomeryInt<uint64_t>>> ntt_polys;
  ntt_polys.reserve(ctx.moduli.size());
  for (int i = 0; i < ctx.moduli.size(); ++i) {
    auto mod_params_i = ctx.moduli[i]->ModParams();
    std::vector<::rlwe::MontgomeryInt<uint64_t>> ntt_poly;
    ntt_poly.reserve(Len());
    for (int j = 0; j < Len(); ++j) {
      RLWE_ASSIGN_OR_RETURN(
          auto val, ::rlwe::MontgomeryInt<uint64_t>::ImportInt(
                        static_cast<uint64_t>(coeffs_[j]), mod_params_i));
      ntt_poly.push_back(std::move(val));
    }
    RLWE_RETURN_IF_ERROR(::rlwe::ForwardNumberTheoreticTransform(
        ntt_poly, ctx.ntt_params[i], *mod_params_i));
    ntt_polys.push_back(std::move(ntt_poly));
  }
  return NttPolynomial(std::move(ntt_polys));
}

template <typename CoeffType>
absl::StatusOr<Polynomial<CoeffType>> Polynomial<CoeffType>::FromNtt(
    const NttPolynomial& ntt_poly, const Context& ctx) {
  int num_primes = ctx.moduli.size();
  if (ntt_poly.ntt_polys_.size() != num_primes) {
    return absl::InvalidArgumentError(
        "NTT polynomial must have the same number of CRT components as the "
        "context.");
  }
  int len = ntt_poly.ntt_polys_[0].size();

  std::vector<std::vector<::rlwe::MontgomeryInt<uint64_t>>> crt_components;
  crt_components.reserve(num_primes);

  std::vector<const ::rlwe::PrimeModulus<::rlwe::MontgomeryInt<uint64_t>>*>
      moduli_ptrs(num_primes);
  for (int i = 0; i < num_primes; ++i) {
    moduli_ptrs[i] = ctx.moduli[i].get();
  }

  for (int i = 0; i < num_primes; ++i) {
    auto mod_params_i = ctx.moduli[i]->ModParams();
    std::vector<::rlwe::MontgomeryInt<uint64_t>> inv_ntt =
        ntt_poly.ntt_polys_[i];
    RLWE_RETURN_IF_ERROR(::rlwe::InverseNumberTheoreticTransform(
        inv_ntt, ctx.ntt_params[i], *mod_params_i));
    crt_components.push_back(std::move(inv_ntt));
  }

  RLWE_ASSIGN_OR_RETURN(
      auto crt_result,
      (::rlwe::CrtInterpolation<::rlwe::MontgomeryInt<uint64_t>,
                                ::rlwe::uint256>(crt_components, moduli_ptrs,
                                                 ctx.modulus_hats,
                                                 ctx.modulus_hat_invs)));

  ::rlwe::uint256 q = ctx.q;
  ::rlwe::uint256 q_half = ctx.q_half;

  std::vector<CoeffType> result(len);
  for (int j = 0; j < len; ++j) {
    ::rlwe::uint256 val = crt_result[j];
    if (val > q_half) {
      // It's a negative value. We want to convert (val - q) mod 2^64 to
      // CoeffType. Since val - q is negative, we can just do (val mod 2^64) -
      // (q mod 2^64)
      uint64_t val_low = absl::Uint128Low64(::rlwe::Uint256Low128(val));
      uint64_t q_low = absl::Uint128Low64(::rlwe::Uint256Low128(q));
      result[j] = static_cast<CoeffType>(val_low - q_low);
    } else {
      result[j] = static_cast<CoeffType>(
          absl::Uint128Low64(::rlwe::Uint256Low128(val)));
    }
  }

  return Polynomial::Create(std::move(result));
}

template <typename CoeffType>
absl::StatusOr<proto::Polynomial> Polynomial<CoeffType>::ToProto(
    int log_modulus) const {
  if (log_modulus <= 0 || log_modulus > 8 * sizeof(CoeffType)) {
    return absl::InvalidArgumentError("Invalid log_modulus.");
  }
  std::string bytes;
  uint8_t current_byte = 0;
  int bits_in_byte = 0;
  for (CoeffType coeff : coeffs_) {
    for (int i = 0; i < log_modulus; ++i) {
      uint8_t bit = (coeff >> i) & 1;
      current_byte |= (bit << bits_in_byte);
      bits_in_byte++;
      if (bits_in_byte == 8) {
        bytes.push_back(current_byte);
        current_byte = 0;
        bits_in_byte = 0;
      }
    }
  }
  if (bits_in_byte > 0) {
    bytes.push_back(current_byte);
  }

  proto::Polynomial result;
  result.set_encoded_coeffs(std::move(bytes));
  return result;
}


template class Polynomial<uint32_t>;
template class Polynomial<uint64_t>;

}  // namespace v2
}  // namespace rlwe
}  // namespace private_membership
