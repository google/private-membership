#ifndef PRIVACY_PRIVATE_MEMBERSHIP_RLWE_V2_CRYPTO_ENCRYPTION_H_
#define PRIVACY_PRIVATE_MEMBERSHIP_RLWE_V2_CRYPTO_ENCRYPTION_H_

#include <vector>

#include "private_membership/rlwe/v2/crypto/polynomial.h"
#include "absl/status/statusor.h"
#include "shell_encryption/prng/prng.h"

namespace private_membership {
namespace rlwe {
namespace v2 {

// Parameters for Ring Learning with Errors (RLWE).
template <typename CoeffType>
class RlweParams {
 public:
  static absl::StatusOr<RlweParams<CoeffType>> Create(int degree,
                                                      CoeffType modulus,
                                                      int variance);
  static absl::StatusOr<RlweParams<CoeffType>> Create(
      int degree, CoeffType modulus, CoeffType modulus1_after_switch,
      CoeffType modulus2_after_switch, int variance);

  // Returns the ring degree.
  int Degree() const { return degree_; }

  // Returns the RLWE modulus.
  CoeffType Modulus() const { return modulus_; }

  // Returns the bit size of the RLWE modulus.
  int LogModulus() const {
    return modulus_ == 0 ? sizeof(CoeffType) * 8
                         : absl::bit_width(modulus_ - 1);
  }

  CoeffType Modulus1AfterSwitch() const { return modulus_1_after_switch_; }
  CoeffType Modulus2AfterSwitch() const { return modulus_2_after_switch_; }

  int LogModulus1AfterSwitch() const {
    return modulus_1_after_switch_ == 0
               ? sizeof(CoeffType) * 8
               : absl::bit_width(modulus_1_after_switch_ - 1);
  }

  int LogModulus2AfterSwitch() const {
    return modulus_2_after_switch_ == 0
               ? sizeof(CoeffType) * 8
               : absl::bit_width(modulus_2_after_switch_ - 1);
  }

  // Variance of the error distribution.
  int Variance() const { return variance_; }

 private:
  RlweParams(int degree, CoeffType modulus, int variance)
      : degree_(degree), modulus_(modulus), modulus_1_after_switch_(modulus),
        modulus_2_after_switch_(modulus), variance_(variance) {}

  RlweParams(int degree, CoeffType modulus, CoeffType modulus1_after_switch,
             CoeffType modulus2_after_switch, int variance)
      : degree_(degree),
        modulus_(modulus),
        modulus_1_after_switch_(modulus1_after_switch),
        modulus_2_after_switch_(modulus2_after_switch),
        variance_(variance) {}

  // Ring degree.
  const int degree_;
  // RLWE modulus
  // A modulus of 0 indicates that we operate in Z[X] / (X^d + 1).
  const CoeffType modulus_;

  const CoeffType modulus_1_after_switch_;
  const CoeffType modulus_2_after_switch_;

  const int variance_;
};

// Parameters for the Gadget decomposition.
struct GadgetParams {
  int log_digit;
  int num_digits;
};

// Samples a secret key from the uniform ternary distribution.
template <typename CoeffType>
absl::StatusOr<Polynomial<CoeffType>> SampleSecretKey(
    const RlweParams<CoeffType>& params, ::rlwe::SecurePrng& prng);

// Samples a vector of a components uniformly from Z_q.
//
// `num_a_components` is the number of polynomials to generate.
template <typename CoeffType>
absl::StatusOr<std::vector<Polynomial<CoeffType>>> SampleAComponents(
    const RlweParams<CoeffType>& params, int num_a_components,
    ::rlwe::SecurePrng& prng);

// Represents a single Ring Learning with Errors (RLWE) sample (a, b).
//
// `a` is a uniformly random polynomial
// `b` is computed as `-a * s + e` where `s` is the secret key and `e` is the
//   small error polynomial.
template <typename CoeffType>
struct RlweSample {
  Polynomial<CoeffType> a;
  Polynomial<CoeffType> b;
};

// Generates RLWE samples: (a_i, -a_i * s + e_i) where a_i are from
// `a_components`, s is `secret_key`, and e_i are sampled from the error
// distribution with `params.variance`.
//
// The returned vector contains pairs of (a_i, b_i), where b_i = -a_i * s + e_i.
template <typename CoeffType>
absl::StatusOr<std::vector<RlweSample<CoeffType>>> GenerateRlweSamples(
    const RlweParams<CoeffType>& params,
    const Polynomial<CoeffType>& secret_key,
    const std::vector<Polynomial<CoeffType>>& a_components,
    ::rlwe::SecurePrng& prng);

// Represents an RLWE Ciphertext (a, b).
//
// `a` is the uniformly random polynomial from the RLWE sample.
// `b` is computed as `-a * s + e + m * \Delta` where `s` is the secret key,
// `e` is the small error polynomial, `m` is the message, and `\Delta` is the
//   scaling factor.
template <typename CoeffType>
struct RlweCiphertext {
  Polynomial<CoeffType> a;
  Polynomial<CoeffType> b;
};

// Encrypts a message using a generated RLWE sample.
//
// `plaintext_modulus` is the plaintext modulus. It must be positive.
// `rlwe_sample` is an RLWE sample (a, b). Its polynomials must have length
//   `params.Degree()`.
// `message` is a vector of scalars representing the message. It must have
//   length `params.Degree()` and its elements must be in
//   [0, plaintext_modulus-1].
//
// The returned RlweCiphertext (A, B) is a ciphertext
// under secret key `s` for `message`.
template <typename CoeffType>
absl::StatusOr<RlweCiphertext<CoeffType>> EncryptFromRlweSample(
    const RlweParams<CoeffType>& params, CoeffType plaintext_modulus,
    const RlweSample<CoeffType>& rlwe_sample,
    const std::vector<CoeffType>& message);

// Encrypts a message.
//
// `plaintext_modulus` is the plaintext modulus. It must be positive.
// `message` is a vector of scalars representing the message. It must have
//   length `params.Degree()` and its elements must be in
//   [0, plaintext_modulus-1].
//
// The returned RlweCiphertext (A, B) is a ciphertext
// under secret key `s` for `message`.
template <typename CoeffType>
absl::StatusOr<RlweCiphertext<CoeffType>> EncryptRlwe(
    const RlweParams<CoeffType>& params, CoeffType plaintext_modulus,
    const Polynomial<CoeffType>& secret_key,
    const std::vector<CoeffType>& message, ::rlwe::SecurePrng& prng);

// Decrypts a ciphertext (A, B) using the secret key `secret_key`.
//
// `plaintext_modulus` is the plaintext modulus. It must be positive.
// `ciphertext` is the ciphertext polynomials (a, b). They must have length
//   `params.Degree()`.
// `secret_key` is the secret key polynomial. It must have length
//   `params.Degree()`.
//
// The returned vector of scalars has length `params.Degree()` and its elements
// are in [0, plaintext_modulus-1].
template <typename CoeffType>
absl::StatusOr<std::vector<CoeffType>> Decrypt(
    const RlweParams<CoeffType>& params, CoeffType plaintext_modulus,
    const RlweCiphertext<CoeffType>& ciphertext,
    const Polynomial<CoeffType>& secret_key);

template <typename CoeffType>
absl::StatusOr<std::vector<CoeffType>> DecryptAfterModulusSwitch(
    const RlweParams<CoeffType>& params, CoeffType plaintext_modulus,
    const RlweCiphertext<CoeffType>& ciphertext,
    const Polynomial<CoeffType>& secret_key);

// Encrypts a gadget ciphertext based on a set of independent RLWE samples.
//
// This function constructs a gadget ciphertext essentially by forming
// pseudo-RLWE ciphertexts encrypting the base-b decomposed components
// of `message`, combined with the RLWE samples.
//
// The output is a gadget ciphertext containing `num_components * num_digits`
// generated RLWE ciphertexts, where each has the source material `message`
// encoded in its message component.
template <typename CoeffType>
absl::StatusOr<std::vector<RlweCiphertext<CoeffType>>>
EncryptGadgetCiphertextFromRlweSamples(
    const RlweParams<CoeffType>& params, const GadgetParams& gadget_params,
    const std::vector<RlweSample<CoeffType>>& rlwe_samples,
    const std::vector<CoeffType>& message);

// RGSW Ciphertext composed of two gadget ciphertext vectors.
//
// `u` contains gadget message encoding for `m * s`.
// `v` contains gadget message encoding for `m`.
template <typename CoeffType>
struct RgswCiphertext {
  std::vector<RlweCiphertext<CoeffType>> u;
  std::vector<RlweCiphertext<CoeffType>> v;
};

// Encrypts an RGSW ciphertext from a set of RLWE samples.
//
// `u_samples` and `v_samples` are vectors of RlweSample (w_i, y'_i) where
//   y'_i = -w_i * s + e_i. They must each have length `num_digits`.
// `secret_key` is the secret key polynomial.
// `message` is the message vector. It must have length params.Degree().
//
// The returned RGSW ciphertext contains `u` and `v`.
template <typename CoeffType>
absl::StatusOr<RgswCiphertext<CoeffType>> EncryptRgswCiphertextFromRlweSamples(
    const RlweParams<CoeffType>& params, const GadgetParams& gadget_params,
    const std::vector<RlweSample<CoeffType>>& u_samples,
    const std::vector<RlweSample<CoeffType>>& v_samples,
    const Polynomial<CoeffType>& secret_key,
    const std::vector<CoeffType>& message);


extern template absl::StatusOr<RlweParams<uint32_t>>
RlweParams<uint32_t>::Create(int, uint32_t, int);
extern template absl::StatusOr<RlweParams<uint64_t>>
RlweParams<uint64_t>::Create(int, uint64_t, int);

extern template absl::StatusOr<Polynomial<uint32_t>> SampleSecretKey<uint32_t>(
    const RlweParams<uint32_t>&, ::rlwe::SecurePrng&);
extern template absl::StatusOr<Polynomial<uint64_t>> SampleSecretKey<uint64_t>(
    const RlweParams<uint64_t>&, ::rlwe::SecurePrng&);

extern template absl::StatusOr<std::vector<Polynomial<uint32_t>>>
SampleAComponents<uint32_t>(const RlweParams<uint32_t>&, int,
                            ::rlwe::SecurePrng&);
extern template absl::StatusOr<std::vector<Polynomial<uint64_t>>>
SampleAComponents<uint64_t>(const RlweParams<uint64_t>&, int,
                            ::rlwe::SecurePrng&);

extern template absl::StatusOr<std::vector<RlweSample<uint32_t>>>
GenerateRlweSamples<uint32_t>(const RlweParams<uint32_t>&,
                              const Polynomial<uint32_t>&,
                              const std::vector<Polynomial<uint32_t>>&,
                              ::rlwe::SecurePrng&);
extern template absl::StatusOr<std::vector<RlweSample<uint64_t>>>
GenerateRlweSamples<uint64_t>(const RlweParams<uint64_t>&,
                              const Polynomial<uint64_t>&,
                              const std::vector<Polynomial<uint64_t>>&,
                              ::rlwe::SecurePrng&);

extern template absl::StatusOr<RlweCiphertext<uint32_t>>
EncryptFromRlweSample<uint32_t>(const RlweParams<uint32_t>&, uint32_t,
                                const RlweSample<uint32_t>&,
                                const std::vector<uint32_t>&);
extern template absl::StatusOr<RlweCiphertext<uint64_t>>
EncryptFromRlweSample<uint64_t>(const RlweParams<uint64_t>&, uint64_t,
                                const RlweSample<uint64_t>&,
                                const std::vector<uint64_t>&);

extern template absl::StatusOr<std::vector<uint32_t>> Decrypt<uint32_t>(
    const RlweParams<uint32_t>&, uint32_t, const RlweCiphertext<uint32_t>&,
    const Polynomial<uint32_t>&);
extern template absl::StatusOr<std::vector<uint64_t>> Decrypt<uint64_t>(
    const RlweParams<uint64_t>&, uint64_t, const RlweCiphertext<uint64_t>&,
    const Polynomial<uint64_t>&);

extern template absl::StatusOr<std::vector<uint32_t>>
DecryptAfterModulusSwitch<uint32_t>(const RlweParams<uint32_t>&, uint32_t,
                                    const RlweCiphertext<uint32_t>&,
                                    const Polynomial<uint32_t>&);
extern template absl::StatusOr<std::vector<uint64_t>>
DecryptAfterModulusSwitch<uint64_t>(const RlweParams<uint64_t>&, uint64_t,
                                    const RlweCiphertext<uint64_t>&,
                                    const Polynomial<uint64_t>&);

extern template absl::StatusOr<std::vector<RlweCiphertext<uint32_t>>>
EncryptGadgetCiphertextFromRlweSamples<uint32_t>(
    const RlweParams<uint32_t>&, const GadgetParams&,
    const std::vector<RlweSample<uint32_t>>&, const std::vector<uint32_t>&);
extern template absl::StatusOr<std::vector<RlweCiphertext<uint64_t>>>
EncryptGadgetCiphertextFromRlweSamples<uint64_t>(
    const RlweParams<uint64_t>&, const GadgetParams&,
    const std::vector<RlweSample<uint64_t>>&, const std::vector<uint64_t>&);

extern template absl::StatusOr<RgswCiphertext<uint32_t>>
EncryptRgswCiphertextFromRlweSamples<uint32_t>(
    const RlweParams<uint32_t>&, const GadgetParams&,
    const std::vector<RlweSample<uint32_t>>&,
    const std::vector<RlweSample<uint32_t>>&, const Polynomial<uint32_t>&,
    const std::vector<uint32_t>&);
extern template absl::StatusOr<RgswCiphertext<uint64_t>>
EncryptRgswCiphertextFromRlweSamples<uint64_t>(
    const RlweParams<uint64_t>&, const GadgetParams&,
    const std::vector<RlweSample<uint64_t>>&,
    const std::vector<RlweSample<uint64_t>>&, const Polynomial<uint64_t>&,
    const std::vector<uint64_t>&);


}  // namespace v2
}  // namespace rlwe
}  // namespace private_membership

#endif  // PRIVACY_PRIVATE_MEMBERSHIP_RLWE_V2_CRYPTO_ENCRYPTION_H_
