#ifndef PRIVACY_PRIVATE_MEMBERSHIP_RLWE_V2_CRYPTO_POLYNOMIAL_H_
#define PRIVACY_PRIVATE_MEMBERSHIP_RLWE_V2_CRYPTO_POLYNOMIAL_H_

#include <cstdint>
#include <memory>
#include <vector>

#include "private_membership/rlwe/v2/proto/rlwe.pb.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "shell_encryption/int256.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/rns/rns_modulus.h"

namespace private_membership {
namespace rlwe {
namespace v2 {

struct FftContext;

// Context holds the precomputed NTT parameters and constants needed for
// polynomial multiplication. It is independent of the coefficient type and
// depends only on the polynomial degree (length).
struct Context {
  // Creates a new context for a given polynomial length 2^log_n.
  static absl::StatusOr<Context> Create(int log_n);

  std::vector<std::unique_ptr<
      const ::rlwe::PrimeModulus<::rlwe::MontgomeryInt<uint64_t>>>>
      moduli;
  std::vector<::rlwe::NttParameters<::rlwe::MontgomeryInt<uint64_t>>>
      ntt_params;
  std::vector<::rlwe::uint256> modulus_hats;
  std::vector<::rlwe::MontgomeryInt<uint64_t>> modulus_hat_invs;
  ::rlwe::uint256 q;
  ::rlwe::uint256 q_half;
};

// Represents a polynomial in NTT form, suitable for efficient multiplication
// and addition.
class NttPolynomial {
 public:
  // Creates a zero polynomial in NTT form.
  static absl::StatusOr<NttPolynomial> CreateZero(int len, const Context& ctx);

  NttPolynomial() = default;

  // Copy constructor.
  NttPolynomial(const NttPolynomial& p) = default;
  NttPolynomial& operator=(const NttPolynomial& that) = default;

  // Move constructor.
  NttPolynomial(NttPolynomial&& p) = default;
  NttPolynomial& operator=(NttPolynomial&& that) = default;

  // Adds two polynomials in NTT form.
  absl::StatusOr<NttPolynomial> Add(const NttPolynomial& that,
                                    const Context& ctx) const;

  // Multiplies two polynomials in NTT form.
  absl::StatusOr<NttPolynomial> Mult(const NttPolynomial& that,
                                     const Context& ctx) const;

  // Accessor for NTT polynomials.
  const std::vector<std::vector<::rlwe::MontgomeryInt<uint64_t>>>& NttPolys()
      const {
    return ntt_polys_;
  }

 private:
  explicit NttPolynomial(
      std::vector<std::vector<::rlwe::MontgomeryInt<uint64_t>>> ntt_polys)
      : ntt_polys_(std::move(ntt_polys)) {}

  // A vector of polynomials in NTT form, one for each RNS modulus.
  std::vector<std::vector<::rlwe::MontgomeryInt<uint64_t>>> ntt_polys_;

  template <typename CoeffType>
  friend class Polynomial;
};

// A polynomial in coefficient representation with operations implicitly
// modulo a power of two defined by the `CoeffType`.
// For `CoeffType` = `uint32_t`, operations are implicitly modulo 2^32.
// For `CoeffType` = `uint64_t`, operations are implicitly modulo 2^64.
template <typename CoeffType>
class Polynomial {
 public:
  // Creates a polynomial with the given coefficients.
  // The length of the coefficient vector must be a power of two.
  static absl::StatusOr<Polynomial> Create(std::vector<CoeffType> coeffs);

  // Creates a polynomial with the given length filled with zeros.
  // The length must be a power of two.
  static absl::StatusOr<Polynomial> CreateZero(int len);

  // Creates a polynomial from a proto::Polynomial for a given degree
  // (length `len`) and `log_modulus`.
  static absl::StatusOr<Polynomial> CreateFromProto(
      const proto::Polynomial& proto, int len, int log_modulus);

  // Default constructor.
  Polynomial() = default;

  // Copy constructor.
  Polynomial(const Polynomial& p) = default;
  Polynomial& operator=(const Polynomial& that) = default;

  // Move constructor.
  Polynomial(Polynomial&& p) = default;
  Polynomial& operator=(Polynomial&& that) = default;

  // Coordinate-wise addition.
  absl::StatusOr<Polynomial> Add(const Polynomial& that) const;

  // Coordinate-wise addition in place.
  absl::Status AddInPlace(const Polynomial& that);

  // Coordinate-wise negation.
  absl::StatusOr<Polynomial> Negate() const;

  // Rescales coefficients from current_modulus to new_modulus.
  // This function computes round(c * new_modulus / current_modulus) for each
  // coefficient c.
  absl::StatusOr<Polynomial> Rescale(int log_current_modulus,
                                     int log_new_modulus) const;

  // Right shifts the coefficients by `amount`.
  // amount must be in [0, 8 * sizeof(CoeffType)].
  absl::StatusOr<Polynomial> RightShift(int amount) const;

  // Keeps the lower `amount` bits of the coefficients.
  // amount must be in [0, 8 * sizeof(CoeffType)].
  absl::StatusOr<Polynomial> LowBits(int amount) const;

  // Performs gadget decomposition of the polynomial. Returns a vector of
  // `num_digits` polynomials, each representing a digit of length `log_digit`.
  absl::StatusOr<std::vector<Polynomial>> GadgetInv(int log_modulus,
                                                    int log_digit,
                                                    int num_digits) const;

  // Automorphism: Given a polynomial representing p(x), returns a
  // polynomial representing p(x^power). Power must be an odd non-negative
  // integer.
  absl::StatusOr<Polynomial> Automorph(int power) const;

  // Multiplication using Number Theoretic Transform.
  // Performs multiplication modulo (X^N + 1) where N is the length of the
  // coefficient vector, and the result's coefficients are modulo 2^32 or 2^64
  // depending on CoeffType.
  //
  // Requires a precomputed Context object configured for the same length N.
  absl::StatusOr<Polynomial> Mult(const Polynomial& that,
                                  const Context& ctx) const;

  // Multiplies using FFT with a reusable precomputed context to avoid
  // plan creation and buffer allocation overhead.
  absl::StatusOr<Polynomial> MultFft(const Polynomial& that, FftContext& ctx,
                                     int this_bits = 8 * sizeof(CoeffType),
                                     int that_bits = 8 * sizeof(CoeffType),
                                     int chunk_bits = 20) const;

  // Computes the inner product of two vectors of polynomials efficiently
  // by performing all accumulations in the FFT domain and doing exactly one
  // set of inverse FFTs at the end.
  static absl::StatusOr<Polynomial> InnerProductFft(
      const std::vector<Polynomial>& u, const std::vector<Polynomial>& v,
      FftContext& ctx,
      int u_bits = 8 * sizeof(CoeffType),
      int v_bits = 8 * sizeof(CoeffType),
      int chunk_bits = 20);

  // Converts the polynomial to NTT form.
  absl::StatusOr<NttPolynomial> ToNtt(const Context& ctx) const;

  // Converts from NTT form to coefficient form.
  static absl::StatusOr<Polynomial> FromNtt(const NttPolynomial& ntt_poly,
                                            const Context& ctx);

  // Accessor for coefficients.
  const std::vector<CoeffType>& Coeffs() const { return coeffs_; }

  int Len() const { return coeffs_.size(); }

  // Serializes the polynomial by packing the lowest `log_modulus` bits of each
  // coefficient into a sequence of bits represented as a byte string.
  absl::StatusOr<proto::Polynomial> ToProto(int log_modulus) const;

 private:
  explicit Polynomial(std::vector<CoeffType> coeffs)
      : coeffs_(std::move(coeffs)) {}

  std::vector<CoeffType> coeffs_;
};

extern template class Polynomial<uint32_t>;
extern template class Polynomial<uint64_t>;

}  // namespace v2
}  // namespace rlwe
}  // namespace private_membership

#endif  // PRIVACY_PRIVATE_MEMBERSHIP_RLWE_V2_CRYPTO_POLYNOMIAL_H_
