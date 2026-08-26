#ifndef PRIVACY_PRIVATE_MEMBERSHIP_RLWE_V2_CRYPTO_PIR_CLIENT_PIR_CLIENT_H_
#define PRIVACY_PRIVATE_MEMBERSHIP_RLWE_V2_CRYPTO_PIR_CLIENT_PIR_CLIENT_H_

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "private_membership/rlwe/v2/crypto/encryption.h"
#include "private_membership/rlwe/v2/crypto/pir/pir_params.h"
#include "private_membership/rlwe/v2/crypto/polynomial.h"
#include "absl/status/statusor.h"
#include "shell_encryption/prng/prng.h"

namespace private_membership {
namespace rlwe {
namespace v2 {

// A PIR client that generates PIR requests and processes PIR responses.
template <typename CoeffType>
class PirClient {
 public:
  // Creates and initializes a PirClient.
  static absl::StatusOr<std::unique_ptr<PirClient<CoeffType>>> Create(
      const PirParams<CoeffType>& params,
      std::unique_ptr<::rlwe::SecurePrng> a_component_prng,
      std::unique_ptr<::rlwe::SecurePrng> prng);

  // Generates a PIR request for the given database index.
  absl::StatusOr<PirRequest<CoeffType>> CreateRequest(int index);

  // Generates a PIR request encrypting an arbitrary vector of length a multiple
  // of ring degree d.
  // Note that this does not support second dimension PIR query.
  absl::StatusOr<PirRequest<CoeffType>> CreateRequestForVector(
      const std::vector<CoeffType>& vec);

  // Returns the internally generated secret key.
  // Exposed primarily for testing and debugging purposes.
  const Polynomial<CoeffType>& GetSecretKey() const { return secret_key_; }

  // Returns PIR parameters.
  const PirParams<CoeffType>& GetParams() const { return params_; }

  // Returns mutable pointer to the internally stored PRNG used to generate "a"
  // component of RLWE ciphertexts.
  ::rlwe::SecurePrng* GetAComponentPrng() const {
    return a_component_prng_.get();
  }

  // Returns mutable pointer to the internally stored PRNG used for other
  // randomness.
  ::rlwe::SecurePrng* GetPrng() const { return prng_.get(); }

  // Processes a PirResponse to extract the resulting message using the inner
  // secret key.
  absl::StatusOr<std::vector<CoeffType>> ProcessResponse(
      const PirResponse<CoeffType>& response) const;

  // Creates the part of PIR query corresponding to the second dimension of
  // database items.
  absl::StatusOr<std::vector<Polynomial<CoeffType>>> CreateSecondDimensionQuery(
      int index) const;

  // Creates the key used for packing..
  absl::StatusOr<std::vector<Polynomial<CoeffType>>> CreatePackingKey() const;

 private:
  // Creates the part of PIR query corresponding to the first dimension of
  // database items.
  absl::StatusOr<std::vector<Polynomial<CoeffType>>> CreateFirstDimensionQuery(
      int index) const;

  absl::StatusOr<std::vector<Polynomial<CoeffType>>>
  CreateFirstDimensionQueryForVector(const std::vector<CoeffType>& vec) const;

  PirClient(const PirParams<CoeffType>& params,
            std::unique_ptr<::rlwe::SecurePrng> a_component_prng,
            std::unique_ptr<::rlwe::SecurePrng> prng,
            Polynomial<CoeffType> secret_key)
      : params_(params),
        a_component_prng_(std::move(a_component_prng)),
        prng_(std::move(prng)),
        secret_key_(std::move(secret_key)) {}

  // Parameters for the PIR scheme.
  const PirParams<CoeffType> params_;

  // PRNG used to generate the "a" component of RLWE ciphertexts.
  const std::unique_ptr<::rlwe::SecurePrng> a_component_prng_;

  // PRNG used for other randomness, such as error generation.
  const std::unique_ptr<::rlwe::SecurePrng> prng_;

  // The secret key of used in the query generation..
  const Polynomial<CoeffType> secret_key_;
};

extern template class PirClient<uint32_t>;
extern template class PirClient<uint64_t>;

}  // namespace v2
}  // namespace rlwe
}  // namespace private_membership

#endif  // PRIVACY_PRIVATE_MEMBERSHIP_RLWE_V2_CRYPTO_PIR_CLIENT_PIR_CLIENT_H_
