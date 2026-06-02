#ifndef PRIVACY_PRIVATE_MEMBERSHIP_OPRF_CLIENT_CPP_OPRF_CLIENT_H_
#define PRIVACY_PRIVATE_MEMBERSHIP_OPRF_CLIENT_CPP_OPRF_CLIENT_H_

#include <memory>
#include <string>
#include <vector>

#include "private_join_and_compute/crypto/context.h"
#include "private_join_and_compute/crypto/ec_commutative_cipher.h"
#include "private_membership/oprf/proto/oprf.pb.h"
#include "private_membership/oprf/shared/oprf_utils.h"
#include "absl/container/flat_hash_map.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"

namespace private_membership::oprf {

// This class generates encrypted requests and processes encrypted responses to
// execute the oblivious PRF (OPRF) protocol.
//
// Each execution of the OPRF protocol requires a new instance of the
// OprfClient.
//
// This class is not thread-safe.
class OprfClient {
 public:
  // Create client with randomly generated key.
  static absl::StatusOr<std::unique_ptr<OprfClient>> Create(
      const UseCase& use_case, const OprfParameters& params);

  // Create client with specified key.
  static absl::StatusOr<std::unique_ptr<OprfClient>> CreateFromKey(
      const UseCase& use_case, const OprfParameters& params,
      absl::string_view key);

  // Decrypt given encrypted value using given encrypted key.
  absl::StatusOr<std::string> DecryptValue(
      const EncryptedIdAndKey& encrypted_id_and_key,
      absl::string_view encrypted_value);

  // Create OPRF re-encrypt request given plaintext IDs and metadata.
  absl::StatusOr<OprfReencryptRequest> CreateRequestWithMetadata(
      int key_version, absl::Span<const PlaintextIdAndMetadata> ids);

  // Create OPRF re-encrypt request given plaintext IDs and no metadata.
  absl::StatusOr<OprfReencryptRequest> CreateRequest(
      int key_version, absl::Span<const std::string> ids);

  // Process response to get mapping from plaintext IDs to encrypted IDs and
  // keys.
  absl::StatusOr<absl::flat_hash_map<std::string, EncryptedIdAndKey>>
  RetrieveEncryptedIdsAndKeys(const OprfReencryptResponse& response);

 private:
  OprfClient(const UseCase& use_case, const OprfParameters& params_,
             std::unique_ptr<::private_join_and_compute::ECCommutativeCipher> cipher);

  const UseCase use_case_;
  const OprfParameters params_;
  const std::unique_ptr<::private_join_and_compute::ECCommutativeCipher> cipher_;
  ::private_join_and_compute::Context context_;
  absl::flat_hash_map<std::string, std::string> encrypted_to_plaintext_id_;
};

}  // namespace private_membership::oprf

#endif  // PRIVACY_PRIVATE_MEMBERSHIP_OPRF_CLIENT_CPP_OPRF_CLIENT_H_
