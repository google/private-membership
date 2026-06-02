#include "private_membership/oprf/client/cpp/oprf_client.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "private_join_and_compute/crypto/ec_commutative_cipher.h"
#include "private_membership/oprf/proto/oprf.pb.h"
#include "private_membership/oprf/shared/oprf_utils.h"
#include "absl/container/flat_hash_map.h"
#include "absl/memory/memory.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "private_join_and_compute/util/status_macros.h"

namespace private_membership::oprf {

OprfClient::OprfClient(const UseCase& use_case, const OprfParameters& params,
                       std::unique_ptr<::private_join_and_compute::ECCommutativeCipher> cipher)
    : use_case_(use_case), params_(params), cipher_(std::move(cipher)) {}

absl::StatusOr<std::unique_ptr<OprfClient>> OprfClient::Create(
    const UseCase& use_case, const OprfParameters& params) {
  PJC_ASSIGN_OR_RETURN(auto cipher,
                   ::private_join_and_compute::ECCommutativeCipher::CreateWithNewKey(
                       params.curve_id, params.hash_type));
  return absl::WrapUnique<OprfClient>(
      new OprfClient(use_case, params, std::move(cipher)));
}

absl::StatusOr<std::unique_ptr<OprfClient>> OprfClient::CreateFromKey(
    const UseCase& use_case, const OprfParameters& params,
    absl::string_view key) {
  PJC_ASSIGN_OR_RETURN(auto cipher, ::private_join_and_compute::ECCommutativeCipher::CreateFromKey(
                                    params.curve_id, key, params.hash_type));
  return absl::WrapUnique<OprfClient>(
      new OprfClient(use_case, params, std::move(cipher)));
}

absl::StatusOr<OprfReencryptRequest> OprfClient::CreateRequest(
    int key_version, absl::Span<const std::string> ids) {
  std::vector<PlaintextIdAndMetadata> ids_and_metadata(ids.size());
  for (int i = 0; i < ids.size(); ++i) {
    ids_and_metadata[i].set_plaintext_id(ids[i]);
  }
  return CreateRequestWithMetadata(key_version, ids_and_metadata);
}

absl::StatusOr<OprfReencryptRequest> OprfClient::CreateRequestWithMetadata(
    int key_version, absl::Span<const PlaintextIdAndMetadata> ids) {
  if (!encrypted_to_plaintext_id_.empty()) {
    return absl::FailedPreconditionError("Cannot create a second request.");
  } else if (ids.empty()) {
    return absl::InvalidArgumentError("Cannot create an empty request.");
  } else if (key_version <= 0) {
    return absl::InvalidArgumentError("Key versions must be positive.");
  }

  OprfReencryptRequest request;
  *request.mutable_use_case() = use_case_;
  request.set_key_version(key_version);

  for (const PlaintextIdAndMetadata& id : ids) {
    PJC_ASSIGN_OR_RETURN(std::string client_encrypted_id,
                     cipher_->Encrypt(id.plaintext_id()));
    auto* encrypted_id = request.add_encrypted_ids();
    encrypted_id->set_encrypted_id(client_encrypted_id);
    encrypted_id->set_public_metadata(id.public_metadata());
    encrypted_to_plaintext_id_[client_encrypted_id] = id.plaintext_id();
  }

  return request;
}

absl::StatusOr<absl::flat_hash_map<std::string, EncryptedIdAndKey>>
OprfClient::RetrieveEncryptedIdsAndKeys(const OprfReencryptResponse& response) {
  if (encrypted_to_plaintext_id_.empty()) {
    return absl::FailedPreconditionError("Must create a request first.");
  }
  absl::flat_hash_map<std::string, EncryptedIdAndKey> result;
  for (const OprfReencryptResponse::DoublyEncryptedId& double_encrypted_id :
       response.doubly_encrypted_ids()) {
    auto it = encrypted_to_plaintext_id_.find(
        double_encrypted_id.queried_encrypted_id());
    if (it == encrypted_to_plaintext_id_.end()) {
      return absl::InvalidArgumentError(
          "Response is malformed containing a non-queried ID.");
    }

    PJC_ASSIGN_OR_RETURN(
        std::string decrypted,
        cipher_->Decrypt(double_encrypted_id.doubly_encrypted_id()));

    EncryptedIdAndKey encrypted_id_and_key;
    PJC_ASSIGN_OR_RETURN(*encrypted_id_and_key.mutable_encrypted_id(),
                     GenerateHashedId(decrypted, params_, context_));
    PJC_ASSIGN_OR_RETURN(*encrypted_id_and_key.mutable_encryption_key(),
                     GenerateEncryptionKey(decrypted, params_));

    result[it->second] = std::move(encrypted_id_and_key);
  }

  return result;
}

absl::StatusOr<std::string> OprfClient::DecryptValue(
    const EncryptedIdAndKey& encrypted_id_and_key,
    absl::string_view encrypted_value) {
  return DecryptValueWithKey(encrypted_id_and_key.encryption_key(),
                             encrypted_value, params_);
}

}  // namespace private_membership::oprf
