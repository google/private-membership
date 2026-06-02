#ifndef PRIVACY_PRIVATE_MEMBERSHIP_OPRF_SHARED_OPRF_UTILS_H_
#define PRIVACY_PRIVATE_MEMBERSHIP_OPRF_SHARED_OPRF_UTILS_H_

#include <memory>
#include <string>

#include "private_join_and_compute/crypto/context.h"
#include "private_join_and_compute/crypto/ec_commutative_cipher.h"
#include "private_membership/rlwe/client/cpp/internal/cipher.h"
#include "private_membership/oprf/proto/oprf.pb.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"

namespace private_membership::oprf {

enum OprfEncryptionInternal {
  OPRF_ENCRYPTION_AES_SIV = 0,
  OPRF_ENCRYPTION_AES_CTR_256_WITH_FIXED_IV = 1,
  OPRF_ENCRYPTION_AES_CTR_256_NO_AUTH = 2,
  OPRF_ENCRYPTION_AES_GCM_256 = 3,
};

// Parameters for the OPRF protocol.
struct OprfParameters {
  int curve_id;
  ::private_join_and_compute::ECCommutativeCipher::HashType hash_type;
  OprfEncryptionInternal encryption_type;
  std::string hashed_id_salt;
  std::string hkdf_salt;
  std::string hkdf_info;
  int encryption_key_length;
  int encrypted_value_byte_length;
};

// Converts OprfEncryptionInternal to Cipher::Scheme.
// Exposed for test use only.
Cipher::Scheme ToCustomCipherScheme(OprfEncryptionInternal encryption_type);

// Generate an ID from an encrypted ID using parameters.
absl::StatusOr<std::string> GenerateHashedId(absl::string_view encrypted_id,
                                             const OprfParameters& params,
                                             ::private_join_and_compute::Context& context);

// Generate encryption key using parameters.
absl::StatusOr<std::string> GenerateEncryptionKey(
    absl::string_view encrypted_id, const OprfParameters& params);

// Encrypt a value using an encrypted ID using parameters.
absl::StatusOr<std::string> EncryptValue(absl::string_view encrypted_id,
                                         absl::string_view value,
                                         const OprfParameters& params);

// Decrypt a value using an encryption key.
absl::StatusOr<std::string> DecryptValueWithKey(
    absl::string_view encryption_key, absl::string_view encrypted_value,
    const OprfParameters& params);

}  // namespace private_membership::oprf

#endif  // PRIVACY_PRIVATE_MEMBERSHIP_OPRF_SHARED_OPRF_UTILS_H_
