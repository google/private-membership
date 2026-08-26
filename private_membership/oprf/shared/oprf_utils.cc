#include "private_membership/oprf/shared/oprf_utils.h"

#include <memory>
#include <string>
#include <utility>

#include "private_join_and_compute/crypto/big_num.h"
#include "private_join_and_compute/crypto/context.h"
#include "private_join_and_compute/crypto/ec_group.h"
#include "private_join_and_compute/crypto/ec_commutative_cipher.h"
#include "private_membership/rlwe/client/cpp/internal/aes_ctr_256_no_auth.h"
#include "private_membership/rlwe/client/cpp/internal/aes_ctr_256_with_fixed_iv.h"
#include "private_membership/rlwe/client/cpp/internal/aes_gcm_256.h"
#include "private_membership/rlwe/client/cpp/internal/cipher.h"
#include "private_membership/rlwe/client/cpp/internal/crypto_utils.h"
#include "private_membership/oprf/proto/oprf.pb.h"
#include "private_membership/oprf/shared/oprf_constants.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "private_join_and_compute/util/status_macros.h"
#include <openssl/nid.h>
#include <openssl/obj_mac.h>
#include "tink/deterministic_aead.h"
#include "tink/subtle/aes_siv_boringssl.h"
#include "tink/subtle/common_enums.h"
#include "tink/subtle/hkdf.h"
#include "tink/util/secret_data.h"

namespace private_membership::oprf {
namespace {

absl::StatusOr<std::unique_ptr<crypto::tink::DeterministicAead>>
GenerateAesSivCipher(absl::string_view encryption_key) {
  // Create AES-SIV cipher to encrypt value.
  crypto::tink::util::SecretData aes_key =
      crypto::tink::util::SecretDataFromStringView(encryption_key);
  absl::StatusOr<std::unique_ptr<crypto::tink::DeterministicAead>> aes_cipher =
      crypto::tink::subtle::AesSivBoringSsl::New(aes_key);
  return aes_cipher;
}

}  // namespace

Cipher::Scheme ToCustomCipherScheme(OprfEncryptionInternal encryption_type) {
  switch (encryption_type) {
    case OPRF_ENCRYPTION_AES_CTR_256_WITH_FIXED_IV:
      return Cipher::Scheme::kAesCtr256WithFixedIv;
    case OPRF_ENCRYPTION_AES_GCM_256:
      return Cipher::Scheme::kAesGcm256;
    case OPRF_ENCRYPTION_AES_CTR_256_NO_AUTH:
      return Cipher::Scheme::kAesCtr256NoAuth;
    case OPRF_ENCRYPTION_AES_SIV:
      return Cipher::Scheme::kAesSiv;
    default:
      return Cipher::Scheme::kEncryptionSchemeUndefined;
  }
}

absl::StatusOr<std::string> GenerateHashedId(absl::string_view encrypted_id,
                                             const OprfParameters& params,
                                             ::private_join_and_compute::Context& context) {
  switch (params.encryption_type) {
    case OPRF_ENCRYPTION_AES_SIV:
      return context.Sha256String(
          absl::StrCat(encrypted_id, params.hashed_id_salt));
    case OPRF_ENCRYPTION_AES_CTR_256_WITH_FIXED_IV:
    case OPRF_ENCRYPTION_AES_GCM_256:
    case OPRF_ENCRYPTION_AES_CTR_256_NO_AUTH:
      return HashEncryptedId(encrypted_id, &context);
    default:
      return absl::InvalidArgumentError("Invalid encryption type.");
  }
}

absl::StatusOr<std::string> GenerateEncryptionKey(
    absl::string_view encrypted_id, const OprfParameters& params) {
  switch (params.encryption_type) {
    case OPRF_ENCRYPTION_AES_SIV:
      return crypto::tink::subtle::Hkdf::ComputeHkdf(
          crypto::tink::subtle::HashType::SHA256, encrypted_id,
          params.hkdf_salt, params.hkdf_info, params.encryption_key_length);
    case OPRF_ENCRYPTION_AES_CTR_256_WITH_FIXED_IV:
    case OPRF_ENCRYPTION_AES_GCM_256:
    case OPRF_ENCRYPTION_AES_CTR_256_NO_AUTH: {
      ::private_join_and_compute::Context ctx;
      return GetValueEncryptionKey(encrypted_id, &ctx);
    }
    default:
      return absl::InvalidArgumentError("Invalid encryption type.");
  }
}

absl::StatusOr<std::string> EncryptValue(const absl::string_view encrypted_id,
                                         const absl::string_view value,
                                         const OprfParameters& params) {
  PJC_ASSIGN_OR_RETURN(std::string encryption_key,
                   GenerateEncryptionKey(encrypted_id, params));

  switch (params.encryption_type) {
    case OPRF_ENCRYPTION_AES_SIV: {
      // Create AES-GCM-SIV cipher to encrypt value.
      PJC_ASSIGN_OR_RETURN(
          std::unique_ptr<crypto::tink::DeterministicAead> aes_cipher,
          GenerateAesSivCipher(encryption_key));
      if (params.encrypted_value_byte_length > 0) {
        PJC_ASSIGN_OR_RETURN(const std::string padded_value,
                         Pad(value, params.encrypted_value_byte_length));
        return aes_cipher->EncryptDeterministically(padded_value, "");
      } else {
        return aes_cipher->EncryptDeterministically(value, "");
      }
    }
    case OPRF_ENCRYPTION_AES_CTR_256_WITH_FIXED_IV:
    case OPRF_ENCRYPTION_AES_GCM_256:
    case OPRF_ENCRYPTION_AES_CTR_256_NO_AUTH: {
      ::private_join_and_compute::Context ctx;
      if (params.encrypted_value_byte_length > 0) {
        return ::private_membership::EncryptValue(
            encrypted_id, value, params.encrypted_value_byte_length,
            ToCustomCipherScheme(params.encryption_type), &ctx);
      } else {
        return ::private_membership::EncryptValue(
            encrypted_id, value, value.length(),
            ToCustomCipherScheme(params.encryption_type), &ctx);
      }
    }
    default:
      return absl::InvalidArgumentError("Invalid encryption type specified.");
  }
}

absl::StatusOr<std::string> DecryptValueWithKey(
    const absl::string_view encryption_key,
    const absl::string_view encrypted_value, const OprfParameters& params) {
  switch (params.encryption_type) {
    case OPRF_ENCRYPTION_AES_SIV: {
      PJC_ASSIGN_OR_RETURN(
          std::unique_ptr<crypto::tink::DeterministicAead> aes_cipher,
          GenerateAesSivCipher(encryption_key));
      auto decrypted_value_or = aes_cipher->DecryptDeterministically(encrypted_value, "");
      if (!decrypted_value_or.ok()) {
        return absl::InternalError("Decryption failed");
      }
      const std::string decrypted_value = decrypted_value_or.ValueOrDie();

      if (params.encrypted_value_byte_length > 0) {
        return Unpad(decrypted_value);
      } else {
        return decrypted_value;
      }
    }
    case OPRF_ENCRYPTION_AES_CTR_256_WITH_FIXED_IV:
    case OPRF_ENCRYPTION_AES_GCM_256:
    case OPRF_ENCRYPTION_AES_CTR_256_NO_AUTH: {
      ::private_join_and_compute::Context ctx;
      return ::private_membership::DecryptValueWithKey(
          encryption_key, encrypted_value,
          ToCustomCipherScheme(params.encryption_type), &ctx);
    }
    default:
      return absl::InvalidArgumentError("Invalid encryption type specified.");
  }
}

}  // namespace private_membership::oprf
