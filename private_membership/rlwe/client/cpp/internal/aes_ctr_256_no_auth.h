#ifndef PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_AES_CTR_256_NO_AUTH_H_
#define PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_AES_CTR_256_NO_AUTH_H_

#include <memory>
#include <string>

#include "private_membership/rlwe/client/cpp/internal/cipher.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"

namespace private_membership {

// Wrapper around OpenSSL AES-CTR with 256 bit keys with no authentication. As a
// result, the decryption will not fail even if the ciphertext is invalid.
//
// This cipher is intended to be used for Targeted User Abuse Prevention
// modification of blinders based PSM protocol. We use unauthenticated
// encryption to hide the membership result from the client See the design doc
// go/psm-preventing-targeted-user-attacks for more detail.
//
// DO NOT USE this cipher for any other use case.
//
// We have our own implementation because we would like to avoid introducing
// dependencies other than OpenSSL to keep the size of the binary small as
// requested by our customer.
class AesCtr256NoAuth : public Cipher {
 public:
  // The size of the key in bytes.
  static constexpr int kKeySize = 32;
  // The size of the block in bytes.
  static constexpr int kBlockSize = 16;
  // The size of the IV in bytes.
  static constexpr int kIVSize = 16;

  // Create a cipher with the input key.
  //
  // Returns INVALID_ARGUMENT if key is an invalid size.
  static absl::StatusOr<std::unique_ptr<AesCtr256NoAuth>> Create(
      absl::string_view key);

  // Encrypts the plaintext.
  //
  // Returns INTERNAL on internal cryptographic errors.
  absl::StatusOr<std::string> Encrypt(
      absl::string_view plaintext) const override;

  // Decrypts the ciphertext.
  //
  // This function does not fail even if the ciphertext is invalid.
  absl::StatusOr<std::string> Decrypt(
      absl::string_view ciphertext) const override;

 private:
  explicit AesCtr256NoAuth(absl::string_view key);

  const std::string key_;
};

}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_AES_CTR_256_NO_AUTH_H_
