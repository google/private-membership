#ifndef PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_AES_GCM_256_H_
#define PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_AES_GCM_256_H_

#include <cstddef>
#include <memory>
#include <string>
#include <utility>

#include "private_membership/rlwe/client/cpp/internal/cipher.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include <openssl/aead.h>
#include <openssl/base.h>

// Implements a wrapper around OpenSSL AES-GCM with 256 bit keys.
// This class is not thread-safe.
namespace private_membership {

class AesGcm256 : public Cipher {
 public:
  static const size_t kKeySize = 32;    // Key length = 256 bits
  static const size_t kNonceSize = 12;  // Nonce length = 96 bits

  // AesGcm256 cannot be copied or copy-assigned
  AesGcm256(const AesGcm256&) = delete;
  AesGcm256& operator=(const AesGcm256&) = delete;

  // AesGcm256 can be moved and move-assigned
  AesGcm256(AesGcm256&&) = default;
  AesGcm256& operator=(AesGcm256&&) = default;

  // Returns an instance of an AesGcm256 object, which can be used to seal or
  // open messages.
  // Returns INVALID_ARGUMENT if the key is not 256 bits.
  // Returns INTERNAL if the creation fails.
  static absl::StatusOr<std::unique_ptr<AesGcm256>> Create(
      absl::string_view key);

  // Encrypts a given message with AES-256, and tags the encryption using GCM.
  // Internally generates a random nonce for the encryption.
  // Returns INTERNAL if the underlying crypto fails.
  absl::StatusOr<std::string> Encrypt(
      absl::string_view plaintext) const override;

  // Decrypts a message, additionally checking that the GCM tag verifies
  // correctly.
  // Returns INVALID_ARGUMENT if ciphertext does not have a nonce of the right
  // size.
  // Returns INTERNAL if the underlying crypto fails. This may indicate that the
  // ciphertext or nonce has been modified.
  absl::StatusOr<std::string> Decrypt(
      absl::string_view ciphertext) const override;

 private:
  explicit AesGcm256(bssl::UniquePtr<EVP_AEAD_CTX> ctx)
      : ctx_(std::move(ctx)) {};

  bssl::UniquePtr<EVP_AEAD_CTX> ctx_;
};

}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_AES_GCM_256_H_
