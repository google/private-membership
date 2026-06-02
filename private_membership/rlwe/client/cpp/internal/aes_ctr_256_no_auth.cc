#include "private_membership/rlwe/client/cpp/internal/aes_ctr_256_no_auth.h"

#include <cstdint>
#include <memory>
#include <string>

#include "absl/memory/memory.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include <openssl/base.h>
#include <openssl/cipher.h>
#include <openssl/rand.h>

namespace private_membership {

namespace {

absl::string_view EnsureNotNull(absl::string_view str) {
  if (str.empty() && str.data() == nullptr) {
    return absl::string_view("");
  }
  return str;
}

}  // namespace

absl::StatusOr<std::unique_ptr<AesCtr256NoAuth>> AesCtr256NoAuth::Create(
    absl::string_view key) {
  if (key.size() != kKeySize) {
    return absl::InvalidArgumentError("Key size is invalid.");
  }
  return absl::WrapUnique<AesCtr256NoAuth>(new AesCtr256NoAuth(key));
}

AesCtr256NoAuth::AesCtr256NoAuth(absl::string_view key) : key_(key) {}

absl::StatusOr<std::string> AesCtr256NoAuth::Encrypt(
    absl::string_view plaintext) const {
  plaintext = EnsureNotNull(plaintext);

  // Setup a context.
  bssl::UniquePtr<EVP_CIPHER_CTX> ctx(EVP_CIPHER_CTX_new());
  if (ctx.get() == nullptr) {
    return absl::InternalError("Could not initialize EVP_CIPHER_CTX.");
  }

  // Generate a random IV.
  std::string iv(kIVSize, '\0');
  if (RAND_bytes(reinterpret_cast<uint8_t*>(iv.data()), kIVSize) != 1) {
    return absl::InternalError("Could not generate an IV.");
  }

  // Initialize the encryption.
  if (EVP_EncryptInit_ex(ctx.get(), EVP_aes_256_ctr(),
                         /*impl=*/nullptr,
                         reinterpret_cast<const uint8_t*>(key_.data()),
                         reinterpret_cast<const uint8_t*>(iv.data())) != 1) {
    return absl::InternalError("Could not initialize encryption.");
  }

  // Encrypt.
  std::string ciphertext(iv);
  ciphertext.resize(kIVSize + plaintext.size(), '\0');
  int len = 0;
  if (EVP_EncryptUpdate(
          ctx.get(), reinterpret_cast<uint8_t*>(ciphertext.data() + kIVSize),
          &len, reinterpret_cast<const uint8_t*>(plaintext.data()),
          plaintext.size()) != 1) {
    return absl::InternalError("Encryption failed.");
  }
  ciphertext.resize(kIVSize + len);

  // Finalize the encryption.
  if (len != plaintext.size()) {
    return absl::InternalError("Ciphertext is incorrect size.");
  }
  return ciphertext;
}

absl::StatusOr<std::string> AesCtr256NoAuth::Decrypt(
    absl::string_view ciphertext) const {
  // Ensure the ciphertext is valid.
  if (ciphertext.data() == nullptr) {
    return absl::InternalError("Invalid null ciphertext.");
  }
  if (ciphertext.size() < kIVSize) {
    return absl::InternalError("Ciphertext is too short.");
  }

  // Setup a context.
  bssl::UniquePtr<EVP_CIPHER_CTX> ctx(EVP_CIPHER_CTX_new());
  if (ctx.get() == nullptr) {
    return absl::InternalError("Could not initialize EVP_CIPHER_CTX.");
  }

  // Split the IV and the ciphertext.
  const absl::string_view iv = ciphertext.substr(0, kIVSize);
  ciphertext = ciphertext.substr(kIVSize);

  // Initialize the decryption.
  if (EVP_DecryptInit_ex(ctx.get(), EVP_aes_256_ctr(),
                         /*impl=*/nullptr,
                         reinterpret_cast<const uint8_t*>(key_.data()),
                         reinterpret_cast<const uint8_t*>(iv.data())) != 1) {
    return absl::InternalError("Could not initialize decryption.");
  }

  // Decrypt.
  std::string plaintext(ciphertext.size(), '\0');
  int len = 0;
  if (EVP_DecryptUpdate(ctx.get(), reinterpret_cast<uint8_t*>(plaintext.data()),
                        &len,
                        reinterpret_cast<const uint8_t*>(ciphertext.data()),
                        ciphertext.size()) != 1) {
    return absl::InternalError("Decryption failed.");
  }
  plaintext.resize(len);

  // Finalize the decryption.
  if (len != ciphertext.size()) {
    return absl::InternalError("Plaintext is incorrect size.");
  }
  return plaintext;
}

}  // namespace private_membership
