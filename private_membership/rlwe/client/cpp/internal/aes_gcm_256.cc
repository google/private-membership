#include "private_membership/rlwe/client/cpp/internal/aes_gcm_256.h"

#include <cstddef>
#include <cstring>
#include <memory>
#include <string>
#include <utility>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include <openssl/aead.h>
#include <openssl/base.h>
#include <openssl/cipher.h>
#include <openssl/err.h>
#include <openssl/rand.h>

namespace private_membership {

namespace {

std::string OpenSSLErrorString() {
  char buf[256];
  ERR_error_string_n(ERR_get_error(), buf, sizeof(buf));
  return buf;
}

}  // namespace

absl::StatusOr<std::unique_ptr<AesGcm256>> AesGcm256::Create(
    const absl::string_view key) {
  if (key.size() != kKeySize) {
    return absl::InvalidArgumentError(absl::StrCat(
        "AesGcm256::Create : Key size must be ", kKeySize, " bytes."));
  }

  bssl::UniquePtr<EVP_AEAD_CTX> ctx(EVP_AEAD_CTX_new(
      EVP_aead_aes_256_gcm(), reinterpret_cast<const uint8_t*>(key.data()),
      kKeySize, /*tag_len=*/0));
  if (!ctx) {
    return absl::InternalError(
        absl::StrCat("EVP_AEAD_CTX_new failed:", OpenSSLErrorString()));
  }

  std::unique_ptr<AesGcm256> cipher(new AesGcm256(std::move(ctx)));
  return cipher;
}

absl::StatusOr<std::string> AesGcm256::Encrypt(
    absl::string_view plaintext) const {
  const size_t max_output_length =
      plaintext.length() + EVP_AEAD_max_overhead(ctx_->aead);
  // The ciphertext is a concatenation of the nonce and the encrypted plaintext.
  std::string ciphertext(kNonceSize + max_output_length, '\0');

  // Generate a random nonce.
  if (1 !=
      RAND_bytes(reinterpret_cast<uint8_t*>(ciphertext.data()), kNonceSize)) {
    return absl::InternalError(
        absl::StrCat("RAND_bytes failed:", OpenSSLErrorString()));
  }

  // Encrypt plaintext
  size_t output_length = 0;
  if (1 != EVP_AEAD_CTX_seal(
               ctx_.get(),
               reinterpret_cast<uint8_t*>(ciphertext.data() + kNonceSize),
               &output_length, max_output_length,
               reinterpret_cast<uint8_t*>(ciphertext.data()), kNonceSize,
               reinterpret_cast<const uint8_t*>(plaintext.data()),
               plaintext.length(),
               /*ad=*/nullptr, /*ad_len=*/0)) {
    return absl::InternalError(
        absl::StrCat("EVP_AEAD_CTX_seal failed:", OpenSSLErrorString()));
  }

  return ciphertext;
}

absl::StatusOr<std::string> AesGcm256::Decrypt(
    absl::string_view ciphertext) const {
  if (ciphertext.length() < kNonceSize) {
    return absl::InvalidArgumentError("Ciphertext is too short.");
  }

  // Split the nonce and the encrypted plaintext.
  const absl::string_view nonce = ciphertext.substr(0, kNonceSize);
  ciphertext = ciphertext.substr(kNonceSize);

  std::string plaintext(ciphertext.length(), '\0');
  size_t output_length;
  // Decrypt ciphertext.
  if (1 != EVP_AEAD_CTX_open(
               ctx_.get(), reinterpret_cast<uint8_t*>(plaintext.data()),
               &output_length, plaintext.length(),
               reinterpret_cast<const uint8_t*>(nonce.data()), kNonceSize,
               reinterpret_cast<const uint8_t*>(ciphertext.data()),
               ciphertext.length(), nullptr, 0)) {
    return absl::InternalError(
        absl::StrCat("EVP_AEAD_CTX_open failed:", OpenSSLErrorString()));
  }

  // Finalize the plaintext.
  plaintext.resize(output_length);
  return plaintext;
}

}  // namespace private_membership
