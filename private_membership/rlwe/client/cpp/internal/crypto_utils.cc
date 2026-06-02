// Copyright 2020 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "private_membership/rlwe/client/cpp/internal/crypto_utils.h"

#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <utility>

#include "private_join_and_compute/crypto/context.h"
#include "private_join_and_compute/crypto/ec_commutative_cipher.h"
#include "private_membership/rlwe/client/cpp/internal/aes_ctr_256_no_auth.h"
#include "private_membership/rlwe/client/cpp/internal/aes_ctr_256_with_fixed_iv.h"
#include "private_membership/rlwe/client/cpp/internal/aes_gcm_256.h"
#include "private_membership/rlwe/client/cpp/internal/cipher.h"
#include "private_membership/rlwe/client/cpp/internal/id_utils.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"

namespace private_membership {
namespace {

std::string LengthToBytes(uint32_t length) {
  uint8_t ret[4];
  ret[0] = length & 0xFF;
  ret[1] = (length >> 8) & 0xFF;
  ret[2] = (length >> 16) & 0xFF;
  ret[3] = (length >> 24) & 0xFF;
  return std::string(ret, ret + 4);
}

absl::StatusOr<uint32_t> BytesToLength(absl::string_view bytes) {
  if (bytes.size() != 4) {
    return absl::InvalidArgumentError("Invalid byte size of length encoding.");
  }
  uint32_t ret = 0;
  ret += static_cast<uint8_t>(bytes[0]);
  ret += (static_cast<uint8_t>(bytes[1]) << 8);
  ret += (static_cast<uint8_t>(bytes[2]) << 16);
  ret += (static_cast<uint8_t>(bytes[3]) << 24);
  return ret;
}

absl::StatusOr<std::unique_ptr<Cipher>> CreateCipher(
    const Cipher::Scheme encryption_scheme, const absl::string_view key) {
  switch (encryption_scheme) {
    case Cipher::Scheme::kAesCtr256WithFixedIv:
      return AesCtr256WithFixedIV::Create(key);
    case Cipher::Scheme::kAesCtr256NoAuth:
      return AesCtr256NoAuth::Create(key);
    case Cipher::Scheme::kAesGcm256:
      return AesGcm256::Create(key);
    default:
      return absl::InvalidArgumentError(
          absl::StrCat("Unsupported encryption scheme: ", encryption_scheme));
  }
}

}  // namespace

absl::StatusOr<std::string> Pad(absl::string_view bytes, int max_byte_length) {
  if (max_byte_length < bytes.length()) {
    return absl::InvalidArgumentError(
        "max_byte_length smaller than the input bytes length.");
  }
  // Prepend byte length of value into the first four bytes. This limits the
  // maximum length of the value to at most 2^32 bytes or ~4.3 gigabytes.
  std::string length_in_byte = LengthToBytes(bytes.length());
  std::string padded_bytes =
      std::string(bytes).append(max_byte_length - bytes.length(), 0);
  return absl::StrCat(length_in_byte, padded_bytes);
}

absl::StatusOr<std::string> Unpad(absl::string_view bytes) {
  if (bytes.length() < 4) {
    return absl::InvalidArgumentError("Invalid bytes does not encode length.");
  }
  absl::StatusOr<uint32_t> bytes_length = BytesToLength(bytes.substr(0, 4));
  if (!bytes_length.ok()) {
    return bytes_length.status();
  }
  if (*bytes_length + 4 > bytes.length()) {
    return absl::InvalidArgumentError("Incorrect bytes length.");
  }
  return std::string(bytes.substr(4, *bytes_length));
}

std::string HashEncryptedId(absl::string_view encrypted_id,
                            ::private_join_and_compute::Context* ctx) {
  // Salt used in SHA256 hash function for creating sensitive matching
  // identifier. Randomly generated salt forcing adversaries to re-compute
  // rainbow tables.
  static constexpr unsigned char kSensitiveIdHashSalt[] = {
      0x3C, 0xD1, 0xF3, 0x69, 0x2B, 0x57, 0x40, 0xEA, 0xD8, 0xE4, 0xF4,
      0x4A, 0xB2, 0x5F, 0x7B, 0xAD, 0xC8, 0x10, 0xAA, 0x3D, 0x4C, 0x6E,
      0xCA, 0x57, 0x78, 0x5C, 0x5A, 0xED, 0x06, 0x81, 0x14, 0x7C};
  return ctx->Sha256String(absl::StrCat(
      std::string(reinterpret_cast<const char*>(kSensitiveIdHashSalt), 32),
      encrypted_id));
}

absl::StatusOr<std::string> GetValueEncryptionKey(
    absl::string_view encrypted_id, ::private_join_and_compute::Context* ctx) {
  // Salt user in SHA256 hash function to generate value encryption key.
  // Randomly generated salt forcing adversaries to re-compute rainbow tables.
  static constexpr unsigned char kValueEncryptionKeyHashSalt[] = {
      0x89, 0xC3, 0x67, 0xA7, 0x8A, 0x68, 0x2B, 0xE8, 0xC6, 0xB2, 0x22,
      0xB7, 0xE0, 0xB7, 0x4A, 0x37, 0x63, 0x8F, 0x10, 0x79, 0x98, 0x91,
      0x31, 0x94, 0x44, 0x03, 0xB6, 0x76, 0x8F, 0x70, 0xEB, 0xBF};
  // AES key length in bytes.
  static constexpr int kRequiredAesKeyLength = 32;
  return PadOrTruncate(
      ctx->Sha256String(absl::StrCat(
          std::string(
              reinterpret_cast<const char*>(kValueEncryptionKeyHashSalt), 32),
          encrypted_id)),
      kRequiredAesKeyLength);
}

absl::StatusOr<std::string> EncryptValue(const absl::string_view encrypted_id,
                                         const absl::string_view value,
                                         const uint32_t max_value_byte_length,
                                         const Cipher::Scheme encryption_scheme,
                                         ::private_join_and_compute::Context* ctx) {
  if (value.length() > max_value_byte_length) {
    return absl::InvalidArgumentError(absl::StrCat(
        "Length of value ", value.length(), " larger than maximum value byte ",
        "length ", max_value_byte_length, "."));
  }
  auto value_encryption_key = GetValueEncryptionKey(encrypted_id, ctx);
  if (!value_encryption_key.ok()) {
    return value_encryption_key.status();
  }

  absl::StatusOr<std::unique_ptr<Cipher>> cipher =
      CreateCipher(encryption_scheme, *value_encryption_key);
  if (!cipher.ok()) {
    return cipher.status();
  }

  absl::StatusOr<std::string> plaintext;
  // Skip padding for unauthenticated encryption because the for failed
  // decryption, the output will not be in the correct format.
  if (encryption_scheme == Cipher::Scheme::kAesCtr256NoAuth) {
    plaintext = std::string(value);
  } else {
    plaintext = Pad(value, max_value_byte_length);
  }
  if (!plaintext.ok()) {
    return plaintext.status();
  }

  return (*cipher)->Encrypt(*std::move(plaintext));
}

absl::StatusOr<std::string> DecryptValue(
    const absl::string_view encrypted_id,
    const absl::string_view encrypted_value,
    const Cipher::Scheme encryption_scheme, ::private_join_and_compute::Context* ctx) {
  auto value_encryption_key = GetValueEncryptionKey(encrypted_id, ctx);
  if (!value_encryption_key.ok()) {
    return value_encryption_key.status();
  }
  return DecryptValueWithKey(*value_encryption_key, encrypted_value,
                             encryption_scheme, ctx);
}

absl::StatusOr<std::string> DecryptValueWithKey(
    const absl::string_view encryption_key,
    const absl::string_view encrypted_value,
    const Cipher::Scheme encryption_scheme, ::private_join_and_compute::Context* ctx) {
  absl::StatusOr<std::unique_ptr<Cipher>> cipher =
      CreateCipher(encryption_scheme, encryption_key);
  if (!cipher.ok()) {
    return cipher.status();
  }
  auto plaintext = (*cipher)->Decrypt(encrypted_value);
  if (!plaintext.ok()) {
    return plaintext.status();
  }
  // Skip unpadding for unauthenticated encryption because the for failed
  // decryption, the output will not be in the correct format.
  if (encryption_scheme == Cipher::Scheme::kAesCtr256NoAuth) {
    return plaintext;
  }
  return Unpad(*std::move(plaintext));
}

absl::StatusOr<std::string> Truncate(std::string_view in, int bit_length) {
  if (bit_length < 0 || bit_length > in.size() * 8) {
    return absl::InvalidArgumentError("Truncation bit length out of bounds.");
  }
  if (bit_length == 0) {
    return std::string();
  }
  int ceiled_byte_length = (bit_length - 1) / 8 + 1;
  std::string res(in.begin(), in.begin() + ceiled_byte_length);
  if (int remainder_bit_count = bit_length % 8; remainder_bit_count > 0) {
    int mask = (1 << remainder_bit_count) - 1;
    res[res.size() - 1] &= (mask << (8 - remainder_bit_count));
  }
  return res;
}
}  // namespace private_membership
