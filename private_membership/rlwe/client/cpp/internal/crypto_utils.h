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

#ifndef PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_CRYPTO_UTILS_H_
#define PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_CRYPTO_UTILS_H_

#include <cstdint>
#include <string>
#include <string_view>

#include "private_join_and_compute/crypto/context.h"
#include "private_join_and_compute/crypto/ec_commutative_cipher.h"
#include "private_membership/rlwe/client/cpp/internal/cipher.h"
#include "private_membership/rlwe/client/cpp/open_source/chromium/private_membership_export.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"

namespace private_membership {


// Pad to `max_byte_length` bytes.
//
// The returned bytes will be max_byte_length + 4 bytes long, where the first
// four bytes represents the length of the original bytes.
PRIVATE_MEMBERSHIP_EXPORT absl::StatusOr<std::string> Pad(
    absl::string_view bytes, int max_byte_length);

// Unpad bytes padded with the `Pad` function.
PRIVATE_MEMBERSHIP_EXPORT absl::StatusOr<std::string> Unpad(
    absl::string_view bytes);

// Hash encrypted ID.
PRIVATE_MEMBERSHIP_EXPORT std::string HashEncryptedId(
    absl::string_view encrypted_id, ::private_join_and_compute::Context* ctx);

// Compute the value encryption key.
PRIVATE_MEMBERSHIP_EXPORT absl::StatusOr<std::string> GetValueEncryptionKey(
    absl::string_view encrypted_id, ::private_join_and_compute::Context* ctx);

// Encrypts the value, first padding it to. The result will be
// max_value_byte_length+4 bytes long, the first 4 bytes being an encryption of
// the length of the original value.
//
// Returns an error if the value is longer than max_value_byte_length.
PRIVATE_MEMBERSHIP_EXPORT absl::StatusOr<std::string> EncryptValue(
    absl::string_view encrypted_id, absl::string_view value,
    uint32_t max_value_byte_length, Cipher::Scheme encryption_scheme,
    ::private_join_and_compute::Context* ctx);

// Decrypt the encrypted_value using a key derived from the encrypted_id.
PRIVATE_MEMBERSHIP_EXPORT absl::StatusOr<std::string> DecryptValue(
    absl::string_view encrypted_id, absl::string_view encrypted_value,
    Cipher::Scheme encryption_scheme, ::private_join_and_compute::Context* ctx);

// Decrypt the encrypted_value using the supplied key.
PRIVATE_MEMBERSHIP_EXPORT absl::StatusOr<std::string> DecryptValueWithKey(
    absl::string_view encryption_key, absl::string_view encrypted_value,
    Cipher::Scheme encryption_scheme, ::private_join_and_compute::Context* ctx);

// Truncates the given bytes by returning the first `bit_length` leftmost bits.
//
// The input bytes `in` is assumed to be in big endian representation. The size
// of the result will be rounded up to the nearest byte. Thus, the last byte of
// the result may contain some extra bits, which will all be set to 0.
//
// For example, given 10111011 11111011 in big endian with `bit_length` == 12,
// the result will be 10111011 11110000.
//
// Returns InvalidArgumentError if the truncation bit length is larger than the
// bit length of the input string.
PRIVATE_MEMBERSHIP_EXPORT absl::StatusOr<std::string> Truncate(
    std::string_view in, int bit_length);

}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_CRYPTO_UTILS_H_
