// Copyright 2026 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_CRYPTO_UTILS_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_CRYPTO_UTILS_H_

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

// 4-byte magic identifier for authenticating encrypted real blocks: "ORAM"
constexpr uint32_t kBlockMagic = 0x4F52414D;

// Header size: magic (4B) + block_id (8B) + version (8B) + payload_len (4B) +
// checksum (4B) = 28B.
constexpr size_t kBlockHeaderSize = sizeof(uint32_t) + sizeof(uint64_t) +
                                    sizeof(uint64_t) + sizeof(uint32_t) +
                                    sizeof(uint32_t);

// Decrypted block payload and metadata.
struct DecryptedBlock {
  BlockId block_id;
  uint64_t version = 0;
  std::vector<uint8_t> payload;
};

// Generates cryptographically secure random bytes of specified length.
std::string GenerateRandomBytes(size_t length = 16);

// Generates a deterministic pseudo-random key derived from master_key and
// index.
std::string DeriveSubKey(absl::string_view master_key, uint64_t index);

// Generates a deterministic pseudo-random key derived from master_key and
// string tag.
std::string DeriveSubKey(absl::string_view master_key, absl::string_view tag);

// Generates a public block label for a real block with block_id.
std::string GenerateBlockLabel(absl::string_view key_l, BlockId block_id);

// Generates a public block label for a dummy block with dummy_index.
std::string GenerateDummyLabel(absl::string_view key_l, uint64_t dummy_index);

// Generates a pseudo-random keystream of the given length.
std::vector<uint8_t> GeneratePseudoRandomStream(absl::string_view key,
                                                absl::string_view label,
                                                size_t length);

// Generates pseudo-random dummy block contents of size `block_size`.
std::vector<uint8_t> GenerateDummyBlock(absl::string_view key_d,
                                        absl::string_view label,
                                        size_t block_size);

// Computes the total ciphertext size for a block with user payload
// `payload_size`.
constexpr size_t CiphertextBlockSize(size_t payload_size) {
  return payload_size + kBlockHeaderSize;
}

// Encrypts a real block with `block_id`, `payload`, and optional `version`.
std::vector<uint8_t> EncryptBlock(absl::string_view key_e,
                                  absl::string_view label, BlockId block_id,
                                  absl::Span<const uint8_t> payload,
                                  size_t payload_size, uint64_t version = 0);

// Decrypts a ciphertext block, validates magic & checksum, and extracts the
// block_id and payload.
absl::StatusOr<DecryptedBlock> DecryptBlock(
    absl::string_view key_e, absl::string_view label,
    absl::Span<const uint8_t> ciphertext, size_t payload_size);

// Generates a uniform random permutation of {0, 1, ..., n-1}.
std::vector<size_t> GenerateRandomPermutation(size_t n, uint64_t seed = 0);

// Helper functions to serialize and deserialize OKVS tables into byte strings.
std::string SerializeOkvsTable(const std::vector<std::vector<uint8_t>>& table);

std::vector<std::vector<uint8_t>> DeserializeOkvsTable(absl::string_view bytes);

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_CRYPTO_UTILS_H_
