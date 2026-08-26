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

#include "research/oram_okvs/oram/crypto_utils.h"

#include <openssl/rand.h>
#include <openssl/sha.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "absl/hash/hash.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

namespace {

uint32_t ComputeChecksum(BlockId block_id, uint64_t version,
                         absl::Span<const uint8_t> payload) {
  size_t h = absl::Hash<std::tuple<BlockId, uint64_t, absl::string_view>>{}(
      {block_id, version,
       absl::string_view(reinterpret_cast<const char*>(payload.data()),
                         payload.size())});
  return static_cast<uint32_t>(h ^ (h >> 32));
}

}  // namespace

std::string GenerateRandomBytes(size_t length) {
  std::string bytes(length, '\0');
  RAND_bytes(reinterpret_cast<uint8_t*>(bytes.data()),
             static_cast<int>(length));
  return bytes;
}

std::string DeriveSubKey(absl::string_view master_key, uint64_t index) {
  uint8_t digest[SHA256_DIGEST_LENGTH];
  SHA256_CTX ctx;
  SHA256_Init(&ctx);
  SHA256_Update(&ctx, master_key.data(), master_key.size());
  SHA256_Update(&ctx, &index, sizeof(index));
  SHA256_Final(digest, &ctx);
  return std::string(reinterpret_cast<const char*>(digest), 16);
}

std::string DeriveSubKey(absl::string_view master_key, absl::string_view tag) {
  uint8_t digest[SHA256_DIGEST_LENGTH];
  SHA256_CTX ctx;
  SHA256_Init(&ctx);
  SHA256_Update(&ctx, master_key.data(), master_key.size());
  SHA256_Update(&ctx, tag.data(), tag.size());
  SHA256_Final(digest, &ctx);
  return std::string(reinterpret_cast<const char*>(digest), 16);
}

std::string GenerateBlockLabel(absl::string_view key_l, BlockId block_id) {
  uint8_t digest[SHA256_DIGEST_LENGTH];
  SHA256_CTX ctx;
  SHA256_Init(&ctx);
  SHA256_Update(&ctx, key_l.data(), key_l.size());
  uint8_t tag = 0x01;  // Tag for real blocks
  SHA256_Update(&ctx, &tag, sizeof(tag));
  SHA256_Update(&ctx, &block_id, sizeof(block_id));
  SHA256_Final(digest, &ctx);
  return std::string(reinterpret_cast<const char*>(digest), 16);
}

std::string GenerateDummyLabel(absl::string_view key_l, uint64_t dummy_index) {
  uint8_t digest[SHA256_DIGEST_LENGTH];
  SHA256_CTX ctx;
  SHA256_Init(&ctx);
  SHA256_Update(&ctx, key_l.data(), key_l.size());
  uint8_t tag = 0x02;  // Tag for dummy blocks
  SHA256_Update(&ctx, &tag, sizeof(tag));
  SHA256_Update(&ctx, &dummy_index, sizeof(dummy_index));
  SHA256_Final(digest, &ctx);
  return std::string(reinterpret_cast<const char*>(digest), 16);
}

std::vector<uint8_t> GeneratePseudoRandomStream(absl::string_view key,
                                                absl::string_view label,
                                                size_t length) {
  std::vector<uint8_t> stream(length);
  size_t generated = 0;
  uint32_t counter = 0;

  while (generated < length) {
    SHA256_CTX ctx;
    SHA256_Init(&ctx);
    SHA256_Update(&ctx, key.data(), key.size());
    SHA256_Update(&ctx, label.data(), label.size());
    SHA256_Update(&ctx, &counter, sizeof(counter));
    uint8_t digest[SHA256_DIGEST_LENGTH];
    SHA256_Final(digest, &ctx);

    size_t to_copy = std::min<size_t>(sizeof(digest), length - generated);
    std::memcpy(stream.data() + generated, digest, to_copy);
    generated += to_copy;
    ++counter;
  }

  return stream;
}

std::vector<uint8_t> GenerateDummyBlock(absl::string_view key_d,
                                        absl::string_view label,
                                        size_t block_size) {
  return GeneratePseudoRandomStream(key_d, label, block_size);
}

std::vector<uint8_t> EncryptBlock(absl::string_view key_e,
                                  absl::string_view label, BlockId block_id,
                                  absl::Span<const uint8_t> payload,
                                  size_t payload_size, uint64_t version) {
  size_t total_size = CiphertextBlockSize(payload_size);
  std::vector<uint8_t> plaintext(total_size, 0);

  // 1. Pack Header: magic (4B) + block_id (8B) + version (8B) + payload_len
  // (4B) + checksum (4B) = 28B
  uint32_t magic = kBlockMagic;
  std::memcpy(plaintext.data(), &magic, sizeof(magic));
  std::memcpy(plaintext.data() + 4, &block_id, sizeof(block_id));
  std::memcpy(plaintext.data() + 12, &version, sizeof(version));

  uint32_t actual_payload_len =
      static_cast<uint32_t>(std::min(payload.size(), payload_size));
  std::memcpy(plaintext.data() + 20, &actual_payload_len,
              sizeof(actual_payload_len));

  // 2. Pack Payload
  if (actual_payload_len > 0) {
    std::memcpy(plaintext.data() + kBlockHeaderSize, payload.data(),
                actual_payload_len);
  }

  uint32_t checksum = ComputeChecksum(
      block_id, version,
      absl::MakeSpan(plaintext.data() + kBlockHeaderSize, payload_size));
  std::memcpy(plaintext.data() + 24, &checksum, sizeof(checksum));

  // 3. Encrypt via Keystream XOR
  std::vector<uint8_t> keystream =
      GeneratePseudoRandomStream(key_e, label, total_size);
  XorInPlace(absl::MakeSpan(plaintext), keystream);

  return plaintext;
}

absl::StatusOr<DecryptedBlock> DecryptBlock(
    absl::string_view key_e, absl::string_view label,
    absl::Span<const uint8_t> ciphertext, size_t payload_size) {
  size_t expected_size = CiphertextBlockSize(payload_size);
  if (ciphertext.size() < expected_size) {
    return absl::InvalidArgumentError(absl::StrFormat(
        "Ciphertext size mismatch: expected at least %zu, got %zu",
        expected_size, ciphertext.size()));
  }

  std::vector<uint8_t> plaintext(ciphertext.begin(),
                                 ciphertext.begin() + expected_size);
  std::vector<uint8_t> keystream =
      GeneratePseudoRandomStream(key_e, label, expected_size);
  XorInPlace(absl::MakeSpan(plaintext), keystream);

  // Validate magic header
  uint32_t magic = 0;
  std::memcpy(&magic, plaintext.data(), sizeof(magic));
  if (magic != kBlockMagic) {
    return absl::InvalidArgumentError(
        "Not a valid real block (magic mismatch).");
  }

  // Unpack Header
  BlockId block_id = 0;
  std::memcpy(&block_id, plaintext.data() + 4, sizeof(block_id));

  uint64_t version = 0;
  std::memcpy(&version, plaintext.data() + 12, sizeof(version));

  uint32_t payload_len = 0;
  std::memcpy(&payload_len, plaintext.data() + 20, sizeof(payload_len));

  if (payload_len > payload_size) {
    return absl::InternalError(absl::StrFormat(
        "Corrupt decrypted block: payload_len %u exceeds payload_size %zu",
        payload_len, payload_size));
  }

  uint32_t checksum = 0;
  std::memcpy(&checksum, plaintext.data() + 24, sizeof(checksum));

  std::vector<uint8_t> payload(payload_size, 0);
  if (payload_len > 0) {
    std::memcpy(payload.data(), plaintext.data() + kBlockHeaderSize,
                payload_len);
  }

  uint32_t computed_checksum =
      ComputeChecksum(block_id, version, absl::MakeSpan(payload));
  if (checksum != computed_checksum) {
    return absl::InvalidArgumentError(
        "Corrupt decrypted block: checksum mismatch.");
  }

  return DecryptedBlock{block_id, version, std::move(payload)};
}

std::vector<size_t> GenerateRandomPermutation(size_t n, uint64_t seed) {
  std::vector<size_t> perm(n);
  for (size_t i = 0; i < n; ++i) {
    perm[i] = i;
  }
  if (n <= 1) return perm;

  if (seed != 0) {
    std::mt19937_64 rng(seed);
    std::shuffle(perm.begin(), perm.end(), rng);
  } else {
    std::random_device rd;
    std::mt19937_64 rng(rd());
    std::shuffle(perm.begin(), perm.end(), rng);
  }
  return perm;
}

std::string SerializeOkvsTable(const std::vector<std::vector<uint8_t>>& table) {
  if (table.empty()) return "";
  uint32_t num_rows = static_cast<uint32_t>(table.size());
  uint32_t row_size = static_cast<uint32_t>(table[0].size());
  std::string out;
  out.resize(sizeof(num_rows) + sizeof(row_size) + num_rows * row_size);
  uint8_t* ptr = reinterpret_cast<uint8_t*>(out.data());
  std::memcpy(ptr, &num_rows, sizeof(num_rows));
  ptr += sizeof(num_rows);
  std::memcpy(ptr, &row_size, sizeof(row_size));
  ptr += sizeof(row_size);
  for (const auto& row : table) {
    std::memcpy(ptr, row.data(), row_size);
    ptr += row_size;
  }
  return out;
}

std::vector<std::vector<uint8_t>> DeserializeOkvsTable(
    absl::string_view bytes) {
  if (bytes.size() < sizeof(uint32_t) * 2) return {};
  uint32_t num_rows = 0;
  uint32_t row_size = 0;
  const uint8_t* ptr = reinterpret_cast<const uint8_t*>(bytes.data());
  std::memcpy(&num_rows, ptr, sizeof(num_rows));
  ptr += sizeof(num_rows);
  std::memcpy(&row_size, ptr, sizeof(row_size));
  ptr += sizeof(row_size);
  if (bytes.size() < sizeof(uint32_t) * 2 + num_rows * row_size) return {};

  std::vector<std::vector<uint8_t>> table(num_rows,
                                          std::vector<uint8_t>(row_size));
  for (uint32_t i = 0; i < num_rows; ++i) {
    std::memcpy(table[i].data(), ptr, row_size);
    ptr += row_size;
  }
  return table;
}

}  // namespace oram_okvs
