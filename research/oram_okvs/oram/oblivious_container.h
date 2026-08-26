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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_OBLIVIOUS_CONTAINER_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_OBLIVIOUS_CONTAINER_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "research/oram_okvs/okvs/polynomial_okvs.h"
#include "research/oram_okvs/okvs/random_band_okvs.h"
#include "research/oram_okvs/oram/crypto_utils.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

// ============================================================================
// Abstract Server-Side Oblivious Container Interface (Definition 1)
// ============================================================================
class ObliviousContainer {
 public:
  virtual ~ObliviousContainer() = default;

  // Evaluates query on this container and returns the retrieved/decoded block.
  virtual absl::StatusOr<std::vector<uint8_t>> ServerResponse(
      absl::string_view label) = 0;

  virtual size_t NumRealBlocks() const = 0;
  virtual size_t NumTotalBlocks() const = 0;
  virtual size_t StorageSizeBytes() const = 0;
};

// ============================================================================
// Keyed Oblivious Container (Appendix A.2)
// Used for Main container in LatORAM / OKLat (N real + MP dummies)
// ============================================================================
class KeyedObliviousContainer : public ObliviousContainer {
 public:
  KeyedObliviousContainer(
      absl::flat_hash_map<std::string, std::vector<uint8_t>> blocks,
      size_t num_real, size_t ciphertext_size, std::string metadata = "")
      : blocks_(std::move(blocks)),
        num_real_(num_real),
        ciphertext_size_(ciphertext_size),
        metadata_(std::move(metadata)) {}

  ~KeyedObliviousContainer() override = default;

  // ServerResponse: O(1) hash table lookup; if not found, returns zeros.
  absl::StatusOr<std::vector<uint8_t>> ServerResponse(
      absl::string_view label) override {
    auto it = blocks_.find(label);
    if (it != blocks_.end()) {
      std::vector<uint8_t> val = std::move(it->second);
      blocks_.erase(it);
      return val;
    }
    return std::vector<uint8_t>(ciphertext_size_, 0);
  }

  size_t NumRealBlocks() const override { return num_real_; }
  size_t NumTotalBlocks() const override { return blocks_.size(); }
  size_t StorageSizeBytes() const override {
    return blocks_.size() * (16 + ciphertext_size_);
  }

  const absl::flat_hash_map<std::string, std::vector<uint8_t>>& blocks() const {
    return blocks_;
  }
  absl::flat_hash_map<std::string, std::vector<uint8_t>>& mutable_blocks() {
    return blocks_;
  }

  // --- Client-Side Helpers (Appendix A.2) ---

  static absl::StatusOr<std::unique_ptr<KeyedObliviousContainer>> Build(
      absl::string_view seed_l, absl::string_view seed_d,
      absl::string_view key_enc,
      absl::Span<const std::pair<BlockId, std::vector<uint8_t>>> real_blocks,
      size_t num_dummies, size_t payload_size, uint64_t seed = 0) {
    size_t num_real = real_blocks.size();
    size_t total = num_real + num_dummies;
    size_t cipher_size = CiphertextBlockSize(payload_size);

    absl::flat_hash_map<std::string, std::vector<uint8_t>> blocks;
    blocks.reserve(total);

    // Real blocks
    for (size_t i = 0; i < num_real; ++i) {
      BlockId id = real_blocks[i].first;
      std::string label = GenerateBlockLabel(seed_l, id);
      std::vector<uint8_t> ct =
          EncryptBlock(key_enc, label, id, real_blocks[i].second, payload_size);
      blocks.emplace(std::move(label), std::move(ct));
    }

    // Dummy blocks
    for (size_t d = 1; d <= num_dummies; ++d) {
      std::string label = GenerateDummyLabel(seed_l, d);
      std::vector<uint8_t> dummy_blk =
          GenerateDummyBlock(seed_d, label, cipher_size);
      blocks.emplace(std::move(label), std::move(dummy_blk));
    }

    return std::make_unique<KeyedObliviousContainer>(std::move(blocks),
                                                     num_real, cipher_size);
  }

  static std::string ClientQuery(absl::string_view seed_l,
                                 std::optional<BlockId> query_id,
                                 uint64_t dummy_index) {
    if (query_id.has_value()) {
      return GenerateBlockLabel(seed_l, *query_id);
    }
    return GenerateDummyLabel(seed_l, dummy_index);
  }

  static absl::StatusOr<std::optional<DecryptedBlock>> ClientRecon(
      absl::string_view key_enc, absl::string_view label,
      absl::Span<const uint8_t> eBlk, size_t payload_size, bool is_dummy) {
    if (is_dummy) {
      return std::nullopt;
    }
    auto dec_or = DecryptBlock(key_enc, label, eBlk, payload_size);
    if (!dec_or.ok()) return dec_or.status();
    return *dec_or;
  }

 private:
  absl::flat_hash_map<std::string, std::vector<uint8_t>> blocks_;
  size_t num_real_;
  size_t ciphertext_size_;
  std::string metadata_;
};

// ============================================================================
// OKVS-Based Oblivious Container (Section 3.2)
// Used for Ring buckets (OKPRing) and lateral stashes (OKLat)
// ============================================================================
template <typename OkvsBackend =
              okvs::RandomBandOkvs<std::string, std::vector<uint8_t>>>
class OkvsObliviousContainer : public ObliviousContainer {
 public:
  OkvsObliviousContainer(std::vector<std::vector<uint8_t>> table, uint64_t seed,
                         size_t num_real, size_t ciphertext_size,
                         std::string metadata = "")
      : table_(std::move(table)),
        seed_(seed),
        num_real_(num_real),
        ciphertext_size_(ciphertext_size),
        metadata_(std::move(metadata)) {}

  ~OkvsObliviousContainer() override = default;

  // ServerResponse: Evaluates OKVS Decode on the fly.
  absl::StatusOr<std::vector<uint8_t>> ServerResponse(
      absl::string_view label) override {
    auto val_or =
        OkvsBackend::DecodeWithTable(std::string(label), table_, seed_);
    if (!val_or.ok()) return val_or.status();
    if (val_or->size() > ciphertext_size_) {
      val_or->resize(ciphertext_size_);
    }
    return *val_or;
  }

  size_t NumRealBlocks() const override { return num_real_; }
  size_t NumTotalBlocks() const override { return table_.size(); }
  size_t StorageSizeBytes() const override {
    if (table_.empty()) return 0;
    return table_.size() * table_[0].size();
  }

  const std::vector<std::vector<uint8_t>>& table() const { return table_; }
  uint64_t seed() const { return seed_; }
  const std::string& metadata() const { return metadata_; }

  // --- Client-Side Helpers (Section 3.2) ---

  static absl::StatusOr<std::unique_ptr<OkvsObliviousContainer<OkvsBackend>>>
  Build(absl::string_view seed_k, absl::string_view key_enc,
        absl::Span<const std::pair<BlockId, std::vector<uint8_t>>> real_blocks,
        size_t payload_size, uint64_t seed, bool with_metadata = false,
        double epsilon = 0.15, size_t fixed_size = 0) {
    size_t num_real = real_blocks.size();
    size_t cipher_size = CiphertextBlockSize(payload_size);

    std::vector<std::string> keys;
    std::vector<std::vector<uint8_t>> values;
    keys.reserve(num_real);
    values.reserve(num_real);

    std::vector<BlockId> initial_ids;
    initial_ids.reserve(num_real);

    for (size_t j = 0; j < num_real; ++j) {
      BlockId id = real_blocks[j].first;
      initial_ids.push_back(id);
      std::string label = GenerateBlockLabel(seed_k, id);
      std::vector<uint8_t> ct =
          EncryptBlock(key_enc, label, id, real_blocks[j].second, payload_size);
      keys.push_back(std::move(label));
      values.push_back(std::move(ct));
    }

    std::vector<std::vector<uint8_t>> table;
    uint64_t okvs_seed = seed;

    if constexpr (std::is_same_v<OkvsBackend,
                                 okvs::PolynomialOkvs<std::string,
                                                      std::vector<uint8_t>>>) {
      size_t target_size = (fixed_size > 0) ? fixed_size : num_real;
      okvs::PolynomialOkvs<std::string, std::vector<uint8_t>> okvs(
          target_size, seed, cipher_size);
      auto status = okvs.Encode(keys, values);
      if (!status.ok()) return status;
      table = okvs.Table();
      okvs_seed = okvs.seed();
    } else {
      okvs::RandomBandOkvs<std::string, std::vector<uint8_t>> okvs(epsilon,
                                                                   seed);
      auto status = okvs.EncodeWithRetries(keys, values);
      if (!status.ok()) return status;
      table = okvs.Table();
      okvs_seed = okvs.seed();
    }

    std::string md_plain;
    if (with_metadata) {
      uint32_t count = static_cast<uint32_t>(initial_ids.size());
      md_plain.resize(sizeof(count) + count * sizeof(BlockId));
      std::memcpy(md_plain.data(), &count, sizeof(count));
      if (count > 0) {
        std::memcpy(md_plain.data() + sizeof(count), initial_ids.data(),
                    count * sizeof(BlockId));
      }
    }

    return std::make_unique<OkvsObliviousContainer<OkvsBackend>>(
        std::move(table), okvs_seed, num_real, cipher_size,
        std::move(md_plain));
  }

  static std::string ClientQuery(absl::string_view seed_k,
                                 std::optional<BlockId> query_id) {
    if (query_id.has_value()) {
      return GenerateBlockLabel(seed_k, *query_id);
    }
    return GenerateRandomBytes(16);
  }

  static absl::StatusOr<std::optional<DecryptedBlock>> ClientRecon(
      absl::string_view key_enc, absl::string_view label,
      absl::Span<const uint8_t> eBlk, size_t payload_size, bool is_dummy) {
    if (is_dummy) {
      return std::nullopt;
    }
    auto dec_or = DecryptBlock(key_enc, label, eBlk, payload_size);
    if (!dec_or.ok()) return dec_or.status();
    return *dec_or;
  }

 private:
  std::vector<std::vector<uint8_t>> table_;
  uint64_t seed_;
  size_t num_real_;
  size_t ciphertext_size_;
  std::string metadata_;
};

// ============================================================================
// Factory Helpers from Protobuf Payloads
// ============================================================================

inline absl::StatusOr<std::unique_ptr<ObliviousContainer>>
CreateObliviousContainerFromStashRequest(
    const LatOramAddLateralStashRequest& request, size_t entry_size) {
  if (!request.okvs_table().empty()) {
    auto table = DeserializeOkvsTable(request.okvs_table());
    size_t actual_entry_size = table.empty() ? entry_size : table[0].size();
    if (request.okvs_backend() == 0) {
      return std::make_unique<OkvsObliviousContainer<
          okvs::PolynomialOkvs<std::string, std::vector<uint8_t>>>>(
          std::move(table), request.okvs_seed(), request.num_real(),
          actual_entry_size);
    } else {
      return std::make_unique<OkvsObliviousContainer<
          okvs::RandomBandOkvs<std::string, std::vector<uint8_t>>>>(
          std::move(table), request.okvs_seed(), request.num_real(),
          actual_entry_size);
    }
  }

  absl::flat_hash_map<std::string, std::vector<uint8_t>> blocks;
  blocks.reserve(request.blocks_size());
  for (const auto& blk : request.blocks()) {
    blocks.emplace(blk.label(), std::vector<uint8_t>(blk.ciphertext().begin(),
                                                     blk.ciphertext().end()));
  }
  return std::make_unique<KeyedObliviousContainer>(
      std::move(blocks), request.blocks_size(), entry_size);
}

inline absl::StatusOr<std::unique_ptr<ObliviousContainer>>
CreateObliviousContainerFromBucketProto(const RingOramBucketProto& bucket,
                                        size_t entry_size) {
  if (!bucket.okvs_table().empty()) {
    auto table = DeserializeOkvsTable(bucket.okvs_table());
    size_t actual_entry_size = table.empty() ? entry_size : table[0].size();
    if (bucket.okvs_backend() == 0) {
      return std::make_unique<OkvsObliviousContainer<
          okvs::PolynomialOkvs<std::string, std::vector<uint8_t>>>>(
          std::move(table), bucket.okvs_seed(), bucket.num_real(),
          actual_entry_size, bucket.metadata());
    } else {
      return std::make_unique<OkvsObliviousContainer<
          okvs::RandomBandOkvs<std::string, std::vector<uint8_t>>>>(
          std::move(table), bucket.okvs_seed(), bucket.num_real(),
          actual_entry_size, bucket.metadata());
    }
  }

  absl::flat_hash_map<std::string, std::vector<uint8_t>> blocks;
  blocks.reserve(bucket.blocks_size());
  for (int i = 0; i < bucket.blocks_size(); ++i) {
    const std::string& blk = bucket.blocks(i);
    blocks.emplace(std::to_string(i),
                   std::vector<uint8_t>(blk.begin(), blk.end()));
  }
  return std::make_unique<KeyedObliviousContainer>(
      std::move(blocks), bucket.blocks_size(), entry_size, bucket.metadata());
}

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_OBLIVIOUS_CONTAINER_H_
