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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_RING_RINGORAM_CLIENT_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_RING_RINGORAM_CLIENT_H_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/block_retrieval.h"
#include "research/oram_okvs/oram/crypto_utils.h"
#include "research/oram_okvs/oram/oram.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/ring/ringoram_server.h"
#include "research/oram_okvs/oram/ring/ringoram_storage.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

// Abstract base client implementation for Ring ORAM protocols.
//
// Implements common ORAM access, remapping, background eviction loop, and
// stash management while delegating scheme-specific query generation,
// bucket serialization, and block extraction to subclasses.
class RingOramClient : public Oram {
 public:
  RingOramClient(const RingOramConfig& config,
                 std::shared_ptr<RingOramServerInterface> server);

  ~RingOramClient() override = default;

  // Factory method returning the appropriate concrete Ring ORAM client
  // based on config.storage_type.
  static std::unique_ptr<RingOramClient> Create(
      const RingOramConfig& config,
      std::shared_ptr<RingOramServerInterface> server);

  // Initializes the Ring ORAM tree on the server.
  virtual absl::Status Init() = 0;

  // Builds and encrypts a complete bucket for the tree.
  virtual RingOramBucketProto BuildBucket(
      uint32_t bucket_id,
      const std::vector<std::pair<BlockId, std::vector<uint8_t>>>&
          real_blocks) = 0;

  // Reads the target block from the specified path on the server.
  virtual absl::StatusOr<std::vector<uint8_t>> ReadBlockFromPath(
      BlockId block_id, uint32_t old_leaf,
      const std::vector<uint32_t>& path) = 0;

  // Extracts all valid decrypted blocks from the path full response and pushes
  // them into the client stash during eviction.
  virtual absl::Status ExtractEvictionBlocks(
      const std::vector<uint32_t>& path,
      const RingOramReadPathFullResponse& full_resp) = 0;

  // Hook called after path buckets have been rebuilt and written to the server.
  virtual void OnEviction(const std::vector<uint32_t>& path,
                          const std::vector<std::vector<BlockId>>& placements) {
  }

  // Hook called when a block is accessed and pulled into client stash.
  virtual void OnBlockRead(BlockId block_id) {}

  // Oram interface methods.
  absl::StatusOr<std::vector<uint8_t>> Read(BlockId block_id) override;
  absl::Status Write(BlockId block_id, absl::Span<const uint8_t> data) override;
  absl::StatusOr<std::vector<uint8_t>> Access(
      BlockId block_id, AccessOp op,
      absl::Span<const uint8_t> data = {}) override;

  absl::StatusOr<std::vector<uint8_t>> Access(
      AccessOp op, BlockId block_id, absl::Span<const uint8_t> new_data = {});

  // Executes one eviction operation on the next path in reverse-lexicographic
  // order.
  absl::Status Evict();

  size_t num_blocks() const override { return config_.num_blocks; }
  size_t block_size() const override { return config_.block_size; }
  OramStats GetStats() const override { return stats_; }
  void ResetStats() override { stats_.Reset(); }

  const RingOramConfig& config() const { return config_; }
  size_t stash_size() const { return stash_.size(); }
  uint32_t eviction_counter() const { return eviction_counter_; }
  uint32_t GetLeaf(BlockId block_id) const { return position_map_[block_id]; }

  // Testing inspection hooks.
  virtual int32_t GetDistance(BlockId /*block_id*/) const { return -1; }
  virtual int32_t GetBucket(BlockId /*block_id*/) const { return -1; }
  virtual int32_t GetSlot(BlockId /*block_id*/) const { return -1; }

 protected:
  RingOramConfig config_;
  std::shared_ptr<RingOramServerInterface> server_;

  // Master cryptographic keys.
  std::string key_e_;  // Block payload encryption
  std::string key_l_;  // Block label / PRF key
  std::string key_d_;  // Dummy block encryption

  size_t tree_depth_;
  size_t num_leaves_;
  size_t bucket_size_;
  size_t ciphertext_block_size_;

  std::vector<uint32_t> position_map_;
  absl::flat_hash_map<BlockId, RingStashEntry> stash_;
  std::vector<uint64_t> block_versions_;
  uint64_t access_seq_ = 0;
  uint32_t eviction_counter_ = 0;

  std::mt19937_64 rng_;
  OramStats stats_;
};

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_RINGORAM_CLIENT_H_
