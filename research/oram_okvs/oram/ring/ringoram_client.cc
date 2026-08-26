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

#include "research/oram_okvs/oram/ring/ringoram_client.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/crypto_utils.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/ring/ok_ringoram_client.h"
#include "research/oram_okvs/oram/ring/okp_ringoram_client.h"
#include "research/oram_okvs/oram/ring/ringoram_server.h"
#include "research/oram_okvs/oram/ring/tree_ringoram_client.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

RingOramClient::RingOramClient(const RingOramConfig& config,
                               std::shared_ptr<RingOramServerInterface> server)
    : config_(config),
      server_(std::move(server)),
      tree_depth_(config.TreeDepth()),
      num_leaves_(config.NumLeaves()),
      bucket_size_(config.BucketSize()),
      ciphertext_block_size_(CiphertextBlockSize(config.block_size)),
      rng_(config.seed != 0 ? config.seed : std::random_device{}()) {
  key_e_ = GenerateRandomBytes(32);
  key_d_ = GenerateRandomBytes(32);
  key_l_ = GenerateRandomBytes(32);
  position_map_.resize(config_.num_blocks);
  block_versions_.assign(config_.num_blocks, 0);
}

std::unique_ptr<RingOramClient> RingOramClient::Create(
    const RingOramConfig& config,
    std::shared_ptr<RingOramServerInterface> server) {
  if (config.storage_type == RingOramStorageType::kOkRing) {
    return std::make_unique<OkRingOramClient>(config, std::move(server));
  } else if (config.storage_type == RingOramStorageType::kOkPRing) {
    return std::make_unique<OkPRingOramClient>(config, std::move(server));
  } else {
    return std::make_unique<TreeRingOramClient>(config, std::move(server));
  }
}

absl::StatusOr<std::vector<uint8_t>> RingOramClient::Read(BlockId block_id) {
  return Access(AccessOp::kRead, block_id);
}

absl::Status RingOramClient::Write(BlockId block_id,
                                   absl::Span<const uint8_t> data) {
  auto res_or = Access(AccessOp::kWrite, block_id, data);
  if (!res_or.ok()) return res_or.status();
  return absl::OkStatus();
}

absl::StatusOr<std::vector<uint8_t>> RingOramClient::Access(
    BlockId block_id, AccessOp op, absl::Span<const uint8_t> data) {
  return Access(op, block_id, data);
}

absl::StatusOr<std::vector<uint8_t>> RingOramClient::Access(
    AccessOp op, BlockId block_id, absl::Span<const uint8_t> new_data) {
  if (block_id >= config_.num_blocks) {
    return absl::InvalidArgumentError(absl::StrFormat(
        "Block ID %u out of bounds (max: %u)", block_id, config_.num_blocks));
  }

  stats_.total_accesses++;
  if (op == AccessOp::kRead) stats_.reads_count++;
  if (op == AccessOp::kWrite || op == AccessOp::kReadAndWrite) {
    stats_.writes_count++;
  }

  // 1. Lookup and Remap in position map
  uint32_t old_leaf = position_map_[block_id];
  uint32_t new_leaf = static_cast<uint32_t>(rng_() % num_leaves_);
  position_map_[block_id] = new_leaf;

  std::vector<uint32_t> path = GetPathBucketIds(old_leaf, tree_depth_);

  // 2. Delegate scheme-specific block retrieval along path
  std::vector<uint8_t> current_payload;
  auto payload_or = ReadBlockFromPath(block_id, old_leaf, path);
  if (payload_or.ok()) {
    current_payload = std::move(*payload_or);
  } else if (stash_.contains(block_id)) {
    current_payload = stash_[block_id].payload;
  } else {
    current_payload = std::vector<uint8_t>(config_.block_size, 0);
  }

  access_seq_++;
  if (block_id < block_versions_.size()) {
    block_versions_[block_id] = access_seq_;
  }

  // 3. Update Stash with target block mapped to new_leaf
  std::vector<uint8_t> stored_payload = current_payload;
  if (op == AccessOp::kWrite || op == AccessOp::kReadAndWrite) {
    stored_payload.assign(new_data.begin(), new_data.end());
    if (stored_payload.size() < config_.block_size) {
      stored_payload.resize(config_.block_size, 0);
    }
  }
  stash_[block_id] = {block_id, new_leaf, std::move(stored_payload)};
  OnBlockRead(block_id);

  // 4. Background Eviction: Evict once every A accesses
  if (config_.eviction_rate > 0 && (access_seq_ % config_.eviction_rate == 0)) {
    auto evict_status = Evict();
    if (!evict_status.ok()) return evict_status;
  }

  return current_payload;
}

absl::Status RingOramClient::Evict() {
  stats_.shuffles_count++;

  // 1. Determine eviction path using reverse-lexicographic order
  uint32_t evict_leaf = ReverseBits(eviction_counter_, tree_depth_);
  eviction_counter_ = (eviction_counter_ + 1) % num_leaves_;

  // 2. Read entire path from server
  RingOramReadPathFullRequest full_req;
  full_req.set_leaf_id(evict_leaf);
  size_t full_req_sz = full_req.ByteSizeLong();
  stats_.bytes_written += full_req_sz;
  stats_.amortized_bytes_uploaded += full_req_sz;
  stats_.amortized_bytes_transferred += full_req_sz;

  auto full_resp_or = server_->HandleReadPathFull(full_req);
  if (!full_resp_or.ok()) return full_resp_or.status();
  size_t full_resp_sz = full_resp_or->ByteSizeLong();
  stats_.bytes_read += full_resp_sz;
  stats_.amortized_bytes_downloaded += full_resp_sz;
  stats_.amortized_bytes_transferred += full_resp_sz;

  std::vector<uint32_t> path = GetPathBucketIds(evict_leaf, tree_depth_);

  // 3. Extract and decrypt all valid blocks in path into Stash
  auto extract_status = ExtractEvictionBlocks(path, *full_resp_or);
  if (!extract_status.ok()) return extract_status;

  // 4. Rebuild path bottom-up (level L down to 0)
  std::vector<RingOramBucketProto> rebuilt_buckets(tree_depth_ + 1);
  std::vector<std::vector<BlockId>> placements(tree_depth_ + 1);

  for (int d = static_cast<int>(tree_depth_); d >= 0; --d) {
    uint32_t b_id = path[d];
    std::vector<std::pair<BlockId, std::vector<uint8_t>>> bucket_real_blocks;

    // Find all blocks in stash that can reside in bucket b_id at level d
    std::vector<BlockId> candidate_block_ids;
    for (const auto& [b_id_val, entry] : stash_) {
      uint32_t subtree_prefix = entry.leaf_id >> (tree_depth_ - d);
      uint32_t bucket_prefix = b_id - ((1 << d) - 1);
      if (subtree_prefix == bucket_prefix) {
        candidate_block_ids.push_back(b_id_val);
      }
    }

    // Pop up to Z blocks from Stash for this bucket
    size_t take_count = std::min(candidate_block_ids.size(), config_.z_real);
    for (size_t i = 0; i < take_count; ++i) {
      BlockId id = candidate_block_ids[i];
      placements[d].push_back(id);
      bucket_real_blocks.push_back({id, stash_[id].payload});
      stash_.erase(id);
    }

    rebuilt_buckets[d] = BuildBucket(b_id, bucket_real_blocks);
  }

  OnEviction(path, placements);

  // 5. Write rebuilt path back to server
  RingOramWritePathFullRequest write_req;
  write_req.set_leaf_id(evict_leaf);
  for (size_t d = 0; d <= tree_depth_; ++d) {
    *write_req.add_buckets() = std::move(rebuilt_buckets[d]);
  }

  size_t write_req_sz = write_req.ByteSizeLong();
  stats_.bytes_written += write_req_sz;
  stats_.amortized_bytes_uploaded += write_req_sz;
  stats_.amortized_bytes_transferred += write_req_sz;

  auto write_resp_or = server_->HandleWritePathFull(write_req);
  if (!write_resp_or.ok()) return write_resp_or.status();
  size_t write_resp_sz = write_resp_or->ByteSizeLong();
  stats_.bytes_read += write_resp_sz;
  stats_.amortized_bytes_downloaded += write_resp_sz;
  stats_.amortized_bytes_transferred += write_resp_sz;

  return absl::OkStatus();
}

}  // namespace oram_okvs
