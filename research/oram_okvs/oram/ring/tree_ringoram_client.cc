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

#include "research/oram_okvs/oram/ring/tree_ringoram_client.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/crypto_utils.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/ring/ringoram_client.h"
#include "research/oram_okvs/oram/ring/ringoram_server.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

TreeRingOramClient::TreeRingOramClient(
    const RingOramConfig& config,
    std::shared_ptr<RingOramServerInterface> server)
    : RingOramClient(config, std::move(server)) {
  block_bucket_.assign(config_.num_blocks, -1);
  block_slot_.assign(config_.num_blocks, -1);
  dummy_counts_.assign(config_.TotalBuckets(), 0);
  bucket_versions_.assign(config_.TotalBuckets(), 0);
  bucket_num_real_.assign(config_.TotalBuckets(), 0);
  bucket_blocks_.resize(config_.TotalBuckets());
}

uint64_t TreeRingOramClient::GetBucketSeed(uint32_t bucket_id,
                                           uint64_t version) const {
  uint64_t base = config_.seed != 0 ? config_.seed : 123456789ULL;
  return base ^ (static_cast<uint64_t>(bucket_id) * 0x9e3779b97f4a7c15ULL) ^
         (version * 0x517cc1b727220a95ULL);
}

absl::Status TreeRingOramClient::Init() {
  stats_.Reset();
  stash_.clear();
  eviction_counter_ = 0;
  block_versions_.assign(config_.num_blocks, 0);
  access_seq_ = 0;

  block_bucket_.assign(config_.num_blocks, -1);
  block_slot_.assign(config_.num_blocks, -1);
  dummy_counts_.assign(config_.TotalBuckets(), 0);
  bucket_versions_.assign(config_.TotalBuckets(), 0);
  bucket_num_real_.assign(config_.TotalBuckets(), 0);
  for (auto& blocks : bucket_blocks_) {
    blocks.clear();
  }

  for (BlockId i = 0; i < config_.num_blocks; ++i) {
    position_map_[i] = static_cast<uint32_t>(rng_() % num_leaves_);
  }

  size_t total_buckets = config_.TotalBuckets();
  RingOramInitTreeRequest init_req;
  init_req.set_tree_depth(static_cast<uint32_t>(tree_depth_));
  init_req.set_bucket_size(static_cast<uint32_t>(bucket_size_));

  for (uint32_t b_id = 0; b_id < total_buckets; ++b_id) {
    std::vector<std::pair<BlockId, std::vector<uint8_t>>> empty_real_blocks;
    *init_req.add_initial_buckets() = BuildBucket(b_id, empty_real_blocks);
  }

  auto init_resp_or = server_->HandleInitTree(init_req);
  if (!init_resp_or.ok()) return init_resp_or.status();
  if (!init_resp_or->success()) {
    return absl::InternalError(
        absl::StrFormat("Server failed to initialize Ring ORAM tree: %s",
                        init_resp_or->error_message()));
  }

  stats_.Reset();
  return absl::OkStatus();
}

RingOramBucketProto TreeRingOramClient::BuildBucket(
    uint32_t bucket_id,
    const std::vector<std::pair<BlockId, std::vector<uint8_t>>>& real_blocks) {
  RingOramBucketProto bucket;
  bucket.set_bucket_id(bucket_id);
  bucket.set_access_count(0);
  dummy_counts_[bucket_id] = 0;
  bucket_num_real_[bucket_id] = real_blocks.size();
  bucket_blocks_[bucket_id].clear();

  uint64_t b_seed = GetBucketSeed(bucket_id, bucket_versions_[bucket_id]);
  std::vector<size_t> perm = GenerateRandomPermutation(bucket_size_, b_seed);

  bucket_versions_[bucket_id]++;

  std::vector<std::string> slot_blocks(bucket_size_);
  size_t real_idx = 0;
  uint32_t dummy_idx = 0;

  for (size_t slot = 0; slot < bucket_size_; ++slot) {
    size_t item_type = perm[slot];
    std::string label = absl::StrFormat("b_%u_%zu", bucket_id, slot);

    if (item_type < real_blocks.size()) {
      BlockId b_id = real_blocks[real_idx].first;
      const auto& payload = real_blocks[real_idx].second;
      uint64_t ver =
          (b_id < block_versions_.size()) ? block_versions_[b_id] : 0;
      std::vector<uint8_t> ct =
          EncryptBlock(key_e_, label, b_id, payload, config_.block_size, ver);
      slot_blocks[slot] =
          std::string(reinterpret_cast<const char*>(ct.data()), ct.size());

      block_bucket_[b_id] = bucket_id;
      block_slot_[b_id] = static_cast<int32_t>(slot);
      bucket_blocks_[bucket_id].push_back(b_id);
      real_idx++;
    } else {
      std::vector<uint8_t> dummy_ct =
          GenerateDummyBlock(key_d_, label, ciphertext_block_size_);
      slot_blocks[slot] = std::string(
          reinterpret_cast<const char*>(dummy_ct.data()), dummy_ct.size());
      dummy_idx++;
    }
  }

  for (auto& sb : slot_blocks) {
    *bucket.add_blocks() = std::move(sb);
  }

  return bucket;
}

absl::StatusOr<std::vector<uint8_t>> TreeRingOramClient::ReadBlockFromPath(
    BlockId block_id, uint32_t old_leaf, const std::vector<uint32_t>& path) {
  int target_depth = -1;
  int32_t target_b =
      (block_id < block_bucket_.size()) ? block_bucket_[block_id] : -1;
  if (target_b >= 0) {
    for (size_t d = 0; d < path.size(); ++d) {
      if (path[d] == static_cast<uint32_t>(target_b)) {
        target_depth = static_cast<int>(d);
        break;
      }
    }
  }

  RingOramReadPathBlocksRequest blk_req;
  blk_req.set_leaf_id(old_leaf);

  for (size_t d = 0; d < path.size(); ++d) {
    uint32_t b_id = path[d];
    uint32_t req_slot = 0;
    if (static_cast<int>(d) == target_depth) {
      req_slot = static_cast<uint32_t>(block_slot_[block_id]);
    } else {
      uint64_t b_seed = GetBucketSeed(b_id, bucket_versions_[b_id] - 1);
      std::vector<size_t> perm =
          GenerateRandomPermutation(bucket_size_, b_seed);
      size_t dummy_offset = bucket_num_real_[b_id] + dummy_counts_[b_id];
      dummy_counts_[b_id]++;

      for (size_t slot = 0; slot < bucket_size_; ++slot) {
        if (perm[slot] == dummy_offset) {
          req_slot = static_cast<uint32_t>(slot);
          break;
        }
      }
    }
    blk_req.add_slot_indices(req_slot);
  }

  size_t req_sz = blk_req.ByteSizeLong();
  stats_.bytes_written += req_sz;
  stats_.online_bytes_uploaded += req_sz;
  stats_.online_bytes_transferred += req_sz;

  auto blk_resp_or = server_->HandleReadPathBlocks(blk_req);
  if (!blk_resp_or.ok()) return blk_resp_or.status();
  size_t resp_sz = blk_resp_or->ByteSizeLong();
  stats_.bytes_read += resp_sz;
  stats_.online_bytes_downloaded += resp_sz;
  stats_.online_bytes_transferred += resp_sz;

  if (target_depth >= 0 && blk_resp_or->has_retrieval_response() &&
      target_depth <
          blk_resp_or->retrieval_response().candidate_blocks_size()) {
    const std::string& target_ct =
        blk_resp_or->retrieval_response().candidate_blocks(target_depth);
    uint32_t b_id = path[target_depth];
    int32_t slot = block_slot_[block_id];
    std::string label = absl::StrFormat("b_%u_%d", b_id, slot);

    auto dec_or =
        DecryptBlock(key_e_, label,
                     absl::Span<const uint8_t>(
                         reinterpret_cast<const uint8_t*>(target_ct.data()),
                         target_ct.size()),
                     config_.block_size);
    if (dec_or.ok()) {
      return dec_or->payload;
    }
  }

  return absl::NotFoundError("Block not in server path");
}

absl::Status TreeRingOramClient::ExtractEvictionBlocks(
    const std::vector<uint32_t>& path,
    const RingOramReadPathFullResponse& full_resp) {
  for (size_t d = 0; d < path.size(); ++d) {
    uint32_t b_id = path[d];
    if (static_cast<int>(d) >= full_resp.buckets_size()) continue;
    const auto& bucket = full_resp.buckets(d);

    for (BlockId b_id_val : bucket_blocks_[b_id]) {
      if (b_id_val >= block_bucket_.size() ||
          block_bucket_[b_id_val] != static_cast<int32_t>(b_id)) {
        continue;
      }

      int32_t slot = block_slot_[b_id_val];
      if (slot < 0 || slot >= bucket.blocks_size()) continue;

      std::string label = absl::StrFormat("b_%u_%d", b_id, slot);
      const std::string& ct = bucket.blocks(slot);
      auto dec_or = DecryptBlock(
          key_e_, label,
          absl::Span<const uint8_t>(reinterpret_cast<const uint8_t*>(ct.data()),
                                    ct.size()),
          config_.block_size);
      if (dec_or.ok()) {
        uint32_t cur_leaf = position_map_[b_id_val];
        stash_[b_id_val] = {b_id_val, cur_leaf, std::move(dec_or->payload)};
        block_bucket_[b_id_val] = -1;
        block_slot_[b_id_val] = -1;
      }
    }
    bucket_blocks_[b_id].clear();
  }
  return absl::OkStatus();
}

void TreeRingOramClient::OnBlockRead(BlockId block_id) {
  if (block_id < block_bucket_.size()) {
    block_bucket_[block_id] = -1;
    block_slot_[block_id] = -1;
  }
}

}  // namespace oram_okvs
