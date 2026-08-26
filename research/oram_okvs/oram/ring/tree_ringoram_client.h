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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_RING_TREE_RINGORAM_CLIENT_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_RING_TREE_RINGORAM_CLIENT_H_

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/ring/ringoram_client.h"
#include "research/oram_okvs/oram/ring/ringoram_server.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

// Tree Ring ORAM client (Baseline - Algorithm 3).
//
// Maintains slot-level metadata and dummy blocks per bucket. On read, queries
// the exact slot index if the block is at that level, or queries a
// pseudo-random dummy slot otherwise.
class TreeRingOramClient : public RingOramClient {
 public:
  TreeRingOramClient(const RingOramConfig& config,
                     std::shared_ptr<RingOramServerInterface> server);

  absl::Status Init() override;

  RingOramBucketProto BuildBucket(
      uint32_t bucket_id,
      const std::vector<std::pair<BlockId, std::vector<uint8_t>>>& real_blocks)
      override;

  absl::StatusOr<std::vector<uint8_t>> ReadBlockFromPath(
      BlockId block_id, uint32_t old_leaf,
      const std::vector<uint32_t>& path) override;

  absl::Status ExtractEvictionBlocks(
      const std::vector<uint32_t>& path,
      const RingOramReadPathFullResponse& full_resp) override;

  void OnBlockRead(BlockId block_id) override;

  int32_t GetBucket(BlockId block_id) const override {
    return (block_id < block_bucket_.size()) ? block_bucket_[block_id] : -1;
  }

  int32_t GetSlot(BlockId block_id) const override {
    return (block_id < block_slot_.size()) ? block_slot_[block_id] : -1;
  }

 private:
  uint64_t GetBucketSeed(uint32_t bucket_id, uint64_t version) const;

  std::vector<int32_t> block_bucket_;
  std::vector<int32_t> block_slot_;
  std::vector<uint32_t> dummy_counts_;
  std::vector<uint64_t> bucket_versions_;
  std::vector<size_t> bucket_num_real_;
  std::vector<std::vector<BlockId>> bucket_blocks_;
};

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_TREE_RINGORAM_CLIENT_H_
