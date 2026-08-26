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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_RING_RINGORAM_STORAGE_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_RING_RINGORAM_STORAGE_H_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "research/oram_okvs/oram/oblivious_container.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

// Unified concrete Ring ORAM server-side bucket tree storage.
// Supports both classical slot-index buckets and OKVS-encoded buckets
// (OKRing/OKPRing).
class RingOramStorage {
 public:
  RingOramStorage() = default;
  virtual ~RingOramStorage() = default;

  // Initializes the tree with tree_depth (L) and bucket_size (B).
  virtual absl::Status InitTree(
      size_t tree_depth, size_t bucket_size,
      std::vector<RingOramBucketProto> initial_buckets);

  // Reads metadata for all buckets along the path to leaf_id.
  virtual absl::StatusOr<RingOramReadPathMetadataResponse> ReadPathMetadata(
      uint32_t leaf_id);

  // Reads selected slot indices (standard) or OKVS queries (OKRing/OKPRing)
  // from path.
  virtual absl::StatusOr<RingOramReadPathBlocksResponse> ReadPathBlocks(
      uint32_t leaf_id, const std::vector<uint32_t>& slot_indices,
      const std::vector<std::string>& okvs_queries);

  // Reads full buckets along an eviction path.
  virtual absl::StatusOr<std::vector<RingOramBucketProto>> ReadPathFull(
      uint32_t leaf_id);

  // Writes back rebuilt buckets along an eviction path.
  virtual absl::Status WritePathFull(uint32_t leaf_id,
                                     std::vector<RingOramBucketProto> buckets);

  // Reshuffles a single bucket.
  virtual absl::Status ReshuffleBucket(uint32_t bucket_id,
                                       RingOramBucketProto bucket);

  size_t tree_depth() const { return tree_depth_; }
  size_t bucket_size() const { return bucket_size_; }
  size_t total_buckets() const { return buckets_.size(); }

 private:
  struct BucketEntry {
    RingOramBucketProto proto;
    std::unique_ptr<ObliviousContainer> container;
  };

  size_t tree_depth_ = 0;
  size_t bucket_size_ = 0;
  size_t ciphertext_size_ = 64;
  std::vector<BucketEntry> buckets_;
};

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_RINGORAM_STORAGE_H_
