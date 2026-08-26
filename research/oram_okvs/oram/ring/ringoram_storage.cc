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

#include "research/oram_okvs/oram/ring/ringoram_storage.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "research/oram_okvs/oram/crypto_utils.h"
#include "research/oram_okvs/oram/oblivious_container.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

absl::Status RingOramStorage::InitTree(
    size_t tree_depth, size_t bucket_size,
    std::vector<RingOramBucketProto> initial_buckets) {
  tree_depth_ = tree_depth;
  bucket_size_ = bucket_size;
  size_t expected_buckets = (size_t{1} << (tree_depth + 1)) - 1;

  if (initial_buckets.size() != expected_buckets) {
    return absl::InvalidArgumentError(
        absl::StrFormat("Expected %zu buckets for depth %zu, got %zu",
                        expected_buckets, tree_depth, initial_buckets.size()));
  }

  // Determine ciphertext size
  ciphertext_size_ = 64;
  for (const auto& b : initial_buckets) {
    if (b.blocks_size() > 0 && !b.blocks(0).empty()) {
      ciphertext_size_ = b.blocks(0).size();
      break;
    }
    if (!b.okvs_table().empty()) {
      auto table = DeserializeOkvsTable(b.okvs_table());
      if (!table.empty() && !table[0].empty()) {
        ciphertext_size_ = table[0].size();
        break;
      }
    }
  }

  buckets_.clear();
  buckets_.reserve(expected_buckets);

  for (auto& b : initial_buckets) {
    BucketEntry entry;
    auto container_or =
        CreateObliviousContainerFromBucketProto(b, ciphertext_size_);
    if (container_or.ok()) {
      entry.container = std::move(*container_or);
    }
    entry.proto = std::move(b);
    buckets_.push_back(std::move(entry));
  }

  return absl::OkStatus();
}

absl::StatusOr<RingOramReadPathMetadataResponse>
RingOramStorage::ReadPathMetadata(uint32_t leaf_id) {
  size_t num_leaves = size_t{1} << tree_depth_;
  if (leaf_id >= num_leaves) {
    return absl::InvalidArgumentError(absl::StrFormat(
        "Leaf ID %u out of bounds (num_leaves: %zu)", leaf_id, num_leaves));
  }

  RingOramReadPathMetadataResponse response;
  std::vector<uint32_t> path = GetPathBucketIds(leaf_id, tree_depth_);

  for (uint32_t b_id : path) {
    if (b_id >= buckets_.size()) {
      response.set_success(false);
      response.set_error_message(
          absl::StrFormat("Bucket ID %u out of range", b_id));
      return response;
    }
    response.add_bucket_ids(b_id);
    response.add_access_counts(buckets_[b_id].proto.access_count());
    response.add_bucket_metadata(buckets_[b_id].proto.metadata());
  }

  response.set_success(true);
  return response;
}

absl::StatusOr<RingOramReadPathBlocksResponse> RingOramStorage::ReadPathBlocks(
    uint32_t leaf_id, const std::vector<uint32_t>& slot_indices,
    const std::vector<std::string>& okvs_queries) {
  size_t num_leaves = size_t{1} << tree_depth_;
  if (leaf_id >= num_leaves) {
    return absl::InvalidArgumentError(absl::StrFormat(
        "Leaf ID %u out of bounds (num_leaves: %zu)", leaf_id, num_leaves));
  }

  std::vector<uint32_t> path = GetPathBucketIds(leaf_id, tree_depth_);
  RingOramReadPathBlocksResponse response;

  if (!okvs_queries.empty()) {
    // OKVS Bucket Evaluation
    response.mutable_retrieval_response()->set_mode(
        RetrievalModeProto::RETRIEVAL_MODE_DIRECT);
    for (size_t d = 0; d < path.size(); ++d) {
      uint32_t b_id = path[d];
      if (b_id >= buckets_.size()) {
        response.set_success(false);
        response.set_error_message("Bucket index out of range");
        return response;
      }

      std::string query_label =
          (okvs_queries.size() == 1) ? okvs_queries[0] : okvs_queries[d];

      std::vector<uint8_t> block_bytes;
      if (buckets_[b_id].container) {
        auto val_or = buckets_[b_id].container->ServerResponse(query_label);
        if (!val_or.ok()) {
          response.set_success(false);
          response.set_error_message(std::string(val_or.status().message()));
          return response;
        }
        block_bytes = std::move(*val_or);
      } else {
        block_bytes.assign(ciphertext_size_, 0);
      }

      buckets_[b_id].proto.set_access_count(
          buckets_[b_id].proto.access_count() + 1);
      response.mutable_retrieval_response()->add_candidate_blocks(
          block_bytes.data(), block_bytes.size());
    }
  } else if (!slot_indices.empty()) {
    // Classical slot index retrieval
    if (slot_indices.size() != path.size()) {
      return absl::InvalidArgumentError(absl::StrFormat(
          "Slot indices count (%zu) must match path length (%zu)",
          slot_indices.size(), path.size()));
    }

    response.mutable_retrieval_response()->set_mode(
        RetrievalModeProto::RETRIEVAL_MODE_DIRECT);
    for (size_t d = 0; d < path.size(); ++d) {
      uint32_t b_id = path[d];
      uint32_t slot = slot_indices[d];

      if (b_id >= buckets_.size()) {
        response.set_success(false);
        response.set_error_message("Bucket index out of range");
        return response;
      }

      if (static_cast<int>(slot) >= buckets_[b_id].proto.blocks_size()) {
        response.set_success(false);
        response.set_error_message("Slot index out of range");
        return response;
      }

      buckets_[b_id].proto.set_access_count(
          buckets_[b_id].proto.access_count() + 1);
      const auto& blk = buckets_[b_id].proto.blocks(slot);
      response.mutable_retrieval_response()->add_candidate_blocks(blk.data(),
                                                                  blk.size());
    }
  }

  response.set_success(true);
  return response;
}

absl::StatusOr<std::vector<RingOramBucketProto>> RingOramStorage::ReadPathFull(
    uint32_t leaf_id) {
  std::vector<uint32_t> path = GetPathBucketIds(leaf_id, tree_depth_);
  std::vector<RingOramBucketProto> result;
  result.reserve(path.size());

  for (uint32_t b_id : path) {
    if (b_id >= buckets_.size()) {
      return absl::InvalidArgumentError(
          absl::StrFormat("Bucket ID %u out of range", b_id));
    }
    result.push_back(buckets_[b_id].proto);
  }

  return result;
}

absl::Status RingOramStorage::WritePathFull(
    uint32_t leaf_id, std::vector<RingOramBucketProto> path_buckets) {
  std::vector<uint32_t> path = GetPathBucketIds(leaf_id, tree_depth_);
  if (path_buckets.size() != path.size()) {
    return absl::InvalidArgumentError(
        absl::StrFormat("WritePathFull: expected %zu buckets, got %zu",
                        path.size(), path_buckets.size()));
  }

  for (size_t d = 0; d < path.size(); ++d) {
    uint32_t b_id = path[d];
    if (b_id >= buckets_.size()) {
      return absl::InvalidArgumentError("Bucket index out of range");
    }
    auto container_or = CreateObliviousContainerFromBucketProto(
        path_buckets[d], ciphertext_size_);
    if (container_or.ok()) {
      buckets_[b_id].container = std::move(*container_or);
    }
    buckets_[b_id].proto = std::move(path_buckets[d]);
  }

  return absl::OkStatus();
}

absl::Status RingOramStorage::ReshuffleBucket(uint32_t bucket_id,
                                              RingOramBucketProto bucket) {
  if (bucket_id >= buckets_.size()) {
    return absl::InvalidArgumentError(
        absl::StrFormat("Bucket ID %u out of range", bucket_id));
  }
  auto container_or =
      CreateObliviousContainerFromBucketProto(bucket, ciphertext_size_);
  if (container_or.ok()) {
    buckets_[bucket_id].container = std::move(*container_or);
  }
  buckets_[bucket_id].proto = std::move(bucket);
  return absl::OkStatus();
}

}  // namespace oram_okvs
