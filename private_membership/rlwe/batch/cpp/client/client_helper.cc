// Copyright 2021 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "private_membership/rlwe/batch/cpp/client/client_helper.h"

#include <cmath>
#include <cstdlib>

#include "absl/status/status.h"
#include "absl/strings/str_cat.h"
#include "private_membership/rlwe/batch/cpp/shared.pb.h"

namespace private_membership {
namespace batch {

// Arbitrary, realistic maximum values to avoid integer overflows.
constexpr int kMaximumLevelsOfRecursion = 30;
constexpr int kMaximumNumShards = (1 << 30);
constexpr int kMaximumNumBucketsPerShard = (1 << 30);

absl::StatusOr<BitVector> GenerateBitVector(
    const EncryptQueriesRequest& request) {
  const Parameters& parameters = request.parameters();

  // Retrieve all relevant parameters.
  int num_queries = request.plaintext_queries_size();
  int num_shards = parameters.shard_parameters().number_of_shards();
  int num_buckets_per_query =
      parameters.shard_parameters().number_of_buckets_per_shard();
  int levels_of_recursion =
      parameters.crypto_parameters().levels_of_recursion();
  if (num_shards <= 0 || num_buckets_per_query <= 0 ||
      levels_of_recursion <= 0) {
    return absl::InvalidArgumentError(
        "The input parameters are invalid because one or more parameters are "
        "non-negative.");
  } else if (num_shards >= kMaximumNumShards ||
             num_buckets_per_query > kMaximumNumBucketsPerShard ||
             levels_of_recursion > kMaximumLevelsOfRecursion) {
    return absl::InvalidArgumentError(
        "The input parameters are invalid because one or more parameters are "
        "too large.");
  }
  int num_slots_per_level = static_cast<int>(
      std::ceil(std::pow(num_buckets_per_query, 1.0 / levels_of_recursion)));
  int num_slots_per_query = num_slots_per_level * levels_of_recursion;
  int total_size = num_queries * num_slots_per_query;

  // Generate all indices that should have a 1-entry.
  BitVector bit_vector;
  bit_vector.indices.reserve(num_queries * levels_of_recursion);
  bit_vector.total_size = total_size;

  // Generate uncompressed vector of indices.
  for (int i = 0; i < request.plaintext_queries_size(); ++i) {
    const auto& plaintext_query = request.plaintext_queries(i);
    int shard_id = plaintext_query.query_metadata().shard_id();
    if (shard_id < 0 || shard_id >= num_shards) {
      return absl::InvalidArgumentError(absl::StrCat(
          "Invalid shard number. ", i, "-th plaintext query has shard ID of ",
          shard_id, " that is negative or larger than the maximum shard ID of ",
          num_shards - 1, "."));
    }
    int bucket_id = plaintext_query.bucket_id();
    if (bucket_id < 0 || bucket_id >= num_buckets_per_query) {
      return absl::InvalidArgumentError(absl::StrCat(
          "Invalid shard number. ", i, "-th plaintext query has bucket ID of ",
          bucket_id,
          " that is negative or larger than the maximum bucket ID of ",
          num_buckets_per_query - 1, "."));
    }
    int range_start = i * num_slots_per_query;
    int remaining_bucket_id = bucket_id;
    for (int j = 0; j < levels_of_recursion; ++j) {
      int index_at_level = remaining_bucket_id % num_slots_per_level;
      bit_vector.indices.push_back(range_start + index_at_level);
      remaining_bucket_id /= num_slots_per_level;
      range_start += num_slots_per_level;
    }
  }

  return bit_vector;
}

}  // namespace batch
}  // namespace private_membership
