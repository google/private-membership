// Copyright 2021 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "private_membership/rlwe/batch/cpp/parameters.h"

#include "absl/status/status.h"
#include "private_membership/rlwe/batch/proto/shared.pb.h"

namespace private_membership {
namespace batch {
namespace {
// Arbitrary, realistic maximum values to avoid integer overflows.
constexpr int kMaxLevelsOfRecursion = 16;
constexpr int kMaxLogCompressionFactor = 32;
constexpr int kMaxLogDegree = 32;
constexpr int kMaxLogDecompositionModulus = 32;
constexpr int kMaxLogT = 8;
constexpr int kMaxNumBucketsPerShard = (1 << 30);
constexpr int kMaxNumShards = (1 << 30);
constexpr int kMaxRequestModulusSize = 8;
constexpr int kMaxRequiredFakeQueries = (1 << 30);
constexpr int kMaxResponseModulusSize = 8;
constexpr int kMaxQueriesPerShard = (1 << 30);
constexpr int kMaxVariance = 256;

absl::Status ValidateCryptoParameters(
    const Parameters::CryptoParameters& parameters) {
  if (parameters.levels_of_recursion() <= 0 ||
      parameters.levels_of_recursion() > kMaxLevelsOfRecursion) {
    return absl::InvalidArgumentError(
        "levels_of_recursion is out of the valid range");
  }
  if (parameters.log_compression_factor() <= 0 ||
      parameters.log_compression_factor() >= kMaxLogCompressionFactor) {
    return absl::InvalidArgumentError(
        "log_compression_factor is out of the valid range");
  }
  if (parameters.log_decomposition_modulus() <= 0 ||
      parameters.log_decomposition_modulus() > kMaxLogDecompositionModulus) {
    return absl::InvalidArgumentError(
        "log_decomposition_modulus is out of the valid range");
  }
  if (parameters.log_degree() <= 0 || parameters.log_degree() > kMaxLogDegree) {
    return absl::InvalidArgumentError("log_degree is out of the valid range");
  }
  if (parameters.log_t() <= 0 || parameters.log_t() > kMaxLogT) {
    return absl::InvalidArgumentError("log_t is out of the valid range");
  }
  if (parameters.request_modulus_size() <= 0 ||
      parameters.request_modulus_size() > kMaxRequestModulusSize) {
    return absl::InvalidArgumentError(
        "Invalid number of request_modulus entries");
  }
  if (parameters.response_modulus_size() <= 0 ||
      parameters.response_modulus_size() > kMaxResponseModulusSize) {
    return absl::InvalidArgumentError(
        "Invalid number of response_modulus entries");
  }
  if (parameters.variance() <= 0 || parameters.variance() > kMaxVariance) {
    return absl::InvalidArgumentError("variance is out of the valid range");
  }
  return absl::OkStatus();
}

absl::Status ValidateShardParameters(
    const Parameters::ShardParameters& parameters) {
  if (parameters.number_of_shards() <= 0 ||
      parameters.number_of_shards() > kMaxNumShards) {
    return absl::InvalidArgumentError(
        "number_of_shards is out of the valid range");
  }
  if (parameters.number_of_buckets_per_shard() <= 0 ||
      parameters.number_of_buckets_per_shard() > kMaxNumBucketsPerShard) {
    return absl::InvalidArgumentError(
        "number_of_shards is out of the valid range");
  }
  if (parameters.required_fake_queries() < 0 ||
      parameters.required_fake_queries() > kMaxRequiredFakeQueries) {
    return absl::InvalidArgumentError(
        "required_fake_queries is out of the valid range");
  }
  if (parameters.required_queries_per_shard() < 0 ||
      parameters.required_queries_per_shard() > kMaxQueriesPerShard) {
    return absl::InvalidArgumentError(
        "required_queries_per_shard is out of the valid range");
  }
  if (std::numeric_limits<int>::max() / parameters.number_of_shards() <=
      parameters.required_queries_per_shard()) {
    return absl::InvalidArgumentError(
        "The total number of required queries is too large");
  }
  if (parameters.required_queries_per_shard() != 0 &&
      parameters.number_of_shards() * parameters.required_queries_per_shard() <=
          parameters.required_fake_queries()) {
    return absl::InvalidArgumentError(
        "The number of required fake queries is too large");
  }
  return absl::OkStatus();
}
}  // namespace

absl::Status ValidateParameters(const Parameters& parameters) {
  if (auto status = ValidateCryptoParameters(parameters.crypto_parameters());
      !status.ok()) {
    return status;
  }
  if (auto status = ValidateShardParameters(parameters.shard_parameters());
      !status.ok()) {
    return status;
  }

  return absl::OkStatus();
}

}  // namespace batch
}  // namespace private_membership
