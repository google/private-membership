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

#include "research/oram_okvs/oram/ring/ringoram_server.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/ring/ringoram_storage.h"
#include "research/oram_okvs/pir/onion_pir_params.h"
#include "research/oram_okvs/pir/onion_pir_server.h"

namespace oram_okvs {

RingOramServer::RingOramServer(std::unique_ptr<RingOramStorage> storage)
    : storage_(std::move(storage)),
      retrieval_server_(BlockRetrievalServer::Create()) {
  if (!storage_) {
    storage_ = std::make_unique<RingOramStorage>();
  }
}

absl::StatusOr<RingOramInitTreeResponse> RingOramServer::HandleInitTree(
    const RingOramInitTreeRequest& request) {
  RingOramInitTreeResponse response;
  std::vector<RingOramBucketProto> buckets;
  buckets.reserve(request.initial_buckets_size());
  for (const auto& b : request.initial_buckets()) {
    buckets.push_back(b);
  }

  auto status = storage_->InitTree(request.tree_depth(), request.bucket_size(),
                                   std::move(buckets));
  if (!status.ok()) {
    response.set_success(false);
    response.set_error_message(std::string(status.message()));
    return response;
  }

  // Initialize retrieval server if evaluation keys are provided
  if (request.has_pir_eval_keys()) {
    auto set_keys_status =
        retrieval_server_->SetEvaluationKeys(request.pir_eval_keys());
    if (!set_keys_status.ok()) {
      response.set_success(false);
      response.set_error_message(std::string(set_keys_status.message()));
      return response;
    }
  }

  response.set_success(true);
  return response;
}

absl::StatusOr<RingOramReadPathMetadataResponse>
RingOramServer::HandleReadPathMetadata(
    const RingOramReadPathMetadataRequest& request) {
  return storage_->ReadPathMetadata(request.leaf_id());
}

absl::StatusOr<RingOramReadPathBlocksResponse>
RingOramServer::HandleReadPathBlocks(
    const RingOramReadPathBlocksRequest& request) {
  auto t0 = std::chrono::steady_clock::now();
  std::vector<uint32_t> slot_indices;
  slot_indices.reserve(request.slot_indices_size());
  for (uint32_t s : request.slot_indices()) {
    slot_indices.push_back(s);
  }

  std::vector<std::string> okvs_queries;
  okvs_queries.reserve(request.okvs_queries_size());
  for (const auto& q : request.okvs_queries()) {
    okvs_queries.push_back(q);
  }

  auto resp_or =
      storage_->ReadPathBlocks(request.leaf_id(), slot_indices, okvs_queries);
  if (!resp_or.ok()) return resp_or;

  if (resp_or->success()) {
    total_accesses_served_++;

    if (request.has_retrieval_query() && retrieval_server_ &&
        resp_or->has_retrieval_response()) {
      std::vector<std::vector<uint8_t>> path_blocks;
      for (const auto& blk : resp_or->retrieval_response().candidate_blocks()) {
        path_blocks.emplace_back(blk.begin(), blk.end());
      }
      auto ret_resp_or =
          retrieval_server_->Evaluate(path_blocks, request.retrieval_query());
      if (!ret_resp_or.ok()) {
        resp_or->set_success(false);
        resp_or->set_error_message(std::string(ret_resp_or.status().message()));
        return resp_or;
      }
      *resp_or->mutable_retrieval_response() = std::move(*ret_resp_or);
    }
  }
  auto t1 = std::chrono::steady_clock::now();
  total_server_time_ms_ +=
      std::chrono::duration<double, std::milli>(t1 - t0).count();
  return resp_or;
}

absl::StatusOr<RingOramReadPathFullResponse> RingOramServer::HandleReadPathFull(
    const RingOramReadPathFullRequest& request) {
  auto t0 = std::chrono::steady_clock::now();
  RingOramReadPathFullResponse response;
  auto buckets_or = storage_->ReadPathFull(request.leaf_id());
  if (!buckets_or.ok()) {
    response.set_success(false);
    response.set_error_message(std::string(buckets_or.status().message()));
    return response;
  }

  for (auto& b : *buckets_or) {
    *response.add_buckets() = std::move(b);
  }

  response.set_success(true);
  auto t1 = std::chrono::steady_clock::now();
  total_server_time_ms_ +=
      std::chrono::duration<double, std::milli>(t1 - t0).count();
  return response;
}

absl::StatusOr<RingOramWritePathFullResponse>
RingOramServer::HandleWritePathFull(
    const RingOramWritePathFullRequest& request) {
  auto t0 = std::chrono::steady_clock::now();
  RingOramWritePathFullResponse response;
  std::vector<RingOramBucketProto> buckets;
  buckets.reserve(request.buckets_size());
  for (const auto& b : request.buckets()) {
    buckets.push_back(b);
  }

  auto status = storage_->WritePathFull(request.leaf_id(), std::move(buckets));
  if (!status.ok()) {
    response.set_success(false);
    response.set_error_message(std::string(status.message()));
    return response;
  }

  total_evictions_served_++;
  response.set_success(true);
  auto t1 = std::chrono::steady_clock::now();
  total_server_time_ms_ +=
      std::chrono::duration<double, std::milli>(t1 - t0).count();
  return response;
}

absl::StatusOr<RingOramReshuffleBucketResponse>
RingOramServer::HandleReshuffleBucket(
    const RingOramReshuffleBucketRequest& request) {
  RingOramReshuffleBucketResponse response;
  auto status =
      storage_->ReshuffleBucket(request.bucket_id(), request.bucket());
  if (!status.ok()) {
    response.set_success(false);
    response.set_error_message(std::string(status.message()));
    return response;
  }

  response.set_success(true);
  return response;
}

}  // namespace oram_okvs
