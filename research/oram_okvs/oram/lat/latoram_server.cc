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

#include "research/oram_okvs/oram/lat/latoram_server.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/crypto_utils.h"
#include "research/oram_okvs/oram/lat/latoram_storage.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/types.h"
#include "research/oram_okvs/pir/onion_pir_params.h"
#include "research/oram_okvs/pir/onion_pir_server.h"

namespace oram_okvs {

LatOramServer::LatOramServer(std::unique_ptr<LatOramStorage> storage)
    : storage_(std::move(storage)),
      retrieval_server_(BlockRetrievalServer::Create()) {
  if (!storage_) {
    storage_ = std::make_unique<LatOramStorage>();
  }
}

absl::StatusOr<LatOramInitMainResponse> LatOramServer::HandleInitMain(
    const LatOramInitMainRequest& request) {
  LatOramInitMainResponse response;
  std::vector<LabeledBlock> blocks;
  blocks.reserve(request.blocks_size());
  for (const auto& blk_proto : request.blocks()) {
    blocks.emplace_back(blk_proto.label(),
                        std::vector<uint8_t>(blk_proto.ciphertext().begin(),
                                             blk_proto.ciphertext().end()));
  }

  auto status = storage_->InitMain(std::move(blocks));
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

absl::StatusOr<LatOramAddLateralStashResponse>
LatOramServer::HandleAddLateralStash(
    const LatOramAddLateralStashRequest& request) {
  auto t0 = std::chrono::steady_clock::now();
  LatOramAddLateralStashResponse response;
  std::vector<LabeledBlock> blocks;
  blocks.reserve(request.blocks_size());
  for (const auto& blk_proto : request.blocks()) {
    blocks.emplace_back(blk_proto.label(),
                        std::vector<uint8_t>(blk_proto.ciphertext().begin(),
                                             blk_proto.ciphertext().end()));
  }

  if (!request.okvs_table().empty()) {
    auto table = DeserializeOkvsTable(request.okvs_table());
    uint64_t seed = request.okvs_seed();
    uint32_t num_real = request.num_real();
    uint32_t backend = request.okvs_backend();
    auto status = storage_->AddLateralStashOkvs(
        std::move(table), seed, num_real, backend, std::move(blocks));
    if (!status.ok()) {
      response.set_success(false);
      response.set_error_message(std::string(status.message()));
      return response;
    }
    response.set_success(true);
    auto t1 = std::chrono::steady_clock::now();
    total_server_time_ms_ +=
        std::chrono::duration<double, std::milli>(t1 - t0).count();
    return response;
  }

  auto status = storage_->AddLateralStash(std::move(blocks));
  if (!status.ok()) {
    response.set_success(false);
    response.set_error_message(std::string(status.message()));
    return response;
  }

  response.set_success(true);
  auto t1 = std::chrono::steady_clock::now();
  total_server_time_ms_ +=
      std::chrono::duration<double, std::milli>(t1 - t0).count();
  return response;
}

absl::StatusOr<LatOramQueryResponse> LatOramServer::HandleQuery(
    const LatOramQueryRequest& request) {
  auto t0 = std::chrono::steady_clock::now();
  LatOramQueryResponse response;
  std::vector<std::string> labels;
  labels.reserve(request.labels_size());
  for (const auto& l : request.labels()) {
    labels.push_back(l);
  }

  auto containers_or = storage_->QueryContainers(labels);
  if (!containers_or.ok()) {
    response.set_success(false);
    response.set_error_message(std::string(containers_or.status().message()));
    return response;
  }

  // Process unified retrieval query
  if (request.has_retrieval_query() && retrieval_server_) {
    auto ret_resp_or =
        retrieval_server_->Evaluate(*containers_or, request.retrieval_query());
    if (!ret_resp_or.ok()) {
      response.set_success(false);
      response.set_error_message(std::string(ret_resp_or.status().message()));
      return response;
    }
    *response.mutable_retrieval_response() = std::move(*ret_resp_or);
  } else if (retrieval_server_) {
    ObliviousRetrievalQueryProto q;
    q.set_mode(labels.size() > 1 ? RetrievalModeProto::RETRIEVAL_MODE_XOR
                                 : RetrievalModeProto::RETRIEVAL_MODE_DIRECT);
    auto ret_resp_or = retrieval_server_->Evaluate(*containers_or, q);
    if (!ret_resp_or.ok()) {
      response.set_success(false);
      response.set_error_message(std::string(ret_resp_or.status().message()));
      return response;
    }
    *response.mutable_retrieval_response() = std::move(*ret_resp_or);
  }

  total_queries_served_++;
  if (labels.size() > 1) {
    total_xor_ops_ += labels.size() - 1;
  }

  response.set_success(true);
  auto t1 = std::chrono::steady_clock::now();
  total_server_time_ms_ +=
      std::chrono::duration<double, std::milli>(t1 - t0).count();
  return response;
}

absl::StatusOr<CacheShuffleInitResponse> LatOramServer::HandleCacheShuffleInit(
    const CacheShuffleInitRequest& request) {
  CacheShuffleInitResponse response;
  auto count_or = storage_->InitCacheShuffle(
      request.source_container(), request.num_rows(), request.num_cols());
  if (!count_or.ok()) {
    response.set_success(false);
    response.set_error_message(std::string(count_or.status().message()));
    return response;
  }

  total_shuffles_served_++;
  response.set_total_blocks(static_cast<uint32_t>(*count_or));
  response.set_success(true);
  return response;
}

absl::StatusOr<CacheShuffleSprayStepResponse>
LatOramServer::HandleCacheShuffleSprayStep(
    const CacheShuffleSprayStepRequest& request) {
  CacheShuffleSprayStepResponse response;

  // 1. If client uploaded a column, save it to matrix T
  if (request.uploaded_column_size() > 0) {
    std::vector<LabeledBlock> col;
    col.reserve(request.uploaded_column_size());
    for (const auto& blk_proto : request.uploaded_column()) {
      col.emplace_back(blk_proto.label(),
                       std::vector<uint8_t>(blk_proto.ciphertext().begin(),
                                            blk_proto.ciphertext().end()));
    }
    auto status =
        storage_->CacheShuffleWriteColumn(request.col_idx(), std::move(col));
    if (!status.ok()) {
      response.set_success(false);
      response.set_error_message(std::string(status.message()));
      return response;
    }
  }

  // 2. If client requested blocks to download for the next spray round
  if (request.download_count() > 0) {
    auto blocks_or = storage_->CacheShuffleReadSource(request.download_count());
    if (!blocks_or.ok()) {
      response.set_success(false);
      response.set_error_message(std::string(blocks_or.status().message()));
      return response;
    }
    for (const auto& [label, blk] : *blocks_or) {
      auto* p = response.add_downloaded_blocks();
      p->set_label(label);
      p->set_ciphertext(blk.data(), blk.size());
    }
  }

  response.set_success(true);
  return response;
}

absl::StatusOr<CacheShuffleRecalibrateStepResponse>
LatOramServer::HandleCacheShuffleRecalibrateStep(
    const CacheShuffleRecalibrateStepRequest& request) {
  CacheShuffleRecalibrateStepResponse response;

  // 1. If client uploaded recalibrated blocks from previous row, save to output
  // container
  if (request.recalibrated_blocks_size() > 0) {
    std::vector<LabeledBlock> blocks;
    blocks.reserve(request.recalibrated_blocks_size());
    for (const auto& blk_proto : request.recalibrated_blocks()) {
      blocks.emplace_back(blk_proto.label(),
                          std::vector<uint8_t>(blk_proto.ciphertext().begin(),
                                               blk_proto.ciphertext().end()));
    }
    auto status = storage_->CacheShuffleAppendOutput(std::move(blocks));
    if (!status.ok()) {
      response.set_success(false);
      response.set_error_message(std::string(status.message()));
      return response;
    }
  }

  // 2. Fetch row from matrix T
  auto row_or = storage_->CacheShuffleReadRow(request.row_idx());
  if (!row_or.ok()) {
    response.set_success(false);
    response.set_error_message(std::string(row_or.status().message()));
    return response;
  }

  for (const auto& [label, blk] : *row_or) {
    auto* p = response.add_row_blocks();
    p->set_label(label);
    p->set_ciphertext(blk.data(), blk.size());
  }

  response.set_success(true);
  return response;
}

absl::StatusOr<TwoMergeInitResponse> LatOramServer::HandleTwoMergeInit(
    const TwoMergeInitRequest& request) {
  TwoMergeInitResponse response;
  auto sizes_or =
      storage_->InitTwoMerge(request.container1(), request.container2());
  if (!sizes_or.ok()) {
    response.set_success(false);
    response.set_error_message(std::string(sizes_or.status().message()));
    return response;
  }

  response.set_size1(static_cast<uint32_t>(sizes_or->first));
  response.set_size2(static_cast<uint32_t>(sizes_or->second));
  response.set_success(true);
  return response;
}

absl::StatusOr<TwoMergeStepResponse> LatOramServer::HandleTwoMergeStep(
    const TwoMergeStepRequest& request) {
  TwoMergeStepResponse response;

  // 1. If client uploaded a merged batch, append to output container
  if (request.merged_batch_size() > 0) {
    std::vector<LabeledBlock> blocks;
    blocks.reserve(request.merged_batch_size());
    for (const auto& blk_proto : request.merged_batch()) {
      blocks.emplace_back(blk_proto.label(),
                          std::vector<uint8_t>(blk_proto.ciphertext().begin(),
                                               blk_proto.ciphertext().end()));
    }
    auto status = storage_->TwoMergeAppendOutput(std::move(blocks));
    if (!status.ok()) {
      response.set_success(false);
      response.set_error_message(std::string(status.message()));
      return response;
    }
  }

  // 2. Read chunks from C1 and C2
  if (request.fetch_c1_count() > 0 || request.fetch_c2_count() > 0) {
    auto chunks_or = storage_->TwoMergeReadChunks(request.fetch_c1_count(),
                                                  request.fetch_c2_count());
    if (!chunks_or.ok()) {
      response.set_success(false);
      response.set_error_message(std::string(chunks_or.status().message()));
      return response;
    }
    for (const auto& [label, blk] : chunks_or->first) {
      auto* p = response.add_c1_blocks();
      p->set_label(label);
      p->set_ciphertext(blk.data(), blk.size());
    }
    for (const auto& [label, blk] : chunks_or->second) {
      auto* p = response.add_c2_blocks();
      p->set_label(label);
      p->set_ciphertext(blk.data(), blk.size());
    }
  }

  response.set_success(true);
  return response;
}

absl::StatusOr<TwoMergeFinalizeResponse> LatOramServer::HandleTwoMergeFinalize(
    const TwoMergeFinalizeRequest& request) {
  TwoMergeFinalizeResponse response;
  std::vector<LabeledBlock> blocks;
  blocks.reserve(request.final_batch_size());
  for (const auto& blk_proto : request.final_batch()) {
    blocks.emplace_back(blk_proto.label(),
                        std::vector<uint8_t>(blk_proto.ciphertext().begin(),
                                             blk_proto.ciphertext().end()));
  }

  auto status = storage_->TwoMergeFinalize(std::move(blocks));
  if (!status.ok()) {
    response.set_success(false);
    response.set_error_message(std::string(status.message()));
    return response;
  }

  response.set_success(true);
  return response;
}

}  // namespace oram_okvs
