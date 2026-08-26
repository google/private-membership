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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LAT_LATORAM_SERVER_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LAT_LATORAM_SERVER_H_

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "research/oram_okvs/oram/block_retrieval.h"
#include "research/oram_okvs/oram/lat/latoram_storage.h"
#include "research/oram_okvs/oram/oram.pb.h"

namespace oram_okvs {

// Abstract interface for LatORAM Server RPC handlers.
class LatOramServerInterface {
 public:
  virtual ~LatOramServerInterface() = default;

  // 1. Basic container RPCs
  virtual absl::StatusOr<LatOramInitMainResponse> HandleInitMain(
      const LatOramInitMainRequest& request) = 0;

  virtual absl::StatusOr<LatOramAddLateralStashResponse> HandleAddLateralStash(
      const LatOramAddLateralStashRequest& request) = 0;

  virtual absl::StatusOr<LatOramQueryResponse> HandleQuery(
      const LatOramQueryRequest& request) = 0;

  // 2. Multi-round CacheShuffle RPCs
  virtual absl::StatusOr<CacheShuffleInitResponse> HandleCacheShuffleInit(
      const CacheShuffleInitRequest& request) = 0;

  virtual absl::StatusOr<CacheShuffleSprayStepResponse>
  HandleCacheShuffleSprayStep(const CacheShuffleSprayStepRequest& request) = 0;

  virtual absl::StatusOr<CacheShuffleRecalibrateStepResponse>
  HandleCacheShuffleRecalibrateStep(
      const CacheShuffleRecalibrateStepRequest& request) = 0;

  // 3. Multi-round TwoMerge RPCs
  virtual absl::StatusOr<TwoMergeInitResponse> HandleTwoMergeInit(
      const TwoMergeInitRequest& request) = 0;

  virtual absl::StatusOr<TwoMergeStepResponse> HandleTwoMergeStep(
      const TwoMergeStepRequest& request) = 0;

  virtual absl::StatusOr<TwoMergeFinalizeResponse> HandleTwoMergeFinalize(
      const TwoMergeFinalizeRequest& request) = 0;

  virtual double total_server_time_ms() const { return 0.0; }
  virtual void ResetStats() {}
};

// LatORAM Server implementation.
// Manages untrusted server-side container storage, executes online XOR queries,
// and manages multi-round server-side CacheShuffle and TwoMerge sessions.
class LatOramServer : public LatOramServerInterface {
 public:
  explicit LatOramServer(std::unique_ptr<LatOramStorage> storage = nullptr);

  ~LatOramServer() override = default;

  absl::StatusOr<LatOramInitMainResponse> HandleInitMain(
      const LatOramInitMainRequest& request) override;

  absl::StatusOr<LatOramAddLateralStashResponse> HandleAddLateralStash(
      const LatOramAddLateralStashRequest& request) override;

  absl::StatusOr<LatOramQueryResponse> HandleQuery(
      const LatOramQueryRequest& request) override;

  absl::StatusOr<CacheShuffleInitResponse> HandleCacheShuffleInit(
      const CacheShuffleInitRequest& request) override;

  absl::StatusOr<CacheShuffleSprayStepResponse> HandleCacheShuffleSprayStep(
      const CacheShuffleSprayStepRequest& request) override;

  absl::StatusOr<CacheShuffleRecalibrateStepResponse>
  HandleCacheShuffleRecalibrateStep(
      const CacheShuffleRecalibrateStepRequest& request) override;

  absl::StatusOr<TwoMergeInitResponse> HandleTwoMergeInit(
      const TwoMergeInitRequest& request) override;

  absl::StatusOr<TwoMergeStepResponse> HandleTwoMergeStep(
      const TwoMergeStepRequest& request) override;

  absl::StatusOr<TwoMergeFinalizeResponse> HandleTwoMergeFinalize(
      const TwoMergeFinalizeRequest& request) override;

  // Server state accessors.
  size_t MainSize() const { return storage_->MainSize(); }
  size_t NumLateralStashes() const { return storage_->NumLateralStashes(); }
  size_t total_queries_served() const { return total_queries_served_; }
  size_t total_xor_ops() const { return total_xor_ops_; }
  size_t total_shuffles_served() const { return total_shuffles_served_; }
  double total_server_time_ms() const override { return total_server_time_ms_; }
  void ResetStats() override {
    total_queries_served_ = 0;
    total_xor_ops_ = 0;
    total_shuffles_served_ = 0;
    total_server_time_ms_ = 0.0;
  }

 private:
  std::unique_ptr<LatOramStorage> storage_;
  std::unique_ptr<BlockRetrievalServer> retrieval_server_;
  size_t total_queries_served_ = 0;

  size_t total_xor_ops_ = 0;
  size_t total_shuffles_served_ = 0;
  double total_server_time_ms_ = 0.0;
};

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LAT_LATORAM_SERVER_H_
