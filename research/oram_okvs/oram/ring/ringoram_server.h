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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_RING_RINGORAM_SERVER_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_RING_RINGORAM_SERVER_H_

#include <cstddef>
#include <memory>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "research/oram_okvs/oram/block_retrieval.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/ring/ringoram_storage.h"

namespace oram_okvs {

// Abstract interface for Ring ORAM Server RPC handlers.
class RingOramServerInterface {
 public:
  virtual ~RingOramServerInterface() = default;

  virtual absl::StatusOr<RingOramInitTreeResponse> HandleInitTree(
      const RingOramInitTreeRequest& request) = 0;

  virtual absl::StatusOr<RingOramReadPathMetadataResponse>
  HandleReadPathMetadata(const RingOramReadPathMetadataRequest& request) = 0;

  virtual absl::StatusOr<RingOramReadPathBlocksResponse> HandleReadPathBlocks(
      const RingOramReadPathBlocksRequest& request) = 0;

  virtual absl::StatusOr<RingOramReadPathFullResponse> HandleReadPathFull(
      const RingOramReadPathFullRequest& request) = 0;

  virtual absl::StatusOr<RingOramWritePathFullResponse> HandleWritePathFull(
      const RingOramWritePathFullRequest& request) = 0;

  virtual absl::StatusOr<RingOramReshuffleBucketResponse> HandleReshuffleBucket(
      const RingOramReshuffleBucketRequest& request) = 0;

  virtual double total_server_time_ms() const { return 0.0; }
  virtual void ResetStats() {}
};

// Concrete in-process Ring ORAM server implementation.
class RingOramServer : public RingOramServerInterface {
 public:
  explicit RingOramServer(std::unique_ptr<RingOramStorage> storage = nullptr);
  ~RingOramServer() override = default;

  absl::StatusOr<RingOramInitTreeResponse> HandleInitTree(
      const RingOramInitTreeRequest& request) override;

  absl::StatusOr<RingOramReadPathMetadataResponse> HandleReadPathMetadata(
      const RingOramReadPathMetadataRequest& request) override;

  absl::StatusOr<RingOramReadPathBlocksResponse> HandleReadPathBlocks(
      const RingOramReadPathBlocksRequest& request) override;

  absl::StatusOr<RingOramReadPathFullResponse> HandleReadPathFull(
      const RingOramReadPathFullRequest& request) override;

  absl::StatusOr<RingOramWritePathFullResponse> HandleWritePathFull(
      const RingOramWritePathFullRequest& request) override;

  absl::StatusOr<RingOramReshuffleBucketResponse> HandleReshuffleBucket(
      const RingOramReshuffleBucketRequest& request) override;

  RingOramStorage* storage() { return storage_.get(); }
  size_t total_accesses_served() const { return total_accesses_served_; }
  size_t total_evictions_served() const { return total_evictions_served_; }
  double total_server_time_ms() const override { return total_server_time_ms_; }
  void ResetStats() override {
    total_accesses_served_ = 0;
    total_evictions_served_ = 0;
    total_server_time_ms_ = 0.0;
  }

 private:
  std::unique_ptr<RingOramStorage> storage_;
  std::unique_ptr<BlockRetrievalServer> retrieval_server_;
  size_t total_accesses_served_ = 0;
  size_t total_evictions_served_ = 0;
  double total_server_time_ms_ = 0.0;
};

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_RING_RINGORAM_SERVER_H_
