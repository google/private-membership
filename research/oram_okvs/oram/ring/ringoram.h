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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_RING_RINGORAM_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_RING_RINGORAM_H_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/oram.h"
#include "research/oram_okvs/oram/ring/ringoram_client.h"
#include "research/oram_okvs/oram/ring/ringoram_server.h"
#include "research/oram_okvs/oram/ring/ringoram_storage.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

// Standalone in-process Ring ORAM instance combining Client and Server.
class RingOram : public Oram {
 public:
  explicit RingOram(const RingOramConfig& config,
                    std::unique_ptr<RingOramStorage> storage = nullptr);

  RingOram(const RingOramConfig& config,
           std::shared_ptr<RingOramServerInterface> server);

  ~RingOram() override = default;

  absl::Status Init() { return client_->Init(); }

  absl::StatusOr<std::vector<uint8_t>> Read(BlockId block_id) override {
    return client_->Read(block_id);
  }

  absl::Status Write(BlockId block_id,
                     absl::Span<const uint8_t> data) override {
    return client_->Write(block_id, data);
  }

  absl::StatusOr<std::vector<uint8_t>> Access(
      BlockId block_id, AccessOp op,
      absl::Span<const uint8_t> data = {}) override {
    return client_->Access(block_id, op, data);
  }

  size_t num_blocks() const override { return client_->num_blocks(); }
  size_t block_size() const override { return client_->block_size(); }
  OramStats GetStats() const override { return client_->GetStats(); }
  void ResetStats() override { client_->ResetStats(); }

  RingOramClient* client() { return client_.get(); }
  const RingOramClient* client() const { return client_.get(); }

 private:
  std::shared_ptr<RingOramServerInterface> server_;
  std::unique_ptr<RingOramClient> client_;
};

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_RINGORAM_H_
