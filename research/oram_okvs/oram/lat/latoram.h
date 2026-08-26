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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LAT_LATORAM_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LAT_LATORAM_H_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/lat/latoram_client.h"
#include "research/oram_okvs/oram/lat/latoram_server.h"
#include "research/oram_okvs/oram/lat/latoram_storage.h"
#include "research/oram_okvs/oram/oram.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

// LatOram: High-level convenient wrapper combining LatOramClient and
// LatOramServer.
class LatOram : public Oram {
 public:
  explicit LatOram(const LatOramConfig& config,
                   std::unique_ptr<LatOramStorage> storage = nullptr);

  LatOram(const LatOramConfig& config,
          std::shared_ptr<LatOramServerInterface> server);

  ~LatOram() override = default;

  absl::Status Init(absl::Span<const std::vector<uint8_t>> initial_data = {}) {
    return client_->Init(initial_data);
  }

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

  LatOramClient* client() { return client_.get(); }
  const LatOramClient* client() const { return client_.get(); }
  LatOramServerInterface* server() { return server_.get(); }

  static std::vector<LabeledBlock> SmallShuffle(
      absl::Span<const std::pair<BlockId, std::vector<uint8_t>>> stash_blocks,
      size_t num_dummies, absl::string_view key_e, absl::string_view key_l,
      absl::string_view key_d, size_t payload_size, uint64_t seed = 0,
      absl::Span<const uint64_t> block_versions = {}) {
    return LatOramClient::SmallShuffle(stash_blocks, num_dummies, key_e, key_l,
                                       key_d, payload_size, seed,
                                       block_versions);
  }

 private:
  std::shared_ptr<LatOramServerInterface> server_;
  std::unique_ptr<LatOramClient> client_;
};

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LATORAM_H_
