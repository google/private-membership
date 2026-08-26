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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LAT_OK_LATORAM_CLIENT_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LAT_OK_LATORAM_CLIENT_H_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/block_retrieval.h"
#include "research/oram_okvs/oram/lat/latoram_client.h"
#include "research/oram_okvs/oram/lat/latoram_server.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

// OKLat ORAM client (OKVS-based LatORAM without PIR - Algorithm 6).
class OkLatOramClient : public LatOramClient {
 public:
  OkLatOramClient(const LatOramConfig& config,
                  std::shared_ptr<LatOramServerInterface> server);

  absl::Status Init(
      absl::Span<const std::vector<uint8_t>> initial_data = {}) override;

  absl::StatusOr<std::vector<uint8_t>> QueryTargetBlock(
      BlockId block_id, size_t target_container) override;

  absl::Status AddLateralStash(
      size_t cp,
      const std::vector<std::pair<BlockId, std::vector<uint8_t>>>& stash_items)
      override;

  bool IsOkvs() const override { return true; }

 protected:
  std::unique_ptr<BlockRetrievalClient> retrieval_client_;
};

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_OK_LATORAM_CLIENT_H_
