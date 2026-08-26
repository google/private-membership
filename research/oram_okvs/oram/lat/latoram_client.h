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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LAT_LATORAM_CLIENT_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LAT_LATORAM_CLIENT_H_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/block_retrieval.h"
#include "research/oram_okvs/oram/crypto_utils.h"
#include "research/oram_okvs/oram/lat/latoram_server.h"
#include "research/oram_okvs/oram/lat/latoram_storage.h"
#include "research/oram_okvs/oram/oram.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

// Base class for LatORAM client implementations.
//
// Encapsulates the phase progression, client position map ($map$), stash ($S$),
// and streaming oblivious shuffles (CacheShuffle and TwoMerge) executed jointly
// with the server.
class LatOramClient : public Oram {
 public:
  LatOramClient(const LatOramConfig& config,
                std::shared_ptr<LatOramServerInterface> server);

  ~LatOramClient() override = default;

  // Factory method returning the appropriate concrete LatORAM client.
  static std::unique_ptr<LatOramClient> Create(
      const LatOramConfig& config,
      std::shared_ptr<LatOramServerInterface> server);

  // Initializes the ORAM storage on the server with optional initial data.
  virtual absl::Status Init(
      absl::Span<const std::vector<uint8_t>> initial_data = {}) = 0;

  // Queries the target block from active containers (Main + Lateral Stashes).
  virtual absl::StatusOr<std::vector<uint8_t>> QueryTargetBlock(
      BlockId block_id, size_t target_container) = 0;

  // Uploads stash items to create a new lateral stash for phase cp.
  virtual absl::Status AddLateralStash(
      size_t cp, const std::vector<std::pair<BlockId, std::vector<uint8_t>>>&
                     stash_items) = 0;

  virtual bool IsOkvs() const = 0;

  // Oram interface implementation.
  absl::StatusOr<std::vector<uint8_t>> Read(BlockId block_id) override;
  absl::Status Write(BlockId block_id, absl::Span<const uint8_t> data) override;
  absl::StatusOr<std::vector<uint8_t>> Access(
      BlockId block_id, AccessOp op,
      absl::Span<const uint8_t> data = {}) override;

  size_t num_blocks() const override { return config_.num_blocks; }
  size_t block_size() const override { return config_.block_size; }
  OramStats GetStats() const override { return stats_; }
  void ResetStats() override { stats_.Reset(); }

  // Protocol state accessors.
  size_t current_phase() const { return cp_; }
  size_t current_query_in_phase() const { return cq_; }
  size_t stash_size() const { return stash_.size(); }
  size_t position_map_size() const { return map_.size(); }
  const LatOramConfig& config() const { return config_; }

  // Subroutines executed jointly with the server.
  static std::vector<LabeledBlock> SmallShuffle(
      absl::Span<const std::pair<BlockId, std::vector<uint8_t>>> stash_blocks,
      size_t num_dummies, absl::string_view key_e, absl::string_view key_l,
      absl::string_view key_d, size_t payload_size, uint64_t seed = 0,
      absl::Span<const uint64_t> block_versions = {});

  absl::Status CacheShuffle(absl::string_view source_container,
                            absl::string_view old_key_e,
                            absl::string_view old_key_l,
                            absl::string_view new_key_e,
                            absl::string_view new_key_l,
                            absl::string_view key_d, size_t payload_size,
                            uint64_t seed = 0);

  absl::Status TwoMerge(
      absl::string_view container1, absl::string_view container2,
      absl::string_view old_key_e, absl::string_view old_key_l_c1,
      absl::string_view old_key_l_c2, absl::string_view new_key_e,
      absl::string_view new_key_l, absl::string_view key_d, size_t payload_size,
      size_t chunk_size, uint64_t seed = 0);

 protected:
  absl::Status Postprocess();
  absl::Status MainReinit();

  LatOramConfig config_;
  size_t effective_m_;
  size_t effective_p_;
  size_t ciphertext_block_size_;

  std::shared_ptr<LatOramServerInterface> server_;

  // Master cryptographic keys.
  std::string key_e_;
  std::string key_l_;
  std::string key_d_;

  // Client state.
  absl::flat_hash_map<BlockId, size_t> map_;
  absl::flat_hash_map<BlockId, std::vector<uint8_t>> stash_;
  std::vector<uint64_t> block_versions_;
  uint64_t access_seq_ = 0;
  size_t cp_;
  size_t cq_;
  std::vector<uint64_t> cd_;

  OramStats stats_;
  bool initialized_ = false;
};

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LATORAM_CLIENT_H_
