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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LAT_LATORAM_STORAGE_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LAT_LATORAM_STORAGE_H_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/oblivious_container.h"

namespace oram_okvs {

// Unified concrete LatORAM server-side storage.

//
// Manages:
//  - Main storage container (KeyedObliviousContainer).
//  - Lateral stash containers (OkvsObliviousContainer or
//  KeyedObliviousContainer).
//  - Multi-round CacheShuffle matrix T and output container.
//  - Multi-round TwoMerge streaming and output container.
class LatOramStorage {
 public:
  LatOramStorage() = default;
  virtual ~LatOramStorage() = default;

  // 1. Basic container operations
  virtual absl::Status InitMain(std::vector<LabeledBlock> blocks);
  virtual absl::Status InitMain(std::unique_ptr<ObliviousContainer> container);

  virtual absl::Status AddLateralStash(std::vector<LabeledBlock> blocks);
  virtual absl::Status AddLateralStash(
      std::unique_ptr<ObliviousContainer> container);
  virtual absl::Status AddLateralStashOkvs(
      std::vector<std::vector<uint8_t>> okvs_table, uint64_t okvs_seed,
      uint32_t num_real, uint32_t okvs_backend,
      std::vector<LabeledBlock> fallback_blocks = {});

  // Evaluates query across active containers and returns the XOR sum.
  virtual absl::StatusOr<std::vector<uint8_t>> QueryAndXor(
      absl::Span<const std::string> labels);

  // Evaluates query on each active container and returns individual container
  // blocks.
  virtual absl::StatusOr<std::vector<std::vector<uint8_t>>> QueryContainers(
      absl::Span<const std::string> labels);

  virtual size_t MainSize() const;
  virtual size_t NumLateralStashes() const { return lateral_stashes_.size(); }
  virtual size_t LateralStashSize(size_t index) const;

  // 2. Server-side CacheShuffle operations
  virtual absl::StatusOr<size_t> InitCacheShuffle(
      absl::string_view source_container, size_t num_rows, size_t num_cols);

  virtual absl::StatusOr<std::vector<LabeledBlock>> CacheShuffleReadSource(
      size_t count);

  virtual absl::Status CacheShuffleWriteColumn(
      size_t col_idx, std::vector<LabeledBlock> column_blocks);

  virtual absl::StatusOr<std::vector<LabeledBlock>> CacheShuffleReadRow(
      size_t row_idx);

  virtual absl::Status CacheShuffleAppendOutput(
      std::vector<LabeledBlock> blocks);

  virtual size_t CacheShuffleOutputSize() const { return cs_output_.size(); }

  // 3. Server-side TwoMerge operations
  virtual absl::StatusOr<std::pair<size_t, size_t>> InitTwoMerge(
      absl::string_view container1, absl::string_view container2);

  virtual absl::StatusOr<
      std::pair<std::vector<LabeledBlock>, std::vector<LabeledBlock>>>
  TwoMergeReadChunks(size_t c1_count, size_t c2_count);

  virtual absl::Status TwoMergeAppendOutput(
      std::vector<LabeledBlock> merged_blocks);

  virtual absl::Status TwoMergeFinalize(std::vector<LabeledBlock> final_blocks);

 private:
  std::unique_ptr<ObliviousContainer> main_container_;
  std::vector<std::unique_ptr<ObliviousContainer>> lateral_stashes_;
  std::vector<std::vector<LabeledBlock>> lateral_stashes_raw_;

  // CacheShuffle state

  std::vector<LabeledBlock> cs_source_blocks_;
  size_t cs_source_read_idx_ = 0;
  size_t cs_num_rows_ = 0;
  size_t cs_num_cols_ = 0;
  std::vector<std::vector<LabeledBlock>> cs_matrix_;
  std::vector<LabeledBlock> cs_output_;

  // TwoMerge state
  std::vector<LabeledBlock> tm_c1_source_;
  std::vector<LabeledBlock> tm_c2_source_;
  size_t tm_c1_read_idx_ = 0;
  size_t tm_c2_read_idx_ = 0;
  std::vector<LabeledBlock> tm_output_;
};

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LATORAM_STORAGE_H_
