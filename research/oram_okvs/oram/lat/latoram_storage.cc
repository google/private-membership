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

#include "research/oram_okvs/oram/lat/latoram_storage.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/match.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "research/oram_okvs/okvs/polynomial_okvs.h"
#include "research/oram_okvs/okvs/random_band_okvs.h"
#include "research/oram_okvs/oram/oblivious_container.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

absl::Status LatOramStorage::InitMain(std::vector<LabeledBlock> blocks) {
  size_t num_real = blocks.size();
  size_t cipher_size = blocks.empty() ? 64 : blocks[0].ciphertext.size();

  absl::flat_hash_map<std::string, std::vector<uint8_t>> blocks_map;
  blocks_map.reserve(blocks.size());
  for (auto& [label, ct] : blocks) {
    blocks_map.emplace(std::move(label), std::move(ct));
  }

  main_container_ = std::make_unique<KeyedObliviousContainer>(
      std::move(blocks_map), num_real, cipher_size);

  lateral_stashes_.clear();
  lateral_stashes_raw_.clear();
  return absl::OkStatus();
}

absl::Status LatOramStorage::InitMain(
    std::unique_ptr<ObliviousContainer> container) {
  main_container_ = std::move(container);
  lateral_stashes_.clear();
  lateral_stashes_raw_.clear();
  return absl::OkStatus();
}

absl::Status LatOramStorage::AddLateralStash(std::vector<LabeledBlock> blocks) {
  size_t num_real = blocks.size();
  size_t cipher_size = blocks.empty() ? 64 : blocks[0].ciphertext.size();

  absl::flat_hash_map<std::string, std::vector<uint8_t>> stash_map;
  stash_map.reserve(blocks.size());
  for (auto& [label, ct] : blocks) {
    stash_map.emplace(label, ct);
  }

  lateral_stashes_raw_.push_back(std::move(blocks));
  lateral_stashes_.push_back(std::make_unique<KeyedObliviousContainer>(
      std::move(stash_map), num_real, cipher_size));
  return absl::OkStatus();
}

absl::Status LatOramStorage::AddLateralStash(
    std::unique_ptr<ObliviousContainer> container) {
  lateral_stashes_.push_back(std::move(container));
  lateral_stashes_raw_.emplace_back();
  return absl::OkStatus();
}

absl::Status LatOramStorage::AddLateralStashOkvs(
    std::vector<std::vector<uint8_t>> okvs_table, uint64_t okvs_seed,
    uint32_t num_real, uint32_t okvs_backend,
    std::vector<LabeledBlock> fallback_blocks) {
  size_t cipher_size = okvs_table.empty() ? 64 : okvs_table[0].size();

  if (okvs_backend == 0) {
    // Polynomial OKVS (Interpolation-OC)
    lateral_stashes_.push_back(
        std::make_unique<OkvsObliviousContainer<
            okvs::PolynomialOkvs<std::string, std::vector<uint8_t>>>>(
            std::move(okvs_table), okvs_seed, num_real, cipher_size));
  } else {
    // Random Band OKVS (RB-OC)
    lateral_stashes_.push_back(
        std::make_unique<OkvsObliviousContainer<
            okvs::RandomBandOkvs<std::string, std::vector<uint8_t>>>>(
            std::move(okvs_table), okvs_seed, num_real, cipher_size));
  }
  lateral_stashes_raw_.push_back(std::move(fallback_blocks));
  return absl::OkStatus();
}

absl::StatusOr<std::vector<std::vector<uint8_t>>>
LatOramStorage::QueryContainers(absl::Span<const std::string> labels) {
  if (labels.empty()) {
    return absl::InvalidArgumentError("Labels span cannot be empty");
  }

  std::vector<std::vector<uint8_t>> results;
  size_t total_containers = 1 + lateral_stashes_.size();
  results.reserve(total_containers);

  // 1. Query Main container (index 0)
  std::string main_label = labels[0];
  if (main_container_) {
    auto val_or = main_container_->ServerResponse(main_label);
    if (!val_or.ok()) return val_or.status();
    results.push_back(std::move(*val_or));
  } else {
    results.push_back(std::vector<uint8_t>(64, 0));
  }

  // 2. Query Lateral Stashes
  for (size_t s = 0; s < lateral_stashes_.size(); ++s) {
    std::string stash_label = (labels.size() > 1) ? labels[s + 1] : labels[0];
    if (lateral_stashes_[s]) {
      auto val_or = lateral_stashes_[s]->ServerResponse(stash_label);
      if (!val_or.ok()) return val_or.status();
      results.push_back(std::move(*val_or));
    } else {
      results.push_back(std::vector<uint8_t>(64, 0));
    }
  }

  return results;
}

absl::StatusOr<std::vector<uint8_t>> LatOramStorage::QueryAndXor(
    absl::Span<const std::string> labels) {
  auto containers_or = QueryContainers(labels);
  if (!containers_or.ok()) return containers_or.status();
  if (containers_or->empty()) {
    return std::vector<uint8_t>();
  }

  std::vector<uint8_t> xor_sum = (*containers_or)[0];
  for (size_t i = 1; i < containers_or->size(); ++i) {
    XorInPlace(absl::MakeSpan(xor_sum), (*containers_or)[i]);
  }
  return xor_sum;
}

size_t LatOramStorage::MainSize() const {
  if (main_container_) return main_container_->NumTotalBlocks();
  return 0;
}

size_t LatOramStorage::LateralStashSize(size_t index) const {
  if (index >= lateral_stashes_.size()) return 0;
  if (lateral_stashes_[index]) return lateral_stashes_[index]->NumRealBlocks();
  return 0;
}

// ============================================================================
// CacheShuffle & TwoMerge Helper
// ============================================================================

namespace {

absl::Status LoadContainerBlocks(
    absl::string_view name,
    const std::unique_ptr<ObliviousContainer>& main_container,
    const std::vector<std::unique_ptr<ObliviousContainer>>& lateral_stashes,
    const std::vector<std::vector<LabeledBlock>>& lateral_stashes_raw,
    const std::vector<LabeledBlock>& cs_output,
    std::vector<LabeledBlock>& dest) {
  if (name == "Main" || name == "main") {
    if (auto* keyed =
            dynamic_cast<KeyedObliviousContainer*>(main_container.get())) {
      dest.reserve(keyed->blocks().size());
      for (const auto& [l, c] : keyed->blocks()) dest.emplace_back(l, c);
    }
  } else if (name == "lateral_stashes" || name == "LateralStashes") {
    for (size_t s = 0; s < lateral_stashes.size(); ++s) {
      if (auto* keyed = dynamic_cast<KeyedObliviousContainer*>(
              lateral_stashes[s].get())) {
        for (const auto& [l, c] : keyed->blocks()) dest.emplace_back(l, c);
      } else if (s < lateral_stashes_raw.size()) {
        for (const auto& [l, c] : lateral_stashes_raw[s])
          dest.emplace_back(l, c);
      }
    }
  } else if (absl::StartsWith(name, "LatStash_") ||
             absl::StartsWith(name, "lat_stash_")) {
    int idx = std::stoi(std::string(name.substr(name.find('_') + 1))) - 1;
    if (idx < 0 || static_cast<size_t>(idx) >= lateral_stashes.size()) {
      return absl::InvalidArgumentError(
          absl::StrFormat("Invalid container: %s", name));
    }
    if (auto* keyed = dynamic_cast<KeyedObliviousContainer*>(
            lateral_stashes[idx].get())) {
      dest.reserve(keyed->blocks().size());
      for (const auto& [l, c] : keyed->blocks()) dest.emplace_back(l, c);
    } else if (static_cast<size_t>(idx) < lateral_stashes_raw.size()) {
      dest.reserve(lateral_stashes_raw[idx].size());
      for (const auto& [l, c] : lateral_stashes_raw[idx])
        dest.emplace_back(l, c);
    }
  } else if (name == "CacheShuffleOutput" || name == "cache_shuffle_output") {
    dest = cs_output;
  } else {
    return absl::InvalidArgumentError(
        absl::StrFormat("Unknown container: %s", name));
  }
  return absl::OkStatus();
}

}  // namespace

absl::StatusOr<size_t> LatOramStorage::InitCacheShuffle(
    absl::string_view source_container, size_t num_rows, size_t num_cols) {
  cs_source_blocks_.clear();
  cs_source_read_idx_ = 0;
  cs_num_rows_ = num_rows;
  cs_num_cols_ = num_cols;
  cs_matrix_.assign(num_cols, std::vector<LabeledBlock>());
  cs_output_.clear();

  auto status =
      LoadContainerBlocks(source_container, main_container_, lateral_stashes_,
                          lateral_stashes_raw_, cs_output_, cs_source_blocks_);
  if (!status.ok()) return status;

  return cs_source_blocks_.size();
}

absl::StatusOr<std::vector<LabeledBlock>>
LatOramStorage::CacheShuffleReadSource(size_t count) {
  std::vector<LabeledBlock> chunk;
  size_t available = cs_source_blocks_.size() - cs_source_read_idx_;
  size_t take = std::min(count, available);
  chunk.reserve(take);
  for (size_t i = 0; i < take; ++i) {
    chunk.push_back(std::move(cs_source_blocks_[cs_source_read_idx_++]));
  }
  return chunk;
}

absl::Status LatOramStorage::CacheShuffleWriteColumn(
    size_t col_idx, std::vector<LabeledBlock> column_blocks) {
  if (col_idx >= cs_matrix_.size()) {
    cs_matrix_.resize(col_idx + 1);
  }
  cs_matrix_[col_idx] = std::move(column_blocks);
  return absl::OkStatus();
}

absl::StatusOr<std::vector<LabeledBlock>> LatOramStorage::CacheShuffleReadRow(
    size_t row_idx) {
  std::vector<LabeledBlock> row_blocks;
  row_blocks.reserve(cs_matrix_.size());
  for (size_t c = 0; c < cs_matrix_.size(); ++c) {
    if (row_idx < cs_matrix_[c].size()) {
      row_blocks.push_back(cs_matrix_[c][row_idx]);
    }
  }
  return row_blocks;
}

absl::Status LatOramStorage::CacheShuffleAppendOutput(
    std::vector<LabeledBlock> blocks) {
  cs_output_.reserve(cs_output_.size() + blocks.size());
  for (auto& blk : blocks) {
    cs_output_.push_back(std::move(blk));
  }
  return absl::OkStatus();
}

// ============================================================================
// TwoMerge Operations
// ============================================================================

absl::StatusOr<std::pair<size_t, size_t>> LatOramStorage::InitTwoMerge(
    absl::string_view container1, absl::string_view container2) {
  tm_c1_source_.clear();
  tm_c2_source_.clear();
  tm_c1_read_idx_ = 0;
  tm_c2_read_idx_ = 0;
  tm_output_.clear();

  auto s1 =
      LoadContainerBlocks(container1, main_container_, lateral_stashes_,
                          lateral_stashes_raw_, cs_output_, tm_c1_source_);
  if (!s1.ok()) return s1;
  auto s2 =
      LoadContainerBlocks(container2, main_container_, lateral_stashes_,
                          lateral_stashes_raw_, cs_output_, tm_c2_source_);
  if (!s2.ok()) return s2;

  return std::make_pair(tm_c1_source_.size(), tm_c2_source_.size());
}

absl::StatusOr<std::pair<std::vector<LabeledBlock>, std::vector<LabeledBlock>>>
LatOramStorage::TwoMergeReadChunks(size_t c1_count, size_t c2_count) {
  std::vector<LabeledBlock> c1_chunk;
  size_t c1_avail = tm_c1_source_.size() - tm_c1_read_idx_;
  size_t c1_take = std::min(c1_count, c1_avail);
  c1_chunk.reserve(c1_take);
  for (size_t i = 0; i < c1_take; ++i) {
    c1_chunk.push_back(std::move(tm_c1_source_[tm_c1_read_idx_++]));
  }

  std::vector<LabeledBlock> c2_chunk;
  size_t c2_avail = tm_c2_source_.size() - tm_c2_read_idx_;
  size_t c2_take = std::min(c2_count, c2_avail);
  c2_chunk.reserve(c2_take);
  for (size_t i = 0; i < c2_take; ++i) {
    c2_chunk.push_back(std::move(tm_c2_source_[tm_c2_read_idx_++]));
  }

  return std::make_pair(std::move(c1_chunk), std::move(c2_chunk));
}

absl::Status LatOramStorage::TwoMergeAppendOutput(
    std::vector<LabeledBlock> merged_blocks) {
  tm_output_.reserve(tm_output_.size() + merged_blocks.size());
  for (auto& blk : merged_blocks) {
    tm_output_.push_back(std::move(blk));
  }
  return absl::OkStatus();
}

absl::Status LatOramStorage::TwoMergeFinalize(
    std::vector<LabeledBlock> final_blocks) {
  tm_output_.reserve(tm_output_.size() + final_blocks.size());
  for (auto& blk : final_blocks) {
    tm_output_.push_back(std::move(blk));
  }
  return InitMain(std::move(tm_output_));
}

}  // namespace oram_okvs
