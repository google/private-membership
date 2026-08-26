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

#include "research/oram_okvs/oram/lat/latoram_client.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/block_retrieval.h"
#include "research/oram_okvs/oram/crypto_utils.h"
#include "research/oram_okvs/oram/lat/latoram_server.h"
#include "research/oram_okvs/oram/lat/latoram_storage.h"
#include "research/oram_okvs/oram/lat/map_latoram_client.h"
#include "research/oram_okvs/oram/lat/multi_perm.h"
#include "research/oram_okvs/oram/lat/ok_latoram_client.h"
#include "research/oram_okvs/oram/lat/okp_latoram_client.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

namespace {

constexpr size_t kStashTarget = static_cast<size_t>(-1);

}  // namespace

LatOramClient::LatOramClient(const LatOramConfig& config,
                             std::shared_ptr<LatOramServerInterface> server)
    : config_(config),
      effective_m_(config.EffectiveStashSize()),
      effective_p_(config.EffectiveNumPhases()),
      ciphertext_block_size_(CiphertextBlockSize(config.block_size)),
      server_(std::move(server)),
      cp_(0),
      cq_(0),
      cd_(effective_p_ + 1, 1),
      initialized_(false) {}

std::unique_ptr<LatOramClient> LatOramClient::Create(
    const LatOramConfig& config,
    std::shared_ptr<LatOramServerInterface> server) {
  if (config.storage_type == LatOramStorageType::kOkLat) {
    return std::make_unique<OkLatOramClient>(config, std::move(server));
  } else if (config.storage_type == LatOramStorageType::kOkPLat) {
    return std::make_unique<OkPLatOramClient>(config, std::move(server));
  } else {
    return std::make_unique<MapLatOramClient>(config, std::move(server));
  }
}

absl::StatusOr<std::vector<uint8_t>> LatOramClient::Read(BlockId block_id) {
  return Access(block_id, AccessOp::kRead);
}

absl::Status LatOramClient::Write(BlockId block_id,
                                  absl::Span<const uint8_t> data) {
  auto result_or = Access(block_id, AccessOp::kWrite, data);
  if (!result_or.ok()) {
    return result_or.status();
  }
  return absl::OkStatus();
}

absl::StatusOr<std::vector<uint8_t>> LatOramClient::Access(
    BlockId block_id, AccessOp op, absl::Span<const uint8_t> data) {
  if (!initialized_) {
    return absl::FailedPreconditionError("LatOram is not initialized.");
  }
  if (block_id >= config_.num_blocks) {
    return absl::InvalidArgumentError(
        absl::StrFormat("Invalid block_id %lu: out of bounds [0, %zu)",
                        block_id, config_.num_blocks));
  }

  stats_.total_accesses++;
  if (op == AccessOp::kRead) {
    stats_.reads_count++;
  } else {
    stats_.writes_count++;
  }

  size_t target_container = 0;
  if (stash_.contains(block_id)) {
    target_container = kStashTarget;
  } else if (auto it = map_.find(block_id); it != map_.end()) {
    target_container = it->second;
    map_.erase(it);
  } else {
    target_container = 0;  // Main
  }

  std::vector<uint8_t> current_payload;
  if (target_container != kStashTarget) {
    auto target_payload_or = QueryTargetBlock(block_id, target_container);
    if (!target_payload_or.ok()) return target_payload_or.status();
    current_payload = std::move(*target_payload_or);
  } else {
    current_payload = stash_[block_id];
  }

  // Update data and client stash
  std::vector<uint8_t> return_payload;
  if (op == AccessOp::kRead) {
    return_payload = current_payload;
    stash_[block_id] = std::move(current_payload);
  } else if (op == AccessOp::kWrite) {
    std::vector<uint8_t> new_payload(config_.block_size, 0);
    size_t bytes_to_copy = std::min(data.size(), config_.block_size);
    if (bytes_to_copy > 0) {
      std::memcpy(new_payload.data(), data.data(), bytes_to_copy);
    }
    return_payload = new_payload;
    stash_[block_id] = std::move(new_payload);
  } else if (op == AccessOp::kReadAndWrite) {
    return_payload = current_payload;
    std::vector<uint8_t> new_payload(config_.block_size, 0);
    size_t bytes_to_copy = std::min(data.size(), config_.block_size);
    if (bytes_to_copy > 0) {
      std::memcpy(new_payload.data(), data.data(), bytes_to_copy);
    }
    stash_[block_id] = std::move(new_payload);
  }

  access_seq_++;
  if (block_id < block_versions_.size()) {
    block_versions_[block_id] = access_seq_;
  }

  // Postprocess
  auto postprocess_status = Postprocess();
  if (!postprocess_status.ok()) {
    return postprocess_status;
  }

  return return_payload;
}

absl::Status LatOramClient::Postprocess() {
  cq_++;
  if (cq_ < effective_m_) {
    return absl::OkStatus();
  }

  // Stash is full: create Lateral Stash via SmallShuffle
  cq_ = 0;
  cp_++;
  stats_.shuffles_count++;

  // Update map references for blocks in stash
  std::vector<std::pair<BlockId, std::vector<uint8_t>>> stash_items;
  stash_items.reserve(stash_.size());
  for (auto& [id, payload] : stash_) {
    map_[id] = cp_;
    stash_items.emplace_back(id, std::move(payload));
  }

  auto add_status = AddLateralStash(cp_, stash_items);
  if (!add_status.ok()) return add_status;

  stash_.clear();

  // If cp == P, trigger Main Reinitialization
  if (cp_ == effective_p_) {
    return MainReinit();
  }

  return absl::OkStatus();
}

absl::Status LatOramClient::CacheShuffle(absl::string_view source_container,
                                         absl::string_view old_key_e,
                                         absl::string_view old_key_l,
                                         absl::string_view new_key_e,
                                         absl::string_view new_key_l,
                                         absl::string_view key_d,
                                         size_t payload_size, uint64_t seed) {
  // 1. Query total count from server
  CacheShuffleInitRequest init_req;
  init_req.set_source_container(std::string(source_container));
  size_t init_req_sz = init_req.ByteSizeLong();
  stats_.bytes_written += init_req_sz;
  stats_.amortized_bytes_uploaded += init_req_sz;
  stats_.amortized_bytes_transferred += init_req_sz;

  auto init_resp_or = server_->HandleCacheShuffleInit(init_req);
  if (!init_resp_or.ok()) return init_resp_or.status();
  size_t init_resp_sz = init_resp_or->ByteSizeLong();
  stats_.bytes_read += init_resp_sz;
  stats_.amortized_bytes_downloaded += init_resp_sz;
  stats_.amortized_bytes_transferred += init_resp_sz;

  size_t total_elements = init_resp_or->total_blocks();
  if (total_elements == 0) {
    return absl::OkStatus();
  }

  size_t c = (effective_m_ > 0) ? effective_m_ : 4;
  size_t r = (total_elements + c - 1) / c;

  std::vector<size_t> freqs(r, c);
  if (r * c > total_elements) {
    size_t extra = r * c - total_elements;
    for (size_t i = 0; i < extra; ++i) {
      freqs[r - 1 - (i % r)]--;
    }
  }
  MultiPermDS multi_perm(freqs, seed);

  std::vector<std::vector<DecryptedBlock>> queues(r);
  std::string temp_key_e(old_key_e);
  std::string temp_key_l = GenerateRandomBytes(16);

  // Phase 1: Spray
  CacheShuffleSprayStepRequest spray_req;
  spray_req.set_col_idx(0);
  spray_req.set_download_count(static_cast<uint32_t>(r));

  size_t spray_req_0_sz = spray_req.ByteSizeLong();
  stats_.bytes_written += spray_req_0_sz;
  stats_.amortized_bytes_uploaded += spray_req_0_sz;
  stats_.amortized_bytes_transferred += spray_req_0_sz;

  auto spray_resp_or = server_->HandleCacheShuffleSprayStep(spray_req);
  if (!spray_resp_or.ok()) return spray_resp_or.status();
  size_t spray_resp_0_sz = spray_resp_or->ByteSizeLong();
  stats_.bytes_read += spray_resp_0_sz;
  stats_.amortized_bytes_downloaded += spray_resp_0_sz;
  stats_.amortized_bytes_transferred += spray_resp_0_sz;

  std::vector<LabeledBlock> downloaded_batch;
  for (const auto& p : spray_resp_or->downloaded_blocks()) {
    downloaded_batch.push_back(
        {p.label(),
         std::vector<uint8_t>(p.ciphertext().begin(), p.ciphertext().end())});
  }

  for (size_t col = 0; col < c; ++col) {
    // Decrypt downloaded blocks immediately and push to queues
    for (size_t b = 0; b < downloaded_batch.size(); ++b) {
      const auto& blk = downloaded_batch[b];
      auto dec_or =
          DecryptBlock(old_key_e, blk.label, blk.ciphertext, payload_size);
      auto z_or = multi_perm.Next();
      if (z_or.ok()) {
        if (dec_or.ok() && dec_or->block_id < config_.num_blocks) {
          queues[*z_or].push_back(std::move(*dec_or));
        }
      }
    }
    downloaded_batch.clear();

    // Pop 1 block from each queue Q_0..Q_{r-1} to form column
    spray_req.Clear();
    spray_req.set_col_idx(static_cast<uint32_t>(col));
    for (size_t k = 0; k < r; ++k) {
      auto* p = spray_req.add_uploaded_column();
      if (!queues[k].empty()) {
        DecryptedBlock item = std::move(queues[k].front());
        queues[k].erase(queues[k].begin());

        std::string temp_label = GenerateBlockLabel(temp_key_l, item.block_id);
        std::vector<uint8_t> temp_blk =
            EncryptBlock(temp_key_e, temp_label, item.block_id, item.payload,
                         payload_size, item.version);
        p->set_label(std::move(temp_label));
        p->set_ciphertext(temp_blk.data(), temp_blk.size());
      } else {
        std::string temp_label =
            GenerateDummyLabel(temp_key_l, col * r + k + 1);
        std::vector<uint8_t> temp_blk = GenerateDummyBlock(
            key_d, temp_label, CiphertextBlockSize(payload_size));
        p->set_label(std::move(temp_label));
        p->set_ciphertext(temp_blk.data(), temp_blk.size());
      }
    }

    if (col + 1 < c) {
      spray_req.set_download_count(static_cast<uint32_t>(r));
    } else {
      spray_req.set_download_count(0);
    }

    size_t spray_req_sz = spray_req.ByteSizeLong();
    stats_.bytes_written += spray_req_sz;
    stats_.amortized_bytes_uploaded += spray_req_sz;
    stats_.amortized_bytes_transferred += spray_req_sz;

    spray_resp_or = server_->HandleCacheShuffleSprayStep(spray_req);
    if (!spray_resp_or.ok()) return spray_resp_or.status();
    size_t spray_resp_sz = spray_resp_or->ByteSizeLong();
    stats_.bytes_read += spray_resp_sz;
    stats_.amortized_bytes_downloaded += spray_resp_sz;
    stats_.amortized_bytes_transferred += spray_resp_sz;

    for (const auto& p : spray_resp_or->downloaded_blocks()) {
      downloaded_batch.push_back(
          {p.label(),
           std::vector<uint8_t>(p.ciphertext().begin(), p.ciphertext().end())});
    }
  }

  // Phase 2: Recalibrate
  std::vector<LabeledBlockProto> recalibrated_upload;

  for (size_t row = 0; row < r; ++row) {
    CacheShuffleRecalibrateStepRequest recal_req;
    recal_req.set_row_idx(static_cast<uint32_t>(row));
    for (auto& proto_blk : recalibrated_upload) {
      *recal_req.add_recalibrated_blocks() = std::move(proto_blk);
    }
    recalibrated_upload.clear();

    size_t recal_req_sz = recal_req.ByteSizeLong();
    stats_.bytes_written += recal_req_sz;
    stats_.amortized_bytes_uploaded += recal_req_sz;
    stats_.amortized_bytes_transferred += recal_req_sz;

    auto recal_resp_or = server_->HandleCacheShuffleRecalibrateStep(recal_req);
    if (!recal_resp_or.ok()) return recal_resp_or.status();
    size_t recal_resp_sz = recal_resp_or->ByteSizeLong();
    stats_.bytes_read += recal_resp_sz;
    stats_.amortized_bytes_downloaded += recal_resp_sz;
    stats_.amortized_bytes_transferred += recal_resp_sz;

    std::vector<DecryptedBlock> valid_blocks;
    for (const auto& p : recal_resp_or->row_blocks()) {
      auto dec_or = DecryptBlock(
          temp_key_e, p.label(),
          absl::Span<const uint8_t>(
              reinterpret_cast<const uint8_t*>(p.ciphertext().data()),
              p.ciphertext().size()),
          payload_size);
      if (dec_or.ok()) {
        valid_blocks.push_back(std::move(*dec_or));
      }
    }

    // Append remaining blocks from queue Q_row
    for (auto& item : queues[row]) {
      valid_blocks.push_back(std::move(item));
    }
    queues[row].clear();

    // Uniformly random shuffle valid blocks locally in O(c) memory
    std::vector<size_t> local_perm =
        GenerateRandomPermutation(valid_blocks.size(), seed + row);
    for (size_t idx : local_perm) {
      const auto& item = valid_blocks[idx];
      std::string new_label = GenerateBlockLabel(new_key_l, item.block_id);
      std::vector<uint8_t> new_blk =
          EncryptBlock(new_key_e, new_label, item.block_id, item.payload,
                       payload_size, item.version);
      LabeledBlockProto p;
      p.set_label(std::move(new_label));
      p.set_ciphertext(new_blk.data(), new_blk.size());
      recalibrated_upload.push_back(std::move(p));
    }
  }

  // Upload final row to finish CacheShuffle
  if (!recalibrated_upload.empty()) {
    CacheShuffleRecalibrateStepRequest recal_req;
    recal_req.set_row_idx(static_cast<uint32_t>(r));
    for (auto& proto_blk : recalibrated_upload) {
      *recal_req.add_recalibrated_blocks() = std::move(proto_blk);
    }
    size_t final_recal_sz = recal_req.ByteSizeLong();
    stats_.bytes_written += final_recal_sz;
    stats_.amortized_bytes_uploaded += final_recal_sz;
    stats_.amortized_bytes_transferred += final_recal_sz;

    auto recal_resp_or = server_->HandleCacheShuffleRecalibrateStep(recal_req);
    if (!recal_resp_or.ok()) return recal_resp_or.status();
    size_t final_resp_sz = recal_resp_or->ByteSizeLong();
    stats_.bytes_read += final_resp_sz;
    stats_.amortized_bytes_downloaded += final_resp_sz;
    stats_.amortized_bytes_transferred += final_resp_sz;
  }

  return absl::OkStatus();
}

absl::Status LatOramClient::TwoMerge(
    absl::string_view container1, absl::string_view container2,
    absl::string_view old_key_e, absl::string_view old_key_l_c1,
    absl::string_view old_key_l_c2, absl::string_view new_key_e,
    absl::string_view new_key_l, absl::string_view key_d, size_t payload_size,
    size_t chunk_size, uint64_t seed) {
  TwoMergeInitRequest init_req;
  init_req.set_container1(std::string(container1));
  init_req.set_container2(std::string(container2));
  init_req.set_chunk_size(static_cast<uint32_t>(chunk_size));
  size_t init_req_sz = init_req.ByteSizeLong();
  stats_.bytes_written += init_req_sz;
  stats_.amortized_bytes_uploaded += init_req_sz;
  stats_.amortized_bytes_transferred += init_req_sz;

  auto init_resp_or = server_->HandleTwoMergeInit(init_req);
  if (!init_resp_or.ok()) return init_resp_or.status();
  size_t init_resp_sz = init_resp_or->ByteSizeLong();
  stats_.bytes_read += init_resp_sz;
  stats_.amortized_bytes_downloaded += init_resp_sz;
  stats_.amortized_bytes_transferred += init_resp_sz;

  size_t size1 = init_resp_or->size1();
  size_t size2 = init_resp_or->size2();
  size_t remaining_c1 = size1;
  size_t remaining_c2 = size2;

  std::vector<std::vector<uint8_t>> final_real_payloads(
      config_.num_blocks, std::vector<uint8_t>(config_.block_size, 0));
  std::vector<uint64_t> merged_versions(config_.num_blocks, 0);

  // Stream blocks from C1 and C2 in chunks of size chunk_size
  while (remaining_c1 > 0 || remaining_c2 > 0) {
    size_t fetch1 = std::min(chunk_size, remaining_c1);
    size_t fetch2 = std::min(chunk_size, remaining_c2);
    remaining_c1 -= fetch1;
    remaining_c2 -= fetch2;

    TwoMergeStepRequest step_req;
    step_req.set_fetch_c1_count(static_cast<uint32_t>(fetch1));
    step_req.set_fetch_c2_count(static_cast<uint32_t>(fetch2));

    size_t step_req_sz = step_req.ByteSizeLong();
    stats_.bytes_written += step_req_sz;
    stats_.amortized_bytes_uploaded += step_req_sz;
    stats_.amortized_bytes_transferred += step_req_sz;

    auto step_resp_or = server_->HandleTwoMergeStep(step_req);
    if (!step_resp_or.ok()) return step_resp_or.status();
    size_t step_resp_sz = step_resp_or->ByteSizeLong();
    stats_.bytes_read += step_resp_sz;
    stats_.amortized_bytes_downloaded += step_resp_sz;
    stats_.amortized_bytes_transferred += step_resp_sz;

    for (const auto& p : step_resp_or->c1_blocks()) {
      auto dec_or = DecryptBlock(
          old_key_e, p.label(),
          absl::Span<const uint8_t>(
              reinterpret_cast<const uint8_t*>(p.ciphertext().data()),
              p.ciphertext().size()),
          payload_size);
      if (dec_or.ok() && dec_or->block_id < config_.num_blocks) {
        if (dec_or->version >= merged_versions[dec_or->block_id]) {
          merged_versions[dec_or->block_id] = dec_or->version;
          final_real_payloads[dec_or->block_id] = std::move(dec_or->payload);
        }
      }
    }
    for (const auto& p : step_resp_or->c2_blocks()) {
      auto dec_or = DecryptBlock(
          old_key_e, p.label(),
          absl::Span<const uint8_t>(
              reinterpret_cast<const uint8_t*>(p.ciphertext().data()),
              p.ciphertext().size()),
          payload_size);
      if (dec_or.ok() && dec_or->block_id < config_.num_blocks) {
        if (dec_or->version >= merged_versions[dec_or->block_id]) {
          merged_versions[dec_or->block_id] = dec_or->version;
          final_real_payloads[dec_or->block_id] = std::move(dec_or->payload);
        }
      }
    }
  }

  // Stream merged and newly encrypted Main container back to server in chunks
  // of size M
  size_t total_dummies = effective_p_ * effective_m_;
  size_t total_main_size = config_.num_blocks + total_dummies;
  std::vector<size_t> final_perm =
      GenerateRandomPermutation(total_main_size, seed);

  size_t uploaded_count = 0;
  while (uploaded_count < total_main_size) {
    size_t batch_size =
        std::min(effective_m_, total_main_size - uploaded_count);
    std::vector<LabeledBlockProto> batch_protos;
    batch_protos.reserve(batch_size);

    for (size_t b = 0; b < batch_size; ++b) {
      size_t item_idx = final_perm[uploaded_count + b];
      LabeledBlockProto p;
      if (item_idx < config_.num_blocks) {
        BlockId b_id = static_cast<BlockId>(item_idx);
        std::string label = GenerateBlockLabel(new_key_l, b_id);
        std::vector<uint8_t> ct =
            EncryptBlock(new_key_e, label, b_id, final_real_payloads[b_id],
                         payload_size, merged_versions[b_id]);
        p.set_label(std::move(label));
        p.set_ciphertext(ct.data(), ct.size());
      } else {
        uint64_t dummy_idx =
            static_cast<uint64_t>(item_idx - config_.num_blocks + 1);
        std::string label = GenerateDummyLabel(new_key_l, dummy_idx);
        std::vector<uint8_t> ct =
            GenerateDummyBlock(key_d, label, CiphertextBlockSize(payload_size));
        p.set_label(std::move(label));
        p.set_ciphertext(ct.data(), ct.size());
      }
      batch_protos.push_back(std::move(p));
    }
    uploaded_count += batch_size;

    if (uploaded_count < total_main_size) {
      TwoMergeStepRequest step_req;
      for (auto& proto_blk : batch_protos) {
        *step_req.add_merged_batch() = std::move(proto_blk);
      }
      size_t step_batch_sz = step_req.ByteSizeLong();
      stats_.bytes_written += step_batch_sz;
      stats_.amortized_bytes_uploaded += step_batch_sz;
      stats_.amortized_bytes_transferred += step_batch_sz;

      auto step_resp_or = server_->HandleTwoMergeStep(step_req);
      if (!step_resp_or.ok()) return step_resp_or.status();
      size_t step_resp_sz = step_resp_or->ByteSizeLong();
      stats_.bytes_read += step_resp_sz;
      stats_.amortized_bytes_downloaded += step_resp_sz;
      stats_.amortized_bytes_transferred += step_resp_sz;
    } else {
      TwoMergeFinalizeRequest fin_req;
      for (auto& proto_blk : batch_protos) {
        *fin_req.add_final_batch() = std::move(proto_blk);
      }
      size_t fin_req_sz = fin_req.ByteSizeLong();
      stats_.bytes_written += fin_req_sz;
      stats_.amortized_bytes_uploaded += fin_req_sz;
      stats_.amortized_bytes_transferred += fin_req_sz;

      auto fin_resp_or = server_->HandleTwoMergeFinalize(fin_req);
      if (!fin_resp_or.ok()) return fin_resp_or.status();
      size_t fin_resp_sz = fin_resp_or->ByteSizeLong();
      stats_.bytes_read += fin_resp_sz;
      stats_.amortized_bytes_downloaded += fin_resp_sz;
      stats_.amortized_bytes_transferred += fin_resp_sz;

      if (!fin_resp_or->success()) {
        return absl::InternalError(
            absl::StrFormat("Server failed TwoMergeFinalize: %s",
                            fin_resp_or->error_message()));
      }
    }
  }

  return absl::OkStatus();
}

absl::Status LatOramClient::MainReinit() {
  stats_.main_reinits_count++;

  std::string new_master_key_l = GenerateRandomBytes(16);
  std::string temp_key_l = GenerateRandomBytes(16);

  // 1. Joint CacheShuffle on server's lateral stashes
  auto shuffle_status =
      CacheShuffle("lateral_stashes", key_e_, key_l_, key_e_, temp_key_l,
                   key_d_, config_.block_size, config_.seed);
  if (!shuffle_status.ok()) return shuffle_status;

  // 2. Joint TwoMerge on server's "main" and "cache_shuffle_output"
  std::string new_key_l_0 = DeriveSubKey(new_master_key_l, uint64_t{0});
  std::string old_key_l_0 = DeriveSubKey(key_l_, uint64_t{0});

  auto merge_status = TwoMerge(
      "main", "cache_shuffle_output", key_e_, old_key_l_0, temp_key_l, key_e_,
      new_key_l_0, key_d_, config_.block_size, effective_m_, config_.seed);
  if (!merge_status.ok()) return merge_status;

  // 3. Reset client state
  key_l_ = new_master_key_l;
  map_.clear();
  stash_.clear();
  cp_ = 0;
  cq_ = 0;
  cd_.assign(effective_p_ + 1, 1);

  return absl::OkStatus();
}

std::vector<LabeledBlock> LatOramClient::SmallShuffle(
    absl::Span<const std::pair<BlockId, std::vector<uint8_t>>> stash_blocks,
    size_t num_dummies, absl::string_view key_e, absl::string_view key_l,
    absl::string_view key_d, size_t payload_size, uint64_t seed,
    absl::Span<const uint64_t> block_versions) {
  size_t num_real = stash_blocks.size();
  size_t total = num_real + num_dummies;
  std::vector<LabeledBlock> container(total);

  std::vector<size_t> perm = GenerateRandomPermutation(total, seed);

  size_t real_cursor = 0;
  uint64_t dummy_cursor = 1;
  size_t cipher_size = CiphertextBlockSize(payload_size);

  for (size_t i = 0; i < total; ++i) {
    if (perm[i] < num_real) {
      const auto& [block_id, payload] = stash_blocks[real_cursor++];
      std::string label = GenerateBlockLabel(key_l, block_id);
      uint64_t ver =
          (block_id < block_versions.size()) ? block_versions[block_id] : 0;
      std::vector<uint8_t> blk =
          EncryptBlock(key_e, label, block_id, payload, payload_size, ver);
      container[i] = {std::move(label), std::move(blk)};
    } else {
      uint64_t dummy_idx = dummy_cursor++;
      std::string label = GenerateDummyLabel(key_l, dummy_idx);
      std::vector<uint8_t> blk = GenerateDummyBlock(key_d, label, cipher_size);
      container[i] = {std::move(label), std::move(blk)};
    }
  }

  return container;
}

}  // namespace oram_okvs
