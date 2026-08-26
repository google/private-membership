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

#include "research/oram_okvs/oram/lat/map_latoram_client.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/block_retrieval.h"
#include "research/oram_okvs/oram/crypto_utils.h"
#include "research/oram_okvs/oram/lat/latoram_client.h"
#include "research/oram_okvs/oram/lat/latoram_server.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

MapLatOramClient::MapLatOramClient(
    const LatOramConfig& config, std::shared_ptr<LatOramServerInterface> server)
    : LatOramClient(config, std::move(server)) {
  retrieval_client_ = BlockRetrievalClient::CreateXor();
}

absl::Status MapLatOramClient::Init(
    absl::Span<const std::vector<uint8_t>> initial_data) {
  if (config_.num_blocks == 0) {
    return absl::InvalidArgumentError("num_blocks must be greater than 0.");
  }
  if (config_.block_size == 0) {
    return absl::InvalidArgumentError("block_size must be greater than 0.");
  }
  if (!server_) {
    return absl::FailedPreconditionError("Server interface is not set.");
  }

  // 1. Generate master cryptographic keys
  key_e_ = GenerateRandomBytes(16);
  key_l_ = GenerateRandomBytes(16);
  key_d_ = GenerateRandomBytes(16);

  // 2. Initialize client state
  map_.clear();
  stash_.clear();
  block_versions_.assign(config_.num_blocks, 0);
  access_seq_ = 0;
  cp_ = 0;
  cq_ = 0;
  cd_.assign(effective_p_ + 1, 1);
  stats_.Reset();

  // 3. Prepare real and dummy blocks for Main container
  size_t total_dummies = effective_p_ * effective_m_;
  size_t total_blocks = config_.num_blocks + total_dummies;

  std::string key_l_0 = DeriveSubKey(key_l_, uint64_t{0});
  std::vector<size_t> perm =
      GenerateRandomPermutation(total_blocks, config_.seed);

  LatOramInitMainRequest request;

  for (size_t i = 0; i < total_blocks; ++i) {
    size_t item_idx = perm[i];
    auto* blk_proto = request.add_blocks();
    if (item_idx < config_.num_blocks) {
      BlockId block_id = static_cast<BlockId>(item_idx);
      absl::Span<const uint8_t> payload;
      if (block_id < initial_data.size() && !initial_data[block_id].empty()) {
        payload = absl::MakeSpan(initial_data[block_id]);
      }
      std::string label = GenerateBlockLabel(key_l_0, block_id);
      std::vector<uint8_t> blk =
          EncryptBlock(key_e_, label, block_id, payload, config_.block_size);
      blk_proto->set_label(std::move(label));
      blk_proto->set_ciphertext(blk.data(), blk.size());
    } else {
      uint64_t dummy_idx =
          static_cast<uint64_t>(item_idx - config_.num_blocks + 1);
      std::string label = GenerateDummyLabel(key_l_0, dummy_idx);
      std::vector<uint8_t> blk =
          GenerateDummyBlock(key_d_, label, ciphertext_block_size_);
      blk_proto->set_label(std::move(label));
      blk_proto->set_ciphertext(blk.data(), blk.size());
    }
  }

  auto response_or = server_->HandleInitMain(request);
  if (!response_or.ok()) return response_or.status();
  if (!response_or->success()) {
    return absl::InternalError(absl::StrFormat(
        "Server failed to initialize Main: %s", response_or->error_message()));
  }

  initialized_ = true;
  stats_.Reset();
  return absl::OkStatus();
}

absl::StatusOr<std::vector<uint8_t>> MapLatOramClient::QueryTargetBlock(
    BlockId block_id, size_t target_container) {
  LatOramQueryRequest request;
  std::vector<std::string> query_labels(cp_ + 1);
  std::vector<std::vector<uint8_t>> dummy_pads(cp_ + 1);

  for (size_t j = 0; j <= cp_; ++j) {
    std::string key_l_j = DeriveSubKey(key_l_, static_cast<uint64_t>(j));
    if (j == target_container) {
      query_labels[j] = GenerateBlockLabel(key_l_j, block_id);
    } else {
      query_labels[j] = GenerateDummyLabel(key_l_j, cd_[j]++);
      dummy_pads[j] =
          GenerateDummyBlock(key_d_, query_labels[j], ciphertext_block_size_);
    }
    request.add_labels(query_labels[j]);
  }

  if (retrieval_client_) {
    size_t target_idx =
        (target_container != static_cast<size_t>(-1)) ? target_container : 0;
    auto query_or = retrieval_client_->CreateQuery(target_idx, cp_ + 1,
                                                   ciphertext_block_size_);
    if (!query_or.ok()) return query_or.status();
    *request.mutable_retrieval_query() = std::move(*query_or);
  }

  size_t req_sz = request.ByteSizeLong();
  stats_.bytes_written += req_sz;
  stats_.online_bytes_uploaded += req_sz;

  auto response_or = server_->HandleQuery(request);
  if (!response_or.ok()) return response_or.status();
  if (!response_or->success()) {
    return absl::InternalError(absl::StrFormat(
        "Server failed to execute query: %s", response_or->error_message()));
  }

  size_t resp_sz = response_or->ByteSizeLong();
  stats_.bytes_read += resp_sz;
  stats_.online_bytes_downloaded += resp_sz;
  stats_.online_bytes_transferred += req_sz + resp_sz;
  stats_.server_xor_operations += cp_;

  std::vector<uint8_t> z;
  if (response_or->has_retrieval_response() && retrieval_client_) {
    auto dec_or = retrieval_client_->DecodeResponse(
        response_or->retrieval_response(), target_container,
        ciphertext_block_size_, dummy_pads);
    if (!dec_or.ok()) return dec_or.status();
    z = std::move(*dec_or);
  }

  std::string dec_label = query_labels[target_container];
  auto decrypted_or = DecryptBlock(key_e_, dec_label, z, config_.block_size);
  if (!decrypted_or.ok()) {
    return decrypted_or.status();
  }

  if (decrypted_or->block_id != block_id) {
    return absl::InternalError(absl::StrFormat(
        "Decrypted block ID %lu does not match requested block_id %lu",
        decrypted_or->block_id, block_id));
  }

  return decrypted_or->payload;
}

absl::Status MapLatOramClient::AddLateralStash(
    size_t cp,
    const std::vector<std::pair<BlockId, std::vector<uint8_t>>>& stash_items) {
  size_t num_dummies = (effective_p_ - cp) * effective_m_;
  std::string key_l_cp = DeriveSubKey(key_l_, static_cast<uint64_t>(cp));

  LatOramAddLateralStashRequest request;
  request.set_phase(static_cast<uint32_t>(cp));

  uint64_t shuffle_seed = config_.seed + cp + 1;
  std::vector<LabeledBlock> lat_blocks =
      SmallShuffle(stash_items, num_dummies, key_e_, key_l_cp, key_d_,
                   config_.block_size, shuffle_seed, block_versions_);
  for (auto& blk : lat_blocks) {
    auto* p = request.add_blocks();
    p->set_label(std::move(blk.label));
    p->set_ciphertext(blk.ciphertext.data(), blk.ciphertext.size());
  }

  size_t add_req_sz = request.ByteSizeLong();
  stats_.bytes_written += add_req_sz;
  stats_.amortized_bytes_uploaded += add_req_sz;
  stats_.amortized_bytes_transferred += add_req_sz;

  auto response_or = server_->HandleAddLateralStash(request);
  if (!response_or.ok()) return response_or.status();
  if (!response_or->success()) {
    return absl::InternalError(
        absl::StrFormat("Server failed to add lateral stash: %s",
                        response_or->error_message()));
  }
  size_t add_resp_sz = response_or->ByteSizeLong();
  stats_.bytes_read += add_resp_sz;
  stats_.amortized_bytes_downloaded += add_resp_sz;
  stats_.amortized_bytes_transferred += add_resp_sz;

  return absl::OkStatus();
}

}  // namespace oram_okvs
