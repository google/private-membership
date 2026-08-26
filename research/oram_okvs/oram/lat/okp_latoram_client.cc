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

#include "research/oram_okvs/oram/lat/okp_latoram_client.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/strings/str_format.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/block_retrieval.h"
#include "research/oram_okvs/oram/crypto_utils.h"
#include "research/oram_okvs/oram/lat/latoram_server.h"
#include "research/oram_okvs/oram/lat/ok_latoram_client.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/types.h"
#include "research/oram_okvs/pir/onion_pir_params.h"

namespace oram_okvs {

OkPLatOramClient::OkPLatOramClient(
    const LatOramConfig& config, std::shared_ptr<LatOramServerInterface> server)
    : OkLatOramClient(config, std::move(server)) {}

absl::Status OkPLatOramClient::Init(
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

  // 4. Initialize OnionPIR client
  size_t num_pir_entries = effective_p_ + 1;
  if (config_.custom_pir_params.has_value()) {
    pir_params_ = *config_.custom_pir_params;
    pir_params_.num_entries = num_pir_entries;
    pir_params_.entry_size_bytes = ciphertext_block_size_;
  } else {
    pir_params_ = pir::OnionPirParams::CreateFastForTesting(
        num_pir_entries, ciphertext_block_size_);
    if (num_pir_entries > pir_params_.n0) {
      pir_params_.n0 = 8;
    }
  }
  auto ret_client_or = BlockRetrievalClient::CreateOnionPir(pir_params_);
  if (!ret_client_or.ok()) return ret_client_or.status();
  retrieval_client_ = std::move(*ret_client_or);

  auto eval_keys_or = retrieval_client_->GetSerializedEvaluationKeys();
  if (!eval_keys_or.ok()) return eval_keys_or.status();
  *request.mutable_pir_eval_keys() = std::move(*eval_keys_or);

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

}  // namespace oram_okvs
