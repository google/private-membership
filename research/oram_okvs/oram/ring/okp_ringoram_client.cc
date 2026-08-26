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

#include "research/oram_okvs/oram/ring/okp_ringoram_client.h"

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
#include "research/oram_okvs/okvs/polynomial_okvs.h"
#include "research/oram_okvs/okvs/random_band_okvs.h"
#include "research/oram_okvs/oram/block_retrieval.h"
#include "research/oram_okvs/oram/crypto_utils.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/ring/ringoram_client.h"
#include "research/oram_okvs/oram/ring/ringoram_server.h"
#include "research/oram_okvs/oram/types.h"
#include "research/oram_okvs/pir/onion_pir_params.h"

namespace oram_okvs {

OkPRingOramClient::OkPRingOramClient(
    const RingOramConfig& config,
    std::shared_ptr<RingOramServerInterface> server)
    : RingOramClient(config, std::move(server)) {
  distance_map_.assign(config_.num_blocks, -1);
  bucket_blocks_.resize(config_.TotalBuckets());
}

absl::Status OkPRingOramClient::Init() {
  stats_.Reset();
  stash_.clear();
  eviction_counter_ = 0;
  distance_map_.assign(config_.num_blocks, -1);
  block_versions_.assign(config_.num_blocks, 0);
  access_seq_ = 0;
  for (auto& blocks : bucket_blocks_) {
    blocks.clear();
  }

  for (BlockId i = 0; i < config_.num_blocks; ++i) {
    position_map_[i] = static_cast<uint32_t>(rng_() % num_leaves_);
  }

  size_t total_buckets = config_.TotalBuckets();
  RingOramInitTreeRequest init_req;
  init_req.set_tree_depth(static_cast<uint32_t>(tree_depth_));
  init_req.set_bucket_size(static_cast<uint32_t>(bucket_size_));

  for (uint32_t b_id = 0; b_id < total_buckets; ++b_id) {
    std::vector<std::pair<BlockId, std::vector<uint8_t>>> empty_real_blocks;
    *init_req.add_initial_buckets() = BuildBucket(b_id, empty_real_blocks);
  }

  // Initialize OnionPIR client
  size_t path_len = tree_depth_ + 1;
  size_t pir_entry_size = ((ciphertext_block_size_ + 15) / 16) * 16;
  if (config_.custom_pir_params.has_value()) {
    pir_params_ = *config_.custom_pir_params;
    pir_params_.num_entries = path_len;
    pir_params_.entry_size_bytes = pir_entry_size;
  } else {
    pir_params_ =
        pir::OnionPirParams::CreateFastForTesting(path_len, pir_entry_size);
    if (path_len > pir_params_.n0) {
      pir_params_.n0 = 8;
    }
  }
  auto ret_client_or = BlockRetrievalClient::CreateOnionPir(pir_params_);
  if (!ret_client_or.ok()) return ret_client_or.status();
  retrieval_client_ = std::move(*ret_client_or);

  auto eval_keys_or = retrieval_client_->GetSerializedEvaluationKeys();
  if (!eval_keys_or.ok()) return eval_keys_or.status();
  *init_req.mutable_pir_eval_keys() = std::move(*eval_keys_or);

  auto init_resp_or = server_->HandleInitTree(init_req);
  if (!init_resp_or.ok()) return init_resp_or.status();
  if (!init_resp_or->success()) {
    return absl::InternalError(
        absl::StrFormat("Server failed to initialize Ring ORAM tree: %s",
                        init_resp_or->error_message()));
  }

  stats_.Reset();
  return absl::OkStatus();
}

RingOramBucketProto OkPRingOramClient::BuildBucket(
    uint32_t bucket_id,
    const std::vector<std::pair<BlockId, std::vector<uint8_t>>>& real_blocks) {
  RingOramBucketProto bucket;
  bucket.set_bucket_id(bucket_id);
  bucket.set_access_count(0);
  bucket_blocks_[bucket_id].clear();

  std::vector<std::string> keys;
  std::vector<std::vector<uint8_t>> values;
  keys.reserve(real_blocks.size());
  values.reserve(real_blocks.size());

  for (const auto& [b_id, payload] : real_blocks) {
    std::string label = GenerateBlockLabel(key_l_, b_id);
    uint64_t ver = (b_id < block_versions_.size()) ? block_versions_[b_id] : 0;
    std::vector<uint8_t> ct =
        EncryptBlock(key_e_, label, b_id, payload, config_.block_size, ver);

    keys.push_back(std::move(label));
    values.push_back(std::move(ct));
    bucket_blocks_[bucket_id].push_back(b_id);
  }

  uint64_t okvs_seed = config_.seed + bucket_id + 1;
  if (config_.okvs_backend == OkvsBackendType::kPolynomial) {
    bucket.set_okvs_backend(0);
    okvs::PolynomialOkvs<std::string, std::vector<uint8_t>> okvs(
        /*fixed_size=*/config_.z_real, okvs_seed,
        /*default_val_len=*/ciphertext_block_size_);
    auto enc_status = okvs.Encode(keys, values);
    if (enc_status.ok()) {
      bucket.set_okvs_table(SerializeOkvsTable(okvs.Table()));
      bucket.set_okvs_seed(okvs.seed());
      bucket.set_num_real(static_cast<uint32_t>(real_blocks.size()));
    }
  } else {
    bucket.set_okvs_backend(1);
    okvs::RandomBandOkvs<std::string, std::vector<uint8_t>> okvs(0.15,
                                                                 okvs_seed);
    auto enc_status = okvs.EncodeWithRetries(keys, values);
    if (enc_status.ok()) {
      bucket.set_okvs_table(SerializeOkvsTable(okvs.Table()));
      bucket.set_okvs_seed(okvs.seed());
      bucket.set_num_real(static_cast<uint32_t>(real_blocks.size()));
    }
  }

  return bucket;
}

absl::StatusOr<std::vector<uint8_t>> OkPRingOramClient::ReadBlockFromPath(
    BlockId block_id, uint32_t old_leaf, const std::vector<uint32_t>& path) {
  int target_depth =
      (block_id < distance_map_.size()) ? distance_map_[block_id] : -1;
  std::string target_label;
  if (target_depth >= 0) {
    target_label = GenerateBlockLabel(key_l_, block_id);
  } else {
    target_label = GenerateRandomBytes(16);
  }

  RingOramReadPathBlocksRequest blk_req;
  blk_req.set_leaf_id(old_leaf);
  blk_req.add_okvs_queries(target_label);

  size_t path_len = tree_depth_ + 1;
  size_t pir_target = (target_depth >= 0) ? target_depth : 0;
  auto query_or = retrieval_client_->CreateQuery(pir_target, path_len,
                                                 ciphertext_block_size_);
  if (!query_or.ok()) return query_or.status();
  *blk_req.mutable_retrieval_query() = std::move(*query_or);

  size_t req_sz = blk_req.ByteSizeLong();
  stats_.bytes_written += req_sz;
  stats_.online_bytes_uploaded += req_sz;
  stats_.online_bytes_transferred += req_sz;

  auto blk_resp_or = server_->HandleReadPathBlocks(blk_req);
  if (!blk_resp_or.ok()) return blk_resp_or.status();
  size_t resp_sz = blk_resp_or->ByteSizeLong();
  stats_.bytes_read += resp_sz;
  stats_.online_bytes_downloaded += resp_sz;
  stats_.online_bytes_transferred += resp_sz;

  if (target_depth >= 0 && blk_resp_or->has_retrieval_response()) {
    auto dec_resp_or =
        retrieval_client_->DecodeResponse(blk_resp_or->retrieval_response(),
                                          target_depth, ciphertext_block_size_);
    if (dec_resp_or.ok()) {
      auto dec_or =
          DecryptBlock(key_e_, target_label, *dec_resp_or, config_.block_size);
      if (dec_or.ok() && dec_or->block_id == block_id) {
        return dec_or->payload;
      }
    }
  }

  return absl::NotFoundError("Block not in server path");
}

absl::Status OkPRingOramClient::ExtractEvictionBlocks(
    const std::vector<uint32_t>& path,
    const RingOramReadPathFullResponse& full_resp) {
  for (size_t d = 0; d < path.size(); ++d) {
    uint32_t b_id = path[d];
    if (static_cast<int>(d) >= full_resp.buckets_size()) continue;
    const auto& bucket = full_resp.buckets(d);

    if (!bucket.okvs_table().empty()) {
      auto table = DeserializeOkvsTable(bucket.okvs_table());
      if (bucket.okvs_backend() == 0) {
        okvs::PolynomialOkvs<std::string, std::vector<uint8_t>> okvs(
            /*fixed_size=*/config_.z_real, bucket.okvs_seed(),
            /*default_val_len=*/ciphertext_block_size_);
        okvs.SetTable(table);
        for (BlockId b_id_val : bucket_blocks_[b_id]) {
          std::string label = GenerateBlockLabel(key_l_, b_id_val);
          auto val_or = okvs.Decode(label);
          if (val_or.ok()) {
            auto dec_or =
                DecryptBlock(key_e_, label, *val_or, config_.block_size);
            if (dec_or.ok() && dec_or->block_id < config_.num_blocks) {
              uint64_t prev_ver = (b_id_val < block_versions_.size())
                                      ? block_versions_[b_id_val]
                                      : 0;
              if (dec_or->version >= prev_ver) {
                uint32_t cur_leaf = position_map_[b_id_val];
                stash_[b_id_val] = {b_id_val, cur_leaf,
                                    std::move(dec_or->payload)};
                distance_map_[b_id_val] = -1;
              }
            }
          }
        }
      } else {
        okvs::RandomBandOkvs<std::string, std::vector<uint8_t>> okvs(
            0.15, bucket.okvs_seed());
        okvs.SetTable(table);
        for (BlockId b_id_val : bucket_blocks_[b_id]) {
          std::string label = GenerateBlockLabel(key_l_, b_id_val);
          auto val_or = okvs.Decode(label);
          if (val_or.ok()) {
            auto dec_or =
                DecryptBlock(key_e_, label, *val_or, config_.block_size);
            if (dec_or.ok() && dec_or->block_id < config_.num_blocks) {
              uint64_t prev_ver = (b_id_val < block_versions_.size())
                                      ? block_versions_[b_id_val]
                                      : 0;
              if (dec_or->version >= prev_ver) {
                uint32_t cur_leaf = position_map_[b_id_val];
                stash_[b_id_val] = {b_id_val, cur_leaf,
                                    std::move(dec_or->payload)};
                distance_map_[b_id_val] = -1;
              }
            }
          }
        }
      }
    }
    bucket_blocks_[b_id].clear();
  }
  return absl::OkStatus();
}

void OkPRingOramClient::OnEviction(
    const std::vector<uint32_t>& path,
    const std::vector<std::vector<BlockId>>& placements) {
  for (size_t d = 0; d < placements.size(); ++d) {
    for (BlockId b_id : placements[d]) {
      if (b_id < distance_map_.size()) {
        distance_map_[b_id] = static_cast<int32_t>(d);
      }
    }
  }
}

void OkPRingOramClient::OnBlockRead(BlockId block_id) {
  if (block_id < distance_map_.size()) {
    distance_map_[block_id] = -1;
  }
}

}  // namespace oram_okvs
