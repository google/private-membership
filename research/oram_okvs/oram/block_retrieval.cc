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

#include "research/oram_okvs/oram/block_retrieval.h"

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
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/types.h"
#include "research/oram_okvs/pir/onion_pir_client.h"
#include "research/oram_okvs/pir/onion_pir_params.h"
#include "research/oram_okvs/pir/onion_pir_server.h"

namespace oram_okvs {

namespace {

// ============================================================================
// 1. Direct / Multi-Block Retrieval (Trivial PIR)
// ============================================================================

class DirectRetrievalClient : public BlockRetrievalClient {
 public:
  absl::StatusOr<ObliviousRetrievalQueryProto> CreateQuery(
      size_t target_index, size_t num_candidates, size_t entry_size) override {
    ObliviousRetrievalQueryProto query;
    query.set_mode(RetrievalModeProto::RETRIEVAL_MODE_DIRECT);
    return query;
  }

  absl::StatusOr<std::vector<uint8_t>> DecodeResponse(
      const ObliviousRetrievalResponseProto& response, size_t target_index,
      size_t entry_size,
      absl::Span<const std::vector<uint8_t>> known_dummy_pads) override {
    if (static_cast<int>(target_index) >= response.candidate_blocks_size()) {
      return absl::InvalidArgumentError(absl::StrFormat(
          "Target index %zu out of range in DirectRetrieval response (size %d)",
          target_index, response.candidate_blocks_size()));
    }
    const std::string& blk =
        response.candidate_blocks(static_cast<int>(target_index));
    std::vector<uint8_t> res(blk.begin(), blk.end());
    if (res.size() > entry_size) {
      res.resize(entry_size);
    }
    return res;
  }
};

class DirectRetrievalServer : public BlockRetrievalServer {
 public:
  absl::StatusOr<ObliviousRetrievalResponseProto> Evaluate(
      absl::Span<const std::vector<uint8_t>> candidate_blocks,
      const ObliviousRetrievalQueryProto& query) override {
    ObliviousRetrievalResponseProto response;
    response.set_mode(RetrievalModeProto::RETRIEVAL_MODE_DIRECT);
    for (const auto& blk : candidate_blocks) {
      response.add_candidate_blocks(blk.data(), blk.size());
    }
    return response;
  }
};

// ============================================================================
// 2. XOR-Trick Dummy Cancellation Retrieval
// ============================================================================

class XorRetrievalClient : public BlockRetrievalClient {
 public:
  absl::StatusOr<ObliviousRetrievalQueryProto> CreateQuery(
      size_t target_index, size_t num_candidates, size_t entry_size) override {
    ObliviousRetrievalQueryProto query;
    query.set_mode(RetrievalModeProto::RETRIEVAL_MODE_XOR);
    return query;
  }

  absl::StatusOr<std::vector<uint8_t>> DecodeResponse(
      const ObliviousRetrievalResponseProto& response, size_t target_index,
      size_t entry_size,
      absl::Span<const std::vector<uint8_t>> known_dummy_pads) override {
    const std::string& xor_bytes = response.xor_block();
    std::vector<uint8_t> z(xor_bytes.begin(), xor_bytes.end());
    if (z.size() < entry_size) {
      z.resize(entry_size, 0);
    }

    for (size_t j = 0; j < known_dummy_pads.size(); ++j) {
      if (j != target_index && !known_dummy_pads[j].empty()) {
        XorInPlace(absl::MakeSpan(z), known_dummy_pads[j]);
      }
    }
    if (z.size() > entry_size) {
      z.resize(entry_size);
    }
    return z;
  }
};

class XorRetrievalServer : public BlockRetrievalServer {
 public:
  absl::StatusOr<ObliviousRetrievalResponseProto> Evaluate(
      absl::Span<const std::vector<uint8_t>> candidate_blocks,
      const ObliviousRetrievalQueryProto& query) override {
    ObliviousRetrievalResponseProto response;
    response.set_mode(RetrievalModeProto::RETRIEVAL_MODE_XOR);
    if (candidate_blocks.empty()) {
      return response;
    }

    std::vector<uint8_t> xor_sum = candidate_blocks[0];
    for (size_t i = 1; i < candidate_blocks.size(); ++i) {
      XorInPlace(absl::MakeSpan(xor_sum), candidate_blocks[i]);
    }
    response.set_xor_block(xor_sum.data(), xor_sum.size());
    return response;
  }
};

// ============================================================================
// 3. Single-Query OnionPIR Retrieval
// ============================================================================

class OnionPirRetrievalClient : public BlockRetrievalClient {
 public:
  explicit OnionPirRetrievalClient(std::unique_ptr<pir::PirClient> pir_client,
                                   pir::OnionPirParams params)
      : pir_client_(std::move(pir_client)), params_(std::move(params)) {}

  absl::StatusOr<ObliviousRetrievalQueryProto> CreateQuery(
      size_t target_index, size_t num_candidates, size_t entry_size) override {
    ObliviousRetrievalQueryProto query;
    query.set_mode(RetrievalModeProto::RETRIEVAL_MODE_PIR);
    auto pir_q_or = pir_client_->CreateQuery(target_index);
    if (!pir_q_or.ok()) return pir_q_or.status();
    *query.mutable_pir_query() = std::move(*pir_q_or);
    return query;
  }

  absl::StatusOr<std::vector<uint8_t>> DecodeResponse(
      const ObliviousRetrievalResponseProto& response, size_t target_index,
      size_t entry_size,
      absl::Span<const std::vector<uint8_t>> known_dummy_pads) override {
    if (!response.has_pir_response()) {
      return absl::InvalidArgumentError(
          "Missing pir_response in PIR retrieval mode.");
    }
    auto decoded_or = pir_client_->DecodeResponse(response.pir_response(),
                                                  params_.entry_size_bytes);
    if (!decoded_or.ok()) return decoded_or.status();
    std::vector<uint8_t> res = std::move(*decoded_or);
    if (res.size() > entry_size) {
      res.resize(entry_size);
    }
    return res;
  }

  absl::StatusOr<pir::SerializedEvaluationKeys> GetSerializedEvaluationKeys()
      override {
    return pir_client_->GetSerializedEvaluationKeys();
  }

 private:
  std::unique_ptr<pir::PirClient> pir_client_;
  pir::OnionPirParams params_;
};

class UnifiedBlockRetrievalServer : public BlockRetrievalServer {
 public:
  UnifiedBlockRetrievalServer() = default;

  explicit UnifiedBlockRetrievalServer(
      std::unique_ptr<pir::PirServer> pir_server, pir::OnionPirParams params)
      : pir_server_(std::move(pir_server)), params_(std::move(params)) {}

  absl::Status SetEvaluationKeys(
      const pir::SerializedEvaluationKeys& keys) override {
    pir::OnionPirParams pir_params =
        pir::OnionPirParams::FromProto(keys.params());
    auto pir_server_or = pir::OnionPirServer::Create(pir_params);
    if (!pir_server_or.ok()) return pir_server_or.status();
    pir_server_ = std::move(*pir_server_or);
    params_ = pir_params;
    return pir_server_->SetEvaluationKeys(keys);
  }

  absl::StatusOr<ObliviousRetrievalResponseProto> Evaluate(
      absl::Span<const std::vector<uint8_t>> candidate_blocks,
      const ObliviousRetrievalQueryProto& query) override {
    ObliviousRetrievalResponseProto response;

    if (query.mode() == RetrievalModeProto::RETRIEVAL_MODE_PIR ||
        query.has_pir_query()) {
      response.set_mode(RetrievalModeProto::RETRIEVAL_MODE_PIR);
      if (!pir_server_) {
        return absl::FailedPreconditionError(
            "PIR evaluation keys not initialized on server.");
      }
      auto set_db_status = pir_server_->SetDatabase(candidate_blocks);
      if (!set_db_status.ok()) return set_db_status;

      auto pir_resp_or = pir_server_->ProcessQuery(query.pir_query());
      if (!pir_resp_or.ok()) return pir_resp_or.status();

      *response.mutable_pir_response() = std::move(*pir_resp_or);
      return response;
    }

    if (query.mode() == RetrievalModeProto::RETRIEVAL_MODE_XOR) {
      response.set_mode(RetrievalModeProto::RETRIEVAL_MODE_XOR);
      if (candidate_blocks.empty()) return response;
      std::vector<uint8_t> xor_sum = candidate_blocks[0];
      for (size_t i = 1; i < candidate_blocks.size(); ++i) {
        XorInPlace(absl::MakeSpan(xor_sum), candidate_blocks[i]);
      }
      response.set_xor_block(xor_sum.data(), xor_sum.size());
      return response;
    }

    // Default: RETRIEVAL_MODE_DIRECT
    response.set_mode(RetrievalModeProto::RETRIEVAL_MODE_DIRECT);
    for (const auto& blk : candidate_blocks) {
      response.add_candidate_blocks(blk.data(), blk.size());
    }
    return response;
  }

 private:
  std::unique_ptr<pir::PirServer> pir_server_;
  pir::OnionPirParams params_;
};

}  // namespace

std::unique_ptr<BlockRetrievalClient> BlockRetrievalClient::CreateDirect() {
  return std::make_unique<DirectRetrievalClient>();
}

std::unique_ptr<BlockRetrievalClient> BlockRetrievalClient::CreateXor() {
  return std::make_unique<XorRetrievalClient>();
}

absl::StatusOr<std::unique_ptr<BlockRetrievalClient>>
BlockRetrievalClient::CreateOnionPir(const pir::OnionPirParams& params) {
  auto pir_client_or = pir::OnionPirClient::Create(params);
  if (!pir_client_or.ok()) return pir_client_or.status();
  return std::make_unique<OnionPirRetrievalClient>(std::move(*pir_client_or),
                                                   params);
}

std::unique_ptr<BlockRetrievalServer> BlockRetrievalServer::Create() {
  return std::make_unique<UnifiedBlockRetrievalServer>();
}

std::unique_ptr<BlockRetrievalServer> BlockRetrievalServer::CreateDirect() {
  return std::make_unique<DirectRetrievalServer>();
}

std::unique_ptr<BlockRetrievalServer> BlockRetrievalServer::CreateXor() {
  return std::make_unique<XorRetrievalServer>();
}

absl::StatusOr<std::unique_ptr<BlockRetrievalServer>>
BlockRetrievalServer::CreateOnionPir(const pir::OnionPirParams& params) {
  auto pir_server_or = pir::OnionPirServer::Create(params);
  if (!pir_server_or.ok()) return pir_server_or.status();
  return std::make_unique<UnifiedBlockRetrievalServer>(
      std::move(*pir_server_or), params);
}

}  // namespace oram_okvs
