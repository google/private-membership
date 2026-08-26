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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_BLOCK_RETRIEVAL_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_BLOCK_RETRIEVAL_H_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/pir/onion_pir_client.h"
#include "research/oram_okvs/pir/onion_pir_params.h"
#include "research/oram_okvs/pir/onion_pir_server.h"

namespace oram_okvs {

// Abstract base class for Client-side oblivious block retrieval.
class BlockRetrievalClient {
 public:
  virtual ~BlockRetrievalClient() = default;

  // Generates a retrieval query for the target candidate index.
  virtual absl::StatusOr<ObliviousRetrievalQueryProto> CreateQuery(
      size_t target_index, size_t num_candidates, size_t entry_size) = 0;

  // Reconstructs the target block ciphertext from the server's response.
  // `known_dummy_pads`: for XOR-trick dummy un-XORing (if applicable).
  virtual absl::StatusOr<std::vector<uint8_t>> DecodeResponse(
      const ObliviousRetrievalResponseProto& response, size_t target_index,
      size_t entry_size,
      absl::Span<const std::vector<uint8_t>> known_dummy_pads = {}) = 0;

  // Returns evaluation keys to initialize the server (for homomorphic PIR).
  virtual absl::StatusOr<pir::SerializedEvaluationKeys>
  GetSerializedEvaluationKeys() {
    return pir::SerializedEvaluationKeys{};
  }

  // Factory methods
  static std::unique_ptr<BlockRetrievalClient> CreateDirect();
  static std::unique_ptr<BlockRetrievalClient> CreateXor();
  static absl::StatusOr<std::unique_ptr<BlockRetrievalClient>> CreateOnionPir(
      const pir::OnionPirParams& params);
};

// Abstract base class for Server-side oblivious block retrieval.
class BlockRetrievalServer {
 public:
  virtual ~BlockRetrievalServer() = default;

  // Sets evaluation keys (for homomorphic PIR).
  virtual absl::Status SetEvaluationKeys(
      const pir::SerializedEvaluationKeys& keys) {
    return absl::OkStatus();
  }

  // Evaluates retrieval across the candidate blocks using the client's query.
  virtual absl::StatusOr<ObliviousRetrievalResponseProto> Evaluate(
      absl::Span<const std::vector<uint8_t>> candidate_blocks,
      const ObliviousRetrievalQueryProto& query) = 0;

  // Factory methods
  static std::unique_ptr<BlockRetrievalServer> Create();
  static std::unique_ptr<BlockRetrievalServer> CreateDirect();
  static std::unique_ptr<BlockRetrievalServer> CreateXor();
  static absl::StatusOr<std::unique_ptr<BlockRetrievalServer>> CreateOnionPir(
      const pir::OnionPirParams& params);
};

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_BLOCK_RETRIEVAL_H_
