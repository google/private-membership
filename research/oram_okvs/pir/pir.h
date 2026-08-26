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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_PIR_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_PIR_H_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "research/oram_okvs/pir/pir.pb.h"

namespace oram_okvs::pir {

// Common block / entry representation.
using Entry = std::vector<uint8_t>;

// Generic PIR client interface.
//
// A PIR client can generate evaluation keys (sent to the server once during
// initialization), construct encrypted queries for target database indices,
// and decrypt the responses returned by the server.
class PirClient {
 public:
  virtual ~PirClient() = default;

  // Generates evaluation keys (Galois keys, RGSW secret keys, etc.) to be
  // stored or registered with the PIR server.
  virtual absl::StatusOr<SerializedEvaluationKeys> GetSerializedEvaluationKeys()
      const = 0;

  // Constructs an encrypted query targeting the specified database `index`.
  // `index` must be in [0, num_entries() - 1].
  virtual absl::StatusOr<SerializedPirQuery> CreateQuery(
      size_t index) const = 0;

  // Decrypts and decodes the server's PIR response to retrieve the original
  // entry payload.
  virtual absl::StatusOr<std::vector<uint8_t>> DecodeResponse(
      const SerializedPirResponse& response, size_t entry_size_bytes) const = 0;

  // Total number of database entries supported by this client configuration.
  virtual size_t num_entries() const = 0;

  // Expected database entry size in bytes.
  virtual size_t entry_size_bytes() const = 0;
};

// Generic PIR server interface.
//
// A PIR server stores or processes a database of entries, holds client
// evaluation keys, and evaluates encrypted queries over the database without
// learning which index is being retrieved.
class PirServer {
 public:
  virtual ~PirServer() = default;

  // Sets or updates the client evaluation keys needed to evaluate queries.
  virtual absl::Status SetEvaluationKeys(
      const SerializedEvaluationKeys& keys) = 0;

  // Loads or updates the server database.
  // Entries are represented as byte arrays.
  virtual absl::Status SetDatabase(
      absl::Span<const std::vector<uint8_t>> db) = 0;

  // Evaluates an encrypted PIR query over the currently loaded database.
  virtual absl::StatusOr<SerializedPirResponse> ProcessQuery(
      const SerializedPirQuery& query) const = 0;

  // Evaluates an encrypted PIR query over a provided database on the fly.
  // This allows dynamic databases (such as ORAM buckets or OKVS tables) to be
  // queried directly without an explicit SetDatabase() call.
  virtual absl::StatusOr<SerializedPirResponse> ProcessQuery(
      const SerializedPirQuery& query,
      absl::Span<const std::vector<uint8_t>> db) const = 0;

  // Total number of database entries configured.
  virtual size_t num_entries() const = 0;

  // Entry size in bytes.
  virtual size_t entry_size_bytes() const = 0;
};

}  // namespace oram_okvs::pir

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_PIR_H_
