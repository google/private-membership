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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_ORAM_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_ORAM_H_

#include <cstddef>
#include <cstdint>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {

// Abstract interface for Oblivious RAM (ORAM).
//
// An ORAM provides an array abstraction of `num_blocks()` logical blocks, each
// of size `block_size()` bytes. Logical block addresses are in the range
// [0, num_blocks() - 1].
//
// All operations ensure access pattern obliviousness: an observer watching the
// storage accesses and network communication cannot distinguish which logical
// block is being accessed, or whether the operation is a Read or Write.
class Oram {
 public:
  virtual ~Oram() = default;

  // Reads the data stored at `block_id`.
  // Returns the block content, or an error status if `block_id` is invalid or
  // an internal protocol failure occurs.
  virtual absl::StatusOr<std::vector<uint8_t>> Read(BlockId block_id) = 0;

  // Writes `data` to the block at `block_id`.
  // If `data` is smaller than `block_size()`, it is zero-padded.
  virtual absl::Status Write(BlockId block_id,
                             absl::Span<const uint8_t> data) = 0;

  // Performs a general Access operation (Read, Write, or ReadAndWrite).
  //  - kRead: returns the current data stored at `block_id`.
  //  - kWrite: overwrites `block_id` with `data` and returns the new data.
  //  - kReadAndWrite: overwrites `block_id` with `data` and returns the
  //  previous data.
  virtual absl::StatusOr<std::vector<uint8_t>> Access(
      BlockId block_id, AccessOp op, absl::Span<const uint8_t> data = {}) = 0;

  // Returns the total number of logical blocks in this ORAM instance.
  virtual size_t num_blocks() const = 0;

  // Returns the fixed block payload size in bytes.
  virtual size_t block_size() const = 0;

  // Returns communication and access metrics.
  virtual OramStats GetStats() const = 0;

  // Resets communication and access metrics.
  virtual void ResetStats() = 0;
};

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_ORAM_H_
