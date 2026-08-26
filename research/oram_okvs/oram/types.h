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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_TYPES_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_TYPES_H_

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "research/oram_okvs/pir/onion_pir_params.h"

namespace oram_okvs {

// Block identifier type.
using BlockId = uint64_t;

// Block content type.
using Block = std::vector<uint8_t>;

// Supported ORAM Access operations.
enum class AccessOp {
  kRead,          // Read existing block content.
  kWrite,         // Write new block content.
  kReadAndWrite,  // Write new content and return previous content.
};

// Returns a human-readable name for an AccessOp.
inline absl::string_view AccessOpToString(AccessOp op) {
  switch (op) {
    case AccessOp::kRead:
      return "Read";
    case AccessOp::kWrite:
      return "Write";
    case AccessOp::kReadAndWrite:
      return "ReadAndWrite";
  }
  return "Unknown";
}

// Performance and bandwidth statistics.
struct OramStats {
  size_t total_accesses = 0;
  size_t reads_count = 0;
  size_t writes_count = 0;
  size_t bytes_read = 0;
  size_t bytes_written = 0;
  size_t online_bytes_uploaded = 0;
  size_t online_bytes_downloaded = 0;
  size_t online_bytes_transferred = 0;
  size_t amortized_bytes_uploaded = 0;
  size_t amortized_bytes_downloaded = 0;
  size_t amortized_bytes_transferred = 0;
  size_t shuffles_count = 0;
  size_t main_reinits_count = 0;
  size_t server_xor_operations = 0;
  double client_time_ms = 0.0;
  double server_time_ms = 0.0;

  void Reset() { *this = OramStats{}; }

  double OnlineBytesPerAccess() const {
    return total_accesses == 0
               ? 0.0
               : static_cast<double>(online_bytes_transferred) / total_accesses;
  }

  double AmortizedBytesPerAccess() const {
    return total_accesses == 0
               ? 0.0
               : static_cast<double>(amortized_bytes_transferred) /
                     total_accesses;
  }

  double TotalBytesPerAccess() const {
    return total_accesses == 0
               ? 0.0
               : static_cast<double>(online_bytes_transferred +
                                     amortized_bytes_transferred) /
                     total_accesses;
  }

  double BandwidthInBlocks(size_t block_size) const {
    return (block_size == 0 || total_accesses == 0)
               ? 0.0
               : TotalBytesPerAccess() / static_cast<double>(block_size);
  }

  double ClientTimeMsPerAccess() const {
    return total_accesses == 0 ? 0.0 : client_time_ms / total_accesses;
  }

  double ServerTimeMsPerAccess() const {
    return total_accesses == 0 ? 0.0 : server_time_ms / total_accesses;
  }

  double CostPer1MAccesses(double network_egress_per_gb = 0.09,
                           double server_cpu_per_core_sec = 1.5e-5) const {
    double network_gb = (TotalBytesPerAccess() * 1e6) / 1e9;
    double network_cost = network_gb * network_egress_per_gb;
    double server_sec = (ServerTimeMsPerAccess() / 1000.0) * 1e6;
    double compute_cost = server_sec * server_cpu_per_core_sec;
    return network_cost + compute_cost;
  }
};

// OKVS Backend Type for oblivious containers.
enum class OkvsBackendType {
  kPolynomial =
      0,  // Interpolation-OC (exact size Z, 0% overhead, default for small Z).
  kRandomBand = 1,  // Random Band OC (linear encoding time).
};

// Storage backend type for LatORAM.
enum class LatOramStorageType {
  kMap = 0,    // Hash map with XOR-trick dummy cancellation (baseline).
  kOkLat = 1,  // OKLat (Algorithm 6 without PIR): KeyedOC Main + OkvsOC lateral
               // stashes.
  kOkPLat = 2,  // OKPLat (Algorithm 6 with Single-Query OnionPIR): KeyedOC Main
                // + OkvsOC lateral stashes.
};

// Configuration parameters for LatORAM.
struct LatOramConfig {
  size_t num_blocks = 1024;  // N: Total number of logical blocks.
  size_t block_size = 64;    // b: Block size in bytes.
  size_t stash_size = 0;  // M: Client stash capacity (0 = auto: ceil(sqrt(N))).
  size_t num_phases =
      0;  // P: Lateral phases before reinit (0 = auto: 2 * ceil(sqrt(N/M))).
  uint64_t seed = 0;  // Seed for deterministic PRNG (0 = random).
  LatOramStorageType storage_type = LatOramStorageType::kMap;
  OkvsBackendType okvs_backend = OkvsBackendType::kRandomBand;
  std::optional<pir::OnionPirParams> custom_pir_params;

  // Computes effective M and P if set to 0 (auto).
  size_t EffectiveStashSize() const {
    if (stash_size > 0) return stash_size;
    size_t m = static_cast<size_t>(
        std::ceil(std::sqrt(static_cast<double>(num_blocks))));
    return std::max<size_t>(m, 4);
  }

  size_t EffectiveNumPhases() const {
    if (num_phases > 0) return num_phases;
    size_t m = EffectiveStashSize();
    if (storage_type == LatOramStorageType::kOkLat ||
        storage_type == LatOramStorageType::kOkPLat) {
      // Dummy-free OKVS lateral stashes eliminate dummy explosion, unlocking
      // P = O(N/M) = O(sqrt(N)) and O(1) amortized background communication.
      size_t p = static_cast<size_t>(
          std::ceil(static_cast<double>(num_blocks) / static_cast<double>(m)));
      return std::max<size_t>(p, 2);
    }
    // Baseline LatORAM with dummy penalty is bounded by P = O(N^{1/4})
    double ratio = static_cast<double>(num_blocks) / static_cast<double>(m);
    size_t p = static_cast<size_t>(std::ceil(2.0 * std::sqrt(ratio)));
    return std::max<size_t>(p, 2);
  }
};

// Storage backend type for Ring ORAM.

enum class RingOramStorageType {
  kTree = 0,    // Standard Tree-based Ring ORAM with metadata & slot retrieval.
  kOkRing = 1,  // OKRing (Algorithm 3): (D+1) labels and candidate block
                // downloads without PIR.
  kOkPRing =
      2,  // OKPRing (Algorithm 4): Single-query OnionPIR with (pos, dist) map.
};

// Configuration parameters for Ring ORAM.
struct RingOramConfig {
  size_t num_blocks = 1024;  // N: Total number of logical blocks.
  size_t block_size = 64;    // Block size in bytes.
  size_t z_real = 2;         // Z: Real blocks per bucket.
  size_t s_dummy = 2;        // S: Dummy blocks per bucket.
  size_t eviction_rate = 1;  // A: Evictions per access operation.
  uint64_t seed = 0;         // PRNG seed (0 = random).
  RingOramStorageType storage_type = RingOramStorageType::kTree;
  OkvsBackendType okvs_backend = OkvsBackendType::kPolynomial;
  bool use_xor_trick = false;
  std::optional<pir::OnionPirParams> custom_pir_params;

  size_t BucketSize() const { return z_real + s_dummy; }

  size_t TreeDepth() const {
    if (num_blocks <= 1) return 1;
    if (eviction_rate <= 1) {
      return static_cast<size_t>(
          std::ceil(std::log2(static_cast<double>(num_blocks))));
    }
    double ratio = (2.0 * static_cast<double>(num_blocks)) /
                   static_cast<double>(eviction_rate);
    if (ratio <= 1.0) return 1;
    return static_cast<size_t>(std::ceil(std::log2(ratio)));
  }

  size_t NumLeaves() const { return static_cast<size_t>(1) << TreeDepth(); }

  size_t TotalBuckets() const {
    return (static_cast<size_t>(1) << (TreeDepth() + 1)) - 1;
  }
};

// Metadata for a single slot in a Ring ORAM bucket.
struct BucketSlotMetadata {
  BucketSlotMetadata() = default;
  BucketSlotMetadata(BlockId id, bool dummy, uint32_t dummy_idx)
      : block_id(id), is_dummy(dummy), dummy_index(dummy_idx) {}

  BlockId block_id = 0;
  bool is_dummy = true;
  uint32_t dummy_index = 0;
};

// Labeled block consisting of a lookup label and encrypted block ciphertext.
struct LabeledBlock {
  LabeledBlock() = default;
  LabeledBlock(std::string l, std::vector<uint8_t> c)
      : label(std::move(l)), ciphertext(std::move(c)) {}

  std::string label;
  std::vector<uint8_t> ciphertext;
};

// Block entry stored in the Ring ORAM client stash.
struct RingStashEntry {
  RingStashEntry() = default;
  RingStashEntry(BlockId id, uint32_t leaf, std::vector<uint8_t> p)
      : block_id(id), leaf_id(leaf), payload(std::move(p)) {}

  BlockId block_id = 0;
  uint32_t leaf_id = 0;
  std::vector<uint8_t> payload;
};

// Computes the sequence of bucket IDs on the path from root (level 0) to
// leaf_id (level depth).
inline std::vector<uint32_t> GetPathBucketIds(uint32_t leaf_id, size_t depth) {
  std::vector<uint32_t> path(depth + 1);
  for (size_t d = 0; d <= depth; ++d) {
    path[d] = (static_cast<uint32_t>(1) << d) - 1 + (leaf_id >> (depth - d));
  }
  return path;
}

// Reverses the lowest `num_bits` bits of `value` (for reverse-lexicographical
// eviction).
inline uint32_t ReverseBits(uint32_t value, size_t num_bits) {
  uint32_t res = 0;
  for (size_t i = 0; i < num_bits; ++i) {
    if ((value >> i) & 1) {
      res |= (static_cast<uint32_t>(1) << (num_bits - 1 - i));
    }
  }
  return res;
}

// Bitwise XOR of two byte spans.
inline std::vector<uint8_t> XorBlocks(absl::Span<const uint8_t> a,
                                      absl::Span<const uint8_t> b) {
  size_t n = std::min(a.size(), b.size());
  std::vector<uint8_t> res(n);
  for (size_t i = 0; i < n; ++i) {
    res[i] = a[i] ^ b[i];
  }
  return res;
}

// In-place bitwise XOR: dst[i] ^= src[i].
inline void XorInPlace(absl::Span<uint8_t> dst, absl::Span<const uint8_t> src) {
  size_t n = std::min(dst.size(), src.size());
  for (size_t i = 0; i < n; ++i) {
    dst[i] ^= src[i];
  }
}

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_TYPES_H_
