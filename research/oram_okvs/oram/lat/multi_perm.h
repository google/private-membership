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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LAT_MULTI_PERM_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LAT_MULTI_PERM_H_

#include <cstddef>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include "absl/hash/hash.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"

namespace oram_okvs {

// MultiPermDS: Client-side prefix-sum data structure used to generate a
// uniform random multi-permutation on the fly without storing an entire
// O(N) array in memory.
//
// Uses a Fenwick Tree (Binary Indexed Tree) for O(log c) updates and queries.
class MultiPermDS {
 public:
  MultiPermDS() : total_remaining_(0), seed_(0), counter_(0) {}

  // Initializes with frequency array L = (L_0, L_1, ..., L_{c-1})
  // and an optional PRNG seed.
  explicit MultiPermDS(absl::Span<const size_t> frequencies,
                       uint64_t seed = 0) {
    Init(frequencies, seed);
  }

  // Resets and initializes the structure with new frequencies.
  void Init(absl::Span<const size_t> frequencies, uint64_t seed = 0) {
    size_t c = frequencies.size();
    frequencies_.assign(frequencies.begin(), frequencies.end());
    tree_.assign(c + 1, 0);
    total_remaining_ = 0;
    seed_ = seed;
    counter_ = 0;

    for (size_t i = 0; i < c; ++i) {
      if (frequencies[i] > 0) {
        Update(i + 1, static_cast<int64_t>(frequencies[i]));
        total_remaining_ += frequencies[i];
      }
    }
  }

  // Returns the total number of items remaining across all bins.
  size_t TotalRemaining() const { return total_remaining_; }

  // Returns the number of bins.
  size_t NumBins() const { return frequencies_.size(); }

  // Returns the remaining frequency for a specific bin.
  size_t Frequency(size_t index) const {
    if (index >= frequencies_.size()) return 0;
    return frequencies_[index];
  }

  // Samples the next bin index z in [0, NumBins() - 1], decrements its
  // frequency, and returns z.
  absl::StatusOr<size_t> Next() {
    if (total_remaining_ == 0) {
      return absl::FailedPreconditionError("MultiPermDS is empty.");
    }

    uint64_t r_raw = NextRandom();
    size_t r = static_cast<size_t>(r_raw % total_remaining_);

    // Binary lifting on Fenwick tree:
    // Finds the largest 1-based index `idx` such that prefix_sum(idx) <= r.
    size_t idx = 0;
    size_t bit_mask = 1;
    while ((bit_mask << 1) < tree_.size()) {
      bit_mask <<= 1;
    }

    int64_t current_sum = 0;
    for (; bit_mask > 0; bit_mask >>= 1) {
      size_t next_idx = idx + bit_mask;
      if (next_idx < tree_.size() &&
          current_sum + tree_[next_idx] <= static_cast<int64_t>(r)) {
        idx = next_idx;
        current_sum += tree_[next_idx];
      }
    }

    // The target bin is at 1-based index (idx + 1), which is 0-based `idx`.
    size_t z = idx;
    if (z >= frequencies_.size() || frequencies_[z] == 0) {
      // Fallback in case of edge cases
      for (size_t i = 0; i < frequencies_.size(); ++i) {
        if (frequencies_[i] > 0) {
          z = i;
          break;
        }
      }
    }

    // Decrement frequency for bin z
    frequencies_[z]--;
    total_remaining_--;
    Update(z + 1, -1);

    return z;
  }

 private:
  void Update(size_t idx, int64_t delta) {
    for (; idx < tree_.size(); idx += idx & -idx) {
      tree_[idx] += delta;
    }
  }

  int64_t QueryPrefix(size_t idx) const {
    int64_t sum = 0;
    for (; idx > 0; idx -= idx & -idx) {
      sum += tree_[idx];
    }
    return sum;
  }

  uint64_t NextRandom() {
    if (seed_ != 0) {
      return absl::HashOf(seed_, counter_++);
    }
    // Fallback pseudo-random using std::random_device / minstd
    static thread_local std::mt19937_64 default_rng(std::random_device{}());
    return default_rng();
  }

  std::vector<size_t> frequencies_;
  std::vector<int64_t> tree_;
  size_t total_remaining_;
  uint64_t seed_;
  uint64_t counter_;
};

}  // namespace oram_okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_ORAM_LAT_MULTI_PERM_H_
