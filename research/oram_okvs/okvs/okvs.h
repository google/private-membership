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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_OKVS_OKVS_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_OKVS_OKVS_H_

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"

namespace oram_okvs::okvs {

// SFINAE helper to detect container types with .data() and .size()
template <typename T, typename = void>
struct is_byte_container : std::false_type {};

template <typename T>
struct is_byte_container<T, std::void_t<decltype(std::declval<T>().data()),
                                        decltype(std::declval<T>().size())>>
    : std::true_type {};

template <typename T>
inline constexpr bool is_byte_container_v = is_byte_container<T>::value;

// Abstract base class template for Oblivious Key-Value Stores (OKVS).
//
// An OKVS encodes a mapping from a set of n key-value pairs {(k_1, v_1), ...,
// (k_n, v_n)} into a compact storage table S of size m = O(n).
//
// Core Properties:
//  - Correctness: Decode(k_i) = v_i for all encoded keys.
//  - Obliviousness: For any key k not in the encoded set, Decode(k) produces a
//    pseudo-random value that reveals no information about the encoded keys or
//    values.
//
// Template Parameters:
//  - Key: The type of keys to encode/decode (e.g. std::string, uint64_t,
//  Gf128).
//  - Value: The type of values stored in the OKVS (e.g. Gf128, absl::uint128,
//  uint64_t).
template <typename Key, typename Value>
class Okvs {
 public:
  virtual ~Okvs() = default;

  // Encodes a mapping from `keys` to `values` into the internal OKVS storage.
  // `keys` and `values` must have identical sizes.
  // Returns absl::OkStatus() on success, or an error status if encoding fails
  // (e.g. due to linear dependencies in the random system).
  virtual absl::Status Encode(absl::Span<const Key> keys,
                              absl::Span<const Value> values) = 0;

  // Convenience overload for vector/span of (Key, Value) pairs.
  virtual absl::Status Encode(
      absl::Span<const std::pair<Key, Value>> key_values) {
    std::vector<Key> keys;
    std::vector<Value> values;
    keys.reserve(key_values.size());
    values.reserve(key_values.size());
    for (const auto& [k, v] : key_values) {
      keys.push_back(k);
      values.push_back(v);
    }
    return Encode(keys, values);
  }

  // Decodes the value associated with `key` from the internal OKVS table.
  // Returns the value if key was encoded, or a pseudo-random value otherwise.
  virtual absl::StatusOr<Value> Decode(const Key& key) const = 0;

  // Returns the total number of variables (table size m) in the OKVS.
  virtual size_t TableSize() const = 0;

  // Returns the number of keys encoded in the most recent Encode call.
  virtual size_t NumKeys() const = 0;

  // Returns a const reference to the underlying encoded storage table.
  virtual const std::vector<Value>& Table() const = 0;

  // Returns a span over the underlying encoded storage table.
  virtual absl::Span<const Value> TableSpan() const = 0;

  // Sets/replaces the internal storage table (e.g. when received from a peer).
  virtual void SetTable(std::vector<Value> table) = 0;
};

}  // namespace oram_okvs::okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_OKVS_OKVS_H_
