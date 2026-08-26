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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_OKVS_RANDOM_BAND_OKVS_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_OKVS_RANDOM_BAND_OKVS_H_

#include <immintrin.h>
#include <openssl/sha.h>
#include <wmmintrin.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "absl/hash/hash.h"
#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "research/oram_okvs/okvs/gf128.h"
#include "research/oram_okvs/okvs/okvs.h"

namespace oram_okvs::okvs {

// ============================================================================
// Fixed-size Bit Set for Random Bands
// ============================================================================

template <size_t NumBits>
class BandBits {
 public:
  static constexpr size_t kNumBits = NumBits;
  static constexpr size_t kNumWords = (NumBits + 63) / 64;

  constexpr BandBits() : words_() {}

  void Set(size_t index, bool val = true) {
    if (index >= NumBits) return;
    if (val) {
      words_[index / 64] |= (1ULL << (index % 64));
    } else {
      words_[index / 64] &= ~(1ULL << (index % 64));
    }
  }

  bool Test(size_t index) const {
    if (index >= NumBits) return false;
    return (words_[index / 64] & (1ULL << (index % 64))) != 0;
  }

  bool IsZero() const {
    for (size_t i = 0; i < kNumWords; ++i) {
      if (words_[i] != 0) return false;
    }
    return true;
  }

  // Returns index of the lowest set bit, or NumBits if all bits are 0.
  int CountTrailingZeros() const {
    for (size_t i = 0; i < kNumWords; ++i) {
      if (words_[i] != 0) {
        return static_cast<int>(i * 64 + __builtin_ctzll(words_[i]));
      }
    }
    return static_cast<int>(NumBits);
  }

  // Shifts bits to the right by `shift` bit positions.
  BandBits& ShiftRight(size_t shift) {
    if (shift >= NumBits) {
      words_.fill(0);
      return *this;
    }
    if (shift == 0) return *this;

    size_t word_shift = shift / 64;
    size_t bit_shift = shift % 64;

    for (size_t i = 0; i < kNumWords; ++i) {
      size_t src = i + word_shift;
      if (src >= kNumWords) {
        words_[i] = 0;
      } else if (bit_shift == 0) {
        words_[i] = words_[src];
      } else {
        uint64_t low = words_[src] >> bit_shift;
        uint64_t high =
            (src + 1 < kNumWords) ? (words_[src + 1] << (64 - bit_shift)) : 0;
        words_[i] = low | high;
      }
    }
    MaskTopWord();
    return *this;
  }

  BandBits& operator^=(const BandBits& other) {
    for (size_t i = 0; i < kNumWords; ++i) {
      words_[i] ^= other.words_[i];
    }
    return *this;
  }

  BandBits operator^(const BandBits& other) const {
    BandBits res = *this;
    res ^= other;
    return res;
  }

  bool operator==(const BandBits& other) const {
    return words_ == other.words_;
  }

  bool operator!=(const BandBits& other) const {
    return words_ != other.words_;
  }

  uint64_t word(size_t i) const { return words_[i]; }
  uint64_t& word(size_t i) { return words_[i]; }

  uint8_t* data() { return reinterpret_cast<uint8_t*>(words_.data()); }
  const uint8_t* data() const {
    return reinterpret_cast<const uint8_t*>(words_.data());
  }

  void MaskTopWord() {
    constexpr size_t kExtraBits = NumBits % 64;
    if constexpr (kExtraBits != 0) {
      constexpr uint64_t kMask = (1ULL << kExtraBits) - 1;
      words_[kNumWords - 1] &= kMask;
    }
  }

 private:
  std::array<uint64_t, kNumWords> words_;
};

// ============================================================================
// Value Type Traits
// ============================================================================

template <typename Value>
struct ValueTraits {
  static Value Zero() { return Value{}; }

  static bool IsZero(const Value& v) {
    if constexpr (std::is_same_v<Value, Gf128>) {
      return v.IsZero();
    } else if constexpr (std::is_integral_v<Value> ||
                         std::is_same_v<Value, absl::uint128>) {
      return v == 0;
    } else {
      return v == Value{};
    }
  }

  static void Accumulate(Value& acc, const Value& b) {
    if constexpr (std::is_same_v<Value, Gf128>) {
      acc = acc + b;
    } else if constexpr (std::is_integral_v<Value> ||
                         std::is_same_v<Value, absl::uint128>) {
      acc ^= b;
    } else {
      acc = acc ^ b;
    }
  }

  static Value Add(const Value& a, const Value& b) {
    Value res = a;
    Accumulate(res, b);
    return res;
  }
};

// Specialization for std::array<uint8_t, N>
template <size_t N>
struct ValueTraits<std::array<uint8_t, N>> {
  using Value = std::array<uint8_t, N>;
  static Value Zero() {
    Value arr;
    arr.fill(0);
    return arr;
  }

  static bool IsZero(const Value& v) {
    for (size_t i = 0; i < N; ++i) {
      if (v[i] != 0) return false;
    }
    return true;
  }

  static void Accumulate(Value& acc, const Value& b) {
    uint8_t* dst = acc.data();
    const uint8_t* src = b.data();
    size_t i = 0;

#if defined(__AVX2__)
    for (; i + 32 <= N; i += 32) {
      __m256i v_dst =
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(dst + i));
      __m256i v_src =
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src + i));
      _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst + i),
                          _mm256_xor_si256(v_dst, v_src));
    }
#elif defined(__SSE2__)
    for (; i + 16 <= N; i += 16) {
      __m128i v_dst =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(dst + i));
      __m128i v_src =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(src + i));
      _mm_storeu_si128(reinterpret_cast<__m128i*>(dst + i),
                       _mm_xor_si128(v_dst, v_src));
    }
#endif
    for (; i + 8 <= N; i += 8) {
      *reinterpret_cast<uint64_t*>(dst + i) ^=
          *reinterpret_cast<const uint64_t*>(src + i);
    }
    for (; i < N; ++i) {
      dst[i] ^= src[i];
    }
  }

  static Value Add(const Value& a, const Value& b) {
    Value res = a;
    Accumulate(res, b);
    return res;
  }
};

// Specialization for std::vector<uint8_t>
template <>
struct ValueTraits<std::vector<uint8_t>> {
  using Value = std::vector<uint8_t>;
  static Value Zero() { return {}; }

  static bool IsZero(const Value& v) {
    for (uint8_t b : v) {
      if (b != 0) return false;
    }
    return true;
  }

  static void Accumulate(Value& acc, const Value& b) {
    if (b.empty()) return;
    if (acc.empty()) {
      acc = b;
      return;
    }
    if (acc.size() < b.size()) {
      acc.resize(b.size(), 0);
    }
    uint8_t* dst = acc.data();
    const uint8_t* src = b.data();
    size_t sz = b.size();
    size_t i = 0;

#if defined(__AVX2__)
    for (; i + 32 <= sz; i += 32) {
      __m256i v_dst =
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(dst + i));
      __m256i v_src =
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src + i));
      _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst + i),
                          _mm256_xor_si256(v_dst, v_src));
    }
#elif defined(__SSE2__)
    for (; i + 16 <= sz; i += 16) {
      __m128i v_dst =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(dst + i));
      __m128i v_src =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(src + i));
      _mm_storeu_si128(reinterpret_cast<__m128i*>(dst + i),
                       _mm_xor_si128(v_dst, v_src));
    }
#endif
    for (; i + 8 <= sz; i += 8) {
      *reinterpret_cast<uint64_t*>(dst + i) ^=
          *reinterpret_cast<const uint64_t*>(src + i);
    }
    for (; i < sz; ++i) {
      dst[i] ^= src[i];
    }
  }

  static Value Add(const Value& a, const Value& b) {
    Value res = a;
    Accumulate(res, b);
    return res;
  }
};

template <typename Value>
Value GenerateRandomValue(uint64_t seed, size_t index) {
  uint8_t digest[SHA256_DIGEST_LENGTH];
  SHA256_CTX ctx;
  SHA256_Init(&ctx);
  SHA256_Update(&ctx, &seed, sizeof(seed));
  SHA256_Update(&ctx, &index, sizeof(index));
  SHA256_Final(digest, &ctx);

  if constexpr (std::is_same_v<Value, Gf128>) {
    uint64_t h1, h2;
    std::memcpy(&h1, digest, 8);
    std::memcpy(&h2, digest + 8, 8);
    return Gf128(h1, h2);
  } else if constexpr (std::is_same_v<Value, absl::uint128>) {
    uint64_t h1, h2;
    std::memcpy(&h1, digest, 8);
    std::memcpy(&h2, digest + 8, 8);
    return absl::MakeUint128(h1, h2);
  } else if constexpr (std::is_same_v<Value, std::vector<uint8_t>>) {
    return std::vector<uint8_t>(digest, digest + 16);
  } else if constexpr (std::is_integral_v<Value>) {
    Value res{};
    std::memcpy(&res, digest, std::min(sizeof(Value), sizeof(digest)));
    return res;
  } else {
    return ValueTraits<Value>::Zero();
  }
}

// ============================================================================
// Key Serializer
// ============================================================================

template <typename Key>
struct KeySerializer {
  static absl::string_view AsBytes(const Key& key, std::string* scratch) {
    if constexpr (std::is_same_v<Key, std::string> ||
                  std::is_same_v<Key, absl::string_view>) {
      return key;
    } else if constexpr (std::is_same_v<Key, const char*>) {
      return absl::string_view(key);
    } else if constexpr (std::is_trivially_copyable_v<Key>) {
      return absl::string_view(reinterpret_cast<const char*>(&key),
                               sizeof(Key));
    } else {
      size_t h = absl::Hash<Key>{}(key);
      scratch->assign(reinterpret_cast<const char*>(&h), sizeof(h));
      return *scratch;
    }
  }
};

// ============================================================================
// Key-to-Band Mapping (Cryptographically Secure SHA-512 PRF)
// ============================================================================

template <size_t BandBitWidth, typename Key>
inline std::pair<size_t, BandBits<BandBitWidth>> HashKeyToBand(const Key& key,
                                                               size_t num_vars,
                                                               uint64_t seed) {
  std::string scratch;
  absl::string_view key_bytes = KeySerializer<Key>::AsBytes(key, &scratch);

  constexpr size_t kBandBytes = (BandBitWidth + 7) / 8;
  constexpr size_t kTotalBytesNeeded = 8 + kBandBytes;

  uint8_t rand_bytes[kTotalBytesNeeded];

  if constexpr (kTotalBytesNeeded <= SHA512_DIGEST_LENGTH) {
    uint8_t digest[SHA512_DIGEST_LENGTH];
    SHA512_CTX ctx;
    SHA512_Init(&ctx);
    SHA512_Update(&ctx, &seed, sizeof(seed));
    uint32_t counter = 0;
    SHA512_Update(&ctx, &counter, sizeof(counter));
    SHA512_Update(&ctx, key_bytes.data(), key_bytes.size());
    SHA512_Final(digest, &ctx);
    std::memcpy(rand_bytes, digest, kTotalBytesNeeded);
  } else {
    size_t bytes_generated = 0;
    uint32_t counter = 0;
    while (bytes_generated < kTotalBytesNeeded) {
      SHA512_CTX ctx;
      SHA512_Init(&ctx);
      SHA512_Update(&ctx, &seed, sizeof(seed));
      SHA512_Update(&ctx, &counter, sizeof(counter));
      SHA512_Update(&ctx, key_bytes.data(), key_bytes.size());
      uint8_t digest[SHA512_DIGEST_LENGTH];
      SHA512_Final(digest, &ctx);

      size_t to_copy =
          std::min<size_t>(sizeof(digest), kTotalBytesNeeded - bytes_generated);
      std::memcpy(rand_bytes + bytes_generated, digest, to_copy);
      bytes_generated += to_copy;
      ++counter;
    }
  }

  uint64_t raw_pos = 0;
  std::memcpy(&raw_pos, rand_bytes, sizeof(uint64_t));

  size_t max_start =
      (num_vars >= BandBitWidth) ? (num_vars - BandBitWidth + 1) : 1;
  size_t start_pos = raw_pos % max_start;

  BandBits<BandBitWidth> band;
  std::memcpy(band.data(), rand_bytes + 8, kBandBytes);
  band.Set(0, true);
  band.MaskTopWord();

  return {start_pos, band};
}

// ============================================================================
// Random Band OKVS Class
// ============================================================================

// RandomBandOkvs implements the Random Band Oblivious Key-Value Store (RB-OKVS
// / PaXoS).
//
// Parameters:
//  - Key: The key type (e.g. std::string, uint64_t, Gf128).
//  - Value: The value type (e.g. Gf128, absl::uint128, uint64_t).
//  - BandBitWidth: Width w of the random band in bits (default: 350).
//    A band width of 350 achieves failure probability < 2^-40 for epsilon =
//    0.05.
template <typename Key, typename Value, size_t BandBitWidth = 350>
class RandomBandOkvs : public Okvs<Key, Value> {
 public:
  static constexpr size_t kBandBitWidth = BandBitWidth;

  // Constructs a RandomBandOkvs with a default expansion ratio epsilon
  // (e.g., 0.10 means m = (1 + 0.10) * n) and a PRF hash seed.
  explicit RandomBandOkvs(double epsilon = 0.10, uint64_t seed = 0)
      : epsilon_(epsilon), seed_(seed), num_keys_(0), num_vars_(0) {}

  // Constructs a RandomBandOkvs with a pre-configured table size.
  RandomBandOkvs(size_t num_variables, uint64_t seed = 0)
      : epsilon_(0.10),
        seed_(seed),
        num_keys_(0),
        num_vars_(num_variables),
        table_(num_variables, ValueTraits<Value>::Zero()) {}

  // Constructs a RandomBandOkvs from an existing encoded table.
  RandomBandOkvs(std::vector<Value> table, uint64_t seed = 0)
      : epsilon_(0.10),
        seed_(seed),
        num_keys_(0),
        num_vars_(table.size()),
        table_(std::move(table)) {}

  ~RandomBandOkvs() override = default;

  using Okvs<Key, Value>::Encode;

  // Encodes a mapping from `keys` to `values`. If `keys` is empty and
  // `default_value_len > 0`, initializes an unconstrained random table with
  // rows of size `default_value_len`.
  absl::Status Encode(absl::Span<const Key> keys,
                      absl::Span<const Value> values) override {
    return Encode(keys, values, 0);
  }

  absl::Status Encode(absl::Span<const Key> keys,
                      absl::Span<const Value> values,
                      size_t default_value_len) {
    if (keys.size() != values.size()) {
      return absl::InvalidArgumentError(
          "Keys and values must have the same length for OKVS encoding.");
    }

    size_t n = keys.size();
    num_keys_ = n;

    if (n == 0 && default_value_len == 0) {
      table_.clear();
      num_vars_ = 0;
      return absl::OkStatus();
    }

    // Determine storage size m = max(n + w, ceil(n * (1 + epsilon)))
    size_t required_vars = std::max(
        n + BandBitWidth, static_cast<size_t>(std::ceil(static_cast<double>(n) *
                                                        (1.0 + epsilon_))));
    if (num_vars_ < required_vars) {
      num_vars_ = required_vars;
    }

    return EncodeInternal(keys, values, num_vars_, seed_, default_value_len);
  }

  // Encodes keys and values, automatically retrying with different PRF seeds
  // up to `max_retries` if a rare linear dependency is encountered.
  absl::Status EncodeWithRetries(absl::Span<const Key> keys,
                                 absl::Span<const Value> values,
                                 int max_retries = 5,
                                 size_t default_value_len = 0) {
    uint64_t original_seed = seed_;
    for (int attempt = 0; attempt < max_retries; ++attempt) {
      seed_ = original_seed + attempt;
      absl::Status status = Encode(keys, values, default_value_len);
      if (status.ok()) {
        return absl::OkStatus();
      }
    }
    return absl::FailedPreconditionError(
        "OKVS encoding failed after maximum retry attempts.");
  }

  // Decodes the value associated with `key` from the internal table.
  absl::StatusOr<Value> Decode(const Key& key) const override {
    if (table_.empty()) {
      return absl::FailedPreconditionError(
          "OKVS table is empty or has not been encoded.");
    }
    return DecodeWithTable(key, table_, seed_);
  }

  static void AccumulateActiveBand(Value& dst, absl::Span<const Value> table,
                                   const uint32_t* active_indices,
                                   size_t num_active) {
    if constexpr (is_byte_container_v<Value>) {
      size_t sz = dst.size();
      uint8_t* dst_ptr = reinterpret_cast<uint8_t*>(dst.data());
      size_t offset = 0;

#if defined(__AVX2__)
      // 128-byte tiles (4x 256-bit YMM registers kept in registers across inner
      // loop)
      for (; offset + 128 <= sz; offset += 128) {
        __m256i acc0 = _mm256_setzero_si256();
        __m256i acc1 = _mm256_setzero_si256();
        __m256i acc2 = _mm256_setzero_si256();
        __m256i acc3 = _mm256_setzero_si256();

        for (size_t i = 0; i < num_active; ++i) {
          const uint8_t* src = reinterpret_cast<const uint8_t*>(
                                   table[active_indices[i]].data()) +
                               offset;
          acc0 = _mm256_xor_si256(
              acc0, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src)));
          acc1 = _mm256_xor_si256(
              acc1,
              _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src + 32)));
          acc2 = _mm256_xor_si256(
              acc2,
              _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src + 64)));
          acc3 = _mm256_xor_si256(
              acc3,
              _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src + 96)));
        }

        _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst_ptr + offset), acc0);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst_ptr + offset + 32),
                            acc1);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst_ptr + offset + 64),
                            acc2);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst_ptr + offset + 96),
                            acc3);
      }
      // 32-byte chunks
      for (; offset + 32 <= sz; offset += 32) {
        __m256i acc = _mm256_setzero_si256();
        for (size_t i = 0; i < num_active; ++i) {
          const uint8_t* src = reinterpret_cast<const uint8_t*>(
                                   table[active_indices[i]].data()) +
                               offset;
          acc = _mm256_xor_si256(
              acc, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src)));
        }
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst_ptr + offset), acc);
      }
#elif defined(__SSE2__)
      for (; offset + 16 <= sz; offset += 16) {
        __m128i acc = _mm_setzero_si128();
        for (size_t i = 0; i < num_active; ++i) {
          const uint8_t* src = reinterpret_cast<const uint8_t*>(
                                   table[active_indices[i]].data()) +
                               offset;
          acc = _mm_xor_si128(
              acc, _mm_loadu_si128(reinterpret_cast<const __m128i*>(src)));
        }
        _mm_storeu_si128(reinterpret_cast<__m128i*>(dst_ptr + offset), acc);
      }
#endif
      for (; offset + 8 <= sz; offset += 8) {
        uint64_t acc = 0;
        for (size_t i = 0; i < num_active; ++i) {
          acc ^= *reinterpret_cast<const uint64_t*>(
              reinterpret_cast<const uint8_t*>(
                  table[active_indices[i]].data()) +
              offset);
        }
        *reinterpret_cast<uint64_t*>(dst_ptr + offset) = acc;
      }
      for (; offset < sz; ++offset) {
        uint8_t acc = 0;
        for (size_t i = 0; i < num_active; ++i) {
          acc ^= reinterpret_cast<const uint8_t*>(
              table[active_indices[i]].data())[offset];
        }
        dst_ptr[offset] = acc;
      }
    } else {
      for (size_t i = 0; i < num_active; ++i) {
        ValueTraits<Value>::Accumulate(dst, table[active_indices[i]]);
      }
    }
  }

  // Static helper to decode a key given an external table span.
  static absl::StatusOr<Value> DecodeWithTable(const Key& key,
                                               absl::Span<const Value> table,
                                               uint64_t seed = 0) {
    if (table.empty()) {
      return absl::InvalidArgumentError("Table is empty.");
    }

    auto [start_pos, band] =
        HashKeyToBand<BandBitWidth>(key, table.size(), seed);

    Value res = ValueTraits<Value>::Zero();
    if constexpr (std::is_same_v<Value, std::vector<uint8_t>>) {
      if (!table.empty() && !table[0].empty()) {
        res.resize(table[0].size(), 0);
      }
    }

    uint32_t active_indices[BandBitWidth];
    size_t num_active = 0;
    constexpr size_t kNumWords = BandBits<BandBitWidth>::kNumWords;

    for (size_t w = 0; w < kNumWords; ++w) {
      uint64_t word = band.word(w);
      size_t base_idx = start_pos + w * 64;
      while (word != 0) {
        int bit = __builtin_ctzll(word);
        size_t var_idx = base_idx + bit;
        if (var_idx < table.size()) {
          active_indices[num_active++] = static_cast<uint32_t>(var_idx);
        }
        word &= word - 1;
      }
    }

    AccumulateActiveBand(res, table, active_indices, num_active);
    return res;
  }

  // Accessors.
  size_t TableSize() const override { return table_.size(); }
  size_t NumKeys() const override { return num_keys_; }
  const std::vector<Value>& Table() const override { return table_; }
  absl::Span<const Value> TableSpan() const override { return table_; }

  void SetTable(std::vector<Value> table) override {
    num_vars_ = table.size();
    table_ = std::move(table);
  }

  double epsilon() const { return epsilon_; }
  void set_epsilon(double epsilon) { epsilon_ = epsilon; }

  uint64_t seed() const { return seed_; }
  void set_seed(uint64_t seed) { seed_ = seed; }

 private:
  struct Equation {
    size_t start_pos;
    BandBits<BandBitWidth> band;
    Value value;
  };

  absl::Status EncodeInternal(absl::Span<const Key> keys,
                              absl::Span<const Value> values, size_t num_vars,
                              uint64_t seed, size_t default_value_len = 0) {
    size_t n = keys.size();
    std::vector<Equation> eqns(n);

    for (size_t i = 0; i < n; ++i) {
      auto [pos, band] = HashKeyToBand<BandBitWidth>(keys[i], num_vars, seed);
      eqns[i] = {pos, band, values[i]};
    }

    // Sort equations by starting position to achieve band-triangular form.
    std::sort(eqns.begin(), eqns.end(),
              [](const Equation& a, const Equation& b) {
                return a.start_pos < b.start_pos;
              });

    // 1. Forward Elimination (Row Echelon Reduction)
    std::vector<BandBits<BandBitWidth>> pivot_rows(num_vars);
    std::vector<Value> pivot_values(num_vars, ValueTraits<Value>::Zero());

    for (size_t i = 0; i < n; ++i) {
      size_t offset = eqns[i].start_pos;
      BandBits<BandBitWidth> band = eqns[i].band;
      Value val = eqns[i].value;

      while (!band.IsZero()) {
        int ctz = band.CountTrailingZeros();
        offset += ctz;
        band.ShiftRight(ctz);

        if (offset >= num_vars) {
          return absl::InternalError(
              "Row offset out of bounds during Gaussian elimination.");
        }

        if (pivot_rows[offset].IsZero()) {
          // Found pivot for column 'offset'.
          pivot_rows[offset] = band;
          pivot_values[offset] = val;
          break;
        }

        // Eliminate leading 1 with existing pivot row.
        band ^= pivot_rows[offset];
        ValueTraits<Value>::Accumulate(val, pivot_values[offset]);
      }

      if (band.IsZero()) {
        if (!ValueTraits<Value>::IsZero(val)) {
          return absl::FailedPreconditionError(
              "OKVS encoding failed: linear dependency with inconsistent "
              "values detected.");
        }
      }
    }

    // 2. Back Substitution (Solving for Table S)
    table_.assign(num_vars, ValueTraits<Value>::Zero());
    constexpr size_t kNumWords = BandBits<BandBitWidth>::kNumWords;

    size_t value_len = default_value_len;
    if constexpr (is_byte_container_v<Value>) {
      if (!values.empty()) {
        value_len = values[0].size();
      }
    }

    for (int64_t i = static_cast<int64_t>(num_vars) - 1; i >= 0; --i) {
      if (pivot_rows[i].IsZero()) {
        if constexpr (std::is_same_v<Value, std::vector<uint8_t>>) {
          std::vector<uint8_t> v(value_len);
          size_t bytes_done = 0;
          uint32_t ctr = 0;
          while (bytes_done < value_len) {
            uint8_t digest[SHA256_DIGEST_LENGTH];
            SHA256_CTX ctx;
            SHA256_Init(&ctx);
            SHA256_Update(&ctx, &seed, sizeof(seed));
            SHA256_Update(&ctx, &i, sizeof(i));
            SHA256_Update(&ctx, &ctr, sizeof(ctr));
            SHA256_Final(digest, &ctx);
            size_t to_copy =
                std::min<size_t>(sizeof(digest), value_len - bytes_done);
            std::memcpy(v.data() + bytes_done, digest, to_copy);
            bytes_done += to_copy;
            ++ctr;
          }
          table_[i] = std::move(v);
        } else {
          table_[i] = GenerateRandomValue<Value>(seed, static_cast<size_t>(i));
        }
        continue;
      }

      Value sum = pivot_values[i];
      const auto& band = pivot_rows[i];

      for (size_t w = 0; w < kNumWords; ++w) {
        uint64_t word = band.word(w);
        if (w == 0) {
          word &= ~1ULL;  // Exclude bit 0 (variable i itself)
        }
        size_t base_idx = static_cast<size_t>(i) + w * 64;
        while (word != 0) {
          int bit = __builtin_ctzll(word);
          size_t var_idx = base_idx + bit;
          if (var_idx < num_vars) {
            ValueTraits<Value>::Accumulate(sum, table_[var_idx]);
          }
          word &= word - 1;
        }
      }
      table_[i] = std::move(sum);
    }

    return absl::OkStatus();
  }

  double epsilon_;
  uint64_t seed_;
  size_t num_keys_;
  size_t num_vars_;
  std::vector<Value> table_;
};

}  // namespace oram_okvs::okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_OKVS_RANDOM_BAND_OKVS_H_
