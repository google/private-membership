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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_OKVS_POLYNOMIAL_OKVS_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_OKVS_POLYNOMIAL_OKVS_H_

#include <immintrin.h>
#include <openssl/sha.h>
#include <wmmintrin.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
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
// Value Arithmetic Helpers for GF(2^128) Scalar Multiplication and Addition
// ============================================================================

namespace internal {

// In-place bitwise addition (XOR in characteristic 2)
template <typename Value>
inline void AddInPlace(Value& dst, const Value& src) {
  if constexpr (std::is_same_v<Value, Gf128>) {
    dst += src;
  } else if constexpr (std::is_same_v<Value, absl::uint128> ||
                       std::is_integral_v<Value>) {
    dst ^= src;
  } else if constexpr (is_byte_container_v<Value>) {
    if (dst.size() < src.size()) {
      dst.resize(src.size(), 0);
    }
    size_t sz = std::min(dst.size(), src.size());
    uint8_t* dst_ptr = reinterpret_cast<uint8_t*>(dst.data());
    const uint8_t* src_ptr = reinterpret_cast<const uint8_t*>(src.data());
    size_t i = 0;
#if defined(__AVX2__)
    for (; i + 32 <= sz; i += 32) {
      __m256i v_dst =
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(dst_ptr + i));
      __m256i v_src =
          _mm256_loadu_si256(reinterpret_cast<const __m256i*>(src_ptr + i));
      _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst_ptr + i),
                          _mm256_xor_si256(v_dst, v_src));
    }
#elif defined(__SSE2__)
    for (; i + 16 <= sz; i += 16) {
      __m128i v_dst =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(dst_ptr + i));
      __m128i v_src =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(src_ptr + i));
      _mm_storeu_si128(reinterpret_cast<__m128i*>(dst_ptr + i),
                       _mm_xor_si128(v_dst, v_src));
    }
#endif
    for (; i + 8 <= sz; i += 8) {
      *reinterpret_cast<uint64_t*>(dst_ptr + i) ^=
          *reinterpret_cast<const uint64_t*>(src_ptr + i);
    }
    for (; i < sz; ++i) {
      dst_ptr[i] ^= src_ptr[i];
    }
  } else {
    dst = dst ^ src;
  }
}

template <typename Value>
inline Value AddValues(const Value& a, const Value& b) {
  Value res = a;
  AddInPlace(res, b);
  return res;
}

// Multiplies a value by scalar in GF(2^128) in-place with SIMD
template <typename Value>
inline Value MultiplyByScalar(const Value& v, Gf128 scalar) {
  if constexpr (std::is_same_v<Value, Gf128>) {
    return v * scalar;
  } else if constexpr (std::is_same_v<Value, absl::uint128>) {
    return (Gf128(v) * scalar).to_uint128();
  } else if constexpr (std::is_same_v<Value, uint64_t>) {
    return (Gf128(v) * scalar).low64();
  } else if constexpr (is_byte_container_v<Value>) {
    if (v.empty()) return {};
    if (scalar.IsZero()) {
      Value zero_res(v.size(), 0);
      return zero_res;
    }
    if (scalar.IsOne()) return v;

    size_t orig_size = v.size();
    size_t num_chunks = (orig_size + 15) / 16;
    size_t padded_size = num_chunks * 16;
    std::vector<uint8_t> result(padded_size, 0);
    std::memcpy(result.data(), v.data(), orig_size);

    uint8_t* ptr = result.data();
    for (size_t c = 0; c < num_chunks; ++c) {
      __m128i word =
          _mm_loadu_si128(reinterpret_cast<const __m128i*>(ptr + c * 16));
      __m128i prod = (Gf128(word) * scalar).to_m128i();
      _mm_storeu_si128(reinterpret_cast<__m128i*>(ptr + c * 16), prod);
    }
    return result;
  } else {
    return v * scalar;
  }
}

// Generates a deterministic pseudo-random value for padding using SHA-256
// (CSPRNG)
template <typename Value>
inline Value GenerateRandomValue(uint64_t seed, size_t index,
                                 size_t byte_len = 16) {
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
    std::vector<uint8_t> v(byte_len, 0);
    size_t bytes_done = 0;
    uint32_t ctr = 0;
    while (bytes_done < byte_len) {
      uint8_t d[SHA256_DIGEST_LENGTH];
      SHA256_CTX c;
      SHA256_Init(&c);
      SHA256_Update(&c, &seed, sizeof(seed));
      SHA256_Update(&c, &index, sizeof(index));
      SHA256_Update(&c, &ctr, sizeof(ctr));
      SHA256_Final(d, &c);
      size_t to_copy = std::min<size_t>(sizeof(d), byte_len - bytes_done);
      std::memcpy(v.data() + bytes_done, d, to_copy);
      bytes_done += to_copy;
      ++ctr;
    }
    return v;
  } else if constexpr (std::is_integral_v<Value>) {
    Value res{};
    std::memcpy(&res, digest, std::min(sizeof(Value), sizeof(digest)));
    return res;
  } else {
    return Value{};
  }
}

}  // namespace internal

// ============================================================================
// Polynomial OKVS (Lagrange Interpolation over GF(2^128))
// ============================================================================
//
// Encodes n key-value pairs {(k_1, v_1), ..., (k_n, v_n)} into a polynomial
// P(X) = \sum_{j=0}^{m-1} C_j X^j of degree < m over GF(2^128).
//
// Properties:
//  - Table size: Exactly m = fixed_size coefficients (zero storage expansion!).
//  - Deterministic: Zero encoding failures (Lagrange interpolation always
//    succeeds over distinct field evaluation points).
//  - Fast: Horner's polynomial evaluation with carryless multiplication
//    (PCLMULQDQ intrinsics).
template <typename Key = std::string, typename Value = std::vector<uint8_t>>
class PolynomialOkvs : public Okvs<Key, Value> {
 public:
  // Constructs with an optional fixed table size (e.g. Z for Ring ORAM
  // buckets), PRNG seed, and default value length.
  explicit PolynomialOkvs(size_t fixed_size = 0, uint64_t seed = 0,
                          size_t default_val_len = 16)
      : fixed_size_(fixed_size),
        seed_(seed),
        default_val_len_(default_val_len) {}

  ~PolynomialOkvs() override = default;

  // Hashes a key to a unique non-zero point in GF(2^128).
  static Gf128 KeyToPoint(const Key& key, uint64_t seed) {
    if constexpr (std::is_same_v<Key, Gf128>) {
      return key.IsZero() ? Gf128(1ULL) : key;
    } else {
      std::string key_str;
      absl::string_view key_view;
      if constexpr (std::is_same_v<Key, std::string> ||
                    std::is_same_v<Key, absl::string_view>) {
        key_view = key;
      } else if constexpr (std::is_trivially_copyable_v<Key>) {
        key_view =
            absl::string_view(reinterpret_cast<const char*>(&key), sizeof(Key));
      } else {
        size_t h = absl::Hash<Key>{}(key);
        key_str.assign(reinterpret_cast<const char*>(&h), sizeof(h));
        key_view = key_str;
      }

      SHA256_CTX ctx;
      SHA256_Init(&ctx);
      SHA256_Update(&ctx, &seed, sizeof(seed));
      SHA256_Update(&ctx, key_view.data(), key_view.size());
      uint8_t hash[SHA256_DIGEST_LENGTH];
      SHA256_Final(hash, &ctx);

      Gf128 pt(absl::string_view(reinterpret_cast<const char*>(hash), 16));
      return pt.IsZero() ? Gf128(1ULL) : pt;
    }
  }

  // Encodes keys and values into polynomial coefficients table.
  absl::Status Encode(absl::Span<const Key> keys,
                      absl::Span<const Value> values) override {
    if (keys.size() != values.size()) {
      return absl::InvalidArgumentError("Keys and values size mismatch");
    }

    num_keys_ = keys.size();
    size_t target_size = (fixed_size_ > 0) ? fixed_size_ : num_keys_;
    if (target_size == 0) {
      table_.clear();
      return absl::OkStatus();
    }

    // Determine value element length (if vector)
    size_t val_len = default_val_len_;
    for (const auto& v : values) {
      if constexpr (std::is_same_v<Value, std::vector<uint8_t>>) {
        if (!v.empty()) {
          val_len = v.size();
          break;
        }
      }
    }
    size_t padded_val_len = ((val_len + 15) / 16) * 16;

    // 1. Collect points (x_i, y_i)
    std::vector<Gf128> xs;
    std::vector<Value> ys;
    xs.reserve(target_size);
    ys.reserve(target_size);

    for (size_t i = 0; i < num_keys_; ++i) {
      Gf128 x = KeyToPoint(keys[i], seed_);
      xs.push_back(x);
      if constexpr (std::is_same_v<Value, std::vector<uint8_t>>) {
        std::vector<uint8_t> y_padded(padded_val_len, 0);
        std::memcpy(y_padded.data(), values[i].data(),
                    std::min(values[i].size(), padded_val_len));
        ys.push_back(std::move(y_padded));
      } else {
        ys.push_back(values[i]);
      }
    }

    // 2. Pad up to target_size with deterministic pseudo-random points if
    // needed
    for (size_t i = num_keys_; i < target_size; ++i) {
      uint64_t dummy_id = (i + 1) * 0x9e3779b97f4a7c15ULL;
      Gf128 x = KeyToPoint(std::to_string(dummy_id), seed_ ^ 0xfeedbeefULL);
      // Ensure distinct x coordinates
      while (std::find(xs.begin(), xs.end(), x) != xs.end()) {
        x += Gf128(1ULL);
      }
      xs.push_back(x);
      ys.push_back(
          internal::GenerateRandomValue<Value>(seed_, i, padded_val_len));
    }

    // 3. Lagrange Interpolation over GF(2^128)
    table_ = InterpolatePolynomial(xs, ys);
    return absl::OkStatus();
  }

  // Decodes value associated with `key` by evaluating polynomial at x =
  // KeyToPoint(key).
  absl::StatusOr<Value> Decode(const Key& key) const override {
    if (table_.empty()) {
      return Value{};
    }
    Gf128 x = KeyToPoint(key, seed_);
    return EvaluatePolynomial(x, table_);
  }

  // Decodes from an external table of polynomial coefficients using Horner's
  // rule.
  static absl::StatusOr<Value> DecodeWithTable(const Key& key,
                                               const std::vector<Value>& table,
                                               uint64_t seed = 0) {
    if (table.empty()) return Value{};
    Gf128 x = KeyToPoint(key, seed);
    return EvaluatePolynomial(x, table);
  }

  // Direct evaluation at field point x.
  static Value EvaluatePolynomial(Gf128 x, const std::vector<Value>& table) {
    if (table.empty()) return Value{};
    if (table.size() == 1) return table[0];

    Value acc = table.back();
    if constexpr (is_byte_container_v<Value>) {
      size_t sz = acc.size();
      size_t num_chunks = (sz + 15) / 16;
      uint8_t* acc_ptr = reinterpret_cast<uint8_t*>(acc.data());

      for (int k = static_cast<int>(table.size()) - 2; k >= 0; --k) {
        const uint8_t* coeff_ptr =
            reinterpret_cast<const uint8_t*>(table[k].data());
        for (size_t c = 0; c < num_chunks; ++c) {
          __m128i word = _mm_loadu_si128(
              reinterpret_cast<const __m128i*>(acc_ptr + c * 16));
          __m128i prod = (Gf128(word) * x).to_m128i();
          __m128i coeff = _mm_loadu_si128(
              reinterpret_cast<const __m128i*>(coeff_ptr + c * 16));
          _mm_storeu_si128(reinterpret_cast<__m128i*>(acc_ptr + c * 16),
                           _mm_xor_si128(prod, coeff));
        }
      }
      return acc;
    } else {
      // Horner's method: P(x) = C_0 + x * (C_1 + x * (C_2 + ...))
      for (int k = static_cast<int>(table.size()) - 2; k >= 0; --k) {
        acc = internal::AddValues(table[k], internal::MultiplyByScalar(acc, x));
      }
      return acc;
    }
  }

  // Table accessors
  size_t TableSize() const override { return table_.size(); }
  size_t NumKeys() const override { return num_keys_; }
  const std::vector<Value>& Table() const override { return table_; }
  absl::Span<const Value> TableSpan() const override { return table_; }
  void SetTable(std::vector<Value> table) override {
    table_ = std::move(table);
  }

  uint64_t seed() const { return seed_; }
  void set_seed(uint64_t seed) { seed_ = seed; }

 private:
  // Lagrange polynomial interpolation in GF(2^128)
  static std::vector<Value> InterpolatePolynomial(absl::Span<const Gf128> xs,
                                                  absl::Span<const Value> ys) {
    size_t m = xs.size();
    if (m == 0) return {};
    if (m == 1) return {ys[0]};

    // 1. Compute master polynomial A(X) = \prod_{j=0}^{m-1} (X + x_j)
    // In GF(2^128), (X - x_j) == (X + x_j) == (X ^ x_j)
    std::vector<Gf128> A = {xs[0], Gf128::One()};  // (X + xs[0])
    for (size_t j = 1; j < m; ++j) {
      std::vector<Gf128> new_A(A.size() + 1, Gf128::Zero());
      Gf128 xj = xs[j];
      // new_A = A * (X + xj) = A * X + A * xj
      for (size_t k = 0; k < A.size(); ++k) {
        new_A[k] += A[k] * xj;
        new_A[k + 1] += A[k];
      }
      A = std::move(new_A);
    }

    // 2. For each point i, compute basis polynomial \ell_i(X) and accumulate
    // coefficients
    std::vector<Value> coefficients(m, Value{});
    if constexpr (std::is_same_v<Value, std::vector<uint8_t>>) {
      size_t byte_len = ys[0].size();
      for (auto& c : coefficients) {
        c.assign(byte_len, 0);
      }
    }

    for (size_t i = 0; i < m; ++i) {
      Gf128 xi = xs[i];

      // Synthetic division: A_i(X) = A(X) / (X + xi)
      std::vector<Gf128> b(m, Gf128::Zero());
      b[m - 1] = A[m];  // Always 1
      for (int k = static_cast<int>(m) - 2; k >= 0; --k) {
        b[k] = A[k + 1] + b[k + 1] * xi;
      }

      // Denominator D_i = A_i(x_i) = \prod_{j != i} (x_i + x_j)
      Gf128 Di = Gf128::Zero();
      for (int k = static_cast<int>(m) - 1; k >= 0; --k) {
        Di = b[k] + Di * xi;
      }

      Gf128 wi = Di.Inverse();

      // Accumulate: coefficients[k] += (wi * b[k]) * ys[i]
      for (size_t k = 0; k < m; ++k) {
        Gf128 weight = wi * b[k];
        if (weight.IsZero()) continue;

        if constexpr (is_byte_container_v<Value>) {
          size_t sz = coefficients[k].size();
          size_t num_chunks = (sz + 15) / 16;
          uint8_t* dst_ptr = reinterpret_cast<uint8_t*>(coefficients[k].data());
          const uint8_t* src_ptr =
              reinterpret_cast<const uint8_t*>(ys[i].data());
          for (size_t c = 0; c < num_chunks; ++c) {
            __m128i word = _mm_loadu_si128(
                reinterpret_cast<const __m128i*>(src_ptr + c * 16));
            __m128i prod = (Gf128(word) * weight).to_m128i();
            __m128i existing = _mm_loadu_si128(
                reinterpret_cast<const __m128i*>(dst_ptr + c * 16));
            _mm_storeu_si128(reinterpret_cast<__m128i*>(dst_ptr + c * 16),
                             _mm_xor_si128(existing, prod));
          }
        } else {
          Value term = internal::MultiplyByScalar(ys[i], weight);
          internal::AddInPlace(coefficients[k], term);
        }
      }
    }

    return coefficients;
  }

  size_t fixed_size_ = 0;
  uint64_t seed_ = 0;
  size_t default_val_len_ = 16;
  size_t num_keys_ = 0;
  std::vector<Value> table_;
};

}  // namespace oram_okvs::okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_OKVS_POLYNOMIAL_OKVS_H_
