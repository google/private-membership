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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_OKVS_GF128_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_OKVS_GF128_H_

#include <emmintrin.h>
#include <smmintrin.h>
#include <stdint.h>
#include <wmmintrin.h>

#include <array>
#include <ostream>
#include <random>
#include <string>
#include <utility>

#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"

namespace oram_okvs::okvs {

// GF(2^128) represents elements of the Galois Field GF(2^128) with the
// standard AES-GCM / cryptographic irreducible reduction polynomial:
//   P(x) = x^128 + x^7 + x^2 + x + 1 (low-order coefficient byte 0x87).
//
// Arithmetic:
// - Addition / Subtraction: Bitwise XOR (characteristic 2).
// - Multiplication: Carryless polynomial multiplication with PCLMULQDQ
//   intrinsics, reduced modulo P(x).
// - Inversion: Constant-time field inversion via Fermat's Little Theorem /
//   Itoh-Tsujii addition chain (a^(2^128 - 2)).
// - Division: a * b^(-1).
class Gf128 {
 public:
  // Default constructor: creates field element 0.
  Gf128() : value_(_mm_setzero_si128()) {}

  // Constructs from a 64-bit unsigned integer (low 64 bits).
  explicit Gf128(uint64_t low)
      : value_(_mm_set_epi64x(0, static_cast<int64_t>(low))) {}

  // Constructs from high and low 64-bit words.
  Gf128(uint64_t high, uint64_t low)
      : value_(_mm_set_epi64x(static_cast<int64_t>(high),
                              static_cast<int64_t>(low))) {}

  // Constructs from __m128i intrinsic register.
  explicit Gf128(__m128i val) : value_(val) {}

  // Constructs from absl::uint128.
  explicit Gf128(absl::uint128 val)
      : Gf128(absl::Uint128High64(val), absl::Uint128Low64(val)) {}

  // Constructs from raw byte array or string_view (at least 16 bytes, or
  // 0-padded).
  explicit Gf128(absl::string_view bytes);

  // Constructs from 16-byte std::array.
  explicit Gf128(const std::array<uint8_t, 16>& bytes) {
    value_ = _mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes.data()));
  }

  // Factory methods.
  static Gf128 Zero() { return Gf128(); }
  static Gf128 One() { return Gf128(1ULL); }

  // Parses a 32-character hexadecimal string into a Gf128 element.
  static absl::StatusOr<Gf128> FromHexString(absl::string_view hex);

  // Creates from a 16-byte raw string.
  static absl::StatusOr<Gf128> FromBytes(absl::string_view bytes);

  // Generates a uniform random field element using a given PRNG.
  template <typename URBG>
  static Gf128 Random(URBG& gen) {
    std::uniform_int_distribution<uint64_t> dist;
    uint64_t low = dist(gen);
    uint64_t high = dist(gen);
    return Gf128(high, low);
  }

  // Element state checks.
  bool IsZero() const { return _mm_testz_si128(value_, value_) != 0; }

  bool IsOne() const { return *this == One(); }

  // Accessors and conversions.
  __m128i to_m128i() const { return value_; }
  explicit operator __m128i() const { return value_; }

  uint64_t low64() const {
    return static_cast<uint64_t>(_mm_cvtsi128_si64(value_));
  }

  uint64_t high64() const {
    return static_cast<uint64_t>(_mm_cvtsi128_si64(_mm_srli_si128(value_, 8)));
  }

  absl::uint128 to_uint128() const {
    return absl::MakeUint128(high64(), low64());
  }

  explicit operator absl::uint128() const { return to_uint128(); }

  std::string ToHexString() const {
    return absl::StrFormat("%016x%016x", high64(), low64());
  }

  std::string ToBytes() const {
    std::string s(16, '\0');
    _mm_storeu_si128(reinterpret_cast<__m128i*>(s.data()), value_);
    return s;
  }

  std::array<uint8_t, 16> ToByteArray() const {
    std::array<uint8_t, 16> arr;
    _mm_storeu_si128(reinterpret_cast<__m128i*>(arr.data()), value_);
    return arr;
  }

  // Field addition and subtraction (both are bitwise XOR in characteristic 2).
  Gf128 operator+(const Gf128& other) const {
    return Gf128(_mm_xor_si128(value_, other.value_));
  }

  Gf128 operator-(const Gf128& other) const {
    return Gf128(_mm_xor_si128(value_, other.value_));
  }

  Gf128 operator^(const Gf128& other) const {
    return Gf128(_mm_xor_si128(value_, other.value_));
  }

  Gf128& operator+=(const Gf128& other) {
    value_ = _mm_xor_si128(value_, other.value_);
    return *this;
  }

  Gf128& operator-=(const Gf128& other) {
    value_ = _mm_xor_si128(value_, other.value_);
    return *this;
  }

  Gf128& operator^=(const Gf128& other) {
    value_ = _mm_xor_si128(value_, other.value_);
    return *this;
  }

  // Field multiplication modulo P(x).
  Gf128 operator*(const Gf128& other) const;

  Gf128& operator*=(const Gf128& other) {
    *this = *this * other;
    return *this;
  }

  // Fast squaring in GF(2^128) using Frobenius endomorphism / no cross terms.
  Gf128 Square() const;

  // Multiplicative inverse in GF(2^128). Returns 0 for 0.
  Gf128 Inverse() const;

  // Field division: *this * other.Inverse().
  Gf128 operator/(const Gf128& other) const { return *this * other.Inverse(); }

  Gf128& operator/=(const Gf128& other) {
    *this = *this / other;
    return *this;
  }

  // Exponentiation by unsigned 128-bit integer exponent.
  Gf128 Pow(absl::uint128 exp) const;

  // Equality and ordering.
  bool operator==(const Gf128& other) const {
    __m128i neq = _mm_xor_si128(value_, other.value_);
    return _mm_testz_si128(neq, neq) != 0;
  }

  bool operator!=(const Gf128& other) const { return !(*this == other); }

  bool operator<(const Gf128& other) const {
    uint64_t h1 = high64(), h2 = other.high64();
    if (h1 != h2) return h1 < h2;
    return low64() < other.low64();
  }

  // AbslHashValue support.
  template <typename H>
  friend H AbslHashValue(H h, const Gf128& v) {
    return H::combine(std::move(h), v.high64(), v.low64());
  }

  // Stream output operator for logging and test assertions.
  friend std::ostream& operator<<(std::ostream& os, const Gf128& v) {
    return os << v.ToHexString();
  }

 private:
  __m128i value_;

  static constexpr uint64_t kModulus = 0x87ULL;
};

// ================= Inline implementations =================

inline Gf128 Gf128::operator*(const Gf128& other) const {
  const __m128i a = value_;
  const __m128i b = other.value_;
  const __m128i modulus = _mm_set_epi64x(0, static_cast<int64_t>(kModulus));

  __m128i mul256_high = _mm_clmulepi64_si128(a, b, 0x11);
  __m128i mul256_low = _mm_clmulepi64_si128(a, b, 0x00);
  __m128i mul256_mid1 = _mm_clmulepi64_si128(a, b, 0x01);
  __m128i mul256_mid2 = _mm_clmulepi64_si128(a, b, 0x10);

  __m128i mul256_mid = _mm_xor_si128(mul256_mid1, mul256_mid2);
  mul256_high = _mm_xor_si128(mul256_high, _mm_srli_si128(mul256_mid, 8));
  mul256_low = _mm_xor_si128(mul256_low, _mm_slli_si128(mul256_mid, 8));

  // Reduce w.r.t. high 64 bits of mul256_high
  __m128i tmp = _mm_clmulepi64_si128(mul256_high, modulus, 0x01);
  mul256_low = _mm_xor_si128(mul256_low, _mm_slli_si128(tmp, 8));
  mul256_high = _mm_xor_si128(mul256_high, _mm_srli_si128(tmp, 8));

  // Reduce w.r.t. low 64 bits of mul256_high
  tmp = _mm_clmulepi64_si128(mul256_high, modulus, 0x00);
  mul256_low = _mm_xor_si128(mul256_low, tmp);

  return Gf128(mul256_low);
}

inline Gf128 Gf128::Square() const {
  const __m128i a = value_;
  const __m128i modulus = _mm_set_epi64x(0, static_cast<int64_t>(kModulus));

  // In GF(2)[x], (a1*x^64 + a0)^2 = a1^2 * x^128 + a0^2 (no middle cross
  // terms).
  __m128i mul256_high = _mm_clmulepi64_si128(a, a, 0x11);
  __m128i mul256_low = _mm_clmulepi64_si128(a, a, 0x00);

  // Reduce w.r.t. high 64 bits of mul256_high
  __m128i tmp = _mm_clmulepi64_si128(mul256_high, modulus, 0x01);
  mul256_low = _mm_xor_si128(mul256_low, _mm_slli_si128(tmp, 8));
  mul256_high = _mm_xor_si128(mul256_high, _mm_srli_si128(tmp, 8));

  // Reduce w.r.t. low 64 bits of mul256_high
  tmp = _mm_clmulepi64_si128(mul256_high, modulus, 0x00);
  mul256_low = _mm_xor_si128(mul256_low, tmp);

  return Gf128(mul256_low);
}

}  // namespace oram_okvs::okvs

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_OKVS_GF128_H_
