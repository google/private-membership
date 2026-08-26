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

#include "research/oram_okvs/okvs/gf128.h"

#include <cstring>
#include <string>

#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/numbers.h"
#include "absl/strings/string_view.h"

namespace oram_okvs::okvs {

namespace {

// Repeated squaring helper for Itoh-Tsujii inversion.
inline Gf128 SquareN(Gf128 x, int n) {
  for (int i = 0; i < n; ++i) {
    x = x.Square();
  }
  return x;
}

}  // namespace

Gf128::Gf128(absl::string_view bytes) {
  if (bytes.size() >= 16) {
    value_ = _mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes.data()));
  } else {
    std::array<uint8_t, 16> buf = {0};
    std::memcpy(buf.data(), bytes.data(), bytes.size());
    value_ = _mm_loadu_si128(reinterpret_cast<const __m128i*>(buf.data()));
  }
}

absl::StatusOr<Gf128> Gf128::FromHexString(absl::string_view hex) {
  if (hex.size() != 32) {
    return absl::InvalidArgumentError(
        "Hex string for Gf128 must be exactly 32 hex digits (128 bits).");
  }
  uint64_t high = 0;
  uint64_t low = 0;
  if (!absl::SimpleHexAtoi(hex.substr(0, 16), &high) ||
      !absl::SimpleHexAtoi(hex.substr(16, 16), &low)) {
    return absl::InvalidArgumentError(
        "Invalid hexadecimal character in string.");
  }
  return Gf128(high, low);
}

absl::StatusOr<Gf128> Gf128::FromBytes(absl::string_view bytes) {
  if (bytes.size() != 16) {
    return absl::InvalidArgumentError(
        "Byte string for Gf128 must be exactly 16 bytes.");
  }
  return Gf128(bytes);
}

Gf128 Gf128::Inverse() const {
  if (IsZero()) {
    return Gf128::Zero();
  }

  // Constant-time inversion via Itoh-Tsujii addition chain for 2^128 - 2.
  // x^(2^128 - 2) = (x^(2^127 - 1))^2
  // 127 = 1 + 2 + 4 + 8 + 16 + 32 + 64 (1111111 in binary)
  // Let beta_k = x^(2^k - 1)
  // beta_1 = x
  Gf128 b1 = *this;

  // beta_2 = (beta_1)^2 * beta_1 = x^3
  Gf128 b2 = b1.Square() * b1;

  // beta_3 = (beta_2)^2 * beta_1 = x^7
  Gf128 b3 = b2.Square() * b1;

  // beta_6 = (beta_3)^(2^3) * beta_3 = x^(2^6 - 1) = x^63
  Gf128 b6 = SquareN(b3, 3) * b3;

  // beta_7 = (beta_6)^2 * beta_1 = x^(2^7 - 1) = x^127
  Gf128 b7 = b6.Square() * b1;

  // beta_14 = (beta_7)^(2^7) * beta_7 = x^(2^14 - 1)
  Gf128 b14 = SquareN(b7, 7) * b7;

  // beta_15 = (beta_14)^2 * beta_1 = x^(2^15 - 1)
  Gf128 b15 = b14.Square() * b1;

  // beta_30 = (beta_15)^(2^15) * beta_15 = x^(2^30 - 1)
  Gf128 b30 = SquareN(b15, 15) * b15;

  // beta_31 = (beta_30)^2 * beta_1 = x^(2^31 - 1)
  Gf128 b31 = b30.Square() * b1;

  // beta_62 = (beta_31)^(2^31) * beta_31 = x^(2^62 - 1)
  Gf128 b62 = SquareN(b31, 31) * b31;

  // beta_63 = (beta_62)^2 * beta_1 = x^(2^63 - 1)
  Gf128 b63 = b62.Square() * b1;

  // beta_126 = (beta_63)^(2^63) * beta_63 = x^(2^126 - 1)
  Gf128 b126 = SquareN(b63, 63) * b63;

  // beta_127 = (beta_126)^2 * beta_1 = x^(2^127 - 1)
  Gf128 b127 = b126.Square() * b1;

  // x^(2^128 - 2) = (x^(2^127 - 1))^2
  return b127.Square();
}

Gf128 Gf128::Pow(absl::uint128 exp) const {
  Gf128 result = Gf128::One();
  Gf128 base = *this;
  while (exp > 0) {
    if (absl::Uint128Low64(exp) & 1) {
      result *= base;
    }
    base = base.Square();
    exp >>= 1;
  }
  return result;
}

}  // namespace oram_okvs::okvs
