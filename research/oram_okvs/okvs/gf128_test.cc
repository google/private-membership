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

#include <random>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/hash/hash.h"
#include "absl/numeric/int128.h"

namespace oram_okvs::okvs {
namespace {

TEST(Gf128Test, ZeroAndOne) {
  Gf128 zero = Gf128::Zero();
  EXPECT_TRUE(zero.IsZero());
  EXPECT_FALSE(zero.IsOne());
  EXPECT_EQ(zero.high64(), 0);
  EXPECT_EQ(zero.low64(), 0);

  Gf128 one = Gf128::One();
  EXPECT_FALSE(one.IsZero());
  EXPECT_TRUE(one.IsOne());
  EXPECT_EQ(one.high64(), 0);
  EXPECT_EQ(one.low64(), 1);

  EXPECT_NE(zero, one);
}

TEST(Gf128Test, ConstructorsAndConversions) {
  Gf128 a(0x1234567890ABCDEFULL, 0xFEDCBA0987654321ULL);
  EXPECT_EQ(a.high64(), 0x1234567890ABCDEFULL);
  EXPECT_EQ(a.low64(), 0xFEDCBA0987654321ULL);

  absl::uint128 u = a.to_uint128();
  EXPECT_EQ(absl::Uint128High64(u), 0x1234567890ABCDEFULL);
  EXPECT_EQ(absl::Uint128Low64(u), 0xFEDCBA0987654321ULL);

  Gf128 from_u(u);
  EXPECT_EQ(a, from_u);

  std::string hex = a.ToHexString();
  EXPECT_EQ(hex, "1234567890abcdeffedcba0987654321");

  auto parsed = Gf128::FromHexString(hex);
  ASSERT_TRUE(parsed.ok());
  EXPECT_EQ(*parsed, a);

  std::string bytes = a.ToBytes();
  EXPECT_EQ(bytes.size(), 16);
  auto from_bytes = Gf128::FromBytes(bytes);
  ASSERT_TRUE(from_bytes.ok());
  EXPECT_EQ(*from_bytes, a);
}

TEST(Gf128Test, AdditionAndSubtraction) {
  Gf128 a(0xAAAAAAAAAAAAAAAAULL, 0x5555555555555555ULL);
  Gf128 b(0x5555555555555555ULL, 0xAAAAAAAAAAAAAAAAULL);

  Gf128 sum = a + b;
  EXPECT_EQ(sum.high64(), 0xFFFFFFFFFFFFFFFFULL);
  EXPECT_EQ(sum.low64(), 0xFFFFFFFFFFFFFFFFULL);

  // In GF(2^128), a + b == a - b == a ^ b
  EXPECT_EQ(a - b, sum);
  EXPECT_EQ(a ^ b, sum);

  // a + a == 0
  EXPECT_TRUE((a + a).IsZero());
  EXPECT_TRUE((a - a).IsZero());
  EXPECT_TRUE((a ^ a).IsZero());

  // a + 0 == a
  EXPECT_EQ(a + Gf128::Zero(), a);
}

TEST(Gf128Test, MultiplicationIdentityAndZero) {
  Gf128 a(0x12345678ULL, 0x9ABCDEF0ULL);
  EXPECT_TRUE((a * Gf128::Zero()).IsZero());
  EXPECT_EQ(a * Gf128::One(), a);
  EXPECT_EQ(Gf128::One() * a, a);
}

TEST(Gf128Test, ReductionPolynomial) {
  // In GF(2^128) with P(x) = x^128 + x^7 + x^2 + x + 1 (modulus 0x87):
  // x^127 * x = x^128 = 0x87.
  Gf128 x_127(1ULL << 63, 0);  // high bit 127 is set
  Gf128 x(0, 2ULL);            // x^1 (bit 1)

  Gf128 prod = x_127 * x;
  EXPECT_EQ(prod.high64(), 0);
  EXPECT_EQ(prod.low64(), 0x87ULL);
}

TEST(Gf128Test, MultiplicationProperties) {
  std::mt19937_64 rng(42);
  for (int i = 0; i < 50; ++i) {
    Gf128 a = Gf128::Random(rng);
    Gf128 b = Gf128::Random(rng);
    Gf128 c = Gf128::Random(rng);

    // Commutativity: a * b == b * a
    EXPECT_EQ(a * b, b * a);

    // Distributivity: a * (b + c) == (a * b) + (a * c)
    EXPECT_EQ(a * (b + c), (a * b) + (a * c));

    // Associativity: (a * b) * c == a * (b * c)
    EXPECT_EQ((a * b) * c, a * (b * c));

    // Squaring: a.Square() == a * a
    EXPECT_EQ(a.Square(), a * a);
  }
}

TEST(Gf128Test, InversionAndDivision) {
  std::mt19937_64 rng(12345);

  // Inversion of 1 is 1
  EXPECT_EQ(Gf128::One().Inverse(), Gf128::One());
  // Inversion of 0 is 0
  EXPECT_EQ(Gf128::Zero().Inverse(), Gf128::Zero());

  for (int i = 0; i < 50; ++i) {
    Gf128 a = Gf128::Random(rng);
    if (a.IsZero()) continue;

    Gf128 inv = a.Inverse();
    EXPECT_EQ(a * inv, Gf128::One());
    EXPECT_EQ(inv * a, Gf128::One());

    // Double inversion: (a^-1)^-1 == a
    EXPECT_EQ(inv.Inverse(), a);

    Gf128 b = Gf128::Random(rng);
    if (b.IsZero()) continue;

    Gf128 div = a / b;
    EXPECT_EQ(div * b, a);
  }
}

TEST(Gf128Test, Exponentiation) {
  std::mt19937_64 rng(999);
  Gf128 a = Gf128::Random(rng);

  EXPECT_EQ(a.Pow(0), Gf128::One());
  EXPECT_EQ(a.Pow(1), a);
  EXPECT_EQ(a.Pow(2), a.Square());
  EXPECT_EQ(a.Pow(3), a.Square() * a);
  EXPECT_EQ(a.Pow(4), a.Square().Square());
}

TEST(Gf128Test, ComparisonAndHashing) {
  Gf128 a(1, 2);
  Gf128 b(1, 2);
  Gf128 c(1, 3);
  Gf128 d(2, 1);

  EXPECT_EQ(a, b);
  EXPECT_NE(a, c);
  EXPECT_TRUE(a < c);
  EXPECT_TRUE(c < d);

  EXPECT_EQ(absl::Hash<Gf128>{}(a), absl::Hash<Gf128>{}(b));
}

}  // namespace
}  // namespace oram_okvs::okvs
