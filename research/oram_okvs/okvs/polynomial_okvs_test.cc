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

#include "research/oram_okvs/okvs/polynomial_okvs.h"

#include <cstdint>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "absl/strings/str_format.h"
#include "research/oram_okvs/okvs/gf128.h"

namespace oram_okvs::okvs {
namespace {

TEST(PolynomialOkvsTest, BasicEncodeDecodeGf128) {
  PolynomialOkvs<std::string, Gf128> okvs(/*fixed_size=*/0, /*seed=*/42);

  std::vector<std::string> keys = {"alpha", "beta", "gamma"};
  std::vector<Gf128> values = {
      Gf128(0x1234, 0x5678),
      Gf128(0xabcd, 0xef01),
      Gf128(0x9999, 0x8888),
  };

  ASSERT_TRUE(okvs.Encode(keys, values).ok());
  EXPECT_EQ(okvs.TableSize(), 3);

  for (size_t i = 0; i < keys.size(); ++i) {
    auto dec_or = okvs.Decode(keys[i]);
    ASSERT_TRUE(dec_or.ok());
    EXPECT_EQ(*dec_or, values[i]) << "Mismatch for key: " << keys[i];
  }

  // Non-encoded key should evaluate to non-zero pseudo-random value
  auto dummy_or = okvs.Decode("unknown_key");
  ASSERT_TRUE(dummy_or.ok());
  EXPECT_NE(*dummy_or, values[0]);
  EXPECT_NE(*dummy_or, values[1]);
  EXPECT_NE(*dummy_or, values[2]);
}

TEST(PolynomialOkvsTest, BasicEncodeDecodeVectorBytes) {
  PolynomialOkvs<std::string, std::vector<uint8_t>> okvs(/*fixed_size=*/0,
                                                         /*seed=*/1234);

  std::vector<std::string> keys = {"block_0", "block_1", "block_2"};
  std::vector<std::vector<uint8_t>> values = {
      std::vector<uint8_t>(64, 0xAA),
      std::vector<uint8_t>(64, 0xBB),
      std::vector<uint8_t>(64, 0xCC),
  };
  values[0][0] = 1;
  values[1][0] = 2;
  values[2][0] = 3;

  ASSERT_TRUE(okvs.Encode(keys, values).ok());
  EXPECT_EQ(okvs.TableSize(), 3);

  for (size_t i = 0; i < keys.size(); ++i) {
    auto dec_or = okvs.Decode(keys[i]);
    ASSERT_TRUE(dec_or.ok());
    EXPECT_EQ(*dec_or, values[i]) << "Mismatch for key: " << keys[i];
  }
}

TEST(PolynomialOkvsTest, FixedSizePadding) {
  // Fixed size = 4, but encoding only 2 items
  PolynomialOkvs<std::string, std::vector<uint8_t>> okvs(/*fixed_size=*/4,
                                                         /*seed=*/777);

  std::vector<std::string> keys = {"key_a", "key_b"};
  std::vector<std::vector<uint8_t>> values = {
      std::vector<uint8_t>(32, 0x11),
      std::vector<uint8_t>(32, 0x22),
  };

  ASSERT_TRUE(okvs.Encode(keys, values).ok());
  EXPECT_EQ(okvs.TableSize(), 4);  // Exactly fixed_size

  for (size_t i = 0; i < keys.size(); ++i) {
    auto dec_or = okvs.Decode(keys[i]);
    ASSERT_TRUE(dec_or.ok());
    EXPECT_EQ(*dec_or, values[i]);
  }
}

TEST(PolynomialOkvsTest, SingleElementAndEmpty) {
  PolynomialOkvs<std::string, std::vector<uint8_t>> okvs(/*fixed_size=*/2,
                                                         /*seed=*/999);

  // Single element with fixed size 2
  std::vector<std::string> keys = {"single_key"};
  std::vector<std::vector<uint8_t>> values = {std::vector<uint8_t>(16, 0x55)};

  ASSERT_TRUE(okvs.Encode(keys, values).ok());
  EXPECT_EQ(okvs.TableSize(), 2);

  auto dec_or = okvs.Decode("single_key");
  ASSERT_TRUE(dec_or.ok());
  EXPECT_EQ(*dec_or, values[0]);

  // Empty encode with fixed size 2
  ASSERT_TRUE(okvs.Encode({}, {}).ok());
  EXPECT_EQ(okvs.TableSize(), 2);
  auto dec_empty_or = okvs.Decode("any_key");
  ASSERT_TRUE(dec_empty_or.ok());
}

TEST(PolynomialOkvsTest, DecodeWithTableStaticMethod) {
  PolynomialOkvs<std::string, std::vector<uint8_t>> okvs(/*fixed_size=*/3,
                                                         /*seed=*/42);

  std::vector<std::string> keys = {"k1", "k2"};
  std::vector<std::vector<uint8_t>> values = {
      std::vector<uint8_t>(32, 0x77),
      std::vector<uint8_t>(32, 0x88),
  };

  ASSERT_TRUE(okvs.Encode(keys, values).ok());

  auto static_val1_or =
      PolynomialOkvs<std::string, std::vector<uint8_t>>::DecodeWithTable(
          "k1", okvs.Table(), /*seed=*/42);
  ASSERT_TRUE(static_val1_or.ok());
  EXPECT_EQ(*static_val1_or, values[0]);

  auto static_val2_or =
      PolynomialOkvs<std::string, std::vector<uint8_t>>::DecodeWithTable(
          "k2", okvs.Table(), /*seed=*/42);
  ASSERT_TRUE(static_val2_or.ok());
  EXPECT_EQ(*static_val2_or, values[1]);
}

}  // namespace
}  // namespace oram_okvs::okvs
