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

#include "research/oram_okvs/okvs/random_band_okvs.h"

#include <array>
#include <cstdint>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/strings/str_cat.h"
#include "research/oram_okvs/okvs/gf128.h"

namespace oram_okvs::okvs {
namespace {

TEST(RandomBandOkvsTest, EmptySet) {
  RandomBandOkvs<std::string, Gf128> okvs;
  EXPECT_TRUE(okvs.Encode({}, {}).ok());
  EXPECT_EQ(okvs.NumKeys(), 0);
  EXPECT_EQ(okvs.TableSize(), 0);
}

TEST(RandomBandOkvsTest, SingleElement) {
  RandomBandOkvs<std::string, Gf128> okvs;
  std::vector<std::string> keys = {"hello"};
  std::vector<Gf128> values = {Gf128(0x12345678, 0x9ABCDEF0)};

  ASSERT_TRUE(okvs.Encode(keys, values).ok());
  EXPECT_EQ(okvs.NumKeys(), 1);
  EXPECT_GT(okvs.TableSize(), 0);

  auto decoded = okvs.Decode("hello");
  ASSERT_TRUE(decoded.ok());
  EXPECT_EQ(*decoded, values[0]);
}

TEST(RandomBandOkvsTest, SmallKeySetGf128) {
  constexpr int kNumKeys = 100;
  std::mt19937_64 rng(1234);

  std::vector<std::string> keys;
  std::vector<Gf128> values;
  keys.reserve(kNumKeys);
  values.reserve(kNumKeys);

  for (int i = 0; i < kNumKeys; ++i) {
    keys.push_back(absl::StrCat("key_", i));
    values.push_back(Gf128::Random(rng));
  }

  RandomBandOkvs<std::string, Gf128> okvs(/*epsilon=*/0.15, /*seed=*/42);
  ASSERT_TRUE(okvs.Encode(keys, values).ok());
  EXPECT_EQ(okvs.NumKeys(), kNumKeys);

  for (int i = 0; i < kNumKeys; ++i) {
    auto decoded = okvs.Decode(keys[i]);
    ASSERT_TRUE(decoded.ok());
    EXPECT_EQ(*decoded, values[i]);
  }
}

TEST(RandomBandOkvsTest, MediumKeySetGf128) {
  constexpr int kNumKeys = 1000;
  std::mt19937_64 rng(5678);

  std::vector<std::string> keys;
  std::vector<Gf128> values;
  keys.reserve(kNumKeys);
  values.reserve(kNumKeys);

  for (int i = 0; i < kNumKeys; ++i) {
    keys.push_back(absl::StrCat("key_", rng()));
    values.push_back(Gf128::Random(rng));
  }

  RandomBandOkvs<std::string, Gf128> okvs(/*epsilon=*/0.10, /*seed=*/100);
  ASSERT_TRUE(okvs.Encode(keys, values).ok());

  for (int i = 0; i < kNumKeys; ++i) {
    auto decoded = okvs.Decode(keys[i]);
    ASSERT_TRUE(decoded.ok());
    EXPECT_EQ(*decoded, values[i]);
  }
}

TEST(RandomBandOkvsTest, NumericKeysAndValues) {
  constexpr int kNumKeys = 200;
  std::mt19937_64 rng(999);

  std::vector<uint64_t> keys;
  std::vector<uint64_t> values;
  for (int i = 0; i < kNumKeys; ++i) {
    keys.push_back(rng());
    values.push_back(rng());
  }

  RandomBandOkvs<uint64_t, uint64_t> okvs(/*epsilon=*/0.12, /*seed=*/1);
  ASSERT_TRUE(okvs.Encode(keys, values).ok());

  for (int i = 0; i < kNumKeys; ++i) {
    auto decoded = okvs.Decode(keys[i]);
    ASSERT_TRUE(decoded.ok());
    EXPECT_EQ(*decoded, values[i]);
  }
}

TEST(RandomBandOkvsTest, Uint128KeysAndValues) {
  constexpr int kNumKeys = 150;
  std::mt19937_64 rng(777);

  std::vector<absl::uint128> keys;
  std::vector<absl::uint128> values;
  for (int i = 0; i < kNumKeys; ++i) {
    keys.push_back(absl::MakeUint128(rng(), rng()));
    values.push_back(absl::MakeUint128(rng(), rng()));
  }

  RandomBandOkvs<absl::uint128, absl::uint128> okvs;
  ASSERT_TRUE(okvs.Encode(keys, values).ok());

  for (int i = 0; i < kNumKeys; ++i) {
    auto decoded = okvs.Decode(keys[i]);
    ASSERT_TRUE(decoded.ok());
    EXPECT_EQ(*decoded, values[i]);
  }
}

TEST(RandomBandOkvsTest, Gf128KeysAndValues) {
  constexpr int kNumKeys = 150;
  std::mt19937_64 rng(333);

  std::vector<Gf128> keys;
  std::vector<Gf128> values;
  for (int i = 0; i < kNumKeys; ++i) {
    keys.push_back(Gf128::Random(rng));
    values.push_back(Gf128::Random(rng));
  }

  RandomBandOkvs<Gf128, Gf128> okvs;
  ASSERT_TRUE(okvs.Encode(keys, values).ok());

  for (int i = 0; i < kNumKeys; ++i) {
    auto decoded = okvs.Decode(keys[i]);
    ASSERT_TRUE(decoded.ok());
    EXPECT_EQ(*decoded, values[i]);
  }
}

TEST(RandomBandOkvsTest, ByteArrayValue) {
  constexpr int kNumKeys = 50;
  std::mt19937_64 rng(111);

  std::vector<std::string> keys;
  std::vector<std::array<uint8_t, 32>> values;

  for (int i = 0; i < kNumKeys; ++i) {
    keys.push_back(absl::StrCat("key_", i));
    std::array<uint8_t, 32> val;
    for (size_t b = 0; b < 32; ++b) {
      val[b] = static_cast<uint8_t>(rng());
    }
    values.push_back(val);
  }

  RandomBandOkvs<std::string, std::array<uint8_t, 32>> okvs;
  ASSERT_TRUE(okvs.Encode(keys, values).ok());

  for (int i = 0; i < kNumKeys; ++i) {
    auto decoded = okvs.Decode(keys[i]);
    ASSERT_TRUE(decoded.ok());
    EXPECT_EQ(*decoded, values[i]);
  }
}

TEST(RandomBandOkvsTest, DynamicByteVectorValue) {
  constexpr int kNumKeys = 50;
  constexpr size_t kBlockSize = 84;  // Arbitrary non-standard size
  std::mt19937_64 rng(222);

  std::vector<std::string> keys;
  std::vector<std::vector<uint8_t>> values;

  for (int i = 0; i < kNumKeys; ++i) {
    keys.push_back(absl::StrCat("key_", i));
    std::vector<uint8_t> val(kBlockSize);
    for (size_t b = 0; b < kBlockSize; ++b) {
      val[b] = static_cast<uint8_t>(rng());
    }
    values.push_back(val);
  }

  RandomBandOkvs<std::string, std::vector<uint8_t>> okvs;
  ASSERT_TRUE(okvs.Encode(keys, values).ok());

  for (int i = 0; i < kNumKeys; ++i) {
    auto decoded = okvs.Decode(keys[i]);
    ASSERT_TRUE(decoded.ok());
    EXPECT_EQ(*decoded, values[i]);
  }
}

TEST(RandomBandOkvsTest, DifferentBandBitWidths) {
  constexpr int kNumKeys = 100;
  std::mt19937_64 rng(444);

  std::vector<std::string> keys;
  std::vector<Gf128> values;
  for (int i = 0; i < kNumKeys; ++i) {
    keys.push_back(absl::StrCat("k_", i));
    values.push_back(Gf128::Random(rng));
  }

  // Band width 128
  RandomBandOkvs<std::string, Gf128, 128> okvs_128(/*epsilon=*/0.20,
                                                   /*seed=*/1);
  ASSERT_TRUE(okvs_128.Encode(keys, values).ok());
  for (int i = 0; i < kNumKeys; ++i) {
    EXPECT_EQ(*okvs_128.Decode(keys[i]), values[i]);
  }

  // Band width 256
  RandomBandOkvs<std::string, Gf128, 256> okvs_256(/*epsilon=*/0.15,
                                                   /*seed=*/2);
  ASSERT_TRUE(okvs_256.Encode(keys, values).ok());
  for (int i = 0; i < kNumKeys; ++i) {
    EXPECT_EQ(*okvs_256.Decode(keys[i]), values[i]);
  }

  // Band width 512
  RandomBandOkvs<std::string, Gf128, 512> okvs_512(/*epsilon=*/0.10,
                                                   /*seed=*/3);
  ASSERT_TRUE(okvs_512.Encode(keys, values).ok());
  for (int i = 0; i < kNumKeys; ++i) {
    EXPECT_EQ(*okvs_512.Decode(keys[i]), values[i]);
  }
}

TEST(RandomBandOkvsTest, EncodePairsOverload) {
  std::vector<std::pair<std::string, Gf128>> pairs = {
      {"apple", Gf128(1, 2)},
      {"banana", Gf128(3, 4)},
      {"cherry", Gf128(5, 6)},
  };

  RandomBandOkvs<std::string, Gf128> okvs;
  ASSERT_TRUE(okvs.Encode(pairs).ok());

  for (const auto& [k, v] : pairs) {
    auto decoded = okvs.Decode(k);
    ASSERT_TRUE(decoded.ok());
    EXPECT_EQ(*decoded, v);
  }
}

TEST(RandomBandOkvsTest, ObliviousDecodingNonExistentKeys) {
  std::vector<std::string> keys = {"alpha", "beta", "gamma"};
  std::vector<Gf128> values = {Gf128(10), Gf128(20), Gf128(30)};

  RandomBandOkvs<std::string, Gf128> okvs;
  ASSERT_TRUE(okvs.Encode(keys, values).ok());

  // Non-existent key must decode to a deterministic value without crashing.
  auto val1 = okvs.Decode("delta");
  auto val2 = okvs.Decode("delta");
  ASSERT_TRUE(val1.ok());
  ASSERT_TRUE(val2.ok());
  EXPECT_EQ(*val1, *val2);

  // Different non-existent key gets another pseudo-random value.
  auto val3 = okvs.Decode("epsilon");
  ASSERT_TRUE(val3.ok());
  EXPECT_NE(*val1, *val3);
}

TEST(RandomBandOkvsTest, TableSerializationAndRemoteDecoding) {
  std::vector<std::string> keys = {"k1", "k2", "k3"};
  std::vector<Gf128> values = {Gf128(1, 1), Gf128(2, 2), Gf128(3, 3)};

  RandomBandOkvs<std::string, Gf128> encoder(/*epsilon=*/0.15, /*seed=*/123);
  ASSERT_TRUE(encoder.Encode(keys, values).ok());

  // Extract table and pass to a new OKVS instance.
  std::vector<Gf128> table = encoder.Table();
  RandomBandOkvs<std::string, Gf128> decoder(table, /*seed=*/123);

  for (size_t i = 0; i < keys.size(); ++i) {
    auto decoded = decoder.Decode(keys[i]);
    ASSERT_TRUE(decoded.ok());
    EXPECT_EQ(*decoded, values[i]);
  }

  // Also test static DecodeWithTable.
  for (size_t i = 0; i < keys.size(); ++i) {
    auto decoded = RandomBandOkvs<std::string, Gf128>::DecodeWithTable(
        keys[i], table, /*seed=*/123);
    ASSERT_TRUE(decoded.ok());
    EXPECT_EQ(*decoded, values[i]);
  }
}

TEST(RandomBandOkvsTest, ErrorHandling) {
  RandomBandOkvs<std::string, Gf128> okvs;

  // Decoding before encoding should fail.
  EXPECT_FALSE(okvs.Decode("test").ok());

  // Mismatched keys and values sizes.
  std::vector<std::string> keys = {"a", "b"};
  std::vector<Gf128> values = {Gf128(1)};
  EXPECT_FALSE(okvs.Encode(keys, values).ok());
}

}  // namespace
}  // namespace oram_okvs::okvs
