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

#include "research/oram_okvs/oram/oblivious_container.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "research/oram_okvs/okvs/polynomial_okvs.h"
#include "research/oram_okvs/okvs/random_band_okvs.h"
#include "research/oram_okvs/oram/crypto_utils.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {
namespace {

TEST(ObliviousContainerTest, KeyedObliviousContainerBuildAndQuery) {
  std::string seed_l = GenerateRandomBytes(16);
  std::string seed_d = GenerateRandomBytes(16);
  std::string key_enc = GenerateRandomBytes(16);

  constexpr size_t kBlockSize = 32;

  std::vector<std::pair<BlockId, std::vector<uint8_t>>> real_blocks = {
      {10, std::vector<uint8_t>(kBlockSize, 0x11)},
      {20, std::vector<uint8_t>(kBlockSize, 0x22)},
      {30, std::vector<uint8_t>(kBlockSize, 0x33)},
  };

  auto container_or = KeyedObliviousContainer::Build(
      seed_l, seed_d, key_enc, real_blocks, /*num_dummies=*/5, kBlockSize);
  ASSERT_TRUE(container_or.ok()) << container_or.status();
  auto container = std::move(*container_or);

  EXPECT_EQ(container->NumRealBlocks(), 3);
  EXPECT_EQ(container->NumTotalBlocks(), 8);

  // Query block 20
  std::string query_label =
      KeyedObliviousContainer::ClientQuery(seed_l, 20, /*dummy_index=*/0);
  auto resp_or = container->ServerResponse(query_label);
  ASSERT_TRUE(resp_or.ok());

  auto recon_or = KeyedObliviousContainer::ClientRecon(
      key_enc, query_label, *resp_or, kBlockSize, /*is_dummy=*/false);
  ASSERT_TRUE(recon_or.ok());
  ASSERT_TRUE(recon_or->has_value());
  EXPECT_EQ((*recon_or)->block_id, 20);
  EXPECT_EQ((*recon_or)->payload, std::vector<uint8_t>(kBlockSize, 0x22));
}

TEST(ObliviousContainerTest, OkvsObliviousContainerPolynomialBackend) {
  std::string seed_k = GenerateRandomBytes(16);
  std::string key_enc = GenerateRandomBytes(16);

  constexpr size_t kBlockSize = 32;

  std::vector<std::pair<BlockId, std::vector<uint8_t>>> real_blocks = {
      {10, std::vector<uint8_t>(kBlockSize, 0x11)},
      {20, std::vector<uint8_t>(kBlockSize, 0x22)},
      {30, std::vector<uint8_t>(kBlockSize, 0x33)},
  };

  using PolyContainer = OkvsObliviousContainer<
      okvs::PolynomialOkvs<std::string, std::vector<uint8_t>>>;
  auto container_or = PolyContainer::Build(
      seed_k, key_enc, real_blocks, kBlockSize, /*seed=*/42,
      /*with_metadata=*/true, /*epsilon=*/0.15, /*fixed_size=*/3);
  ASSERT_TRUE(container_or.ok()) << container_or.status();
  auto container = std::move(*container_or);

  // Query block 20
  std::string query_label = PolyContainer::ClientQuery(seed_k, 20);
  auto resp_or = container->ServerResponse(query_label);
  ASSERT_TRUE(resp_or.ok());

  auto recon_or = PolyContainer::ClientRecon(key_enc, query_label, *resp_or,
                                             kBlockSize, /*is_dummy=*/false);
  ASSERT_TRUE(recon_or.ok());
  ASSERT_TRUE(recon_or->has_value());
  EXPECT_EQ((*recon_or)->block_id, 20);
  EXPECT_EQ((*recon_or)->payload, std::vector<uint8_t>(kBlockSize, 0x22));
}

TEST(ObliviousContainerTest, OkvsObliviousContainerRandomBandBackend) {
  std::string seed_k = GenerateRandomBytes(16);
  std::string key_enc = GenerateRandomBytes(16);

  constexpr size_t kBlockSize = 32;

  std::vector<std::pair<BlockId, std::vector<uint8_t>>> real_blocks;
  for (BlockId i = 0; i < 20; ++i) {
    real_blocks.emplace_back(
        i + 1, std::vector<uint8_t>(kBlockSize, static_cast<uint8_t>(i + 1)));
  }

  using BandContainer = OkvsObliviousContainer<
      okvs::RandomBandOkvs<std::string, std::vector<uint8_t>>>;
  auto container_or = BandContainer::Build(
      seed_k, key_enc, real_blocks, kBlockSize, /*seed=*/1001,
      /*with_metadata=*/true, /*epsilon=*/0.15);
  ASSERT_TRUE(container_or.ok()) << container_or.status();
  auto container = std::move(*container_or);

  // Query block 5
  std::string query_label = BandContainer::ClientQuery(seed_k, 5);
  auto resp_or = container->ServerResponse(query_label);
  ASSERT_TRUE(resp_or.ok());

  auto recon_or = BandContainer::ClientRecon(key_enc, query_label, *resp_or,
                                             kBlockSize, /*is_dummy=*/false);
  ASSERT_TRUE(recon_or.ok());
  ASSERT_TRUE(recon_or->has_value());
  EXPECT_EQ((*recon_or)->block_id, 5);
  EXPECT_EQ((*recon_or)->payload, std::vector<uint8_t>(kBlockSize, 5));

  // Dummy query
  std::string dummy_label = BandContainer::ClientQuery(seed_k, std::nullopt);
  auto dummy_resp_or = container->ServerResponse(dummy_label);
  ASSERT_TRUE(dummy_resp_or.ok());

  auto dummy_recon_or = BandContainer::ClientRecon(
      key_enc, dummy_label, *dummy_resp_or, kBlockSize, /*is_dummy=*/true);
  ASSERT_TRUE(dummy_recon_or.ok());
  EXPECT_FALSE(dummy_recon_or->has_value());
}

}  // namespace
}  // namespace oram_okvs
