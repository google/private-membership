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

#include "research/oram_okvs/oram/ring/ringoram.h"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "research/oram_okvs/oram/oram.h"
#include "research/oram_okvs/oram/ring/ok_ringoram_client.h"
#include "research/oram_okvs/oram/ring/okp_ringoram_client.h"
#include "research/oram_okvs/oram/ring/ringoram_client.h"
#include "research/oram_okvs/oram/ring/ringoram_server.h"
#include "research/oram_okvs/oram/ring/ringoram_storage.h"
#include "research/oram_okvs/oram/ring/tree_ringoram_client.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {
namespace {

TEST(RingOramTest, ConfigCalculations) {
  RingOramConfig config;
  config.num_blocks = 16;
  config.block_size = 64;
  config.z_real = 2;
  config.s_dummy = 2;

  EXPECT_EQ(config.TreeDepth(), 4);
  EXPECT_EQ(config.NumLeaves(), 16);
  EXPECT_EQ(config.TotalBuckets(), 31);
  EXPECT_EQ(config.BucketSize(), 4);
}

TEST(RingOramTest, ReverseBits) {
  EXPECT_EQ(ReverseBits(0, 4), 0);
  EXPECT_EQ(ReverseBits(1, 4), 8);    // 0001 -> 1000
  EXPECT_EQ(ReverseBits(2, 4), 4);    // 0010 -> 0100
  EXPECT_EQ(ReverseBits(3, 4), 12);   // 0011 -> 1100
  EXPECT_EQ(ReverseBits(7, 4), 14);   // 0111 -> 1110
  EXPECT_EQ(ReverseBits(15, 4), 15);  // 1111 -> 1111
}

TEST(RingOramTest, GetPathBucketIds) {
  // Depth = 3: root=0, level 1: 1,2; level 2: 3,4,5,6; level 3 (leaves): 7..14
  // Leaf 0 -> leaf bucket 7. Path from root: 0, 1, 3, 7
  std::vector<uint32_t> path0 = GetPathBucketIds(0, 3);
  std::vector<uint32_t> expected0 = {0, 1, 3, 7};
  EXPECT_EQ(path0, expected0);

  // Leaf 7 -> leaf bucket 14. Path from root: 0, 2, 6, 14
  std::vector<uint32_t> path7 = GetPathBucketIds(7, 3);
  std::vector<uint32_t> expected7 = {0, 2, 6, 14};
  EXPECT_EQ(path7, expected7);
}

TEST(RingOramTest, BasicReadAndWrite) {
  RingOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.z_real = 2;
  config.s_dummy = 2;
  config.eviction_rate = 1;
  config.seed = 42;

  RingOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  std::vector<uint8_t> data1(32, 0xAA);
  data1[0] = 5;
  data1[31] = 99;

  ASSERT_TRUE(oram.Write(3, data1).ok());

  auto read_or = oram.Read(3);
  ASSERT_TRUE(read_or.ok());
  EXPECT_EQ(*read_or, data1);

  // Overwrite block 3
  std::vector<uint8_t> data2(32, 0xBB);
  data2[0] = 42;
  ASSERT_TRUE(oram.Write(3, data2).ok());

  auto read2_or = oram.Read(3);
  ASSERT_TRUE(read2_or.ok());
  EXPECT_EQ(*read2_or, data2);

  // Read unwritten block (should return zeroes)
  auto read_unwritten_or = oram.Read(5);
  ASSERT_TRUE(read_unwritten_or.ok());
  EXPECT_EQ(*read_unwritten_or, std::vector<uint8_t>(32, 0));
}

TEST(RingOramTest, SequentialAndRandomAccessWithEvictions) {
  RingOramConfig config;
  config.num_blocks = 32;
  config.block_size = 16;
  config.z_real = 2;
  config.s_dummy = 2;
  config.eviction_rate = 1;
  config.seed = 789;

  auto server = std::make_shared<RingOramServer>();
  auto client = RingOramClient::Create(config, server);
  ASSERT_TRUE(client->Init().ok());

  std::vector<std::vector<uint8_t>> oracle(
      config.num_blocks, std::vector<uint8_t>(config.block_size, 0));

  // Initialize all blocks with distinct values
  for (BlockId i = 0; i < config.num_blocks; ++i) {
    oracle[i] =
        std::vector<uint8_t>(config.block_size, static_cast<uint8_t>(i + 1));
    ASSERT_TRUE(client->Write(i, oracle[i]).ok());
  }

  // Perform mixed random reads and writes
  std::mt19937_64 rng(99999);
  for (int step = 0; step < 60; ++step) {
    BlockId block_id = rng() % config.num_blocks;
    bool is_write = (rng() % 2 == 0);

    if (is_write) {
      std::vector<uint8_t> new_val(config.block_size,
                                   static_cast<uint8_t>(rng() & 0xFF));
      oracle[block_id] = new_val;
      ASSERT_TRUE(client->Write(block_id, new_val).ok());
    } else {
      auto read_or = client->Read(block_id);
      ASSERT_TRUE(read_or.ok())
          << "Read failed at step " << step << " for block " << block_id;
      EXPECT_EQ(*read_or, oracle[block_id])
          << "Mismatch at step " << step << " for block " << block_id;
    }
  }

  // Final verification of all blocks
  for (BlockId i = 0; i < config.num_blocks; ++i) {
    auto read_or = client->Read(i);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, oracle[i]) << "Final mismatch at block " << i;
  }

  OramStats stats = client->GetStats();
  EXPECT_GT(stats.total_accesses, 60);
  EXPECT_GT(stats.shuffles_count, 0);  // Evictions occurred
  EXPECT_GT(server->total_accesses_served(), 0);
  EXPECT_GT(server->total_evictions_served(), 0);
}

TEST(RingOramTest, OutOfBoundsAccessReturnsError) {
  RingOramConfig config;
  config.num_blocks = 10;
  config.block_size = 16;

  RingOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  EXPECT_FALSE(oram.Read(10).ok());
  EXPECT_FALSE(oram.Read(100).ok());
  EXPECT_FALSE(oram.Write(10, std::vector<uint8_t>(16, 1)).ok());
}

TEST(RingOramTest, SmallReproducer) {
  RingOramConfig config;
  config.num_blocks = 8;
  config.block_size = 16;
  config.z_real = 2;
  config.s_dummy = 2;
  config.eviction_rate = 1;
  config.seed = 123;

  auto server = std::make_shared<RingOramServer>();
  auto client = RingOramClient::Create(config, server);
  ASSERT_TRUE(client->Init().ok());

  std::vector<std::vector<uint8_t>> oracle(
      config.num_blocks, std::vector<uint8_t>(config.block_size, 0));

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    oracle[i] =
        std::vector<uint8_t>(config.block_size, static_cast<uint8_t>(i + 1));
    ASSERT_TRUE(client->Write(i, oracle[i]).ok());
  }

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    auto read_or = client->Read(i);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, oracle[i]) << "Initial read mismatch at " << i;
  }
}

TEST(OkvsRingOramTest, BasicReadAndWriteWithOkvs) {
  RingOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.storage_type = RingOramStorageType::kOkRing;
  config.seed = 42;

  auto storage = std::make_unique<RingOramStorage>();
  auto server = std::make_shared<RingOramServer>(std::move(storage));

  auto client = RingOramClient::Create(config, server);
  ASSERT_TRUE(client->Init().ok());

  std::vector<uint8_t> data1(32, 0xAA);
  data1[0] = 5;
  data1[31] = 99;

  ASSERT_TRUE(client->Write(3, data1).ok());

  auto read_or = client->Read(3);
  ASSERT_TRUE(read_or.ok());
  EXPECT_EQ(*read_or, data1);
}

TEST(OkvsRingOramTest, SequentialAndRandomAccessWithEvictions) {
  RingOramConfig config;
  config.num_blocks = 32;
  config.block_size = 32;
  config.z_real = 2;
  config.s_dummy = 2;
  config.eviction_rate = 1;
  config.storage_type = RingOramStorageType::kOkRing;
  config.seed = 9999;

  RingOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  std::vector<std::vector<uint8_t>> oracle(
      config.num_blocks, std::vector<uint8_t>(config.block_size, 0));

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    oracle[i] =
        std::vector<uint8_t>(config.block_size, static_cast<uint8_t>(i + 1));
    ASSERT_TRUE(oram.Write(i, oracle[i]).ok());
  }

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    auto read_or = oram.Read(i);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, oracle[i]) << "Initial read mismatch at " << i;
  }

  std::mt19937_64 rng(88888);
  for (int step = 0; step < 40; ++step) {
    BlockId block_id = rng() % config.num_blocks;
    bool is_write = (rng() % 2 == 0);

    if (is_write) {
      std::vector<uint8_t> new_val(config.block_size,
                                   static_cast<uint8_t>(rng() & 0xFF));
      oracle[block_id] = new_val;
      ASSERT_TRUE(oram.Write(block_id, new_val).ok());
    } else {
      auto read_or = oram.Read(block_id);
      ASSERT_TRUE(read_or.ok())
          << "Read failed at step " << step << " for block " << block_id;
      EXPECT_EQ(*read_or, oracle[block_id])
          << "Mismatch at step " << step << " for block " << block_id;
    }
  }

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    auto read_or = oram.Read(i);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, oracle[i]) << "Final mismatch at block " << i;
  }

  OramStats stats = oram.GetStats();
  EXPECT_GT(stats.total_accesses, 40);
  EXPECT_GT(stats.shuffles_count, 0);
}

TEST(OkvsRingOramTest, RandomBandOkvs) {
  RingOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.z_real = 4;
  config.s_dummy = 2;
  config.storage_type = RingOramStorageType::kOkRing;
  config.okvs_backend = OkvsBackendType::kRandomBand;
  config.seed = 54321;

  RingOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    std::vector<uint8_t> data(32, static_cast<uint8_t>(i + 10));
    ASSERT_TRUE(oram.Write(i, data).ok());
  }

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    auto read_or = oram.Read(i);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, std::vector<uint8_t>(32, static_cast<uint8_t>(i + 10)));
  }
}

TEST(OkvsPirRingOramTest, BasicReadAndWriteWithOnionPir) {
  RingOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.storage_type = RingOramStorageType::kOkPRing;
  config.seed = 42;

  auto storage = std::make_unique<RingOramStorage>();
  auto server = std::make_shared<RingOramServer>(std::move(storage));
  auto client = RingOramClient::Create(config, server);
  ASSERT_TRUE(client->Init().ok());

  std::vector<uint8_t> data1(32, 0xAA);
  data1[0] = 5;
  data1[31] = 99;

  ASSERT_TRUE(client->Write(3, data1).ok());

  auto read_or = client->Read(3);
  ASSERT_TRUE(read_or.ok());
  EXPECT_EQ(*read_or, data1);
}

TEST(OkvsPirRingOramTest, RandomBandOkvsWithOnionPir) {
  RingOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.z_real = 4;
  config.s_dummy = 2;
  config.storage_type = RingOramStorageType::kOkPRing;
  config.okvs_backend = OkvsBackendType::kRandomBand;
  config.seed = 1234;

  RingOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    std::vector<uint8_t> data(32, static_cast<uint8_t>(i + 20));
    ASSERT_TRUE(oram.Write(i, data).ok());
  }

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    auto read_or = oram.Read(i);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, std::vector<uint8_t>(32, static_cast<uint8_t>(i + 20)));
  }
}

TEST(OkvsPirRingOramTest, SequentialAndRandomAccessWithEvictions) {
  RingOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.z_real = 2;
  config.s_dummy = 2;
  config.eviction_rate = 1;
  config.storage_type = RingOramStorageType::kOkPRing;
  config.seed = 9999;

  RingOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  std::vector<std::vector<uint8_t>> oracle(
      config.num_blocks, std::vector<uint8_t>(config.block_size, 0));

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    oracle[i] =
        std::vector<uint8_t>(config.block_size, static_cast<uint8_t>(i + 1));
    ASSERT_TRUE(oram.Write(i, oracle[i]).ok());
  }

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    auto read_or = oram.Read(i);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, oracle[i]) << "Initial read mismatch at " << i;
  }

  std::mt19937_64 rng(77777);
  for (int step = 0; step < 20; ++step) {
    BlockId block_id = rng() % config.num_blocks;
    bool is_write = (rng() % 2 == 0);

    if (is_write) {
      std::vector<uint8_t> new_val(config.block_size,
                                   static_cast<uint8_t>(rng() & 0xFF));
      oracle[block_id] = new_val;
      ASSERT_TRUE(oram.Write(block_id, new_val).ok());
    } else {
      auto read_or = oram.Read(block_id);
      ASSERT_TRUE(read_or.ok());
      EXPECT_EQ(*read_or, oracle[block_id])
          << "Mismatch at step " << step << " for block " << block_id;
    }
  }

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    auto read_or = oram.Read(i);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, oracle[i]) << "Final mismatch at block " << i;
  }

  OramStats stats = oram.GetStats();
  EXPECT_GT(stats.total_accesses, 20);
  EXPECT_GT(stats.shuffles_count, 0);
}

TEST(RingOramTest, LargeUint64BlockIdSupport) {
  RingOramConfig config;
  config.num_blocks = 1000000;  // 1 million blocks (Depth = 20)
  config.block_size = 64;
  config.z_real = 2;
  config.s_dummy = 2;
  config.seed = 8888;
  config.storage_type = RingOramStorageType::kTree;

  auto server = std::make_shared<RingOramServer>();
  auto client = RingOramClient::Create(config, server);
  ASSERT_TRUE(client->Init().ok());

  BlockId large_id1 = 0x0000000000000005ULL;
  BlockId large_id2 = 999999ULL;

  std::vector<uint8_t> val1(64, 0xAA);
  std::vector<uint8_t> val2(64, 0xBB);

  ASSERT_TRUE(client->Write(large_id1, val1).ok());
  ASSERT_TRUE(client->Write(large_id2, val2).ok());

  auto r1 = client->Read(large_id1);
  ASSERT_TRUE(r1.ok());
  EXPECT_EQ(*r1, val1);

  auto r2 = client->Read(large_id2);
  ASSERT_TRUE(r2.ok());
  EXPECT_EQ(*r2, val2);
}

TEST(OkvsRingOramTest, NoDummiesInOkvsBuckets) {
  RingOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.z_real = 2;
  config.s_dummy = 2;
  config.storage_type = RingOramStorageType::kOkRing;
  config.seed = 12345;

  auto storage = std::make_unique<RingOramStorage>();
  auto server = std::make_shared<RingOramServer>(std::move(storage));

  auto client = RingOramClient::Create(config, server);
  ASSERT_TRUE(client->Init().ok());

  // Initial empty bucket built via BuildBucket has 0 blocks and non-empty OKVS
  // table
  auto empty_bucket = client->BuildBucket(0, {});
  EXPECT_EQ(empty_bucket.blocks_size(), 0);
  EXPECT_FALSE(empty_bucket.okvs_table().empty());

  // Populated bucket has 0 raw blocks (no dummies, all encoded in okvs_table)
  std::vector<uint8_t> payload(32, 0xAB);
  auto populated_bucket = client->BuildBucket(0, {{1, payload}});
  EXPECT_EQ(populated_bucket.blocks_size(), 0);
  EXPECT_FALSE(populated_bucket.okvs_table().empty());
}

TEST(RingOramTest, ClientStoredPositionMapTracking) {
  RingOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.z_real = 2;
  config.s_dummy = 2;
  config.eviction_rate = 1;
  config.seed = 42;

  auto server = std::make_shared<RingOramServer>();
  auto client = RingOramClient::Create(config, server);
  ASSERT_TRUE(client->Init().ok());

  // Initially, all blocks are marked in stash or unassigned
  EXPECT_EQ(client->GetBucket(0), -1);
  EXPECT_EQ(client->GetSlot(0), -1);

  // Write block 0
  std::vector<uint8_t> data(32, 0xAA);
  ASSERT_TRUE(client->Write(0, data).ok());

  // Eviction occurred during Write(0), placing block 0 into a bucket in the
  // tree or remaining in stash
  int32_t b0_bucket = client->GetBucket(0);
  int32_t b0_slot = client->GetSlot(0);
  if (b0_bucket >= 0) {
    EXPECT_GE(b0_slot, 0);
    EXPECT_LT(b0_slot, static_cast<int32_t>(config.BucketSize()));
  }

  // Read block 0
  auto read_or = client->Read(0);
  ASSERT_TRUE(read_or.ok());
  EXPECT_EQ(*read_or, data);
}

TEST(RingOramTest, OramInterfacePolymorphicUsage) {
  RingOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.z_real = 2;
  config.s_dummy = 2;
  config.seed = 99;

  std::unique_ptr<RingOram> oram = std::make_unique<RingOram>(config);
  ASSERT_TRUE(oram->Init().ok());

  Oram* oram_ptr = oram.get();
  EXPECT_EQ(oram_ptr->num_blocks(), 16);
  EXPECT_EQ(oram_ptr->block_size(), 32);

  std::vector<uint8_t> data(32, 0x42);
  ASSERT_TRUE(oram_ptr->Write(3, data).ok());

  auto read_or = oram_ptr->Read(3);
  ASSERT_TRUE(read_or.ok());
  EXPECT_EQ(*read_or, data);

  // Test Access with kReadAndWrite
  std::vector<uint8_t> new_data(32, 0x99);
  auto access_or = oram_ptr->Access(3, AccessOp::kReadAndWrite, new_data);
  ASSERT_TRUE(access_or.ok());
  EXPECT_EQ(*access_or, data);  // Returns previous value

  auto read_new_or = oram_ptr->Read(3);
  ASSERT_TRUE(read_new_or.ok());
  EXPECT_EQ(*read_new_or, new_data);

  OramStats stats = oram_ptr->GetStats();
  EXPECT_GT(stats.total_accesses, 0);
  oram_ptr->ResetStats();
  EXPECT_EQ(oram_ptr->GetStats().total_accesses, 0);
}

TEST(RingOramTest, DirectSlotAccessWithoutMetadataTransfer) {
  RingOramConfig config;
  config.num_blocks = 16;
  config.block_size = 64;
  config.z_real = 2;
  config.s_dummy = 2;
  config.eviction_rate = 0;  // Disable evictions to isolate online access cost
  config.seed = 1234;

  auto server = std::make_shared<RingOramServer>();
  auto client = RingOramClient::Create(config, server);
  ASSERT_TRUE(client->Init().ok());

  client->ResetStats();
  std::vector<uint8_t> data(64, 0x55);
  ASSERT_TRUE(client->Write(2, data).ok());

  // In 1-RTT access: online bytes transferred contains only the single
  // ReadPathBlocks request and response
  OramStats stats = client->GetStats();
  EXPECT_EQ(stats.total_accesses, 1);
  EXPECT_GT(stats.online_bytes_transferred, 0);
  EXPECT_EQ(stats.shuffles_count, 0);
}

}  // namespace
}  // namespace oram_okvs
