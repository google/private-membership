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

#include "research/oram_okvs/oram/lat/latoram.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "research/oram_okvs/oram/crypto_utils.h"
#include "research/oram_okvs/oram/lat/latoram_client.h"
#include "research/oram_okvs/oram/lat/latoram_server.h"
#include "research/oram_okvs/oram/lat/latoram_storage.h"
#include "research/oram_okvs/oram/lat/multi_perm.h"
#include "research/oram_okvs/oram/oram.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {
namespace {

// ============================================================================
// MultiPermDS Tests
// ============================================================================

TEST(MultiPermDSTest, BasicSampling) {
  std::vector<size_t> freqs = {3, 2, 1};
  MultiPermDS ds(freqs, 42);

  EXPECT_EQ(ds.TotalRemaining(), 6);
  EXPECT_EQ(ds.NumBins(), 3);

  std::vector<size_t> counts(3, 0);
  while (ds.TotalRemaining() > 0) {
    auto z_or = ds.Next();
    ASSERT_TRUE(z_or.ok());
    ASSERT_LT(*z_or, 3);
    counts[*z_or]++;
  }

  EXPECT_EQ(counts[0], 3);
  EXPECT_EQ(counts[1], 2);
  EXPECT_EQ(counts[2], 1);
  EXPECT_EQ(ds.TotalRemaining(), 0);
  EXPECT_FALSE(ds.Next().ok());
}

TEST(MultiPermDSTest, LargeBinSampling) {
  constexpr size_t kBins = 10;
  constexpr size_t kFreqPerBin = 20;
  std::vector<size_t> freqs(kBins, kFreqPerBin);
  MultiPermDS ds(freqs, 12345);

  EXPECT_EQ(ds.TotalRemaining(), kBins * kFreqPerBin);

  std::vector<size_t> counts(kBins, 0);
  while (ds.TotalRemaining() > 0) {
    auto z_or = ds.Next();
    ASSERT_TRUE(z_or.ok());
    counts[*z_or]++;
  }

  for (size_t i = 0; i < kBins; ++i) {
    EXPECT_EQ(counts[i], kFreqPerBin);
  }
}

// ============================================================================
// Crypto Utils Tests
// ============================================================================

TEST(CryptoUtilsTest, RandomBytesUniqueness) {
  std::string k1 = GenerateRandomBytes(16);
  std::string k2 = GenerateRandomBytes(16);
  EXPECT_EQ(k1.size(), 16);
  EXPECT_EQ(k2.size(), 16);
  EXPECT_NE(k1, k2);
}

TEST(CryptoUtilsTest, SubKeyDerivation) {
  std::string master = GenerateRandomBytes(16);
  std::string sub0 = DeriveSubKey(master, uint64_t{0});
  std::string sub1 = DeriveSubKey(master, uint64_t{1});
  std::string sub0_repeat = DeriveSubKey(master, uint64_t{0});

  EXPECT_EQ(sub0, sub0_repeat);
  EXPECT_NE(sub0, sub1);
}

TEST(CryptoUtilsTest, LabelGeneration) {
  std::string key = GenerateRandomBytes(16);
  std::string l_real0 = GenerateBlockLabel(key, 0);
  std::string l_real1 = GenerateBlockLabel(key, 1);
  std::string l_dummy0 = GenerateDummyLabel(key, 0);

  EXPECT_NE(l_real0, l_real1);
  EXPECT_NE(l_real0, l_dummy0);
}

TEST(CryptoUtilsTest, EncryptDecryptRoundtrip) {
  std::string key = GenerateRandomBytes(16);
  std::string label = GenerateBlockLabel(key, 42);

  std::vector<uint8_t> payload = {1, 2, 3, 4, 5, 6, 7, 8};
  constexpr size_t kBlockSize = 32;

  std::vector<uint8_t> cipher =
      EncryptBlock(key, label, 42, payload, kBlockSize);
  EXPECT_EQ(cipher.size(), CiphertextBlockSize(kBlockSize));

  auto dec_or = DecryptBlock(key, label, cipher, kBlockSize);
  ASSERT_TRUE(dec_or.ok());
  EXPECT_EQ(dec_or->block_id, 42);
  EXPECT_EQ(dec_or->payload.size(), kBlockSize);
  EXPECT_EQ(std::vector<uint8_t>(dec_or->payload.begin(),
                                 dec_or->payload.begin() + payload.size()),
            payload);
}

TEST(CryptoUtilsTest, DummyXorIsolation) {
  std::string key_e = GenerateRandomBytes(16);
  std::string key_l = GenerateRandomBytes(16);
  std::string key_d = GenerateRandomBytes(16);

  constexpr size_t kPayloadSize = 64;
  constexpr size_t kCipherSize = CiphertextBlockSize(kPayloadSize);

  BlockId target_id = 7;
  std::vector<uint8_t> target_payload(kPayloadSize, 0xAB);
  std::string real_label = GenerateBlockLabel(key_l, target_id);
  std::vector<uint8_t> real_cipher =
      EncryptBlock(key_e, real_label, target_id, target_payload, kPayloadSize);

  std::vector<std::string> dummy_labels = {
      GenerateDummyLabel(key_l, 1),
      GenerateDummyLabel(key_l, 2),
      GenerateDummyLabel(key_l, 3),
  };

  std::vector<std::vector<uint8_t>> dummy_blocks;
  for (const auto& d_lab : dummy_labels) {
    dummy_blocks.push_back(GenerateDummyBlock(key_d, d_lab, kCipherSize));
  }

  std::vector<uint8_t> z = real_cipher;
  for (const auto& d_blk : dummy_blocks) {
    XorInPlace(absl::MakeSpan(z), d_blk);
  }

  for (const auto& d_lab : dummy_labels) {
    std::vector<uint8_t> d_blk = GenerateDummyBlock(key_d, d_lab, kCipherSize);
    XorInPlace(absl::MakeSpan(z), d_blk);
  }

  auto dec_or = DecryptBlock(key_e, real_label, z, kPayloadSize);
  ASSERT_TRUE(dec_or.ok());
  EXPECT_EQ(dec_or->block_id, target_id);
  EXPECT_EQ(dec_or->payload, target_payload);
}

// ============================================================================
// LatORAM Server Proto Message Tests
// ============================================================================

TEST(LatOramServerTest, HandleInitAndQueryProtos) {
  LatOramServer server;

  // 1. Init Main via proto
  LatOramInitMainRequest init_req;
  auto* blk0 = init_req.add_blocks();
  blk0->set_label("L0");
  blk0->set_ciphertext(std::string({1, 2, 3}));

  auto* blk1 = init_req.add_blocks();
  blk1->set_label("L1");
  blk1->set_ciphertext(std::string({4, 5, 6}));

  auto init_resp_or = server.HandleInitMain(init_req);
  ASSERT_TRUE(init_resp_or.ok());
  EXPECT_TRUE(init_resp_or->success());
  EXPECT_EQ(server.MainSize(), 2);

  // 2. Add Lateral Stash via proto
  LatOramAddLateralStashRequest lat_req;
  lat_req.set_phase(1);
  auto* lat_blk = lat_req.add_blocks();
  lat_blk->set_label("L_lat0");
  lat_blk->set_ciphertext(std::string({10, 20, 30}));

  auto lat_resp_or = server.HandleAddLateralStash(lat_req);
  ASSERT_TRUE(lat_resp_or.ok());
  EXPECT_TRUE(lat_resp_or->success());
  EXPECT_EQ(server.NumLateralStashes(), 1);

  // 3. Query via proto
  LatOramQueryRequest query_req;
  query_req.add_labels("L0");
  query_req.add_labels("L_lat0");

  auto query_resp_or = server.HandleQuery(query_req);
  ASSERT_TRUE(query_resp_or.ok());
  EXPECT_TRUE(query_resp_or->success());

  std::vector<uint8_t> expected = {0x01 ^ 0x0A, 0x02 ^ 0x14, 0x03 ^ 0x1E};
  EXPECT_EQ(std::vector<uint8_t>(
                query_resp_or->retrieval_response().xor_block().begin(),
                query_resp_or->retrieval_response().xor_block().end()),
            expected);

  EXPECT_EQ(server.MainSize(), 1);
  EXPECT_EQ(server.total_queries_served(), 1);
}

// ============================================================================
// LatORAM SmallShuffle Tests
// ============================================================================

TEST(LatOramSmallShuffleTest, SmallShuffleIntegrity) {
  std::string key_e = GenerateRandomBytes(16);
  std::string key_l = GenerateRandomBytes(16);
  std::string key_d = GenerateRandomBytes(16);

  constexpr size_t kPayloadSize = 32;
  std::vector<std::pair<BlockId, std::vector<uint8_t>>> stash = {
      {10, std::vector<uint8_t>(kPayloadSize, 1)},
      {20, std::vector<uint8_t>(kPayloadSize, 2)},
      {30, std::vector<uint8_t>(kPayloadSize, 3)},
  };
  size_t num_dummies = 5;

  auto container = LatOramClient::SmallShuffle(stash, num_dummies, key_e, key_l,
                                               key_d, kPayloadSize, 42);

  EXPECT_EQ(container.size(), stash.size() + num_dummies);

  std::vector<BlockId> recovered_ids;
  for (const auto& [label, blk] : container) {
    auto dec_or = DecryptBlock(key_e, label, blk, kPayloadSize);
    if (dec_or.ok()) {
      recovered_ids.push_back(dec_or->block_id);
    }
  }

  std::sort(recovered_ids.begin(), recovered_ids.end());
  std::vector<BlockId> expected_ids = {10, 20, 30};
  EXPECT_EQ(recovered_ids, expected_ids);
}

// ============================================================================
// LatORAM End-to-End Client-Server Tests
// ============================================================================

TEST(LatOramClientServerTest, SeparateClientServerDeployment) {
  LatOramConfig config;
  config.num_blocks = 32;
  config.block_size = 64;
  config.stash_size = 4;
  config.num_phases = 2;
  config.seed = 123;

  auto server = std::make_shared<LatOramServer>();
  auto client = LatOramClient::Create(config, server);
  ASSERT_TRUE(client->Init().ok());

  std::vector<uint8_t> data(64, 0x55);
  data[0] = 12;
  data[63] = 34;

  ASSERT_TRUE(client->Write(7, data).ok());

  auto read_or = client->Read(7);
  ASSERT_TRUE(read_or.ok());
  EXPECT_EQ(*read_or, data);

  OramStats stats = client->GetStats();
  EXPECT_GT(stats.bytes_written, 0);
  EXPECT_GT(stats.bytes_read, 0);
  EXPECT_GT(stats.online_bytes_transferred, 0);
}

TEST(LatOramTest, UninitializedReadReturnsZeros) {
  LatOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.stash_size = 4;
  config.num_phases = 2;
  config.seed = 42;

  LatOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    auto read_or = oram.Read(i);
    ASSERT_TRUE(read_or.ok()) << "Failed to read block " << i;
    EXPECT_EQ(read_or->size(), config.block_size);
    for (uint8_t byte : *read_or) {
      EXPECT_EQ(byte, 0);
    }
  }
}

TEST(LatOramTest, WriteAndReadSingleBlock) {
  LatOramConfig config;
  config.num_blocks = 32;
  config.block_size = 64;
  config.stash_size = 4;
  config.num_phases = 2;
  config.seed = 99;

  LatOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  std::vector<uint8_t> data(64, 0xFE);
  data[0] = 42;
  data[63] = 99;

  ASSERT_TRUE(oram.Write(5, data).ok());

  auto read_or = oram.Read(5);
  ASSERT_TRUE(read_or.ok());
  EXPECT_EQ(*read_or, data);
}

TEST(LatOramTest, OverwriteBlockMultipleTimes) {
  LatOramConfig config;
  config.num_blocks = 16;
  config.block_size = 16;
  config.stash_size = 4;
  config.num_phases = 2;

  LatOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  for (int iter = 0; iter < 10; ++iter) {
    std::vector<uint8_t> data(16, static_cast<uint8_t>(iter * 11));
    ASSERT_TRUE(oram.Write(3, data).ok());
    auto read_or = oram.Read(3);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, data);
  }
}

TEST(LatOramTest, ReadAndWriteAtomicity) {
  LatOramConfig config;
  config.num_blocks = 16;
  config.block_size = 16;
  config.stash_size = 4;
  config.num_phases = 2;

  LatOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  std::vector<uint8_t> initial(16, 0x11);
  ASSERT_TRUE(oram.Write(2, initial).ok());

  std::vector<uint8_t> updated(16, 0x22);
  auto prev_or = oram.Access(2, AccessOp::kReadAndWrite, updated);
  ASSERT_TRUE(prev_or.ok());
  EXPECT_EQ(*prev_or, initial);

  auto curr_or = oram.Read(2);
  ASSERT_TRUE(curr_or.ok());
  EXPECT_EQ(*curr_or, updated);
}

TEST(LatOramTest, SequentialAllBlocksReadWrite) {
  LatOramConfig config;
  config.num_blocks = 64;
  config.block_size = 32;
  config.stash_size = 8;
  config.num_phases = 4;
  config.seed = 12345;

  LatOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    std::vector<uint8_t> data(config.block_size, static_cast<uint8_t>(i));
    ASSERT_TRUE(oram.Write(i, data).ok());
  }

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    std::vector<uint8_t> expected(config.block_size, static_cast<uint8_t>(i));
    auto read_or = oram.Read(i);
    ASSERT_TRUE(read_or.ok()) << "Failed to read block " << i;
    EXPECT_EQ(*read_or, expected) << "Mismatch at block " << i;
  }
}

TEST(LatOramTest, CrossPhaseAndMainReinitCycleWithServerShuffling) {
  LatOramConfig config;
  config.num_blocks = 32;
  config.block_size = 16;
  config.stash_size = 4;  // M = 4
  config.num_phases = 3;  // P = 3 (Reinit every P * M = 12 accesses)
  config.seed = 777;

  auto server = std::make_shared<LatOramServer>();
  auto client = LatOramClient::Create(config, server);
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
    EXPECT_EQ(*read_or, oracle[i])
        << "Mismatch after initial 32 writes for block " << i;
  }

  std::mt19937_64 rng(54321);
  for (int step = 0; step < 50; ++step) {
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

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    auto read_or = client->Read(i);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, oracle[i]);
  }

  OramStats stats = client->GetStats();
  EXPECT_GT(stats.total_accesses, 50);
  EXPECT_GT(stats.shuffles_count, 0);
  EXPECT_GT(stats.main_reinits_count, 0);
  EXPECT_GT(server->total_shuffles_served(), 0);
}

TEST(LatOramTest, OutOfBoundsAccessReturnsError) {
  LatOramConfig config;
  config.num_blocks = 10;
  config.block_size = 16;

  LatOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  EXPECT_FALSE(oram.Read(10).ok());
  EXPECT_FALSE(oram.Read(100).ok());
  EXPECT_FALSE(oram.Write(10, std::vector<uint8_t>(16, 1)).ok());
}

// ============================================================================
// LatORAM with OKVS Storage (Trivial PIR Evaluation) Tests
// ============================================================================

TEST(OkvsLatOramTest, BasicReadAndWrite) {
  LatOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.stash_size = 4;
  config.num_phases = 2;
  config.seed = 101;
  config.storage_type = LatOramStorageType::kOkLat;

  LatOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  std::vector<uint8_t> data(32, 0xAB);
  data[0] = 7;
  data[31] = 42;

  ASSERT_TRUE(oram.Write(5, data).ok());

  auto read_or = oram.Read(5);
  ASSERT_TRUE(read_or.ok());
  EXPECT_EQ(*read_or, data);
}

TEST(OkvsLatOramTest, OverwriteBlockMultipleTimes) {
  LatOramConfig config;
  config.num_blocks = 16;
  config.block_size = 16;
  config.stash_size = 4;
  config.num_phases = 2;
  config.seed = 202;
  config.storage_type = LatOramStorageType::kOkLat;

  LatOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  for (int iter = 0; iter < 10; ++iter) {
    std::vector<uint8_t> data(16, static_cast<uint8_t>(iter * 13));
    ASSERT_TRUE(oram.Write(2, data).ok());
    auto read_or = oram.Read(2);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, data);
  }
}

TEST(OkvsLatOramTest, CrossPhaseAndMainReinitCycleWithOkvsStorage) {
  LatOramConfig config;
  config.num_blocks = 32;
  config.block_size = 16;
  config.stash_size = 4;  // M = 4
  config.num_phases = 3;  // P = 3 (Reinit every P * M = 12 accesses)
  config.seed = 777;
  config.storage_type = LatOramStorageType::kOkLat;

  auto storage = std::make_unique<LatOramStorage>();
  auto server = std::make_shared<LatOramServer>(std::move(storage));

  auto client = LatOramClient::Create(config, server);
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
    EXPECT_EQ(*read_or, oracle[i])
        << "Mismatch after initial 32 writes for block " << i;
  }

  std::mt19937_64 rng(54321);
  for (int step = 0; step < 40; ++step) {
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

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    auto read_or = client->Read(i);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, oracle[i]);
  }

  OramStats stats = client->GetStats();
  EXPECT_GT(stats.total_accesses, 40);
  EXPECT_GT(stats.shuffles_count, 0);
  EXPECT_GT(stats.main_reinits_count, 0);
}

TEST(OkvsLatOramTest, PolynomialOkvsBackend) {
  LatOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.stash_size = 4;
  config.num_phases = 2;
  config.storage_type = LatOramStorageType::kOkLat;
  config.okvs_backend = OkvsBackendType::kPolynomial;
  config.seed = 42;

  LatOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    std::vector<uint8_t> data(32, static_cast<uint8_t>(i + 15));
    ASSERT_TRUE(oram.Write(i, data).ok());
  }

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    auto read_or = oram.Read(i);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, std::vector<uint8_t>(32, static_cast<uint8_t>(i + 15)));
  }
}

// ============================================================================
// OkvsPirLatOramTest (LatORAM with OKVS + OnionPIR protocol)
// ============================================================================

TEST(OkvsPirLatOramTest, BasicReadAndWriteWithOnionPir) {
  LatOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.stash_size = 4;
  config.num_phases = 2;
  config.seed = 101;
  config.storage_type = LatOramStorageType::kOkPLat;

  auto storage = std::make_unique<LatOramStorage>();
  auto server = std::make_shared<LatOramServer>(std::move(storage));
  auto client = LatOramClient::Create(config, server);
  auto init_status = client->Init();
  ASSERT_TRUE(init_status.ok()) << init_status;

  std::vector<uint8_t> data(32, 0xAB);
  data[0] = 7;
  data[31] = 42;

  ASSERT_TRUE(client->Write(5, data).ok());

  auto read_or = client->Read(5);
  ASSERT_TRUE(read_or.ok());
  EXPECT_EQ(*read_or, data);
}

TEST(OkvsPirLatOramTest, OverwriteBlockMultipleTimesWithOnionPir) {
  LatOramConfig config;
  config.num_blocks = 16;
  config.block_size = 16;
  config.stash_size = 4;
  config.num_phases = 2;
  config.seed = 202;
  config.storage_type = LatOramStorageType::kOkPLat;

  auto storage = std::make_unique<LatOramStorage>();
  auto server = std::make_shared<LatOramServer>(std::move(storage));
  auto client = LatOramClient::Create(config, server);
  auto init_status = client->Init();
  ASSERT_TRUE(init_status.ok()) << init_status;

  for (int iter = 0; iter < 10; ++iter) {
    std::vector<uint8_t> data(16, static_cast<uint8_t>(iter * 13 + 1));
    ASSERT_TRUE(client->Write(2, data).ok());
    auto read_or = client->Read(2);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, data);
  }
}

TEST(OkvsPirLatOramTest, CrossPhaseAndMainReinitCycleWithOnionPir) {
  LatOramConfig config;
  config.num_blocks = 16;
  config.block_size = 16;
  config.stash_size = 4;  // M = 4
  config.num_phases = 2;  // P = 2 (Reinit every P * M = 8 accesses)
  config.seed = 777;
  config.storage_type = LatOramStorageType::kOkPLat;

  auto storage = std::make_unique<LatOramStorage>();
  auto server = std::make_shared<LatOramServer>(std::move(storage));

  auto client = LatOramClient::Create(config, server);
  auto init_status = client->Init();
  ASSERT_TRUE(init_status.ok()) << init_status;

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
    EXPECT_EQ(*read_or, oracle[i])
        << "Mismatch after initial writes for block " << i;
  }

  std::mt19937_64 rng(54321);
  for (int step = 0; step < 20; ++step) {
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

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    auto read_or = client->Read(i);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, oracle[i]);
  }

  OramStats stats = client->GetStats();
  EXPECT_GT(stats.total_accesses, 20);
  EXPECT_GT(stats.shuffles_count, 0);
  EXPECT_GT(stats.main_reinits_count, 0);
}

TEST(OkvsPirLatOramTest, PolynomialOkvsBackendWithOnionPir) {
  LatOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.stash_size = 4;
  config.num_phases = 2;
  config.storage_type = LatOramStorageType::kOkPLat;
  config.okvs_backend = OkvsBackendType::kPolynomial;
  config.seed = 999;

  LatOram oram(config);
  ASSERT_TRUE(oram.Init().ok());

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    std::vector<uint8_t> data(32, static_cast<uint8_t>(i + 50));
    ASSERT_TRUE(oram.Write(i, data).ok());
  }

  for (BlockId i = 0; i < config.num_blocks; ++i) {
    auto read_or = oram.Read(i);
    ASSERT_TRUE(read_or.ok());
    EXPECT_EQ(*read_or, std::vector<uint8_t>(32, static_cast<uint8_t>(i + 50)));
  }
}

TEST(LatOramTest, LargeUint64BlockIdSupport) {
  LatOramConfig config;
  config.num_blocks = 1000000;  // 1 million blocks
  config.block_size = 64;
  config.stash_size = 8;
  config.num_phases = 2;
  config.seed = 9999;
  config.storage_type = LatOramStorageType::kMap;

  auto storage = std::make_unique<LatOramStorage>();
  auto server = std::make_shared<LatOramServer>(std::move(storage));
  auto client = LatOramClient::Create(config, server);
  ASSERT_TRUE(client->Init().ok());

  // Test accessing large 64-bit Block IDs
  BlockId large_id1 = 0x0000000000000001ULL;
  BlockId large_id2 = 999999ULL;

  std::vector<uint8_t> val1(64, 0x11);
  std::vector<uint8_t> val2(64, 0x22);

  ASSERT_TRUE(client->Write(large_id1, val1).ok());
  ASSERT_TRUE(client->Write(large_id2, val2).ok());

  auto r1 = client->Read(large_id1);
  ASSERT_TRUE(r1.ok());
  EXPECT_EQ(*r1, val1);

  auto r2 = client->Read(large_id2);
  ASSERT_TRUE(r2.ok());
  EXPECT_EQ(*r2, val2);
}

TEST(OkvsLatOramTest, NoDummiesInOkvsContainers) {
  LatOramConfig config;
  config.num_blocks = 16;
  config.block_size = 32;
  config.stash_size = 4;
  config.num_phases = 2;
  config.storage_type = LatOramStorageType::kOkLat;
  config.seed = 12345;

  auto storage = std::make_unique<LatOramStorage>();
  auto* raw_storage = storage.get();
  auto server = std::make_shared<LatOramServer>(std::move(storage));

  auto client = LatOramClient::Create(config, server);
  ASSERT_TRUE(client->Init().ok());

  // Main container retains P * M dummies (16 + 2*4 = 24) to support streaming
  // in O(M) client RAM.
  EXPECT_EQ(raw_storage->MainSize(), 24);

  // Fill stash to trigger lateral stash creation
  for (BlockId i = 0; i < 4; ++i) {
    std::vector<uint8_t> data(32, static_cast<uint8_t>(i + 1));
    ASSERT_TRUE(client->Write(i, data).ok());
  }

  // Lateral stash is dummy-free (exactly 4 blocks, not 4 + dummies)
  EXPECT_EQ(raw_storage->NumLateralStashes(), 1);
  EXPECT_EQ(raw_storage->LateralStashSize(0), 4);
}

}  // namespace
}  // namespace oram_okvs
