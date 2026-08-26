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

#include "research/oram_okvs/oram/block_retrieval.h"

#include <cstdint>
#include <memory>
#include <vector>

#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "research/oram_okvs/oram/oram.pb.h"
#include "research/oram_okvs/pir/onion_pir_params.h"

namespace oram_okvs {
namespace {

TEST(BlockRetrievalTest, DirectRetrievalRoundTrip) {
  auto client = BlockRetrievalClient::CreateDirect();
  auto server = BlockRetrievalServer::CreateDirect();

  std::vector<std::vector<uint8_t>> candidates = {
      std::vector<uint8_t>(32, 0x11),
      std::vector<uint8_t>(32, 0x22),
      std::vector<uint8_t>(32, 0x33),
      std::vector<uint8_t>(32, 0x44),
  };

  size_t target_index = 2;
  auto query_or = client->CreateQuery(target_index, candidates.size(), 32);
  ASSERT_TRUE(query_or.ok());

  auto response_or = server->Evaluate(candidates, *query_or);
  ASSERT_TRUE(response_or.ok());

  auto decoded_or = client->DecodeResponse(*response_or, target_index, 32);
  ASSERT_TRUE(decoded_or.ok());
  EXPECT_EQ(*decoded_or, candidates[target_index]);
}

TEST(BlockRetrievalTest, XorRetrievalRoundTripWithDummyCancellation) {
  auto client = BlockRetrievalClient::CreateXor();
  auto server = BlockRetrievalServer::CreateXor();

  std::vector<std::vector<uint8_t>> dummy_pads = {
      std::vector<uint8_t>(32, 0xAA),
      std::vector<uint8_t>(32, 0xBB),
      std::vector<uint8_t>(32, 0x00),  // Real block slot (no dummy)
      std::vector<uint8_t>(32, 0xDD),
  };

  std::vector<uint8_t> real_block(32, 0x77);

  std::vector<std::vector<uint8_t>> server_blocks = {
      dummy_pads[0],
      dummy_pads[1],
      real_block,
      dummy_pads[3],
  };

  size_t target_index = 2;
  auto query_or = client->CreateQuery(target_index, server_blocks.size(), 32);
  ASSERT_TRUE(query_or.ok());

  auto response_or = server->Evaluate(server_blocks, *query_or);
  ASSERT_TRUE(response_or.ok());

  auto decoded_or =
      client->DecodeResponse(*response_or, target_index, 32, dummy_pads);
  ASSERT_TRUE(decoded_or.ok());
  EXPECT_EQ(*decoded_or, real_block);
}

TEST(BlockRetrievalTest, OnionPirRetrievalRoundTrip) {
  size_t num_candidates = 8;
  size_t entry_size = 32;
  auto params =
      pir::OnionPirParams::CreateFastForTesting(num_candidates, entry_size);
  params.n0 = 8;

  auto client_or = BlockRetrievalClient::CreateOnionPir(params);
  ASSERT_TRUE(client_or.ok());
  auto client = std::move(*client_or);

  auto server_or = BlockRetrievalServer::CreateOnionPir(params);
  ASSERT_TRUE(server_or.ok());
  auto server = std::move(*server_or);

  auto eval_keys_or = client->GetSerializedEvaluationKeys();
  ASSERT_TRUE(eval_keys_or.ok());
  ASSERT_TRUE(server->SetEvaluationKeys(*eval_keys_or).ok());

  std::vector<std::vector<uint8_t>> candidates(num_candidates);
  for (size_t i = 0; i < num_candidates; ++i) {
    candidates[i] =
        std::vector<uint8_t>(entry_size, static_cast<uint8_t>(i + 1));
  }

  for (size_t target_index = 0; target_index < num_candidates; ++target_index) {
    auto query_or =
        client->CreateQuery(target_index, num_candidates, entry_size);
    ASSERT_TRUE(query_or.ok());

    auto response_or = server->Evaluate(candidates, *query_or);
    ASSERT_TRUE(response_or.ok());

    auto decoded_or =
        client->DecodeResponse(*response_or, target_index, entry_size);
    ASSERT_TRUE(decoded_or.ok());
    EXPECT_EQ(*decoded_or, candidates[target_index]);
  }
}

}  // namespace
}  // namespace oram_okvs
