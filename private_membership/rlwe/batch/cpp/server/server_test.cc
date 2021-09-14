// Copyright 2021 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "private_membership/rlwe/batch/cpp/server/server.h"

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <google/protobuf/repeated_field.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/memory/memory.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "private_membership/rlwe/batch/cpp/client/client.h"
#include "private_membership/rlwe/batch/cpp/constants.h"
#include "private_membership/rlwe/batch/cpp/context.h"
#include "private_membership/rlwe/batch/cpp/encoding.h"
#include "private_membership/rlwe/batch/cpp/test_helper.h"
#include "private_membership/rlwe/batch/proto/client.pb.h"
#include "private_membership/rlwe/batch/proto/server.pb.h"
#include "private_membership/rlwe/batch/proto/shared.pb.h"
#include "shell_encryption/context.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/serialization.pb.h"
#include "shell_encryption/symmetric_encryption.h"

namespace private_membership {
namespace batch {
namespace {

using ::testing::UnorderedElementsAreArray;

const int kEmptyBucketId = kNumberOfBucketsPerShard - 2;

// Define struct to enable initialization.
struct RawShard {
  int shard_id;

  // List of bucket IDs and contents.
  std::vector<std::pair<int, std::string>> buckets;
};

bool operator==(const RawShard& lhs, const RawShard& rhs) {
  return (lhs.shard_id == rhs.shard_id) && (lhs.buckets == rhs.buckets);
}

std::string CreateTestBucketContent(int shard_id, int bucket_id) {
  // Second last bucket will be empty.
  if (bucket_id == kEmptyBucketId) {
    return "";
  }
  return absl::StrCat(shard_id, ",", bucket_id);
}

std::vector<RawShard> CreateTestDatabase() {
  std::vector<RawShard> database(kNumberOfShards);
  for (int i = 0; i < kNumberOfShards; ++i) {
    database[i].shard_id = i;
    for (int j = 0; j < kNumberOfBucketsPerShard; ++j) {
      database[i].buckets.push_back({j, CreateTestBucketContent(i, j)});
    }
  }
  return database;
}

EncodeDatabaseRequest CreateEncodeDatabaseRequest(
    const std::vector<RawShard>& raw_shards) {
  EncodeDatabaseRequest request;
  *request.mutable_parameters() = CreateTestParameters();
  for (const RawShard& raw_shard : raw_shards) {
    RawDatabaseShard* shard = request.add_shards();
    shard->set_shard_index(raw_shard.shard_id);
    for (const auto& [bucket_id, bucket_contents] : raw_shard.buckets) {
      RawDatabaseShard::Bucket* request_bucket = shard->add_buckets();
      request_bucket->set_bucket_id(bucket_id);
      request_bucket->set_bucket_contents(bucket_contents);
    }
  }
  return request;
}

absl::StatusOr<SumCiphertextsRequest> CreateSumCiphertextsRequest(
    const std::vector<std::vector<Int>>& addends, const Context& context,
    const rlwe::SymmetricRlweKey<ModularInt>& private_key) {
  SumCiphertextsRequest request;
  *request.mutable_parameters() = CreateTestParameters();
  for (const std::vector<Int>& sublist : addends) {
    SumCiphertextsRequest::Summation* summation = request.add_summations();
    for (Int addend : sublist) {
      absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> ciphertext =
          EncryptSingleInt(context, private_key, addend);
      if (!ciphertext.ok()) {
        return ciphertext.status();
      }
      auto serialized = ciphertext->Serialize();
      if (!serialized.ok()) {
        return serialized.status();
      }
      *summation->add_ciphertexts() = *std::move(serialized);
    }
  }
  return request;
}

absl::StatusOr<FinalizeResultsRequest> CreateFinalizeResultsRequest(
    const std::vector<std::vector<Int>>& results, const Context& context,
    const rlwe::SymmetricRlweKey<ModularInt>& private_key) {
  FinalizeResultsRequest request;
  *request.mutable_parameters() = CreateTestParameters();
  for (const std::vector<Int>& result_sublist : results) {
    EncryptedQueryResult* query_result = request.add_encrypted_results();
    for (Int partial_result : result_sublist) {
      absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> ciphertext =
          EncryptSingleInt(context, private_key, partial_result);
      if (!ciphertext.ok()) {
        return ciphertext.status();
      }
      auto serialized = ciphertext->Serialize();
      if (!serialized.ok()) {
        return serialized.status();
      }
      *query_result->add_ciphertexts() = *std::move(serialized);
    }
  }
  return request;
}

class ServerTest : public ::testing::Test {
 protected:
  void SetUp() override {
    auto request_context = CreateRlweRequestContext(CreateTestParameters());
    ASSERT_TRUE(request_context.ok());
    ASSERT_NE(request_context->get(), nullptr);
    request_context_ = *std::move(request_context);

    auto response_context = CreateRlweResponseContext(CreateTestParameters());
    ASSERT_TRUE(response_context.ok());
    ASSERT_NE(response_context->get(), nullptr);
    response_context_ = *std::move(response_context);

    absl::StatusOr<GenerateKeysResponse> keys = CreateTestKeys();
    ASSERT_TRUE(keys.ok());
    keys_ = *std::move(keys);

    absl::StatusOr<rlwe::SymmetricRlweKey<ModularInt>> request_private_key =
        DeserializePrivateKey(keys_.private_key().request_key(),
                              *request_context_);
    ASSERT_TRUE(request_private_key.ok());
    request_private_key_ =
        absl::make_unique<rlwe::SymmetricRlweKey<ModularInt>>(
            *std::move(request_private_key));

    absl::StatusOr<rlwe::SymmetricRlweKey<ModularInt>> response_private_key =
        DeserializePrivateKey(keys_.private_key().response_key(),
                              *response_context_);
    ASSERT_TRUE(response_private_key.ok());
    response_private_key_ =
        absl::make_unique<rlwe::SymmetricRlweKey<ModularInt>>(
            *std::move(response_private_key));
  }

  absl::StatusOr<ApplyQueriesRequest> CreateApplyQueriesRequest(
      const std::vector<PlaintextQuery>& plaintext_queries,
      bool encode_database) {
    return CreateApplyQueriesRequest(plaintext_queries, CreateTestDatabase(),
                                     encode_database);
  }

  absl::StatusOr<ApplyQueriesRequest> CreateApplyQueriesRequest(
      const std::vector<PlaintextQuery>& plaintext_queries,
      const std::vector<RawShard>& database, bool encode_database) {
    absl::StatusOr<EncryptedQueries> encrypted_queries =
        CreateEncryptedQueries(keys_, plaintext_queries);
    if (!encrypted_queries.ok()) {
      return encrypted_queries.status();
    }

    ApplyQueriesRequest request;
    *request.mutable_parameters() = CreateTestParameters();
    *request.mutable_public_key() = keys_.public_key();
    *request.add_queries() = *std::move(encrypted_queries);

    if (encode_database) {
      EncodeDatabaseRequest encode_request =
          CreateEncodeDatabaseRequest(database);
      absl::StatusOr<EncodeDatabaseResponse> encode_response =
          EncodeDatabase(encode_request);
      if (!encode_response.ok()) {
        return encode_response.status();
      }
      *request.mutable_encoded_database()->mutable_shards() =
          std::move(encode_response)->shards();
    } else {
      for (const RawShard& raw_shard : database) {
        RawDatabaseShard* added_shard =
            request.mutable_raw_database()->add_shards();
        added_shard->set_shard_index(raw_shard.shard_id);
        for (const auto& [bucket_id, contents] : raw_shard.buckets) {
          RawDatabaseShard::Bucket* bucket = added_shard->add_buckets();
          bucket->set_bucket_id(bucket_id);
          bucket->set_bucket_contents(contents);
        }
      }
    }

    return request;
  }

  std::unique_ptr<const Context> request_context_;
  std::unique_ptr<const Context> response_context_;
  std::unique_ptr<const rlwe::SymmetricRlweKey<ModularInt>>
      request_private_key_;
  std::unique_ptr<const rlwe::SymmetricRlweKey<ModularInt>>
      response_private_key_;
  GenerateKeysResponse keys_;
};

TEST_F(ServerTest, EncodingDatabases) {
  // List of raw shards.
  const std::vector<RawShard> raw_shards = {
      {0, {{0, "hello"}, {17, "goodbye"}, {2, "testing"}}},
      {12, {{39, std::string(6000, 'B')}}},
      {7, {{7, ""}}},
  };

  EncodeDatabaseRequest request = CreateEncodeDatabaseRequest(raw_shards);
  absl::StatusOr<EncodeDatabaseResponse> response = EncodeDatabase(request);
  ASSERT_TRUE(response.ok());

  ASSERT_EQ(response->shards_size(), raw_shards.size());

  std::vector<RawShard> decoded_shards;
  decoded_shards.reserve(response->shards_size());
  int index = 0;
  for (const EncodedDatabaseShard& shard : response->shards()) {
    ASSERT_EQ(shard.buckets_size(), raw_shards[index].buckets.size())
        << " at index " << index;
    int shard_index = shard.shard_index();
    std::vector<std::pair<int, std::string>> decoded_buckets;
    decoded_buckets.reserve(shard.buckets_size());
    for (const EncodedDatabaseShard::Bucket& encoded_bucket : shard.buckets()) {
      std::vector<rlwe::Polynomial<ModularInt>> encoded_polynomials;
      encoded_polynomials.reserve(encoded_bucket.chunks_size());
      for (const auto& chunk : encoded_bucket.chunks()) {
        absl::StatusOr<rlwe::Polynomial<ModularInt>> polynomial =
            rlwe::Polynomial<ModularInt>::Deserialize(
                chunk, request_context_->GetModulusParams());
        ASSERT_TRUE(polynomial.ok()) << polynomial.status();
        encoded_polynomials.push_back(*std::move(polynomial));
      }
      absl::StatusOr<std::string> decoded_bucket =
          DecodePolynomial(*request_context_, encoded_polynomials);
      ASSERT_TRUE(decoded_bucket.ok()) << decoded_bucket.status();
      decoded_buckets.push_back({
          encoded_bucket.bucket_id(),
          *std::move(decoded_bucket),
      });
    }
    decoded_shards.push_back({
        .shard_id = shard_index,
        .buckets = decoded_buckets,
    });
    ++index;
  }

  EXPECT_EQ(decoded_shards, raw_shards);
}

TEST_F(ServerTest, EncodingEmptyDatabaseList) {
  // Empty list of shards.
  const std::vector<RawShard> raw_shards = {};

  EncodeDatabaseRequest request = CreateEncodeDatabaseRequest(raw_shards);
  absl::StatusOr<EncodeDatabaseResponse> response = EncodeDatabase(request);
  ASSERT_TRUE(response.ok());

  EXPECT_EQ(response->shards_size(), 0);
}

TEST_F(ServerTest, EncodingEmptyShard) {
  // Shard list with empty shard.
  const std::vector<RawShard> raw_shards = {{0, {}}};

  EncodeDatabaseRequest request = CreateEncodeDatabaseRequest(raw_shards);
  absl::StatusOr<EncodeDatabaseResponse> response = EncodeDatabase(request);
  EXPECT_FALSE(response.ok());
}

class ApplyQueriesTest : public ServerTest,
                         public testing::WithParamInterface<bool> {
 protected:
  void RunTest(const std::vector<PlaintextQuery>& queries,
               bool finalize_results) {
    RunTest(queries, CreateTestDatabase(), finalize_results);
  }

  void RunTest(const std::vector<PlaintextQuery>& queries,
               const std::vector<RawShard>& database, bool finalize_results) {
    absl::flat_hash_map<int, std::string> expected_results;
    for (const PlaintextQuery& query : queries) {
      int query_id = query.query_metadata().query_id();
      int shard_id = query.query_metadata().shard_id();
      int bucket_id = query.bucket_id();
      for (const RawShard& shard : database) {
        if (shard.shard_id == shard_id) {
          for (const auto& bucket : shard.buckets) {
            if (bucket.first == bucket_id) {
              expected_results[query_id] = bucket.second;
            }
          }
        }
      }
    }

    absl::StatusOr<ApplyQueriesRequest> request = CreateApplyQueriesRequest(
        queries, database, /*encode_database=*/GetParam());
    ASSERT_TRUE(request.ok());
    request->set_finalize_results(finalize_results);

    absl::StatusOr<ApplyQueriesResponse> response = ApplyQueries(*request);
    ASSERT_TRUE(response.ok());
    ASSERT_EQ(response->query_results_size(), expected_results.size());

    DecryptQueriesRequest decrypt_request;
    *decrypt_request.mutable_parameters() = CreateTestParameters();
    *decrypt_request.mutable_private_key() = keys_.private_key();
    *decrypt_request.mutable_public_key() = keys_.public_key();

    if (finalize_results) {
      *decrypt_request.mutable_encrypted_queries() = response->query_results();
    } else {
      FinalizeResultsRequest finalize_request;
      *finalize_request.mutable_parameters() = CreateTestParameters();
      *finalize_request.mutable_encrypted_results() = response->query_results();
      absl::StatusOr<FinalizeResultsResponse> finalize_response =
          FinalizeResults(finalize_request);
      ASSERT_TRUE(finalize_response.ok());
      ASSERT_EQ(finalize_response->finalized_results_size(), queries.size());
      *decrypt_request.mutable_encrypted_queries() =
          finalize_response->finalized_results();
    }

    absl::StatusOr<DecryptQueriesResponse> decrypt_response =
        DecryptQueries(decrypt_request);
    ASSERT_TRUE(decrypt_response.ok());
    ASSERT_EQ(decrypt_response->result_size(), queries.size());

    absl::flat_hash_map<int, std::string> results;
    for (const auto& result : decrypt_response->result()) {
      int query_id = result.query_metadata().query_id();
      results[query_id] = result.result();
    }

    EXPECT_THAT(results, UnorderedElementsAreArray(expected_results));
  }
};

TEST_P(ApplyQueriesTest, ApplyQueries) {
  const std::vector<PlaintextQuery> queries = {
      CreatePlaintextQuery(/*query_id=*/0, /*shard_id=*/0, /*bucket_index=*/0),
      CreatePlaintextQuery(/*query_id=*/1, /*shard_id=*/1, /*bucket_index=*/1),
      CreatePlaintextQuery(/*query_id=*/2, /*shard_id=*/1,
                           /*bucket_index=*/kEmptyBucketId),
  };

  RunTest(queries, /*finalize_results=*/false);
}

TEST_P(ApplyQueriesTest, ApplyQueriesWithFinalizing) {
  const std::vector<PlaintextQuery> queries = {
      CreatePlaintextQuery(/*query_id=*/10, /*shard_id=*/0,
                           /*bucket_index=*/0),
      CreatePlaintextQuery(/*query_id=*/110,
                           /*shard_id=*/1, /*bucket_index=*/1),
      CreatePlaintextQuery(/*query_id=*/102, /*shard_id=*/1,
                           /*bucket_index=*/kEmptyBucketId),
  };

  RunTest(queries, /*finalize_results=*/true);
}

TEST_P(ApplyQueriesTest, ApplyQueriesWithNoPublicKeys) {
  const std::vector<PlaintextQuery> queries = {
      CreatePlaintextQuery(/*query_id=*/0, /*shard_id=*/0, /*bucket_index=*/0),
  };
  absl::StatusOr<ApplyQueriesRequest> request =
      CreateApplyQueriesRequest(queries, /*encode_database=*/GetParam());
  ASSERT_TRUE(request.ok());

  request->clear_public_key();

  absl::StatusOr<ApplyQueriesResponse> response = ApplyQueries(*request);
  EXPECT_EQ(response.status().code(), absl::StatusCode::kInvalidArgument);
}

TEST_P(ApplyQueriesTest, ApplyQueriesWithNoQueries) {
  const std::vector<PlaintextQuery> queries = {
      CreatePlaintextQuery(/*query_id=*/0, /*shard_id=*/0, /*bucket_index=*/0),
  };
  absl::StatusOr<ApplyQueriesRequest> request =
      CreateApplyQueriesRequest(queries, /*encode_database=*/GetParam());
  ASSERT_TRUE(request.ok());

  request->clear_queries();

  absl::StatusOr<ApplyQueriesResponse> response = ApplyQueries(*request);
  EXPECT_EQ(response.status().code(), absl::StatusCode::kInvalidArgument);
}

TEST_P(ApplyQueriesTest, ApplyQueriesWithNoShards) {
  const std::vector<PlaintextQuery> queries = {
      CreatePlaintextQuery(/*query_id=*/0, /*shard_id=*/0, /*bucket_index=*/0),
  };
  absl::StatusOr<ApplyQueriesRequest> request =
      CreateApplyQueriesRequest(queries, /*encode_database=*/GetParam());
  ASSERT_TRUE(request.ok());

  request->clear_encoded_database();
  request->clear_raw_database();

  absl::StatusOr<ApplyQueriesResponse> response = ApplyQueries(*request);
  EXPECT_EQ(response.status().code(), absl::StatusCode::kInvalidArgument);
}

TEST_P(ApplyQueriesTest, ApplyQueriesWithMissingShard) {
  const std::vector<RawShard> raw_shards = {
      {0, {{0, "hello"}, {1, "goodbye"}}},
  };
  // Send a query for shard index 1 that is not specified in the input.
  const std::vector<PlaintextQuery> queries = {
      CreatePlaintextQuery(/*query_id=*/0, /*shard_id=*/0, /*bucket_index=*/0),
      CreatePlaintextQuery(/*query_id=*/1, /*shard_id=*/1, /*bucket_index=*/0),
  };
  absl::StatusOr<ApplyQueriesRequest> request = CreateApplyQueriesRequest(
      queries, raw_shards, /*encode_database=*/GetParam());
  ASSERT_TRUE(request.ok());

  absl::StatusOr<ApplyQueriesResponse> response = ApplyQueries(*request);
  ASSERT_TRUE(response.ok());

  EXPECT_EQ(response->query_results_size(), 1);
}

TEST_P(ApplyQueriesTest, ApplyQueriesWithDuplicateBucket) {
  // Specify a database with multiple definitions of a single bucket.
  const std::vector<RawShard> raw_shards = {
      {0, {{0, "hello"}, {0, "goodbye"}}},
  };
  const std::vector<PlaintextQuery> queries = {
      CreatePlaintextQuery(/*query_id=*/0, /*shard_id=*/0, /*bucket_index=*/0),
  };
  absl::StatusOr<ApplyQueriesRequest> request = CreateApplyQueriesRequest(
      queries, raw_shards, /*encode_database=*/GetParam());
  ASSERT_TRUE(request.ok());

  absl::StatusOr<ApplyQueriesResponse> response = ApplyQueries(*request);
  EXPECT_EQ(response.status().code(), absl::StatusCode::kInvalidArgument);
}

TEST_P(ApplyQueriesTest, ApplyQueriesWithShardSubset) {
  // Specify a database with one bucket per shard.
  const std::vector<RawShard> raw_shards = {
      {0, {{17, "hello"}}},
      {1, {{97, "goodbye"}}},
  };

  const std::vector<PlaintextQuery> queries = {
      CreatePlaintextQuery(/*query_id=*/0, /*shard_id=*/0, /*bucket_index=*/17),
      CreatePlaintextQuery(/*query_id=*/1, /*shard_id=*/1, /*bucket_index=*/97),
  };

  RunTest(queries, raw_shards, /*finalize_results=*/true);
}

INSTANTIATE_TEST_SUITE_P(ApplyQueriesTest, ApplyQueriesTest, ::testing::Bool());

TEST_F(ServerTest, SumCiphertexts) {
  // Note each coefficient enables homomorphic additions mod 3 with the test
  // parameters.
  const std::vector<std::vector<Int>> addends = {
      {1, 1, 0}, {0, 1}, {1}, {0}, {0, 1, 0, 0}, {1, 1, 1}, {2, 2}};
  const std::vector<Int> expected_sums = {2, 1, 1, 0, 1, 0, 1};

  absl::StatusOr<SumCiphertextsRequest> request = CreateSumCiphertextsRequest(
      addends, *request_context_, *request_private_key_);
  ASSERT_TRUE(request.ok());

  absl::StatusOr<SumCiphertextsResponse> response = SumCiphertexts(*request);
  ASSERT_TRUE(response.ok());

  ASSERT_EQ(response->ciphertexts_size(), expected_sums.size());

  int index = 0;
  for (const rlwe::SerializedSymmetricRlweCiphertext& ciphertext :
       response->ciphertexts()) {
    absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> deserialized =
        rlwe::SymmetricRlweCiphertext<ModularInt>::Deserialize(
            ciphertext, request_context_->GetModulusParams(),
            request_context_->GetErrorParams());
    ASSERT_TRUE(deserialized.ok()) << deserialized.status();
    absl::StatusOr<Int> decrypted = DecryptSingleInt(
        *request_context_, *request_private_key_, *deserialized);
    ASSERT_TRUE(decrypted.ok()) << decrypted.status();
    EXPECT_EQ(*decrypted, expected_sums[index]) << "at index " << index;
    ++index;
  }
}

TEST_F(ServerTest, SumEmptyList) {
  // An empty list input should return an empty list.
  const std::vector<std::vector<Int>> addends = {};

  absl::StatusOr<SumCiphertextsRequest> request = CreateSumCiphertextsRequest(
      addends, *request_context_, *request_private_key_);
  ASSERT_TRUE(request.ok());

  absl::StatusOr<SumCiphertextsResponse> response = SumCiphertexts(*request);
  ASSERT_TRUE(response.ok());

  EXPECT_EQ(response->ciphertexts_size(), 0);
}

TEST_F(ServerTest, SumEmptySublist) {
  // An empty sublist will return an error.
  const std::vector<std::vector<Int>> addends = {{}, {1}};

  absl::StatusOr<SumCiphertextsRequest> request = CreateSumCiphertextsRequest(
      addends, *request_context_, *request_private_key_);
  ASSERT_TRUE(request.ok());

  absl::StatusOr<SumCiphertextsResponse> response = SumCiphertexts(*request);
  EXPECT_FALSE(response.ok());
}

TEST_F(ServerTest, FinalizeResults) {
  // List of list of results. Each result will correspond to a sublist.
  const std::vector<std::vector<Int>> results = {
      {2, 1, 2, 0, 1, 2}, {1, 2}, {0}, {0, 0}};

  absl::StatusOr<FinalizeResultsRequest> request = CreateFinalizeResultsRequest(
      results, *request_context_, *request_private_key_);
  ASSERT_TRUE(request.ok());

  absl::StatusOr<FinalizeResultsResponse> response = FinalizeResults(*request);
  ASSERT_TRUE(response.ok());

  ASSERT_EQ(response->finalized_results_size(), results.size());

  int result_index = 0;
  for (const EncryptedQueryResult& result : response->finalized_results()) {
    int result_sublist_index = 0;
    ASSERT_EQ(result.ciphertexts_size(), results[result_index].size())
        << " at result list " << result_index;
    for (const rlwe::SerializedSymmetricRlweCiphertext& ciphertext :
         result.ciphertexts()) {
      absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> deserialized =
          rlwe::SymmetricRlweCiphertext<ModularInt>::Deserialize(
              ciphertext, response_context_->GetModulusParams(),
              response_context_->GetErrorParams());
      ASSERT_TRUE(deserialized.ok()) << deserialized.status();

      absl::StatusOr<Int> decrypted = DecryptSingleInt(
          *response_context_, *response_private_key_, *deserialized);
      ASSERT_TRUE(decrypted.ok()) << decrypted.status();
      EXPECT_EQ(*decrypted, results[result_index][result_sublist_index])
          << "at ciphertext " << result_sublist_index << " in result sub-list "
          << result_index;
      ++result_sublist_index;
    }
    ++result_index;
  }
}

TEST_F(ServerTest, FinalizeEmptyResultList) {
  // Empty list of results.
  const std::vector<std::vector<Int>> results = {};

  absl::StatusOr<FinalizeResultsRequest> request = CreateFinalizeResultsRequest(
      results, *request_context_, *request_private_key_);
  ASSERT_TRUE(request.ok());

  absl::StatusOr<FinalizeResultsResponse> response = FinalizeResults(*request);
  ASSERT_TRUE(response.ok());

  EXPECT_EQ(response->finalized_results_size(), 0);
}

TEST_F(ServerTest, FinalizeEmptyResultSublist) {
  // List of results with an invalid empty sublist.
  const std::vector<std::vector<Int>> results = {{1}, {}};

  absl::StatusOr<FinalizeResultsRequest> request = CreateFinalizeResultsRequest(
      results, *request_context_, *request_private_key_);
  ASSERT_TRUE(request.ok());

  absl::StatusOr<FinalizeResultsResponse> response = FinalizeResults(*request);
  EXPECT_FALSE(response.ok());
}

}  // namespace
}  // namespace batch
}  // namespace private_membership
