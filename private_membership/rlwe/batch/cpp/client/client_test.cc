// Copyright 2021 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "private_membership/rlwe/batch/cpp/client/client.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <google/protobuf/repeated_field.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "private_membership/rlwe/batch/proto/client.pb.h"
#include "private_membership/rlwe/batch/cpp/client/client_helper.h"
#include "private_membership/rlwe/batch/cpp/constants.h"
#include "private_membership/rlwe/batch/cpp/context.h"
#include "private_membership/rlwe/batch/cpp/encoding.h"
#include "private_membership/rlwe/batch/cpp/padding.h"
#include "private_membership/rlwe/batch/cpp/prng.h"
#include "private_membership/rlwe/batch/proto/shared.pb.h"
#include "private_membership/rlwe/batch/cpp/test_helper.h"
#include "shell_encryption/context.h"
#include "shell_encryption/galois_key.h"
#include "shell_encryption/int256.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/oblivious_expand.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/serialization.pb.h"
#include "shell_encryption/symmetric_encryption.h"
#include "shell_encryption/symmetric_encryption_with_prng.h"
#include "shell_encryption/transcription.h"

namespace private_membership {
namespace batch {
namespace {

// With 2 levels of recursion, the number of slots per level is
// ceil(200^{1/2}) = 15.
constexpr int kSlotsPerLevel = 15;

constexpr std::array<int, 3> kTestQueryIds = {0, 1, 2};
constexpr std::array<int, 3> kTestShardIds = {1, 3, 2};
constexpr std::array<int, 3> kTestBucketIds = {173, 14, 82};
constexpr std::array<absl::string_view, 3> kTestBuckets = {"abc", "defghijk",
                                                           "lmnop"};

// Note that the following are true:
// - 143 = 8 + (15 * 11)
// - 14 = 14 + (15 * 0)
// - 82 = 7 + (15 * 5)
constexpr std::array<int, 6> kExpectedQueryIndices = {8,
                                                      kSlotsPerLevel + 11,
                                                      kSlotsPerLevel * 2 + 14,
                                                      kSlotsPerLevel * 3 + 0,
                                                      kSlotsPerLevel * 4 + 7,
                                                      kSlotsPerLevel * 5 + 5};

// Each of the 3 queries requires 15 slots per level of recursion.
constexpr int kExpectedNumberOfSlots = 3 * kSlotsPerLevel * kLevelsOfRecursion;

// Avoid EqualsProto that is not available with lite protos.
MATCHER_P(IsQueryMetadata, other, "") {
  return testing::ExplainMatchResult(
      testing::AllOf(testing::Property("query_id", &QueryMetadata::query_id,
                                       other.query_id()),
                     testing::Property("shard_id", &QueryMetadata::shard_id,
                                       other.shard_id())),
      arg, result_listener);
}

GenerateKeysRequest CreateTestGenerateKeysRequest() {
  GenerateKeysRequest request;
  *(request.mutable_parameters()) = CreateTestParameters();
  return request;
}

std::vector<PlaintextQuery> CreateTestPlaintextQueries() {
  std::vector<PlaintextQuery> test_plaintext_queries(kTestQueryIds.size());

  for (int i = 0; i < kTestQueryIds.size(); ++i) {
    // Set query metadata.
    QueryMetadata* query_metadata =
        test_plaintext_queries[i].mutable_query_metadata();
    query_metadata->set_query_id(kTestQueryIds.at(i));
    query_metadata->set_shard_id(kTestShardIds.at(i));

    // Set bucket ID.
    test_plaintext_queries[i].set_bucket_id(kTestBucketIds.at(i));
  }

  return test_plaintext_queries;
}

absl::StatusOr<EncryptQueriesRequest> CreateTestEncryptQueriesRequest() {
  GenerateKeysRequest keys_request = CreateTestGenerateKeysRequest();

  absl::StatusOr<GenerateKeysResponse> keys_response =
      GenerateKeys(keys_request);
  if (!keys_response.ok()) {
    return keys_response.status();
  }

  EncryptQueriesRequest request;
  *request.mutable_parameters() = CreateTestParameters();
  *request.mutable_private_key() = keys_response->private_key();
  *request.mutable_public_key() = keys_response->public_key();

  std::vector<PlaintextQuery> test_plaintext_queries =
      CreateTestPlaintextQueries();
  *request.mutable_plaintext_queries() = {test_plaintext_queries.begin(),
                                          test_plaintext_queries.end()};

  return request;
}

absl::StatusOr<DecryptQueriesRequest> CreateTestDecryptQueriesRequest() {
  GenerateKeysRequest keys_request = CreateTestGenerateKeysRequest();

  absl::StatusOr<GenerateKeysResponse> keys_response =
      GenerateKeys(keys_request);
  if (!keys_response.ok()) {
    return keys_response.status();
  }

  DecryptQueriesRequest request;
  *request.mutable_parameters() = CreateTestParameters();
  *request.mutable_private_key() = keys_response->private_key();
  *request.mutable_public_key() = keys_response->public_key();

  auto response_context = CreateRlweResponseContext(request.parameters());
  if (!response_context.ok()) {
    return response_context.status();
  }
  auto private_key = DeserializePrivateKey(request.private_key().response_key(),
                                           **response_context);
  if (!private_key.ok()) {
    return private_key.status();
  }

  std::vector<PlaintextQuery> test_plaintext_queries =
      CreateTestPlaintextQueries();
  for (int i = 0; i < test_plaintext_queries.size(); ++i) {
    EncryptedQueryResult* encrypted_result = request.add_encrypted_queries();
    *encrypted_result->mutable_query_metadata() =
        test_plaintext_queries[i].query_metadata();
    absl::StatusOr<std::vector<rlwe::Polynomial<ModularInt>>> ntt_polynomials =
        EncodeBytes(**response_context, kTestBuckets[i]);
    if (!ntt_polynomials.ok()) {
      return ntt_polynomials.status();
    }

    for (const rlwe::Polynomial<ModularInt>& ntt_polynomial :
         *ntt_polynomials) {
      auto prng = CreatePrng();
      if (!prng.ok()) {
        return prng.status();
      }
      absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> encryption =
          rlwe::Encrypt(*private_key, ntt_polynomial,
                        (*response_context)->GetErrorParams(), prng->get());
      if (!encryption.ok()) {
        return encryption.status();
      }
      absl::StatusOr<rlwe::SerializedSymmetricRlweCiphertext>
          serialized_encryption = encryption->Serialize();
      if (!serialized_encryption.ok()) {
        return serialized_encryption.status();
      }
      *encrypted_result->add_ciphertexts() = *std::move(serialized_encryption);
    }
  }

  return request;
}

TEST(ClientTest, GeneratesKeys) {
  GenerateKeysRequest request = CreateTestGenerateKeysRequest();

  absl::StatusOr<GenerateKeysResponse> response = GenerateKeys(request);
  ASSERT_TRUE(response.ok()) << response.status();

  auto request_context = CreateRlweRequestContext(request.parameters());
  ASSERT_TRUE(request_context.ok());
  ASSERT_NE(request_context->get(), nullptr);
  auto response_context = CreateRlweResponseContext(request.parameters());
  ASSERT_TRUE(response_context.ok());
  ASSERT_NE(response_context->get(), nullptr);

  EXPECT_TRUE(DeserializePrivateKey(response->private_key().request_key(),
                                    **request_context)
                  .ok());
  EXPECT_TRUE(DeserializePrivateKey(response->private_key().response_key(),
                                    **response_context)
                  .ok());
  for (int i = 0; i < response->public_key().galois_key_size(); ++i) {
    EXPECT_TRUE(DeserializePublicGaloisKey(response->public_key().galois_key(i),
                                           **request_context)
                    .ok());
  }
}

TEST(ClientTest, RandomKeys) {
  GenerateKeysRequest request = CreateTestGenerateKeysRequest();

  int num_responses = 10;
  std::vector<GenerateKeysResponse> responses(num_responses);
  for (int i = 0; i < num_responses; ++i) {
    absl::StatusOr<GenerateKeysResponse> response = GenerateKeys(request);
    ASSERT_TRUE(response.ok());
    responses[i] = *std::move(response);
  }

  auto request_context = CreateRlweRequestContext(request.parameters());
  ASSERT_TRUE(request_context.ok());
  ASSERT_NE(request_context->get(), nullptr);
  auto response_context = CreateRlweResponseContext(request.parameters());
  ASSERT_TRUE(response_context.ok());
  ASSERT_NE(response_context->get(), nullptr);

  for (int i = 0; i < num_responses; ++i) {
    EXPECT_NE(responses[i].private_key().request_key().coeffs(),
              responses[i].private_key().response_key().coeffs());
  }

  // No need to check public keys as long as private keys are random.
  for (int i = 0; i < num_responses; ++i) {
    for (int j = 0; j < i; ++j) {
      EXPECT_NE(responses[i].private_key().request_key().coeffs(),
                responses[j].private_key().request_key().coeffs());
      EXPECT_NE(responses[i].private_key().response_key().coeffs(),
                responses[j].private_key().response_key().coeffs());
    }
  }
}

TEST(ClientTest, CheckPrivateKeyRelation) {
  GenerateKeysRequest request = CreateTestGenerateKeysRequest();

  absl::StatusOr<GenerateKeysResponse> response = GenerateKeys(request);
  ASSERT_TRUE(response.ok());

  auto request_context = CreateRlweRequestContext(request.parameters());
  ASSERT_TRUE(request_context.ok());
  ASSERT_NE(request_context->get(), nullptr);
  auto response_context = CreateRlweResponseContext(request.parameters());
  ASSERT_TRUE(response_context.ok());
  ASSERT_NE(response_context->get(), nullptr);

  auto request_key = DeserializePrivateKey(
      response->private_key().request_key(), **request_context);
  ASSERT_TRUE(request_key.ok());
  auto response_key = DeserializePrivateKey(
      response->private_key().response_key(), **response_context);
  ASSERT_TRUE(response_key.ok());

  std::vector<Int> values = {0, 1, 2, 1, 0, 2};
  for (int i = 0; i < values.size(); ++i) {
    absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> ciphertext =
        EncryptSingleInt(**request_context, *request_key, values[i]);
    ASSERT_TRUE(ciphertext.ok()) << ciphertext.status();
    absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>>
        switched_ciphertext = ciphertext->SwitchModulus(
            (*request_context)->GetNttParams(),
            (*response_context)->GetModulusParams(),
            (*response_context)->GetNttParams(),
            (*response_context)->GetErrorParams(), (*response_context)->GetT());
    ASSERT_TRUE(switched_ciphertext.ok()) << switched_ciphertext.status();
    absl::StatusOr<Int> decrypted = DecryptSingleInt(
        **response_context, *response_key, *switched_ciphertext);
    ASSERT_TRUE(decrypted.ok()) << decrypted.status();
    EXPECT_EQ(*decrypted, values[i]) << " at index i";
  }
}

TEST(ClientTest, GenerateBitVector) {
  absl::StatusOr<EncryptQueriesRequest> encrypt_request =
      CreateTestEncryptQueriesRequest();
  ASSERT_TRUE(encrypt_request.ok());

  absl::StatusOr<BitVector> bit_vector = GenerateBitVector(*encrypt_request);
  ASSERT_TRUE(bit_vector.ok());

  // With 2 levels of recursion, ceil(200^{0.5}) = 15. Each level uses 15 slots.
  // So, the total number of slots is 15 * 2 * number of queries.
  EXPECT_EQ(bit_vector->total_size,
            kLevelsOfRecursion * kSlotsPerLevel * kTestQueryIds.size());

  EXPECT_EQ(bit_vector->indices, std::vector<int>(kExpectedQueryIndices.begin(),
                                                  kExpectedQueryIndices.end()));
}

TEST(ClientTest, EncryptQueries) {
  absl::StatusOr<EncryptQueriesRequest> encrypt_request =
      CreateTestEncryptQueriesRequest();
  ASSERT_TRUE(encrypt_request.ok());

  absl::StatusOr<EncryptQueriesResponse> encrypt_response =
      EncryptQueries(*encrypt_request);
  ASSERT_TRUE(encrypt_response.ok());

  // Check that metadata is populated correctly.
  auto test_plaintext_queries = CreateTestPlaintextQueries();
  auto encrypted_queries = encrypt_response->encrypted_queries();
  ASSERT_EQ(encrypted_queries.query_metadata_size(),
            test_plaintext_queries.size());
  for (int i = 0; i < test_plaintext_queries.size(); ++i) {
    EXPECT_THAT(test_plaintext_queries[i].query_metadata(),
                IsQueryMetadata(encrypted_queries.query_metadata(i)));
  }

  // Generate necessary contexts, PRNGs and private keys.
  auto prng = CreatePrngFromSeed(encrypted_queries.prng_seed());
  ASSERT_TRUE(prng.ok());
  auto request_context =
      CreateRlweRequestContext(encrypt_request->parameters());
  ASSERT_TRUE(request_context.ok());
  ASSERT_NE(request_context->get(), nullptr);
  auto private_key = DeserializePrivateKey(
      encrypt_request->private_key().request_key(), **request_context);
  ASSERT_TRUE(private_key.ok());

  // Deserialize encrypted and compressed queries.
  int num_encrypted_queries = encrypted_queries.encrypted_queries_size();
  std::vector<rlwe::Polynomial<ModularInt>> deserialized_ciphertexts;
  deserialized_ciphertexts.reserve(num_encrypted_queries);
  for (int i = 0; i < num_encrypted_queries; ++i) {
    auto deserialized = rlwe::Polynomial<ModularInt>::Deserialize(
        encrypted_queries.encrypted_queries(i),
        (*request_context)->GetModulusParams());
    ASSERT_TRUE(deserialized.ok());
    deserialized_ciphertexts.push_back(*std::move(deserialized));
  }

  // Expand queries using PRNG.
  auto expanded_ciphertexts = rlwe::ExpandFromPrng<ModularInt>(
      std::move(deserialized_ciphertexts),
      (*request_context)->GetModulusParams(),
      (*request_context)->GetNttParams(), (*request_context)->GetErrorParams(),
      prng->get());
  ASSERT_TRUE(expanded_ciphertexts.ok());

  // Deserialize public key for oblivious expansion.
  int num_galois_keys = encrypt_request->public_key().galois_key_size();
  std::vector<rlwe::GaloisKey<ModularInt>> galois_keys;
  galois_keys.reserve(num_galois_keys);
  for (const auto& serialized_galois_key :
       encrypt_request->public_key().galois_key()) {
    auto galois_key = rlwe::GaloisKey<ModularInt>::Deserialize(
        serialized_galois_key, (*request_context)->GetModulusParams(),
        (*request_context)->GetNttParams());
    ASSERT_TRUE(galois_key.ok());
    galois_keys.push_back(*std::move(galois_key));
  }
  auto expander = rlwe::GaloisKeysObliviousExpander<ModularInt>::Create(
      std::move(galois_keys), (*request_context)->GetLogT(),
      (*request_context)->GetModulusParams(),
      (*request_context)->GetNttParams());
  ASSERT_TRUE(expander.ok());

  // Obliviously expand ciphertexts.
  auto oblivious_expanded_ciphertexts =
      (*expander)->ObliviousExpand(*std::move(expanded_ciphertexts),
                                   num_galois_keys, kExpectedNumberOfSlots);
  ASSERT_EQ(oblivious_expanded_ciphertexts->size(), kExpectedNumberOfSlots);

  // Decrypt and transcribe plaintexts.
  std::vector<ModularInt::Int> decrypted_queries;
  for (int i = 0; i < oblivious_expanded_ciphertexts->size(); ++i) {
    auto decrypted_query =
        rlwe::Decrypt(*private_key, oblivious_expanded_ciphertexts->at(i));
    ASSERT_TRUE(decrypted_query.ok());
    absl::StatusOr<std::vector<uint8_t>> bytes =
        rlwe::TranscribeBits<ModularInt::Int, uint8_t>(
            *std::move(decrypted_query),
            private_key->Len() * private_key->BitsPerCoeff(),
            private_key->BitsPerCoeff(), 8);
    ASSERT_TRUE(bytes.ok());
    ASSERT_FALSE(bytes->empty());
    ASSERT_TRUE(bytes->at(0) == 0 || bytes->at(0) == 1);
    // All bytes beyond the lowest order byte should be 0.
    for (int i = 1; i < bytes->size(); ++i) {
      EXPECT_EQ(bytes->at(i), 0);
    }
    decrypted_queries.push_back(bytes->at(0));
  }

  // Check that all 1-indices exist as expected after decryption.
  absl::StatusOr<BitVector> bit_vector = GenerateBitVector(*encrypt_request);
  ASSERT_TRUE(bit_vector.ok());
  const std::vector<int>& one_indices = bit_vector->indices;
  for (int i = 0; i < kExpectedNumberOfSlots; ++i) {
    if (std::find(one_indices.begin(), one_indices.end(), i) !=
        one_indices.end()) {
      EXPECT_EQ(decrypted_queries[i], 1);
    } else {
      EXPECT_EQ(decrypted_queries[i], 0);
    }
  }
}

TEST(ClientTest, RandomEncryptedQueries) {
  int num_responses = 10;
  std::vector<EncryptQueriesResponse> responses(num_responses);
  for (int i = 0; i < num_responses; ++i) {
    absl::StatusOr<EncryptQueriesRequest> encrypt_request =
        CreateTestEncryptQueriesRequest();
    ASSERT_TRUE(encrypt_request.ok());

    absl::StatusOr<EncryptQueriesResponse> encrypt_response =
        EncryptQueries(*encrypt_request);
    ASSERT_TRUE(encrypt_response.ok());
    responses[i] = *std::move(encrypt_response);
  }

  for (int i = 0; i < num_responses; ++i) {
    for (int j = 0; j < i; ++j) {
      auto encrypted_query1 = responses[i].encrypted_queries();
      auto encrypted_query2 = responses[j].encrypted_queries();
      ASSERT_EQ(encrypted_query1.encrypted_queries_size(),
                encrypted_query2.encrypted_queries_size());
      int num_encrypted_queries = encrypted_query1.encrypted_queries_size();
      for (int k = 0; k < num_encrypted_queries; ++k) {
        EXPECT_NE(
            responses[i].encrypted_queries().encrypted_queries(k).coeffs(),
            responses[j].encrypted_queries().encrypted_queries(k).coeffs());
      }
    }
  }
}

TEST(ClientTest, InvalidBucketId) {
  absl::StatusOr<EncryptQueriesRequest> encrypt_request =
      CreateTestEncryptQueriesRequest();
  ASSERT_TRUE(encrypt_request.ok());

  EncryptQueriesRequest negative_bucket_id = *encrypt_request;
  negative_bucket_id.mutable_plaintext_queries(0)->set_bucket_id(-1);
  EXPECT_FALSE(EncryptQueries(negative_bucket_id).ok());

  EncryptQueriesRequest large_bucket_id = *encrypt_request;
  large_bucket_id.mutable_plaintext_queries(0)->set_bucket_id(
      kNumberOfBucketsPerShard);
  EXPECT_FALSE(EncryptQueries(large_bucket_id).ok());
}

TEST(ClientTest, InvalidShardId) {
  absl::StatusOr<EncryptQueriesRequest> encrypt_request =
      CreateTestEncryptQueriesRequest();
  ASSERT_TRUE(encrypt_request.ok());

  EncryptQueriesRequest negative_shard_id = *encrypt_request;
  negative_shard_id.mutable_plaintext_queries(0)
      ->mutable_query_metadata()
      ->set_shard_id(-1);
  EXPECT_FALSE(EncryptQueries(negative_shard_id).ok());

  EncryptQueriesRequest large_shard_id = *encrypt_request;
  large_shard_id.mutable_plaintext_queries(0)
      ->mutable_query_metadata()
      ->set_shard_id(kNumberOfShards);
  EXPECT_FALSE(EncryptQueries(large_shard_id).ok());
}

TEST(ClientTest, DecryptQueries) {
  absl::StatusOr<DecryptQueriesRequest> request =
      CreateTestDecryptQueriesRequest();
  ASSERT_TRUE(request.ok());

  absl::StatusOr<DecryptQueriesResponse> response = DecryptQueries(*request);
  ASSERT_TRUE(response.ok());

  auto test_plaintext_queries = CreateTestPlaintextQueries();
  ASSERT_EQ(response->result_size(), test_plaintext_queries.size());
  for (int i = 0; i < test_plaintext_queries.size(); ++i) {
    // Check that metadata is populated correctly.
    EXPECT_THAT(test_plaintext_queries[i].query_metadata(),
                IsQueryMetadata(response->result(i).query_metadata()));

    // Check the result is correct.
    EXPECT_EQ(response->result(i).result(), kTestBuckets[i]);
  }
}

}  // namespace
}  // namespace batch
}  // namespace private_membership
