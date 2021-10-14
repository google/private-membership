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

#include <iostream>
#include <string>
#include <vector>

#include <google/protobuf/repeated_field.h>
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/time/clock.h"
#include "absl/time/time.h"
#include "private_membership/rlwe/batch/cpp/client/client.h"
#include "private_membership/rlwe/batch/cpp/server/server.h"
#include "private_membership/rlwe/batch/proto/client.pb.h"
#include "private_membership/rlwe/batch/proto/server.pb.h"
#include "private_membership/rlwe/batch/proto/shared.pb.h"

namespace private_membership {
namespace batch {
namespace {

// A 20,000 bucket database split into 5 shards. Each bucket is 400 bytes.
constexpr int kNumberOfShards = 5;
constexpr int kNumberOfBucketsPerShard = 4000;
constexpr int kContentByteLength = 400;

// Create parameters for example.
Parameters CreateParameters() {
  Parameters parameters;

  Parameters::ShardParameters* shard_parameters =
      parameters.mutable_shard_parameters();
  shard_parameters->set_number_of_shards(kNumberOfShards);
  shard_parameters->set_number_of_buckets_per_shard(kNumberOfBucketsPerShard);

  // Example test parameters for the underlying RLWE cryptosystem. The security
  // of these parameters may be calculated using the code found at:
  // https://bitbucket.org/malb/lwe-estimator/src/master/
  Parameters::CryptoParameters* crypto_parameters =
      parameters.mutable_crypto_parameters();
  crypto_parameters->add_request_modulus(18446744073708380161ULL);
  crypto_parameters->add_request_modulus(137438953471ULL);
  crypto_parameters->add_response_modulus(2056193ULL);
  crypto_parameters->set_log_degree(12);
  crypto_parameters->set_log_t(1);
  crypto_parameters->set_variance(8);
  crypto_parameters->set_levels_of_recursion(2);
  crypto_parameters->set_log_compression_factor(4);
  crypto_parameters->set_log_decomposition_modulus(10);

  return parameters;
}

std::string CreateDatabaseContents(int shard_id, int bucket_id) {
  // Every 700-th bucket will be empty.
  if (bucket_id % 700 == 6) {
    return std::string(kContentByteLength, '0');
  }

  // All odd buckets will be "shard_id,bucket_id".
  if (bucket_id % 2 == 1) {
    std::string bucket = absl::StrCat(shard_id, ",", bucket_id);
    bucket.append(kContentByteLength - bucket.size(), '0');
    return bucket;
  }

  // All even buckets will be "bucket_id:shard_id".
  std::string bucket = absl::StrCat(bucket_id, ":", shard_id);
  bucket.append(kContentByteLength - bucket.size(), '0');
  return bucket;
}

std::vector<RawDatabaseShard> CreateRawDatabase() {
  std::vector<RawDatabaseShard> database(kNumberOfShards);
  for (int shard_id = 0; shard_id < kNumberOfShards; ++shard_id) {
    database[shard_id].set_shard_index(shard_id);
    for (int bucket_id = 0; bucket_id < kNumberOfBucketsPerShard; ++bucket_id) {
      RawDatabaseShard::Bucket* bucket = database[shard_id].add_buckets();
      bucket->set_bucket_id(bucket_id);
      bucket->set_bucket_contents(CreateDatabaseContents(shard_id, bucket_id));
    }
  }
  return database;
}

PlaintextQuery CreatePlaintextQuery(int query_id, int shard_id, int bucket_id) {
  PlaintextQuery query;
  query.mutable_query_metadata()->set_query_id(query_id);
  query.mutable_query_metadata()->set_shard_id(shard_id);
  query.set_bucket_id(bucket_id);
  return query;
}

std::vector<PlaintextQuery> CreatePlaintextQueries() {
  std::vector<PlaintextQuery> queries = {
      // Empty buckets (bucket_id divisible by 7).
      CreatePlaintextQuery(/*query_id=*/0, /*shard_id=*/0, /*bucket_id=*/28),
      CreatePlaintextQuery(/*query_id=*/6, /*shard_id=*/3,
                           /*bucket_id=*/707),

      // Odd buckets.
      CreatePlaintextQuery(/*query_id=*/10, /*shard_id=*/4,
                           /*bucket_id=*/17),
      CreatePlaintextQuery(/*query_id=*/12, /*shard_id=*/3,
                           /*bucket_id=*/1971),
      CreatePlaintextQuery(/*query_id=*/13, /*shard_id=*/2,
                           /*bucket_id=*/333),
      CreatePlaintextQuery(/*query_id=*/14, /*shard_id=*/0,
                           /*bucket_id=*/971),
      CreatePlaintextQuery(/*query_id=*/19, /*shard_id=*/1,
                           /*bucket_id=*/1003),

      // Even buckets.
      CreatePlaintextQuery(/*query_id=*/21, /*shard_id=*/3,
                           /*bucket_id=*/20),
      CreatePlaintextQuery(/*query_id=*/22, /*shard_id=*/3,
                           /*bucket_id=*/20),
      CreatePlaintextQuery(/*query_id=*/23, /*shard_id=*/1,
                           /*bucket_id=*/20),
  };
  return queries;
}

absl::Status RunExample() {
  // Used to time operations.
  absl::Time start = absl::Now();
  absl::Time end = absl::Now();

  // Parameters for the system.
  Parameters parameters = CreateParameters();

  // Client generates public and private key pairs.
  std::cout << "Generating keys." << std::endl;
  start = absl::Now();
  GenerateKeysRequest keys_request;
  *keys_request.mutable_parameters() = parameters;
  absl::StatusOr<GenerateKeysResponse> response = GenerateKeys(keys_request);
  if (!response.ok()) {
    return response.status();
  }
  const PrivateKey& private_key = response->private_key();
  const PublicKey& public_key = response->public_key();
  end = absl::Now();
  std::cout << "Key generation: " << absl::ToInt64Milliseconds(end - start)
            << " ms" << std::endl;

  // Client generates request.
  std::cout << "Generating client request." << std::endl;
  start = absl::Now();
  EncryptQueriesRequest encrypt_request;
  *encrypt_request.mutable_parameters() = parameters;
  *encrypt_request.mutable_public_key() = public_key;
  *encrypt_request.mutable_private_key() = private_key;
  std::vector<PlaintextQuery> plaintext_queries = CreatePlaintextQueries();
  *encrypt_request.mutable_plaintext_queries() = {plaintext_queries.begin(),
                                                  plaintext_queries.end()};
  absl::StatusOr<EncryptQueriesResponse> encrypt_response =
      EncryptQueries(encrypt_request);
  if (!encrypt_response.ok()) {
    return encrypt_response.status();
  }
  end = absl::Now();

  std::cout << "Generate client request: "
            << absl::ToInt64Milliseconds(end - start) << " ms" << std::endl;
  std::cout << "Client request size: " << encrypt_response->ByteSizeLong()
            << " bytes" << std::endl;

  // Server encodes database.
  std::cout << "Encoding server database." << std::endl;
  start = absl::Now();
  EncodeDatabaseRequest encode_request;
  *encode_request.mutable_parameters() = parameters;
  std::vector<RawDatabaseShard> raw_database = CreateRawDatabase();
  *encode_request.mutable_shards() = {raw_database.begin(), raw_database.end()};
  absl::StatusOr<EncodeDatabaseResponse> encode_response =
      EncodeDatabase(encode_request);
  if (!encode_response.ok()) {
    return encode_response.status();
  }
  end = absl::Now();

  std::cout << "Encode server database: "
            << absl::ToInt64Milliseconds(end - start) << " ms" << std::endl;

  // Server processes request and produces response.
  std::cout << "Applying queries to server to produce response." << std::endl;
  start = absl::Now();
  ApplyQueriesRequest apply_request;
  *apply_request.mutable_parameters() = parameters;
  *apply_request.mutable_public_key() = public_key;
  *apply_request.mutable_encoded_database()->mutable_shards() =
      encode_response->shards();
  *apply_request.add_queries() = encrypt_response->encrypted_queries();
  apply_request.set_finalize_results(true);
  absl::StatusOr<ApplyQueriesResponse> apply_response =
      ApplyQueries(apply_request);
  if (!apply_response.ok()) {
    return apply_response.status();
  }
  end = absl::Now();

  std::cout << "Apply queries: " << absl::ToInt64Milliseconds(end - start)
            << " ms" << std::endl;
  std::cout << "Server response size: " << apply_response->ByteSizeLong()
            << " bytes" << std::endl;

  // Client processes response to obtain final result.
  std::cout << "Decrypting response." << std::endl;
  start = absl::Now();
  DecryptQueriesRequest decrypt_request;
  *decrypt_request.mutable_parameters() = parameters;
  *decrypt_request.mutable_public_key() = public_key;
  *decrypt_request.mutable_private_key() = private_key;
  *decrypt_request.mutable_encrypted_queries() =
      apply_response->query_results();
  absl::StatusOr<DecryptQueriesResponse> decrypt_response =
      DecryptQueries(decrypt_request);
  if (!decrypt_response.ok()) {
    return decrypt_response.status();
  }
  end = absl::Now();

  std::cout << "Decrypting response: " << absl::ToInt64Milliseconds(end - start)
            << " ms" << std::endl;

  // Check expected results vs. actual results.
  std::cout << "Checking results." << std::endl;

  if (decrypt_response->result_size() < plaintext_queries.size()) {
    return absl::InternalError("Encrypted responses are missing.");
  } else if (decrypt_response->result_size() > plaintext_queries.size()) {
    return absl::InternalError("Too many encrypted results returned.");
  }

  for (const PlaintextQuery& query : plaintext_queries) {
    int query_id = query.query_metadata().query_id();
    bool found = false;
    for (const DecryptedQueryResult& result : decrypt_response->result()) {
      if (query_id == result.query_metadata().query_id()) {
        found = true;
        std::string expected = CreateDatabaseContents(
            query.query_metadata().shard_id(), query.bucket_id());
        std::string actual = result.result();
        if (expected != actual) {
          return absl::InternalError(absl::StrCat(
              "Expected and actual results differed for query with ID: ",
              query_id, ". Expected: ", expected, " vs. Actual: ", actual));
        }
        break;
      }
    }
    if (!found) {
      return absl::InternalError(
          absl::StrCat("Missing result for query with ID: ", query_id));
    }
  }

  return absl::OkStatus();
}

}  // namespace
}  // namespace batch
}  // namespace private_membership

int main(int argc, char** argv) {
  absl::Status status = private_membership::batch::RunExample();
  if (!status.ok()) {
    std::cout << status << std::endl;
    return -1;
  } else {
    std::cout << "Simple example completed successfully!" << std::endl;
  }
  return 0;
}
