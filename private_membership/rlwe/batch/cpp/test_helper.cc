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

#include "private_membership/rlwe/batch/cpp/test_helper.h"

#include <memory>
#include <utility>
#include <vector>

#include <google/protobuf/repeated_field.h>
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "private_membership/rlwe/batch/cpp/client/client.h"
#include "private_membership/rlwe/batch/cpp/constants.h"
#include "private_membership/rlwe/batch/cpp/prng.h"
#include "private_membership/rlwe/batch/proto/client.pb.h"
#include "private_membership/rlwe/batch/proto/shared.pb.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/serialization.pb.h"
#include "shell_encryption/symmetric_encryption.h"

namespace private_membership {
namespace batch {

Parameters CreateTestParameters() {
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
  crypto_parameters->set_levels_of_recursion(kLevelsOfRecursion);
  crypto_parameters->set_log_compression_factor(4);
  crypto_parameters->set_log_decomposition_modulus(10);

  return parameters;
}

absl::StatusOr<GenerateKeysResponse> CreateTestKeys() {
  GenerateKeysRequest request;
  *(request.mutable_parameters()) = CreateTestParameters();
  return GenerateKeys(request);
}

PlaintextQuery CreatePlaintextQuery(int query_id, int shard_id,
                                    int bucket_index) {
  PlaintextQuery query;
  query.mutable_query_metadata()->set_query_id(query_id);
  query.mutable_query_metadata()->set_shard_id(shard_id);
  query.set_bucket_id(bucket_index);
  return query;
}

absl::StatusOr<EncryptedQueries> CreateEncryptedQueries(
    const GenerateKeysResponse& keys,
    const std::vector<PlaintextQuery>& plaintext_queries) {
  EncryptQueriesRequest request;
  *request.mutable_parameters() = CreateTestParameters();
  *request.mutable_private_key() = keys.private_key();
  *request.mutable_public_key() = keys.public_key();
  *request.mutable_plaintext_queries() = {plaintext_queries.begin(),
                                          plaintext_queries.end()};

  absl::StatusOr<EncryptQueriesResponse> encrypted_queries =
      EncryptQueries(request);
  if (!encrypted_queries.ok()) {
    return encrypted_queries.status();
  }
  return encrypted_queries->encrypted_queries();
}

absl::StatusOr<rlwe::SymmetricRlweKey<ModularInt>> DeserializePrivateKey(
    const rlwe::SerializedNttPolynomial& serialized_key,
    const Context& context) {
  return rlwe::SymmetricRlweKey<ModularInt>::Deserialize(
      context.GetVariance(), context.GetLogT(), serialized_key,
      context.GetModulusParams(), context.GetNttParams());
}

absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> EncryptSingleInt(
    const Context& context,
    const rlwe::SymmetricRlweKey<ModularInt>& private_key, Int plaintext) {
  std::vector<ModularInt> coeffs;
  coeffs.reserve(context.GetN());
  for (int i = 0; i < context.GetN(); ++i) {
    coeffs.push_back(ModularInt::ImportZero(context.GetModulusParams()));
  }
  absl::StatusOr<ModularInt> imported_int =
      ModularInt::ImportInt(plaintext, context.GetModulusParams());
  if (!imported_int.ok()) {
    return imported_int.status();
  }
  coeffs[0] = *std::move(imported_int);
  absl::StatusOr<rlwe::Polynomial<ModularInt>> encoded =
      rlwe::Polynomial<ModularInt>::ConvertToNtt(coeffs, context.GetNttParams(),
                                                 context.GetModulusParams());
  if (!encoded.ok()) {
    return encoded.status();
  }
  auto prng = CreatePrng();
  if (!prng.ok()) {
    return prng.status();
  }
  return rlwe::Encrypt(private_key, *encoded, context.GetErrorParams(),
                       prng->get());
}

absl::StatusOr<Int> DecryptSingleInt(
    const Context& context,
    const rlwe::SymmetricRlweKey<ModularInt>& private_key,
    const rlwe::SymmetricRlweCiphertext<ModularInt>& ciphertext) {
  absl::StatusOr<std::vector<typename ModularInt::Int>> decryption =
      rlwe::Decrypt(private_key, ciphertext);
  if (!decryption.ok()) {
    return decryption.status();
  }

  if (decryption->empty()) {
    return absl::InternalError("Decryption should not be empty.");
  }

  return decryption->at(0);
}

}  // namespace batch
}  // namespace private_membership
