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

#include <cmath>
#include <iterator>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <google/protobuf/repeated_field.h>
#include "absl/algorithm/container.h"
#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/meta/type_traits.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "private_membership/rlwe/batch/cpp/constants.h"
#include "private_membership/rlwe/batch/cpp/context.h"
#include "private_membership/rlwe/batch/cpp/encoding.h"
#include "private_membership/rlwe/batch/cpp/prng.h"
#include "private_membership/rlwe/batch/proto/server.pb.h"
#include "private_membership/rlwe/batch/proto/shared.pb.h"
#include "shell_encryption/galois_key.h"
#include "shell_encryption/int256.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/oblivious_expand.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/serialization.pb.h"
#include "shell_encryption/symmetric_encryption.h"
#include "shell_encryption/symmetric_encryption_with_prng.h"
#include "shell_encryption/transcription.h"

namespace private_membership {
namespace batch {

namespace {

// Represents a query deserialized in RAM.
struct DeserializedQuery {
  int query_id = 0;
  int shard_id = 0;
  std::vector<rlwe::SymmetricRlweCiphertext<ModularInt>> ciphertexts;
};

// Represents a single bucket.
struct Bucket {
  int id = 0;
  std::vector<rlwe::Polynomial<ModularInt>> contents;
};

// Represents a single shard.
struct Shard {
  int width = 0;
  std::vector<Bucket> buckets;
};

// Context keeping track of necessary information for the ApplyQueries function.
// This context is only applicable for a single column and must be refreshed for
// each column.
class ApplyQueriesContext {
 public:
  // Create context from parameters.
  static absl::StatusOr<ApplyQueriesContext> Create(
      const Parameters& parameters) {
    int levels_of_recursion =
        parameters.crypto_parameters().levels_of_recursion();
    if (levels_of_recursion <= 0) {
      return absl::InvalidArgumentError(
          "Invalid parameters: levels_of_recursion must be positive.");
    }
    int buckets_per_shard =
        parameters.shard_parameters().number_of_buckets_per_shard();
    if (buckets_per_shard <= 0) {
      return absl::InvalidArgumentError(
          "Invalid parameters: buckets_per_shard must be positive.");
    }
    int queries_per_level =
        std::ceil(std::pow(buckets_per_shard, 1.0 / levels_of_recursion));
    if (queries_per_level <= 0) {
      return absl::InvalidArgumentError(
          "Invalid parameters: queries_per_level must be positive.");
    }
    return ApplyQueriesContext(
        levels_of_recursion, buckets_per_shard, queries_per_level,
        (buckets_per_shard + queries_per_level - 1) / queries_per_level);
  }

  // Compute information for next level of applying queries.
  void AdvanceToNextLevel() {
    next_database_size_ =
        (next_database_size_ + queries_per_level_ - 1) / queries_per_level_;
    next_query_range_start_ += queries_per_level_;
  }

  // Check that queries and database size are compatible at this level.
  absl::Status VerifyState(const DeserializedQuery& deserialized_query,
                           int current_database_size) const {
    if (next_query_range_start_ + queries_per_level_ - 1 >=
        deserialized_query.ciphertexts.size()) {
      return absl::InvalidArgumentError("Invalid parameters: too few queries.");
    }
    if ((current_database_size - 1) / queries_per_level_ >=
        next_database_size_) {
      return absl::InvalidArgumentError(
          "Invalid parameters: too many buckets for given queries.");
    }
    return absl::OkStatus();
  }

  int levels_of_recursion() const { return levels_of_recursion_; }
  int buckets_per_shard() const { return buckets_per_shard_; }
  int queries_per_level() const { return queries_per_level_; }
  int next_database_size() const { return next_database_size_; }
  int next_query_range_start() const { return next_query_range_start_; }

 private:
  ApplyQueriesContext(int levels_of_recursion, int buckets_per_shard,
                      int queries_per_level, int next_database_size)
      : levels_of_recursion_(levels_of_recursion),
        buckets_per_shard_(buckets_per_shard),
        queries_per_level_(queries_per_level),
        next_database_size_(next_database_size),
        next_query_range_start_(0) {}

  const int levels_of_recursion_;
  const int buckets_per_shard_;
  const int queries_per_level_;
  int next_database_size_;
  int next_query_range_start_;
};

absl::StatusOr<std::vector<rlwe::Polynomial<ModularInt>>>
DeserializePolynomials(
    const Context& context,
    const ::google::protobuf::RepeatedPtrField<rlwe::SerializedNttPolynomial>&
        serialized_polynomials) {
  std::vector<rlwe::Polynomial<ModularInt>> decoded;
  decoded.reserve(serialized_polynomials.size());
  for (const rlwe::SerializedNttPolynomial& serialized :
       serialized_polynomials) {
    absl::StatusOr<rlwe::Polynomial<ModularInt>> deserialized =
        rlwe::Polynomial<ModularInt>::Deserialize(serialized,
                                                  context.GetModulusParams());
    if (!deserialized.ok()) {
      return deserialized.status();
    }
    decoded.push_back(*std::move(deserialized));
  }
  return decoded;
}

absl::StatusOr<Bucket> DeserializeRawBucket(
    const RawDatabaseShard::Bucket& raw_bucket, const Context& context) {
  // Encode bucket contents to polynomials.
  absl::StatusOr<std::vector<rlwe::Polynomial<ModularInt>>>
      encoded_polynomials = EncodeBytes(context, raw_bucket.bucket_contents());
  if (!encoded_polynomials.ok()) {
    return encoded_polynomials.status();
  }

  Bucket bucket = {
      .id = raw_bucket.bucket_id(),
      .contents = *std::move(encoded_polynomials),
  };
  return bucket;
}

absl::StatusOr<EncodedDatabaseShard::Bucket> EncodeRawBucket(
    const RawDatabaseShard::Bucket& raw_bucket, const Context& context) {
  // Encode bucket contents to polynomials.
  absl::StatusOr<std::vector<rlwe::Polynomial<ModularInt>>>
      encoded_polynomials = EncodeBytes(context, raw_bucket.bucket_contents());
  if (!encoded_polynomials.ok()) {
    return encoded_polynomials.status();
  }

  EncodedDatabaseShard::Bucket encoded_bucket;
  encoded_bucket.set_bucket_id(raw_bucket.bucket_id());

  // Serialize polynomials in response.
  for (const rlwe::Polynomial<ModularInt>& encoded_polynomial :
       *encoded_polynomials) {
    absl::StatusOr<rlwe::SerializedNttPolynomial> serialized_polynomial =
        encoded_polynomial.Serialize(context.GetModulusParams());
    if (!serialized_polynomial.ok()) {
      return serialized_polynomial.status();
    }
    *encoded_bucket.add_chunks() = *std::move(serialized_polynomial);
  }

  return encoded_bucket;
}

absl::StatusOr<rlwe::SerializedSymmetricRlweCiphertext> SumCiphertextList(
    const Context& context,
    const ::google::protobuf::RepeatedPtrField<rlwe::SerializedSymmetricRlweCiphertext>&
        list) {
  if (list.empty()) {
    return absl::InvalidArgumentError(
        "Each summation must contain at least one ciphertext.");
  }

  absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> sum =
      rlwe::SymmetricRlweCiphertext<ModularInt>::Deserialize(
          list.Get(0), context.GetModulusParams(), context.GetErrorParams());
  if (!sum.ok()) {
    return sum.status();
  }

  for (int i = 1; i < list.size(); ++i) {
    absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> addend =
        rlwe::SymmetricRlweCiphertext<ModularInt>::Deserialize(
            list.Get(i), context.GetModulusParams(), context.GetErrorParams());
    if (!addend.ok()) {
      return addend.status();
    }
    absl::Status add_status = sum->AddInPlace(*addend);
    if (!add_status.ok()) {
      return add_status;
    }
  }

  return sum->Serialize();
}

// Finalize results when given deserialized ciphertexts. The input result is
// passed by value because modulus switching requires the object to be
// non-const.
absl::StatusOr<rlwe::SerializedSymmetricRlweCiphertext> FinalizeResult(
    const Context& request_context, const Context& response_context,
    rlwe::SymmetricRlweCiphertext<ModularInt> result) {
  absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> finalized =
      result.SwitchModulus(
          request_context.GetNttParams(), response_context.GetModulusParams(),
          response_context.GetNttParams(), response_context.GetErrorParams(),
          response_context.GetT());
  if (!finalized.ok()) {
    return finalized.status();
  }
  return finalized->Serialize();
}

// Finalize results when given serialized ciphertexts.
absl::StatusOr<rlwe::SerializedSymmetricRlweCiphertext> FinalizeResult(
    const Context& request_context, const Context& response_context,
    const rlwe::SerializedSymmetricRlweCiphertext& result) {
  absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> deserialized =
      rlwe::SymmetricRlweCiphertext<ModularInt>::Deserialize(
          result, request_context.GetModulusParams(),
          request_context.GetErrorParams());
  if (!deserialized.ok()) {
    return deserialized.status();
  }
  return FinalizeResult(request_context, response_context,
                        *std::move(deserialized));
}

absl::StatusOr<absl::flat_hash_map<int, std::vector<DeserializedQuery>>>
SeparateQueries(
    std::vector<rlwe::SymmetricRlweCiphertext<ModularInt>> query_list,
    const ::google::protobuf::RepeatedPtrField<QueryMetadata>& query_metadata_list,
    int ciphertexts_per_query) {
  if (query_list.size() != query_metadata_list.size() * ciphertexts_per_query) {
    return absl::InvalidArgumentError(
        "Invalid number of ciphertexts provided.");
  }
  absl::flat_hash_map<int, std::vector<DeserializedQuery>> queries;
  auto it = query_list.begin();
  for (const QueryMetadata& query_metadata : query_metadata_list) {
    std::vector<rlwe::SymmetricRlweCiphertext<ModularInt>> query_ciphertexts;
    query_ciphertexts.reserve(ciphertexts_per_query);
    std::move(it, it + ciphertexts_per_query,
              std::back_inserter(query_ciphertexts));
    int shard_id = query_metadata.shard_id();
    queries[shard_id].push_back({
        .query_id = query_metadata.query_id(),
        .shard_id = shard_id,
        .ciphertexts = std::move(query_ciphertexts),
    });
    it += ciphertexts_per_query;
  }
  return queries;
}

// Create the Expander object used to expand ciphertexts from the serialized
// public Galois keys.
absl::StatusOr<std::unique_ptr<rlwe::GaloisKeysObliviousExpander<ModularInt>>>
CreateExpander(
    const Context& context,
    const ::google::protobuf::RepeatedPtrField<rlwe::SerializedGaloisKey>& public_keys) {
  const int levels_of_expansion = public_keys.size();
  std::vector<rlwe::GaloisKey<ModularInt>> public_galois_keys;
  public_galois_keys.reserve(levels_of_expansion);
  for (const rlwe::SerializedGaloisKey& serialized_key : public_keys) {
    absl::StatusOr<rlwe::GaloisKey<ModularInt>> public_galois_key =
        DeserializePublicGaloisKey(serialized_key, context);
    if (!public_galois_key.ok()) {
      return public_galois_key.status();
    }
    public_galois_keys.push_back(*std::move(public_galois_key));
  }
  return rlwe::GaloisKeysObliviousExpander<ModularInt>::Create(
      std::move(public_galois_keys), context.GetLogT(),
      context.GetModulusParams(), context.GetNttParams());
}

absl::StatusOr<absl::flat_hash_map<int, Shard>> DeserializeDatabase(
    const Context& context,
    const ::google::protobuf::RepeatedPtrField<EncodedDatabaseShard>& encoded_shards) {
  absl::flat_hash_map<int, Shard> database;
  database.reserve(encoded_shards.size());

  // Keep track of seen (shard_id, bucket_id) pairs.
  absl::flat_hash_set<std::pair<int, int>> seen_buckets;
  seen_buckets.reserve(encoded_shards.size());

  for (const EncodedDatabaseShard& shard : encoded_shards) {
    int shard_index = shard.shard_index();
    Shard& deserialized_shard = database[shard_index];
    for (const EncodedDatabaseShard::Bucket& bucket : shard.buckets()) {
      int bucket_id = bucket.bucket_id();
      absl::StatusOr<std::vector<rlwe::Polynomial<ModularInt>>> entire_bucket =
          DeserializePolynomials(context, bucket.chunks());
      if (!entire_bucket.ok()) {
        return entire_bucket.status();
      }
      auto [unused, inserted] = seen_buckets.insert({shard_index, bucket_id});
      if (!inserted) {
        return absl::InvalidArgumentError("Database contains bucket twice.");
      }
      if (entire_bucket->size() > deserialized_shard.width) {
        deserialized_shard.width = entire_bucket->size();
      }
      deserialized_shard.buckets.push_back(
          {.id = bucket_id, .contents = *std::move(entire_bucket)});
    }
  }
  return database;
}

absl::StatusOr<absl::flat_hash_map<int, Shard>> DeserializeDatabase(
    const Context& context,
    const ::google::protobuf::RepeatedPtrField<RawDatabaseShard>& raw_shards) {
  absl::flat_hash_map<int, Shard> database;
  database.reserve(raw_shards.size());

  // Keep track of seen (shard_id, bucket_id) pairs.
  absl::flat_hash_set<std::pair<int, int>> seen_buckets;
  seen_buckets.reserve(raw_shards.size());

  for (const RawDatabaseShard& shard : raw_shards) {
    int shard_index = shard.shard_index();
    Shard& deserialized_shard = database[shard_index];
    for (const RawDatabaseShard::Bucket& bucket : shard.buckets()) {
      int bucket_id = bucket.bucket_id();
      absl::StatusOr<Bucket> entire_bucket =
          DeserializeRawBucket(bucket, context);
      if (!entire_bucket.ok()) {
        return entire_bucket.status();
      }
      auto [unused, inserted] = seen_buckets.insert({shard_index, bucket_id});
      if (!inserted) {
        return absl::InvalidArgumentError("Database contains bucket twice.");
      }
      if (entire_bucket->contents.size() > deserialized_shard.width) {
        deserialized_shard.width = entire_bucket->contents.size();
      }
      deserialized_shard.buckets.push_back(*std::move(entire_bucket));
    }
  }
  return database;
}

absl::StatusOr<absl::flat_hash_map<int, std::vector<DeserializedQuery>>>
DeserializeQueries(const Context& context, const Parameters& parameters,
                   EncryptedQueries encrypted_queries, int levels_of_expansion,
                   rlwe::GaloisKeysObliviousExpander<ModularInt>* expander) {
  const int levels_of_recursion =
      parameters.crypto_parameters().levels_of_recursion();
  if (levels_of_recursion <= 0) {
    return absl::InvalidArgumentError(
        "Invalid parameters. Levels of recursion must be positive.");
  }
  const int buckets_per_query =
      parameters.shard_parameters().number_of_buckets_per_shard();
  const int ciphertexts_per_level = static_cast<int>(
      std::ceil(std::pow(buckets_per_query, 1.0 / levels_of_recursion)));
  const int ciphertexts_per_query = ciphertexts_per_level * levels_of_recursion;

  // Deserialize to compressed queries.
  absl::StatusOr<std::vector<rlwe::Polynomial<ModularInt>>> compressed_queries =
      DeserializePolynomials(context, encrypted_queries.encrypted_queries());
  if (!compressed_queries.ok()) {
    return compressed_queries.status();
  }

  // Deserialize PRNG.
  absl::StatusOr<std::unique_ptr<rlwe::SecurePrng>> prng =
      CreatePrngFromSeed(encrypted_queries.prng_seed());
  if (!prng.ok()) {
    return prng.status();
  }

  // Expand queries using PRNG.
  absl::StatusOr<std::vector<rlwe::SymmetricRlweCiphertext<ModularInt>>>
      expanded_queries = rlwe::ExpandFromPrng(
          *std::move(compressed_queries), context.GetModulusParams(),
          context.GetNttParams(), context.GetErrorParams(), prng->get());
  if (!expanded_queries.ok()) {
    return expanded_queries.status();
  }

  // Expand queries using public key.
  const int number_of_ciphertexts =
      ciphertexts_per_query * encrypted_queries.query_metadata_size();
  absl::StatusOr<std::vector<rlwe::SymmetricRlweCiphertext<ModularInt>>>
      final_queries =
          expander->ObliviousExpand(*std::move(expanded_queries),
                                    levels_of_expansion, number_of_ciphertexts);
  if (!final_queries.ok()) {
    return final_queries.status();
  }

  // Map expanded ciphertexts back to individual queries.
  return SeparateQueries(*std::move(final_queries),
                         encrypted_queries.query_metadata(),
                         ciphertexts_per_query);
}

// Applies first level of queries to database. The first level consists of the
// original encoded database being applied to the deserialized queries.
absl::StatusOr<
    std::vector<absl::optional<rlwe::SymmetricRlweCiphertext<ModularInt>>>>
ApplyFirstLevel(const Context& context,
                const ApplyQueriesContext& apply_context, const Shard& shard,
                const DeserializedQuery& deserialized_query, int column_index) {
  std::vector<absl::optional<rlwe::SymmetricRlweCiphertext<ModularInt>>> result(
      apply_context.next_database_size());

  // Verify the current state of apply context is valid for this level.
  absl::Status verify_status = apply_context.VerifyState(
      deserialized_query, apply_context.buckets_per_shard());
  if (!verify_status.ok()) {
    return verify_status;
  }

  for (const Bucket& bucket : shard.buckets) {
    // Skip buckets without any columns left.
    if (bucket.contents.size() <= column_index) {
      continue;
    }
    int query_index = bucket.id % apply_context.queries_per_level();
    int result_index = bucket.id / apply_context.queries_per_level();

    absl::optional<rlwe::SymmetricRlweCiphertext<ModularInt>>& ciphertext =
        result[result_index];
    if (!ciphertext.has_value()) {
      ciphertext = rlwe::SymmetricRlweCiphertext<ModularInt>(
          context.GetModulusParams(), context.GetErrorParams());
    }
    absl::Status status = ciphertext->FusedAbsorbAddInPlaceLazily(
        deserialized_query.ciphertexts[query_index],
        bucket.contents[column_index]);
    if (!status.ok()) {
      return status;
    }
  }

  for (int i = 0; i < result.size(); ++i) {
    if (result[i].has_value()) {
      absl::Status status = result[i]->MergeLazyOperations();
      if (!status.ok()) {
        return status;
      }
    }
  }

  return result;
}

// Apply the next level of queries to the intermediate database. Unlike the
// first level, this applies queries to ciphertexts from prior levels.
absl::StatusOr<
    std::vector<absl::optional<rlwe::SymmetricRlweCiphertext<ModularInt>>>>
ApplyNextLevel(
    const Context& context, const ApplyQueriesContext& apply_context,
    const DeserializedQuery& deserialized_query,
    std::vector<absl::optional<rlwe::SymmetricRlweCiphertext<ModularInt>>>&
        intermediate_database) {
  std::vector<absl::optional<rlwe::SymmetricRlweCiphertext<ModularInt>>> result(
      apply_context.next_database_size());

  // Verify the current state of apply context is valid for this level.
  absl::Status verify_status = apply_context.VerifyState(
      deserialized_query, intermediate_database.size());
  if (!verify_status.ok()) {
    return verify_status;
  }

  for (int i = 0; i < intermediate_database.size(); ++i) {
    absl::optional<rlwe::SymmetricRlweCiphertext<ModularInt>>& database_entry =
        intermediate_database[i];
    if (!database_entry.has_value()) {
      continue;
    }

    int query_index = apply_context.next_query_range_start() +
                      (i % apply_context.queries_per_level());
    int result_index = i / apply_context.queries_per_level();

    absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> product =
        database_entry.value() * deserialized_query.ciphertexts[query_index];
    if (!product.ok()) {
      return product.status();
    }

    absl::optional<rlwe::SymmetricRlweCiphertext<ModularInt>>& ciphertext =
        result[result_index];
    if (!ciphertext.has_value()) {
      ciphertext = *std::move(product);
    } else {
      absl::Status add_status = ciphertext->AddInPlace(*std::move(product));
      if (!add_status.ok()) {
        return add_status;
      }
    }
  }

  return result;
}

// Apply queries to a single column of the database.
absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> ApplyQueryToColumn(
    const Context& context, const Parameters& parameters, const Shard& shard,
    const DeserializedQuery& deserialized_query, int column_index) {
  // Store all necessary information for applying queries.
  absl::StatusOr<ApplyQueriesContext> apply_context =
      ApplyQueriesContext::Create(parameters);
  if (!apply_context.ok()) {
    return apply_context.status();
  }

  absl::StatusOr<
      std::vector<absl::optional<rlwe::SymmetricRlweCiphertext<ModularInt>>>>
      result = ApplyFirstLevel(context, *apply_context, shard,
                               deserialized_query, column_index);
  if (!result.ok()) {
    return result.status();
  }
  apply_context->AdvanceToNextLevel();

  for (int level = 1; level < apply_context->levels_of_recursion(); ++level) {
    absl::StatusOr<
        std::vector<absl::optional<rlwe::SymmetricRlweCiphertext<ModularInt>>>>
        intermediate_result = ApplyNextLevel(context, *apply_context,
                                             deserialized_query, *result);
    if (!intermediate_result.ok()) {
      return intermediate_result.status();
    }
    result = *std::move(intermediate_result);
    apply_context->AdvanceToNextLevel();
  }

  if (result->size() != 1 || !(*result)[0].has_value()) {
    return absl::InvalidArgumentError(
        "Parameters are invalid to apply queries correctly.");
  }

  return std::move((*result)[0].value());
}

absl::StatusOr<std::vector<rlwe::SymmetricRlweCiphertext<ModularInt>>>
ApplyQuery(const Context& context, const Parameters& parameters,
           const Shard& shard, const DeserializedQuery& deserialized_query) {
  std::vector<rlwe::SymmetricRlweCiphertext<ModularInt>> results;
  results.reserve(shard.width);
  for (int i = 0; i < shard.width; ++i) {
    absl::StatusOr<rlwe::SymmetricRlweCiphertext<ModularInt>> result =
        ApplyQueryToColumn(context, parameters, shard, deserialized_query,
                           /*column_index=*/i);
    if (!result.ok()) {
      return result.status();
    }
    results.push_back(*std::move(result));
  }
  return results;
}

absl::StatusOr<std::vector<rlwe::SerializedSymmetricRlweCiphertext>>
SerializeOrFinalizeQueries(
    const Context& context, const Context& finalize_context,
    const std::vector<rlwe::SymmetricRlweCiphertext<ModularInt>>& ciphertexts,
    bool finalize_results) {
  std::vector<rlwe::SerializedSymmetricRlweCiphertext> serialized_result;
  serialized_result.reserve(ciphertexts.size());
  for (int i = 0; i < ciphertexts.size(); ++i) {
    absl::StatusOr<rlwe::SerializedSymmetricRlweCiphertext> serialized;
    if (finalize_results) {
      serialized =
          FinalizeResult(context, finalize_context, std::move(ciphertexts[i]));
    } else {
      serialized = ciphertexts[i].Serialize();
    }
    if (!serialized.ok()) {
      return serialized.status();
    }
    serialized_result.push_back(*std::move(serialized));
  }
  return serialized_result;
}

absl::Status ApplyQueryToShard(
    const Context& context, const Context& finalize_context,
    const Parameters& parameters,
    const std::vector<DeserializedQuery>& shard_queries, const Shard& shard,
    bool finalize_results, ApplyQueriesResponse& response) {
  for (const DeserializedQuery& shard_query : shard_queries) {
    // Apply query to the relevant shard.
    absl::StatusOr<std::vector<rlwe::SymmetricRlweCiphertext<ModularInt>>>
        result_ciphertexts =
            ApplyQuery(context, parameters, shard, shard_query);
    if (!result_ciphertexts.ok()) {
      return result_ciphertexts.status();
    }

    // Serialize (and potentially finalize) resulting ciphertexts.
    absl::StatusOr<std::vector<rlwe::SerializedSymmetricRlweCiphertext>>
        serialized_result = SerializeOrFinalizeQueries(
            context, finalize_context, *result_ciphertexts, finalize_results);
    if (!serialized_result.ok()) {
      return serialized_result.status();
    }

    // Populate metadata and move ciphertexts into response.
    EncryptedQueryResult* query_result = response.add_query_results();

    QueryMetadata* metadata = query_result->mutable_query_metadata();
    metadata->set_shard_id(shard_query.shard_id);
    metadata->set_query_id(shard_query.query_id);

    absl::c_move(*serialized_result, google::protobuf::RepeatedPtrFieldBackInserter(
                                         query_result->mutable_ciphertexts()));
  }
  return absl::OkStatus();
}

}  // namespace

absl::StatusOr<EncodeDatabaseResponse> EncodeDatabase(
    const EncodeDatabaseRequest& request) {
  auto request_context = CreateRlweRequestContext(request.parameters());
  if (!request_context.ok()) {
    return request_context.status();
  }

  EncodeDatabaseResponse response;

  for (const RawDatabaseShard& raw_shard : request.shards()) {
    if (raw_shard.buckets_size() == 0) {
      return absl::InvalidArgumentError("Raw shard should not be empty.");
    }
    EncodedDatabaseShard* encoded_shard = response.add_shards();
    encoded_shard->set_shard_index(raw_shard.shard_index());
    for (const RawDatabaseShard::Bucket& raw_bucket : raw_shard.buckets()) {
      absl::StatusOr<EncodedDatabaseShard::Bucket> encoded_bucket =
          EncodeRawBucket(raw_bucket, **request_context);
      if (!encoded_bucket.ok()) {
        return encoded_bucket.status();
      }
      *encoded_shard->add_buckets() = *std::move(encoded_bucket);
    }
  }

  return response;
}

absl::StatusOr<ApplyQueriesResponse> ApplyQueries(
    const ApplyQueriesRequest& request) {
  if (request.public_key().galois_key().empty()) {
    return absl::InvalidArgumentError("Expect at least one public key");
  } else if (request.has_encoded_database() &&
             request.encoded_database().shards().empty()) {
    return absl::InvalidArgumentError("Expect at least one encoded shard");
  } else if (request.has_raw_database() &&
             request.raw_database().shards().empty()) {
    return absl::InvalidArgumentError("Expect at least one encoded shard");
  } else if (request.queries().empty()) {
    return absl::InvalidArgumentError("Expect at least one query");
  }

  // Compute necessary information from parameters.
  auto context = CreateRlweRequestContext(request.parameters());
  if (!context.ok()) {
    return context.status();
  }
  // Compute necessary information from parameters.
  auto finalize_context = CreateRlweResponseContext(request.parameters());
  if (!finalize_context.ok()) {
    return finalize_context.status();
  }

  // Deserialize public keys.
  const int levels_of_expansion = request.public_key().galois_key_size();
  absl::StatusOr<std::unique_ptr<rlwe::GaloisKeysObliviousExpander<ModularInt>>>
      expander = CreateExpander(**context, request.public_key().galois_key());
  if (!expander.ok()) {
    return expander.status();
  }

  // Deserialize database shards.
  absl::StatusOr<absl::flat_hash_map<int, Shard>> database;
  if (request.has_encoded_database()) {
    database =
        DeserializeDatabase(**context, request.encoded_database().shards());
  } else if (request.has_raw_database()) {
    database = DeserializeDatabase(**context, request.raw_database().shards());
  } else {
    return absl::InvalidArgumentError(
        "One of raw or encoded databases must be specified.");
  }
  if (!database.ok()) {
    return database.status();
  }

  ApplyQueriesResponse response;
  for (const EncryptedQueries& encrypted_queries : request.queries()) {
    // Deserialize all queries.
    absl::StatusOr<absl::flat_hash_map<int, std::vector<DeserializedQuery>>>
        queries = DeserializeQueries(**context, request.parameters(),
                                     encrypted_queries, levels_of_expansion,
                                     expander->get());
    if (!queries.ok()) {
      return queries.status();
    }

    for (const auto& [shard_id, shard_queries] : *queries) {
      // Skip queries without any corresponding database.
      auto database_shard = database->find(shard_id);
      if (database_shard == database->end()) {
        continue;
      }

      absl::Status apply_status = ApplyQueryToShard(
          **context, **finalize_context, request.parameters(), shard_queries,
          database_shard->second, request.finalize_results(), response);
      if (!apply_status.ok()) {
        return apply_status;
      }
    }
  }

  return response;
}

absl::StatusOr<SumCiphertextsResponse> SumCiphertexts(
    const SumCiphertextsRequest& request) {
  auto context = CreateRlweRequestContext(request.parameters());
  if (!context.ok()) {
    return context.status();
  }

  SumCiphertextsResponse response;
  for (const SumCiphertextsRequest::Summation& summation :
       request.summations()) {
    if (summation.ciphertexts_size() == 1) {
      // If only one ciphertext appears in a sub-list, return it without
      // deserializing.
      *response.add_ciphertexts() = summation.ciphertexts(0);
    } else {
      absl::StatusOr<rlwe::SerializedSymmetricRlweCiphertext> sum =
          SumCiphertextList(**context, summation.ciphertexts());
      if (!sum.ok()) {
        return sum.status();
      }
      *response.add_ciphertexts() = *std::move(sum);
    }
  }

  return response;
}

absl::StatusOr<FinalizeResultsResponse> FinalizeResults(
    const FinalizeResultsRequest& request) {
  auto request_context = CreateRlweRequestContext(request.parameters());
  if (!request_context.ok()) {
    return request_context.status();
  }
  auto response_context = CreateRlweResponseContext(request.parameters());
  if (!response_context.ok()) {
    return response_context.status();
  }

  FinalizeResultsResponse response;
  for (const EncryptedQueryResult& result : request.encrypted_results()) {
    EncryptedQueryResult* finalized_result = response.add_finalized_results();
    *finalized_result->mutable_query_metadata() = result.query_metadata();
    if (result.ciphertexts_size() <= 0) {
      return absl::InvalidArgumentError(
          "Each result should encode a non-empty string.");
    }
    for (const rlwe::SerializedSymmetricRlweCiphertext& ciphertext :
         result.ciphertexts()) {
      absl::StatusOr<rlwe::SerializedSymmetricRlweCiphertext>
          finalized_ciphertext =
              FinalizeResult(**request_context, **response_context, ciphertext);
      if (!finalized_ciphertext.ok()) {
        return finalized_ciphertext.status();
      }
      *finalized_result->add_ciphertexts() = *std::move(finalized_ciphertext);
    }
  }
  return response;
}

}  // namespace batch
}  // namespace private_membership
