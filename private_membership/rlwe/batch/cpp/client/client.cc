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

#include "private_membership/rlwe/batch/cpp/client/client.h"

#include <cstdint>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "private_membership/rlwe/batch/proto/client.pb.h"
#include "private_membership/rlwe/batch/cpp/client/client_helper.h"
#include "private_membership/rlwe/batch/cpp/constants.h"
#include "private_membership/rlwe/batch/cpp/context.h"
#include "private_membership/rlwe/batch/cpp/padding.h"
#include "private_membership/rlwe/batch/cpp/prng.h"
#include "private_membership/rlwe/batch/proto/shared.pb.h"
#include "shell_encryption/context.h"
#include "shell_encryption/galois_key.h"
#include "shell_encryption/integral_types.h"
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

absl::StatusOr<rlwe::SymmetricRlweKey<ModularInt>> GeneratePrivateKey(
    const Context& context) {
  auto prng = CreatePrng();
  if (!prng.ok()) {
    return prng.status();
  }

  return rlwe::SymmetricRlweKey<ModularInt>::Sample(
      context.GetLogN(), context.GetVariance(), context.GetLogT(),
      context.GetModulusParams(), context.GetNttParams(), prng->get());
}

absl::StatusOr<rlwe::SymmetricRlweKey<ModularInt>> DeserializePrivateKey(
    const rlwe::SerializedNttPolynomial& serialization,
    const Context& context) {
  return rlwe::SymmetricRlweKey<ModularInt>::Deserialize(
      context.GetVariance(), context.GetLogT(), serialization,
      context.GetModulusParams(), context.GetNttParams());
}

absl::StatusOr<ModularInt> CreateNormalizer(int log_compression_factor,
                                            const Context& context) {
  // Create normalizer to normalize compressed plaintexts.
  rlwe::Uint64 normalizer =
      rlwe::ObliviousExpander<ModularInt>::ComputeNormalizer(
          log_compression_factor, context.GetLogT());
  return ModularInt::ImportInt(normalizer, context.GetModulusParams());
}

absl::StatusOr<std::vector<rlwe::Polynomial<ModularInt>>> CompressBitVector(
    const BitVector& bit_vector, const Parameters& parameters,
    const Context& context) {
  // Create compressed vector of plaintext values.
  int log_compression_factor =
      parameters.crypto_parameters().log_compression_factor();
  absl::StatusOr<std::vector<rlwe::Polynomial<ModularInt>>>
      compressed_plaintexts = rlwe::MakeCompressedVector(
          bit_vector.total_size, bit_vector.indices, log_compression_factor,
          context.GetModulusParams(), context.GetNttParams());
  if (!compressed_plaintexts.ok()) {
    return compressed_plaintexts.status();
  }

  // Create normalizer to normalize compressed plaintexts.
  absl::StatusOr<ModularInt> normalizer =
      CreateNormalizer(log_compression_factor, context);
  if (!normalizer.ok()) {
    return normalizer.status();
  }

  // Normalize compressed plaintexts.
  for (rlwe::Polynomial<ModularInt>& plaintext : *compressed_plaintexts) {
    absl::Status status =
        plaintext.MulInPlace(*normalizer, context.GetModulusParams());
    if (!status.ok()) {
      return status;
    }
  }
  return *compressed_plaintexts;
}

// Two different PRNGs used during encryption.
struct EncryptionPrng {
  // PRNG used for compression. The seed needs to be used to decompress.
  PrngWithSeed compression_prng;

  // PRNG used for randomness. The seed may be thrown out afterwards.
  std::unique_ptr<rlwe::SecurePrng> encryption_prng;
};

absl::StatusOr<EncryptionPrng> CreateEncryptionPrng() {
  EncryptionPrng encryption_prng;
  absl::StatusOr<PrngWithSeed> compression_prng_with_seed =
      CreatePrngWithSeed();
  if (!compression_prng_with_seed.ok()) {
    return compression_prng_with_seed.status();
  }
  auto prng = CreatePrng();
  if (!prng.ok()) {
    return prng.status();
  }

  encryption_prng.compression_prng = *std::move(compression_prng_with_seed);
  encryption_prng.encryption_prng = *std::move(prng);

  return encryption_prng;
}

absl::StatusOr<std::string> DecryptQueryResult(
    const EncryptedQueryResult& result,
    const rlwe::SymmetricRlweKey<ModularInt>& private_key,
    const Context& context) {
  std::vector<uint8_t> byte_result;
  for (const rlwe::SerializedSymmetricRlweCiphertext& serialized_ciphertext :
       result.ciphertexts()) {
    // Deserialize the ciphertext.
    auto ciphertext = rlwe::SymmetricRlweCiphertext<ModularInt>::Deserialize(
        serialized_ciphertext, context.GetModulusParams(),
        context.GetErrorParams());
    if (!ciphertext.ok()) {
      return ciphertext.status();
    }

    // Decrypt the ciphertext.
    absl::StatusOr<std::vector<ModularInt::Int>> decryption =
        rlwe::Decrypt(private_key, *ciphertext);
    if (!decryption.ok()) {
      return decryption.status();
    }

    // Transcribe the plaintext polynomial back to bytes.
    absl::StatusOr<std::vector<uint8_t>> bytes =
        rlwe::TranscribeBits<ModularInt::Int, uint8_t>(
            *std::move(decryption),
            private_key.Len() * private_key.BitsPerCoeff(),
            private_key.BitsPerCoeff(), 8);
    if (!bytes.ok()) {
      return bytes.status();
    }
    byte_result.insert(byte_result.end(),
                       std::make_move_iterator(bytes->begin()),
                       std::make_move_iterator(bytes->end()));
  }

  // Unpad the result.
  absl::string_view padded_result(
      reinterpret_cast<const char*>(byte_result.data()), byte_result.size());
  return Unpad(padded_result);
}

}  // namespace

absl::StatusOr<GenerateKeysResponse> GenerateKeys(
    const GenerateKeysRequest& request) {
  auto request_context = CreateRlweRequestContext(request.parameters());
  if (!request_context.ok()) {
    return request_context.status();
  }

  auto response_context = CreateRlweResponseContext(request.parameters());
  if (!response_context.ok()) {
    return response_context.status();
  }

  GenerateKeysResponse response;
  auto request_key = GeneratePrivateKey(**request_context);
  if (!request_key.ok()) {
    return request_key.status();
  }
  auto serialized_request_key = request_key->Serialize();
  if (!serialized_request_key.ok()) {
    return serialized_request_key.status();
  }
  *response.mutable_private_key()->mutable_request_key() =
      *std::move(serialized_request_key);

  auto response_key =
      request_key->SwitchModulus((*response_context)->GetModulusParams(),
                                 (*response_context)->GetNttParams());
  if (!response_key.ok()) {
    return response_key.status();
  }
  auto serialized_response_key = response_key->Serialize();
  if (!serialized_response_key.ok()) {
    return serialized_response_key.status();
  }
  *response.mutable_private_key()->mutable_response_key() =
      *std::move(serialized_response_key);

  const int log_compression_factor =
      request.parameters().crypto_parameters().log_compression_factor();
  const int log_decomposition_modulus =
      request.parameters().crypto_parameters().log_decomposition_modulus();
  const int request_key_length = request_key->Len();
  for (int i = 0; i < log_compression_factor; ++i) {
    int substitution_power = (request_key_length >> i) + 1;
    auto galois_key = rlwe::GaloisKey<ModularInt>::Create(
        *request_key, rlwe::PRNG_TYPE_HKDF, substitution_power,
        log_decomposition_modulus);
    if (!galois_key.ok()) {
      return galois_key.status();
    }
    auto serialized_galois_key = galois_key->Serialize();
    if (!serialized_galois_key.ok()) {
      return serialized_galois_key.status();
    }
    *response.mutable_public_key()->add_galois_key() =
        *std::move(serialized_galois_key);
  }

  return response;
}

absl::StatusOr<EncryptQueriesResponse> EncryptQueries(
    const EncryptQueriesRequest& request) {
  // Generate context from parameters.
  const Parameters& parameters = request.parameters();
  auto request_context = CreateRlweRequestContext(parameters);
  if (!request_context.ok()) {
    return request_context.status();
  }

  // Deserialize private key.
  auto private_key = DeserializePrivateKey(request.private_key().request_key(),
                                           **request_context);
  if (!private_key.ok()) {
    return private_key.status();
  }

  // Generate bit vector to be compressed and encrypted.
  absl::StatusOr<BitVector> bit_vector = GenerateBitVector(request);
  if (!bit_vector.ok()) {
    return bit_vector.status();
  }

  // Compress bit vector.
  auto compressed_bit_vector =
      CompressBitVector(*bit_vector, parameters, **request_context);
  if (!compressed_bit_vector.ok()) {
    return compressed_bit_vector.status();
  }

  // Generate necessary PRNGs.
  absl::StatusOr<EncryptionPrng> encryption_prng = CreateEncryptionPrng();
  if (!encryption_prng.ok()) {
    return encryption_prng.status();
  }

  // Encrypt compressed plaintexts.
  auto ciphertexts =
      rlwe::EncryptWithPrng(*private_key, *compressed_bit_vector,
                            encryption_prng->compression_prng.prng.get(),
                            encryption_prng->encryption_prng.get());
  if (!ciphertexts.ok()) {
    return ciphertexts.status();
  }

  // Serialize encryptions to response.
  EncryptQueriesResponse response;
  EncryptedQueries* encrypted_query = response.mutable_encrypted_queries();
  for (const auto& plaintext_query : request.plaintext_queries()) {
    *(encrypted_query->add_query_metadata()) = plaintext_query.query_metadata();
  }
  for (auto& ciphertext : *ciphertexts) {
    auto serialized_ciphertext =
        ciphertext.Serialize((*request_context)->GetModulusParams());
    if (!serialized_ciphertext.ok()) {
      return serialized_ciphertext.status();
    }
    (*encrypted_query->add_encrypted_queries()) =
        *std::move(serialized_ciphertext);
  }
  encrypted_query->set_prng_seed(encryption_prng->compression_prng.seed);
  return response;
}

absl::StatusOr<DecryptQueriesResponse> DecryptQueries(
    const DecryptQueriesRequest& request) {
  // Generate context from parameters.
  const Parameters& parameters = request.parameters();
  auto response_context = CreateRlweResponseContext(parameters);
  if (!response_context.ok()) {
    return response_context.status();
  }

  // Deserialize private key.
  auto private_key = DeserializePrivateKey(request.private_key().response_key(),
                                           **response_context);
  if (!private_key.ok()) {
    return private_key.status();
  }

  // Decrypt and transcribe plaintexts.
  DecryptQueriesResponse response;
  for (const EncryptedQueryResult& result : request.encrypted_queries()) {
    // Decrypt result.
    absl::StatusOr<std::string> plaintext_result =
        DecryptQueryResult(result, *private_key, **response_context);
    if (!plaintext_result.ok()) {
      return plaintext_result.status();
    }

    // Serialize result.
    DecryptedQueryResult* query_result = response.add_result();
    *query_result->mutable_query_metadata() = result.query_metadata();
    query_result->set_result(*std::move(plaintext_result));
  }

  return response;
}

}  // namespace batch
}  // namespace private_membership
