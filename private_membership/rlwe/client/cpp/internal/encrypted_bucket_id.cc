// Copyright 2020 Google LLC
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

#include "private_membership/rlwe/client/cpp/internal/encrypted_bucket_id.h"

#include "private_membership/rlwe/client/cpp/internal/crypto_utils.h"
#include "private_membership/rlwe/client/cpp/internal/rlwe_id_utils.h"
#include "private_membership/rlwe/client/cpp/internal/utils.h"
#include "shell_encryption/status_macros.h"

namespace private_membership {
namespace rlwe {

::rlwe::StatusOr<EncryptedBucketId> EncryptedBucketId::Create(
    const std::string& encrypted_bucket_id, int bit_length) {
  if (!IsValid(encrypted_bucket_id, bit_length)) {
    return absl::InvalidArgumentError("Invalid bit_length.");
  }
  return EncryptedBucketId(encrypted_bucket_id, bit_length);
}

::rlwe::StatusOr<EncryptedBucketId> EncryptedBucketId::Create(
    const RlwePlaintextId& id, const EncryptedBucketsParameters& params,
    ::private_join_and_compute::ECCommutativeCipher* ec_cipher, ::private_join_and_compute::Context* ctx) {
  if (ec_cipher == nullptr || ctx == nullptr) {
    return absl::InvalidArgumentError(
        "ECCipher and Context must both be non-null.");
  }
  std::string full_id = HashRlwePlaintextId(id);
  RLWE_ASSIGN_OR_RETURN(std::string encrypted_id, ec_cipher->Encrypt(full_id));
  return EncryptedBucketId::Create(encrypted_id, params, ctx);
}

::rlwe::StatusOr<EncryptedBucketId> EncryptedBucketId::Create(
    const std::string& encrypted_id, const EncryptedBucketsParameters& params,
    ::private_join_and_compute::Context* ctx) {
  if (ctx == nullptr) {
    return absl::InvalidArgumentError("Context must be non-null.");
  }
  std::string hashed_encrypted_id = HashEncryptedId(encrypted_id, ctx);
  RLWE_ASSIGN_OR_RETURN(
      std::string encrypted_bucket_id_data,
      Truncate(hashed_encrypted_id, params.encrypted_bucket_id_length()));
  return EncryptedBucketId::Create(encrypted_bucket_id_data,
                                   params.encrypted_bucket_id_length());
}

::rlwe::StatusOr<uint32_t> EncryptedBucketId::ToUint32() const {
  if (bit_length_ > 32) {
    return absl::InternalError("Bit length exceeds 32 bits.");
  }
  return TruncateAsUint32(encrypted_bucket_id_bytes_, bit_length_);
}

EncryptedBucketId::EncryptedBucketId(
    absl::string_view encrypted_bucket_id_bytes, int bit_length)
    : encrypted_bucket_id_bytes_(encrypted_bucket_id_bytes),
      bit_length_(bit_length) {}

}  // namespace rlwe
}  // namespace private_membership
