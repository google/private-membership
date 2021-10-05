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

#ifndef PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_RLWE_ID_UTILS_H_
#define PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_RLWE_ID_UTILS_H_

#include <string>

#include "private_join_and_compute/crypto/ec_commutative_cipher.h"
#include "private_membership/rlwe/client/proto/private_membership.pb.h"
#include "private_membership/rlwe/client/proto/private_membership_rlwe.pb.h"
#include "shell_encryption/statusor.h"

namespace private_membership {
namespace rlwe {

// Computes the representation of an encrypted ID stored within buckets starting
// from the plaintext.
//
// Returns an error when the cryptographic functions fails or the parameters are
// invalid.
::rlwe::StatusOr<std::string> ComputeBucketStoredEncryptedId(
    const RlwePlaintextId& id, const EncryptedBucketsParameters& params,
    ::private_join_and_compute::ECCommutativeCipher* ec_cipher, ::private_join_and_compute::Context* ctx);

// Computes the representation of an encrypted ID stored within buckets using
// the encrypted ID.
//
// Returns an error when the cryptographic functions fails or the parameters are
// invalid.
::rlwe::StatusOr<std::string> ComputeBucketStoredEncryptedId(
    absl::string_view encrypted_id, const EncryptedBucketsParameters& params,
    ::private_join_and_compute::Context* ctx);

// Function used to injectively hash RlwePlaintextId proto to string. This hash
// is not cryptographically secure, nor very compact.
std::string HashRlwePlaintextId(const RlwePlaintextId& id);

// Function used to hash the nonsensitive portion of a RlwePlaintextId, using a
// salt to force adversaries to recompute rainbow tables.
::rlwe::StatusOr<std::string> HashNonsensitiveIdWithSalt(
    absl::string_view nsid, HashType hash_type, ::private_join_and_compute::Context* ctx);

}  // namespace rlwe
}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_RLWE_ID_UTILS_H_
