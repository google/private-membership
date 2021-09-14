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

#ifndef PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_ENCODING_H_
#define PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_ENCODING_H_

#include <string>
#include <vector>

#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "private_membership/rlwe/batch/cpp/constants.h"
#include "shell_encryption/context.h"
#include "shell_encryption/galois_key.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/serialization.pb.h"

namespace private_membership {
namespace batch {

// Encode bytes into polynomial.
absl::StatusOr<std::vector<rlwe::Polynomial<ModularInt>>> EncodeBytes(
    const Context& context, absl::string_view bytes);

// Decode polynomial into bytes.
absl::StatusOr<std::string> DecodePolynomial(
    const Context& context,
    const std::vector<rlwe::Polynomial<ModularInt>>& polynomial);

// Decodes public Galois key.
absl::StatusOr<rlwe::GaloisKey<ModularInt>> DeserializePublicGaloisKey(
    const rlwe::SerializedGaloisKey& serialized_key, const Context& context);

}  // namespace batch
}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_ENCODING_H_
