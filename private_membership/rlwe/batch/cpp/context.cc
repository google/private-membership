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

#include "private_membership/rlwe/batch/cpp/context.h"

#include <stddef.h>
#include <stdint.h>

#include <memory>
#include <utility>

#include <google/protobuf/repeated_field.h>
#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "private_membership/rlwe/batch/cpp/constants.h"
#include "private_membership/rlwe/batch/proto/shared.pb.h"
#include "shell_encryption/context.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"


namespace private_membership {
namespace batch {

namespace {

absl::StatusOr<Int> CreateModulus(
    const google::protobuf::RepeatedField<uint64_t>& serialized_modulus) {
  if (serialized_modulus.empty()) {
    return absl::InvalidArgumentError("No modulus serialized.");
  } else if (serialized_modulus.size() > 2) {
    return absl::InvalidArgumentError("Modulus does not fit into uint128.");
  } else if (serialized_modulus.size() == 2) {
    return static_cast<rlwe::Uint128>(
        absl::MakeUint128(serialized_modulus.at(1), serialized_modulus.at(0)));
  } else {
    return static_cast<rlwe::Uint128>(
        absl::MakeUint128(0, serialized_modulus.at(0)));
  }
}

}  // namespace

absl::StatusOr<std::unique_ptr<const rlwe::RlweContext<ModularInt>>>
CreateRlweRequestContext(const Parameters& parameters) {
  const auto& crypto_parameters = parameters.crypto_parameters();
  absl::StatusOr<Int> modulus =
      CreateModulus(crypto_parameters.request_modulus());
  if (!modulus.ok()) {
    return modulus.status();
  }
  return rlwe::RlweContext<ModularInt>::Create(
      {.modulus = *std::move(modulus),
       .log_n = static_cast<size_t>(crypto_parameters.log_degree()),
       .log_t = static_cast<size_t>(crypto_parameters.log_t()),
       .variance = static_cast<size_t>(crypto_parameters.variance())});
}

absl::StatusOr<std::unique_ptr<const rlwe::RlweContext<ModularInt>>>
CreateRlweResponseContext(const Parameters& parameters) {
  const auto& crypto_parameters = parameters.crypto_parameters();
  absl::StatusOr<Int> modulus =
      CreateModulus(crypto_parameters.response_modulus());
  if (!modulus.ok()) {
    return modulus.status();
  }
  return rlwe::RlweContext<ModularInt>::Create(
      {.modulus = *std::move(modulus),
       .log_n = static_cast<size_t>(crypto_parameters.log_degree()),
       .log_t = static_cast<size_t>(crypto_parameters.log_t()),
       .variance = static_cast<size_t>(crypto_parameters.variance())});
}

}  // namespace batch
}  // namespace private_membership
