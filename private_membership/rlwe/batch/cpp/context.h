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

#ifndef PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_CONTEXT_H_
#define PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_CONTEXT_H_

#include <memory>

#include "absl/status/statusor.h"
#include "private_membership/rlwe/batch/cpp/constants.h"
#include "private_membership/rlwe/batch/proto/shared.pb.h"
#include "shell_encryption/context.h"

namespace private_membership {
namespace batch {

// Creates context for generating and processing requests.
absl::StatusOr<std::unique_ptr<const Context>> CreateRlweRequestContext(
    const Parameters& parameters);

// Creates context for generating and processing responses.
absl::StatusOr<std::unique_ptr<const Context>> CreateRlweResponseContext(
    const Parameters& parameters);

}  // namespace batch
}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_CONTEXT_H_
