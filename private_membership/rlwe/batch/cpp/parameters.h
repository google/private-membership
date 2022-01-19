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

#ifndef PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_PARAMETERS_H_
#define PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_PARAMETERS_H_

#include "absl/status/status.h"
#include "private_membership/rlwe/batch/proto/shared.pb.h"

namespace private_membership {
namespace batch {

// Ensures values in `parameters` are within acceptable ranges.
absl::Status ValidateParameters(const Parameters& parameters);

}  // namespace batch
}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_PARAMETERS_H_
