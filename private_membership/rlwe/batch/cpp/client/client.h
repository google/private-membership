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

#ifndef PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_CLIENT_CLIENT_H_
#define PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_CLIENT_CLIENT_H_

#include "absl/status/statusor.h"
#include "private_membership/rlwe/batch/proto/client.pb.h"

namespace private_membership {
namespace batch {

absl::StatusOr<GenerateKeysResponse> GenerateKeys(
    const GenerateKeysRequest& request);

absl::StatusOr<EncryptQueriesResponse> EncryptQueries(
    const EncryptQueriesRequest& request);

absl::StatusOr<DecryptQueriesResponse> DecryptQueries(
    const DecryptQueriesRequest& request);

}  // namespace batch
}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_CLIENT_CLIENT_H_
