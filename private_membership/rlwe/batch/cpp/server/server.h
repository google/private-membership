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

#ifndef PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_SERVER_SERVER_H_
#define PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_SERVER_SERVER_H_

#include "absl/status/statusor.h"
#include "private_membership/rlwe/batch/proto/server.pb.h"

namespace private_membership {
namespace batch {

// Encodes raw databases to be ready for computation.
absl::StatusOr<EncodeDatabaseResponse> EncodeDatabase(
    const EncodeDatabaseRequest& request);

// Applies provided queries against provided encoded database shards.
absl::StatusOr<ApplyQueriesResponse> ApplyQueries(
    const ApplyQueriesRequest& request);

// Receives lists of lists of ciphertexts. Within each sublist, the function
// will sum all ciphertexts. The response will contain a list of all summed
// ciphertexts.
absl::StatusOr<SumCiphertextsResponse> SumCiphertexts(
    const SumCiphertextsRequest& request);

// Receives a list of results. The function finalizes the results so they can be
// returned to the querier. During finalization, the server performs modulus
// switching to compress the results.
absl::StatusOr<FinalizeResultsResponse> FinalizeResults(
    const FinalizeResultsRequest& request);

}  // namespace batch
}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_SERVER_SERVER_H_
