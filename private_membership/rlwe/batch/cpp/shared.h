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

#ifndef PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_CONSTANTS_H_
#define PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_CONSTANTS_H_

#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "private_membership/rlwe/batch/cpp/shared.pb.h"
#include "context.h"
#include "integral_types.h"
#include "montgomery.h"
#include "prng/prng.h"

namespace private_membership {
namespace batch {

using Int = rlwe::Uint128;
using ModularInt = rlwe::MontgomeryInt<Int>;
using Context = rlwe::RlweContext<ModularInt>;

// Creates context for generating and processing requests.
absl::StatusOr<std::unique_ptr<const Context>> CreateRlweRequestContext(
    const Parameters& parameters);

// Creates context for generating and processing responses.
absl::StatusOr<std::unique_ptr<const Context>> CreateRlweResponseContext(
    const Parameters& parameters);

// Create a PRNG with a random seed.
absl::StatusOr<std::unique_ptr<rlwe::SecurePrng>> CreatePrng();

// Struct holding a pair of PRNG and seed.
struct PrngWithSeed {
  std::unique_ptr<rlwe::SecurePrng> prng;
  std::string seed;
};

// Create a PRNG with a random seed and store seed in parameter.
absl::StatusOr<PrngWithSeed> CreatePrngWithSeed();

// Create a PRNG using the given seed.
absl::StatusOr<std::unique_ptr<rlwe::SecurePrng>> CreatePrngFromSeed(
    absl::string_view seed);

// Prepend length of input using the first 4 bytes in little-endian.
std::string PrependLength(absl::string_view input);

// Unpad string that has its real length prepended using the above function and
// is potentially padded with arbitrary characters afterwards.
absl::StatusOr<std::string> Unpad(absl::string_view input);

}  // namespace batch
}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_CONSTANTS_H_
