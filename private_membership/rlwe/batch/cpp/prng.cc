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

#include "private_membership/rlwe/batch/cpp/prng.h"

#include <memory>
#include <string>
#include <utility>

#include <google/protobuf/repeated_field.h>
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "private_membership/rlwe/batch/proto/shared.pb.h"
#include "shell_encryption/context.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/transcription.h"

namespace private_membership {
namespace batch {

absl::StatusOr<PrngWithSeed> CreatePrngWithSeed() {
  absl::StatusOr<std::string> seed = rlwe::SingleThreadHkdfPrng::GenerateSeed();
  if (!seed.ok()) {
    return seed.status();
  }
  auto prng = CreatePrngFromSeed(*seed);
  if (!prng.ok()) {
    return prng.status();
  }
  return PrngWithSeed{
      .prng = *std::move(prng),
      .seed = *std::move(seed),
  };
}

absl::StatusOr<std::unique_ptr<rlwe::SecurePrng>> CreatePrng() {
  auto prng_with_seed = CreatePrngWithSeed();
  if (!prng_with_seed.ok()) {
    return prng_with_seed.status();
  }
  return std::move(prng_with_seed->prng);
}

absl::StatusOr<std::unique_ptr<rlwe::SecurePrng>> CreatePrngFromSeed(
    absl::string_view seed) {
  return rlwe::SingleThreadHkdfPrng::Create(seed);
}

}  // namespace batch
}  // namespace private_membership
