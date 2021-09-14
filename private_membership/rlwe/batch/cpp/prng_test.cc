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

#include <stdint.h>

#include <memory>
#include <string>

#include <gtest/gtest.h>
#include "absl/status/statusor.h"
#include "shell_encryption/prng/prng.h"

namespace private_membership {
namespace batch {
namespace {

TEST(PrngTest, DifferentOutput) {
  auto prng = CreatePrng();
  ASSERT_TRUE(prng.ok());

  absl::StatusOr<uint64_t> rand1 = (*prng)->Rand64();
  ASSERT_TRUE(rand1.ok());

  absl::StatusOr<uint64_t> rand2 = (*prng)->Rand64();
  ASSERT_TRUE(rand2.ok());

  EXPECT_NE(*rand1, *rand2);
}

TEST(PrngTest, RandomPrngs) {
  auto prng1 = CreatePrng();
  ASSERT_TRUE(prng1.ok());

  auto prng2 = CreatePrng();
  ASSERT_TRUE(prng2.ok());

  absl::StatusOr<uint64_t> rand1 = (*prng1)->Rand64();
  ASSERT_TRUE(rand1.ok());

  absl::StatusOr<uint64_t> rand2 = (*prng2)->Rand64();
  ASSERT_TRUE(rand2.ok());

  EXPECT_NE(*rand1, *rand2);
}

TEST(PrngTest, SameSeedPrngs) {
  auto prng1 = CreatePrngWithSeed();
  ASSERT_TRUE(prng1.ok());

  auto prng2 = CreatePrngFromSeed(prng1->seed);
  ASSERT_TRUE(prng2.ok());

  absl::StatusOr<uint64_t> rand1 = prng1->prng->Rand64();
  ASSERT_TRUE(rand1.ok());

  absl::StatusOr<uint64_t> rand2 = (*prng2)->Rand64();
  ASSERT_TRUE(rand2.ok());

  EXPECT_EQ(*rand1, *rand2);
}

}  // namespace
}  // namespace batch
}  // namespace private_membership
