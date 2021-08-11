// Copyright 2021 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "private_membership/rlwe/batch/cpp/shared.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace private_membership {
namespace batch {
namespace {

TEST(SharedTest, PrependLengthSuccess) {
  std::string expected_output("\x05\x00\x00\x00input", 9);
  EXPECT_EQ(PrependLength("input"), expected_output);
}

TEST(SharedTest, UnpadByteLengthSmallerThan4) {
  EXPECT_FALSE(Unpad("abc").ok());
}

TEST(SharedTest, UnpadIncorrectByteLength) {
  std::string input("\x06\x00\x00\x00input", 9);
  EXPECT_FALSE(Unpad(input).ok());
}

TEST(SharedTest, UnpadSuccess) {
  std::string input("\x05\x00\x00\x00input\x00\x00", 11);
  absl::StatusOr<std::string> x = Unpad(input);
  ASSERT_TRUE(x.ok());
  EXPECT_EQ(*x, "input");
}

TEST(SharedTest, UnpadSuccessWithPrependLength) {
  std::string input = absl::StrCat(PrependLength("input"), "abc");
  absl::StatusOr<std::string> x = Unpad(input);
  ASSERT_TRUE(x.ok());
  EXPECT_EQ(*x, "input");
}

}  // namespace
}  // namespace batch
}  // namespace private_membership
