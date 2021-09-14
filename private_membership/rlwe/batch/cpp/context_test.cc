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

#include <gtest/gtest.h>
#include "absl/status/statusor.h"
#include "private_membership/rlwe/batch/proto/shared.pb.h"
#include "private_membership/rlwe/batch/cpp/test_helper.h"

namespace private_membership {
namespace batch {
namespace {

TEST(ContextTest, CreateRequestContext) {
  auto context = CreateRlweRequestContext(CreateTestParameters());
  EXPECT_TRUE(context.ok());
}

TEST(ContextTest, CreateResponseContext) {
  auto context = CreateRlweResponseContext(CreateTestParameters());
  EXPECT_TRUE(context.ok());
}

TEST(ContextTest, EmptyParameters) {
  Parameters parameters;

  auto request_context = CreateRlweRequestContext(parameters);
  EXPECT_FALSE(request_context.ok());

  auto response_context = CreateRlweResponseContext(parameters);
  EXPECT_FALSE(response_context.ok());
}

}  // namespace
}  // namespace batch
}  // namespace private_membership
