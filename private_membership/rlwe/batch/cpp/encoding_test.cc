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

#include "private_membership/rlwe/batch/cpp/encoding.h"

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include "absl/status/statusor.h"
#include "private_membership/rlwe/batch/proto/client.pb.h"
#include "private_membership/rlwe/batch/cpp/constants.h"
#include "private_membership/rlwe/batch/cpp/context.h"
#include "private_membership/rlwe/batch/cpp/padding.h"
#include "private_membership/rlwe/batch/cpp/prng.h"
#include "private_membership/rlwe/batch/proto/shared.pb.h"
#include "private_membership/rlwe/batch/cpp/test_helper.h"
#include "shell_encryption/galois_key.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/serialization.pb.h"

namespace private_membership {
namespace batch {
namespace {

TEST(SharedTest, EncodeAndDecodeBytes) {
  auto context = CreateRlweRequestContext(CreateTestParameters());
  ASSERT_TRUE(context.ok());

  std::string bytes = "Test Byte Contents!";
  absl::StatusOr<std::vector<rlwe::Polynomial<ModularInt>>> encoded =
      EncodeBytes(**context, bytes);
  ASSERT_TRUE(encoded.ok());

  absl::StatusOr<std::string> decoded = DecodePolynomial(**context, *encoded);
  ASSERT_TRUE(decoded.ok());
  EXPECT_EQ(bytes, *decoded);
}

TEST(SharedTest, EncodeAndDecodeLongBytes) {
  auto context = CreateRlweRequestContext(CreateTestParameters());
  ASSERT_TRUE(context.ok());

  std::string bytes = std::string(6000, 'B');
  absl::StatusOr<std::vector<rlwe::Polynomial<ModularInt>>> encoded =
      EncodeBytes(**context, bytes);
  ASSERT_TRUE(encoded.ok());

  absl::StatusOr<std::string> decoded = DecodePolynomial(**context, *encoded);
  ASSERT_TRUE(decoded.ok());
  EXPECT_EQ(bytes, *decoded);
}

TEST(SharedTest, EncodeAndDecodeEmptyBytes) {
  auto context = CreateRlweRequestContext(CreateTestParameters());
  ASSERT_TRUE(context.ok());

  std::string bytes = "";
  absl::StatusOr<std::vector<rlwe::Polynomial<ModularInt>>> encoded =
      EncodeBytes(**context, bytes);
  ASSERT_TRUE(encoded.ok());

  absl::StatusOr<std::string> decoded = DecodePolynomial(**context, *encoded);
  ASSERT_TRUE(decoded.ok());
  EXPECT_EQ(bytes, *decoded);
}

TEST(SharedTest, DeserializePublicGaloisKey) {
  absl::StatusOr<GenerateKeysResponse> keys = CreateTestKeys();
  ASSERT_TRUE(keys.ok());

  auto context = CreateRlweRequestContext(CreateTestParameters());
  ASSERT_TRUE(context.ok());

  for (const rlwe::SerializedGaloisKey& serialized_key :
       keys->public_key().galois_key()) {
    absl::StatusOr<rlwe::GaloisKey<ModularInt>> public_key =
        DeserializePublicGaloisKey(serialized_key, **context);
    ASSERT_TRUE(public_key.ok());
  }
}

}  // namespace
}  // namespace batch
}  // namespace private_membership
