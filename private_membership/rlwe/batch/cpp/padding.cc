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

#include "private_membership/rlwe/batch/cpp/padding.h"

#include <cstdint>
#include <string>

#include <google/protobuf/repeated_field.h>
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
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

namespace {

std::string LengthToBytes(uint32_t length) {
  uint8_t ret[4];
  ret[0] = length & 0xFF;
  ret[1] = (length >> 8) & 0xFF;
  ret[2] = (length >> 16) & 0xFF;
  ret[3] = (length >> 24) & 0xFF;
  return std::string(ret, ret + 4);
}

absl::StatusOr<uint32_t> BytesToLength(absl::string_view bytes) {
  if (bytes.size() != 4) {
    return absl::InvalidArgumentError("Invalid byte size of length encoding.");
  }
  uint32_t ret = 0;
  ret += static_cast<uint8_t>(bytes[0]);
  ret += (static_cast<uint8_t>(bytes[1]) << 8);
  ret += (static_cast<uint8_t>(bytes[2]) << 16);
  ret += (static_cast<uint8_t>(bytes[3]) << 24);
  return ret;
}

}  // namespace

std::string PrependLength(absl::string_view input) {
  // Prepend byte length of value into the first four bytes. This limits the
  // maximum length of the value to at most 2^32 bytes or ~4.3 gigabytes.
  std::string length_in_bytes = LengthToBytes(input.length());
  return absl::StrCat(length_in_bytes, input);
}

absl::StatusOr<std::string> Unpad(absl::string_view input) {
  if (input.length() < 4) {
    return absl::InvalidArgumentError("Invalid input does not encode length.");
  }
  absl::StatusOr<uint32_t> input_length = BytesToLength(input.substr(0, 4));
  if (!input_length.ok()) {
    return input_length.status();
  }
  if (*input_length + 4 > input.length()) {
    return absl::InvalidArgumentError("Incorrect encoded length.");
  }
  return std::string(input.substr(4, *input_length));
}

}  // namespace batch
}  // namespace private_membership
