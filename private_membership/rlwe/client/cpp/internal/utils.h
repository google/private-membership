// Copyright 2020 Google LLC
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

#ifndef PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_UTILS_H_
#define PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_UTILS_H_

#include <cstdint>

#include "private_membership/rlwe/client/cpp/open_source/chromium/private_membership_export.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"

namespace private_membership {
namespace rlwe {

// Returns the first `bit_length` leftmost bits of `in` as an unsigned 32 bit
// integer.
//
// The input bytes `in` is assumed to be in big endian representation.
//
// For example, given 10111011 11111011 in big endian with `bit_length` == 12,
// the result will be 10111011 1111 which is equal to 3007 in decimal.
//
// Returns InvalidArgumentError if `bit_length` is larger than the bit length of
// the input string or larger than 32.
PRIVATE_MEMBERSHIP_EXPORT absl::StatusOr<uint32_t> TruncateAsUint32(
    absl::string_view in, int bit_length);

// Returns true if `in` is a valid byte representation of a value truncated to
// `bit_length`.
//
// That is, it returns true if and only if all bits after `bit_length`th bit is
// set to 0 and `bit_length` is in the range
// ((in.size() - 1) * 8, in.size() * 8].
PRIVATE_MEMBERSHIP_EXPORT bool IsValid(absl::string_view in, int bit_length);

}  // namespace rlwe
}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_UTILS_H_
