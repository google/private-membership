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

#ifndef PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_PADDING_H_
#define PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_PADDING_H_

#include <string>

#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"

namespace private_membership {
namespace batch {

// Prepend length of input using the first 4 bytes in little-endian.
std::string PrependLength(absl::string_view input);

// Unpad string that has its real length prepended using the above function and
// is potentially padded with arbitrary characters afterwards.
absl::StatusOr<std::string> Unpad(absl::string_view input);

}  // namespace batch
}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_BATCH_CPP_PADDING_H_
