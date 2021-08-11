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

#include "private_membership/rlwe/batch/cpp/shared.h"

#include <memory>

#include <google/protobuf/repeated_field.h>
#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "private_membership/rlwe/batch/cpp/shared.pb.h"
#include "context.h"
#include "integral_types.h"
#include "montgomery.h"
#include "prng/prng.h"
#include "prng/single_thread_hkdf_prng.h"

namespace private_membership {
namespace batch {
namespace {

absl::StatusOr<Int> CreateModulus(
    const google::protobuf::RepeatedField<uint64_t>& serialized_modulus) {
  if (serialized_modulus.empty()) {
    return absl::InvalidArgumentError("No modulus serialized.");
  } else if (serialized_modulus.size() > 2) {
    return absl::InvalidArgumentError("Modulus does not fit into uint128.");
  } else if (serialized_modulus.size() == 2) {
    return static_cast<rlwe::Uint128>(
        absl::MakeUint128(serialized_modulus.at(1), serialized_modulus.at(0)));
  } else {
    return static_cast<rlwe::Uint128>(
        absl::MakeUint128(0, serialized_modulus.at(0)));
  }
}

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

absl::StatusOr<std::unique_ptr<const rlwe::RlweContext<ModularInt>>>
CreateRlweRequestContext(const Parameters& parameters) {
  const auto& crypto_parameters = parameters.crypto_parameters();
  absl::StatusOr<Int> modulus =
      CreateModulus(crypto_parameters.request_modulus());
  if (!modulus.ok()) {
    return modulus.status();
  }
  return rlwe::RlweContext<ModularInt>::Create(
      {.modulus = *std::move(modulus),
       .log_n = static_cast<size_t>(crypto_parameters.log_degree()),
       .log_t = static_cast<size_t>(crypto_parameters.log_t()),
       .variance = static_cast<size_t>(crypto_parameters.variance())});
}

absl::StatusOr<std::unique_ptr<const rlwe::RlweContext<ModularInt>>>
CreateRlweResponseContext(const Parameters& parameters) {
  const auto& crypto_parameters = parameters.crypto_parameters();
  absl::StatusOr<Int> modulus =
      CreateModulus(crypto_parameters.response_modulus());
  if (!modulus.ok()) {
    return modulus.status();
  }
  return rlwe::RlweContext<ModularInt>::Create(
      {.modulus = *std::move(modulus),
       .log_n = static_cast<size_t>(crypto_parameters.log_degree()),
       .log_t = static_cast<size_t>(crypto_parameters.log_t()),
       .variance = static_cast<size_t>(crypto_parameters.variance())});
}

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
