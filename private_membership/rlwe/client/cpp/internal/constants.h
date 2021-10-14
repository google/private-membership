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

#ifndef PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_CONSTANTS_H_
#define PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_CONSTANTS_H_

#include <cstdint>

#include <openssl/obj_mac.h>
#include "shell_encryption/montgomery.h"
#include "shell_encryption/prng/chacha_prng.h"
#include "shell_encryption/prng/single_thread_chacha_prng.h"

namespace private_membership {
namespace rlwe {

// Identifier of the elliptic curve.
constexpr int kCurveId = NID_X9_62_prime256v1;

// Byte length of the encrypted id stored inside buckets.
constexpr int kStoredEncryptedIdByteLength = 13;

// A PRNG that is not thread safe.
typedef ::rlwe::SingleThreadChaChaPrng SingleThreadPrng;

// A PRNG that is thread safe.
typedef ::rlwe::ChaChaPrng Prng;

// Modular int type for 64-bit moduli.
typedef ::rlwe::MontgomeryInt<uint64_t> ModularInt64;

// Modular int type for 128-bit moduli.
typedef ::rlwe::MontgomeryInt<absl::uint128> ModularInt128;

}  // namespace rlwe
}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_CONSTANTS_H_
