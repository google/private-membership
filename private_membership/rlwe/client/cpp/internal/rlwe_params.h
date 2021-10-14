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

#ifndef PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_RLWE_PARAMS_H_
#define PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_RLWE_PARAMS_H_

#include "private_membership/rlwe/client/proto/private_membership_rlwe.pb.h"
#include "shell_encryption/context.h"
#include "shell_encryption/error_params.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/statusor.h"

namespace private_membership {
namespace rlwe {

template <typename ModularInt>
::rlwe::StatusOr<
    std::vector<std::unique_ptr<const ::rlwe::RlweContext<ModularInt>>>>
CreateContexts(const RlweParameters& rlwe_params);

template <typename ModularInt>
::rlwe::StatusOr<std::unique_ptr<const typename ModularInt::Params>>
CreateModulusParams(const Uint128& modulus);

template <typename ModularInt>
::rlwe::StatusOr<std::unique_ptr<const ::rlwe::NttParameters<ModularInt>>>
CreateNttParams(const RlweParameters& rlwe_params,
                const typename ModularInt::Params* modulus_params);

template <typename ModularInt>
::rlwe::StatusOr<std::unique_ptr<const ::rlwe::ErrorParams<ModularInt>>>
CreateErrorParams(const RlweParameters& rlwe_params,
                  const typename ModularInt::Params* modulus_params,
                  const ::rlwe::NttParameters<ModularInt>* ntt_params);

}  // namespace rlwe
}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_RLWE_PARAMS_H_
