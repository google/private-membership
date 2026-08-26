// Copyright 2026 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "absl/strings/str_cat.h"
#include "benchmark/benchmark.h"
#include "research/oram_okvs/pir/fast_ntt.h"
#include "shell_encryption/dft_transformations.h"
#include "shell_encryption/dft_transformations_hwy.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"

namespace oram_okvs::pir {
namespace {

using ModularInt = rlwe::MontgomeryInt<uint64_t>;
using PrimeModulus = rlwe::PrimeModulus<ModularInt>;
using RnsContext = rlwe::RnsContext<ModularInt>;

std::vector<ModularInt> GenerateTestCoeffs(
    size_t n, const ModularInt::Params* mod_params) {
  std::vector<ModularInt> coeffs;
  coeffs.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    uint64_t val = (i * 37 + 101) % mod_params->modulus;
    coeffs.push_back(ModularInt::ImportInt(val, mod_params).value());
  }
  return coeffs;
}

// ----------------------------------------------------------------------------
// Forward NTT Benchmarks
// ----------------------------------------------------------------------------

// 1. RLWE Baseline Scalar Forward NTT
void BM_Forward_RlweScalar(benchmark::State& state) {
  int log_n = state.range(0);
  int n = 1 << log_n;
  auto context =
      RnsContext::CreateForBfv(log_n, {1152921504606830593ULL}, {}, 65537)
          .value();
  const auto* prime_mod = context.MainPrimeModuli()[0];
  const auto* ntt_params = prime_mod->NttParams();
  const auto* mod_params = prime_mod->ModParams();

  auto init_coeffs = GenerateTestCoeffs(n, mod_params);
  std::vector<ModularInt> coeffs = init_coeffs;

  for (auto _ : state) {
    coeffs = init_coeffs;
    benchmark::DoNotOptimize(rlwe::ForwardNumberTheoreticTransform(
        coeffs, *ntt_params, *mod_params));
  }
  state.SetBytesProcessed(state.iterations() * n * sizeof(uint64_t));
  state.SetItemsProcessed(state.iterations() * n);
}
BENCHMARK(BM_Forward_RlweScalar)->Arg(9)->Arg(10)->Arg(11)->Arg(12);

// 2. RLWE Existing Highway Forward NTT
void BM_Forward_RlweHighway(benchmark::State& state) {
  int log_n = state.range(0);
  int n = 1 << log_n;
  auto context =
      RnsContext::CreateForBfv(log_n, {1152921504606830593ULL}, {}, 65537)
          .value();
  const auto* prime_mod = context.MainPrimeModuli()[0];
  const auto* ntt_params = prime_mod->NttParams();
  const auto* mod_params = prime_mod->ModParams();

  auto init_coeffs = GenerateTestCoeffs(n, mod_params);
  std::vector<ModularInt> coeffs = init_coeffs;

  for (auto _ : state) {
    coeffs = init_coeffs;
    benchmark::DoNotOptimize(rlwe::internal::ForwardNumberTheoreticTransformHwy(
        coeffs, *ntt_params, *mod_params));
  }
  state.SetBytesProcessed(state.iterations() * n * sizeof(uint64_t));
  state.SetItemsProcessed(state.iterations() * n);
}
BENCHMARK(BM_Forward_RlweHighway)->Arg(9)->Arg(10)->Arg(11)->Arg(12);

// 3. Our Custom Optimized Forward NTT (4x Unrolled + Harvey Lazy Reductions)
void BM_Forward_CustomFast(benchmark::State& state) {
  int log_n = state.range(0);
  int n = 1 << log_n;
  auto context =
      RnsContext::CreateForBfv(log_n, {1152921504606830593ULL}, {}, 65537)
          .value();
  const auto* prime_mod = context.MainPrimeModuli()[0];
  const auto* ntt_params = prime_mod->NttParams();
  const auto* mod_params = prime_mod->ModParams();

  auto init_coeffs = GenerateTestCoeffs(n, mod_params);
  std::vector<ModularInt> coeffs = init_coeffs;

  for (auto _ : state) {
    coeffs = init_coeffs;
    benchmark::DoNotOptimize(FastForwardNtt(coeffs, *ntt_params, *mod_params));
  }
  state.SetBytesProcessed(state.iterations() * n * sizeof(uint64_t));
  state.SetItemsProcessed(state.iterations() * n);
}
BENCHMARK(BM_Forward_CustomFast)->Arg(9)->Arg(10)->Arg(11)->Arg(12);

// ----------------------------------------------------------------------------
// Inverse NTT Benchmarks
// ----------------------------------------------------------------------------

// 1. RLWE Baseline Scalar Inverse NTT
void BM_Inverse_RlweScalar(benchmark::State& state) {
  int log_n = state.range(0);
  int n = 1 << log_n;
  auto context =
      RnsContext::CreateForBfv(log_n, {1152921504606830593ULL}, {}, 65537)
          .value();
  const auto* prime_mod = context.MainPrimeModuli()[0];
  const auto* ntt_params = prime_mod->NttParams();
  const auto* mod_params = prime_mod->ModParams();

  auto init_coeffs = GenerateTestCoeffs(n, mod_params);
  std::vector<ModularInt> coeffs = init_coeffs;

  for (auto _ : state) {
    coeffs = init_coeffs;
    benchmark::DoNotOptimize(rlwe::InverseNumberTheoreticTransform(
        coeffs, *ntt_params, *mod_params));
  }
  state.SetBytesProcessed(state.iterations() * n * sizeof(uint64_t));
  state.SetItemsProcessed(state.iterations() * n);
}
BENCHMARK(BM_Inverse_RlweScalar)->Arg(9)->Arg(10)->Arg(11)->Arg(12);

// 2. Our Custom Optimized Inverse NTT (Gentleman-Sande + Highway SIMD)
void BM_Inverse_CustomFast(benchmark::State& state) {
  int log_n = state.range(0);
  int n = 1 << log_n;
  auto context =
      RnsContext::CreateForBfv(log_n, {1152921504606830593ULL}, {}, 65537)
          .value();
  const auto* prime_mod = context.MainPrimeModuli()[0];
  const auto* ntt_params = prime_mod->NttParams();
  const auto* mod_params = prime_mod->ModParams();

  auto init_coeffs = GenerateTestCoeffs(n, mod_params);
  std::vector<ModularInt> coeffs = init_coeffs;

  for (auto _ : state) {
    coeffs = init_coeffs;
    benchmark::DoNotOptimize(FastInverseNtt(coeffs, *ntt_params, *mod_params));
  }
  state.SetBytesProcessed(state.iterations() * n * sizeof(uint64_t));
  state.SetItemsProcessed(state.iterations() * n);
}
BENCHMARK(BM_Inverse_CustomFast)->Arg(9)->Arg(10)->Arg(11)->Arg(12);

}  // namespace
}  // namespace oram_okvs::pir
