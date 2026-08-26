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

#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "absl/strings/str_cat.h"
#include "benchmark/benchmark.h"
#include "research/oram_okvs/okvs/gf128.h"
#include "research/oram_okvs/okvs/polynomial_okvs.h"

namespace oram_okvs::okvs {
namespace {

void BM_PolynomialOkvs_Decode(benchmark::State& state) {
  size_t Z = state.range(0);
  size_t B = state.range(1);

  std::mt19937_64 rng(42 + Z);
  std::vector<std::string> keys(Z);
  std::vector<std::vector<uint8_t>> values(Z, std::vector<uint8_t>(B));
  for (size_t i = 0; i < Z; ++i) {
    keys[i] = absl::StrCat("key_", i, "_", rng());
    for (size_t b = 0; b < B; ++b) {
      values[i][b] = static_cast<uint8_t>(rng());
    }
  }

  PolynomialOkvs<std::string, std::vector<uint8_t>> okvs(Z, 42, B);
  (void)okvs.Encode(keys, values);

  size_t idx = 0;
  for (auto _ : state) {
    auto val_or = okvs.Decode(keys[idx % Z]);
    benchmark::DoNotOptimize(val_or);
    ++idx;
  }
  state.SetBytesProcessed(state.iterations() * B);
}

BENCHMARK(BM_PolynomialOkvs_Decode)
    ->Args({100, 1024})
    ->Args({100, 2048})
    ->Args({100, 4096})
    ->Args({435, 1024})
    ->Args({435, 2048})
    ->Args({435, 4096});

void BM_PolynomialOkvs_Encode(benchmark::State& state) {
  size_t Z = state.range(0);
  size_t B = state.range(1);

  std::mt19937_64 rng(42 + Z);
  std::vector<std::string> keys(Z);
  std::vector<std::vector<uint8_t>> values(Z, std::vector<uint8_t>(B));
  for (size_t i = 0; i < Z; ++i) {
    keys[i] = absl::StrCat("key_", i, "_", rng());
    for (size_t b = 0; b < B; ++b) {
      values[i][b] = static_cast<uint8_t>(rng());
    }
  }

  PolynomialOkvs<std::string, std::vector<uint8_t>> okvs(Z, 42, B);

  for (auto _ : state) {
    auto status = okvs.Encode(keys, values);
    benchmark::DoNotOptimize(status);
  }
  state.SetBytesProcessed(state.iterations() * Z * B);
}

BENCHMARK(BM_PolynomialOkvs_Encode)
    ->Args({100, 1024})
    ->Args({100, 2048})
    ->Args({100, 4096})
    ->Args({435, 1024})
    ->Args({435, 2048})
    ->Args({435, 4096});

}  // namespace
}  // namespace oram_okvs::okvs
