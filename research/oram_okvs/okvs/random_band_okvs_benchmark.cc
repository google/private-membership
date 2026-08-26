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

#include <cstdint>
#include <random>
#include <string>
#include <vector>

#include "benchmark/benchmark.h"
#include "absl/strings/str_cat.h"
#include "research/oram_okvs/okvs/gf128.h"
#include "research/oram_okvs/okvs/random_band_okvs.h"

namespace oram_okvs::okvs {
namespace {

void BM_Encode(benchmark::State& state) {
  const size_t num_keys = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<std::string> keys(num_keys);
  std::vector<Gf128> values(num_keys);
  for (size_t i = 0; i < num_keys; ++i) {
    keys[i] = absl::StrCat("key_", i, "_", rng());
    values[i] = Gf128::Random(rng);
  }

  for (auto s : state) {
    RandomBandOkvs<std::string, Gf128, 384> okvs(/*epsilon=*/0.10, /*seed=*/1);
    auto status = okvs.Encode(keys, values);
    if (!status.ok()) {
      state.SkipWithError("Encoding failed");
      break;
    }
  }

  state.SetItemsProcessed(state.iterations() * num_keys);
}
BENCHMARK(BM_Encode)->RangeMultiplier(4)->Range(1024, 65536);

void BM_Decode(benchmark::State& state) {
  const size_t num_keys = state.range(0);
  std::mt19937_64 rng(42);

  std::vector<std::string> keys(num_keys);
  std::vector<Gf128> values(num_keys);
  for (size_t i = 0; i < num_keys; ++i) {
    keys[i] = absl::StrCat("key_", i, "_", rng());
    values[i] = Gf128::Random(rng);
  }

  RandomBandOkvs<std::string, Gf128, 384> okvs(/*epsilon=*/0.10, /*seed=*/1);
  auto status = okvs.Encode(keys, values);
  if (!status.ok()) {
    state.SkipWithError("Encoding failed in setup");
    return;
  }

  size_t idx = 0;
  for (auto s : state) {
    auto val = okvs.Decode(keys[idx]);
    benchmark::DoNotOptimize(val);
    idx = (idx + 1) % num_keys;
  }

  state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Decode)->RangeMultiplier(4)->Range(1024, 65536);

}  // namespace
}  // namespace oram_okvs::okvs
