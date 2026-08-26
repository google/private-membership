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
#include <vector>

#include "benchmark/benchmark.h"
#include "research/oram_okvs/oram/lat/latoram.h"
#include "research/oram_okvs/oram/lat/multi_perm.h"
#include "research/oram_okvs/oram/types.h"

namespace oram_okvs {
namespace {

void BM_MultiPermNext(benchmark::State& state) {
  size_t num_bins = state.range(0);
  std::vector<size_t> freqs(num_bins, 1000000);
  MultiPermDS ds(freqs, 42);

  for (auto _ : state) {
    auto z_or = ds.Next();
    benchmark::DoNotOptimize(z_or);
  }
}
BENCHMARK(BM_MultiPermNext)->Arg(2)->Arg(8)->Arg(32)->Arg(128);

void BM_LatOramRead(benchmark::State& state) {
  size_t n = state.range(0);
  size_t block_size = state.range(1);

  LatOramConfig config;
  config.num_blocks = n;
  config.block_size = block_size;
  config.seed = 42;

  LatOram oram(config);
  auto init_status = oram.Init();
  if (!init_status.ok()) {
    state.SkipWithError("LatOram Init failed.");
    return;
  }

  std::mt19937_64 rng(12345);
  for (auto _ : state) {
    BlockId id = rng() % n;
    auto res = oram.Read(id);
    benchmark::DoNotOptimize(res);
  }
}
BENCHMARK(BM_LatOramRead)
    ->Args({256, 64})
    ->Args({1024, 64})
    ->Args({1024, 256})
    ->Args({4096, 64});

void BM_LatOramWrite(benchmark::State& state) {
  size_t n = state.range(0);
  size_t block_size = state.range(1);

  LatOramConfig config;
  config.num_blocks = n;
  config.block_size = block_size;
  config.seed = 42;

  LatOram oram(config);
  auto init_status = oram.Init();
  if (!init_status.ok()) {
    state.SkipWithError("LatOram Init failed.");
    return;
  }

  std::vector<uint8_t> data(block_size, 0xAB);
  std::mt19937_64 rng(12345);
  for (auto _ : state) {
    BlockId id = rng() % n;
    auto res = oram.Write(id, data);
    benchmark::DoNotOptimize(res);
  }
}
void BM_OkvsLatOramRead(benchmark::State& state) {
  size_t n = state.range(0);
  size_t block_size = state.range(1);

  LatOramConfig config;
  config.num_blocks = n;
  config.block_size = block_size;
  config.storage_type = LatOramStorageType::kOkLat;
  config.seed = 42;

  LatOram oram(config);
  auto init_status = oram.Init();
  if (!init_status.ok()) {
    state.SkipWithError("OkvsLatOram Init failed.");
    return;
  }

  std::mt19937_64 rng(12345);
  for (auto _ : state) {
    BlockId id = rng() % n;
    auto res = oram.Read(id);
    benchmark::DoNotOptimize(res);
  }
}
BENCHMARK(BM_OkvsLatOramRead)
    ->Args({256, 64})
    ->Args({1024, 64})
    ->Args({1024, 256});

void BM_OkvsLatOramWrite(benchmark::State& state) {
  size_t n = state.range(0);
  size_t block_size = state.range(1);

  LatOramConfig config;
  config.num_blocks = n;
  config.block_size = block_size;
  config.storage_type = LatOramStorageType::kOkLat;
  config.seed = 42;

  LatOram oram(config);
  auto init_status = oram.Init();
  if (!init_status.ok()) {
    state.SkipWithError("OkvsLatOram Init failed.");
    return;
  }

  std::vector<uint8_t> data(block_size, 0xAB);
  std::mt19937_64 rng(12345);
  for (auto _ : state) {
    BlockId id = rng() % n;
    auto res = oram.Write(id, data);
    benchmark::DoNotOptimize(res);
  }
}
BENCHMARK(BM_OkvsLatOramWrite)
    ->Args({256, 64})
    ->Args({1024, 64})
    ->Args({1024, 256});

}  // namespace
}  // namespace oram_okvs
