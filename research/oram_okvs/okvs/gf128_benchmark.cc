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

#include <random>
#include <vector>

#include "benchmark/benchmark.h"
#include "absl/base/optimization.h"
#include "research/oram_okvs/okvs/gf128.h"

namespace oram_okvs::okvs {
namespace {

constexpr size_t kBatchSize = 1024;

void BM_Gf128Add(benchmark::State& state) {
  std::mt19937_64 rng(42);
  std::vector<Gf128> a(kBatchSize), b(kBatchSize), c(kBatchSize);
  for (size_t i = 0; i < kBatchSize; ++i) {
    a[i] = Gf128::Random(rng);
    b[i] = Gf128::Random(rng);
  }

  size_t idx = 0;
  for (auto s : state) {
    c[idx] = a[idx] + b[idx];
    benchmark::DoNotOptimize(c[idx]);
    idx = (idx + 1) % kBatchSize;
  }
  state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Gf128Add);

void BM_Gf128Multiply(benchmark::State& state) {
  std::mt19937_64 rng(42);
  std::vector<Gf128> a(kBatchSize), b(kBatchSize), c(kBatchSize);
  for (size_t i = 0; i < kBatchSize; ++i) {
    a[i] = Gf128::Random(rng);
    b[i] = Gf128::Random(rng);
  }

  size_t idx = 0;
  for (auto s : state) {
    c[idx] = a[idx] * b[idx];
    benchmark::DoNotOptimize(c[idx]);
    idx = (idx + 1) % kBatchSize;
  }
  state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Gf128Multiply);

void BM_Gf128Square(benchmark::State& state) {
  std::mt19937_64 rng(42);
  std::vector<Gf128> a(kBatchSize), c(kBatchSize);
  for (size_t i = 0; i < kBatchSize; ++i) {
    a[i] = Gf128::Random(rng);
  }

  size_t idx = 0;
  for (auto s : state) {
    c[idx] = a[idx].Square();
    benchmark::DoNotOptimize(c[idx]);
    idx = (idx + 1) % kBatchSize;
  }
  state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Gf128Square);

void BM_Gf128Inverse(benchmark::State& state) {
  std::mt19937_64 rng(42);
  std::vector<Gf128> a(kBatchSize), c(kBatchSize);
  for (size_t i = 0; i < kBatchSize; ++i) {
    a[i] = Gf128::Random(rng);
  }

  size_t idx = 0;
  for (auto s : state) {
    c[idx] = a[idx].Inverse();
    benchmark::DoNotOptimize(c[idx]);
    idx = (idx + 1) % kBatchSize;
  }
  state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Gf128Inverse);

void BM_Gf128Divide(benchmark::State& state) {
  std::mt19937_64 rng(42);
  std::vector<Gf128> a(kBatchSize), b(kBatchSize), c(kBatchSize);
  for (size_t i = 0; i < kBatchSize; ++i) {
    a[i] = Gf128::Random(rng);
    b[i] = Gf128::Random(rng);
  }

  size_t idx = 0;
  for (auto s : state) {
    c[idx] = a[idx] / b[idx];
    benchmark::DoNotOptimize(c[idx]);
    idx = (idx + 1) % kBatchSize;
  }
  state.SetItemsProcessed(state.iterations());
}
BENCHMARK(BM_Gf128Divide);

}  // namespace
}  // namespace oram_okvs::okvs
