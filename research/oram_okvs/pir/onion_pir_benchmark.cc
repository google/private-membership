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
#include "research/oram_okvs/pir/onion_pir_client.h"
#include "research/oram_okvs/pir/onion_pir_params.h"
#include "research/oram_okvs/pir/onion_pir_server.h"
#include "research/oram_okvs/pir/onion_pir_types.h"
#include "shell_encryption/dft_transformations.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"

namespace oram_okvs::pir {
namespace {

std::vector<std::vector<uint8_t>> GenerateBenchmarkDb(size_t num_entries,
                                                      size_t entry_size) {
  std::vector<std::vector<uint8_t>> db(num_entries);
  for (size_t i = 0; i < num_entries; ++i) {
    db[i].resize(entry_size);
    for (size_t j = 0; j < entry_size; ++j) {
      db[i][j] = static_cast<uint8_t>((i * 37 + j * 17) & 0xFF);
    }
  }
  return db;
}

void BM_ClientQueryPack(benchmark::State& state) {
  size_t num_entries = state.range(0);
  size_t entry_size = 64;
  OnionPirParams params =
      OnionPirParams::CreateForOram(num_entries, entry_size);
  auto client = OnionPirClient::Create(params).value();

  size_t query_idx = num_entries / 2;
  for (auto _ : state) {
    auto query = client->CreateQuery(query_idx).value();
    benchmark::DoNotOptimize(query);
  }
  state.SetBytesProcessed(state.iterations() * params.poly_degree *
                          sizeof(uint64_t));
}

BENCHMARK(BM_ClientQueryPack)->Arg(64)->Arg(256)->Arg(1024);

void BM_ServerQueryUnpack(benchmark::State& state) {
  size_t num_entries = state.range(0);
  size_t entry_size = 64;
  OnionPirParams params =
      OnionPirParams::CreateForOram(num_entries, entry_size);
  auto client = OnionPirClient::Create(params).value();
  auto server = OnionPirServer::Create(params).value();

  auto eval_keys = client->GetSerializedEvaluationKeys().value();
  server->SetEvaluationKeys(eval_keys).IgnoreError();

  auto query = client->CreatePackedQueryCiphertext(num_entries / 2).value();

  for (auto _ : state) {
    auto unpacked = server->UnpackQuery(query).value();
    benchmark::DoNotOptimize(unpacked);
  }
}

BENCHMARK(BM_ServerQueryUnpack)->Arg(64)->Arg(256)->Arg(1024);

void BM_ServerProcessQueryPrecomputed(benchmark::State& state) {
  size_t num_entries = state.range(0);
  size_t entry_size = 64;
  OnionPirParams params =
      OnionPirParams::CreateForOram(num_entries, entry_size);
  auto client = OnionPirClient::Create(params).value();
  auto server = OnionPirServer::Create(params).value();

  auto eval_keys = client->GetSerializedEvaluationKeys().value();
  server->SetEvaluationKeys(eval_keys).IgnoreError();

  auto db = GenerateBenchmarkDb(num_entries, entry_size);
  server->SetDatabase(db).IgnoreError();

  auto query = client->CreatePackedQueryCiphertext(num_entries / 2).value();

  for (auto _ : state) {
    auto resp = server->ProcessQueryDirect(query, /*db=*/{}).value();
    benchmark::DoNotOptimize(resp);
  }

  int64_t total_db_bytes = static_cast<int64_t>(num_entries * entry_size);
  state.SetBytesProcessed(state.iterations() * total_db_bytes);
}

BENCHMARK(BM_ServerProcessQueryPrecomputed)->Arg(64)->Arg(256)->Arg(1024);

void BM_ServerProcessQueryOnTheFly(benchmark::State& state) {
  size_t num_entries = state.range(0);
  size_t entry_size = 64;
  OnionPirParams params =
      OnionPirParams::CreateForOram(num_entries, entry_size);
  auto client = OnionPirClient::Create(params).value();
  auto server = OnionPirServer::Create(params).value();

  auto eval_keys = client->GetSerializedEvaluationKeys().value();
  server->SetEvaluationKeys(eval_keys).IgnoreError();

  auto db = GenerateBenchmarkDb(num_entries, entry_size);
  auto query = client->CreatePackedQueryCiphertext(num_entries / 2).value();

  for (auto _ : state) {
    auto resp = server->ProcessQueryDirect(query, db).value();
    benchmark::DoNotOptimize(resp);
  }

  int64_t total_db_bytes = static_cast<int64_t>(num_entries * entry_size);
  state.SetBytesProcessed(state.iterations() * total_db_bytes);
}

BENCHMARK(BM_ServerProcessQueryOnTheFly)->Arg(64)->Arg(256)->Arg(1024);

void BM_ClientDecodeResponse(benchmark::State& state) {
  size_t num_entries = 1024;
  size_t entry_size = state.range(0);
  OnionPirParams params =
      OnionPirParams::CreateForOram(num_entries, entry_size);
  auto client = OnionPirClient::Create(params).value();
  auto server = OnionPirServer::Create(params).value();

  auto eval_keys = client->GetSerializedEvaluationKeys().value();
  server->SetEvaluationKeys(eval_keys).IgnoreError();

  auto db = GenerateBenchmarkDb(num_entries, entry_size);
  auto query = client->CreatePackedQueryCiphertext(num_entries / 2).value();
  auto resp_cts = server->ProcessQueryDirect(query, db).value();

  for (auto _ : state) {
    auto decoded =
        client->DecryptResponseCiphertexts(resp_cts, entry_size).value();
    benchmark::DoNotOptimize(decoded);
  }
  state.SetBytesProcessed(state.iterations() * entry_size);
}

BENCHMARK(BM_ClientDecodeResponse)->Arg(64)->Arg(256)->Arg(1024)->Arg(3072);

// ============================================================================
// Benchmarks for Varying N0 on N = 1024 (4 MB Database)
// ============================================================================

void BM_VaryingN0_QueryUnpack(benchmark::State& state) {
  size_t num_entries = 1024;
  size_t n0_val = state.range(0);
  size_t entry_size = 4096;  // 4 KB

  OnionPirParams params =
      OnionPirParams::CreateForOram(num_entries, entry_size);
  params.n0 = n0_val;
  auto client = OnionPirClient::Create(params).value();
  auto server = OnionPirServer::Create(params).value();

  auto eval_keys = client->GetSerializedEvaluationKeys().value();
  server->SetEvaluationKeys(eval_keys).IgnoreError();

  auto query = client->CreatePackedQueryCiphertext(num_entries / 2).value();

  for (auto _ : state) {
    auto unpacked = server->UnpackQuery(query).value();
    benchmark::DoNotOptimize(unpacked);
  }
}

BENCHMARK(BM_VaryingN0_QueryUnpack)
    ->Arg(32)
    ->Arg(64)
    ->Arg(128)
    ->Arg(256)
    ->Arg(512)
    ->Arg(1024);

void BM_VaryingN0_ProcessQueryOnTheFly(benchmark::State& state) {
  size_t num_entries = 1024;
  size_t n0_val = state.range(0);
  size_t entry_size = 4096;  // 4 KB per entry -> 4.19 MB DB

  OnionPirParams params =
      OnionPirParams::CreateForOram(num_entries, entry_size);
  params.n0 = n0_val;
  auto client = OnionPirClient::Create(params).value();
  auto server = OnionPirServer::Create(params).value();

  auto eval_keys = client->GetSerializedEvaluationKeys().value();
  server->SetEvaluationKeys(eval_keys).IgnoreError();

  auto db = GenerateBenchmarkDb(num_entries, entry_size);
  auto query = client->CreatePackedQueryCiphertext(num_entries / 2).value();

  for (auto _ : state) {
    auto resp = server->ProcessQueryDirect(query, db).value();
    benchmark::DoNotOptimize(resp);
  }

  int64_t total_db_bytes = static_cast<int64_t>(num_entries * entry_size);
  state.SetBytesProcessed(state.iterations() * total_db_bytes);
}

BENCHMARK(BM_VaryingN0_ProcessQueryOnTheFly)
    ->Arg(32)
    ->Arg(64)
    ->Arg(128)
    ->Arg(256)
    ->Arg(512)
    ->Arg(1024);

// ============================================================================
// Benchmarks for Database Scaling with Optimal N0
// ============================================================================

void BM_ScalingDbSize_OptimalN0(benchmark::State& state) {
  size_t num_entries = state.range(0);
  size_t entry_size = 4096;  // 4 KB

  OnionPirParams params =
      OnionPirParams::CreateForOram(num_entries, entry_size);
  params.n0 = std::min<size_t>(num_entries, 1024);
  auto client = OnionPirClient::Create(params).value();
  auto server = OnionPirServer::Create(params).value();

  auto eval_keys = client->GetSerializedEvaluationKeys().value();
  server->SetEvaluationKeys(eval_keys).IgnoreError();

  auto db = GenerateBenchmarkDb(num_entries, entry_size);
  auto query = client->CreatePackedQueryCiphertext(num_entries / 2).value();

  for (auto _ : state) {
    auto resp = server->ProcessQueryDirect(query, db).value();
    benchmark::DoNotOptimize(resp);
  }

  int64_t total_db_bytes = static_cast<int64_t>(num_entries * entry_size);
  state.SetBytesProcessed(state.iterations() * total_db_bytes);
}

BENCHMARK(BM_ScalingDbSize_OptimalN0)
    ->Arg(64)
    ->Arg(256)
    ->Arg(1024)
    ->Arg(2048)
    ->Arg(4096);

// ============================================================================
// NTT Microbenchmarks: Baseline RLWE vs Fast Vectorized NTT
// ============================================================================

void BM_Ntt_BaselineForward(benchmark::State& state) {
  int log_n = state.range(0);
  int n = 1 << log_n;
  auto context = rlwe::RnsContext<ModularInt>::CreateForBfv(
                     log_n, {1152921504606830593ULL}, {}, 65537)
                     .value();
  const auto* prime_mod = context.MainPrimeModuli()[0];
  const auto* ntt_params = prime_mod->NttParams();
  const auto* mod_params = prime_mod->ModParams();

  std::vector<ModularInt> coeffs(n, ModularInt::ImportZero(mod_params));
  for (int i = 0; i < n; ++i) {
    coeffs[i] =
        ModularInt::ImportInt(i % mod_params->modulus, mod_params).value();
  }

  for (auto _ : state) {
    auto tmp = coeffs;
    benchmark::DoNotOptimize(
        rlwe::ForwardNumberTheoreticTransform(tmp, *ntt_params, *mod_params));
  }
  state.SetBytesProcessed(state.iterations() * n * sizeof(uint64_t));
}
BENCHMARK(BM_Ntt_BaselineForward)->Arg(10)->Arg(11)->Arg(12);

void BM_Ntt_FastForward(benchmark::State& state) {
  int log_n = state.range(0);
  int n = 1 << log_n;
  auto context = rlwe::RnsContext<ModularInt>::CreateForBfv(
                     log_n, {1152921504606830593ULL}, {}, 65537)
                     .value();
  const auto* prime_mod = context.MainPrimeModuli()[0];
  const auto* ntt_params = prime_mod->NttParams();
  const auto* mod_params = prime_mod->ModParams();

  std::vector<ModularInt> coeffs(n, ModularInt::ImportZero(mod_params));
  for (int i = 0; i < n; ++i) {
    coeffs[i] =
        ModularInt::ImportInt(i % mod_params->modulus, mod_params).value();
  }

  for (auto _ : state) {
    auto tmp = coeffs;
    benchmark::DoNotOptimize(FastForwardNtt(tmp, *ntt_params, *mod_params));
  }
  state.SetBytesProcessed(state.iterations() * n * sizeof(uint64_t));
}
BENCHMARK(BM_Ntt_FastForward)->Arg(10)->Arg(11)->Arg(12);

void BM_Ntt_BaselineInverse(benchmark::State& state) {
  int log_n = state.range(0);
  int n = 1 << log_n;
  auto context = rlwe::RnsContext<ModularInt>::CreateForBfv(
                     log_n, {1152921504606830593ULL}, {}, 65537)
                     .value();
  const auto* prime_mod = context.MainPrimeModuli()[0];
  const auto* ntt_params = prime_mod->NttParams();
  const auto* mod_params = prime_mod->ModParams();

  std::vector<ModularInt> coeffs(n, ModularInt::ImportZero(mod_params));
  for (int i = 0; i < n; ++i) {
    coeffs[i] =
        ModularInt::ImportInt(i % mod_params->modulus, mod_params).value();
  }

  for (auto _ : state) {
    auto tmp = coeffs;
    benchmark::DoNotOptimize(
        rlwe::InverseNumberTheoreticTransform(tmp, *ntt_params, *mod_params));
  }
  state.SetBytesProcessed(state.iterations() * n * sizeof(uint64_t));
}
BENCHMARK(BM_Ntt_BaselineInverse)->Arg(10)->Arg(11)->Arg(12);

void BM_Ntt_FastInverse(benchmark::State& state) {
  int log_n = state.range(0);
  int n = 1 << log_n;
  auto context = rlwe::RnsContext<ModularInt>::CreateForBfv(
                     log_n, {1152921504606830593ULL}, {}, 65537)
                     .value();
  const auto* prime_mod = context.MainPrimeModuli()[0];
  const auto* ntt_params = prime_mod->NttParams();
  const auto* mod_params = prime_mod->ModParams();

  std::vector<ModularInt> coeffs(n, ModularInt::ImportZero(mod_params));
  for (int i = 0; i < n; ++i) {
    coeffs[i] =
        ModularInt::ImportInt(i % mod_params->modulus, mod_params).value();
  }

  for (auto _ : state) {
    auto tmp = coeffs;
    benchmark::DoNotOptimize(FastInverseNtt(tmp, *ntt_params, *mod_params));
  }
  state.SetBytesProcessed(state.iterations() * n * sizeof(uint64_t));
}
BENCHMARK(BM_Ntt_FastInverse)->Arg(10)->Arg(11)->Arg(12);

}  // namespace
}  // namespace oram_okvs::pir
