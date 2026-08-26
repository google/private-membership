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

#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_format.h"
#include "research/oram_okvs/okvs/gf128.h"
#include "research/oram_okvs/okvs/random_band_okvs.h"

ABSL_FLAG(size_t, block_size, 2048,
          "Payload block size B in bytes (default: 2048 B = 2 KB)");
ABSL_FLAG(int, trials, 10,
          "Number of trials per N for accurate timing and success rate");

namespace oram_okvs::okvs {
namespace {

constexpr size_t kBandWidth = 350;
constexpr double kEpsilon = 0.05;

struct BenchmarkResult {
  size_t N;
  size_t m_table_size;
  double encode_time_ms;
  double decode_time_us;
  double encode_throughput_mb_s;
  double decode_throughput_mb_s;
  double success_rate;
};

template <size_t B>
BenchmarkResult RunBenchmarkForN(size_t N, int num_trials) {
  std::mt19937_64 rng(1337 + N);

  // Generate N random keys and B-byte values
  std::vector<std::string> keys(N);
  std::vector<std::array<uint8_t, B>> values(N);

  for (size_t i = 0; i < N; ++i) {
    keys[i] = absl::StrCat("key_", i, "_", rng());
    for (size_t b = 0; b < B; ++b) {
      values[i][b] = static_cast<uint8_t>(rng() & 0xFF);
    }
  }

  int success_count = 0;
  double total_encode_ms = 0.0;
  size_t actual_table_size = 0;

  // 1. Benchmark Encode across trials
  for (int trial = 0; trial < num_trials; ++trial) {
    RandomBandOkvs<std::string, std::array<uint8_t, B>, kBandWidth> okvs(
        kEpsilon, /*seed=*/rng());

    auto t0 = std::chrono::steady_clock::now();
    auto status = okvs.Encode(keys, values);
    auto t1 = std::chrono::steady_clock::now();

    if (status.ok()) {
      ++success_count;
      total_encode_ms +=
          std::chrono::duration<double, std::milli>(t1 - t0).count();
      actual_table_size = okvs.TableSize();
    }
  }

  double avg_encode_ms =
      (success_count > 0) ? (total_encode_ms / success_count) : 0.0;
  double success_rate = static_cast<double>(success_count) / num_trials * 100.0;

  // 2. Benchmark Decode on a successful OKVS instance
  double avg_decode_us = 0.0;
  if (success_count > 0) {
    RandomBandOkvs<std::string, std::array<uint8_t, B>, kBandWidth> okvs(
        kEpsilon, /*seed=*/42);
    auto status = okvs.Encode(keys, values);
    if (status.ok()) {
      const size_t num_queries =
          std::max<size_t>(1000, std::min<size_t>(N * 5, 20000));
      auto t0 = std::chrono::steady_clock::now();
      for (size_t q = 0; q < num_queries; ++q) {
        size_t idx = q % N;
        auto val_or = okvs.Decode(keys[idx]);
        if (!val_or.ok()) {
          std::cerr << "Decode error: " << val_or.status() << "\n";
        }
      }
      auto t1 = std::chrono::steady_clock::now();
      double total_decode_us =
          std::chrono::duration<double, std::micro>(t1 - t0).count();
      avg_decode_us = total_decode_us / num_queries;
    }
  }

  double total_bytes_mb = (static_cast<double>(N * B)) / (1024.0 * 1024.0);
  double enc_throughput =
      (avg_encode_ms > 0) ? (total_bytes_mb / (avg_encode_ms / 1000.0)) : 0.0;
  double dec_throughput = (avg_decode_us > 0)
                              ? ((static_cast<double>(B) / (1024.0 * 1024.0)) /
                                 (avg_decode_us / 1e6))
                              : 0.0;

  return BenchmarkResult{
      .N = N,
      .m_table_size = actual_table_size,
      .encode_time_ms = avg_encode_ms,
      .decode_time_us = avg_decode_us,
      .encode_throughput_mb_s = enc_throughput,
      .decode_throughput_mb_s = dec_throughput,
      .success_rate = success_rate,
  };
}

}  // namespace
}  // namespace oram_okvs::okvs

int main(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);

  size_t B = absl::GetFlag(FLAGS_block_size);
  int trials = absl::GetFlag(FLAGS_trials);

  std::cout << "\n============================================================="
               "============================================\n";
  std::cout
      << "  RANDOM BAND OKVS (RB-OKVS / RIBBON) EMPIRICAL BENCHMARK SUITE\n";
  std::cout << absl::StrFormat(
      "  Parameters: Ribbon Width w = %zu, Expansion epsilon = %.2f (Overhead "
      "factor %.2fx)\n",
      oram_okvs::okvs::kBandWidth, oram_okvs::okvs::kEpsilon,
      1.0 + oram_okvs::okvs::kEpsilon);
  std::cout << absl::StrFormat(
      "  Block Size B = %zu B (%zu KB), Trials per N = %d\n", B, B / 1024,
      trials);
  std::cout << "  N Sweep: 2^9 (512) up to 2^16 (65536), multiplying by 2\n";
  std::cout << "==============================================================="
               "==========================================\n";
  std::cout << absl::StrFormat(
      "%-8s | %-8s | %-12s | %-16s | %-16s | %-16s | %-12s\n", "N (Exp)",
      "N (blks)", "Table Size m", "Encode Time (ms)", "Encode (MB/s)",
      "Decode Time (µs)", "Success Rate");
  std::cout << "---------------------------------------------------------------"
               "------------------------------------------\n";

  std::vector<size_t> exponents = {9, 10, 11, 12, 13, 14, 15, 16};

  for (size_t exp : exponents) {
    size_t N = 1ULL << exp;
    oram_okvs::okvs::BenchmarkResult res;
    if (B == 2048) {
      res = oram_okvs::okvs::RunBenchmarkForN<2048>(N, trials);
    } else if (B == 1024) {
      res = oram_okvs::okvs::RunBenchmarkForN<1024>(N, trials);
    } else {
      res = oram_okvs::okvs::RunBenchmarkForN<128>(N, trials);
    }

    std::string n_exp_str = absl::StrFormat("2^%zu", exp);
    std::cout << absl::StrFormat(
        "%-8s | %-8zu | %-12zu | %-16.3f | %-16.2f | %-16.3f | %-11.1f%%\n",
        n_exp_str, res.N, res.m_table_size, res.encode_time_ms,
        res.encode_throughput_mb_s, res.decode_time_us, res.success_rate);
  }

  std::cout << "==============================================================="
               "==========================================\n\n";

  return 0;
}
