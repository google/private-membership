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

#include "research/oram_okvs/pir/fast_ntt.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <random>
#include <vector>

#include <gtest/gtest.h>
#include "absl/status/status.h"
#include "shell_encryption/dft_transformations.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"

namespace oram_okvs::pir {
namespace {

using ModularInt = rlwe::MontgomeryInt<uint64_t>;
using PrimeModulus = rlwe::PrimeModulus<ModularInt>;
using RnsContext = rlwe::RnsContext<ModularInt>;
using RnsPolynomial = rlwe::RnsPolynomial<ModularInt>;

std::vector<ModularInt> GenerateRandomCoeffs(
    size_t n, const ModularInt::Params* mod_params, uint64_t seed = 42) {
  std::mt19937_64 rng(seed);
  std::vector<ModularInt> coeffs;
  coeffs.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    uint64_t val = rng() % mod_params->modulus;
    coeffs.push_back(ModularInt::ImportInt(val, mod_params).value());
  }
  return coeffs;
}

TEST(FastNttTest, MatchesRlweBaselineForwardNtt) {
  for (int log_n : {6, 8, 10, 11}) {
    int n = 1 << log_n;
    auto context_or = RnsContext::CreateForBfv(
        log_n, /*qs=*/{1152921504606830593ULL}, /*ps=*/{}, /*t=*/65537);
    ASSERT_TRUE(context_or.ok());
    auto context = std::move(*context_or);

    const auto* prime_mod = context.MainPrimeModuli()[0];
    const auto* ntt_params = prime_mod->NttParams();
    const auto* mod_params = prime_mod->ModParams();

    auto coeffs_baseline = GenerateRandomCoeffs(n, mod_params, 12345);
    auto coeffs_fast = coeffs_baseline;

    ASSERT_TRUE(rlwe::ForwardNumberTheoreticTransform(coeffs_baseline,
                                                      *ntt_params, *mod_params)
                    .ok());
    ASSERT_TRUE(FastForwardNtt(coeffs_fast, *ntt_params, *mod_params).ok());

    for (size_t i = 0; i < n; ++i) {
      EXPECT_EQ(coeffs_fast[i].ExportInt(mod_params),
                coeffs_baseline[i].ExportInt(mod_params))
          << "Mismatch at index " << i << " for log_n = " << log_n;
    }
  }
}

TEST(FastNttTest, MatchesRlweBaselineInverseNtt) {
  for (int log_n : {6, 8, 10, 11}) {
    int n = 1 << log_n;
    auto context_or = RnsContext::CreateForBfv(
        log_n, /*qs=*/{1152921504606830593ULL}, /*ps=*/{}, /*t=*/65537);
    ASSERT_TRUE(context_or.ok());
    auto context = std::move(*context_or);

    const auto* prime_mod = context.MainPrimeModuli()[0];
    const auto* ntt_params = prime_mod->NttParams();
    const auto* mod_params = prime_mod->ModParams();

    auto coeffs_baseline = GenerateRandomCoeffs(n, mod_params, 54321);
    auto coeffs_fast = coeffs_baseline;

    ASSERT_TRUE(rlwe::InverseNumberTheoreticTransform(coeffs_baseline,
                                                      *ntt_params, *mod_params)
                    .ok());
    ASSERT_TRUE(FastInverseNtt(coeffs_fast, *ntt_params, *mod_params).ok());

    for (size_t i = 0; i < n; ++i) {
      EXPECT_EQ(coeffs_fast[i].ExportInt(mod_params),
                coeffs_baseline[i].ExportInt(mod_params))
          << "Mismatch at index " << i << " for log_n = " << log_n;
    }
  }
}

TEST(FastNttTest, ForwardAndInverseRoundtrip) {
  for (int log_n : {6, 8, 10, 11, 12}) {
    int n = 1 << log_n;
    auto context_or = RnsContext::CreateForBfv(
        log_n, /*qs=*/{1152921504606830593ULL, 1152921504606584833ULL},
        /*ps=*/{}, /*t=*/65537);
    ASSERT_TRUE(context_or.ok());
    auto context = std::move(*context_or);

    for (const auto* prime_mod : context.MainPrimeModuli()) {
      const auto* ntt_params = prime_mod->NttParams();
      const auto* mod_params = prime_mod->ModParams();

      auto orig_coeffs = GenerateRandomCoeffs(n, mod_params, 9999 + log_n);
      auto coeffs = orig_coeffs;

      ASSERT_TRUE(FastForwardNtt(coeffs, *ntt_params, *mod_params).ok());
      ASSERT_TRUE(FastInverseNtt(coeffs, *ntt_params, *mod_params).ok());

      for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(coeffs[i].ExportInt(mod_params),
                  orig_coeffs[i].ExportInt(mod_params))
            << "Mismatch at index " << i << " for log_n = " << log_n;
      }
    }
  }
}

TEST(FastNttTest, RnsPolynomialConvertToNttFormAndBack) {
  int log_n = 10;
  auto context_or = RnsContext::CreateForBfv(
      log_n, /*qs=*/{1152921504606830593ULL, 1152921504606584833ULL},
      /*ps=*/{}, /*t=*/65537);
  ASSERT_TRUE(context_or.ok());
  auto context = std::move(*context_or);
  auto moduli = context.MainPrimeModuli();

  std::vector<std::vector<ModularInt>> rns_coeffs;
  for (const auto* prime_mod : moduli) {
    rns_coeffs.push_back(
        GenerateRandomCoeffs(1 << log_n, prime_mod->ModParams(), 777));
  }

  auto orig_poly_or =
      RnsPolynomial::Create(std::vector<std::vector<ModularInt>>(rns_coeffs),
                            /*is_ntt=*/false);
  ASSERT_TRUE(orig_poly_or.ok());
  RnsPolynomial orig_poly = std::move(*orig_poly_or);

  RnsPolynomial poly = orig_poly.Clone();
  ASSERT_FALSE(poly.IsNttForm());

  ASSERT_TRUE(FastConvertToNttForm(poly, absl::MakeSpan(moduli)).ok());
  EXPECT_TRUE(poly.IsNttForm());

  ASSERT_TRUE(FastConvertToCoeffForm(poly, absl::MakeSpan(moduli)).ok());
  EXPECT_FALSE(poly.IsNttForm());

  EXPECT_EQ(poly, orig_poly);
}

}  // namespace
}  // namespace oram_okvs::pir
