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
#include <string>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/numeric/int128.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "research/oram_okvs/pir/fast_ntt.h"
#include "research/oram_okvs/pir/onion_pir_client.h"
#include "research/oram_okvs/pir/onion_pir_params.h"
#include "research/oram_okvs/pir/onion_pir_server.h"
#include "research/oram_okvs/pir/onion_pir_types.h"
#include "research/oram_okvs/pir/pir.h"
#include "research/oram_okvs/pir/pir.pb.h"
#include "shell_encryption/constants.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/prng/hkdf_prng.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/rns/rns_bfv_ciphertext.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/statusor.h"
#include "shell_encryption/testing/status_testing.h"

namespace oram_okvs::pir {
namespace {

// Helper to generate deterministic test data entries.
std::vector<std::vector<uint8_t>> GenerateTestDatabase(size_t num_entries,
                                                       size_t entry_size) {
  std::vector<std::vector<uint8_t>> db(num_entries);
  for (size_t i = 0; i < num_entries; ++i) {
    db[i].resize(entry_size);
    for (size_t j = 0; j < entry_size; ++j) {
      db[i][j] = static_cast<uint8_t>((i * 37 + j * 17 + 11) & 0xFF);
    }
  }
  return db;
}

TEST(OnionPirTypesTest, ModularInverseCorrectness) {
  rlwe::Uint64 m = 8193;  // (1 << 13) + 1
  for (size_t h = 0; h < 10; ++h) {
    rlwe::Uint64 w = static_cast<rlwe::Uint64>(1) << h;
    rlwe::Uint64 inv_w = ModularInverse(w, m);
    EXPECT_EQ((inv_w * w) % m, 1ULL);
  }
}

TEST(OnionPirTypesTest, BitReverseCorrectness) {
  // 3-bit reversal
  EXPECT_EQ(BitReverse(0b000, 3), 0b000);
  EXPECT_EQ(BitReverse(0b001, 3), 0b100);
  EXPECT_EQ(BitReverse(0b010, 3), 0b010);
  EXPECT_EQ(BitReverse(0b011, 3), 0b110);
  EXPECT_EQ(BitReverse(0b100, 3), 0b001);
  EXPECT_EQ(BitReverse(0b101, 3), 0b101);
  EXPECT_EQ(BitReverse(0b110, 3), 0b011);
  EXPECT_EQ(BitReverse(0b111, 3), 0b111);
}

TEST(OnionPirTypesTest, BitPackAndUnpackRoundTripAcrossBitWidths) {
  std::vector<size_t> test_bit_widths = {1,  3,  5,  8,  9,  14,
                                         16, 17, 18, 24, 27, 32};
  for (size_t num_bits : test_bit_widths) {
    uint64_t mask = (num_bits >= 64) ? ~0ULL : ((1ULL << num_bits) - 1ULL);
    size_t num_coeffs = 100;
    std::vector<uint64_t> coeffs(num_coeffs);
    for (size_t i = 0; i < num_coeffs; ++i) {
      coeffs[i] = (i * 1234567891ULL + 7ULL) & mask;
    }

    std::string packed = BitPack(coeffs, num_bits);
    size_t expected_bytes = (num_coeffs * num_bits + 7) / 8;
    EXPECT_EQ(packed.size(), expected_bytes)
        << "Failed size for num_bits = " << num_bits;

    ASSERT_OK_AND_ASSIGN(std::vector<uint64_t> unpacked,
                         BitUnpack(packed, num_coeffs, num_bits));
    EXPECT_EQ(unpacked, coeffs)
        << "Failed roundtrip for num_bits = " << num_bits;
  }
}

TEST(OnionPirTypesTest, GadgetDecompositionReconstructsOriginal) {
  size_t poly_degree = 1024;
  size_t log_n = 10;
  rlwe::Uint64 q = rlwe::kModulus59;
  rlwe::Uint64 t = 8193;

  ASSERT_OK_AND_ASSIGN(auto rns_context,
                       RnsContext::CreateForBfv(log_n, {q}, /*ps=*/{}, t));
  auto main_moduli = rns_context.MainPrimeModuli();

  int level = main_moduli.size() - 1;
  ASSERT_OK_AND_ASSIGN(auto q_hats,
                       rns_context.MainPrimeModulusComplements(level));
  ASSERT_OK_AND_ASSIGN(auto q_hat_invs,
                       rns_context.MainPrimeModulusCrtFactors(level));

  size_t log_b = 6;
  std::vector<size_t> log_bs = {log_b};
  ASSERT_OK_AND_ASSIGN(auto gadget_val,
                       rlwe::RnsGadget<ModularInt>::Create(
                           log_n, log_bs, q_hats, q_hat_invs, main_moduli));
  auto gadget = std::make_unique<const rlwe::RnsGadget<ModularInt>>(
      std::move(gadget_val));
  size_t ell = gadget->Dimension();

  std::string prng_seed;
  ASSERT_OK_AND_ASSIGN(prng_seed, rlwe::HkdfPrng::GenerateSeed());
  std::unique_ptr<rlwe::HkdfPrng> prng;
  ASSERT_OK_AND_ASSIGN(prng, rlwe::HkdfPrng::Create(prng_seed));

  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial poly_ntt,
      RnsPolynomial::SampleUniform(log_n, prng.get(), main_moduli));

  uint64_t base_b = static_cast<uint64_t>(1) << log_b;

  std::vector<RnsPolynomial> decomposed_ntt;
  ASSERT_OK_AND_ASSIGN(decomposed_ntt,
                       GadgetDecompose(gadget.get(), poly_ntt, main_moduli));
  EXPECT_EQ(decomposed_ntt.size(), ell);

  // Reconstruct: sum_{k=0}^{ell-1} decomposed[k] * B^k
  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial reconstructed,
      RnsPolynomial::CreateZero(log_n, main_moduli, /*is_ntt=*/true));
  uint64_t g_k = 1;
  for (size_t k = 0; k < ell; ++k) {
    ASSERT_OK_AND_ASSIGN(RnsPolynomial term,
                         decomposed_ntt[k].Mul(g_k % q, main_moduli));
    ASSERT_OK(reconstructed.AddInPlace(term, main_moduli));
    g_k = absl::Uint128Low64((static_cast<absl::uint128>(g_k) * base_b) % q);
  }

  // Compare reconstructed against original in coefficient form
  RnsPolynomial orig_coeff = poly_ntt;
  ASSERT_OK(FastConvertToCoeffForm(orig_coeff, main_moduli));
  RnsPolynomial recon_coeff = reconstructed;
  ASSERT_OK(FastConvertToCoeffForm(recon_coeff, main_moduli));

  const auto& orig_coeffs = orig_coeff.Coeffs()[0];
  const auto& recon_coeffs = recon_coeff.Coeffs()[0];
  const ModularIntParams* mod_params = main_moduli[0]->ModParams();

  for (size_t j = 0; j < poly_degree; ++j) {
    EXPECT_EQ(recon_coeffs[j].ExportInt(mod_params),
              orig_coeffs[j].ExportInt(mod_params))
        << "Mismatch at coefficient " << j;
  }
}

TEST(OnionPirTest, ExternalProductWithRgsw) {
  OnionPirParams params = OnionPirParams::CreateFastForTesting(16, 64);
  params.n0 = 8;
  ASSERT_OK_AND_ASSIGN(auto client, OnionPirClient::Create(params));
  ASSERT_OK_AND_ASSIGN(auto server, OnionPirServer::Create(params));

  ASSERT_OK_AND_ASSIGN(auto eval_keys, client->GetSerializedEvaluationKeys());
  ASSERT_OK(server->SetEvaluationKeys(eval_keys));

  // Client creates packed query for index 8 (m0 = 0, m1 = 1)
  ASSERT_OK_AND_ASSIGN(auto query, client->CreatePackedQueryCiphertext(8));
  ASSERT_OK_AND_ASSIGN(auto traced_query, server->HomomorphicTrace(query));
  ASSERT_OK_AND_ASSIGN(auto leaves, server->ExpandBfv(traced_query));

  auto main_moduli = client->main_moduli();
  const ModularIntParams* mod_params = client->mod_params();
  rlwe::Uint64 q = params.primary_ciphertext_modulus();

  // Verify GadgetDecompose reconstructs leaves[8] c0
  ASSERT_OK_AND_ASSIGN(auto l8_c0, leaves[8].Component(0));
  ASSERT_OK_AND_ASSIGN(auto l8_c1, leaves[8].Component(1));

  ASSERT_OK_AND_ASSIGN(auto d0_ntt,
                       GadgetDecompose(client->gadget(), l8_c0, main_moduli));

  ASSERT_OK_AND_ASSIGN(
      RnsPolynomial recon0,
      RnsPolynomial::CreateZero(client->rns_context()->LogN(), main_moduli,
                                /*is_ntt=*/true));
  for (size_t k = 0; k < client->gadget()->Dimension(); ++k) {
    RnsPolynomial term = d0_ntt[k];
    ASSERT_OK(term.MulInPlace(client->gadget()->Component(k), main_moduli));
    ASSERT_OK(recon0.AddInPlace(term, main_moduli));
  }
  RnsPolynomial c0_inv = l8_c0;
  ASSERT_OK(FastConvertToCoeffForm(c0_inv, main_moduli));
  RnsPolynomial recon0_inv = recon0;
  ASSERT_OK(FastConvertToCoeffForm(recon0_inv, main_moduli));

  for (size_t m_idx = 0; m_idx < main_moduli.size(); ++m_idx) {
    const ModularIntParams* m_params = main_moduli[m_idx]->ModParams();
    for (size_t j = 0; j < params.poly_degree; ++j) {
      EXPECT_EQ(recon0_inv.Coeffs()[m_idx][j].ExportInt(m_params),
                c0_inv.Coeffs()[m_idx][j].ExportInt(m_params));
    }
  }

  // ExternalProduct(RGSW(s), leaves[8]) should decrypt to s * leaves[8]
  ASSERT_OK_AND_ASSIGN(
      auto bot_row,
      ExternalProduct(client->rgsw_secret_key(), leaves[8], main_moduli,
                      client->gadget(), client->rns_context(),
                      client->error_params()));

  ASSERT_OK_AND_ASSIGN(auto bot_c0, bot_row.Component(0));
  ASSERT_OK_AND_ASSIGN(auto bot_c1, bot_row.Component(1));

  ASSERT_OK_AND_ASSIGN(RnsPolynomial bot_c1_s,
                       bot_c1.Mul(client->secret_key().Key(), main_moduli));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial bot_noisy,
                       bot_c0.Add(bot_c1_s, main_moduli));
  RnsPolynomial bot_noisy_coeff = bot_noisy;
  ASSERT_OK(FastConvertToCoeffForm(bot_noisy_coeff, main_moduli));
  const auto& bot_coeffs = bot_noisy_coeff.Coeffs()[0];

  ASSERT_OK_AND_ASSIGN(RnsPolynomial l8_c1_s,
                       l8_c1.Mul(client->secret_key().Key(), main_moduli));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial l8_noisy, l8_c0.Add(l8_c1_s, main_moduli));
  ASSERT_OK_AND_ASSIGN(RnsPolynomial expected_bot,
                       l8_noisy.Mul(client->secret_key().Key(), main_moduli));
  RnsPolynomial exp_noisy_coeff = expected_bot;
  ASSERT_OK(FastConvertToCoeffForm(exp_noisy_coeff, main_moduli));
  const auto& exp_coeffs = exp_noisy_coeff.Coeffs()[0];

  // Check that the error between bot_noisy and leaves[8] * s is within
  // gadget noise bound
  double base_b = static_cast<double>(static_cast<uint64_t>(1)
                                      << params.conversion_log_base());
  double expected_std_dev =
      std::sqrt(2.0 * params.conversion_ell() * params.poly_degree *
                (base_b * base_b / 12.0) * params.variance);
  uint64_t noise_bound = static_cast<uint64_t>(10.0 * expected_std_dev);

  for (size_t j = 0; j < params.poly_degree; ++j) {
    rlwe::Uint64 diff = (bot_coeffs[j].ExportInt(mod_params) + q -
                         exp_coeffs[j].ExportInt(mod_params)) %
                        q;
    if (diff > q / 2) diff = q - diff;
    EXPECT_LT(diff, noise_bound);
  }
}

TEST(OnionPirTest, EndToEndSingleDimensionRetrieval) {
  OnionPirParams params = OnionPirParams::CreateFastForTesting(8, 64);
  params.n0 = 8;  // Single dimension (d = 1)
  EXPECT_EQ(params.num_dimensions(), 1);

  std::unique_ptr<OnionPirClient> client;
  ASSERT_OK_AND_ASSIGN(client, OnionPirClient::Create(params));
  std::unique_ptr<OnionPirServer> server;
  ASSERT_OK_AND_ASSIGN(server, OnionPirServer::Create(params));

  SerializedEvaluationKeys eval_keys;
  ASSERT_OK_AND_ASSIGN(eval_keys, client->GetSerializedEvaluationKeys());
  ASSERT_OK(server->SetEvaluationKeys(eval_keys));

  auto db = GenerateTestDatabase(params.num_entries, params.entry_size_bytes);
  ASSERT_OK(server->SetDatabase(db));

  for (size_t idx = 0; idx < params.num_entries; ++idx) {
    SerializedPirQuery query;
    ASSERT_OK_AND_ASSIGN(query, client->CreateQuery(idx));
    auto resp_status = server->ProcessQuery(query);
    ASSERT_TRUE(resp_status.ok())
        << "ProcessQuery failed: " << resp_status.status();
    SerializedPirResponse response = std::move(resp_status).value();
    std::vector<uint8_t> retrieved;
    ASSERT_OK_AND_ASSIGN(
        retrieved, client->DecodeResponse(response, params.entry_size_bytes));
    EXPECT_EQ(retrieved, db[idx])
        << "Failed to retrieve entry at index " << idx;
  }
}

TEST(OnionPirTest, EndToEndMultiDimension2D) {
  OnionPirParams params = OnionPirParams::CreateFastForTesting(16, 64);
  params.n0 = 8;  // 2 dimensions (d = 2, 8 x 2)
  EXPECT_EQ(params.num_dimensions(), 2);

  std::unique_ptr<OnionPirClient> client;
  ASSERT_OK_AND_ASSIGN(client, OnionPirClient::Create(params));
  std::unique_ptr<OnionPirServer> server;
  ASSERT_OK_AND_ASSIGN(server, OnionPirServer::Create(params));

  SerializedEvaluationKeys eval_keys;
  ASSERT_OK_AND_ASSIGN(eval_keys, client->GetSerializedEvaluationKeys());
  ASSERT_OK(server->SetEvaluationKeys(eval_keys));

  auto db = GenerateTestDatabase(params.num_entries, params.entry_size_bytes);
  ASSERT_OK(server->SetDatabase(db));

  for (size_t idx = 0; idx < params.num_entries; ++idx) {
    SerializedPirQuery query;
    ASSERT_OK_AND_ASSIGN(query, client->CreateQuery(idx));
    SerializedPirResponse response;
    ASSERT_OK_AND_ASSIGN(response, server->ProcessQuery(query));
    std::vector<uint8_t> retrieved;
    ASSERT_OK_AND_ASSIGN(
        retrieved, client->DecodeResponse(response, params.entry_size_bytes));
    EXPECT_EQ(retrieved, db[idx])
        << "Failed to retrieve entry at index " << idx;
  }
}

TEST(OnionPirTest, EndToEndMultiDimension3DAnd4D) {
  OnionPirParams params = OnionPirParams::CreateFastForTesting(32, 64);
  params.n0 = 8;  // 3 dimensions (d = 3, 8 x 2 x 2 = 32)
  EXPECT_EQ(params.num_dimensions(), 3);

  std::unique_ptr<OnionPirClient> client;
  ASSERT_OK_AND_ASSIGN(client, OnionPirClient::Create(params));
  std::unique_ptr<OnionPirServer> server;
  ASSERT_OK_AND_ASSIGN(server, OnionPirServer::Create(params));

  SerializedEvaluationKeys eval_keys;
  ASSERT_OK_AND_ASSIGN(eval_keys, client->GetSerializedEvaluationKeys());
  ASSERT_OK(server->SetEvaluationKeys(eval_keys));

  auto db = GenerateTestDatabase(params.num_entries, params.entry_size_bytes);

  // Test spot checks across the database
  std::vector<size_t> test_indices = {0, 1, 7, 8, 15, 16, 23, 24, 31};
  for (size_t idx : test_indices) {
    SerializedPirQuery query;
    ASSERT_OK_AND_ASSIGN(query, client->CreateQuery(idx));
    // Test on-the-fly database evaluation
    SerializedPirResponse response;
    ASSERT_OK_AND_ASSIGN(response, server->ProcessQuery(query, db));
    std::vector<uint8_t> retrieved;
    ASSERT_OK_AND_ASSIGN(
        retrieved, client->DecodeResponse(response, params.entry_size_bytes));
    EXPECT_EQ(retrieved, db[idx])
        << "Failed on-the-fly retrieval at index " << idx;
  }
}

TEST(OnionPirTest, OnTheFlyNttMatchesPrecomputedDb) {
  OnionPirParams params = OnionPirParams::CreateFastForTesting(16, 128);
  params.n0 = 8;

  std::unique_ptr<OnionPirClient> client;
  ASSERT_OK_AND_ASSIGN(client, OnionPirClient::Create(params));
  std::unique_ptr<OnionPirServer> server;
  ASSERT_OK_AND_ASSIGN(server, OnionPirServer::Create(params));

  SerializedEvaluationKeys eval_keys;
  ASSERT_OK_AND_ASSIGN(eval_keys, client->GetSerializedEvaluationKeys());
  ASSERT_OK(server->SetEvaluationKeys(eval_keys));

  auto db = GenerateTestDatabase(params.num_entries, params.entry_size_bytes);
  ASSERT_OK(server->SetDatabase(db));

  for (size_t idx : {0, 3, 8, 15}) {
    SerializedPirQuery query;
    ASSERT_OK_AND_ASSIGN(query, client->CreateQuery(idx));

    // Precomputed NTT query
    SerializedPirResponse resp_precomputed;
    ASSERT_OK_AND_ASSIGN(resp_precomputed, server->ProcessQuery(query));
    std::vector<uint8_t> res_precomputed;
    ASSERT_OK_AND_ASSIGN(
        res_precomputed,
        client->DecodeResponse(resp_precomputed, params.entry_size_bytes));

    // On-the-fly NTT query
    SerializedPirResponse resp_onthefly;
    ASSERT_OK_AND_ASSIGN(resp_onthefly, server->ProcessQuery(query, db));
    std::vector<uint8_t> res_onthefly;
    ASSERT_OK_AND_ASSIGN(
        res_onthefly,
        client->DecodeResponse(resp_onthefly, params.entry_size_bytes));

    EXPECT_EQ(res_precomputed, db[idx])
        << "Precomputed failed at index " << idx;
    EXPECT_EQ(res_onthefly, db[idx]) << "Onthefly failed at index " << idx;
  }
}

TEST(OnionPirTest, MultiPolynomialLargeEntrySize) {
  // Entry size larger than poly_degree (poly_degree = 1024 bytes, entry = 2500
  // bytes -> 3 polynomials)
  OnionPirParams params = OnionPirParams::CreateFastForTesting(8, 2500);
  params.n0 = 8;
  EXPECT_EQ(params.num_polynomials_per_entry(), 3);

  std::unique_ptr<OnionPirClient> client;
  ASSERT_OK_AND_ASSIGN(client, OnionPirClient::Create(params));
  std::unique_ptr<OnionPirServer> server;
  ASSERT_OK_AND_ASSIGN(server, OnionPirServer::Create(params));

  SerializedEvaluationKeys eval_keys;
  ASSERT_OK_AND_ASSIGN(eval_keys, client->GetSerializedEvaluationKeys());
  ASSERT_OK(server->SetEvaluationKeys(eval_keys));

  auto db = GenerateTestDatabase(params.num_entries, params.entry_size_bytes);
  ASSERT_OK(server->SetDatabase(db));

  for (size_t idx : {0, 2, 5, 7}) {
    SerializedPirQuery query;
    ASSERT_OK_AND_ASSIGN(query, client->CreateQuery(idx));
    SerializedPirResponse response;
    ASSERT_OK_AND_ASSIGN(response, server->ProcessQuery(query, db));
    std::vector<uint8_t> retrieved;
    ASSERT_OK_AND_ASSIGN(
        retrieved, client->DecodeResponse(response, params.entry_size_bytes));
    EXPECT_EQ(retrieved, db[idx])
        << "Multi-polynomial retrieval failed at index " << idx;
  }
}

TEST(OnionPirTest, GenericPirInterfacePolymorphism) {
  OnionPirParams params = OnionPirParams::CreateFastForTesting(16, 64);
  params.n0 = 8;

  std::unique_ptr<OnionPirClient> client_impl;
  ASSERT_OK_AND_ASSIGN(client_impl, OnionPirClient::Create(params));
  std::unique_ptr<OnionPirServer> server_impl;
  ASSERT_OK_AND_ASSIGN(server_impl, OnionPirServer::Create(params));

  std::unique_ptr<PirClient> client = std::move(client_impl);
  std::unique_ptr<PirServer> server = std::move(server_impl);

  EXPECT_EQ(client->num_entries(), 16);
  EXPECT_EQ(server->num_entries(), 16);

  SerializedEvaluationKeys eval_keys;
  ASSERT_OK_AND_ASSIGN(eval_keys, client->GetSerializedEvaluationKeys());
  ASSERT_OK(server->SetEvaluationKeys(eval_keys));

  auto db = GenerateTestDatabase(params.num_entries, params.entry_size_bytes);
  ASSERT_OK(server->SetDatabase(db));

  size_t target_index = 5;
  SerializedPirQuery query;
  ASSERT_OK_AND_ASSIGN(query, client->CreateQuery(target_index));
  SerializedPirResponse response;
  ASSERT_OK_AND_ASSIGN(response, server->ProcessQuery(query));
  std::vector<uint8_t> retrieved;
  ASSERT_OK_AND_ASSIGN(
      retrieved, client->DecodeResponse(response, params.entry_size_bytes));

  EXPECT_EQ(retrieved, db[target_index]);
}

TEST(OnionPirTest, ModulusSwitchingReducesResponseSizeAndDecryptsCorrectly) {
  OnionPirParams params_scaled = OnionPirParams::CreateFastForTesting(16, 64);
  params_scaled.response_modulus = 262144ULL;

  auto db = GenerateTestDatabase(params_scaled.num_entries,
                                 params_scaled.entry_size_bytes);

  ASSERT_OK_AND_ASSIGN(auto client, OnionPirClient::Create(params_scaled));
  ASSERT_OK_AND_ASSIGN(auto server, OnionPirServer::Create(params_scaled));
  ASSERT_OK_AND_ASSIGN(auto keys, client->GetSerializedEvaluationKeys());
  ASSERT_OK(server->SetEvaluationKeys(keys));
  ASSERT_OK_AND_ASSIGN(auto query, client->CreateQuery(7));
  ASSERT_OK_AND_ASSIGN(auto resp, server->ProcessQuery(query, db));

  ASSERT_OK_AND_ASSIGN(auto dec, client->DecodeResponse(resp, 64));
  EXPECT_EQ(dec, db[7]);

  // Verify response format and bit-packed size (18 bits per coeff -> 1024 * 18
  // / 8 = 2304 bytes)
  EXPECT_GT(resp.scaled_ct_size(), 0);
  EXPECT_EQ(resp.scaled_ct(0).modulus(), 262144ULL);
  size_t expected_bytes_per_poly = (params_scaled.poly_degree * 18 + 7) / 8;
  EXPECT_EQ(resp.scaled_ct(0).c0().size(), expected_bytes_per_poly);
  EXPECT_EQ(resp.scaled_ct(0).c1().size(), expected_bytes_per_poly);
}

TEST(OnionPirTest, IndependenceHeuristicNoiseEstimation) {
  // Test FastForTesting parameters
  OnionPirParams fast_params = OnionPirParams::CreateFastForTesting(16, 64);
  NoiseEstimate fast_est = fast_params.EstimateNoise(/*confidence_factor=*/6.0);

  EXPECT_TRUE(fast_est.is_decryption_correct);
  EXPECT_GT(fast_est.noise_budget_bits, 0.0);
  EXPECT_GT(fast_est.fresh_noise_variance, 0.0);
  EXPECT_GT(fast_est.keyswitch_noise_variance, 0.0);
  EXPECT_GT(fast_est.unpacked_query_noise_variance,
            fast_est.fresh_noise_variance);
  EXPECT_GT(fast_est.dim0_dot_product_noise_variance, 0.0);
  EXPECT_GT(fast_est.final_noise_std_dev, 0.0);

  // Test 128BitSecure parameters
  OnionPirParams secure_params =
      OnionPirParams::Create128BitSecure(/*num_blocks=*/1024, /*block_size=*/64,
                                         /*n0=*/32);
  NoiseEstimate secure_est =
      secure_params.EstimateNoise(/*confidence_factor=*/6.0);

  EXPECT_TRUE(secure_est.is_decryption_correct);
  EXPECT_GT(secure_est.noise_budget_bits, 0.0);
  EXPECT_GT(secure_est.higher_dims_noise_variance, 0.0);

  // Test single dimension (d = 1) has zero higher dimensions noise
  OnionPirParams single_dim_params =
      OnionPirParams::CreateFastForTesting(8, 64);
  single_dim_params.n0 = 8;
  EXPECT_EQ(single_dim_params.num_dimensions(), 1);
  NoiseEstimate single_est = single_dim_params.EstimateNoise();
  EXPECT_EQ(single_est.higher_dims_noise_variance, 0.0);
  EXPECT_TRUE(single_est.is_decryption_correct);
}

TEST(OnionPirTest,
     HomomorphicTraceZeroesGarbageCoefficientsAndScalesCorrectly) {
  OnionPirParams params = OnionPirParams::CreateFastForTesting(16, 64);
  params.n0 = 8;
  ASSERT_OK_AND_ASSIGN(auto client, OnionPirClient::Create(params));
  ASSERT_OK_AND_ASSIGN(auto server, OnionPirServer::Create(params));

  SerializedEvaluationKeys eval_keys;
  ASSERT_OK_AND_ASSIGN(eval_keys, client->GetSerializedEvaluationKeys());
  ASSERT_OK(server->SetEvaluationKeys(eval_keys));

  // Generate a compressed query and deserialize it (creates zeros at non-stride
  // positions)
  ASSERT_OK_AND_ASSIGN(SerializedPirQuery query, client->CreateQuery(5));
  ASSERT_OK_AND_ASSIGN(auto reconstructed_ct,
                       server->DeserializeQueryCiphertext(query));

  // HomomorphicTrace should eliminate garbage from a * s at non-stride
  // positions
  ASSERT_OK_AND_ASSIGN(auto traced_ct,
                       server->HomomorphicTrace(reconstructed_ct));

  // Expand the traced query
  ASSERT_OK_AND_ASSIGN(auto leaves, server->ExpandBfv(traced_ct));
  EXPECT_EQ(leaves.size(), params.total_packed_values());
}

TEST(OnionPirTest, CompressedQueryOptimizationCutsSizeDrastically) {
  OnionPirParams params = OnionPirParams::CreateFastForTesting(16, 64);
  ASSERT_OK_AND_ASSIGN(auto client, OnionPirClient::Create(params));
  ASSERT_OK_AND_ASSIGN(auto server, OnionPirServer::Create(params));

  SerializedEvaluationKeys eval_keys;
  ASSERT_OK_AND_ASSIGN(eval_keys, client->GetSerializedEvaluationKeys());
  ASSERT_OK(server->SetEvaluationKeys(eval_keys));

  auto db = GenerateTestDatabase(params.num_entries, params.entry_size_bytes);
  ASSERT_OK(server->SetDatabase(db));

  // 1. Generate optimized query with compressed coefficients + PRG seed
  ASSERT_OK_AND_ASSIGN(SerializedPirQuery query, client->CreateQuery(5));
  EXPECT_FALSE(query.prng_seed().empty());
  EXPECT_EQ(query.prng_seed().size(), 64);
  EXPECT_EQ(query.num_coeffs(), params.query_w());
  EXPECT_FALSE(query.compressed_coeffs().empty());

  // 2. Generate 2-polynomial full ciphertext for size comparison
  ASSERT_OK_AND_ASSIGN(auto full_ct, client->CreatePackedQueryCiphertext(5));
  ASSERT_OK_AND_ASSIGN(auto full_serialized, full_ct.Serialize());

  // 3. Verify query size is dramatically reduced
  size_t compressed_query_bytes = query.ByteSizeLong();
  size_t full_query_bytes = full_serialized.ByteSizeLong();
  EXPECT_GT(full_query_bytes, 0);
  EXPECT_LT(compressed_query_bytes, full_query_bytes);
  // Compressed query should be a small fraction of the full ciphertext size
  EXPECT_LT(compressed_query_bytes, full_query_bytes / 10);

  // 4. Verify server correctly processes compressed query and decrypts
  ASSERT_OK_AND_ASSIGN(auto resp, server->ProcessQuery(query));
  ASSERT_OK_AND_ASSIGN(auto dec,
                       client->DecodeResponse(resp, params.entry_size_bytes));
  EXPECT_EQ(dec, db[5]);
}

TEST(OnionPirTest, ThreeIndependentGadgetDecompositions) {
  OnionPirParams params = OnionPirParams::CreateFastForTesting(16, 64);
  params.n0 = 8;  // 2D: 8 x 2
  // Explicitly distinct gadget bases for conversion, expansion, and trace
  params.conversion_gadget_log_base = 7;  // ell = 8
  params.expansion_gadget_log_base = 4;   // ell = 14
  params.trace_gadget_log_base = 3;       // ell = 18

  EXPECT_EQ(params.conversion_ell(), 8);
  EXPECT_EQ(params.expansion_ell(), 14);
  EXPECT_EQ(params.trace_ell(), 18);

  ASSERT_OK_AND_ASSIGN(auto client, OnionPirClient::Create(params));
  ASSERT_OK_AND_ASSIGN(auto server, OnionPirServer::Create(params));

  EXPECT_EQ(client->conversion_gadget()->Dimension(), 8);
  EXPECT_EQ(client->expansion_gadget()->Dimension(), 14);
  EXPECT_EQ(client->trace_gadget()->Dimension(), 18);

  EXPECT_EQ(server->conversion_gadget()->Dimension(), 8);
  EXPECT_EQ(server->expansion_gadget()->Dimension(), 14);
  EXPECT_EQ(server->trace_gadget()->Dimension(), 18);

  SerializedEvaluationKeys eval_keys;
  ASSERT_OK_AND_ASSIGN(eval_keys, client->GetSerializedEvaluationKeys());
  ASSERT_OK(server->SetEvaluationKeys(eval_keys));

  auto db = GenerateTestDatabase(params.num_entries, params.entry_size_bytes);
  ASSERT_OK(server->SetDatabase(db));

  for (size_t idx : {0, 3, 7, 8, 12, 15}) {
    SerializedPirQuery query;
    ASSERT_OK_AND_ASSIGN(query, client->CreateQuery(idx));
    SerializedPirResponse response;
    ASSERT_OK_AND_ASSIGN(response, server->ProcessQuery(query, db));
    std::vector<uint8_t> retrieved;
    ASSERT_OK_AND_ASSIGN(
        retrieved, client->DecodeResponse(response, params.entry_size_bytes));
    EXPECT_EQ(retrieved, db[idx])
        << "Failed retrieval with custom 3-gadget params at index " << idx;
  }
}

TEST(OnionPirTest, ExtraCoefficientsAndTraceReductionCorrectness) {
  // Test across multiple w_sent values from minimal w to N
  std::vector<size_t> test_w_sent = {0, 32, 64, 128, 256, 512, 1024};

  for (size_t w_sent : test_w_sent) {
    OnionPirParams params = OnionPirParams::CreateFastForTesting(16, 64);
    params.n0 = 8;  // 2D: 8 x 2 -> u = 8 + 8 = 16 -> w = 16
    params.num_query_coeffs = w_sent;

    ASSERT_OK_AND_ASSIGN(auto client, OnionPirClient::Create(params));
    ASSERT_OK_AND_ASSIGN(auto server, OnionPirServer::Create(params));

    SerializedEvaluationKeys eval_keys;
    ASSERT_OK_AND_ASSIGN(eval_keys, client->GetSerializedEvaluationKeys());
    ASSERT_OK(server->SetEvaluationKeys(eval_keys));

    auto db = GenerateTestDatabase(params.num_entries, params.entry_size_bytes);
    ASSERT_OK(server->SetDatabase(db));

    for (size_t idx : {0, 5, 11, 15}) {
      ASSERT_OK_AND_ASSIGN(SerializedPirQuery query, client->CreateQuery(idx));
      EXPECT_EQ(query.num_coeffs(), params.query_w());

      ASSERT_OK_AND_ASSIGN(SerializedPirResponse response,
                           server->ProcessQuery(query, db));
      ASSERT_OK_AND_ASSIGN(
          std::vector<uint8_t> retrieved,
          client->DecodeResponse(response, params.entry_size_bytes));
      EXPECT_EQ(retrieved, db[idx])
          << "Failed retrieval for w_sent = " << w_sent << " at index " << idx;
    }
  }
}

TEST(OnionPirTest, SetMaxQueryBytesCalculatesOptimalCoefficients) {
  OnionPirParams params =
      OnionPirParams::Create128BitSecure(1024, 64, /*n0=*/32);
  // Default w() = 64 (n0=32, d=6 -> u = 32 + 8*5 = 72 -> w = 128)
  size_t min_w = params.w();

  params.SetMaxQueryBytes(1024);
  EXPECT_GE(params.query_w(), min_w);
  EXPECT_LE(params.query_w(), params.poly_degree);

  // Bit-packed query with query_w() coefficients should fit in <= 1024 bytes
  ASSERT_OK_AND_ASSIGN(auto client, OnionPirClient::Create(params));
  ASSERT_OK_AND_ASSIGN(SerializedPirQuery query, client->CreateQuery(0));
  EXPECT_LE(query.ByteSizeLong(), 1024);
}

}  // namespace
}  // namespace oram_okvs::pir
