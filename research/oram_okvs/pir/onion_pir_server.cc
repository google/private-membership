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

#include "research/oram_okvs/pir/onion_pir_server.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "absl/memory/memory.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/types/span.h"
#include "research/oram_okvs/pir/fast_ntt.h"
#include "research/oram_okvs/pir/onion_pir_params.h"
#include "research/oram_okvs/pir/onion_pir_types.h"
#include "research/oram_okvs/pir/pir.pb.h"
#include "shell_encryption/constants.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/prng/hkdf_prng.h"
#include "shell_encryption/rns/rns_bfv_ciphertext.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_gadget.h"
#include "shell_encryption/rns/rns_galois_key.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/status_macros.h"

namespace oram_okvs::pir {

OnionPirServer::OnionPirServer(
    OnionPirParams params, std::unique_ptr<const RnsContext> rns_context,
    std::unique_ptr<const RnsErrorParams> error_params,
    std::unique_ptr<const rlwe::RnsGadget<ModularInt>> conversion_gadget,
    std::unique_ptr<const rlwe::RnsGadget<ModularInt>> expansion_gadget,
    std::unique_ptr<const rlwe::RnsGadget<ModularInt>> trace_gadget,
    std::vector<RnsPolynomial> neg_monomials_ntt)
    : params_(std::move(params)),
      rns_context_(std::move(rns_context)),
      error_params_(std::move(error_params)),
      conversion_gadget_(std::move(conversion_gadget)),
      expansion_gadget_(std::move(expansion_gadget)),
      trace_gadget_(std::move(trace_gadget)),
      neg_monomials_ntt_(std::move(neg_monomials_ntt)) {}

absl::StatusOr<std::unique_ptr<OnionPirServer>> OnionPirServer::Create(
    OnionPirParams params) {
  if (params.poly_degree == 0 ||
      (params.poly_degree & (params.poly_degree - 1)) != 0) {
    return absl::InvalidArgumentError(
        "Polynomial degree must be a power of 2.");
  }
  size_t log_n = 0;
  while ((static_cast<size_t>(1) << log_n) < params.poly_degree) {
    ++log_n;
  }

  RLWE_ASSIGN_OR_RETURN(
      auto rns_context_val,
      RnsContext::CreateForBfv(log_n, params.ciphertext_moduli, /*ps=*/{},
                               params.plaintext_modulus()));
  auto rns_context =
      std::make_unique<const RnsContext>(std::move(rns_context_val));
  auto main_moduli = rns_context->MainPrimeModuli();

  RLWE_ASSIGN_OR_RETURN(
      auto error_params_val,
      RnsErrorParams::Create(log_n, main_moduli, /*aux_moduli=*/{},
                             params.log_t, sqrt(params.variance)));
  auto error_params =
      std::make_unique<const RnsErrorParams>(std::move(error_params_val));

  int level = main_moduli.size() - 1;
  RLWE_ASSIGN_OR_RETURN(auto q_hats,
                        rns_context->MainPrimeModulusComplements(level));
  RLWE_ASSIGN_OR_RETURN(auto q_hat_invs,
                        rns_context->MainPrimeModulusCrtFactors(level));

  std::vector<size_t> conv_log_bs = params.conversion_gadget_log_bs;
  if (conv_log_bs.empty()) {
    conv_log_bs.assign(main_moduli.size(), params.conversion_log_base());
  }
  RLWE_ASSIGN_OR_RETURN(
      auto conv_gadget_val,
      rlwe::RnsGadget<ModularInt>::Create(log_n, conv_log_bs, q_hats,
                                          q_hat_invs, main_moduli));
  auto conversion_gadget = std::make_unique<const rlwe::RnsGadget<ModularInt>>(
      std::move(conv_gadget_val));

  std::vector<size_t> exp_log_bs = params.expansion_gadget_log_bs;
  if (exp_log_bs.empty()) {
    exp_log_bs.assign(main_moduli.size(), params.expansion_log_base());
  }
  RLWE_ASSIGN_OR_RETURN(
      auto exp_gadget_val,
      rlwe::RnsGadget<ModularInt>::Create(log_n, exp_log_bs, q_hats, q_hat_invs,
                                          main_moduli));
  auto expansion_gadget = std::make_unique<const rlwe::RnsGadget<ModularInt>>(
      std::move(exp_gadget_val));

  std::vector<size_t> trace_log_bs = params.trace_gadget_log_bs;
  if (trace_log_bs.empty()) {
    trace_log_bs.assign(main_moduli.size(), params.trace_log_base());
  }
  RLWE_ASSIGN_OR_RETURN(
      auto trace_gadget_val,
      rlwe::RnsGadget<ModularInt>::Create(log_n, trace_log_bs, q_hats,
                                          q_hat_invs, main_moduli));
  auto trace_gadget = std::make_unique<const rlwe::RnsGadget<ModularInt>>(
      std::move(trace_gadget_val));

  // Precompute negative monomial polynomials -x^{n - k_shift} in NTT form for
  // each level L in [0, log_n - 1], where k_shift = 1 << L.
  std::vector<RnsPolynomial> neg_monomials_ntt;
  neg_monomials_ntt.reserve(log_n);

  for (size_t L = 0; L < log_n; ++L) {
    size_t k_shift = static_cast<size_t>(1) << L;
    size_t deg = params.poly_degree - k_shift;
    std::vector<std::vector<ModularInt>> rns_coeffs;
    rns_coeffs.reserve(main_moduli.size());
    for (const auto* prime_mod : main_moduli) {
      const auto* mod_params = prime_mod->ModParams();
      std::vector<ModularInt> coeffs(params.poly_degree,
                                     ModularInt::ImportZero(mod_params));
      RLWE_ASSIGN_OR_RETURN(
          coeffs[deg],
          ModularInt::ImportInt(prime_mod->Modulus() - 1, mod_params));
      rns_coeffs.push_back(std::move(coeffs));
    }
    RLWE_ASSIGN_OR_RETURN(
        RnsPolynomial poly_ntt,
        RnsPolynomial::Create(std::move(rns_coeffs), /*is_ntt=*/false));
    RLWE_RETURN_IF_ERROR(FastConvertToNttForm(poly_ntt, main_moduli));
    neg_monomials_ntt.push_back(std::move(poly_ntt));
  }

  return absl::WrapUnique(new OnionPirServer(
      std::move(params), std::move(rns_context), std::move(error_params),
      std::move(conversion_gadget), std::move(expansion_gadget),
      std::move(trace_gadget), std::move(neg_monomials_ntt)));
}

absl::Status OnionPirServer::SetEvaluationKeysDirect(
    std::vector<rlwe::RnsGaloisKey<ModularInt>> galois_keys,
    RgswCiphertext rgsw_secret_key) {
  galois_keys_ = std::move(galois_keys);
  rgsw_secret_key_ = std::move(rgsw_secret_key);
  keys_loaded_ = true;
  return absl::OkStatus();
}

absl::Status OnionPirServer::SetEvaluationKeys(
    const SerializedEvaluationKeys& keys) {
  auto main_moduli = rns_context_->MainPrimeModuli();
  size_t trace_steps = params_.query_trace_steps();
  std::vector<rlwe::RnsGaloisKey<ModularInt>> gkeys;
  gkeys.reserve(keys.galois_keys_size());
  for (int i = 0; i < keys.galois_keys_size(); ++i) {
    const auto* g_ptr = (static_cast<size_t>(i) < trace_steps)
                            ? trace_gadget_.get()
                            : expansion_gadget_.get();
    RLWE_ASSIGN_OR_RETURN(auto gkey,
                          rlwe::RnsGaloisKey<ModularInt>::Deserialize(
                              keys.galois_keys(i), g_ptr, main_moduli));
    gkeys.push_back(std::move(gkey));
  }

  RLWE_ASSIGN_OR_RETURN(
      auto rgsw_key,
      RgswCiphertext::Deserialize(keys.rgsw_secret_key(), main_moduli,
                                  error_params_.get(), rns_context_.get()));

  return SetEvaluationKeysDirect(std::move(gkeys), std::move(rgsw_key));
}

absl::Status OnionPirServer::SetDatabase(
    absl::Span<const std::vector<uint8_t>> db) {
  size_t num_poly_per_entry = params_.num_polynomials_per_entry();
  size_t bytes_per_poly = params_.bytes_per_polynomial();
  auto main_moduli = rns_context_->MainPrimeModuli();

  precomputed_db_ntt_.clear();
  precomputed_db_ntt_.resize(db.size());

  for (size_t i = 0; i < db.size(); ++i) {
    precomputed_db_ntt_[i].reserve(num_poly_per_entry);
    for (size_t p = 0; p < num_poly_per_entry; ++p) {
      size_t byte_offset = p * bytes_per_poly;
      absl::Span<const uint8_t> chunk_span;
      if (byte_offset < db[i].size()) {
        size_t len = std::min(db[i].size() - byte_offset, bytes_per_poly);
        chunk_span = absl::MakeConstSpan(db[i].data() + byte_offset, len);
      }
      RLWE_ASSIGN_OR_RETURN(RnsPolynomial poly_ntt,
                            EncodeBytesToPolynomialNtt(
                                chunk_span, params_.poly_degree, main_moduli));
      precomputed_db_ntt_[i].push_back(std::move(poly_ntt));
    }
  }

  db_loaded_ = true;
  return absl::OkStatus();
}

absl::StatusOr<RnsBfvCiphertext> OnionPirServer::HomomorphicTrace(
    RnsBfvCiphertext ct) const {
  if (!keys_loaded_) {
    return absl::FailedPreconditionError("Evaluation keys not loaded.");
  }
  size_t log_n = rns_context_->LogN();
  size_t trace_steps = params_.query_trace_steps();

  if (galois_keys_.size() < log_n) {
    return absl::InvalidArgumentError(
        absl::StrCat("Insufficient Galois keys: expected at least ", log_n,
                     ", got ", galois_keys_.size()));
  }

  for (size_t s = 0; s < trace_steps; ++s) {
    size_t k = static_cast<size_t>(1) << s;
    rlwe::Uint64 sub_power = (params_.poly_degree / k) + 1;
    const auto& gkey = galois_keys_[s];

    RLWE_ASSIGN_OR_RETURN(auto sub_ct,
                          ct.Substitute(static_cast<int>(sub_power)));
    RLWE_ASSIGN_OR_RETURN(auto c_prime, gkey.ApplyTo(sub_ct));
    RLWE_RETURN_IF_ERROR(ct.AddInPlace(c_prime));
  }
  return ct;
}

absl::StatusOr<std::vector<RnsBfvCiphertext>> OnionPirServer::ExpandBfv(
    const RnsBfvCiphertext& packed_query, size_t target_u) const {
  if (!keys_loaded_) {
    return absl::FailedPreconditionError("Evaluation keys not loaded.");
  }

  size_t log_n = rns_context_->LogN();
  size_t query_w = params_.query_w();
  size_t u = (target_u > 0 && target_u < params_.total_packed_values())
                 ? target_u
                 : params_.total_packed_values();
  size_t tree_height = params_.query_expansion_tree_height();
  size_t trace_steps = params_.query_trace_steps();
  auto moduli = rns_context_->MainPrimeModuli();

  if (galois_keys_.size() < log_n) {
    return absl::InvalidArgumentError(
        absl::StrCat("Insufficient Galois keys: expected at least ", log_n,
                     ", got ", galois_keys_.size()));
  }

  // Array representing the expansion tree: indices 1 to 2*query_w - 1.
  std::vector<RnsBfvCiphertext> c(
      2 * query_w, RnsBfvCiphertext::CreateZero(moduli, error_params_.get(),
                                                rns_context_.get()));
  c[1] = packed_query;

  for (size_t h = 0; h < tree_height; ++h) {
    size_t k = static_cast<size_t>(1) << h;
    rlwe::Uint64 sub_power = (query_w / k) + 1;
    size_t gkey_idx = trace_steps + h;
    const auto& gkey = galois_keys_[gkey_idx];
    const auto& x_inv_ntt = neg_monomials_ntt_[gkey_idx];

    for (size_t i = k; i < 2 * k; ++i) {
      // Tree pruning: if the subtree only leads to leaf indices >= u, skip it.
      size_t min_leaf = i * (query_w / k) - query_w;
      if (min_leaf >= u) {
        continue;
      }

      // Substitute x -> x^{sub_power}
      RLWE_ASSIGN_OR_RETURN(auto sub_ct,
                            c[i].Substitute(static_cast<int>(sub_power)));

      // Apply GaloisKey to relinearize back to s(x)
      RLWE_ASSIGN_OR_RETURN(auto c_prime, gkey.ApplyTo(sub_ct));

      // c[2i] = c[i] + c'
      RLWE_ASSIGN_OR_RETURN(c[2 * i], c[i] + c_prime);

      // c[2i + 1] = (c[i] - c') * x^{-k_shift}
      RLWE_ASSIGN_OR_RETURN(auto diff, c[i] - c_prime);
      RLWE_ASSIGN_OR_RETURN(c[2 * i + 1], diff.AbsorbSimple(x_inv_ntt));
    }
  }

  // Extract the first u leaf ciphertexts: c[query_w ... query_w + u - 1]
  std::vector<RnsBfvCiphertext> unpacked;
  unpacked.reserve(u);
  for (size_t idx = 0; idx < u; ++idx) {
    unpacked.push_back(std::move(c[query_w + idx]));
  }
  return unpacked;
}

absl::StatusOr<
    std::pair<std::vector<RnsBfvCiphertext>, std::vector<RgswCiphertext>>>
OnionPirServer::UnpackQuery(const RnsBfvCiphertext& packed_query,
                            size_t active_db_size) const {
  size_t n0 = params_.n0;
  size_t d = params_.num_dimensions();
  size_t ell = conversion_gadget_->Dimension();
  auto moduli = rns_context_->MainPrimeModuli();

  size_t active_n0 =
      (active_db_size > 0 && active_db_size <= n0) ? active_db_size : n0;
  size_t active_d = (active_db_size > 0 && active_db_size <= n0) ? 1 : d;
  size_t target_u = (active_d == 1) ? active_n0 : (n0 + (active_d - 1) * ell);

  // Apply Homomorphic Trace to zero out non-stride garbage coefficients
  RLWE_ASSIGN_OR_RETURN(auto traced_query, HomomorphicTrace(packed_query));

  // Expand the traced query into individual selection ciphertexts
  RLWE_ASSIGN_OR_RETURN(auto leaves, ExpandBfv(traced_query, target_u));

  // Initial dimension BFV query vector: first active_n0 leaves
  std::vector<RnsBfvCiphertext> c_bfv;
  c_bfv.reserve(active_n0);
  for (size_t k = 0; k < active_n0; ++k) {
    c_bfv.push_back(std::move(leaves[k]));
  }

  // Higher dimension RGSW selection ciphertexts
  std::vector<RgswCiphertext> c_rgsw;
  if (active_d > 1) {
    c_rgsw.reserve(active_d - 1);
    for (size_t i = 1; i < active_d; ++i) {
      size_t base = n0 + ell * (i - 1);
      std::vector<RnsBfvCiphertext> rgsw_rows;
      rgsw_rows.reserve(2 * ell);

      // Top ell rows: c[base + k]
      for (size_t k = 0; k < ell; ++k) {
        rgsw_rows.push_back(leaves[base + k]);
      }

      // Bottom ell rows: ExternalProduct(A, c[base + k])
      for (size_t k = 0; k < ell; ++k) {
        RLWE_ASSIGN_OR_RETURN(
            auto bot_row,
            ExternalProduct(rgsw_secret_key_, leaves[base + k], moduli,
                            conversion_gadget_.get(), rns_context_.get(),
                            error_params_.get()));
        rgsw_rows.push_back(std::move(bot_row));
      }

      c_rgsw.emplace_back(std::move(rgsw_rows));
    }
  }

  return std::make_pair(std::move(c_bfv), std::move(c_rgsw));
}

absl::StatusOr<std::vector<RnsBfvCiphertext>>
OnionPirServer::ProcessQueryDirect(
    const RnsBfvCiphertext& packed_query,
    absl::Span<const std::vector<uint8_t>> db) const {
  size_t active_db_size = !db.empty() ? db.size()
                                      : (db_loaded_ ? precomputed_db_ntt_.size()
                                                    : params_.num_entries);

  RLWE_ASSIGN_OR_RETURN(auto unpacked,
                        UnpackQuery(packed_query, active_db_size));
  return EvaluateDatabaseDirect(unpacked.first, unpacked.second, db);
}

absl::StatusOr<std::vector<RnsBfvCiphertext>>
OnionPirServer::EvaluateDatabaseDirect(
    absl::Span<const RnsBfvCiphertext> c_bfv,
    absl::Span<const RgswCiphertext> c_rgsw,
    absl::Span<const std::vector<uint8_t>> db) const {
  size_t active_db_size = !db.empty() ? db.size()
                                      : (db_loaded_ ? precomputed_db_ntt_.size()
                                                    : params_.num_entries);

  size_t n0 = c_bfv.size();
  size_t d = 1 + c_rgsw.size();
  size_t n_cols = (d > 1) ? (static_cast<size_t>(1) << (d - 1)) : 1;
  size_t num_poly_per_entry = params_.num_polynomials_per_entry();
  size_t bytes_per_poly = params_.bytes_per_polynomial();
  auto moduli = rns_context_->MainPrimeModuli();

  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial zero_poly,
      RnsPolynomial::CreateZero(rns_context_->LogN(), moduli, /*is_ntt=*/true));

  std::vector<RnsBfvCiphertext> response_cts;
  response_cts.reserve(num_poly_per_entry);

  for (size_t p = 0; p < num_poly_per_entry; ++p) {
    // Initial dimension evaluation (dimension 0 dot products)
    std::vector<RnsBfvCiphertext> db_prime;
    db_prime.reserve(n_cols);

    for (size_t j = 0; j < n_cols; ++j) {
      if (j * n0 >= active_db_size) {
        RLWE_ASSIGN_OR_RETURN(
            RnsPolynomial zero_c0,
            RnsPolynomial::CreateZero(rns_context_->LogN(), moduli,
                                      /*is_ntt=*/true));
        RLWE_ASSIGN_OR_RETURN(
            RnsPolynomial zero_c1,
            RnsPolynomial::CreateZero(rns_context_->LogN(), moduli,
                                      /*is_ntt=*/true));
        std::vector<RnsPolynomial> comps;
        comps.reserve(2);
        comps.push_back(std::move(zero_c0));
        comps.push_back(std::move(zero_c1));
        db_prime.push_back(RnsBfvCiphertext(
            std::move(comps), moduli, /*power_of_s=*/1,
            /*error=*/0.0, error_params_.get(), rns_context_.get()));
        continue;
      }

      std::optional<RnsBfvCiphertext> accum;
      for (size_t k = 0; k < n0; ++k) {
        size_t entry_idx = j * n0 + k;
        if (entry_idx >= active_db_size) {
          break;
        }
        RnsPolynomial entry_ntt = zero_poly;

        if (!db.empty()) {
          // On-the-fly NTT conversion of database entry
          if (entry_idx < db.size()) {
            size_t byte_offset = p * bytes_per_poly;
            absl::Span<const uint8_t> chunk_span;
            if (byte_offset < db[entry_idx].size()) {
              size_t len =
                  std::min(db[entry_idx].size() - byte_offset, bytes_per_poly);
              chunk_span =
                  absl::MakeConstSpan(db[entry_idx].data() + byte_offset, len);
            }
            RLWE_ASSIGN_OR_RETURN(entry_ntt,
                                  EncodeBytesToPolynomialNtt(
                                      chunk_span, params_.poly_degree, moduli));
          }
        } else if (db_loaded_) {
          // Use precomputed NTT database
          if (entry_idx < precomputed_db_ntt_.size() &&
              p < precomputed_db_ntt_[entry_idx].size()) {
            entry_ntt = precomputed_db_ntt_[entry_idx][p];
          }
        } else {
          return absl::FailedPreconditionError(
              "No database provided or loaded on server.");
        }

        RLWE_ASSIGN_OR_RETURN(auto term, c_bfv[k].AbsorbSimple(entry_ntt));
        if (!accum.has_value()) {
          accum = std::move(term);
        } else {
          RLWE_RETURN_IF_ERROR(accum->AddInPlace(term));
        }
      }
      if (!accum.has_value()) {
        RLWE_ASSIGN_OR_RETURN(
            RnsPolynomial zero_c0,
            RnsPolynomial::CreateZero(rns_context_->LogN(), moduli,
                                      /*is_ntt=*/true));
        RLWE_ASSIGN_OR_RETURN(
            RnsPolynomial zero_c1,
            RnsPolynomial::CreateZero(rns_context_->LogN(), moduli,
                                      /*is_ntt=*/true));
        std::vector<RnsPolynomial> comps;
        comps.reserve(2);
        comps.push_back(std::move(zero_c0));
        comps.push_back(std::move(zero_c1));
        accum = RnsBfvCiphertext(std::move(comps), moduli, /*power_of_s=*/1,
                                 /*error=*/0.0, error_params_.get(),
                                 rns_context_.get());
      }
      db_prime.push_back(std::move(*accum));
    }

    // Higher dimensions evaluation (Algorithm 4 lines 7-14)
    for (size_t i = 1; i < d; ++i) {
      size_t half = db_prime.size() / 2;
      for (size_t j = 0; j < half; ++j) {
        RLWE_ASSIGN_OR_RETURN(RnsBfvCiphertext db_diff,
                              db_prime[half + j] - db_prime[j]);
        RLWE_ASSIGN_OR_RETURN(
            RnsBfvCiphertext ext_res,
            ExternalProduct(c_rgsw[i - 1], db_diff, moduli,
                            conversion_gadget_.get(), rns_context_.get(),
                            error_params_.get()));
        RLWE_RETURN_IF_ERROR(db_prime[j].AddInPlace(ext_res));
      }
      db_prime.erase(db_prime.begin() + half, db_prime.end());
    }

    if (db_prime.empty()) {
      return absl::InternalError("DB evaluation produced empty result.");
    }
    auto final_ct = std::move(db_prime[0]);
    response_cts.push_back(std::move(final_ct));
  }

  return response_cts;
}

absl::StatusOr<RnsBfvCiphertext> OnionPirServer::DeserializeQueryCiphertext(
    const SerializedPirQuery& query) const {
  if (query.prng_seed().empty() || query.compressed_coeffs().empty()) {
    return absl::InvalidArgumentError(
        "Invalid or empty SerializedPirQuery: missing compressed_coeffs or "
        "prng_seed.");
  }

  size_t n = params_.poly_degree;
  size_t query_w =
      query.num_coeffs() > 0 ? query.num_coeffs() : params_.query_w();
  size_t stride = n / query_w;
  auto moduli = rns_context_->MainPrimeModuli();

  RLWE_ASSIGN_OR_RETURN(auto prng, rlwe::HkdfPrng::Create(query.prng_seed()));
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial a,
      RnsPolynomial::SampleUniform(rns_context_->LogN(), prng.get(), moduli));

  std::vector<std::vector<ModularInt>> c0_coeffs;
  c0_coeffs.reserve(moduli.size());

  absl::string_view packed_data = query.compressed_coeffs();
  size_t byte_offset = 0;

  for (size_t m_idx = 0; m_idx < moduli.size(); ++m_idx) {
    const auto* prime_mod = moduli[m_idx];
    const auto* mod_params = prime_mod->ModParams();
    uint64_t q_i = prime_mod->Modulus();
    size_t num_bits = NumBitsForModulus(q_i);
    size_t limb_bytes = (query_w * num_bits + 7) / 8;

    if (byte_offset + limb_bytes > packed_data.size()) {
      return absl::InvalidArgumentError(
          absl::StrCat("Malformed compressed_coeffs buffer: offset ",
                       byte_offset + limb_bytes, " > ", packed_data.size()));
    }

    absl::string_view limb_packed = packed_data.substr(byte_offset, limb_bytes);
    byte_offset += limb_bytes;

    RLWE_ASSIGN_OR_RETURN(std::vector<uint64_t> limb_vals,
                          BitUnpack(limb_packed, query_w, num_bits));

    std::vector<ModularInt> coeffs(n, ModularInt::ImportZero(mod_params));
    for (size_t l = 0; l < query_w; ++l) {
      RLWE_ASSIGN_OR_RETURN(coeffs[l * stride],
                            ModularInt::ImportInt(limb_vals[l], mod_params));
    }
    c0_coeffs.push_back(std::move(coeffs));
  }

  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial c0,
      RnsPolynomial::Create(std::move(c0_coeffs), /*is_ntt=*/false));
  RLWE_RETURN_IF_ERROR(FastConvertToNttForm(c0, moduli));

  std::vector<RnsPolynomial> comps;
  comps.reserve(2);
  comps.push_back(std::move(c0));
  comps.push_back(std::move(a));

  return RnsBfvCiphertext(std::move(comps), moduli, /*power_of_s=*/1,
                          /*error=*/0.0, error_params_.get(),
                          rns_context_.get());
}

absl::StatusOr<SerializedScaledCiphertext>
OnionPirServer::ScaleAndSwitchCiphertext(RnsBfvCiphertext ct) const {
  RLWE_RETURN_IF_ERROR(ct.ConvertToCoeffForm());

  RLWE_ASSIGN_OR_RETURN(auto c0_poly, ct.Component(0));
  RLWE_ASSIGN_OR_RETURN(auto c1_poly, ct.Component(1));

  auto main_moduli = rns_context_->MainPrimeModuli();
  size_t num_moduli = ct.NumModuli();
  uint64_t p =
      params_.response_modulus > 0 ? params_.response_modulus : 262144ULL;
  size_t n = params_.poly_degree;

  SerializedScaledCiphertext scaled_ct;
  scaled_ct.set_modulus(p);

  std::vector<uint64_t> c0_scaled(n);
  std::vector<uint64_t> c1_scaled(n);

  if (num_moduli == 1) {
    const auto* mod_params_q0 = main_moduli[0]->ModParams();
    uint64_t q0 = main_moduli[0]->Modulus();
    const auto& c0_coeffs = c0_poly.Coeffs()[0];
    const auto& c1_coeffs = c1_poly.Coeffs()[0];

    for (size_t i = 0; i < n; ++i) {
      uint64_t val0 = c0_coeffs[i].ExportInt(mod_params_q0);
      uint64_t val1 = c1_coeffs[i].ExportInt(mod_params_q0);
      uint64_t s0 = (static_cast<unsigned __int128>(val0) * p + (q0 / 2)) / q0;
      uint64_t s1 = (static_cast<unsigned __int128>(val1) * p + (q0 / 2)) / q0;
      c0_scaled[i] = s0 % p;
      c1_scaled[i] = s1 % p;
    }
  } else {
    // 2-limb Garner CRT reconstruction over Q = q0 * q1 and scale to p
    const auto* mod_params_q0 = main_moduli[0]->ModParams();
    const auto* mod_params_q1 = main_moduli[1]->ModParams();
    uint64_t q0 = main_moduli[0]->Modulus();
    uint64_t q1 = main_moduli[1]->Modulus();
    absl::uint128 total_q = static_cast<absl::uint128>(q0) * q1;
    uint64_t q0_inv_q1 = ModularInverse(q0, q1);

    const auto& c0_l0 = c0_poly.Coeffs()[0];
    const auto& c0_l1 = c0_poly.Coeffs()[1];
    const auto& c1_l0 = c1_poly.Coeffs()[0];
    const auto& c1_l1 = c1_poly.Coeffs()[1];

    for (size_t i = 0; i < n; ++i) {
      uint64_t val0_l0 = c0_l0[i].ExportInt(mod_params_q0);
      uint64_t val0_l1 = c0_l1[i].ExportInt(mod_params_q1);

      uint64_t val1_l0 = c1_l0[i].ExportInt(mod_params_q0);
      uint64_t val1_l1 = c1_l1[i].ExportInt(mod_params_q1);

      // Garner CRT for c0
      uint64_t diff0 = (val0_l1 >= (val0_l0 % q1))
                           ? (val0_l1 - (val0_l0 % q1))
                           : (q1 - ((val0_l0 % q1) - val0_l1));
      uint64_t x1_prime0 = absl::Uint128Low64(
          (static_cast<absl::uint128>(diff0) * q0_inv_q1) % q1);
      absl::uint128 v0_Q = static_cast<absl::uint128>(x1_prime0) * q0 + val0_l0;

      // Garner CRT for c1
      uint64_t diff1 = (val1_l1 >= (val1_l0 % q1))
                           ? (val1_l1 - (val1_l0 % q1))
                           : (q1 - ((val1_l0 % q1) - val1_l1));
      uint64_t x1_prime1 = absl::Uint128Low64(
          (static_cast<absl::uint128>(diff1) * q0_inv_q1) % q1);
      absl::uint128 v1_Q = static_cast<absl::uint128>(x1_prime1) * q0 + val1_l0;

      // Scale to p
      uint64_t s0 = static_cast<uint64_t>((v0_Q * p + (total_q / 2)) / total_q);
      uint64_t s1 = static_cast<uint64_t>((v1_Q * p + (total_q / 2)) / total_q);

      c0_scaled[i] = s0 % p;
      c1_scaled[i] = s1 % p;
    }
  }

  size_t num_bits = NumBitsForModulus(p);
  scaled_ct.set_c0(BitPack(c0_scaled, num_bits));
  scaled_ct.set_c1(BitPack(c1_scaled, num_bits));

  return scaled_ct;
}

absl::StatusOr<SerializedPirResponse> OnionPirServer::ProcessQuery(
    const SerializedPirQuery& query) const {
  RLWE_ASSIGN_OR_RETURN(auto packed_query, DeserializeQueryCiphertext(query));
  RLWE_ASSIGN_OR_RETURN(auto response_cts,
                        ProcessQueryDirect(packed_query, /*db=*/{}));

  SerializedPirResponse response;
  for (auto& ct : response_cts) {
    RLWE_ASSIGN_OR_RETURN(auto scaled_ct,
                          ScaleAndSwitchCiphertext(std::move(ct)));
    *response.add_scaled_ct() = std::move(scaled_ct);
  }
  return response;
}

absl::StatusOr<SerializedPirResponse> OnionPirServer::ProcessQuery(
    const SerializedPirQuery& query,
    absl::Span<const std::vector<uint8_t>> db) const {
  RLWE_ASSIGN_OR_RETURN(auto packed_query, DeserializeQueryCiphertext(query));
  RLWE_ASSIGN_OR_RETURN(auto response_cts,
                        ProcessQueryDirect(packed_query, db));

  SerializedPirResponse response;
  for (auto& ct : response_cts) {
    RLWE_ASSIGN_OR_RETURN(auto scaled_ct,
                          ScaleAndSwitchCiphertext(std::move(ct)));
    *response.add_scaled_ct() = std::move(scaled_ct);
  }
  return response;
}

}  // namespace oram_okvs::pir
