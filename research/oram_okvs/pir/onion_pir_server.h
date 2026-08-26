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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_ONION_PIR_SERVER_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_ONION_PIR_SERVER_H_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "research/oram_okvs/pir/onion_pir_params.h"
#include "research/oram_okvs/pir/onion_pir_types.h"
#include "research/oram_okvs/pir/pir.h"
#include "research/oram_okvs/pir/pir.pb.h"
#include "shell_encryption/rns/rns_bfv_ciphertext.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_gadget.h"
#include "shell_encryption/rns/rns_galois_key.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"

namespace oram_okvs::pir {

// OnionPIR Server implementation using RNS RLWE primitives.
class OnionPirServer : public PirServer {
 public:
  // Factory method to create and initialize an OnionPIR server.
  static absl::StatusOr<std::unique_ptr<OnionPirServer>> Create(
      OnionPirParams params);

  ~OnionPirServer() override = default;

  // PirServer interface methods.
  absl::Status SetEvaluationKeys(const SerializedEvaluationKeys& keys) override;

  absl::Status SetDatabase(absl::Span<const std::vector<uint8_t>> db) override;

  absl::StatusOr<SerializedPirResponse> ProcessQuery(
      const SerializedPirQuery& query) const override;

  absl::StatusOr<SerializedPirResponse> ProcessQuery(
      const SerializedPirQuery& query,
      absl::Span<const std::vector<uint8_t>> db) const override;

  size_t num_entries() const override { return params_.num_entries; }
  size_t entry_size_bytes() const override { return params_.entry_size_bytes; }

  // Sets in-memory evaluation keys directly.
  absl::Status SetEvaluationKeysDirect(
      std::vector<rlwe::RnsGaloisKey<ModularInt>> galois_keys,
      RgswCiphertext rgsw_secret_key);

  // Directly unpacks a packed query into CBFV and CRGSW ciphertexts.
  absl::StatusOr<
      std::pair<std::vector<RnsBfvCiphertext>, std::vector<RgswCiphertext>>>
  UnpackQuery(const RnsBfvCiphertext& packed_query,
              size_t active_db_size = 0) const;

  // Expands the packed query using the BFV expansion tree.
  absl::StatusOr<std::vector<RnsBfvCiphertext>> ExpandBfv(
      const RnsBfvCiphertext& packed_query, size_t target_u = 0) const;

  // Evaluates Homomorphic Trace to clear non-stride garbage coefficients.
  absl::StatusOr<RnsBfvCiphertext> HomomorphicTrace(RnsBfvCiphertext ct) const;

  // Directly evaluates query ciphertexts on an in-memory or on-the-fly
  // database.
  absl::StatusOr<std::vector<RnsBfvCiphertext>> ProcessQueryDirect(
      const RnsBfvCiphertext& packed_query,
      absl::Span<const std::vector<uint8_t>> db) const;

  // Evaluates unpacked query vectors (CBFV and CRGSW) directly on database.
  absl::StatusOr<std::vector<RnsBfvCiphertext>> EvaluateDatabaseDirect(
      absl::Span<const RnsBfvCiphertext> c_bfv,
      absl::Span<const RgswCiphertext> c_rgsw,
      absl::Span<const std::vector<uint8_t>> db) const;

  // Internal accessors.
  const OnionPirParams& params() const { return params_; }
  const RnsContext* rns_context() const { return rns_context_.get(); }
  const RnsErrorParams* error_params() const { return error_params_.get(); }
  const rlwe::RnsGadget<ModularInt>* conversion_gadget() const {
    return conversion_gadget_.get();
  }
  const rlwe::RnsGadget<ModularInt>* expansion_gadget() const {
    return expansion_gadget_.get();
  }
  const rlwe::RnsGadget<ModularInt>* trace_gadget() const {
    return trace_gadget_.get();
  }
  const rlwe::RnsGadget<ModularInt>* gadget() const {
    return conversion_gadget_.get();
  }
  std::vector<const PrimeModulus*> main_moduli() const {
    return rns_context_->MainPrimeModuli();
  }
  const ModularIntParams* mod_params() const {
    return rns_context_->MainPrimeModuli()[0]->ModParams();
  }

  // Deserializes a SerializedPirQuery with compressed coefficients.
  absl::StatusOr<RnsBfvCiphertext> DeserializeQueryCiphertext(
      const SerializedPirQuery& query) const;

  // Scales and rounds a Level 0 BFV ciphertext to an arbitrary response modulus
  // p.
  absl::StatusOr<SerializedScaledCiphertext> ScaleAndSwitchCiphertext(
      RnsBfvCiphertext ct) const;

 private:
  explicit OnionPirServer(
      OnionPirParams params, std::unique_ptr<const RnsContext> rns_context,
      std::unique_ptr<const RnsErrorParams> error_params,
      std::unique_ptr<const rlwe::RnsGadget<ModularInt>> conversion_gadget,
      std::unique_ptr<const rlwe::RnsGadget<ModularInt>> expansion_gadget,
      std::unique_ptr<const rlwe::RnsGadget<ModularInt>> trace_gadget,
      std::vector<RnsPolynomial> neg_monomials_ntt);

  OnionPirParams params_;
  std::unique_ptr<const RnsContext> rns_context_;
  std::unique_ptr<const RnsErrorParams> error_params_;
  std::unique_ptr<const rlwe::RnsGadget<ModularInt>> conversion_gadget_;
  std::unique_ptr<const rlwe::RnsGadget<ModularInt>> expansion_gadget_;
  std::unique_ptr<const rlwe::RnsGadget<ModularInt>> trace_gadget_;
  std::vector<RnsPolynomial>
      neg_monomials_ntt_;  // Precomputed -x^{n-k} in NTT form

  std::vector<rlwe::RnsGaloisKey<ModularInt>> galois_keys_;
  RgswCiphertext rgsw_secret_key_;
  bool keys_loaded_ = false;

  // Optional precomputed NTT representation of loaded database:
  // precomputed_db_ntt_[entry_index][poly_chunk_index]
  std::vector<std::vector<RnsPolynomial>> precomputed_db_ntt_;
  bool db_loaded_ = false;
};

}  // namespace oram_okvs::pir

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_ONION_PIR_SERVER_H_
