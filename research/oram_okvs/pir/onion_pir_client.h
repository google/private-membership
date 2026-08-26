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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_ONION_PIR_CLIENT_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_ONION_PIR_CLIENT_H_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/types/span.h"
#include "research/oram_okvs/pir/onion_pir_params.h"
#include "research/oram_okvs/pir/onion_pir_types.h"
#include "research/oram_okvs/pir/pir.h"
#include "research/oram_okvs/pir/pir.pb.h"
#include "shell_encryption/rns/coefficient_encoder.h"
#include "shell_encryption/rns/rns_bfv_ciphertext.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_gadget.h"
#include "shell_encryption/rns/rns_galois_key.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/rns_secret_key.h"

namespace oram_okvs::pir {

// OnionPIR Client implementation using RNS RLWE primitives.
class OnionPirClient : public PirClient {
 public:
  // Factory method to create and initialize an OnionPIR client.
  static absl::StatusOr<std::unique_ptr<OnionPirClient>> Create(
      OnionPirParams params);

  ~OnionPirClient() override = default;

  // PirClient interface methods.
  absl::StatusOr<SerializedEvaluationKeys> GetSerializedEvaluationKeys()
      const override;

  absl::StatusOr<SerializedPirQuery> CreateQuery(size_t index) const override;

  absl::StatusOr<std::vector<uint8_t>> DecodeResponse(
      const SerializedPirResponse& response,
      size_t entry_size_bytes) const override;

  size_t num_entries() const override { return params_.num_entries; }
  size_t entry_size_bytes() const override { return params_.entry_size_bytes; }

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

  const rlwe::RnsRlweSecretKey<ModularInt>& secret_key() const {
    return secret_key_;
  }

  const RgswCiphertext& rgsw_secret_key() const { return rgsw_secret_key_; }

  // Directly creates in-memory packed query ciphertext without serialization.
  absl::StatusOr<RnsBfvCiphertext> CreatePackedQueryCiphertext(
      size_t index) const;

  // Creates the (c1, prng_seed) components of the packed query ciphertext.
  absl::StatusOr<std::pair<RnsPolynomial, std::string>>
  CreatePackedQuerySeedComponents(size_t index) const;

  // Directly decrypts in-memory BFV response ciphertext to raw bytes.
  absl::StatusOr<std::vector<uint8_t>> DecryptResponseCiphertexts(
      absl::Span<const RnsBfvCiphertext> response_cts,
      size_t entry_size_bytes) const;

 private:
  explicit OnionPirClient(
      OnionPirParams params, std::unique_ptr<const RnsContext> rns_context,
      std::unique_ptr<const RnsErrorParams> error_params,
      std::unique_ptr<const rlwe::RnsGadget<ModularInt>> conversion_gadget,
      std::unique_ptr<const rlwe::RnsGadget<ModularInt>> expansion_gadget,
      std::unique_ptr<const rlwe::RnsGadget<ModularInt>> trace_gadget,
      std::unique_ptr<const rlwe::CoefficientEncoder<ModularInt>> encoder,
      rlwe::RnsRlweSecretKey<ModularInt> secret_key,
      std::vector<rlwe::RnsGaloisKey<ModularInt>> galois_keys,
      RgswCiphertext rgsw_secret_key, std::vector<int64_t> sk_coeffs);

  OnionPirParams params_;
  std::unique_ptr<const RnsContext> rns_context_;
  std::unique_ptr<const RnsErrorParams> error_params_;
  std::unique_ptr<const rlwe::RnsGadget<ModularInt>> conversion_gadget_;
  std::unique_ptr<const rlwe::RnsGadget<ModularInt>> expansion_gadget_;
  std::unique_ptr<const rlwe::RnsGadget<ModularInt>> trace_gadget_;
  std::unique_ptr<const rlwe::CoefficientEncoder<ModularInt>> encoder_;
  rlwe::RnsRlweSecretKey<ModularInt> secret_key_;
  std::vector<rlwe::RnsGaloisKey<ModularInt>> galois_keys_;
  RgswCiphertext rgsw_secret_key_;
  std::vector<int64_t> sk_coeffs_;
};

}  // namespace oram_okvs::pir

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_ONION_PIR_CLIENT_H_
