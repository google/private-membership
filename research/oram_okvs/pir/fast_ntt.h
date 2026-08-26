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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_FAST_NTT_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_FAST_NTT_H_

#include <cstddef>
#include <cstdint>
#include <vector>

#include "absl/status/status.h"
#include "absl/types/span.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"

namespace oram_okvs::pir {

// ============================================================================
// High-Performance Vectorized Number Theoretic Transform (NTT)
// ============================================================================
//
// Optimizations over baseline third_party/rlwe NTT:
// 1. Fully vectorized Cooley-Tukey (Forward) and Gentleman-Sande (Inverse)
//    butterflies using Highway SIMD (AVX2 / AVX-512 / NEON).
// 2. Vectorized intra-register butterflies for bottom layers (m < num_lanes),
//    eliminating all scalar fallbacks.
// 3. Harvey's lazy reduction: keeps intermediate values in [0, 4q) and defers
//    modular reduction to every 2 layers, reducing branch/mask overhead by
//    >75%.
// 4. 4x unrolled vector loops for outer layers to saturate CPU execution
// pipelines.
// 5. Zero-copy in-place transform on RnsPolynomial and raw coefficient spans.

// Forward NTT on a vector of ModularInt coefficients in-place.
template <typename ModularInt>
absl::Status FastForwardNtt(std::vector<ModularInt>& coeffs,
                            const rlwe::NttParameters<ModularInt>& ntt_params,
                            const typename ModularInt::Params& mod_params);

// Inverse NTT on a vector of ModularInt coefficients in-place (includes n^-1
// scaling).
template <typename ModularInt>
absl::Status FastInverseNtt(std::vector<ModularInt>& coeffs,
                            const rlwe::NttParameters<ModularInt>& ntt_params,
                            const typename ModularInt::Params& mod_params);

// Fast forward NTT conversion on RnsPolynomial in-place across all RNS moduli.
template <typename ModularInt>
absl::Status FastConvertToNttForm(
    rlwe::RnsPolynomial<ModularInt>& poly,
    absl::Span<const rlwe::PrimeModulus<ModularInt>* const> moduli);

template <typename ModularInt>
inline absl::Status FastConvertToNttForm(
    rlwe::RnsPolynomial<ModularInt>& poly,
    absl::Span<const rlwe::PrimeModulus<ModularInt>*> moduli) {
  return FastConvertToNttForm(
      poly, absl::Span<const rlwe::PrimeModulus<ModularInt>* const>(
                moduli.data(), moduli.size()));
}

template <typename ModularInt>
inline absl::Status FastConvertToNttForm(
    rlwe::RnsPolynomial<ModularInt>& poly,
    absl::Span<rlwe::PrimeModulus<ModularInt>*> moduli) {
  return FastConvertToNttForm(
      poly, absl::Span<const rlwe::PrimeModulus<ModularInt>* const>(
                reinterpret_cast<const rlwe::PrimeModulus<ModularInt>* const*>(
                    moduli.data()),
                moduli.size()));
}

template <typename ModularInt>
inline absl::Status FastConvertToNttForm(
    rlwe::RnsPolynomial<ModularInt>& poly,
    const std::vector<const rlwe::PrimeModulus<ModularInt>*>& moduli) {
  return FastConvertToNttForm(
      poly, absl::Span<const rlwe::PrimeModulus<ModularInt>* const>(
                moduli.data(), moduli.size()));
}

template <typename ModularInt>
inline absl::Status FastConvertToNttForm(
    rlwe::RnsPolynomial<ModularInt>& poly,
    const std::vector<rlwe::PrimeModulus<ModularInt>*>& moduli) {
  return FastConvertToNttForm(
      poly, absl::Span<const rlwe::PrimeModulus<ModularInt>* const>(
                reinterpret_cast<const rlwe::PrimeModulus<ModularInt>* const*>(
                    moduli.data()),
                moduli.size()));
}

// Fast inverse NTT conversion on RnsPolynomial in-place across all RNS moduli.
template <typename ModularInt>
absl::Status FastConvertToCoeffForm(
    rlwe::RnsPolynomial<ModularInt>& poly,
    absl::Span<const rlwe::PrimeModulus<ModularInt>* const> moduli);

template <typename ModularInt>
inline absl::Status FastConvertToCoeffForm(
    rlwe::RnsPolynomial<ModularInt>& poly,
    absl::Span<const rlwe::PrimeModulus<ModularInt>*> moduli) {
  return FastConvertToCoeffForm(
      poly, absl::Span<const rlwe::PrimeModulus<ModularInt>* const>(
                moduli.data(), moduli.size()));
}

template <typename ModularInt>
inline absl::Status FastConvertToCoeffForm(
    rlwe::RnsPolynomial<ModularInt>& poly,
    absl::Span<rlwe::PrimeModulus<ModularInt>*> moduli) {
  return FastConvertToCoeffForm(
      poly, absl::Span<const rlwe::PrimeModulus<ModularInt>* const>(
                reinterpret_cast<const rlwe::PrimeModulus<ModularInt>* const*>(
                    moduli.data()),
                moduli.size()));
}

template <typename ModularInt>
inline absl::Status FastConvertToCoeffForm(
    rlwe::RnsPolynomial<ModularInt>& poly,
    const std::vector<const rlwe::PrimeModulus<ModularInt>*>& moduli) {
  return FastConvertToCoeffForm(
      poly, absl::Span<const rlwe::PrimeModulus<ModularInt>* const>(
                moduli.data(), moduli.size()));
}

template <typename ModularInt>
inline absl::Status FastConvertToCoeffForm(
    rlwe::RnsPolynomial<ModularInt>& poly,
    const std::vector<rlwe::PrimeModulus<ModularInt>*>& moduli) {
  return FastConvertToCoeffForm(
      poly, absl::Span<const rlwe::PrimeModulus<ModularInt>* const>(
                reinterpret_cast<const rlwe::PrimeModulus<ModularInt>* const*>(
                    moduli.data()),
                moduli.size()));
}

}  // namespace oram_okvs::pir

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_FAST_NTT_H_
