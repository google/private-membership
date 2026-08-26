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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_ONION_PIR_TYPES_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_ONION_PIR_TYPES_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "research/oram_okvs/pir/fast_ntt.h"
#include "research/oram_okvs/pir/onion_pir_params.h"
#include "research/oram_okvs/pir/pir.pb.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/ntt_parameters.h"
#include "shell_encryption/rns/rns_bfv_ciphertext.h"
#include "shell_encryption/rns/rns_ciphertext.h"
#include "shell_encryption/rns/rns_context.h"
#include "shell_encryption/rns/rns_error_params.h"
#include "shell_encryption/rns/rns_gadget.h"
#include "shell_encryption/rns/rns_modulus.h"
#include "shell_encryption/rns/rns_polynomial.h"
#include "shell_encryption/rns/rns_serialization.pb.h"
#include "shell_encryption/status_macros.h"

namespace oram_okvs::pir {

using ModularInt = rlwe::MontgomeryInt<rlwe::Uint64>;
using ModularIntParams = ModularInt::Params;
using PrimeModulus = rlwe::PrimeModulus<ModularInt>;
using RnsPolynomial = rlwe::RnsPolynomial<ModularInt>;
using RnsBfvCiphertext = rlwe::RnsBfvCiphertext<ModularInt>;
using RnsContext = rlwe::RnsContext<ModularInt>;
using RnsErrorParams = rlwe::RnsErrorParams<ModularInt>;

// Type aliases for seamless compatibility.
using Polynomial = RnsPolynomial;
using BfvCiphertext = RnsBfvCiphertext;

// Computes the modular inverse of `a` modulo `m` using the Extended Euclidean
// Algorithm. Requires gcd(a, m) == 1.
inline rlwe::Uint64 ModularInverse(rlwe::Uint64 a, rlwe::Uint64 m) {
  int64_t t = 0, newt = 1;
  int64_t r = static_cast<int64_t>(m);
  int64_t newr = static_cast<int64_t>(a % m);
  while (newr != 0) {
    int64_t quotient = r / newr;
    int64_t temp_t = t - quotient * newt;
    t = newt;
    newt = temp_t;
    int64_t temp_r = r - quotient * newr;
    r = newr;
    newr = temp_r;
  }
  if (t < 0) {
    t += m;
  }
  return static_cast<rlwe::Uint64>(t);
}

// Computes the bit reversal of `val` using `num_bits` bits.
inline uint32_t BitReverse(uint32_t val, size_t num_bits) {
  uint32_t res = 0;
  for (size_t i = 0; i < num_bits; ++i) {
    if ((val >> i) & 1) {
      res |= (1u << (num_bits - 1 - i));
    }
  }
  return res;
}

// Returns the minimum number of bits needed to represent any value in [0, p -
// 1].
inline size_t NumBitsForModulus(uint64_t p) {
  if (p <= 1) return 1;
  size_t bits = 0;
  while ((static_cast<uint64_t>(1) << bits) < p) {
    ++bits;
  }
  return bits;
}

// Bit-packs a slice of uint64 coefficients (each requiring `num_bits`) into a
// compact byte string.
inline std::string BitPack(absl::Span<const uint64_t> coeffs, size_t num_bits) {
  if (num_bits == 0 || coeffs.empty()) return "";
  size_t total_bits = coeffs.size() * num_bits;
  size_t total_bytes = (total_bits + 7) / 8;
  std::string buffer(total_bytes, '\0');
  uint8_t* out = reinterpret_cast<uint8_t*>(buffer.data());

  uint64_t bit_buf = 0;
  size_t bits_in_buf = 0;
  size_t byte_idx = 0;

  for (uint64_t val : coeffs) {
    bit_buf |= (val << bits_in_buf);
    bits_in_buf += num_bits;
    while (bits_in_buf >= 8) {
      out[byte_idx++] = static_cast<uint8_t>(bit_buf & 0xFF);
      bit_buf >>= 8;
      bits_in_buf -= 8;
    }
  }
  if (bits_in_buf > 0) {
    out[byte_idx++] = static_cast<uint8_t>(bit_buf & 0xFF);
  }
  return buffer;
}

// Unpacks `num_coeffs` integers of `num_bits` each from a bit-packed byte
// string.
inline absl::StatusOr<std::vector<uint64_t>> BitUnpack(absl::string_view bytes,
                                                       size_t num_coeffs,
                                                       size_t num_bits) {
  std::vector<uint64_t> coeffs(num_coeffs, 0);
  if (num_bits == 0 || num_coeffs == 0) return coeffs;

  size_t expected_bytes = (num_coeffs * num_bits + 7) / 8;
  if (bytes.size() < expected_bytes) {
    return absl::InvalidArgumentError(
        "Malformed bit-packed buffer: insufficient bytes.");
  }

  const uint8_t* in = reinterpret_cast<const uint8_t*>(bytes.data());
  size_t in_size = bytes.size();
  size_t byte_idx = 0;
  uint64_t bit_buf = 0;
  size_t bits_in_buf = 0;
  uint64_t mask = (num_bits >= 64) ? ~0ULL : ((1ULL << num_bits) - 1ULL);

  for (size_t i = 0; i < num_coeffs; ++i) {
    while (bits_in_buf < num_bits && byte_idx < in_size) {
      bit_buf |= (static_cast<uint64_t>(in[byte_idx++]) << bits_in_buf);
      bits_in_buf += 8;
    }
    coeffs[i] = bit_buf & mask;
    bit_buf >>= num_bits;
    bits_in_buf = (bits_in_buf >= num_bits) ? (bits_in_buf - num_bits) : 0;
  }
  return coeffs;
}

// Represents an RGSW ciphertext consisting of 2*ell BFV ciphertexts in NTT
// form.
struct RgswCiphertext {
  std::vector<RnsBfvCiphertext> rows;  // 2*ell rows

  RgswCiphertext() = default;
  explicit RgswCiphertext(std::vector<RnsBfvCiphertext> r)
      : rows(std::move(r)) {}

  absl::StatusOr<SerializedRgswCiphertext> Serialize() const {
    SerializedRgswCiphertext proto;
    for (const auto& row : rows) {
      RLWE_ASSIGN_OR_RETURN(*proto.add_rows(), row.Serialize());
    }
    return proto;
  }

  static absl::StatusOr<RgswCiphertext> Deserialize(
      const SerializedRgswCiphertext& proto,
      std::vector<const PrimeModulus*> moduli,
      const RnsErrorParams* error_params, const RnsContext* context) {
    std::vector<RnsBfvCiphertext> r;
    r.reserve(proto.rows_size());
    for (int i = 0; i < proto.rows_size(); ++i) {
      RLWE_ASSIGN_OR_RETURN(auto rlwe_ct,
                            rlwe::RnsRlweCiphertext<ModularInt>::Deserialize(
                                proto.rows(i), moduli, error_params));
      RLWE_ASSIGN_OR_RETURN(auto c0, rlwe_ct.Component(0));
      RLWE_ASSIGN_OR_RETURN(auto c1, rlwe_ct.Component(1));
      std::vector<RnsPolynomial> comps;
      comps.reserve(2);
      comps.push_back(std::move(c0));
      comps.push_back(std::move(c1));
      r.emplace_back(std::move(comps), moduli, rlwe_ct.PowerOfS(),
                     rlwe_ct.Error(), error_params, context);
    }
    return RgswCiphertext(std::move(r));
  }
};

// Decomposes a polynomial in NTT form into `gadget->Dimension()` polynomials in
// NTT form using rlwe::RnsGadget.
inline absl::StatusOr<std::vector<RnsPolynomial>> GadgetDecompose(
    const rlwe::RnsGadget<ModularInt>* gadget, const RnsPolynomial& poly_ntt,
    absl::Span<const PrimeModulus* const> moduli) {
  if (gadget == nullptr) {
    return absl::InvalidArgumentError("`gadget` must not be null.");
  }
  RnsPolynomial poly_coeff = poly_ntt;
  if (poly_coeff.IsNttForm()) {
    RLWE_RETURN_IF_ERROR(FastConvertToCoeffForm(poly_coeff, moduli));
  }
  RLWE_ASSIGN_OR_RETURN(std::vector<RnsPolynomial> digits,
                        gadget->Decompose(poly_coeff, moduli));
  for (auto& digit : digits) {
    RLWE_RETURN_IF_ERROR(FastConvertToNttForm(digit, moduli));
  }
  return digits;
}

// Evaluates the External Product between an RGSW ciphertext and a BFV
// ciphertext: ExternalProduct(RGSW, BFV) -> BFV ciphertext encrypting the
// homomorphic product.
inline absl::StatusOr<RnsBfvCiphertext> ExternalProduct(
    const RgswCiphertext& rgsw, const RnsBfvCiphertext& bfv,
    absl::Span<const PrimeModulus* const> moduli,
    const rlwe::RnsGadget<ModularInt>* gadget, const RnsContext* context,
    const RnsErrorParams* error_params) {
  if (gadget == nullptr) {
    return absl::InvalidArgumentError("`gadget` must not be null.");
  }
  size_t gadget_dim = gadget->Dimension();
  if (rgsw.rows.size() != 2 * gadget_dim) {
    return absl::InvalidArgumentError(
        absl::StrCat("RGSW row count does not match 2 * gadget_dim (",
                     rgsw.rows.size(), " vs ", 2 * gadget_dim, ")"));
  }

  RLWE_ASSIGN_OR_RETURN(auto bfv_c0, bfv.Component(0));
  RLWE_ASSIGN_OR_RETURN(auto bfv_c1, bfv.Component(1));

  RLWE_ASSIGN_OR_RETURN(std::vector<RnsPolynomial> d0,
                        GadgetDecompose(gadget, bfv_c0, moduli));
  RLWE_ASSIGN_OR_RETURN(std::vector<RnsPolynomial> d1,
                        GadgetDecompose(gadget, bfv_c1, moduli));

  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial res_c0,
      RnsPolynomial::CreateZero(bfv.LogN(), moduli, /*is_ntt=*/true));
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial res_c1,
      RnsPolynomial::CreateZero(bfv.LogN(), moduli, /*is_ntt=*/true));

  for (size_t k = 0; k < gadget_dim; ++k) {
    RLWE_ASSIGN_OR_RETURN(auto row_k_c0, rgsw.rows[k].Component(0));
    RLWE_ASSIGN_OR_RETURN(auto row_k_c1, rgsw.rows[k].Component(1));
    RLWE_RETURN_IF_ERROR(res_c0.FusedMulAddInPlace(d0[k], row_k_c0, moduli));
    RLWE_RETURN_IF_ERROR(res_c1.FusedMulAddInPlace(d0[k], row_k_c1, moduli));

    RLWE_ASSIGN_OR_RETURN(auto row_ell_k_c0,
                          rgsw.rows[gadget_dim + k].Component(0));
    RLWE_ASSIGN_OR_RETURN(auto row_ell_k_c1,
                          rgsw.rows[gadget_dim + k].Component(1));
    RLWE_RETURN_IF_ERROR(
        res_c0.FusedMulAddInPlace(d1[k], row_ell_k_c0, moduli));
    RLWE_RETURN_IF_ERROR(
        res_c1.FusedMulAddInPlace(d1[k], row_ell_k_c1, moduli));
  }

  std::vector<RnsPolynomial> res_comps;
  res_comps.reserve(2);
  res_comps.push_back(std::move(res_c0));
  res_comps.push_back(std::move(res_c1));

  std::vector<const PrimeModulus*> moduli_vec(moduli.begin(), moduli.end());
  return RnsBfvCiphertext(std::move(res_comps), std::move(moduli_vec),
                          bfv.PowerOfS(), bfv.Error(), error_params, context);
}

// Encodes raw byte span into an NTT polynomial for database evaluation.
// Byte values are placed directly into coefficients in [0, 255] < t.
inline absl::StatusOr<RnsPolynomial> EncodeBytesToPolynomialNtt(
    absl::Span<const uint8_t> bytes, size_t poly_degree,
    absl::Span<const PrimeModulus* const> moduli) {
  std::vector<std::vector<ModularInt>> rns_coeffs;
  rns_coeffs.reserve(moduli.size());
  size_t count = std::min(bytes.size(), poly_degree);
  for (const auto* prime_mod : moduli) {
    const auto* mod_params = prime_mod->ModParams();
    std::vector<ModularInt> coeffs(poly_degree,
                                   ModularInt::ImportZero(mod_params));
    for (size_t i = 0; i < count; ++i) {
      RLWE_ASSIGN_OR_RETURN(
          coeffs[i], ModularInt::ImportInt(static_cast<rlwe::Uint64>(bytes[i]),
                                           mod_params));
    }
    rns_coeffs.push_back(std::move(coeffs));
  }
  RLWE_ASSIGN_OR_RETURN(
      RnsPolynomial poly,
      RnsPolynomial::Create(std::move(rns_coeffs), /*is_ntt=*/false));
  RLWE_RETURN_IF_ERROR(FastConvertToNttForm(poly, moduli));
  return poly;
}

// Decodes coefficients from decrypted polynomial into raw bytes.
inline std::vector<uint8_t> DecodePolynomialCoeffsToBytes(
    absl::Span<const rlwe::Uint64> coeffs, size_t expected_bytes) {
  std::vector<uint8_t> bytes(expected_bytes, 0);
  size_t count = std::min(coeffs.size(), expected_bytes);
  for (size_t i = 0; i < count; ++i) {
    bytes[i] = static_cast<uint8_t>(coeffs[i] & 0xFF);
  }
  return bytes;
}

}  // namespace oram_okvs::pir

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_ONION_PIR_TYPES_H_
