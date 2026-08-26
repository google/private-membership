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

#ifndef PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_ONION_PIR_PARAMS_H_
#define PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_ONION_PIR_PARAMS_H_

#include <cmath>
#include <cstddef>
#include <cstdint>

#include "research/oram_okvs/pir/pir.pb.h"
#include "shell_encryption/constants.h"
#include "shell_encryption/integral_types.h"

namespace oram_okvs::pir {

// Forward declaration
struct OnionPirParams;

// Holds noise variance and standard deviation components computed under the
// independence heuristic for OnionPIR.
struct NoiseEstimate {
  double fresh_noise_variance = 0.0;
  double keyswitch_noise_variance = 0.0;
  double unpacked_query_noise_variance = 0.0;
  double dim0_dot_product_noise_variance = 0.0;
  double higher_dims_noise_variance = 0.0;
  double server_eval_noise_variance = 0.0;
  double final_noise_variance = 0.0;
  double final_noise_std_dev = 0.0;
  double max_noise_capacity = 0.0;  // Q_dec / (2 * t)
  double noise_budget_bits = 0.0;   // log2(max_noise_capacity / (C * std_dev))
  double confidence_factor = 6.0;   // C (e.g. 6.0 for 6-sigma tail bound)
  bool is_decryption_correct = false;
};

NoiseEstimate EstimateOnionPirNoise(const OnionPirParams& params,
                                    double confidence_factor = 6.0);

// Configuration parameters for the OnionPIR protocol.
struct OnionPirParams {
  size_t num_entries = 1024;     // Total number of logical entries N in DB.
  size_t entry_size_bytes = 64;  // Entry payload size in bytes.
  size_t poly_degree = 2048;     // RLWE polynomial degree n (power of 2).
  size_t n0 = 32;                // Initial dimension size N0 (e.g. 16, 32, 64).
  std::vector<rlwe::Uint64> ciphertext_moduli = {
      134176769ULL,
      134111233ULL};  // Two 27-bit NTT primes (chain product ~ 2^54).
  size_t log_t =
      8;  // Plaintext modulus bit-length: t = (1 << log_t) + 1 = 257.
  size_t variance = 10;  // Error distribution variance (sigma^2 = 10).
  // 1. RLWE -> RGSW conversion key gadget parameter (for higher dimensions d >
  // 1).
  size_t conversion_gadget_log_base = 7;         // log2(B_conv)
  std::vector<size_t> conversion_gadget_log_bs;  // Optional per-modulus log_bs.

  // 2. Expansion tree Galois key gadget parameter (for slot expansion
  // ExpandBfv).
  size_t expansion_gadget_log_base = 6;         // log2(B_exp)
  std::vector<size_t> expansion_gadget_log_bs;  // Optional per-modulus log_bs.

  // 3. Homomorphic trace Galois key gadget parameter (for query trace
  // reduction).
  size_t trace_gadget_log_base = 5;         // log2(B_trace)
  std::vector<size_t> trace_gadget_log_bs;  // Optional per-modulus log_bs.

  rlwe::Uint64 response_modulus =
      262144ULL;  // Response modulus p (e.g. 262144 = 2^18).

  size_t num_query_coeffs =
      0;  // Number of query coefficients sent (w_sent). 0 means minimal w.

  // Helper to compute gadget dimension ell across ciphertext_moduli for given
  // log_bs.
  size_t GadgetEllForLogBs(const std::vector<size_t>& log_bs,
                           size_t default_log_b) const {
    if (ciphertext_moduli.empty() || default_log_b == 0) return 0;
    size_t total_ell = 0;
    for (size_t i = 0; i < ciphertext_moduli.size(); ++i) {
      size_t log_b =
          (i < log_bs.size() && log_bs[i] > 0) ? log_bs[i] : default_log_b;
      if (log_b == 0) continue;
      size_t log2_q = 0;
      while ((static_cast<rlwe::Uint64>(1) << log2_q) < ciphertext_moduli[i]) {
        ++log2_q;
      }
      total_ell += (log2_q + log_b - 1) / log_b;
    }
    return total_ell;
  }

  // Computes the gadget dimension ell for a uniform log_base across
  // ciphertext_moduli.
  size_t GadgetEllForLogBase(size_t log_b) const {
    return GadgetEllForLogBs(/*log_bs=*/{}, log_b);
  }

  // Returns effective conversion gadget log base.
  size_t conversion_log_base() const {
    return conversion_gadget_log_base > 0 ? conversion_gadget_log_base : 7;
  }

  // Returns effective expansion gadget log base.
  size_t expansion_log_base() const {
    if (expansion_gadget_log_base > 0) return expansion_gadget_log_base;
    if (trace_gadget_log_base > 0) return trace_gadget_log_base;
    return 5;
  }

  // Returns effective trace gadget log base.
  size_t trace_log_base() const {
    return trace_gadget_log_base > 0 ? trace_gadget_log_base : 5;
  }

  // Returns effective conversion gadget dimension ell derived from moduli and
  // log_base.
  size_t conversion_ell() const {
    return GadgetEllForLogBs(conversion_gadget_log_bs, conversion_log_base());
  }

  // Returns effective expansion gadget dimension ell derived from moduli and
  // log_base.
  size_t expansion_ell() const {
    return GadgetEllForLogBs(expansion_gadget_log_bs, expansion_log_base());
  }

  // Returns effective trace gadget dimension ell derived from moduli and
  // log_base.
  size_t trace_ell() const {
    return GadgetEllForLogBs(trace_gadget_log_bs, trace_log_base());
  }

  // Backward compatibility aliases for conversion gadget dimension.
  size_t gadget_ell() const { return conversion_ell(); }
  size_t ell() const { return conversion_ell(); }

  // Computes the estimated noise under the independence heuristic.
  NoiseEstimate EstimateNoise(double confidence_factor = 6.0) const {
    return EstimateOnionPirNoise(*this, confidence_factor);
  }

  // Returns the noise budget in bits (positive means correct with high
  // probability).
  double NoiseBudgetBits(double confidence_factor = 6.0) const {
    return EstimateNoise(confidence_factor).noise_budget_bits;
  }

  // Primary (first) ciphertext prime modulus q_0.
  rlwe::Uint64 primary_ciphertext_modulus() const {
    return ciphertext_moduli.empty() ? 134176769ULL : ciphertext_moduli[0];
  }

  // Computes the number of dimensions d in the hypercube.
  // d = 1 when num_entries <= n0, otherwise 1 + ceil(log2(ceil(num_entries /
  // n0))).
  size_t num_dimensions() const {
    if (num_entries <= n0) return 1;
    size_t high_dim_size = (num_entries + n0 - 1) / n0;
    size_t log2_high = 0;
    while ((static_cast<size_t>(1) << log2_high) < high_dim_size) {
      ++log2_high;
    }
    return 1 + log2_high;
  }

  // Returns the total number of higher dimensions (d - 1).
  size_t num_higher_dimensions() const { return num_dimensions() - 1; }

  // Total number of packed values u = N0 + ell * (d - 1).
  size_t total_packed_values() const {
    return n0 + conversion_ell() * num_higher_dimensions();
  }

  // Smallest power of two w >= u.
  size_t w() const {
    size_t u = total_packed_values();
    size_t w_val = 1;
    while (w_val < u) {
      w_val <<= 1;
    }
    return w_val;
  }

  // Tree height = log2(w).
  size_t expansion_tree_height() const {
    size_t w_val = w();
    size_t h = 0;
    while ((static_cast<size_t>(1) << h) < w_val) {
      ++h;
    }
    return h;
  }

  // Number of query coefficients sent (w_sent), power of two in [w(),
  // poly_degree].
  size_t query_w() const {
    size_t min_w = w();
    if (num_query_coeffs <= min_w) {
      return min_w;
    }
    size_t w_val = 1;
    while (w_val < num_query_coeffs) {
      w_val <<= 1;
    }
    return std::min(w_val, poly_degree);
  }

  // Expansion tree height starting from query_w: log2(query_w).
  size_t query_expansion_tree_height() const {
    size_t w_val = query_w();
    size_t h = 0;
    while ((static_cast<size_t>(1) << h) < w_val) {
      ++h;
    }
    return h;
  }

  // Homomorphic trace steps needed on server: log2(N / query_w).
  size_t query_trace_steps() const {
    size_t log_n = 0;
    while ((static_cast<size_t>(1) << log_n) < poly_degree) {
      ++log_n;
    }
    size_t tree_h = query_expansion_tree_height();
    return log_n > tree_h ? (log_n - tree_h) : 0;
  }

  // Configures num_query_coeffs to send as many coefficients as fit within
  // max_bytes.
  void SetMaxQueryBytes(size_t max_bytes) {
    size_t bits_per_coeff = 0;
    for (auto q_val : ciphertext_moduli) {
      size_t b = 0;
      while ((static_cast<rlwe::Uint64>(1) << b) < q_val) {
        ++b;
      }
      bits_per_coeff += b;
    }
    if (bits_per_coeff == 0) bits_per_coeff = 54;
    size_t seed_bytes = 32;
    if (max_bytes <= seed_bytes) {
      num_query_coeffs = w();
      return;
    }
    size_t payload_bytes = max_bytes - seed_bytes;
    size_t max_coeffs = (payload_bytes * 8) / bits_per_coeff;
    size_t min_w = w();
    size_t cand = poly_degree;
    while (cand > min_w && cand > max_coeffs) {
      cand >>= 1;
    }
    num_query_coeffs = cand;
  }

  // Effective plaintext modulus t = (1 << log_t) + 1.
  rlwe::Uint64 plaintext_modulus() const {
    return (static_cast<rlwe::Uint64>(1) << log_t) + 1;
  }

  // Scaling factor Delta = floor(q_0 / t).
  rlwe::Uint64 delta() const {
    return primary_ciphertext_modulus() / plaintext_modulus();
  }

  // Usable payload bytes that fit within a single polynomial of degree n.
  // Using floor(log_t / 8) bytes per coefficient for direct byte mapping,
  // or bit-packing across all n coefficients.
  // With log_t = 8, 1 byte (8 bits) per coefficient is mapped directly
  // with coefficient in [0, 255] < t = 257. Thus, degree n provides n bytes
  // capacity.
  size_t bytes_per_polynomial() const { return poly_degree; }

  // Number of polynomials required to store an entry of size
  // `entry_size_bytes`.
  size_t num_polynomials_per_entry() const {
    size_t bytes_per_poly = bytes_per_polynomial();
    return (entry_size_bytes + bytes_per_poly - 1) / bytes_per_poly;
  }

  // Returns default parameters suitable for ORAM (e.g. 64-byte blocks, 1024
  // blocks).
  static OnionPirParams CreateForOram(size_t num_blocks, size_t block_size) {
    OnionPirParams params;
    params.num_entries = num_blocks;
    params.entry_size_bytes = block_size;
    params.poly_degree = 2048;
    params.n0 = 32;
    params.ciphertext_moduli = {134176769ULL, 134111233ULL};
    params.log_t = 8;
    params.variance = 10;
    params.conversion_gadget_log_base = 7;
    params.expansion_gadget_log_base = 6;
    params.trace_gadget_log_base = 5;
    params.response_modulus = 262144ULL;
    return params;
  }

  // Returns 128-bit secure RLWE parameters for OnionPIR (n = 2048).
  static OnionPirParams Create128BitSecure(size_t num_blocks, size_t block_size,
                                           size_t n0 = 32) {
    OnionPirParams params;
    params.num_entries = num_blocks;
    params.entry_size_bytes = block_size;
    params.poly_degree = 2048;
    params.n0 = n0;
    params.ciphertext_moduli = {134176769ULL, 134111233ULL};
    params.log_t = 8;
    params.variance = 10;
    params.conversion_gadget_log_base = 7;
    params.expansion_gadget_log_base = 6;
    params.trace_gadget_log_base = 5;
    params.response_modulus = 262144ULL;
    return params;
  }

  // Serializes parameters to protobuf.
  SerializedPirParams ToProto() const {
    SerializedPirParams proto;
    proto.set_num_entries(num_entries);
    proto.set_entry_size_bytes(entry_size_bytes);
    proto.set_poly_degree(poly_degree);
    proto.set_n0(n0);
    for (auto q_val : ciphertext_moduli) {
      proto.add_ciphertext_moduli(q_val);
    }
    proto.set_log_t(log_t);
    proto.set_variance(variance);
    proto.set_conversion_gadget_log_base(conversion_log_base());
    for (auto b_val : conversion_gadget_log_bs) {
      proto.add_conversion_gadget_log_bs(b_val);
    }
    proto.set_expansion_gadget_log_base(expansion_log_base());
    for (auto b_val : expansion_gadget_log_bs) {
      proto.add_expansion_gadget_log_bs(b_val);
    }
    proto.set_trace_gadget_log_base(trace_log_base());
    for (auto b_val : trace_gadget_log_bs) {
      proto.add_trace_gadget_log_bs(b_val);
    }
    proto.set_response_modulus(response_modulus);
    proto.set_num_query_coeffs(num_query_coeffs);
    return proto;
  }

  // Deserializes parameters from protobuf.
  static OnionPirParams FromProto(const SerializedPirParams& proto) {
    OnionPirParams params;
    params.num_entries = proto.num_entries();
    params.entry_size_bytes = proto.entry_size_bytes();
    params.poly_degree = proto.poly_degree();
    params.n0 = proto.n0();
    params.ciphertext_moduli.assign(proto.ciphertext_moduli().begin(),
                                    proto.ciphertext_moduli().end());
    if (params.ciphertext_moduli.empty()) {
      params.ciphertext_moduli = {134176769ULL, 134111233ULL};
    }
    params.log_t = proto.log_t();
    params.variance = proto.variance();

    params.conversion_gadget_log_base = proto.conversion_gadget_log_base() > 0
                                            ? proto.conversion_gadget_log_base()
                                            : 7;
    params.conversion_gadget_log_bs.assign(
        proto.conversion_gadget_log_bs().begin(),
        proto.conversion_gadget_log_bs().end());

    params.expansion_gadget_log_base = proto.expansion_gadget_log_base() > 0
                                           ? proto.expansion_gadget_log_base()
                                           : 5;
    params.expansion_gadget_log_bs.assign(
        proto.expansion_gadget_log_bs().begin(),
        proto.expansion_gadget_log_bs().end());

    params.trace_gadget_log_base =
        proto.trace_gadget_log_base() > 0 ? proto.trace_gadget_log_base() : 5;
    params.trace_gadget_log_bs.assign(proto.trace_gadget_log_bs().begin(),
                                      proto.trace_gadget_log_bs().end());

    params.response_modulus =
        proto.response_modulus() > 0 ? proto.response_modulus() : 262144ULL;
    params.num_query_coeffs = proto.num_query_coeffs();
    return params;
  }

  // Returns fast lightweight parameters for fast testing.
  static OnionPirParams CreateFastForTesting(size_t num_blocks,
                                             size_t block_size) {
    OnionPirParams params;
    params.num_entries = num_blocks;
    params.entry_size_bytes = block_size;
    params.poly_degree = 1024;
    params.n0 = 8;
    params.ciphertext_moduli = {134176769ULL, 134111233ULL};
    params.log_t = 8;
    params.variance = 10;
    params.conversion_gadget_log_base = 7;
    params.expansion_gadget_log_base = 5;
    params.trace_gadget_log_base = 5;
    params.response_modulus = 262144ULL;
    return params;
  }
};

// Computes the expected noise variance and remaining noise budget for OnionPIR
// operations under the independence heuristic (Gaussian approximation).
inline NoiseEstimate EstimateOnionPirNoise(const OnionPirParams& params,
                                           double confidence_factor) {
  NoiseEstimate est;
  est.confidence_factor = confidence_factor;

  double n = static_cast<double>(params.poly_degree);
  double sigma_sq = static_cast<double>(params.variance);
  double t = static_cast<double>(params.plaintext_modulus());
  double b_conv = static_cast<double>(static_cast<uint64_t>(1)
                                      << params.conversion_log_base());
  double ell_conv = static_cast<double>(params.conversion_ell());
  double b_exp = static_cast<double>(static_cast<uint64_t>(1)
                                     << params.expansion_log_base());
  double ell_exp = static_cast<double>(params.expansion_ell());
  double b_trace =
      static_cast<double>(static_cast<uint64_t>(1) << params.trace_log_base());
  double ell_trace = static_cast<double>(params.trace_ell());
  double n0 = static_cast<double>(params.n0);
  double d = static_cast<double>(params.num_dimensions());

  // 1. Fresh encryption noise variance:
  est.fresh_noise_variance = sigma_sq;

  // 2. Key-switching noise variance per automorphism / gadget multiplication:
  double v_digit_trace = (b_trace * b_trace - 1.0) / 12.0;
  double v_ks_trace = ell_trace * n * v_digit_trace * sigma_sq;

  double v_digit_exp = (b_exp * b_exp - 1.0) / 12.0;
  double v_ks_exp = ell_exp * n * v_digit_exp * sigma_sq;

  // Representative key-switching noise variance (trace)
  est.keyswitch_noise_variance = v_ks_trace;

  // 3. Unpacked query noise after Homomorphic Trace (L_trace levels) and Tree
  // Expansion (L_exp levels): Total levels = log2(N). Fresh noise is scaled by
  // N, contributing N^2 * V_fresh variance. Homomorphic trace: (N^2 - w_sent^2)
  // / 3 * V_ks_trace Tree expansion: (w_sent^2 - 1) / 3 * V_ks_exp
  double w_sent = static_cast<double>(params.query_w());
  double trace_ks_contrib = ((n * n - w_sent * w_sent) / 3.0) * v_ks_trace;
  double exp_ks_contrib = ((w_sent * w_sent - 1.0) / 3.0) * v_ks_exp;

  est.unpacked_query_noise_variance =
      n * n * est.fresh_noise_variance + trace_ks_contrib + exp_ks_contrib;

  // 4. Initial dimension 0 dot product noise:
  double v_pt = (t * t - 1.0) / 12.0;
  est.dim0_dot_product_noise_variance =
      n0 * n * v_pt * est.unpacked_query_noise_variance;

  // 5. Higher dimension (d - 1) folding via External Product:
  double v_digit_conv = (b_conv * b_conv - 1.0) / 12.0;
  double num_high_dims = (d > 1.0) ? (d - 1.0) : 0.0;
  double v_rgsw_row = est.unpacked_query_noise_variance +
                      2.0 * ell_conv * n * v_digit_conv * sigma_sq;
  double v_ext_prod = 2.0 * ell_conv * n * v_digit_conv * v_rgsw_row;
  est.higher_dims_noise_variance = num_high_dims * v_ext_prod;

  // Total server evaluation noise variance before modulus reduction
  est.server_eval_noise_variance =
      est.dim0_dot_product_noise_variance + est.higher_dims_noise_variance;

  // 6. Ciphertext modulus calculations:
  double q0 = static_cast<double>(params.primary_ciphertext_modulus());
  double total_q = 1.0;
  for (auto q_i : params.ciphertext_moduli) {
    total_q *= static_cast<double>(q_i);
  }
  if (params.ciphertext_moduli.empty()) {
    total_q = q0;
  }

  // 7. Modulus switching to response modulus p:
  double p = static_cast<double>(
      params.response_modulus > 0 ? params.response_modulus : 262144ULL);
  double q_dec = p;
  double scale = p / total_q;
  // With Gaussian secret distribution (sigma_sq = params.variance),
  // Rounding noise: V_round = (1 + N * sigma_sq) / 12
  double v_round = (1.0 + n * sigma_sq) / 12.0;
  est.final_noise_variance =
      (scale * scale) * est.server_eval_noise_variance + v_round;

  est.final_noise_std_dev = std::sqrt(est.final_noise_variance);
  est.max_noise_capacity = q_dec / (2.0 * t);

  double bound = confidence_factor * est.final_noise_std_dev;
  if (bound > 0.0 && est.max_noise_capacity > bound) {
    est.noise_budget_bits = std::log2(est.max_noise_capacity / bound);
    est.is_decryption_correct = true;
  } else if (bound > 0.0) {
    est.noise_budget_bits = -std::log2(bound / est.max_noise_capacity);
    est.is_decryption_correct = false;
  } else {
    est.noise_budget_bits = 0.0;
    est.is_decryption_correct = false;
  }

  return est;
}

}  // namespace oram_okvs::pir

#endif  // PRIVATE_MEMBERSHIP_RESEARCH_ORAM_OKVS_PIR_ONION_PIR_PARAMS_H_
