// Copyright 2021 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "private_membership/rlwe/batch/cpp/encoding.h"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <google/protobuf/repeated_field.h>
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"
#include "private_membership/rlwe/batch/cpp/constants.h"
#include "private_membership/rlwe/batch/cpp/padding.h"
#include "private_membership/rlwe/batch/proto/shared.pb.h"
#include "shell_encryption/context.h"
#include "shell_encryption/galois_key.h"
#include "shell_encryption/integral_types.h"
#include "shell_encryption/montgomery.h"
#include "shell_encryption/polynomial.h"
#include "shell_encryption/prng/prng.h"
#include "shell_encryption/prng/single_thread_hkdf_prng.h"
#include "shell_encryption/serialization.pb.h"
#include "shell_encryption/transcription.h"

namespace private_membership {
namespace batch {

namespace {

absl::StatusOr<std::vector<ModularInt>> GetCoefficients(
    const std::vector<ModularInt::Int>& polynomial, const Context& context) {
  std::vector<ModularInt> transcribed_coeffs;
  transcribed_coeffs.reserve(polynomial.size());
  for (const ModularInt::Int& coeff : polynomial) {
    absl::StatusOr<ModularInt> coeff_modular_int =
        ModularInt::ImportInt(coeff, context.GetModulusParams());
    if (!coeff_modular_int.ok()) {
      return coeff_modular_int.status();
    }
    transcribed_coeffs.push_back(*std::move(coeff_modular_int));
  }
  return transcribed_coeffs;
}

absl::StatusOr<rlwe::Polynomial<ModularInt>> EncodeSinglePolynomial(
    const Context& context, const std::vector<uint8_t>& bytes) {
  absl::StatusOr<std::vector<ModularInt::Int>> transcribed =
      rlwe::TranscribeBits<uint8_t, ModularInt::Int>(
          bytes, /*input_bit_length=*/bytes.size() * 8,
          /*input_bits_per_int=*/8,
          /*output_bits_per_int=*/context.GetLogT());
  if (!transcribed.ok()) {
    return transcribed.status();
  }
  absl::StatusOr<std::vector<ModularInt>> transcribed_coeffs =
      GetCoefficients(*transcribed, context);
  if (!transcribed_coeffs.ok()) {
    return transcribed_coeffs.status();
  }
  return rlwe::Polynomial<ModularInt>::ConvertToNtt(
      *std::move(transcribed_coeffs), context.GetNttParams(),
      context.GetModulusParams());
}

}  // namespace

absl::StatusOr<std::vector<rlwe::Polynomial<ModularInt>>> EncodeBytes(
    const Context& context, absl::string_view bytes) {
  std::vector<rlwe::Polynomial<ModularInt>> encoded_polynomials;
  if (bytes.empty()) {
    return encoded_polynomials;
  }

  int bytes_per_polynomial = (context.GetN() * context.GetLogT()) / 8;
  if (bytes_per_polynomial <= 0) {
    return absl::InvalidArgumentError("Invalid parameters used for encoding.");
  }
  std::string length_prepended_bucket = PrependLength(bytes);
  // Ceiling of length_prepended_bucket.size() / bytes_per_polynomial.
  int number_of_ciphertexts =
      (length_prepended_bucket.size() + bytes_per_polynomial - 1) /
      bytes_per_polynomial;

  encoded_polynomials.reserve(number_of_ciphertexts);

  std::vector<uint8_t> padded_plaintext(bytes_per_polynomial, 0);
  for (int i = 0; i < number_of_ciphertexts; ++i) {
    int start = i * bytes_per_polynomial;
    int end = (i + 1) * bytes_per_polynomial;

    // On only the last iteration, the plaintext may be shorter.
    if (end > length_prepended_bucket.size()) {
      end = length_prepended_bucket.size();
    }

    std::copy(length_prepended_bucket.begin() + start,
              length_prepended_bucket.begin() + end, padded_plaintext.begin());

    absl::StatusOr<rlwe::Polynomial<ModularInt>> encoded_polynomial =
        EncodeSinglePolynomial(context, padded_plaintext);
    encoded_polynomials.push_back(*std::move(encoded_polynomial));
  }
  return encoded_polynomials;
}

absl::StatusOr<std::string> DecodePolynomial(
    const Context& context,
    const std::vector<rlwe::Polynomial<ModularInt>>& polynomials) {
  if (polynomials.empty()) {
    return "";
  }

  std::string padded_result;
  padded_result.reserve(polynomials.size() * context.GetN() *
                        context.GetLogT());
  for (const auto& polynomial : polynomials) {
    std::vector<ModularInt> poly_coeffs = polynomial.InverseNtt(
        context.GetNttParams(), context.GetModulusParams());
    std::vector<ModularInt::Int> int_coeffs;
    int_coeffs.reserve(polynomial.Len());
    for (const ModularInt& coeff : poly_coeffs) {
      int_coeffs.push_back(coeff.ExportInt(context.GetModulusParams()));
    }

    // Transcribe the plaintext polynomial back to bytes.
    absl::StatusOr<std::vector<uint8_t>> bytes =
        rlwe::TranscribeBits<ModularInt::Int, uint8_t>(
            int_coeffs,
            /*input_bit_length=*/int_coeffs.size() * context.GetLogT(),
            /*input_bits_per_int=*/context.GetLogT(),
            /*output_bits_per_int=*/8);
    if (!bytes.ok()) {
      return bytes.status();
    }
    padded_result.append(reinterpret_cast<const char*>(bytes->data()),
                         bytes->size());
  }

  return Unpad(padded_result);
}

absl::StatusOr<rlwe::GaloisKey<ModularInt>> DeserializePublicGaloisKey(
    const rlwe::SerializedGaloisKey& serialized_key, const Context& context) {
  return rlwe::GaloisKey<ModularInt>::Deserialize(
      serialized_key, context.GetModulusParams(), context.GetNttParams());
}

}  // namespace batch
}  // namespace private_membership
