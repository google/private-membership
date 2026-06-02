#ifndef PLAINTEXT_VALUE_H
#define PLAINTEXT_VALUE_H

#include <cassert>
#include <vector>

#include "seal/seal.h"
#include "seal/util/polyarithsmallmod.h"

#include "int_utils.h"
#include "prime_field.h"

using namespace std;

class PlaintextValue {
public:
  PlaintextValue() = default;

  explicit PlaintextValue(uint64_t field_size) : modulus_(field_size) {}

  PlaintextValue(const vector<uint64_t> &coeffs, uint64_t field_size) : plaintext_(coeffs.size()),
                                                                        modulus_(field_size) {
    copy(coeffs.begin(), coeffs.end(), &plaintext_[0]);
  }

  vector<uint64_t> GetValue() const {
    vector<uint64_t> res;
    for (int i = 0; i < plaintext_.coeff_count(); ++i) {
      res.push_back(plaintext_[i]);
    }
    return res;
  }

  uint64_t Modulus() const {
    return modulus_.value();
  }

  inline void ScalarMultSubtractVecInPlace(uint64_t y, __m512i y_batch, uint64_t y_barrett, __m512i y_barrett_batch,
                                           const PlaintextValue &other, __m512i q_batch) {
    ScalarMulSubVecMod(plaintext_.data(), other.plaintext_.data(), y, y_batch, y_barrett, y_barrett_batch,
                       modulus_.value(), q_batch,
                       plaintext_.coeff_count());
  }

  inline void SubtractVectorScalarMultInPlace(const PlaintextValue &other, uint64_t z) {
    if (other.plaintext_.coeff_count() == 0) {
      return;
    }
    __m512i z_batch = _mm512_set1_epi64(z);
    uint64_t z_barrett = (static_cast<__uint128_t>(z) << 64) / modulus_.value();
    __m512i z_barrett_batch = _mm512_set1_epi64(z_barrett);
    __m512i q_batch = _mm512_set1_epi64(modulus_.value());
    SubVecScalarMulMod(plaintext_.data(), other.plaintext_.data(), z, z_batch, z_barrett, z_barrett_batch,
                       modulus_.value(), q_batch,
                       plaintext_.coeff_count());
  }

  inline void ScalarMultInPlace(uint64_t y) {
    if (plaintext_.is_zero()) {
      return;
    }
    __m512i y_batch = _mm512_set1_epi64(y);
    uint64_t y_barrett = (static_cast<__uint128_t>(y) << 64) / modulus_.value();
    __m512i y_barrett_batch = _mm512_set1_epi64(y_barrett);
    __m512i q_batch = _mm512_set1_epi64(modulus_.value());
    ScalarMulMod(plaintext_.data(), y, y_batch, y_barrett, y_barrett_batch, modulus_.value(), q_batch,
                 plaintext_.coeff_count());
    //seal::util::multiply_poly_scalar_coeffmod(
    //        plaintext_.data(), plaintext_.coeff_count(), scalar, modulus_, plaintext_.data());
  }

  inline void AddInPlace(const PlaintextValue &other) {
    if (other.plaintext_.coeff_count() == 0) {
      return;
    }
    if (plaintext_.coeff_count() == 0) {
      plaintext_ = other.plaintext_;
      return;
    }
    assert(plaintext_.coeff_count() == other.plaintext_.coeff_count());
    seal::util::add_poly_poly_coeffmod(
            plaintext_.data(), other.plaintext_.data(), plaintext_.coeff_count(), modulus_, plaintext_.data());
  }

  inline void SubtractInPlace(const PlaintextValue &other) {
    if (other.plaintext_.coeff_count() == 0) {
      return;
    }
    if (plaintext_.coeff_count() == 0) {
      plaintext_ = other.plaintext_;
      return;
    }
    assert(plaintext_.coeff_count() == other.plaintext_.coeff_count());
    seal::util::sub_poly_poly_coeffmod(
            plaintext_.data(), other.plaintext_.data(), plaintext_.coeff_count(), modulus_, plaintext_.data());
  }

  bool operator==(const PlaintextValue &other) {
    return modulus_ == other.modulus_ &&
           plaintext_ == other.plaintext_;
  }

  bool operator!=(const PlaintextValue &other) {
    return !(*this == other);
  }

  bool IsZero() const {
    return plaintext_.is_zero();
  }

  int CoeffCount() {
    return plaintext_.coeff_count();
  }

  void Print() {
    for (int i = 0; i < plaintext_.coeff_count(); ++i) {
      cout << plaintext_[i] << " ";
    }
    cout << endl;
  }

private:
  seal::Modulus modulus_;
  seal::Plaintext plaintext_;
};

uint64_t CoefficientsPerElement(uint32_t logtp, uint64_t element_size) {
  return ceil(8 * element_size / (double) logtp);
}

vector<uint64_t> ValueToCoeffs(uint64_t t, uint64_t logt, const string &value) {
  uint64_t size_out = CoefficientsPerElement(logt, value.size());
  vector<uint64_t> output(size_out);
  uint32_t room = logt;
  uint64_t *target = &output[0];
  for (int i = 0; i < value.size(); ++i) {
    uint8_t src = value[i];
    uint32_t rest = 8;
    while (rest > 0) {
      if (room == 0) {
        target++;
        room = logt;
      }
      uint32_t shift = rest;
      if (room < rest) {
        shift = room;
      }
      *target = *target << shift;
      *target = *target | (src >> (8 - shift));
      src = src << shift;
      room -= shift;
      rest -= shift;
    }
  }
  return output;
}

seal::Plaintext CoeffsToPlaintext(const vector<uint64_t> &coeffs, int poly_deg, int coeff_length) {
  assert(coeffs.size() <= poly_deg);
  seal::Plaintext plaintext(poly_deg);
  vector<uint64_t> raw_coeffs(poly_deg);
  for (int i = 0; i < coeffs.size(); ++i) {
    raw_coeffs[i] = coeffs[i];
  }
  for (int i = coeff_length; i < poly_deg; ++i) {
    raw_coeffs[i] = 1;
  }
  seal::util::set_uint_uint(raw_coeffs.data(), raw_coeffs.size(), plaintext.data());
  return plaintext;
}

vector<uint64_t> PlaintextToCoeffs(uint64_t t, const seal::Plaintext &plaintext) {
  vector<uint64_t> res;
  for (uint32_t i = 0; i < plaintext.coeff_count(); i++) {
    res.push_back(plaintext[i]);
  }
  return res;
}

vector<PlaintextValue> ToPlaintextValues(uint64_t t, uint64_t logt, const vector<string> &raw_values) {
  vector<PlaintextValue> res;
  res.reserve(raw_values.size());
  for (const string &raw_value: raw_values) {
    vector<uint64_t> coeffs = ValueToCoeffs(t, logt, raw_value);
    res.push_back(PlaintextValue(coeffs, t));
  }
  return res;
}

string CoeffsToValue(uint64_t t, uint64_t logt, int offset, const vector<uint64_t> &coeffs, uint64_t element_size) {
  uint64_t coeff_size = t;
  vector<uint64_t> raw_coeffs;
  raw_coeffs.reserve(coeffs.size());
  for (int i = 0; i < coeffs.size(); ++i) {
    raw_coeffs.push_back(coeffs[i] % coeff_size);
  }

  uint32_t room = 8;
  uint32_t j = 0;
  string output(element_size, ' ');
  char *target = &output[0];

  for (uint32_t i = offset; i < coeffs.size(); ++i) {
    uint64_t src = raw_coeffs[i];
    uint32_t rest = logt;
    while (rest && j < element_size) {
      uint32_t shift = rest;
      if (room < rest) {
        shift = room;
      }
      target[j] = target[j] << shift;
      target[j] = target[j] | (src >> (logt - shift));
      src = src << shift;
      room -= shift;
      rest -= shift;
      if (room == 0) {
        j++;
        room = 8;
      }
    }
  }
  return output;
}

void VectorToPlaintext(const vector<uint64_t> &coeffs, seal::Plaintext &plain) {
  uint32_t coeff_count = coeffs.size();
  plain.resize(coeff_count);
  seal::util::set_uint_uint(coeffs.data(), coeff_count, plain.data());
}

#endif //PLAINTEXT_VALUE_H
