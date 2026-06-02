#ifndef LINEAR_SYSTEM_ROW_H
#define LINEAR_SYSTEM_ROW_H

#include <assert.h>
#include <stdint.h>
#include <vector>
#include <algorithm>

#include "band_row_vector.h"
#include "prime_field.h"
#include "utils.h"

using namespace std;

class TestValue {
public:
  TestValue(uint64_t field_size) : value_(0, field_size) {}

  explicit TestValue(const PrimeField& value) : value_(value) {}

  PrimeField GetValue() const {
    return value_;
  }

  void ScalarMultInPlace(const PrimeField& scalar) {
    value_.MultInPlace(scalar);
  }

  void AddInPlace(const TestValue& value) {
    value_.AddInPlace(value.value_);
  }

  void SubtractInPlace(const TestValue& value) {
    value_.SubtractInPlace(value.value_);
  }

private:
  PrimeField value_;
};

template <typename ValueType>
class LinearSystemRow {
public:
  LinearSystemRow(const BandRowVector& band_row_vector, const ValueType& value) :
    band_row_vector_(band_row_vector), ref_value_(&value), value_(value.Modulus()) {}

  uint64_t GetFirstNonZeroElement() const {
    return band_row_vector_.GetRowElement(band_row_vector_.GetOffset());
  }

  void FinalizeValue() {
    value_ = *ref_value_;
  }

  void Reduce(const LinearSystemRow& other, bool reduce_value = true) {
    uint64_t scalar = MultInverse(GetFirstNonZeroElement(), band_row_vector_.Modulus());
    scalar = (static_cast<uint128_t>(scalar) * other.GetFirstNonZeroElement()) % band_row_vector_.Modulus();
    if (reduce_value) {
      value_.ScalarMultInPlace(scalar);
      value_.SubtractInPlace(other.value_);
    }
  }

  int GetOffset() const {
    return band_row_vector_.GetOffset();
  }

  bool IsZero() const {
    return band_row_vector_.IsZero();
  }

  uint64_t GetRowElement(int index) const {
    return band_row_vector_.GetRowElement(index);
  }

  ValueType GetValue() const {
    return value_;
  }

  void Print() {
    band_row_vector_.Print();
  }

private:
  BandRowVector band_row_vector_;
  const ValueType* ref_value_;
  ValueType value_;
};

#endif //LINEAR_SYSTEM_ROW_H
