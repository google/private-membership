#ifndef PRIME_FIELD_H
#define PRIME_FIELD_H

#include <stdint.h>
#include <algorithm>

using namespace std;

class PrimeField {
public:
  PrimeField(uint64_t value, uint64_t field_size) : value_(value), field_size_(field_size) {}

  uint64_t GetRawValue() const {
    return value_;
  }

  uint64_t GetFieldSize() const {
    return field_size_;
  }

  bool IsZero() const {
    return value_ == 0;
  }

  void AddInPlace(const PrimeField& other) {
    value_ = (value_ + other.GetRawValue()) % field_size_;
  }

  PrimeField Add(const PrimeField& other) const {
    return PrimeField((value_ + other.GetRawValue()) % field_size_, field_size_);
  }

  void SubtractInPlace(const PrimeField& other) {
    AddInPlace(other.Negation());
  }

  PrimeField Subtract(const PrimeField& other) {
    return Add(other.Negation());
  }

  void MultInPlace(const PrimeField& other) {
    value_ = (static_cast<uint128_t>(value_) * other.GetRawValue()) % field_size_;
  }

  PrimeField Mult(const PrimeField& other) {
    return PrimeField((static_cast<uint128_t>(value_) * other.GetRawValue()) % field_size_, field_size_);
  }

  PrimeField MultInverse() const {
    uint64_t inv = BinPow(value_, field_size_ - 2, field_size_);
    return PrimeField(inv, field_size_);
  }

  PrimeField Negation() const {
    return PrimeField(field_size_ - value_, field_size_);
  }

private:
  uint64_t BinPow(uint64_t a, uint64_t b, uint64_t mod) const {
    uint64_t res = 1;
    while (b > 0) {
      if (b & 1) {
        res = (static_cast<uint128_t>(res) * a) % mod;
      }
      a = (static_cast<uint128_t>(a) * a) % mod;
      b >>= 1;
    }
    return res;
  }

  uint64_t value_;
  uint64_t field_size_;
};

#endif //PRIME_FIELD_H
