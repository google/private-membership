#ifndef BAND_ROW_VECTOR_H
#define BAND_ROW_VECTOR_H

#include <assert.h>
#include <deque>
#include <vector>

#include "seal/seal.h"
#include "seal/util/polyarithsmallmod.h"
#include "cryptoTools/Common/block.h"

#include "prime_field.h"
#include "int_utils.h"

namespace oc = osuCrypto;

class BandRowVector {
public:
  BandRowVector() = default;

  BandRowVector(uint64_t field_size) : band_length_(0), offset_(0),
    modulus_(field_size) {
  }

  BandRowVector(int offset, const vector<uint64_t>& band, uint64_t field_size) :
    band_length_(band.size()), offset_(offset),
    //band_(band.size()),
    modulus_(field_size) {
    //band_.resize(band.size());
    copy(band.begin(), band.end(), band_);
  }

  inline uint64_t GetFirstElement() const {
    return band_[0];
  }

  uint64_t* RawBand() {
    return band_;
  }

  void Reduce(const BandRowVector& other, uint64_t y, __m512i y_batch, uint64_t y_barrett,  __m512i y_barrett_batch,
              const __m512i& q_batch) {
    assert(offset_ == other.offset_);

    int band_length = max(band_length_, other.band_length_);
    ScalarMulSubVecMod(band_, other.band_, y, y_batch, y_barrett, y_barrett_batch,
                           modulus_.value(),
                           q_batch, band_length);

    int leading_zero_cnt = 0;
    for (int i = 0; i < band_length; ++i) {
      //int j = band_offset_ + i;
      if (band_[i] != 0) {
        copy(band_ + leading_zero_cnt, band_ + band_length, band_);
        band_length_ = band_length - i;
        fill(band_ + band_length_, band_ + BAND_CAPACITY, 0);
        break;
      }
      ++leading_zero_cnt;
    }
    if (leading_zero_cnt == band_length) {
      band_length_ = 0;
    }
    offset_ += leading_zero_cnt;
  }


  int GetLength() const {
    return band_length_;
  }

  int GetOffset() const {
    return offset_;
  }

  uint64_t GetRowElement(int index) const {
    if (index < offset_ || index >= offset_ + band_length_) {
      return 0;
    }
    return band_[index - offset_];
  }

  bool IsZero() const {
    return band_length_ == 0;
  }

  void Print() {
    cout << offset_ << " ";
    for (int i = 0; i < band_length_; ++i) {
      cout << band_[i] << " ";
    }
    cout << endl;
  }

  uint64_t Modulus() const {
    return modulus_.value();
  }

private:
  static constexpr int BAND_CAPACITY = 500;
  int band_length_ = 0;
  int offset_ = 0;

  uint64_t band_[BAND_CAPACITY] = {0};
  seal::Modulus modulus_;
};

template <typename T>
class Band {
 public:
  Band() : band_start_(0), band_(0) {}

  Band(int band_start, T band)
      : band_start_(band_start), band_(band) {}

  void Set(int band_start, T band) {
    band_start_ = band_start;
    band_ = band;
  }

  int BandStart() const {
    return band_start_;
  }

  T RawBand() const {
    return band_;
  }

 private:
  int band_start_;
  T band_;
  oc::block value_;
};


template <typename T, typename V>
class BandAndValue {
 public:
  BandAndValue() : band_start_(0), band_(0), value_(0) {}

  BandAndValue(const Band<T>& band, V value)
      : band_start_(band.BandStart()), band_(band.RawBand()), value_(value) {}

  BandAndValue(int band_start, T band, V value)
      : band_start_(band_start), band_(band), value_(value) {}

  void Set(int band_start, T band, V value) {
    band_start_ = band_start;
    band_ = band;
    value_ = value;
  }

  int BandStart() const {
    return band_start_;
  }

  T RawBand() const {
    return band_;
  }

  V RawValue() const {
    return value_;
  }

 private:
  int band_start_;
  T band_;
  V value_;
};

template<typename V>
class HashAndValue {
 public:
  HashAndValue() = default;

  HashAndValue(oc::block hash, V value) : hash_(hash), value_(value) {}

  oc::block RawHash() const {
    return hash_;
  }

  V RawValue() const {
    return value_;
  }

 private:
  oc::block hash_;
  V value_;
};

struct lessthan {
  int B_;

  lessthan() = default;

  lessthan(int B) : B_(B) {}

  template <typename T>
  inline int operator()(const Band<T>& x, const Band<T>& y) const {
    return x.BandStart() < y.BandStart();
  }

  template <typename T, typename V>
  inline int operator()(const BandAndValue<T, V>& x, const BandAndValue<T, V>& y) const {
    return x.BandStart() < y.BandStart();
  }

  template <typename V>
  inline int operator()(const HashAndValue<V>& x, const HashAndValue<V>& y) const {
    oc::block hash1 = x.RawHash();
    oc::block hash2 = y.RawHash();
    uint32_t idx1 = hash1.get<uint32_t>(0) % B_;
    uint32_t idx2 = hash2.get<uint32_t>(0) % B_;
    return idx1 < idx2;
  }

  inline int operator()(const oc::block& x, const oc::block& y) const {
    uint32_t idx1 = x.get<uint32_t>(0) % B_;
    uint32_t idx2 = y.get<uint32_t>(0) % B_;
    return idx1 < idx2;
  }
};

struct rightshift {
  int B_;

  rightshift() = default;

  rightshift(int B) : B_(B) {}

  template <typename T>
  inline int operator()(const Band<T>& x, const unsigned int offset) const {
    return x.BandStart() >> offset;
  }

  template <typename T, typename V>
  inline int operator()(const BandAndValue<T, V>& x, const unsigned int offset) const {
    return x.BandStart() >> offset;
  }

  template <typename V>
  inline int operator()(const HashAndValue<V>& x, const unsigned int offset) const {
    oc::block hash = x.RawHash();
    return (hash.get<uint32_t>(0) % B_) >> offset;
  }

  inline int operator()(const oc::block& x, const unsigned int offset) const {
    return (x.get<uint32_t>(0) % B_) >> offset;
  }
};

class uint256;

constexpr inline uint256 MakeUint256(__uint128_t high, __uint128_t low);
constexpr inline uint256 operator|(uint256 lhs, uint256 rhs);
constexpr inline uint256 operator&(uint256 lhs, uint256 rhs);
constexpr inline __uint128_t operator&(uint256 lhs, int rhs);
constexpr inline uint256 operator^(uint256 lhs, uint256 rhs);
constexpr inline uint256 operator<<(uint256 lhs, int amount);

constexpr inline uint256 operator>>(uint256 lhs, int amount);
constexpr inline bool operator==(uint256 lhs, uint256 rhs);
constexpr inline bool operator!=(uint256 lhs, uint256 rhs);

class uint256 {
 public:
  uint256() = default;

  constexpr uint256(__uint128_t high, __uint128_t low)
      : low_(low), high_(high) {}

  constexpr uint256(int v) : low_(static_cast<__uint128_t>(v)), high_(0) {}

  /*constexpr uint256(__uint128_t high, __uint128_t low)
      : high_(high), low_(low) {}

  constexpr uint256(int v) : high_(0), low_(static_cast<__uint128_t>(v)) {}*/

  inline uint256& operator=(int v) { return *this = uint256(v); }

  constexpr explicit operator int() const { return static_cast<int>(low_); }

  inline uint256& operator&=(uint256 other) {
    *this = *this & other;
    return *this;
  }

  inline uint256& operator|=(uint256 other) {
    *this = *this | other;
    return *this;
  }

  inline uint256& operator^=(uint256 other) {
    //_mm256_store_epi64(this, _mm256_xor_epi64(_mm256_load_epi64(this), _mm256_load_epi64(&other)));
    *this = *this ^ other;
    return *this;
  }

  inline uint256& operator<<=(int amount) {
    *this = *this << amount;
    return *this;
  }

  inline uint256& operator>>=(int amount) {
    *this = *this >> amount;
    return *this;
  }

  friend constexpr __uint128_t Uint256Low128(uint256 v) { return v.low_; }

  friend constexpr __uint128_t Uint256High128(uint256 v) { return v.high_; }

 private:
  //__uint128_t high_ = 0;
  __uint128_t low_ = 0;
  __uint128_t high_ = 0;
};

constexpr inline uint256 MakeUint256(__uint128_t high, __uint128_t low) {
  return uint256(high, low);
}

constexpr inline uint256 operator|(uint256 lhs, uint256 rhs) {
  return MakeUint256(Uint256High128(lhs) | Uint256High128(rhs),
                     Uint256Low128(lhs) | Uint256Low128(rhs));
}

constexpr inline uint256 operator&(uint256 lhs, uint256 rhs) {
  return MakeUint256(Uint256High128(lhs) & Uint256High128(rhs),
                     Uint256Low128(lhs) & Uint256Low128(rhs));
}

constexpr inline __uint128_t operator&(uint256 lhs, int rhs) {
  return Uint256Low128(lhs) & rhs;
}

constexpr inline __uint128_t operator&(uint256 lhs, uint64_t rhs) {
  return Uint256Low128(lhs) & rhs;
}

constexpr inline __uint128_t operator&(uint256 lhs, __uint128_t rhs) {
  return Uint256Low128(lhs) & rhs;
}

constexpr inline uint256 operator^(uint256 lhs, uint256 rhs) {
  return MakeUint256(Uint256High128(lhs) ^ Uint256High128(rhs),
                     Uint256Low128(lhs) ^ Uint256Low128(rhs));
}

constexpr inline uint256 operator<<(uint256 lhs, int amount) {
  return amount >= 128 ? MakeUint256(Uint256Low128(lhs) << (amount - 128), 0)
         : amount == 0 ? lhs
                       : MakeUint256((Uint256High128(lhs) << amount) |
                                         (Uint256Low128(lhs) >> (128 - amount)),
                                     Uint256Low128(lhs) << amount);
}

constexpr inline uint256 operator>>(uint256 lhs, int amount) {
  return amount >= 128 ? MakeUint256(0, Uint256High128(lhs) >> (amount - 128))
         : amount == 0
             ? lhs
             : MakeUint256(Uint256High128(lhs) >> amount,
                           (Uint256Low128(lhs) >> amount) |
                               (Uint256High128(lhs) << (128 - amount)));
}

constexpr inline bool operator==(uint256 lhs, uint256 rhs) {
  return Uint256Low128(lhs) == Uint256Low128(rhs) &&
         Uint256High128(lhs) == Uint256High128(rhs);
}

constexpr inline bool operator!=(uint256 lhs, uint256 rhs) {
  return Uint256Low128(lhs) != Uint256Low128(rhs) ||
         Uint256High128(lhs) != Uint256High128(rhs);
}

constexpr inline bool operator==(uint256 lhs, int rhs) {
  return Uint256Low128(lhs) == rhs && Uint256High128(lhs) == 0;
}

constexpr inline bool operator!=(uint256 lhs, int rhs) {
  return Uint256Low128(lhs) != rhs || Uint256High128(lhs) != 0;
}

constexpr inline int ctz(uint256 v) {
  uint64_t a = static_cast<uint64_t>(Uint256Low128(v));
  if (a != 0) {
    return __builtin_ctzll(a);
  }
  uint64_t b = static_cast<uint64_t>(Uint256Low128(v) >> 64);
  if (b != 0) {
    return __builtin_ctzll(b) + 64;
  }
  uint64_t c = static_cast<uint64_t>(Uint256High128(v));
  if (c != 0) {
    return __builtin_ctzll(b) + 128;
  }
  uint64_t d = static_cast<uint64_t>(Uint256High128(v) >> 64);
  return __builtin_ctzll(d) + 192;
}

#endif //BAND_ROW_VECTOR_H