#ifndef BANDOKVS_BAND_H_
#define BANDOKVS_BAND_H_

#include "cryptoTools/Crypto/PRNG.h"
#include <iostream>

namespace band_okvs {

template<typename T>
class Band {
 public:
  Band() : band_start_(0), band_(0), idx_(0) {}

  Band(int band_start, T band, int idx)
      : band_start_(band_start), band_(band), idx_(idx) {}

  void Set(int band_start, T band, int idx) {
    band_start_ = band_start;
    band_ = band;
    idx_ = idx;
  }

  inline int BandStart() const {
    return band_start_;
  }

  inline T RawBand() const {
    return band_;
  }

  inline int Index() const {
    return idx_;
  }

  int GetBit(int index) const {
    if (index < band_start_ || index >= band_start_ + 8 * sizeof(T)) {
      return 0;
    }
    return band_.GetBit(index - band_start_);
  }

  size_t band_start_;
  T band_;
  size_t idx_;
};

template<typename T, typename V>
class BandAndValue {
 public:
  BandAndValue() : band_start_(0), value_(0), band_(0) {}

  BandAndValue(const Band<T>& band, V value)
      : band_start_(band.BandStart()), value_(value), band_(band.RawBand()) {}

  BandAndValue(int band_start, T band, V value)
      : band_start_(band_start), value_(value), band_(band) {}

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

  size_t band_start_;
  V value_;
  T band_;
};

struct lessthan {
  int B_;

  lessthan() = default;

  lessthan(int B) : B_(B) {}

  template<typename T>
  inline int operator()(const Band<T>& x, const Band<T>& y) const {
    return x.BandStart() < y.BandStart();
  }

  template<typename T, typename V>
  inline int operator()(const BandAndValue<T, V>& x,
                        const BandAndValue<T, V>& y) const {
    return x.BandStart() < y.BandStart();
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

  template<typename T>
  inline int operator()(const Band<T>& x, const unsigned int offset) const {
    return x.BandStart() >> offset;
  }

  template<typename T, typename V>
  inline int operator()(const BandAndValue<T, V>& x,
                        const unsigned int offset) const {
    return x.BandStart() >> offset;
  }

  inline int operator()(const oc::block& x, const unsigned int offset) const {
    return (x.get<uint32_t>(0) % B_) >> offset;
  }
};

}  // namespace band_okvs

#endif //BANDOKVS_BAND_H_
