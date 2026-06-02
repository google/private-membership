#ifndef BAND_ROW_VECTOR_HASHER_H
#define BAND_ROW_VECTOR_HASHER_H

#include <string>

#include "band_row_vector.h"
#include "bin_allocator.h"

using namespace std;

class BandRowVectorHasher {
public:
  BandRowVectorHasher() {}

  BandRowVectorHasher(const BinAllocator &bin_allocator, int dimension, int band_width, uint64_t field_size,
                      uint64_t seed) :
          bin_allocator_(&bin_allocator), band_width_(band_width),
          field_size_(field_size), seed_(seed), aes_(oc::block(seed)) {}

  Band<uint256> Hash(const oc::block &key) const {
    int offset = ComputeBinOffset(key, band_width_);
    uint256 band(0, ComputeBand(key, band_width_));
    return Band<uint256>(offset, band);
  }

  BandRowVector HashToBandRowVector(oc::block key) const {
    int offset = ComputeBinOffset(key, band_width_);
    __uint128_t band = ComputeBand(key, band_width_);

    vector<uint64_t> band_vec;
    band_vec.reserve(band_width_);
    while (band != 0) {
      band_vec.push_back(band & 1);
      band >>= 1;
    }
    reverse(band_vec.begin(), band_vec.end());
    return BandRowVector(offset, band_vec, field_size_);
  }

private:
  int ComputeBinOffset(const oc::block &key, int band_width) const {
    oc::block hash = aes_.hashBlock(key);
    return static_cast<int>(hash.get<uint64_t>(0) % (bin_allocator_->BinSize() - band_width + 1));
  }

  __uint128_t ComputeBand(const oc::block &key, int band_width) const {
    if (band_width == 0) {
      return 0;
    }
    __uint128_t mask = band_width == 128 ? -1 : ((static_cast<__uint128_t>(1) << band_width) - 1);
    __uint128_t salted_key = aes_.hashBlock(key).get<__uint128_t>(0);
    return (salted_key & mask) | 1;
  }

  const BinAllocator *bin_allocator_;
  int dimension_;
  int band_width_;
  uint64_t field_size_;
  uint64_t seed_;
  oc::AES aes_;
};

#endif //BAND_ROW_VECTOR_HASHER_H
