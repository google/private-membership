#ifndef BANDOKVS_BAND_OKVS_H_
#define BANDOKVS_BAND_OKVS_H_

#define OC_ENALBE_AESNI

#include "cryptoTools/Crypto/PRNG.h"
#include <immintrin.h>

#include <cassert>
#include <string>
#include <vector>

#include "band.h"
#include "uint.h"
#include "utils.h"
#include <boost/sort/spreadsort/spreadsort.hpp>

namespace band_okvs {

namespace {

template<typename T>
inline T CreateBand(oc::span<oc::block> blocks,
                    __uint128_t high_mask) {
  return T();
}

template<>
inline uint<1> CreateBand(oc::span<oc::block> blocks,
                          __uint128_t high_mask) {
  return uint < 1 > (((blocks[0]).get<__uint128_t>(0) & high_mask) | 1);
}

template<>
inline uint<2> CreateBand(oc::span<oc::block> blocks,
                          __uint128_t high_mask) {
  return uint < 2 > ((blocks[0]).get<__uint128_t>(0) | 1,
      (blocks[1]).get<__uint128_t>(0) & high_mask);
}

template<>
inline uint<3> CreateBand(oc::span<oc::block> blocks,
                          __uint128_t high_mask) {
  return uint < 3 > ((blocks[0]).get<__uint128_t>(0) | 1,
      (blocks[1]).get<__uint128_t>(0),
      (blocks[2]).get<__uint128_t>(0) & high_mask);
}

template<>
inline uint<4> CreateBand(oc::span<oc::block> blocks,
                          __uint128_t high_mask) {
  return uint < 4 > ((blocks[0]).get<__uint128_t>(0) | 1,
      (blocks[1]).get<__uint128_t>(0),
      (blocks[2]).get<__uint128_t>(0),
      (blocks[3]).get<__uint128_t>(0) & high_mask);
}

template<>
inline uint<5> CreateBand(oc::span<oc::block> blocks,
                          __uint128_t high_mask) {
  return uint < 5 > ((blocks[0]).get<__uint128_t>(0) | 1,
      (blocks[1]).get<__uint128_t>(0),
      (blocks[2]).get<__uint128_t>(0),
      (blocks[3]).get<__uint128_t>(0),
      (blocks[4]).get<__uint128_t>(0) & high_mask);
}

template<>
inline uint<6> CreateBand(oc::span<oc::block> blocks,
                          __uint128_t high_mask) {
  return uint < 6 > ((blocks[0]).get<__uint128_t>(0) | 1,
      (blocks[1]).get<__uint128_t>(0),
      (blocks[2]).get<__uint128_t>(0),
      (blocks[3]).get<__uint128_t>(0),
      (blocks[4]).get<__uint128_t>(0),
      (blocks[5]).get<__uint128_t>(0) & high_mask);
}

template<typename T>
inline void GenBands(int n, const oc::block* keys,
                     int okvs_length, int band_length,
                     Band<T>* out) {
  int num_blocks = (band_length + 127) / 128;
  std::vector<oc::block> expanded_keys(num_blocks);
  std::vector<oc::block> blocks(num_blocks);
  std::vector<oc::block> xor_blocks(num_blocks);
  for (int i = 0; i < num_blocks; i++) {
    xor_blocks[i] = oc::block(i);
  }

  __uint128_t
      high_mask = band_length % 128 == 0 ? static_cast<__uint128_t>(-1) :
                  (static_cast<__uint128_t>(1) << (band_length % 128)) - 1;
  uint32_t divisor = okvs_length - band_length + 1;

  oc::AES aes(oc::ZeroBlock);
#pragma GCC unroll 8
  for (int i = 0; i < n; ++i) {
    oc::block block = aes.hashBlock(keys[i]);
    uint32_t p = block.get<uint32_t>(0);
    uint32_t start_pos = p % divisor;

    for (int k = 0; k < blocks.size(); k++) {
      expanded_keys[k] = block ^ xor_blocks[k];
    }
    aes.hashBlocks(expanded_keys, blocks);
    out[i].band_start_ = start_pos;
    out[i].band_ = CreateBand<T>(blocks, high_mask);
    out[i].idx_ = i;
  }
}

template<typename T, typename V>
inline void GenBandsAndValues(int n, const oc::block* keys,
                              const V* values,
                              int okvs_length, int band_length,
                              BandAndValue<T, V>* out) {
  int num_blocks = (band_length + 127) / 128;
  std::vector<oc::block> expanded_keys(num_blocks);
  std::vector<oc::block> blocks(num_blocks);
  std::vector<oc::block> xor_blocks(num_blocks);
  for (int i = 0; i < num_blocks; i++) {
    xor_blocks[i] = oc::block(i);
  }

  __uint128_t
      high_mask = band_length % 128 == 0 ? static_cast<__uint128_t>(-1) :
                  (static_cast<__uint128_t>(1) << (band_length % 128)) - 1;
  uint32_t divisor = okvs_length - band_length + 1;
  oc::AES aes(oc::ZeroBlock);
#pragma GCC unroll 16
  for (int i = 0; i < n; ++i) {

    oc::block block = aes.hashBlock(keys[i]);
    uint32_t p = block.get<uint32_t>(0);
    uint32_t start_pos = p % divisor;

    for (int k = 0; k < num_blocks; k++) {
      expanded_keys[k] = block ^ xor_blocks[k];
    }
    aes.hashBlocks(expanded_keys, blocks);
    out[i].band_start_ = start_pos;
    out[i].value_ = values[i];
    out[i].band_ = CreateBand<T>(blocks, high_mask);
  }
}

}

std::vector<__uint128_t> GenRandomValues(int n, uint64_t seed1 = 0,
                                         uint64_t seed2 = 0);
std::vector<oc::block> GenRandomValuesBlocks(int n, uint64_t seed1 = 0,
                                             uint64_t seed2 = 0);

void GenRandomValuesBlocks(oc::span<oc::block> out, uint64_t seed1 = 0,
                           uint64_t seed2 = 0);

class BandOkvs {
 public:
  void Init(int n,
            int m,
            int band_length,
            oc::block prng_seed = oc::ZeroBlock) {
    num_eqns_ = n;
    num_vars_ = m;
    band_length_ = band_length;
    prng_seed_ = prng_seed;
  }

  // Add additional space
  int Size() const {
    return num_vars_ + 64;
  }

  int NumVars() const {
    return num_vars_;
  }

  int NumEqns() const {
    return num_eqns_;
  }

  template<typename V>
  bool Encode(const oc::block* keys, const V* values, V* out) {
    if (band_length_ <= 128) {
      return EncodeHelper<uint<1>, V>(keys, values, out);
    } else if (band_length_ <= 256) {
      return EncodeHelper<uint<2>, V>(keys, values, out);
    } else if (band_length_ <= 384) {
      return EncodeHelper<uint<3>, V>(keys, values, out);
    } else if (band_length_ <= 512) {
      return EncodeHelper<uint<4>, V>(keys, values, out);
    } else if (band_length_ <= 640) {
      return EncodeHelper<uint<5>, V>(keys, values, out);
    } else if (band_length_ <= 768) {
      return EncodeHelper<uint<6>, V>(keys, values, out);
    }
    std::cout << "Encode failed! Invalid band length: " << band_length_ << "\n";
    return false;
  }

  template<typename V>
  void Decode(const oc::block* keys, const V* okvs, V* out) {
    if (band_length_ <= 128) {
      return DecodeHelper<uint<1>, V>(keys, okvs, out);
    } else if (band_length_ <= 256) {
      return DecodeHelper<uint<2>, V>(keys, okvs, out);
    } else if (band_length_ <= 384) {
      return DecodeHelper<uint<3>, V>(keys, okvs, out);
    } else if (band_length_ <= 512) {
      return DecodeHelper<uint<4>, V>(keys, okvs, out);
    } else if (band_length_ <= 640) {
      return DecodeHelper<uint<5>, V>(keys, okvs, out);
    } else if (band_length_ <= 768) {
      return DecodeHelper<uint<6>, V>(keys, okvs, out);
    }
    std::cout << "Decode failed! Invalid band length: " << band_length_ << "\n";
  }

 private:
  template<typename T, typename V>
  inline bool ReduceToRowEchelonHelper(const BandAndValue<T, V>* bands,
                                       T* reduced_matrix,
                                       V* reduced_values) {
    int n = num_eqns_;
    for (int i = 0; i < n; ++i) {
      int offset = bands[i].BandStart();
      T raw_band = bands[i].RawBand();
      V value = bands[i].RawValue();

      while (reduced_matrix[offset] != 0) {
        raw_band ^= reduced_matrix[offset];
        value ^= reduced_values[offset];
        while ((raw_band & 1) == 0) {
          raw_band >>= 1;
          ++offset;
        }
      }
      if (raw_band == 0) {
        continue;
      }
      reduced_matrix[offset] = raw_band;
      reduced_values[offset] = value;
    }
    return true;
  }

  // oc::block type
  template<typename T, typename V>
  inline typename std::enable_if<std::is_same<V, oc::block>::value, bool>::type
  ReduceToRowEchelon(const BandAndValue<T, V>* bands,
                     T* reduced_matrix,
                     V* reduced_values) {
    return ReduceToRowEchelonHelper<T, V>
        (bands, reduced_matrix, reduced_values);
  }

  // all other types
  template<typename T, typename V>
  inline typename std::enable_if<not std::is_same<V, oc::block>::value,
                                 bool>::type
  ReduceToRowEchelon(const BandAndValue<T, V>* bands,
                     T* reduced_matrix,
                     V* reduced_values) {
    return ReduceToRowEchelonHelper<T, V>
        (bands, reduced_matrix, reduced_values, 0);
  }

  template<typename T, typename V>
  bool EncodeHelper(const BandAndValue<T, V>* bands,
              V* out) {
    std::vector<T> reduced_matrix(num_vars_);

    if (!ReduceToRowEchelon<T, V>(bands, reduced_matrix.data(), out)) {
      return false;
    }

    Solve<T, V>(reduced_matrix.data(), out);
    return true;
  }

  template<typename T, typename V>
  bool EncodeHelper(const oc::block* keys, const V* values, V* out) {
    int n = num_eqns_;
    oc::PRNG prng(prng_seed_);
    prng.get(out, num_vars_);
    BandAndValue<T, V>* bands =
        (BandAndValue<T, V>*) malloc(n * sizeof(BandAndValue<T, V>));
    GenBandsAndValues<T, V>(n, keys, values, num_vars_,
                            band_length_, bands);
    //std::sort(bands, bands + n, lessthan());
    boost::sort::spreadsort::integer_sort(bands, bands + n,
                                          rightshift(), lessthan());

    bool encoded = EncodeHelper<T, V>(bands, out);
    //free(bands);
    return encoded;
  }

  template<typename T, typename V>
  inline void DecodeHelper(const Band<T>* bands, const V* okvs, V* out) {
    int n = num_eqns_;
    for (int i = 0; i < n; ++i) {
      int band_start = bands[i].BandStart();
      __m512i res = DoXor(band_start, bands[i].RawBand(), okvs, band_length_);
      __m128i res2 = DoXor512(res);

      _mm_storeu_si128((__m128i * ) & out[bands[i].Index()], res2);
    }
  }

  template<typename T, typename V>
  void DecodeHelper(const oc::block* keys, const V* okvs, V* out) {
    int n = num_eqns_;
    Band<T>* bands = (Band<T>*) malloc(n * sizeof(Band<T>));

    GenBands<T>(n, keys, num_vars_, band_length_,
                bands);

    boost::sort::spreadsort::integer_sort(bands, bands + n,
                                          rightshift(), lessthan());

    DecodeHelper(bands, okvs, out);
    //free(bands);
  }

  inline static constexpr __mmask8
  table_[16] = {
    0,
        3,
        12,
        15,
        48,
        51,
        60,
        63,
        192,
        195,
        204,
        207,
        240,
        243,
        252,
        255
  };

  // Xors values in range [i, i + length) according to the mask defined by
  // raw_band. Assumes V is 128 bit type.
  template<typename T, typename V>
  inline __m512i DoXor(int i, const T& raw_band, const V* reduced_values,
                       int length) {
    __m512i res = _mm512_setzero_si512();
    if (raw_band == 0) {
      return res;
    }

#pragma GCC unroll 8
    for (int j = 0, k = 0; j < length; j += 8, k++) {
      uint8_t raw_band_mask = *((uint8_t * )(&raw_band) + k);

      __mmask8 mask = table_[raw_band_mask & 0b1111];
      res = _mm512_xor_epi64(
          res, _mm512_maskz_loadu_epi64(mask, &reduced_values[i + j]));

      __mmask8 mask2 = table_[(raw_band_mask >> 4) & 0b1111];
      res = _mm512_xor_epi64(res, _mm512_maskz_loadu_epi64(
          mask2, &reduced_values[i + j + 4]));
    }
    return res;
  }

  // Xors the 4 128 bit blocks in the 512 bit register.
  static inline __m128i DoXor512(__m512i n) {
    __m128i a = _mm512_extracti64x2_epi64(n, 0);
    __m128i b = _mm512_extracti64x2_epi64(n, 1);
    __m128i c = _mm512_extracti64x2_epi64(n, 2);
    __m128i d = _mm512_extracti64x2_epi64(n, 3);
    __m128i ab = _mm_xor_si128(a, b);
    __m128i cd = _mm_xor_si128(c, d);
    __m128i res = _mm_xor_si128(ab, cd);
    return res;
  }

  template<typename T, typename V>
  void Solve(const T* reduced_matrix, V* reduced_values) {
    for (int i = num_vars_ - 1; i >= 0; --i) {
      if (reduced_matrix[i] == 0) {
        continue;
      }
      __m512i res = DoXor(i, reduced_matrix[i], reduced_values, band_length_);
      __m128i res2 = DoXor512(res);

      _mm_storeu_si128((__m128i * ) & reduced_values[i], res2);
    }
  }

  int num_eqns_ = 0;
  int num_vars_ = 0;
  int band_length_ = 0;
  oc::block prng_seed_;
};

}

#endif //BANDOKVS_BAND_OKVS_H_
