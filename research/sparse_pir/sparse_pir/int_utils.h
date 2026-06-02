#ifndef SPARSE_PIR_INT_UTILS_H
#define SPARSE_PIR_INT_UTILS_H

#include <immintrin.h>
#include <cstdint>

using namespace std;

template <int InputModFactor = 2>
inline __m512i _mm512_small_mod_epu64(__m512i x, __m512i q,
                                           __m512i* q_times_2 = nullptr,
                                           __m512i* q_times_4 = nullptr) {
  if (InputModFactor == 1) {
    return x;
  }
  if (InputModFactor == 2) {
    return _mm512_min_epu64(x, _mm512_sub_epi64(x, q));
  }
  if (InputModFactor == 4) {
    x = _mm512_min_epu64(x, _mm512_sub_epi64(x, *q_times_2));
    return _mm512_min_epu64(x, _mm512_sub_epi64(x, q));
  }
  if (InputModFactor == 8) {
    x = _mm512_min_epu64(x, _mm512_sub_epi64(x, *q_times_4));
    x = _mm512_min_epu64(x, _mm512_sub_epi64(x, *q_times_2));
    return _mm512_min_epu64(x, _mm512_sub_epi64(x, q));
  }
  return x;  // Return dummy value
}

inline __m512i _mm512_mulhi_epi64(__m512i x, __m512i y) {
  __m512i lo_mask = _mm512_set1_epi64(0x00000000ffffffff);
  __m512i x_hi = _mm512_shuffle_epi32(x, (_MM_PERM_ENUM)0xB1);
  __m512i y_hi = _mm512_shuffle_epi32(y, (_MM_PERM_ENUM)0xB1);
  __m512i z_lo_lo = _mm512_mul_epu32(x, y);        // x_lo * y_lo
  __m512i z_lo_hi = _mm512_mul_epu32(x, y_hi);     // x_lo * y_hi
  __m512i z_hi_lo = _mm512_mul_epu32(x_hi, y);     // x_hi * y_lo
  __m512i z_hi_hi = _mm512_mul_epu32(x_hi, y_hi);  // x_hi * y_hi

  __m512i z_lo_lo_shift = _mm512_srli_epi64(z_lo_lo, 32);

  __m512i sum_tmp = _mm512_add_epi64(z_lo_hi, z_lo_lo_shift);
  __m512i sum_lo = _mm512_and_si512(sum_tmp, lo_mask);
  __m512i sum_mid = _mm512_srli_epi64(sum_tmp, 32);

  __m512i sum_mid2 = _mm512_add_epi64(z_hi_lo, sum_lo);
  __m512i sum_mid2_hi = _mm512_srli_epi64(sum_mid2, 32);
  __m512i sum_hi = _mm512_add_epi64(z_hi_hi, sum_mid);
  return _mm512_add_epi64(sum_hi, sum_mid2_hi);
}

inline uint64_t MulMod(uint64_t x, uint64_t y, uint64_t y_barrett, uint64_t q) {
  uint64_t xy = x * y;
  uint64_t r = (static_cast<__uint128_t>(x) * y_barrett) >> 64;
  uint64_t rq = r * q;
  r = xy - rq;
  r = r < q ? r : r - q;
  return r;
}

inline __m512i MulMod(__m512i x, __m512i y, __m512i y_barrett, __m512i q) {
  __m512i xy = _mm512_mullo_epi64(x, y);
  __m512i r = _mm512_mulhi_epi64(x, y_barrett);
  __m512i rq = _mm512_mullo_epi64(r, q);
  r = _mm512_sub_epi64(xy, rq);
  r = _mm512_small_mod_epu64(r, q);
  return r;
}

// x * y - z
inline uint64_t MulSubMod(uint64_t x, uint64_t y, uint64_t z, uint64_t y_barrett, uint64_t q) {
  uint64_t xy = x * y;
  uint64_t r = (static_cast<__uint128_t>(x) * y_barrett) >> 64;
  uint64_t rq = r * q;
  r = xy - rq;
  r = r < q ? r : r - q;
  r -= z;
  r = r < q ? r : r + q;
  return r;
}

// x * y - z
inline __m512i MulSubMod(__m512i x, __m512i y, __m512i z, __m512i y_barrett, __m512i q) {
  __m512i xy = _mm512_mullo_epi64(x, y);
  __m512i r = _mm512_mulhi_epi64(x, y_barrett);
  __m512i rq = _mm512_mullo_epi64(r, q);
  r = _mm512_sub_epi64(xy, rq);
  r = _mm512_small_mod_epu64(r, q);
  r = _mm512_sub_epi64(r, z);
  r = _mm512_min_epu64(r, _mm512_add_epi64(r, q));
  return r;
}

// x - y * z
inline uint64_t SubMulMod(uint64_t x, uint64_t y, uint64_t z, uint64_t z_barrett, uint64_t q) {
  uint64_t yz = y * z;
  uint64_t r = (static_cast<__uint128_t>(y) * z_barrett) >> 64;
  uint64_t rq = r * q;
  r = yz - rq;
  r = r < q ? r : r - q;
  r = x - r;
  r = r < q ? r : r + q;
  return r;
}

// x - y * z
inline __m512i SubMulMod(__m512i x, __m512i y, __m512i z, __m512i z_barrett, __m512i q) {
  __m512i yz = _mm512_mullo_epi64(y, z);
  __m512i r = _mm512_mulhi_epi64(y, z_barrett);
  __m512i rq = _mm512_mullo_epi64(r, q);
  r = _mm512_sub_epi64(yz, rq);
  r = _mm512_small_mod_epu64(r, q);
  r = _mm512_sub_epi64(x, r);
  r = _mm512_min_epu64(r, _mm512_add_epi64(r, q));
  return r;
}

inline void ScalarMulMod(uint64_t* v1, uint64_t y, const __m512i& y_batch, uint64_t y_barrett,
                                  const __m512i& y_barrett_batch,
                                  uint64_t q,
                                  const __m512i& q_batch, int length) {
  static constexpr int batch_size = sizeof(__m512i) / sizeof(uint64_t);
  int i = 0;
  for (; i + batch_size < length; i += batch_size) {
    __m512i x_batch = _mm512_loadu_epi64(v1 + i);
    __m512i res = MulMod(x_batch, y_batch, y_barrett_batch, q_batch);
    _mm512_storeu_epi64(v1 + i, res);
  }
  for (; i < length; ++i) {
    v1[i] = MulMod(v1[i], y, y_barrett, q);
  }
}

inline void ScalarMulSubVecMod(uint64_t* v1, const uint64_t* v2, uint64_t y, const __m512i& y_batch, uint64_t y_barrett,
                                  const __m512i& y_barrett_batch, uint64_t q, const __m512i& q_batch, int length) {
  if (length == 0) {
    return;
  }
  static constexpr int batch_size = sizeof(__m512i) / sizeof(uint64_t);
  int i = 0;
  for (; i + batch_size < length; i += batch_size) {
    __m512i x_batch = _mm512_loadu_epi64(v1 + i);
    __m512i z_batch = _mm512_loadu_epi64(v2 + i);
    __m512i res = MulSubMod(x_batch, y_batch, z_batch, y_barrett_batch, q_batch);
    _mm512_storeu_epi64(v1 + i, res);
  }
  for (; i < length; ++i) {
    v1[i] = MulSubMod(v1[i], y, v2[i], y_barrett, q);
  }
}

inline void SubVecScalarMulMod(uint64_t* v1, const uint64_t* v2, uint64_t z, const __m512i& z_batch, uint64_t z_barrett,
                                  const __m512i& z_barrett_batch, uint64_t q, const __m512i& q_batch, int length) {
  if (length == 0) {
    return;
  }
  static constexpr int batch_size = sizeof(__m512i) / sizeof(uint64_t);
  int i = 0;
  for (; i + batch_size < length; i += batch_size) {
    __m512i x_batch = _mm512_loadu_epi64(v1 + i);
    __m512i y_batch = _mm512_loadu_epi64(v2 + i);
    __m512i res = SubMulMod(x_batch, y_batch, z_batch, z_barrett_batch, q_batch);
    _mm512_storeu_epi64(v1 + i, res);
  }
  for (; i < length; ++i) {
    v1[i] = SubMulMod(v1[i], v2[i], z, z_barrett, q);
  }
}

uint64_t BinPow(uint64_t a, uint64_t b, uint64_t mod) {
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

uint64_t MultInverse(uint64_t a, uint64_t mod) {
  return BinPow(a, mod - 2, mod);
}

#endif //SPARSE_PIR_INT_UTILS_H
