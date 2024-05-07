#ifndef BANDOKVS_UINT_H_
#define BANDOKVS_UINT_H_

namespace band_okvs {

template<int Size>
class uint;

constexpr inline uint<1> MakeUint(__uint128_t v0);
constexpr inline uint<2> MakeUint(__uint128_t v0, __uint128_t v1);
constexpr inline uint<3> MakeUint(__uint128_t v0, __uint128_t v1,
                                  __uint128_t v2);
constexpr inline uint<4> MakeUint(__uint128_t v0, __uint128_t v1,
                                  __uint128_t v2, __uint128_t v3);
constexpr inline uint<5> MakeUint(__uint128_t v0,
                                  __uint128_t v1,
                                  __uint128_t v2,
                                  __uint128_t v3,
                                  __uint128_t v4);
constexpr inline uint<6> MakeUint(__uint128_t v0,
                                  __uint128_t v1,
                                  __uint128_t v2,
                                  __uint128_t v3,
                                  __uint128_t v4,
                                  __uint128_t v5);

template<int Size>
constexpr inline uint<Size> operator|(uint<Size> lhs, uint<Size> rhs);

template<int Size>
constexpr inline uint<Size> operator&(uint<Size> lhs, uint<Size> rhs);

template<int Size>
constexpr inline __uint128_t operator&(uint<Size> lhs, int rhs);

template<int Size>
constexpr inline uint<Size> operator^(uint<Size> lhs, uint<Size> rhs);

template<int Size>
constexpr inline uint<Size> operator>>(uint<Size> lhs, int amount);

template<int Size>
constexpr inline bool operator==(uint<Size> lhs, uint<Size> rhs);

template<int Size>
constexpr inline bool operator!=(uint<Size> lhs, uint<Size> rhs);

template<int Size>
class uint {
 public:
  uint() = default;

  void set(oc::block* blocks, int len, __uint128_t mask) {
    std::memcpy(values_, blocks, len * sizeof(oc::block));
    values_[0] |= 1;
    values_[Size - 1] &= mask;
  }

  constexpr uint(__uint128_t v0) {
    values_[0] = v0;
  }

  constexpr uint(__uint128_t v0, __uint128_t v1) {
    values_[0] = v0;
    values_[1] = v1;
  }

  constexpr uint(__uint128_t v0, __uint128_t v1, __uint128_t v2) {
    values_[0] = v0;
    values_[1] = v1;
    values_[2] = v2;
  }

  constexpr uint(__uint128_t v0, __uint128_t v1, __uint128_t v2,
                 __uint128_t v3) {
    values_[0] = v0;
    values_[1] = v1;
    values_[2] = v2;
    values_[3] = v3;
  }

  constexpr uint(__uint128_t v0, __uint128_t v1, __uint128_t v2,
                 __uint128_t v3, __uint128_t v4) {
    values_[0] = v0;
    values_[1] = v1;
    values_[2] = v2;
    values_[3] = v3;
    values_[4] = v4;
  }

  constexpr uint(__uint128_t v0, __uint128_t v1, __uint128_t v2,
                 __uint128_t v3, __uint128_t v4, __uint128_t v5) {
    values_[0] = v0;
    values_[1] = v1;
    values_[2] = v2;
    values_[3] = v3;
    values_[4] = v4;
    values_[5] = v5;
  }

  inline uint& operator=(int v) { return *this = uint(v); }

  constexpr explicit operator int() const { return static_cast<int>(values_[0]); }

  inline uint& operator&=(uint<Size> other) {
    *this = *this & other;
    return *this;
  }

  inline uint& operator|=(uint<Size> other) {
    *this = *this | other;
    return *this;
  }

  inline uint& operator^=(uint<Size> other) {
    *this = *this ^ other;
    return *this;
  }

  inline uint& operator>>=(int amount) {
    *this = *this >> amount;
    return *this;
  }

  inline void ToggleBit(int bit_pos) {
    int int_index = bit_pos / 128;
    int bit_index = bit_pos % 128;
    values_[int_index] ^= (static_cast<__uint128_t>(1) << bit_index);
  }

  inline int GetBit(int bit_pos) const {
    int int_index = bit_pos / 128;
    int bit_index = bit_pos % 128;
    return (values_[int_index] & (static_cast<__uint128_t>(1) << bit_index))
        >> bit_index;
  }

  friend constexpr __uint128_t Uint(uint<Size> v,
                                    int index) { return v.values_[index]; }

 private:
  __uint128_t values_[Size] = {0};
};

constexpr inline uint<1> MakeUint(__uint128_t v0) {
  return uint<1>(v0);
}

constexpr inline uint<2> MakeUint(__uint128_t v0, __uint128_t v1) {
  return uint<2>(v0, v1);
}

constexpr inline uint<3> MakeUint(__uint128_t v0, __uint128_t v1,
                                  __uint128_t v2) {
  return uint<3>(v0, v1, v2);
}

constexpr inline uint<4> MakeUint(__uint128_t v0, __uint128_t v1,
                                  __uint128_t v2, __uint128_t v3) {
  return uint<4>(v0, v1, v2, v3);
}

constexpr inline uint<5> MakeUint(__uint128_t v0, __uint128_t v1,
                                  __uint128_t v2, __uint128_t v3,
                                  __uint128_t v4) {
  return uint<5>(v0, v1, v2, v3, v4);
}

constexpr inline uint<6> MakeUint(__uint128_t v0, __uint128_t v1,
                                  __uint128_t v2, __uint128_t v3,
                                  __uint128_t v4, __uint128_t v5) {
  return uint<6>(v0, v1, v2, v3, v4, v5);
}

template<>
constexpr inline uint<1> operator|(uint<1> lhs, uint<1> rhs) {
  return MakeUint(Uint(lhs, 0) | Uint(rhs, 0));
}

template<>
constexpr inline uint<2> operator|(uint<2> lhs, uint<2> rhs) {
  return MakeUint(Uint(lhs, 0) | Uint(rhs, 0),
                  Uint(lhs, 1) | Uint(rhs, 1));
}

template<>
constexpr inline uint<3> operator|(uint<3> lhs, uint<3> rhs) {
  return MakeUint(Uint(lhs, 0) | Uint(rhs, 0),
                  Uint(lhs, 1) | Uint(rhs, 1),
                  Uint(lhs, 2) | Uint(rhs, 2));
}

template<>
constexpr inline uint<4> operator|(uint<4> lhs, uint<4> rhs) {
  return MakeUint(Uint(lhs, 0) | Uint(rhs, 0),
                  Uint(lhs, 1) | Uint(rhs, 1),
                  Uint(lhs, 2) | Uint(rhs, 2),
                  Uint(lhs, 3) | Uint(rhs, 3));
}

template<>
constexpr inline uint<5> operator|(uint<5> lhs, uint<5> rhs) {
  return MakeUint(Uint(lhs, 0) | Uint(rhs, 0),
                  Uint(lhs, 1) | Uint(rhs, 1),
                  Uint(lhs, 2) | Uint(rhs, 2),
                  Uint(lhs, 3) | Uint(rhs, 3),
                  Uint(lhs, 4) | Uint(rhs, 4));
}

template<>
constexpr inline uint<6> operator|(uint<6> lhs, uint<6> rhs) {
  return MakeUint(Uint(lhs, 0) | Uint(rhs, 0),
                  Uint(lhs, 1) | Uint(rhs, 1),
                  Uint(lhs, 2) | Uint(rhs, 2),
                  Uint(lhs, 3) | Uint(rhs, 3),
                  Uint(lhs, 4) | Uint(rhs, 4),
                  Uint(lhs, 5) | Uint(rhs, 5));
}

template<>
constexpr inline uint<1> operator&(uint<1> lhs, uint<1> rhs) {
  return MakeUint(Uint(lhs, 0) & Uint(rhs, 0));
}

template<>
constexpr inline uint<2> operator&(uint<2> lhs, uint<2> rhs) {
  return MakeUint(Uint(lhs, 0) & Uint(rhs, 0),
                  Uint(lhs, 1) & Uint(rhs, 1));
}

template<>
constexpr inline uint<3> operator&(uint<3> lhs, uint<3> rhs) {
  return MakeUint(Uint(lhs, 0) & Uint(rhs, 0),
                  Uint(lhs, 1) & Uint(rhs, 1),
                  Uint(lhs, 2) & Uint(rhs, 2));
}

template<>
constexpr inline uint<4> operator&(uint<4> lhs, uint<4> rhs) {
  return MakeUint(Uint(lhs, 0) & Uint(rhs, 0),
                  Uint(lhs, 1) & Uint(rhs, 1),
                  Uint(lhs, 2) & Uint(rhs, 2),
                  Uint(lhs, 3) & Uint(rhs, 3));
}

template<>
constexpr inline uint<5> operator&(uint<5> lhs, uint<5> rhs) {
  return MakeUint(Uint(lhs, 0) & Uint(rhs, 0),
                  Uint(lhs, 1) & Uint(rhs, 1),
                  Uint(lhs, 2) & Uint(rhs, 2),
                  Uint(lhs, 3) & Uint(rhs, 3),
                  Uint(lhs, 4) & Uint(rhs, 4));
}

template<>
constexpr inline uint<6> operator&(uint<6> lhs, uint<6> rhs) {
  return MakeUint(Uint(lhs, 0) & Uint(rhs, 0),
                  Uint(lhs, 1) & Uint(rhs, 1),
                  Uint(lhs, 2) & Uint(rhs, 2),
                  Uint(lhs, 3) & Uint(rhs, 3),
                  Uint(lhs, 4) & Uint(rhs, 4),
                  Uint(lhs, 5) & Uint(rhs, 5));
}

template<int Size>
constexpr inline __uint128_t operator&(uint<Size> lhs, int rhs) {
  return Uint(lhs, 0) & rhs;
}

template<int Size>
constexpr inline __uint128_t operator&(uint<Size> lhs, uint64_t rhs) {
  return Uint(lhs, 0) & rhs;
}

template<int Size>
constexpr inline __uint128_t operator&(uint<Size> lhs, __uint128_t rhs) {
  return Uint(lhs, 0) & rhs;
}

template<>
constexpr inline uint<1> operator^(uint<1> lhs, uint<1> rhs) {
  return MakeUint(Uint(lhs, 0) ^ Uint(rhs, 0));
}

template<>
constexpr inline uint<2> operator^(uint<2> lhs, uint<2> rhs) {
  return MakeUint(Uint(lhs, 0) ^ Uint(rhs, 0),
                  Uint(lhs, 1) ^ Uint(rhs, 1));
}

template<>
constexpr inline uint<3> operator^(uint<3> lhs, uint<3> rhs) {
  return MakeUint(Uint(lhs, 0) ^ Uint(rhs, 0),
                  Uint(lhs, 1) ^ Uint(rhs, 1),
                  Uint(lhs, 2) ^ Uint(rhs, 2));
}

template<>
constexpr inline uint<4> operator^(uint<4> lhs, uint<4> rhs) {
  return MakeUint(Uint(lhs, 0) ^ Uint(rhs, 0),
                  Uint(lhs, 1) ^ Uint(rhs, 1),
                  Uint(lhs, 2) ^ Uint(rhs, 2),
                  Uint(lhs, 3) ^ Uint(rhs, 3));
}

template<>
constexpr inline uint<5> operator^(uint<5> lhs, uint<5> rhs) {
  return MakeUint(Uint(lhs, 0) ^ Uint(rhs, 0),
                  Uint(lhs, 1) ^ Uint(rhs, 1),
                  Uint(lhs, 2) ^ Uint(rhs, 2),
                  Uint(lhs, 3) ^ Uint(rhs, 3),
                  Uint(lhs, 4) ^ Uint(rhs, 4));
}

template<>
constexpr inline uint<6> operator^(uint<6> lhs, uint<6> rhs) {
  return MakeUint(Uint(lhs, 0) ^ Uint(rhs, 0),
                  Uint(lhs, 1) ^ Uint(rhs, 1),
                  Uint(lhs, 2) ^ Uint(rhs, 2),
                  Uint(lhs, 3) ^ Uint(rhs, 3),
                  Uint(lhs, 4) ^ Uint(rhs, 4),
                  Uint(lhs, 5) ^ Uint(rhs, 5));
}

template<int Size>
constexpr inline uint<Size> operator>>(uint<Size> lhs, int amount);

template<>
constexpr inline uint<1> operator>>(uint<1> lhs, int amount) {
  // Only for amount < 128
  return MakeUint(Uint(lhs, 0) >> amount);
}

template<>
constexpr inline uint<2> operator>>(uint<2> lhs, int amount) {
  // Only for amount < 128
  return MakeUint(
      Uint(lhs, 0) >> amount | (Uint(lhs, 1) << (128 - amount)),
      Uint(lhs, 1) >> amount);
}

template<>
constexpr inline uint<3> operator>>(uint<3> lhs, int amount) {
  // Only for amount < 128
  return MakeUint(
      Uint(lhs, 0) >> amount | (Uint(lhs, 1) << (128 - amount)),
      Uint(lhs, 1) >> amount | (Uint(lhs, 2) << (128 - amount)),
      Uint(lhs, 2) >> amount);
}

template<>
constexpr inline uint<4> operator>>(uint<4> lhs, int amount) {
  // Only for amount < 128
  return MakeUint(
      Uint(lhs, 0) >> amount | (Uint(lhs, 1) << (128 - amount)),
      Uint(lhs, 1) >> amount | (Uint(lhs, 2) << (128 - amount)),
      Uint(lhs, 2) >> amount | (Uint(lhs, 3) << (128 - amount)),
      Uint(lhs, 3) >> amount);
}

template<>
constexpr inline uint<5> operator>>(uint<5> lhs, int amount) {
  // Only for amount < 128
  return MakeUint(
      Uint(lhs, 0) >> amount | (Uint(lhs, 1) << (128 - amount)),
      Uint(lhs, 1) >> amount | (Uint(lhs, 2) << (128 - amount)),
      Uint(lhs, 2) >> amount | (Uint(lhs, 3) << (128 - amount)),
      Uint(lhs, 3) >> amount | (Uint(lhs, 4) << (128 - amount)),
      Uint(lhs, 4) >> amount);
}

template<>
constexpr inline uint<6> operator>>(uint<6> lhs, int amount) {
  // Only for amount < 128
  return MakeUint(
      Uint(lhs, 0) >> amount | (Uint(lhs, 1) << (128 - amount)),
      Uint(lhs, 1) >> amount | (Uint(lhs, 2) << (128 - amount)),
      Uint(lhs, 2) >> amount | (Uint(lhs, 3) << (128 - amount)),
      Uint(lhs, 3) >> amount | (Uint(lhs, 4) << (128 - amount)),
      Uint(lhs, 4) >> amount | (Uint(lhs, 5) << (128 - amount)),
      Uint(lhs, 5) >> amount);
}

template<int Size>
constexpr inline bool operator==(uint<Size> lhs, uint<Size> rhs) {
  bool equals = true;
  for (int i = 0; i < Size; i++) {
    equals = equals && (Uint(lhs, i) == Uint(rhs, i));
  }
  return equals;
}

template<int Size>
constexpr inline bool operator!=(uint<Size> lhs, uint<Size> rhs) {
  bool not_equals = false;
  for (int i = 0; i < Size; i++) {
    not_equals = not_equals || (Uint(lhs, i) != Uint(rhs, i));
  }
  return not_equals;
}

template<int Size>
constexpr inline bool operator==(uint<Size> lhs, int rhs) {
  bool equals = Uint(lhs, 0) == rhs;
  for (int i = 1; i < Size; i++) {
    equals = equals && (Uint(lhs, i) == 0);
  }
  return equals;
}

template<int Size>
constexpr inline bool operator!=(uint<Size> lhs, int rhs) {
  bool not_equals = Uint(lhs, 0) != rhs;
  for (int i = 1; i < Size; i++) {
    not_equals = not_equals || (Uint(lhs, i) != 0);
  }
  return not_equals;
}

}

#endif //BANDOKVS_UINT_H_
