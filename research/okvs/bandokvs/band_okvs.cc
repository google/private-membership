#ifndef BANDOKVS_BAND_OKVS_H_
#define BANDOKVS_BAND_OKVS_H_

#include <vector>

#include "band.h"
#include "band_okvs.h"
#include "uint.h"
#include "utils.h"

#include "cryptoTools/Crypto/PRNG.h"

namespace band_okvs {

std::vector<__uint128_t> GenRandomValues(int n, uint64_t seed1,
                                         uint64_t seed2) {
  std::vector<oc::block> blocks(n);
  oc::PRNG prng(oc::toBlock(seed2, seed1));
  prng.get<oc::block>(blocks);

  std::vector<__uint128_t> res(n);
  for (int i = 0; i < n; ++i) {
    res[i] = blocks[i].get<__uint128_t>(0);
  }
  return res;
}

std::vector<oc::block> GenRandomValuesBlocks(int n,
                                             uint64_t seed1,
                                             uint64_t seed2) {
  std::vector<oc::block> blocks(n);
  oc::PRNG prng(oc::toBlock(seed2, seed1));
  prng.get<oc::block>(blocks);
  return blocks;
}

void GenRandomValuesBlocks(oc::span<oc::block> out, uint64_t seed1,
                           uint64_t seed2) {
  oc::PRNG prng(oc::toBlock(seed2, seed1));
  prng.get<oc::block>(out);
}

}

#endif //BANDOKVS_BAND_OKVS_H_