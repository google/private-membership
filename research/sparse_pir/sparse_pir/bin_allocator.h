#ifndef BIN_ALLOCATOR_H
#define BIN_ALLOCATOR_H

#include <map>
#include <string>
#include <vector>

#include "cryptoTools/Common/block.h"
#include "cryptoTools/Crypto/AES.h"

using namespace std;

namespace oc = osuCrypto;

class Bin {
public:
  explicit Bin(const vector<oc::block> &keys) : keys_(keys) {}

  const vector<oc::block> &Keys() const {
    return keys_;
  }

private:
  vector<oc::block> keys_;
};

class BinAllocator {
public:
  BinAllocator(int num_bins, int bin_size, int num_chunks = 1) : hash_seed_(0),
                                                                 num_bins_(num_bins), bin_size_(bin_size),
                                                                 num_chunks_(num_chunks), aes_(oc::mAesFixedKey) {}

  bool SetKeys(const vector<oc::block> &keys, int bin_buffer = 0, int max_tries = 1000) {
    auto time_begin = chrono::steady_clock::now();
    hash_seed_ = 0;
    for (int j = 0; j < max_tries; ++j) {
      vector<int> bins(num_bins_);
      bool exceeded = false;
      for (const auto &key: keys) {
        int allocated_bin = Allocate(key);
        if (bins[allocated_bin] > num_chunks_ * (bin_size_ - bin_buffer)) {
          exceeded = true;
          break;
        } else {
          ++bins[allocated_bin];
        }
      }
      if (!exceeded) {
        num_flattened_bins_ = num_bins_;
        keys_to_bin_indices_.clear();
        for (const auto &key: keys) {
          keys_to_bin_indices_[key] = Allocate(key);
        }
        auto time_end = chrono::steady_clock::now();
        cout << "Found allocation in " << chrono::duration_cast<chrono::milliseconds>(time_end - time_begin).count()
             << " ms" << endl;
        return true;
      }
      ++hash_seed_;
    }
    return false;
  }

  int BinSize() const {
    return bin_size_;
  }

  int NumBins() const {
    return num_flattened_bins_;
  }

  int TotalCapacity() const {
    return bin_size_ * num_flattened_bins_;
  }

  int NumChunks() const {
    return num_chunks_;
  }

  int AllocateFlattened(const oc::block &key) const {
    return keys_to_bin_indices_.at(key);
  }

private:
  int Allocate(const oc::block &key) const {
    oc::block level_hash = aes_.hashBlock(oc::block(0));
    uint64_t hash = (aes_.hashBlock(key) ^ level_hash).get<uint64_t>(0);
    return hash % num_bins_;
  }

  uint64_t hash_seed_;
  int num_bins_;
  int bin_size_;
  int num_chunks_;
  oc::AES aes_;

  int num_flattened_bins_;
  map<oc::block, int> keys_to_bin_indices_;
};

#endif //BIN_ALLOCATOR_H