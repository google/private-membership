#ifndef LINEAR_SYSTEM_H
#define LINEAR_SYSTEM_H

#include <algorithm>
#include <chrono>

#include "band_row_vector.h"
#include "bin_allocator.h"
#include "band_row_vector_hasher.h"
#include "linear_system_row.h"
#include "utils.h"
#include "cryptoTools/Common/block.h"
#include "cryptoTools/Crypto/AES.h"
#include "cryptoTools/Crypto/PRNG.h"

using namespace std;

namespace oc = osuCrypto;

template<typename ValueType>
class LinearSystem {
public:
  LinearSystem() {}

  void Init(int n, int m, int band_width, uint64_t field_size) {
    num_equations_ = n;
    num_vars_ = m;
    band_width_ = band_width;
    field_size_ = field_size;
  }

  void
  SetKeysAndValues(const vector<oc::block> &keys, const vector<int> &values_indices, const vector<ValueType> *values) {
    keys_ = keys;
    values_indices_ = values_indices;
    values_ptr_ = values;
  }

  struct Key {
    int index = 0;
    int offset = 0;
    oc::block key;
  };

  template<typename Hasher>
  bool ReduceToRowEchelonForm(bool reduce_values, const Hasher &hasher) {
    reduced_matrix_ = vector<BandRowVector>(num_vars_);
    if (reduce_values) {
      reduced_values_ = vector<ValueType>(num_vars_);
    }

    std::vector<Key> sorted_keys(keys_.size());
    for (int i = 0; i < keys_.size(); ++i) {
      Band<uint256> band = hasher.Hash(keys_[i]);
      sorted_keys[i] = {i, band.BandStart(), keys_[i]};
    }
    sort(sorted_keys.begin(), sorted_keys.end(), [](const auto &a, const auto &b) {
      return a.offset < b.offset;
    });

    __m512i q_batch = _mm512_set1_epi64(field_size_);
    for (int i = 0; i < sorted_keys.size(); ++i) {
      BandRowVector band = hasher.HashToBandRowVector(sorted_keys[i].key);
      //band.Print();
      ValueType value;
      if (reduce_values) {
        value = (*values_ptr_)[values_indices_[sorted_keys[i].index]];
      }

      int offset = band.GetOffset();
      while (!reduced_matrix_[offset].IsZero()) {
        const BandRowVector &other = reduced_matrix_[offset];
        uint64_t y = MultInverse(band.GetFirstElement(), band.Modulus());
        y = MulMod(y, other.GetFirstElement(),
                   (static_cast<uint128_t>(other.GetFirstElement()) << 64) / band.Modulus(),
                   band.Modulus());
        __m512i y_batch = _mm512_set1_epi64(y);
        uint64_t y_barrett = (static_cast<__uint128_t>(y) << 64) / field_size_;
        __m512i y_barrett_batch = _mm512_set1_epi64(y_barrett);
        band.Reduce(reduced_matrix_[offset], y, y_batch, y_barrett, y_barrett_batch, q_batch);
        if (reduce_values) {
          value.ScalarMultSubtractVecInPlace(y, y_batch, y_barrett, y_barrett_batch, reduced_values_[offset], q_batch);
        }
        if (band.IsZero()) {
          cout << "Failed at index: " << i << endl;
          return false;
        }
        offset = band.GetOffset();
      }
      reduced_matrix_[offset] = std::move(band);
      if (reduce_values) {
        reduced_values_[offset] = std::move(value);
      }
    }
    return true;
  }

  void SolveLinearSystem() {
    for (int i = num_vars_ - 1; i >= 0; i--) {
      if (reduced_matrix_[i].IsZero()) {
        continue;
      }
      int band_length = reduced_matrix_[i].GetLength();
      uint64_t *raw_band = reduced_matrix_[i].RawBand();
      for (int j = 1; j < band_length; ++j) {
        if (raw_band[j] != 0) {
          reduced_values_[i].SubtractVectorScalarMultInPlace(reduced_values_[i + j], raw_band[j]);
        }
      }
      uint64_t inv = MultInverse(raw_band[0], field_size_);
      reduced_values_[i].ScalarMultInPlace(inv);
    }
  }

  vector<ValueType> GetSolution() {
    return reduced_values_;
  }

private:
  int num_vars_;
  int num_equations_;
  int band_width_;
  uint64_t field_size_;

  std::vector<oc::block> keys_;
  std::vector<int> values_indices_;
  const std::vector<ValueType> *values_ptr_;

  vector<BandRowVector> reduced_matrix_;
  vector<ValueType> reduced_values_;
};

template<typename ValueType, typename Hasher>
class BucketedLinearSystems {
public:
  bool
  SetKeysAndValues(const BinAllocator &bin_allocator, const vector<oc::block> &keys, const vector<ValueType> &values,
                   int dimension, int band_width, int band_width_lower, uint64_t field_size,
                   int max_iter = 30, bool print_result = false) {
    assert(keys.size() == values.size());

    total_elapsed_time_ = chrono::duration<int, milli>(0);
    print_result_ = print_result;
    if (print_result) {
      cout << "Preprocessing..." << endl;
    }
    auto time_begin = chrono::steady_clock::now();
    num_bins_ = bin_allocator.NumBins();
    vector<vector<oc::block>> allocated_keys(num_bins_);
    vector<vector<int>> allocated_values_indices(num_bins_);
    allocated_keys[0].reserve(keys.size());
    allocated_values_indices[0].reserve(values.size());
    for (int i = 0; i < keys.size(); ++i) {
      int index = bin_allocator.AllocateFlattened(keys[i]);
      allocated_keys[index].push_back(keys[i]);
      allocated_values_indices[index].push_back(i);
    }
    auto time_end = chrono::steady_clock::now();
    if (print_result) {
      cout << "Done! Elapsed time: " << chrono::duration_cast<chrono::milliseconds>(time_end - time_begin).count()
           << " ms" << endl;
    }

    uint64_t seed = 0;
    int iteration = 0;
    if (print_result) {
      cout << "Transforming keys and values to linear systems..." << endl;
    }
    chrono::duration<int, milli> total_elapsed_time(0);
    for (int iter = 0; iter < max_iter; iter++) {
      if (print_result) {
        cout << "Iteration " << iteration << endl;
      }
      vector<vector<LinearSystem<ValueType>>>
              linear_systems(num_bins_,
                             vector<LinearSystem<ValueType>>(bin_allocator.NumChunks()));

      for (int i = 0; i < linear_systems.size(); ++i) {
        for (int j = 0; j < linear_systems[i].size(); ++j) {
          int eqn_per_chunk = ceil((allocated_keys[i].size() + 0.0) / bin_allocator.NumChunks());
          linear_systems[i][j].Init(eqn_per_chunk, bin_allocator.BinSize(), band_width, field_size);
          vector<oc::block> keys_per_chunk(allocated_keys[i].begin() + j * eqn_per_chunk,
                                           allocated_keys[i].begin() + min(static_cast<int>(allocated_keys[i].size()),
                                                                           (j + 1) * eqn_per_chunk));
          vector<int> values_indices_per_chunk(allocated_values_indices[i].begin() + j * eqn_per_chunk,
                                               allocated_values_indices[i].begin() +
                                               min(static_cast<int>(allocated_values_indices[i].size()),
                                                   (j + 1) * eqn_per_chunk));
          linear_systems[i][j].SetKeysAndValues(keys_per_chunk, values_indices_per_chunk, &values);
        }
      }

      Hasher hasher(bin_allocator, dimension, band_width, field_size, seed);
      bool ok = true;
      if (print_result) {
        cout << "Reducing to row echelon form..." << endl;
      }
      chrono::duration<int, milli> elapsed_time(0);
      int failed_index = 0;
      int failed_chunk = 0;
      for (int i = 0; i < num_bins_ && ok; ++i) {
        for (int j = 0; j < bin_allocator.NumChunks(); ++j) {
          auto time_begin = chrono::steady_clock::now();
          bool reduced = linear_systems[i][j].ReduceToRowEchelonForm(false, hasher);
          auto time_end = chrono::steady_clock::now();
          elapsed_time += chrono::duration_cast<chrono::milliseconds>(time_end - time_begin);
          total_elapsed_time_ += chrono::duration_cast<chrono::milliseconds>(time_end - time_begin);
          if (!reduced) {
            ok = false;
            failed_index = i;
            failed_chunk = j;
            break;
          } else {
            total_elapsed_time_ -= chrono::duration_cast<chrono::milliseconds>(time_end - time_begin);
          }
        }
      }
      total_elapsed_time += elapsed_time;
      if (ok) {
        if (print_result) {
          cout << "Row reduction complete!" << endl;
          cout << "Elapsed time : " << elapsed_time.count() << " ms" << endl;
          cout << endl;
          cout << "Total elapsed time: " << total_elapsed_time.count() << " ms" << endl;
          cout << endl;
        }
        bucketed_systems_ = move(linear_systems);
        /*for (auto& bucket : bucketed_systems_) {
          for (auto& chunk : bucket) {
            chunk.FinalizeValues();
          }
        }*/
        hasher_ = hasher;
        return true;
      }
      if (print_result) {
        cout << "Row reduction failed for bucket at index " << failed_index << " and chunk " << failed_chunk << endl;
        cout << "Elapsed time : " << elapsed_time.count() << " ms" << endl;
        cout << endl;
      }
      ++seed;
      ++iteration;
    }
    return false;
  }

  Hasher &GetHasher() {
    return hasher_;
  }

  vector<vector<vector<ValueType>>> SolveLinearSystems() {
    if (print_result_) {
      cout << "Solving linear systems..." << endl;
    }
    auto time_begin = chrono::steady_clock::now();
    vector<vector<vector<ValueType>>> solutions(num_bins_);
    for (int i = 0; i < num_bins_; ++i) {
      for (int j = 0; j < bucketed_systems_[i].size(); ++j) {
        auto time_begin = chrono::steady_clock::now();
        bool reduced = bucketed_systems_[i][j].ReduceToRowEchelonForm(true, hasher_);
        auto time_end = chrono::steady_clock::now();
        total_elapsed_time_ += chrono::duration_cast<chrono::milliseconds>(time_end - time_begin);
        if (print_result_) {
          cout << "Row reduction done" << endl;
          cout << chrono::duration_cast<chrono::milliseconds>(time_end - time_begin).count() << " ms" << endl;
        }
        assert(reduced);
        time_begin = chrono::steady_clock::now();
        bucketed_systems_[i][j].SolveLinearSystem();
        time_end = chrono::steady_clock::now();
        total_elapsed_time_ += chrono::duration_cast<chrono::milliseconds>(time_end - time_begin);
        if (print_result_) {
          cout << "Back substitution done" << endl;
          cout << chrono::duration_cast<chrono::milliseconds>(time_end - time_begin).count() << " ms" << endl;
        }
        solutions[i].push_back(bucketed_systems_[i][j].GetSolution());
      }
    }
    auto time_end = chrono::steady_clock::now();
    if (print_result_) {
      cout << "Done! Elapsed time: " <<
           chrono::duration_cast<chrono::milliseconds>(time_end - time_begin).count() << " ms" << endl;
    }
    return solutions;
  }

  int TotalElapsedTime() const {
    return total_elapsed_time_.count();
  }

private:
  int num_bins_;
  vector<vector<LinearSystem<ValueType>>> bucketed_systems_;
  Hasher hasher_;
  chrono::duration<int, milli> total_elapsed_time_;
  bool print_result_;
};

#endif //LINEAR_SYSTEM_H
