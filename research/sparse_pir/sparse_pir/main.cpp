#include <bitset>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "seal/seal.h"

#include "bin_allocator.h"
#include "band_row_vector_hasher.h"
#include "int_utils.h"
#include "linear_system.h"
#include "linear_system_row.h"
#include "plaintext_value.h"
#include "prime_field.h"
#include "utils.h"

#include "cryptoTools/Common/block.h"

using namespace std;

namespace oc = osuCrypto;

/*vector<string> GenerateKeys(int n) {
  vector<string> res;
  res.reserve(n);
  for (int i = 0; i < n; ++i) {
    res.push_back("key-" + to_string(i));
  }
  return res;
}

vector<TestValue> GenerateValues(int n, uint64_t field_size) {
  vector<TestValue> res;
  res.reserve(n);
  for (int i = 0; i < n; ++i) {
    res.push_back(TestValue(PrimeField(i % field_size, field_size)));
  }
  return res;
}
*/

/*std::vector<uint64_t> Encode(const std::vector<unsigned char>& str, size_t plain_modulus_bit_count) {
  std::vector<uint64_t> res;
  std::string bit_str;
  int plain_data_bits = plain_modulus_bit_count - 1;
  int n = str.size();
  int remain = ((n / 2) * 8) % plain_data_bits;
  for (int iter = 0; iter < 2; iter++) {
    int start_byte = iter * (n / 2);
    for (int i = 0; i < n / 2; i++) {
      bit_str += std::bitset<8>(str[start_byte + i]).to_string();
    }
    if (remain != 0) {
      for (int i = 0; i < (plain_data_bits - remain); i++)
        bit_str += "1";
    }
  }
  for (int i = 0; i < bit_str.length(); i += plain_data_bits)
    res.push_back((uint64_t)std::stoi(bit_str.substr(i, plain_data_bits), nullptr, 2));

  return res;
}

std::vector<seal::Plaintext> EncodedValues(const std::vector<string>& values, int N) {
  std::vector<std::vector<uint64_t> > extended_db(values.size(), std::vector<uint64_t>(N, 1));
  int row_size = N / 2;

  for (int i = 0; i < num_obj; i++) {
    std::vector<uint64_t> temp = encode(db[i]);

    int row = (i / row_size);
    int col = (i % row_size);
    for (int j = 0; j < num_columns_per_obj / 2; j++)
    {
      extended_db[row][col] = temp[j];
      extended_db[row][col + row_size] = temp[j + (num_columns_per_obj / 2)];
      row += num_query_ciphertext;
    }

  }

  std::vector<seal::Plaintext> encoded_values(values.size());
}
*/

/*vector<PlaintextValue> GeneratePlaintextValues(int n, uint64_t field_size, seal::BatchEncoder& batch_encoder) {
  vector<PlaintextValue> res;
  res.reserve(n);
  for (int i = 0; i < n; ++i) {
    string value = "value-" + to_string(i);
    batch_encoder.encode()
    res.push_back(PlaintextValue(PrimeField(i % field_size, field_size)));
  }
  return res;
}*/

int Test1(int n, double epsilon, int dimension, int band_width, int band_width_lower, int max_iter,
          const vector<PlaintextValue2>& values, bool test_queries, bool print_result) {

  int num_equations = n;
  int num_chunks = 1;

  //uint64_t field_size = 257;
  //uint64_t field_size = 8388617;
  uint64_t field_size = 1152921504606830593ULL;
  uint64_t t = field_size;
  uint64_t logt = 60;
  uint64_t element_size = 256;

  int bin_size = ceil((1 + epsilon) * n);
  while (bin_size % dimension != 0) {
    ++bin_size;
  }
  if (print_result) {
    cout << "New epsilon: " << (bin_size + 0.0) / n - 1 << endl;
  }
  BinAllocator2 allocator(1, 2, bin_size , num_chunks);
  vector<oc::block> keys = GenerateKeysBlocks(n);
  if (!allocator.SetKeysWithoutOverflow(keys, 0, 1000)) {
    cout << "Failed to find appropriate allocations!" << endl;
    return -1;
  }

  //int max_overflow_bins = ceil(initial_num_bins * 0.01);
  /*int max_overflow_bins = 14;
  cout << "Overflow bins: " << max_overflow_bins << endl;
  cout << "New epsilon: " << ((dimension + 1.0) * (initial_num_bins + max_overflow_bins) / n - 1.0) << endl;
  if (!allocator.SetKeys(keys, 39, max_overflow_bins, 1000)) {
    cout << "Failed to find appropriate allocations!" << endl;
    return 0;
  }*/

  if (print_result) {
    cout << "Num bins: " << allocator.NumBins() << endl;
    cout << "Generating values..." << endl;
  }
  /*vector<PlaintextValue> values = ToPlaintextValues(t, logt, raw_values);

  BucketedLinearSystems<PlaintextValue> linear_systems;
  linear_systems.SetKeysAndValues(allocator, keys, values, band_width, field_size);

  vector<vector<vector<PlaintextValue>>> res = linear_systems.SolveLinearSystems();*/

  if (print_result) {
    cout << "Generating raw values done" << endl;
  }

  if (print_result) {
    cout << "Done!" << endl;
    cout << "Num coeffs: " << values[0].GetValue().size() << endl;
  }

  BucketedLinearSystems<PlaintextValue2, BandHasher> linear_systems;
  if (!linear_systems.SetKeysAndValues(allocator, keys, values, dimension, band_width, band_width_lower, field_size,
                                       max_iter, print_result)) {
    return -1;
  }

  vector<vector<vector<PlaintextValue2>>> res = linear_systems.SolveLinearSystems();

  if (!test_queries) {
    return linear_systems.TotalElapsedTime();
  }

  cout << "Testing queries..." << endl;
  for (int i = 0; i < num_equations; ++i) {
    int index = allocator.AllocateFlattened(keys[i]);
    auto row = linear_systems.GetHasher().HashToBandRowVector(keys[i]);
    //seal::Plaintext& raw_band = row.RawBand();
    int band_length = row.GetLength();
    //cout << band_length << endl;
    uint64_t* raw_band = row.RawBand();
    //band_row_vector.Print();
    bool is_equal = false;
    for (int k = 0; k < num_chunks; ++k) {
      PlaintextValue2 v(field_size);
      for (int j = 0; j < band_length; ++j) {
        if (raw_band[j] != 0) {
          v.AddInPlace(res[index][k][row.GetOffset() + j]);
        }
      }

      /*string raw_value = CoeffsToValue(t, logt, v.GetValue(), element_size);
      if (raw_value == raw_values[i]) {
        is_equal = true;
        break;
      }*/
      if (v == values[i]) {
        //values[i].Print();
        //v.Print();
        is_equal = true;
        break;
      } else {
        //values[i].Print();
        //v.Print();
        break;
      }
      /*for (int j =  0; j < 288; j++) {
        if (raw_value[j] != raw_values[i][j]) {
          cout << j << " " << (int)raw_value[j] << " " << (int)raw_values[i][j] << endl;
        }
      }*/
      //cout << raw_value << endl;
    }

    if (!is_equal) {
      cout << "Different, unexpected result! " << i << endl;
      exit(0);
    }
  }
  cout << "All queries correct!" << endl;
  return linear_systems.TotalElapsedTime();
}

/*void Test2() {
  double eps = 0.02;
  int n = (1 << 18);

  int num_equations = n;
  int num_chunks = 1;
  int dimension = 128;
  int band_width = 128;
  //uint64_t field_size = 257;
  //uint64_t field_size = 8388617;
  uint64_t field_size = 1152921504606830593ULL;
  uint64_t t = field_size;
  uint64_t logt = 60;
  uint64_t element_size = 300;
  cout << field_size << endl;

  int bin_size = ceil((1 + eps) * n);
  while (bin_size % 128 != 0) {
    ++bin_size;
  }
  cout << "New epsilon: " << (bin_size + 0.0) / n - 1 << endl;
  BinAllocator2 allocator(1, 2, bin_size , num_chunks);
  vector<oc::block> keys = GenerateKeysBlocks(n);
  if (!allocator.SetKeysWithoutOverflow(keys, 0, 1000)) {
    cout << "Failed to find appropriate allocations!" << endl;
    return;
  }

  //int max_overflow_bins = ceil(initial_num_bins * 0.01);

  cout << "Num bins: " << allocator.NumBins() << endl;

  cout << "Generating values..." << endl;
  vector<string> raw_values = GenerateRawValues(n, element_size);
  vector<BitVector> values = RawValuesToBitVectors(raw_values);

  BucketedBandLinearSystems<BitVector, BandHasher> linear_systems;
  linear_systems.SetKeysAndValues(allocator, keys, values, band_width);

  auto solution = linear_systems.SolveLinearSystems();
}*/

void Test3() {
  uint64_t y = 12345678912313123LL;
  uint64_t q = 12346763463631414LL;
  uint64_t y_barrett = (static_cast<__uint128_t>(y) << 64) / q;
  for (int i = 0; i < 1000; ++i) {
    cout << static_cast<uint64_t>(static_cast<__uint128_t>(y * i) % q) << " " << MulMod(i, y, y_barrett, q) << endl;
  }
}

void Experiment() {
  vector<double> epsilons = {0.04};
  vector<int> band_widths = {20};

  int n = 1 << 10;
  int num_runs = 5;
  int dimension = 32;
  uint64_t field_size = 1152921504606830593ULL;
  uint64_t t = field_size;
  uint64_t logt = 60;
  uint64_t element_size = 256;
  vector<string> raw_values = GenerateRawValues(n, element_size);
  vector<PlaintextValue2> values = ToPlaintextValues2(t, logt, raw_values);
  for (double epsilon : epsilons) {
    int min_total_time = 1 << 30;
    int best_band_width = -1;
    cout << "Epsilon: " << epsilon << endl;
    for (int band_width: band_widths) {
      cout << "Band width: " << band_width << endl;
      int total_time = 0;
      bool failed = false;
      for (int num_run = 0; num_run < num_runs; num_run++) {
        int elapsed_time = Test1(n, epsilon, dimension, band_width, 20, 1000, values, true, true);
        cout << "Run " << num_run << " done " << elapsed_time << endl;
        if (elapsed_time == -1) {
          failed = true;
          break;
        }
        total_time += elapsed_time;
      }
      if (!failed) {
        if (min_total_time > total_time) {
          min_total_time = total_time;
          best_band_width = band_width;
        }
      }
    }
    if (min_total_time != 1 << 30) {
      cout << "Epsilon: " << epsilon << ", Time: " << static_cast<double>(min_total_time) / num_runs
           << " ms, Band width: " << best_band_width << endl;
    } else {
      cout << "Epsilon " << epsilon << " failed" << endl;
    }
  }
}

void Test(int n, int dimension) {
  vector<uint64_t> best_keys(n);
  vector<uint64_t> best_boundaries(n);
  vector<uint32_t> best_boundaries_lengths(n);
  int min_max_len = 1 << 30;
  for (int k = 0; k < 500; k++) {
    vector<oc::block> blocks = GenerateKeysBlocks(n);
    vector<uint64_t> keys(n);
    for (int i = 0; i < n; i++) {
      keys[i] = blocks[i].get<uint64_t>(0);
    }
    sort(keys.begin(), keys.end());
    int max_len = 0;
    vector<uint64_t> boundaries;
    boundaries.reserve(n / dimension);
    vector<uint32_t> lengths;
    lengths.reserve(n / dimension);
    for (int i = dimension; i < n; i += dimension) {
      uint64_t a = keys[i - 1];
      uint64_t b = keys[i];
      uint64_t boundary = 1;
      for (int j = 0; j < 64; j++) {
        uint64_t mask = static_cast<uint64_t>(1) << (63 - j);
        uint64_t b1 = a & mask;
        uint64_t b2 = b & mask;
        if (b1 != b2) {
          max_len = max(max_len, j + 1);
          boundaries.push_back(boundary);
          lengths.push_back(j + 1);
          break;
        }
        boundary = (boundary << 1) | (b2 >> (63 - j));
      }
    }
    if (max_len < min_max_len) {
      min_max_len = max_len;
      best_keys = keys;
      best_boundaries = boundaries;
      best_boundaries_lengths = lengths;
    }
  }

  int len_per_key = ceil(log2(min_max_len));

  vector<uint64_t> values;
  values.reserve(best_boundaries.size());

  for (int i = 0; i < best_boundaries.size(); i++) {
    uint64_t value = best_boundaries[i];
    values.push_back(value);
  }
  sort(values.begin(), values.end());
  uint64_t max_diff = 0;
  ofstream fout("/usr/local/google/home/jyseo/Desktop/output.txt", ofstream::out);
  for (int i = 1; i < values.size(); i++) {
    uint64_t diff = values[i] - values[i - 1];
    fout << static_cast<uint32_t>(diff);
    max_diff = max(max_diff, diff);
  }
  cout << max_diff << endl;
  fout.close();
}

int BallsIntoBins(int balls, int bins) {
  random_device rd;
  uniform_int_distribution<int> dist(0, bins);
  vector<int> count(bins);
  for (int i = 0; i < balls; i++) {
    count[dist(rd)]++;
  }
  int max_cnt = 0;
  for (int i = 0; i < bins; i++) {
    max_cnt = max(max_cnt, count[i]);
  }
  return max_cnt;
}

void Test4() {
  double epsilon = 0.05;
  int n = 1 << 14;
  int dimension = n;
  uint64_t field_size = 1152921504606830593ULL;
  uint64_t t = field_size;
  uint64_t logt = 60;
  uint64_t element_size = 64;
  int band_width = 100;
  int m = ceil((1 + epsilon) * n);
  BinAllocator2 bin_allocator(1, 2, m, 1);
  uint64_t seed = 0;
  vector<oc::block> keys = GenerateKeysBlocks(n);
  vector<string> raw_values = GenerateRawValues(n, element_size);
  vector<PlaintextValue2> values = ToPlaintextValues2(t, logt, raw_values);
  int max_coeff_count = 0;
  for (int i = 0; i < values.size(); i++) {
    max_coeff_count = max(max_coeff_count, values[i].CoeffCount());
  }

  vector<PlaintextValue2> solution;
  while (true) {
    seed++;
    BandHasher hasher(bin_allocator, m, band_width, band_width, field_size, seed);
    LinearSystem2<PlaintextValue2> linear_system;
    linear_system.Init(m, n, field_size);
    vector<BandRowVector> column_bands;
    column_bands.reserve(n);
    for (const auto& key : keys) {
      column_bands.push_back(hasher.HashToBandRowVector(key));
    }
    vector<PlaintextValue2> new_values(m, PlaintextValue2(max_coeff_count, field_size));
    for (int i = 0; i < column_bands.size(); i++) {
      for (int j = column_bands[i].GetOffset(); j < column_bands[i].GetOffset() + band_width; j++) {
        if (column_bands[i].GetRowElement(j) != 0) {
          new_values[j].AddInPlace(values[i]);
        }
      }
    }

    linear_system.SetColumnVectorsAndValues(column_bands, new_values, m, band_width);

    auto time_begin = chrono::steady_clock::now();
    cout << "Run " << seed << "\n";
    if (linear_system.ReduceToRowEchelonForm(true)) {
      //linear_system.ReduceToRowEchelonForm(true);
      linear_system.SolveLinearSystem();
      solution = linear_system.GetSolution();
      auto time_end = chrono::steady_clock::now();
      std::cout << chrono::duration_cast<chrono::milliseconds>(time_end - time_begin).count() << "\n";
      break;
    }
  }

  for (int i = 0; i <values.size(); i++) {
    if (values[0] != solution[0]) {
      cout << "Wrong answer " << i << endl;
      break;
    }
  }
  cout << "All correct!" << endl;
}

int main() {
  Test4();
  //Experiment();
  //Test3();

  //Test(1 << 20, 1010);

  return 0;
}