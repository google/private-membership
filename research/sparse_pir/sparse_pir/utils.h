#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>

#include "linear_system_row.h"
#include "plaintext_value.h"

#include "cryptoTools/Common/block.h"
#include "cryptoTools/Crypto/PRNG.h"

using namespace std;

namespace oc = osuCrypto;

vector<oc::block> GenerateKeys(int n) {
  random_device rd;
  uniform_int_distribution<uint64_t> dist;
  uint64_t seed = dist(rd);
  oc::PRNG prng(oc::block(seed, seed));
  vector<oc::block> res(n);
  prng.get<oc::block>(res);
  return res;
}

vector<string> GenerateRawValues(int n, int value_size) {
  vector<string> res(n);
  for (int i = 0; i < n; ++i) {
    string value(value_size, 'v');
    string v = "value-value-value-value-" + to_string(i);
    copy(v.begin(), v.end(), value.begin());
    res[i] = move(value);
  }
  return res;
}

PlaintextValue MergePlaintextValues(const vector<PlaintextValue> &plaintexts, int coeff_length) {
  if (plaintexts.size() == 1) {
    return plaintexts[0];
  }
  uint64_t modulus = plaintexts[0].Modulus();
  vector<uint64_t> res(plaintexts.size() * coeff_length);
  //cout << "Merge size: " << res.size() << endl;
  for (int j = 0; j < plaintexts.size(); ++j) {
    vector<uint64_t> coeffs = plaintexts[j].GetValue();
    if (coeffs.empty()) {
      continue;
    }
    assert(coeffs.size() == coeff_length);
    for (int i = 0; i < coeff_length; ++i) {
      res[coeff_length * j + i] = coeffs[i];
    }
  }
  return PlaintextValue(res, modulus);
}

vector<seal::Plaintext>
ToFlattenedDb(const vector<vector<vector<PlaintextValue>>> &bucketed_db, int poly_deg, const vector<uint64_t> &nvec) {
  int coeff_length = 0;
  for (const auto &a: bucketed_db) {
    for (const auto &b: a) {
      for (const auto &c: b) {
        coeff_length = max(coeff_length, static_cast<int>(c.GetValue().size()));
      }
    }
  }

  vector<vector<PlaintextValue>> bucketed_db2(bucketed_db.size());
  for (int k = 0; k < bucketed_db.size(); ++k) {
    for (int j = 0; j < bucketed_db[k][0].size(); ++j) {
      vector<PlaintextValue> plaintexts;
      for (int i = 0; i < bucketed_db[k].size(); ++i) {
        plaintexts.push_back(bucketed_db[k][i][j]);
      }
      PlaintextValue merged = MergePlaintextValues(plaintexts, coeff_length);
      bucketed_db2[k].push_back(merged);
    }
  }
  coeff_length *= bucketed_db[0].size();

  vector<seal::Plaintext> res;
  /*for (const auto& values : bucketed_db) {
    for (const auto& value : values) {
      seal::Plaintext plaintext = CoeffsToPlaintext(value.GetValue(), poly_deg);
      res.push_back(plaintext);
    }
  }*/
  uint64_t prod = 1;
  for (int i = 1; i < nvec.size(); ++i) {
    prod *= nvec[i];
  }
  vector<uint64_t> padding(poly_deg, 20);
  int cnt = 0;
  for (int i = 0; i < bucketed_db2[0].size(); ++i) {
    for (int j = 0; j < prod; j++) {
      seal::Plaintext plaintext;
      if (j < bucketed_db.size()) {
        plaintext = CoeffsToPlaintext(bucketed_db2[j][i].GetValue(), poly_deg, coeff_length);
      } else {
        VectorToPlaintext(padding, plaintext);
        ++cnt;
      }
      res.push_back(plaintext);
    }
  }
  cout << "Padding: " << cnt << endl;
  return res;
}

int Contains(const vector<uint64_t> &a, const vector<uint64_t> &b) {
  for (int i = 0; i < a.size() - b.size() + 1; ++i) {
    bool eq = true;
    for (int j = 0; j < b.size(); ++j) {
      if (b[j] != a[i + j]) {
        eq = false;
        break;
      }
    }
    if (eq) {
      return i;
    }
  }
  return -1;
}

class GrayCode {
public:
  GrayCode(vector<uint64_t> a) : dimensions_(a) {
    int N = a.size();
    for (int i = 0; i < N; i++) {
      a[i]--;
    }
    vector<int> a1(N);
    vector<int> a2(N);
    vector<int> d(N);
    vector<int> f(N + 1);
    vector<int> x(N);
    int j;

    int count = 0;

    for (int i = 0; i < N; i++) {
      if (a[i] >= 1) {
        x[count] = i;
        count++;
      }
      a1[i] = 0;
      a2[i] = a[i];
      d[i] = 1;
      f[i] = i;
    }
    f[N] = N;
    j = 0;

    value_to_gray_[0] = a1;
    gray_to_value_[a1] = 0;

    int k = 1;
    if (count > 0) {
      while (true) {
        j = f[0];
        f[0] = 0;
        if (j == count) break;
        else {
          a1[x[j]] = a1[x[j]] + d[j];
          a2[x[j]] = a[x[j]] - a1[x[j]];
        }
        if ((a1[x[j]] == 0) || (a1[x[j]] == a[x[j]])) {
          d[j] = -d[j];
          f[j] = f[j + 1];
          f[j + 1] = j + 1;
        }

        value_to_gray_[k] = a1;
        gray_to_value_[a1] = k;
        k++;
      }
    }
  }

  int Size() {
    return value_to_gray_.size();
  }

  vector<int> ToGrayCode(int value) {
    return value_to_gray_[value];
  }

  int ToValue(const vector<int> &gray_code) {
    return gray_to_value_[gray_code];
  }

  int ConvertToInt(const vector<int> &gray_code) {
    int res = 0;
    int d = 1;
    for (int dim: dimensions_) {
      d *= dim;
    }
    for (int i = 0; i < gray_code.size(); i++) {
      d /= dimensions_[i];
      res += d * gray_code[i];
    }
    return res;
  }

private:
  vector<uint64_t> dimensions_;
  map<int, vector<int>> value_to_gray_;
  map<vector<int>, int> gray_to_value_;
};

vector<int> EmbedToHypercube(int n, const vector<uint64_t> &dimensions, bool use_gray_code = true) {
  int d1 = dimensions[0];
  if (n % d1 != 0) {
    cout << "Invalid n, n must be divisible by the first dimension: " << n << " " << d1 << endl;
    exit(0);
  }

  int m = 1;
  for (int d: dimensions) {
    m *= d;
  }

  vector<int> p(m);
  for (int i = 0; i < m; i++) {
    p[i] = i;
  }

  int num_partitions = m / d1;
  vector<bool> should_reverse(num_partitions);
  if (use_gray_code) {
    vector<uint64_t> dimensions2(dimensions.begin() + 1, dimensions.end());
    GrayCode gray_code(dimensions2);
    for (int i = 0; i < num_partitions; i++) {
      vector<int> code = gray_code.ToGrayCode(i);
      int original_index = gray_code.ConvertToInt(code);
      for (int j = 0; j < d1; j++) {
        p[d1 * i + j] = d1 * original_index + j;
      }
      if (i % 2 == 1) {
        should_reverse[original_index] = true;
      }
    }
    /*for (int i = 0; i < num_partitions; i++) {
      if (should_reverse[i]) {
        int last = p[d1 * (i + 1) - 1];
        for (int j = d1 * i, k = 0; j < d1 * (i + 1); j++, k++) {
          p[j] = last - k;
        }
      }
    }*/
  }

  vector<int> p3(n);
  int d = m / d1;
  for (int i = 0; i < n; i++) {
    int j = p[i];
    int partition = j / d1;
    int offset = j % d1;
    if (should_reverse[partition]) {
      offset = d1 - offset - 1;
    }
    p3[i] = d * offset + partition;
  }
  return p3;
}

vector<int> ToGray(int base, int digits, int value) {
  vector<int> baseN(digits);
  int i;
  int val = value;
  for (i = 0; i < digits; i++) {
    baseN[i] = val % base;
    val = val / base;
  }
  int shift = 0;
  vector<int> res(digits);
  while (i--) {
    res[i] = (baseN[i] + shift) % base;
    shift = shift + base - res[i];
  }
  return res;
}

vector<vector<int>> ToQueryVector(const BandRowVector &band_vector, const vector<uint64_t> &nvec) {
  int d1 = nvec[0];

  int offset = band_vector.GetOffset();
  int offset2 = offset % d1;

  int index = offset / d1;
  int len = band_vector.GetLength();
  bool overflow = false;
  if (offset2 + band_vector.GetLength() > d1) {
    len = d1 - offset2;
    overflow = true;
  }
  vector<int> first_dim(d1);
  for (int i = 0; i < len; ++i) {
    first_dim[i + offset2] = band_vector.GetRowElement(offset + i);
  }
  if (overflow && 2 * len != band_vector.GetLength()) {
    exit(0);
  }
  if (index % 2 == 1) {
    reverse(first_dim.begin(), first_dim.end());
  }
  cout << "Overflow: " << overflow << endl;


  vector<uint64_t> dimensions2(nvec.begin() + 1, nvec.end());
  GrayCode gray_code(dimensions2);
  vector<int> gray_code1 = gray_code.ToGrayCode(index);
  vector<int> gray_code2 = gray_code.ToGrayCode(index + 1);

  vector<vector<int>> res;
  res.push_back(first_dim);
  for (int i = 0; i < gray_code1.size(); ++i) {
    vector<int> dim(nvec[i + 1]);
    dim[gray_code1[i]] = 1;
    if (overflow) {
      dim[gray_code2[i]] = 1;
    }
    res.push_back(dim);
  }
  return res;
}

int ComputeDbIndex(int num_vars, int initial_index, const vector<uint64_t> &nvec) {
  int prod = 1;
  for (int i = 1; i < nvec.size(); ++i) {
    prod *= nvec[i];
  }

  int d1 = nvec[0];
  int index = initial_index / d1;
  int offset2 = initial_index % d1;
  int offset = prod * offset2;
  if (index % 2 == 1) {
    offset = prod * (d1 - offset2 - 1);
  }

  vector<int> gray_code = ToGray(4, nvec.size() - 1, index);
  int index2 = 0;
  for (int i = 0; i < gray_code.size(); ++i) {
    index2 = 4 * index2 + gray_code[i];
    //cout << gray_code[i];
  }
  //cout << endl;
  return offset + index2;
}

vector<seal::Plaintext> ToFlattenedDb3(const vector<PlaintextValue> &db, int poly_deg, const vector<uint64_t> &nvec) {
  int coeff_length = 0;
  for (const auto &a: db) {
    coeff_length = max(coeff_length, static_cast<int>(a.GetValue().size()));
  }

  /*for (const auto& values : bucketed_db) {
    for (const auto& value : values) {
      seal::Plaintext plaintext = CoeffsToPlaintext(value.GetValue(), poly_deg);
      res.push_back(plaintext);
    }
  }*/
  uint64_t prod = 1;
  for (int i = 1; i < nvec.size(); ++i) {
    prod *= nvec[i];
  }
  vector<uint64_t> padding(poly_deg, 1);
  seal::Plaintext plaintext;
  VectorToPlaintext(padding, plaintext);
  vector<seal::Plaintext> res(nvec[0] * prod, plaintext);
  vector<bool> filled(nvec[0] * prod);

  cout << "Coeff length: " << coeff_length << endl;
  vector<int> p = EmbedToHypercube(db.size(), nvec, true);
  for (int i = 0; i < db.size(); ++i) {
    int index = p[i];
    if (filled[index]) {
      cout << "Filled!" << endl;
    }
    res[index] = CoeffsToPlaintext(db[i].GetValue(), poly_deg, coeff_length);
    filled[index] = true;
  }
  /*for (int i = 0; i < filled.size(); ++i) {
    if (!filled[i]) {
      cout << "Not filled!" << endl;
    }
  }*/
  return res;
}

#endif //UTILS_H
