#include <random>

#include "../bandokvs/band_okvs.h"
#include "../bandokvs/uint.h"

using namespace band_okvs;

void BandOkvsTest() {
  double epsilon = 0.05;
  int n = 1 << 20;
  int m = static_cast<int>((1 + epsilon) * n);
  int band_length = 383;

  std::random_device rd;
  std::uniform_int_distribution<uint64_t> dist;

  oc::PRNG prng(oc::block(dist(rd), dist(rd)));
  std::vector<oc::block> keys(n);
  prng.get<oc::block>(keys);

  using band_type = band_okvs::uint<3>;

  std::vector<oc::block> values(n);
  GenRandomValuesBlocks(values, dist(rd), dist(rd));

  BandOkvs okvs;
  okvs.Init(n, m, band_length, oc::block(dist(rd), dist(rd)));

  std::vector<oc::block> out(okvs.Size());

  std::cout << "Encoding..." << std::endl;
  StartTimer();
  if (!okvs.Encode<band_type, oc::block>(keys.data(),
                                         values.data(),
                                         out.data())) {
    std::cout << "Failed to encode!" << std::endl;
    exit(0);
  }
  EndTimer();
  std::cout << "Encoding done." << std::endl;
  std::cout << "Elapsed time: " << GetElapsedTime().count() << std::endl;

  std::vector<oc::block> decoded(n);
  std::cout << "Decoding..." << std::endl;
  StartTimer();
  okvs.Decode<band_type>(keys.data(), out.data(), decoded.data());
  EndTimer();
  std::cout << "Decoding done" << std::endl;
  std::cout << "Elapsed time: " << GetElapsedTime().count() << std::endl;

  for (int i = 0; i < n; i++) {
    if (values[i] != decoded[i]) {
      std::cout << "Wrong result " << i << " " << values[i] << " " << decoded[i]
                << std::endl;
      exit(0);
    }
  }
  std::cout << "All correct!" << std::endl;
}

int main() {
  BandOkvsTest();
  return 0;
}
