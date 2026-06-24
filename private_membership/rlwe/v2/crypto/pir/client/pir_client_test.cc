#include "private_membership/rlwe/v2/crypto/pir/client/pir_client.h"

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "private_membership/rlwe/v2/crypto/encryption.h"
#include "private_membership/rlwe/v2/crypto/pir/pir_params.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "shell_encryption/prng/chacha_prng.h"
#include "shell_encryption/prng/prng.h"

namespace private_membership {
namespace rlwe {
namespace v2 {
namespace {

using ::testing::HasSubstr;
using ::testing::Test;
using ::rlwe::testing::StatusIs;

constexpr int kDegree = 1024;
constexpr int kVariance = 8;
constexpr int kLogDigit0 = 2;
constexpr int kNumDigits0 = 4;
constexpr int kLogDigit1 = 3;
constexpr int kNumDigits1 = 3;

template <typename CoeffType>
class PirClientTest : public Test {
 protected:
  void SetUp() override {
    ASSERT_OK_AND_ASSIGN(
        auto rlwe_params,
        RlweParams<CoeffType>::Create(kDegree, GetModulus(), kVariance));
    rlwe_params_ =
        std::make_unique<RlweParams<CoeffType>>(std::move(rlwe_params));
    ASSERT_OK_AND_ASSIGN(a_component_prng_, this->CreatePrng());
    ASSERT_OK_AND_ASSIGN(prng_, this->CreatePrng());
  }

  CoeffType GetModulus() const {
    if (sizeof(CoeffType) == 4) {
      return 0;  // uint32_t uses implicit modulus 2^32 (0)
    } else {
      return 0;  // uint64_t uses implicit modulus 2^64 (0)
    }
  }

  CoeffType GetPlaintextModulus() const { return 256; }

  GadgetParams GetGadgetParams0() const {
    return GadgetParams{/*log_digit=*/kLogDigit0, /*num_digits=*/kNumDigits0};
  }

  GadgetParams GetGadgetParams1() const {
    return GadgetParams{/*log_digit=*/kLogDigit1, /*num_digits=*/kNumDigits1};
  }

  absl::StatusOr<PirParams<CoeffType>> CreateValidPirParams() const {
    int num_entries = 4096;
    int interpolation_degree = 2;
    return PirParams<CoeffType>::Create(*rlwe_params_, GetPlaintextModulus(),
                                        GetGadgetParams0(), GetGadgetParams1(),
                                        num_entries, interpolation_degree);
  }

  absl::StatusOr<std::unique_ptr<::rlwe::SecurePrng>> CreatePrng() const {
    RLWE_ASSIGN_OR_RETURN(auto prng_seed, ::rlwe::ChaChaPrng::GenerateSeed());
    return ::rlwe::ChaChaPrng::Create(prng_seed);
  }

  std::unique_ptr<RlweParams<CoeffType>> rlwe_params_;
  std::unique_ptr<::rlwe::SecurePrng> a_component_prng_;
  std::unique_ptr<::rlwe::SecurePrng> prng_;
};

using CoeffTypes = ::testing::Types<uint32_t, uint64_t>;
TYPED_TEST_SUITE(PirClientTest, CoeffTypes);

TYPED_TEST(PirClientTest, CreateFailsWithNullPrng) {
  ASSERT_OK_AND_ASSIGN(auto params, this->CreateValidPirParams());
  EXPECT_THAT(PirClient<TypeParam>::Create(params, nullptr, nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("a_component_prng must not be null")));

  EXPECT_THAT(PirClient<TypeParam>::Create(
                  params, std::move(this->a_component_prng_), nullptr),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("prng must not be null")));
}

TYPED_TEST(PirClientTest, CreateRequestFailsWithIndexOutOfBounds) {
  ASSERT_OK_AND_ASSIGN(auto params, this->CreateValidPirParams());
  ASSERT_OK_AND_ASSIGN(
      auto client,
      PirClient<TypeParam>::Create(params, std::move(this->a_component_prng_),
                                   std::move(this->prng_)));

  EXPECT_THAT(client->CreateRequest(-1),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Index out of bounds")));

  EXPECT_THAT(client->CreateRequest(params.NumEntries()),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Index out of bounds")));
}

TYPED_TEST(PirClientTest, CreateRequestGeneratesValidCiphertexts) {
  ASSERT_OK_AND_ASSIGN(auto params, this->CreateValidPirParams());

  std::vector<int> test_indices = {0, 1, 35, 100, params.NumEntries() - 1};

  for (int index : test_indices) {
    ASSERT_OK_AND_ASSIGN(auto prng_seed, ::rlwe::ChaChaPrng::GenerateSeed());
    ASSERT_OK_AND_ASSIGN(auto a_component_prng,
                         ::rlwe::ChaChaPrng::Create(prng_seed));
    ASSERT_OK_AND_ASSIGN(auto a_component_prng_clone,
                         ::rlwe::ChaChaPrng::Create(prng_seed));
    ASSERT_OK_AND_ASSIGN(auto prng, this->CreatePrng());

    ASSERT_OK_AND_ASSIGN(
        auto client, PirClient<TypeParam>::Create(
                         params, std::move(a_component_prng), std::move(prng)));

    ASSERT_OK_AND_ASSIGN(auto request, client->CreateRequest(index));

    const auto& rlwe_params = params.RlweParameters();
    int d = rlwe_params.Degree();
    int t = params.InterpolationDegree();

    int expected_first_dim = params.NumEntries() / t / d;
    EXPECT_EQ(request.first_dimension_query.size(), expected_first_dim);

    int expected_second_dim = 2 * params.PolyEvalGadgetParams().num_digits;
    EXPECT_EQ(request.second_dimension_query.size(), expected_second_dim);

    int expected_packing = 2 * params.PackGadgetParams().num_digits;
    EXPECT_EQ(request.packing_key.size(), expected_packing);

    int num_samples = params.NumEntries() / t / d;
    ASSERT_OK_AND_ASSIGN(
        auto a_first,
        SampleAComponents(rlwe_params, num_samples, *a_component_prng_clone));

    int entries_per_row = t;
    int row_index = index / entries_per_row;
    int expected_i = row_index / d;
    int expected_j = row_index % d;

    for (int idx = 0; idx < request.first_dimension_query.size(); ++idx) {
      RlweCiphertext<TypeParam> ciphertext{
          std::move(a_first[idx]),
          std::move(request.first_dimension_query[idx])};
      ASSERT_OK_AND_ASSIGN(auto plaintext,
                           Decrypt(rlwe_params, params.PlaintextModulus(),
                                   ciphertext, client->GetSecretKey()));

      std::vector<TypeParam> expected_message(d, 0);
      if (idx == expected_i) {
        expected_message[expected_j] = 1;
      }
      EXPECT_EQ(plaintext, expected_message);
    }
  }
}

}  // namespace
}  // namespace v2
}  // namespace rlwe
}  // namespace private_membership
