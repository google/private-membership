#include "private_membership/rlwe/client/cpp/internal/oprf_utils.h"

#include <sys/stat.h>

#include <memory>
#include <string>

#include "private_join_and_compute/crypto/ec_commutative_cipher.h"
#include "private_membership/rlwe/client/proto/private_membership_rlwe.pb.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "absl/status/status.h"

namespace private_membership {
namespace rlwe {
namespace {


using ::rlwe::testing::StatusIs;

TEST(OprfUtilsTest, CreateOprfCipher) {
  EXPECT_OK(
      CreateOprfCipher(OprfCurveAndHashType::kPrime256v1CurveIdWithSha256));
  EXPECT_OK(
      CreateOprfCipher(OprfCurveAndHashType::kPrime256v1CurveIdWithSswuRo));
}

TEST(OprfUtilsTest, CreateOprfCipherWithUndefinedCurveAndHashType) {
  EXPECT_THAT(
      CreateOprfCipher(OprfCurveAndHashType::kUndefinedOprfCurveAndHashType),
      StatusIs(absl::StatusCode::kInvalidArgument, ::testing::_));
}

TEST(OprfUtilsTest, CreateOprfCipherWithCorrectKeys) {
  ASSERT_OK_AND_ASSIGN(
      std::unique_ptr<::private_join_and_compute::ECCommutativeCipher> cipher1,
      CreateOprfCipher(OprfCurveAndHashType::kPrime256v1CurveIdWithSswuRo));
  ASSERT_OK_AND_ASSIGN(
      std::unique_ptr<::private_join_and_compute::ECCommutativeCipher> cipher2,
      CreateOprfCipher(OprfCurveAndHashType::kPrime256v1CurveIdWithSswuRo));
  ASSERT_OK_AND_ASSIGN(
      std::unique_ptr<::private_join_and_compute::ECCommutativeCipher> cipher3,
      CreateOprfCipher(OprfCurveAndHashType::kPrime256v1CurveIdWithSswuRo,
                       cipher1->GetPrivateKeyBytes()));
  ASSERT_OK_AND_ASSIGN(std::string encrypted_id1, cipher1->Encrypt("test"));
  ASSERT_OK_AND_ASSIGN(std::string encrypted_id2, cipher2->Encrypt("test"));
  ASSERT_OK_AND_ASSIGN(std::string encrypted_id3, cipher3->Encrypt("test"));

  EXPECT_EQ(encrypted_id1, encrypted_id3);
  EXPECT_NE(encrypted_id1, encrypted_id2);
}

TEST(OprfUtilsTest, CreateOprfCipherWithCorrectCurveAndHashType) {
  ASSERT_OK_AND_ASSIGN(
      std::unique_ptr<::private_join_and_compute::ECCommutativeCipher> cipher1,
      CreateOprfCipher(OprfCurveAndHashType::kPrime256v1CurveIdWithSswuRo));
  ASSERT_OK_AND_ASSIGN(
      std::unique_ptr<::private_join_and_compute::ECCommutativeCipher> cipher2,
      CreateOprfCipher(OprfCurveAndHashType::kPrime256v1CurveIdWithSha256,
                       cipher1->GetPrivateKeyBytes()));
  ASSERT_OK_AND_ASSIGN(
      std::unique_ptr<::private_join_and_compute::ECCommutativeCipher> cipher3,
      CreateOprfCipher(OprfCurveAndHashType::kPrime256v1CurveIdWithSswuRo,
                       cipher1->GetPrivateKeyBytes()));
  ASSERT_OK_AND_ASSIGN(std::string encrypted_id1, cipher1->Encrypt("test"));
  ASSERT_OK_AND_ASSIGN(std::string encrypted_id2, cipher2->Encrypt("test"));
  ASSERT_OK_AND_ASSIGN(std::string encrypted_id3, cipher3->Encrypt("test"));

  EXPECT_EQ(encrypted_id1, encrypted_id3);
  EXPECT_NE(encrypted_id1, encrypted_id2);
}

TEST(OprfUtilsTest, ToInternalEnum) {
  ASSERT_OK_AND_ASSIGN(auto val1, ToInternalEnum(OprfParams::OPRF_CONFIG_NID_X9_62_PRIME256V1_WITH_SHA256));
  EXPECT_EQ(val1, OprfCurveAndHashType::kPrime256v1CurveIdWithSha256);
  ASSERT_OK_AND_ASSIGN(auto val2, ToInternalEnum(OprfParams::OPRF_CONFIG_NID_X9_62_PRIME256V1_WITH_SSWU_RO));
  EXPECT_EQ(val2, OprfCurveAndHashType::kPrime256v1CurveIdWithSswuRo);
  ASSERT_OK_AND_ASSIGN(auto val3, ToInternalEnum(OprfParams::OPRF_CONFIG_UNSPECIFIED));
  EXPECT_EQ(val3, OprfCurveAndHashType::kUndefinedOprfCurveAndHashType);
  EXPECT_THAT(ToInternalEnum(static_cast<OprfParams>(-1)),
              StatusIs(absl::StatusCode::kInvalidArgument, ::testing::_));
}

}  // namespace
}  // namespace rlwe
}  // namespace private_membership
