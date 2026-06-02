#include "private_membership/rlwe/client/cpp/internal/aes_ctr_256_no_auth.h"

#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "absl/status/status.h"
#include "absl/strings/string_view.h"

namespace private_membership {
namespace {

using ::testing::IsEmpty;
using ::testing::Ne;

using ::rlwe::testing::StatusIs;

std::string GetRandomKey() {
  // https://xkcd.com/221/
  return "12345678901234567890123456789012";
}

constexpr absl::string_view kPlaintext =
    R"(" !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~")";

TEST(AesCtr256NoAuthTest, EncryptDecryptEmptyString) {
  ASSERT_OK_AND_ASSIGN(auto cipher, AesCtr256NoAuth::Create(GetRandomKey()));

  ASSERT_OK_AND_ASSIGN(const std::string ciphertext, cipher->Encrypt(""));
  ASSERT_OK_AND_ASSIGN(auto plaintext, cipher->Decrypt(ciphertext));
  EXPECT_THAT(plaintext, IsEmpty());
}

TEST(AesCtr256NoAuthTest, EncryptDecryptPartialBlock) {
  ASSERT_OK_AND_ASSIGN(auto cipher, AesCtr256NoAuth::Create(GetRandomKey()));

  for (int extra_bytes = 0; extra_bytes < kPlaintext.size(); ++extra_bytes) {
    if (extra_bytes % AesCtr256NoAuth::kBlockSize == 0) continue;
    ASSERT_OK_AND_ASSIGN(const std::string ciphertext,
                         cipher->Encrypt(kPlaintext.substr(0, extra_bytes)));
    ASSERT_OK_AND_ASSIGN(auto plaintext, cipher->Decrypt(ciphertext));
    EXPECT_THAT(plaintext, kPlaintext.substr(0, extra_bytes));
  }
}

TEST(AesCtr256NoAuthTest, EncryptDecryptFullBlock) {
  ASSERT_OK_AND_ASSIGN(auto cipher, AesCtr256NoAuth::Create(GetRandomKey()));

  for (int extra_bytes = 0; extra_bytes < kPlaintext.size(); ++extra_bytes) {
    if (extra_bytes % AesCtr256NoAuth::kBlockSize != 0) continue;
    ASSERT_OK_AND_ASSIGN(const std::string ciphertext,
                         cipher->Encrypt(kPlaintext.substr(0, extra_bytes)));
    ASSERT_OK_AND_ASSIGN(auto plaintext, cipher->Decrypt(ciphertext));
    EXPECT_THAT(plaintext, kPlaintext.substr(0, extra_bytes));
  }
}

TEST(AesCtr256NoAuthTest, RandomIV) {
  ASSERT_OK_AND_ASSIGN(auto cipher, AesCtr256NoAuth::Create(GetRandomKey()));

  ASSERT_OK_AND_ASSIGN(const std::string ciphertext1,
                       cipher->Encrypt(kPlaintext));
  ASSERT_OK_AND_ASSIGN(const std::string ciphertext2,
                       cipher->Encrypt(kPlaintext));
  EXPECT_NE(ciphertext1, ciphertext2);
}

TEST(AesCtr256NoAuthTest, InvalidKey) {
  EXPECT_THAT(
      AesCtr256NoAuth::Create(""),
      StatusIs(absl::StatusCode::kInvalidArgument, "Key size is invalid."));
  EXPECT_THAT(
      AesCtr256NoAuth::Create("1234567890123456789012345678901"),
      StatusIs(absl::StatusCode::kInvalidArgument, "Key size is invalid."));
  EXPECT_THAT(
      AesCtr256NoAuth::Create("123456789012345678901234567890123"),
      StatusIs(absl::StatusCode::kInvalidArgument, "Key size is invalid."));
}

TEST(AesCtr256NoAuthTest,
     DecryptingWithDifferentKeyDoesNotFailButReturnsDifferentPlaintext) {
  ASSERT_OK_AND_ASSIGN(auto cipher1, AesCtr256NoAuth::Create(
                                         "12345678901234567890123456789012"));
  ASSERT_OK_AND_ASSIGN(auto cipher2, AesCtr256NoAuth::Create(
                                         "abcdefghijklmnopqrstuvwxyzabcdef"));

  ASSERT_OK_AND_ASSIGN(const std::string ciphertext,
                       cipher1->Encrypt(kPlaintext));
  ASSERT_OK_AND_ASSIGN(auto plaintext, cipher2->Decrypt(ciphertext));
  EXPECT_THAT(plaintext, Ne(kPlaintext));
}

TEST(AesCtr256NoAuthTest, EncryptingNullPlaintextDoesNotFail) {
  ASSERT_OK_AND_ASSIGN(auto cipher, AesCtr256NoAuth::Create(GetRandomKey()));

  EXPECT_OK(cipher->Encrypt(
      absl::NullSafeStringView(static_cast<const char*>(nullptr))));
}

TEST(AesCtr256NoAuthTest, DecryptingInvalidCiphertextDoesNotFail) {
  ASSERT_OK_AND_ASSIGN(auto cipher, AesCtr256NoAuth::Create(GetRandomKey()));

  EXPECT_OK(cipher->Decrypt("invalid_ciphertext"));
}

TEST(AesCtr256NoAuthTest, DecryptingTooShortCiphertextDoesFail) {
  ASSERT_OK_AND_ASSIGN(auto cipher, AesCtr256NoAuth::Create(GetRandomKey()));

  EXPECT_THAT(cipher->Decrypt("c"), StatusIs(absl::StatusCode::kInternal,
                                             "Ciphertext is too short."));
  EXPECT_THAT(cipher->Decrypt(""), StatusIs(absl::StatusCode::kInternal,
                                            "Ciphertext is too short."));
}

}  // namespace
}  // namespace private_membership
