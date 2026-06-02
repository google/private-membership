// Copyright 2020 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "private_membership/rlwe/client/cpp/internal/crypto_utils.h"

#include <cstdint>
#include <memory>
#include <string>

#include "private_join_and_compute/crypto/context.h"
#include "private_join_and_compute/crypto/ec_commutative_cipher.h"
#include "private_membership/rlwe/client/cpp/internal/aes_ctr_256_no_auth.h"
#include "private_membership/rlwe/client/cpp/internal/aes_ctr_256_with_fixed_iv.h"
#include "private_membership/rlwe/client/cpp/internal/aes_gcm_256.h"
#include "private_membership/rlwe/client/cpp/internal/cipher.h"
#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include <openssl/nid.h>
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"

namespace private_membership {
namespace {

using ::rlwe::testing::StatusIs;
using ::testing::Gt;
using ::testing::HasSubstr;
using ::testing::SizeIs;

constexpr int kCurveId = NID_X9_62_prime256v1;

TEST(CryptoUtilsTest, PadMaxByteLengthSmallerThanInputByteLength) {
  EXPECT_THAT(
      Pad("input", 4),
      StatusIs(
          absl::StatusCode::kInvalidArgument,
          HasSubstr("max_byte_length smaller than the input bytes length.")));
}

TEST(CryptoUtilsTest, PadSuccess) {
  std::string expected_output("\x05\x00\x00\x00input\x00\x00", 11);
  ASSERT_OK_AND_ASSIGN(auto x, Pad("input", 7));
  EXPECT_EQ(x, expected_output);
}

TEST(CryptoUtilsTest, UnpadByteLengthSmallerThan4) {
  EXPECT_THAT(Unpad("abc"),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Invalid bytes does not encode length.")));
}

TEST(CryptoUtilsTest, UnpadIncorrectByteLength) {
  std::string input("\x06\x00\x00\x00input", 9);
  EXPECT_THAT(Unpad(input), StatusIs(absl::StatusCode::kInvalidArgument,
                                     HasSubstr("Incorrect bytes length.")));
}

TEST(CryptoUtilsTest, UnpadSuccess) {
  std::string input("\x05\x00\x00\x00input\x00\x00", 11);
  ASSERT_OK_AND_ASSIGN(auto x, Unpad(input));
  EXPECT_EQ(x, "input");
}

TEST(CryptoUtilsTest, HashEncryptedId) {
  ::private_join_and_compute::Context ctx;
  std::string id = "id";
  auto hash1 = HashEncryptedId(id, &ctx);
  auto hash2 = HashEncryptedId(id, &ctx);
  EXPECT_EQ(hash1, hash2);
  EXPECT_NE(hash1, id);
}

class EncryptDecryptValueTest
    : public ::testing::TestWithParam<Cipher::Scheme> {
 protected:
  absl::StatusOr<std::unique_ptr<Cipher>> CreateCipher(
      Cipher::Scheme encryption_scheme, absl::string_view key) {
    switch (encryption_scheme) {
      case Cipher::Scheme::kAesCtr256WithFixedIv:
        return AesCtr256WithFixedIV::Create(key);
      case Cipher::Scheme::kAesCtr256NoAuth:
        return AesCtr256NoAuth::Create(key);
      case Cipher::Scheme::kAesGcm256:
        return AesGcm256::Create(key);
        default:
        return absl::InvalidArgumentError(
            absl::StrCat("Unsupported encryption scheme: ", encryption_scheme));
    }
  }
};

INSTANTIATE_TEST_SUITE_P(
    EncryptValueTest, EncryptDecryptValueTest,
    ::testing::Values(Cipher::Scheme::kAesCtr256WithFixedIv,
                      Cipher::Scheme::kAesCtr256NoAuth,
                      Cipher::Scheme::kAesGcm256));

TEST_P(EncryptDecryptValueTest, EncryptValue) {
  ::private_join_and_compute::Context ctx;
  std::string encrypted_id = "encrypt_value_id";
  std::string value = "value";
  uint32_t max_value_byte_length = 32;
  Cipher::Scheme encryption_scheme = GetParam();
  ASSERT_OK_AND_ASSIGN(std::string ciphertext,
                       EncryptValue(encrypted_id, value, max_value_byte_length,
                                    encryption_scheme, &ctx));
  EXPECT_NE(ciphertext, value);

  ASSERT_OK_AND_ASSIGN(auto value_encryption_key,
                       GetValueEncryptionKey(encrypted_id, &ctx));

  ASSERT_OK_AND_ASSIGN(
      auto x, DecryptValue(encrypted_id, ciphertext, encryption_scheme, &ctx));
  EXPECT_EQ(x, value);

  ASSERT_OK_AND_ASSIGN(x, DecryptValueWithKey(value_encryption_key, ciphertext,
                                              encryption_scheme, &ctx));
  EXPECT_EQ(x, value);
}

TEST_P(EncryptDecryptValueTest, EncryptValueCiphertextLength) {
  ::private_join_and_compute::Context ctx;
  std::string encrypted_id = "encrypt_value_id_length";
  uint32_t max_value_byte_length = 16;
  Cipher::Scheme encryption_scheme = GetParam();
  ASSERT_OK_AND_ASSIGN(
      std::string ciphertext_empty_string,
      EncryptValue(encrypted_id, /*value=*/"", max_value_byte_length,
                   encryption_scheme, &ctx));
  for (int i = 1; i <= max_value_byte_length; ++i) {
    std::string value(i, 0);
    ASSERT_OK_AND_ASSIGN(
        std::string ciphertext,
        EncryptValue(encrypted_id, value, max_value_byte_length,
                     encryption_scheme, &ctx));
    if (encryption_scheme == Cipher::Scheme::kAesCtr256NoAuth) {
      // Because there is no padding, the size should be strictly greater
      // than the size of the empty string.
      EXPECT_THAT(ciphertext, SizeIs(Gt(ciphertext_empty_string.length())));
    } else {
      EXPECT_THAT(ciphertext, SizeIs(ciphertext_empty_string.length()));
    }
  }
}

TEST_P(EncryptDecryptValueTest, EncryptInvalidTooLongValue) {
  ::private_join_and_compute::Context ctx;
  std::string encrypted_id = "encrypt_invalid_value_id";
  uint32_t max_value_byte_length = 8;
  std::string value(max_value_byte_length + 1, 'a');
  Cipher::Scheme encryption_scheme = GetParam();
  EXPECT_THAT(EncryptValue(encrypted_id, value, max_value_byte_length,
                           encryption_scheme, &ctx),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("larger than maximum value byte length")));
}

TEST_P(EncryptDecryptValueTest, DecryptInvalidCiphertextWithoutLengthEncoding) {
  ::private_join_and_compute::Context ctx;
  std::string encrypted_id = "no_length_encoding";
  std::string invalid_plaintext = "abc";
  Cipher::Scheme encryption_scheme = GetParam();
  ASSERT_OK_AND_ASSIGN(auto value_encryption_key,
                       GetValueEncryptionKey(encrypted_id, &ctx));
  ASSERT_OK_AND_ASSIGN(auto cipher,
                       CreateCipher(encryption_scheme, value_encryption_key));
  ASSERT_OK_AND_ASSIGN(auto encrypted_value,
                       cipher->Encrypt(invalid_plaintext));
  absl::StatusOr<std::string> decrypted_value =
      DecryptValue(encrypted_id, encrypted_value, encryption_scheme, &ctx);
  if (encryption_scheme == Cipher::Scheme::kAesCtr256NoAuth) {
    // Because no auth scheme is used, there is no length encoding at all,
    // so it will decrypt successfully to the "invalid plaintext".
    ASSERT_OK(decrypted_value);
    EXPECT_THAT(*decrypted_value, testing::Eq("abc"));
  } else {
    EXPECT_THAT(decrypted_value, StatusIs(absl::StatusCode::kInvalidArgument,
                                          HasSubstr("does not encode length")));
  }
}

TEST_P(EncryptDecryptValueTest, DecryptInvalidCiphertextWithIncorrectLength) {
  ::private_join_and_compute::Context ctx;
  std::string encrypted_id = "incorrect_length";
  std::string invalid_plaintext = "aaaa";
  Cipher::Scheme encryption_scheme = GetParam();
  ASSERT_OK_AND_ASSIGN(auto value_encryption_key,
                       GetValueEncryptionKey(encrypted_id, &ctx));
  ASSERT_OK_AND_ASSIGN(auto cipher,
                       CreateCipher(encryption_scheme, value_encryption_key));
  ASSERT_OK_AND_ASSIGN(auto encrypted_value,
                       cipher->Encrypt(invalid_plaintext));
  absl::StatusOr<std::string> decrypted_value =
      DecryptValue(encrypted_id, encrypted_value, encryption_scheme, &ctx);
  if (encryption_scheme == Cipher::Scheme::kAesCtr256NoAuth) {
    // Because no auth scheme is used, there is no length encoding at all,
    // so it will decrypt successfully to the "invalid plaintext".
    ASSERT_OK(decrypted_value);
    EXPECT_THAT(*decrypted_value, testing::Eq("aaaa"));
  } else {
    EXPECT_THAT(decrypted_value, StatusIs(absl::StatusCode::kInvalidArgument,
                                          HasSubstr("Incorrect bytes length")));
  }
}

TEST(NoAuth, AesCtr256NoAuthDecryptionDoesNotFail) {
  ::private_join_and_compute::Context ctx;
  EXPECT_OK(DecryptValue("random_encrypted_id", "random_encrypted_value",
                         Cipher::Scheme::kAesCtr256NoAuth, &ctx));
}

TEST(CryptoUtilsTest, HashEncryptedIdRegression) {
  ::private_join_and_compute::Context ctx;
  std::string id = "id";
  std::string hash = HashEncryptedId(id, &ctx);

  constexpr unsigned char expected_hash[] = {
      0xd0, 0xd5, 0x50, 0xc1, 0x88, 0xa8, 0xcd, 0x8a, 0x4,  0x97, 0x26,
      0x86, 0x90, 0x4,  0xe1, 0xbc, 0xbd, 0x9d, 0xba, 0x19, 0x3f, 0x2d,
      0x36, 0xa8, 0xb5, 0x87, 0x95, 0xe8, 0xe2, 0x57, 0x1,  0x76,
  };
  EXPECT_EQ(hash, std::string(reinterpret_cast<const char*>(expected_hash),
                              hash.length()));
}

TEST(CryptoUtilsTest, GetValueEncryptionKeyRegression) {
  ::private_join_and_compute::Context ctx;
  std::string id = "id";
  ASSERT_OK_AND_ASSIGN(std::string hash, GetValueEncryptionKey(id, &ctx));

  constexpr unsigned char expected_hash[] = {
      0x2a, 0x3b, 0xe5, 0xe4, 0x33, 0xd0, 0x85, 0x73, 0xdd, 0x76, 0xad,
      0xd1, 0x86, 0x33, 0x30, 0xac, 0x5b, 0x7f, 0x81, 0xbd, 0x4e, 0xdd,
      0xdb, 0x80, 0xc3, 0x51, 0xb5, 0xad, 0x0,  0x38, 0xca, 0xbd,
  };
  EXPECT_EQ(hash, std::string(reinterpret_cast<const char*>(expected_hash),
                              hash.length()));
}

TEST(CryptoUtilsTest, CorrectlyTruncate) {
  // \xa5 = b10100101
  EXPECT_THAT(Truncate("\xa5\xff", -1),
              StatusIs(absl::StatusCode::kInvalidArgument, ::testing::_));
  ASSERT_OK_AND_ASSIGN(std::string zero, Truncate("\xa5\xff", 0));
  EXPECT_EQ(zero, "");
  ASSERT_OK_AND_ASSIGN(std::string one, Truncate("\xa5\xff", 1));
  EXPECT_EQ(one, "\x80");
  ASSERT_OK_AND_ASSIGN(std::string two, Truncate("\xa5\xff", 2));
  EXPECT_EQ(two, "\x80");
  ASSERT_OK_AND_ASSIGN(std::string three, Truncate("\xa5\xff", 3));
  EXPECT_EQ(three, "\xa0");
  ASSERT_OK_AND_ASSIGN(std::string four, Truncate("\xa5\xff", 4));
  EXPECT_EQ(four, "\xa0");
  ASSERT_OK_AND_ASSIGN(std::string five, Truncate("\xa5\xff", 5));
  EXPECT_EQ(five, "\xa0");
  ASSERT_OK_AND_ASSIGN(std::string six, Truncate("\xa5\xff", 6));
  EXPECT_EQ(six, "\xa4");
  ASSERT_OK_AND_ASSIGN(std::string seven, Truncate("\xa5\xff", 7));
  EXPECT_EQ(seven, "\xa4");
  ASSERT_OK_AND_ASSIGN(std::string eight, Truncate("\xa5\xff", 8));
  EXPECT_EQ(eight, "\xa5");
  ASSERT_OK_AND_ASSIGN(std::string nine, Truncate("\xa5\xff", 9));
  EXPECT_EQ(nine, "\xa5\x80");
  ASSERT_OK_AND_ASSIGN(std::string ten, Truncate("\xa5\xff", 10));
  EXPECT_EQ(ten, "\xa5\xc0");
  ASSERT_OK_AND_ASSIGN(std::string eleven, Truncate("\xa5\xff", 11));
  EXPECT_EQ(eleven, "\xa5\xe0");
  ASSERT_OK_AND_ASSIGN(std::string twelve, Truncate("\xa5\xff", 12));
  EXPECT_EQ(twelve, "\xa5\xf0");
  ASSERT_OK_AND_ASSIGN(std::string thirteen, Truncate("\xa5\xff", 13));
  EXPECT_EQ(thirteen, "\xa5\xf8");
  ASSERT_OK_AND_ASSIGN(std::string fourteen, Truncate("\xa5\xff", 14));
  EXPECT_EQ(fourteen, "\xa5\xfc");
  ASSERT_OK_AND_ASSIGN(std::string fifteen, Truncate("\xa5\xff", 15));
  EXPECT_EQ(fifteen, "\xa5\xfe");
  ASSERT_OK_AND_ASSIGN(std::string sixteen, Truncate("\xa5\xff", 16));
  EXPECT_EQ(sixteen, "\xa5\xff");
  EXPECT_THAT(Truncate("\xa5\xff", 17),
              StatusIs(absl::StatusCode::kInvalidArgument, ::testing::_));
}

}  // namespace
}  // namespace private_membership
