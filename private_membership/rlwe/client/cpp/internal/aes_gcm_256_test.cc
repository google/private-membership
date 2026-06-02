#include "private_membership/rlwe/client/cpp/internal/aes_gcm_256.h"

#include <memory>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "shell_encryption/testing/status_matchers.h"
#include "shell_encryption/testing/status_testing.h"
#include "absl/status/status.h"

namespace private_membership {
namespace {

using ::testing::HasSubstr;
using ::testing::Lt;
using ::testing::SizeIs;
using ::rlwe::testing::StatusIs;

const char kTestKey[] = "thequickbrownfoxjumpsoverthelazy";  // 32 bytes
const char kTooLongTestKey[] =
    "thequickbrownfoxjumpsoverthelazydog";  // 35 bytes

TEST(AesGcm256Test, CreateFailsOnKeyOfWrongSize) {
  // Too short.
  EXPECT_THAT(AesGcm256::Create(""),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Key size must be 32 bytes.")));

  // Too long.
  EXPECT_THAT(AesGcm256::Create(kTooLongTestKey),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Key size must be 32 bytes.")));
}

TEST(AesGcm256Test, EncryptProducesCiphertextOfSufficientSize) {
  std::string message = "Hello, World!";

  ASSERT_OK_AND_ASSIGN(std::unique_ptr<AesGcm256> aes_gcm_256,
                       AesGcm256::Create(kTestKey));
  ASSERT_OK_AND_ASSIGN(std::string ciphertext, aes_gcm_256->Encrypt(message));
  EXPECT_THAT(message, SizeIs(Lt(ciphertext.size() - AesGcm256::kNonceSize)));
}

TEST(AesGcm256Test, EncryptAndDecryptRecoversOriginalMessage) {
  std::string message = "Hello, World!";

  ASSERT_OK_AND_ASSIGN(std::unique_ptr<AesGcm256> aes_gcm_256,
                       AesGcm256::Create(kTestKey));
  ASSERT_OK_AND_ASSIGN(std::string ciphertext, aes_gcm_256->Encrypt(message));
  ASSERT_OK_AND_ASSIGN(std::string decrypted, aes_gcm_256->Decrypt(ciphertext));
  EXPECT_EQ(message, decrypted);
}

TEST(AesGcm256Test, DifferentMessagesProduceDifferentCiphertext) {
  std::string message1 = "Hello, World!";
  std::string message2 = "Hi, World!";

  ASSERT_OK_AND_ASSIGN(std::unique_ptr<AesGcm256> aes_gcm_256,
                       AesGcm256::Create(kTestKey));
  ASSERT_OK_AND_ASSIGN(std::string ciphertext1, aes_gcm_256->Encrypt(message1));
  ASSERT_OK_AND_ASSIGN(std::string ciphertext2, aes_gcm_256->Encrypt(message2));
  EXPECT_NE(ciphertext1, ciphertext2);
}

TEST(AesGcm256Test, DifferentKeysProduceDifferentCiphertext) {
  std::string message = "Hello, World!";

  ASSERT_OK_AND_ASSIGN(std::unique_ptr<AesGcm256> aes_gcm_256_1,
                       AesGcm256::Create(kTestKey));
  ASSERT_OK_AND_ASSIGN(std::unique_ptr<AesGcm256> aes_gcm_256_2,
                       AesGcm256::Create(kTestKey));
  ASSERT_OK_AND_ASSIGN(std::string ciphertext1,
                       aes_gcm_256_1->Encrypt(message));
  ASSERT_OK_AND_ASSIGN(std::string ciphertext2,
                       aes_gcm_256_2->Encrypt(message));
  EXPECT_NE(ciphertext1, ciphertext2);
}

TEST(AesGcm256Test, CiphertextWithBadNonceDoesNotDecrypt) {
  std::string message = "Hello, World!";

  ASSERT_OK_AND_ASSIGN(std::unique_ptr<AesGcm256> aes_gcm_256,
                       AesGcm256::Create(kTestKey));
  ASSERT_OK_AND_ASSIGN(std::string ciphertext, aes_gcm_256->Encrypt(message));
  for (int i = 0; i < AesGcm256::kNonceSize; ++i) {
    ciphertext[i] = '\0';
  }

  EXPECT_THAT(aes_gcm_256->Decrypt(ciphertext),
              StatusIs(absl::StatusCode::kInternal, ::testing::_));
}

TEST(AesGcm256Test, MauledCiphertextDoesNotDecrypt) {
  std::string message = "Hello, World!";

  ASSERT_OK_AND_ASSIGN(std::unique_ptr<AesGcm256> aes_gcm_256,
                       AesGcm256::Create(kTestKey));
  ASSERT_OK_AND_ASSIGN(std::string ciphertext, aes_gcm_256->Encrypt(message));
  for (int i = AesGcm256::kNonceSize; i < ciphertext.size(); ++i) {
    ciphertext[i] = '\0';
  }

  EXPECT_THAT(aes_gcm_256->Decrypt(ciphertext),
              StatusIs(absl::StatusCode::kInternal, ::testing::_));
}

}  // namespace
}  // namespace private_membership
