#ifndef PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_CIPHER_H_
#define PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_CIPHER_H_

#include <string>

#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"

namespace private_membership {

// Abstract interface for a cipher.
class Cipher {
 public:
  // LINT.IfChange
  // **IMPORTANT**: Keep the integer value in sync with the enum in
  // google3/privacy/private_membership/proto/private_membership.proto
  enum class Scheme {
    kEncryptionSchemeUndefined = 0,
    kAesCtr256WithFixedIv = 1,
    kAesCtr256NoAuth = 2,
    kAesGcm256 = 3,
    kAesSiv = 4,
  };
  // LINT.ThenChange(
  //    //depot/google3/privacy/private_membership/proto/private_membership.proto,
  //    //depot/google3/googlemac/iPhone/Firebase/Analytics/AppMeasurement/Nano/PrivateSetMembership/Protos/private_membership.proto)

  virtual ~Cipher() = default;

  virtual absl::StatusOr<std::string> Encrypt(
      absl::string_view plaintext) const = 0;
  virtual absl::StatusOr<std::string> Decrypt(
      absl::string_view ciphertext) const = 0;
};

}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_CIPHER_H_
