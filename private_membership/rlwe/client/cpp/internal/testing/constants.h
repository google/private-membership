#ifndef PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_TESTING_CONSTANTS_H_
#define PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_TESTING_CONSTANTS_H_

#include <openssl/obj_mac.h>

namespace private_membership {
namespace rlwe {

// Identifier of the elliptic curve used in tests.
constexpr int kTestCurveId = NID_X9_62_prime256v1;

}  // namespace rlwe
}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_TESTING_CONSTANTS_H_
