#ifndef PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_OPRF_UTILS_H_
#define PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_OPRF_UTILS_H_

#include <memory>
#include <optional>
#include <string>

#include "private_join_and_compute/crypto/ec_commutative_cipher.h"
#include "private_membership/rlwe/client/proto/private_membership_rlwe.pb.h"
#include "absl/status/statusor.h"

namespace private_membership {
namespace rlwe {

enum class OprfCurveAndHashType {
  kUndefinedOprfCurveAndHashType = 0,
  kPrime256v1CurveIdWithSha256 = 1,
  kPrime256v1CurveIdWithSswuRo = 2,
};

absl::StatusOr<OprfCurveAndHashType> ToInternalEnum(OprfParams oprf_params);

absl::StatusOr<std::unique_ptr<::private_join_and_compute::ECCommutativeCipher>>
CreateOprfCipher(OprfCurveAndHashType oprf_curve_and_hash_type,
                 std::optional<std::string> ec_cipher_key = std::nullopt);

}  // namespace rlwe
}  // namespace private_membership

#endif  // PRIVATE_MEMBERSHIP_RLWE_CLIENT_CPP_INTERNAL_OPRF_UTILS_H_
