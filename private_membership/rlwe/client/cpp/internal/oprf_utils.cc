#include "private_membership/rlwe/client/cpp/internal/oprf_utils.h"

#include <memory>
#include <optional>
#include <string>

#include "private_join_and_compute/crypto/ec_commutative_cipher.h"
#include "private_membership/rlwe/client/proto/private_membership_rlwe.pb.h"
#include "private_membership/rlwe/client/cpp/internal/constants.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"

namespace private_membership {
namespace rlwe {

absl::StatusOr<OprfCurveAndHashType> ToInternalEnum(OprfParams oprf_params) {
  switch (oprf_params) {
    case OprfParams::OPRF_CONFIG_NID_X9_62_PRIME256V1_WITH_SHA256:
      return OprfCurveAndHashType::kPrime256v1CurveIdWithSha256;
    case OprfParams::OPRF_CONFIG_NID_X9_62_PRIME256V1_WITH_SSWU_RO:
      return OprfCurveAndHashType::kPrime256v1CurveIdWithSswuRo;
    case OprfParams::OPRF_CONFIG_UNSPECIFIED:
      return OprfCurveAndHashType::kUndefinedOprfCurveAndHashType;
    default:
      return absl::InvalidArgumentError("Invalid OPRF params.");
  }
}

absl::StatusOr<std::unique_ptr<::private_join_and_compute::ECCommutativeCipher>>
CreateOprfCipher(OprfCurveAndHashType oprf_curve_and_hash_type,
                 std::optional<std::string> ec_cipher_key) {
  int curve_id;
  ::private_join_and_compute::ECCommutativeCipher::HashType hash_type;
  switch (oprf_curve_and_hash_type) {
    case OprfCurveAndHashType::kPrime256v1CurveIdWithSha256:
      curve_id = kPrime256v1CurveId;
      hash_type = ::private_join_and_compute::ECCommutativeCipher::HashType::SHA256;
      break;
    case OprfCurveAndHashType::kPrime256v1CurveIdWithSswuRo:
      curve_id = kPrime256v1CurveId;
      hash_type = ::private_join_and_compute::ECCommutativeCipher::HashType::SSWU_RO;
      break;
    case OprfCurveAndHashType::kUndefinedOprfCurveAndHashType:
    default:
      return absl::InvalidArgumentError(
          "Unsupported OprfCurveAndHashType specified.");
  }
  return ec_cipher_key.has_value()
             ? ::private_join_and_compute::ECCommutativeCipher::CreateFromKey(
                   curve_id, ec_cipher_key.value(), hash_type)
             : ::private_join_and_compute::ECCommutativeCipher::CreateWithNewKey(curve_id,
                                                                 hash_type);
}

}  // namespace rlwe
}  // namespace private_membership
