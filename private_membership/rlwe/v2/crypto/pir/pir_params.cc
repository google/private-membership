#include "private_membership/rlwe/v2/crypto/pir/pir_params.h"

#include <cstdint>

#include "private_membership/rlwe/v2/crypto/encryption.h"
#include "absl/numeric/bits.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/str_cat.h"

namespace private_membership {
namespace rlwe {
namespace v2 {

template <typename CoeffType>
absl::StatusOr<PirParams<CoeffType>> PirParams<CoeffType>::Create(
    const RlweParams<CoeffType>& rlwe_params, const CoeffType plaintext_modulus,
    const GadgetParams& pack_gadget_params,
    const GadgetParams& poly_eval_gadget_params, const int num_entries,
    const int interpolation_degree, const int entry_size_multiple) {
  if (num_entries <= 0) {
    return absl::InvalidArgumentError("num_entries must be strictly positive.");
  }

  if (entry_size_multiple <= 0) {
    return absl::InvalidArgumentError(
        "entry_size_multiple must be strictly positive.");
  }

  if (interpolation_degree <= 0 ||
      absl::has_single_bit(static_cast<unsigned int>(interpolation_degree)) !=
          1) {
    return absl::InvalidArgumentError(
        "interpolation_degree must be a power of two.");
  }

  if (num_entries % interpolation_degree != 0) {
    return absl::InvalidArgumentError(
        absl::StrCat("interpolation_degree (", interpolation_degree,
                     ") must divide num_entries (", num_entries, ")."));
  }

  const int rows = num_entries / interpolation_degree;
  if (rows % rlwe_params.Degree() != 0) {
    return absl::InvalidArgumentError(
        absl::StrCat("(num_entries / interpolation_degree) (", rows,
                     ") must be divisible by the ring degree (",
                     rlwe_params.Degree(), ")."));
  }

  if (interpolation_degree > 2 * rlwe_params.Degree()) {
    return absl::InvalidArgumentError(
        absl::StrCat("interpolation_degree (", interpolation_degree,
                     ") must be less than or equal to 2 * ring degree (",
                     2 * rlwe_params.Degree(), ")."));
  }

  return PirParams<CoeffType>(rlwe_params, plaintext_modulus,
                              pack_gadget_params, poly_eval_gadget_params,
                              num_entries, interpolation_degree,
                              entry_size_multiple);
}

template class PirParams<uint32_t>;
template class PirParams<uint64_t>;

}  // namespace v2
}  // namespace rlwe
}  // namespace private_membership
