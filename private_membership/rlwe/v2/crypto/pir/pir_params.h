#ifndef PRIVACY_PRIVATE_MEMBERSHIP_RLWE_V2_CRYPTO_PIR_PIR_PARAMS_H_
#define PRIVACY_PRIVATE_MEMBERSHIP_RLWE_V2_CRYPTO_PIR_PIR_PARAMS_H_

#include <cstdint>
#include <vector>

#include "private_membership/rlwe/v2/crypto/encryption.h"
#include "private_membership/rlwe/v2/proto/pir.pb.h"
#include "absl/numeric/bits.h"
#include "shell_encryption/status_macros.h"
#include "absl/status/statusor.h"

namespace private_membership {
namespace rlwe {
namespace v2 {

// Parameters for PIR.
//
// We treat each database entry as Z_p^d where d is the ring degree and p is the
// plaintext modulus.
template <typename CoeffType>
class PirParams {
 public:
  // Creates and validates a new PirParams instance.
  // entry_size_multiple is the number of RLWE plaintexts needed to encode each
  // database entry.
  // Returns InvalidArgumentError if:
  // - `num_entries` <= 0
  // - `interpolation_degree` is not a power of 2
  // - `interpolation_degree` does not divide num_entries
  // - `num_entries` / `interpolation_degree` is not divisible by
  //     rlwe_params.Degree()
  static absl::StatusOr<PirParams<CoeffType>> Create(
      const RlweParams<CoeffType>& rlwe_params, CoeffType plaintext_modulus,
      const GadgetParams& pack_gadget_params,
      const GadgetParams& poly_eval_gadget_params, int num_entries,
      int interpolation_degree, int entry_size_multiple = 1);

  // Accessors
  const RlweParams<CoeffType>& RlweParameters() const { return rlwe_params_; }

  CoeffType PlaintextModulus() const { return plaintext_modulus_; }

  const GadgetParams& PackGadgetParams() const { return pack_gadget_params_; }

  const GadgetParams& PolyEvalGadgetParams() const {
    return poly_eval_gadget_params_;
  }

  int NumEntries() const { return num_entries_; }

  int InterpolationDegree() const { return interpolation_degree_; }

  int EntrySizeMultiple() const { return entry_size_multiple_; }

 private:
  PirParams(const RlweParams<CoeffType>& rlwe_params,
            CoeffType plaintext_modulus, const GadgetParams& pack_gadget_params,
            const GadgetParams& poly_eval_gadget_params, int num_entries,
            int interpolation_degree, int entry_size_multiple)
      : rlwe_params_(rlwe_params),
        plaintext_modulus_(plaintext_modulus),
        pack_gadget_params_(pack_gadget_params),
        poly_eval_gadget_params_(poly_eval_gadget_params),
        num_entries_(num_entries),
        interpolation_degree_(interpolation_degree),
        entry_size_multiple_(entry_size_multiple) {}

  const RlweParams<CoeffType> rlwe_params_;
  const CoeffType plaintext_modulus_;
  const GadgetParams pack_gadget_params_;
  const GadgetParams poly_eval_gadget_params_;
  const int num_entries_;
  const int interpolation_degree_;
  const int entry_size_multiple_;
};

extern template class PirParams<uint32_t>;
extern template class PirParams<uint64_t>;

template <typename CoeffType>
absl::StatusOr<proto::PirParams> SerializePirParams(
    const PirParams<CoeffType>& params,
    proto::PirParams::DbValueType db_value_type,
    proto::PirParams::MatCoeffType mat_coeff_type) {
  proto::PirParams proto;

  // Serialize RlweParams
  auto* rlwe_proto = proto.mutable_rlwe_params();
  const auto& rlwe_params = params.RlweParameters();
  rlwe_proto->set_log_degree(
      absl::bit_width(static_cast<uint32_t>(rlwe_params.Degree() - 1)));
  rlwe_proto->set_log_modulus(rlwe_params.LogModulus());
  rlwe_proto->set_log_modulus1_after_switch(
      rlwe_params.LogModulus1AfterSwitch());
  rlwe_proto->set_log_modulus2_after_switch(
      rlwe_params.LogModulus2AfterSwitch());
  rlwe_proto->set_variance(rlwe_params.Variance());

  proto.set_plaintext_modulus(static_cast<int64_t>(params.PlaintextModulus()));

  // Serialize GadgetParams
  auto* pack_proto = proto.mutable_pack_gadget_params();
  pack_proto->set_log_digit(params.PackGadgetParams().log_digit);
  pack_proto->set_num_digits(params.PackGadgetParams().num_digits);

  auto* eval_proto = proto.mutable_poly_eval_gadget_params();
  eval_proto->set_log_digit(params.PolyEvalGadgetParams().log_digit);
  eval_proto->set_num_digits(params.PolyEvalGadgetParams().num_digits);

  proto.set_num_entries(params.NumEntries());
  proto.set_interpolation_degree(params.InterpolationDegree());
  proto.set_entry_size_multiple(params.EntrySizeMultiple());

  // Set type parameters
  proto.set_db_value_type(db_value_type);
  if constexpr (std::is_same_v<CoeffType, uint32_t>) {
    proto.set_coeff_type(proto::PirParams::COEFF_TYPE_UINT32);
  } else if constexpr (std::is_same_v<CoeffType, uint64_t>) {
    proto.set_coeff_type(proto::PirParams::COEFF_TYPE_UINT64);
  }
  proto.set_mat_coeff_type(mat_coeff_type);

  return proto;
}

template <typename CoeffType>
absl::StatusOr<PirParams<CoeffType>> DeserializePirParams(
    const proto::PirParams& proto) {
  // Reconstruct RlweParams
  const auto& rlwe_proto = proto.rlwe_params();
  int degree = 1 << rlwe_proto.log_degree();
  CoeffType modulus = rlwe_proto.log_modulus() == 8 * sizeof(CoeffType)
                          ? 0
                          : (CoeffType{1} << rlwe_proto.log_modulus());
  CoeffType modulus1 =
      rlwe_proto.log_modulus1_after_switch() == 8 * sizeof(CoeffType)
          ? 0
          : (CoeffType{1} << rlwe_proto.log_modulus1_after_switch());
  CoeffType modulus2 =
      rlwe_proto.log_modulus2_after_switch() == 8 * sizeof(CoeffType)
          ? 0
          : (CoeffType{1} << rlwe_proto.log_modulus2_after_switch());

  RLWE_ASSIGN_OR_RETURN(auto rlwe_params, RlweParams<CoeffType>::Create(
                                         degree, modulus, modulus1, modulus2,
                                         rlwe_proto.variance()));

  GadgetParams pack_gadget_params{
      .log_digit = proto.pack_gadget_params().log_digit(),
      .num_digits = proto.pack_gadget_params().num_digits(),
  };

  GadgetParams poly_eval_gadget_params{
      .log_digit = proto.poly_eval_gadget_params().log_digit(),
      .num_digits = proto.poly_eval_gadget_params().num_digits(),
  };

  return PirParams<CoeffType>::Create(
      rlwe_params, static_cast<CoeffType>(proto.plaintext_modulus()),
      pack_gadget_params, poly_eval_gadget_params, proto.num_entries(),
      proto.interpolation_degree(), proto.entry_size_multiple());
}

// A request object generated by the client, consisting of "b" components of the
// RLWE ciphertexts.
template <typename CoeffType>
struct PirRequest {
  // The part of PIR query corresponding to the first dimension of database
  // items.
  std::vector<Polynomial<CoeffType>> first_dimension_query;

  // The part of PIR query corresponding to the second dimension of database
  // items.
  std::vector<Polynomial<CoeffType>> second_dimension_query;

  // The key used for packing..
  std::vector<Polynomial<CoeffType>> packing_key;
};

// A response object containing the result of the PIR query.
template <typename CoeffType>
struct PirResponse {
  // The ciphertexts containing the queried item shards, encrypted under
  // client's secret key.
  std::vector<RlweCiphertext<CoeffType>> ciphertexts;
};

}  // namespace v2
}  // namespace rlwe
}  // namespace private_membership

#endif  // PRIVACY_PRIVATE_MEMBERSHIP_RLWE_V2_CRYPTO_PIR_PIR_PARAMS_H_
