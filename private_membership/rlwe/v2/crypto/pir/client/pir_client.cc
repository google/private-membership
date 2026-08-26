#include "private_membership/rlwe/v2/crypto/pir/client/pir_client.h"

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "private_membership/rlwe/v2/crypto/encryption.h"
#include "private_membership/rlwe/v2/crypto/pir/pir_params.h"
#include "private_membership/rlwe/v2/crypto/polynomial.h"
#include "absl/memory/memory.h"
#include "absl/status/status.h"
#include "shell_encryption/status_macros.h"
#include "absl/status/statusor.h"
#include "shell_encryption/prng/prng.h"

namespace private_membership {
namespace rlwe {
namespace v2 {

template <typename CoeffType>
absl::StatusOr<std::unique_ptr<PirClient<CoeffType>>>
PirClient<CoeffType>::Create(
    const PirParams<CoeffType>& params,
    std::unique_ptr<::rlwe::SecurePrng> a_component_prng,
    std::unique_ptr<::rlwe::SecurePrng> prng) {
  if (a_component_prng == nullptr) {
    return absl::InvalidArgumentError("a_component_prng must not be null.");
  }
  if (prng == nullptr) {
    return absl::InvalidArgumentError("prng must not be null.");
  }
  RLWE_ASSIGN_OR_RETURN(auto secret_key,
                   SampleSecretKey(params.RlweParameters(), *prng));
  return absl::WrapUnique(
      new PirClient<CoeffType>(params, std::move(a_component_prng),
                               std::move(prng), std::move(secret_key)));
}

template <typename CoeffType>
absl::StatusOr<PirRequest<CoeffType>> PirClient<CoeffType>::CreateRequest(
    const int index) {
  if (index < 0 || index >= params_.NumEntries()) {
    return absl::InvalidArgumentError("Index out of bounds.");
  }

  PirRequest<CoeffType> request;
  RLWE_ASSIGN_OR_RETURN(request.first_dimension_query,
                   CreateFirstDimensionQuery(index));
  if (params_.InterpolationDegree() > 1) {
    RLWE_ASSIGN_OR_RETURN(request.second_dimension_query,
                     CreateSecondDimensionQuery(index));
  } else {
    request.second_dimension_query = {};
  }
  RLWE_ASSIGN_OR_RETURN(request.packing_key, CreatePackingKey());

  return request;
}

template <typename CoeffType>
absl::StatusOr<PirRequest<CoeffType>>
PirClient<CoeffType>::CreateRequestForVector(
    const std::vector<CoeffType>& vec) {
  if (params_.InterpolationDegree() != 1) {
    return absl::InvalidArgumentError(
        "CreateRequestForVector is only supported for interpolation degree 1.");
  }
  PirRequest<CoeffType> request;
  RLWE_ASSIGN_OR_RETURN(request.first_dimension_query,
                   CreateFirstDimensionQueryForVector(vec));
  request.second_dimension_query = {};
  RLWE_ASSIGN_OR_RETURN(request.packing_key, CreatePackingKey());
  return request;
}

template <typename CoeffType>
absl::StatusOr<std::vector<Polynomial<CoeffType>>>
PirClient<CoeffType>::CreateFirstDimensionQueryForVector(
    const std::vector<CoeffType>& vec) const {
  const auto& rlwe_params = params_.RlweParameters();
  const int d = rlwe_params.Degree();
  if (vec.empty() || vec.size() % d != 0) {
    return absl::InvalidArgumentError(
        "Vector size must be a non-zero multiple of ring degree d.");
  }
  const int num_samples = vec.size() / d;
  RLWE_ASSIGN_OR_RETURN(auto a_first, SampleAComponents(rlwe_params, num_samples,
                                                   *a_component_prng_));
  RLWE_ASSIGN_OR_RETURN(
      auto rlwe_samples,
      GenerateRlweSamples(rlwe_params, secret_key_, a_first, *prng_));

  std::vector<Polynomial<CoeffType>> query;
  query.reserve(num_samples);
  for (int idx = 0; idx < num_samples; ++idx) {
    std::vector<CoeffType> chunk(vec.begin() + idx * d,
                                 vec.begin() + (idx + 1) * d);
    RLWE_ASSIGN_OR_RETURN(
        auto ct, EncryptFromRlweSample(rlwe_params, params_.PlaintextModulus(),
                                       rlwe_samples[idx], chunk));
    query.push_back(std::move(ct.b));
  }
  return query;
}

template <typename CoeffType>
absl::StatusOr<std::vector<Polynomial<CoeffType>>>
PirClient<CoeffType>::CreateFirstDimensionQuery(const int index) const {
  const auto& rlwe_params = params_.RlweParameters();
  const int d = rlwe_params.Degree();
  const int t = params_.InterpolationDegree();

  const int entries_per_row = t;
  const int row_index = index / entries_per_row;
  const int i = row_index / d;
  const int j = row_index % d;

  const int num_samples = params_.NumEntries() / t / d;
  RLWE_ASSIGN_OR_RETURN(auto a_first, SampleAComponents(rlwe_params, num_samples,
                                                   *a_component_prng_));
  RLWE_ASSIGN_OR_RETURN(
      auto rlwe_samples,
      GenerateRlweSamples(rlwe_params, secret_key_, a_first, *prng_));

  std::vector<Polynomial<CoeffType>> query;
  query.reserve(rlwe_samples.size());

  for (int idx = 0; idx < rlwe_samples.size(); ++idx) {
    std::vector<CoeffType> message(d, 0);
    if (idx == i) {
      message[j] = 1;
    }
    RLWE_ASSIGN_OR_RETURN(
        auto ct, EncryptFromRlweSample(rlwe_params, params_.PlaintextModulus(),
                                       rlwe_samples[idx], message));
    query.push_back(std::move(ct.b));
  }

  return query;
}

template <typename CoeffType>
absl::StatusOr<std::vector<Polynomial<CoeffType>>>
PirClient<CoeffType>::CreateSecondDimensionQuery(const int index) const {
  const auto& rlwe_params = params_.RlweParameters();
  const int d = rlwe_params.Degree();
  const int t = params_.InterpolationDegree();

  const int entries_per_row = t;
  const int col_index = index % entries_per_row;

  const int num_samples = 2 * params_.PolyEvalGadgetParams().num_digits;

  RLWE_ASSIGN_OR_RETURN(
      std::vector<Polynomial<CoeffType>> a_second,
      SampleAComponents(rlwe_params, num_samples, *a_component_prng_));
  RLWE_ASSIGN_OR_RETURN(
      std::vector<RlweSample<CoeffType>> rlwe_samples,
      GenerateRlweSamples(rlwe_params, secret_key_, a_second, *prng_));

  std::vector<CoeffType> message(d, 0);
  const int second_dim_shift = (2 * d / t) * col_index;
  if (second_dim_shift >= d) {
    message[second_dim_shift - d] = rlwe_params.Modulus() - 1;
  } else {
    message[second_dim_shift] = 1;
  }

  const int half = rlwe_samples.size() / 2;
  std::vector<RlweSample<CoeffType>> u_samples(rlwe_samples.begin(),
                                               rlwe_samples.begin() + half);
  std::vector<RlweSample<CoeffType>> v_samples(rlwe_samples.begin() + half,
                                               rlwe_samples.end());

  RLWE_ASSIGN_OR_RETURN(auto rgsw_ct,
                   EncryptRgswCiphertextFromRlweSamples(
                       rlwe_params, params_.PolyEvalGadgetParams(), u_samples,
                       v_samples, secret_key_, message));

  std::vector<Polynomial<CoeffType>> query;
  query.reserve(rlwe_samples.size());
  for (auto& ct : rgsw_ct.u) {
    query.push_back(std::move(ct.b));
  }
  for (auto& ct : rgsw_ct.v) {
    query.push_back(std::move(ct.b));
  }
  return query;
}

template <typename CoeffType>
absl::StatusOr<std::vector<Polynomial<CoeffType>>>
PirClient<CoeffType>::CreatePackingKey() const {
  const auto& rlwe_params = params_.RlweParameters();
  const int d = rlwe_params.Degree();

  const int num_samples = 2 * params_.PackGadgetParams().num_digits;
  RLWE_ASSIGN_OR_RETURN(
      std::vector<Polynomial<CoeffType>> a_pack,
      SampleAComponents(rlwe_params, num_samples, *a_component_prng_));
  RLWE_ASSIGN_OR_RETURN(
      std::vector<RlweSample<CoeffType>> rlwe_samples,
      GenerateRlweSamples(rlwe_params, secret_key_, a_pack, *prng_));

  RLWE_ASSIGN_OR_RETURN(Polynomial<CoeffType> sk_g, secret_key_.Automorph(5));
  RLWE_ASSIGN_OR_RETURN(Polynomial<CoeffType> sk_h,
                   secret_key_.Automorph(2 * d - 1));

  const int half = rlwe_samples.size() / 2;
  std::vector<RlweSample<CoeffType>> pack_samples_g(
      rlwe_samples.begin(), rlwe_samples.begin() + half);
  std::vector<RlweSample<CoeffType>> pack_samples_h(rlwe_samples.begin() + half,
                                                    rlwe_samples.end());

  RLWE_ASSIGN_OR_RETURN(std::vector<RlweCiphertext<CoeffType>> ksm_g,
                   EncryptGadgetCiphertextFromRlweSamples(
                       rlwe_params, params_.PackGadgetParams(), pack_samples_g,
                       sk_g.Coeffs()));
  RLWE_ASSIGN_OR_RETURN(std::vector<RlweCiphertext<CoeffType>> ksm_h,
                   EncryptGadgetCiphertextFromRlweSamples(
                       rlwe_params, params_.PackGadgetParams(), pack_samples_h,
                       sk_h.Coeffs()));

  std::vector<Polynomial<CoeffType>> packing_key;
  packing_key.reserve(rlwe_samples.size());
  for (auto& ct : ksm_g) {
    packing_key.push_back(std::move(ct.b));
  }
  for (auto& ct : ksm_h) {
    packing_key.push_back(std::move(ct.b));
  }
  return packing_key;
}

template <typename CoeffType>
absl::StatusOr<std::vector<CoeffType>> PirClient<CoeffType>::ProcessResponse(
    const PirResponse<CoeffType>& response) const {
  const int entry_size_multiple = params_.EntrySizeMultiple();
  if (response.ciphertexts.size() != entry_size_multiple) {
    return absl::InvalidArgumentError(
        "Response must contain exactly entry_size_multiple ciphertexts.");
  }
  std::vector<CoeffType> decrypted_message;
  decrypted_message.reserve(entry_size_multiple *
                            params_.RlweParameters().Degree());

  for (int k = 0; k < entry_size_multiple; ++k) {
    RLWE_ASSIGN_OR_RETURN(std::vector<CoeffType> shard_decrypted,
                     DecryptAfterModulusSwitch(
                         params_.RlweParameters(), params_.PlaintextModulus(),
                         response.ciphertexts[k], secret_key_));
    decrypted_message.insert(decrypted_message.end(), shard_decrypted.begin(),
                             shard_decrypted.end());
  }
  return decrypted_message;
}

template class PirClient<uint32_t>;
template class PirClient<uint64_t>;

}  // namespace v2
}  // namespace rlwe
}  // namespace private_membership
