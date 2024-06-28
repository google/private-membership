//
// Created by Haris Mughees on 4/22/21.
//

#ifndef EXTERNAL_PROD_PIR_CLIENT_H
#define EXTERNAL_PROD_PIR_CLIENT_H

#include "external_prod.h"
#include "pir.h"

using namespace std;

class pir_client {
public:
    pir_client(const seal::EncryptionParameters &parms,
              const PirParams &pirparms);

    GaloisKeys generate_galois_keys();
    PirQuery generate_query(std::uint64_t desiredIndex);

    PirQuery generate_query_combined(std::uint64_t desiredIndex);
    PirQuery generate_query_combined2(const vector<int>& first_dim_indices, uint64_t index);
    PirQuery generate_query_combined3(const vector<vector<int>>& indices);

    PirQuery generate_perm_query(std::uint64_t desiredIndex);
    void compute_inverse_scales();
    // Index and offset of an element in an FV plaintext
    uint64_t get_fv_index(uint64_t element_idx, uint64_t ele_size);
    uint64_t get_fv_offset(uint64_t element_idx, uint64_t ele_size);
    void decrypt_results(std::vector<seal::Ciphertext> reply);

    Plaintext decrypt_result(std::vector<seal::Ciphertext> reply);

    GSWCiphertext get_enc_sk();

    void test_query_expansion(PirQuery query,GaloisKeys galkey);

    SecretKey get_decryptor();

private:
    seal::EncryptionParameters params_;
    PirParams pir_params_;

    std::unique_ptr<seal::Encryptor> encryptor_;
    std::unique_ptr<seal::Decryptor> decryptor_;
    std::unique_ptr<seal::Evaluator> evaluator_;
    std::unique_ptr<seal::KeyGenerator> keygen_;
    std::shared_ptr<seal::SEALContext> newcontext_;


    vector<uint64_t> indices_; // the indices for retrieval.

    vector<int> first_dim_indices_;
    vector<uint64_t> indices2_;

    vector<uint64_t> inverse_scales_;

    seal::Ciphertext compose_to_ciphertext(std::vector<seal::Plaintext> plains);


    friend class PIRServer;
};



#endif //EXTERNAL_PROD_PIR_CLIENT_H
