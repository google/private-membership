//
// Created by Haris Mughees on 4/21/21.
//

#include "external_prod.h"
#include "pir.h"
#include "pir_client.h"


#ifndef EXTERNAL_PROD_PIR_SERVER_H
#define EXTERNAL_PROD_PIR_SERVER_H

#endif //EXTERNAL_PROD_PIR_SERVER_H

typedef std::vector<seal::Plaintext> Database;
typedef std::vector<seal::Ciphertext> PirReply;



class pir_server {
public:
    pir_server(const seal::EncryptionParameters &params, const PirParams &pir_params);
    void set_galois_key(std::uint32_t client_id, seal::GaloisKeys galkey);
    void set_database(const std::unique_ptr<const std::uint8_t[]> &bytes, std::uint64_t ele_num, std::uint64_t ele_size);
    void set_database(std::unique_ptr<std::vector<seal::Plaintext>> &&db, bool use_padding = false);

    // NOTE: server takes over ownership of db and frees it when it exits.
    // Caller cannot free db
    PirReply generate_reply(PirQuery query, uint32_t client_id, SecretKey sk);
    PirReply generate_reply_combined(PirQuery query, uint32_t client_id, SecretKey sk);
    void set_enc_sk(GSWCiphertext sk_enc);

    void preprocess_database();
private:
    seal::EncryptionParameters params_; // SEAL parameters
    PirParams pir_params_;// PIR parameters
    parms_id_type parms_id_;
    std::unique_ptr<Database> db_;
    bool is_db_preprocessed_;
    std::shared_ptr<seal::SEALContext> newcontext_;
    seal::GaloisKeys galoisKeys_;
    std::unique_ptr<seal::Evaluator> evaluator_;
    vector<uint64_t *> plain_decom;
    vector<vector<uint64_t *>> split_db;
    GSWCiphertext sk_enc_;



};




