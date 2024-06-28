//
// Created by Haris Mughees on 4/14/21.
//

#pragma once

#ifndef EXTERNAL_PROD_EXTERNAL_PROD_H
#define EXTERNAL_PROD_EXTERNAL_PROD_H

 //EXTERNAL_PROD_EXTERNAL_PROD_H


#include <iostream>
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/polycore.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include <pthread.h>
#include "nfl.hpp"
#include "tools.h"
#include "seal/seal.h"
#include "waksman.h"

using namespace std;
using namespace std::chrono;
using namespace std;
using namespace seal;
using namespace seal::util;
using poly_t = nfl::poly_from_modulus<uint64_t, 4096, 128>;

#endif

typedef vector<Ciphertext> GSWCiphertext;


void multiply_add_plain_no_scaling_variant(
        const Plaintext &plain,
        const SEALContext::ContextData &context_data,
        const int shift_amount,
        uint64_t *destination, seal::util::MemoryPool &pool, uint64_t inv

);

void poly_nfllib_prod_with_no_ntt(std::uint64_t *p1, std::uint64_t *p2, std::uint64_t *res, const size_t coeff_count,
const std::uint64_t coeff_mod_count);

void poly_nfllib_add(std::uint64_t *p1, std::uint64_t *p2, std::uint64_t *res);

void poc_nfllib_ntt_rlwe_decomp(vector<uint64_t *> &rlwe_expansion);

void poc_nfllib_ntt_gsw(vector<Ciphertext> &gsw_enc, shared_ptr<SEALContext> &context);

void poc_nfllib_ntt_ct(Ciphertext &ct , shared_ptr<SEALContext> &context);
void poc_nfllib_intt_ct(Ciphertext &ct , shared_ptr<SEALContext> &context);

void poc_nfllib_add_ct(Ciphertext &ct1 ,Ciphertext &ct2 , shared_ptr<SEALContext> &context);

void poc_nfllib_plain_ct_prod(Ciphertext &ct , Plaintext &pt,
                              shared_ptr<SEALContext> &context, Ciphertext &res_ct);

void
plain_decompositions(Plaintext &pt, shared_ptr<SEALContext> &context, const uint64_t decomp_size, const uint64_t base_bit,
                     vector<uint64_t *> &plain_decom);

void poc_decomp_plain(Plaintext pt, const uint64_t decomp_size, shared_ptr<SEALContext> context,
                      vector<uint64_t *> &vec_ciphertexts, int base_bit, seal::util::MemoryPool &pool);

void poc_plain_gsw_enc128(const uint64_t decomp_size, const uint64_t base_bit, shared_ptr<SEALContext> context,
                          const SecretKey sk, vector<Ciphertext> &gsw_ciphertext, Plaintext gsw_plain,
                          seal::util::MemoryPool &pool, uint64_t inv);

void poc_plain_gsw_enc128_combined(const uint64_t decomp_size, const uint64_t base_bit, shared_ptr<SEALContext> context,
                          const SecretKey sk, vector<Ciphertext> &gsw_ciphertext, Plaintext gsw_plain,
                          seal::util::MemoryPool &pool, uint64_t inv, uint64_t gap);

vector<Ciphertext> poc_rlwe_expand(Ciphertext packedquery, shared_ptr<SEALContext> context, seal::GaloisKeys galkey, uint64_t size);

void poc_rlwe_expand_threaded(Ciphertext packedquery, shared_ptr<SEALContext> context, seal::GaloisKeys galkey, uint64_t size, vector<Ciphertext> &result);

void poc_decompose_array(uint64_t *value, size_t count, std::vector<Modulus> coeff_modulus, size_t coeff_mod_count,
                         MemoryPoolHandle pool);

void poc_expand_flat(vector<GSWCiphertext>::iterator &result, vector<Ciphertext> &packed_swap_bits,
                     shared_ptr<SEALContext> context, int size, seal::GaloisKeys &galkey);


void poc_expand_flat_combined(vector<GSWCiphertext>::iterator &result, vector<Ciphertext> &packed_swap_bits,
                     shared_ptr<SEALContext> context, int size, seal::GaloisKeys &galkey);



void poc_expand_flat_threaded(vector<GSWCiphertext>::iterator &result, vector<Ciphertext> &packed_swap_bits,
                     shared_ptr<SEALContext> context, int size, seal::GaloisKeys &galkey);

void rwle_decompositions(Ciphertext rlwe_ct_1, shared_ptr<SEALContext> context, const uint64_t l, const uint64_t base_bit,
                         vector<uint64_t *> &rlwe_decom);

void multiply_power_of_X(const Ciphertext &encrypted, Ciphertext &destination,
                         uint32_t index, shared_ptr<SEALContext> context);

GaloisKeys generate_galois_keys(shared_ptr<SEALContext> context, KeyGenerator &keygen);
/*!
 * This function sets parameters for BFV encryption.
 * @param parms
 */
void set_bfv_parms(EncryptionParameters &parms);

/*!
 * GSW encryption coeff modulus of 128 bits
 * @param l - number of rows in gsw ciphertexts
 * @param base_bit - size of base for decomposition
 * @param context - BFV public context information
 * @param sk- Secret key of the user
 * @param gsw_ciphertext - Output ciphertext
 * @param gsw_plain - input plain text to be encrypted
 * @param pool - memory pool
 * @param inv - this encrypts msg* 1/ inv
 */
void poc_gsw_enc128(const uint64_t l, const uint64_t base_bit, shared_ptr<SEALContext> context,
                    const SecretKey sk, vector<Ciphertext> &gsw_ciphertext, Plaintext gsw_plain,
                    seal::util::MemoryPool &pool, uint64_t inv);

void poc_enc_sk_gsw(SecretKey sk, shared_ptr<SEALContext> context,  const int base_bit, vector<Ciphertext>& sk_gsw_ciphertext);

void poc_half_gsw_enc128(const uint64_t l, const uint64_t base_bit, shared_ptr<SEALContext> context,
                         const SecretKey sk, vector<Ciphertext> &gsw_ciphertext, Plaintext gsw_plain,
                         seal::util::MemoryPool &pool, uint64_t inv);

void poc_half_gsw_enc128_combined(const uint64_t l, const uint64_t base_bit, shared_ptr<SEALContext> context,
                         const SecretKey sk, vector<Ciphertext> &gsw_ciphertext, Plaintext gsw_plain,
                         seal::util::MemoryPool &pool, uint64_t total_expand_size, uint64_t gap, uint64_t dimension_size);

void mymultiply_add_plain_without_scaling_variant(
        const Plaintext &plain,
        const SEALContext::ContextData &context_data,
        const int shift_amount,
        uint64_t *destination, seal::util::MemoryPool &pool, uint64_t inv = 0

);

void mymultiply_add_plain_without_scaling_variant_combined(
        const Plaintext &plain,
        const SEALContext::ContextData &context_data,
        const int shift_amount,
        uint64_t *destination, seal::util::MemoryPool &pool, uint64_t total_expand_size, uint64_t gap, uint64_t dimension_size
        , uint64_t curr_component
);

void poc_gsw_enc128_sk(const uint64_t l, const uint64_t base_bit, shared_ptr<SEALContext> context,
                       const SecretKey sk, vector<Ciphertext> &gsw_ciphertext, Plaintext gsw_plain,
                       seal::util::MemoryPool &pool);

void thread_server_expand(vector<GSWCiphertext>::iterator &result, vector<Ciphertext> packed_swap_bits,
                          shared_ptr<SEALContext> context, int begin, int end, int size, seal::GaloisKeys galkey,
                          const int l, const int base_bit, const int lsk, const int bsk,vector<Ciphertext> sk_gsw_ciphertext2);

void gsw_server_expand_combined(vector<GSWCiphertext>::iterator &result, vector<Ciphertext> packed_swap_bits,
                                                                shared_ptr<SEALContext> context, int begin, int end, int total_dim_size, seal::GaloisKeys galkey,
                                                                const int l, const int base_bit,const int lsk, const int bsk, vector<Ciphertext> sk_gsw_ciphertext2,
                                                                uint64_t total_expand_size);
void mymultiply_add_plain_without_scaling_variant_sk(
        const Plaintext &plain,
        const SEALContext::ContextData &context_data,
        const int shift_amount,
        uint64_t *destination,
        seal::util::MemoryPool &pool
);


void my_poc_gsw_enc128_sk(const uint64_t l, const uint64_t base_bit, shared_ptr<SEALContext> context,
                          const SecretKey sk, vector<Ciphertext> &gsw_ciphertext, Plaintext gsw_plain,
                          seal::util::MemoryPool &pool);
void my_decompose_array(uint64_t *value, size_t count, std::vector<Modulus> coeff_modulus, size_t coeff_mod_count,
                        MemoryPoolHandle pool);

void my_poc_decomp_rlwe128(Ciphertext ct, const uint64_t l, shared_ptr<SEALContext> context,
                        vector<uint64_t *> &vec_ciphertexts, int base_bit, seal::util::MemoryPool &pool);
void
my_rwle_decompositions(Ciphertext rlwe_ct_1, shared_ptr<SEALContext> context, const uint64_t l, const uint64_t base_bit,
                    vector<uint64_t *> &rlwe_decom);

void set_ciphertext(Ciphertext &ct, shared_ptr<SEALContext> context);
vector<Ciphertext>
rlweExpand(Ciphertext packedquery, shared_ptr<SEALContext> context, seal::GaloisKeys galkey, uint64_t size);

void my_poc_nfllib_external_product(vector<Ciphertext> gsw_enc, vector<uint64_t *> rlwe_expansion,
                                 shared_ptr<SEALContext> context,
                                 int l, Ciphertext &res_ct, int is_reusable=1);

/*!
 *  A helper function used in gsw encryption. As GSW is not a scale invariant scheme therefore plain text is encrypted without multiplying to (delta)
 * @param plain - input plaintext
 * @param context_data - public context data
 * @param shift_amount - gadget matrix entry
 * @param destination - output (polynomial pointer)
 * @param pool
 * @param inv
 */
void poc_multiply_add_plain_without_scaling_variant(
        const Plaintext &plain,
        const SEALContext::ContextData &context_data,
        const int shift_amount,
        uint64_t *destination, seal::util::MemoryPool &pool, uint64_t inv

);

void poc_multiply_add_plain_without_scaling_variant_sk(
        const Plaintext &plain,
        const SEALContext::ContextData &context_data,
        const int shift_amount,
        uint64_t *destination,
        seal::util::MemoryPool &pool
);



/*!
 * main function to decompose bfv ciphertext into 2l parts, parts are created in size of base_bit
 * @param ct - input bfv ciphertext
 * @param l - number of rows in GSW
 * @param context -
 * @param vec_ciphertexts -output decomposed
 * @param base_bit - size of each part
 * @param pool
 */
void poc_decomp_rlwe128(Ciphertext ct, const uint64_t l, shared_ptr<SEALContext> context,
                        vector<uint64_t *> &vec_ciphertexts, int base_bit, seal::util::MemoryPool &pool);



/*!
 * Main function to perform external product between GSW ciphertext and decomposed vector of BFV ciphertext
 * @param gsw_enc -GSW ciphertext, due to noise growth this should only encrypt {0,1}
 * @param rlwe_expansion - decomposed vector of BFV ciphertext
 * @param context
 * @param l
 * @param res_ct - output ciphertext
 * @param is_reusable - reuse nTT conversion
 */
void poc_nfllib_external_product(vector<Ciphertext> &gsw_enc, vector<uint64_t *> &rlwe_expansion,
                                 shared_ptr<SEALContext> &context,
                                 int l, Ciphertext &res_ct, int is_reusable);


/*!
 * Fast function based on NFLlib for polynomial multiplications
 * @param p1
 * @param p2
 * @param res
 * @param coeff_count
 * @param coeff_mod_count
 * @param is_reusable
 */
void poly_nfllib_mul(std::uint64_t *p1, std::uint64_t *p2, std::uint64_t *res, const size_t coeff_count,
                     const std::uint64_t coeff_mod_count, int is_reusable);