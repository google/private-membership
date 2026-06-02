#include <iostream>
#include "seal/util/polyarithsmallmod.h"
#include "seal/util/polycore.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <random>
#include <pthread.h>
#include "nfl.hpp"
#include "tools.h"
#include "seal/seal.h"
#include "external_prod.h"
#include "util.h"
#include "pir.h"
#include "pir_server.h"
#include "pir_client.h"

#include "../sparse_pir/linear_system.h"
#include "../sparse_pir/plaintext_value.h"
#include "../sparse_pir/utils.h"

using namespace std;
using namespace std::chrono;
using namespace std;
using namespace seal;
using namespace seal::util;

typedef vector<Ciphertext> GSWCiphertext;

void
test_external_prod_with_sk(Evaluator &evaluator1, Encryptor &encryptor1, Decryptor &decryptor1, KeyGenerator &keygen,
                           shared_ptr<SEALContext> context, SecretKey sk) {

    const auto &context_data2 = context->first_context_data();
    auto &parms = context_data2->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);
    const int base_bits = 31;


    Plaintext gsw_plain(to_string(1));
    Plaintext msg;
    msg.resize(coeff_count);
    for (int h = 0; h < coeff_count; h++) {
        if (h == 0) {

            msg.data()[h] = 1;
            //cout<< msg.data()[h ]<<endl;
        } else
            msg.data()[h] = 1;
    }

    vector<Ciphertext> sk_gsw_ciphertext;




    ///test ct of rlwe
    Plaintext test_rlwe_pt("1");

    Ciphertext test_rlwe_ct;
    encryptor1.encrypt_symmetric(test_rlwe_pt, test_rlwe_ct);

    cout << "-----------------------------------------------" << endl;
    cout << "Noise budget before external product=" << decryptor1.invariant_noise_budget(test_rlwe_ct) << endl;

    vector<uint64_t *> rlwe_decom;


    int duration = 0;
    int interations = 0;
    for (int i = base_bits; i > 1; i = ceil(i / 2)) {
        interations++;
        const int lvl = context_data2->total_coeff_modulus_bit_count() / i;

        sk_gsw_ciphertext.clear();
        poc_enc_sk_gsw(sk, context, i, sk_gsw_ciphertext);
        //poc_gsw_enc128(lvl, i, context, sk, choice_bit, msg, pool, 0);


        //poc_nfllib_ntt_gsw(choice_bit, context);

        auto gsw_enc_time_start = std::chrono::steady_clock::now();

        rwle_decompositions(test_rlwe_ct, context, lvl, i, rlwe_decom);


        poc_nfllib_ntt_rlwe_decomp(rlwe_decom);


        /// steps for external product. Both rlwe and gsw must be crt-decomposed
        Ciphertext res_ct;
        res_ct.resize(context, context->first_context_data()->parms_id(), 2);


        poc_nfllib_external_product(sk_gsw_ciphertext, rlwe_decom, context, lvl, res_ct, 1);
        poc_nfllib_intt_ct(res_ct, context);
        auto gsw_enc_time_end = std::chrono::steady_clock::now();

        for (auto p : rlwe_decom) {
            free(p);
        }
        rlwe_decom.clear();


        duration = duration_cast<std::chrono::microseconds>(gsw_enc_time_end - gsw_enc_time_start).count();

        cout << "---------------------------------" << endl;
        cout << "For Base bits=" << i << endl;
        cout << "Noise budget after external product= "
             << decryptor1.invariant_noise_budget(res_ct) << endl;

        cout << "External prod duration= " << duration << "us" << endl;
    }


}

void test_external_prod(Evaluator &evaluator1, Encryptor &encryptor1, Decryptor &decryptor1, KeyGenerator &keygen,
                        shared_ptr<SEALContext> context, SecretKey sk) {

    const auto &context_data2 = context->first_context_data();
    auto &parms = context_data2->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);
    const int base_bits = 31;


    Plaintext gsw_plain(to_string(1));
    Plaintext msg;
    msg.resize(coeff_count);
    for (int h = 0; h < coeff_count; h++) {
        if (h == 0) {

            msg.data()[h] = 1;
            //cout<< msg.data()[h ]<<endl;
        } else
            msg.data()[h] = 1;
    }


    GSWCiphertext choice_bit;


    ///test ct of rlwe
    Plaintext test_rlwe_pt("1");

    Ciphertext test_rlwe_ct;
    encryptor1.encrypt_symmetric(test_rlwe_pt, test_rlwe_ct);

    cout << "-----------------------------------------------" << endl;
    cout << "Noise budget before external product=" << decryptor1.invariant_noise_budget(test_rlwe_ct) << endl;

    vector<uint64_t *> rlwe_decom;


    int duration = 0;
    int interations = 0;
    for (int i = base_bits; i > 1; i = ceil(i / 2)) {
        interations++;
        const int lvl = context_data2->total_coeff_modulus_bit_count() / i;
        choice_bit.clear();


        poc_gsw_enc128(lvl, i, context, sk, choice_bit, msg, pool, 0);


        poc_nfllib_ntt_gsw(choice_bit, context);

        auto gsw_enc_time_start = std::chrono::steady_clock::now();

        rwle_decompositions(test_rlwe_ct, context, lvl, i, rlwe_decom);


        poc_nfllib_ntt_rlwe_decomp(rlwe_decom);


        /// steps for external product. Both rlwe and gsw must be crt-decomposed
        Ciphertext res_ct;
        res_ct.resize(context, context->first_context_data()->parms_id(), 2);


        poc_nfllib_external_product(choice_bit, rlwe_decom, context, lvl, res_ct, 1);
        poc_nfllib_intt_ct(res_ct, context);
        auto gsw_enc_time_end = std::chrono::steady_clock::now();

        for (auto p : rlwe_decom) {
            free(p);
        }
        rlwe_decom.clear();


        duration = duration_cast<std::chrono::microseconds>(gsw_enc_time_end - gsw_enc_time_start).count();

        cout << "---------------------------------" << endl;
        cout << "For Base bits=" << i << endl;
        cout << "Noise budget after external product= "
             << decryptor1.invariant_noise_budget(res_ct) << endl;

        cout << "External prod duration= " << duration << "us" << endl;
    }


}

void test_nfllib_ct_add(Evaluator &evaluator1, Encryptor &encryptor1, Decryptor &decryptor1, KeyGenerator &keygen,
                        shared_ptr<SEALContext> context, SecretKey sk) {

    const auto &context_data2 = context->first_context_data();
    auto &parms = context_data2->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);
    Plaintext msg;
    msg.resize(coeff_count);
    for (int h = 0; h < coeff_count; h++) {
        if (h == 0) {

            msg.data()[h] = 1;
//cout<< msg.data()[h ]<<endl;
        } else
            msg.data()[h] = 0;
    }
    Ciphertext ct1, ct2;

    encryptor1.encrypt_symmetric(msg, ct1);
    encryptor1.encrypt_symmetric(msg, ct2);

    auto gsw_enc_time_start = std::chrono::high_resolution_clock::now();
    //poc_nfllib_ntt_ct(test_rlwe_ct, context);

    //poc_nfllib_plain_ct_prod(test_rlwe_ct , msg, context, res_ct);

    //poc_nfllib_add_ct(ct1,ct2,context);
    evaluator1.add_inplace(ct1, ct2);

    auto gsw_enc_time_end = std::chrono::high_resolution_clock::now();
    //poc_nfllib_intt_ct(test_rlwe_ct, context);

    int duration = duration_cast<std::chrono::microseconds>(gsw_enc_time_end - gsw_enc_time_start).count();
    Plaintext ppt;
    decryptor1.decrypt(ct1, ppt);
    cout << duration << endl;


}

void test_gsw_expansion(Evaluator &evaluator1, Encryptor &encryptor1, Decryptor &decryptor1, KeyGenerator &keygen,
                        shared_ptr<SEALContext> context, SecretKey sk) {

    GaloisKeys galois_keys = generate_galois_keys(context, keygen);
    const auto &context_data2 = context->first_context_data();
    auto &parms = context_data2->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);
    uint64_t base_bit = 16;
    const int l = 5;
    context_data2->total_coeff_modulus_bit_count() / base_bit;

    //int64_t swapbitsSize = get_swapbits_size(32);
    uint64_t swapbitsSize = 4096;
    int logsize = ceil(log2(swapbitsSize));
    int gap = ceil(coeff_count / (1 << logsize));


    cout << "swapbits size = " << swapbitsSize << endl;

    Plaintext msg;
    msg.resize(coeff_count);
    for (int h = 0; h < swapbitsSize; h++) {
        if (h == 0)
            msg.data()[h * gap] = ((int64_t) 1);
    }

    // get upper half (l) part of gsw where b*B^i is added to C_0
    vector<Ciphertext> half_gsw_ciphertext;
    //poc_l_pack_enc128(l, base_bit, context, sk, half_gsw_ciphertext, msg, decryptor1,  pool);


    poc_half_gsw_enc128(l, base_bit, context, sk, half_gsw_ciphertext, msg, pool, (1 << logsize));


    int bsk = 16;
    int lsk = context_data2->total_coeff_modulus_bit_count() / base_bit;;
    vector<Ciphertext> sk_gsw_ciphertext;
    poc_enc_sk_gsw(sk, context, bsk, sk_gsw_ciphertext);


    vector<GSWCiphertext> CtMuxBits;
    CtMuxBits.resize((1 << logsize), GSWCiphertext(2 * l));


    int size = (1 << logsize);
    vector<GSWCiphertext>::iterator gswCiphers_ptr = CtMuxBits.begin();
    auto expand_start = std::chrono::high_resolution_clock::now();
    thread_server_expand(gswCiphers_ptr, half_gsw_ciphertext, context, 0, l, size, galois_keys, l, base_bit, lsk, bsk,
                         sk_gsw_ciphertext);

    cout << "client gap = " << logsize << endl;
    poc_nfllib_ntt_gsw(gswCiphers_ptr[0], context);


    Plaintext test_rlwe_pt("1");
    Ciphertext test_rlwe_ct;
    encryptor1.encrypt_symmetric(test_rlwe_pt, test_rlwe_ct);

    ///steps to crt-compose -> baseB-decompose -> crt-decompose
    vector<uint64_t *> rlwe_decom;
    rwle_decompositions(test_rlwe_ct, context, l, base_bit, rlwe_decom);
    poc_nfllib_ntt_rlwe_decomp(rlwe_decom);



    /// steps for external product. Both rlwe and gsw must be crt-decomposed
    Ciphertext res_ct;
    res_ct.resize(context, context->first_context_data()->parms_id(), 2);
    set_ciphertext(res_ct, context);

    //poc_external_product(gswCiphers_ptr[0], rlwe_decom, context, l, res_ct);
    poc_nfllib_external_product(gswCiphers_ptr[0], rlwe_decom, context, l, res_ct, 1);
    poc_nfllib_intt_ct(res_ct, context);

    int i = 1;
    int duration = 0;
    Plaintext pp;
    decryptor1.decrypt(res_ct, pp);
    cout << "decrypted=" << pp.to_string() << endl;
    cout << "noise budget=" << decryptor1.invariant_noise_budget(res_ct) << endl;


}


void
test_homomorphic_permutation(Evaluator &evaluator1, Encryptor &encryptor1, Decryptor &decryptor1, KeyGenerator &keygen,
                             shared_ptr<SEALContext> context, SecretKey sk) {
    GaloisKeys galois_keys = generate_galois_keys(context, keygen);
    const auto &context_data = context->first_context_data();
    auto &parms = context_data->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data->small_ntt_tables();
    size_t coeff_count = parms.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);
    uint64_t base_bit = 16;
    const int l = 5;
//    context_data->total_coeff_modulus_bit_count() / base_bit;

    //DATA NEEDED  FOR gsw expansion
    int bsk = 16;
    int lsk = context_data->total_coeff_modulus_bit_count() / base_bit;;
    vector<Ciphertext> sk_gsw_ciphertext;
    poc_enc_sk_gsw(sk, context, bsk, sk_gsw_ciphertext);




    //total elements in database and their relevant elements
    uint64_t total_elements = 32;
    uint64_t swapbitsSize = get_swapbits_size(total_elements);

    //spacing of gap needed for packing
    int logsize = ceil(log2(swapbitsSize));
    int gap = ceil(coeff_count / (1 << logsize));
    cout << "swapbits size = " << swapbitsSize << endl;


    //setting up server data
    vector<Ciphertext> server;
    server.resize(total_elements);
    fill_server_bkt(server, total_elements, encryptor1);
    vector<Ciphertext>::iterator input = server.begin();


    //decide the permutation
    vector<uint32_t> permutation;
    permutation.resize(total_elements);
    iota(permutation.begin(), permutation.end(), 0);
    int ttp = total_elements - 1;
    for (int i = 0; i < total_elements; i++) {
        permutation[i] = ttp;
        ttp--;
    }


    int *inverse = computeInversePermutation((int *) permutation.data(), total_elements);
    vector<int> swapbits = sortingNetworkBits(inverse, total_elements);

    assert(swapbits.size() == swapbitsSize);
    Plaintext msg;
    msg.resize(coeff_count);
    for (int h = 0; h < swapbitsSize; h++) {
        msg.data()[h * gap] = ((int64_t) swapbits[h]);
        cout<< "swap bit number=" << h << "="<< swapbits[h]<<endl;
    }

    //encrypted packked permutation cxtx
    vector<Ciphertext> packed_perm_cxtx;
    poc_half_gsw_enc128(l, base_bit, context, sk, packed_perm_cxtx, msg, pool, (1 << logsize));



    vector<GSWCiphertext> CtMuxBits;
    CtMuxBits.resize((1 << logsize), GSWCiphertext(2 * l));


    int size = (1 << logsize);
    vector<GSWCiphertext>::iterator gswCiphers_ptr = CtMuxBits.begin();

    thread_server_expand(gswCiphers_ptr, packed_perm_cxtx, context, 0, l, size, galois_keys, l, base_bit, lsk, bsk,
                         sk_gsw_ciphertext);

    gswCiphers_ptr = CtMuxBits.begin();




    Plaintext test_rlwe_pt("1");
    Ciphertext test_rlwe_ct;
    encryptor1.encrypt_symmetric(test_rlwe_pt, test_rlwe_ct);

    ///steps to crt-compose -> baseB-decompose -> crt-decompose
    vector<uint64_t *> rlwe_decom;
    rwle_decompositions(test_rlwe_ct, context, l, base_bit, rlwe_decom);
    poc_nfllib_ntt_rlwe_decomp(rlwe_decom);


    /// steps for external product. Both rlwe and gsw must be crt-decomposed

    Ciphertext res_ct;
    res_ct.resize(context, context->first_context_data()->parms_id(), 2);


    //poc_external_product(gswCiphers_ptr[0], rlwe_decom, context, l, res_ct);

    for (int i=0;i<swapbitsSize;i++) {
        set_ciphertext(res_ct, context);
        poc_nfllib_ntt_gsw(gswCiphers_ptr[i], context);
        poc_nfllib_external_product(gswCiphers_ptr[i], rlwe_decom, context, l, res_ct, 1);
        poc_nfllib_intt_ct(res_ct, context);

        Plaintext pp;
        decryptor1.decrypt(res_ct, pp);
        cout << "decrypted=" << pp.to_string() << endl;
        cout << "noise budget=" << decryptor1.invariant_noise_budget(res_ct) << endl;
    }


}

void test_plain_expansion(Evaluator &evaluator1, Encryptor &encryptor1, Decryptor &decryptor1, KeyGenerator &keygen,
                          shared_ptr<SEALContext> context, SecretKey sk) {

    GaloisKeys galois_keys = generate_galois_keys(context, keygen);
    const auto &context_data2 = context->first_context_data();
    auto &parms = context_data2->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);
    int logsize = ceil(log2(64));
    int gap = ceil(coeff_count / (1 << logsize));

    int idx = 10;
    Plaintext msg(coeff_count);
    msg.set_zero();
    msg[idx * gap] = 1;

    Plaintext pt;
    pt.resize(coeff_count);
    pt.set_zero();
    pt[0] = 1234567;


    const int base_bits = 20;

    const int decomp_size = parms.plain_modulus().bit_count() / base_bits;

    //gen gsw ct
    GSWCiphertext packed_ct;

    uint64_t dimension_size = 64;

    poc_plain_gsw_enc128(decomp_size, base_bits, context, sk, packed_ct, msg, pool, dimension_size);

    //evaluator1.add_inplace(test,choice_bit[1]);




    vector<GSWCiphertext> list_enc;
    list_enc.resize(dimension_size, GSWCiphertext(3));

    vector<GSWCiphertext>::iterator list_enc_ptr = list_enc.begin();

    auto gsw_enc_time_start = std::chrono::steady_clock::now();

    //poc_expand_flat_threaded(list_enc_ptr, packed_ct, context, dimension_size, galois_keys);
    poc_expand_flat(list_enc_ptr, packed_ct, context, dimension_size, galois_keys);
    auto gsw_enc_time_end = std::chrono::steady_clock::now();


    vector<uint64_t *> plain_decom;
    plain_decompositions(pt, context, decomp_size, base_bits, plain_decom);



    //evaluator1.transform_to_ntt_inplace(pt,context->first_parms_id());
    //evaluator1.transform_to_ntt_inplace(choice_bit[0]);

    poc_nfllib_ntt_rlwe_decomp(plain_decom);

    poc_nfllib_ntt_gsw(list_enc[idx], context);


    Ciphertext res_ct;
    res_ct.resize(context, context->first_context_data()->parms_id(), 2);
    poc_nfllib_external_product(list_enc[idx], plain_decom, context, decomp_size, res_ct, 1);


    poc_nfllib_intt_ct(res_ct, context);


    int duration = duration_cast<std::chrono::microseconds>(gsw_enc_time_end - gsw_enc_time_start).count();
    cout << "Plain prod duration= " << duration << "us" << endl;

    Plaintext ppt;
    decryptor1.decrypt(res_ct, ppt);
    cout << ppt.to_string() << endl;
    cout << "Noise budget after external product " << decryptor1.invariant_noise_budget(res_ct) << endl;

}

void test_rlwe_expansion(Evaluator &evaluator1, Encryptor &encryptor1, Decryptor &decryptor1, KeyGenerator &keygen,
                         shared_ptr<SEALContext> context, SecretKey sk) {

    GaloisKeys galois_keys = generate_galois_keys(context, keygen);
    const auto &context_data2 = context->first_context_data();
    auto &parms = context_data2->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);

    Plaintext msg(coeff_count);


    const int base_bits = 20;

    const int decomp_size = parms.plain_modulus().bit_count() / base_bits;

    //gen gsw ct
    GSWCiphertext packed_ct;

    uint64_t dimension_size = 48;

    msg.set_zero();

    int logsize = ceil(log2(dimension_size));
    int gap = ceil(coeff_count / (1 << logsize));
    msg[1 * gap] = 1;
    //msg[1]=1;

    Ciphertext ct;
    encryptor1.encrypt_symmetric(msg, ct);
    poc_plain_gsw_enc128(decomp_size, base_bits, context, sk, packed_ct, msg, pool, (1 << logsize));

    //evaluator1.add_inplace(test,choice_bit[1]);

    for (int i = 0; i < decomp_size; i++) {
        vector<Ciphertext> list_enc;
        auto gsw_enc_time_start = std::chrono::steady_clock::now();
        list_enc = poc_rlwe_expand(packed_ct[2], context, galois_keys, (1 << logsize));
        auto gsw_enc_time_end = std::chrono::steady_clock::now();
        int duration = duration_cast<std::chrono::milliseconds>(gsw_enc_time_end - gsw_enc_time_start).count();
        cout << "Plain prod duration= " << duration << " ms" << endl;


        //2^40=1099511627776 = hex= 10000000000
        //2^20=1048576 = hex= 100000
        //2^0=1 = hex= 1

        Plaintext ppt;
        decryptor1.decrypt(list_enc[1], ppt);
        cout << ppt.to_string() << endl;
        cout << "Noise budget after external product " << decryptor1.invariant_noise_budget(list_enc[0]) << endl;
    }

}


void test_plain_flatening(Evaluator &evaluator1, Encryptor &encryptor1, Decryptor &decryptor1, KeyGenerator &keygen,
                          shared_ptr<SEALContext> context, SecretKey sk) {

    const auto &context_data2 = context->first_context_data();
    auto &parms = context_data2->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);

    Plaintext msg;
    msg.resize(coeff_count);
    msg.set_zero();
    msg[0] = 1;
//    for (int h = 0; h < coeff_count; h++) {
//        if (h == 0) {
//
//            msg.data()[h] = 1;
//            //cout<< msg.data()[h ]<<endl;
//        } else
//            msg.data()[h] = 0;
//    }

    Plaintext pt;
    pt.resize(coeff_count);
    pt.set_zero();
    pt[0] = 1234567891234567890;
//    for (int h = 0; h < coeff_count; h++) {
//        if (h == 0) {
//
//            pt.data()[h] = 1234567891234567890;
//            //cout<< msg.data()[h ]<<endl;
//        } else
//            pt.data()[h] = 0;
//    }

    const int base_bits = 30;

    const int decomp_size = parms.plain_modulus().bit_count() / base_bits;

    //gen gsw ct
    GSWCiphertext choice_bit;


    poc_plain_gsw_enc128(decomp_size, base_bits, context, sk, choice_bit, msg, pool, 0);

    //evaluator1.add_inplace(test,choice_bit[1]);



    vector<uint64_t *> plain_decom;
    plain_decompositions(pt, context, decomp_size, base_bits, plain_decom);



    //evaluator1.transform_to_ntt_inplace(pt,context->first_parms_id());
    //evaluator1.transform_to_ntt_inplace(choice_bit[0]);

    poc_nfllib_ntt_rlwe_decomp(plain_decom);

    poc_nfllib_ntt_gsw(choice_bit, context);


    auto gsw_enc_time_start = std::chrono::steady_clock::now();

    Ciphertext res_ct;
    res_ct.resize(context, context->first_context_data()->parms_id(), 2);
    poc_nfllib_external_product(choice_bit, plain_decom, context, decomp_size, res_ct, 1);

    auto gsw_enc_time_end = std::chrono::steady_clock::now();

    poc_nfllib_intt_ct(res_ct, context);


    int duration = duration_cast<std::chrono::microseconds>(gsw_enc_time_end - gsw_enc_time_start).count();
    cout << "Plain prod duration= " << duration << "us" << endl;

    Plaintext ppt;
    decryptor1.decrypt(res_ct, ppt);
    cout << ppt.to_string() << endl;
    cout << "Noise budget after external product " << decryptor1.invariant_noise_budget(res_ct) << endl;

}

void test_external_prod_chain(Evaluator &evaluator1, Encryptor &encryptor1, Decryptor &decryptor1, KeyGenerator &keygen,
                              shared_ptr<SEALContext> context, SecretKey sk) {

    const auto &context_data2 = context->first_context_data();
    auto &parms = context_data2->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);


    const int base_bits = 40;
    int iterations = 1000;


    Plaintext gsw_plain(to_string(1));
    Plaintext msg;
    msg.resize(coeff_count);
    msg.set_zero();
    msg[0] = 1;
//    for (int h = 0; h < coeff_count; h++) {
//        if (h == 0) {
//
//            msg.data()[h] = 1;
//            //cout<< msg.data()[h ]<<endl;
//        } else
//            msg.data()[h] = 0;
//    }


    const int lvl = context_data2->total_coeff_modulus_bit_count() / base_bits;

    //gen gsw ct
    GSWCiphertext choice_bit;
    poc_gsw_enc128(lvl, base_bits, context, sk, choice_bit, msg, pool, 0);


    ///gen ct of rlwe
    Plaintext test_rlwe_pt("12345678");
    Ciphertext test_rlwe_ct;
    encryptor1.encrypt_symmetric(test_rlwe_pt, test_rlwe_ct);

    ///steps to crt-compose -> baseB-decompose -> crt-decompose
    vector<uint64_t *> rlwe_decom;
    rwle_decompositions(test_rlwe_ct, context, lvl, base_bits, rlwe_decom);

    Ciphertext res_ct;
    res_ct.resize(context, context->first_context_data()->parms_id(), 2);
    poc_nfllib_external_product(choice_bit, rlwe_decom, context, lvl, res_ct, 1);


    cout << "-----------------------------------------------" << endl;
    cout << " Testing external product chain " << endl;
    cout << "-----------------------------------------------" << endl;


    GSWCiphertext chain_gsw;
    poc_gsw_enc128(lvl, base_bits, context, sk, chain_gsw, msg, pool, 0);

    int i = 1;
    while (decryptor1.invariant_noise_budget(res_ct) > 0 && i < iterations) {
        i++;

        Plaintext pp;
        decryptor1.decrypt(res_ct, pp);
        cout << "Noise budget after " << i << " external product " << decryptor1.invariant_noise_budget(res_ct) << endl;

        vector<uint64_t *> rlwe_decom;
        rwle_decompositions(res_ct, context, lvl, base_bits, rlwe_decom);

        poc_nfllib_external_product(choice_bit, rlwe_decom, context, lvl, res_ct, 1);
        poc_nfllib_intt_ct(res_ct, context);
        for (auto p : rlwe_decom) {
            free(p);
        }
        rlwe_decom.clear();

    }

    cout << "-----------------------------------------------" << endl;
    cout << " Testing bfv product chain " << endl;
    cout << "-----------------------------------------------" << endl;


    Plaintext left_rlwe_pt("1");
    Ciphertext left_rlwe_ct;
    encryptor1.encrypt_symmetric(left_rlwe_pt, left_rlwe_ct);


    Ciphertext res_ct_;


    i = 0;
    while (decryptor1.invariant_noise_budget(test_rlwe_ct) > 0 && i < iterations) {


        evaluator1.multiply_inplace(test_rlwe_ct, left_rlwe_ct);
        cout << "Noise budget after " << i << " product " << decryptor1.invariant_noise_budget(test_rlwe_ct) << endl;
        i++;
    }
    if (i == 0)
        cout << "Noise budget after " << i << " product " << decryptor1.invariant_noise_budget(test_rlwe_ct) << endl;

}


void test_seal(Evaluator &evaluator1, Encryptor &encryptor1, Decryptor &decryptor1, KeyGenerator &keygen,
               shared_ptr<SEALContext> context, SecretKey sk) {

    const auto &context_data2 = context->first_context_data();
    const seal::EncryptionParameters &parms = context_data2->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);


    Plaintext gsw_plain(to_string(1));
    Plaintext msg;
    msg.resize(coeff_count);
    for (int h = 0; h < coeff_count; h++) {

        msg.data()[h] = 1;

    }

    Ciphertext test_rlwe_ct, temp;
    encryptor1.encrypt_symmetric(msg, test_rlwe_ct);


    evaluator1.transform_to_ntt_inplace(msg, context->first_parms_id());
    evaluator1.transform_to_ntt_inplace(test_rlwe_ct);
    //evaluator1.transform_to_ntt_inplace(test_rlwe_ct);

    int duration = 0;
    int iterations = 1000;
    Ciphertext res_ct;
    for (int i = 0; i < iterations; i++) {


        auto gsw_enc_time_start = std::chrono::high_resolution_clock::now();
        //poc_nfllib_ntt_ct(test_rlwe_ct, context);
        evaluator1.multiply_plain_inplace(test_rlwe_ct, msg);
        //poc_nfllib_plain_ct_prod(test_rlwe_ct , msg, context, res_ct);
        auto gsw_enc_time_end = std::chrono::high_resolution_clock::now();
        //poc_nfllib_intt_ct(test_rlwe_ct, context);

        duration = duration + duration_cast<std::chrono::microseconds>(gsw_enc_time_end - gsw_enc_time_start).count();
    }

    //evaluator1.transform_from_ntt_inplace(res_ct);
    Plaintext ppt;
    //decryptor1.decrypt(res_ct, ppt);
    //cout<< ppt.to_string()<<endl;

    //evaluator1.transform_to_ntt_inplace(msg,context->first_parms_id());



    //poc_nfllib_plain_ct_prod2(test_rlwe_ct,msg,context,temp);
    //poc_nfllib_plain_ct_prod(test_rlwe_ct,msg,context,temp);
    //evaluator1.multiply_plain(test_rlwe_ct,msg,temp);

    //evaluator1.multiply_plain(test_rlwe_ct,msg,temp);



    cout << "prod duration= " << duration / iterations << " us" << endl;
    //evaluator1.add_inplace(test_rlwe_ct,temp);

}

//int main() {
//
//    EncryptionParameters parms(scheme_type::BFV);
//    set_bfv_parms(parms);
//    auto context = SEALContext::Create(parms);
//    print_line(__LINE__);
//    print_parameters(context);
//
//    KeyGenerator keygen(context);
//
//    //generating secret key
//    Plaintext secret_key_pt;
//    SecretKey secret_key = keygen.secret_key();
//
//
//    /// generating encryptor, decryptor and evaluator
//    Encryptor encryptor(context, secret_key);
//    Decryptor decryptor(context, secret_key);
//    Evaluator evaluator(context);
//
////    test_gsw_expansion(evaluator, encryptor, decryptor, keygen, context, secret_key);
//    test_homomorphic_permutation(evaluator, encryptor, decryptor, keygen, context, secret_key);
//    //test_rlwe_expansion(evaluator, encryptor, decryptor, keygen,  context, secret_key);
//    //test_plain_expansion(evaluator, encryptor, decryptor, keygen,  context, secret_key);
//    //test_nfllib_ct_add(evaluator, encryptor, decryptor, keygen,  context, secret_key);
//    //test_plain_flatening(evaluator, encryptor, decryptor, keygen,  context, secret_key);
//    //test_external_prod_with_sk(evaluator, encryptor, decryptor, keygen,  context, secret_key);
//    //test_external_prod(evaluator, encryptor, decryptor, keygen,  context, secret_key);
//    //test_external_prod_chain(evaluator, encryptor, decryptor, keygen,  context, secret_key);
//    //test_seal(evaluator, encryptor, decryptor, keygen,  context, secret_key);
//    return 0;
//}

int main(){

    uint64_t number_of_items = 1<<17;
    uint64_t size_per_item = 256; // in bytes
    uint32_t N = 4096;

    // Recommended values: (logt, d) = (12, 2) or (8, 1).
    //uint32_t logt = 60;
    PirParams pir_params;


    EncryptionParameters parms(scheme_type::BFV);
    set_bfv_parms(parms);

    int n = 1 << 16;

    int num_equations = n;
    int dimension = 128;
    double eps = 0.45;
    int band_width = 60;
    int num_chunks = 1;
    //uint64_t field_size = 7;
    uint64_t field_size = 1152921504606830593ULL;
    uint64_t t = field_size;
    uint64_t logt = 60;
    uint64_t element_size = size_per_item;
    cout << field_size << endl;

    int num_bins = ceil(((1 + eps) * n / dimension / num_chunks));
    cout << "Num bins for key allocation: " << num_bins << "\n";
    BinAllocator allocator(num_bins, dimension, num_chunks);
    // We can't afford to assign up to dimension keys to each bin, since it's likely that the associated matrix is not
    // full rank. We give some buffer, e.g. up to dimension - bin_buffer keys can be assigned to each bin. Thus, the
    // linear system has more variables than equations.
    int bin_buffer = static_cast<int>(dimension * 0.05);
    vector<oc::block> keys = GenerateKeys(n);

    if (!allocator.SetKeys(keys, bin_buffer, 5000)) {
      cout << "Failed to find appropriate allocations!" << endl;
      return 0;
    }

    cout << "Found bin allocation.\n";

    cout << "Num bins: " << allocator.NumBins() << endl;

    vector<string> raw_values = GenerateRawValues(n, element_size);
    vector<PlaintextValue> values = ToPlaintextValues(t, logt, raw_values);
    cout << "Num Coeffs: " << values[0].GetValue().size() << endl;

    BucketedLinearSystems<PlaintextValue, BandRowVectorHasher> linear_systems;
    if (!linear_systems.SetKeysAndValues(allocator, keys, values, dimension, band_width, band_width, field_size)) {
      cout << "fail" << endl;
      exit(0);
    }

    vector<vector<vector<PlaintextValue>>> res = linear_systems.SolveLinearSystems();

    gen_params( number_of_items,  size_per_item, N, logt,
                pir_params, res.size() * res[0][0].size());

    for (int d : pir_params.nvec) {
      cout << d << " ";
    }
    cout << endl;

    vector<Plaintext> flattened = ToFlattenedDb(res, N, pir_params.nvec);
    assert(res.size() == 1);
    assert(res[0].size() == 1);
    //vector<Plaintext> flattened = ToFlattenedDb3(res[0][0], N, pir_params.nvec);
    cout << "Flattened size: " << flattened.size() << endl;

    cout << "Main: Initializing the database (this may take some time) ..." << endl;

    // Initialize PIR Server
    cout << "Main: Initializing server and client" << endl;
    pir_server server(parms, pir_params);

    // Initialize PIR client....
    pir_client client(parms, pir_params);
    GaloisKeys galois_keys = client.generate_galois_keys();

    cout << "Main: Setting Galois keys...";
    server.set_galois_key(0, galois_keys);

    auto time_pre_s = high_resolution_clock::now();
    //server.set_database(move(db), number_of_items, size_per_item);
    cout << "Main: Preprocessing database..." << endl;
    unique_ptr<vector<Plaintext>> flattened_db = make_unique<vector<Plaintext>>(std::move(flattened));
    server.set_database(std::move(flattened_db), true);
    server.preprocess_database();
    cout << "Main: Database preprocessing done" << endl;
    auto time_pre_e = high_resolution_clock::now();
    auto time_pre_us = duration_cast<microseconds>(time_pre_e - time_pre_s).count();

    random_device rd;
    uniform_int_distribution<int> dist(0, n);
    int key_index = dist(rd);
    cout << "Key index: " << key_index << endl;
    int index = allocator.AllocateFlattened(keys[key_index]);
    auto row = linear_systems.GetHasher().HashToBandRowVector(keys[key_index]);
    uint64_t* raw_band = row.RawBand();
    vector<int> row_vector_indices;
    for (int i = 0; i < row.GetLength(); i++) {
      if (raw_band[i] != 0) {
        row_vector_indices.push_back(row.GetOffset() + i);
      }
    }

    PirQuery query = client.generate_query_combined2(row_vector_indices, index);

    cout<<"Main: query size = "<< query.size()<< endl;

    SecretKey sk = client.get_decryptor();

    GSWCiphertext enc_sk=client.get_enc_sk();
    server.set_enc_sk(enc_sk);

    auto time_server_s = high_resolution_clock::now();
    PirReply reply = server.generate_reply_combined(query, 0, sk);
    auto time_server_e = high_resolution_clock::now();
    auto time_server_us = duration_cast<microseconds>(time_server_e - time_server_s).count();

    Plaintext rep = client.decrypt_result(reply);

    cout << "Original coeffs: ";
    for (uint64_t field : values[key_index].GetValue()) {
      cout << field << " ";
    }

    cout << endl;
    vector<uint64_t> coeffs = PlaintextToCoeffs(t, rep);
    cout << "Coeff size: " << coeffs.size() << endl;
    cout << "Retrieved coeffs: ";
    for (int i = 0; i < N; ++i) {
      cout << coeffs[i] % field_size << " ";
    }
    cout << endl;
    cout << "Diff: ";
    for (int i = 0; i < values[key_index].GetValue().size(); ++i) {
      auto a = coeffs[i] % field_size;
      auto b = values[key_index].GetValue()[i];
      if (a > b) {
        cout << a - b << " ";
      } else {
        cout << b - a << " ";
      }
    }
    cout << endl;

    int idx = Contains(coeffs, values[key_index].GetValue());
    if (idx != -1) {
      cout << "Correct!" << endl;
      string value = CoeffsToValue(t, logt, idx, coeffs, element_size);
      cout << value << endl;
    } else {
      cout << "Incorrect!" << endl;
      //string value = CoeffsToValue(t, logt, idx, coeffs, element_size);
      //cout << value << endl;
    }

    return 0;
}