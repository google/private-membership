//
// Created by Haris Mughees on 4/14/21.
//

#include "external_prod.h"
poly_t *resa = alloc_aligned<poly_t, 32>(1),
        *resb = alloc_aligned<poly_t, 32>(1),
        *resc = alloc_aligned<poly_t, 32>(1);


void set_bfv_parms(EncryptionParameters &parms) {

    // size of polynomial degree
    size_t poly_modulus_degree = 4096;

    //number of bits in plaintext modulu, seal does not allow plaintexts to be more than 62 bits
    uint64_t plaintxt_mod_bits= 60;


    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, plaintxt_mod_bits));
    //parms.set_plain_modulus((1LL << 60) + 33);

    cout << "Plaintext Modulus: " << parms.plain_modulus().value() << endl;

    //we use special CRT components of coeff mod that are compatible with nfllib
    parms.set_coeff_modulus({4611686018326724609, 4611686018309947393, 4611686018282684417});
    //parms.set_coeff_modulus({ 1073479681, 1072496641,1071513601});
}



void poc_gsw_enc128(const uint64_t l, const uint64_t base_bit, shared_ptr<SEALContext> context, const SecretKey sk,
                    vector<Ciphertext> &gsw_ciphertext, Plaintext gsw_plain, seal::util::MemoryPool &pool,
                    uint64_t inv) {
    Encryptor encryptor(context, sk);
    Decryptor decryptor(context, sk);
    const auto &context_data = context->first_context_data();
    auto &parms = context_data->parms();
    auto &coeff_modulus = parms.coeff_modulus();


    size_t ct_poly_count = context_data->parms().coeff_modulus().size();/// find good way of getting it
    int total_bits;
    uint64_t r_l = l;

    Ciphertext t;
    for (int j = 0; j < 2; j++) {// c0, c1
        total_bits = (context_data->total_coeff_modulus_bit_count());
        for (int p = 0; p < r_l; p++) {
            const int shift_amount = ((total_bits) - ((p + 1) * base_bit));
            Ciphertext res;
            encryptor.encrypt_zero_symmetric(res);
            //set_ciphertext(res, context);

            poc_multiply_add_plain_without_scaling_variant(gsw_plain, *context->first_context_data(), shift_amount,
                                                           res.data(j), pool, inv);


            gsw_ciphertext.push_back(res);
        }

    }

}

void poc_multiply_add_plain_without_scaling_variant(const Plaintext &plain, const SEALContext::ContextData &context_data,
                                                    const int shift_amount, uint64_t *destination,
                                                    seal::util::MemoryPool &pool, uint64_t inv = 0) {
    auto &parms = context_data.parms();
    size_t coeff_count = parms.poly_modulus_degree();
    size_t plain_coeff_count = plain.coeff_count();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_mod_count = coeff_modulus.size();
    auto plain_modulus = context_data.parms().plain_modulus();
    auto coeff_div_plain_modulus = context_data.coeff_div_plain_modulus();
    uint64_t h;


    //cout<< shift_amount << endl;

    for (size_t i = 0; i < plain_coeff_count; i++) {
        // Add to ciphertext: h * m
        for (size_t j = 0; j < coeff_mod_count; j++) {
            //init empty 128 bit integers
            auto ptr(allocate_uint(coeff_mod_count, pool));
            auto ptr2(allocate_uint(coeff_mod_count, pool));
            auto ptr3(allocate_uint(coeff_mod_count, pool));
            //set 1 in lsb (it will be used for bit shifts)

            uint64_t poly_inv;
            uint64_t plain_coeff;
            if (inv > 0) {
                seal::util::try_invert_uint_mod(inv, coeff_modulus[j], poly_inv);
                plain_coeff = seal::util::multiply_uint_uint_mod(plain.data()[i], poly_inv, coeff_modulus[j]);

            } else {
                plain_coeff = plain.data()[i];
            }


            ptr2[0] = 0;
            ptr2[1] = 0;
            ptr[0] = 1;
            ptr[1] = 0;
            //use 128 bit implementation for left shifts 1<<shiftamount
            util::left_shift_uint128(ptr.get(), shift_amount, ptr2.get());
            h = seal::util::barrett_reduce_128(ptr2.get(), coeff_modulus[j]);

            //seal::util::multiply_uint_uint_mod(plain_coeff, context_data.coeff_div_plain_modulus(),coeff_modulus[j]);

            //barret reduction is used for converting 128 bit interger to mod q1, q2 where q1, q2 are max 64 bits

            h = seal::util::multiply_uint_uint_mod(h, plain_coeff, coeff_modulus[j]);
            //cout<<h<<",";
            destination[i + (j * coeff_count)] = seal::util::add_uint_uint_mod(
                    destination[i + (j * coeff_count)], h, coeff_modulus[j]);
        }
    }

}

void
rwle_decompositions(Ciphertext rlwe_ct_1, shared_ptr<SEALContext> context, const uint64_t l, const uint64_t base_bit,
                    vector<uint64_t *> &rlwe_decom) {
    const auto &context_data2 = context->first_context_data();
    auto &parms2 = context_data2->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms2.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);


    //compose ciphertext for all q_i's
    context_data2->rns_tool()->base_q()->compose_array(rlwe_ct_1.data(0), coeff_count, pool);
    context_data2->rns_tool()->base_q()->compose_array(rlwe_ct_1.data(1), coeff_count, pool);


    //128 bits decomp as given in external product
    poc_decomp_rlwe128(rlwe_ct_1, l, context, rlwe_decom, base_bit, pool);

    //auto rlwe_start = std::chrono::high_resolution_clock::now();
    int ssize = rlwe_decom.size();
    for (int i = 0; i < ssize; i++) {
        //rwle_crt_decompose and poc_decompose_array does same thing but rwle_crt_decompose is slower
        //rwle_crt_decompose(rlwe_decom[i], context, pool);
        //cout<<i<<endl;
        poc_decompose_array(rlwe_decom[i], coeff_count, coeff_modulus, coeff_modulus_size, pool);

    }

}

void poc_decomp_plain(Plaintext pt, const uint64_t decomp_size, shared_ptr<SEALContext> context,
                        vector<uint64_t *> &vec_ciphertexts, int base_bit, seal::util::MemoryPool &pool) {

    assert(vec_ciphertexts.size() == 0);
    const uint64_t base = UINT64_C(1) << base_bit;
    const uint64_t mask = base - 1;

    const auto &context_data = context->get_context_data(context->first_parms_id());
    auto &parms = context_data->parms();
    auto &plain_modulus = parms.plain_modulus();
    size_t coeff_count = parms.poly_modulus_degree();
    size_t coeff_mod_count = parms.coeff_modulus().size();
    size_t pt_poly_count = 1;// there is only one poilynomial in plaintext


    uint64_t r_l = decomp_size;
    int total_bits;
    std::uint64_t *res;

    for (int j = 0; j < pt_poly_count; j++) {// c0, c1
        total_bits = (plain_modulus.bit_count()); // in normal rlwe decomp we use total bits of q, here total bits of t is required
        uint64_t *raw_ptr = pt.data();


        for (int p = 0; p < r_l; p++) {
            vector<uint64_t *> results;
            res = (std::uint64_t *) calloc((coeff_count * coeff_mod_count), sizeof(uint64_t)); //we are allocating larger space to cater for ct modulus later
            const int shift_amount = ((total_bits) - ((p + 1) * base_bit));


            for (size_t k = 0; k < coeff_count; k++) {
                auto ptr(allocate_uint(1, pool));
                ptr[0] = 0;


                //seal::util::right_shift_uint(&raw_ptr[k],shift_amount,1,ptr.get());
                seal::util::right_shift_uint128(&raw_ptr[k], shift_amount, ptr.get());
                uint64_t temp1 = ptr[0] & mask;
                res[k*coeff_mod_count] = temp1;

            }
            //results.push_back(res);
            vec_ciphertexts.push_back(res);
        }

    }

}


void
plain_decompositions(Plaintext &pt, shared_ptr<SEALContext> &context, const uint64_t decomp_size, const uint64_t base_bit,
                     vector<uint64_t *> &plain_decom) {
    const auto &context_data2 = context->first_context_data();
    auto &parms2 = context_data2->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms2.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);


    //compose ciphertext for all q_i's
//s    context_data2->rns_tool()->base_q()->compose_array(rlwe_ct_1.data(0), coeff_count, pool);
//    context_data2->rns_tool()->base_q()->compose_array(rlwe_ct_1.data(1), coeff_count, pool);


    //128 bits decomp as given in external product


    poc_decomp_plain(pt, decomp_size, context, plain_decom, base_bit, pool);

    //auto rlwe_start = std::chrono::high_resolution_clock::now();
    int ssize = plain_decom.size();
    for (int i = 0; i < ssize; i++) {
        //rwle_crt_decompose and poc_decompose_array does same thing but rwle_crt_decompose is slower
        //rwle_crt_decompose(rlwe_decom[i], context, pool);
        //cout<<i<<endl;
        poc_decompose_array(plain_decom[i], coeff_count, coeff_modulus, coeff_modulus_size, pool);

    }

}

void poc_decomp_rlwe128(Ciphertext ct, const uint64_t l, shared_ptr<SEALContext> context,
                        vector<uint64_t *> &vec_ciphertexts, int base_bit, seal::util::MemoryPool &pool) {

    assert(vec_ciphertexts.size() == 0);
    const uint64_t base = UINT64_C(1) << base_bit;
    const uint64_t mask = base - 1;

    const auto &context_data = context->get_context_data(ct.parms_id());
    auto &parms = context_data->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_count = parms.poly_modulus_degree();
    size_t coeff_mod_count = coeff_modulus.size();
    size_t ct_poly_count = ct.size();

    uint64_t r_l = l;
    int total_bits;
    std::uint64_t *res;

    for (int j = 0; j < ct_poly_count; j++) {// c0, c1
        total_bits = (context_data->total_coeff_modulus_bit_count());
        uint64_t *encrypted1_ptr = ct.data(j);

        for (int p = 0; p < r_l; p++) {
            vector<uint64_t *> results;
            res = (std::uint64_t *) calloc((coeff_count * coeff_mod_count), sizeof(uint64_t));
            const int shift_amount = ((total_bits) - ((p + 1) * base_bit));
            for (size_t k = 0; k < coeff_mod_count * coeff_count; k = k + 2) {
                auto ptr(allocate_uint(2, pool));
                ptr[0] = 0;
                ptr[1] = 0;
                seal::util::right_shift_uint128(&encrypted1_ptr[k], shift_amount, ptr.get());
                uint64_t temp1 = ptr[0] & mask;
                res[k] = temp1;

            }
            //results.push_back(res);
            vec_ciphertexts.push_back(res);
        }

    }

}

void poc_decompose_array(uint64_t *value, size_t count, std::vector<Modulus> coeff_modulus, size_t coeff_mod_count,
                         MemoryPoolHandle pool) {
    if (!value) {
        throw invalid_argument("value cannot be null");
    }
    if (!pool) {
        throw invalid_argument("pool is uninitialized");
    }

    if (coeff_mod_count > 1) {
        if (!product_fits_in(count, coeff_mod_count)) {
            throw logic_error("invalid parameters");
        }

        // Decompose an array of multi-precision integers into an array of arrays,
        // one per each base element
        auto value_copy(allocate_uint(count * coeff_mod_count, pool));

        auto temp_array(allocate_uint(count * coeff_mod_count, pool));

        // Merge the coefficients first
        for (size_t i = 0; i < count; i++) {
            for (size_t j = 0; j < coeff_mod_count; j++) {
                temp_array[j + (i * coeff_mod_count)] = value[j + (i * coeff_mod_count)];
            }
        }

        set_zero_uint(count * coeff_mod_count, value);

        for (size_t i = 0; i < count; i++) {
            //set_uint_uint(value, size_, value_copy.get());

            // Temporary space for 128-bit reductions
            for (size_t j = 0; j < coeff_mod_count; j++) {
                // Reduce in blocks
                uint64_t temp[2]{0, temp_array[(i * coeff_mod_count) + coeff_mod_count - 1]};
                for (size_t k = coeff_mod_count - 1; k--;) {
                    temp[0] = temp_array[(i * coeff_mod_count) + k];
                    temp[1] = barrett_reduce_128(temp, coeff_modulus[j]);
                }

                // Save the result modulo i-th base element
                //value[i] = temp[1];
                value[(j * count) + i] = temp[1];
            }
        }
    }

}

GaloisKeys generate_galois_keys(shared_ptr<SEALContext> context, KeyGenerator &keygen) {
    // Generate the Galois keys needed for coeff_select.
    //check https://github.com/microsoft/SealPIR/blob/dd06e2ff10d966177dc9892e92852e60add3f8b6/pir_client.cpp#L137
    //In the polynomial view (not batching), a Galois automorphism by a Galois element p changes
    //        Enc(plain(x)) to Enc(plain(x^p)).

    std::vector<uint32_t> galois_elts;
    const auto &context_data = context->first_context_data();
    auto &parms = context_data->parms();
    int N = parms.poly_modulus_degree();
    int logN = seal::util::get_power_of_two(N);
    for (int i = 0; i < logN; i++) {
        galois_elts.push_back((N + seal::util::exponentiate_uint64(2, i)) / seal::util::exponentiate_uint64(2, i));

    }

    return keygen.galois_keys_local(galois_elts);
}


void poc_plain_gsw_enc128(const uint64_t decomp_size, const uint64_t base_bit, shared_ptr<SEALContext> context,
                         const SecretKey sk, vector<Ciphertext> &gsw_ciphertext, Plaintext gsw_plain,
                         seal::util::MemoryPool &pool, uint64_t inv) {

    Encryptor encryptor(context, sk);
    Decryptor decryptor(context, sk);
    const auto &context_data = context->first_context_data();
    auto &parms = context_data->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    auto &plain_modulus = parms.plain_modulus();

    //size_t ct_poly_count = context_data->parms().coeff_modulus().size();/// find good way of getting it
    int total_bits;
    uint64_t r_l = decomp_size;
    total_bits = (plain_modulus.bit_count()); // in normal gsw we use total bits of q, here total bits of t is required
    auto ptr(allocate_uint(2, pool));
    auto ptr2(allocate_uint(2, pool));
    Plaintext   ppt(gsw_plain.coeff_count());
    ppt.set_zero();

    uint64_t h=1;
    for (int p = 0; p < r_l; p++) {

        const int shift_amount = ((total_bits) - ((p + 1) * base_bit));
        ptr2[0] = 0;
        ptr2[1] = 0;
        ptr[0] = 1;
        ptr[1] = 0;
        util::left_shift_uint128(ptr.get(), shift_amount, ptr2.get());
        h = seal::util::barrett_reduce_128(ptr2.get(), plain_modulus.value());
        if(inv>0){
            uint64_t tt;
            util::try_invert_uint_mod(inv, parms.plain_modulus().value(),tt);
            h = util::multiply_uint_uint_mod(h ,tt , parms.plain_modulus());
        }

        for (int j=0;j<gsw_plain.coeff_count();j++){

            if(gsw_plain.data()[j]==1) {
                //cout <<"here i am "<< j << endl;
                ppt[j] = h;
            }


        }
        Ciphertext res;
        encryptor.encrypt_symmetric(ppt,res);
        //encryptor.encrypt_zero_symmetric(res);
        //set_ciphertext(res, context);
        //poc_multiply_add_plain_without_scaling_variant(gsw_plain, *context->first_context_data(), shift_amount,res.data(0), pool, inv);
        gsw_ciphertext.push_back(res);



    }



}


void poc_plain_gsw_enc128_combined(const uint64_t decomp_size, const uint64_t base_bit, shared_ptr<SEALContext> context,
                          const SecretKey sk, vector<Ciphertext> &gsw_ciphertext, Plaintext gsw_plain,
                          seal::util::MemoryPool &pool, uint64_t inv, uint64_t gap) {

    Encryptor encryptor(context, sk);
    Decryptor decryptor(context, sk);
    const auto &context_data = context->first_context_data();
    auto &parms = context_data->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    auto &plain_modulus = parms.plain_modulus();

    //size_t ct_poly_count = context_data->parms().coeff_modulus().size();/// find good way of getting it
    int total_bits;
    uint64_t r_l = decomp_size;
    total_bits = (plain_modulus.bit_count()); // in normal gsw we use total bits of q, here total bits of t is required
    auto ptr(allocate_uint(2, pool));
    auto ptr2(allocate_uint(2, pool));
    Plaintext   ppt(gsw_plain.coeff_count());
    ppt.set_zero();

    uint64_t h=1;
    for (int p = 0; p < r_l; p++) {

        const int shift_amount = ((total_bits) - ((p + 1) * base_bit));
        ptr2[0] = 0;
        ptr2[1] = 0;
        ptr[0] = 1;
        ptr[1] = 0;
        util::left_shift_uint128(ptr.get(), shift_amount, ptr2.get());
        h = seal::util::barrett_reduce_128(ptr2.get(), plain_modulus.value());
        if(inv>0){
            uint64_t tt;
            util::try_invert_uint_mod(2* inv, parms.plain_modulus().value(),tt);
            h = util::multiply_uint_uint_mod(h ,tt , parms.plain_modulus());
        }

        uint64_t  total_dim_with_gap = inv*gap;
        for (int j=0;j<total_dim_with_gap;j++){

            if(gsw_plain.data()[j]==1) {

                ppt[j+p*(total_dim_with_gap)] = h;

            }


        }




    }
    Ciphertext res;
    encryptor.encrypt_symmetric(ppt,res);
    //encryptor.encrypt_zero_symmetric(res);
    //set_ciphertext(res, context);
    //poc_multiply_add_plain_without_scaling_variant(gsw_plain, *context->first_context_data(), shift_amount,res.data(0), pool, inv);
    gsw_ciphertext.push_back(res);



}



void poc_expand_flat(vector<GSWCiphertext>::iterator &result, vector<Ciphertext> &packed_swap_bits,
                          shared_ptr<SEALContext> context, int size, seal::GaloisKeys &galkey ) {


    const auto &context_data2 = context->first_context_data();
    auto &parms2 = context_data2->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms2.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);
    int logN = seal::util::get_power_of_two(coeff_count);
    vector<Ciphertext> expanded_ciphers(coeff_count);

    //outloop is from 0-to-(l-1)
    for (int i = 0; i < packed_swap_bits.size(); i++) {
        auto rlwe_start = std::chrono::high_resolution_clock::now();
        expanded_ciphers = poc_rlwe_expand(packed_swap_bits[i], context, galkey, size);
        auto rlwe_end = std::chrono::high_resolution_clock::now();

        vector<uint64_t *> rlwe_decom;
        for (int j = 0; j < size; j++) {
            ///put jth expanded ct in ith idx slot  of jt gswct
            result[j][i] = expanded_ciphers[j];


        }



        //cout <<"---rlwe---"<< duration_cast<std::chrono::milliseconds >(rlwe_end - rlwe_start).count() <<"ms"<< endl;

    }


}

void poc_expand_flat_combined(vector<GSWCiphertext>::iterator &result, vector<Ciphertext> &packed_swap_bits,
                     shared_ptr<SEALContext> context, int size, seal::GaloisKeys &galkey ) {


    const auto &context_data2 = context->first_context_data();
    auto &parms2 = context_data2->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms2.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);
    int logN = seal::util::get_power_of_two(coeff_count);
    vector<Ciphertext> expanded_ciphers(coeff_count);

    //outloop is from 0-to-(l-1)


    assert(packed_swap_bits.size()==1);
        auto rlwe_start = std::chrono::high_resolution_clock::now();
        expanded_ciphers = poc_rlwe_expand(packed_swap_bits[0], context, galkey, 2* size);
        auto rlwe_end = std::chrono::high_resolution_clock::now();


    for (int i = 0; i < 2; i++) {

        //components of first dimension
        vector<uint64_t *> rlwe_decom;
        for (int j = 0; j < size; j++) {
            ///put jth expanded ct in ith idx slot  of jt gswct
            result[j][i] = expanded_ciphers[j+i*size];


        }



        //cout <<"---rlwe---"<< duration_cast<std::chrono::milliseconds >(rlwe_end - rlwe_start).count() <<"ms"<< endl;

    }


}

void poc_rlwe_expand_threaded(Ciphertext packedquery, shared_ptr<SEALContext> context, seal::GaloisKeys galkey, uint64_t size, vector<Ciphertext> &result) {
    /// this function return size  vector of RLWE ciphertexts
    /// it takes a single RLWE packed ciphertext

    //result.clear();
    //result.resize(size);

    cout<<result.size()<<endl;



    Evaluator evaluator1(context);

    const auto &context_data = context->first_context_data();

    auto &parms = context_data->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    int N2 = parms.poly_modulus_degree();


    Ciphertext tempctxt;
    Ciphertext tempctxt_rotated;
    Ciphertext tempctxt_shifted;


    vector<Ciphertext> temp;
    Ciphertext tmp;
    temp.push_back(packedquery);
    int numIters = ceil(log2(size));   // size is a ---rlwe---power of 2.


    if (numIters > ceil(log2(N2))) {
        throw logic_error("m > coeff_count is not allowed.");
    }

    int startIndex = static_cast<int>(log2(N2) - numIters);



    auto time_server_s = high_resolution_clock::now();
    for (long i = 0; i < numIters; i++) {
        vector<Ciphertext> newtemp(temp.size() << 1);
        int index = startIndex + i;
        int power = (N2 >> index) + 1;//k
        int ai = (1 << index);
        for (int j = 0; j < (1 << i); j++) {

            //check this power
            // tempctxt_rotated = subs(result[j])
            //evaluator1.apply_galois(temp[j], galois_elts[i], galkey, tempctxt_rotated);

            evaluator1.apply_galois(temp[j], power, galkey, tempctxt_rotated);

            // result[j+ 2**i] = result[j] - tempctxt_rotated;

            evaluator1.sub(temp[j], tempctxt_rotated, newtemp[j + (1 << i)]);

            // divide by x^ai = multiply by x^(2N - ai).

            multiply_power_of_X(newtemp[j + (1 << i)], tempctxt_shifted, (N2 << 1) - ai, context);

            newtemp[j + (1 << i)] = tempctxt_shifted;

            evaluator1.add(tempctxt_rotated, temp[j], newtemp[j]);

        }

        temp = newtemp;
    }
    cout<<"here"<<endl;
 //   auto time_server_e = high_resolution_clock::now();
//    int dur =  duration_cast<microseconds>(time_server_e - time_server_s).count();

    //cout<<"inner loop=========="<<dur<<endl;



   // result.resize(temp.size() );

    //std::move(std::begin(temp), std::end(temp), std::back_inserter(result));
    //assert(temp.size()==result.size());
    //result=temp;

}
vector<Ciphertext> poc_rlwe_expand(Ciphertext packedquery, shared_ptr<SEALContext> context, seal::GaloisKeys galkey, uint64_t size) {
    /// this function return size  vector of RLWE ciphertexts
    /// it takes a single RLWE packed ciphertext

    Evaluator evaluator1(context);

    const auto &context_data = context->first_context_data();

    auto &parms = context_data->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    int N2 = parms.poly_modulus_degree();


    Ciphertext tempctxt;
    Ciphertext tempctxt_rotated;
    Ciphertext tempctxt_shifted;


    vector<Ciphertext> temp;
    Ciphertext tmp;
    temp.push_back(packedquery);
    int numIters = ceil(log2(size));   // size is a ---rlwe---power of 2.


    if (numIters > ceil(log2(N2))) {
        throw logic_error("m > coeff_count is not allowed.");
    }

    int startIndex = static_cast<int>(log2(N2) - numIters);



    auto time_server_s = high_resolution_clock::now();
    for (long i = 0; i < numIters; i++) {
        vector<Ciphertext> newtemp(temp.size() << 1);
        int index = startIndex + i;
        int power = (N2 >> index) + 1;//k
        int ai = (1 << index);
        for (int j = 0; j < (1 << i); j++) {

            //check this power
            // tempctxt_rotated = subs(result[j])
            //evaluator1.apply_galois(temp[j], galois_elts[i], galkey, tempctxt_rotated);

            evaluator1.apply_galois(temp[j], power, galkey, tempctxt_rotated);

            // result[j+ 2**i] = result[j] - tempctxt_rotated;

            evaluator1.sub(temp[j], tempctxt_rotated, newtemp[j + (1 << i)]);

            // divide by x^ai = multiply by x^(2N - ai).

            multiply_power_of_X(newtemp[j + (1 << i)], tempctxt_shifted, (N2 << 1) - ai, context);

            newtemp[j + (1 << i)] = tempctxt_shifted;

            evaluator1.add(tempctxt_rotated, temp[j], newtemp[j]);

        }
        temp = newtemp;
    }

    auto time_server_e = high_resolution_clock::now();
    int dur =  duration_cast<microseconds>(time_server_e - time_server_s).count();

    //cout<<"inner loop=========="<<dur<<endl;

    return temp;

}


void multiply_power_of_X(const Ciphertext &encrypted, Ciphertext &destination,
                         uint32_t index, shared_ptr<SEALContext> context) {

    const auto &context_data = context->first_context_data();
    auto &parms = context_data->parms();
    auto coeff_mod_count = parms.coeff_modulus().size();
    auto coeff_count = parms.poly_modulus_degree();
    auto encrypted_count = encrypted.size();


    destination = encrypted;


    for (int i = 0; i < encrypted_count; i++) {
        for (int j = 0; j < coeff_mod_count; j++) {
            seal::util::negacyclic_shift_poly_coeffmod(encrypted.data(i) + (j * coeff_count),
                                                       coeff_count, index,
                                                       parms.coeff_modulus()[j],
                                                       destination.data(i) + (j * coeff_count));
        }
    }
}
void poc_encrypt_gsw_sk(vector<Ciphertext> &sk_gsw_ciphertext, shared_ptr<SEALContext> context, SecretKey sk, const int base_bit){
    const auto &context_data2 = context->first_context_data();
    auto &parms2 = context_data2->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms2.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);
    int logN = seal::util::get_power_of_two(coeff_count);
    const int lvl = context_data2->total_coeff_modulus_bit_count() / base_bit;

    Plaintext secret_key_pt;


    secret_key_pt.resize(coeff_count * coeff_modulus_size);
    for (int i = 0; i < coeff_modulus_size * coeff_count; i++) {
        secret_key_pt.data()[i] = sk.data().data()[i];

    }

    for (int i = 0; i < coeff_modulus_size; i++) {
        inverse_ntt_negacyclic_harvey(secret_key_pt.data() + i * coeff_count, small_ntt_tables[i]);
    }

    //intt version of gsw encryption of sk
    poc_gsw_enc128_sk(lvl, base_bit, context, sk, sk_gsw_ciphertext, secret_key_pt, pool);


}




void poc_multiply_add_plain_without_scaling_variant_sk(
        const Plaintext &plain,
        const SEALContext::ContextData &context_data,
        const int shift_amount,
        uint64_t *destination,
        seal::util::MemoryPool &pool
) {
    auto &parms = context_data.parms();
    size_t coeff_count = parms.poly_modulus_degree();
    size_t plain_coeff_count = plain.coeff_count() / 2;
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_mod_count = coeff_modulus.size();
    auto plain_modulus = context_data.parms().plain_modulus();
    auto coeff_div_plain_modulus = context_data.coeff_div_plain_modulus();
    uint64_t h;


    for (size_t i = 0; i < plain_coeff_count; i++) {

        // Add to ciphertext: h * m
        for (size_t j = 0; j < coeff_mod_count; j++) {

            //init empty 128 bit integers
            auto ptr(allocate_uint(coeff_mod_count, pool));
            auto ptr2(allocate_uint(coeff_mod_count, pool));

            //set 1 in lsb (it will be used for bit shifts)
            //ptr[0]=plain.data()[i];

            uint64_t plain_coeff;
            plain_coeff = plain.data()[i + (j * coeff_count)];
            ptr[0] = 1;
            ptr[1] = 0;
            ptr2[0] = 0;
            ptr2[1] = 0;
            //use 128 bit implementation for left shifts 1<<shiftamount
            util::left_shift_uint128(ptr.get(), shift_amount, ptr2.get());


//            cout<<"------------------"<<endl;
//            cout<< plain.data()[i + (j * coeff_count)]<<endl;
//            cout<< ptr2[0]<<ptr2[1]<<endl;
//            cout<< shift_amount<<endl;

            h = seal::util::barrett_reduce_128(ptr2.get(), coeff_modulus[j]);

            h = seal::util::multiply_uint_uint_mod(h, plain_coeff, coeff_modulus[j]);

            destination[i + (j * coeff_count)] = seal::util::add_uint_uint_mod(
                    destination[i + (j * coeff_count)], h, coeff_modulus[j]);


        }


    }


}




void poc_nfllib_ntt_ct(Ciphertext &ct , shared_ptr<SEALContext> &context){
    const auto &context_data = context->get_context_data(context->first_parms_id());
    auto &parms2 = context_data->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_count = parms2.poly_modulus_degree();
    size_t coeff_mod_count = coeff_modulus.size();
    size_t encrypted1_size = 2;
    for (size_t j = 0; j < encrypted1_size; j++) {
        uint64_t *ct_ptr = ct.data(j);

        for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
            for (size_t i = 0; i < poly_t::degree; i++) {
                resa[0](cm, i) = ct_ptr[(cm * poly_t::degree) + i];
            }
        }


       resa[0].ntt_pow_phi();



        for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
            for (size_t i = 0; i < poly_t::degree; i++) {
                ct_ptr[cm * poly_t::degree + i] = resa[0](cm, i);
            }
        }


    }
}

void poc_nfllib_ntt_rlwe_decomp(vector<uint64_t *> &rlwe_expansion){



    for (int k = 0; k < rlwe_expansion.size(); k++) {
        uint64_t *encrypted_rlwe_ptr = rlwe_expansion[k];

        for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
            for (size_t i = 0; i < poly_t::degree; i++) {
                resa[0](cm, i) = encrypted_rlwe_ptr[(cm * poly_t::degree) + i];
            }
        }

        resa[0].ntt_pow_phi();



        for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
            for (size_t i = 0; i < poly_t::degree; i++) {
                encrypted_rlwe_ptr[cm * poly_t::degree + i] = resa[0](cm, i);
            }
        }

    }



}

void poc_nfllib_intt_ct(Ciphertext &ct , shared_ptr<SEALContext> &context){
    const auto &context_data = context->get_context_data(context->first_parms_id());
    auto &parms2 = context_data->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_count = parms2.poly_modulus_degree();
    size_t coeff_mod_count = coeff_modulus.size();
    size_t encrypted1_size = 2;
    for (size_t j = 0; j < encrypted1_size; j++) {
        uint64_t *ct_ptr = ct.data(j);

        for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
            for (size_t i = 0; i < poly_t::degree; i++) {
                resa[0](cm, i) = ct_ptr[(cm * poly_t::degree) + i];
            }
        }


        resa[0].invntt_pow_invphi();



        for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
            for (size_t i = 0; i < poly_t::degree; i++) {
                ct_ptr[cm * poly_t::degree + i] = resa[0](cm, i);
            }
        }


    }
}

void poc_nfllib_plain_ct_prod(Ciphertext &ct , Plaintext &pt,
                                shared_ptr<SEALContext> &context, Ciphertext &res_ct){
   // auto start = std::chrono::steady_clock::now();
//return;
    const auto &context_data = context->get_context_data(context->first_parms_id());
    auto &parms2 = context_data->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_count = parms2.poly_modulus_degree();
    size_t coeff_mod_count = coeff_modulus.size();
    size_t encrypted1_size = 2;


//    poly_t *resa = alloc_aligned<poly_t, 32>(1),
//            *resb = alloc_aligned<poly_t, 32>(1),
//            *resc = alloc_aligned<poly_t, 32>(1);





    for (size_t j = 0; j < encrypted1_size; j++) {
        uint64_t *ct_ptr = ct.data(j);
        uint64_t *pt_ptr = pt.data();

        for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
            for (size_t i = 0; i < poly_t::degree; i++) {
                resa[0](cm, i) = ct_ptr[(cm * poly_t::degree) + i];
                resb[0](cm, i) = pt_ptr[(cm * poly_t::degree) + i];
            }
        }


        nfl::mul(resc[0], resa[0], resb[0]);



        for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
            for (size_t i = 0; i < poly_t::degree; i++) {
                ct_ptr[cm * poly_t::degree + i] = resc[0](cm, i);
            }
        }


    }

//    auto end = std::chrono::steady_clock::now();
//
//    int duration = duration_cast<std::chrono::microseconds>(start - end).count();
//    cout << "===================== " << duration <<"us" << endl;


//    free_aligned(1, resa);
//    free_aligned(1, resb);
//    free_aligned(1, resc);


    //end = std::chrono::steady_clock::now();





}

void poc_nfllib_external_product(vector<Ciphertext> &gsw_enc, vector<uint64_t *> &rlwe_expansion,
                                 shared_ptr<SEALContext> &context, int l, Ciphertext &res_ct, int is_reusable=1) {

    //assert(gsw_enc.size()==l);

    const auto &context_data = context->get_context_data(gsw_enc[0].parms_id());
    auto &parms2 = context_data->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_count = parms2.poly_modulus_degree();
    size_t coeff_mod_count = coeff_modulus.size();
    size_t encrypted1_size = 2;
    //auto small_ntt_tables = context_data->small_ntt_tables();

    //assert(gsw_enc.size() == encrypted1_size * l);
    //assert(rlwe_expansion.size() == encrypted1_size * l);



    std::uint64_t *result;
    Plaintext pt_tmp;

//    std::uint64_t *c00;
//    c00 = (std::uint64_t *) calloc((coeff_count), sizeof(uint64_t));
//    std::uint64_t *c01;
//    c01 = (std::uint64_t *) calloc((coeff_count), sizeof(uint64_t));
//    std::uint64_t *c10;
//    c10 = (std::uint64_t *) calloc((coeff_count), sizeof(uint64_t));
//    std::uint64_t *c11;
//    c11 = (std::uint64_t *) calloc((coeff_count), sizeof(uint64_t));
    //auto expand_start = std::chrono::high_resolution_clock::now();





    for (int k = 0; k < gsw_enc.size(); k++) {



        for (size_t j = 0; j < encrypted1_size; j++) {
            //j==0,j=1
            uint64_t *encrypted_gsw_ptr = gsw_enc[k].data(j);
            uint64_t *encrypted_rlwe_ptr = rlwe_expansion[k];
            result = (std::uint64_t *) calloc((coeff_count * coeff_mod_count), sizeof(uint64_t));

            auto start = std::chrono::steady_clock::now();
            poly_nfllib_prod_with_no_ntt(encrypted_gsw_ptr, encrypted_rlwe_ptr, result, coeff_count, coeff_mod_count);


            //poly_nfllib_mul(encrypted_gsw_ptr, encrypted_rlwe_ptr, result, coeff_count, coeff_mod_count, is_reusable);
            //poly_nfllib_mul_preprocessed(encrypted_gsw_ptr, encrypted_rlwe_ptr, result, coeff_count, coeff_mod_count);



            poly_nfllib_add(result,res_ct.data(j),res_ct.data(j));

            auto end = std::chrono::steady_clock::now();
            int duration= duration_cast<std::chrono::microseconds >(end - start).count();
            //cout<<"external~"<< duration<<endl;

            //cout<<"myduration="<<duration<<endl;
//            for (size_t i = 0; i < coeff_mod_count; i++) {
//
////                seal::util::add_poly_poly_coeffmod(res_ct.data(j) + (i*coeff_count), result+ (i*coeff_count), coeff_count, coeff_modulus[i].value(), res_ct.data(j) + (i*coeff_count));
//
//                if (j == 0 && i == 0) {
//                    seal::util::add_poly_poly_coeffmod(c00, result, coeff_count, coeff_modulus[i].value(), c00);
//
//                } else if (j == 0 && i == 1) {
//                    seal::util::add_poly_poly_coeffmod(c01, result + coeff_count, coeff_count, coeff_modulus[i].value(),
//                                                       c01);
//                    //myadd_poly_poly_coeffmod(c01, result, coeff_count, coeff_modulus[i].value(), c01);
//                } else if (j == 1 && i == 0) {
//                    seal::util::add_poly_poly_coeffmod(c10, result, coeff_count, coeff_modulus[i].value(), c10);
//                    // myadd_poly_poly_coeffmod(c10, result, coeff_count, coeff_modulus[i].value(), c10);
//                } else if (j == 1 && i == 1) {
//                    seal::util::add_poly_poly_coeffmod(c11, result + coeff_count, coeff_count, coeff_modulus[i].value(),
//                                                       c11);
//                    //myadd_poly_poly_coeffmod(c11, result, coeff_count, coeff_modulus[i].value(), c11);
//                }
//
//
//            }


            free(result);


        }






    }





//    seal::util::set_poly_poly(c00, coeff_count, 1, res_ct.data(0));
//    seal::util::set_poly_poly(c01, coeff_count, 1, res_ct.data(0) + coeff_count);
//    seal::util::set_poly_poly(c10, coeff_count, 1, res_ct.data(1));
//    seal::util::set_poly_poly(c11, coeff_count, 1, res_ct.data(1) + coeff_count);

    //poc_nfllib_intt_ct(res_ct,context);

//    free(c00);
//    free(c01);
//    free(c10);
//    free(c11);

}



void poly_nfllib_add(std::uint64_t *p1, std::uint64_t *p2, std::uint64_t *res){
    for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
        for (size_t i = 0; i < poly_t::degree; i++) {
            resa[0](cm, i) = p1[(cm * poly_t::degree) + i];
            resb[0](cm, i) = p2[(cm * poly_t::degree) + i];
        }
    }
    nfl::add(resc[0], resa[0], resb[0]);
    for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
        for (size_t i = 0; i < poly_t::degree; i++) {
            res[cm * poly_t::degree + i] = resc[0](cm, i);
        }
    }

}
void poly_nfllib_mul(std::uint64_t *p1, std::uint64_t *p2, std::uint64_t *res, const size_t coeff_count,
                     const std::uint64_t coeff_mod_count, int is_reusable=1) {

    //ternary is_reusable = 1 gsw is not reused in future, is_reusable=2 it will be reused in future so save p1 in ntt, is_reusable=3 it is reused in current mul so p1 is already in ntt

    //using poly_t = nfl::poly_from_modulus<uint32_t , 4096, 64>;

    //using poly_t = nfl::poly_from_modulus<uint64_t, 2048, 128>;
    //using poly_t = nfl::poly_from_modulus<uint32_t , 8192, 64>;
    //4096, 124, uint64_t





//    if(is_reusable>1)
//    cout<<"is_reusable="<< is_reusable<<endl;
    //start = std::chrono::steady_clock::now();
//    poly_t *resaa = alloc_aligned<poly_t, 32>(1),
//            *resbb = alloc_aligned<poly_t, 32>(1),
//            *rescc = alloc_aligned<poly_t, 32>(1);

    //end = std::chrono::steady_clock::now();



    for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
        for (size_t i = 0; i < poly_t::degree; i++) {
            resa[0](cm, i) = p1[(cm * poly_t::degree) + i];
            resb[0](cm, i) = p2[(cm * poly_t::degree) + i];
        }
    }


//    std::fill(resa, resa + 1, 123234);
//    std::fill(resb, resb + 1, 123244);



    if(is_reusable!=3) {
        resa[0].ntt_pow_phi();
    }
    resb[0].ntt_pow_phi();

    auto start = std::chrono::steady_clock::now();

    nfl::mul(resc[0], resa[0], resb[0]);
    auto end = std::chrono::steady_clock::now();

    int duration = duration_cast<std::chrono::microseconds>(start - end).count();
    //cout << "===================== " << duration <<"us" << endl;

    resc[0].invntt_pow_invphi();


    for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
        for (size_t i = 0; i < poly_t::degree; i++) {
            if(is_reusable==2){
                p1[cm * poly_t::degree + i] = resa[0](cm, i);
            }
            res[cm * poly_t::degree + i] = resc[0](cm, i);
        }
    }


//    free_aligned(1, resa);
//    free_aligned(1, resb);
//    free_aligned(1, resc);


    //std::cout << "Time per polynomial NTT: " << get_time_us(start, end, 1) << " us" << std::endl;

}


void poly_nfllib_prod_with_no_ntt(std::uint64_t *p1, std::uint64_t *p2, std::uint64_t *res, const size_t coeff_count,
                     const std::uint64_t coeff_mod_count) {

    //ternary is_reusable = 1 gsw is not reused in future, is_reusable=2 it will be reused in future so save p1 in ntt, is_reusable=3 it is reused in current mul so p1 is already in ntt

    //using poly_t = nfl::poly_from_modulus<uint32_t , 4096, 64>;

    //using poly_t = nfl::poly_from_modulus<uint64_t, 2048, 128>;
    //using poly_t = nfl::poly_from_modulus<uint32_t , 8192, 64>;
    //4096, 124, uint64_t





//    if(is_reusable>1)
//    cout<<"is_reusable="<< is_reusable<<endl;
    //start = std::chrono::steady_clock::now();
//    poly_t *resaa = alloc_aligned<poly_t, 32>(1),
//            *resbb = alloc_aligned<poly_t, 32>(1),
//            *rescc = alloc_aligned<poly_t, 32>(1);

    //end = std::chrono::steady_clock::now();



    for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
        for (size_t i = 0; i < poly_t::degree; i++) {
            resa[0](cm, i) = p1[(cm * poly_t::degree) + i];
            resb[0](cm, i) = p2[(cm * poly_t::degree) + i];
        }
    }


//    std::fill(resa, resa + 1, 123234);
//    std::fill(resb, resb + 1, 123244);


    auto start = std::chrono::steady_clock::now();
    nfl::mul(resc[0], resa[0], resb[0]);
    auto end = std::chrono::steady_clock::now();

    //resc[0].invntt_pow_invphi();




    for (size_t cm = 0; cm < poly_t::nmoduli; cm++) {
        for (size_t i = 0; i < poly_t::degree; i++) {

            res[cm * poly_t::degree + i] = resc[0](cm, i);
        }
    }


//    free_aligned(1, resa);
//    free_aligned(1, resb);
//    free_aligned(1, resc);


    //std::cout << "Time per polynomial NTT: " << duration_cast<microseconds>(start - end).count() << " us" << std::endl;

}

void poc_nfllib_ntt_gsw(vector<Ciphertext> &gsw_enc, shared_ptr<SEALContext> &context) {


    for (int k = 0; k < gsw_enc.size(); k++) {

//        auto start = std::chrono::steady_clock::now();
        poc_nfllib_ntt_ct(gsw_enc[k],context);
//        auto end = std::chrono::steady_clock::now();
//        cout<< duration_cast<std::chrono::microseconds>(end - start).count()<<endl;
    }

}

void poc_nfllib_add_ct(Ciphertext &ct1, Ciphertext &ct2, shared_ptr<SEALContext> &context) {
    const auto &context_data = context->get_context_data(context->first_parms_id());
    auto &parms2 = context_data->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_count = parms2.poly_modulus_degree();
    size_t coeff_mod_count = coeff_modulus.size();
    size_t encrypted1_size = ct1.size();

    for (int j=0;j<encrypted1_size;j++) {
        uint64_t *encrypted_gsw_ptr1 = ct1.data(j);
        uint64_t *encrypted_gsw_ptr2 = ct2.data(j);
        poly_nfllib_add(encrypted_gsw_ptr1, encrypted_gsw_ptr2, encrypted_gsw_ptr1);
    }

}

void poc_expand_flat_threaded(vector<GSWCiphertext>::iterator &result, vector<Ciphertext> &packed_swap_bits,
                              shared_ptr<SEALContext> context, int size, seal::GaloisKeys &galkey) {

    const auto &context_data2 = context->first_context_data();
    auto &parms2 = context_data2->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms2.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);
    int logN = seal::util::get_power_of_two(coeff_count);

    vector<vector<Ciphertext>> threaded_expanded_ciphers;

    std::vector<std::thread> threads;

    vector<Ciphertext> expanded_t1(size);
    vector<Ciphertext> expanded_t2(size);
    vector<Ciphertext> expanded_t3(size);
    //outloop is from 0-to-(l-1)
    threaded_expanded_ciphers.emplace_back(expanded_t1);
    threaded_expanded_ciphers.emplace_back(expanded_t2);
    threaded_expanded_ciphers.emplace_back(expanded_t3);

    threads.emplace_back(std::thread(poc_rlwe_expand_threaded, packed_swap_bits[0], context,galkey,size,std::ref(expanded_t1)));
    threads.emplace_back(std::thread(poc_rlwe_expand_threaded, packed_swap_bits[1], context,galkey,size,std::ref(expanded_t2)));
    threads.emplace_back(std::thread(poc_rlwe_expand_threaded, packed_swap_bits[2], context,galkey,size,std::ref(expanded_t3)));

        //expanded_ciphers = poc_rlwe_expand(packed_swap_bits[i], context, galkey, size);

        //cout <<"---rlwe---"<< duration_cast<std::chrono::milliseconds >(rlwe_end - rlwe_start).count() <<"ms"<< endl;

    int k=0;
    for(auto& th: threads){
        th.join();
        vector<uint64_t *> rlwe_decom;
        for (int j = 0; j < size; j++) {
            ///put jth expanded ct in ith idx slot  of jt gswct
            result[j][k] = threaded_expanded_ciphers[k][j];
        }
    }

}

void poc_half_gsw_enc128(const uint64_t l, const uint64_t base_bit, shared_ptr<SEALContext> context, const SecretKey sk,
                         vector<Ciphertext> &gsw_ciphertext, Plaintext gsw_plain, seal::util::MemoryPool &pool,
                         uint64_t inv) {

    Encryptor encryptor(context, sk);
    Decryptor decryptor(context, sk);
    const auto &context_data = context->first_context_data();
    auto &parms = context_data->parms();
    auto &coeff_modulus = parms.coeff_modulus();


    //size_t ct_poly_count = context_data->parms().coeff_modulus().size();/// find good way of getting it
    int total_bits;
    uint64_t r_l = l;
    total_bits = (context_data->total_coeff_modulus_bit_count());
    for (int p = 0; p < r_l; p++) {
        const int shift_amount = ((total_bits) - ((p + 1) * base_bit));
        Ciphertext res;
        encryptor.encrypt_zero_symmetric(res);
        //set_ciphertext(res, context);
        mymultiply_add_plain_without_scaling_variant(gsw_plain, *context->first_context_data(), shift_amount,
                                                     res.data(0), pool, inv);
        gsw_ciphertext.push_back(res);
    }


}

void poc_half_gsw_enc128_combined(const uint64_t l, const uint64_t base_bit, shared_ptr<SEALContext> context, const SecretKey sk,
                         vector<Ciphertext> &gsw_ciphertext, Plaintext gsw_plain, seal::util::MemoryPool &pool,
                         uint64_t total_expand_size, uint64_t gap, uint64_t dimension_size) {

    Encryptor encryptor(context, sk);
    Decryptor decryptor(context, sk);
    const auto &context_data = context->first_context_data();
    auto &parms = context_data->parms();
    auto &coeff_modulus = parms.coeff_modulus();

    Ciphertext res;
    encryptor.encrypt_zero_symmetric(res);
    //size_t ct_poly_count = context_data->parms().coeff_modulus().size();/// find good way of getting it
    int total_bits;
    uint64_t r_l = l;
    total_bits = (context_data->total_coeff_modulus_bit_count());
    for (int p = 0; p < r_l; p++) {
        const int shift_amount = ((total_bits) - ((p + 1) * base_bit));

        //set_ciphertext(res, context);
        mymultiply_add_plain_without_scaling_variant_combined(gsw_plain, *context->first_context_data(), shift_amount,
                                                     res.data(0), pool, total_expand_size,  gap,  dimension_size, p);

    }
    gsw_ciphertext.push_back(res);
}

void mymultiply_add_plain_without_scaling_variant_combined(const Plaintext &plain, const SEALContext::ContextData &context_data,
                                                  const int shift_amount, uint64_t *destination,
                                                  seal::util::MemoryPool &pool, uint64_t total_expand_size,
                                                  uint64_t gap, uint64_t dimension_size,
                                                  uint64_t curr_component) {
    auto &parms = context_data.parms();
    size_t coeff_count = parms.poly_modulus_degree();
    size_t plain_coeff_count = plain.coeff_count();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_mod_count = coeff_modulus.size();
    auto plain_modulus = context_data.parms().plain_modulus();
    auto coeff_div_plain_modulus = context_data.coeff_div_plain_modulus();
    uint64_t h;


    //cout<< shift_amount << endl;

    uint64_t  total_dim_with_gap = dimension_size*gap;

    for (size_t i = 0; i < total_dim_with_gap; i++) {
        // Add to ciphertext: h * m
        for (size_t j = 0; j < coeff_mod_count; j++) {
            //init empty 128 bit integers
            auto ptr(allocate_uint(coeff_mod_count, pool));
            auto ptr2(allocate_uint(coeff_mod_count, pool));
            auto ptr3(allocate_uint(coeff_mod_count, pool));
            //set 1 in lsb (it will be used for bit shifts)

            uint64_t poly_inv;
            uint64_t plain_coeff;
            if (total_expand_size > 0) {
                seal::util::try_invert_uint_mod(total_expand_size, coeff_modulus[j], poly_inv);
                plain_coeff = seal::util::multiply_uint_uint_mod(plain.data()[i], poly_inv, coeff_modulus[j]);

            } else {
                plain_coeff = plain.data()[i];
            }


            ptr2[0] = 0;
            ptr2[1] = 0;
            ptr[0] = 1;
            ptr[1] = 0;
            //use 128 bit implementation for left shifts 1<<shiftamount
            util::left_shift_uint128(ptr.get(), shift_amount, ptr2.get());
            h = seal::util::barrett_reduce_128(ptr2.get(), coeff_modulus[j]);

            //barret reduction is used for converting 128 bit interger to mod q1, q2 where q1, q2 are max 64 bits

            h = seal::util::multiply_uint_uint_mod(h, plain_coeff, coeff_modulus[j]);

            uint64_t index= i+ curr_component*(total_dim_with_gap) + (j * coeff_count);
            destination[index] = seal::util::add_uint_uint_mod(
                    destination[index], h, coeff_modulus[j]);
        }
    }

}



void mymultiply_add_plain_without_scaling_variant(const Plaintext &plain, const SEALContext::ContextData &context_data,
                                                  const int shift_amount, uint64_t *destination,
                                                  seal::util::MemoryPool &pool, uint64_t inv) {
    auto &parms = context_data.parms();
    size_t coeff_count = parms.poly_modulus_degree();
    size_t plain_coeff_count = plain.coeff_count();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_mod_count = coeff_modulus.size();
    auto plain_modulus = context_data.parms().plain_modulus();
    auto coeff_div_plain_modulus = context_data.coeff_div_plain_modulus();
    uint64_t h;


    //cout<< shift_amount << endl;

    for (size_t i = 0; i < plain_coeff_count; i++) {
        // Add to ciphertext: h * m
        for (size_t j = 0; j < coeff_mod_count; j++) {
            //init empty 128 bit integers
            auto ptr(allocate_uint(coeff_mod_count, pool));
            auto ptr2(allocate_uint(coeff_mod_count, pool));
            auto ptr3(allocate_uint(coeff_mod_count, pool));
            //set 1 in lsb (it will be used for bit shifts)

            uint64_t poly_inv;
            uint64_t plain_coeff;
            if (inv > 0) {
                seal::util::try_invert_uint_mod(inv, coeff_modulus[j], poly_inv);
                plain_coeff = seal::util::multiply_uint_uint_mod(plain.data()[i], poly_inv, coeff_modulus[j]);

            } else {
                plain_coeff = plain.data()[i];
            }


            ptr2[0] = 0;
            ptr2[1] = 0;
            ptr[0] = 1;
            ptr[1] = 0;
            //use 128 bit implementation for left shifts 1<<shiftamount
            util::left_shift_uint128(ptr.get(), shift_amount, ptr2.get());
            h = seal::util::barrett_reduce_128(ptr2.get(), coeff_modulus[j]);

            //barret reduction is used for converting 128 bit interger to mod q1, q2 where q1, q2 are max 64 bits

            h = seal::util::multiply_uint_uint_mod(h, plain_coeff, coeff_modulus[j]);
            //cout<<h<<",";
            destination[i + (j * coeff_count)] = seal::util::add_uint_uint_mod(
                    destination[i + (j * coeff_count)], h, coeff_modulus[j]);
        }
    }

}

void poc_enc_sk_gsw(SecretKey sk, shared_ptr<SEALContext> context,  const int base_bit, vector<Ciphertext>& sk_gsw_ciphertext){

    const auto &context_data2 = context->first_context_data();
    auto &parms2 = context_data2->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms2.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);
    const int l = context_data2->total_coeff_modulus_bit_count() / base_bit;


    int logN = seal::util::get_power_of_two(coeff_count);

    vector<Ciphertext> expanded_ciphers(coeff_count);

    //we need to get secret key data --> inverse NTT
    Plaintext secret_key_pt;


    secret_key_pt.resize(coeff_count * coeff_modulus_size);
    for (int i = 0; i < coeff_modulus_size * coeff_count; i++) {
        secret_key_pt.data()[i] = sk.data().data()[i];

    }

    for (int i = 0; i < coeff_modulus_size; i++) {
        inverse_ntt_negacyclic_harvey(secret_key_pt.data() + i * coeff_count, small_ntt_tables[i]);
    }



    //get gsw encrypting seckey
    poc_gsw_enc128_sk(l, base_bit, context, sk, sk_gsw_ciphertext, secret_key_pt, pool);

    poc_nfllib_ntt_gsw(sk_gsw_ciphertext, context);
}

void thread_server_expand(vector<GSWCiphertext>::iterator &result, vector<Ciphertext> packed_swap_bits,
                          shared_ptr<SEALContext> context, int begin, int end, int size, seal::GaloisKeys galkey,
                         const int l, const int base_bit,const int lsk, const int bsk, vector<Ciphertext> sk_gsw_ciphertext2) {

    const auto &context_data2 = context->first_context_data();
    auto &parms2 = context_data2->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms2.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);



    int logN = seal::util::get_power_of_two(coeff_count);

    vector<Ciphertext> expanded_ciphers(coeff_count);
    int duration=0;
    auto rlwe_start = std::chrono::high_resolution_clock::now();
    for (int i = begin; i < end; i++) {

        //expanded_ciphers = rlweExpand(packed_swap_bits[i], context, galkey, size);
        expanded_ciphers = poc_rlwe_expand(packed_swap_bits[i], context, galkey, size);


        vector<uint64_t *> rlwe_decom;

        for (int j = 0; j < size; j++) {

            result[j][i] = expanded_ciphers[j];
            Ciphertext res_ct;

            rlwe_decom.clear();
            rwle_decompositions(expanded_ciphers[j], context, lsk, bsk, rlwe_decom);


            poc_nfllib_ntt_rlwe_decomp(rlwe_decom);


            res_ct.resize(context, context->first_context_data()->parms_id(), 2);
            poc_nfllib_external_product(sk_gsw_ciphertext2, rlwe_decom, context, lsk, res_ct);
            for (auto p : rlwe_decom) {
                free(p);
            }

            poc_nfllib_intt_ct(res_ct,context);
            result[j][i + l] = res_ct;
        }



    }
    auto rlwe_end = std::chrono::high_resolution_clock::now();
    duration = duration+ duration_cast<std::chrono::milliseconds >(rlwe_end - rlwe_start).count();
    //cout <<"GSW exansion time="<< duration <<" ms"<< endl;



}


void gsw_server_expand_combined(vector<GSWCiphertext>::iterator &result, vector<Ciphertext> packed_swap_bits,
                          shared_ptr<SEALContext> context, int begin, int end, int total_dim_size, seal::GaloisKeys galkey,
                          const int l, const int base_bit,const int lsk, const int bsk, vector<Ciphertext> sk_gsw_ciphertext2,
                          uint64_t total_expand_size) {

    const auto &context_data2 = context->first_context_data();
    auto &parms2 = context_data2->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms2.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);



    int logN = seal::util::get_power_of_two(coeff_count);

    vector<Ciphertext> expanded_ciphers(coeff_count);
    int duration=0;

    expanded_ciphers = poc_rlwe_expand(packed_swap_bits[0], context, galkey, total_expand_size);

    auto rlwe_start = std::chrono::high_resolution_clock::now();
    for (int i = begin; i < end; i++) {

        //expanded_ciphers = rlweExpand(packed_swap_bits[i], context, galkey, size);
        //expanded_ciphers = poc_rlwe_expand(packed_swap_bits[i], context, galkey, size);


        vector<uint64_t *> rlwe_decom;

        for (int j = 0; j < total_dim_size; j++) {

            result[j][i] = expanded_ciphers[j + i*total_dim_size];
            Ciphertext res_ct;

            rlwe_decom.clear();
            rwle_decompositions(expanded_ciphers[j+ i*total_dim_size], context, lsk, bsk, rlwe_decom);


            poc_nfllib_ntt_rlwe_decomp(rlwe_decom);


            res_ct.resize(context, context->first_context_data()->parms_id(), 2);
            poc_nfllib_external_product(sk_gsw_ciphertext2, rlwe_decom, context, lsk, res_ct);
            for (auto p : rlwe_decom) {
                free(p);
            }

            poc_nfllib_intt_ct(res_ct,context);
            result[j][i + l] = res_ct;
        }



    }
    auto rlwe_end = std::chrono::high_resolution_clock::now();
    duration = duration+ duration_cast<std::chrono::milliseconds >(rlwe_end - rlwe_start).count();
   // cout <<"GSW exansion time="<< duration <<" ms"<< endl;



}

void my_poc_gsw_enc128_sk(const uint64_t l, const uint64_t base_bit, shared_ptr<SEALContext> context,
                       const SecretKey sk, vector<Ciphertext> &gsw_ciphertext, Plaintext gsw_plain,
                       seal::util::MemoryPool &pool) {
    Encryptor encryptor(context, sk);
    //Decryptor decryptor(context, sk);
    const auto &context_data = context->first_context_data();
    auto &parms = context_data->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t ct_poly_count = context_data->parms().coeff_modulus().size();/// find good way of getting it
    int total_bits;
    uint64_t r_l = l;
    Ciphertext t;

    total_bits = (context_data->total_coeff_modulus_bit_count());

    for (int j = 0; j < 2; j++) {// c0, c1
        for (int p = 0; p < r_l; p++) {
            const int shift_amount = ((total_bits) - ((p + 1) * base_bit));
            Ciphertext res;
            encryptor.encrypt_zero_symmetric(res);
            //set_ciphertext(res, context);
            mymultiply_add_plain_without_scaling_variant_sk(gsw_plain, *context->first_context_data(), shift_amount,
                                                            res.data(j), pool);
            gsw_ciphertext.push_back(res);
        }
    }
}


void poc_gsw_enc128_sk(const uint64_t l, const uint64_t base_bit, shared_ptr<SEALContext> context, const SecretKey sk,
                       vector<Ciphertext> &gsw_ciphertext, Plaintext gsw_plain, seal::util::MemoryPool &pool) {
    Encryptor encryptor(context, sk);
    //Decryptor decryptor(context, sk);
    const auto &context_data = context->first_context_data();
    auto &parms = context_data->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t ct_poly_count = context_data->parms().coeff_modulus().size();/// find good way of getting it
    int total_bits;
    uint64_t r_l = l;
    Ciphertext t;

    total_bits = (context_data->total_coeff_modulus_bit_count());

    for (int j = 0; j < 2; j++) {// c0, c1
        for (int p = 0; p < r_l; p++) {
            const int shift_amount = ((total_bits) - ((p + 1) * base_bit));
            Ciphertext res;
            encryptor.encrypt_zero_symmetric(res);
            //set_ciphertext(res, context);
            mymultiply_add_plain_without_scaling_variant_sk(gsw_plain, *context->first_context_data(), shift_amount,
                                                            res.data(j), pool);
            gsw_ciphertext.push_back(res);
        }
    }

}

void
mymultiply_add_plain_without_scaling_variant_sk(const Plaintext &plain, const SEALContext::ContextData &context_data,
                                                const int shift_amount, uint64_t *destination,
                                                seal::util::MemoryPool &pool) {

    auto &parms = context_data.parms();
    size_t coeff_count = parms.poly_modulus_degree();
    size_t plain_coeff_count = plain.coeff_count() / 2;
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_mod_count = coeff_modulus.size();
    auto plain_modulus = context_data.parms().plain_modulus();
    auto coeff_div_plain_modulus = context_data.coeff_div_plain_modulus();
    uint64_t h;


    for (size_t i = 0; i < plain_coeff_count; i++) {

        // Add to ciphertext: h * m
        for (size_t j = 0; j < coeff_mod_count; j++) {

            //init empty 128 bit integers
            auto ptr(allocate_uint(coeff_mod_count, pool));
            auto ptr2(allocate_uint(coeff_mod_count, pool));

            //set 1 in lsb (it will be used for bit shifts)
            //ptr[0]=plain.data()[i];

            uint64_t plain_coeff;
            plain_coeff = plain.data()[i + (j * coeff_count)];
            ptr[0] = 1;
            ptr[1] = 0;
            ptr2[0] = 0;
            ptr2[1] = 0;
            //use 128 bit implementation for left shifts 1<<shiftamount
            util::left_shift_uint128(ptr.get(), shift_amount, ptr2.get());


//            cout<<"------------------"<<endl;
//            cout<< plain.data()[i + (j * coeff_count)]<<endl;
//            cout<< ptr2[0]<<ptr2[1]<<endl;
//            cout<< shift_amount<<endl;
            //barret reduction is used for converting 128 bit interger to mod q1, q2 where q1, q2 are max 64 bits
            h = seal::util::barrett_reduce_128(ptr2.get(), coeff_modulus[j]);

            h = seal::util::multiply_uint_uint_mod(h, plain_coeff, coeff_modulus[j]);

            destination[i + (j * coeff_count)] = seal::util::add_uint_uint_mod(
                    destination[i + (j * coeff_count)], h, coeff_modulus[j]);


        }


    }

}

vector<Ciphertext>
rlweExpand(Ciphertext packedquery, shared_ptr<SEALContext> context, seal::GaloisKeys galkey, uint64_t size) {
    Evaluator evaluator1(context);

    const auto &context_data = context->first_context_data();

    auto &parms = context_data->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    int N2 = parms.poly_modulus_degree();


    Ciphertext tempctxt;
    Ciphertext tempctxt_rotated;
    Ciphertext tempctxt_shifted;


//    //std::vector<std::uint64_t> galois_elts;
//    for (int i = 0; i < ceil(log2(N2)); i++) {
//        galois_elts.push_back((N2 + exponentiate_uint64(2, i)) / exponentiate_uint64(2, i));
//    }

    vector<Ciphertext> temp;
    Ciphertext tmp;
    temp.push_back(packedquery);
    int numIters = ceil(log2(size));   // size is a ---rlwe---power of 2.


    if (numIters > ceil(log2(N2))) {
        throw logic_error("m > coeff_count is not allowed.");
    }

    int startIndex = static_cast<int>(log2(N2) - numIters);



    auto time_server_s = high_resolution_clock::now();
    for (long i = 0; i < numIters; i++) {
        vector<Ciphertext> newtemp(temp.size() << 1);
        int index = startIndex + i;
        int power = (N2 >> index) + 1;//k
        int ai = (1 << index);
        for (int j = 0; j < (1 << i); j++) {

            //check this power
            // tempctxt_rotated = subs(result[j])
            //evaluator1.apply_galois(temp[j], galois_elts[i], galkey, tempctxt_rotated);

            evaluator1.apply_galois(temp[j], power, galkey, tempctxt_rotated);

            // result[j+ 2**i] = result[j] - tempctxt_rotated;

            evaluator1.sub(temp[j], tempctxt_rotated, newtemp[j + (1 << i)]);

            // divide by x^ai = multiply by x^(2N - ai).

            multiply_power_of_X(newtemp[j + (1 << i)], tempctxt_shifted, (N2 << 1) - ai, context);

            newtemp[j + (1 << i)] = tempctxt_shifted;

            evaluator1.add(tempctxt_rotated, temp[j], newtemp[j]);

        }
        temp = newtemp;
    }

    auto time_server_e = high_resolution_clock::now();
    int dur =  duration_cast<microseconds>(time_server_e - time_server_s).count();

    cout<<"inner loop=========="<<dur<<endl;

    return temp;
}

void my_poc_nfllib_external_product(vector<Ciphertext> gsw_enc, vector<uint64_t *> rlwe_expansion,
                                 shared_ptr<SEALContext> context, int l, Ciphertext &res_ct, int is_reusable) {

    const auto &context_data = context->get_context_data(gsw_enc[0].parms_id());
    auto &parms2 = context_data->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_count = parms2.poly_modulus_degree();
    size_t coeff_mod_count = coeff_modulus.size();
    size_t encrypted1_size = 2;
    //auto small_ntt_tables = context_data->small_ntt_tables();

    //assert(gsw_enc.size() == encrypted1_size * l);
    //assert(rlwe_expansion.size() == encrypted1_size * l);


    int duration = 0;

    std::uint64_t *result;
    Plaintext pt_tmp;

    std::uint64_t *c00;
    c00 = (std::uint64_t *) calloc((coeff_count), sizeof(uint64_t));
    std::uint64_t *c01;
    c01 = (std::uint64_t *) calloc((coeff_count), sizeof(uint64_t));
    std::uint64_t *c10;
    c10 = (std::uint64_t *) calloc((coeff_count), sizeof(uint64_t));
    std::uint64_t *c11;
    c11 = (std::uint64_t *) calloc((coeff_count), sizeof(uint64_t));
    //auto expand_start = std::chrono::high_resolution_clock::now();

    //start = std::chrono::steady_clock::now();

    for (int k = 0; k < encrypted1_size * l; k++) {

        for (size_t j = 0; j < encrypted1_size; j++) {
            //j==0,j=1
            uint64_t *encrypted_gsw_ptr = gsw_enc[k].data(j);
            uint64_t *encrypted_rlwe_ptr = rlwe_expansion[k];
            result = (std::uint64_t *) calloc((coeff_count * coeff_mod_count), sizeof(uint64_t));


            poly_nfllib_mul(encrypted_gsw_ptr, encrypted_rlwe_ptr, result, coeff_count, coeff_mod_count, is_reusable);
            //poly_nfllib_mul_preprocessed(encrypted_gsw_ptr, encrypted_rlwe_ptr, result, coeff_count, coeff_mod_count);

            for (size_t i = 0; i < coeff_mod_count; i++) {

//                seal::util::add_poly_poly_coeffmod(res_ct.data(j) + (i*coeff_count), result+ (i*coeff_count), coeff_count, coeff_modulus[i].value(), res_ct.data(j) + (i*coeff_count));

                if (j == 0 && i == 0) {
                    seal::util::add_poly_poly_coeffmod(c00, result, coeff_count, coeff_modulus[i].value(), c00);

                } else if (j == 0 && i == 1) {
                    seal::util::add_poly_poly_coeffmod(c01, result + coeff_count, coeff_count, coeff_modulus[i].value(),
                                                       c01);
                    //myadd_poly_poly_coeffmod(c01, result, coeff_count, coeff_modulus[i].value(), c01);
                } else if (j == 1 && i == 0) {
                    seal::util::add_poly_poly_coeffmod(c10, result, coeff_count, coeff_modulus[i].value(), c10);
                    // myadd_poly_poly_coeffmod(c10, result, coeff_count, coeff_modulus[i].value(), c10);
                } else if (j == 1 && i == 1) {
                    seal::util::add_poly_poly_coeffmod(c11, result + coeff_count, coeff_count, coeff_modulus[i].value(),
                                                       c11);
                    //myadd_poly_poly_coeffmod(c11, result, coeff_count, coeff_modulus[i].value(), c11);
                }


            }

            free(result);


        }


    }


//    end = std::chrono::steady_clock::now();
//    duration= duration_cast<std::chrono::milliseconds >(end - start).count();
//    cout<<"myduration="<<duration<<endl;


    seal::util::set_poly_poly(c00, coeff_count, 1, res_ct.data(0));
    seal::util::set_poly_poly(c01, coeff_count, 1, res_ct.data(0) + coeff_count);
    seal::util::set_poly_poly(c10, coeff_count, 1, res_ct.data(1));
    seal::util::set_poly_poly(c11, coeff_count, 1, res_ct.data(1) + coeff_count);


    free(c00);
    free(c01);
    free(c10);
    free(c11);

}

void set_ciphertext(Ciphertext &ct, shared_ptr<SEALContext> context) {
    const auto &context_data = context->get_context_data(ct.parms_id());
    auto &parms2 = context_data->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_count = parms2.poly_modulus_degree();
    size_t coeff_mod_count = coeff_modulus.size();
    size_t encrypted1_size = ct.size();


    for (size_t j = 0; j < encrypted1_size; j++) {
        uint64_t *encrypted1_ptr = ct.data(j);
        for (size_t i = 0; i < coeff_mod_count; i++) {
            uint64_t *tmp = encrypted1_ptr + (i * coeff_count);
            for (size_t k = 0; k < coeff_count; k++) {
                tmp[k] = 0;
            }
            //if(j==b)tmp[0]=h;
        }

    }

}

void
my_rwle_decompositions(Ciphertext rlwe_ct_1, shared_ptr<SEALContext> context, const uint64_t l, const uint64_t base_bit,
                       vector<uint64_t *> &rlwe_decom) {
    const auto &context_data2 = context->first_context_data();
    auto &parms2 = context_data2->parms();
    auto &coeff_modulus = parms2.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms2.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);


    //compose ciphertext for all q_i's
    context_data2->rns_tool()->base_q()->compose_array(rlwe_ct_1.data(0), coeff_count, pool);
    context_data2->rns_tool()->base_q()->compose_array(rlwe_ct_1.data(1), coeff_count, pool);


    //128 bits decomp as given in external product
    my_poc_decomp_rlwe128(rlwe_ct_1, l, context, rlwe_decom, base_bit, pool);

    //auto rlwe_start = std::chrono::high_resolution_clock::now();
    int ssize = rlwe_decom.size();
    for (int i = 0; i < ssize; i++) {
        //rwle_crt_decompose and my_decompose_array does same thing but rwle_crt_decompose is slower
        //rwle_crt_decompose(rlwe_decom[i], context, pool);
        //cout<<i<<endl;
        my_decompose_array(rlwe_decom[i], coeff_count, coeff_modulus, coeff_modulus_size, pool);

    }
}

void my_poc_decomp_rlwe128(Ciphertext ct, const uint64_t l, shared_ptr<SEALContext> context,
                           vector<uint64_t *> &vec_ciphertexts, int base_bit, seal::util::MemoryPool &pool) {
    assert(vec_ciphertexts.size() == 0);
    const uint64_t base = UINT64_C(1) << base_bit;
    const uint64_t mask = base - 1;

    const auto &context_data = context->get_context_data(ct.parms_id());
    auto &parms = context_data->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_count = parms.poly_modulus_degree();
    size_t coeff_mod_count = coeff_modulus.size();
    size_t ct_poly_count = ct.size();

    uint64_t r_l = l;
    int total_bits;
    std::uint64_t *res;

    for (int j = 0; j < ct_poly_count; j++) {// c0, c1
        total_bits = (context_data->total_coeff_modulus_bit_count());
        uint64_t *encrypted1_ptr = ct.data(j);

        for (int p = 0; p < r_l; p++) {
            vector<uint64_t *> results;
            res = (std::uint64_t *) calloc((coeff_count * coeff_mod_count), sizeof(uint64_t));
            const int shift_amount = ((total_bits) - ((p + 1) * base_bit));
            for (size_t k = 0; k < coeff_mod_count * coeff_count; k = k + 2) {
                auto ptr(allocate_uint(2, pool));
                ptr[0] = 0;
                ptr[1] = 0;
                seal::util::right_shift_uint128(&encrypted1_ptr[k], shift_amount, ptr.get());
                uint64_t temp1 = ptr[0] & mask;
                res[k] = temp1;

            }
            //results.push_back(res);
            vec_ciphertexts.push_back(res);
        }

    }
}

void my_decompose_array(uint64_t *value, size_t count, std::vector<Modulus> coeff_modulus, size_t coeff_mod_count,
                        MemoryPoolHandle pool) {
    if (!value) {
        throw invalid_argument("value cannot be null");
    }
    if (!pool) {
        throw invalid_argument("pool is uninitialized");
    }

    if (coeff_mod_count > 1) {
        if (!product_fits_in(count, coeff_mod_count)) {
            throw logic_error("invalid parameters");
        }

        // Decompose an array of multi-precision integers into an array of arrays,
        // one per each base element
        auto value_copy(allocate_uint(count * coeff_mod_count, pool));

        auto temp_array(allocate_uint(count * coeff_mod_count, pool));

        // Merge the coefficients first
        for (size_t i = 0; i < count; i++) {
            for (size_t j = 0; j < coeff_mod_count; j++) {
                temp_array[j + (i * coeff_mod_count)] = value[j + (i * coeff_mod_count)];
            }
        }

        set_zero_uint(count * coeff_mod_count, value);

        for (size_t i = 0; i < count; i++) {
            //set_uint_uint(value, size_, value_copy.get());

            // Temporary space for 128-bit reductions
            for (size_t j = 0; j < coeff_mod_count; j++) {
                // Reduce in blocks
                uint64_t temp[2]{0, temp_array[(i * coeff_mod_count) + coeff_mod_count - 1]};
                for (size_t k = coeff_mod_count - 1; k--;) {
                    temp[0] = temp_array[(i * coeff_mod_count) + k];
                    temp[1] = barrett_reduce_128(temp, coeff_modulus[j]);
                }

                // Save the result modulo i-th base element
                //value[i] = temp[1];
                value[(j * count) + i] = temp[1];
            }
        }
    }
}








