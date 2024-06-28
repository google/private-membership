//
// Created by Haris Mughees on 4/22/21.
//

#include "pir_client.h"

void print(const vector<int>& v) {
  for (int n : v) {
    cout << n << " ";
  }
  cout << endl;
}

pir_client::pir_client(const EncryptionParameters &params,
                                const PirParams &pir_parms): params_(params){

    newcontext_ = SEALContext::Create(params_);

    pir_params_ = pir_parms;

    keygen_ = make_unique<KeyGenerator>(newcontext_);


    SecretKey secret_key = keygen_->secret_key();
    encryptor_ = make_unique<Encryptor>(newcontext_, secret_key);


    decryptor_ = make_unique<Decryptor>(newcontext_, secret_key);
    evaluator_ = make_unique<Evaluator>(newcontext_);

}

GaloisKeys pir_client::generate_galois_keys() {
    // Generate the Galois keys needed for coeff_select.
    //check https://github.com/microsoft/SealPIR/blob/dd06e2ff10d966177dc9892e92852e60add3f8b6/pir_client.cpp#L137
    //In the polynomial view (not batching), a Galois automorphism by a Galois element p changes
    //        Enc(plain(x)) to Enc(plain(x^p)).

    std::vector<uint32_t> galois_elts;
    const auto &context_data = newcontext_->first_context_data();
    auto &parms = context_data->parms();
    int N = parms.poly_modulus_degree();
    int logN = seal::util::get_power_of_two(N);
    for (int i = 0; i < logN; i++) {
        galois_elts.push_back((N + seal::util::exponentiate_uint64(2, i)) / seal::util::exponentiate_uint64(2, i));

    }

    return keygen_->galois_keys_local(galois_elts);
}

uint64_t pir_client::get_fv_index(uint64_t element_idx, uint64_t ele_size) {
    auto N = params_.poly_modulus_degree();
    auto logt = params_.plain_modulus().bit_count();
    auto ele_per_ptxt = elements_per_ptxt(logt, N, ele_size);
    return static_cast<uint64_t>(element_idx / ele_per_ptxt);
}

uint64_t pir_client::get_fv_offset(uint64_t element_idx, uint64_t ele_size) {
    uint32_t N = params_.poly_modulus_degree();
    auto logt = params_.plain_modulus().bit_count();

    uint64_t ele_per_ptxt = elements_per_ptxt(logt, N, ele_size);
    return element_idx % ele_per_ptxt;
}


void pir_client::compute_inverse_scales(){
    if (indices_.size() != pir_params_.nvec.size()){
        throw invalid_argument("size mismatch");
    }
    int logt = params_.plain_modulus().bit_count();

    uint64_t N = params_.poly_modulus_degree();
    uint64_t t = params_.plain_modulus().value();
    int logN = log2(N);
    int logm = logN;

    inverse_scales_.clear();

    for(int i = 0; i < pir_params_.nvec.size(); i++){
        uint64_t index_modN = indices_[i] % N;
        uint64_t numCtxt = ceil ( (pir_params_.nvec[i] + 0.0) / N);  // number of query ciphertexts.
        uint64_t batchId = indices_[i] / N;
        if (batchId == numCtxt - 1) {
            cout << "Client: adjusting the logm value..." << endl;
            logm = ceil(log2((pir_params_.nvec[i] % N)));
        }

        uint64_t inverse_scale;


        int quo = logm / logt;
        int mod = logm % logt;
        inverse_scale = pow(2, logt - mod);
        if ((quo +1) %2 != 0){
            inverse_scale =  params_.plain_modulus().value() - pow(2, logt - mod);
        }
        inverse_scales_.push_back(inverse_scale);
        if ( (inverse_scale << logm)  % t != 1){
            throw logic_error("something wrong");
        }
        cout << "Client: logm, inverse scale, t = " << logm << ", " << inverse_scale << ", " << t << endl;
    }
}
PirQuery pir_client::generate_query(uint64_t desiredIndex) {

    indices_ = compute_indices(desiredIndex, pir_params_.nvec);

    //compute_inverse_scales();
    GSWCiphertext packed_ct;
    uint64_t dimension_size=FIRST_DIM;

    int logsize;
    int gap ;

     int base_bits = pir_params_.plain_base;
     int decomp_size = params_.plain_modulus().bit_count() / base_bits;
    SecretKey sk= keygen_->secret_key();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);

    vector<vector<Ciphertext> > result(2);
    int N = params_.poly_modulus_degree();

    Plaintext pt(params_.poly_modulus_degree());

    //handling first dimension first
    for (uint32_t i = 0; i < 1; i++) {
        uint32_t num_ptxts = ceil( (pir_params_.nvec[i] + 0.0) / N);//mostly 1
        dimension_size= pir_params_.nvec[i];// first dimension is always 64
        logsize = ceil(log2(dimension_size));//6
        gap = ceil(N / (1 << logsize));//64
        // initialize result.
        cout << "Client: index " << i +1  <<  "/ " <<  indices_.size() << " = " << indices_[i] << endl;
        cout << "Client: number of ctxts needed for query = " << num_ptxts << endl;
        for (uint32_t j =0; j < num_ptxts; j++){
            pt.set_zero();
            if (indices_[i] > N*(j+1) || indices_[i] < N*j){
#ifdef DEBUG
                cout << "Client: coming here: so just encrypt zero." << endl;
#endif
                // just encrypt zero
            } else{
#ifdef DEBUG
                cout << "Client: encrypting a real thing " << endl;
#endif
                uint64_t real_index = indices_[i] - N*j;
                pt[real_index*gap] = 1;
            }
            //Ciphertext dest;

            packed_ct.clear();
            poc_plain_gsw_enc128(decomp_size, base_bits, newcontext_, sk, packed_ct, pt, pool, dimension_size);
            //encryptor_->encrypt(pt, dest);
            //dest.parms_id() = params_.parms_id();



            result[i]=(packed_ct);

//            Plaintext ppt;
//            decryptor_->decrypt(result[i][2],ppt);
//            cout<<ppt.to_string()<<endl;
        }
    }

    int previous_dimension_size=0;
    int new_dimension_size=0;
    vector<uint64_t> new_indices;
    //compressing all the remaining dimensions into one dimension of size equal to sum of remaining dimensions

    if(indices_.size()>1) {

        for (uint32_t i = 1; i < indices_.size(); i++) {

            dimension_size = pir_params_.nvec[i];
            uint64_t real_index = indices_[i]+new_dimension_size ;
            new_indices.push_back(real_index);
            //previous_dimension_size=dimension_size;
            new_dimension_size=new_dimension_size+dimension_size;
        }


    pt.set_zero();
    dimension_size= new_dimension_size;
    logsize = ceil(log2(dimension_size));
    gap = ceil(N / (1 << logsize));

    for (uint32_t i = 0; i < new_indices.size(); i++) {
        uint32_t num_ptxts = ceil((pir_params_.nvec[i] + 0.0) / N);

        // initialize result.
        cout << "Client: index " << i + 2 << "/ " << indices_.size() << " = " << indices_[i] << endl;
//        cout << "Client: number of ctxts needed for query = " << num_ptxts << endl;
        uint64_t real_index = new_indices[i];
        pt[real_index*gap] = 1;

    }

    vector<Ciphertext> half_gsw_ciphertext;
    //poc_l_pack_enc128(l, base_bit, context, sk, half_gsw_ciphertext, msg, decryptor1,  pool);

    base_bits = pir_params_.gsw_base;
    //decomp_size = newcontext_->first_context_data()->total_coeff_modulus_bit_count() / base_bits;
    decomp_size = pir_params_.gsw_decomp_size;
    poc_half_gsw_enc128(decomp_size, base_bits, newcontext_, sk, half_gsw_ciphertext, pt, pool, (1 << logsize));

    result[1]=(half_gsw_ciphertext);
//    Plaintext ppt;
//   decryptor_->decrypt(result[1][0],ppt);
//   cout<<ppt.to_string()<<endl;
    }

    return result;
}



PirQuery pir_client::generate_query_combined(uint64_t desiredIndex) {

    indices_ = compute_indices(desiredIndex, pir_params_.nvec);

    //compute_inverse_scales();
    GSWCiphertext packed_ct;
    uint64_t dimension_size=FIRST_DIM;
    int new_dimension_size=0;

    int logsize;
    int gap ;

    int base_bits = pir_params_.plain_base;
    int decomp_size = params_.plain_modulus().bit_count() / base_bits;
    SecretKey sk= keygen_->secret_key();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);

    vector<vector<Ciphertext> > result(2);
    int N = params_.poly_modulus_degree();//4096

    Plaintext pt(params_.poly_modulus_degree());

    //handling first dimension first
    for (uint32_t i = 0; i < 1; i++) {
        uint32_t num_ptxts = ceil( (pir_params_.nvec[i] + 0.0) / N);//mostly 1

        cout << "Client: index " << i+1  << "/ " << indices_.size() << " = " << indices_[i] << endl;

        dimension_size= 2*pir_params_.nvec[i];// first dimension is always 64
        logsize = ceil(log2(dimension_size));//8
        gap = ceil(N / (1 << logsize));//16
        // initialize result.


        for (uint32_t j =0; j < num_ptxts; j++){
            pt.set_zero();
            if (indices_[i] > N*(j+1) || indices_[i] < N*j){

            } else{
                uint64_t real_index = indices_[i] - N*j;// N*j =0 mostly
                pt[real_index*gap] = 1;
            }
            //Ciphertext dest;

            packed_ct.clear();
            poc_plain_gsw_enc128_combined(decomp_size, base_bits, newcontext_, sk, packed_ct, pt, pool, pir_params_.nvec[i], gap);




            result[i]=(packed_ct);

//            Plaintext ppt;
//            decryptor_->decrypt(result[i][2],ppt);
//            cout<<ppt.to_string()<<endl;
        }
    }

    int previous_dimension_size=0;
    new_dimension_size=0;
    vector<uint64_t> new_indices;
    //compressing all the remaining dimensions into one dimension of size equal to sum of remaining dimensions

    if(indices_.size()>1) {

        for (uint32_t i = 1; i < indices_.size(); i++) {

            cout << "Client: index " << i+1  << "/ " << indices_.size() << " = " << indices_[i] << endl;

            dimension_size = pir_params_.nvec[i];
            uint64_t real_index = indices_[i]+new_dimension_size ;
            new_indices.push_back(real_index);
            //previous_dimension_size=dimension_size;
            new_dimension_size=new_dimension_size+dimension_size;
        }


        pt.set_zero();
        dimension_size= new_dimension_size;
        logsize = ceil(log2(dimension_size*pir_params_.gsw_decomp_size));
        gap = ceil(N / (1 << logsize));

        for (uint32_t i = 0; i < new_indices.size(); i++) {
            uint32_t num_ptxts = ceil((pir_params_.nvec[i] + 0.0) / N);


            uint64_t real_index = new_indices[i];
            pt[real_index*gap] = 1;

        }

        vector<Ciphertext> half_gsw_ciphertext;

        base_bits = pir_params_.gsw_base;
        decomp_size = pir_params_.gsw_decomp_size;
        poc_half_gsw_enc128_combined(decomp_size, base_bits, newcontext_, sk, half_gsw_ciphertext, pt, pool, (1 << logsize), gap, dimension_size);

        result[1]=(half_gsw_ciphertext);

    }

    return result;
}

PirQuery pir_client::generate_query_combined2(const vector<int>& first_dim_indices, uint64_t index) {
    cout << "Bucket index: " << index << endl;

    auto p = compute_indices2(index, first_dim_indices, pir_params_.nvec);

    first_dim_indices_ = p.first;
    indices2_ = p.second;

    //compute_inverse_scales();
    GSWCiphertext packed_ct;
    uint64_t dimension_size=FIRST_DIM;
    int new_dimension_size=0;

    int logsize;
    int gap ;

    int base_bits = pir_params_.plain_base;
    int decomp_size = params_.plain_modulus().bit_count() / base_bits;
    SecretKey sk= keygen_->secret_key();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);

    vector<vector<Ciphertext> > result(2);
    int N = params_.poly_modulus_degree();//4096

    Plaintext pt(params_.poly_modulus_degree());

    cout << endl;

    //handling first dimension first
    for (uint32_t i = 0; i < 1; i++) {
        uint32_t num_ptxts = ceil( (pir_params_.nvec[i] + 0.0) / N);//mostly 1
        assert(num_ptxts == 1);

        cout << "Client: First dimension / Size: " << first_dim_indices.size() << " / ";
        print(first_dim_indices_);

        dimension_size= 2*pir_params_.nvec[i];// first dimension is always 64
        logsize = ceil(log2(dimension_size));//8
        gap = ceil(N / (1 << logsize));//16
        // initialize result.


        for (uint32_t j =0; j < num_ptxts; j++){
            pt.set_zero();
            for (int k : first_dim_indices_) {
              if (k > N * (j + 1) || k < N * j) {

              } else {
                uint64_t real_index = k - N * j;// N*j =0 mostly
                pt[real_index * gap] = 1;
              }
              //Ciphertext dest;

//            Plaintext ppt;
//            decryptor_->decrypt(result[i][2],ppt);
//            cout<<ppt.to_string()<<endl;
              //break;
            }
            packed_ct.clear();
            poc_plain_gsw_enc128_combined(decomp_size, base_bits, newcontext_, sk, packed_ct, pt, pool,
                                          pir_params_.nvec[i], gap);


            result[i] = (packed_ct);
        }
    }

    int previous_dimension_size=0;
    new_dimension_size=0;
    vector<uint64_t> new_indices;
    //compressing all the remaining dimensions into one dimension of size equal to sum of remaining dimensions

    if(indices2_.size()>0) {

        for (uint32_t i = 0; i < indices2_.size(); i++) {

            cout << "Client: index " << i+1  << "/ " << indices2_.size() << " = " << indices2_[i] << endl;

            dimension_size = pir_params_.nvec[i + 1];
            uint64_t real_index = indices2_[i]+new_dimension_size ;
            new_indices.push_back(real_index);
            //previous_dimension_size=dimension_size;
            new_dimension_size=new_dimension_size+dimension_size;
        }


        pt.set_zero();
        dimension_size= new_dimension_size;
        logsize = ceil(log2(dimension_size*pir_params_.gsw_decomp_size));
        gap = ceil(N / (1 << logsize));

        for (uint32_t i = 0; i < new_indices.size(); i++) {
            uint32_t num_ptxts = ceil((pir_params_.nvec[i] + 0.0) / N);


            uint64_t real_index = new_indices[i];
            pt[real_index*gap] = 1;

        }

        vector<Ciphertext> half_gsw_ciphertext;

        base_bits = pir_params_.gsw_base;
        decomp_size = pir_params_.gsw_decomp_size;
        poc_half_gsw_enc128_combined(decomp_size, base_bits, newcontext_, sk, half_gsw_ciphertext, pt, pool, (1 << logsize), gap, dimension_size);

        result[1]=(half_gsw_ciphertext);

    }

    return result;
}

PirQuery pir_client::generate_query_combined3(const vector<vector<int>>& indices) {

  /*auto p = compute_indices2(index, first_dim_indices, pir_params_.nvec);

  first_dim_indices_ = p.first;
  indices2_ = p.second;*/

  //compute_inverse_scales();
  GSWCiphertext packed_ct;
  uint64_t dimension_size = FIRST_DIM;
  int new_dimension_size = 0;

  int logsize;
  int gap;

  int base_bits = pir_params_.plain_base;
  int decomp_size = params_.plain_modulus().bit_count() / base_bits;
  SecretKey sk = keygen_->secret_key();
  auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);

  vector<vector<Ciphertext> > result(2);
  int N = params_.poly_modulus_degree();//4096

  Plaintext pt(params_.poly_modulus_degree());

  cout << endl;

  //handling first dimension first
  for (uint32_t i = 0; i < 1; i++) {
    uint32_t num_ptxts = ceil((pir_params_.nvec[i] + 0.0) / N);//mostly 1
    assert(num_ptxts == 1);

    auto first_dim_indices = indices[0];
    cout << "Client: First dimension / Size: " << first_dim_indices.size() << " / ";
    print(first_dim_indices);

    dimension_size = 2 * pir_params_.nvec[i];// first dimension is always 64
    logsize = ceil(log2(dimension_size));//8
    gap = ceil(N / (1 << logsize));//16
    // initialize result.


    for (uint32_t j = 0; j < num_ptxts; j++) {
      pt.set_zero();
      for (int k = 0; k < first_dim_indices.size(); ++k) {
        if (first_dim_indices[k] == 0) {
          continue;
        }
        if (k > N * (j + 1) || k < N * j) {

        } else {
          uint64_t real_index = k - N * j;// N*j =0 mostly
          pt[real_index * gap] = 1;
        }
        //Ciphertext dest;

//            Plaintext ppt;
//            decryptor_->decrypt(result[i][2],ppt);
//            cout<<ppt.to_string()<<endl;
        //break;
      }
      packed_ct.clear();
      poc_plain_gsw_enc128_combined(decomp_size, base_bits, newcontext_, sk, packed_ct, pt, pool,
                                    pir_params_.nvec[i], gap);


      result[i] = (packed_ct);
    }
  }

  int previous_dimension_size = 0;
  new_dimension_size = 0;
  vector<uint64_t> new_indices;
  //compressing all the remaining dimensions into one dimension of size equal to sum of remaining dimensions

  if (indices.size() > 1) {

    for (uint32_t i = 1; i < indices.size(); i++) {

      cout << "Client: index " << i << "/ " << indices[i].size() << " = ";
      for (int j = 0; j < indices[i].size(); j++) {
        cout << indices[i][j] << " ";
      }
      cout << endl;

      dimension_size = pir_params_.nvec[i];
      for (int j = 0; j < indices[i].size(); j++) {
        if (indices[i][j] == 1) {
          uint64_t real_index = j + new_dimension_size;
          new_indices.push_back(real_index);
          //break;
        }
      }
      //previous_dimension_size=dimension_size;
      new_dimension_size = new_dimension_size + dimension_size;
    }


    pt.set_zero();
    dimension_size = new_dimension_size;
    logsize = ceil(log2(dimension_size * pir_params_.gsw_decomp_size));
    gap = ceil(N / (1 << logsize));

    for (uint32_t i = 0; i < new_indices.size(); i++) {
      uint32_t num_ptxts = ceil((pir_params_.nvec[i] + 0.0) / N);


      uint64_t real_index = new_indices[i];
      pt[real_index * gap] = 1;

    }

    vector<Ciphertext> half_gsw_ciphertext;

    base_bits = pir_params_.gsw_base;
    decomp_size = pir_params_.gsw_decomp_size;
    poc_half_gsw_enc128_combined(decomp_size, base_bits, newcontext_, sk, half_gsw_ciphertext, pt, pool, (1 << logsize),
                                 gap, dimension_size);

    result[1] = (half_gsw_ciphertext);

  }

  return result;
}

void pir_client::decrypt_results(std::vector<seal::Ciphertext> reply) {
    for (int i=0; i< reply.size();i++){
        Plaintext ppt;
        decryptor_->decrypt(reply[i],ppt);
        cout<< ppt.to_string()<<endl;
        cout<< decryptor_->invariant_noise_budget(reply[i])<<endl;
    }

}

Plaintext pir_client::decrypt_result(std::vector<seal::Ciphertext> reply) {

        Plaintext ppt;
        decryptor_->decrypt(reply[0],ppt);

        return ppt;

}


SecretKey pir_client::get_decryptor() {
    return keygen_->secret_key();
}

void pir_client::test_query_expansion(PirQuery query, GaloisKeys galkey) {

    const auto &context_data2 = newcontext_->first_context_data();
    auto &parms = context_data2->parms();
    auto &coeff_modulus = parms.coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    auto small_ntt_tables = context_data2->small_ntt_tables();
    size_t coeff_count = parms.poly_modulus_degree();
    auto pool = MemoryManager::GetPool(mm_prof_opt::FORCE_NEW);

    Plaintext pt;
    pt.resize(coeff_count);
    pt.set_zero();
    pt[0]=1;

    Plaintext msg(coeff_count);
    msg.set_zero();
    msg[1]=1;


    int logsize = ceil(log2(64));
    int gap = ceil(coeff_count / (1 << logsize));
    const int base_bits = 20;

    const int decomp_size = parms.plain_modulus().bit_count() / base_bits;

    vector<GSWCiphertext> list_enc;
    list_enc.resize(64, GSWCiphertext(3));

    vector<GSWCiphertext>::iterator list_enc_ptr = list_enc.begin();

    auto gsw_enc_time_start = std::chrono::steady_clock::now();

    GSWCiphertext packed_ct;

    uint64_t dimension_size=64;

    poc_plain_gsw_enc128(decomp_size, base_bits, newcontext_, keygen_->secret_key(), packed_ct, msg, pool, dimension_size);

    //poc_expand_flat_threaded(list_enc_ptr, packed_ct, context, dimension_size, galois_keys);
    poc_expand_flat(list_enc_ptr, query[0], newcontext_, 64, galkey);
    auto gsw_enc_time_end = std::chrono::steady_clock::now();


    vector<uint64_t *> plain_decom;
    plain_decompositions(pt, newcontext_, decomp_size, base_bits, plain_decom);

    poc_nfllib_ntt_rlwe_decomp(plain_decom);

    for(int ii=0;ii<list_enc.size();ii++) {
        Plaintext ppt;
        decryptor_->decrypt(list_enc[ii][2], ppt);
        cout << ppt.to_string() << endl;

//        poc_nfllib_ntt_gsw(list_enc[ii], newcontext_);
//        Ciphertext res_ct;
//        res_ct.resize(newcontext_, newcontext_->first_context_data()->parms_id(), 2);
//        poc_nfllib_external_product(list_enc[0], plain_decom, newcontext_, decomp_size, res_ct, 1);
//        poc_nfllib_intt_ct(res_ct, newcontext_);



    }

}

GSWCiphertext pir_client::get_enc_sk() {
    vector<Ciphertext> sk_gsw_ciphertext;
    poc_enc_sk_gsw(keygen_->secret_key(),newcontext_,pir_params_.secret_base, sk_gsw_ciphertext);
    return sk_gsw_ciphertext;
}

PirQuery pir_client::generate_perm_query(std::uint64_t desiredIndex) {
    uint64_t swapbitsSize = 16;
    int logsize = ceil(log2(swapbitsSize));
    int gap = ceil(pir_params_.n / (1 << logsize));


}
