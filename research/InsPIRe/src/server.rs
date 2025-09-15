#[cfg(target_feature = "avx2")]
use std::arch::x86_64::*;
use std::cmp::min;
use std::{marker::PhantomData, ops::Range, time::Instant};
use std::collections::HashMap;

use log::debug;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

use spiral_rs::aligned_memory::AlignedMemory64;
use spiral_rs::{
    arith::*,
    client::*,
    params::*,
    poly::*,
};
use rayon::prelude::*;

use crate::convolution::naive_multiply_matrices;
use crate::measurement::{get_size_bytes_one, Measurement};

use super::{
    bits::*,
    client::*,
    convolution::{negacyclic_perm_u32, Convolution},
    kernel::*,
    lwe::*,
    matmul::matmul_vec_packed,
    modulus_switch::ModulusSwitch,
    packing::*,
    params::*,
    scheme::*,
    transpose::*,
    util::*,
};

pub fn generate_y_constants<'a>(
    params: &'a Params,
) -> (Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>) {
    let mut y_constants = Vec::new();
    let mut neg_y_constants = Vec::new();
    for num_cts_log2 in 1..params.poly_len_log2 + 1 {
        let num_cts = 1 << num_cts_log2;

        // Y = X^(poly_len / num_cts)
        let mut y_raw = PolyMatrixRaw::zero(params, 1, 1);
        y_raw.data[params.poly_len / num_cts] = 1;
        let y = y_raw.ntt();

        let mut neg_y_raw = PolyMatrixRaw::zero(params, 1, 1);
        neg_y_raw.data[params.poly_len / num_cts] = params.modulus - 1;
        let neg_y = neg_y_raw.ntt();

        y_constants.push(y);
        neg_y_constants.push(neg_y);
    }

    (y_constants, neg_y_constants)
}

/// Takes a matrix of u64s and returns a matrix of T's.
///
/// Input is row x cols u64's.
/// Output is out_rows x cols T's.
pub fn split_alloc(
    buf: &[u64],
    special_bit_offs: usize,
    rows: usize,
    crossing_point: usize,
    cols: usize,
    out_rows: usize,
    inp_mod_bits: usize,
    pt_bits: usize,
    // inp_transposed: bool,
) -> Vec<u16> {
    let mut out = vec![0u16; out_rows * cols];

    assert!(out_rows >= rows);
    let blowup_factor = inp_mod_bits as f64 / pt_bits as f64;
    assert!(
        out_rows as f64
            >= (crossing_point as f64 * blowup_factor).ceil()
                + ((rows - crossing_point) as f64 * blowup_factor).ceil()
    );
    assert!(inp_mod_bits >= pt_bits);

    for j in 0..cols {
        let mut bytes_tmp = vec![0u8; out_rows * inp_mod_bits / 8];

        // read this column
        let mut bit_offs = 0;
        for i in 0..crossing_point {
            let inp = buf[i * cols + j];
            // if i >= crossing_point {
            //     bit_offs = special_bit_offs;
            // }
            write_bits(&mut bytes_tmp, inp, bit_offs, inp_mod_bits);
            bit_offs += inp_mod_bits;
        }
        bit_offs = special_bit_offs;
        for i in crossing_point..rows {
            let inp = buf[i * cols + j];
            write_bits(&mut bytes_tmp, inp, bit_offs, inp_mod_bits);
            bit_offs += inp_mod_bits;
        }

        // now, 'stretch' the column vertically
        let mut bit_offs = 0;
        for i in 0..out_rows {
            let out_val = read_bits(&bytes_tmp, bit_offs, pt_bits);
            out[i * cols + j] = out_val as u16;
            bit_offs += pt_bits;
            if bit_offs >= out_rows * inp_mod_bits {
                break;
            }
        }
    }

    out
}

pub fn generate_fake_pack_pub_params<'a>(params: &'a Params) -> Vec<PolyMatrixNTT<'a>> {
    let pack_pub_params = raw_generate_expansion_params(
        &params,
        &PolyMatrixRaw::zero(&params, 1, 1),
        params.poly_len_log2,
        params.t_exp_left,
        &mut ChaCha20Rng::from_entropy(),
        &mut ChaCha20Rng::from_seed(STATIC_SEED_2),
    );
    pack_pub_params
}

pub type Precomp<'a> = Vec<(PolyMatrixNTT<'a>, Vec<PolyMatrixNTT<'a>>, Vec<Vec<usize>>)>;

#[derive(Clone)]
pub struct OfflinePrecomputedValues<'a> {
    // pub hint_0: Vec<u64>,
    pub hint_1: Vec<u64>,
    pub pseudorandom_query_1: Vec<PolyMatrixNTT<'a>>,
    pub y_constants: (Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
    pub smaller_server: Option<YServer<'a, u16>>,
    // pub vec_prepacked_lwe: Vec<Vec<Vec<PolyMatrixNTT<'a>>>>,
    pub fake_pack_pub_params: Vec<PolyMatrixNTT<'a>>,
    pub precomp: Precomp<'a>,
    pub precomp_inspir_vec_first_layer: Vec<PrecompInsPIR<'a>>,
    pub precomp_inspir_vec: Vec<PrecompInsPIR<'a>>,
    // pub w_all: PolyMatrixNTT<'a>,
    // pub w_bar_all: PolyMatrixNTT<'a>,
    // pub v_mask: PolyMatrixNTT<'a>,
    pub offline_packing_keys: OfflinePackingKeys<'a>,
}

pub struct Response {
    pub packed_mask_mod_switched: Vec<Vec<u8>>,
    pub packed_body_mod_switched: Vec<Vec<u8>>,
}

impl Response {
    pub fn get_size_bytes(&self) -> usize {
        get_size_bytes_one(&self.packed_mask_mod_switched) + get_size_bytes_one(&self.packed_body_mod_switched)
    }
}

#[derive(Clone)]
pub struct YServer<'a, T: Sync> {
    pub params: &'a Params,
    pub packing_params_set: HashMap<usize, PackParams<'a>>,
    pub half_packing_params_set: HashMap<usize, PackParams<'a>>,
    pub smaller_params: Params,
    pub db_buf_aligned: AlignedMemory64, // db_buf: Vec<u8>, // stored transposed
    pub phantom: PhantomData<T>,
    pub protocol_type: ProtocolType,
    pub second_level_packing_mask: PackingType,
    pub second_level_packing_body: PackingType,
}

impl<'a, T: Sync> YServer<'a, T>
where
    T: Sized + Copy + ToU64 + Default + std::marker::Sync,
    *const T: ToM512,
    u64: From<T>,
{
    pub fn new<'b, I>(
        params: &'a Params,
        mut db: I,
        protocol_type: ProtocolType,
        second_level_packing_mask: PackingType,
        second_level_packing_body: PackingType,
        inp_transposed: bool,
        // pad_rows: bool,
        gammas: Vec<usize>,
    ) -> Self
    where
        I: Iterator<Item = T>,
    {
        let bytes_per_pt_el = std::mem::size_of::<T>(); //1; //((lwe_params.pt_modulus as f64).log2() / 8.).ceil() as usize;
        debug!("bytes_per_pt_el: {}", bytes_per_pt_el);

        let db_rows = 1 << (params.db_dim_1 + params.poly_len_log2);
        let db_cols = match protocol_type {
            ProtocolType::SimplePIR => params.instances * params.poly_len,
            ProtocolType::DoublePIR => 1 << (params.db_dim_2 + params.poly_len_log2),
            ProtocolType::InsPIRe => params.instances * params.poly_len,
        };

        let sz_bytes = db_rows * db_cols * bytes_per_pt_el;

        let mut db_buf_aligned = AlignedMemory64::new(sz_bytes / 8);
        let db_buf_mut = as_bytes_mut(&mut db_buf_aligned);
        let db_buf_ptr = db_buf_mut.as_mut_ptr() as *mut T;

        for i in 0..db_rows {
            for j in 0..db_cols {
                let idx = if inp_transposed { i * db_cols + j } else { j * db_rows + i };

                unsafe {
                    *db_buf_ptr.add(idx) = db.next().unwrap();
                }
            }
        }

        // Parameters for the second round (the "DoublePIR" round)
        let smaller_params = match protocol_type {
            ProtocolType::SimplePIR => params.clone(),
            ProtocolType::DoublePIR => {
                let lwe_params = LWEParams::default();
                let pt_bits = (params.pt_modulus as f64).log2().floor() as usize;
                let blowup_factor = lwe_params.q2_bits as f64 / pt_bits as f64;
                let mut smaller_params = params.clone();
                smaller_params.db_dim_1 = params.db_dim_2;
                smaller_params.db_dim_2 = ((blowup_factor * (lwe_params.n + 1) as f64)
                    / params.poly_len as f64)
                    .log2()
                    .ceil() as usize;

                let out_rows = 1 << (smaller_params.db_dim_2 + params.poly_len_log2);
                assert_eq!(smaller_params.db_dim_1, params.db_dim_2);
                assert!(out_rows as f64 >= (blowup_factor * (lwe_params.n + 1) as f64));
                smaller_params
            }
            ProtocolType::InsPIRe => {
                let this_gamma = gammas[0];
                let pt_bits = (params.pt_modulus as f64).log2().floor() as usize;
                let blowup_factor = params.q2_bits as f64 / pt_bits as f64;
                let mut smaller_params = params.clone();
                smaller_params.db_dim_1 = (db_cols as f64 / this_gamma as f64).ceil().log2().ceil()
                    as usize
                    - params.poly_len_log2;
                // println!("smaller_params.db_dim_1:{}", smaller_params.db_dim_1);
                smaller_params.db_dim_2 = (((blowup_factor * params.poly_len as f64).ceil() + (blowup_factor * this_gamma as f64).ceil())
                    / params.poly_len as f64)
                    .log2()
                    .ceil() as usize;
                smaller_params
            }
        };

        let mut packing_params_set: HashMap<usize, PackParams> = HashMap::new(); 
        let mut half_packing_params_set: HashMap<usize, PackParams> = HashMap::new(); 
        for gamma in gammas.clone() {
            if !packing_params_set.contains_key(&gamma) {
                let packing_params = PackParams::new(&params, gamma);
                let half_packing_params = PackParams::new(&params, gamma >> 1);
                packing_params_set.insert(gamma, packing_params);
                half_packing_params_set.insert(gamma, half_packing_params);
            }

        }

        Self {
            params,
            packing_params_set: packing_params_set,
            half_packing_params_set: half_packing_params_set,
            smaller_params,
            db_buf_aligned,
            phantom: PhantomData,
            // pad_rows,
            protocol_type,
            second_level_packing_mask,
            second_level_packing_body,
        }
    }

    pub fn new_small<'b, I>(
        params: &'a Params,
        mut db: I,
        protocol_type: ProtocolType,
        inp_transposed: bool,
        // pad_rows: bool
    ) -> Self
    where
        I: Iterator<Item = T>,
    {
        let bytes_per_pt_el = std::mem::size_of::<T>(); //1; //((lwe_params.pt_modulus as f64).log2() / 8.).ceil() as usize;

        let db_rows = 1 << (params.db_dim_1 + params.poly_len_log2);
        let db_cols = match protocol_type {
            ProtocolType::SimplePIR => panic!("Shouldn't be here!"),
            ProtocolType::DoublePIR => 1 << (params.db_dim_2 + params.poly_len_log2),
            ProtocolType::InsPIRe => 1 << (params.db_dim_2 + params.poly_len_log2),
        };

        let sz_bytes = db_rows * db_cols * bytes_per_pt_el;

        let mut db_buf_aligned = AlignedMemory64::new(sz_bytes / 8);
        let db_buf_mut = as_bytes_mut(&mut db_buf_aligned);
        let db_buf_ptr = db_buf_mut.as_mut_ptr() as *mut T;

        for i in 0..db_rows {
            for j in 0..db_cols {
                let idx = if inp_transposed { i * db_cols + j } else { j * db_rows + i };

                unsafe {
                    *db_buf_ptr.add(idx) = db.next().unwrap();
                    // *db_buf_ptr.add(idx) = if i < db_rows {
                    //     db.next().unwrap()
                    // } else {
                    //     T::default()
                    // };
                }
            }
        }

        Self {
            params,
            packing_params_set: HashMap::new(),
            half_packing_params_set: HashMap::new(),
            smaller_params: params.clone(),
            db_buf_aligned,
            phantom: PhantomData,
            protocol_type,
            second_level_packing_mask: PackingType::NoPacking,
            second_level_packing_body: PackingType::NoPacking,
        }
    }

    pub fn db_rows(&self) -> usize {
        1 << (self.params.db_dim_1 + self.params.poly_len_log2)
    }

    pub fn db_cols(&self) -> usize {
        match self.protocol_type {
            ProtocolType::SimplePIR => self.params.instances * self.params.poly_len,
            ProtocolType::DoublePIR => 1 << (self.params.db_dim_2 + self.params.poly_len_log2),
            ProtocolType::InsPIRe => 1 << (self.params.db_dim_2 + self.params.poly_len_log2), /* self.params.instances * self.params.poly_len // TODO: Fix this!!! */
        }
    }

    pub fn get_smaller_params(&self) -> &Params {
        &self.smaller_params
    }

    pub fn multiply_batched_with_db_packed<const K: usize>(
        &self,
        aligned_query_packed: &[u64],
        query_rows: usize,
    ) -> AlignedMemory64 {
        // let db_rows = 1 << (self.params.db_dim_1 + self.params.poly_len_log2);
        let db_rows_padded = self.db_rows();
        let db_cols = self.db_cols();
        // assert_eq!(aligned_query_packed.len(), K * query_rows * db_rows_padded);
        assert_eq!(query_rows, 1);

        let now = Instant::now();
        let mut result = AlignedMemory64::new(db_cols);
        fast_batched_dot_product_generic::<_>(
            self.params,
            result.as_mut_slice(),
            &aligned_query_packed[..query_rows * db_rows_padded],
            db_rows_padded,
            &self.db(),
            db_rows_padded,
            db_cols,
        );
        debug!("Fast dot product in {} us", now.elapsed().as_micros());

        result
    }

    pub fn lwe_multiply_batched_with_db_packed<const K: usize>(
        &self,
        aligned_query_packed: &[u32],
    ) -> Vec<u32> {
        let _db_rows = 1 << (self.params.db_dim_1 + self.params.poly_len_log2);
        let db_cols = self.db_cols();
        let db_rows_padded = self.db_rows();
        assert_eq!(aligned_query_packed.len(), K * db_rows_padded);
        // assert_eq!(aligned_query_packed[db_rows + 1], 0);

        let mut result = vec![0u32; (db_cols + 8) * K];
        let now = Instant::now();
        
        let a_rows = db_cols;
        let a_true_cols = db_rows_padded;
        let a_cols = a_true_cols / 4; // order is inverted on purpose, because db is transposed
        let b_rows = a_true_cols;
        let b_cols = K;
        matmul_vec_packed(
            result.as_mut_slice(),
            self.db_u32(),
            aligned_query_packed,
            a_rows,
            a_cols,
            b_rows,
            b_cols,
        );
        let t = Instant::now();
        let result = transpose_generic(&result, db_cols, K);
        debug!("Transpose in {} us", t.elapsed().as_micros());
        debug!("Fast dot product in {} us", now.elapsed().as_micros());

        result
    }

    pub fn multiply_with_db_ring(
        &self,
        preprocessed_query: &[PolyMatrixNTT],
        col_range: Range<usize>,
        seed_idx: u8,
    ) -> Vec<u64> {
        // let db_rows_poly = 1 << (self.params.db_dim_1);
        // println!("{}", self.params.db_dim_1);
        let db_rows_poly =
            if (self.params.db_dim_1 as isize) < 0 { 1 } else { 1 << (self.params.db_dim_1) };

        let db_rows = 1 << (self.params.db_dim_1 + self.params.poly_len_log2);
        assert_eq!(preprocessed_query.len(), db_rows_poly);

        let mut result = Vec::new();
        let db = self.db();

        let mut prod = PolyMatrixNTT::zero(self.params, 1, 1);
        let mut db_elem_poly = PolyMatrixRaw::zero(self.params, 1, 1);
        let mut db_elem_ntt = PolyMatrixNTT::zero(self.params, 1, 1);

        for col in col_range.clone() {
            let mut sum = PolyMatrixNTT::zero(self.params, 1, 1);

            for row in 0..db_rows_poly {
                // let step = self.params.poly_len / 2;
                let how_many =
                    self.params.poly_len / (1 << (-min(0 as isize, self.params.db_dim_1 as isize)));

                for z in 0..how_many {
                    db_elem_poly.data[z] = db[col * db_rows + row * how_many + z].to_u64();
                }
                to_ntt(&mut db_elem_ntt, &db_elem_poly);

                multiply(&mut prod, &preprocessed_query[row], &db_elem_ntt);

                if row == db_rows_poly - 1 {
                    add_into(&mut sum, &prod);
                } else {
                    add_into_no_reduce(&mut sum, &prod);
                }
            }

            let sum_raw = sum.raw();

            // do negacyclic permutation (for first mul only)
            if seed_idx == SEED_0 && (self.protocol_type == ProtocolType::DoublePIR) {
                let sum_raw_transformed =
                    negacyclic_perm(sum_raw.get_poly(0, 0), 0, self.params.modulus);
                result.extend(&sum_raw_transformed);
            } else {
                result.extend(sum_raw.as_slice());
            }
        }

        // result
        let now = Instant::now();
        let res = transpose_generic(&result, col_range.len(), self.params.poly_len);
        debug!("transpose in {} us", now.elapsed().as_micros());
        res
    }

    pub fn generate_pseudorandom_query(&self, public_seed_idx: u8) -> Vec<PolyMatrixNTT<'a>> {
        let mut client = Client::init(&self.params);
        client.generate_secret_keys();
        let y_client = YClient::new(&mut client, &self.params);
        let query =
            y_client.generate_query_impl(public_seed_idx, self.params.db_dim_1, PackingType::InspiRING, 0); // TODO: is this correct type of packing?
        let query_mapped = query.iter().map(|x| x.submatrix(0, 0, 1, 1)).collect::<Vec<_>>();

        let mut preprocessed_query = Vec::new();
        for query_raw in query_mapped {
            // let query_raw_transformed =
            //     negacyclic_perm(query_raw.get_poly(0, 0), 0, self.params.modulus);
            // let query_raw_transformed = query_raw.get_poly(0, 0);
            let query_raw_transformed = if public_seed_idx == SEED_0 {
                negacyclic_perm(query_raw.get_poly(0, 0), 0, self.params.modulus)
                // query_raw.get_poly(0, 0).to_owned()
            } else {
                negacyclic_perm(query_raw.get_poly(0, 0), 0, self.params.modulus)
            };
            let mut query_transformed_pol = PolyMatrixRaw::zero(self.params, 1, 1);
            query_transformed_pol.as_mut_slice().copy_from_slice(&query_raw_transformed);
            preprocessed_query.push(query_transformed_pol.ntt());
        }

        preprocessed_query
    }

    pub fn answer_hint_ring(&self, public_seed_idx: u8, cols: usize) -> Vec<u64> {
        let preprocessed_query = self.generate_pseudorandom_query(public_seed_idx);
        let res = self.multiply_with_db_ring(&preprocessed_query, 0..cols, public_seed_idx);
        res
    }

    pub fn generate_hint_0(&self) -> Vec<u64> {
        let _db_rows = 1 << (self.params.db_dim_1 + self.params.poly_len_log2);
        let db_cols = self.db_cols();

        let mut rng_pub = ChaCha20Rng::from_seed(get_seed(SEED_0));
        let lwe_params = LWEParams::default();

        // pseudorandom LWE query is n x db_rows
        let psuedorandom_query =
            generate_matrix_ring(&mut rng_pub, lwe_params.n, lwe_params.n, db_cols);

        // db is db_cols x db_rows (!!!)
        // hint_0 is n x db_cols
        let hint_0 = naive_multiply_matrices(
            &psuedorandom_query,
            lwe_params.n,
            db_cols,
            &self.db(),
            self.db_rows(), // TODO: doesn't quite work
            db_cols,
            true,
        );
        hint_0.iter().map(|&x| x as u64).collect::<Vec<_>>()
    }

    pub fn generate_hint_0_ring(&self) -> Vec<u64> {
        let db_rows = 1 << (self.params.db_dim_1 + self.params.poly_len_log2);
        let db_cols = self.db_cols();

        let lwe_params = LWEParams::default();
        let n = lwe_params.n;
        let conv = Convolution::new(n);

        let mut hint_0 = vec![0u64; n * db_cols];

        let convd_len = conv.params().crt_count * conv.params().poly_len;

        let mut rng_pub = ChaCha20Rng::from_seed(get_seed(SEED_0));

        let mut v_nega_perm_a = Vec::new();
        for _ in 0..db_rows / n {
            let mut a = vec![0u32; n];
            for idx in 0..n {
                a[idx] = rng_pub.sample::<u32, _>(rand::distributions::Standard);
            }
            let nega_perm_a = negacyclic_perm_u32(&a);
            let nega_perm_a_ntt = conv.ntt(&nega_perm_a);
            v_nega_perm_a.push(nega_perm_a_ntt);
        }

        // limit on the number of times we can add results modulo M before we wrap
        let log2_conv_output =
            log2(lwe_params.modulus) + log2(lwe_params.n as u64) + log2(lwe_params.pt_modulus);
        let log2_modulus = log2(conv.params().modulus);
        let log2_max_adds = log2_modulus - log2_conv_output - 1;
        assert!(log2_max_adds > 0);
        let max_adds = 1 << log2_max_adds;

        for col in 0..db_cols {
            let mut tmp_col = vec![0u64; convd_len];
            for outer_row in 0..db_rows / n {
                let start_idx = col * self.db_rows() + outer_row * n;
                let pt_col = &self.db()[start_idx..start_idx + n];
                let pt_col_u32 = pt_col.iter().map(|&x| x.to_u64() as u32).collect::<Vec<_>>();
                assert_eq!(pt_col_u32.len(), n);
                let pt_ntt = conv.ntt(&pt_col_u32);

                let convolved_ntt = conv.pointwise_mul(&v_nega_perm_a[outer_row], &pt_ntt);

                for r in 0..convd_len {
                    tmp_col[r] += convolved_ntt[r] as u64;
                }

                if outer_row % max_adds == max_adds - 1 || outer_row == db_rows / n - 1 {
                    let mut col_poly_u32 = vec![0u32; convd_len];
                    for i in 0..conv.params().crt_count {
                        for j in 0..conv.params().poly_len {
                            let val = barrett_coeff_u64(
                                conv.params(),
                                tmp_col[i * conv.params().poly_len + j],
                                i,
                            );
                            col_poly_u32[i * conv.params().poly_len + j] = val as u32;
                        }
                    }
                    let col_poly_raw = conv.raw(&col_poly_u32);
                    for i in 0..n {
                        hint_0[i * db_cols + col] += col_poly_raw[i] as u64;
                        hint_0[i * db_cols + col] %= 1u64 << 32;
                    }
                    tmp_col.fill(0);
                }
            }
        }

        hint_0
    }

    pub fn answer_query(&self, aligned_query_packed: &[u64]) -> AlignedMemory64 {
        self.multiply_batched_with_db_packed::<1>(aligned_query_packed, 1)
    }

    pub fn answer_batched_queries<const K: usize>(
        &self,
        aligned_queries_packed: &[u64],
    ) -> AlignedMemory64 {
        self.multiply_batched_with_db_packed::<K>(aligned_queries_packed, 1)
    }

    pub fn perform_offline_precomputation_simplepir(
        &self,
        gamma: usize,
        measurement: Option<&mut Measurement>,
        online_only: bool
    ) -> OfflinePrecomputedValues {
        // Set up some parameters

        println!(" -- Offline 1");
        let params = self.params;
        assert_eq!(self.protocol_type, ProtocolType::SimplePIR);
        let db_cols = params.instances * params.poly_len;
        let db_cols_prime = (db_cols as f64 / gamma as f64).ceil() as usize;

        if online_only {
            assert_eq!(self.second_level_packing_mask, PackingType::InspiRING);
            let mut precomp_inspir_vec = Vec::with_capacity(db_cols_prime);
            let thing_size = if gamma == params.poly_len { params.poly_len / 2 } else {gamma};
            for _ in 0..db_cols_prime {
                precomp_inspir_vec.push(PrecompInsPIR {
                    a_hat: PolyMatrixRaw::random(&params, 1, 1),
                    bold_t_condensed: PolyMatrixNTT::random(&params, thing_size - 1, params.t_exp_left),
                    bold_t_bar_condensed: PolyMatrixNTT::random(&params, thing_size - 1, params.t_exp_left),
                    bold_t_hat_condensed: PolyMatrixNTT::random(&params, 1, params.t_exp_left)
                });
            }

            return  OfflinePrecomputedValues {
                hint_1: vec![],
                pseudorandom_query_1: vec![],
                y_constants: (vec![], vec![]),
                smaller_server: None,
                fake_pack_pub_params: vec![],
                precomp: vec![], //todo: put some garbage here
                precomp_inspir_vec_first_layer: vec![],
                precomp_inspir_vec,
                offline_packing_keys: OfflinePackingKeys::init_empty()
            }
        } 
        
        // Begin offline precomputation

        let now = Instant::now();
        let hint_0: Vec<u64> = self.answer_hint_ring(SEED_0, db_cols);
        println!(" -- Offline 2");

        // hint_0 is poly_len x db_cols
        if let Some(measurement) = measurement {
            measurement.offline.simplepir_prep_time_ms = now.elapsed().as_millis() as usize;
        }

        let now = Instant::now();

        let combined = [&hint_0[..], &vec![0u64; db_cols]].concat();
        assert_eq!(combined.len(), db_cols * (params.poly_len + 1)); // TODO: get rid of this +1 by replacing it in the underlying functions too

        let prepacked_lwe = prep_pack_many_lwes(&params, &combined, db_cols_prime, gamma);
        println!(" -- Offline 3");

        let offline_packing_keys = if self.second_level_packing_mask == PackingType::CDKS {
            OfflinePackingKeys::init_empty()
        } else {
            if gamma <= params.poly_len / 2 {
                OfflinePackingKeys::init(&self.packing_params_set[&gamma], W_SEED)
            } else {
                OfflinePackingKeys::init_full(&self.packing_params_set[&gamma], W_SEED, V_SEED)
            }
        };
        println!(" -- Offline 4");

        let mut y_constants: (Vec<PolyMatrixNTT<'_>>, Vec<PolyMatrixNTT<'_>>) =
            (Vec::new(), Vec::new());
        let mut fake_pack_pub_params: Vec<PolyMatrixNTT<'_>> = Vec::new();

        if self.second_level_packing_mask == PackingType::CDKS {
            y_constants = generate_y_constants(&params);
            fake_pack_pub_params = generate_fake_pack_pub_params(&params);
        }

        let mut precomp: Precomp = Vec::with_capacity(db_cols_prime);
        for _ in 0..db_cols_prime {
            precomp.push((PolyMatrixNTT::zero(&params, 1, 1), Vec::new(), Vec::new()))
        }

        let mut precomp_inspir_vec = Vec::with_capacity(db_cols_prime);
        for _ in 0..db_cols_prime {
            precomp_inspir_vec.push(PrecompInsPIR{
                a_hat: PolyMatrixRaw::zero(&params, 1, 1),
                bold_t_condensed: PolyMatrixNTT::zero(&params, 1, 1),
                bold_t_bar_condensed: PolyMatrixNTT::zero(&params, 1, 1),
                bold_t_hat_condensed: PolyMatrixNTT::zero(&params, 1, 1),
            })
        }

        match self.second_level_packing_mask {
            PackingType::InspiRING => {
                let results: Vec<_> = (0..db_cols_prime)
                .into_par_iter()
                .map(|i| {
                    let mut a_ct_tilde = Vec::new();
                    for j in 0..gamma {
                        a_ct_tilde.push(prepacked_lwe[i][j].submatrix(0, 0, 1, 1))
                    }
            
                    if gamma <= params.poly_len / 2 {
                        let w_all = offline_packing_keys.w_all.as_ref().unwrap();
                        packing_with_preprocessing_offline(
                            &self.packing_params_set[&gamma],
                            &w_all,
                            &a_ct_tilde,
                        )
                    } else {
                        let w_all = offline_packing_keys.w_all.as_ref().unwrap();
                        let w_bar_all = offline_packing_keys.w_bar_all.as_ref().unwrap();
                        let v_mask = offline_packing_keys.v_mask.as_ref().unwrap();
            
                        full_packing_with_preprocessing_offline(
                            &self.packing_params_set[&gamma],
                            &w_all,
                            &w_bar_all,
                            &v_mask,
                            &a_ct_tilde,
                        )
                    }
                })
                .collect();
    
                precomp_inspir_vec = results;    
        
            },
            PackingType::CDKS => {
                for i in 0..db_cols_prime {
                    let tup = precompute_pack(
                        params,
                        params.poly_len_log2,
                        &prepacked_lwe[i],
                        &fake_pack_pub_params,
                        &y_constants,
                    );
                    precomp.push(tup);
                }
            }
            PackingType::NoPacking => {
                panic!("Shouldn't be here!");
            }
        }
        debug!("Precomp in {} us", now.elapsed().as_micros());
        println!(" -- Offline 5");

        OfflinePrecomputedValues {
            hint_1: vec![],
            pseudorandom_query_1: vec![],
            y_constants,
            smaller_server: None,
            fake_pack_pub_params,
            precomp,
            precomp_inspir_vec_first_layer: vec![],
            precomp_inspir_vec,
            offline_packing_keys,
        }
    }

    pub fn perform_offline_precomputation(
        &self,
        gamma: usize,
        measurement: Option<&mut Measurement>,
        online_only: bool,
    ) -> OfflinePrecomputedValues {
        // Set up some parameters

        let params = self.params;

        let lwe_params = LWEParams::default();
        let db_cols = 1 << (params.db_dim_2 + params.poly_len_log2);

        // LWE reduced moduli
        let lwe_q_prime_bits = lwe_params.q2_bits as usize;
        let lwe_q_prime = lwe_params.get_q_prime_2();

        // The number of bits represented by a plaintext RLWE coefficient
        let pt_bits = (params.pt_modulus as f64).log2().floor() as usize;
        // assert_eq!(pt_bits, 16);

        // The factor by which ciphertext values are bigger than plaintext values
        let blowup_factor = lwe_q_prime_bits as f64 / pt_bits as f64;
        debug!("blowup_factor: {}", blowup_factor);

        // The starting index of the final value (the '1' in lwe_params.n + 1)
        // This is rounded to start on a pt_bits boundary
        let special_offs =
            ((lwe_params.n * lwe_q_prime_bits) as f64 / pt_bits as f64).ceil() as usize;
        let special_bit_offs = special_offs * pt_bits;

        let out_rows = 1 << (self.smaller_params.db_dim_2 + params.poly_len_log2); // = (tau * (n+1) ) rounded up to multiple of d
        // let rho = 1 << self.smaller_params.db_dim_2; // number of rlwe cts
        assert_eq!(self.smaller_params.db_dim_1, params.db_dim_2);
        assert!(out_rows as f64 >= (blowup_factor * (lwe_params.n + 1) as f64));

        let num_rlwes = (out_rows as f64 / gamma as f64).ceil() as usize;

        debug!(
            "the first {} LWE output ciphertexts of the DoublePIR round (out of {} total) are query-indepednent",
            special_offs, out_rows
        );
        debug!(
            "the next {} LWE output ciphertexts are query-dependent",
            blowup_factor.ceil() as usize
        );
        debug!("the rest are zero");

        // Begin offline precomputation

        if online_only {
            let mut rng_pub = ChaCha20Rng::from_seed(get_seed(SEED_1));
            let mut smaller_db : Vec<u16> = Vec::with_capacity(db_cols * out_rows);
            for _ in 0..db_cols * out_rows {
                smaller_db.push((rng_pub.sample::<u64, _>(rand::distributions::Standard) % params.pt_modulus) as u16);
            }

            assert_eq!(smaller_db.len(), db_cols * out_rows);
            let smaller_server: YServer<u16> = YServer::<u16>::new_small(
                &self.smaller_params,
                smaller_db.into_iter(),
                ProtocolType::InsPIRe,
                true,
            );

            let mut hint_1 : Vec<u64> = Vec::with_capacity(params.poly_len * out_rows);
            for _ in 0..params.poly_len * out_rows {
                hint_1.push(rng_pub.sample::<u64, _>(rand::distributions::Standard) % params.modulus);
            }

            let vec_size = 1 << smaller_server.params.db_dim_1;
            let mut pseudorandom_query_1 = Vec::with_capacity(vec_size);
            for _ in 0..vec_size {
                pseudorandom_query_1.push(PolyMatrixNTT::random(&params, 1, 1))
            }
            let mut y_constants = (Vec::new(), Vec::new());
            let mut fake_pack_pub_params = Vec::new();
    
            let offline_packing_keys = if self.second_level_packing_mask == PackingType::CDKS {
                y_constants = generate_y_constants(&params);
                fake_pack_pub_params = generate_fake_pack_pub_params(&params);
                OfflinePackingKeys::init_empty()
            } else {
                if gamma <= params.poly_len / 2 {
                    OfflinePackingKeys::init(&self.packing_params_set[&gamma], W_SEED)
                } else {
                    OfflinePackingKeys::init_full(&self.packing_params_set[&gamma], W_SEED, V_SEED)
                }
            };

            let precomp: Precomp = Vec::new();    
            let mut precomp_inspir_vec = Vec::new();

            let thing_size = if gamma == params.poly_len { params.poly_len / 2 } else {gamma};
            assert_eq!(self.second_level_packing_mask, PackingType::InspiRING);
            for _ in 0..num_rlwes {
                if self.second_level_packing_mask == PackingType::InspiRING {
                    precomp_inspir_vec.push(PrecompInsPIR {
                        a_hat: PolyMatrixRaw::random(&params, 1, 1),
                        bold_t_condensed: PolyMatrixNTT::random(&params, thing_size - 1, params.t_exp_left),
                        bold_t_bar_condensed: PolyMatrixNTT::random(&params, thing_size - 1, params.t_exp_left),
                        bold_t_hat_condensed: PolyMatrixNTT::random(&params, 1, params.t_exp_left)
                    });
                } else if self.second_level_packing_mask == PackingType::CDKS {
                    println!("No online-only CDKS");
                }
            }            

            return OfflinePrecomputedValues {
                // hint_0,
                hint_1,
                pseudorandom_query_1,
                y_constants,
                smaller_server: Some(smaller_server),
                fake_pack_pub_params,
                precomp,
                precomp_inspir_vec_first_layer: vec![],
                precomp_inspir_vec,
                offline_packing_keys,
            }            
        }


        let now = Instant::now();
        let hint_0: Vec<u64> = self.generate_hint_0_ring();
        // hint_0 is n x db_cols
        if let Some(measurement) = measurement {
            measurement.offline.simplepir_prep_time_ms = now.elapsed().as_millis() as usize;
        }
        debug!("Answered hint (ring) in {} us", now.elapsed().as_micros());

        // compute (most of) the secondary hint
        let intermediate_cts = [&hint_0[..], &vec![0u64; db_cols]].concat();
        let intermediate_cts_rescaled = intermediate_cts
            .iter()
            .map(|x| rescale(*x, lwe_params.modulus, lwe_q_prime))
            .collect::<Vec<_>>();

        // split and do a second PIR over intermediate_cts
        // split into blowup_factor=q/p instances (so that all values are now mod p)
        // the second PIR is over a database of db_cols x (blowup_factor * (lwe_params.n
        // + 1)) values mod p

        // inp: (lwe_params.n + 1, db_cols)
        // out: (out_rows >= (lwe_params.n + 1) * blowup_factor, db_cols)
        //      we are 'stretching' the columns (and padding)

        debug!("Splitting intermediate cts...");

        let smaller_db: Vec<u16> = split_alloc(
            &intermediate_cts_rescaled,
            special_bit_offs,
            lwe_params.n + 1,
            lwe_params.n,
            db_cols,
            out_rows,
            lwe_q_prime_bits,
            pt_bits,
            // false,
        );
        assert_eq!(smaller_db.len(), db_cols * out_rows);

        debug!("Done splitting intermediate cts.");

        // This is the 'intermediate' db after the first pass of PIR and expansion
        let smaller_server: YServer<u16> = YServer::<u16>::new_small(
            &self.smaller_params,
            smaller_db.into_iter(),
            ProtocolType::DoublePIR,
            true,
        );
        debug!("gen'd smaller server.");

        let hint_1 = smaller_server.answer_hint_ring(
            SEED_1,
            1 << (smaller_server.params.db_dim_2 + smaller_server.params.poly_len_log2),
        );
        assert_eq!(hint_1.len(), params.poly_len * out_rows);

        let pseudorandom_query_1 = smaller_server.generate_pseudorandom_query(SEED_1);

        let combined = [&hint_1[..], &vec![0u64; out_rows]].concat();
        assert_eq!(combined.len(), out_rows * (params.poly_len + 1));
        let prepacked_lwe = prep_pack_many_lwes(&params, &combined, num_rlwes, gamma);

        let now = Instant::now();

        let mut y_constants = (Vec::new(), Vec::new());
        let mut fake_pack_pub_params = Vec::new();

        let offline_packing_keys = if self.second_level_packing_mask == PackingType::CDKS {
            y_constants = generate_y_constants(&params);
            fake_pack_pub_params = generate_fake_pack_pub_params(&params);
            OfflinePackingKeys::init_empty()
        } else {
            if gamma <= params.poly_len / 2 {
                OfflinePackingKeys::init(&self.packing_params_set[&gamma], W_SEED)
            } else {
                OfflinePackingKeys::init_full(&self.packing_params_set[&gamma], W_SEED, V_SEED)
            }
        };

        let mut precomp: Precomp = Vec::new();
        let mut precomp_inspir_vec = Vec::new();
        for i in 0..num_rlwes {
            match self.second_level_packing_mask {
                PackingType::InspiRING => {
                    // Our version
                    let mut a_ct_tilde = Vec::with_capacity(prepacked_lwe[i].len());
                    for j in 0..gamma {
                        a_ct_tilde.push(prepacked_lwe[i][j].submatrix(0, 0, 1, 1))
                    }

                    let tup_inspir = if gamma <= params.poly_len/2 {
                        let w_all = offline_packing_keys.w_all.as_ref().unwrap();
                        packing_with_preprocessing_offline_without_rotations(
                            &self.packing_params_set[&gamma],
                            &w_all,
                            &a_ct_tilde,
                        )
                    } else {
                        let w_all = offline_packing_keys.w_all.as_ref().unwrap();
                        let w_bar_all = offline_packing_keys.w_bar_all.as_ref().unwrap();
                        let v_mask = offline_packing_keys.v_mask.as_ref().unwrap();

                        full_packing_with_preprocessing_offline_without_rotations(
                            &self.packing_params_set[&gamma],
                            &w_all,
                            &w_bar_all,
                            &v_mask,
                            &a_ct_tilde,
                        )
                    };
                    precomp_inspir_vec.push(tup_inspir)
                }
                PackingType::CDKS => {
                    let tup = precompute_pack(
                        params,
                        params.poly_len_log2,
                        &prepacked_lwe[i],
                        &fake_pack_pub_params,
                        &y_constants,
                    );
                    precomp.push(tup);
                }
                PackingType::NoPacking => {
                    panic!("Shouldn't be here!");
                }
            }
        }
        debug!("Precomp in {} us", now.elapsed().as_micros());

        OfflinePrecomputedValues {
            // hint_0,
            hint_1,
            pseudorandom_query_1,
            y_constants,
            smaller_server: Some(smaller_server),
            fake_pack_pub_params,
            precomp,
            precomp_inspir_vec_first_layer: vec![],
            precomp_inspir_vec,
            offline_packing_keys,
        }
    }

    pub fn perform_offline_precomputation_medium_payload(
        &self,
        gammas: &Vec<usize>,
        // gamma[0] is used to pack the first layer
        // gamma[1] is used to pack the mask in the second layer
        // gamma[2] is used to pack the body in the second layer
        measurement: Option<&mut Measurement>,
        online_only: bool,
    ) -> OfflinePrecomputedValues {
        // Set up some parameters

        assert_eq!(self.second_level_packing_mask, PackingType::InspiRING);
        assert!(gammas.len() >= 2);

        let params = self.params;
        let db_cols = params.instances * params.poly_len;
        let db_cols_prime = (db_cols as f64 / gammas[0] as f64).ceil() as usize;
        let rlwe_q_prime_2 = params.get_q_prime_2();
        let q2_bits = (rlwe_q_prime_2 as f64).log2().ceil() as usize;
        let pt_bits = (params.pt_modulus as f64).log2().floor() as usize;
        let blowup_factor = q2_bits as f64 / pt_bits as f64;
        let mask_decomposed_num_coeffs = (params.poly_len as f64 * blowup_factor).ceil() as usize;
        // let body_decomposed_num_coeffs = (gamma as f64 * blowup_factor).ceil() as usize;
        // let smaller_db_total_num_cols = mask_decomposed_num_coeffs + body_decomposed_num_coeffs;
        let smaller_db_total_num_cols_rounded_up =
            1 << (self.smaller_params.db_dim_2 + params.poly_len_log2); // out_columns after first multiplication and splitting
        // assert_eq!(mask_decomposed_num_coeffs % gamma, 0); -> still necessary?
        let mask_num_packed_rlwes = (mask_decomposed_num_coeffs as f64 / gammas[1] as f64).ceil() as usize;

        if online_only {
            let mut rng_pub = ChaCha20Rng::from_seed(get_seed(SEED_0));

            let mut hint_1 : Vec<u64> = Vec::with_capacity(params.poly_len * smaller_db_total_num_cols_rounded_up);
            for _ in 0..params.poly_len * smaller_db_total_num_cols_rounded_up {
                hint_1.push(rng_pub.sample::<u64, _>(rand::distributions::Standard) % params.modulus);
            }

            let mut rng_pub = ChaCha20Rng::from_seed(get_seed(SEED_1));
            let mut smaller_db : Vec<u16> = Vec::with_capacity(smaller_db_total_num_cols_rounded_up * db_cols_prime);
            for _ in 0..smaller_db_total_num_cols_rounded_up * db_cols_prime {
                smaller_db.push((rng_pub.sample::<u64, _>(rand::distributions::Standard) % params.pt_modulus) as u16);
            }

            assert_eq!(smaller_db.len(), smaller_db_total_num_cols_rounded_up * db_cols_prime);
            let smaller_server: YServer<u16> = YServer::<u16>::new_small(
                &self.smaller_params,
                smaller_db.into_iter(),
                ProtocolType::InsPIRe,
                true,
            );

            let vec_size = 1 << smaller_server.params.db_dim_1;
            let mut pseudorandom_query_1 = Vec::with_capacity(vec_size);
            for _ in 0..vec_size {
                pseudorandom_query_1.push(PolyMatrixNTT::random(&params, 1, 1))
            }

            let mut precomp_inspir_vec_first_layer = Vec::new();
            let thing_size = if gammas[0] == params.poly_len { params.poly_len / 2 } else {gammas[0]};
            for _ in 0..db_cols_prime {
                precomp_inspir_vec_first_layer.push(PrecompInsPIR {
                    a_hat: PolyMatrixRaw::random(&params, 1, 1),
                    bold_t_condensed: PolyMatrixNTT::random(&params, thing_size - 1, params.t_exp_left),
                    bold_t_bar_condensed: PolyMatrixNTT::random(&params, thing_size - 1, params.t_exp_left),
                    bold_t_hat_condensed: PolyMatrixNTT::random(&params, 1, params.t_exp_left)
                });
            }

            let mut precomp_inspir_vec = Vec::new();
            let thing_size = if gammas[1] == params.poly_len { params.poly_len / 2 } else {gammas[1]};
            for _ in 0..mask_num_packed_rlwes {
                precomp_inspir_vec.push(PrecompInsPIR {
                    a_hat: PolyMatrixRaw::random(&params, 1, 1),
                    bold_t_condensed: PolyMatrixNTT::random(&params, thing_size - 1, params.t_exp_left),
                    bold_t_bar_condensed: PolyMatrixNTT::random(&params, thing_size - 1, params.t_exp_left),
                    bold_t_hat_condensed: PolyMatrixNTT::random(&params, 1, params.t_exp_left)
                });
            }
            let offline_packing_keys = if gammas.len() > 2 {
                if gammas[2] <= params.poly_len / 2 {
                    OfflinePackingKeys::init(&self.packing_params_set[&gammas[2]], W_SEED)
                } else {
                    OfflinePackingKeys::init_full(&self.packing_params_set[&gammas[2]], W_SEED, V_SEED)
                }
            } else {
                OfflinePackingKeys::init_empty()
            };

            return OfflinePrecomputedValues {
                hint_1,
                pseudorandom_query_1,
                y_constants: (vec![], vec![]),
                smaller_server: Some(smaller_server),
                fake_pack_pub_params: vec![],
                precomp: vec![],
                precomp_inspir_vec_first_layer,
                precomp_inspir_vec,
                offline_packing_keys: offline_packing_keys
            }            
        }

        // Begin offline precomputation

        let now: Instant = Instant::now();
        let hint_0: Vec<u64> = self.answer_hint_ring(SEED_0, db_cols);

        // hint_0 is poly_len x db_cols
        if let Some(measurement) = measurement {
            measurement.offline.simplepir_prep_time_ms = now.elapsed().as_millis() as usize;
        }
        debug!("Answered hint (ring) in {} us", now.elapsed().as_micros());

        let combined = [&hint_0[..], &vec![0u64; db_cols]].concat();

        assert_eq!(combined.len(), (params.poly_len + 1) * db_cols);

        let gamma_first_layer = gammas[0];

        let prepacked_lwe = prep_pack_many_lwes(&params, &combined, db_cols_prime, gamma_first_layer);

        let now = Instant::now();

        let mut offline_packing_keys_set: HashMap<usize, OfflinePackingKeys> = HashMap::new(); // Assuming gamma is f64
        for gamma in gammas {
            if !offline_packing_keys_set.contains_key(&gamma) {
                let offline_packing_keys = if *gamma <= params.poly_len / 2 {
                    OfflinePackingKeys::init(&self.packing_params_set[&gamma], W_SEED)
                } else {
                    OfflinePackingKeys::init_full(&self.packing_params_set[&gamma], W_SEED, V_SEED)
                };
                offline_packing_keys_set.insert(*gamma, offline_packing_keys);
            } 
        }

        let mut precomp_inspir_vec_first_layer = Vec::new();

        let offline_packing_key_first_layer = offline_packing_keys_set.get(&gamma_first_layer).unwrap();
        for i in 0..db_cols_prime {
            let mut a_ct_tilde = Vec::new();
            for j in 0..gamma_first_layer {
                a_ct_tilde.push(prepacked_lwe[i][j].submatrix(0, 0, 1, 1))
            }

            let tup_inspir = if gamma_first_layer <= params.poly_len / 2 {
                let w_all = offline_packing_key_first_layer.w_all.as_ref().unwrap();
                packing_with_preprocessing_offline(&self.packing_params_set[&gamma_first_layer], &w_all, &a_ct_tilde)
            } else {
                let w_all = offline_packing_key_first_layer.w_all.as_ref().unwrap();
                let w_bar_all = offline_packing_key_first_layer.w_bar_all.as_ref().unwrap();
                let v_mask = offline_packing_key_first_layer.v_mask.as_ref().unwrap();

                full_packing_with_preprocessing_offline(
                    &self.packing_params_set[&gamma_first_layer],
                    &w_all,
                    &w_bar_all,
                    &v_mask,
                    &a_ct_tilde,
                )
            };
            precomp_inspir_vec_first_layer.push(tup_inspir);
        }
        debug!("Precomp in {} us", now.elapsed().as_micros());

        // compute (most of) the secondary hint
        // let intermediate_cts = [&hint_0[..], &vec![0u64; db_cols]].concat();

        let new_hint = precomp_inspir_vec_first_layer
            .iter()
            .map(|x| [x.a_hat.get_poly(0, 0).to_vec(), vec![0u64; gamma_first_layer]].concat())
            .collect::<Vec<_>>()
            .concat();

        // new_hint is num_rlwe_outputs * (poly_len + packed_per_rlwe)

        let intermediate_cts_rescaled = new_hint
            .iter()
            .map(|x| rescale(*x, params.modulus, rlwe_q_prime_2))
            .collect::<Vec<_>>();
        // let intermediate_cts_rescaled = [&intermediate_cts_rescaled[..], &vec![0u64;
        // num_rlwe_outputs*packed_per_rlwe]].concat();

        // split and do a second PIR over intermediate_cts
        // split into blowup_factor=q/p instances (so that all values are now mod p)
        // the second PIR is over a database of db_cols x (blowup_factor * (lwe_params.n
        // + 1)) values mod p

        // inp: (params.poly_len + 1, db_cols)
        // out: (out_rows >= (params.poly_len + \gamma) * blowup_factor, db_cols)
        //      we are 'stretching' the columns (and padding)

        debug!("Splitting intermediate cts...");

        let intermediate_cts_rescaled_transposed  = transpose_generic(&intermediate_cts_rescaled, db_cols_prime, params.poly_len + gamma_first_layer);
        let smaller_db: Vec<u16> = split_alloc(
            &intermediate_cts_rescaled_transposed,
            mask_decomposed_num_coeffs * pt_bits,
            params.poly_len + gamma_first_layer,
            params.poly_len,
            db_cols_prime,
            smaller_db_total_num_cols_rounded_up,
            q2_bits,
            pt_bits,
        );
        assert_eq!(smaller_db.len(), smaller_db_total_num_cols_rounded_up * db_cols_prime);

        debug!("Done splitting intermediate cts.");

        let gamma_second_layer_mask = gammas[1];
        // This is the 'intermediate' db after the first pass of PIR and expansion
        let smaller_server: YServer<u16> = YServer::<u16>::new_small(
            &self.smaller_params,
            smaller_db.into_iter(),
            ProtocolType::InsPIRe,
            true,
        );
        debug!("gen'd smaller server.");

        let hint_1 = smaller_server.answer_hint_ring(SEED_1, smaller_db_total_num_cols_rounded_up);

        assert_eq!(hint_1.len(), params.poly_len * smaller_db_total_num_cols_rounded_up);

        let pseudorandom_query_1 = smaller_server.generate_pseudorandom_query(SEED_1);

        let combined = [&hint_1[..], &vec![0u64; smaller_db_total_num_cols_rounded_up]].concat();
        assert_eq!(combined.len(), smaller_db_total_num_cols_rounded_up * (params.poly_len + 1));
        // Instead of rho, it should probably only be blowup_factor * poly_len
        let prepacked_lwe = prep_pack_many_lwes(
            &params,
            &combined,
            smaller_db_total_num_cols_rounded_up / gamma_second_layer_mask,
            gamma_second_layer_mask,
        );
        debug!("prepacked lwes");

        let now = Instant::now();

        let offline_packing_key_second_layer_mask = offline_packing_keys_set.get(&gammas[1]).unwrap();
        let mut precomp_inspir_vec = Vec::new();
        for i in 0..mask_num_packed_rlwes {
            // Our version
            let mut a_ct_tilde = Vec::new();
            for j in 0..gamma_second_layer_mask {
                a_ct_tilde.push(prepacked_lwe[i][j].submatrix(0, 0, 1, 1))
            }

            let tup_inspir = if gamma_second_layer_mask <= params.poly_len / 2 {
                let w_all = offline_packing_key_second_layer_mask.w_all.as_ref().unwrap();
                packing_with_preprocessing_offline(
                    &self.packing_params_set[&gamma_second_layer_mask],
                    &w_all,
                    &a_ct_tilde,
                )
            } else {
                let w_all = offline_packing_key_second_layer_mask.w_all.as_ref().unwrap();
                let w_bar_all = offline_packing_key_second_layer_mask.w_bar_all.as_ref().unwrap();
                let v_mask = offline_packing_key_second_layer_mask.v_mask.as_ref().unwrap();
                full_packing_with_preprocessing_offline(
                    &self.packing_params_set[&gamma_second_layer_mask],
                    &w_all,
                    &w_bar_all,
                    &v_mask,
                    &a_ct_tilde,
                )
            };
            precomp_inspir_vec.push(tup_inspir);
        }
        println!("Precomp in {} us", now.elapsed().as_micros());

        let offline_packing_keys = if self.second_level_packing_body == PackingType::InspiRING {
            offline_packing_keys_set.get(&gammas[2]).unwrap().clone()
        } else {
            OfflinePackingKeys::init_empty()
        };

        OfflinePrecomputedValues {
            // hint_0,
            hint_1,
            pseudorandom_query_1,
            y_constants: (vec![], vec![]),
            smaller_server: Some(smaller_server),
            fake_pack_pub_params: vec![],
            precomp: vec![],
            precomp_inspir_vec_first_layer,
            precomp_inspir_vec,
            offline_packing_keys,
        }
    }

    /// Perform SimplePIR-style YPIR
    pub fn perform_online_computation_simplepir(
        &self,
        gamma: usize,
        first_dim_queries_packed: &[u64],
        offline_vals: &OfflinePrecomputedValues<'a>,
        packing_keys: &mut PackingKeys<'a>,
        measurement: Option<&mut Measurement>,
    ) -> Vec<Vec<u8>> {
        assert_eq!(self.protocol_type, ProtocolType::SimplePIR);

        let params = self.params;

        // RLWE reduced moduli
        let rlwe_q_prime_1 = params.get_q_prime_1();
        let rlwe_q_prime_2 = params.get_q_prime_2();

        let db_rows = 1 << (params.db_dim_1 + params.poly_len_log2);
        let db_cols = params.instances * params.poly_len;
        let db_cols_prime = (db_cols as f64 / gamma as f64).ceil() as usize;

        assert_eq!(first_dim_queries_packed.len(), db_rows);

        // Begin online computation

        let first_pass = Instant::now();
        debug!("Performing mul...");
        let mut intermediate = AlignedMemory64::new(db_cols);
        fast_batched_dot_product_generic::<T>(
            &params,
            intermediate.as_mut_slice(),
            first_dim_queries_packed,
            db_rows,
            self.db(),
            db_rows,
            db_cols,
        );
        debug!("Done w mul...");

        let first_pass_time_us = first_pass.elapsed().as_micros();

        let mut packing_key_rotations_time_us = 0;
        let ring_packing = Instant::now();

        let packed = match self.second_level_packing_mask {
            PackingType::InspiRING => {
                let precomp_inspir_vec = &offline_vals.precomp_inspir_vec;
                let packing_key_rotations_time = Instant::now();
                packing_keys.expand();
                packing_key_rotations_time_us = packing_key_rotations_time.elapsed().as_micros();
                pack_many_lwes_inspir(
                    &self.packing_params_set[&gamma],
                    &precomp_inspir_vec,
                    intermediate.as_slice(),
                    &packing_keys,
                    gamma,
                )                
            },
            PackingType::CDKS => {
                let y_constants = &offline_vals.y_constants;
                let precomp = &offline_vals.precomp;
                pack_many_lwes(
                    &params,
                    &precomp,
                    intermediate.as_slice(),
                    db_cols_prime,
                    &packing_keys.pack_pub_params_row_1s,
                    &y_constants,
                )
            },
            PackingType::NoPacking => {
                panic!("Shouldn't be here!");
            },
        };

        let total_ring_packing_time_us = ring_packing.elapsed().as_micros();

        debug!("Packed...");
        if let Some(m) = measurement {
            m.online.first_pass_time_us = first_pass_time_us as usize;
            m.online.packing_key_rotations_time_us = packing_key_rotations_time_us as usize;
            m.online.first_pack_time_us = total_ring_packing_time_us as usize;
        }

        let mut packed_mod_switched = Vec::with_capacity(packed.len());
        for ct in packed.iter() {
            let res_switched = ct.switch_and_keep(rlwe_q_prime_1, rlwe_q_prime_2, gamma);
            packed_mod_switched.push(res_switched);
        }

        packed_mod_switched
    }

    pub fn perform_online_computation(
        &self,
        gamma: usize,
        packing_keys: &mut PackingKeys<'a>,
        offline_vals: &mut OfflinePrecomputedValues<'a>,
        first_dim_queries_packed: &[u32],
        second_dim_queries: &[&[u64]],
        mut measurement: Option<&mut Measurement>,
    ) -> Vec<Vec<Vec<u8>>> {
        // Set up some parameters

        let params = self.params;
        let lwe_params = LWEParams::default();

        let db_cols = self.db_cols();

        // RLWE reduced moduli
        let rlwe_q_prime_1 = params.get_q_prime_1();
        let rlwe_q_prime_2 = params.get_q_prime_2();

        // LWE reduced moduli
        let lwe_q_prime_bits = lwe_params.q2_bits as usize;
        let lwe_q_prime = lwe_params.get_q_prime_2();

        // The number of bits represented by a plaintext RLWE coefficient
        let pt_bits = (params.pt_modulus as f64).log2().floor() as usize;
        // assert_eq!(pt_bits, 16);

        // The factor by which ciphertext values are bigger than plaintext values
        let blowup_factor = lwe_q_prime_bits as f64 / pt_bits as f64;
        let blowup_factor_ceil = blowup_factor.ceil() as usize;
        debug!("blowup_factor: {}", blowup_factor);

        // The starting index of the final value (the '1' in lwe_params.n + 1)
        // This is rounded to start on a pt_bits boundary
        let special_offs =
            ((lwe_params.n * lwe_q_prime_bits) as f64 / pt_bits as f64).ceil() as usize;

        let out_rows = 1 << (self.smaller_params.db_dim_2 + params.poly_len_log2);
        // let num_rlwes = (out_rows as f64 / gamma as f64).ceil() as usize;
        let rho = 1 << self.smaller_params.db_dim_2;
        
        assert_eq!(self.smaller_params.db_dim_1, params.db_dim_2);
        assert!(out_rows as f64 >= (blowup_factor * (lwe_params.n + 1) as f64));

        // Load offline precomputed values
        let hint_1_combined = &mut offline_vals.hint_1;
        let pseudorandom_query_1 = &offline_vals.pseudorandom_query_1;
        let y_constants = &offline_vals.y_constants;
        let smaller_server = offline_vals.smaller_server.as_mut().unwrap();
        // let prepacked_lwe = &offline_vals.prepacked_lwe;
        let fake_pack_pub_params = &offline_vals.fake_pack_pub_params;
        let precomp = &offline_vals.precomp;
        let precomp_inspir_vec = &offline_vals.precomp_inspir_vec;

        // Begin online computation

        let online_phase = Instant::now();
        let first_pass = Instant::now();
        let intermediate = self.lwe_multiply_batched_with_db_packed::<1>(first_dim_queries_packed);
        let simplepir_resp_bytes = intermediate.len() / (lwe_q_prime_bits as usize) / 8;
        debug!("simplepir_resp_bytes {} bytes", simplepir_resp_bytes);
        let first_pass_time_us = first_pass.elapsed().as_micros();
        debug!("First pass took {} us", first_pass_time_us);

        debug!("intermediate.len(): {}", intermediate.len());
        let mut second_pass_time_us = 0;
        let mut total_ring_packing_time_ms = 0;

        let packing_key_rotations_time_us = 0;
        let mut packing_mask_online_time_us = 0;
        let mut packing_body_online_time_us = 0;

        let mut responses = Vec::new();
        for (intermediate_chunk, packed_query_col) in
            intermediate.as_slice().chunks(db_cols).zip(second_dim_queries.iter())
        {
            let second_pass = Instant::now();
            let intermediate_cts_rescaled = intermediate_chunk
                .iter()
                .map(|x| rescale(*x as u64, lwe_params.modulus, lwe_q_prime))
                .collect::<Vec<_>>();
            assert_eq!(intermediate_cts_rescaled.len(), db_cols);
            debug!("intermediate_cts_rescaled[0] = {}", intermediate_cts_rescaled[0]);

            let now = Instant::now();
            // modify the smaller_server db to include the intermediate values
            // let mut smaller_server_clone = smaller_server.clone();
            {
                // remember, this is stored in 'transposed' form
                // so it is out_cols x db_cols
                let smaller_db_mut: &mut [u16] = smaller_server.db_mut();
                for j in 0..db_cols {
                    // new value to write into the db
                    let val = intermediate_cts_rescaled[j];

                    for m in 0..blowup_factor_ceil as usize {
                        // index in the transposed db
                        let out_idx = (special_offs + m) * db_cols + j;

                        // part of the value to write into the db
                        let val_part = ((val >> (m * pt_bits)) & ((1 << pt_bits) - 1)) as u16;

                        // assert_eq!(smaller_db_mut[out_idx], DoubleType::default());
                        smaller_db_mut[out_idx] = val_part;
                    }
                }
            }
            debug!("load secondary hint {} us", now.elapsed().as_micros());

            let now = Instant::now();
            {

                let phase = Instant::now();
                let secondary_hint = smaller_server.multiply_with_db_ring(
                    &pseudorandom_query_1,
                    special_offs..special_offs + blowup_factor_ceil,
                    SEED_1,
                );
                debug!("multiply_with_db_ring took: {} us", phase.elapsed().as_micros());
                // let phase = Instant::now();
                // let secondary_hint =
                //     smaller_server_clone.answer_hint(SEED_1, special_offs..special_offs +
                // blowup_factor_ceil); debug!(
                //     "traditional answer_hint took: {} us",
                //     phase.elapsed().as_micros()
                // );

                assert_eq!(secondary_hint.len(), params.poly_len * blowup_factor_ceil);

                for i in 0..params.poly_len {
                    for j in 0..blowup_factor_ceil {
                        let inp_idx = i * blowup_factor_ceil + j;
                        let out_idx = i * out_rows + special_offs + j;

                        // assert_eq!(hint_1_combined[out_idx], 0); // we no longer clone for each
                        // query, just overwrite
                        hint_1_combined[out_idx] = secondary_hint[inp_idx];
                    }
                }
            }
            debug!("compute secondary hint in {} us", now.elapsed().as_micros());

            assert_eq!(hint_1_combined.len(), params.poly_len * out_rows);

            let response: AlignedMemory64 = smaller_server.answer_query(packed_query_col);

            second_pass_time_us += second_pass.elapsed().as_micros();
            let total_ring_packing = Instant::now();
            let now = Instant::now();
            assert_eq!(response.len(), 1 * out_rows);

            // combined is now (poly_len + 1) * (out_rows)
            // let combined = [&hint_1_combined[..], response.as_slice()].concat();

            let mut excess_cts = Vec::with_capacity(blowup_factor_ceil as usize);
            for j in special_offs..special_offs + blowup_factor_ceil as usize {
                let mut rlwe_ct = PolyMatrixRaw::zero(&params, 2, 1);

                // 'a' vector
                // put this in negacyclic order
                let mut poly = Vec::new();
                for k in 0..params.poly_len {
                    poly.push(hint_1_combined[k * out_rows + j]);
                }
                let nega = negacyclic_perm(&poly, 0, params.modulus);

                rlwe_ct.get_poly_mut(0, 0).copy_from_slice(&nega);

                if self.second_level_packing_body == PackingType::InspiRING {
                    rlwe_ct.get_poly_mut(1, 0)[0] = response[j];
                } else if self.second_level_packing_body == PackingType::NoPacking {
                    rlwe_ct.get_poly_mut(1, 0)[0] = response[j];
                }

                // let j_within_last = j % params.poly_len;
                // prepacked_lwe_mut.last_mut().unwrap()[j_within_last] = rlwe_ct.ntt();
                excess_cts.push(rlwe_ct);
            }
            debug!("in between: {} us", now.elapsed().as_micros());

            // assert_eq!(pack_pub_params_row_1s[0].rows, 1);
            let packing_mask = Instant::now();

            let mut packed = match self.second_level_packing_mask {
                PackingType::InspiRING => pack_many_lwes_inspir_without_rotations(
                    &self.packing_params_set[&gamma],
                    &precomp_inspir_vec,
                    &response.as_slice()[..special_offs],
                    &packing_keys,
                    gamma,
                ),
                PackingType::CDKS => {
                    let pack_pub_params_row_1s = &packing_keys.pack_pub_params_row_1s;
                    pack_many_lwes(
                        &params,
                        // &prepacked_lwe,
                        &precomp,
                        response.as_slice(),
                        rho,
                        pack_pub_params_row_1s,
                        &y_constants,
                    )
                }
                PackingType::NoPacking => {
                    println!("This shouldn't happen!");
                    Vec::new()
                }
            };
            // println!("Packing with preprocessing took {} us",
            // packing_mask_online_time_ms.elapsed().as_micros());
            packing_mask_online_time_us = packing_mask.elapsed().as_micros();

            let packing_body = Instant::now();

            match self.second_level_packing_body {
                PackingType::NoPacking => {
                    for t in 0..excess_cts.len() {
                        packed.push(excess_cts[t].clone());
                    }
                }
                PackingType::InspiRING => {
                    let mut b_poly = PolyMatrixRaw::zero(&params, 1, 1);
                    let non_zeros = excess_cts.len() as usize;
                    assert!(non_zeros <= gamma);
                    for j in 0..non_zeros {
                        b_poly.get_poly_mut(0, 0)[j] = excess_cts[j].get_poly(1, 0)[0];
                    }
                    let mut a_ct_tilde = Vec::new();
                    for j in 0..non_zeros {
                        a_ct_tilde.push(excess_cts[j].submatrix(0, 0, 1, 1).ntt());
                    }

                    let other_packed = if gamma <= params.poly_len/2 {
                        packing_fully_online_without_rotations( 
                            &self.packing_params_set[&gamma],
                            &offline_vals.offline_packing_keys,
                            &a_ct_tilde,
                            packing_keys,
                            &b_poly,
                        )
                    } else {
                        // half packing is enough because there's very few items to pack
                        half_packing_fully_online_without_rotations( 
                            &self.half_packing_params_set[&gamma],
                            &offline_vals.offline_packing_keys.w_all.as_ref().unwrap(),
                            &a_ct_tilde,
                            packing_keys,
                            &b_poly,
                        )
                    };

                    packed.push(other_packed);
                }
                PackingType::CDKS => {
                    let now = Instant::now();
                    let pack_pub_params_row_1s = &packing_keys.pack_pub_params_row_1s;
                    let mut pack_pub_params = fake_pack_pub_params.clone();
                    for i in 0..pack_pub_params.len() {
                        let uncondensed = uncondense_matrix(params, &pack_pub_params_row_1s[i]);
                        pack_pub_params[i].copy_into(&uncondensed, 1, 0);
                    }
                    debug!("uncondense pub params: {} us", now.elapsed().as_micros());

                    let other_packed = pack_using_single_with_offset(
                        &params,
                        &pack_pub_params,
                        &excess_cts.iter().map(|x| x.ntt()).collect::<Vec<_>>(),
                        special_offs,
                    );
                    let packed_clone = packed[0].clone();
                    add_raw(&mut packed[0], &packed_clone, &other_packed);
                }
            }

            debug!("pack_using_single_with_offset: {} us", now.elapsed().as_micros());

            // println!("packing body online took: {} us",
            // packing_body_online_time_ms.elapsed().as_micros());

            packing_body_online_time_us = packing_body.elapsed().as_micros();

            let mut packed_mod_switched = Vec::with_capacity(packed.len());
            for ct in packed.iter() {
                let res_switched = ct.switch_and_keep(rlwe_q_prime_1, rlwe_q_prime_2, gamma);
                packed_mod_switched.push(res_switched);
            }
            // debug!("Preprocessing pack in {} us", now.elapsed().as_micros());
            // debug!("");
            total_ring_packing_time_ms += total_ring_packing.elapsed().as_millis();

            // packed is blowup_factor ring ct's
            // these encode, contiguously [poly_len + 1, blowup_factor]
            // (and some padding)
            // match self.second_level_packing_body {
            //     PackingType::NoPacking => {
            //         assert_eq!(packed.len(), rho + blowup_factor_ceil as usize);
            //     }
            //     PackingType::InspiRING => {
            //         assert_eq!(packed.len(), rho + 1);
            //     }
            //     PackingType::CDKS => {
            //         assert_eq!(packed.len(), rho);
            //     }
            // }

            responses.push(packed_mod_switched);
        }
        debug!("Total online time: {} us", online_phase.elapsed().as_micros());
        debug!("");

        if let Some(ref mut m) = measurement {
            m.online.first_pass_time_us = first_pass_time_us as usize;
            m.online.second_pass_time_us = second_pass_time_us as usize;
            m.online.packing_key_rotations_time_us = packing_key_rotations_time_us as usize;
            m.online.packing_mask_online_time_us = packing_mask_online_time_us as usize;
            m.online.packing_body_online_time_us = packing_body_online_time_us as usize;
            m.online.total_ring_packing_time_ms = total_ring_packing_time_ms as usize;
        }

        responses
    }

    pub fn perform_online_computation_medium_payload(
        &self,
        offline_vals: &mut OfflinePrecomputedValues<'a>,
        first_dim_queries_packed: &[u64],
        second_dim_queries: &[&[u64]],
        mut packing_keys_set: HashMap<usize, PackingKeys<'a>>,
        gammas: Vec<usize>,
        mut measurement: Option<&mut Measurement>,
    // ) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    ) -> Response {
        // Set up some parameters

        let params = self.params;
        let precomp_inspir_vec_first_layer = &offline_vals.precomp_inspir_vec_first_layer;
        let precomp_inspir_vec = &offline_vals.precomp_inspir_vec;
        let hint_1_combined = &mut offline_vals.hint_1;
        let smaller_server = offline_vals.smaller_server.as_mut().unwrap();
        let pseudorandom_query_1 = &offline_vals.pseudorandom_query_1;

        let gamma_first_layer = gammas[0];
        let gamma_second_layer_mask = gammas[1];
        let gamma_second_layer_body = if self.second_level_packing_body == PackingType::InspiRING {
            gammas[2]
        } else { 1 };

        // RLWE reduced moduli
        let rlwe_q_prime_1 = params.get_q_prime_1(); // 20 bits
        let rlwe_q_prime_2 = params.get_q_prime_2(); // 28 bits

        let db_rows = 1 << (params.db_dim_1 + params.poly_len_log2);
        let db_cols = params.instances * params.poly_len;

        let db_cols_prime = (db_cols as f64 / gamma_first_layer as f64).ceil() as usize;

        let q2_bits = (rlwe_q_prime_2 as f64).log2().ceil() as usize;

        // The number of bits represented by a plaintext RLWE coefficient
        let pt_bits = (params.pt_modulus as f64).log2().floor() as usize;
        let blowup_factor = q2_bits as f64 / pt_bits as f64;

        let mask_decomposed_num_coeffs = (params.poly_len as f64 * blowup_factor).ceil() as usize;
        let body_decomposed_num_coeffs = (gamma_first_layer as f64 * blowup_factor).ceil() as usize;

        // let smaller_db_total_num_cols = mask_decomposed_num_coeffs + body_decomposed_num_coeffs;
        let smaller_db_total_num_cols_rounded_up =
            1 << (self.smaller_params.db_dim_2 + params.poly_len_log2); // out_columns after first multiplication and splitting

        // assert_eq!(mask_decomposed_num_coeffs % gamma, 0);
        let mask_num_packed_rlwes = (mask_decomposed_num_coeffs as f64 / gamma_second_layer_mask as f64).ceil() as usize;

        assert_eq!(first_dim_queries_packed.len(), db_rows);

        let online_phase = Instant::now();

        let packing_key_rotations_time = Instant::now();
        packing_keys_set.get_mut(&gamma_first_layer).unwrap().expand();
        packing_keys_set.get_mut(&gamma_second_layer_mask).unwrap().expand();
        if self.second_level_packing_body == PackingType::InspiRING {
            packing_keys_set.get_mut(&gammas[2]).unwrap().expand();
        }
        let packing_key_rotations_time_us = packing_key_rotations_time.elapsed().as_micros();
        debug!("KSK gen took {} us", packing_key_rotations_time_us);

        let first_pass = Instant::now();

        // Begin online computation

        debug!("Performing mul...");
        let mut intermediate = AlignedMemory64::new(db_cols);

        fast_batched_dot_product_avx512::<1, T>(
            &params,
            intermediate.as_mut_slice(),
            first_dim_queries_packed,
            db_rows,
            self.db(),
            db_rows,
            db_cols,
        );
        debug!("Done w mul...");

        // let simplepir_resp_bytes = intermediate.len() / (lwe_q_prime_bits as usize) / 8;
        let first_pass_time_us = first_pass.elapsed().as_micros();
        debug!("First pass took {} us", first_pass_time_us);

        let first_pack = Instant::now();

        let packed_0 = pack_many_lwes_inspir(
            &self.packing_params_set[&gamma_first_layer],
            &precomp_inspir_vec_first_layer,
            intermediate.as_slice(),
            &packing_keys_set[&gamma_first_layer],
            gamma_first_layer,
        );
        debug!("Packed...");

        let first_pack_time_us = first_pack.elapsed().as_micros();
        debug!("First pack took {} us", first_pack_time_us);

        let intermediate_cts_rescaled = packed_0
            .iter()
            .map(|x| x.get_poly(1, 0)[..gamma_first_layer].to_vec())
            .collect::<Vec<_>>()
            .concat()
            .iter()
            .map(|x| rescale(*x as u64, params.modulus, rlwe_q_prime_2))
            .collect::<Vec<_>>();

        assert_eq!(intermediate_cts_rescaled.len(), db_cols_prime * gamma_first_layer);

        let intermediate_cts_rescaled_transposed = transpose_generic(&intermediate_cts_rescaled, db_cols_prime, gamma_first_layer);

        let part_to_add = split_alloc(
            &intermediate_cts_rescaled_transposed,
            0,
            gamma_first_layer,
            gamma_first_layer,
            db_cols_prime,
            body_decomposed_num_coeffs,
            q2_bits,
            pt_bits,
        );

        assert_eq!(part_to_add.len(), body_decomposed_num_coeffs*db_cols_prime);

        let smaller_db_mut: &mut [u16] = smaller_server.db_mut();
        for j0 in 0..db_cols_prime {
            for j1 in 0..body_decomposed_num_coeffs {
                smaller_db_mut[(mask_decomposed_num_coeffs + j1) * db_cols_prime + j0] =
                    part_to_add[j1 * db_cols_prime + j0];
            }
        }
        
        let secondary_hint = smaller_server.multiply_with_db_ring(
            &pseudorandom_query_1,
            mask_decomposed_num_coeffs..mask_decomposed_num_coeffs + body_decomposed_num_coeffs,
            SEED_1,
        );

        assert_eq!(secondary_hint.len(), params.poly_len * body_decomposed_num_coeffs);

        for i in 0..params.poly_len {
            for j in 0..body_decomposed_num_coeffs {
                let inp_idx = i * body_decomposed_num_coeffs + j;
                let out_idx =
                    i * smaller_db_total_num_cols_rounded_up + mask_decomposed_num_coeffs + j;

                hint_1_combined[out_idx] = secondary_hint[inp_idx];
            }
        }

        assert_eq!(hint_1_combined.len(), params.poly_len * smaller_db_total_num_cols_rounded_up);

        let second_pass = Instant::now();

        let packed_query_col = second_dim_queries[0];

        let response: AlignedMemory64 = smaller_server.answer_query(packed_query_col);

        let second_pass_time_us = second_pass.elapsed().as_micros();

        assert_eq!(response.len(), 1 * smaller_db_total_num_cols_rounded_up);

        let packing_mask = Instant::now();
        let packed_mask = pack_many_lwes_inspir(
            &self.packing_params_set[&gamma_second_layer_mask],
            &precomp_inspir_vec,
            &response.as_slice()[..mask_decomposed_num_coeffs],
            &packing_keys_set[&gamma_second_layer_mask],
            gamma_second_layer_mask,
        );
        assert_eq!(packed_mask.len(), mask_num_packed_rlwes);

        let packing_mask_time_us = packing_mask.elapsed().as_micros();
        debug!("Second pack took {} us", packing_mask_time_us);

        let packing_body_online_time = Instant::now();

        // let part1 = Instant::now();

        let mut excess_cts = Vec::with_capacity(body_decomposed_num_coeffs);
        for j in mask_decomposed_num_coeffs
            ..mask_decomposed_num_coeffs + body_decomposed_num_coeffs as usize
        {
            let mut rlwe_ct = PolyMatrixRaw::zero(&params, 2, 1);

            let mut poly = Vec::new();
            for k in 0..params.poly_len {
                poly.push(hint_1_combined[k * smaller_db_total_num_cols_rounded_up + j]);
            }
            let nega = negacyclic_perm(&poly, 0, params.modulus);

            rlwe_ct.get_poly_mut(0, 0).copy_from_slice(&nega);
            rlwe_ct.get_poly_mut(1, 0)[0] = response[j];

            excess_cts.push(rlwe_ct);
        }

        // println!("part1: {} us", part1.elapsed().as_micros());

        let mut packed_body = Vec::new();

        match self.second_level_packing_body {
            PackingType::NoPacking => {
                for t in 0..body_decomposed_num_coeffs {
                    packed_body.push(excess_cts[t].clone());
                }
            }
            PackingType::InspiRING => {
                
                let packed_online_num = (body_decomposed_num_coeffs as f64 / gamma_second_layer_body as f64).ceil() as usize;
                for i in 0..packed_online_num {

                    // let part2 = Instant::now();
                    let mut a_ct_tilde = Vec::new();
                    let mut b_poly = PolyMatrixRaw::zero(&params, 1, 1);
                    for j in 0..gamma_second_layer_body {
                        if i * gamma_second_layer_body + j < body_decomposed_num_coeffs {
                            a_ct_tilde.push(excess_cts[i * gamma_second_layer_body + j].submatrix(0, 0, 1, 1).ntt());
                            b_poly.get_poly_mut(0, 0)[j] = excess_cts[i * gamma_second_layer_body + j].get_poly(1, 0)[0];
                        }
                    }
                    // println!("part2: {} us", part2.elapsed().as_micros());
                    
                    // let part3 = Instant::now();
                    let other_packed = packing_fully_online(
                        &self.packing_params_set[&gamma_second_layer_body],
                        &offline_vals.offline_packing_keys,
                        &a_ct_tilde,
                        &packing_keys_set[&gamma_second_layer_body],
                        &b_poly,
                    );
                    // println!("part3: {} us", part3.elapsed().as_micros());
                    packed_body.push(other_packed);
                }
            },
            PackingType::CDKS => {
                panic!("Not implemented");
            }
        };

        let packing_body_online_time_us = packing_body_online_time.elapsed().as_micros();

        let mut packed_mask_mod_switched = Vec::with_capacity(packed_mask.len());
        for ct in packed_mask.iter() {
            let res_switched = ct.switch_and_keep(rlwe_q_prime_1, rlwe_q_prime_2, gamma_second_layer_mask);
            packed_mask_mod_switched.push(res_switched);
        }

        let mut packed_body_mod_switched = Vec::with_capacity(packed_body.len());
        for ct in packed_body.iter() {
            let res_switched = ct.switch_and_keep(rlwe_q_prime_1, rlwe_q_prime_2, gamma_second_layer_body);
            packed_body_mod_switched.push(res_switched);
        }

        debug!("Total online time: {} us", online_phase.elapsed().as_micros());
        debug!("");

        if let Some(ref mut m) = measurement {
            m.online.packing_key_rotations_time_us = packing_key_rotations_time_us as usize;
            m.online.first_pass_time_us = first_pass_time_us as usize;
            m.online.second_pass_time_us = second_pass_time_us as usize;
            m.online.first_pack_time_us = first_pack_time_us as usize;
            // m.online.second_pack_time_us = second_pack_time_us as usize;
            m.online.packing_mask_online_time_us = packing_mask_time_us as usize;
            m.online.packing_body_online_time_us = packing_body_online_time_us as usize;
            // m.online.total_ring_packing_time_ms = total_ring_packing_time_ms as usize;
        }

        // assert_eq!(packed_mod_switched.len(), mask_num_packed_rlwes + body_decomposed_num_coeffs);

        // packed_mod_switched
        Response {
            packed_mask_mod_switched,
            packed_body_mod_switched,
        }

    }

    // generic function that returns a u8 or u16:
    pub fn db(&self) -> &[T] {
        unsafe {
            std::slice::from_raw_parts(
                self.db_buf_aligned.as_ptr() as *const T,
                self.db_buf_aligned.len() * 8 / std::mem::size_of::<T>(),
            )
        }
    }

    pub fn db_mut(&mut self) -> &mut [T] {
        unsafe {
            std::slice::from_raw_parts_mut(
                self.db_buf_aligned.as_ptr() as *mut T,
                self.db_buf_aligned.len() * 8 / std::mem::size_of::<T>(),
            )
        }
    }

    pub fn db_u16(&self) -> &[u16] {
        unsafe {
            std::slice::from_raw_parts(
                self.db_buf_aligned.as_ptr() as *const u16,
                self.db_buf_aligned.len() * 8 / std::mem::size_of::<u16>(),
            )
        }
    }

    pub fn db_u32(&self) -> &[u32] {
        unsafe {
            std::slice::from_raw_parts(
                self.db_buf_aligned.as_ptr() as *const u32,
                self.db_buf_aligned.len() * 8 / std::mem::size_of::<u32>(),
            )
        }
    }

    pub fn get_elem(&self, row: usize, col: usize) -> T {
        self.db()[col * self.db_rows() + row] // stored transposed
    }

    pub fn get_row(&self, row: usize) -> Vec<T> {
        let db_cols = self.db_cols();
        let mut res = Vec::with_capacity(db_cols);
        for col in 0..db_cols {
            res.push(self.get_elem(row, col));
        }
        res
        // // convert to u8 contiguously
        // let mut res_u8 = Vec::with_capacity(db_cols *
        // std::mem::size_of::<T>()); for &x in res.iter() {
        //     res_u8.extend_from_slice(&x.to_u64().to_le_bytes()[..
        // std::mem::size_of::<T>()]); }
        // res_u8
    }
}

#[cfg(not(target_feature = "avx2"))]
#[allow(non_camel_case_types)]
type __m512i = u64;

pub trait ToM512 {
    fn to_m512(self) -> __m512i;
}

#[cfg(target_feature = "avx512f")]
mod m512_impl {
    use super::*;

    impl ToM512 for *const u8 {
        #[inline(always)]
        fn to_m512(self) -> __m512i {
            unsafe {
                let x = _mm_loadl_epi64(self as *const _);
                _mm512_cvtepu8_epi64(x)
            }
        }
    }

    impl ToM512 for *const u16 {
        #[inline(always)]
        fn to_m512(self) -> __m512i {
            unsafe { _mm512_cvtepu16_epi64(_mm_load_si128(self as *const _)) }
        }
    }

    impl ToM512 for *const u32 {
        #[inline(always)]
        fn to_m512(self) -> __m512i {
            unsafe { _mm512_cvtepu32_epi64(_mm256_load_si256(self as *const _)) }
        }
    }
}

#[cfg(not(target_feature = "avx512f"))]
mod m512_impl {
    use super::*;

    impl ToM512 for *const u8 {
        #[inline(always)]
        fn to_m512(self) -> __m512i {
            self as __m512i
        }
    }

    impl ToM512 for *const u16 {
        #[inline(always)]
        fn to_m512(self) -> __m512i {
            self as __m512i
        }
    }

    impl ToM512 for *const u32 {
        #[inline(always)]
        fn to_m512(self) -> __m512i {
            self as __m512i
        }
    }
}

pub trait ToU64 {
    fn to_u64(self) -> u64;
}

impl ToU64 for u8 {
    fn to_u64(self) -> u64 {
        self as u64
    }
}

impl ToU64 for u16 {
    fn to_u64(self) -> u64 {
        self as u64
    }
}

impl ToU64 for u32 {
    fn to_u64(self) -> u64 {
        self as u64
    }
}

impl ToU64 for u64 {
    fn to_u64(self) -> u64 {
        self
    }
}
