use inspire::packing::ToStr;

use clap::Parser;
use inspire::packing::PackingType;
use inspire::scheme::ProtocolType;

use log::debug;
use rand::{thread_rng, Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;
use spiral_rs::arith::multiply_uint_mod;
use spiral_rs::arith;

use spiral_rs::aligned_memory::AlignedMemory64;
use spiral_rs::{
    client::*,
    params::*,
    poly::*,
    gadget::*,
};

use spiral_rs::poly::{PolyMatrix, PolyMatrixRaw, PolyMatrixNTT};
use spiral_rs::number_theory::invert_uint_mod;

use inspire::{bits::*, client::*, kernel::*, measurement::*, modulus_switch::*, packing::*, params::*, scheme::*, server::*};

use inspire::interpolate::*;

use std::{marker::PhantomData, time::Instant};
use std::collections::HashMap;

pub const RGSW_SEEDS: [[u8; 32]; 6] = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,],
];

pub trait RGSWPIR<'a> {
    fn new_rgswpir<'b>(
        params: &'a Params,
        db: &[u16],
        interpolate_degree: usize,
    ) -> Self;

    fn perform_online_computation_simplepir_and_rgsw(
        &self,
        gamma: usize,
        first_dim_queries_packed: &[u64],
        ct_gsw: PolyMatrixNTT<'a>,
        offline_vals: &OfflinePrecomputedValues<'a>,
        packing_keys: &mut PackingKeys<'a>,
        measurement: Option<&mut Measurement>,
    ) -> Vec<Vec<u8>>;
    
}

impl<'a, T: Sync> RGSWPIR<'a> for YServer<'a, T> where
T: Sized + Copy + ToU64 + Default + std::marker::Sync,
*const T: ToM512,
u64: From<T>,
{
    fn new_rgswpir<'b>(
        params: &'a Params,
        db: &[u16],
        interpolate_degree: usize,
    ) -> Self {
        let second_level_packing_mask = PackingType::InspiRING;
        let bytes_per_pt_el = std::mem::size_of::<u16>();
        debug!("bytes_per_pt_el: {}", bytes_per_pt_el);

        let poly_len = params.poly_len;

        let db_rows = 1 << (params.db_dim_1 + params.poly_len_log2);
        let db_cols = params.instances * params.poly_len;

        let sz_bytes = db_rows * db_cols * bytes_per_pt_el;

        let gamma = poly_len;
        let db_cols_prime = db_cols / gamma;
        let c = db_cols_prime / interpolate_degree;

        let mut db_buf_aligned = AlignedMemory64::new(sz_bytes / 8);
        let db_buf_mut = as_bytes_mut(&mut db_buf_aligned);
        let db_buf_ptr = db_buf_mut.as_mut_ptr() as *mut u16;

        let can_use_pt_ntt = params.pt_modulus % (2*poly_len as u64) == 1;

        let mod_str = format!("[{}]", params.pt_modulus);
        let plaintext_params = if can_use_pt_ntt {
            internal_params_for(poly_len, 0, 0, 1, 28, 3, &mod_str)
        } else {
            let modulus = params.pt_modulus;
            let moduli = [params.pt_modulus; 4];
            let modulus_log2 = arith::log2_ceil(modulus);
            let (barrett_cr_0, barrett_cr_1) = arith::get_barrett(&moduli);
            let (barrett_cr_0_modulus, barrett_cr_1_modulus) = arith::get_barrett_crs(modulus);
            Params {
                poly_len: poly_len,
                poly_len_log2: params.poly_len_log2,
                ntt_tables: Vec::new(),
                scratch: Vec::new(),
                crt_count: 1,
                barrett_cr_0: barrett_cr_0,
                barrett_cr_1: barrett_cr_1,
                barrett_cr_0_modulus: barrett_cr_0_modulus,
                barrett_cr_1_modulus: barrett_cr_1_modulus,
                mod0_inv_mod1: 0,
                mod1_inv_mod0: 0,
                moduli: moduli,
                modulus: params.pt_modulus,
                modulus_log2: modulus_log2 ,
                noise_width: 0.,
                n: 0,
                pt_modulus: 0,
                q2_bits: 0,
                t_conv: 0,
                t_exp_left: 0,
                t_exp_right: 0,
                t_gsw: 0,
                expand_queries: false,
                db_dim_1: 0,
                db_dim_2: 0,
                instances: 0,
                db_item_size: 0,
                version: 0
            }
        };


        let mut monomials = Vec::with_capacity(2*poly_len);

        if can_use_pt_ntt {
            for i in 0..poly_len {
                let mut monomial = PolyMatrixRaw::zero(&plaintext_params, 1, 1);
                monomial.get_poly_mut(0, 0)[i] = 1;
                monomials.push(monomial.ntt());
            }
            for i in 0..poly_len {
                monomials.push(-&monomials[i]);
            }
        }

        let num_inv = invert_uint_mod(interpolate_degree as u64, plaintext_params.modulus).unwrap();

        let step = db_rows / 8;
        for i in 0..db_rows {
            if i % step == 0 {
                print!("i={} -> ",i);
            }

            for which_poly in 0..c {

                let mut points = Vec::new();
    
                for j_prime in 0..interpolate_degree {
                    let mut point = PolyMatrixRaw::zero(&plaintext_params, 1, 1);
                    for k in 0..gamma {
                        let j = which_poly * interpolate_degree * gamma + j_prime * gamma + k;
                        let idx = j * db_rows + i;
                        point.get_poly_mut(0, 0)[k] = db[idx] as u64;
                    }
                    points.push(point);
                }
                let coeffs = if can_use_pt_ntt {
                    cooley_tukey(
                        &plaintext_params,
                        points.iter().map(|x| x.ntt()).collect::<Vec<_>>(),
                        monomials.as_slice()).iter().map(|x| x.raw()).collect::<Vec<_>>()
                } else {
                    cooley_tukey_without_ntt(&plaintext_params, points)
                } ;
    
                for j_prime in 0..interpolate_degree {
                    let raw_coeff = &coeffs[j_prime];
                    let temp = raw_coeff.get_poly(0, 0);
                    for k in 0..gamma {
                        let j = which_poly * interpolate_degree * gamma + j_prime * gamma + k;
                        let idx = j * db_rows + i;
                        
                        unsafe {
                            *db_buf_ptr.add(idx) = multiply_uint_mod(temp[k] as u64, num_inv, plaintext_params.modulus) as u16;
                        }
                        
                    }
                }

            }

        }

        // Parameters for the second round (the "DoublePIR" round)
        let smaller_params = params.clone();

        let mut packing_params_set: HashMap<usize, PackParams> = HashMap::new(); 
        let mut half_packing_params_set: HashMap<usize, PackParams> = HashMap::new(); 
        for gamma in vec![poly_len] {
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
            protocol_type: ProtocolType::SimplePIR,
            second_level_packing_mask,
            second_level_packing_body:PackingType::NoPacking,
        }
    }

    /// Perform SimplePIR-style YPIR
    fn perform_online_computation_simplepir_and_rgsw(
        &self,
        interpolate_degree: usize,
        first_dim_queries_packed: &[u64],
        ct_gsw_body: PolyMatrixNTT<'a>,
        offline_vals: &OfflinePrecomputedValues<'a>,
        packing_keys: &mut PackingKeys<'a>,
        measurement: Option<&mut Measurement>,
    ) -> Vec<Vec<u8>> {
        assert_eq!(self.protocol_type, ProtocolType::SimplePIR);

        let params = self.params;
        let gamma = params.poly_len;

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

        let rgsw_time = Instant::now();


        let mut ct_gsw = ct_gsw_body.pad_top(1);

        for i in 0..ct_gsw.cols {
            let a = PolyMatrixRaw::random_rng(params, 1, 1, &mut ChaCha20Rng::from_seed(RGSW_SEEDS[i]));
            ct_gsw.copy_into(&(-&a).ntt(), 0, i);
        }

        let ell = ct_gsw.cols / 2;
        let mut ginv_c = PolyMatrixRaw::zero(&params, 2 * ell, 1);
        let mut ginv_c_ntt = PolyMatrixNTT::zero(&params, 2 * ell, 1);
        
        let mut results = Vec::new();

        assert_eq!(db_cols_prime, packed.len());
        assert!(db_cols_prime % interpolate_degree == 0);
        let c = db_cols_prime / interpolate_degree;

        for which_poly in 0..c {
            let mut sum = PolyMatrixRaw::zero(&params, 2, 1);
            for i in (0..interpolate_degree).rev() {
                let mut prod = PolyMatrixNTT::zero(&params, 2, 1);
                gadget_invert(&mut ginv_c, &sum);
                to_ntt(&mut ginv_c_ntt, &ginv_c);
                multiply(
                    &mut prod,
                    &ct_gsw,
                    &ginv_c_ntt,
                );
                sum = &prod.raw() + &packed[which_poly * interpolate_degree + i];
                sum.reduce_mod(params.modulus);
            }
            results.push(sum);
        }
        
        let total_rgsw_time_us = rgsw_time.elapsed().as_micros();
        debug!("RGSW Time: {} us", total_rgsw_time_us);

        debug!("Packed...");
        if let Some(m) = measurement {
            m.online.first_pass_time_us = first_pass_time_us as usize;
            m.online.packing_key_rotations_time_us = packing_key_rotations_time_us as usize;
            m.online.first_pack_time_us = total_ring_packing_time_us as usize;
            m.online.rgsw_time_us = total_rgsw_time_us as usize;
        }

        let mut packed_mod_switched = Vec::with_capacity(results.len());
        for ct in results.iter() {
            let res_switched = ct.switch_and_keep(rlwe_q_prime_1, rlwe_q_prime_2, gamma);
            packed_mod_switched.push(res_switched);
        }

        packed_mod_switched

    }
}

pub fn serialize_everything(
    params: &Params,
    packing_keys: &mut PackingKeys<'_>,
    packed_query_row: AlignedMemory64,
    ct_gsw_body: PolyMatrixNTT<'_>,
) -> Vec<u8> {
    let mut all_u64 = Vec::new();
    let db_rows = 1 << (params.db_dim_1 + params.poly_len_log2);
    for i in 0..params.t_exp_left {
        for k in 0..params.poly_len {
            all_u64.push(packing_keys.y_body_condensed.as_ref().unwrap().get_poly(0, i)[k]);
            all_u64.push(packing_keys.z_body_condensed.as_ref().unwrap().get_poly(0, i)[k]);
        }
    }
    for i in 0..db_rows {
        all_u64.push(packed_query_row[i]);
    }
    // let ct_gsw_body_raw = ct_gsw_body.raw();
    for i in 0..2*params.t_gsw {
        let poly = ct_gsw_body.get_poly(0, i);
        for k in 0..2*params.poly_len {
            all_u64.push(poly[k]);
        }
    }
    all_u64.iter().flat_map(|&n| n.to_ne_bytes()).collect()
}

pub fn deserialize_everything<'a>(
    params: &'a Params,
    all_u8: Vec<u8>,
) -> (PackingKeys<'a>, AlignedMemory64, PolyMatrixNTT<'a>) {

    let all_u64: Vec<_> = all_u8.as_slice()
     .chunks_exact(8)
     .map(|chunk| u64::from_ne_bytes(chunk.try_into().unwrap()))
     .collect();
    

    let mut cnt = 0; 
    let db_rows = 1 << (params.db_dim_1 + params.poly_len_log2);
    let packing_params = PackParams::new(&params, params.poly_len);

    let mut y_body_condensed = PolyMatrixNTT::zero(&params, 1, params.t_exp_left);
    let mut z_body_condensed = PolyMatrixNTT::zero(&params, 1, params.t_exp_left);

    for i in 0..params.t_exp_left {
        for k in 0..params.poly_len {
            y_body_condensed.get_poly_mut(0, i)[k] = all_u64[cnt];
            cnt += 1;
            z_body_condensed.get_poly_mut(0, i)[k] = all_u64[cnt];
            cnt += 1;
        }
    }

    let packing_keys = PackingKeys {
        packing_type: PackingType::InspiRING,
        packing_params: Some(packing_params),
        full_key: true,
        y_body: Some(y_body_condensed.clone()),
        z_body: Some(z_body_condensed.clone()),
        y_body_condensed: Some(y_body_condensed),
        z_body_condensed: Some(z_body_condensed),
        expanded: false,
        y_all_condensed: None,
        y_bar_all_condensed: None,
        params: None,
        pack_pub_params_row_1s: vec![],
        fake_pack_pub_params: vec![],
    };
    let mut packed_query_row = AlignedMemory64::new(db_rows);
    packed_query_row
        .as_mut_slice()
        .copy_from_slice(&all_u64[cnt..cnt+db_rows]);

    cnt += db_rows;

    let mut ct_gsw_body_raw = PolyMatrixNTT::zero(&params, 1, 2*params.t_gsw); // before_gsw.clone();
    for i in 0..2*params.t_gsw {
        let poly = ct_gsw_body_raw.get_poly_mut(0, i);
        for k in 0..2*params.poly_len {
            poly[k] = all_u64[cnt];
            cnt += 1;
        }
    }
    (packing_keys, packed_query_row, ct_gsw_body_raw)
}


pub fn run_simple_ypir_rgsw_on_params(
    params: Params,
    interpolate_degree: usize,
    trials: usize,
    online_only: bool
) -> Measurement {
    let gamma = params.poly_len;
    let packing_type = PackingType::InspiRING;

    let db_rows = 1 << (params.db_dim_1 + params.poly_len_log2);
    let db_cols = params.instances * params.poly_len;

    let mut rng = thread_rng();

    // RLWE reduced moduli
    let rlwe_q_prime_1 = params.get_q_prime_1();
    let rlwe_q_prime_2 = params.get_q_prime_2();

    let db_cols_prime = db_cols / gamma;
    let c = db_cols_prime / interpolate_degree;

    let now = Instant::now();
    type T = u16;
    // let usable_modulus = (2.0 as f64).powi((params.pt_modulus as f64).log2().floor() as i32) as u64;
    // let usable_modulus = 256;
    // println!("usable_modulus: {}", usable_modulus);
    let mut pt_iter = std::iter::repeat_with(|| (T::sample() as u64 % params.pt_modulus) as T);
    let mut actual_db = Vec::with_capacity(db_cols*db_rows);
    for _ in 0..db_cols*db_rows {
        actual_db.push(pt_iter.next().unwrap());
    }
    let start_encode_comp = Instant::now();
    let y_server = YServer::<T>::new_rgswpir(
        &params,
        &actual_db,
        interpolate_degree,
    );
    let encode_server_time_ms = start_encode_comp.elapsed().as_millis();
    println!("Done encode...");

    println!("Created server in {} us", now.elapsed().as_micros());
    debug!(
        "Database of {} bytes",
        y_server.db().len() * (params.pt_modulus as f64).log2().ceil() as usize / 8
    );
    // assert_eq!(
    //     y_server.db().len() * std::mem::size_of::<T>(),
    //     db_rows_padded * db_cols * (params.pt_modulus as f64).log2().ceil() as
    // usize / 8 );
    assert_eq!(y_server.db().len(), db_rows * db_cols);

    // ================================================================
    // OFFLINE PHASE
    // ================================================================
    let mut measurements = vec![Measurement::default(); trials + 1];

    println!("Starting offline...");
    let start_offline_comp = Instant::now();
    let offline_values = y_server
        .perform_offline_precomputation_simplepir(gamma, Some(&mut measurements[0]), online_only);
    let offline_server_time_ms = start_offline_comp.elapsed().as_millis();
    println!("Done offline...");

    for trial in 0..trials + 1 {
        debug!("trial: {}", trial);
        let mut measurement = &mut measurements[trial];
        measurement.offline.encode_time_ms = encode_server_time_ms as usize;
        measurement.offline.server_time_ms = offline_server_time_ms as usize;

        // ================================================================
        // QUERY GENERATION PHASE
        // ================================================================
        // let mut queries = Vec::new();

        let mut client = Client::init(&params);

        let target_idx: usize = rng.gen::<usize>() % (db_rows * db_cols);
        let target_row = target_idx / db_cols;
        let target_col = target_idx % db_cols;

        let target_sub_col = (target_col % (interpolate_degree * gamma)) / gamma;
        debug!("Target item: {} ({}, {} ({}))", target_idx, target_row, target_col, target_sub_col);

        let start = Instant::now();
        client.generate_secret_keys();
        let sk_reg = client.get_sk_reg();

        let packing_keys = match packing_type {
            PackingType::InspiRING => {
                if gamma <= params.poly_len / 2 {
                    PackingKeys::init(&y_server.packing_params_set[&gamma], sk_reg, W_SEED)
                } else {
                    PackingKeys::init_full(&y_server.packing_params_set[&gamma], sk_reg, W_SEED, V_SEED)
                }
            },
            PackingType::CDKS => {
                PackingKeys::init_cdks(&params, sk_reg, STATIC_SEED_2)
            },
            PackingType::NoPacking => {
                panic!("Shouldn't be here");
            }
        };

        let y_client = YClient::new(&mut client, &params);
        let mut ct_gsw_body = PolyMatrixNTT::zero(&params, 1, 2 * params.t_gsw);
        debug!("t_gsw: {}", params.t_gsw);

        let bits_per = get_bits_per(&params, params.t_gsw);
        for j in 0..params.t_gsw {
            let mut sigma = PolyMatrixRaw::zero(&params, 1, 1);
            let exponent = (2 * params.poly_len * target_sub_col / interpolate_degree) % (2 * params.poly_len);
            sigma.get_poly_mut(0, 0)[exponent % params.poly_len] = if exponent < params.poly_len {
                1u64 << (bits_per * j)
            } else {
                params.modulus - (1u64 << (bits_per * j))
            };
            let sigma_ntt = sigma.ntt();
            let ct = y_client.client().encrypt_matrix_reg(
                &sigma_ntt,
                &mut ChaCha20Rng::from_entropy(),
                &mut ChaCha20Rng::from_seed(RGSW_SEEDS[2*j+1]),
            );
            ct_gsw_body.copy_into(&ct.submatrix(1, 0, 1, 1), 0, 2 * j + 1);
            
            let prod = &y_client.client().get_sk_reg().ntt() * &sigma_ntt;
            let ct = &y_client.client().encrypt_matrix_reg(
                &prod,
                &mut ChaCha20Rng::from_entropy(),
                &mut ChaCha20Rng::from_seed(RGSW_SEEDS[2*j]),
            );
            ct_gsw_body.copy_into(&ct.submatrix(1, 0, 1, 1), 0, 2 * j);
        }        

        let pub_params_size = packing_keys.get_size_bytes();
        debug!("RGSW Size: {} B", get_vec_pm_size_bytes(&[ct_gsw_body.clone()]));
        debug!("pub params size: {} bytes", pub_params_size);

        let query_row = y_client.generate_query_over_prime(
            SEED_0,
            params.db_dim_1,
            packing_type,
            target_row,
        );

        assert_eq!(query_row.len(), (params.poly_len + 1) * db_rows);
        let query_row_last_row: &[u64] = &query_row[params.poly_len * db_rows..];
        
        assert_eq!(query_row_last_row.len(), db_rows);
        let packed_query_row = pack_query(&params, query_row_last_row);

        let query_size =
            ((packed_query_row.len() as f64 * params.modulus_log2 as f64) / 8.0).ceil() as usize
            + get_vec_pm_size_bytes(&[ct_gsw_body.clone()]);

        measurement.online.client_query_gen_time_ms = start.elapsed().as_millis() as usize;
        debug!("Generated query in {} us", start.elapsed().as_micros());

        let online_upload_bytes = query_size + pub_params_size;
        debug!("Query size: {} bytes", online_upload_bytes);

        let serialized = serialize_everything(&params, &mut packing_keys.clone(), packed_query_row.clone(), ct_gsw_body.clone());
        let (mut packing_keys, packed_query_row, ct_gsw_body) = deserialize_everything(&params, serialized);

        // ================================================================
        // ONLINE PHASE
        // ================================================================

        println!("Starting online...");
        let start_online_comp = Instant::now();

        let sum_switched = y_server.perform_online_computation_simplepir_and_rgsw(
            interpolate_degree,
            packed_query_row.as_slice(),
            ct_gsw_body,
            &offline_values,
            &mut packing_keys,
            Some(&mut measurement),
        );

        let online_server_time_ms = start_online_comp.elapsed().as_millis();
        let online_download_bytes = get_size_bytes(&vec![sum_switched.clone()]);
        println!("Done online...");

        let (target_row, target_col) = (target_idx / db_cols, target_idx % db_cols);
        let target_sub_col = (target_col % (interpolate_degree * gamma)) / gamma;
        let mut correct_from_db = Vec::with_capacity(db_cols);
        for i in 0..db_cols {
            correct_from_db.push(actual_db[i * db_rows + target_row].to_u64());
        }
        let mut sub_corr_result = Vec::new();
        for which_poly in 0..c {
            let sub_index = which_poly * interpolate_degree * gamma + target_sub_col * gamma;
            for k in 0..gamma {
                sub_corr_result.push(correct_from_db[sub_index+k]);
            }
        }
        
        let start_decode = Instant::now();

        let mut results = Vec::new();
        for which_poly in 0..c {
            let sum = PolyMatrixRaw::recover_how_many(&params, rlwe_q_prime_1, rlwe_q_prime_2, gamma, sum_switched[which_poly].as_slice());
            results.push(sum);
        }

        let rgsw_ans = results.iter().flat_map(|ct| {
            decrypt_ct_reg_measured(y_client.client(), &params, &ct.ntt(), params.poly_len).as_slice()[..gamma].to_vec()
        }).collect::<Vec<_>>();
        
        let final_result = rgsw_ans.as_slice();
        measurement.online.client_decode_time_ms = start_decode.elapsed().as_millis() as usize;

        if online_only {
            println!("Warning!")
        } else {
            assert_eq!(final_result, sub_corr_result);
        }

        measurement.online.upload_keys = pub_params_size;
        measurement.online.upload_query = query_size;

        measurement.online.upload_bytes = online_upload_bytes;
        measurement.online.download_bytes = online_download_bytes;
        measurement.online.total_bytes = online_upload_bytes + online_download_bytes;
        measurement.online.server_time_ms = online_server_time_ms as usize;
    }

    // discard the first measurement (if there were multiple trials)
    // copy offline values from the first measurement to the second measurement
    if trials > 1 {
        measurements[1].offline = measurements[0].offline.clone();
        measurements.remove(0);
    }

    let mut final_measurement = measurements[0].clone();
    final_measurement.online.server_time_ms =
        mean(&measurements.iter().map(|m| m.online.server_time_ms).collect::<Vec<_>>()).round()
            as usize;
    final_measurement.online.all_server_times_ms =
        measurements.iter().map(|m| m.online.server_time_ms).collect::<Vec<_>>();
    final_measurement.online.std_dev_server_time_ms =
        std_dev(&final_measurement.online.all_server_times_ms);

    final_measurement.specs.protocol_type = ProtocolType::SimplePIR.to_str();
    final_measurement.specs.second_level_packing_mask = PackingType::NoPacking.to_str();
    final_measurement.specs.second_level_packing_body = PackingType::NoPacking.to_str();
    final_measurement.specs.poly_len = params.poly_len;
    final_measurement.specs.modulus_bits = params.modulus_log2 as usize;
    final_measurement.specs.pt_modulus = params.pt_modulus as usize;

    final_measurement.specs.gamma_0 = gamma;

    final_measurement
}

fn max_interpolate_degree(modulus: f64, d0: f64, p: f64, t_exp: f64, t_rgsw: f64, poly_len: f64) -> usize {
    // let sigma_x = 16.042421 as f64;
    let sigma_x = 6.4 as f64;
    let modulus_len = modulus.log2();
    let z1 = (2.0 as f64).powi((modulus_len / t_exp).ceil() as i32);
    let z2 = (2.0 as f64).powi((modulus_len / t_rgsw).ceil() as i32);
    let term1_variance = d0 * p.powi(2) * sigma_x.powi(2);
    let term2_variance = t_exp * poly_len.powi(2) * z1.powi(2) * sigma_x.powi(2) / 4.0;
    let term3_variance = t_rgsw * poly_len * z2.powi(2) * sigma_x.powi(2) / 2.0;

    debug!("d0 = {}", d0);
    debug!("-- term1 std: {}", term1_variance.log2()/2.);
    debug!("-- term2 std: {}", term2_variance.log2()/2.);
    debug!("-- term3 std: {}", term3_variance.log2()/2.);
    
    let total_log2_std_before_poly_eval = (term1_variance + term2_variance + term3_variance).sqrt().log2();
    // let error_rate = 31.0 * (2 as f64).ln();

    let log2_std_upper_bound = (modulus / (2.*2.*p)).log2() - (2.0 * 41.0 * (2 as f64).ln()).sqrt().log2();
    debug!("log2_std_upper_bound: {}", log2_std_upper_bound);
    debug!("total_log2_std_before_poly_eval: {}", total_log2_std_before_poly_eval);
    let max_log_interpolate_degree = 2. * (log2_std_upper_bound - total_log2_std_before_poly_eval);
    assert!(max_log_interpolate_degree >= 0.);
    let max_interpolate_degree = (2.0 as f64).powf(max_log_interpolate_degree).floor() as usize;
    std::cmp::min(max_interpolate_degree as usize, poly_len as usize)

}

pub fn params_rgswpir_given_interpolate_degree(input_num_items: usize, interpolate_degree: usize, input_item_size_bits: usize) -> (Params, usize, (usize, usize, usize)) {

    let poly_len_log2 = 11;
    let poly_len = 1 << poly_len_log2; 
    let p = 12289;
    let log_p = (p as f64).log2().floor() as usize;

    let input_item_num_pts = (input_item_size_bits as f64 / log_p as f64).ceil() as usize; // -> number of pt modului required to represent one item

    let total_db_num_pts = (((input_num_items * input_item_size_bits) as f64) / log_p as f64).ceil() as usize; // -> total number of pt moduli required to represent the database
    // let sqrt_num_db_pts = (total_db_num_pts as f64).sqrt().floor() as usize; // -> sqrt of the previous line

    let mut c = ((total_db_num_pts * 56) as f64 / (interpolate_degree * poly_len * poly_len * (20 + 28)) as f64).sqrt().floor() as usize;
    if c == 0 {
        c = 1;
    }

    let col_size = c * interpolate_degree * poly_len;
    let new_item_size_num_pts = round_up_to_multiple_of(col_size, input_item_num_pts); // -> rounded up 
    let new_num_items = (total_db_num_pts as f64 / new_item_size_num_pts as f64).ceil() as usize;
    
    // needs to be at least 8 for code to work correctly (becuase of AlignedMemory)
    let db_rows = if new_num_items >= 8 {
        new_num_items
    } else {
        8
    };
    
    // number of required rlwe ciphertexts (in first layer)
    let db_cols_poly = (new_item_size_num_pts as f64 / (poly_len as f64)).ceil() as usize;

    let nu_1 = db_rows.next_power_of_two().trailing_zeros() as usize - poly_len_log2;
    debug!("chose nu_1: {}", nu_1);

    let q2_bits = 28; // modulus after packing
    let t_exp_left = 3;

    let mut params = internal_params_for(poly_len, nu_1, 0, p, q2_bits, t_exp_left, DEF_MOD_STR);
    params.instances = db_cols_poly;
    (params, interpolate_degree, (db_rows, 1, (db_cols_poly*poly_len)*log_p))
}

pub fn params_rgswpir_given_input_size_and_dim0(input_num_items: usize, input_item_size_bits: usize, dim0: usize) -> (Params, usize, (usize, usize, usize)) {

    let poly_len_log2 = 11;
    let poly_len = 1 << poly_len_log2; 
    // let p = 12289;
    // let p = 65537;
    let p = 65535;
    // let p = 61441;
    // let log_p = (p as f64).log2().floor() as usize;
    let log_p = 16;
    let q2_bits = 28; // modulus after packing
    let t_exp_left = 3;
    let modulus = (268369921u64*249561089u64) as f64;

    let max_interpolate_degree = max_interpolate_degree(modulus, dim0 as f64, p as f64, t_exp_left as f64, t_exp_left as f64, poly_len as f64) as usize;
    assert!(max_interpolate_degree > 1);

    let bits_per_poly = (poly_len * log_p) as f64;

    let mut factor = (bits_per_poly as f64 / input_item_size_bits as f64).floor() as usize;
    if factor == 0 {
        factor = 1;
    }
    let input_num_items = (input_num_items as f64 / factor as f64).ceil() as usize;
    let input_item_size_bits = input_item_size_bits * factor;

    let padded_item_size_num_bits = ((input_item_size_bits as f64) / bits_per_poly).ceil() as usize * bits_per_poly as usize; 
    let padded_item_num_pts = (padded_item_size_num_bits as f64 / log_p as f64).ceil() as usize; // -> number of pt modului required to represent one item

    let dim1_lower_bound = padded_item_num_pts * (input_num_items as f64 / dim0 as f64).ceil() as usize; 
    let mut current_dim1 = padded_item_num_pts;
    let mut interpolate_degree = 1;
    while current_dim1 < dim1_lower_bound {
        interpolate_degree *= 2;
        current_dim1 *= 2;
        if 2*interpolate_degree > max_interpolate_degree {
            break;
        }
    };

    let new_item_size_num_pts = round_up_to_multiple_of(dim1_lower_bound, current_dim1);

    // needs to be at least 8 for code to work correctly (becuase of AlignedMemory)
    let db_rows = if dim0 >= 8 {
        dim0
    } else {
        8
    };
    
    // number of required rlwe ciphertexts (in first layer)
    let db_cols_poly = (new_item_size_num_pts as f64 / (poly_len as f64)).ceil() as usize;

    let nu_1 = db_rows.next_power_of_two().trailing_zeros() as usize - poly_len_log2;

    let mut params = internal_params_for(poly_len, nu_1, 0, p, q2_bits, t_exp_left, DEF_MOD_STR);
    params.instances = db_cols_poly;
    (params, interpolate_degree, (db_rows, 1, new_item_size_num_pts*log_p))
}

pub fn params_rgswpir_given_input_size_and_dim0_small(input_num_items: usize, input_item_size_bits: usize, dim0: usize) -> (Params, usize, (usize, usize, usize)) {

    let poly_len_log2 = 10;
    let poly_len = 1 << poly_len_log2; 
    // let p = 12289;
    // let p = 17;
    let p = 65;
    // let p = 65537;
    let log_p = (p as f64).log2().floor() as usize;
    let q2_bits = 28; // modulus after packing
    let t_exp_left = 12;
    pub const CUSTOM_DEF_MOD_STR: &str = "[\"268369921\"]";
    let modulus = 268369921.;

    let max_interpolate_degree = max_interpolate_degree(modulus, dim0 as f64, p as f64, t_exp_left as f64, t_exp_left as f64, poly_len as f64) as usize;
    assert!(max_interpolate_degree > 1);

    let bits_per_poly = (poly_len * log_p) as f64;

    let mut factor = (bits_per_poly as f64 / input_item_size_bits as f64).floor() as usize;
    if factor == 0 {
        factor = 1;
    }
    let input_num_items = (input_num_items as f64 / factor as f64).ceil() as usize;
    let input_item_size_bits = input_item_size_bits * factor;


    let padded_item_size_num_bits = ((input_item_size_bits as f64) / bits_per_poly).ceil() as usize * bits_per_poly as usize; 
    let padded_item_num_pts = (padded_item_size_num_bits as f64 / log_p as f64).ceil() as usize; // -> number of pt modului required to represent one item

    let dim1_lower_bound = padded_item_num_pts * (input_num_items as f64 / dim0 as f64).ceil() as usize; 
    let mut current_dim1 = padded_item_num_pts;
    let mut interpolate_degree = 1;
    while current_dim1 < dim1_lower_bound {
        interpolate_degree *= 2;
        current_dim1 *= 2;
        if 2*interpolate_degree > max_interpolate_degree {
            break;
        }
    };

    let new_item_size_num_pts = round_up_to_multiple_of(dim1_lower_bound, current_dim1);

    // needs to be at least 8 for code to work correctly (becuase of AlignedMemory)
    let db_rows = if dim0 >= 8 {
        dim0
    } else {
        8
    };
    
    // number of required rlwe ciphertexts (in first layer)
    let db_cols_poly = (new_item_size_num_pts as f64 / (poly_len as f64)).ceil() as usize;

    let nu_1 = db_rows.next_power_of_two().trailing_zeros() as usize - poly_len_log2;


    let mut params = internal_params_for(poly_len, nu_1, 0, p, q2_bits, t_exp_left, CUSTOM_DEF_MOD_STR);
    params.instances = db_cols_poly;
    (params, interpolate_degree, (db_rows, 1, new_item_size_num_pts*log_p))
}



/// Run the YPIR scheme with the given parameters
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Number of items in the database
    #[clap(long)]
    num_items: Option<usize>,

    /// Interpolation degree
    #[clap(long)]
    dim0: Option<usize>,

    /// Interpolation degree
    #[clap(long)]
    interpolate_degree: Option<usize>,

    #[clap(long)]
    item_size_bits: Option<usize>,

    /// Small params mode (optional)
    /// if set, small paramseters will be used
    #[clap(long, short, action)]
    small_params: bool,

    #[clap(long)]
    trials: Option<usize>,

    /// Output report file (optional)
    /// where results will be written in JSON.
    #[clap(long)]
    out_report_json: Option<String>,

    /// Just a label
    #[clap(long)]
    label: Option<String>,

    /// Verbose mode (optional)
    /// if set, the program will print debug logs to stderr.
    #[clap(long, short, action)]
    verbose: bool,

    #[clap(long, short, action)]
    online_only: bool,

}

fn main() {
    let args = Args::parse();
    let Args {
        num_items,
        dim0,
        interpolate_degree,
        item_size_bits,
        small_params,
        trials,
        out_report_json,
        label,
        verbose,
        online_only,
    } = args;

    if verbose {
        println!("Running in verbose mode.");
        env_logger::Builder::new()
            .filter_level(log::LevelFilter::Debug)
            .write_style(env_logger::WriteStyle::Always)
            .init();
    } else {
        env_logger::init();
    }

    let num_items = num_items.unwrap_or(4096);
    let dim0 = dim0.unwrap_or(2048);
    let interpolate_degree = interpolate_degree.unwrap_or(0);
    let item_size_bits = item_size_bits.unwrap_or(4096 * 16);// 262144);
    let small_params = small_params;
    let trials = trials.unwrap_or(1);
    let label = label.unwrap_or("".to_string());

    println!(
        "Protocol=InsPIRe, DB={} KB, trials={}",
        (num_items * item_size_bits) / 8192,
        trials
    );

    let given_d0 = interpolate_degree == 0;
    let (params, interpolate_degree, (resized_db_first_dim, resized_db_second_dim, resized_item_size_bits)) = if given_d0 {
        if small_params {
            params_rgswpir_given_input_size_and_dim0_small(num_items, item_size_bits, dim0)
        } else {
            params_rgswpir_given_input_size_and_dim0(num_items, item_size_bits, dim0)
        }
    } else {
        params_rgswpir_given_interpolate_degree(num_items, interpolate_degree, item_size_bits)
    };

    println!("DB Size={} MB", (resized_db_first_dim * resized_db_second_dim * resized_item_size_bits) as f64 / (8. * 1024. * 1024.));

    let mut measurement = run_simple_ypir_rgsw_on_params(params, interpolate_degree, trials, online_only);
        
    measurement.specs.resized_db_first_dim = resized_db_first_dim;
    measurement.specs.resized_db_second_dim = resized_db_second_dim;
    measurement.specs.resized_item_size_bits = resized_item_size_bits;
    measurement.specs.resized_database_size_mb = (resized_db_first_dim * resized_db_second_dim * resized_item_size_bits) as f64 / (8. * 1024. * 1024.);
    measurement.specs.online_only = online_only;

    measurement.specs.input_num_items = num_items;
    measurement.specs.input_item_size_bits = item_size_bits;
    measurement.specs.input_database_size_mb = (num_items * item_size_bits) as f64 / (8. * 1024. * 1024.);
    measurement.specs.interpolate_degree = interpolate_degree;
    measurement.specs.label = label;

    if let Some(out_report_json) = out_report_json {
        println!("Writing report to {}", out_report_json);
        let mut file = std::fs::File::create(out_report_json).unwrap();
        serde_json::to_writer_pretty(&mut file, &measurement).unwrap();
        println!("Report written.");
    } else {
        println!(
            "Measurement completed. See the README for details on what the
            following fields mean."
        );
        println!("Result:");
        println!("{}", serde_json::to_string_pretty(&measurement).unwrap());    
    }
}
