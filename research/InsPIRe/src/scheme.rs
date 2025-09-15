use std::collections::HashMap;
use std::time::Instant;

use log::debug;
use rand::{thread_rng, Rng};

use spiral_rs::aligned_memory::AlignedMemory64;
use spiral_rs::arith::rescale;
use spiral_rs::poly::{PolyMatrix, PolyMatrixRaw};
use spiral_rs::{client::*, params::*};

use crate::bits::{read_bits, u64s_to_contiguous_bytes};
use crate::modulus_switch::ModulusSwitch;
use crate::noise_analysis::YPIRSchemeParams;
use crate::packing::{PackingKeys, PackingType, ToStr};

use super::{client::*, lwe::LWEParams, measurement::*, params::*, server::*};
use strum_macros::{Display, EnumString};

pub const STATIC_PUBLIC_SEED: [u8; 32] = [0u8; 32];
pub const SEED_0: u8 = 0;
pub const SEED_1: u8 = 1;

pub const STATIC_SEED_2: [u8; 32] = [
    2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

pub const W_SEED: [u8; 32] = [
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
];
pub const V_SEED: [u8; 32] = [
    8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
];

// have three protocol types called "simplepir", "doublepir", and "InsPIRe"
// implement a display function too
#[derive(Debug, Clone, Copy, PartialEq, EnumString, Display)]
pub enum ProtocolType {
    #[strum(serialize = "SimplePIR")]
    SimplePIR,
    #[strum(serialize = "DoublePIR")]
    DoublePIR,
    #[strum(serialize = "InsPIRe")]
    InsPIRe,
}

impl Default for ProtocolType {
    fn default() -> Self {
        ProtocolType::DoublePIR // Specify the default variant
    }
}

impl ToStr for ProtocolType {
    fn to_str(&self) -> String {
        match self {
            ProtocolType::SimplePIR => "SimplePIR",
            ProtocolType::DoublePIR => "DoublePIR",
            ProtocolType::InsPIRe => "InsPIRe",            
        }.to_string()
    } 
}

pub fn run_ypir_batched(
    num_items: usize,
    item_size_bits: usize,
    protocol_type: ProtocolType,
    second_level_packing_mask: PackingType,
    second_level_packing_body: PackingType,
    performance_factor: usize,
    small_params: bool,
    gammas: Vec<usize>,
    trials: usize,
    online_only: bool
) -> Measurement {
    let (params, (resized_db_first_dim, resized_db_second_dim, resized_item_size_bits)) = match protocol_type {
        ProtocolType::SimplePIR => {
            if small_params {
                params_for_scenario_simplepir_small(num_items, item_size_bits)
            } else {
                params_for_scenario_simplepir(num_items, item_size_bits)
            }
        },
        ProtocolType::DoublePIR => params_for_scenario(num_items, item_size_bits),
        ProtocolType::InsPIRe => params_for_scenario_medium_payload(num_items, item_size_bits, gammas.clone(), performance_factor),
    };

    // println!("size: {} MB", (resized_db_first_dim * resized_db_second_dim * resized_item_size_bits) as f64 / (8. * 1024. * 1024.));

    let mut measurement = match protocol_type {
        ProtocolType::SimplePIR => {
            run_simple_ypir_on_params(params, second_level_packing_mask, gammas[0], trials, online_only)
        }
        ProtocolType::DoublePIR => {
            run_ypir_on_params(params, second_level_packing_mask, second_level_packing_body, gammas[0], trials, online_only)
        }
        ProtocolType::InsPIRe => run_ypir_with_medium_payload_on_params(
            params,
            second_level_packing_mask,
            second_level_packing_body,
            gammas,
            trials,
            online_only
        ),
    };
        
    measurement.specs.resized_db_first_dim = resized_db_first_dim;
    measurement.specs.resized_db_second_dim = resized_db_second_dim;
    measurement.specs.resized_item_size_bits = resized_item_size_bits;
    measurement.specs.resized_database_size_mb = (resized_db_first_dim * resized_db_second_dim * resized_item_size_bits) as f64 / (8. * 1024. * 1024.);
    measurement.specs.online_only = online_only;
    debug!("{:#?}", measurement);
    measurement
}

pub trait Sample {
    fn sample() -> Self;
}

impl Sample for u8 {
    fn sample() -> Self {
        fastrand::u8(..)
    }
}

impl Sample for u16 {
    fn sample() -> Self {
        fastrand::u16(..)
    }
}

impl Sample for u32 {
    fn sample() -> Self {
        fastrand::u32(..)
    }
}

pub fn run_simple_ypir_on_params(
    params: Params,
    packing_type: PackingType,
    gamma: usize,
    trials: usize,
    online_only: bool
) -> Measurement {
    println!("Starting... ");
    // let is_simplepir = true;
    let db_rows = 1 << (params.db_dim_1 + params.poly_len_log2);
    let db_cols = params.instances * params.poly_len;

    let mut rng = thread_rng();

    // RLWE reduced moduli
    let rlwe_q_prime_1 = params.get_q_prime_1();
    let rlwe_q_prime_2 = params.get_q_prime_2();

    let db_cols_prime = db_cols / gamma;

    let now = Instant::now();
    type T = u16;
    let pt_iter = std::iter::repeat_with(|| (T::sample() as u64 % params.pt_modulus) as T);
    // let pt_iter = std::iter::repeat_with(|| (5 as u64 % params.pt_modulus) as T);
    let y_server = YServer::<T>::new(
        &params,
        pt_iter,
        ProtocolType::SimplePIR,
        packing_type,
        PackingType::NoPacking,
        false,
        vec![gamma],
    );

    debug!("Created server in {} us", now.elapsed().as_micros());
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
    println!("Starting Offline... ");

    let start_offline_comp = Instant::now();
    let offline_values = y_server
        .perform_offline_precomputation_simplepir(gamma, Some(&mut measurements[0]), online_only);
    let offline_server_time_ms = start_offline_comp.elapsed().as_millis();
    println!("Starting after online... ");

    let packed_query_row_sz = db_rows;
    // let mut all_queries_packed = AlignedMemory64::new(K * packed_query_row_sz);

    for trial in 0..trials + 1 {
        debug!("trial: {}", trial);
        let mut measurement = &mut measurements[trial];
        measurement.offline.server_time_ms = offline_server_time_ms as usize;

        // ================================================================
        // QUERY GENERATION PHASE
        // ================================================================
        let mut queries = Vec::new();

        let mut client = Client::init(&params);

        let target_idx: usize = rng.gen::<usize>() % (db_rows * db_cols);
        let target_row = target_idx / db_cols;
        let target_col = target_idx % db_cols;
        debug!("Target item: {} ({}, {})", target_idx, target_row, target_col);

        let start = Instant::now();
        client.generate_secret_keys();
        let sk_reg = client.get_sk_reg();

        let mut packing_keys = match packing_type {
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
        let pub_params_size = packing_keys.get_size_bytes();                
        debug!("pub params size: {} bytes", pub_params_size);

        let y_client = YClient::new(&mut client, &params);
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
            ((packed_query_row.len() as f64 * params.modulus_log2 as f64) / 8.0).ceil() as usize;

        measurement.online.client_query_gen_time_ms = start.elapsed().as_millis() as usize;
        debug!("Generated query in {} us", start.elapsed().as_micros());

        let online_upload_bytes = query_size + pub_params_size;
        debug!("Query size: {} bytes", online_upload_bytes);

        queries.push((y_client, target_idx, packed_query_row));

        let mut all_queries_packed = AlignedMemory64::new(packed_query_row_sz);
        for (i, chunk_mut) in
            all_queries_packed.as_mut_slice().chunks_mut(packed_query_row_sz).enumerate()
        {
            (&mut chunk_mut[..db_rows]).copy_from_slice(queries[i].2.as_slice());
        }

        // ================================================================
        // ONLINE PHASE
        // ================================================================

        println!("Starting Online... ");

        let start_online_comp = Instant::now();

        let responses = vec![y_server.perform_online_computation_simplepir(
            gamma,
            all_queries_packed.as_slice(),
            &offline_values,
            &mut packing_keys,
            Some(&mut measurement),
        )];
        let online_server_time_ms = start_online_comp.elapsed().as_millis();
        let online_download_bytes = get_size_bytes(&responses); // TODO: this is not quite right for multiple clients

        // check correctness
        for (response_switched, (y_client, target_idx, _)) in
            responses.iter().zip(queries.iter())
        {
            let (target_row, _target_col) = (target_idx / db_cols, target_idx % db_cols);
            let corr_result =
                y_server.get_row(target_row).iter().map(|x| x.to_u64()).collect::<Vec<_>>();

            // let scheme_params = YPIRSchemeParams::from_params(&params, &lwe_params);
            // let log2_corr_err = scheme_params.delta().log2();
            // let log2_expected_outer_noise = scheme_params.expected_outer_noise().log2();
            // debug!("log2_correctness_err: {}", log2_corr_err);
            // debug!("log2_expected_outer_noise: {}", log2_expected_outer_noise);

            let start_decode = Instant::now();

            debug!("rescaling response...");
            let mut response = Vec::new();
            for ct_bytes in response_switched.iter() {
                let ct = PolyMatrixRaw::recover_how_many(&params, rlwe_q_prime_1, rlwe_q_prime_2, gamma, ct_bytes);
                response.push(ct);
            }

            debug!("decrypting outer cts...");
            let outer_ct = response
                .iter()
                .flat_map(|ct| {
                    decrypt_ct_reg_measured(y_client.client(), &params, &ct.ntt(), params.poly_len)
                        .as_slice()[..gamma]
                        .to_vec()
                })
                .collect::<Vec<_>>();
            assert_eq!(outer_ct.len(), db_cols_prime * gamma);
            // debug!("outer_ct: {:?}", &outer_ct[..]);
            // let outer_ct_t_u8 = u64s_to_contiguous_bytes(&outer_ct, pt_bits);
            let final_result = outer_ct.as_slice();
            measurement.online.client_decode_time_ms = start_decode.elapsed().as_millis() as usize;

            // debug!("got      {:?}", &final_result[..256]);
            // debug!("expected {:?}", &corr_result[..256]);
            // debug!("got {:?}, expected {:?}", &final_result[..256], &corr_result[..256]);
            if online_only {
                println!("Warning!")
            } else {
                assert_eq!(final_result, corr_result);
            }
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
    final_measurement.specs.second_level_packing_mask = packing_type.to_str();
    final_measurement.specs.second_level_packing_body = PackingType::NoPacking.to_str();
    final_measurement.specs.poly_len = params.poly_len;
    final_measurement.specs.modulus_bits = params.modulus_log2 as usize;
    final_measurement.specs.pt_modulus = params.pt_modulus as usize;

    final_measurement.specs.gamma_0 = gamma;

    final_measurement
}

pub fn run_ypir_on_params(
    params: Params,
    second_level_packing_mask: PackingType,
    second_level_packing_body: PackingType,
    gamma: usize,
    trials: usize,
    online_only: bool
) -> Measurement {

    let lwe_params = LWEParams::default();

    // let gamma = params.poly_len;

    let db_rows = 1 << (params.db_dim_1 + params.poly_len_log2);
    let db_cols = 1 << (params.db_dim_2 + params.poly_len_log2);

    let sqrt_n_bytes = db_cols * (lwe_params.pt_modulus as f64).log2().floor() as usize / 8;

    let mut rng = thread_rng();

    // RLWE reduced moduli
    let rlwe_q_prime_1 = params.get_q_prime_1();
    let rlwe_q_prime_2 = params.get_q_prime_2();

    // LWE reduced moduli
    let lwe_q_prime_bits = lwe_params.q2_bits as usize;

    // The number of bits represented by a plaintext RLWE coefficient
    let pt_bits = (params.pt_modulus as f64).log2().floor() as usize;
    // assert_eq!(pt_bits, 16);

    // The factor by which ciphertext values are bigger than plaintext values
    let blowup_factor = lwe_q_prime_bits as f64 / pt_bits as f64;
    debug!("blowup_factor: {}", blowup_factor);

    let mut smaller_params = params.clone();
    smaller_params.db_dim_1 = params.db_dim_2;
    smaller_params.db_dim_2 = ((blowup_factor * (lwe_params.n + 1) as f64) / params.poly_len as f64)
        .log2()
        .ceil() as usize;

    let out_rows = 1 << (smaller_params.db_dim_2 + params.poly_len_log2);
    let rho: usize = 1 << smaller_params.db_dim_2; // rho

    debug!("rho: {}", rho);

    assert_eq!(smaller_params.db_dim_1, params.db_dim_2);
    assert!(out_rows as f64 >= (blowup_factor * (lwe_params.n + 1) as f64));

    // --

    let lwe_q_bits = (lwe_params.modulus as f64).log2().ceil() as usize;

    let rlwe_q_prime_1_bits = (rlwe_q_prime_1 as f64).log2().ceil() as usize;
    let rlwe_q_prime_2_bits = (rlwe_q_prime_2 as f64).log2().ceil() as usize;
    let simplepir_hint_bytes = (lwe_params.n * db_cols * lwe_q_prime_bits) / 8;
    let doublepir_hint_bytes = (params.poly_len * out_rows * rlwe_q_prime_2_bits) / 8;
    let simplepir_query_bytes = db_rows * lwe_q_bits / 8;
    let doublepir_query_bytes = db_cols * params.modulus_log2 as usize / 8;
    let simplepir_resp_bytes = (db_cols * lwe_q_prime_bits) / 8;
    let doublepir_resp_bytes = ((rho * params.poly_len) * rlwe_q_prime_2_bits
        + (rho * params.poly_len) * rlwe_q_prime_1_bits)
        / 8;
    debug!("          \"simplepirHintBytes\": {},", simplepir_hint_bytes);
    debug!("          \"doublepirHintBytes\": {}", doublepir_hint_bytes);
    debug!("          \"simplepirQueryBytes\": {},", simplepir_query_bytes);
    debug!("          \"doublepirQueryBytes\": {},", doublepir_query_bytes);
    debug!("          \"simplepirRespBytes\": {},", simplepir_resp_bytes);
    debug!("          \"doublepirRespBytes\": {},", doublepir_resp_bytes);

    // --

    let now = Instant::now();
    let pt_iter = std::iter::repeat_with(|| u8::sample());
    let y_server = YServer::<u8>::new(
        &params,
        pt_iter,
        ProtocolType::DoublePIR,
        second_level_packing_mask,
        second_level_packing_body,
        false,
        vec![gamma],
    );
    debug!("Created server in {} us", now.elapsed().as_micros());
    debug!("Database of {} bytes", y_server.db().len() * std::mem::size_of::<u8>());
    // let db_pt_modulus = lwe_params.pt_modulus;
    // assert_eq!(
    //     y_server.db().len() * std::mem::size_of::<u8>(),
    //     db_rows * db_cols * (db_pt_modulus as f64).log2().ceil() as usize / 8
    // );

    // ================================================================
    // OFFLINE PHASE
    // ================================================================
    let mut measurements = vec![Measurement::default(); trials + 1];

    let start_offline_comp = Instant::now();
    let offline_values = y_server.perform_offline_precomputation(
        gamma,
        Some(&mut measurements[0]),
        online_only,
    );
    let offline_server_time_ms = start_offline_comp.elapsed().as_millis();

    let packed_query_row_sz = db_rows;
    // let mut all_queries_packed = AlignedMemory64::new(K * packed_query_row_sz);

    for trial in 0..trials + 1 {
        debug!("trial: {}", trial);
        let mut measurement = &mut measurements[trial];
        measurement.offline.server_time_ms = offline_server_time_ms as usize;
        measurement.offline.simplepir_hint_bytes = simplepir_hint_bytes;
        measurement.offline.doublepir_hint_bytes = doublepir_hint_bytes;
        measurement.online.simplepir_query_bytes = simplepir_query_bytes;
        measurement.online.doublepir_query_bytes = doublepir_query_bytes;
        measurement.online.simplepir_resp_bytes = simplepir_resp_bytes;
        measurement.online.doublepir_resp_bytes = doublepir_resp_bytes;

        // ================================================================
        // QUERY GENERATION PHASE
        // ================================================================
        let online_upload_bytes: usize;
        let mut queries = Vec::new();

        let mut client = Client::init(&params);

        let target_idx: usize = rng.gen::<usize>() % (db_rows * db_cols);
        let target_row = target_idx / db_cols;
        let target_col = target_idx % db_cols;
        debug!("Target item: {} ({}, {})", target_idx, target_row, target_col);

        let start = Instant::now();
        client.generate_secret_keys();
        let sk_reg = &client.get_sk_reg();

        let mut packing_keys = match second_level_packing_mask {
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
        let pub_params_size = packing_keys.get_size_bytes();                
        debug!("pub params size: {} bytes", pub_params_size);

        let y_client = YClient::new(&mut client, &params);
        let query_row = y_client.generate_query_over_u32(SEED_0, params.db_dim_1, target_row);
        let query_row_last_row: &[u64] = &query_row[lwe_params.n * db_rows..];
        let mut aligned_query_packed = AlignedMemory64::new(query_row_last_row.len());
        aligned_query_packed.as_mut_slice().copy_from_slice(&query_row_last_row);
        let packed_query_row = aligned_query_packed;
        let packed_query_row_u32 =
            packed_query_row.as_slice().iter().map(|x| *x as u32).collect::<Vec<_>>();

        let query_col = y_client.generate_query_over_prime(
            SEED_1,
            params.db_dim_2,
            second_level_packing_mask,
            target_col,
        );
        let query_col_last_row = &query_col[params.poly_len * db_cols..];
        let packed_query_col = pack_query(&params, query_col_last_row);

        let query_size = query_row_last_row.len() * 4 + query_col_last_row.len() * 8;

        measurement.online.client_query_gen_time_ms = start.elapsed().as_millis() as usize;
        debug!("Generated query in {} us", start.elapsed().as_micros());

        online_upload_bytes = query_size + pub_params_size;
        debug!("Query size: {} bytes", online_upload_bytes);

        queries.push((
            y_client,
            target_idx,
            packed_query_row_u32,
            packed_query_col,
        ));

        let mut all_queries_packed = vec![0u32; packed_query_row_sz];
        for (i, chunk_mut) in
            all_queries_packed.as_mut_slice().chunks_mut(packed_query_row_sz).enumerate()
        {
            (&mut chunk_mut[..db_rows]).copy_from_slice(queries[i].2.as_slice());
        }

        let mut offline_values = offline_values.clone();

        // ================================================================
        // ONLINE PHASE
        // ================================================================

        let start_online_comp = Instant::now();
        let responses = y_server.perform_online_computation(
            gamma,
            &mut packing_keys,
            &mut offline_values,
            &all_queries_packed,
            &queries.iter().map(|x| x.3.as_slice()).collect::<Vec<_>>(),
            Some(&mut measurement),
        );
        let online_server_time_ms = start_online_comp.elapsed().as_millis();
        let online_download_bytes = get_size_bytes(&responses);
        // check correctness
        for (response_switched, (y_client, target_idx, _, _)) in
            responses.iter().zip(queries.iter())
        {
            let (target_row, target_col) = (target_idx / db_cols, target_idx % db_cols);
            let corr_result = y_server.get_elem(target_row, target_col).to_u64();

            let scheme_params = YPIRSchemeParams::from_params(&params, &lwe_params);
            let log2_corr_err = scheme_params.delta().log2();
            let log2_expected_outer_noise = scheme_params.expected_outer_noise().log2();
            debug!("log2_correctness_err: {}", log2_corr_err);
            debug!("log2_expected_outer_noise: {}", log2_expected_outer_noise);

            let start_decode = Instant::now();

            debug!("rescaling response...");
            let mut response = Vec::new();
            for ct_bytes in response_switched.iter() {
                let ct = PolyMatrixRaw::recover_how_many(&params, rlwe_q_prime_1, rlwe_q_prime_2, gamma, ct_bytes);
                response.push(ct);
            }

            let s_2 = ct_reg_measure(y_client.client(), &params, &response[0].ntt(), gamma);
            measurement.online.noise_width = s_2.log2() / 2.0;
            // println!("log2(measured noise): {}", s_2.log2());

            debug!("decrypting outer cts...");
            let outer_ct = response
                .iter()
                .flat_map(|ct| {
                    decrypt_ct_reg_measured(y_client.client(), &params, &ct.ntt(), gamma)
                        .as_slice()[..gamma]
                        .to_vec()
                })
                .collect::<Vec<_>>();
            // assert_eq!(outer_ct.len(), response.len() * out_rows);
            // debug!("outer_ct: {:?}", &outer_ct[..]);
            let outer_ct_t_u8 = u64s_to_contiguous_bytes(&outer_ct, pt_bits);

            let mut inner_ct = PolyMatrixRaw::zero(&params, 2, 1);
            let mut bit_offs = 0;
            let lwe_q_prime = lwe_params.get_q_prime_2();
            let special_offs =
                ((lwe_params.n * lwe_q_prime_bits) as f64 / pt_bits as f64).ceil() as usize;
            for z in 0..lwe_params.n {
                let val = read_bits(&outer_ct_t_u8, bit_offs, lwe_q_prime_bits);
                bit_offs += lwe_q_prime_bits;
                if !online_only {
                    assert!(val < lwe_q_prime, "val: {}, lwe_q_prime: {}", val, lwe_q_prime);
                }
                inner_ct.data[z] = rescale(val, lwe_q_prime, lwe_params.modulus);
            }

            let num_rlwes_for_mask = (special_offs as f64 / gamma as f64).ceil() as usize;
            let mut val = 0;
            for i in 0..blowup_factor.ceil() as usize {
                match second_level_packing_body {
                    PackingType::InspiRING => {
                        val |= (outer_ct[params.poly_len + i] % params.pt_modulus) << (i * pt_bits);
                    }
                    PackingType::CDKS => {
                        val |= outer_ct[special_offs + i] << (i * pt_bits);
                    }
                    PackingType::NoPacking => {
                        val |= (outer_ct[num_rlwes_for_mask * gamma +  i * gamma] % params.pt_modulus)
                            << (i * pt_bits);
                    }
                }
            }
            if !online_only {
                assert!(val < lwe_q_prime, "val: {}, lwe_q_prime: {}", val, lwe_q_prime);
            }
            debug!("got b_val of: {}", val);
            inner_ct.data[lwe_params.n] = rescale(val, lwe_q_prime, lwe_params.modulus);

            debug!("decrypting inner ct...");
            // let plaintext = decrypt_ct_reg_measured(y_client.client(), &params,
            // &inner_ct.ntt(), 1); let final_result = plaintext.data[0];
            let inner_ct_as_u32 = inner_ct
                .as_slice()
                .iter()
                .take(lwe_params.n + 1)
                .map(|x| *x as u32)
                .collect::<Vec<_>>();
            let decrypted = y_client.lwe_client().decrypt(&inner_ct_as_u32);
            let final_result = rescale(decrypted as u64, lwe_params.modulus, lwe_params.pt_modulus);

            measurement.online.client_decode_time_ms = start_decode.elapsed().as_millis() as usize;

            debug!("got {}, expected {}", final_result, corr_result);
            // debug!("was correct? {}", final_result == corr_result);
            if !online_only {
                assert_eq!(final_result, corr_result);
            }
        }

        measurement.online.upload_query = query_size;
        measurement.online.upload_keys = pub_params_size;
        measurement.online.upload_bytes = online_upload_bytes;
        measurement.online.download_bytes = online_download_bytes;
        measurement.online.total_bytes = online_upload_bytes + online_download_bytes;
        measurement.online.server_time_ms = online_server_time_ms as usize;
        measurement.online.sqrt_n_bytes = sqrt_n_bytes;
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
    
    final_measurement.online.noise_width =
        mean_f64(&measurements.iter().map(|m| m.online.noise_width).collect::<Vec<_>>());    
    let all_noise_width =
        measurements.iter().map(|m| m.online.noise_width).collect::<Vec<_>>();
    final_measurement.online.std_dev_server_time_ms =
        std_dev_f64(&all_noise_width);

    final_measurement.specs.protocol_type = ProtocolType::DoublePIR.to_str();
    final_measurement.specs.second_level_packing_mask = second_level_packing_mask.to_str();
    final_measurement.specs.second_level_packing_body = second_level_packing_body.to_str();
    final_measurement.specs.poly_len = params.poly_len;
    final_measurement.specs.modulus_bits = params.modulus_log2 as usize;
    final_measurement.specs.pt_modulus = params.pt_modulus as usize;

    final_measurement.specs.gamma_0 = gamma;

    final_measurement
}

pub fn run_ypir_with_medium_payload_on_params(
    params: Params,
    second_level_packing_mask: PackingType,
    second_level_packing_body: PackingType,
    gammas: Vec<usize>,
    trials: usize,
    online_only: bool
) -> Measurement {

    if second_level_packing_body == PackingType::InspiRING {
        assert_eq!(gammas.len(), 3);
    } else {
        assert!(gammas.len() >= 2);
    }

    let gamma_first_layer = gammas[0];
    let gamma_second_layer_mask = gammas[1];
    let gamma_second_layer_body = if second_level_packing_body == PackingType::InspiRING {
        gammas[2]
    } else { 1 };

    let db_rows = 1 << (params.db_dim_1 + params.poly_len_log2);
    let db_cols = params.instances * params.poly_len;

    let db_cols_prime = db_cols / gamma_first_layer;

    let mut rng = thread_rng();

    // RLWE reduced moduli
    let rlwe_q_prime_1 = params.get_q_prime_1();
    let rlwe_q_prime_2 = params.get_q_prime_2();
    // let rlwe_q_prime_1_bits = (rlwe_q_prime_1 as f64).log2().ceil() as usize;
    // let rlwe_q_prime_2_bits = (rlwe_q_prime_2 as f64).log2().ceil() as usize;

    // LWE reduced moduli
    let q_prime_bits = params.q2_bits as usize;

    // The number of bits represented by a plaintext RLWE coefficient
    let pt_bits = (params.pt_modulus as f64).log2().floor() as usize;

    // The factor by which ciphertext values are bigger than plaintext values
    let blowup_factor = q_prime_bits as f64 / pt_bits as f64;
    let mask_decomposed_num_coeffs = (params.poly_len as f64 * blowup_factor).ceil() as usize;
    let body_decomposed_num_coeffs = (gamma_first_layer as f64 * blowup_factor).ceil() as usize;

    // assert_eq!(mask_decomposed_num_coeffs % gamma, 0);
    let mask_num_packed_rlwes = (mask_decomposed_num_coeffs as f64 / gamma_second_layer_mask as f64).ceil() as usize;

    // let simplepir_hint_bytes = (params.poly_len * db_cols * q_prime_bits) / 8;
    // let doublepir_hint_bytes = (params.poly_len * out_rows * rlwe_q_prime_2_bits) / 8;
    let simplepir_query_bytes = db_rows * params.modulus_log2 as usize / 8;
    let doublepir_query_bytes = db_cols_prime * params.modulus_log2 as usize / 8;
    // let simplepir_resp_bytes = (db_cols * q_prime_bits) / 8;
    // let doublepir_resp_bytes = (mask_num_packed_rlwes + body_decomposed_num_coeffs)
    //                         * params.poly_len
    //                         * (rlwe_q_prime_2_bits + rlwe_q_prime_1_bits) / 8;
    // ((rho * params.poly_len) * rlwe_q_prime_2_bits
    //     + (rho * params.poly_len) * rlwe_q_prime_1_bits) / 8;

    // debug!("          \"simplepirHintBytes\": {},", simplepir_hint_bytes);
    // debug!("          \"doublepirHintBytes\": {}", doublepir_hint_bytes);
    debug!("          \"simplepirQueryBytes\": {},", simplepir_query_bytes);
    debug!("          \"doublepirQueryBytes\": {},", doublepir_query_bytes);
    // debug!("          \"simplepirRespBytes\": {},", simplepir_resp_bytes);
    // debug!("          \"doublepirRespBytes\": {},", doublepir_resp_bytes);

    // --

    let now = Instant::now();
    // let pt_iter = std::iter::repeat_with(|| u8::sample());
    type T = u16;
    // assert!(params.pt_modulus <= 256);
    let pt_iter = std::iter::repeat_with(|| (T::sample() as u64 % params.pt_modulus) as T);

    let y_server = YServer::<T>::new(
        &params,
        pt_iter,
        ProtocolType::InsPIRe,
        second_level_packing_mask,
        second_level_packing_body,
        false,
        gammas.clone(),
    );
    debug!("Created server in {} us", now.elapsed().as_micros());
    debug!("Database of {} bytes", y_server.db().len() * std::mem::size_of::<u8>());
    // let db_pt_modulus = params.pt_modulus;
    // assert_eq!(
    //     y_server.db().len() * std::mem::size_of::<u8>(),
    //     db_rows * db_cols * (db_pt_modulus as f64).log2().ceil() as usize / 8
    // );

    // ================================================================
    // OFFLINE PHASE
    // ================================================================
    let mut measurements = vec![Measurement::default(); trials + 1];

    let start_offline_comp = Instant::now();
    let offline_values =
        y_server.perform_offline_precomputation_medium_payload(&gammas, Some(&mut measurements[0]), online_only);
    let offline_server_time_ms = start_offline_comp.elapsed().as_millis();

    debug!("Done offline!");

    let packed_query_row_sz = db_rows;
    // let mut all_queries_packed = AlignedMemory64::new(K * packed_query_row_sz);

    for trial in 0..trials + 1 {
        debug!("trial: {}", trial);
        let mut measurement = &mut measurements[trial];
        measurement.offline.server_time_ms = offline_server_time_ms as usize;
        // measurement.offline.simplepir_hint_bytes = simplepir_hint_bytes;
        // measurement.offline.doublepir_hint_bytes = doublepir_hint_bytes;
        measurement.online.simplepir_query_bytes = simplepir_query_bytes;
        measurement.online.doublepir_query_bytes = doublepir_query_bytes;
        // measurement.online.simplepir_resp_bytes = simplepir_resp_bytes;
        // measurement.online.doublepir_resp_bytes = doublepir_resp_bytes;

        // ================================================================
        // QUERY GENERATION PHASE
        // ================================================================
        let mut queries = Vec::new();

        let mut client = Client::init(&params);

        let target_idx: usize = rng.gen::<usize>() % (db_rows * db_cols);
        let target_row = target_idx / db_cols;
        let target_col = target_idx % db_cols;
        let target_rlwe_col = target_col / gamma_first_layer;
        debug!("Target item: {} ({}, {})", target_idx, target_row, target_col);

        let start = Instant::now();
        client.generate_secret_keys();
        let sk_reg = &client.get_sk_reg();

        let mut packing_keys_set : HashMap<usize, PackingKeys> = HashMap::new();
        for gamma in gammas.clone() {
            if !packing_keys_set.contains_key(&gamma) {
                let packing_keys = if gamma <= params.poly_len / 2 { 
                    PackingKeys::init(&y_server.packing_params_set[&gamma], sk_reg, W_SEED)
                } else {
                    PackingKeys::init_full(&y_server.packing_params_set[&gamma], sk_reg, W_SEED, V_SEED)
                };
                packing_keys_set.insert(gamma, packing_keys);
            }
        }
        let mut pub_params_size = 0;
        for packing_key in packing_keys_set.values() {
            pub_params_size += packing_key.get_size_bytes();
        }

        let y_client = YClient::new(&mut client, &params);
        let query_row = y_client.generate_query_over_prime(
            SEED_0,
            params.db_dim_1,
            PackingType::InspiRING,
            target_row,
        );

        assert_eq!(query_row.len(), (params.poly_len + 1) * db_rows);
        let query_row_last_row: &[u64] = &query_row[params.poly_len * db_rows..];

        assert_eq!(query_row_last_row.len(), db_rows);
        let packed_query_row = pack_query(&params, query_row_last_row);

        let query_col = y_client.generate_query_over_prime(
            SEED_1,
            y_server.get_smaller_params().db_dim_1,
            PackingType::InspiRING,
            target_rlwe_col,
        );
        assert_eq!(query_col.len(), (params.poly_len + 1) * db_cols_prime);
        let query_col_last_row = &query_col[params.poly_len * db_cols_prime..];
        debug!("query_col_last_row len:{}", query_col_last_row.len());
        let packed_query_col = pack_query(&params, query_col_last_row);

        let query_size =
            (query_row_last_row.len() * params.modulus_log2 as usize + query_col_last_row.len() * params.modulus_log2 as usize) / 8;

        measurement.online.client_query_gen_time_ms = start.elapsed().as_millis() as usize;
        debug!("Generated query in {} us", start.elapsed().as_micros());

        let online_upload_bytes = query_size + pub_params_size;
        debug!("Query size: {} bytes", online_upload_bytes);

        queries.push((y_client, target_idx, packed_query_row, packed_query_col));

        let mut all_queries_packed = AlignedMemory64::new(packed_query_row_sz);
        for (i, chunk_mut) in
            all_queries_packed.as_mut_slice().chunks_mut(packed_query_row_sz).enumerate()
        {
            (&mut chunk_mut[..db_rows]).copy_from_slice(queries[i].2.as_slice());
        }

        let mut offline_values = offline_values.clone();

        // ================================================================
        // ONLINE PHASE
        // ================================================================

        let start_online_comp = Instant::now();
        let responses = vec![y_server.perform_online_computation_medium_payload(
            &mut offline_values,
            all_queries_packed.as_slice(),
            &queries.iter().map(|x| x.3.as_slice()).collect::<Vec<_>>(),
            packing_keys_set,
            gammas.clone(),
            Some(&mut measurement),
        )];
        let online_server_time_ms = start_online_comp.elapsed().as_millis();
        // let online_download_bytes = get_size_bytes(&responses);
        let online_download_bytes = responses[0].get_size_bytes();
        // check correctness
        for (response_object, (y_client, target_idx, _, _)) in
            responses.iter().zip(queries.iter())
        {

            let (target_row, target_col) = (target_idx / db_cols, target_idx % db_cols);
            let corr_result = y_server.get_elem(target_row, target_col).to_u64();

            let start_decode = Instant::now();

            debug!("rescaling response...");
            let mut response_mask = Vec::new();
            for ct_bytes in response_object.packed_mask_mod_switched.iter() {
                let ct = PolyMatrixRaw::recover_how_many(&params, rlwe_q_prime_1, rlwe_q_prime_2, gamma_second_layer_mask, ct_bytes);
                response_mask.push(ct);
            }
            let mut response_body = Vec::new();
            for ct_bytes in response_object.packed_body_mod_switched.iter() {
                let ct = PolyMatrixRaw::recover_how_many(&params, rlwe_q_prime_1, rlwe_q_prime_2, gamma_second_layer_body, ct_bytes);
                response_body.push(ct);
            }            

            let mut inner_ct = PolyMatrixRaw::zero(&params, 2, 1);

            debug!("decrypting outer cts...");
            let outer_ct = response_mask
                .iter()
                .flat_map(|ct| {
                    decrypt_ct_reg_measured(y_client.client(), &params, &ct.ntt(), params.poly_len)
                        .as_slice()[..gamma_second_layer_mask]
                        .to_vec()
                })
                .collect::<Vec<_>>();
            assert_eq!(outer_ct.len(), mask_num_packed_rlwes * gamma_second_layer_mask);
            
            let outer_ct_t_u8 = u64s_to_contiguous_bytes(&outer_ct, pt_bits);

            let mut bit_offs = 0;
            for z in 0..params.poly_len {
                let val = read_bits(&outer_ct_t_u8, bit_offs, q_prime_bits);
                bit_offs += q_prime_bits;
                if !online_only {
                    assert!(val < rlwe_q_prime_2, "val: {}, rlwe_q_prime_2: {}", val, rlwe_q_prime_2);
                }
                inner_ct.data[z] = rescale(val, rlwe_q_prime_2, params.modulus);
            }

            let outer_ct : Vec<u64>;
            match second_level_packing_body {
                PackingType::NoPacking => {
                    outer_ct = response_body//[mask_num_packed_rlwes..mask_num_packed_rlwes + body_decomposed_num_coeffs]
                    .iter()
                    .map(|ct| {
                        decrypt_ct_reg_measured(y_client.client(), &params, &ct.ntt(), params.poly_len)
                            .as_slice()[0]
                    })
                    .collect::<Vec<_>>();
                    assert_eq!(outer_ct.len(), body_decomposed_num_coeffs);
                },
                PackingType::InspiRING => {
                    let gamma_second_layer_body = gammas[2];

                    let packed_online_num = (body_decomposed_num_coeffs as f64 / gamma_second_layer_body as f64).ceil() as usize;
                    outer_ct = response_body//[mask_num_packed_rlwes..mask_num_packed_rlwes + packed_online_num]
                    .iter()
                    .flat_map(|ct| {
                        decrypt_ct_reg_measured(y_client.client(), &params, &ct.ntt(), params.poly_len)
                            .as_slice()[0..gamma_second_layer_body].to_vec()
                    })
                    .collect::<Vec<_>>();
                    assert_eq!(outer_ct.len(), packed_online_num * gamma_second_layer_body);
                },
                PackingType::CDKS => {
                    panic!("Not implemented");
                },
            };
            let outer_ct_t_u8 = u64s_to_contiguous_bytes(&outer_ct, pt_bits);

            let mut bit_offs = 0;
            for z in 0..gamma_first_layer {
                let val = read_bits(&outer_ct_t_u8, bit_offs, q_prime_bits);
                bit_offs += q_prime_bits;
                if !online_only {
                    assert!(val < rlwe_q_prime_2, "val: {}, lwe_q_prime: {}", val, rlwe_q_prime_2);
                }
                inner_ct.data[params.poly_len + z] = rescale(val, rlwe_q_prime_2, params.modulus);
            }

            debug!("decrypting inner ct...");
            let decrypted_rlwe = decrypt_ct_reg_measured(
                y_client.client(),
                &params,
                &inner_ct.ntt(),
                params.poly_len,
            )
            .as_slice()[..gamma_first_layer]
            .to_vec();
            let final_result = decrypted_rlwe[target_col % gamma_first_layer];
            debug!("decrypted_rlwe: {}", final_result);

            measurement.online.client_decode_time_ms = start_decode.elapsed().as_millis() as usize;

            debug!("got {}, expected {}", final_result, corr_result);
            if online_only {
                println!("Warning! It's not correct");
            } else {
                assert_eq!(final_result, corr_result);
            }
        }

        measurement.online.upload_query = query_size;
        measurement.online.upload_keys = pub_params_size;
        measurement.online.upload_bytes = online_upload_bytes;
        measurement.online.download_bytes = online_download_bytes;
        measurement.online.total_bytes = online_upload_bytes + online_download_bytes;
        measurement.online.server_time_ms = online_server_time_ms as usize;
        // measurement.online.sqrt_n_bytes = sqrt_n_bytes;
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

    final_measurement.specs.protocol_type = ProtocolType::InsPIRe.to_str();
    final_measurement.specs.second_level_packing_mask = second_level_packing_mask.to_str();
    final_measurement.specs.second_level_packing_body = second_level_packing_body.to_str();
    final_measurement.specs.poly_len = params.poly_len;
    final_measurement.specs.modulus_bits = params.modulus_log2 as usize;
    final_measurement.specs.pt_modulus = params.pt_modulus as usize;

    final_measurement.specs.gamma_0 = gammas[0];
    final_measurement.specs.gamma_1 = gammas[1];
    if second_level_packing_body == PackingType::InspiRING {
        final_measurement.specs.gamma_2 = gammas[2];
    }
    
    final_measurement
}

pub fn mean(xs: &[usize]) -> f64 {
    xs.iter().map(|x| *x as f64).sum::<f64>() / xs.len() as f64
}

pub fn mean_f64(xs: &[f64]) -> f64 {
    xs.iter().map(|x| *x as f64).sum::<f64>() / xs.len() as f64
}

pub fn std_dev(xs: &[usize]) -> f64 {
    let mean = mean(xs);
    let mut variance = 0.;
    for x in xs {
        variance += (*x as f64 - mean).powi(2);
    }
    (variance / xs.len() as f64).sqrt()
}

pub fn std_dev_f64(xs: &[f64]) -> f64 {
    let mean = mean_f64(xs);
    let mut variance = 0.;
    for x in xs {
        variance += (*x as f64 - mean).powi(2);
    }
    (variance / xs.len() as f64).sqrt()
}

#[cfg(test)]
mod test {
    // use super::*;
    // use test_log::test;

    // #[test]
    // fn test_ypir_basic() {
    //     run_ypir_batched(1 << 30, 1, 1, false, false, 1);
    // }

    // #[test]
    // fn test_ypir_simplepir_basic() {
    //     run_ypir_batched(1 << 17, 65536 * 8, 1, true, false, 1);
    // }

    // #[test]
    // fn test_ypir_simplepir_rectangle() {
    //     run_ypir_batched(1 << 16, 16384 * 8, 1, true, false, 1);
    // }

    // #[test]
    // fn test_ypir_simplepir_rectangle_8gb() {
    //     run_ypir_batched(1 << 17, 65536 * 8, 1, true, false, 1);
    // }

    // #[test]
    // fn test_ypir_many_clients() {
    //     run_ypir_batched(1 << 30, 1, 2, false, false, 1);
    // }

    // #[test]
    // fn test_ypir_many_clients_and_trials() {
    //     run_ypir_batched(1 << 30, 1, 2, false, false, 5);
    // }

    // #[test]
    // #[ignore]
    // fn test_ypir_1gb() {
    //     run_ypir_batched(1 << 33, 1, 1, false, false, 5);
    // }

    // #[test]
    // #[ignore]
    // fn test_ypir_2gb() {
    //     run_ypir_batched(1 << 34, 1, 1, false, false, 5);
    // }

    // #[test]
    // #[ignore]
    // fn test_ypir_4gb() {
    //     run_ypir_batched(1 << 35, 1, 1, false, false, 5);
    // }

    // #[test]
    // #[ignore]
    // fn test_ypir_8gb() {
    //     run_ypir_batched(1 << 36, 1, 1, false, false, 5);
    // }

    // #[test]
    // #[ignore]
    // fn test_ypir_16gb() {
    //     run_ypir_batched(1 << 37, 1, 1, false, false, 5);
    // }

    // #[test]
    // #[ignore]
    // fn test_ypir_32gb() {
    //     run_ypir_batched(1 << 38, 1, 1, false, false, 5);
    // }

    // #[test]
    // #[ignore]
    // fn test_batched_4_ypir() {
    //     run_ypir_batched(1 << 30, 1, 4, false, false, 5);
    // }
}
