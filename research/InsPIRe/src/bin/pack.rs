use core::panic;
use serde::{Deserialize, Serialize};
use spiral_rs::{arith::*, poly::*};
use std::time::Instant;
use log::debug;

use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use spiral_rs::{client::Client, number_theory::invert_uint_mod, params::*};

use inspire::{packing::*, params::*, server::*, scheme::*};

fn file_name(
    packing_type: PackingType,
    pre_rotate: bool,
    gamma: usize,
    total_num_to_pack: usize,
) -> String {
    match packing_type {
        PackingType::CDKS => format!("{}-{}.json", packing_type.to_str(), total_num_to_pack),
        PackingType::InspiRING => format!(
            "{}-{}-gamma={}-prerotate={}.json",
            packing_type.to_str(),
            total_num_to_pack,
            gamma,
            pre_rotate
        ),
        PackingType::NoPacking => panic!(),
    }
}

fn mean_u128(vec: Vec<u128>) -> u128 {
    if vec.is_empty() {
        return 0;
    }
    let sum: u128 = vec.iter().sum();
    sum / vec.len() as u128
}

fn std_dev_u128(vec: Vec<u128>) -> f64 {
    if vec.is_empty() {
        return 0.0;
    }
    let mean = mean_u128(vec.clone());
    let mut sum_of_squares = 0;
    for x in &vec {
        sum_of_squares += (x - mean) * (x - mean);
    }
    let sum_of_squares = sum_of_squares as f64;
    let len = vec.len() as f64;
    (sum_of_squares / len).sqrt()
}

pub fn params_for_scenario_packing(input_num_items: usize, input_item_size_bits: usize) -> (Params, (usize, usize, usize)) {

    let poly_len_log2 = 11;
    let poly_len = 1 << poly_len_log2; 
    let log_p = 15;

    let input_item_num_pts = (input_item_size_bits as f64 / log_p as f64).ceil() as usize; // -> number of pt modului required to represent one item

    let total_db_num_pts = (((input_num_items * input_item_size_bits) as f64) / log_p as f64).ceil() as usize; // -> total number of pt moduli required to represent the database
    let sqrt_num_db_pts = (total_db_num_pts as f64).sqrt().floor() as usize; // -> sqrt of the previous line

    let new_item_size_num_pts = round_up_to_multiple_of(sqrt_num_db_pts, input_item_num_pts); // -> rounded up 
    let new_num_items = (total_db_num_pts as f64 / new_item_size_num_pts as f64).ceil() as usize;

    let db_rows = new_num_items;
    
    // number of required rlwe ciphertexts (in first layer)
    let db_cols_poly = (new_item_size_num_pts as f64 / (poly_len as f64)).ceil() as usize;

    debug!("db_rows: {}, db_cols_poly: {}", db_rows, db_cols_poly);

    let nu_1 = db_rows.next_power_of_two().trailing_zeros() as usize - poly_len_log2;
    debug!("chose nu_1: {}", nu_1);

    let p = 1 << log_p;
    let q2_bits = 28; // modulus after packing
    let t_exp_left = 3;

    let mut params = internal_params_for(poly_len, nu_1, 0, p, q2_bits, t_exp_left, DEF_MOD_STR);
    params.instances = db_cols_poly;
    (params, (new_num_items, 1, new_item_size_num_pts*log_p))
}

pub fn params_for_scenario_packing_small(input_num_items: usize, input_item_size_bits: usize) -> (Params, (usize, usize, usize)) {

    let poly_len_log2 = 10;
    let poly_len = 1 << poly_len_log2; 
    let log_p = 8;

    let input_item_num_pts = (input_item_size_bits as f64 / log_p as f64).ceil() as usize; // -> number of pt modului required to represent one item

    let total_db_num_pts = (((input_num_items * input_item_size_bits) as f64) / log_p as f64).ceil() as usize; // -> total number of pt moduli required to represent the database
    let sqrt_num_db_pts = (total_db_num_pts as f64).sqrt().floor() as usize; // -> sqrt of the previous line

    let new_item_size_num_pts = round_up_to_multiple_of(sqrt_num_db_pts, input_item_num_pts); // -> rounded up 
    let new_num_items = (total_db_num_pts as f64 / new_item_size_num_pts as f64).ceil() as usize;

    let db_rows = new_num_items;
    
    // number of required rlwe ciphertexts (in first layer)
    let db_cols_poly = (new_item_size_num_pts as f64 / (poly_len as f64)).ceil() as usize;

    debug!("db_rows: {}, db_cols_poly: {}", db_rows, db_cols_poly);

    let nu_1 = db_rows.next_power_of_two().trailing_zeros() as usize - poly_len_log2;
    debug!("chose nu_1: {}", nu_1);

    let p = 1 << log_p;
    let q2_bits = 30; // modulus after packing
    let t_exp_left = 10;

    // pub const CUSTOM_DEF_MOD_STR: &str = "[\"4293918721\"]";
    // pub const CUSTOM_DEF_MOD_STR: &str = "[\"268369921\"]";
    pub const CUSTOM_DEF_MOD_STR: &str = "[\"1073479681\"]";
    // let mut params = internal_params_for(poly_len, nu_1, 0, p, q2_bits, t_exp_left, CUSTOM_DEF_MOD_STR);

    let mut params = ext_params_from_json(&format!(
        r#"
        {{
            "poly_len": {},
            "n": 1,
            "nu_1": {},
            "nu_2": {},
            "p": {},
            "q2_bits": {},
            "t_gsw": 3,
            "t_conv": 4,
            "t_exp_left": {},
            "t_exp_right": 2,
            "instances": 1,
            "db_item_size": 0,
            "moduli": {},
            "noise_width": 16.042421
        }}
        "#,
        poly_len, nu_1, 0, p, q2_bits, t_exp_left, CUSTOM_DEF_MOD_STR
    ));

    params.instances = db_cols_poly;
    (params, (new_num_items, 1, new_item_size_num_pts*log_p))
}

fn benchmark_packing(
    packing_type: PackingType,
    pre_rotate: bool,
    gamma: usize,
    total_num_to_pack: usize,
    num_trials: usize,
    version: &str,
) {
    let filename = file_name(packing_type, pre_rotate, gamma, total_num_to_pack);

    let (params, _) = params_for_scenario_packing(1 << 30, 1);
    let mut measurement = PackingMeasurement::default();

    measurement.packing_type = packing_type.to_str();
    measurement.num_trials = num_trials;
    match packing_type {
        PackingType::CDKS => {
            measurement.total_num_to_pack = total_num_to_pack;
            measurement.gamma = params.poly_len;
            assert!(total_num_to_pack % params.poly_len == 0);
        }
        PackingType::InspiRING => {
            measurement.total_num_to_pack = total_num_to_pack;
            measurement.gamma = gamma;
            measurement.pre_rotate = pre_rotate;
        }
        PackingType::NoPacking => panic!("Shouldn't be here"),
    };

    let gamma = measurement.gamma;

    let num_sets = (total_num_to_pack as f64 / gamma as f64).ceil() as usize;

    let packing_params = PackParams::new(&params, gamma);

    let pack_seed = [1u8; 32];
    let cts_seed = [2u8; 32];
    let mut ct_pub_rng = ChaCha20Rng::from_seed(cts_seed);

    let mut client = Client::init(&params);
    client.generate_secret_keys();
    let sk_reg = client.get_sk_reg();

    let y_constants = generate_y_constants(&params);

    let mut packing_keys = match packing_type {
        PackingType::CDKS => PackingKeys::init_cdks(&params, sk_reg, pack_seed),
        PackingType::InspiRING => {
            if gamma <= params.poly_len / 2 {
                PackingKeys::init(&packing_params, sk_reg, W_SEED)
            } else {
                PackingKeys::init_full(&packing_params, sk_reg, W_SEED, V_SEED)
            }
        },
        PackingType::NoPacking => panic!()
    };
    measurement.keys_size_bytes = packing_keys.get_size_bytes();

    let mut gold = Vec::new();
    for i in 0..total_num_to_pack {
        gold.push(i as u64 % params.pt_modulus);
    }

    // generate poly_len ciphertexts
    let mut v_ct = Vec::new();
    let mut b_values = Vec::new();
    let scale_k = params.modulus / params.pt_modulus;
    for i in 0..total_num_to_pack {
        let mut pt = PolyMatrixRaw::zero(&params, 1, 1);
        let val = gold[i] as u64 % params.pt_modulus;
        let mod_inv = match packing_type {
            PackingType::CDKS => invert_uint_mod(gamma as u64, params.modulus).unwrap(),
            PackingType::InspiRING => 1,
            PackingType::NoPacking => panic!(),
        };
        let val_to_enc = multiply_uint_mod(val * scale_k, mod_inv, params.modulus);
        pt.data[0] = val_to_enc;
        let ct =
            client.encrypt_matrix_reg(&pt.ntt(), &mut ChaCha20Rng::from_entropy(), &mut ct_pub_rng);
        let mut ct_raw = ct.raw();

        // get the b value
        b_values.push(ct_raw.get_poly(1, 0)[0]);

        // zero out all of the second poly
        ct_raw.get_poly_mut(1, 0).fill(0);
        v_ct.push(ct_raw.ntt());
    }

    let mut a_ct_tilde = Vec::new();
    for i in 0..total_num_to_pack {
        a_ct_tilde.push(v_ct[i].submatrix(0, 0, 1, 1));
    }

    let mut precomp_cdks_list = Vec::new();
    let mut precomp_inspiring_list = Vec::new();

    let w_mask = PolyMatrixNTT::random_rng(
        &params,
        1,
        params.t_exp_left,
        &mut ChaCha20Rng::from_seed(W_SEED),
    );
    let v_mask = PolyMatrixNTT::random_rng(
        &params,
        1,
        params.t_exp_left,
        &mut ChaCha20Rng::from_seed(V_SEED),
    );

    let now = Instant::now();
    for i in 0..num_sets {
        match packing_type {
            PackingType::CDKS => {
                let v_ct_slice = v_ct[i * gamma..(i + 1) * gamma].to_vec();
                precomp_cdks_list.push(precompute_pack(
                    &params,
                    params.poly_len_log2,
                    &v_ct_slice,
                    &packing_keys.fake_pack_pub_params,
                    &y_constants,
                ));
            }
            PackingType::InspiRING => {
                let a_ct_tilde_select = a_ct_tilde[i * gamma..(i + 1) * gamma].to_vec();
                if gamma <= params.poly_len / 2 {
                    let w_all = generate_rotations(&packing_params, &w_mask);
                    if pre_rotate {
                        precomp_inspiring_list.push(packing_with_preprocessing_offline(
                            &packing_params,
                            &w_all,
                            &a_ct_tilde_select,
                        ));
                    } else {
                        precomp_inspiring_list.push(
                            packing_with_preprocessing_offline_without_rotations(
                                &packing_params,
                                &w_all,
                                &a_ct_tilde_select,
                            ),
                        );
                    }
                } else {
                    let (w_all, w_bar_all) = generate_rotations_double(&packing_params, &w_mask);
                    if pre_rotate {
                        precomp_inspiring_list.push(full_packing_with_preprocessing_offline(
                            &packing_params,
                            &w_all,
                            &w_bar_all,
                            &v_mask,
                            &a_ct_tilde_select,
                        ));
                    } else {
                        precomp_inspiring_list.push(
                            full_packing_with_preprocessing_offline_without_rotations(
                                &packing_params,
                                &w_all,
                                &w_bar_all,
                                &v_mask,
                                &a_ct_tilde_select,
                            ),
                        );
                    }
                }
            }
            PackingType::NoPacking => panic!(),
        }
    }

    let offline_time = now.elapsed().as_micros();
    // println!("Precomputing for packing took {} us", offline_time);
    measurement.time_offline_us = offline_time;

    let mut time_key_expansion_us_list = Vec::new();
    let mut packing_time_us_list = Vec::new();
    let mut online_time_us_list = Vec::new();
    let mut noise_bits = Vec::new();

    for _ in 0..num_trials {
        let online_time = Instant::now();

        let mut packing_time = Instant::now();
        let packed = match packing_type {
            PackingType::CDKS => {
                let pack_pub_params_row_1s = &packing_keys.pack_pub_params_row_1s;
                pack_many_lwes(
                    &params,
                    // &prepacked_lwe,
                    &precomp_cdks_list,
                    &b_values,
                    num_sets,
                    &pack_pub_params_row_1s,
                    &y_constants,
                )
            }
            PackingType::InspiRING => {
                if pre_rotate {
                    let time_key_expansion = Instant::now();
                    packing_keys.expand();
                    time_key_expansion_us_list.push(time_key_expansion.elapsed().as_micros());
                    packing_time = Instant::now();
                    pack_many_lwes_inspir(
                        &packing_params,
                        &precomp_inspiring_list,
                        &b_values,
                        &packing_keys,
                        gamma,
                    )
                } else {
                    packing_time = Instant::now();
                    pack_many_lwes_inspir_without_rotations(
                        &packing_params,
                        &precomp_inspiring_list,
                        &b_values,
                        &packing_keys,
                        gamma,
                    )
                }
            }
            PackingType::NoPacking => panic!(),
        };
        packing_time_us_list.push(packing_time.elapsed().as_micros());

        // let q1 = params.get_q_prime_1();
        // let q2 = params.get_q_prime_2();

        // let mut packed_mod_switched = Vec::new();
        // for ct in packed.iter() {
        //     let ct_switched = ct.switch(q1, q2);
        //     packed_mod_switched.push(ct_switched);
        // }

        online_time_us_list.push(online_time.elapsed().as_micros());
        // println!("Packing took {} us", online_time);

        // measurement.output_size_bytes = get_size_bytes(&vec![packed_mod_switched]);

        measurement.output_size_bytes = ((packed.len() * params.poly_len * 2 * params.modulus_log2 as usize) + 7) / 8;

        // decrypt + decode
        let mut rescaled = Vec::new();
        for i in 0..num_sets {
            let dec = client.decrypt_matrix_reg(&packed[i].ntt());
            let dec_raw = dec.raw();

            // noise_bits.push(measure_noise_width_bits(&params, &dec, gamma)as f64);

            // rescale
            // let mut errors = Vec::new();
            let mut max_noise = 0;
            for i in 0..gamma {
                let mu = dec_raw.data[i];
                let mu_i64 = mu as i64;
                let q_i64 = params.modulus as i64;
                // println!("mod bits: {}", (params.modulus as f64).log2());
                let q_half_i64 = (params.modulus / 2) as i64;
                let message = rescale(mu, params.modulus, params.pt_modulus);
                rescaled.push(message);

                let scaled_message = rescale(message, params.pt_modulus, params.modulus) as i64;                

                let mut diff = if mu_i64 > scaled_message { 
                    mu_i64 - scaled_message
                } else {
                    q_i64 + mu_i64 - scaled_message
                };

                if diff >= q_half_i64 {
                    diff = q_i64 - diff
                }

                if diff > max_noise {
                    max_noise = diff;
                }

            }
            noise_bits.push((max_noise as f64).log2());
        }
        // println!("modulus bitlength: {}", params.modulus_log2);
        // let max_error = errors.iter().max().unwrap();
        // println!("error bitlength: {}", (*max_error as f64).log2());


        // println!("rescaled: {:?}", &rescaled.as_slice()[..50]);
        assert_eq!(rescaled.as_slice()[..total_num_to_pack], gold.as_slice()[..total_num_to_pack]);
    }

    measurement.time_key_expansion_us = mean_u128(time_key_expansion_us_list.clone());
    measurement.packing_time_us = mean_u128(packing_time_us_list.clone());
    measurement.time_online_us = mean_u128(online_time_us_list.clone());
    measurement.time_key_expansion_std_dev = std_dev_u128(time_key_expansion_us_list.clone()) as f64;
    measurement.packing_time_std_dev = std_dev_u128(packing_time_us_list.clone()) as f64;
    measurement.time_online_std_dev = std_dev_u128(online_time_us_list.clone()) as f64;

    measurement.noise_bits_mean = mean_f64(noise_bits.as_slice());
    measurement.noise_bits_std = std_dev_f64(noise_bits.as_slice());
    measurement.noise_bits_max = *(noise_bits.iter().max_by(|a, b| a.total_cmp(b)).unwrap());

    measurement.d = params.poly_len;
    measurement.log_p = params.pt_modulus.trailing_zeros() as usize;
    measurement.log_q = params.modulus_log2 as usize;


    // (measurement, filename)
    let dirname = format!("evaluation/results/packing/{}", version);
    std::fs::create_dir_all(&dirname).unwrap();
    let mut file = std::fs::File::create(format!("{}/{}", dirname, filename)).unwrap();
    serde_json::to_writer_pretty(&mut file, &measurement).unwrap();
}

#[derive(Serialize, Deserialize, Debug, Default, Clone)]
#[serde(rename_all = "camelCase")]
pub struct PackingMeasurement {
    pub log_p: usize,
    pub log_q: usize,
    pub d: usize,

    pub packing_type: String,
    pub pre_rotate: bool,
    pub total_num_to_pack: usize,
    pub gamma: usize,
    pub num_trials: usize,
    
    pub time_offline_us: u128,
    pub time_key_expansion_us: u128,
    pub packing_time_us: u128,
    pub time_online_us: u128,
    pub time_key_expansion_std_dev: f64,
    pub packing_time_std_dev: f64,
    pub time_online_std_dev: f64,

    pub keys_size_bytes: usize,
    pub output_size_bytes: usize,
    pub noise_bits_mean: f64,
    pub noise_bits_std: f64,
    pub noise_bits_max: f64,

}


fn main() {
    let num_trials = 3;
    let version = "";

    for i in 1..2 {
        let total_num_to_pack = 4096 * i;
        println!("total_num_to_pack={}", total_num_to_pack);

        benchmark_packing(PackingType::CDKS, true, 0, total_num_to_pack, num_trials, version);
        for gamma in vec![2048, 1024, 512, ] {
            println!("gamma={}", gamma);
            benchmark_packing(
                PackingType::InspiRING,
                true,
                gamma,
                total_num_to_pack,
                num_trials,
                version,
            );
            benchmark_packing(
                PackingType::InspiRING,
                false,
                gamma,
                total_num_to_pack,
                num_trials,
                version,
            );
        }
    }

}
