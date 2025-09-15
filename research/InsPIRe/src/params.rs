use log::debug;
use serde_json::Value;

use spiral_rs::{arith::*, params::*};

use super::lwe::LWEParams;

static DEFAULT_MODULI: [u64; 2] = [268369921u64, 249561089u64];
pub const DEF_MOD_STR: &str = "[\"268369921\", \"249561089\"]";

pub fn ext_params_from_json(json_str: &str) -> Params {
    let v: Value = serde_json::from_str(json_str).unwrap();

    let n = v["n"].as_u64().unwrap() as usize;
    let db_dim_1 = v["nu_1"].as_u64().unwrap() as usize;
    let db_dim_2 = v["nu_2"].as_u64().unwrap() as usize;
    let instances = v["instances"].as_u64().unwrap_or(1) as usize;
    let p: u64 = v["p"].as_u64().unwrap();
    let q2_bits = u64::max(v["q2_bits"].as_u64().unwrap(), MIN_Q2_BITS);
    let t_gsw = v["t_gsw"].as_u64().unwrap() as usize;
    let t_conv = v["t_conv"].as_u64().unwrap() as usize;
    let t_exp_left = v["t_exp_left"].as_u64().unwrap() as usize;
    let t_exp_right = v["t_exp_right"].as_u64().unwrap() as usize;
    let do_expansion = v.get("direct_upload").is_none();

    let mut db_item_size = v["db_item_size"].as_u64().unwrap_or(0) as usize;
    if db_item_size == 0 {
        db_item_size = instances * n * n;
        db_item_size = db_item_size * 2048 * log2_ceil(p) as usize / 8;
    }

    let version = v["version"].as_u64().unwrap_or(0) as usize;

    let poly_len = v["poly_len"].as_u64().unwrap_or(2048) as usize;
    let moduli = v["moduli"]
        .as_array()
        .map(|x| {
            x.as_slice()
                .iter()
                .map(|y| {
                    y.as_u64()
                        .unwrap_or_else(|| y.as_str().unwrap().parse().unwrap())
                })
                .collect::<Vec<_>>()
        })
        .unwrap_or(DEFAULT_MODULI.to_vec());
    let noise_width = v["noise_width"].as_f64().unwrap_or(6.4);

    Params::init(
        poly_len,
        &moduli,
        noise_width,
        n,
        p,
        q2_bits,
        t_conv,
        t_exp_left,
        t_exp_right,
        t_gsw,
        do_expansion,
        db_dim_1,
        db_dim_2,
        instances,
        db_item_size,
        version,
    )
}

pub fn internal_params_for(
    poly_len : usize,
    nu_1: usize,
    nu_2: usize,
    p: u64,
    q2_bits: usize,
    t_exp_left: usize,
    moduli: &str,
) -> Params {
    ext_params_from_json(&format!(
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
        poly_len, nu_1, nu_2, p, q2_bits, t_exp_left, moduli
    ))
}

pub fn round_up_to_multiple_of(x : usize, multiplicand: usize) -> usize {
    (x as f64 / multiplicand as f64).ceil() as usize * multiplicand
}

pub fn get_variance(dim: f64, p: f64, sigma_x: f64, ell_ks: f64, z: f64, poly_len: f64, gamma: f64) -> f64 {
    (dim * p.powi(2) * sigma_x.powi(2) + ell_ks * gamma * poly_len * z.powi(2) * sigma_x.powi(2) / 4.0).log2() / 2.0
}

pub fn params_for_scenario_medium_payload(input_num_items: usize, input_item_size_bits: usize, gammas: Vec<usize>, performance_factor: usize) -> (Params, (usize, usize, usize)) {
    let gamma_0 = gammas[0];
    let gamma_1 = gammas[1];
    let gamma_2 = gammas[2];

    let poly_len_log2 = 11;
    let poly_len = 1 << poly_len_log2; 

    let log_p = 16;
    assert!(log_p <= 16);

    let p = 1 << log_p;
    let q2_bits = 28;
    let t_exp_left = 3;

    let working_num_large_items = (input_num_items as f64 * input_item_size_bits as f64 / (log_p * gamma_0) as f64).ceil();

    let num_tiles = working_num_large_items as f64 / (poly_len * poly_len) as f64;
    let num_tiles_log2 = num_tiles.ceil().log2().ceil() as usize;

    let log_factor = performance_factor.trailing_zeros() as usize;
    let (nu_1, nu_2) = if num_tiles_log2 % 2 == 0 {
        ((num_tiles_log2 / 2) + log_factor, (num_tiles_log2 / 2) - log_factor)
    } else {
        (((num_tiles_log2 + 1) / 2) + log_factor, ((num_tiles_log2 - 1) / 2) - log_factor)
    };

    let db_rows = nu_1;
    let db_cols = 1 << (nu_2 + (gamma_0.trailing_zeros() as usize));

    debug!("db_rows: {}, instances: {}", db_rows, db_cols);

    let size_over_t = 1 << (nu_1 + poly_len_log2);
    let t = 1 << (nu_2 + poly_len_log2); 
    let sigma_x = 6.4 as f64;
    let z = (1 << 19) as f64;

    let term_0_variance = get_variance(size_over_t as f64, p as f64, sigma_x, t_exp_left as f64, z as f64, poly_len as f64, gamma_0 as f64);
    let term_1_variance = get_variance(t as f64, p as f64, sigma_x, t_exp_left as f64, z as f64, poly_len as f64, gamma_1 as f64);
    let term_2_variance = get_variance(t as f64, p as f64, sigma_x, t_exp_left as f64, z as f64, poly_len as f64, gamma_2 as f64);

    // // print the three terms
    // println!("term_0_variance: {}", term_0_variance);
    // println!("term_1_variance: {}", term_1_variance);
    // println!("term_2_variance: {}", term_2_variance);

    let max_variance = term_0_variance.max(term_1_variance).max(term_2_variance);
    // println!("max_variance: {}", max_variance);

    let modulus_log2 = (2.*2.* (p as f64)).log2() + (2.0 * 41.0 * (2 as f64).ln()).sqrt().log2() + max_variance;
    // println!("best modulus_log2: {}", modulus_log2);

    pub const CUSTOM_MOD_STR: &str = "[\"67043329\", \"132120577\"]";
    let modulus: f64 = 67043329. * 132120577.;
    assert!(modulus.log2() >= modulus_log2);
    let mut params: Params = internal_params_for(poly_len, nu_1, 0, p, q2_bits, t_exp_left, CUSTOM_MOD_STR);
    params.instances = db_cols;
    (params, (1 << (nu_1 + poly_len_log2), 1 << (nu_2 + poly_len_log2), log_p * gamma_0))
}

pub fn params_for_scenario(input_num_items: usize, input_item_size_bits: usize) -> (Params, (usize, usize, usize)) {
    let poly_len_log2 = 11;
    let poly_len = 1 << poly_len_log2;

    // The plaintext moodulus of the first layer is 8 bits
    assert!(input_item_size_bits <= 8);

    let total_db_bytes = input_num_items * input_item_size_bits / 8;
    let lwe_pt_word_bytes = 1;
    let num_items = total_db_bytes / lwe_pt_word_bytes;
    let num_tiles = num_items as f64 / (2048. * 2048.); // shouldn't this be d1 and d2?
    let num_tiles_usize = num_tiles.ceil() as usize;
    let num_tiles_log2 = (num_tiles_usize as f64).log2().ceil() as usize;

    let (nu_1, nu_2) = if num_tiles_log2 % 2 == 0 {
        (num_tiles_log2 / 2, num_tiles_log2 / 2)
    } else {
        ((num_tiles_log2 + 1) / 2, (num_tiles_log2 - 1) / 2)
    };

    debug!("chose nu_1: {}, nu_2: {}", nu_1, nu_2);


    let logp = 15; // This is the plaintext modulus of the second layer ONLY
    let q2_bits = 28;
    let t_exp_left = 3;
    let p = 1 << logp;

    let params = internal_params_for(poly_len, nu_1, nu_2, p, q2_bits, t_exp_left, DEF_MOD_STR);
    (params, (1 << (nu_1 + poly_len_log2), 1 << (nu_2 + poly_len_log2), 8))

}

pub fn params_for_scenario_simplepir(input_num_items: usize, input_item_size_bits: usize) -> (Params, (usize, usize, usize)) {

    let poly_len_log2 = 11;
    let poly_len = 1 << poly_len_log2; 
    let log_p = 14;

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

    // let CUSTOM_DEF_MOD_STR: &str = "[\"61441\", \"12289\"]";    

    let mut params = internal_params_for(poly_len, nu_1, 0, p, q2_bits, t_exp_left, DEF_MOD_STR);
    params.instances = db_cols_poly;
    (params, (new_num_items, 1, new_item_size_num_pts*log_p))
}

pub fn params_for_scenario_simplepir_small(input_num_items: usize, input_item_size_bits: usize) -> (Params, (usize, usize, usize)) {

    let poly_len_log2 = 10;
    let poly_len = 1 << poly_len_log2; 
    let log_p = 6;

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
    let t_exp_left = 8;

    // let CUSTOM_DEF_MOD_STR: &str = "[\"16515073\", \"33292289\"]";    
    pub const CUSTOM_DEF_MOD_STR: &str = "[\"268369921\"]";

    let mut params = internal_params_for(poly_len, nu_1, 0, p, q2_bits, t_exp_left, CUSTOM_DEF_MOD_STR);
    params.instances = db_cols_poly;
    (params, (new_num_items, 1, new_item_size_num_pts*log_p))
}

pub trait GetQPrime {
    /// The smaller reduced modulus, used on the second row of the encoding
    fn get_q_prime_1(&self) -> u64;

    /// The larger reduced modulus, used on the first row of the encoding
    fn get_q_prime_2(&self) -> u64;
}

impl GetQPrime for Params {
    fn get_q_prime_1(&self) -> u64 {
        1 << 20
    }

    fn get_q_prime_2(&self) -> u64 {
        if self.q2_bits == self.modulus_log2 {
            self.modulus
        } else {
            Q2_VALUES[self.q2_bits as usize]
        }
    }
}

impl GetQPrime for LWEParams {
    fn get_q_prime_1(&self) -> u64 {
        u64::MAX // unsupported
    }

    fn get_q_prime_2(&self) -> u64 {
        if self.q2_bits == (self.modulus as f64).log2().ceil() as usize {
            self.modulus
        } else {
            Q2_VALUES[self.q2_bits as usize]
        }
    }
}