use serde::{Deserialize, Serialize};

use spiral_rs::poly::PolyMatrixNTT;

#[derive(Serialize, Deserialize, Debug, Default, Clone)]
#[serde(rename_all = "camelCase")]
pub struct Measurement {
    pub specs: Specs, 
    pub offline: Offline,
    pub online: Online,
}

// do serde rename
#[derive(Serialize, Deserialize, Debug, Default, Clone)]
#[serde(rename_all = "camelCase")]
pub struct Specs {
    pub label: String,
    pub input_database_size_mb: f64,
    pub input_num_items: usize,
    pub input_item_size_bits: usize,
    pub poly_len: usize,
    pub modulus_bits: usize,
    pub pt_modulus: usize,
    pub protocol_type: String,
    pub second_level_packing_mask: String,
    pub second_level_packing_body: String,
    pub gamma_0: usize,
    pub gamma_1: usize,
    pub gamma_2: usize,
    pub interpolate_degree: usize,
    pub resized_db_first_dim: usize,
    pub resized_db_second_dim: usize,
    pub resized_item_size_bits: usize,
    pub resized_database_size_mb: f64,
    pub online_only: bool,
}

// do serde rename
#[derive(Serialize, Deserialize, Debug, Default, Clone)]
#[serde(rename_all = "camelCase")]
pub struct Offline {
    pub upload_bytes: usize,
    pub download_bytes: usize,
    pub encode_time_ms: usize,
    pub server_time_ms: usize,
    pub client_time_ms: usize,
    pub simplepir_prep_time_ms: usize,
    pub simplepir_hint_bytes: usize,
    pub doublepir_hint_bytes: usize,
}

#[derive(Serialize, Deserialize, Debug, Default, Clone)]
#[serde(rename_all = "camelCase")]
pub struct Online {
    pub upload_bytes: usize,
    pub upload_keys: usize,
    pub upload_query: usize,
    pub download_bytes: usize,
    pub total_bytes: usize,
    pub simplepir_query_bytes: usize,
    pub doublepir_query_bytes: usize,
    pub simplepir_resp_bytes: usize,
    pub doublepir_resp_bytes: usize,
    pub server_time_ms: usize,
    pub client_query_gen_time_ms: usize,
    pub client_decode_time_ms: usize,

    pub first_pass_time_us: usize,
    pub second_pass_time_us: usize,

    pub first_pack_time_us: usize,
    pub rgsw_time_us: usize,

    pub packing_key_rotations_time_us: usize,
    pub packing_mask_online_time_us: usize,
    pub packing_body_online_time_us: usize,
    pub total_ring_packing_time_ms: usize,
    pub sqrt_n_bytes: usize,
    pub all_server_times_ms: Vec<usize>,
    pub std_dev_server_time_ms: f64,

    pub noise_width: f64,

}

pub fn get_vec_pm_size_bytes(v_p: &[PolyMatrixNTT]) -> usize {
    v_p.len()
        * v_p[0].rows
        * v_p[0].cols
        * v_p[0].params.poly_len
        * v_p[0].params.modulus_log2 as usize
        / 8
}

pub fn get_size_bytes(response: &[Vec<Vec<u8>>]) -> usize {
    let mut size_bytes = 0;
    for i in 0..response.len() {
        for j in 0..response[i].len() {
            size_bytes += response[i][j].len();
        }
    }
    size_bytes
}

pub fn get_size_bytes_one(response: &Vec<Vec<u8>>) -> usize {
    let mut size_bytes = 0;
        for j in 0..response.len() {
            size_bytes += response[j].len();
        }
    size_bytes
}
