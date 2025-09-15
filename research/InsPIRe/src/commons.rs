use super::packing::PackingType;
use spiral_rs::aligned_memory::AlignedMemory64;
use spiral_rs::{
    params::*,
};
use spiral_rs::poly::{PolyMatrix, PolyMatrixNTT};
use super::{packing::*, params::*};
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use pad::PadStr;

pub const RGSW_SEEDS: [[u8; 32]; 6] = [
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,],
    [7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,],
];

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
    packing_params: &'a PackParams,
    all_u8: Vec<u8>,
) -> (PackingKeys<'a>, AlignedMemory64, PolyMatrixNTT<'a>) {

    let all_u64: Vec<_> = all_u8.as_slice()
     .chunks_exact(8)
     .map(|chunk| u64::from_ne_bytes(chunk.try_into().unwrap()))
     .collect();
    

    let mut cnt = 0; 
    let db_rows = 1 << (params.db_dim_1 + params.poly_len_log2);
    // let packing_params = PackParams::new(&params, params.poly_len);

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
        packing_params: Some(packing_params.clone()),
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

fn max_interpolate_degree(modulus: f64, d0: f64, p: f64, t_exp: f64, t_rgsw: f64, poly_len: f64) -> usize {
    let sigma_x = 6.4 as f64;
    let modulus_len = modulus.log2();
    let z1 = (2.0 as f64).powi((modulus_len / t_exp).ceil() as i32);
    let z2 = (2.0 as f64).powi((modulus_len / t_rgsw).ceil() as i32);
    let term1_variance = d0 * p.powi(2) * sigma_x.powi(2);
    let term2_variance = t_exp * poly_len.powi(2) * z1.powi(2) * sigma_x.powi(2) / 4.0;
    let term3_variance = t_rgsw * poly_len * z2.powi(2) * sigma_x.powi(2) / 2.0;

    let total_log2_std_before_poly_eval = (term1_variance + term2_variance + term3_variance).sqrt().log2();

    let log2_std_upper_bound = (modulus / (2.*2.*p)).log2() - (2.0 * 41.0 * (2 as f64).ln()).sqrt().log2();
    let max_log_interpolate_degree = 2. * (log2_std_upper_bound - total_log2_std_before_poly_eval);
    assert!(max_log_interpolate_degree >= 0.);
    let max_interpolate_degree = (2.0 as f64).powf(max_log_interpolate_degree).floor() as usize;
    std::cmp::min(max_interpolate_degree as usize, poly_len as usize)

}

pub fn params_rgswpir_given_input_size_and_dim0<'a>(input_num_items: usize, input_item_size_bits: usize, dim0: usize) -> (Params, usize, (usize, usize, usize)) {

    let poly_len_log2 = 11;
    let poly_len = 1 << poly_len_log2; 
    let p = 65535;
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


pub fn read_file_into_matrix(path: &Path, num_items: usize, item_size_bits: usize) -> io::Result<Vec<Vec<u16>>> {
    let file = File::open(path)?;
    // Use a BufReader for efficient, buffered I/O.
    let reader = BufReader::new(file);

    // Pre-allocate the outer vector since we know the number of lines.
    let mut matrix: Vec<Vec<u16>> = Vec::with_capacity(num_items);

    for line_result in reader.lines() {
        let line = line_result?; // Propagate potential I/O errors.

        let line = line.pad_to_width(item_size_bits / 8);

        // Process the line as a byte slice for performance with ASCII.
        let row: Vec<u16> = line.as_bytes()
            // Group the bytes into non-overlapping chunks of 2.
            .chunks_exact(2)
            // Map each pair of bytes [c1, c2] to a single u16.
            .map(|chunk| {
                // Pack two 8-bit bytes into one 16-bit integer.
                // The first byte becomes the "high" byte, the second is the "low" byte.
                (chunk[0] as u16) << 8 | (chunk[1] as u16)
            })
            .collect();
        
        matrix.push(row);
    }

    while matrix.len() < num_items {
        matrix.push(vec![0u16; item_size_bits / 16]);
    }

    Ok(matrix)
}