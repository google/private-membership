use std::arch::x86_64::*;

use spiral_rs::{arith::*, params::*};

use super::server::ToM512;

pub fn fast_batched_dot_product_avx512<const K: usize, T: Copy>(
    params: &Params,
    c: &mut [u64],
    a: &[u64],
    a_elems: usize,
    b_t: &[T], // transposed
    b_rows: usize,
    b_cols: usize,
) where
    *const T: ToM512,
{
    assert_eq!(a_elems, b_rows);

    // debug!("Multiplying {}x{} by {}x{}", K, a_elems, b_rows, b_cols);

    let mut simd_width = 8;
    while simd_width > a_elems {
        simd_width = simd_width / 2;   
    }

    let chunk_size = (8192 / K.next_power_of_two()).min(a_elems / simd_width);
    let num_chunks = (a_elems / simd_width) / chunk_size;
    // debug!("k_chunk_size: {}, k_num_chunks: {}", chunk_size, num_chunks);

    let j_chunk_size = 1;
    let j_num_chunks = b_cols / j_chunk_size;
    // debug!(
    //     "j_chunk_size: {}, j_num_chunks: {}",
    //     j_chunk_size, j_num_chunks
    // );

    // let mut result = AlignedMemory64::new(b_cols);
    // let res_mut_slc = result.as_mut_slice();
    let res_mut_slc = c;

    unsafe {
        let mut a_slcs: [&[u64]; K] = [&[]; K];
        for (slc_mut, chunk) in a_slcs.iter_mut().zip(a.chunks_exact(a.len() / K)) {
            *slc_mut = chunk;
        }

        // let a_ptr = a.as_ptr();
        let b_ptr = b_t.as_ptr();

        for k_outer in 0..num_chunks {
            for j_outer in 0..j_num_chunks {
                for j_inner in 0..j_chunk_size {
                    let j = j_outer * j_chunk_size + j_inner;

                    let mut total_sum_lo = [_mm512_setzero_si512(); K];
                    let mut total_sum_hi = [_mm512_setzero_si512(); K];
                    let mut tmp = [_mm512_setzero_si512(); K];

                    for k_inner in 0..chunk_size {
                        let k = simd_width * (k_outer * chunk_size + k_inner);

                        let a_idx = k;
                        let b_idx = j * b_rows + k;
                        let b_val_simd = b_ptr.add(b_idx).to_m512();

                        for batch in 0..K {
                            tmp[batch] =
                                _mm512_load_si512(a_slcs[batch].as_ptr().add(a_idx) as *const _);
                        }


                        for batch in 0..K {
                            let a_val_lo = tmp[batch];
                            let a_val_hi = _mm512_srli_epi64(tmp[batch], 32);

                            total_sum_lo[batch] = _mm512_add_epi64(
                                total_sum_lo[batch],
                                _mm512_mul_epu32(a_val_lo, b_val_simd),
                            );
                            total_sum_hi[batch] = _mm512_add_epi64(
                                total_sum_hi[batch],
                                _mm512_mul_epu32(a_val_hi, b_val_simd),
                            );
                        }
                    }

                    let res_mut_slcs = res_mut_slc.chunks_exact_mut(res_mut_slc.len() / K);
                    for (batch, res_mut_slc) in (0..K).zip(res_mut_slcs) {
                        let mut values_lo = [0u64; 8];
                        let mut values_hi = [0u64; 8];
                        _mm512_store_si512(
                            (&mut values_lo).as_mut_ptr() as *mut _,
                            total_sum_lo[batch],
                        );
                        _mm512_store_si512(
                            (&mut values_hi).as_mut_ptr() as *mut _,
                            total_sum_hi[batch],
                        );

                        let res_lo = values_lo[0]
                            + values_lo[1]
                            + values_lo[2]
                            + values_lo[3]
                            + values_lo[4]
                            + values_lo[5]
                            + values_lo[6]
                            + values_lo[7];

                        let res_hi = values_hi[0]
                            + values_hi[1]
                            + values_hi[2]
                            + values_hi[3]
                            + values_hi[4]
                            + values_hi[5]
                            + values_hi[6]
                            + values_hi[7];

                        // res_mut_slc[j] = res_lo + res_hi;

                        let (lo, hi) = (
                            barrett_coeff_u64(params, res_lo as u64, 0),
                            barrett_coeff_u64(params, res_hi as u64, 1),
                        );

                        // res_mut_slc[j] = lo | (hi << 32);

                        let res = params.crt_compose_2(lo, hi);
                        res_mut_slc[j] = barrett_u64(params, res_mut_slc[j] + res);
                    }
                }
            }
        }
    }
}

#[repr(C)]
#[repr(align(64))]
struct AlignedU64Array8([u64; 8]);
pub fn fast_batched_dot_product_avx512_k1<T: Copy>(
    params: &Params,
    c: &mut [u64], // Output slice
    a: &[u64],     // Input slice 'a'
    a_elems: usize,
    b_t: &[T], // Transposed input slice 'b'
    b_rows: usize,
    b_cols: usize,
) where
    *const T: ToM512,
{
    assert_eq!(a_elems, b_rows);
    assert_eq!(a.len(), a_elems, "For K=1, 'a' should have a_elems elements.");
    assert_eq!(c.len(), b_cols, "For K=1, 'c' should have b_cols elements.");

    let mut simd_width = 8; 
    while simd_width > a_elems && simd_width > 1 {
        simd_width /= 2;
    }

    let k_chunk_size_factor = if simd_width > 0 { a_elems / simd_width } else { 0 };

    let chunk_size = (8192).min(k_chunk_size_factor);
    let num_chunks = if chunk_size > 0 { k_chunk_size_factor / chunk_size } else { 0 };

    let j_num_chunks = b_cols;

    let res_mut_slc = c;

    unsafe {
        let b_ptr = b_t.as_ptr();

        for k_outer in 0..num_chunks {
            for j_outer in 0..j_num_chunks {

                let mut total_sum_lo_val = _mm512_setzero_si512();
                let mut total_sum_hi_val = _mm512_setzero_si512();

                for k_inner in 0..chunk_size {
                    let k_base = simd_width * (k_outer * chunk_size + k_inner);

                    let a_idx = k_base;
                    let b_idx = j_outer * b_rows + k_base; 
                                            
                    let b_val_simd = b_ptr.add(b_idx).to_m512();

                    let a_m512_val = _mm512_load_si512(a.as_ptr().add(a_idx) as *const _);

                    let a_val_lo = a_m512_val;
                    let a_val_hi = _mm512_srli_epi64(a_m512_val, 32);

                    total_sum_lo_val = _mm512_add_epi64(
                        total_sum_lo_val,
                        _mm512_mul_epu32(a_val_lo, b_val_simd),
                    );
                    total_sum_hi_val = _mm512_add_epi64(
                        total_sum_hi_val,
                        _mm512_mul_epu32(a_val_hi, b_val_simd),
                    );
                } 

                let mut values_lo = AlignedU64Array8([0u64; 8]);
                let mut values_hi = AlignedU64Array8([0u64; 8]);
                _mm512_store_si512(values_lo.0.as_mut_ptr() as *mut _, total_sum_lo_val);
                _mm512_store_si512(values_hi.0.as_mut_ptr() as *mut _, total_sum_hi_val);

                let res_lo = values_lo.0.iter().sum::<u64>();
                let res_hi = values_hi.0.iter().sum::<u64>();

                // _mm512_store_si512(
                //     values_lo.as_mut_ptr() as *mut _,
                //     total_sum_lo_val,
                // );
                // _mm512_store_si512(
                //     values_hi.as_mut_ptr() as *mut _,
                //     total_sum_hi_val,
                // );

                // let res_lo = values_lo.iter().sum::<u64>();
                // let res_hi = values_hi.iter().sum::<u64>();
                
                let (lo_coeff, hi_coeff) = (
                    barrett_coeff_u64(params, res_lo, 0),
                    barrett_coeff_u64(params, res_hi, 1),
                );

                let composed_res = params.crt_compose_2(lo_coeff, hi_coeff);
                res_mut_slc[j_outer] = barrett_u64(params, res_mut_slc[j_outer] + composed_res);
            } 
        } 

    }
}

pub fn fast_batched_dot_product_avx512_k1_crt1_tester(
    params: &Params,
    c: &mut [u64], // Output slice
    a: &[u64],     // Input slice 'a'
    a_elems: usize,
    b_t: &[u16], // Transposed input slice 'b'
    b_rows: usize,
    b_cols: usize
){
    assert_eq!(a_elems, b_rows);
    for i in 0..b_cols {
        c[i] = 0;
        for j in 0..a_elems {
            let b_thing = b_t[i * b_rows + j];
            // println!("b_thing: {:?}", b_thing);
            // println!("a[j]: {:?}", a[j]);
            c[i] += multiply_uint_mod(b_thing as u64,  a[j], params.modulus);
            
            // let cast_b = value_as_u16 as u64;
            // c[i] += (cast_b * a[j]) % params.modulus;
            c[i] %= params.modulus;
        }
    }
}

#[repr(C)]
#[repr(align(64))]
struct AlignedU64Array8888([u64; 8]);
pub fn fast_batched_dot_product_avx512_k1_crt1<T: Copy>(
    params: &Params,
    c: &mut [u64], // Output slice
    a: &[u64],     // Input slice 'a'
    a_elems: usize,
    b_t: &[T], // Transposed input slice 'b'
    b_rows: usize,
    b_cols: usize,
) where
    *const T: ToM512,
{
    assert_eq!(a_elems, b_rows);
    assert_eq!(a.len(), a_elems, "For K=1, 'a' should have a_elems elements.");
    assert_eq!(c.len(), b_cols, "For K=1, 'c' should have b_cols elements.");

    let mut simd_width = 8; 
    while simd_width > a_elems && simd_width > 1 {
        simd_width /= 2;
    }

    let k_chunk_size_factor = if simd_width > 0 { a_elems / simd_width } else { 0 };

    let chunk_size = (8192).min(k_chunk_size_factor);
    let num_chunks = if chunk_size > 0 { k_chunk_size_factor / chunk_size } else { 0 };

    let j_num_chunks = b_cols;

    let res_mut_slc = c;

    unsafe {
        let b_ptr = b_t.as_ptr();

        for k_outer in 0..num_chunks {
            for j_outer in 0..j_num_chunks {

                let mut total_sum_lo_val = _mm512_setzero_si512();

                for k_inner in 0..chunk_size {
                    let k_base = simd_width * (k_outer * chunk_size + k_inner);

                    let a_idx = k_base;
                    let b_idx = j_outer * b_rows + k_base; 
                                            
                    let b_val_simd = b_ptr.add(b_idx).to_m512();

                    let a_m512_val = _mm512_load_si512(a.as_ptr().add(a_idx) as *const _);

                    let a_val_lo = a_m512_val;

                    total_sum_lo_val = _mm512_add_epi64(
                        total_sum_lo_val,
                        _mm512_mul_epu32(a_val_lo, b_val_simd),
                    );
                } 

                let mut values_lo = AlignedU64Array8888([0u64; 8]);
                _mm512_store_si512(values_lo.0.as_mut_ptr() as *mut _, total_sum_lo_val);

                let res_lo = values_lo.0.iter().sum::<u64>();

                // let mut values_lo = [0u64; 8];
                // _mm512_store_si512(
                //     values_lo.as_mut_ptr() as *mut _,
                //     total_sum_lo_val,
                // );

                // let res_lo = values_lo.iter().sum::<u64>();
                
                let lo_coeff = barrett_coeff_u64(params, res_lo, 0);

                let composed_res = params.crt_compose_1(lo_coeff);
                res_mut_slc[j_outer] = barrett_u64(params, res_mut_slc[j_outer] + composed_res);
            } 
        } 

    }
}

pub fn fast_batched_dot_product_generic<T:Copy>(
    params: &Params,
    c: &mut [u64], // Output slice
    a: &[u64],     // Input slice 'a'
    a_elems: usize,
    b_t: &[T], // Transposed input slice 'b'
    b_rows: usize,
    b_cols: usize,
) where
*const T: ToM512, {
    if params.crt_count == 1 {
        fast_batched_dot_product_avx512_k1_crt1(params, c, a, a_elems, b_t, b_rows, b_cols)
    } else {
        fast_batched_dot_product_avx512_k1(params, c, a, a_elems, b_t, b_rows, b_cols)    
    }
}

// pub fn fast_batched_dot_product_generic_u16(
//     params: &Params,
//     c: &mut [u64], // Output slice
//     a: &[u64],     // Input slice 'a'
//     a_elems: usize,
//     b_t: &[u16], // Transposed input slice 'b'
//     b_rows: usize,
//     b_cols: usize,
// ) {
//     // if params.crt_count == 1 {
//         fast_batched_dot_product_avx512_k1_crt1_tester(params, c, a, a_elems, b_t, b_rows, b_cols)
//     // } else {
//     //     fast_batched_dot_product_avx512_k1(params, c, a, a_elems, b_t, b_rows, b_cols)    
//     // }
// }

#[cfg(test)]
mod test {
    use std::time::Instant;

    use log::debug;
    use spiral_rs::aligned_memory::AlignedMemory64;
    use spiral_rs::poly::*;

    use super::super::util::test_params;
    use super::*;
    use test_log::test;

    #[test]
    fn test_fast_batched_dot_product_avx512() {
        let params = test_params();

        let a_cols = 4096;
        let b_rows = a_cols;
        let b_cols = 16;

        let mut a = PolyMatrixRaw::zero(&params, 1, a_cols / params.poly_len);
        for i in 0..a.cols {
            for j in 0..params.poly_len {
                a.get_poly_mut(0, i)[j] = 1u64;
            }
        }
        let mut b_u64 = AlignedMemory64::new(b_rows * b_cols);
        let b_u16: [u16; 65536] = [1u16; 4096 * 16];
        let mut c = AlignedMemory64::new(1 * b_cols);
        let mut c_other = AlignedMemory64::new(b_cols);
        let trials = 1;
        let mut sum = 0u64;
        let mut sum_time = 0;
        for _ in 0..trials {
            for i in 0..b_u64.len() {
                // b[i] = fastrand::u64(..);
                b_u64[i] = 1u64;
            }
            let b_u16_slc =
                unsafe { std::slice::from_raw_parts(b_u64.as_ptr() as *const u16, b_u64.len() * 4) };

            let now = Instant::now();
            // fast_dot_product_avx512(&params, a.as_slice(), rows, b_as_t_slice, rows,
            // cols); fast_dot_product_avx512_fastest(&params, a.as_slice(),
            // rows, b_as_t_slice, rows, cols);
            fast_batched_dot_product_avx512_k1::<_>(
                &params,
                c.as_mut_slice(),
                a.as_slice(),
                a_cols,
                b_u16_slc,
                b_rows,
                b_cols,
            );
            sum_time += now.elapsed().as_micros();
            sum += c.as_slice()[fastrand::usize(..c.len())];

            fast_batched_dot_product_avx512_k1_crt1_tester(
                &params,
                c_other.as_mut_slice(),
                a.as_slice(),
                a_cols,
                b_u16.as_slice(),
                b_rows,
                b_cols
            );
            assert_eq!(c.as_slice(), c_other.as_slice());
        }
        debug!("fast_matmul_avx512 in {} us ({}: {} x {})", sum_time, trials, a_cols, b_cols);
        debug!("");
        debug!("{}", sum);
    }

    use crate::params::params_for_scenario_simplepir_small;
    
    #[test]
    fn test_fast_batched_dot_product_generic() {
        // let params = test_params();
        let (params, _) = params_for_scenario_simplepir_small(1 << 30, 1);

        let a_cols = 4096;
        let b_rows = a_cols;
        let b_cols = 16;

        let a = PolyMatrixRaw::random(&params, 1, a_cols / params.poly_len);
        // for i in 0..a.cols {
        //     for j in 0..params.poly_len {
        //         a.get_poly_mut(0, i)[j] = 1u64;
        //     }
        // }

        let mut b = AlignedMemory64::new(b_rows * b_cols);
        // let mut b_u16: [u16; 65536] = [1u16; 4096 * 16];

        let mut c = AlignedMemory64::new(b_cols);
        let mut c_other = AlignedMemory64::new(b_cols);
        let trials = 1;
        let mut sum = 0u64;
        let mut sum_time = 0;
        for _ in 0..trials {
            for i in 0..b.len() {
                b[i] = fastrand::u64(..);
            }
            let b_u16_slc: &[u16] =
                unsafe { std::slice::from_raw_parts(b.as_ptr() as *const u16, b.len() * 4) };

            let now = Instant::now();
            fast_batched_dot_product_avx512_k1_crt1::<_>(
                &params,
                c.as_mut_slice(),
                a.as_slice(),
                a_cols,
                b_u16_slc,
                b_rows,
                b_cols,
            );
            sum_time += now.elapsed().as_micros();
            sum += c.as_slice()[fastrand::usize(..c.len())];

            fast_batched_dot_product_avx512_k1_crt1_tester(
                &params,
                c_other.as_mut_slice(),
                a.as_slice(),
                a_cols,
                b_u16_slc,
                b_rows,
                b_cols
            );

            assert_eq!(c.as_slice(), c_other.as_slice());
        }
        debug!("fast_matmul_avx512 in {} us ({}: {} x {})", sum_time, trials, a_cols, b_cols);
        debug!("");
        debug!("{}", sum);
    }    
}
