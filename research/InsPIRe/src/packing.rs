use std::{arch::x86_64::*, time::Instant};

use log::debug;

// use spiral_rs::poly::multiply_add_poly_avx;
use spiral_rs::{
    arith::*, discrete_gaussian::*, gadget::*, ntt::*, number_theory::invert_uint_mod, params::*,
    poly::*,
};

use crate::server::Precomp;
use rand::SeedableRng;
use rand_chacha::ChaCha20Rng;
use std::str::FromStr;

use super::{client::*, measurement::*};

use super::util::*;

pub fn gadget_invert_transposed_alloc<'a>(
    inp: &PolyMatrixRaw<'a>,
    num_digits: usize,
) -> PolyMatrixRaw<'a> {
    assert_eq!(inp.cols, 1);

    let params = inp.params;
    let mut out = PolyMatrixRaw::zero(&params, inp.rows, num_digits);

    let num_elems = num_digits;
    let bits_per = get_bits_per(params, num_elems);
    let mask = (1u64 << bits_per) - 1;

    for i in 0..inp.rows {
        for z in 0..params.poly_len {
            let val = inp.get_poly(i, 0)[z];
            for k in 0..num_elems {
                let bit_offs = k * bits_per;
                let piece = if bit_offs >= 64 { 0 } else { (val >> bit_offs) & mask };
                out.get_poly_mut(i, k)[z] = piece;
            }
        }
    }
    out
}

fn homomorphic_automorph<'a>(
    params: &'a Params,
    t: usize,
    t_exp: usize,
    ct: &PolyMatrixNTT<'a>,
    pub_param: &PolyMatrixNTT<'a>,
) -> PolyMatrixNTT<'a> {
    assert_eq!(ct.rows, 2);
    assert_eq!(ct.cols, 1);

    let ct_raw = ct.raw();
    let ct_auto = automorph_alloc(&ct_raw, t);

    let mut ginv_ct = PolyMatrixRaw::zero(params, t_exp, 1);
    gadget_invert_rdim(&mut ginv_ct, &ct_auto, 1);
    let mut ginv_ct_ntt = PolyMatrixNTT::zero(params, t_exp, 1);
    for i in 1..t_exp {
        let pol_src: &[u64] = ginv_ct.get_poly(i, 0);
        let pol_dst = ginv_ct_ntt.get_poly_mut(i, 0);
        reduce_copy(params, pol_dst, pol_src);
        ntt_forward(params, pol_dst);
    }
    // let ginv_ct_ntt = ginv_ct.ntt();
    let w_times_ginv_ct = pub_param * &ginv_ct_ntt;

    let mut ct_auto_1 = PolyMatrixRaw::zero(params, 1, 1);
    ct_auto_1.data.as_mut_slice().copy_from_slice(ct_auto.get_poly(1, 0));
    let ct_auto_1_ntt = ct_auto_1.ntt();

    &ct_auto_1_ntt.pad_top(1) + &w_times_ginv_ct
}

pub fn pack_lwes_inner<'a>(
    params: &'a Params,
    ell: usize,
    start_idx: usize,
    rlwe_cts: &[PolyMatrixNTT<'a>],
    pub_params: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> PolyMatrixNTT<'a> {
    assert_eq!(pub_params.len(), params.poly_len_log2);

    if ell == 0 {
        return rlwe_cts[start_idx].clone();
    }

    let step = 1 << (params.poly_len_log2 - ell);
    let even = start_idx;
    let odd = start_idx + step;

    let mut ct_even = pack_lwes_inner(params, ell - 1, even, rlwe_cts, pub_params, y_constants);
    let ct_odd = pack_lwes_inner(params, ell - 1, odd, rlwe_cts, pub_params, y_constants);

    let (y, neg_y) = (&y_constants.0[ell - 1], &y_constants.1[ell - 1]);

    let y_times_ct_odd = scalar_multiply_alloc(&y, &ct_odd);
    let neg_y_times_ct_odd = scalar_multiply_alloc(&neg_y, &ct_odd);

    let mut ct_sum_1 = ct_even.clone();
    add_into(&mut ct_sum_1, &neg_y_times_ct_odd);
    add_into(&mut ct_even, &y_times_ct_odd);

    // let now = Instant::now();
    let ct_sum_1_automorphed = homomorphic_automorph(
        params,
        (1 << ell) + 1,
        params.t_exp_left,
        &ct_sum_1,
        &pub_params[params.poly_len_log2 - 1 - (ell - 1)],
    );
    // debug!("Homomorphic automorph in {} us", now.elapsed().as_micros());

    &ct_even + &ct_sum_1_automorphed
}

fn pack_lwes_inner_non_recursive<'a>(
    params: &'a Params,
    ell: usize,
    _start_idx: usize,
    rlwe_cts: &[PolyMatrixNTT<'a>],
    pub_params: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
    prepared_vals: Option<&[PolyMatrixNTT<'a>]>,
    mut output_prepared_vals: Option<&mut Vec<PolyMatrixNTT<'a>>>,
) -> PolyMatrixNTT<'a> {
    assert!(pub_params.len() == params.poly_len_log2 || pub_params.len() == 0);
    assert_eq!(params.crt_count, 2);

    // let mut working_set = Vec::with_capacity(1 << (ell - 1));
    // let num_out = 1 << (ell - 1);
    // for i in 0..num_out {
    //     let combined = combine(
    //         params,
    //         1,
    //         &rlwe_cts[i],
    //         &rlwe_cts[num_out + i],
    //         pub_params,
    //         y_constants,
    //     );
    //     working_set.push(combined);
    // }

    let mut working_set = rlwe_cts.to_vec();

    let mut y_times_ct_odd = PolyMatrixNTT::zero(params, 2, 1);
    let mut neg_y_times_ct_odd = PolyMatrixNTT::zero(params, 2, 1);
    let mut ct_sum_1 = PolyMatrixNTT::zero(params, 2, 1);

    let mut ct_raw = PolyMatrixRaw::zero(params, 1, 1);
    let mut ct_auto = PolyMatrixRaw::zero(params, 1, 1);
    let mut ginv_ct = PolyMatrixRaw::zero(params, params.t_exp_left, 1);
    let mut ginv_ct_ntt = PolyMatrixNTT::zero(params, params.t_exp_left, 1);
    let mut ct_auto_1_ntt = PolyMatrixNTT::zero(params, 1, 1);
    let mut w_times_ginv_ct = PolyMatrixNTT::zero(params, 2, 1);
    let mut scratch = PolyMatrixNTT::zero(params, 2, 1);
    let scratch_mut_slc = scratch.as_mut_slice();

    let mut total_0 = 0;
    let mut total_1 = 0;
    let mut total_2 = 0;
    let mut total_3 = 0;
    let mut total_4 = 0;

    let mut num_ntts = 0;

    let mut ct_raw_1_auto_ntt;

    for cur_ell in 1..=ell {
        let num_in = 1 << (ell - cur_ell + 1);
        let num_out = num_in >> 1;

        let (first_half, second_half) = (&mut working_set[..num_in]).split_at_mut(num_out);

        for i in 0..num_out {
            let now = Instant::now();
            let ct_even = &mut first_half[i];
            let ct_odd = &second_half[i];

            let (y, neg_y) = (&y_constants.0[cur_ell - 1], &y_constants.1[cur_ell - 1]);

            // if i == 5 {
            //     debug!("neg_y: {:?}", neg_y.as_slice());
            // }

            scalar_multiply_avx(&mut y_times_ct_odd, &y, &ct_odd);
            scalar_multiply_avx(&mut neg_y_times_ct_odd, &neg_y, &ct_odd);

            // if i == 5 && cur_ell == 1 && output_prepared_vals.is_none() {
            //     debug!(
            //         "ct_even[0]: {:?}",
            //         params.crt_compose(ct_even.get_poly(1, 0), 0)
            //     );
            //     debug!(
            //         "ct_odd[0]: {:?}",
            //         params.crt_compose(ct_odd.get_poly(1, 0), 0)
            //     );
            // }

            ct_sum_1.as_mut_slice().copy_from_slice(ct_even.as_slice());
            add_into(&mut ct_sum_1, &neg_y_times_ct_odd);
            fast_add_into_no_reduce(ct_even, &y_times_ct_odd);
            total_3 += now.elapsed().as_micros();

            {
                let ct: &PolyMatrixNTT<'_> = &ct_sum_1;
                let t = (1 << cur_ell) + 1;
                let t_exp = params.t_exp_left;
                let (cur_ginv_ct_ntt, cur_ct_auto_1_ntt) = if cur_ell == 1
                    && prepared_vals.is_some()
                {
                    let ginv_ct_ntt = &prepared_vals.unwrap()[i];

                    // In this first round, this value is always zero
                    ct_raw_1_auto_ntt = PolyMatrixNTT::zero(params, 1, 1);
                    (ginv_ct_ntt, &ct_raw_1_auto_ntt)
                } else {
                    let now = Instant::now();
                    // let ct_raw = ct.raw();
                    // nb: scratch has 2nd row of ct in uncrtd form,
                    //     ct_raw has only first row
                    from_ntt_scratch(&mut ct_raw, scratch_mut_slc, ct);
                    if cur_ell == 1 {
                        num_ntts += 2;
                    }
                    total_0 += now.elapsed().as_micros();
                    let now = Instant::now();
                    automorph(&mut ct_auto, &ct_raw, t);
                    total_1 += now.elapsed().as_micros();

                    gadget_invert_rdim(&mut ginv_ct, &ct_auto, 1);

                    let skip_first_gadget_dim = true;
                    if skip_first_gadget_dim {
                        for i in 1..t_exp {
                            let pol_src = ginv_ct.get_poly(i, 0);
                            let pol_dst = ginv_ct_ntt.get_poly_mut(i, 0);
                            pol_dst[..params.poly_len].copy_from_slice(pol_src);
                            pol_dst[params.poly_len..].copy_from_slice(pol_src);

                            ntt_forward(params, pol_dst);
                            if cur_ell == 1 {
                                num_ntts += 1;
                            }
                        }
                    } else {
                        to_ntt(&mut ginv_ct_ntt, &ginv_ct);
                        // num_ntts += ginv_ct_ntt.rows * ginv_ct_ntt.cols;
                    }

                    let now = Instant::now();
                    automorph_poly_uncrtd(params, ct_auto_1_ntt.as_mut_slice(), scratch_mut_slc, t);
                    ntt_forward(params, ct_auto_1_ntt.as_mut_slice());
                    // num_ntts += 1;

                    total_4 += now.elapsed().as_micros();

                    (&ginv_ct_ntt, &ct_auto_1_ntt)
                };

                if output_prepared_vals.is_some() {
                    let opv_mut = output_prepared_vals.as_deref_mut();
                    opv_mut.unwrap().push(cur_ginv_ct_ntt.clone());
                    continue;
                }

                let pub_param = &pub_params[params.poly_len_log2 - 1 - (cur_ell - 1)];
                // let ginv_ct_ntt = ginv_ct.ntt();
                // let w_times_ginv_ct = pub_param * &ginv_ct_ntt;
                w_times_ginv_ct.as_mut_slice().fill(0);
                multiply_no_reduce(&mut w_times_ginv_ct, &pub_param, &cur_ginv_ct_ntt, 1);

                // &ct_auto_1_ntt.pad_top(1) + &w_times_ginv_ct
                let now = Instant::now();
                add_into_at_no_reduce(ct_even, &cur_ct_auto_1_ntt, 1, 0);
                add_into(ct_even, &w_times_ginv_ct);
                total_2 += now.elapsed().as_micros();
            };
        }

        if output_prepared_vals.is_some() {
            return PolyMatrixNTT::zero(params, 2, 1);
        }
    }

    if false {
        debug!("num_ntts: {}", num_ntts);
        debug!("total_0: {} us", total_0);
        debug!("total_1: {} us", total_1);
        debug!("total_2: {} us", total_2);
        debug!("total_3: {} us", total_3);
        debug!("total_4: {} us", total_4);
    }

    // let mut res = PolyMatrixNTT::zero(params, 2, 1);
    // // let mut res = working_set[0].clone();
    // add_into(&mut res, &working_set[0]);
    // res

    working_set[0].clone()
}

pub fn precompute_pack<'a>(
    params: &'a Params,
    ell: usize,
    rlwe_cts: &[PolyMatrixNTT<'a>],
    fake_pub_params: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> (PolyMatrixNTT<'a>, Vec<PolyMatrixNTT<'a>>, Vec<Vec<usize>>) {
    assert!(fake_pub_params.len() == params.poly_len_log2);
    assert_eq!(params.crt_count, 2);

    let mut working_set = rlwe_cts.to_vec();

    let mut y_times_ct_odd = PolyMatrixNTT::zero(params, 2, 1);
    let mut neg_y_times_ct_odd = PolyMatrixNTT::zero(params, 2, 1);
    let mut ct_sum_1 = PolyMatrixNTT::zero(params, 2, 1);

    let mut ct_raw = PolyMatrixRaw::zero(params, 1, 1);
    let mut ct_auto = PolyMatrixRaw::zero(params, 1, 1);
    let mut ginv_ct = PolyMatrixRaw::zero(params, params.t_exp_left, 1);
    let mut ginv_ct_ntt = PolyMatrixNTT::zero(params, params.t_exp_left, 1);
    let mut ct_auto_1_ntt = PolyMatrixNTT::zero(params, 1, 1);
    let mut w_times_ginv_ct = PolyMatrixNTT::zero(params, 2, 1);
    let mut scratch = PolyMatrixNTT::zero(params, 2, 1);
    let scratch_mut_slc = scratch.as_mut_slice();

    let mut total_0 = 0;
    let mut total_1 = 0;
    let mut total_2 = 0;
    let mut total_3 = 0;
    let mut total_4 = 0;

    let mut num_ntts = 0;

    let mut res = Vec::new();

    for cur_ell in 1..=ell {
        let num_in = 1 << (ell - cur_ell + 1);
        let num_out = num_in >> 1;

        let (first_half, second_half) = (&mut working_set[..num_in]).split_at_mut(num_out);

        for i in 0..num_out {
            let now = Instant::now();
            let ct_even = &mut first_half[i];
            let ct_odd = &second_half[i];

            let (y, neg_y) = (&y_constants.0[cur_ell - 1], &y_constants.1[cur_ell - 1]);

            scalar_multiply_avx(&mut y_times_ct_odd, &y, &ct_odd);
            scalar_multiply_avx(&mut neg_y_times_ct_odd, &neg_y, &ct_odd);

            ct_sum_1.as_mut_slice().copy_from_slice(ct_even.as_slice());
            add_into(&mut ct_sum_1, &neg_y_times_ct_odd);
            fast_add_into_no_reduce(ct_even, &y_times_ct_odd);
            total_3 += now.elapsed().as_micros();

            {
                let ct: &PolyMatrixNTT<'_> = &ct_sum_1;
                let t = (1 << cur_ell) + 1;
                let t_exp = params.t_exp_left;
                let (cur_ginv_ct_ntt, cur_ct_auto_1_ntt) = {
                    let now = Instant::now();
                    // let ct_raw = ct.raw();

                    // nb: scratch has 2nd row of ct in uncrtd form,
                    //     ct_raw has only first row
                    from_ntt_scratch(&mut ct_raw, scratch_mut_slc, ct);
                    if cur_ell == 1 {
                        num_ntts += 2;
                    }
                    total_0 += now.elapsed().as_micros();
                    let now = Instant::now();
                    automorph(&mut ct_auto, &ct_raw, t);
                    total_1 += now.elapsed().as_micros();

                    gadget_invert_rdim(&mut ginv_ct, &ct_auto, 1);

                    let skip_first_gadget_dim = false;
                    if skip_first_gadget_dim {
                        for i in 1..t_exp {
                            let pol_src = ginv_ct.get_poly(i, 0);
                            let pol_dst = ginv_ct_ntt.get_poly_mut(i, 0);
                            pol_dst[..params.poly_len].copy_from_slice(pol_src);
                            pol_dst[params.poly_len..].copy_from_slice(pol_src);

                            ntt_forward(params, pol_dst);
                            if cur_ell == 1 {
                                num_ntts += 1;
                            }
                        }
                    } else {
                        to_ntt(&mut ginv_ct_ntt, &ginv_ct);
                        // num_ntts += ginv_ct_ntt.rows * ginv_ct_ntt.cols;
                    }

                    let now = Instant::now();
                    automorph_poly_uncrtd(params, ct_auto_1_ntt.as_mut_slice(), scratch_mut_slc, t);
                    ntt_forward(params, ct_auto_1_ntt.as_mut_slice());
                    // num_ntts += 1;

                    total_4 += now.elapsed().as_micros();

                    (&ginv_ct_ntt, &ct_auto_1_ntt)
                };

                // println!(
                //     "ct_auto_1_ntt.raw(): {:?}",
                //     &ct_auto_1_ntt.raw().as_slice()[..30]
                // );

                res.push(condense_matrix(params, cur_ginv_ct_ntt));

                let pub_param = &fake_pub_params[params.poly_len_log2 - 1 - (cur_ell - 1)];
                // let ginv_ct_ntt = ginv_ct.ntt();
                // let w_times_ginv_ct = pub_param * &ginv_ct_ntt;
                w_times_ginv_ct.as_mut_slice().fill(0);
                multiply_no_reduce(&mut w_times_ginv_ct, &pub_param, &cur_ginv_ct_ntt, 0);

                // &ct_auto_1_ntt.pad_top(1) + &w_times_ginv_ct
                let now = Instant::now();
                add_into_at_no_reduce(ct_even, &cur_ct_auto_1_ntt, 1, 0);
                add_into(ct_even, &w_times_ginv_ct);
                total_2 += now.elapsed().as_micros();
            };
        }
    }

    if false {
        debug!("num_ntts: {}", num_ntts);
        debug!("total_0: {} us", total_0);
        debug!("total_1: {} us", total_1);
        debug!("total_2: {} us", total_2);
        debug!("total_3: {} us", total_3);
        debug!("total_4: {} us", total_4);
    }

    let tables = generate_automorph_tables_brute_force(&params);

    (working_set[0].clone(), res, tables)
}

// pub fn fake_precompute_pack<'a>(
//     params: &'a Params,
//     ell: usize,
//     rlwe_cts: &[PolyMatrixNTT<'a>],
//     fake_pub_params: &[PolyMatrixNTT<'a>],
//     y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
// ) {
//     let tables = generate_automorph_tables_brute_force(&params);

//     (working_set[0].clone(), res, tables)    
// }

pub fn pack_using_precomp_vals<'a>(
    params: &'a Params,
    ell: usize,
    pub_params: &[PolyMatrixNTT<'a>],
    b_values: &[u64],
    precomp_res: &PolyMatrixNTT<'a>,
    precomp_vals: &[PolyMatrixNTT<'a>],
    precomp_tables: &[Vec<usize>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> PolyMatrixRaw<'a> {
    // let now = Instant::now();
    // let mut working_set = vec![PolyMatrixNTT::zero(params, 1, 1); 1 << ell];
    let mut working_set = Vec::with_capacity(1 << (ell - 1));
    for _ in 0..(1 << (ell - 1)) {
        working_set.push(PolyMatrixNTT::zero(params, 1, 1));
    }

    let mut y_times_ct_odd = PolyMatrixNTT::zero(params, 1, 1);
    let mut neg_y_times_ct_odd = PolyMatrixNTT::zero(params, 1, 1);
    let mut ct_sum_1 = PolyMatrixNTT::zero(params, 1, 1);
    let mut w_times_ginv_ct = PolyMatrixNTT::zero(params, 1, 1);

    // println!("time_-1: {} us", now.elapsed().as_micros());

    let mut time_0 = 0;
    let mut time_1 = 0;
    let mut time_2 = 0;
    let mut time_3 = 0;
    let mut time_4 = 0;

    let mut idx_precomp = 0;
    let mut num_muls = 0;
    for cur_ell in 1..=ell {
        let mut num_in = 1 << (ell - cur_ell + 1);
        let num_out = num_in >> 1;

        if num_in == params.poly_len {
            num_in = num_out;
        }

        let (first_half, second_half) = (&mut working_set[..num_in]).split_at_mut(num_out);

        for i in 0..num_out {
            let now = Instant::now();
            let ct_even = &mut first_half[i];

            let (y, neg_y) = (&y_constants.0[cur_ell - 1], &y_constants.1[cur_ell - 1]);

            if cur_ell > 1 {
                let ct_odd = &mut second_half[i];
                scalar_multiply_avx(&mut y_times_ct_odd, &y, &ct_odd);
                scalar_multiply_avx(&mut neg_y_times_ct_odd, &neg_y, &ct_odd);
            }

            time_0 += now.elapsed().as_micros();

            let now = Instant::now();
            if cur_ell > 1 {
                ct_sum_1.as_mut_slice().copy_from_slice(ct_even.as_slice());
                fast_add_into_no_reduce(&mut ct_sum_1, &neg_y_times_ct_odd);
                fast_add_into_no_reduce(ct_even, &y_times_ct_odd);
            }
            time_1 += now.elapsed().as_micros();

            // --

            let now = Instant::now();
            let ct: &PolyMatrixNTT<'_> = &ct_sum_1;
            let t = (1 << cur_ell) + 1;

            let cur_ginv_ct_ntt = &precomp_vals[idx_precomp];
            idx_precomp += 1;

            let w = &pub_params[params.poly_len_log2 - 1 - (cur_ell - 1)];
            // let w = pub_param.submatrix(1, 0, 1, pub_param.cols);
            // let w_times_ginv_ct = &w * cur_ginv_ct_ntt;
            // multiply(&mut w_times_ginv_ct, &w, &cur_ginv_ct_ntt);
            // w_times_ginv_ct.as_mut_slice().fill(0);
            fast_multiply_no_reduce(params, &mut w_times_ginv_ct, &w, &cur_ginv_ct_ntt);
            num_muls += 1;
            time_2 += now.elapsed().as_micros();

            if cur_ell > 1 {
                let now = Instant::now();
                apply_automorph_ntt(params, &precomp_tables, &ct, ct_even, t);

                // fast_add_into_no_reduce(ct_even, &ct_auto_1_ntt);
                time_3 += now.elapsed().as_micros();
                let now = Instant::now();

                // second condition prevents overflow
                if i < num_out / 2 && ((cur_ell - 1) % 5 != 0) {
                    fast_add_into_no_reduce(ct_even, &w_times_ginv_ct);
                } else {
                    // reduction right before or after addition is much faster than at
                    // multiplication time
                    fast_add_into(ct_even, &w_times_ginv_ct);
                }
                time_4 += now.elapsed().as_micros();
            } else {
                let now = Instant::now();
                if i < num_out / 2 {
                    fast_add_into_no_reduce(ct_even, &w_times_ginv_ct);
                } else {
                    fast_add_into(ct_even, &w_times_ginv_ct);
                }
                time_4 += now.elapsed().as_micros();
            }
        }
    }
    // let now = Instant::now();

    if false {
        println!("time_0: {} us", time_0);
        println!("time_1: {} us", time_1);
        println!("time_2: {} us", time_2);
        println!("time_3: {} us", time_3);
        println!("time_4: {} us", time_4);
        println!("idx_precomp: {}", idx_precomp);
        println!("num_muls: {}", num_muls);
    }

    assert_eq!(idx_precomp, precomp_vals.len());

    let mut resulting_row_1 = working_set[0].clone();
    fast_reduce(&mut resulting_row_1);

    let resulting_row_1 = resulting_row_1.as_slice();

    let mut res = precomp_res.clone();
    // {
    //     let r = res.raw();
    //     println!("res row 0: {:?}", &r.get_poly(0, 0)[..30]);
    //     println!("res row 1: {:?}", &r.get_poly(1, 0)[..30]);
    // }
    res.get_poly_mut(1, 0).copy_from_slice(resulting_row_1);

    // println!(
    //     "precomp_res     row 1: {:?}",
    //     &precomp_res.raw().get_poly(1, 0)[..30]
    // );
    // println!(
    //     "resulting_row_1 row 1: {:?}",
    //     &working_set[0].raw().get_poly(1, 0)[..30]
    // );

    // let mut res = precomp_res.clone();

    let mut out_raw = res.raw();
    for z in 0..params.poly_len {
        let val = barrett_reduction_u128(params, b_values[z] as u128 * params.poly_len as u128);
        let idx = params.poly_len + z;
        out_raw.data[idx] += val;
        if out_raw.data[idx] >= params.modulus {
            out_raw.data[idx] -= params.modulus;
        }
    }

    // println!("time_5: {} us", now.elapsed().as_micros());

    // for z in 0..params.poly_len {
    //     let b_value = b_values[z];
    //     let val = barrett_u64(params, res.get_poly(1, 0)[z] + b_value);
    //     res.get_poly_mut(1, 0)[z] = val;
    // }

    out_raw
}

pub fn pack_single_lwe<'a>(
    params: &'a Params,
    pub_params: &[PolyMatrixNTT<'a>],
    lwe_ct: &PolyMatrixNTT<'a>,
) -> PolyMatrixNTT<'a> {
    // computing:
    // r0 = f
    // r1 = r0 + automorph(r0, ts[0])
    // r2 = r1 + automorph(r1, ts[1])
    // ...
    // r_\log d = ...

    let mut cur_r = lwe_ct.clone();
    for i in 0..params.poly_len_log2 {
        let t = (params.poly_len / (1 << i)) + 1;
        let pub_param = &pub_params[i];
        let tau_of_r = homomorphic_automorph(params, t, params.t_exp_left, &cur_r, pub_param);
        add_into(&mut cur_r, &tau_of_r);
    }
    cur_r
}

pub fn add_all_rotations<'a>(
    params: &'a Params,
    pub_params: &[PolyMatrixNTT<'a>],
    lwe_ct: &PolyMatrixNTT<'a>,
) -> PolyMatrixNTT<'a> {
    // computing:
    // r0 = f
    // r1 = r0 + automorph(r0, ts[0])
    // r2 = r1 + automorph(r1, ts[1])
    // ...
    // r_\log d = ...

    let mut cur_r = lwe_ct.clone();
    for i in 0..params.poly_len_log2 {
        let t = (params.poly_len / (1 << i)) + 1;
        let pub_param = &pub_params[i];
        let tau_of_r = homomorphic_automorph(params, t, params.t_exp_left, &cur_r, pub_param);
        add_into(&mut cur_r, &tau_of_r);
    }
    cur_r
}

pub fn fast_barrett_raw_u64(input: u64, const_ratio_1: u64, modulus: u64) -> u64 {
    let tmp = (((input as u128) * (const_ratio_1 as u128)) >> 64) as u64;

    // Barrett subtraction
    let res = input - tmp * modulus;

    res
}

pub fn fast_add_into(res: &mut PolyMatrixNTT, a: &PolyMatrixNTT) {
    assert!(res.rows == a.rows);
    assert!(res.cols == a.cols);

    let params = res.params;
    for i in 0..res.rows {
        for j in 0..res.cols {
            let res_poly = res.get_poly_mut(i, j);
            let a_poly = a.get_poly(i, j);
            for c in 0..params.crt_count {
                for i in 0..params.poly_len {
                    let idx = c * params.poly_len + i;
                    unsafe {
                        let p_res = res_poly.as_mut_ptr().add(idx);
                        let p_a = a_poly.as_ptr().add(idx);
                        let val = *p_res + *p_a;
                        let reduced =
                            fast_barrett_raw_u64(val, params.barrett_cr_1[c], params.moduli[c]);
                        *p_res = reduced;
                    }
                }
            }
        }
    }
}

pub fn fast_multiply_no_reduce(
    params: &Params,
    res: &mut PolyMatrixNTT,
    a: &PolyMatrixNTT,
    b: &PolyMatrixNTT,
) {
    // assert_eq!(res.rows, a.rows);
    // assert_eq!(res.cols, b.cols);
    // assert_eq!(res.rows, 1);
    // assert_eq!(res.cols, 1);

    // assert_eq!(a.cols, b.rows);
    assert!(a.cols == 1 || a.rows == 1);
    assert!(b.cols == 1 || b.rows == 1);
    assert_eq!(a.cols * a.rows, b.cols * b.rows);

    // assert_eq!(params.crt_count * params.poly_len, 2 * 2048);
    // assert_eq!((params.crt_count * params.poly_len) % 1024, 0);
    assert_eq!(params.crt_count, 2);
    assert_eq!(params.poly_len, 2048);

    unsafe {
        let a_ptr = a.as_slice().as_ptr();
        let b_ptr = b.as_slice().as_ptr();
        let res_ptr = res.as_mut_slice().as_mut_ptr();
        let pol_sz = params.poly_len;

        for idx in (0..pol_sz).step_by(8) {
            let mut sum_lo = _mm512_setzero_si512();
            let mut sum_hi = _mm512_setzero_si512();
            for k in 0..a.cols {
                let p_x = a_ptr.add(k * 2 * pol_sz + idx);
                let p_y = b_ptr.add(k * 2 * pol_sz + idx);

                let x = _mm512_load_si512(p_x as *const _);
                let x_lo = x;
                let x_hi = _mm512_srli_epi64(x, 32);
                let y = _mm512_load_si512(p_y as *const _);
                let y_lo = y;
                let y_hi = _mm512_srli_epi64(y, 32);

                let product_lo = _mm512_mul_epu32(x_lo, y_lo);
                let product_hi = _mm512_mul_epu32(x_hi, y_hi);

                sum_lo = _mm512_add_epi64(sum_lo, product_lo);
                sum_hi = _mm512_add_epi64(sum_hi, product_hi);
            }

            let p_z = res_ptr.add(idx);
            _mm512_store_si512(p_z as *mut _, sum_lo);
            let p_z = res_ptr.add(pol_sz + idx);
            _mm512_store_si512(p_z as *mut _, sum_hi);
        }
    }
}

pub fn fast_multiply_no_reduce_in_range(
    params: &Params,
    res: &mut PolyMatrixNTT,
    a: &PolyMatrixNTT,
    b: &PolyMatrixNTT,
    begin_a: usize,
    begin_b: usize,
    length: usize,
) {
    // assert_eq!(res.rows, a.rows);
    // assert_eq!(res.cols, b.cols);
    // assert_eq!(res.rows, 1);
    // assert_eq!(res.cols, 1);

    // assert_eq!(a.cols, b.rows);
    // assert!(a.cols == 1 || a.rows == 1);
    // assert!(b.cols == 1 || b.rows == 1);
    //assert_eq!(a.cols * a.rows, b.cols * b.rows);

    assert_eq!(params.crt_count, 2);
    assert_eq!(params.poly_len % 1024, 0);

    unsafe {
        let a_ptr = a.as_slice().as_ptr();
        let b_ptr = b.as_slice().as_ptr();
        let res_ptr = res.as_mut_slice().as_mut_ptr();
        let pol_sz = params.poly_len;

        for idx in (0..pol_sz).step_by(8) {
            let mut sum_lo = _mm512_setzero_si512();
            let mut sum_hi = _mm512_setzero_si512();
            // for k in 0..(a.cols * a.rows) {
            for k in 0..length {
                let p_x = a_ptr.add((begin_a + k) * 2 * pol_sz + idx);
                let p_y = b_ptr.add((begin_b + k) * 2 * pol_sz + idx);

                let x = _mm512_load_si512(p_x as *const _);
                let x_lo = x;
                let x_hi = _mm512_srli_epi64(x, 32);
                let y = _mm512_load_si512(p_y as *const _);
                let y_lo = y;
                let y_hi = _mm512_srli_epi64(y, 32);

                let product_lo = _mm512_mul_epu32(x_lo, y_lo);
                let product_hi = _mm512_mul_epu32(x_hi, y_hi);

                sum_lo = _mm512_add_epi64(sum_lo, product_lo);
                sum_hi = _mm512_add_epi64(sum_hi, product_hi);
            }

            let p_z = res_ptr.add(idx);
            _mm512_store_si512(p_z as *mut _, sum_lo);
            let p_z = res_ptr.add(pol_sz + idx);
            _mm512_store_si512(p_z as *mut _, sum_hi);
        }
    }
}

pub fn fast_multiply_no_reduce_in_range_generic(
    params: &Params,
    res: &mut PolyMatrixNTT,
    a: &PolyMatrixNTT,
    b: &PolyMatrixNTT,
    begin_a: usize,
    begin_b: usize,
    length: usize,
){
    if params.crt_count == 1 {
        let a_ptr = a.as_slice();
        let b_ptr = b.as_slice();
        let res_ptr = res.as_mut_slice();
        multiply_poly_avx_in_range(params, res_ptr, a_ptr, b_ptr, begin_a, begin_b, length);
        
    } else  { //if params.crt_count == 2 {
        fast_multiply_no_reduce_in_range(params, res, a, b, begin_a, begin_b, length);
    }
}

pub fn condense_matrix<'a>(params: &'a Params, a: &PolyMatrixNTT<'a>) -> PolyMatrixNTT<'a> {
    if params.crt_count == 2 {
        let mut res = PolyMatrixNTT::zero(params, a.rows, a.cols);
        for i in 0..a.rows {
            for j in 0..a.cols {
                let res_poly = &mut res.get_poly_mut(i, j);
                let a_poly = a.get_poly(i, j);
                for z in 0..params.poly_len {
                    res_poly[z] = a_poly[z] | (a_poly[z + params.poly_len] << 32); // TODO: Does this assume that the digit decomposition is exactly 2?
                }
            }
        }
        res
    } else {
        a.clone()
    }
}

pub fn uncondense_matrix<'a>(params: &'a Params, a: &PolyMatrixNTT<'a>) -> PolyMatrixNTT<'a> {
    let mut res = PolyMatrixNTT::zero(params, a.rows, a.cols);
    for i in 0..a.rows {
        for j in 0..a.cols {
            let res_poly = &mut res.get_poly_mut(i, j);
            let a_poly = a.get_poly(i, j);
            for z in 0..params.poly_len {
                res_poly[z] = a_poly[z] & ((1u64 << 32) - 1);
                res_poly[z + params.poly_len] = a_poly[z] >> 32;
            }
        }
    }
    res
}

pub fn multiply_poly_avx(_params: &Params, res: &mut [u64], a: &[u64], b: &[u64]) {
    unsafe {
        let a_ptr = a.as_ptr();
        let b_ptr = b.as_ptr();
        let res_ptr = res.as_mut_ptr();

        for i in (0..res.len()).step_by(8) {
            let p_x = a_ptr.add(i);
            let p_y = b_ptr.add(i);
            let p_z = res_ptr.add(i);

            let x = _mm512_load_si512(p_x as *const _);
            let y = _mm512_load_si512(p_y as *const _);

            let product = _mm512_mul_epu32(x, y);

            _mm512_store_si512(p_z as *mut _, product);
        }
    }
}

pub fn multiply_poly_avx_in_range(_params: &Params, res: &mut [u64], a: &[u64], b: &[u64], begin_a: usize, begin_b: usize, length: usize) {
    unsafe {
        let a_ptr = a.as_ptr();
        let b_ptr = b.as_ptr();
        let res_ptr = res.as_mut_ptr();
        let poly_sz = _params.poly_len;

        for i in (0..poly_sz).step_by(8) {
            let p_z = res_ptr.add(i);
            let mut sum = _mm512_setzero_si512();
            for k in 0..length {
                let p_x = a_ptr.add((begin_a + k) * poly_sz + i);
                let p_y = b_ptr.add((begin_b + k) * poly_sz + i);

                let x = _mm512_load_si512(p_x as *const _);
                let y = _mm512_load_si512(p_y as *const _);

                let product = _mm512_mul_epu32(x, y);
                sum = _mm512_add_epi64(sum, product);
            }
            _mm512_store_si512(p_z as *mut _, sum);
        }
    }
}


pub fn fast_add_into_no_reduce(res: &mut PolyMatrixNTT, a: &PolyMatrixNTT) {
    assert!(res.rows == a.rows);
    assert!(res.cols == a.cols);

    let a_slc = a.as_slice();
    let res_slc = res.as_mut_slice();
    for (res_chunk, a_chunk) in res_slc.chunks_exact_mut(8).zip(a_slc.chunks_exact(8)) {
        for i in 0..8 {
            res_chunk[i] += a_chunk[i];
        }
    }
}

pub fn fast_reduce(res: &mut PolyMatrixNTT) {
    let params = res.params;
    let res_slc = res.as_mut_slice();
    for m in 0..params.crt_count {
        for i in 0..params.poly_len {
            let idx = m * params.poly_len + i;
            // res_slc[idx] = barrett_coeff_u64(params, res_slc[idx], m);
            unsafe {
                let p = res_slc.as_mut_ptr().add(idx);
                *p = barrett_coeff_u64(params, *p, m);
            }
        }
    }
}

pub fn combine<'a>(
    params: &'a Params,
    cur_ell: usize,
    ct_even: &mut PolyMatrixNTT<'a>,
    ct_odd: &PolyMatrixNTT<'a>,
    pub_params: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) {
    let (y, neg_y) = (&y_constants.0[cur_ell - 1], &y_constants.1[cur_ell - 1]);

    let y_times_ct_odd = scalar_multiply_alloc(&y, &ct_odd);
    let neg_y_times_ct_odd = scalar_multiply_alloc(&neg_y, &ct_odd);

    let mut ct_sum_1 = ct_even.clone();
    add_into(&mut ct_sum_1, &neg_y_times_ct_odd);
    add_into(ct_even, &y_times_ct_odd);

    let ct_sum_1_automorphed = homomorphic_automorph(
        params,
        (1 << cur_ell) + 1,
        params.t_exp_left,
        &ct_sum_1,
        &pub_params[params.poly_len_log2 - 1 - (cur_ell - 1)],
    );

    add_into(ct_even, &ct_sum_1_automorphed);
}

pub fn prep_pack_lwes<'a>(
    params: &'a Params,
    lwe_cts: &[u64],
    // packed_per_rlwe: usize,
) -> Vec<PolyMatrixNTT<'a>> {
    let lwe_cts_size = params.poly_len * (params.poly_len + 1);
    assert_eq!(lwe_cts.len(), lwe_cts_size);

    // assert!(cols_to_do == params.poly_len);

    let mut rlwe_cts: Vec<PolyMatrixNTT<'_>> = Vec::new();
    for i in 0..params.poly_len {
        let mut rlwe_ct = PolyMatrixRaw::zero(params, 2, 1);

        // 'a' vector
        // put this in negacyclic order
        let mut poly = Vec::new();
        for j in 0..params.poly_len {
            poly.push(lwe_cts[j * params.poly_len + i])
        }
        let nega = negacyclic_perm(&poly, 0, params.modulus);

        for j in 0..params.poly_len {
            rlwe_ct.get_poly_mut(0, 0)[j] = nega[j];
        }
        // 'b' scalar (skip)

        rlwe_cts.push(rlwe_ct.ntt());
    }

    rlwe_cts
}

// (d+1) * t * d
// v = (0 -> d) + (td -> td+d) + (2td -> 2td+d) ... + ()

pub fn prep_pack_many_lwes<'a>(
    params: &'a Params,
    lwe_cts: &[u64],
    num_rlwe_outputs: usize, // num_rlwe_outputs = db_cols / gamma
    gamma: usize,
) -> Vec<Vec<PolyMatrixNTT<'a>>> {
    let lwe_cts_size = (params.poly_len + 1) * (num_rlwe_outputs * gamma); // = (params.poly_len + 1) * db_cols
    assert_eq!(lwe_cts.len(), lwe_cts_size);

    let mut res = Vec::new();
    let ratio = params.poly_len / gamma;
    for i in 0..(num_rlwe_outputs / ratio) {
        // db_cols / params.poly_len

        let mut v = Vec::new();
        for j in 0..params.poly_len + 1 {
            let index = j * ((num_rlwe_outputs / ratio) * params.poly_len) + i * params.poly_len; // = j * db_cols + i * params.poly_len
            v.extend(&lwe_cts[index..index + params.poly_len]);
        }
        let to_push = prep_pack_lwes(params, &v);
        for k in 0..ratio {
            res.push(to_push[k * gamma..(k + 1) * gamma].to_vec());
        }
    }

    res
}

pub fn prepare_packed_vals_pack_lwes<'a>(
    params: &'a Params,
    preped_rlwe_cts: &[PolyMatrixNTT<'a>],
    _cols_to_do: usize,
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> Vec<PolyMatrixNTT<'a>> {
    let now = Instant::now();
    let mut output_preped_packed_vals = Vec::new();
    pack_lwes_inner_non_recursive(
        params,
        params.poly_len_log2,
        0,
        &preped_rlwe_cts,
        &[],
        y_constants,
        None,
        Some(&mut output_preped_packed_vals),
    );
    debug!("prepack: {} us", now.elapsed().as_micros());
    output_preped_packed_vals
}

/// Returns the `prep_packed_vals` value.
pub fn prep_pack_many_lwes_packed_vals<'a>(
    params: &'a Params,
    prep_rlwe_cts: &[Vec<PolyMatrixNTT<'a>>],
    num_rlwe_outputs: usize,
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> Vec<Vec<PolyMatrixNTT<'a>>> {
    assert_eq!(prep_rlwe_cts.len(), num_rlwe_outputs);
    assert_eq!(prep_rlwe_cts[0].len(), params.poly_len);

    let mut res = Vec::new();
    for i in 0..num_rlwe_outputs {
        res.push(prepare_packed_vals_pack_lwes(
            params,
            &prep_rlwe_cts[i],
            params.poly_len,
            y_constants,
        ));
    }

    res
}

pub fn pack_lwes<'a>(
    params: &'a Params,
    b_values: &[u64],
    preped_rlwe_cts: &[PolyMatrixNTT<'a>],
    preped_packed_vals: &[PolyMatrixNTT<'a>],
    cols_to_do: usize,
    pub_params: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> PolyMatrixNTT<'a> {
    assert_eq!(preped_rlwe_cts.len(), cols_to_do);
    assert_eq!(cols_to_do, params.poly_len);
    assert_eq!(b_values.len(), params.poly_len);

    let now = Instant::now();
    let preped_packed_val_opt =
        if preped_packed_vals.len() == 0 { None } else { Some(preped_packed_vals) };
    let out = pack_lwes_inner_non_recursive(
        params,
        params.poly_len_log2,
        0,
        &preped_rlwe_cts,
        pub_params,
        y_constants,
        preped_packed_val_opt,
        None,
    );
    let mut out_raw = out.raw();
    for z in 0..params.poly_len {
        let val = barrett_reduction_u128(params, b_values[z] as u128 * params.poly_len as u128);
        out_raw.get_poly_mut(1, 0)[z] = barrett_u64(params, out_raw.get_poly(1, 0)[z] + val);
    }
    let res = out_raw.ntt();
    debug!("True pack took {} us", now.elapsed().as_micros());

    res
}

pub fn pack_many_lwes<'a>(
    params: &'a Params,
    // prep_rlwe_cts: &[Vec<PolyMatrixNTT<'a>>],
    precomp: &Precomp<'a>,
    b_values: &[u64],
    num_rlwe_outputs: usize,
    pack_pub_params_row_1s: &[PolyMatrixNTT<'a>],
    y_constants: &(Vec<PolyMatrixNTT<'a>>, Vec<PolyMatrixNTT<'a>>),
) -> Vec<PolyMatrixRaw<'a>> {
    // assert_eq!(prep_rlwe_cts.len(), num_rlwe_outputs);
    // assert_eq!(prep_rlwe_cts[0].len(), params.poly_len);
    assert_eq!(b_values.len(), num_rlwe_outputs * params.poly_len);

    let mut res = Vec::new();
    for i in 0..num_rlwe_outputs {
        let (precomp_res, precomp_vals, precomp_tables) = &precomp[i];

        let packed = pack_using_precomp_vals(
            &params,
            params.poly_len_log2,
            &pack_pub_params_row_1s,
            &b_values[i * params.poly_len..(i + 1) * params.poly_len],
            &precomp_res,
            &precomp_vals,
            &precomp_tables,
            &y_constants,
        );

        res.push(packed);
    }

    res
}

pub fn pack_many_lwes_inspir<'a>(
    packing_params: &'a PackParams,
    precomp_inspir_vec: &Vec<PrecompInsPIR<'a>>,
    b_values: &[u64],
    packing_keys: &PackingKeys<'a>,
    gamma: usize,
) -> Vec<PolyMatrixRaw<'a>> {
    // assert!(b_values.len() >= num_rlwe_outputs * packed_per_rlwe);
    let num_rlwe_outputs = (b_values.len() as f64 / gamma as f64).ceil() as usize;

    let mut res = Vec::with_capacity(num_rlwe_outputs);
    let group_size = packing_params.params.poly_len / gamma;
    for i in 0..num_rlwe_outputs {
        let which_group = i / group_size;
        let within_group = i % group_size;

        let mut b_poly = PolyMatrixRaw::zero(&packing_params.params, 1, 1);
        // let start = which_group * packing_params.params.poly_len + within_group * gamma;
        // let end = min(start+gamma, b_values.len());
        // b_poly.data.copy_from_slice(&b_values[start..end]); 

        for j in 0..gamma {
            let index = which_group * packing_params.params.poly_len + within_group * gamma + j;
            if index >= b_values.len() {
                continue;
            }
            b_poly.get_poly_mut(0, 0)[j] = b_values[index];
        }

        // let now = Instant::now();
        let packed = if gamma <= packing_params.params.poly_len / 2 {
            
            let y_all_condensed = packing_keys.y_all_condensed.as_ref().unwrap();
            packing_with_preprocessing_online(
                &packing_params,
                &precomp_inspir_vec[i],
                &b_poly,
                &y_all_condensed,
            )
        } else {

            let y_all_condensed = packing_keys.y_all_condensed.as_ref().unwrap();
            let y_bar_all_condensed = packing_keys.y_bar_all_condensed.as_ref().unwrap();
            let z_body_condensed = packing_keys.z_body_condensed.as_ref().unwrap();            

            full_packing_with_preprocessing_online(
                &packing_params,
                &precomp_inspir_vec[i],
                &b_poly,
                &y_all_condensed,
                &y_bar_all_condensed,
                &z_body_condensed,
            )
        };
        // debug!("Packing took {} us", now.elapsed().as_micros());

        res.push(packed);
    }

    res
}


pub fn pack_many_lwes_inspir_without_rotations<'a>(
    packing_params: &'a PackParams,
    precomp_inspir_vec: &Vec<PrecompInsPIR<'a>>,
    b_values: &[u64],
    packing_keys: &PackingKeys<'a>,
    gamma: usize,
) -> Vec<PolyMatrixRaw<'a>> {
    // assert!(b_values.len() >= num_rlwe_outputs * packed_per_rlwe);
    let num_rlwe_outputs = (b_values.len() as f64 / gamma as f64).ceil() as usize;

    let mut res = Vec::with_capacity(num_rlwe_outputs);
    let group_size = packing_params.params.poly_len / gamma;
    for i in 0..num_rlwe_outputs {
        let which_group = i / group_size;
        let within_group = i % group_size;

        let mut b_poly = PolyMatrixRaw::zero(&packing_params.params, 1, 1);

        for j in 0..gamma {
            let index = which_group * packing_params.params.poly_len + within_group * gamma + j;
            if index >= b_values.len() {
                continue;
            }
            b_poly.get_poly_mut(0, 0)[j] = b_values[index];
        }

        // let now = Instant::now();
        let packed = if gamma <= packing_params.params.poly_len / 2 {
            let y_condensed = packing_keys.y_body_condensed.as_ref().unwrap();
            
            packing_with_preprocessing_online_without_rotations(
                &packing_params,
                &precomp_inspir_vec[i],
                &b_poly,
                &y_condensed,
            )
        } else {
            let y_condensed = packing_keys.y_body_condensed.as_ref().unwrap();
            let z_body_condensed = packing_keys.z_body_condensed.as_ref().unwrap();
            full_packing_with_preprocessing_online_without_rotations(
                &packing_params,
                &precomp_inspir_vec[i],
                &b_poly,
                &y_condensed,
                &z_body_condensed,
            )
        };
        // debug!("Packing took {} us", now.elapsed().as_micros());

        res.push(packed);
    }

    res
}

fn rotation_poly<'a>(params: &'a Params, amount: usize) -> PolyMatrixNTT<'a> {
    let mut res = PolyMatrixRaw::zero(params, 1, 1);
    res.data[amount] = 1;
    res.ntt()
}

pub fn pack_using_single_with_offset<'a>(
    params: &'a Params,
    pub_params: &[PolyMatrixNTT<'a>],
    cts: &[PolyMatrixNTT<'a>],
    offset: usize,
) -> PolyMatrixRaw<'a> {
    let mut res = PolyMatrixNTT::zero(params, 2, 1);
    for i in 0..cts.len() {
        let packed_single = pack_single_lwe(params, pub_params, &cts[i]);
        let rotation = rotation_poly(params, offset + i);
        let rotated = scalar_multiply_alloc(&rotation, &packed_single);
        add_into(&mut res, &rotated);
    }
    let res_raw = res.raw();
    return res_raw;
}

fn swap_midpoint<T>(a: &mut [T]) {
    let len = a.len();
    let (a, b) = a.split_at_mut(len / 2);
    a.swap_with_slice(b);
}

pub fn produce_table(poly_len: usize, chunk_size: usize) -> Vec<usize> {
    let mut cur = (0..poly_len).collect::<Vec<_>>();

    let outer_chunk_size = poly_len / (chunk_size / 2);
    println!("outer_chunk_size {}", outer_chunk_size);

    let mut do_it = true;
    for outer_chunk in cur.chunks_mut(outer_chunk_size) {
        if !do_it {
            do_it = true;
            continue;
        }
        do_it = false;

        for chunk in outer_chunk.chunks_mut(chunk_size) {
            let mut offs = 0;
            let mut to_add_to_offs = (chunk_size / 2).min(chunk.len() / 2); // weird hack
            while to_add_to_offs > 0 {
                swap_midpoint(&mut chunk[offs..]);
                offs += to_add_to_offs;
                to_add_to_offs /= 2;
            }
        }
    }

    cur
}

pub fn automorph_ntt_tables(poly_len: usize, log2_poly_len: usize) -> Vec<Vec<usize>> {
    let mut tables = Vec::new();
    for i in 0..log2_poly_len {
        let chunk_size = 1 << i;
        println!("table {}", i);
        let table = produce_table(poly_len, 2 * chunk_size);
        println!("table {:?}", &table.as_slice()[..32]);
        tables.push(table);
    }

    tables
}

pub fn bit_reversal(x: usize, log_len: usize) -> usize {
    let mut reversed = 0;
    let mut i = 0;
    let mut j = log_len - 1;

    while i < j {
        let bit_i = (x >> i) & 1;
        let bit_j = (x >> j) & 1;

        reversed |= (bit_i << j) | (bit_j << i);
        i += 1;
        j -= 1;
    }

    // If the length is odd, the middle bit remains unchanged
    if log_len % 2 == 1 {
        reversed |= (x >> (log_len / 2)) << (log_len / 2);
    }

    reversed
}
pub fn generate_automorph_tables_brute_force(params: &Params) -> Vec<Vec<usize>> {
    let mut tables = Vec::new();
    // for i in (1..=params.poly_len_log2) {
    for t in (1..2 * params.poly_len).step_by(2) {
        let mut table_candidate = vec![0usize; params.poly_len];

        // for 2048 balls and 2^28 bins, we will have a collision ~1% of the time
        // so, we redo if necessary
        loop {
            // let t = (1 << i) + 1;

            let poly = PolyMatrixRaw::random(&params, 1, 1);
            let poly_ntt = poly.ntt();

            let poly_auto = automorph_alloc(&poly, t);
            let poly_auto_ntt = poly_auto.ntt();

            let pol_orig = (&poly_ntt.get_poly(0, 0)[..params.poly_len]).to_vec();
            let pol_auto = (&poly_auto_ntt.get_poly(0, 0)[..params.poly_len]).to_vec();

            let mut must_redo = false;

            for i in 0..params.poly_len {
                let mut total = 0;
                let mut found = None;
                for j in 0..params.poly_len {
                    if pol_orig[i] == pol_auto[j] {
                        total += 1;
                        found = Some(j);
                    }
                }
                table_candidate[found.unwrap()] = i;
                if total != 1 {
                    must_redo = true;
                    break;
                }
            }

            if !must_redo {
                break;
            }
        }
        tables.push(table_candidate);
    }
    tables
}

// pub fn generate_automorph_tables_brute_force(params: &Params) -> Vec<Vec<usize>> {
//     let mut tables = Vec::new();
//     let n = params.poly_len as usize;
//     let m = 2 * n;

//     // Iterate through the same odd automorphism factors
//     for t in (1..m).step_by(2) {
//         let mut table_candidate = vec![0usize; n as usize];

//         // We need the modular inverse of t to reverse the permutation
//         let t_inv = invert_uint_mod(t as u64, m as u64).unwrap() as usize;

//         // For each new index 'j', find the original index 'i' it came from
//         for j in 0..n {
//             let j_exp = 2 * j + 1;

//             // Apply the inverse permutation on the exponent
//             let i_exp = (t_inv * j_exp) % m;
            
//             // Convert the old exponent back to an index
//             let i = (i_exp - 1) / 2;

//             // Store the mapping: table[new_index] = old_index
//             table_candidate[j as usize] = i as usize;
//         }
//         tables.push(table_candidate);
//     }
//     tables
// }


pub fn apply_automorph_ntt_raw<'a>(
    params: &Params,
    poly: &[u64],
    out: &mut [u64],
    t: usize,
    tables: &[Vec<usize>],
) {
    let poly_len = params.poly_len;
    // table_idx = log2(poly_len / (t - 1))
    // let table_idx = (poly_len / (t - 1)).trailing_zeros() as usize;
    let table_idx = (t - 1) / 2;
    let table = &tables[table_idx];

    for i in 0..poly_len {
        out[i] += poly[table[i]];
    }
}

pub fn apply_automorph_ntt<'a>(
    params: &'a Params,
    tables: &[Vec<usize>],
    mat: &PolyMatrixNTT<'a>,
    res: &mut PolyMatrixNTT<'a>,
    t: usize,
) {
    // run apply_automorph_ntt on each poly in the matrix
    // let mut res = PolyMatrixNTT::zero(params, mat.rows, mat.cols);
    for i in 0..mat.rows {
        for j in 0..mat.cols {
            let poly = mat.get_poly(i, j);
            let mut res_poly: Vec<&mut [u64]> = res.get_poly_mut(i, j).chunks_exact_mut(params.poly_len).collect();
            // for (chunk, res_chunk) in
                // poly.chunks_exact(params.poly_len).zip(res_poly.chunks_exact_mut(params.poly_len))
            for (index, chunk) in poly.chunks_exact(params.poly_len).enumerate()
            {
                apply_automorph_ntt_raw(params, chunk, res_poly[index], t, tables);
            }
        }
    }
    // res
}

pub fn apply_automorph_ntt_double<'a>(
    params: &'a Params,
    tables: &[Vec<usize>],
    mat: &PolyMatrixNTT<'a>,
    res_1: &mut PolyMatrixNTT<'a>,
    res_2: &mut PolyMatrixNTT<'a>,
    t: usize,
) {
    // run apply_automorph_ntt on each poly in the matrix
    // let mut res = PolyMatrixNTT::zero(params, mat.rows, mat.cols);
    for i in 0..mat.rows {
        for j in 0..mat.cols {
            let poly = mat.get_poly(i, j);
            let mut res_1_poly: Vec<&mut [u64]> = res_1.get_poly_mut(i, j).chunks_exact_mut(params.poly_len).collect();
            let mut res_2_poly: Vec<&mut [u64]> = res_2.get_poly_mut(i, j).chunks_exact_mut(params.poly_len).collect();
            // for (chunk, res_1_chunk, res_2_chunk) in
            //     izip!(poly.chunks_exact(params.poly_len), res_1_poly.chunks_exact_mut(params.poly_len), res_2_poly.chunks_exact_mut(params.poly_len))
            for (index, chunk) in poly.chunks_exact(params.poly_len).enumerate()
            {
                apply_automorph_ntt_raw(params, chunk, res_1_poly[index], t, tables);
                apply_automorph_ntt_raw(params, chunk, res_2_poly[index], 2*params.poly_len - t, tables);
            }
        }
    }
    // res
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PackingType {
    NoPacking,
    CDKS,
    InspiRING,
}

impl FromStr for PackingType {
    type Err = String;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        match input.to_lowercase().as_str() {
            "nopacking" => Ok(PackingType::NoPacking),
            "cdks" => Ok(PackingType::CDKS),
            "inspiring" => Ok(PackingType::InspiRING),
            _ => Err(format!("Invalid input: {}", input)),
        }
    }
}
/// Trait to convert an enum to a string representation
pub trait ToStr {
    fn to_str(&self) -> String;
}

impl ToStr for PackingType {
    fn to_str(&self) -> String {
        match self {
            PackingType::NoPacking => "NoPacking",
            PackingType::CDKS => "CDKS",
            PackingType::InspiRING => "InspiRING",
        }
        .to_string()
    }
}

impl Default for PackingType {
    fn default() -> Self {
        PackingType::CDKS // Specify the default variant
    }
}

#[derive(Clone)]
pub struct PackParams<'a> {
    pub params: &'a Params,
    pub num_to_pack: usize,
    pub tables: Vec<Vec<usize>>,
    pub gen_pows: Vec<usize>,
    pub mod_inv_poly: PolyMatrixNTT<'a>,
    pub monomial_ntts: Vec<PolyMatrixNTT<'a>>,
    pub neg_monomial_ntts: Vec<PolyMatrixNTT<'a>>,
}

impl PackParams<'_> {
    pub fn new<'a>(params: &'a Params, num_to_pack: usize) -> PackParams<'a> {
        debug!("Starting tables");
        let tables = generate_automorph_tables_brute_force(params);
        debug!("Got tables");
        let gen: usize =
            if num_to_pack < params.poly_len { (2 * params.poly_len / num_to_pack) + 1 } else { 5 };

        let mut gen_pows = Vec::new();
        for i in 0..params.poly_len {
            gen_pows
                .push(exponentiate_uint_mod(gen as u64, i as u64, 2 * params.poly_len as u64)
                    as usize);
        }
        let mod_inv = invert_uint_mod(num_to_pack as u64, params.modulus).unwrap();
        let mod_inv_poly = single_poly(params, mod_inv).ntt();

        let mut monomial_ntts = Vec::new();
        let mut neg_monomial_ntts = Vec::new();
        for j in 0..params.poly_len {
            let mut monomial = PolyMatrixRaw::zero(params, 1, 1);
            monomial.get_poly_mut(0, 0)[j] = 1;
            let mono_ntt = monomial.ntt();
            monomial_ntts.push(mono_ntt.clone());
            neg_monomial_ntts.push(-&mono_ntt);
        }
        PackParams {
            params,
            num_to_pack,
            tables,
            gen_pows,
            mod_inv_poly,
            monomial_ntts,
            neg_monomial_ntts,
        }
    }

    // An incomplete version of the packing parameters, used by the client
    pub fn new_fast<'a>(params: &'a Params, num_to_pack: usize) -> PackParams<'a> {
        let gen: usize =
            if num_to_pack < params.poly_len { (2 * params.poly_len / num_to_pack) + 1 } else { 5 };

        let mut gen_pows = Vec::new();
        for i in 0..params.poly_len {
            gen_pows
                .push(exponentiate_uint_mod(gen as u64, i as u64, 2 * params.poly_len as u64)
                    as usize);
        }

        PackParams {
            params,
            num_to_pack,
            tables: vec![],
            gen_pows,
            mod_inv_poly:PolyMatrixNTT::zero(&params, 1, 1),
            monomial_ntts: vec![],
            neg_monomial_ntts: vec![],
        }
    }
    
    pub fn empty<'a>(params: &'a Params) -> PackParams<'_> {
        PackParams {
            params,
            num_to_pack: 0,
            tables: vec![],
            gen_pows: vec![],
            mod_inv_poly: PolyMatrixNTT::zero(&params, 1, 1),
            monomial_ntts: vec![],
            neg_monomial_ntts: vec![],
        }
    }
}

#[derive(Clone)]
pub struct PrecompInsPIR<'a> {
    pub a_hat: PolyMatrixRaw<'a>,
    pub bold_t_condensed: PolyMatrixNTT<'a>,
    pub bold_t_bar_condensed: PolyMatrixNTT<'a>,
    pub bold_t_hat_condensed: PolyMatrixNTT<'a>,
}

pub fn generate_rotations<'a>(
    packing_params: &PackParams<'a>,
    to_rotate: &PolyMatrixNTT<'a>,
) -> PolyMatrixNTT<'a> {
    let params = packing_params.params;
    let num_to_pack = packing_params.num_to_pack;
    let tables = &packing_params.tables;
    let gen_pows = &packing_params.gen_pows;

    let mut rotations_all = PolyMatrixNTT::zero(&params, num_to_pack - 1, params.t_exp_left);
    for i in 0..num_to_pack - 1 {
        let mut rotated_w = PolyMatrixNTT::zero(&params, 1, params.t_exp_left);
        apply_automorph_ntt(&params, &tables, &to_rotate, &mut rotated_w, gen_pows[i]);
        rotations_all.copy_into(&rotated_w, i, 0);
    }
    rotations_all
}

pub fn generate_rotations_double<'a>(
    packing_params: &PackParams<'a>,
    to_rotate: &PolyMatrixNTT<'a>,
) -> (PolyMatrixNTT<'a>, PolyMatrixNTT<'a>) {
    let params = packing_params.params;
    let num_to_pack = packing_params.num_to_pack;
    let tables = &packing_params.tables;
    let gen_pows = &packing_params.gen_pows;
    assert_eq!(num_to_pack, params.poly_len);

    let num_rotations = params.poly_len / 2 - 1;

    let mut rotations_all = PolyMatrixNTT::zero(&params, num_rotations, params.t_exp_left);
    let mut rotations_bar_all = PolyMatrixNTT::zero(&params, num_rotations, params.t_exp_left);
    for i in 0..num_rotations {
        let mut rotated_w_1 = PolyMatrixNTT::zero(&params, 1, params.t_exp_left);
        let mut rotated_w_2 = PolyMatrixNTT::zero(&params, 1, params.t_exp_left);
        apply_automorph_ntt_double(&params, &tables, &to_rotate, &mut rotated_w_1, &mut rotated_w_2, gen_pows[i]);
        rotations_all.copy_into(&rotated_w_1, i, 0);
        rotations_bar_all.copy_into(&rotated_w_2, i, 0);
    }
    (rotations_all, rotations_bar_all)
}

#[derive(Clone)]
pub struct OfflinePackingKeys<'a> {

    pub packing_params: Option<&'a PackParams<'a>>,
    pub full_key: bool,

    pub w_seed: [u8; 32],
    pub v_seed: [u8; 32],

    pub w_mask: Option<PolyMatrixNTT<'a>>,
    pub v_mask: Option<PolyMatrixNTT<'a>>,

    pub w_all: Option<PolyMatrixNTT<'a>>,
    pub w_bar_all: Option<PolyMatrixNTT<'a>>,

}

impl OfflinePackingKeys<'_> {

    pub fn init_empty<'a>() -> OfflinePackingKeys<'a> {
        OfflinePackingKeys {
            packing_params: None,
            full_key: false,
            w_seed: [0; 32],
            v_seed: [0; 32],
            w_mask: None,
            w_all: None,
            w_bar_all: None,
            v_mask: None,
        }
    }

    pub fn init<'a>(packing_params: &'a PackParams<'a>, w_seed: [u8; 32]) -> OfflinePackingKeys<'_> {
        let w_mask = PolyMatrixNTT::random_rng(
            &packing_params.params,
            1,
            packing_params.params.t_exp_left,
            &mut ChaCha20Rng::from_seed(w_seed),
        );
        // let w_mask = PolyMatrixNTT::zero(&packing_params.params, 1, packing_params.params.t_exp_left);
        let x = generate_rotations(&packing_params, &w_mask);
        let w_all = Some(x);
        OfflinePackingKeys{
            packing_params: Some(packing_params),
            full_key: false,
            w_seed,
            v_seed: [0; 32],
            w_mask: Some(w_mask),
            w_all: w_all,
            w_bar_all: None,
            v_mask: None,
        }
    }

    pub fn init_full<'a>(packing_params: &'a PackParams<'a>, w_seed: [u8; 32], v_seed: [u8; 32]) -> OfflinePackingKeys<'a> {
        let w_mask = PolyMatrixNTT::random_rng(
            &packing_params.params,
            1,
            packing_params.params.t_exp_left,
            &mut ChaCha20Rng::from_seed(w_seed),
        );
        let v_mask = PolyMatrixNTT::random_rng(
            &packing_params.params,
            1,
            packing_params.params.t_exp_left,
            &mut ChaCha20Rng::from_seed(v_seed),
        );

        let (x,y) = generate_rotations_double(&packing_params, &w_mask);
        let (w_all, w_bar_all) = (Some(x), Some(y));

        OfflinePackingKeys{
            packing_params: Some(packing_params),
            full_key: true,
            w_seed,
            v_seed,
            w_mask: Some(w_mask),
            w_all: w_all,
            w_bar_all: w_bar_all,
            v_mask: Some(v_mask),
        }
    }

}

#[derive(Clone)]
pub struct PackingKeys<'a> {

    pub packing_type: PackingType,

    // Inspiring Stuff
    pub full_key: bool,
    pub packing_params: Option<PackParams<'a>>,
    pub y_body: Option<PolyMatrixNTT<'a>>,
    pub z_body: Option<PolyMatrixNTT<'a>>,
    pub y_body_condensed: Option<PolyMatrixNTT<'a>>,
    pub z_body_condensed: Option<PolyMatrixNTT<'a>>,    

    pub expanded: bool,
    pub y_all_condensed: Option<PolyMatrixNTT<'a>>,
    pub y_bar_all_condensed: Option<PolyMatrixNTT<'a>>,    

    // CDKS stuff
    pub params: Option<Params>,
    pub pack_pub_params_row_1s: Vec<PolyMatrixNTT<'a>>,    
    pub fake_pack_pub_params: Vec<PolyMatrixNTT<'a>>,
}

impl PackingKeys<'_> {

    pub fn init_full<'a>(
        packing_params: &PackParams<'a>,
        sk_reg: &PolyMatrixRaw<'a>, 
        w_seed :[u8; 32],
        v_seed :[u8; 32],
    ) -> PackingKeys<'a> {
        let w_mask: PolyMatrixNTT<'_> = PolyMatrixNTT::random_rng(
            &packing_params.params,
            1,
            packing_params.params.t_exp_left,
            &mut ChaCha20Rng::from_seed(w_seed),
        );
        let v_mask: PolyMatrixNTT<'_> = PolyMatrixNTT::random_rng(
            &packing_params.params,
            1,
            packing_params.params.t_exp_left,
            &mut ChaCha20Rng::from_seed(v_seed),
        );                
        let y_body = generate_ksk_body(
            &packing_params.params,
            &sk_reg,
            packing_params.gen_pows[1],
            &w_mask,
            &mut ChaCha20Rng::from_entropy(),
        );
        let z_body = generate_ksk_body(
            &packing_params.params,
            &sk_reg,
            2 * packing_params.params.poly_len - 1,
            &v_mask,
            &mut ChaCha20Rng::from_entropy(),
        );
        let y_body_condensed = condense_matrix(&packing_params.params, &y_body);
        let z_body_condensed = condense_matrix(&packing_params.params, &z_body);
        PackingKeys {
            packing_type: PackingType::InspiRING,
            packing_params: Some(packing_params.clone()),
            full_key: true,
            y_body: Some(y_body),
            z_body: Some(z_body),
            y_body_condensed: Some(y_body_condensed),
            z_body_condensed: Some(z_body_condensed),
            expanded: false,
            y_all_condensed: None,
            y_bar_all_condensed: None,
            params: None,
            pack_pub_params_row_1s: vec![],
            fake_pack_pub_params: vec![],
        }
    }

    pub fn init<'a> (
        packing_params: &PackParams<'a>,
        sk_reg: &PolyMatrixRaw<'a>, 
        w_seed :[u8; 32],
    ) -> PackingKeys<'a> {
        let w_mask: PolyMatrixNTT<'_> = PolyMatrixNTT::random_rng(
            &packing_params.params,
            1,
            packing_params.params.t_exp_left,
            &mut ChaCha20Rng::from_seed(w_seed),
        );
        // let w_mask = PolyMatrixNTT::zero(&packing_params.params, 1, packing_params.params.t_exp_left);
        let y_body = generate_ksk_body(
            &packing_params.params,
            &sk_reg,
            packing_params.gen_pows[1],
            &w_mask,
            &mut ChaCha20Rng::from_entropy(),
        );
        let y_body_condensed = condense_matrix(&packing_params.params, &y_body);
        PackingKeys {
            packing_type: PackingType::InspiRING,
            packing_params: Some(packing_params.clone()),
            full_key: false,
            y_body : Some(y_body),
            z_body : None,
            y_body_condensed : Some(y_body_condensed),
            z_body_condensed : None,
            expanded: false,
            y_all_condensed: None,
            y_bar_all_condensed: None,
            params: None,
            pack_pub_params_row_1s: vec![],
            fake_pack_pub_params: vec![],
        }        
    }

    pub fn init_cdks<'a>(
        params: &'a Params,
        sk_reg: &PolyMatrixRaw<'a>,
        static_seed_2: [u8; 32],  
    ) -> PackingKeys<'a> {
        let pack_pub_params = raw_generate_expansion_params(
            &params,
            &sk_reg,
            params.poly_len_log2,
            params.t_exp_left,
            &mut ChaCha20Rng::from_entropy(),
            &mut ChaCha20Rng::from_seed(static_seed_2),
        );
        // let pub_params_size = get_vec_pm_size_bytes(&pack_pub_params) / 2;
        let mut pack_pub_params_row_1s = pack_pub_params.to_vec();
        for i in 0..pack_pub_params.len() {
            pack_pub_params_row_1s[i] =
                pack_pub_params[i].submatrix(1, 0, 1, pack_pub_params[i].cols);
            pack_pub_params_row_1s[i] =
                condense_matrix(&params, &pack_pub_params_row_1s[i]);
        }

        let mut fake_pack_pub_params = pack_pub_params.clone();
        // zero out all of the second rows
        for i in 0..pack_pub_params.len() {
            for col in 0..pack_pub_params[i].cols {
                fake_pack_pub_params[i].get_poly_mut(1, col).fill(0);
            }
        }

        PackingKeys {
            packing_type: PackingType::CDKS,
            packing_params: None,
            full_key: false,
            y_body : None,
            z_body : None,
            y_body_condensed : None,
            z_body_condensed : None,
            expanded: false,
            y_all_condensed: None,
            y_bar_all_condensed: None,
            params: Some(params.clone()),
            pack_pub_params_row_1s : pack_pub_params_row_1s,
            fake_pack_pub_params: fake_pack_pub_params,
        }
    }

    pub fn get_gamma(&self) -> usize {
        if self.packing_type == PackingType::InspiRING {
            self.packing_params.as_ref().unwrap().num_to_pack
        } else {
            self.params.as_ref().unwrap().poly_len
        }
    }

    pub fn get_size_bytes(&self) -> usize {
        if self.packing_type == PackingType::InspiRING {
            if self.full_key {
                get_vec_pm_size_bytes(&vec![self.y_body.as_ref().unwrap().clone()])
                + get_vec_pm_size_bytes(&vec![self.z_body.as_ref().unwrap().clone()])
            } else {
                get_vec_pm_size_bytes(&vec![self.y_body.as_ref().unwrap().clone()])
            }
        } else {
            get_vec_pm_size_bytes(&self.pack_pub_params_row_1s)
        }
    }

    pub fn expand<'a>(&mut self) {
        assert_eq!(self.packing_type, PackingType::InspiRING);
        if !self.expanded {
            let packing_params = self.packing_params.as_ref().unwrap();
            if self.full_key {
                let (x,y) = generate_rotations_double(&packing_params, &self.y_body_condensed.as_ref().unwrap());
                (self.y_all_condensed, self.y_bar_all_condensed) = (Some(x), Some(y));
            } else {
                let x = generate_rotations(&packing_params, &self.y_body_condensed.as_ref().unwrap());
                self.y_all_condensed = Some(x);
            }
            self.expanded = true;
        }
    }

}

pub fn packing_with_preprocessing_offline<'a>(
    packing_params: &PackParams<'a>,
    w_all: &PolyMatrixNTT<'a>,
    a_ct_tilde: &Vec<PolyMatrixNTT<'a>>
) -> PrecompInsPIR<'a> {
    let params = packing_params.params;
    let num_to_pack = packing_params.num_to_pack;
    let tables = &packing_params.tables;
    let gen_pows = &packing_params.gen_pows;
    let monomial_ntts = &packing_params.monomial_ntts;
    let neg_monomial_ntts = &packing_params.neg_monomial_ntts;
    let mod_inv_poly = &packing_params.mod_inv_poly;

    let now = Instant::now();

    let non_zeros = a_ct_tilde.len();

    let mut r_all = Vec::with_capacity(num_to_pack);
    for i in 0..num_to_pack {
        let mut r_pow_i = PolyMatrixNTT::zero(&params, 1, 1);
        for j in 0..non_zeros {
            let res_poly = r_pow_i.get_poly_mut(0, 0);
            let index = (j * packing_params.gen_pows[(params.poly_len - i) % params.poly_len])
                % (2 * params.poly_len);
            let pol1 = if index < params.poly_len {
                monomial_ntts[index % params.poly_len].get_poly(0, 0)
            } else {
                neg_monomial_ntts[index % params.poly_len].get_poly(0, 0)
            };
            let pol2 = a_ct_tilde[j].get_poly(0, 0);
            multiply_add_poly_avx(&params, res_poly, pol1, pol2);
            let reduction_steps = 1 << (64 - 2 * params.q2_bits);
            if (j + 1) % reduction_steps == 0 {
                fast_reduce(&mut r_pow_i);
            }
        }
        fast_reduce(&mut r_pow_i);

        let mut r_pow_i_mul = PolyMatrixNTT::zero(&params, 1, 1);
        scalar_multiply(&mut r_pow_i_mul, &mod_inv_poly, &r_pow_i);

        let mut r_pow_i_mul_rotated = PolyMatrixNTT::zero(&params, 1, 1);
        apply_automorph_ntt(&params, &tables, &r_pow_i_mul, &mut r_pow_i_mul_rotated, gen_pows[i]);
        r_all.push(r_pow_i_mul_rotated);
    }
    debug!("Inner product for R took {} us", now.elapsed().as_micros());

    let now = Instant::now();

    let mut bold_t = PolyMatrixNTT::zero(&params, num_to_pack - 1, params.t_exp_left);

    for i in (0..num_to_pack - 1).rev() {
        let gadget = gadget_invert_transposed_alloc(&r_all[i + 1].raw(), params.t_exp_left).ntt();
        bold_t.copy_into(&gadget, i, 0);

        let mut res = PolyMatrixNTT::zero(&params, 1, 1);
        let poly_res = res.get_poly_mut(0, 0);
        for k in 0..params.t_exp_left {
            let poly1 = w_all.get_poly(i, k);
            let poly2 = bold_t.get_poly(i, k);
            multiply_add_poly_avx(params, poly_res, poly1, poly2);
        }
        fast_add_into(&mut r_all[i], &res);
    }
    debug!("Gadget inversion for R took {} us", now.elapsed().as_micros());

    PrecompInsPIR {
        a_hat: r_all[0].raw(),
        bold_t_condensed: condense_matrix(&params, &bold_t),
        bold_t_bar_condensed: PolyMatrixNTT::zero(&params, 1, 1),
        bold_t_hat_condensed: PolyMatrixNTT::zero(&params, 1, 1),
    }
}


pub fn packing_with_preprocessing_offline_without_rotations<'a>(
    packing_params: &PackParams<'a>,
    w_all: &PolyMatrixNTT<'a>,
    a_ct_tilde: &Vec<PolyMatrixNTT<'a>>
) -> PrecompInsPIR<'a> {
    let params = packing_params.params;
    let num_to_pack = packing_params.num_to_pack;
    let tables = &packing_params.tables;
    let gen_pows = &packing_params.gen_pows;
    let monomial_ntts = &packing_params.monomial_ntts;
    let neg_monomial_ntts = &packing_params.neg_monomial_ntts;
    let mod_inv_poly = &packing_params.mod_inv_poly;

    let now = Instant::now();

    let non_zeros = a_ct_tilde.len();

    let mut r_all = Vec::with_capacity(num_to_pack);
    for i in 0..num_to_pack {
        let mut r_pow_i = PolyMatrixNTT::zero(&params, 1, 1);
        for j in 0..non_zeros {
            let res_poly = r_pow_i.get_poly_mut(0, 0);
            let index = (j * packing_params.gen_pows[(params.poly_len - i) % params.poly_len])
                % (2 * params.poly_len);
            let pol1 = if index < params.poly_len {
                monomial_ntts[index % params.poly_len].get_poly(0, 0)
            } else {
                neg_monomial_ntts[index % params.poly_len].get_poly(0, 0)
            };
            let pol2 = a_ct_tilde[j].get_poly(0, 0);
            multiply_add_poly_avx(&params, res_poly, pol1, pol2);
            let reduction_steps = 1 << (64 - 2 * params.q2_bits);
            if (j + 1) % reduction_steps == 0 {
                fast_reduce(&mut r_pow_i);
            }
        }
        fast_reduce(&mut r_pow_i);

        let mut r_pow_i_mul = PolyMatrixNTT::zero(&params, 1, 1);
        scalar_multiply(&mut r_pow_i_mul, &mod_inv_poly, &r_pow_i);

        let mut r_pow_i_mul_rotated = PolyMatrixNTT::zero(&params, 1, 1);
        apply_automorph_ntt(&params, &tables, &r_pow_i_mul, &mut r_pow_i_mul_rotated, gen_pows[i]);
        r_all.push(r_pow_i_mul_rotated.raw());
    }
    debug!("Inner product for R took {} us", now.elapsed().as_micros());

    let now = Instant::now();

    let mut bold_t = PolyMatrixNTT::zero(&params, num_to_pack - 1, params.t_exp_left);
    let mut bold_t_prime = PolyMatrixNTT::zero(&params, num_to_pack - 1, params.t_exp_left);

    for i in (0..num_to_pack - 1).rev() {
        let gadget = gadget_invert_transposed_alloc(&r_all[i + 1], params.t_exp_left).ntt();
        bold_t.copy_into(&gadget, i, 0);

        let mut rotated_gadget = PolyMatrixNTT::zero(&params, 1, params.t_exp_left);
        apply_automorph_ntt(&params, &packing_params.tables, &gadget, &mut rotated_gadget, gen_pows[(params.poly_len - i) % params.poly_len]);        
        bold_t_prime.copy_into(&rotated_gadget, i, 0);

        let mut res = PolyMatrixNTT::zero(&params, 1, 1);
        let poly_res = res.get_poly_mut(0, 0);
        for k in 0..params.t_exp_left {
            let poly1 = w_all.get_poly(i, k);
            let poly2 = bold_t.get_poly(i, k);
            multiply_add_poly_avx(params, poly_res, poly1, poly2);
        }
        fast_reduce(&mut res);
        r_all[i] = &r_all[i] + &res.raw();
    }
    debug!("Gadget inversion for R took {} us", now.elapsed().as_micros());


    PrecompInsPIR {
        a_hat: r_all[0].clone(),
        bold_t_condensed: condense_matrix(&params, &bold_t_prime),
        bold_t_bar_condensed: PolyMatrixNTT::zero(&params, 1, 1),
        bold_t_hat_condensed: PolyMatrixNTT::zero(&params, 1, 1),
    }
}


pub fn packing_with_preprocessing_online<'a>(
    packing_params: &PackParams<'a>,
    precomp_inspiring: &PrecompInsPIR<'a>,
    b_poly: &PolyMatrixRaw<'a>,
    y_all_condensed: &PolyMatrixNTT<'a>,
) -> PolyMatrixRaw<'a> {
    let params = packing_params.params;
    let num_to_pack = packing_params.num_to_pack;

    let a_hat = &precomp_inspiring.a_hat;
    let bold_t_condensed = &precomp_inspiring.bold_t_condensed;

    let mut sum_poly = PolyMatrixNTT::zero(&params, 1, 1);

    let addition_capacity = 1 << (64 - 2 * params.q2_bits - 1);
    let mut num_added = 0;

    let now = Instant::now();
    let mut temp_poly = PolyMatrixNTT::zero(&params, 1, 1);
    for i in 0..(num_to_pack - 1) {

        fast_multiply_no_reduce_in_range_generic(
            &params,
            &mut temp_poly,
            &y_all_condensed,
            &bold_t_condensed,
            i*params.t_exp_left,
            i*params.t_exp_left,
            params.t_exp_left,
        );
        fast_add_into_no_reduce(&mut sum_poly, &temp_poly);
        num_added += params.t_exp_left;
        if (num_added >= addition_capacity) || (i==num_to_pack-2) {
            fast_reduce(&mut sum_poly);
            num_added = 0;
        }
    }

    debug!("Multiplication and addition for R took {} us", now.elapsed().as_micros());

    let now = Instant::now();
    let mut final_b_poly = PolyMatrixRaw::zero(&params, 1, 1);
    add_raw(&mut final_b_poly, &b_poly, &sum_poly.raw());
    debug!("Addition for B took {} us", now.elapsed().as_micros());

    let mut packed_raw = PolyMatrixRaw::zero(&params, 2, 1);
    packed_raw.copy_into(a_hat, 0, 0);
    packed_raw.copy_into(&final_b_poly, 1, 0);

    debug!("Packing online took {} us", now.elapsed().as_micros());

    packed_raw
}

pub fn packing_with_preprocessing_online_without_rotations<'a>(
    packing_params: &PackParams<'a>,
    precomp_inspiring: &PrecompInsPIR<'a>,
    b_poly: &PolyMatrixRaw<'a>,
    y_condensed: &PolyMatrixNTT<'a>,
) -> PolyMatrixRaw<'a> {
    let params = packing_params.params;
    let num_to_pack = packing_params.num_to_pack;

    let a_hat = &precomp_inspiring.a_hat;
    let bold_t_condensed = &precomp_inspiring.bold_t_condensed;

    let mut sum_poly = PolyMatrixNTT::zero(&params, 1, 1);

    let addition_capacity = 1 << (64 - 2 * params.q2_bits - 1);
    let mut num_added = 0;

    let now = Instant::now();
    for i in 0..(num_to_pack - 1) {

        let mut temp_poly = PolyMatrixNTT::zero(&params, 1, 1);
        fast_multiply_no_reduce_in_range_generic(
            &params,
            &mut temp_poly,
            &y_condensed,
            &bold_t_condensed,
	        0,
	        i*params.t_exp_left,
	        params.t_exp_left,
        );
        
        let mut rotated_temp = PolyMatrixNTT::zero(&params, 1, 1);
        apply_automorph_ntt(&params, &packing_params.tables, &temp_poly, &mut rotated_temp, packing_params.gen_pows[i]);        
        
        fast_add_into_no_reduce(&mut sum_poly, &rotated_temp);
        num_added += params.t_exp_left;
        if (num_added >= addition_capacity) || (i==num_to_pack-2) {
            fast_reduce(&mut sum_poly);
            num_added = 0;
        }        
    }

    debug!("Multiplication and addition for R took {} us", now.elapsed().as_micros());

    let now = Instant::now();
    fast_reduce(&mut sum_poly);
    debug!("Reduction for R took {} us", now.elapsed().as_micros());

    let now = Instant::now();
    let mut final_b_poly = PolyMatrixRaw::zero(&params, 1, 1);
    add_raw(&mut final_b_poly, &b_poly, &sum_poly.raw());
    debug!("Addition for B took {} us", now.elapsed().as_micros());

    let mut packed_raw = PolyMatrixRaw::zero(&params, 2, 1);
    packed_raw.copy_into(a_hat, 0, 0);
    packed_raw.copy_into(&final_b_poly, 1, 0);

    debug!("Packing online took {} us", now.elapsed().as_micros());

    packed_raw
}


pub fn full_packing_with_preprocessing_offline<'a>(
    packing_params: &PackParams<'a>,
    w_all: &PolyMatrixNTT<'a>,
    w_bar_all: &PolyMatrixNTT<'a>,
    v_mask: &PolyMatrixNTT<'a>,
    a_ct_tilde: &Vec<PolyMatrixNTT<'a>>,
) -> PrecompInsPIR<'a> {
    let params = packing_params.params;
    let tables = &packing_params.tables;
    let gen_pows = &packing_params.gen_pows;
    let monomial_ntts = &packing_params.monomial_ntts;
    let neg_monomial_ntts = &packing_params.neg_monomial_ntts;
    let mod_inv_poly = &packing_params.mod_inv_poly;

    let num_to_pack = packing_params.num_to_pack;
    assert_eq!(num_to_pack, params.poly_len);
    let num_to_pack_half = num_to_pack >> 1;

    let non_zeros = a_ct_tilde.len();

    let mut r_all = Vec::with_capacity(num_to_pack_half);
    let mut r_bar_all = Vec::with_capacity(num_to_pack_half);
    for i in 0..num_to_pack_half {
        let mut r_pow_i = PolyMatrixNTT::zero(&params, 1, 1);
        let mut r_bar_pow_i = PolyMatrixNTT::zero(&params, 1, 1);
        for j in 0..non_zeros {
            let pol2 = a_ct_tilde[j].get_poly(0, 0);

            let res_poly = r_pow_i.get_poly_mut(0, 0);
            let index = (j * packing_params.gen_pows[(params.poly_len - i) % params.poly_len])
                % (2 * params.poly_len);
            let pol1 = if index < params.poly_len {
                monomial_ntts[index % params.poly_len].get_poly(0, 0)
            } else {
                neg_monomial_ntts[index % params.poly_len].get_poly(0, 0)
            };
            multiply_add_poly_avx(&params, res_poly, pol1, pol2);

            let res_poly = r_bar_pow_i.get_poly_mut(0, 0);
            let index = (2 * params.poly_len
                - (j * packing_params.gen_pows[(params.poly_len - i) % params.poly_len]))
                % (2 * params.poly_len);
            let pol1 = if index < params.poly_len {
                monomial_ntts[index % params.poly_len].get_poly(0, 0)
            } else {
                neg_monomial_ntts[index % params.poly_len].get_poly(0, 0)
            };
            multiply_add_poly_avx(&params, res_poly, pol1, pol2);
            let reduction_steps = 1 << (64 - 2 * params.q2_bits - 1);
            if (j + 1) % reduction_steps == 0 {
                fast_reduce(&mut r_pow_i);
                fast_reduce(&mut r_bar_pow_i);
            }
        }
        fast_reduce(&mut r_pow_i);
        fast_reduce(&mut r_bar_pow_i);

        let mut r_pow_i_mul = PolyMatrixNTT::zero(&params, 1, 1);
        let mut r_bar_pow_i_mul = PolyMatrixNTT::zero(&params, 1, 1);
        scalar_multiply(&mut r_pow_i_mul, &mod_inv_poly, &r_pow_i);
        scalar_multiply(&mut r_bar_pow_i_mul, &mod_inv_poly, &r_bar_pow_i);

        let mut r_pow_i_mul_rotated = PolyMatrixNTT::zero(&params, 1, 1);
        let mut r_bar_pow_i_mul_rotated = PolyMatrixNTT::zero(&params, 1, 1);
        apply_automorph_ntt(&params, &tables, &r_pow_i_mul, &mut r_pow_i_mul_rotated, gen_pows[i]);
        apply_automorph_ntt(
            &params,
            &tables,
            &r_bar_pow_i_mul,
            &mut r_bar_pow_i_mul_rotated,
            2 * params.poly_len - gen_pows[i],
        );

        r_all.push(r_pow_i_mul_rotated);
        r_bar_all.push(r_bar_pow_i_mul_rotated);
    }

    let mut bold_t_g = PolyMatrixNTT::zero(&params, num_to_pack_half - 1, params.t_exp_left);
    let mut bold_t_bar_g = PolyMatrixNTT::zero(&params, num_to_pack_half - 1, params.t_exp_left);

    for i in (0..num_to_pack_half - 1).rev() {
        let gadget = gadget_invert_transposed_alloc(&r_all[i + 1].raw(), params.t_exp_left).ntt();
        bold_t_g.copy_into(&gadget, i, 0);

        let mut res = PolyMatrixNTT::zero(&params, 1, 1);
        let poly_res = res.get_poly_mut(0, 0);
        for k in 0..params.t_exp_left {
            let poly1 = w_all.get_poly(i, k);
            let poly2 = bold_t_g.get_poly(i, k);
            multiply_add_poly_avx(params, poly_res, poly1, poly2);
        }
        fast_add_into(&mut r_all[i] ,&res);

        let gadget = gadget_invert_transposed_alloc(&r_bar_all[i + 1].raw(), params.t_exp_left).ntt();
        bold_t_bar_g.copy_into(&gadget, i, 0);

        let mut res = PolyMatrixNTT::zero(&params, 1, 1);
        let poly_res = res.get_poly_mut(0, 0);
        for k in 0..params.t_exp_left {
            let poly1 = w_bar_all.get_poly(i, k);
            let poly2 = bold_t_bar_g.get_poly(i, k);
            multiply_add_poly_avx(params, poly_res, poly1, poly2);
        }
        fast_add_into(&mut r_bar_all[i] ,&res);
    }

    let bold_t_hat = gadget_invert_transposed_alloc(&r_bar_all[0].raw(), params.t_exp_left).ntt();

    let mut res = PolyMatrixNTT::zero(&params, 1, 1);
    let poly_res = res.get_poly_mut(0, 0);
    for k in 0..params.t_exp_left {
        let poly1 = v_mask.get_poly(0, k);
        let poly2 = bold_t_hat.get_poly(0, k);
        multiply_add_poly_avx(params, poly_res, poly1, poly2);
    }
    fast_add_into(&mut r_all[0] ,&res);

    PrecompInsPIR {
        a_hat: r_all[0].raw(),
        bold_t_condensed: condense_matrix(&params, &bold_t_g),
        bold_t_bar_condensed: condense_matrix(&params, &bold_t_bar_g),
        bold_t_hat_condensed: condense_matrix(&params, &bold_t_hat),
    }
}

pub fn full_packing_with_preprocessing_offline_without_rotations<'a>(
    packing_params: &PackParams<'a>,
    w_all: &PolyMatrixNTT<'a>,
    w_bar_all: &PolyMatrixNTT<'a>,
    v_mask: &PolyMatrixNTT<'a>,
    a_ct_tilde: &Vec<PolyMatrixNTT<'a>>,
) -> PrecompInsPIR<'a> {
    let params = packing_params.params;
    let tables = &packing_params.tables;
    let gen_pows = &packing_params.gen_pows;
    let monomial_ntts = &packing_params.monomial_ntts;
    let neg_monomial_ntts = &packing_params.neg_monomial_ntts;
    let mod_inv_poly = &packing_params.mod_inv_poly;

    let num_to_pack = packing_params.num_to_pack;
    assert_eq!(num_to_pack, params.poly_len);
    let num_to_pack_half = num_to_pack >> 1;

    let non_zeros = a_ct_tilde.len();

    let mut r_all = Vec::new();
    let mut r_bar_all = Vec::new();
    for i in 0..num_to_pack_half {
        let mut r_pow_i = PolyMatrixNTT::zero(&params, 1, 1);
        let mut r_bar_pow_i = PolyMatrixNTT::zero(&params, 1, 1);
        for j in 0..non_zeros {
            let pol2 = a_ct_tilde[j].get_poly(0, 0);

            let res_poly = r_pow_i.get_poly_mut(0, 0);
            let index = (j * packing_params.gen_pows[(params.poly_len - i) % params.poly_len])
                % (2 * params.poly_len);
            let pol1 = if index < params.poly_len {
                monomial_ntts[index % params.poly_len].get_poly(0, 0)
            } else {
                neg_monomial_ntts[index % params.poly_len].get_poly(0, 0)
            };
            multiply_add_poly_avx(&params, res_poly, pol1, pol2);

            let res_poly = r_bar_pow_i.get_poly_mut(0, 0);
            let index = (2 * params.poly_len
                - (j * packing_params.gen_pows[(params.poly_len - i) % params.poly_len]))
                % (2 * params.poly_len);
            let pol1 = if index < params.poly_len {
                monomial_ntts[index % params.poly_len].get_poly(0, 0)
            } else {
                neg_monomial_ntts[index % params.poly_len].get_poly(0, 0)
            };
            multiply_add_poly_avx(&params, res_poly, pol1, pol2);
            let reduction_steps = 1 << (64 - 2 * params.q2_bits - 1);
            if (j + 1) % reduction_steps == 0 {
                fast_reduce(&mut r_pow_i);
                fast_reduce(&mut r_bar_pow_i);
            }
        }
        fast_reduce(&mut r_pow_i);
        fast_reduce(&mut r_bar_pow_i);

        let mut r_pow_i_mul = PolyMatrixNTT::zero(&params, 1, 1);
        let mut r_bar_pow_i_mul = PolyMatrixNTT::zero(&params, 1, 1);
        scalar_multiply(&mut r_pow_i_mul, &mod_inv_poly, &r_pow_i);
        scalar_multiply(&mut r_bar_pow_i_mul, &mod_inv_poly, &r_bar_pow_i);

        let mut r_pow_i_mul_rotated = PolyMatrixNTT::zero(&params, 1, 1);
        let mut r_bar_pow_i_mul_rotated = PolyMatrixNTT::zero(&params, 1, 1);
        apply_automorph_ntt(&params, &tables, &r_pow_i_mul, &mut r_pow_i_mul_rotated, gen_pows[i]);
        apply_automorph_ntt(
            &params,
            &tables,
            &r_bar_pow_i_mul,
            &mut r_bar_pow_i_mul_rotated,
            2 * params.poly_len - gen_pows[i],
        );

        r_all.push(r_pow_i_mul_rotated.raw());
        r_bar_all.push(r_bar_pow_i_mul_rotated.raw());
    }

    let mut bold_t_g = PolyMatrixNTT::zero(&params, num_to_pack_half - 1, params.t_exp_left);
    let mut bold_t_bar_g = PolyMatrixNTT::zero(&params, num_to_pack_half - 1, params.t_exp_left);
    let mut bold_t_g_prime = PolyMatrixNTT::zero(&params, num_to_pack_half - 1, params.t_exp_left);
    let mut bold_t_bar_g_prime = PolyMatrixNTT::zero(&params, num_to_pack_half - 1, params.t_exp_left);

    for i in (0..num_to_pack_half - 1).rev() {
        let gadget = gadget_invert_transposed_alloc(&r_all[i + 1], params.t_exp_left).ntt();
        bold_t_g.copy_into(&gadget, i, 0);

        let mut rotated_gadget = PolyMatrixNTT::zero(&params, 1, params.t_exp_left);
        apply_automorph_ntt(&params, &packing_params.tables, &gadget, &mut rotated_gadget, gen_pows[(params.poly_len - i) % params.poly_len]);        
        bold_t_g_prime.copy_into(&rotated_gadget, i, 0);

        let mut res = PolyMatrixNTT::zero(&params, 1, 1);
        let poly_res = res.get_poly_mut(0, 0);
        for k in 0..params.t_exp_left {
            let poly1 = w_all.get_poly(i, k);
            let poly2 = bold_t_g.get_poly(i, k);
            multiply_add_poly_avx(params, poly_res, poly1, poly2);
        }
        fast_reduce(&mut res);
        r_all[i] = &r_all[i] + &res.raw();

        let gadget = gadget_invert_transposed_alloc(&r_bar_all[i + 1], params.t_exp_left).ntt();
        bold_t_bar_g.copy_into(&gadget, i, 0);

        let mut rotated_gadget = PolyMatrixNTT::zero(&params, 1, params.t_exp_left);
        apply_automorph_ntt(&params, &packing_params.tables, &gadget, &mut rotated_gadget, 2*params.poly_len - gen_pows[(params.poly_len - i) % params.poly_len]);        
        bold_t_bar_g_prime.copy_into(&rotated_gadget, i, 0);

        let mut res = PolyMatrixNTT::zero(&params, 1, 1);
        let poly_res = res.get_poly_mut(0, 0);
        for k in 0..params.t_exp_left {
            let poly1 = w_bar_all.get_poly(i, k);
            let poly2 = bold_t_bar_g.get_poly(i, k);
            multiply_add_poly_avx(params, poly_res, poly1, poly2);
        }
        fast_reduce(&mut res);
        r_bar_all[i] = &r_bar_all[i] + &res.raw();
    }

    let bold_t_hat = gadget_invert_transposed_alloc(&r_bar_all[0], params.t_exp_left).ntt();

    let mut res = PolyMatrixNTT::zero(&params, 1, 1);
    let poly_res = res.get_poly_mut(0, 0);
    for k in 0..params.t_exp_left {
        let poly1 = v_mask.get_poly(0, k);
        let poly2 = bold_t_hat.get_poly(0, k);
        multiply_add_poly_avx(params, poly_res, poly1, poly2);
    }
    fast_reduce(&mut res);
    r_all[0] = &r_all[0] + &res.raw();

    PrecompInsPIR {
        a_hat: r_all[0].clone(),
        bold_t_condensed: condense_matrix(&params, &bold_t_g_prime),
        bold_t_bar_condensed: condense_matrix(&params, &bold_t_bar_g_prime),
        bold_t_hat_condensed: condense_matrix(&params, &bold_t_hat),
    }
}

pub fn full_packing_with_preprocessing_online<'a>(
    packing_params: &'a PackParams<'a>,
    precomp_inspiring: &PrecompInsPIR<'a>,
    b_poly: &PolyMatrixRaw<'a>,
    y_all_condensed: &PolyMatrixNTT<'a>,
    y_bar_all_condensed: &PolyMatrixNTT<'a>,
    z_body_condensed: &PolyMatrixNTT<'a>,
) -> PolyMatrixRaw<'a> {
    let params = packing_params.params;

    let a_hat = &precomp_inspiring.a_hat;
    let bold_t_condensed = &precomp_inspiring.bold_t_condensed;
    let bold_t_bar_condensed = &precomp_inspiring.bold_t_bar_condensed;
    let bold_t_hat_condensed = &precomp_inspiring.bold_t_hat_condensed;

    let now = Instant::now();

    let num_to_pack = params.poly_len;
    let num_to_pack_half = num_to_pack >> 1;

    let mut sum_poly = PolyMatrixNTT::zero(&params, 1, 1);

    let addition_capacity = 1 << (64 - 2 * params.q2_bits - 1);
    let mut num_added = 0;

    for i in 0..num_to_pack_half - 1 {
        let mut temp_poly = PolyMatrixNTT::zero(&params, 1, 1);
        fast_multiply_no_reduce_in_range_generic(
            &params,
            &mut temp_poly,
            &y_all_condensed,
            &bold_t_condensed,
        i*params.t_exp_left,
            i*params.t_exp_left,
            params.t_exp_left,
        );
        fast_add_into_no_reduce(&mut sum_poly, &temp_poly);
        num_added += params.t_exp_left;

        let mut temp_poly = PolyMatrixNTT::zero(&params, 1, 1);
        fast_multiply_no_reduce_in_range_generic(
            &params,
            &mut temp_poly,
            &y_bar_all_condensed,
            &bold_t_bar_condensed,
            i*params.t_exp_left,
        i*params.t_exp_left,
            params.t_exp_left,
        );
        fast_add_into_no_reduce(&mut sum_poly, &temp_poly);
        num_added += params.t_exp_left;

        if (num_added >= addition_capacity) || (i == num_to_pack_half - 2) {
            fast_reduce(&mut sum_poly);
            num_added = 0;
        }
    }

    let mut temp_poly = PolyMatrixNTT::zero(&params, 1, 1);
    fast_multiply_no_reduce_in_range_generic(
        &params,
        &mut temp_poly,
        &z_body_condensed,
        &bold_t_hat_condensed,
        0,
    0,
        params.t_exp_left,
    );
    fast_add_into(&mut sum_poly, &temp_poly);    

    let mut final_b_poly = PolyMatrixRaw::zero(&params, 1, 1);
    add_raw(&mut final_b_poly, &b_poly, &sum_poly.raw());

    let mut packed_raw = PolyMatrixRaw::zero(&params, 2, 1);
    packed_raw.copy_into(&a_hat, 0, 0);
    packed_raw.copy_into(&final_b_poly, 1, 0);

    debug!("Packing online took {} us", now.elapsed().as_micros());

    packed_raw
}

pub fn full_packing_with_preprocessing_online_without_rotations<'a>(
    packing_params: &'a PackParams<'a>,
    precomp_inspiring: &PrecompInsPIR<'a>,
    b_poly: &PolyMatrixRaw<'a>,
    y_condensed: &PolyMatrixNTT<'a>,
    z_body_condensed: &PolyMatrixNTT<'a>,
) -> PolyMatrixRaw<'a> {
    let params = packing_params.params;

    let a_hat = &precomp_inspiring.a_hat;
    let bold_t_condensed = &precomp_inspiring.bold_t_condensed;
    let bold_t_bar_condensed = &precomp_inspiring.bold_t_bar_condensed;
    let bold_t_hat_condensed = &precomp_inspiring.bold_t_hat_condensed;

    let now = Instant::now();

    let num_to_pack = params.poly_len;
    let num_to_pack_half = num_to_pack >> 1;

    let mut sum_poly = PolyMatrixNTT::zero(&params, 1, 1);

    let addition_capacity = 1 << (64 - 2 * params.q2_bits - 1);
    let mut num_added = 0;

    for i in 0..num_to_pack_half - 1 {
        let mut temp_poly = PolyMatrixNTT::zero(&params, 1, 1);
        fast_multiply_no_reduce_in_range_generic(
            &params,
            &mut temp_poly,
            &y_condensed,
            &bold_t_condensed,
	    0,
	    i*params.t_exp_left,
            params.t_exp_left,
        );
        let mut rotated_temp = PolyMatrixNTT::zero(&params, 1, 1);
        apply_automorph_ntt(&params, &packing_params.tables, &temp_poly, &mut rotated_temp, packing_params.gen_pows[i]);        
        
        fast_add_into_no_reduce(&mut sum_poly, &rotated_temp);        
        num_added += params.t_exp_left;

        let mut temp_poly = PolyMatrixNTT::zero(&params, 1, 1);
        fast_multiply_no_reduce_in_range_generic(
            &params,
            &mut temp_poly,
            &y_condensed,
            &bold_t_bar_condensed,
	    0,
	    i*params.t_exp_left,
	    params.t_exp_left,
        );
        let mut rotated_temp = PolyMatrixNTT::zero(&params, 1, 1);
        apply_automorph_ntt(&params, &packing_params.tables, &temp_poly, &mut rotated_temp, 2 * params.poly_len - packing_params.gen_pows[i]);        
        
        fast_add_into_no_reduce(&mut sum_poly, &rotated_temp);        
        num_added += params.t_exp_left;

        if (num_added >= addition_capacity) || (i == num_to_pack_half - 2) {
            fast_reduce(&mut sum_poly);
            num_added = 0;
        }

    }

    let mut temp_poly = PolyMatrixNTT::zero(&params, 1, 1);
    fast_multiply_no_reduce_in_range_generic(
        &params,
        &mut temp_poly,
        &z_body_condensed,
        &bold_t_hat_condensed,
        0,
	0,
        params.t_exp_left,
    );
    fast_add_into(&mut sum_poly, &temp_poly);    

    let mut final_b_poly = PolyMatrixRaw::zero(&params, 1, 1);
    add_raw(&mut final_b_poly, &b_poly, &sum_poly.raw());

    let mut packed_raw = PolyMatrixRaw::zero(&params, 2, 1);
    packed_raw.copy_into(&a_hat, 0, 0);
    packed_raw.copy_into(&final_b_poly, 1, 0);

    debug!("Packing online took {} us", now.elapsed().as_micros());

    packed_raw
}

pub fn packing_fully_online<'a>(
    packing_params: &'a PackParams<'a>,
    offline_packing_keys: &OfflinePackingKeys<'a>,
    a_ct_tilde: &Vec<PolyMatrixNTT<'a>>,
    packing_keys: &PackingKeys<'a>,
    b_poly: &PolyMatrixRaw<'a>,
) -> PolyMatrixRaw<'a> {
    let poly_len = packing_params.params.poly_len;

    // let part5 : Instant;

    let packed = if (packing_params.num_to_pack == poly_len) && (a_ct_tilde.len() > poly_len / 2) {
        // let part4 = Instant::now();

        let w_all = offline_packing_keys.w_all.as_ref().unwrap();
        let w_bar_all = offline_packing_keys.w_bar_all.as_ref().unwrap();
        let v_mask = offline_packing_keys.v_mask.as_ref().unwrap();

        let precomp_inspiring = full_packing_with_preprocessing_offline(
            &packing_params,
            &w_all,
            &w_bar_all,
            &v_mask,
            &a_ct_tilde,
        );

        let y_all_condensed = &packing_keys.y_all_condensed.as_ref().unwrap();
        let y_bar_all_condensed = &packing_keys.y_bar_all_condensed.as_ref().unwrap();
        let z_body_condensed = &packing_keys.z_body_condensed.as_ref().unwrap();        

        // println!("part4: {} us", part4.elapsed().as_micros());

        // part5 = Instant::now();
        full_packing_with_preprocessing_online(
            &packing_params,
            &precomp_inspiring, // parts generated offline
            b_poly,
            y_all_condensed,
            y_bar_all_condensed,
            z_body_condensed, // parts generated online
        )

    } else {
        // let part4 = Instant::now();
        
        let w_all = offline_packing_keys.w_all.as_ref().unwrap();
        let y_all_condensed = &packing_keys.y_all_condensed.as_ref().unwrap();

        let precomp_inspiring =
            packing_with_preprocessing_offline(&packing_params, &w_all, &a_ct_tilde);
            // println!("part4: {} us", part4.elapsed().as_micros());
        
        // part5 = Instant::now();
        packing_with_preprocessing_online(
            &packing_params, &precomp_inspiring, // parts generated offline
            b_poly,
            &y_all_condensed, // parts generated online
        )
    };
    // println!("part5: {} us", part5.elapsed().as_micros());

    packed
}

pub fn packing_fully_online_without_rotations<'a>(
    packing_params: &'a PackParams<'a>,
    offline_packing_keys: &OfflinePackingKeys<'a>,
    a_ct_tilde: &Vec<PolyMatrixNTT<'a>>,
    packing_keys: &PackingKeys<'a>,
    b_poly: &PolyMatrixRaw<'a>,
) -> PolyMatrixRaw<'a> {
    let poly_len = packing_params.params.poly_len;
    let packed = if (packing_params.num_to_pack == poly_len) && (a_ct_tilde.len() > poly_len / 2) {

        let w_all = offline_packing_keys.w_all.as_ref().unwrap();
        let w_bar_all = offline_packing_keys.w_bar_all.as_ref().unwrap();
        let v_mask = offline_packing_keys.v_mask.as_ref().unwrap();
        
        let precomp_inspiring = full_packing_with_preprocessing_offline(
            &packing_params,
            &w_all,
            &w_bar_all,
            &v_mask,
            &a_ct_tilde,
        );
        let y_condensed = packing_keys.y_body_condensed.as_ref().unwrap(); 
        let z_body_condensed = packing_keys.z_body_condensed.as_ref().unwrap(); 
        full_packing_with_preprocessing_online_without_rotations(
            &packing_params,
            &precomp_inspiring, // parts generated offline
            b_poly,
            y_condensed,
            z_body_condensed, // parts generated online
        )
    } else {
        let w_all = offline_packing_keys.w_all.as_ref().unwrap();
        let precomp_inspiring =
            packing_with_preprocessing_offline_without_rotations(&packing_params, &w_all, &a_ct_tilde);
        let y_condensed = packing_keys.y_body_condensed.as_ref().unwrap(); 
        packing_with_preprocessing_online_without_rotations(
            &packing_params, &precomp_inspiring, // parts generated offline
            b_poly,
            y_condensed, // parts generated online
        )
    };
    packed
}


pub fn half_packing_fully_online<'a>(
    half_packing_params: &'a PackParams<'a>,
    w_all: &PolyMatrixNTT<'a>,
    a_ct_tilde: &Vec<PolyMatrixNTT<'a>>,
    y_all_condensed: &PolyMatrixNTT<'a>,
    b_poly: &PolyMatrixRaw<'a>,
) -> PolyMatrixRaw<'a> {
    let precomp_inspiring =
        packing_with_preprocessing_offline(&half_packing_params, &w_all, &a_ct_tilde);
    packing_with_preprocessing_online(
        &half_packing_params, &precomp_inspiring, // parts generated offline
        b_poly, y_all_condensed, // parts generated online
    )
}

pub fn half_packing_fully_online_without_rotations<'a>(
    half_packing_params: &'a PackParams<'a>,
    w_all: &PolyMatrixNTT<'a>,
    a_ct_tilde: &Vec<PolyMatrixNTT<'a>>,
    // y_condensed: &PolyMatrixNTT<'a>,
    packing_keys: &PackingKeys<'a>,
    b_poly: &PolyMatrixRaw<'a>,
) -> PolyMatrixRaw<'a> {
    let precomp_inspiring =
        packing_with_preprocessing_offline_without_rotations(&half_packing_params, &w_all, &a_ct_tilde);
    let y_condensed = packing_keys.y_body_condensed.as_ref().unwrap();
    packing_with_preprocessing_online_without_rotations(
        &half_packing_params, &precomp_inspiring, // parts generated offline
        b_poly, y_condensed, // parts generated online
    )
}

pub fn generate_ksk_body<'a>(
    params: &'a Params,
    sk_reg: &PolyMatrixRaw<'a>,
    gen: usize,
    mask: &PolyMatrixNTT<'a>,
    rng: &mut ChaCha20Rng,
) -> PolyMatrixNTT<'a> {
    let tau_sk_reg = automorph_alloc(&sk_reg, gen);
    let minus_s_times_mask = &sk_reg.ntt() * &(-mask);
    let error_poly = PolyMatrixRaw::noise(
        &params,
        1,
        params.t_exp_left,
        &DiscreteGaussian::init(params.noise_width),
        rng,
    );
    // let error_poly = PolyMatrixRaw::zero(&params, 1, params.t_exp_left);
    let g_exp_ntt = build_gadget(&params, 1, params.t_exp_left).ntt();
    let ksk = &tau_sk_reg.ntt() * &g_exp_ntt;
    let body = &minus_s_times_mask + &error_poly.ntt();
    let result = &body + &ksk;
    result
}


pub fn query_gen<'a>(
    packing_params: &'a PackParams,
    sk_reg: &PolyMatrixRaw<'a>,
    w_mask: &PolyMatrixNTT<'a>,
    v_mask: &PolyMatrixNTT<'a>,
    messages: &Vec<u64>,
    a_ct_tilde: &Vec<PolyMatrixNTT<'a>>,
    rng_y: &mut ChaCha20Rng,
    rng_z: &mut ChaCha20Rng,
) -> (PolyMatrixRaw<'a>, PolyMatrixNTT<'a>, PolyMatrixNTT<'a>) {
    let params = packing_params.params;
    let gen = packing_params.gen_pows[1];
    // let num_to_pack = packing_params.num_to_pack;
    let gamma = a_ct_tilde.len();
    let y_body = generate_ksk_body(&params, &sk_reg, gen, &w_mask, rng_y);
    let z_body = generate_ksk_body(&params, &sk_reg, 2 * params.poly_len - 1, &v_mask, rng_z);

    let mut b_poly = PolyMatrixRaw::zero(&params, 1, 1);
    for i in 0..gamma {
        let mut pt: PolyMatrixRaw<'_> = PolyMatrixRaw::zero(&params, 1, 1);
        pt.get_poly_mut(0, 0)[0] = rescale(messages[i], params.pt_modulus, params.modulus);

        let mut b_ntt = &sk_reg.ntt() * &(-&a_ct_tilde[i]);
        fast_add_into(&mut b_ntt, &pt.ntt());

        b_poly.get_poly_mut(0, 0)[i] = b_ntt.raw().get_poly(0, 0)[0];
    }
    (b_poly, y_body, z_body)
}

#[cfg(test)]
mod test {

    use rand::SeedableRng;
    use spiral_rs::{client::Client, number_theory::invert_uint_mod};

    use crate::{
        client::raw_generate_expansion_params, params::params_for_scenario,
        server::generate_y_constants,
    };
    use rand_chacha::ChaCha20Rng;

    use crate::scheme::{V_SEED, W_SEED};

    use super::*;

    #[test]
    fn test_packing() {
        // let params = get_test_params();
        let (params, _) = params_for_scenario(1 << 30, 1);
        let mut client = Client::init(&params);
        client.generate_secret_keys();
        let y_constants = generate_y_constants(&params);

        let pack_seed = [1u8; 32];
        let cts_seed = [2u8; 32];
        let mut ct_pub_rng = ChaCha20Rng::from_seed(cts_seed);

        let pack_pub_params = raw_generate_expansion_params(
            &params,
            client.get_sk_reg(),
            params.poly_len_log2,
            params.t_exp_left,
            &mut ChaCha20Rng::from_entropy(),
            &mut ChaCha20Rng::from_seed(pack_seed),
        );

        let num_to_pack = params.poly_len; // can't change this
        let mut gold = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..num_to_pack {
            gold.data[i] = i as u64 % params.pt_modulus;
        }

        // generate poly_len ciphertexts
        let mut v_ct = Vec::new();
        let mut b_values = Vec::new();
        for i in 0..num_to_pack {
            let mut pt: PolyMatrixRaw<'_> = PolyMatrixRaw::zero(&params, 1, 1);
            let val = gold.data[i];
            let scale_k = params.modulus / params.pt_modulus;
            let mod_inv = invert_uint_mod(num_to_pack as u64, params.modulus).unwrap();
            let val_to_enc = multiply_uint_mod(val * scale_k, mod_inv, params.modulus);
            pt.data[0] = val_to_enc;
            let ct = client.encrypt_matrix_reg(
                &pt.ntt(),
                &mut ChaCha20Rng::from_entropy(),
                &mut ct_pub_rng,
            );
            let mut ct_raw = ct.raw();

            // get the b value
            b_values.push(ct_raw.get_poly(1, 0)[0]);

            // zero out all of the second poly
            ct_raw.get_poly_mut(1, 0).fill(0);
            v_ct.push(ct_raw.ntt());
        }

        let now = Instant::now();
        let packed = pack_lwes(
            &params,
            &b_values,
            &v_ct,
            &[],
            params.poly_len,
            &pack_pub_params,
            &y_constants,
        );

        println!("Packing took {} us", now.elapsed().as_micros());

        // let packed_raw = packed.raw();
        // println!("packed_0: {:?}", &packed_raw.get_poly(0, 0)[..10]);
        // assert_eq!(packed_raw.get_poly(0, 0)[0], 47649720264253743u64);

        // decrypt + decode
        let dec = client.decrypt_matrix_reg(&packed);
        let dec_raw = dec.raw();

        // rescale
        let mut rescaled = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..params.poly_len {
            rescaled.data[i] = rescale(dec_raw.data[i], params.modulus, params.pt_modulus);
        }

        println!("rescaled: {:?}", &rescaled.as_slice()[..num_to_pack.min(50)]);
        assert_eq!(rescaled.as_slice()[..num_to_pack], gold.as_slice()[..num_to_pack]);
    }

    #[test]
    fn test_packing_inspir() {
        let (params, _) = params_for_scenario(1 << 30, 1);
        let mut client = Client::init(&params);
        client.generate_secret_keys();
        let cts_seed = [2u8; 32];
        let mut ct_pub_rng = ChaCha20Rng::from_seed(cts_seed);

        let num_to_pack = 16;

        let mut gold = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..num_to_pack {
            gold.data[i] = i as u64 % params.pt_modulus;
        }

        // generate poly_len ciphertexts
        let mut v_ct = Vec::new();
        let mut b_values = Vec::new();
        for i in 0..num_to_pack {
            let mut pt: PolyMatrixRaw<'_> = PolyMatrixRaw::zero(&params, 1, 1);
            let val = gold.data[i];
            let scale_k = params.modulus / params.pt_modulus;
            // let mod_inv = invert_uint_mod(num_to_pack as u64, params.modulus).unwrap();
            let mod_inv = 1;
            let val_to_enc = multiply_uint_mod(val * scale_k, mod_inv, params.modulus);
            pt.data[0] = val_to_enc;
            let ct = client.encrypt_matrix_reg(
                &pt.ntt(),
                &mut ChaCha20Rng::from_entropy(),
                &mut ct_pub_rng,
            );
            let mut ct_raw = ct.raw();

            // get the b value
            b_values.push(ct_raw.get_poly(1, 0)[0]);

            // zero out all of the second poly
            ct_raw.get_poly_mut(1, 0).fill(0);
            v_ct.push(ct_raw.ntt());
        }
        let mut a_ct_tilde = Vec::new();
        for i in 0..num_to_pack {
            a_ct_tilde.push(v_ct[i].submatrix(0, 0, 1, 1));
        }
        let mut b_poly = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..num_to_pack {
            b_poly.get_poly_mut(0, 0)[i] = b_values[i];
        }

        let packing_params = if num_to_pack <= (params.poly_len >> 1) {
            PackParams::new(&params, params.poly_len >> 1)
        } else {
            PackParams::new(&params, params.poly_len)
        };
        let now = Instant::now();

        let offline_packing_keys = if num_to_pack <= params.poly_len / 2 {
            OfflinePackingKeys::init(&packing_params, W_SEED)
        } else {
            OfflinePackingKeys::init_full(&packing_params, W_SEED, V_SEED)
        };

        let packing_keys = if num_to_pack <= (params.poly_len >> 1) {
            PackingKeys::init(&packing_params, &client.get_sk_reg(), W_SEED)
        } else {
            PackingKeys::init_full(&packing_params, &client.get_sk_reg(), W_SEED, V_SEED)
        };

        let packed = packing_fully_online(
            &packing_params,
            &offline_packing_keys,
            &a_ct_tilde,
            &packing_keys,
            &b_poly,
        )
        .ntt();

        println!("Packing took {} us", now.elapsed().as_micros());

        // decrypt + decode
        let dec = client.decrypt_matrix_reg(&packed);
        let dec_raw = dec.raw();

        // rescale
        let mut rescaled = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..params.poly_len {
            rescaled.data[i] = rescale(dec_raw.data[i], params.modulus, params.pt_modulus);
        }

        println!("rescaled: {:?}", &rescaled.as_slice()[..num_to_pack.min(50)]);
        assert_eq!(rescaled.as_slice()[..num_to_pack], gold.as_slice()[..num_to_pack]);
    }

    #[test]
    fn test_precompute_packing() {
        let (params, _) = params_for_scenario(1 << 30, 1);
        let mut client = Client::init(&params);
        client.generate_secret_keys();
        let y_constants = generate_y_constants(&params);

        // let log2_num_to_pack = 11;
        // let num_to_pack = 1 << log2_num_to_pack;

        // let total_num_to_pack = 2048;
        let num_to_pack = params.poly_len as usize;

        let pack_seed = [1u8; 32];
        let cts_seed = [2u8; 32];
        let mut ct_pub_rng = ChaCha20Rng::from_seed(cts_seed);

        let pack_pub_params = raw_generate_expansion_params(
            &params,
            client.get_sk_reg(),
            params.poly_len_log2,
            params.t_exp_left,
            &mut ChaCha20Rng::from_entropy(),
            &mut ChaCha20Rng::from_seed(pack_seed),
        );

        let mut fake_pack_pub_params = pack_pub_params.clone();
        // zero out all of the second rows
        for i in 0..pack_pub_params.len() {
            for col in 0..pack_pub_params[i].cols {
                fake_pack_pub_params[i].get_poly_mut(1, col).fill(0);
            }
        }

        let mut pack_pub_params_row_1s = pack_pub_params.clone();
        for i in 0..pack_pub_params.len() {
            pack_pub_params_row_1s[i] =
                pack_pub_params[i].submatrix(1, 0, 1, pack_pub_params[i].cols);
            pack_pub_params_row_1s[i] = condense_matrix(&params, &pack_pub_params_row_1s[i]);
        }

        let mut gold = Vec::new();
        for i in 0..params.poly_len {
            gold.push(i as u64 % params.pt_modulus);
        }

        // generate poly_len ciphertexts
        let mut v_ct = Vec::new();
        let mut b_values = Vec::new();
        let scale_k = params.modulus / params.pt_modulus;
        for i in 0..num_to_pack {
            let mut pt = PolyMatrixRaw::zero(&params, 1, 1);
            let val = gold[i] as u64 % params.pt_modulus;
            let mod_inv = invert_uint_mod(num_to_pack as u64, params.modulus).unwrap();
            // let mod_inv = 1;
            let val_to_enc = multiply_uint_mod(val * scale_k, mod_inv, params.modulus);
            pt.data[0] = val_to_enc;
            let ct = client.encrypt_matrix_reg(
                &pt.ntt(),
                &mut ChaCha20Rng::from_entropy(),
                &mut ct_pub_rng,
            );
            let mut ct_raw = ct.raw();

            // get the b value
            b_values.push(ct_raw.get_poly(1, 0)[0]);

            // zero out all of the second poly
            ct_raw.get_poly_mut(1, 0).fill(0);
            v_ct.push(ct_raw.ntt());
        }

        let mut a_ct_tilde = Vec::new();
        for i in 0..num_to_pack {
            a_ct_tilde.push(v_ct[i].submatrix(0, 0, 1, 1));
        }

        let now = Instant::now();
        let (precomp_res, precomp_vals, precomp_tables) = precompute_pack(
            &params,
            params.poly_len_log2,
            &v_ct,
            &fake_pack_pub_params,
            &y_constants,
        );

        println!("Precomputing for packing took {} us", now.elapsed().as_micros());

        // let sk_reg = client.get_sk_reg();
        // let y_body = generate_ksk_body(&params, sk_reg, packing_params.gen_pows[1],
        // &w_mask, &mut ChaCha20Rng::from_entropy()); let z_body =
        // generate_ksk_body(&params, sk_reg, 2*params.poly_len - 1, &v_mask, &mut
        // ChaCha20Rng::from_entropy());
        let mut b_poly = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..params.poly_len {
            // b_poly.get_poly_mut(0, 0)[i] = multiply_uint_mod(b_values[i], params.poly_len
            // as u64, params.modulus);
            b_poly.get_poly_mut(0, 0)[i] = b_values[i];
        }

        let now = Instant::now();
        let packed = pack_using_precomp_vals(
            &params,
            params.poly_len_log2,
            &pack_pub_params_row_1s,
            &b_values,
            &precomp_res,
            &precomp_vals,
            &precomp_tables,
            &y_constants,
        );

        println!("Packing took {} us", now.elapsed().as_micros());

        // decrypt + decode
        let dec = client.decrypt_matrix_reg(&packed.ntt());
        let dec_raw = dec.raw();

        // rescale
        let mut rescaled = Vec::new();
        let mut errors = Vec::new();
        for i in 0..params.poly_len {
            let mu = dec_raw.data[i];
            let message = rescale(mu, params.modulus, params.pt_modulus);
            rescaled.push(message);
            // let err = mu - message * scale_k;
            // let error = min(err, params.modulus-err);
            // errors.push(error);

            let mut err = mu - message * scale_k;
            if err > params.modulus / 2 {
                err = (params.modulus / 2) - err;
            }
            errors.push(err);
        }
        println!("modulus bitlength: {}", params.modulus_log2);
        let max_error = errors.iter().max().unwrap();
        println!("error bitlength: {}", (*max_error as f64).log2());

        println!("rescaled: {:?}", &rescaled.as_slice()[..50]);
        assert_eq!(rescaled.as_slice()[..num_to_pack], gold.as_slice()[..num_to_pack]);
    }

    use rayon::current_num_threads;
    
    #[test]
    fn test_packing_with_preprocessing() {
        println!("Number of threads in the thread pool: {}", current_num_threads());

        let (params, _) = params_for_scenario(1 << 30, 1);
        
        let max_num_to_pack = 2048 as usize;
        let gamma = max_num_to_pack as usize;

        let packing_params = PackParams::new(&params, max_num_to_pack);

        let w_mask = &PolyMatrixNTT::random(&params, 1, params.t_exp_left);
        let v_mask = &PolyMatrixNTT::random(&params, 1, params.t_exp_left);

        let w_mask_condensed = condense_matrix(&params, &w_mask); 
        let v_mask_condensed = condense_matrix(&params, &v_mask); 

        // generate poly_len ciphertexts
        let mut a_ct_tilde = Vec::new();
        for _ in 0..gamma {
            a_ct_tilde.push(PolyMatrixNTT::random(&params, 1, 1));
        }

        let now_offline = Instant::now();

        let precomp_inspiring = if max_num_to_pack < params.poly_len {
            let w_all = generate_rotations(&packing_params, &w_mask);
                packing_with_preprocessing_offline(&packing_params, &w_all, &a_ct_tilde)
            
                // packing_with_preprocessing_offline_without_rotations(&packing_params, &w_all, &a_ct_tilde)
            
        } else {
            let (w_all, w_bar_all) = generate_rotations_double(&packing_params, &w_mask_condensed);
            full_packing_with_preprocessing_offline(
                &packing_params,
                &w_all,
                &w_bar_all,
                &v_mask_condensed,
                &a_ct_tilde,
            )
            // full_packing_with_preprocessing_offline_without_rotations(
            //     &packing_params,
            //     &w_all,
            //     &w_bar_all,
            //     &v_mask,
            //     &a_ct_tilde,
            // )
        };

        println!("Offline packing took {} us", now_offline.elapsed().as_micros());

        ///////////////////////// Client Query

        let mut client = Client::init(&params);
        client.generate_secret_keys();
        let sk_reg = client.get_sk_reg();
        let mut messages = vec![0; gamma];
        for i in 0..gamma {
            messages[i] = (i + 1) as u64 % params.pt_modulus;
        }

        let (b_poly, y_body, z_body) = query_gen(
            &packing_params,
            sk_reg,
            &w_mask,
            &v_mask,
            &messages,
            &a_ct_tilde,
            &mut ChaCha20Rng::from_entropy(),
            &mut ChaCha20Rng::from_entropy(),
        );
        // let mut y_bar = PolyMatrixNTT::zero(&params, 1, params.t_exp_left);
        // apply_automorph_ntt(
        //     &params,
        //     &packing_params.tables,
        //     &y_body,
        //     &mut y_bar,
        //     2 * params.poly_len - 1,
        // );        
        let y_condensed = condense_matrix(&params, &y_body);
        // let y_bar_condensed = condense_matrix(&params, &y_bar);
        let z_body_condensed = condense_matrix(&params, &z_body);

        ///////////////////////// Online server

        let now_online = Instant::now();

        let packed = if max_num_to_pack <= params.poly_len/2 {
            let now_prerotate = Instant::now();
            let y_all_condensed = generate_rotations(&packing_params, &y_condensed);
            println!("Prerotations took {} us", now_prerotate.elapsed().as_micros());
            packing_with_preprocessing_online(&packing_params, &precomp_inspiring, &b_poly, &y_all_condensed)
            // packing_with_preprocessing_online_without_rotations(&packing_params, &precomp_inspiring, &b_poly, &y_condensed)
        } else {
            let now_prerotate = Instant::now();
            let (y_all_condensed, y_bar_all_condensed) = generate_rotations_double(&packing_params, &y_condensed);
            println!("Prerotations took {} us", now_prerotate.elapsed().as_micros());
            full_packing_with_preprocessing_online(
                &packing_params,
                &precomp_inspiring,
                &b_poly,
                &y_all_condensed,
                &y_bar_all_condensed,
                &z_body_condensed,
            )
            // full_packing_with_preprocessing_online_without_rotations(
            //     &packing_params,
            //     &precomp_inspiring,
            //     &b_poly,
            //     &y_condensed,
            //     // &y_bar_condensed,
            //     &z_body_condensed,
            // )
        };
        println!("Online packing took {} us", now_online.elapsed().as_micros());

        ///////////////////////// Result Extract

        let dec: PolyMatrixNTT<'_> = client.decrypt_matrix_reg(&packed.ntt());
        let dec_raw = dec.raw();
        // let scale = params.modulus / params.pt_modulus;

        let mut rescaled = Vec::new();
        // let mut error = Vec::new();
        for i in 0..gamma {
            let mu = dec_raw.get_poly(0, 0)[i];
            let message = rescale(mu, params.modulus, params.pt_modulus);
            rescaled.push(message);

            // // let scaled_message = rescale(message, params.pt_modulus,
            // params.modulus); let mut err = mu - message * scale;
            // if err > params.modulus / 2 {
            //     err = params.modulus - err;
            // }
            // error.push(err);
        }
        // let max_value = error.iter().max().unwrap();
        // let inf_norm_error = *max_value;
        assert_eq!(rescaled.as_slice()[..gamma], messages.as_slice()[..gamma]);
        // print!("error bitlength: {} ", (inf_norm_error as f64).log2());
    }

    #[test]
    fn test_polynomial_mult() {
        let (params, _) = params_for_scenario(1 << 30, 8);
        println!("cnt: {}", params.crt_count);
        let t = 1;

        let p1 = PolyMatrixNTT::random(&params, 1, t);
        let p2 = PolyMatrixNTT::random(&params, t, 1);

        let mut p3 = PolyMatrixNTT::zero(&params, t, 1);

        let p1_condensed = condense_matrix(&params, &p1);
        let p2_condensed = condense_matrix(&params, &p2);
        for j in 0..t {
            p3.get_poly_mut(j, 0).copy_from_slice(&p2.get_poly(j, 0)[..]);
        }

        let now = Instant::now();
        let mut res1 = PolyMatrixNTT::zero(&params, 1, 1);
        for _ in 0..1000 {
            res1 = PolyMatrixNTT::zero(&params, 1, 1);
            fast_multiply_no_reduce(&params, &mut res1, &p1_condensed, &p2_condensed);
            // fast_reduce(&mut res1);
        }
        println!("Method1: {}", now.elapsed().as_micros());

        let now = Instant::now();
        let mut res2 = PolyMatrixNTT::zero(&params, 1, 1);
        for _ in 0..1000 {
            res2 = PolyMatrixNTT::zero(&params, 1, 1);
            let res_poly = res2.get_poly_mut(0, 0);
            for k in 0..t {
                let poly1 = p1.get_poly(0, k);
                let poly2 = p3.get_poly(k, 0);
                multiply_add_poly_avx(&params, res_poly, poly1, poly2);
            }
            // fast_reduce(&mut res2);
        }
        println!("Method2: {}", now.elapsed().as_micros());

        assert_eq!(res1.get_poly(0, 0)[0], res2.get_poly(0, 0)[0]);
    }

    #[test]
    fn test_single_packing() {
        let (params, _) = params_for_scenario(1 << 30, 1);
        let mut client = Client::init(&params);
        client.generate_secret_keys();

        let pack_seed = [1u8; 32];
        let cts_seed = [2u8; 32];
        let mut ct_pub_rng = ChaCha20Rng::from_seed(cts_seed);

        let pack_pub_params = raw_generate_expansion_params(
            &params,
            client.get_sk_reg(),
            params.poly_len_log2,
            params.t_exp_left,
            &mut ChaCha20Rng::from_entropy(),
            &mut ChaCha20Rng::from_seed(pack_seed),
        );

        // generate 1 ciphertext
        let sentinel_val = 99;
        let mut v_ct = Vec::new();
        for _i in 0..1 {
            let mut pt = PolyMatrixRaw::zero(&params, 1, 1);
            let val = sentinel_val % params.pt_modulus;
            let scale_k = params.modulus / params.pt_modulus;
            let mod_inv = invert_uint_mod(params.poly_len as u64, params.modulus).unwrap();
            let val_to_enc = multiply_uint_mod(val * scale_k, mod_inv, params.modulus);
            pt.data[0] = val_to_enc;
            let ct = client.encrypt_matrix_reg(
                &pt.ntt(),
                &mut ChaCha20Rng::from_entropy(),
                &mut ct_pub_rng,
            );

            v_ct.push(ct);
        }
        let now = Instant::now();
        let packed = pack_single_lwe(&params, &pack_pub_params, &v_ct[0]);
        println!("Packing took {} us", now.elapsed().as_micros());

        let packed_raw = packed.raw();
        println!("packed_0: {:?}", &packed_raw.get_poly(0, 0)[..10]);

        // decrypt + decode
        let dec = client.decrypt_matrix_reg(&packed);
        let dec_raw = dec.raw();

        // rescale
        let mut rescaled = PolyMatrixRaw::zero(&params, 1, 1);
        for i in 0..params.poly_len {
            rescaled.data[i] = rescale(dec_raw.data[i], params.modulus, params.pt_modulus);
        }

        println!("rescaled: {:?}", &rescaled.as_slice()[..50]);
        assert_eq!(rescaled.as_slice()[0], sentinel_val);
    }

    #[test]
    fn test_automorph_tables() {
        let (params, _) = params_for_scenario(1 << 30, 1);

        let now = Instant::now();
        let tables = generate_automorph_tables_brute_force(&params);
        println!("Generating tables took {} us", now.elapsed().as_micros());

        let mut rng = ChaCha20Rng::from_seed([7u8; 32]);

        for t in ([3, 9, 17, 33, 65, 129, 257, 513, 1025, 2049]).into_iter().rev() {
            let poly = PolyMatrixRaw::random_rng(&params, 1, 1, &mut rng);
            let poly_ntt = poly.ntt();

            let poly_auto = automorph_alloc(&poly, t);
            // let poly_auto_ntt = poly_auto.ntt();
            let mut poly_auto_ntt_using_tables = PolyMatrixNTT::zero(&params, 1, 1);
            apply_automorph_ntt(&params, &tables, &poly_ntt, &mut poly_auto_ntt_using_tables, t);

            let poly_auto_ntt_using_tables_back = poly_auto_ntt_using_tables.raw();
            println!("poly_ntt: {:?}", &poly_auto.as_slice()[..10]);
            println!(
                "poly_auto_ntt_using_tables: {:?}",
                &poly_auto_ntt_using_tables_back.as_slice()[..10]
            );

            assert_eq!(
                &poly_auto.as_slice(),
                &poly_auto_ntt_using_tables_back.as_slice(),
                "t: {}",
                t
            );
        }
    }

    #[test]
    fn test_random() {
        let (params, _) = params_for_scenario(1 << 30, 1);
        let mut rng = ChaCha20Rng::from_seed([7u8; 32]);
        let pt = PolyMatrixRaw::random_rng(&params, 1, 1, &mut rng);
        println!("{:?}", &pt.as_slice()[..10]);
    }
    use crate::params::params_for_scenario_simplepir_small;

    #[test]
    fn test_new_mult(){
        let (params, _) = params_for_scenario_simplepir_small(1<<30, 1);
        let t = 3;
        let ell = 3;
        // let p1 = PolyMatrixNTT::random(&params, 1, ell*t);
        // let p2 = PolyMatrixNTT::random(&params, 1, ell*t);

        // let p1 = PolyMatrixRaw::single_value(&params, 5).ntt();
        // let p2 = PolyMatrixRaw::single_value(&params, 3).ntt();

        let p1_raw = PolyMatrixRaw::random(&params, 1, ell*t);
        // for i in 0..ell*t {
        //     p1_raw.data[i*params.poly_len] = (i+1) as u64;
        // }
        let p1 = p1_raw.ntt();

        let p2_raw = PolyMatrixRaw::random(&params, 1, ell*t);
        // for i in 0..ell*t {
        //     p2_raw.data[i*params.poly_len] = 3;
        // }
        let p2 = p2_raw.ntt();

        let p1_condensed = condense_matrix(&params, &p1);
        let p2_condensed = condense_matrix(&params, &p2);

        let now = Instant::now();
        let mut res1 = PolyMatrixNTT::zero(&params, 1, 1);
        for l in 0..ell {
            let mut temp = PolyMatrixNTT::zero(&params, 1, 1);
            fast_multiply_no_reduce_in_range_generic(&params, &mut temp, &p1_condensed, &p2_condensed, t*l, t*l, t);
            fast_add_into(&mut res1, &temp);
            // fast_reduce(&mut res1);
            // println!("{}", res1.get_poly(0, 0)[0]);
            // println!("{}", res1.raw().get_poly(0, 0)[0]);
        }
        println!("Method1: {}", now.elapsed().as_micros());

        let now = Instant::now();
        let mut res2 = PolyMatrixNTT::zero(&params, 1, 1);
        let res_poly = res2.get_poly_mut(0, 0);
        
        for l in 0..ell {
            for k in 0..t {
                let poly1 = p1.get_poly(0, l*t+k);
                let poly2 = p2.get_poly(0, l*t+k);
                multiply_add_poly_avx(&params, res_poly, poly1, poly2);
            }
            // println!("{}", res_poly[0]);
            // println!("{}", res2.clone().raw().get_poly(0, 0)[0]);            
        }
        fast_reduce(&mut res2);
        println!("Method2: {}", now.elapsed().as_micros());

        // println!("{:?}\n{:?}", res1.raw().get_poly(0, 0)[0], res2.raw().get_poly(0, 0)[0]);
        assert_eq!(res1.raw().get_poly(0, 0)[..10], res2.raw().get_poly(0, 0)[..10]);
    }
}
