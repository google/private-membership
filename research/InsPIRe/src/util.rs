use spiral_rs::{arith::*, params::*, poly::*, util};

use super::server::ToU64;

pub fn negacyclic_perm(a: &[u64], shift: usize, modulus: u64) -> Vec<u64> {
    let n = a.len();
    let mut out = vec![0u64; n];

    for i in 0..shift + 1 {
        out[i] = a[shift - i];
    }

    for i in shift + 1..n {
        out[i] = modulus - (a[n - (i - shift)] % modulus);
        if out[i] == modulus {
            out[i] = 0;
        }
    }

    out
}

pub fn negacyclic_matrix(a: &[u64], modulus: u64, how_many: usize) -> Vec<u64> {
    let n = a.len();
    let mut out = vec![0u64; n * how_many];

    for i in 0..how_many {
        let perm = negacyclic_perm(a, i, modulus);
        for j in 0..n {
            out[j * how_many + i] = perm[j];
        }
    }

    out
}

pub fn get_negacylic<'a>(poly: &PolyMatrixRaw<'a>) -> PolyMatrixRaw<'a> {
    let mut out = poly.clone();

    (&mut out.as_mut_slice()[1..]).reverse();

    for z in 1..out.data.len() {
        out.data[z] = out.params.modulus - out.data[z];
    }

    out
}

pub fn reduce_copy(params: &Params, out: &mut [u64], inp: &[u64]) {
    for n in 0..params.crt_count {
        for z in 0..params.poly_len {
            out[n * params.poly_len + z] = barrett_coeff_u64(params, inp[z], n);
        }
    }
}

pub fn add_into_no_reduce(res: &mut PolyMatrixNTT, a: &PolyMatrixNTT) {
    assert!(res.rows == a.rows);
    assert!(res.cols == a.cols);

    for i in 0..res.rows {
        for j in 0..res.cols {
            let res_poly = res.get_poly_mut(i, j);
            let pol2 = a.get_poly(i, j);
            for (res_poly, pol2) in res_poly.iter_mut().zip(pol2.iter()) {
                *res_poly += *pol2;
            }
        }
    }
}

pub fn add_into_at_no_reduce(
    res: &mut PolyMatrixNTT,
    a: &PolyMatrixNTT,
    t_row: usize,
    t_col: usize,
) {
    let params = res.params;
    for i in 0..a.rows {
        for j in 0..a.cols {
            let res_poly = res.get_poly_mut(t_row + i, t_col + j);
            let pol2 = a.get_poly(i, j);
            for z in 0..params.crt_count * params.poly_len {
                res_poly[z] += pol2[z];
            }
        }
    }
}

pub fn modular_reduce_poly<'a>(a: &mut PolyMatrixNTT<'a>) {
    let params = a.params;
    for i in 0..a.rows {
        for j in 0..a.cols {
            let pol = a.get_poly_mut(i, j);
            for z in 0..params.crt_count * params.poly_len {
                pol[z] = barrett_coeff_u64(params, pol[z], z / params.poly_len);
            }
        }
    }
}

pub fn scalar_multiply_avx(res: &mut PolyMatrixNTT, a: &PolyMatrixNTT, b: &PolyMatrixNTT) {
    assert_eq!(a.rows, 1);
    assert_eq!(a.cols, 1);

    let params = res.params;
    let pol2 = a.get_poly(0, 0);
    for i in 0..b.rows {
        for j in 0..b.cols {
            let res_poly = res.get_poly_mut(i, j);
            let pol1 = b.get_poly(i, j);
            crate::packing::multiply_poly_avx(params, res_poly, pol1, pol2);
        }
    }
}

pub fn concat_horizontal(v_a: &[Vec<u64>], a_rows: usize, a_cols: usize) -> Vec<u64> {
    let mut out = vec![0u64; a_rows * a_cols * v_a.len()];

    for i in 0..a_rows {
        for j in 0..a_cols {
            for k in 0..v_a.len() {
                let idx = i * a_cols + j;
                let out_idx = i * a_cols * v_a.len() + k * a_cols + j;
                out[out_idx] = v_a[k][idx];
            }
        }
    }

    out
}

#[cfg(target_feature = "avx2")]
pub fn is_avx() -> bool {
    true
}

#[cfg(not(target_feature = "avx2"))]
pub fn is_avx() -> bool {
    false
}

pub fn test_params() -> Params {
    let params_str = r#"{
        "n": 1,
        "nu_1": 0,
        "nu_2": 0,
        "p": 256,
        "q2_bits": 22,
        "t_gsw": 3,
        "t_conv": 2,
        "t_exp_left": 2,
        "t_exp_right": 2,
        "instances": 1,
        "db_item_size": 0,
        "version": 2
    }"#;
    let params = util::params_from_json(&params_str);
    params
}

pub fn multiply_matrices_raw_not_transposed<T>(
    params: &Params,
    a: &[u64],
    a_rows: usize,
    a_cols: usize,
    b: &[T], // NOT transposed
    b_rows: usize,
    b_cols: usize,
) -> Vec<u64>
where
    T: ToU64 + Copy,
{
    assert_eq!(a_cols, b_rows);

    let mut result = vec![0u128; a_rows * b_cols];

    for i in 0..a_rows {
        for k in 0..a_cols {
            for j in 0..b_cols {
                let a_idx = i * a_cols + k;
                let b_idx = k * b_cols + j;
                let res_idx = i * b_cols + j;

                unsafe {
                    let a_val = *a.get_unchecked(a_idx);
                    let b_val = (*b.get_unchecked(b_idx)).to_u64();

                    let prod = a_val as u128 * b_val as u128;
                    result[res_idx] += prod;
                }
            }
        }
    }

    let mut result_u64 = vec![0u64; a_rows * b_cols];
    for i in 0..result.len() {
        result_u64[i] = barrett_reduction_u128(params, result[i]);
    }

    result_u64
}

#[cfg(test)]
mod test {
    use spiral_rs::poly::*;

    use super::super::transpose::transpose_generic;
    use super::*;

    #[test]
    fn test_negacyclic_mul_db_col() {
        let params = test_params();
        let pol_a = PolyMatrixRaw::random(&params, 1, 1);
        let pol_b: PolyMatrixRaw<'_> = PolyMatrixRaw::random(&params, 1, 1);
        let a = pol_a.get_poly(0, 0);
        let b = pol_b.get_poly(0, 0);
        let negacylic_a = negacyclic_matrix(&a, params.modulus, params.poly_len);
        let negacyclic_a_t = transpose_generic(&negacylic_a, params.poly_len, params.poly_len);

        // the twist is that b is a 'column of the db'
        // we want to compute (b as a row vector) * (transpose(negacylic_a)))

        assert_eq!(negacylic_a[0], a[0]);
        assert_eq!(
            negacylic_a[params.poly_len],
            (params.modulus - a[params.poly_len - 1]) % params.modulus
        );
        let prod = multiply_matrices_raw_not_transposed(
            &params,
            b,
            1,
            params.poly_len,
            &negacyclic_a_t,
            params.poly_len,
            params.poly_len,
        );

        // we think this is equivalent to:
        // poly_b * poly([a_0 -a_d-1 -a_d-2 ... -a_1])
        // = poly_b * poly(negacyclic_perm(a, 0))

        let transformed_a = negacyclic_perm(a, 0, params.modulus);
        let mut pol_a_transformed = PolyMatrixRaw::zero(&params, 1, 1);
        pol_a_transformed
            .data
            .as_mut_slice()
            .copy_from_slice(&transformed_a);
        let pol_c = (&pol_a_transformed.ntt() * &pol_b.ntt()).raw();
        let c = pol_c.get_poly(0, 0);

        for i in 0..params.poly_len {
            assert_eq!(prod[i] % params.modulus, c[i] % params.modulus, "i = {}", i);
        }
    }

    #[test]
    fn test_negacyclic_perm() {
        let a = vec![1, 2, 3, 4, 5, 6, 7, 8];
        let modulus = 11;
        let shift = 3;

        let result = negacyclic_perm(&a, shift, modulus);
        assert_eq!(result, vec![4, 3, 2, 1, 11 - 8, 11 - 7, 11 - 6, 11 - 5]);
    }

    #[test]
    fn test_crt() {
        let params = test_params();

        let a = (1u64 << 45) + 77;

        let a0 = a % params.moduli[0];
        let a1 = a % params.moduli[1];

        let a_goal = params.crt_compose_2(a0, a1) % params.modulus;
        assert_eq!(a_goal, a);

        let a0_mod = (a0 * 5123) % params.moduli[0];
        let a1_mod = (a1 * 5123) % params.moduli[1];

        let a_modified = params.crt_compose_2(a0_mod, a1_mod) % params.modulus;
        assert_eq!(a_modified, (a * 5123) % params.modulus);
    }

    #[test]
    fn test_negacyclic_mul() {
        let params = test_params();
        let pol_a = PolyMatrixRaw::random(&params, 1, 1);
        let pol_b = PolyMatrixRaw::random(&params, 1, 1);
        let a = pol_a.get_poly(0, 0);
        let b = pol_b.get_poly(0, 0);
        let negacylic_a = negacyclic_matrix(&a, params.modulus, params.poly_len);
        // let negacylic_a_t = transpose_generic(&negacylic_a, params.poly_len, params.poly_len);
        assert_eq!(negacylic_a[0], a[0]);
        assert_eq!(
            negacylic_a[params.poly_len],
            (params.modulus - a[params.poly_len - 1]) % params.modulus
        );
        let prod = multiply_matrices_raw_not_transposed(
            &params,
            b,
            1,
            params.poly_len,
            &negacylic_a,
            params.poly_len,
            params.poly_len,
        );

        let pol_c = (&pol_a.ntt() * &pol_b.ntt()).raw();
        let c = pol_c.get_poly(0, 0);

        for i in 0..params.poly_len {
            assert_eq!(prod[i] % params.modulus, c[i] % params.modulus, "i = {}", i);
        }
    }
}
