use std::cmp::min;

use log::debug;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

use spiral_rs::aligned_memory::AlignedMemory64;
use spiral_rs::{
    arith::*, client::*, discrete_gaussian::*, gadget::*, number_theory::*, params::*, poly::*,
};

use crate::packing::PackingType;

use super::convolution::negacyclic_matrix_u32;
use super::{lwe::*, noise_analysis::measure_noise_width_squared, scheme::*, util::*};

pub fn rlwe_to_lwe<'a>(params: &'a Params, ct: &PolyMatrixRaw<'a>, how_many: usize) -> Vec<u64> {
    let a = ct.get_poly(0, 0);
    let mut negacylic_a = negacyclic_matrix(&a, params.modulus, how_many);
    negacylic_a.extend(&ct.get_poly(1, 0)[0..how_many]);

    negacylic_a
}

pub fn pack_query(params: &Params, query: &[u64]) -> AlignedMemory64 {
    let query_packed = if params.crt_count == 2 {
        query
            .iter()
            .enumerate()
            .map(|(_i, x)| {
                let crt0 = (*x) % params.moduli[0];
                let crt1 = (*x) % params.moduli[1];
                crt0 | (crt1 << 32)
            })
            .collect::<Vec<_>>()
    } else {
        query
            .iter()
            .enumerate()
            .map(|(_i, x)| {
                let crt0 = (*x) % params.moduli[0];
                crt0
            })
            .collect::<Vec<_>>()
    };
    let mut aligned_query_packed = AlignedMemory64::new(query_packed.len());
    aligned_query_packed
        .as_mut_slice()
        .copy_from_slice(&query_packed);
    aligned_query_packed
}

pub fn get_reg_sample<'a>(
    params: &'a Params,
    sk_reg: &PolyMatrixRaw<'a>,
    rng: &mut ChaCha20Rng,
    rng_pub: &mut ChaCha20Rng,
) -> PolyMatrixNTT<'a> {
    let a = PolyMatrixRaw::random_rng(params, 1, 1, rng_pub);
    let e = PolyMatrixRaw::noise(
        params,
        1,
        1,
        &DiscreteGaussian::init(params.noise_width),
        rng,
    );
    let b_p = &sk_reg.ntt() * &a.ntt();
    let b = &e.ntt() + &b_p;
    let mut p = PolyMatrixNTT::zero(params, 2, 1);
    p.copy_into(&(-&a).ntt(), 0, 0);
    p.copy_into(&b, 1, 0);
    p
}

pub fn get_fresh_reg_public_key<'a>(
    params: &'a Params,
    sk_reg: &PolyMatrixRaw<'a>,
    m: usize,
    rng: &mut ChaCha20Rng,
    rng_pub: &mut ChaCha20Rng,
) -> PolyMatrixNTT<'a> {
    let mut p = PolyMatrixNTT::zero(params, 2, m);

    for i in 0..m {
        p.copy_into(&get_reg_sample(params, sk_reg, rng, rng_pub), 0, i);
    }
    p
}

pub fn raw_generate_expansion_params<'a>(
    params: &'a Params,
    sk_reg: &PolyMatrixRaw<'a>,
    num_exp: usize,
    m_exp: usize,
    rng: &mut ChaCha20Rng,
    rng_pub: &mut ChaCha20Rng,
) -> Vec<PolyMatrixNTT<'a>> {
    let g_exp = build_gadget(params, 1, m_exp);
    debug!("using gadget base {}", g_exp.get_poly(0, 1)[0]);
    let g_exp_ntt = g_exp.ntt();
    let mut res = Vec::new();

    for i in 0..num_exp {
        let t = (params.poly_len / (1 << i)) + 1;
        let tau_sk_reg = automorph_alloc(&sk_reg, t);
        let prod = &tau_sk_reg.ntt() * &g_exp_ntt;

        // let w_exp_i = client.encrypt_matrix_reg(&prod, rng, rng_pub);
        let sample = get_fresh_reg_public_key(params, &sk_reg, m_exp, rng, rng_pub);
        let w_exp_i = &sample + &prod.pad_top(1);
        res.push(w_exp_i);
    }

    res
}

pub fn decrypt_ct_reg_measured<'a>(
    client: &Client<'a>,
    params: &'a Params,
    ct: &PolyMatrixNTT<'a>,
    coeffs_to_measure: usize,
) -> PolyMatrixRaw<'a> {
    let dec_result = client.decrypt_matrix_reg(ct).raw();

    let mut dec_rescaled = PolyMatrixRaw::zero(&params, dec_result.rows, dec_result.cols);
    for z in 0..dec_rescaled.data.len() {
        dec_rescaled.data[z] = rescale(dec_result.data[z], params.modulus, params.pt_modulus);
    }

    // measure noise width
    let _s_2 = measure_noise_width_squared(params, client, ct, &dec_rescaled, coeffs_to_measure);
    // debug!("log2(measured noise): {}", s_2.log2());

    dec_rescaled
}

pub fn ct_reg_measure<'a>(
    client: &Client<'a>,
    params: &'a Params,
    ct: &PolyMatrixNTT<'a>,
    coeffs_to_measure: usize,
) -> f64 {
    let dec_result = client.decrypt_matrix_reg(ct).raw();

    let mut dec_rescaled = PolyMatrixRaw::zero(&params, dec_result.rows, dec_result.cols);
    for z in 0..dec_rescaled.data.len() {
        dec_rescaled.data[z] = rescale(dec_result.data[z], params.modulus, params.pt_modulus);
    }

    // measure noise width
    // let _s_2 = 
    measure_noise_width_squared(params, client, ct, &dec_rescaled, coeffs_to_measure)
    // debug!("log2(measured noise): {}", s_2.log2());

    // dec_rescaled
}

pub struct YClient<'a> {
    inner: &'a mut Client<'a>,
    params: &'a Params,
    lwe_client: LWEClient,
}

pub fn get_seed(public_seed_idx: u8) -> [u8; 32] {
    let mut seed = STATIC_PUBLIC_SEED;
    seed[0] = public_seed_idx;
    seed
}

pub fn generate_matrix_ring(
    rng_pub: &mut ChaCha20Rng,
    n: usize,
    rows: usize,
    cols: usize,
) -> Vec<u32> {
    assert_eq!(rows % n, 0);
    assert_eq!(cols % n, 0);
    let rows_outer = rows / n;
    let cols_outer = cols / n;

    let mut out = vec![0u32; rows * cols];
    for i in 0..rows_outer {
        for j in 0..cols_outer {
            let mut a = vec![0u32; n];
            for idx in 0..n {
                a[idx] = rng_pub.sample::<u32, _>(rand::distributions::Standard);
            }

            let mat = negacyclic_matrix_u32(&a);
            for k in 0..n {
                for l in 0..n {
                    let idx = (i * n + k) * cols + (j * n + l);
                    out[idx] = mat[k * n + l];
                }
            }
        }
    }

    out
}

impl<'a> YClient<'a> {
    pub fn new(inner: &'a mut Client<'a>, params: &'a Params) -> Self {
        Self {
            inner,
            params,
            lwe_client: LWEClient::new(LWEParams::default()),
        }
    }

    pub fn lwe_client(&self) -> &LWEClient {
        &self.lwe_client
    }

    fn rlwes_to_lwes(&self, ct: &[PolyMatrixRaw<'a>], how_many : usize) -> Vec<u64> {
        let v = ct
            .iter()
            .map(|ct| rlwe_to_lwe(self.params, ct, how_many))
            .collect::<Vec<_>>();
        concat_horizontal(&v, self.params.poly_len + 1, how_many) //self.params.poly_len)
    }

    pub fn generate_query_impl(
        &self,
        public_seed_idx: u8,
        dim_log2: usize,
        packing_type: PackingType,
        index: usize,
    ) -> Vec<PolyMatrixRaw<'a>> {

        let dim_log2 = if (dim_log2 as isize) < 0 { 0 } else {dim_log2};  // checks if dim_log2 went below zero

        // let db_cols = 1 << (self.params.db_dim_2 + self.params.poly_len_log2);
        // let idx_dim1 = index / db_cols;

        let mut rng_pub = ChaCha20Rng::from_seed(get_seed(public_seed_idx));

        // Generate dim1_bits LWE samples under public randomness
        let mut out = Vec::new();

        let scale_k = self.params.modulus / self.params.pt_modulus;

        for i in 0..(1 << dim_log2) {
            let mut scalar = PolyMatrixRaw::zero(self.params, 1, 1);
            let is_nonzero = i == (index / self.params.poly_len);

            if is_nonzero {
                scalar.data[index % self.params.poly_len] =
                    if packing_type == PackingType::CDKS {
                        multiply_uint_mod(
                            invert_uint_mod(self.params.poly_len as u64, self.params.modulus).unwrap(),
                            scale_k,
                            self.params.modulus
                        )
                    } else {
                        scale_k
                    }
            }

            // if public_seed_idx == SEED_0 {
            //     out.push(scalar.pad_top(1));
            //     continue;
            // }

            let ct = {
                let factor = if packing_type == PackingType::InspiRING {
                    1
                } else {
                    invert_uint_mod(self.params.poly_len as u64, self.params.modulus).unwrap()
                };
                self.inner.encrypt_matrix_scaled_reg(
                    &scalar.ntt(),
                    &mut ChaCha20Rng::from_entropy(),
                    &mut rng_pub,
                    factor,
                )
            };

            let ct_raw = ct.raw();

            out.push(ct_raw);
        }

        out
    }

    pub fn generate_query_over_u32 (
        &self,
        public_seed_idx: u8,
        dim_log2: usize,
        index_row: usize,
    ) -> Vec<u64> {
        let lwe_params = LWEParams::default();
        let dim = 1 << (dim_log2 + self.params.poly_len_log2);

        // lwes must be (n + 1) x (dim) matrix
        let mut lwes = vec![0u64; (lwe_params.n + 1) * dim];

        let scale_k = lwe_params.scale_k() as u32;
        let mut vals_to_encrypt = vec![0u32; dim];
        vals_to_encrypt[index_row] = scale_k;

        let mut rng_pub = ChaCha20Rng::from_seed(get_seed(public_seed_idx));

        for i in (0..dim).step_by(lwe_params.n) {
            let out = self
                .lwe_client
                .encrypt_many(&mut rng_pub, &vals_to_encrypt[i..i + lwe_params.n])
                .iter()
                .map(|x| *x as u64)
                .collect::<Vec<_>>();
            assert_eq!(out.len(), (lwe_params.n + 1) * lwe_params.n);
            for r in 0..lwe_params.n + 1 {
                for c in 0..lwe_params.n {
                    lwes[r * dim + i + c] = out[r * lwe_params.n + c];
                }
            }
        }

        lwes
    }

    pub fn generate_query_over_prime (
        &self,
        public_seed_idx: u8,
        dim_log2: usize,
        packing_type: PackingType,
        index_row: usize,
    ) -> Vec<u64> {
        let out = self.generate_query_impl(public_seed_idx, dim_log2, packing_type, index_row);
        let how_many = self.params.poly_len / (1 << (-min(0 as isize, dim_log2 as isize)));
        let lwes = self.rlwes_to_lwes(&out, how_many);
        lwes
    }

    pub fn decode_response(&self, response: &[u64]) -> Vec<u64> {
        debug!("Decoding response: {:?}", &response[..16]);
        let db_cols = 1 << (self.params.db_dim_2 + self.params.poly_len_log2);

        let sk = self.inner.get_sk_reg().as_slice().to_vec();

        let mut out = Vec::new();
        for col in 0..db_cols {
            let mut sum = 0u128;
            for i in 0..self.params.poly_len {
                let v1 = response[i * db_cols + col];
                let v2 = sk[i];
                sum += v1 as u128 * v2 as u128;
            }

            sum += response[self.params.poly_len * db_cols + col] as u128;

            let result = (sum % self.params.modulus as u128) as u64;
            let result_rescaled = rescale(result, self.params.modulus, self.params.pt_modulus);
            out.push(result_rescaled);
        }

        out
    }

    pub fn client(&self) -> &Client<'a> {
        self.inner
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_lwe() {
        let lwe_params = LWEParams::default();
        let client = LWEClient::new(lwe_params.clone());
        let pt = fastrand::u32(0..lwe_params.pt_modulus as u32);
        let scaled_pt = pt.wrapping_mul(lwe_params.scale_k() as u32);
        let ct = client.encrypt(&mut ChaCha20Rng::from_entropy(), scaled_pt);
        let pt_dec = client.decrypt(&ct);
        let result = rescale(pt_dec as u64, lwe_params.modulus, lwe_params.pt_modulus) as u32;
        assert_eq!(result, pt);
    }
}
