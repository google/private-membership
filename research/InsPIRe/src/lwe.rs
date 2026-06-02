use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

use spiral_rs::discrete_gaussian::*;

use super::convolution::negacyclic_matrix_u32;

#[derive(Clone, Debug)]
pub struct LWEParams {
    pub n: usize,
    pub modulus: u64,
    pub pt_modulus: u64,
    pub q2_bits: usize,
    pub noise_width: f64,
}

impl Default for LWEParams {
    fn default() -> Self {
        Self {
            n: 1024,
            modulus: 1u64 << 32,
            pt_modulus: 256,
            q2_bits: 28,
            noise_width: 27.57291103, // 11 * sqrt(2*pi)
        }
    }
}

impl LWEParams {
    pub fn scale_k(&self) -> u64 {
        self.modulus / self.pt_modulus
    }
}

pub struct LWEClient {
    lwe_params: LWEParams,
    sk: Vec<u32>,
}

impl LWEClient {
    pub fn new(lwe_params: LWEParams) -> Self {
        let mut rng = ChaCha20Rng::from_entropy();
        let dg = DiscreteGaussian::init(lwe_params.noise_width);
        let sk = (0..lwe_params.n)
            .map(|_| dg.sample(lwe_params.modulus, &mut rng) as u32)
            .collect::<Vec<_>>();
        Self { lwe_params, sk }
    }

    pub fn get_sk(&self) -> &[u32] {
        &self.sk
    }

    pub fn encrypt(&self, rng_pub: &mut ChaCha20Rng, pt: u32) -> Vec<u32> {
        let dg = DiscreteGaussian::init(self.lwe_params.noise_width);
        let mut rng = ChaCha20Rng::from_entropy();
        let e = dg.sample(self.lwe_params.modulus, &mut rng) as u32;
        let mut ct = Vec::new();
        let mut sum = 0u32;
        for i in 0..self.lwe_params.n {
            let v = rng_pub.sample::<u32, _>(rand::distributions::Standard);
            ct.push(v);
            sum = sum.wrapping_add(v.wrapping_mul(self.sk[i]));
        }
        // neg_sum = 2^32 - sum (mod 2^32)
        let neg_sum = sum.wrapping_neg();
        let b = neg_sum.wrapping_add(pt.wrapping_add(e));
        ct.push(b);

        ct
    }

    pub fn encrypt_many(&self, rng_pub: &mut ChaCha20Rng, v_pt: &[u32]) -> Vec<u32> {
        assert_eq!(v_pt.len(), self.lwe_params.n);
        let mut rng = ChaCha20Rng::from_entropy();
        let dg = DiscreteGaussian::init(self.lwe_params.noise_width);

        let mut a = Vec::new();
        for _ in 0..self.lwe_params.n {
            let v = rng_pub.sample::<u32, _>(rand::distributions::Standard);
            a.push(v);
        }

        let nega_a = negacyclic_matrix_u32(&a);
        let mut last_row = vec![0u32; self.lwe_params.n];

        // correctness test for convolution
        // let conv = Convolution::new(self.lwe_params.n);
        // let a_ntt = conv.ntt(&a);
        // let sk_ntt = conv.ntt(&self.sk);
        // let b_ntt = conv.pointwise_mul(&a_ntt, &sk_ntt);
        // let test_vals = conv.raw(&b_ntt);

        let n = self.lwe_params.n;
        for col in 0..n {
            let mut sum = 0u32;
            for row in 0..n {
                let idx = row * n + col;
                sum = sum.wrapping_add(nega_a[idx].wrapping_mul(self.sk[row]));
            }
            // assert_eq!(sum, test_vals[col]);

            let e = dg.sample(self.lwe_params.modulus, &mut rng) as u32;
            let val = sum.wrapping_neg().wrapping_add(v_pt[col].wrapping_add(e));
            last_row[col] = val;
        }

        // // neg_sum = 2^32 - sum (mod 2^32)
        // let neg_sum = sum.wrapping_neg();
        // let b = neg_sum.wrapping_add(pt.wrapping_add(e));
        // ct.push(b);

        let ct = [nega_a, last_row].concat();
        ct
    }

    pub fn decrypt(&self, ct: &[u32]) -> u32 {
        let mut sum = 0u32;
        for i in 0..self.lwe_params.n {
            let v1 = ct[i];
            let v2 = self.sk[i];
            sum = sum.wrapping_add(v1.wrapping_mul(v2));
        }
        sum = sum.wrapping_add(ct[self.lwe_params.n]);

        sum
    }
}
