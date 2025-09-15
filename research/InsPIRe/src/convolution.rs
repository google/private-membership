//! Implements convolutions between two integer sequences.

// Say n is the sequence length, and the inputs are mod p and mod q. Take FFT of
// both sequences mod some NTT-friendly modulus that is (log p + log q + log n)
// bits long. Then convolve using standard FFT, multiply pointwise, and take
// inverse FFT.

use spiral_rs::{arith::barrett_coeff_u64, params::Params, poly::*};

use crate::server::ToU64;

static DEFAULT_MODULI: [u64; 2] = [268369921u64, 249561089u64];

pub struct Convolution {
    n: usize,
    params: Params,
}

impl Convolution {
    pub fn params_for(n: usize) -> Params {
        Params::init(
            n,
            &DEFAULT_MODULI,
            6.4,
            1,
            2,
            28,
            0,
            0,
            0,
            0,
            true,
            1,
            1,
            1,
            1,
            0,
        )
    }

    pub fn new(n: usize) -> Self {
        let params = Self::params_for(n);

        Self { n, params }
    }

    pub fn params(&self) -> &Params {
        &self.params
    }

    pub fn ntt(&self, a: &[u32]) -> Vec<u32> {
        assert_eq!(a.len(), self.n);
        let mut inp = PolyMatrixRaw::zero(&self.params, 1, 1);
        for i in 0..self.n {
            inp.data[i] = a[i] as u64;
        }
        let ntt = inp.ntt();
        let mut out = vec![0u32; self.params.crt_count * self.n];
        for i in 0..self.params.crt_count * self.n {
            out[i] = ntt.data[i] as u32;
        }
        out
    }

    pub fn raw(&self, a: &[u32]) -> Vec<u32> {
        assert_eq!(a.len(), self.params.crt_count * self.n);
        let mut inp = PolyMatrixNTT::zero(&self.params, 1, 1);
        for i in 0..self.params.crt_count * self.n {
            inp.data[i] = a[i] as u64;
        }
        let raw = inp.raw();
        let mut out = vec![0u32; self.n];
        for i in 0..self.n {
            let mut val = raw.data[i] as i64;
            assert!(val < self.params.modulus as i64);
            if val > self.params.modulus as i64 / 2 {
                val -= self.params.modulus as i64;
            }
            if val < 0 {
                val += ((self.params.modulus as i64) / (1i64 << 32)) * (1i64 << 32);
                val += 1i64 << 32;
            }
            assert!(val >= 0);
            out[i] = (val % (1i64 << 32)) as u32;
        }
        out
    }

    pub fn pointwise_mul(&self, a: &[u32], b: &[u32]) -> Vec<u32> {
        assert_eq!(a.len(), self.params.crt_count * self.n);
        assert_eq!(b.len(), self.params.crt_count * self.n);

        let mut out = vec![0u32; self.params.crt_count * self.n];
        for m in 0..self.params.crt_count {
            for i in 0..self.n {
                let idx = m * self.n + i;
                let val = a[idx] as u64 * b[idx] as u64;
                out[idx] = barrett_coeff_u64(&self.params, val as u64, m) as u32;
            }
        }
        out
    }

    pub fn convolve(&self, a: &[u32], b: &[u32]) -> Vec<u32> {
        let a_ntt = self.ntt(a);
        let b_ntt = self.ntt(b);
        let res = self.pointwise_mul(&a_ntt, &b_ntt);
        self.raw(&res)
    }
}

pub fn naive_negacyclic_convolve(a: &[u32], b: &[u32]) -> Vec<u32> {
    assert_eq!(a.len(), b.len());
    let n = a.len();
    let mut res = vec![0u32; n];
    for i in 0..n {
        for j in 0..n {
            let mut b_val = b[(n + i - j) % n];
            if i < j {
                b_val = b_val.wrapping_neg();
            }
            res[i] += a[j] * b_val;
        }
    }
    res
}

pub fn negacyclic_matrix_u32(b: &[u32]) -> Vec<u32> {
    let n = b.len();
    let mut res = vec![0u32; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut b_val = b[(n + i - j) % n];
            if i < j {
                b_val = b_val.wrapping_neg();
            }
            res[j * n + i] = b_val; // nb: transposed
        }
    }
    res
}

pub fn negacyclic_perm_u32(a: &[u32]) -> Vec<u32> {
    let n = a.len();
    let mut res = vec![0u32; n];
    res[0] = a[0];
    for i in 1..n {
        res[i] = a[(n - i) % n].wrapping_neg();
    }
    res
}

pub fn naive_multiply_matrices<T: ToU64 + Copy>(
    a: &[u32],
    a_rows: usize,
    a_cols: usize,
    b_t: &[T], // transposed
    b_rows: usize,
    b_cols: usize,
    is_b_transposd: bool,
) -> Vec<u32> {
    // performs wrapping arithmetic

    assert_eq!(a_cols, b_rows);

    // debug!("Multiplying {}x{} by {}x{}", a_rows, a_cols, b_rows, b_cols);

    let mut result = vec![0u32; a_rows * b_cols];
    for i in 0..a_rows {
        for j in 0..b_cols {
            for k in 0..a_cols {
                let a_idx = i * a_cols + k;
                let b_idx = if is_b_transposd {
                    j * b_rows + k // on purpose, since transposed
                } else {
                    k * b_cols + j
                };
                let res_idx = i * b_cols + j;

                unsafe {
                    let a_val = *a.get_unchecked(a_idx);
                    let b_val = (b_t.get_unchecked(b_idx)).to_u64() as u32;

                    result[res_idx] = result[res_idx].wrapping_add(a_val.wrapping_mul(b_val));
                }
            }
        }
    }

    result
}

#[cfg(test)]
mod test {
    use std::time::Instant;

    use log::debug;

    use super::*;

    #[test]
    fn test_ntt_raw() {
        let n = 1 << 10;
        let conv = Convolution::new(n);
        let a = (0..n).map(|_| fastrand::u32(..)).collect::<Vec<_>>();
        let ntt = conv.ntt(&a);
        let raw = conv.raw(&ntt);
        assert_eq!(a, raw);
    }

    #[test]
    fn test_convolve_simple() {
        let n = 1 << 10;
        let conv = Convolution::new(n);
        let mut a = vec![0u32; n];
        (&mut a[..4]).copy_from_slice(&[1, 2, 3, 4]);
        let mut b = vec![0u32; n];
        (&mut b[..4]).copy_from_slice(&[5, 6, 7, 8]);

        let res = conv.convolve(&a, &b);
        // debug!("{:?}", &res[..16]);
        let mut expected = vec![0u32; n];
        (&mut expected[..7]).copy_from_slice(&[5, 16, 34, 60, 61, 52, 32]);
        assert_eq!(res, expected);

        let naive = naive_negacyclic_convolve(&a, &b);
        assert_eq!(expected, naive);

        let nega_b = negacyclic_matrix_u32(&b);
        let naive_matmul = naive_multiply_matrices(&a, 1, n, &nega_b, n, n, false);
        assert_eq!(expected, naive_matmul);
    }

    #[test]
    fn test_convolve_from_nega_matrix() {
        let n = 1 << 10;
        let conv = Convolution::new(n);
        let mut a = vec![0u32; n];
        a.iter_mut().for_each(|x| *x = fastrand::u32(..));
        let mut b = vec![0u32; n];
        b.iter_mut().for_each(|x| *x = fastrand::u32(..256));

        let nega_a = negacyclic_matrix_u32(&a);
        let naive_matmul = naive_multiply_matrices(&nega_a, n, n, &b, n, 1, false);

        let nega_perm_a = negacyclic_perm_u32(&a);

        let now = Instant::now();
        let res = conv.convolve(&nega_perm_a, &b);
        debug!("convolve: {:?}", now.elapsed());

        debug!("{:?}", &res[..16]);
        debug!("{:?}", &naive_matmul[..16]);

        assert_eq!(res, naive_matmul);
    }

    #[test]
    fn test_convolve_negacyclic() {
        let n = 1 << 10;
        let conv = Convolution::new(n);
        let mut a = vec![0u32; n];
        (&mut a[..2]).copy_from_slice(&[1, 2]);
        (&mut a[(n - 2)..]).copy_from_slice(&[3, 4]);
        let mut b = vec![0u32; n];
        (&mut b[..2]).copy_from_slice(&[5, 6]);

        let res = conv.convolve(&a, &b);
        debug!("{:?}...{:?}", &res[..16], &res[(n - 16)..]);
        let mut expected = vec![0u32; n];
        (&mut expected[..3]).copy_from_slice(&[((1u64 << 32) - 19) as u32, 16, 12]);
        (&mut expected[(n - 2)..]).copy_from_slice(&[15, 38]);
        assert_eq!(res, expected);

        let naive = naive_negacyclic_convolve(&a, &b);
        assert_eq!(expected, naive);

        let nega_b = negacyclic_matrix_u32(&b);
        let naive_matmul = naive_multiply_matrices(&a, 1, n, &nega_b, n, n, false);
        assert_eq!(expected, naive_matmul);
    }
}
