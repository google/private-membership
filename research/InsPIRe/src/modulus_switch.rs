use spiral_rs::{arith::*, params::*, poly::*, util::*};

pub trait ModulusSwitch<'a> {
    fn switch(&self, q_1: u64, q_2: u64) -> Vec<u8>;
    fn switch_and_keep(&self, q_1: u64, q_2: u64, num_to_keep: usize) -> Vec<u8>;
    fn recover(params: &'a Params, q_1: u64, q_2: u64, ciphertext: &[u8]) -> Self;
    fn recover_how_many(params: &'a Params, q_1: u64, q_2: u64, how_many: usize, ciphertext: &[u8]) -> Self;
}

impl<'a> ModulusSwitch<'a> for PolyMatrixRaw<'a> {
    fn switch(&self, q_1: u64, q_2: u64) -> Vec<u8> {
        assert_eq!(self.rows, 2);
        assert_eq!(self.cols, 1);
        let q_1_bits = (q_1 as f64).log2().ceil() as usize;
        let q_2_bits = (q_2 as f64).log2().ceil() as usize;
        let total_sz_bits = (q_2_bits + q_1_bits) * self.params.poly_len;
        let total_sz_bytes = (total_sz_bits + 7) / 8;

        let mut res = vec![0u8; total_sz_bytes];
        let mut bit_offs = 0;
        let (row_0, row_1) = (self.get_poly(0, 0), self.get_poly(1, 0));
        for z in 0..self.params.poly_len {
            let val = row_0[z];
            let val_rescaled = rescale(val, self.params.modulus, q_2);
            write_arbitrary_bits(&mut res, val_rescaled, bit_offs, q_2_bits);
            bit_offs += q_2_bits;
        }
        for z in 0..self.params.poly_len {
            let val = row_1[z];
            let val_rescaled = rescale(val, self.params.modulus, q_1);
            write_arbitrary_bits(&mut res, val_rescaled, bit_offs, q_1_bits);
            bit_offs += q_1_bits;
        }

        res
        // self.as_slice().to_vec()
    }

    fn switch_and_keep(&self, q_1: u64, q_2: u64, num_to_keep: usize) -> Vec<u8> {
        assert_eq!(self.rows, 2);
        assert_eq!(self.cols, 1);
        assert!(num_to_keep <= self.params.poly_len);
        let q_1_bits = (q_1 as f64).log2().ceil() as usize;
        let q_2_bits = (q_2 as f64).log2().ceil() as usize;
        let total_sz_bits = ((q_2_bits * self.params.poly_len + q_1_bits * num_to_keep + 63) / 64) * 64; // rounding up to multiple of 64
        let total_sz_bytes = (total_sz_bits + 7) / 8;

        let mut res = vec![0u8; total_sz_bytes];
        let mut bit_offs = 0;
        let (row_0, row_1) = (self.get_poly(0, 0), self.get_poly(1, 0));
        for z in 0..self.params.poly_len {
            let val = row_0[z];
            let val_rescaled = rescale(val, self.params.modulus, q_2);
            write_arbitrary_bits(&mut res, val_rescaled, bit_offs, q_2_bits);
            bit_offs += q_2_bits;
        }
        for z in 0..num_to_keep {
            let val = row_1[z];
            let val_rescaled = rescale(val, self.params.modulus, q_1);
            write_arbitrary_bits(&mut res, val_rescaled, bit_offs, q_1_bits);
            bit_offs += q_1_bits;
        }

        res
        // self.as_slice().to_vec()
    }    

    fn recover(params: &'a Params, q_1: u64, q_2: u64, ciphertext: &[u8]) -> Self {
        let q_1_bits = (q_1 as f64).log2().ceil() as usize;
        let q_2_bits = (q_2 as f64).log2().ceil() as usize;
        let total_sz_bits = (q_2_bits + q_1_bits) * params.poly_len;
        let total_sz_bytes = (total_sz_bits + 7) / 8;
        assert_eq!(ciphertext.len(), total_sz_bytes);

        let mut res = PolyMatrixRaw::zero(params, 2, 1);
        let mut bit_offs = 0;
        let (row_0, row_1) = res.data.as_mut_slice().split_at_mut(params.poly_len);
        for z in 0..params.poly_len {
            let val = read_arbitrary_bits(&ciphertext, bit_offs, q_2_bits);
            row_0[z] = rescale(val, q_2, params.modulus);
            bit_offs += q_2_bits;
        }
        for z in 0..params.poly_len {
            let val = read_arbitrary_bits(&ciphertext, bit_offs, q_1_bits);
            row_1[z] = rescale(val, q_1, params.modulus);
            bit_offs += q_1_bits;
        }

        res
        // let mut res = PolyMatrixRaw::zero(params, 2, 1);
        // res.as_mut_slice().copy_from_slice(ciphertext);
        // res
    }

    fn recover_how_many(params: &'a Params, q_1: u64, q_2: u64, how_many: usize, ciphertext: &[u8]) -> Self {
        let q_1_bits = (q_1 as f64).log2().ceil() as usize;
        let q_2_bits = (q_2 as f64).log2().ceil() as usize;
        let total_sz_bits = ((q_2_bits * params.poly_len + q_1_bits * how_many + 63) / 64) * 64;
        let total_sz_bytes = (total_sz_bits + 7) / 8;
        assert_eq!(ciphertext.len(), total_sz_bytes);
        assert!(how_many <= params.poly_len);

        let mut res = PolyMatrixRaw::zero(params, 2, 1);
        let mut bit_offs = 0;
        let (row_0, row_1) = res.data.as_mut_slice().split_at_mut(params.poly_len);
        for z in 0..params.poly_len {
            let val = read_arbitrary_bits(&ciphertext, bit_offs, q_2_bits);
            row_0[z] = rescale(val, q_2, params.modulus);
            bit_offs += q_2_bits;
        }
        for z in 0..how_many {
            let val = read_arbitrary_bits(&ciphertext, bit_offs, q_1_bits);
            row_1[z] = rescale(val, q_1, params.modulus);
            bit_offs += q_1_bits;
        }

        res
        // let mut res = PolyMatrixRaw::zero(params, 2, 1);
        // res.as_mut_slice().copy_from_slice(ciphertext);
        // res
    }    
}
