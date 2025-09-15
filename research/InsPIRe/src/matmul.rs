use log::debug;

extern "C" {
    fn matMulVecPacked(out: *mut u32, a: *const u32, b: *const u32, a_rows: usize, a_cols: usize);
    fn matMulVecPacked2(
        out: *mut u32,
        a: *const u32,
        b_full: *const u32,
        a_rows: usize,
        a_cols: usize,
    );
    fn matMulVecPacked4(
        out: *mut u32,
        a: *const u32,
        b_full: *const u32,
        a_rows: usize,
        a_cols: usize,
    );
    fn matMulVecPacked8(
        out: *mut u32,
        a: *const u32,
        b_full: *const u32,
        a_rows: usize,
        a_cols: usize,
    );
}

/// Returns `out = a * b` where `a` is a matrix of shape `(a_rows, a_cols)` and
/// `b` is a vector of shape `(b_rows, b_cols)`. `out` is a vector of shape
/// `(a_rows, b_cols)`.
///
/// `a_cols * 4` must equal `b_rows`. `b_cols` must be 1. `out` must have length
/// greater than or equal to `a_rows + 8`.
pub fn matmul_vec_packed(
    out: &mut [u32],
    a: &[u32],
    b: &[u32],
    a_rows: usize,
    a_cols: usize,
    b_rows: usize,
    b_cols: usize,
) {
    debug!(
        "matmul_vec_packed: a_rows: {}, a_cols: {}, b_rows: {}, b_cols: {}",
        a_rows, a_cols, b_rows, b_cols
    );
    assert_eq!(a.len(), a_rows * a_cols);
    assert_eq!(b.len(), b_rows * b_cols);
    assert_eq!(a_cols * 4, b_rows);
    assert!(out.len() >= a_rows + 8);

    unsafe {
        if b_cols == 1 {
            matMulVecPacked(out.as_mut_ptr(), a.as_ptr(), b.as_ptr(), a_rows, a_cols);
        } else if b_cols == 2 {
            matMulVecPacked2(out.as_mut_ptr(), a.as_ptr(), b.as_ptr(), a_rows, a_cols);
        } else if b_cols == 4 {
            matMulVecPacked4(out.as_mut_ptr(), a.as_ptr(), b.as_ptr(), a_rows, a_cols);
        // } else if b_cols == 6 {
        //     matMulVecPacked6(out.as_mut_ptr(), a.as_ptr(), b.as_ptr(), a_rows, a_cols);
        } else if b_cols == 8 {
            matMulVecPacked8(out.as_mut_ptr(), a.as_ptr(), b.as_ptr(), a_rows, a_cols);
        } else {
            panic!("b_cols must be 1, 2, 4, or 8");
        }
    }
}

#[cfg(test)]
mod test {
    use std::time::Instant;

    use super::*;

    #[test]
    #[ignore]
    fn test_matmul_vec_packed() {
        let a_rows = 32768;
        let a_cols = 8192;
        let b_rows = a_cols * 4;

        let mut a = vec![0u32; a_rows * a_cols];
        for i in 0..a.len() {
            a[i] = fastrand::u32(..);
        }

        let mut b = vec![0u32; b_rows];
        for i in 0..b.len() {
            b[i] = fastrand::u32(..);
        }

        let mut out = vec![0u32; a_rows + 8];

        let now = Instant::now();
        matmul_vec_packed(&mut out, &a, &b, a_rows, a_cols, b_rows, 1);
        debug!("matmul_vec_packed: {:?}", now.elapsed());

        let now = Instant::now();
        matmul_vec_packed(&mut out, &a, &b, a_rows, a_cols, b_rows, 1);
        debug!("matmul_vec_packed: {:?}", now.elapsed());
    }

    #[test]
    #[ignore]
    fn test_matmul_vec_packed_8() {
        // 32 GB:
        // matmul_vec_packed: a_rows: 131072, a_cols: 65536, b_rows: 262144, b_cols: 8
        // let a_rows = 131072;
        // let a_cols = 65536;

        // 10922 x 3642

        // let a_rows = 131072;
        // let a_cols = 131072 / 4;

        //185363, m=185365

        let a_rows = 184320; //131072;
        let a_true_cols = 184320 + 1000; //262144;
        let b_cols = 8;

        let a_cols = (a_true_cols) / 4;
        let b_rows = a_cols * 4;

        // let mut a_rows = 131072; //(185363 / 8) * 8;
        // let mut a_cols = (131072 + 131072 / 16) / 4; //(131072 + 511) / 4; //(185363 / 4 / 8) * 8;
        // debug!("a_rows: {}, a_cols: {}", a_rows, a_cols);
        // a_rows = a_rows / 8 * 8;
        // a_cols = a_cols / 8 * 8;
        debug!("a_rows: {}, a_cols: {}", a_rows, a_cols);

        // let a_rows = ((131072.0 * (4.0 / 3.0)) as usize / 8) * 8;
        // let a_cols = ((131072.0 / 4.0 * (4.0 / 3.0)) as usize / 8) * 8;

        let mut a = vec![0u32; a_rows * a_cols];
        for i in 0..a.len() {
            a[i] = fastrand::u32(..);
        }

        let mut b = vec![0u32; b_rows * b_cols];
        for i in 0..b.len() {
            b[i] = fastrand::u32(..);
        }

        let mut out = vec![0u32; (a_rows + 8) * b_cols];

        let now = Instant::now();
        matmul_vec_packed(&mut out, &a, &b, a_rows, a_cols, b_rows, b_cols);
        debug!("matmul_vec_packed: {:?}", now.elapsed());

        // let now = Instant::now();
        // matmul_vec_packed(&mut out, &a, &b, a_rows, a_cols, b_rows, b_cols);
        // debug!("matmul_vec_packed: {:?}", now.elapsed());
    }
}
