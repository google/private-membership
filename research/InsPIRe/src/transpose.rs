pub fn transpose(buf: &[u8], rows: usize, cols: usize, bytes_per_pt_el: usize) -> Vec<u8> {
    let mut out = vec![0u8; buf.len()];

    for i in 0..rows {
        for j in 0..cols {
            for k in 0..bytes_per_pt_el {
                out[j * rows * bytes_per_pt_el + i * bytes_per_pt_el + k] =
                    buf[i * cols * bytes_per_pt_el + j * bytes_per_pt_el + k];
            }
        }
    }

    out
}

pub fn transpose_elems<T: Sized + Copy>(
    buf: &[&[T]],
    inp_rows: usize,
    inp_cols: usize,
) -> Vec<Vec<T>> {
    let mut out = Vec::new();

    for j in 0..inp_cols {
        let mut row = Vec::new();
        for i in 0..inp_rows {
            row.push(buf[i][j]);
        }
        out.push(row);
    }

    out
}

static TRANSPOSE_TILE_SIZE: usize = 32;

pub fn transpose_generic<T>(a: &[T], a_rows: usize, a_cols: usize) -> Vec<T>
where
    T: Sized + Copy + Default,
{
    // assert!(a_rows >= TRANSPOSE_TILE_SIZE);
    // assert!(a_cols >= TRANSPOSE_TILE_SIZE);
    // assert_eq!(a_rows % TRANSPOSE_TILE_SIZE, 0);
    // assert_eq!(a_cols % TRANSPOSE_TILE_SIZE, 0);
    let mut transpose_tile_size = 32.min(a_rows).min(a_cols);
    if a_rows % transpose_tile_size != 0 || a_cols % transpose_tile_size != 0 {
        transpose_tile_size = 1;
    }
    let mut out = vec![T::default(); a_rows * a_cols];

    for i_outer in (0..a_rows).step_by(transpose_tile_size) {
        for j_outer in (0..a_cols).step_by(transpose_tile_size) {
            for i_inner in 0..transpose_tile_size {
                for j_inner in 0..transpose_tile_size {
                    let i = i_outer + i_inner;
                    let j = j_outer + j_inner;
                    out[j * a_rows + i] = a[i * a_cols + j];
                }
            }
        }
    }

    out
}

pub fn transpose_f64(out: &mut [f64], a: &[f64], a_rows: usize, a_cols: usize) {
    assert!(a_rows >= TRANSPOSE_TILE_SIZE);
    assert!(a_cols >= TRANSPOSE_TILE_SIZE);
    assert_eq!(a_rows % TRANSPOSE_TILE_SIZE, 0);
    assert_eq!(a_cols % TRANSPOSE_TILE_SIZE, 0);
    let transpose_tile_size = TRANSPOSE_TILE_SIZE; //32.min(a_rows).min(a_cols);

    for i_outer in (0..a_rows).step_by(transpose_tile_size) {
        for j_outer in (0..a_cols).step_by(transpose_tile_size) {
            for i_inner in 0..transpose_tile_size {
                for j_inner in 0..transpose_tile_size {
                    unsafe {
                        let i = i_outer + i_inner;
                        let j = j_outer + j_inner;

                        // out[j * a_rows + i] = a[i * a_cols + j];
                        *out.get_unchecked_mut(j * a_rows + i) = *a.get_unchecked(i * a_cols + j);
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_transpose() {
        let mut buf = vec![0u8; 16];
        for i in 0..16 {
            buf[i] = i as u8;
        }

        let out = transpose(&buf, 2, 8, 1);

        assert_eq!(
            out,
            vec![0, 8, 1, 9, 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15,]
        );
    }
}
