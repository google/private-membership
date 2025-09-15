use spiral_rs::{
    params::*,
    poly::*,
};

use spiral_rs::poly::{PolyMatrix, PolyMatrixNTT};


pub fn cooley_tukey<'a>(params: &Params, y_points: Vec<PolyMatrixNTT<'a>>, monomials: &[PolyMatrixNTT<'a>]) -> Vec<PolyMatrixNTT<'a>> {
    let num_points = y_points.len();
    let poly_len = params.poly_len;
    assert_eq!(num_points & (num_points-1), 0);
    if num_points == 1 {
        let mut res: Vec<PolyMatrixNTT<'a>> = Vec::new();
        res.push(y_points[0].clone());
        return res;
    }
    let num_points_half = num_points / 2; 

    let y_points_even = y_points.iter().enumerate().filter(|(i, _)| i % 2 == 0).map(|(_, item)| item.clone()).collect();
    let y_points_odd = y_points.iter().enumerate().filter(|(i, _)| i % 2 == 1).map(|(_, item)| item.clone()).collect();

    let a_even = cooley_tukey(&params, y_points_even, monomials);
    let a_odd = cooley_tukey(&params, y_points_odd, monomials);
    let mut result : Vec<PolyMatrixNTT> = Vec::with_capacity(num_points);
    for i in 0..num_points_half {
        result.push(a_even[i].clone());
    }
    for i in 0..num_points_half {
        result.push(a_even[i].clone());
    }

    let mut temp = PolyMatrixNTT::zero(&params, 1, 1);
    for i in 0..num_points_half {
        let exponent = (2 * poly_len - 2 * poly_len * i / num_points) % (2 * poly_len);
        let mono = &monomials[exponent];
        
        multiply(&mut temp, &mono, &a_odd[i]);
        add_into(&mut result[i], &temp);
        sub_into(&mut result[i + num_points_half], &temp);
    }
    
    result
}

pub fn mono_mult(modulus: u64, poly_len: usize, c: &mut [u64], a: &[u64], exponent: usize){
    for i in 0..poly_len {
        let new_place = (i+exponent) % (2*poly_len);
        if new_place < poly_len {
            c[new_place] = a[i];
        } else {
            c[new_place - poly_len] = modulus - a[i];
        }
    }
}

pub fn cooley_tukey_without_ntt<'a>(params: &Params, y_points: Vec<PolyMatrixRaw<'a>>) -> Vec<PolyMatrixRaw<'a>> {
    let num_points = y_points.len();
    let poly_len = params.poly_len;
    assert_eq!(num_points & (num_points-1), 0);
    if num_points == 1 {
        let mut res: Vec<PolyMatrixRaw<'a>> = Vec::new();
        res.push(y_points[0].clone());
        return res;
    }
    let num_points_half = num_points / 2; 

    let y_points_even = y_points.iter().enumerate().filter(|(i, _)| i % 2 == 0).map(|(_, item)| item.clone()).collect();
    let y_points_odd = y_points.iter().enumerate().filter(|(i, _)| i % 2 == 1).map(|(_, item)| item.clone()).collect();

    let a_even = cooley_tukey_without_ntt(&params, y_points_even);
    let a_odd = cooley_tukey_without_ntt(&params, y_points_odd);
    let mut result : Vec<PolyMatrixRaw> = Vec::with_capacity(num_points);
    for i in 0..num_points_half {
        result.push(a_even[i].clone());
    }
    for i in 0..num_points_half {
        result.push(a_even[i].clone());
    }

    let mut temp = PolyMatrixRaw::zero(&params, 1, 1);
    let mut temp_poly = temp.get_poly_mut(0, 0);
    for i in 0..num_points_half {
        let exponent = (2 * poly_len - 2 * poly_len * i / num_points) % (2 * poly_len);
        // let mono = &monomials[exponent];
        
        // multiply(&mut temp, &mono, &a_odd[i]);
        // add_into(&mut result[i], &temp);
        // sub_into(&mut result[i + num_points_half], &temp);

        mono_mult(params.modulus, params.poly_len, &mut temp_poly, &a_odd[i].get_poly(0, 0), exponent);

        let mut this_result_poly = result[i].get_poly_mut(0, 0);
        add_poly_into(&params, &mut this_result_poly, temp_poly);

        let mut this_result_poly = result[i+num_points_half].get_poly_mut(0, 0);
        sub_poly_into(&params, &mut this_result_poly, temp_poly);

    }
    
    result
}


// If coeffs represents a polynomial, evaluate that polynomial at P(x) = X^{monomial_exponent}
pub fn evalauate_at_monomial_poly<'a>(params : &'a Params, coeffs : &[PolyMatrixNTT<'a>], monomial_exponent: usize, monomials: &[PolyMatrixNTT<'a>]) -> PolyMatrixNTT<'a>{
    let poly_eval_degree = coeffs.len();
    let mut result = PolyMatrixNTT::zero(&params, 1, 1);
    for i in 0..poly_eval_degree {
        let exponent = (monomial_exponent * i) % (2 * params.poly_len);
        let mut res = PolyMatrixNTT::zero(&params, 1, 1);
        multiply(&mut res, &coeffs[i], &monomials[exponent]);
        add_into(&mut result, &res);
    }
    return result
} 