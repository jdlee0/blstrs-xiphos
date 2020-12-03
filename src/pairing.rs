use crate::{Fp12, G1Affine, G2Affine};

use core::borrow::Borrow;
use fff::Field;

use blst::*;

/// Execute a complete pairing operation `(p, q)`.
pub fn pairing(p: G1Affine, q: G2Affine) -> Fp12 {
    let mut tmp = blst_fp12::default();
    unsafe { blst_miller_loop(&mut tmp, &q.0, &p.0) };

    let mut out = blst_fp12::default();
    unsafe { blst_final_exp(&mut out, &tmp) };

    out.into()
}

/// Execute a complete multi-pairing operation Prod_i (p_i, q_i)
pub fn multi_pairing<I, J>(ps: I, qs: J) -> Fp12
where
  I: IntoIterator,
  I::Item: Borrow<G1Affine>,
  J: IntoIterator,
  J::Item: Borrow<G2Affine>,
{
    let mut ps = ps.into_iter();
    let mut qs = qs.into_iter();

    // Get sizes
    let (p_lo, p_hi) = ps.by_ref().size_hint();
    let (q_lo, q_hi) = qs.by_ref().size_hint();

    // They should all be equal
    assert_eq!(p_lo, q_lo);
    assert_eq!(p_hi, Some(p_lo));
    assert_eq!(q_hi, Some(q_lo));

    let mut acc = Fp12::one();
    // NB: blst internally has a function miller_loop_n, which implements this
    // loop with additional efficiency savings.
    for (p, q) in ps.zip(qs) {
        let mut tmp = Fp12::default();

        unsafe { blst_miller_loop(&mut tmp.0, &q.borrow().0, &p.borrow().0) };
        acc *= tmp;
    }

    unsafe { blst_final_exp(&mut acc.0, &acc.0) };   
    acc
}
