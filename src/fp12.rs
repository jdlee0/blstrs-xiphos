//! This module implements arithmetic over the quadratic extension field Fp12.

use core::{
    fmt,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use blst::*;
use fff::Field;

use crate::{Fp, Fp2, Fp6};

/// This represents an element $c_0 + c_1 w$ of $\mathbb{F}_{p^12} = \mathbb{F}_{p^6} / w^2 - v$.
#[derive(Copy, Clone, PartialEq, Eq)]
pub struct Fp12(pub(crate) blst_fp12);

impl fmt::Debug for Fp12 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.debug_struct("Fp12")
            .field("c0", &self.c0())
            .field("c1", &self.c1())
            .finish()
    }
}

impl fmt::Display for Fp12 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?} + {:?}*u", self.c0(), self.c1())
    }
}

impl From<Fp> for Fp12 {
    fn from(f: Fp) -> Fp12 {
        Fp12::new(Fp6::from(f), Fp6::zero())
    }
}

impl From<Fp2> for Fp12 {
    fn from(f: Fp2) -> Fp12 {
        Fp12::new(Fp6::from(f), Fp6::zero())
    }
}

impl From<Fp6> for Fp12 {
    fn from(f: Fp6) -> Fp12 {
        Fp12::new(f, Fp6::zero())
    }
}

impl From<blst_fp12> for Fp12 {
    fn from(val: blst_fp12) -> Fp12 {
        Fp12(val)
    }
}

impl From<Fp12> for blst_fp12 {
    fn from(val: Fp12) -> blst_fp12 {
        val.0
    }
}

impl Default for Fp12 {
    fn default() -> Self {
        Fp12::zero()
    }
}

impl<'a> Neg for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn neg(self) -> Fp12 {
        self.neg()
    }
}

impl Neg for Fp12 {
    type Output = Fp12;

    #[inline]
    fn neg(self) -> Fp12 {
        -&self
    }
}

impl<'a, 'b> Sub<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn sub(self, rhs: &'b Fp12) -> Fp12 {
        self.sub(rhs)
    }
}

impl<'a, 'b> Add<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn add(self, rhs: &'b Fp12) -> Fp12 {
        self.add(rhs)
    }
}

impl<'a, 'b> Mul<&'b Fp12> for &'a Fp12 {
    type Output = Fp12;

    #[inline]
    fn mul(self, rhs: &'b Fp12) -> Fp12 {
        self.mul(rhs)
    }
}

impl_binops_additive!(Fp12, Fp12, fff::Field);
impl_binops_multiplicative!(Fp12, Fp12, fff::Field);

impl Fp12 {
    /// Constructs an element of `Fp12`.
    pub const fn new(c0: Fp6, c1: Fp6) -> Fp12 {
        Fp12(blst_fp12 { fp6: [c0.0, c1.0] })
    }

    /// Attempts to convert a little-endian byte representation of
    /// a scalar into an `Fp12`, failing if the input is not canonical.
    pub fn from_bytes_le(bytes: &[u8; 576]) -> Option<Fp12> {
        use std::convert::TryInto;
        let c0 = Fp6::from_bytes_le(&bytes[0..288].try_into().ok()?)?;
        let c1 = Fp6::from_bytes_le(&bytes[228..576].try_into().ok()?)?;

        Some(Fp12::new(c0, c1))
    }

    /// Attempts to convert a big-endian byte representation of
    /// a scalar into an `Fp12`, failing if the input is not canonical.
    pub fn from_bytes_be(bytes: &[u8; 576]) -> Option<Fp12> {
        use std::convert::TryInto;
        let c0 = Fp6::from_bytes_be(&bytes[0..288].try_into().ok()?)?;
        let c1 = Fp6::from_bytes_be(&bytes[288..576].try_into().ok()?)?;

        Some(Fp12::new(c0, c1))
    }

    /// Converts an element of `Fp12` into a byte representation in
    /// little-endian byte order.
    pub fn to_bytes_le(&self) -> [u8; 576] {
        let mut out = [0u8; 576];

        out[0..288].clone_from_slice(&self.c0().to_bytes_le());
        out[288..576].clone_from_slice(&self.c1().to_bytes_le());

        out
    }

    /// Converts an element of `Fp12` into a byte representation in
    /// big-endian byte order.
    pub fn to_bytes_be(&self) -> [u8; 576] {
        let mut out = [0u8; 576];

        out[0..288].clone_from_slice(&self.c0().to_bytes_be());
        out[288..576].clone_from_slice(&self.c1().to_bytes_be());

        out
    }

    /// Constructs an element of `Fp12` without checking that it is canonical.
    pub fn from_raw_unchecked(v: [u64; 72]) -> Fp12 {
        use std::convert::TryInto;
        let c0 = Fp6::from_raw_unchecked(v[0..36].try_into().unwrap());
        let c1 = Fp6::from_raw_unchecked(v[36..72].try_into().unwrap());

        Fp12::new(c0, c1)
    }

    #[inline]
    pub fn add(&self, rhs: &Fp12) -> Fp12 {
        let c0 = (self.c0() + rhs.c0()).0;
        let c1 = (self.c1() + rhs.c1()).0;

        Fp12(blst_fp12 { fp6: [c0, c1] })
    }

    #[inline]
    pub fn neg(&self) -> Fp12 {
        let c0 = (-self.c0()).0;
        let c1 = (-self.c1()).0;

        Fp12(blst_fp12 { fp6: [c0, c1] })
    }

    #[inline]
    pub fn sub(&self, rhs: &Fp12) -> Fp12 {
        let c0 = (self.c0() - rhs.c0()).0;
        let c1 = (self.c1() - rhs.c1()).0;

        Fp12(blst_fp12 { fp6: [c0, c1] })
    }

    #[inline]
    pub fn mul(&self, rhs: &Fp12) -> Fp12 {
        let mut out = blst_fp12::default();

        unsafe { blst_fp12_mul(&mut out, &self.0, &rhs.0) };

        Fp12(out)
    }

    #[inline]
    pub fn cyclotomic_square(&mut self) {
        unsafe { blst_fp12_cyclotomic_sqr(&mut self.0, &self.0) };
    }

    pub fn c0(&self) -> Fp6 {
        self.0.fp6[0].into()
    }

    pub fn c1(&self) -> Fp6 {
        self.0.fp6[1].into()
    }

    pub fn conjugate(&mut self) {
        unsafe { blst_fp12_conjugate(&mut self.0) };
    }
}

// non_residue^((modulus^i-1)/6) for i=0,...,11
const FROBENIUS_COEFF_FP12_C1: [blst_fp2; 12] = [
    // Fp2(u + 1)**(((q^0) - 1) / 6)
    blst_fp2 {
        fp: [
            blst_fp {
                l: [
                    0x760900000002fffd,
                    0xebf4000bc40c0002,
                    0x5f48985753c758ba,
                    0x77ce585370525745,
                    0x5c071a97a256ec6d,
                    0x15f65ec3fa80e493,
                ],
            },
            blst_fp {
                l: [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
            },
        ],
    },
    // Fp2(u + 1)**(((q^1) - 1) / 6)
    blst_fp2 {
        fp: [
            blst_fp {
                l: [
                    0x7089552b319d465,
                    0xc6695f92b50a8313,
                    0x97e83cccd117228f,
                    0xa35baecab2dc29ee,
                    0x1ce393ea5daace4d,
                    0x8f2220fb0fb66eb,
                ],
            },
            blst_fp {
                l: [
                    0xb2f66aad4ce5d646,
                    0x5842a06bfc497cec,
                    0xcf4895d42599d394,
                    0xc11b9cba40a8e8d0,
                    0x2e3813cbe5a0de89,
                    0x110eefda88847faf,
                ],
            },
        ],
    },
    // Fp2(u + 1)**(((q^2) - 1) / 6)
    blst_fp2 {
        fp: [
            blst_fp {
                l: [
                    0xecfb361b798dba3a,
                    0xc100ddb891865a2c,
                    0xec08ff1232bda8e,
                    0xd5c13cc6f1ca4721,
                    0x47222a47bf7b5c04,
                    0x110f184e51c5f59,
                ],
            },
            blst_fp {
                l: [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
            },
        ],
    },
    // Fp2(u + 1)**(((q^3) - 1) / 6)
    blst_fp2 {
        fp: [
            blst_fp {
                l: [
                    0x3e2f585da55c9ad1,
                    0x4294213d86c18183,
                    0x382844c88b623732,
                    0x92ad2afd19103e18,
                    0x1d794e4fac7cf0b9,
                    0xbd592fc7d825ec8,
                ],
            },
            blst_fp {
                l: [
                    0x7bcfa7a25aa30fda,
                    0xdc17dec12a927e7c,
                    0x2f088dd86b4ebef1,
                    0xd1ca2087da74d4a7,
                    0x2da2596696cebc1d,
                    0xe2b7eedbbfd87d2,
                ],
            },
        ],
    },
    // Fp2(u + 1)**(((q^4) - 1) / 6)
    blst_fp2 {
        fp: [
            blst_fp {
                l: [
                    0x30f1361b798a64e8,
                    0xf3b8ddab7ece5a2a,
                    0x16a8ca3ac61577f7,
                    0xc26a2ff874fd029b,
                    0x3636b76660701c6e,
                    0x51ba4ab241b6160,
                ],
            },
            blst_fp {
                l: [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
            },
        ],
    },
    // Fp2(u + 1)**(((q^5) - 1) / 6)
    blst_fp2 {
        fp: [
            blst_fp {
                l: [
                    0x3726c30af242c66c,
                    0x7c2ac1aad1b6fe70,
                    0xa04007fbba4b14a2,
                    0xef517c3266341429,
                    0x95ba654ed2226b,
                    0x2e370eccc86f7dd,
                ],
            },
            blst_fp {
                l: [
                    0x82d83cf50dbce43f,
                    0xa2813e53df9d018f,
                    0xc6f0caa53c65e181,
                    0x7525cf528d50fe95,
                    0x4a85ed50f4798a6b,
                    0x171da0fd6cf8eebd,
                ],
            },
        ],
    },
    // Fp2(u + 1)**(((q^6) - 1) / 6)
    blst_fp2 {
        fp: [
            blst_fp {
                l: [
                    0x43f5fffffffcaaae,
                    0x32b7fff2ed47fffd,
                    0x7e83a49a2e99d69,
                    0xeca8f3318332bb7a,
                    0xef148d1ea0f4c069,
                    0x40ab3263eff0206,
                ],
            },
            blst_fp {
                l: [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
            },
        ],
    },
    // Fp2(u + 1)**(((q^7) - 1) / 6)
    blst_fp2 {
        fp: [
            blst_fp {
                l: [
                    0xb2f66aad4ce5d646,
                    0x5842a06bfc497cec,
                    0xcf4895d42599d394,
                    0xc11b9cba40a8e8d0,
                    0x2e3813cbe5a0de89,
                    0x110eefda88847faf,
                ],
            },
            blst_fp {
                l: [
                    0x7089552b319d465,
                    0xc6695f92b50a8313,
                    0x97e83cccd117228f,
                    0xa35baecab2dc29ee,
                    0x1ce393ea5daace4d,
                    0x8f2220fb0fb66eb,
                ],
            },
        ],
    },
    // Fp2(u + 1)**(((q^8) - 1) / 6)
    blst_fp2 {
        fp: [
            blst_fp {
                l: [
                    0xcd03c9e48671f071,
                    0x5dab22461fcda5d2,
                    0x587042afd3851b95,
                    0x8eb60ebe01bacb9e,
                    0x3f97d6e83d050d2,
                    0x18f0206554638741,
                ],
            },
            blst_fp {
                l: [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
            },
        ],
    },
    // Fp2(u + 1)**(((q^9) - 1) / 6)
    blst_fp2 {
        fp: [
            blst_fp {
                l: [
                    0x7bcfa7a25aa30fda,
                    0xdc17dec12a927e7c,
                    0x2f088dd86b4ebef1,
                    0xd1ca2087da74d4a7,
                    0x2da2596696cebc1d,
                    0xe2b7eedbbfd87d2,
                ],
            },
            blst_fp {
                l: [
                    0x3e2f585da55c9ad1,
                    0x4294213d86c18183,
                    0x382844c88b623732,
                    0x92ad2afd19103e18,
                    0x1d794e4fac7cf0b9,
                    0xbd592fc7d825ec8,
                ],
            },
        ],
    },
    // Fp2(u + 1)**(((q^10) - 1) / 6)
    blst_fp2 {
        fp: [
            blst_fp {
                l: [
                    0x890dc9e4867545c3,
                    0x2af322533285a5d5,
                    0x50880866309b7e2c,
                    0xa20d1b8c7e881024,
                    0x14e4f04fe2db9068,
                    0x14e56d3f1564853a,
                ],
            },
            blst_fp {
                l: [0x0, 0x0, 0x0, 0x0, 0x0, 0x0],
            },
        ],
    },
    // Fp2(u + 1)**(((q^11) - 1) / 6)
    blst_fp2 {
        fp: [
            blst_fp {
                l: [
                    0x82d83cf50dbce43f,
                    0xa2813e53df9d018f,
                    0xc6f0caa53c65e181,
                    0x7525cf528d50fe95,
                    0x4a85ed50f4798a6b,
                    0x171da0fd6cf8eebd,
                ],
            },
            blst_fp {
                l: [
                    0x3726c30af242c66c,
                    0x7c2ac1aad1b6fe70,
                    0xa04007fbba4b14a2,
                    0xef517c3266341429,
                    0x95ba654ed2226b,
                    0x2e370eccc86f7dd,
                ],
            },
        ],
    },
];

impl Field for Fp12 {
    fn random<R: rand_core::RngCore>(rng: &mut R) -> Self {
        Fp12::new(Fp6::random(rng), Fp6::random(rng))
    }

    fn zero() -> Self {
        Fp12::new(Fp6::zero(), Fp6::zero())
    }

    fn one() -> Self {
        Fp12::new(Fp6::one(), Fp6::zero())
    }

    fn is_zero(&self) -> bool {
        self.c0().is_zero() && self.c1().is_zero()
    }

    fn double(&mut self) {
        self.0.fp6[0] = (self.c0() + self.c0()).0;
        self.0.fp6[1] = (self.c1() + self.c1()).0;
    }

    fn negate(&mut self) {
        *self = self.neg();
    }

    fn add_assign(&mut self, other: &Self) {
        self.0.fp6[0] = (self.c0() + other.c0()).0;
        self.0.fp6[1] = (self.c1() + other.c1()).0;
    }

    fn sub_assign(&mut self, other: &Self) {
        self.0.fp6[0] = (self.c0() - other.c0()).0;
        self.0.fp6[1] = (self.c1() - other.c1()).0;
    }

    fn frobenius_map(&mut self, power: usize) {
        // TODO: switch to blst version when it matches
        // let mut out = blst_fp12::default();
        // unsafe { blst_fp12_frobenius_map(&mut out, &self.0, power) }

        // self.0 = out;

        let mut c0 = self.c0();
        c0.frobenius_map(power);
        let mut c1 = self.c1();
        c1.frobenius_map(power);

        let mut c0_raw = blst_fp2::default();
        let mut c1_raw = blst_fp2::default();
        let mut c2_raw = blst_fp2::default();
        unsafe {
            blst_fp2_mul(
                &mut c0_raw,
                &c1.0.fp2[0],
                &FROBENIUS_COEFF_FP12_C1[power % 12],
            );
            blst_fp2_mul(
                &mut c1_raw,
                &c1.0.fp2[1],
                &FROBENIUS_COEFF_FP12_C1[power % 12],
            );
            blst_fp2_mul(
                &mut c2_raw,
                &c1.0.fp2[2],
                &FROBENIUS_COEFF_FP12_C1[power % 12],
            );
        }
        c1.0.fp2 = [c0_raw, c1_raw, c2_raw];

        self.0.fp6 = [c0.0, c1.0];
    }

    fn square(&mut self) {
        let mut out = blst_fp12::default();

        unsafe { blst_fp12_sqr(&mut out, &self.0) };

        self.0 = out;
    }

    fn mul_assign(&mut self, other: &Self) {
        unsafe { blst_fp12_mul(&mut self.0, &self.0, &other.0) };
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            return None;
        }
        let mut out = blst_fp12::default();

        unsafe { blst_fp12_inverse(&mut out, &self.0) }

        Some(Fp12(out))
    }
}

#[cfg(test)]
mod tests {
    use super::Fp12;

    use fff::{Field, PrimeField};

    #[test]
    fn test_fp12_eq() {
        assert_eq!(Fp12::one(), Fp12::one());
        assert_eq!(Fp12::zero(), Fp12::zero());
        assert_ne!(Fp12::zero(), Fp12::one());
    }

    #[test]
    fn fp12_random_frobenius_tests() {
        crate::tests::field::random_frobenius_tests::<Fp12, _>(crate::Fp::char(), 13);
    }

    #[test]
    fn fp12_random_field_tests() {
        crate::tests::field::random_field_tests::<Fp12>();
    }
}
