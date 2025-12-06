// Code developed by FARAOUN Kamel Mohamed.
// EEDIS-Laboratory. UDL-University. Algeria
// During May 2024.

use crate::{curves::{curve_arithmetics::EcPoint, g2::G2Element, g2_primitives::phi::phi_bls48}, 
            fields::prime_fields::FieldElement, tools::{arithmetic_interface::ArithmeticOperations, 
            recoders::{recod_scalar_gls16, recod_scalar_gls4, recod_scalar_gls8}}};
use num_bigint::BigUint;
use num_traits::{Num, One,ToPrimitive};
use std::{ops::BitAnd, str::FromStr};

use super::phi::phi_bls24;

fn compute_res<const PRAMASIZE: usize, const R: usize, const N: usize,
    const MAX_COEFS_COUNT: usize>(sign: i8, lookup: &[G2Element<PRAMASIZE, R, N, MAX_COEFS_COUNT>; 8],
                                  idx: usize) -> G2Element<PRAMASIZE, R, N, MAX_COEFS_COUNT> {
    match sign {
        1 => lookup[idx],
        _ => lookup[idx].negate()
    }
}

fn create_infinit<const PRAMASIZE:usize, const R: usize, const N: usize, const MAX_COEFS_COUNT: usize>(
    input: &G2Element<PRAMASIZE, R, N, MAX_COEFS_COUNT>) -> G2Element<PRAMASIZE, R, N, MAX_COEFS_COUNT> {
    G2Element {
        point: EcPoint {
            x: input.point.x.one(),
            y: input.point.x.one(),
            z: input.point.x.zero()
        },
        consts: input.consts
    }
}

pub fn gls4_multiply_bls12<const PRAMASIZE:usize, const R: usize, const N: usize, const MAX_COEFS_COUNT: usize>
            (input: &G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>, scalar :&FieldElement<R>) -> G2Element<PRAMASIZE, R, N, MAX_COEFS_COUNT>
{   
    // Constant-Time multiplication for elements in G2 (4-GLS): GLS Implementation of points multiplication on G2, for the BLS12 Curve
    // Joppe W. Bos, Craig Costello, and Michael Naehrig https://eprint.iacr.org/2013/458.pdf
    
    let ff: BigUint = BigUint::from_str("255").unwrap();
    let mut code = recod_scalar_gls4(&scalar.to_big_uint(),&scalar.negate().to_big_uint(),input.consts.u.abs() as u128);

    let infinit = create_infinit(input);
    let mut lookup = [infinit;8];    
    populate_lookup1(input, &mut lookup);

    let mut limb : u8 = (&code).bitand(ff.clone()).to_u8().unwrap();
    let out_sig : i8 =  1 - ((limb & 1) << 1) as i8;
    code = code >> 1;
    limb = (&code).bitand(ff.clone()).to_u8().unwrap();
    let mut sign : i8 = ((limb & 1) << 1) as i8 - 1;
    let mut idx = ((limb & 15) >> 1) as usize;
    let mut result = compute_res(sign, &lookup, idx);
    code = code >> 4;
    while code != BigUint::one() {
        result = result.double();
        limb = (&code).bitand(ff.clone()).to_u8().unwrap();
        sign = ((limb & 1) << 1) as i8 - 1;
        idx =  ((limb & 15) >> 1) as usize;
        result = result.addto(&compute_res(sign, &lookup, idx));
        code = code >> 4;
    }
    if out_sig == 1 {
        result
    } else {
        result.negate()
    }
}

fn populate_lookup1<const PRAMASIZE: usize, const R: usize,
    const N: usize, const MAX_COEFS_COUNT: usize>(input: &G2Element<PRAMASIZE, R, N, MAX_COEFS_COUNT>,
                                                  lookup1: &mut [G2Element<PRAMASIZE, R, N, MAX_COEFS_COUNT>; 8]) {
    let mut p1 = input.phi();
    let p2 = p1.phi();
    let mut p3 = p2.phi();

    if input.consts.u < 0 {
        p1 = p1.negate();
        p3 = p3.negate()
    };

    lookup1[0] = input.clone();
    lookup1[1] = p1.addto(&lookup1[0]);
    lookup1[2] = p2.addto(&lookup1[0]);
    lookup1[3] = p2.addto(&lookup1[1]);
    lookup1[4] = p3.addto(&lookup1[0]);
    lookup1[5] = p3.addto(&lookup1[1]);
    lookup1[6] = p3.addto(&lookup1[2]);
    lookup1[7] = p3.addto(&lookup1[3]);
}

fn create_lookup<const PRAMASIZE: usize, const R: usize,
    const N: usize, const MAX_COEFS_COUNT: usize>(infinit: G2Element<PRAMASIZE, R, N, MAX_COEFS_COUNT>,
                                                  lookup1: &[G2Element<PRAMASIZE, R, N, MAX_COEFS_COUNT>; 8],
                                                  degree: u8) -> [G2Element<PRAMASIZE, R, N, MAX_COEFS_COUNT>; 8] {
    let mut arg = [infinit; 8];

    match degree {
        24 => {
            arg.iter_mut()
                .zip(lookup1.iter().map(|x| phi_bls24(x,4)))
                .for_each(|(dest, src)| *dest = src);
        },
        48 => {
            arg.iter_mut()
                .zip(lookup1.iter().map(|x| phi_bls48(x,4)))
                .for_each(|(dest, src)| *dest = src);
        },
        _ => panic!("Invalid degree value")
    }

    arg
}

pub fn gls8_multiply_bls24<const PRAMASIZE:usize, const R: usize, const N: usize, const MAX_COEFS_COUNT: usize>
            (input :&G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>, scalar :&FieldElement<R>) -> G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>
{   
    //  Constant-Time multiplication for elements in G2 (8-GLS): GLS Implementation of points multiplication on G2, for the BLS24 Curve
    //  Joppe W. Bos, Craig Costello, and Michael Naehrig https://eprint.iacr.org/2013/458.pdf    
    let ff: BigUint = BigUint::from_str("255").unwrap();
    let mut code = recod_scalar_gls8(&scalar.to_big_uint(),
                                               input.consts.u.abs() as u128,
                                               &BigUint::from_str_radix(&input.consts.lambda.fieldparams.modulo_as_strhex[2..],16).unwrap());
    let infinit = create_infinit(input);
    let mut lookup1: [G2Element<PRAMASIZE, R, N, MAX_COEFS_COUNT>; 8] = [infinit;8];

    populate_lookup1(input, &mut lookup1);

    let lookup2  = create_lookup(infinit, &lookup1, 24);

    let mut limb : u8 = (&code).bitand(ff.clone()).to_u8().unwrap();
    let out_sig : i8 =  1 - ((limb & 1) << 1) as i8;    
    code = code >> 1;
    limb = (&code).bitand(ff.clone()).to_u8().unwrap();
    let (mut c1, mut c2) = ((limb & 15) as usize , ((limb >> 4) & 15) as usize);
    let (mut sign1 , mut sign2) = (((c1 & 1) << 1) as i8 - 1, ((c2 & 1) << 1) as i8 - 1);
    let mut result = compute_res(sign1, &lookup1, c1 >> 1).
        addto(&compute_res(sign2, &lookup2, c2 >> 1));
    code = code >> 8;
    while code != BigUint::one() {
        result = result.double();
        limb = (&code).bitand(ff.clone()).to_u8().unwrap();
        (c1, c2) = ((limb & 15) as usize , ((limb >> 4) & 15) as usize);
        (sign1, sign2) = (((c1 & 1) << 1) as i8 - 1, ((c2 & 1) << 1) as i8 - 1);
        result = result.
            addto(&compute_res(sign1, &lookup1, c1 >> 1)).
            addto(&compute_res(sign2, &lookup2, c2 >> 1));
        code = code >> 8;
    }
    if out_sig == 1 {result} else {result.negate()}
}

fn get_sign_tuple(c1: &usize, c2: &usize, c3: &usize, c4: &usize) -> (i8, i8, i8, i8) {
    (((c1 & 1) << 1) as i8 - 1, ((c2 & 1) << 1) as i8 - 1, ((c3 & 1) << 1) as i8 - 1, ((c4 & 1) << 1) as i8 - 1)
}

pub fn gls16_multiply_bls48<const PRAMASIZE:usize, const R: usize, const N: usize, const MAX_COEFS_COUNT: usize>
            (input :&G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>, scalar :&FieldElement<R>) -> G2Element<PRAMASIZE, R,N,MAX_COEFS_COUNT>
{   
    //  Constant-Time multiplication for elements in G2 (16-GLS): GLS Implementation of points multiplication on G2 for the BLS48 Curve    //  Inspired from Joppe W. Bos, Craig Costello, and Michael Naehrig https://eprint.iacr.org/2013/458.pdf
    let ffff: BigUint = BigUint::from_str("65535").unwrap();
    let mut code = recod_scalar_gls16(&scalar.to_big_uint(),
                                               input.consts.u.abs() as u64,
                                               &BigUint::from_str_radix(&input.consts.lambda.fieldparams.modulo_as_strhex[2..],16).unwrap());

    let infinit = create_infinit(input);

    let mut lookup1: [G2Element<PRAMASIZE, R, N, MAX_COEFS_COUNT>; 8] = [infinit;8];

    populate_lookup1(input, &mut lookup1);

    let lookup2 = create_lookup(infinit, &lookup1, 48);
    let lookup3 = create_lookup(infinit, &lookup2, 48);
    let lookup4 = create_lookup(infinit, &lookup3, 48);

    let mut limb : u16 = (&code).bitand(ffff.clone()).to_u16().unwrap();
    let out_sig : i8 =  1 - ((limb & 1) << 1) as i8;    
    code = code >> 1;
    limb = (&code).bitand(ffff.clone()).to_u16().unwrap();
    let (mut c1, mut c2,mut c3,mut c4) = ((limb & 15) as usize , ((limb >> 4) & 15) as usize,
                                                                      ((limb >> 8) & 15) as usize,((limb >> 12) & 15) as usize);
    let (mut sign1 , mut sign2, mut sign3, mut sign4) = get_sign_tuple(&c1, &c2, &c3, &c4);
    let mut result1 = compute_res(sign1, &lookup1, c1 >> 1).
        addto(&compute_res(sign2, &lookup2, c2 >> 1));
    let mut result2 = compute_res(sign3, &lookup3, c3 >> 1).
        addto(&compute_res(sign4, &lookup4, c4 >> 1));
    let mut result = result1.addto(&result2);
    code = code >> 16;
    while code != BigUint::one() {
        result = result.double();
        limb = (&code).bitand(ffff.clone()).to_u16().unwrap();
        (c1,c2,c3,c4) = ((limb & 15) as usize , ((limb >> 4) & 15) as usize,((limb >> 8) & 15) as usize,((limb >> 12) & 15) as usize);
        (sign1,sign2,sign3,sign4) = get_sign_tuple(&c1, &c2, &c3, &c4);
        result1 = compute_res(sign1, &lookup1, c1 >> 1).addto(&compute_res(sign2, &lookup2, c2 >> 1));
        result2 = compute_res(sign3, &lookup3, c3 >> 1).addto(&compute_res(sign4, &lookup4, c4 >> 1));
        result = result.addto(&result1.addto(&result2));
        code = code >> 16;
    }
    if out_sig == 1 {result} else {result.negate()}
}