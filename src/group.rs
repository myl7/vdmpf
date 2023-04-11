// Copyright (C) myl7
// SPDX-License-Identifier: Apache-2.0

use std::ops::{Add, AddAssign};

/// Bytes with xor operations as a mathmatical group
#[derive(Clone)]
pub struct BGroup(BGroupInner);

#[derive(Clone)]
enum BGroupInner {
    Bytes(Vec<u8>),
    Bit(bool),
}

impl Into<Vec<u8>> for BGroup {
    fn into(self) -> Vec<u8> {
        match self.0 {
            BGroupInner::Bytes(b) => b,
            BGroupInner::Bit(_) => panic!("BGroup contains bytes"),
        }
    }
}

impl Into<bool> for BGroup {
    fn into(self) -> bool {
        match self.0 {
            BGroupInner::Bytes(_) => panic!("BGroup contains a bit"),
            BGroupInner::Bit(b) => b,
        }
    }
}

impl From<Vec<u8>> for BGroup {
    fn from(v: Vec<u8>) -> Self {
        Self(BGroupInner::Bytes(v))
    }
}

impl From<bool> for BGroup {
    fn from(b: bool) -> Self {
        Self(BGroupInner::Bit(b))
    }
}

impl AddAssign<&[u8]> for BGroup {
    fn add_assign(&mut self, rhs: &[u8]) {
        match &mut self.0 {
            BGroupInner::Bytes(b) => {
                b.iter_mut().zip(rhs.iter()).for_each(|(a, b)| *a ^= b);
            }
            BGroupInner::Bit(_) => panic!("BGroup inner type mismatching"),
        }
    }
}

impl AddAssign<bool> for BGroup {
    fn add_assign(&mut self, rhs: bool) {
        fn xor_bool(a: bool, b: bool) -> bool {
            if a == b {
                false
            } else {
                true
            }
        }

        match &mut self.0 {
            BGroupInner::Bytes(_) => panic!("BGroup inner type mismatching"),
            BGroupInner::Bit(b) => *b = xor_bool(*b, rhs),
        }
    }
}

impl AddAssign for BGroup {
    fn add_assign(&mut self, rhs: Self) {
        match rhs.0 {
            BGroupInner::Bytes(b) => *self += b.as_ref(),
            BGroupInner::Bit(b) => *self += b,
        }
    }
}

impl Add<&[u8]> for BGroup {
    type Output = Self;

    fn add(mut self, rhs: &[u8]) -> Self::Output {
        self += rhs;
        self
    }
}

impl Add<bool> for BGroup {
    type Output = Self;

    fn add(mut self, rhs: bool) -> Self::Output {
        self += rhs;
        self
    }
}

impl Add for BGroup {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}
