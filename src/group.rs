// Copyright (C) myl7
// SPDX-License-Identifier: Apache-2.0

pub fn xor_bytes(a: &mut [u8], b: &[u8]) {
    assert!(a.len() == b.len());
    a.iter_mut().zip(b.iter()).for_each(|(a, b)| *a ^= b);
}

pub fn xor_bool(a: bool, b: bool) -> bool {
    a != b
}
