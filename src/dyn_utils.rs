use std::cell::RefCell;

use aes::cipher::generic_array::GenericArray;
use aes::cipher::{BlockEncrypt, KeyInit};
use aes::Aes128;
use bitvec::prelude::*;
use hex_literal::hex;
use lazy_static::lazy_static;
use rand::prelude::*;
use rand_chacha::ChaChaRng;

// use crate::dmpf::{ISampler, Permu};
use crate::dmpf::{ISampler, Permu};
use crate::dpf::{BSampler, Hash, HashPrime, PRG};

#[derive(Default)]
pub struct Aes128PRP {}

impl Permu for Aes128PRP {
    fn domains(&self) -> (u32, u32) {
        (126, 128)
    }

    fn permu(&self, seed: &[u8], x: &mut [u8]) {
        assert_eq!(seed.len(), 16);
        assert_eq!(x.len(), 16);
        let key = GenericArray::from_slice(seed);
        let block = GenericArray::from_mut_slice(x);
        let cipher = Aes128::new(key);
        cipher.encrypt_block(block);
    }
}

pub struct ChaChaISampler {
    rng: RefCell<ChaChaRng>,
}

impl Default for ChaChaISampler {
    fn default() -> Self {
        Self {
            rng: RefCell::new(ChaChaRng::seed_from_u64(thread_rng().gen())),
        }
    }
}

impl ChaChaISampler {
    pub fn new(seed: u64) -> Self {
        Self {
            rng: RefCell::new(ChaChaRng::seed_from_u64(seed)),
        }
    }
}

impl ISampler for ChaChaISampler {
    fn sample(&self, n: usize) -> usize {
        self.rng.borrow_mut().gen_range(0..n - 1)
    }
}

lazy_static! {
    static ref AES128_CIPHERS: [Aes128; 4] = [
        Aes128::new(&GenericArray::from(hex!(
            "445299ccf9df04b5fad62dd09271ccf5"
        ))),
        Aes128::new(&GenericArray::from(hex!(
            "d83473e081a00060dbb57a6257122581"
        ))),
        Aes128::new(&GenericArray::from(hex!(
            "d068d6471f97b590c96de02b6399dd73"
        ))),
        Aes128::new(&GenericArray::from(hex!(
            "36f0e28cc8e5766b66cc411b5975e68f"
        )))
    ];
}

#[derive(Default)]
pub struct Aes128Hash {}

impl Hash for Aes128Hash {
    fn gen(&self, x0: &mut [u8], x1: &[u8]) {
        assert!(x0.len() == 64 && x1.len() == 16);
        let mut buf = [0; 16];
        for i in 0..16 {
            x0[16 + i] = x1[i];
        }

        buf.copy_from_slice(&x0[0..16]);
        let mut block = GenericArray::from_mut_slice(&mut buf);
        AES128_CIPHERS[0].encrypt_block(&mut block);
        for i in 0..16 {
            x0[i] ^= block[i];
        }

        buf.copy_from_slice(&x0[16..32]);
        let mut block = GenericArray::from_mut_slice(&mut buf);
        AES128_CIPHERS[1].encrypt_block(&mut block);
        for i in 0..16 {
            x0[16 + i] ^= block[i];
        }

        buf.copy_from_slice(&x1);
        let mut block = GenericArray::from_mut_slice(&mut buf);
        AES128_CIPHERS[0].encrypt_block(&mut block);
        for i in 0..16 {
            x0[32 + i] = x1[i] ^ block[i];
        }

        buf.copy_from_slice(&x1);
        let mut block = GenericArray::from_mut_slice(&mut buf);
        AES128_CIPHERS[1].encrypt_block(&mut block);
        for i in 0..16 {
            x0[48 + i] = x1[i] ^ block[i];
        }
    }
}

impl HashPrime for Aes128Hash {
    fn gen(&self, x: &mut [u8]) {
        assert_eq!(x.len(), 64);
        let mut buf = [0; 16];
        let mut block = GenericArray::from_mut_slice(&mut buf);

        for i in 0..4 {
            let block_in = GenericArray::from_slice(&x[i * 16..(i + 1) * 16]);
            AES128_CIPHERS[i].encrypt_block_b2b(block_in, &mut block);
            for j in 0..16 {
                x[i * 16 + j] ^= block[j]
            }
        }

        for i in 0..32 {
            x[i] ^= x[i + 32];
            x[i + 32] = 0;
        }
    }
}

#[derive(Default)]
pub struct Aes128PRG {}

impl PRG for Aes128PRG {
    fn gen(&self, x: &mut [u8], x2: &mut [u8]) -> [bool; 2] {
        assert_eq!(x.len(), 16);
        assert!(x.len() == x2.len());

        x2.copy_from_slice(x);
        let block = GenericArray::from_mut_slice(x2);
        AES128_CIPHERS[0].encrypt_block(block);
        for i in 0..16 {
            x2[i] ^= x[i];
        }
        let ts = [x2[0].view_bits::<Msb0>()[0], x2[0].view_bits::<Msb0>()[1]];

        x2.copy_from_slice(x);
        let block = GenericArray::from_mut_slice(x2);
        AES128_CIPHERS[2].encrypt_block(block);
        for i in 0..16 {
            x2[i] ^= x[i];
        }

        let mut buf = [0; 16];
        buf.copy_from_slice(x);
        let block = GenericArray::from_mut_slice(&mut buf);
        AES128_CIPHERS[1].encrypt_block(block);
        for i in 0..16 {
            x[i] ^= buf[i];
        }

        ts
    }
}

pub struct ChaChaBSampler {
    pub rng: RefCell<ChaChaRng>,
}

impl Default for ChaChaBSampler {
    fn default() -> Self {
        Self {
            rng: RefCell::new(ChaChaRng::seed_from_u64(thread_rng().gen())),
        }
    }
}

impl ChaChaBSampler {
    pub fn new(seed: u64) -> Self {
        Self {
            rng: RefCell::new(ChaChaRng::seed_from_u64(seed)),
        }
    }
}

impl BSampler for ChaChaBSampler {
    fn sample(&self, len: usize) -> Vec<u8> {
        let mut buf = vec![0u8; len];
        self.rng.borrow_mut().fill_bytes(&mut buf);
        buf
    }
}
