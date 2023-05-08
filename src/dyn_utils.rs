use std::cell::RefCell;

use aes::cipher::generic_array::GenericArray;
use aes::cipher::{BlockEncrypt, KeyInit};
use aes::Aes128;
use bitvec::prelude::*;
use hex_literal::hex;
use rand::prelude::*;
use rand_chacha::ChaChaRng;

// use crate::dmpf::{ISampler, Permu};
use crate::dpf::{BSampler, Hash, HashPrime, PRG};
use crate::dmpf::{ISampler, Permu};

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

#[derive(Default)]
pub struct Aes128Hash {}

impl Aes128Hash {
    const BLOCKS: [u8; 64] = hex!("445299ccf9df04b5fad62dd09271ccf55d068d6471f97b590c96de02b6399dd7336f0e28cc8e5766b66cc411b5975e68336f0e28cc8e5766b66cc411b5975e68");
}

impl Hash for Aes128Hash {
    fn gen(&self, x0: &mut [u8], x1: &[u8]) {
        assert!(x0.len() == 64 && x1.len() == 16, "x0.len() == {} && x1.len() == {}", x0.len(), x1.len());
        for i in 0..4 {
            let key = GenericArray::clone_from_slice(if i == 0{
                x1
            } else {
                &x0[16 * (i - 1)..16 * i]
            });
            x0[16 * i..16 * (i + 1)].copy_from_slice(&Self::BLOCKS[16 * i..16 * (i + 1)]);
            let block = GenericArray::from_mut_slice(&mut x0[16 * i..16 * (i + 1)]);
            let cipher = Aes128::new(&key);
            cipher.encrypt_block(block);
            for j in 0..16 {
                x0[16 * i + j] ^= Self::BLOCKS[16 * i + j];
            }
        }
    }
}

impl HashPrime for Aes128Hash {
    fn gen(&self, x: &mut [u8]) {
        todo!()
    }
}

#[derive(Default)]
pub struct Aes128PRG {}

impl Aes128PRG {
    const BLOCKS: [u8; 48] = hex!("445299ccf9df04b5fad62dd09271ccf55d068d6471f97b590c96de02b6399dd7336f0e28cc8e5766b66cc411b5975e68");
}

impl PRG for Aes128PRG {
    fn gen(&self, x: &mut [u8], x2: &mut [u8]) -> [bool; 2] {
        assert_eq!(x.len(), 16);
        assert!(x.len() == x2.len());

        let key = GenericArray::from_slice(x);
        x2.copy_from_slice(&Self::BLOCKS[0..16]);
        let block = GenericArray::from_mut_slice(x2);
        let cipher = Aes128::new(key);
        cipher.encrypt_block(block);
        for i in 0..16 {
            x2[i] ^= Self::BLOCKS[i];
        }
        let ts = [x2[0].view_bits::<Msb0>()[0], x2[0].view_bits::<Msb0>()[1]];

        let key = GenericArray::from_slice(x2);
        x.copy_from_slice(&Self::BLOCKS[16..32]);
        let block = GenericArray::from_mut_slice(x);
        let cipher = Aes128::new(key);
        cipher.encrypt_block(block);
        for i in 0..16 {
            x[i] ^= Self::BLOCKS[16 + i];
        }

        let key = GenericArray::from_slice(x);
        x2.copy_from_slice(&Self::BLOCKS[32..48]);
        let block = GenericArray::from_mut_slice(x2);
        let cipher = Aes128::new(key);
        cipher.encrypt_block(block);
        for i in 0..16 {
            x2[i] ^= Self::BLOCKS[32 + i];
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
