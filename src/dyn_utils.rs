use std::cell::RefCell;

use aes::cipher::generic_array::GenericArray;
use aes::cipher::{BlockEncrypt, KeyInit};
use aes::Aes256;
use rand::prelude::*;
use rand_chacha::ChaChaRng;
use rand_pcg::Pcg64Mcg;
use rand_seeder::Seeder;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

use crate::dmpf::{ISampler, Permu};
use crate::dpf::{BSampler, Gen};

#[derive(Default)]
pub struct Aes256PRP {}

impl Permu for Aes256PRP {
    fn domains(&self) -> (u32, u32) {
        (126, 128)
    }

    fn permu(&self, seed: &[u8], x: &mut [u8]) {
        assert_eq!(x.len(), 16);
        let mut rng: ChaChaRng = Seeder::from(seed).make_rng();
        let key = GenericArray::from(rng.gen::<[u8; 32]>());
        let mut block = GenericArray::from_mut_slice(&mut x[0..16]);
        let cipher = Aes256::new(&key);
        cipher.encrypt_block(&mut block);
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
pub struct Shake256Hash {}

impl Gen for Shake256Hash {
    fn gen(&self, input: &[u8], output_len: usize) -> Vec<u8> {
        let mut hasher = Shake256::default();
        hasher.update(input);
        let mut reader = hasher.finalize_xof();
        let mut output = vec![0u8; output_len];
        reader.read(&mut output);
        output
    }
}

#[derive(Default)]
pub struct ChaChaPRG {}

impl Gen for ChaChaPRG {
    fn gen(&self, input: &[u8], output_len: usize) -> Vec<u8> {
        let mut rng: ChaChaRng = Seeder::from(input).make_rng();
        let mut output = vec![0u8; output_len];
        rng.fill_bytes(&mut output);
        output
    }
}

#[derive(Default)]
pub struct Pcg64McgPRG {}

impl Gen for Pcg64McgPRG {
    fn gen(&self, input: &[u8], output_len: usize) -> Vec<u8> {
        let mut rng: Pcg64Mcg = Seeder::from(input).make_rng();
        let mut output = vec![0u8; output_len];
        rng.fill_bytes(&mut output);
        output
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
