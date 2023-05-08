use std::cell::RefCell;

use aes::cipher::generic_array::GenericArray;
use aes::cipher::{BlockEncrypt, KeyInit};
use aes::{Aes128, Aes256};
use hex_literal::hex;
use rand::prelude::*;
use rand_chacha::ChaChaRng;
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
        let mut seed = [0u8; 32];
        seed[0..16].copy_from_slice(input);
        seed[16..32].copy_from_slice(input);
        let mut rng = ChaChaRng::from_seed(seed);
        // let mut rng: ChaChaRng = Seeder::from(input).make_rng();
        let mut output = vec![0u8; output_len];
        rng.fill_bytes(&mut output);
        output
    }
}

#[derive(Default)]
pub struct Aes128PRG {}

impl Aes128PRG {
    const BLOCKS: [u8; 48] = hex!("445299ccf9df04b5fad62dd09271ccf55d068d6471f97b590c96de02b6399dd7336f0e28cc8e5766b66cc411b5975e68");
}

impl Gen for Aes128PRG {
    fn gen(&self, input: &[u8], output_len: usize) -> Vec<u8> {
        assert_eq!(input.len(), 16);
        assert!(output_len <= 16 * 3);
        let mut output = vec![0u8; 16 * 3];
        let mut key = GenericArray::clone_from_slice(input);
        for i in 0..3 {
            output[i * 16..(i + 1) * 16].copy_from_slice(&Self::BLOCKS[i * 16..(i + 1) * 16]);
            let mut block = GenericArray::from_mut_slice(&mut output[i * 16..(i + 1) * 16]);
            let cipher = Aes128::new(&key);
            cipher.encrypt_block(&mut block);
            for j in 0..16 {
                output[i * 16 + j] ^= Self::BLOCKS[i * 16 + j];
            }
            key = GenericArray::clone_from_slice(&output[i * 16..(i + 1) * 16]);
        }
        output.truncate(output_len);
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
