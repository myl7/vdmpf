#![feature(test)]

extern crate test;

use test::Bencher;

use std::cell::RefCell;

use aes::cipher::generic_array::GenericArray;
use aes::cipher::{BlockEncrypt, KeyInit};
use aes::Aes256;
use hex_literal::hex;
use rand::prelude::*;
use rand_chacha::ChaChaRng;
use rand_seeder::Seeder;
use sha3::digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;

use vdmpf::dmpf::{ISampler, Permu, VDMPF};
use vdmpf::dpf::{BSampler, Gen, PointFn};
use vdmpf::group::BGroup;

#[derive(Default)]
pub struct PRP {}

impl Permu for PRP {
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

pub struct CHSampler {
    rng: RefCell<ChaChaRng>,
}

impl Default for CHSampler {
    fn default() -> Self {
        Self {
            rng: RefCell::new(ChaChaRng::seed_from_u64(thread_rng().gen())),
        }
    }
}

impl ISampler for CHSampler {
    fn sample(&self, n: usize) -> usize {
        self.rng.borrow_mut().gen_range(0..n - 1)
    }
}

#[derive(Default)]
pub struct Hash {}

impl Gen for Hash {
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
pub struct PRG {}

impl Gen for PRG {
    fn gen(&self, input: &[u8], output_len: usize) -> Vec<u8> {
        let mut rng: ChaChaRng = Seeder::from(input).make_rng();
        let mut output = vec![0u8; output_len];
        rng.fill_bytes(&mut output);
        output
    }
}

pub struct Sampler {
    pub rng: RefCell<ChaChaRng>,
}

impl Default for Sampler {
    fn default() -> Self {
        Self {
            rng: RefCell::new(ChaChaRng::seed_from_u64(thread_rng().gen())),
        }
    }
}

impl BSampler for Sampler {
    fn sample(&self, len: usize) -> Vec<u8> {
        let mut buf = vec![0u8; len];
        self.rng.borrow_mut().fill_bytes(&mut buf);
        buf
    }
}

#[bench]
fn main(b: &mut Bencher) {
    b.iter(|| {
        let vdmpf = VDMPF::new(
            80f64,
            1000,
            Box::new(PRP::default()),
            Box::new(CHSampler::default()),
            Box::new(Sampler::default()),
            1000,
            16,
            Box::new(PRG::default()),
            Box::new(Hash::default()),
            Box::new(Hash::default()),
            Box::new(Hash::default()),
            Box::new(Sampler::default()),
            1000,
        );
        let mut fs = Vec::with_capacity(100);
        let mut f_rng = thread_rng();
        for _ in 0..100 {
            let mut a = vec![0; 16];
            f_rng.fill_bytes(&mut a);
            a[0] = 0x07;
            let mut b = vec![0; 16];
            f_rng.fill_bytes(&mut b);
            fs.push(PointFn { a, b });
        }
        let fs_ref = fs.iter().collect::<Vec<_>>();
        let mut mshare0 = vdmpf.gen(&fs_ref).unwrap();
        let mut mshare1 = mshare0.clone();
        mshare0
            .ks
            .iter_mut()
            .for_each(|k| k.s0s = vec![k.s0s[0].clone()]);
        mshare1
            .ks
            .iter_mut()
            .for_each(|k| k.s0s = vec![k.s0s[1].clone()]);

        let mut xs = fs.iter().map(|f| f.a.clone()).collect::<Vec<_>>();
        for _ in 0..100 {
            let mut a = vec![0; 16];
            f_rng.fill_bytes(&mut a);
            a[0] = 0x0f;
            xs.push(a);
        }
        let xs_ref = xs.iter().map(|x| x.as_ref()).collect::<Vec<_>>();
        let (y0s, pi0) = vdmpf.eval(false, &mshare0, &xs_ref, fs.len());
        let (y1s, pi1) = vdmpf.eval(true, &mshare1, &xs_ref, fs.len());
        assert_eq!(vdmpf.verify(&[&pi0, &pi1]), true);
        for ((x, y0), y1) in xs.iter().zip(y0s.iter()).zip(y1s.iter()) {
            let y: Vec<u8> = (BGroup::from(y0.to_owned()) + y1.as_ref()).into();
            match fs.iter().find(|f| f.a == *x) {
                Some(f) => assert_eq!(y, f.b),
                None => assert_eq!(y, hex!("00000000000000000000000000000000")),
            }
        }
    });
}
