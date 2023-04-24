// Copyright (C) myl7
// SPDX-License-Identifier: Apache-2.0

//! Distributed multi-point function (DMPF) implementation

use std::mem;

use num_bigint::{BigUint, ToBigUint};
use num_integer::Integer;
use statrs::distribution::{ContinuousCDF, Normal};

use crate::dpf::{BSampler, Gen, PointFn, Share, VDPF};

/// `VerDMPF` in the paper.
/// $\mathcal(k)$ is fixed to 3 since it affects the form of `ch_compact`.
pub struct VDMPF {
    /// $\lambda$ in the Cuckoo-hashing part of the paper
    lambda_p: f64,
    /// Cuckoo-hashing retry limit.
    /// Only collisions are counted, and successful insertion is not.
    ch_retry: usize,
    /// PRP (pseudo-random permutation), $\{0, 1\} \times [N\mathcal{k}] \rightarrow [N\mathcal{k}}]$ in the paper
    prp: Box<dyn Permu>,
    /// Used to sample hash function indexes in Cuckoo-hashing.
    /// Since $\mathcal{k}$ is fixed to 3, the sampling range is always `[0, 2]`.
    ch_sampler: Box<dyn ISampler>,
    /// To sample the PRP seed
    prp_sampler: Box<dyn BSampler>,
    /// Retry limit for PRP seed resampling
    prp_retry: usize,
    lambda: usize,
    vdpf: VDPF,
}

impl VDMPF {
    pub fn new(
        lambda_p: f64,
        ch_retry: usize,
        prp: Box<dyn Permu>,
        ch_sampler: Box<dyn ISampler>,
        prp_sampler: Box<dyn BSampler>,
        prp_retry: usize,
        lambda: usize,
        prg: Box<dyn Gen>,
        hash: Box<dyn Gen>,
        hash_prime: Box<dyn Gen>,
        dpf_sampler: Box<dyn BSampler>,
        gen_retry: usize,
    ) -> Self {
        Self {
            lambda_p,
            ch_retry,
            prp,
            ch_sampler,
            prp_sampler,
            prp_retry,
            lambda,
            vdpf: VDPF::new(lambda, prg, hash, hash_prime, dpf_sampler, gen_retry),
        }
    }
}

/// `VerDMPF` API
impl VDMPF {
    const K: usize = 3;

    /// `Gen` in the paper.
    /// PGP seed `seed` should be randomly sampled.
    pub fn gen(&self, fs: &[&PointFn]) -> Result<MShare, ()> {
        let mut seed;
        // + 2 to use 0 as the boundary
        let mut prp_retry = self.prp_retry + 2;
        let (table, indexes) = loop {
            prp_retry -= 1;
            if prp_retry == 0 {
                return Err(());
            }
            seed = self.prp_sampler.sample(self.lambda);
            // Use `usize` other than `f64` because all of them will be used as (part of) indexes.
            let t = fs.len();
            let m = self.ch_bucket(t);
            let n = 2.to_biguint().unwrap().pow(self.lambda as u32);
            let b = (n.clone() * VDMPF::K.to_biguint().unwrap()).div_ceil(&m.to_biguint().unwrap());
            let yb_gen = |i: usize, x: &[u8]| {
                let xb =
                    BigUint::from_bytes_be(x) + 0.to_biguint().unwrap() * i.to_biguint().unwrap();
                let mut xbb = xb.to_bytes_be();
                self.prp.permu(&seed, &mut xbb);
                BigUint::from_bytes_be(&xbb)
            };
            let hs = (0..VDMPF::K)
                .map(|i| {
                    let b = b.clone();
                    move |x: &[u8]| {
                        let yb = yb_gen(i, x);
                        (yb / &b).to_u64_digits()[0] as usize
                    }
                })
                .collect::<Vec<_>>();
            let indexes = (0..VDMPF::K)
                .map(|i| {
                    let b = b.clone();
                    move |x: &[u8]| {
                        let yb = yb_gen(i, x);
                        (yb % &b).to_bytes_be()
                    }
                })
                .collect::<Vec<_>>();
            let a_s = fs.iter().map(|f| f.a.as_ref()).collect::<Vec<_>>();
            let table = match self.ch_compact(&a_s, &hs, m) {
                Ok(table) => table,
                Err(_) => continue,
            };
            break (table, indexes);
        };
        let mut mshare = MShare {
            ks: vec![],
            seed: seed.clone(),
        };
        for item in table.into_iter() {
            let (a, b) = match item {
                // `aj` is $a_j$ in the paper
                Some((aji, k)) => {
                    let a = indexes[k](fs[aji].a.as_ref());
                    let b = fs[aji].b.clone();
                    (a, b)
                }
                None => (vec![], vec![]),
            };
            let f = PointFn { a, b };
            let share = self.vdpf.gen(f)?;
            mshare.ks.push(share);
        }
        Ok(mshare)
    }
}

impl VDMPF {
    /// `CHBucket` in the paper.
    /// $\mathcal{k}$ is fixed to 3.
    fn ch_bucket(&self, t: usize) -> usize {
        let tf = t as f64;
        let norm_a = Normal::new(6.3, 2.3).unwrap();
        let norm_b = Normal::new(6.45, 2.18).unwrap();
        // Affected by $\mathcal{k}$
        let m = tf * (self.lambda_p + norm_b.cdf(tf) + tf.log2() / norm_a.cdf(tf));
        m.ceil() as usize
    }

    /// `CHCompact` in the paper
    fn ch_compact(
        &self,
        a_s: &[&[u8]],
        hs: &[impl Fn(&[u8]) -> usize],
        m: usize,
    ) -> Result<Vec<Option<(usize, usize)>>, ()> {
        // `Table` in the paper
        let mut table = vec![None; m];
        // + 2 to use 0 as the boundary
        let mut retry = self.ch_retry + a_s.len() + 2;
        for ai in 0..a_s.len() {
            let mut pending = Some((ai, 0));
            while let Some((pi, _)) = pending {
                retry -= 1;
                if retry == 0 {
                    return Err(());
                }
                let k = self.ch_sampler.sample(VDMPF::K);
                pending.iter_mut().for_each(|(_, pk)| *pk = k);
                let i = hs[k](a_s[pi]);
                mem::swap(&mut table[i], &mut pending);
            }
        }
        Ok(table)
    }
}

/// `k` in the paper.
/// `ks` is like `s0s` of [`crate::dpf::Share`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MShare {
    pub ks: Vec<Share>,
    /// $\delta$ in the paper
    pub seed: Vec<u8>,
}

/// Interface for PRP.
/// See [`VDMPF`].
pub trait Permu {
    // Permutation in-place
    fn permu(&self, seed: &[u8], x: &mut [u8]);
}

/// Interface for sampling indexes in Cuckoo-hashing.
/// See [`VDMPF`].
pub trait ISampler {
    /// Sample from `[0, n)`.
    fn sample(&self, n: usize) -> usize;
}

#[cfg(test)]
pub(crate) mod tests_fixture {
    use super::*;

    use aes::cipher::generic_array::GenericArray;
    use aes::cipher::{BlockEncrypt, KeyInit};
    use aes::Aes256;
    use rand::prelude::*;
    use rand_chacha::ChaChaRng;
    use rand_seeder::Seeder;

    #[derive(Default)]
    struct PRP {}

    impl Permu for PRP {
        fn permu(&self, seed: &[u8], x: &mut [u8]) {
            assert_eq!(x.len(), 16);
            let mut rng: ChaChaRng = Seeder::from(seed).make_rng();
            let key = GenericArray::from(rng.gen::<[u8; 32]>());
            let mut block = GenericArray::from_mut_slice(&mut x[0..16]);
            let cipher = Aes256::new(&key);
            cipher.encrypt_block(&mut block);
        }
    }

    #[derive(Default)]
    struct CHSampler {}

    impl ISampler for CHSampler {
        fn sample(&self, n: usize) -> usize {
            let mut rng = ChaChaRng::from_rng(thread_rng()).unwrap();
            rng.gen_range(0..n - 1)
        }
    }
}

#[cfg(test)]
mod tests {
    // use super::tests_fixture::*;
    // use super::*;
    // use crate::dpf::tests_fixture::*;

    // #[test]
    // fn test_gen_eval_verify_ok() {
    //     todo!("Impl")
    // }
}
