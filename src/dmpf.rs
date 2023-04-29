// Copyright (C) myl7
// SPDX-License-Identifier: Apache-2.0

//! Distributed multi-point function (DMPF) implementation

use std::collections::HashMap;
use std::mem;
use std::rc::Rc;

use num_bigint::{BigUint, ToBigUint};
use num_integer::Integer;
use statrs::distribution::{ContinuousCDF, Normal};

use crate::dpf::{BSampler, Gen, PointFn, Share, VDPF};
use crate::group::BGroup;

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
    hash_prime: Box<dyn Gen>,
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
        hash_prime2: Box<dyn Gen>,
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
            hash_prime: hash_prime2,
        }
    }
}

/// `VerDMPF` API
impl VDMPF {
    const K: usize = 3;

    /// `Gen` in the paper
    pub fn gen(&self, fs: &[&PointFn]) -> Result<MShare, ()> {
        // Ensure $\alpha < 2^126$ so that $n\mathcal{k} < 2^128$ and we can use AES as PRP
        fs.iter()
            .for_each(|f| assert!(BigUint::from_bytes_be(&f.a) < 2.to_biguint().unwrap().pow(126)));

        // + 2 to use 0 as the boundary
        let mut prp_retry = self.prp_retry + 2;
        let prp_rc = Rc::new(&self.prp);
        let n = 2.to_biguint().unwrap().pow(126);
        let (table, seed, n, b, n_prime) = loop {
            prp_retry -= 1;
            if prp_retry == 0 {
                return Err(());
            }
            let seed = self.prp_sampler.sample(self.lambda);
            // Use `usize` other than `f64` because all of them will be used as (part of) indexes.
            let t = fs.len();
            let m = self.ch_bucket(t);
            let b = 2
                .to_biguint()
                .unwrap()
                .pow(128)
                .div_ceil(&m.to_biguint().unwrap());
            // $n' / 8$ in the paper
            let n_prime = (b.bits() as f64 / 8 as f64).ceil() as usize;
            let hs = (0..VDMPF::K)
                .map(|i| {
                    let n = n.clone();
                    let prp = prp_rc.clone();
                    let b = b.clone();
                    let seed = seed.clone();
                    move |x: &[u8]| {
                        let xb = BigUint::from_bytes_be(x) + &n * i.to_biguint().unwrap();
                        let mut xbb = xb.to_bytes_be();
                        prp.permu(&seed, &mut xbb);
                        let yb = BigUint::from_bytes_be(&xbb);
                        (yb / &b).to_u64_digits()[0] as usize
                    }
                })
                .collect::<Vec<_>>();
            let a_s = fs.iter().map(|f| f.a.as_ref()).collect::<Vec<_>>();
            let table = match self.ch_compact(&a_s, &hs, m) {
                Ok(table) => table,
                Err(_) => continue,
            };
            break (table, seed, n, b, n_prime);
        };
        let mut mshare = MShare {
            ks: vec![],
            seed: seed.clone(),
        };
        for item in table.into_iter() {
            let (a, b) = match item {
                // `aj` is $a_j$ in the paper
                Some((aji, k)) => {
                    let xb = BigUint::from_bytes_be(fs[aji].a.as_ref())
                        + n.clone() * k.to_biguint().unwrap();
                    let mut xbb = xb.to_bytes_be();
                    prp_rc.permu(&seed, &mut xbb);
                    let yb = BigUint::from_bytes_be(&xbb);
                    let mut index = (yb % &b).to_bytes_be();
                    for _ in 0..(n_prime - index.len()) {
                        index.insert(0, 0);
                    }
                    let a = index;
                    let b = fs[aji].b.clone();
                    (a, b)
                }
                None => (vec![0; n_prime], vec![0; self.lambda]),
            };
            let f = PointFn { a, b };
            let share = self.vdpf.gen(f)?;
            mshare.ks.push(share);
        }
        Ok(mshare)
    }

    /// `BVEval` in the paper
    pub fn eval(
        &self,
        b_party: bool,
        mshare: &MShare,
        xs: &[&[u8]],
        t: usize,
    ) -> (Vec<Vec<u8>>, Vec<u8>) {
        assert_eq!(mshare.ks.first().map(|k| k.s0s.len()), Some(1));
        let n = 2.to_biguint().unwrap().pow(126);
        let m = self.ch_bucket(t);
        let b = 2
            .to_biguint()
            .unwrap()
            .pow(128)
            .div_ceil(&m.to_biguint().unwrap());
        let n_prime = (b.bits() as f64 / 8 as f64).ceil() as usize;
        let mut inputs = vec![vec![]; m];
        let mut dedup = HashMap::new();
        xs.iter().enumerate().for_each(|(eta, x)| {
            let is = (0..VDMPF::K)
                .map(|i| {
                    let xb = BigUint::from_bytes_be(x) + &n * i.to_biguint().unwrap();
                    let mut xbb = xb.to_bytes_be();
                    self.prp.permu(&mshare.seed, &mut xbb);
                    let yb = BigUint::from_bytes_be(&xbb);
                    (yb / &b).to_u64_digits()[0] as usize
                })
                .collect::<Vec<_>>();
            let js = (0..VDMPF::K)
                .map(|i| {
                    let xb = BigUint::from_bytes_be(x) + n.clone() * i.to_biguint().unwrap();
                    let mut xbb = xb.to_bytes_be();
                    self.prp.permu(&mshare.seed, &mut xbb);
                    let yb = BigUint::from_bytes_be(&xbb);
                    let mut index = (yb % &b).to_bytes_be();
                    for _ in 0..(n_prime - index.len()) {
                        index.insert(0, 0);
                    }
                    index
                })
                .collect::<Vec<_>>();
            js.into_iter().zip(is.iter()).for_each(|(j, i)| {
                if let None = dedup.get(&(j.clone(), eta)) {
                    dedup.insert((j.clone(), eta), ());
                    inputs[*i].push((j, eta));
                }
            });
        });
        let mut outputs = vec![BGroup::from(vec![0; self.lambda]); xs.len()];
        let mut pi = BGroup::from(vec![0; self.lambda]);
        inputs.iter().enumerate().for_each(|(i, input)| {
            let js = input.iter().map(|(j, _)| j.as_ref()).collect::<Vec<_>>();
            let (ys, pii) = self.vdpf.eval(b_party, &mshare.ks[i], &js);
            ys.into_iter().zip(input.iter()).for_each(|(y, (_, eta))| {
                outputs[*eta] += y.as_ref();
                let pi_tmp: Vec<u8> = (pi.clone() + pii.as_ref()).into();
                pi += self.hash_prime.gen(&pi_tmp, self.lambda).as_ref();
            });
        });
        let output_vals: Vec<Vec<u8>> = outputs.into_iter().map(|output| output.into()).collect();
        (output_vals, pi.into())
    }

    /// `Verify` in the paper
    pub fn verify(&self, pis: &[&[u8]; 2]) -> bool {
        pis[0] == pis[1]
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
    /// Permutation $[n\mathcal{k} \rightarrow n\mathcal{k}]$ in-place
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
    use std::cell::RefCell;

    use super::*;

    use aes::cipher::generic_array::GenericArray;
    use aes::cipher::{BlockEncrypt, KeyInit};
    use aes::Aes256;
    use rand::prelude::*;
    use rand_chacha::ChaChaRng;
    use rand_seeder::Seeder;

    #[derive(Default)]
    pub struct PRP {}

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

    pub struct CHSampler {
        rng: RefCell<ChaChaRng>,
    }

    impl CHSampler {
        pub fn new(seed: u64) -> Self {
            Self {
                rng: RefCell::new(ChaChaRng::seed_from_u64(seed)),
            }
        }
    }

    impl ISampler for CHSampler {
        fn sample(&self, n: usize) -> usize {
            self.rng.borrow_mut().gen_range(0..n - 1)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::tests_fixture::*;
    use super::*;
    use crate::dpf::tests_fixture::*;

    use hex_literal::hex;

    #[test]
    fn test_gen_eval_verify_ok() {
        let fs = [
            PointFn {
                a: hex!("01b2c3d4a1b2c3d4a1b2c3d4a1b2c3d4").to_vec(),
                b: hex!("e5f67890e5f67890e5f67890e5f67890").to_vec(),
            },
            PointFn {
                a: hex!("01b2c3d4a1b2c3d4a1b2c3d4a1b2c3d5").to_vec(),
                b: hex!("e5f67890e5f67890e5f67890e5f67891").to_vec(),
            },
        ];
        let ch_seed = 7;
        let prp_seed = 8;
        let dpf_seed = 9;
        let vdmpf = VDMPF::new(
            80f64,
            1000,
            Box::new(PRP::default()),
            Box::new(CHSampler::new(ch_seed)),
            Box::new(Sampler::new(prp_seed)),
            1000,
            16,
            Box::new(PRG::default()),
            Box::new(Hash::default()),
            Box::new(Hash::default()),
            Box::new(Hash::default()),
            Box::new(Sampler::new(dpf_seed)),
            1000,
        );
        let gen_res = vdmpf.gen(&fs.iter().collect::<Vec<_>>());
        assert!(gen_res.is_ok(), "VDMPF gen failed");
        let mut mshare0 = gen_res.unwrap();
        let mut mshare1 = mshare0.clone();
        mshare0
            .ks
            .iter_mut()
            .for_each(|k| k.s0s = vec![k.s0s[0].clone()]);
        mshare1
            .ks
            .iter_mut()
            .for_each(|k| k.s0s = vec![k.s0s[1].clone()]);

        let xs = &[
            hex!("01b2c3d4a1b2c3d4a1b2c3d4a1b2c3d4").as_ref(),
            hex!("01b2c3d4a1b2c3d4a1b2c3d4a1b2c3d5").as_ref(),
            hex!("0d3c2b1a4d3c2b1a4d3c2b1a4d3c2b1a").as_ref(),
            hex!("0d3c2b1a4d3c2b1a4d3c2b1a4d3c2b1b").as_ref(),
        ];
        let (y0s, pi0) = vdmpf.eval(false, &mshare0, xs, 2);
        let (y1s, pi1) = vdmpf.eval(true, &mshare1, xs, 2);
        assert_eq!(vdmpf.verify(&[&pi0, &pi1]), true);
        for ((x, y0), y1) in xs.iter().zip(y0s.iter()).zip(y1s.iter()) {
            let y: Vec<u8> = (BGroup::from(y0.to_owned()) + y1.as_ref()).into();
            if x == &fs[0].a {
                assert_eq!(y, fs[0].b);
            } else if x == &fs[1].a {
                assert_eq!(y, fs[1].b);
            } else {
                assert_eq!(y, hex!("00000000000000000000000000000000"));
            }
        }
    }
}
