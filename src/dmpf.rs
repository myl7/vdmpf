// Copyright (C) myl7
// SPDX-License-Identifier: Apache-2.0

//! Distributed multi-point function (DMPF) implementation

use std::collections::HashMap;
use std::mem;

use anyhow::anyhow;
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
    /// $H'$ used by VDMPF particularly
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
        hash_prime_dmpf: Box<dyn Gen>,
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
            hash_prime: hash_prime_dmpf,
        }
    }
}

/// `VerDMPF` API
impl VDMPF {
    const K: usize = 3;

    /// `Gen` in the paper
    pub fn gen(&self, fs: &[&PointFn]) -> anyhow::Result<MShare> {
        // Ensure $\alpha < 2^126$ so that $n\mathcal{k} < 2^128$ and we can use AES as PRP
        // TODO: Guard according to the permu interface.
        // fs.iter().for_each(|f| assert!(f.a[0] & 0xc0 == 0));

        // pow for generating `n` and `b`
        let (n_pow, bn_pow) = self.prp.domains();
        let n = 2.to_biguint().unwrap().pow(n_pow);
        let t = fs.len();
        // Use `usize` other than `f64` because all of them will be used as (part of) indexes
        let m = self.ch_bucket(t);
        let b = 2
            .to_biguint()
            .unwrap()
            .pow(bn_pow)
            .div_ceil(&m.to_biguint().unwrap());
        // $n' / 8$ in the paper
        let n_prime = b.bits() as usize;
        let n_prime_byte = (b.bits() as f64 / 8f64).ceil() as usize;

        // + 2 to use 0 as the boundary
        let mut prp_retry = self.prp_retry + 2;
        let (table, seed) = loop {
            prp_retry -= 1;
            if prp_retry == 0 {
                return Err(anyhow!("VDMPF fails"));
            }
            let seed = self.prp_sampler.sample(self.lambda);
            let hs = (0..VDMPF::K)
                .map(|i| {
                    let n = n.clone();
                    let b = b.clone();
                    let seed = seed.clone();
                    move |x: &[u8]| {
                        let yb = self.prp_result(x, &n, i, &seed);
                        usize_biguint(yb / &b)
                    }
                })
                .collect::<Vec<_>>();
            let a_s = fs.iter().map(|f| f.a.as_ref()).collect::<Vec<_>>();
            let table = match self.ch_compact(&a_s, &hs, m) {
                Ok(table) => table,
                Err(_) => continue,
            };
            break (table, seed);
        };
        let mut mshare = MShare { ks: vec![], seed };
        for item in table.into_iter() {
            let (a, b) = match item {
                // `aj` is $a_j$ in the paper
                Some((aji, k)) => {
                    let yb = self.prp_result(fs[aji].a.as_ref(), &n, k, &mshare.seed);
                    let mut index = (yb % &b).to_bytes_be();
                    pad_biguint(&mut index, n_prime_byte);
                    let a = index;
                    let b = fs[aji].b.clone();
                    (a, b)
                }
                None => (vec![0; n_prime_byte], vec![0; self.lambda]),
            };
            let f = PointFn {
                a,
                b,
                a_leap: Some(n_prime_byte * 8 - n_prime),
            };
            let share = self.vdpf.gen(f)?;
            mshare.ks.push(share);
        }
        Ok(mshare)
    }

    /// `BVEval` in the paper.
    /// `t` presents possible points of the MPF, and it should work if `t` is larger than the actual one generally,
    /// but now if `t` is larger/smaller, evaluation fails.
    pub fn eval(
        &self,
        b_party: bool,
        mshare: &MShare,
        xs: &[&[u8]],
        t: usize,
    ) -> (Vec<Vec<u8>>, Vec<u8>) {
        assert_eq!(mshare.ks.first().map(|k| k.s0s.len()), Some(1));

        let (n_pow, bn_pow) = self.prp.domains();
        let n = 2.to_biguint().unwrap().pow(n_pow);
        // Use `usize` other than `f64` because all of them will be used as (part of) indexes
        let m = self.ch_bucket(t);
        let b = 2
            .to_biguint()
            .unwrap()
            .pow(bn_pow)
            .div_ceil(&m.to_biguint().unwrap());
        // $n' / 8$ in the paper
        let n_prime = b.bits() as usize;
        let n_prime_byte = (b.bits() as f64 / 8f64).ceil() as usize;

        let mut inputs = vec![vec![]; m];
        let mut dedup = HashMap::new();
        xs.iter().enumerate().for_each(|(eta, x)| {
            let is = (0..VDMPF::K)
                .map(|i| {
                    let yb = self.prp_result(x, &n, i, &mshare.seed);
                    usize_biguint(yb / &b)
                })
                .collect::<Vec<_>>();
            let js = (0..VDMPF::K)
                .map(|i| {
                    let yb = self.prp_result(x, &n, i, &mshare.seed);
                    let mut index = (yb % &b).to_bytes_be();
                    pad_biguint(&mut index, n_prime_byte);
                    index
                })
                .collect::<Vec<_>>();
            js.into_iter().zip(is.iter()).for_each(|(j, i)| {
                if dedup.get(&(j.clone(), eta)).is_none() {
                    dedup.insert((j.clone(), eta), ());
                    inputs[*i].push((j, eta));
                }
            });
        });
        let mut outputs = vec![BGroup::from(vec![0; self.lambda]); xs.len()];
        let mut pi = BGroup::from(vec![0; self.lambda]);
        inputs.iter().enumerate().for_each(|(i, input)| {
            let js = input.iter().map(|(j, _)| j.as_ref()).collect::<Vec<_>>();
            let (ys, pii) = self.vdpf.eval(
                b_party,
                &mshare.ks[i],
                &js,
                Some(n_prime_byte * 8 - n_prime),
            );
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
    /// Use `usize` other than `f64` because it will be used as (part of) indexes.
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

    fn prp_result(&self, x: &[u8], n: &BigUint, i: usize, seed: &[u8]) -> BigUint {
        let xb = BigUint::from_bytes_be(x) + n * i.to_biguint().unwrap();
        let mut xbb = xb.to_bytes_be();
        pad_biguint(&mut xbb, 16);
        self.prp.permu(seed, &mut xbb);
        BigUint::from_bytes_be(&xbb)
    }
}

fn pad_biguint(b: &mut Vec<u8>, n: usize) {
    for _ in 0..(n - b.len()) {
        b.insert(0, 0);
    }
}

fn usize_biguint(b: BigUint) -> usize {
    let bs = b.to_u64_digits();
    if !bs.is_empty() {
        bs[0] as usize
    } else {
        0
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
    /// It returns the domains of the input $x / \mathcal{k}$ and output.
    /// $\mathcal{k}$ is fixed to 3. See [`VDMPF`].
    /// The output domain should be larger than the input one.
    /// Both the domains are presented as the power of 2, which limitted the possible domain size values,
    /// but it should be reasonable for most use cases.
    fn domains(&self) -> (u32, u32);

    /// Permutation $[n\mathcal{k} \rightarrow n\mathcal{k}]$ in the paper.
    /// In-place because `x` has a fixed size.
    /// However, the output domain can be a little larger than the input one as long as the space is enough.
    /// We suggest using an AES-based PRP as the tests for convenience.
    fn permu(&self, seed: &[u8], x: &mut [u8]);
}

/// Interface for sampling indexes in Cuckoo-hashing.
/// See [`VDMPF`].
pub trait ISampler {
    /// Sample from `[0, n)`.
    /// Use `&self`. See [`crate::dpf::BSampler::sample`] for the reason.
    fn sample(&self, n: usize) -> usize;
}

#[cfg(all(test, feature = "dyn_utils"))]
mod tests {
    use super::*;

    use hex_literal::hex;

    use crate::dyn_utils::*;

    #[test]
    fn run_ok() {
        let fs = [
            PointFn {
                a: hex!("01b2c3d4a1b2c3d4a1b2c3d4a1b2c3d4").to_vec(),
                b: hex!("e5f67890e5f67890e5f67890e5f67890").to_vec(),
                a_leap: None,
            },
            PointFn {
                a: hex!("01b2c3d4a1b2c3d4a1b2c3d4a1b2c3d5").to_vec(),
                b: hex!("e5f67890e5f67890e5f67890e5f67891").to_vec(),
                a_leap: None,
            },
        ];
        let ch_seed = 7;
        let prp_seed = 8;
        let dpf_seed = 9;
        let vdmpf = VDMPF::new(
            80f64,
            1000,
            Box::new(Aes256PRP::default()),
            Box::new(ChaChaISampler::new(ch_seed)),
            Box::new(ChaChaBSampler::new(prp_seed)),
            1000,
            16,
            Box::new(ChaChaPRG::default()),
            Box::new(Shake256Hash::default()),
            Box::new(Shake256Hash::default()),
            Box::new(Shake256Hash::default()),
            Box::new(ChaChaBSampler::new(dpf_seed)),
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
        let (y0s, pi0) = vdmpf.eval(false, &mshare0, xs, fs.len());
        let (y1s, pi1) = vdmpf.eval(true, &mshare1, xs, fs.len());
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
