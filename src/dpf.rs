// Copyright (C) myl7
// SPDX-License-Identifier: Apache-2.0

//! Distributed point function (DPF) implementation

// `s` & `t` as single letter var follow the definition in the paper.
// Left (L) = 0, right (R) = 1.
// For bits, true = 1, false = 0.
// When both L / R and 0 / 1 exist, put 0 / 1 first in indexing.

use std::ops::Add;

use anyhow::anyhow;
use bitvec::prelude::*;

use crate::group::{xor_bool, xor_bytes, BGroup};

/// `VerDPF` in the paper
pub struct VDPF {
    /// $\lambda$ in the paper.
    /// The byte len of point function a and b, and sampled seeds should be equal to it.
    lambda: usize,
    /// PRG (pseudo-random generator), $\mathcal{G}$: $\{0, 1\}^{\lambda} \rightarrow \{0, 1\}^{2\lambda + 2}$ in the paper
    prg: Box<dyn PRG>,
    /// Hash function that is both collision-resistant and xor-collision-resistant, $H$: $\{0, 1\}^{n + \lambda} \rightarrow \{0, 1\}^{4\lambda}$ in the paper
    hash: Box<dyn Hash>,
    /// Hash function that is collision-resistant, $H'$: $\{0, 1\}^{4\lambda} \rightarrow \{0, 1\}^{2\lambda}$ in the paper
    hash_prime: Box<dyn HashPrime>,
    /// To sample starting seeds
    sampler: Box<dyn BSampler>,
    /// Retry limit for starting seed resampling
    gen_retry: usize,
}

impl VDPF {
    pub fn new(
        lambda: usize,
        prg: Box<dyn PRG>,
        hash: Box<dyn Hash>,
        hash_prime: Box<dyn HashPrime>,
        sampler: Box<dyn BSampler>,
        gen_retry: usize,
    ) -> Self {
        Self {
            lambda,
            prg,
            hash,
            hash_prime,
            sampler,
            gen_retry,
        }
    }
}

/// `VerDPF` API
impl VDPF {
    /// `Gen` in the paper
    pub fn gen(&self, f: PointFn) -> anyhow::Result<Share> {
        for _ in 0..self.gen_retry {
            let mut s0_buf = self.sampler.sample(self.lambda * 2);
            let s01 = s0_buf.split_off(self.lambda);
            let s00 = s0_buf;
            let s0s = vec![s00.clone(), s01.clone()];

            let n = f.a.view_bits::<Msb0>().len();
            let mut nodes = [(s00, false), (s01, true)];
            let mut nodes_buf = nodes.clone();
            let mut cws = vec![
                CW {
                    s: vec![0; self.lambda],
                    ts: [false, false]
                };
                n
            ];
            for i in 0..n {
                self.cw_gen(
                    &mut nodes,
                    f.a.view_bits::<Msb0>()[i],
                    &mut nodes_buf,
                    &mut cws[i],
                );
            }
            let mut pi0 = f.a.clone();
            pi0.extend_from_slice(&vec![0; 3 * self.lambda]);
            self.hash.gen(&mut pi0, &nodes[0].0);
            let mut pi1 = f.a.clone();
            pi1.extend_from_slice(&vec![0; 3 * self.lambda]);
            self.hash.gen(&mut pi1, &nodes[1].0);
            let cs: Vec<u8> = (BGroup::from(pi0) + pi1.as_ref()).into();
            nodes[0].1 = nodes[0].0.view_bits::<Lsb0>()[0];
            nodes[1].1 = nodes[1].0.view_bits::<Lsb0>()[0];
            if nodes[0].1 == nodes[1].1 {
                continue;
            }
            // Since we use xor as plus, -a == a in the group
            let ocw: Vec<u8> =
                (BGroup::from(f.b) + nodes[0].0.as_ref() + nodes[1].0.as_ref()).into();
            return Ok(Share { s0s, cws, cs, ocw });
        }
        Err(anyhow!("VDPF fails"))
    }

    /// `BVEval` in the paper.
    /// `b` is the party num, which is 0 / 1.
    pub fn eval(&self, b: bool, share: &Share, xs: &[&[u8]]) -> (Vec<Vec<u8>>, Vec<u8>) {
        assert_eq!(share.s0s.len(), 1);

        let mut ys: Vec<Vec<u8>> = vec![];
        let mut pi = share.cs.clone();
        let mut hash_prime_buf = vec![0; self.lambda * 4];
        for x in xs {
            let mut node = (share.s0s[0].clone(), b);
            let mut node_buf = node.clone();
            for i in 0..x.view_bits::<Msb0>().len() {
                self.node_expand(&mut node, &mut node_buf, &share.cws[i]);
                if x.view_bits::<Msb0>()[i] {
                    node = node_buf.clone();
                }
            }
            hash_prime_buf[0..self.lambda].copy_from_slice(x);
            self.hash.gen(&mut hash_prime_buf, &node.0);
            node.1 = node.0.view_bits::<Lsb0>()[0];
            // Since we use xor as plus, -a == a in the group
            ys.push(correct(BGroup::from(node.0), share.ocw.as_ref(), node.1).into());
            correct_bytes(&mut hash_prime_buf, &share.cs, node.1);
            xor_bytes(&mut hash_prime_buf, &pi);
            xor_bytes(&mut pi, &mut hash_prime_buf)
        }
        (ys, pi.into())
    }

    /// `Verify` in the paper
    pub fn verify(&self, pis: &[&[u8]; 2]) -> bool {
        pis[0] == pis[1]
    }
}

impl VDPF {
    /// `NodeExpand` in the paper
    fn node_expand(&self, (sl, tl): &mut (Vec<u8>, bool), (sr, tr): &mut (Vec<u8>, bool), cw: &CW) {
        let t = *tl;
        let [t0, t1] = self.prg.gen(sl, sr);
        correct_bytes(sl, &cw.s, t);
        *tl = correct_bit(t0, cw.ts[0], t);
        correct_bytes(sr, &cw.s, t);
        *tr = correct_bit(t1, cw.ts[1], t);
    }

    /// `CWGen` in the paper
    fn cw_gen(
        &self,
        [(s0l, t0), (s0r, t1)]: &mut [(Vec<u8>, bool); 2],
        x: bool,
        [(s1l, _), (s1r, _)]: &mut [(Vec<u8>, bool); 2],
        cw: &mut CW,
    ) {
        let [t1l, t1r] = self.prg.gen(s0r, s1r);
        *s1l = s0r.clone();
        let [t0l, t0r] = self.prg.gen(s0l, s0r);
        let ss = [[s0l, s0r], [s1l, s1r]];
        let ts = [[t0l, t0r], [t1l, t1r]];
        let (diff, same) = if x { (1, 0) } else { (0, 1) };
        cw.s.copy_from_slice(ss[0][same]);
        xor_bytes(&mut cw.s, ss[1][same]);
        let tcs: [bool; 2] = [
            xor_bool(xor_bool(xor_bool(ts[0][0], ts[1][0]), true), x),
            xor_bool(xor_bool(ts[0][1], ts[1][1]), x),
        ];
        correct_bytes(ss[0][diff], &cw.s, *t0);
        if diff == 1 {
            *ss[0][0] = ss[0][diff].clone();
        }
        correct_bytes(ss[1][diff], &cw.s, *t1);
        *ss[0][1] = ss[1][diff].clone();
        *t0 = correct_bit(ts[0][diff], tcs[diff], *t0);
        *t1 = correct_bit(ts[1][diff], tcs[diff], *t1);
        cw.ts = tcs;
    }
}

/// `correct` in the paper
fn correct_bytes(a: &mut [u8], b: &[u8], t: bool) {
    if t {
        xor_bytes(a, b)
    }
}

/// `correct` in the paper
fn correct_bit(a: bool, b: bool, t: bool) -> bool {
    if t {
        xor_bool(a, b)
    } else {
        a
    }
}

/// `correct` in the paper
fn correct<Lhs, Rhs>(a: Lhs, b: Rhs, t: bool) -> Lhs
where
    Lhs: Add<Rhs, Output = Lhs>,
{
    if t {
        a + b
    } else {
        a
    }
}

/// PRG trait. See [`VDPF`].
pub trait PRG {
    fn gen(&self, x: &mut [u8], x2: &mut [u8]) -> [bool; 2];
}

/// Hash trait. See [`VDPF`].
pub trait Hash {
    fn gen(&self, x0: &mut [u8], x1: &[u8]);
}

/// Hash prime trait. See [`VDPF`].
pub trait HashPrime {
    fn gen(&self, x: &mut [u8]);
}

/// Sample random bytes.
pub trait BSampler {
    /// Use `&self` because it does not need to remember any states
    fn sample(&self, len: usize) -> Vec<u8>;
}

/// Point function, $f_{\alpha, \beta}$ in the paper
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PointFn {
    /// $\alpha$ in the paper
    pub a: Vec<u8>,
    /// $\beta$ in the paper
    pub b: Vec<u8>,
}

/// `k` in the paper.
/// Since $k_0$, $k_1$ share all other fields except `s0s`, when generating, `s0s` contains all 2 start seeds, and when evaluating, `s0s` contains only 1.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Share {
    pub s0s: Vec<Vec<u8>>,
    pub cws: Vec<CW>,
    pub cs: Vec<u8>,
    pub ocw: Vec<u8>,
}

/// Correction word, `cw` in the paper
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CW {
    pub s: Vec<u8>,
    pub ts: [bool; 2],
}

#[cfg(all(test, feature = "dyn_utils"))]
mod tests {
    use super::*;

    use hex_literal::hex;

    use crate::dyn_utils::*;

    #[test]
    fn run_ok() {
        let f = PointFn {
            a: hex!("a1b2c3d4a1b2c3d4a1b2c3d4a1b2c3d4").to_vec(),
            b: hex!("e5f67890e5f67890e5f67890e5f67890").to_vec(),
        };
        let seed = 7;
        let vdpf = VDPF::new(
            16,
            Box::new(Aes128PRG::default()),
            Box::new(Aes128Hash::default()),
            Box::new(Aes128Hash::default()),
            Box::new(ChaChaBSampler::new(seed)),
            1000,
        );
        let gen_res = vdpf.gen(f.clone());
        assert!(gen_res.is_ok(), "VDPF gen failed");
        let mut share = gen_res.unwrap();

        let s00 = share.s0s[0].clone();
        let s01 = share.s0s[1].clone();
        let xs = &[
            hex!("a1b2c3d4a1b2c3d4a1b2c3d4a1b2c3d5").as_ref(),
            hex!("a1b2c3d4a1b2c3d4a1b2c3d4a1b2c3d4").as_ref(),
            hex!("4d3c2b1a4d3c2b1a4d3c2b1a4d3c2b1a").as_ref(),
        ];
        share.s0s = vec![s01];
        let (y1s, pi1) = vdpf.eval(true, &share, xs);
        share.s0s = vec![s00];
        let (y0s, pi0) = vdpf.eval(false, &share, xs);
        assert_eq!(vdpf.verify(&[&pi0, &pi1]), true);
        for ((x, y0), y1) in xs.iter().zip(y0s.iter()).zip(y1s.iter()) {
            let y: Vec<u8> = (BGroup::from(y0.to_owned()) + y1.as_ref()).into();
            if x == &f.a {
                assert_eq!(y, f.b);
            } else {
                assert_eq!(y, hex!("00000000000000000000000000000000"));
            }
        }
    }

    // #[test]
    // fn ab_diff_len() {
    //     let f = PointFn {
    //         a: hex!("a1b2c3d4a1b2c3d4").to_vec(),
    //         b: hex!("e5f67890e5f67890e5f67890e5f67890").to_vec(),
    //     };
    //     let seed = 7;
    //     let vdpf = VDPF::new(
    //         16,
    //         Box::new(ChaChaPRG::default()),
    //         Box::new(Shake256Hash::default()),
    //         Box::new(Shake256Hash::default()),
    //         Box::new(ChaChaBSampler::new(seed)),
    //         1000,
    //     );
    //     let gen_res = vdpf.gen(f.clone());
    //     assert!(gen_res.is_ok(), "VDPF gen failed");
    //     let mut share = gen_res.unwrap();

    //     let s00 = share.s0s[0].clone();
    //     let s01 = share.s0s[1].clone();
    //     let xs = &[
    //         hex!("a1b2c3d4a1b2c3d4").as_ref(),
    //         hex!("a1b2c3d4a1b2c3d4").as_ref(),
    //         hex!("4d3c2b1a4d3c2b1a").as_ref(),
    //     ];
    //     share.s0s = vec![s01];
    //     let (y1s, pi1) = vdpf.eval(true, &share, xs);
    //     share.s0s = vec![s00];
    //     let (y0s, pi0) = vdpf.eval(false, &share, xs);
    //     assert_eq!(vdpf.verify(&[&pi0, &pi1]), true);
    //     for ((x, y0), y1) in xs.iter().zip(y0s.iter()).zip(y1s.iter()) {
    //         let y: Vec<u8> = (BGroup::from(y0.to_owned()) + y1.as_ref()).into();
    //         if x == &f.a {
    //             assert_eq!(y, f.b);
    //         } else {
    //             assert_eq!(y, hex!("00000000000000000000000000000000"));
    //         }
    //     }
    // }

    // #[test]
    // fn ab_zeros() {
    //     let f = PointFn {
    //         a: hex!("00000000000000000000000000000000").to_vec(),
    //         b: hex!("00000000000000000000000000000000").to_vec(),
    //     };
    //     let seed = 7;
    //     let vdpf = VDPF::new(
    //         16,
    //         Box::new(ChaChaPRG::default()),
    //         Box::new(Shake256Hash::default()),
    //         Box::new(Shake256Hash::default()),
    //         Box::new(ChaChaBSampler::new(seed)),
    //         1000,
    //     );
    //     let gen_res = vdpf.gen(f.clone());
    //     assert!(gen_res.is_ok(), "VDPF gen failed");
    // }
}
