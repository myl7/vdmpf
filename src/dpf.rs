// Copyright (C) myl7
// SPDX-License-Identifier: Apache-2.0

//! Distributed point function (DPF) implementation

// `s` & `t` as single letter var follow the definition in the paper.
// Left (L) = 0, right (R) = 1.
// For bits, true = 1, false = 0.
// When both L / R and 0 / 1 exist, put 0 / 1 first in indexing.

use std::ops::Add;

use bitvec::prelude::*;

use crate::group::BGroup;

/// `VerDPF` in the paper
pub struct VDPF {
    /// $\lambda$ in the paper.
    /// The byte len of point function a and b, and sampled seeds should be equal to it.
    lambda: usize,
    /// PRG (pseudo-random generator), $\mathcal{G}$: $\{0, 1\}^{\lambda} \rightarrow \{0, 1\}^{2\lambda + 2}$ in the paper
    prg: Box<dyn Gen>,
    /// Hash function that is both collision-resistant and xor-collision-resistant, $H$: $\{0, 1\}^{n + \lambda} \rightarrow \{0, 1\}^{4\lambda}$ in the paper
    hash: Box<dyn Gen>,
    /// Hash function that is collision-resistant, $H'$: $\{0, 1\}^{4\lambda} \rightarrow \{0, 1\}^{2\lambda}$ in the paper
    hash_prime: Box<dyn Gen>,
    /// To sample starting seeds
    sampler: Box<dyn BSampler>,
    /// Retry limit for starting seed resampling
    gen_retry: usize,
}

impl VDPF {
    pub fn new(
        lambda: usize,
        prg: Box<dyn Gen>,
        hash: Box<dyn Gen>,
        hash_prime: Box<dyn Gen>,
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
    /// `Gen` in the paper.
    /// Starting seeds `s0s` should be randomly sampled, but not required to be different.
    /// They should be both $\lambda$ bytes, which is the `lambda` field in the struct, otherwise panic.
    pub fn gen(&self, f: PointFn) -> Result<Share, ()> {
        for _ in 0..self.gen_retry {
            let mut s0_buf = self.sampler.sample(self.lambda * 2);
            let s01 = s0_buf.split_off(self.lambda);
            let s00 = s0_buf.split_off(self.lambda);
            let s0s = vec![s00.clone(), s01.clone()];

            let mut nodes = [vec![(s00, false)], vec![(s01, true)]];
            let mut cws = vec![];
            let n = f.a.view_bits::<Msb0>().len();
            for i in 0..n {
                let (cw, [node0, node1]) =
                    self.cw_gen(&[&nodes[0][i], &nodes[1][i]], f.a.view_bits::<Msb0>()[i]);
                nodes[0].push(node0);
                nodes[1].push(node1);
                cws.push(cw);
            }
            let pi0 = self.hash.gen(
                &[f.a.clone(), nodes[0][n].0.clone()].concat(),
                self.lambda * 4,
            );
            let pi1 = self.hash.gen(
                &[f.a.clone(), nodes[1][n].0.clone()].concat(),
                self.lambda * 4,
            );
            let cs: Vec<u8> = (BGroup::from(pi0) + pi1.as_ref()).into();
            nodes[0].push((nodes[0][n].0.clone(), nodes[0][n].0.view_bits::<Lsb0>()[0]));
            nodes[1].push((nodes[1][n].0.clone(), nodes[1][n].0.view_bits::<Lsb0>()[0]));
            if nodes[0][n + 1].1 == nodes[1][n + 1].1 {
                continue;
            }
            // Since we use xor as plus, -a == a in the group
            let ocw: Vec<u8> =
                (BGroup::from(f.b) + nodes[0][n + 1].0.as_ref() + nodes[1][n + 1].0.as_ref())
                    .into();
            return Ok(Share { s0s, cws, cs, ocw });
        }
        Err(())
    }

    /// `BVEval` in the paper.
    /// `b` is the party num, which is 0 / 1.
    pub fn eval(&self, b: bool, share: &Share, xs: &[&[u8]]) -> (Vec<Vec<u8>>, Vec<u8>) {
        assert_eq!(share.s0s.len(), 1);

        let mut ys: Vec<Vec<u8>> = vec![];
        let mut pi = BGroup::from(share.cs.clone());
        for x in xs {
            let mut node = (share.s0s[0].clone(), b);
            for i in 0..x.view_bits::<Msb0>().len() {
                let [node0, node1] = self.node_expand(&node, share.cws[i].clone());
                if x.view_bits::<Msb0>()[i] {
                    node = node1;
                } else {
                    node = node0;
                }
            }
            let pi_tmp = self
                .hash
                .gen(&[x.to_vec(), node.0.clone()].concat(), self.lambda * 4);
            node.1 = node.0.view_bits::<Lsb0>()[0];
            // Since we use xor as plus, -a == a in the group
            ys.push(correct(BGroup::from(node.0), share.ocw.as_ref(), node.1).into());
            let hash_prime_buf: Vec<u8> =
                (pi.clone() + correct(BGroup::from(pi_tmp), share.cs.as_ref(), node.1)).into();
            pi += self
                .hash_prime
                .gen(&hash_prime_buf, self.lambda * 2)
                .as_ref();
        }
        (ys, pi.into())
    }

    /// `Verify` in the paper.
    pub fn verify(&self, pis: &[&[u8]; 2]) -> bool {
        pis[0] == pis[1]
    }
}

impl VDPF {
    /// `NodeExpand` in the paper
    fn node_expand(&self, (s, t): &(Vec<u8>, bool), cw: CW) -> [(Vec<u8>, bool); 2] {
        let [(sl, tl), (sr, tr)] = self.prg.prg_gen(self.lambda, s);
        [
            (
                correct(BGroup::from(sl), cw.s.as_ref(), *t).into(),
                correct(BGroup::from(tl), cw.ts[0], *t).into(),
            ),
            (
                correct(BGroup::from(sr), cw.s.as_ref(), *t).into(),
                correct(BGroup::from(tr), cw.ts[1], *t).into(),
            ),
        ]
    }

    /// `CWGen` in the paper
    fn cw_gen(
        &self,
        [(s0, t0), (s1, t1)]: &[&(Vec<u8>, bool); 2],
        x: bool,
    ) -> (CW, [(Vec<u8>, bool); 2]) {
        let [(s0l, t0l), (s0r, t0r)] = self.prg.prg_gen(self.lambda, s0);
        let [(s1l, t1l), (s1r, t1r)] = self.prg.prg_gen(self.lambda, s1);
        let ss = [[s0l, s0r], [s1l, s1r]];
        let ts = [[t0l, t0r], [t1l, t1r]];
        let (diff, same) = if x { (1, 0) } else { (0, 1) };
        let sc: Vec<u8> = (BGroup::from(ss[0][same].clone()) + ss[1][same].as_ref()).into();
        let tcs: [bool; 2] = [
            (BGroup::from(ts[0][0]) + ts[1][0] + true + x).into(),
            (BGroup::from(ts[0][1]) + ts[1][1] + x).into(),
        ];
        let nodes = [
            (
                correct(BGroup::from(ss[0][diff].clone()), sc.as_ref(), *t0).into(),
                correct(BGroup::from(ts[0][diff]), tcs[diff], *t0).into(),
            ),
            (
                correct(BGroup::from(ss[1][diff].clone()), sc.as_ref(), *t1).into(),
                correct(BGroup::from(ts[1][diff]), tcs[diff], *t1).into(),
            ),
        ];
        let cw = CW { s: sc, ts: tcs };
        (cw, nodes)
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

/// Interface for PRG and Hash functions.
/// For PRG, `input` is the seed.
/// See [`VDPF`].
pub trait Gen {
    fn gen(&self, input: &[u8], output_len: usize) -> Vec<u8>;
}

impl dyn Gen {
    fn prg_gen(&self, lambda: usize, input: &[u8]) -> [(Vec<u8>, bool); 2] {
        // Other than (s, t, s, t) in the paper, here we take (s, s, t, t) for convenience
        let mut buf = self.gen(input, lambda * 2 + 1);
        let t_buf = buf.split_off(2 * lambda);
        let ts = [t_buf.view_bits::<Msb0>()[0], t_buf.view_bits::<Msb0>()[1]];
        let mut s_buf = buf.split_off(lambda);
        let sr = buf;
        s_buf.split_off(lambda).truncate(0);
        let sl = s_buf;
        [(sl, ts[0]), (sr, ts[1])]
    }
}

/// Sample random bytes.
pub trait BSampler {
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

#[cfg(test)]
pub(crate) mod tests_fixture {
    use std::cell::RefCell;

    use super::*;

    use rand::prelude::*;
    use rand_chacha::ChaChaRng;
    use rand_seeder::Seeder;
    use sha3::digest::{ExtendableOutput, Update, XofReader};
    use sha3::Shake256;

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

    impl Sampler {
        pub fn new(seed: u64) -> Self {
            Self {
                rng: RefCell::new(ChaChaRng::seed_from_u64(seed)),
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
}

#[cfg(test)]
mod tests {
    use super::tests_fixture::*;
    use super::*;

    use hex_literal::hex;

    #[test]
    fn test_gen_eval_verify_ok() {
        let f = PointFn {
            a: hex!("a1b2c3d4a1b2c3d4a1b2c3d4a1b2c3d4").to_vec(),
            b: hex!("e5f67890e5f67890e5f67890e5f67890").to_vec(),
        };
        let seed = 7;
        let vdpf = VDPF::new(
            16,
            Box::new(PRG::default()),
            Box::new(Hash::default()),
            Box::new(Hash::default()),
            Box::new(Sampler::new(seed)),
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

    #[test]
    fn test_gen_eval_verify_ok_with_diff_len() {
        let f = PointFn {
            a: hex!("a1b2c3d4a1b2c3d4").to_vec(),
            b: hex!("e5f67890e5f67890e5f67890e5f67890").to_vec(),
        };
        let seed = 7;
        let vdpf = VDPF::new(
            16,
            Box::new(PRG::default()),
            Box::new(Hash::default()),
            Box::new(Hash::default()),
            Box::new(Sampler::new(seed)),
            1000,
        );
        let gen_res = vdpf.gen(f.clone());
        assert!(gen_res.is_ok(), "VDPF gen failed");
        let mut share = gen_res.unwrap();

        let s00 = share.s0s[0].clone();
        let s01 = share.s0s[1].clone();
        let xs = &[
            hex!("a1b2c3d4a1b2c3d4").as_ref(),
            hex!("a1b2c3d4a1b2c3d4").as_ref(),
            hex!("4d3c2b1a4d3c2b1a").as_ref(),
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

    #[test]
    fn test_gen_zeros() {
        let f = PointFn {
            a: hex!("00000000000000000000000000000000").to_vec(),
            b: hex!("00000000000000000000000000000000").to_vec(),
        };
        let seed = 7;
        let vdpf = VDPF::new(
            16,
            Box::new(PRG::default()),
            Box::new(Hash::default()),
            Box::new(Hash::default()),
            Box::new(Sampler::new(seed)),
            1000,
        );
        let gen_res = vdpf.gen(f.clone());
        assert!(gen_res.is_ok(), "VDPF gen failed");
    }
}
