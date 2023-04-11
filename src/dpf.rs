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
    lambda: usize,
    /// PRG (pseudo-random generator), $\mathcal{G}$: $\{0, 1\}^{\lambda} \rightarrow \{0, 1\}^{2\lambda + 2}$ in the paper
    prg: Box<dyn Gen>,
    /// Hash function that is both collision-resistant and xor-collision-resistant, $H$: $\{0, 1\}^{n + \lambda} \rightarrow \{0, 1\}^{4\lambda}$ in the paper
    hash: Box<dyn Gen>,
    /// Hash function that is collision-resistant, $H'$: $\{0, 1\}^{4\lambda} \rightarrow \{0, 1\}^{2\lambda}$ in the paper
    hash_prime: Box<dyn Gen>,
}

impl VDPF {
    pub fn new(
        lambda: usize,
        prg: Box<dyn Gen>,
        hash: Box<dyn Gen>,
        hash_prime: Box<dyn Gen>,
    ) -> Self {
        Self {
            lambda,
            prg,
            hash,
            hash_prime,
        }
    }
}

/// `VerDPF` API
impl VDPF {
    /// `Gen` in the paper.
    /// Starting seeds `s0s` should be randomly sampled, but not required to be different.
    /// They should be both $\lambda$ bytes, which is the `lambda` field in the struct, otherwise panic.
    pub fn gen(&self, s0s: [Vec<u8>; 2], f: PointFn) -> Result<Share, ()> {
        for s0 in s0s.iter() {
            assert_eq!(s0.len(), self.lambda);
        }

        let mut nodes = [vec![(s0s[0].clone(), false)], vec![(s0s[1].clone(), true)]];
        let mut cws = vec![];
        let n = f.a.view_bits::<Lsb0>().len();
        for i in 0..n {
            let (cw, [node0, node1]) =
                self.cw_gen(&[&nodes[0][i], &nodes[1][i]], f.a.view_bits::<Lsb0>()[i]);
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
            return Err(());
        }
        // Since we use xor as plus, -a == a in the group
        let ocw: Vec<u8> =
            (BGroup::from(f.b) + nodes[0][n + 1].0.as_ref() + nodes[0][n + 1].0.as_ref()).into();
        Ok(Share {
            s0s: s0s.to_vec(),
            cws,
            cs,
            ocw,
        })
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
    pub fn verify(pis: [Vec<u8>; 2]) -> bool {
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
/// See [`VDPF`].
pub trait Gen {
    fn gen(&self, input: &[u8], output_len: usize) -> Vec<u8>;
}

impl dyn Gen {
    fn prg_gen(&self, lambda: usize, input: &[u8]) -> [(Vec<u8>, bool); 2] {
        // Other than (s, t, s, t) in the paper, here we take (s, s, t, t) for convenience
        let mut buf = self.gen(input, lambda * 2 + 1);
        let t_buf = buf.split_off(2 * lambda);
        let ts = [t_buf.view_bits::<Lsb0>()[0], t_buf.view_bits::<Lsb0>()[1]];
        let mut s_buf = buf.split_off(lambda);
        let sr = buf;
        s_buf.split_off(lambda).truncate(0);
        let sl = s_buf;
        [(sl, ts[0]), (sr, ts[1])]
    }
}

/// Point function, $f_{\alpha, \beta}$ in the paper
pub struct PointFn {
    /// $\alpha$ in the paper
    pub a: Vec<u8>,
    /// $\beta$ in the paper
    pub b: Vec<u8>,
}

/// `k` in the paper.
/// Since $k_0$, $k_1$ share all other fields except `s0s`, when generating, `s0s` contains all 2 start seeds, and when evaluating, `s0s` contains only 1.
pub struct Share {
    pub s0s: Vec<Vec<u8>>,
    pub cws: Vec<CW>,
    pub cs: Vec<u8>,
    pub ocw: Vec<u8>,
}

/// Correction word, `cw` in the paper
#[derive(Clone)]
pub struct CW {
    pub s: Vec<u8>,
    pub ts: [bool; 2],
}
