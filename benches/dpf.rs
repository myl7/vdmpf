use criterion::{black_box, criterion_group, criterion_main, Criterion};

use rand::prelude::*;

use vdmpf::dpf::{PointFn, Share, VDPF};
use vdmpf::dyn_utils::*;

fn gen(f: PointFn) -> [Share; 2] {
    let vdpf = VDPF::new(
        16,
        Box::new(Aes128PRG::default()),
        Box::new(Aes128Hash::default()),
        Box::new(Aes128Hash::default()),
        Box::new(ChaChaBSampler::default()),
        1000,
    );
    let mut mshare0 = vdpf.gen(f).unwrap();
    let mut mshare1 = mshare0.clone();
    mshare0.s0s = vec![mshare0.s0s[0].clone()];
    mshare1.s0s = vec![mshare1.s0s[1].clone()];
    [mshare0, mshare1]
}

fn eval_xs100(share: Share, xs: Vec<Vec<u8>>) -> (Vec<Vec<u8>>, Vec<u8>) {
    let vdpf = VDPF::new(
        16,
        Box::new(Aes128PRG::default()),
        Box::new(Aes128Hash::default()),
        Box::new(Aes128Hash::default()),
        Box::new(ChaChaBSampler::default()),
        1000,
    );
    let xs_ref = xs.iter().map(|x| x.as_ref()).collect::<Vec<_>>();
    let (y0s, pi0) = vdpf.eval(false, &share, &xs_ref);
    (y0s, pi0)
}

fn criterion_benches(c: &mut Criterion) {
    let mut f_rng = thread_rng();
    let mut a = vec![0; 16];
    f_rng.fill_bytes(&mut a);
    let mut b = vec![0; 16];
    f_rng.fill_bytes(&mut b);
    let f = PointFn { a, b };

    c.bench_function("gen", |b| b.iter(|| gen(black_box(f.clone()))));

    let xs_n = 100;
    let mut xs = vec![f.a.clone()];
    for _ in 0..(xs_n - xs.len()) {
        let mut a = vec![0; 16];
        f_rng.fill_bytes(&mut a);
        xs.push(a);
    }

    let [share0, _] = gen(f);
    c.bench_function("eval_xs100", |b| {
        b.iter(|| eval_xs100(black_box(share0.clone()), black_box(xs.clone())))
    });
}

criterion_group!(benches, criterion_benches);
criterion_main!(benches);
