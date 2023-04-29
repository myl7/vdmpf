use criterion::{black_box, criterion_group, criterion_main, Criterion};

use hex_literal::hex;
use rand::prelude::*;

use vdmpf::dmpf::VDMPF;
use vdmpf::dpf::PointFn;
use vdmpf::dyn_utils::{Aes256PRP, ChaChaBSampler, ChaChaISampler, ChaChaPRG, Shake256Hash};
use vdmpf::group::BGroup;

fn var_t(point_n: usize) {
    let vdmpf = VDMPF::new(
        80f64,
        1000,
        Box::new(Aes256PRP::default()),
        Box::new(ChaChaISampler::default()),
        Box::new(ChaChaBSampler::default()),
        1000,
        16,
        Box::new(ChaChaPRG::default()),
        Box::new(Shake256Hash::default()),
        Box::new(Shake256Hash::default()),
        Box::new(Shake256Hash::default()),
        Box::new(ChaChaBSampler::default()),
        1000,
    );
    let mut fs = Vec::with_capacity(point_n);
    let mut f_rng = thread_rng();
    for _ in 0..point_n {
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
    for _ in 0..point_n {
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
}

fn criterion_bench(c: &mut Criterion) {
    c.bench_function("t 1000", |b| b.iter(|| var_t(black_box(1000))));
}

criterion_group! {
    name = benches;
    // This can be any expression that returns a `Criterion` object.
    config = Criterion::default().sample_size(10);
    targets = criterion_bench
}
criterion_main!(benches);
