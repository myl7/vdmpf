use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion};

use rand::prelude::*;

use vdmpf::dmpf::{MShare, VDMPF};
use vdmpf::dpf::{Gen, PointFn};
use vdmpf::dyn_utils::*;

fn gen_t1k(fs: Vec<PointFn>, prg: Box<dyn Gen>) -> [MShare; 2] {
    let vdmpf = VDMPF::new(
        80f64,
        1000,
        Box::new(Aes256PRP::default()),
        Box::new(ChaChaISampler::default()),
        Box::new(ChaChaBSampler::default()),
        1000,
        16,
        prg,
        Box::new(Shake256Hash::default()),
        Box::new(Shake256Hash::default()),
        Box::new(Shake256Hash::default()),
        Box::new(ChaChaBSampler::default()),
        1000,
    );
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
    [mshare0, mshare1]
}

fn eval_xs100k(
    mshare: MShare,
    xs: Vec<Vec<u8>>,
    t: usize,
    prg: Box<dyn Gen>,
) -> (Vec<Vec<u8>>, Vec<u8>) {
    let vdmpf = VDMPF::new(
        80f64,
        1000,
        Box::new(Aes256PRP::default()),
        Box::new(ChaChaISampler::default()),
        Box::new(ChaChaBSampler::default()),
        1000,
        16,
        prg,
        Box::new(Shake256Hash::default()),
        Box::new(Shake256Hash::default()),
        Box::new(Shake256Hash::default()),
        Box::new(ChaChaBSampler::default()),
        1000,
    );
    let xs_ref = xs.iter().map(|x| x.as_ref()).collect::<Vec<_>>();
    let (y0s, pi0) = vdmpf.eval(false, &mshare, &xs_ref, t);
    (y0s, pi0)
}

fn criterion_benches(c: &mut Criterion) {
    // Multi-point function with 1000 points
    let point_n = 1000;
    let mut fs = Vec::with_capacity(point_n);
    let mut f_rng = thread_rng();
    for _ in 0..point_n {
        let mut a = vec![0; 16];
        f_rng.fill_bytes(&mut a);
        a[0] = 0x07;
        let mut b = vec![0; 16];
        f_rng.fill_bytes(&mut b);
        fs.push(PointFn { a, b, a_leap: None });
    }

    c.bench_function("aes128prg_gen_t1k", |b| {
        b.iter(|| {
            gen_t1k(
                black_box(fs.clone()),
                black_box(Box::new(Aes128PRG::default())),
            )
        })
    });
    c.bench_function("chachaprg_gen_t1k", |b| {
        b.iter(|| {
            gen_t1k(
                black_box(fs.clone()),
                black_box(Box::new(ChaChaPRG::default())),
            )
        })
    });

    let xs_n = 100_000;
    let mut xs = fs.iter().map(|f| f.a.clone()).collect::<Vec<_>>();
    for _ in 0..(xs_n - fs.len()) {
        let mut a = vec![0; 16];
        f_rng.fill_bytes(&mut a);
        a[0] = 0x0f;
        xs.push(a);
    }

    let [mshare0_aes128, _] = gen_t1k(fs.clone(), Box::new(Aes128PRG::default()));
    c.bench_function("aes128prg_eval_xs100k", |b| {
        b.iter(|| {
            eval_xs100k(
                black_box(mshare0_aes128.clone()),
                black_box(xs.clone()),
                black_box(fs.len()),
                black_box(Box::new(Aes128PRG::default())),
            )
        })
    });

    let [mshare0_chacha, _] = gen_t1k(fs.clone(), Box::new(ChaChaPRG::default()));
    c.bench_function("chachaprg_eval_xs100k", |b| {
        b.iter(|| {
            eval_xs100k(
                black_box(mshare0_chacha.clone()),
                black_box(xs.clone()),
                black_box(fs.len()),
                black_box(Box::new(ChaChaPRG::default())),
            )
        })
    });
}

criterion_group! {
    name = benches;
    config = Criterion::default().sample_size(10).measurement_time(Duration::from_secs(90));
    targets = criterion_benches
}
criterion_main!(benches);
