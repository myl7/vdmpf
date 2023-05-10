use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion};

use rand::prelude::*;

use vdmpf::dmpf::{MShare, VDMPF};
use vdmpf::dpf::PointFn;
use vdmpf::dyn_utils::*;

fn gen_t1k(fs: &[&PointFn]) -> MShare {
    let vdmpf = VDMPF::new(
        80f64,
        1000,
        Box::new(Aes128PRP::default()),
        Box::new(ChaChaISampler::default()),
        Box::new(ChaChaBSampler::default()),
        1000,
        16,
        Box::new(Aes128PRG::default()),
        Box::new(Aes128Hash::default()),
        Box::new(Aes128Hash::default()),
        Box::new(Aes128Hash::default()),
        Box::new(ChaChaBSampler::default()),
        1000,
    );
    let mshare = vdmpf.gen(fs).unwrap();
    mshare
}

fn eval_xs100k(mshare: &MShare, xs: &[&[u8]], ys: &mut [&mut [u8]], t: usize) -> Vec<u8> {
    let vdmpf = VDMPF::new(
        80f64,
        1000,
        Box::new(Aes128PRP::default()),
        Box::new(ChaChaISampler::default()),
        Box::new(ChaChaBSampler::default()),
        1000,
        16,
        Box::new(Aes128PRG::default()),
        Box::new(Aes128Hash::default()),
        Box::new(Aes128Hash::default()),
        Box::new(Aes128Hash::default()),
        Box::new(ChaChaBSampler::default()),
        1000,
    );
    let pi = vdmpf.eval(false, mshare, xs, ys, t);
    pi
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
        fs.push(PointFn { a, b });
    }
    let fs_ref = fs.iter().collect::<Vec<_>>();

    c.bench_function("gen_t1k", |b| b.iter(|| gen_t1k(black_box(&fs_ref))));

    let xs_n = 100_000;
    let mut xs = fs.iter().map(|f| f.a.clone()).collect::<Vec<_>>();
    for _ in 0..(xs_n - fs.len()) {
        let mut a = vec![0; 16];
        f_rng.fill_bytes(&mut a);
        a[0] = 0x0f;
        xs.push(a);
    }

    let mut mshare0 = gen_t1k(&fs_ref);
    mshare0
        .ks
        .iter_mut()
        .for_each(|k| k.s0s = vec![k.s0s[0].clone()]);
    let xs_ref = xs.iter().map(|x| x.as_ref()).collect::<Vec<_>>();
    let mut ys = vec![vec![0; 16]; xs.len()];
    let mut ys_ref = ys.iter_mut().map(|y| y.as_mut()).collect::<Vec<_>>();
    c.bench_function("eval_xs100k", |b| {
        b.iter(|| {
            eval_xs100k(
                black_box(&mshare0),
                black_box(&xs_ref),
                black_box(&mut ys_ref),
                black_box(fs.len()),
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
