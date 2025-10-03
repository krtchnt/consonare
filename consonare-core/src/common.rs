#[inline]
pub fn cents_to_ratio(c: f32) -> f32 {
    (c / 1200.0).exp2() // 2^(c/1200)
}

#[inline]
pub fn ratio_to_cents(r: f32) -> f32 {
    1200.0 * r.log2()
}

#[inline]
pub fn cbw_hz_plomp_levelt(f_hz: f32) -> f32 {
    // A common analytical fit for critical bandwidth around f (Hz).
    // This specific shape isn't too critical when used only as a pruning limit.
    // Here we reuse the well-known Glasberg–Moore ERB conversion as a proxy.
    // (It scales similarly in the range we care about and is cheap to compute.)
    let fk = f_hz / 1000.0;
    24.7 * (4.37 * fk + 1.0) // Hz
}

#[inline]
pub fn sethares_pair_roughness(delta_f: f32, f_mean: f32) -> f32 {
    // Sethares/Plomp–Levelt kernel parameters
    let a = 3.5_f32;
    let b = 5.75_f32;
    // Critical-band scaling term
    let s = 0.24_f32 / (0.021_f32 * f_mean + 19.0_f32);
    let x = s * delta_f.abs();
    (-(a * x)).exp() - (-(b * x)).exp()
}

#[inline]
pub fn sethares_pair_roughness_no_erb(delta_f: f32) -> f32 {
    // Same a,b as usual, but no ERB/CBW scaling (s = 1).
    let a = 3.5_f32;
    let b = 5.75_f32;
    let x = delta_f.abs();
    (-(a * x)).exp() - (-(b * x)).exp()
}

/// Compute median of a list.
pub fn median(xs: &mut [f32]) -> f32 {
    if xs.is_empty() {
        return 0.0;
    }
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = xs.len() / 2;
    if xs.len() % 2 == 1 {
        xs[m]
    } else {
        0.5 * (xs[m - 1] + xs[m])
    }
}

/// Cents difference between two positive frequencies.#
#[inline]
pub fn cents_diff(f1: f32, f2: f32) -> f32 {
    if f1 <= 0.0 || f2 <= 0.0 {
        return f32::INFINITY;
    }
    1200.0 * ((f2 / f1) as f64).log2() as f32
}
