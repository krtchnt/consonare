//! Step 7: Dissonance profile for two-note chords
//!
//! Given perceptually weighted partials S_A and S_B (Step 6), sweep interval
//! ratios r over a cents grid and compute total roughness
//!
//!     D(r) = sum_{i in A} sum_{j in B} w_i^A w_j^B psi(|f_j^B(r)-f_i^A|, f_m)
//!
//! where f_m = 0.5*(f_i^A + f_j^B(r)) and psi() is a pairwise roughness kernel.
//! We use the widely referenced Sethares/Plomp–Levelt-style kernel:
//!
//!     s = 0.24 / (0.021*f_m + 19)
//!     x = s * |Δf|
//!     psi = exp(-a*x) - exp(-b*x)     with a=3.5, b=5.75
//!
//! This places the roughness maximum near Δf ≈ 0.25*CBW, consistent with
//! Plomp & Levelt’s empirical curve. See:
//! - Sethares, “Local consonance and the relationship between timbre and scale” (1993/1998).
//! - Plomp & Levelt, “Tonal Consonance and Critical Bandwidth” (1965).
//!
//! If B is None, we synthesize B by scaling A (i.e., f_i^B(r) = r * f_i^A).
//!
//! The output includes a lightly smoothed profile and detected local minima.
//! Minima “depth” is estimated against neighboring shoulders to help rank them.
//!
//! API mirrors prior steps: a Config, a Result, and a single `run_*` entry point.

use crate::{
    common::{
        cbw_hz_plomp_levelt, cents_to_ratio, sethares_pair_roughness,
        sethares_pair_roughness_no_erb,
    },
    weighting::{WeightedPartial, WeightingResult},
};

#[derive(Clone, Debug)]
pub struct DissonanceConfig {
    /// Min and max cents to sweep (e.g., -1200..2400).
    pub min_cents: f32,
    pub max_cents: f32,
    /// Cents step (1–2 recommended).
    pub step_cents: f32,
    /// Moving-average window (odd; >=1) for light smoothing of D(r).
    pub smooth_window: usize,
    /// If Some(k), skip pair contributions when |Δf| > k * CBW(f_m).
    /// This speeds things up without affecting the curve shape much.
    pub max_deltaf_over_cbw: Option<f32>,
    /// Number of minima to keep (by increasing D value).
    pub top_k_minima: usize,
    /// If false, use a fixed roughness kernel (no ERB/CBW scaling).
    pub normalize_by_cbw: bool,
}

impl Default for DissonanceConfig {
    fn default() -> Self {
        Self {
            min_cents: -1200.0,
            max_cents: 2400.0,
            step_cents: 2.0,
            smooth_window: 5,
            max_deltaf_over_cbw: Some(2.0),
            top_k_minima: 16,
            normalize_by_cbw: true,
        }
    }
}

#[derive(Clone, Debug)]
pub struct DissonancePoint {
    pub cents: f32,
    pub value: f32,
}

#[derive(Clone, Debug)]
pub struct LocalMinimum {
    pub cents: f32,
    pub value: f32,
    /// Simple depth estimate against shoulders (higher = “stronger valley”).
    pub depth: f32,
}

#[derive(Clone, Debug)]
pub struct DissonanceResult {
    pub grid: Vec<DissonancePoint>,
    pub smoothed: Vec<DissonancePoint>,
    pub minima: Vec<LocalMinimum>,
}

/// Compute D(r) over a cents grid using the specified kernel.
pub fn run_dissonance_step(
    a: &WeightingResult,
    b_opt: Option<&WeightingResult>,
    cfg: &DissonanceConfig,
) -> DissonanceResult {
    let a_parts = &a.partials;
    let b_parts: &[WeightedPartial] = if let Some(b) = b_opt {
        &b.partials
    } else {
        // If B is absent, we'll reuse A but scaled by r each step;
        // keep a copy to avoid realloc each iteration.
        &a.partials
    };

    let mut grid: Vec<DissonancePoint> = Vec::new();
    if a_parts.is_empty() || b_parts.is_empty() {
        return DissonanceResult {
            grid,
            smoothed: Vec::new(),
            minima: Vec::new(),
        };
    }

    // Pre-extract weights and freqs for speed.
    let a_freqs: Vec<f32> = a_parts.iter().map(|p| p.freq_hz).collect();
    let a_ws: Vec<f32> = a_parts.iter().map(|p| p.weight).collect();
    let b_freqs0: Vec<f32> = b_parts.iter().map(|p| p.freq_hz).collect();
    let b_ws: Vec<f32> = b_parts.iter().map(|p| p.weight).collect();

    let mut c = cfg.min_cents;
    while c <= cfg.max_cents + 1e-6 {
        let r = cents_to_ratio(c);
        let mut d_total = 0.0_f32;

        // Scale B by r
        for (j, &f_b0) in b_freqs0.iter().enumerate() {
            let f_b = r * f_b0;
            let w_b = b_ws[j];
            if w_b <= 0.0 {
                continue;
            }

            for (i, &f_a) in a_freqs.iter().enumerate() {
                let w = w_b * a_ws[i];
                if w <= 0.0 {
                    continue;
                }

                let f_m = 0.5 * (f_a + f_b);
                let delta = (f_b - f_a).abs();

                if let Some(k) = cfg.max_deltaf_over_cbw {
                    let cbw = cbw_hz_plomp_levelt(f_m);
                    if delta > k * cbw {
                        continue;
                    }
                }

                let pr = if cfg.normalize_by_cbw {
                    sethares_pair_roughness(delta, f_m)
                } else {
                    sethares_pair_roughness_no_erb(delta)
                };
                d_total += w * pr.max(0.0);
            }
        }

        grid.push(DissonancePoint {
            cents: c,
            value: d_total,
        });
        c += cfg.step_cents;
    }

    // Light smoothing to suppress numerical roughness
    let smoothed = moving_average(&grid, cfg.smooth_window);

    // Locate local minima on smoothed curve
    let minima = find_local_minima(&smoothed, cfg.top_k_minima);

    DissonanceResult {
        grid,
        smoothed,
        minima,
    }
}

fn moving_average(xs: &[DissonancePoint], win: usize) -> Vec<DissonancePoint> {
    if xs.is_empty() || win <= 1 || win.is_multiple_of(2) {
        return xs.to_vec();
    }
    let n = xs.len();
    let h = win / 2;
    let mut out = Vec::with_capacity(n);
    let mut acc = 0.0_f32;
    // initialize first window
    for xs_i in xs.iter().take(win.min(n)) {
        acc += xs_i.value;
    }
    for i in 0..n {
        let l = i.saturating_sub(h);
        let r = (i + h).min(n - 1);
        // recompute on edges; otherwise slide
        let v = if i == 0 || i + h >= n || win > (r - l + 1) {
            let mut sum = 0.0;
            let mut cnt = 0.0;
            for xs_k in xs.iter().take(r + 1).skip(l) {
                sum += xs_k.value;
                cnt += 1.0;
            }
            sum / cnt
        } else {
            acc / (win as f32)
        };
        out.push(DissonancePoint {
            cents: xs[i].cents,
            value: v,
        });

        // slide accumulator if fully inside bounds
        if i + h + 1 < n && i >= h {
            acc += xs[i + h + 1].value;
            acc -= xs[i - h].value;
        }
    }
    out
}

fn find_local_minima(xs: &[DissonancePoint], top_k: usize) -> Vec<LocalMinimum> {
    let n = xs.len();
    if n < 3 {
        return Vec::new();
    }
    let mut mins: Vec<LocalMinimum> = Vec::new();

    // First pass: collect candidate minima
    for i in 1..n - 1 {
        let (l, m, r) = (xs[i - 1].value, xs[i].value, xs[i + 1].value);
        if m <= l && m <= r && (m < l || m < r) {
            // Estimate “depth” against nearest shoulders (local maxima)
            let left_peak = find_left_peak(xs, i);
            let right_peak = find_right_peak(xs, i);
            let shoulder = 0.5 * (left_peak + right_peak);
            let depth = (shoulder - m).max(0.0);
            mins.push(LocalMinimum {
                cents: xs[i].cents,
                value: m,
                depth,
            });
        }
    }

    // Rank by value ascending (deeper minima first); tie-break by depth
    mins.sort_by(|a, b| {
        a.value
            .partial_cmp(&b.value)
            .unwrap()
            .then(b.depth.partial_cmp(&a.depth).unwrap())
    });
    if mins.len() > top_k {
        mins.truncate(top_k);
    }
    mins
}

fn find_left_peak(xs: &[DissonancePoint], idx: usize) -> f32 {
    let mut i = idx;
    let mut best = xs[idx].value;
    while i > 0 {
        if xs[i - 1].value > best {
            best = xs[i - 1].value;
        } else {
            break;
        }
        i -= 1;
    }
    best
}

fn find_right_peak(xs: &[DissonancePoint], idx: usize) -> f32 {
    let mut i = idx;
    let mut best = xs[idx].value;
    while i + 1 < xs.len() {
        if xs[i + 1].value > best {
            best = xs[i + 1].value;
        } else {
            break;
        }
        i += 1;
    }
    best
}
