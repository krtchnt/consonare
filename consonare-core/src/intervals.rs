//! Step 8: Interval selection and naming
//!
//! Inputs: the smoothed dissonance profile and detected local minima from Step 7.
//! For each minimum (at r* ~ c* cents), generate rational candidates p/q within a
//! tolerance using a denominator cap, score with a weighted objective, pick best,
//! and assign a human-friendly label where possible.
//!
//! Notes (theory):
//! - Cents conversion: c = 1200 * log2(r) and r = 2^(c/1200).
//! - Complexity metric: Tenney height ~ log2(p) + log2(q) ~ ln(p) + ln(q).
//! - Local minima of the dissonance curve correlate with consonant intervals.
//!
//! Public API mirrors other steps: a Config, a Result, and `run_*` entry point.

use std::cmp::Ordering;

use crate::{
    common::{cents_to_ratio, ratio_to_cents},
    dissonance::{DissonancePoint, DissonanceResult},
};

/// Configuration for Step 8.
#[derive(Clone, Debug)]
pub struct IntervalNamingConfig {
    /// Candidate acceptance window (± cents around a local minimum).
    pub tolerance_cents: f32, // e.g., 5.0
    /// Maximum denominator to search for p/q candidates.
    pub max_denominator: u32, // e.g., 64
    /// Scoring weights: minimize J = α * D_norm + β * C + γ * |Δcents|.
    pub alpha: f32,
    pub beta: f32,
    pub gamma: f32,
    /// When true, compute a crude sharpness (2nd difference) around each min for info/confidence.
    pub estimate_sharpness: bool,
    /// Maximum number of alternative candidates (besides the best) to keep.
    pub max_alternatives: usize,
}

impl Default for IntervalNamingConfig {
    fn default() -> Self {
        Self {
            tolerance_cents: 5.0,
            max_denominator: 64,
            alpha: 1.0,
            beta: 0.5,
            gamma: 1.0,
            estimate_sharpness: true,
            max_alternatives: 5,
        }
    }
}

/// A candidate rational approximation near a local minimum.
#[derive(Clone, Debug)]
pub struct IntervalCandidate {
    pub p: u32,
    pub q: u32,
    pub ratio: f32,
    pub cents: f32,
    pub cents_error: f32, // cents(candidate) - cents(min)
    pub complexity: f32,  // ~ Tenney height (natural log base scaling)
    pub score: f32,       // lower is better
    pub label: String,
    pub confidence: f32, // 0..1 heuristic from depth, error, and sharpness
}

/// Annotated local minimum with its chosen label and alternatives.
#[derive(Clone, Debug)]
pub struct NamedMinimum {
    pub at_cents: f32,
    pub at_value: f32,
    pub depth: f32,
    pub sharpness: f32,
    pub best: Option<IntervalCandidate>,
    pub alternatives: Vec<IntervalCandidate>,
}

/// Final result for Step 8.
#[derive(Clone, Debug)]
pub struct IntervalNamingResult {
    pub named: Vec<NamedMinimum>,
}

/// Entry point for Step 8.
pub fn run_interval_naming_step(
    d: &DissonanceResult,
    cfg: &IntervalNamingConfig,
) -> IntervalNamingResult {
    // Precompute normalization for D: map [min(grid.value), max(grid.value)] -> [0,1]
    let (gmin, gmax) = min_max_values(&d.smoothed);
    let d_span = (gmax - gmin).max(1e-12);
    let norm_d = |v: f32| (v - gmin) / d_span;

    // Precompute a crude sharpness map via 2nd difference on smoothed grid
    let sharpness_series = if cfg.estimate_sharpness {
        second_diff_magnitude(&d.smoothed)
    } else {
        vec![0.0; d.smoothed.len()]
    };

    let mut named = Vec::with_capacity(d.minima.len());
    for m in &d.minima {
        // candidate generation near this minimum
        let mut candidates =
            generate_candidates_for_min(m.cents, cfg.tolerance_cents, cfg.max_denominator);
        // score and label
        for cand in &mut candidates {
            let d_norm = norm_d(m.value);
            let j = cfg.alpha * d_norm
                + cfg.beta * cand.complexity
                + cfg.gamma * cand.cents_error.abs();
            cand.score = j;
            // heuristic confidence (0..1)
            let s = estimate_sharpness_at(&d.smoothed, &sharpness_series, m.cents);
            let depth_factor = 1.0 - d_norm; // deeper valley => closer to 1
            let error_factor =
                (1.0 - (cand.cents_error.abs() / cfg.tolerance_cents)).clamp(0.0, 1.0);
            let sharp_factor = (s / (s + 1.0)).clamp(0.0, 1.0); // monotone map
            cand.confidence = 0.5 * depth_factor + 0.3 * error_factor + 0.2 * sharp_factor;
            // label
            cand.label = suggest_label(cand.p, cand.q);
        }

        // Rank by score (lower is better), ties by |Δcents| then complexity
        candidates.sort_by(|a, b| {
            a.score
                .partial_cmp(&b.score)
                .unwrap_or(Ordering::Equal)
                .then_with(|| {
                    a.cents_error
                        .abs()
                        .partial_cmp(&b.cents_error.abs())
                        .unwrap_or(Ordering::Equal)
                })
                .then_with(|| {
                    a.complexity
                        .partial_cmp(&b.complexity)
                        .unwrap_or(Ordering::Equal)
                })
        });

        let best = candidates.first().cloned();
        let alternatives = candidates
            .into_iter()
            .skip(1)
            .take(cfg.max_alternatives)
            .collect();

        // Compute sharpness once per minimum for reporting
        let sharp = estimate_sharpness_at(&d.smoothed, &sharpness_series, m.cents);

        named.push(NamedMinimum {
            at_cents: m.cents,
            at_value: m.value,
            depth: m.depth,
            sharpness: sharp,
            best,
            alternatives,
        });
    }

    IntervalNamingResult { named }
}

/// Generate all reduced fractions p/q with 1 <= q <= max_q such that the cents
/// error from target_cents is within tolerance.
fn generate_candidates_for_min(
    target_cents: f32,
    tol_cents: f32,
    max_q: u32,
) -> Vec<IntervalCandidate> {
    let target_ratio = cents_to_ratio(target_cents);
    let mut out = Vec::new();
    for q in 1..=max_q {
        let p_float = target_ratio * (q as f32);
        let p = p_float.round().max(1.0) as u32;
        let g = gcd(p, q);
        let (pn, qn) = (p / g, q / g);
        let ratio = (pn as f32) / (qn as f32);
        let cents = ratio_to_cents(ratio);
        let err = cents - target_cents;
        if err.abs() <= tol_cents {
            let complexity = (pn as f32).ln() + (qn as f32).ln(); // ~ Tenney height (natural log base)
            out.push(IntervalCandidate {
                p: pn,
                q: qn,
                ratio,
                cents,
                cents_error: err,
                complexity,
                score: 0.0, // fill later
                label: String::new(),
                confidence: 0.0,
            });
        }
    }
    // Deduplicate (p,q) in case multiple q produced same reduced form
    out.sort_by(|a, b| (a.p, a.q).cmp(&(b.p, b.q)));
    out.dedup_by(|a, b| a.p == b.p && a.q == b.q);
    out
}

/// GCD for u32
fn gcd(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// Return (min, max) over smoothed values.
fn min_max_values(grid: &[DissonancePoint]) -> (f32, f32) {
    let mut min_v = f32::INFINITY;
    let mut max_v = f32::NEG_INFINITY;
    for p in grid.iter() {
        if p.value < min_v {
            min_v = p.value;
        }
        if p.value > max_v {
            max_v = p.value;
        }
    }
    if !min_v.is_finite() || !max_v.is_finite() {
        (0.0, 1.0)
    } else {
        (min_v, max_v.max(min_v + 1e-9))
    }
}

/// Second difference magnitude series for a uniform cents grid (assumes near-uniform spacing).
fn second_diff_magnitude(grid: &[DissonancePoint]) -> Vec<f32> {
    if grid.len() < 3 {
        return vec![0.0; grid.len()];
    }
    let mut out = vec![0.0; grid.len()];
    for k in 1..grid.len() - 1 {
        let y_ = grid[k - 1].value - 2.0 * grid[k].value + grid[k + 1].value;
        out[k] = y_.abs();
    }
    out
}

/// Estimate sharpness near a cents location by sampling nearest index in the series.
fn estimate_sharpness_at(grid: &[DissonancePoint], sharp: &[f32], cents: f32) -> f32 {
    if grid.is_empty() || sharp.is_empty() {
        return 0.0;
    }
    // find nearest index (linear scan; grids are typically a few thousand points)
    let mut best_k = 0usize;
    let mut best_d = f32::INFINITY;
    for (k, p) in grid.iter().enumerate() {
        let d = (p.cents - cents).abs();
        if d < best_d {
            best_d = d;
            best_k = k;
        }
    }
    sharp[best_k]
}

/// Provide a friendly label for common small-ratio intervals. Fallback to "p:q".
fn suggest_label(p: u32, q: u32) -> String {
    const MAP: &[(u32, u32, &str)] = &[
        (1, 1, "unison"),
        (16, 15, "minor second"),
        (9, 8, "major second"),
        (6, 5, "minor third"),
        (5, 4, "major third"),
        (4, 3, "perfect fourth"),
        (7, 5, "septimal tritone"),
        (45, 32, "augmented fourth"),
        (64, 45, "diminished fifth"),
        (10, 7, "septimal tritone"),
        (3, 2, "perfect fifth"),
        (8, 5, "minor sixth"),
        (5, 3, "major sixth"),
        (7, 4, "harmonic seventh"),
        (16, 9, "small minor seventh"),
        (9, 5, "large minor seventh"),
        (15, 8, "major seventh"),
        (2, 1, "octave"),
    ];
    for (pn, qn, name) in MAP.iter() {
        if *pn == p && *qn == q {
            return format!("{}/{} ({})", p, q, name);
        }
    }
    format!("{}/{}", p, q)
}
