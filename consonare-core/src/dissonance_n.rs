//! Step 9: Extension to N-note chords (multi-note roughness minimisation)
//!
//! We generalise Step 7 (two-note dissonance profile) to chords with N notes
//! that share the *same timbre* (so every note's spectrum is a scaled copy of
//! a reference note’s weighted spectrum).
//!
//! Approach (runtime-friendly):
//! 1) Pick a reference note (index 0, fixed at ratio r_0 = 1).
//! 2) For each other note k ∈ {1..N-1}, compute the two-note dissonance profile
//!    against the reference using Step 7 and extract its top-K local minima as
//!    *seeds* for r_k (in cents).
//! 3) Perform a *beam search* over these per-note candidate sets (width B)
//!    expanding one note at a time. At partial depth d we score the partial
//!    assignment using the total roughness over the d+1 notes chosen so far.
//!    This avoids the K^(N-1) combinatorial explosion while exploring good
//!    combinations. B and K are small (e.g. 5) so the search is fast.
//! 4) Locally refine each full assignment via 1D coordinate descent around each
//!    r_k (small ±window) while holding others fixed.
//! 5) Project each ratio to a nearby simple rational p/q (continued fractions)
//!    within a tolerance in cents and provide a readable label.
//!
//! Complexity: O(B * (N-1) * cost_eval), with cost_eval ≈ O(P^2 * M^2), where
//! M is the number of notes in the *partial* chord and P is partials per note.
//! This scales linearly in N (for fixed B, K) rather than exponentially.
//!
//! The pairwise roughness kernel and pruning are identical to Step 7.
//!
//! Public API mirrors other steps: a Config, a Result, and a `run_*` entry point.

use std::cmp::Ordering;

use crate::{
    // helper functions are below
    common::{cbw_hz_plomp_levelt, cents_to_ratio, ratio_to_cents, sethares_pair_roughness},
    dissonance::{DissonanceConfig, run_dissonance_step},
    weighting::WeightingResult,
};

/// Configuration for Step 9 multi-note search.
#[derive(Clone, Debug)]
pub struct MultiDissonanceConfig {
    /// Number of total notes in the chord, including the fixed reference.
    pub n_notes: usize,
    /// Sweep bounds used for per-note 2-note seeding (cents).
    pub min_cents: f32,
    pub max_cents: f32,
    pub step_cents: f32,
    /// Max |Δf|/CBW pruning (as in Step 7).
    pub max_deltaf_over_cbw: Option<f32>,
    /// How many local minima to keep per coordinate when seeding.
    pub seeds_per_note: usize,
    /// Beam width for combinatorial search (larger explores more).
    pub beam_width: usize,
    /// Local coordinate-descent refinement window (± cents).
    pub refine_window_cents: f32,
    /// Refinement step size in cents.
    pub refine_step_cents: f32,
    /// Number of refinement sweeps.
    pub refine_sweeps: usize,
    /// Max denominator for p/q projection.
    pub approx_max_den: u32,
    /// Accept p/q only if within this cents tolerance.
    pub approx_tol_cents: f32,
    /// Return top-K multi-note solutions (by total roughness).
    pub top_k_solutions: usize,
}

impl Default for MultiDissonanceConfig {
    fn default() -> Self {
        Self {
            n_notes: 3,
            min_cents: -1200.0,
            max_cents: 2400.0,
            step_cents: 2.0,
            max_deltaf_over_cbw: Some(3.0),
            seeds_per_note: 5,
            beam_width: 5,
            refine_window_cents: 6.0,
            refine_step_cents: 1.0,
            refine_sweeps: 2,
            approx_max_den: 32,
            approx_tol_cents: 6.0,
            top_k_solutions: 6,
        }
    }
}

/// A single ratio relative to the reference, with an optional rational name.
#[derive(Clone, Debug)]
pub struct NamedRatio {
    pub ratio: f32,
    pub cents: f32,
    pub approx: Option<ApproxPQ>,
}

/// Rational approximation and naming, à la Step 8.
#[derive(Clone, Debug)]
pub struct ApproxPQ {
    pub p: u32,
    pub q: u32,
    pub label: String,
    pub cents_error: f32,
    pub complexity: f32,
}

/// A full chord hypothesis.
#[derive(Clone, Debug)]
pub struct MultiSolution {
    /// Ratios r_0..r_{N-1} with r_0 fixed to 1.0.
    pub ratios: Vec<NamedRatio>,
    /// Total pairwise roughness over all note pairs (smaller is better).
    pub total_roughness: f32,
}

/// Output for Step 9.
#[derive(Clone, Debug)]
pub struct MultiDissonanceResult {
    pub solutions: Vec<MultiSolution>,
}

/// Entry point for Step 9.
///
/// `a_ref` is the reference note’s weighted spectrum (Step 6). If you have
/// multiple *measured* weighted spectra for the other notes with different
/// timbres, pass them in `others`; otherwise pass `None` to synthesise notes
/// by scaling `a_ref` (same timbre case).
pub fn run_multi_dissonance_step(
    a_ref: &WeightingResult,
    others: Option<&[WeightingResult]>,
    cfg: &MultiDissonanceConfig,
) -> MultiDissonanceResult {
    let n = cfg.n_notes.max(2);
    // 1) Seeds from two-note dissonance minima against the reference
    let seed_cfg = DissonanceConfig {
        min_cents: cfg.min_cents,
        max_cents: cfg.max_cents,
        step_cents: cfg.step_cents,
        smooth_window: 5,
        max_deltaf_over_cbw: cfg.max_deltaf_over_cbw,
        top_k_minima: cfg.seeds_per_note,
        ..Default::default()
    };
    let seeds_per_coord: Vec<Vec<f32>> = (1..n)
        .map(|k| {
            // Use the other note's spectrum if provided; else synthesise from A.
            let b_opt = others.and_then(|xs| xs.get(k - 1));
            let d2 = run_dissonance_step(a_ref, b_opt, &seed_cfg);
            d2.minima.iter().map(|m| m.cents).collect::<Vec<f32>>()
        })
        .collect();

    // 2) Beam search over seeds (coordinates are cents r_k)
    #[derive(Clone)]
    struct BeamState {
        cents: Vec<f32>, // length d (number of chosen non-ref notes)
        score: f32,      // total roughness on partial chord [ref + chosen]
    }

    let mut beam: Vec<BeamState> = vec![BeamState {
        cents: vec![],
        score: 0.0,
    }];

    for coord_seeds in seeds_per_coord.iter() {
        let mut next: Vec<BeamState> = Vec::new();
        for st in &beam {
            let last_key = st.cents.last().map(|c| round_key(*c));
            for &c_k in coord_seeds.iter() {
                // Enforce nondecreasing order in rounded cents to avoid permutations.
                if let Some(lk) = last_key
                    && round_key(c_k) < lk
                {
                    continue;
                }
                let mut cents = st.cents.clone();
                cents.push(c_k);
                // Score partial assignment using roughness among currently-chosen notes.
                let score =
                    total_roughness_for_cents(a_ref, others, &cents, cfg.max_deltaf_over_cbw);
                next.push(BeamState { cents, score });
            }
        }
        // Keep the best `beam_width` states
        next.sort_by(|a, b| a.score.partial_cmp(&b.score).unwrap_or(Ordering::Equal));
        next.truncate(cfg.beam_width.max(1));
        beam = next;
    }

    // 3) Local coordinate-descent refinement around each completed assignment
    let mut full_solutions: Vec<(Vec<f32>, f32)> = Vec::new();
    for st in beam {
        let mut cents_vec = st.cents.clone();
        for _ in 0..cfg.refine_sweeps {
            for i in 0..cents_vec.len() {
                let base = cents_vec[i];
                let mut best_c = base;
                let mut best_s =
                    total_roughness_for_cents(a_ref, others, &cents_vec, cfg.max_deltaf_over_cbw);
                let step = cfg.refine_step_cents.max(0.5);
                let w = cfg.refine_window_cents.max(step);
                let mut c = base - w;
                while c <= base + w + 1e-6 {
                    cents_vec[i] = c;
                    let s = total_roughness_for_cents(
                        a_ref,
                        others,
                        &cents_vec,
                        cfg.max_deltaf_over_cbw,
                    );
                    if s < best_s {
                        best_s = s;
                        best_c = c;
                    }
                    c += step;
                }
                cents_vec[i] = best_c;
            }
        }
        let final_score =
            total_roughness_for_cents(a_ref, others, &cents_vec, cfg.max_deltaf_over_cbw);
        full_solutions.push((cents_vec, final_score));
    }

    // 4) Sort, de-duplicate near-identical solutions (order-invariant), and label.
    full_solutions.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));

    let mut uniq: Vec<(Vec<f32>, f32)> = Vec::new();
    'outer: for (cs, s) in full_solutions {
        // Deduplicate by sorted rounded cents (order-invariant).
        let mut key: Vec<i32> = cs.iter().map(|c| round_key(*c)).collect();
        key.sort_unstable();

        for (ocs, _) in &uniq {
            let mut okey: Vec<i32> = ocs.iter().map(|c| round_key(*c)).collect();
            okey.sort_unstable();
            if key == okey {
                continue 'outer;
            }
        }
        uniq.push((cs, s));
        if uniq.len() >= cfg.top_k_solutions {
            break;
        }
    }

    let mut solutions: Vec<MultiSolution> = Vec::new();
    for (cents_nonref, score) in uniq {
        // Canonicalize order for display: sort non-ref cents ascending.
        let mut cents_sorted = cents_nonref.clone();
        cents_sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));

        // Assemble r_0 = 1 plus non-ref ratios
        let mut ratios = Vec::with_capacity(cents_sorted.len() + 1);
        ratios.push(NamedRatio {
            ratio: 1.0,
            cents: 0.0,
            approx: Some(ApproxPQ {
                p: 1,
                q: 1,
                label: "unison".to_string(),
                cents_error: 0.0,
                complexity: 0.0,
            }),
        });
        for &c in &cents_sorted {
            let r = cents_to_ratio(c);
            let approx = approximate_ratio_with_label(r, cfg.approx_max_den, cfg.approx_tol_cents);
            ratios.push(NamedRatio {
                ratio: r,
                cents: c,
                approx,
            });
        }
        solutions.push(MultiSolution {
            ratios,
            total_roughness: score,
        });
    }

    MultiDissonanceResult { solutions }
}

/// Compute total roughness for the chord specified by *non-reference* cents
/// values (reference is fixed at 0 cents). If `others` is provided, each
/// non-reference note uses its measured spectrum; otherwise we scale `a_ref`.
///
/// Quantize cents for order-invariant comparisons and pruning.
#[inline]
fn round_key(cents: f32) -> i32 {
    // Round to nearest 0.5 cent to suppress tiny jitter.
    (2.0 * cents).round() as i32
}

pub fn total_roughness_for_cents(
    a_ref: &WeightingResult,
    others: Option<&[WeightingResult]>,
    cents_nonref: &[f32],
    max_deltaf_over_cbw: Option<f32>,
) -> f32 {
    // Build spectra for all notes (ratios r_0..r_{M-1})
    let mut notes: Vec<Vec<(f32, f32)>> = Vec::new(); // (freq_hz, weight)
    notes.push(
        a_ref
            .partials
            .iter()
            .map(|p| (p.freq_hz, p.weight))
            .collect(),
    );
    for (k, &c) in cents_nonref.iter().enumerate() {
        let r = cents_to_ratio(c);
        match others.and_then(|xs| xs.get(k)) {
            Some(wb) => {
                notes.push(
                    wb.partials
                        .iter()
                        .map(|p| (r * p.freq_hz, p.weight))
                        .collect(),
                );
            }
            None => {
                notes.push(
                    a_ref
                        .partials
                        .iter()
                        .map(|p| (r * p.freq_hz, p.weight))
                        .collect(),
                );
            }
        }
    }
    // Sum pairwise roughness over all a < b
    let mut total = 0.0f32;
    for a in 0..notes.len() {
        for b in (a + 1)..notes.len() {
            let sa = &notes[a];
            let sb = &notes[b];
            for &(fa, wa) in sa.iter() {
                for &(fb, wb) in sb.iter() {
                    let fm = 0.5 * (fa + fb);
                    let df = (fb - fa).abs();
                    if let Some(kmax) = max_deltaf_over_cbw {
                        let cbw = cbw_hz_plomp_levelt(fm);
                        if df > kmax * cbw {
                            continue;
                        }
                    }
                    let psi = sethares_pair_roughness(df, fm);
                    total += wa * wb * psi;
                }
            }
        }
    }
    total
}

/// Continued fractions p/q approximation of a ratio r with denominator cap.
fn best_rationals_near(r: f32, max_den: u32) -> Vec<(u32, u32)> {
    // Generate convergents up to max_den, include mediants near the end.
    let mut res: Vec<(u32, u32)> = Vec::new();
    let mut x = r as f64;
    let mut a0 = x.floor();
    let mut p0: u128 = 1;
    let mut q0: u128 = 0;
    let mut p1: u128 = a0 as u128;
    let mut q1: u128 = 1;
    res.push((p1 as u32, q1 as u32));
    let mut iter = 0;
    while q1 <= max_den as u128 && iter < 32 {
        iter += 1;
        x = 1.0 / (x - a0);
        a0 = x.floor();
        let p2 = (a0 as u128).saturating_mul(p1).saturating_add(p0);
        let q2 = (a0 as u128).saturating_mul(q1).saturating_add(q0);
        if q2 == 0 {
            break;
        }
        // Stop if denominator would exceed cap (no push), or if values got too large to cast.
        if q2 <= max_den as u128 && p2 <= u32::MAX as u128 {
            res.push((p2 as u32, q2 as u32));
        } else {
            break;
        }
        p0 = p1;
        q0 = q1;
        p1 = p2;
        q1 = q2;
    }
    // Add a few nearby mediants to densify candidates
    let mut extras: Vec<(u32, u32)> = Vec::new();
    for i in 0..res.len().saturating_sub(1) {
        let (p_a, q_a) = res[i];
        let (p_b, q_b) = res[i + 1];
        let p_m = (p_a as u64) + (p_b as u64);
        let q_m = (q_a as u64) + (q_b as u64);
        if q_m <= max_den as u64 && p_m <= u32::MAX as u64 {
            extras.push((p_m as u32, q_m as u32));
        }
    }
    res.extend(extras);
    // Unify
    res.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    res.dedup();
    res
}

fn suggest_label(p: u32, q: u32) -> &'static str {
    // A small map consistent with Step 8.
    const MAP: &[(u32, u32, &str)] = &[
        (1, 1, "unison"),
        (16, 15, "diatonic semitone"),
        (15, 14, "diatonic semitone"),
        (9, 8, "major second"),
        (6, 5, "minor third"),
        (5, 4, "major third"),
        (4, 3, "perfect fourth"),
        (7, 5, "septimal tritone"),
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
            return name;
        }
    }
    ""
}

fn approximate_ratio_with_label(r: f32, max_den: u32, tol_cents: f32) -> Option<ApproxPQ> {
    let cents = ratio_to_cents(r);
    let mut best: Option<ApproxPQ> = None;
    for (p, q) in best_rationals_near(r, max_den) {
        if q == 0 {
            continue;
        }
        let rr = (p as f32) / (q as f32);
        let c = ratio_to_cents(rr);
        let err = (c - cents).abs();
        if err <= tol_cents {
            let label = suggest_label(p, q);
            let complexity = (p as f32).ln() + (q as f32).ln();
            let apq = ApproxPQ {
                p,
                q,
                label: if label.is_empty() {
                    format!("{}/{}", p, q)
                } else {
                    format!("{}/{} ({})", p, q, label)
                },
                cents_error: err,
                complexity,
            };
            match &best {
                None => best = Some(apq),
                Some(b) => {
                    // Prefer lower complexity; tie-break by |err|
                    if apq.complexity < b.complexity - 1e-6
                        || ((apq.complexity - b.complexity).abs() < 1e-6
                            && apq.cents_error < b.cents_error)
                    {
                        best = Some(apq);
                    }
                }
            }
        }
    }
    best
}
