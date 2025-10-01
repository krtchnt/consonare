//! Step 10: Quality control and diagnostics
//!
//! Flags low-confidence results and prints a compact diagnostic report.
//!
//! Heuristics (configurable):
//! - Few stable partials: N_stable < min_partials
//! - Poor f0 agreement: median per-frame method disagreement (MAD, cents) > max_f0_mad_cents
//! - Flat dissonance: max(minimum depth) / dynamic_range < min_norm_min_depth
//!
//! Additionally prints:
//! - Extracted partial summaries (median Hz/dB/SNR, track counts)
//! - Inharmonicity coefficient B if available
//! - Top named intervals table with p/q, cents, cents error, and confidence
//!
//! References for metrics/thresholds:
//! - YIN robustness and pitch confidence via method agreement (de Cheveigné & Kawahara, 2002).
//! - Dissonance profile minima depth as a proxy for consonance salience (Sethares, 1993; Vassilakis, 2001).

use crate::common::{cents_diff, median};
use crate::dissonance::DissonanceResult;
use crate::f0::F0Result;
use crate::inharmonicity::InharmonicityResult;
use crate::intervals::{IntervalCandidate, IntervalNamingResult};
use crate::partials::PartialsResult;

/// Configuration for Step 10 quality checks.
#[derive(Clone, Debug)]
pub struct QualityConfig {
    /// Minimum count of stable partials required.
    pub min_partials: usize,
    /// Maximum acceptable median absolute disagreement (cents) between f0 methods per frame.
    pub max_f0_mad_cents: f32,
    /// Minimum acceptable normalized minimum depth (max depth / (max-min)) for dissonance curve.
    pub min_norm_min_depth: f32,
    /// Only count partials with median SNR (dB) above this value.
    pub min_partial_snr_db: f32,
    /// Limit how many partials to print in the report (0 = all).
    pub max_partials_to_print: usize,
}

impl Default for QualityConfig {
    fn default() -> Self {
        Self {
            min_partials: 6,
            max_f0_mad_cents: 5.0,
            min_norm_min_depth: 0.05,
            min_partial_snr_db: 0.0,
            max_partials_to_print: 20,
        }
    }
}

/// Flags indicating potential quality issues.
#[derive(Clone, Debug, Default)]
pub struct QualityFlags {
    pub few_partials: bool,
    pub poor_f0_agreement: bool,
    pub flat_dissonance_curve: bool,
}

/// Aggregated diagnostics.
#[derive(Clone, Debug, Default)]
pub struct QualityResult {
    pub flags: QualityFlags,
    pub n_partials_kept: usize,
    pub f0_mad_cents: f32,
    pub dissonance_norm_depth: f32,
}

/// Run Step 10.
pub fn run_quality_step(
    parts: &PartialsResult,
    f0: &F0Result,
    diss: &DissonanceResult,
    names: &IntervalNamingResult,
    inh: &InharmonicityResult,
    cfg: &QualityConfig,
) -> QualityResult {
    // 1) Count stable partials above SNR threshold
    let n_partials = parts
        .summaries
        .iter()
        .filter(|s| s.median_snr_db >= cfg.min_partial_snr_db)
        .count();

    // 2) f0 method agreement MAD (cents), aggregated over frames
    let mut per_frame_mad: Vec<f32> = Vec::new();
    for fr in &f0.frames {
        let mut vals_hz: Vec<f32> = Vec::with_capacity(3);
        if let Some(y) = fr.yin_hz
            && y > 0.0
        {
            vals_hz.push(y);
        }
        if let Some(c) = fr.cep_hz
            && c > 0.0
        {
            vals_hz.push(c);
        }
        if let Some(s) = fr.spec_hz
            && s > 0.0
        {
            vals_hz.push(s);
        }
        if vals_hz.len() >= 2 {
            let med_hz = median(&mut vals_hz.clone());
            let mut cents_abs: Vec<f32> = vals_hz
                .iter()
                .map(|v| cents_diff(med_hz, *v).abs())
                .collect();
            per_frame_mad.push(median(&mut cents_abs));
        }
    }
    let f0_mad_cents = if per_frame_mad.is_empty() {
        0.0
    } else {
        median(&mut per_frame_mad.clone())
    };

    // 3) Dissonance minima normalized depth
    let mut norm_depth = 0.0;
    if !diss.smoothed.is_empty() {
        let min_v = diss
            .smoothed
            .iter()
            .map(|p| p.value)
            .fold(f32::INFINITY, f32::min);
        let max_v = diss
            .smoothed
            .iter()
            .map(|p| p.value)
            .fold(f32::NEG_INFINITY, f32::max);
        let dr = (max_v - min_v).max(1e-9);
        let max_depth = diss.minima.iter().map(|m| m.depth).fold(0.0, f32::max);
        norm_depth = (max_depth / dr).clamp(0.0, 1.0);
    }

    let flags = QualityFlags {
        few_partials: n_partials < cfg.min_partials,
        poor_f0_agreement: f0_mad_cents > cfg.max_f0_mad_cents,
        flat_dissonance_curve: norm_depth < cfg.min_norm_min_depth,
    };

    // 4) Print diagnostics report
    println!("\n== Step 10: Quality & Diagnostics ==");
    println!(
        "Stable partials kept: {} (SNR≥{:.1} dB){}",
        n_partials,
        cfg.min_partial_snr_db,
        if flags.few_partials { "  [LOW]" } else { "" }
    );
    println!(
        "f0 method disagreement (MAD): {:.2} cents{}",
        f0_mad_cents,
        if flags.poor_f0_agreement {
            "  [HIGH]"
        } else {
            ""
        }
    );
    if diss.smoothed.is_empty() {
        println!("Dissonance profile: (empty)");
    } else {
        println!(
            "Dissonance curve normalized max-min depth: {:.3}{}",
            norm_depth,
            if flags.flat_dissonance_curve {
                "  [FLAT]"
            } else {
                ""
            }
        );
    }

    if let Some(b) = inh.b {
        println!(
            "Inharmonicity: B ≈ {:.6e}  (R²={:.3}, used {})",
            b, inh.r2, inh.used_points
        );
        if let Some(med_c) = inh.median_abs_cents {
            println!("Model median abs cents error ≈ {:.2} cents", med_c);
        }
    } else {
        println!("Inharmonicity: no reliable stiff-string fit (keeping empirical partials).");
    }

    // Partials table
    if !parts.summaries.is_empty() {
        println!("Extracted partials (median over time):");
        println!(
            "  {:>4}  {:>10}  {:>8}  {:>8}  {:>6}",
            "id", "Hz", "dB", "SNRdB", "count"
        );
        for (i, s) in parts
            .summaries
            .iter()
            .take(cfg.max_partials_to_print)
            .enumerate()
        {
            println!(
                "  {:>4}  {:>10.3}  {:>8.2}  {:>8.2}  {:>6}",
                i, s.median_hz, s.median_db, s.median_snr_db, s.count
            );
        }
        if cfg.max_partials_to_print > 0 && parts.summaries.len() > cfg.max_partials_to_print {
            println!(
                "  ... ({} more)",
                parts.summaries.len() - cfg.max_partials_to_print
            );
        }
    }

    // Intervals table
    if !names.named.is_empty() {
        println!("Recommended intervals:");
        println!(
            "  {:>7}  {:>10}  {:>9}  {:>9}  {:>10}  label [conf]",
            "cents", "D", "depth", "sharpness", "p/q"
        );
        for nm in &names.named {
            let pq_label = match &nm.best {
                Some(IntervalCandidate { p, q, label, .. }) => format!("{}/{}  {}", p, q, label),
                None => "-".to_string(),
            };
            let conf = match &nm.best {
                Some(c) => format!("{:.2}", c.confidence),
                None => "-".to_string(),
            };
            println!(
                "  {:>7.2}  {:>10.6}  {:>9.6}  {:>9.6}  {:>10}  [{}]",
                nm.at_cents, nm.at_value, nm.depth, nm.sharpness, pq_label, conf
            );
        }
    }

    // Final summary of flags
    if flags.few_partials || flags.poor_f0_agreement || flags.flat_dissonance_curve {
        println!("QUALITY WARNINGS:");
        if flags.few_partials {
            println!(" - Few stable partials detected (< {}).", cfg.min_partials);
        }
        if flags.poor_f0_agreement {
            println!(
                " - f0 estimators disagree beyond {:.1} cents MAD.",
                cfg.max_f0_mad_cents
            );
        }
        if flags.flat_dissonance_curve {
            println!(
                " - Dissonance curve is flat (normalized minima depth < {:.3}).",
                cfg.min_norm_min_depth
            );
        }
    } else {
        println!("All quality checks passed.");
    }

    QualityResult {
        flags,
        n_partials_kept: n_partials,
        f0_mad_cents,
        dissonance_norm_depth: norm_depth,
    }
}
