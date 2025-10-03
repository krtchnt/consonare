//! Step 6: Perceptual weighting
//!
//! - Convert partial magnitudes (dB) to SPL proxies and apply A-weighting (dB).
//! - Within-note masking: down-weight partials > `mask_drop_db` below the
//!   strongest *A-weighted* component and attenuate components that fall within
//!   the ERB neighborhood of a stronger masker by at least `mask_margin_db`.
//! - Output per note: (f_i, w_i), where w_i is a unit-less weight in [0, 1].
//!
//! References:
//! - IEC 61672 A-weighting approximation
//! - Glasberg & Moore ERB: ERB(f) ≈ 24.7 * (4.37 * f_kHz + 1) [Hz]

use crate::partials::PartialsResult;

/// Configuration for Step 6 perceptual weighting.
#[derive(Clone, Debug)]
pub struct WeightingConfig {
    /// Drop partials more than this many dB below the strongest A-weighted partial.
    pub mask_drop_db: f32, // e.g., 40.0
    /// If a partial lies within ERB of a stronger neighbor by at least this margin,
    /// scale it by `erb_mask_atten` instead of dropping completely.
    pub mask_margin_db: f32, // e.g., 30.0
    /// Attenuation factor (linear, 0..1) when ERB-masked by a stronger local component.
    pub erb_mask_atten: f32, // e.g., 0.25
    /// Maximum number of partials to consider from `summaries` (sorted by freq).
    pub max_partials: usize, // e.g., 20
    /// If false, skip ERB masking (treat all partials as unmasked).
    pub enable_masking: bool,
    /// If false, set weights to 1.0 and A-weighting contribution to 0 dB.
    pub enable_perceptual_weighting: bool,
}

impl Default for WeightingConfig {
    fn default() -> Self {
        Self {
            mask_drop_db: 40.0,
            mask_margin_db: 30.0,
            erb_mask_atten: 0.25,
            max_partials: 20,
            enable_masking: true,
            enable_perceptual_weighting: true,
        }
    }
}

/// Output structure: weighted partials and some diagnostics.
#[derive(Clone, Debug)]
pub struct WeightedPartial {
    pub freq_hz: f32,
    /// Weight in [0,1], normalized to the strongest partial after A-weighting and masking.
    pub weight: f32,
    /// A-weighting applied (dB).
    pub a_weight_db: f32,
    /// Original median level (dB) from Step 4 summaries.
    pub median_db: f32,
    /// Whether ERB masking attenuation was applied.
    pub erb_masked: bool,
}

#[derive(Clone, Debug)]
pub struct WeightingResult {
    pub partials: Vec<WeightedPartial>,
    /// Index (into `partials`) of the strongest component used for normalization, if any.
    pub anchor_idx: Option<usize>,
}

/// Main entry point for Step 6.
pub fn run_weighting_step(parts: &PartialsResult, cfg: &WeightingConfig) -> WeightingResult {
    // 1) Collect top-N partials from summaries (already sorted by frequency).
    let candidates: Vec<PartialCandidate> = parts
        .summaries
        .iter()
        .take(cfg.max_partials)
        .map(|s| PartialCandidate {
            freq_hz: s.median_hz,
            median_db: s.median_db,
            a_db: a_weight_db(s.median_hz),
            level_aw_db: s.median_db + a_weight_db(s.median_hz),
            erb_hz: erb_bw_hz(s.median_hz),
        })
        .collect();

    if candidates.is_empty() {
        return WeightingResult {
            partials: vec![],
            anchor_idx: None,
        };
    }

    // 2) Find the strongest A-weighted component as the anchor.
    let (anchor_idx, max_level) = candidates
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.level_aw_db.partial_cmp(&b.1.level_aw_db).unwrap())
        .map(|(i, c)| (Some(i), c.level_aw_db))
        .unwrap_or((None, f32::NEG_INFINITY));

    // 3) Apply within-note masking rules.
    let mut out: Vec<WeightedPartial> = Vec::with_capacity(candidates.len());
    for i in 0..candidates.len() {
        let c = &candidates[i];

        // Global drop vs. anchor
        let rel_to_anchor = max_level - c.level_aw_db;
        if rel_to_anchor > cfg.mask_drop_db {
            // Fully drop
            out.push(WeightedPartial {
                freq_hz: c.freq_hz,
                weight: 0.0,
                a_weight_db: c.a_db,
                median_db: c.median_db,
                erb_masked: false,
            });
            continue;
        }

        // ERB masking by any stronger neighbor
        let mut masked = false;
        if cfg.enable_masking {
            for (j, candidates_j) in candidates.iter().enumerate() {
                if j == i {
                    continue;
                }
                let cj = &candidates_j;
                if cj.level_aw_db >= c.level_aw_db + cfg.mask_margin_db
                    && (c.freq_hz - cj.freq_hz).abs() <= cj.erb_hz
                {
                    masked = true;
                    break;
                }
            }
        }

        // Linear weight relative to anchor using dB difference after A-weighting.
        let mut w_lin = if cfg.enable_perceptual_weighting {
            // 0 dB -> 1.0, -6 dB -> ~0.5, etc.
            db_to_linear(-(rel_to_anchor)) // A-weighted, relative to anchor
        } else {
            1.0 // flat weighting
        };
        if masked {
            w_lin *= cfg.erb_mask_atten;
        }

        out.push(WeightedPartial {
            freq_hz: c.freq_hz,
            weight: w_lin,
            a_weight_db: if cfg.enable_perceptual_weighting {
                c.a_db
            } else {
                0.0
            },
            median_db: c.median_db,
            erb_masked: masked,
        });
    }

    // 4) Normalize to max(weight)=1 (safety: anchor should already be 1.0, but normalize anyway).
    if let Some(max_w) = out
        .iter()
        .map(|p| p.weight)
        .fold(None, |acc: Option<f32>, w| {
            Some(acc.map_or(w, |m| if w > m { w } else { m }))
        })
        && max_w > 0.0
    {
        for p in &mut out {
            p.weight /= max_w;
        }
    }

    WeightingResult {
        partials: out,
        anchor_idx,
    }
}

/// Internal representation to compute A-weighted levels and ERB bandwidths.
#[derive(Clone, Debug)]
struct PartialCandidate {
    freq_hz: f32,
    median_db: f32,
    a_db: f32,
    level_aw_db: f32,
    erb_hz: f32,
}

// IEC 61672 A-weighting (dB), f > 0 (Hz)
pub fn a_weight_db(f_hz: f32) -> f32 {
    let f = f_hz.max(1e-3);
    let f2 = f * f;
    let ra = (12200.0_f32.powi(2) * f2 * f2)
        / ((f2 + 20.6_f32.powi(2))
            * ((f2 + 107.7_f32.powi(2)) * (f2 + 737.9_f32.powi(2))).sqrt()
            * (f2 + 12200.0_f32.powi(2)));
    20.0 * ra.log10() + 2.0 // A(1 kHz) ≈ 0 dB
}

/// Glasberg & Moore ERB bandwidth in Hz:  ERB(f) ≈ 24.7 * (4.37 * f_kHz + 1).
pub fn erb_bw_hz(f_hz: f32) -> f32 {
    let fk = f_hz / 1000.0;
    24.7 * (4.37 * fk + 1.0)
}

#[inline]
fn db_to_linear(db: f32) -> f32 {
    10f32.powf(db / 20.0)
}
