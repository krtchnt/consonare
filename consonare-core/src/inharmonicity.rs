//! Step 5: Inharmonicity modelling (optional)
//!
//! Fit the stiff-string inharmonicity coefficient `B` from measured partials.
//!
//! Model (Fletcher & Rossing):
//!     f_n = n * f0 * sqrt(1 + B * n^2)
//!
//! Using an f0 estimate \hat{f0}, define r_n = f_n / (n * \hat{f0}).
//! Then r_n^2 = s^2 (1 + B n^2) with s = f0_true / \hat{f0}.
//! This is linear in n^2:  y = a + b x, with y=r_n^2, x=n^2, a=s^2, b=s^2 B.
//! Hence B = b/a and a provides a refined f0: f0_refined = sqrt(a) * \hat{f0}.
//!
//! We perform a (optionally weighted) least-squares fit over reliable partials
//! chosen from Step 4 summaries. If the fit quality is poor (few points,
//! negative/degenerate `a`, or R^2 below a threshold) we skip and return None.

use crate::{common::median, partials::PartialsResult};

#[derive(Clone, Debug)]
pub struct InharmonicityConfig {
    /// Minimum number of distinct partial indices (n) required to attempt a fit.
    pub min_points: usize,
    /// Maximum harmonic order to consider (helps reject wild outliers).
    pub max_n: usize,
    /// Cents tolerance when mapping a partial frequency to an integer n.
    pub map_tol_cents: f32,
    /// Minimum median SNR (dB) for a summary to be used.
    pub min_snr_db: f32,
    /// Require this minimum coefficient of determination to accept the model.
    pub min_r2_accept: f32,
    /// If true, weight the regression by exp(snr_db / 20).
    pub weight_by_snr: bool,
}

impl Default for InharmonicityConfig {
    fn default() -> Self {
        Self {
            min_points: 5,
            max_n: 32,
            map_tol_cents: 30.0,
            min_snr_db: 10.0,
            min_r2_accept: 0.85,
            weight_by_snr: true,
        }
    }
}

#[derive(Clone, Debug)]
pub struct InharmonicityResult {
    /// Estimated coefficient B if accepted, otherwise None.
    pub b: Option<f32>,
    /// Refined f0 (Hz) if fit accepted; otherwise the input f0_hat is echoed.
    pub f0_refined_hz: f32,
    /// Number of points used in the regression.
    pub used_points: usize,
    /// Coefficient of determination of the linear fit.
    pub r2: f32,
    /// Median absolute cents error of the model over used points.
    pub median_abs_cents: Option<f32>,
    /// The (n, measured_f_hz) tuples that were used.
    pub used_pairs: Vec<(usize, f32)>,
}

/// Map measured partial frequency to nearest integer harmonic index n given f0.
fn map_to_harmonic_index(f_hz: f32, f0_hz: f32) -> (usize, f32) {
    if f0_hz <= 0.0 || f_hz <= 0.0 {
        return (0, f_hz);
    }
    let n = (f_hz / f0_hz).round().max(1.0) as usize;
    (n, f_hz)
}

/// Weighted linear regression for y = a + b x.
fn weighted_linreg(x: &[f32], y: &[f32], w: &[f32]) -> Option<(f32, f32, f32)> {
    if x.len() != y.len() || x.len() != w.len() || x.is_empty() {
        return None;
    }
    let mut sw = 0.0f64;
    let mut swx = 0.0f64;
    let mut swy = 0.0f64;
    let mut swxx = 0.0f64;
    let mut swxy = 0.0f64;

    for i in 0..x.len() {
        let wi = w[i].max(1e-12) as f64;
        let xi = x[i] as f64;
        let yi = y[i] as f64;
        sw += wi;
        swx += wi * xi;
        swy += wi * yi;
        swxx += wi * xi * xi;
        swxy += wi * xi * yi;
    }
    let denom = sw * swxx - swx * swx;
    if denom.abs() < 1e-12 {
        return None;
    }
    let b = (sw * swxy - swx * swy) / denom;
    let a = (swy - b * swx) / sw;

    // Compute weighted R^2
    let mut ss_tot = 0.0f64;
    let y_bar = swy / sw;
    let mut ss_res = 0.0f64;
    for i in 0..x.len() {
        let wi = w[i].max(1e-12) as f64;
        let yi = y[i] as f64;
        let yhat = a + b * (x[i] as f64);
        ss_tot += wi * (yi - y_bar) * (yi - y_bar);
        ss_res += wi * (yi - yhat) * (yi - yhat);
    }
    let r2 = if ss_tot <= 0.0 {
        1.0
    } else {
        1.0 - (ss_res / ss_tot)
    };
    Some((a as f32, b as f32, r2 as f32))
}

/// Cents error between two positive frequencies.
fn cents_err(f_meas: f32, f_pred: f32) -> f32 {
    if f_meas <= 0.0 || f_pred <= 0.0 {
        return f32::INFINITY;
    }
    1200.0 * ((f_meas / f_pred) as f64).log2() as f32
}

/// Main entry point for Step 5.
pub fn run_inharmonicity_step(
    parts: &PartialsResult,
    f0_hat_hz: Option<f32>,
    cfg: &InharmonicityConfig,
) -> InharmonicityResult {
    let f0_hat = f0_hat_hz.unwrap_or(0.0);
    // Gather candidate (n, f) from partial summaries.
    let mut pairs: Vec<(usize, f32, f32)> = Vec::new(); // (n, f_hz, weight)
    for s in &parts.summaries {
        if s.median_snr_db < cfg.min_snr_db || s.median_hz <= 0.0 {
            continue;
        }
        if f0_hat > 0.0 {
            let (n, _f) = map_to_harmonic_index(s.median_hz, f0_hat);
            if n == 0 || n > cfg.max_n {
                continue;
            }
            // Check mapping tolerance in cents against the nearest harmonic of f0_hat.
            let pred = (n as f32) * f0_hat;
            let cents = cents_err(s.median_hz, pred).abs();
            if cents > cfg.map_tol_cents {
                continue;
            }
            let w = if cfg.weight_by_snr {
                (10f32).powf(s.median_snr_db / 20.0) // proportional to linear SNR
            } else {
                1.0
            };
            pairs.push((n, s.median_hz, w));
        }
    }
    // Keep only the strongest occurrence per n (if multiple tracks map to same n).
    pairs.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| b.2.total_cmp(&a.2)));
    pairs.dedup_by(|a, b| a.0 == b.0); // assume sorted by n then weight
    let used_points = pairs.len();

    if used_points < cfg.min_points || f0_hat <= 0.0 {
        return InharmonicityResult {
            b: None,
            f0_refined_hz: f0_hat.max(0.0),
            used_points,
            r2: 0.0,
            median_abs_cents: None,
            used_pairs: pairs.into_iter().map(|(n, f, _)| (n, f)).collect(),
        };
    }

    // Build regression data: y = r_n^2 = (f_n / (n f0_hat))^2 ; x = n^2
    let mut x: Vec<f32> = Vec::with_capacity(used_points);
    let mut y: Vec<f32> = Vec::with_capacity(used_points);
    let mut w: Vec<f32> = Vec::with_capacity(used_points);
    for (n, f, wi) in &pairs {
        let nn = (*n as f32) * (*n as f32);
        let rn = f / ((*n as f32) * f0_hat);
        x.push(nn);
        y.push(rn * rn);
        w.push(*wi);
    }

    let (a, b, r2) = match weighted_linreg(&x, &y, &w) {
        Some(t) => t,
        None => {
            return InharmonicityResult {
                b: None,
                f0_refined_hz: f0_hat,
                used_points,
                r2: 0.0,
                median_abs_cents: None,
                used_pairs: pairs.into_iter().map(|(n, f, _)| (n, f)).collect(),
            };
        }
    };

    // Recover s^2=a and B=b/a. Accept only if positive and fit quality good.
    if a <= 0.0 || !a.is_finite() || !b.is_finite() || r2 < cfg.min_r2_accept {
        return InharmonicityResult {
            b: None,
            f0_refined_hz: f0_hat,
            used_points,
            r2,
            median_abs_cents: None,
            used_pairs: pairs.into_iter().map(|(n, f, _)| (n, f)).collect(),
        };
    }
    let s = a.sqrt();
    let b_est = (b / a).max(0.0);

    // Compute a quick median cents deviation against the accepted model.
    let mut cents_errors: Vec<f32> = Vec::with_capacity(used_points);
    for (n, f, _) in &pairs {
        let pred = (*n as f32) * (s * f0_hat) * (1.0 + b_est * (*n as f32) * (*n as f32)).sqrt();
        cents_errors.push(cents_err(*f, pred).abs());
    }
    let med_cents = median(&mut cents_errors);

    InharmonicityResult {
        b: Some(b_est),
        f0_refined_hz: s * f0_hat,
        used_points,
        r2,
        median_abs_cents: Some(med_cents),
        used_pairs: pairs.into_iter().map(|(n, f, _)| (n, f)).collect(),
    }
}
