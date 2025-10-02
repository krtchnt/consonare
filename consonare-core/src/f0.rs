//! Step #3: Fundamental (f0) estimation via fast YIN (FFT autocorr) + log-spectral autocorrelation,
//! then consensus & refinement with a harmonic comb score, plus confidence reporting.
//!
//! This version replaces naive DFTs with realfft/rustfft FFTs and computes YIN using an
//! FFT-based autocorrelation + O(N) cumulative sums (a common "fast YIN" approach).
//!
//! Expected integration remains the same:
//! ```
//! use crate::preprocess::SpectralResult;
//! let f0 = fundamental::run_f0_step(&result.samples, sr, Some(&spec), &F0Config::default())?;
//! ```

use std::f32::consts::PI;

use crate::preprocess::SpectralResult;

use pitch_detection::detector::PitchDetector;
use pitch_detection::detector::yin::YINDetector;
use realfft::RealFftPlanner;

#[derive(Clone, Debug)]
pub struct F0Config {
    /// Minimum and maximum f0 to consider (Hz)
    pub fmin_hz: f32,
    pub fmax_hz: f32,
    /// Frame length for YIN in samples (should roughly match the STFT window)
    pub yin_frame: usize,
    /// Hop between frames (samples). If `None`, it will try to infer from SpectralResult.times_s
    /// and fall back to yin_frame/4.
    pub yin_hop: Option<usize>,
    /// YIN power threshold (reject low-energy frames). Typical ~5.0.
    pub yin_power_threshold: f32,
    /// YIN clarity threshold ∈ [0,1] (confidence). Typical 0.6–0.9.
    pub yin_clarity_threshold: f32,
    /// Number of harmonics to score in the spectral comb
    pub harm_peaks: usize,
    /// Bin tolerance in cents when checking harmonic peaks
    pub cents_tolerance: f32,
}

impl Default for F0Config {
    fn default() -> Self {
        Self {
            fmin_hz: 50.0,
            fmax_hz: 2000.0,
            yin_frame: 4096,
            yin_hop: None,
            yin_power_threshold: 5.0,
            yin_clarity_threshold: 0.7,
            harm_peaks: 6,
            cents_tolerance: 35.0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct F0Frame {
    pub time_s: f32,
    pub yin_hz: Option<f32>,
    pub yin_strength: f32, // 0..1 (1 is best)
    pub cep_hz: Option<f32>,
    pub cep_strength: f32,    // 0..1
    pub spec_hz: Option<f32>, // harmonic comb refined candidate from either yin/cep
    pub spec_score: f32,      // 0..1
    pub fused_hz: Option<f32>,
    /// Agreement proxy: higher when methods agree and are stable
    pub agreement: f32, // 0..1
}

#[derive(Clone, Debug)]
pub struct F0Result {
    pub frames: Vec<F0Frame>,
    /// Robust consensus over frames (RANSAC-ish median within an inlier band)
    pub consensus_f0_hz: Option<f32>,
    /// Global confidence in 0..1 based on agreement and inlier fraction
    pub consensus_confidence: f32,
}

pub fn run_f0_step(
    samples: &[f32],
    sr: u32,
    spec_opt: Option<&SpectralResult>,
    cfg: &F0Config,
) -> Result<F0Result, String> {
    if samples.is_empty() {
        return Ok(F0Result {
            frames: vec![],
            consensus_f0_hz: None,
            consensus_confidence: 0.0,
        });
    }

    let sr_f = sr as f32;
    let hop = cfg.yin_hop.unwrap_or(cfg.yin_frame / 4);

    // If we have STFT times, align YIN frame centers to those for better cross-method agreement.
    let frame_centers: Vec<usize> = if let Some(spec) = spec_opt {
        spec.times_s
            .iter()
            .map(|&t| {
                let c = (t * sr_f).round() as isize;
                c.clamp(0, (samples.len() as isize) - 1) as usize
            })
            .collect()
    } else {
        // Build our own grid
        let mut centers = Vec::new();
        let mut start = cfg.yin_frame as isize / 2;
        while (start as usize) + cfg.yin_frame / 2 < samples.len() {
            centers.push(start as usize);
            start += hop as isize;
        }
        centers
    };

    // Precompute spectral meta if provided
    let (spec_freqs_hz, spec_mag_db_row, spec_supp_row, _df) = if let Some(spec) = spec_opt {
        let df = if spec.freqs_hz.len() > 1 {
            (spec.freqs_hz[1] - spec.freqs_hz[0]).abs().max(1e-6)
        } else {
            // Fallback, approximate
            (sr_f / cfg.yin_frame as f32).max(1e-6)
        };
        (
            Some(&spec.freqs_hz),
            Some(&spec.mag_db),
            Some(&spec.suppressed),
            df,
        )
    } else {
        (None, None, None, (sr_f / cfg.yin_frame as f32).max(1e-6))
    };

    let mut frames = Vec::with_capacity(frame_centers.len());
    // A planner caches FFTs by size; reusing it avoids re-alloc and twiddle recompute.
    let mut fft_planner = RealFftPlanner::<f32>::new();
    // Reuse one YIN detector with fixed window size.
    let mut yin_detector = YINDetector::<f32>::new(cfg.yin_frame, cfg.yin_frame / 2);

    for (fi, &center) in frame_centers.iter().enumerate() {
        // Extract a Hann-windowed frame centered at `center`
        let half = cfg.yin_frame / 2;
        let start = center.saturating_sub(half);
        let end = (start + cfg.yin_frame).min(samples.len());
        if end <= start + 8 {
            continue;
        }
        // Build a fixed-size, Hann-windowed frame for YINDetector (zero-padded beyond bounds).
        let mut frame = vec![0.0_f32; cfg.yin_frame];
        for (j, frame_j) in frame.iter_mut().enumerate().take(cfg.yin_frame) {
            let idx = start + j;
            let s = if idx < samples.len() {
                samples[idx]
            } else {
                0.0
            };
            *frame_j = s * hann(cfg.yin_frame, j);
        }

        // --- YIN via `pitch_detection` crate ---
        let pitch = yin_detector.get_pitch(
            &frame,
            sr as usize,
            cfg.yin_power_threshold,
            cfg.yin_clarity_threshold,
        );
        let mut yin_hz = pitch.as_ref().map(|p| p.frequency);
        if let Some(f) = yin_hz
            && !(f.is_finite() && f >= cfg.fmin_hz && f <= cfg.fmax_hz)
        {
            yin_hz = None;
        }
        let yin_strength = pitch.map(|p| p.clarity).unwrap_or(0.0).clamp(0.0, 1.0);

        // --- Cepstral-like from spectrum ---
        let (cep_hz, cep_strength) = if let (Some(freqs), Some(mag_db_rows)) =
            (spec_freqs_hz, spec_mag_db_row)
        {
            if let Some(row) = mag_db_rows.get(fi) {
                let sup_row = spec_supp_row.and_then(|s| s.get(fi));
                spectral_autocorr_f0(row, freqs, sup_row.map(|v| &**v), cfg.fmin_hz, cfg.fmax_hz)
            } else {
                (None, 0.0)
            }
        } else {
            // If we don't have a spectrum row aligned, compute one directly from the frame via FFT
            spectral_autocorr_f0_from_frame(
                &frame,
                sr_f,
                cfg.fmin_hz,
                cfg.fmax_hz,
                &mut fft_planner,
            )
        };

        // --- Build candidate list ---
        let mut cand = vec![];
        if let Some(hz) = yin_hz {
            cand.push((hz, yin_strength, "yin"));
        }
        if let Some(hz) = cep_hz {
            cand.push((hz, cep_strength, "cep"));
        }

        // Evaluate each candidate with a harmonic comb over the magnitude spectrum (if available)
        let (spec_hz, spec_score) = if let (Some(freqs), Some(mag_db_rows)) =
            (spec_freqs_hz, spec_mag_db_row)
        {
            if let Some(row) = mag_db_rows.get(fi) {
                let best = cand
                    .iter()
                    .filter_map(|(hz, s, _)| {
                        let (refined, score) =
                            harmonic_refine(*hz, row, freqs, cfg.harm_peaks, cfg.cents_tolerance);
                        refined.map(|r| (r, score * *s))
                    })
                    .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

                match best {
                    Some((hz, score)) => (Some(hz), score),
                    None => (None, 0.0),
                }
            } else {
                (None, 0.0)
            }
        } else {
            (None, 0.0)
        };

        // Fuse: prefer agreement; otherwise best-by-score
        let (fused_hz, agree) = fuse_candidates(
            yin_hz,
            yin_strength,
            cep_hz,
            cep_strength,
            spec_hz,
            spec_score,
        );

        frames.push(F0Frame {
            time_s: if let Some(spec) = spec_opt {
                spec.times_s.get(fi).copied().unwrap_or(0.0)
            } else {
                center as f32 / sr_f
            },
            yin_hz,
            yin_strength,
            cep_hz,
            cep_strength,
            spec_hz,
            spec_score,
            fused_hz,
            agreement: agree,
        });
    }

    // RANSAC-ish consensus: pick the densest band within ±50 cents around the median,
    // then re-median inliers for final f0.
    let (consensus_hz, conf) = consensus_over_time(&frames);

    Ok(F0Result {
        frames,
        consensus_f0_hz: consensus_hz,
        consensus_confidence: conf,
    })
}

// ------------------------------ Helpers ------------------------------

fn hann(n: usize, i: usize) -> f32 {
    if n <= 1 {
        return 1.0;
    }
    (PI * 2.0 * i as f32 / (n as f32 - 1.0))
        .cos()
        .mul_add(-0.5, 0.5)
}

fn parabolic_interpolate(y_m1: f32, y0: f32, y_p1: f32, x0: f32) -> f32 {
    // Estimate vertex location of parabola passing through (-1,y_m1),(0,y0),(1,y_p1) shifted by x0
    let denom = y_m1 - 2.0 * y0 + y_p1;
    if denom.abs() < 1e-12 {
        return x0;
    }
    let delta = 0.5 * (y_m1 - y_p1) / denom;
    x0 + delta
}

// --------- Cepstral-like: autocorr over (log) spectrum to find harmonic spacing ----------

fn spectral_autocorr_f0(
    mag_db_row: &[f32],
    freqs_hz: &[f32],
    suppressed_row: Option<&[bool]>,
    fmin: f32,
    fmax: f32,
) -> (Option<f32>, f32) {
    if freqs_hz.len() < 3 || mag_db_row.len() != freqs_hz.len() {
        return (None, 0.0);
    }
    let df = (freqs_hz[1] - freqs_hz[0]).abs().max(1e-6);
    let mut lspec = Vec::with_capacity(mag_db_row.len());
    for (k, &db) in mag_db_row.iter().enumerate() {
        if let Some(sup) = suppressed_row
            && sup.get(k).copied().unwrap_or(false)
        {
            // De-emphasize suppressed bins
            lspec.push(db * 0.2);
            continue;
        }
        lspec.push(db);
    }
    // Mean-normalize to remove bias
    let mean = lspec.iter().copied().sum::<f32>() / (lspec.len().max(1) as f32);
    for v in &mut lspec {
        *v -= mean;
    }

    let lag_min = (fmin / df).ceil() as usize;
    let lag_max = (fmax / df).floor() as usize;
    if lag_max <= lag_min + 1 {
        return (None, 0.0);
    }

    // Raw autocorr; normalize by energy
    let energy = lspec.iter().map(|v| v * v).sum::<f32>().max(1e-12);
    let mut best_lag = 0usize;
    let mut best_val = f32::MIN;
    for lag in lag_min..=lag_max {
        let mut s = 0.0;
        let lim = lspec.len().saturating_sub(lag);
        for k in 0..lim {
            s += lspec[k] * lspec[k + lag];
        }
        let norm = s / energy;
        if norm > best_val {
            best_val = norm;
            best_lag = lag;
        }
    }

    if best_lag == 0 {
        return (None, 0.0);
    }

    // Parabolic interp on the autocorr curve
    let left = if best_lag > lag_min {
        autocorr_at_lag(&lspec, best_lag - 1) / energy
    } else {
        best_val
    };
    let right = if best_lag < lag_max {
        autocorr_at_lag(&lspec, best_lag + 1) / energy
    } else {
        best_val
    };
    let lag_f = parabolic_interpolate(left, best_val, right, best_lag as f32);
    let hz = lag_f * df;
    let strength = best_val.clamp(0.0, 1.0);
    (Some(hz), strength)
}

// Compute spectrum from frame via FFT (fast) and then reuse spectral_autocorr_f0
fn spectral_autocorr_f0_from_frame(
    frame: &[f32],
    sr: f32,
    fmin: f32,
    fmax: f32,
    planner: &mut RealFftPlanner<f32>,
) -> (Option<f32>, f32) {
    let n = frame.len().max(16);
    let n_fft = n.next_power_of_two().max(1024);
    let r2c = planner.plan_fft_forward(n_fft);

    // zero-padded input
    let mut x = r2c.make_input_vec();
    x[..frame.len()].copy_from_slice(frame);
    for x_i in x.iter_mut().take(n_fft).skip(frame.len()) {
        *x_i = 0.0;
    }

    // FFT
    let mut spec = r2c.make_output_vec();
    let _ = r2c.process(&mut x, &mut spec);

    // Magnitude -> dB
    let mut mag = vec![0.0_f32; spec.len()];
    for (k, s) in spec.iter().enumerate() {
        let m = s.norm(); // sqrt(re^2 + im^2)
        mag[k] = 20.0 * (m.max(1e-12)).log10();
    }

    // Frequency axis for real FFT half-spectrum
    let mut freqs = vec![0.0_f32; mag.len()];
    let df = sr / (n_fft as f32);
    for (k, f) in freqs.iter_mut().enumerate() {
        *f = k as f32 * df;
    }

    spectral_autocorr_f0(&mag, &freqs, None, fmin, fmax)
}

fn autocorr_at_lag(xs: &[f32], lag: usize) -> f32 {
    let lim = xs.len().saturating_sub(lag);
    let mut s = 0.0;
    for i in 0..lim {
        s += xs[i] * xs[i + lag];
    }
    s
}

fn harmonic_refine(
    f0_hz: f32,
    mag_db_row: &[f32],
    freqs_hz: &[f32],
    max_harm: usize,
    cents_tol: f32,
) -> (Option<f32>, f32) {
    if freqs_hz.len() < 3 || mag_db_row.len() != freqs_hz.len() {
        return (None, 0.0);
    }
    if f0_hz <= 0.0 {
        return (None, 0.0);
    }

    // Evaluate comb score near the candidate by sampling three offsets (for parabolic refine)
    let offs = [-0.015, 0.0, 0.015]; // ~±1.5% (~±26 cents) coarse search
    let mut vals = [0.0_f32; 3];
    for (j, off) in offs.iter().enumerate() {
        let f = f0_hz * (1.0 + *off);
        vals[j] = comb_score(f, mag_db_row, freqs_hz, max_harm, cents_tol);
    }
    // Parabolic refine in (offset, score) space
    let v_m1 = vals[0];
    let v0 = vals[1];
    let v_p1 = vals[2];
    let delta = parabolic_peak_location(v_m1, v0, v_p1).clamp(-0.03, 0.03);
    let refined = f0_hz * (1.0 + delta);
    let score = comb_score(refined, mag_db_row, freqs_hz, max_harm, cents_tol);
    (Some(refined), score)
}

fn parabolic_peak_location(y_m1: f32, y0: f32, y_p1: f32) -> f32 {
    let denom = y_m1 - 2.0 * y0 + y_p1;
    if denom.abs() < 1e-12 {
        return 0.0;
    }
    0.5 * (y_m1 - y_p1) / denom
}

fn comb_score(
    f0: f32,
    mag_db_row: &[f32],
    freqs_hz: &[f32],
    max_harm: usize,
    cents_tol: f32,
) -> f32 {
    if freqs_hz.len() < 2 {
        return 0.0;
    }
    let mut acc = 0.0_f32;
    let mut count = 0usize;
    let df = (freqs_hz[1] - freqs_hz[0]).abs().max(1e-6);
    let tol_ratio = 2f32.powf(cents_tol / 1200.0);
    let fmax = *freqs_hz.last().unwrap_or(&0.0);
    for h in 1..=max_harm {
        let target = f0 * (h as f32);
        if target > fmax {
            break;
        }
        // Find nearest bin within cents tolerance
        let k = (target / df).round() as isize;
        let mut best = None::<f32>;
        for dk in -2..=2 {
            let idx = k + dk;
            if idx < 0 || (idx as usize) >= freqs_hz.len() {
                continue;
            }
            let f = freqs_hz[idx as usize];
            let ratio = (f / target).max(1e-6);
            if ratio < 1.0 / tol_ratio || ratio > tol_ratio {
                continue;
            }
            let val = mag_db_row[idx as usize];
            if best.is_none_or(|b| val > b) {
                best = Some(val);
            }
        }
        if let Some(v) = best {
            acc += v;
            count += 1;
        }
    }
    if count == 0 {
        return 0.0;
    }
    // Normalize by count and dynamic range
    let avg = acc / count as f32;
    // Map dB-ish average to 0..1 via logistic-ish squash around 0 dB
    1.0 / (1.0 + (-avg / 12.0).exp())
}

fn fuse_candidates(
    yin_hz: Option<f32>,
    yin_s: f32,
    cep_hz: Option<f32>,
    cep_s: f32,
    spec_hz: Option<f32>,
    spec_score: f32,
) -> (Option<f32>, f32) {
    let mut agree = 0.0f32;
    let mut fused: Option<f32> = None;

    // Quick hum test
    let is_hum = |x: f32| (48.0..=52.5).contains(&x) || (59.0..=61.5).contains(&x);

    match (yin_hz, cep_hz) {
        (Some(a), Some(b)) => {
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };

            // Prefer the higher candidate when they form a near-integer 2–4× relation
            // (typical subharmonic confusion; bells/trombone often trigger this)
            let near_integer = (2..=4).any(|k| {
                let cents = 1200.0 * ((hi / (k as f32 * lo)).abs().max(1e-9)).log2().abs();
                cents < 30.0
            });

            if near_integer {
                fused = Some(hi);
                agree = (yin_s.max(cep_s) * 0.8).min(1.0);
            } else {
                // Original weighting when not an integer relation
                let cents = 1200.0 * ((a / b).max(1e-9)).log2().abs();
                let w = ((1.0 - (cents / 100.0)).clamp(0.0, 1.0)) * 0.7 + 0.3 * (yin_s.min(cep_s));
                agree = w;
                let wa = yin_s.max(1e-3);
                let wb = cep_s.max(1e-3);
                fused = Some((a * wa + b * wb) / (wa + wb));
            }

            // Hum guard: if one candidate sits at 50/60 Hz and the other is well above 120 Hz,
            // keep the higher one.
            if fused.is_some() && (is_hum(a) ^ is_hum(b)) && (a.max(b) > 120.0) {
                fused = Some(a.max(b));
                agree = agree.max(0.6);
            }
        }
        (Some(a), None) => {
            fused = Some(a);
            agree = yin_s * 0.6;
        }
        (None, Some(b)) => {
            fused = Some(b);
            agree = cep_s * 0.6;
        }
        (None, None) => { /* keep None */ }
    }

    // If spectral comb gives a strong refinement, prefer it
    if let (Some(f), Some(s)) = (spec_hz, fused) {
        let cents = 1200.0 * ((f / s).max(1e-9)).log2().abs();
        if spec_score > 0.5 && cents > 15.0 {
            fused = Some(f);
            agree = (agree * 0.5 + spec_score * 0.5).clamp(0.0, 1.0);
        } else {
            // small tweak toward spectral
            fused = Some((f * spec_score + s * (1.0 - spec_score)).max(1e-6));
            agree = (agree * 0.7 + spec_score * 0.3).clamp(0.0, 1.0);
        }
    } else if spec_hz.is_some() && (yin_hz.is_none() && cep_hz.is_none()) {
        fused = spec_hz;
        agree = spec_score * 0.6;
    }

    (fused, agree.clamp(0.0, 1.0))
}

fn consensus_over_time(frames: &[F0Frame]) -> (Option<f32>, f32) {
    let mut xs: Vec<f32> = frames.iter().filter_map(|f| f.fused_hz).collect();
    if xs.is_empty() {
        return (None, 0.0);
    }
    xs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = xs[xs.len() / 2];

    // Inliers within ±50 cents
    let tol = 50.0;
    let mut inliers = Vec::new();
    for &x in &xs {
        let cents = 1200.0 * ((x / med).max(1e-9)).log2().abs();
        if cents <= tol {
            inliers.push(x);
        }
    }
    if inliers.is_empty() {
        return (Some(med), 0.2);
    }
    inliers.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let consensus = inliers[inliers.len() / 2];
    let frac = inliers.len() as f32 / xs.len() as f32;

    // Confidence combines inlier fraction and average per-frame agreement
    let avg_agree = frames.iter().map(|f| f.agreement).sum::<f32>() / (frames.len() as f32);
    let conf = (0.6 * frac + 0.4 * avg_agree).clamp(0.0, 1.0);
    (Some(consensus), conf)
}
