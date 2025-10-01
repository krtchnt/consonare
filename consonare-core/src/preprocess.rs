//! Spectral preprocessing (Step #2): band-limit, STFT, time-median, noise-floor suppression.
//!
//! Implements the spec:
//! - High-pass at 20 Hz; low-pass at min(Nyquist, 12 kHz).
//! - STFT with Hann window (4096–8192 samples recommended), ~75% overlap.
//! - Magnitude spectra; median-filter across time to emphasize stable partials.
//! - Noise floor estimate via percentile; suppress bins below `floor + k` dB.
//!
//! The API is intentionally standalone, but mirrors Step #1 style (Config/Result,
//! thiserror, etc.). Feed it mono f32 samples ([-1,1]) and the sample rate.

use realfft::{RealFftPlanner, num_complex::Complex32};
use thiserror::Error;

use crate::common::median;

#[derive(Debug, Error)]
pub enum SpectralError {
    #[error("invalid-arg: {0}")]
    InvalidArg(String),
}

#[derive(Clone, Debug)]
pub struct SpectralConfig {
    /// High-pass cutoff in Hz (applied as spectral mask).
    pub highpass_hz: f32, // e.g., 20.0
    /// Low-pass cutoff in Hz (applied as spectral mask); final cutoff = min(lowpass_hz, Nyquist).
    pub lowpass_hz: f32, // e.g., 12_000.0
    /// STFT window size in samples (e.g., 4096 or 8192).
    pub win_size: usize,
    /// Hop size in samples (win_size/4 ≈ 75% overlap if None); if 0, computed as win_size/4.
    pub hop_size: usize,
    /// Time-median filter length in frames (odd, e.g., 5 or 7).
    pub median_len: usize,
    /// Percentile (0–100) for frame-wise noise floor on dB spectrogram (e.g., 20.0).
    pub floor_percentile: f32,
    /// Suppress bins below floor + k dB (e.g., k in 6–12).
    pub suppress_k_db: f32,
    /// Use power (mag^2) before dB conversion if true; else magnitude.
    pub use_power: bool,
}

impl Default for SpectralConfig {
    fn default() -> Self {
        Self {
            highpass_hz: 20.0,
            lowpass_hz: 12_000.0,
            win_size: 4096,
            hop_size: 0, // interpret as win/4
            median_len: 5,
            floor_percentile: 20.0,
            suppress_k_db: 9.0,
            use_power: false,
        }
    }
}

#[derive(Clone, Debug)]
pub struct SpectralResult {
    /// Center time (s) for each frame.
    pub times_s: Vec<f32>,
    /// Frequency axis for bins (Hz), length = n_bins = win_size/2 + 1.
    pub freqs_hz: Vec<f32>,
    /// Unsuppressed dB spectrogram *after* time-median smoothing; [frames][bins].
    pub mag_db: Vec<Vec<f32>>,
    /// Per-frame noise floor (dB) from percentile.
    pub floor_db: Vec<f32>,
    /// Mask of bins suppressed by (mag_db < floor_db + k_db).
    pub suppressed: Vec<Vec<bool>>,
}

/// Run Step #2 on mono time-domain samples.
pub fn run_spectral_step(
    samples: &[f32],
    sample_rate: u32,
    cfg: &SpectralConfig,
) -> Result<SpectralResult, SpectralError> {
    if samples.is_empty() {
        return Err(SpectralError::InvalidArg("empty signal".into()));
    }
    if sample_rate == 0 {
        return Err(SpectralError::InvalidArg("sample_rate must be > 0".into()));
    }
    let sr = sample_rate as f32;

    let n = cfg.win_size.max(8);
    let hop = if cfg.hop_size == 0 {
        (cfg.win_size / 4).max(1)
    } else {
        cfg.hop_size
    };
    let n_bins = n / 2 + 1;
    if n_bins < 3 {
        return Err(SpectralError::InvalidArg("win_size too small".into()));
    }
    if cfg.median_len == 0 || cfg.median_len.is_multiple_of(2) {
        return Err(SpectralError::InvalidArg(
            "median_len must be positive and odd".into(),
        ));
    }
    if !(0.0..=100.0).contains(&cfg.floor_percentile) {
        return Err(SpectralError::InvalidArg(
            "floor_percentile must be within [0,100]".into(),
        ));
    }

    // Precompute Hann window and frequency axis.
    let hann: Vec<f32> = (0..n)
        .map(|i| 0.5 - 0.5 * ((2.0 * std::f32::consts::PI * i as f32) / (n as f32 - 1.0)).cos())
        .collect();
    let freqs_hz: Vec<f32> = (0..n_bins).map(|k| (k as f32) * sr / (n as f32)).collect();

    // Band mask derived from HP/LP spec (applied in frequency domain).
    let nyq = sr * 0.5;
    let hp = cfg.highpass_hz.max(0.0);
    let lp = cfg.lowpass_hz.min(nyq).max(hp);
    let band_mask: Vec<bool> = freqs_hz.iter().map(|&f| f >= hp && f <= lp).collect();

    // STFT using realfft.
    let mut planner = RealFftPlanner::<f32>::new();
    let r2c = planner.plan_fft_forward(n);
    let mut in_buf = r2c.make_input_vec();
    let mut spec_buf = r2c.make_output_vec();

    // Collect magnitude spectra per frame (pre-smoothing), linear magnitude.
    let mut mags: Vec<Vec<f32>> = Vec::new();
    let mut times: Vec<f32> = Vec::new();

    let mut i = 0usize;
    while i + n <= samples.len() {
        // Windowed frame
        for j in 0..n {
            in_buf[j] = samples[i + j] * hann[j];
        }
        r2c.process(&mut in_buf, &mut spec_buf)
            .map_err(|e| SpectralError::InvalidArg(format!("FFT error: {e}")))?;

        let mut mag = vec![0.0f32; n_bins];
        for k in 0..n_bins {
            let Complex32 { re, im } = spec_buf[k];
            let m = (re * re + im * im).sqrt();
            mag[k] = if band_mask[k] { m } else { 0.0 };
        }
        mags.push(mag);
        let center = (i + n / 2) as f32 / sr;
        times.push(center);
        i += hop;
    }

    // Time-median smoothing across frames per bin.
    let smoothed = time_median_per_bin(&mags, cfg.median_len);

    // Convert to dB (magnitude or power), estimate per-frame floor, and mask.
    let mut mag_db = vec![vec![0.0f32; n_bins]; smoothed.len()];
    let mut floor_db = vec![0.0f32; smoothed.len()];
    let mut suppressed = vec![vec![false; n_bins]; smoothed.len()];

    for (t, row) in smoothed.iter().enumerate() {
        // Convert to dB
        if cfg.use_power {
            // power = mag^2
            for (k, row_k) in row.iter().enumerate().take(n_bins) {
                let p = row_k * row_k;
                mag_db[t][k] = pow_to_db(p);
            }
        } else {
            for (k, &row_k) in row.iter().enumerate().take(n_bins) {
                mag_db[t][k] = amp_to_db(row_k);
            }
        }
        // Frame noise floor via percentile across bins (on dB domain).
        let floor = percentile(&mag_db[t], cfg.floor_percentile);
        floor_db[t] = floor;
        let th = floor + cfg.suppress_k_db;
        for (k, band_mask_k) in band_mask.iter().enumerate().take(n_bins) {
            suppressed[t][k] = mag_db[t][k] < th || !band_mask_k;
        }
    }

    Ok(SpectralResult {
        times_s: times,
        freqs_hz,
        mag_db,
        floor_db,
        suppressed,
    })
}

fn time_median_per_bin(mags: &[Vec<f32>], win_len: usize) -> Vec<Vec<f32>> {
    if mags.is_empty() {
        return vec![];
    }
    let t_len = mags.len();
    let n_bins = mags[0].len();
    let half = win_len / 2;
    let mut out = vec![vec![0.0f32; n_bins]; t_len];

    // Naive O(T * B * W log W) median; fine for offline analysis.
    for (t, out_t) in out.iter_mut().enumerate().take(t_len) {
        let t0 = t.saturating_sub(half);
        let t1 = (t + half + 1).min(t_len);
        let span = t1 - t0;
        for k in 0..n_bins {
            let mut buf = Vec::with_capacity(span);
            for mags_tt in mags.iter().take(t1).skip(t0) {
                buf.push(mags_tt[k]);
            }
            out_t[k] = median(&mut buf);
        }
    }
    out
}

#[inline]
fn amp_to_db(a: f32) -> f32 {
    // 20*log10(a); floor at -160 dB to be safe.
    if a <= 0.0 {
        return -160.0;
    }
    20.0 * a.log10().max(-8.0) // -8 decades ~ -160 dB
}

#[inline]
fn pow_to_db(p: f32) -> f32 {
    // 10*log10(p); floor at -160 dB
    if p <= 0.0 {
        return -160.0;
    }
    10.0 * p.log10().max(-16.0)
}

/// Percentile on a slice; p in [0,100]. Linear interpolation between ranks.
fn percentile(xs: &[f32], p: f32) -> f32 {
    assert!(!xs.is_empty());
    let mut v: Vec<f32> = xs.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if p <= 0.0 {
        return v[0];
    }
    if p >= 100.0 {
        return v[v.len() - 1];
    }
    let pos = (p / 100.0) * ((v.len() - 1) as f32);
    let i = pos.floor() as usize;
    let frac = pos - (i as f32);
    if i + 1 < v.len() {
        v[i] * (1.0 - frac) + v[i + 1] * frac
    } else {
        v[i]
    }
}

// --------------------------- Tests -----------------------------

#[cfg(test)]
mod tests {
    use crate::common::median;

    use super::*;

    #[test]
    fn median_basic() {
        let mut v = vec![3.0, 1.0, 4.0, 1.5, 2.0];
        assert!((median(&mut v) - 2.0).abs() < 1e-6);
    }

    #[test]
    fn percentile_basic() {
        let xs = [0.0, 1.0, 2.0, 3.0, 4.0];
        assert!((percentile(&xs, 0.0) - 0.0).abs() < 1e-6);
        assert!((percentile(&xs, 50.0) - 2.0).abs() < 1e-6);
        assert!((percentile(&xs, 100.0) - 4.0).abs() < 1e-6);
        assert!((percentile(&xs, 25.0) - 1.0).abs() < 1e-6);
        assert!((percentile(&xs, 75.0) - 3.0).abs() < 1e-6);
    }
}
