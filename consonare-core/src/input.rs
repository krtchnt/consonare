//! Steady-pitch preprocessing and checks (using `pitch_detection::YINDetector`).
//!
//! Pipeline:
//! 1) Decode mono WAV/FLAC (sr ≥ 44_100, bits ≥ 16).
//! 2) Trim leading/trailing silence (RMS gate).
//! 3) Normalize peak to −1 dBFS.
//! 4) Track f0(t) with `pitch_detection::detector::yin::YINDetector` over sliding windows;
//!    flag if max |1200*log2(f0/median)| > 100 cents.

use std::fs::File;
use std::path::Path;

use symphonia::core::{
    audio::{SampleBuffer, SignalSpec},
    codecs::DecoderOptions,
    formats::FormatOptions,
    io::MediaSourceStream,
    meta::MetadataOptions,
    probe::Hint,
};
use thiserror::Error;

// NEW: bring in YIN from the pitch_detection crate
use pitch_detection::detector::PitchDetector;
use pitch_detection::detector::yin::YINDetector;

#[derive(Debug, Error)]
pub enum PitchError {
    #[error("decode error: {0}")]
    Decode(String),
    #[error("unsupported: {0}")]
    Unsupported(String),
    #[error("invalid-data: {0}")]
    InvalidData(String),
}

#[derive(Clone, Debug)]
pub struct Config {
    /// Silence trim RMS threshold in dBFS (e.g., -50.0)
    pub trim_rms_dbfs: f32,
    /// Min continuous non-silence required to keep (ms) – guards against choppy tails.
    pub min_keep_ms: u32,
    /// Peak normalization target in dBFS (e.g., -1.0)
    pub target_peak_dbfs: f32,
    /// Analysis window in ms (e.g., 50)
    pub yin_window_ms: u32,
    /// Hop in ms (e.g., 10)
    pub hop_ms: u32,
    /// YIN power threshold (rejects low-energy frames). Typical ~5.0 per crate examples.
    pub yin_power_threshold: f32,
    /// YIN clarity threshold ∈ [0,1] (confidence). Typical 0.6–0.9.
    pub yin_clarity_threshold: f32,
    /// Min/Max pitch search band (Hz)
    pub fmin_hz: f32,
    pub fmax_hz: f32,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            trim_rms_dbfs: -70.0, // was -50.0
            min_keep_ms: 80,      // >= yin_window_ms
            target_peak_dbfs: -1.0,
            yin_window_ms: 50,
            hop_ms: 10,
            yin_power_threshold: 5.0,
            yin_clarity_threshold: 0.7,
            fmin_hz: 40.0, // give room for very low pedal tones
            fmax_hz: 1000.0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct AnalysisResult {
    pub sample_rate: u32,
    pub bits_per_sample: Option<u32>,
    pub trimmed_samples: usize,
    pub original_samples: usize,
    pub applied_gain_db: f32,
    pub f0_hz: Vec<Option<f32>>, // per frame
    pub times_s: Vec<f32>,       // center time per frame
    pub median_f0_hz: Option<f32>,
    pub max_abs_dev_cents: Option<f32>,
    pub steady_pitch: bool,
}

#[derive(Clone, Debug)]
pub struct Processed {
    /// Audio after trim + normalization, mono f32 in [-1, 1].
    pub samples: Vec<f32>,
    pub analysis: AnalysisResult,
}

/// Main entry: load, check, trim, normalize, analyze.
pub fn process_file<P: AsRef<Path>>(path: P, cfg: &Config) -> Result<Processed, PitchError> {
    let Decoded {
        mut samples,
        sr,
        bits_per_sample,
        channels,
    } = decode_to_mono_f32(path)?;

    if channels != 1 {
        return Err(PitchError::Unsupported("Only mono is accepted".into()));
    }
    if sr < 44_100 {
        return Err(PitchError::Unsupported(format!(
            "Sample rate must be ≥ 44.1 kHz (got {})",
            sr
        )));
    }
    if let Some(bits) = bits_per_sample
        && bits < 16
    {
        return Err(PitchError::Unsupported(format!(
            "Bit depth must be ≥ 16-bit (got {}-bit)",
            bits
        )));
    }

    let original_samples = samples.len();

    // 1) Trim silence.
    let (start, end) = trim_silence_bounds_rms(&samples, sr, cfg.trim_rms_dbfs, cfg.min_keep_ms);
    samples = samples[start..end].to_vec();

    // 2) Normalize to target peak.
    let (gain_db, gain_lin) = gain_to_target_peak(&samples, cfg.target_peak_dbfs);
    if gain_lin.is_finite() && gain_lin > 0.0 {
        for s in &mut samples {
            *s *= gain_lin;
        }
    }

    // 3) Pitch tracking (YIN via crate).
    let (f0_hz, times_s) = track_f0_yin_pitch_detection(&samples, sr, cfg);

    let median_f0_hz = median_nonempty(&f0_hz);
    let (max_abs_dev_cents, steady_pitch) = if let Some(med) = median_f0_hz {
        let max_cents = f0_hz.iter().flatten().fold(0.0_f32, |acc, &f| {
            let cents = 1200.0 * (f / med).log2();
            acc.max(cents.abs())
        });
        (Some(max_cents), max_cents <= 100.0)
    } else {
        (None, false)
    };

    Ok(Processed {
        samples: samples.clone(),
        analysis: AnalysisResult {
            sample_rate: sr,
            bits_per_sample,
            trimmed_samples: samples.len(),
            original_samples,
            applied_gain_db: gain_db,
            f0_hz,
            times_s,
            median_f0_hz,
            max_abs_dev_cents,
            steady_pitch,
        },
    })
}

// --------------------------- Decoding -----------------------------

struct Decoded {
    samples: Vec<f32>,
    sr: u32,
    bits_per_sample: Option<u32>,
    channels: u16,
}

fn decode_to_mono_f32<P: AsRef<Path>>(path: P) -> Result<Decoded, PitchError> {
    let file = File::open(&path).map_err(|e| PitchError::Decode(e.to_string()))?;
    let mss = MediaSourceStream::new(Box::new(file), Default::default());

    let mut hint = Hint::new();
    if let Some(ext) = path.as_ref().extension().and_then(|s| s.to_str()) {
        hint.with_extension(ext);
    }

    let probed = symphonia::default::get_probe()
        .format(
            &hint,
            mss,
            &FormatOptions::default(),
            &MetadataOptions::default(),
        )
        .map_err(|e| PitchError::Decode(e.to_string()))?;
    let mut format = probed.format;

    let track = format
        .default_track()
        .ok_or_else(|| PitchError::InvalidData("No default audio track".into()))?;

    let mut decoder = symphonia::default::get_codecs()
        .make(&track.codec_params, &DecoderOptions::default())
        .map_err(|e| PitchError::Decode(e.to_string()))?;

    let sr = track
        .codec_params
        .sample_rate
        .ok_or_else(|| PitchError::InvalidData("Missing sample rate".into()))?;

    let bits_per_sample = track
        .codec_params
        .bits_per_coded_sample
        .or(track.codec_params.bits_per_sample);
    let channels = track
        .codec_params
        .channels
        .map(|c| c.count() as u16)
        .unwrap_or(1);

    let mut out = Vec::<f32>::new();
    let mut _spec: Option<SignalSpec> = None;

    // Correct, single-pass decode of all packets → mono f32
    loop {
        let packet = match format.next_packet() {
            Ok(p) => p,
            Err(err) => {
                if matches!(err, symphonia::core::errors::Error::ResetRequired) {
                    return Err(PitchError::Decode(
                        "Decoder requires a reset (unsupported midstream change)".into(),
                    ));
                }
                break; // End of stream
            }
        };

        let decoded = decoder
            .decode(&packet)
            .map_err(|e| PitchError::Decode(e.to_string()))?;

        if _spec.is_none() {
            _spec = Some(*decoded.spec());
        }

        let spec = *decoded.spec();
        let ch = spec.channels.count();
        let mut sbuf = SampleBuffer::<f32>::new(decoded.capacity() as u64, spec);
        sbuf.copy_interleaved_ref(decoded);
        push_planar_as_mono(sbuf.samples(), ch, &mut out);
    }

    Ok(Decoded {
        samples: out,
        sr,
        bits_per_sample,
        channels,
    })
}

fn push_planar_as_mono(samples: &[f32], channels: usize, out: &mut Vec<f32>) {
    if channels == 1 {
        out.extend_from_slice(samples);
        return;
    }
    let frames = samples.len() / channels;
    for i in 0..frames {
        let mut acc = 0.0f32;
        for c in 0..channels {
            acc += samples[i * channels + c];
        }
        out.push(acc / channels as f32);
    }
}

// --------------------------- Trimming -----------------------------

fn db_to_lin(db: f32) -> f32 {
    10f32.powf(db / 20.0)
}

/// Return [start, end) bounds to keep after trimming silence using short-time RMS.
fn trim_silence_bounds_rms(
    samples: &[f32],
    sr: u32,
    thresh_dbfs: f32,
    min_keep_ms: u32,
) -> (usize, usize) {
    if samples.is_empty() {
        return (0, 0);
    }
    // Short-time window for RMS (25 ms).
    let win = ((0.025f32 * sr as f32).round() as usize).max(8);
    let half = win / 2;
    let thr = db_to_lin(thresh_dbfs);

    let rms = moving_rms(samples, win);
    let mut start = 0usize;
    while start < rms.len() && rms[start] < thr {
        start += 1;
    }
    let mut end = rms.len();
    while end > 0 && rms[end - 1] < thr {
        end -= 1;
    }
    // If nothing passed the gate, keep a centered minimum chunk instead of returning (0,0).
    if start >= end {
        let min_keep = ((min_keep_ms as f32 / 1000.0) * sr as f32).round() as usize;
        if samples.len() <= min_keep {
            return (0, samples.len());
        }
        let center = samples.len() / 2;
        let half_keep = min_keep / 2;
        let s = center.saturating_sub(half_keep);
        let e = (center + half_keep).min(samples.len());
        return (s, e);
    }
    // Expand by half-window to map RMS centers back to sample indices.
    let mut s = start.saturating_mul(1).saturating_sub(half);
    let mut e = (end.saturating_mul(1)).saturating_add(half);
    s = s.min(samples.len());
    e = e.min(samples.len());

    // Keep a tiny minimum non-silent chunk to avoid over-trimming.
    let min_keep = ((min_keep_ms as f32 / 1000.0) * sr as f32).round() as usize;
    if e - s < min_keep && samples.len() >= min_keep {
        let center = samples.len() / 2;
        let half_keep = min_keep / 2;
        s = center.saturating_sub(half_keep);
        e = (center + half_keep).min(samples.len());
    }
    (s, e)
}

fn moving_rms(x: &[f32], win: usize) -> Vec<f32> {
    if x.is_empty() {
        return vec![];
    }
    let w = win.max(1);
    let mut out = Vec::with_capacity(x.len());
    let mut acc = 0.0f64;
    let mut q = std::collections::VecDeque::<f32>::new();
    for &s in x {
        let v = (s as f64) * (s as f64);
        acc += v;
        q.push_back(s);
        if q.len() > w {
            let u = q.pop_front().unwrap();
            acc -= (u as f64) * (u as f64);
        }
        let denom = q.len().max(1) as f64;
        out.push((acc / denom).sqrt() as f32);
    }
    out
}

// --------------------------- Normalization -----------------------------

/// Returns (gain_db, gain_lin)
fn gain_to_target_peak(samples: &[f32], target_dbfs: f32) -> (f32, f32) {
    let peak = samples
        .iter()
        .fold(0.0f32, |m, &s| m.max(s.abs()))
        .max(1e-12);
    let target_lin = db_to_lin(target_dbfs);
    let gain_lin = (target_lin / peak).min(1000.0);
    let gain_db = 20.0 * gain_lin.log10();
    (gain_db, gain_lin)
}

// --------------------------- YIN f0 via crate -----------------------------

fn track_f0_yin_pitch_detection(
    x: &[f32],
    sr: u32,
    Config {
        yin_window_ms: win_ms,
        hop_ms,
        yin_power_threshold: power_threshold,
        yin_clarity_threshold: clarity_threshold,
        fmin_hz,
        fmax_hz,
        ..
    }: &Config,
) -> (Vec<Option<f32>>, Vec<f32>) {
    let sr_f = sr as f32;
    let w = ((*win_ms as f32 / 1000.0) * sr_f).round() as usize;
    let hop = ((*hop_ms as f32 / 1000.0) * sr_f).round() as usize;
    let w = w.max(32);
    let hop = hop.max(1);

    // The detector can be reused across frames for efficiency.
    let padding = w / 2; // matches crate examples (SIZE/2)
    let mut detector = YINDetector::<f32>::new(w, padding);

    let mut f0s = Vec::new();
    let mut times = Vec::new();

    let mut i = 0usize;
    while i + w <= x.len() {
        let frame = &x[i..i + w];
        let pitch_opt =
            detector.get_pitch(frame, sr as usize, *power_threshold, *clarity_threshold);
        let f0 = pitch_opt
            .map(|p| p.frequency)
            .filter(|&f| f.is_finite() && f >= *fmin_hz && f <= *fmax_hz);
        let center_time = (i + w / 2) as f32 / sr_f;
        f0s.push(f0);
        times.push(center_time);
        i += hop;
    }
    (f0s, times)
}

// --------------------------- Stats helpers -----------------------------

fn median_nonempty(xs: &[Option<f32>]) -> Option<f32> {
    let mut v: Vec<f32> = xs.iter().copied().flatten().collect();
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let m = v.len() / 2;
    if v.len() % 2 == 1 {
        Some(v[m])
    } else {
        Some(0.5 * (v[m - 1] + v[m]))
    }
}

// --------------------------- Tests (quick sanity) -----------------------------

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn db_to_lin_roundtrip() {
        let db = -1.0;
        let lin = db_to_lin(db);
        assert!((20.0 * lin.log10() - db).abs() < 1e-5);
    }
}
