//! Step 4: Partial extraction & tracking
//!
//! - Peak-pick per frame with SNR gate against the estimated noise floor
//! - Refine each peak via quadratic interpolation on log-magnitude (dB)
//!   around the bin maximum (parabolic three-point fit).
//! - Greedy track continuation across frames with a cents jump limit
//! - Report time-medians (Hz, dB, SNR) per track (robust to outliers)
//!
//! References:
//! - Parabolic (quadratic) peak interpolation: J.O. Smith, "Quadratic Interpolation of Spectral Peaks"
//! - Sinusoidal modeling & peak continuation: McAulay–Quatieri (1986); Serra SMS (1990)

use std::error::Error;

use crate::preprocess::SpectralResult;

/// Configuration for partial extraction & tracking.
#[derive(Debug, Clone)]
pub struct PartialsConfig {
    /// Max number of partials to keep (per frame and overall).
    pub max_partials: usize,
    /// Minimum per-frame SNR (dB) for a peak to be considered.
    pub min_snr_db: f32,
    /// Neighborhood half-width (in bins) for local-maximum test.
    pub local_max_half_width: usize,
    /// Maximum allowed inter-frame jump (in cents) when linking peaks.
    pub track_max_jump_cents: f32,
    /// Minimum number of observations required for a track to be reported.
    pub min_track_len: usize,
}

impl Default for PartialsConfig {
    fn default() -> Self {
        Self {
            max_partials: 20,
            min_snr_db: 15.0,           // "until SNR < 15 dB" (spec)
            local_max_half_width: 1,    // compare immediate neighbors by default
            track_max_jump_cents: 35.0, // conservative continuity
            min_track_len: 3,
        }
    }
}

/// One refined peak in a frame.
#[derive(Debug, Clone)]
pub struct Peak {
    pub bin: usize,       // integer FFT bin index of the local maximum
    pub offset_bins: f32, // sub-bin offset from quadratic fit ([-0.5, 0.5] typically)
    pub freq_hz: f32,     // refined frequency (Hz)
    pub mag_db: f32,      // refined log-magnitude (dB)
    pub snr_db: f32,      // mag_db - noise_floor_db(frame)
}

#[derive(Debug, Clone)]
pub struct FramePartials {
    pub frame_idx: usize,
    pub peaks: Vec<Peak>, // sorted descending by snr_db
}

/// One partial track across time.
#[derive(Debug, Clone)]
pub struct PartialTrack {
    pub id: usize,
    /// (frame_idx, freq_hz, mag_db, snr_db)
    pub points: Vec<(usize, f32, f32, f32)>,
}

#[derive(Debug, Clone)]
pub struct PartialSummary {
    pub track_id: usize,
    pub median_hz: f32,
    pub median_db: f32,
    pub median_snr_db: f32,
    pub count: usize,
}

#[derive(Debug, Clone)]
pub struct PartialsResult {
    pub frames: Vec<FramePartials>,
    pub tracks: Vec<PartialTrack>,
    pub summaries: Vec<PartialSummary>, // sorted by median_hz
}

/// Main entry point for Step 4.
pub fn run_partials_step(
    spec: &SpectralResult,
    _consensus_f0_hz: Option<f32>, // not required, but can be handy for downstream labeling
    cfg: &PartialsConfig,
) -> Result<PartialsResult, Box<dyn Error>> {
    // 1) Per-frame peak picking and refinement
    let frames = extract_refined_peaks(spec, cfg);

    // 2) Track continuation across frames (greedy nearest-neighbor in cents space)
    let tracks = link_tracks_greedy(&frames, cfg.track_max_jump_cents, cfg.max_partials);

    // 3) Compute time-medians for each track (robust to outliers)
    let mut summaries: Vec<PartialSummary> = tracks
        .iter()
        .filter(|t| t.points.len() >= cfg.min_track_len)
        .map(|t| PartialSummary {
            track_id: t.id,
            median_hz: median_f32(t.points.iter().map(|(_, f, _, _)| *f)),
            median_db: median_f32(t.points.iter().map(|(_, _, d, _)| *d)),
            median_snr_db: median_f32(t.points.iter().map(|(_, _, _, s)| *s)),
            count: t.points.len(),
        })
        .collect();

    // Stable ordering by frequency for readability.
    summaries.sort_by(|a, b| a.median_hz.partial_cmp(&b.median_hz).unwrap());

    Ok(PartialsResult {
        frames,
        tracks,
        summaries,
    })
}

/// Peak picking + refinement for all frames.
fn extract_refined_peaks(spec: &SpectralResult, cfg: &PartialsConfig) -> Vec<FramePartials> {
    let n_frames = spec.mag_db.len();
    let n_bins = spec.freqs_hz.len();
    let global_floor_med = if !spec.floor_db.is_empty() {
        median_f32(spec.floor_db.iter().copied())
    } else {
        -120.0
    };

    let mut out = Vec::with_capacity(n_frames);

    for (t, mags_db) in spec.mag_db.iter().enumerate() {
        // Per-frame noise floor (dB)
        let floor_db = spec.floor_db.get(t).copied().unwrap_or(global_floor_med);

        let suppressed = spec
            .suppressed
            .get(t)
            .map(|row| row.as_slice())
            .unwrap_or(&[]);

        // 1) Find local maxima with SNR >= threshold and not suppressed.
        let mut candidates: Vec<(usize, f32)> = Vec::new(); // (bin, snr_db)

        // avoid very first/last bins for 3-point interpolation
        for k in 1..n_bins.saturating_sub(1) {
            if !suppressed.is_empty() && suppressed[k] {
                continue;
            }

            // local-maximum test within neighborhood
            if is_local_max(mags_db, k, cfg.local_max_half_width) {
                let snr_db = mags_db[k] - floor_db;
                if snr_db >= cfg.min_snr_db {
                    candidates.push((k, snr_db));
                }
            }
        }

        // Sort by SNR desc, keep top-N
        candidates.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        if candidates.len() > cfg.max_partials {
            candidates.truncate(cfg.max_partials);
        }

        // 2) Refine each candidate via quadratic interpolation (in dB domain).
        let mut peaks: Vec<Peak> = Vec::with_capacity(candidates.len());
        for (k, _snr_db) in candidates {
            let (ym1, y0, yp1) = (mags_db[k - 1], mags_db[k], mags_db[k + 1]);
            let (offset, y_refined) = quad_interp_three_point(ym1, y0, yp1);

            // Convert sub-bin offset to Hz. Assume (near-)uniform bin spacing.
            let f0 = spec.freqs_hz[k];
            let bw = if k + 1 < n_bins {
                spec.freqs_hz[k + 1] - spec.freqs_hz[k]
            } else {
                spec.freqs_hz[k] - spec.freqs_hz[k - 1]
            };
            let f_ref = f0 + offset * bw;

            peaks.push(Peak {
                bin: k,
                offset_bins: offset,
                freq_hz: f_ref.max(0.0),
                mag_db: y_refined,
                snr_db: y_refined - floor_db,
            });
        }

        // Sort refined peaks by SNR desc for deterministic ordering.
        peaks.sort_by(|a, b| b.snr_db.partial_cmp(&a.snr_db).unwrap());

        out.push(FramePartials {
            frame_idx: t,
            peaks,
        });
    }

    out
}

/// Simple local-maximum test within a half-width neighborhood.
fn is_local_max(xs: &[f32], k: usize, half_w: usize) -> bool {
    let n = xs.len();
    let start = k.saturating_sub(half_w);
    let end = (k + half_w + 1).min(n);
    let center = xs[k];

    // strict greater-or-equal at center vs neighbors, and strictly greater than at least one side
    let mut ge_all = true;
    let mut gt_any = false;
    for (i, &v) in xs[start..end].iter().enumerate() {
        let idx = start + i;
        if idx == k {
            continue;
        }
        if center < v {
            ge_all = false;
            break;
        }
        if center > v {
            gt_any = true;
        }
    }
    ge_all && gt_any
}

/// Three-point quadratic interpolation on log-magnitude (dB).
///
/// Returns (offset_in_bins, refined_value_db).
/// Offset formula follows standard parabola fit through y[-1], y[0], y[+1]:
///   d = 0.5*(y[-1] - y[+1]) / (y[-1] - 2*y[0] + y[+1])
/// See J.O. Smith, "Quadratic Interpolation of Spectral Peaks".
fn quad_interp_three_point(ym1: f32, y0: f32, yp1: f32) -> (f32, f32) {
    let denom = ym1 - 2.0 * y0 + yp1;
    if denom.abs() < 1e-12 {
        return (0.0, y0);
    }
    let d = 0.5 * (ym1 - yp1) / denom; // sub-bin offset
    // Peak value of the parabola at x=d:
    // y(d) = y0 - 0.25*(ym1 - yp1)*d
    let y = y0 - 0.25 * (ym1 - yp1) * d;
    (d.clamp(-0.5, 0.5), y)
}

/// Greedy nearest-neighbor continuation with a cents jump cap.
fn link_tracks_greedy(
    frames: &[FramePartials],
    max_jump_cents: f32,
    max_tracks: usize,
) -> Vec<PartialTrack> {
    let mut next_id = 0usize;
    let mut tracks: Vec<PartialTrack> = Vec::new();

    if frames.is_empty() {
        return tracks;
    }

    // Initialize with first frame’s peaks
    {
        let init = &frames[0];
        let n = init.peaks.len().min(max_tracks);
        for p in init.peaks.iter().take(n) {
            tracks.push(PartialTrack {
                id: {
                    let id = next_id;
                    next_id += 1;
                    id
                },
                points: vec![(init.frame_idx, p.freq_hz, p.mag_db, p.snr_db)],
            });
        }
    }

    // Continue across remaining frames
    for frame in frames.iter().skip(1) {
        // Prepare availability flags
        let mut peak_used = vec![false; frame.peaks.len()];
        let mut track_used = vec![false; tracks.len()];

        // Build all candidate (track, peak) pairs with cost = |cents|
        let mut pairs: Vec<(usize, usize, f32)> = Vec::new();
        for (ti, tr) in tracks.iter().enumerate() {
            let &(last_frame, last_f, _, _) = tr.points.last().unwrap();
            let _ = last_frame; // unused, but kept for clarity

            for (pi, pk) in frame.peaks.iter().enumerate() {
                let cost_cents = cents_diff(last_f, pk.freq_hz).abs();
                if cost_cents <= max_jump_cents {
                    pairs.push((ti, pi, cost_cents));
                }
            }
        }

        // Sort by ascending cost (closest links first)
        pairs.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());

        // Greedy assignment
        for (ti, pi, _cost) in pairs {
            if track_used[ti] || peak_used[pi] {
                continue;
            }
            let pk = &frame.peaks[pi];
            tracks[ti]
                .points
                .push((frame.frame_idx, pk.freq_hz, pk.mag_db, pk.snr_db));
            track_used[ti] = true;
            peak_used[pi] = true;
        }

        // Start new tracks from any strong unused peaks (until max_tracks)
        for (pi, pk) in frame.peaks.iter().enumerate() {
            if peak_used[pi] {
                continue;
            }
            if tracks.len() >= max_tracks {
                break;
            }
            tracks.push(PartialTrack {
                id: {
                    let id = next_id;
                    next_id += 1;
                    id
                },
                points: vec![(frame.frame_idx, pk.freq_hz, pk.mag_db, pk.snr_db)],
            });
        }
    }

    tracks
}

/// Median utility over an iterator of f32 values.
/// Returns NaN if empty.
fn median_f32<I: Iterator<Item = f32>>(iter: I) -> f32 {
    let mut v: Vec<f32> = iter.collect();
    if v.is_empty() {
        return f32::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = v.len() / 2;
    if v.len() % 2 == 1 {
        v[mid]
    } else {
        0.5 * (v[mid - 1] + v[mid])
    }
}

/// Cents difference between two positive frequencies.
fn cents_diff(f1: f32, f2: f32) -> f32 {
    if f1 <= 0.0 || f2 <= 0.0 {
        return f32::INFINITY;
    }
    1200.0 * ((f2 / f1) as f64).log2() as f32
}
