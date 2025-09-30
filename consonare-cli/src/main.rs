use consonare_core::dissonance::{DissonanceConfig, run_dissonance_step};
use consonare_core::{
    f0::{F0Config, run_f0_step},
    input::{Config, process_file},
    intervals::{IntervalNamingConfig, run_interval_naming_step},
    partials::{PartialsConfig, run_partials_step},
    preprocess::{SpectralConfig, run_spectral_step},
    weighting::{WeightingConfig, run_weighting_step},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // ---- Step #1: steady-pitch preprocessing & checks ----
    let cfg = Config::default(); // tweak thresholds/window/hop here if needed
    let result = process_file("samples/piano-c4-1s.wav", &cfg)?;
    //let result = process_file("samples/violin-c5-1s.wav", &cfg)?;

    println!("== Step 1 ==");
    println!("Sample rate: {} Hz", result.analysis.sample_rate);
    println!(
        "Trimmed {} -> {} samples",
        result.analysis.original_samples, result.analysis.trimmed_samples
    );
    println!("Applied gain: {:.2} dB", result.analysis.applied_gain_db);

    if let (Some(med), Some(max_dev)) = (
        result.analysis.median_f0_hz,
        result.analysis.max_abs_dev_cents,
    ) {
        println!(
            "Median f0: {:.2} Hz, Max deviation: {:.1} cents",
            med, max_dev
        );
        println!("Steady pitch? {}", result.analysis.steady_pitch);
    } else {
        println!(
            "No reliable f0 detected; steady pitch? {}",
            result.analysis.steady_pitch
        );
    }

    if result.samples.is_empty() {
        eprintln!("No samples remaining after trim; skipping Step #2.");
        return Ok(());
    }

    // ---- Step #2: spectral preprocessing (band-limit, STFT, median, floor/mask) ----
    let sr = result.analysis.sample_rate;
    let spec_cfg = SpectralConfig {
        ..Default::default()
    };

    let spec = run_spectral_step(&result.samples, sr, &spec_cfg)?;

    println!("\n== Step 2 ==");
    println!(
        "Frames: {} | Bins: {}",
        spec.times_s.len(),
        spec.freqs_hz.len()
    );
    if let (Some(t0), Some(t1)) = (spec.times_s.first(), spec.times_s.last()) {
        println!("Time span: {:.3} … {:.3} s", t0, t1);
    }
    if let (Some(f0), Some(f1)) = (spec.freqs_hz.first(), spec.freqs_hz.last()) {
        println!("Freq span: {:.1} … {:.1} Hz", f0, f1);
    }

    // Suppression stats: how much of the spectrogram is considered sub-floor
    let (mut suppressed, mut total) = (0usize, 0usize);
    for row in &spec.suppressed {
        total += row.len();
        suppressed += row.iter().filter(|&&b| b).count();
    }
    if total > 0 {
        let pct = 100.0 * suppressed as f64 / total as f64;
        println!("Suppressed bins: {suppressed}/{total} ({pct:.1}%)");
    }

    // Noise-floor summary (dB)
    let (min_floor, med_floor, max_floor) = summary_min_median_max(&spec.floor_db);
    println!(
        "Noise floor (dB): min/median/max = {:.1} / {:.1} / {:.1}",
        min_floor, med_floor, max_floor
    );

    // Example: peek at a couple of frame values (smoothed dB after median)
    for &probe_idx in &[0usize, spec.mag_db.len().saturating_sub(1)] {
        if let Some(row) = spec.mag_db.get(probe_idx) {
            let (min_db, med_db, max_db) = summary_min_median_max(row);
            println!(
                "Frame {probe_idx}: smoothed dB min/median/max = {:.1} / {:.1} / {:.1}",
                min_db, med_db, max_db
            );
        }
    }

    // ---- Step #3: fundamental estimation (YIN + cepstral + harmonic consensus) ----
    let f0_cfg = F0Config {
        ..Default::default()
    };

    let f0 = run_f0_step(&result.samples, sr, Some(&spec), &f0_cfg)?;

    println!("\n== Step 3 ==");
    println!("Analyzed {} frames", f0.frames.len());
    if let Some(hz) = f0.consensus_f0_hz {
        println!("Consensus f0: {:.3} Hz", hz);
        println!("Confidence: {:.2}", f0.consensus_confidence);
    } else {
        println!("Consensus f0: (none)");
    }
    if let (Some(first), Some(last)) = (f0.frames.first(), f0.frames.last()) {
        println!(
            "First frame: yin={:.3?} Hz, cep={:.3?} Hz, fused={:.3?} Hz, agree={:.2}",
            first.yin_hz, first.cep_hz, first.fused_hz, first.agreement
        );
        println!(
            "Last frame:  yin={:.3?} Hz, cep={:.3?} Hz, fused={:.3?} Hz, agree={:.2}",
            last.yin_hz, last.cep_hz, last.fused_hz, last.agreement
        );
    }

    // ---- Step #4: partial extraction (peaks -> tracks -> medians) ----
    let p_cfg = PartialsConfig {
        // You can tweak these to be stricter/looser:
        // max_partials: 20,
        // min_snr_db: 15.0,
        // local_max_half_width: 1,
        // track_max_jump_cents: 35.0,
        // min_track_len: 3,
        ..Default::default()
    };

    let consensus_f0_hz_f32 = f0.consensus_f0_hz;
    let parts = run_partials_step(&spec, consensus_f0_hz_f32, &p_cfg)?;

    println!("\n== Step 4 ==");
    println!("Tracks formed: {}", parts.tracks.len());
    println!(
        "Reported summaries (len >= {}): {}",
        p_cfg.min_track_len,
        parts.summaries.len()
    );

    // Print a concise table of top partials (by frequency)
    let mut shown = 0usize;
    for s in &parts.summaries {
        let harm = if let Some(f0hz) = consensus_f0_hz_f32 {
            if f0hz > 0.0 {
                (s.median_hz / f0hz).round() as i32
            } else {
                0
            }
        } else {
            0
        };

        println!(
            "k={:<2}  f≈{:>9.3} Hz  level≈{:>6.1} dB  SNR≈{:>5.1} dB  n={:<3}{}",
            s.track_id,
            s.median_hz,
            s.median_db,
            s.median_snr_db,
            s.count,
            if harm > 0 {
                format!("  (≈ harmonic {})", harm)
            } else {
                String::new()
            }
        );

        shown += 1;
        if shown >= 12 {
            break;
        }
    }

    // ---- Step #6: perceptual weighting (A-weighting + within-note masking) ----

    let w_cfg = WeightingConfig {
        // Tune these if desired:
        // - `mask_drop_db` (default 40 dB): drop very weak partials relative to the strongest (A-weighted) partial.
        // - `mask_margin_db` (default 30 dB) & `erb_mask_atten` (default 0.25): partials within the ERB of a much stronger neighbor are attenuated.
        // - `max_partials` caps how many partials from Step 4 summaries are considered.
        // - `erb_mask_atten`: 0.25,
        ..Default::default()
    };
    let weighted = run_weighting_step(&parts, &w_cfg);

    println!("\n== Step 6 ==");
    if weighted.partials.is_empty() {
        println!("No partials to weight.");
    } else {
        println!("Anchor partial: {:?}", weighted.anchor_idx);
        println!("f (Hz)     A(dB)   med(dB)   masked   weight");
        for wp in &weighted.partials {
            println!(
                "{:>8.3}  {:>6.2}   {:>7.2}   {:>5}   {:>6.3}",
                wp.freq_hz,
                wp.a_weight_db,
                wp.median_db,
                if wp.erb_masked { "yes" } else { " no " },
                wp.weight
            );
        }
    }

    // ---- Step #7: dissonance profile sweep over interval (cents) ----
    let d_cfg = DissonanceConfig {
        ..Default::default()
    };
    let d_result = run_dissonance_step(&weighted, None, &d_cfg);

    println!("\n== Step 7 ==");
    if d_result.smoothed.is_empty() {
        println!("No dissonance profile (empty partial sets).");
    } else {
        println!("Top local minima (by D value):");
        for m in &d_result.minima {
            let ratio = (2.0_f32).powf(m.cents / 1200.0);
            println!(
                "  {:>7.1} cents  (ratio ≈ {:.5})   D ≈ {:>8.6}   depth ≈ {:>7.6}",
                m.cents, ratio, m.value, m.depth
            );
        }
        let vals: Vec<f32> = d_result.smoothed.iter().map(|p| p.value).collect();
        let (dmin, dmed, dmax) = summary_min_median_max(&vals);
        println!(
            "D(r) stats (smoothed): min={:.6}  med={:.6}  max={:.6}",
            dmin, dmed, dmax
        );
    }

    // ---- Step #8: interval selection & naming near D(r) minima ----
    let i_cfg = IntervalNamingConfig {
        ..Default::default()
    };
    let i_result = run_interval_naming_step(&d_result, &i_cfg);

    println!("\n== Step 8 ==");
    for (k, nm) in i_result.named.iter().enumerate() {
        println!(
            "Min #{} at {:.2} cents: D={:.6} (depth {:.6}, sharpness {:.6})",
            k + 1,
            nm.at_cents,
            nm.at_value,
            nm.depth,
            nm.sharpness
        );
        if let Some(best) = &nm.best {
            println!(
                "  → best: {}  (~{:.5})  err={:+.2} cents  J={:.4}  conf={:.2}",
                best.label, best.ratio, best.cents_error, best.score, best.confidence
            );
        } else {
            println!("  → no rational within ±{:.1} cents", i_cfg.tolerance_cents);
        }
        if !nm.alternatives.is_empty() {
            println!("    alternatives:");
            for alt in &nm.alternatives {
                println!(
                    "      - {}/{}  (~{:.5})  err={:+.2}c  J={:.4}",
                    alt.p, alt.q, alt.ratio, alt.cents_error, alt.score
                );
            }
        }
    }

    Ok(())
}

fn summary_min_median_max(xs: &[f32]) -> (f32, f32, f32) {
    if xs.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let mut v = xs.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let min = v[0];
    let max = *v.last().unwrap();
    let mid = v.len() / 2;
    let med = if v.len() % 2 == 1 {
        v[mid]
    } else {
        0.5 * (v[mid - 1] + v[mid])
    };
    (min, med, max)
}
