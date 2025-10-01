use consonare_core::{
    diagnostics::{QualityConfig, run_quality_step},
    dissonance::{DissonanceConfig, run_dissonance_step},
    dissonance_n::{MultiDissonanceConfig, run_multi_dissonance_step},
    f0::{F0Config, run_f0_step},
    inharmonicity::{InharmonicityConfig, run_inharmonicity_step},
    input::{Config, process_file},
    intervals::{IntervalNamingConfig, run_interval_naming_step},
    musician::{
        NamedInterval, export_dissonance_csv, export_overtone_csv, overtone_table,
        tuning_suggestion_multi_note, tuning_suggestion_two_note,
    },
    partials::{PartialsConfig, run_partials_step},
    preprocess::{SpectralConfig, run_spectral_step},
    weighting::{WeightingConfig, run_weighting_step},
};

use crate::common::gen_target_fn;

mod common;
#[cfg(feature = "visualise")]
mod heatmap;
#[cfg(feature = "visualise")]
mod plot;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let fname = "keyboards_upright-piano_e4.wav";

    // ---- Step #1: steady-pitch preprocessing & checks ----
    let cfg = Config::default(); // tweak thresholds/window/hop here if needed
    let result = process_file(format!("samples/{fname}"), &cfg)?;

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

    // ---- Step #5: inharmonicity modelling (optional) ----
    //
    // Fit the stiff-string coefficient B from partials (if applicable) using
    // the linearized model r_n^2 = a + b n^2 with r_n = f_n / (n * f0_hat).
    // Accept only if the fit quality is good; otherwise skip.
    let ih_cfg = InharmonicityConfig {
        ..Default::default()
    };
    let f0_hat_opt = f0.consensus_f0_hz;
    let ih = run_inharmonicity_step(&parts, f0_hat_opt, &ih_cfg);

    println!("\n== Step 5 ==");
    if ih.b.is_none() {
        if let Some(f0h) = f0_hat_opt {
            println!(
                "No reliable stiff-string fit (used {} points, R^2={:.3}). Keeping empirical partials. f0≈{:.2} Hz",
                ih.used_points, ih.r2, f0h
            );
        } else {
            println!("No f0 available → cannot attempt inharmonicity fit.");
        }
    } else {
        let b_est = ih.b.unwrap();
        println!(
            "Accepted stiff-string model: B ≈ {:.6e}, R^2={:.3}, points={}",
            b_est, ih.r2, ih.used_points
        );
        println!(
            "Refined f0 ≈ {:.2} Hz (median abs error ≈ {:.2} cents)",
            ih.f0_refined_hz,
            ih.median_abs_cents.unwrap_or(0.0)
        );
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

        #[cfg(feature = "visualise")]
        {
            if let Err(e) = plot::plot_dissonance(&d_result, fname) {
                eprintln!("Plotting failed: {e}");
            } else {
                println!("Saved dissonance profile plot as dissonance.png");
            }
        }
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

    // ---- Step #9: multi-note optimisation (N ≥ 3) ----
    let m_cfg = MultiDissonanceConfig {
        ..Default::default()
    };
    let m_res = run_multi_dissonance_step(&weighted, None, &m_cfg);
    println!("\n== Step 9 ({}-note) ==", m_cfg.n_notes);
    for (idx, sol) in m_res.solutions.iter().enumerate() {
        println!("Solution #{idx}:  total D ≈ {:.6}", sol.total_roughness);
        for (k, r) in sol.ratios.iter().enumerate() {
            if k == 0 {
                println!("  Note {k}:  r=1.00000 (0.0 cents)  [reference]");
            } else if let Some(ap) = &r.approx {
                println!(
                    "  Note {k}:  r≈{:.5} ({:.1} cents)  ~ {}  |Δc|≈{:.2}c",
                    r.ratio, r.cents, ap.label, ap.cents_error
                );
            } else {
                println!("  Note {k}:  r≈{:.5} ({:.1} cents)", r.ratio, r.cents);
            }
        }
    }

    #[cfg(feature = "visualise")]
    {
        if m_cfg.n_notes == 3 {
            // Use the same sweep range as Step 7 by default.

            use crate::heatmap::{compute_dissonance_grid_3, plot_dissonance_surface_3};
            let grid = compute_dissonance_grid_3(
                &weighted,
                None, // or Some(&[b1, b2]) if you have measured spectra for the other 2 notes
                m_cfg.min_cents,
                m_cfg.max_cents,
                m_cfg.step_cents,
                m_cfg.max_deltaf_over_cbw,
            );
            if let Err(e) = plot_dissonance_surface_3(
                &grid,
                format!("Dissonance Surface of `{}` (N=3)", fname),
                fname,
            ) {
                eprintln!("Surface plotting failed: {e}");
            } else {
                println!("Saved dissonance surface as dissonance_3.png");
            }
        }
    }

    // ---- Step #10: quality control & diagnostics ----
    let q_cfg = QualityConfig {
        ..Default::default()
    };
    let _qres = run_quality_step(&parts, &f0, &d_result, &i_result, &ih, &q_cfg);

    // ---- Step #11: musician‑facing outputs ----
    println!("\n== Step 11 (Musician‑facing outputs) ==");
    // Choose the best f0 to reference: prefer refined if available.
    let f0_for_table: Option<f32> = if ih.b.is_some() {
        Some(ih.f0_refined_hz)
    } else {
        f0.consensus_f0_hz
    };
    if let Some(f_ref) = f0_for_table {
        println!("Reference f0 for overtone table: {:.3} Hz", f_ref);
    } else {
        println!("Reference f0 for overtone table: (none)");
    }

    // Collect weighted partials (freq_hz, median_db, weight)
    let wtriples: Vec<(f32, f32, f32)> = weighted
        .partials
        .iter()
        .map(|wp| (wp.freq_hz, wp.median_db, wp.weight))
        .collect();

    // Build overtone table (top 20 by frequency)
    let overtones = overtone_table(f0_for_table, &wtriples, 20);
    if overtones.is_empty() {
        println!("No overtones available to print.");
    } else {
        println!("Overtone table (first {}):", overtones.len());
        println!(" n   f_n (Hz)    Δcents    level(dB)  weight");
        for r in &overtones {
            let dc = r
                .cents_from_harm
                .map(|x| format!("{:+6.2}", x))
                .unwrap_or_else(|| "   --  ".into());
            println!(
                " {:>2}  {:>9.3}  {}   {:>8.2}  {:>6.3}",
                r.n, r.freq_hz, dc, r.amplitude_db, r.weight
            );
        }
    }

    // Dissonance CSV export (cents, D)
    // We don't assume d_result.smoothed carries the cents; reconstruct from cfg
    let mut dsamples: Vec<(f32, f32)> = Vec::with_capacity(d_result.smoothed.len());
    for (i, p) in d_result.smoothed.iter().enumerate() {
        let cents = d_cfg.min_cents + (i as f32) * d_cfg.step_cents;
        dsamples.push((cents, p.value));
    }
    if let Err(e) = export_dissonance_csv(
        gen_target_fn(format!("{}_dissonance-curve.csv", fname)),
        &dsamples,
    ) {
        eprintln!("Could not write dissonance_curve.csv: {e}");
    } else {
        println!("Saved dissonance_curve.csv (cents,D).");
    }
    if let Err(e) = export_overtone_csv(
        gen_target_fn(format!("{}_overtone-table.csv", fname)),
        &overtones,
    ) {
        eprintln!("Could not write overtone_table.csv: {e}");
    } else {
        println!("Saved overtone_table.csv.");
    }

    // Two‑note tuning suggestion from Step 8
    let mut named_for_step11 = Vec::new();
    for nm in &i_result.named {
        if let Some(best) = &nm.best {
            named_for_step11.push(NamedInterval {
                label: best.label.clone(),
                ratio: best.ratio,
                at_cents: nm.at_cents,
                cents_error: best.cents_error,
                confidence: best.confidence,
            });
        }
    }
    if let Some(msg) = tuning_suggestion_two_note(&named_for_step11) {
        println!("\n{}", msg);
    } else {
        println!("\nNo two‑note tuning suggestion (no named minima).");
    }

    // N‑note suggestion (Step 9): use the first solution, if any
    let multi_for_step11 = m_res.solutions.first().map(|sol| {
        sol.ratios
            .iter()
            .map(|r| {
                if let Some(ap) = &r.approx {
                    (
                        r.ratio,
                        r.cents,
                        Some(ap.label.clone()),
                        Some(ap.cents_error),
                    )
                } else {
                    (r.ratio, r.cents, None, None)
                }
            })
            .collect::<Vec<_>>()
    });

    if let Some(summary) = tuning_suggestion_multi_note(multi_for_step11.as_deref()) {
        println!("{}", summary);
    } else {
        println!("No N‑note tuning summary available.");
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
