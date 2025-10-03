use clap::{ArgAction, Parser, ValueEnum};
use std::path::{Path, PathBuf};

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

use crate::{
    common::gen_target_fn,
    util::{mean, nearest_shift, summary_min_median_max, write_csv},
};

mod common;
#[cfg(feature = "visualise")]
mod heatmap;
#[cfg(feature = "visualise")]
mod plot;
pub mod util;

/// Print only when `verbose` is true.
macro_rules! vprintln {
    ($v:expr, $($arg:tt)*) => {
        if $v { println!($($arg)*); }
    };
}

/// CLI for consonare analysis.
#[derive(Parser, Debug)]
#[command(
    name = "consonare-cli",
    about = "Analyze a sustained note from an audio file (F0, partials, dissonance, intervals) and print musician-facing guidance.",
    version
)]
struct Args {
    /// Path to the audio file (e.g., samples/upright-piano_e4.wav)
    #[arg(value_name = "FILE_PATH")]
    file_path: PathBuf,

    /// Increase verbosity (show diagnostics for Steps 1–10). `-v` or `--verbose`.
    #[arg(short, long, action = ArgAction::Count)]
    verbose: u8,

    /// Output directory for CSV exports.
    #[arg(long, value_name = "DIR", default_value = ".")]
    out_dir: PathBuf,

    /// Which CSVs to export.
    #[arg(long, value_enum, default_value_t = Export::All)]
    export: Export,

    /// Ablations to run (repeatable). Example: --ablate weighting --ablate masking
    #[arg(long, value_enum, action = ArgAction::Append)]
    ablate: Vec<Ablate>,

    /// Export ablation summaries as CSV alongside other outputs.
    #[arg(long, default_value_t = true)]
    ablation_report: bool,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
enum Export {
    None,
    Dissonance,
    Overtones,
    All,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
enum Ablate {
    /// Disable perceptual weighting (weights=1, A-weight=0 dB)
    Weighting,
    /// Disable ERB masking (no partials hidden by masking)
    Masking,
    /// Disable ERB/CBW normalization in roughness kernel
    ErbNorm,
    /// Apply all of the above
    All,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // ---- CLI ----
    let args = Args::parse();
    let verbose = args.verbose > 0;
    let audio_path: &Path = &args.file_path;
    let display_name = audio_path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("audio");
    let out_dir = args.out_dir;

    // Use file stem for output file prefixes (nicer than including .wav in the name).
    let out_prefix = audio_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("output");

    // ---- Step #1: steady-pitch preprocessing & checks ----
    let cfg = Config::default(); // tweak thresholds/window/hop here if needed
    let result = process_file(audio_path.to_string_lossy().to_string(), &cfg)?;

    vprintln!(verbose, "== Step 1 ==");
    vprintln!(verbose, "Sample rate: {} Hz", result.analysis.sample_rate);
    vprintln!(
        verbose,
        "Trimmed {} -> {} samples",
        result.analysis.original_samples,
        result.analysis.trimmed_samples
    );
    vprintln!(
        verbose,
        "Applied gain: {:.2} dB",
        result.analysis.applied_gain_db
    );

    if let (Some(med), Some(max_dev)) = (
        result.analysis.median_f0_hz,
        result.analysis.max_abs_dev_cents,
    ) {
        vprintln!(
            verbose,
            "Median f0: {:.2} Hz, Max deviation: {:.1} cents",
            med,
            max_dev
        );
        vprintln!(verbose, "Steady pitch? {}", result.analysis.steady_pitch);
        if !verbose {
            println!(
                "Loaded `{}` • median f0 ≈ {:.2} Hz • steady_pitch={}",
                display_name, med, result.analysis.steady_pitch
            );
        }
    } else {
        vprintln!(
            verbose,
            "No reliable f0 detected; steady pitch? {}",
            result.analysis.steady_pitch
        );
        if !verbose {
            println!(
                "Loaded `{}` • f0 not reliable • steady_pitch={}",
                display_name, result.analysis.steady_pitch
            );
        }
    }

    if result.samples.is_empty() {
        eprintln!("No samples remaining after trim; aborting.");
        return Ok(());
    }

    // ---- Step #2: spectral preprocessing (band-limit, STFT, median, floor/mask) ----
    let sr = result.analysis.sample_rate;

    // Try larger → smaller windows until we get at least 3 frames.
    let mut spec = None;
    for (win, hop, medlen) in [
        (16384usize, 4096usize, 5usize), // ~0.37 s window @44.1k
        (8192, 2048, 5),                 // ~0.19 s
        (4096, 1024, 3),                 // ~0.09 s
    ] {
        let cfg_try = SpectralConfig {
            win_size: win,
            hop_size: hop,
            median_len: medlen,
            ..Default::default()
        };
        let s = run_spectral_step(&result.samples, sr, &cfg_try)?;
        if s.times_s.len() >= 3 {
            spec = Some(s);
            break;
        }
    }
    // Fallback if even the smallest window failed (extremely short / heavy trim)
    let spec = spec.unwrap_or_else(|| {
        run_spectral_step(
            &result.samples,
            sr,
            &SpectralConfig {
                win_size: 4096,
                hop_size: 512,
                median_len: 3,
                ..Default::default()
            },
        )
        .expect("spectral step (fallback)")
    });

    vprintln!(verbose, "\n== Step 2 ==");
    vprintln!(
        verbose,
        "Frames: {} | Bins: {}",
        spec.times_s.len(),
        spec.freqs_hz.len()
    );
    if let (Some(t0), Some(t1)) = (spec.times_s.first(), spec.times_s.last()) {
        vprintln!(verbose, "Time span: {:.3} … {:.3} s", t0, t1);
    }
    if let (Some(f0), Some(f1)) = (spec.freqs_hz.first(), spec.freqs_hz.last()) {
        vprintln!(verbose, "Freq span: {:.1} … {:.1} Hz", f0, f1);
    }

    // Suppression stats
    let (mut suppressed, mut total) = (0usize, 0usize);
    for row in &spec.suppressed {
        total += row.len();
        suppressed += row.iter().filter(|&&b| b).count();
    }
    if total > 0 {
        let pct = 100.0 * suppressed as f64 / total as f64;
        vprintln!(verbose, "Suppressed bins: {suppressed}/{total} ({pct:.1}%)");
    }

    // Noise-floor summary (dB)
    let (min_floor, med_floor, max_floor) = summary_min_median_max(&spec.floor_db);
    vprintln!(
        verbose,
        "Noise floor (dB): min/median/max = {:.1} / {:.1} / {:.1}",
        min_floor,
        med_floor,
        max_floor
    );

    // Example frames
    for &probe_idx in &[0usize, spec.mag_db.len().saturating_sub(1)] {
        if let Some(row) = spec.mag_db.get(probe_idx) {
            let (min_db, med_db, max_db) = summary_min_median_max(row);
            vprintln!(
                verbose,
                "Frame {probe_idx}: smoothed dB min/median/max = {:.1} / {:.1} / {:.1}",
                min_db,
                med_db,
                max_db
            );
        }
    }

    // ---- Step #3: fundamental estimation ----
    let mut f0_cfg = F0Config {
        // align with the STFT grid (helps agreement/consensus too)
        yin_frame: spec.freqs_hz.len().saturating_sub(1) * 2, // ≈ win_size
        yin_hop: Some(if spec.times_s.len() > 1 {
            ((spec.times_s[1] - spec.times_s[0]) * sr as f32) as usize
        } else {
            1024
        }),
        ..Default::default()
    };

    if let Some(med) = result.analysis.median_f0_hz {
        // A conservative but effective band that avoids 50/60 Hz subharmonics:
        f0_cfg.fmin_hz = (med * 0.6).max(50.0);
        f0_cfg.fmax_hz = med * 1.8;
    }

    let f0 = run_f0_step(&result.samples, sr, Some(&spec), &f0_cfg)?;

    vprintln!(verbose, "\n== Step 3 ==");
    vprintln!(verbose, "Analyzed {} frames", f0.frames.len());
    if let Some(hz) = f0.consensus_f0_hz {
        vprintln!(verbose, "Consensus f0: {:.3} Hz", hz);
        vprintln!(verbose, "Confidence: {:.2}", f0.consensus_confidence);
        if !verbose {
            println!("f0 ≈ {:.3} Hz (conf {:.2})", hz, f0.consensus_confidence);
        }
    } else {
        vprintln!(verbose, "Consensus f0: (none)");
        if !verbose {
            println!("f0: (none)");
        }
    }
    if let (Some(first), Some(last)) = (f0.frames.first(), f0.frames.last()) {
        vprintln!(
            verbose,
            "First frame: yin={:.3?} Hz, cep={:.3?} Hz, fused={:.3?} Hz, agree={:.2}",
            first.yin_hz,
            first.cep_hz,
            first.fused_hz,
            first.agreement
        );
        vprintln!(
            verbose,
            "Last frame:  yin={:.3?} Hz, cep={:.3?} Hz, fused={:.3?} Hz, agree={:.2}",
            last.yin_hz,
            last.cep_hz,
            last.fused_hz,
            last.agreement
        );
    }

    // ---- Step #4: partial extraction ----
    let p_cfg = {
        let mut cfg = PartialsConfig {
            ..Default::default()
        };
        if spec.times_s.len() < 4 {
            cfg.min_track_len = 2; // accept 2-frame tracks when the signal is short
            cfg.track_max_jump_cents = 50.0; // slightly easier linking across coarse hops
        }
        cfg
    };
    let consensus_f0_hz_f32 = f0.consensus_f0_hz;
    let parts = run_partials_step(&spec, consensus_f0_hz_f32, &p_cfg)?;

    vprintln!(verbose, "\n== Step 4 ==");
    vprintln!(verbose, "Tracks formed: {}", parts.tracks.len());
    vprintln!(
        verbose,
        "Reported summaries (len >= {}): {}",
        p_cfg.min_track_len,
        parts.summaries.len()
    );

    // Concise table of top partials (by frequency)
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
        vprintln!(
            verbose,
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

    // ---- Step #5: inharmonicity modelling ----
    let ih_cfg = InharmonicityConfig {
        ..Default::default()
    };
    let f0_hat_opt = f0.consensus_f0_hz;
    let ih = run_inharmonicity_step(&parts, f0_hat_opt, &ih_cfg);

    vprintln!(verbose, "\n== Step 5 ==");
    if ih.b.is_none() {
        if let Some(f0h) = f0_hat_opt {
            vprintln!(
                verbose,
                "No reliable stiff-string fit (used {} points, R^2={:.3}). Keeping empirical partials. f0≈{:.2} Hz",
                ih.used_points,
                ih.r2,
                f0h
            );
        } else {
            vprintln!(
                verbose,
                "No f0 available → cannot attempt inharmonicity fit."
            );
        }
    } else {
        let b_est = ih.b.unwrap();
        vprintln!(
            verbose,
            "Accepted stiff-string model: B ≈ {:.6e}, R^2={:.3}, points={}",
            b_est,
            ih.r2,
            ih.used_points
        );
        vprintln!(
            verbose,
            "Refined f0 ≈ {:.2} Hz (median abs error ≈ {:.2} cents)",
            ih.f0_refined_hz,
            ih.median_abs_cents.unwrap_or(0.0)
        );
    }

    // ---- Step #6: perceptual weighting ----
    let w_cfg = WeightingConfig {
        ..Default::default()
    };
    let weighted = run_weighting_step(&parts, &w_cfg);

    vprintln!(verbose, "\n== Step 6 ==");
    if weighted.partials.is_empty() {
        vprintln!(verbose, "No partials to weight.");
    } else {
        vprintln!(verbose, "Anchor partial: {:?}", weighted.anchor_idx);
        vprintln!(verbose, "f (Hz)     A(dB)   med(dB)   masked   weight");
        for wp in &weighted.partials {
            vprintln!(
                verbose,
                "{:>8.3}  {:>6.2}   {:>7.2}   {:>5}   {:>6.3}",
                wp.freq_hz,
                wp.a_weight_db,
                wp.median_db,
                if wp.erb_masked { "yes" } else { " no " },
                wp.weight
            );
        }
    }

    // ---- Step #7: dissonance profile sweep ----
    let d_cfg = DissonanceConfig {
        ..Default::default()
    };
    let d_result = run_dissonance_step(&weighted, None, &d_cfg);

    vprintln!(verbose, "\n== Step 7 ==");
    if d_result.smoothed.is_empty() {
        vprintln!(verbose, "No dissonance profile (empty partial sets).");
    } else {
        vprintln!(verbose, "Top local minima (by D value):");
        for m in &d_result.minima {
            let ratio = (2.0_f32).powf(m.cents / 1200.0);
            vprintln!(
                verbose,
                "  {:>7.1} cents  (ratio ≈ {:.5})   D ≈ {:>8.6}   depth ≈ {:>7.6}",
                m.cents,
                ratio,
                m.value,
                m.depth
            );
        }
        let vals: Vec<f32> = d_result.smoothed.iter().map(|p| p.value).collect();
        let (dmin, dmed, dmax) = summary_min_median_max(&vals);
        vprintln!(
            verbose,
            "D(r) stats (smoothed): min={:.6}  med={:.6}  max={:.6}",
            dmin,
            dmed,
            dmax
        );

        #[cfg(feature = "visualise")]
        {
            if let Err(e) = plot::plot_dissonance(&d_result, out_prefix) {
                eprintln!("Plotting failed: {e}");
            } else {
                vprintln!(verbose, "Saved dissonance profile plot as dissonance.png");
            }
        }
    }

    // ---- Step #8: interval selection & naming ----
    let i_cfg = IntervalNamingConfig {
        ..Default::default()
    };
    let i_result = run_interval_naming_step(&d_result, &i_cfg);

    vprintln!(verbose, "\n== Step 8 ==");
    for (k, nm) in i_result.named.iter().enumerate() {
        vprintln!(
            verbose,
            "Min #{} at {:.2} cents: D={:.6} (depth {:.6}, sharpness {:.6})",
            k + 1,
            nm.at_cents,
            nm.at_value,
            nm.depth,
            nm.sharpness
        );
        if let Some(best) = &nm.best {
            vprintln!(
                verbose,
                "  → best: {}  (~{:.5})  err={:+.2} cents  J={:.4}  conf={:.2}",
                best.label,
                best.ratio,
                best.cents_error,
                best.score,
                best.confidence
            );
        } else {
            vprintln!(
                verbose,
                "  → no rational within ±{:.1} cents",
                i_cfg.tolerance_cents
            );
        }
        if !nm.alternatives.is_empty() {
            vprintln!(verbose, "    alternatives:");
            for alt in &nm.alternatives {
                vprintln!(
                    verbose,
                    "      - {}/{}  (~{:.5})  err={:+.2}c  J={:.4}",
                    alt.p,
                    alt.q,
                    alt.ratio,
                    alt.cents_error,
                    alt.score
                );
            }
        }
    }

    // ---- Step #9: multi-note optimisation ----
    let m_cfg = MultiDissonanceConfig {
        ..Default::default()
    };
    let m_res = run_multi_dissonance_step(&weighted, None, &m_cfg);

    vprintln!(verbose, "\n== Step 9 ({}-note) ==", m_cfg.n_notes);
    for (idx, sol) in m_res.solutions.iter().enumerate() {
        vprintln!(
            verbose,
            "Solution #{idx}:  total D ≈ {:.6}",
            sol.total_roughness
        );
        for (k, r) in sol.ratios.iter().enumerate() {
            if k == 0 {
                vprintln!(verbose, "  Note {k}:  r=1.00000 (0.0 cents)  [reference]");
            } else if let Some(ap) = &r.approx {
                vprintln!(
                    verbose,
                    "  Note {k}:  r≈{:.5} ({:.1} cents)  ~ {}  |Δc|≈{:.2}c",
                    r.ratio,
                    r.cents,
                    ap.label,
                    ap.cents_error
                );
            } else {
                vprintln!(
                    verbose,
                    "  Note {k}:  r≈{:.5} ({:.1} cents)",
                    r.ratio,
                    r.cents
                );
            }
        }
    }

    #[cfg(feature = "visualise")]
    {
        if m_cfg.n_notes == 3 {
            use crate::heatmap::{compute_dissonance_grid_3, plot_dissonance_surface_3};
            let grid = compute_dissonance_grid_3(
                &weighted,
                None,
                m_cfg.min_cents,
                m_cfg.max_cents,
                m_cfg.step_cents,
                m_cfg.max_deltaf_over_cbw,
            );
            if let Err(e) = plot_dissonance_surface_3(&grid, out_prefix) {
                eprintln!("Surface plotting failed: {e}");
            } else {
                vprintln!(verbose, "Saved dissonance surface as dissonance_3.png");
            }
        }
    }

    // ---- Step #10: quality control & diagnostics ----
    let q_cfg = QualityConfig {
        ..Default::default()
    };
    let _qres = run_quality_step(&parts, &f0, &d_result, &i_result, &ih, &q_cfg);

    // ---- Step #11: musician-facing outputs (ALWAYS printed) ----
    println!("\n== Step 11 (Musician-facing outputs) ==");
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

    // Prepare dissonance CSV samples (cents, D)
    let mut dsamples: Vec<(f32, f32)> = Vec::with_capacity(d_result.smoothed.len());
    for (i, p) in d_result.smoothed.iter().enumerate() {
        let cents = d_cfg.min_cents + (i as f32) * d_cfg.step_cents;
        dsamples.push((cents, p.value));
    }

    // CSV export control
    if matches!(args.export, Export::All | Export::Dissonance) {
        let path = format!("{}/{}_dissonance-curve.csv", out_dir.display(), out_prefix);
        if let Err(e) = export_dissonance_csv(gen_target_fn(path), &dsamples) {
            eprintln!("Could not write dissonance_curve.csv: {e}");
        } else {
            println!("Saved dissonance_curve.csv.");
        }
    }
    if matches!(args.export, Export::All | Export::Overtones) {
        let path = format!("{}/{}_overtone-table.csv", out_dir.display(), out_prefix);
        if let Err(e) = export_overtone_csv(gen_target_fn(path), &overtones) {
            eprintln!("Could not write overtone_table.csv: {e}");
        } else {
            println!("Saved overtone_table.csv.");
        }
    }

    // Two-note tuning suggestion from Step 8
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
        println!("\nNo two-note tuning suggestion (no named minima).");
    }

    // N-note suggestion (Step 9): use the first solution, if any
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
        println!("No N-note tuning summary available.");
    }

    // ---- Step 12: Ablation runs & metrics (Outcome #9) ----
    let do_any_ablate = !args.ablate.is_empty();
    if do_any_ablate {
        println!("\n== Step 12 (Ablation study) ==");

        // Baseline artifacts we’ll compare against
        let baseline_named = i_result.named.clone();
        let baseline_cents: Vec<f32> = baseline_named.iter().map(|nm| nm.at_cents).collect();
        let baseline_depths: Vec<f32> = baseline_named.iter().map(|nm| nm.depth).collect();
        let baseline_sharp: Vec<f32> = baseline_named.iter().map(|nm| nm.sharpness).collect();
        let baseline_conf: Vec<f32> = baseline_named
            .iter()
            .filter_map(|nm| nm.best.as_ref().map(|b| b.confidence))
            .collect();

        // Utility to run one variant with toggles applied
        let run_variant = |label: &str,
                           ablate_weighting: bool,
                           ablate_masking: bool,
                           ablate_erb: bool| {
            // Step 6 (recompute weighting with ablation flags)
            let mut w_cfg_ab = WeightingConfig {
                ..Default::default()
            };
            if ablate_masking {
                w_cfg_ab.enable_masking = false;
            }
            if ablate_weighting {
                w_cfg_ab.enable_perceptual_weighting = false;
            }
            let weighted_ab = run_weighting_step(&parts, &w_cfg_ab);

            // Step 7 (dissonance with or without ERB normalization)
            let mut d_cfg_ab = DissonanceConfig {
                ..Default::default()
            };
            d_cfg_ab.normalize_by_cbw = !ablate_erb;
            let d_res_ab = run_dissonance_step(&weighted_ab, None, &d_cfg_ab);

            // Step 8 (naming)
            let i_res_ab = run_interval_naming_step(&d_res_ab, &i_cfg);

            // Metrics
            let vals_ab: Vec<f32> = d_res_ab.smoothed.iter().map(|p| p.value).collect();
            let (dmin, dmed, dmax) = summary_min_median_max(&vals_ab);
            let named_ab = &i_res_ab.named;
            let cents_ab: Vec<f32> = named_ab.iter().map(|nm| nm.at_cents).collect();
            let depths_ab: Vec<f32> = named_ab.iter().map(|nm| nm.depth).collect();
            let sharp_ab: Vec<f32> = named_ab.iter().map(|nm| nm.sharpness).collect();
            let conf_ab: Vec<f32> = named_ab
                .iter()
                .filter_map(|nm| nm.best.as_ref().map(|b| b.confidence))
                .collect();

            let (mean_shift, shifts) = nearest_shift(&baseline_cents, &cents_ab);

            // Print one-line summary
            println!(
                "{:<12} | Dmin={:.6}  Dmed={:.6}  mins={}  depth_top={:.6}  sharp_mean={:.6}  conf_mean={:.3}  Δcents_vs_base={:.2}",
                label,
                dmin,
                dmed,
                named_ab.len(),
                depths_ab.iter().cloned().fold(0.0, f32::max),
                mean(&sharp_ab),
                mean(&conf_ab),
                mean_shift
            );

            (
                label.to_string(),
                dmin,
                dmed,
                dmax,
                named_ab.len() as i32,
                depths_ab.into_iter().fold(0.0, f32::max),
                mean(&sharp_ab),
                mean(&conf_ab),
                mean_shift,
                shifts,
                cents_ab,
            )
        };

        // Determine which variants to run
        let ablate_all = args.ablate.contains(&Ablate::All);
        let want_w = ablate_all || args.ablate.contains(&Ablate::Weighting);
        let want_m = ablate_all || args.ablate.contains(&Ablate::Masking);
        let want_e = ablate_all || args.ablate.contains(&Ablate::ErbNorm);

        // Always include a baseline row at top of CSV (already printed in Steps 7-8)
        let vals_base: Vec<f32> = d_result.smoothed.iter().map(|p| p.value).collect();
        let (bdmin, bdmed, bdmax) = summary_min_median_max(&vals_base);
        println!(
            "{:<12} | Dmin={:.6}  Dmed={:.6}  mins={}  depth_top={:.6}  sharp_mean={:.6}  conf_mean={:.3}  Δcents_vs_base={:.2}",
            "baseline",
            bdmin,
            bdmed,
            baseline_named.len(),
            baseline_depths.iter().cloned().fold(0.0, f32::max),
            mean(&baseline_sharp),
            mean(&baseline_conf),
            0.0
        );

        // Collect rows for CSV export
        let mut summary_rows = Vec::<String>::new();
        summary_rows.push("variant,Dmin,Dmed,Dmax,num_minima,top_depth,mean_sharpness,mean_confidence,mean_shift_vs_baseline".into());
        summary_rows.push(format!(
            "baseline,{:.6},{:.6},{:.6},{}, {:.6},{:.6},{:.6},{:.2}",
            bdmin,
            bdmed,
            bdmax,
            baseline_named.len(),
            baseline_depths.iter().cloned().fold(0.0, f32::max),
            mean(&baseline_sharp),
            mean(&baseline_conf),
            0.0
        ));

        let mut minima_rows = Vec::<String>::new();
        minima_rows.push("variant,baseline_min_index,baseline_at_cents,nearest_shift_cents".into());

        // Run requested variants (independent toggles)
        if want_w {
            let (name, dmin, dmed, dmax, nmins, topd, msh, mc, mshift, shifts, _cents) =
                run_variant("no-weight", true, false, false);
            summary_rows.push(format!(
                "{},{:.6},{:.6},{:.6},{},{:.6},{:.6},{:.6},{:.2}",
                name, dmin, dmed, dmax, nmins, topd, msh, mc, mshift
            ));
            for (k, s) in shifts.into_iter().enumerate() {
                minima_rows.push(format!("{},{},{:.3},{:.3}", name, k, baseline_cents[k], s));
            }
        }
        if want_m {
            let (name, dmin, dmed, dmax, nmins, topd, msh, mc, mshift, shifts, _cents) =
                run_variant("no-mask", false, true, false);
            summary_rows.push(format!(
                "{},{:.6},{:.6},{:.6},{},{:.6},{:.6},{:.6},{:.2}",
                name, dmin, dmed, dmax, nmins, topd, msh, mc, mshift
            ));
            for (k, s) in shifts.into_iter().enumerate() {
                minima_rows.push(format!("{},{},{:.3},{:.3}", name, k, baseline_cents[k], s));
            }
        }
        if want_e {
            let (name, dmin, dmed, dmax, nmins, topd, msh, mc, mshift, shifts, _cents) =
                run_variant("no-erb", false, false, true);
            summary_rows.push(format!(
                "{},{:.6},{:.6},{:.6},{},{:.6},{:.6},{:.6},{:.2}",
                name, dmin, dmed, dmax, nmins, topd, msh, mc, mshift
            ));
            for (k, s) in shifts.into_iter().enumerate() {
                minima_rows.push(format!("{},{},{:.3},{:.3}", name, k, baseline_cents[k], s));
            }
        }

        // CSV export
        if args.ablation_report {
            let summary_path = format!("{}/{}_ablation_summary.csv", out_dir.display(), out_prefix);
            let minima_path = format!(
                "{}/{}_ablation_minima_shift.csv",
                out_dir.display(),
                out_prefix
            );
            if let Err(e) = write_csv(&summary_path, &summary_rows) {
                eprintln!("Could not write ablation_summary.csv: {e}");
            } else {
                println!("Saved ablation_summary.csv.");
            }
            if let Err(e) = write_csv(&minima_path, &minima_rows) {
                eprintln!("Could not write ablation_minima_shift.csv: {e}");
            } else {
                println!("Saved ablation_minima_shift.csv.");
            }
        }
    }

    Ok(())
}
