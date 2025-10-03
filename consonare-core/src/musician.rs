use std::collections::HashMap;
// Step 11: Musician-facing outputs (overtone table, dissonance export, tuning suggestions)
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone, Copy)]
pub struct OvertoneRow {
    pub n: usize,
    pub freq_hz: f32,
    /// Cents offset from the nearest harmonic n * f0 (None if f0 unknown/invalid).
    pub cents_from_harm: Option<f32>,
    /// Median (spectral) level in dB (as provided by Step 4/6 summary).
    pub amplitude_db: f32,
    /// Linear amplitude proxy derived from amplitude_db (20*log10).
    pub lin_amp: f32,
    /// Perceptual weight from Step 6.
    pub weight: f32,
}

#[derive(Debug, Clone)]
pub struct NamedInterval {
    /// e.g., "3:2 (perfect fifth)" or similar label from Step 8
    pub label: String,
    /// Best-fit ratio for the minimum (e.g., 1.5 for 3:2)
    pub ratio: f32,
    /// The D(r) minimum location in cents
    pub at_cents: f32,
    /// |c(min) - c(p/q)| in cents (signed error where available from Step 8)
    pub cents_error: f32,
    /// Confidence from Step 8
    pub confidence: f32,
}

#[inline]
fn ratio_to_cents(r: f32) -> f32 {
    1200.0 * (r as f64).log2() as f32
}

/// Build an overtone table from perceptually-weighted partials.
/// `weighted` is (freq_hz, amplitude_db, weight).
pub fn overtone_table(
    f0_hz: Option<f32>,
    weighted: &[(f32, f32, f32)],
    max_rows: usize,
) -> Vec<OvertoneRow> {
    // Build raw rows with lin_amp precomputed.
    let mut rows: Vec<OvertoneRow> = weighted
        .iter()
        .map(|&(f_hz, amp_db, w)| OvertoneRow {
            n: 0,
            freq_hz: f_hz,
            cents_from_harm: None,
            amplitude_db: amp_db,
            lin_amp: 10f32.powf(amp_db / 20.0),
            weight: w,
        })
        .collect();

    // If we have a valid f0, compute harmonic index n and Δcents.
    let f0_opt = f0_hz.filter(|f| *f > 0.0);
    if let Some(f0) = f0_opt {
        for row in rows.iter_mut() {
            let n = (row.freq_hz / f0).round().max(1.0) as usize;
            row.n = n;
            let ideal = n as f32 * f0;
            let cents = if ideal > 0.0 {
                1200.0 * ((row.freq_hz / ideal) as f64).log2() as f32
            } else {
                0.0
            };
            row.cents_from_harm = Some(cents);
        }

        // Keep the strongest (by salience) per harmonic n.
        let mut best_by_n = HashMap::<usize, OvertoneRow>::new();
        for row in rows.into_iter() {
            let salience = row.lin_amp * row.weight;
            best_by_n
                .entry(row.n)
                .and_modify(|best| {
                    let best_sal = best.lin_amp * best.weight;
                    let better_sal = salience > best_sal + 1e-9;
                    let tie_sal = (salience - best_sal).abs() <= 1e-9;
                    let row_dc = row.cents_from_harm.unwrap_or(f32::INFINITY).abs();
                    let best_dc = best.cents_from_harm.unwrap_or(f32::INFINITY).abs();
                    let better_dc = row_dc < best_dc - 1e-6;
                    let tie_dc = (row_dc - best_dc).abs() <= 1e-6;
                    let better_freq = row.freq_hz < best.freq_hz;
                    if better_sal || (tie_sal && (better_dc || (tie_dc && better_freq))) {
                        *best = row;
                    }
                })
                .or_insert(row);
        }

        // Collect winners and sort: salience desc, then frequency asc.
        let mut dedup: Vec<OvertoneRow> = best_by_n.into_values().collect();
        dedup.sort_by(|a, b| {
            let sa = a.lin_amp * a.weight;
            let sb = b.lin_amp * b.weight;
            sb.partial_cmp(&sa)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    a.freq_hz
                        .partial_cmp(&b.freq_hz)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        });
        if dedup.len() > max_rows {
            dedup.truncate(max_rows);
        }
        dedup
    } else {
        // No f0: rank purely by salience, then frequency.
        rows.sort_by(|a, b| {
            let sa = a.lin_amp * a.weight;
            let sb = b.lin_amp * b.weight;
            sb.partial_cmp(&sa)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    a.freq_hz
                        .partial_cmp(&b.freq_hz)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        });
        if rows.len() > max_rows {
            rows.truncate(max_rows);
        }
        rows
    }
}

/// Create a succinct two‑note tuning suggestion from Step 8 named minima.
/// Picks the candidate with the highest confidence (breaking ties by smallest |cents_error|).
pub fn tuning_suggestion_two_note(named: &[NamedInterval]) -> Option<String> {
    if named.is_empty() {
        return None;
    }
    let mut best = named[0].clone();
    for c in named.iter().skip(1) {
        let better_conf = c.confidence > best.confidence + 1e-6;
        let same_conf_better_err = (c.confidence - best.confidence).abs() <= 1e-6
            && c.cents_error.abs() < best.cents_error.abs();
        if better_conf || same_conf_better_err {
            best = c.clone();
        }
    }
    // Signed offset of D-minimum from the exact rational position.
    let rational_cents = ratio_to_cents(best.ratio);
    let offset = best.at_cents - rational_cents;
    let dir = if offset < 0.0 { "flat of" } else { "sharp of" };
    let msg = format!(
        "Best interval: {} (~{:.5}). D-minimum at {:.1} cents → about {:.1} cents {} {}.\n\
         Suggestion: tune the upper note ≈ {:.1} cents {} the {} minimum to align prominent partials.",
        best.label,
        best.ratio,
        best.at_cents,
        offset.abs(),
        dir,
        best.label,
        offset.abs(),
        dir,
        best.label
    );
    Some(msg)
}

type MultiNote<'a> = Option<&'a [(f32, f32, Option<String>, Option<f32>)]>;

/// Build a compact textual summary for the N‑note case (Step 9).
/// `multi` is per-note: (ratio r_k, cents, optional label, optional |cents_error|).
pub fn tuning_suggestion_multi_note(multi: MultiNote) -> Option<String> {
    let list = match multi {
        Some(l) if !l.is_empty() => l,
        _ => return None,
    };
    // Note 0 is the reference by construction.
    let mut parts: Vec<String> = Vec::new();
    for (idx, (ratio, cents, label_opt, _err_opt)) in list.iter().enumerate() {
        if idx == 0 {
            parts.push("Note 0: r=1.00000 (0.0 c) [reference]".to_string());
        } else if let Some(label) = label_opt {
            parts.push(format!(
                "Note {}: ~ {} ({:+.1} c, r≈{:.5})",
                idx, label, cents, ratio
            ));
        } else {
            parts.push(format!("Note {}: ({:+.1} c, r≈{:.5})", idx, cents, ratio));
        }
    }
    Some(format!(
        "Suggested set relative to reference → {}.",
        parts.join(" · ")
    ))
}

pub fn export_dissonance_csv(
    path: impl AsRef<Path>,
    samples: &[(f32, f32)],
) -> std::io::Result<()> {
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    writeln!(w, "cents,D")?;
    for (c, d) in samples {
        writeln!(w, "{:.6},{:.9}", c, d)?;
    }
    Ok(())
}

pub fn export_overtone_csv(path: impl AsRef<Path>, rows: &[OvertoneRow]) -> std::io::Result<()> {
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    writeln!(w, "n,freq_hz,cents_from_harm,amplitude_db,lin_amp,weight")?;
    for r in rows {
        let cents_str = r
            .cents_from_harm
            .map(|x| format!("{:.3}", x))
            .unwrap_or_else(|| "".into());
        writeln!(
            w,
            "{},{:.6},{},{:.3},{:.6},{:.6}",
            r.n, r.freq_hz, cents_str, r.amplitude_db, r.lin_amp, r.weight
        )?;
    }
    Ok(())
}
