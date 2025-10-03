pub(super) fn summary_min_median_max(xs: &[f32]) -> (f32, f32, f32) {
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

pub(super) fn mean(xs: &[f32]) -> f32 {
    if xs.is_empty() {
        0.0
    } else {
        xs.iter().sum::<f32>() / xs.len() as f32
    }
}

pub(super) fn nearest_shift(baseline: &[f32], variant: &[f32]) -> (f32, Vec<f32>) {
    // For each baseline min, find nearest variant min and report |Δcents|.
    let mut shifts = Vec::with_capacity(baseline.len());
    for &c0 in baseline {
        let d = variant
            .iter()
            .map(|&c1| (c1 - c0).abs())
            .fold(f32::INFINITY, f32::min);
        shifts.push(d);
    }
    (mean(&shifts), shifts)
}

pub(super) fn write_csv(path: &str, rows: &[String]) -> std::io::Result<()> {
    use std::io::Write;
    let f = std::fs::File::create(path)?;
    let mut w = std::io::BufWriter::new(f);
    for line in rows {
        writeln!(w, "{line}")?;
    }
    Ok(())
}
