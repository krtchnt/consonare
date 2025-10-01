/// A dense 2D grid of total roughness for N=3 (two free ratios vs a fixed reference).
#[cfg(feature = "visualise")]
#[derive(Clone, Debug)]
pub struct DissonanceGrid2D {
    pub x_cents: Vec<f32>,     // r1 cents grid (relative to reference)
    pub y_cents: Vec<f32>,     // r2 cents grid (relative to reference)
    pub values: Vec<Vec<f32>>, // values[j][i] is D at (x_cents[i], y_cents[j])
    pub vmin: f32,
    pub vmax: f32,
}

/// Compute a 2D roughness grid for N=3 by sweeping r1 and r2 cents.
/// If `others` is Some(&[b1, b2]), their spectra are used; else both notes reuse `a_ref`’s timbre.
#[cfg(feature = "visualise")]
pub fn compute_dissonance_grid_3(
    a_ref: &consonare_core::weighting::WeightingResult,
    others: Option<&[consonare_core::weighting::WeightingResult]>,
    min_cents: f32,
    max_cents: f32,
    step_cents: f32,
    max_deltaf_over_cbw: Option<f32>,
) -> DissonanceGrid2D {
    let mut x_cents = Vec::new();
    let mut c = min_cents;
    while c <= max_cents + 1e-9 {
        x_cents.push(c);
        c += step_cents;
    }
    let y_cents = x_cents.clone();

    let w = x_cents.len();
    let h = y_cents.len();
    let mut values = vec![vec![0.0f32; w]; h];
    let mut vmin = f32::INFINITY;
    let mut vmax = f32::NEG_INFINITY;

    for (j, &cy) in y_cents.iter().enumerate() {
        for (i, &cx) in x_cents.iter().enumerate() {
            // Reuse the pairwise kernel & pruning from Step 7 via our multi-note evaluator.
            let v = consonare_core::dissonance_n::total_roughness_for_cents(
                a_ref,
                others,
                &[cx, cy],
                max_deltaf_over_cbw,
            );
            values[j][i] = v;
            if v < vmin {
                vmin = v;
            }
            if v > vmax {
                vmax = v;
            }
        }
    }

    DissonanceGrid2D {
        x_cents,
        y_cents,
        values,
        vmin,
        vmax,
    }
}

#[cfg(feature = "visualise")]
pub fn plot_dissonance_surface_3(
    grid: &DissonanceGrid2D,
    title: impl AsRef<str>,
    fname: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    use plotters::prelude::*;

    use crate::common::gen_target_fn;

    let out_file = gen_target_fn(format!("{}_dissonance-3.png", fname));

    // Canvas
    let img_w = 1200;
    let img_h = 1000;
    let root = BitMapBackend::new(&out_file, (img_w, img_h)).into_drawing_area();
    root.fill(&WHITE)?;

    let (xmin, xmax) = (
        *grid.x_cents.first().unwrap_or(&0.0),
        *grid.x_cents.last().unwrap_or(&1.0),
    );
    let (ymin, ymax) = (
        *grid.y_cents.first().unwrap_or(&0.0),
        *grid.y_cents.last().unwrap_or(&1.0),
    );

    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 24))
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(60)
        .build_cartesian_2d(xmin..xmax, ymin..ymax)?;

    chart
        .configure_mesh()
        .x_desc("Note 1 ratio (cents)")
        .y_desc("Note 2 ratio (cents)")
        .x_labels(((xmax - xmin) / 100.0).ceil() as usize + 1)
        .y_labels(((ymax - ymin) / 100.0).ceil() as usize + 1)
        .x_label_formatter(&|x| format!("{:.0}", x))
        .y_label_formatter(&|y| format!("{:.0}", y))
        .draw()?;

    // Simple heat colormap (blue -> yellow -> red)
    let map_color = |v: f32| -> RGBColor {
        let vmin = grid.vmin;
        let vmax = grid.vmax.max(vmin + 1e-9);
        let t = ((v - vmin) / (vmax - vmin)).clamp(0.0, 1.0);
        if t < 0.5 {
            let u = t / 0.5;
            RGBColor(
                (255.0 * u) as u8,
                (255.0 * u) as u8,
                (255.0 * (1.0 - u)) as u8,
            ) // blue->yellow
        } else {
            let u = (t - 0.5) / 0.5;
            RGBColor(255, (255.0 * (1.0 - u)) as u8, 0) // yellow->red
        }
    };

    let nx = grid.x_cents.len();
    let ny = grid.y_cents.len();

    // Draw each cell as a filled rectangle in data coords.
    chart.draw_series((0..ny.saturating_sub(1)).flat_map(|j| {
        (0..nx.saturating_sub(1)).map(move |i| {
            let x0 = grid.x_cents[i];
            let x1 = grid.x_cents[i + 1];
            let y0 = grid.y_cents[j];
            let y1 = grid.y_cents[j + 1];
            let v = grid.values[j][i];
            let color = map_color(v);
            Rectangle::new([(x0, y0), (x1, y1)], ShapeStyle::from(&color).filled())
        })
    }))?;

    Ok(())
}
