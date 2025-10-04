#[cfg(feature = "visualise")]
use std::fmt::Display;

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

#[cfg(feature = "visualise")]
pub fn compute_dissonance_grid_3(
    a_ref: &consonare_core::weighting::WeightingResult,
    others: Option<&[consonare_core::weighting::WeightingResult]>,
    min_cents: f32,
    max_cents: f32,
    step_cents: f32,
    max_deltaf_over_cbw: Option<f32>,
) -> DissonanceGrid2D {
    use rayon::prelude::*;

    // 1) Axes with integer-derived cents to avoid drift
    let nx = ((max_cents - min_cents) / step_cents).round() as usize + 1;
    let x_cents: Vec<f32> = (0..nx).map(|k| min_cents + k as f32 * step_cents).collect();
    let y_cents = x_cents.clone();

    let delta_min = min_cents - max_cents;
    let delta_max = max_cents - min_cents;
    let ndelta = ((delta_max - delta_min) / step_cents).round() as usize + 1;
    let delta_cents: Vec<f32> = (0..ndelta)
        .map(|k| delta_min + k as f32 * step_cents)
        .collect();

    // 2) Resolve timbres per pair
    let (b1, b2) = match others {
        Some(arr) if arr.len() >= 2 => (&arr[0], &arr[1]),
        _ => (a_ref, a_ref),
    };

    // Helper: pair roughness at one interval (prefer a dedicated pair fn if available)
    let pair = |left: &_, right: &_, cents: f32| -> f32 {
        // If you add a dedicated pair API, call it here instead.
        consonare_core::dissonance_n::total_roughness_for_cents(
            left,
            Some(std::slice::from_ref(right)),
            std::slice::from_ref(&cents),
            max_deltaf_over_cbw,
        )
    };

    // 3) Precompute the three 1-D curves in parallel
    let (lx, (ly, ldelta)) = rayon::join(
        || {
            x_cents
                .par_iter()
                .map(|&cx| pair(a_ref, b1, cx))
                .collect::<Vec<_>>()
        },
        || {
            rayon::join(
                || {
                    y_cents
                        .par_iter()
                        .map(|&cy| pair(a_ref, b2, cy))
                        .collect::<Vec<_>>()
                },
                || {
                    delta_cents
                        .par_iter()
                        .map(|&d| pair(b1, b2, d))
                        .collect::<Vec<_>>()
                },
            )
        },
    );

    // 4) Assemble the grid: v[i,j] = Lx[i] + Ly[j] + LΔ[j-i]
    let w = x_cents.len();
    let h = y_cents.len();
    let mut values = vec![vec![0.0f32; w]; h];

    let mut vmin = f32::INFINITY;
    let mut vmax = f32::NEG_INFINITY;

    values
        .par_iter_mut()
        .enumerate()
        .map(|(j, row)| {
            let cy = y_cents[j];
            let base_y = ly[j];
            let mut local_min = f32::INFINITY;
            let mut local_max = f32::NEG_INFINITY;

            for i in 0..w {
                let cx = x_cents[i];
                // Index Δ robustly (all points are multiples of step_cents)
                let d = cy - cx;
                let k = (((d - delta_min) / step_cents).round() as isize)
                    .clamp(0, (ndelta - 1) as isize) as usize;

                let v = lx[i] + base_y + ldelta[k];
                row[i] = v;
                if v < local_min {
                    local_min = v;
                }
                if v > local_max {
                    local_max = v;
                }
            }
            (local_min, local_max)
        })
        .collect::<Vec<_>>()
        .into_iter()
        .for_each(|(lo, hi)| {
            if lo < vmin {
                vmin = lo;
            }
            if hi > vmax {
                vmax = hi;
            }
        });

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
    fname: &str,
    out_dir: impl Display,
) -> Result<(), Box<dyn std::error::Error>> {
    use image::{DynamicImage, imageops::FilterType};
    use plotters::{backend::RGBPixel, prelude::*};

    let out_file = format!("{}/{}_dissonance-3.png", out_dir, fname);
    let (img_w, img_h) = (1200, 1000);
    let root = BitMapBackend::new(&out_file, (img_w, img_h)).into_drawing_area();
    root.fill(&WHITE)?;

    let xmin = *grid.x_cents.first().unwrap_or(&0.0);
    let xmax = *grid.x_cents.last().unwrap_or(&1.0);
    let ymin = *grid.y_cents.first().unwrap_or(&0.0);
    let ymax = *grid.y_cents.last().unwrap_or(&1.0);

    let caption = format!("Dissonance Surface of `{}` (N=3)", fname);
    let mut chart = ChartBuilder::on(&root)
        .caption(caption, ("sans-serif", 24))
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

    // Build a pixel-sized image (one pixel per grid cell).
    let w = grid.x_cents.len().saturating_sub(1);
    let h = grid.y_cents.len().saturating_sub(1);
    let mut raw: Vec<u8> = vec![0; w * h * 3];

    let vmin = grid.vmin;
    let vmax = (grid.vmax).max(vmin + 1e-9);
    let inv = 1.0 / (vmax - vmin);

    // Precompute a 1025-entry color LUT
    let lut: Vec<[u8; 3]> = (0..=1024)
        .map(|i| {
            let t = i as f64 / 1024.0;
            let c = colorous::CUBEHELIX.eval_continuous(t);
            [c.r, c.g, c.b]
        })
        .collect();

    for j in 0..h {
        for i in 0..w {
            let v = grid.values[j][i];
            let t = ((v - vmin) * inv).clamp(0.0, 1.0);
            let idx = (t * 1024.0) as usize;
            // Flip vertically so chart axes align
            let jj = h - 1 - j;
            let offset = (jj * w + i) * 3;
            raw[offset..offset + 3].copy_from_slice(&lut[idx]);
        }
    }

    // Scale the raw grid to the plotting area
    let (pw, ph) = chart.plotting_area().dim_in_pixel();
    let dyn_scaled =
        DynamicImage::ImageRgb8(image::RgbImage::from_raw(w as u32, h as u32, raw).unwrap())
            .resize_exact(pw, ph, FilterType::Nearest);

    let rgb = dyn_scaled.to_rgb8();
    let (w, h) = rgb.dimensions();
    let buf = rgb.into_raw();

    let bitmap = BitMapElement::<_, RGBPixel>::with_owned_buffer((xmin, ymax), (w, h), buf)
        .expect("buffer size mismatch");

    chart.draw_series(std::iter::once(bitmap))?;
    root.present()?;
    Ok(())
}
