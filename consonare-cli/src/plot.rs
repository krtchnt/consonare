#[cfg(feature = "visualise")]
pub fn plot_dissonance(
    d_result: &consonare_core::dissonance::DissonanceResult,
    fname: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    use plotters::prelude::*;

    use crate::common::gen_target_fn;

    let out_file = gen_target_fn(format!("{}_dissonance.png", fname));

    let root = BitMapBackend::new(&out_file, (1200, 900)).into_drawing_area();
    root.fill(&WHITE)?;

    let (min_cents, max_cents) = (
        d_result.smoothed.first().map(|p| p.cents).unwrap_or(0.0),
        d_result.smoothed.last().map(|p| p.cents).unwrap_or(1.0),
    );
    let (min_val, _, max_val) = super::summary_min_median_max(
        &d_result
            .smoothed
            .iter()
            .map(|p| p.value)
            .collect::<Vec<_>>(),
    );

    let mut chart = ChartBuilder::on(&root)
        .caption(
            format!("Dissonance Profile of `{}`", fname),
            ("sans-serif", 24),
        )
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(min_cents..max_cents, min_val..max_val)?;

    chart
        .configure_mesh()
        .x_desc("Interval (cents)")
        .y_desc("Dissonance")
        .x_labels(((max_cents - min_cents) / 100.0).ceil() as usize + 1)
        .x_label_formatter(&|x| format!("{:.0}", x))
        .draw()?;

    chart
        .draw_series(LineSeries::new(
            d_result.smoothed.iter().map(|p| (p.cents, p.value)),
            &RED,
        ))?
        .label("Smoothed")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], RED));

    chart
        .draw_series(
            d_result
                .minima
                .iter()
                .map(|m| Circle::new((m.cents, m.value), 5, ShapeStyle::from(&BLUE).filled())),
        )?
        .label("Minima")
        .legend(|(x, y)| Circle::new((x + 10, y), 5, ShapeStyle::from(&BLUE).filled()));

    chart.configure_series_labels().border_style(BLACK).draw()?;

    Ok(())
}
