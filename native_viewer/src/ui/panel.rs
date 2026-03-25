//! Individual genome panel widget for rendering pileup data.

use egui::{Rect, Sense, Ui, Vec2};

use crate::filters::DisplayFilters;
use crate::genome::region::GenomicRegion;
use crate::pileup::PileupData;
use crate::ui::colors;

/// Render a single genome panel with header and pileup area.
pub fn render_panel(
    ui: &mut Ui,
    title: &str,
    region: &Option<GenomicRegion>,
    pileup: &Option<PileupData>,
    filters: &DisplayFilters,
) {
    // Panel header
    let header_height = 28.0;
    let (header_rect, _) = ui.allocate_exact_size(
        Vec2::new(ui.available_width(), header_height),
        Sense::hover(),
    );

    ui.painter()
        .rect_filled(header_rect, 0.0, colors::PANEL_HEADER_BG);

    // Draw header text
    let mut text_pos = header_rect.left_center();
    text_pos.x += 8.0;
    ui.painter().text(
        text_pos,
        egui::Align2::LEFT_CENTER,
        title,
        egui::FontId::proportional(14.0),
        egui::Color32::WHITE,
    );

    if let Some(region) = region {
        let region_text_pos = egui::pos2(text_pos.x + 120.0, text_pos.y);
        ui.painter().text(
            region_text_pos,
            egui::Align2::LEFT_CENTER,
            format!("{region}"),
            egui::FontId::proportional(11.0),
            egui::Color32::from_rgb(180, 180, 180),
        );
    }

    // Pileup rendering area
    let available = ui.available_rect_before_wrap();
    let pileup_rect =
        Rect::from_min_size(available.min, Vec2::new(available.width(), available.height()));

    // Background
    ui.painter()
        .rect_filled(pileup_rect, 0.0, colors::PANEL_BG);

    let (response, painter) = ui.allocate_painter(pileup_rect.size(), Sense::hover());
    let rect = response.rect;

    match (region, pileup) {
        (Some(region), Some(pileup)) => {
            render_pileup(&painter, rect, region, pileup, filters);
        }
        (Some(_region), None) => {
            // Region set but no data loaded yet
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "Loading...",
                egui::FontId::proportional(14.0),
                egui::Color32::GRAY,
            );
        }
        _ => {
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "No region selected",
                egui::FontId::proportional(14.0),
                egui::Color32::GRAY,
            );
        }
    }
}

/// Render pileup reads within a panel.
fn render_pileup(
    painter: &egui::Painter,
    rect: Rect,
    region: &GenomicRegion,
    pileup: &PileupData,
    filters: &DisplayFilters,
) {
    let region_span = region.span().max(1) as f32;
    let row_height = filters.row_height();
    let row_gap = filters.row_gap();

    // Coordinate axis at the top
    let axis_height = 20.0;
    render_axis(painter, rect, region, axis_height);

    let pileup_top = rect.min.y + axis_height;

    // Render each row of reads
    for (row_idx, row) in pileup.rows.iter().enumerate() {
        let y = pileup_top + row_idx as f32 * (row_height + row_gap);

        // Stop if we've exceeded the panel height
        if y + row_height > rect.max.y {
            break;
        }

        for read in &row.reads {
            // Calculate x positions from genomic coordinates
            let x_start = rect.min.x
                + (read.start.saturating_sub(region.start) as f32 / region_span) * rect.width();
            let x_end = rect.min.x
                + (read.end.saturating_sub(region.start) as f32 / region_span) * rect.width();

            // Clamp to panel bounds
            let x_start = x_start.max(rect.min.x);
            let x_end = x_end.min(rect.max.x);

            if x_end <= x_start {
                continue;
            }

            let read_color = colors::read_color(read.hp_tag, filters.squished);
            let read_rect = Rect::from_min_max(
                egui::pos2(x_start, y),
                egui::pos2(x_end, y + row_height),
            );
            painter.rect_filled(read_rect, 1.0, read_color);
        }
    }

    // Draw read count indicator
    let count_text = format!(
        "{} reads ({} rows)",
        pileup.total_reads,
        pileup.rows.len()
    );
    painter.text(
        egui::pos2(rect.max.x - 5.0, rect.min.y + axis_height + 2.0),
        egui::Align2::RIGHT_TOP,
        count_text,
        egui::FontId::proportional(10.0),
        egui::Color32::from_rgba_premultiplied(100, 100, 100, 180),
    );
}

/// Render the coordinate axis at the top of a panel.
fn render_axis(painter: &egui::Painter, rect: Rect, region: &GenomicRegion, height: f32) {
    let axis_rect = Rect::from_min_size(rect.min, Vec2::new(rect.width(), height));

    // Axis background
    painter.rect_filled(axis_rect, 0.0, egui::Color32::from_rgb(245, 245, 245));

    // Axis line
    painter.line_segment(
        [
            egui::pos2(rect.min.x, rect.min.y + height - 1.0),
            egui::pos2(rect.max.x, rect.min.y + height - 1.0),
        ],
        egui::Stroke::new(1.0, egui::Color32::from_rgb(180, 180, 180)),
    );

    // Tick marks and labels
    let span = region.span();
    if span == 0 {
        return;
    }

    // Calculate nice tick interval
    let approx_ticks = 5;
    let raw_interval = span as f64 / approx_ticks as f64;
    let magnitude = 10f64.powf(raw_interval.log10().floor());
    let normalized = raw_interval / magnitude;
    let nice_interval = if normalized <= 1.5 {
        1.0
    } else if normalized <= 3.5 {
        2.0
    } else if normalized <= 7.5 {
        5.0
    } else {
        10.0
    } * magnitude;

    let tick_interval = nice_interval.max(1.0) as u64;
    let first_tick = (region.start / tick_interval) * tick_interval;

    let mut tick = first_tick;
    while tick <= region.end {
        if tick >= region.start {
            let x_frac = (tick - region.start) as f32 / span as f32;
            let x = rect.min.x + x_frac * rect.width();

            // Tick line
            painter.line_segment(
                [
                    egui::pos2(x, rect.min.y + height - 5.0),
                    egui::pos2(x, rect.min.y + height - 1.0),
                ],
                egui::Stroke::new(1.0, egui::Color32::from_rgb(150, 150, 150)),
            );

            // Label
            let label = format_position(tick);
            painter.text(
                egui::pos2(x, rect.min.y + height - 7.0),
                egui::Align2::CENTER_BOTTOM,
                label,
                egui::FontId::proportional(9.0),
                egui::Color32::from_rgb(80, 80, 80),
            );
        }
        tick += tick_interval;
    }
}

/// Format a genomic position for display (e.g. "1.5 Mb", "250 kb").
fn format_position(pos: u64) -> String {
    if pos >= 1_000_000 {
        format!("{:.1} Mb", pos as f64 / 1_000_000.0)
    } else if pos >= 1000 {
        format!("{:.1} kb", pos as f64 / 1000.0)
    } else {
        format!("{pos}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_position() {
        assert_eq!(format_position(500), "500");
        assert_eq!(format_position(1500), "1.5 kb");
        assert_eq!(format_position(1_500_000), "1.5 Mb");
    }
}
