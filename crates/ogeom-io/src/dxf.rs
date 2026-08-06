//! Writing DXF: 2D polylines, the drawing interchange the field expects.
//!
//! R12 ASCII, the most widely readable dialect: a TABLES section declaring
//! the two linetypes and two layers a technical drawing needs, then one
//! `POLYLINE` per curve. Visible curves go on the `VISIBLE` layer with
//! continuous lines; hidden curves on `HIDDEN`, dashed. The writer takes
//! bare polylines rather than a drawing type, so anything that produces 2D
//! curves — the hidden-line projector, a section outline, a sketch — writes
//! without this crate knowing where they came from.

use ogeom_math::Point2;
use std::fmt::Write as _;

/// Write polylines as an R12 DXF document.
///
/// `visible` draws continuous on layer `VISIBLE`; `hidden` draws dashed on
/// layer `HIDDEN`. Polylines with fewer than two points are skipped — a
/// point is not a line in a drawing, here as everywhere.
#[must_use]
pub fn write_dxf(visible: &[Vec<Point2>], hidden: &[Vec<Point2>]) -> String {
    let mut out = String::new();
    // Header: R12 says almost nothing and needs almost nothing.
    push(
        &mut out,
        &[("0", "SECTION"), ("2", "HEADER"), ("0", "ENDSEC")],
    );

    // Tables: the linetypes first, because the layers name them.
    push(&mut out, &[("0", "SECTION"), ("2", "TABLES")]);
    push(&mut out, &[("0", "TABLE"), ("2", "LTYPE"), ("70", "2")]);
    push(
        &mut out,
        &[
            ("0", "LTYPE"),
            ("2", "CONTINUOUS"),
            ("70", "0"),
            ("3", "Solid line"),
            ("72", "65"),
            ("73", "0"),
            ("40", "0.0"),
        ],
    );
    push(
        &mut out,
        &[
            ("0", "LTYPE"),
            ("2", "DASHED"),
            ("70", "0"),
            ("3", "Dashed line"),
            ("72", "65"),
            ("73", "2"),
            ("40", "0.75"),
            ("49", "0.5"),
            ("49", "-0.25"),
        ],
    );
    push(&mut out, &[("0", "ENDTAB")]);
    push(&mut out, &[("0", "TABLE"), ("2", "LAYER"), ("70", "2")]);
    push(
        &mut out,
        &[
            ("0", "LAYER"),
            ("2", "VISIBLE"),
            ("70", "0"),
            ("62", "7"),
            ("6", "CONTINUOUS"),
        ],
    );
    push(
        &mut out,
        &[
            ("0", "LAYER"),
            ("2", "HIDDEN"),
            ("70", "0"),
            ("62", "8"),
            ("6", "DASHED"),
        ],
    );
    push(&mut out, &[("0", "ENDTAB"), ("0", "ENDSEC")]);

    push(&mut out, &[("0", "SECTION"), ("2", "ENTITIES")]);
    for (layer, curves) in [("VISIBLE", visible), ("HIDDEN", hidden)] {
        for curve in curves {
            if curve.len() < 2 {
                continue;
            }
            push(
                &mut out,
                &[("0", "POLYLINE"), ("8", layer), ("66", "1"), ("70", "0")],
            );
            for p in curve {
                push(&mut out, &[("0", "VERTEX"), ("8", layer)]);
                let _ = writeln!(out, "10\n{}\n20\n{}\n30\n0.0", real(p.x), real(p.y));
            }
            push(&mut out, &[("0", "SEQEND")]);
        }
    }
    push(&mut out, &[("0", "ENDSEC"), ("0", "EOF")]);
    out
}

/// A DXF real: shortest exact form, decimal point guaranteed.
fn real(v: f64) -> String {
    let s = format!("{v:?}");
    if s.contains('.') || s.contains('e') {
        s
    } else {
        format!("{s}.0")
    }
}

/// Append group-code/value pairs, one per line each.
fn push(out: &mut String, pairs: &[(&str, &str)]) {
    for (code, value) in pairs {
        let _ = writeln!(out, "{code}\n{value}");
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn a_drawing_writes_layers_polylines_and_exact_coordinates() {
        let visible = vec![vec![
            Point2::new(0.0, 0.0),
            Point2::new(10.0, 0.0),
            Point2::new(10.0, 5.0),
        ]];
        let hidden = vec![vec![Point2::new(1.5, 2.25), Point2::new(3.0, 2.25)]];
        let text = write_dxf(&visible, &hidden);

        assert!(text.starts_with("0\nSECTION"));
        assert!(text.trim_end().ends_with("EOF"));
        assert_eq!(text.matches("POLYLINE").count(), 2);
        assert_eq!(text.matches("VERTEX").count(), 5);
        assert_eq!(text.matches("SEQEND").count(), 2);
        // Both layers declared and used, hidden dashed.
        assert!(text.contains("VISIBLE"));
        assert!(text.contains("HIDDEN"));
        assert!(text.contains("DASHED"));
        // Coordinates exact and point-carrying.
        assert!(text.contains("10\n10.0\n20\n5.0"));
        assert!(text.contains("10\n1.5\n20\n2.25"));
    }

    #[test]
    fn degenerate_polylines_are_dropped() {
        let text = write_dxf(&[vec![Point2::new(1.0, 1.0)]], &[vec![]]);
        assert_eq!(text.matches("POLYLINE").count(), 0);
    }
}

/// The polylines a DXF carries, by the layer they were drawn on.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct DxfDrawing {
    /// Curves on the `VISIBLE` layer, or on no named layer at all.
    pub visible: Vec<Vec<Point2>>,
    /// Curves on the `HIDDEN` layer.
    pub hidden: Vec<Vec<Point2>>,
}

/// Read the polylines out of an ASCII DXF.
///
/// `POLYLINE`/`VERTEX`/`SEQEND` and `LWPOLYLINE` become polylines; `LINE`
/// becomes a polyline of two points. Layers are read by name: `HIDDEN`
/// separates, everything else is visible — which is the convention this
/// crate's own writer uses and the one a drawing's reader can act on
/// without guessing at linetypes.
///
/// Everything else in a DXF — blocks, text, dimensions, splines, hatches —
/// is skipped. This reads *drawings as curves*, which is what a kernel has
/// use for; it does not pretend to be a DXF application.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// file is not group-coded in pairs, which is the one thing every DXF is.
pub fn read_dxf(text: &str) -> ogeom_core::OgeomResult<DxfDrawing> {
    // A DXF is a stream of (code, value) pairs, one per line each.
    let lines: Vec<&str> = text.lines().map(str::trim).collect();
    if !lines.len().is_multiple_of(2) {
        // A trailing blank line is ordinary; anything else is not pairs.
        if !lines.last().is_some_and(|l| l.is_empty()) {
            ogeom_core::ogeom_bail!(
                Construction,
                "a DXF is group codes and values in pairs; this has an odd number of lines"
            );
        }
    }

    let mut out = DxfDrawing::default();
    let mut entity = String::new();
    let mut layer = String::new();
    let mut points: Vec<Point2> = Vec::new();
    let mut pending: Option<(Option<f64>, Option<f64>)> = None;
    let mut line_ends: [Option<Point2>; 2] = [None, None];

    let flush = |entity: &str, layer: &str, points: &mut Vec<Point2>, out: &mut DxfDrawing| {
        if points.len() >= 2 {
            let drawn = core::mem::take(points);
            if layer.eq_ignore_ascii_case("HIDDEN") {
                out.hidden.push(drawn);
            } else {
                out.visible.push(drawn);
            }
        } else {
            points.clear();
        }
        let _ = entity;
    };

    let mut i = 0;
    while i + 1 < lines.len() {
        let code = lines[i];
        let value = lines[i + 1];
        i += 2;
        match code {
            "0" => {
                // A new entity ends whatever was being gathered.
                match entity.as_str() {
                    "POLYLINE" | "LWPOLYLINE" => {
                        if let Some((Some(x), Some(y))) = pending.take() {
                            points.push(Point2::new(x, y));
                        }
                        if value != "VERTEX" && value != "SEQEND" {
                            flush(&entity, &layer, &mut points, &mut out);
                        }
                    }
                    "LINE" => {
                        if let [Some(a), Some(b)] = line_ends {
                            let mut segment = vec![a, b];
                            if layer.eq_ignore_ascii_case("HIDDEN") {
                                out.hidden.push(core::mem::take(&mut segment));
                            } else {
                                out.visible.push(core::mem::take(&mut segment));
                            }
                        }
                        line_ends = [None, None];
                    }
                    _ => {}
                }
                if value == "VERTEX" {
                    if let Some((Some(x), Some(y))) = pending.take() {
                        points.push(Point2::new(x, y));
                    }
                    pending = Some((None, None));
                    continue;
                }
                if value == "SEQEND" {
                    flush("POLYLINE", &layer, &mut points, &mut out);
                    entity.clear();
                    continue;
                }
                entity = value.to_string();
                if entity == "POLYLINE" || entity == "LWPOLYLINE" {
                    points.clear();
                    pending = Some((None, None));
                } else if entity == "LINE" {
                    line_ends = [None, None];
                }
                layer.clear();
            }
            "8" => layer = value.to_string(),
            "10" | "20" | "11" | "21" => {
                let Ok(number) = value.parse::<f64>() else {
                    continue;
                };
                match (entity.as_str(), code) {
                    ("LINE", "10") => {
                        line_ends[0] = Some(Point2::new(number, line_ends[0].map_or(0.0, |p| p.y)));
                    }
                    ("LINE", "20") => {
                        let x = line_ends[0].map_or(0.0, |p| p.x);
                        line_ends[0] = Some(Point2::new(x, number));
                    }
                    ("LINE", "11") => {
                        line_ends[1] = Some(Point2::new(number, line_ends[1].map_or(0.0, |p| p.y)));
                    }
                    ("LINE", "21") => {
                        let x = line_ends[1].map_or(0.0, |p| p.x);
                        line_ends[1] = Some(Point2::new(x, number));
                    }
                    (_, "10") => {
                        // An LWPOLYLINE gives its vertices as repeated
                        // 10/20 pairs without a VERTEX entity between them,
                        // so a fresh 10 closes the one before it.
                        if let Some((Some(x), Some(y))) = pending {
                            points.push(Point2::new(x, y));
                        }
                        pending = Some((Some(number), None));
                    }
                    (_, "20") => {
                        if let Some((x, _)) = pending {
                            pending = Some((x, Some(number)));
                        }
                    }
                    _ => {}
                }
            }
            _ => {}
        }
    }
    // Whatever the last entity was gathering.
    if let Some((Some(x), Some(y))) = pending.take() {
        points.push(Point2::new(x, y));
    }
    flush(&entity, &layer, &mut points, &mut out);
    if entity == "LINE"
        && let [Some(a), Some(b)] = line_ends
    {
        if layer.eq_ignore_ascii_case("HIDDEN") {
            out.hidden.push(vec![a, b]);
        } else {
            out.visible.push(vec![a, b]);
        }
    }
    Ok(out)
}
