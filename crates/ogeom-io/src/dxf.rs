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
