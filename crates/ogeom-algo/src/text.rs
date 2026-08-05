//! Text to wires: a built-in single-stroke engraving font.
//!
//! The glyphs are the kernel's own, defined here as polyline strokes on a
//! six-by-nine grid — the shapes a CNC engraver or a drawing title block
//! wants, dependency-free and deterministic. Outline fonts are a file
//! format problem: when a TTF reader arrives with the exchange work, its
//! glyph outlines feed the same wire builder. Until then, what this speaks
//! it speaks exactly, and any character it does not carry is refused by
//! name rather than dropped.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::LineCurve;
use ogeom_math::{Frame, Point};
use ogeom_topo::{Model, Shape};

use crate::build::{make_compound, make_edge_between, make_vertex, make_wire};
use crate::history::{Built, History};

/// One glyph: strokes on a 6-wide grid with cap height 9, baseline 0.
type Glyph = &'static [&'static [(f64, f64)]];

/// The advance from one character origin to the next, in grid units.
const ADVANCE: f64 = 8.0;
/// The cap height the requested text height maps onto.
const CAP: f64 = 9.0;

fn glyph_of(c: char) -> Option<Glyph> {
    let glyph: Glyph = match c.to_ascii_uppercase() {
        'A' => &[
            &[(0.0, 0.0), (3.0, 9.0), (6.0, 0.0)],
            &[(1.0, 3.0), (5.0, 3.0)],
        ],
        'B' => &[
            &[
                (0.0, 0.0),
                (0.0, 9.0),
                (4.0, 9.0),
                (5.0, 8.0),
                (5.0, 6.0),
                (4.0, 5.0),
                (0.0, 5.0),
            ],
            &[(4.0, 5.0), (5.0, 4.0), (5.0, 1.0), (4.0, 0.0), (0.0, 0.0)],
        ],
        'C' => &[&[
            (6.0, 7.0),
            (4.0, 9.0),
            (2.0, 9.0),
            (0.0, 7.0),
            (0.0, 2.0),
            (2.0, 0.0),
            (4.0, 0.0),
            (6.0, 2.0),
        ]],
        'D' => &[&[
            (0.0, 0.0),
            (0.0, 9.0),
            (3.0, 9.0),
            (6.0, 7.0),
            (6.0, 2.0),
            (3.0, 0.0),
            (0.0, 0.0),
        ]],
        'E' => &[
            &[(6.0, 0.0), (0.0, 0.0), (0.0, 9.0), (6.0, 9.0)],
            &[(0.0, 5.0), (4.0, 5.0)],
        ],
        'F' => &[
            &[(0.0, 0.0), (0.0, 9.0), (6.0, 9.0)],
            &[(0.0, 5.0), (4.0, 5.0)],
        ],
        'G' => &[&[
            (6.0, 7.0),
            (4.0, 9.0),
            (2.0, 9.0),
            (0.0, 7.0),
            (0.0, 2.0),
            (2.0, 0.0),
            (4.0, 0.0),
            (6.0, 2.0),
            (6.0, 4.0),
            (3.0, 4.0),
        ]],
        'H' => &[
            &[(0.0, 0.0), (0.0, 9.0)],
            &[(6.0, 0.0), (6.0, 9.0)],
            &[(0.0, 5.0), (6.0, 5.0)],
        ],
        'I' => &[
            &[(3.0, 0.0), (3.0, 9.0)],
            &[(1.0, 9.0), (5.0, 9.0)],
            &[(1.0, 0.0), (5.0, 0.0)],
        ],
        'J' => &[&[(5.0, 9.0), (5.0, 2.0), (3.0, 0.0), (1.0, 0.0), (0.0, 2.0)]],
        'K' => &[
            &[(0.0, 0.0), (0.0, 9.0)],
            &[(6.0, 9.0), (0.0, 4.0)],
            &[(2.0, 5.5), (6.0, 0.0)],
        ],
        'L' => &[&[(0.0, 9.0), (0.0, 0.0), (6.0, 0.0)]],
        'M' => &[&[(0.0, 0.0), (0.0, 9.0), (3.0, 4.0), (6.0, 9.0), (6.0, 0.0)]],
        'N' => &[&[(0.0, 0.0), (0.0, 9.0), (6.0, 0.0), (6.0, 9.0)]],
        'O' => &[&[
            (2.0, 0.0),
            (0.0, 2.0),
            (0.0, 7.0),
            (2.0, 9.0),
            (4.0, 9.0),
            (6.0, 7.0),
            (6.0, 2.0),
            (4.0, 0.0),
            (2.0, 0.0),
        ]],
        'P' => &[&[
            (0.0, 0.0),
            (0.0, 9.0),
            (5.0, 9.0),
            (6.0, 8.0),
            (6.0, 6.0),
            (5.0, 5.0),
            (0.0, 5.0),
        ]],
        'Q' => &[
            &[
                (2.0, 0.0),
                (0.0, 2.0),
                (0.0, 7.0),
                (2.0, 9.0),
                (4.0, 9.0),
                (6.0, 7.0),
                (6.0, 2.0),
                (4.0, 0.0),
                (2.0, 0.0),
            ],
            &[(4.0, 2.0), (6.0, 0.0)],
        ],
        'R' => &[
            &[
                (0.0, 0.0),
                (0.0, 9.0),
                (5.0, 9.0),
                (6.0, 8.0),
                (6.0, 6.0),
                (5.0, 5.0),
                (0.0, 5.0),
            ],
            &[(2.0, 5.0), (6.0, 0.0)],
        ],
        'S' => &[&[
            (6.0, 7.0),
            (4.0, 9.0),
            (2.0, 9.0),
            (0.0, 7.5),
            (0.0, 6.0),
            (2.0, 5.0),
            (4.0, 5.0),
            (6.0, 3.5),
            (6.0, 2.0),
            (4.0, 0.0),
            (2.0, 0.0),
            (0.0, 2.0),
        ]],
        'T' => &[&[(3.0, 0.0), (3.0, 9.0)], &[(0.0, 9.0), (6.0, 9.0)]],
        'U' => &[&[
            (0.0, 9.0),
            (0.0, 2.0),
            (2.0, 0.0),
            (4.0, 0.0),
            (6.0, 2.0),
            (6.0, 9.0),
        ]],
        'V' => &[&[(0.0, 9.0), (3.0, 0.0), (6.0, 9.0)]],
        'W' => &[&[(0.0, 9.0), (1.5, 0.0), (3.0, 6.0), (4.5, 0.0), (6.0, 9.0)]],
        'X' => &[&[(0.0, 0.0), (6.0, 9.0)], &[(0.0, 9.0), (6.0, 0.0)]],
        'Y' => &[
            &[(0.0, 9.0), (3.0, 4.0), (6.0, 9.0)],
            &[(3.0, 4.0), (3.0, 0.0)],
        ],
        'Z' => &[&[(0.0, 9.0), (6.0, 9.0), (0.0, 0.0), (6.0, 0.0)]],
        '0' => &[
            &[
                (2.0, 0.0),
                (0.0, 2.0),
                (0.0, 7.0),
                (2.0, 9.0),
                (4.0, 9.0),
                (6.0, 7.0),
                (6.0, 2.0),
                (4.0, 0.0),
                (2.0, 0.0),
            ],
            &[(1.0, 1.0), (5.0, 8.0)],
        ],
        '1' => &[
            &[(1.0, 7.0), (3.0, 9.0), (3.0, 0.0)],
            &[(1.0, 0.0), (5.0, 0.0)],
        ],
        '2' => &[&[
            (0.0, 7.0),
            (2.0, 9.0),
            (4.0, 9.0),
            (6.0, 7.0),
            (6.0, 5.0),
            (0.0, 0.0),
            (6.0, 0.0),
        ]],
        '3' => &[&[
            (0.0, 8.0),
            (2.0, 9.0),
            (4.0, 9.0),
            (6.0, 7.0),
            (6.0, 6.0),
            (4.0, 5.0),
            (6.0, 4.0),
            (6.0, 2.0),
            (4.0, 0.0),
            (2.0, 0.0),
            (0.0, 1.0),
        ]],
        '4' => &[&[(4.0, 0.0), (4.0, 9.0), (0.0, 3.0), (6.0, 3.0)]],
        '5' => &[&[
            (6.0, 9.0),
            (0.0, 9.0),
            (0.0, 5.0),
            (4.0, 5.0),
            (6.0, 3.0),
            (6.0, 2.0),
            (4.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0),
        ]],
        '6' => &[&[
            (5.0, 9.0),
            (2.0, 9.0),
            (0.0, 6.0),
            (0.0, 2.0),
            (2.0, 0.0),
            (4.0, 0.0),
            (6.0, 2.0),
            (6.0, 3.0),
            (4.0, 5.0),
            (0.0, 5.0),
        ]],
        '7' => &[&[(0.0, 9.0), (6.0, 9.0), (2.0, 0.0)]],
        '8' => &[&[
            (2.0, 5.0),
            (0.0, 6.5),
            (0.0, 7.5),
            (2.0, 9.0),
            (4.0, 9.0),
            (6.0, 7.5),
            (6.0, 6.5),
            (4.0, 5.0),
            (2.0, 5.0),
            (0.0, 3.5),
            (0.0, 1.5),
            (2.0, 0.0),
            (4.0, 0.0),
            (6.0, 1.5),
            (6.0, 3.5),
            (4.0, 5.0),
        ]],
        '9' => &[&[
            (1.0, 0.0),
            (4.0, 0.0),
            (6.0, 3.0),
            (6.0, 7.0),
            (4.0, 9.0),
            (2.0, 9.0),
            (0.0, 7.0),
            (0.0, 6.0),
            (2.0, 4.0),
            (6.0, 4.0),
        ]],
        '-' => &[&[(1.0, 4.0), (5.0, 4.0)]],
        '.' => &[&[(2.6, 0.0), (3.4, 0.0), (3.4, 0.8), (2.6, 0.8), (2.6, 0.0)]],
        _ => return None,
    };
    Some(glyph)
}

/// Build `text` as engraved wires in the `xy` plane of `frame`, `height`
/// tall at the capitals, reading along `x`.
///
/// The result is a compound of open wires — one per stroke — with history
/// generating each from nothing. A space advances; any other character the
/// font does not carry refuses by name.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
/// height is not finite and positive, or a character is not in the font.
pub fn make_text(
    model: &mut Model,
    text: &str,
    frame: &Frame,
    height: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if !height.is_finite() || height <= tol.confusion() {
        ogeom_bail!(Construction, "text {height} tall is not text");
    }
    let scale = height / CAP;
    let place = |pen: f64, x: f64, y: f64| -> Point {
        frame.origin() + frame.x().vector() * ((pen + x) * scale) + frame.y().vector() * (y * scale)
    };

    let mut wires: Vec<Shape> = Vec::new();
    let mut pen = 0.0;
    for c in text.chars() {
        if c == ' ' {
            pen += ADVANCE;
            continue;
        }
        let Some(glyph) = glyph_of(c) else {
            ogeom_bail!(
                Construction,
                "the built-in font does not carry {c:?}; the supported set is A-Z, 0-9, \
                 space, dash and dot"
            );
        };
        for stroke in glyph {
            if stroke.len() < 2 {
                continue;
            }
            let vertices: Vec<Shape> = stroke
                .iter()
                .map(|(x, y)| make_vertex(model, place(pen, *x, *y)).shape)
                .collect();
            let mut edges = Vec::with_capacity(stroke.len() - 1);
            for i in 0..stroke.len() - 1 {
                let a = place(pen, stroke[i].0, stroke[i].1);
                let b = place(pen, stroke[i + 1].0, stroke[i + 1].1);
                edges.push(
                    make_edge_between(
                        model,
                        LineCurve::segment(a, b, tol)?.into(),
                        (0.0, a.distance(b)),
                        &vertices[i],
                        &vertices[i + 1],
                        tol,
                    )?
                    .shape,
                );
            }
            wires.push(make_wire(model, &edges, tol)?.shape);
        }
        pen += ADVANCE;
    }

    let compound = make_compound(model, &wires)?;
    let mut history = History::new();
    for wire in &wires {
        history.generate(wire, compound.shape.clone());
    }
    Ok(Built::new(compound.shape, history))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use ogeom_topo::{ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn text_engraves_at_its_height_and_refuses_what_it_lacks() {
        let mut model = Model::new();
        let built = make_text(&mut model, "OG-42.0", &Frame::WORLD, 12.0, T).unwrap();
        let wires = explore_unique(&model, &built.shape, ShapeType::Wire).unwrap();
        // O(1) G(1) -(1) 4(1) 2(1) .(1) 0(2) strokes.
        assert_eq!(wires.len(), 8, "one wire per stroke");

        // The capitals reach exactly the requested height, and nothing
        // rises above it.
        let mut top = f64::NEG_INFINITY;
        for v in explore_unique(&model, &built.shape, ShapeType::Vertex).unwrap() {
            let Some(data) = model.node(&v).and_then(|n| n.data().as_vertex()) else {
                continue;
            };
            top = top.max(data.point.z.mul_add(0.0, data.point.y));
        }
        assert!((top - 12.0).abs() < 1e-9, "cap height {top}");

        let unsupported = make_text(&mut model, "π", &Frame::WORLD, 12.0, T);
        assert!(unsupported.is_err(), "an honest refusal, not a blank");
    }
}
