//! Holes drilled beside the corner of a box rounded on every edge. A hole
//! crosses a fillet's tangent line on a flat face, where the section on
//! the flat and the section on the fillet meet tangentially, and near the
//! corner it crosses fillets and the corner ball together.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_box, make_cylinder, volume_properties};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

/// Cut and common with each hole are valid and share the box's volume.
#[test]
fn holes_beside_a_rounded_corner_cut_and_fill_valid() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let rounded = ogeom::fillet::fillet_edges(&mut model, &block, &edges, 2.0, T)
        .unwrap()
        .shape;
    let whole = volume(&model, &rounded);
    for (x, y, r) in [
        (4.3, 2.3, 1.4),
        (1.1, 4.45, 1.0),
        (0.65, 4.5, 0.75),
        (1.1, 0.6, 1.4),
        (1.35, 1.9, 1.4),
        (4.6, 3.0, 1.3),
    ] {
        let at = Frame::new(Point::new(x, y, -5.0), Direction::Z, Direction::X, T).unwrap();
        let drill = make_cylinder(&mut model, at, r, 20.0, T).unwrap().shape;
        let mut shares = 0.0;
        for (name, made) in [
            ("cut", ogeom::boolean::cut(&mut model, &rounded, &drill, T)),
            (
                "common",
                ogeom::boolean::common(&mut model, &rounded, &drill, T),
            ),
        ] {
            let made = made.unwrap_or_else(|e| panic!("{name} at ({x}, {y}): {e}"));
            let diagnosis = check(&model, &made.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{name} at ({x}, {y}): {diagnosis}");
            shares += volume(&model, &made.shape);
        }
        assert!(
            (shares - whole).abs() < whole * 1e-4,
            "at ({x}, {y}): {shares} against {whole}"
        );
    }
}
