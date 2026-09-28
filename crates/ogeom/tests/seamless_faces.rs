//! Periodic faces imported without the seam or pole edge their chart needs:
//! bands between rings that each go round once, caps bounded by one ring,
//! faces reaching a pole with no edge along the pole's row. Healed, they
//! close in the chart, face the way their loops say, and a boolean cuts
//! them.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::with_chord(0.05).unwrap(), T)
        .unwrap()
        .mass
}

/// A part whose bands, caps and pole corners all arrived seamless: healed,
/// it is valid and measures the same, and a drill through a spherical
/// pocket and the bands around it cuts, the cut and the common adding up to
/// the part.
#[test]
fn seamless_faces_heal_and_a_drill_cuts_through_them() {
    let mut import = ogeom::io::read_step(&corpus("nist_ftc_06_asme1_rd.stp"), T).unwrap();
    let solid = import.solids[0].clone();
    let model = import.document.model_mut();
    let before = volume(model, &solid);
    let (healed, count) = ogeom::heal::reanchor_periodic_rings(model, &solid, T).unwrap();
    assert!(count > 0);
    let part = healed.shape;
    assert!(ogeom::algo::check(model, &part, T).unwrap().is_valid());
    let whole = volume(model, &part);
    assert!(
        (whole - before).abs() < before * 1e-4,
        "{whole} against {before}"
    );

    // The pocket's floor faces into the sphere: the material is outside it.
    let beside_pole = Point::new(114.3, 46.625, -107.95);
    let within_pocket = Point::new(114.3, 48.625, -107.95);
    let fine = Deflection::with_chord(0.01).unwrap();
    for (at, want) in [
        (beside_pole, ogeom::algo::Containment::In),
        (within_pocket, ogeom::algo::Containment::Out),
    ] {
        let found = ogeom::algo::classify_in_solid(model, &part, at, fine, T).unwrap();
        assert_eq!(found, want, "at {at:?}");
    }

    let frame = Frame::new(
        Point::new(
            81.848_489_699_818_77,
            80.088_274_165_849_1,
            -453.028_484_344_653,
        ),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let drill = ogeom::algo::make_cylinder(model, frame, 42.991_498, 493.500_258, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(model, &part, &drill, T).unwrap().shape;
    let common = ogeom::boolean::common(model, &part, &drill, T)
        .unwrap()
        .shape;
    for result in [&cut, &common] {
        assert!(ogeom::algo::check(model, result, T).unwrap().is_valid());
    }
    let (a, b) = (volume(model, &cut), volume(model, &common));
    assert!(
        (a + b - whole).abs() < whole * 1e-4,
        "{a} + {b} against {whole}"
    );
}
