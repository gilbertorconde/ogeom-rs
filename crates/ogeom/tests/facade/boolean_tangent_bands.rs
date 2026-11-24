//! A straight crease's band carried on past its end into a rim's band, on a
//! cube a drum shaves: the rim's torus and the crease's cylinder are both
//! tangent to the cube's top, and their contact curves on it cross, so the
//! two bands touch each other tangentially there, and the top with them.
//! The section between the bands is marched up to that point and stops a
//! few microns off it, while the two bands' edges on the top cross there in
//! closed form; each edge splits at that crossing.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"))
}

/// The part in common with the band's tool, and the part less it: both
/// valid, and adding up to the part.
#[test]
fn a_band_touching_a_rim_band_tangentially_cuts_the_part() {
    let (mut model, roots) =
        ogeom::io::native::read(&corpus("straight_band_past_a_rim_band.ogeom")).unwrap();
    let (part, tool) = (&roots[0], &roots[1]);
    let volume = |model: &ogeom::topo::Model, shape: &ogeom::topo::Shape| {
        ogeom::algo::volume_properties(model, shape, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass
    };
    let taken = ogeom::boolean::common(&mut model, tool, part, T)
        .unwrap_or_else(|e| panic!("common: {e}"))
        .shape;
    assert!(ogeom::algo::check(&model, &taken, T).unwrap().is_valid());
    let kept = ogeom::boolean::cut(&mut model, part, tool, T)
        .unwrap_or_else(|e| panic!("cut: {e}"))
        .shape;
    assert!(ogeom::algo::check(&model, &kept, T).unwrap().is_valid());
    let (whole, a, b) = (
        volume(&model, part),
        volume(&model, &taken),
        volume(&model, &kept),
    );
    assert!(
        (whole - a - b).abs() <= whole * 1e-9,
        "the part {whole} against {a} taken and {b} kept"
    );
}
