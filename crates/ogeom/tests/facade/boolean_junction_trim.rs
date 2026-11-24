//! The junction trim of two rounds meeting at a bore's drafted chord: the
//! chord through the bore's axis, its wall drafted by 0.1. The arc's
//! band, run on past the corner, is kept within either side of the
//! junction (the region the other band's tool and its end plane bound).
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"))
}

/// The band in common with the junction's region: valid, and the two
/// pieces either side of the region (common and cut) adding up to the band.
#[test]
fn a_band_in_common_with_a_junction_through_the_bore_s_axis_closes() {
    let (mut model, roots) =
        ogeom::io::native::read(&corpus("bore_chord_junction_trim.ogeom")).unwrap();
    let (band, region) = (&roots[0], &roots[1]);
    let volume = |model: &ogeom::topo::Model, shape: &ogeom::topo::Shape| {
        ogeom::algo::volume_properties(model, shape, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass
    };
    let kept = ogeom::boolean::common(&mut model, band, region, T)
        .unwrap_or_else(|e| panic!("common: {e}"))
        .shape;
    let diagnosis = ogeom::algo::check(&model, &kept, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let off = ogeom::boolean::cut(&mut model, band, region, T)
        .unwrap_or_else(|e| panic!("cut: {e}"))
        .shape;
    let (whole, inside, outside) = (
        volume(&model, band),
        volume(&model, &kept),
        volume(&model, &off),
    );
    assert!(
        (whole - inside - outside).abs() <= whole * 1e-9,
        "the band {whole} against {inside} in common and {outside} cut"
    );
}
