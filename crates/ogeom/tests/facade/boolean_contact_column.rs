//! Cuts through an imported inside corner, a fraction of a micron out of
//! square, by a band whose cylinder the part already carries from an
//! earlier cut. The contact along the corner's wall runs level with the
//! cylinder face's edge in its chart; its image is folded onto the face's
//! own column, and the face it bounds closes.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"))
}

/// The part less the band: valid, and as much smaller as the band's
/// volume inside the part, measured by the common of the two.
fn cut_closes(name: &str) {
    let (mut model, roots) = ogeom::io::native::read(&corpus(name)).unwrap();
    let (part, band) = (&roots[0], &roots[1]);
    let volume = |model: &ogeom::topo::Model, shape: &ogeom::topo::Shape| {
        ogeom::algo::volume_properties(model, shape, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass
    };
    let cut = ogeom::boolean::cut(&mut model, part, band, T)
        .unwrap_or_else(|e| panic!("{name}: {e}"))
        .shape;
    let diagnosis = ogeom::algo::check(&model, &cut, T).unwrap();
    assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
    let inside = ogeom::boolean::common(&mut model, part, band, T)
        .unwrap_or_else(|e| panic!("{name}, common: {e}"))
        .shape;
    let (whole, kept, taken) = (
        volume(&model, part),
        volume(&model, &cut),
        volume(&model, &inside),
    );
    assert!(
        (whole - kept - taken).abs() <= whole * 1e-9,
        "{name}: the cut kept {kept} of {whole}, the band takes {taken}"
    );
}

#[test]
fn a_cut_through_a_corner_whose_vertex_stands_off_the_wall_closes() {
    cut_closes("inner_corner_off_wall_second_cut.ogeom");
}

#[test]
fn a_cut_through_a_corner_whose_edge_runs_past_its_vertex_closes() {
    cut_closes("inner_corner_overrun_second_cut.ogeom");
}
