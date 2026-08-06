//! Draft against a pull direction and least material thickness — the two
//! analyses that ride on the pick structure, each answering what it measures
//! and stating how it sampled.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame};
use ogeom_mesh::Deflection;
use ogeom_select::Pickable;
use ogeom_topo::Model;

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_slab_reads_its_draft_and_its_thickness() {
    let mut model = Model::new();
    let slab = ogeom_algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 5.0), T)
        .unwrap()
        .shape;
    let scene = Pickable::build(&model, &slab, Deflection::default(), T).unwrap();

    // Draft against +z: the top face reads a quarter turn, the bottom its
    // negative, and every wall reads zero — straight, undrafted.
    let draft = scene.draft_analysis(Direction::Z);
    assert_eq!(draft.len(), 6);
    let quarter = core::f64::consts::FRAC_PI_2;
    let tops = draft
        .iter()
        .filter(|d| (d.min - quarter).abs() < 1e-9 && (d.max - quarter).abs() < 1e-9)
        .count();
    let bottoms = draft
        .iter()
        .filter(|d| (d.min + quarter).abs() < 1e-9 && (d.max + quarter).abs() < 1e-9)
        .count();
    let walls = draft
        .iter()
        .filter(|d| d.min.abs() < 1e-9 && d.max.abs() < 1e-9)
        .count();
    assert_eq!((tops, bottoms, walls), (1, 1, 4), "{draft:?}");

    // Thickness: the big faces see the opposite wall five away; the thin
    // walls see across the slab's own footprint or the five, whichever
    // their rays strike first — every reading is one of the slab's spans.
    let thickness = scene.thickness_analysis();
    let least = thickness
        .iter()
        .map(|t| t.least)
        .fold(f64::INFINITY, f64::min);
    assert!(
        (least - 5.0).abs() < 1e-6,
        "the slab is five thick: {least}"
    );
}
