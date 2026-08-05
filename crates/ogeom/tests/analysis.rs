//! The §6 analysis tail: draft against a pull direction, least material
//! thickness, and self-intersection — each answering what it measures and
//! stating how it sampled.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::{LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Plane, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_slab_reads_its_draft_and_its_thickness() {
    let mut model = Model::new();
    let slab = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 5.0), T)
        .unwrap()
        .shape;
    let scene = ogeom::select::Pickable::build(&model, &slab, Deflection::default(), T).unwrap();

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

#[test]
fn crossing_sheets_confess_and_a_box_does_not() {
    let mut model = Model::new();
    let box_shape = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    assert!(
        ogeom::algo::check_self_intersection(&model, &box_shape, T)
            .unwrap()
            .is_empty(),
        "a valid solid interferes with nothing"
    );

    // Two square sheets crossing like an X, sharing no topology.
    let sheet = |model: &mut Model, frame: Frame| -> Shape {
        let corners = [
            frame.origin() + frame.x().vector() * -5.0 + frame.y().vector() * -5.0,
            frame.origin() + frame.x().vector() * 5.0 + frame.y().vector() * -5.0,
            frame.origin() + frame.x().vector() * 5.0 + frame.y().vector() * 5.0,
            frame.origin() + frame.x().vector() * -5.0 + frame.y().vector() * 5.0,
        ];
        let vertices: Vec<Shape> = corners
            .iter()
            .map(|c| ogeom::algo::make_vertex(model, *c).shape)
            .collect();
        let edges: Vec<Shape> = (0..4)
            .map(|i| {
                let (a, b) = (corners[i], corners[(i + 1) % 4]);
                ogeom::algo::make_edge_between(
                    model,
                    LineCurve::segment(a, b, T).unwrap().into(),
                    (0.0, a.distance(b)),
                    &vertices[i],
                    &vertices[(i + 1) % 4],
                    T,
                )
                .unwrap()
                .shape
            })
            .collect();
        ogeom::algo::make_face_with_pcurves(
            model,
            SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(frame))),
            &[edges],
            T,
        )
        .unwrap()
        .shape
    };
    let flat = sheet(&mut model, Frame::WORLD);
    let tilted_frame = Frame::new(
        Point::new(0.0, 0.0, -2.0),
        Direction::from_coords(0.0, 1.0, 0.2, T).unwrap(),
        Direction::X,
        T,
    )
    .unwrap();
    let tilted = sheet(&mut model, tilted_frame);
    let pair = ogeom::algo::make_compound(&mut model, &[flat, tilted])
        .unwrap()
        .shape;

    let crossings = ogeom::algo::check_self_intersection(&model, &pair, T).unwrap();
    assert_eq!(crossings.len(), 1, "the sheets cross once: {crossings:?}");
}
