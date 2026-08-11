//! IGES, held to the same standard as every other format here: what goes out
//! comes back, and the measurement says so.
//!
//! Each case builds a solid with the kernel's own operations, writes it as an
//! IGES manifold solid B-rep, reads the deck back, and measures the recovered
//! solid's volume against the original's — not against a hope, against the
//! number. The corpus covers the surface vocabulary a real file exercises:
//! planes, a periodic cylinder wall with its seam, a torus, a boolean result
//! carrying both, and a B-spline loft.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// Write one solid, read it back, and return (original volume, recovered
/// volume, the import) for the assertions each case cares about.
fn round_trip(mut model: Model, solid: Shape) -> (f64, f64, ogeom::io::IgesImport) {
    let original = volume(&model, &solid);
    let mut document = ogeom::doc::Document::over(std::mem::take(&mut model));
    document.add_part("part", solid);
    let text = ogeom::io::write_iges(&document, T).unwrap();
    let import = ogeom::io::read_iges(&text, T).unwrap();
    assert_eq!(import.solids.len(), 1, "one solid out, one back");
    let recovered = volume(import.document.model(), &import.solids[0]);
    (original, recovered, import)
}

#[test]
fn a_box_round_trips_to_its_own_volume() {
    let mut model = Model::new();
    let solid = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 12.0, 8.0), T)
        .unwrap()
        .shape;
    let (original, recovered, import) = round_trip(model, solid);
    assert!((original - 20.0 * 12.0 * 8.0).abs() < 1e-9);
    assert!(
        (recovered - original).abs() < 1e-9,
        "{recovered} against {original}"
    );
    assert!(
        import.report.warnings.is_empty(),
        "{:?}",
        import.report.warnings
    );
}

#[test]
fn a_cylinder_round_trips_with_its_seam() {
    let mut model = Model::new();
    let solid = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 6.0, 15.0, T)
        .unwrap()
        .shape;
    let (original, recovered, _) = round_trip(model, solid);
    // The wall wraps its whole period; a reader that lost the seam would
    // refuse to tessellate at all, and one that miscounted the turn would
    // miss by a factor, not an epsilon.
    assert!(
        (recovered - original).abs() < original * 1e-6,
        "{recovered} against {original}"
    );
}

#[test]
fn a_torus_round_trips_doubly_periodic() {
    let mut model = Model::new();
    let solid = ogeom::algo::make_torus(&mut model, Frame::WORLD, 10.0, 3.0, T)
        .unwrap()
        .shape;
    let (original, recovered, _) = round_trip(model, solid);
    assert!(
        (recovered - original).abs() < original * 1e-6,
        "{recovered} against {original}"
    );
}

#[test]
fn a_drilled_block_round_trips_through_its_boolean_faces() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (30.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(15.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, frame, 4.0, 12.0, T)
        .unwrap()
        .shape;
    let solid = ogeom::boolean::cut(&mut model, &block, &drill, T)
        .unwrap()
        .shape;
    let (original, recovered, _) = round_trip(model, solid);
    assert!(
        (recovered - original).abs() < original * 1e-6,
        "{recovered} against {original}"
    );
}

#[test]
fn a_spline_walled_prism_round_trips_through_126_and_128() {
    // A lens-shaped profile from two fitted open B-spline arcs, extruded:
    // the walls are extrusion surfaces over B-spline curves, which the
    // writer spells as rational B-spline surfaces (128) bounded by rational
    // B-spline curves (126), and the fitted-pcurve path carries back in.
    //
    // The walls are open on purpose. A wall that wraps a whole period — a
    // closed skinned loft — does not survive *either* exchange format today,
    // and docs/PLAN.md F5 owns that; a test here could only hold IGES to a
    // standard STEP does not meet.
    let mut model = Model::new();
    let arc = |sign: f64| -> ogeom::geom::Curve {
        let pts: Vec<Point> = (0..=8)
            .map(|i| {
                let t = f64::from(i) / 8.0;
                let x = -10.0 + 20.0 * t;
                let y = sign * 6.0 * (core::f64::consts::PI * t).sin();
                Point::new(x, y, 0.0)
            })
            .collect();
        ogeom::geom::Curve::from(
            ogeom::geom::fit::fit_points(&pts, 3, 1e-9, T)
                .unwrap()
                .curve,
        )
    };
    let upper = arc(1.0);
    let lower = arc(-1.0);
    let va = ogeom::algo::make_vertex(&mut model, Point::new(-10.0, 0.0, 0.0)).shape;
    let vb = ogeom::algo::make_vertex(&mut model, Point::new(10.0, 0.0, 0.0)).shape;
    let (u0, u1) = ogeom::geom::Curve3d::domain(&upper);
    let (l0, l1) = ogeom::geom::Curve3d::domain(&lower);
    let e_up = ogeom::algo::make_edge_between(&mut model, upper, (u0, u1), &va, &vb, T)
        .unwrap()
        .shape;
    let e_dn = ogeom::algo::make_edge_between(&mut model, lower, (l0, l1), &va, &vb, T)
        .unwrap()
        .shape;
    let face = ogeom::algo::make_face_with_pcurves(
        &mut model,
        ogeom::geom::SurfaceGeometry::Plane(ogeom::geom::PlaneSurface::new(
            ogeom::math::Plane::new(Frame::WORLD),
        )),
        // Lower arc first, then the upper against its sense: counter-
        // clockwise about +z, so the prism winds outward.
        &[vec![e_dn, e_up.reversed()]],
        T,
    )
    .unwrap()
    .shape;
    let solid = ogeom::algo::make_prism(
        &mut model,
        &face,
        ogeom::math::Vector::new(0.0, 0.0, 12.0),
        T,
    )
    .unwrap()
    .shape;
    let (original, recovered, _) = round_trip(model, solid);
    assert!(
        (recovered - original).abs() < original * 1e-6,
        "{recovered} against {original}"
    );
}

#[test]
fn a_sphere_round_trips_through_its_seam_only_boundary() {
    // A full sphere's only written boundary is its seam, listed twice; the
    // poles have no curve to write. The reader recognises a seam-only loop
    // as the whole chart and rebuilds the natural face, degenerate boundary
    // and all. The two tessellations sample the chart differently, so the
    // comparison is against the closed form, at a chord fine enough that
    // the inscribed deficit stays inside the stated bound.
    let mut model = Model::new();
    let solid = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 7.0, T)
        .unwrap()
        .shape;
    let mut document = ogeom::doc::Document::over(std::mem::take(&mut model));
    document.add_part("ball", solid);
    let text = ogeom::io::write_iges(&document, T).unwrap();
    let import = ogeom::io::read_iges(&text, T).unwrap();
    let fine = Deflection::with_chord(0.01).unwrap();
    let recovered =
        ogeom::algo::volume_properties(import.document.model(), &import.solids[0], fine, T)
            .unwrap()
            .mass;
    let exact = 4.0 / 3.0 * core::f64::consts::PI * 7.0_f64.powi(3);
    // An inscribed tessellation at chord ε on radius r reads low by about
    // 3ε/r — every triangle sits up to ε inside — which at 0.01 on 7 is
    // 0.43%. The bound is that, with a third again for the sampling's
    // unevenness, and it is a statement about tessellation rather than
    // about the exchange: the same solid measured natively reads the same.
    let budget = 3.0 * 0.01 / 7.0 * (4.0 / 3.0);
    assert!(
        (recovered - exact).abs() < exact * budget,
        "{recovered} against {exact}, budget {budget}"
    );
}

#[test]
fn an_inch_file_scales_into_millimetres() {
    let mut model = Model::new();
    let solid = ogeom::algo::make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
        .unwrap()
        .shape;
    let mut document = ogeom::doc::Document::over(std::mem::take(&mut model));
    document.add_part("cube", solid);
    let text = ogeom::io::write_iges(&document, T).unwrap();
    // The same deck, restated in inches: the global units flag and name are
    // the only difference, and every coordinate must arrive 25.4 times as
    // large. The substitution is same-length so the fixed columns survive.
    let inch = text.replace(",2,2HMM,", ",1,2HIN,");
    assert_ne!(inch, text, "the substitution must strike");
    let import = ogeom::io::read_iges(&inch, T).unwrap();
    assert!((import.report.scale_mm - 25.4).abs() < 1e-12);
    let recovered = volume(import.document.model(), &import.solids[0]);
    let expected = 25.4_f64.powi(3);
    assert!(
        (recovered - expected).abs() < expected * 1e-9,
        "{recovered} against {expected}"
    );
}

#[test]
fn an_empty_deck_is_refused_by_name() {
    let err = ogeom::io::read_iges("not iges at all\n", T).unwrap_err();
    assert!(
        err.to_string().contains("section letter"),
        "unexpected message: {err}"
    );
}
