//! Measuring what the operations leave behind: a result whose faces are all
//! analytic is measured in closed form, however the arrangement split its
//! rims on the way.
//!
//! The exact path reports a deflection of zero, and that is the pin here —
//! a number that comes out right at one chord and differently at another
//! was read off a mesh, and a part whose faces are a cylinder and two
//! planes should not have to be meshed to be weighed.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

/// Measured twice at chords three orders apart, and held to agreeing.
///
/// `closed_form` asks as well that the measurement was not read off a mesh
/// at all — which planar faces do not need, since a mesh of a flat thing is
/// the flat thing, but a curved one does.
fn weigh(model: &Model, shape: &Shape, closed_form: bool) -> f64 {
    let coarse = ogeom::algo::volume_properties(model, shape, Deflection::default(), T).unwrap();
    let fine =
        ogeom::algo::volume_properties(model, shape, Deflection::with_chord(1e-4).unwrap(), T)
            .unwrap();
    if closed_form {
        assert_eq!(coarse.deflection, 0.0, "the exact path was taken");
    }
    assert!(
        (coarse.mass - fine.mass).abs() < coarse.mass * 1e-12,
        "the same at either chord: {} and {}",
        coarse.mass,
        fine.mass
    );
    coarse.mass
}

/// A drum cut down to size still weighs what a drum weighs.
///
/// The arrangement splits a closed rim to give its walker somewhere to
/// start — including the rim at the far end, which the cut never reached —
/// so a drum that arrives with one full-turn circle on each disc leaves
/// with two arcs on each. The disc those arcs bound is the same disc. And
/// the wall's seam keeps the column it was born with, so its chart's hull
/// stood a whole millimetre above the rim the cut left: the seam is read
/// over the edge's own range now.
///
/// It was out by 1.3% at the default chord, which is the kind of error
/// that looks like a tolerance and is not one.
#[test]
fn a_boolean_result_is_still_weighed_in_closed_form() {
    let pi = core::f64::consts::PI;
    let mut model = Model::new();

    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T)
        .unwrap()
        .shape;
    let whole = weigh(&model, &drum, true);
    assert!(
        (whole - pi * 250.0).abs() < 1e-9,
        "the drum itself: {whole}"
    );

    // A cut that takes the top off, and reaches nothing at the bottom.
    let seat = Frame::new(Point::new(0.0, 0.0, 9.0), Direction::Z, Direction::X, T).unwrap();
    let lid = ogeom::algo::make_cylinder(&mut model, seat, 6.0, 2.0, T)
        .unwrap()
        .shape;
    let shortened = ogeom::boolean::cut(&mut model, &drum, &lid, T)
        .unwrap()
        .shape;
    let shorter = weigh(&model, &shortened, true);
    assert!(
        (shorter - pi * 25.0 * 9.0).abs() < 1e-9,
        "a drum nine deep: {shorter}"
    );

    // A planar pair, where the same splitting happens to the box's rims.
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let corner = Frame::new(Point::new(4.0, 4.0, 4.0), Direction::Z, Direction::X, T).unwrap();
    let bite = ogeom::algo::make_box(&mut model, corner, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let notched = ogeom::boolean::cut(&mut model, &block, &bite, T)
        .unwrap()
        .shape;
    // Flat faces, and an L-shaped one among them, which is no chart
    // rectangle: this one is meshed, and a mesh of flat faces is exact.
    let bitten = weigh(&model, &notched, false);
    assert!(
        (bitten - (1000.0 - 216.0)).abs() < 1e-9,
        "a box less its corner: {bitten}"
    );
}

/// And a blend's own faces are analytic, so a blended drum is weighed the
/// same way — the torus included.
#[test]
fn a_rim_blend_leaves_a_solid_that_is_still_weighed_exactly() {
    use ogeom::geom::Curve3d as _;
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T)
        .unwrap()
        .shape;
    // The top rim, asked for away from the seam so the seam cannot answer.
    let rim = ogeom::topo::explore_unique(&model, &drum, ogeom::topo::ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = model
                .node(e)
                .and_then(|n| n.data().as_edge())
                .and_then(ogeom::topo::EdgeData::curve3d)
            else {
                return false;
            };
            let Some(geometry) = model.geometry().curve(*curve) else {
                return false;
            };
            (0..=8).all(|i| {
                let t = range.0 + (range.1 - range.0) * f64::from(i) / 8.0;
                geometry
                    .point_at(t, T)
                    .is_ok_and(|p| (p.z - 10.0).abs() < 1e-9 && (p.x.hypot(p.y) - 5.0).abs() < 1e-6)
            })
        })
        .expect("the drum has its top rim");
    let blended = ogeom::fillet::fillet_edge(&mut model, &drum, &rim, 1.0, T)
        .unwrap()
        .shape;
    let mass = weigh(&model, &blended, true);
    // Pappus over the corner the ball rolled out of: the square less the
    // quarter disc, swept about the axis. Its first moment about the
    // cylinder of ball centres is the square's less the quarter disc's,
    // r³/2 − r³/3, so the region rides at (R − r) + r³/6 over its own area.
    let pi = core::f64::consts::PI;
    let (radius, r) = (5.0_f64, 1.0_f64);
    let area = (1.0 - pi / 4.0) * r * r;
    let removed = 2.0 * pi * ((radius - r) * area + r * r * r / 6.0);
    let want = pi * radius * radius * 10.0 - removed;
    assert!(
        (mass - want).abs() < want * 1e-9,
        "the drum less the wedge its rim shed: {mass} against {want}"
    );
}

/// A face with a hole in it is left to the mesh, for now.
///
/// The integral would be a region less a region — a plate with a bore is a
/// rectangle less a disc — and the arithmetic is easy. What is not is
/// knowing which way each face looks: the closed-form integral reads each
/// face's own flag, and a part in the corpus has flags that disagree with
/// each other, which the tessellator repairs and this cannot. Taking them
/// at their word made that part a third heavy, its bore counted as
/// material, so a holed face waits on a way to tell a wrong flag from a
/// right one. Issue #39 carries it.
#[test]
fn a_face_with_a_hole_waits_for_the_mesh() {
    let pi = core::f64::consts::PI;
    let mut model = Model::new();
    let plate = ogeom::algo::make_box(
        &mut model,
        Frame::new(Point::new(-10.0, -10.0, 0.0), Direction::Z, Direction::X, T).unwrap(),
        (20.0, 20.0, 10.0),
        T,
    )
    .unwrap()
    .shape;
    let seat = Frame::new(Point::new(0.0, 0.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, seat, 4.0, 12.0, T)
        .unwrap()
        .shape;
    let bored = ogeom::boolean::cut(&mut model, &plate, &drill, T)
        .unwrap()
        .shape;
    let measured =
        ogeom::algo::volume_properties(&model, &bored, Deflection::with_chord(1e-4).unwrap(), T)
            .unwrap();
    assert!(measured.deflection > 0.0, "the mesh was asked");
    let want = 20.0 * 20.0 * 10.0 - pi * 16.0 * 10.0;
    assert!(
        (measured.mass - want).abs() < want * 1e-4,
        "and answers within its chord: {} against {want}",
        measured.mass
    );
}
