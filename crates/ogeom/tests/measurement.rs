//! Measuring what the operations leave behind: a result whose faces are all
//! analytic is measured in closed form, however the arrangement split its
//! rims on the way.
//!
//! The exact path reports a deflection of zero, and that is the pin here:
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
/// at all, which planar faces do not need, since a mesh of a flat thing is
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
/// start (including the rim at the far end, which the cut never reached),
/// so a drum that arrives with one full-turn circle on each disc leaves
/// with two arcs on each. The disc those arcs bound is the same disc. And
/// the wall's seam keeps the column it was born with, so its chart's hull
/// stands a whole millimetre above the rim the cut left: the seam is read
/// over the edge's own range, not the hull's.
///
/// Read off a mesh instead, the drum is out by 1.3% at the default chord,
/// which is the kind of error that looks like a tolerance and is not one.
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
/// same way, the torus included.
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

/// A face with a hole in it is a region less a region, and the integral is
/// their sum: a plate with a bore is a rectangle less a disc, a tube's end
/// face a disc less a disc. A mesh of a circle is a polygon inscribed in
/// it, so measured off a mesh both come out heavy.
#[test]
fn a_face_with_a_hole_is_weighed_in_closed_form() {
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
    let measured = weigh(&model, &bored, true);
    let want = 20.0 * 20.0 * 10.0 - pi * 16.0 * 10.0;
    assert!(
        (measured - want).abs() < 1e-8,
        "a plate with a bore: {measured} against {want}"
    );

    let outer = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 9.0, 4.0, T)
        .unwrap()
        .shape;
    let bore = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 4.0, T)
        .unwrap()
        .shape;
    let tube = ogeom::boolean::cut(&mut model, &outer, &bore, T)
        .unwrap()
        .shape;
    let measured = weigh(&model, &tube, true);
    let want = pi * (81.0 - 25.0) * 4.0;
    assert!(
        (measured - want).abs() < 1e-8,
        "a tube: {measured} against {want}"
    );
}

/// A shell whose faces disagree about which way is out is not weighed in
/// closed form, however analytic its surfaces are.
///
/// The flag on a face is the only thing that says which side of its surface
/// the material is on, and nothing in a shell makes the flags agree. They
/// can be asked about each other, though: an edge between two faces is
/// walked by each with its own material on its left, so the two walks run
/// opposite ways along it. Where they do not, this hands the shape to the
/// tessellator, which repairs such a shell by flipping whichever side of
/// the disagreement is in the minority.
///
/// An imported part can arrive exactly like this (a bore wall's flag
/// pointing into the solid), and taking the flags at their word counts the
/// bore as material and makes the part a third heavy.
#[test]
fn a_shell_whose_faces_disagree_is_left_to_the_mesh() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 6.0, 4.0), T)
        .unwrap()
        .shape;
    let honest = ogeom::algo::volume_properties(&model, &block, Deflection::default(), T).unwrap();
    assert_eq!(honest.deflection, 0.0, "a box is weighed exactly");
    assert!((honest.mass - 240.0).abs() < 1e-9);

    // Every face in turn, since which one is turned over decides how wrong
    // believing it would be, and one of them is the plane the moments are
    // measured from, where believing it costs nothing at all.
    for which in 0..6 {
        let mut faces =
            ogeom::topo::explore_unique(&model, &block, ogeom::topo::ShapeType::Face).unwrap();
        faces[which] = faces[which].clone().reversed();
        let shell = ogeom::algo::make_shell(&mut model, &faces).unwrap().shape;
        let turned = ogeom::algo::make_solid(&mut model, std::slice::from_ref(&shell))
            .unwrap()
            .shape;
        let measured =
            ogeom::algo::volume_properties(&model, &turned, Deflection::default(), T).unwrap();
        assert!(
            measured.deflection > 0.0,
            "face {which} turned over: the mesh was asked instead"
        );
        assert!(
            (measured.mass - 240.0).abs() < 1e-9,
            "face {which} turned over: and it mends the flag, {}",
            measured.mass
        );
    }
}

/// A face the edge walks cannot compare (a cone cap whose apex has no
/// normal) is asked by probing the solid off both its sides, and turned
/// over it is left to the mesh like any other disagreement, never weighed
/// in closed form with its share of the volume counted backwards.
///
/// The cap is all but flat, so counting it backwards is only three parts
/// in a thousand wrong: the mesh is asked at a chord fine enough to land
/// within one.
#[test]
fn a_face_the_walks_cannot_compare_is_probed_before_it_is_believed() {
    use ogeom::topo::{ShapeType, explore_unique};
    let mut model = Model::new();
    // A drum of radius 5 and height 10 along `y`, its top a cone rising
    // 0.05 to the axis.
    let (radius, height, rise) = (5.0, 10.0, 0.05);
    let pts = [
        (0.0, 0.0),
        (radius, 0.0),
        (radius, height),
        (0.0, height + rise),
    ]
    .map(|(x, y)| Point::new(x, y, 0.0));
    let wire = ogeom::algo::make_polygon(&mut model, &pts, true, T)
        .unwrap()
        .shape;
    let plane = ogeom::geom::PlaneSurface::new(ogeom::math::Plane::new(Frame::WORLD)).into();
    let profile = ogeom::algo::make_face(&mut model, plane, &[wire], T)
        .unwrap()
        .shape;
    let axis = ogeom::math::Axis {
        location: Point::ORIGIN,
        direction: Direction::Y,
    };
    let drum = ogeom::algo::make_revolution(&mut model, &profile, axis, std::f64::consts::TAU, T)
        .unwrap()
        .shape;
    let pi = std::f64::consts::PI;
    let want = pi * radius * radius * (height + rise / 3.0);
    let honest = ogeom::algo::volume_properties(&model, &drum, Deflection::default(), T).unwrap();
    assert_eq!(honest.deflection, 0.0, "the drum is weighed exactly");
    assert!(
        (honest.mass - want).abs() < 1e-9,
        "{} against {want}",
        honest.mass
    );

    let faces = explore_unique(&model, &drum, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 3);
    for which in 0..faces.len() {
        let mut held = faces.clone();
        held[which] = held[which].reversed();
        let shell = ogeom::algo::make_shell(&mut model, &held).unwrap().shape;
        let turned = ogeom::algo::make_solid(&mut model, &[shell]).unwrap().shape;
        let measured = ogeom::algo::volume_properties(
            &model,
            &turned,
            Deflection::with_chord(1e-3).unwrap(),
            T,
        )
        .unwrap();
        assert!(
            measured.deflection > 0.0,
            "face {which} turned over: the mesh was asked, not the closed form ({})",
            measured.mass
        );
        assert!(
            (measured.mass - want).abs() < want * 1e-3,
            "face {which} turned over: and it mends the flag, {} against {want}",
            measured.mass
        );
    }
}

/// A lump is one set of faces that agree, and any face of it says which way
/// the set faces. A sheet thinner than the probe's shortest step has only
/// outside either side of its broad faces, so where one of those is asked
/// first and cannot tell, its rim is asked next, and two lumps that face
/// out are still weighed in closed form.
#[test]
fn a_sheet_whose_broad_faces_cannot_tell_is_settled_by_its_rim() {
    let mut model = Model::new();
    let thick = 5e-4;
    // The sheet thin along `z`, so its first face is a broad one.
    let sheet = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, thick), T)
        .unwrap()
        .shape;
    let drum = ogeom::algo::make_cylinder(
        &mut model,
        Frame::new(Point::new(20.0, 5.0, 0.0), Direction::Z, Direction::X, T).unwrap(),
        5.0,
        10.0,
        T,
    )
    .unwrap()
    .shape;
    let both = ogeom::algo::make_compound(&mut model, &[sheet, drum])
        .unwrap()
        .shape;
    let want = thick * 100.0 + std::f64::consts::PI * 250.0;
    let measured = ogeom::algo::volume_properties(&model, &both, Deflection::default(), T).unwrap();
    assert_eq!(measured.deflection, 0.0, "the closed form was taken");
    assert!(
        (measured.mass - want).abs() < 1e-9,
        "{} against {want}",
        measured.mass
    );
}

#[test]
fn a_metre_cube_and_a_long_drum_measure_exactly_and_at_once() {
    let mut model = Model::new();
    let cube = ogeom::algo::make_box(&mut model, Frame::WORLD, (1000.0, 1000.0, 1000.0), T)
        .unwrap()
        .shape;
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 3000.0, T)
        .unwrap()
        .shape;
    let started = std::time::Instant::now();
    for (shape, volume, area) in [
        (&cube, 1e9, 6e6),
        (
            &drum,
            core::f64::consts::PI * 25.0 * 3000.0,
            core::f64::consts::TAU * 5.0 * (5.0 + 3000.0),
        ),
    ] {
        let v = ogeom::algo::volume_properties(&model, shape, Deflection::default(), T).unwrap();
        let a = ogeom::algo::surface_properties(&model, shape, Deflection::default(), T).unwrap();
        assert_eq!((v.deflection, a.deflection), (0.0, 0.0));
        assert!(
            (v.mass - volume).abs() <= volume * 1e-12,
            "{} against {volume}",
            v.mass
        );
        assert!(
            (a.mass - area).abs() <= area * 1e-12,
            "{} against {area}",
            a.mass
        );
    }
    // A face a metre across is a handful of panels, not thousands.
    assert!(
        started.elapsed().as_secs_f64() < 2.0,
        "{:?}",
        started.elapsed()
    );
}

/// A void is a shell of its own, so turned inside out as a whole (every
/// face of it reversed) no edge notices: the faces agree with each other
/// and the mesh, which mends a minority of turned faces within each piece,
/// would weigh the cavity as material. Each shell of a solid with several
/// is asked which way it faces: the whole void turned over is turned back
/// and weighed in closed form, and a block with a face of its own turned as
/// well is meshed with every face asked, the void still subtracted.
#[test]
fn a_void_turned_inside_out_is_still_a_void() {
    use ogeom::topo::{ShapeType, explore_unique};
    let chord = Deflection::with_chord(1e-3).unwrap();
    for spherical in [false, true] {
        let mut model = Model::new();
        let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let (cavity, void) = if spherical {
            let centre =
                Frame::new(Point::new(5.0, 5.0, 5.0), Direction::Z, Direction::X, T).unwrap();
            (
                ogeom::algo::make_sphere(&mut model, centre, 2.0, T)
                    .unwrap()
                    .shape,
                4.0 / 3.0 * std::f64::consts::PI * 8.0,
            )
        } else {
            let corner =
                Frame::new(Point::new(3.0, 3.0, 3.0), Direction::Z, Direction::X, T).unwrap();
            (
                ogeom::algo::make_box(&mut model, corner, (4.0, 4.0, 4.0), T)
                    .unwrap()
                    .shape,
                64.0,
            )
        };
        let want = 1000.0 - void;
        let outer = explore_unique(&model, &block, ShapeType::Face).unwrap();
        // The cavity's own faces face out of it, into the block's material.
        let turned_void = explore_unique(&model, &cavity, ShapeType::Face).unwrap();
        let honest_void: Vec<Shape> = turned_void.iter().map(Shape::reversed).collect();
        let mut turned_outer = outer.clone();
        turned_outer[0] = turned_outer[0].reversed();
        for (name, outer, inner, exact) in [
            ("as built", &outer, &honest_void, true),
            ("void turned", &outer, &turned_void, true),
            ("outer face turned", &turned_outer, &honest_void, false),
            ("both turned", &turned_outer, &turned_void, false),
        ] {
            let o = ogeom::algo::make_shell(&mut model, outer).unwrap().shape;
            let i = ogeom::algo::make_shell(&mut model, inner).unwrap().shape;
            let solid = ogeom::algo::make_solid(&mut model, &[o, i]).unwrap().shape;
            let measured = ogeom::algo::volume_properties(&model, &solid, chord, T).unwrap();
            assert_eq!(
                measured.deflection == 0.0,
                exact,
                "sphere {spherical}, {name}: closed form {exact}"
            );
            assert!(
                (measured.mass - want).abs() < want * 1e-4,
                "sphere {spherical}, {name}: {} against {want}",
                measured.mass
            );
        }
    }
}
