//! Contact configurations the boolean resolves rather than refuses: curved
//! same-domain pairs unify, and contact confined to lines or points is
//! classified off the contact, so it shares no volume.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn fine() -> Deflection {
    Deflection::with_chord(1e-3).unwrap()
}

fn volume(model: &Model, shape: &ogeom::topo::Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, fine(), T)
        .unwrap()
        .mass
}

/// Coaxial cylinders sharing one surface: flush stack, partial overlap,
/// and the cut: the curved same-domain family, held to closed forms.
#[test]
fn curved_same_domain_pairs_unify() {
    let pi = core::f64::consts::PI;

    // Flush stack: walls meet rim to rim on one infinite cylinder.
    {
        let mut model = Model::new();
        let bottom = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T)
            .unwrap()
            .shape;
        let up = Frame::new(Point::new(0.0, 0.0, 10.0), Direction::Z, Direction::X, T).unwrap();
        let top = ogeom::algo::make_cylinder(&mut model, up, 5.0, 10.0, T)
            .unwrap()
            .shape;
        let fused = ogeom::boolean::fuse(&mut model, &bottom, &top, T)
            .unwrap()
            .shape;
        let v = volume(&model, &fused);
        assert!(
            (v - pi * 25.0 * 20.0).abs() / (pi * 25.0 * 20.0) < 1e-3,
            "the stack is one drum: {v}"
        );
    }

    // Partial overlap: the walls overlap in a band and split each other at
    // the other's rims.
    {
        let mut model = Model::new();
        let bottom = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T)
            .unwrap()
            .shape;
        let up = Frame::new(Point::new(0.0, 0.0, 5.0), Direction::Z, Direction::X, T).unwrap();
        let top = ogeom::algo::make_cylinder(&mut model, up, 5.0, 10.0, T)
            .unwrap()
            .shape;
        let fused = ogeom::boolean::fuse(&mut model, &bottom, &top, T)
            .unwrap()
            .shape;
        let v = volume(&model, &fused);
        assert!(
            (v - pi * 25.0 * 15.0).abs() / (pi * 25.0 * 15.0) < 1e-3,
            "the overlap fuses to one taller drum: {v}"
        );
    }
    {
        let mut model = Model::new();
        let bottom = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T)
            .unwrap()
            .shape;
        let up = Frame::new(Point::new(0.0, 0.0, 5.0), Direction::Z, Direction::X, T).unwrap();
        let top = ogeom::algo::make_cylinder(&mut model, up, 5.0, 10.0, T)
            .unwrap()
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &bottom, &top, T)
            .unwrap()
            .shape;
        let v = volume(&model, &cut);
        assert!(
            (v - pi * 25.0 * 5.0).abs() / (pi * 25.0 * 5.0) < 1e-3,
            "the cut keeps the un-overlapped stub: {v}"
        );
    }
}

/// What a boolean gave, measured: volume, validity, face count and solid
/// count.
struct Measured {
    volume: f64,
    valid: bool,
    faces: usize,
    solids: usize,
}

fn measured(model: &Model, shape: &ogeom::topo::Shape) -> Measured {
    Measured {
        volume: volume(model, shape),
        valid: ogeom::algo::check(model, shape, T).unwrap().is_valid(),
        faces: explore_unique(model, shape, ShapeType::Face).unwrap().len(),
        solids: explore_unique(model, shape, ShapeType::Solid)
            .unwrap()
            .len(),
    }
}

/// Fuse, common and cut of the pair `make` builds, each on a fresh model,
/// held to the expected volume, face count and solid count, and to `check`.
fn three_ways(
    make: &dyn Fn(&mut Model) -> (ogeom::topo::Shape, ogeom::topo::Shape),
    want: [(f64, usize, usize); 3],
) {
    for (op, (volume, faces, solids)) in ["fuse", "common", "cut"].into_iter().zip(want) {
        let mut model = Model::new();
        let (a, b) = make(&mut model);
        let got = match op {
            "fuse" => ogeom::boolean::fuse(&mut model, &a, &b, T),
            "common" => ogeom::boolean::common(&mut model, &a, &b, T),
            _ => ogeom::boolean::cut(&mut model, &a, &b, T),
        }
        .unwrap_or_else(|e| panic!("{op}: {e}"))
        .shape;
        let m = measured(&model, &got);
        assert!(
            (m.volume - volume).abs() <= 1e-6 * volume.max(1.0),
            "{op}: volume {} against {volume}",
            m.volume
        );
        assert!(m.valid, "{op}: check reports the result invalid");
        assert_eq!((m.faces, m.solids), (faces, solids), "{op}: faces, solids");
    }
}

fn cube(model: &mut Model, at: (f64, f64, f64), side: f64) -> ogeom::topo::Shape {
    let frame = Frame::new(Point::new(at.0, at.1, at.2), Direction::Z, Direction::X, T).unwrap();
    ogeom::algo::make_box(model, frame, (side, side, side), T)
        .unwrap()
        .shape
}

fn direction(x: f64, y: f64, z: f64) -> Direction {
    Direction::new(ogeom::math::Vector::new(x, y, z), T).unwrap()
}

/// Boxes touching along a whole edge, along part of one, and at a corner.
/// Touching shares no volume: the fuse holds both boxes (one solid on the
/// non-manifold edge the model permits where the edge is shared, two
/// solids where only a corner is), the common is empty, and the cut leaves
/// the first box as it was.
#[test]
fn boxes_touching_along_an_edge_or_at_a_corner_share_no_volume() {
    for at in [(10.0, 0.0, 10.0), (10.0, 5.0, 10.0)] {
        three_ways(
            &|model| (cube(model, (0.0, 0.0, 0.0), 10.0), cube(model, at, 10.0)),
            [(2000.0, 12, 1), (0.0, 0, 0), (1000.0, 6, 1)],
        );
    }
    three_ways(
        &|model| {
            (
                cube(model, (0.0, 0.0, 0.0), 10.0),
                cube(model, (10.0, 10.0, 10.0), 10.0),
            )
        },
        [(2000.0, 12, 2), (0.0, 0, 0), (1000.0, 6, 1)],
    );
}

/// Solids resting on a box's top face: a cylinder on its side touches it
/// along a line, a box stood on one corner and a cone stood on its apex
/// touch it at a point. None shares volume with the box.
#[test]
fn solids_resting_on_a_face_along_a_line_or_at_a_point_share_no_volume() {
    let pi = core::f64::consts::PI;
    // Radius 2, length 6, axis along x at height 12: the line y = 5, z = 10.
    three_ways(
        &|model| {
            let frame =
                Frame::new(Point::new(2.0, 5.0, 12.0), Direction::X, Direction::Y, T).unwrap();
            let roller = ogeom::algo::make_cylinder(model, frame, 2.0, 6.0, T)
                .unwrap()
                .shape;
            (cube(model, (0.0, 0.0, 0.0), 10.0), roller)
        },
        [(1000.0 + 24.0 * pi, 9, 2), (0.0, 0, 0), (1000.0, 6, 1)],
    );
    // A cube of side 2 whose three edges from one corner all climb at the
    // same angle: that corner is its lowest point, set on (5, 5, 10).
    let c = 1.0 / 3f64.sqrt();
    three_ways(
        &|model| {
            let frame = Frame::new(
                Point::new(5.0, 5.0, 10.0),
                direction(-1.0 / 6f64.sqrt(), -1.0 / 2f64.sqrt(), c),
                direction((2.0f64 / 3.0).sqrt(), 0.0, c),
                T,
            )
            .unwrap();
            let tilted = ogeom::algo::make_box(model, frame, (2.0, 2.0, 2.0), T)
                .unwrap()
                .shape;
            (cube(model, (0.0, 0.0, 0.0), 10.0), tilted)
        },
        [(1008.0, 12, 2), (0.0, 0, 0), (1000.0, 6, 1)],
    );
    // Base radius 3 at height 14, apex on (5, 5, 10).
    three_ways(
        &|model| {
            let frame = Frame::new(
                Point::new(5.0, 5.0, 14.0),
                direction(0.0, 0.0, -1.0),
                Direction::X,
                T,
            )
            .unwrap();
            let cone = ogeom::algo::make_cone(model, frame, 3.0, 0.0, 4.0, T)
                .unwrap()
                .shape;
            (cube(model, (0.0, 0.0, 0.0), 10.0), cone)
        },
        [(1000.0 + 12.0 * pi, 8, 2), (0.0, 0, 0), (1000.0, 6, 1)],
    );
}

/// A cylinder inscribed in a box touches its four walls along the lines
/// where its quarter turns stand, which are exactly the halves and quarters
/// of its wall's chart: every regular probe of the unsplit wall lands on a
/// contact line. Off those lines the wall reads plainly inside the box.
/// The cut leaves the four corners, pairwise touching along those lines;
/// with the caps flush, the box's top and bottom keep the disc as a hole
/// that touches their outline at four points.
#[test]
fn a_cylinder_inscribed_in_a_box_is_read_off_its_contact_lines() {
    let pi = core::f64::consts::PI;
    let disc = 25.0 * pi;
    for (low, height, fuse_volume, fuse_faces, cut_faces) in [
        (-2.0, 14.0, 1000.0 + 4.0 * disc, 16, 13),
        (0.0, 10.0, 1000.0, 8, 7),
    ] {
        three_ways(
            &|model| {
                let frame =
                    Frame::new(Point::new(5.0, 5.0, low), Direction::Z, Direction::X, T).unwrap();
                let drum = ogeom::algo::make_cylinder(model, frame, 5.0, height, T)
                    .unwrap()
                    .shape;
                (cube(model, (0.0, 0.0, 0.0), 10.0), drum)
            },
            [
                (fuse_volume, fuse_faces, 1),
                (10.0 * disc, 3, 1),
                (1000.0 - 10.0 * disc, cut_faces, 1),
            ],
        );
    }
}

/// A ball seated in a bore touches it along the equator and crosses it
/// nowhere. The boolean keeps that curve out of its arithmetic (a contact
/// carries no parity, so nothing is inside on one side of it), and the
/// section still reports it, because the curve is there.
#[test]
fn a_tangential_contact_is_sectioned_but_not_classified() {
    let mut model = Model::new();
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 2.0, T)
        .unwrap()
        .shape;
    let bore = ogeom::algo::make_cylinder(
        &mut model,
        Frame::new(Point::new(0.0, 0.0, -3.0), Direction::Z, Direction::X, T).unwrap(),
        2.0,
        6.0,
        T,
    )
    .unwrap()
    .shape;

    // Classification is untouched by the touch: the ball is inside the bore,
    // so their fuse is the bore and their common is the ball.
    let pi = core::f64::consts::PI;
    let fused = ogeom::boolean::fuse(&mut model, &ball, &bore, T)
        .unwrap()
        .shape;
    let cylinder_volume = pi * 4.0 * 6.0;
    assert!(
        (volume(&model, &fused) - cylinder_volume).abs() < cylinder_volume * 1e-3,
        "the ball adds nothing outside the bore: {}",
        volume(&model, &fused)
    );

    // The section is the contact circle: radius 2 in the plane z = 0.
    let cut = ogeom::boolean::section(&mut model, &ball, &bore, T)
        .unwrap()
        .shape;
    let edges = ogeom::topo::explore_unique(&model, &cut, ogeom::topo::ShapeType::Edge).unwrap();
    assert_eq!(edges.len(), 1, "one contact curve, once");
    let length = ogeom::algo::linear_properties(&model, &edges[0], fine(), T)
        .unwrap()
        .mass;
    let circle = 2.0 * pi * 2.0;
    assert!(
        (length - circle).abs() < circle * 1e-3,
        "the whole equator: {length}"
    );
}

/// A sphere's chart is bounded above and below by *poles* (edges that are
/// points in space) and by a seam it meets twice. Leave the poles out of
/// the arrangement and the chart has no top or bottom, so nothing can be
/// arranged inside it and every boolean over a ball fails, whether or not
/// the ball is anywhere near the other solid. They are in it.
#[test]
fn a_ball_is_boolean_material_like_anything_else() {
    let corner = Frame::new(Point::new(-5.0, -5.0, -5.0), Direction::Z, Direction::X, T).unwrap();
    let at = |p: Point| Frame::new(p, Direction::Z, Direction::X, T).unwrap();
    let pi = core::f64::consts::PI;
    let ball_volume = 4.0 / 3.0 * pi * 8.0;

    // Apart: the fuse is both, and the common is nothing.
    let mut model = Model::new();
    let far = ogeom::algo::make_sphere(&mut model, at(Point::new(20.0, 0.0, 0.0)), 2.0, T)
        .unwrap()
        .shape;
    let brick = ogeom::algo::make_box(&mut model, corner, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let apart = ogeom::boolean::fuse(&mut model, &far, &brick, T)
        .unwrap()
        .shape;
    assert!(
        (volume(&model, &apart) - (1000.0 + ball_volume)).abs() < 1.0,
        "both lumps, untouched: {}",
        volume(&model, &apart)
    );

    // Swallowed: the ball is inside, so the fuse is the brick and the cut
    // hollows a spherical void out of it.
    let mut model = Model::new();
    let inside = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 2.0, T)
        .unwrap()
        .shape;
    let brick = ogeom::algo::make_box(&mut model, corner, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let swallowed = ogeom::boolean::fuse(&mut model, &inside, &brick, T)
        .unwrap()
        .shape;
    assert!(
        (volume(&model, &swallowed) - 1000.0).abs() < 1.0,
        "the brick already held it: {}",
        volume(&model, &swallowed)
    );
    let hollow = ogeom::boolean::cut(&mut model, &brick, &inside, T)
        .unwrap()
        .shape;
    assert!(
        (volume(&model, &hollow) - (1000.0 - ball_volume)).abs() < 1.0,
        "a spherical void: {}",
        volume(&model, &hollow)
    );

    // Sitting on the lid: the section circle wraps the sphere's seam, and
    // the halves classify either side of it.
    let mut model = Model::new();
    let dome = ogeom::algo::make_sphere(&mut model, at(Point::new(0.0, 0.0, 5.0)), 2.0, T)
        .unwrap()
        .shape;
    let brick = ogeom::algo::make_box(&mut model, corner, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let capped = ogeom::boolean::fuse(&mut model, &dome, &brick, T)
        .unwrap()
        .shape;
    assert!(
        (volume(&model, &capped) - (1000.0 + ball_volume / 2.0)).abs() < 1.0,
        "the brick and the half that stands proud: {}",
        volume(&model, &capped)
    );
}

/// A ball resting exactly on a lid touches it at one point. A point bounds
/// no material (there is no side of it that is inside on one hand and
/// outside on the other), so the boolean carries the touch instead of
/// refusing it, and what comes back is both volumes joined at that point.
#[test]
fn a_point_touch_is_carried_rather_than_refused() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let ball = ogeom::algo::make_sphere(
        &mut model,
        Frame::new(Point::new(10.0, 10.0, 13.0), Direction::Z, Direction::X, T).unwrap(),
        3.0,
        T,
    )
    .unwrap()
    .shape;

    let pi = core::f64::consts::PI;
    let ball_volume = 4.0 / 3.0 * pi * 27.0;
    let fused = ogeom::boolean::fuse(&mut model, &block, &ball, T)
        .unwrap()
        .shape;
    assert!(
        (volume(&model, &fused) - (4000.0 + ball_volume)).abs() < 1.0,
        "both volumes, joined at the point they share: {}",
        volume(&model, &fused)
    );
    // And the cut takes nothing: the ball meets the block in a point, and a
    // point has no volume to remove.
    let cut = ogeom::boolean::cut(&mut model, &block, &ball, T)
        .unwrap()
        .shape;
    assert!(
        (volume(&model, &cut) - 4000.0).abs() < 1e-6,
        "a point removes nothing: {}",
        volume(&model, &cut)
    );
}

#[test]
fn a_fitted_edge_on_a_shared_cylinder_still_melts_the_same_domain_contact() {
    // Two drums on the *identical* cylinder chart, one of them scooped at
    // the top by a crossing cylinder: its wall's upper edges are marched,
    // fitted curves with fitted pcurves, which no closed-form projection can
    // carry into the other wall's chart. On the identical chart the stored
    // pcurve already is the projection, so the same-domain melt takes it in
    // place of one, and the fuse of a contained solid comes out as the
    // container.
    let mut model = Model::new();

    let tall = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 6.0, T)
        .unwrap()
        .shape;
    let scoop_frame = Frame::new(
        Point::new(-5.0, 0.0, 6.7),
        ogeom::math::Direction::X,
        ogeom::math::Direction::Y,
        T,
    )
    .unwrap();
    let scoop = ogeom::algo::make_cylinder(&mut model, scoop_frame, 2.5, 10.0, T)
        .unwrap()
        .shape;
    let wavy = ogeom::boolean::cut(&mut model, &tall, &scoop, T)
        .unwrap()
        .shape;
    let short = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 4.0, T)
        .unwrap()
        .shape;

    let before = volume(&model, &wavy);
    let fused = ogeom::boolean::fuse(&mut model, &short, &wavy, T).unwrap();
    assert!(
        ogeom::algo::check(&model, &fused.shape, T)
            .unwrap()
            .is_valid(),
        "the fused drum is a valid solid"
    );
    // The container integrates on its exact surfaces. The fused drum's
    // fitted trims leave it to the mesh, measured a decade finer so its
    // deficit stays inside the budget. The budget is the mesh, not the
    // melt: the fuse either resolves the contact or refuses by name.
    let finer = Deflection::with_chord(1e-4).unwrap();
    let measured = ogeom::algo::volume_properties(&model, &fused.shape, finer, T)
        .unwrap()
        .mass;
    assert!(
        (measured - before).abs() < 2e-2,
        "fusing a contained drum should give the container: {measured} vs {before}"
    );
    // Every face of the contained drum melted rather than surviving as a
    // skin. The wall stays split where the contained drum's rim touched it:
    // one face more than the container had, none of them inside.
    assert_eq!(
        explore_unique(&model, &fused.shape, ShapeType::Face)
            .unwrap()
            .len(),
        explore_unique(&model, &wavy, ShapeType::Face)
            .unwrap()
            .len()
            + 1
    );
}

/// A scaled copy shares whole planes with its original, and the shared
/// regions nest rather than match: the small box's faces at the origin lie
/// strictly inside the big box's. The contact is real same-domain contact,
/// so it has to reach the melt, which it only does while the scale is
/// carried as the placement it is, since a restated plane does not say it
/// is one and no closed form recognizes the pair.
#[test]
fn a_box_and_its_doubled_copy_fuse_into_the_bigger_box() {
    let mut model = Model::with_tolerances(T);
    let a = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let doubled = ogeom::math::GeneralTransform {
        linear: ogeom::math::Matrix3::from_columns(
            ogeom::math::Vector::new(2.0, 0.0, 0.0),
            ogeom::math::Vector::new(0.0, 2.0, 0.0),
            ogeom::math::Vector::new(0.0, 0.0, 2.0),
        ),
        translation: ogeom::math::Vector::ZERO,
    };
    let b = ogeom::algo::general_transformed_shape(&mut model, &a, &doubled, T)
        .unwrap()
        .shape;

    let fused = ogeom::boolean::fuse(&mut model, &a, &b, T).unwrap();

    let diagnosis = ogeom::algo::check(&model, &fused.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    // The small box is inside the big one, so the union is the big one.
    let measured = volume(&model, &fused.shape);
    assert!(
        (measured - 8000.0).abs() < 1e-6,
        "fused volume {measured} against 8000"
    );
}

/// The same pair slid along the plane they share, so neither shared face
/// contains the other and the overlap is partial in both directions.
#[test]
fn a_doubled_copy_slid_along_its_shared_plane_fuses_over_a_partial_contact() {
    let build = |model: &mut Model| {
        let a = ogeom::algo::make_box(model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let slid = ogeom::math::GeneralTransform {
            linear: ogeom::math::Matrix3::from_columns(
                ogeom::math::Vector::new(2.0, 0.0, 0.0),
                ogeom::math::Vector::new(0.0, 2.0, 0.0),
                ogeom::math::Vector::new(0.0, 0.0, 2.0),
            ),
            translation: ogeom::math::Vector::new(5.0, 5.0, 0.0),
        };
        let b = ogeom::algo::general_transformed_shape(model, &a, &slid, T)
            .unwrap()
            .shape;
        (a, b)
    };

    // a = [0,10]³, b = [5,25]×[5,25]×[0,20]; they share the z = 0 plane and
    // overlap there on [5,10]², which neither face contains.
    let mut model = Model::with_tolerances(T);
    let (a, b) = build(&mut model);
    let fused = ogeom::boolean::fuse(&mut model, &a, &b, T).unwrap();
    assert!(
        ogeom::algo::check(&model, &fused.shape, T)
            .unwrap()
            .is_valid(),
        "the fused body is not valid"
    );
    let measured = volume(&model, &fused.shape);
    assert!(
        (measured - 8750.0).abs() < 1e-6,
        "fused volume {measured} against 8750"
    );

    let mut model = Model::with_tolerances(T);
    let (a, b) = build(&mut model);
    let shared = ogeom::boolean::common(&mut model, &a, &b, T).unwrap();
    let measured = volume(&model, &shared.shape);
    assert!(
        (measured - 250.0).abs() < 1e-6,
        "common volume {measured} against 250"
    );

    let mut model = Model::with_tolerances(T);
    let (a, b) = build(&mut model);
    let rest = ogeom::boolean::cut(&mut model, &a, &b, T).unwrap();
    let measured = volume(&model, &rest.shape);
    assert!(
        (measured - 750.0).abs() < 1e-6,
        "cut volume {measured} against 750"
    );
}

/// A cylinder seated on the face it pierces: the two solids share the plane
/// they both stand on, and the cylinder's cap lies strictly inside the
/// block's bottom face. Both arguments therefore describe that one disk, and
/// exactly one of the two descriptions may survive: a question `cut` never
/// has to ask and `common` must.
#[test]
fn a_cylinder_seated_on_the_face_it_pierces_shares_only_the_segment_between_them() {
    let build = |model: &mut Model| {
        let frame = Frame::new(Point::new(10.0, 10.0, 0.0), Direction::Z, Direction::X, T).unwrap();
        let circle = ogeom::math::Circle::new(frame, 4.0, T).unwrap();
        let curve = ogeom::geom::Curve::Circle(ogeom::geom::CircleCurve::new(circle));
        let range = <ogeom::geom::Curve as ogeom::geom::Curve3d>::domain(&curve);
        let edge = ogeom::algo::make_edge(model, curve, range, T)
            .unwrap()
            .shape;
        let round =
            ogeom::geom::PlaneSurface::over(ogeom::math::Plane::XY, (5.0, 15.0), (5.0, 15.0))
                .unwrap();
        let cap = ogeom::algo::make_face_with_pcurves(model, round.into(), &[vec![edge]], T)
            .unwrap()
            .shape;
        let cylinder =
            ogeom::algo::make_prism(model, &cap, ogeom::math::Vector::new(0.0, 0.0, 20.0), T)
                .unwrap()
                .shape;

        let corners = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(20.0, 0.0, 0.0),
            Point::new(20.0, 20.0, 0.0),
            Point::new(0.0, 20.0, 0.0),
        ];
        let wire = ogeom::algo::make_polygon(model, &corners, true, T)
            .unwrap()
            .shape;
        let flat =
            ogeom::geom::PlaneSurface::over(ogeom::math::Plane::XY, (-1.0, 21.0), (-1.0, 21.0))
                .unwrap();
        let edges = model.children_of(&wire).unwrap();
        let base = ogeom::algo::make_face_with_pcurves(model, flat.into(), &[edges], T)
            .unwrap()
            .shape;
        let block =
            ogeom::algo::make_prism(model, &base, ogeom::math::Vector::new(0.0, 0.0, 10.0), T)
                .unwrap()
                .shape;
        (block, cylinder)
    };
    let pi = core::f64::consts::PI;

    let mut model = Model::with_tolerances(T);
    let (block, cylinder) = build(&mut model);
    let shared = ogeom::boolean::common(&mut model, &block, &cylinder, T).unwrap();
    assert!(
        ogeom::algo::check(&model, &shared.shape, T)
            .unwrap()
            .is_valid(),
        "the shared post is not valid"
    );
    // The post the block's height cuts out of the cylinder: two disks and the
    // wall between them. Three faces. A fourth would be the shared disk
    // described twice.
    assert_eq!(
        explore_unique(&model, &shared.shape, ShapeType::Face)
            .unwrap()
            .len(),
        3,
        "the disk the two solids share is described more than once"
    );
    let expected = pi * 16.0 * 10.0;
    let measured = volume(&model, &shared.shape);
    assert!(
        (measured - expected).abs() < expected * 1e-3,
        "common volume {measured} against {expected}"
    );

    // The cut of the same pair, which never describes the disk twice.
    let mut model = Model::with_tolerances(T);
    let (block, cylinder) = build(&mut model);
    let bored = ogeom::boolean::cut(&mut model, &block, &cylinder, T).unwrap();
    let expected = 20.0_f64.mul_add(20.0 * 10.0, -(pi * 16.0 * 10.0));
    let measured = volume(&model, &bored.shape);
    assert!(
        (measured - expected).abs() < expected * 1e-3,
        "cut volume {measured} against {expected}"
    );
}

/// A shear is the transform a placement cannot express, so the body is
/// restated as patches: there is no other way to carry a box's planes
/// under one. What the patches lose is the *word* plane, and
/// coincidence is decided on what the geometry says: the two faces on the
/// z = 0 plane are one surface, and each's edges are carried into the
/// other's chart over the stretch that lies on its window. The union, the
/// difference and the intersection come out at their exact volumes.
#[test]
fn a_sheared_copy_sharing_a_plane_combines_to_its_exact_volumes() {
    let mut model = Model::with_tolerances(T);
    let a = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    // Unit determinant, so the copy keeps its volume, and slid along the
    // z = 0 plane the two go on sharing, overlapping by a wedge.
    let shear = ogeom::math::GeneralTransform {
        linear: ogeom::math::Matrix3 {
            rows: [[1.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        },
        translation: ogeom::math::Vector::new(5.0, 0.0, 0.0),
    };
    let b = ogeom::algo::general_transformed_shape(&mut model, &a, &shear, T)
        .unwrap()
        .shape;

    // The overlap: over y in 0..10 the copy starts at x = 5 + y / 2, so the
    // shared section is 25 and the wedge 250.
    for (name, op, expected) in [
        (
            "fuse",
            ogeom::boolean::fuse
                as fn(&mut Model, &ogeom::topo::Shape, &ogeom::topo::Shape, Tolerances) -> _,
            1750.0,
        ),
        ("cut", ogeom::boolean::cut, 750.0),
        ("common", ogeom::boolean::common, 250.0),
    ] {
        let out = op(&mut model, &a, &b, T).unwrap_or_else(|e| panic!("{name}: {e}"));
        let diagnosis = ogeom::algo::check(&model, &out.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
        let measured = volume(&model, &out.shape);
        assert!(
            (measured - expected).abs() < expected * 1e-6,
            "{name} volume {measured} against {expected}"
        );
    }
}

/// A sheared bar laid across a sheared drum on the plane both stand on:
/// the drum's rim leaves the bar's floor window and comes back onto it,
/// so the rim is carried into the bar's chart as one stretch for each
/// time it lies on the window, and the bar's floor is split by every one.
/// The rim's seam sits once on the window (three stretches, two of them
/// ending at the seam) and once off it (two stretches, both ends of each
/// on the window's rim). The shear has unit determinant, so each result
/// keeps the volume of the unsheared bar and drum: the bar `|y| <= 1`
/// meets the disk of radius 5 in a band of area `2·sqrt(24) + 50·asin(1/5)`.
#[test]
#[ignore = "heavy"]
fn a_rim_weaving_across_a_sheared_bar_s_window_combines_to_its_exact_volumes() {
    let shear = ogeom::math::GeneralTransform {
        linear: ogeom::math::Matrix3 {
            rows: [[1.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        },
        translation: ogeom::math::Vector::new(0.0, 0.0, 0.0),
    };
    let build = |model: &mut Model, seam: Direction| {
        let corner =
            Frame::new(Point::new(-10.0, -1.0, 0.0), Direction::Z, Direction::X, T).unwrap();
        let bar = ogeom::algo::make_box(model, corner, (20.0, 2.0, 2.0), T)
            .unwrap()
            .shape;
        let axis = Frame::new(Point::ORIGIN, Direction::Z, seam, T).unwrap();
        let drum = ogeom::algo::make_cylinder(model, axis, 5.0, 3.0, T)
            .unwrap()
            .shape;
        let bar = ogeom::algo::general_transformed_shape(model, &bar, &shear, T)
            .unwrap()
            .shape;
        let drum = ogeom::algo::general_transformed_shape(model, &drum, &shear, T)
            .unwrap()
            .shape;
        (bar, drum)
    };
    let band = 2.0f64.mul_add(24.0f64.sqrt(), 50.0 * 0.2f64.asin());
    let (bar, drum) = (80.0, core::f64::consts::PI * 25.0 * 3.0);
    let shared = band * 2.0;
    for seam in [Direction::X, Direction::Y] {
        for (name, op, expected) in [
            (
                "fuse",
                ogeom::boolean::fuse
                    as fn(&mut Model, &ogeom::topo::Shape, &ogeom::topo::Shape, Tolerances) -> _,
                bar + drum - shared,
            ),
            ("cut", ogeom::boolean::cut, bar - shared),
            ("common", ogeom::boolean::common, shared),
        ] {
            let mut model = Model::with_tolerances(T);
            let (a, b) = build(&mut model, seam);
            let out = op(&mut model, &a, &b, T).unwrap_or_else(|e| panic!("{name}: {e}"));
            let diagnosis = ogeom::algo::check(&model, &out.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{name}, seam {seam:?}: {diagnosis}");
            let measured = volume(&model, &out.shape);
            assert!(
                (measured - expected).abs() < expected * 1e-6,
                "{name}, seam {seam:?}: volume {measured} against {expected}"
            );
        }
    }
}

/// A box cut from an L-bracket flush with the bracket's wall: the box's
/// wall-side face lies on the wall's plane below the wall, and its edge on
/// the end face runs along the line of that face's own edge up the wall,
/// a length below it. The contact carried onto the end face is clipped to
/// the edge it runs along, not read over the whole line: read over the
/// line, the strip's side is never paved and the end face keeps the
/// strip. Three placements against the wall and one short of it, each to
/// its exact volume.
#[test]
fn a_box_cut_flush_with_a_bracket_s_wall_paves_the_end_face() {
    let bracket = |model: &mut Model| {
        let block = ogeom::algo::make_box(model, Frame::WORLD, (2.0, 2.0, 2.0), T)
            .unwrap()
            .shape;
        let seat = Frame::new(Point::new(1.0, -0.5, 1.0), Direction::Z, Direction::X, T).unwrap();
        let notch = ogeom::algo::make_box(model, seat, (2.0, 3.0, 2.0), T)
            .unwrap()
            .shape;
        ogeom::boolean::cut(model, &block, &notch, T).unwrap().shape
    };
    for (label, low, high) in [
        (
            "on the wall, the end and the leg's top",
            [1.0, 0.0, 0.5],
            [2.0, 0.5, 1.0],
        ),
        (
            "on the wall and the end, below the top",
            [1.0, 0.0, 0.3],
            [2.0, 0.5, 0.8],
        ),
        ("through the wall", [0.8, 0.0, 0.5], [2.0, 0.5, 1.0]),
        ("short of the wall", [1.2, 0.0, 0.5], [2.0, 0.5, 1.0]),
    ] {
        let mut model = Model::with_tolerances(T);
        let l = bracket(&mut model);
        let seat = Frame::new(
            Point::new(low[0], low[1], low[2]),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        let size = (high[0] - low[0], high[1] - low[1], high[2] - low[2]);
        let tool = ogeom::algo::make_box(&mut model, seat, size, T)
            .unwrap()
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &l, &tool, T)
            .unwrap_or_else(|e| panic!("{label}: {e}"))
            .shape;
        let diagnosis = ogeom::algo::check(&model, &cut, T).unwrap();
        assert!(diagnosis.is_valid(), "{label}: {:?}", diagnosis.problems);
        // Every placement lies wholly within the bracket (through the
        // wall too, the wall's material starting where the leg's top ends),
        // so the tool takes exactly its own volume.
        let want = 6.0 - size.0 * size.1 * size.2;
        let got = volume(&model, &cut);
        assert!((got - want).abs() < 1e-3, "{label}: {got} against {want}");
    }
}

/// A countersink whose bore is exactly its cone's narrow end: the cone's
/// rim lies in the bore's wall, so the one circle is both an edge the cone
/// keeps and a section the bore's wall is split along. Cut in either order
/// or as one tool, the plate comes out valid and short by exactly the
/// frustum and the bore below it.
#[test]
fn a_countersink_whose_bore_is_its_narrow_end_cuts_clean() {
    let mut model = Model::new();
    let at = |z: f64| Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X, T).unwrap();
    let corner = Frame::new(Point::new(-10.0, -10.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let plate = ogeom::algo::make_box(&mut model, corner, (20.0, 20.0, 5.0), T)
        .unwrap()
        .shape;
    let sink = ogeom::algo::make_cone(&mut model, at(1.0), 1.0, 5.0, 4.0, T)
        .unwrap()
        .shape;
    let bore = ogeom::algo::make_cylinder(&mut model, at(-1.0), 1.0, 7.0, T)
        .unwrap()
        .shape;
    let pi = core::f64::consts::PI;
    let want = 2000.0 - pi * 4.0 / 3.0 * (1.0 + 5.0 + 25.0) - pi;

    let sunk = ogeom::boolean::cut(&mut model, &plate, &sink, T)
        .unwrap()
        .shape;
    let one = ogeom::boolean::cut(&mut model, &sunk, &bore, T)
        .unwrap()
        .shape;
    let bored = ogeom::boolean::cut(&mut model, &plate, &bore, T)
        .unwrap()
        .shape;
    let other = ogeom::boolean::cut(&mut model, &bored, &sink, T)
        .unwrap()
        .shape;
    let tool = ogeom::boolean::fuse(&mut model, &sink, &bore, T)
        .unwrap()
        .shape;
    let whole = ogeom::boolean::cut(&mut model, &plate, &tool, T)
        .unwrap()
        .shape;
    for result in [&one, &other, &whole] {
        assert!(ogeom::algo::check(&model, result, T).unwrap().is_valid());
        let got = volume(&model, result);
        assert!((got - want).abs() < want * 1e-9, "{got} against {want}");
    }
}

/// A plate with a through slot chamfered at both openings, and a prism of
/// its top face pushed down through it. The prism's slot walls run through
/// the chamfers' outer rims, and at the bottom rim two of the plate's faces
/// (the underside and the chamfer) meet the wall along one line: the wall
/// is split there once, whichever face's section says so. Short of the
/// bottom, at it and through it, the three operations come out exact.
#[test]
fn a_prism_through_a_chamfered_slot_splits_its_walls_once() {
    let at = |x: f64, y: f64, z: f64| {
        Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
    };
    let mut model = Model::with_tolerances(T);
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 5.0), T)
        .unwrap()
        .shape;
    let slot = ogeom::algo::make_box(&mut model, at(3.0, 4.0, -1.0), (4.0, 2.0, 7.0), T)
        .unwrap()
        .shape;
    let plate = ogeom::boolean::cut(&mut model, &block, &slot, T)
        .unwrap()
        .shape;
    let rims: Vec<_> = explore_unique(&model, &plate, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let b = ogeom::algo::tight_bounds(&model, e, T).unwrap();
            let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
            (hi.z - lo.z).abs() < 1e-9
                && (lo.z.abs() < 1e-9 || (lo.z - 5.0).abs() < 1e-9)
                && lo.x > 2.9
                && hi.x < 7.1
                && lo.y > 3.9
                && hi.y < 6.1
        })
        .collect();
    assert_eq!(rims.len(), 8);
    let part = ogeom::fillet::chamfer_edges(&mut model, &plate, &rims, 0.8, T)
        .unwrap()
        .shape;
    let top = explore_unique(&model, &part, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let (p, n) = ogeom::algo::face_normal(&model, f, T).unwrap();
            n.z > 0.999 && (p.z - 5.0).abs() < 1e-9
        })
        .unwrap();
    let v_part = volume(&model, &part);
    for depth in [3.0, 5.0, 10.0] {
        let pad = ogeom::algo::make_prism(
            &mut model,
            &top,
            ogeom::math::Vector::new(0.0, 0.0, -depth),
            T,
        )
        .unwrap()
        .shape;
        let v_pad = volume(&model, &pad);
        // Inside the part: the top face's area down to the bottom at most.
        let v_common = v_pad * depth.min(5.0) / depth;
        for (name, op, expected) in [
            (
                "fuse",
                ogeom::boolean::fuse
                    as fn(&mut Model, &ogeom::topo::Shape, &ogeom::topo::Shape, Tolerances) -> _,
                v_part + v_pad - v_common,
            ),
            ("cut", ogeom::boolean::cut, v_part - v_common),
            ("common", ogeom::boolean::common, v_common),
        ] {
            let out = op(&mut model, &part, &pad, T)
                .unwrap_or_else(|e| panic!("{name} at depth {depth}: {e}"));
            let diagnosis = ogeom::algo::check(&model, &out.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{name} at depth {depth}: {diagnosis}");
            let measured = volume(&model, &out.shape);
            assert!(
                (measured - expected).abs() < expected.max(1.0) * 1e-6,
                "{name} at depth {depth}: {measured} against {expected}"
            );
        }
    }
}

/// A block whose one wall leans a few microns off the vertical, and a
/// prism of its top face pushed down into it and through it: the prism's
/// side meets the wall along their shared top edge and parts from it by the
/// lean. The side's edges run beside the wall's, not along them, and the
/// three operations come out sound at a lean of a micron, ten, and a
/// hundred, either way.
#[test]
fn a_prism_beside_a_leaning_wall_parts_from_it() {
    let p = Point::new;
    for delta in [1e-6, -1e-6, 1e-5, -1e-5, 1e-4, -1e-4] {
        let mut model = Model::with_tolerances(T);
        let corners = [
            p(0.0, 0.0, 0.0),
            p(10.0 + delta, 0.0, 0.0),
            p(10.0 + delta, 10.0, 0.0),
            p(0.0, 10.0, 0.0),
            p(0.0, 0.0, 5.0),
            p(10.0, 0.0, 5.0),
            p(10.0, 10.0, 5.0),
            p(0.0, 10.0, 5.0),
        ];
        let part = ogeom::algo::make_hexahedron(&mut model, corners, T)
            .unwrap()
            .shape;
        let top = explore_unique(&model, &part, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                let (q, n) = ogeom::algo::face_normal(&model, f, T).unwrap();
                n.z > 0.999 && (q.z - 5.0).abs() < 1e-9
            })
            .unwrap();
        for depth in [3.0, 10.0] {
            let pad = ogeom::algo::make_prism(
                &mut model,
                &top,
                ogeom::math::Vector::new(0.0, 0.0, -depth),
                T,
            )
            .unwrap()
            .shape;
            let (v_part, v_pad) = (volume(&model, &part), volume(&model, &pad));
            let mut got = Vec::new();
            for (name, op) in [
                (
                    "fuse",
                    ogeom::boolean::fuse
                        as fn(
                            &mut Model,
                            &ogeom::topo::Shape,
                            &ogeom::topo::Shape,
                            Tolerances,
                        ) -> _,
                ),
                ("cut", ogeom::boolean::cut),
                ("common", ogeom::boolean::common),
            ] {
                let out = op(&mut model, &part, &pad, T)
                    .unwrap_or_else(|e| panic!("{name}, lean {delta}, depth {depth}: {e}"));
                let diagnosis = ogeom::algo::check(&model, &out.shape, T).unwrap();
                assert!(
                    diagnosis.is_valid(),
                    "{name}, lean {delta}, depth {depth}: {diagnosis}"
                );
                got.push(volume(&model, &out.shape));
            }
            let [fused, cut, common] = got[..] else {
                unreachable!()
            };
            let scale = v_part + v_pad;
            assert!(
                (fused + common - scale).abs() < 1e-6 * scale,
                "lean {delta}, depth {depth}: fuse {fused} + common {common} against {scale}"
            );
            assert!(
                (cut + common - v_part).abs() < 1e-6 * scale,
                "lean {delta}, depth {depth}: cut {cut} + common {common} against {v_part}"
            );
        }
    }
}

/// A block with rounded upright edges, meshed, its points moved a couple
/// of microns and stored in single precision a hundred millimetres out (as
/// an exported mesh holds them), and converted as it is: its rounds are
/// strips of facets, each leaning its own few microns. A prism of its top
/// face pushed down into it and through it runs along those facets, and
/// at a fuzz wider than the lean the three operations come out sound and
/// hold what the two inputs do.
#[test]
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    reason = "the rounding to single precision is the point"
)]
fn a_prism_along_a_faceted_round_combines_at_a_fuzz() {
    let mut source = Model::with_tolerances(T);
    let block = ogeom::algo::make_box(&mut source, Frame::WORLD, (2.0, 14.0, 8.5), T)
        .unwrap()
        .shape;
    let upright: Vec<_> = explore_unique(&source, &block, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let b = ogeom::algo::tight_bounds(&source, e, T).unwrap();
            b.high().unwrap().z - b.low().unwrap().z > 1.0
        })
        .collect();
    let rounded = ogeom::fillet::fillet_edges(&mut source, &block, &upright, 0.5, T)
        .unwrap()
        .shape;
    let mut mesh =
        ogeom::mesh::triangulate(&source, &rounded, Deflection::with_chord(0.05).unwrap(), T)
            .unwrap();
    let mut state: u64 = 1;
    let mut nudge = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((state >> 11) as f64 / (1u64 << 53) as f64).mul_add(2.0, -1.0) * 2e-6
    };
    for p in &mut mesh.positions {
        let (dx, dy, dz) = (nudge(), nudge(), nudge());
        *p = Point::new(
            f64::from((p.x + 100.0 + dx) as f32),
            f64::from((p.y + 100.0 + dy) as f32),
            f64::from((p.z + dz) as f32),
        );
    }
    let mut model = Model::with_tolerances(T);
    let options = ogeom::algo::MeshSolidOptions {
        recognize: false,
        keep_vertices: true,
        quantum: Some(ogeom::algo::single_precision_quantum(&mesh)),
        ..ogeom::algo::MeshSolidOptions::default()
    };
    let part = ogeom::algo::solid_from_mesh(&mut model, &mesh, &options, T)
        .unwrap()
        .shape;
    let top = explore_unique(&model, &part, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let (q, n) = ogeom::algo::face_normal(&model, f, T).unwrap();
            n.z > 0.999_999 && (q.z - 8.5).abs() < 1e-5
        })
        .unwrap();
    let fuzz = 1e-5;
    let loose = Tolerances::with_scale(1e-7 / fuzz).unwrap();
    for depth in [3.0, 10.0] {
        let pad = ogeom::algo::make_prism(
            &mut model,
            &top,
            ogeom::math::Vector::new(0.0, 0.0, -depth),
            T,
        )
        .unwrap()
        .shape;
        let (v_part, v_pad) = (volume(&model, &part), volume(&model, &pad));
        let fused = ogeom::boolean::fuse_fuzzy(&mut model, &part, &pad, fuzz, T)
            .unwrap_or_else(|e| panic!("fuse at depth {depth}: {e}"));
        let cut = ogeom::boolean::cut_fuzzy(&mut model, &part, &pad, fuzz, T)
            .unwrap_or_else(|e| panic!("cut at depth {depth}: {e}"));
        let common = ogeom::boolean::common(&mut model, &part, &pad, loose)
            .unwrap_or_else(|e| panic!("common at depth {depth}: {e}"));
        let mut got = Vec::new();
        for (name, out) in [("fuse", &fused), ("cut", &cut), ("common", &common)] {
            let diagnosis = ogeom::algo::check(&model, &out.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{name} at depth {depth}: {diagnosis}");
            got.push(volume(&model, &out.shape));
        }
        let scale = v_part + v_pad;
        assert!(
            (got[0] + got[2] - scale).abs() < 1e-5 * scale,
            "depth {depth}: fuse {} + common {} against {scale}",
            got[0],
            got[2]
        );
        assert!(
            (got[1] + got[2] - v_part).abs() < 1e-5 * scale,
            "depth {depth}: cut {} + common {} against {v_part}",
            got[1],
            got[2]
        );
    }
}

/// A drum two microns wider than the box it stands in pokes out of each
/// wall along a strip too narrow for a mesh's samples to land on. Whatever
/// the boolean makes of it, it does not answer as if one lay within the
/// other.
#[test]
fn a_solid_poking_microns_through_another_is_not_taken_as_nested() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (4.0, 4.0, 4.0), T)
        .unwrap()
        .shape;
    let turn: f64 = 0.05;
    let seat = Frame::new(
        Point::new(2.0, 2.0, 1.0),
        Direction::Z,
        Direction::new(ogeom::math::Vector::new(turn.cos(), turn.sin(), 0.0), T).unwrap(),
        T,
    )
    .unwrap();
    let drum = ogeom::algo::make_cylinder(&mut model, seat, 2.002, 2.0, T)
        .unwrap()
        .shape;
    if let Ok(fused) = ogeom::boolean::fuse(&mut model, &block, &drum, T) {
        assert!(!fused.shape.is_same(&block), "the union is not the box");
    }
    if let Ok(common) = ogeom::boolean::common(&mut model, &block, &drum, T) {
        assert!(
            !common.shape.is_same(&drum),
            "the intersection is not the drum"
        );
    }
}

/// Two boxes stacked a few microns apart, closer than the boolean welds:
/// their facing walls are not kept as a wall of no thickness inside one
/// shell. Below the weld the union is the one tall box, as touching boxes
/// fuse. Past it, the two boxes apart.
#[test]
fn boxes_stacked_microns_apart_fuse_without_a_membrane() {
    for gap in [1e-6, 5e-6, 1e-5, 2e-5] {
        let mut model = Model::new();
        let low = ogeom::algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T)
            .unwrap()
            .shape;
        let seat = Frame::new(
            Point::new(0.0, 0.0, 2.0 + gap),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        let high = ogeom::algo::make_box(&mut model, seat, (2.0, 2.0, 2.0), T)
            .unwrap()
            .shape;
        let fused = ogeom::boolean::fuse(&mut model, &low, &high, T)
            .unwrap()
            .shape;
        assert!(ogeom::algo::check(&model, &fused, T).unwrap().is_valid());
        let faces = explore_unique(&model, &fused, ShapeType::Face).unwrap();
        let solids = explore_unique(&model, &fused, ShapeType::Solid).unwrap();
        // Every edge of every shell bounds two faces.
        for shell in explore_unique(&model, &fused, ShapeType::Shell).unwrap() {
            let shell_faces = explore_unique(&model, &shell, ShapeType::Face).unwrap();
            for edge in explore_unique(&model, &shell, ShapeType::Edge).unwrap() {
                let users = shell_faces
                    .iter()
                    .filter(|f| {
                        explore_unique(&model, f, ShapeType::Edge)
                            .unwrap()
                            .iter()
                            .any(|e| e.node() == edge.node())
                    })
                    .count();
                assert_eq!(users, 2, "gap {gap}");
            }
        }
        match solids.len() {
            1 => assert_eq!(faces.len(), 10, "gap {gap}"),
            2 => assert_eq!(faces.len(), 12, "gap {gap}"),
            n => panic!("gap {gap}: {n} solids"),
        }
    }
}

/// A solid cut by itself, or by a copy described another way, is cut away
/// entirely.
#[test]
fn a_solid_cut_by_its_own_copy_leaves_nothing() {
    let mut model = Model::new();
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 1.0, T)
        .unwrap()
        .shape;
    let same = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 1.0, T)
        .unwrap()
        .shape;
    let turned = Frame::new(Point::ORIGIN, Direction::X, Direction::Y, T).unwrap();
    let other = ogeom::algo::make_sphere(&mut model, turned, 1.0, T)
        .unwrap()
        .shape;
    for tool in [same, other] {
        let cut = ogeom::boolean::cut(&mut model, &ball, &tool, T).unwrap();
        assert!(
            explore_unique(&model, &cut.shape, ShapeType::Face)
                .unwrap()
                .is_empty()
        );
    }
}

/// A drum a hair wider than the box it stands in pokes out of every side
/// by a sliver. Its caps' rims bow past their sampled chords by far more
/// than that, and whichever way the drum's seam is turned the caps meet the
/// box sides and every operation keeps its volumes.
#[test]
fn a_drum_poking_a_sliver_out_of_a_box_works_at_any_turn() {
    let mut reference: Option<(f64, f64, f64)> = None;
    for turn in [0.0_f64, 0.05, 0.123, 0.3] {
        let mut model = Model::new();
        let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (4.0, 4.0, 4.0), T)
            .unwrap()
            .shape;
        let frame = Frame::new(
            Point::new(2.0, 2.0, 1.0),
            Direction::Z,
            Direction::new(ogeom::math::Vector::new(turn.cos(), turn.sin(), 0.0), T).unwrap(),
            T,
        )
        .unwrap();
        let drum = ogeom::algo::make_cylinder(&mut model, frame, 2.002, 2.0, T)
            .unwrap()
            .shape;
        let (va, vb) = (volume(&model, &block), volume(&model, &drum));
        let fuse = ogeom::boolean::fuse(&mut model, &block, &drum, T).unwrap();
        let common = ogeom::boolean::common(&mut model, &block, &drum, T).unwrap();
        let cut = ogeom::boolean::cut(&mut model, &block, &drum, T).unwrap();
        let (f, c, k) = (
            volume(&model, &fuse.shape),
            volume(&model, &common.shape),
            volume(&model, &cut.shape),
        );
        assert!((f + c - va - vb).abs() < 1e-6, "turn {turn}");
        assert!((k + c - va).abs() < 1e-6, "turn {turn}");
        match reference {
            None => reference = Some((f, c, k)),
            Some((f0, c0, k0)) => {
                assert!((f - f0).abs() < 1e-6 && (c - c0).abs() < 1e-6 && (k - k0).abs() < 1e-6);
            }
        }
    }
}

/// Two drums a thousandth to a hundred-thousandth apart, in radius or in
/// axis, overlapping half their height: the wider one's cap leaves an
/// annulus or a crescent that thin, narrower than its outline's chords bow,
/// whose ends meet the rim at a shallow angle, and each operation keeps its
/// volumes to the closed forms.
#[test]
fn drums_a_thousandth_apart_keep_their_volumes() {
    let at =
        |x: f64, z: f64| Frame::new(Point::new(x, 0.0, z), Direction::Z, Direction::X, T).unwrap();
    for (radius, offset) in [
        (2.001, 0.0),
        (2.0001, 0.0),
        (2.00003, 0.0),
        (2.00002, 0.0),
        (2.0, 1e-3),
        (2.0, 1e-4),
    ] {
        let mut model = Model::new();
        let low = ogeom::algo::make_cylinder(&mut model, at(0.0, 0.0), 2.0, 4.0, T)
            .unwrap()
            .shape;
        let high = ogeom::algo::make_cylinder(&mut model, at(offset, 2.0), radius, 4.0, T)
            .unwrap()
            .shape;
        let (va, vb) = (volume(&model, &low), volume(&model, &high));
        let fuse = ogeom::boolean::fuse(&mut model, &low, &high, T).unwrap();
        let common = ogeom::boolean::common(&mut model, &low, &high, T).unwrap();
        let cut = ogeom::boolean::cut(&mut model, &low, &high, T).unwrap();
        let (f, c, k) = (
            volume(&model, &fuse.shape),
            volume(&model, &common.shape),
            volume(&model, &cut.shape),
        );
        assert!(
            (f + c - va - vb).abs() < 1e-6,
            "{radius} {offset}: {f} + {c}"
        );
        assert!((k + c - va).abs() < 1e-6, "{radius} {offset}: {k} + {c}");
        // The common part is the narrower drum's upper half, less the
        // crescent the offset leaves outside the other.
        let pi = core::f64::consts::PI;
        let lens = 8.0 * (offset / 4.0).acos() - offset / 2.0 * (16.0 - offset * offset).sqrt();
        assert!(
            (c - 2.0 * lens.min(4.0 * pi)).abs() < 1e-6,
            "{radius} {offset}: {c}"
        );
    }
}

/// Two drums whose axes stand the weld distance apart, or a hair more,
/// overlapping half their height. Where the wider cap meets the other
/// drum it leaves a crescent that thin, and the rim's arc across it lies
/// inside the cap by the full offset only at its middle angle: asked
/// where it lies inside by less than the boundary's own doubt, the arc
/// still counts as inside, and each drum is split along the two lines the
/// other crosses it on.
#[test]
fn drums_the_weld_distance_apart_in_axis_are_split_apart() {
    let at =
        |x: f64, z: f64| Frame::new(Point::new(x, 0.0, z), Direction::Z, Direction::X, T).unwrap();
    for offset in [1e-5, 1.03e-5, 1.06e-5] {
        let mut model = Model::new();
        let low = ogeom::algo::make_cylinder(&mut model, at(0.0, 0.0), 2.0, 4.0, T)
            .unwrap()
            .shape;
        let high = ogeom::algo::make_cylinder(&mut model, at(offset, 2.0), 2.0, 4.0, T)
            .unwrap()
            .shape;
        let va = volume(&model, &low);
        let common = ogeom::boolean::common(&mut model, &low, &high, T)
            .unwrap_or_else(|e| panic!("{offset}: common: {e}"));
        let cut = ogeom::boolean::cut(&mut model, &low, &high, T)
            .unwrap_or_else(|e| panic!("{offset}: cut: {e}"));
        let fuse = ogeom::boolean::fuse(&mut model, &low, &high, T)
            .unwrap_or_else(|e| panic!("{offset}: fuse: {e}"));
        for (name, made) in [("common", &common), ("cut", &cut), ("fuse", &fuse)] {
            let diagnosis = ogeom::algo::check(&model, &made.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{offset}: {name}: {diagnosis}");
        }
        let (c, k) = (volume(&model, &common.shape), volume(&model, &cut.shape));
        assert!((k + c - va).abs() < 1e-6, "{offset}: {k} + {c}");
        let lens = 8.0 * (offset / 4.0).acos() - offset / 2.0 * (16.0 - offset * offset).sqrt();
        assert!((c - 2.0 * lens).abs() < 1e-6, "{offset}: {c}");
        // The fuse keeps both crescents as faces of their own.
        let faces = explore_unique(&model, &fuse.shape, ShapeType::Face)
            .unwrap()
            .len();
        assert_eq!(faces, 6, "{offset}");
    }
}

/// Faces that meet the other solid's within tolerance without being its
/// coincident partners: a drill whose wall pokes a tenth of a micron out of
/// a box side, a ball resting a tenth of a micron into a box top, and a box
/// stacked on another and tilted a microradian. Each piece read on the
/// other's boundary is settled by asking just off it on both sides, and
/// every operation keeps its volumes.
#[test]
fn faces_within_tolerance_of_the_other_solid_are_settled_from_both_sides() {
    let at = |x: f64, y: f64, z: f64| {
        Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
    };
    let pi = core::f64::consts::PI;
    let cases: Vec<(&str, f64, f64, f64)> = vec![
        // (case, fuse, common, cut)
        (
            "drill",
            64.0 + 6.0 * pi - 4.0 * pi,
            4.0 * pi,
            64.0 - 4.0 * pi,
        ),
        ("ball", 64.0 + 4.0 / 3.0 * pi, 0.0, 64.0),
        ("tilt", 16.0 - 1e-6, 0.0, 8.0 - 1e-6),
    ];
    for (case, fuse_want, common_want, cut_want) in cases {
        let mut model = Model::new();
        let (a, b) = match case {
            "drill" => (
                ogeom::algo::make_box(&mut model, Frame::WORLD, (4.0, 4.0, 4.0), T)
                    .unwrap()
                    .shape,
                ogeom::algo::make_cylinder(&mut model, at(1.0 - 1e-7, 2.0, -1.0), 1.0, 6.0, T)
                    .unwrap()
                    .shape,
            ),
            "ball" => (
                ogeom::algo::make_box(&mut model, Frame::WORLD, (4.0, 4.0, 4.0), T)
                    .unwrap()
                    .shape,
                ogeom::algo::make_sphere(&mut model, at(2.0, 2.0, 5.0 - 1e-7), 1.0, T)
                    .unwrap()
                    .shape,
            ),
            _ => {
                let low = ogeom::algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T)
                    .unwrap()
                    .shape;
                let high = ogeom::algo::make_box(&mut model, at(0.0, 0.0, 2.0), (2.0, 2.0, 2.0), T)
                    .unwrap()
                    .shape;
                let axis = ogeom::math::Axis {
                    location: Point::new(1.0, 1.0, 2.0),
                    direction: Direction::X,
                };
                let tilted = ogeom::algo::transformed(
                    &mut model,
                    &high,
                    ogeom::math::Transform::rotation(axis, 1e-6),
                )
                .unwrap()
                .shape;
                (low, tilted)
            }
        };
        for (name, built, want) in [
            (
                "fuse",
                ogeom::boolean::fuse(&mut model, &a, &b, T),
                fuse_want,
            ),
            (
                "common",
                ogeom::boolean::common(&mut model, &a, &b, T),
                common_want,
            ),
            ("cut", ogeom::boolean::cut(&mut model, &a, &b, T), cut_want),
        ] {
            let shape = built.unwrap_or_else(|e| panic!("{case} {name}: {e}")).shape;
            let v = if explore_unique(&model, &shape, ShapeType::Face)
                .unwrap()
                .is_empty()
            {
                0.0
            } else {
                assert!(
                    ogeom::algo::check(&model, &shape, T).unwrap().is_valid(),
                    "{case} {name}"
                );
                volume(&model, &shape)
            };
            assert!((v - want).abs() < 1e-5, "{case} {name}: {v} against {want}");
        }
    }
}

/// A ball cut from a box it touches from inside, at the top wall or a tenth
/// of a micron through it: one solid with the ball as its void, not the
/// box and an inside-out ball beside it.
#[test]
fn a_cavity_touching_its_wall_is_the_solid_s_void() {
    for poke in [0.0, 1e-7] {
        let mut model = Model::new();
        let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (4.0, 4.0, 4.0), T)
            .unwrap()
            .shape;
        let centre = Frame::new(
            Point::new(2.0, 2.0, 3.0 + poke),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        let ball = ogeom::algo::make_sphere(&mut model, centre, 1.0, T)
            .unwrap()
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &block, &ball, T)
            .unwrap()
            .shape;
        let solids = explore_unique(&model, &cut, ShapeType::Solid).unwrap();
        assert_eq!(solids.len(), 1, "poke {poke}");
        assert_eq!(
            explore_unique(&model, &solids[0], ShapeType::Shell)
                .unwrap()
                .len(),
            2
        );
        assert!(
            ogeom::algo::check(&model, &cut, T).unwrap().is_valid(),
            "poke {poke}"
        );
        let want = 64.0 - 4.0 / 3.0 * core::f64::consts::PI;
        assert!((volume(&model, &cut) - want).abs() < 1e-6, "poke {poke}");
        let inside =
            ogeom::algo::classify_in_solid_exact(&model, &solids[0], Point::new(2.0, 2.0, 3.0), T)
                .unwrap();
        assert_eq!(inside, ogeom::algo::Containment::Out, "poke {poke}");
    }
}

/// Slivers thinner than the weld distance (a hundred confusion distances)
/// are welded: drums a few microns apart in radius, axis or end, and a box
/// tilted five microradians on another. Every operation answers with a
/// valid solid whose volume is the exact one to the sliver it welded, and
/// nothing in it states a tolerance past the weld distance.
#[test]
fn slivers_under_the_weld_distance_are_welded() {
    use ogeom::topo::Shape;
    fn at(x: f64, z: f64) -> Frame {
        Frame::new(Point::new(x, 0.0, z), Direction::Z, Direction::X, T).unwrap()
    }
    let pi = core::f64::consts::PI;
    let drum = |model: &mut Model, frame: Frame, r: f64, h: f64| {
        ogeom::algo::make_cylinder(model, frame, r, h, T)
            .unwrap()
            .shape
    };
    type Build = Box<dyn Fn(&mut Model) -> (Shape, Shape)>;
    let cases: Vec<(&str, Build, [f64; 3])> = vec![
        (
            "radius 5e-6",
            Box::new(move |m: &mut Model| {
                (
                    drum(m, at(0.0, 0.0), 2.0, 4.0),
                    drum(m, at(0.0, 2.0), 2.000005, 4.0),
                )
            }),
            [
                16.0 * pi + 4.0 * pi * 2.000005f64.powi(2) - 8.0 * pi,
                8.0 * pi,
                8.0 * pi,
            ],
        ),
        (
            "axis 1e-6",
            Box::new(move |m: &mut Model| {
                (
                    drum(m, at(0.0, 0.0), 2.0, 4.0),
                    drum(m, at(1e-6, 2.0), 2.0, 4.0),
                )
            }),
            [24.0 * pi, 8.0 * pi, 8.0 * pi],
        ),
        (
            "end 5e-6",
            Box::new(move |m: &mut Model| {
                (
                    drum(m, at(0.0, 0.0), 2.0, 4.0),
                    drum(m, at(0.0, 4.0 - 5e-6), 1.0, 4.0),
                )
            }),
            [20.0 * pi, 0.0, 16.0 * pi],
        ),
        (
            "tilt 5e-6",
            Box::new(|m: &mut Model| {
                let low = ogeom::algo::make_box(m, Frame::WORLD, (2.0, 2.0, 2.0), T)
                    .unwrap()
                    .shape;
                let high = ogeom::algo::make_box(m, at(0.0, 2.0), (2.0, 2.0, 2.0), T)
                    .unwrap()
                    .shape;
                let axis = ogeom::math::Axis {
                    location: Point::new(1.0, 1.0, 2.0),
                    direction: Direction::X,
                };
                let tilted = ogeom::algo::transformed(
                    m,
                    &high,
                    ogeom::math::Transform::rotation(axis, 5e-6),
                )
                .unwrap()
                .shape;
                (low, tilted)
            }),
            [16.0, 0.0, 8.0],
        ),
    ];
    for (case, build, wants) in cases {
        let mut model = Model::new();
        let (a, b) = build(&mut model);
        for (name, built, want) in [
            (
                "fuse",
                ogeom::boolean::fuse(&mut model, &a, &b, T),
                wants[0],
            ),
            (
                "common",
                ogeom::boolean::common(&mut model, &a, &b, T),
                wants[1],
            ),
            ("cut", ogeom::boolean::cut(&mut model, &a, &b, T), wants[2]),
        ] {
            let shape = built.unwrap_or_else(|e| panic!("{case} {name}: {e}")).shape;
            if explore_unique(&model, &shape, ShapeType::Face)
                .unwrap()
                .is_empty()
            {
                assert!(want.abs() < 1e-3, "{case} {name}: nothing, against {want}");
                continue;
            }
            assert!(
                ogeom::algo::check(&model, &shape, T).unwrap().is_valid(),
                "{case} {name}"
            );
            let v = volume(&model, &shape);
            assert!((v - want).abs() < 1e-3, "{case} {name}: {v} against {want}");
            for kind in [ShapeType::Edge, ShapeType::Vertex] {
                for sub in explore_unique(&model, &shape, kind).unwrap() {
                    let stated = model.tolerance_of(&sub).unwrap().map_or(0.0, |t| t.get());
                    assert!(stated <= 1e-5, "{case} {name}: a tolerance of {stated}");
                }
            }
        }
    }
}

/// Blocks tilted by a few hundredths of a milliradian, a sliver of real
/// geometry: stacked and hinged about a horizontal axis through the middle
/// of the joint (the wedge each side is the tilt's own volume, tan t), and
/// standing in one another turned about the vertical (four corner wedges,
/// 4 tan t between them). The two faces near-parallel along the hinge cross on
/// their true line, and every operation keeps its volume.
#[test]
fn blocks_tilted_by_a_hair_keep_their_volumes() {
    for (turn, about_x) in [(1e-4_f64, true), (3e-5, true), (3e-5, false)] {
        let mut model = Model::new();
        let low = ogeom::algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T)
            .unwrap()
            .shape;
        let (base, axis) = if about_x {
            (
                Frame::new(Point::new(0.0, 0.0, 2.0), Direction::Z, Direction::X, T).unwrap(),
                ogeom::math::Axis {
                    location: Point::new(1.0, 1.0, 2.0),
                    direction: Direction::X,
                },
            )
        } else {
            (
                Frame::WORLD,
                ogeom::math::Axis {
                    location: Point::new(1.0, 1.0, 0.0),
                    direction: Direction::Z,
                },
            )
        };
        let other = ogeom::algo::make_box(&mut model, base, (2.0, 2.0, 2.0), T)
            .unwrap()
            .shape;
        let turned = ogeom::algo::transformed(
            &mut model,
            &other,
            ogeom::math::Transform::rotation(axis, turn),
        )
        .unwrap()
        .shape;
        let (fuse_want, common_want, cut_want) = if about_x {
            (16.0 - turn.tan(), turn.tan(), 8.0 - turn.tan())
        } else {
            (
                8.0 + 4.0 * turn.tan(),
                8.0 - 4.0 * turn.tan(),
                4.0 * turn.tan(),
            )
        };
        for (name, built, want) in [
            (
                "fuse",
                ogeom::boolean::fuse(&mut model, &low, &turned, T),
                fuse_want,
            ),
            (
                "common",
                ogeom::boolean::common(&mut model, &low, &turned, T),
                common_want,
            ),
            (
                "cut",
                ogeom::boolean::cut(&mut model, &low, &turned, T),
                cut_want,
            ),
        ] {
            let shape = built
                .unwrap_or_else(|e| panic!("{turn} {about_x} {name}: {e}"))
                .shape;
            assert!(
                ogeom::algo::check(&model, &shape, T).unwrap().is_valid(),
                "{turn} {about_x} {name}"
            );
            let v = volume(&model, &shape);
            assert!(
                (v - want).abs() < 2e-6,
                "{turn} {about_x} {name}: {v} against {want}"
            );
        }
    }
}

/// A compound of solids that do not meet is a tool like any solid: a slab
/// cut by two cubes standing apart through it loses both overlaps, and a
/// plate cut by a grid of drums at once is the plate cut by each in turn.
#[test]
fn a_compound_of_disjoint_solids_cuts_as_one_tool() {
    let at = |x: f64, y: f64, z: f64| {
        Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
    };
    let mut model = Model::new();
    let a = ogeom::algo::make_box(&mut model, at(0.0, 0.0, 0.0), (1.0, 1.0, 1.0), T)
        .unwrap()
        .shape;
    let b = ogeom::algo::make_box(&mut model, at(3.0, 0.0, 0.0), (1.0, 1.0, 1.0), T)
        .unwrap()
        .shape;
    let pair = model.add_compound(&[a, b]).unwrap();
    let slab = ogeom::algo::make_box(&mut model, at(-1.0, -1.0, 0.5), (6.0, 3.0, 1.0), T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(&mut model, &slab, &pair, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &cut, T).unwrap().is_valid());
    assert!((volume(&model, &cut) - 17.0).abs() < 1e-9);

    let plate = ogeom::algo::make_box(&mut model, at(0.0, 0.0, 0.0), (15.0, 15.0, 4.0), T)
        .unwrap()
        .shape;
    let mut drums = Vec::new();
    for i in 0..3 {
        for j in 0..3 {
            let (x, y) = (2.5 + 5.0 * f64::from(i), 2.5 + 5.0 * f64::from(j));
            drums.push(
                ogeom::algo::make_cylinder(&mut model, at(x, y, -1.0), 1.0, 6.0, T)
                    .unwrap()
                    .shape,
            );
        }
    }
    let grid = model.add_compound(&drums).unwrap();
    let bored = ogeom::boolean::cut(&mut model, &plate, &grid, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &bored, T).unwrap().is_valid());
    let want = 15.0 * 15.0 * 4.0 - 9.0 * core::f64::consts::PI * 4.0;
    assert!(
        (volume(&model, &bored) - want).abs() < 1e-6,
        "{}",
        volume(&model, &bored)
    );
}

/// A drill lying tangent to the underside of a thin plate, its lowest line
/// running through the corner where a slot's straight wall meets the floor
/// and its rounded end. The drill's circle on that wall touches the wall's
/// corner at its own start, which the section crossing finds a hair short
/// of a full turn: the arc from there up to the plate's top is read from
/// the circle's start, not shrunk to the hair. Cut and common add up to the
/// plate, and the common is the drill's segment below the plate's top along
/// the plate, less what lies in the slot: half the segment across the
/// straight walls and, past the corner, the segment under the rounded end,
/// integrated here by the midpoint rule.
#[test]
fn a_drill_tangent_under_a_plate_through_a_slot_corner_cuts() {
    let mut model = Model::new();
    let m = &mut model;
    let at = |x, y, z| Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap();
    let (floor, top) = (76.2, 78.107_54);
    let (length, slot_x, slot_width, slot_y) = (241.3, -136.042_4, 19.05, 210.972_400_167_151_8);
    let end = slot_width / 2.0;
    let plate = ogeom::algo::make_box(
        m,
        at(-length, 150.0, floor),
        (length, 120.0, top - floor),
        T,
    )
    .unwrap()
    .shape;
    let straight = ogeom::algo::make_box(
        m,
        at(slot_x, 140.0, floor - 1.0),
        (slot_width, slot_y - 140.0, 4.0),
        T,
    )
    .unwrap()
    .shape;
    let round = ogeom::algo::make_cylinder(m, at(slot_x + end, slot_y, floor - 1.0), end, 4.0, T)
        .unwrap()
        .shape;
    let slot = ogeom::boolean::fuse(m, &straight, &round, T).unwrap().shape;
    let part = ogeom::boolean::cut(m, &plate, &slot, T).unwrap().shape;

    let (r, axis_y) = (23.616_669_978_788_245, 210.972_399_999_999_96);
    let axis_z = floor + r;
    let frame = Frame::new(
        Point::new(-473.877_236_225_169_9, axis_y, axis_z),
        Direction::X,
        Direction::Z,
        T,
    )
    .unwrap();
    let drill = ogeom::algo::make_cylinder(m, frame, r, 586.769_672_450_339_8, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(m, &part, &drill, T).unwrap().shape;
    let common = ogeom::boolean::common(m, &part, &drill, T).unwrap().shape;
    for (name, shape) in [("cut", &cut), ("common", &common)] {
        assert!(
            ogeom::algo::check(m, shape, T).unwrap().is_valid(),
            "{name}: check reports the result invalid"
        );
    }

    // The drill's depth below the plate's top at `y`, and its segment.
    let depth = |y: f64| (top - axis_z + (r * r - (y - axis_y).powi(2)).max(0.0).sqrt()).max(0.0);
    let h = top - floor;
    let segment = r * r * ((r - h) / r).acos() - (r - h) * (2.0 * r * h - h * h).sqrt();
    let steps = 200_000;
    let dy = end / f64::from(steps);
    let under_round: f64 = (0..steps)
        .map(|i| {
            let s = (f64::from(i) + 0.5) * dy;
            2.0 * (end * end - s * s).sqrt() * depth(slot_y + s) * dy
        })
        .sum();
    let want = segment * length - segment / 2.0 * slot_width - under_round;
    let (v_part, v_cut, v_common) = (volume(m, &part), volume(m, &cut), volume(m, &common));
    assert!(
        (v_cut + v_common - v_part).abs() <= 1e-6 * v_part,
        "cut {v_cut} + common {v_common} against the part {v_part}"
    );
    assert!(
        (v_common - want).abs() <= 1e-6 * want,
        "common {v_common} against {want}"
    );
}
