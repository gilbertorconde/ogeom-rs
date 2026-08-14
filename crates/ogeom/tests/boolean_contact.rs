//! The contact configurations an earlier plan still called refusals,
//! pinned as the working behaviour they have become: curved same-domain
//! pairs unify, and contact confined to an edge or a vertex passes through
//! the boolean without harm.
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
/// and the cut — the curved same-domain family, held to closed forms.
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

/// Contact confined to one edge or one vertex: the fuse of edge-touching
/// boxes is one valid solid of both volumes — the shared edge is the
/// non-manifold seam the model permits — and cutting a corner-touching
/// tool removes nothing.
#[test]
fn edge_and_vertex_contact_pass_through_the_boolean() {
    {
        let mut model = Model::new();
        let a = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let f = Frame::new(Point::new(10.0, 0.0, 10.0), Direction::Z, Direction::X, T).unwrap();
        let b = ogeom::algo::make_box(&mut model, f, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let fused = ogeom::boolean::fuse(&mut model, &a, &b, T).unwrap().shape;
        let v = volume(&model, &fused);
        assert!((v - 2000.0).abs() < 1e-6, "both boxes survive: {v}");
        assert!(ogeom::algo::check(&model, &fused, T).unwrap().is_valid());
    }
    {
        let mut model = Model::new();
        let a = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let f = Frame::new(Point::new(10.0, 10.0, 10.0), Direction::Z, Direction::X, T).unwrap();
        let b = ogeom::algo::make_box(&mut model, f, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &a, &b, T).unwrap().shape;
        let v = volume(&model, &cut);
        assert!(
            (v - 1000.0).abs() < 1e-6,
            "a corner touch removes nothing: {v}"
        );
    }
}

/// A ball seated in a bore touches it along the equator and crosses it
/// nowhere. The boolean keeps that curve out of its arithmetic — a contact
/// carries no parity, so nothing is inside on one side of it — and the
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

/// A sphere's chart is bounded above and below by *poles* — edges that are
/// points in space — and by a seam it meets twice. Leave the poles out of
/// the arrangement and the chart has no top or bottom, so nothing can be
/// arranged inside it and every boolean over a ball fails, whether or not
/// the ball is anywhere near the other solid. They are in it now.
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
/// no material — there is no side of it that is inside on one hand and
/// outside on the other — so the boolean carries the touch instead of
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
    // carry into the other wall's chart. The same-domain melt used to refuse
    // exactly here; now the stored pcurve — which on the identical chart
    // already is the projection — stands in, and the fuse of a contained
    // solid comes out as the container.
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
    // The budget is two separately meshed fitted trims at this chord, not
    // the melt: the fuse either resolves the contact or refuses by name.
    let measured = volume(&model, &fused.shape);
    assert!(
        (measured - before).abs() < 2e-2,
        "fusing a contained drum should give the container: {measured} vs {before}"
    );
    // Every face of the contained drum melted rather than surviving as a
    // skin. The wall stays split where the contained drum's rim touched it —
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
/// so it has to reach the melt — which it only does while the scale is
/// carried as the placement it is, since a restated plane no longer says it
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
/// exactly one of the two descriptions may survive — the question `cut` never
/// has to ask, which is why it closed while `common` did not.
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
    // wall between them. Three faces — a fourth would be the shared disk
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

    // The neighbour that always worked, kept working.
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
/// restated as patches — correctly; there is no other way to carry a box's
/// planes under one. What the patches lose is the *word* plane, and
/// coincidence is decided on what the geometry says: nothing answers `Same`
/// for two patches, so the pair went to the marcher, which documents that it
/// is not a coincidence detector and has no crossing to trace here. It came
/// back with a section made of noise, and the boolean refused some stages
/// later for edge or vertex contact, which this is not.
#[test]
fn a_sheared_copy_sharing_a_plane_is_refused_as_the_same_domain_contact_it_is() {
    let mut model = Model::with_tolerances(T);
    let a = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    // Unit determinant, so the copy keeps its volume, and slid clear along
    // the z = 0 plane the two go on sharing.
    let shear = ogeom::math::GeneralTransform {
        linear: ogeom::math::Matrix3 {
            rows: [[1.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        },
        translation: ogeom::math::Vector::new(5.0, 0.0, 0.0),
    };
    let b = ogeom::algo::general_transformed_shape(&mut model, &a, &shear, T)
        .unwrap()
        .shape;

    let refused = ogeom::boolean::fuse(&mut model, &a, &b, T)
        .expect_err("a same-domain pair with no closed-form pcurve is not resolved yet");

    let said = refused.to_string();
    assert!(
        said.contains("same-domain contact"),
        "refused as something else: {said}"
    );
    assert!(
        !said.contains("edge or vertex contact"),
        "still refusing by the name of a different configuration: {said}"
    );
}
