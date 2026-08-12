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
