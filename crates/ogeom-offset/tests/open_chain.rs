//! An open chain thickened to one side or both, and a rib grown from one,
//! measured against closed-form areas and volumes and the band's own
//! definition: every point of a side stands at its offset from the chain.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Plane, Point};
use ogeom_mesh::Deflection;
use ogeom_offset::{Cap, Join, make_rib, thicken_open_wire};
use ogeom_topo::{Filter, Model, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

fn p(x: f64, y: f64) -> Point {
    Point::new(x, y, 0.0)
}

fn polyline(model: &mut Model, corners: &[Point]) -> Shape {
    ogeom_algo::make_polygon(model, corners, false, T)
        .unwrap()
        .shape
}

/// An arc of `radius` about `centre` in the xy plane, counter-clockwise
/// from angle `a0` to `a1`, between the given vertices.
fn arc(
    model: &mut Model,
    centre: Point,
    radius: f64,
    (a0, a1): (f64, f64),
    ends: (&Shape, &Shape),
) -> Shape {
    let frame = Frame::new(centre, Direction::Z, Direction::X, T).unwrap();
    let circle = ogeom_math::Circle::new(frame, radius, T).unwrap();
    let curve = ogeom_geom::Curve::Circle(ogeom_geom::CircleCurve::new(circle));
    ogeom_algo::make_edge_between(model, curve, (a0, a1), ends.0, ends.1, T)
        .unwrap()
        .shape
}

fn segment(model: &mut Model, from: &Shape, to: &Shape, a: Point, b: Point) -> Shape {
    let curve = ogeom_geom::Curve::Line(ogeom_geom::LineCurve::segment(a, b, T).unwrap());
    let domain = ogeom_geom::Curve3d::domain(&curve);
    ogeom_algo::make_edge_between(model, curve, domain, from, to, T)
        .unwrap()
        .shape
}

fn vertex(model: &mut Model, at: Point) -> Shape {
    ogeom_algo::make_vertex(model, at).shape
}

/// A quarter arc of radius 5 about the origin, counter-clockwise from
/// (5, 0) to (0, 5).
fn quarter_arc(model: &mut Model) -> Shape {
    let (a, b) = (vertex(model, p(5.0, 0.0)), vertex(model, p(0.0, 5.0)));
    let edge = arc(model, Point::ORIGIN, 5.0, (0.0, PI / 2.0), (&a, &b));
    ogeom_algo::make_wire(model, &[edge], T).unwrap().shape
}

/// A line along x from the origin to (10, 0), then an arc of `radius`
/// about (10 + radius, 0) turning down from (10, 0) through a quarter.
fn line_then_arc(model: &mut Model, radius: f64) -> Shape {
    let v0 = vertex(model, p(0.0, 0.0));
    let v1 = vertex(model, p(10.0, 0.0));
    let v2 = vertex(model, p(10.0 + radius, -radius));
    let line = segment(model, &v0, &v1, p(0.0, 0.0), p(10.0, 0.0));
    let turn = arc(
        model,
        p(10.0 + radius, 0.0),
        radius,
        (PI, 1.5 * PI),
        (&v1, &v2),
    );
    ogeom_algo::make_wire(model, &[line, turn], T)
        .unwrap()
        .shape
}

/// The area the outline bounds, from the planar face it bounds.
fn area(model: &mut Model, outline: &Shape) -> f64 {
    let faces =
        ogeom_algo::faces_from_planar_wires(model, Plane::XY, std::slice::from_ref(outline), T)
            .unwrap()
            .faces;
    assert_eq!(faces.len(), 1);
    ogeom_algo::surface_properties(model, &faces[0], Deflection::default(), T)
        .unwrap()
        .mass
}

/// Points along every edge of `shape`, 65 per edge, ends included.
fn samples(model: &Model, shape: &Shape) -> Vec<Point> {
    let mut out = Vec::new();
    for edge in explore(model, shape, Filter::OfType(ShapeType::Edge)).unwrap() {
        let (curve, range) = ogeom_algo::placed_edge_curve(model, &edge, T)
            .unwrap()
            .unwrap();
        for k in 0..=64 {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 64.0;
            out.push(curve.point_at(t, T).unwrap());
        }
    }
    out
}

fn vertices(model: &Model, shape: &Shape) -> Vec<Point> {
    samples(model, shape)
        .into_iter()
        .enumerate()
        .filter(|(i, _)| i % 65 == 0 || i % 65 == 64)
        .map(|(_, q)| q)
        .collect()
}

fn has_vertex(model: &Model, shape: &Shape, at: Point) -> bool {
    vertices(model, shape).iter().any(|v| v.distance(at) < 1e-9)
}

/// Distance from `q` to the segment `a`..`b`.
fn to_segment(q: Point, a: Point, b: Point) -> f64 {
    let d = b - a;
    let t = ((q - a).dot(d) / d.dot(d)).clamp(0.0, 1.0);
    q.distance(a + d * t)
}

fn to_polyline(q: Point, corners: &[Point]) -> f64 {
    corners
        .windows(2)
        .map(|w| to_segment(q, w[0], w[1]))
        .fold(f64::INFINITY, f64::min)
}

#[test]
fn a_straight_chain_thickened_one_side_is_the_rectangle() {
    let mut model = Model::new();
    let chain = polyline(&mut model, &[p(0.0, 0.0), p(10.0, 0.0)]);
    let band = thicken_open_wire(
        &mut model,
        &chain,
        Plane::XY,
        2.0,
        0.0,
        Cap::Flat,
        Join::Intersection,
        T,
    )
    .unwrap()
    .shape;
    for corner in [p(0.0, 0.0), p(10.0, 0.0), p(10.0, 2.0), p(0.0, 2.0)] {
        assert!(has_vertex(&model, &band, corner), "no corner at {corner:?}");
    }
    let measured = area(&mut model, &band);
    assert!((measured - 20.0).abs() < 1e-9, "area {measured}, want 20");
    assert!(ogeom_algo::check(&model, &band, T).unwrap().is_valid());
}

#[test]
fn the_caps_close_a_straddling_band_flat_round_and_extended() {
    // Width 2 about the chain; a round cap adds a half disc of radius 1 at
    // each end, an extended one a 1 by 2 rectangle.
    for (cap, want) in [
        (Cap::Flat, 20.0),
        (Cap::Round, 20.0 + PI),
        (Cap::Extended, 24.0),
    ] {
        let mut model = Model::new();
        let chain = polyline(&mut model, &[p(0.0, 0.0), p(10.0, 0.0)]);
        let band = thicken_open_wire(
            &mut model,
            &chain,
            Plane::XY,
            -1.0,
            1.0,
            cap,
            Join::Intersection,
            T,
        )
        .unwrap()
        .shape;
        let measured = area(&mut model, &band);
        assert!(
            (measured - want).abs() < 1e-6,
            "{cap:?}: area {measured}, want {want}"
        );
    }
}

#[test]
fn an_l_chain_thickened_inside_its_turn_is_trimmed_and_flat_ended() {
    // The chain turns left at (10, 0); its left side 2 away is on the
    // inside of the turn: two 2-wide strips sharing their 2 by 2 corner.
    let mut model = Model::new();
    let chain = polyline(&mut model, &[p(0.0, 0.0), p(10.0, 0.0), p(10.0, 8.0)]);
    let band = thicken_open_wire(
        &mut model,
        &chain,
        Plane::XY,
        2.0,
        0.0,
        Cap::Flat,
        Join::Intersection,
        T,
    )
    .unwrap()
    .shape;
    let measured = area(&mut model, &band);
    assert!((measured - 32.0).abs() < 1e-9, "area {measured}, want 32");
    // The end edges are square to the chain through its end points.
    for corner in [
        p(0.0, 0.0),
        p(0.0, 2.0),
        p(10.0, 8.0),
        p(8.0, 8.0),
        p(8.0, 2.0),
    ] {
        assert!(has_vertex(&model, &band, corner), "no corner at {corner:?}");
    }
    assert!(ogeom_algo::check(&model, &band, T).unwrap().is_valid());
}

#[test]
fn an_l_chain_thickened_outside_its_turn_gains_the_mitre() {
    // Right side 2 away, outside the turn: 2 x the chain's length 18, plus
    // the corner's 2 by 2 mitre, or its quarter disc with the arc join.
    for (join, want) in [(Join::Intersection, 40.0), (Join::Arc, 36.0 + PI)] {
        let mut model = Model::new();
        let corners = [p(0.0, 0.0), p(10.0, 0.0), p(10.0, 8.0)];
        let chain = polyline(&mut model, &corners);
        let band = thicken_open_wire(&mut model, &chain, Plane::XY, 0.0, -2.0, Cap::Flat, join, T)
            .unwrap()
            .shape;
        let measured = area(&mut model, &band);
        assert!(
            (measured - want).abs() < 1e-6,
            "{join:?}: area {measured}, want {want}"
        );
        if join == Join::Intersection {
            assert!(has_vertex(&model, &band, p(12.0, -2.0)), "no mitre corner");
        } else {
            // Every point of the outline is on the chain, 2 from it, or on
            // a flat cap (x = 0 at the start, y = 8 at the end).
            let mut worst: f64 = 0.0;
            for q in samples(&model, &band) {
                let d = to_polyline(q, &corners);
                let on_cap = q.x.abs() < 1e-12 || (q.y - 8.0).abs() < 1e-12;
                if !on_cap {
                    worst = worst.max(d.min((d - 2.0).abs()));
                }
            }
            assert!(worst < 1e-9, "a side point stands {worst} off its offset");
        }
    }
}

#[test]
fn an_arc_chain_thickens_into_an_annular_sector() {
    // Counter-clockwise, the chain's left is toward the centre.
    for (left, right, inner, outer) in [(1.0, 0.0, 4.0, 5.0), (0.0, -1.0, 5.0, 6.0)] {
        let mut model = Model::new();
        let chain = quarter_arc(&mut model);
        let band = thicken_open_wire(
            &mut model,
            &chain,
            Plane::XY,
            left,
            right,
            Cap::Flat,
            Join::Intersection,
            T,
        )
        .unwrap()
        .shape;
        let want = PI / 4.0 * (outer * outer - inner * inner);
        let measured = area(&mut model, &band);
        assert!(
            (measured - want).abs() < 1e-6,
            "area {measured}, want {want}"
        );
        // Both sides concentric: every point is at the inner or outer
        // radius, or on a cap along an axis.
        let mut worst: f64 = 0.0;
        for q in samples(&model, &band) {
            let r = q.distance(Point::ORIGIN);
            let on_cap = q.x.abs() < 1e-12 || q.y.abs() < 1e-12;
            if !on_cap {
                worst = worst.max((r - inner).abs().min((r - outer).abs()));
            }
        }
        assert!(worst < 1e-9, "a side point stands {worst} off its radius");
        assert!(ogeom_algo::check(&model, &band, T).unwrap().is_valid());
    }
}

#[test]
fn a_thickness_past_the_arcs_radius_is_refused_by_name() {
    for left in [5.0, 6.0] {
        let mut model = Model::new();
        let chain = quarter_arc(&mut model);
        let err = thicken_open_wire(
            &mut model,
            &chain,
            Plane::XY,
            left,
            0.0,
            Cap::Flat,
            Join::Intersection,
            T,
        )
        .unwrap_err()
        .to_string();
        assert!(err.contains("past the centre"), "{err}");
    }
}

#[test]
fn an_intersection_join_extends_an_arc_to_its_meeting_or_refuses() {
    // A line along x, then an arc about (15, 0) turning down from (10, 0):
    // a right turn, so the left side gaps. The line's offset y = 2 meets
    // the arc's inner offset (radius 3) at x = 15 - sqrt(5).
    let mut model = Model::new();
    let chain = line_then_arc(&mut model, 5.0);
    let band = thicken_open_wire(
        &mut model,
        &chain,
        Plane::XY,
        2.0,
        0.0,
        Cap::Flat,
        Join::Intersection,
        T,
    )
    .unwrap()
    .shape;
    let met = p(15.0 - 5.0_f64.sqrt(), 2.0);
    assert!(has_vertex(&model, &band, met), "no meeting at {met:?}");
    assert!(ogeom_algo::check(&model, &band, T).unwrap().is_valid());

    // A tighter arc: radius 1, its inner side 0.6 in has radius 0.4 and
    // the line's offset y = 0.6 passes clear of it.
    let mut model = Model::new();
    let chain = line_then_arc(&mut model, 1.0);
    let err = thicken_open_wire(
        &mut model,
        &chain,
        Plane::XY,
        0.6,
        0.0,
        Cap::Flat,
        Join::Intersection,
        T,
    )
    .unwrap_err()
    .to_string();
    assert!(err.contains("no meeting"), "{err}");
}

#[test]
fn a_closed_wire_or_a_band_of_no_width_is_refused() {
    let mut model = Model::new();
    let closed = ogeom_algo::make_polygon(
        &mut model,
        &[p(0.0, 0.0), p(4.0, 0.0), p(4.0, 3.0)],
        true,
        T,
    )
    .unwrap()
    .shape;
    assert!(
        thicken_open_wire(
            &mut model,
            &closed,
            Plane::XY,
            1.0,
            0.0,
            Cap::Flat,
            Join::Arc,
            T
        )
        .is_err()
    );
    let chain = polyline(&mut model, &[p(0.0, 0.0), p(10.0, 0.0)]);
    assert!(
        thicken_open_wire(
            &mut model,
            &chain,
            Plane::XY,
            1.0,
            1.0,
            Cap::Flat,
            Join::Arc,
            T
        )
        .is_err()
    );
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom_algo::volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

fn cuboid(model: &mut Model, at: Point, size: (f64, f64, f64)) -> Shape {
    ogeom_algo::make_box(model, Frame::WORLD.with_origin(at), size, T)
        .unwrap()
        .shape
}

/// A 30 x 10 x 20 block with a slot 10 wide (x from 10 to 20) cut down
/// from the top to a floor at z = 5, through in y.
fn channel(model: &mut Model) -> Shape {
    let block = cuboid(model, Point::ORIGIN, (30.0, 10.0, 20.0));
    let slot = cuboid(model, Point::new(10.0, -1.0, 5.0), (10.0, 12.0, 16.0));
    ogeom_bool::cut(model, &block, &slot, T).unwrap().shape
}

const CHANNEL: f64 = 30.0 * 10.0 * 20.0 - 10.0 * 10.0 * 15.0;

#[test]
fn a_rib_across_a_channel_fills_the_slot_below_its_chain() {
    // A chain across the slot at height z, grown down to the floor at
    // z = 5: the web is 2 thick, 10 wide (the slot) and z - 5 deep; upward
    // it runs out of the slot and meets nothing. At z = 15 the chain stops
    // short of both walls and is carried on into them; at z = 20 it runs
    // from wall to wall along the slot's top edges.
    for (z, x0, x1) in [(15.0, 12.0, 18.0), (20.0, 10.0, 20.0)] {
        let web = 2.0 * 10.0 * (z - 5.0);
        for along_normal in [false, true] {
            let mut model = Model::new();
            let solid = channel(&mut model);
            let plane = if along_normal {
                // In the plane of the chain's height, the band 1 each side
                // of y = 5.
                Plane::through(Point::new(0.0, 0.0, z), Direction::Z)
            } else {
                // In the plane y = 5, thickened 1 each way across it.
                Plane::through(Point::new(0.0, 5.0, 0.0), Direction::Y)
            };
            let chain = ogeom_algo::make_polygon(
                &mut model,
                &[Point::new(x0, 5.0, z), Point::new(x1, 5.0, z)],
                false,
                T,
            )
            .unwrap()
            .shape;
            let rib = make_rib(
                &mut model,
                &solid,
                &chain,
                plane,
                (1.0, 1.0),
                along_normal,
                T,
            )
            .unwrap()
            .shape;
            let measured = volume(&model, &rib);
            let want = CHANNEL + web;
            assert!(
                (measured - want).abs() < want * 1e-9,
                "z {z}, along_normal {along_normal}: volume {measured}, want {want}"
            );
            let diagnosis = ogeom_algo::check(&model, &rib, T).unwrap();
            assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        }
    }
}

#[test]
fn a_gusset_past_a_short_wall_stays_within_the_chains_end() {
    // A plate (z 0 to 2) with a wall at its left end only 10 high (x 0 to
    // 2). The chain runs from (2, 14), above the wall, down to the plate at
    // (28, 2). Its start points at nothing, so the web ends square to the
    // growth through it: the line from (2, 14) along the growth meets the
    // wall's top, and the web is the polygon chain, plate, wall face, wall
    // top, that line. It does not run round the wall's top.
    let mut model = Model::new();
    let plate = cuboid(&mut model, Point::ORIGIN, (30.0, 10.0, 2.0));
    let wall = cuboid(&mut model, Point::ORIGIN, (2.0, 10.0, 10.0));
    let solid = ogeom_bool::fuse(&mut model, &plate, &wall, T)
        .unwrap()
        .shape;
    let solids = explore(&model, &solid, Filter::OfType(ShapeType::Solid)).unwrap();
    let solid = solids[0].clone();
    let base = volume(&model, &solid);
    assert!((base - (600.0 + 160.0)).abs() < 1e-9, "base volume {base}");

    let (a, e) = ((2.0, 14.0), (28.0, 2.0));
    let chain = ogeom_algo::make_polygon(
        &mut model,
        &[Point::new(a.0, 5.0, a.1), Point::new(e.0, 5.0, e.1)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let plane = Plane::through(Point::new(0.0, 5.0, 0.0), Direction::Y);
    let rib = make_rib(&mut model, &solid, &chain, plane, (0.5, 0.5), false, T)
        .unwrap()
        .shape;

    // The growth is square to the chord, toward the corner.
    let (ux, uz) = (e.0 - a.0, e.1 - a.1);
    let len = ux.hypot(uz);
    let (dx, dz) = (uz / len, -ux / len);
    let (dx, dz) = if dx < 0.0 { (dx, dz) } else { (-dx, -dz) };
    let s = (10.0 - a.1) / dz;
    let b = (a.0 + dx * s, 10.0);
    let polygon = [a, e, (2.0, 2.0), (2.0, 10.0), b];
    let mut twice = 0.0;
    for i in 0..polygon.len() {
        let (p0, p1) = (polygon[i], polygon[(i + 1) % polygon.len()]);
        twice += p0.0 * p1.1 - p1.0 * p0.1;
    }
    let web = twice.abs() / 2.0 * 1.0;
    let measured = volume(&model, &rib);
    let want = base + web;
    assert!(
        (measured - want).abs() < want * 1e-9,
        "volume {measured}, want {want} (web {web})"
    );
    let diagnosis = ogeom_algo::check(&model, &rib, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
}

#[test]
fn a_rib_whose_chain_meets_no_wall_or_lies_inside_is_refused_by_name() {
    let mut model = Model::new();
    let solid = cuboid(&mut model, Point::ORIGIN, (10.0, 10.0, 10.0));
    // In the plane y = 50, far beside the block.
    let chain = ogeom_algo::make_polygon(
        &mut model,
        &[Point::new(0.0, 50.0, 5.0), Point::new(10.0, 50.0, 5.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let plane = Plane::through(Point::new(0.0, 50.0, 0.0), Direction::Y);
    let err = make_rib(&mut model, &solid, &chain, plane, (1.0, 1.0), false, T)
        .unwrap_err()
        .to_string();
    assert!(err.contains("meets no wall"), "{err}");

    // Wholly inside the block, in the plane y = 5: nothing to add.
    let inside = ogeom_algo::make_polygon(
        &mut model,
        &[Point::new(2.0, 5.0, 5.0), Point::new(8.0, 5.0, 5.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let plane = Plane::through(Point::new(0.0, 5.0, 0.0), Direction::Y);
    let err = make_rib(&mut model, &solid, &inside, plane, (1.0, 1.0), false, T)
        .unwrap_err()
        .to_string();
    assert!(err.contains("inside the solid"), "{err}");
}
