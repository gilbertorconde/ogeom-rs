//! Pipe shells that have closed forms: a flat profile edge down a straight
//! leg is a plane, a circle down a line a drum; a ring walked by reversed
//! edges sweeps as the same ring walked forward; a profile square to a
//! curved spine's start sweeps along it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/walks.rs"]
mod walks;

use ogeom::algo::{make_edge, make_face, make_polygon, make_wire, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve, Curve3d as _, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Circle, Direction, Frame, Plane, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn line(model: &mut Model, a: Point, b: Point) -> Shape {
    let curve = LineCurve::segment(a, b, T).unwrap();
    let range = curve.domain();
    make_edge(model, Curve::Line(curve), range, T)
        .unwrap()
        .shape
}

fn xy_face(model: &mut Model, wire: Shape) -> Shape {
    let plane = PlaneSurface::new(Plane::new(Frame::WORLD));
    make_face(model, plane.into(), &[wire], T).unwrap().shape
}

fn square(model: &mut Model, half: f64) -> Shape {
    let pts = [(-half, -half), (half, -half), (half, half), (-half, half)]
        .map(|(x, y)| Point::new(x, y, 0.0));
    let wire = make_polygon(model, &pts, true, T).unwrap().shape;
    xy_face(model, wire)
}

/// Volume, and whether it was integrated exactly.
fn measure(model: &Model, shape: &Shape) -> (f64, bool) {
    let p = volume_properties(model, shape, Deflection::default(), T).unwrap();
    (p.mass, p.deflection == 0.0)
}

fn planes_only(model: &Model, shape: &Shape) -> bool {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .iter()
        .all(|f| {
            let Some(NodeData::Face(d)) = model.node(f).map(|n| n.data()) else {
                return false;
            };
            matches!(
                model.geometry().surface(d.surface),
                Some(SurfaceGeometry::Plane(_))
            )
        })
}

#[test]
fn a_square_down_an_l_is_planes_and_measures_exactly() {
    let mut model = Model::new();
    let profile = square(&mut model, 2.0);
    let spine = make_polygon(
        &mut model,
        &[
            Point::ORIGIN,
            Point::new(0.0, 0.0, 20.0),
            Point::new(20.0, 0.0, 20.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let pipe = ogeom::offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    assert!(planes_only(&model, &pipe));
    let (v, exact) = measure(&model, &pipe);
    assert!(exact, "measured from a mesh");
    assert!((v - 640.0).abs() < 640.0 * 1e-9, "{v}");
}

#[test]
fn a_ring_of_reversed_edges_sweeps_as_the_forward_ring() {
    let mut model = Model::new();
    let pts =
        [(-2.0, -2.0), (2.0, -2.0), (2.0, 2.0), (-2.0, 2.0)].map(|(x, y)| Point::new(x, y, 0.0));
    let forward = make_polygon(&mut model, &pts, true, T).unwrap().shape;
    let edges = model.ordered_children_of(&forward).unwrap();
    let reversed: Vec<Shape> = edges.iter().rev().map(Shape::reversed).collect();
    let wire = make_wire(&mut model, &reversed, T).unwrap().shape;
    let profile = xy_face(&mut model, wire);
    let e = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, 20.0));
    let spine = make_wire(&mut model, &[e], T).unwrap().shape;
    let pipe = ogeom::offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    let (v, _) = measure(&model, &pipe);
    assert!((v - 320.0).abs() < 320.0 * 1e-9, "{v}");
}

#[test]
fn a_circle_down_a_line_is_a_drum() {
    let mut model = Model::new();
    let circle = Circle::new(Frame::WORLD, 2.0, T).unwrap();
    let e = make_edge(
        &mut model,
        CircleCurve::new(circle).into(),
        (0.0, core::f64::consts::TAU),
        T,
    )
    .unwrap()
    .shape;
    let wire = make_wire(&mut model, &[e], T).unwrap().shape;
    let profile = xy_face(&mut model, wire);
    let s = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, 20.0));
    let spine = make_wire(&mut model, &[s], T).unwrap().shape;
    let pipe = ogeom::offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    let (v, exact) = measure(&model, &pipe);
    let want = core::f64::consts::PI * 4.0 * 20.0;
    assert!(exact, "measured from a mesh");
    assert!((v - want).abs() < want * 1e-9, "{v} against {want}");
}

#[test]
fn a_square_square_to_an_arcs_start_sweeps_along_it() {
    let quarter = core::f64::consts::FRAC_PI_2;
    // A quarter circle from the origin to (10, 0, 10), centred on
    // (10, 0, 0) in the XZ plane: its start tangent is +Z. Walked forward,
    // and as the same arc stored the other way round and used reversed.
    for backward in [false, true] {
        let mut model = Model::new();
        let profile = square(&mut model, 2.0);
        let (x, range) = if backward {
            (Direction::Z, (0.0, quarter))
        } else {
            (-Direction::X, (0.0, quarter))
        };
        let frame = Frame::new(
            Point::new(10.0, 0.0, 0.0),
            if backward {
                -Direction::Y
            } else {
                Direction::Y
            },
            x,
            T,
        )
        .unwrap();
        let circle = Circle::new(frame, 10.0, T).unwrap();
        let edge = make_edge(&mut model, CircleCurve::new(circle).into(), range, T)
            .unwrap()
            .shape;
        let edge = if backward { edge.reversed() } else { edge };
        let spine = make_wire(&mut model, &[edge], T).unwrap().shape;
        let pipe = ogeom::offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
            .unwrap_or_else(|e| panic!("backward {backward}: {e}"))
            .shape;
        let (v, exact) = measure(&model, &pipe);
        let want = 16.0 * quarter * 10.0;
        assert!(exact, "measured from a mesh");
        assert!((v - want).abs() < want * 1e-9, "{v} against {want}");
    }
}

/// A square centred on a conical helix and square to it, walked either
/// way: valid, every edge walked once each way, caps included, and the
/// volume its area times the helix's length.
#[test]
fn a_profile_square_to_a_conical_helix_sweeps() {
    for back in [false, true] {
        let mut model = Model::new();
        let frame = Frame::WORLD;
        let taper = 5.0 * 10f64.to_radians().tan();
        let helix = ogeom::geom::HelixCurve::conical(
            frame,
            10.0,
            5.0,
            taper,
            0.0,
            2.0 * core::f64::consts::TAU,
        )
        .unwrap();
        let range = helix.domain();
        let start = helix.point_at(range.0, T).unwrap();
        let tangent = helix.d1_at(range.0, T).unwrap();
        let length = simpson(range.0, range.1, 4000, |t| {
            helix.d1_at(t, T).unwrap().magnitude()
        });
        let edge = make_edge(&mut model, helix.into(), range, T).unwrap().shape;
        let spine = make_wire(&mut model, &[edge], T).unwrap().shape;
        // A small square square to the start tangent, centred on the start.
        let normal = Direction::new(tangent, T).unwrap();
        let plane = Plane::through(start, normal);
        let (u, w) = (plane.frame().x().vector(), plane.frame().y().vector());
        let mut pts: Vec<Point> = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]
            .iter()
            .map(|(a, b)| start + u * *a + w * *b)
            .collect();
        if back {
            pts.reverse();
        }
        let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
        let face = make_face(&mut model, PlaneSurface::new(plane).into(), &[wire], T)
            .unwrap()
            .shape;
        let pipe = ogeom::offset::make_pipe_shell(&mut model, &face, &spine, false, 1e-3, T)
            .unwrap()
            .shape;
        assert_eq!(walks::edges_walked_one_way(&model, &pipe), 0, "back {back}");
        let found = ogeom::algo::check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "back {back}: {found}");
        let fine = Deflection {
            chord: 1e-3,
            ..Deflection::default()
        };
        let v = volume_properties(&model, &pipe, fine, T).unwrap().mass;
        assert!(
            (v - length).abs() < length * 1e-3,
            "back {back}: {v} against {length}"
        );
    }
}

/// Simpson's rule over `[a, b]` in `n` (even) steps.
fn simpson(a: f64, b: f64, n: usize, f: impl Fn(f64) -> f64) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let h = (b - a) / (n as f64);
    let mut sum = f(a) + f(b);
    for i in 1..n {
        #[allow(clippy::cast_precision_loss)]
        let x = a + h * (i as f64);
        sum += f(x) * if i % 2 == 1 { 4.0 } else { 2.0 };
    }
    sum * h / 3.0
}

#[test]
fn a_placed_profile_square_to_an_arc_sweeps_along_it() {
    let mut model = Model::new();
    let profile = square(&mut model, 2.0);
    let profile = ogeom::algo::transformed(&mut model, &profile, ogeom::math::Transform::IDENTITY)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(10.0, 0.0, 0.0), Direction::Y, -Direction::X, T).unwrap();
    let circle = Circle::new(frame, 10.0, T).unwrap();
    let quarter = core::f64::consts::FRAC_PI_2;
    let edge = make_edge(
        &mut model,
        CircleCurve::new(circle).into(),
        (0.0, quarter),
        T,
    )
    .unwrap()
    .shape;
    let spine = make_wire(&mut model, &[edge], T).unwrap().shape;
    let pipe = ogeom::offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    let (v, exact) = measure(&model, &pipe);
    let want = 16.0 * quarter * 10.0;
    assert!(exact);
    assert!((v - want).abs() < want * 1e-9, "{v} against {want}");
}

#[test]
/// A square, holed or not, round a closed ellipse under either frame law:
/// centred on the spine, it sweeps its area times the ellipse's length, the
/// hole a void tunnel whichever way its ring is wound.
fn a_holed_profile_rounds_an_elliptic_ring() {
    let length = {
        let n = 20_000;
        (0..n)
            .map(|i| {
                let t = core::f64::consts::TAU * (f64::from(i) + 0.5) / f64::from(n);
                (20.0 * t.sin()).hypot(15.0 * t.cos())
            })
            .sum::<f64>()
            * core::f64::consts::TAU
            / f64::from(n)
    };
    for frenet in [false, true] {
        for holed in [false, true] {
            let mut model = Model::new();
            // An ellipse 20 by 15 in the XY plane, closed.
            let ellipse = ogeom::math::Ellipse::new(Frame::WORLD, 20.0, 15.0, T).unwrap();
            let e = make_edge(
                &mut model,
                ogeom::geom::EllipseCurve::new(ellipse).into(),
                (0.0, core::f64::consts::TAU),
                T,
            )
            .unwrap()
            .shape;
            let spine = make_wire(&mut model, &[e], T).unwrap().shape;
            // A 4 x 4 square in the XZ plane at (20, 0, 0), square to +Y.
            let pts = [(-2.0, -2.0), (2.0, -2.0), (2.0, 2.0), (-2.0, 2.0)]
                .map(|(a, b)| Point::new(20.0 + a, 0.0, b));
            let outer = make_polygon(&mut model, &pts, true, T).unwrap().shape;
            let frame =
                Frame::new(Point::new(20.0, 0.0, 0.0), Direction::Y, Direction::X, T).unwrap();
            let mut wires = vec![outer];
            if holed {
                let c = Circle::new(frame, 1.0, T).unwrap();
                let he = make_edge(
                    &mut model,
                    CircleCurve::new(c).into(),
                    (0.0, core::f64::consts::TAU),
                    T,
                )
                .unwrap()
                .shape;
                wires.push(make_wire(&mut model, &[he], T).unwrap().shape);
            }
            let face = make_face(
                &mut model,
                PlaneSurface::new(Plane::new(frame)).into(),
                &wires,
                T,
            )
            .unwrap()
            .shape;
            let r = ogeom::offset::make_pipe_shell(&mut model, &face, &spine, frenet, 1e-3, T);
            let pipe = r
                .unwrap_or_else(|e| panic!("frenet {frenet} holed {holed}: {e}"))
                .shape;
            let (v, _) = measure(&model, &pipe);
            let area = if holed {
                16.0 - core::f64::consts::PI
            } else {
                16.0
            };
            let want = area * length;
            assert!(
                (v - want).abs() < want * 2e-3,
                "frenet {frenet} holed {holed}: {v} against {want}"
            );
        }
    }
}

/// A disc piped round an L of straight legs is exact whatever the circle's
/// seam: each leg a drum run on past the corner and trimmed at the mitre,
/// so the tube measures its area times the path's length.
#[test]
fn a_disc_round_an_l_is_exact_wherever_its_seam() {
    use ogeom::math::{Transform, Vector};
    let disc = |model: &mut Model, frame: Frame| {
        let circle = Circle::new(frame, 1.0, T).unwrap();
        let e = make_edge(
            model,
            CircleCurve::new(circle).into(),
            (0.0, core::f64::consts::TAU),
            T,
        )
        .unwrap()
        .shape;
        let w = make_wire(model, &[e], T).unwrap().shape;
        make_face(model, PlaneSurface::new(Plane::new(frame)).into(), &[w], T)
            .unwrap()
            .shape
    };
    for placed in [true, false] {
        let mut model = Model::new();
        let spine = make_polygon(
            &mut model,
            &[
                Point::new(10.0, 10.0, 0.0),
                Point::new(0.0, 10.0, 0.0),
                Point::ORIGIN,
            ],
            false,
            T,
        )
        .unwrap()
        .shape;
        let profile = if placed {
            // The disc in the XZ plane, turned a quarter about Z and moved to
            // the spine's start: its seam on the far side.
            let f = disc(
                &mut model,
                Frame::new(Point::ORIGIN, -Direction::Y, Direction::X, T).unwrap(),
            );
            let turn = Transform::rotation(
                ogeom::math::Axis::new(Point::ORIGIN, Direction::Z),
                -core::f64::consts::FRAC_PI_2,
            );
            let moved = Transform::translation(Vector::new(10.0, 10.0, 0.0)) * turn;
            ogeom::algo::transformed(&mut model, &f, moved)
                .unwrap()
                .shape
        } else {
            disc(
                &mut model,
                Frame::new(Point::new(10.0, 10.0, 0.0), Direction::X, Direction::Y, T).unwrap(),
            )
        };
        let pipe = ogeom::offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
            .unwrap_or_else(|e| panic!("placed {placed}: {e}"))
            .shape;
        let (v, exact) = measure(&model, &pipe);
        let want = 20.0 * core::f64::consts::PI;
        assert!(exact, "measured from a mesh");
        assert!(
            (v - want).abs() < want * 1e-9,
            "placed {placed}: {v} against {want}"
        );
    }
}

/// A square down a spine of three straight legs turning out of plane:
/// mitred exactly at both corners, its volume the area times the length.
#[test]
fn a_square_down_a_skew_polyline_is_mitred_exactly() {
    let mut model = Model::new();
    let profile = square(&mut model, 1.0);
    let spine = make_polygon(
        &mut model,
        &[
            Point::ORIGIN,
            Point::new(0.0, 0.0, 12.0),
            Point::new(10.0, 0.0, 12.0),
            Point::new(10.0, 8.0, 12.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let pipe = ogeom::offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &pipe, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let (v, exact) = measure(&model, &pipe);
    let want = 4.0 * 30.0;
    assert!(exact, "measured from a mesh");
    assert!((v - want).abs() < want * 1e-9, "{v} against {want}");
}

/// Two quarter arcs turning opposite ways meet at a corner, a square
/// swept down them: the flat strips of the second leg start on the
/// corner's join row, which lies along the leg's own chord, and still walk
/// every edge once each way. Valid, and the volume the union of the two
/// annuli and the corner fill makes.
#[test]
fn a_square_round_a_corner_between_curved_legs_walks_each_edge_once_each_way() {
    let mut model = Model::new();
    let r = 20.0;
    let vertex = |model: &mut Model, p: Point| ogeom::algo::make_vertex(model, p).shape;
    let (va, vb, vd) = (
        vertex(&mut model, Point::new(r, 0.0, 0.0)),
        vertex(&mut model, Point::new(0.0, r, 0.0)),
        vertex(&mut model, Point::new(r, 2.0 * r, 0.0)),
    );
    let arc = |model: &mut Model, centre: Point, range: (f64, f64), from: &Shape, to: &Shape| {
        let frame = Frame::new(centre, Direction::Z, Direction::X, T).unwrap();
        let curve = CircleCurve::new(Circle::new(frame, r, T).unwrap()).into();
        ogeom::algo::make_edge_between(model, curve, range, from, to, T)
            .unwrap()
            .shape
    };
    let half_pi = core::f64::consts::FRAC_PI_2;
    let first = arc(&mut model, Point::ORIGIN, (0.0, half_pi), &va, &vb);
    // About (r, r) from b back to d: walked against its parameter, turning
    // the other way.
    let second = arc(
        &mut model,
        Point::new(r, r, 0.0),
        (half_pi, core::f64::consts::PI),
        &vd,
        &vb,
    );
    let spine = make_wire(&mut model, &[first, second.reversed()], T)
        .unwrap()
        .shape;
    // A 4 x 4 square square to the spine's start, centred on it.
    let w = 2.0;
    let pts = [(-w, -w), (w, -w), (w, w), (-w, w)].map(|(a, b)| Point::new(r + a, 0.0, b));
    let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
    let plane =
        Plane::new(Frame::new(Point::new(r, 0.0, 0.0), Direction::Y, Direction::X, T).unwrap());
    let profile = make_face(&mut model, PlaneSurface::new(plane).into(), &[wire], T)
        .unwrap()
        .shape;
    let pipe = ogeom::offset::make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    assert_eq!(walks::edges_walked_one_way(&model, &pipe), 0);
    let found = ogeom::algo::check(&model, &pipe, T).unwrap();
    assert!(found.is_valid(), "{found}");
    // The plane slice's width at height y: the first annulus's span
    // (x >= 0), the second's (x <= r), and the corner fill, merged.
    let width = |y: f64| -> f64 {
        let (ro, ri) = (r + w, r - w);
        let mut spans: Vec<(f64, f64)> = Vec::new();
        if (0.0..=ro).contains(&y) {
            spans.push(((ri * ri - y * y).max(0.0).sqrt(), (ro * ro - y * y).sqrt()));
        }
        let e = y - r;
        if (0.0..=ro).contains(&e) {
            spans.push((
                r - (ro * ro - e * e).sqrt(),
                r - (ri * ri - e * e).max(0.0).sqrt(),
            ));
        }
        if (r - w..=r).contains(&y) {
            spans.push((-w, 0.0));
        }
        spans.sort_by(|p, q| p.0.total_cmp(&q.0));
        let (mut total, mut reach) = (0.0, f64::NEG_INFINITY);
        for (lo, hi) in spans {
            let lo = lo.max(reach);
            if hi > lo {
                total += hi - lo;
                reach = hi;
            }
        }
        total
    };
    let want = simpson(0.0, 2.0 * r + w, 40000, width) * 2.0 * w;
    let fine = Deflection {
        chord: 1e-3,
        ..Deflection::default()
    };
    let v = volume_properties(&model, &pipe, fine, T).unwrap().mass;
    assert!((v - want).abs() < want * 5e-3, "{v} against {want}");
}
