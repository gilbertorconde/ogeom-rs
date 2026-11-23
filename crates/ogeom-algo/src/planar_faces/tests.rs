use std::f64::consts::{FRAC_1_SQRT_2, PI, TAU};

use ogeom_geom::{BSplineCurve, CircleCurve, Curve, EllipseCurve, LineCurve};
use ogeom_math::{Circle, Direction, Ellipse, KnotVector, Vector, Weighted};
use ogeom_mesh::Deflection;

use super::*;
use crate::build::{make_edge, make_edge_between, make_polygon, make_vertex, make_wire};

fn tol() -> Tolerances {
    Tolerances::default()
}

fn p(x: f64, y: f64) -> Point {
    Point::new(x, y, 0.0)
}

/// A full circle in the `xy` plane, counterclockwise about `+z` when `ccw`.
fn circle_wire(m: &mut Model, centre: Point, r: f64, ccw: bool) -> Shape {
    let z = if ccw { Direction::Z } else { -Direction::Z };
    let circle = Circle::new(Frame::about(centre, z), r, tol()).unwrap();
    let edge = make_edge(m, CircleCurve::new(circle).into(), (0.0, TAU), tol())
        .unwrap()
        .shape;
    make_wire(m, &[edge], tol()).unwrap().shape
}

/// A closed polygon through `points` in order.
fn polygon(m: &mut Model, points: &[Point]) -> Shape {
    make_polygon(m, points, true, tol()).unwrap().shape
}

/// An axis-aligned square of side `side` about `centre`.
fn square_wire(m: &mut Model, centre: Point, side: f64, ccw: bool) -> Shape {
    let h = side / 2.0;
    let mut corners = vec![
        Point::new(centre.x - h, centre.y - h, centre.z),
        Point::new(centre.x + h, centre.y - h, centre.z),
        Point::new(centre.x + h, centre.y + h, centre.z),
        Point::new(centre.x - h, centre.y + h, centre.z),
    ];
    if !ccw {
        corners.reverse();
    }
    polygon(m, &corners)
}

/// A wire through the given curves in turn, each over its range, each
/// ending where the next starts and the last where the first starts.
fn chained(m: &mut Model, pieces: Vec<(Curve, (f64, f64))>) -> Shape {
    let starts: Vec<Shape> = pieces
        .iter()
        .map(|(c, r)| make_vertex(m, c.point_at(r.0, tol()).unwrap()).shape)
        .collect();
    let n = pieces.len();
    let edges: Vec<Shape> = pieces
        .into_iter()
        .enumerate()
        .map(|(k, (c, r))| {
            make_edge_between(m, c, r, &starts[k], &starts[(k + 1) % n], tol())
                .unwrap()
                .shape
        })
        .collect();
    make_wire(m, &edges, tol()).unwrap().shape
}

fn segment(a: Point, b: Point) -> (Curve, (f64, f64)) {
    (
        LineCurve::segment(a, b, tol()).unwrap().into(),
        (0.0, a.distance(b)),
    )
}

fn area(m: &Model, face: &Shape) -> f64 {
    crate::surface_properties(m, face, Deflection::default(), tol())
        .unwrap()
        .mass
}

/// Every face valid, facing `normal`, its outer wire counterclockwise
/// about it and its holes clockwise; returns the faces' areas.
fn checked(m: &Model, built: &PlanarFaces, normal: Vector) -> Vec<f64> {
    assert_eq!(built.faces.len(), built.wire_of.len());
    built
        .faces
        .iter()
        .map(|face| {
            let diagnosis = crate::check(m, face, tol()).unwrap();
            assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
            let (_, n) = crate::face_normal(m, face, tol()).unwrap();
            assert!(
                (n.normalized(tol()).unwrap() - normal).magnitude() < 1e-12,
                "face normal {n:?} against {normal:?}"
            );
            for (k, wire) in m.ordered_children_of(face).unwrap().iter().enumerate() {
                let turn = crate::build::wire_turn(m, wire, normal, tol()).unwrap();
                assert_eq!(turn > 0.0, k == 0, "wire {k} turns {turn}");
            }
            area(m, face)
        })
        .collect()
}

fn close(measured: f64, expected: f64) {
    let bound = 1e-9 * expected.abs().max(1.0);
    assert!(
        (measured - expected).abs() <= bound,
        "measured {measured}, expected {expected} (off {:.3e}, bound {bound:.1e})",
        (measured - expected).abs()
    );
}

fn refusal(result: OgeomResult<PlanarFaces>) -> String {
    match result {
        Ok(built) => panic!("expected a refusal, got {:?}", built.wire_of),
        Err(e) => e.to_string(),
    }
}

#[test]
fn an_annulus_with_a_square_hole_and_an_island_disc() {
    let mut m = Model::new();
    let outer = circle_wire(&mut m, p(0.0, 0.0), 20.0, true);
    let hole = square_wire(&mut m, p(0.0, 0.0), 10.0, true);
    let island = circle_wire(&mut m, p(0.0, 0.0), 3.0, true);
    let built = faces_from_planar_wires(&mut m, Plane::XY, &[outer, hole, island], tol()).unwrap();
    assert_eq!(built.wire_of, vec![vec![0, 1], vec![2]]);
    let areas = checked(&m, &built, Vector::Z);
    close(areas[0], PI * 400.0 - 100.0);
    close(areas[1], PI * 9.0);
}

#[test]
fn wires_given_clockwise_and_shuffled_bound_the_same_faces() {
    // Each draw: a permutation and which wires run clockwise.
    let mut seed: u64 = 0x2545_f491_4f6c_dd1d;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    for _ in 0..8 {
        let mut order = [0_usize, 1, 2, 3];
        for k in (1..order.len()).rev() {
            order.swap(
                k,
                usize::try_from(next() % (u64::try_from(k).unwrap() + 1)).unwrap(),
            );
        }
        let flips = next();
        let mut m = Model::new();
        let wires: Vec<Shape> = order
            .iter()
            .enumerate()
            .map(|(slot, &which)| {
                let ccw = flips >> slot & 1 == 0;
                match which {
                    0 => square_wire(&mut m, p(0.0, 0.0), 40.0, ccw),
                    1 => circle_wire(&mut m, p(0.0, 0.0), 15.0, ccw),
                    2 => square_wire(&mut m, p(0.0, 0.0), 10.0, ccw),
                    _ => circle_wire(&mut m, p(50.0, 0.0), 4.0, ccw),
                }
            })
            .collect();
        let built = faces_from_planar_wires(&mut m, Plane::XY, &wires, tol()).unwrap();
        let slot = |which: usize| order.iter().position(|&w| w == which).unwrap();
        let mut expected = vec![vec![slot(0), slot(1)], vec![slot(2)], vec![slot(3)]];
        expected.sort();
        assert_eq!(built.wire_of, expected, "order {order:?}");
        let areas = checked(&m, &built, Vector::Z);
        for (bounding, measured) in built.wire_of.iter().zip(areas) {
            let want = match order[bounding[0]] {
                0 => 1600.0 - PI * 225.0,
                2 => 100.0,
                _ => PI * 16.0,
            };
            close(measured, want);
        }
    }
}

#[test]
fn a_square_with_a_square_hole() {
    let mut m = Model::new();
    let hole = square_wire(&mut m, p(1.0, 2.0), 4.0, true);
    let outer = square_wire(&mut m, p(0.0, 0.0), 10.0, false);
    let built = faces_from_planar_wires(&mut m, Plane::XY, &[hole, outer], tol()).unwrap();
    assert_eq!(built.wire_of, vec![vec![1, 0]]);
    close(checked(&m, &built, Vector::Z)[0], 84.0);
}

#[test]
fn islands_alternate_with_holes_down_a_nest() {
    let mut m = Model::new();
    let radii = [10.0, 8.0, 6.0, 4.0, 2.0];
    let wires: Vec<Shape> = radii
        .iter()
        .rev()
        .map(|&r| circle_wire(&mut m, p(0.0, 0.0), r, r > 5.0))
        .collect();
    let built = faces_from_planar_wires(&mut m, Plane::XY, &wires, tol()).unwrap();
    // Given innermost first: 2, 4, 6, 8, 10.
    assert_eq!(built.wire_of, vec![vec![0], vec![2, 1], vec![4, 3]]);
    let areas = checked(&m, &built, Vector::Z);
    close(areas[0], PI * 4.0);
    close(areas[1], PI * (36.0 - 16.0));
    close(areas[2], PI * (100.0 - 64.0));
}

#[test]
fn a_ring_of_circles_is_classified_in_one_sweep() {
    const N: u32 = 200;
    let mut m = Model::new();
    let mut wires = vec![circle_wire(&mut m, p(0.0, 0.0), 120.0, true)];
    for k in 0..N {
        let a = TAU * f64::from(k) / f64::from(N);
        wires.push(circle_wire(
            &mut m,
            p(100.0 * a.cos(), 100.0 * a.sin()),
            1.0,
            k % 2 == 0,
        ));
    }
    let (built, tests) = arranged(&mut m, Plane::XY, &wires, tol()).unwrap();
    let all: Vec<usize> = (0..=N).map(|k| k as usize).collect();
    assert_eq!(built.wire_of, vec![all]);
    close(
        checked(&m, &built, Vector::Z)[0],
        PI * (120.0 * 120.0 - f64::from(N)),
    );
    // One point per small circle against the big one; neighbours 3.14
    // apart never share a box.
    assert!(tests <= N as usize, "{tests} points classified");

    // The same ring with no circle around it: a face per circle.
    let mut m = Model::new();
    let ring: Vec<Shape> = (0..N)
        .map(|k| {
            let a = TAU * f64::from(k) / f64::from(N);
            circle_wire(&mut m, p(100.0 * a.cos(), 100.0 * a.sin()), 1.0, true)
        })
        .collect();
    let (built, tests) = arranged(&mut m, Plane::XY, &ring, tol()).unwrap();
    assert_eq!(built.faces.len(), N as usize);
    assert_eq!(tests, 0);
}

#[test]
fn arcs_ellipses_and_b_splines_bound_exact_areas() {
    let mut m = Model::new();
    // A slot: two half circles of radius 2 joined by segments 6 long.
    let half = |centre: Point, from: f64| -> (Curve, (f64, f64)) {
        let frame =
            Frame::from_axes(centre, Direction::X, Direction::Y, Direction::Z, tol()).unwrap();
        let circle = Circle::new(frame, 2.0, tol()).unwrap();
        (CircleCurve::new(circle).into(), (from, from + PI))
    };
    let slot = chained(
        &mut m,
        vec![
            segment(p(0.0, -2.0), p(6.0, -2.0)),
            half(p(6.0, 0.0), -PI / 2.0),
            segment(p(6.0, 2.0), p(0.0, 2.0)),
            half(p(0.0, 0.0), PI / 2.0),
        ],
    );
    // A parabolic cap y = h(1 - x^2) on [-1, 1] as a quadratic Bezier,
    // closed by its chord: 4h/3.
    let h = 1.5;
    let (dx, dy) = (20.0, 0.0);
    let cap_curve: Curve = BSplineCurve::new(
        KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap(),
        vec![p(dx - 1.0, dy), p(dx, dy + 2.0 * h), p(dx + 1.0, dy)],
        tol(),
    )
    .unwrap()
    .into();
    let cap = chained(
        &mut m,
        vec![
            (cap_curve, (0.0, 1.0)),
            segment(p(dx + 1.0, dy), p(dx - 1.0, dy)),
        ],
    );
    // A quarter disc of radius 3 whose arc is a rational quadratic: 9 pi/4.
    let r = 3.0;
    let o = p(30.0, 0.0);
    let arc: Curve = BSplineCurve::rational(
        KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap(),
        vec![
            Weighted::new(p(o.x + r, o.y), 1.0, tol()).unwrap(),
            Weighted::new(p(o.x + r, o.y + r), FRAC_1_SQRT_2, tol()).unwrap(),
            Weighted::new(p(o.x, o.y + r), 1.0, tol()).unwrap(),
        ],
    )
    .unwrap()
    .into();
    let quarter = chained(
        &mut m,
        vec![
            (arc, (0.0, 1.0)),
            segment(p(o.x, o.y + r), o),
            segment(o, p(o.x + r, o.y)),
        ],
    );
    // An ellipse 5 by 4 around the quarter disc, clockwise: a hole.
    let clockwise = Frame::from_axes(
        p(31.0, 1.0),
        Direction::X,
        -Direction::Y,
        -Direction::Z,
        tol(),
    )
    .unwrap();
    let ellipse = Ellipse::new(clockwise, 5.0, 4.0, tol()).unwrap();
    let ring_edge = make_edge(&mut m, EllipseCurve::new(ellipse).into(), (0.0, TAU), tol())
        .unwrap()
        .shape;
    let around = make_wire(&mut m, &[ring_edge], tol()).unwrap().shape;
    // A square round the cap.
    let frame = square_wire(&mut m, p(dx, dy + 0.5), 4.0, true);

    let built = faces_from_planar_wires(
        &mut m,
        Plane::XY,
        &[slot, cap, quarter, around, frame],
        tol(),
    )
    .unwrap();
    // The quarter disc is a hole in the ellipse, the cap one in the square.
    assert_eq!(built.wire_of, vec![vec![0], vec![3, 2], vec![4, 1]]);
    let areas = checked(&m, &built, Vector::Z);
    close(areas[0], 4.0 * 6.0 + PI * 4.0);
    close(areas[1], PI * 20.0 - PI * r * r / 4.0);
    close(areas[2], 16.0 - 4.0 * h / 3.0);
}

#[test]
fn wires_of_different_faces_may_touch() {
    let mut m = Model::new();
    // Two squares meeting at one corner.
    let a = square_wire(&mut m, p(0.0, 0.0), 2.0, true);
    let b = square_wire(&mut m, p(2.0, 2.0), 2.0, false);
    // An annulus, and an island in its hole tangent to the hole at (30, 0).
    let outer = circle_wire(&mut m, p(20.0, 0.0), 15.0, true);
    let hole = circle_wire(&mut m, p(20.0, 0.0), 10.0, true);
    let island = circle_wire(&mut m, p(25.0, 0.0), 5.0, false);
    // A disc outside the annulus, tangent to it at (35, 0).
    let beside = circle_wire(&mut m, p(38.0, 0.0), 3.0, false);
    let built = faces_from_planar_wires(
        &mut m,
        Plane::XY,
        &[a, b, outer, hole, island, beside],
        tol(),
    )
    .unwrap();
    assert_eq!(
        built.wire_of,
        vec![vec![0], vec![1], vec![2, 3], vec![4], vec![5]]
    );
    let areas = checked(&m, &built, Vector::Z);
    close(areas[0], 4.0);
    close(areas[1], 4.0);
    close(areas[2], PI * 125.0);
    close(areas[3], PI * 25.0);
    close(areas[4], PI * 9.0);
}

#[test]
fn wires_of_one_face_touching_are_refused_by_name() {
    // A hole tangent to its outer boundary.
    let mut m = Model::new();
    let lone = square_wire(&mut m, p(-40.0, 0.0), 1.0, true);
    let outer = circle_wire(&mut m, p(0.0, 0.0), 10.0, true);
    let hole = circle_wire(&mut m, p(5.0, 0.0), 5.0, true);
    let message = refusal(faces_from_planar_wires(
        &mut m,
        Plane::XY,
        &[lone, outer, hole],
        tol(),
    ));
    assert!(message.contains("wires 1 and 2 touch"), "{message}");

    // A diamond whose corners sit on its square's sides.
    let mut m = Model::new();
    let square = square_wire(&mut m, p(0.0, 0.0), 2.0, true);
    let diamond = polygon(
        &mut m,
        &[p(1.0, 0.0), p(0.0, 1.0), p(-1.0, 0.0), p(0.0, -1.0)],
    );
    let message = refusal(faces_from_planar_wires(
        &mut m,
        Plane::XY,
        &[diamond, square],
        tol(),
    ));
    assert!(message.contains("wires 0 and 1 touch"), "{message}");

    // Two holes of one face touching each other.
    let mut m = Model::new();
    let outer = square_wire(&mut m, p(0.0, 0.0), 20.0, true);
    let left = circle_wire(&mut m, p(-2.0, 0.0), 2.0, true);
    let right = circle_wire(&mut m, p(2.0, 0.0), 2.0, true);
    let message = refusal(faces_from_planar_wires(
        &mut m,
        Plane::XY,
        &[outer, left, right],
        tol(),
    ));
    assert!(message.contains("wires 1 and 2 touch"), "{message}");
}

#[test]
fn crossing_wires_are_refused_by_name() {
    let mut m = Model::new();
    let lone = circle_wire(&mut m, p(-40.0, 0.0), 1.0, true);
    let a = circle_wire(&mut m, p(0.0, 0.0), 5.0, true);
    let b = square_wire(&mut m, p(6.0, 0.0), 6.0, true);
    let message = refusal(faces_from_planar_wires(
        &mut m,
        Plane::XY,
        &[lone, a, b],
        tol(),
    ));
    assert!(message.contains("wires 1 and 2 cross"), "{message}");

    // A circle round a square's corner.
    let mut m = Model::new();
    let s = square_wire(&mut m, p(0.0, 0.0), 2.0, true);
    let c = circle_wire(&mut m, p(1.0, 1.0), 1.0, true);
    let message = refusal(faces_from_planar_wires(&mut m, Plane::XY, &[s, c], tol()));
    assert!(message.contains("wires 0 and 1 cross"), "{message}");
}

#[test]
fn wires_sharing_a_stretch_are_refused_by_name() {
    let mut m = Model::new();
    let a = square_wire(&mut m, p(0.0, 0.0), 2.0, true);
    let b = square_wire(&mut m, p(2.0, 0.0), 2.0, true);
    let message = refusal(faces_from_planar_wires(&mut m, Plane::XY, &[a, b], tol()));
    assert!(
        message.contains("wires 0 and 1 run along each other"),
        "{message}"
    );

    let mut m = Model::new();
    let a = circle_wire(&mut m, p(0.0, 0.0), 2.0, true);
    let b = circle_wire(&mut m, p(0.0, 0.0), 2.0, false);
    let message = refusal(faces_from_planar_wires(&mut m, Plane::XY, &[a, b], tol()));
    assert!(
        message.contains("wires 0 and 1 run along each other"),
        "{message}"
    );
}

#[test]
fn open_crossed_and_off_plane_wires_are_refused_by_name() {
    let mut m = Model::new();
    let fine = square_wire(&mut m, p(0.0, 0.0), 2.0, true);
    let open = make_polygon(
        &mut m,
        &[p(5.0, 0.0), p(6.0, 0.0), p(6.0, 1.0)],
        false,
        tol(),
    )
    .unwrap()
    .shape;
    let message = refusal(faces_from_planar_wires(
        &mut m,
        Plane::XY,
        &[fine.clone(), open],
        tol(),
    ));
    assert!(message.contains("wire 1 is open"), "{message}");

    let bow_tie = polygon(
        &mut m,
        &[p(5.0, 0.0), p(7.0, 2.0), p(7.0, 0.0), p(5.0, 2.0)],
    );
    let message = refusal(faces_from_planar_wires(
        &mut m,
        Plane::XY,
        &[fine.clone(), bow_tie],
        tol(),
    ));
    assert!(message.contains("wire 1 crosses itself"), "{message}");

    let lifted = circle_wire(&mut m, Point::new(5.0, 0.0, 1.0), 1.0, true);
    let message = refusal(faces_from_planar_wires(
        &mut m,
        Plane::XY,
        &[fine.clone(), lifted],
        tol(),
    ));
    assert!(
        message.contains("wire 1 does not lie in the plane"),
        "{message}"
    );

    let tilted = polygon(
        &mut m,
        &[
            Point::new(5.0, 0.0, 0.0),
            Point::new(6.0, 0.0, 0.0),
            Point::new(6.0, 1.0, 0.1),
        ],
    );
    let message = refusal(faces_from_planar_wires(
        &mut m,
        Plane::XY,
        &[fine, tilted],
        tol(),
    ));
    assert!(
        message.contains("wire 1 does not lie in the plane"),
        "{message}"
    );
}

fn left_handed(plane: Plane) -> Plane {
    let f = plane.frame();
    Plane::new(Frame::from_axes(f.origin(), f.x(), -f.y(), f.z(), tol()).unwrap())
}

#[test]
fn faces_face_the_plane_s_normal_however_it_points() {
    let normal = Vector::new(1.0, -2.0, 0.5).normalized(tol()).unwrap();
    let origin = Point::new(3.0, -1.0, 7.0);
    for plane in [
        Plane::through(origin, Direction::new(normal, tol()).unwrap()),
        Plane::through(origin, Direction::new(-normal, tol()).unwrap()),
        left_handed(Plane::through(
            origin,
            Direction::new(normal, tol()).unwrap(),
        )),
    ] {
        let frame = plane.frame();
        let at = |x: f64, y: f64| frame.origin() + frame.x().vector() * x + frame.y().vector() * y;
        let mut m = Model::new();
        let outer = polygon(
            &mut m,
            &[at(0.0, 0.0), at(6.0, 0.0), at(6.0, 4.0), at(0.0, 4.0)],
        );
        let hole = polygon(
            &mut m,
            &[at(1.0, 1.0), at(2.0, 1.0), at(2.0, 3.0), at(1.0, 3.0)],
        );
        let built = faces_from_planar_wires(&mut m, plane, &[outer, hole], tol()).unwrap();
        assert_eq!(built.wire_of, vec![vec![0, 1]]);
        close(checked(&m, &built, frame.z().vector())[0], 22.0);
    }
}

#[test]
fn no_wires_make_no_faces() {
    let mut m = Model::new();
    let built = faces_from_planar_wires(&mut m, Plane::XY, &[], tol()).unwrap();
    assert!(built.faces.is_empty() && built.wire_of.is_empty());
}
