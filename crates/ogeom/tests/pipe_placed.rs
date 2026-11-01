//! Circular pipes along wires of several edges, and profiles brought onto
//! their spine before a pipe shell sweeps them, measured against closed
//! forms.
//!
//! A disc of radius `r` square to a spine, its centre on it, holds by
//! Pappus `π r²` times the spine's length, mitred corners included: each
//! leg is a cylinder cut on planes through its axis's ends.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, find_plane, make_edge, make_edge_between, make_face, make_polygon, make_vertex,
    make_wire, surface_properties, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{BSplineCurve, CircleCurve, Curve3d as _, LineCurve, PlaneSurface};
use ogeom::math::{Circle, Direction, Frame, KnotVector, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::offset::{
    PipeContact, PipeCorners, PipeLaw, PipeLaws, PipePlacement, make_pipe, make_pipe_shell,
    make_pipe_shell_laws, make_pipe_skinned, place_on_spine,
};
use ogeom::topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// The volume and centroid, and whether they were integrated exactly.
fn measure(model: &Model, shape: &Shape) -> (f64, Point, bool) {
    let p = volume_properties(model, shape, Deflection::default(), T).unwrap();
    (p.mass, p.centre, p.deflection == 0.0)
}

/// Up `z` to `(0, 0, 5)`, a quarter arc of radius 3 over to `(3, 0, 8)`
/// heading `+x`, and on along `x` to `(8, 0, 8)`, each edge tangent to the
/// next; then, where `turn`, a sharp quarter turn along `y` to
/// `(8, 6, 8)`.
fn line_arc_line(model: &mut Model, turn: bool) -> Shape {
    let points = [
        Point::ORIGIN,
        Point::new(0.0, 0.0, 5.0),
        Point::new(3.0, 0.0, 8.0),
        Point::new(8.0, 0.0, 8.0),
        Point::new(8.0, 6.0, 8.0),
    ];
    let v: Vec<Shape> = points
        .iter()
        .map(|p| make_vertex(model, *p).shape)
        .collect();
    let segment = |model: &mut Model, i: usize| {
        let s = LineCurve::segment(points[i], points[i + 1], T).unwrap();
        let range = s.domain();
        make_edge_between(model, s.into(), range, &v[i], &v[i + 1], T)
            .unwrap()
            .shape
    };
    let frame = Frame::new(Point::new(3.0, 0.0, 5.0), Direction::Y, -Direction::X, T).unwrap();
    let arc = CircleCurve::new(Circle::new(frame, 3.0, T).unwrap());
    let mut edges = vec![segment(model, 0)];
    edges.push(
        make_edge_between(model, arc.into(), (0.0, PI / 2.0), &v[1], &v[2], T)
            .unwrap()
            .shape,
    );
    edges.push(segment(model, 2));
    if turn {
        edges.push(segment(model, 3));
    }
    make_wire(model, &edges, T).unwrap().shape
}

fn analytic_only(model: &Model, shape: &Shape) -> bool {
    use ogeom::geom::SurfaceGeometry as S;
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .iter()
        .all(|f| {
            let Some(NodeData::Face(d)) = model.node(f).map(|n| n.data()) else {
                return false;
            };
            matches!(
                model.geometry().surface(d.surface),
                Some(S::Plane(_) | S::Cylinder(_) | S::Torus(_))
            )
        })
}

/// A pipe of radius `r` along three straight legs, 10, 10 and 7 long,
/// mitred at both corners: `π r² · 27`, exactly. Along a line, a quarter
/// arc of radius 3 and a line, meeting tangent: `π r² (10 + 3π/2)`,
/// exactly, by Pappus.
#[test]
fn a_pipe_along_a_wire_holds_pappus_s_volume_exactly() {
    let r = 0.5;
    let mut model = Model::new();
    let legs = make_polygon(
        &mut model,
        &[
            Point::ORIGIN,
            Point::new(0.0, 0.0, 10.0),
            Point::new(10.0, 0.0, 10.0),
            Point::new(10.0, 7.0, 10.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let smooth = line_arc_line(&mut model, false);
    for (spine, length) in [(legs, 27.0), (smooth, 10.0 + 1.5 * PI)] {
        let pipe = make_pipe(&mut model, &spine, r, T).unwrap().shape;
        let found = check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "{found}");
        assert!(analytic_only(&model, &pipe));
        let (v, _, exact) = measure(&model, &pipe);
        assert!(exact);
        let want = PI * r * r * length;
        assert!((v - want).abs() < want * 1e-12, "{v} against {want}");
    }
}

/// The rounded rectangle in the `xy` plane about the origin: sides `2 a`
/// along `x` and `2 b` along `y` between quarter arcs of radius `rho`
/// centred at `(±a, ±b, 0)`, each meeting the next tangent to it, from
/// `(a + rho, -b, 0)` heading `+y`. With `a = 0` it is the stadium, its
/// two half circles each two quarter arcs. Its length is
/// `4 a + 4 b + 2 π rho`.
fn rounded(model: &mut Model, a: f64, b: f64, rho: f64) -> Shape {
    let centres = [(a, b), (-a, b), (-a, -b), (a, -b)].map(|(x, y)| Point::new(x, y, 0.0));
    let mut points = Vec::new();
    for (k, centre) in centres.iter().enumerate() {
        let from = f64::from(u8::try_from(k).unwrap()) * PI / 2.0;
        let to = from + PI / 2.0;
        points.push(*centre + Vector::new(rho * from.cos(), rho * from.sin(), 0.0));
        points.push(*centre + Vector::new(rho * to.cos(), rho * to.sin(), 0.0));
    }
    let v: Vec<Shape> = points
        .iter()
        .map(|p| make_vertex(model, *p).shape)
        .collect();
    let mut edges = Vec::new();
    for (k, centre) in centres.iter().enumerate() {
        let before = (2 * k + 7) % 8;
        // Where the side has no length the arc starts where the last ends.
        let start = if points[before].distance(points[2 * k]) > T.confusion() {
            let s = LineCurve::segment(points[before], points[2 * k], T).unwrap();
            let range = s.domain();
            edges.push(
                make_edge_between(model, s.into(), range, &v[before], &v[2 * k], T)
                    .unwrap()
                    .shape,
            );
            &v[2 * k]
        } else {
            &v[before]
        };
        let frame = Frame::new(*centre, Direction::Z, Direction::X, T).unwrap();
        let arc = CircleCurve::new(Circle::new(frame, rho, T).unwrap());
        let from = f64::from(u8::try_from(k).unwrap()) * PI / 2.0;
        edges.push(
            make_edge_between(
                model,
                arc.into(),
                (from, from + PI / 2.0),
                start,
                &v[(2 * k + 1) % 8],
                T,
            )
            .unwrap()
            .shape,
        );
    }
    make_wire(model, &edges, T).unwrap().shape
}

/// A pipe of radius `r` round a closed spine of lines and arcs meeting
/// tangent is a closed ring: no cap, every wall a cylinder or a torus,
/// holding `π r²` times the spine's length exactly by Pappus, its centroid
/// the spine's centre of symmetry. The stadium of two lines 4 long and
/// half circles of radius 3 is `8 + 6π` long; the rounded rectangle with
/// sides 6 and 4 and quarter arcs of radius 1.5 is `20 + 3π`.
#[test]
fn a_pipe_round_a_closed_wire_of_lines_and_arcs_is_a_ring_of_pappus_s_volume() {
    let r = 0.5;
    let mut model = Model::new();
    let stadium = rounded(&mut model, 0.0, 2.0, 3.0);
    let rectangle = rounded(&mut model, 3.0, 2.0, 1.5);
    for (spine, length) in [(stadium, 8.0 + 6.0 * PI), (rectangle, 20.0 + 3.0 * PI)] {
        let pipe = make_pipe(&mut model, &spine, r, T).unwrap().shape;
        let found = check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "{found}");
        assert!(analytic_only(&model, &pipe));
        let faces = explore_unique(&model, &pipe, ShapeType::Face).unwrap();
        assert!(
            faces.iter().all(|f| {
                let Some(NodeData::Face(d)) = model.node(f).map(|n| n.data()) else {
                    return false;
                };
                !matches!(
                    model.geometry().surface(d.surface),
                    Some(ogeom::geom::SurfaceGeometry::Plane(_))
                )
            }),
            "a closed ring has no cap"
        );
        let (v, centre, exact) = measure(&model, &pipe);
        assert!(exact);
        // Integrated exactly: the error is the arithmetic's.
        let want = PI * r * r * length;
        assert!((v - want).abs() < want * 1e-12, "{v} against {want}");
        assert!(centre.distance(Point::ORIGIN) < 1e-9, "{centre:?}");
    }
}

/// Up `z` to `(0, 0, 5)`, then sharply onto a quarter arc of radius 3 in
/// the plane `z = 5` about `(0, 3, 5)`, heading `+x` at its start: a
/// corner that turns the arc out of its own plane.
fn line_corner_flat_arc(model: &mut Model) -> Shape {
    let (a, b) = (Point::ORIGIN, Point::new(0.0, 0.0, 5.0));
    let (va, vb) = (make_vertex(model, a).shape, make_vertex(model, b).shape);
    let segment = LineCurve::segment(a, b, T).unwrap();
    let range = segment.domain();
    let up = make_edge_between(model, segment.into(), range, &va, &vb, T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(0.0, 3.0, 5.0), Direction::Z, -Direction::Y, T).unwrap();
    let arc = CircleCurve::new(Circle::new(frame, 3.0, T).unwrap());
    let vc = make_vertex(model, arc.point_at(PI / 2.0, T).unwrap()).shape;
    let round = make_edge_between(model, arc.into(), (0.0, PI / 2.0), &vb, &vc, T)
        .unwrap()
        .shape;
    make_wire(model, &[up, round], T).unwrap().shape
}

/// The disc of radius `r` about `at`, square to `heading`.
fn disc(model: &mut Model, at: Point, heading: Direction, r: f64) -> Shape {
    let frame = Frame::about(at, heading);
    let circle = CircleCurve::new(Circle::new(frame, r, T).unwrap());
    let rim = make_edge(model, circle.into(), (0.0, 2.0 * PI), T)
        .unwrap()
        .shape;
    let ring = make_wire(model, &[rim], T).unwrap().shape;
    make_face(
        model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[ring],
        T,
    )
    .unwrap()
    .shape
}

/// A sharp corner between a line and an arc is mitred as the pipe shell
/// of the disc mitres it, the default corners: here exactly, each side run
/// on straight past the corner and cut on the plane halving it. The pipe
/// is that solid: the same walls, the same volume.
#[test]
fn a_pipe_round_a_sharp_corner_onto_an_arc_is_the_pipe_shell_of_its_disc() {
    let r = 0.5;
    let mut model = Model::new();
    let spine = line_corner_flat_arc(&mut model);
    let pipe = make_pipe(&mut model, &spine, r, T).unwrap().shape;
    let found = check(&model, &pipe, T).unwrap();
    assert!(found.is_valid(), "{found}");
    assert!(analytic_only(&model, &pipe));
    let profile = disc(&mut model, Point::ORIGIN, Direction::Z, r);
    let shell = make_pipe_shell(&mut model, &profile, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    let count = |shape: &Shape| {
        explore_unique(&model, shape, ShapeType::Face)
            .unwrap()
            .len()
    };
    assert_eq!(count(&pipe), count(&shell));
    let ((v, c, exact), (w, d, exact_shell)) = (measure(&model, &pipe), measure(&model, &shell));
    assert!(exact && exact_shell);
    // Both integrated exactly: the error is the arithmetic's.
    assert!((v - w).abs() < w * 1e-12, "{v} against {w}");
    assert!(c.distance(d) < 1e-9, "{c:?} against {d:?}");
}

/// What the pipe shell sweeps only as a skin, or does not sweep, is
/// refused by name, with the function that takes a tolerance for it: a
/// free-form edge. A tube wider than an arc bends, and a leg shorter than
/// the reach of the crossing inside its corner, are refused as before.
#[test]
fn a_pipe_along_a_wire_refuses_what_it_does_not_sweep_exactly_by_name() {
    let mut model = Model::new();
    let (a, b) = (Point::ORIGIN, Point::new(0.0, 0.0, 5.0));
    let (va, vb) = (
        make_vertex(&mut model, a).shape,
        make_vertex(&mut model, b).shape,
    );
    let segment = LineCurve::segment(a, b, T).unwrap();
    let range = segment.domain();
    let up = make_edge_between(&mut model, segment.into(), range, &va, &vb, T)
        .unwrap()
        .shape;
    let free = free_form_after_line(&mut model, &up, &vb);
    let err = make_pipe(&mut model, &free, 0.5, T).unwrap_err();
    assert!(err.to_string().contains("make_pipe_skinned"), "{err}");
    let tight = line_arc_line(&mut model, false);
    let err = make_pipe(&mut model, &tight, 3.0, T).unwrap_err();
    assert!(err.to_string().contains("swallows its spine"), "{err}");
    let short = planar_path(
        &mut model,
        &[
            Step::Line(0.3),
            Step::Turn(PI / 2.0),
            Step::Arc(3.0, PI / 2.0),
        ],
    );
    let err = make_pipe(&mut model, &short, 0.5, T).unwrap_err();
    assert!(
        err.to_string().contains("shorter than its corner's reach"),
        "{err}"
    );
}

/// One step of a path drawn in the `xy` plane from the origin heading `+x`:
/// a line of a length, an arc of a radius through a signed angle (positive
/// to the left), tangent to what comes before, or a sharp turn through a
/// signed angle on the spot.
#[derive(Debug, Clone, Copy)]
enum Step {
    Line(f64),
    Arc(f64, f64),
    Turn(f64),
}

/// The wire `steps` draw.
fn planar_path(model: &mut Model, steps: &[Step]) -> Shape {
    let mut at = Point::ORIGIN;
    let mut heading = 0.0_f64;
    let mut from = make_vertex(model, at).shape;
    let mut edges = Vec::new();
    for step in steps {
        match *step {
            Step::Turn(angle) => heading += angle,
            Step::Line(length) => {
                let to = at + Vector::new(heading.cos(), heading.sin(), 0.0) * length;
                let end = make_vertex(model, to).shape;
                let segment = LineCurve::segment(at, to, T).unwrap();
                let range = segment.domain();
                edges.push(
                    make_edge_between(model, segment.into(), range, &from, &end, T)
                        .unwrap()
                        .shape,
                );
                (at, from) = (to, end);
            }
            Step::Arc(radius, angle) => {
                let left = Vector::new(-heading.sin(), heading.cos(), 0.0);
                let centre = at + left * (radius * angle.signum());
                let axis = if angle > 0.0 {
                    Direction::Z
                } else {
                    -Direction::Z
                };
                let frame =
                    Frame::new(centre, axis, Direction::new(at - centre, T).unwrap(), T).unwrap();
                let arc = CircleCurve::new(Circle::new(frame, radius, T).unwrap());
                let to = arc.point_at(angle.abs(), T).unwrap();
                let end = make_vertex(model, to).shape;
                edges.push(
                    make_edge_between(model, arc.into(), (0.0, angle.abs()), &from, &end, T)
                        .unwrap()
                        .shape,
                );
                (at, from, heading) = (to, end, heading + angle);
            }
        }
    }
    make_wire(model, &edges, T).unwrap().shape
}

/// `p` turned about `c` through `angle`, in the plane.
fn turned(p: (f64, f64), c: (f64, f64), angle: f64) -> (f64, f64) {
    let (x, y) = (p.0 - c.0, p.1 - c.1);
    let (s, k) = angle.sin_cos();
    (c.0 + x * k - y * s, c.1 + x * s + y * k)
}

/// Inside a corner turning left through `turn` from the leg `into` onto the
/// leg `out_of`, where the paths of the section's point offset `inside`
/// into the turn cross: how far back along the first leg and on along the
/// second, each along the point's own path. The corner stands at the
/// origin heading `+x`; Newton's method from the crossing two lines would
/// have, on the two path lengths.
fn crossing_inside(into: Step, out_of: Step, turn: f64, inside: f64) -> (f64, f64) {
    let (ahead, left) = ((turn.cos(), turn.sin()), (-turn.sin(), turn.cos()));
    let start = (inside * left.0, inside * left.1);
    let back = |a: f64| match into {
        Step::Arc(radius, angle) => {
            let s = angle.signum();
            turned(
                (0.0, inside),
                (0.0, s * radius),
                -s * a / (radius - s * inside),
            )
        }
        _ => (-a, inside),
    };
    let on = |b: f64| match out_of {
        Step::Arc(radius, angle) => {
            let s = angle.signum();
            let centre = (s * radius * left.0, s * radius * left.1);
            turned(start, centre, s * b / (radius - s * inside))
        }
        _ => (start.0 + b * ahead.0, start.1 + b * ahead.1),
    };
    let miss = |a: f64, b: f64| {
        let (p, q) = (back(a), on(b));
        (p.0 - q.0, p.1 - q.1)
    };
    let guess = inside * (turn / 2.0).tan();
    let (mut a, mut b) = (guess, guess);
    for _ in 0..100 {
        let h = 1e-7;
        let f = miss(a, b);
        let (fa, fb) = (miss(a + h, b), miss(a, b + h));
        let (j11, j21) = ((fa.0 - f.0) / h, (fa.1 - f.1) / h);
        let (j12, j22) = ((fb.0 - f.0) / h, (fb.1 - f.1) / h);
        let det = j11 * j22 - j12 * j21;
        let (da, db) = ((f.0 * j22 - j12 * f.1) / det, (j11 * f.1 - j21 * f.0) / det);
        (a, b) = (a - da, b - db);
        if da.abs().max(db.abs()) < 1e-15 {
            break;
        }
    }
    (a, b)
}

/// The volume of the wedge a disc of radius `r` reaches past the plane
/// square to its spine, cut on a plane leaning `alpha` from that one:
/// each point at distance `ρ` outside the axis runs `ρ tan α` on, over the
/// half disc, `tan α ∫₀ʳ 2ρ √(r² - ρ²) dρ`.
fn mitre_wedge(r: f64, alpha: f64) -> f64 {
    2.0 / 3.0 * r * r * r * alpha.tan()
}

/// Simpson's rule for `f` over `[lo, hi]` in `n` (even) intervals.
fn simpson(f: impl Fn(f64) -> f64, lo: f64, hi: f64, n: u32) -> f64 {
    let h = (hi - lo) / f64::from(n);
    let inner: f64 = (1..n)
        .map(|k| f(lo + h * f64::from(k)) * if k % 2 == 1 { 4.0 } else { 2.0 })
        .sum();
    (f(lo) + f(hi) + inner) * h / 3.0
}

/// The volume of the mitred pipe of radius `r` down `steps`, as the pipe
/// shell sweeps it. Under the rotation-minimizing frame of a path in one
/// plane a point of the disc at offset `ρ` to the left of the spine runs
/// the path offset by `ρ` (a line's length, an arc's `(R ∓ ρ)` times its
/// angle), so a slice of the pipe at height `w` is the strip the offsets
/// `|ρ| ≤ √(r² - w²)` sweep, and the volume is `∫ ℓ(ρ) 2√(r² - ρ²) dρ`
/// over `[-r, r]`, `ℓ(ρ)` the length of the offset path. At a corner
/// turning through `θ` the offset path outside the turn runs straight on
/// to the plane halving it, `|ρ| tan(θ/2)` along either leg; inside, it
/// stops where its two legs' paths cross. Integrated with `ρ = r sin u`,
/// smooth on either side of the spine, by Simpson's rule.
fn generators_volume(steps: &[Step], r: f64) -> f64 {
    let length = |rho: f64| -> f64 {
        let mut total = 0.0;
        for (i, step) in steps.iter().enumerate() {
            total += match *step {
                Step::Line(length) => length,
                Step::Arc(radius, angle) => (radius - angle.signum() * rho) * angle.abs(),
                Step::Turn(turn) => {
                    // Mirrored, where it turns right, to turn left.
                    let m = turn.signum();
                    let mirrored = |s: Step| match s {
                        Step::Arc(radius, angle) => Step::Arc(radius, angle * m),
                        other => other,
                    };
                    let inside = rho * m;
                    if inside <= 0.0 {
                        2.0 * inside.abs() * (turn.abs() / 2.0).tan()
                    } else {
                        let (a, b) = crossing_inside(
                            mirrored(steps[i - 1]),
                            mirrored(steps[i + 1]),
                            turn.abs(),
                            inside,
                        );
                        -(a + b)
                    }
                }
            };
        }
        total
    };
    let slice = |u: f64| {
        let c = u.cos();
        length(r * u.sin()) * 2.0 * r * r * c * c
    };
    simpson(slice, -PI / 2.0, 0.0, 20_000) + simpson(slice, 0.0, PI / 2.0, 20_000)
}

/// A sharp corner turning in an arc's plane, between a line and an arc
/// either way round, or between two arcs, and a sharp corner past an arc,
/// are swept exactly: each leg the cylinder or torus segment it sweeps,
/// the two legs' cuts on the plane halving each corner one curve (no face
/// of that plane is left), the only planes the two ends. The volume is the
/// generators' (see [`generators_volume`]): Pappus's down each leg, the
/// wedges outside each turn, `⅔ r³ tan(θ/2)` each, added, and the legs'
/// overlap inside it taken off. The pipe shell of the disc is that solid.
#[test]
fn a_pipe_round_a_corner_in_an_arc_s_plane_is_its_tori_and_cylinders() {
    use Step::{Arc, Line, Turn};
    let r = 0.5;
    // The wedge's closed form against its integral, and a corner between
    // two lines under the generators' volume: the wedges outside the turn
    // and the overlap inside it cancel, leaving `π r²` times the legs'
    // length.
    for alpha in [0.1, PI / 6.0, PI / 4.0, 1.2] {
        let integral = simpson(
            |u: f64| {
                let (s, c) = u.sin_cos();
                r * s * alpha.tan() * 2.0 * r * c * r * c
            },
            0.0,
            PI / 2.0,
            20_000,
        );
        let wedge = mitre_wedge(r, alpha);
        assert!(
            (integral - wedge).abs() < wedge * 1e-13,
            "{integral} against {wedge}"
        );
        let lines = [Line(4.0), Turn(2.0 * alpha), Line(3.0)];
        let v = generators_volume(&lines, r);
        let want = PI * r * r * 7.0;
        assert!((v - want).abs() < want * 1e-13, "{v} against {want}");
    }
    let paths: [&[Step]; 7] = [
        &[Line(5.0), Turn(PI / 3.0), Arc(3.0, PI / 2.0)],
        &[Line(5.0), Turn(PI / 3.0), Arc(3.0, -PI / 2.0)],
        &[Arc(3.0, PI / 2.0), Turn(PI / 3.0), Line(5.0)],
        &[Arc(3.0, PI / 2.0), Turn(-PI / 3.0), Line(5.0)],
        &[Arc(3.0, PI / 2.0), Turn(-PI / 2.0), Arc(4.0, PI / 2.0)],
        &[
            Line(5.0),
            Arc(3.0, PI / 2.0),
            Line(5.0),
            Turn(PI / 2.0),
            Line(5.0),
        ],
        &[
            Line(4.0),
            Turn(PI / 4.0),
            Arc(2.0, PI / 2.0),
            Turn(-PI / 3.0),
            Line(3.0),
            Arc(1.5, -PI / 2.0),
            Turn(PI / 2.0),
            Arc(2.5, PI / 3.0),
            Line(2.0),
        ],
    ];
    for steps in paths {
        let mut model = Model::new();
        let spine = planar_path(&mut model, steps);
        let pipe = make_pipe(&mut model, &spine, r, T)
            .unwrap_or_else(|e| panic!("{steps:?}: {e}"))
            .shape;
        let found = check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "{steps:?}: {found}");
        assert!(analytic_only(&model, &pipe), "{steps:?}");
        let radii: Vec<f64> = steps
            .iter()
            .filter_map(|s| match s {
                Arc(radius, _) => Some(*radius),
                _ => None,
            })
            .collect();
        let mut planes = 0;
        for face in explore_unique(&model, &pipe, ShapeType::Face).unwrap() {
            let Some(NodeData::Face(d)) = model.node(&face).map(|n| n.data()) else {
                panic!("a face holds no data");
            };
            match model.geometry().surface(d.surface) {
                Some(ogeom::geom::SurfaceGeometry::Torus(t)) => {
                    let t = t.torus();
                    assert!((t.minor_radius() - r).abs() < 1e-12, "{steps:?}");
                    assert!(
                        radii
                            .iter()
                            .any(|radius| (t.major_radius() - radius).abs() < 1e-12),
                        "{steps:?}: a torus of radius {}",
                        t.major_radius()
                    );
                }
                Some(ogeom::geom::SurfaceGeometry::Cylinder(c)) => {
                    assert!((c.cylinder().radius() - r).abs() < 1e-12, "{steps:?}");
                }
                _ => planes += 1,
            }
        }
        assert_eq!(planes, 2, "{steps:?}: a cut on a mitre plane left a face");
        let (v, _, exact) = measure(&model, &pipe);
        assert!(exact, "{steps:?}");
        let want = generators_volume(steps, r);
        assert!(
            (v - want).abs() < want * 1e-9,
            "{steps:?}: {v} against {want}"
        );

        let profile = disc(&mut model, Point::ORIGIN, Direction::X, r);
        let shell = make_pipe_shell(&mut model, &profile, &spine, false, 1e-2, T)
            .unwrap()
            .shape;
        assert!(analytic_only(&model, &shell), "{steps:?}");
        let (w, _, exact_shell) = measure(&model, &shell);
        assert!(exact_shell);
        assert!((v - w).abs() < want * 1e-12, "{steps:?}: {v} against {w}");
    }

    // Past an arc in the `xz` plane, a corner between two lines turning
    // out of it: Pappus's `π r²` times the spine's length, `16 + 3π/2`.
    let mut model = Model::new();
    let spine = line_arc_line(&mut model, true);
    let pipe = make_pipe(&mut model, &spine, r, T).unwrap().shape;
    let found = check(&model, &pipe, T).unwrap();
    assert!(found.is_valid(), "{found}");
    assert!(analytic_only(&model, &pipe));
    let (v, _, exact) = measure(&model, &pipe);
    assert!(exact);
    let want = PI * r * r * (16.0 + 1.5 * PI);
    assert!((v - want).abs() < want * 1e-12, "{v} against {want}");
}

/// The line `up` to `(0, 0, 5)` and on from its end `top`, tangent to
/// it, the cubic Bézier through `(0, 0, 7)` and `(2, 0, 9)` to `(5, 1, 10)`.
fn free_form_after_line(model: &mut Model, up: &Shape, top: &Shape) -> Shape {
    let control = vec![
        Point::new(0.0, 0.0, 5.0),
        Point::new(0.0, 0.0, 7.0),
        Point::new(2.0, 0.0, 9.0),
        Point::new(5.0, 1.0, 10.0),
    ];
    let knots = KnotVector::clamped_uniform(3, 4).unwrap();
    let spline = BSplineCurve::new(knots, control, T).unwrap();
    let end = make_vertex(model, Point::new(5.0, 1.0, 10.0)).shape;
    let bend = make_edge_between(model, spline.into(), (0.0, 1.0), top, &end, T)
        .unwrap()
        .shape;
    make_wire(model, &[up.clone(), bend], T).unwrap().shape
}

/// The volume of `solid` against the tube `π r²` times `length`, to within
/// the skin's tolerance and the mesh's chord over the walls' area, at most
/// the circle's perimeter times the spine's length.
fn assert_tube_volume(model: &Model, solid: &Shape, r: f64, length: f64, skin: f64) {
    const CHORD: f64 = 1e-4;
    let found = check(model, solid, T).unwrap();
    assert!(found.is_valid(), "at {skin}: {found}");
    let v = volume_properties(model, solid, Deflection::with_chord(CHORD).unwrap(), T)
        .unwrap()
        .mass;
    let want = PI * r * r * length;
    let bound = (skin + CHORD) * 2.0 * PI * r * length;
    assert!(
        (v - want).abs() < bound,
        "at {skin}: {v} against {want}, bound {bound}"
    );
}

/// The pipe shell of a disc down a line, round a tangent arc, along a
/// line and round a sharp turn onto a last line closes at every skin
/// tolerance: the arc's skinned run and the straight runs either side of it
/// meet on one section, whichever fit each run settles on. Its volume is
/// `π r²` times the spine's length, `16 + 3π/2`, by Pappus.
#[test]
fn a_disc_piped_past_an_arc_and_a_sharp_turn_closes_at_every_tolerance() {
    let r = 0.5;
    for skin in [1e-2, 3e-3, 1.5e-3, 1e-3, 3e-4, 2e-4, 1e-4] {
        let mut model = Model::new();
        let spine = line_arc_line(&mut model, true);
        let profile = disc(&mut model, Point::ORIGIN, Direction::Z, r);
        let shell = make_pipe_shell(&mut model, &profile, &spine, false, skin, T)
            .unwrap_or_else(|e| panic!("at {skin}: {e}"))
            .shape;
        assert_tube_volume(&model, &shell, r, 16.0 + 1.5 * PI, skin);
    }
}

/// The skinned pipe takes a wire of free-form edges too, as the pipe shell
/// of its disc: down a line and on along a Bézier bend, `π r²` times the
/// spine's length by Pappus, at every skin tolerance: the straight run and
/// the bend's skinned run meet on one section, whichever fit each settles
/// on.
#[test]
fn a_skinned_pipe_runs_along_a_wire_with_a_free_form_edge() {
    let r = 0.5;
    // The bend's length by Simpson's rule over its speed, far finer than
    // the bound below.
    let control = [
        Vector::new(0.0, 0.0, 5.0),
        Vector::new(0.0, 0.0, 7.0),
        Vector::new(2.0, 0.0, 9.0),
        Vector::new(5.0, 1.0, 10.0),
    ];
    let speed = |t: f64| {
        let s = 1.0 - t;
        ((control[1] - control[0]) * (3.0 * s * s)
            + (control[2] - control[1]) * (6.0 * s * t)
            + (control[3] - control[2]) * (3.0 * t * t))
            .magnitude()
    };
    let n = 2000;
    let h = 1.0 / f64::from(n);
    let bend: f64 = (0..=n)
        .map(|i| {
            let w = if i == 0 || i == n {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            w * speed(f64::from(i) * h)
        })
        .sum::<f64>()
        * h
        / 3.0;
    let length = 5.0 + bend;
    for skin in [1e-3, 3e-4, 1e-4, 7e-5, 3e-5] {
        let mut model = Model::new();
        let (a, b) = (Point::ORIGIN, Point::new(0.0, 0.0, 5.0));
        let (va, vb) = (
            make_vertex(&mut model, a).shape,
            make_vertex(&mut model, b).shape,
        );
        let segment = LineCurve::segment(a, b, T).unwrap();
        let range = segment.domain();
        let up = make_edge_between(&mut model, segment.into(), range, &va, &vb, T)
            .unwrap()
            .shape;
        let spine = free_form_after_line(&mut model, &up, &vb);
        let pipe = make_pipe_skinned(&mut model, &spine, r, skin, T)
            .unwrap_or_else(|e| panic!("at {skin}: {e}"))
            .shape;
        assert_tube_volume(&model, &pipe, r, length, skin);
    }
}

/// A `2 × 2` square drawn about `(5, 3, 1)` on the plane facing
/// `(1, 1, 1)`, brought onto a spine up `z` from the origin by contact and
/// correction: the prism `4 · 4`, its centroid half way up the spine, at
/// `(0, 0, 2)`, exactly.
#[test]
fn a_leaning_square_off_the_spine_is_brought_onto_it_and_swept() {
    let mut model = Model::new();
    let at = Point::new(5.0, 3.0, 1.0);
    let normal = Direction::new(Vector::new(1.0, 1.0, 1.0), T).unwrap();
    let frame = Frame::about(at, normal);
    let (x, y) = (frame.x().vector(), frame.y().vector());
    let corners =
        [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)].map(|(a, b)| at + x * a + y * b);
    let ring = make_polygon(&mut model, &corners, true, T).unwrap().shape;
    let square = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[ring],
        T,
    )
    .unwrap()
    .shape;
    let line = LineCurve::segment(Point::ORIGIN, Point::new(0.0, 0.0, 4.0), T).unwrap();
    let range = line.domain();
    let spine = make_edge(&mut model, line.into(), range, T).unwrap().shape;

    let leaning = make_pipe_shell(&mut model, &square, &spine, false, 1e-3, T);
    assert!(leaning.is_err(), "a leaning square off the spine swept");
    let placed = place_on_spine(
        &mut model,
        &square,
        &spine,
        PipePlacement {
            contact: true,
            correction: true,
        },
        T,
    )
    .unwrap()
    .shape;
    let pipe = make_pipe_shell(&mut model, &placed, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    let found = check(&model, &pipe, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let (v, centre, exact) = measure(&model, &pipe);
    assert!(exact);
    assert!((v - 16.0).abs() < 16.0 * 1e-12, "{v}");
    assert!(
        centre.distance(Point::new(0.0, 0.0, 2.0)) < 10.0 * T.confusion(),
        "{centre:?}"
    );
}

/// Contact alone moves a wire unturned: a `2 × 4` rectangle about
/// `(3, 4, 0)` in the `xy` plane lands with its corners at `(±1, ±2, 0)`.
/// Correction alone turns a square leaning a sixth of a turn about `x`
/// through the spine's start, its centroid there already, square to the
/// spine up `z`: the prism `4 · 3`, exactly.
#[test]
fn contact_and_correction_each_move_the_profile_as_named() {
    let mut model = Model::new();
    let line = LineCurve::segment(Point::ORIGIN, Point::new(0.0, 0.0, 3.0), T).unwrap();
    let range = line.domain();
    let spine = make_edge(&mut model, line.into(), range, T).unwrap().shape;

    let rectangle = make_polygon(
        &mut model,
        &[(2.0, 2.0), (4.0, 2.0), (4.0, 6.0), (2.0, 6.0)].map(|(x, y)| Point::new(x, y, 0.0)),
        true,
        T,
    )
    .unwrap()
    .shape;
    let moved = place_on_spine(
        &mut model,
        &rectangle,
        &spine,
        PipePlacement {
            contact: true,
            correction: false,
        },
        T,
    )
    .unwrap()
    .shape;
    assert_eq!(model.kind_of(&moved).unwrap(), ShapeType::Wire);
    for v in explore_unique(&model, &moved, ShapeType::Vertex).unwrap() {
        let Some(NodeData::Vertex(d)) = model.node(&v).map(|n| n.data()) else {
            panic!("a vertex holds no point");
        };
        let p = d.point;
        assert!(
            (p.x.abs() - 1.0).abs() < 1e-12 && (p.y.abs() - 2.0).abs() < 1e-12 && p.z.abs() < 1e-12,
            "{p:?}"
        );
    }

    let tilt = PI / 6.0;
    let normal = Direction::new(Vector::new(0.0, -tilt.sin(), tilt.cos()), T).unwrap();
    let frame = Frame::new(Point::ORIGIN, normal, Direction::X, T).unwrap();
    let (x, y) = (frame.x().vector(), frame.y().vector());
    let corners = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
        .map(|(a, b)| Point::ORIGIN + x * a + y * b);
    let ring = make_polygon(&mut model, &corners, true, T).unwrap().shape;
    let square = make_face(
        &mut model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[ring],
        T,
    )
    .unwrap()
    .shape;
    let turned = place_on_spine(
        &mut model,
        &square,
        &spine,
        PipePlacement {
            contact: false,
            correction: true,
        },
        T,
    )
    .unwrap()
    .shape;
    let pipe = make_pipe_shell(&mut model, &turned, &spine, false, 1e-3, T)
        .unwrap()
        .shape;
    let (v, _, exact) = measure(&model, &pipe);
    assert!(exact);
    assert!((v - 12.0).abs() < 12.0 * 1e-12, "{v}");
}

/// The pipe shell of `profile` down `spine` under the rotation-minimizing
/// frame, mitred, its profile brought on as `contact` asks.
fn piped(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    contact: PipeContact,
) -> ogeom::core::OgeomResult<Shape> {
    let laws = PipeLaws {
        contact,
        ..PipeLaws::default()
    };
    Ok(make_pipe_shell_laws(
        model,
        profile,
        spine,
        &PipeLaw::RotationMinimizing,
        &laws,
        PipeCorners::Mitre,
        1e-3,
        T,
    )?
    .shape)
}

/// A spine straight up `z` from the origin to `(0, 0, length)`.
fn up_z(model: &mut Model, length: f64) -> Shape {
    let line = LineCurve::segment(Point::ORIGIN, Point::new(0.0, 0.0, length), T).unwrap();
    let range = line.domain();
    make_edge(model, line.into(), range, T).unwrap().shape
}

/// Contact moves a profile drawn square to the spine but 5 back from its
/// start along its plane's normal onto the start: down a line, a quarter
/// arc and a line, a disc of radius `r` so moved holds Pappus's
/// `π r² (10 + 3π/2)` as the disc drawn at the start does, and an open half
/// circle so moved sweeps the sheet of the half circle drawn at the start,
/// to the same area. Swept as placed, the profiles off the start are
/// refused, as before.
#[test]
fn contact_moves_a_profile_off_the_start_onto_it() {
    let r = 0.5;
    let mut model = Model::new();
    let spine = line_arc_line(&mut model, false);
    let back = Point::new(0.0, 0.0, -5.0);

    let off = disc(&mut model, back, Direction::Z, r);
    let err = piped(&mut model, &off, &spine, PipeContact::AsPlaced).unwrap_err();
    assert!(err.to_string().contains("does not sit"), "{err}");
    let moved = piped(&mut model, &off, &spine, PipeContact::Contact).unwrap();
    let found = check(&model, &moved, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let at_start = disc(&mut model, Point::ORIGIN, Direction::Z, r);
    let placed = piped(&mut model, &at_start, &spine, PipeContact::AsPlaced).unwrap();
    let want = PI * r * r * (10.0 + 1.5 * PI);
    let (v_moved, c_moved, exact) = measure(&model, &moved);
    let (v_placed, c_placed, _) = measure(&model, &placed);
    assert!(exact);
    assert!((v_moved - want).abs() < want * 1e-12, "{v_moved} vs {want}");
    assert!(
        (v_moved - v_placed).abs() < want * 1e-12,
        "{v_moved} vs {v_placed}"
    );
    assert!(c_moved.distance(c_placed) < 10.0 * T.confusion());

    let half = |model: &mut Model, at: Point| -> Shape {
        let frame = Frame::about(at, Direction::Z);
        let circle = CircleCurve::new(Circle::new(frame, r, T).unwrap());
        let arc = make_edge(model, circle.into(), (0.0, PI), T).unwrap().shape;
        make_wire(model, &[arc], T).unwrap().shape
    };
    let off = half(&mut model, back);
    let err = piped(&mut model, &off, &spine, PipeContact::AsPlaced).unwrap_err();
    assert!(err.to_string().contains("does not sit"), "{err}");
    let moved = piped(&mut model, &off, &spine, PipeContact::Contact).unwrap();
    let at_start = half(&mut model, Point::ORIGIN);
    let placed = piped(&mut model, &at_start, &spine, PipeContact::AsPlaced).unwrap();
    let area = |model: &Model, s: &Shape| {
        surface_properties(model, s, Deflection::default(), T)
            .unwrap()
            .mass
    };
    let (a_moved, a_placed) = (area(&model, &moved), area(&model, &placed));
    assert!(
        (a_moved - a_placed).abs() < a_placed * 1e-12,
        "{a_moved} vs {a_placed}"
    );
}

/// The normals of a solid's planar faces.
fn plane_normals(model: &Model, solid: &Shape) -> Vec<Vector> {
    explore_unique(model, solid, ShapeType::Face)
        .unwrap()
        .iter()
        .filter_map(|f| find_plane(model, f, T).unwrap())
        .map(|p| p.normal().vector())
        .collect()
}

/// Correction turns a disc of radius `r` centred on the start of a
/// straight spine of length `L`, leaning `30°` off square, about the start
/// until it stands square: the cylinder `π r² L`, centred half way up, its
/// two end faces square to the spine. A disc facing back down the spine
/// is square already and is not turned. As placed, or moved into contact
/// without a turn, the leaning disc is refused, as before.
#[test]
fn correction_turns_a_leaning_disc_square_to_the_spine() {
    let (r, length) = (0.75, 4.0);
    let mut model = Model::new();
    let spine = up_z(&mut model, length);
    let tilt = PI / 6.0;
    let leaning = Direction::new(Vector::new(0.0, -tilt.sin(), tilt.cos()), T).unwrap();
    let profile = disc(&mut model, Point::ORIGIN, leaning, r);
    for contact in [PipeContact::AsPlaced, PipeContact::Contact] {
        let err = piped(&mut model, &profile, &spine, contact).unwrap_err();
        assert!(err.to_string().contains("leans"), "{err}");
    }

    let want = PI * r * r * length;
    let back = disc(&mut model, Point::ORIGIN, -Direction::Z, r);
    for profile in [profile, back] {
        let pipe = piped(&mut model, &profile, &spine, PipeContact::Correction).unwrap();
        let found = check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "{found}");
        let (v, centre, exact) = measure(&model, &pipe);
        assert!(exact);
        assert!((v - want).abs() < want * 1e-12, "{v} vs {want}");
        assert!(
            centre.distance(Point::new(0.0, 0.0, length / 2.0)) < 10.0 * T.confusion(),
            "{centre:?}"
        );
        let normals = plane_normals(&model, &pipe);
        assert_eq!(normals.len(), 2);
        for n in normals {
            assert!(n.cross(Vector::Z).magnitude() < 1e-12, "{n:?}");
        }
    }
}

/// Correction is refused by name where it is not one turn: a profile with
/// no plane to stand square, and a disc whose plane holds the spine's
/// start tangent, a quarter turn from square either way round.
#[test]
fn correction_refuses_a_profile_it_cannot_turn_by_name() {
    let mut model = Model::new();
    let spine = up_z(&mut model, 4.0);
    let bent = make_polygon(
        &mut model,
        &[
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.5),
            Point::new(-1.0, 0.0, 0.0),
            Point::new(0.0, -1.0, 0.5),
        ],
        true,
        T,
    )
    .unwrap()
    .shape;
    let err = piped(&mut model, &bent, &spine, PipeContact::Correction).unwrap_err();
    assert!(err.to_string().contains("not planar"), "{err}");

    let along = disc(&mut model, Point::ORIGIN, Direction::X, 0.5);
    let err = piped(&mut model, &along, &spine, PipeContact::Correction).unwrap_err();
    assert!(err.to_string().contains("quarter turn"), "{err}");
}
