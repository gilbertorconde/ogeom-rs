//! Circular pipes along wires of several edges, and profiles brought onto
//! their spine before a pipe shell sweeps them, measured against closed
//! forms.
//!
//! A disc of radius `r` square to a spine, its centre on it, holds by
//! Pappus `π r²` times the spine's length, mitred corners included: each
//! leg is a cylinder cut on planes through its axis's ends.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, make_edge, make_edge_between, make_face, make_polygon, make_vertex, make_wire,
    volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{BSplineCurve, CircleCurve, Curve3d as _, LineCurve, PlaneSurface};
use ogeom::math::{Circle, Direction, Frame, KnotVector, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::offset::{PipePlacement, make_pipe, make_pipe_shell, make_pipe_skinned, place_on_spine};
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

/// What the pipe shell sweeps only as a skin is refused by name, with
/// the function that takes a tolerance for it: a corner turning in an
/// arc's plane, a sharp corner past an arc (the straight run then mitred
/// against the skinned one), a free-form edge. A tube wider than an arc
/// bends is refused as before.
#[test]
fn a_pipe_along_a_wire_refuses_what_it_does_not_sweep_exactly_by_name() {
    let mut model = Model::new();
    let turned = line_arc_line(&mut model, true);
    let err = make_pipe(&mut model, &turned, 0.5, T).unwrap_err();
    assert!(err.to_string().contains("not swept exactly"), "{err}");
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
    let frame = Frame::new(Point::new(0.0, 0.0, 8.0), Direction::Y, -Direction::Z, T).unwrap();
    let arc = CircleCurve::new(Circle::new(frame, 3.0, T).unwrap());
    let vc = make_vertex(&mut model, arc.point_at(PI / 2.0, T).unwrap()).shape;
    let round = make_edge_between(&mut model, arc.into(), (0.0, PI / 2.0), &vb, &vc, T)
        .unwrap()
        .shape;
    let in_plane = make_wire(&mut model, &[up.clone(), round], T)
        .unwrap()
        .shape;
    let err = make_pipe(&mut model, &in_plane, 0.5, T).unwrap_err();
    assert!(err.to_string().contains("not swept exactly"), "{err}");
    let free = free_form_after_line(&mut model, &up, &vb);
    let err = make_pipe(&mut model, &free, 0.5, T).unwrap_err();
    assert!(err.to_string().contains("make_pipe_skinned"), "{err}");
    let tight = line_arc_line(&mut model, false);
    let err = make_pipe(&mut model, &tight, 3.0, T).unwrap_err();
    assert!(err.to_string().contains("swallows its spine"), "{err}");
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
