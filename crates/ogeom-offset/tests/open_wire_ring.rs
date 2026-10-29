//! Open wires swept round closed spines of lines and arcs: a ring of the
//! sheet an open spine sweeps, closed on itself where the spine closes,
//! against closed forms.
//!
//! A segment round a circle sweeps a plane annulus, a cylinder band or a
//! cone band; any profile round a circle sweeps Pappus' area, its length
//! times the path of its centroid. Round a closed spine's straight legs and
//! corners each point's path is as down an open spine's: `α tan(θ/2)` on
//! along each leg for a point `α` outside a mitred corner turning `θ`, or
//! the arc `α θ` round a rounded one.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{FRAC_PI_2, PI};

use ogeom_algo::{check, curve_length, make_polygon, surface_properties};
use ogeom_core::{FastMap, Tolerances};
use ogeom_math::{Circle, Direction, Frame, Point, Transform, Vector};
use ogeom_mesh::Deflection;
use ogeom_offset::{PipeCorners, PipeLaw, make_pipe_shell_with};
use ogeom_topo::{EdgeRepr, Filter, Location, Model, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

/// The surfaces are planes, cylinders, cones and the surfaces a line turns
/// about an axis, measured exactly; the bound covers rounding.
const EXACT: f64 = 1e-9;

fn close(got: f64, want: f64) -> bool {
    (got - want).abs() <= want.abs() * EXACT
}

fn sweep(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    corners: PipeCorners,
) -> Result<Shape, String> {
    make_pipe_shell_with(model, profile, spine, law, corners, 1e-3, T)
        .map(|b| b.shape)
        .map_err(|e| e.to_string())
}

/// A sheet's area and the length of its boundary after checking it is one
/// shell with no finding but its being open, every edge of one face or
/// two, bounded by two loops.
fn measured(model: &Model, sheet: &Shape) -> (f64, f64) {
    assert_eq!(model.kind_of(sheet).unwrap(), ShapeType::Shell);
    let diagnosis = check(model, sheet, T).unwrap();
    assert!(
        diagnosis
            .problems
            .iter()
            .all(|p| p.kind == ShapeType::Shell),
        "{diagnosis}"
    );
    let area = surface_properties(model, sheet, Deflection::with_chord(1e-4).unwrap(), T)
        .unwrap()
        .mass;
    let mut uses: FastMap<ogeom_topo::TShapeId, (Shape, usize)> = FastMap::default();
    for face in explore(model, sheet, Filter::OfType(ShapeType::Face)).unwrap() {
        for edge in explore(model, &face, Filter::OfType(ShapeType::Edge)).unwrap() {
            uses.entry(edge.node()).or_insert((edge, 0)).1 += 1;
        }
    }
    assert!(
        uses.values().all(|(_, n)| *n <= 2),
        "an edge of three faces"
    );
    let mut boundary = 0.0;
    for (edge, _) in uses.values().filter(|(_, n)| *n == 1) {
        let data = model.node(edge).unwrap().data().as_edge().unwrap();
        if data.degenerate {
            continue;
        }
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            panic!("an edge with no curve");
        };
        let curve = model.geometry().curve(*curve).unwrap();
        boundary += curve_length(curve, *range, T).unwrap();
    }
    (area, boundary)
}

/// An open polyline through `points`.
fn polyline(model: &mut Model, points: &[Point]) -> Shape {
    make_polygon(model, points, false, T).unwrap().shape
}

/// The circle of radius `r` about the z axis, from `(r, 0, 0)` heading
/// `+y`, as a wire of one edge.
fn circle(model: &mut Model, r: f64) -> Shape {
    let circle = Circle::new(Frame::WORLD, r, T).unwrap();
    let curve: ogeom_geom::Curve = ogeom_geom::CircleCurve::new(circle).into();
    let edge = ogeom_algo::make_edge(model, curve, (0.0, 2.0 * PI), T)
        .unwrap()
        .shape;
    ogeom_algo::make_wire(model, &[edge], T).unwrap().shape
}

/// The arc of radius `r` about `centre` in the plane `z = 0` between the
/// angles `range`, on the vertices `from` and `to`.
fn arc(
    model: &mut Model,
    centre: Point,
    r: f64,
    range: (f64, f64),
    (from, to): (&Shape, &Shape),
) -> Shape {
    let frame = Frame::new(centre, Direction::Z, Direction::X, T).unwrap();
    let curve: ogeom_geom::Curve =
        ogeom_geom::CircleCurve::new(Circle::new(frame, r, T).unwrap()).into();
    ogeom_algo::make_edge_between(model, curve, range, from, to, T)
        .unwrap()
        .shape
}

/// The stadium of two legs of 10 and two half circles of radius 5 about
/// `(0, 0, 0)` and `(10, 0, 0)`, from `(0, -5, 0)` heading `+x`.
fn stadium(model: &mut Model) -> Shape {
    let at = [(0.0, -5.0), (10.0, -5.0), (10.0, 5.0), (0.0, 5.0)];
    let v: Vec<Shape> = at
        .iter()
        .map(|&(x, y)| ogeom_algo::make_vertex(model, Point::new(x, y, 0.0)).shape)
        .collect();
    let line = |model: &mut Model, a: usize, b: usize| {
        let (p, q) = (
            Point::new(at[a].0, at[a].1, 0.0),
            Point::new(at[b].0, at[b].1, 0.0),
        );
        let curve = ogeom_geom::Curve::Line(ogeom_geom::LineCurve::segment(p, q, T).unwrap());
        ogeom_algo::make_edge_between(model, curve, (0.0, p.distance(q)), &v[a], &v[b], T)
            .unwrap()
            .shape
    };
    let edges = [
        line(model, 0, 1),
        arc(
            model,
            Point::new(10.0, 0.0, 0.0),
            5.0,
            (-FRAC_PI_2, FRAC_PI_2),
            (&v[1], &v[2]),
        ),
        line(model, 2, 3),
        arc(
            model,
            Point::ORIGIN,
            5.0,
            (FRAC_PI_2, 3.0 * FRAC_PI_2),
            (&v[3], &v[0]),
        ),
    ];
    ogeom_algo::make_wire(model, &edges, T).unwrap().shape
}

/// The square of side 10 in the plane `z = 0`, from the origin heading
/// `+x` and turning left at each corner.
fn square(model: &mut Model) -> Shape {
    let at = [
        Point::ORIGIN,
        Point::new(10.0, 0.0, 0.0),
        Point::new(10.0, 10.0, 0.0),
        Point::new(0.0, 10.0, 0.0),
    ];
    make_polygon(model, &at, true, T).unwrap().shape
}

/// A segment round a circle in the plane square to it at the start sweeps
/// the surface of revolution of the segment about the circle's axis: a
/// plane annulus lying in the circle's plane, a cylinder band standing
/// along its axis, a cone band leaning; an L, Pappus' area.
#[test]
fn a_profile_round_a_circle_sweeps_its_revolution() {
    let r = 5.0;
    let at = |x: f64, z: f64| Point::new(x, 0.0, z);
    let annulus = PI * ((r + 2.0) * (r + 2.0) - (r + 1.0) * (r + 1.0));
    let cylinder = 2.0 * PI * (r + 1.0) * 2.0;
    let cone = PI * ((r + 1.0) + (r + 2.0)) * 2.0_f64.sqrt();
    let cases: [(&str, Vec<Point>, f64, f64); 4] = [
        (
            "annulus",
            vec![at(r + 1.0, 0.0), at(r + 2.0, 0.0)],
            annulus,
            2.0 * PI * (2.0 * r + 3.0),
        ),
        (
            "cylinder",
            vec![at(r + 1.0, 0.0), at(r + 1.0, 2.0)],
            cylinder,
            2.0 * 2.0 * PI * (r + 1.0),
        ),
        (
            "cone",
            vec![at(r + 1.0, 0.0), at(r + 2.0, 1.0)],
            cone,
            2.0 * PI * (2.0 * r + 3.0),
        ),
        (
            "L, Pappus",
            vec![at(r - 1.0, 0.0), at(r + 1.0, 0.0), at(r + 1.0, 1.0)],
            2.0 * 2.0 * PI * r + 1.0 * 2.0 * PI * (r + 1.0),
            2.0 * PI * (r - 1.0) + 2.0 * PI * (r + 1.0),
        ),
    ];
    for (name, points, want, around) in cases {
        for law in [PipeLaw::RotationMinimizing, PipeLaw::Frenet] {
            let mut model = Model::new();
            let spine = circle(&mut model, r);
            let profile = polyline(&mut model, &points);
            let sheet = sweep(&mut model, &profile, &spine, &law, PipeCorners::Mitre)
                .unwrap_or_else(|e| panic!("{name} {law:?}: {e}"));
            let (area, boundary) = measured(&model, &sheet);
            assert!(close(area, want), "{name} {law:?}: {area} against {want}");
            assert!(
                close(boundary, around),
                "{name} {law:?}: {boundary} against {around}"
            );
        }
    }
}

/// Round a stadium the legs sweep prisms and the half circles half turns,
/// all joined tangent: a segment outside it sweeps `20 + 13π`; turned to
/// stand along the axis, its band the stadium's length times its height.
#[test]
fn a_segment_round_a_stadium_sweeps_its_prisms_and_turns() {
    let at = |y: f64, z: f64| Point::new(0.0, y, z);
    let cases: [(&str, [Point; 2], f64, f64); 2] = [
        (
            "flat",
            [at(-6.0, 0.0), at(-7.0, 0.0)],
            20.0 + 13.0 * PI,
            (20.0 + 12.0 * PI) + (20.0 + 14.0 * PI),
        ),
        (
            "standing",
            [at(-6.0, 0.0), at(-6.0, 2.0)],
            2.0 * (20.0 + 12.0 * PI),
            2.0 * (20.0 + 12.0 * PI),
        ),
    ];
    for (name, points, want, around) in cases {
        let mut model = Model::new();
        let spine = stadium(&mut model);
        let profile = polyline(&mut model, &points);
        let sheet = sweep(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            PipeCorners::Mitre,
        )
        .unwrap_or_else(|e| panic!("{name}: {e}"));
        let (area, boundary) = measured(&model, &sheet);
        assert!(close(area, want), "{name}: {area} against {want}");
        assert!(
            close(boundary, around),
            "{name}: {boundary} against {around}"
        );
    }
}

/// Round a closed square, its spine starting at a corner, each corner is
/// turned as any other: a segment `1` to `2` outside the legs sweeps
/// `40 + 4 ∫ 2α dα = 52` mitred and `40 + 4 ∫ α π/2 dα = 40 + 3π` rounded;
/// one standing `2` high at `1` outside, `2 (40 + 8)` and `2 (40 + 2π)`;
/// one crossing the spine from `1` inside to `1` outside, its inner half
/// cut on the planes halving the corners, `80` mitred and `80 + 4 (π/4 -
/// 1)` rounded. The binormal law square to the square's plane sweeps the
/// same sheet.
#[test]
fn a_segment_round_a_closed_square_turns_every_corner() {
    let at = |y: f64, z: f64| Point::new(0.0, y, z);
    let cases: [(&str, [Point; 2], PipeCorners, f64, f64); 6] = [
        (
            "crossing mitred",
            [at(1.0, 0.0), at(-1.0, 0.0)],
            PipeCorners::Mitre,
            80.0,
            48.0 + 32.0,
        ),
        (
            "crossing rounded",
            [at(1.0, 0.0), at(-1.0, 0.0)],
            PipeCorners::Round,
            76.0 + PI,
            (40.0 + 2.0 * PI) + 32.0,
        ),
        (
            "flat mitred",
            [at(-1.0, 0.0), at(-2.0, 0.0)],
            PipeCorners::Mitre,
            52.0,
            (40.0 + 8.0) + (40.0 + 16.0),
        ),
        (
            "flat rounded",
            [at(-1.0, 0.0), at(-2.0, 0.0)],
            PipeCorners::Round,
            40.0 + 3.0 * PI,
            (40.0 + 2.0 * PI) + (40.0 + 4.0 * PI),
        ),
        (
            "standing mitred",
            [at(-1.0, 0.0), at(-1.0, 2.0)],
            PipeCorners::Mitre,
            2.0 * 48.0,
            2.0 * 48.0,
        ),
        (
            "standing rounded",
            [at(-1.0, 0.0), at(-1.0, 2.0)],
            PipeCorners::Round,
            2.0 * (40.0 + 2.0 * PI),
            2.0 * (40.0 + 2.0 * PI),
        ),
    ];
    for (name, points, corners, want, around) in cases {
        for law in [PipeLaw::RotationMinimizing, PipeLaw::Binormal(Direction::Z)] {
            let mut model = Model::new();
            let spine = square(&mut model);
            let profile = polyline(&mut model, &points);
            let sheet = sweep(&mut model, &profile, &spine, &law, corners)
                .unwrap_or_else(|e| panic!("{name} {law:?}: {e}"));
            let (area, boundary) = measured(&model, &sheet);
            assert!(close(area, want), "{name} {law:?}: {area} against {want}");
            assert!(
                close(boundary, around),
                "{name} {law:?}: {boundary} against {around}"
            );
        }
    }
}

/// Round a D, a half circle of radius 5 closed by its diameter, its spine
/// starting at a corner on the arc, a segment `1` to `2` outside sweeps
/// the half annulus, the leg's band and at each corner `∫ 2α dα = 3`
/// mitred, the arc run on straight past it: `16 + 6.5π`.
#[test]
fn a_segment_round_a_d_mitres_its_corners_against_the_arc() {
    let mut model = Model::new();
    let (a, b) = (Point::new(0.0, -5.0, 0.0), Point::new(0.0, 5.0, 0.0));
    let (va, vb) = (
        ogeom_algo::make_vertex(&mut model, a).shape,
        ogeom_algo::make_vertex(&mut model, b).shape,
    );
    let half = arc(
        &mut model,
        Point::ORIGIN,
        5.0,
        (-FRAC_PI_2, FRAC_PI_2),
        (&va, &vb),
    );
    let curve = ogeom_geom::Curve::Line(ogeom_geom::LineCurve::segment(b, a, T).unwrap());
    let diameter = ogeom_algo::make_edge_between(&mut model, curve, (0.0, 10.0), &vb, &va, T)
        .unwrap()
        .shape;
    let spine = ogeom_algo::make_wire(&mut model, &[half, diameter], T)
        .unwrap()
        .shape;
    let profile = polyline(
        &mut model,
        &[Point::new(0.0, -6.0, 0.0), Point::new(0.0, -7.0, 0.0)],
    );
    let sheet = sweep(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        PipeCorners::Mitre,
    )
    .unwrap();
    let (area, boundary) = measured(&model, &sheet);
    let (want, around) = (16.0 + 6.5 * PI, (14.0 + 6.0 * PI) + (18.0 + 7.0 * PI));
    assert!(close(area, want), "{area} against {want}");
    assert!(close(boundary, around), "{boundary} against {around}");
}

/// A closed spine placed by a location sweeps the ring the spine where it
/// stands sweeps.
#[test]
fn a_ring_round_a_placed_closed_spine_is_the_ring_where_it_stands() {
    let motion = Transform::rotation(
        ogeom_math::Axis {
            location: Point::new(1.0, 2.0, 3.0),
            direction: Direction::new(Vector::new(1.0, 2.0, -0.5), T).unwrap(),
        },
        1.1,
    ) * Transform::translation(Vector::new(5.0, -3.0, 2.0));
    let points = [Point::new(0.0, -1.0, 0.0), Point::new(0.0, -2.0, 1.0)];
    let moved: Vec<Point> = points.iter().map(|p| motion.apply(*p)).collect();
    let mut model = Model::new();
    let place = Location::of(model.add_datum(motion));
    let spine = square(&mut model).moved(&place);
    let profile = polyline(&mut model, &moved);
    let sheet = sweep(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        PipeCorners::Round,
    )
    .unwrap();
    let mut drawn = Model::new();
    let spine = square(&mut drawn);
    let profile = polyline(&mut drawn, &points);
    let in_place = sweep(
        &mut drawn,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        PipeCorners::Round,
    )
    .unwrap();
    let (area, boundary) = measured(&model, &sheet);
    let (want, around) = measured(&drawn, &in_place);
    assert!(close(area, want), "{area} against {want}");
    assert!(close(boundary, around), "{boundary} against {around}");
}

/// Round a closed spine of lines that does not lie in one plane the
/// rotation-minimizing frame comes home turned (by about 19 degrees round
/// this pentagon), so the sheet would not
/// close onto itself; and round a closed spine with a free-form edge no
/// ring is built. Each is refused by name.
#[test]
fn a_ring_that_would_not_close_is_refused_by_name() {
    let mut model = Model::new();
    let skew = make_polygon(
        &mut model,
        &[
            Point::ORIGIN,
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0, 10.0, 5.0),
            Point::new(0.0, 10.0, 0.0),
            Point::new(-3.0, 4.0, -4.0),
        ],
        true,
        T,
    )
    .unwrap()
    .shape;
    let profile = polyline(
        &mut model,
        &[Point::new(0.0, -1.0, 0.0), Point::new(0.0, -2.0, 0.0)],
    );
    let refusal = sweep(
        &mut model,
        &profile,
        &skew,
        &PipeLaw::RotationMinimizing,
        PipeCorners::Mitre,
    )
    .unwrap_err();
    assert!(refusal.contains("comes home turned"), "{refusal}");

    let mut model = Model::new();
    let control = vec![
        Point::new(5.0, 0.0, 0.0),
        Point::new(5.0, 8.0, 0.0),
        Point::new(-5.0, 8.0, 0.0),
        Point::new(-5.0, -8.0, 0.0),
        Point::new(5.0, -8.0, 0.0),
        Point::new(5.0, 0.0, 0.0),
    ];
    let knots = vec![0.0, 0.0, 0.0, 0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0, 1.0, 1.0, 1.0];
    let curve = ogeom_geom::Curve::BSpline(
        ogeom_geom::BSplineCurve::new(ogeom_math::KnotVector::new(knots, 3).unwrap(), control, T)
            .unwrap(),
    );
    let edge = ogeom_algo::make_edge(&mut model, curve, (0.0, 1.0), T)
        .unwrap()
        .shape;
    let loop_ = ogeom_algo::make_wire(&mut model, &[edge], T).unwrap().shape;
    let profile = polyline(
        &mut model,
        &[Point::new(6.0, 0.0, 0.0), Point::new(7.0, 0.0, 0.0)],
    );
    let refusal = sweep(
        &mut model,
        &profile,
        &loop_,
        &PipeLaw::RotationMinimizing,
        PipeCorners::Mitre,
    )
    .unwrap_err();
    assert!(
        refusal.contains("closed spine of lines and arcs"),
        "{refusal}"
    );
}
