//! Pipe shells round smooth closed spines whose section turns against the
//! rotation-minimizing frame: under a leaning binormal, toward an
//! auxiliary spine, or twisted whole turns, measured against closed forms.
//!
//! A section standing square to the spine with its centroid on it holds,
//! by Pappus, its area times the spine's length, however it turns about
//! the spine, wherever it never meets itself.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, edge_vertices, make_edge, make_edge_between, make_face, make_polygon, make_wire,
    volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve3d as _, LineCurve, PlaneSurface};
use ogeom::math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::offset::{
    PipeCorners, PipeEvolution, PipeLaw, make_pipe_shell_evolving, make_pipe_shell_with,
};
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// The skin's tolerance.
const SKIN: f64 = 1e-3;

/// The mesh's chord the volume and the reach are measured at.
const CHORD: f64 = 1e-4;

/// A rectangle `2 a` across (along `-x` turned with `heading`) and `2 b`
/// along `z` about `at`, square to `heading` in the `xy` plane.
fn rectangle(model: &mut Model, at: Point, heading: Vector, a: f64, b: f64) -> Shape {
    let up = Vector::new(0.0, 0.0, 1.0);
    let across = up.cross(heading);
    let across = across / across.magnitude();
    let pts = [(-a, -b), (a, -b), (a, b), (-a, b)].map(|(s, t)| at + across * s + up * t);
    let wire = make_polygon(model, &pts, true, T).unwrap().shape;
    let frame = Frame::new(
        at,
        Direction::new(heading, T).unwrap(),
        Direction::new(across, T).unwrap(),
        T,
    )
    .unwrap();
    make_face(
        model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape
}

/// The circle of `radius` about the origin in the `xy` plane, from
/// `(radius, 0, 0)` heading `+y`.
fn ring(model: &mut Model, radius: f64) -> Shape {
    let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap();
    let circle = Circle::new(frame, radius, T).unwrap();
    let edge = make_edge(model, CircleCurve::new(circle).into(), (0.0, 2.0 * PI), T)
        .unwrap()
        .shape;
    make_wire(model, &[edge], T).unwrap().shape
}

/// The stadium in the `xy` plane: the line `x = r` from `y = -l/2` up to
/// `l/2`, the half circle of radius `r` over the top, the line back down
/// `x = -r` and the half circle under the bottom, each meeting the next
/// tangent to it.
fn stadium(model: &mut Model, r: f64, l: f64) -> Shape {
    let h = l / 2.0;
    let arc = |model: &mut Model, centre: Point, from: f64| -> Shape {
        let frame = Frame::new(centre, Direction::Z, Direction::X, T).unwrap();
        let circle = Circle::new(frame, r, T).unwrap();
        make_edge(model, CircleCurve::new(circle).into(), (from, from + PI), T)
            .unwrap()
            .shape
    };
    let top = arc(model, Point::new(0.0, h, 0.0), 0.0);
    let bottom = arc(model, Point::new(0.0, -h, 0.0), PI);
    let (top_start, top_end) = edge_vertices(model, &top).unwrap().unwrap();
    let (bottom_start, bottom_end) = edge_vertices(model, &bottom).unwrap().unwrap();
    let line = |model: &mut Model, a: Point, b: Point, va, vb| -> Shape {
        let segment = LineCurve::segment(a, b, T).unwrap();
        let range = segment.domain();
        make_edge_between(model, segment.into(), range, va, vb, T)
            .unwrap()
            .shape
    };
    let right = line(
        model,
        Point::new(r, -h, 0.0),
        Point::new(r, h, 0.0),
        &bottom_end,
        &top_start,
    );
    let left = line(
        model,
        Point::new(-r, h, 0.0),
        Point::new(-r, -h, 0.0),
        &top_end,
        &bottom_start,
    );
    make_wire(model, &[right, top, left, bottom], T)
        .unwrap()
        .shape
}

/// The volume, integrated over the mesh at [`CHORD`].
fn measure(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(CHORD).unwrap(), T)
        .unwrap()
        .mass
}

/// The highest `z` the solid's mesh reaches.
fn highest(model: &Model, shape: &Shape) -> f64 {
    let mesh =
        ogeom::mesh::triangulate(model, shape, Deflection::with_chord(CHORD).unwrap(), T).unwrap();
    mesh.positions
        .iter()
        .map(|p| p.z)
        .fold(f64::NEG_INFINITY, f64::max)
}

/// The measured volume of a skinned sweep strays from the swept solid's by
/// at most the skin's tolerance and the mesh's chord over the walls' area:
/// their sum times a bound on the area, the perimeter times the longest
/// path a profile point runs.
fn skin_bound(perimeter: f64, longest_path: f64) -> f64 {
    (SKIN + CHORD) * perimeter * longest_path
}

/// A `2a x 2b` rectangle round a circle of radius `big`, under a binormal
/// leaning off the circle's axis and toward an auxiliary circle off the
/// spine's centre. Each frame turns the rectangle about the spine as it
/// goes, against the rotation-minimizing frame, which keeps its long side
/// level; each comes home to itself, so the ring closes. Pappus gives
/// `2 π big · 4 a b`.
///
/// The leaning binormal `(1, 1, 1)` keeps its part square to the tangent
/// as the section's binormal: at angle `θ` round the circle it leans
/// `ψ = atan2(1, cos θ + sin θ)` off the outward radius, `45°` at the
/// start and `atan2(1, −√2)` at `θ = 5π/4`.
/// The guide, the circle of radius `big` about `(d, 0, h)`, crosses the
/// plane square to the spine at angle `θ` at `ρ = d cos θ + √(big² − d²
/// sin² θ)`, so the frame's normal leans toward `(ρ − big, h)`: from
/// `atan2(h, d)` off level at the start to `atan2(h, −d)` half way round.
/// Either way the rectangle turns past its diagonal's angle `atan(a / b)`
/// from its start, so its corner reaches the half diagonal `√(a² + b²)` in
/// `z`, which no level rectangle does.
#[test]
fn a_rectangle_round_a_circle_turns_with_its_law_and_holds_pappus_s_volume() {
    let (big, a, b) = (10.0_f64, 1.0_f64, 0.5_f64);
    let (d, h) = (2.0_f64, 3.0_f64);
    let leaning = PipeLaw::Binormal(Direction::new(Vector::new(1.0, 1.0, 1.0), T).unwrap());
    assert!(1.0_f64.atan2(-core::f64::consts::SQRT_2) - PI / 4.0 > (a / b).atan());
    // The guide turns it from atan2(h, d) to atan2(h, -d).
    assert!(h.atan2(-d) - h.atan2(d) > (a / b).atan());
    let mut model = Model::new();
    let guide = {
        let frame = Frame::new(Point::new(d, 0.0, h), Direction::Z, Direction::X, T).unwrap();
        let circle = Circle::new(frame, big, T).unwrap();
        make_edge(
            &mut model,
            CircleCurve::new(circle).into(),
            (0.0, 2.0 * PI),
            T,
        )
        .unwrap()
        .shape
    };
    let auxiliary = PipeLaw::Auxiliary { guide: &guide };
    for (name, law) in [("binormal", leaning), ("auxiliary", auxiliary)] {
        let spine = ring(&mut model, big);
        let profile = rectangle(
            &mut model,
            Point::new(big, 0.0, 0.0),
            Vector::new(0.0, 1.0, 0.0),
            a,
            b,
        );
        let pipe = make_pipe_shell_with(
            &mut model,
            &profile,
            &spine,
            &law,
            PipeCorners::Mitre,
            SKIN,
            T,
        )
        .unwrap_or_else(|e| panic!("{name}: {e}"))
        .shape;
        let found = check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "{name}: {found}");
        let v = measure(&model, &pipe);
        let want = 2.0 * PI * big * 4.0 * a * b;
        let reach = a.hypot(b);
        let bound = skin_bound(4.0 * (a + b), 2.0 * PI * (big + reach));
        assert!(
            (v - want).abs() < bound,
            "{name}: {v} against {want}, bound {bound}"
        );
        let top = highest(&model, &pipe);
        assert!(
            (top - reach).abs() < SKIN + CHORD,
            "{name}: reaches z {top}, not {reach}"
        );
    }
}

/// A rectangle round the stadium of two half circles and two lines under
/// the leaning binormal: the frame is the law's own on each edge and
/// continuous across their joins, and the loop closes. Pappus gives the
/// rectangle's area times the stadium's length `2 π r + 2 l`.
#[test]
fn a_rectangle_round_a_stadium_under_a_leaning_binormal_holds_pappus_s_volume() {
    let (r, l, a, b) = (8.0, 12.0, 1.0, 0.5);
    let mut model = Model::new();
    let spine = stadium(&mut model, r, l);
    let profile = rectangle(
        &mut model,
        Point::new(r, -l / 2.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
        a,
        b,
    );
    let law = PipeLaw::Binormal(Direction::new(Vector::new(1.0, 1.0, 1.0), T).unwrap());
    let pipe = make_pipe_shell_with(
        &mut model,
        &profile,
        &spine,
        &law,
        PipeCorners::Mitre,
        SKIN,
        T,
    )
    .unwrap()
    .shape;
    let found = check(&model, &pipe, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let v = measure(&model, &pipe);
    let length = 2.0 * PI * r + 2.0 * l;
    let want = length * 4.0 * a * b;
    let reach = a.hypot(b);
    let bound = skin_bound(4.0 * (a + b), length + 2.0 * PI * reach);
    assert!(
        (v - want).abs() < bound,
        "{v} against {want}, bound {bound}"
    );
}

/// A rectangle twisted one whole turn round a circle under the
/// rotation-minimizing frame comes home to itself. It turns through its
/// diagonal's angle, so its corner reaches the half diagonal in `z`;
/// Pappus gives `2 π big · 4 a b`. A point at `ρ` from the spine runs at
/// most `√((1 + ρ / big)² + (2 π ρ / L)²)` per unit of the spine's length
/// `L`, which bounds the walls' area.
#[test]
fn a_rectangle_twisted_whole_turns_round_a_circle_holds_pappus_s_volume() {
    let (big, a, b) = (10.0, 1.0, 0.5);
    let mut model = Model::new();
    let spine = ring(&mut model, big);
    let profile = rectangle(
        &mut model,
        Point::new(big, 0.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
        a,
        b,
    );
    for turns in [1.0, -2.0] {
        let pipe = make_pipe_shell_evolving(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            PipeEvolution {
                scale: 1.0,
                twist: 2.0 * PI * turns,
            },
            SKIN,
            T,
        )
        .unwrap_or_else(|e| panic!("{turns} turns: {e}"))
        .shape;
        let found = check(&model, &pipe, T).unwrap();
        assert!(found.is_valid(), "{turns} turns: {found}");
        let v = measure(&model, &pipe);
        let want = 2.0 * PI * big * 4.0 * a * b;
        let reach = a.hypot(b);
        let length = 2.0 * PI * big;
        let path = (length + 2.0 * PI * reach).hypot(2.0 * PI * turns.abs() * reach);
        let bound = skin_bound(4.0 * (a + b), path);
        assert!(
            (v - want).abs() < bound,
            "{turns} turns: {v} against {want}, bound {bound}"
        );
        let top = highest(&model, &pipe);
        assert!(
            (top - reach).abs() < SKIN + CHORD,
            "{turns} turns: reaches z {top}, not {reach}"
        );
    }
}

#[test]
fn a_closed_spine_refuses_a_section_that_would_not_come_home_by_name() {
    let mut model = Model::new();
    let spine = ring(&mut model, 10.0);
    let profile = rectangle(
        &mut model,
        Point::new(10.0, 0.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
        1.0,
        0.5,
    );
    let evolve = |scale: f64, twist: f64| PipeEvolution { scale, twist };
    let law = PipeLaw::RotationMinimizing;
    let err = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &law,
        evolve(1.0, PI / 2.0),
        SKIN,
        T,
    )
    .unwrap_err();
    assert!(err.to_string().contains("whole number of turns"), "{err}");
    let err = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &law,
        evolve(0.5, 0.0),
        SKIN,
        T,
    )
    .unwrap_err();
    assert!(err.to_string().contains("would not come home"), "{err}");
}
