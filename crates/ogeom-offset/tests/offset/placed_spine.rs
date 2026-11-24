//! A pipe along a spine placed by a location, on its edges and on its
//! wire, is the pipe along the same spine built where it stands: every
//! route of the pipe shell, the evolving and law sweeps, `make_pipe` and
//! `make_pipe_skinned` read the spine where it is placed, and the pipe
//! shell reads an auxiliary guide so too.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{FRAC_PI_2, PI};

use ogeom_core::Tolerances;
use ogeom_geom::{Curve, Curve3d as _, Transformable as _};
use ogeom_math::{Axis, Circle, Direction, Frame, Point, Transform, Vector};
use ogeom_offset::{
    PipeCorners, PipeEvolution, PipeLaw, PipeLaws, SpineFunction, SpineLaw, make_pipe,
    make_pipe_shell_evolving, make_pipe_shell_laws, make_pipe_shell_with, make_pipe_skinned,
};
use ogeom_topo::{Filter, Location, Model, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

/// A piece of a spine where it is drawn.
#[derive(Clone)]
enum Piece {
    Line(Point, Point),
    /// An arc about the z axis through a centre: its radius and angles.
    Arc(Point, f64, (f64, f64)),
    /// A cubic through control points on uniform knots.
    Spline(Vec<Point>),
}

fn curve_of(piece: &Piece) -> (Curve, (f64, f64)) {
    match piece {
        Piece::Line(a, b) => {
            let curve = Curve::Line(ogeom_geom::LineCurve::segment(*a, *b, T).unwrap());
            let domain = curve.domain();
            (curve, domain)
        }
        Piece::Arc(centre, r, range) => {
            let frame = Frame::new(*centre, Direction::Z, Direction::X, T).unwrap();
            let circle = Circle::new(frame, *r, T).unwrap();
            (ogeom_geom::CircleCurve::new(circle).into(), *range)
        }
        Piece::Spline(control) => {
            let n = control.len();
            let mut knots = vec![0.0; 4];
            #[allow(clippy::cast_precision_loss)]
            knots.extend((1..n - 3).map(|k| k as f64 / (n - 3) as f64));
            knots.extend([1.0; 4]);
            let curve = Curve::BSpline(
                ogeom_geom::BSplineCurve::new(
                    ogeom_math::KnotVector::new(knots, 3).unwrap(),
                    control.clone(),
                    T,
                )
                .unwrap(),
            );
            let domain = curve.domain();
            (curve, domain)
        }
    }
}

/// The spine of `pieces` with `motion` baked into its curves: a lone edge,
/// or a wire whose last edge closes on its first vertex where it comes
/// back there.
fn baked(model: &mut Model, pieces: &[Piece], motion: Transform, as_edge: bool) -> Shape {
    let curves: Vec<(Curve, (f64, f64))> = pieces
        .iter()
        .map(|p| {
            let (curve, range) = curve_of(p);
            (curve.transformed(&motion, T).unwrap(), range)
        })
        .collect();
    if as_edge {
        let (curve, range) = curves[0].clone();
        return ogeom_algo::make_edge(model, curve, range, T).unwrap().shape;
    }
    let first = curves[0].0.point_at(curves[0].1.0, T).unwrap();
    let start = ogeom_algo::make_vertex(model, first).shape;
    let mut from = start.clone();
    let mut edges = Vec::new();
    for (i, (curve, range)) in curves.iter().enumerate() {
        let end = curve.point_at(range.1, T).unwrap();
        let to = if i + 1 == curves.len() && end.distance(first) <= T.confusion() {
            start.clone()
        } else {
            ogeom_algo::make_vertex(model, end).shape
        };
        let edge = ogeom_algo::make_edge_between(model, curve.clone(), *range, &from, &to, T);
        edges.push(edge.unwrap().shape);
        from = to;
    }
    ogeom_algo::make_wire(model, &edges, T).unwrap().shape
}

/// The spine of `pieces` built where it is drawn and placed: a lone edge
/// by `outer` over `inner`, a wire's edges by `inner` and the wire by
/// `outer`.
fn placed(
    model: &mut Model,
    pieces: &[Piece],
    (outer, inner): (Transform, Transform),
    as_edge: bool,
) -> Shape {
    let outer = Location::of(model.add_datum(outer));
    let inner = Location::of(model.add_datum(inner));
    let drawn = baked(model, pieces, Transform::IDENTITY, as_edge);
    if as_edge {
        return drawn.moved(&inner).moved(&outer);
    }
    let edges: Vec<Shape> = model
        .ordered_children_of(&drawn)
        .unwrap()
        .iter()
        .map(|e| e.moved(&inner))
        .collect();
    ogeom_algo::make_wire(model, &edges, T)
        .unwrap()
        .shape
        .moved(&outer)
}

/// A unit square square to `tangent` at `centre`: its face, or its wire.
fn square(model: &mut Model, centre: Point, tangent: Vector, face: bool) -> Shape {
    let plane = ogeom_math::Plane::through(centre, Direction::new(tangent, T).unwrap());
    let frame = plane.frame();
    let corners: Vec<Point> = [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]
        .iter()
        .map(|(a, b)| centre + frame.x().vector() * *a + frame.y().vector() * *b)
        .collect();
    let wire = ogeom_algo::make_polygon(model, &corners, true, T)
        .unwrap()
        .shape;
    if !face {
        return wire;
    }
    let surface: ogeom_geom::SurfaceGeometry =
        ogeom_geom::PlaneSurface::over(plane, (-2.0, 2.0), (-2.0, 2.0))
            .unwrap()
            .into();
    ogeom_algo::make_face(model, surface, std::slice::from_ref(&wire), T)
        .unwrap()
        .shape
}

/// A solid's volume or a shell's area, the shape's bound, and whether it
/// checks clean.
fn measured(model: &Model, shape: &Shape) -> (f64, Point, Point, bool) {
    let solid = !explore(model, shape, Filter::OfType(ShapeType::Solid))
        .unwrap()
        .is_empty();
    let deflection = ogeom_mesh::Deflection::default();
    let mass = if solid {
        ogeom_algo::volume_properties(model, shape, deflection, T)
    } else {
        ogeom_algo::surface_properties(model, shape, deflection, T)
    };
    let bound = ogeom_algo::shape_bounds(model, shape, T).unwrap();
    let valid = ogeom_algo::check(model, shape, T).unwrap().is_valid();
    (
        mass.unwrap().mass,
        bound.low().unwrap(),
        bound.high().unwrap(),
        valid,
    )
}

/// What is swept along the spine, and how.
#[derive(Clone, Copy, Debug)]
enum Route {
    /// The pipe shell of the unit square's face or wire, under a law
    /// (rotation-minimizing, Frenet, auxiliary with the guide placed as the
    /// spine is, binormal) and a way of turning corners.
    Shell {
        law: u8,
        corners: PipeCorners,
        face: bool,
    },
    /// The square scaled to a half and twisted by a radian.
    Evolving,
    /// The square twisted by a law given as a function.
    TwistFunction,
    Pipe,
    PipeSkinned,
}

/// Sweep along `spine` as `route` asks, the profile square to the spine's
/// start, where `motion` places the spine's pieces; an auxiliary guide is
/// the spine lifted three along the drawing's z, placed as `place` asks.
fn swept(
    model: &mut Model,
    route: Route,
    pieces: &[Piece],
    spine: &Shape,
    motion: Transform,
    place: Option<(Transform, Transform)>,
) -> Result<Shape, String> {
    let (at, tangent) = {
        let (curve, range) = curve_of(&pieces[0]);
        let curve = curve.transformed(&motion, T).unwrap();
        (
            curve.point_at(range.0, T).unwrap(),
            curve.d1_at(range.0, T).unwrap(),
        )
    };
    let up = Direction::new(motion.apply_vector(Vector::new(0.0, 0.0, 1.0)), T).unwrap();
    let lifted: Vec<Piece> = pieces
        .iter()
        .map(|p| {
            let lift = |q: &Point| *q + Vector::new(0.0, 0.0, 3.0);
            match p {
                Piece::Line(a, b) => Piece::Line(lift(a), lift(b)),
                Piece::Arc(c, r, range) => Piece::Arc(lift(c), *r, *range),
                Piece::Spline(points) => Piece::Spline(points.iter().map(lift).collect()),
            }
        })
        .collect();
    let guide = match place {
        Some(place) => placed(model, &lifted, place, false),
        None => baked(model, &lifted, motion, false),
    };
    let built = match route {
        Route::Shell { law, corners, face } => {
            let profile = square(model, at, tangent, face);
            let law = match law {
                0 => PipeLaw::RotationMinimizing,
                1 => PipeLaw::Frenet,
                2 => PipeLaw::Auxiliary { guide: &guide },
                _ => PipeLaw::Binormal(up),
            };
            make_pipe_shell_with(model, &profile, spine, &law, corners, 1e-3, T)
        }
        Route::Evolving => {
            let profile = square(model, at, tangent, true);
            let evolution = PipeEvolution {
                scale: 0.5,
                twist: 1.0,
            };
            let law = PipeLaw::RotationMinimizing;
            make_pipe_shell_evolving(model, &profile, spine, &law, evolution, 1e-3, T)
        }
        Route::TwistFunction => {
            let profile = square(model, at, tangent, true);
            let laws = PipeLaws {
                twist: SpineLaw::Function(SpineFunction::new(|s| 2.0 * PI * s * s)),
                ..PipeLaws::default()
            };
            let law = PipeLaw::RotationMinimizing;
            make_pipe_shell_laws(
                model,
                &profile,
                spine,
                &law,
                &laws,
                PipeCorners::Mitre,
                1e-3,
                T,
            )
        }
        Route::Pipe => make_pipe(model, spine, 0.5, T),
        Route::PipeSkinned => make_pipe_skinned(model, spine, 0.5, 1e-3, T),
    };
    built.map(|b| b.shape).map_err(|e| e.to_string())
}

/// The placement: a turn about a skew axis over a shift.
fn placement() -> (Transform, Transform) {
    let outer = Transform::rotation(
        Axis {
            location: Point::new(1.0, 2.0, 3.0),
            direction: Direction::new(Vector::new(1.0, 2.0, -0.5), T).unwrap(),
        },
        1.1,
    );
    (outer, Transform::translation(Vector::new(5.0, -3.0, 2.0)))
}

/// Sweep `route` along `pieces` with the spine or the guide placed as
/// `spine_placed` and `guide_placed` ask, and with both baked in place, and
/// hold the two to one volume or area, one bound and one verdict of the
/// checker.
fn same_as_baked(
    name: &str,
    route: Route,
    pieces: &[Piece],
    as_edge: bool,
    (spine_placed, guide_placed): (bool, bool),
) {
    let (outer, inner) = placement();
    let motion = outer * inner;
    let case = format!("{name} {route:?}");
    let mut model = Model::new();
    let spine = if spine_placed {
        placed(&mut model, pieces, (outer, inner), as_edge)
    } else {
        baked(&mut model, pieces, motion, as_edge)
    };
    let place = guide_placed.then_some((outer, inner));
    let got = swept(&mut model, route, pieces, &spine, motion, place)
        .unwrap_or_else(|e| panic!("{case}, placed: {e}"));
    let mut in_place = Model::new();
    let spine = baked(&mut in_place, pieces, motion, as_edge);
    let want = swept(&mut in_place, route, pieces, &spine, motion, None)
        .unwrap_or_else(|e| panic!("{case}, baked: {e}"));
    let (m, low, high, valid) = measured(&model, &got);
    let (m0, low0, high0, valid0) = measured(&in_place, &want);
    assert!(valid0, "{case}: the sweep in place checks clean");
    assert!(valid, "{case}: the placed sweep checks clean");
    assert!(
        (m - m0).abs() <= 1e-9 * m0.abs().max(1.0),
        "{case}: {m} against {m0}"
    );
    let off = low.distance(low0).max(high.distance(high0));
    assert!(off <= 1e-9, "{case}: the bounds stand {off} apart");
}

fn arc() -> Vec<Piece> {
    vec![Piece::Arc(Point::ORIGIN, 5.0, (0.0, FRAC_PI_2))]
}

fn ell() -> Vec<Piece> {
    vec![
        Piece::Line(Point::ORIGIN, Point::new(10.0, 0.0, 0.0)),
        Piece::Line(Point::new(10.0, 0.0, 0.0), Point::new(10.0, 10.0, 0.0)),
    ]
}

fn circle() -> Vec<Piece> {
    vec![Piece::Arc(Point::ORIGIN, 5.0, (0.0, 2.0 * PI))]
}

fn square_loop() -> Vec<Piece> {
    let at = [
        (0.0, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (0.0, 10.0),
        (0.0, 0.0),
    ];
    at.windows(2)
        .map(|w| {
            Piece::Line(
                Point::new(w[0].0, w[0].1, 0.0),
                Point::new(w[1].0, w[1].1, 0.0),
            )
        })
        .collect()
}

fn stadium() -> Vec<Piece> {
    vec![
        Piece::Line(Point::new(0.0, -5.0, 0.0), Point::new(10.0, -5.0, 0.0)),
        Piece::Arc(Point::new(10.0, 0.0, 0.0), 5.0, (-FRAC_PI_2, FRAC_PI_2)),
        Piece::Line(Point::new(10.0, 5.0, 0.0), Point::new(0.0, 5.0, 0.0)),
        Piece::Arc(Point::ORIGIN, 5.0, (FRAC_PI_2, 3.0 * FRAC_PI_2)),
    ]
}

fn wave() -> Vec<Piece> {
    vec![Piece::Spline(vec![
        Point::ORIGIN,
        Point::new(3.0, 2.0, 0.0),
        Point::new(6.0, -2.0, 0.0),
        Point::new(9.0, 0.0, 0.0),
    ])]
}

/// The pipe shell of a face or a closed wire, by each route it takes:
/// skinned along a smooth spine, exact round straight legs' corners, run
/// by run round a law's corners, closed round a smooth or a cornered loop.
#[test]
fn a_pipe_shell_along_a_placed_spine_is_the_shell_along_it_in_place() {
    let shell = |law: u8, corners: PipeCorners, face: bool| Route::Shell { law, corners, face };
    let (mitre, extended, round) = (
        PipeCorners::Mitre,
        PipeCorners::Extended,
        PipeCorners::Round,
    );
    let cases: [(&str, Vec<Piece>, bool, Route); 9] = [
        ("arc edge", arc(), true, shell(0, mitre, true)),
        ("arc wire", arc(), false, shell(1, round, false)),
        ("wave edge", wave(), true, shell(3, mitre, true)),
        ("ell", ell(), false, shell(0, mitre, true)),
        ("ell", ell(), false, shell(0, extended, true)),
        ("ell", ell(), false, shell(3, round, true)),
        ("circle", circle(), false, shell(0, mitre, true)),
        ("stadium", stadium(), false, shell(1, mitre, false)),
        ("square", square_loop(), false, shell(0, extended, true)),
    ];
    for (name, pieces, as_edge, route) in &cases {
        same_as_baked(name, *route, pieces, *as_edge, (true, false));
    }
}

/// An auxiliary guide placed by a location steers the section where it
/// stands.
#[test]
fn a_placed_auxiliary_guide_steers_the_section_where_it_stands() {
    let route = Route::Shell {
        law: 2,
        corners: PipeCorners::Mitre,
        face: true,
    };
    for (name, pieces, as_edge) in [("arc", arc(), true), ("circle", circle(), false)] {
        same_as_baked(name, route, &pieces, as_edge, (false, true));
    }
}

/// The evolving sweep, a twist law given as a function, and the circular
/// pipe exact or skinned, along placed spines.
#[test]
fn evolving_sweeps_and_pipes_along_a_placed_spine_are_those_in_place() {
    let cases: [(&str, Vec<Piece>, bool, Route); 8] = [
        ("arc", arc(), true, Route::Evolving),
        ("square", square_loop(), false, Route::Evolving),
        ("arc", arc(), false, Route::TwistFunction),
        ("arc", arc(), true, Route::Pipe),
        ("ell", ell(), false, Route::Pipe),
        ("circle", circle(), true, Route::Pipe),
        ("wave", wave(), true, Route::PipeSkinned),
        ("stadium", stadium(), false, Route::PipeSkinned),
    ];
    for (name, pieces, as_edge, route) in &cases {
        same_as_baked(name, *route, pieces, *as_edge, (true, false));
    }
}
