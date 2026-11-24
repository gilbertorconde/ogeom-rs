//! An operation on inputs placed by a location, on their edges and on the
//! wire or face that holds them, builds what it builds from the same
//! inputs drawn where they stand: the helical sweep, the loft to a point,
//! the multisection pipe, the evolved sweep, bringing a profile onto a
//! spine, the four-sided filling, the wire offset, the tapered prism and a
//! face given its pcurves all read edges where they are placed. A prism of
//! a mirrored profile measures the volume it encloses.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::FRAC_PI_2;

use ogeom_core::Tolerances;
use ogeom_geom::{Curve, Curve3d as _, Transformable as _};
use ogeom_math::{Axis, Circle, Direction, Frame, Plane, Point, Transform, Vector};
use ogeom_offset::{Join, PipePlacement};
use ogeom_topo::{Filter, Location, Model, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

/// A piece of a wire where it is drawn.
#[derive(Clone)]
enum Piece {
    Line(Point, Point),
    /// An arc about the z axis through a centre: its radius and angles.
    Arc(Point, f64, (f64, f64)),
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
    }
}

/// The straight pieces round `points`, closed back to the first.
fn polygon(points: &[Point]) -> Vec<Piece> {
    (0..points.len())
        .map(|i| Piece::Line(points[i], points[(i + 1) % points.len()]))
        .collect()
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

fn motion() -> Transform {
    let (outer, inner) = placement();
    outer * inner
}

/// How an input is given: drawn where it stands, placed whole by the
/// motion, or its edges placed by the shift inside a wire or face placed
/// by the turn.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Given {
    Baked,
    Placed,
    Nested,
}

fn location(model: &mut Model, transform: Transform) -> Location {
    Location::of(model.add_datum(transform))
}

/// A wire of `pieces` with `motion` baked into its curves, or a lone edge.
fn drawn(model: &mut Model, pieces: &[Piece], motion: Transform, as_edge: bool) -> Shape {
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

/// The edges of a wire drawn at the origin, each placed by `inner`.
fn shifted_edges(model: &mut Model, pieces: &[Piece], inner: &Location) -> Vec<Shape> {
    let wire = drawn(model, pieces, Transform::IDENTITY, false);
    model
        .ordered_children_of(&wire)
        .unwrap()
        .iter()
        .map(|e| e.moved(inner))
        .collect()
}

fn wire(model: &mut Model, pieces: &[Piece], given: Given, as_edge: bool) -> Shape {
    match given {
        Given::Baked => drawn(model, pieces, motion(), as_edge),
        Given::Placed => {
            let at = location(model, motion());
            drawn(model, pieces, Transform::IDENTITY, as_edge).moved(&at)
        }
        Given::Nested => {
            let (outer, inner) = placement();
            let (outer, inner) = (location(model, outer), location(model, inner));
            if as_edge {
                let edge = drawn(model, pieces, Transform::IDENTITY, true);
                return edge.moved(&inner).moved(&outer);
            }
            let edges = shifted_edges(model, pieces, &inner);
            let wire = ogeom_algo::make_wire(model, &edges, T).unwrap().shape;
            wire.moved(&outer)
        }
    }
}

/// The plane through `origin` square to `normal`, carried by `motion`.
fn plane(origin: Point, normal: Vector, motion: Transform) -> Plane {
    let x = if normal.cross(Vector::new(1.0, 0.0, 0.0)).magnitude() > 0.5 {
        Vector::new(1.0, 0.0, 0.0)
    } else {
        Vector::new(0.0, 1.0, 0.0)
    };
    let frame = Frame::new(
        motion.apply(origin),
        Direction::new(motion.apply_vector(normal), T).unwrap(),
        Direction::new(motion.apply_vector(x), T).unwrap(),
        T,
    )
    .unwrap();
    Plane::new(frame)
}

/// A planar face bounded by `pieces` on the plane through `origin` square
/// to `normal`, with its pcurves on it where `charted`.
fn face(
    model: &mut Model,
    pieces: &[Piece],
    (origin, normal): (Point, Vector),
    given: Given,
    charted: bool,
) -> Shape {
    let on = |model: &mut Model, plane: Plane, edges: Vec<Shape>| {
        let surface: ogeom_geom::SurfaceGeometry =
            ogeom_geom::PlaneSurface::over(plane, (-100.0, 100.0), (-100.0, 100.0))
                .unwrap()
                .into();
        if charted {
            return ogeom_algo::make_face_with_pcurves(model, surface, &[edges], T)
                .unwrap()
                .shape;
        }
        let ring = ogeom_algo::make_wire(model, &edges, T).unwrap().shape;
        ogeom_algo::make_face(model, surface, &[ring], T)
            .unwrap()
            .shape
    };
    let (outer, inner) = placement();
    match given {
        Given::Baked => {
            let ring = drawn(model, pieces, motion(), false);
            let edges = model.ordered_children_of(&ring).unwrap();
            on(model, plane(origin, normal, motion()), edges)
        }
        Given::Placed => {
            let ring = drawn(model, pieces, Transform::IDENTITY, false);
            let edges = model.ordered_children_of(&ring).unwrap();
            let face = on(model, plane(origin, normal, Transform::IDENTITY), edges);
            let at = location(model, motion());
            face.moved(&at)
        }
        Given::Nested => {
            let at = location(model, inner);
            let edges = shifted_edges(model, pieces, &at);
            let face = on(model, plane(origin, normal, inner), edges);
            let at = location(model, outer);
            face.moved(&at)
        }
    }
}

/// Volume, area or length; the bounds; and the checker's verdict.
fn measured(model: &Model, shape: &Shape) -> (f64, Point, Point, bool) {
    let has = |kind| {
        !explore(model, shape, Filter::OfType(kind))
            .unwrap()
            .is_empty()
    };
    let deflection = ogeom_mesh::Deflection::default();
    let mass = if has(ShapeType::Solid) {
        ogeom_algo::volume_properties(model, shape, deflection, T)
    } else if has(ShapeType::Face) {
        ogeom_algo::surface_properties(model, shape, deflection, T)
    } else {
        ogeom_algo::linear_properties(model, shape, deflection, T)
    };
    let bound = ogeom_algo::tight_bounds(model, shape, T).unwrap();
    let valid = ogeom_algo::check(model, shape, T).unwrap().is_valid();
    (
        mass.unwrap().mass,
        bound.low().unwrap(),
        bound.high().unwrap(),
        valid,
    )
}

type Build = dyn Fn(&mut Model, Given) -> ogeom_core::OgeomResult<Shape>;

/// Build `case` from placed inputs, nested or whole, and from the same
/// inputs baked in place, and hold the results to one volume, area or
/// length, one bound and a clean check. Both are the same exact geometry,
/// apart from where the motion is applied, so they agree to rounding: a
/// part in a billion.
fn same_as_baked(case: &str, build: &Build) {
    let mut in_place = Model::new();
    let want = build(&mut in_place, Given::Baked).unwrap_or_else(|e| panic!("{case}, baked: {e}"));
    let (m0, low0, high0, valid0) = measured(&in_place, &want);
    assert!(valid0, "{case}: built in place, it checks clean");
    for given in [Given::Placed, Given::Nested] {
        let mut model = Model::new();
        let got = build(&mut model, given).unwrap_or_else(|e| panic!("{case} {given:?}: {e}"));
        let (m, low, high, valid) = measured(&model, &got);
        assert!(valid, "{case} {given:?}: it checks clean");
        assert!(
            (m - m0).abs() <= 1e-9 * m0.abs().max(1.0),
            "{case} {given:?}: {m} against {m0}"
        );
        let off = low.distance(low0).max(high.distance(high0));
        assert!(
            off <= 1e-9,
            "{case} {given:?}: the bounds stand {off} apart"
        );
    }
}

fn arc() -> Vec<Piece> {
    vec![Piece::Arc(Point::ORIGIN, 5.0, (0.0, FRAC_PI_2))]
}

/// A square of side `side` in the xz plane, `x` from `at`.
fn standing_square(at: f64, side: f64) -> Vec<Piece> {
    polygon(&[
        Point::new(at, 0.0, 0.0),
        Point::new(at + side, 0.0, 0.0),
        Point::new(at + side, 0.0, side),
        Point::new(at, 0.0, side),
    ])
}

/// A square of side `side` in the xy plane, centred on the z axis at `z`.
fn level_square(z: f64, side: f64) -> Vec<Piece> {
    let h = side / 2.0;
    polygon(&[
        Point::new(-h, -h, z),
        Point::new(h, -h, z),
        Point::new(h, h, z),
        Point::new(-h, h, z),
    ])
}

fn up() -> (Point, Vector) {
    (Point::ORIGIN, Vector::new(0.0, 0.0, 1.0))
}

fn across(x: f64) -> (Point, Vector) {
    (Point::new(x, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0))
}

#[test]
fn profiles_and_paths_are_read_where_they_are_placed() {
    same_as_baked("a face given its pcurves", &|m, given| {
        Ok(face(m, &level_square(0.0, 2.0), up(), given, true))
    });
    same_as_baked("a helical sweep", &|m, given| {
        let profile = face(m, &standing_square(2.0, 1.0), across(2.0), given, false);
        let axis = Axis {
            location: motion().apply(Point::ORIGIN),
            direction: Direction::new(motion().apply_vector(Vector::new(0.0, 0.0, 1.0)), T)
                .unwrap(),
        };
        ogeom_offset::make_helical_sweep(m, &profile, axis, 3.0, 1.5, false, 0.0, T)
            .map(|b| b.shape)
    });
    same_as_baked("a loft to a point", &|m, given| {
        let section = wire(m, &level_square(0.0, 2.0), given, false);
        let apex = Point::new(0.0, 0.0, 3.0);
        let apex = if given == Given::Baked {
            ogeom_algo::make_vertex(m, motion().apply(apex)).shape
        } else {
            let at = location(m, motion());
            ogeom_algo::make_vertex(m, apex).shape.moved(&at)
        };
        ogeom_offset::make_loft(m, &section, &apex, T).map(|b| b.shape)
    });
    same_as_baked("a multisection pipe", &|m, given| {
        let spine = wire(m, &arc(), given, true);
        let start = polygon(&[
            Point::new(4.5, 0.0, -0.5),
            Point::new(5.5, 0.0, -0.5),
            Point::new(5.5, 0.0, 0.5),
            Point::new(4.5, 0.0, 0.5),
        ]);
        let end = polygon(&[
            Point::new(0.0, 4.75, -0.25),
            Point::new(0.0, 5.25, -0.25),
            Point::new(0.0, 5.25, 0.25),
            Point::new(0.0, 4.75, 0.25),
        ]);
        let sections = [wire(m, &start, given, false), wire(m, &end, given, false)];
        ogeom_offset::make_pipe_sections(m, &sections, &spine, false, 1e-3, T).map(|b| b.shape)
    });
    same_as_baked("an evolved sweep", &|m, given| {
        let spine = wire(m, &level_square(0.0, 10.0), given, false);
        let profile = polygon(&[
            Point::new(-5.0, -6.5, 0.0),
            Point::new(-5.0, -5.5, 0.0),
            Point::new(-5.0, -5.5, 1.0),
            Point::new(-5.0, -6.5, 1.0),
        ]);
        let profile = wire(m, &profile, given, false);
        ogeom_offset::make_evolved(m, &spine, &profile, T).map(|b| b.shape)
    });
    same_as_baked("a profile brought onto a spine", &|m, given| {
        let spine = wire(m, &arc(), given, false);
        let profile = face(m, &level_square(0.0, 1.0), up(), given, false);
        let placement = PipePlacement {
            contact: true,
            correction: true,
        };
        ogeom_offset::place_on_spine(m, &profile, &spine, placement, T).map(|b| b.shape)
    });
    same_as_baked("a filling", &|m, given| {
        let sides = polygon(&[
            Point::ORIGIN,
            Point::new(4.0, 0.0, 0.5),
            Point::new(4.0, 4.0, 0.0),
            Point::new(0.0, 4.0, 1.0),
        ]);
        let ring = wire(m, &sides, given, false);
        let edges: [Shape; 4] = model_edges(m, &ring).try_into().unwrap();
        ogeom_offset::make_filling(m, &edges, 16, 1e-3, T).map(|b| b.shape)
    });
    same_as_baked("a wire offset", &|m, given| {
        let ring = vec![
            Piece::Line(Point::ORIGIN, Point::new(4.0, 0.0, 0.0)),
            Piece::Arc(Point::new(4.0, 2.0, 0.0), 2.0, (-FRAC_PI_2, FRAC_PI_2)),
            Piece::Line(Point::new(4.0, 4.0, 0.0), Point::new(0.0, 4.0, 0.0)),
            Piece::Line(Point::new(0.0, 4.0, 0.0), Point::ORIGIN),
        ];
        let ring = wire(m, &ring, given, false);
        ogeom_offset::offset_wire(m, &ring, 0.5, Join::Arc, T).map(|b| b.shape)
    });
    same_as_baked("a tapered prism", &|m, given| {
        let profile = face(m, &level_square(0.0, 2.0), up(), given, false);
        let travel = motion().apply_vector(Vector::new(0.0, 0.0, 3.0));
        ogeom_algo::make_prism_tapered(m, &profile, travel, 0.1, T).map(|b| b.shape)
    });
}

fn model_edges(model: &Model, wire: &Shape) -> Vec<Shape> {
    model.ordered_children_of(wire).unwrap()
}

#[test]
fn a_prism_of_a_mirrored_profile_holds_its_volume() {
    // A 2 by 2 square mirrored in a skew plane, swept 3 along its normal:
    // 12, whichever way the mirror turns its faces.
    let mut model = Model::new();
    let square = face(
        &mut model,
        &level_square(0.0, 2.0),
        up(),
        Given::Baked,
        true,
    );
    let mirror = Transform::plane_mirror(
        Point::ORIGIN,
        Direction::new(Vector::new(1.0, 0.3, 0.2), T).unwrap(),
    );
    let mirrored = ogeom_algo::transformed(&mut model, &square, mirror)
        .unwrap()
        .shape;
    let travel = (mirror * motion()).apply_vector(Vector::new(0.0, 0.0, 3.0));
    let prism = ogeom_algo::make_prism(&mut model, &mirrored, travel, T)
        .unwrap()
        .shape;
    let (volume, _, _, valid) = measured(&model, &prism);
    assert!(valid, "the mirrored prism checks clean");
    // The exact integral, to rounding.
    assert!((volume - 12.0).abs() <= 1e-9 * 12.0, "volume {volume}");
}
