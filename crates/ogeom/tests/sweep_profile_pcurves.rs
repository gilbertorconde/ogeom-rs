//! A planar face whose edges are offset curves or parabolas swept into a
//! solid: the prism and the partial turn keep the profile as their ends, so
//! every edge needs its pcurve on the profile's plane.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, make_edge_between, make_face, make_prism, make_revolution, make_vertex, make_wire,
    volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{
    BSplineCurve, CircleCurve, Curve, LineCurve, OffsetCurve, ParabolaCurve, PlaneSurface,
    TrimmedCurve,
};
use ogeom::math::{Axis, Circle, Direction, Frame, KnotVector, Parabola, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// The diamond's half diagonal: its sides are 2 long, its area `2 R²`, 4.
const R: f64 = core::f64::consts::SQRT_2;

/// The diamond's centre, 3 from the `z` axis it is turned about.
const CENTRE: Point = Point::new(3.0, 0.0, 0.0);

/// How far each offset side stands from its basis.
const D: f64 = 0.25;

/// How far a parabola side bulges out of the diamond.
const BULGE: f64 = 0.4;

/// The profile plane's normal; the diamond runs anticlockwise about it.
fn normal() -> Direction {
    Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap()
}

/// The diamond's corners in the plane `y = 0`.
fn corners() -> [Point; 4] {
    [0.0, 1.0, 2.0, 3.0].map(|k: f64| {
        let a = PI / 2.0 * k;
        CENTRE + Vector::new(R * a.cos(), 0.0, R * a.sin())
    })
}

/// The clamped cubic B-spline whose four control points lie evenly from
/// `a` to `b`: the segment, on `[0, 1]`.
fn straight_spline(a: Point, b: Point) -> Curve {
    let control: Vec<Point> = (0..4).map(|k| a + (b - a) * (f64::from(k) / 3.0)).collect();
    let knots = KnotVector::clamped_uniform(3, 4).unwrap();
    Curve::BSpline(BSplineCurve::new(knots, control, T).unwrap())
}

/// How a side from one corner to the next is spelled.
#[derive(Debug, Clone, Copy)]
enum Side {
    /// A B-spline line `D` inside the side, offset onto it.
    OffsetSpline,
    /// The same, its reference and distance both turned round.
    OffsetSplineTurned,
    /// A line `D` inside the side, offset onto it.
    OffsetLine,
    /// The middle half of a B-spline line twice the side's length, `D`
    /// inside it, trimmed, then offset onto the side.
    OffsetOfTrimmed,
    /// The same B-spline line offset onto the side, then trimmed.
    TrimmedOffset,
    /// Half a circle of radius `1 - D` on the side, offset outward to the
    /// half disc on the side.
    OffsetHalfCircle,
    /// A parabola through both corners, bulging `BULGE` outward.
    Parabola,
}

/// The side from `a` to `b`, its range, and the area it adds outside the
/// straight side.
fn side(a: Point, b: Point, kind: Side) -> (Curve, (f64, f64), f64) {
    let n = normal();
    let d = b - a;
    let half = d.magnitude() / 2.0;
    let middle = a + d * 0.5;
    // Anticlockwise about the normal, the tangent crossed with it points
    // out of the diamond: an offset along `+D` with the normal as its
    // reference moves a basis `D` inside the side onto it.
    let inset = d.cross(n.vector()) * (D / d.magnitude());
    let offset = |basis: Curve| Curve::Offset(Box::new(OffsetCurve::new(basis, D, n).unwrap()));
    let long_spline = || straight_spline(a - inset - d * 0.5, b - inset + d * 0.5);
    match kind {
        Side::OffsetSpline => (
            offset(straight_spline(a - inset, b - inset)),
            (0.0, 1.0),
            0.0,
        ),
        Side::OffsetSplineTurned => (
            Curve::Offset(Box::new(
                OffsetCurve::new(straight_spline(a - inset, b - inset), -D, -n).unwrap(),
            )),
            (0.0, 1.0),
            0.0,
        ),
        Side::OffsetLine => (
            offset(LineCurve::segment(a - inset, b - inset, T).unwrap().into()),
            (0.0, 2.0 * half),
            0.0,
        ),
        Side::OffsetOfTrimmed => (
            offset(
                TrimmedCurve::new(long_spline(), 0.25, 0.75, T)
                    .unwrap()
                    .into(),
            ),
            (0.25, 0.75),
            0.0,
        ),
        Side::TrimmedOffset => (
            TrimmedCurve::new(offset(long_spline()), 0.25, 0.75, T)
                .unwrap()
                .into(),
            (0.25, 0.75),
            0.0,
        ),
        Side::OffsetHalfCircle => {
            // From `a` at angle 0 to `b` at a half turn, through the
            // frame's `y`, the normal crossed with the side's reverse:
            // outward.
            let frame = Frame::new(middle, n, Direction::new(a - b, T).unwrap(), T).unwrap();
            let circle = CircleCurve::new(Circle::new(frame, half - D, T).unwrap());
            (offset(circle.into()), (0.0, PI), PI * half * half / 2.0)
        }
        Side::Parabola => {
            // Opening inward from its apex `BULGE` outside the side's
            // middle, through both corners at `±half`; the segment it
            // closes off with the side is two thirds of the box on it.
            let inward = Direction::new(n.vector().cross(d), T).unwrap();
            let apex = middle - inward.vector() * BULGE;
            let frame = Frame::new(apex, -n, inward, T).unwrap();
            let parabola = Parabola::new(frame, half * half / (4.0 * BULGE), T).unwrap();
            (
                ParabolaCurve::new(parabola, half).unwrap().into(),
                (-half, half),
                4.0 / 3.0 * half * BULGE,
            )
        }
    }
}

/// The diamond face with sides of `kind`, and its area.
fn diamond(model: &mut Model, kind: Side) -> (Shape, f64) {
    let at = corners();
    let v: Vec<Shape> = at.iter().map(|p| make_vertex(model, *p).shape).collect();
    let mut edges = Vec::new();
    let mut area = 2.0 * R * R;
    for i in 0..4 {
        let (curve, range, outside) = side(at[i], at[(i + 1) % 4], kind);
        area += outside;
        let edge = make_edge_between(model, curve, range, &v[i], &v[(i + 1) % 4], T).unwrap();
        edges.push(edge.shape);
    }
    let wire = make_wire(model, &edges, T).unwrap().shape;
    let plane = Plane::new(Frame::new(Point::ORIGIN, normal(), Direction::X, T).unwrap());
    let surface = PlaneSurface::over(plane, (-10.0, 10.0), (-10.0, 10.0)).unwrap();
    let face = make_face(model, surface.into(), &[wire], T).unwrap().shape;
    (face, area)
}

/// The diamond with offset or parabola sides, swept 3 along its normal and
/// turned a full turn and three quarters about `z`, is a valid solid of
/// the closed-form volume: the area times 3, and by Pappus the area times
/// the turn times the centre's distance 3 from the axis. The four sides
/// share the diamond's symmetry, so the centroid stays at its centre.
#[test]
fn a_profile_with_offset_or_parabola_sides_sweeps_into_a_valid_solid() {
    for kind in [
        Side::OffsetSpline,
        Side::OffsetSplineTurned,
        Side::OffsetLine,
        Side::OffsetOfTrimmed,
        Side::TrimmedOffset,
        Side::OffsetHalfCircle,
        Side::Parabola,
    ] {
        for turn in [None, Some(2.0 * PI), Some(1.5 * PI)] {
            let at = format!("{kind:?} sides, turn {turn:?}");
            let mut model = Model::new();
            let (face, area) = diamond(&mut model, kind);
            let (built, want) = match turn {
                None => (
                    make_prism(&mut model, &face, normal().vector() * 3.0, T),
                    area * 3.0,
                ),
                Some(angle) => {
                    let axis = Axis {
                        location: Point::ORIGIN,
                        direction: Direction::Z,
                    };
                    (
                        make_revolution(&mut model, &face, axis, angle, T),
                        area * angle * 3.0,
                    )
                }
            };
            let built = built.unwrap_or_else(|e| panic!("{at}: {e}")).shape;
            let found = check(&model, &built, T).unwrap();
            assert!(found.is_valid(), "{at}: {found}");
            let v = volume_properties(&model, &built, Deflection::with_chord(1e-3).unwrap(), T)
                .unwrap_or_else(|e| panic!("{at}: {e}"));
            assert_eq!(v.deflection, 0.0, "{at}: integrated on the exact surfaces");
            // Integrated on the exact surfaces, the volume differs from the
            // closed form by quadrature rounding alone.
            assert!(
                (v.mass - want).abs() <= 1e-9 * want,
                "{at}: {} against {want}",
                v.mass
            );
        }
    }
}
