//! A planar face whose edges are any curve lying in its plane swept into a
//! solid: the prism and the partial turn keep the profile as their ends, so
//! every edge needs its pcurve on the profile's plane.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, make_edge_between, make_face, make_prism, make_revolution, make_vertex, make_wire,
    volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{
    BSplineCurve, Circle2d, CircleCurve, Curve, Curve3d as _, CurveOnSurface, EllipseCurve,
    HyperbolaCurve, Line2d, LineCurve, OffsetCurve, ParabolaCurve, PlanarCurve, PlaneSurface,
    Reversible as _, SurfaceGeometry, TrimmedCurve, TrimmedSurface,
};
use ogeom::math::{
    Axis, Circle, Circle2, Direction, Direction2, Ellipse, Frame, Frame2, Hyperbola, KnotVector,
    Parabola, Plane, Point, Point2, Vector, Vector2, Weighted,
};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

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
    /// The side itself, a line.
    Line,
    /// A circular arc through both corners, bulging `BULGE` outward.
    Arc,
    /// The same arc, built from the far corner and reversed.
    ReversedArc,
    /// The same arc as a rational quadratic B-spline.
    RationalArc,
    /// The same arc, trimmed from a whole circle and reversed.
    ReversedTrimmedArc,
    /// Half an ellipse on the side, `BULGE` deep.
    HalfEllipse,
    /// A hyperbola arc through both corners, bulging `BULGE` outward, its
    /// branch opening into the diamond.
    Hyperbola,
    /// The same arc, built from the far corner and reversed.
    ReversedHyperbola,
    /// The same arc, trimmed from a longer stretch of the branch.
    TrimmedHyperbola,
    /// The same arc on a branch whose domain runs past it on one side,
    /// built from the far corner and reversed.
    ReversedLongHyperbola,
    /// The straight side as a line drawn on the profile's own plane.
    OnProfilePlane,
    /// The straight side as a line drawn on the profile's plane under a
    /// chart turned and moved against the profile's.
    OnTurnedChart,
    /// The arc drawn as a circle on the profile's plane under the turned
    /// chart, trimmed to the side.
    ArcOnTurnedChart,
    /// The arc drawn as a circle on the profile's plane under a chart whose
    /// normal is turned round, a mirror image of the profile's, and the
    /// plane trimmed.
    ArcOnMirroredChart,
}

/// The hyperbola's half parameter range: the arc is `[-H, H]`.
const H: f64 = 1.0;

/// The profile's plane, which the diamond lies in.
fn profile_plane() -> Plane {
    Plane::new(Frame::new(Point::ORIGIN, normal(), Direction::X, T).unwrap())
}

/// A frame at `origin` whose `x` is `x` and whose `y` is `y`, `x` and `y`
/// perpendicular.
fn frame_xy(origin: Point, x: Vector, y: Vector) -> Frame {
    Frame::new(
        origin,
        Direction::new(x.cross(y), T).unwrap(),
        Direction::new(x, T).unwrap(),
        T,
    )
    .unwrap()
}

/// A conic bulge from `a` to `b`, bulging `BULGE` towards `outward`, of
/// `kind` (an arc, a rational arc, a half ellipse or a hyperbola), its
/// range, and the area between it and the side.
fn bulge(a: Point, b: Point, outward: Vector, kind: Side) -> (Curve, (f64, f64), f64) {
    let d = b - a;
    let half = d.magnitude() / 2.0;
    let along = d / d.magnitude();
    let middle = a + d * 0.5;
    match kind {
        Side::Arc | Side::RationalArc => {
            // The circle through both corners and the apex `BULGE` out:
            // its centre `r - BULGE` inside the side, the arc `2α` wide.
            let r = half.mul_add(half, BULGE * BULGE) / (2.0 * BULGE);
            let alpha = (half / r).asin();
            let centre = middle - outward * (r - BULGE);
            let segment = r * r / 2.0 * (2.0 * alpha - (2.0 * alpha).sin());
            if matches!(kind, Side::RationalArc) {
                // The end tangents meet `r / cos α` out from the centre;
                // the middle weight `cos α` makes the conic a circle.
                let control = vec![
                    Weighted::new(a, 1.0, T).unwrap(),
                    Weighted::new(centre + outward * (r / alpha.cos()), alpha.cos(), T).unwrap(),
                    Weighted::new(b, 1.0, T).unwrap(),
                ];
                let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap();
                let spline = BSplineCurve::rational(knots, control).unwrap();
                return (spline.into(), (0.0, 1.0), segment);
            }
            // Angle `π/2` points out; `x` against the side puts `a` at
            // `π/2 - α`.
            let frame = frame_xy(centre, -along, outward);
            let circle = CircleCurve::new(Circle::new(frame, r, T).unwrap());
            let range = (PI / 2.0 - alpha, PI / 2.0 + alpha);
            (circle.into(), range, segment)
        }
        Side::HalfEllipse => {
            let frame = frame_xy(middle, -along, outward);
            let ellipse = Ellipse::new(frame, half, BULGE, T).unwrap();
            (
                EllipseCurve::new(ellipse).into(),
                (0.0, PI),
                PI * half * BULGE / 2.0,
            )
        }
        Side::Hyperbola => {
            // `x = A cosh t`, `y = B sinh t` with `x` pointing in: the
            // vertex `A` from the centre is the apex, the chord at
            // `A cosh H` is the side. The segment between them is
            // `2 ∫₀ᴴ (A cosh H - A cosh t) B cosh t dt = A B (sinh 2H / 2 - H)`.
            let major = BULGE / (H.cosh() - 1.0);
            let minor = half / H.sinh();
            let centre = middle + outward * (BULGE + major);
            let frame = frame_xy(centre, -outward, along);
            let hyperbola = Hyperbola::new(frame, major, minor, T).unwrap();
            (
                HyperbolaCurve::new(hyperbola, H).unwrap().into(),
                (-H, H),
                major * minor * ((2.0 * H).sinh() / 2.0 - H),
            )
        }
        _ => unreachable!("{kind:?} is not a conic bulge"),
    }
}

/// `curve` over `range` run backwards, over the mirrored range.
fn reversed(curve: &Curve, range: (f64, f64)) -> (Curve, (f64, f64)) {
    let (lo, hi) = curve.domain();
    (curve.reversed(), (lo + hi - range.1, lo + hi - range.0))
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
    let outward = d.cross(n.vector()) / d.magnitude();
    let inset = outward * D;
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
        Side::Line => (
            LineCurve::segment(a, b, T).unwrap().into(),
            (0.0, 2.0 * half),
            0.0,
        ),
        Side::Arc | Side::RationalArc | Side::HalfEllipse | Side::Hyperbola => {
            bulge(a, b, outward, kind)
        }
        Side::ReversedArc | Side::ReversedHyperbola => {
            let forward = if matches!(kind, Side::ReversedArc) {
                Side::Arc
            } else {
                Side::Hyperbola
            };
            let (curve, range, area) = bulge(b, a, outward, forward);
            let (curve, range) = reversed(&curve, range);
            (curve, range, area)
        }
        Side::ReversedTrimmedArc => {
            let (circle, range, area) = bulge(b, a, outward, Side::Arc);
            let trimmed: Curve = TrimmedCurve::new(circle, range.0, range.1, T)
                .unwrap()
                .into();
            let (curve, range) = reversed(&trimmed, range);
            (curve, range, area)
        }
        Side::TrimmedHyperbola => {
            let (curve, _, area) = bulge(a, b, outward, Side::Hyperbola);
            let Curve::Hyperbola(h) = curve else {
                unreachable!()
            };
            let long = HyperbolaCurve::new(h.hyperbola(), 2.0 * H).unwrap();
            (
                TrimmedCurve::new(long.into(), -H, H, T).unwrap().into(),
                (-H, H),
                area,
            )
        }
        Side::ReversedLongHyperbola => {
            let (curve, range, area) = bulge(b, a, outward, Side::Hyperbola);
            let Curve::Hyperbola(h) = curve else {
                unreachable!()
            };
            let long: Curve = HyperbolaCurve::over(h.hyperbola(), -H, 2.0 * H)
                .unwrap()
                .into();
            let (curve, range) = reversed(&long, range);
            (curve, range, area)
        }
        Side::OnProfilePlane | Side::OnTurnedChart => {
            // The chart is the profile's own, or one turned an eighth
            // about the normal and moved within the plane.
            let chart = if matches!(kind, Side::OnProfilePlane) {
                profile_plane().frame()
            } else {
                turned_chart(normal())
            };
            let flat = flat_in(chart);
            let line = PlanarCurve::from(Line2d::segment(flat(a), flat(b), T).unwrap());
            (on_chart(line, chart, false), (0.0, 2.0 * half), 0.0)
        }
        Side::ArcOnTurnedChart | Side::ArcOnMirroredChart => {
            // The arc's circle read in the chart: its centre and axes
            // there, the chart's own `y` turned round with its normal.
            let mirrored = matches!(kind, Side::ArcOnMirroredChart);
            let chart = turned_chart(if mirrored { -normal() } else { normal() });
            let (circle, range, area) = bulge(a, b, outward, Side::Arc);
            let Curve::Circle(circle) = circle else {
                unreachable!()
            };
            let circle = circle.circle();
            let flat = flat_in(chart);
            let flat_direction = |d: Direction| {
                let l = chart.vector_to_local(d.vector());
                Direction2::new(Vector2::new(l.x, l.y), T).unwrap()
            };
            let frame = Frame2::from_axes(
                flat(circle.centre()),
                flat_direction(circle.frame().x()),
                flat_direction(circle.frame().y()),
                T,
            )
            .unwrap();
            let drawn = PlanarCurve::from(Circle2d::new(
                Circle2::new(frame, circle.radius(), T).unwrap(),
            ));
            (on_chart(drawn, chart, mirrored), range, area)
        }
    }
}

/// A chart on the profile's plane with normal `normal`, turned an eighth
/// against the profile's and moved within the plane.
fn turned_chart(normal: Direction) -> Frame {
    Frame::new(
        Point::new(1.0, 0.0, 2.0),
        normal,
        Direction::new(Vector::new(1.0, 0.0, 1.0), T).unwrap(),
        T,
    )
    .unwrap()
}

/// A point's coordinates in `chart`'s `x` and `y`.
fn flat_in(chart: Frame) -> impl Fn(Point) -> Point2 {
    move |p: Point| {
        let l = chart.to_local(p);
        Point2::new(l.x, l.y)
    }
}

/// `drawn` on the plane with `chart`'s frame, that plane trimmed or not.
fn on_chart(drawn: PlanarCurve, chart: Frame, trimmed: bool) -> Curve {
    let plane: SurfaceGeometry =
        PlaneSurface::over(Plane::new(chart), (-10.0, 10.0), (-10.0, 10.0))
            .unwrap()
            .into();
    let surface = if trimmed {
        TrimmedSurface::new(plane, (-8.0, 8.0), (-8.0, 8.0), T)
            .unwrap()
            .into()
    } else {
        plane
    };
    Curve::OnSurface(Box::new(CurveOnSurface::new(drawn, surface)))
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
    let surface = PlaneSurface::over(profile_plane(), (-10.0, 10.0), (-10.0, 10.0)).unwrap();
    let face = make_face(model, surface.into(), &[wire], T).unwrap().shape;
    (face, area)
}

/// The diamond with sides of every curve that lies in a plane, swept 3 along its normal and
/// turned a full turn and three quarters about `z`, is a valid solid of
/// the closed-form volume: the area times 3, and by Pappus the area times
/// the turn times the centre's distance 3 from the axis. The four sides
/// share the diamond's symmetry, so the centroid stays at its centre.
#[test]
fn a_profile_with_sides_of_any_planar_curve_sweeps_into_a_valid_solid() {
    let mut failures = Vec::new();
    for kind in [
        Side::OffsetSpline,
        Side::OffsetSplineTurned,
        Side::OffsetLine,
        Side::OffsetOfTrimmed,
        Side::TrimmedOffset,
        Side::OffsetHalfCircle,
        Side::Parabola,
        Side::Line,
        Side::Arc,
        Side::ReversedArc,
        Side::RationalArc,
        Side::ReversedTrimmedArc,
        Side::HalfEllipse,
        Side::Hyperbola,
        Side::ReversedHyperbola,
        Side::TrimmedHyperbola,
        Side::ReversedLongHyperbola,
        Side::OnProfilePlane,
        Side::OnTurnedChart,
        Side::ArcOnTurnedChart,
        Side::ArcOnMirroredChart,
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
            let built = match built {
                Ok(built) => built.shape,
                Err(e) => {
                    failures.push(format!("{at}: refused: {e}"));
                    continue;
                }
            };
            let found = check(&model, &built, T).unwrap();
            if !found.is_valid() {
                failures.push(format!(
                    "{at}: invalid: {}",
                    found.to_string().lines().next().unwrap_or("")
                ));
                continue;
            }
            let v =
                match volume_properties(&model, &built, Deflection::with_chord(1e-3).unwrap(), T) {
                    Ok(v) => v,
                    Err(e) => {
                        failures.push(format!("{at}: no volume: {e}"));
                        continue;
                    }
                };
            // Every end trim is exact, so no edge or vertex widens past
            // the confusion distance it was built with.
            let widest = [ShapeType::Edge, ShapeType::Vertex]
                .into_iter()
                .flat_map(|kind| explore_unique(&model, &built, kind).unwrap())
                .filter_map(|shape| model.tolerance_of(&shape).unwrap())
                .map(|t| t.get())
                .fold(0.0, f64::max);
            if widest > T.confusion() {
                failures.push(format!("{at}: widened to {widest:.3e}"));
            }
            // Integrated on the exact surfaces, the volume differs from the
            // closed form by quadrature rounding alone.
            if v.deflection != 0.0 || (v.mass - want).abs() > 1e-9 * want {
                failures.push(format!(
                    "{at}: {} against {want}, deflection {}",
                    v.mass, v.deflection
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// A side that leaves the profile's plane is not the profile's boundary,
/// whether or not the plane's closed forms would project it there: sweeping
/// the profile into a solid that keeps it as an end refuses, naming how far
/// off the side sits, rather than building an end whose trim misses the
/// side.
#[test]
fn a_side_leaving_the_profile_plane_refuses_to_sweep_into_ends() {
    let at = corners();
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Z,
    };
    for reversed_trim in [true, false] {
        let mut model = Model::new();
        let v: Vec<Shape> = at
            .iter()
            .map(|p| make_vertex(&mut model, *p).shape)
            .collect();
        let mut edges = Vec::new();
        for i in 0..4 {
            let (a, b) = (at[i], at[(i + 1) % 4]);
            let (curve, range) = if i == 0 {
                // The arc, its plane turned about the side so its middle
                // stands `BULGE sin 0.1` off the profile; as it is, or
                // built from the far corner, trimmed and reversed.
                let d = b - a;
                let outward = d.cross(normal().vector()) / d.magnitude();
                let tilted = outward * 0.1_f64.cos() + normal().vector() * 0.1_f64.sin();
                if reversed_trim {
                    let (circle, range, _) = bulge(b, a, tilted, Side::Arc);
                    let trimmed: Curve = TrimmedCurve::new(circle, range.0, range.1, T)
                        .unwrap()
                        .into();
                    reversed(&trimmed, range)
                } else {
                    let (circle, range, _) = bulge(a, b, tilted, Side::Arc);
                    (circle, range)
                }
            } else {
                let (curve, range, _) = side(a, b, Side::Line);
                (curve, range)
            };
            let edge =
                make_edge_between(&mut model, curve, range, &v[i], &v[(i + 1) % 4], T).unwrap();
            edges.push(edge.shape);
        }
        let wire = make_wire(&mut model, &edges, T).unwrap().shape;
        let surface = PlaneSurface::over(profile_plane(), (-10.0, 10.0), (-10.0, 10.0)).unwrap();
        let face = make_face(&mut model, surface.into(), &[wire], T)
            .unwrap()
            .shape;
        let prism = make_prism(&mut model, &face, normal().vector() * 3.0, T);
        let turned = make_revolution(&mut model, &face, axis, 1.5 * PI, T);
        for built in [prism, turned] {
            let e = built.expect_err("the sweep refuses");
            assert!(
                e.to_string().contains("sits 3.99e-2 off its surface"),
                "reversed trim {reversed_trim}: {e}"
            );
        }
    }
}

/// The hyperbola-sided prism and partial turn, their ends trimmed by the
/// branch in the plane's chart, cross the native, `.brep` and STEP formats
/// and come back valid, of the same closed-form volume; the native text
/// reads back to itself. STEP states no pcurve for the branch, and the
/// reader charts the edge again.
#[test]
fn a_hyperbola_sided_sweep_crosses_the_exchange_formats() {
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Z,
    };
    for turn in [None, Some(1.5 * PI)] {
        let mut model = Model::new();
        let (face, area) = diamond(&mut model, Side::Hyperbola);
        let (built, want) = match turn {
            None => (
                make_prism(&mut model, &face, normal().vector() * 3.0, T),
                area * 3.0,
            ),
            Some(angle) => (
                make_revolution(&mut model, &face, axis, angle, T),
                area * angle * 3.0,
            ),
        };
        let built = built.unwrap().shape;
        let text =
            ogeom::io::native::write(&model, std::slice::from_ref(&built), Default::default())
                .unwrap();
        assert!(
            text.contains("hyperbola2"),
            "{turn:?}: the ends carry the branch"
        );
        let (native, roots) = ogeom::io::native::read(&text).unwrap();
        let again = ogeom::io::native::write(&native, &roots, Default::default()).unwrap();
        assert_eq!(
            text, again,
            "{turn:?}: the native text reads back to itself"
        );
        let brep = ogeom::io::brep::write(&model, &built, T).unwrap();
        let (from_brep, shape) = ogeom::io::brep::read(&brep, T).unwrap();
        let mut document = ogeom::doc::Document::over(model);
        document.add_part("swept", built);
        let step = ogeom::io::write_step(&document, T).unwrap();
        let import = ogeom::io::read_step(&step, T).unwrap();
        assert_eq!(import.solids.len(), 1, "{turn:?}");
        for (back, shape) in [
            (&native, &roots[0]),
            (&from_brep, &shape),
            (import.document.model(), &import.solids[0]),
        ] {
            let found = check(back, shape, T).unwrap();
            assert!(found.is_valid(), "{turn:?}: {found}");
            let v =
                volume_properties(back, shape, Deflection::with_chord(1e-3).unwrap(), T).unwrap();
            assert!(
                (v.mass - want).abs() <= 1e-9 * want,
                "{turn:?}: {} against {want}",
                v.mass
            );
        }
    }
}
