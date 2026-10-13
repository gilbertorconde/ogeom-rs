//! An edge projected orthogonally onto a plane, as an exact curve in the
//! plane's own coordinates.
//!
//! Orthogonal projection onto a plane is affine, so it keeps every curve
//! whose family is closed under affine maps: a line stays a line (or
//! collapses to a point), a circle or ellipse becomes an ellipse (a circle
//! when its plane is parallel, a segment when it stands square), a parabola
//! a parabola and a hyperbola a hyperbola (each a segment when its plane
//! stands square), and a B-spline becomes the B-spline on its projected
//! control points, weights unchanged, clamped first where its knots are
//! not. An offset of a line is a line, and an offset of a circle in its
//! own plane a circle, and each projects as one. Any other curve is
//! fitted, and the fit's error is returned with it.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    BSpline2d, BSplineCurve, CircleCurve, Curve, Curve3d as _, OffsetCurve, Transformable as _,
};
use ogeom_math::{
    Circle, Frame, KnotVector, Plane, Point, Point2, Transform, Vector, Vector2, Weighted,
};
use ogeom_topo::{EdgeRepr, Model, NodeData, Shape};

/// An edge's orthogonal projection onto a plane, in the plane frame's `x`
/// and `y` coordinates.
#[derive(Debug, Clone, PartialEq)]
pub enum ProjectedCurve {
    /// The edge projects to one point: a line along the plane's normal.
    Point(Point2),
    /// A segment.
    Line {
        /// One end.
        start: Point2,
        /// The other end.
        end: Point2,
    },
    /// A circle or circular arc, `centre + radius (cos t, sin t)`, turning
    /// counter-clockwise about the plane's normal from `range.0` to
    /// `range.1`; a full turn where the range spans one.
    Circle {
        /// The centre.
        centre: Point2,
        /// The radius.
        radius: f64,
        /// The arc's angles, `range.0 < range.1`.
        range: (f64, f64),
    },
    /// An ellipse or elliptical arc, `centre + major cos t + minor sin t`
    /// with `minor` the major axis turned a quarter counter-clockwise and
    /// scaled by `ratio`, from `range.0` to `range.1`.
    Ellipse {
        /// The centre.
        centre: Point2,
        /// The major semi-axis, as a vector.
        major: Vector2,
        /// The minor semi-axis over the major, in `(0, 1)`.
        ratio: f64,
        /// The arc's eccentric angles, `range.0 < range.1`.
        range: (f64, f64),
    },
    /// A parabola or a piece of one, `vertex + axis s^2 / (4 focal) + side s`
    /// with `axis` a unit vector and `side` it turned a quarter
    /// counter-clockwise, over `s` from `range.0` to `range.1`. `s` is the
    /// edge's own parameter scaled and shifted, increasing with it or
    /// against it.
    Parabola {
        /// The vertex.
        vertex: Point2,
        /// The unit direction the parabola opens in.
        axis: Vector2,
        /// The distance from the vertex to the focus.
        focal: f64,
        /// The piece's parameters, `range.0 < range.1`.
        range: (f64, f64),
    },
    /// A hyperbola's branch or a piece of one, `centre + major cosh s +
    /// minor sinh s` with `major` and `minor` square, over `s` from
    /// `range.0` to `range.1`: the edge's own parameter shifted.
    Hyperbola {
        /// The centre, off the branch.
        centre: Point2,
        /// The major semi-axis: centre to vertex.
        major: Vector2,
        /// The minor semi-axis, square to the major.
        minor: Vector2,
        /// The piece's parameters, `range.0 < range.1`.
        range: (f64, f64),
    },
    /// A B-spline: exact for a B-spline edge, fitted for any other curve at
    /// the edge's own parameter, with the fit's largest miss in
    /// `fit_error`.
    BSpline {
        /// The curve.
        curve: BSpline2d,
        /// `None` where the projection is exact.
        fit_error: Option<f64>,
    },
}

/// The orthogonal projection of an edge onto a plane, in the plane's own
/// 2D coordinates (its frame's `x` and `y` axes). The edge's placement is
/// applied first.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `edge` is not an edge with a curve of its own, or if a curve with no
/// closed-form projection cannot be fitted; plus whatever evaluating the
/// curve refuses.
pub fn project_edge_onto_plane(
    model: &Model,
    edge: &Shape,
    plane: &Plane,
    tol: Tolerances,
) -> OgeomResult<ProjectedCurve> {
    let Some(NodeData::Edge(data)) = model.node(edge).map(|n| n.data()) else {
        ogeom_bail!(Construction, "only an edge projects onto a plane");
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        ogeom_bail!(
            Construction,
            "the edge has no curve in space (a degenerate edge), so there is nothing to project"
        );
    };
    let Some(curve) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "the edge's curve is not in this model");
    };
    let curve = curve
        .clone()
        .transformed(&edge.transform(model.datums())?, tol)?;
    project_curve(&curve, *range, &Onto::new(plane), tol)
}

/// The plane, as the map it is.
struct Onto {
    origin: Point,
    x: Vector,
    y: Vector,
}

impl Onto {
    fn new(plane: &Plane) -> Self {
        let frame = plane.frame();
        Self {
            origin: frame.origin(),
            x: frame.x().vector(),
            y: frame.y().vector(),
        }
    }

    fn point(&self, p: Point) -> Point2 {
        let d = p - self.origin;
        Point2::new(d.dot(self.x), d.dot(self.y))
    }

    fn vector(&self, v: Vector) -> Vector2 {
        Vector2::new(v.dot(self.x), v.dot(self.y))
    }
}

fn project_curve(
    curve: &Curve,
    range: (f64, f64),
    onto: &Onto,
    tol: Tolerances,
) -> OgeomResult<ProjectedCurve> {
    match curve {
        Curve::Line(_) => {
            let start = onto.point(curve.point_at(range.0, tol)?);
            let end = onto.point(curve.point_at(range.1, tol)?);
            Ok(if start.distance(end) <= tol.confusion() {
                ProjectedCurve::Point(start.midpoint(end))
            } else {
                ProjectedCurve::Line { start, end }
            })
        }
        Curve::Circle(c) => {
            let circle = c.circle();
            let frame = circle.frame();
            let r = circle.radius();
            conic(
                onto.point(frame.origin()),
                onto.vector(frame.x().vector() * r),
                onto.vector(frame.y().vector() * r),
                range,
                tol,
            )
        }
        Curve::Ellipse(e) => {
            let ellipse = e.ellipse();
            let frame = ellipse.frame();
            conic(
                onto.point(frame.origin()),
                onto.vector(frame.x().vector() * ellipse.major_radius()),
                onto.vector(frame.y().vector() * ellipse.minor_radius()),
                range,
                tol,
            )
        }
        Curve::Parabola(c) if !c.is_reversed() => {
            let parabola = c.parabola();
            let frame = parabola.frame();
            Ok(parabolic(
                onto.point(parabola.apex()),
                onto.vector(frame.x().vector() / (4.0 * parabola.focal())),
                onto.vector(frame.y().vector()),
                range,
                tol,
            ))
        }
        Curve::Hyperbola(c) if !c.is_reversed() => {
            let hyperbola = c.hyperbola();
            let frame = hyperbola.frame();
            Ok(hyperbolic(
                onto.point(hyperbola.centre()),
                onto.vector(frame.x().vector() * hyperbola.major_radius()),
                onto.vector(frame.y().vector() * hyperbola.minor_radius()),
                range,
                tol,
            ))
        }
        Curve::BSpline(spline) => {
            let piece = clamped(&restricted(spline, range, tol)?, tol)?;
            Ok(ProjectedCurve::BSpline {
                curve: projected_spline(&piece, onto, tol)?,
                fit_error: None,
            })
        }
        Curve::Trimmed(trimmed) => {
            // A trimmed curve's parameter is its basis's, run backwards
            // where it is reversed; the projection takes no side.
            let (s, e) = trimmed.domain();
            let on_basis = if trimmed.is_reversed() {
                (s + e - range.1, s + e - range.0)
            } else {
                range
            };
            project_curve(trimmed.basis(), on_basis, onto, tol)
        }
        Curve::Offset(offset) => match exact_offset(offset, range, tol)? {
            Some(exact) => project_curve(&exact, range, onto, tol),
            None => fitted(curve, range, onto, tol),
        },
        _ => fitted(curve, range, onto, tol),
    }
}

/// An offset curve as the curve it is, at the same parameter, where that
/// has a closed form: a line's offset is the line moved square to it, and
/// a circle's, offset in its own plane, the circle about the same centre
/// with the offset added to its radius. `None` for any other.
fn exact_offset(
    offset: &OffsetCurve,
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<Option<Curve>> {
    match offset.basis() {
        Curve::Line(line) => {
            // The offset direction is the same all along a line.
            let at = f64::midpoint(range.0, range.1);
            let shift = offset.point_at(at, tol)? - line.point_at(at, tol)?;
            Ok(Some(
                Curve::from(*line).transformed(&Transform::translation(shift), tol)?,
            ))
        }
        Curve::Circle(circle) => {
            let circle = circle.circle();
            let frame = circle.frame();
            let normal = frame.z().vector();
            let reference = offset.reference().vector();
            if normal.cross(reference).magnitude() > tol.angular() {
                return Ok(None);
            }
            // Square to the tangent and the circle's normal: outward from
            // the centre with the reference along the normal, inward
            // against it.
            let radius = circle.radius() + offset.distance() * normal.dot(reference).signum();
            if radius.abs() <= tol.confusion() {
                return Ok(None);
            }
            // A radius turned negative runs the circle half a turn round,
            // at the same parameter.
            let frame = if radius > 0.0 {
                frame
            } else {
                Frame::new(frame.origin(), frame.z(), -frame.x(), tol)?
            };
            Ok(Some(
                CircleCurve::new(Circle::new(frame, radius.abs(), tol)?).into(),
            ))
        }
        _ => Ok(None),
    }
}

/// The projection of `centre + u cos t + v sin t` over `range`: an
/// ellipse in general, a circle where `u` and `v` stay square and equal,
/// a segment where the ellipse closes to no width.
fn conic(
    centre: Point2,
    u: Vector2,
    v: Vector2,
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<ProjectedCurve> {
    // Principal axes of the conjugate pair: at `t0` the point is furthest
    // from the centre, and the pair turned by `t0` is square.
    let t0 = 0.5 * (2.0 * u.dot(v)).atan2(u.dot(u) - v.dot(v));
    let major = u * t0.cos() + v * t0.sin();
    let minor = v * t0.cos() - u * t0.sin();
    let (major, minor, t0) = if minor.magnitude() > major.magnitude() {
        (minor, -major, t0 + core::f64::consts::FRAC_PI_2)
    } else {
        (major, minor, t0)
    };
    // The point is `centre + major cos(t - t0) + minor sin(t - t0)`.
    let length = major.magnitude();
    if length <= tol.confusion() {
        return Ok(ProjectedCurve::Point(centre));
    }
    let width = minor.magnitude();
    if width <= tol.confusion() {
        return Ok(segment(centre, major, (range.0 - t0, range.1 - t0)));
    }
    // Counter-clockwise where the minor axis leads the major by a quarter
    // turn; otherwise the parameter runs the other way round.
    let (from, to) = if major.cross(minor) > 0.0 {
        (range.0 - t0, range.1 - t0)
    } else {
        (t0 - range.1, t0 - range.0)
    };
    let ratio = width / length;
    if (1.0 - ratio) * length <= tol.confusion() {
        let angle = major.y.atan2(major.x);
        return Ok(ProjectedCurve::Circle {
            centre,
            radius: length,
            range: normalized(from + angle, to + angle),
        });
    }
    Ok(ProjectedCurve::Ellipse {
        centre,
        major,
        ratio,
        range: normalized(from, to),
    })
}

/// The projection of `origin + p t^2 + q t` over `range`: a parabola, a
/// segment where `q` lies along `p` or `p` vanishes (the parabola seen
/// edge-on), a point where both do.
fn parabolic(
    origin: Point2,
    p: Vector2,
    q: Vector2,
    range: (f64, f64),
    tol: Tolerances,
) -> ProjectedCurve {
    let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
    let reach = (lo.abs().max(hi.abs())).max(1.0);
    // `q` split along `p` and square to it.
    let pp = p.dot(p);
    let alpha = if pp > 0.0 { q.dot(p) / pp } else { 0.0 };
    let across = q - p * alpha;
    let at = |t: f64| origin + p * (t * t) + q * t;
    if p.magnitude() * reach * reach <= tol.confusion()
        || across.magnitude() * reach <= tol.confusion()
    {
        // On one line: the ends and, where it lies inside, the turn.
        let mut points = vec![at(lo), at(hi)];
        if pp > 0.0 {
            let turn = -0.5 * alpha;
            if turn > lo && turn < hi {
                points.push(at(turn));
            }
        }
        return line_through(&points, tol);
    }
    // With `s = w (t + alpha / 2)`, `w = |across|`, the point is
    // `vertex + p s^2 / w^2 + across / w s`.
    let w = across.magnitude();
    let vertex = origin - p * (alpha * alpha / 4.0) - across * (alpha / 2.0);
    let axis = p / p.magnitude();
    let focal = w * w / (4.0 * p.magnitude());
    // `side` the axis turned counter-clockwise; where the parabola runs the
    // other way, `s` runs against `t`.
    let sense = if axis.cross(across) > 0.0 { 1.0 } else { -1.0 };
    let (a, b) = (
        sense * w * (lo + alpha / 2.0),
        sense * w * (hi + alpha / 2.0),
    );
    ProjectedCurve::Parabola {
        vertex,
        axis,
        focal,
        range: (a.min(b), a.max(b)),
    }
}

/// The projection of `centre + u cosh t + v sinh t` over `range`: a
/// hyperbola on its principal axes, a segment where `u` and `v` lie along
/// one line (the hyperbola seen edge-on).
fn hyperbolic(
    centre: Point2,
    u: Vector2,
    v: Vector2,
    range: (f64, f64),
    tol: Tolerances,
) -> ProjectedCurve {
    let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
    let at = |t: f64| centre + u * t.cosh() + v * t.sinh();
    let scale = u.magnitude().max(v.magnitude()) * lo.abs().max(hi.abs()).cosh();
    let squared = u.dot(u) + v.dot(v);
    if u.cross(v).abs() * scale <= tol.confusion() * squared.max(f64::MIN_POSITIVE).sqrt()
        || 2.0 * u.dot(v).abs() >= squared
    {
        // On one line through the centre: the ends and, where the point
        // turns back along it, the turn, where `tanh t` is `-b / a`.
        let mut points = vec![at(lo), at(hi)];
        let along = if u.magnitude() >= v.magnitude() { u } else { v };
        if along.magnitude() > 0.0 {
            let (a, b) = (u.dot(along), v.dot(along));
            if b.abs() < a.abs() {
                let turn = (-b / a).atanh();
                if turn > lo && turn < hi {
                    points.push(at(turn));
                }
            }
        }
        return line_through(&points, tol);
    }
    // Shifted by `t0`, the conjugate pair turns square:
    // `tanh 2 t0 = -2 u.v / (|u|^2 + |v|^2)`.
    let t0 = 0.5 * (-2.0 * u.dot(v) / squared).atanh();
    let major = u * t0.cosh() + v * t0.sinh();
    let minor = u * t0.sinh() + v * t0.cosh();
    ProjectedCurve::Hyperbola {
        centre,
        major,
        minor,
        range: (lo - t0, hi - t0),
    }
}

/// The segment holding points on one line, or the point they all are.
fn line_through(points: &[Point2], tol: Tolerances) -> ProjectedCurve {
    let mut start = points[0];
    let mut end = points[0];
    for &p in points {
        if p.distance(start) > end.distance(start) {
            end = p;
        }
    }
    for &p in points {
        if p.distance(end) > start.distance(end) {
            start = p;
        }
    }
    if start.distance(end) <= tol.confusion() {
        ProjectedCurve::Point(start.midpoint(end))
    } else {
        ProjectedCurve::Line { start, end }
    }
}

/// A range moved by whole turns to start in `[0, 2 pi)`.
fn normalized(from: f64, to: f64) -> (f64, f64) {
    let start = from.rem_euclid(core::f64::consts::TAU);
    (start, start + (to - from))
}

/// The segment `centre + major cos s` covers for `s` over `range`: its
/// ends are the extremes of `cos s`, which an arc through a multiple of pi
/// reaches inside the range.
fn segment(centre: Point2, major: Vector2, range: (f64, f64)) -> ProjectedCurve {
    let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
    let mut values = vec![lo.cos(), hi.cos()];
    let pi = core::f64::consts::PI;
    let mut k = (lo / pi).ceil();
    while k * pi <= hi && values.len() < 4 {
        values.push((k * pi).cos());
        k += 1.0;
    }
    let least = values.iter().copied().fold(f64::INFINITY, f64::min);
    let most = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    ProjectedCurve::Line {
        start: centre + major * least,
        end: centre + major * most,
    }
}

/// The part of a spline an edge covers, split out where the edge covers
/// less than the whole.
fn restricted(
    spline: &BSplineCurve,
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<BSplineCurve> {
    let (lo, hi) = spline.domain();
    let reach = tol.parametric();
    let mut piece = spline.clone();
    if range.0 > lo + reach {
        piece = piece.split_at(range.0, tol)?.1;
    }
    if range.1 < hi - reach {
        piece = piece.split_at(range.1, tol)?.0;
    }
    Ok(piece)
}

/// A spline clamped at both ends of its domain, the curve unchanged; a
/// clamped spline as it is.
fn clamped(spline: &BSplineCurve, tol: Tolerances) -> OgeomResult<BSplineCurve> {
    let knots = spline.knots();
    if knots.is_clamped() {
        return Ok(spline.clone());
    }
    let control: Vec<Weighted<Point>> = spline.control_points().to_vec();
    let (knots, control) = clamped_start(knots, &control, tol)?;
    // The end clamped as the start of the spline run backwards.
    let (knots, control) = reversed(&knots, &control)?;
    let (knots, control) = clamped_start(&knots, &control, tol)?;
    let (knots, control) = reversed(&knots, &control)?;
    if spline.is_rational() {
        BSplineCurve::rational(knots, control)
    } else {
        BSplineCurve::new(knots, control.iter().map(|c| c.scaled).collect(), tol)
    }
}

/// The domain's first knot taken to the degree's multiplicity by knot
/// insertion, and the knots before its first copy dropped with the control
/// points they alone held.
fn clamped_start(
    knots: &KnotVector,
    control: &[Weighted<Point>],
    tol: Tolerances,
) -> OgeomResult<(KnotVector, Vec<Weighted<Point>>)> {
    let p = knots.degree();
    let (lo, _) = knots.domain();
    let (refined, points) = ogeom_math::bspline::insert_knot(
        knots,
        control,
        lo,
        p.saturating_sub(knots.multiplicity_of(lo)),
        tol,
    )?;
    let cut = refined.knots().partition_point(|k| *k < lo);
    let mut start = vec![lo];
    start.extend_from_slice(&refined.knots()[cut..]);
    Ok((
        KnotVector::new(start, p)?,
        points[cut.saturating_sub(1)..].to_vec(),
    ))
}

/// The same spline run backwards: knots mirrored, control points reversed.
fn reversed(
    knots: &KnotVector,
    control: &[Weighted<Point>],
) -> OgeomResult<(KnotVector, Vec<Weighted<Point>>)> {
    let mirrored: Vec<f64> = knots.knots().iter().rev().map(|k| -k).collect();
    Ok((
        KnotVector::new(mirrored, knots.degree())?,
        control.iter().rev().copied().collect(),
    ))
}

/// A spline's projection: the same knots, each control point projected,
/// each weight kept.
fn projected_spline(spline: &BSplineCurve, onto: &Onto, tol: Tolerances) -> OgeomResult<BSpline2d> {
    let knots: KnotVector = spline.knots().clone();
    if spline.is_rational() {
        // The homogeneous point `w p` maps to `w` times the projected point,
        // which is the projection of `w p` less `w` times the origin's.
        let control = spline
            .control_points()
            .iter()
            .map(|c| {
                let d = c.scaled.to_vector() - onto.origin.to_vector() * c.weight;
                Weighted {
                    scaled: Point2::new(d.dot(onto.x), d.dot(onto.y)),
                    weight: c.weight,
                }
            })
            .collect();
        BSpline2d::rational(knots, control)
    } else {
        let control = spline
            .control_points()
            .iter()
            .map(|c| onto.point(c.scaled))
            .collect();
        BSpline2d::new(knots, control, tol)
    }
}

/// A curve with no closed-form projection, fitted at the curve's own
/// parameters through its projected points to a hundred times the
/// confusion distance, sampled at every span of the curve and held to
/// that between the samples as well as at them.
fn fitted(
    curve: &Curve,
    range: (f64, f64),
    onto: &Onto,
    tol: Tolerances,
) -> OgeomResult<ProjectedCurve> {
    const SAMPLES: usize = 200;
    let target = tol.confusion() * 100.0;
    let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
    let stations = crate::traced::stations(curve, (lo, hi), SAMPLES);
    let fit = crate::traced::fit_traced_2d(
        |t| Ok(onto.point(curve.point_at(t, tol)?)),
        &stations,
        3,
        target,
        tol,
    )?;
    if !fit.met {
        ogeom_bail!(
            Construction,
            "the projected curve could not be fitted within {target}; its best fit misses by {}",
            fit.error
        );
    }
    Ok(ProjectedCurve::BSpline {
        curve: fit.curve,
        fit_error: Some(fit.error),
    })
}
