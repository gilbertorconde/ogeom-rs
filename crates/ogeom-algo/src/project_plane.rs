//! An edge projected orthogonally onto a plane, as an exact curve in the
//! plane's own coordinates.
//!
//! Orthogonal projection onto a plane is affine, so it keeps every curve
//! whose family is closed under affine maps: a line stays a line (or
//! collapses to a point), a circle or ellipse becomes an ellipse (a circle
//! when its plane is parallel, a segment when it stands square), and a
//! B-spline becomes the B-spline on its projected control points, weights
//! unchanged. Any other curve is fitted, and the fit's error is returned
//! with it.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{BSpline2d, BSplineCurve, Curve, Curve3d as _, Transformable as _};
use ogeom_math::{KnotVector, Plane, Point, Point2, Vector, Vector2, Weighted};
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
    /// A B-spline: exact for a B-spline edge, fitted for any other curve,
    /// with the fit's largest miss in `fit_error`.
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
        Curve::BSpline(spline) => {
            let piece = restricted(spline, range, tol)?;
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
        _ => fitted(curve, range, onto, tol),
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

/// A curve with no closed-form projection, fitted through its projected
/// points to a hundred times the confusion distance.
fn fitted(
    curve: &Curve,
    range: (f64, f64),
    onto: &Onto,
    tol: Tolerances,
) -> OgeomResult<ProjectedCurve> {
    const SAMPLES: u32 = 200;
    let mut points = Vec::with_capacity(SAMPLES as usize + 1);
    for k in 0..=SAMPLES {
        let t = range.0 + (range.1 - range.0) * f64::from(k) / f64::from(SAMPLES);
        let p = onto.point(curve.point_at(t, tol)?);
        points.push(Point::new(p.x, p.y, 0.0));
    }
    let fit = ogeom_geom::fit::fit_points(&points, 3, tol.confusion() * 100.0, tol)?;
    if !fit.met {
        ogeom_bail!(
            Construction,
            "the projected curve could not be fitted within {}; its best fit misses by {}",
            tol.confusion() * 100.0,
            fit.error
        );
    }
    let flat = BSplineCurve::new(
        fit.curve.knots().clone(),
        fit.curve
            .control_points()
            .iter()
            .map(|c| c.scaled)
            .collect(),
        tol,
    )?;
    let onto_xy = Onto {
        origin: Point::ORIGIN,
        x: Vector::X,
        y: Vector::Y,
    };
    Ok(ProjectedCurve::BSpline {
        curve: projected_spline(&flat, &onto_xy, tol)?,
        fit_error: Some(fit.error),
    })
}
