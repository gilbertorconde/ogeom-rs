//! Curves in the plane — pcurves.
//!
//! These carry a curve through a surface's `(u, v)` parameter space. An edge
//! holds one per adjacent face (`docs/DATA_MODEL.md` §6), and boolean face
//! splitting happens entirely in this space: without a pcurve on each face
//! there is nothing to split *with*.
//!
//! A separate type hierarchy from [`crate::curve`] rather than a generic
//! parameter, because a pcurve is used differently from a spatial curve.
//! Distance in parameter space is not distance in space — the same parametric
//! step covers a metre near a cylinder's equator and nothing at all near a
//! sphere's pole — so a function that treats the two alike is wrong, and
//! separate types keep that from compiling.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{
    Axis2, Circle2, Direction2, Ellipse2, KnotVector, Point2, Transform2, Vector2, Weighted,
    bspline, elementary,
};

use crate::traits::{Curve2d, CurveKind, Reversible};

const TAU: f64 = core::f64::consts::TAU;

/// A curve in the plane.
#[derive(Debug, Clone, PartialEq)]
pub enum PlanarCurve {
    /// A straight line.
    Line(Line2d),
    /// A circle or arc.
    Circle(Circle2d),
    /// An ellipse or arc.
    Ellipse(Ellipse2d),
    /// A polynomial or rational B-spline.
    BSpline(BSpline2d),
    /// Another planar curve restricted to a sub-interval.
    Trimmed(Box<Trimmed2d>),
}

/// A straight line in the plane, parameterized by length.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Line2d {
    axis: Axis2,
    domain: (f64, f64),
}

/// A circle in the plane, parameterized by angle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circle2d {
    circle: Circle2,
    reversed: bool,
}

/// An ellipse in the plane, parameterized by eccentric angle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipse2d {
    ellipse: Ellipse2,
    reversed: bool,
}

/// A B-spline in the plane, polynomial or rational.
#[derive(Debug, Clone, PartialEq)]
pub struct BSpline2d {
    knots: KnotVector,
    control: Vec<Weighted<Point2>>,
    rational: bool,
}

/// Another planar curve restricted to a sub-interval.
#[derive(Debug, Clone, PartialEq)]
pub struct Trimmed2d {
    basis: PlanarCurve,
    domain: (f64, f64),
    reversed: bool,
}

/// Reverse a parameter within `[a, b]`, preserving the interval.
fn mirror(u: f64, a: f64, b: f64) -> f64 {
    a + b - u
}

impl Line2d {
    /// A segment between two distinct points, parameterized by arc length.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the points
    /// coincide.
    pub fn segment(from: Point2, to: Point2, tol: Tolerances) -> OgeomResult<Self> {
        Ok(Self {
            axis: Axis2::through(from, to, tol)?,
            domain: (0.0, from.distance(to)),
        })
    }

    /// A line over an explicit parameter range.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the range
    /// is empty or non-finite.
    pub fn over(axis: Axis2, start: f64, end: f64) -> OgeomResult<Self> {
        if !start.is_finite() || !end.is_finite() || end <= start {
            ogeom_bail!(Construction, "line range [{start}, {end}] is empty");
        }
        Ok(Self {
            axis,
            domain: (start, end),
        })
    }

    /// The underlying axis.
    #[must_use]
    pub const fn axis(&self) -> Axis2 {
        self.axis
    }
}

impl Circle2d {
    /// A full circle.
    #[must_use]
    pub const fn new(circle: Circle2) -> Self {
        Self {
            circle,
            reversed: false,
        }
    }

    /// The underlying circle.
    #[must_use]
    pub const fn circle(&self) -> Circle2 {
        self.circle
    }

    /// Whether the curve runs backwards along its underlying circle.
    ///
    /// Part of the curve's state and not derivable from its circle, so
    /// anything that has to reproduce this curve exactly — the native format
    /// above all — needs to be able to read it.
    #[must_use]
    pub const fn is_reversed(&self) -> bool {
        self.reversed
    }
}

impl Ellipse2d {
    /// A full ellipse.
    #[must_use]
    pub const fn new(ellipse: Ellipse2) -> Self {
        Self {
            ellipse,
            reversed: false,
        }
    }

    /// The underlying ellipse.
    #[must_use]
    pub const fn ellipse(&self) -> Ellipse2 {
        self.ellipse
    }

    /// Whether the curve runs backwards along its underlying ellipse.
    ///
    /// Part of the curve's state and not derivable from its ellipse, so
    /// anything that has to reproduce this curve exactly — the native format
    /// above all — needs to be able to read it.
    #[must_use]
    pub const fn is_reversed(&self) -> bool {
        self.reversed
    }
}

impl BSpline2d {
    /// A polynomial B-spline.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) on a shape mismatch.
    pub fn new(knots: KnotVector, control: Vec<Point2>, tol: Tolerances) -> OgeomResult<Self> {
        let weighted = control
            .into_iter()
            .map(|p| Weighted::new(p, 1.0, tol))
            .collect::<OgeomResult<Vec<_>>>()?;
        Self::rational(knots, weighted)
    }

    /// A rational B-spline.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) on a shape mismatch.
    pub fn rational(knots: KnotVector, control: Vec<Weighted<Point2>>) -> OgeomResult<Self> {
        if control.len() != knots.control_point_count() {
            ogeom_bail!(
                Dimension,
                "knot vector describes {} control points, got {}",
                knots.control_point_count(),
                control.len()
            );
        }
        let first = control[0].weight;
        let rational = control
            .iter()
            .any(|w| (w.weight - first).abs() > 1e-12 * first.abs());
        Ok(Self {
            knots,
            control,
            rational,
        })
    }

    /// The knot vector.
    #[must_use]
    pub const fn knots(&self) -> &KnotVector {
        &self.knots
    }

    /// The weighted control points.
    #[must_use]
    pub fn control_points(&self) -> &[Weighted<Point2>] {
        &self.control
    }

    /// Whether the weights differ.
    #[must_use]
    pub const fn is_rational(&self) -> bool {
        self.rational
    }
}

impl Trimmed2d {
    /// Restrict `basis` to `[start, end]`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if the range is empty or
    /// leaves the basis curve's domain.
    pub fn new(basis: PlanarCurve, start: f64, end: f64, tol: Tolerances) -> OgeomResult<Self> {
        let (a, b) = basis.domain();
        if !start.is_finite() || !end.is_finite() || end <= start + tol.parametric() {
            ogeom_bail!(Domain, "trim range [{start}, {end}] is empty");
        }
        if !basis.is_periodic() && (start < a - tol.parametric() || end > b + tol.parametric()) {
            ogeom_bail!(Domain, "trim range [{start}, {end}] leaves [{a}, {b}]");
        }
        Ok(Self {
            basis,
            domain: (start, end),
            reversed: false,
        })
    }

    /// The curve being trimmed.
    #[must_use]
    pub const fn basis(&self) -> &PlanarCurve {
        &self.basis
    }

    /// Whether the curve runs backwards along its underlying curve.
    ///
    /// Part of the curve's state and not derivable from its basis curve, so
    /// anything that has to reproduce this curve exactly — the native format
    /// above all — needs to be able to read it.
    #[must_use]
    pub const fn is_reversed(&self) -> bool {
        self.reversed
    }

    /// This curve's parameter mapped onto the basis curve's.
    fn basis_parameter(&self, u: f64, tol: Tolerances) -> OgeomResult<f64> {
        let u = clamp_to_domain(u, self.domain, false, tol)?;
        Ok(if self.reversed {
            mirror(u, self.domain.0, self.domain.1)
        } else {
            u
        })
    }
}

/// Bring `u` into `domain`, wrapping if periodic.
fn clamp_to_domain(
    u: f64,
    domain: (f64, f64),
    periodic: bool,
    tol: Tolerances,
) -> OgeomResult<f64> {
    let (a, b) = domain;
    if periodic {
        return Ok(a + (u - a).rem_euclid(b - a));
    }
    if !u.is_finite() || u < a - tol.parametric() || u > b + tol.parametric() {
        ogeom_bail!(Domain, "parameter {u} outside curve domain [{a}, {b}]");
    }
    Ok(u.clamp(a, b))
}

/// Fill a derivative list to `n + 1` entries, padding with zeros.
fn pad(mut out: Vec<Vector2>, n: usize) -> Vec<Vector2> {
    out.resize(n.max(out.len().saturating_sub(1)) + 1, Vector2::ZERO);
    out.truncate(n + 1);
    out
}

impl Curve2d for Line2d {
    fn domain(&self) -> (f64, f64) {
        self.domain
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point2> {
        Ok(self
            .axis
            .point_at(clamp_to_domain(u, self.domain, false, tol)?))
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector2> {
        clamp_to_domain(u, self.domain, false, tol)?;
        Ok(self.axis.direction.vector())
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector2>> {
        let p = self.point_at(u, tol)?;
        Ok(pad(vec![p.to_vector(), self.axis.direction.vector()], n))
    }

    fn kind(&self) -> CurveKind {
        CurveKind::Line
    }

    fn is_closed(&self, _tol: Tolerances) -> bool {
        false
    }

    fn is_periodic(&self) -> bool {
        false
    }
}

impl Curve2d for Circle2d {
    fn domain(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point2> {
        Ok(Point2::from_vector(self.derivatives_at(u, 0, tol)?[0]))
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector2> {
        Ok(self.derivatives_at(u, 1, tol)?[1])
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector2>> {
        let u = clamp_to_domain(u, self.domain(), true, tol)?;
        let angle = if self.reversed { -u } else { u };
        let f = self.circle.frame();
        let r = self.circle.radius();
        let (sin, cos) = angle.sin_cos();
        let (x, y) = (f.x().vector(), f.y().vector());
        let point = self.circle.centre() + x * (r * cos) + y * (r * sin);
        // Each order of the reversal picks up a factor of -1, so odd orders flip.
        let sign = if self.reversed { -1.0 } else { 1.0 };
        Ok(pad(
            vec![
                point.to_vector(),
                (x * (-r * sin) + y * (r * cos)) * sign,
                x * (-r * cos) + y * (-r * sin),
            ],
            n,
        ))
    }

    fn kind(&self) -> CurveKind {
        CurveKind::Circle
    }

    fn is_closed(&self, _tol: Tolerances) -> bool {
        true
    }

    fn is_periodic(&self) -> bool {
        true
    }
}

impl Curve2d for Ellipse2d {
    fn domain(&self) -> (f64, f64) {
        (0.0, TAU)
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point2> {
        Ok(Point2::from_vector(self.derivatives_at(u, 0, tol)?[0]))
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector2> {
        Ok(self.derivatives_at(u, 1, tol)?[1])
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector2>> {
        let u = clamp_to_domain(u, self.domain(), true, tol)?;
        let angle = if self.reversed { -u } else { u };
        let f = self.ellipse.frame();
        let (a, b) = (self.ellipse.major_radius(), self.ellipse.minor_radius());
        let (sin, cos) = angle.sin_cos();
        let (x, y) = (f.x().vector(), f.y().vector());
        let point = self.ellipse.centre() + x * (a * cos) + y * (b * sin);
        let sign = if self.reversed { -1.0 } else { 1.0 };
        Ok(pad(
            vec![
                point.to_vector(),
                (x * (-a * sin) + y * (b * cos)) * sign,
                x * (-a * cos) + y * (-b * sin),
            ],
            n,
        ))
    }

    fn kind(&self) -> CurveKind {
        CurveKind::Ellipse
    }

    fn is_closed(&self, _tol: Tolerances) -> bool {
        true
    }

    fn is_periodic(&self) -> bool {
        true
    }
}

impl Curve2d for BSpline2d {
    fn domain(&self) -> (f64, f64) {
        self.knots.domain()
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point2> {
        let u = clamp_to_domain(u, self.domain(), false, tol)?;
        if self.rational {
            bspline::evaluate_rational(&self.knots, &self.control, u, tol)
        } else {
            Ok(bspline::evaluate(&self.knots, &self.control, u, tol)?.point())
        }
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector2> {
        Ok(self.derivatives_at(u, 1, tol)?[1])
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector2>> {
        let u = clamp_to_domain(u, self.domain(), false, tol)?;
        let points = bspline::rational_derivatives(&self.knots, &self.control, u, n, tol)?;
        Ok(points.into_iter().map(Point2::to_vector).collect())
    }

    fn kind(&self) -> CurveKind {
        CurveKind::BSpline
    }

    fn is_closed(&self, tol: Tolerances) -> bool {
        self.control[0]
            .point()
            .is_equal(self.control[self.control.len() - 1].point(), tol)
    }

    fn is_periodic(&self) -> bool {
        false
    }
}

impl Curve2d for Trimmed2d {
    fn domain(&self) -> (f64, f64) {
        self.domain
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point2> {
        self.basis.point_at(self.basis_parameter(u, tol)?, tol)
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector2> {
        let d = self.basis.d1_at(self.basis_parameter(u, tol)?, tol)?;
        Ok(if self.reversed { -d } else { d })
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector2>> {
        let t = self.basis_parameter(u, tol)?;
        let mut out = self.basis.derivatives_at(t, n, tol)?;
        if self.reversed {
            for (order, d) in out.iter_mut().enumerate() {
                if order % 2 == 1 {
                    *d = -*d;
                }
            }
        }
        Ok(out)
    }

    fn kind(&self) -> CurveKind {
        CurveKind::Trimmed
    }

    fn is_closed(&self, tol: Tolerances) -> bool {
        match (self.start(tol), self.end(tol)) {
            (Ok(a), Ok(b)) => a.is_equal(b, tol),
            _ => false,
        }
    }

    fn is_periodic(&self) -> bool {
        false
    }
}

macro_rules! dispatch {
    ($self:ident, $c:ident => $body:expr) => {
        match $self {
            Self::Line($c) => $body,
            Self::Circle($c) => $body,
            Self::Ellipse($c) => $body,
            Self::BSpline($c) => $body,
            Self::Trimmed($c) => $body,
        }
    };
}

impl Curve2d for PlanarCurve {
    fn domain(&self) -> (f64, f64) {
        dispatch!(self, c => c.domain())
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point2> {
        dispatch!(self, c => c.point_at(u, tol))
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector2> {
        dispatch!(self, c => c.d1_at(u, tol))
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector2>> {
        dispatch!(self, c => c.derivatives_at(u, n, tol))
    }

    fn kind(&self) -> CurveKind {
        dispatch!(self, c => c.kind())
    }

    fn is_closed(&self, tol: Tolerances) -> bool {
        dispatch!(self, c => c.is_closed(tol))
    }

    fn is_periodic(&self) -> bool {
        dispatch!(self, c => c.is_periodic())
    }
}

impl PlanarCurve {
    /// This curve moved by a planar similarity.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the result
    /// would be degenerate.
    pub fn transformed(&self, t: &Transform2, tol: Tolerances) -> OgeomResult<Self> {
        let scale = t.scale_factor().abs();
        Ok(match self {
            Self::Line(c) => Self::Line(Line2d {
                axis: Axis2::new(
                    t.apply(c.axis.location),
                    Direction2::new(t.apply_vector(c.axis.direction.vector()), tol)?,
                ),
                // A line's parameter is a length, so the domain rescales with it.
                domain: (c.domain.0 * scale, c.domain.1 * scale),
            }),
            Self::Circle(c) => Self::Circle(Circle2d {
                circle: c.circle.transformed(t, tol)?,
                ..*c
            }),
            Self::Ellipse(c) => Self::Ellipse(Ellipse2d {
                ellipse: c.ellipse.transformed(t, tol)?,
                ..*c
            }),
            Self::BSpline(c) => Self::BSpline(BSpline2d {
                control: c
                    .control
                    .iter()
                    .map(|w| Weighted::new(t.apply(w.point()), w.weight, tol))
                    .collect::<OgeomResult<Vec<_>>>()?,
                ..c.clone()
            }),
            Self::Trimmed(c) => Self::Trimmed(Box::new(Trimmed2d {
                basis: c.basis.transformed(t, tol)?,
                domain: if matches!(c.basis, Self::Line(_)) {
                    (c.domain.0 * scale, c.domain.1 * scale)
                } else {
                    c.domain
                },
                reversed: c.reversed,
            })),
        })
    }
}

impl Reversible for PlanarCurve {
    fn reversed(&self) -> Self {
        match self {
            Self::Line(c) => Self::Line(Line2d {
                axis: Axis2::new(
                    c.axis.point_at(c.domain.0 + c.domain.1),
                    c.axis.direction.reversed(),
                ),
                domain: c.domain,
            }),
            Self::Circle(c) => Self::Circle(Circle2d {
                reversed: !c.reversed,
                ..*c
            }),
            Self::Ellipse(c) => Self::Ellipse(Ellipse2d {
                reversed: !c.reversed,
                ..*c
            }),
            Self::BSpline(c) => {
                let (knots, control) = bspline::reverse(&c.knots, &c.control);
                Self::BSpline(BSpline2d {
                    knots,
                    control,
                    ..c.clone()
                })
            }
            Self::Trimmed(c) => Self::Trimmed(Box::new(Trimmed2d {
                reversed: !c.reversed,
                ..(**c).clone()
            })),
        }
    }
}

/// The angle a planar curve's tangent makes with the `u` axis at `t`.
///
/// The natural way to ask which way a pcurve is heading, which is what wire
/// ordering and outer/inner classification are built from.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) at a cusp, where
/// the derivative vanishes.
pub fn tangent_angle(curve: &PlanarCurve, t: f64, tol: Tolerances) -> OgeomResult<f64> {
    let d = curve.d1_at(t, tol)?;
    if d.is_zero(tol) {
        ogeom_bail!(Construction, "curve has no tangent at {t}");
    }
    Ok(elementary::wrap_signed_angle(d.y.atan2(d.x)))
}

impl From<Line2d> for PlanarCurve {
    fn from(c: Line2d) -> Self {
        Self::Line(c)
    }
}
impl From<Circle2d> for PlanarCurve {
    fn from(c: Circle2d) -> Self {
        Self::Circle(c)
    }
}
impl From<Ellipse2d> for PlanarCurve {
    fn from(c: Ellipse2d) -> Self {
        Self::Ellipse(c)
    }
}
impl From<BSpline2d> for PlanarCurve {
    fn from(c: BSpline2d) -> Self {
        Self::BSpline(c)
    }
}
impl From<Trimmed2d> for PlanarCurve {
    fn from(c: Trimmed2d) -> Self {
        Self::Trimmed(Box::new(c))
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ogeom_math::Frame2;

    const T: Tolerances = Tolerances::millimetres();

    fn frame() -> Frame2 {
        Frame2::new(Point2::new(2.0, -1.0), Direction2::from_angle(0.4))
    }

    fn every_curve() -> Vec<PlanarCurve> {
        let spline = {
            let control = vec![
                Point2::new(0.0, 0.0),
                Point2::new(1.0, 2.0),
                Point2::new(3.0, 1.0),
                Point2::new(5.0, 0.0),
                Point2::new(6.0, -1.0),
            ];
            BSpline2d::new(
                KnotVector::clamped_uniform(3, control.len()).unwrap(),
                control,
                T,
            )
            .unwrap()
        };
        vec![
            Line2d::segment(Point2::ORIGIN, Point2::new(3.0, 4.0), T)
                .unwrap()
                .into(),
            Circle2d::new(Circle2::new(frame(), 2.0, T).unwrap()).into(),
            Ellipse2d::new(Ellipse2::new(frame(), 5.0, 3.0, T).unwrap()).into(),
            spline.clone().into(),
            Trimmed2d::new(spline.into(), 0.2, 0.8, T).unwrap().into(),
        ]
    }

    fn interior(c: &PlanarCurve, n: usize) -> Vec<f64> {
        let (a, b) = c.domain();
        (1..n)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let t = i as f64 / n as f64 + 0.0413;
                a + (b - a) * t
            })
            .collect()
    }

    #[test]
    fn derivatives_agree_with_finite_differences() {
        let h = 1e-6;
        for c in every_curve() {
            for u in interior(&c, 8) {
                let d1 = c.d1_at(u, T).unwrap();
                let numeric = (c.point_at(u + h, T).unwrap() - c.point_at(u - h, T).unwrap())
                    * (1.0 / (2.0 * h));
                assert!(
                    (d1 - numeric).magnitude() <= 1e-5 * numeric.magnitude().max(1.0),
                    "{:?} at {u}",
                    c.kind()
                );
            }
        }
    }

    #[test]
    fn derivatives_at_zero_returns_the_point() {
        for c in every_curve() {
            for u in interior(&c, 4) {
                let d = c.derivatives_at(u, 0, T).unwrap();
                assert_eq!(d.len(), 1);
                assert!(Point2::from_vector(d[0]).is_equal(c.point_at(u, T).unwrap(), T));
            }
        }
    }

    #[test]
    fn out_of_domain_parameters_follow_periodicity() {
        for c in every_curve() {
            let (a, b) = c.domain();
            if c.is_periodic() {
                assert!(c.point_at(b + 1.0, T).is_ok(), "{:?}", c.kind());
            } else {
                assert!(c.point_at(b + 1.0, T).is_err(), "{:?}", c.kind());
                assert!(c.point_at(a - 1.0, T).is_err(), "{:?}", c.kind());
            }
        }
    }

    #[test]
    fn reversal_traverses_the_same_points_backwards() {
        for c in every_curve() {
            let r = c.reversed();
            let (a, b) = c.domain();
            assert_eq!(r.domain(), (a, b), "{:?} changed its domain", c.kind());
            for i in 0..=8 {
                let t = f64::from(i) / 8.0;
                let forward = c.point_at(a + (b - a) * t, T).unwrap();
                let backward = r.point_at(a + (b - a) * (1.0 - t), T).unwrap();
                assert!(forward.is_equal(backward, T), "{:?} at {t}", c.kind());
            }
        }
    }

    #[test]
    fn reversing_twice_is_the_identity() {
        for c in every_curve() {
            let twice = c.reversed().reversed();
            for u in interior(&c, 8) {
                assert!(
                    c.point_at(u, T)
                        .unwrap()
                        .is_equal(twice.point_at(u, T).unwrap(), T),
                    "{:?}",
                    c.kind()
                );
            }
        }
    }

    #[test]
    fn a_reversed_curve_heads_the_other_way() {
        for c in every_curve() {
            let r = c.reversed();
            let (a, b) = c.domain();
            let u = a + (b - a) * 0.4;
            let forward = tangent_angle(&c, u, T).unwrap();
            let backward = tangent_angle(&r, mirror(u, a, b), T).unwrap();
            let difference = (forward - backward).abs();
            assert!(
                (difference - core::f64::consts::PI).abs() < 1e-9,
                "{:?}: {forward} vs {backward}",
                c.kind()
            );
        }
    }

    #[test]
    fn transforms_move_curves() {
        let t = Transform2::rotation(Point2::new(1.0, 1.0), 0.7);
        for c in every_curve() {
            let moved = c.transformed(&t, T).unwrap();
            assert_eq!(moved.kind(), c.kind());
            for u in interior(&c, 8) {
                let expected = t.apply(c.point_at(u, T).unwrap());
                assert!(
                    moved.point_at(u, T).unwrap().is_equal(expected, T),
                    "{:?} at {u}",
                    c.kind()
                );
            }
        }
    }

    #[test]
    fn a_line_segments_parameter_is_arc_length() {
        let l = Line2d::segment(Point2::ORIGIN, Point2::new(3.0, 4.0), T).unwrap();
        assert_eq!(l.domain(), (0.0, 5.0));
        assert!(
            l.point_at(2.5, T)
                .unwrap()
                .is_equal(Point2::new(1.5, 2.0), T)
        );
        assert_relative_eq!(l.d1_at(1.0, T).unwrap().magnitude(), 1.0, epsilon = 1e-15);
        assert!(Line2d::segment(Point2::ORIGIN, Point2::ORIGIN, T).is_err());
        assert!(Line2d::over(Axis2::X, 1.0, 1.0).is_err());
    }

    #[test]
    fn a_rational_quadratic_traces_an_exact_arc() {
        let w = core::f64::consts::FRAC_1_SQRT_2;
        let control: Vec<_> = [
            (Point2::new(1.0, 0.0), 1.0),
            (Point2::new(1.0, 1.0), w),
            (Point2::new(0.0, 1.0), 1.0),
        ]
        .iter()
        .map(|(p, w)| Weighted::new(*p, *w, T).unwrap())
        .collect();
        let c = BSpline2d::rational(KnotVector::clamped_uniform(2, 3).unwrap(), control).unwrap();
        assert!(c.is_rational());
        for i in 0..=20 {
            let u = f64::from(i) / 20.0;
            assert_relative_eq!(
                c.point_at(u, T).unwrap().to_vector().magnitude(),
                1.0,
                epsilon = 1e-14
            );
        }
    }

    #[test]
    fn tangent_angle_reports_the_heading() {
        // Along +x, then +y after a quarter turn of the circle.
        let l: PlanarCurve = Line2d::segment(Point2::ORIGIN, Point2::new(5.0, 0.0), T)
            .unwrap()
            .into();
        assert_relative_eq!(tangent_angle(&l, 1.0, T).unwrap(), 0.0, epsilon = 1e-15);

        let c: PlanarCurve =
            Circle2d::new(Circle2::centred(Point2::ORIGIN, 1.0, T).unwrap()).into();
        assert_relative_eq!(
            tangent_angle(&c, 0.0, T).unwrap(),
            core::f64::consts::FRAC_PI_2,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            tangent_angle(&c, core::f64::consts::FRAC_PI_2, T).unwrap(),
            core::f64::consts::PI,
            epsilon = 1e-12
        );
    }

    #[test]
    fn trimming_is_bounds_checked_and_agrees_with_its_basis() {
        let base: PlanarCurve = Line2d::over(Axis2::X, 0.0, 10.0).unwrap().into();
        assert!(Trimmed2d::new(base.clone(), 2.0, 8.0, T).is_ok());
        assert!(Trimmed2d::new(base.clone(), 8.0, 2.0, T).is_err());
        assert!(Trimmed2d::new(base.clone(), -1.0, 5.0, T).is_err());

        let trimmed = Trimmed2d::new(base.clone(), 2.0, 8.0, T).unwrap();
        assert_eq!(trimmed.domain(), (2.0, 8.0));
        for i in 0..=6 {
            let u = 2.0 + 6.0 * f64::from(i) / 6.0;
            assert!(
                trimmed
                    .point_at(u, T)
                    .unwrap()
                    .is_equal(base.point_at(u, T).unwrap(), T)
            );
        }
        assert!(trimmed.point_at(1.0, T).is_err());
    }

    #[test]
    fn a_circle_in_parameter_space_closes_on_itself() {
        let c: PlanarCurve = Circle2d::new(Circle2::new(frame(), 2.0, T).unwrap()).into();
        assert!(c.is_closed(T) && c.is_periodic());
        let base = c.point_at(0.7, T).unwrap();
        for k in [-2.0_f64, 1.0, 3.0] {
            assert!(base.is_equal(c.point_at(k.mul_add(TAU, 0.7), T).unwrap(), T));
        }
    }
}
