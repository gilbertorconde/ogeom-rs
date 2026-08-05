//! Concrete space curves.
//!
//! Each analytic curve is a thin parameterization over the shape descriptions
//! in `ogeom-math`; the spline curve carries its own control points. All of them
//! are reachable through [`Curve`], an enum rather than a boxed trait object.
//!
//! # Why an enum
//!
//! Curves are stored in their millions in a real model, compared constantly,
//! and eventually serialized. An enum makes each of those cheap: no allocation,
//! no vtable, `Clone` and `PartialEq` derived, and — most usefully — exhaustive
//! matching, so adding a curve type produces a compile error at every site that
//! needs to know rather than a silent fallthrough.
//!
//! Deliberately *not* `#[non_exhaustive]`. Marking it so would force every
//! match outside this crate to carry a wildcard arm, which is exactly the
//! silent fallthrough the enum exists to prevent: a new curve type would then
//! compile everywhere and be mishandled everywhere. The cost is that adding a
//! variant is a breaking change, which for a kernel this size is the right
//! trade — a curve type nobody handles is worse than a version bump.
//!
//! [`CurveKind`] *is* non-exhaustive, because matching on it is for opting into
//! an analytic shortcut and a caller that does not recognise a kind should fall
//! back to the general path rather than fail to compile.
//!
//! The [`Curve3d`] trait is still the interface algorithms are written against;
//! [`Curve`] implements it and forwards.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{
    Axis, Circle, Ellipse, Frame, Hyperbola, KnotVector, Parabola, Point, Transform, Vector,
    Weighted, bspline, elementary,
};

use crate::traits::{Continuity, Curve3d, CurveKind, Reversible, Transformable};

/// How far along a line the default domain reaches either side of its origin.
///
/// A line is unbounded, but every interface here works on a finite interval, so
/// an unbounded curve needs *some* domain. This is far beyond any real model
/// while staying well short of the range where `f64` spacing becomes coarse.
pub const LINE_EXTENT: f64 = 1.0e9;

/// A curve in space.
#[derive(Debug, Clone, PartialEq)]
pub enum Curve {
    /// A straight line.
    Line(LineCurve),
    /// A circle or arc.
    Circle(CircleCurve),
    /// An ellipse or arc.
    Ellipse(EllipseCurve),
    /// One branch of a hyperbola.
    Hyperbola(HyperbolaCurve),
    /// A parabola.
    Parabola(ParabolaCurve),
    /// A polynomial or rational B-spline.
    BSpline(BSplineCurve),
    /// A helix about an axis, the one transcendental the vocabulary keeps.
    Helix(HelixCurve),
    /// Another curve restricted to a sub-interval.
    Trimmed(Box<TrimmedCurve>),
}

/// A straight line, parameterized by length from its origin.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LineCurve {
    axis: Axis,
    domain: (f64, f64),
}

/// A circle, parameterized by angle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CircleCurve {
    circle: Circle,
    reversed: bool,
}

/// An ellipse, parameterized by eccentric angle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EllipseCurve {
    ellipse: Ellipse,
    reversed: bool,
}

/// One branch of a hyperbola.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HyperbolaCurve {
    hyperbola: Hyperbola,
    domain: (f64, f64),
    reversed: bool,
}

/// A parabola.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParabolaCurve {
    parabola: Parabola,
    domain: (f64, f64),
    reversed: bool,
}

/// A helix about its frame's `z`, parameterized by turn angle.
///
/// The point at `t` sits at angle `t` around the axis, radius out along the
/// turned `x`, risen by `pitch·t/2π` along `z` — so one full turn advances
/// exactly one pitch, and a negative pitch winds the other hand. The speed
/// `√(r² + (pitch/2π)²)` is constant, which gives the arc length a closed
/// form. A helix is transcendental: no rational B-spline states it exactly,
/// which is why it is its own type rather than a conversion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HelixCurve {
    frame: Frame,
    radius: f64,
    pitch: f64,
    domain: (f64, f64),
    reversed: bool,
}

/// A B-spline curve, polynomial or rational.
#[derive(Debug, Clone, PartialEq)]
pub struct BSplineCurve {
    knots: KnotVector,
    control: Vec<Weighted<Point>>,
    rational: bool,
    periodic: bool,
}

/// Another curve restricted to a sub-interval of its domain.
#[derive(Debug, Clone, PartialEq)]
pub struct TrimmedCurve {
    basis: Curve,
    domain: (f64, f64),
    reversed: bool,
}

impl LineCurve {
    /// A line along `axis`, spanning [`LINE_EXTENT`] either side of its origin.
    #[must_use]
    pub const fn new(axis: Axis) -> Self {
        Self {
            axis,
            domain: (-LINE_EXTENT, LINE_EXTENT),
        }
    }

    /// A line segment between two distinct points.
    ///
    /// The domain runs from zero to the distance between them, so the parameter
    /// is arc length — which makes every length query along the segment exact.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the points
    /// coincide.
    pub fn segment(from: Point, to: Point, tol: Tolerances) -> OgeomResult<Self> {
        let axis = Axis::through(from, to, tol)?;
        Ok(Self {
            axis,
            domain: (0.0, from.distance(to)),
        })
    }

    /// A line over an explicit parameter range.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the range
    /// is empty or non-finite.
    pub fn over(axis: Axis, start: f64, end: f64) -> OgeomResult<Self> {
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
    pub const fn axis(&self) -> Axis {
        self.axis
    }
}

impl CircleCurve {
    /// A full circle, running counter-clockwise about its frame's `z`.
    #[must_use]
    pub const fn new(circle: Circle) -> Self {
        Self {
            circle,
            reversed: false,
        }
    }

    /// The underlying circle.
    #[must_use]
    pub const fn circle(&self) -> Circle {
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

impl EllipseCurve {
    /// A full ellipse.
    #[must_use]
    pub const fn new(ellipse: Ellipse) -> Self {
        Self {
            ellipse,
            reversed: false,
        }
    }

    /// The underlying ellipse.
    #[must_use]
    pub const fn ellipse(&self) -> Ellipse {
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

impl HyperbolaCurve {
    /// A hyperbola branch over `[-extent, extent]` in its natural parameter.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `extent` is
    /// not finite and positive.
    pub fn new(hyperbola: Hyperbola, extent: f64) -> OgeomResult<Self> {
        if !extent.is_finite() || extent <= 0.0 {
            ogeom_bail!(
                Construction,
                "hyperbola extent {extent} must be finite and positive"
            );
        }
        Ok(Self {
            hyperbola,
            domain: (-extent, extent),
            reversed: false,
        })
    }

    /// A hyperbola branch over an arbitrary `[start, end]` in its natural
    /// parameter.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// interval is not finite and increasing.
    pub fn over(hyperbola: Hyperbola, start: f64, end: f64) -> OgeomResult<Self> {
        if !start.is_finite() || !end.is_finite() || start >= end {
            ogeom_bail!(
                Construction,
                "hyperbola domain [{start}, {end}] must be finite and increasing"
            );
        }
        Ok(Self {
            hyperbola,
            domain: (start, end),
            reversed: false,
        })
    }

    /// The underlying hyperbola.
    #[must_use]
    pub const fn hyperbola(&self) -> Hyperbola {
        self.hyperbola
    }

    /// Whether the curve runs backwards along its underlying hyperbola.
    ///
    /// Part of the curve's state and not derivable from its hyperbola, so
    /// anything that has to reproduce this curve exactly — the native format
    /// above all — needs to be able to read it.
    #[must_use]
    pub const fn is_reversed(&self) -> bool {
        self.reversed
    }
}

impl ParabolaCurve {
    /// A parabola over `[-extent, extent]`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `extent` is
    /// not finite and positive.
    pub fn new(parabola: Parabola, extent: f64) -> OgeomResult<Self> {
        if !extent.is_finite() || extent <= 0.0 {
            ogeom_bail!(
                Construction,
                "parabola extent {extent} must be finite and positive"
            );
        }
        Ok(Self {
            parabola,
            domain: (-extent, extent),
            reversed: false,
        })
    }

    /// A parabola over an arbitrary `[start, end]`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// interval is not finite and increasing.
    pub fn over(parabola: Parabola, start: f64, end: f64) -> OgeomResult<Self> {
        if !start.is_finite() || !end.is_finite() || start >= end {
            ogeom_bail!(
                Construction,
                "parabola domain [{start}, {end}] must be finite and increasing"
            );
        }
        Ok(Self {
            parabola,
            domain: (start, end),
            reversed: false,
        })
    }

    /// The underlying parabola.
    #[must_use]
    pub const fn parabola(&self) -> Parabola {
        self.parabola
    }

    /// Whether the curve runs backwards along its underlying parabola.
    ///
    /// Part of the curve's state and not derivable from its parabola, so
    /// anything that has to reproduce this curve exactly — the native format
    /// above all — needs to be able to read it.
    #[must_use]
    pub const fn is_reversed(&self) -> bool {
        self.reversed
    }
}

impl HelixCurve {
    /// A helix over `turns` full revolutions from angle zero.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
    /// the radius is not finite and positive, the pitch is not finite and
    /// non-zero — a zero pitch is a circle, and there is a type for that —
    /// or `turns` is not finite and positive.
    pub fn new(frame: Frame, radius: f64, pitch: f64, turns: f64) -> OgeomResult<Self> {
        if !turns.is_finite() || turns <= 0.0 {
            ogeom_bail!(Construction, "a helix over {turns} turns is not a curve");
        }
        Self::over(frame, radius, pitch, 0.0, core::f64::consts::TAU * turns)
    }

    /// A helix over an arbitrary increasing angle interval.
    ///
    /// # Errors
    ///
    /// As [`HelixCurve::new`], with the interval checked instead of the turn
    /// count.
    pub fn over(frame: Frame, radius: f64, pitch: f64, start: f64, end: f64) -> OgeomResult<Self> {
        if !radius.is_finite() || radius <= 0.0 {
            ogeom_bail!(
                Construction,
                "helix radius {radius} must be finite and positive"
            );
        }
        if !pitch.is_finite() || pitch == 0.0 {
            ogeom_bail!(
                Construction,
                "helix pitch {pitch} must be finite and non-zero; a zero pitch is a circle"
            );
        }
        if !start.is_finite() || !end.is_finite() || start >= end {
            ogeom_bail!(
                Construction,
                "helix domain [{start}, {end}] must be finite and increasing"
            );
        }
        Ok(Self {
            frame,
            radius,
            pitch,
            domain: (start, end),
            reversed: false,
        })
    }

    /// The frame the helix turns about.
    #[must_use]
    pub const fn frame(&self) -> &Frame {
        &self.frame
    }

    /// The radius.
    #[must_use]
    pub const fn radius(&self) -> f64 {
        self.radius
    }

    /// The advance along the axis per full turn; negative winds left-handed.
    #[must_use]
    pub const fn pitch(&self) -> f64 {
        self.pitch
    }

    /// Whether evaluation runs the domain backwards.
    #[must_use]
    pub const fn is_reversed(&self) -> bool {
        self.reversed
    }

    /// The exact arc length between two parameters: constant speed times the
    /// swept angle.
    #[must_use]
    pub fn arc_length(&self, from: f64, to: f64) -> f64 {
        (to - from).abs() * self.radius.hypot(self.pitch / core::f64::consts::TAU)
    }

    /// Point and first three derivatives at the raw (unreversed) angle.
    fn at(&self, t: f64) -> (Point, Vector, Vector, Vector) {
        let (sin, cos) = t.sin_cos();
        let x = self.frame.x().vector();
        let y = self.frame.y().vector();
        let z = self.frame.z().vector();
        let rise = self.pitch / core::f64::consts::TAU;
        let point = self.frame.origin()
            + x * (self.radius * cos)
            + y * (self.radius * sin)
            + z * (rise * t);
        let d1 = x * (-self.radius * sin) + y * (self.radius * cos) + z * rise;
        let d2 = x * (-self.radius * cos) + y * (-self.radius * sin);
        let d3 = x * (self.radius * sin) + y * (-self.radius * cos);
        (point, d1, d2, d3)
    }
}

impl Curve3d for HelixCurve {
    fn domain(&self) -> (f64, f64) {
        self.domain
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point> {
        let u = self.normalize_parameter(u, tol)?;
        let t = if self.reversed {
            mirror(u, self.domain.0, self.domain.1)
        } else {
            u
        };
        Ok(self.at(t).0)
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector> {
        Ok(self.derivatives_at(u, 1, tol)?[1])
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector>> {
        let u = self.normalize_parameter(u, tol)?;
        let t = if self.reversed {
            mirror(u, self.domain.0, self.domain.1)
        } else {
            u
        };
        let (point, d1, d2, d3) = self.at(t);
        // The chain rule for the reversal: odd orders flip sign.
        let sign = if self.reversed { -1.0 } else { 1.0 };
        let mut out = vec![point.to_vector(), d1 * sign, d2, d3 * sign];
        out.resize(n.max(3) + 1, Vector::ZERO);
        out.truncate(n + 1);
        Ok(out)
    }

    fn kind(&self) -> CurveKind {
        CurveKind::Helix
    }

    fn continuity(&self) -> Continuity {
        Continuity::CInfinity
    }

    fn is_closed(&self, _tol: Tolerances) -> bool {
        false
    }

    fn is_periodic(&self) -> bool {
        false
    }
}

impl BSplineCurve {
    /// A polynomial B-spline from a knot vector and control points.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) if the control point
    /// count disagrees with the knot vector.
    pub fn new(knots: KnotVector, control: Vec<Point>, tol: Tolerances) -> OgeomResult<Self> {
        let weighted = control
            .into_iter()
            .map(|p| Weighted::new(p, 1.0, tol))
            .collect::<OgeomResult<Vec<_>>>()?;
        Self::rational(knots, weighted)
    }

    /// A rational B-spline from a knot vector and weighted control points.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) if the control point
    /// count disagrees with the knot vector.
    pub fn rational(knots: KnotVector, control: Vec<Weighted<Point>>) -> OgeomResult<Self> {
        if control.len() != knots.control_point_count() {
            ogeom_bail!(
                Dimension,
                "knot vector describes {} control points, got {}",
                knots.control_point_count(),
                control.len()
            );
        }
        // A curve whose weights are all equal is polynomial regardless of what
        // that common value is, and saying so lets evaluation skip the divide.
        let first = control[0].weight;
        let rational = control
            .iter()
            .any(|w| (w.weight - first).abs() > 1e-12 * first.abs());
        Ok(Self {
            knots,
            control,
            rational,
            periodic: false,
        })
    }

    /// The knot vector.
    #[must_use]
    pub const fn knots(&self) -> &KnotVector {
        &self.knots
    }

    /// The weighted control points.
    #[must_use]
    pub fn control_points(&self) -> &[Weighted<Point>] {
        &self.control
    }

    /// Whether the weights differ, so the curve is genuinely rational.
    #[must_use]
    pub const fn is_rational(&self) -> bool {
        self.rational
    }

    /// The degree.
    #[must_use]
    pub const fn degree(&self) -> usize {
        self.knots.degree()
    }

    /// Insert a knot without moving the curve.
    ///
    /// # Errors
    ///
    /// As [`bspline::insert_knot`].
    pub fn with_knot_inserted(&self, u: f64, count: usize, tol: Tolerances) -> OgeomResult<Self> {
        let (knots, control) = bspline::insert_knot(&self.knots, &self.control, u, count, tol)?;
        Ok(Self {
            knots,
            control,
            ..self.clone()
        })
    }

    /// Raise the degree without moving the curve.
    ///
    /// # Errors
    ///
    /// As [`bspline::elevate_degree`].
    pub fn elevated(&self, tol: Tolerances) -> OgeomResult<Self> {
        let (knots, control) = bspline::elevate_degree(&self.knots, &self.control, tol)?;
        Ok(Self {
            knots,
            control,
            ..self.clone()
        })
    }

    /// Split into two curves meeting at `u`.
    ///
    /// # Errors
    ///
    /// As [`bspline::split`].
    pub fn split_at(&self, u: f64, tol: Tolerances) -> OgeomResult<(Self, Self)> {
        let ((lk, lc), (rk, rc)) = bspline::split(&self.knots, &self.control, u, tol)?;
        Ok((
            Self {
                knots: lk,
                control: lc,
                ..self.clone()
            },
            Self {
                knots: rk,
                control: rc,
                ..self.clone()
            },
        ))
    }
}

impl TrimmedCurve {
    /// Restrict `basis` to `[start, end]`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if the range is empty or
    /// falls outside the basis curve's own domain.
    pub fn new(basis: Curve, start: f64, end: f64, tol: Tolerances) -> OgeomResult<Self> {
        let (a, b) = basis.domain();
        if !start.is_finite() || !end.is_finite() || end <= start + tol.parametric() {
            ogeom_bail!(Domain, "trim range [{start}, {end}] is empty");
        }
        if !basis.is_periodic() && (start < a - tol.parametric() || end > b + tol.parametric()) {
            ogeom_bail!(
                Domain,
                "trim range [{start}, {end}] leaves the basis domain [{a}, {b}]"
            );
        }
        Ok(Self {
            basis,
            domain: (start, end),
            reversed: false,
        })
    }

    /// The curve being trimmed.
    #[must_use]
    pub const fn basis(&self) -> &Curve {
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
        let u = self.normalize_parameter(u, tol)?;
        Ok(if self.reversed {
            mirror(u, self.domain.0, self.domain.1)
        } else {
            u
        })
    }
}

/// Reverse a parameter within `[a, b]`, so the curve runs the other way over
/// the same interval.
///
/// Preserving the domain matters: trimming ranges elsewhere refer to it, and a
/// reversal that also renumbered the parameters would invalidate them.
fn mirror(u: f64, a: f64, b: f64) -> f64 {
    a + b - u
}

impl Curve3d for LineCurve {
    fn domain(&self) -> (f64, f64) {
        self.domain
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point> {
        let u = self.normalize_parameter(u, tol)?;
        Ok(self.axis.point_at(u))
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector> {
        self.normalize_parameter(u, tol)?;
        Ok(self.axis.direction.vector())
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector>> {
        let p = self.point_at(u, tol)?;
        let mut out = vec![p.to_vector(), self.axis.direction.vector()];
        out.resize(n + 1, Vector::ZERO);
        out.truncate(n + 1);
        Ok(out)
    }

    fn kind(&self) -> CurveKind {
        CurveKind::Line
    }

    fn continuity(&self) -> Continuity {
        Continuity::CInfinity
    }

    fn is_closed(&self, _tol: Tolerances) -> bool {
        false
    }

    fn is_periodic(&self) -> bool {
        false
    }
}

impl Curve3d for CircleCurve {
    fn domain(&self) -> (f64, f64) {
        (0.0, core::f64::consts::TAU)
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point> {
        let u = self.normalize_parameter(u, tol)?;
        let angle = if self.reversed { -u } else { u };
        Ok(elementary::circle_at(&self.circle, angle).point)
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector> {
        Ok(self.derivatives_at(u, 1, tol)?[1])
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector>> {
        let u = self.normalize_parameter(u, tol)?;
        let angle = if self.reversed { -u } else { u };
        let c = elementary::circle_at(&self.circle, angle);
        // The chain rule for the reversal: each derivative picks up a factor of
        // -1 per order, so odd orders flip sign.
        let sign = if self.reversed { -1.0 } else { 1.0 };
        let mut out = vec![c.point.to_vector(), c.d1 * sign, c.d2];
        out.resize(n.max(2) + 1, Vector::ZERO);
        out.truncate(n + 1);
        Ok(out)
    }

    fn kind(&self) -> CurveKind {
        CurveKind::Circle
    }

    fn continuity(&self) -> Continuity {
        Continuity::CInfinity
    }

    fn is_closed(&self, _tol: Tolerances) -> bool {
        true
    }

    fn is_periodic(&self) -> bool {
        true
    }
}

impl Curve3d for EllipseCurve {
    fn domain(&self) -> (f64, f64) {
        (0.0, core::f64::consts::TAU)
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point> {
        let u = self.normalize_parameter(u, tol)?;
        let angle = if self.reversed { -u } else { u };
        Ok(elementary::ellipse_at(&self.ellipse, angle).point)
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector> {
        Ok(self.derivatives_at(u, 1, tol)?[1])
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector>> {
        let u = self.normalize_parameter(u, tol)?;
        let angle = if self.reversed { -u } else { u };
        let c = elementary::ellipse_at(&self.ellipse, angle);
        let sign = if self.reversed { -1.0 } else { 1.0 };
        let mut out = vec![c.point.to_vector(), c.d1 * sign, c.d2];
        out.resize(n.max(2) + 1, Vector::ZERO);
        out.truncate(n + 1);
        Ok(out)
    }

    fn kind(&self) -> CurveKind {
        CurveKind::Ellipse
    }

    fn continuity(&self) -> Continuity {
        Continuity::CInfinity
    }

    fn is_closed(&self, _tol: Tolerances) -> bool {
        true
    }

    fn is_periodic(&self) -> bool {
        true
    }
}

impl Curve3d for HyperbolaCurve {
    fn domain(&self) -> (f64, f64) {
        self.domain
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point> {
        let u = self.normalize_parameter(u, tol)?;
        let t = if self.reversed {
            mirror(u, self.domain.0, self.domain.1)
        } else {
            u
        };
        Ok(elementary::hyperbola_at(&self.hyperbola, t).point)
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector> {
        Ok(self.derivatives_at(u, 1, tol)?[1])
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector>> {
        let u = self.normalize_parameter(u, tol)?;
        let t = if self.reversed {
            mirror(u, self.domain.0, self.domain.1)
        } else {
            u
        };
        let c = elementary::hyperbola_at(&self.hyperbola, t);
        let sign = if self.reversed { -1.0 } else { 1.0 };
        let mut out = vec![c.point.to_vector(), c.d1 * sign, c.d2];
        out.resize(n.max(2) + 1, Vector::ZERO);
        out.truncate(n + 1);
        Ok(out)
    }

    fn kind(&self) -> CurveKind {
        CurveKind::Hyperbola
    }

    fn continuity(&self) -> Continuity {
        Continuity::CInfinity
    }

    fn is_closed(&self, _tol: Tolerances) -> bool {
        false
    }

    fn is_periodic(&self) -> bool {
        false
    }
}

impl Curve3d for ParabolaCurve {
    fn domain(&self) -> (f64, f64) {
        self.domain
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point> {
        let u = self.normalize_parameter(u, tol)?;
        let t = if self.reversed {
            mirror(u, self.domain.0, self.domain.1)
        } else {
            u
        };
        Ok(elementary::parabola_at(&self.parabola, t).point)
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector> {
        Ok(self.derivatives_at(u, 1, tol)?[1])
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector>> {
        let u = self.normalize_parameter(u, tol)?;
        let t = if self.reversed {
            mirror(u, self.domain.0, self.domain.1)
        } else {
            u
        };
        let c = elementary::parabola_at(&self.parabola, t);
        let sign = if self.reversed { -1.0 } else { 1.0 };
        let mut out = vec![c.point.to_vector(), c.d1 * sign, c.d2];
        out.resize(n.max(2) + 1, Vector::ZERO);
        out.truncate(n + 1);
        Ok(out)
    }

    fn kind(&self) -> CurveKind {
        CurveKind::Parabola
    }

    fn continuity(&self) -> Continuity {
        Continuity::CInfinity
    }

    fn is_closed(&self, _tol: Tolerances) -> bool {
        false
    }

    fn is_periodic(&self) -> bool {
        false
    }
}

impl Curve3d for BSplineCurve {
    fn domain(&self) -> (f64, f64) {
        self.knots.domain()
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point> {
        let u = self.normalize_parameter(u, tol)?;
        if self.rational {
            bspline::evaluate_rational(&self.knots, &self.control, u, tol)
        } else {
            // All weights equal: the projection is a no-op up to that common
            // factor, so the cheaper polynomial path is exact here.
            Ok(bspline::evaluate(&self.knots, &self.control, u, tol)?.point())
        }
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector> {
        Ok(self.derivatives_at(u, 1, tol)?[1])
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector>> {
        let u = self.normalize_parameter(u, tol)?;
        let points = bspline::rational_derivatives(&self.knots, &self.control, u, n, tol)?;
        Ok(points.into_iter().map(Point::to_vector).collect())
    }

    fn kind(&self) -> CurveKind {
        CurveKind::BSpline
    }

    /// Continuity across the whole curve.
    ///
    /// A degree-`p` B-spline is `C^(p - m)` at an interior knot of multiplicity
    /// `m`, and the worst interior knot governs the curve. With no interior
    /// knots at all the curve is a single polynomial piece and so genuinely
    /// smooth to every order.
    ///
    /// Higher orders than `C2` report as `C2`, which is the highest
    /// [`Continuity`] names short of `CInfinity`. Reporting `CInfinity` for a
    /// merely-`C3` curve would be a claim that is false, and the distinction
    /// above `C2` is not one any algorithm here asks about.
    fn continuity(&self) -> Continuity {
        let degree = self.knots.degree();
        let (a, b) = self.knots.domain();
        let worst = self
            .knots
            .distinct()
            .into_iter()
            .filter(|(v, _)| *v > a && *v < b)
            .map(|(_, m)| m)
            .max();
        match worst {
            None => Continuity::CInfinity,
            Some(m) => match degree.saturating_sub(m) {
                0 => Continuity::C0,
                1 => Continuity::C1,
                _ => Continuity::C2,
            },
        }
    }

    fn is_closed(&self, tol: Tolerances) -> bool {
        let (first, last) = (self.control[0], self.control[self.control.len() - 1]);
        first.point().is_equal(last.point(), tol)
    }

    fn is_periodic(&self) -> bool {
        self.periodic
    }
}

impl Curve3d for TrimmedCurve {
    fn domain(&self) -> (f64, f64) {
        self.domain
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point> {
        self.basis.point_at(self.basis_parameter(u, tol)?, tol)
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector> {
        let d = self.basis.d1_at(self.basis_parameter(u, tol)?, tol)?;
        Ok(if self.reversed { -d } else { d })
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector>> {
        let t = self.basis_parameter(u, tol)?;
        let mut out = self.basis.derivatives_at(t, n, tol)?;
        if self.reversed {
            // Chain rule for u -> (s + e - u): each order picks up a factor of
            // -1, so odd orders flip.
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

    fn continuity(&self) -> Continuity {
        self.basis.continuity()
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

/// Dispatch a method across every curve variant.
macro_rules! dispatch {
    ($self:ident, $c:ident => $body:expr) => {
        match $self {
            Self::Line($c) => $body,
            Self::Circle($c) => $body,
            Self::Ellipse($c) => $body,
            Self::Hyperbola($c) => $body,
            Self::Parabola($c) => $body,
            Self::BSpline($c) => $body,
            Self::Helix($c) => $body,
            Self::Trimmed($c) => $body,
        }
    };
}

impl Curve3d for Curve {
    fn domain(&self) -> (f64, f64) {
        dispatch!(self, c => c.domain())
    }

    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point> {
        dispatch!(self, c => c.point_at(u, tol))
    }

    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector> {
        dispatch!(self, c => c.d1_at(u, tol))
    }

    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector>> {
        dispatch!(self, c => c.derivatives_at(u, n, tol))
    }

    fn kind(&self) -> CurveKind {
        dispatch!(self, c => c.kind())
    }

    fn continuity(&self) -> Continuity {
        dispatch!(self, c => c.continuity())
    }

    fn is_closed(&self, tol: Tolerances) -> bool {
        dispatch!(self, c => c.is_closed(tol))
    }

    fn is_periodic(&self) -> bool {
        dispatch!(self, c => c.is_periodic())
    }
}

impl Transformable for Curve {
    fn transformed(&self, t: &Transform, tol: Tolerances) -> OgeomResult<Self> {
        Ok(match self {
            Self::Line(c) => Self::Line(LineCurve {
                axis: Axis::new(
                    t.apply(c.axis.location),
                    t.apply_direction(c.axis.direction, tol)?,
                ),
                // The parameter is a length, so a scaling rescales the domain
                // with it — otherwise the trimmed extent would silently change.
                domain: (
                    c.domain.0 * t.scale_factor().abs(),
                    c.domain.1 * t.scale_factor().abs(),
                ),
            }),
            // No parameter flip for a mirror: transforming the frame already
            // carries its x and y axes with it, so evaluating at the same angle
            // lands on the transformed point. Flipping as well would reverse
            // the arc twice.
            Self::Circle(c) => Self::Circle(CircleCurve {
                circle: c.circle.transformed(t, tol)?,
                ..*c
            }),
            Self::Ellipse(c) => Self::Ellipse(EllipseCurve {
                ellipse: c.ellipse.transformed(t, tol)?,
                ..*c
            }),
            Self::Hyperbola(c) => Self::Hyperbola(HyperbolaCurve {
                hyperbola: c.hyperbola.transformed(t, tol)?,
                ..*c
            }),
            Self::Parabola(c) => Self::Parabola(ParabolaCurve {
                parabola: c.parabola.transformed(t, tol)?,
                ..*c
            }),
            Self::BSpline(c) => {
                let control = c
                    .control
                    .iter()
                    .map(|w| Weighted::new(t.apply(w.point()), w.weight, tol))
                    .collect::<OgeomResult<Vec<_>>>()?;
                Self::BSpline(BSplineCurve {
                    control,
                    ..c.clone()
                })
            }
            Self::Helix(c) => Self::Helix(HelixCurve {
                frame: t.apply_frame(&c.frame, tol)?,
                radius: c.radius * t.scale_factor().abs(),
                pitch: c.pitch * t.scale_factor().abs(),
                ..*c
            }),
            Self::Trimmed(c) => Self::Trimmed(Box::new(TrimmedCurve {
                basis: c.basis.transformed(t, tol)?,
                // A line's parameter is a length and rescales; every other
                // curve's is an angle or a spline parameter and does not.
                domain: if matches!(c.basis, Self::Line(_)) {
                    let s = t.scale_factor().abs();
                    (c.domain.0 * s, c.domain.1 * s)
                } else {
                    c.domain
                },
                reversed: c.reversed,
            })),
        })
    }
}

impl Reversible for Curve {
    fn reversed(&self) -> Self {
        match self {
            Self::Line(c) => Self::Line(LineCurve {
                axis: Axis::new(
                    c.axis.point_at(c.domain.0 + c.domain.1),
                    c.axis.direction.reversed(),
                ),
                domain: c.domain,
            }),
            Self::Circle(c) => Self::Circle(CircleCurve {
                reversed: !c.reversed,
                ..*c
            }),
            Self::Ellipse(c) => Self::Ellipse(EllipseCurve {
                reversed: !c.reversed,
                ..*c
            }),
            Self::Hyperbola(c) => Self::Hyperbola(HyperbolaCurve {
                reversed: !c.reversed,
                ..*c
            }),
            Self::Parabola(c) => Self::Parabola(ParabolaCurve {
                reversed: !c.reversed,
                ..*c
            }),
            Self::BSpline(c) => {
                let (knots, control) = bspline::reverse(&c.knots, &c.control);
                Self::BSpline(BSplineCurve {
                    knots,
                    control,
                    ..c.clone()
                })
            }
            Self::Helix(c) => Self::Helix(HelixCurve {
                reversed: !c.reversed,
                ..*c
            }),
            // A flag rather than reversing the basis: mirroring the trim range
            // within the basis domain would move this curve's own domain, and
            // trimming ranges held elsewhere refer to it.
            Self::Trimmed(c) => Self::Trimmed(Box::new(TrimmedCurve {
                reversed: !c.reversed,
                ..(**c).clone()
            })),
        }
    }
}

impl From<LineCurve> for Curve {
    fn from(c: LineCurve) -> Self {
        Self::Line(c)
    }
}
impl From<CircleCurve> for Curve {
    fn from(c: CircleCurve) -> Self {
        Self::Circle(c)
    }
}
impl From<EllipseCurve> for Curve {
    fn from(c: EllipseCurve) -> Self {
        Self::Ellipse(c)
    }
}
impl From<HelixCurve> for Curve {
    fn from(c: HelixCurve) -> Self {
        Self::Helix(c)
    }
}

impl From<HyperbolaCurve> for Curve {
    fn from(c: HyperbolaCurve) -> Self {
        Self::Hyperbola(c)
    }
}
impl From<ParabolaCurve> for Curve {
    fn from(c: ParabolaCurve) -> Self {
        Self::Parabola(c)
    }
}
impl From<BSplineCurve> for Curve {
    fn from(c: BSplineCurve) -> Self {
        Self::BSpline(c)
    }
}
impl From<TrimmedCurve> for Curve {
    fn from(c: TrimmedCurve) -> Self {
        Self::Trimmed(Box::new(c))
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ogeom_math::{Direction, Frame};

    const T: Tolerances = Tolerances::millimetres();

    fn tilted() -> Frame {
        Frame::new(
            Point::new(1.0, -2.0, 3.0),
            Direction::from_coords(1.0, 2.0, 3.0, T).unwrap(),
            Direction::X,
            T,
        )
        .unwrap()
    }

    fn every_curve() -> Vec<Curve> {
        let spline = {
            let control = vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 2.0, 0.0),
                Point::new(3.0, 1.0, 1.0),
                Point::new(5.0, 0.0, 2.0),
                Point::new(6.0, -1.0, 0.0),
            ];
            let knots = KnotVector::clamped_uniform(3, control.len()).unwrap();
            BSplineCurve::new(knots, control, T).unwrap()
        };
        vec![
            LineCurve::segment(Point::ORIGIN, Point::new(3.0, 4.0, 0.0), T)
                .unwrap()
                .into(),
            CircleCurve::new(Circle::new(tilted(), 2.0, T).unwrap()).into(),
            EllipseCurve::new(Ellipse::new(tilted(), 5.0, 3.0, T).unwrap()).into(),
            HyperbolaCurve::new(Hyperbola::new(tilted(), 3.0, 4.0, T).unwrap(), 1.5)
                .unwrap()
                .into(),
            ParabolaCurve::new(Parabola::new(tilted(), 2.0, T).unwrap(), 4.0)
                .unwrap()
                .into(),
            spline.clone().into(),
            HelixCurve::new(tilted(), 2.5, 1.25, 2.0).unwrap().into(),
            TrimmedCurve::new(spline.into(), 0.2, 0.8, T)
                .unwrap()
                .into(),
        ]
    }

    #[test]
    fn a_helix_rises_one_pitch_per_turn_and_knows_its_length() {
        let helix = HelixCurve::new(Frame::WORLD, 3.0, 2.0, 2.0).unwrap();
        let tau = core::f64::consts::TAU;
        let start = helix.point_at(0.0, T).unwrap();
        let after_one_turn = helix.point_at(tau, T).unwrap();
        assert_relative_eq!(start.x, 3.0);
        assert_relative_eq!(after_one_turn.x, 3.0, epsilon = 1e-12);
        assert_relative_eq!(after_one_turn.y, 0.0, epsilon = 1e-12);
        assert_relative_eq!(after_one_turn.z - start.z, 2.0, epsilon = 1e-12);

        // Closed-form length: constant speed times swept angle — checked
        // against a fine chordal sum.
        let exact = helix.arc_length(0.0, 2.0 * tau);
        assert_relative_eq!(exact, 2.0 * tau * 3.0f64.hypot(2.0 / tau), epsilon = 1e-12);
        let mut chords = 0.0;
        let n = 20_000;
        for i in 0..n {
            let a = 2.0 * tau * f64::from(i) / f64::from(n);
            let b = 2.0 * tau * f64::from(i + 1) / f64::from(n);
            chords += helix
                .point_at(a, T)
                .unwrap()
                .distance(helix.point_at(b, T).unwrap());
        }
        assert!((exact - chords) / exact < 1e-6, "{exact} vs {chords}");

        // A negative pitch winds the other hand: same rise magnitude, the
        // quarter-turn point mirrored through the xz-plane... the y stays,
        // the z descends.
        let left = HelixCurve::new(Frame::WORLD, 3.0, -2.0, 2.0).unwrap();
        let q = left.point_at(tau / 4.0, T).unwrap();
        assert_relative_eq!(q.y, 3.0, epsilon = 1e-12);
        assert!(q.z < 0.0);
    }

    #[test]
    fn a_reversed_helix_swaps_its_ends_and_flips_its_tangent() {
        let helix: Curve = HelixCurve::new(tilted(), 2.0, 1.0, 1.5).unwrap().into();
        let (lo, hi) = helix.domain();
        let back = helix.reversed();
        assert_relative_eq!(
            helix
                .point_at(lo, T)
                .unwrap()
                .distance(back.point_at(hi, T).unwrap()),
            0.0,
            epsilon = 1e-12
        );
        let d_fwd = helix.d1_at(f64::midpoint(lo, hi), T).unwrap();
        let d_back = back.d1_at(f64::midpoint(lo, hi), T).unwrap();
        assert_relative_eq!((d_fwd + d_back).magnitude(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn a_helix_has_no_exact_spline_and_says_so() {
        let helix: Curve = HelixCurve::new(Frame::WORLD, 1.0, 1.0, 1.0).unwrap().into();
        assert!(helix.to_bspline(T).is_err());
    }

    /// Sample a curve evenly across its domain, avoiding the exact ends.
    fn interior(c: &Curve, n: usize) -> Vec<f64> {
        let (a, b) = c.domain();
        (1..n)
            .map(|i| {
                #[allow(clippy::cast_precision_loss)]
                let t = i as f64 / n as f64;
                a + (b - a) * t
            })
            .collect()
    }

    #[test]
    fn every_curves_derivative_agrees_with_finite_differences() {
        let h = 1e-6;
        for c in every_curve() {
            for u in interior(&c, 8) {
                let d1 = c.d1_at(u, T).unwrap();
                let numeric = (c.point_at(u + h, T).unwrap() - c.point_at(u - h, T).unwrap())
                    * (1.0 / (2.0 * h));
                let scale = numeric.magnitude().max(1.0);
                assert!(
                    (d1 - numeric).magnitude() <= 1e-5 * scale,
                    "{:?} at {u}: {d1:?} vs {numeric:?}",
                    c.kind()
                );
            }
        }
    }

    #[test]
    fn derivatives_at_zero_returns_the_point_itself() {
        for c in every_curve() {
            for u in interior(&c, 4) {
                let d = c.derivatives_at(u, 0, T).unwrap();
                assert_eq!(d.len(), 1);
                assert!(Point::from_vector(d[0]).is_equal(c.point_at(u, T).unwrap(), T));
            }
        }
    }

    #[test]
    fn out_of_domain_parameters_are_refused_for_non_periodic_curves() {
        for c in every_curve() {
            let (a, b) = c.domain();
            if c.is_periodic() {
                // A periodic curve accepts anything and wraps it.
                assert!(c.point_at(b + 1.0, T).is_ok());
                assert!(c.point_at(a - 1.0, T).is_ok());
            } else {
                assert!(c.point_at(b + 1.0, T).is_err(), "{:?}", c.kind());
                assert!(c.point_at(a - 1.0, T).is_err(), "{:?}", c.kind());
            }
        }
    }

    #[test]
    fn a_periodic_curve_wraps_to_the_same_point() {
        let c: Curve = CircleCurve::new(Circle::new(tilted(), 2.0, T).unwrap()).into();
        assert!(c.is_periodic() && c.is_closed(T));
        let base = c.point_at(0.7, T).unwrap();
        for k in [-2.0_f64, -1.0, 1.0, 3.0] {
            let wrapped = c
                .point_at(k.mul_add(core::f64::consts::TAU, 0.7), T)
                .unwrap();
            assert!(base.is_equal(wrapped, T), "wrap by {k} moved the point");
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
                assert!(
                    forward.is_equal(backward, T),
                    "{:?} at t = {t}: {forward:?} vs {backward:?}",
                    c.kind()
                );
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
    fn a_reversed_curves_tangent_points_the_other_way() {
        for c in every_curve() {
            let r = c.reversed();
            let (a, b) = c.domain();
            let u = a + (b - a) * 0.4;
            let forward = c.tangent_at(u, T).unwrap();
            let backward = r.tangent_at(mirror(u, a, b), T).unwrap();
            assert!(forward.is_opposite(backward, T), "{:?}", c.kind());
        }
    }

    #[test]
    fn transforms_move_curves_and_preserve_their_shape() {
        let t =
            Transform::rotation(Axis::X, 0.7) * Transform::translation(Vector::new(1.0, 2.0, 3.0));
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
    fn a_scaling_rescales_a_lines_arc_length_domain() {
        // A line's parameter is a length, so scaling must rescale the domain or
        // the segment silently changes extent.
        let line = LineCurve::segment(Point::ORIGIN, Point::new(3.0, 4.0, 0.0), T).unwrap();
        let c: Curve = line.into();
        assert_eq!(c.domain(), (0.0, 5.0));
        let scaled = c
            .transformed(&Transform::scaling(Point::ORIGIN, 2.0, T).unwrap(), T)
            .unwrap();
        assert_eq!(scaled.domain(), (0.0, 10.0));
        assert!(
            scaled
                .end(T)
                .unwrap()
                .is_equal(Point::new(6.0, 8.0, 0.0), T)
        );
    }

    #[test]
    fn mirroring_a_circle_moves_every_point_by_the_mirror() {
        // The frame carries its own axes through the transform, so evaluating
        // the mirrored circle at a parameter lands exactly where the mirror
        // sends the original point. No extra parameter flip is involved.
        let c: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.0, T).unwrap()).into();
        let m = Transform::plane_mirror(Point::ORIGIN, Direction::X);
        let mirrored = c.transformed(&m, T).unwrap();
        for u in interior(&c, 8) {
            let expected = m.apply(c.point_at(u, T).unwrap());
            assert!(
                mirrored.point_at(u, T).unwrap().is_equal(expected, T),
                "at {u}"
            );
        }
    }

    #[test]
    fn a_line_segments_parameter_is_arc_length() {
        let line = LineCurve::segment(Point::ORIGIN, Point::new(3.0, 4.0, 0.0), T).unwrap();
        assert_eq!(line.domain(), (0.0, 5.0));
        assert!(line.point_at(0.0, T).unwrap().is_equal(Point::ORIGIN, T));
        assert!(
            line.point_at(5.0, T)
                .unwrap()
                .is_equal(Point::new(3.0, 4.0, 0.0), T)
        );
        assert!(
            line.point_at(2.5, T)
                .unwrap()
                .is_equal(Point::new(1.5, 2.0, 0.0), T)
        );
        assert_relative_eq!(
            line.d1_at(1.0, T).unwrap().magnitude(),
            1.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn degenerate_constructions_are_refused() {
        assert!(LineCurve::segment(Point::ORIGIN, Point::ORIGIN, T).is_err());
        assert!(LineCurve::over(Axis::X, 1.0, 1.0).is_err());
        assert!(LineCurve::over(Axis::X, 0.0, f64::NAN).is_err());
        assert!(HyperbolaCurve::new(Hyperbola::new(tilted(), 1.0, 1.0, T).unwrap(), 0.0).is_err());
        assert!(ParabolaCurve::new(Parabola::new(tilted(), 1.0, T).unwrap(), -1.0).is_err());
    }

    #[test]
    fn trimming_is_bounds_checked() {
        let base: Curve = LineCurve::over(Axis::X, 0.0, 10.0).unwrap().into();
        assert!(TrimmedCurve::new(base.clone(), 2.0, 8.0, T).is_ok());
        assert!(
            TrimmedCurve::new(base.clone(), 8.0, 2.0, T).is_err(),
            "empty"
        );
        assert!(
            TrimmedCurve::new(base.clone(), 5.0, 5.0, T).is_err(),
            "empty"
        );
        assert!(
            TrimmedCurve::new(base, -1.0, 5.0, T).is_err(),
            "outside the basis"
        );
    }

    #[test]
    fn a_trimmed_curve_agrees_with_its_basis() {
        let base: Curve = CircleCurve::new(Circle::new(tilted(), 2.0, T).unwrap()).into();
        let trimmed = TrimmedCurve::new(base.clone(), 0.5, 2.0, T).unwrap();
        assert_eq!(trimmed.domain(), (0.5, 2.0));
        for i in 0..=8 {
            let u = 0.5 + 1.5 * (f64::from(i) / 8.0);
            assert!(
                trimmed
                    .point_at(u, T)
                    .unwrap()
                    .is_equal(base.point_at(u, T).unwrap(), T)
            );
        }
        assert!(trimmed.point_at(0.4, T).is_err());
        assert!(trimmed.point_at(2.1, T).is_err());
    }

    #[test]
    fn a_rational_curve_is_recognized_and_a_uniformly_weighted_one_is_not() {
        let knots = KnotVector::clamped_uniform(2, 3).unwrap();
        let points = [
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];

        let uniform: Vec<_> = points
            .iter()
            .map(|p| Weighted::new(*p, 3.0, T).unwrap())
            .collect();
        assert!(
            !BSplineCurve::rational(knots.clone(), uniform)
                .unwrap()
                .is_rational(),
            "equal weights are polynomial whatever their value"
        );

        let w = core::f64::consts::FRAC_1_SQRT_2;
        let arc: Vec<_> = points
            .iter()
            .zip([1.0, w, 1.0])
            .map(|(p, w)| Weighted::new(*p, w, T).unwrap())
            .collect();
        let c = BSplineCurve::rational(knots, arc).unwrap();
        assert!(c.is_rational());
        // And it is an exact circular arc, which no polynomial curve is.
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
    fn spline_continuity_follows_interior_knot_multiplicity() {
        let control = vec![
            Point::ORIGIN,
            Point::new(1.0, 1.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
            Point::new(3.0, 1.0, 0.0),
            Point::new(4.0, 0.0, 0.0),
        ];
        // Degree 3 with simple interior knots is C2 there — not C-infinity,
        // which would be a false claim about a piecewise polynomial.
        let smooth = BSplineCurve::new(
            KnotVector::clamped_uniform(3, control.len()).unwrap(),
            control.clone(),
            T,
        )
        .unwrap();
        assert_eq!(smooth.continuity(), Continuity::C2);

        // A single Bezier piece has no interior knots and is smooth to every
        // order.
        let bezier = BSplineCurve::new(
            KnotVector::clamped_uniform(4, control.len()).unwrap(),
            control.clone(),
            T,
        )
        .unwrap();
        assert_eq!(bezier.continuity(), Continuity::CInfinity);

        // Degree 3, interior multiplicity 2: C1.
        let kinked = BSplineCurve::new(
            KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0], 3).unwrap(),
            {
                let mut c = control.clone();
                c.push(Point::new(5.0, 1.0, 0.0));
                c
            },
            T,
        )
        .unwrap();
        assert_eq!(kinked.continuity(), Continuity::C1);

        // An interior knot at full multiplicity is a corner.
        let corner = BSplineCurve::new(
            KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0], 2).unwrap(),
            control,
            T,
        )
        .unwrap();
        assert_eq!(corner.continuity(), Continuity::C0);
    }

    #[test]
    fn spline_refinement_does_not_move_the_curve() {
        let control = vec![
            Point::ORIGIN,
            Point::new(1.0, 2.0, 0.0),
            Point::new(3.0, 1.0, 1.0),
            Point::new(5.0, 0.0, 2.0),
        ];
        let c = BSplineCurve::new(
            KnotVector::clamped_uniform(3, control.len()).unwrap(),
            control,
            T,
        )
        .unwrap();
        let refined = c.with_knot_inserted(0.5, 1, T).unwrap();
        let elevated = c.elevated(T).unwrap();
        assert_eq!(elevated.degree(), 4);
        for i in 0..=20 {
            let u = f64::from(i) / 20.0;
            let base = c.point_at(u, T).unwrap();
            assert!(refined.point_at(u, T).unwrap().is_equal(base, T));
            assert!(elevated.point_at(u, T).unwrap().is_equal(base, T));
        }
    }

    #[test]
    fn curvature_of_a_circle_is_the_reciprocal_of_its_radius() {
        for r in [0.5_f64, 2.0, 50.0] {
            let c: Curve = CircleCurve::new(Circle::new(tilted(), r, T).unwrap()).into();
            assert_relative_eq!(
                c.curvature_at(1.1, T).unwrap(),
                1.0 / r,
                max_relative = 1e-12
            );
        }
        // A line has no curvature.
        let line: Curve = LineCurve::segment(Point::ORIGIN, Point::new(1.0, 1.0, 1.0), T)
            .unwrap()
            .into();
        assert_relative_eq!(line.curvature_at(0.5, T).unwrap(), 0.0);
    }

    #[test]
    fn kinds_are_reported_for_dispatch() {
        let kinds: Vec<_> = every_curve().iter().map(Curve3d::kind).collect();
        assert_eq!(
            kinds,
            vec![
                CurveKind::Line,
                CurveKind::Circle,
                CurveKind::Ellipse,
                CurveKind::Hyperbola,
                CurveKind::Parabola,
                CurveKind::BSpline,
                CurveKind::Helix,
                CurveKind::Trimmed,
            ]
        );
    }
}
