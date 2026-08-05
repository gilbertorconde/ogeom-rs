//! The adaptor traits: what every algorithm sees geometry through.
//!
//! Intersection, projection, extrema and tessellation are written against
//! [`Curve3d`], [`Curve2d`] and [`Surface`], never against a concrete type.
//! Nothing downstream needs to know whether it holds an analytic cylinder or a
//! NURBS patch, which is what keeps one implementation of each algorithm rather
//! than one per surface type.
//!
//! # Fast paths, opted into
//!
//! [`Curve3d::kind`] and [`Surface::kind`] let an algorithm *ask*. Plane/plane
//! intersection is two lines of algebra and should never go through a marching
//! intersector; a caller that wants that shortcut matches on the kind and takes
//! it. The general path never has to know what it is looking at, so adding a
//! surface type does not break existing algorithms — it only forgoes a
//! shortcut until someone writes one.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_math::{Direction, Point, Point2, Transform, Vector, Vector2};

/// How smooth a curve or surface is.
///
/// Ordered from least to most smooth, so `>=` is a meaningful test.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Continuity {
    /// Positions agree; tangents may not. A corner.
    C0,
    /// Tangent *directions* agree, but their magnitudes need not. Enough for a
    /// visually smooth join, and the usual requirement for a wire.
    G1,
    /// First derivatives agree exactly.
    C1,
    /// Curvature agrees.
    G2,
    /// Second derivatives agree exactly.
    C2,
    /// Differentiable to any order — an analytic surface away from its
    /// degeneracies.
    CInfinity,
}

/// What kind of curve this is, for analytic fast paths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum CurveKind {
    /// A straight line.
    Line,
    /// A circle or circular arc.
    Circle,
    /// An ellipse or elliptical arc.
    Ellipse,
    /// One branch of a hyperbola.
    Hyperbola,
    /// A parabola.
    Parabola,
    /// A polynomial or rational Bézier.
    Bezier,
    /// A polynomial or rational B-spline.
    BSpline,
    /// A helix about an axis.
    Helix,
    /// A restriction of another curve to a sub-interval.
    Trimmed,
    /// A curve offset from another.
    Offset,
    /// A pcurve composed with the surface it is drawn on.
    OnSurface,
}

/// What kind of surface this is, for analytic fast paths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SurfaceKind {
    /// A plane.
    Plane,
    /// A circular cylinder.
    Cylinder,
    /// A circular cone.
    Cone,
    /// A sphere.
    Sphere,
    /// A torus.
    Torus,
    /// A polynomial or rational Bézier patch.
    Bezier,
    /// A polynomial or rational B-spline patch.
    BSpline,
    /// A surface of revolution.
    Revolution,
    /// A surface swept by translating a curve.
    Extrusion,
    /// A restriction of another surface to a sub-rectangle.
    Trimmed,
    /// A surface offset from another.
    Offset,
}

impl SurfaceKind {
    /// Whether this surface is a quadric — a plane, cylinder, cone or sphere.
    ///
    /// Quadric pairs have closed-form intersections, which is worth a great
    /// deal: it is the difference between an exact conic and a marched
    /// approximation of one.
    #[must_use]
    pub const fn is_quadric(self) -> bool {
        matches!(
            self,
            Self::Plane | Self::Cylinder | Self::Cone | Self::Sphere
        )
    }

    /// Whether this surface has an exact analytic form, as opposed to being
    /// defined by control points.
    #[must_use]
    pub const fn is_analytic(self) -> bool {
        matches!(
            self,
            Self::Plane | Self::Cylinder | Self::Cone | Self::Sphere | Self::Torus
        )
    }
}

/// A parametric curve in space.
///
/// Implementors must guarantee:
///
/// - [`Curve3d::domain`] returns `(a, b)` with `a < b`, both finite;
/// - [`Curve3d::point_at`] and the derivative methods agree — `d1_at` is the
///   derivative of `point_at`, and so on;
/// - a periodic curve's period is exactly its domain width.
pub trait Curve3d {
    /// The parameter interval over which the curve is defined.
    fn domain(&self) -> (f64, f64);

    /// The point at `u`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if `u` is outside the
    /// domain and the curve is not periodic.
    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point>;

    /// The first derivative at `u`.
    ///
    /// # Errors
    ///
    /// As [`Curve3d::point_at`].
    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector>;

    /// Derivatives up to order `n`, with `result[0]` the point itself.
    ///
    /// # Errors
    ///
    /// As [`Curve3d::point_at`].
    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector>>;

    /// What kind of curve this is.
    fn kind(&self) -> CurveKind;

    /// How smooth the curve is across its domain.
    fn continuity(&self) -> Continuity;

    /// Whether the curve's ends meet.
    fn is_closed(&self, tol: Tolerances) -> bool;

    /// Whether the curve continues past its domain by repeating.
    ///
    /// Distinct from being closed: a full circle is both, a closed B-spline
    /// whose ends merely coincide is closed but not periodic, and evaluating
    /// the latter outside its domain is an error.
    fn is_periodic(&self) -> bool;

    /// The unit tangent at `u`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) at a cusp,
    /// where the derivative vanishes.
    fn tangent_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Direction> {
        Direction::new(self.d1_at(u, tol)?, tol)
    }

    /// The curvature at `u`.
    ///
    /// # Errors
    ///
    /// As [`Curve3d::point_at`].
    fn curvature_at(&self, u: f64, tol: Tolerances) -> OgeomResult<f64> {
        let d = self.derivatives_at(u, 2, tol)?;
        let speed = d[1].magnitude();
        if speed == 0.0 {
            return Ok(0.0);
        }
        Ok(d[1].cross(d[2]).magnitude() / (speed * speed * speed))
    }

    /// The start point.
    ///
    /// # Errors
    ///
    /// As [`Curve3d::point_at`].
    fn start(&self, tol: Tolerances) -> OgeomResult<Point> {
        self.point_at(self.domain().0, tol)
    }

    /// The end point.
    ///
    /// # Errors
    ///
    /// As [`Curve3d::point_at`].
    fn end(&self, tol: Tolerances) -> OgeomResult<Point> {
        self.point_at(self.domain().1, tol)
    }

    /// Bring `u` into the domain, wrapping if the curve is periodic.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if `u` is outside the
    /// domain of a non-periodic curve by more than `tol.parametric()`.
    fn normalize_parameter(&self, u: f64, tol: Tolerances) -> OgeomResult<f64> {
        let (a, b) = self.domain();
        if self.is_periodic() {
            let period = b - a;
            return Ok(a + (u - a).rem_euclid(period));
        }
        if !u.is_finite() || u < a - tol.parametric() || u > b + tol.parametric() {
            return Err(ogeom_core::ogeom_err!(
                Domain,
                "parameter {u} outside curve domain [{a}, {b}]"
            ));
        }
        Ok(u.clamp(a, b))
    }
}

/// A parametric curve in the plane.
///
/// The same contract as [`Curve3d`], one dimension down. Kept as a separate
/// trait rather than a generic parameter because pcurves — curves in a
/// surface's parameter space — are used differently enough from spatial curves
/// that conflating them invites mistakes.
pub trait Curve2d {
    /// The parameter interval over which the curve is defined.
    fn domain(&self) -> (f64, f64);

    /// The point at `u`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if `u` is outside the
    /// domain and the curve is not periodic.
    fn point_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Point2>;

    /// The first derivative at `u`.
    ///
    /// # Errors
    ///
    /// As [`Curve2d::point_at`].
    fn d1_at(&self, u: f64, tol: Tolerances) -> OgeomResult<Vector2>;

    /// Derivatives up to order `n`, with `result[0]` the point itself.
    ///
    /// # Errors
    ///
    /// As [`Curve2d::point_at`].
    fn derivatives_at(&self, u: f64, n: usize, tol: Tolerances) -> OgeomResult<Vec<Vector2>>;

    /// What kind of curve this is.
    fn kind(&self) -> CurveKind;

    /// Whether the curve's ends meet.
    fn is_closed(&self, tol: Tolerances) -> bool;

    /// Whether the curve continues past its domain by repeating.
    fn is_periodic(&self) -> bool;

    /// The start point.
    ///
    /// # Errors
    ///
    /// As [`Curve2d::point_at`].
    fn start(&self, tol: Tolerances) -> OgeomResult<Point2> {
        self.point_at(self.domain().0, tol)
    }

    /// The end point.
    ///
    /// # Errors
    ///
    /// As [`Curve2d::point_at`].
    fn end(&self, tol: Tolerances) -> OgeomResult<Point2> {
        self.point_at(self.domain().1, tol)
    }
}

/// A parametric surface.
///
/// Implementors must guarantee:
///
/// - both domain intervals are non-empty and finite;
/// - the derivative methods agree with [`Surface::point_at`];
/// - `du` and `dv` are the derivatives along the first and second parameters
///   respectively, in that order, for every surface.
pub trait Surface {
    /// The `u` and `v` parameter intervals.
    fn domain(&self) -> ((f64, f64), (f64, f64));

    /// The point at `(u, v)`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if the parameters lie
    /// outside the domain in a direction that is not periodic.
    fn point_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Point>;

    /// The two first derivatives at `(u, v)`.
    ///
    /// # Errors
    ///
    /// As [`Surface::point_at`].
    fn d1_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector)>;

    /// The three second derivatives at `(u, v)`: `d2u`, `duv`, `d2v`.
    ///
    /// # Errors
    ///
    /// As [`Surface::point_at`].
    fn d2_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector, Vector)>;

    /// What kind of surface this is.
    fn kind(&self) -> SurfaceKind;

    /// How smooth the surface is across its domain.
    fn continuity(&self) -> Continuity;

    /// Whether the surface closes on itself along `u`.
    fn is_closed_u(&self, tol: Tolerances) -> bool;

    /// Whether the surface closes on itself along `v`.
    fn is_closed_v(&self, tol: Tolerances) -> bool;

    /// Whether `u` repeats past the domain.
    fn is_periodic_u(&self) -> bool;

    /// Whether `v` repeats past the domain.
    fn is_periodic_v(&self) -> bool;

    /// The unit normal at `(u, v)`, following `du x dv`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) at a
    /// degeneracy — a pole or an apex — where the tangents determine no normal.
    fn normal_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Direction> {
        let (du, dv) = self.d1_at(u, v, tol)?;
        // Compared against the square of the larger tangent: at a pole one
        // tangent vanishes, and a test relative to the *product* would find the
        // cross product large by comparison and call the point healthy.
        let scale = du.magnitude().max(dv.magnitude());
        if du.cross(dv).magnitude() <= tol.angular() * scale * scale {
            return Err(ogeom_core::ogeom_err!(
                Construction,
                "surface is degenerate at ({u}, {v}); no normal is determined"
            ));
        }
        Direction::new(du.cross(dv), tol)
    }

    /// Whether the surface degenerates at `(u, v)`.
    ///
    /// # Errors
    ///
    /// As [`Surface::point_at`].
    fn is_degenerate_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<bool> {
        let (du, dv) = self.d1_at(u, v, tol)?;
        let scale = du.magnitude().max(dv.magnitude());
        Ok(du.cross(dv).magnitude() <= tol.angular() * scale * scale)
    }

    /// Bring `(u, v)` into the domain, wrapping in whichever directions are
    /// periodic.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if a parameter is outside
    /// a non-periodic direction's domain by more than `tol.parametric()`.
    fn normalize_parameters(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(f64, f64)> {
        let ((ua, ub), (va, vb)) = self.domain();
        let fix = |t: f64, a: f64, b: f64, periodic: bool, name: &str| -> OgeomResult<f64> {
            if periodic {
                return Ok(a + (t - a).rem_euclid(b - a));
            }
            if !t.is_finite() || t < a - tol.parametric() || t > b + tol.parametric() {
                return Err(ogeom_core::ogeom_err!(
                    Domain,
                    "{name} parameter {t} outside [{a}, {b}]"
                ));
            }
            Ok(t.clamp(a, b))
        };
        Ok((
            fix(u, ua, ub, self.is_periodic_u(), "u")?,
            fix(v, va, vb, self.is_periodic_v(), "v")?,
        ))
    }
}

/// Geometry that can be moved by a similarity.
///
/// Separate from the evaluation traits because a *view* of geometry — a
/// borrowed adaptor over someone else's data — can be evaluated but not moved.
pub trait Transformable: Sized {
    /// This geometry moved by `t`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// transformed geometry would be degenerate.
    fn transformed(&self, t: &Transform, tol: Tolerances) -> OgeomResult<Self>;
}

/// Reversible geometry: the same point set, traversed the other way.
pub trait Reversible: Sized {
    /// This geometry with its parameter direction reversed.
    ///
    /// The domain is preserved, so a curve reversed still runs over the same
    /// interval — only the direction of travel changes. Preserving the domain
    /// matters because trimming ranges elsewhere refer to it.
    #[must_use]
    fn reversed(&self) -> Self;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn continuity_is_ordered_from_least_to_most_smooth() {
        assert!(Continuity::C0 < Continuity::G1);
        assert!(Continuity::G1 < Continuity::C1);
        assert!(Continuity::C1 < Continuity::G2);
        assert!(Continuity::G2 < Continuity::C2);
        assert!(Continuity::C2 < Continuity::CInfinity);
        // The ordering is what makes a requirement expressible as a comparison.
        assert!(Continuity::CInfinity >= Continuity::C1);
    }

    #[test]
    fn quadrics_and_analytics_are_classified_correctly() {
        for k in [
            SurfaceKind::Plane,
            SurfaceKind::Cylinder,
            SurfaceKind::Cone,
            SurfaceKind::Sphere,
        ] {
            assert!(k.is_quadric(), "{k:?}");
            assert!(k.is_analytic(), "{k:?}");
        }
        // A torus is analytic but quartic, not quadric — a distinction that
        // matters, since quadric pairs have closed-form intersections and
        // torus pairs do not.
        assert!(!SurfaceKind::Torus.is_quadric());
        assert!(SurfaceKind::Torus.is_analytic());

        for k in [
            SurfaceKind::BSpline,
            SurfaceKind::Bezier,
            SurfaceKind::Revolution,
            SurfaceKind::Extrusion,
            SurfaceKind::Trimmed,
            SurfaceKind::Offset,
        ] {
            assert!(!k.is_quadric(), "{k:?}");
            assert!(!k.is_analytic(), "{k:?}");
        }
    }
}
