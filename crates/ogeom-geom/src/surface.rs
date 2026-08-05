//! Concrete surfaces.
//!
//! The five analytic surfaces plus a NURBS patch, a surface of revolution, a
//! surface of extrusion, and a trimmed restriction of any of them. All reachable
//! through [`SurfaceGeometry`], an enum for the same reasons curves are one — see
//! [`crate::curve`].
//!
//! # Keeping analytic surfaces analytic
//!
//! A cylinder could be written as a NURBS patch, and some kernels do exactly
//! that. Keeping it a cylinder is worth the extra types: intersection can take a
//! closed-form path, measurement can report a radius rather than a fit, files
//! stay small, and a fillet knows it is filleting a cylinder. The cost is that
//! every algorithm must handle a handful of cases — which the `kind` method
//! makes an explicit choice rather than a hidden one.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{
    Axis, Cone, ControlGrid, Cylinder, Direction, KnotVector, Plane, Point, Sphere, Torus,
    Transform, Vector, Weighted, bspline, elementary,
};

use crate::curve::Curve;
use crate::traits::{Continuity, Curve3d, Surface, SurfaceKind, Transformable};

/// How far an unbounded surface's default domain reaches.
///
/// A plane and a cylinder are unbounded, but every interface here works on a
/// finite rectangle. Far past any real model, and well short of where `f64`
/// spacing turns coarse.
pub const SURFACE_EXTENT: f64 = 1.0e9;

const TAU: f64 = core::f64::consts::TAU;

/// A surface.
#[derive(Debug, Clone, PartialEq)]
pub enum SurfaceGeometry {
    /// A plane.
    Plane(PlaneSurface),
    /// A circular cylinder.
    Cylinder(CylinderSurface),
    /// A circular cone.
    Cone(ConeSurface),
    /// A sphere.
    Sphere(SphereSurface),
    /// A torus.
    Torus(TorusSurface),
    /// A polynomial or rational B-spline patch.
    BSpline(BSplineSurface),
    /// A curve revolved about an axis.
    Revolution(Box<RevolutionSurface>),
    /// A curve swept along a direction.
    Extrusion(Box<ExtrusionSurface>),
    /// Another surface restricted to a sub-rectangle.
    Trimmed(Box<TrimmedSurface>),
    /// A surface at a constant signed distance along another's normal.
    Offset(Box<OffsetSurface>),
}

/// A plane, parameterized by its frame's `x` and `y` axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlaneSurface {
    plane: Plane,
    domain: ((f64, f64), (f64, f64)),
}

/// A cylinder, parameterized by `(angle, height)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CylinderSurface {
    cylinder: Cylinder,
    height: (f64, f64),
}

/// A cone, parameterized by `(angle, height)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConeSurface {
    cone: Cone,
    height: (f64, f64),
}

/// A sphere, parameterized by `(longitude, latitude)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SphereSurface {
    sphere: Sphere,
}

/// A torus, parameterized by `(angle about the axis, angle around the tube)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TorusSurface {
    torus: Torus,
}

/// A tensor-product B-spline patch, polynomial or rational.
#[derive(Debug, Clone, PartialEq)]
pub struct BSplineSurface {
    u_knots: KnotVector,
    v_knots: KnotVector,
    grid: ControlGrid<Weighted<Point>>,
    rational: bool,
}

/// A curve revolved about an axis.
///
/// `u` is the angle of revolution; `v` is the generating curve's own parameter.
#[derive(Debug, Clone, PartialEq)]
pub struct RevolutionSurface {
    curve: Curve,
    axis: Axis,
    angle: (f64, f64),
}

/// A curve swept along a direction.
///
/// `u` is the generating curve's parameter; `v` is the distance swept.
#[derive(Debug, Clone, PartialEq)]
pub struct ExtrusionSurface {
    curve: Curve,
    direction: Direction,
    extent: (f64, f64),
}

/// Another surface restricted to a sub-rectangle of its domain.
#[derive(Debug, Clone, PartialEq)]
pub struct TrimmedSurface {
    basis: SurfaceGeometry,
    domain: ((f64, f64), (f64, f64)),
}

/// A surface displaced a constant signed distance along its basis's normal,
/// sharing the basis's parameterization.
///
/// Point and first derivatives are exact — the normal's derivative is the
/// projection formula over the basis's second derivatives. The *second*
/// derivative would need the basis's third, which the vocabulary does not
/// carry, so `d2_at` refuses by name rather than differencing quietly. For
/// an analytic basis the offset is itself analytic and the direct type is
/// the better spelling; this type exists for the free-form bases that have
/// no such spelling.
#[derive(Debug, Clone, PartialEq)]
pub struct OffsetSurface {
    basis: SurfaceGeometry,
    distance: f64,
}

/// Reject an empty or non-finite parameter range.
fn check_range(name: &str, lo: f64, hi: f64) -> OgeomResult<()> {
    if !lo.is_finite() || !hi.is_finite() || hi <= lo {
        ogeom_bail!(
            Construction,
            "{name} range [{lo}, {hi}] is empty or non-finite"
        );
    }
    Ok(())
}

impl PlaneSurface {
    /// A plane spanning [`SURFACE_EXTENT`] in both directions.
    #[must_use]
    pub const fn new(plane: Plane) -> Self {
        Self {
            plane,
            domain: (
                (-SURFACE_EXTENT, SURFACE_EXTENT),
                (-SURFACE_EXTENT, SURFACE_EXTENT),
            ),
        }
    }

    /// A plane over an explicit rectangle of its own coordinates.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if either
    /// range is empty.
    pub fn over(plane: Plane, u: (f64, f64), v: (f64, f64)) -> OgeomResult<Self> {
        check_range("u", u.0, u.1)?;
        check_range("v", v.0, v.1)?;
        Ok(Self {
            plane,
            domain: (u, v),
        })
    }

    /// The underlying plane.
    #[must_use]
    pub const fn plane(&self) -> Plane {
        self.plane
    }
}

impl CylinderSurface {
    /// A cylinder over a height range.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the range
    /// is empty.
    pub fn new(cylinder: Cylinder, height: (f64, f64)) -> OgeomResult<Self> {
        check_range("height", height.0, height.1)?;
        Ok(Self { cylinder, height })
    }

    /// The underlying cylinder.
    #[must_use]
    pub const fn cylinder(&self) -> Cylinder {
        self.cylinder
    }
}

impl ConeSurface {
    /// A cone over a height range.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the range
    /// is empty.
    pub fn new(cone: Cone, height: (f64, f64)) -> OgeomResult<Self> {
        check_range("height", height.0, height.1)?;
        Ok(Self { cone, height })
    }

    /// The underlying cone.
    #[must_use]
    pub const fn cone(&self) -> Cone {
        self.cone
    }

    /// The height at which the radius vanishes — the apex.
    #[must_use]
    pub fn apex_height(&self) -> f64 {
        -self.cone.reference_radius() / self.cone.half_angle().tan()
    }
}

impl SphereSurface {
    /// A whole sphere.
    #[must_use]
    pub const fn new(sphere: Sphere) -> Self {
        Self { sphere }
    }

    /// The underlying sphere.
    #[must_use]
    pub const fn sphere(&self) -> Sphere {
        self.sphere
    }
}

impl TorusSurface {
    /// A whole torus.
    #[must_use]
    pub const fn new(torus: Torus) -> Self {
        Self { torus }
    }

    /// The underlying torus.
    #[must_use]
    pub const fn torus(&self) -> Torus {
        self.torus
    }
}

impl BSplineSurface {
    /// A polynomial patch.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) on a shape mismatch.
    pub fn new(
        u_knots: KnotVector,
        v_knots: KnotVector,
        grid: &ControlGrid<Point>,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        let weighted = ControlGrid::new(
            grid.points()
                .iter()
                .map(|p| Weighted::new(*p, 1.0, tol))
                .collect::<OgeomResult<Vec<_>>>()?,
            grid.u_count(),
            grid.v_count(),
        )?;
        Self::rational(u_knots, v_knots, weighted)
    }

    /// A rational patch.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) on a shape mismatch.
    pub fn rational(
        u_knots: KnotVector,
        v_knots: KnotVector,
        grid: ControlGrid<Weighted<Point>>,
    ) -> OgeomResult<Self> {
        if grid.u_count() != u_knots.control_point_count()
            || grid.v_count() != v_knots.control_point_count()
        {
            ogeom_bail!(
                Dimension,
                "knot vectors describe a {}x{} grid, got {}x{}",
                u_knots.control_point_count(),
                v_knots.control_point_count(),
                grid.u_count(),
                grid.v_count()
            );
        }
        let first = grid.points()[0].weight;
        let rational = grid
            .points()
            .iter()
            .any(|w| (w.weight - first).abs() > 1e-12 * first.abs());
        Ok(Self {
            u_knots,
            v_knots,
            grid,
            rational,
        })
    }

    /// The `u` knot vector.
    #[must_use]
    pub const fn u_knots(&self) -> &KnotVector {
        &self.u_knots
    }

    /// The `v` knot vector.
    #[must_use]
    pub const fn v_knots(&self) -> &KnotVector {
        &self.v_knots
    }

    /// The control grid.
    #[must_use]
    pub const fn grid(&self) -> &ControlGrid<Weighted<Point>> {
        &self.grid
    }

    /// Whether the weights differ, so the patch is genuinely rational.
    #[must_use]
    pub const fn is_rational(&self) -> bool {
        self.rational
    }
}

impl RevolutionSurface {
    /// Revolve `curve` about `axis` through `angle` radians.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the angle
    /// is not finite and positive, or exceeds a full turn.
    pub fn new(curve: Curve, axis: Axis, angle: f64) -> OgeomResult<Self> {
        if !angle.is_finite() || angle <= 0.0 || angle > TAU + 1e-12 {
            ogeom_bail!(
                Construction,
                "revolution angle {angle} must be in (0, 2*pi]"
            );
        }
        Ok(Self {
            curve,
            axis,
            angle: (0.0, angle.min(TAU)),
        })
    }

    /// The generating curve.
    #[must_use]
    pub const fn curve(&self) -> &Curve {
        &self.curve
    }

    /// The axis of revolution.
    #[must_use]
    pub const fn axis(&self) -> Axis {
        self.axis
    }
}

impl ExtrusionSurface {
    /// Sweep `curve` along `direction` over `[0, distance]`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `distance`
    /// is not finite and positive.
    pub fn new(curve: Curve, direction: Direction, distance: f64) -> OgeomResult<Self> {
        if !distance.is_finite() || distance <= 0.0 {
            ogeom_bail!(
                Construction,
                "extrusion distance {distance} must be positive"
            );
        }
        Ok(Self {
            curve,
            direction,
            extent: (0.0, distance),
        })
    }

    /// The generating curve.
    #[must_use]
    pub const fn curve(&self) -> &Curve {
        &self.curve
    }

    /// The sweep direction.
    #[must_use]
    pub const fn direction(&self) -> Direction {
        self.direction
    }
}

impl TrimmedSurface {
    /// Restrict `basis` to a sub-rectangle.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) if either range is empty
    /// or leaves the basis surface's domain in a non-periodic direction.
    pub fn new(
        basis: SurfaceGeometry,
        u: (f64, f64),
        v: (f64, f64),
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        check_range("u", u.0, u.1)?;
        check_range("v", v.0, v.1)?;
        let ((ua, ub), (va, vb)) = basis.domain();
        let eps = tol.parametric();
        if !basis.is_periodic_u() && (u.0 < ua - eps || u.1 > ub + eps) {
            ogeom_bail!(Domain, "u range [{}, {}] leaves [{ua}, {ub}]", u.0, u.1);
        }
        if !basis.is_periodic_v() && (v.0 < va - eps || v.1 > vb + eps) {
            ogeom_bail!(Domain, "v range [{}, {}] leaves [{va}, {vb}]", v.0, v.1);
        }
        Ok(Self {
            basis,
            domain: (u, v),
        })
    }

    /// The surface being trimmed.
    #[must_use]
    pub const fn basis(&self) -> &SurfaceGeometry {
        &self.basis
    }
}

impl OffsetSurface {
    /// Offset `basis` by a signed `distance` along its own normal.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// distance is not finite and non-zero.
    pub fn new(basis: SurfaceGeometry, distance: f64) -> OgeomResult<Self> {
        if !distance.is_finite() || distance == 0.0 {
            ogeom_bail!(
                Construction,
                "an offset of {distance} is not a displacement"
            );
        }
        Ok(Self { basis, distance })
    }

    /// The surface being offset.
    #[must_use]
    pub const fn basis(&self) -> &SurfaceGeometry {
        &self.basis
    }

    /// The signed displacement along the basis normal.
    #[must_use]
    pub const fn distance(&self) -> f64 {
        self.distance
    }
}

impl Surface for OffsetSurface {
    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        self.basis.domain()
    }

    fn point_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Point> {
        let base = self.basis.point_at(u, v, tol)?;
        let normal = self.basis.normal_at(u, v, tol)?;
        Ok(base + normal.vector() * self.distance)
    }

    fn d1_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector)> {
        // Differentiate S + d·c/|c| with c = Su x Sv: the unit normal's
        // derivative is the tangential projection of c's derivative, scaled
        // by the magnitude — exact from the basis's first and second
        // derivatives.
        let (su, sv) = self.basis.d1_at(u, v, tol)?;
        let (suu, suv, svv) = self.basis.d2_at(u, v, tol)?;
        let c = su.cross(sv);
        let m = c.magnitude();
        let scale = su.magnitude().max(sv.magnitude());
        if m <= tol.angular() * scale * scale {
            ogeom_bail!(
                Construction,
                "the basis is degenerate at ({u}, {v}); the offset has no tangent plane there"
            );
        }
        let n = c / m;
        let cu = suu.cross(sv) + su.cross(suv);
        let cv = suv.cross(sv) + su.cross(svv);
        let nu = (cu - n * n.dot(cu)) / m;
        let nv = (cv - n * n.dot(cv)) / m;
        Ok((su + nu * self.distance, sv + nv * self.distance))
    }

    fn d2_at(&self, _u: f64, _v: f64, _tol: Tolerances) -> OgeomResult<(Vector, Vector, Vector)> {
        ogeom_bail!(
            Construction,
            "an offset surface's second derivative needs its basis's third, which the              vocabulary does not carry; offset the basis analytically or fit at a stated              tolerance instead"
        )
    }

    fn kind(&self) -> SurfaceKind {
        SurfaceKind::Offset
    }

    fn continuity(&self) -> Continuity {
        // Offsetting spends one order of smoothness.
        match self.basis.continuity() {
            Continuity::CInfinity => Continuity::CInfinity,
            Continuity::C2 | Continuity::G2 => Continuity::C1,
            Continuity::C1 | Continuity::G1 | Continuity::C0 => Continuity::C0,
        }
    }

    fn is_closed_u(&self, tol: Tolerances) -> bool {
        self.basis.is_closed_u(tol)
    }

    fn is_closed_v(&self, tol: Tolerances) -> bool {
        self.basis.is_closed_v(tol)
    }

    fn is_periodic_u(&self) -> bool {
        self.basis.is_periodic_u()
    }

    fn is_periodic_v(&self) -> bool {
        self.basis.is_periodic_v()
    }
}

impl Surface for PlaneSurface {
    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        self.domain
    }

    fn point_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Point> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        Ok(elementary::plane_at(&self.plane, u, v).point)
    }

    fn d1_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector)> {
        self.normalize_parameters(u, v, tol)?;
        let f = self.plane.frame();
        Ok((f.x().vector(), f.y().vector()))
    }

    fn d2_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector, Vector)> {
        self.normalize_parameters(u, v, tol)?;
        Ok((Vector::ZERO, Vector::ZERO, Vector::ZERO))
    }

    fn kind(&self) -> SurfaceKind {
        SurfaceKind::Plane
    }

    fn continuity(&self) -> Continuity {
        Continuity::CInfinity
    }

    fn is_closed_u(&self, _tol: Tolerances) -> bool {
        false
    }

    fn is_closed_v(&self, _tol: Tolerances) -> bool {
        false
    }

    fn is_periodic_u(&self) -> bool {
        false
    }

    fn is_periodic_v(&self) -> bool {
        false
    }
}

impl Surface for CylinderSurface {
    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        ((0.0, TAU), self.height)
    }

    fn point_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Point> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        Ok(elementary::cylinder_at(&self.cylinder, u, v).point)
    }

    fn d1_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector)> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        let p = elementary::cylinder_at(&self.cylinder, u, v);
        Ok((p.du, p.dv))
    }

    fn d2_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector, Vector)> {
        let (u, _) = self.normalize_parameters(u, v, tol)?;
        let f = self.cylinder.frame();
        let r = self.cylinder.radius();
        let (sin, cos) = u.sin_cos();
        // Second derivative in u points back at the axis; the surface is ruled
        // along v, so everything involving v vanishes.
        Ok((
            f.x() * (-r * cos) + f.y() * (-r * sin),
            Vector::ZERO,
            Vector::ZERO,
        ))
    }

    fn kind(&self) -> SurfaceKind {
        SurfaceKind::Cylinder
    }

    fn continuity(&self) -> Continuity {
        Continuity::CInfinity
    }

    fn is_closed_u(&self, _tol: Tolerances) -> bool {
        true
    }

    fn is_closed_v(&self, _tol: Tolerances) -> bool {
        false
    }

    fn is_periodic_u(&self) -> bool {
        true
    }

    fn is_periodic_v(&self) -> bool {
        false
    }
}

impl Surface for ConeSurface {
    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        ((0.0, TAU), self.height)
    }

    fn point_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Point> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        Ok(elementary::cone_at(&self.cone, u, v).point)
    }

    fn d1_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector)> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        let p = elementary::cone_at(&self.cone, u, v);
        Ok((p.du, p.dv))
    }

    fn d2_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector, Vector)> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        let f = self.cone.frame();
        let r = self.cone.radius_at(v);
        let slope = self.cone.half_angle().tan();
        let (sin, cos) = u.sin_cos();
        Ok((
            // Radial, scaled by the radius at this height.
            f.x() * (-r * cos) + f.y() * (-r * sin),
            // The radius grows linearly along v, so the mixed partial is the
            // u-tangent's rate of growth.
            f.x() * (-slope * sin) + f.y() * (slope * cos),
            // Straight along the ruling.
            Vector::ZERO,
        ))
    }

    fn kind(&self) -> SurfaceKind {
        SurfaceKind::Cone
    }

    fn continuity(&self) -> Continuity {
        Continuity::CInfinity
    }

    fn is_closed_u(&self, _tol: Tolerances) -> bool {
        true
    }

    fn is_closed_v(&self, _tol: Tolerances) -> bool {
        false
    }

    fn is_periodic_u(&self) -> bool {
        true
    }

    fn is_periodic_v(&self) -> bool {
        false
    }
}

impl Surface for SphereSurface {
    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        let half = core::f64::consts::FRAC_PI_2;
        ((0.0, TAU), (-half, half))
    }

    fn point_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Point> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        Ok(elementary::sphere_at(&self.sphere, u, v).point)
    }

    fn d1_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector)> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        let p = elementary::sphere_at(&self.sphere, u, v);
        Ok((p.du, p.dv))
    }

    fn d2_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector, Vector)> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        let f = self.sphere.frame();
        let r = self.sphere.radius();
        let (sin_u, cos_u) = u.sin_cos();
        let (sin_v, cos_v) = v.sin_cos();
        let (x, y, z) = (f.x().vector(), f.y().vector(), f.z().vector());
        let ring = r * cos_v;
        Ok((
            x * (-ring * cos_u) + y * (-ring * sin_u),
            x * (r * sin_v * sin_u) + y * (-r * sin_v * cos_u),
            x * (-ring * cos_u) + y * (-ring * sin_u) + z * (-r * sin_v),
        ))
    }

    fn kind(&self) -> SurfaceKind {
        SurfaceKind::Sphere
    }

    fn continuity(&self) -> Continuity {
        Continuity::CInfinity
    }

    fn is_closed_u(&self, _tol: Tolerances) -> bool {
        true
    }

    fn is_closed_v(&self, _tol: Tolerances) -> bool {
        false
    }

    fn is_periodic_u(&self) -> bool {
        true
    }

    fn is_periodic_v(&self) -> bool {
        // Latitude runs pole to pole and stops. Treating it as periodic would
        // let a parameter past the pole wrap round to the far side, which is a
        // different point.
        false
    }
}

impl Surface for TorusSurface {
    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        ((0.0, TAU), (0.0, TAU))
    }

    fn point_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Point> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        Ok(elementary::torus_at(&self.torus, u, v).point)
    }

    fn d1_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector)> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        let p = elementary::torus_at(&self.torus, u, v);
        Ok((p.du, p.dv))
    }

    fn d2_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector, Vector)> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        let f = self.torus.frame();
        let (major, minor) = (self.torus.major_radius(), self.torus.minor_radius());
        let (sin_u, cos_u) = u.sin_cos();
        let (sin_v, cos_v) = v.sin_cos();
        let (x, y, z) = (f.x().vector(), f.y().vector(), f.z().vector());
        let out = x * cos_u + y * sin_u;
        let side = x * -sin_u + y * cos_u;
        let radius = minor.mul_add(cos_v, major);
        Ok((
            out * -radius,
            side * (-minor * sin_v),
            out * (-minor * cos_v) + z * (-minor * sin_v),
        ))
    }

    fn kind(&self) -> SurfaceKind {
        SurfaceKind::Torus
    }

    fn continuity(&self) -> Continuity {
        Continuity::CInfinity
    }

    fn is_closed_u(&self, _tol: Tolerances) -> bool {
        true
    }

    fn is_closed_v(&self, _tol: Tolerances) -> bool {
        true
    }

    fn is_periodic_u(&self) -> bool {
        true
    }

    fn is_periodic_v(&self) -> bool {
        true
    }
}

impl Surface for BSplineSurface {
    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        (self.u_knots.domain(), self.v_knots.domain())
    }

    fn point_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Point> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        if self.rational {
            bspline::evaluate_rational_surface(&self.u_knots, &self.v_knots, &self.grid, u, v, tol)
        } else {
            Ok(
                bspline::evaluate_surface(&self.u_knots, &self.v_knots, &self.grid, u, v, tol)?
                    .point(),
            )
        }
    }

    fn d1_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector)> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        let d = bspline::rational_surface_derivatives(
            &self.u_knots,
            &self.v_knots,
            &self.grid,
            u,
            v,
            1,
            tol,
        )?;
        Ok((d[1][0].to_vector(), d[0][1].to_vector()))
    }

    fn d2_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector, Vector)> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        let d = bspline::rational_surface_derivatives(
            &self.u_knots,
            &self.v_knots,
            &self.grid,
            u,
            v,
            2,
            tol,
        )?;
        Ok((
            d[2][0].to_vector(),
            d[1][1].to_vector(),
            d[0][2].to_vector(),
        ))
    }

    fn kind(&self) -> SurfaceKind {
        SurfaceKind::BSpline
    }

    fn continuity(&self) -> Continuity {
        // The worse of the two directions governs the patch.
        let worst = |k: &KnotVector| {
            let (a, b) = k.domain();
            k.distinct()
                .into_iter()
                .filter(|(x, _)| *x > a && *x < b)
                .map(|(_, m)| m)
                .max()
                .map_or(Continuity::CInfinity, |m| {
                    match k.degree().saturating_sub(m) {
                        0 => Continuity::C0,
                        1 => Continuity::C1,
                        _ => Continuity::C2,
                    }
                })
        };
        worst(&self.u_knots).min(worst(&self.v_knots))
    }

    fn is_closed_u(&self, tol: Tolerances) -> bool {
        let last = self.grid.u_count() - 1;
        (0..self.grid.v_count()).all(|j| match (self.grid.get(0, j), self.grid.get(last, j)) {
            (Some(a), Some(b)) => a.point().is_equal(b.point(), tol),
            _ => false,
        })
    }

    fn is_closed_v(&self, tol: Tolerances) -> bool {
        let last = self.grid.v_count() - 1;
        (0..self.grid.u_count()).all(|i| match (self.grid.get(i, 0), self.grid.get(i, last)) {
            (Some(a), Some(b)) => a.point().is_equal(b.point(), tol),
            _ => false,
        })
    }

    fn is_periodic_u(&self) -> bool {
        false
    }

    fn is_periodic_v(&self) -> bool {
        false
    }
}

impl Surface for RevolutionSurface {
    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        (self.angle, self.curve.domain())
    }

    fn point_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Point> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        let p = self.curve.point_at(v, tol)?;
        Ok(Transform::rotation(self.axis, u).apply(p))
    }

    fn d1_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector)> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        let rotate = Transform::rotation(self.axis, u);
        let p = rotate.apply(self.curve.point_at(v, tol)?);
        // Rotating about an axis moves a point along a circle centred on the
        // axis, so the u-tangent is the axis direction crossed with the radius.
        let radius = p - self.axis.project(p);
        Ok((
            self.axis.direction.cross_with(radius),
            rotate.apply_vector(self.curve.d1_at(v, tol)?),
        ))
    }

    fn d2_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector, Vector)> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        let rotate = Transform::rotation(self.axis, u);
        let p = rotate.apply(self.curve.point_at(v, tol)?);
        let radius = p - self.axis.project(p);
        let d = self.curve.derivatives_at(v, 2, tol)?;
        let curve_d1 = rotate.apply_vector(d[1]);
        // Circular motion: the second derivative in u points back at the axis.
        let d2u = -radius;
        // The mixed partial is the u-derivative of the rotated curve tangent,
        // which is the same circular relation applied to that tangent.
        let duv = self.axis.direction.cross_with(curve_d1);
        Ok((d2u, duv, rotate.apply_vector(d[2])))
    }

    fn kind(&self) -> SurfaceKind {
        SurfaceKind::Revolution
    }

    fn continuity(&self) -> Continuity {
        // Rotation is smooth, so the generating curve governs.
        self.curve.continuity()
    }

    fn is_closed_u(&self, _tol: Tolerances) -> bool {
        (self.angle.1 - self.angle.0 - TAU).abs() <= 1e-12
    }

    fn is_closed_v(&self, tol: Tolerances) -> bool {
        self.curve.is_closed(tol)
    }

    fn is_periodic_u(&self) -> bool {
        (self.angle.1 - self.angle.0 - TAU).abs() <= 1e-12
    }

    fn is_periodic_v(&self) -> bool {
        self.curve.is_periodic()
    }
}

impl Surface for ExtrusionSurface {
    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        (self.curve.domain(), self.extent)
    }

    fn point_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Point> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        Ok(self.curve.point_at(u, tol)? + self.direction * v)
    }

    fn d1_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector)> {
        let (u, _) = self.normalize_parameters(u, v, tol)?;
        Ok((self.curve.d1_at(u, tol)?, self.direction.vector()))
    }

    fn d2_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector, Vector)> {
        let (u, _) = self.normalize_parameters(u, v, tol)?;
        // Ruled along v, so only the curve's own curvature survives.
        Ok((
            self.curve.derivatives_at(u, 2, tol)?[2],
            Vector::ZERO,
            Vector::ZERO,
        ))
    }

    fn kind(&self) -> SurfaceKind {
        SurfaceKind::Extrusion
    }

    fn continuity(&self) -> Continuity {
        self.curve.continuity()
    }

    fn is_closed_u(&self, tol: Tolerances) -> bool {
        self.curve.is_closed(tol)
    }

    fn is_closed_v(&self, _tol: Tolerances) -> bool {
        false
    }

    fn is_periodic_u(&self) -> bool {
        self.curve.is_periodic()
    }

    fn is_periodic_v(&self) -> bool {
        false
    }
}

impl Surface for TrimmedSurface {
    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        self.domain
    }

    fn point_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Point> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        self.basis.point_at(u, v, tol)
    }

    fn d1_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector)> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        self.basis.d1_at(u, v, tol)
    }

    fn d2_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector, Vector)> {
        let (u, v) = self.normalize_parameters(u, v, tol)?;
        self.basis.d2_at(u, v, tol)
    }

    fn kind(&self) -> SurfaceKind {
        SurfaceKind::Trimmed
    }

    fn continuity(&self) -> Continuity {
        self.basis.continuity()
    }

    fn is_closed_u(&self, _tol: Tolerances) -> bool {
        // A trim narrower than the basis cannot close, whatever the basis does.
        false
    }

    fn is_closed_v(&self, _tol: Tolerances) -> bool {
        false
    }

    fn is_periodic_u(&self) -> bool {
        false
    }

    fn is_periodic_v(&self) -> bool {
        false
    }
}

/// Dispatch a method across every surface variant.
macro_rules! dispatch {
    ($self:ident, $s:ident => $body:expr) => {
        match $self {
            Self::Plane($s) => $body,
            Self::Cylinder($s) => $body,
            Self::Cone($s) => $body,
            Self::Sphere($s) => $body,
            Self::Torus($s) => $body,
            Self::BSpline($s) => $body,
            Self::Revolution($s) => $body,
            Self::Extrusion($s) => $body,
            Self::Trimmed($s) => $body,
            Self::Offset($s) => $body,
        }
    };
}

impl Surface for SurfaceGeometry {
    fn domain(&self) -> ((f64, f64), (f64, f64)) {
        dispatch!(self, s => s.domain())
    }

    fn point_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<Point> {
        dispatch!(self, s => s.point_at(u, v, tol))
    }

    fn d1_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector)> {
        dispatch!(self, s => s.d1_at(u, v, tol))
    }

    fn d2_at(&self, u: f64, v: f64, tol: Tolerances) -> OgeomResult<(Vector, Vector, Vector)> {
        dispatch!(self, s => s.d2_at(u, v, tol))
    }

    fn kind(&self) -> SurfaceKind {
        dispatch!(self, s => s.kind())
    }

    fn continuity(&self) -> Continuity {
        dispatch!(self, s => s.continuity())
    }

    fn is_closed_u(&self, tol: Tolerances) -> bool {
        dispatch!(self, s => s.is_closed_u(tol))
    }

    fn is_closed_v(&self, tol: Tolerances) -> bool {
        dispatch!(self, s => s.is_closed_v(tol))
    }

    fn is_periodic_u(&self) -> bool {
        dispatch!(self, s => s.is_periodic_u())
    }

    fn is_periodic_v(&self) -> bool {
        dispatch!(self, s => s.is_periodic_v())
    }
}

impl Transformable for SurfaceGeometry {
    fn transformed(&self, t: &Transform, tol: Tolerances) -> OgeomResult<Self> {
        let scale = t.scale_factor().abs();
        Ok(match self {
            Self::Plane(s) => Self::Plane(PlaneSurface {
                plane: s.plane.transformed(t, tol)?,
                // The parameters are lengths along the frame axes, so a scaling
                // rescales the extent with them.
                domain: (
                    (s.domain.0.0 * scale, s.domain.0.1 * scale),
                    (s.domain.1.0 * scale, s.domain.1.1 * scale),
                ),
            }),
            Self::Cylinder(s) => Self::Cylinder(CylinderSurface {
                cylinder: s.cylinder.transformed(t, tol)?,
                height: (s.height.0 * scale, s.height.1 * scale),
            }),
            Self::Offset(s) => Self::Offset(Box::new(OffsetSurface {
                basis: s.basis.transformed(t, tol)?,
                distance: s.distance * scale,
            })),
            Self::Cone(s) => Self::Cone(ConeSurface {
                cone: s.cone.transformed(t, tol)?,
                height: (s.height.0 * scale, s.height.1 * scale),
            }),
            Self::Sphere(s) => Self::Sphere(SphereSurface {
                sphere: s.sphere.transformed(t, tol)?,
            }),
            Self::Torus(s) => Self::Torus(TorusSurface {
                torus: s.torus.transformed(t, tol)?,
            }),
            Self::BSpline(s) => {
                let grid = ControlGrid::new(
                    s.grid
                        .points()
                        .iter()
                        .map(|w| Weighted::new(t.apply(w.point()), w.weight, tol))
                        .collect::<OgeomResult<Vec<_>>>()?,
                    s.grid.u_count(),
                    s.grid.v_count(),
                )?;
                Self::BSpline(BSplineSurface { grid, ..s.clone() })
            }
            Self::Revolution(s) => Self::Revolution(Box::new(RevolutionSurface {
                curve: s.curve.transformed(t, tol)?,
                axis: Axis::new(
                    t.apply(s.axis.location),
                    t.apply_direction(s.axis.direction, tol)?,
                ),
                angle: s.angle,
            })),
            Self::Extrusion(s) => Self::Extrusion(Box::new(ExtrusionSurface {
                curve: s.curve.transformed(t, tol)?,
                direction: t.apply_direction(s.direction, tol)?,
                extent: (s.extent.0 * scale, s.extent.1 * scale),
            })),
            Self::Trimmed(s) => {
                let basis = s.basis.transformed(t, tol)?;
                // The trim range lives in the basis surface's parameters, and
                // those rescale exactly when the basis domain does.
                let ((oa, _), (ob, _)) = s.basis.domain();
                let ((na, _), (nb, _)) = basis.domain();
                let ur = if oa == 0.0 { 1.0 } else { na / oa };
                let vr = if ob == 0.0 { 1.0 } else { nb / ob };
                Self::Trimmed(Box::new(TrimmedSurface {
                    basis,
                    domain: (
                        (s.domain.0.0 * ur, s.domain.0.1 * ur),
                        (s.domain.1.0 * vr, s.domain.1.1 * vr),
                    ),
                }))
            }
        })
    }
}

impl From<PlaneSurface> for SurfaceGeometry {
    fn from(s: PlaneSurface) -> Self {
        Self::Plane(s)
    }
}
impl From<CylinderSurface> for SurfaceGeometry {
    fn from(s: CylinderSurface) -> Self {
        Self::Cylinder(s)
    }
}
impl From<ConeSurface> for SurfaceGeometry {
    fn from(s: ConeSurface) -> Self {
        Self::Cone(s)
    }
}
impl From<SphereSurface> for SurfaceGeometry {
    fn from(s: SphereSurface) -> Self {
        Self::Sphere(s)
    }
}
impl From<TorusSurface> for SurfaceGeometry {
    fn from(s: TorusSurface) -> Self {
        Self::Torus(s)
    }
}
impl From<BSplineSurface> for SurfaceGeometry {
    fn from(s: BSplineSurface) -> Self {
        Self::BSpline(s)
    }
}
impl From<RevolutionSurface> for SurfaceGeometry {
    fn from(s: RevolutionSurface) -> Self {
        Self::Revolution(Box::new(s))
    }
}
impl From<ExtrusionSurface> for SurfaceGeometry {
    fn from(s: ExtrusionSurface) -> Self {
        Self::Extrusion(Box::new(s))
    }
}
impl From<TrimmedSurface> for SurfaceGeometry {
    fn from(s: TrimmedSurface) -> Self {
        Self::Trimmed(Box::new(s))
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::curve::{CircleCurve, LineCurve};
    use approx::assert_relative_eq;
    use ogeom_math::{Circle, Frame};

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn an_offset_cylinder_is_the_larger_cylinder() {
        use ogeom_math::Cylinder;
        let basis = SurfaceGeometry::Cylinder(
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 2.0, T).unwrap(), (0.0, 5.0)).unwrap(),
        );
        let offset = OffsetSurface::new(basis, 1.5).unwrap();
        // The cylinder's outward normal is radial, so every offset point
        // stands at the grown radius, at the same height.
        for (u, v) in [(0.0, 0.0), (1.0, 2.0), (3.5, 4.5)] {
            let p = offset.point_at(u, v, T).unwrap();
            let radial = (p.x * p.x + p.y * p.y).sqrt();
            assert_relative_eq!(radial, 3.5, epsilon = 1e-12);
        }
        // Exact first derivatives against differencing the exact points.
        let h = 1e-6;
        let (du, dv) = offset.d1_at(1.0, 2.0, T).unwrap();
        let fdu = (offset.point_at(1.0 + h, 2.0, T).unwrap()
            - offset.point_at(1.0 - h, 2.0, T).unwrap())
            / (2.0 * h);
        let fdv = (offset.point_at(1.0, 2.0 + h, T).unwrap()
            - offset.point_at(1.0, 2.0 - h, T).unwrap())
            / (2.0 * h);
        assert_relative_eq!((du - fdu).magnitude(), 0.0, epsilon = 1e-5);
        assert_relative_eq!((dv - fdv).magnitude(), 0.0, epsilon = 1e-5);
        // The second derivative refuses by name.
        assert!(offset.d2_at(1.0, 2.0, T).is_err());
    }

    fn tilted() -> Frame {
        Frame::new(
            Point::new(1.0, -2.0, 3.0),
            Direction::from_coords(1.0, 2.0, 3.0, T).unwrap(),
            Direction::X,
            T,
        )
        .unwrap()
    }

    fn patch() -> BSplineSurface {
        let (nu, nv) = (5, 4);
        let mut points = Vec::with_capacity(nu * nv);
        for i in 0..nu {
            for j in 0..nv {
                #[allow(clippy::cast_precision_loss)]
                let (x, y) = (i as f64, j as f64);
                points.push(Point::new(x, y, (x * 0.7).sin() * (y * 0.5).cos()));
            }
        }
        BSplineSurface::new(
            KnotVector::clamped_uniform(3, nu).unwrap(),
            KnotVector::clamped_uniform(2, nv).unwrap(),
            &ControlGrid::new(points, nu, nv).unwrap(),
            T,
        )
        .unwrap()
    }

    fn every_surface() -> Vec<SurfaceGeometry> {
        let circle: Curve = CircleCurve::new(
            Circle::new(
                Frame::new(Point::new(5.0, 0.0, 0.0), Direction::Y, Direction::X, T).unwrap(),
                1.0,
                T,
            )
            .unwrap(),
        )
        .into();
        let line: Curve =
            LineCurve::segment(Point::new(2.0, 0.0, 0.0), Point::new(2.0, 0.0, 4.0), T)
                .unwrap()
                .into();
        vec![
            PlaneSurface::over(Plane::new(tilted()), (-5.0, 5.0), (-3.0, 3.0))
                .unwrap()
                .into(),
            CylinderSurface::new(Cylinder::new(tilted(), 2.0, T).unwrap(), (-4.0, 4.0))
                .unwrap()
                .into(),
            ConeSurface::new(Cone::new(tilted(), 3.0, 0.6, T).unwrap(), (-1.0, 5.0))
                .unwrap()
                .into(),
            SphereSurface::new(Sphere::new(tilted(), 4.0, T).unwrap()).into(),
            TorusSurface::new(Torus::new(tilted(), 5.0, 2.0, T).unwrap()).into(),
            patch().into(),
            RevolutionSurface::new(circle, Axis::Z, TAU).unwrap().into(),
            ExtrusionSurface::new(line, Direction::X, 3.0)
                .unwrap()
                .into(),
            TrimmedSurface::new(patch().into(), (0.2, 0.8), (0.3, 0.7), T)
                .unwrap()
                .into(),
        ]
    }

    /// Interior sample parameters, avoiding the domain edges *and* the tidy
    /// fractions where a spline's interior knots sit.
    ///
    /// At a knot a spline's second derivative is one-sided, while a central
    /// difference straddles two different polynomial pieces — so a comparison
    /// there measures the discontinuity rather than the derivative. The offset
    /// keeps samples clear of knots at halves, thirds and quarters.
    fn interior(s: &SurfaceGeometry, n: usize) -> Vec<(f64, f64)> {
        const OFF: f64 = 0.0413;
        let ((ua, ub), (va, vb)) = s.domain();
        let mut out = Vec::new();
        for i in 1..n {
            for j in 1..n {
                #[allow(clippy::cast_precision_loss)]
                let (tu, tv) = (i as f64 / n as f64 + OFF, j as f64 / n as f64 + OFF);
                out.push((ua + (ub - ua) * tu, va + (vb - va) * tv));
            }
        }
        out
    }

    #[test]
    fn every_surfaces_first_partials_agree_with_finite_differences() {
        let h = 1e-6;
        for s in every_surface() {
            for (u, v) in interior(&s, 5) {
                let (du, dv) = s.d1_at(u, v, T).unwrap();
                let nu = (s.point_at(u + h, v, T).unwrap() - s.point_at(u - h, v, T).unwrap())
                    * (1.0 / (2.0 * h));
                let nv = (s.point_at(u, v + h, T).unwrap() - s.point_at(u, v - h, T).unwrap())
                    * (1.0 / (2.0 * h));
                assert!(
                    (du - nu).magnitude() <= 1e-5 * nu.magnitude().max(1.0),
                    "{:?} du at ({u}, {v}): {du:?} vs {nu:?}",
                    s.kind()
                );
                assert!(
                    (dv - nv).magnitude() <= 1e-5 * nv.magnitude().max(1.0),
                    "{:?} dv at ({u}, {v}): {dv:?} vs {nv:?}",
                    s.kind()
                );
            }
        }
    }

    #[test]
    fn every_surfaces_second_partials_agree_with_finite_differences() {
        let h = 1e-5;
        for s in every_surface() {
            for (u, v) in interior(&s, 4) {
                let (d2u, duv, d2v) = s.d2_at(u, v, T).unwrap();

                let nuu = (s.point_at(u + h, v, T).unwrap().to_vector()
                    - s.point_at(u, v, T).unwrap().to_vector() * 2.0
                    + s.point_at(u - h, v, T).unwrap().to_vector())
                    * (1.0 / (h * h));
                let nvv = (s.point_at(u, v + h, T).unwrap().to_vector()
                    - s.point_at(u, v, T).unwrap().to_vector() * 2.0
                    + s.point_at(u, v - h, T).unwrap().to_vector())
                    * (1.0 / (h * h));
                let nuv = (s.point_at(u + h, v + h, T).unwrap()
                    - s.point_at(u + h, v - h, T).unwrap()
                    - (s.point_at(u - h, v + h, T).unwrap()
                        - s.point_at(u - h, v - h, T).unwrap()))
                    * (1.0 / (4.0 * h * h));

                for (analytic, numeric, name) in
                    [(d2u, nuu, "d2u"), (duv, nuv, "duv"), (d2v, nvv, "d2v")]
                {
                    assert!(
                        (analytic - numeric).magnitude() <= 1e-3 * numeric.magnitude().max(1.0),
                        "{:?} {name} at ({u}, {v}): {analytic:?} vs {numeric:?}",
                        s.kind()
                    );
                }
            }
        }
    }

    #[test]
    fn out_of_domain_parameters_are_refused_where_the_surface_is_not_periodic() {
        for s in every_surface() {
            let ((ua, ub), (va, vb)) = s.domain();
            let inside = ((ua + ub) / 2.0, (va + vb) / 2.0);
            if s.is_periodic_u() {
                assert!(s.point_at(ub + 1.0, inside.1, T).is_ok(), "{:?}", s.kind());
            } else {
                assert!(s.point_at(ub + 1.0, inside.1, T).is_err(), "{:?}", s.kind());
            }
            if s.is_periodic_v() {
                assert!(s.point_at(inside.0, vb + 1.0, T).is_ok(), "{:?}", s.kind());
            } else {
                assert!(s.point_at(inside.0, vb + 1.0, T).is_err(), "{:?}", s.kind());
            }
        }
    }

    #[test]
    fn periodic_parameters_wrap_to_the_same_point() {
        for s in every_surface() {
            let ((ua, ub), (va, vb)) = s.domain();
            let (u, v) = ((ua + ub) * 0.4, (va + vb) * 0.4);
            let base = s.point_at(u, v, T).unwrap();
            if s.is_periodic_u() {
                let wrapped = s.point_at(u + (ub - ua), v, T).unwrap();
                assert!(base.is_equal(wrapped, T), "{:?} u wrap", s.kind());
            }
            if s.is_periodic_v() {
                let wrapped = s.point_at(u, v + (vb - va), T).unwrap();
                assert!(base.is_equal(wrapped, T), "{:?} v wrap", s.kind());
            }
        }
    }

    #[test]
    fn analytic_surfaces_contain_their_own_points() {
        // Cross-check against the independent distance functions in ogeom-math.
        let cyl = Cylinder::new(tilted(), 2.0, T).unwrap();
        let cone = Cone::new(tilted(), 3.0, 0.6, T).unwrap();
        let sph = Sphere::new(tilted(), 4.0, T).unwrap();
        let tor = Torus::new(tilted(), 5.0, 2.0, T).unwrap();
        let plane = Plane::new(tilted());

        let cyl_s = CylinderSurface::new(cyl, (-4.0, 4.0)).unwrap();
        let cone_s = ConeSurface::new(cone, (-1.0, 5.0)).unwrap();
        let sph_s = SphereSurface::new(sph);
        let tor_s = TorusSurface::new(tor);
        let plane_s = PlaneSurface::over(plane, (-5.0, 5.0), (-3.0, 3.0)).unwrap();

        for i in 1..6 {
            for j in 1..6 {
                let (tu, tv) = (f64::from(i) / 6.0, f64::from(j) / 6.0);
                assert!(
                    plane.contains(
                        plane_s
                            .point_at(-5.0 + 10.0 * tu, -3.0 + 6.0 * tv, T)
                            .unwrap(),
                        T
                    )
                );
                assert!(cyl.contains(cyl_s.point_at(TAU * tu, -4.0 + 8.0 * tv, T).unwrap(), T));
                assert!(cone.contains(cone_s.point_at(TAU * tu, -1.0 + 6.0 * tv, T).unwrap(), T));
                let half = core::f64::consts::FRAC_PI_2;
                assert!(
                    sph.contains(
                        sph_s
                            .point_at(TAU * tu, -half + core::f64::consts::PI * tv, T)
                            .unwrap(),
                        T
                    )
                );
                assert!(tor.contains(tor_s.point_at(TAU * tu, TAU * tv, T).unwrap(), T));
            }
        }
    }

    #[test]
    fn normals_of_analytic_surfaces_match_their_own_definitions() {
        let sph = Sphere::new(tilted(), 4.0, T).unwrap();
        let s = SphereSurface::new(sph);
        for i in 1..6 {
            for j in 1..6 {
                let u = TAU * f64::from(i) / 6.0;
                let v = -1.2 + 2.4 * f64::from(j) / 6.0;
                let p = s.point_at(u, v, T).unwrap();
                assert!(
                    s.normal_at(u, v, T)
                        .unwrap()
                        .is_equal(sph.normal_at(p, T).unwrap(), T)
                );
            }
        }
    }

    #[test]
    fn a_sphere_degenerates_at_its_poles_and_says_so() {
        let s = SphereSurface::new(Sphere::new(tilted(), 4.0, T).unwrap());
        let half = core::f64::consts::FRAC_PI_2;
        for pole in [-half, half] {
            assert!(s.is_degenerate_at(1.0, pole, T).unwrap());
            assert!(s.normal_at(1.0, pole, T).is_err());
        }
        assert!(!s.is_degenerate_at(1.0, 0.0, T).unwrap());
        assert!(s.normal_at(1.0, 0.0, T).is_ok());
    }

    #[test]
    fn a_cone_degenerates_at_its_apex() {
        let cone = Cone::new(tilted(), 3.0, 0.6, T).unwrap();
        let apex_height = -3.0 / 0.6_f64.tan();
        let s = ConeSurface::new(cone, (apex_height, 5.0)).unwrap();
        assert!(s.is_degenerate_at(1.0, apex_height, T).unwrap());
        assert!(s.normal_at(1.0, apex_height, T).is_err());
        assert!(
            s.point_at(1.0, apex_height, T)
                .unwrap()
                .is_equal(cone.apex(), T)
        );
    }

    #[test]
    fn revolving_a_circle_about_an_offset_axis_gives_a_torus() {
        // The strongest cross-check available here: two independent
        // constructions of the same surface must agree pointwise.
        let major = 5.0;
        let minor = 1.0;
        // The generator's normal is -Y, not +Y. A frame's second axis is
        // `z x x`, so a normal of +Y with an x of +X gives a y of -Z and winds
        // the circle backwards relative to the torus's v parameter. Getting
        // this wrong produces a surface that is the right shape and the wrong
        // parameterization, which is exactly the sort of mistake that survives
        // a visual check.
        let generator: Curve = CircleCurve::new(
            Circle::new(
                Frame::new(Point::new(major, 0.0, 0.0), -Direction::Y, Direction::X, T).unwrap(),
                minor,
                T,
            )
            .unwrap(),
        )
        .into();
        let revolved = RevolutionSurface::new(generator, Axis::Z, TAU).unwrap();
        let torus = TorusSurface::new(Torus::new(Frame::WORLD, major, minor, T).unwrap());

        for i in 0..12 {
            for j in 0..12 {
                let u = TAU * f64::from(i) / 12.0;
                let v = TAU * f64::from(j) / 12.0;
                let a = revolved.point_at(u, v, T).unwrap();
                let b = torus.point_at(u, v, T).unwrap();
                assert!(a.is_equal(b, T), "at ({u}, {v}): {a:?} vs {b:?}");
            }
        }
    }

    #[test]
    fn extruding_a_line_gives_a_plane_and_a_circle_gives_a_cylinder() {
        let line: Curve = LineCurve::segment(Point::ORIGIN, Point::new(0.0, 4.0, 0.0), T)
            .unwrap()
            .into();
        let sheet = ExtrusionSurface::new(line, Direction::Z, 3.0).unwrap();
        let plane = Plane::new(Frame::WORLD);
        for i in 0..=4 {
            for j in 0..=4 {
                let p = sheet
                    .point_at(4.0 * f64::from(i) / 4.0, 3.0 * f64::from(j) / 4.0, T)
                    .unwrap();
                assert!(plane.distance_to(p) < 1e-12 || p.x.abs() < 1e-12);
            }
        }

        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.0, T).unwrap()).into();
        let tube = ExtrusionSurface::new(circle, Direction::Z, 5.0).unwrap();
        let cylinder = Cylinder::new(Frame::WORLD, 2.0, T).unwrap();
        for i in 0..8 {
            for j in 0..=4 {
                let p = tube
                    .point_at(TAU * f64::from(i) / 8.0, 5.0 * f64::from(j) / 4.0, T)
                    .unwrap();
                assert!(cylinder.contains(p, T));
            }
        }
    }

    #[test]
    fn degenerate_constructions_are_refused() {
        assert!(PlaneSurface::over(Plane::new(tilted()), (1.0, 1.0), (0.0, 1.0)).is_err());
        assert!(
            CylinderSurface::new(Cylinder::new(tilted(), 1.0, T).unwrap(), (2.0, 1.0)).is_err()
        );
        let circle: Curve = CircleCurve::new(Circle::new(tilted(), 1.0, T).unwrap()).into();
        assert!(RevolutionSurface::new(circle.clone(), Axis::Z, 0.0).is_err());
        assert!(RevolutionSurface::new(circle.clone(), Axis::Z, 7.0).is_err());
        assert!(RevolutionSurface::new(circle.clone(), Axis::Z, TAU).is_ok());
        assert!(ExtrusionSurface::new(circle.clone(), Direction::Z, 0.0).is_err());
        assert!(ExtrusionSurface::new(circle, Direction::Z, f64::NAN).is_err());
    }

    #[test]
    fn trimming_is_bounds_checked_and_agrees_with_its_basis() {
        let basis: SurfaceGeometry = patch().into();
        assert!(TrimmedSurface::new(basis.clone(), (0.2, 0.8), (0.3, 0.7), T).is_ok());
        assert!(TrimmedSurface::new(basis.clone(), (0.8, 0.2), (0.3, 0.7), T).is_err());
        assert!(TrimmedSurface::new(basis.clone(), (-0.5, 0.8), (0.3, 0.7), T).is_err());

        let trimmed = TrimmedSurface::new(basis.clone(), (0.2, 0.8), (0.3, 0.7), T).unwrap();
        assert_eq!(trimmed.domain(), ((0.2, 0.8), (0.3, 0.7)));
        for i in 0..=4 {
            for j in 0..=4 {
                let u = 0.2 + 0.6 * f64::from(i) / 4.0;
                let v = 0.3 + 0.4 * f64::from(j) / 4.0;
                assert!(
                    trimmed
                        .point_at(u, v, T)
                        .unwrap()
                        .is_equal(basis.point_at(u, v, T).unwrap(), T)
                );
            }
        }
        assert!(trimmed.point_at(0.1, 0.5, T).is_err());
    }

    #[test]
    fn transforms_move_surfaces_and_preserve_their_kind() {
        let t =
            Transform::rotation(Axis::X, 0.7) * Transform::translation(Vector::new(1.0, 2.0, 3.0));
        for s in every_surface() {
            let moved = s.transformed(&t, T).unwrap();
            assert_eq!(moved.kind(), s.kind());
            for (u, v) in interior(&s, 4) {
                let expected = t.apply(s.point_at(u, v, T).unwrap());
                assert!(
                    moved.point_at(u, v, T).unwrap().is_equal(expected, T),
                    "{:?} at ({u}, {v})",
                    s.kind()
                );
            }
        }
    }

    #[test]
    fn a_scaling_rescales_length_valued_parameters() {
        // A cylinder's v is a height, so it must rescale; its u is an angle and
        // must not.
        let s: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 2.0, T).unwrap(), (0.0, 4.0))
                .unwrap()
                .into();
        let scaled = s
            .transformed(&Transform::scaling(Point::ORIGIN, 3.0, T).unwrap(), T)
            .unwrap();
        assert_eq!(scaled.domain(), ((0.0, TAU), (0.0, 12.0)));
        assert!(
            scaled
                .point_at(0.0, 12.0, T)
                .unwrap()
                .is_equal(Point::new(6.0, 0.0, 12.0), T)
        );
    }

    #[test]
    fn a_rational_patch_is_recognized_and_a_uniformly_weighted_one_is_not() {
        let g = patch();
        assert!(!g.is_rational());

        let w = core::f64::consts::FRAC_1_SQRT_2;
        let points: Vec<_> = [
            (Point::new(1.0, 0.0, 0.0), 1.0),
            (Point::new(1.0, 1.0, 0.0), w),
            (Point::new(0.0, 1.0, 0.0), 1.0),
            (Point::new(1.0, 0.0, 1.0), 1.0),
            (Point::new(1.0, 1.0, 1.0), w),
            (Point::new(0.0, 1.0, 1.0), 1.0),
        ]
        .iter()
        .map(|(p, w)| Weighted::new(*p, *w, T).unwrap())
        .collect();
        let arc = BSplineSurface::rational(
            KnotVector::clamped_uniform(1, 2).unwrap(),
            KnotVector::clamped_uniform(2, 3).unwrap(),
            ControlGrid::new(points, 2, 3).unwrap(),
        )
        .unwrap();
        assert!(arc.is_rational());
        // Each isoparametric line is an exact unit arc in the xy plane.
        for i in 0..=4 {
            for j in 0..=8 {
                let p = arc
                    .point_at(f64::from(i) / 4.0, f64::from(j) / 8.0, T)
                    .unwrap();
                assert_relative_eq!(p.x.hypot(p.y), 1.0, epsilon = 1e-14);
            }
        }
    }

    #[test]
    fn closure_and_periodicity_are_reported_correctly() {
        let s = every_surface();
        let expect = [
            (SurfaceKind::Plane, false, false, false, false),
            (SurfaceKind::Cylinder, true, false, true, false),
            (SurfaceKind::Cone, true, false, true, false),
            (SurfaceKind::Sphere, true, false, true, false),
            (SurfaceKind::Torus, true, true, true, true),
            (SurfaceKind::BSpline, false, false, false, false),
            (SurfaceKind::Revolution, true, true, true, true),
            (SurfaceKind::Extrusion, false, false, false, false),
            (SurfaceKind::Trimmed, false, false, false, false),
        ];
        for (surface, (kind, cu, cv, pu, pv)) in s.iter().zip(expect) {
            assert_eq!(surface.kind(), kind);
            assert_eq!(surface.is_closed_u(T), cu, "{kind:?} closed u");
            assert_eq!(surface.is_closed_v(T), cv, "{kind:?} closed v");
            assert_eq!(surface.is_periodic_u(), pu, "{kind:?} periodic u");
            assert_eq!(surface.is_periodic_v(), pv, "{kind:?} periodic v");
        }
    }

    #[test]
    fn a_sphere_is_not_periodic_in_latitude() {
        // Wrapping past a pole would land on the far side of the sphere, which
        // is a different point — so latitude stops rather than repeating.
        let s = SphereSurface::new(Sphere::new(Frame::WORLD, 1.0, T).unwrap());
        assert!(!s.is_periodic_v());
        assert!(s.point_at(0.0, core::f64::consts::PI, T).is_err());
    }
}
