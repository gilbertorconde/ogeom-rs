//! Elementary surfaces: plane, cylinder, cone, sphere, torus.
//!
//! Each is described by a [`Frame`] and its size parameters. As with the conics
//! the frame is not decoration: it fixes the parameterization, and therefore
//! fixes where a cylinder's seam falls and which way its normal points.
//!
//! These five plus the plane cover the overwhelming majority of real mechanical
//! geometry. Keeping them as exact analytic descriptions rather than converting
//! everything to NURBS is what lets intersection take analytic shortcuts, lets
//! measurement report a radius rather than a fitted approximation of one, and
//! keeps files small.
//!
//! Evaluation and derivatives are in [`crate::elementary`]; this module holds
//! the descriptions and the queries that follow directly from them.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};

use crate::{Axis, Direction, Frame, Point, Transform, Vector};

/// Reject a size parameter that cannot describe a real shape.
fn check_positive(name: &str, value: f64, tol: Tolerances) -> OgeomResult<()> {
    if !value.is_finite() || value <= tol.confusion() {
        ogeom_bail!(Construction, "{name} {value} must be finite and positive");
    }
    Ok(())
}

/// An unbounded plane.
///
/// The frame's `z` is the normal; `x` and `y` span the surface and fix its
/// parameterization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    frame: Frame,
}

/// An unbounded circular cylinder, with the frame's `z` as its axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cylinder {
    frame: Frame,
    radius: f64,
}

/// An unbounded circular cone.
///
/// The frame's `z` is the axis and its origin sits on the reference circle, of
/// radius [`Cone::reference_radius`]. The radius grows in `+z` for a positive
/// half angle. The apex is the point where it reaches zero.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cone {
    frame: Frame,
    reference_radius: f64,
    half_angle: f64,
}

/// A sphere.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sphere {
    frame: Frame,
    radius: f64,
}

/// A torus, with the frame's `z` as its axis of revolution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Torus {
    frame: Frame,
    major_radius: f64,
    minor_radius: f64,
}

/// How a torus's minor radius compares with its major radius.
///
/// The three cases are genuinely different surfaces, and an algorithm that
/// assumes the first will produce nonsense on the others: a spindle torus
/// self-intersects, and a horn torus is tangent to itself at the poles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TorusKind {
    /// `minor < major`: a ring, with a hole.
    Ring,
    /// `minor == major`: the hole closes to a single point at each pole.
    Horn,
    /// `minor > major`: the surface passes through itself.
    Spindle,
}

impl Plane {
    /// The `xy` plane.
    pub const XY: Self = Self {
        frame: Frame::WORLD,
    };

    /// The plane of `frame`, with `frame`'s `z` as its normal.
    #[must_use]
    pub const fn new(frame: Frame) -> Self {
        Self { frame }
    }

    /// The plane through `origin` with the given `normal`, parameterized
    /// arbitrarily but deterministically.
    #[must_use]
    pub fn through(origin: Point, normal: Direction) -> Self {
        Self {
            frame: Frame::about(origin, normal),
        }
    }

    /// The plane through three points, with the normal following the right-hand
    /// rule around `a`, `b`, `c`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the points
    /// are collinear.
    pub fn through_points(a: Point, b: Point, c: Point, tol: Tolerances) -> OgeomResult<Self> {
        let normal = Direction::from_cross(b - a, c - a, tol)?;
        Ok(Self::through(a, normal))
    }

    /// The frame positioning this plane.
    #[must_use]
    pub const fn frame(&self) -> Frame {
        self.frame
    }

    /// A point on the plane.
    #[must_use]
    pub const fn origin(&self) -> Point {
        self.frame.origin()
    }

    /// The normal.
    #[must_use]
    pub const fn normal(&self) -> Direction {
        self.frame.z()
    }

    /// The signed distance from `p`, positive on the side the normal points to.
    #[must_use]
    pub fn signed_distance_to(&self, p: Point) -> f64 {
        self.frame.signed_distance_to_plane(p)
    }

    /// The distance from `p`.
    #[must_use]
    pub fn distance_to(&self, p: Point) -> f64 {
        self.signed_distance_to(p).abs()
    }

    /// The closest point on the plane to `p`.
    #[must_use]
    pub fn project(&self, p: Point) -> Point {
        p - self.normal() * self.signed_distance_to(p)
    }

    /// Whether `p` lies on the plane within `tol.confusion()`.
    #[must_use]
    pub fn contains(&self, p: Point, tol: Tolerances) -> bool {
        self.distance_to(p) <= tol.confusion()
    }

    /// This plane with its normal reversed.
    #[must_use]
    pub const fn reversed(&self) -> Self {
        Self {
            frame: self.frame.with_z_reversed(),
        }
    }

    /// This plane moved by `t`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// transformed frame is degenerate.
    pub fn transformed(&self, t: &Transform, tol: Tolerances) -> OgeomResult<Self> {
        Ok(Self::new(t.apply_frame(&self.frame, tol)?))
    }
}

impl Cylinder {
    /// A cylinder of `radius` about `frame`'s `z` axis.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `radius` is
    /// not finite and positive.
    pub fn new(frame: Frame, radius: f64, tol: Tolerances) -> OgeomResult<Self> {
        check_positive("cylinder radius", radius, tol)?;
        Ok(Self { frame, radius })
    }

    /// A cylinder about an axis, parameterized arbitrarily but
    /// deterministically.
    ///
    /// # Errors
    ///
    /// As [`Cylinder::new`].
    pub fn about(axis: Axis, radius: f64, tol: Tolerances) -> OgeomResult<Self> {
        Self::new(Frame::about(axis.location, axis.direction), radius, tol)
    }

    /// The frame positioning this cylinder.
    #[must_use]
    pub const fn frame(&self) -> Frame {
        self.frame
    }

    /// The axis of revolution.
    #[must_use]
    pub const fn axis(&self) -> Axis {
        self.frame.axis()
    }

    /// The radius.
    #[must_use]
    pub const fn radius(&self) -> f64 {
        self.radius
    }

    /// The signed distance from `p`, negative inside.
    #[must_use]
    pub fn signed_distance_to(&self, p: Point) -> f64 {
        self.axis().distance_to(p) - self.radius
    }

    /// The distance from `p` to the surface.
    #[must_use]
    pub fn distance_to(&self, p: Point) -> f64 {
        self.signed_distance_to(p).abs()
    }

    /// Whether `p` lies on the surface within `tol.confusion()`.
    #[must_use]
    pub fn contains(&self, p: Point, tol: Tolerances) -> bool {
        self.distance_to(p) <= tol.confusion()
    }

    /// The outward unit normal at the point of the surface nearest `p`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `p` lies on
    /// the axis, where the nearest point — and so the normal — is not unique.
    pub fn normal_at(&self, p: Point, tol: Tolerances) -> OgeomResult<Direction> {
        let axis = self.axis();
        Direction::new(p - axis.project(p), tol)
    }

    /// The area of a section of this cylinder `height` long.
    #[must_use]
    pub fn lateral_area(&self, height: f64) -> f64 {
        core::f64::consts::TAU * self.radius * height.abs()
    }

    /// The volume enclosed by a section `height` long.
    #[must_use]
    pub fn volume(&self, height: f64) -> f64 {
        core::f64::consts::PI * self.radius * self.radius * height.abs()
    }

    /// This cylinder moved by `t`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// transformed frame is degenerate.
    pub fn transformed(&self, t: &Transform, tol: Tolerances) -> OgeomResult<Self> {
        Self::new(
            t.apply_frame(&self.frame, tol)?,
            self.radius * t.scale_factor().abs(),
            tol,
        )
    }
}

impl Cone {
    /// A cone with the given reference radius and half angle.
    ///
    /// The reference circle lies in `frame`'s `xy` plane. A positive half angle
    /// widens the cone in `+z`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
    /// `reference_radius` is negative or non-finite, or if `half_angle` is not
    /// strictly between `0` and `pi/2`. At zero the cone is a cylinder and at
    /// `pi/2` it is a plane; both are different surfaces with their own types,
    /// and admitting them here would produce a cone whose apex is at infinity.
    pub fn new(
        frame: Frame,
        reference_radius: f64,
        half_angle: f64,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        if !reference_radius.is_finite() || reference_radius < 0.0 {
            ogeom_bail!(
                Construction,
                "cone reference radius {reference_radius} must be finite and non-negative"
            );
        }
        if !half_angle.is_finite()
            || half_angle.abs() <= tol.angular()
            || half_angle.abs() >= core::f64::consts::FRAC_PI_2 - tol.angular()
        {
            ogeom_bail!(
                Construction,
                "cone half angle {half_angle} must lie strictly between 0 and pi/2"
            );
        }
        Ok(Self {
            frame,
            reference_radius,
            half_angle,
        })
    }

    /// The frame positioning this cone.
    #[must_use]
    pub const fn frame(&self) -> Frame {
        self.frame
    }

    /// The axis of revolution.
    #[must_use]
    pub const fn axis(&self) -> Axis {
        self.frame.axis()
    }

    /// The radius of the circle in the frame's `xy` plane.
    #[must_use]
    pub const fn reference_radius(&self) -> f64 {
        self.reference_radius
    }

    /// The half angle at the apex, in `(0, pi/2)`.
    #[must_use]
    pub const fn half_angle(&self) -> f64 {
        self.half_angle
    }

    /// The apex.
    #[must_use]
    pub fn apex(&self) -> Point {
        // The radius shrinks at `tan(half_angle)` per unit along the axis, so
        // the apex is that many units back from the reference circle.
        self.frame.origin() - self.frame.z() * (self.reference_radius / self.half_angle.tan())
    }

    /// The radius at signed distance `z` along the axis from the frame origin.
    #[must_use]
    pub fn radius_at(&self, z: f64) -> f64 {
        self.half_angle.tan().mul_add(z, self.reference_radius)
    }

    /// The distance from `p` to the surface, ignoring the far nappe.
    ///
    /// A double cone extends both sides of its apex; this measures to the
    /// surface as a whole, which is what a surface query means.
    #[must_use]
    pub fn distance_to(&self, p: Point) -> f64 {
        let local = self.frame.to_local(p);
        let radial = local.xy().to_vector().magnitude();
        let apex_z = -self.reference_radius / self.half_angle.tan();
        // A cone is a *double* cone: the quadric has two nappes meeting at the
        // apex, and the surface type built on this parameterizes both — its
        // height range may cross the apex, exactly as the conventional
        // kernel's conical surface does. An earlier version measured one nappe
        // and clamped everything past the apex to the apex, which reported a
        // point *on* the second nappe as almost a unit away — and it was the
        // intersection benchmark that caught it, by flagging a correctly
        // traced curve as off the surface.
        //
        // In the (radial, axial) half-plane each nappe is a ray from the apex;
        // the distance is the nearer of the two, each clamped to its own ray
        // so a point in the wedge beyond the apex measures to the apex.
        let (sin, cos) = self.half_angle.sin_cos();
        let height = local.z - apex_z;
        let apex_distance = radial.hypot(height);
        let nappe = |along: f64, across: f64| {
            if along <= 0.0 {
                apex_distance
            } else {
                across.abs()
            }
        };
        let up = nappe(
            height.mul_add(cos, radial * sin),
            height.mul_add(sin, -(radial * cos)),
        );
        let down = nappe(
            height.mul_add(-cos, radial * sin),
            height.mul_add(sin, radial * cos),
        );
        up.min(down)
    }

    /// Whether `p` lies on the surface within `tol.confusion()`.
    #[must_use]
    pub fn contains(&self, p: Point, tol: Tolerances) -> bool {
        self.distance_to(p) <= tol.confusion()
    }

    /// This cone moved by `t`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// transformed frame is degenerate.
    pub fn transformed(&self, t: &Transform, tol: Tolerances) -> OgeomResult<Self> {
        // A similarity scales lengths uniformly, so the half angle survives it
        // unchanged. That is exactly why the transform type is restricted to
        // similarities: a non-uniform scale would leave a surface that is no
        // longer a circular cone at all.
        Self::new(
            t.apply_frame(&self.frame, tol)?,
            self.reference_radius * t.scale_factor().abs(),
            self.half_angle,
            tol,
        )
    }
}

impl Sphere {
    /// A sphere of `radius` centred on `frame`'s origin.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `radius` is
    /// not finite and positive.
    pub fn new(frame: Frame, radius: f64, tol: Tolerances) -> OgeomResult<Self> {
        check_positive("sphere radius", radius, tol)?;
        Ok(Self { frame, radius })
    }

    /// A sphere from a centre and a radius, parameterized arbitrarily but
    /// deterministically.
    ///
    /// # Errors
    ///
    /// As [`Sphere::new`].
    pub fn centred(centre: Point, radius: f64, tol: Tolerances) -> OgeomResult<Self> {
        Self::new(Frame::about(centre, Direction::Z), radius, tol)
    }

    /// The frame positioning this sphere.
    #[must_use]
    pub const fn frame(&self) -> Frame {
        self.frame
    }

    /// The centre.
    #[must_use]
    pub const fn centre(&self) -> Point {
        self.frame.origin()
    }

    /// The radius.
    #[must_use]
    pub const fn radius(&self) -> f64 {
        self.radius
    }

    /// The signed distance from `p`, negative inside.
    #[must_use]
    pub fn signed_distance_to(&self, p: Point) -> f64 {
        self.centre().distance(p) - self.radius
    }

    /// The distance from `p` to the surface.
    #[must_use]
    pub fn distance_to(&self, p: Point) -> f64 {
        self.signed_distance_to(p).abs()
    }

    /// Whether `p` lies on the surface within `tol.confusion()`.
    #[must_use]
    pub fn contains(&self, p: Point, tol: Tolerances) -> bool {
        self.distance_to(p) <= tol.confusion()
    }

    /// Whether `p` lies strictly inside.
    #[must_use]
    pub fn encloses(&self, p: Point, tol: Tolerances) -> bool {
        self.signed_distance_to(p) < -tol.confusion()
    }

    /// The outward unit normal at the point nearest `p`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `p` is the
    /// centre, where no nearest point is unique.
    pub fn normal_at(&self, p: Point, tol: Tolerances) -> OgeomResult<Direction> {
        Direction::new(p - self.centre(), tol)
    }

    /// The surface area.
    #[must_use]
    pub fn area(&self) -> f64 {
        4.0 * core::f64::consts::PI * self.radius * self.radius
    }

    /// The volume enclosed.
    #[must_use]
    pub fn volume(&self) -> f64 {
        4.0 / 3.0 * core::f64::consts::PI * self.radius.powi(3)
    }

    /// This sphere moved by `t`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// transformed frame is degenerate.
    pub fn transformed(&self, t: &Transform, tol: Tolerances) -> OgeomResult<Self> {
        Self::new(
            t.apply_frame(&self.frame, tol)?,
            self.radius * t.scale_factor().abs(),
            tol,
        )
    }
}

impl Torus {
    /// A torus about `frame`'s `z` axis.
    ///
    /// `major_radius` is the distance from the axis to the centre of the tube;
    /// `minor_radius` is the tube's own radius.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if either
    /// radius is not finite and positive. A minor radius exceeding the major is
    /// *allowed*: that is a spindle torus, a real self-intersecting surface, and
    /// [`Torus::kind`] reports which case this is.
    pub fn new(
        frame: Frame,
        major_radius: f64,
        minor_radius: f64,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        check_positive("major radius", major_radius, tol)?;
        check_positive("minor radius", minor_radius, tol)?;
        Ok(Self {
            frame,
            major_radius,
            minor_radius,
        })
    }

    /// The frame positioning this torus.
    #[must_use]
    pub const fn frame(&self) -> Frame {
        self.frame
    }

    /// The centre.
    #[must_use]
    pub const fn centre(&self) -> Point {
        self.frame.origin()
    }

    /// The axis of revolution.
    #[must_use]
    pub const fn axis(&self) -> Axis {
        self.frame.axis()
    }

    /// The distance from the axis to the centre of the tube.
    #[must_use]
    pub const fn major_radius(&self) -> f64 {
        self.major_radius
    }

    /// The radius of the tube.
    #[must_use]
    pub const fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// Which of the three topological cases this torus falls into.
    #[must_use]
    pub fn kind(&self, tol: Tolerances) -> TorusKind {
        let difference = self.major_radius - self.minor_radius;
        if difference.abs() <= tol.confusion() {
            TorusKind::Horn
        } else if difference > 0.0 {
            TorusKind::Ring
        } else {
            TorusKind::Spindle
        }
    }

    /// Whether the surface passes through itself.
    #[must_use]
    pub fn self_intersects(&self, tol: Tolerances) -> bool {
        self.kind(tol) == TorusKind::Spindle
    }

    /// The signed distance from `p`, negative inside the tube.
    ///
    /// Meaningful for a [`TorusKind::Ring`] or [`TorusKind::Horn`] torus, where
    /// inside and outside are well defined. A [`TorusKind::Spindle`] torus
    /// passes through itself and has no consistent interior, so use
    /// [`Torus::distance_to`] there.
    #[must_use]
    pub fn signed_distance_to(&self, p: Point) -> f64 {
        let local = self.frame.to_local(p);
        // Distance to the tube's centre circle, then out by the tube radius.
        let radial = local.xy().to_vector().magnitude() - self.major_radius;
        radial.hypot(local.z) - self.minor_radius
    }

    /// The distance from `p` to the surface.
    ///
    /// Correct for all three kinds. In the half-plane at a fixed azimuth, the
    /// surface's profile is the generating circle *folded* about the axis: when
    /// the minor radius exceeds the major, part of that circle lies on the far
    /// side of the axis and sweeps to the near side. Measuring only to the
    /// unfolded branch — which is what the signed form does — then reports a
    /// point on the surface as being some distance off it.
    #[must_use]
    pub fn distance_to(&self, p: Point) -> f64 {
        let local = self.frame.to_local(p);
        let radial = local.xy().to_vector().magnitude();
        let near = (radial - self.major_radius).hypot(local.z) - self.minor_radius;
        let folded = (radial + self.major_radius).hypot(local.z) - self.minor_radius;
        near.abs().min(folded.abs())
    }

    /// Whether `p` lies on the surface within `tol.confusion()`.
    #[must_use]
    pub fn contains(&self, p: Point, tol: Tolerances) -> bool {
        self.distance_to(p) <= tol.confusion()
    }

    /// The surface area, for a ring torus.
    ///
    /// Meaningless for a spindle torus, whose surface overlaps itself.
    #[must_use]
    pub fn area(&self) -> f64 {
        4.0 * core::f64::consts::PI * core::f64::consts::PI * self.major_radius * self.minor_radius
    }

    /// The volume enclosed, for a ring torus.
    #[must_use]
    pub fn volume(&self) -> f64 {
        2.0 * core::f64::consts::PI
            * core::f64::consts::PI
            * self.major_radius
            * self.minor_radius
            * self.minor_radius
    }

    /// This torus moved by `t`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// transformed frame is degenerate.
    pub fn transformed(&self, t: &Transform, tol: Tolerances) -> OgeomResult<Self> {
        let s = t.scale_factor().abs();
        Self::new(
            t.apply_frame(&self.frame, tol)?,
            self.major_radius * s,
            self.minor_radius * s,
            tol,
        )
    }
}

/// The unit normal to a plane through `origin` containing `a` and `b`.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `a` and `b` are
/// collinear.
pub fn plane_normal(a: Vector, b: Vector, tol: Tolerances) -> OgeomResult<Direction> {
    Direction::from_cross(a, b, tol)
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn plane_distance_is_signed_by_the_normal() {
        let p = Plane::XY;
        assert_relative_eq!(p.signed_distance_to(Point::new(1.0, 2.0, 3.0)), 3.0);
        assert_relative_eq!(p.signed_distance_to(Point::new(1.0, 2.0, -3.0)), -3.0);
        assert_relative_eq!(p.distance_to(Point::new(0.0, 0.0, -3.0)), 3.0);
        assert!(p.contains(Point::new(9.0, -9.0, 0.0), T));
        assert!(
            p.project(Point::new(1.0, 2.0, 3.0))
                .is_equal(Point::new(1.0, 2.0, 0.0), T)
        );
    }

    #[test]
    fn reversing_a_plane_flips_the_sign_but_not_the_surface() {
        let p = Plane::XY;
        let r = p.reversed();
        let q = Point::new(0.0, 0.0, 5.0);
        assert_relative_eq!(r.signed_distance_to(q), -p.signed_distance_to(q));
        assert_relative_eq!(r.distance_to(q), p.distance_to(q));
        assert!(r.contains(Point::new(1.0, 1.0, 0.0), T));
    }

    #[test]
    fn plane_through_three_points_follows_the_right_hand_rule() {
        let p = Plane::through_points(
            Point::ORIGIN,
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            T,
        )
        .unwrap();
        assert!(p.normal().is_equal(Direction::Z, T));
        // Reversing the winding reverses the normal.
        let q = Plane::through_points(
            Point::ORIGIN,
            Point::new(0.0, 1.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            T,
        )
        .unwrap();
        assert!(q.normal().is_equal(-Direction::Z, T));
        assert!(
            Plane::through_points(
                Point::ORIGIN,
                Point::new(1.0, 1.0, 1.0),
                Point::new(2.0, 2.0, 2.0),
                T
            )
            .is_err()
        );
    }

    #[test]
    fn plane_through_tiny_triangles_still_works() {
        // Same trap as the circumcircle: the cross product scales as the square
        // of the triangle size.
        let s = 1e-6;
        let p = Plane::through_points(
            Point::ORIGIN,
            Point::new(s, 0.0, 0.0),
            Point::new(0.0, s, 0.0),
            T,
        )
        .unwrap();
        assert!(p.normal().is_equal(Direction::Z, T));
    }

    #[test]
    fn cylinder_distance_and_normal() {
        let c = Cylinder::new(Frame::WORLD, 2.0, T).unwrap();
        assert_relative_eq!(c.signed_distance_to(Point::new(3.0, 0.0, 100.0)), 1.0);
        assert_relative_eq!(c.signed_distance_to(Point::new(1.0, 0.0, -50.0)), -1.0);
        assert_relative_eq!(c.signed_distance_to(Point::ORIGIN), -2.0);
        assert!(c.contains(Point::new(0.0, 2.0, 7.0), T));
        assert!(
            c.normal_at(Point::new(3.0, 0.0, 5.0), T)
                .unwrap()
                .is_equal(Direction::X, T)
        );
        // A point on the axis has no unique normal.
        assert!(c.normal_at(Point::new(0.0, 0.0, 5.0), T).is_err());
    }

    #[test]
    fn cylinder_measurements() {
        let c = Cylinder::new(Frame::WORLD, 2.0, T).unwrap();
        assert_relative_eq!(c.lateral_area(5.0), core::f64::consts::TAU * 10.0);
        assert_relative_eq!(c.volume(5.0), core::f64::consts::PI * 20.0);
        // Height is a magnitude; a negative one is the same section.
        assert_relative_eq!(c.volume(-5.0), c.volume(5.0));
    }

    #[test]
    fn cone_degenerate_angles_are_refused() {
        let f = Frame::WORLD;
        assert!(
            Cone::new(f, 1.0, 0.0, T).is_err(),
            "zero angle is a cylinder"
        );
        assert!(
            Cone::new(f, 1.0, core::f64::consts::FRAC_PI_2, T).is_err(),
            "a right angle is a plane"
        );
        assert!(Cone::new(f, 1.0, f64::NAN, T).is_err());
        assert!(Cone::new(f, -1.0, 0.5, T).is_err());
        // A zero reference radius is fine: the frame origin is then the apex.
        assert!(Cone::new(f, 0.0, 0.5, T).is_ok());
        assert!(Cone::new(f, 1.0, 0.5, T).is_ok());
    }

    #[test]
    fn cone_apex_and_radius_profile() {
        // Half angle of 45 degrees: radius grows one unit per unit of height.
        let quarter = core::f64::consts::FRAC_PI_4;
        let c = Cone::new(Frame::WORLD, 3.0, quarter, T).unwrap();
        assert!(c.apex().is_equal(Point::new(0.0, 0.0, -3.0), T));
        assert_relative_eq!(c.radius_at(0.0), 3.0, epsilon = 1e-12);
        assert_relative_eq!(c.radius_at(2.0), 5.0, epsilon = 1e-12);
        assert_relative_eq!(c.radius_at(-3.0), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn cone_distance_is_zero_on_the_surface() {
        let quarter = core::f64::consts::FRAC_PI_4;
        let c = Cone::new(Frame::WORLD, 3.0, quarter, T).unwrap();
        for z in [-3.0_f64, -1.0, 0.0, 2.0, 10.0] {
            let r = c.radius_at(z);
            for angle in [0.0_f64, 1.0, 2.5] {
                let p = Point::new(r * angle.cos(), r * angle.sin(), z);
                assert!(c.contains(p, T), "z = {z}, distance = {}", c.distance_to(p));
            }
        }
    }

    #[test]
    fn cone_distance_off_the_surface_is_perpendicular() {
        let quarter = core::f64::consts::FRAC_PI_4;
        let c = Cone::new(Frame::WORLD, 0.0, quarter, T).unwrap();
        // Apex at the origin, opening along +z at 45 degrees. The point (1,0,0)
        // sits at perpendicular distance sin(45) from the surface line.
        assert_relative_eq!(
            c.distance_to(Point::new(1.0, 0.0, 0.0)),
            core::f64::consts::FRAC_1_SQRT_2,
            epsilon = 1e-12
        );
        // A cone is a double cone: behind the apex is the second nappe, and a
        // point on the axis there measures perpendicular to it, not to the
        // apex. The earlier claim here — apex distance, 5.0 — encoded a
        // single-nappe convention that disagreed with the surface type built
        // on this, and the intersection benchmark caught the disagreement by
        // flagging a correctly traced second-nappe curve as off the surface.
        assert_relative_eq!(
            c.distance_to(Point::new(0.0, 0.0, -5.0)),
            5.0 * core::f64::consts::FRAC_1_SQRT_2,
            epsilon = 1e-12
        );
        // A point *on* the second nappe is on the cone.
        assert_relative_eq!(
            c.distance_to(Point::new(2.0, 0.0, -2.0)),
            0.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn sphere_queries() {
        let s = Sphere::centred(Point::new(1.0, 2.0, 3.0), 5.0, T).unwrap();
        assert_relative_eq!(s.signed_distance_to(Point::new(1.0, 2.0, 3.0)), -5.0);
        assert_relative_eq!(s.signed_distance_to(Point::new(6.0, 2.0, 3.0)), 0.0);
        assert_relative_eq!(s.signed_distance_to(Point::new(11.0, 2.0, 3.0)), 5.0);
        assert!(s.encloses(s.centre(), T));
        assert!(!s.encloses(Point::new(6.0, 2.0, 3.0), T));
        assert!(s.contains(Point::new(6.0, 2.0, 3.0), T));
        assert!(
            s.normal_at(Point::new(6.0, 2.0, 3.0), T)
                .unwrap()
                .is_equal(Direction::X, T)
        );
        assert!(s.normal_at(s.centre(), T).is_err());
    }

    #[test]
    fn sphere_measurements() {
        let s = Sphere::centred(Point::ORIGIN, 3.0, T).unwrap();
        assert_relative_eq!(s.area(), 4.0 * core::f64::consts::PI * 9.0);
        assert_relative_eq!(s.volume(), 4.0 / 3.0 * core::f64::consts::PI * 27.0);
    }

    #[test]
    fn torus_kinds_are_distinguished() {
        let f = Frame::WORLD;
        assert_eq!(Torus::new(f, 5.0, 1.0, T).unwrap().kind(T), TorusKind::Ring);
        assert_eq!(Torus::new(f, 5.0, 5.0, T).unwrap().kind(T), TorusKind::Horn);
        assert_eq!(
            Torus::new(f, 5.0, 8.0, T).unwrap().kind(T),
            TorusKind::Spindle
        );
        assert!(Torus::new(f, 5.0, 8.0, T).unwrap().self_intersects(T));
        assert!(!Torus::new(f, 5.0, 1.0, T).unwrap().self_intersects(T));
        // A spindle torus is a real surface and must be constructible.
        assert!(Torus::new(f, 1.0, 2.0, T).is_ok());
        assert!(Torus::new(f, 0.0, 1.0, T).is_err());
    }

    #[test]
    fn spindle_torus_distance_accounts_for_the_folded_branch() {
        // Minor radius exceeds major: the generating circle crosses the axis,
        // so part of the surface comes from the far side of it. A point there
        // is on the surface, and the unfolded formula alone would not say so.
        let t = Torus::new(Frame::WORLD, 1.0, 3.0, T).unwrap();
        assert_eq!(t.kind(T), TorusKind::Spindle);

        // Parametric point with `major + minor*cos(v) < 0`, which lands on the
        // folded branch.
        let v = core::f64::consts::PI;
        let radial = 1.0 + 3.0 * v.cos(); // = -2
        let p = Point::new(radial.abs(), 0.0, 3.0 * v.sin());
        assert!(t.contains(p, T), "distance was {}", t.distance_to(p));

        // The ring case is unaffected: near branch still wins everywhere.
        let ring = Torus::new(Frame::WORLD, 5.0, 2.0, T).unwrap();
        for p in [
            Point::new(7.0, 0.0, 0.0),
            Point::new(3.0, 0.0, 0.0),
            Point::new(5.0, 0.0, 2.0),
        ] {
            assert_relative_eq!(ring.distance_to(p), 0.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn torus_distance_and_measurements() {
        let t = Torus::new(Frame::WORLD, 5.0, 2.0, T).unwrap();
        // Outer equator, inner equator, and the top of the tube.
        assert!(t.contains(Point::new(7.0, 0.0, 0.0), T));
        assert!(t.contains(Point::new(3.0, 0.0, 0.0), T));
        assert!(t.contains(Point::new(5.0, 0.0, 2.0), T));
        // The centre of the tube is one tube radius inside.
        assert_relative_eq!(t.signed_distance_to(Point::new(5.0, 0.0, 0.0)), -2.0);
        // The centre of the hole is far outside the surface.
        assert_relative_eq!(t.signed_distance_to(Point::ORIGIN), 3.0);

        let pi2 = core::f64::consts::PI * core::f64::consts::PI;
        assert_relative_eq!(t.area(), 4.0 * pi2 * 10.0);
        assert_relative_eq!(t.volume(), 2.0 * pi2 * 5.0 * 4.0);
    }

    #[test]
    fn transforms_scale_sizes_and_preserve_shape() {
        let scale = Transform::scaling(Point::ORIGIN, 3.0, T).unwrap();

        let c = Cylinder::new(Frame::WORLD, 2.0, T).unwrap();
        assert_relative_eq!(
            c.transformed(&scale, T).unwrap().radius(),
            6.0,
            epsilon = 1e-12
        );

        let s = Sphere::centred(Point::new(1.0, 0.0, 0.0), 2.0, T).unwrap();
        let moved = s.transformed(&scale, T).unwrap();
        assert_relative_eq!(moved.radius(), 6.0, epsilon = 1e-12);
        assert!(moved.centre().is_equal(Point::new(3.0, 0.0, 0.0), T));

        let t = Torus::new(Frame::WORLD, 5.0, 2.0, T).unwrap();
        let scaled = t.transformed(&scale, T).unwrap();
        assert_relative_eq!(scaled.major_radius(), 15.0, epsilon = 1e-12);
        assert_relative_eq!(scaled.minor_radius(), 6.0, epsilon = 1e-12);

        // A similarity leaves a cone's half angle alone — the reason transforms
        // are restricted to similarities in the first place.
        let cone = Cone::new(Frame::WORLD, 3.0, 0.5, T).unwrap();
        let big = cone.transformed(&scale, T).unwrap();
        assert_relative_eq!(big.half_angle(), 0.5);
        assert_relative_eq!(big.reference_radius(), 9.0, epsilon = 1e-12);
    }

    #[test]
    fn transformed_surfaces_still_contain_their_transformed_points() {
        let t =
            Transform::rotation(Axis::X, 0.7) * Transform::translation(Vector::new(1.0, 2.0, 3.0));
        let s = Sphere::centred(Point::ORIGIN, 4.0, T).unwrap();
        let moved = s.transformed(&t, T).unwrap();
        let on_surface = Point::new(4.0, 0.0, 0.0);
        assert!(s.contains(on_surface, T));
        assert!(moved.contains(t.apply(on_surface), T));
    }
}
