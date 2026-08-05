//! Axes and coordinate frames.
//!
//! An [`Axis`] is a point and a direction. A [`Frame`] is a full local
//! coordinate system: an origin plus three mutually perpendicular directions
//! plus a handedness.
//!
//! Frames are how every piece of analytic geometry in the kernel is positioned.
//! A cylinder is a radius and a frame; a circle is a radius and a frame; the
//! parameterization of each is defined *relative to* its frame, which is what
//! makes "the seam of this cylinder" a well-defined place rather than an
//! accident of how the surface was built.
//!
//! Unlike the conventional design, which splits right-handed and
//! possibly-left-handed frames into two separate types, there is one [`Frame`]
//! carrying a [`Handedness`]. The split buys nothing and costs a conversion at
//! every boundary between them.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};

use crate::{Direction, Direction2, Matrix3, Point, Point2, Vector, Vector2};

/// Whether a frame's third direction follows the right-hand rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Handedness {
    /// `x × y = z`. The usual case.
    #[default]
    Right,
    /// `x × y = -z`. Arises from mirroring.
    Left,
}

impl Handedness {
    /// `+1` for right-handed, `-1` for left-handed.
    #[must_use]
    pub const fn sign(self) -> f64 {
        match self {
            Self::Right => 1.0,
            Self::Left => -1.0,
        }
    }

    /// The opposite handedness.
    #[must_use]
    pub const fn flipped(self) -> Self {
        match self {
            Self::Right => Self::Left,
            Self::Left => Self::Right,
        }
    }
}

/// A point and a direction: an oriented line through space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Axis {
    /// A point on the axis.
    pub location: Point,
    /// The axis direction.
    pub direction: Direction,
}

/// A point and a direction in the plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Axis2 {
    /// A point on the axis.
    pub location: Point2,
    /// The axis direction.
    pub direction: Direction2,
}

/// A local coordinate system in space.
///
/// The `z` direction is primary — it is the axis of revolution for a cylinder,
/// the normal of a plane, the axis of a circle. `x` fixes where parameterization
/// starts. `y` is derived and always consistent with the handedness.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame {
    origin: Point,
    z: Direction,
    x: Direction,
    y: Direction,
    handedness: Handedness,
}

/// A local coordinate system in the plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Frame2 {
    origin: Point2,
    x: Direction2,
    y: Direction2,
    handedness: Handedness,
}

impl Axis {
    /// The X axis through the origin.
    pub const X: Self = Self {
        location: Point::ORIGIN,
        direction: Direction::X,
    };
    /// The Y axis through the origin.
    pub const Y: Self = Self {
        location: Point::ORIGIN,
        direction: Direction::Y,
    };
    /// The Z axis through the origin.
    pub const Z: Self = Self {
        location: Point::ORIGIN,
        direction: Direction::Z,
    };

    /// An axis from a point and a direction.
    #[must_use]
    pub const fn new(location: Point, direction: Direction) -> Self {
        Self {
            location,
            direction,
        }
    }

    /// The axis through two distinct points, directed from `from` to `to`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the points
    /// coincide.
    pub fn through(from: Point, to: Point, tol: Tolerances) -> OgeomResult<Self> {
        Ok(Self::new(from, Direction::new(to - from, tol)?))
    }

    /// This axis with its direction reversed.
    #[must_use]
    pub const fn reversed(self) -> Self {
        Self::new(self.location, self.direction.reversed())
    }

    /// The point at parameter `t`, measured from [`Axis::location`] in units of
    /// length along the direction.
    #[must_use]
    pub fn point_at(self, t: f64) -> Point {
        self.location + self.direction * t
    }

    /// The parameter of the projection of `p` onto this axis.
    #[must_use]
    pub fn parameter_of(self, p: Point) -> f64 {
        self.direction.dot_vector(p - self.location)
    }

    /// The closest point on this axis to `p`.
    #[must_use]
    pub fn project(self, p: Point) -> Point {
        self.point_at(self.parameter_of(p))
    }

    /// The perpendicular distance from `p` to this axis.
    #[must_use]
    pub fn distance_to(self, p: Point) -> f64 {
        // The cross product with a unit direction gives the perpendicular
        // component directly, without the cancellation that subtracting the
        // projection would introduce for a point far along the axis.
        self.direction.cross_with(p - self.location).magnitude()
    }

    /// Whether `p` lies on this axis within `tol.confusion()`.
    #[must_use]
    pub fn contains(self, p: Point, tol: Tolerances) -> bool {
        self.distance_to(p) <= tol.confusion()
    }

    /// Whether two axes are the same line with the same sense.
    #[must_use]
    pub fn is_coaxial(self, other: Self, tol: Tolerances) -> bool {
        self.direction.is_equal(other.direction, tol)
            && self.contains(other.location, tol)
            && other.contains(self.location, tol)
    }

    /// Whether two axes lie on the same line, ignoring sense.
    #[must_use]
    pub fn is_collinear(self, other: Self, tol: Tolerances) -> bool {
        self.direction.is_parallel(other.direction, tol)
            && self.contains(other.location, tol)
            && other.contains(self.location, tol)
    }
}

impl Axis2 {
    /// The X axis through the origin.
    pub const X: Self = Self {
        location: Point2::ORIGIN,
        direction: Direction2::X,
    };
    /// The Y axis through the origin.
    pub const Y: Self = Self {
        location: Point2::ORIGIN,
        direction: Direction2::Y,
    };

    /// An axis from a point and a direction.
    #[must_use]
    pub const fn new(location: Point2, direction: Direction2) -> Self {
        Self {
            location,
            direction,
        }
    }

    /// The axis through two distinct points.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the points
    /// coincide.
    pub fn through(from: Point2, to: Point2, tol: Tolerances) -> OgeomResult<Self> {
        Ok(Self::new(from, Direction2::new(to - from, tol)?))
    }

    /// This axis with its direction reversed.
    #[must_use]
    pub const fn reversed(self) -> Self {
        Self::new(self.location, self.direction.reversed())
    }

    /// The point at parameter `t`.
    #[must_use]
    pub fn point_at(self, t: f64) -> Point2 {
        self.location + self.direction * t
    }

    /// The parameter of the projection of `p` onto this axis.
    #[must_use]
    pub fn parameter_of(self, p: Point2) -> f64 {
        self.direction.vector().dot(p - self.location)
    }

    /// The closest point on this axis to `p`.
    #[must_use]
    pub fn project(self, p: Point2) -> Point2 {
        self.point_at(self.parameter_of(p))
    }

    /// The signed distance from `p` to this axis, positive on the left.
    #[must_use]
    pub fn signed_distance_to(self, p: Point2) -> f64 {
        self.direction.vector().cross(p - self.location)
    }

    /// The distance from `p` to this axis.
    #[must_use]
    pub fn distance_to(self, p: Point2) -> f64 {
        self.signed_distance_to(p).abs()
    }
}

impl Default for Frame {
    fn default() -> Self {
        Self::WORLD
    }
}

impl Frame {
    /// The identity frame: origin at the origin, axes along X, Y and Z.
    pub const WORLD: Self = Self {
        origin: Point::ORIGIN,
        z: Direction::Z,
        x: Direction::X,
        y: Direction::Y,
        handedness: Handedness::Right,
    };

    /// A right-handed frame from an origin, a primary direction and a reference
    /// for the first axis.
    ///
    /// `x_reference` need not be perpendicular to `z`: its component along `z`
    /// is removed. It must not be parallel to `z`, since then there is nothing
    /// left to orient by.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `z` and
    /// `x_reference` are parallel.
    pub fn new(
        origin: Point,
        z: Direction,
        x_reference: Direction,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        // Gram-Schmidt: remove the component of the reference along z, then
        // renormalize. Fails cleanly when nothing is left to normalize.
        let v = x_reference.vector() - z.vector() * z.dot(x_reference);
        let Ok(x) = Direction::new(v, tol) else {
            ogeom_bail!(
                Construction,
                "frame reference direction is parallel to the primary direction"
            );
        };
        let y = Direction::new(z.cross_vector(x), tol)?;
        Ok(Self {
            origin,
            z,
            x,
            y,
            handedness: Handedness::Right,
        })
    }

    /// A right-handed frame with an arbitrary but deterministic first axis.
    ///
    /// For geometry with rotational symmetry — a sphere, a full circle — the
    /// choice of `x` is immaterial, and requiring the caller to invent one is
    /// noise.
    #[must_use]
    pub fn about(origin: Point, z: Direction) -> Self {
        let x = z.any_perpendicular();
        // z and x are perpendicular unit vectors, so their cross product is
        // already unit length.
        let y =
            Direction::new(z.cross_vector(x), Tolerances::millimetres()).unwrap_or(Direction::Y);
        Self {
            origin,
            z,
            x,
            y,
            handedness: Handedness::Right,
        }
    }

    /// A frame from three directions given explicitly.
    ///
    /// Handedness is inferred from the triple product rather than asserted.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the three
    /// are not mutually perpendicular within `tol.angular()`, or if they are
    /// coplanar.
    pub fn from_axes(
        origin: Point,
        x: Direction,
        y: Direction,
        z: Direction,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        for (a, b, names) in [(x, y, "x/y"), (y, z, "y/z"), (z, x, "z/x")] {
            if !a.is_normal(b, tol) {
                ogeom_bail!(Construction, "frame axes {names} are not perpendicular");
            }
        }
        let triple = x.vector().triple(y.vector(), z.vector());
        if triple.abs() <= tol.angular() {
            ogeom_bail!(Construction, "frame axes are coplanar");
        }
        let handedness = if triple > 0.0 {
            Handedness::Right
        } else {
            Handedness::Left
        };
        Ok(Self {
            origin,
            z,
            x,
            y,
            handedness,
        })
    }

    /// The frame's origin.
    #[must_use]
    pub const fn origin(&self) -> Point {
        self.origin
    }

    /// The primary direction — the normal of a plane, the axis of a cylinder.
    #[must_use]
    pub const fn z(&self) -> Direction {
        self.z
    }

    /// The first axis, fixing where parameterization starts.
    #[must_use]
    pub const fn x(&self) -> Direction {
        self.x
    }

    /// The second axis.
    #[must_use]
    pub const fn y(&self) -> Direction {
        self.y
    }

    /// This frame's handedness.
    #[must_use]
    pub const fn handedness(&self) -> Handedness {
        self.handedness
    }

    /// The axis along the primary direction.
    #[must_use]
    pub const fn axis(&self) -> Axis {
        Axis::new(self.origin, self.z)
    }

    /// This frame moved to a new origin.
    #[must_use]
    pub const fn with_origin(&self, origin: Point) -> Self {
        Self { origin, ..*self }
    }

    /// This frame with its handedness flipped, by reversing `y`.
    #[must_use]
    pub const fn mirrored(&self) -> Self {
        Self {
            y: self.y.reversed(),
            handedness: self.handedness.flipped(),
            ..*self
        }
    }

    /// This frame with the primary direction reversed.
    ///
    /// `x` is kept, so `y` must flip to preserve handedness — reversing a
    /// plane's normal should not silently turn its parameterization inside out.
    #[must_use]
    pub const fn with_z_reversed(&self) -> Self {
        Self {
            z: self.z.reversed(),
            y: self.y.reversed(),
            ..*self
        }
    }

    /// Local coordinates of a point given in world coordinates.
    #[must_use]
    pub fn to_local(&self, p: Point) -> Point {
        let v = p - self.origin;
        Point::new(
            self.x.dot_vector(v),
            self.y.dot_vector(v),
            self.z.dot_vector(v),
        )
    }

    /// World coordinates of a point given in this frame's local coordinates.
    #[must_use]
    pub fn to_world(&self, p: Point) -> Point {
        self.origin + self.x * p.x + self.y * p.y + self.z * p.z
    }

    /// Local components of a world-space vector. Unaffected by the origin.
    #[must_use]
    pub fn vector_to_local(&self, v: Vector) -> Vector {
        Vector::new(
            self.x.dot_vector(v),
            self.y.dot_vector(v),
            self.z.dot_vector(v),
        )
    }

    /// World components of a vector given in local coordinates.
    #[must_use]
    pub fn vector_to_world(&self, v: Vector) -> Vector {
        self.x * v.x + self.y * v.y + self.z * v.z
    }

    /// The rotation taking local coordinates to world coordinates.
    #[must_use]
    pub fn to_matrix(&self) -> Matrix3 {
        Matrix3::from_columns(self.x.vector(), self.y.vector(), self.z.vector())
    }

    /// Whether two frames agree in origin and all three directions.
    #[must_use]
    pub fn is_equal(&self, other: &Self, tol: Tolerances) -> bool {
        self.origin.is_equal(other.origin, tol)
            && self.x.is_equal(other.x, tol)
            && self.y.is_equal(other.y, tol)
            && self.z.is_equal(other.z, tol)
    }

    /// The signed distance from `p` to this frame's XY plane, positive on the
    /// side the primary direction points to.
    #[must_use]
    pub fn signed_distance_to_plane(&self, p: Point) -> f64 {
        self.z.dot_vector(p - self.origin)
    }
}

impl Default for Frame2 {
    fn default() -> Self {
        Self::WORLD
    }
}

impl Frame2 {
    /// The identity frame.
    pub const WORLD: Self = Self {
        origin: Point2::ORIGIN,
        x: Direction2::X,
        y: Direction2::Y,
        handedness: Handedness::Right,
    };

    /// A right-handed frame from an origin and a first axis.
    #[must_use]
    pub const fn new(origin: Point2, x: Direction2) -> Self {
        Self {
            origin,
            x,
            y: x.perpendicular(),
            handedness: Handedness::Right,
        }
    }

    /// A frame with an explicit second axis, whose handedness is inferred.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the two
    /// axes are not perpendicular within `tol.angular()`.
    pub fn from_axes(
        origin: Point2,
        x: Direction2,
        y: Direction2,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        if !x.is_normal(y, tol) {
            ogeom_bail!(Construction, "frame axes are not perpendicular");
        }
        let handedness = if x.cross(y) > 0.0 {
            Handedness::Right
        } else {
            Handedness::Left
        };
        Ok(Self {
            origin,
            x,
            y,
            handedness,
        })
    }

    /// The frame's origin.
    #[must_use]
    pub const fn origin(&self) -> Point2 {
        self.origin
    }

    /// The first axis.
    #[must_use]
    pub const fn x(&self) -> Direction2 {
        self.x
    }

    /// The second axis.
    #[must_use]
    pub const fn y(&self) -> Direction2 {
        self.y
    }

    /// This frame's handedness.
    #[must_use]
    pub const fn handedness(&self) -> Handedness {
        self.handedness
    }

    /// This frame with its handedness flipped.
    #[must_use]
    pub const fn mirrored(&self) -> Self {
        Self {
            y: self.y.reversed(),
            handedness: self.handedness.flipped(),
            ..*self
        }
    }

    /// Local coordinates of a point given in world coordinates.
    #[must_use]
    pub fn to_local(&self, p: Point2) -> Point2 {
        let v = p - self.origin;
        Point2::new(self.x.vector().dot(v), self.y.vector().dot(v))
    }

    /// World coordinates of a point given in local coordinates.
    #[must_use]
    pub fn to_world(&self, p: Point2) -> Point2 {
        self.origin + self.x * p.x + self.y * p.y
    }

    /// Local components of a world-space vector.
    #[must_use]
    pub fn vector_to_local(&self, v: Vector2) -> Vector2 {
        Vector2::new(self.x.vector().dot(v), self.y.vector().dot(v))
    }

    /// World components of a vector given in local coordinates.
    #[must_use]
    pub fn vector_to_world(&self, v: Vector2) -> Vector2 {
        self.x * v.x + self.y * v.y
    }

    /// Whether two frames agree in origin and both directions.
    #[must_use]
    pub fn is_equal(&self, other: &Self, tol: Tolerances) -> bool {
        self.origin.is_equal(other.origin, tol)
            && self.x.is_equal(other.x, tol)
            && self.y.is_equal(other.y, tol)
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn axis_projection_and_distance() {
        let a = Axis::new(Point::new(1.0, 0.0, 0.0), Direction::Z);
        let p = Point::new(4.0, 0.0, 7.0);
        assert_relative_eq!(a.parameter_of(p), 7.0);
        assert_eq!(a.project(p), Point::new(1.0, 0.0, 7.0));
        assert_relative_eq!(a.distance_to(p), 3.0);
        assert!(a.contains(Point::new(1.0, 0.0, -5.0), T));
        assert!(!a.contains(p, T));
    }

    #[test]
    fn axis_distance_stays_accurate_far_along_the_axis() {
        // Subtracting the projection would cancel two numbers around 1e9 to
        // recover a distance of 3. The cross-product form does not.
        let a = Axis::Z;
        let p = Point::new(3.0, 0.0, 1.0e9);
        assert_relative_eq!(a.distance_to(p), 3.0, epsilon = 1e-9);
    }

    #[test]
    fn axis_through_coincident_points_is_refused() {
        let p = Point::new(1.0, 2.0, 3.0);
        assert!(Axis::through(p, p, T).is_err());
        assert!(Axis::through(p, Point::new(1.0, 2.0, 4.0), T).is_ok());
    }

    #[test]
    fn coaxial_and_collinear_differ_by_sense() {
        let a = Axis::Z;
        let b = Axis::new(Point::new(0.0, 0.0, 5.0), Direction::Z);
        let c = b.reversed();
        assert!(a.is_coaxial(b, T));
        assert!(!a.is_coaxial(c, T), "opposite sense is not coaxial");
        assert!(a.is_collinear(c, T), "but it is collinear");
        assert!(!a.is_collinear(Axis::X, T));
    }

    #[test]
    fn frame_orthonormalizes_a_non_perpendicular_reference() {
        // The reference leans heavily into z; only its perpendicular part
        // should survive.
        let reference = Direction::from_coords(1.0, 0.0, 10.0, T).unwrap();
        let f = Frame::new(Point::ORIGIN, Direction::Z, reference, T).unwrap();
        assert!(f.x().is_equal(Direction::X, T));
        assert!(f.y().is_equal(Direction::Y, T));
        assert!(f.to_matrix().is_orthonormal(1e-14));
    }

    #[test]
    fn frame_refuses_a_parallel_reference() {
        assert!(Frame::new(Point::ORIGIN, Direction::Z, Direction::Z, T).is_err());
        assert!(Frame::new(Point::ORIGIN, Direction::Z, -Direction::Z, T).is_err());
    }

    #[test]
    fn frame_about_works_for_every_primary_direction() {
        for z in [
            Direction::X,
            Direction::Y,
            Direction::Z,
            -Direction::Y,
            Direction::from_coords(1.0, 1.0, 1.0, T).unwrap(),
        ] {
            let f = Frame::about(Point::new(1.0, 2.0, 3.0), z);
            assert!(f.z().is_equal(z, T));
            assert!(f.to_matrix().is_orthonormal(1e-14));
            assert_eq!(f.handedness(), Handedness::Right);
        }
    }

    #[test]
    fn local_and_world_coordinates_round_trip() {
        let f = Frame::new(
            Point::new(10.0, -5.0, 2.0),
            Direction::from_coords(1.0, 1.0, 1.0, T).unwrap(),
            Direction::X,
            T,
        )
        .unwrap();
        for p in [
            Point::ORIGIN,
            Point::new(1.0, 2.0, 3.0),
            Point::new(-100.0, 0.5, 7.0),
        ] {
            assert!(f.to_world(f.to_local(p)).is_equal(p, T));
        }
        // The origin maps to local zero, and the axes to the unit vectors.
        assert!(f.to_local(f.origin()).is_equal(Point::ORIGIN, T));
        assert!(
            f.to_local(f.origin() + f.x() * 1.0)
                .is_equal(Point::new(1.0, 0.0, 0.0), T)
        );
    }

    #[test]
    fn vectors_ignore_the_origin_but_points_do_not() {
        let f = Frame::new(Point::new(100.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
        let v = Vector::new(1.0, 2.0, 3.0);
        assert!(
            f.vector_to_local(v).is_equal(v, T),
            "aligned frame, offset origin"
        );
        assert!(
            !f.to_local(Point::from_vector(v))
                .is_equal(Point::from_vector(v), T)
        );
    }

    #[test]
    fn handedness_is_inferred_not_asserted() {
        let right =
            Frame::from_axes(Point::ORIGIN, Direction::X, Direction::Y, Direction::Z, T).unwrap();
        assert_eq!(right.handedness(), Handedness::Right);

        let left =
            Frame::from_axes(Point::ORIGIN, Direction::X, Direction::Y, -Direction::Z, T).unwrap();
        assert_eq!(left.handedness(), Handedness::Left);
        assert_relative_eq!(left.handedness().sign(), -1.0);
    }

    #[test]
    fn from_axes_rejects_non_orthogonal_and_coplanar_input() {
        let skew = Direction::from_coords(1.0, 1.0, 0.0, T).unwrap();
        assert!(Frame::from_axes(Point::ORIGIN, Direction::X, skew, Direction::Z, T).is_err());
        assert!(
            Frame::from_axes(Point::ORIGIN, Direction::X, Direction::Y, Direction::X, T).is_err()
        );
    }

    #[test]
    fn reversing_the_primary_direction_preserves_handedness() {
        let f = Frame::WORLD;
        let r = f.with_z_reversed();
        assert!(r.z().is_equal(-Direction::Z, T));
        assert!(r.x().is_equal(Direction::X, T), "x is kept");
        assert!(r.y().is_equal(-Direction::Y, T), "y flips to compensate");
        assert_eq!(r.handedness(), Handedness::Right);
        assert_relative_eq!(
            r.x().vector().triple(r.y().vector(), r.z().vector()),
            1.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn mirroring_flips_handedness() {
        let m = Frame::WORLD.mirrored();
        assert_eq!(m.handedness(), Handedness::Left);
        assert_eq!(m.mirrored().handedness(), Handedness::Right);
        assert_relative_eq!(
            m.x().vector().triple(m.y().vector(), m.z().vector()),
            -1.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn signed_distance_to_the_frame_plane() {
        let f = Frame::WORLD;
        assert_relative_eq!(f.signed_distance_to_plane(Point::new(1.0, 2.0, 3.0)), 3.0);
        assert_relative_eq!(f.signed_distance_to_plane(Point::new(1.0, 2.0, -3.0)), -3.0);
        assert_relative_eq!(
            f.with_z_reversed()
                .signed_distance_to_plane(Point::new(0.0, 0.0, 3.0)),
            -3.0
        );
    }

    #[test]
    fn frame2_round_trips_and_infers_handedness() {
        let f = Frame2::new(Point2::new(3.0, 4.0), Direction2::from_angle(0.6));
        assert_eq!(f.handedness(), Handedness::Right);
        for p in [Point2::ORIGIN, Point2::new(-2.0, 7.0)] {
            assert!(f.to_world(f.to_local(p)).is_equal(p, T));
        }
        let left = Frame2::from_axes(Point2::ORIGIN, Direction2::X, -Direction2::Y, T).unwrap();
        assert_eq!(left.handedness(), Handedness::Left);
        assert!(Frame2::from_axes(Point2::ORIGIN, Direction2::X, Direction2::X, T).is_err());
    }

    #[test]
    fn axis2_signed_distance_is_positive_on_the_left() {
        let a = Axis2::X;
        assert_relative_eq!(a.signed_distance_to(Point2::new(5.0, 2.0)), 2.0);
        assert_relative_eq!(a.signed_distance_to(Point2::new(5.0, -2.0)), -2.0);
        assert_relative_eq!(a.distance_to(Point2::new(5.0, -2.0)), 2.0);
        assert_eq!(a.project(Point2::new(5.0, 2.0)), Point2::new(5.0, 0.0));
    }
}
