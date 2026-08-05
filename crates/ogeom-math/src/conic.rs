//! Conic sections: circle, ellipse, hyperbola, parabola.
//!
//! Each is a shape described by a [`Frame`] and one or two size parameters. The
//! frame is not decoration — it fixes the parameterization. A circle's angular
//! parameter is measured from its frame's `x` axis towards its `y` axis, so
//! "the point at 0" is a specific place that survives the shape being stored,
//! reloaded and transformed.
//!
//! Conics live in their frame's `xy` plane, with the frame's `z` as the normal.
//!
//! Evaluation and derivatives are in [`crate::elementary`]; this module holds
//! the descriptions and the queries that follow directly from them.
//!
//! A straight line needs no type of its own: it is exactly an [`Axis`](crate::Axis),
//! and the conventional design's separate line type carries no information the
//! axis does not.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};

use crate::{Frame, Frame2, Point, Point2, Transform, Transform2};

/// The complete elliptic integral of the second kind, `E(m)`, for `m` in
/// `[0, 1]`.
///
/// Computed by the arithmetic-geometric mean, which converges quadratically —
/// the number of correct digits doubles each iteration, so `f64` accuracy costs
/// about seven steps regardless of `m`.
///
/// `E(0) = pi/2` (a circle) and `E(1) = 1` (a degenerate ellipse, a doubled
/// line segment).
#[must_use]
pub fn complete_elliptic_e(m: f64) -> f64 {
    let m = m.clamp(0.0, 1.0);
    if m >= 1.0 {
        return 1.0;
    }
    let mut a = 1.0_f64;
    let mut b = (1.0 - m).sqrt();
    // The running sum of `2^(n-1) * c_n^2`, starting with the n = 0 term.
    let mut sum = m * 0.5;
    let mut power = 1.0_f64;
    // Quadratic convergence reaches the f64 floor well inside this bound; the
    // limit is a backstop, not the expected exit.
    for _ in 0..20 {
        let c = (a - b) * 0.5;
        let next_a = f64::midpoint(a, b);
        b = (a * b).sqrt();
        a = next_a;
        sum += power * c * c;
        power *= 2.0;
        if c.abs() <= f64::EPSILON * a {
            break;
        }
    }
    // K(m) = pi / (2 * AGM), and E(m) = K(m) * (1 - sum).
    core::f64::consts::FRAC_PI_2 / a * (1.0 - sum)
}

/// Reject a size parameter that cannot describe a real shape.
fn check_positive(name: &str, value: f64, tol: Tolerances) -> OgeomResult<()> {
    if !value.is_finite() || value <= tol.confusion() {
        ogeom_bail!(Construction, "{name} {value} must be finite and positive");
    }
    Ok(())
}

/// A circle in space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circle {
    frame: Frame,
    radius: f64,
}

/// An ellipse in space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipse {
    frame: Frame,
    major_radius: f64,
    minor_radius: f64,
}

/// A hyperbola in space.
///
/// Only the branch on the positive `x` side of its frame is described; the
/// other branch is the same hyperbola with the frame's `x` reversed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hyperbola {
    frame: Frame,
    major_radius: f64,
    minor_radius: f64,
}

/// A parabola in space, opening along its frame's `x` axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Parabola {
    frame: Frame,
    focal: f64,
}

/// A circle in the plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Circle2 {
    frame: Frame2,
    radius: f64,
}

/// An ellipse in the plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipse2 {
    frame: Frame2,
    major_radius: f64,
    minor_radius: f64,
}

/// A hyperbola in the plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hyperbola2 {
    frame: Frame2,
    major_radius: f64,
    minor_radius: f64,
}

/// A parabola in the plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Parabola2 {
    frame: Frame2,
    focal: f64,
}

impl Circle {
    /// A circle of `radius` in the `xy` plane of `frame`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `radius` is
    /// not finite and positive.
    pub fn new(frame: Frame, radius: f64, tol: Tolerances) -> OgeomResult<Self> {
        check_positive("circle radius", radius, tol)?;
        Ok(Self { frame, radius })
    }

    /// The circle through three points.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the points
    /// are collinear or any two coincide — no circle passes through them.
    pub fn through(a: Point, b: Point, c: Point, tol: Tolerances) -> OgeomResult<Self> {
        let (ab, ac) = (b - a, c - a);

        // `from_cross` doubles as the collinearity test, and it is the scale-free
        // one: it asks whether `|ab x ac|` is small *relative to* `|ab| |ac|`,
        // so three points a micron apart are judged by the same standard as
        // three a metre apart. Comparing the cross product against a length
        // tolerance instead would reject small triangles whose plane is
        // perfectly well determined.
        let z = crate::Direction::from_cross(ab, ac, tol)?;

        let normal = ab.cross(ac);
        let to_centre = (normal.cross(ab) * ac.square_magnitude()
            + ac.cross(normal) * ab.square_magnitude())
            / (2.0 * normal.square_magnitude());
        let centre = a + to_centre;

        // A circle smaller than the confusion tolerance is degenerate, and
        // `Direction::new` reports that.
        let x = crate::Direction::new(a - centre, tol)?;
        Self::new(Frame::new(centre, z, x, tol)?, to_centre.magnitude(), tol)
    }

    /// The frame positioning this circle.
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

    /// The circumference.
    #[must_use]
    pub fn length(&self) -> f64 {
        core::f64::consts::TAU * self.radius
    }

    /// The area enclosed.
    #[must_use]
    pub fn area(&self) -> f64 {
        core::f64::consts::PI * self.radius * self.radius
    }

    /// The shortest distance from `p` to the circle.
    ///
    /// Zero on the circle, and positive both inside and outside — this is the
    /// distance to the curve, not to the disc it bounds.
    #[must_use]
    pub fn distance_to(&self, p: Point) -> f64 {
        let local = self.frame.to_local(p);
        // Distance from a point to a circle, in the cylindrical coordinates of
        // the circle's own frame: radially out to the circle, then along the
        // axis.
        let radial = local.xy().to_vector().magnitude() - self.radius;
        radial.hypot(local.z)
    }

    /// Whether `p` lies on the circle within `tol.confusion()`.
    #[must_use]
    pub fn contains(&self, p: Point, tol: Tolerances) -> bool {
        self.distance_to(p) <= tol.confusion()
    }

    /// This circle moved by `t`.
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

impl Ellipse {
    /// An ellipse with the given radii in the `xy` plane of `frame`, the major
    /// radius along `x`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if either
    /// radius is not finite and positive, or if the minor radius exceeds the
    /// major.
    pub fn new(
        frame: Frame,
        major_radius: f64,
        minor_radius: f64,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        check_positive("major radius", major_radius, tol)?;
        check_positive("minor radius", minor_radius, tol)?;
        if minor_radius > major_radius + tol.confusion() {
            ogeom_bail!(
                Construction,
                "minor radius {minor_radius} exceeds major radius {major_radius}"
            );
        }
        Ok(Self {
            frame,
            major_radius,
            minor_radius,
        })
    }

    /// The frame positioning this ellipse.
    #[must_use]
    pub const fn frame(&self) -> Frame {
        self.frame
    }

    /// The centre.
    #[must_use]
    pub const fn centre(&self) -> Point {
        self.frame.origin()
    }

    /// The major radius.
    #[must_use]
    pub const fn major_radius(&self) -> f64 {
        self.major_radius
    }

    /// The minor radius.
    #[must_use]
    pub const fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// The distance from the centre to either focus.
    #[must_use]
    pub fn focal_distance(&self) -> f64 {
        // Written as a difference of squares factored into a product, which
        // keeps precision for a nearly circular ellipse where the two radii
        // almost cancel.
        ((self.major_radius - self.minor_radius) * (self.major_radius + self.minor_radius)).sqrt()
    }

    /// The eccentricity, in `[0, 1)`. Zero for a circle.
    #[must_use]
    pub fn eccentricity(&self) -> f64 {
        self.focal_distance() / self.major_radius
    }

    /// The two foci, on the `x` axis either side of the centre.
    #[must_use]
    pub fn foci(&self) -> (Point, Point) {
        let offset = self.frame.x() * self.focal_distance();
        (self.centre() + offset, self.centre() - offset)
    }

    /// The area enclosed.
    #[must_use]
    pub fn area(&self) -> f64 {
        core::f64::consts::PI * self.major_radius * self.minor_radius
    }

    /// The circumference.
    ///
    /// Exact to machine precision, via the complete elliptic integral of the
    /// second kind — see [`complete_elliptic_e`].
    ///
    /// Ramanujan's well-known approximation was the obvious alternative and is
    /// not good enough: it is excellent for a nearly circular ellipse but its
    /// relative error reaches `1.2e-5` by an axis ratio of 10:1. Perimeter
    /// feeds arc-length parameterization and measurement, where that is a
    /// visible error rather than a rounding detail.
    #[must_use]
    pub fn length(&self) -> f64 {
        // m = e^2 = 1 - (b/a)^2, written factored to avoid cancellation when
        // the ellipse is nearly circular.
        let ratio = self.minor_radius / self.major_radius;
        let m = (1.0 - ratio) * (1.0 + ratio);
        4.0 * self.major_radius * complete_elliptic_e(m)
    }

    /// This ellipse moved by `t`.
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

impl Hyperbola {
    /// A hyperbola with the given radii in the `xy` plane of `frame`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if either
    /// radius is not finite and positive. Unlike an ellipse, the minor radius
    /// may exceed the major.
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

    /// The frame positioning this hyperbola.
    #[must_use]
    pub const fn frame(&self) -> Frame {
        self.frame
    }

    /// The centre — the midpoint of the two vertices, not a point on the curve.
    #[must_use]
    pub const fn centre(&self) -> Point {
        self.frame.origin()
    }

    /// The major radius: centre to vertex.
    #[must_use]
    pub const fn major_radius(&self) -> f64 {
        self.major_radius
    }

    /// The minor radius, governing how fast the branches open.
    #[must_use]
    pub const fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// The distance from the centre to either focus.
    #[must_use]
    pub fn focal_distance(&self) -> f64 {
        self.major_radius.hypot(self.minor_radius)
    }

    /// The eccentricity, always greater than 1.
    #[must_use]
    pub fn eccentricity(&self) -> f64 {
        self.focal_distance() / self.major_radius
    }

    /// The two foci.
    #[must_use]
    pub fn foci(&self) -> (Point, Point) {
        let offset = self.frame.x() * self.focal_distance();
        (self.centre() + offset, self.centre() - offset)
    }

    /// The vertex of the described branch.
    #[must_use]
    pub fn vertex(&self) -> Point {
        self.centre() + self.frame.x() * self.major_radius
    }

    /// This hyperbola moved by `t`.
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

impl Parabola {
    /// A parabola with the given focal length in the `xy` plane of `frame`.
    ///
    /// The apex is at the frame origin and the curve opens along `+x`; the
    /// focus sits at distance `focal` from the apex along `+x`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `focal` is
    /// not finite and positive. A zero focal length degenerates to a ray.
    pub fn new(frame: Frame, focal: f64, tol: Tolerances) -> OgeomResult<Self> {
        check_positive("focal length", focal, tol)?;
        Ok(Self { frame, focal })
    }

    /// The frame positioning this parabola.
    #[must_use]
    pub const fn frame(&self) -> Frame {
        self.frame
    }

    /// The apex.
    #[must_use]
    pub const fn apex(&self) -> Point {
        self.frame.origin()
    }

    /// The focal length: apex to focus.
    #[must_use]
    pub const fn focal(&self) -> f64 {
        self.focal
    }

    /// The focus.
    #[must_use]
    pub fn focus(&self) -> Point {
        self.apex() + self.frame.x() * self.focal
    }

    /// The eccentricity, which is `1` for every parabola.
    #[must_use]
    pub const fn eccentricity(&self) -> f64 {
        1.0
    }

    /// This parabola moved by `t`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the
    /// transformed frame is degenerate.
    pub fn transformed(&self, t: &Transform, tol: Tolerances) -> OgeomResult<Self> {
        Self::new(
            t.apply_frame(&self.frame, tol)?,
            self.focal * t.scale_factor().abs(),
            tol,
        )
    }
}

impl Circle2 {
    /// A circle of `radius` centred on `frame`'s origin.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `radius` is
    /// not finite and positive.
    pub fn new(frame: Frame2, radius: f64, tol: Tolerances) -> OgeomResult<Self> {
        check_positive("circle radius", radius, tol)?;
        Ok(Self { frame, radius })
    }

    /// A circle from a centre and a radius, with the default orientation.
    ///
    /// # Errors
    ///
    /// As [`Circle2::new`].
    pub fn centred(centre: Point2, radius: f64, tol: Tolerances) -> OgeomResult<Self> {
        Self::new(Frame2::new(centre, crate::Direction2::X), radius, tol)
    }

    /// The frame positioning this circle.
    #[must_use]
    pub const fn frame(&self) -> Frame2 {
        self.frame
    }

    /// The centre.
    #[must_use]
    pub const fn centre(&self) -> Point2 {
        self.frame.origin()
    }

    /// The radius.
    #[must_use]
    pub const fn radius(&self) -> f64 {
        self.radius
    }

    /// The circumference.
    #[must_use]
    pub fn length(&self) -> f64 {
        core::f64::consts::TAU * self.radius
    }

    /// The area enclosed.
    #[must_use]
    pub fn area(&self) -> f64 {
        core::f64::consts::PI * self.radius * self.radius
    }

    /// The signed distance from `p` to the circle, negative inside.
    #[must_use]
    pub fn signed_distance_to(&self, p: Point2) -> f64 {
        self.centre().distance(p) - self.radius
    }

    /// Whether `p` lies on the circle within `tol.confusion()`.
    #[must_use]
    pub fn contains(&self, p: Point2, tol: Tolerances) -> bool {
        self.signed_distance_to(p).abs() <= tol.confusion()
    }

    /// Whether `p` lies strictly inside.
    #[must_use]
    pub fn encloses(&self, p: Point2, tol: Tolerances) -> bool {
        self.signed_distance_to(p) < -tol.confusion()
    }

    /// This circle moved by `t`.
    ///
    /// The frame goes through the transform intact rather than being rebuilt
    /// from the centre. The frame fixes where the angular parameter starts, so
    /// discarding its orientation would keep the shape and silently renumber
    /// every point on it — invisible to a distance check, and wrong for
    /// anything holding a parameter.
    ///
    /// # Errors
    ///
    /// As [`Circle2::new`].
    pub fn transformed(&self, t: &Transform2, tol: Tolerances) -> OgeomResult<Self> {
        Self::new(
            t.apply_frame(&self.frame, tol)?,
            self.radius * t.scale_factor().abs(),
            tol,
        )
    }
}

impl Ellipse2 {
    /// An ellipse with the given radii, the major radius along `frame`'s `x`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a radius is
    /// not finite and positive, or the minor exceeds the major.
    pub fn new(
        frame: Frame2,
        major_radius: f64,
        minor_radius: f64,
        tol: Tolerances,
    ) -> OgeomResult<Self> {
        check_positive("major radius", major_radius, tol)?;
        check_positive("minor radius", minor_radius, tol)?;
        if minor_radius > major_radius + tol.confusion() {
            ogeom_bail!(
                Construction,
                "minor radius {minor_radius} exceeds major radius {major_radius}"
            );
        }
        Ok(Self {
            frame,
            major_radius,
            minor_radius,
        })
    }

    /// The frame positioning this ellipse.
    #[must_use]
    pub const fn frame(&self) -> Frame2 {
        self.frame
    }

    /// The centre.
    #[must_use]
    pub const fn centre(&self) -> Point2 {
        self.frame.origin()
    }

    /// The major radius.
    #[must_use]
    pub const fn major_radius(&self) -> f64 {
        self.major_radius
    }

    /// The minor radius.
    #[must_use]
    pub const fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// The distance from the centre to either focus.
    #[must_use]
    pub fn focal_distance(&self) -> f64 {
        ((self.major_radius - self.minor_radius) * (self.major_radius + self.minor_radius)).sqrt()
    }

    /// The eccentricity, in `[0, 1)`.
    #[must_use]
    pub fn eccentricity(&self) -> f64 {
        self.focal_distance() / self.major_radius
    }

    /// The area enclosed.
    #[must_use]
    pub fn area(&self) -> f64 {
        core::f64::consts::PI * self.major_radius * self.minor_radius
    }

    /// This ellipse moved by `t`.
    ///
    /// # Errors
    ///
    /// As [`Ellipse2::new`].
    pub fn transformed(&self, t: &Transform2, tol: Tolerances) -> OgeomResult<Self> {
        let s = t.scale_factor().abs();
        Self::new(
            t.apply_frame(&self.frame, tol)?,
            self.major_radius * s,
            self.minor_radius * s,
            tol,
        )
    }
}

impl Hyperbola2 {
    /// A hyperbola with the given radii.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a radius is
    /// not finite and positive.
    pub fn new(
        frame: Frame2,
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

    /// The frame positioning this hyperbola.
    #[must_use]
    pub const fn frame(&self) -> Frame2 {
        self.frame
    }

    /// The centre.
    #[must_use]
    pub const fn centre(&self) -> Point2 {
        self.frame.origin()
    }

    /// The major radius.
    #[must_use]
    pub const fn major_radius(&self) -> f64 {
        self.major_radius
    }

    /// The minor radius.
    #[must_use]
    pub const fn minor_radius(&self) -> f64 {
        self.minor_radius
    }

    /// The eccentricity, always greater than 1.
    #[must_use]
    pub fn eccentricity(&self) -> f64 {
        self.major_radius.hypot(self.minor_radius) / self.major_radius
    }
}

impl Parabola2 {
    /// A parabola with the given focal length, opening along `frame`'s `x`.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `focal` is
    /// not finite and positive.
    pub fn new(frame: Frame2, focal: f64, tol: Tolerances) -> OgeomResult<Self> {
        check_positive("focal length", focal, tol)?;
        Ok(Self { frame, focal })
    }

    /// The frame positioning this parabola.
    #[must_use]
    pub const fn frame(&self) -> Frame2 {
        self.frame
    }

    /// The apex.
    #[must_use]
    pub const fn apex(&self) -> Point2 {
        self.frame.origin()
    }

    /// The focal length.
    #[must_use]
    pub const fn focal(&self) -> f64 {
        self.focal
    }

    /// The focus.
    #[must_use]
    pub fn focus(&self) -> Point2 {
        self.apex() + self.frame.x() * self.focal
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::{Axis, Direction, Vector};
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();

    #[test]
    fn degenerate_sizes_are_refused() {
        let f = Frame::WORLD;
        assert!(Circle::new(f, 0.0, T).is_err());
        assert!(Circle::new(f, -1.0, T).is_err());
        assert!(Circle::new(f, f64::NAN, T).is_err());
        assert!(Circle::new(f, f64::INFINITY, T).is_err());
        assert!(Parabola::new(f, 0.0, T).is_err());
        // An ellipse whose minor radius exceeds its major is not an ellipse with
        // the axes swapped, it is a construction error.
        assert!(Ellipse::new(f, 1.0, 2.0, T).is_err());
        assert!(Ellipse::new(f, 2.0, 1.0, T).is_ok());
        // A hyperbola has no such constraint.
        assert!(Hyperbola::new(f, 1.0, 2.0, T).is_ok());
    }

    #[test]
    fn circle_measurements() {
        let c = Circle::new(Frame::WORLD, 2.0, T).unwrap();
        assert_relative_eq!(c.length(), core::f64::consts::TAU * 2.0);
        assert_relative_eq!(c.area(), core::f64::consts::PI * 4.0);
        assert_eq!(c.centre(), Point::ORIGIN);
    }

    #[test]
    fn circle_distance_is_to_the_curve_not_the_disc() {
        let c = Circle::new(Frame::WORLD, 5.0, T).unwrap();
        // The centre is 5 from the circle, not 0.
        assert_relative_eq!(c.distance_to(Point::ORIGIN), 5.0);
        assert_relative_eq!(c.distance_to(Point::new(5.0, 0.0, 0.0)), 0.0);
        assert_relative_eq!(c.distance_to(Point::new(7.0, 0.0, 0.0)), 2.0);
        assert_relative_eq!(c.distance_to(Point::new(3.0, 0.0, 0.0)), 2.0);
        // Off the plane, the distance combines radial and axial parts.
        assert_relative_eq!(c.distance_to(Point::new(5.0, 0.0, 3.0)), 3.0);
        assert_relative_eq!(c.distance_to(Point::new(1.0, 0.0, 3.0)), 5.0);
        assert!(c.contains(Point::new(0.0, 5.0, 0.0), T));
    }

    #[test]
    fn circle_through_three_points() {
        // Three points on the unit circle in the xy plane.
        let c = Circle::through(
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(-1.0, 0.0, 0.0),
            T,
        )
        .unwrap();
        assert_relative_eq!(c.radius(), 1.0, epsilon = 1e-12);
        assert!(c.centre().is_equal(Point::ORIGIN, T));
        assert!(c.contains(Point::new(0.0, -1.0, 0.0), T));
    }

    #[test]
    fn circle_through_three_points_works_off_axis_and_at_scale() {
        for scale in [1e-3_f64, 1.0, 1e3] {
            let a = Point::new(3.0 * scale, scale, 2.0 * scale);
            let b = Point::new(scale, 5.0 * scale, -scale);
            let c = Point::new(-2.0 * scale, 0.0, 4.0 * scale);
            let circle = Circle::through(a, b, c, T).unwrap();
            // Defining property: all three are on it, equidistant from centre.
            for p in [a, b, c] {
                assert_relative_eq!(
                    circle.centre().distance(p),
                    circle.radius(),
                    max_relative = 1e-10
                );
            }
        }
    }

    #[test]
    fn collinear_points_admit_no_circle() {
        let a = Point::ORIGIN;
        let b = Point::new(1.0, 1.0, 1.0);
        let c = Point::new(2.0, 2.0, 2.0);
        assert!(Circle::through(a, b, c, T).is_err());
        assert!(Circle::through(a, b, b, T).is_err());
        assert!(Circle::through(a, a, a, T).is_err());
    }

    #[test]
    fn collinearity_test_is_scale_free() {
        // Three points a micron apart are not collinear just because they are
        // close together. An absolute threshold on the cross product would
        // wrongly reject this.
        let s = 1e-6;
        let circle = Circle::through(
            Point::new(s, 0.0, 0.0),
            Point::new(0.0, s, 0.0),
            Point::new(-s, 0.0, 0.0),
            T,
        )
        .unwrap();
        assert_relative_eq!(circle.radius(), s, max_relative = 1e-9);
    }

    #[test]
    fn ellipse_focal_geometry() {
        let e = Ellipse::new(Frame::WORLD, 5.0, 3.0, T).unwrap();
        assert_relative_eq!(e.focal_distance(), 4.0, epsilon = 1e-12);
        assert_relative_eq!(e.eccentricity(), 0.8, epsilon = 1e-12);
        let (f1, f2) = e.foci();
        assert!(f1.is_equal(Point::new(4.0, 0.0, 0.0), T));
        assert!(f2.is_equal(Point::new(-4.0, 0.0, 0.0), T));
        // The defining property: distances to the two foci sum to 2a.
        let on_curve = Point::new(0.0, 3.0, 0.0);
        assert_relative_eq!(
            on_curve.distance(f1) + on_curve.distance(f2),
            10.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn nearly_circular_ellipse_keeps_focal_precision() {
        // a^2 - b^2 with a and b nearly equal cancels catastrophically; the
        // factored form does not.
        let a = 1.0;
        let b = 1.0 - 1e-12;
        let e = Ellipse::new(Frame::WORLD, a, b, T).unwrap();
        let expected = ((a - b) * (a + b)).sqrt();
        assert_relative_eq!(e.focal_distance(), expected, max_relative = 1e-12);
        assert!(e.focal_distance() > 0.0, "must not collapse to zero");
    }

    #[test]
    fn ellipse_length_matches_the_circle_it_degenerates_to() {
        let r = 3.0;
        let e = Ellipse::new(Frame::WORLD, r, r, T).unwrap();
        assert_relative_eq!(e.length(), core::f64::consts::TAU * r, max_relative = 1e-12);
        assert_relative_eq!(e.area(), core::f64::consts::PI * r * r);
        assert_relative_eq!(e.eccentricity(), 0.0);
    }

    #[test]
    fn ellipse_length_is_accurate_at_every_eccentricity() {
        // Reference values from numerical quadrature of the arc-length
        // integral, independent of the implementation under test.
        let cases = [
            (10.0, 1.0, 40.639_741_801_0),
            (2.0, 1.0, 9.688_448_220_5),
            (1.0, 1.0, core::f64::consts::TAU),
            (100.0, 1.0, 400.109_832_972_2),
        ];
        for (a, b, expected) in cases {
            let e = Ellipse::new(Frame::WORLD, a, b, T).unwrap();
            // The references come from quadrature, which itself carries error,
            // so the bound reflects the reference rather than the method.
            assert_relative_eq!(e.length(), expected, max_relative = 1e-9);
        }
    }

    #[test]
    fn complete_elliptic_integral_endpoints() {
        assert_relative_eq!(complete_elliptic_e(0.0), core::f64::consts::FRAC_PI_2);
        assert_relative_eq!(complete_elliptic_e(1.0), 1.0);
        // Monotonically decreasing in m.
        let mut previous = complete_elliptic_e(0.0);
        for i in 1..=20 {
            let e = complete_elliptic_e(f64::from(i) / 20.0);
            assert!(
                e < previous,
                "not decreasing at m = {}",
                f64::from(i) / 20.0
            );
            previous = e;
        }
    }

    #[test]
    fn hyperbola_focal_geometry() {
        let h = Hyperbola::new(Frame::WORLD, 3.0, 4.0, T).unwrap();
        assert_relative_eq!(h.focal_distance(), 5.0, epsilon = 1e-12);
        assert_relative_eq!(h.eccentricity(), 5.0 / 3.0, epsilon = 1e-12);
        assert!(h.eccentricity() > 1.0);
        assert!(h.vertex().is_equal(Point::new(3.0, 0.0, 0.0), T));
    }

    #[test]
    fn parabola_focal_geometry() {
        let p = Parabola::new(Frame::WORLD, 2.0, T).unwrap();
        assert!(p.apex().is_equal(Point::ORIGIN, T));
        assert!(p.focus().is_equal(Point::new(2.0, 0.0, 0.0), T));
        assert_relative_eq!(p.eccentricity(), 1.0);
    }

    #[test]
    fn transforms_scale_sizes_and_move_frames() {
        let c = Circle::new(Frame::WORLD, 2.0, T).unwrap();
        let t = Transform::scaling(Point::ORIGIN, 3.0, T).unwrap()
            * Transform::translation(Vector::new(1.0, 0.0, 0.0));
        let moved = c.transformed(&t, T).unwrap();
        assert_relative_eq!(moved.radius(), 6.0, epsilon = 1e-12);
        assert!(moved.centre().is_equal(Point::new(3.0, 0.0, 0.0), T));

        // A rotation leaves the radius alone but moves the frame.
        let r = Transform::rotation(Axis::X, core::f64::consts::FRAC_PI_2);
        let rotated = c.transformed(&r, T).unwrap();
        assert_relative_eq!(rotated.radius(), 2.0, epsilon = 1e-12);
        assert!(rotated.frame().z().is_equal(-Direction::Y, T));
    }

    #[test]
    fn mirroring_a_circle_keeps_a_positive_radius() {
        // A negative scale factor must not produce a negative radius; the shape
        // is mirrored through its frame, not inverted.
        let c = Circle::new(Frame::WORLD, 2.0, T).unwrap();
        let m = Transform::point_mirror(Point::ORIGIN);
        let mirrored = c.transformed(&m, T).unwrap();
        assert_relative_eq!(mirrored.radius(), 2.0);
    }

    #[test]
    fn circle2_signed_distance_and_containment() {
        let c = Circle2::centred(Point2::new(1.0, 1.0), 2.0, T).unwrap();
        assert_relative_eq!(c.signed_distance_to(Point2::new(1.0, 1.0)), -2.0);
        assert_relative_eq!(c.signed_distance_to(Point2::new(3.0, 1.0)), 0.0);
        assert_relative_eq!(c.signed_distance_to(Point2::new(5.0, 1.0)), 2.0);
        assert!(c.encloses(Point2::new(1.0, 1.0), T));
        assert!(
            !c.encloses(Point2::new(3.0, 1.0), T),
            "on the boundary is not inside"
        );
        assert!(c.contains(Point2::new(3.0, 1.0), T));
    }

    #[test]
    fn a_planar_circles_frame_survives_a_transform() {
        // The frame fixes where the angular parameter starts. Rebuilding it
        // from the centre would keep the shape and renumber every point on it,
        // which is invisible to a distance check and wrong for anything that
        // refers to a parameter.
        let f = Frame2::new(Point2::new(1.0, 2.0), crate::Direction2::from_angle(0.9));
        let c = Circle2::new(f, 2.0, T).unwrap();
        let rot = Transform2::rotation(Point2::ORIGIN, 0.5);
        let moved = c.transformed(&rot, T).unwrap();

        assert!(moved.centre().is_equal(rot.apply(c.centre()), T));
        let expected = rot.apply_direction(f.x(), T).unwrap();
        assert!(moved.frame().x().is_equal(expected, T));

        // The point at angle zero moves with the transform, rather than jumping
        // to wherever a rebuilt frame would have put it.
        let at_zero = c.centre() + f.x() * 2.0;
        let moved_at_zero = moved.centre() + moved.frame().x() * 2.0;
        assert!(moved_at_zero.is_equal(rot.apply(at_zero), T));
    }

    #[test]
    fn planar_conics_mirror_their_spatial_counterparts() {
        let e = Ellipse2::new(Frame2::WORLD, 5.0, 3.0, T).unwrap();
        assert_relative_eq!(e.focal_distance(), 4.0, epsilon = 1e-12);
        assert_relative_eq!(e.eccentricity(), 0.8, epsilon = 1e-12);
        assert_relative_eq!(e.area(), core::f64::consts::PI * 15.0);

        let h = Hyperbola2::new(Frame2::WORLD, 3.0, 4.0, T).unwrap();
        assert_relative_eq!(h.eccentricity(), 5.0 / 3.0, epsilon = 1e-12);

        let p = Parabola2::new(Frame2::WORLD, 2.0, T).unwrap();
        assert!(p.focus().is_equal(Point2::new(2.0, 0.0), T));

        assert!(Ellipse2::new(Frame2::WORLD, 1.0, 2.0, T).is_err());
    }
}
