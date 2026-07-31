//! Parameterization of the analytic primitives.
//!
//! Evaluation, derivatives and parameter inversion for every shape in
//! [`crate::conic`] and [`crate::quadric`]. Each shape's parameterization is
//! defined relative to its own [`Frame`](crate::Frame), which is what makes a
//! parameter value mean the same place across a save, a reload and a transform.
//!
//! # Conventions
//!
//! - **Circle, ellipse:** angle from the frame's `x` axis towards `y`.
//! - **Hyperbola:** `(a cosh t, b sinh t)`, describing the `+x` branch.
//! - **Parabola:** `(t^2 / (4 f), t)` with the apex at `t = 0`.
//! - **Plane:** `(u, v)` are the coordinates along `x` and `y`.
//! - **Cylinder:** `(angle, height)`.
//! - **Cone:** `(angle, height)`, with the radius varying along the height.
//! - **Sphere:** `(longitude, latitude)`, latitude in `[-pi/2, pi/2]`.
//! - **Torus:** `(angle about the axis, angle around the tube)`.
//!
//! Every surface here follows the same `(u, v)` order, `u` going around and `v`
//! going along. Getting that consistent matters: a caller that has to remember
//! which surface reverses the convention will eventually forget.

use og_core::{OgResult, Tolerances, og_bail};

use crate::{
    Axis, Circle, Cone, Cylinder, Direction, Ellipse, Hyperbola, Parabola, Plane, Point, Sphere,
    Torus, Vector,
};

/// A point on a curve with its first two derivatives.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurvePoint {
    /// The position.
    pub point: Point,
    /// First derivative with respect to the parameter.
    pub d1: Vector,
    /// Second derivative.
    pub d2: Vector,
}

impl CurvePoint {
    /// The unit tangent.
    ///
    /// # Errors
    ///
    /// [`OgError::Construction`](og_core::OgError::Construction) at a cusp,
    /// where the first derivative vanishes and there is no tangent direction.
    pub fn tangent(&self, tol: Tolerances) -> OgResult<Direction> {
        Direction::new(self.d1, tol)
    }

    /// The curvature.
    ///
    /// `|d1 x d2| / |d1|^3`. Zero on a straight section, and the reciprocal of
    /// the radius on a circle.
    #[must_use]
    pub fn curvature(&self) -> f64 {
        let speed = self.d1.magnitude();
        if speed == 0.0 {
            return 0.0;
        }
        self.d1.cross(self.d2).magnitude() / (speed * speed * speed)
    }
}

/// A point on a surface with its first derivatives.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfacePoint {
    /// The position.
    pub point: Point,
    /// Derivative along `u`.
    pub du: Vector,
    /// Derivative along `v`.
    pub dv: Vector,
}

impl SurfacePoint {
    /// The unit normal, `du x dv` normalized.
    ///
    /// # Errors
    ///
    /// [`OgError::Construction`](og_core::OgError::Construction) at a
    /// degeneracy — a sphere's pole or a cone's apex — where no normal is
    /// determined by the tangents.
    pub fn normal(&self, tol: Tolerances) -> OgResult<Direction> {
        if self.is_degenerate(tol) {
            og_bail!(
                Construction,
                "surface point is degenerate; the tangents determine no normal"
            );
        }
        Direction::new(self.du.cross(self.dv), tol)
    }

    /// Whether the tangents fail to determine a normal.
    ///
    /// Compared against the *square of the larger tangent*, not against the
    /// product of the two. The product test asks whether the tangents are
    /// collinear, which is only one of the two ways a surface degenerates: at a
    /// sphere's pole or a cone's apex one tangent *vanishes*, and there the
    /// product is itself near zero, so a relative-to-product test finds the
    /// cross product respectably large by comparison and declares the point
    /// healthy. Squaring the larger tangent keeps a fixed scale to judge
    /// against and catches both cases.
    #[must_use]
    pub fn is_degenerate(&self, tol: Tolerances) -> bool {
        let scale = self.du.magnitude().max(self.dv.magnitude());
        self.du.cross(self.dv).magnitude() <= tol.angular() * scale * scale
    }
}

/// Wrap `angle` into `[0, 2*pi)`.
#[must_use]
pub fn wrap_angle(angle: f64) -> f64 {
    let wrapped = angle.rem_euclid(core::f64::consts::TAU);
    // `rem_euclid` can return exactly TAU for a tiny negative input, which
    // would put a supposedly-normalized angle outside its own range.
    if wrapped >= core::f64::consts::TAU {
        0.0
    } else {
        wrapped
    }
}

/// Wrap `angle` into `(-pi, pi]`.
#[must_use]
pub fn wrap_signed_angle(angle: f64) -> f64 {
    let wrapped = wrap_angle(angle);
    if wrapped > core::f64::consts::PI {
        wrapped - core::f64::consts::TAU
    } else {
        wrapped
    }
}

/// Evaluate a line at `t`, measured in length along its direction.
#[must_use]
pub fn line_at(axis: Axis, t: f64) -> CurvePoint {
    CurvePoint {
        point: axis.point_at(t),
        d1: axis.direction.vector(),
        d2: Vector::ZERO,
    }
}

/// The parameter of the projection of `p` onto a line.
#[must_use]
pub fn line_parameter(axis: Axis, p: Point) -> f64 {
    axis.parameter_of(p)
}

/// Evaluate a circle at `angle`.
#[must_use]
pub fn circle_at(circle: &Circle, angle: f64) -> CurvePoint {
    let f = circle.frame();
    let r = circle.radius();
    let (sin, cos) = angle.sin_cos();
    let (x, y) = (f.x().vector(), f.y().vector());
    CurvePoint {
        point: circle.centre() + x * (r * cos) + y * (r * sin),
        d1: x * (-r * sin) + y * (r * cos),
        d2: x * (-r * cos) + y * (-r * sin),
    }
}

/// The angle of the point on a circle nearest `p`, in `[0, 2*pi)`.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if `p` lies on the
/// circle's axis, where every angle is equally near.
pub fn circle_parameter(circle: &Circle, p: Point, tol: Tolerances) -> OgResult<f64> {
    let local = circle.frame().to_local(p);
    if local.x.hypot(local.y) <= tol.confusion() {
        og_bail!(
            Construction,
            "point is on the circle's axis; no nearest angle"
        );
    }
    Ok(wrap_angle(local.y.atan2(local.x)))
}

/// Evaluate an ellipse at `angle`.
///
/// The parameter is the *eccentric* angle, not the polar one: the point is
/// `(a cos t, b sin t)`. That keeps evaluation free of trigonometric inversion
/// and matches every exchange format.
#[must_use]
pub fn ellipse_at(ellipse: &Ellipse, angle: f64) -> CurvePoint {
    let f = ellipse.frame();
    let (a, b) = (ellipse.major_radius(), ellipse.minor_radius());
    let (sin, cos) = angle.sin_cos();
    let (x, y) = (f.x().vector(), f.y().vector());
    CurvePoint {
        point: ellipse.centre() + x * (a * cos) + y * (b * sin),
        d1: x * (-a * sin) + y * (b * cos),
        d2: x * (-a * cos) + y * (-b * sin),
    }
}

/// The eccentric angle of the point on an ellipse in the direction of `p`.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if `p` projects to
/// the centre.
pub fn ellipse_parameter(ellipse: &Ellipse, p: Point, tol: Tolerances) -> OgResult<f64> {
    let local = ellipse.frame().to_local(p);
    if local.x.hypot(local.y) <= tol.confusion() {
        og_bail!(Construction, "point projects to the ellipse's centre");
    }
    // Undo the axis scaling before taking the angle, or the result is the polar
    // angle rather than the eccentric one.
    Ok(wrap_angle(
        (local.y / ellipse.minor_radius()).atan2(local.x / ellipse.major_radius()),
    ))
}

/// Evaluate a hyperbola at `t`, on the `+x` branch.
#[must_use]
pub fn hyperbola_at(hyperbola: &Hyperbola, t: f64) -> CurvePoint {
    let f = hyperbola.frame();
    let (a, b) = (hyperbola.major_radius(), hyperbola.minor_radius());
    let (cosh, sinh) = (t.cosh(), t.sinh());
    let (x, y) = (f.x().vector(), f.y().vector());
    CurvePoint {
        point: hyperbola.centre() + x * (a * cosh) + y * (b * sinh),
        d1: x * (a * sinh) + y * (b * cosh),
        // The second derivative of cosh is cosh, and of sinh is sinh, so this
        // is the position vector again — a hyperbola's acceleration points
        // away from its centre.
        d2: x * (a * cosh) + y * (b * sinh),
    }
}

/// The parameter of the point on a hyperbola in the direction of `p`.
///
/// # Errors
///
/// [`OgError::Domain`](og_core::OgError::Domain) if `p` projects to the far
/// branch, which this hyperbola does not describe.
pub fn hyperbola_parameter(hyperbola: &Hyperbola, p: Point, tol: Tolerances) -> OgResult<f64> {
    let local = hyperbola.frame().to_local(p);
    let a = local.x / hyperbola.major_radius();
    if a < 1.0 - tol.parametric() {
        og_bail!(
            Domain,
            "point is not on the described branch of the hyperbola"
        );
    }
    // asinh rather than acosh: acosh loses precision near t = 0, where its
    // argument approaches 1 and its derivative is unbounded.
    Ok((local.y / hyperbola.minor_radius()).asinh())
}

/// Evaluate a parabola at `t`.
#[must_use]
pub fn parabola_at(parabola: &Parabola, t: f64) -> CurvePoint {
    let f = parabola.frame();
    let focal = parabola.focal();
    let (x, y) = (f.x().vector(), f.y().vector());
    let scale = 1.0 / (4.0 * focal);
    CurvePoint {
        point: parabola.apex() + x * (t * t * scale) + y * t,
        d1: x * (2.0 * t * scale) + y,
        d2: x * (2.0 * scale),
    }
}

/// The parameter of the point on a parabola level with `p`.
#[must_use]
pub fn parabola_parameter(parabola: &Parabola, p: Point) -> f64 {
    parabola.frame().to_local(p).y
}

/// Evaluate a plane at `(u, v)`, the coordinates along its `x` and `y` axes.
#[must_use]
pub fn plane_at(plane: &Plane, u: f64, v: f64) -> SurfacePoint {
    let f = plane.frame();
    SurfacePoint {
        point: f.origin() + f.x() * u + f.y() * v,
        du: f.x().vector(),
        dv: f.y().vector(),
    }
}

/// The `(u, v)` of the projection of `p` onto a plane.
#[must_use]
pub fn plane_parameters(plane: &Plane, p: Point) -> (f64, f64) {
    let local = plane.frame().to_local(p);
    (local.x, local.y)
}

/// Evaluate a cylinder at `(angle, height)`.
#[must_use]
pub fn cylinder_at(cylinder: &Cylinder, angle: f64, height: f64) -> SurfacePoint {
    let f = cylinder.frame();
    let r = cylinder.radius();
    let (sin, cos) = angle.sin_cos();
    let (x, y, z) = (f.x().vector(), f.y().vector(), f.z().vector());
    SurfacePoint {
        point: f.origin() + x * (r * cos) + y * (r * sin) + z * height,
        du: x * (-r * sin) + y * (r * cos),
        dv: z,
    }
}

/// The `(angle, height)` of the point on a cylinder nearest `p`.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if `p` lies on the
/// axis.
pub fn cylinder_parameters(cylinder: &Cylinder, p: Point, tol: Tolerances) -> OgResult<(f64, f64)> {
    let local = cylinder.frame().to_local(p);
    if local.x.hypot(local.y) <= tol.confusion() {
        og_bail!(
            Construction,
            "point is on the cylinder's axis; no nearest angle"
        );
    }
    Ok((wrap_angle(local.y.atan2(local.x)), local.z))
}

/// Evaluate a cone at `(angle, height)`.
///
/// At the apex the two tangents are collinear and the surface has no normal;
/// [`SurfacePoint::normal`] reports that rather than returning a made-up
/// direction.
#[must_use]
pub fn cone_at(cone: &Cone, angle: f64, height: f64) -> SurfacePoint {
    let f = cone.frame();
    let r = cone.radius_at(height);
    let slope = cone.half_angle().tan();
    let (sin, cos) = angle.sin_cos();
    let (x, y, z) = (f.x().vector(), f.y().vector(), f.z().vector());
    SurfacePoint {
        point: f.origin() + x * (r * cos) + y * (r * sin) + z * height,
        du: x * (-r * sin) + y * (r * cos),
        dv: x * (slope * cos) + y * (slope * sin) + z,
    }
}

/// The `(angle, height)` of the point on a cone nearest `p`.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if `p` lies on the
/// axis.
pub fn cone_parameters(cone: &Cone, p: Point, tol: Tolerances) -> OgResult<(f64, f64)> {
    let local = cone.frame().to_local(p);
    if local.x.hypot(local.y) <= tol.confusion() {
        og_bail!(
            Construction,
            "point is on the cone's axis; no nearest angle"
        );
    }
    Ok((wrap_angle(local.y.atan2(local.x)), local.z))
}

/// Evaluate a sphere at `(longitude, latitude)`.
///
/// Latitude runs from `-pi/2` at the `-z` pole to `+pi/2` at `+z`.
#[must_use]
pub fn sphere_at(sphere: &Sphere, longitude: f64, latitude: f64) -> SurfacePoint {
    let f = sphere.frame();
    let r = sphere.radius();
    let (sin_lon, cos_lon) = longitude.sin_cos();
    let (sin_lat, cos_lat) = latitude.sin_cos();
    let (x, y, z) = (f.x().vector(), f.y().vector(), f.z().vector());
    let ring = r * cos_lat;
    SurfacePoint {
        point: sphere.centre() + x * (ring * cos_lon) + y * (ring * sin_lon) + z * (r * sin_lat),
        du: x * (-ring * sin_lon) + y * (ring * cos_lon),
        dv: x * (-r * sin_lat * cos_lon) + y * (-r * sin_lat * sin_lon) + z * (r * cos_lat),
    }
}

/// The `(longitude, latitude)` of the point on a sphere nearest `p`.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if `p` is the
/// centre. At a pole the latitude is well defined but the longitude is
/// arbitrary; zero is returned rather than an error, since the *position* is
/// unambiguous and callers overwhelmingly want it.
pub fn sphere_parameters(sphere: &Sphere, p: Point, tol: Tolerances) -> OgResult<(f64, f64)> {
    let local = sphere.frame().to_local(p);
    let ring = local.x.hypot(local.y);
    if ring.hypot(local.z) <= tol.confusion() {
        og_bail!(Construction, "point is the sphere's centre");
    }
    let longitude = if ring <= tol.confusion() {
        0.0
    } else {
        wrap_angle(local.y.atan2(local.x))
    };
    // atan2 of z against the ring radius, not asin of z/r: the point need not
    // be exactly on the sphere, and this stays correct and accurate when it is
    // not.
    Ok((longitude, local.z.atan2(ring)))
}

/// Evaluate a torus at `(around the axis, around the tube)`.
#[must_use]
pub fn torus_at(torus: &Torus, u: f64, v: f64) -> SurfacePoint {
    let f = torus.frame();
    let (major, minor) = (torus.major_radius(), torus.minor_radius());
    let (sin_u, cos_u) = u.sin_cos();
    let (sin_v, cos_v) = v.sin_cos();
    let (x, y, z) = (f.x().vector(), f.y().vector(), f.z().vector());
    let out = x * cos_u + y * sin_u;
    let radius = minor.mul_add(cos_v, major);
    SurfacePoint {
        point: f.origin() + out * radius + z * (minor * sin_v),
        du: (x * -sin_u + y * cos_u) * radius,
        dv: out * (-minor * sin_v) + z * (minor * cos_v),
    }
}

/// The `(u, v)` of the point on a torus nearest `p`.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if `p` lies on the
/// axis, or on the tube's centre circle where every `v` is equally near.
pub fn torus_parameters(torus: &Torus, p: Point, tol: Tolerances) -> OgResult<(f64, f64)> {
    let local = torus.frame().to_local(p);
    let ring = local.x.hypot(local.y);
    if ring <= tol.confusion() {
        og_bail!(Construction, "point is on the torus's axis");
    }
    let u = wrap_angle(local.y.atan2(local.x));
    let radial = ring - torus.major_radius();
    if radial.hypot(local.z) <= tol.confusion() {
        og_bail!(Construction, "point is on the tube's centre circle");
    }
    Ok((u, wrap_angle(local.z.atan2(radial))))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::Frame;
    use approx::assert_relative_eq;

    const T: Tolerances = Tolerances::millimetres();
    const PI: f64 = core::f64::consts::PI;

    fn tilted() -> Frame {
        Frame::new(
            Point::new(1.0, -2.0, 3.0),
            Direction::from_coords(1.0, 2.0, 3.0, T).unwrap(),
            Direction::from_coords(1.0, 0.0, 0.0, T).unwrap(),
            T,
        )
        .unwrap()
    }

    /// Check a curve's analytic derivatives against central differences.
    fn check_curve_derivatives(f: impl Fn(f64) -> CurvePoint, t: f64) {
        let h = 1e-6;
        let (a, b, c) = (f(t - h), f(t), f(t + h));
        let d1 = (c.point - a.point) * (1.0 / (2.0 * h));
        let d2 = (c.d1 - a.d1) * (1.0 / (2.0 * h));
        let scale = d1.magnitude().max(1.0);
        assert!((b.d1 - d1).magnitude() <= 1e-5 * scale, "d1 wrong at {t}");
        assert!(
            (b.d2 - d2).magnitude() <= 1e-5 * d2.magnitude().max(1.0),
            "d2 wrong at {t}"
        );
    }

    /// Check a surface's analytic derivatives against central differences.
    fn check_surface_derivatives(f: impl Fn(f64, f64) -> SurfacePoint, u: f64, v: f64) {
        let h = 1e-6;
        let du = (f(u + h, v).point - f(u - h, v).point) * (1.0 / (2.0 * h));
        let dv = (f(u, v + h).point - f(u, v - h).point) * (1.0 / (2.0 * h));
        let p = f(u, v);
        assert!(
            (p.du - du).magnitude() <= 1e-5 * du.magnitude().max(1.0),
            "du wrong"
        );
        assert!(
            (p.dv - dv).magnitude() <= 1e-5 * dv.magnitude().max(1.0),
            "dv wrong"
        );
    }

    #[test]
    fn angle_wrapping_stays_inside_its_range() {
        for a in [-10.0_f64, -PI, -1e-18, 0.0, 1.0, PI, 7.0, 100.0] {
            let w = wrap_angle(a);
            assert!(
                (0.0..core::f64::consts::TAU).contains(&w),
                "{a} wrapped to {w}"
            );
            let s = wrap_signed_angle(a);
            assert!(
                s > -PI - 1e-15 && s <= PI + 1e-15,
                "{a} signed-wrapped to {s}"
            );
        }
        // A tiny negative input is the case rem_euclid can round up to exactly
        // tau, which would put a normalized angle outside its own range.
        assert_eq!(wrap_angle(-1e-300), 0.0);
    }

    #[test]
    fn line_evaluation_and_inversion() {
        let axis = Axis::new(Point::new(1.0, 2.0, 3.0), Direction::Z);
        let c = line_at(axis, 5.0);
        assert!(c.point.is_equal(Point::new(1.0, 2.0, 8.0), T));
        assert!(c.d1.is_equal(Vector::Z, T));
        assert_relative_eq!(c.curvature(), 0.0);
        assert_relative_eq!(line_parameter(axis, c.point), 5.0, epsilon = 1e-12);
    }

    #[test]
    fn circle_evaluation_derivatives_and_inversion() {
        let c = Circle::new(tilted(), 3.0, T).unwrap();
        for i in 0..12 {
            let angle = f64::from(i) * PI / 6.0;
            let p = circle_at(&c, angle);
            assert!(c.contains(p.point, T));
            assert_relative_eq!(
                circle_parameter(&c, p.point, T).unwrap(),
                wrap_angle(angle),
                epsilon = 1e-12
            );
            check_curve_derivatives(|t| circle_at(&c, t), angle);
            // Curvature is the reciprocal of the radius, everywhere.
            assert_relative_eq!(p.curvature(), 1.0 / 3.0, epsilon = 1e-12);
        }
        assert!(circle_parameter(&c, c.centre(), T).is_err());
    }

    #[test]
    fn circle_starts_on_its_frames_x_axis() {
        // This is the whole point of carrying a frame: parameter zero is a
        // specific, reproducible place.
        let f = tilted();
        let c = Circle::new(f, 2.0, T).unwrap();
        assert!(
            circle_at(&c, 0.0)
                .point
                .is_equal(c.centre() + f.x() * 2.0, T)
        );
        assert!(
            circle_at(&c, PI / 2.0)
                .point
                .is_equal(c.centre() + f.y() * 2.0, T)
        );
    }

    #[test]
    fn ellipse_evaluation_derivatives_and_inversion() {
        let e = Ellipse::new(tilted(), 5.0, 3.0, T).unwrap();
        for i in 0..12 {
            let angle = f64::from(i) * PI / 6.0;
            let p = ellipse_at(&e, angle);
            assert_relative_eq!(
                ellipse_parameter(&e, p.point, T).unwrap(),
                wrap_angle(angle),
                epsilon = 1e-12
            );
            check_curve_derivatives(|t| ellipse_at(&e, t), angle);
        }
        // The eccentric angle is not the polar one; at 45 degrees eccentric the
        // point is not at 45 degrees polar.
        let p = ellipse_at(&e, PI / 4.0);
        let local = e.frame().to_local(p.point);
        assert_relative_eq!(
            local.x,
            5.0 * core::f64::consts::FRAC_1_SQRT_2,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            local.y,
            3.0 * core::f64::consts::FRAC_1_SQRT_2,
            epsilon = 1e-12
        );
    }

    #[test]
    fn ellipse_curvature_is_extreme_at_the_ends_of_its_axes() {
        let e = Ellipse::new(Frame::WORLD, 5.0, 3.0, T).unwrap();
        // At the end of the major axis, curvature is b/a^2 * a... = a/b^2 form:
        // kappa = a / b^2 at the minor-axis end, b / a^2 at the major-axis end.
        assert_relative_eq!(ellipse_at(&e, 0.0).curvature(), 5.0 / 9.0, epsilon = 1e-12);
        assert_relative_eq!(
            ellipse_at(&e, PI / 2.0).curvature(),
            3.0 / 25.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn hyperbola_evaluation_and_inversion() {
        let h = Hyperbola::new(tilted(), 3.0, 4.0, T).unwrap();
        for t in [-2.0_f64, -0.5, 0.0, 0.5, 2.0] {
            let p = hyperbola_at(&h, t);
            assert_relative_eq!(
                hyperbola_parameter(&h, p.point, T).unwrap(),
                t,
                epsilon = 1e-11
            );
            check_curve_derivatives(|s| hyperbola_at(&h, s), t);
        }
        assert!(hyperbola_at(&h, 0.0).point.is_equal(h.vertex(), T));
        // The far branch is a different curve.
        let far = h.centre() - h.frame().x() * 5.0;
        assert!(hyperbola_parameter(&h, far, T).is_err());
    }

    #[test]
    fn hyperbola_inversion_is_accurate_near_the_vertex() {
        // The reason for asinh rather than acosh: at t near zero, acosh's
        // argument approaches 1 where its derivative is unbounded, so it loses
        // most of its precision exactly where curves are usually trimmed.
        let h = Hyperbola::new(Frame::WORLD, 3.0, 4.0, T).unwrap();
        for t in [1e-8_f64, 1e-5, 1e-3] {
            let p = hyperbola_at(&h, t);
            assert_relative_eq!(
                hyperbola_parameter(&h, p.point, T).unwrap(),
                t,
                max_relative = 1e-9
            );
        }
    }

    #[test]
    fn parabola_evaluation_and_inversion() {
        let p = Parabola::new(tilted(), 2.0, T).unwrap();
        for t in [-4.0_f64, -1.0, 0.0, 1.0, 4.0] {
            let c = parabola_at(&p, t);
            assert_relative_eq!(parabola_parameter(&p, c.point), t, epsilon = 1e-11);
            check_curve_derivatives(|s| parabola_at(&p, s), t);
        }
        assert!(parabola_at(&p, 0.0).point.is_equal(p.apex(), T));
        // Every point is equidistant from the focus and the directrix.
        let focus = p.focus();
        for t in [-3.0_f64, 1.0, 5.0] {
            let point = parabola_at(&p, t).point;
            let local = p.frame().to_local(point);
            let to_directrix = local.x + p.focal();
            assert_relative_eq!(point.distance(focus), to_directrix, epsilon = 1e-11);
        }
    }

    #[test]
    fn plane_evaluation_and_inversion() {
        let plane = Plane::new(tilted());
        for (u, v) in [(0.0, 0.0), (3.0, -2.0), (-100.0, 50.0)] {
            let p = plane_at(&plane, u, v);
            assert!(plane.contains(p.point, T));
            let (bu, bv) = plane_parameters(&plane, p.point);
            assert_relative_eq!(bu, u, epsilon = 1e-11);
            assert_relative_eq!(bv, v, epsilon = 1e-11);
            assert!(p.normal(T).unwrap().is_equal(plane.normal(), T));
        }
        check_surface_derivatives(|u, v| plane_at(&plane, u, v), 1.0, 2.0);
    }

    #[test]
    fn cylinder_evaluation_inversion_and_normal() {
        let c = Cylinder::new(tilted(), 2.0, T).unwrap();
        for i in 0..8 {
            let angle = f64::from(i) * PI / 4.0;
            for h in [-5.0_f64, 0.0, 7.0] {
                let p = cylinder_at(&c, angle, h);
                assert!(c.contains(p.point, T));
                let (ba, bh) = cylinder_parameters(&c, p.point, T).unwrap();
                assert_relative_eq!(ba, wrap_angle(angle), epsilon = 1e-11);
                assert_relative_eq!(bh, h, epsilon = 1e-11);
                // The normal is radial, so perpendicular to the axis.
                assert!(p.normal(T).unwrap().dot(c.frame().z()).abs() < 1e-12);
            }
        }
        check_surface_derivatives(|u, v| cylinder_at(&c, u, v), 0.7, 3.0);
        assert!(cylinder_parameters(&c, c.frame().origin(), T).is_err());
    }

    #[test]
    fn cone_evaluation_inversion_and_apex_degeneracy() {
        let c = Cone::new(tilted(), 3.0, 0.6, T).unwrap();
        for i in 0..8 {
            let angle = f64::from(i) * PI / 4.0;
            for h in [-1.0_f64, 0.0, 4.0] {
                let p = cone_at(&c, angle, h);
                assert!(
                    c.contains(p.point, T),
                    "distance {}",
                    c.distance_to(p.point)
                );
                let (ba, bh) = cone_parameters(&c, p.point, T).unwrap();
                assert_relative_eq!(ba, wrap_angle(angle), epsilon = 1e-11);
                assert_relative_eq!(bh, h, epsilon = 1e-11);
            }
        }
        check_surface_derivatives(|u, v| cone_at(&c, u, v), 0.7, 2.0);

        // At the apex the radius is zero, so the u-tangent vanishes and there
        // is no normal. Reporting that beats inventing a direction.
        let apex_height = -3.0 / 0.6_f64.tan();
        let at_apex = cone_at(&c, 1.0, apex_height);
        assert!(at_apex.point.is_equal(c.apex(), T));
        assert!(at_apex.is_degenerate(T));
        assert!(at_apex.normal(T).is_err());
    }

    #[test]
    fn sphere_evaluation_inversion_and_poles() {
        let s = Sphere::new(tilted(), 4.0, T).unwrap();
        for i in 0..8 {
            let lon = f64::from(i) * PI / 4.0;
            for lat in [-1.2_f64, -0.4, 0.0, 0.9] {
                let p = sphere_at(&s, lon, lat);
                assert!(s.contains(p.point, T));
                let (blon, blat) = sphere_parameters(&s, p.point, T).unwrap();
                assert_relative_eq!(blon, wrap_angle(lon), epsilon = 1e-11);
                assert_relative_eq!(blat, lat, epsilon = 1e-11);
                // The normal is radial.
                assert!(
                    p.normal(T)
                        .unwrap()
                        .is_equal(s.normal_at(p.point, T).unwrap(), T)
                );
            }
        }
        check_surface_derivatives(|u, v| sphere_at(&s, u, v), 1.1, 0.3);

        // At a pole the position is unambiguous even though the longitude is
        // not, so inversion succeeds and picks zero.
        let north = sphere_at(&s, 2.0, PI / 2.0);
        assert!(north.point.is_equal(s.centre() + s.frame().z() * 4.0, T));
        assert!(north.is_degenerate(T));
        let (lon, lat) = sphere_parameters(&s, north.point, T).unwrap();
        assert_relative_eq!(lon, 0.0);
        assert_relative_eq!(lat, PI / 2.0, epsilon = 1e-8);
        assert!(sphere_parameters(&s, s.centre(), T).is_err());
    }

    #[test]
    fn torus_evaluation_inversion_and_normal() {
        let t = Torus::new(tilted(), 5.0, 2.0, T).unwrap();
        for i in 0..6 {
            let u = f64::from(i) * PI / 3.0;
            for j in 0..6 {
                let v = f64::from(j) * PI / 3.0;
                let p = torus_at(&t, u, v);
                assert!(
                    t.contains(p.point, T),
                    "distance {}",
                    t.distance_to(p.point)
                );
                let (bu, bv) = torus_parameters(&t, p.point, T).unwrap();
                assert_relative_eq!(bu, wrap_angle(u), epsilon = 1e-10);
                assert_relative_eq!(bv, wrap_angle(v), epsilon = 1e-10);
                assert!(p.normal(T).is_ok());
            }
        }
        check_surface_derivatives(|u, v| torus_at(&t, u, v), 0.7, 2.0);
        assert!(torus_parameters(&t, t.centre(), T).is_err());
    }

    #[test]
    fn curvature_of_a_circle_is_the_reciprocal_of_its_radius() {
        for r in [0.1_f64, 1.0, 100.0] {
            let c = Circle::new(Frame::WORLD, r, T).unwrap();
            assert_relative_eq!(
                circle_at(&c, 1.3).curvature(),
                1.0 / r,
                max_relative = 1e-12
            );
        }
    }

    #[test]
    fn tangent_of_a_circle_is_perpendicular_to_its_radius() {
        let c = Circle::new(tilted(), 3.0, T).unwrap();
        for i in 0..8 {
            let angle = f64::from(i) * PI / 4.0;
            let p = circle_at(&c, angle);
            let radius = p.point - c.centre();
            assert!(p.tangent(T).unwrap().dot_vector(radius).abs() < 1e-12);
        }
    }
}
