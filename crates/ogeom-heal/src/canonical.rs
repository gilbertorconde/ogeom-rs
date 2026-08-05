//! Canonical recognition: deciding that a free-form patch *is* a plane, a
//! cylinder, a cone, a sphere or a torus — not resembles one.
//!
//! The input is samples with normals; the output is the canonical surface
//! and the worst deviation actually measured, or nothing. The deciding rule
//! is the deferred entry's: a fit is easy, the *decision* is the product,
//! and a wrong yes gives a solid that measures nearly right with the wrong
//! surface under every later operation. So every candidate is verified
//! against all the samples at the caller's stated tolerance, and the
//! reported deviation is the certificate.
//!
//! The estimators are classical. A plane is the point covariance's smallest
//! direction. A sphere is linear least squares through the `|c|² − r²`
//! substitution. A cylinder's axis is the direction the *normals* avoid —
//! their covariance's smallest eigenvector — and its section is a 2D circle
//! fit. A cone adds the linear taper of radius against height. A torus is
//! the one genuinely nonlinear case and takes a few Gauss–Newton steps from
//! the cylinder-style seed.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Surface as _, SurfaceGeometry};
use ogeom_math::{Cone, Cylinder, Direction, Frame, Plane, Point, Sphere, Torus, Vector};

/// A canonical surface a patch was recognized as.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Canonical {
    /// A plane.
    Plane(Plane),
    /// A cylinder.
    Cylinder(Cylinder),
    /// A cone.
    Cone(Cone),
    /// A sphere.
    Sphere(Sphere),
    /// A torus.
    Torus(Torus),
}

/// A recognition with its certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Recognized {
    /// What the samples are.
    pub surface: Canonical,
    /// The worst distance from any sample to it — measured, not promised.
    pub deviation: f64,
}

/// Recognize a canonical surface from samples with unit normals.
///
/// Candidates are tried simplest first — plane, sphere, cylinder, cone,
/// torus — and the first whose *measured* worst deviation meets `tolerance`
/// wins; `None` says the samples are genuinely free-form at that tolerance,
/// which is an answer, not a failure.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// fewer than six samples arrive or the tolerance is not positive.
pub fn recognize_points(
    points: &[Point],
    normals: &[Vector],
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Recognized>> {
    if points.len() < 6 || points.len() != normals.len() {
        ogeom_bail!(
            Construction,
            "recognition needs at least six samples with matching normals"
        );
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        ogeom_bail!(Construction, "a tolerance of {tolerance} is not a distance");
    }

    let candidates = [
        fit_plane(points, tol),
        fit_sphere(points, tol),
        fit_cylinder(points, normals, tol),
        fit_cone(points, normals, tol),
        fit_torus(points, normals, tol),
    ];
    for candidate in candidates.into_iter().flatten() {
        let deviation = worst_deviation(&candidate, points);
        if deviation <= tolerance {
            return Ok(Some(Recognized {
                surface: candidate,
                deviation,
            }));
        }
    }
    Ok(None)
}

/// Recognize a surface by sampling it over its own domain.
///
/// # Errors
///
/// As [`recognize_points`], plus whatever evaluation refuses.
pub fn recognize_surface(
    surface: &SurfaceGeometry,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Recognized>> {
    let ((ua, ub), (va, vb)) = surface.domain();
    let mut points = Vec::new();
    let mut normals = Vec::new();
    let n = 9;
    for i in 0..=n {
        for j in 0..=n {
            let u = ua + (ub - ua) * f64::from(i) / f64::from(n);
            let v = va + (vb - va) * f64::from(j) / f64::from(n);
            let Ok(p) = surface.point_at(u, v, tol) else {
                continue;
            };
            let Ok(normal) = surface.normal_at(u, v, tol) else {
                continue;
            };
            points.push(p);
            normals.push(normal.vector());
        }
    }
    recognize_points(&points, &normals, tolerance, tol)
}

/// The worst distance from any sample to the candidate.
fn worst_deviation(candidate: &Canonical, points: &[Point]) -> f64 {
    points
        .iter()
        .map(|p| match candidate {
            Canonical::Plane(plane) => plane.signed_distance_to(*p).abs(),
            Canonical::Cylinder(c) => c.distance_to(*p),
            Canonical::Cone(c) => c.distance_to(*p),
            Canonical::Sphere(s) => (p.distance(s.centre()) - s.radius()).abs(),
            Canonical::Torus(t) => t.distance_to(*p),
        })
        .fold(0.0, f64::max)
}

fn centroid(points: &[Point]) -> Vector {
    let mut sum = Vector::ZERO;
    for p in points {
        sum += p.to_vector();
    }
    #[allow(clippy::cast_precision_loss, reason = "sample counts are small")]
    let n = points.len() as f64;
    sum / n
}

/// The eigenvector of a 3×3 symmetric matrix for its smallest eigenvalue.
fn smallest_direction(m: nalgebra::Matrix3<f64>) -> Option<Vector> {
    let eigen = nalgebra::SymmetricEigen::new(m);
    let mut best = 0;
    for i in 1..3 {
        if eigen.eigenvalues[i] < eigen.eigenvalues[best] {
            best = i;
        }
    }
    let v = eigen.eigenvectors.column(best);
    Some(Vector::new(v[0], v[1], v[2]))
}

fn covariance(vectors: impl Iterator<Item = Vector>) -> nalgebra::Matrix3<f64> {
    let mut m = nalgebra::Matrix3::zeros();
    for v in vectors {
        let n = nalgebra::Vector3::new(v.x, v.y, v.z);
        m += n * n.transpose();
    }
    m
}

fn fit_plane(points: &[Point], tol: Tolerances) -> Option<Canonical> {
    let c = centroid(points);
    let m = covariance(points.iter().map(|p| p.to_vector() - c));
    let normal = Direction::new(smallest_direction(m)?, tol).ok()?;
    let frame = Frame::about(Point::from_vector(c), normal);
    Some(Canonical::Plane(Plane::new(frame)))
}

fn fit_sphere(points: &[Point], tol: Tolerances) -> Option<Canonical> {
    // |p|² − 2p·c + Q = 0 with Q = |c|² − r²: linear in (c, Q).
    let mut a = nalgebra::Matrix4::zeros();
    let mut b = nalgebra::Vector4::zeros();
    for p in points {
        let row = nalgebra::Vector4::new(-2.0 * p.x, -2.0 * p.y, -2.0 * p.z, 1.0);
        let rhs = -(p.to_vector().dot(p.to_vector()));
        a += row * row.transpose();
        b += row * rhs;
    }
    let solved = a.lu().solve(&b)?;
    let centre = Point::new(solved[0], solved[1], solved[2]);
    let r2 = centre.to_vector().dot(centre.to_vector()) - solved[3];
    if r2 <= tol.confusion() {
        return None;
    }
    Some(Canonical::Sphere(
        Sphere::centred(centre, r2.sqrt(), tol).ok()?,
    ))
}

/// Axis, then the section circle in the axis plane.
fn fit_cylinder(points: &[Point], normals: &[Vector], tol: Tolerances) -> Option<Canonical> {
    // A cylinder's normals all lie perpendicular to the axis.
    let mean = {
        let mut sum = Vector::ZERO;
        for n in normals {
            sum += *n;
        }
        #[allow(clippy::cast_precision_loss, reason = "sample counts are small")]
        let count = normals.len() as f64;
        sum / count
    };
    let axis = Direction::new(
        smallest_direction(covariance(normals.iter().map(|n| *n - mean)))?,
        tol,
    )
    .ok()?;
    let (centre2, radius) = section_circle(points, axis)?;
    let frame = Frame::about(centre2, axis);
    Some(Canonical::Cylinder(Cylinder::new(frame, radius, tol).ok()?))
}

/// Project points into the plane perpendicular to `axis` through the
/// centroid and fit a circle there — shared by the cylinder and the seed of
/// the torus.
fn section_circle(points: &[Point], axis: Direction) -> Option<(Point, f64)> {
    let c = centroid(points);
    let a = axis.vector();
    // An in-plane orthonormal basis: the 2D circle problem solved as one.
    let seed = if a.x.abs() < 0.9 {
        Vector::new(1.0, 0.0, 0.0)
    } else {
        Vector::new(0.0, 1.0, 0.0)
    };
    let e1 = {
        let v = seed - a * seed.dot(a);
        v / v.magnitude()
    };
    let e2 = a.cross(e1);
    let mut m = nalgebra::Matrix3::zeros();
    let mut b = nalgebra::Vector3::zeros();
    for p in points {
        let d = p.to_vector() - c;
        let flat = d - a * d.dot(a);
        let (x, y) = (flat.dot(e1), flat.dot(e2));
        let row = nalgebra::Vector3::new(-2.0 * x, -2.0 * y, 1.0);
        let rhs = -x.mul_add(x, y * y);
        m += row * row.transpose();
        b += row * rhs;
    }
    let solved = m.lu().solve(&b)?;
    let r2 = solved[0].mul_add(solved[0], solved[1] * solved[1]) - solved[2];
    if r2 <= 0.0 {
        return None;
    }
    Some((
        Point::from_vector(c + e1 * solved[0] + e2 * solved[1]),
        r2.sqrt(),
    ))
}

fn fit_cone(points: &[Point], normals: &[Vector], tol: Tolerances) -> Option<Canonical> {
    // The axis: same estimator as the cylinder — a cone's normals keep a
    // constant angle with the axis, so their variance along it is least.
    let mean = {
        let mut sum = Vector::ZERO;
        for n in normals {
            sum += *n;
        }
        #[allow(clippy::cast_precision_loss, reason = "sample counts are small")]
        let count = normals.len() as f64;
        sum / count
    };
    let axis = Direction::new(
        smallest_direction(covariance(normals.iter().map(|n| *n - mean)))?,
        tol,
    )
    .ok()?;
    let a = axis.vector();
    let c = centroid(points);
    // Radius against height is a line: ρ = k·h + ρ₀.
    let (mut sh, mut shh, mut sr, mut shr, mut count) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    for p in points {
        let d = p.to_vector() - c;
        let h = d.dot(a);
        let rho = (d - a * h).magnitude();
        sh += h;
        shh += h * h;
        sr += rho;
        shr += h * rho;
        count += 1.0;
    }
    let det = count.mul_add(shh, -(sh * sh));
    if det.abs() <= f64::MIN_POSITIVE {
        return None;
    }
    let k = count.mul_add(shr, -(sh * sr)) / det;
    let rho0 = shh.mul_add(sr, -(sh * shr)) / det;
    if k.abs() <= tol.angular() {
        // No taper is a cylinder, which was already tried.
        return None;
    }
    let half_angle = k.atan();
    // Also centre the axis transversally through the section fit.
    let (centre, _) = section_circle(points, axis)?;
    let through =
        Point::from_vector(centre.to_vector() + a * (c.dot(a) - centre.to_vector().dot(a)));
    let frame = Frame::about(through, axis);
    Some(Canonical::Cone(
        Cone::new(frame, half_angle.abs(), rho0.max(tol.confusion()), tol).ok()?,
    ))
}

fn fit_torus(points: &[Point], normals: &[Vector], tol: Tolerances) -> Option<Canonical> {
    // The axis, from the normals' own geometry: every torus normal line
    // passes through the tube-centre circle, so the points `p − r·n` become
    // *coplanar* exactly at the true tube radius. A one-dimensional search
    // over r — both signs, since the sampled normals may face either way —
    // finds the radius as the flattest configuration, and the flat plane's
    // normal is the axis. No sampling distribution can bias it.
    let spread = {
        let c = centroid(points);
        points
            .iter()
            .map(|p| (p.to_vector() - c).magnitude())
            .fold(0.0, f64::max)
    };
    let flatness = |r: f64| -> Option<(f64, Vector, Vector)> {
        let shifted: Vec<Point> = points
            .iter()
            .zip(normals)
            .map(|(p, n)| *p - *n * r)
            .collect();
        let c = centroid(&shifted);
        let m = covariance(shifted.iter().map(|p| p.to_vector() - c));
        let eigen = nalgebra::SymmetricEigen::new(m);
        let mut best = 0;
        for i in 1..3 {
            if eigen.eigenvalues[i] < eigen.eigenvalues[best] {
                best = i;
            }
        }
        let v = eigen.eigenvectors.column(best);
        Some((
            eigen.eigenvalues[best].max(0.0),
            Vector::new(v[0], v[1], v[2]),
            c,
        ))
    };
    let mut best: Option<(f64, f64, Vector, Vector)> = None;
    let steps = 400;
    for k in 0..=steps {
        let r = -spread + 2.0 * spread * f64::from(k) / f64::from(steps);
        if r.abs() <= tol.confusion() {
            continue;
        }
        let (flat, axis_v, centre_v) = flatness(r)?;
        if best.is_none_or(|(held, ..)| flat < held) {
            best = Some((flat, r, axis_v, centre_v));
        }
    }
    let (_, tube_signed, axis_v, _) = best?;
    // Refine the radius by golden-ratio narrowing around the winner.
    let mut lo = tube_signed - 2.0 * spread / f64::from(steps);
    let mut hi = tube_signed + 2.0 * spread / f64::from(steps);
    for _ in 0..48 {
        let m1 = lo + (hi - lo) * 0.382;
        let m2 = lo + (hi - lo) * 0.618;
        let (f1, ..) = flatness(m1)?;
        let (f2, ..) = flatness(m2)?;
        if f1 < f2 {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let tube_signed = f64::midpoint(lo, hi);
    let (_, axis_refined, _) = flatness(tube_signed)?;
    let _ = axis_v;
    let axis = Direction::new(axis_refined, tol).ok()?;

    // The shifted points sit *exactly* on the mid-circle, so fitting that
    // circle reads the centre and major radius without any distribution
    // bias — an algebraic circle fit is exact for on-circle points however
    // they cluster.
    let mid: Vec<Point> = points
        .iter()
        .zip(normals)
        .map(|(p, n)| *p - *n * tube_signed)
        .collect();
    let (centre, major) = section_circle(&mid, axis)?;
    let minor = tube_signed.abs();
    if minor <= tol.confusion() || minor >= major {
        return None;
    }
    Some(Canonical::Torus(
        Torus::new(Frame::about(centre, axis), major, minor, tol).ok()?,
    ))
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_geom::{CylinderSurface, SphereSurface, TorusSurface};

    const T: Tolerances = Tolerances::millimetres();

    fn recognized(surface: &SurfaceGeometry) -> Recognized {
        recognize_surface(surface, 1e-6, T)
            .unwrap()
            .expect("the canonical should be recognized")
    }

    #[test]
    fn splined_canonicals_confess_what_they_are() {
        // Each canonical, converted to its NURBS restatement, must be
        // recognized back exactly — the round trip the healing entry wants.
        let cylinder = SurfaceGeometry::Cylinder(
            CylinderSurface::new(Cylinder::new(Frame::WORLD, 4.0, T).unwrap(), (0.0, 10.0))
                .unwrap(),
        );
        let cyl_spline: SurfaceGeometry = cylinder.to_bspline(T).unwrap().into();
        let found = recognized(&cyl_spline);
        let Canonical::Cylinder(c) = found.surface else {
            panic!("a splined cylinder is a cylinder: {found:?}");
        };
        assert!((c.radius() - 4.0).abs() < 1e-6, "{c:?}");

        let sphere = SurfaceGeometry::Sphere(SphereSurface::new(
            Sphere::centred(Point::new(1.0, 2.0, 3.0), 5.0, T).unwrap(),
        ));
        let sphere_spline: SurfaceGeometry = sphere.to_bspline(T).unwrap().into();
        let found = recognized(&sphere_spline);
        let Canonical::Sphere(s) = found.surface else {
            panic!("a splined sphere is a sphere: {found:?}");
        };
        assert!((s.radius() - 5.0).abs() < 1e-6);
        assert!(s.centre().distance(Point::new(1.0, 2.0, 3.0)) < 1e-6);

        let torus = SurfaceGeometry::Torus(TorusSurface::new(
            Torus::new(Frame::WORLD, 10.0, 2.0, T).unwrap(),
        ));
        let torus_spline: SurfaceGeometry = torus.to_bspline(T).unwrap().into();
        let found = recognize_surface(&torus_spline, 1e-4, T).unwrap().unwrap();
        let Canonical::Torus(t) = found.surface else {
            panic!("a splined torus is a torus: {found:?}");
        };
        assert!((t.major_radius() - 10.0).abs() < 1e-3, "{t:?}");
        assert!((t.minor_radius() - 2.0).abs() < 1e-3, "{t:?}");
    }

    #[test]
    fn a_free_form_patch_refuses_every_canonical() {
        // The saddle z = x·y over a square: nothing canonical fits it.
        let mut points = Vec::new();
        let mut normals = Vec::new();
        for i in 0..=9 {
            for j in 0..=9 {
                let (x, y) = (f64::from(i) - 4.5, f64::from(j) - 4.5);
                points.push(Point::new(x, y, x * y));
                let n = Vector::new(-y, -x, 1.0);
                normals.push(n / n.magnitude());
            }
        }
        assert!(
            recognize_points(&points, &normals, 1e-3, T)
                .unwrap()
                .is_none(),
            "a saddle is not canonical"
        );
    }
}
