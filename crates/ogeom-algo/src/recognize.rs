//! Canonical recognition: deciding that a set of points *is* a plane, a
//! cylinder, a cone, a sphere or a torus, not that it resembles one.
//!
//! The input is samples with normals; the output is the canonical surface
//! and the worst deviation actually measured, or nothing. A fit is easy;
//! the decision is the product, and a wrong yes gives a solid that measures
//! nearly right with the wrong surface under every later operation. So
//! every candidate is verified against all the samples at the caller's
//! stated tolerance, and the reported deviation is the certificate.
//!
//! The first estimates are closed forms. A plane is the point covariance's
//! smallest direction, and a sphere linear least squares through the
//! `|c|² − r²` substitution. A cylinder, a cone and a torus are surfaces of
//! revolution, every normal line of which meets the axis: the axis is the
//! line that best meets them all, found linearly in Plücker coordinates and
//! reweighted against the few samples off the surface, and the kind's
//! profile (a line, a slanted line) is then fitted in the plane through
//! it. A torus's tube radius is read from how fast its normals turn, and its
//! spine as the circle the samples land on when shifted back along their
//! normals by that radius. Normals estimated from a mesh's facets are only
//! good to a fraction of the facet angle, so each estimate is refined by
//! least squares on the samples' own distances to the surface, which is
//! what the verification measures.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Cone, Cylinder, Direction, Frame, Plane, Point, Sphere, Torus, Vector};

/// A canonical surface a set of samples was recognized as.
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

impl Canonical {
    /// The distance from `p` to the surface.
    #[must_use]
    pub fn distance_to(&self, p: Point) -> f64 {
        match self {
            Self::Plane(plane) => plane.signed_distance_to(p).abs(),
            Self::Cylinder(c) => c.distance_to(p),
            Self::Cone(c) => c.distance_to(p),
            Self::Sphere(s) => (p.distance(s.centre()) - s.radius()).abs(),
            Self::Torus(t) => t.distance_to(p),
        }
    }

    /// The distance from `p` to the surface, positive on the side its
    /// radius grows toward (outside a cylinder, sphere or torus, off the
    /// axis of a cone) and along a plane's normal. Near the surface it
    /// grows as the distance does, which is what a solve onto the surface
    /// asks of it; a cone's is taken to its nearer nappe.
    #[must_use]
    pub fn signed_distance_to(&self, p: Point) -> f64 {
        let radial = |frame: ogeom_math::Frame| {
            let w = p - frame.origin();
            let h = w.dot(frame.z().vector());
            ((w - frame.z().vector() * h).magnitude(), h)
        };
        match self {
            Self::Plane(plane) => plane.signed_distance_to(p),
            Self::Cylinder(c) => radial(c.frame()).0 - c.radius(),
            Self::Cone(c) => {
                let (r, h) = radial(c.frame());
                let (sin, cos) = c.half_angle().sin_cos();
                // In the half-plane through the axis the nappe is the line
                // through (reference radius, 0) at the half angle.
                (r - c.reference_radius()).mul_add(cos, -(h * sin))
            }
            Self::Sphere(s) => p.distance(s.centre()) - s.radius(),
            Self::Torus(t) => {
                let (r, h) = radial(t.frame());
                (r - t.major_radius()).hypot(h) - t.minor_radius()
            }
        }
    }
}

/// A recognition with its certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Recognized {
    /// What the samples are.
    pub surface: Canonical,
    /// The worst distance from any sample to it: measured, not promised.
    pub deviation: f64,
}

/// Recognize a canonical surface from samples with unit normals.
///
/// A plane is tried first; otherwise every curved kind is fitted, and the
/// closest whose *measured* worst deviation meets `tolerance` is taken, a
/// simpler kind winning a tie. `None` says the samples are free-form at
/// that tolerance, which is an answer, not a failure. Each kind has its own
/// sample floor, below which the fit is underdetermined and verification
/// would rubber-stamp whatever the algebra produced.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// fewer than three samples arrive, the normals do not match the points,
/// or the tolerance is not a positive distance.
pub fn recognize_points(
    points: &[Point],
    normals: &[Vector],
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Recognized>> {
    if points.len() < 3 || points.len() != normals.len() {
        ogeom_bail!(
            Construction,
            "recognition needs at least three samples with matching normals"
        );
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        ogeom_bail!(Construction, "a tolerance of {tolerance} is not a distance");
    }
    if let Some(plane) = fit_plane(points, tol) {
        let deviation = worst_deviation(&plane, points);
        if deviation <= tolerance {
            return Ok(Some(Recognized {
                surface: plane,
                deviation,
            }));
        }
    }
    Ok(recognize_curved(points, normals, &[], tolerance, tol))
}

/// Whether the samples lie within `tolerance` of one plane.
pub(crate) fn is_flat(points: &[Point], tolerance: f64, tol: Tolerances) -> bool {
    fit_plane(points, tol).is_some_and(|plane| worst_deviation(&plane, points) <= tolerance)
}

/// As [`recognize_points`], among the curved kinds only.
///
/// `chords` are straight segments between samples (a mesh's edges),
/// which break ties among fits: samples on two equal coaxial circles lie
/// as exactly on a sphere as on a cylinder, but the cylinder's rulings run
/// between them and the mesh draws those as edges on the surface, midpoint
/// and all, where no chord of a sphere lies on it.
pub(crate) fn recognize_curved(
    points: &[Point],
    normals: &[Vector],
    chords: &[(Point, Point)],
    tolerance: f64,
    tol: Tolerances,
) -> Option<Recognized> {
    choose(
        &curved_fits(points, normals, tolerance, tol),
        chords,
        tolerance,
    )
}

/// Every curved kind fitted and refined, with its worst deviation over all
/// the samples. Fitted on an even subsample of a large set and verified on
/// all of it: the fit's cost grows with every point, the certificate's
/// only by a distance each.
fn curved_fits(
    points: &[Point],
    normals: &[Vector],
    hopeless: f64,
    tol: Tolerances,
) -> Vec<Recognized> {
    const FIT_SAMPLES: usize = 1500;
    let stride = points.len().div_ceil(FIT_SAMPLES).max(1);
    let (sub_points, sub_normals): (Vec<Point>, Vec<Vector>) = points
        .iter()
        .zip(normals)
        .step_by(stride)
        .map(|(p, n)| (*p, *n))
        .unzip();
    let n = points.len();
    type Fit = fn(&[Point], &[Vector], Tolerances) -> Option<Canonical>;
    // Each kind is fitted only to more samples than it has parameters, by
    // half as many again: a torus's seven fit eight points on a sphere and
    // a cylinder's end together exactly, and says nothing by it.
    let attempts: [(usize, Fit); 4] = [
        (6, |p, _, tol| fit_sphere(p, tol)),
        (8, fit_cylinder),
        (9, fit_cone),
        (11, fit_torus),
    ];
    let mut fits = Vec::with_capacity(4);
    for (floor, fit) in attempts {
        if n < floor {
            continue;
        }
        let Some(seed) = fit(&sub_points, &sub_normals, tol) else {
            continue;
        };
        let refined = refine(seed, &sub_points, hopeless, tol).unwrap_or(seed);
        fits.push(Recognized {
            surface: refined,
            deviation: worst_deviation(&refined, points),
        });
    }
    fits
}

/// The fit to take among those within the tolerance: the closest (a
/// small patch of a torus may sit within the tolerance of a sphere too,
/// and the sphere would not extrapolate) unless a simpler kind ties it,
/// within twice the best or at the noise floor of a thousandth of the
/// tolerance, since a torus can mimic a cylinder as closely as it likes.
///
/// Before either, the fit that lays the most of the chords' midpoints on
/// itself: a mesh's edges along a surface's rulings are on the surface.
fn choose(fits: &[Recognized], chords: &[(Point, Point)], tolerance: f64) -> Option<Recognized> {
    let within: Vec<&Recognized> = fits.iter().filter(|f| f.deviation <= tolerance).collect();
    let on = |f: &Recognized| {
        chords
            .iter()
            .filter(|(a, b)| {
                let mid = Point::from_vector((a.to_vector() + b.to_vector()) * 0.5);
                f.surface.distance_to(mid) <= tolerance
            })
            .count()
    };
    let most = within.iter().map(|f| on(f)).max()?;
    let within: Vec<&Recognized> = within.into_iter().filter(|f| on(f) == most).collect();
    let best = within
        .iter()
        .map(|f| f.deviation)
        .fold(f64::INFINITY, f64::min);
    within
        .into_iter()
        .find(|f| f.deviation <= (2.0 * best).max(tolerance * 1e-3))
        .copied()
}

/// A curved surface through samples of which a few are not on it: fitted,
/// the tenth that miss the closest fit farthest dropped (a few corners
/// far off the surface pull a fit a long way, so a median cut would not
/// isolate them), and fitted again, until what is left fits within the
/// tolerance or too little is left. Returns the surface and which samples
/// it keeps, or, when nothing fits, how close the closest fit came.
pub(crate) fn recognize_trimmed(
    points: &[Point],
    normals: &[Vector],
    chords: &[(Point, Point)],
    tolerance: f64,
    tol: Tolerances,
) -> Result<(Recognized, Vec<bool>), f64> {
    let mut keep = vec![true; points.len()];
    let mut closest = f64::INFINITY;
    for _ in 0..3 {
        let (kept_points, kept_normals): (Vec<Point>, Vec<Vector>) = points
            .iter()
            .zip(normals)
            .zip(&keep)
            .filter(|(_, k)| **k)
            .map(|((p, n), _)| (*p, *n))
            .unzip();
        if kept_points.len() < 8 {
            break;
        }
        let fits = curved_fits(&kept_points, &kept_normals, tolerance, tol);
        if let Some(found) = choose(&fits, chords, tolerance) {
            return Ok((found, keep));
        }
        let Some(candidate) = fits
            .iter()
            .min_by(|a, b| a.deviation.total_cmp(&b.deviation))
        else {
            break;
        };
        closest = closest.min(candidate.deviation);
        let mut misses: Vec<f64> = kept_points
            .iter()
            .map(|p| candidate.surface.distance_to(*p))
            .collect();
        misses.sort_by(f64::total_cmp);
        // Trimming drops a few samples off a surface the rest are on; when
        // the typical sample misses too, the rest are on no surface.
        if misses[misses.len() / 2] > tolerance * 10.0 {
            break;
        }
        let bound = misses[misses.len() * 9 / 10].max(tolerance);
        let before = keep.iter().filter(|k| **k).count();
        for (k, p) in keep.iter_mut().zip(points) {
            *k = *k && candidate.surface.distance_to(*p) <= bound;
        }
        if keep.iter().filter(|k| **k).count() == before {
            break;
        }
    }
    Err(closest)
}

/// The worst distance from any sample to the candidate.
pub(crate) fn worst_deviation(candidate: &Canonical, points: &[Point]) -> f64 {
    points
        .iter()
        .map(|p| candidate.distance_to(*p))
        .fold(0.0, f64::max)
}

fn centroid(points: &[Point]) -> Vector {
    let mut sum = Vector::ZERO;
    for p in points {
        sum += p.to_vector();
    }
    #[allow(
        clippy::cast_precision_loss,
        reason = "sample counts are far below 2^52"
    )]
    let n = points.len() as f64;
    sum / n
}

/// The eigenvector of a 3×3 symmetric matrix for its smallest eigenvalue.
fn smallest_direction(m: nalgebra::Matrix3<f64>) -> Vector {
    let eigen = nalgebra::SymmetricEigen::new(m);
    let mut best = 0;
    for i in 1..3 {
        if eigen.eigenvalues[i] < eigen.eigenvalues[best] {
            best = i;
        }
    }
    let v = eigen.eigenvectors.column(best);
    Vector::new(v[0], v[1], v[2])
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
    let normal = Direction::new(smallest_direction(m), tol).ok()?;
    Some(Canonical::Plane(Plane::new(Frame::about(
        Point::from_vector(c),
        normal,
    ))))
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

/// The axis of a surface of revolution, from its normals: every normal
/// line of such a surface lies in a plane through the axis, so meets it.
/// In Plücker coordinates a line `(d, m)` meets the normal line `(n, p×n)`
/// exactly when `d·(p×n) + m·n = 0`, which is linear; the axis is the
/// least-squares null vector of those conditions, put back on the Klein
/// quadric. A sphere's normals all meet at one point and determine no
/// axis, which the fits after this one catch.
fn revolution_axis(
    points: &[Point],
    normals: &[Vector],
    tol: Tolerances,
) -> Option<(Point, Direction)> {
    // Reweighted: a few samples off the surface (a flat face's corner met
    // along a tangent line) have normal lines far from the axis, and in
    // plain least squares those few decide it. Each round weighs a line by
    // how far it passes from the last round's axis, relative to the
    // typical miss.
    let mut weights = vec![1.0; points.len()];
    let mut axis = None;
    for _ in 0..8 {
        let found = weighted_axis(points, normals, &weights, tol)?;
        let (through, direction) = found;
        axis = Some(found);
        let a = direction.vector();
        let misses: Vec<f64> = points
            .iter()
            .zip(normals)
            .map(|(p, n)| {
                let w = *p - through;
                let across = a.cross(*n);
                let m = across.magnitude();
                if m > 1e-9 {
                    w.dot(across).abs() / m
                } else {
                    w.cross(a).magnitude()
                }
            })
            .collect();
        let mut sorted = misses.clone();
        sorted.sort_by(f64::total_cmp);
        let typical = sorted[sorted.len() / 2].max(tol.confusion());
        for (w, miss) in weights.iter_mut().zip(&misses) {
            let r = miss / (2.0 * typical);
            *w = 1.0 / r.mul_add(r, 1.0);
        }
    }
    axis
}

fn weighted_axis(
    points: &[Point],
    normals: &[Vector],
    weights: &[f64],
    tol: Tolerances,
) -> Option<(Point, Direction)> {
    // Minimized with the direction held to unit length, not the whole
    // six-vector: a cylinder's normals are all perpendicular to its axis,
    // and the direction-free "line at infinity" meets every one of them.
    // For a fixed direction the best moment is linear in it, so the
    // moment is eliminated and the direction is the Schur complement's
    // smallest eigenvector.
    let mut dd = nalgebra::Matrix3::<f64>::zeros();
    let mut dm = nalgebra::Matrix3::<f64>::zeros();
    let mut mm = nalgebra::Matrix3::<f64>::zeros();
    for ((p, n), w) in points.iter().zip(normals).zip(weights) {
        let moment = p.to_vector().cross(*n);
        let a = nalgebra::Vector3::new(moment.x, moment.y, moment.z);
        let b = nalgebra::Vector3::new(n.x, n.y, n.z);
        dd += a * a.transpose() * *w;
        dm += a * b.transpose() * *w;
        mm += b * b.transpose() * *w;
    }
    let scale = mm.norm().max(f64::MIN_POSITIVE);
    let pinv = mm.pseudo_inverse(scale * 1e-9).ok()?;
    let schur = dd - dm * pinv * dm.transpose();
    let schur = (schur + schur.transpose()) * 0.5;
    let d = smallest_direction(schur);
    let dv = nalgebra::Vector3::new(d.x, d.y, d.z);
    let mv = -(pinv * dm.transpose() * dv);
    let moment = Vector::new(mv[0], mv[1], mv[2]);
    let length = d.magnitude();
    if length <= f64::MIN_POSITIVE {
        return None;
    }
    let d = d / length;
    let moment = (moment - d * d.dot(moment)) / length;
    let through = Point::from_vector(d.cross(moment));
    Some((through, Direction::new(d, tol).ok()?))
}

/// Each point's height along the axis and distance from it: the profile
/// the surface of revolution sweeps.
fn profile(points: &[Point], through: Point, axis: Direction) -> Vec<(f64, f64)> {
    let a = axis.vector();
    points
        .iter()
        .map(|p| {
            let w = *p - through;
            let h = w.dot(a);
            ((w - a * h).magnitude(), h)
        })
        .collect()
}

/// The axis point at the samples' mean height, so the frame sits among
/// them rather than wherever the axis estimate happened to pass.
fn centred_on(points: &[Point], through: Point, axis: Direction) -> Point {
    let c = Point::from_vector(centroid(points));
    through + axis.vector() * (c - through).dot(axis.vector())
}

fn fit_cylinder(points: &[Point], normals: &[Vector], tol: Tolerances) -> Option<Canonical> {
    let (through, axis) = revolution_axis(points, normals, tol)?;
    let rows = profile(points, through, axis);
    #[allow(
        clippy::cast_precision_loss,
        reason = "sample counts are far below 2^52"
    )]
    let radius = rows.iter().map(|(rho, _)| rho).sum::<f64>() / rows.len() as f64;
    Some(Canonical::Cylinder(
        Cylinder::new(
            Frame::about(centred_on(points, through, axis), axis),
            radius,
            tol,
        )
        .ok()?,
    ))
}

fn fit_cone(points: &[Point], normals: &[Vector], tol: Tolerances) -> Option<Canonical> {
    let (through, axis) = revolution_axis(points, normals, tol)?;
    let origin = centred_on(points, through, axis);
    // Radius against height is a line: ρ = k·h + ρ₀.
    let (mut sh, mut shh, mut sr, mut shr, mut count) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for (rho, h) in profile(points, origin, axis) {
        sh += h;
        shh += h * h;
        sr += rho;
        shr += h * rho;
        count += 1.0;
    }
    let det = f64::mul_add(count, shh, -(sh * sh));
    if det.abs() <= f64::MIN_POSITIVE {
        return None;
    }
    let k = f64::mul_add(count, shr, -(sh * sr)) / det;
    let rho0 = f64::mul_add(shh, sr, -(sh * shr)) / det;
    if k.abs() <= tol.angular() {
        // No taper is a cylinder.
        return None;
    }
    // A negative taper is the same cone seen from its other end.
    let (axis, k) = if k < 0.0 { (-axis, -k) } else { (axis, k) };
    Some(Canonical::Cone(
        Cone::new(
            Frame::about(origin, axis),
            rho0.max(tol.confusion()),
            k.atan(),
            tol,
        )
        .ok()?,
    ))
}

/// A circle through points in space: its plane from their covariance, its
/// centre and radius from an algebraic fit in that plane, and the rms of
/// how far the points miss it.
fn circle_through(points: &[Point], tol: Tolerances) -> Option<(Point, Direction, f64, f64)> {
    let c = centroid(points);
    let normal = Direction::new(
        smallest_direction(covariance(points.iter().map(|p| p.to_vector() - c))),
        tol,
    )
    .ok()?;
    let a = normal.vector();
    let e1 = normal.any_perpendicular().vector();
    let e2 = a.cross(e1);
    let mut m = nalgebra::Matrix3::zeros();
    let mut b = nalgebra::Vector3::zeros();
    for p in points {
        let d = p.to_vector() - c;
        let (x, y) = (d.dot(e1), d.dot(e2));
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
    let radius = r2.sqrt();
    let centre = Point::from_vector(c + e1 * solved[0] + e2 * solved[1]);
    let mut miss = 0.0;
    for p in points {
        let d = *p - centre;
        let off = d.dot(a);
        let within = (d - a * off).magnitude() - radius;
        miss += off.mul_add(off, within * within);
    }
    #[allow(
        clippy::cast_precision_loss,
        reason = "sample counts are far below 2^52"
    )]
    let rms = (miss / points.len() as f64).sqrt();
    Some((centre, normal, radius, rms))
}

fn fit_torus(points: &[Point], normals: &[Vector], tol: Tolerances) -> Option<Canonical> {
    // Every normal line of a torus passes through its spine circle, a tube
    // radius from the surface: at the right radius the points shifted back
    // along their normals lie on one circle, whose plane gives the axis.
    // Searched log-spaced over both signs, since the normals may face
    // either way and the tube may be many times the patch's size.
    // The tube is the tighter of the torus's two curvatures, so its radius
    // is read from how fast the normals turn between samples: the largest
    // turn per unit distance, taken robustly as a high quantile over pairs.
    // The search then runs finely within a factor of two of it: the
    // spine's fit sharpens to a narrow valley at the true radius, which a
    // coarse sweep over every scale steps across.
    let stride = points.len().div_ceil(160).max(1);
    let picked: Vec<usize> = (0..points.len()).step_by(stride).collect();
    let mut turns: Vec<f64> = Vec::with_capacity(picked.len() * picked.len() / 2);
    for (a, &i) in picked.iter().enumerate() {
        for &j in &picked[a + 1..] {
            let gap = points[i].distance(points[j]);
            if gap > 0.0 {
                turns.push((normals[i] - normals[j]).magnitude() / gap);
            }
        }
    }
    if turns.is_empty() {
        return None;
    }
    turns.sort_by(f64::total_cmp);
    let sharpest = turns[turns.len() * 9 / 10];
    if sharpest <= 0.0 {
        return None;
    }
    let estimate = 1.0 / sharpest;
    let miss = |r: f64| -> Option<(Point, Direction, f64, f64)> {
        let shifted: Vec<Point> = points
            .iter()
            .zip(normals)
            .map(|(p, n)| *p - *n * r)
            .collect();
        circle_through(&shifted, tol)
    };
    // Only a spine wider than the tube makes a torus: shifted far enough,
    // the points crowd onto a small circle near the axis, which fits well
    // and describes nothing.
    let score = |r: f64| {
        miss(r)
            .filter(|m| m.2 > r.abs())
            .map_or(f64::INFINITY, |m| m.3)
    };
    let steps = 60_i32;
    let ratio = 4.0_f64.powf(1.0 / f64::from(steps));
    let mut best: Option<(f64, f64)> = None;
    for sign in [-1.0, 1.0] {
        for k in 0..=steps {
            let r = sign * estimate * 0.5 * ratio.powi(k);
            let rms = score(r);
            if rms.is_finite() && best.is_none_or(|(held, _)| rms < held) {
                best = Some((rms, r));
            }
        }
    }
    let (_, tube) = best?;
    let (mut lo, mut hi) = (tube / ratio, tube * ratio);
    if lo > hi {
        core::mem::swap(&mut lo, &mut hi);
    }
    for _ in 0..40 {
        let m1 = lo + (hi - lo) * 0.382;
        let m2 = lo + (hi - lo) * 0.618;
        if score(m1) < score(m2) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let tube = f64::midpoint(lo, hi);
    let (centre, axis, major, _) = miss(tube)?;
    let minor = tube.abs();
    if minor <= tol.confusion() || minor >= major {
        return None;
    }
    Some(Canonical::Torus(
        Torus::new(Frame::about(centre, axis), major, minor, tol).ok()?,
    ))
}

/// A surface's defining numbers, as the refinement moves them: a point, a
/// direction (not kept unit while it moves), and the kind's radii.
fn parameters(surface: &Canonical) -> Option<Vec<f64>> {
    let flat = |o: Point, d: Direction, rest: &[f64]| {
        let v = d.vector();
        let mut out = vec![o.x, o.y, o.z, v.x, v.y, v.z];
        out.extend_from_slice(rest);
        out
    };
    Some(match surface {
        Canonical::Plane(_) => return None,
        Canonical::Cylinder(c) => flat(c.frame().origin(), c.frame().z(), &[c.radius()]),
        Canonical::Cone(c) => flat(
            c.frame().origin(),
            c.frame().z(),
            &[c.reference_radius(), c.half_angle()],
        ),
        Canonical::Sphere(s) => {
            let o = s.centre();
            vec![o.x, o.y, o.z, s.radius()]
        }
        Canonical::Torus(t) => flat(
            t.frame().origin(),
            t.frame().z(),
            &[t.major_radius(), t.minor_radius()],
        ),
    })
}

/// A point's signed distance to the surface `x` describes, of the kind
/// `like` is.
fn residual(like: &Canonical, x: &[f64], p: Point) -> f64 {
    let o = Point::new(x[0], x[1], x[2]);
    if let Canonical::Sphere(_) = like {
        return p.distance(o) - x[3];
    }
    let d = Vector::new(x[3], x[4], x[5]);
    let m = d.magnitude();
    let a = if m > 0.0 { d / m } else { Vector::Z };
    let w = p - o;
    let h = w.dot(a);
    let rho = (w - a * h).magnitude();
    match like {
        Canonical::Cylinder(_) => rho - x[6],
        Canonical::Cone(_) => {
            let (r0, angle) = (x[6], x[7]);
            (rho - h.mul_add(angle.tan(), r0)) * angle.cos()
        }
        Canonical::Torus(_) => (rho - x[6]).hypot(h) - x[7],
        Canonical::Plane(_) | Canonical::Sphere(_) => 0.0,
    }
}

fn rebuild(like: &Canonical, x: &[f64], tol: Tolerances) -> Option<Canonical> {
    let o = Point::new(x[0], x[1], x[2]);
    if let Canonical::Sphere(_) = like {
        return Some(Canonical::Sphere(Sphere::centred(o, x[3], tol).ok()?));
    }
    let axis = Direction::new(Vector::new(x[3], x[4], x[5]), tol).ok()?;
    let frame = Frame::about(o, axis);
    Some(match like {
        Canonical::Cylinder(_) => Canonical::Cylinder(Cylinder::new(frame, x[6], tol).ok()?),
        Canonical::Cone(_) => {
            if x[7] <= 0.0 || x[7] >= core::f64::consts::FRAC_PI_2 {
                return None;
            }
            Canonical::Cone(Cone::new(frame, x[6].max(tol.confusion()), x[7], tol).ok()?)
        }
        Canonical::Torus(_) => {
            if x[7] <= 0.0 || x[7] >= x[6] {
                return None;
            }
            Canonical::Torus(Torus::new(frame, x[6], x[7], tol).ok()?)
        }
        Canonical::Plane(_) | Canonical::Sphere(_) => return None,
    })
}

/// The same surface with its axis direction unit and its origin the axis
/// point nearest the samples' centroid: the two freedoms the residual does
/// not see, pinned so the solve cannot wander along them.
fn regauged(like: &Canonical, mut x: Vec<f64>, points: &[Point]) -> Vec<f64> {
    if let Canonical::Sphere(_) = like {
        return x;
    }
    let d = Vector::new(x[3], x[4], x[5]);
    let m = d.magnitude();
    if m == 0.0 {
        return x;
    }
    let a = d / m;
    x[3..6].copy_from_slice(&[a.x, a.y, a.z]);
    // A torus's centre is a point, not a place along a line.
    if let Canonical::Torus(_) = like {
        return x;
    }
    let o = Point::new(x[0], x[1], x[2]);
    let c = Point::from_vector(centroid(points));
    let shift = (c - o).dot(a);
    let moved = o + a * shift;
    if let Canonical::Cone(_) = like {
        // The reference radius is the radius at the origin.
        x[6] = shift.mul_add(x[7].tan(), x[6]);
    }
    x[..3].copy_from_slice(&[moved.x, moved.y, moved.z]);
    x
}

/// Least squares on the points' distances to the surface, from the
/// estimate: Gauss–Newton steps through the Jacobian's singular value
/// decomposition, the directions the samples barely determine truncated
/// rather than amplified, and each step halved until it lowers the cost.
/// The Jacobian is by central differences: a handful of numbers, and a
/// residual cheap enough that exactness in it buys nothing.
fn refine(seed: Canonical, points: &[Point], hopeless: f64, tol: Tolerances) -> Option<Canonical> {
    let mut x = parameters(&seed)?;
    let n = x.len();
    let cost = |x: &[f64]| -> f64 { points.iter().map(|p| residual(&seed, x, *p).powi(2)).sum() };
    let scale = points
        .iter()
        .map(|p| p.distance(points[0]))
        .fold(tol.confusion(), f64::max);
    let mut current = cost(&x);
    #[allow(
        clippy::cast_precision_loss,
        reason = "sample counts are far below 2^52"
    )]
    let count = points.len() as f64;
    let mut slowed = false;
    // A torus on a small patch is ill-conditioned (its tube and its sweep
    // trade off against each other) and closes in slowly before it closes
    // in fast: it gets the steps, and is never given up for slowing.
    let patient = matches!(seed, Canonical::Torus(_));
    let (steps_allowed, least_gain) = if patient { (60, 1e-10) } else { (30, 1e-4) };
    for iteration in 0..steps_allowed {
        // Converging onto samples that are this surface, a step divides
        // the cost many times over; onto samples that are some other
        // surface, the steps stall at a floor. A fit that has stopped
        // halving its cost while still missing by more than the tolerance
        // will not get there.
        if !patient && iteration >= 3 && slowed && (current / count).sqrt() > hopeless {
            break;
        }
        let mut jacobian = nalgebra::DMatrix::<f64>::zeros(points.len(), n);
        let mut r = nalgebra::DVector::<f64>::zeros(points.len());
        let steps: Vec<f64> = x.iter().map(|v| 1e-7 * v.abs().max(scale)).collect();
        for (i, p) in points.iter().enumerate() {
            r[i] = residual(&seed, &x, *p);
            for k in 0..n {
                let mut up = x.clone();
                let mut down = x.clone();
                up[k] += steps[k];
                down[k] -= steps[k];
                jacobian[(i, k)] =
                    (residual(&seed, &up, *p) - residual(&seed, &down, *p)) / (2.0 * steps[k]);
            }
        }
        let svd = jacobian.svd(true, true);
        let largest = svd.singular_values.max();
        let Ok(step) = svd.solve(&(-r), largest * 1e-8) else {
            break;
        };
        let mut scale_step = 1.0;
        let mut accepted = false;
        for _ in 0..12 {
            let trial: Vec<f64> = x
                .iter()
                .zip(step.iter())
                .map(|(a, b)| b.mul_add(scale_step, *a))
                .collect();
            let trial = regauged(&seed, trial, points);
            let c = cost(&trial);
            if c.is_finite() && c < current {
                let gain = current - c;
                slowed = c > current * 0.5;
                x = trial;
                current = c;
                // Converged once a step gains less than a ten-thousandth:
                // exact samples fall by orders of magnitude a step until
                // they reach the rounding floor, and samples on no such
                // surface creep toward a miss the tolerance will refuse.
                accepted = gain > current * least_gain;
                break;
            }
            scale_step *= 0.5;
        }
        if !accepted {
            break;
        }
    }
    rebuild(&seed, &x, tol)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;

    const T: Tolerances = Tolerances::millimetres();

    /// Points on a cylinder with normals tipped a few degrees off true, as
    /// a mesh's averaged facet normals are: the refinement recovers the
    /// cylinder from the points to well within a micron.
    #[test]
    fn a_cylinder_is_recognized_from_rough_normals() {
        let axis = Direction::new(Vector::new(0.2, -0.1, 1.0), T).unwrap();
        let frame = Frame::about(Point::new(3.0, -2.0, 1.0), axis);
        let (mut points, mut normals) = (Vec::new(), Vec::new());
        for i in 0..24 {
            for j in 0..5 {
                let u = f64::from(i) * 0.2;
                let h = f64::from(j) * 2.0;
                let radial = frame.x().vector() * u.cos() + frame.y().vector() * u.sin();
                points.push(frame.origin() + radial * 7.5 + axis.vector() * h);
                let tipped = radial + axis.vector() * 0.05 * f64::from(j % 2);
                normals.push(tipped / tipped.magnitude());
            }
        }
        let found = recognize_points(&points, &normals, 1e-6, T)
            .unwrap()
            .unwrap();
        let Canonical::Cylinder(c) = found.surface else {
            panic!("a cylinder: {found:?}");
        };
        assert!((c.radius() - 7.5).abs() < 1e-7, "{c:?}");
    }

    /// A small patch of a thick torus (a few millimetres of a tube ten
    /// millimetres in radius, its normals tipped off true) is recognized
    /// as that torus, exactly: the tube radius lies far outside the patch's
    /// own size, and the fit has to find it anyway.
    #[test]
    fn a_small_patch_of_a_thick_torus_is_that_torus() {
        let torus = Torus::new(Frame::WORLD, 40.0, 10.0, T).unwrap();
        let (mut points, mut normals) = (Vec::new(), Vec::new());
        for i in 0..6 {
            for j in 0..5 {
                let (u, v) = (0.3 + f64::from(i) * 0.0126, 1.0 + f64::from(j) * 0.0314);
                let at = ogeom_math::elementary::torus_at(&torus, u, v);
                points.push(at.point);
                let n = at.du.cross(at.dv);
                normals.push(n / n.magnitude() + Vector::new(0.003, 0.0, 0.0));
            }
        }
        let found = recognize_points(&points, &normals, 1e-9, T)
            .unwrap()
            .unwrap();
        let Canonical::Torus(t) = found.surface else {
            panic!("a torus: {found:?}");
        };
        assert!((t.minor_radius() - 10.0).abs() < 1e-6, "{t:?}");
        assert!((t.major_radius() - 40.0).abs() < 1e-6, "{t:?}");
    }

    #[test]
    fn a_free_form_patch_refuses_every_canonical() {
        let (mut points, mut normals) = (Vec::new(), Vec::new());
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
                .is_none()
        );
    }
}
