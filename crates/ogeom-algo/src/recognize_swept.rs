//! Surfaces of revolution and extrusions read from a smooth region's
//! samples: a profile of any shape swept round an axis or along a
//! direction, for the regions none of the canonical surfaces fits.
//!
//! Each is held to the converter's standard: the profile is fitted to
//! within the tolerance, and every sample is then measured against the
//! swept surface itself. A profile that is straight sweeps a plane, a
//! cylinder or a cone, which the canonical recognition answers better, and
//! is refused here.

use ogeom_core::Tolerances;
use ogeom_geom::{Curve, ExtrusionSurface, RevolutionSurface, SurfaceGeometry};
use ogeom_math::{Axis, Direction, Point, Vector};

/// A swept surface and the worst distance from any sample to it.
#[derive(Debug, Clone)]
pub(crate) struct Swept {
    pub(crate) surface: SurfaceGeometry,
    pub(crate) deviation: f64,
}

/// The fewest distinct profile points a sweep is fitted through.
const PROFILE_POINTS: usize = 8;

/// An extrusion through the samples: their normals all square to one
/// direction, and their shadows along it on one curve.
pub(crate) fn fit_extrusion(
    points: &[Point],
    normals: &[Vector],
    tolerance: f64,
    tol: Tolerances,
) -> Option<Swept> {
    if points.len() < PROFILE_POINTS || points.len() != normals.len() {
        return None;
    }
    // The direction the normals are all square to: the smallest
    // eigenvector of their scatter.
    let mut scatter = nalgebra::Matrix3::<f64>::zeros();
    for n in normals {
        let v = nalgebra::Vector3::new(n.x, n.y, n.z);
        scatter += v * v.transpose();
    }
    let eigen = nalgebra::SymmetricEigen::new(scatter);
    let (smallest, _) = eigen
        .eigenvalues
        .iter()
        .enumerate()
        .min_by(|a, b| a.1.total_cmp(b.1))?;
    let d = eigen.eigenvectors.column(smallest);
    let direction = Direction::new(Vector::new(d[0], d[1], d[2]), tol).ok()?;
    let a = direction.vector();
    // Along the direction the samples must spread, or there is nothing
    // swept; across it the normals must agree with it to within a few
    // degrees, the rest left to the measurement.
    if normals.iter().any(|n| n.dot(a).abs() > 0.1) {
        return None;
    }
    let origin = points[0];
    let frame_x = {
        let any = if a.x.abs() < 0.9 {
            Vector::X
        } else {
            Vector::Y
        };
        let x = any - a * any.dot(a);
        x / x.magnitude()
    };
    let frame_y = a.cross(frame_x);
    let heights: Vec<f64> = points.iter().map(|p| (*p - origin).dot(a)).collect();
    let (lo, hi) = heights
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), &x| {
            (l.min(x), h.max(x))
        });
    if hi - lo <= tolerance * 10.0 {
        return None;
    }
    let shadows: Vec<(f64, f64)> = points
        .iter()
        .map(|p| {
            let w = *p - origin;
            (w.dot(frame_x), w.dot(frame_y))
        })
        .collect();
    let chain = ordered_profile(&shadows, tolerance)?;
    let lift = |(x, y): (f64, f64)| origin + frame_x * x + frame_y * y;
    let on: Vec<Point> = chain.iter().map(|&q| lift(q)).collect();
    if straight(&on, tolerance) {
        return None;
    }
    let curve = profile_curve(&on, tolerance, tol)?;
    let surface = ExtrusionSurface::over(curve.clone(), direction, (lo, hi)).ok()?;
    // Every sample against the curve, its shadow on the profile's plane.
    let mut deviation = 0.0_f64;
    for (p, h) in points.iter().zip(&heights) {
        let shadow = *p - a * *h;
        let foot = crate::project_on_curve(&curve, shadow, 64, tol).ok()?;
        deviation = deviation.max(foot.distance);
        if deviation > tolerance {
            return None;
        }
    }
    Some(Swept {
        surface: SurfaceGeometry::Extrusion(Box::new(surface)),
        deviation,
    })
}

/// A surface of revolution through the samples: their normals all meet one
/// axis, and their distances from it against their heights along it lie on
/// one profile.
pub(crate) fn fit_revolution(
    points: &[Point],
    normals: &[Vector],
    tolerance: f64,
    tol: Tolerances,
) -> Option<Swept> {
    if points.len() < PROFILE_POINTS || points.len() != normals.len() {
        return None;
    }
    let (through, axis) = crate::recognize::revolution_axis(points, normals, tol)?;
    let (through, axis) = refined_axis(points, through, axis, tol);
    let (through, axis) =
        settled_axis(points, through, axis, tolerance, tol).unwrap_or((through, axis));
    let a = axis.vector();
    // The profile lies in the half-plane through the axis and the first
    // sample, every sample turned into it.
    let first = points[0] - through;
    let radial = first - a * first.dot(a);
    let m = radial.magnitude();
    if m <= tolerance {
        return None;
    }
    let out = radial / m;
    let profile: Vec<(f64, f64)> = points
        .iter()
        .map(|p| {
            let w = *p - through;
            let h = w.dot(a);
            ((w - a * h).magnitude(), h)
        })
        .collect();
    let chain = ordered_profile(&profile, tolerance)?;
    let lift = |(r, h): (f64, f64)| through + out * r + a * h;
    let on: Vec<Point> = chain.iter().map(|&q| lift(q)).collect();
    if straight(&on, tolerance) {
        return None;
    }
    let curve = profile_curve(&on, tolerance, tol)?;
    // The turn the samples sweep, so the surface stands where they do: a
    // whole turn where they go round.
    let surface = RevolutionSurface::new(
        curve.clone(),
        Axis::new(through, axis),
        core::f64::consts::TAU,
    )
    .ok()?;
    let mut deviation = 0.0_f64;
    for &(r, h) in &profile {
        let foot = crate::project_on_curve(&curve, lift((r, h)), 64, tol).ok()?;
        deviation = deviation.max(foot.distance);
        if deviation > tolerance {
            return None;
        }
    }
    Some(Swept {
        surface: SurfaceGeometry::Revolution(Box::new(surface)),
        deviation,
    })
}

/// A curve through the profile's points: fitted to within half the
/// tolerance where that many control points suffice, and through every one
/// of them where they are too few to approximate. Either way the samples
/// are measured against the sweep afterwards.
fn profile_curve(points: &[Point], tolerance: f64, tol: Tolerances) -> Option<Curve> {
    if let Ok(fitted) = crate::fit::approximate_within(points, tolerance * 0.5, tol)
        && fitted.met
    {
        return Some(fitted.curve.into());
    }
    crate::fit::interpolate(points, 3, crate::fit::Spacing::Centripetal, tol)
        .ok()
        .map(Into::into)
}

/// The axis moved to where the points' distances from it agree best with
/// their heights along it, roughly.
///
/// The normals a mesh gives lean by up to its facets' turn, and an axis
/// read from them can be off by a good part of a facet's size on a coarse
/// mesh. On the true axis points at one height stand at one distance, so
/// the points are binned by height, each bin's distances fitted by a line
/// in height, and the axis's two tilts and two offsets moved by a simplex
/// search to make the misfit least. That brings the axis close from far
/// off, but the profile's own bend within a bin outweighs a small tilt, and
/// over a whole turn the tilt left is past the tolerance; [`settled_axis`]
/// takes it from there.
fn refined_axis(
    points: &[Point],
    through: Point,
    axis: Direction,
    tol: Tolerances,
) -> (Point, Direction) {
    let a = axis.vector();
    let e1 = {
        let any = if a.x.abs() < 0.9 {
            Vector::X
        } else {
            Vector::Y
        };
        let x = any - a * any.dot(a);
        x / x.magnitude()
    };
    let e2 = a.cross(e1);
    let centre = {
        let mut sum = Vector::ZERO;
        for p in points {
            sum += p.to_vector();
        }
        #[allow(clippy::cast_precision_loss, reason = "sample counts are small")]
        let n = points.len() as f64;
        Point::from_vector(sum / n)
    };
    let foot = through + a * (centre - through).dot(a);
    let reach = points
        .iter()
        .map(|p| p.distance(foot))
        .fold(0.0_f64, f64::max)
        .max(tol.confusion());
    // An axis from four numbers: tilts toward e1 and e2 (radians), and
    // offsets along them (as fractions of the reach).
    let place = |x: &[f64; 4]| -> Option<(Point, Direction)> {
        let d = a + e1 * x[0] + e2 * x[1];
        let origin = foot + (e1 * x[2] + e2 * x[3]) * reach;
        Some((origin, Direction::new(d, tol).ok()?))
    };
    let misfit = |x: &[f64; 4]| -> f64 {
        let Some((o, d)) = place(x) else {
            return f64::INFINITY;
        };
        let d = d.vector();
        let mut rh: Vec<(f64, f64)> = points
            .iter()
            .map(|p| {
                let w = *p - o;
                let h = w.dot(d);
                (h, (w - d * h).magnitude())
            })
            .collect();
        rh.sort_by(|p, q| p.0.total_cmp(&q.0));
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            clippy::cast_precision_loss,
            reason = "a bin count"
        )]
        let bins = ((rh.len() as f64).sqrt() as usize).max(1);
        let size = rh.len().div_ceil(bins).max(3);
        let mut sum = 0.0;
        for chunk in rh.chunks(size) {
            if chunk.len() < 3 {
                continue;
            }
            // Least squares r = c + s h over the bin.
            #[allow(clippy::cast_precision_loss, reason = "bin sizes are small")]
            let n = chunk.len() as f64;
            let (mh, mr) = chunk
                .iter()
                .fold((0.0, 0.0), |(x, y), &(h, r)| (x + h / n, y + r / n));
            let (mut shh, mut shr) = (0.0, 0.0);
            for &(h, r) in chunk {
                shh += (h - mh) * (h - mh);
                shr += (h - mh) * (r - mr);
            }
            let slope = if shh > 0.0 { shr / shh } else { 0.0 };
            for &(h, r) in chunk {
                let e = r - (mr + slope * (h - mh));
                sum += e * e;
            }
        }
        sum
    };
    let best = nelder_mead(&misfit, [0.0; 4], 1e-3, 400);
    place(&best).unwrap_or((through, axis))
}

/// A simplex search for the least of `f` from `start`, its first steps
/// `step` along each axis.
fn nelder_mead(
    f: &dyn Fn(&[f64; 4]) -> f64,
    start: [f64; 4],
    step: f64,
    rounds: usize,
) -> [f64; 4] {
    let mut simplex: Vec<([f64; 4], f64)> = Vec::with_capacity(5);
    simplex.push((start, f(&start)));
    for i in 0..4 {
        let mut x = start;
        x[i] += step;
        simplex.push((x, f(&x)));
    }
    let combine = |a: &[f64; 4], b: &[f64; 4], t: f64| -> [f64; 4] {
        let mut out = [0.0; 4];
        for k in 0..4 {
            out[k] = a[k] + (b[k] - a[k]) * t;
        }
        out
    };
    for _ in 0..rounds {
        simplex.sort_by(|p, q| p.1.total_cmp(&q.1));
        let mut centroid = [0.0; 4];
        for (x, _) in &simplex[..4] {
            for k in 0..4 {
                centroid[k] += x[k] / 4.0;
            }
        }
        let worst = simplex[4];
        let reflected = combine(&centroid, &worst.0, -1.0);
        let fr = f(&reflected);
        if fr < simplex[0].1 {
            let expanded = combine(&centroid, &worst.0, -2.0);
            let fe = f(&expanded);
            simplex[4] = if fe < fr {
                (expanded, fe)
            } else {
                (reflected, fr)
            };
        } else if fr < simplex[3].1 {
            simplex[4] = (reflected, fr);
        } else {
            let contracted = combine(&centroid, &worst.0, 0.5);
            let fc = f(&contracted);
            if fc < worst.1 {
                simplex[4] = (contracted, fc);
            } else {
                let best = simplex[0].0;
                for entry in simplex.iter_mut().skip(1) {
                    let x = combine(&best, &entry.0, 0.5);
                    *entry = (x, f(&x));
                }
            }
        }
        let spread = simplex[4].1 - simplex[0].1;
        if spread.abs() <= simplex[0].1.abs() * 1e-12 + f64::MIN_POSITIVE {
            break;
        }
    }
    simplex.sort_by(|p, q| p.1.total_cmp(&q.1));
    simplex[0].0
}

/// The axis moved until the samples lie on one profile turned about it.
///
/// The samples lie on the surface, so near the axis it is settled from
/// them. Each round fits the profile through the samples about the current
/// axis and measures each sample's offset from it, across the profile. One
/// linear least-squares system then solves for the axis's two offsets and
/// two tilts together with a correction to the profile, piecewise linear
/// along it. A tilt or an offset moves a sample across the profile by an
/// amount that turns with its angle about the axis, which no change of the
/// profile follows, so the system tells the two apart over part of a turn
/// as well as over a whole one. The rounds stop when the step is
/// negligible.
///
/// It needs a start near the axis. About an axis that is off, the samples
/// spread into a band, the profile is threaded back and forth through it,
/// and the offsets measured from it are noise. `None` where the profile
/// cannot be fitted or the rounds do not settle.
fn settled_axis(
    points: &[Point],
    through: Point,
    axis: Direction,
    tolerance: f64,
    tol: Tolerances,
) -> Option<(Point, Direction)> {
    let centre = {
        let mut sum = Vector::ZERO;
        for p in points {
            sum += p.to_vector();
        }
        #[allow(clippy::cast_precision_loss, reason = "sample counts are small")]
        let n = points.len() as f64;
        Point::from_vector(sum / n)
    };
    let (mut origin, mut direction) = (through, axis);
    for _ in 0..8 {
        let a = direction.vector();
        let e1 = {
            let any = if a.x.abs() < 0.9 {
                Vector::X
            } else {
                Vector::Y
            };
            let x = any - a * any.dot(a);
            x / x.magnitude()
        };
        let e2 = a.cross(e1);
        let foot = origin + a * (centre - origin).dot(a);
        // Each sample's distance from the axis, height along it and angle
        // round it.
        let polar: Vec<(f64, f64, f64)> = points
            .iter()
            .map(|p| {
                let w = *p - foot;
                let (x, y) = (w.dot(e1), w.dot(e2));
                (x.hypot(y), w.dot(a), y.atan2(x))
            })
            .collect();
        let profile: Vec<(f64, f64)> = polar.iter().map(|&(r, h, _)| (r, h)).collect();
        // About an axis still a little off, each parallel's samples spread
        // over a little more than the tolerance; they are merged coarser
        // until they chain.
        let chain = [1.0, 4.0, 16.0, 64.0]
            .iter()
            .find_map(|&k| ordered_profile(&profile, tolerance * k))?;
        let on: Vec<Point> = chain.iter().map(|&(r, h)| Point::new(r, h, 0.0)).collect();
        let curve = profile_curve(&on, tolerance, tol)?;
        let (u0, u1) = ogeom_geom::Curve3d::domain(&curve);
        let knots = (chain.len() / 2).clamp(4, 40);
        let unknowns = 4 + knots;
        let mut ata = nalgebra::DMatrix::<f64>::zeros(unknowns, unknowns);
        let mut atb = nalgebra::DVector::<f64>::zeros(unknowns);
        for &(r, h, angle) in &polar {
            let foot = crate::project_on_curve(&curve, Point::new(r, h, 0.0), 64, tol).ok()?;
            let tangent = ogeom_geom::Curve3d::d1_at(&curve, foot.parameter, tol).ok()?;
            let length = tangent.magnitude();
            if length <= 0.0 {
                continue;
            }
            // The profile's normal, and the sample's offset along it.
            let (nr, nh) = (-tangent.y / length, tangent.x / length);
            let offset = (r - foot.point.x) * nr + (h - foot.point.y) * nh;
            // Moving the axis by (d1, d2) and tilting it by (t1, t2) moves
            // the sample's distance by -(d + h t) along its own direction
            // round the axis, and its height by r t along it.
            let (c, s) = (angle.cos(), angle.sin());
            let mut row = nalgebra::DVector::<f64>::zeros(unknowns);
            row[0] = -c * nr;
            row[1] = -s * nr;
            row[2] = c * (r * nh - h * nr);
            row[3] = s * (r * nh - h * nr);
            #[allow(clippy::cast_precision_loss, reason = "a knot count")]
            let at = (foot.parameter - u0) / (u1 - u0) * (knots - 1) as f64;
            #[allow(
                clippy::cast_possible_truncation,
                clippy::cast_sign_loss,
                reason = "a knot index, clamped"
            )]
            let k = (at.floor().max(0.0) as usize).min(knots - 2);
            #[allow(clippy::cast_precision_loss, reason = "a knot index")]
            let within = at - k as f64;
            row[4 + k] = 1.0 - within;
            row[5 + k] = within;
            ata += &row * row.transpose();
            atb -= &row * offset;
        }
        // A trace of damping keeps a knot no sample reaches from leaving
        // the system singular.
        for k in 0..unknowns {
            ata[(k, k)] += 1e-12 * (1.0 + ata[(k, k)]);
        }
        let step = ata.lu().solve(&atb)?;
        origin = foot + e1 * step[0] + e2 * step[1];
        direction = Direction::new(a + e1 * step[2] + e2 * step[3], tol).ok()?;
        let reach = points
            .iter()
            .map(|p| p.distance(foot))
            .fold(0.0_f64, f64::max);
        let moved = step[0].hypot(step[1]) + step[2].hypot(step[3]) * reach;
        if moved <= tol.confusion() * 1e-3 {
            return Some((origin, direction));
        }
    }
    None
}

/// Whether the points lie within `tolerance` of the line through the two
/// farthest apart of them.
fn straight(points: &[Point], tolerance: f64) -> bool {
    let (Some(&a), Some(&b)) = (points.first(), points.last()) else {
        return true;
    };
    let along = b - a;
    let length = along.magnitude();
    if length <= tolerance {
        return true;
    }
    let unit = along / length;
    points.iter().all(|p| {
        let w = *p - a;
        (w - unit * w.dot(unit)).magnitude() <= tolerance
    })
}

/// The profile's points in order along it: points within `merge` of one
/// another taken as one, then chained from the point farthest along the
/// points' longest spread, each to its nearest unchained neighbour. `None`
/// where the chain jumps (the points are not on one open curve) or too few
/// points remain to fit through.
fn ordered_profile(points: &[(f64, f64)], merge: f64) -> Option<Vec<(f64, f64)>> {
    // Merged on a grid of the merging distance, each cell's mean.
    let mut cells: std::collections::HashMap<(i64, i64), (f64, f64, f64)> =
        std::collections::HashMap::new();
    #[allow(clippy::cast_possible_truncation, reason = "a grid cell")]
    for &(x, y) in points {
        let key = ((x / merge).floor() as i64, (y / merge).floor() as i64);
        let cell = cells.entry(key).or_insert((0.0, 0.0, 0.0));
        cell.0 += x;
        cell.1 += y;
        cell.2 += 1.0;
    }
    let mut merged: Vec<(f64, f64)> = cells.values().map(|&(x, y, n)| (x / n, y / n)).collect();
    merged.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    if merged.len() < PROFILE_POINTS {
        return None;
    }
    // The longest spread: the pair of extremes along x or y, the wider.
    let spread = |f: &dyn Fn(&(f64, f64)) -> f64| {
        let lo = merged.iter().map(f).fold(f64::INFINITY, f64::min);
        let hi = merged.iter().map(f).fold(f64::NEG_INFINITY, f64::max);
        hi - lo
    };
    let along_x = spread(&|p| p.0) >= spread(&|p| p.1);
    let key = |p: &(f64, f64)| if along_x { p.0 } else { p.1 };
    let start = merged
        .iter()
        .enumerate()
        .min_by(|a, b| key(a.1).total_cmp(&key(b.1)))?
        .0;
    let mut chain = vec![merged.swap_remove(start)];
    let mut steps = Vec::new();
    while !merged.is_empty() {
        let last = *chain.last()?;
        let (next, step) = merged
            .iter()
            .enumerate()
            .map(|(i, p)| (i, (p.0 - last.0).hypot(p.1 - last.1)))
            .min_by(|a, b| a.1.total_cmp(&b.1))?;
        let point = merged.swap_remove(next);
        // Two cells either side of a cell wall can hold one point.
        if step > merge {
            steps.push(step);
            chain.push(point);
        }
    }
    if chain.len() < PROFILE_POINTS {
        return None;
    }
    // A jump far beyond the typical step is two pieces, or a loop chained
    // across.
    let mut sorted = steps.clone();
    sorted.sort_by(f64::total_cmp);
    let typical = sorted[sorted.len() / 2].max(merge);
    if steps.iter().any(|&s| s > typical * 6.0) {
        return None;
    }
    Some(chain)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code")]
mod tests {
    use super::*;
    use ogeom_geom::Surface as _;

    const T: Tolerances = Tolerances::millimetres();

    /// Samples of a surface on a grid of its parameters, with its normals.
    fn sampled(
        surface: &SurfaceGeometry,
        u: (f64, f64),
        v: (f64, f64),
    ) -> (Vec<Point>, Vec<Vector>) {
        let (mut points, mut normals) = (Vec::new(), Vec::new());
        for i in 0..=24 {
            for j in 0..=12 {
                let s = u.0 + (u.1 - u.0) * f64::from(i) / 24.0;
                let t = v.0 + (v.1 - v.0) * f64::from(j) / 12.0;
                points.push(surface.point_at(s, t, T).unwrap());
                normals.push(surface.normal_at(s, t, T).unwrap().vector());
            }
        }
        (points, normals)
    }

    /// A wavy profile: a quarter of a sine over ten.
    fn wave() -> Curve {
        let pts: Vec<Point> = (0..=20)
            .map(|k| {
                let x = f64::from(k) * 0.5;
                Point::new(x, 3.0 + (x * 0.4).sin(), 0.0)
            })
            .collect();
        crate::fit::interpolate(&pts, 3, crate::fit::Spacing::Centripetal, T)
            .unwrap()
            .into()
    }

    /// A wavy profile extruded comes back an extrusion along its direction,
    /// every sample on it.
    #[test]
    fn a_wavy_wall_is_read_as_an_extrusion() {
        let wall = SurfaceGeometry::Extrusion(Box::new(
            ExtrusionSurface::new(wave(), Direction::Z, 5.0).unwrap(),
        ));
        let (u, v) = wall.domain();
        let (points, normals) = sampled(&wall, u, (v.0, v.1));
        let found = fit_extrusion(&points, &normals, 1e-4, T).unwrap();
        assert!(found.deviation <= 1e-4);
        let SurfaceGeometry::Extrusion(e) = &found.surface else {
            panic!("{:?}", found.surface);
        };
        assert!(e.direction().vector().cross(Vector::Z).magnitude() < 1e-6);
        assert!(fit_revolution(&points, &normals, 1e-4, T).is_none());
    }

    /// A wavy profile turned about an axis comes back a surface of
    /// revolution about it, every sample on it.
    #[test]
    fn a_turned_wave_is_read_as_a_revolution() {
        let axis = Axis::new(Point::ORIGIN, Direction::X);
        let turned = SurfaceGeometry::Revolution(Box::new(
            RevolutionSurface::new(wave(), axis, core::f64::consts::PI).unwrap(),
        ));
        let (u, v) = turned.domain();
        let (points, normals) = sampled(&turned, (u.0, u.0 + core::f64::consts::PI), v);
        let found = fit_revolution(&points, &normals, 1e-4, T).unwrap();
        assert!(found.deviation <= 1e-4);
        assert!(fit_extrusion(&points, &normals, 1e-4, T).is_none());
    }

    /// A cylinder's straight profile is the canonical recognition's, and
    /// neither sweep takes it.
    #[test]
    fn a_straight_profile_is_left_to_the_canonical_surfaces() {
        let line: Curve = ogeom_geom::LineCurve::segment(
            Point::new(0.0, 3.0, 0.0),
            Point::new(10.0, 3.0, 0.0),
            T,
        )
        .unwrap()
        .into();
        let drum = SurfaceGeometry::Revolution(Box::new(
            RevolutionSurface::new(
                line,
                Axis::new(Point::ORIGIN, Direction::X),
                core::f64::consts::PI,
            )
            .unwrap(),
        ));
        let (u, v) = drum.domain();
        let (points, normals) = sampled(&drum, (u.0, u.0 + core::f64::consts::PI), v);
        assert!(fit_revolution(&points, &normals, 1e-4, T).is_none());
    }
}
