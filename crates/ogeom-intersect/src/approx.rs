//! The approximation stage: a traced branch becomes curves.
//!
//! A traced branch is a polyline with a stated chord tolerance: honest, and
//! not what anything downstream wants to hold. An edge wants a curve in space;
//! a face wants that curve in its *own parameter space*, because splitting a
//! face happens there and a curve the face cannot express is a curve it cannot
//! be split along (`docs/DATA_MODEL.md` §6).
//!
//! So one branch becomes three fits sharing one tolerance: the 3D curve, and
//! one pcurve per surface, each fitted from the samples the tracer already
//! recorded. The tracer kept the parameters on both surfaces at every point
//! precisely for this moment; re-deriving them here would be a projection per
//! point, solving again what the marcher already solved.
//!
//! # The tolerance story, stated once
//!
//! The result's tolerance is a *sum of stated parts*, not a hope: the trace
//! sits within its chord tolerance of the true intersection, and the fit sits
//! within its own reported error of the trace. Both numbers are carried, and
//! the total is what an edge built on this curve must widen its tolerance to.
//! Nothing here rounds a miss up to a hit; a fit that could not reach its
//! target says so, and the caller decides whether the looser curve is usable.
//!
//! # Seams
//!
//! A branch crossing a periodic surface's seam has parameter samples that jump
//! by a period: the pcurve polyline tears even though the curve in space is
//! smooth. The samples are unwrapped before fitting: each step is folded to
//! the nearest image, so the pcurve runs continuously past the seam and may
//! legitimately leave `[0, 2π)`. That is what a pcurve on a periodic surface
//! is; folding it back would re-tear it.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{BSpline2d, BSplineCurve, Surface, SurfaceGeometry};
use ogeom_math::Point2;

use crate::march::Traced;

/// A branch of an intersection, as curves.
#[derive(Debug, Clone, PartialEq)]
pub struct IntersectionCurve {
    /// The curve in space.
    pub curve: BSplineCurve,
    /// The same curve in the first surface's parameter space.
    pub on_a: BSpline2d,
    /// And in the second's.
    pub on_b: BSpline2d,
    /// How far the *fits* may sit from the traced polyline, in millimetres.
    ///
    /// Measured: the curve against the trace's samples, each pcurve lifted
    /// against the curve, and, where the surfaces meet at a shallow angle,
    /// how far along them their crossing may sit from the curve. The
    /// distance to the true intersection adds the trace's own chord
    /// tolerance on top; both are stated so an edge built on this knows
    /// what to carry.
    pub fit_error: f64,
    /// Whether every fit met the tolerance it was asked for.
    pub met: bool,
    /// Whether the branch is a closed loop.
    pub closed: bool,
}

/// Fit one traced branch to curves, within `tolerance`.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the branch has
/// fewer than two points or the tolerance is not a positive distance.
pub fn approximate_branch(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    branch: &Traced,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<IntersectionCurve> {
    if branch.points.len() < 2 {
        ogeom_bail!(
            Construction,
            "a branch of {} points is not a curve",
            branch.points.len()
        );
    }

    // Marching correction can leave consecutive samples closer than the
    // rounding it converged within, and two samples at one chord-length
    // parameter are a knot span with no data in it: the fitting system
    // reports itself singular where the real defect is the duplicate. Thin
    // them here, where the trace's own step says what "too close" means.
    let mut points: Vec<ogeom_math::Point> = Vec::with_capacity(branch.points.len());
    let mut kept_a = Vec::with_capacity(branch.on_a.len());
    let mut kept_b = Vec::with_capacity(branch.on_b.len());
    // A sample is one point seen three ways, and where the three disagree
    // it is not data: through a point where the surfaces touch, the tracer
    // can report a step's position with its neighbour's parameters, and
    // the joint fit, asked to pass through both descriptions at once,
    // stalls a thousand times above its budget at that one sample.
    let agrees = |i: usize, p: &ogeom_math::Point| -> bool {
        let limit = tolerance.max(tol.confusion());
        let (ua, va) = branch.on_a[i];
        let (ub, vb) = branch.on_b[i];
        a.point_at(ua, va, tol)
            .is_ok_and(|q| q.distance(*p) <= limit)
            && b.point_at(ub, vb, tol)
                .is_ok_and(|q| q.distance(*p) <= limit)
    };
    for (i, p) in branch.points.iter().enumerate() {
        let end = i == 0 || i + 1 == branch.points.len();
        if let Some(last) = points.last()
            && last.distance(*p) <= tol.confusion() * 10.0
            && i + 1 != branch.points.len()
        {
            continue;
        }
        if !end && !agrees(i, p) {
            continue;
        }
        points.push(*p);
        kept_a.push(branch.on_a[i]);
        kept_b.push(branch.on_b[i]);
    }
    if points.len() < 2 {
        ogeom_bail!(Construction, "a branch of coincident points is not a curve");
    }

    // One fit in seven dimensions: the curve and both parameter images
    // together. Fitted separately, each fit's parameter correction drifts
    // its parameterization independently and the three results silently stop
    // being same-parameter: a pcurve claiming 1e-7 can evaluate millimetres
    // from its own curve. Jointly, one
    // parameterization and one knot vector serve all three, and the reported
    // error bounds every coordinate.
    let unwrapped_a = unwrap_periodic(a, &kept_a, tol);
    let unwrapped_b = unwrap_periodic(b, &kept_b, tol);
    // A closed branch takes the loop-smoothing fit: the join's tangents are
    // constrained to agree in all seven coordinates, so the section curve and
    // both pcurves cross their own seam without a crease. A loop winding
    // once round a periodic surface ends its image there a period from where
    // it began, and that fit, closing only where every image does, takes it
    // as an open trace whose ends lie close together, which can stall far
    // from the trace. Where it does, the loop is fitted closed in space with
    // the images' join C1 across the seam, right where the chart runs at one
    // speed across it, as a periodic surface's does (a patch that merely
    // meets itself at its seam need not), and the closer of the two stands.
    let winds_periodically = |surface: &SurfaceGeometry, image: &[Point2]| {
        let (first, last) = (image[0], image[image.len() - 1]);
        ((last.x - first.x).abs() <= tol.parametric() || surface.is_periodic_u())
            && ((last.y - first.y).abs() <= tol.parametric() || surface.is_periodic_v())
    };
    let free = || -> OgeomResult<Joint> {
        if branch.closed() {
            let closed = ogeom_geom::fit::fit_points_joint_closed(
                &points,
                &unwrapped_a,
                &unwrapped_b,
                3,
                tolerance,
                tol,
            )?;
            if !closed.0.met
                && winds_periodically(a, &unwrapped_a)
                && winds_periodically(b, &unwrapped_b)
            {
                let winding = ogeom_geom::fit::fit_points_joint_winding(
                    &points,
                    &unwrapped_a,
                    &unwrapped_b,
                    3,
                    tolerance,
                    tol,
                )?;
                if winding.0.error < closed.0.error {
                    return Ok(winding);
                }
            }
            Ok(closed)
        } else {
            ogeom_geom::fit::fit_points_joint(
                &points,
                &unwrapped_a,
                &unwrapped_b,
                3,
                tolerance,
                tol,
            )
        }
    };
    // The walk steps by how far its chord sags, which for a cubic is far
    // finer than the fit needs: between samples a step `h` apart on a curve
    // turning at `k`, a cubic misses by about `h^4 k^3 / 384` where the
    // chord sags by `h^2 k / 8`. So the knots start where that estimate
    // places a cubic's spans, every sample measures the fit, and only a
    // fit that misses one is made the free way as well, the closer of the
    // two standing.
    let stations = cubic_stations(&points, &unwrapped_a, &unwrapped_b, tolerance);
    let seeded = |closed: bool| {
        ogeom_geom::fit::fit_points_joint_from(
            &points,
            &unwrapped_a,
            &unwrapped_b,
            &stations,
            closed,
            3,
            tolerance,
            tol,
        )
    };
    let meets_itself = {
        let n = points.len();
        let (pa, pb) = (
            unwrapped_a[n - 1] - unwrapped_a[0],
            unwrapped_b[n - 1] - unwrapped_b[0],
        );
        let gap = (points[n - 1] - points[0]).square_magnitude()
            + pa.square_magnitude()
            + pb.square_magnitude();
        gap.sqrt() <= tol.confusion()
    };
    let mut fitted = seeded(branch.closed() && meets_itself).ok();
    if branch.closed()
        && !meets_itself
        && fitted.as_ref().is_none_or(|f| f.0.error > tolerance)
        && winds_periodically(a, &unwrapped_a)
        && winds_periodically(b, &unwrapped_b)
        && let Ok(winding) = seeded(true)
        && fitted.as_ref().is_none_or(|f| winding.0.error < f.0.error)
    {
        fitted = Some(winding);
    }
    let (space, on_a, on_b) = match fitted {
        Some(fitted) if fitted.0.error <= tolerance => fitted,
        Some(fitted) => {
            let other = free()?;
            if other.0.error < fitted.0.error {
                other
            } else {
                fitted
            }
        }
        None => free()?,
    };

    // Measured where it is promised, in millimetres. The joint fit's own
    // residual mixes space and chart coordinates and says little about
    // either alone, so it serves only as a bound on what is measured here.
    let lifted = lift_error(a, b, &on_a, &on_b, &space.curve, tol);
    // Where the surfaces meet at a shallow angle, a curve microns off both
    // can still sit far along them from where they cross. That distance is
    // the surfaces' gap over the sine of their angle, and it is no more than
    // the fit's chart residual carried through each surface's stretch, which
    // bounds how far the lifted pcurves stand from the trace whatever the
    // angle. The smaller of the two stands.
    let charted =
        space_error(a, &on_a, space.error, tol).max(space_error(b, &on_b, space.error, tol));
    let lifted = lifted
        .along
        .max(lifted.across.min(charted.max(space.error)));
    // Each sample's distance from the curve, which the fit's residual
    // already bounds: read only where that bound would decide the error.
    let fit_error = if space.error <= lifted {
        lifted
    } else {
        lifted.max(trace_error(&space.curve, &points, space.error, tol))
    };
    Ok(IntersectionCurve {
        fit_error,
        met: space.error <= tolerance,
        curve: space.curve,
        on_a,
        on_b,
        closed: branch.closed(),
    })
}

/// A joint fit: the curve, and its images on the two surfaces.
type Joint = (ogeom_geom::fit::Fitted<BSplineCurve>, BSpline2d, BSpline2d);

/// The samples a cubic fit's knots start at, by index, the first and last
/// among them: each stretch between two holds a cubic's estimated miss,
/// `L Θ³ / 384` for a stretch `L` long turning through `Θ`, to an eighth
/// of `tolerance`, in space and in either chart, the turn read off the
/// trace's own polylines.
fn cubic_stations(
    points: &[ogeom_math::Point],
    image_a: &[Point2],
    image_b: &[Point2],
    tolerance: f64,
) -> Vec<usize> {
    let n = points.len();
    let target = tolerance / 8.0;
    let traces: [Vec<[f64; 3]>; 3] = [
        points.iter().map(|p| [p.x, p.y, p.z]).collect(),
        image_a.iter().map(|p| [p.x, p.y, 0.0]).collect(),
        image_b.iter().map(|p| [p.x, p.y, 0.0]).collect(),
    ];
    // Per trace: each segment's length, and the turn at each sample
    // between two segments.
    let shape: Vec<(Vec<f64>, Vec<f64>)> = traces
        .iter()
        .map(|trace| {
            let segment =
                |k: usize| -> [f64; 3] { core::array::from_fn(|d| trace[k + 1][d] - trace[k][d]) };
            let dot = |u: [f64; 3], v: [f64; 3]| u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
            let lengths: Vec<f64> = (0..n - 1)
                .map(|k| dot(segment(k), segment(k)).sqrt())
                .collect();
            let mut turns = vec![0.0; n];
            for k in 1..n - 1 {
                let scale = lengths[k - 1] * lengths[k];
                if scale > 0.0 {
                    turns[k] = (dot(segment(k - 1), segment(k)) / scale)
                        .clamp(-1.0, 1.0)
                        .acos();
                }
            }
            (lengths, turns)
        })
        .collect();
    let mut out = vec![0];
    let mut from = 0;
    while from + 1 < n {
        let mut sums: Vec<(f64, f64)> = shape
            .iter()
            .map(|(lengths, _)| (lengths[from], 0.0))
            .collect();
        let mut to = from + 1;
        while to + 1 < n {
            let grown: Vec<(f64, f64)> = sums
                .iter()
                .zip(&shape)
                .map(|(&(length, turn), (lengths, turns))| (length + lengths[to], turn + turns[to]))
                .collect();
            if grown
                .iter()
                .any(|&(length, turn)| length * turn.powi(3) / 384.0 > target)
            {
                break;
            }
            sums = grown;
            to += 1;
        }
        out.push(to);
        from = to;
    }
    out
}

/// What lifting the pcurves onto their surfaces shows.
struct Lifted {
    /// How far either pcurve, lifted, stands from the curve at the same
    /// parameter.
    along: f64,
    /// That gap over the sine of the surfaces' angle there: how far the
    /// surfaces' true crossing may sit from the curve.
    across: f64,
}

/// Each pcurve lifted through its surface against the curve, at four
/// stations in every span of the curve's knots and at least two hundred
/// along it, between the trace's samples as well as at them. Near a cone's
/// apex or a sphere's pole the chart turns fast, and a fit that holds at
/// every sample can wander between them.
fn lift_error(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    on_a: &BSpline2d,
    on_b: &BSpline2d,
    curve: &BSplineCurve,
    tol: Tolerances,
) -> Lifted {
    use ogeom_geom::{Curve2d as _, Curve3d as _};
    let (lo, hi) = curve.knots().domain();
    let spans = curve.knots().distinct().len().saturating_sub(1).max(1);
    let stations = (4 * spans).max(200);
    let mut out = Lifted {
        along: 0.0,
        across: 0.0,
    };
    for k in 0..=stations {
        #[allow(clippy::cast_precision_loss)]
        let t = lo + (hi - lo) * k as f64 / stations as f64;
        let Ok(on) = curve.point_at(t, tol) else {
            continue;
        };
        let mut gap = 0.0_f64;
        let mut normals = Vec::with_capacity(2);
        for (surface, pcurve) in [(a, on_a), (b, on_b)] {
            let Ok(at) = pcurve.point_at(t, tol) else {
                continue;
            };
            let Ok(lifted) = surface.point_at(at.x, at.y, tol) else {
                continue;
            };
            gap = gap.max(lifted.distance(on));
            if let Ok(normal) = surface.normal_at(at.x, at.y, tol) {
                normals.push(normal.vector());
            }
        }
        out.along = out.along.max(gap);
        if let [na, nb] = normals[..] {
            let sine = na.cross(nb).magnitude().max(tol.angular());
            out.across = out.across.max(gap / sine);
        }
    }
    out
}

/// How far the curve stands from the trace it was fitted to: each sample's
/// distance to its nearest point on the curve.
///
/// The samples run in order along the curve, so each one's foot is found by
/// Newton steps from the last one's. Where that does not settle within
/// `bound`, the fit's own residual, which already bounds each sample's
/// distance from the curve at the fit's parameter, the nearest of the
/// curve's stations seeds the search instead, and the result never exceeds
/// `bound`.
fn trace_error(
    curve: &BSplineCurve,
    samples: &[ogeom_math::Point],
    bound: f64,
    tol: Tolerances,
) -> f64 {
    use ogeom_geom::Curve3d as _;
    let (lo, hi) = curve.knots().domain();
    let spans = curve.knots().distinct().len().saturating_sub(1).max(1);
    let count = (4 * spans).max(2 * samples.len()).max(200);
    #[allow(clippy::cast_precision_loss)]
    let at = |k: usize| lo + (hi - lo) * k as f64 / count as f64;
    let mut stations: Option<Vec<(f64, ogeom_math::Point)>> = None;
    let foot = |p: ogeom_math::Point, mut t: f64| -> (f64, f64) {
        let mut best = (t, f64::INFINITY);
        for _ in 0..8 {
            let (Ok(q), Ok(d)) = (curve.point_at(t, tol), curve.d1_at(t, tol)) else {
                break;
            };
            let gap = q.distance(p);
            if gap < best.1 {
                best = (t, gap);
            }
            let speed = d.dot(d);
            if speed <= f64::MIN_POSITIVE {
                break;
            }
            let next = (t + (p - q).dot(d) / speed).clamp(lo, hi);
            if (next - t).abs() <= (hi - lo) * 1e-12 {
                break;
            }
            t = next;
        }
        if let Ok(q) = curve.point_at(t, tol)
            && q.distance(p) < best.1
        {
            best = (t, q.distance(p));
        }
        best
    };
    let mut t = lo;
    let mut worst = 0.0_f64;
    for p in samples {
        let mut found = foot(*p, t);
        if found.1 > bound {
            let stations = stations.get_or_insert_with(|| {
                (0..=count)
                    .filter_map(|k| curve.point_at(at(k), tol).ok().map(|q| (at(k), q)))
                    .collect()
            });
            if let Some(k) = (0..stations.len()).min_by(|&x, &y| {
                stations[x]
                    .1
                    .distance(*p)
                    .total_cmp(&stations[y].1.distance(*p))
            }) {
                // Scanned finely over the stations either side, then
                // narrowed by golden section: where the curve all but stops
                // in space (its chart image swinging round a pole) or kinks
                // between two samples, Newton's steps are blind and the
                // distance has more than one dip between stations.
                let gap = |u: f64| {
                    curve
                        .point_at(u, tol)
                        .map_or(f64::INFINITY, |q| q.distance(*p))
                };
                let (from, to) = (
                    stations[k.saturating_sub(2)].0,
                    stations[(k + 2).min(stations.len() - 1)].0,
                );
                const FINE: u32 = 256;
                let h = (to - from) / f64::from(FINE);
                let start = (0..=FINE)
                    .map(|j| from + h * f64::from(j))
                    .min_by(|&x, &y| gap(x).total_cmp(&gap(y)))
                    .unwrap_or(from);
                let (mut a, mut b) = ((start - h).max(lo), (start + h).min(hi));
                let ratio = 0.5 * (5.0_f64.sqrt() - 1.0);
                for _ in 0..60 {
                    let (x, y) = (b - ratio * (b - a), a + ratio * (b - a));
                    if gap(x) <= gap(y) {
                        b = y;
                    } else {
                        a = x;
                    }
                }
                let again = foot(*p, 0.5 * (a + b));
                if again.1 < found.1 {
                    found = again;
                }
            }
        }
        t = found.0;
        worst = worst.max(found.1.min(bound));
    }
    worst
}

/// The fit's chart residual carried into space through the surface's
/// stretch along the pcurve: a bound on how far the lifted pcurve stands
/// from the trace, loose by the stretch where the residual is mostly in the
/// other coordinates, and so used only to cap an estimate, never stated.
fn space_error(
    surface: &SurfaceGeometry,
    pcurve: &BSpline2d,
    parameter_error: f64,
    tol: Tolerances,
) -> f64 {
    use ogeom_geom::Curve2d;
    // Convert the parameter-space error back through the surface's local
    // stretch at a few places; take the worst.
    let (lo, hi) = pcurve.domain();
    let mut worst = 0.0_f64;
    for i in 0..=16 {
        #[allow(clippy::cast_precision_loss)]
        let u = lo + (hi - lo) * f64::from(i) / 16.0;
        let Ok(at) = pcurve.point_at(u, tol) else {
            continue;
        };
        let Ok((du, dv)) = surface.d1_at(at.x, at.y, tol) else {
            continue;
        };
        let stretch = du.magnitude().max(dv.magnitude());
        worst = worst.max(parameter_error * stretch);
    }
    worst
}

/// Unfold parameter samples across a periodic surface's seam.
///
/// Each step is folded to the nearest image of the next sample, so a branch
/// crossing `u = 0` continues to `-0.1` rather than tearing to `2π - 0.1`. The
/// result may leave the surface's stated domain, which is what a pcurve
/// crossing a seam *is*.
fn unwrap_periodic(
    surface: &SurfaceGeometry,
    samples: &[(f64, f64)],
    tol: Tolerances,
) -> Vec<Point2> {
    let ((ua, ub), (va, vb)) = surface.domain();
    // Closure as well as periodicity: a converted drum is a clamped patch
    // that meets itself at its seam, and a loop walked round it lands on
    // either side of that seam by the walk's own rounding. Folded by the
    // chart's span like a period, the trace is the continuous curve it is.
    // Left as sampled, it jumps a whole span at the seam and the closed
    // fit chases the jump far from the trace.
    let u_period = if surface.is_periodic_u() || surface.is_closed_u(tol) {
        Some(ub - ua)
    } else {
        None
    };
    let v_period = if surface.is_periodic_v() || surface.is_closed_v(tol) {
        Some(vb - va)
    } else {
        None
    };
    let fold = |previous: f64, next: f64, period: Option<f64>| match period {
        None => next,
        Some(period) => {
            let mut candidate = next;
            while candidate - previous > period * 0.5 {
                candidate -= period;
            }
            while previous - candidate > period * 0.5 {
                candidate += period;
            }
            candidate
        }
    };

    let mut out = Vec::with_capacity(samples.len());
    let mut at = Point2::new(samples[0].0, samples[0].1);
    out.push(at);
    for sample in &samples[1..] {
        at = Point2::new(
            fold(at.x, sample.0, u_period),
            fold(at.y, sample.1, v_period),
        );
        out.push(at);
    }
    out
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use crate::march::{Marching, branches};
    use ogeom_geom::{Curve2d, Curve3d, CylinderSurface, PlaneSurface, SphereSurface};
    use ogeom_math::{Cylinder, Direction, Frame, Plane, Point, Sphere, Vector};

    const T: Tolerances = Tolerances::millimetres();

    fn sphere(radius: f64) -> SurfaceGeometry {
        SphereSurface::new(Sphere::centred(Point::ORIGIN, radius, T).unwrap()).into()
    }

    fn cylinder(radius: f64) -> SurfaceGeometry {
        CylinderSurface::new(Cylinder::new(Frame::WORLD, radius, T).unwrap(), (-4.0, 4.0))
            .unwrap()
            .into()
    }

    fn plane(origin: Point, normal: Vector) -> SurfaceGeometry {
        PlaneSurface::over(
            Plane::through(origin, Direction::new(normal, T).unwrap()),
            (-6.0, 6.0),
            (-6.0, 6.0),
        )
        .unwrap()
        .into()
    }

    fn options() -> Marching {
        Marching {
            chord: 1e-5,
            ..Marching::default()
        }
    }

    /// The distance of a fitted curve from both surfaces, sampled densely.
    ///
    /// This is the measure the whole stage exists for: the *fit* (not the
    /// polyline it came from) is what downstream code holds, so the fit is
    /// what must lie on both surfaces.
    fn fitted_deviation(a: &SurfaceGeometry, b: &SurfaceGeometry, curve: &BSplineCurve) -> f64 {
        let off = |surface: &SurfaceGeometry, p: Point| match surface {
            SurfaceGeometry::Plane(x) => x.plane().distance_to(p),
            SurfaceGeometry::Sphere(x) => x.sphere().distance_to(p),
            SurfaceGeometry::Cylinder(x) => x.cylinder().distance_to(p),
            _ => 0.0,
        };
        let (lo, hi) = curve.knots().domain();
        let mut worst = 0.0_f64;
        for i in 0..=800 {
            #[allow(clippy::cast_precision_loss)]
            let u = lo + (hi - lo) * f64::from(i) / 800.0;
            if let Ok(p) = curve.point_at(u, T) {
                worst = worst.max(off(a, p).abs().max(off(b, p).abs()));
            }
        }
        worst
    }

    #[test]
    fn a_fitted_branch_lies_on_both_surfaces_to_the_stated_total() {
        // The tolerance story end to end: trace within 1e-5, fit within 1e-4,
        // so the fitted curve is within the sum of the two of the true
        // intersection, measured against the surfaces, not the polyline.
        let a = sphere(3.0);
        let b = cylinder(1.5);
        let found = branches(&a, &b, options(), T).unwrap();
        assert_eq!(found.len(), 2);

        for branch in &found {
            let fitted = approximate_branch(&a, &b, branch, 1e-4, T).unwrap();
            assert!(fitted.met, "fit error {:e}", fitted.fit_error);
            assert!(fitted.closed);
            let off = fitted_deviation(&a, &b, &fitted.curve);
            assert!(
                off <= 1e-4 + 1e-5,
                "the fitted curve is {off:e} off the surfaces"
            );
            // And it is compact: a curve, not a decorated polyline.
            assert!(
                fitted.curve.control_points().len() * 4 < branch.points.len(),
                "{} control points for {} samples",
                fitted.curve.control_points().len(),
                branch.points.len()
            );
        }
    }

    #[test]
    fn the_pcurves_lift_back_onto_the_curve() {
        // A pcurve is only worth having if evaluating it and lifting through
        // its surface lands on the intersection. Checked through both
        // surfaces at matched ends and sampled interiors.
        let a = sphere(3.0);
        let b = cylinder(1.5);
        let found = branches(&a, &b, options(), T).unwrap();
        let branch = &found[0];
        let fitted = approximate_branch(&a, &b, branch, 1e-4, T).unwrap();

        for (surface, pcurve) in [(&a, &fitted.on_a), (&b, &fitted.on_b)] {
            let (lo, hi) = pcurve.domain();
            for i in 0..=200 {
                #[allow(clippy::cast_precision_loss)]
                let u = lo + (hi - lo) * f64::from(i) / 200.0;
                let at = pcurve.point_at(u, T).unwrap();
                let lifted = surface.point_at(at.x, at.y, T).unwrap();
                // The lifted point is on its own surface by construction; what
                // matters is that it is on the *other* one too, i.e. on the
                // intersection.
                let off = match (surface as &SurfaceGeometry, &a, &b) {
                    _ if core::ptr::eq(surface, &a) => match &b {
                        SurfaceGeometry::Cylinder(c) => c.cylinder().distance_to(lifted),
                        _ => 0.0,
                    },
                    _ => match &a {
                        SurfaceGeometry::Sphere(s) => s.sphere().distance_to(lifted),
                        _ => 0.0,
                    },
                };
                assert!(
                    off.abs() < 5e-4,
                    "a lifted pcurve point is {off:e} off the intersection"
                );
            }
        }
    }

    #[test]
    fn a_branch_across_the_seam_gets_a_continuous_pcurve() {
        // A plane through a cylinder's axis at an angle produces an ellipse
        // whose pcurve crosses the cylinder's u = 0 seam. Folded naively the
        // pcurve tears by 2π; unwrapped it runs smoothly and leaves the stated
        // domain, which is what crossing a seam means.
        let a = cylinder(2.0);
        let b = plane(Point::ORIGIN, Vector::new(0.0, 0.4, 1.0));
        let found = branches(&a, &b, options(), T).unwrap();
        assert_eq!(found.len(), 1, "an oblique plane cuts one ellipse");
        let fitted = approximate_branch(&a, &b, &found[0], 1e-4, T).unwrap();

        // Continuity: no two adjacent samples of the fitted pcurve jump by
        // anything near a period.
        let (lo, hi) = fitted.on_a.domain();
        let mut previous = fitted.on_a.point_at(lo, T).unwrap();
        for i in 1..=400 {
            #[allow(clippy::cast_precision_loss)]
            let u = lo + (hi - lo) * f64::from(i) / 400.0;
            let at = fitted.on_a.point_at(u, T).unwrap();
            assert!(
                (at.x - previous.x).abs() < 1.0,
                "the pcurve tears at the seam: {} to {}",
                previous.x,
                at.x
            );
            previous = at;
        }
    }

    /// A loop walked round a converted drum is closed, seam or no seam.
    ///
    /// A cylinder converted to a patch is clamped, not periodic: it meets
    /// itself at its seam. A plane across it cuts a circle the walk reaches
    /// the seam on from both sides, each half stopping a fraction of a step
    /// short of it, and the joined branch has coincident ends. Left flagged
    /// as having left the domain, the arrangement downstream would hold a
    /// circle with two ends at one point. It is closed, and fitted as a loop
    /// whose chart image runs continuously across the seam.
    #[test]
    fn a_loop_cut_at_a_converted_drum_s_seam_is_closed() {
        let drum: SurfaceGeometry = cylinder(2.0).to_bspline(T).unwrap().into();
        assert!(matches!(drum, SurfaceGeometry::BSpline(_)));
        let cut = plane(Point::new(0.0, 0.0, 1.0), Vector::new(0.0, 0.2, 1.0));
        let found = branches(&drum, &cut, options(), T).unwrap();
        assert_eq!(found.len(), 1, "an oblique plane cuts one loop");
        assert!(found[0].closed(), "the loop closes on the seam");
        let fitted = approximate_branch(&drum, &cut, &found[0], 1e-4, T).unwrap();
        assert!(fitted.closed);
        assert!(
            fitted.fit_error < 1e-3,
            "the loop fits as one: {}",
            fitted.fit_error
        );
        let (lo, hi) = fitted.on_a.domain();
        let mut previous = fitted.on_a.point_at(lo, T).unwrap();
        for i in 1..=400 {
            let u = lo + (hi - lo) * f64::from(i) / 400.0;
            let at = fitted.on_a.point_at(u, T).unwrap();
            assert!(
                (at.x - previous.x).abs() < 0.5,
                "the chart image tears at the seam: {} to {}",
                previous.x,
                at.x
            );
            previous = at;
        }
    }

    /// A bore across a converted drum states the error its curves have.
    ///
    /// The joint fit's residual mixes millimetres with the drum chart's
    /// coordinates, and carried through the drum's stretch it reads hundreds
    /// of times larger than the distance the fitted curves stand from the
    /// two surfaces. The surfaces cross steeply here, so that distance is
    /// what is stated, within a small multiple.
    #[test]
    fn a_section_across_a_wide_drum_states_its_measured_error() {
        let radius = 23.6;
        let wall = CylinderSurface::new(
            Cylinder::new(Frame::WORLD, radius, T).unwrap(),
            (-30.0, 30.0),
        )
        .unwrap();
        let drum: SurfaceGeometry = SurfaceGeometry::from(wall).to_bspline(T).unwrap().into();
        let bore = Cylinder::new(
            Frame::new(
                Point::new(0.0, 0.0, 3.0),
                Direction::new(Vector::X, T).unwrap(),
                Direction::new(Vector::Y, T).unwrap(),
                T,
            )
            .unwrap(),
            9.0,
            T,
        )
        .unwrap();
        let drill: SurfaceGeometry = CylinderSurface::new(bore, (-40.0, 40.0)).unwrap().into();
        let marching = Marching {
            chord: 1e-4,
            ..Marching::default()
        };
        let found = branches(&drum, &drill, marching, T).unwrap();
        assert!(!found.is_empty());
        for branch in &found {
            let fitted = approximate_branch(&drum, &drill, branch, 1e-4, T).unwrap();
            let (lo, hi) = fitted.curve.knots().domain();
            let mut off = 0.0_f64;
            for i in 0..=2000 {
                let t = lo + (hi - lo) * f64::from(i) / 2000.0;
                let p = fitted.curve.point_at(t, T).unwrap();
                off = off.max(
                    wall.cylinder()
                        .distance_to(p)
                        .abs()
                        .max(bore.distance_to(p).abs()),
                );
            }
            // Honest both ways: no less than the curve stands off, and not
            // hundreds of times more.
            assert!(
                fitted.fit_error + marching.chord >= off,
                "states {:e} for a curve {off:e} off its surfaces",
                fitted.fit_error
            );
            assert!(
                fitted.fit_error <= 10.0 * off.max(marching.chord),
                "states {:e} for a curve {off:e} off its surfaces",
                fitted.fit_error
            );
        }
    }

    #[test]
    fn what_cannot_be_fitted_is_refused() {
        let a = sphere(1.0);
        let b = plane(Point::ORIGIN, Vector::Z);
        let found = branches(&a, &b, options(), T).unwrap();
        assert!(approximate_branch(&a, &b, &found[0], 0.0, T).is_err());
        assert!(approximate_branch(&a, &b, &found[0], -1.0, T).is_err());

        let empty = Traced {
            points: vec![],
            on_a: vec![],
            on_b: vec![],
            stopped: crate::march::Stopped::Stalled,
        };
        assert!(approximate_branch(&a, &b, &empty, 1e-4, T).is_err());
    }
}
