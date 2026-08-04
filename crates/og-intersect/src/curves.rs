//! Where two curves meet, in the plane and in space.
//!
//! *Elsewhere* these are `Geom2dAPI_InterCurveCurve` and `IntCurve` for the
//! plane, and extrema-based crossing for space. The planar case is the
//! load-bearing one: boolean face splitting happens in a surface's parameter
//! space, and the curves it splits with are pcurves — so 2D curve/curve is the
//! operation the whole §8 pipeline stands on.
//!
//! # Two curves in space generically miss
//!
//! In the plane, two curves that cross, cross. In space they pass by: a
//! crossing is two points closer than a tolerance, not an exact common point,
//! and pretending otherwise would make every 3D result empty. So the 3D
//! answer reports the *gap* it achieved at each crossing, and the caller's
//! tolerance decides what counts. The 2D answer reports gaps too — a solved
//! crossing is still a pair of floats — but there the gap is rounding, not
//! geometry.
//!
//! # Overlap is an answer, not a failure
//!
//! Two collinear lines, two arcs of one circle: where the supports coincide,
//! "the intersection points" do not exist — the intersection is a stretch of
//! curve. That is reported as an overlap with the parameter ranges involved.
//! Detected for the analytic same-support cases; two B-splines that happen to
//! trace the same path are *not* detected as overlapping, and that limit is
//! recorded rather than discovered.
//!
//! # The general path is honest about resolution
//!
//! Non-analytic pairs are seeded by sampling both curves into segments and
//! testing the pairs, then polished by Newton onto the true crossing. Like the
//! surface seeding it mirrors, it finds what the sampling resolves: two
//! crossings closer together than a sample step can read as one. The sampling
//! density is a stated knob, not a hidden constant.

use og_core::{OgResult, Tolerances, og_bail};
use og_geom::{Curve, Curve2d, Curve3d, PlanarCurve};
use og_math::{Point, Point2, solve};

/// One crossing of two curves.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Crossing<P> {
    /// The parameter on the first curve.
    pub on_a: f64,
    /// The parameter on the second.
    pub on_b: f64,
    /// Where, taken from the first curve.
    pub point: P,
    /// How far apart the two curves are there.
    ///
    /// Rounding for a planar crossing; real geometry for a spatial one, where
    /// two curves generically miss and "crossing" means passing within the
    /// caller's tolerance.
    pub gap: f64,
}

/// A stretch where two curves share their support.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Overlap {
    /// The parameter range on the first curve.
    pub on_a: (f64, f64),
    /// The corresponding range on the second.
    pub on_b: (f64, f64),
}

/// What two curves do to each other.
#[derive(Debug, Clone, PartialEq)]
pub struct CurveIntersection<P> {
    /// Isolated crossings, in order along the first curve.
    pub crossings: Vec<Crossing<P>>,
    /// Stretches of shared support.
    ///
    /// Only the analytic same-support cases are detected — collinear lines,
    /// arcs of one circle. Two free-form curves tracing the same path come
    /// back as whatever isolated crossings the sampling finds, and that limit
    /// is stated here rather than discovered downstream.
    pub overlaps: Vec<Overlap>,
}

impl<P> CurveIntersection<P> {
    /// No contact at all.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.crossings.is_empty() && self.overlaps.is_empty()
    }

    const fn empty() -> Self {
        Self {
            crossings: Vec::new(),
            overlaps: Vec::new(),
        }
    }
}

/// How hard the general path looks.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CurveCurveOptions {
    /// How many segments each curve is sampled into when seeding.
    ///
    /// The resolution knob: two crossings inside one segment read as one.
    pub samples: usize,
    /// The widest gap that still counts as a crossing, in space.
    ///
    /// Meaningful for 3D, where curves generically miss. In 2D a genuine
    /// crossing converges to rounding and this only rejects near-misses.
    pub gap: f64,
}

impl Default for CurveCurveOptions {
    fn default() -> Self {
        Self {
            samples: 128,
            gap: 1e-7,
        }
    }
}

/// Where two planar curves meet.
///
/// Analytic pairs — lines and circles — are answered in closed form, overlaps
/// included. Everything else goes through sampling and Newton.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the options
/// are unusable.
pub fn intersect_curves_2d(
    a: &PlanarCurve,
    b: &PlanarCurve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> OgResult<CurveIntersection<Point2>> {
    check(options)?;
    match (a, b) {
        (PlanarCurve::Line(x), PlanarCurve::Line(y)) => Ok(line_line_2d(x, y, tol)),
        (PlanarCurve::Line(x), PlanarCurve::Circle(y)) => Ok(line_circle_2d(x, y, false, tol)),
        (PlanarCurve::Circle(x), PlanarCurve::Line(y)) => Ok(line_circle_2d(y, x, true, tol)),
        (PlanarCurve::Circle(x), PlanarCurve::Circle(y)) => Ok(circle_circle_2d(x, y, tol)),
        _ => general_2d(a, b, options, tol),
    }
}

/// Where two space curves pass within `options.gap` of each other.
///
/// # Errors
///
/// As [`intersect_curves_2d`].
pub fn intersect_curves(
    a: &Curve,
    b: &Curve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> OgResult<CurveIntersection<Point>> {
    check(options)?;
    if let (Curve::Line(x), Curve::Line(y)) = (a, b) {
        return Ok(line_line_3d(x, y, options, tol));
    }
    general_3d(a, b, options, tol)
}

fn check(options: CurveCurveOptions) -> OgResult<()> {
    if options.samples < 2 {
        og_bail!(Construction, "seeding needs at least two segments");
    }
    if !options.gap.is_finite() || options.gap <= 0.0 {
        og_bail!(Construction, "a gap of {} is not a distance", options.gap);
    }
    Ok(())
}

// --- analytic, planar --------------------------------------------------------

fn line_line_2d(
    a: &og_geom::Line2d,
    b: &og_geom::Line2d,
    tol: Tolerances,
) -> CurveIntersection<Point2> {
    let (oa, da) = (a.axis().location, a.axis().direction.vector());
    let (ob, db) = (b.axis().location, b.axis().direction.vector());
    let cross = da.cross(db);

    if cross.abs() <= tol.angular() {
        // Parallel. Collinear if one origin is on the other line.
        let between = ob - oa;
        if between.cross(da).abs() > tol.confusion() {
            return CurveIntersection::empty();
        }
        // The shared stretch, as each line's own parameter range.
        let (a_lo, a_hi) = a.domain();
        let (b_lo, b_hi) = b.domain();
        // Where b's range lands on a's parameter: t_a = (p - oa)·da.
        let project = |p: Point2| (p - oa).dot(da);
        let (s0, s1) = (project(ob + db * b_lo), project(ob + db * b_hi));
        let (lo, hi) = (s0.min(s1).max(a_lo), s0.max(s1).min(a_hi));
        if lo >= hi {
            return CurveIntersection::empty();
        }
        // And back onto b.
        let back = |t: f64| (oa + da * t - ob).dot(db);
        return CurveIntersection {
            crossings: Vec::new(),
            overlaps: vec![Overlap {
                on_a: (lo, hi),
                on_b: order(back(lo), back(hi)),
            }],
        };
    }

    let between = ob - oa;
    let t = between.cross(db) / cross;
    let s = between.cross(da) / cross;
    let (a_lo, a_hi) = a.domain();
    let (b_lo, b_hi) = b.domain();
    if t < a_lo - tol.parametric()
        || t > a_hi + tol.parametric()
        || s < b_lo - tol.parametric()
        || s > b_hi + tol.parametric()
    {
        return CurveIntersection::empty();
    }
    CurveIntersection {
        crossings: vec![Crossing {
            on_a: t,
            on_b: s,
            point: oa + da * t,
            gap: 0.0,
        }],
        overlaps: Vec::new(),
    }
}

fn line_circle_2d(
    line: &og_geom::Line2d,
    circle: &og_geom::Circle2d,
    swapped: bool,
    tol: Tolerances,
) -> CurveIntersection<Point2> {
    let (o, d) = (line.axis().location, line.axis().direction.vector());
    let c = circle.circle();
    let centre = c.centre();
    let radius = c.radius();

    // Foot of the perpendicular from the centre onto the line.
    let along = (centre - o).dot(d);
    let foot = o + d * along;
    let gap = foot.distance(centre);
    if gap > radius + tol.confusion() {
        return CurveIntersection::empty();
    }
    let half = radius.mul_add(radius, -(gap * gap)).max(0.0).sqrt();
    let candidates = if half <= tol.confusion() {
        vec![along]
    } else {
        vec![along - half, along + half]
    };

    let (l_lo, l_hi) = line.domain();
    let mut crossings = Vec::new();
    for t in candidates {
        if t < l_lo - tol.parametric() || t > l_hi + tol.parametric() {
            continue;
        }
        let p = o + d * t;
        let Some(s) = circle_parameter(circle, p, tol) else {
            continue;
        };
        let (on_a, on_b) = if swapped { (s, t) } else { (t, s) };
        crossings.push(Crossing {
            on_a,
            on_b,
            point: p,
            gap: 0.0,
        });
    }
    sort_crossings(&mut crossings);
    CurveIntersection {
        crossings,
        overlaps: Vec::new(),
    }
}

fn circle_circle_2d(
    a: &og_geom::Circle2d,
    b: &og_geom::Circle2d,
    tol: Tolerances,
) -> CurveIntersection<Point2> {
    let (ca, cb) = (a.circle(), b.circle());
    let between = cb.centre() - ca.centre();
    let distance = between.magnitude();
    let (ra, rb) = (ca.radius(), cb.radius());

    if distance <= tol.confusion() {
        if (ra - rb).abs() <= tol.confusion() {
            // The same circle: the overlap is both whole domains.
            return CurveIntersection {
                crossings: Vec::new(),
                overlaps: vec![Overlap {
                    on_a: a.domain(),
                    on_b: b.domain(),
                }],
            };
        }
        return CurveIntersection::empty();
    }
    if distance > ra + rb + tol.confusion() || distance < (ra - rb).abs() - tol.confusion() {
        return CurveIntersection::empty();
    }

    // The radical line: where the two circles' equations agree.
    let along = distance.mul_add(distance, ra.mul_add(ra, -(rb * rb))) / (2.0 * distance);
    let squared = ra.mul_add(ra, -(along * along));
    let direction = between * (1.0 / distance);
    let foot = ca.centre() + direction * along;
    let mut crossings = Vec::new();
    let mut push = |p: Point2| {
        if let (Some(s), Some(t)) = (circle_parameter(a, p, tol), circle_parameter(b, p, tol)) {
            crossings.push(Crossing {
                on_a: s,
                on_b: t,
                point: p,
                gap: 0.0,
            });
        }
    };
    if squared <= tol.confusion() * tol.confusion() {
        push(foot);
    } else {
        let offset = og_math::Vector2::new(-direction.y, direction.x) * squared.max(0.0).sqrt();
        push(foot + offset);
        push(foot - offset);
    }
    sort_crossings(&mut crossings);
    CurveIntersection {
        crossings,
        overlaps: Vec::new(),
    }
}

/// The parameter at which a circle passes through a point on it.
fn circle_parameter(curve: &og_geom::Circle2d, p: Point2, tol: Tolerances) -> Option<f64> {
    let c = curve.circle();
    let local = p - c.centre();
    let x = local.dot(c.frame().x().vector());
    let y = local.dot(c.frame().y().vector());
    let mut angle = y.atan2(x);
    if curve.is_reversed() {
        angle = -angle;
    }
    let angle = angle.rem_euclid(core::f64::consts::TAU);
    let (lo, hi) = curve.domain();
    // Fold into the arc's own range where the arc covers it.
    if angle >= lo - tol.parametric() && angle <= hi + tol.parametric() {
        return Some(angle.clamp(lo, hi));
    }
    let shifted = angle - core::f64::consts::TAU;
    if shifted >= lo - tol.parametric() && shifted <= hi + tol.parametric() {
        return Some(shifted.clamp(lo, hi));
    }
    None
}

// --- analytic, spatial -------------------------------------------------------

fn line_line_3d(
    a: &og_geom::LineCurve,
    b: &og_geom::LineCurve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> CurveIntersection<Point> {
    let (oa, da) = (a.axis().location, a.axis().direction.vector());
    let (ob, db) = (b.axis().location, b.axis().direction.vector());
    let cross = da.cross(db);
    let denominator = cross.square_magnitude();

    if denominator <= tol.angular() * tol.angular() {
        // Parallel: collinear overlap or nothing.
        let between = ob - oa;
        if between.cross(da).magnitude() > tol.confusion() {
            return CurveIntersection::empty();
        }
        let (a_lo, a_hi) = a.domain();
        let (b_lo, b_hi) = b.domain();
        let project = |p: Point| (p - oa).dot(da);
        let (s0, s1) = (project(ob + db * b_lo), project(ob + db * b_hi));
        let (lo, hi) = (s0.min(s1).max(a_lo), s0.max(s1).min(a_hi));
        if lo >= hi {
            return CurveIntersection::empty();
        }
        let back = |t: f64| (oa + da * t - ob).dot(db);
        return CurveIntersection {
            crossings: Vec::new(),
            overlaps: vec![Overlap {
                on_a: (lo, hi),
                on_b: order(back(lo), back(hi)),
            }],
        };
    }

    // Closest approach of two skew lines, in closed form.
    let between = ob - oa;
    let t = between.cross(db).dot(cross) / denominator;
    let s = between.cross(da).dot(cross) / denominator;
    let pa = oa + da * t;
    let pb = ob + db * s;
    let gap = pa.distance(pb);
    let (a_lo, a_hi) = a.domain();
    let (b_lo, b_hi) = b.domain();
    if gap > options.gap
        || t < a_lo - tol.parametric()
        || t > a_hi + tol.parametric()
        || s < b_lo - tol.parametric()
        || s > b_hi + tol.parametric()
    {
        return CurveIntersection::empty();
    }
    CurveIntersection {
        crossings: vec![Crossing {
            on_a: t,
            on_b: s,
            point: pa,
            gap,
        }],
        overlaps: Vec::new(),
    }
}

// --- the general path --------------------------------------------------------

/// Sampled segments of one curve, with the parameters they span.
struct Sampled<P> {
    points: Vec<P>,
    parameters: Vec<f64>,
}

fn sample_2d(curve: &PlanarCurve, n: usize, tol: Tolerances) -> Sampled<Point2> {
    let (lo, hi) = curve.domain();
    let mut points = Vec::with_capacity(n + 1);
    let mut parameters = Vec::with_capacity(n + 1);
    for i in 0..=n {
        #[allow(clippy::cast_precision_loss)]
        let t = lo + (hi - lo) * i as f64 / n as f64;
        if let Ok(p) = curve.point_at(t, tol) {
            points.push(p);
            parameters.push(t);
        }
    }
    Sampled { points, parameters }
}

fn sample_3d(curve: &Curve, n: usize, tol: Tolerances) -> Sampled<Point> {
    let (lo, hi) = curve.domain();
    let mut points = Vec::with_capacity(n + 1);
    let mut parameters = Vec::with_capacity(n + 1);
    for i in 0..=n {
        #[allow(clippy::cast_precision_loss)]
        let t = lo + (hi - lo) * i as f64 / n as f64;
        if let Ok(p) = curve.point_at(t, tol) {
            points.push(p);
            parameters.push(t);
        }
    }
    Sampled { points, parameters }
}

fn general_2d(
    a: &PlanarCurve,
    b: &PlanarCurve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> OgResult<CurveIntersection<Point2>> {
    let sa = sample_2d(a, options.samples, tol);
    let sb = sample_2d(b, options.samples, tol);

    let mut crossings: Vec<Crossing<Point2>> = Vec::new();
    for i in 1..sa.points.len() {
        for j in 1..sb.points.len() {
            let Some((ta, tb)) = segments_cross_2d(
                (sa.points[i - 1], sa.points[i]),
                (sb.points[j - 1], sb.points[j]),
            ) else {
                continue;
            };
            let seed_a = sa.parameters[i - 1] + (sa.parameters[i] - sa.parameters[i - 1]) * ta;
            let seed_b = sb.parameters[j - 1] + (sb.parameters[j] - sb.parameters[j - 1]) * tb;
            if let Some(found) = polish_2d(a, b, seed_a, seed_b, options, tol) {
                push_unique_2d(&mut crossings, found, tol);
            }
        }
    }
    sort_crossings(&mut crossings);
    Ok(CurveIntersection {
        crossings,
        overlaps: Vec::new(),
    })
}

fn general_3d(
    a: &Curve,
    b: &Curve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> OgResult<CurveIntersection<Point>> {
    let sa = sample_3d(a, options.samples, tol);
    let sb = sample_3d(b, options.samples, tol);

    // Segment pairs whose closest approach is within reach seed the polish.
    // The threshold is the sampling sag plus the acceptable gap: what could
    // converge is seeded, what could not is skipped.
    let mut reach = options.gap;
    for s in [&sa, &sb] {
        let longest = s
            .points
            .windows(2)
            .map(|w| w[0].distance(w[1]))
            .fold(0.0_f64, f64::max);
        reach += longest;
    }

    let mut crossings: Vec<Crossing<Point>> = Vec::new();
    for i in 1..sa.points.len() {
        for j in 1..sb.points.len() {
            let (ta, tb, gap) = segments_approach_3d(
                (sa.points[i - 1], sa.points[i]),
                (sb.points[j - 1], sb.points[j]),
            );
            if gap > reach {
                continue;
            }
            let seed_a = sa.parameters[i - 1] + (sa.parameters[i] - sa.parameters[i - 1]) * ta;
            let seed_b = sb.parameters[j - 1] + (sb.parameters[j] - sb.parameters[j - 1]) * tb;
            if let Some(found) = polish_3d(a, b, seed_a, seed_b, options, tol) {
                push_unique_3d(&mut crossings, found, tol);
            }
        }
    }
    sort_crossings(&mut crossings);
    Ok(CurveIntersection {
        crossings,
        overlaps: Vec::new(),
    })
}

/// Newton on `c1(t) - c2(s) = 0` in the plane.
fn polish_2d(
    a: &PlanarCurve,
    b: &PlanarCurve,
    seed_a: f64,
    seed_b: f64,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Option<Crossing<Point2>> {
    let system = |x: &[f64]| {
        let (t, s) = (clamp_2d(a, x[0]), clamp_2d(b, x[1]));
        let pa = a.point_at(t, tol).unwrap_or(Point2::ORIGIN);
        let pb = b.point_at(s, tol).unwrap_or(Point2::ORIGIN);
        let da = a.d1_at(t, tol).unwrap_or(og_math::Vector2::new(0.0, 0.0));
        let db = b.d1_at(s, tol).unwrap_or(og_math::Vector2::new(0.0, 0.0));
        (
            vec![pa.x - pb.x, pa.y - pb.y],
            vec![vec![da.x, -db.x], vec![da.y, -db.y]],
        )
    };
    let criteria = solve::Criteria {
        residual: tol.confusion() * 0.01,
        step: tol.parametric(),
        max_iterations: 40,
    };
    let found = solve::newton_system(system, &[seed_a, seed_b], criteria).ok()?;
    let (t, s) = (clamp_2d(a, found.value[0]), clamp_2d(b, found.value[1]));
    let pa = a.point_at(t, tol).ok()?;
    let pb = b.point_at(s, tol).ok()?;
    let gap = pa.distance(pb);
    if gap > options.gap {
        return None;
    }
    Some(Crossing {
        on_a: t,
        on_b: s,
        point: pa,
        gap,
    })
}

/// Gauss–Newton on the closest approach of two space curves.
///
/// Three equations would be overdetermined for two unknowns, so the system is
/// the two *stationarity* conditions — the gap vector perpendicular to both
/// tangents — whose solutions are the local closest approaches. The gap test
/// afterwards decides whether the approach found is a crossing.
fn polish_3d(
    a: &Curve,
    b: &Curve,
    seed_a: f64,
    seed_b: f64,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Option<Crossing<Point>> {
    let system = |x: &[f64]| {
        let (t, s) = (clamp_3d(a, x[0]), clamp_3d(b, x[1]));
        let pa = a.point_at(t, tol).unwrap_or(Point::ORIGIN);
        let pb = b.point_at(s, tol).unwrap_or(Point::ORIGIN);
        let da = a.derivatives_at(t, 2, tol).unwrap_or_default();
        let db = b.derivatives_at(s, 2, tol).unwrap_or_default();
        let zero = og_math::Vector::ZERO;
        let (d1a, d2a) = (
            da.get(1).copied().unwrap_or(zero),
            da.get(2).copied().unwrap_or(zero),
        );
        let (d1b, d2b) = (
            db.get(1).copied().unwrap_or(zero),
            db.get(2).copied().unwrap_or(zero),
        );
        let gap = pa - pb;
        (
            vec![gap.dot(d1a), -gap.dot(d1b)],
            vec![
                vec![d1a.dot(d1a) + gap.dot(d2a), -d1a.dot(d1b)],
                vec![-d1a.dot(d1b), d1b.dot(d1b) - gap.dot(d2b)],
            ],
        )
    };
    let criteria = solve::Criteria {
        residual: tol.confusion() * 0.01,
        step: tol.parametric(),
        max_iterations: 40,
    };
    let found = solve::newton_system(system, &[seed_a, seed_b], criteria).ok()?;
    let (t, s) = (clamp_3d(a, found.value[0]), clamp_3d(b, found.value[1]));
    let pa = a.point_at(t, tol).ok()?;
    let pb = b.point_at(s, tol).ok()?;
    let gap = pa.distance(pb);
    if gap > options.gap {
        return None;
    }
    Some(Crossing {
        on_a: t,
        on_b: s,
        point: pa,
        gap,
    })
}

// --- small helpers -----------------------------------------------------------

fn clamp_2d(curve: &PlanarCurve, t: f64) -> f64 {
    let (lo, hi) = curve.domain();
    if curve.is_periodic() {
        let span = hi - lo;
        if span > 0.0 {
            return lo + (t - lo).rem_euclid(span);
        }
    }
    t.clamp(lo, hi)
}

fn clamp_3d(curve: &Curve, t: f64) -> f64 {
    let (lo, hi) = curve.domain();
    if curve.is_periodic() {
        let span = hi - lo;
        if span > 0.0 {
            return lo + (t - lo).rem_euclid(span);
        }
    }
    t.clamp(lo, hi)
}

/// Where two planar segments cross, as fractions along each.
fn segments_cross_2d(a: (Point2, Point2), b: (Point2, Point2)) -> Option<(f64, f64)> {
    let da = a.1 - a.0;
    let db = b.1 - b.0;
    let cross = da.cross(db);
    if cross.abs() <= f64::MIN_POSITIVE {
        return None;
    }
    let between = b.0 - a.0;
    let t = between.cross(db) / cross;
    let s = between.cross(da) / cross;
    if !(0.0..=1.0).contains(&t) || !(0.0..=1.0).contains(&s) {
        return None;
    }
    Some((t, s))
}

/// The closest approach of two spatial segments, as fractions and a distance.
fn segments_approach_3d(a: (Point, Point), b: (Point, Point)) -> (f64, f64, f64) {
    let da = a.1 - a.0;
    let db = b.1 - b.0;
    let between = a.0 - b.0;
    let (aa, bb, ab) = (da.dot(da), db.dot(db), da.dot(db));
    let (ad, bd) = (da.dot(between), db.dot(between));
    let denominator = ab.mul_add(-ab, aa * bb);

    let (mut t, mut s) = if denominator.abs() <= f64::MIN_POSITIVE {
        (
            0.0,
            if bb > 0.0 {
                (bd / bb).clamp(0.0, 1.0)
            } else {
                0.0
            },
        )
    } else {
        (
            (ab.mul_add(bd, -(bb * ad)) / denominator).clamp(0.0, 1.0),
            (aa.mul_add(bd, -(ab * ad)) / denominator).clamp(0.0, 1.0),
        )
    };
    // One clamped end may pull the other; a single re-projection settles it.
    if bb > 0.0 {
        s = ((da.dot(between) + t * aa - 0.0).mul_add(0.0, db.dot(between + da * t)) / bb)
            .clamp(0.0, 1.0);
    }
    if aa > 0.0 {
        t = (da.dot(db * s - between) / aa).clamp(0.0, 1.0);
    }
    let pa = a.0 + da * t;
    let pb = b.0 + db * s;
    (t, s, pa.distance(pb))
}

fn order(a: f64, b: f64) -> (f64, f64) {
    if a <= b { (a, b) } else { (b, a) }
}

fn sort_crossings<P>(crossings: &mut [Crossing<P>]) {
    crossings.sort_by(|x, y| {
        x.on_a
            .partial_cmp(&y.on_a)
            .unwrap_or(core::cmp::Ordering::Equal)
    });
}

fn push_unique_2d(crossings: &mut Vec<Crossing<Point2>>, found: Crossing<Point2>, tol: Tolerances) {
    let reach = tol.confusion() * 100.0;
    if crossings
        .iter()
        .any(|c| c.point.distance(found.point) <= reach)
    {
        return;
    }
    crossings.push(found);
}

fn push_unique_3d(crossings: &mut Vec<Crossing<Point>>, found: Crossing<Point>, tol: Tolerances) {
    let reach = tol.confusion() * 100.0;
    if crossings
        .iter()
        .any(|c| c.point.distance(found.point) <= reach)
    {
        return;
    }
    crossings.push(found);
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;
    use og_geom::{BSpline2d, Circle2d, CircleCurve, Line2d, LineCurve};
    use og_math::{Circle, Circle2, Direction2, Frame, Frame2, KnotVector, Vector2};

    const T: Tolerances = Tolerances::millimetres();

    fn line2(from: Point2, to: Point2) -> PlanarCurve {
        Line2d::segment(from, to, T).unwrap().into()
    }

    fn circle2(centre: Point2, radius: f64) -> PlanarCurve {
        Circle2d::new(
            Circle2::new(
                Frame2::new(centre, Direction2::new(Vector2::new(1.0, 0.0), T).unwrap()),
                radius,
                T,
            )
            .unwrap(),
        )
        .into()
    }

    #[test]
    fn two_lines_cross_where_algebra_says() {
        let a = line2(Point2::new(0.0, 0.0), Point2::new(4.0, 4.0));
        let b = line2(Point2::new(0.0, 4.0), Point2::new(4.0, 0.0));
        let found = intersect_curves_2d(&a, &b, CurveCurveOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 1);
        let hit = &found.crossings[0];
        assert!(hit.point.is_equal(Point2::new(2.0, 2.0), T));
        // Parameters are arc length on a segment.
        approx::assert_relative_eq!(hit.on_a, 8.0_f64.sqrt(), epsilon = 1e-9);

        // Segments that would cross beyond their ends do not.
        let short = line2(Point2::new(0.0, 4.0), Point2::new(1.0, 3.0));
        assert!(
            intersect_curves_2d(&a, &short, CurveCurveOptions::default(), T)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn collinear_lines_overlap_rather_than_crossing_everywhere() {
        let a = line2(Point2::new(0.0, 0.0), Point2::new(10.0, 0.0));
        let b = line2(Point2::new(4.0, 0.0), Point2::new(20.0, 0.0));
        let found = intersect_curves_2d(&a, &b, CurveCurveOptions::default(), T).unwrap();
        assert!(found.crossings.is_empty());
        assert_eq!(found.overlaps.len(), 1);
        let overlap = &found.overlaps[0];
        approx::assert_relative_eq!(overlap.on_a.0, 4.0, epsilon = 1e-9);
        approx::assert_relative_eq!(overlap.on_a.1, 10.0, epsilon = 1e-9);
        approx::assert_relative_eq!(overlap.on_b.0, 0.0, epsilon = 1e-9);
        approx::assert_relative_eq!(overlap.on_b.1, 6.0, epsilon = 1e-9);

        // Parallel but apart: nothing.
        let above = line2(Point2::new(0.0, 1.0), Point2::new(10.0, 1.0));
        assert!(
            intersect_curves_2d(&a, &above, CurveCurveOptions::default(), T)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn a_line_meets_a_circle_in_two_points_one_or_none() {
        let circle = circle2(Point2::new(0.0, 0.0), 2.0);
        let through = line2(Point2::new(-5.0, 0.0), Point2::new(5.0, 0.0));
        let found =
            intersect_curves_2d(&through, &circle, CurveCurveOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 2);
        for hit in &found.crossings {
            approx::assert_relative_eq!(
                hit.point.distance(Point2::new(0.0, 0.0)),
                2.0,
                epsilon = 1e-9
            );
            // The circle parameter really evaluates to the crossing point.
            let PlanarCurve::Circle(_) = &circle else {
                unreachable!()
            };
            let on_circle = circle.point_at(hit.on_b, T).unwrap();
            assert!(on_circle.is_equal(hit.point, T));
        }

        let tangent = line2(Point2::new(-5.0, 2.0), Point2::new(5.0, 2.0));
        assert_eq!(
            intersect_curves_2d(&tangent, &circle, CurveCurveOptions::default(), T)
                .unwrap()
                .crossings
                .len(),
            1
        );
        let missing = line2(Point2::new(-5.0, 3.0), Point2::new(5.0, 3.0));
        assert!(
            intersect_curves_2d(&missing, &circle, CurveCurveOptions::default(), T)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn two_circles_cross_touch_coincide_or_miss() {
        let a = circle2(Point2::new(0.0, 0.0), 2.0);

        let crossing = circle2(Point2::new(3.0, 0.0), 2.0);
        let found = intersect_curves_2d(&a, &crossing, CurveCurveOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 2);
        for hit in &found.crossings {
            let on_a = a.point_at(hit.on_a, T).unwrap();
            let on_b = crossing.point_at(hit.on_b, T).unwrap();
            assert!(on_a.is_equal(hit.point, T));
            assert!(on_b.is_equal(hit.point, T));
        }

        let touching = circle2(Point2::new(4.0, 0.0), 2.0);
        assert_eq!(
            intersect_curves_2d(&a, &touching, CurveCurveOptions::default(), T)
                .unwrap()
                .crossings
                .len(),
            1
        );

        let same = circle2(Point2::new(0.0, 0.0), 2.0);
        let coincident = intersect_curves_2d(&a, &same, CurveCurveOptions::default(), T).unwrap();
        assert!(coincident.crossings.is_empty());
        assert_eq!(coincident.overlaps.len(), 1);

        let apart = circle2(Point2::new(10.0, 0.0), 2.0);
        assert!(
            intersect_curves_2d(&a, &apart, CurveCurveOptions::default(), T)
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn the_general_path_handles_what_has_no_closed_form() {
        // A spline sine-ish wave against a line: three crossings, found by
        // sampling and polished by Newton to rounding.
        let wave: PlanarCurve = BSpline2d::new(
            KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0], 3).unwrap(),
            vec![
                Point2::new(0.0, -1.0),
                Point2::new(1.0, 3.0),
                Point2::new(2.0, -3.0),
                Point2::new(3.0, 3.0),
                Point2::new(4.0, -1.0),
            ],
            T,
        )
        .unwrap()
        .into();
        let axis = line2(Point2::new(-1.0, 0.0), Point2::new(5.0, 0.0));
        let found = intersect_curves_2d(&wave, &axis, CurveCurveOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 3, "a wave crosses its axis thrice");
        for hit in &found.crossings {
            assert!(hit.gap < 1e-9);
            assert!(hit.point.y.abs() < 1e-9);
            let on_wave = wave.point_at(hit.on_a, T).unwrap();
            assert!(on_wave.is_equal(hit.point, T));
        }
    }

    #[test]
    fn space_curves_cross_within_a_gap_and_report_it() {
        // Two circles that would cross in a shared plane, with one lifted a
        // hair out of it: the crossings become passes with a real, small gap
        // that must be reported, not zeroed. (Not chain links — a first draft
        // of this test used linked circles, and linked circles never approach:
        // passing through each other's *disks* is what linked means, and these
        // radii hold the curves a constant two units apart.)
        let a: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.0, T).unwrap()).into();
        let lifted = Frame::new(
            Point::new(3.0, 0.0, 0.001),
            og_math::Direction::Z,
            og_math::Direction::X,
            T,
        )
        .unwrap();
        let b: Curve = CircleCurve::new(Circle::new(lifted, 2.0, T).unwrap()).into();

        let options = CurveCurveOptions {
            gap: 1e-2,
            ..CurveCurveOptions::default()
        };
        let found = intersect_curves(&a, &b, options, T).unwrap();
        assert_eq!(found.crossings.len(), 2, "two near-crossings");
        for hit in &found.crossings {
            assert!(hit.gap > 1e-4, "the gap is real and must not be zeroed");
            assert!(hit.gap < 2e-3, "but small: {}", hit.gap);
        }

        // Tighten the gap below the offset and the crossings vanish.
        let strict = CurveCurveOptions {
            gap: 1e-5,
            ..CurveCurveOptions::default()
        };
        assert!(intersect_curves(&a, &b, strict, T).unwrap().is_empty());
    }

    #[test]
    fn skew_lines_in_space_miss_and_close_ones_meet() {
        let a: Curve = LineCurve::segment(Point::ORIGIN, Point::new(10.0, 0.0, 0.0), T)
            .unwrap()
            .into();
        let skew: Curve =
            LineCurve::segment(Point::new(0.0, -5.0, 1.0), Point::new(0.0, 5.0, 1.0), T)
                .unwrap()
                .into();
        assert!(
            intersect_curves(&a, &skew, CurveCurveOptions::default(), T)
                .unwrap()
                .is_empty(),
            "a unit apart is not a crossing"
        );

        let meeting: Curve =
            LineCurve::segment(Point::new(5.0, -5.0, 0.0), Point::new(5.0, 5.0, 0.0), T)
                .unwrap()
                .into();
        let found = intersect_curves(&a, &meeting, CurveCurveOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 1);
        assert!(
            found.crossings[0]
                .point
                .is_equal(Point::new(5.0, 0.0, 0.0), T)
        );
        assert!(found.crossings[0].gap < 1e-12);

        // Collinear 3D lines overlap.
        let collinear: Curve =
            LineCurve::segment(Point::new(4.0, 0.0, 0.0), Point::new(20.0, 0.0, 0.0), T)
                .unwrap()
                .into();
        let shared = intersect_curves(&a, &collinear, CurveCurveOptions::default(), T).unwrap();
        assert_eq!(shared.overlaps.len(), 1);
    }

    #[test]
    fn unusable_options_are_refused() {
        let a = line2(Point2::new(0.0, 0.0), Point2::new(1.0, 0.0));
        for options in [
            CurveCurveOptions {
                samples: 1,
                ..CurveCurveOptions::default()
            },
            CurveCurveOptions {
                gap: 0.0,
                ..CurveCurveOptions::default()
            },
            CurveCurveOptions {
                gap: f64::NAN,
                ..CurveCurveOptions::default()
            },
        ] {
            assert!(intersect_curves_2d(&a, &a.clone(), options, T).is_err());
        }
    }
}
