//! Where two curves meet, in the plane and in space.
//!
//! *Elsewhere* these are `Geom2dAPI_InterCurveCurve` and `IntCurve` for the
//! plane, and extrema-based crossing for space. The planar case is the
//! load-bearing one: boolean face splitting happens in a surface's parameter
//! space, and the curves it splits with are pcurves, so 2D curve/curve is the
//! operation the whole §8 pipeline stands on.
//!
//! # Two curves in space generically miss
//!
//! In the plane, two curves that cross, cross. In space they pass by: a
//! crossing is two points closer than a tolerance, not an exact common point,
//! and pretending otherwise would make every 3D result empty. So the 3D
//! answer reports the *gap* it achieved at each crossing, and the caller's
//! tolerance decides what counts. The 2D answer reports gaps too (a solved
//! crossing is still a pair of floats), but there the gap is rounding, not
//! geometry.
//!
//! # Overlap is an answer, not a failure
//!
//! Two collinear lines, two arcs of one circle: where the supports coincide,
//! "the intersection points" do not exist; the intersection is a stretch of
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

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve2d, Curve3d, PlanarCurve};
use ogeom_math::{Point, Point2, solve};

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
    /// How far along the curves this contact could honestly sit: zero for a
    /// transversal crossing, the length of the touching run where the curves
    /// meet tangentially: there the closest approach is anywhere in a
    /// valley the width of the gap, and a consumer placing a vertex at it
    /// owns that much doubt.
    pub reach: f64,
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
    /// The analytic same-support cases (collinear lines, arcs of one
    /// circle) come back exactly. In space, the sampling path also reports
    /// a stretch along which the first curve's samples stay within the gap
    /// of the second, its ends bisected to parametric resolution and its
    /// correspondence stated by those ends alone: a fitted section tracing
    /// the arc it was cut along is one overlap, not a row of crossings. A
    /// stretch shorter than two samples of the first curve is still read as
    /// whatever crossings the sampling finds; in the plane, only the
    /// analytic cases are detected.
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
/// Analytic pairs (lines and circles) are answered in closed form, overlaps
/// included. Everything else goes through sampling and Newton.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the options
/// are unusable.
pub fn intersect_curves_2d(
    a: &PlanarCurve,
    b: &PlanarCurve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> OgeomResult<CurveIntersection<Point2>> {
    check(options)?;
    // Through any trim to the basis the closed forms answer for, as in
    // space; the answer is clipped to the trims' windows after.
    let (basis_a, window_a) = through_trim_2d(a);
    let (basis_b, window_b) = through_trim_2d(b);
    let found = match (basis_a, basis_b) {
        (PlanarCurve::Line(x), PlanarCurve::Line(y)) => line_line_2d(x, y, tol),
        (PlanarCurve::Line(x), PlanarCurve::Circle(y)) => line_circle_2d(x, y, false, tol),
        (PlanarCurve::Circle(x), PlanarCurve::Line(y)) => line_circle_2d(y, x, true, tol),
        (PlanarCurve::Circle(x), PlanarCurve::Circle(y)) => circle_circle_2d(x, y, tol),
        _ => return general_2d(a, b, options, tol),
    };
    let side = |curve: &PlanarCurve, window: Option<(f64, f64)>| Side {
        period: {
            let (lo, hi) = curve.domain();
            (curve.is_periodic() && hi > lo).then_some(hi - lo)
        },
        domain: curve.domain(),
        window,
    };
    Ok(clipped_to_windows(
        found,
        side(basis_a, window_a),
        side(basis_b, window_b),
        tol,
    ))
}

/// A planar curve seen through a trim, as [`through_trim`] sees one in
/// space.
fn through_trim_2d(curve: &PlanarCurve) -> (&PlanarCurve, Option<(f64, f64)>) {
    match curve {
        PlanarCurve::Trimmed(t) if !t.is_reversed() => (t.basis(), Some(t.domain())),
        other => (other, None),
    }
}

/// One side of a curve pair as clipping sees it: its basis's period and
/// domain, and the window a trim restricts it to.
#[derive(Debug, Clone, Copy)]
struct Side {
    period: Option<f64>,
    domain: (f64, f64),
    window: Option<(f64, f64)>,
}

impl Side {
    fn of(curve: &Curve, window: Option<(f64, f64)>) -> Self {
        let (lo, hi) = curve.domain();
        Self {
            period: (curve.is_periodic() && hi > lo).then_some(hi - lo),
            domain: (lo, hi),
            window,
        }
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
) -> OgeomResult<CurveIntersection<Point>> {
    check(options)?;
    let (basis_a, window_a) = through_trim(a);
    let (basis_b, window_b) = through_trim(b);
    if let Some(found) = same_curve_3d(basis_a, basis_b) {
        return Ok(clipped_to_windows(
            found,
            Side::of(basis_a, window_a),
            Side::of(basis_b, window_b),
            tol,
        ));
    }
    if let Some(found) = analytic_3d(basis_a, basis_b, options, tol) {
        return Ok(clipped_to_windows(
            found,
            Side::of(basis_a, window_a),
            Side::of(basis_b, window_b),
            tol,
        ));
    }
    general_3d(a, b, options, tol)
}

/// A curve seen through a trim: the curve the analytic path can answer for,
/// and the window it is restricted to.
///
/// The trim shares its basis's parameterization, so the window is stated in
/// the same numbers the analytic answer comes back in and clipping is an
/// interval intersection rather than a change of variable. A *reversed* trim
/// does renumber, so it is left to the sampling path rather than mis-read.
fn through_trim(curve: &Curve) -> (&Curve, Option<(f64, f64)>) {
    match curve {
        Curve::Trimmed(t) if !t.is_reversed() => (t.basis(), Some(Curve3d::domain(&**t))),
        other => (other, None),
    }
}

/// The closed-form answers for a pair of space curves, or `None` where there
/// is none and the sampling path is the honest route.
fn analytic_3d(
    a: &Curve,
    b: &Curve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Option<CurveIntersection<Point>> {
    match (a, b) {
        (Curve::Line(x), Curve::Line(y)) => Some(line_line_3d(x, y, options, tol)),
        (Curve::Circle(x), Curve::Circle(y)) => same_circle_3d(x, y, tol)
            .or_else(|| coplanar_circles_3d(x, y, options, tol))
            .or_else(|| skew_conics_3d(a, b, options, tol)),
        (Curve::Ellipse(x), Curve::Ellipse(y)) => {
            same_ellipse_3d(x, y, tol).or_else(|| skew_conics_3d(a, b, options, tol))
        }
        (Curve::Circle(_), Curve::Ellipse(_)) | (Curve::Ellipse(_), Curve::Circle(_)) => {
            skew_conics_3d(a, b, options, tol)
        }
        (Curve::Line(x), Curve::Circle(_) | Curve::Ellipse(_)) => line_conic_3d(x, b, options, tol)
            .map(|mut found| {
                for c in &mut found.crossings {
                    core::mem::swap(&mut c.on_a, &mut c.on_b);
                }
                found.crossings.sort_by(|x, y| x.on_a.total_cmp(&y.on_a));
                found
            }),
        (Curve::Circle(_) | Curve::Ellipse(_), Curve::Line(y)) => line_conic_3d(y, a, options, tol),
        _ => None,
    }
}

/// A circle or ellipse, as its centre, the two semi-axis vectors its
/// parameter turns between, and its plane's normal; `None` for any other
/// curve, or one whose parameter runs backwards.
fn conic_of(
    curve: &Curve,
) -> Option<(
    Point,
    ogeom_math::Vector,
    ogeom_math::Vector,
    ogeom_math::Vector,
)> {
    match curve {
        Curve::Circle(c) if !c.is_reversed() => {
            let circle = c.circle();
            let f = circle.frame();
            Some((
                f.origin(),
                f.x().vector() * circle.radius(),
                f.y().vector() * circle.radius(),
                f.z().vector(),
            ))
        }
        Curve::Ellipse(e) if !e.is_reversed() => {
            let ellipse = e.ellipse();
            let f = ellipse.frame();
            Some((
                f.origin(),
                f.x().vector() * ellipse.major_radius(),
                f.y().vector() * ellipse.minor_radius(),
                f.z().vector(),
            ))
        }
        _ => None,
    }
}

/// Two circles or ellipses in planes that are not one plane, in closed
/// form: where the first meets the second's plane (`alpha cos t + beta sin
/// t + gamma = 0`, at most twice), kept where the second curve passes
/// within the gap. Parallel planes apart share no point. `None` for a
/// coplanar pair, which the sampling path answers.
fn skew_conics_3d(
    a: &Curve,
    b: &Curve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Option<CurveIntersection<Point>> {
    let (ca, ua, va, na) = conic_of(a)?;
    let (cb, _, _, nb) = conic_of(b)?;
    if na.cross(nb).magnitude() <= tol.angular() {
        return ((ca - cb).dot(nb).abs() > options.gap.max(tol.confusion()))
            .then(CurveIntersection::empty);
    }
    let (alpha, beta, gamma) = (nb.dot(ua), nb.dot(va), nb.dot(ca - cb));
    let size = alpha.hypot(beta);
    // The first conic stands within the gap of the other's plane over a
    // stretch as wide as the gap is against how far it rises from that
    // plane. Planes all but parallel (a rim fitted to one facet group, a
    // section through another) make that stretch wide, and the curves can
    // pass within the gap anywhere along it, not only where one crosses the
    // other's plane; the sampling path measures such a pass.
    if options.gap.max(tol.confusion()) > size * 1e-3 {
        return None;
    }
    let mut crossings: Vec<Crossing<Point>> = Vec::new();
    if size > 0.0 && gamma.abs() <= size * (1.0 + 1e-12) {
        let phase = beta.atan2(alpha);
        let turn = (-gamma / size).clamp(-1.0, 1.0).acos();
        let (a_lo, a_hi) = a.domain();
        let (b_lo, b_hi) = b.domain();
        let tau = core::f64::consts::TAU;
        let into = |t: f64, lo: f64| lo + (t - lo).rem_euclid(tau);
        let roots = if turn <= 1e-12 {
            vec![phase]
        } else {
            vec![phase - turn, phase + turn]
        };
        for root in roots {
            let t = into(root, a_lo);
            if t > a_hi + tol.parametric() {
                continue;
            }
            let point = a.point_at(t, tol).ok()?;
            // The closed form reads the curve as its centre and semi-axes;
            // a point that does not land on the other plane means it read
            // it wrong, and the sampling path answers instead.
            if (point - cb).dot(nb).abs() > tol.confusion() * 10.0 {
                return None;
            }
            let s = match b {
                Curve::Circle(c) => {
                    ogeom_math::elementary::circle_parameter(&c.circle(), point, tol)
                }
                Curve::Ellipse(e) => {
                    ogeom_math::elementary::ellipse_parameter(&e.ellipse(), point, tol)
                }
                _ => return None,
            };
            let Ok(s) = s else {
                continue;
            };
            let s = into(s, b_lo);
            if s > b_hi + tol.parametric() {
                continue;
            }
            let gap = b.point_at(s, tol).ok()?.distance(point);
            if gap > options.gap {
                continue;
            }
            crossings.push(Crossing {
                on_a: t,
                on_b: s,
                point,
                gap,
                reach: 0.0,
            });
        }
    }
    crossings.sort_by(|x, y| x.on_a.total_cmp(&y.on_a));
    Some(CurveIntersection {
        crossings,
        overlaps: Vec::new(),
    })
}

/// A circle or ellipse against a line, in closed form, the conic first.
///
/// A line in the conic's plane meets it where a quadratic along the line
/// vanishes, at most twice; a line through the plane meets it at most where
/// it pierces the plane. Kept where the two pass within the gap. A line
/// running nearly along the plane without lying in it, or a crossing nearly
/// tangent, is left to the sampling path, which measures how far such a
/// touch reaches.
fn line_conic_3d(
    line: &ogeom_geom::LineCurve,
    conic: &Curve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Option<CurveIntersection<Point>> {
    let (centre, u, v, normal) = conic_of(conic)?;
    let (origin, along) = (line.axis().location, line.axis().direction.vector());
    let (a_len, b_len) = (u.magnitude(), v.magnitude());
    if a_len <= tol.confusion() || b_len <= tol.confusion() {
        return None;
    }
    let (ux, vy) = (u / a_len, v / b_len);
    let lean = along.dot(normal);
    let height = (origin - centre).dot(normal);
    let mut ts: Vec<f64> = Vec::new();
    if lean.abs() <= tol.angular() {
        if height.abs() > options.gap.max(tol.confusion()) {
            return Some(CurveIntersection::empty());
        }
        // (x0 + t dx)^2 / a^2 + (y0 + t dy)^2 / b^2 = 1, in the plane.
        let (x0, y0) = ((origin - centre).dot(ux), (origin - centre).dot(vy));
        let (dx, dy) = (along.dot(ux), along.dot(vy));
        let qa = dx * dx / (a_len * a_len) + dy * dy / (b_len * b_len);
        let qb = 2.0 * (x0 * dx / (a_len * a_len) + y0 * dy / (b_len * b_len));
        let qc = x0 * x0 / (a_len * a_len) + y0 * y0 / (b_len * b_len) - 1.0;
        let disc = qb.mul_add(qb, -4.0 * qa * qc);
        if qa <= 0.0 {
            return None;
        }
        if disc < 0.0 {
            let t = -qb / (2.0 * qa);
            let p = origin + along * t;
            let foot = conic_parameter(conic, p, tol)?;
            let gap = conic.point_at(foot, tol).ok()?.distance(p);
            return (gap > options.gap).then(CurveIntersection::empty);
        }
        let root = disc.sqrt();
        ts.push((-qb - root) / (2.0 * qa));
        ts.push((-qb + root) / (2.0 * qa));
    } else if lean.abs() >= 0.1 {
        ts.push(-height / lean);
    } else {
        return None;
    }
    let (lo, hi) = line.domain();
    let (c_lo, c_hi) = conic.domain();
    let tau = core::f64::consts::TAU;
    let mut crossings: Vec<Crossing<Point>> = Vec::new();
    for t in ts {
        if t < lo - tol.parametric() || t > hi + tol.parametric() {
            continue;
        }
        let point = origin + along * t;
        let s = conic_parameter(conic, point, tol)?;
        let s = c_lo + (s - c_lo).rem_euclid(tau);
        if s > c_hi + tol.parametric() {
            continue;
        }
        let on_conic = conic.point_at(s, tol).ok()?;
        let gap = on_conic.distance(point);
        if gap > options.gap {
            continue;
        }
        let tangent = conic.d1_at(s, tol).ok()?;
        if tangent.cross(along).magnitude() <= 1e-3 * tangent.magnitude() {
            return None;
        }
        crossings.push(Crossing {
            on_a: s,
            on_b: t,
            point: on_conic,
            gap,
            reach: 0.0,
        });
    }
    crossings.sort_by(|x, y| x.on_a.total_cmp(&y.on_a));
    crossings.dedup_by(|x, y| (x.on_a - y.on_a).abs() <= tol.parametric());
    Some(CurveIntersection {
        crossings,
        overlaps: Vec::new(),
    })
}

/// Where a point lies along a circle or ellipse, by projection.
fn conic_parameter(conic: &Curve, point: Point, tol: Tolerances) -> Option<f64> {
    match conic {
        Curve::Circle(c) => ogeom_math::elementary::circle_parameter(&c.circle(), point, tol).ok(),
        Curve::Ellipse(e) => {
            ogeom_math::elementary::ellipse_parameter(&e.ellipse(), point, tol).ok()
        }
        _ => None,
    }
}

/// Restrict an answer about two whole curves to the windows their trims
/// actually cover.
///
/// Both parts matter. A crossing is kept only where *both* parameters fall
/// inside their window, on a periodic basis after whichever whole turn
/// brings them there. An overlap is an interval on each side tied by an
/// affine correspondence, so it is clipped on one side, carried across, and
/// clipped again, and what comes back is the stretch both trims really share.
fn clipped_to_windows<P>(
    found: CurveIntersection<P>,
    a: Side,
    b: Side,
    tol: Tolerances,
) -> CurveIntersection<P> {
    let (window_a, window_b) = (a.window, b.window);
    if window_a.is_none() && window_b.is_none() {
        return found;
    }
    let (pa, pb) = (a.period, b.period);
    let slack = tol.parametric();
    let placed = |t: f64, window: Option<(f64, f64)>, period: Option<f64>| -> Option<f64> {
        let Some((lo, hi)) = window else {
            return Some(t);
        };
        for k in [0.0, 1.0, -1.0, 2.0, -2.0] {
            let shifted = period.map_or(t, |p| p.mul_add(k, t));
            if shifted >= lo - slack && shifted <= hi + slack {
                return Some(shifted);
            }
            if period.is_none() {
                break;
            }
        }
        None
    };

    let mut crossings = Vec::with_capacity(found.crossings.len());
    for crossing in found.crossings {
        let (Some(on_a), Some(on_b)) = (
            placed(crossing.on_a, window_a, pa),
            placed(crossing.on_b, window_b, pb),
        ) else {
            continue;
        };
        crossings.push(Crossing {
            on_a,
            on_b,
            ..crossing
        });
    }

    let mut overlaps: Vec<Overlap> = Vec::with_capacity(found.overlaps.len());
    let ordered = |r: (f64, f64)| if r.0 <= r.1 { r } else { (r.1, r.0) };
    let meet = |x: (f64, f64), y: (f64, f64)| -> Option<(f64, f64)> {
        let both = (x.0.max(y.0), x.1.min(y.1));
        (both.1 - both.0 > slack).then_some(both)
    };
    let (domain_a, domain_b) = (a.domain, b.domain);
    let shifts = |period: Option<f64>| -> Vec<f64> {
        period.map_or_else(
            || vec![0.0],
            |p| [0.0, 1.0, -1.0, 2.0, -2.0].iter().map(|k| k * p).collect(),
        )
    };
    for overlap in found.overlaps {
        let span_a = overlap.on_a.1 - overlap.on_a.0;
        let span_b = overlap.on_b.1 - overlap.on_b.0;
        if span_a.abs() <= f64::MIN_POSITIVE || span_b.abs() <= f64::MIN_POSITIVE {
            continue;
        }
        let rate = span_b / span_a;
        let to_b = |t: f64| overlap.on_b.0 + rate * (t - overlap.on_a.0);
        let to_a = |t: f64| overlap.on_a.0 + (t - overlap.on_b.0) / rate;
        // An overlap a whole turn long is the same set traced twice: the
        // correspondence holds round and round, and each side's window is
        // all that limits it. Otherwise it holds on its own stretch, at
        // whichever whole turn each window meets it.
        let whole = pa.is_some_and(|p| span_a.abs() >= p - slack)
            && pb.is_some_and(|p| span_b.abs() >= p - slack);
        let window_a = window_a.map_or(domain_a, ordered);
        let window_b = window_b.map_or(domain_b, ordered);
        let mut pieces_a: Vec<(f64, f64, f64)> = Vec::new();
        if whole {
            pieces_a.push((window_a.0, window_a.1, 0.0));
        } else {
            for shift in shifts(pa) {
                let stretch = ordered(overlap.on_a);
                if let Some(piece) = meet((stretch.0 + shift, stretch.1 + shift), window_a) {
                    pieces_a.push((piece.0, piece.1, shift));
                }
            }
        }
        for (lo, hi, shift_a) in pieces_a {
            // On `b`, at whichever whole turn of it lands in its window.
            let image = ordered((to_b(lo - shift_a), to_b(hi - shift_a)));
            for shift_b in shifts(pb) {
                let Some(on_b) = meet((image.0 + shift_b, image.1 + shift_b), window_b) else {
                    continue;
                };
                let back = |t: f64| to_a(t - shift_b) + shift_a;
                let on_a = ordered((back(on_b.0), back(on_b.1)));
                if on_a.1 - on_a.0 <= slack {
                    continue;
                }
                // Walked as `a` walks: `b`'s ends follow `a`'s.
                let forward = |t: f64| to_b(t - shift_a) + shift_b;
                let (on_a, on_b) = if span_a >= 0.0 {
                    (on_a, (forward(on_a.0), forward(on_a.1)))
                } else {
                    ((on_a.1, on_a.0), (forward(on_a.1), forward(on_a.0)))
                };
                let repeated = overlaps.iter().any(|o| {
                    (o.on_a.0 - on_a.0).abs() <= slack && (o.on_a.1 - on_a.1).abs() <= slack
                });
                if !repeated {
                    overlaps.push(Overlap { on_a, on_b });
                }
            }
        }
    }
    CurveIntersection {
        crossings,
        overlaps,
    }
}

/// Two circles tracing the same point set in space: the circle counterpart of
/// collinear lines, and the one 3D circle pair the sampling path cannot
/// answer: every sample is a hit, and "the crossings" do not exist. Distinct
/// circles return `None` and fall through to the general machinery, which
/// handles genuinely crossing pairs.
/// Two circles in one plane, in closed form: where they cross, from the
/// radical line in the plane. Exact however shallow the crossing: two rims a
/// micron apart meet at a fraction of a milliradian, and sampled, the touch
/// would read as a run millimetres long, which is how far the two stay
/// within tolerance of each other, not where they cross. `None` for circles
/// in different planes, sharing a centre, or within the weld distance of
/// touching, which the other paths answer.
fn coplanar_circles_3d(
    a: &ogeom_geom::CircleCurve,
    b: &ogeom_geom::CircleCurve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Option<CurveIntersection<Point>> {
    let (ca, cb) = (a.circle(), b.circle());
    let normal = ca.frame().z().vector();
    if normal.cross(cb.frame().z().vector()).magnitude() > tol.angular() {
        return None;
    }
    let between = cb.centre() - ca.centre();
    if between.dot(normal).abs() > tol.confusion() {
        return None;
    }
    let in_plane = between - normal * between.dot(normal);
    let distance = in_plane.magnitude();
    let (ra, rb) = (ca.radius(), cb.radius());
    if distance <= tol.confusion() {
        return None;
    }
    if distance > ra + rb + options.gap || distance < (ra - rb).abs() - options.gap {
        return Some(CurveIntersection::empty());
    }
    // Within the weld distance of touching, the two crossings are one touch
    // the root of a rounding error has pulled apart.
    let weld = tol.confusion() * 1e2;
    if (distance - (ra + rb)).abs() <= weld || (distance - (ra - rb).abs()).abs() <= weld {
        return None;
    }
    let along = distance.mul_add(distance, ra.mul_add(ra, -(rb * rb))) / (2.0 * distance);
    let squared = ra.mul_add(ra, -(along * along));
    if squared <= tol.confusion() * tol.confusion() {
        return None;
    }
    let half = squared.sqrt();
    let ux = in_plane / distance;
    let uy = normal.cross(ux);
    let parameter = |curve: &ogeom_geom::CircleCurve, p: Point| -> f64 {
        let local = curve.circle().frame().to_local(p);
        let angle = local.y.atan2(local.x);
        let angle = if curve.is_reversed() { -angle } else { angle };
        let (lo, _) = Curve3d::domain(curve);
        lo + (angle - lo).rem_euclid(core::f64::consts::TAU)
    };
    let mut crossings: Vec<Crossing<Point>> = [half, -half]
        .into_iter()
        .map(|h| {
            let point = ca.centre() + ux * along + uy * h;
            Crossing {
                on_a: parameter(a, point),
                on_b: parameter(b, point),
                point,
                gap: 0.0,
                reach: 0.0,
            }
        })
        .collect();
    sort_crossings(&mut crossings);
    Some(CurveIntersection {
        crossings,
        overlaps: Vec::new(),
    })
}

fn same_circle_3d(
    a: &ogeom_geom::CircleCurve,
    b: &ogeom_geom::CircleCurve,
    tol: Tolerances,
) -> Option<CurveIntersection<Point>> {
    let (ca, cb) = (a.circle(), b.circle());
    if ca.centre().distance(cb.centre()) > tol.confusion() {
        return None;
    }
    if (ca.radius() - cb.radius()).abs() > tol.confusion() {
        return None;
    }
    // Parallel or antiparallel axes both trace the same set, at whatever
    // phase and winding each was written with.
    let (za, zb) = (ca.frame().z().vector(), cb.frame().z().vector());
    if za.cross(zb).magnitude() > tol.angular() {
        return None;
    }
    // The ranges are a *correspondence*, which is what an overlap means and
    // what a caller carrying a split across the pair relies on: `on_b`'s ends
    // are the parameters at which `b` stands where `a`'s own ends do. Phase
    // comes from where `a` starts on `b`, winding from whether the two run
    // the same way there, and a pair written with opposite windings runs
    // `on_b` backwards, which is exactly the truth about them.
    let (lo, hi) = Curve3d::domain(a);
    let start = a.point_at(lo, tol).ok()?;
    let local = cb.frame().to_local(start);
    let angle = local.y.atan2(local.x);
    let phase = if b.is_reversed() { -angle } else { angle }.rem_euclid(core::f64::consts::TAU);
    let along_a = a.d1_at(lo, tol).ok()?;
    let along_b = b.d1_at(phase, tol).ok()?;
    let winding: f64 = if along_a.dot(along_b) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    Some(CurveIntersection {
        crossings: Vec::new(),
        overlaps: vec![Overlap {
            on_a: (lo, hi),
            on_b: (phase, winding.mul_add(hi - lo, phase)),
        }],
    })
}

/// The *same description* twice: one curve object meeting itself, forward
/// or reversed. A fitted seam reused as a wedge's apex ring is exactly this
/// pair, and the sampling path (every sample a hit) cannot answer it, for
/// the same reason it cannot answer coincident circles. Equality here is
/// structural, so two independent fits of one path still fall through to
/// the general machinery, which is the honest place for them.
fn same_curve_3d(a: &Curve, b: &Curve) -> Option<CurveIntersection<Point>> {
    let (lo, hi) = Curve3d::domain(a);
    if a == b {
        return Some(CurveIntersection {
            crossings: Vec::new(),
            overlaps: vec![Overlap {
                on_a: (lo, hi),
                on_b: (lo, hi),
            }],
        });
    }
    use ogeom_geom::Reversible as _;
    if *a == b.clone().reversed() {
        let (blo, bhi) = Curve3d::domain(b);
        return Some(CurveIntersection {
            crossings: Vec::new(),
            overlaps: vec![Overlap {
                on_a: (lo, hi),
                on_b: (bhi, blo),
            }],
        });
    }
    None
}

/// Two ellipses tracing the same point set in space: the ellipse counterpart
/// of [`same_circle_3d`], and just as invisible to the sampling path. Unlike
/// a circle, an ellipse's natural parameter is pinned to its major axis, so
/// the correspondence is affine only when the two `x` axes line up (parallel
/// or antiparallel) as well as the planes and radii; anything else falls
/// through to the general machinery.
fn same_ellipse_3d(
    a: &ogeom_geom::EllipseCurve,
    b: &ogeom_geom::EllipseCurve,
    tol: Tolerances,
) -> Option<CurveIntersection<Point>> {
    let (ea, eb) = (a.ellipse(), b.ellipse());
    if ea.frame().origin().distance(eb.frame().origin()) > tol.confusion() {
        return None;
    }
    if (ea.major_radius() - eb.major_radius()).abs() > tol.confusion()
        || (ea.minor_radius() - eb.minor_radius()).abs() > tol.confusion()
    {
        return None;
    }
    let (za, zb) = (ea.frame().z().vector(), eb.frame().z().vector());
    if za.cross(zb).magnitude() > tol.angular() {
        return None;
    }
    let (xa, xb) = (ea.frame().x().vector(), eb.frame().x().vector());
    if xa.cross(xb).magnitude() > tol.angular() {
        return None;
    }
    // As for circles: phase from where `a` starts on `b`, winding from
    // whether the two run the same way there, and the ranges come back as
    // the correspondence an overlap means.
    let (lo, hi) = Curve3d::domain(a);
    let start = a.point_at(lo, tol).ok()?;
    let angle = ogeom_math::elementary::ellipse_parameter(&eb, start, tol).ok()?;
    let phase = if b.is_reversed() { -angle } else { angle }.rem_euclid(core::f64::consts::TAU);
    let along_a = a.d1_at(lo, tol).ok()?;
    let along_b = b.d1_at(phase, tol).ok()?;
    let winding: f64 = if along_a.dot(along_b) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    Some(CurveIntersection {
        crossings: Vec::new(),
        overlaps: vec![Overlap {
            on_a: (lo, hi),
            on_b: (phase, winding.mul_add(hi - lo, phase)),
        }],
    })
}

fn check(options: CurveCurveOptions) -> OgeomResult<()> {
    if options.samples < 2 {
        ogeom_bail!(Construction, "seeding needs at least two segments");
    }
    if !options.gap.is_finite() || options.gap <= 0.0 {
        ogeom_bail!(Construction, "a gap of {} is not a distance", options.gap);
    }
    Ok(())
}

// --- analytic, planar --------------------------------------------------------

fn line_line_2d(
    a: &ogeom_geom::Line2d,
    b: &ogeom_geom::Line2d,
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
        // And back onto b.
        let back = |t: f64| (oa + da * t - ob).dot(db);
        if hi - lo <= tol.confusion() {
            // Segments meeting end to end share a point, not a stretch.
            if lo - hi > tol.confusion() {
                return CurveIntersection::empty();
            }
            let t = f64::midpoint(lo, hi).clamp(a_lo, a_hi);
            return CurveIntersection {
                crossings: vec![Crossing {
                    on_a: t,
                    on_b: back(t).clamp(b_lo, b_hi),
                    point: oa + da * t,
                    gap: 0.0,
                    reach: 0.0,
                }],
                overlaps: Vec::new(),
            };
        }
        return CurveIntersection {
            crossings: Vec::new(),
            overlaps: vec![Overlap {
                on_a: (lo, hi),
                // Paired end to end with `on_a`, not sorted: two lines written in
                // opposite directions run `on_b` backwards, and a consumer
                // carrying a stretch across by the correspondence (the
                // boolean clipping a contact to the edge it runs along)
                // reads a sorted pair as the reflected stretch.
                on_b: (back(lo), back(hi)),
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
            reach: 0.0,
        }],
        overlaps: Vec::new(),
    }
}

fn line_circle_2d(
    line: &ogeom_geom::Line2d,
    circle: &ogeom_geom::Circle2d,
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
            reach: 0.0,
        });
    }
    sort_crossings(&mut crossings);
    CurveIntersection {
        crossings,
        overlaps: Vec::new(),
    }
}

fn circle_circle_2d(
    a: &ogeom_geom::Circle2d,
    b: &ogeom_geom::Circle2d,
    tol: Tolerances,
) -> CurveIntersection<Point2> {
    let (ca, cb) = (a.circle(), b.circle());
    let between = cb.centre() - ca.centre();
    let distance = between.magnitude();
    let (ra, rb) = (ca.radius(), cb.radius());

    if distance <= tol.confusion() {
        if (ra - rb).abs() <= tol.confusion() {
            // The same circle, the whole turn of each. As in space: phase
            // from where `a` starts on `b`, winding from whether the two run
            // the same way there, so `on_b` names where `b` stands at
            // `a`'s own ends.
            use ogeom_geom::Curve2d as _;
            let (lo, hi) = a.domain();
            let correspondence = (|| {
                let start = a.point_at(lo, tol).ok()?;
                let phase = circle_parameter(b, start, tol)?;
                let along_a = a.d1_at(lo, tol).ok()?;
                let along_b = b.d1_at(phase, tol).ok()?;
                let winding: f64 = if along_a.dot(along_b) >= 0.0 {
                    1.0
                } else {
                    -1.0
                };
                Some((phase, winding.mul_add(hi - lo, phase)))
            })();
            return CurveIntersection {
                crossings: Vec::new(),
                overlaps: vec![Overlap {
                    on_a: (lo, hi),
                    on_b: correspondence.unwrap_or_else(|| b.domain()),
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
                reach: 0.0,
            });
        }
    };
    if squared <= tol.confusion() * tol.confusion() {
        push(foot);
    } else {
        let offset = ogeom_math::Vector2::new(-direction.y, direction.x) * squared.max(0.0).sqrt();
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
fn circle_parameter(curve: &ogeom_geom::Circle2d, p: Point2, tol: Tolerances) -> Option<f64> {
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
    a: &ogeom_geom::LineCurve,
    b: &ogeom_geom::LineCurve,
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
        let back = |t: f64| (oa + da * t - ob).dot(db);
        if hi - lo <= tol.confusion() {
            // Segments meeting end to end share a point, not a stretch.
            if lo - hi > tol.confusion() {
                return CurveIntersection::empty();
            }
            let t = f64::midpoint(lo, hi).clamp(a_lo, a_hi);
            let s = back(t).clamp(b_lo, b_hi);
            let point = oa + da * t;
            return CurveIntersection {
                crossings: vec![Crossing {
                    on_a: t,
                    on_b: s,
                    point,
                    gap: point.distance(ob + db * s),
                    reach: 0.0,
                }],
                overlaps: Vec::new(),
            };
        }
        return CurveIntersection {
            crossings: Vec::new(),
            overlaps: vec![Overlap {
                on_a: (lo, hi),
                // Paired end to end with `on_a`, not sorted: two lines written in
                // opposite directions run `on_b` backwards, and a consumer
                // carrying a stretch across by the correspondence (the
                // boolean clipping a contact to the edge it runs along)
                // reads a sorted pair as the reflected stretch.
                on_b: (back(lo), back(hi)),
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
            reach: 0.0,
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
) -> OgeomResult<CurveIntersection<Point2>> {
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
) -> OgeomResult<CurveIntersection<Point>> {
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

    // Where the first curve runs along the second, every segment pair near
    // the stretch is in reach, and each seeds a polish that lands on the
    // same shared support the overlap pass reports whole and removes the
    // crossings of. So the samples' feet are found first, and a segment of
    // the first curve lying wholly in such a stretch seeds nothing: a
    // section along an edge would otherwise run thousands of polishes for
    // crossings that are all discarded.
    let feet: Vec<Option<(f64, f64)>> = sa
        .points
        .iter()
        .map(|p| foot_via_samples(b, &sb, *p, tol))
        .collect();
    let overlaps = if hugging_runs(&feet, options.gap).contains(&true) {
        shared_support_3d(a, b, &sa, &sb, &feet, options, tol)
    } else {
        Vec::new()
    };
    let shared = |lo: f64, hi: f64| {
        overlaps.iter().any(|o| {
            let (from, to) = order(o.on_a.0, o.on_a.1);
            lo >= from && hi <= to
        })
    };

    // Every segment pair's closest approach, then seeds only where it is a
    // local minimum among the neighbouring pairs. Each crossing sits in
    // such a basin; the pairs round it would all polish to the same place.
    // Two curves running near each other (a section beside the edge it was
    // cut along) put hundreds of pairs in reach, and polishing every one
    // found the same few crossings over and over: five in six polishes
    // were repeats.
    let (na, nb) = (sa.points.len(), sb.points.len());
    // Pairs whose boxes stand further apart than the reach cannot be in
    // it: they are not measured, and read as out of reach, which is what
    // measuring them would have said, for the reach test and for the
    // basin test alike.
    let boxes = |s: &Sampled<Point>| -> Vec<(Point, Point)> {
        s.points
            .windows(2)
            .map(|w| {
                (
                    Point::new(w[0].x.min(w[1].x), w[0].y.min(w[1].y), w[0].z.min(w[1].z)),
                    Point::new(w[0].x.max(w[1].x), w[0].y.max(w[1].y), w[0].z.max(w[1].z)),
                )
            })
            .collect()
    };
    let (boxes_a, boxes_b) = (boxes(&sa), boxes(&sb));
    let apart = |x: &(Point, Point), y: &(Point, Point)| {
        x.0.x - y.1.x > reach
            || y.0.x - x.1.x > reach
            || x.0.y - y.1.y > reach
            || y.0.y - x.1.y > reach
            || x.0.z - y.1.z > reach
            || y.0.z - x.1.z > reach
    };
    let mut approach: Vec<(f64, f64, f64)> =
        Vec::with_capacity(na.saturating_sub(1) * nb.saturating_sub(1));
    for i in 1..na {
        for j in 1..nb {
            approach.push(if apart(&boxes_a[i - 1], &boxes_b[j - 1]) {
                (0.0, 0.0, f64::INFINITY)
            } else {
                segments_approach_3d(
                    (sa.points[i - 1], sa.points[i]),
                    (sb.points[j - 1], sb.points[j]),
                )
            });
        }
    }
    let cols = nb.saturating_sub(1);
    let gap_at = |i: usize, j: usize| approach[(i - 1) * cols + (j - 1)].2;
    // Lowest along either curve: a pair nearer than both its neighbours in
    // the first curve's direction, or in the second's.
    let basin = |i: usize, j: usize| {
        let here = gap_at(i, j);
        let along_a = (i.saturating_sub(1).max(1)..=(i + 1).min(na - 1))
            .all(|p| p == i || here <= gap_at(p, j));
        let along_b = (j.saturating_sub(1).max(1)..=(j + 1).min(nb - 1))
            .all(|q| q == j || here <= gap_at(i, q));
        along_a || along_b
    };
    let mut crossings: Vec<Crossing<Point>> = Vec::new();
    for i in 1..na {
        let (lo, hi) = order(sa.parameters[i - 1], sa.parameters[i]);
        if shared(lo, hi) {
            continue;
        }
        for j in 1..nb {
            let (ta, tb, gap) = approach[(i - 1) * cols + (j - 1)];
            if gap > reach || !basin(i, j) {
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

    // A tangential contact is one crossing, however many the polish
    // returns. Where two curves touch, the stationarity conditions go flat
    // along the contact: every seed converges somewhere in a valley the
    // width of the gap, and an arc ending on the line it is tangent to
    // comes back as thirty crossings inside a micron or two. Consecutive
    // crossings with the first curve staying within the gap of the second
    // all the way between them are the same contact, and the nearest
    // approach among them speaks for it.
    if crossings.len() > 1 {
        let mut merged: Vec<Crossing<Point>> = Vec::with_capacity(crossings.len());
        let mut run_start: Option<Point> = None;
        for c in crossings {
            if let Some(last) = merged.last_mut()
                && contact_between_3d(a, b, last, &c, options, tol)
            {
                let start = run_start.get_or_insert(last.point);
                let reach = start.distance(c.point).max(last.reach);
                if c.gap < last.gap {
                    *last = c;
                }
                last.reach = reach;
                continue;
            }
            run_start = None;
            merged.push(c);
        }
        // A touch astride the first curve's period seam comes back as a
        // crossing at each end of the parameter range; the two are one
        // contact as well.
        if merged.len() > 1 && a.is_periodic() {
            let (lo, hi) = a.domain();
            let (first, last) = (merged[0], merged[merged.len() - 1]);
            let wrapped = Crossing {
                on_a: first.on_a + (hi - lo),
                ..first
            };
            if contact_between_3d(a, b, &last, &wrapped, options, tol) {
                let reach = last
                    .reach
                    .max(first.reach)
                    .max(last.point.distance(first.point));
                let keep = if first.gap <= last.gap {
                    0
                } else {
                    merged.len() - 1
                };
                merged[keep].reach = reach;
                if keep == 0 {
                    merged.pop();
                } else {
                    merged.remove(0);
                }
            }
        }
        crossings = merged;
    }

    // Stretches where the first curve stays within the gap of the second
    // are shared support, not a row of crossings. A fitted section tracing
    // the arc it was cut along wobbles about it by less than the gap and
    // "crosses" it at every wobble; read as crossings, those shatter the
    // curve into hundreds of pieces and pave the edge at each. So every
    // sample of the first curve asks its foot on the second, a run of
    // consecutive samples within the gap is an overlap with its ends
    // bisected to parametric resolution, and the crossings inside it are
    // the overlap's, not the caller's.
    if !overlaps.is_empty() {
        crossings.retain(|c| {
            !overlaps.iter().any(|o| {
                let (lo, hi) = order(o.on_a.0, o.on_a.1);
                c.on_a >= lo - tol.parametric() && c.on_a <= hi + tol.parametric()
            })
        });
    }
    // Every surviving crossing owns the valley it sits in: how far along the
    // first curve the second stays within the caller's gap. A transversal
    // crossing leaves the gap within a gap's length and says nothing; a
    // tangential one (a line touching a fitted rim that wobbles about its
    // circle by the fit's budget) stays inside for the root of gap times
    // radius on either side, and the polish lands on whichever wobble's
    // floor it found. The consumer placing a vertex there owns that much
    // doubt, which the spread of several polished crossings only stated
    // when there were several. A valley longer than a tangency's (the
    // radius being at most the shorter curve's length) is a shared stretch
    // the overlap pass speaks for, not a crossing's to own.
    let gap = options.gap.max(tol.confusion());
    let extent = {
        let along = |s: &Sampled<Point>| {
            s.points
                .windows(2)
                .map(|w| w[0].distance(w[1]))
                .sum::<f64>()
        };
        along(&sa).min(along(&sb))
    };
    let cap = 4.0 * (gap * extent).sqrt();
    for c in &mut crossings {
        let valley = valley_extent_3d(a, b, &sb, c, gap, tol);
        if valley > gap * 8.0 && valley <= cap {
            c.reach = c.reach.max(valley);
        }
    }
    Ok(CurveIntersection {
        crossings,
        overlaps,
    })
}

/// Whether the first curve stays within the gap of the second all the way
/// from one crossing to the next: three stations between them, each foot
/// seeded from the crossings' own parameters.
fn contact_between_3d(
    a: &Curve,
    b: &Curve,
    from: &Crossing<Point>,
    to: &Crossing<Point>,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> bool {
    if (to.on_a - from.on_a).abs() <= tol.parametric() {
        return true;
    }
    (1..=3).all(|k| {
        let f = f64::from(k) / 4.0;
        let t = from.on_a + (to.on_a - from.on_a) * f;
        let seed = from.on_b + (to.on_b - from.on_b) * f;
        a.point_at(t, tol)
            .ok()
            .and_then(|p| foot_on_3d(b, p, seed, tol))
            .is_some_and(|(_, gap)| gap <= options.gap)
    })
}

/// Runs of the first curve's samples whose feet on the second lie within
/// the gap, each bisected to its parametric ends.
/// The foot of a point on the second curve, seeded from its sampled
/// polyline's nearest segment: the parameter and the distance there.
fn foot_via_samples(
    b: &Curve,
    sb: &Sampled<Point>,
    p: Point,
    tol: Tolerances,
) -> Option<(f64, f64)> {
    let mut seed = (f64::INFINITY, 0.0);
    for j in 1..sb.points.len() {
        let (_, tb, gap) = segments_approach_3d((p, p), (sb.points[j - 1], sb.points[j]));
        if gap < seed.0 {
            seed = (
                gap,
                sb.parameters[j - 1] + (sb.parameters[j] - sb.parameters[j - 1]) * tb,
            );
        }
    }
    if !seed.0.is_finite() {
        return None;
    }
    foot_on_3d(b, p, seed.1, tol)
}

/// Which samples of the first curve lie in a run of two or more whose feet
/// on the second are within the gap: the stretches [`shared_support_3d`]
/// reports as overlaps, sample by sample.
fn hugging_runs(feet: &[Option<(f64, f64)>], gap: f64) -> Vec<bool> {
    let within = |i: usize| feet[i].is_some_and(|(_, g)| g <= gap);
    let mut out = vec![false; feet.len()];
    let mut i = 0;
    while i < feet.len() {
        if !within(i) {
            i += 1;
            continue;
        }
        let start = i;
        while i + 1 < feet.len() && within(i + 1) {
            i += 1;
        }
        if i > start {
            out[start..=i].fill(true);
        }
        i += 1;
    }
    out
}

fn shared_support_3d(
    a: &Curve,
    b: &Curve,
    sa: &Sampled<Point>,
    sb: &Sampled<Point>,
    feet: &[Option<(f64, f64)>],
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Vec<Overlap> {
    let hugs = |t: f64| -> Option<(f64, f64)> {
        let p = a.point_at(t, tol).ok()?;
        let (s, gap) = foot_via_samples(b, sb, p, tol)?;
        (gap <= options.gap).then_some((s, gap))
    };
    let within = |i: usize| feet[i].is_some_and(|(_, gap)| gap <= options.gap);

    let mut overlaps = Vec::new();
    let mut i = 0;
    while i < sa.points.len() {
        if !within(i) {
            i += 1;
            continue;
        }
        let start = i;
        while i + 1 < sa.points.len() && within(i + 1) {
            i += 1;
        }
        let end = i;
        i += 1;
        if end == start {
            continue;
        }
        // The run's ends: where the samples stop hugging, bisected between
        // the last inside sample and the first outside one.
        let refine = |inside: usize, outside: Option<usize>| -> (f64, f64) {
            let (mut t_in, s_in) = (sa.parameters[inside], feet[inside].map_or(0.0, |f| f.0));
            let Some(out) = outside else {
                return (t_in, s_in);
            };
            let mut s_at = s_in;
            let mut t_out = sa.parameters[out];
            for _ in 0..48 {
                if (t_out - t_in).abs() <= tol.parametric() {
                    break;
                }
                let mid = f64::midpoint(t_in, t_out);
                match hugs(mid) {
                    Some((s, _)) => {
                        t_in = mid;
                        s_at = s;
                    }
                    None => t_out = mid,
                }
            }
            (t_in, s_at)
        };
        let (lo_a, lo_b) = refine(start, start.checked_sub(1));
        let (hi_a, hi_b) = refine(end, (end + 1 < sa.points.len()).then_some(end + 1));
        if hi_a - lo_a <= tol.parametric() || (hi_b - lo_b).abs() <= tol.parametric() {
            continue;
        }
        overlaps.push(Overlap {
            on_a: (lo_a, hi_a),
            on_b: (lo_b, hi_b),
        });
    }
    overlaps
}

/// Newton on the foot-point condition `(c(s) - p) . c'(s) = 0` from a seed,
/// clamped to the curve's domain; the parameter and the distance there.
fn foot_on_3d(curve: &Curve, p: Point, seed: f64, tol: Tolerances) -> Option<(f64, f64)> {
    let mut s = clamp_3d(curve, seed);
    let mut best = (s, curve.point_at(s, tol).ok()?.distance(p));
    for _ in 0..30 {
        let d = curve.derivatives_at(s, 2, tol).ok()?;
        let zero = ogeom_math::Vector::ZERO;
        let (c, d1, d2) = (
            d.first().copied().unwrap_or(zero),
            d.get(1).copied().unwrap_or(zero),
            d.get(2).copied().unwrap_or(zero),
        );
        let gap = c - (p - Point::ORIGIN);
        let g = gap.dot(d1);
        let dg = d1.dot(d1) + gap.dot(d2);
        if dg.abs() <= f64::MIN_POSITIVE {
            break;
        }
        let next = clamp_3d(curve, s - g / dg);
        let dist = curve.point_at(next, tol).ok()?.distance(p);
        let moved = (next - s).abs();
        s = next;
        if dist < best.1 {
            best = (s, dist);
        }
        if moved <= tol.parametric() {
            break;
        }
    }
    Some(best)
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
    let system = |x: &[f64; 2]| {
        let (t, s) = (clamp_2d(a, x[0]), clamp_2d(b, x[1]));
        let pa = a.point_at(t, tol).unwrap_or(Point2::ORIGIN);
        let pb = b.point_at(s, tol).unwrap_or(Point2::ORIGIN);
        let da = a
            .d1_at(t, tol)
            .unwrap_or(ogeom_math::Vector2::new(0.0, 0.0));
        let db = b
            .d1_at(s, tol)
            .unwrap_or(ogeom_math::Vector2::new(0.0, 0.0));
        ([pa.x - pb.x, pa.y - pb.y], [[da.x, -db.x], [da.y, -db.y]])
    };
    let criteria = solve::Criteria {
        residual: tol.confusion() * 0.01,
        step: tol.parametric(),
        max_iterations: 40,
    };
    let found = solve::newton_system_fixed(system, [seed_a, seed_b], criteria).ok()?;
    let (t, s) = (clamp_2d(a, found.0[0]), clamp_2d(b, found.0[1]));
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
        reach: 0.0,
    })
}

/// Gauss–Newton on the closest approach of two space curves.
///
/// Three equations would be overdetermined for two unknowns, so the system is
/// the two *stationarity* conditions (the gap vector perpendicular to both
/// tangents), whose solutions are the local closest approaches. The gap test
/// afterwards decides whether the approach found is a crossing.
fn polish_3d(
    a: &Curve,
    b: &Curve,
    seed_a: f64,
    seed_b: f64,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Option<Crossing<Point>> {
    // Two unknowns: the allocation-free solver, and each curve's point
    // read from the same derivative table as its derivatives.
    let system = |x: [f64; 2]| {
        let (t, s) = (clamp_3d(a, x[0]), clamp_3d(b, x[1]));
        let (Ok(da), Ok(db)) = (a.derivatives_at(t, 2, tol), b.derivatives_at(s, 2, tol)) else {
            // A zero here would read as a root; infinite, the damped step
            // backs off instead.
            return ([f64::INFINITY; 2], [[0.0; 2]; 2]);
        };
        let zero = ogeom_math::Vector::ZERO;
        let at = |d: &[ogeom_math::Vector], k: usize| d.get(k).copied().unwrap_or(zero);
        let (pa, d1a, d2a) = (at(&da, 0), at(&da, 1), at(&da, 2));
        let (pb, d1b, d2b) = (at(&db, 0), at(&db, 1), at(&db, 2));
        let gap = pa - pb;
        (
            [gap.dot(d1a), -gap.dot(d1b)],
            [
                [d1a.dot(d1a) + gap.dot(d2a), -d1a.dot(d1b)],
                [-d1a.dot(d1b), d1b.dot(d1b) - gap.dot(d2b)],
            ],
        )
    };
    let criteria = solve::Criteria {
        residual: tol.confusion() * 0.01,
        step: tol.parametric(),
        max_iterations: 40,
    };
    let ([t, s], ..) = solve::newton_system_2(system, [seed_a, seed_b], criteria).ok()?;
    let (t, s) = (clamp_3d(a, t), clamp_3d(b, s));
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
        reach: 0.0,
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

/// The distance from `p` to the curve `b`, through its samples and a local
/// polish on the nearest segment's parameter span.
fn distance_to_curve_3d(b: &Curve, sb: &Sampled<Point>, p: Point, tol: Tolerances) -> f64 {
    let mut best = (0_usize, f64::INFINITY);
    for i in 1..sb.points.len() {
        let (q0, q1) = (sb.points[i - 1], sb.points[i]);
        let d = q1 - q0;
        let len2 = d.dot(d);
        let f = if len2 <= f64::MIN_POSITIVE {
            0.0
        } else {
            ((p - q0).dot(d) / len2).clamp(0.0, 1.0)
        };
        let dist = p.distance(q0 + d * f);
        if dist < best.1 {
            best = (i, dist);
        }
    }
    if best.0 == 0 {
        return best.1;
    }
    let (mut lo, mut hi) = (sb.parameters[best.0 - 1], sb.parameters[best.0]);
    let at = |t: f64| -> f64 { b.point_at(t, tol).map_or(f64::INFINITY, |q| q.distance(p)) };
    // Golden-section on the segment's span: the distance is unimodal there
    // at any sampling that resolved the curve at all.
    let phi = 0.5 * (3.0 - 5.0_f64.sqrt());
    let (mut x1, mut x2) = (lo + phi * (hi - lo), hi - phi * (hi - lo));
    let (mut f1, mut f2) = (at(x1), at(x2));
    for _ in 0..48 {
        if f1 < f2 {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = lo + phi * (hi - lo);
            f1 = at(x1);
        } else {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = hi - phi * (hi - lo);
            f2 = at(x2);
        }
    }
    f1.min(f2).min(best.1)
}

/// How far from a crossing, along the first curve, the second curve stays
/// within `gap`: the larger of the two directions, in space.
fn valley_extent_3d(
    a: &Curve,
    b: &Curve,
    sb: &Sampled<Point>,
    crossing: &Crossing<Point>,
    gap: f64,
    tol: Tolerances,
) -> f64 {
    let (lo, hi) = a.domain();
    let span = hi - lo;
    if span <= 0.0 {
        return 0.0;
    }
    let mut extent = 0.0_f64;
    for direction in [-1.0, 1.0] {
        let inside = |t: f64| -> bool {
            if t < lo || t > hi {
                return false;
            }
            a.point_at(t, tol)
                .is_ok_and(|q| distance_to_curve_3d(b, sb, q, tol) <= gap)
        };
        let mut step = span * 1e-6;
        let mut last_in = crossing.on_a;
        let mut first_out: Option<f64> = None;
        while step <= span {
            let t = crossing.on_a + direction * step;
            if inside(t) {
                last_in = t;
                step *= 2.0;
            } else {
                first_out = Some(t);
                break;
            }
        }
        let edge = match first_out {
            Some(mut out) => {
                let mut r#in = last_in;
                for _ in 0..30 {
                    let mid = f64::midpoint(r#in, out);
                    if inside(mid) {
                        r#in = mid;
                    } else {
                        out = mid;
                    }
                }
                r#in
            }
            None => last_in,
        };
        if let Ok(q) = a.point_at(edge, tol) {
            extent = extent.max(q.distance(crossing.point));
        }
    }
    extent
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
    use ogeom_geom::{BSpline2d, Circle2d, CircleCurve, Line2d, LineCurve};
    use ogeom_math::{Circle, Circle2, Direction2, Frame, Frame2, KnotVector, Vector2};

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

    /// One circle written twice in the plane, a quarter turn apart and
    /// wound against each other: the overlap's ranges correspond.
    #[test]
    fn one_planar_circle_written_twice_states_the_correspondence() {
        let a = circle2(Point2::new(1.0, 2.0), 3.0);
        let quarter = Frame2::new(
            Point2::new(1.0, 2.0),
            Direction2::new(Vector2::new(0.0, 1.0), T).unwrap(),
        );
        let b: PlanarCurve = Circle2d::new(Circle2::new(quarter, 3.0, T).unwrap()).into();
        let flipped: PlanarCurve = ogeom_geom::Reversible::reversed(&b);
        for other in [b, flipped] {
            let found = intersect_curves_2d(&a, &other, CurveCurveOptions::default(), T).unwrap();
            assert_eq!(found.overlaps.len(), 1);
            let overlap = found.overlaps[0];
            for i in 0..=8 {
                let t = f64::from(i) / 8.0;
                let ta = (overlap.on_a.1 - overlap.on_a.0).mul_add(t, overlap.on_a.0);
                let tb = (overlap.on_b.1 - overlap.on_b.0).mul_add(t, overlap.on_b.0);
                let pa = a.point_at(ta, T).unwrap();
                let pb = other
                    .point_at(tb.rem_euclid(core::f64::consts::TAU), T)
                    .unwrap();
                assert!(pa.distance(pb) < 1e-9, "at {t}: {pa:?} against {pb:?}");
            }
        }
    }

    /// Collinear segments meeting end to end share their end point, in the
    /// plane and in space, as a perpendicular pair meeting there does.
    #[test]
    fn collinear_segments_meeting_end_to_end_share_a_point() {
        let a = line2(Point2::new(0.0, 0.0), Point2::new(5.0, 0.0));
        let b = line2(Point2::new(5.0, 0.0), Point2::new(9.0, 0.0));
        let found = intersect_curves_2d(&a, &b, CurveCurveOptions::default(), T).unwrap();
        assert!(found.overlaps.is_empty());
        assert_eq!(found.crossings.len(), 1);
        assert!(found.crossings[0].point.distance(Point2::new(5.0, 0.0)) < 1e-12);

        let a: Curve = LineCurve::segment(Point::ORIGIN, Point::new(0.0, 0.0, 5.0), T)
            .unwrap()
            .into();
        let b: Curve = LineCurve::segment(Point::new(0.0, 0.0, 5.0), Point::new(0.0, 0.0, 7.0), T)
            .unwrap()
            .into();
        let found = intersect_curves(&a, &b, CurveCurveOptions::default(), T).unwrap();
        assert!(found.overlaps.is_empty());
        assert_eq!(found.crossings.len(), 1);
        assert!(found.crossings[0].point.distance(Point::new(0.0, 0.0, 5.0)) < 1e-12);
    }

    /// Two descriptions of one circle overlap over the whole turn, and the
    /// overlap's two ranges *correspond*: `on_b`'s ends are where `b` stands
    /// at `a`'s own ends. A caller carrying a split from one to the other
    /// (the boolean, pairing a hole's arcs against the disc that fills them)
    /// reads that correspondence and gets the same point back, whatever phase
    /// and winding the two were written with.
    #[test]
    fn one_circle_written_twice_states_the_correspondence_between_them() {
        use ogeom_math::{Direction, Vector};
        let a: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 3.0, T).unwrap()).into();
        // The same circle seen from underneath, started a third of a turn
        // round: opposite winding, arbitrary phase.
        let third = 2.0 * core::f64::consts::PI / 3.0;
        let flipped = Frame::new(
            Point::ORIGIN,
            -Direction::Z,
            Direction::new(Vector::new(third.cos(), third.sin(), 0.0), T).unwrap(),
            T,
        )
        .unwrap();
        let b: Curve = CircleCurve::new(Circle::new(flipped, 3.0, T).unwrap()).into();
        for pair in [(&a, &b), (&b, &a)] {
            let found = intersect_curves(pair.0, pair.1, CurveCurveOptions::default(), T).unwrap();
            assert!(
                found.crossings.is_empty(),
                "every point is a hit, so none is"
            );
            assert_eq!(found.overlaps.len(), 1);
            let overlap = &found.overlaps[0];
            let span = overlap.on_a.1 - overlap.on_a.0;
            for i in 0..=8 {
                let t = f64::from(i) / 8.0;
                let ta = span.mul_add(t, overlap.on_a.0);
                let tb = (overlap.on_b.1 - overlap.on_b.0).mul_add(t, overlap.on_b.0);
                let pa = pair.0.point_at(ta, T).unwrap();
                let pb = pair
                    .1
                    .point_at(tb.rem_euclid(core::f64::consts::TAU), T)
                    .unwrap();
                assert!(
                    pa.distance(pb) < 1e-9,
                    "at {t}: {pa:?} against {pb:?} (on_a {:?} on_b {:?})",
                    overlap.on_a,
                    overlap.on_b
                );
            }
        }
    }

    #[test]
    fn space_curves_cross_within_a_gap_and_report_it() {
        // Two circles that would cross in a shared plane, with one lifted a
        // hair out of it: the crossings become passes with a real, small gap
        // that must be reported, not zeroed. (Not chain links: a first draft
        // of this test used linked circles, and linked circles never approach:
        // passing through each other's *disks* is what linked means, and these
        // radii hold the curves a constant two units apart.)
        let a: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.0, T).unwrap()).into();
        let lifted = Frame::new(
            Point::new(3.0, 0.0, 0.001),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
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
    fn a_fitted_curve_tracing_an_arc_is_one_overlap_not_a_row_of_crossings() {
        // A spline fitted along a circle's arc sits within its fit budget
        // of the circle everywhere, and "crosses" it at every wobble. The
        // sampling path reports the stretch as one overlap and keeps no
        // crossing inside it; read as crossings, a section tracing the arc
        // it was cut along shattered into hundreds of pieces.
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 4.0, T).unwrap()).into();
        let points: Vec<Point> = (0..=40)
            .map(|i| {
                let a = 0.2 + 1.0 * f64::from(i) / 40.0;
                Point::new(4.0 * a.cos(), 4.0 * a.sin(), 0.0)
            })
            .collect();
        let fitted: Curve = ogeom_geom::fit::fit_points(&points, 3, 1e-7, T)
            .unwrap()
            .curve
            .into();
        let options = CurveCurveOptions {
            gap: 1e-5,
            ..CurveCurveOptions::default()
        };
        let found = intersect_curves(&fitted, &circle, options, T).unwrap();
        assert_eq!(found.overlaps.len(), 1, "one shared stretch: {found:?}");
        let (lo, hi) = found.overlaps[0].on_a;
        let (fa, fb) = fitted.domain();
        assert!(
            lo - fa < 1e-3 && fb - hi < 1e-3,
            "the whole fit runs along the circle"
        );
        assert!(
            found.crossings.is_empty(),
            "no crossing survives inside the overlap: {:?}",
            found.crossings
        );
    }

    #[test]
    fn an_arc_ending_tangent_to_a_line_is_one_crossing_with_its_reach() {
        // A circle and its tangent line touch at one point, but the
        // stationarity conditions go flat along the touch and every seed
        // converges somewhere in a valley the width of the gap. One contact
        // comes back, at the touch, owning the valley's length as its reach.
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 4.0, T).unwrap()).into();
        let line: Curve =
            LineCurve::segment(Point::new(4.0, -3.0, 0.0), Point::new(4.0, 3.0, 0.0), T)
                .unwrap()
                .into();
        let options = CurveCurveOptions {
            gap: 1e-5,
            ..CurveCurveOptions::default()
        };
        let found = intersect_curves(&circle, &line, options, T).unwrap();
        assert_eq!(found.crossings.len(), 1, "one touch: {:?}", found.crossings);
        let touch = found.crossings[0];
        assert!(
            touch.point.distance(Point::new(4.0, 0.0, 0.0)) < 2e-2,
            "{touch:?}"
        );
        assert!(touch.reach < 5e-2, "the valley is short: {touch:?}");
        assert!(found.overlaps.is_empty());
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

    /// Plane sections of one cylinder: a circle and two ellipses tilted
    /// different ways. Any two meet where their planes' common line pierces
    /// the cylinder, twice, and a circle lifted parallel to another meets
    /// it nowhere.
    #[test]
    fn circles_and_ellipses_in_different_planes_meet_where_the_planes_do() {
        use ogeom_geom::EllipseCurve;
        use ogeom_math::{Direction, Ellipse, Vector};
        let radius = 2.0;
        let tilted = |normal: Vector, major: Vector, lean: f64| -> Curve {
            let frame = Frame::new(
                Point::ORIGIN,
                Direction::new(normal, T).unwrap(),
                Direction::new(major, T).unwrap(),
                T,
            )
            .unwrap();
            EllipseCurve::new(Ellipse::new(frame, radius / lean.cos(), radius, T).unwrap()).into()
        };
        let (p, q) = (0.4_f64, 0.7_f64);
        let about_x = tilted(
            Vector::new(0.0, -p.sin(), p.cos()),
            Vector::new(0.0, p.cos(), p.sin()),
            p,
        );
        let about_y = tilted(
            Vector::new(-q.sin(), 0.0, q.cos()),
            Vector::new(q.cos(), 0.0, q.sin()),
            q,
        );
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, radius, T).unwrap()).into();
        for (a, b) in [
            (&about_x, &about_y),
            (&about_y, &about_x),
            (&circle, &about_x),
            (&about_y, &circle),
        ] {
            let found = intersect_curves(a, b, CurveCurveOptions::default(), T).unwrap();
            assert_eq!(found.crossings.len(), 2, "{found:?}");
            for hit in &found.crossings {
                let on_a = a.point_at(hit.on_a, T).unwrap();
                let on_b = b.point_at(hit.on_b, T).unwrap();
                assert!(on_a.distance(on_b) < 1e-9, "{on_a:?} against {on_b:?}");
                assert!((on_a.x.hypot(on_a.y) - radius).abs() < 1e-9);
            }
        }
        let lifted = Frame::new(Point::new(0.0, 0.0, 1.0), Direction::Z, Direction::X, T).unwrap();
        let above: Curve = CircleCurve::new(Circle::new(lifted, radius, T).unwrap()).into();
        assert!(
            intersect_curves(&circle, &above, CurveCurveOptions::default(), T)
                .unwrap()
                .is_empty()
        );
    }

    /// Two circles in planes a millionth of a radian apart, passing within
    /// the gap of each other where their shadows cross: the first crosses
    /// the second's plane a twentieth of a unit away from there, where the
    /// two are far apart. The near pass is found all the same.
    #[test]
    fn circles_in_all_but_one_plane_meet_where_they_pass() {
        use ogeom_math::{Direction, Vector};
        let flat: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.0, T).unwrap()).into();
        let lean = 1e-6;
        let normal = Direction::new(Vector::new(lean, 0.0, 1.0), T).unwrap();
        let centre = Point::new(1.0, 0.0, lean.mul_add(-0.5, 5e-8));
        let tilted_frame = Frame::new(
            centre,
            normal,
            Direction::new(Vector::new(1.0, 0.0, -lean), T).unwrap(),
            T,
        )
        .unwrap();
        let tilted: Curve = CircleCurve::new(Circle::new(tilted_frame, 2.0, T).unwrap()).into();
        let found = intersect_curves(&flat, &tilted, CurveCurveOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 2, "{found:?}");
        for hit in &found.crossings {
            assert!(hit.gap < 1e-7, "{hit:?}");
            let p = flat.point_at(hit.on_a, T).unwrap();
            assert!((p.x - 0.5).abs() < 1e-4, "{p:?}");
        }
    }

    /// A line across an ellipse in its plane meets it twice, a line through
    /// a circle's plane meets it once where it pierces the circle, and a
    /// line through the plane inside the circle misses it. Either order of
    /// the pair answers the same, parameters swapped.
    #[test]
    fn lines_meet_circles_and_ellipses_in_closed_form() {
        use ogeom_geom::EllipseCurve;
        use ogeom_math::Ellipse;
        let ellipse: Curve =
            EllipseCurve::new(Ellipse::new(Frame::WORLD, 3.0, 2.0, T).unwrap()).into();
        let across: Curve =
            LineCurve::segment(Point::new(-5.0, 1.0, 0.0), Point::new(5.0, 1.0, 0.0), T)
                .unwrap()
                .into();
        for (a, b, line_first) in [(&ellipse, &across, false), (&across, &ellipse, true)] {
            let found = intersect_curves(a, b, CurveCurveOptions::default(), T).unwrap();
            assert_eq!(found.crossings.len(), 2, "{found:?}");
            for hit in &found.crossings {
                let p = a.point_at(hit.on_a, T).unwrap();
                let q = b.point_at(hit.on_b, T).unwrap();
                assert!(p.distance(q) < 1e-9, "{p:?} against {q:?}");
                let on_line = if line_first { p } else { q };
                assert!((on_line.y - 1.0).abs() < 1e-12);
            }
        }
        let circle: Curve = CircleCurve::new(Circle::new(Frame::WORLD, 2.0, T).unwrap()).into();
        let through: Curve =
            LineCurve::segment(Point::new(2.0, 0.0, -1.0), Point::new(2.0, 0.0, 1.0), T)
                .unwrap()
                .into();
        let found = intersect_curves(&circle, &through, CurveCurveOptions::default(), T).unwrap();
        assert_eq!(found.crossings.len(), 1);
        assert!(found.crossings[0].on_a.abs() < 1e-9);
        assert!((found.crossings[0].on_b - 1.0).abs() < 1e-9);
        let inside: Curve =
            LineCurve::segment(Point::new(1.0, 0.0, -1.0), Point::new(1.0, 0.0, 1.0), T)
                .unwrap()
                .into();
        assert!(
            intersect_curves(&circle, &inside, CurveCurveOptions::default(), T)
                .unwrap()
                .is_empty()
        );
    }

    /// Two trims of one circle share whatever stretch they share, across the
    /// seam or not, walked the same way or the other: every piece comes
    /// back inside both windows, and each end of it is the same point on
    /// both curves.
    #[test]
    fn arcs_of_one_circle_overlap_across_the_seam() {
        use ogeom_geom::{Curve3d as _, TrimmedCurve};
        let circle = |x: ogeom_math::Direction, z: ogeom_math::Direction| -> Curve {
            let frame = Frame::new(Point::new(10.0, 20.0, 30.0), z, x, T).unwrap();
            CircleCurve::new(Circle::new(frame, 50.0, T).unwrap()).into()
        };
        let trim = |c: &Curve, a: f64, b: f64| -> Curve {
            TrimmedCurve::new(c.clone(), a, b, T).unwrap().into()
        };
        let (x, y) = (ogeom_math::Direction::X, ogeom_math::Direction::Y);
        let (up, down) = (ogeom_math::Direction::Z, -ogeom_math::Direction::Z);
        let plain = circle(x, up);
        let turned = circle(y, up);
        let backwards = circle(x, down);
        // `a`'s share of its own length that `b` covers.
        for (a, b, share) in [
            (trim(&plain, 5.5, 7.0), trim(&plain, 0.2, 1.0), 0.517 / 1.5),
            (trim(&plain, 0.2, 1.0), trim(&plain, 5.5, 7.0), 0.517 / 0.8),
            (
                trim(&plain, 5.8, 6.8),
                trim(&backwards, 5.0, 6.2),
                0.4336 / 1.0,
            ),
            (trim(&plain, 5.0, 7.5), backwards.clone(), 1.0),
            (trim(&plain, 0.5, 6.0), trim(&turned, 4.0, 6.0), 1.218 / 5.5),
        ] {
            let found = intersect_curves(&a, &b, CurveCurveOptions::default(), T).unwrap();
            let (wa, wb) = (a.domain(), b.domain());
            let mut covered = 0.0;
            for o in &found.overlaps {
                for (t, w) in [
                    (o.on_a.0, wa),
                    (o.on_a.1, wa),
                    (o.on_b.0, wb),
                    (o.on_b.1, wb),
                ] {
                    assert!(t >= w.0 - 1e-9 && t <= w.1 + 1e-9, "{t} outside {w:?}");
                }
                for (ta, tb) in [(o.on_a.0, o.on_b.0), (o.on_a.1, o.on_b.1)] {
                    let gap = a
                        .point_at(ta, T)
                        .unwrap()
                        .distance(b.point_at(tb, T).unwrap());
                    assert!(gap < 1e-9, "ends {gap} apart");
                }
                covered += (o.on_a.1 - o.on_a.0).abs();
            }
            let want = share * (wa.1 - wa.0);
            assert!((covered - want).abs() < 2e-3, "{covered} against {want}");
        }
    }

    /// A trimmed arc answers as its circle does inside its window: a line
    /// tangent to it, one crossing it twice a hair below the tangent, and a
    /// circle touching it from outside are all found.
    #[test]
    fn a_trimmed_arc_meets_tangents_as_its_circle_does() {
        use ogeom_geom::Trimmed2d;
        let options = CurveCurveOptions::default();
        let circle = circle2(Point2::new(0.0, 0.0), 10.0);
        let arc: PlanarCurve = Trimmed2d::new(circle, 0.3, 2.9, T).unwrap().into();
        let tangent = line2(Point2::new(-20.0, 10.0), Point2::new(20.0, 10.0));
        let found = intersect_curves_2d(&arc, &tangent, options, T).unwrap();
        assert_eq!(found.crossings.len(), 1);
        assert!(found.crossings[0].point.is_equal(Point2::new(0.0, 10.0), T));
        let grazing = line2(
            Point2::new(-20.0, 10.0 - 1e-4),
            Point2::new(20.0, 10.0 - 1e-4),
        );
        let found = intersect_curves_2d(&arc, &grazing, options, T).unwrap();
        assert_eq!(found.crossings.len(), 2);
        let beside: PlanarCurve =
            Trimmed2d::new(circle2(Point2::new(20.0, 0.0), 10.0), 2.0, 4.5, T)
                .unwrap()
                .into();
        let whole = circle2(Point2::new(0.0, 0.0), 10.0);
        let found = intersect_curves_2d(&whole, &beside, options, T).unwrap();
        assert_eq!(found.crossings.len(), 1);
        assert!(found.crossings[0].point.is_equal(Point2::new(10.0, 0.0), T));
        // Outside the window, nothing.
        let below = line2(Point2::new(-20.0, -10.0), Point2::new(20.0, -10.0));
        assert!(
            intersect_curves_2d(&arc, &below, options, T)
                .unwrap()
                .is_empty()
        );
    }
}
