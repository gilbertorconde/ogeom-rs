//! Where two curves meet, in the plane and in space.
//!
//! Crossings in the plane are solved directly; crossings in space go through
//! extrema. The planar case is the load-bearing one: boolean face splitting
//! happens in a surface's parameter space, and the curves it splits with are
//! pcurves, so 2D curve/curve is the operation the whole boolean pipeline
//! stands on.
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

mod closed_form;
mod conics;
mod general;
mod helpers;

use closed_form::{circle_circle_2d, line_circle_2d, line_line_2d, line_line_3d};
use conics::{analytic_3d, same_curve_3d, tangent_line_conic_3d};
use general::{Sampled, Segments, foot_on_3d, general_2d, general_3d, sample_3d};
use helpers::{
    clamp_2d, clamp_3d, distance_to_curve_3d, one_crossing_radius, order, push_unique_2d,
    push_unique_3d, segments_approach_3d, segments_cross_2d, sort_crossings, valley_ends_3d,
    valley_span_3d,
};

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
    if let Some(found) = tangent_line_conic_3d(a, b, options, tol) {
        return Ok(found);
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

fn check(options: CurveCurveOptions) -> OgeomResult<()> {
    if options.samples < 2 {
        ogeom_bail!(Construction, "seeding needs at least two segments");
    }
    if !options.gap.is_finite() || options.gap <= 0.0 {
        ogeom_bail!(Construction, "a gap of {} is not a distance", options.gap);
    }
    Ok(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests;
