//! Where two surfaces meet: the one call.
//!
//! Everything else in this crate is a stage: closed forms, seeding, tracing,
//! fitting. This is the function an application calls, and the one `ogeom-bool`
//! builds on: give it two surfaces, get back what they do to each other,
//! with the analytic path taken where it exists and the marched-and-fitted
//! path where it does not. The caller does not choose; the pair does.
//!
//! # What a section curve carries
//!
//! Three descriptions, because three consumers: the curve in space for the
//! edge, and a pcurve per surface for the faces; face splitting happens in
//! parameter space, and a curve a face cannot express is one it cannot be
//! split along. Analytic results carry exact pcurves where the projection has
//! a closed form and `None` where it does not; fitted results always carry
//! fitted pcurves, because the tracer recorded the parameters as it walked.
//!
//! A pcurve here is **same-parameter** with its 3D curve: evaluating either at
//! the same `t` lands on the same point of the intersection. That is the claim
//! `docs/DATA_MODEL.md` §6 makes edges carry, and it is arranged here by
//! construction (the 2D curves inherit the 3D curve's own parameterization)
//! rather than asserted and repaired later.

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{
    Circle2d, Curve, Curve2d as _, Curve3d, Ellipse2d, Line2d, Offset2d, PlanarCurve, Surface,
    SurfaceGeometry,
};
use ogeom_math::{Circle2, Ellipse2, Frame2, Point, Point2};

use crate::approx::approximate_branch;
use crate::contact::trace_tangential;
use crate::march::{Marching, TANGENTIAL_SINE, branches_from_touches, leaves_at_an_angle};
use crate::surface::{Meeting, surface_surface};

/// How to intersect, when the general path runs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IntersectOptions {
    /// The tolerance the fitted curves are held to.
    pub tolerance: f64,
    /// The marching settings, for pairs with no closed form.
    pub marching: Marching,
}

impl Default for IntersectOptions {
    fn default() -> Self {
        Self {
            tolerance: 1e-6,
            marching: Marching::default(),
        }
    }
}

/// One curve of a section, with its parameter-space descriptions.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionCurve {
    /// The curve in space.
    pub curve: Curve,
    /// The curve in the first surface's parameter space, where it has one.
    ///
    /// Always present for a fitted curve. For an exact curve, present when the
    /// projection has a closed form (a line on a plane, a circle on the
    /// cylinder it wraps) and `None` where it does not, which is a statement
    /// about the projection rather than about the curve.
    pub on_a: Option<PlanarCurve>,
    /// The same, on the second surface.
    pub on_b: Option<PlanarCurve>,
    /// How far this curve may sit from the true intersection.
    ///
    /// Zero for an exact curve. For a fitted one, the trace's chord tolerance
    /// plus the fit's reported error: the sum of the stated parts.
    pub tolerance: f64,
    /// Whether the curve came from a closed form.
    pub exact: bool,
    /// Whether it is a closed loop.
    pub closed: bool,
    /// Whether the surfaces *touch* along this curve rather than crossing
    /// it.
    ///
    /// A tangential contact is a real curve (the two surfaces meet there,
    /// and a drawing has to show it), but it carries no boundary parity:
    /// neither surface passes through the other, so nothing is inside on
    /// one side and outside on the other. Consumers that classify by
    /// crossing must leave these out of that arithmetic; consumers that
    /// draw or measure contact want them.
    pub tangential: bool,
}

/// What two surfaces do to each other.
#[derive(Debug, Clone, PartialEq)]
pub enum SurfaceIntersection {
    /// They do not meet.
    ///
    /// From the general path this means *no crossing was found at the seeding
    /// resolution*: a branch thinner than the sampling grid is found only from
    /// a border it pierces or from the approach of the surfaces near it, and
    /// the completeness instrument in `tests/support/coverage.rs` is what
    /// checks.
    Apart,
    /// They touch at isolated points without crossing.
    Touching(Vec<Point>),
    /// They meet along these curves.
    Along(Vec<SectionCurve>),
    /// They are the same surface wherever they overlap.
    Same,
}

/// Where two surfaces meet.
///
/// The analytic path answers the pairs with closed forms, exactly, with
/// tolerance zero. Every other pair is seeded, traced and fitted to
/// `options.tolerance`, after a pair without a closed form is first measured
/// for coincidence: two patches lying on one another to within
/// `options.tolerance` wherever they overlap are [`Same`]. One call, and the
/// pair decides the path.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the options
/// are unusable. A pair the marcher finds nothing for is [`Apart`], not an
/// error; see that variant for what it can and cannot claim.
///
/// [`Apart`]: SurfaceIntersection::Apart
/// [`Same`]: SurfaceIntersection::Same
pub fn intersect_surfaces(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: IntersectOptions,
    tol: Tolerances,
) -> OgeomResult<SurfaceIntersection> {
    intersect_surfaces_within(a, b, (None, None), options, tol)
}

/// A rectangle of a surface's parameters, `((u0, u1), (v0, v1))`.
pub type ParameterWindow = ((f64, f64), (f64, f64));

/// [`intersect_surfaces`] with the marched path held to a window of each
/// surface's parameters.
///
/// Where the pair has a closed form the windows change nothing. Where it
/// is marched, seeding and tracing walk only each surface's window, so a
/// face trimmed from a larger patch is not seeded and traced over the
/// whole patch; a branch is followed to the window's edge and stops there.
/// The chart is untouched, so every pcurve returned means on the surface
/// what it would have meant without the window.
///
/// A window is clamped to its surface's domain in a direction that does
/// not repeat. In a periodic direction it is kept as given, so it may
/// straddle the seam and the walk reads it in the window's own period;
/// one that covers a period leaves that direction whole. A drum keeps its
/// turn whole that way and is narrowed along its axis alone; any other
/// surface with a periodic direction left whole is walked whole, since a
/// window does not wrap. `None` leaves that surface whole.
///
/// # Errors
///
/// As [`intersect_surfaces`].
pub fn intersect_surfaces_within(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    windows: (Option<ParameterWindow>, Option<ParameterWindow>),
    options: IntersectOptions,
    tol: Tolerances,
) -> OgeomResult<SurfaceIntersection> {
    if !options.tolerance.is_finite() || options.tolerance <= 0.0 {
        ogeom_bail!(
            Construction,
            "a tolerance of {} is not a distance",
            options.tolerance
        );
    }

    // A plane all but along a drum's axis meets it in an ellipse
    // kilometres long, whose parameter is too coarse a ruler for the few
    // millimetres of it the drum's height holds: a crossing solved on it
    // lands tens of microns off. Over that height it is two lines.
    if let Some(sections) = near_parallel_plane_drum(a, b, tol) {
        return Ok(if sections.is_empty() {
            SurfaceIntersection::Apart
        } else {
            SurfaceIntersection::Along(sections)
        });
    }
    match surface_surface(a, b, tol) {
        Ok(Meeting::Apart) => Ok(SurfaceIntersection::Apart),
        Ok(Meeting::Same) => Ok(SurfaceIntersection::Same),
        Ok(Meeting::Touching(points)) => Ok(SurfaceIntersection::Touching(points)),
        Ok(Meeting::Along(curves)) => {
            let sections: Vec<SectionCurve> = curves
                .into_iter()
                .filter_map(|curve| exact_section(curve, a, b, tol))
                .collect();
            Ok(if sections.is_empty() {
                // Every curve fell outside the surfaces' stated extents: the
                // unbounded geometries meet, the surfaces as given do not.
                SurfaceIntersection::Apart
            } else {
                SurfaceIntersection::Along(sections)
            })
        }
        // No closed form for this pair. Coincidence is measured first: the
        // marcher seeds on sign changes, and a pair that never separates
        // has none, so what it would trace over a coincident pair is noise
        // under a section's name, and costs seconds to produce.
        Err(_) if surfaces_coincide(a, b, options.tolerance, tol) => Ok(SurfaceIntersection::Same),
        // Otherwise the marcher, unless the pair is two drums all but
        // parallel.
        Err(_) => match near_parallel_drums(a, b, tol)
            .or_else(|| ball_through_drum(a, b, tol))
            .or_else(|| axial_plane_revolution(a, b, tol))
            .or_else(|| crate::extrusion::plane_extrusion(a, b, tol))
            .or_else(|| plane_along_spline_lines(a, b, tol))
        {
            Some(sections) if sections.is_empty() => Ok(SurfaceIntersection::Apart),
            Some(sections) => Ok(SurfaceIntersection::Along(sections)),
            None => marched(a, b, windows, options, tol),
        },
    }
}

/// Whether two surfaces are one surface wherever they overlap, measured to
/// `reach`.
///
/// The closed forms answer this for the pairs they know. Where there is no
/// closed form it does not stop being a fair question (two patches restated
/// from one plane are the same surface, and nothing in their control points
/// says so), but it stops being answerable exactly, so it is measured, and
/// only for the pairs the analytic layer has declined.
///
/// Sampled on the *smaller* window, because the answer is about the region
/// the two share and a stated window is not that region: a plane's own
/// extends for a billion units either way, and a grid over it samples
/// nothing. A sample whose foot lands on the rim of the other is skipped
/// rather than counted against: the other patch simply does not reach that
/// far, and a distance measured to its rim is about the window, not the
/// surface.
///
/// The samples that agree must spread over an area, not lie along one
/// line of either chart: samples along one curve prove only that the curve
/// lies in both surfaces, as a wall's border does in the plane it leaves.
/// Where they do lie along one line, each is stepped off it in the sampled
/// chart, by half a grid step and then by halves of that until its foot
/// lands inside the other's window; one such sample must agree, and none
/// may stand off.
///
/// One-sided by construction: a pair that crosses puts interior samples well
/// off the other, so it cannot pass, and a pair this cannot resolve is
/// marched.
fn surfaces_coincide(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    reach: f64,
    tol: Tolerances,
) -> bool {
    /// Samples per direction over the window, and how many must land inside
    /// the other before agreement means anything.
    const GRID: i32 = 6;
    const EVIDENCE: usize = 4;
    /// How many times a step off a line of samples is halved looking for
    /// the other's window: down to a millionth of a grid step.
    const HALVINGS: u32 = 20;

    let span = |s: &SurfaceGeometry| -> f64 {
        let ((ua, ub), (va, vb)) = s.domain();
        (ub - ua).abs().max((vb - va).abs())
    };
    let (sampled, against) = if span(a) <= span(b) { (a, b) } else { (b, a) };
    let ((ua, ub), (va, vb)) = sampled.domain();
    if !(ua.is_finite() && ub.is_finite() && va.is_finite() && vb.is_finite()) {
        return false;
    }
    let ((wu0, wu1), (wv0, wv1)) = against.domain();
    // A twentieth of the window in from each rim: enough that a foot the
    // search pinned to the rim is not read as one the surface truly reaches.
    let (mu, mv) = ((wu1 - wu0) * 0.05, (wv1 - wv0) * 0.05);
    let (du, dv) = ((ub - ua) / f64::from(GRID), (vb - va) / f64::from(GRID));

    // The foot on the other of the sample at `(u, v)`: `Some(None)` where it
    // lands on the window's rim, `None` where the sample stands off the
    // other or cannot be measured.
    let measure = |u: f64, v: f64| -> Option<Option<(f64, f64)>> {
        let p = sampled.point_at(u, v, tol).ok()?;
        let foot = ogeom_geom::project_on_surface(against, p, 16, tol).ok()?;
        let (fu, fv) = foot.parameters;
        if fu <= wu0 + mu || fu >= wu1 - mu || fv <= wv0 + mv || fv >= wv1 - mv {
            return Some(None);
        }
        (foot.distance <= reach).then_some(Some(foot.parameters))
    };
    let mut agreeing: Vec<((i32, i32), (f64, f64))> = Vec::new();
    for i in 0..=GRID {
        for j in 0..=GRID {
            let (u, v) = (du.mul_add(f64::from(i), ua), dv.mul_add(f64::from(j), va));
            match measure(u, v) {
                None => return false,
                Some(None) => {}
                Some(Some(foot)) => agreeing.push(((i, j), foot)),
            }
        }
    }
    if agreeing.len() < EVIDENCE {
        return false;
    }

    // Along one line of the sampled grid, exactly; and of the other's
    // chart, to the parametric confusion of the feet's spread.
    let (first, _) = agreeing[0];
    let lined_up = agreeing
        .iter()
        .find(|(at, _)| *at != first)
        .map(|&(at, _)| (at.0 - first.0, at.1 - first.1))
        .filter(|&(di, dj)| {
            agreeing
                .iter()
                .all(|&((i, j), _)| di * (j - first.1) == dj * (i - first.0))
        });
    let feet_lined_up = {
        let mut far = (agreeing[0].1, agreeing[0].1, 0.0_f64);
        for (x, (_, p)) in agreeing.iter().enumerate() {
            for (_, q) in &agreeing[x + 1..] {
                let d = (q.0 - p.0).hypot(q.1 - p.1);
                if d > far.2 {
                    far = (*p, *q, d);
                }
            }
        }
        let (p, q, d) = far;
        agreeing.iter().all(|(_, r)| {
            ((q.0 - p.0) * (r.1 - p.1) - (q.1 - p.1) * (r.0 - p.0)).abs()
                <= tol.parametric() * d * d
        })
    };
    if lined_up.is_none() && !feet_lined_up {
        return true;
    }

    // Off the line: across it in the sampled chart, or either way where only
    // the feet line up.
    let across: &[(f64, f64)] = match lined_up {
        Some((0, _)) => &[(1.0, 0.0)],
        Some(_) => &[(0.0, 1.0)],
        None => &[(1.0, 0.0), (0.0, 1.0)],
    };
    let mut spread = false;
    for &((i, j), _) in &agreeing {
        for &(ku, kv) in across {
            for sign in [-1.0, 1.0] {
                let mut h = 0.5 * sign;
                for _ in 0..HALVINGS {
                    let (u, v) = (
                        du.mul_add(ku.mul_add(h, f64::from(i)), ua),
                        dv.mul_add(kv.mul_add(h, f64::from(j)), va),
                    );
                    if !(ua..=ub).contains(&u) || !(va..=vb).contains(&v) {
                        break;
                    }
                    match measure(u, v) {
                        None => return false,
                        Some(None) => h *= 0.5,
                        Some(Some(_)) => {
                            spread = true;
                            break;
                        }
                    }
                }
            }
        }
    }
    spread
}

/// A plane through the axis of a surface of revolution whose profile lies
/// in a plane through that axis: the profile turned to each angle that sets
/// its plane on the cut, a column of the revolution's chart.
///
/// Marched, such a section runs through the pole wherever the profile meets
/// the axis, where the chart pinches to a point and the march stalls. Here
/// the profile is cut where it crosses the axis, and each piece is turned to
/// both angles, the one setting its side of the axis on each half of the
/// cut: a revolution whose profile crosses the axis covers the section
/// twice, once from each column, and a face on either column finds its own
/// piece. `None` for any other plane or profile, and where the turns fall
/// outside the sweep the answer is empty.
fn axial_plane_revolution(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<Vec<SectionCurve>> {
    const SAMPLES: u32 = 64;
    let (plane, revolution, plane_first) = match (a, b) {
        (SurfaceGeometry::Plane(p), SurfaceGeometry::Revolution(r)) => (p, r, true),
        (SurfaceGeometry::Revolution(r), SurfaceGeometry::Plane(p)) => (p, r, false),
        _ => return None,
    };
    let cut = plane.plane();
    let axis = revolution.axis();
    let normal = cut.normal();
    if normal.dot(axis.direction).abs() > tol.angular()
        || cut.signed_distance_to(axis.location).abs() > tol.confusion()
    {
        return None;
    }
    let profile = revolution.curve();
    let (v0, v1) = profile.domain();
    let at = |k: u32| v0 + (v1 - v0) * f64::from(k) / f64::from(SAMPLES);
    let radial = |v: f64| {
        let p = profile.point_at(v, tol).ok()?;
        Some(p - axis.project(p))
    };
    // The profile's own side of the axis, from its point furthest off it.
    let mut widest = ogeom_math::Vector::ZERO;
    for k in 0..=SAMPLES {
        let r = radial(at(k))?;
        if r.magnitude() > widest.magnitude() {
            widest = r;
        }
    }
    let side = ogeom_math::Direction::new(widest, tol).ok()?;
    let across = axis.direction.cross_with(side.vector());
    // Every point of the profile in the plane of the axis and that side.
    let offset = |v: f64| radial(v).map(|r| (r.dot(side.vector()), r.dot(across)));
    for k in 0..=SAMPLES {
        let (_, off) = offset(at(k))?;
        if off.abs() > tol.confusion() {
            return None;
        }
    }
    // The pieces between the profile's crossings of the axis, each crossing
    // narrowed by bisection.
    let mut cuts = vec![v0];
    for k in 0..SAMPLES {
        let (mut lo, mut hi) = (at(k), at(k + 1));
        let (s_lo, s_hi) = (offset(lo)?.0, offset(hi)?.0);
        if s_lo.abs() <= tol.confusion() || s_lo * s_hi >= 0.0 {
            continue;
        }
        for _ in 0..80 {
            let mid = f64::midpoint(lo, hi);
            if offset(mid)?.0 * s_lo > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        cuts.push(f64::midpoint(lo, hi));
    }
    cuts.push(v1);
    // The turns setting the profile's side on the cut's two halves.
    let out = axis.direction.cross_with(normal.vector());
    let first = across.dot(out).atan2(side.vector().dot(out));
    let (u0, u1) = revolution.domain().0;
    let turns: Vec<f64> = [first, first + core::f64::consts::PI]
        .into_iter()
        .filter_map(|u| {
            let u = u0 + (u - u0).rem_euclid(core::f64::consts::TAU);
            let u = if (u - u0 - core::f64::consts::TAU).abs() <= tol.angular() {
                u0
            } else {
                u
            };
            (u <= u1 + tol.angular()).then_some(u.min(u1))
        })
        .collect();
    let mut sections = Vec::new();
    for &u in &turns {
        let turned = ogeom_geom::Transformable::transformed(
            profile,
            &ogeom_math::Transform::rotation(axis, u),
            tol,
        )
        .ok()?;
        for piece in cuts.windows(2) {
            let (va, vb) = (piece[0], piece[1]);
            if vb - va <= tol.parametric() {
                continue;
            }
            let curve: Curve = ogeom_geom::TrimmedCurve::new(turned.clone(), va, vb, tol)
                .ok()?
                .into();
            let column: PlanarCurve = Line2d::over(
                ogeom_math::Axis2::new(Point2::new(u, 0.0), ogeom_math::Direction2::Y),
                va,
                vb,
            )
            .ok()?
            .into();
            let flat = exact_pcurve(&curve, (va, vb), a_or_b(plane_first, a, b), tol);
            let (on_a, on_b) = if plane_first {
                (flat, Some(column))
            } else {
                (Some(column), flat)
            };
            sections.push(SectionCurve {
                on_a,
                on_b,
                tolerance: 0.0,
                exact: true,
                closed: false,
                tangential: false,
                curve,
            });
        }
    }
    Some(sections)
}

/// A plane holding whole columns or rows of a spline surface, and meeting
/// it nowhere else: those iso lines, exactly.
///
/// A surface of revolution converted to a spline and scaled keeps its
/// meridians as columns, and a plane through its axis holds two of them,
/// one often the seam along the chart's border. Marched, such a section
/// runs along the chart's edge or through its poles and is not found. An
/// iso line lies in the plane where every control point of it does, which
/// is a root of each control point's weighted distance to the plane, found
/// along the chart. The answer stands only where a grid over the chart
/// finds the surface on one side of the plane between the lines found;
/// anything else meets the plane elsewhere too and is marched.
///
/// A line on an open border of the chart has the surface on one side only.
/// It stands where the surface leaves the plane there at an angle all
/// along it, as [`leaves_at_an_angle`] has it, against the grid's nearest
/// row in and the band its samples ignore: the plane then meets the face
/// along its boundary edge, as a plane through a prism's side edge meets
/// the side faces. A border the surface leaves tangentially is a contact,
/// and the pair is marched.
///
/// An extrusion of a spline profile is read as the spline patch it is, on
/// the same chart: its columns are its straight rulings, and are stated as
/// lines.
fn plane_along_spline_lines(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<Vec<SectionCurve>> {
    const SAMPLES: u32 = 96;
    let (plane, other, plane_first) = match (a, b) {
        (SurfaceGeometry::Plane(p), other) => (p.plane(), other, true),
        (other, SurfaceGeometry::Plane(p)) => (p.plane(), other, false),
        _ => return None,
    };
    let extruded;
    let (spline, ruled) = match other {
        SurfaceGeometry::BSpline(s) => (s, None),
        SurfaceGeometry::Extrusion(e) => {
            extruded = extrusion_patch(e, tol)?;
            (&extruded, Some(e))
        }
        _ => return None,
    };
    let distance = |p: Point| plane.signed_distance_to(p);
    let ((u0, u1), (v0, v1)) = spline.domain();
    // An iso line's control points' weighted distances to the plane, and
    // whether the line has any length.
    let line_at = |along_u: bool, t: f64| -> Option<ogeom_geom::BSplineCurve> {
        if along_u {
            spline.iso_u_curve(t, tol).ok()
        } else {
            spline.iso_v_curve(t, tol).ok()
        }
    };
    let weighted = |curve: &ogeom_geom::BSplineCurve| -> Vec<f64> {
        curve
            .control_points()
            .iter()
            .map(|w| w.weight * distance(w.point()))
            .collect()
    };
    let lies_in = |curve: &ogeom_geom::BSplineCurve| {
        curve
            .control_points()
            .iter()
            .all(|w| distance(w.point()).abs() <= tol.confusion())
    };
    let has_length = |curve: &ogeom_geom::BSplineCurve| {
        let first = curve.control_points()[0].point();
        curve
            .control_points()
            .iter()
            .any(|w| w.point().distance(first) > tol.confusion())
    };
    // The iso lines of one family lying in the plane: the chart's borders,
    // and every root of the control point that strays furthest.
    let found = |along_u: bool| -> Option<Vec<f64>> {
        let (lo, hi) = if along_u { (u0, u1) } else { (v0, v1) };
        let at = |k: u32| lo + (hi - lo) * f64::from(k) / f64::from(SAMPLES);
        let rows: Vec<Vec<f64>> = (0..=SAMPLES)
            .map(|k| line_at(along_u, at(k)).map(|c| weighted(&c)))
            .collect::<Option<_>>()?;
        let count = rows[0].len();
        if rows.iter().any(|r| r.len() != count) {
            return None;
        }
        let widest = (0..count).max_by(|&i, &j| {
            let spread = |i: usize| rows.iter().fold(0.0_f64, |m, r| m.max(r[i].abs()));
            spread(i).total_cmp(&spread(j))
        })?;
        let mut roots = vec![lo, hi];
        for k in 0..SAMPLES {
            let (mut a, mut b) = (at(k), at(k + 1));
            let (da, db) = (rows[k as usize][widest], rows[k as usize + 1][widest]);
            if da == 0.0 {
                roots.push(a);
                continue;
            }
            if da * db > 0.0 {
                continue;
            }
            let sign = da.signum();
            for _ in 0..80 {
                let mid = f64::midpoint(a, b);
                let d = weighted(&line_at(along_u, mid)?)[widest];
                if d * sign > 0.0 {
                    a = mid;
                } else {
                    b = mid;
                }
            }
            roots.push(f64::midpoint(a, b));
        }
        roots.sort_by(f64::total_cmp);
        roots.dedup_by(|x, y| (*x - *y).abs() <= tol.parametric());
        Some(
            roots
                .into_iter()
                .filter(|&t| line_at(along_u, t).is_some_and(|c| lies_in(&c) && has_length(&c)))
                .collect(),
        )
    };
    let columns = found(true)?;
    let rows = found(false)?;
    if columns.is_empty() && rows.is_empty() {
        return None;
    }
    // The chart cut by the lines found into cells: the surface keeps to one
    // side of the plane within each, or it meets the plane elsewhere too.
    let strip = |lines: &[f64], t: f64| lines.iter().filter(|&&x| x < t).count();
    let near_line = |lines: &[f64], t: f64, span: f64| {
        lines
            .iter()
            .any(|&x| (x - t).abs() <= span / f64::from(SAMPLES) * 0.25)
    };
    let mut sides: ogeom_core::FastMap<(usize, usize), f64> = ogeom_core::FastMap::default();
    let band = tol.confusion() * 10.0;
    for i in 0..=SAMPLES {
        let u = u0 + (u1 - u0) * (f64::from(i) + 0.5) / f64::from(SAMPLES + 1);
        if near_line(&columns, u, u1 - u0) {
            continue;
        }
        for j in 0..=SAMPLES {
            let v = v0 + (v1 - v0) * (f64::from(j) + 0.5) / f64::from(SAMPLES + 1);
            if near_line(&rows, v, v1 - v0) {
                continue;
            }
            let d = distance(spline.point_at(u, v, tol).ok()?);
            if d.abs() <= band {
                continue;
            }
            let cell = (strip(&columns, u), strip(&rows, v));
            match sides.get(&cell) {
                Some(side) if side * d < 0.0 => return None,
                Some(_) => {}
                None => {
                    sides.insert(cell, d.signum());
                }
            }
        }
    }
    // Each line must be a crossing, the cells either side of it on
    // opposite sides of the plane (across the border of a closed chart,
    // the cells at its two ends), or lie on an open border the surface
    // leaves at an angle, the one cell beside it on a side. A line the
    // surface only touches, or with no cell beside it to say, is a tangency
    // the marcher and the contact handling answer, and so is the whole pair.
    let closed_u = spline.is_closed_u(tol);
    let closed_v = spline.is_closed_v(tol);
    let plane_normal = plane.normal().vector();
    let leaves_plane = |along_u: bool, t: f64| -> bool {
        let (lo, hi) = if along_u { (v0, v1) } else { (u0, u1) };
        // The grid's nearest row in from the border, half a step in.
        let (from, to) = if along_u { (u0, u1) } else { (v0, v1) };
        let half_step = (to - from) / (2.0 * f64::from(SAMPLES + 1));
        let inside = if (t - from).abs() <= (to - t).abs() {
            t + half_step
        } else {
            t - half_step
        };
        [0.0, 0.07, 0.19, 0.37, 0.53, 0.71, 0.89, 1.0]
            .iter()
            .all(|&f| {
                let s = (hi - lo).mul_add(f, lo);
                let (on, near) = if along_u {
                    ((t, s), (inside, s))
                } else {
                    ((s, t), (s, inside))
                };
                let (Ok(n), Ok(p), Ok(q)) = (
                    spline.normal_at(on.0, on.1, tol),
                    spline.point_at(on.0, on.1, tol),
                    spline.point_at(near.0, near.1, tol),
                ) else {
                    return false;
                };
                let sine = n.vector().cross(plane_normal).magnitude();
                leaves_at_an_angle(sine, p.distance(q), distance(q).abs(), band)
            })
    };
    let on_open_border = |t: f64, (lo, hi): (f64, f64), closed: bool| {
        !closed && ((t - lo).abs() <= tol.parametric() || (t - hi).abs() <= tol.parametric())
    };
    let side_of = |cu: Option<usize>, cv: Option<usize>| -> Option<f64> {
        let mut found = sides
            .iter()
            .filter(|((u, v), _)| cu.is_none_or(|c| c == *u) && cv.is_none_or(|c| c == *v))
            .map(|(_, s)| *s);
        let first = found.next()?;
        found.all(|s| s == first).then_some(first)
    };
    let crosses = |k: usize,
                   count: usize,
                   closed: bool,
                   border: bool,
                   cell: &dyn Fn(usize) -> Option<f64>| {
        let below = cell(k).or_else(|| closed.then(|| (0..=count).rev().find_map(cell)).flatten());
        let above = cell(k + 1).or_else(|| closed.then(|| (0..=count).find_map(cell)).flatten());
        match (below, above) {
            (Some(x), Some(y)) => x * y < 0.0,
            (Some(_), None) | (None, Some(_)) => border,
            (None, None) => false,
        }
    };
    for (k, &u) in columns.iter().enumerate() {
        let border = on_open_border(u, (u0, u1), closed_u) && leaves_plane(true, u);
        if !crosses(k, columns.len(), closed_u, border, &|c| {
            side_of(Some(c), None)
        }) {
            return None;
        }
    }
    for (k, &v) in rows.iter().enumerate() {
        let border = on_open_border(v, (v0, v1), closed_v) && leaves_plane(false, v);
        if !crosses(k, rows.len(), closed_v, border, &|c| side_of(None, Some(c))) {
            return None;
        }
    }
    let mut sections = Vec::new();
    let mut emit = |along_u: bool, t: f64| -> Option<()> {
        let iso = line_at(along_u, t)?;
        let curve: Curve = match ruled {
            Some(extrusion) if along_u => {
                let range = Curve3d::domain(&iso);
                let axis = ogeom_math::Axis {
                    location: extrusion.curve().point_at(t, tol).ok()?,
                    direction: extrusion.direction(),
                };
                ogeom_geom::LineCurve::over(axis, range.0, range.1)
                    .ok()?
                    .into()
            }
            _ => iso.into(),
        };
        let range = curve.domain();
        let chart: PlanarCurve = if along_u {
            Line2d::over(
                ogeom_math::Axis2::new(Point2::new(t, 0.0), ogeom_math::Direction2::Y),
                range.0,
                range.1,
            )
        } else {
            Line2d::over(
                ogeom_math::Axis2::new(Point2::new(0.0, t), ogeom_math::Direction2::X),
                range.0,
                range.1,
            )
        }
        .ok()?
        .into();
        let flat = exact_pcurve(&curve, range, a_or_b(plane_first, a, b), tol)?;
        let (on_a, on_b) = if plane_first {
            (Some(flat), Some(chart))
        } else {
            (Some(chart), Some(flat))
        };
        sections.push(SectionCurve {
            on_a,
            on_b,
            tolerance: 0.0,
            exact: true,
            closed: curve.is_closed(tol),
            tangential: false,
            curve,
        });
        Some(())
    };
    // Every line found: the chart's two borders across a closed surface
    // are one line in space, stated once.
    for (k, &u) in columns.iter().enumerate() {
        if closed_u
            && k + 1 == columns.len()
            && k > 0
            && columns[0] == u0
            && (u - u1).abs() <= tol.parametric()
        {
            continue;
        }
        emit(true, u)?;
    }
    for (k, &v) in rows.iter().enumerate() {
        if closed_v
            && k + 1 == rows.len()
            && k > 0
            && rows[0] == v0
            && (v - v1).abs() <= tol.parametric()
        {
            continue;
        }
        emit(false, v)?;
    }
    Some(sections)
}

/// An extrusion of a spline profile as a spline patch on the extrusion's
/// own chart: the profile's knots along `u`, and a straight run along `v`
/// over the extrusion's window, so each point and each parameter is the
/// extrusion's. `None` for any other profile.
fn extrusion_patch(
    extrusion: &ogeom_geom::ExtrusionSurface,
    tol: Tolerances,
) -> Option<ogeom_geom::BSplineSurface> {
    let Curve::BSpline(profile) = extrusion.curve() else {
        return None;
    };
    let (_, (v0, v1)) = extrusion.domain();
    let step = extrusion.direction().vector();
    let mut grid = Vec::with_capacity(profile.control_points().len() * 2);
    for c in profile.control_points() {
        for v in [v0, v1] {
            grid.push(ogeom_math::Weighted {
                scaled: Point::from_vector(c.scaled.to_vector() + step * (v * c.weight)),
                weight: c.weight,
            });
        }
    }
    let count = profile.control_points().len();
    ogeom_geom::BSplineSurface::rational(
        profile.knots().clone(),
        ogeom_math::KnotVector::new(vec![v0, v0, v1, v1], 1).ok()?,
        ogeom_math::ControlGrid::new(grid, count, 2).ok()?,
    )
    .ok()
    .filter(|patch| {
        let same = |p: (f64, f64), q: (f64, f64)| {
            (p.0 - q.0).abs() <= tol.parametric() && (p.1 - q.1).abs() <= tol.parametric()
        };
        let (pu, pv) = patch.domain();
        let (eu, ev) = extrusion.domain();
        same(pu, eu) && same(pv, ev)
    })
}

/// The first surface where `first` holds, else the second.
fn a_or_b<'s>(first: bool, a: &'s SurfaceGeometry, b: &'s SurfaceGeometry) -> &'s SurfaceGeometry {
    if first { a } else { b }
}

/// Two drums whose axes are all but parallel, over the height they share.
///
/// Parallel drums meet in straight lines along their axes, and drums whose
/// axes lean a ten-thousandth apart (a drilled hole beside a fillet of a
/// converted mesh, each axis fitted to its own facets) meet in a quartic
/// that departs from those lines by less than a micron over any height a
/// part has. Marched, it comes back as fitted curves that cost seconds to
/// cross and wander where the drums nearly touch. Here each is solved in
/// the cross-sections along the shared height and kept as the line through
/// its ends where every station lies near it, that departure stated as the
/// section's tolerance.
///
/// `None` where the axes lean further, where the drums do not cross
/// cleanly at every station (a crossing starting part way up, or a near
/// touch), or where a station strays: the marcher answers those. An empty
/// answer is drums that share no height.
fn near_parallel_drums(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<Vec<SectionCurve>> {
    const LEAN: f64 = 1e-3;
    let (SurfaceGeometry::Cylinder(sa), SurfaceGeometry::Cylinder(sb)) = (a, b) else {
        return None;
    };
    let (ca, cb) = (sa.cylinder(), sb.cylinder());
    let (axis_a, axis_b) = (ca.axis(), cb.axis());
    let (da, db) = (axis_a.direction.vector(), axis_b.direction.vector());
    let (ra, rb) = (ca.radius(), cb.radius());
    let cos = da.dot(db);
    if da.cross(db).magnitude() > LEAN || cos.abs() < 0.5 {
        return None;
    }
    let (pa, pb) = (axis_a.location, axis_b.location);
    // The shared height, measured along the first axis.
    let (_, (a0, a1)) = a.domain();
    let (_, (b0, b1)) = b.domain();
    let along = |v: f64| (pb - pa).dot(da) + v * cos;
    let (lo, hi) = (
        a0.min(a1).max(along(b0).min(along(b1))),
        a0.max(a1).min(along(b0).max(along(b1))),
    );
    if !(lo.is_finite() && hi.is_finite()) {
        return None;
    }
    if hi - lo <= tol.confusion() {
        return Some(Vec::new());
    }
    // Where the two cross-sections at a station meet, left and right of
    // the line of centres: the second drum's section is an ellipse only a
    // square of its lean away from a circle, a stated part of the stray.
    let meet = |z: f64| -> Option<[Point; 2]> {
        let centre_a = pa + da * z;
        let s = (centre_a - pb).dot(da) / cos;
        let centre_b = pb + db * s;
        let mut between = centre_b - centre_a;
        between = between - da * between.dot(da);
        let d = between.magnitude();
        // Axes, or a touch, well inside the weld distance are one: the
        // sliver between the drums is welded rather than sectioned. Half
        // of it, so a sliver at the weld distance itself is sectioned
        // whole rather than lost between the two readings.
        let margin = tol.confusion() * 50.0;
        if d <= margin || d >= ra + rb - margin || d <= (ra - rb).abs() + margin {
            return None;
        }
        let x = (d * d + ra * ra - rb * rb) / (2.0 * d);
        let h = (ra * ra - x * x).max(0.0).sqrt();
        let ex = between / d;
        let ey = da.cross(ex);
        Some([centre_a + ex * x + ey * h, centre_a + ex * x - ey * h])
    };
    lines_through_stations(lo, hi, meet, rb * (1.0 / cos.abs() - 1.0), tol)
}

/// How far a near-parallel pair's sections may stray from the true
/// crossing: what a fitted section typically carries.
const NEAR_PARALLEL_STRAY: f64 = 1e-5;

/// The two curves a near-parallel pair meets in over the height `lo..hi`,
/// from where `meet` puts the crossing at each height: the line through the
/// ends where every station lies within a micron of it, else a cubic
/// through the stations at their heights, checked midway between them.
/// Either is kept within [`NEAR_PARALLEL_STRAY`], the departure stated as
/// its tolerance. `None` where a station has no clean crossing or the
/// curve strays.
fn lines_through_stations(
    lo: f64,
    hi: f64,
    meet: impl Fn(f64) -> Option<[Point; 2]>,
    stated: f64,
    tol: Tolerances,
) -> Option<Vec<SectionCurve>> {
    const STATIONS: u32 = 32;
    const STRAIGHT: f64 = 1e-6;
    let at = |k: f64| (hi - lo).mul_add(k / f64::from(STATIONS), lo);
    let heights: Vec<f64> = (0..=STATIONS).map(|k| at(f64::from(k))).collect();
    let met: Vec<[Point; 2]> = heights.iter().map(|&z| meet(z)).collect::<Option<_>>()?;
    let between: Vec<[Point; 2]> = (0..STATIONS)
        .map(|k| meet(at(f64::from(k) + 0.5)))
        .collect::<Option<_>>()?;
    let mut out = Vec::with_capacity(2);
    for side in 0..2 {
        let (from, to) = (met[0][side], met[met.len() - 1][side]);
        let span = to - from;
        let length = span.magnitude();
        if length <= tol.confusion() {
            return None;
        }
        let off_line = |p: Point| {
            let t = (p - from).dot(span) / (length * length);
            p.distance(from + span * t)
        };
        let stray = met
            .iter()
            .chain(&between)
            .map(|pair| off_line(pair[side]))
            .fold(0.0_f64, f64::max);
        let (curve, stray): (Curve, f64) = if stray <= STRAIGHT {
            (
                ogeom_geom::LineCurve::segment(from, to, tol).ok()?.into(),
                stray,
            )
        } else {
            let points: Vec<Point> = met.iter().map(|pair| pair[side]).collect();
            let fitted =
                ogeom_geom::fit::fit_points_at(&heights, &points, 3, tol.confusion(), tol).ok()?;
            let curve: Curve = fitted.curve.into();
            let mut worst = fitted.error;
            for (k, pair) in (0..STATIONS).zip(&between) {
                let p = curve.point_at(at(f64::from(k) + 0.5), tol).ok()?;
                worst = worst.max(p.distance(pair[side]));
            }
            (curve, worst)
        };
        let tolerance = stray + stated + tol.confusion();
        if tolerance > NEAR_PARALLEL_STRAY {
            return None;
        }
        out.push(SectionCurve {
            curve,
            on_a: None,
            on_b: None,
            tolerance,
            exact: false,
            closed: false,
            tangential: false,
        });
    }
    Some(out)
}

/// A drum passing clean through a ball: every line along the drum meets
/// the ball twice, within the drum's height.
///
/// Then each of the two loops the drum and ball meet in is a function of
/// the angle round the drum: at each angle, where the line along the drum
/// enters and leaves the ball is a quadratic's two roots. The loops are
/// sampled so, exactly, and fitted closed, the fit's error stated as the
/// section's tolerance. Marched instead, a drum that all but grazes the
/// ball's far side leaves loops long and thin, and the trace wanders along
/// them past any bound. `None` where some line misses or grazes the ball,
/// or leaves the drum's height: the marcher answers those.
fn ball_through_drum(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<Vec<SectionCurve>> {
    const SAMPLES: u32 = 256;
    const STRAY: f64 = 1e-5;
    let (ball, drum, ball_first) = match (a, b) {
        (SurfaceGeometry::Sphere(s), SurfaceGeometry::Cylinder(c)) => (s, c, true),
        (SurfaceGeometry::Cylinder(c), SurfaceGeometry::Sphere(s)) => (s, c, false),
        _ => return None,
    };
    let (sphere, cylinder) = (ball.sphere(), drum.cylinder());
    let frame = cylinder.frame();
    let (x, y, d) = (frame.x().vector(), frame.y().vector(), frame.z().vector());
    let (origin, r) = (frame.origin(), cylinder.radius());
    let (centre, big) = (sphere.centre(), sphere.radius());
    let ball_frame = sphere.frame();
    let (_, (h0, h1)) = drum.domain();
    // A line that only just meets the ball leaves the loop turning sharply
    // there; a tenth of the drum's radius of chord inside the ball keeps
    // the loops smooth enough to fit.
    let margin = r * 0.1;
    // Where the line along the drum at `angle` enters and leaves the ball.
    let heights = |angle: f64| -> Option<[f64; 2]> {
        let foot = origin + (x * angle.cos() + y * angle.sin()) * r;
        let w = foot - centre;
        let half = d.dot(w);
        let disc = half.mul_add(half, -(w.dot(w) - big * big));
        if disc <= margin * margin {
            return None;
        }
        let root = disc.sqrt();
        let pair = [-half - root, -half + root];
        pair.iter().all(|v| *v >= h0 && *v <= h1).then_some(pair)
    };
    let at = |angle: f64, v: f64| origin + (x * angle.cos() + y * angle.sin()) * r + d * v;
    // The ball's longitude and latitude of a point, as its chart reads them.
    let on_ball = |p: Point, before: Option<Point2>| -> Point2 {
        let local = ball_frame.to_local(p);
        let lat = local.z.atan2(local.x.hypot(local.y));
        let mut lon = local.y.atan2(local.x).rem_euclid(core::f64::consts::TAU);
        if let Some(prev) = before {
            while lon - prev.x > core::f64::consts::PI {
                lon -= core::f64::consts::TAU;
            }
            while prev.x - lon > core::f64::consts::PI {
                lon += core::f64::consts::TAU;
            }
        }
        Point2::new(lon, lat)
    };
    let angle_of = |k: f64| core::f64::consts::TAU * k / f64::from(SAMPLES);
    let params: Vec<f64> = (0..=SAMPLES).map(|k| angle_of(f64::from(k))).collect();
    let mut sampled: Vec<[f64; 2]> = Vec::with_capacity(params.len());
    for &angle in &params {
        sampled.push(heights(angle)?);
    }
    let mut out = Vec::with_capacity(2);
    for side in 0..2 {
        let points: Vec<Point> = params
            .iter()
            .zip(&sampled)
            .map(|(&angle, pair)| at(angle, pair[side]))
            .collect();
        let on_drum: Vec<Point2> = params
            .iter()
            .zip(&sampled)
            .map(|(&angle, pair)| Point2::new(angle, pair[side]))
            .collect();
        let mut on_sphere: Vec<Point2> = Vec::with_capacity(points.len());
        for p in &points {
            let q = on_ball(*p, on_sphere.last().copied());
            on_sphere.push(q);
        }
        let target = tol.confusion() * 10.0;
        let curve: Curve = ogeom_geom::fit::fit_points_at(&params, &points, 3, target, tol)
            .ok()?
            .curve
            .into();
        let drum_image: PlanarCurve =
            ogeom_geom::fit::fit_points_2d_at(&params, &on_drum, 3, target, tol)
                .ok()?
                .curve
                .into();
        let ball_image: PlanarCurve =
            ogeom_geom::fit::fit_points_2d_at(&params, &on_sphere, 3, target, tol)
                .ok()?
                .curve
                .into();
        // Checked at the samples and midway between them: the curve, and
        // each surface read through its image, against the true meeting.
        let mut stray = 0.0_f64;
        for k in 0..(2 * SAMPLES) {
            let angle = angle_of(f64::from(k) / 2.0);
            let truth = at(angle, heights(angle)?[side]);
            let on_curve = curve.point_at(angle, tol).ok()?;
            let uv = drum_image.point_at(angle, tol).ok()?;
            let through_drum = drum.point_at(uv.x, uv.y, tol).ok()?;
            let uv = ball_image.point_at(angle, tol).ok()?;
            let through_ball = ball.point_at(uv.x, uv.y, tol).ok()?;
            stray = stray
                .max(truth.distance(on_curve))
                .max(truth.distance(through_drum))
                .max(truth.distance(through_ball));
        }
        let tolerance = stray.max(tol.confusion());
        if tolerance > STRAY {
            return None;
        }
        let (on_a, on_b) = if ball_first {
            (ball_image, drum_image)
        } else {
            (drum_image, ball_image)
        };
        out.push(SectionCurve {
            curve,
            on_a: Some(on_a),
            on_b: Some(on_b),
            tolerance,
            exact: false,
            closed: true,
            tangential: false,
        });
    }
    Some(out)
}

/// A plane leaning all but along a drum's axis, over the drum's height.
///
/// The closed form is an ellipse whose long axis is the drum's radius over
/// the lean, kilometres for a facet group fitted a hundred-thousandth off
/// a hole's axis. Its parameter spans the few millimetres the drum holds in
/// a millionth of a turn, and crossings solved on it are only as good as
/// that ruler. The crossing is solved instead in the drum's cross-sections
/// along its height and kept as two lines where they hold, as
/// [`near_parallel_drums`] does. `None` where the lean is exactly nothing
/// (the closed form's lines are exact) or more than a thousandth, or where
/// the plane does not cross the drum cleanly all the way up.
fn near_parallel_plane_drum(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<Vec<SectionCurve>> {
    const LEAN: f64 = 1e-3;
    const SPAN: f64 = 3e4;
    let (plane, drum, surface) = match (a, b) {
        (SurfaceGeometry::Plane(p), SurfaceGeometry::Cylinder(c)) => (p.plane(), c.cylinder(), b),
        (SurfaceGeometry::Cylinder(c), SurfaceGeometry::Plane(p)) => (p.plane(), c.cylinder(), a),
        _ => return None,
    };
    let axis = drum.axis();
    let (d, r) = (axis.direction.vector(), drum.radius());
    let n = plane.normal().vector();
    let lean = n.dot(d).abs();
    // Only where the ellipse is thirty metres or more across: there a
    // parameter solved to its last billionth lands tens of nanometres off in
    // space, past the weld of a face with tight edges. A shorter one is
    // ruler enough, and its closed form crosses faster than a fitted curve.
    if lean <= tol.angular() || lean > LEAN || r / lean < SPAN {
        return None;
    }
    let across = n - d * n.dot(d);
    let k = across.magnitude();
    let e1 = across / k;
    let e2 = d.cross(e1);
    let (_, (lo, hi)) = surface.domain();
    if !(lo.is_finite() && hi.is_finite()) || hi - lo <= tol.confusion() {
        return None;
    }
    let meet = |z: f64| -> Option<[Point; 2]> {
        let centre = axis.location + d * z;
        let u = -plane.signed_distance_to(centre) / k;
        let margin = tol.confusion() * 1e3;
        if u.abs() >= r - margin {
            return None;
        }
        let w = r.mul_add(r, -(u * u)).sqrt();
        Some([centre + e1 * u + e2 * w, centre + e1 * u - e2 * w])
    };
    lines_through_stations(lo, hi, meet, 0.0, tol)
}

/// An exact curve dressed as a section, clipped to the surfaces it lies on.
///
/// The analytic layer works on the unbounded geometry (a plane and a cylinder
/// meet in unbounded lines), but the *surfaces* carry finite extents, and a
/// section running a billion units past both is not something an edge can be
/// built on. A line is clipped to the parameter interval where it is inside
/// both extents, through its exact pcurves; a curve wholly outside either
/// extent is dropped, or the boolean above would see a phantom edge on a
/// region the face does not have.
///
/// A *closed* curve partially outside an extent is kept whole: cutting it into
/// arcs is the restriction problem, and the restriction that matters is the
/// face's trim, which is the boolean's job. The extent here is only the
/// surface's parameterization window.
fn exact_section(
    curve: Curve,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<SectionCurve> {
    let closed = match &curve {
        Curve::Circle(_) | Curve::Ellipse(_) => true,
        _ => curve.is_closed(tol),
    };
    let range = curve.domain();
    let on_a = exact_pcurve(&curve, range, a, tol);
    let on_b = exact_pcurve(&curve, range, b, tol);

    if let Curve::Line(_) = &curve {
        // Clip through whichever pcurves exist; a missing pcurve leaves that
        // surface's extent unenforced, which errs long rather than wrong.
        let mut interval = curve.domain();
        if let Some(p) = &on_a {
            interval = intersect_intervals(interval, inside_box(p, a))?;
        }
        if let Some(p) = &on_b {
            interval = intersect_intervals(interval, inside_box(p, b))?;
        }
        let (lo, hi) = interval;
        let Curve::Line(line) = &curve else {
            unreachable!()
        };
        let clipped: Curve = ogeom_geom::LineCurve::over(line.axis(), lo, hi)
            .ok()?
            .into();
        let clip2 = |p: &PlanarCurve| -> Option<PlanarCurve> {
            let PlanarCurve::Line(l) = p else {
                return Some(p.clone());
            };
            Some(Line2d::over(l.axis(), lo, hi).ok()?.into())
        };
        let (ca, cb) = (on_a.as_ref().and_then(clip2), on_b.as_ref().and_then(clip2));
        let tangential = touching_along(&clipped, ca.as_ref(), cb.as_ref(), a, b, tol);
        return Some(SectionCurve {
            on_a: ca,
            on_b: cb,
            tolerance: 0.0,
            exact: true,
            closed: false,
            tangential,
            curve: clipped,
        });
    }

    // A closed curve: dropped only when wholly outside an extent it has a
    // pcurve to check against.
    for (pcurve, surface) in [(&on_a, a), (&on_b, b)] {
        if let Some(p) = pcurve
            && !touches_box(p, surface, tol)
        {
            return None;
        }
    }
    let tangential = touching_along(&curve, on_a.as_ref(), on_b.as_ref(), a, b, tol);
    Some(SectionCurve {
        on_a,
        on_b,
        tolerance: 0.0,
        exact: true,
        closed,
        tangential,
        curve,
    })
}

/// Whether the surfaces touch along an exact curve rather than crossing it:
/// their normals parallel at stations along its length.
///
/// Decided through the curve's own pcurves, which is where the normals can
/// be read without inverting anything. A curve missing a pcurve on either
/// surface is reported as a crossing, the honest default, since a section
/// nobody can place in a chart is one nothing can classify as contact
/// either.
fn touching_along(
    curve: &Curve,
    on_a: Option<&PlanarCurve>,
    on_b: Option<&PlanarCurve>,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> bool {
    // The chart position of a sample: through the pcurve where one exists,
    // through the surface's own closed-form inversion where not. A meridian
    // through a sphere's poles has no pcurve (its longitude jumps half a
    // turn at each pole), but every *point* of it inverts fine, and a
    // tangency that would be missed for want of a pcurve becomes a crossing
    // section lying along a face's own boundary, which is the worst thing a
    // section can be.
    let sample_uv = |pc: Option<&PlanarCurve>,
                     surface: &SurfaceGeometry,
                     t: f64|
     -> Option<ogeom_math::Point2> {
        if let Some(pc) = pc {
            return pc.point_at(t, tol).ok();
        }
        let p = curve.point_at(t, tol).ok()?;
        chart_inversion(surface, p, tol)
    };
    let (lo, hi) = curve.domain();
    // Offsets chosen off the round fractions, so a curve through a chart
    // degeneracy (a meridian's poles sit at quarters of its turn) is
    // sampled beside the degenerate points rather than on them. A sample
    // whose inversion still fails is skipped: the point says nothing,
    // not that the surfaces cross.
    let mut judged = 0_usize;
    for f in [0.07, 0.19, 0.37, 0.53, 0.71, 0.89] {
        let t = (hi - lo).mul_add(f, lo);
        let (Some(ua), Some(ub)) = (sample_uv(on_a, a, t), sample_uv(on_b, b, t)) else {
            continue;
        };
        let (Ok(na), Ok(nb)) = (a.normal_at(ua.x, ua.y, tol), b.normal_at(ub.x, ub.y, tol)) else {
            continue;
        };
        if na.vector().cross(nb.vector()).magnitude() > 1e-6 {
            return false;
        }
        judged += 1;
    }
    judged >= 3
}

/// A point's chart position on an analytic surface, by closed form.
fn chart_inversion(
    surface: &SurfaceGeometry,
    p: ogeom_math::Point,
    tol: Tolerances,
) -> Option<ogeom_math::Point2> {
    use ogeom_math::elementary;
    let (u, v) = match surface {
        SurfaceGeometry::Plane(s) => elementary::plane_parameters(&s.plane(), p),
        SurfaceGeometry::Cylinder(s) => {
            elementary::cylinder_parameters(&s.cylinder(), p, tol).ok()?
        }
        SurfaceGeometry::Cone(s) => elementary::cone_parameters(&s.cone(), p, tol).ok()?,
        SurfaceGeometry::Sphere(s) => elementary::sphere_parameters(&s.sphere(), p, tol).ok()?,
        SurfaceGeometry::Torus(s) => elementary::torus_parameters(&s.torus(), p, tol).ok()?,
        _ => return None,
    };
    Some(ogeom_math::Point2::new(u, v))
}

/// The parameter interval over which a 2D line stays inside a surface's
/// parameter box. `None` when it never enters.
pub(crate) fn inside_box(pcurve: &PlanarCurve, surface: &SurfaceGeometry) -> Option<(f64, f64)> {
    // The pcurve as a point and a rate along its own parameter: a line, or
    // a degree-one spline of two points (a cone's ruling), linear in it.
    let (o, d) = match pcurve {
        PlanarCurve::Line(line) => {
            let axis = line.axis();
            (axis.location, axis.direction.vector())
        }
        PlanarCurve::BSpline(spline)
            if spline.knots().degree() == 1 && spline.control_points().len() == 2 =>
        {
            let (t0, t1) = spline.knots().domain();
            let (p0, p1) = (
                spline.control_points()[0].point(),
                spline.control_points()[1].point(),
            );
            if t1 <= t0 {
                return None;
            }
            let rate = (p1 - p0) / (t1 - t0);
            (p0 - rate * t0, rate)
        }
        _ => return None,
    };
    let ((ua, ub), (va, vb)) = surface.domain();

    // The slab test, one axis at a time.
    let mut lo = f64::NEG_INFINITY;
    let mut hi = f64::INFINITY;
    for (origin, direction, low, high) in [(o.x, d.x, ua, ub), (o.y, d.y, va, vb)] {
        if direction.abs() <= f64::MIN_POSITIVE {
            if origin < low || origin > high {
                return None;
            }
            continue;
        }
        let (a, b) = ((low - origin) / direction, (high - origin) / direction);
        let (near, far) = if a < b { (a, b) } else { (b, a) };
        lo = lo.max(near);
        hi = hi.min(far);
    }
    if lo >= hi {
        return None;
    }
    Some((lo, hi))
}

/// Whether a closed pcurve may pass through the surface's box.
pub(crate) fn touches_box(
    pcurve: &PlanarCurve,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> bool {
    use ogeom_geom::Curve2d;
    let ((ua, ub), (va, vb)) = surface.domain();
    let (lo, hi) = pcurve.domain();
    // Asked of the spans between samples, not the samples alone: a plane all
    // but parallel to a cylinder's axis meets it in an ellipse kilometres
    // long, whose image on the cylinder's chart sweeps through a window a few
    // millimetres tall in a sliver of its turn, between any two samples.
    // Each span is taken as its chord's box widened by the chord's length,
    // which holds the curve between them wherever it bends no tighter than
    // the samples are apart. Kept wrongly, a curve costs a section the trim
    // then cuts to nothing; dropped wrongly, the faces never split.
    const SPANS: u32 = 64;
    let points: Vec<Option<ogeom_math::Point2>> = (0..=SPANS)
        .map(|i| {
            pcurve
                .point_at(lo + (hi - lo) * f64::from(i) / f64::from(SPANS), tol)
                .ok()
        })
        .collect();
    points.windows(2).any(|pair| {
        let (Some(p), Some(q)) = (pair[0], pair[1]) else {
            return false;
        };
        let pad = p.distance(q);
        // Periodic directions always contain; only a bounded one excludes.
        let u_ok =
            surface.is_periodic_u() || (p.x.max(q.x) + pad >= ua && p.x.min(q.x) - pad <= ub);
        let v_ok =
            surface.is_periodic_v() || (p.y.max(q.y) + pad >= va && p.y.min(q.y) - pad <= vb);
        u_ok && v_ok
    })
}

/// The overlap of two intervals. `None` when they miss.
fn intersect_intervals(a: (f64, f64), b: Option<(f64, f64)>) -> Option<(f64, f64)> {
    let b = b?;
    let (lo, hi) = (a.0.max(b.0), a.1.min(b.1));
    if lo >= hi {
        return None;
    }
    Some((lo, hi))
}

/// The general path: seed, trace, fit.
fn marched(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    windows: (Option<ParameterWindow>, Option<ParameterWindow>),
    options: IntersectOptions,
    tol: Tolerances,
) -> OgeomResult<SurfaceIntersection> {
    // Seeding and tracing walk the windowed views; everything after reads
    // the surfaces as given, whose chart the views share.
    let view_a = windows.0.and_then(|w| windowed(a, w, tol));
    let view_b = windows.1.and_then(|w| windowed(b, w, tol));
    let traced = branches_from_touches(
        view_a.as_ref().unwrap_or(a),
        view_b.as_ref().unwrap_or(b),
        options.marching,
        tol,
    )?;
    if traced.is_empty() {
        return Ok(SurfaceIntersection::Apart);
    }
    let mut out = Vec::with_capacity(traced.len());
    let mut contacts: Vec<crate::march::Traced> = Vec::new();
    // A stretch of a border of either chart lying in the other surface
    // runs on that border, each point at the border's own parameter. It
    // was kept only where the surface leaves the other at an angle there,
    // which it is not asked again.
    let domains = (
        view_a.as_ref().unwrap_or(a).domain(),
        view_b.as_ref().unwrap_or(b).domain(),
    );
    let on_border = |params: &[(f64, f64)], ((u0, u1), (v0, v1)): ((f64, f64), (f64, f64))| {
        let Some(&(u, v)) = params.first() else {
            return false;
        };
        ((u == u0 || u == u1) && params.iter().all(|p| p.0 == u))
            || ((v == v0 || v == v1) && params.iter().all(|p| p.1 == v))
    };
    for (branch, crossing) in &traced {
        // A branch along which the two surfaces share their normal is a
        // tangency, not a crossing: the marcher's seeding cannot tell the
        // noise floor of a tangential valley from a genuine sign change, and
        // what it traces there is a stalled fragment of the valley, not a
        // section. The valley is still a curve, though, and the tangential
        // walker is the one that can follow it, so the fragment becomes a
        // seed rather than a discard, and what comes back is marked as
        // contact so nobody classifies by it. A branch leaving a touch
        // along an arm, where the section crosses itself, is a crossing
        // however shallow: near the touch no sample of it crosses
        // measurably, and the valley walker would follow some other curve
        // within its reach of both.
        let lying = on_border(&branch.on_a, domains.0) || on_border(&branch.on_b, domains.1);
        if !lying && !crossing && branch_is_tangential(a, b, branch, tol)? {
            if let Some(contact) = walk_contact(a, b, branch, &contacts, options.marching, tol)? {
                contacts.push(contact);
            }
            continue;
        }
        if branch.stopped == crate::march::Stopped::RanOut {
            ogeom_bail!(
                NotDone,
                "a marched section ran out of its point budget before \
                 finishing; the seam is longer than the chord affords and \
                 fitting the truncation would state a curve that is not there"
            );
        }
        // A fit past its budget is still honest data: the error it reached
        // is carried on the record and every consumer widens by it: an
        // imported part's ragged pair can trace branches nothing fits, and
        // those sections fall outside every trim downstream. Only a trace
        // cut off by the point budget, refused above, states a curve that
        // is not there. (A boolean marching an *exact* pair whose image has
        // no closed form holds its own marched sections to a budget, in
        // its own fallback, where a miss is a miss.)
        for fitted in fitted_in_pieces(a, b, branch, options.tolerance, tol)? {
            out.push(SectionCurve {
                curve: fitted.curve.into(),
                on_a: Some(fitted.on_a.into()),
                on_b: Some(fitted.on_b.into()),
                // The sum of the stated parts: the trace is within its chord of
                // the truth, the fit within its error of the trace.
                tolerance: options.marching.chord + fitted.fit_error,
                exact: false,
                closed: fitted.closed,
                tangential: false,
            });
        }
    }
    for contact in &contacts {
        let fitted = approximate_branch(a, b, contact, options.tolerance, tol)?;
        out.push(SectionCurve {
            curve: fitted.curve.into(),
            on_a: Some(fitted.on_a.into()),
            on_b: Some(fitted.on_b.into()),
            tolerance: options.marching.chord + fitted.fit_error,
            exact: false,
            closed: fitted.closed,
            tangential: true,
        });
    }
    if out.is_empty() {
        return Ok(SurfaceIntersection::Apart);
    }
    Ok(SurfaceIntersection::Along(out))
}

/// A surface narrowed to a window of its parameters, or `None` where the
/// window leaves it whole: see [`intersect_surfaces_within`].
fn windowed(
    surface: &SurfaceGeometry,
    window: ParameterWindow,
    tol: Tolerances,
) -> Option<SurfaceGeometry> {
    let (du, dv) = surface.domain();
    let eps = tol.parametric();
    // The part of `want` inside `have`, or `None` where a periodic
    // direction's window spans a period. A periodic window may straddle the
    // seam: the trim holds the parameter inside the window and the basis
    // wraps it, so the walk reads the window in the face's own terms.
    let narrow = |want: (f64, f64), have: (f64, f64), periodic: bool| -> Option<(f64, f64)> {
        if periodic {
            return (want.1 - want.0 < have.1 - have.0 - eps).then_some(want);
        }
        Some((want.0.max(have.0), want.1.min(have.1)))
    };
    let finite = |w: (f64, f64)| w.0.is_finite() && w.1.is_finite() && w.1 > w.0;
    if !finite(window.0) || !finite(window.1) {
        return None;
    }
    let narrower = |got: (f64, f64), have: (f64, f64)| {
        got.1 - got.0 < have.1 - have.0 - eps || got.0 > have.0 + eps || got.1 < have.1 - eps
    };
    let u = narrow(window.0, du, surface.is_periodic_u()).filter(|&u| narrower(u, du));
    let v = narrow(window.1, dv, surface.is_periodic_v()).filter(|&v| narrower(v, dv));
    if u.is_some_and(|u| u.1 - u.0 <= eps) || v.is_some_and(|v| v.1 - v.0 <= eps) {
        return None;
    }
    // A drum's height narrows in its own terms, keeping the turn periodic
    // where the window cannot cut it: a trim does not wrap.
    let base: Option<SurfaceGeometry> = match (surface, v) {
        (SurfaceGeometry::Cylinder(c), Some(v)) => {
            ogeom_geom::CylinderSurface::new(c.cylinder(), v)
                .ok()
                .map(Into::into)
        }
        (SurfaceGeometry::Cone(c), Some(v)) => ogeom_geom::ConeSurface::new(c.cone(), v)
            .ok()
            .map(Into::into),
        _ => None,
    };
    let trim = |basis: SurfaceGeometry, u: (f64, f64), v: (f64, f64)| {
        ogeom_geom::TrimmedSurface::new(basis, u, v, tol)
            .ok()
            .map(Into::into)
    };
    match (base, u) {
        (Some(base), Some(u)) => {
            let v = base.domain().1;
            trim(base, u, v)
        }
        (Some(base), None) => Some(base),
        (None, _) => {
            // A periodic direction the window leaves whole would stop
            // wrapping under a trim.
            let wraps_u = surface.is_periodic_u() && u.is_none();
            let wraps_v = surface.is_periodic_v() && v.is_none();
            if (u.is_none() && v.is_none()) || wraps_u || wraps_v {
                return None;
            }
            trim(surface.clone(), u.unwrap_or(du), v.unwrap_or(dv))
        }
    }
}

/// A traced branch fitted, in pieces where whole it will not fit.
///
/// A trace winding several turns round a drum (a thread's flank meeting a
/// bore) is long and turns the same way throughout, and one fit of it can
/// run out of room and come back with an error of the drum's size. An open
/// branch whose fit misses by more than a hundred times its tolerance, or
/// strays farther from the trace than the finest step the marcher took
/// along it, is split at its middle sample and each half fitted the same
/// way, down to a floor of samples and depth; the pieces meet at the shared
/// sample. The marcher shortens its step where the section turns sharply
/// (round a neck a thousandth wide where a drill all but touches a torus),
/// and a fit missing by more than that step there wobbles across the neck.
/// A closed branch is split the same way, its two halves open and meeting
/// at both ends: a loop round a thin drum lying all but tangent inside a
/// wider one turns sharply at its tip, and fitted whole it can come back
/// off the trace by the drum's size. A fit that misses its tolerance by
/// less stands whole, its error stated: a caller takes one curve per branch
/// where it can, and a few microns do not warrant more. So does a branch no
/// split helps.
fn fitted_in_pieces(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    branch: &crate::march::Traced,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<crate::approx::IntersectionCurve>> {
    const DEPTH: u32 = 6;
    const FLOOR: usize = 16;
    fn go(
        a: &SurfaceGeometry,
        b: &SurfaceGeometry,
        branch: &crate::march::Traced,
        tolerance: f64,
        depth: u32,
        tol: Tolerances,
    ) -> OgeomResult<Vec<crate::approx::IntersectionCurve>> {
        let whole = approximate_branch(a, b, branch, tolerance, tol)?;
        // The finest step the marcher took between the branch's ends: the
        // first and last steps land on a patch edge or close the loop, and
        // their length says nothing of how sharply the section turns.
        let points = &branch.points;
        let inner = points.get(1..points.len().saturating_sub(1)).unwrap_or(&[]);
        let step = inner
            .windows(2)
            .map(|w| w[0].distance(w[1]))
            .fold(f64::INFINITY, f64::min);
        if whole.met
            || whole.fit_error <= step.min(tolerance * 1e2)
            || depth == 0
            || branch.points.len() < 2 * FLOOR
        {
            return Ok(vec![whole]);
        }
        let middle = branch.points.len() / 2;
        // A half of a loop is open: it stops where the other half starts.
        let stopped = if branch.closed() {
            crate::march::Stopped::Stalled
        } else {
            branch.stopped
        };
        let half = |range: core::ops::RangeInclusive<usize>| crate::march::Traced {
            points: branch.points[range.clone()].to_vec(),
            on_a: branch.on_a[range.clone()].to_vec(),
            on_b: branch.on_b[range].to_vec(),
            stopped,
        };
        let mut pieces = go(a, b, &half(0..=middle), tolerance, depth - 1, tol)?;
        pieces.extend(go(
            a,
            b,
            &half(middle..=branch.points.len() - 1),
            tolerance,
            depth - 1,
            tol,
        )?);
        // Worse in pieces than whole (a trace that is noise, not length):
        // the whole stands.
        let worst = pieces.iter().map(|p| p.fit_error).fold(0.0_f64, f64::max);
        Ok(if worst < whole.fit_error {
            pieces
        } else {
            vec![whole]
        })
    }
    go(a, b, branch, tolerance, DEPTH, tol)
}

/// Follow the contact a tangential fragment sits on, unless one already
/// traced covers it.
///
/// A tangential valley hands the crossing marcher several stalled fragments
/// (the seeds converge onto the contact from wherever they started and
/// wander there), so the fragments are candidates for *one* curve, not
/// several. A fragment whose middle already lies on a traced contact is one
/// of those repeats.
fn walk_contact(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    fragment: &crate::march::Traced,
    already: &[crate::march::Traced],
    marching: Marching,
    tol: Tolerances,
) -> OgeomResult<Option<crate::march::Traced>> {
    let middle = fragment.points.len() / 2;
    let Some(point) = fragment.points.get(middle).copied() else {
        return Ok(None);
    };
    for traced in already {
        // Traced points sit a step apart, so "on this curve" has to allow
        // half a step of gap to the nearest sample plus the chord budget.
        let spacing = traced
            .points
            .windows(2)
            .map(|w| w[0].distance(w[1]))
            .fold(0.0f64, f64::max);
        let near = traced
            .points
            .iter()
            .map(|p| p.distance(point))
            .fold(f64::INFINITY, f64::min);
        if near <= spacing.mul_add(0.5, marching.chord.max(tol.confusion())) {
            return Ok(None);
        }
    }
    let seed = crate::march::Contact {
        point,
        on_a: fragment.on_a[middle],
        on_b: fragment.on_b[middle],
    };
    // The walker refuses a seed that is not a contact; that refusal is an
    // answer, not a failure: the fragment simply had nothing to follow.
    // A walk that stalls where it started says the same thing in points:
    // too few to fit, so there is no contact curve to report here.
    Ok(trace_tangential(a, b, seed, marching, tol)
        .ok()
        .filter(|traced| traced.points.len() >= 4))
}

/// Whether a traced branch runs along a tangency of the two surfaces:
/// their normals parallel, sampled along its length.
fn branch_is_tangential(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    branch: &crate::march::Traced,
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_geom::Surface as _;
    let count = branch.points.len();
    if count == 0 {
        return Ok(true);
    }
    let exact = exact_surface(a) && exact_surface(b);
    for k in 0..5 {
        let i = ((k * (count - 1)) / 4).min(count - 1);
        let (ua, va) = branch.on_a[i];
        let (ub, vb) = branch.on_b[i];
        let (dau, dav) = a.d1_at(ua, va, tol)?;
        let (dbu, dbv) = b.d1_at(ub, vb, tol)?;
        let na = dau.cross(dav);
        let nb = dbu.cross(dbv);
        let (ma, mb) = (na.magnitude(), nb.magnitude());
        if ma <= tol.confusion() || mb <= tol.confusion() {
            continue;
        }
        // The threshold carries the fitted world: a blend surface within a
        // fit tolerance of true tangency crosses its host at an angle that
        // grows as the square root of that tolerance, and calling such a
        // graze transversal splits faces along slivers no classifier can
        // hold. Genuinely transversal analytic pairs meeting under two
        // degrees are the pathology, not the rule.
        if na.cross(nb).magnitude() / (ma * mb) > TANGENTIAL_SINE {
            return Ok(false);
        }
        // Two surfaces that are their own geometry hold no fit to graze
        // within: under the threshold they are tangent only where they do
        // not measurably cross.
        let at = crate::march::Contact {
            point: branch.points[i],
            on_a: (ua, va),
            on_b: (ub, vb),
        };
        if exact && crate::march::crosses_measurably(a, b, at, tol) {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether a surface is its own geometry rather than what may be a fit of
/// another: an analytic surface, or one swept from an analytic curve. A
/// patch, an offset, or a sweep of a spline may be a fit, tangent to its
/// neighbours only to within it.
fn exact_surface(surface: &SurfaceGeometry) -> bool {
    fn exact_curve(curve: &Curve) -> bool {
        match curve {
            Curve::Line(_)
            | Curve::Circle(_)
            | Curve::Ellipse(_)
            | Curve::Hyperbola(_)
            | Curve::Parabola(_)
            | Curve::Helix(_) => true,
            Curve::Trimmed(t) => exact_curve(t.basis()),
            _ => false,
        }
    }
    match surface {
        SurfaceGeometry::Plane(_)
        | SurfaceGeometry::Cylinder(_)
        | SurfaceGeometry::Cone(_)
        | SurfaceGeometry::Sphere(_)
        | SurfaceGeometry::Torus(_) => true,
        SurfaceGeometry::Extrusion(e) => exact_curve(e.curve()),
        SurfaceGeometry::Revolution(r) => exact_curve(r.curve()),
        SurfaceGeometry::Trimmed(t) => exact_surface(t.basis()),
        _ => false,
    }
}

/// The exact pcurve of a curve lying on a surface, where the projection has
/// a closed form; `None` where it does not.
///
/// Public because the boolean's same-domain handling needs it: two faces on
/// one geometric surface may still carry different charts, and the other
/// face's boundary edges have to be spoken in this face's parameters before
/// they can split it.
#[must_use]
pub fn exact_pcurve_of(
    curve: &Curve,
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    exact_pcurve(curve, curve.domain(), surface, tol)
}

/// As [`exact_pcurve_of`], with the parameter range the caller actually
/// uses.
///
/// A curve's chart image can depend on *which part* of the curve is meant: a
/// ruling on a cone crosses the apex, and its angle on the far nappe is half
/// a turn from its angle on the near one. The curve's own domain may span
/// both (an imported line's usually does), so a caller that knows its edge's
/// range must say so, or the exact projection may answer for the wrong side.
#[must_use]
pub fn exact_pcurve_over(
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    exact_pcurve(curve, range, surface, tol)
}

/// The exact pcurve of an analytic curve on an analytic surface, where the
/// projection has a closed form.
///
/// Same-parameter by construction: each 2D curve inherits the 3D curve's own
/// parameterization, so the two evaluate to the same point of the intersection
/// at the same `t`. The cases are the ones where that inheritance is exact;
/// anything else returns `None` rather than a fit, because an *exact* result
/// with a fitted pcurve would be a curve whose descriptions disagree by an
/// amount nothing on it records.
fn exact_pcurve(
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    // A trim is a statement about *where* on a curve, not about what it is:
    // the basis carries the shape and the trim shares its parameter, so the
    // pcurve is the basis's own pcurve trimmed the same way. Answered here
    // rather than in every surface's own case, because the answer does not
    // depend on the surface at all. A *reversed* trim reads its basis at
    // `s + e - u` over its window `[s, e]`, and so does the planar trim
    // reversed, so the range the caller uses is mirrored the same way
    // before it reaches the basis.
    if let Curve::Trimmed(trimmed) = curve {
        let window = ogeom_geom::Curve3d::domain(&**trimmed);
        let range = if trimmed.is_reversed() {
            let sum = window.0 + window.1;
            (sum - range.1, sum - range.0)
        } else {
            range
        };
        let basis = exact_pcurve(trimmed.basis(), range, surface, tol)?;
        let pcurve: PlanarCurve = ogeom_geom::Trimmed2d::new(basis, window.0, window.1, tol)
            .ok()?
            .into();
        return Some(if trimmed.is_reversed() {
            ogeom_geom::Reversible::reversed(&pcurve)
        } else {
            pcurve
        });
    }
    // The cases below read a circle's or an ellipse's frame as the way its
    // parameter runs, so one running backwards is read as its forward
    // spelling, the same point at every parameter.
    let forward = match curve {
        Curve::Circle(c) if c.is_reversed() => Some(Curve::Circle(c.forward(tol).ok()?)),
        Curve::Ellipse(e) if e.is_reversed() => Some(Curve::Ellipse(e.forward(tol).ok()?)),
        _ => None,
    };
    if let Some(forward) = forward {
        return exact_pcurve(&forward, range, surface, tol);
    }
    // The circle cases read the frame's `z` as the axis the parameter turns
    // about, which holds on a right-handed frame only. A mirrored circle
    // stands on a left-handed one, its parameter turning from `x` towards a
    // `y` that runs clockwise about `z`, and is read as its right-handed
    // spelling: the same `x` and `y`, `z` reversed, the same point at every
    // parameter.
    if let Curve::Circle(c) = curve
        && c.circle().frame().handedness() == ogeom_math::Handedness::Left
    {
        let frame = c.circle().frame();
        let right =
            ogeom_math::Frame::from_axes(frame.origin(), frame.x(), frame.y(), -frame.z(), tol)
                .ok()?;
        let circle = ogeom_math::Circle::new(right, c.circle().radius(), tol).ok()?;
        return exact_pcurve(
            &Curve::Circle(ogeom_geom::CircleCurve::new(circle)),
            range,
            surface,
            tol,
        );
    }
    if let (Curve::Offset(offset), SurfaceGeometry::Plane(p)) = (curve, surface) {
        return offset_on_plane(offset, range, surface, p.plane(), tol);
    }
    match surface {
        SurfaceGeometry::Plane(p) => on_plane(curve, p.plane(), tol),
        SurfaceGeometry::Cylinder(c) => on_cylinder(curve, range, c.cylinder(), tol),
        SurfaceGeometry::Sphere(s) => on_sphere(curve, range, s.sphere(), tol),
        SurfaceGeometry::Torus(t) => on_torus(curve, range, t.torus(), tol),
        SurfaceGeometry::Cone(c) => on_cone(curve, range, c.cone(), tol),
        _ => None,
    }
}

/// The parameter window of a circle's pcurve that inherits the circle's own
/// parameter: the whole turn, widened to the range the edge uses. An arc
/// whose range runs past the turn's end (its start late in the turn, its
/// end past the circle's origin) is read on the same straight chart line
/// continued, not wrapped.
fn circle_window(range: (f64, f64)) -> (f64, f64) {
    let tau = core::f64::consts::TAU;
    if range.0.is_finite() && range.1.is_finite() {
        (range.0.min(range.1).min(0.0), range.0.max(range.1).max(tau))
    } else {
        (0.0, tau)
    }
}

/// The axis an analytic surface's angle `u` turns about: `x × y`, which is
/// the frame's `z` on a right-handed frame and `-z` on a left-handed
/// (mirrored) one. A circle round the axis runs up `u` where it turns about
/// this axis and down `u` where it turns against it.
fn spin(frame: &ogeom_math::Frame) -> ogeom_math::Vector {
    frame.z().vector() * frame.handedness().sign()
}

/// The pcurve of a curve on a cone, for the two straight-line families.
///
/// A ruling (through the apex, on the surface) runs at constant `u`; a
/// circle perpendicular to the axis, centred on it, with the radius the cone
/// has at that height, runs at constant `v`. Both inherit the 3D curve's own
/// parameter, the circle with phase and winding exactly as the cylinder case.
/// The ruling's angle is measured over `range`, because the same line has
/// the opposite angle on the other side of the apex.
fn on_cone(
    curve: &Curve,
    range: (f64, f64),
    cone: ogeom_math::Cone,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    let frame = cone.frame();
    let axis_z = frame.z().vector();
    let tau = core::f64::consts::TAU;
    match curve {
        Curve::Circle(c) => {
            let circle = c.circle();
            if circle.frame().z().vector().cross(axis_z).magnitude() > tol.angular() {
                return None;
            }
            let local = frame.to_local(circle.centre());
            if local.x.hypot(local.y) > tol.confusion() {
                return None;
            }
            // The cone's radius at the circle's height must be the circle's,
            // or (past the apex, on the far nappe, where the radius runs
            // negative) its negative: the same parallel half a turn round.
            let expected = cone
                .half_angle()
                .tan()
                .mul_add(local.z, cone.reference_radius());
            let turned = if (expected - circle.radius()).abs() <= tol.confusion() * 10.0 {
                0.0
            } else if (expected + circle.radius()).abs() <= tol.confusion() * 10.0 {
                core::f64::consts::PI
            } else {
                return None;
            };
            let start = circle.centre() + circle.frame().x().vector() * circle.radius();
            let at = frame.to_local(start);
            let phase = at.y.atan2(at.x) + turned;
            let winding = circle.frame().z().vector().dot(spin(&frame)).signum();
            let towards =
                ogeom_math::Direction2::new(ogeom_math::Vector2::new(winding, 0.0), tol).ok()?;
            Some(
                Line2d::over(
                    ogeom_math::Axis2::new(Point2::new(phase, local.z), towards),
                    circle_window(range).0,
                    circle_window(range).1,
                )
                .ok()?
                .into(),
            )
        }
        Curve::Line(line) => {
            // A ruling: verified by sample, not assumed: three points on
            // the surface pin a line to it.
            let axis = line.axis();
            let on = |t: f64| {
                let p = axis.location + axis.direction.vector() * t;
                cone.distance_to(p) <= tol.confusion() * 10.0
            };
            if !on(0.0) || !on(1.0) || !on(-1.0) {
                return None;
            }
            // A ruling reaching the tip may be *stated* from the apex
            // itself (where the angle is atan2(0, 0), garbage) and its
            // own domain usually spans both nappes, where the angles differ
            // by half a turn. Measure the angle at whichever end of the
            // *used* range stands farthest from the axis: that is the side
            // the caller means.
            let (lo, hi) = if range.0.is_finite() && range.1.is_finite() && range.0 != range.1 {
                range
            } else {
                line.domain()
            };
            // Only the used range votes. The line's own origin is stated
            // wherever the file likes (some writers park it hundreds of
            // kilometres down the infinite line, past the apex on the other
            // nappe), and letting it compete reads the angle half a turn
            // from the side the edge actually uses.
            let mut local: Option<ogeom_math::Point> = None;
            for t in [lo, hi] {
                if !t.is_finite() {
                    continue;
                }
                let candidate = frame.to_local(axis.location + axis.direction.vector() * t);
                if local.is_none_or(|held| candidate.x.hypot(candidate.y) > held.x.hypot(held.y)) {
                    local = Some(candidate);
                }
            }
            let local = local?;
            if local.x.hypot(local.y) <= tol.confusion() {
                return None;
            }
            // The chart's point at angle `u` and height `z` stands at the
            // radius `radius_at(z)` along `u`'s direction, and past the apex
            // that radius runs negative: a point on the far nappe is on the
            // column half a turn from its own bearing. The whole ruling,
            // through the apex or not, is that one column.
            let side = cone.radius_at(local.z).signum();
            let u = (local.y * side).atan2(local.x * side).rem_euclid(tau);
            // Same-parameter exactly: a degree-one spline over the used
            // range maps t linearly onto the chart column, whatever rate
            // the slant climbs at.
            let v_at = |t: f64| {
                frame
                    .to_local(axis.location + axis.direction.vector() * t)
                    .z
            };
            let knots = ogeom_math::KnotVector::new(vec![lo, lo, hi, hi], 1).ok()?;
            Some(
                ogeom_geom::BSpline2d::new(
                    knots,
                    vec![Point2::new(u, v_at(lo)), Point2::new(u, v_at(hi))],
                    tol,
                )
                .ok()?
                .into(),
            )
        }
        _ => None,
    }
}

/// The pcurve of a circle on a torus, for the two families that are straight
/// lines in `(u, v)`.
///
/// A *parallel* (centred on the axis, in a plane perpendicular to it) runs
/// at constant `v`; a *tube circle* (minor radius, centred on the tube's
/// spine, in a plane through the axis) runs at constant `u`. Both inherit
/// the circle's own angle, phase and winding included, exactly as the
/// cylinder case does. Fillet faces are tori more often than not, so the
/// STEP reader is the chief consumer.
fn on_torus(
    curve: &Curve,
    range: (f64, f64),
    torus: ogeom_math::Torus,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    let Curve::Circle(c) = curve else {
        return None;
    };
    let circle = c.circle();
    let frame = torus.frame();
    let axis_z = frame.z().vector();
    let normal = circle.frame().z().vector();
    let local = frame.to_local(circle.centre());

    // A parallel of the sweep.
    if normal.cross(axis_z).magnitude() <= tol.angular()
        && local.x.hypot(local.y) <= tol.confusion()
    {
        let sin_v = local.z / torus.minor_radius();
        // On its own side of the axis, or (on a spindle, whose tube swallows
        // the axis) on the tube's folded half past it, where the sweep's
        // radius runs negative: the same parallel half a turn round.
        let (cos_v, turned) = [
            (circle.radius() - torus.major_radius(), 0.0),
            (
                -circle.radius() - torus.major_radius(),
                core::f64::consts::PI,
            ),
        ]
        .into_iter()
        .map(|(reach, turned)| (reach / torus.minor_radius(), turned))
        .find(|(cos_v, _)| (sin_v.hypot(*cos_v) - 1.0).abs() <= tol.confusion())?;
        let v = sin_v.atan2(cos_v);
        let start = circle.centre() + circle.frame().x().vector() * circle.radius();
        let at = frame.to_local(start);
        let phase = at.y.atan2(at.x) + turned;
        let winding = normal.dot(spin(&frame)).signum();
        let towards =
            ogeom_math::Direction2::new(ogeom_math::Vector2::new(winding, 0.0), tol).ok()?;
        return Some(
            Line2d::over(
                ogeom_math::Axis2::new(Point2::new(phase, v), towards),
                circle_window(range).0,
                circle_window(range).1,
            )
            .ok()?
            .into(),
        );
    }

    // A circle of the tube.
    if (circle.radius() - torus.minor_radius()).abs() <= tol.confusion()
        && normal.dot(axis_z).abs() <= tol.angular()
        && (local.x.hypot(local.y) - torus.major_radius()).abs() <= tol.confusion()
        && local.z.abs() <= tol.confusion()
    {
        let u = local.y.atan2(local.x);
        let radial = frame.x().vector() * u.cos() + frame.y().vector() * u.sin();
        let xc = circle.frame().x().vector();
        let phase = xc.dot(axis_z).atan2(xc.dot(radial));
        let winding = normal.dot(radial.cross(axis_z)).signum();
        let towards =
            ogeom_math::Direction2::new(ogeom_math::Vector2::new(0.0, winding), tol).ok()?;
        return Some(
            Line2d::over(
                ogeom_math::Axis2::new(Point2::new(u, phase), towards),
                circle_window(range).0,
                circle_window(range).1,
            )
            .ok()?
            .into(),
        );
    }
    None
}

/// The pcurve of an offset curve on a plane its reference stands square
/// to: the basis's own pcurve, offset in the chart.
///
/// The offset moves each point along the tangent crossed with the
/// reference. Within the plane, the tangent crossed with the chart's
/// `x × y` is the tangent turned a quarter right, so the chart offset is
/// the distance along the left normal with the sign flipped, and flipped
/// again where the reference runs against `x × y`. A reference leaning `ε`
/// off the normal turns the offset direction by at most `ε`, moving the
/// point by at most `2 |d| sin(ε / 2)`, under `|d| tan ε`; the offset is
/// read in the plane only where that stays within the confusion distance.
/// Same-parameter as the basis's pcurve is.
fn offset_on_plane(
    offset: &ogeom_geom::OffsetCurve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    plane: ogeom_math::Plane,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    let frame = plane.frame();
    let up = frame.x().vector().cross(frame.y().vector());
    let reference = offset.reference().vector();
    let (cos, sin) = (reference.dot(up), reference.cross(up).magnitude());
    if cos == 0.0 || offset.distance().abs() * sin > tol.confusion() * cos.abs() {
        return None;
    }
    let basis = exact_pcurve(offset.basis(), range, surface, tol)?;
    Offset2d::new(basis, -cos.signum() * offset.distance())
        .ok()
        .map(|o| PlanarCurve::Offset(Box::new(o)))
}

/// Project a curve lying in a plane into the plane's own coordinates.
///
/// Exact for a line, a circle, an ellipse, a hyperbola, a B-spline, a
/// parabola and a curve drawn on this plane under another chart: the
/// plane's frame is orthonormal, so lengths and the curves' own
/// parameterizations survive the projection unchanged.
fn on_plane(curve: &Curve, plane: ogeom_math::Plane, tol: Tolerances) -> Option<PlanarCurve> {
    let frame = plane.frame();
    let flat = |p: Point| {
        let local = frame.to_local(p);
        Point2::new(local.x, local.y)
    };
    // A direction is turned into the plane's axes alone: carried to a
    // point off the frame's origin and back, it would pick up the
    // rounding of the origin's coordinates, which far from the world's
    // origin tips it past the angular tolerance.
    let flat_direction = |d: ogeom_math::Direction| {
        let local = frame.vector_to_local(d.vector());
        ogeom_math::Direction2::new(ogeom_math::Vector2::new(local.x, local.y), tol).ok()
    };
    match curve {
        Curve::Line(line) => {
            let axis = line.axis();
            let through = flat(axis.location);
            let direction = flat_direction(axis.direction)?;
            let (lo, hi) = line.domain();
            Some(
                Line2d::over(ogeom_math::Axis2::new(through, direction), lo, hi)
                    .ok()?
                    .into(),
            )
        }
        Curve::Circle(c) => {
            let circle = c.circle();
            let frame2 = Frame2::from_axes(
                flat(circle.centre()),
                flat_direction(circle.frame().x())?,
                flat_direction(circle.frame().y())?,
                tol,
            )
            .ok()?;
            Some(Circle2d::new(Circle2::new(frame2, circle.radius(), tol).ok()?).into())
        }
        Curve::Ellipse(e) => {
            let ellipse = e.ellipse();
            let frame2 = Frame2::from_axes(
                flat(ellipse.centre()),
                flat_direction(ellipse.frame().x())?,
                flat_direction(ellipse.frame().y())?,
                tol,
            )
            .ok()?;
            Some(
                Ellipse2d::new(
                    Ellipse2::new(frame2, ellipse.major_radius(), ellipse.minor_radius(), tol)
                        .ok()?,
                )
                .into(),
            )
        }
        Curve::Hyperbola(h) => {
            // Run backwards over a domain symmetric about zero, the branch
            // is the forward one with `y` turned round: `cosh` is even and
            // `sinh` odd, and the chart image runs forward, as the `.brep`
            // format can state it. Over any other
            // domain the reversal reads the branch at `lo + hi - u`, a
            // shift no frame states, and the chart image is reversed over
            // the same domain, which reads it there too.
            let hyperbola = h.hyperbola();
            let (lo, hi) = h.domain();
            let turned = h.is_reversed() && lo + hi == 0.0;
            let y = if turned {
                -hyperbola.frame().y()
            } else {
                hyperbola.frame().y()
            };
            let frame2 = Frame2::from_axes(
                flat(hyperbola.centre()),
                flat_direction(hyperbola.frame().x())?,
                flat_direction(y)?,
                tol,
            )
            .ok()?;
            let flat_hyperbola = ogeom_math::Hyperbola2::new(
                frame2,
                hyperbola.major_radius(),
                hyperbola.minor_radius(),
                tol,
            )
            .ok()?;
            let image: PlanarCurve = ogeom_geom::Hyperbola2d::over(flat_hyperbola, lo, hi)
                .ok()?
                .into();
            Some(if h.is_reversed() && !turned {
                ogeom_geom::Reversible::reversed(&image)
            } else {
                image
            })
        }
        Curve::BSpline(b) => {
            // Affine invariance: a (rational) B-spline in the plane projects
            // into the plane's own coordinates control point by control
            // point, knots and weights untouched: exact, and same-parameter
            // by construction.
            let control = b
                .control_points()
                .iter()
                .map(|w| ogeom_math::Weighted::new(flat((*w).point()), w.weight, tol))
                .collect::<Result<Vec<_>, _>>()
                .ok()?;
            Some(
                ogeom_geom::BSpline2d::rational(b.knots().clone(), control)
                    .ok()?
                    .into(),
            )
        }
        Curve::Parabola(parabola) => {
            // Quadratic in its own parameter, run either way: the one
            // quadratic Bézier over its domain from end to end, its middle
            // control point where the end tangents meet, is the same curve
            // at every parameter.
            let (lo, hi) = parabola.domain();
            let start = parabola.point_at(lo, tol).ok()?;
            let middle = start + parabola.d1_at(lo, tol).ok()? * ((hi - lo) / 2.0);
            let end = parabola.point_at(hi, tol).ok()?;
            let knots = ogeom_math::KnotVector::new(vec![lo, lo, lo, hi, hi, hi], 2).ok()?;
            Some(
                ogeom_geom::BSpline2d::new(knots, vec![flat(start), flat(middle), flat(end)], tol)
                    .ok()?
                    .into(),
            )
        }
        Curve::OnSurface(on) => {
            // Drawn on this plane under another chart: the other chart's
            // origin and axes, read in this one, are a frame, and the
            // pcurve moved by it is the same point at every parameter. A
            // chart on a plane tipped or standing off this one is not this
            // plane's, and is left alone.
            let chart = chart_plane(on.surface())?.frame();
            if !chart.z().is_parallel(frame.z(), tol)
                || plane.distance_to(chart.origin()) > tol.confusion()
            {
                return None;
            }
            let frame2 = Frame2::from_axes(
                flat(chart.origin()),
                flat_direction(chart.x())?,
                flat_direction(chart.y())?,
                tol,
            )
            .ok()?;
            on.pcurve()
                .transformed(&ogeom_math::Transform2::from_frame(&frame2), tol)
                .ok()
        }
        _ => None,
    }
}

/// The plane whose frame's `x` and `y` are a surface's chart, where it has
/// one: a plane, or a plane trimmed, which keeps its basis's chart.
fn chart_plane(surface: &SurfaceGeometry) -> Option<ogeom_math::Plane> {
    match surface {
        SurfaceGeometry::Plane(p) => Some(p.plane()),
        SurfaceGeometry::Trimmed(t) => chart_plane(t.basis()),
        _ => None,
    }
}

/// The pcurve of a curve on a cylinder, where it is a straight line in
/// parameter space.
///
/// A line along the axis runs at constant `u`; a full circle around it runs at
/// constant `v`. Both are lines in `(u, v)`, exactly, and both inherit the 3D
/// curve's own parameter: height for the line, angle for the circle.
fn on_cylinder(
    curve: &Curve,
    range: (f64, f64),
    cylinder: ogeom_math::Cylinder,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    let axis = cylinder.axis();
    let frame = cylinder.frame();
    match curve {
        Curve::Line(line) => {
            // Parallel to the axis, on the surface.
            let direction = line.axis().direction;
            let along = direction.dot(axis.direction);
            if !direction.is_parallel(axis.direction, tol) {
                return None;
            }
            let through = line.axis().location;
            if (axis.distance_to(through) - cylinder.radius()).abs() > tol.confusion() {
                return None;
            }
            let local = frame.to_local(through);
            let u = local.y.atan2(local.x).rem_euclid(core::f64::consts::TAU);
            // The 3D line's parameter is length from its origin; at constant u
            // the pcurve's `v` runs at the same rate, signed by whether the
            // line runs with the axis or against it.
            let (lo, hi) = line.domain();
            let start = Point2::new(u, local.z);
            let towards =
                ogeom_math::Direction2::new(ogeom_math::Vector2::new(0.0, along.signum()), tol)
                    .ok()?;
            Some(
                Line2d::over(ogeom_math::Axis2::new(start, towards), lo, hi)
                    .ok()?
                    .into(),
            )
        }
        Curve::Circle(c) => {
            let circle = c.circle();
            // Perpendicular to the axis, centred on it, of the same radius.
            if circle
                .frame()
                .z()
                .cross_with(axis.direction.vector())
                .magnitude()
                > tol.angular()
            {
                return None;
            }
            if axis.distance_to(circle.centre()) > tol.confusion() {
                return None;
            }
            if (circle.radius() - cylinder.radius()).abs() > tol.confusion() {
                return None;
            }
            let local = frame.to_local(circle.centre());
            // Where the circle's own angle zero sits in the cylinder's angle,
            // and which way its parameter runs around the axis. A section
            // circle inherits its winding from the pair that made it, and one
            // wound against the cylinder's `u` (a circle cut by a plane whose
            // normal opposes the axis) runs its pcurve in `-u`. Written `+u`
            // unconditionally, the pcurve evaluates half a turn away from the
            // curve, and the face's arrangement tears along a seam that is
            // not there.
            let start = circle.centre() + circle.frame().x().vector() * circle.radius();
            let at = frame.to_local(start);
            let phase = at.y.atan2(at.x);
            let winding = circle.frame().z().vector().dot(spin(&frame)).signum();
            let towards =
                ogeom_math::Direction2::new(ogeom_math::Vector2::new(winding, 0.0), tol).ok()?;
            Some(
                Line2d::over(
                    ogeom_math::Axis2::new(Point2::new(phase, local.z), towards),
                    circle_window(range).0,
                    circle_window(range).1,
                )
                .ok()?
                .into(),
            )
        }
        Curve::Ellipse(_) => {
            // An oblique plane's section: its plan projection is the
            // cylinder's own cross-section circle traced *uniformly*, so
            // the chart trace is u = s·t + φ, v = c₀ + a·cos t + b·sin t:
            // the trig-affine family. Derived from the curve's own
            // evaluations and verified by sample, never assumed.
            use ogeom_geom::Curve3d as _;
            let tau = core::f64::consts::TAU;
            let local = |t: f64| -> Option<ogeom_math::Point> {
                Some(frame.to_local(curve.point_at(t, tol).ok()?))
            };
            let l0 = local(0.0)?;
            let lq = local(tau / 4.0)?;
            let lh = local(tau / 2.0)?;
            // On the surface at all: plan radius must be the cylinder's.
            let r = cylinder.radius();
            for l in [&l0, &lq, &lh] {
                if (l.x.hypot(l.y) - r).abs() > tol.confusion() * 10.0 {
                    return None;
                }
            }
            let phase = l0.y.atan2(l0.x);
            // Winding from the quarter-turn sample: uniform tracing puts it
            // a quarter turn away, one side or the other.
            let uq = lq.y.atan2(lq.x);
            let step = (uq - phase).rem_euclid(tau);
            let winding = if (step - tau / 4.0).abs() < 1e-6 {
                1.0
            } else if (step - 3.0 * tau / 4.0).abs() < 1e-6 {
                -1.0
            } else {
                return None;
            };
            // Height coefficients from three samples.
            let c0 = f64::midpoint(l0.z, lh.z);
            let a = (l0.z - lh.z) / 2.0;
            let b = lq.z - c0;
            // The trig formula is global (cosine wraps, the linear angle
            // unwraps the chart), so the pcurve lives on whatever range the
            // edge actually spans, a loop crossing the period included.
            let candidate = ogeom_geom::Trig2d::new(
                Point2::new(phase, c0),
                ogeom_math::Vector2::new(winding, 0.0),
                ogeom_math::Vector2::new(0.0, a),
                ogeom_math::Vector2::new(0.0, b),
                range,
            )
            .ok()?;
            // The same-parameter law, verified at points the derivation
            // never touched, inside the range the edge will use.
            use ogeom_geom::Curve2d as _;
            for i in 0..7 {
                let t = range.0 + (range.1 - range.0) * (0.09 + 0.13 * f64::from(i)) / 0.91;
                let l = local(t)?;
                let chart = candidate.point_at(t, tol).ok()?;
                let du = (chart.x - l.y.atan2(l.x)).rem_euclid(tau);
                if du.min(tau - du) > 1e-9 {
                    return None;
                }
                if (chart.y - l.z).abs() > tol.confusion() * 10.0 {
                    return None;
                }
            }
            Some(PlanarCurve::Trig(candidate))
        }
        _ => None,
    }
}

/// The pcurve of half a meridian: a great circle through both poles,
/// restricted to one side of them.
///
/// The whole circle has no chart image a single curve can carry (its
/// longitude jumps by half a turn at each pole), but each *half* does, and it
/// is a straight line. Writing the circle's own parameter as `t` and the
/// sphere's axis as `Z = cos α·X + sin α·Y` in the circle's own frame, the
/// point's height above the equator is `r·cos(t − α)`, so the latitude is
/// `asin(cos(t − α))`, which on `t − α ∈ [0, π]` is exactly `π/2 − (t − α)`,
/// affine in `t`, with slope one. The longitude is constant on that half and
/// half a turn away on the other. So the pcurve is a vertical line in the
/// chart, sharing the circle's parameter exactly, and the caller's `range` is
/// what says which half is meant.
///
/// The half is not assumed: the returned line is lifted back through the
/// sphere at stations along the range and compared against the circle, so a
/// misread orientation is caught here rather than downstream.
fn on_meridian(
    curve: &ogeom_geom::CircleCurve,
    range: (f64, f64),
    sphere: ogeom_math::Sphere,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    let circle = curve.circle();
    // A reversed circle runs its own angle backwards, and the shifted angle
    // below is measured in the *curve's* parameter, so the sign travels with
    // it: the sweep flips and so do both the latitude's slope and which half
    // of the circle a range names.
    let sweep = if curve.is_reversed() { -1.0 } else { 1.0 };
    let frame = sphere.frame();
    let z = frame.z().vector();
    // A great circle: the sphere's own centre and radius, in a plane holding
    // the axis. Anything else is not a meridian.
    if circle.centre().distance(sphere.centre()) > tol.confusion() {
        return None;
    }
    if (circle.radius() - sphere.radius()).abs() > tol.confusion() {
        return None;
    }
    let (cx, cy) = (circle.frame().x().vector(), circle.frame().y().vector());
    let (xz, yz) = (cx.dot(z), cy.dot(z));
    // The axis must lie *in* the circle's plane, or the circle is neither a
    // parallel nor a meridian and has no closed-form chart image at all.
    if xz.hypot(yz) < 1.0 - tol.angular() {
        return None;
    }
    let raw_alpha = yz.atan2(xz);
    // `w` is the circle's own horizontal direction: the axis turned a quarter
    // turn within the circle's plane.
    let w = cx * -raw_alpha.sin() + cy * raw_alpha.cos();
    let local = frame.to_local(sphere.centre() + w);
    let longitude = local.y.atan2(local.x);

    let half = core::f64::consts::PI;
    let mid = f64::midpoint(range.0, range.1);
    // Where the range sits relative to the poles, in the shifted angle
    // `x = sweep·t − α` that measures the descent from the north pole.
    let x_mid = (sweep * mid - raw_alpha).rem_euclid(core::f64::consts::TAU);
    let x_mid = if x_mid > half {
        x_mid - core::f64::consts::TAU
    } else {
        x_mid
    };
    let span = sweep * (range.1 - range.0);
    let (mut x0, mut x1) = (x_mid - span / 2.0, x_mid + span / 2.0);
    if x0 > x1 {
        core::mem::swap(&mut x0, &mut x1);
    }
    // The turn count `α` was written with is what decides whether the
    // latitude comes out inside the chart or a whole turn away from it, so
    // the branch the range actually sits on is the one the line is built
    // from.
    let alpha = sweep.mul_add(mid, -x_mid);
    // A range read off points found by projection overshoots a pole by
    // what that projection leaves: an angle the radius turns into the
    // confusion distance is no overshoot, and the lift below measures it.
    let slack = tol
        .parametric()
        .max(1e-9)
        .max(tol.confusion() / sphere.radius());
    let (axis_point, towards) = if x0 >= -slack && x1 <= half + slack {
        // The descending half: latitude π/2 − (sweep·t − α), longitude
        // constant.
        (
            Point2::new(longitude, half.mul_add(0.5, alpha)),
            ogeom_math::Vector2::new(0.0, -sweep),
        )
    } else if x0 >= -half - slack && x1 <= slack {
        // The ascending half, half a turn round the chart.
        (
            Point2::new(longitude + half, half.mul_add(0.5, -alpha)),
            ogeom_math::Vector2::new(0.0, sweep),
        )
    } else {
        // The range straddles a pole: no one line covers it.
        return None;
    };
    let towards = ogeom_math::Direction2::new(towards, tol).ok()?;
    let margin = (range.1 - range.0) * 0.25;
    let line: PlanarCurve = Line2d::over(
        ogeom_math::Axis2::new(axis_point, towards),
        range.0 - margin,
        range.1 + margin,
    )
    .ok()?
    .into();

    // Measured, not assumed: the chart line lifted back through the sphere is
    // the circle it claims to be.
    for k in 0..=4 {
        let t = (range.1 - range.0).mul_add(f64::from(k) / 4.0, range.0);
        let uv = line.point_at(t, tol).ok()?;
        let lifted = ogeom_math::elementary::sphere_at(&sphere, uv.x, uv.y).point;
        let want = curve.point_at(t, tol).ok()?;
        if lifted.distance(want) > tol.confusion() {
            return None;
        }
    }
    Some(line)
}

/// The pcurve of a circle on a sphere: a parallel of latitude, or one half of
/// a meridian.
fn on_sphere(
    curve: &Curve,
    range: (f64, f64),
    sphere: ogeom_math::Sphere,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    let Curve::Circle(c) = curve else {
        return None;
    };
    let circle = c.circle();
    let frame = sphere.frame();
    // Perpendicular to the sphere's axis and centred on it: a parallel of
    // latitude, which is a horizontal line in (longitude, latitude).
    if circle
        .frame()
        .z()
        .cross_with(frame.z().vector())
        .magnitude()
        > tol.angular()
    {
        return on_meridian(c, range, sphere, tol);
    }
    let local = frame.to_local(circle.centre());
    if local.x.abs() > tol.confusion() || local.y.abs() > tol.confusion() {
        return None;
    }
    let latitude = (local.z / sphere.radius()).clamp(-1.0, 1.0).asin();
    // Sanity: the circle's radius must be the parallel's.
    if (circle.radius() - sphere.radius() * latitude.cos()).abs() > tol.confusion() {
        return None;
    }
    let start = circle.centre() + circle.frame().x().vector() * circle.radius();
    let at = frame.to_local(start);
    let phase = at.y.atan2(at.x);
    // Phase and winding exactly as the cylinder case: a parallel whose own
    // axis opposes the sphere's marches its angle *down* the longitude.
    let winding = circle.frame().z().vector().dot(spin(&frame)).signum();
    let towards = ogeom_math::Direction2::new(ogeom_math::Vector2::new(winding, 0.0), tol).ok()?;
    Some(
        Line2d::over(
            ogeom_math::Axis2::new(Point2::new(phase, latitude), towards),
            circle_window(range).0,
            circle_window(range).1,
        )
        .ok()?
        .into(),
    )
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_geom::{Curve2d, Curve3d, CylinderSurface, PlaneSurface, SphereSurface};
    use ogeom_math::{Cylinder, Direction, Frame, Plane, Sphere, Vector};

    const T: Tolerances = Tolerances::millimetres();

    fn sphere(centre: Point, radius: f64) -> SurfaceGeometry {
        SphereSurface::new(Sphere::centred(centre, radius, T).unwrap()).into()
    }

    fn cylinder(axis: Vector, radius: f64) -> SurfaceGeometry {
        let frame = Frame::new(
            Point::ORIGIN,
            Direction::new(axis, T).unwrap(),
            Direction::from_cross(axis, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
            T,
        )
        .unwrap();
        CylinderSurface::new(Cylinder::new(frame, radius, T).unwrap(), (-4.0, 4.0))
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

    /// Same-parameter: pcurve lifted through its surface equals the 3D curve,
    /// at the same parameter, everywhere sampled.
    fn assert_same_parameter(
        section: &SectionCurve,
        surface: &SurfaceGeometry,
        pcurve: &PlanarCurve,
        samples: usize,
    ) {
        let (lo, hi) = section.curve.domain();
        let (plo, phi) = pcurve.domain();
        assert!(
            (lo - plo).abs() < 1e-9 && (hi - phi).abs() < 1e-9,
            "domains disagree: [{lo}, {hi}] against [{plo}, {phi}]"
        );
        for i in 0..=samples {
            #[allow(clippy::cast_precision_loss)]
            let t = lo + (hi - lo) * i as f64 / samples as f64;
            let on_curve = section.curve.point_at(t, T).unwrap();
            let at = pcurve.point_at(t, T).unwrap();
            let lifted = surface.point_at(at.x, at.y, T).unwrap();
            assert!(
                on_curve.is_equal(lifted, T),
                "at t = {t}: curve {on_curve:?}, lifted {lifted:?}"
            );
        }
    }

    /// A straight side from `(√2, 0)` to `(0, √2)` as a clamped cubic
    /// B-spline, its control points even along it, extruded 3 along `z`
    /// from `z = 0` and over `z` in `[1, 4]`; and the same side as a
    /// cubic-by-linear spline patch from `z = 0`.
    fn straight_spline_sides() -> [SurfaceGeometry; 3] {
        let (from, to) = (
            Point::new(2.0_f64.sqrt(), 0.0, 0.0),
            Point::new(0.0, 2.0_f64.sqrt(), 0.0),
        );
        let along = |k: u32| from + (to - from) * (f64::from(k) / 3.0);
        let knots = || ogeom_math::KnotVector::clamped_uniform(3, 4).unwrap();
        let side = ogeom_geom::BSplineCurve::new(knots(), (0..4).map(along).collect(), T).unwrap();
        let lowered =
            ogeom_geom::ExtrusionSurface::over(side.clone().into(), Direction::Z, (1.0, 4.0))
                .unwrap();
        let extruded = ogeom_geom::ExtrusionSurface::new(side.into(), Direction::Z, 3.0).unwrap();
        let rows: Vec<Point> = (0..4)
            .flat_map(|k| [along(k), along(k) + Vector::new(0.0, 0.0, 3.0)])
            .collect();
        let patch = ogeom_geom::BSplineSurface::new(
            knots(),
            ogeom_math::KnotVector::clamped_uniform(1, 2).unwrap(),
            &ogeom_math::ControlGrid::new(rows, 4, 2).unwrap(),
            T,
        )
        .unwrap();
        [extruded.into(), lowered.into(), patch.into()]
    }

    /// A plane holding the extrusion direction through either end of a
    /// spline side meets the side in the straight ruling at that border,
    /// exactly, as a plane through a planar side's edge does; through the
    /// side's inside, in the ruling there. A border the surface leaves
    /// tangentially is not stated as a crossing.
    #[test]
    fn a_plane_through_a_spline_patch_s_border_ruling_meets_it_there() {
        let r = 2.0_f64.sqrt();
        for side in straight_spline_sides() {
            for (through, normal, u) in [
                (Point::ORIGIN, Vector::Y, 0.0),
                (Point::ORIGIN, Vector::X, 1.0),
                (Point::new(0.5, 0.0, 0.0), Vector::X, 1.0 - 0.5 / r),
            ] {
                let cut = plane(through, normal);
                for (a, b, cut_first) in [(&cut, &side, true), (&side, &cut, false)] {
                    let met = intersect_surfaces(a, b, IntersectOptions::default(), T).unwrap();
                    let SurfaceIntersection::Along(sections) = met else {
                        panic!("the plane through {through:?} holds a ruling: {met:?}");
                    };
                    assert_eq!(sections.len(), 1, "{through:?}, {normal:?}");
                    let section = &sections[0];
                    assert!(section.exact && !section.tangential, "{section:?}");
                    let (_, (v0, v1)) = side.domain();
                    let want = [
                        side.point_at(u, v0, T).unwrap(),
                        side.point_at(u, v1, T).unwrap(),
                    ];
                    let (lo, hi) = section.curve.domain();
                    let ends = [
                        section.curve.point_at(lo, T).unwrap(),
                        section.curve.point_at(hi, T).unwrap(),
                    ];
                    assert!(
                        ends[0].is_equal(want[0], T) && ends[1].is_equal(want[1], T),
                        "{through:?}: the ruling at u = {u}, {want:?}, got {ends:?}"
                    );
                    let (on_side, on_cut) = if cut_first {
                        (&section.on_b, &section.on_a)
                    } else {
                        (&section.on_a, &section.on_b)
                    };
                    assert_same_parameter(section, &side, on_side.as_ref().unwrap(), 16);
                    assert_same_parameter(section, &cut, on_cut.as_ref().unwrap(), 16);
                }
            }
        }

        // A side leaving `(0, 0)` along `x` and bending away: the plane
        // `y = 0` holds its first ruling, tangent to it there.
        let bend = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.5), (2.0, 2.0)]
            .map(|(x, y)| Point::new(x, y, 0.0))
            .to_vec();
        let knots = ogeom_math::KnotVector::clamped_uniform(3, 4).unwrap();
        let side = ogeom_geom::BSplineCurve::new(knots, bend, T).unwrap();
        let bent: SurfaceGeometry =
            ogeom_geom::ExtrusionSurface::new(side.into(), Direction::Z, 3.0)
                .unwrap()
                .into();
        let cut = plane(Point::ORIGIN, Vector::Y);
        if let SurfaceIntersection::Along(sections) =
            intersect_surfaces(&cut, &bent, IntersectOptions::default(), T).unwrap()
        {
            assert!(
                sections.iter().all(|s| !s.exact || s.tangential),
                "a tangent border is no crossing: {sections:?}"
            );
        }
    }

    /// A side running from the corner `(√2, 0)` toward `(0, √2)`, stood
    /// up three tall, in each spelling no closed form or iso-line path
    /// reads: an extrusion of half an ellipse bulging outward, trimmed; of
    /// a B-spline line, trimmed; of a B-spline line offset onto the side;
    /// a spline patch; an extrusion of a B-spline line. Each starts on the
    /// ruling through the corner, a border of its chart, and leaves it
    /// along `(1, 1)` or `(-1, 1)`, into `y > 0`.
    fn sides_from_the_corner() -> Vec<SurfaceGeometry> {
        use ogeom_geom::{BSplineCurve, EllipseCurve, ExtrusionSurface, OffsetCurve, TrimmedCurve};
        let r = 2.0_f64.sqrt();
        let (from, to) = (Point::new(r, 0.0, 0.0), Point::new(0.0, r, 0.0));
        let line = |from: Point, to: Point| -> Curve {
            let along = |k: u32| from + (to - from) * (f64::from(k) / 3.0);
            let knots = ogeom_math::KnotVector::clamped_uniform(3, 4).unwrap();
            BSplineCurve::new(knots, (0..4).map(along).collect(), T)
                .unwrap()
                .into()
        };
        let middle = from.midpoint(to);
        let frame = Frame::new(
            middle,
            Direction::Z,
            Direction::new(from - to, T).unwrap(),
            T,
        )
        .unwrap();
        let half_ellipse = TrimmedCurve::new(
            EllipseCurve::new(ogeom_math::Ellipse::new(frame, 1.0, 0.4, T).unwrap()).into(),
            0.0,
            core::f64::consts::PI,
            T,
        )
        .unwrap();
        let run = to - from;
        let trimmed_line =
            TrimmedCurve::new(line(from - run * 0.5, to + run * 0.5), 0.25, 0.75, T).unwrap();
        let inset = Vector::new(1.0, 1.0, 0.0) * (0.25 / r);
        let offset = OffsetCurve::new(line(from - inset, to - inset), 0.25, Direction::Z).unwrap();
        let extruded = |curve: Curve| -> SurfaceGeometry {
            ExtrusionSurface::new(curve, Direction::Z, 3.0)
                .unwrap()
                .into()
        };
        let [_, _, patch] = straight_spline_sides();
        vec![
            extruded(half_ellipse.into()),
            extruded(trimmed_line.into()),
            extruded(Curve::Offset(Box::new(offset))),
            patch,
            extruded(line(from, to)),
        ]
    }

    /// Where a border of one surface's chart lies in the other surface
    /// and the two leave each other there at an angle, they meet along
    /// it, whichever surface comes first: here the ruling `x = √2, y = 0`
    /// of each side, met by the plane `y = 0` and by the drum of radius
    /// one about `(√2, -1)`, which rests on the ruling from below. Every
    /// side lies in `y >= 0`, so the ruling is all they share.
    #[test]
    fn a_border_lying_in_the_other_surface_is_a_section_along_it() {
        let r = 2.0_f64.sqrt();
        let drum: SurfaceGeometry = CylinderSurface::new(
            Cylinder::new(
                Frame::new(Point::new(r, -1.0, -1.0), Direction::Z, Direction::X, T).unwrap(),
                1.0,
                T,
            )
            .unwrap(),
            (0.0, 5.0),
        )
        .unwrap()
        .into();
        let cutters = [
            ("the plane y = 0", plane(Point::ORIGIN, Vector::Y)),
            ("the drum", drum),
        ];
        for (k, side) in sides_from_the_corner().iter().enumerate() {
            for (name, cutter) in &cutters {
                for (a, b) in [(cutter, side), (side, cutter)] {
                    let met = intersect_surfaces(a, b, IntersectOptions::default(), T).unwrap();
                    let SurfaceIntersection::Along(sections) = met else {
                        panic!("side {k}, {name}: the ruling is shared, got {met:?}");
                    };
                    assert_eq!(sections.len(), 1, "side {k}, {name}: {sections:?}");
                    let section = &sections[0];
                    assert!(!section.tangential, "side {k}, {name}: {section:?}");
                    // Each point of the section within its stated tolerance
                    // of the ruling, the ends on its ends.
                    let (lo, hi) = section.curve.domain();
                    let near = section.tolerance + T.confusion();
                    for i in 0..=16 {
                        let p = section
                            .curve
                            .point_at(lo + (hi - lo) * f64::from(i) / 16.0, T)
                            .unwrap();
                        let off = (p.x - r).hypot(p.y);
                        assert!(
                            off <= near,
                            "side {k}, {name}: {p:?} is {off} off the ruling"
                        );
                    }
                    let ends = [
                        section.curve.point_at(lo, T).unwrap(),
                        section.curve.point_at(hi, T).unwrap(),
                    ];
                    let mut heights = [ends[0].z, ends[1].z];
                    heights.sort_by(f64::total_cmp);
                    assert!(
                        heights[0].abs() <= near && (heights[1] - 3.0).abs() <= near,
                        "side {k}, {name}: the section runs {heights:?}"
                    );
                }
            }
        }

        // The drum of radius one about `(√2 - 1/√2, -1/√2)` touches the
        // spline patch's plane along the ruling, from outside: the ruling
        // lies in it, and the two leave each other tangentially there.
        // That is a contact, not a crossing.
        let [_, _, patch] = straight_spline_sides();
        let resting: SurfaceGeometry = CylinderSurface::new(
            Cylinder::new(
                Frame::new(
                    Point::new(r - 1.0 / r, -1.0 / r, -1.0),
                    Direction::Z,
                    Direction::X,
                    T,
                )
                .unwrap(),
                1.0,
                T,
            )
            .unwrap(),
            (0.0, 5.0),
        )
        .unwrap()
        .into();
        for (a, b) in [(&resting, &patch), (&patch, &resting)] {
            if let SurfaceIntersection::Along(sections) =
                intersect_surfaces(a, b, IntersectOptions::default(), T).unwrap()
            {
                assert!(
                    sections.iter().all(|s| s.tangential),
                    "a border left tangentially is no crossing: {sections:?}"
                );
            }
        }
    }

    /// The patch over `0 ≤ x ≤ 1`, `0 ≤ y ≤ 12` standing `z = s y + q y²`
    /// off the plane `z = 0`, linear along `x`, quadratic along `y` with
    /// the middle control point on the start's heading half way: its border
    /// `y = 0` lies in the plane, which it leaves there at the sine
    /// `s / √(1 + s²)`.
    fn wall_leaving(s: f64, q: f64) -> SurfaceGeometry {
        let points: Vec<Point> = [0.0, 1.0]
            .into_iter()
            .flat_map(|x| {
                [(0.0, 0.0), (6.0, 6.0 * s), (12.0, 12.0 * s + 144.0 * q)]
                    .map(|(y, z)| Point::new(x, y, z))
            })
            .collect();
        ogeom_geom::BSplineSurface::new(
            ogeom_math::KnotVector::clamped_uniform(1, 2).unwrap(),
            ogeom_math::KnotVector::clamped_uniform(2, 3).unwrap(),
            &ogeom_math::ControlGrid::new(points, 2, 3).unwrap(),
            T,
        )
        .unwrap()
        .into()
    }

    /// A window of the plane `z = 0` reaching a unit past the wall's border
    /// `y = 0` and no further than `y = 1`, short of the wall's next row
    /// of samples at `y = 2`, meets the wall along that border however
    /// shallowly the wall leaves it: the samples inside the window all lie
    /// on the border, which says only that the border lies in the plane,
    /// and the wall's rise `s d` a distance `d` in is resolved within a
    /// fraction of a millimetre. A wall leaving the plane tangentially there
    /// is a contact, never a crossing.
    #[test]
    fn a_wall_leaving_a_plane_window_shallowly_meets_it_along_its_border() {
        let window: SurfaceGeometry =
            PlaneSurface::over(Plane::new(Frame::WORLD), (-1.0, 2.0), (-1.0, 1.0))
                .unwrap()
                .into();
        for s in [0.005, 0.012, 0.023] {
            let wall = wall_leaving(s, 0.002);
            for (a, b) in [(&window, &wall), (&wall, &window)] {
                let met = intersect_surfaces(a, b, IntersectOptions::default(), T).unwrap();
                let SurfaceIntersection::Along(sections) = met else {
                    panic!("sine {s}: the border is shared, got {met:?}");
                };
                assert_eq!(sections.len(), 1, "sine {s}: {sections:?}");
                let section = &sections[0];
                assert!(!section.tangential, "sine {s}: {section:?}");
                let (lo, hi) = section.curve.domain();
                let near = section.tolerance + T.confusion();
                for i in 0..=16 {
                    let p = section
                        .curve
                        .point_at(lo + (hi - lo) * f64::from(i) / 16.0, T)
                        .unwrap();
                    assert!(
                        p.y.hypot(p.z) <= near && (-near..=1.0 + near).contains(&p.x),
                        "sine {s}: {p:?} is off the border"
                    );
                }
                let mut xs = [lo, hi].map(|t| section.curve.point_at(t, T).unwrap().x);
                xs.sort_by(f64::total_cmp);
                assert!(
                    xs[0].abs() <= near && (xs[1] - 1.0).abs() <= near,
                    "sine {s}: the section runs {xs:?}"
                );
            }
        }
        let resting = wall_leaving(0.0, 0.002);
        for (a, b) in [(&window, &resting), (&resting, &window)] {
            if let SurfaceIntersection::Along(sections) =
                intersect_surfaces(a, b, IntersectOptions::default(), T).unwrap()
            {
                assert!(
                    sections.iter().all(|s| s.tangential),
                    "a border left tangentially is no crossing: {sections:?}"
                );
            }
        }
    }

    #[test]
    fn an_analytic_pair_comes_back_exact_with_matching_pcurves() {
        // A plane through a cylinder's axis: two lines, and every description
        // agrees at the same parameter, which is the claim edges carry and
        // booleans rely on.
        let drum = cylinder(Vector::Z, 2.0);
        let cut = plane(Point::ORIGIN, Vector::X);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&drum, &cut, IntersectOptions::default(), T).unwrap()
        else {
            panic!("a plane through a cylinder meets it along curves");
        };
        assert_eq!(curves.len(), 2);
        for section in &curves {
            assert!(section.exact);
            assert!((section.tolerance - 0.0).abs() < f64::EPSILON);
            let on_a = section.on_a.as_ref().expect("a line has a cylinder pcurve");
            let on_b = section.on_b.as_ref().expect("and a plane pcurve");
            assert_same_parameter(section, &drum, on_a, 50);
            assert_same_parameter(section, &cut, on_b, 50);
        }
    }

    #[test]
    fn an_oblique_cut_gives_the_ellipse_a_trig_pcurve_on_the_drum() {
        // The oblique ellipse's pcurve runs linearly in the chart angle and
        // sinusoidally in height (the trig-affine family), exactly,
        // same-parameter, both sides.
        let drum = cylinder(Vector::Z, 2.0);
        let angle: f64 = 0.5;
        let cut = plane(Point::ORIGIN, Vector::new(0.0, angle.sin(), angle.cos()));
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&drum, &cut, IntersectOptions::default(), T).unwrap()
        else {
            panic!("an oblique plane meets the cylinder along its ellipse");
        };
        assert_eq!(curves.len(), 1);
        let section = &curves[0];
        assert!(section.exact);
        assert!(matches!(section.curve, Curve::Ellipse(_)));
        let on_drum = section
            .on_a
            .as_ref()
            .expect("the oblique ellipse now carries its cylinder pcurve");
        assert!(
            matches!(on_drum, PlanarCurve::Trig(_)),
            "the chart trace is trig-affine: {on_drum:?}"
        );
        assert_same_parameter(section, &drum, on_drum, 60);
        let on_plane = section.on_b.as_ref().expect("and its plane pcurve");
        assert_same_parameter(section, &cut, on_plane, 60);
    }

    #[test]
    fn a_perpendicular_cut_gives_a_circle_with_a_straight_pcurve() {
        let drum = cylinder(Vector::Z, 2.0);
        let cut = plane(Point::new(0.0, 0.0, 1.0), Vector::Z);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&drum, &cut, IntersectOptions::default(), T).unwrap()
        else {
            panic!("expected curves");
        };
        assert_eq!(curves.len(), 1);
        let section = &curves[0];
        assert!(section.closed);
        assert!(matches!(section.curve, Curve::Circle(_)));
        // On the cylinder the circle is a horizontal line in (u, v).
        assert!(matches!(
            section.on_a.as_ref().unwrap(),
            PlanarCurve::Line(_)
        ));
        assert_same_parameter(section, &drum, section.on_a.as_ref().unwrap(), 60);
        assert_same_parameter(section, &cut, section.on_b.as_ref().unwrap(), 60);
    }

    #[test]
    fn coaxial_cylinder_and_sphere_give_circles_with_pcurves_on_both() {
        let drum = cylinder(Vector::Z, 1.5);
        let ball = sphere(Point::ORIGIN, 3.0);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&drum, &ball, IntersectOptions::default(), T).unwrap()
        else {
            panic!("expected curves");
        };
        assert_eq!(curves.len(), 2);
        for section in &curves {
            assert!(section.exact);
            assert_same_parameter(section, &drum, section.on_a.as_ref().unwrap(), 40);
            assert_same_parameter(section, &ball, section.on_b.as_ref().unwrap(), 40);
        }
    }

    fn torus(origin: Point, axis: Vector, major: f64, minor: f64) -> SurfaceGeometry {
        let frame = Frame::new(
            origin,
            Direction::new(axis, T).unwrap(),
            Direction::from_cross(axis, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
            T,
        )
        .unwrap();
        ogeom_geom::TorusSurface::new(ogeom_math::Torus::new(frame, major, minor, T).unwrap())
            .into()
    }

    #[test]
    fn an_axis_normal_plane_meets_a_torus_in_two_parallels_with_pcurves() {
        let ring = torus(Point::ORIGIN, Vector::Z, 2.0, 0.5);
        let cut = plane(Point::new(0.0, 0.0, 0.3), Vector::Z);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&ring, &cut, IntersectOptions::default(), T).unwrap()
        else {
            panic!("an axis-normal plane through the tube meets it along curves");
        };
        assert_eq!(curves.len(), 2);
        let spread = 0.5_f64.mul_add(0.5, -(0.3 * 0.3)).sqrt();
        let mut radii: Vec<f64> = curves
            .iter()
            .map(|s| {
                let Curve::Circle(c) = &s.curve else {
                    panic!("a parallel is a circle");
                };
                c.circle().radius()
            })
            .collect();
        radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((radii[0] - (2.0 - spread)).abs() < 1e-12);
        assert!((radii[1] - (2.0 + spread)).abs() < 1e-12);
        for section in &curves {
            assert!(section.exact);
            assert_same_parameter(section, &ring, section.on_a.as_ref().unwrap(), 48);
            assert_same_parameter(section, &cut, section.on_b.as_ref().unwrap(), 48);
        }
    }

    #[test]
    fn the_plane_a_ball_rolls_on_touches_its_torus_along_the_circle_it_rolled() {
        // Tangency with length is reported as the curve it is (the way a
        // tangent plane reports its line on a cylinder), because the blend
        // machinery builds faces whose boundaries are exactly these circles,
        // and a Touching with no curve in it would read as a refusal upstream.
        let ring = torus(Point::ORIGIN, Vector::Z, 2.0, 0.5);
        let cut = plane(Point::new(0.0, 0.0, 0.5), Vector::Z);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&ring, &cut, IntersectOptions::default(), T).unwrap()
        else {
            panic!("the rolling plane touches along a circle, not at points");
        };
        assert_eq!(curves.len(), 1);
        let Curve::Circle(c) = &curves[0].curve else {
            panic!("the tangency is a circle");
        };
        assert!((c.circle().radius() - 2.0).abs() < 1e-12);
        assert_same_parameter(&curves[0], &ring, curves[0].on_a.as_ref().unwrap(), 48);
        assert_same_parameter(&curves[0], &cut, curves[0].on_b.as_ref().unwrap(), 48);
    }

    #[test]
    fn a_coaxial_cylinder_meets_a_torus_in_two_parallels_and_touches_in_one() {
        let ring = torus(Point::ORIGIN, Vector::Z, 2.0, 0.5);
        let drum = cylinder(Vector::Z, 2.2);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&drum, &ring, IntersectOptions::default(), T).unwrap()
        else {
            panic!("a coaxial cylinder through the tube meets it along curves");
        };
        assert_eq!(curves.len(), 2);
        for section in &curves {
            assert!(section.exact);
            let Curve::Circle(c) = &section.curve else {
                panic!("a parallel is a circle");
            };
            assert!((c.circle().radius() - 2.2).abs() < 1e-12);
            assert_same_parameter(section, &drum, section.on_a.as_ref().unwrap(), 48);
            assert_same_parameter(section, &ring, section.on_b.as_ref().unwrap(), 48);
        }

        // Tangent at the tube's outer equator: one circle, with both pcurves.
        let grazing = cylinder(Vector::Z, 2.5);
        let SurfaceIntersection::Along(touch) =
            intersect_surfaces(&grazing, &ring, IntersectOptions::default(), T).unwrap()
        else {
            panic!("the grazing cylinder touches along the equator");
        };
        assert_eq!(touch.len(), 1);
        assert_same_parameter(&touch[0], &grazing, touch[0].on_a.as_ref().unwrap(), 48);
        assert_same_parameter(&touch[0], &ring, touch[0].on_b.as_ref().unwrap(), 48);
    }

    #[test]
    fn coaxial_tori_are_the_same_or_meet_in_parallels() {
        let ring = torus(Point::ORIGIN, Vector::Z, 2.0, 0.5);
        assert!(matches!(
            intersect_surfaces(&ring, &ring.clone(), IntersectOptions::default(), T).unwrap(),
            SurfaceIntersection::Same
        ));

        // The same tube lifted half a radius: the profile circles cross
        // twice, and each crossing revolves into a parallel shared exactly.
        let lifted = torus(Point::new(0.0, 0.0, 0.5), Vector::Z, 2.0, 0.5);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&ring, &lifted, IntersectOptions::default(), T).unwrap()
        else {
            panic!("lifted coaxial tori meet along curves");
        };
        assert_eq!(curves.len(), 2);
        for section in &curves {
            assert!(section.exact);
            assert_same_parameter(section, &ring, section.on_a.as_ref().unwrap(), 48);
            assert_same_parameter(section, &lifted, section.on_b.as_ref().unwrap(), 48);
        }
    }

    #[test]
    fn a_pair_with_no_closed_form_comes_back_fitted_with_pcurves() {
        // Crossed cylinders: the marched path, end to end through one call.
        let a = cylinder(Vector::Z, 1.0);
        let b = cylinder(Vector::X, 1.6);
        let options = IntersectOptions {
            tolerance: 1e-5,
            marching: Marching {
                chord: 1e-5,
                ..Marching::default()
            },
        };
        let SurfaceIntersection::Along(curves) = intersect_surfaces(&a, &b, options, T).unwrap()
        else {
            panic!("crossed cylinders meet along curves");
        };
        assert_eq!(curves.len(), 2);
        for section in &curves {
            assert!(!section.exact);
            assert!(section.closed);
            assert!(
                section.tolerance <= 1e-5 + 1e-4,
                "got {}",
                section.tolerance
            );
            assert!(section.on_a.is_some() && section.on_b.is_some());

            // The fitted curve lies on both cylinders to its stated tolerance.
            let (lo, hi) = section.curve.domain();
            for i in 0..=200 {
                #[allow(clippy::cast_precision_loss)]
                let t = lo + (hi - lo) * f64::from(i) / 200.0;
                let p = section.curve.point_at(t, T).unwrap();
                let (SurfaceGeometry::Cylinder(x), SurfaceGeometry::Cylinder(y)) = (&a, &b) else {
                    unreachable!()
                };
                let off = x
                    .cylinder()
                    .distance_to(p)
                    .abs()
                    .max(y.cylinder().distance_to(p).abs());
                assert!(
                    off <= section.tolerance * 2.0,
                    "at t = {t} the fitted curve is {off:e} off, tolerance {}",
                    section.tolerance
                );
            }
        }
    }

    /// A plane all but parallel to a drum's axis meets it in an ellipse ten
    /// metres long, which crosses the drum's few units of height only in a
    /// sliver of its turn. It is still a section of the two.
    #[test]
    fn a_plane_all_but_along_the_axis_still_meets_a_short_drum() {
        let drum = cylinder(Vector::Z, 1.0);
        let wall: SurfaceGeometry = PlaneSurface::over(
            Plane::through(
                Point::new(0.0, 0.6, 0.0),
                Direction::new(Vector::new(0.0, 1.0, 1e-4), T).unwrap(),
            ),
            (-1e9, 1e9),
            (-1e9, 1e9),
        )
        .unwrap()
        .into();
        let met = intersect_surfaces(&wall, &drum, IntersectOptions::default(), T).unwrap();
        let SurfaceIntersection::Along(sections) = met else {
            panic!("the wall crosses the drum: {met:?}");
        };
        assert_eq!(sections.len(), 1);
        let curve = &sections[0].curve;
        let (lo, hi) = curve.domain();
        let inside = (0..=100_000).any(|k| {
            let p = curve
                .point_at(lo + (hi - lo) * f64::from(k) / 100_000.0, T)
                .unwrap();
            p.z.abs() <= 4.0
        });
        assert!(inside, "and the section runs through the drum's height");
    }

    /// Every point of a section within its stated tolerance of both
    /// surfaces, sampled along it.
    fn on_both(section: &SectionCurve, a: &SurfaceGeometry, b: &SurfaceGeometry) {
        let (lo, hi) = section.curve.domain();
        for k in 0..=64 {
            let p = section
                .curve
                .point_at(lo + (hi - lo) * f64::from(k) / 64.0, T)
                .unwrap();
            for surface in [a, b] {
                let off = match surface {
                    SurfaceGeometry::Plane(plane) => plane.plane().signed_distance_to(p).abs(),
                    SurfaceGeometry::Cylinder(drum) => {
                        let axis = drum.cylinder().axis();
                        let rel = p - axis.location;
                        let d = axis.direction.vector();
                        ((rel - d * rel.dot(d)).magnitude() - drum.cylinder().radius()).abs()
                    }
                    _ => unreachable!("planes and drums only"),
                };
                assert!(
                    off <= section.tolerance + 1e-9,
                    "{p:?} is {off:e} off, stated {:e}",
                    section.tolerance
                );
            }
        }
    }

    /// A plane leaning two hundred-thousandths off a drum's axis, grazing
    /// it: the closed form's ellipse is fifty metres long, its parameter
    /// too coarse for the drum's eight units of height. The two sections
    /// come back as curves along that height, within their stated
    /// tolerance of both surfaces.
    #[test]
    fn a_plane_all_but_along_a_drums_axis_meets_it_in_two_near_lines() {
        let drum = cylinder(Vector::Z, 1.0);
        let wall: SurfaceGeometry = PlaneSurface::over(
            Plane::through(
                Point::new(0.0, 0.99, 0.0),
                Direction::new(Vector::new(0.0, 1.0, 2e-5), T).unwrap(),
            ),
            (-1e9, 1e9),
            (-1e9, 1e9),
        )
        .unwrap()
        .into();
        let met = intersect_surfaces(&wall, &drum, IntersectOptions::default(), T).unwrap();
        let SurfaceIntersection::Along(sections) = met else {
            panic!("the wall crosses the drum: {met:?}");
        };
        assert_eq!(sections.len(), 2);
        for section in &sections {
            assert!(section.tolerance > 0.0 && section.tolerance <= 1e-5);
            on_both(section, &wall, &drum);
        }
    }

    /// Two drums whose axes lean five hundred-thousandths apart meet in two
    /// curves all but straight, returned as such over the height they share
    /// rather than marched.
    #[test]
    fn drums_all_but_parallel_meet_in_two_near_lines() {
        let drill = cylinder(Vector::Z, 1.0);
        let frame = Frame::new(
            Point::new(1.5, 0.0, 0.0),
            Direction::new(Vector::new(5e-5, 0.0, 1.0), T).unwrap(),
            Direction::X,
            T,
        )
        .unwrap();
        let bore: SurfaceGeometry =
            CylinderSurface::new(Cylinder::new(frame, 1.0, T).unwrap(), (-3.0, 3.0))
                .unwrap()
                .into();
        let met = intersect_surfaces(&drill, &bore, IntersectOptions::default(), T).unwrap();
        let SurfaceIntersection::Along(sections) = met else {
            panic!("the drums cross: {met:?}");
        };
        assert_eq!(sections.len(), 2);
        for section in &sections {
            assert!(!section.exact && section.tolerance <= 1e-5);
            let (lo, hi) = section.curve.domain();
            let (p, q) = (
                section.curve.point_at(lo, T).unwrap(),
                section.curve.point_at(hi, T).unwrap(),
            );
            assert!(
                (p.z - q.z).abs() > 5.9,
                "over the shared height: {p:?} {q:?}"
            );
            on_both(section, &drill, &bore);
        }
    }

    #[test]
    fn exact_lines_are_clipped_to_the_surfaces_extents() {
        // The analytic layer answers for the unbounded geometry; the surfaces
        // are finite. A section line a billion units long is not something an
        // edge can be built on, and one wholly outside the extents is a
        // phantom.
        let drum = cylinder(Vector::Z, 2.0);
        let cut = plane(Point::ORIGIN, Vector::X);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&drum, &cut, IntersectOptions::default(), T).unwrap()
        else {
            panic!("expected curves");
        };
        for section in &curves {
            let (lo, hi) = section.curve.domain();
            // Bounded by the cylinder's height, not by LINE_EXTENT.
            assert!(
                hi - lo <= 8.0 + 1e-9,
                "the line was not clipped: [{lo}, {hi}]"
            );
            let start = section.curve.point_at(lo, T).unwrap();
            let end = section.curve.point_at(hi, T).unwrap();
            assert!(start.z >= -4.0 - 1e-9 && end.z <= 4.0 + 1e-9);
        }

        // A circle at a height the bounded cylinder does not reach is not an
        // intersection of these surfaces, however truly the unbounded ones
        // meet there.
        let high = plane(Point::new(0.0, 0.0, 10.0), Vector::Z);
        assert_eq!(
            intersect_surfaces(&drum, &high, IntersectOptions::default(), T).unwrap(),
            SurfaceIntersection::Apart
        );
    }

    #[test]
    fn the_degenerate_answers_pass_through() {
        assert_eq!(
            intersect_surfaces(
                &sphere(Point::ORIGIN, 1.0),
                &sphere(Point::new(5.0, 0.0, 0.0), 1.0),
                IntersectOptions::default(),
                T
            )
            .unwrap(),
            SurfaceIntersection::Apart
        );
        assert_eq!(
            intersect_surfaces(
                &sphere(Point::ORIGIN, 1.0),
                &sphere(Point::ORIGIN, 1.0),
                IntersectOptions::default(),
                T
            )
            .unwrap(),
            SurfaceIntersection::Same
        );
        assert!(matches!(
            intersect_surfaces(
                &plane(Point::ORIGIN, Vector::Z),
                &sphere(Point::new(0.0, 0.0, 2.0), 2.0),
                IntersectOptions::default(),
                T
            )
            .unwrap(),
            SurfaceIntersection::Touching(ref p) if p.len() == 1
        ));
    }

    #[test]
    fn unusable_options_are_refused() {
        let a = sphere(Point::ORIGIN, 1.0);
        let b = plane(Point::ORIGIN, Vector::Z);
        for tolerance in [0.0, -1.0, f64::NAN] {
            let options = IntersectOptions {
                tolerance,
                ..IntersectOptions::default()
            };
            assert!(intersect_surfaces(&a, &b, options, T).is_err());
        }
    }

    #[test]
    fn a_circle_wound_against_the_axis_keeps_its_pcurve_same_parameter() {
        // A plane whose normal opposes the cylinder's axis cuts a circle
        // wound against the cylinder's `u`, and the pcurve must run in `-u`
        // with it. Written `+u` unconditionally, the pcurve evaluates half a
        // turn away from the curve and every face built on the section tears
        // in parameter space. Both windings are pinned by lifting the pcurve
        // through the surface and demanding the curve's own point back.
        let drum: SurfaceGeometry = CylinderSurface::new(
            Cylinder::new(
                Frame::new(Point::new(2.0, 2.0, -1.0), Direction::Z, Direction::X, T).unwrap(),
                0.5,
                T,
            )
            .unwrap(),
            (0.0, 3.0),
        )
        .unwrap()
        .into();
        for normal in [Direction::Z, -Direction::Z] {
            let frame = Frame::new(Point::ORIGIN, normal, Direction::X, T).unwrap();
            let ground: SurfaceGeometry =
                PlaneSurface::over(Plane::new(frame), (-4.0, 4.0), (-4.0, 4.0))
                    .unwrap()
                    .into();
            let met = intersect_surfaces(&ground, &drum, IntersectOptions::default(), T).unwrap();
            let SurfaceIntersection::Along(curves) = met else {
                panic!("a plane through a cylinder sections it");
            };
            for sc in &curves {
                let pcurve = sc
                    .on_b
                    .as_ref()
                    .expect("a circle on its cylinder has a pcurve");
                let (lo, hi) = sc.curve.domain();
                for i in 0..8 {
                    let t = lo + (hi - lo) * f64::from(i) / 8.0;
                    let p3 = sc.curve.point_at(t, T).unwrap();
                    let uv = pcurve.point_at(t, T).unwrap();
                    let lifted = drum
                        .point_at(uv.x.rem_euclid(core::f64::consts::TAU), uv.y, T)
                        .unwrap();
                    assert!(
                        p3.distance(lifted) < 1e-9,
                        "normal {normal:?}, t {t}: pcurve lifts {lifted:?} against {p3:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn a_circle_on_a_left_handed_frame_keeps_its_pcurve_same_parameter() {
        // A mirrored circle stands on a left-handed frame: its parameter
        // runs from `x` towards a `y` that turns clockwise about `z`. Read
        // by `z` alone, its pcurve runs the other way round the axis from
        // the curve, and lifts to the circle's mirror image across its `x`.
        // Each circle is taken on both frames through its `x`, and its
        // pcurve lifted through the surface must give the curve's own point
        // back at every parameter.
        use ogeom_geom::ConeSurface;
        let cone: SurfaceGeometry = ConeSurface::new(
            ogeom_math::Cone::new(Frame::WORLD, 1.0, core::f64::consts::FRAC_PI_6, T).unwrap(),
            (-1.0, 3.0),
        )
        .unwrap()
        .into();
        let on = |centre: Point, normal: Direction, radius: f64| {
            (Frame::new(centre, normal, Direction::X, T).unwrap(), radius)
        };
        let cases = [
            (
                cylinder(Vector::Z, 2.0),
                on(Point::new(0.0, 0.0, 1.0), Direction::Z, 2.0),
            ),
            (
                sphere(Point::ORIGIN, 2.0),
                on(Point::new(0.0, 0.0, 1.0), Direction::Z, 3.0_f64.sqrt()),
            ),
            (
                cone.clone(),
                on(
                    Point::new(0.0, 0.0, 1.0),
                    Direction::Z,
                    1.0 + core::f64::consts::FRAC_PI_6.tan(),
                ),
            ),
            (
                torus(Point::ORIGIN, Vector::Z, 2.0, 0.5),
                on(Point::ORIGIN, Direction::Z, 2.5),
            ),
            (
                torus(Point::ORIGIN, Vector::Z, 2.0, 0.5),
                on(Point::new(2.0, 0.0, 0.0), Direction::Y, 0.5),
            ),
        ];
        for (surface, (frame, radius)) in &cases {
            for frame in [*frame, frame.mirrored()] {
                let circle: Curve = ogeom_geom::CircleCurve::new(
                    ogeom_math::Circle::new(frame, *radius, T).unwrap(),
                )
                .into();
                let pcurve = exact_pcurve_of(&circle, surface, T)
                    .expect("a circle of the surface has a closed-form pcurve");
                for i in 0..=8 {
                    let t = core::f64::consts::TAU * f64::from(i) / 8.0;
                    let uv = pcurve.point_at(t, T).unwrap();
                    let lifted = surface
                        .point_at(uv.x.rem_euclid(core::f64::consts::TAU), uv.y, T)
                        .or_else(|_| surface.point_at(uv.x, uv.y, T))
                        .unwrap();
                    let p3 = circle.point_at(t, T).unwrap();
                    // Both sides are closed forms of one point: rounding
                    // only.
                    assert!(
                        p3.distance(lifted) < 1e-9,
                        "{:?} frame, t {t}: pcurve lifts {lifted:?} against {p3:?}",
                        frame.handedness()
                    );
                }
            }
        }
    }

    /// A surface on a left-handed (mirrored) frame turns its angle `u`
    /// about `-z`: a circle round its axis that turns about `+z` runs down
    /// `u`, not up it. Each surface is taken on both frames through its own
    /// `x`, each circle on it on both of its own, and the circle's pcurve
    /// lifted through the surface must give the curve's own point back at
    /// every parameter.
    #[test]
    fn a_circle_on_a_left_handed_surface_keeps_its_pcurve_same_parameter() {
        use ogeom_geom::{ConeSurface, TorusSurface};
        let on = |centre: Point, normal: Direction, radius: f64| {
            (Frame::new(centre, normal, Direction::X, T).unwrap(), radius)
        };
        let cylinder = |frame: Frame| -> SurfaceGeometry {
            CylinderSurface::new(Cylinder::new(frame, 2.0, T).unwrap(), (-4.0, 4.0))
                .unwrap()
                .into()
        };
        let cone = |frame: Frame| -> SurfaceGeometry {
            ConeSurface::new(
                ogeom_math::Cone::new(frame, 1.0, core::f64::consts::FRAC_PI_6, T).unwrap(),
                (-1.0, 3.0),
            )
            .unwrap()
            .into()
        };
        let sphere = |frame: Frame| -> SurfaceGeometry {
            SphereSurface::new(Sphere::new(frame, 2.0, T).unwrap()).into()
        };
        let torus = |frame: Frame| -> SurfaceGeometry {
            TorusSurface::new(ogeom_math::Torus::new(frame, 2.0, 0.5, T).unwrap()).into()
        };
        type Make<'a> = &'a dyn Fn(Frame) -> SurfaceGeometry;
        let surfaces: [(Make, Vec<(Frame, f64)>); 4] = [
            (
                &cylinder,
                vec![on(Point::new(0.0, 0.0, 1.0), Direction::Z, 2.0)],
            ),
            (
                &cone,
                vec![on(
                    Point::new(0.0, 0.0, 1.0),
                    Direction::Z,
                    1.0 + core::f64::consts::FRAC_PI_6.tan(),
                )],
            ),
            (
                &sphere,
                vec![
                    on(Point::new(0.0, 0.0, 1.0), Direction::Z, 3.0_f64.sqrt()),
                    on(Point::ORIGIN, Direction::Y, 2.0),
                ],
            ),
            (
                &torus,
                vec![
                    on(Point::ORIGIN, Direction::Z, 2.5),
                    on(Point::new(0.0, 0.0, 0.5), Direction::Z, 2.0),
                    on(Point::new(2.0, 0.0, 0.0), Direction::Y, 0.5),
                ],
            ),
        ];
        let tau = core::f64::consts::TAU;
        for (make, circles) in &surfaces {
            for surface_frame in [Frame::WORLD, Frame::WORLD.mirrored()] {
                let surface = make(surface_frame);
                for (frame, radius) in circles {
                    for frame in [*frame, frame.mirrored()] {
                        let circle: Curve = ogeom_geom::CircleCurve::new(
                            ogeom_math::Circle::new(frame, *radius, T).unwrap(),
                        )
                        .into();
                        // A meridian's half is imaged over its own range.
                        let range = if frame.z().vector().dot(Vector::Z).abs() < 0.5 {
                            let half = core::f64::consts::FRAC_PI_2 - 0.1;
                            (-half, half)
                        } else {
                            (0.0, tau)
                        };
                        let pcurve = exact_pcurve_over(&circle, range, &surface, T)
                            .expect("a circle of the surface has a closed-form pcurve");
                        for i in 0..=8 {
                            let t = (range.1 - range.0).mul_add(f64::from(i) / 8.0, range.0);
                            let uv = pcurve.point_at(t, T).unwrap();
                            let lifted = surface
                                .point_at(uv.x.rem_euclid(tau), uv.y, T)
                                .or_else(|_| {
                                    surface.point_at(uv.x.rem_euclid(tau), uv.y.rem_euclid(tau), T)
                                })
                                .unwrap();
                            let p3 = circle.point_at(t, T).unwrap();
                            // Both sides are closed forms of one point:
                            // rounding only.
                            assert!(
                                p3.distance(lifted) < 1e-9,
                                "{:?} surface, {:?} circle {frame:?}, t {t}: pcurve lifts \
                                 {lifted:?} against {p3:?}",
                                surface_frame.handedness(),
                                frame.handedness()
                            );
                        }
                    }
                }
            }
        }
    }

    /// A plane through a ball's own axis cuts a meridian. The whole circle has
    /// no chart image (its longitude jumps half a turn at each pole), but
    /// each half is a straight line in the chart, exactly, at the circle's own
    /// parameter. Pinned by lifting the line back through the sphere and
    /// demanding the circle's point, on every half of every orientation.
    #[test]
    fn a_meridian_half_has_an_exact_line_for_a_pcurve() {
        use ogeom_geom::Surface as _;
        let half = core::f64::consts::PI;
        for (centre, radius) in [(Point::ORIGIN, 4.0), (Point::new(1.0, -2.0, 0.5), 1.25)] {
            let ball = sphere(centre, radius);
            let SurfaceGeometry::Sphere(s) = &ball else {
                panic!("a sphere surface");
            };
            // Three planes through the axis, at different azimuths, so the
            // constant longitude is not accidentally zero.
            for azimuth in [0.0_f64, 0.7, 2.4] {
                let normal = Vector::new(-azimuth.sin(), azimuth.cos(), 0.0);
                let cut = plane(centre, normal);
                let SurfaceIntersection::Along(curves) =
                    intersect_surfaces(&ball, &cut, IntersectOptions::default(), T).unwrap()
                else {
                    panic!("a plane through the centre meets the ball along a circle");
                };
                assert_eq!(curves.len(), 1, "one great circle");
                let circle = &curves[0].curve;
                assert!(curves[0].exact);
                // The whole circle has no chart image; each half does.
                assert!(
                    exact_pcurve_over(circle, circle.domain(), &ball, T).is_none(),
                    "the whole meridian has no single chart image"
                );
                for (lo, hi) in [(0.0, half), (half, 2.0 * half), (0.3, half - 0.1)] {
                    let pcurve = exact_pcurve_over(circle, (lo, hi), &ball, T)
                        .expect("half a meridian has an exact pcurve");
                    assert!(
                        matches!(pcurve, PlanarCurve::Line(_)),
                        "and it is a straight line in the chart"
                    );
                    for i in 0..=16 {
                        let t = (hi - lo).mul_add(f64::from(i) / 16.0, lo);
                        let want = circle.point_at(t, T).unwrap();
                        let uv = pcurve.point_at(t, T).unwrap();
                        assert!(
                            uv.y >= -half.mul_add(0.5, 1e-12) && uv.y <= half.mul_add(0.5, 1e-12),
                            "the latitude stays inside the chart: {}",
                            uv.y
                        );
                        let lifted = ball
                            .point_at(uv.x.rem_euclid(core::f64::consts::TAU), uv.y, T)
                            .unwrap();
                        assert!(
                            want.distance(lifted) < 1e-9,
                            "azimuth {azimuth}, t {t}: {lifted:?} against {want:?}"
                        );
                    }
                }
                // A range straddling a pole has none, and says so rather than
                // answering for one side.
                assert!(
                    exact_pcurve_over(circle, (half - 0.2, half + 0.2), &ball, T).is_none(),
                    "a range across a pole has no one line"
                );
                let _ = s;
            }
        }
    }

    /// A half meridian whose ends were read off the poles by projection
    /// runs a few nanoradians past them. That is the projection's own
    /// error, far inside the confusion distance, and the half is still the
    /// straight line it is, lifting onto the circle across the whole range.
    #[test]
    fn a_meridian_half_read_slightly_past_its_poles_is_still_a_line() {
        let half = core::f64::consts::PI;
        let ball = sphere(Point::new(1.0, -2.0, 0.5), 0.2);
        let cut = plane(
            Point::new(1.0, -2.0, 0.5),
            Vector::new(-0.7_f64.sin(), 0.7_f64.cos(), 0.0),
        );
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&ball, &cut, IntersectOptions::default(), T).unwrap()
        else {
            panic!("a plane through the centre meets the ball along a circle");
        };
        let circle = &curves[0].curve;
        let SurfaceGeometry::Sphere(s) = &ball else {
            panic!("a sphere surface");
        };
        for over in [2e-9, 1e-8, 1e-7] {
            for (lo, hi) in [(-over, half + over), (half - over, 2.0 * half + over)] {
                let pcurve = exact_pcurve_over(circle, (lo, hi), &ball, T)
                    .unwrap_or_else(|| panic!("{over} past the poles: no line"));
                for i in 0..=16 {
                    let t = (hi - lo).mul_add(f64::from(i) / 16.0, lo);
                    let uv = pcurve.point_at(t, T).unwrap();
                    let lifted = ogeom_math::elementary::sphere_at(&s.sphere(), uv.x, uv.y).point;
                    let want = circle.point_at(t, T).unwrap();
                    assert!(
                        want.distance(lifted) < 1e-9,
                        "{over}, t {t}: {lifted:?} against {want:?}"
                    );
                }
            }
        }
    }

    /// A trim says *where* on a curve, not what it is. The basis carries the
    /// shape and the trim shares its parameter, so a trimmed curve's pcurve is
    /// the basis's own pcurve trimmed the same way, on every surface, since
    /// the answer does not depend on the surface at all.
    ///
    /// A fillet's own end cap is a plane and the edges bounding it are trimmed
    /// curves. Without this, the boolean refuses the coincidence because it
    /// cannot put a trimmed curve into a chart it plainly lies in.
    #[test]
    fn a_trimmed_curve_carries_its_basis_pcurve_trimmed_the_same_way() {
        use ogeom_geom::TrimmedCurve;
        let drum = cylinder(Vector::Z, 2.0);
        let ground = plane(Point::new(0.0, 0.0, 1.0), Vector::Z);
        // The circle where they meet, and a quarter of it.
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&drum, &ground, IntersectOptions::default(), T).unwrap()
        else {
            panic!("a plane across a cylinder meets it in a circle");
        };
        let whole = curves[0].curve.clone();
        let (lo, hi) = whole.domain();
        let quarter: Curve = TrimmedCurve::new(whole.clone(), lo + 0.3, lo + (hi - lo) / 4.0, T)
            .unwrap()
            .into();

        for surface in [&drum, &ground] {
            let full = exact_pcurve_of(&whole, surface, T).expect("the whole circle has one");
            let part = exact_pcurve_of(&quarter, surface, T).expect("and so does a quarter of it");
            // Same parameter, same point: the trim changed the range and
            // nothing else.
            let (a, b) = quarter.domain();
            for i in 0..=8 {
                let t = (b - a).mul_add(f64::from(i) / 8.0, a);
                let (whole_at, part_at) =
                    (full.point_at(t, T).unwrap(), part.point_at(t, T).unwrap());
                assert!(
                    whole_at.distance(part_at) < 1e-12,
                    "the trim carries the basis: {whole_at:?} against {part_at:?}"
                );
                // And it lifts back onto the curve it came from.
                let lifted = surface
                    .point_at(part_at.x.rem_euclid(core::f64::consts::TAU), part_at.y, T)
                    .or_else(|_| surface.point_at(part_at.x, part_at.y, T))
                    .unwrap();
                assert!(
                    lifted.distance(quarter.point_at(t, T).unwrap()) < 1e-9,
                    "same-parameter, still"
                );
            }
        }
    }

    /// An arc whose range runs past the end of its circle's turn (from late
    /// in the turn across the circle's origin) has a pcurve on a surface of
    /// revolution over that whole range, the same point at every parameter
    /// of it.
    #[test]
    fn an_arc_across_its_circles_origin_has_a_pcurve_over_its_range() {
        let drum = cylinder(Vector::Z, 2.0);
        let ball = sphere(Point::ORIGIN, 2.0);
        let ground = plane(Point::ORIGIN, Vector::Z);
        for surface in [&drum, &ball] {
            let SurfaceIntersection::Along(curves) =
                intersect_surfaces(surface, &ground, IntersectOptions::default(), T).unwrap()
            else {
                panic!("a plane through the axis's normal meets it in a circle");
            };
            let circle = curves[0].curve.clone();
            let tau = core::f64::consts::TAU;
            for range in [(4.7, tau + 1.0), (-1.0, 1.5)] {
                let pcurve = exact_pcurve_over(&circle, range, surface, T).expect("a parallel");
                for i in 0..=8 {
                    let t = (range.1 - range.0).mul_add(f64::from(i) / 8.0, range.0);
                    let at = pcurve.point_at(t, T).unwrap();
                    let lifted = surface.point_at(at.x.rem_euclid(tau), at.y, T).unwrap();
                    assert!(
                        lifted.distance(circle.point_at(t, T).unwrap()) < 1e-9,
                        "{range:?} at {t}"
                    );
                }
            }
        }
    }

    /// A plane through a cone's apex: tangent, it touches along one ruling;
    /// steeper, it holds two; shallower, it meets the apex alone. Every
    /// ruling lies on both surfaces and stays within the cone's window.
    #[test]
    fn a_plane_through_a_cones_apex_holds_its_rulings() {
        use ogeom_geom::ConeSurface;
        let frame = Frame::new(
            Point::new(100.0, 200.0, 300.0),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        let cone = ogeom_math::Cone::new(frame, 10.0, core::f64::consts::FRAC_PI_4, T).unwrap();
        let surface: SurfaceGeometry = ConeSurface::new(cone, (-5.0, 50.0)).unwrap().into();
        let apex = Point::new(100.0, 200.0, 290.0);
        let plane = |normal: Vector| -> SurfaceGeometry {
            PlaneSurface::new(Plane::through(apex, Direction::new(normal, T).unwrap())).into()
        };
        let cases = [
            (Vector::new(1.0, 0.0, -1.0), 1, true),
            (Vector::new(1.0, 0.0, 0.0), 2, false),
        ];
        for (normal, count, tangent) in cases {
            let cut = plane(normal);
            let SurfaceIntersection::Along(sections) =
                intersect_surfaces(&surface, &cut, IntersectOptions::default(), T).unwrap()
            else {
                panic!("{normal:?}: rulings");
            };
            assert_eq!(sections.len(), count, "{normal:?}");
            for section in &sections {
                assert_eq!(section.tangential, tangent, "{normal:?}");
                let (lo, hi) = section.curve.domain();
                for k in 0..=4 {
                    let p = section
                        .curve
                        .point_at(lo + (hi - lo) * f64::from(k) / 4.0, T)
                        .unwrap();
                    assert!(cone.distance_to(p) < 1e-9, "{p:?} on the cone");
                    let height = p.z - 300.0;
                    assert!(
                        (-5.0 - 1e-9..=50.0 + 1e-9).contains(&height),
                        "{p:?} in the window"
                    );
                }
            }
        }
        let shallow = plane(Vector::new(0.2, 0.0, 1.0));
        assert!(matches!(
            intersect_surfaces(&surface, &shallow, IntersectOptions::default(), T).unwrap(),
            SurfaceIntersection::Touching(_) | SurfaceIntersection::Apart
        ));
    }

    #[test]
    fn a_far_stated_ruling_reads_its_angle_on_the_used_nappe() {
        use ogeom_geom::ConeSurface;
        // A 45-degree cone opening along +z, reference radius 24 at the
        // frame's origin; a ruling at chart angle 0.01, exactly as a real
        // file states it: the line's own origin parked seven hundred
        // kilometres down the infinite line, past the apex on the other
        // nappe. Only the used range may vote on the angle, or the pcurve
        // lands half a turn away and the face triangulates as a fan across
        // the whole chart.
        let cone =
            ogeom_math::Cone::new(Frame::WORLD, 24.0, core::f64::consts::FRAC_PI_4, T).unwrap();
        let surface: SurfaceGeometry = ConeSurface::new(cone, (-1e5, 1e5)).unwrap().into();
        let u_true = 0.01_f64;
        let radial = Vector::new(u_true.cos(), u_true.sin(), 0.0);
        // The ruling climbs outward at 45 degrees; its stated origin sits
        // far beyond the apex (z = -24 on this cone), on the other nappe.
        let direction =
            Direction::new((radial + Vector::new(0.0, 0.0, 1.0)) / 2f64.sqrt(), T).unwrap();
        let far = -7.0e5;
        let origin = Point::ORIGIN + radial * 24.0 + direction.vector() * far;
        let line = ogeom_geom::LineCurve::over(
            ogeom_math::Axis::new(origin, direction),
            far.abs() - 1.0,
            far.abs() + 1.0,
        )
        .unwrap();
        let curve: Curve = line.into();
        let range = ogeom_geom::Curve3d::domain(&curve);
        let pcurve = exact_pcurve_over(&curve, range, &surface, T).expect("a ruling inverts");
        let at = pcurve.point_at(range.0, T).unwrap();
        let tau = core::f64::consts::TAU;
        let gap = (at.x - u_true)
            .rem_euclid(tau)
            .min(tau - (at.x - u_true).rem_euclid(tau));
        assert!(
            gap < 1e-6,
            "the ruling's chart angle must be the used side's: got u {} against {u_true}",
            at.x
        );
    }

    /// A thin drum crossing a wide one obliquely, its far side passing two
    /// hundredths of a millimetre inside the wide wall: two loops (the
    /// second cut short by the wide drum's end), each turning sharply at
    /// the tip where the drums all but touch. Each comes back as curves
    /// within a few microns of both drums; the whole first loop fitted as
    /// one misses its trace by millimetres.
    #[test]
    fn loops_turning_sharply_where_drums_all_but_touch_fit_in_pieces() {
        let drum = |origin: Point, axis: Vector, x: Vector, radius: f64, height: (f64, f64)| {
            let frame = Frame::new(
                origin,
                Direction::new(axis, T).unwrap(),
                Direction::new(x, T).unwrap(),
                T,
            )
            .unwrap();
            let cylinder = Cylinder::new(frame, radius, T).unwrap();
            (
                cylinder,
                SurfaceGeometry::from(CylinderSurface::new(cylinder, height).unwrap()),
            )
        };
        let radius = 3.175;
        let (thin_drum, thin) = drum(
            Point::new(10.994_218_762_109_735, 53.975, -209.55),
            Vector::new(0.0, 0.0, -1.0),
            Vector::new(-1.0, 0.0, 0.0),
            radius,
            (-1000.0, 1000.0),
        );
        let (wide_drum, wide) = drum(
            Point::new(
                -14.478_610_818_124_423,
                3.999_371_635_181_902,
                -277.138_401_510_994_7,
            ),
            Vector::new(
                -0.565_016_635_381_368_8,
                0.528_780_602_945_684_7,
                0.633_361_883_673_714_3,
            ),
            Vector::new(0.0, 0.767_637_390_304_296_3, -0.640_884_417_821_817),
            57.088_766_757_544_09,
            (0.0, 171.306_147_621_762_6),
        );
        let options = IntersectOptions {
            tolerance: 1e-5,
            marching: crate::Marching {
                chord: 1e-5,
                ..crate::Marching::default()
            },
        };
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&thin, &wide, options, T).unwrap()
        else {
            panic!("the drums cross");
        };
        let mut length = 0.0;
        for section in &curves {
            assert!(!section.tangential);
            assert!(
                section.tolerance < 1e-3,
                "a section states {} of doubt",
                section.tolerance
            );
            let (lo, hi) = section.curve.domain();
            let mut previous = section.curve.point_at(lo, T).unwrap();
            for i in 1..=400 {
                let t = lo + (hi - lo) * f64::from(i) / 400.0;
                let p = section.curve.point_at(t, T).unwrap();
                let off = thin_drum
                    .distance_to(p)
                    .abs()
                    .max(wide_drum.distance_to(p).abs());
                assert!(off < 1e-3, "a section stands {off} off the drums");
                length += previous.distance(p);
                previous = p;
            }
        }
        // Each loop runs the thin drum's girth stretched along the wide
        // wall: the two together far longer than two girths.
        assert!(
            length > 4.0 * core::f64::consts::PI * radius,
            "the sections cover {length}"
        );
    }

    #[test]
    fn coincidence_is_measured_over_the_overlap_and_nowhere_else() {
        // Patches restated from planes: the geometry no longer says "plane",
        // which is the whole reason this measurement exists.
        let patch = |plane: ogeom_math::Plane, u: (f64, f64), v: (f64, f64)| {
            let surface: SurfaceGeometry = PlaneSurface::over(plane, u, v).unwrap().into();
            SurfaceGeometry::from(surface.to_bspline(T).unwrap())
        };
        let reach = T.confusion() * 1e2;

        // Two windows on one plane, overlapping over a quarter of each. They
        // are the same surface exactly where they meet, which is the claim.
        let here = patch(ogeom_math::Plane::XY, (0.0, 10.0), (0.0, 10.0));
        let over = patch(ogeom_math::Plane::XY, (5.0, 15.0), (5.0, 15.0));
        assert!(surfaces_coincide(&here, &over, reach, T));

        // The same plane lifted clear of itself is not the same surface, and
        // a plane square to it crosses rather than coincides: the case that
        // must keep marching, since a crossing has a section to find.
        let above = patch(
            ogeom_math::Plane::new(
                Frame::new(Point::new(0.0, 0.0, 1.0), Direction::Z, Direction::X, T).unwrap(),
            ),
            (0.0, 10.0),
            (0.0, 10.0),
        );
        assert!(!surfaces_coincide(&here, &above, reach, T));
        let across = patch(
            ogeom_math::Plane::new(
                Frame::new(Point::new(5.0, 0.0, 0.0), Direction::X, Direction::Y, T).unwrap(),
            ),
            (0.0, 10.0),
            (0.0, 10.0),
        );
        assert!(!surfaces_coincide(&here, &across, reach, T));
    }

    /// A curve offset square to a plane's normal, and a parabola, lying in
    /// the plane, have exact pcurves on it, the same point at every
    /// parameter: an offset of a B-spline and of a circle, the reference
    /// along the normal or against it, the distance either sign. An
    /// offset whose reference leans off the normal leaves the plane and
    /// has none.
    #[test]
    fn an_offset_square_to_a_plane_and_a_parabola_have_pcurves_on_it() {
        use ogeom_geom::{BSplineCurve, CircleCurve, OffsetCurve, ParabolaCurve};
        let normal = Vector::new(0.3, -0.5, 0.8);
        let up = Direction::new(normal, T).unwrap();
        let origin = Point::new(1.0, 2.0, -0.5);
        let ground = plane(origin, normal);
        let frame = Frame::new(
            origin,
            up,
            Direction::new(Vector::new(1.0, 0.0, -0.375), T).unwrap(),
            T,
        )
        .unwrap();
        let on = |x: f64, y: f64| frame.origin() + frame.x().vector() * x + frame.y().vector() * y;
        let knots = ogeom_math::KnotVector::clamped_uniform(3, 5).unwrap();
        let wiggle: Curve = BSplineCurve::new(
            knots,
            vec![
                on(0.0, 0.0),
                on(1.0, 1.0),
                on(2.0, -0.5),
                on(3.0, 0.5),
                on(4.0, 0.0),
            ],
            T,
        )
        .unwrap()
        .into();
        let round: Curve = CircleCurve::new(ogeom_math::Circle::new(frame, 1.5, T).unwrap()).into();
        let parabola: Curve =
            ParabolaCurve::over(ogeom_math::Parabola::new(frame, 0.7, T).unwrap(), -1.5, 2.0)
                .unwrap()
                .into();
        let mut curves = vec![parabola];
        for basis in [wiggle, round] {
            for reference in [up, -up] {
                for distance in [0.4, -0.3] {
                    let offset = OffsetCurve::new(basis.clone(), distance, reference).unwrap();
                    curves.push(Curve::Offset(Box::new(offset)));
                }
            }
        }
        for curve in &curves {
            let pcurve = exact_pcurve_of(curve, &ground, T).expect("it lies in the plane");
            let (lo, hi) = curve.domain();
            for i in 0..=16 {
                let t = (hi - lo).mul_add(f64::from(i) / 16.0, lo);
                let at = pcurve.point_at(t, T).unwrap();
                let lifted = ground.point_at(at.x, at.y, T).unwrap();
                let point = curve.point_at(t, T).unwrap();
                assert!(
                    lifted.distance(point) < 1e-9,
                    "{curve:?} at {t}: {lifted:?} against {point:?}"
                );
            }
        }
        // Leaned by `ε`, the reference lifts a point offset `d` off the
        // plane by up to `d tan ε`, here past the confusion distance.
        let leaned = Direction::new(normal + frame.x().vector() * 0.05, T).unwrap();
        let knots = ogeom_math::KnotVector::clamped_uniform(3, 4).unwrap();
        let straight: Curve = BSplineCurve::new(
            knots,
            vec![on(0.0, 0.0), on(1.0, 0.0), on(2.0, 0.0), on(3.0, 0.0)],
            T,
        )
        .unwrap()
        .into();
        let off = Curve::Offset(Box::new(OffsetCurve::new(straight, 0.4, leaned).unwrap()));
        assert!(exact_pcurve_of(&off, &ground, T).is_none());
    }

    /// A hyperbola branch in a plane has an exact pcurve there at the same
    /// parameter, with its derivatives: forward, on a frame turned over, run
    /// backwards over a domain symmetric about zero or over any other, and
    /// offset.
    #[test]
    fn a_hyperbola_in_a_plane_has_its_pcurve_there() {
        use ogeom_geom::{HyperbolaCurve, OffsetCurve, Reversible as _};
        let normal = Vector::new(0.3, -0.5, 0.8);
        let up = Direction::new(normal, T).unwrap();
        let origin = Point::new(1.0, 2.0, -0.5);
        let ground = plane(origin, normal);
        let x = Direction::new(Vector::new(1.0, 0.0, -0.375), T).unwrap();
        let base = Frame::new(origin, up, x, T).unwrap();
        let centre = origin + base.x().vector() * 0.2 + base.y().vector() * 0.1;
        let frame = Frame::new(centre, up, x, T).unwrap();
        let turned_over =
            Frame::from_axes(frame.origin(), frame.x(), -frame.y(), -frame.z(), T).unwrap();
        let branch = |frame: Frame, lo: f64, hi: f64| -> Curve {
            HyperbolaCurve::over(
                ogeom_math::Hyperbola::new(frame, 1.5, 0.8, T).unwrap(),
                lo,
                hi,
            )
            .unwrap()
            .into()
        };
        let forward = branch(frame, -1.2, 0.7);
        let mut curves = vec![
            forward.clone(),
            branch(turned_over, -1.2, 0.7),
            branch(frame, -1.2, 1.2).reversed(),
            branch(frame, -1.2, 0.7).reversed(),
        ];
        for reference in [up, -up] {
            let offset = OffsetCurve::new(forward.clone(), 0.3, reference).unwrap();
            curves.push(Curve::Offset(Box::new(offset)));
        }
        for curve in &curves {
            let pcurve = exact_pcurve_of(curve, &ground, T).expect("it lies in the plane");
            let (lo, hi) = curve.domain();
            for i in 0..=16 {
                let t = (hi - lo).mul_add(f64::from(i) / 16.0, lo);
                let at = pcurve.point_at(t, T).unwrap();
                let lifted = ground.point_at(at.x, at.y, T).unwrap();
                let point = curve.point_at(t, T).unwrap();
                assert!(
                    lifted.distance(point) < 1e-9,
                    "{curve:?} at {t}: {lifted:?} against {point:?}"
                );
                if matches!(curve, Curve::Offset(_)) {
                    continue;
                }
                // The plane's chart is an isometry, so each derivative in
                // the chart, lifted along its axes, is the curve's own.
                let flat = pcurve.derivatives_at(t, 2, T).unwrap();
                let space = curve.derivatives_at(t, 2, T).unwrap();
                let axes = ground_axes(&ground);
                for (order, (f, s)) in flat.iter().zip(&space).enumerate().skip(1) {
                    let lifted = axes.0 * f.x + axes.1 * f.y;
                    assert!(
                        (lifted - *s).magnitude() < 1e-9,
                        "{curve:?} at {t}, order {order}: {lifted:?} against {s:?}"
                    );
                }
            }
        }
    }

    /// Each lifted point of `pcurve` on `surface`, and its first derivative
    /// lifted along the chart's tangents, against the curve's own at 65
    /// parameters over `range`, to 1e-12.
    fn lifts_onto(
        curve: &Curve,
        pcurve: &PlanarCurve,
        range: (f64, f64),
        surface: &SurfaceGeometry,
    ) {
        use ogeom_geom::Surface as _;
        for i in 0..=64 {
            let t = (range.1 - range.0).mul_add(f64::from(i) / 64.0, range.0);
            let at = pcurve.derivatives_at(t, 1, T).unwrap();
            let lifted = surface.point_at(at[0].x, at[0].y, T).unwrap();
            let point = curve.point_at(t, T).unwrap();
            assert!(
                lifted.distance(point) < 1e-12,
                "{curve:?} at {t}: {lifted:?} against {point:?}"
            );
            let (su, sv) = surface.d1_at(at[0].x, at[0].y, T).unwrap();
            let tangent = su * at[1].x + sv * at[1].y;
            let own = curve.d1_at(t, T).unwrap();
            assert!(
                (tangent - own).magnitude() < 1e-12,
                "{curve:?} at {t}: tangent {tangent:?} against {own:?}"
            );
        }
    }

    /// A trim run backwards reads its basis at `s + e - u` over its window
    /// `[s, e]`, and its pcurve is the basis's trimmed the same way and
    /// reversed: a circle arc inside its turn and across its seam, a
    /// hyperbola arc, and a trim of a trim, each at the same parameter.
    #[test]
    fn a_reversed_trim_in_a_plane_has_its_pcurve_there() {
        use ogeom_geom::{CircleCurve, HyperbolaCurve, Reversible as _, TrimmedCurve};
        let normal = Vector::new(0.3, -0.5, 0.8);
        let up = Direction::new(normal, T).unwrap();
        let origin = Point::new(1.0, 2.0, -0.5);
        let ground = plane(origin, normal);
        let x = Direction::new(Vector::new(1.0, 0.0, -0.375), T).unwrap();
        let base = Frame::new(origin, up, x, T).unwrap();
        let centre = origin + base.x().vector() * 0.2 + base.y().vector() * 0.1;
        let frame = Frame::new(centre, up, x, T).unwrap();
        let circle: Curve =
            CircleCurve::new(ogeom_math::Circle::new(frame, 1.3, T).unwrap()).into();
        let hyperbola: Curve = HyperbolaCurve::over(
            ogeom_math::Hyperbola::new(frame, 1.5, 0.8, T).unwrap(),
            -1.2,
            0.7,
        )
        .unwrap()
        .into();
        let trim = |basis: &Curve, lo: f64, hi: f64| -> Curve {
            TrimmedCurve::new(basis.clone(), lo, hi, T).unwrap().into()
        };
        let curves = [
            trim(&circle, 0.4, 2.9).reversed(),
            trim(&circle, 5.5, 7.1).reversed(),
            trim(&circle.reversed(), 0.4, 2.9).reversed(),
            trim(&hyperbola, -0.9, 0.5).reversed(),
            trim(&trim(&hyperbola, -1.0, 0.6).reversed(), -0.9, 0.5).reversed(),
        ];
        for curve in &curves {
            let pcurve = exact_pcurve_of(curve, &ground, T).expect("it lies in the plane");
            lifts_onto(curve, &pcurve, curve.domain(), &ground);
        }
    }

    /// A ruling on a cone, trimmed and run backwards: its pcurve is the
    /// ruling's column trimmed the same way and reversed, at the same
    /// parameter.
    #[test]
    fn a_reversed_trim_of_a_ruling_on_a_cone_has_its_pcurve_there() {
        use ogeom_geom::{ConeSurface, Reversible as _, TrimmedCurve};
        let cone =
            ogeom_math::Cone::new(Frame::WORLD, 24.0, core::f64::consts::FRAC_PI_4, T).unwrap();
        let surface: SurfaceGeometry = ConeSurface::new(cone, (-1e3, 1e3)).unwrap().into();
        let radial = Vector::new(0.3_f64.cos(), 0.3_f64.sin(), 0.0);
        let direction =
            Direction::new((radial + Vector::new(0.0, 0.0, 1.0)) / 2f64.sqrt(), T).unwrap();
        let apex = Point::new(0.0, 0.0, -24.0);
        let line: Curve =
            ogeom_geom::LineCurve::over(ogeom_math::Axis::new(apex, direction), 5.0, 50.0)
                .unwrap()
                .into();
        let curve = Curve::from(TrimmedCurve::new(line, 10.0, 30.0, T).unwrap()).reversed();
        let pcurve = exact_pcurve_of(&curve, &surface, T).expect("a ruling inverts");
        lifts_onto(&curve, &pcurve, curve.domain(), &surface);
    }

    /// A ruling past a cone's apex lies on the far nappe, where the chart's
    /// radius runs negative, so its chart column is half a turn from its own
    /// bearing there. The ruling is used wholly on either side of the apex,
    /// through the apex with its farther end on either side, from an origin
    /// stated past the apex on the side it does not use, and run backwards.
    /// Each lifts onto its line at the same parameter, on a cone widening up
    /// its axis and on one narrowing.
    #[test]
    fn a_ruling_past_the_apex_has_its_pcurve_on_its_own_nappe() {
        use ogeom_geom::{ConeSurface, Reversible as _};
        let frame = Frame::new(
            Point::new(1.0, -2.0, 0.5),
            Direction::new(Vector::new(0.2, -0.3, 1.0), T).unwrap(),
            Direction::from_cross(Vector::new(0.2, -0.3, 1.0), Vector::new(0.0, 1.0, 0.0), T)
                .unwrap(),
            T,
        )
        .unwrap();
        for half_angle in [core::f64::consts::FRAC_PI_4, -core::f64::consts::FRAC_PI_6] {
            let cone = ogeom_math::Cone::new(frame, 3.0, half_angle, T).unwrap();
            let surface: SurfaceGeometry = ConeSurface::new(cone, (-200.0, 200.0)).unwrap().into();
            // The column at angle 0.3, climbing the chart's height: through
            // the apex at the height where the radius is zero.
            let u = 0.3_f64;
            let (x, y, z) = (frame.x().vector(), frame.y().vector(), frame.z().vector());
            let bearing = x * u.cos() + y * u.sin();
            let apex = frame.origin() + z * (-3.0 / half_angle.tan());
            let up = Direction::new(bearing * half_angle.tan() + z, T).unwrap();
            let stated = |origin: f64, range: (f64, f64)| -> Curve {
                let at = apex + up.vector() * origin;
                ogeom_geom::LineCurve::over(
                    ogeom_math::Axis::new(at, up),
                    range.0 - origin,
                    range.1 - origin,
                )
                .unwrap()
                .into()
            };
            let lines = [
                stated(0.0, (5.0, 40.0)),
                stated(0.0, (-40.0, -5.0)),
                stated(0.0, (-10.0, 30.0)),
                stated(0.0, (-30.0, 10.0)),
                stated(60.0, (-40.0, -5.0)),
                stated(-60.0, (5.0, 40.0)),
                stated(0.0, (-30.0, 10.0)).reversed(),
                stated(0.0, (-40.0, -5.0)).reversed(),
            ];
            for line in &lines {
                let range = line.domain();
                let pcurve = exact_pcurve_over(line, range, &surface, T).expect("a ruling inverts");
                lifts_onto(line, &pcurve, range, &surface);
            }
        }
    }

    /// A curve drawn on a plane under another chart has its pcurve in the
    /// plane's own, at the same parameter, wherever that chart's plane is
    /// this one: the same chart, one turned and moved within the plane,
    /// one whose normal is turned round (a mirror image in the chart), and
    /// one trimmed. A chart on a plane tipped off this one, or standing off
    /// it, is another plane, and has none.
    #[test]
    fn a_curve_drawn_on_the_plane_under_any_chart_has_its_pcurve_there() {
        use ogeom_geom::{
            BSpline2d, Circle2d, CurveOnSurface, Line2d, PlaneSurface, Reversible as _,
            TrimmedSurface,
        };
        let normal = Vector::new(0.3, -0.5, 0.8);
        let up = Direction::new(normal, T).unwrap();
        let origin = Point::new(1.0, 2.0, -0.5);
        let ground = plane(origin, normal);
        let x = Direction::new(Vector::new(1.0, 0.0, -0.375), T).unwrap();
        let own = Frame::new(origin, up, x, T).unwrap();
        let elsewhere = origin + own.x().vector() * 0.4 + own.y().vector() * 0.5;
        let chart = |frame: Frame| -> SurfaceGeometry {
            PlaneSurface::over(Plane::new(frame), (-6.0, 6.0), (-6.0, 6.0))
                .unwrap()
                .into()
        };
        let turned = Frame::new(
            elsewhere,
            up,
            Direction::from_cross(normal, Vector::new(0.2, 0.9, 0.1), T).unwrap(),
            T,
        )
        .unwrap();
        let mirrored = Frame::new(elsewhere, -up, x, T).unwrap();
        let trimmed: SurfaceGeometry =
            TrimmedSurface::new(chart(turned), (-3.0, 3.0), (-3.0, 3.0), T)
                .unwrap()
                .into();
        let ring = ogeom_math::Circle2::new(
            ogeom_math::Frame2::new(Point2::new(0.3, -0.2), ogeom_math::Direction2::X),
            0.9,
            T,
        )
        .unwrap();
        let knots = ogeom_math::KnotVector::clamped_uniform(3, 5).unwrap();
        let wave = BSpline2d::new(
            knots,
            vec![
                Point2::new(-1.0, 0.0),
                Point2::new(-0.5, 0.7),
                Point2::new(0.0, -0.4),
                Point2::new(0.6, 0.5),
                Point2::new(1.1, 0.1),
            ],
            T,
        )
        .unwrap();
        let pcurves: [PlanarCurve; 4] = [
            Line2d::segment(Point2::new(-0.5, 0.25), Point2::new(1.0, 0.75), T)
                .unwrap()
                .into(),
            Circle2d::new(ring).into(),
            PlanarCurve::from(Circle2d::new(ring)).reversed(),
            wave.into(),
        ];
        for surface in [chart(own), chart(turned), chart(mirrored), trimmed] {
            for drawn in &pcurves {
                let curve = Curve::OnSurface(Box::new(CurveOnSurface::new(
                    drawn.clone(),
                    surface.clone(),
                )));
                let pcurve = exact_pcurve_of(&curve, &ground, T).expect("it lies in the plane");
                lifts_onto(&curve, &pcurve, curve.domain(), &ground);
            }
        }
        let tipped = Frame::new(
            elsewhere,
            Direction::new(normal + Vector::new(1e-3, 0.0, 0.0), T).unwrap(),
            x,
            T,
        );
        let standing = Frame::new(elsewhere + up.vector() * 1e-3, up, x, T).unwrap();
        for frame in [tipped.unwrap(), standing] {
            let curve = Curve::OnSurface(Box::new(CurveOnSurface::new(
                pcurves[0].clone(),
                chart(frame),
            )));
            assert!(exact_pcurve_of(&curve, &ground, T).is_none());
        }
    }

    /// The plane's chart axes in space.
    fn ground_axes(ground: &SurfaceGeometry) -> (Vector, Vector) {
        let SurfaceGeometry::Plane(p) = ground else {
            unreachable!("a plane")
        };
        let frame = p.plane().frame();
        (frame.x().vector(), frame.y().vector())
    }

    /// The cone `ρ = 3.2 − z / 10` about the `z` axis, down from `z = 2`.
    fn drafted_drum() -> SurfaceGeometry {
        let frame = Frame::new(Point::new(0.0, 0.0, 2.0), -Direction::Z, -Direction::X, T).unwrap();
        ogeom_geom::ConeSurface::new(
            ogeom_math::Cone::new(frame, 3.0, 0.1_f64.atan(), T).unwrap(),
            (0.0, 2.0),
        )
        .unwrap()
        .into()
    }

    /// The rim of [`drafted_drum`] rounded at radius ½ under `z = 2`, its
    /// section's arc in the plane `y = 0` from its contact with the cone
    /// over `turn`, swept along `y`, the rim's tangent there, over `sweep`:
    /// a band of the drum of radius ½ about `y` through `(c, 0, 3 / 2)`,
    /// `c = 3.05 − √1.01 / 2`. Returns the band and `c`.
    fn carried_round(turn: f64, sweep: (f64, f64)) -> (SurfaceGeometry, f64) {
        let lean = 1.01_f64.sqrt();
        let c = lean.mul_add(-0.5, 3.05);
        let frame = Frame::from_axes(
            Point::new(c, 0.0, 1.5),
            Direction::new(Vector::new(1.0, 0.0, 0.1), T).unwrap(),
            Direction::new(Vector::new(-0.1, 0.0, 1.0), T).unwrap(),
            -Direction::Y,
            T,
        )
        .unwrap();
        let arc: Curve = ogeom_geom::TrimmedCurve::new(
            ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(frame, 0.5, T).unwrap()).into(),
            0.0,
            turn,
            T,
        )
        .unwrap()
        .into();
        let band = ogeom_geom::ExtrusionSurface::over(arc, Direction::Y, sweep)
            .unwrap()
            .into();
        (band, c)
    }

    /// [`carried_round`] over `turn` and `y` over `sweep`, turned `angle`
    /// about the cone's axis, against the cone in both orders, the band
    /// walked within `window` where one is given. The arc touches the cone
    /// at its start, where the sweep's ruling is the rim's tangent: the band
    /// stands off the cone by about `0.164 y² − u² / 4` near there, so two
    /// branches leave that node on the band's border, one to each end of
    /// the sweep or the window. Each is the point of the band at the turn
    /// `u` solving `√(x(u)² + y²) = 3.2 − z(u) / 10` at each `y`, turned
    /// with the band. Every section must be a crossing on both surfaces,
    /// stating no wider a tolerance than the chord and the fit asked for,
    /// the two branches must end at the node and at both ends, and each end
    /// must end one section only: a branch returned twice ends there twice.
    /// No section may end anywhere else.
    fn crosses_from_the_node_to_both_ends(
        (turn, sweep): (f64, (f64, f64)),
        angle: f64,
        window: Option<ParameterWindow>,
    ) {
        let (band, c) = carried_round(turn, sweep);
        let about = ogeom_math::Transform::rotation(
            ogeom_math::Axis {
                location: Point::ORIGIN,
                direction: Direction::Z,
            },
            angle,
        );
        let band = ogeom_geom::Transformable::transformed(&band, &about, T).unwrap();
        let (sin, cos) = angle.sin_cos();
        let turned = |p: Point| {
            Point::new(
                p.x.mul_add(cos, -p.y * sin),
                p.x.mul_add(sin, p.y * cos),
                p.z,
            )
        };
        let back = |p: Point| {
            Point::new(
                p.x.mul_add(cos, p.y * sin),
                p.y.mul_add(cos, -p.x * sin),
                p.z,
            )
        };
        let drum = drafted_drum();
        let lean = 1.01_f64.sqrt();
        let at = |u: f64, y: f64| {
            let (cos, sin) = (u.cos(), u.sin());
            Point::new(
                0.5f64.mul_add((cos - 0.1 * sin) / lean, c),
                y,
                0.5f64.mul_add((0.1f64.mul_add(cos, sin)) / lean, 1.5),
            )
        };
        let off = |p: Point| p.x.hypot(p.y) - 0.1f64.mul_add(-p.z, 3.2);
        let crossing = |y: f64| {
            let (mut lo, mut hi) = (1e-9, 0.425);
            for _ in 0..200 {
                let mid = f64::midpoint(lo, hi);
                if off(at(mid, y)).signum() == off(at(lo, y)).signum() {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            turned(at(f64::midpoint(lo, hi), y))
        };
        let node = turned(at(0.0, 0.0));
        let ends = window.map_or(sweep, |w| w.1);
        let far = [crossing(ends.0), crossing(ends.1)];
        // A boolean marches at a chord of a ten-millionth of the parts.
        let options = window.map_or_else(IntersectOptions::default, |_| IntersectOptions {
            tolerance: 1e-7,
            marching: Marching {
                chord: 1e-7,
                ..Marching::default()
            },
        });
        for band_first in [true, false] {
            let met = if band_first {
                intersect_surfaces_within(&band, &drum, (window, None), options, T)
            } else {
                intersect_surfaces_within(&drum, &band, (None, window), options, T)
            }
            .unwrap();
            let SurfaceIntersection::Along(sections) = met else {
                panic!("turned {angle}: the band crosses the cone, got {met:?}");
            };
            let mut found: Vec<Point> = Vec::new();
            for section in &sections {
                assert!(
                    !section.tangential,
                    "turned {angle}, band first {band_first}: {section:?}"
                );
                assert!(
                    section.tolerance <= options.marching.chord + options.tolerance,
                    "turned {angle}, band first {band_first}: a tolerance of {:e}",
                    section.tolerance
                );
                let (lo, hi) = section.curve.domain();
                let near = section.tolerance + T.confusion();
                for i in 0..=16 {
                    let p = section
                        .curve
                        .point_at(lo + (hi - lo) * f64::from(i) / 16.0, T)
                        .unwrap();
                    let q = back(p);
                    let from_band = ((q.x - c).hypot(q.z - 1.5) - 0.5).abs();
                    assert!(
                        from_band <= near && off(p).abs() <= near,
                        "{p:?} is off the band by {from_band:e} or the cone by {:e}",
                        off(p)
                    );
                }
                found.extend([lo, hi].map(|t| section.curve.point_at(t, T).unwrap()));
            }
            let near = sections.iter().map(|s| s.tolerance).fold(0.0, f64::max) + T.confusion();
            assert!(
                found.iter().any(|p| p.distance(node) <= near),
                "turned {angle}: no section ends at the node {node:?}: {found:?}"
            );
            for want in far {
                let count = found.iter().filter(|p| p.distance(want) <= near).count();
                assert_eq!(
                    count, 1,
                    "turned {angle}: {count} sections end at {want:?}: {found:?}"
                );
            }
            for p in &found {
                assert!(
                    p.distance(node) <= near || far.iter().any(|q| p.distance(*q) <= near),
                    "turned {angle}, band first {band_first}: a section ends at {p:?}, \
                     neither the node nor an end: {found:?}"
                );
            }
        }
    }

    /// [`crosses_from_the_node_to_both_ends`] with `y` from −0.012 to 0.02:
    /// the two surfaces' normals part by a sine of at most 0.018 along the
    /// branches, under the marcher's tangency threshold. Both are sections
    /// the two cross along, not a contact: they were traced, read as
    /// tangential by their sine alone, and the pair came back apart.
    #[test]
    fn a_band_crossing_a_cone_shallowly_from_a_node_on_its_border_crosses_it() {
        crosses_from_the_node_to_both_ends((0.425, (-0.012, 0.02)), 0.0, None);
    }

    /// [`crosses_from_the_node_to_both_ends`] with the branch to the lower
    /// end of the sweep 0.0021 and 0.0043 long, within the sixty chords
    /// round the node where the march's samples count as on it, over a
    /// short turn and a long one. That branch was all samples on the node
    /// and went; where the walk from it stalled onto the other branch, the
    /// arm it leaves the node along was walked out only from past the
    /// node's radius, outside the band. The section stopped at the node.
    #[test]
    fn a_band_crossing_a_cone_from_a_node_keeps_its_arm_short_of_the_node_s_radius() {
        for turn in [0.425, 1.2] {
            for low in [-0.004, -0.002] {
                crosses_from_the_node_to_both_ends((turn, (low, 0.02)), 0.0, None);
            }
        }
    }

    /// [`crosses_from_the_node_to_both_ends`] with the branch to the lower
    /// end of the sweep 0.0008 and 0.0005 long, under the ten chords an arc
    /// of the march is debris within, over a short turn and a long one.
    /// That branch, the section out to the band's border, was dropped as
    /// debris; kept, every sample of it lies where the surfaces part by
    /// less than they measurably cross at, and it was read as a contact.
    #[test]
    fn a_band_crossing_a_cone_from_a_node_keeps_its_arm_shorter_than_ten_chords() {
        for turn in [0.425, 1.2] {
            for low in [-0.000_74, -0.000_5] {
                crosses_from_the_node_to_both_ends((turn, (low, 0.02)), 0.0, None);
            }
        }
    }

    /// [`crosses_from_the_node_to_both_ends`] over the band a fillet
    /// junction carries on past its cap at a drum drafted foot out:
    /// the round's whole quarter turn swept along `y` from −2.42 to 4.83,
    /// walked within the window of the face a boolean cuts from it, `u` to
    /// 0.42 and `y` from −0.004 to 0.02, turned about the cone's axis every
    /// half degree. Both branches leave that window through its ends in
    /// `y`, crossing the cone there. At the odd multiples of 7.5 degrees the
    /// cells the march seeds from found neither, and the pair came back
    /// apart. Turned 60 degrees, the cone named first, the lower branch was
    /// traced twice, the two tracings standing microns apart across it
    /// where it nears the node: the stretch of one there was kept beside
    /// the other, read as a contact, and the valley walker followed the
    /// band's border, tangent to the cone at the node, 0.0078 either way.
    #[test]
    fn a_band_crossing_a_cone_is_two_crossings_at_every_turn() {
        let window = (
            (0.0, 0.424_721_486_693_660_6),
            (-0.003_979_668_916_835_106, 0.019_898_344_584_175_24),
        );
        for k in 0..720 {
            crosses_from_the_node_to_both_ends(
                (
                    1.471_127_674_303_734_5,
                    (-2.415_259_779_545_657_6, 4.830_519_559_091_315),
                ),
                (f64::from(k) * 0.5).to_radians(),
                Some(window),
            );
        }
    }

    /// [`crosses_from_the_node_to_both_ends`] over the band of
    /// [`a_band_crossing_a_cone_is_two_crossings_at_every_turn`], walked
    /// within a window cut `1e-5` below the node, turned every half degree:
    /// the lower branch is an arm `1.1e-5` long, out to the window's lower
    /// border. Inside the node's snap, `√(2ε / bend)`, about `3.2e-4` here,
    /// every point stands within `ε` of both surfaces, and the walk running
    /// through the node landed on that border `1.7e-4` along it from where
    /// the arm crosses it: at a dozen turns that came back as a third
    /// section beside the arm, ending where no section does.
    #[test]
    fn a_band_crossing_a_cone_in_a_window_cut_by_its_node_has_no_section_along_the_cut() {
        let window = (
            (0.0, 0.424_721_486_693_660_6),
            (-1e-5, 0.019_898_344_584_175_24),
        );
        for k in 0..720 {
            crosses_from_the_node_to_both_ends(
                (
                    1.471_127_674_303_734_5,
                    (-2.415_259_779_545_657_6, 4.830_519_559_091_315),
                ),
                (f64::from(k) * 0.5).to_radians(),
                Some(window),
            );
        }
    }

    /// [`crosses_from_the_node_to_both_ends`] with the branch to the lower
    /// end of the sweep about `3.4e-4`, `2.2e-4` and `1.2e-4` long, at the
    /// default chord of `1e-4`, turned every degree: each arm leaves the
    /// band at or inside the node's snap, `3.2e-4`, where no walk is set
    /// out from, and was missing.
    #[test]
    fn a_band_crossing_a_cone_from_a_node_keeps_its_arm_shorter_than_the_snap() {
        for low in [-0.000_32, -0.000_2, -0.000_11] {
            for k in 0..360 {
                crosses_from_the_node_to_both_ends(
                    (0.425, (low, 0.02)),
                    f64::from(k).to_radians(),
                    None,
                );
            }
        }
    }

    /// [`crosses_from_the_node_to_both_ends`] turned 60 degrees, the cone
    /// named first: the lower crossing's samples by the node stood across
    /// the arm by up to `9e-6`, each within `3e-10` of both surfaces, and
    /// its fit stated a tolerance of `6.5e-4` at a chord of `1e-7`. Every
    /// half degree is weighed by
    /// [`a_band_crossing_a_cone_is_two_crossings_at_every_turn`].
    #[test]
    fn a_band_crossing_a_cone_states_a_tight_tolerance_by_its_node() {
        let window = (
            (0.0, 0.424_721_486_693_660_6),
            (-0.003_979_668_916_835_106, 0.019_898_344_584_175_24),
        );
        crosses_from_the_node_to_both_ends(
            (
                1.471_127_674_303_734_5,
                (-2.415_259_779_545_657_6, 4.830_519_559_091_315),
            ),
            60.0_f64.to_radians(),
            Some(window),
        );
    }

    /// A drum of radius ½ laid along a ruling of [`drafted_drum`], outside
    /// it, its axis that far off the ruling along the cone's normal, less
    /// `into`: the cone's tangent plane holds along the ruling, so resting
    /// there the two touch along it, and pressed in by `ε` they cross along
    /// two curves either side of it with `ε` between them. Resting, or
    /// pressed in by less than the confusion distance a marched point
    /// stands within, they are a contact; pressed in by a micron or ten,
    /// where the two cross at a sine of 0.002 or 0.007, they cross along
    /// two sections, one each side of the ruling, every point on both. The
    /// sine alone read the crossings as a contact.
    #[test]
    fn a_drum_pressed_into_a_cone_along_a_ruling_crosses_it_only_measurably() {
        let drum = drafted_drum();
        let lean = 1.01_f64.sqrt();
        let up = Vector::new(-0.1, 0.0, 1.0) / lean;
        let out = Vector::new(1.0, 0.0, 0.1) / lean;
        let off = |p: Point| p.x.hypot(p.y) - 0.1f64.mul_add(-p.z, 3.2);
        for into in [0.0, 0.4 * T.confusion(), 1e-6, 1e-5] {
            let foot = Point::new(3.0, 0.0, 2.0) + out * (0.5 - into);
            let frame = Frame::new(
                foot,
                Direction::new(up, T).unwrap(),
                Direction::new(out, T).unwrap(),
                T,
            )
            .unwrap();
            let roller: SurfaceGeometry =
                CylinderSurface::new(Cylinder::new(frame, 0.5, T).unwrap(), (-1.8, -0.2))
                    .unwrap()
                    .into();
            let from_axis = |p: Point| {
                let d = p - foot;
                (d - up * d.dot(up)).magnitude()
            };
            for (a, b) in [(&roller, &drum), (&drum, &roller)] {
                let met = intersect_surfaces(a, b, IntersectOptions::default(), T).unwrap();
                let SurfaceIntersection::Along(sections) = met else {
                    panic!("pressed in {into:e}: {met:?}");
                };
                if into <= T.confusion() {
                    assert!(
                        sections.iter().all(|s| s.tangential),
                        "pressed in {into:e}: {sections:?}"
                    );
                    continue;
                }
                assert_eq!(sections.len(), 2, "pressed in {into:e}: {sections:?}");
                let mut sides = Vec::new();
                for section in &sections {
                    assert!(!section.tangential, "pressed in {into:e}: {section:?}");
                    let (lo, hi) = section.curve.domain();
                    let near = section.tolerance + T.confusion();
                    let mut side = 0.0;
                    for i in 0..=16 {
                        let p = section
                            .curve
                            .point_at(lo + (hi - lo) * f64::from(i) / 16.0, T)
                            .unwrap();
                        assert!(
                            (from_axis(p) - 0.5).abs() <= near && off(p).abs() <= near,
                            "pressed in {into:e}: {p:?} is off the drum or the cone"
                        );
                        side += p.y;
                    }
                    sides.push(side);
                }
                assert!(
                    sides[0] * sides[1] < 0.0,
                    "pressed in {into:e}: both on one side {sides:?}"
                );
            }
        }
    }
}
