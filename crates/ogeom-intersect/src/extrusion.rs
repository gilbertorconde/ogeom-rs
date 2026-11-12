//! A plane across a linear extrusion, in closed form.
//!
//! An extrusion sweeps a curve `C` along a unit direction `d`: its point at
//! `(u, v)` is `C(u) + v d`. A plane through `o` with unit normal `n` holds
//! that point where `v = n·(o − C(u)) / n·d`, so where `n·d` is not zero the
//! section is `C` carried along `d` onto the plane, at `C`'s own parameter:
//! the image of `C` under an affine map. An affine map takes a line to a
//! line, a circle or an ellipse to an ellipse, a parabola to a parabola,
//! and a B-spline, polynomial or rational, to the B-spline of its control
//! points carried, on its own knots and weights. Where `n·d` is zero the
//! plane holds `d`, and the section is the rulings through `C`'s crossings
//! with the plane.
//!
//! Each section is cut to the extrusion's own window: to the stretches of
//! `C`'s domain where the height `v` stays within the sweep, found as `C`'s
//! crossings with the plane moved to each end of the sweep, each a root of
//! a closed form. A B-spline's crossings have none, and its section is
//! taken only where its control points' heights all lie within the sweep,
//! which holds every point's. Where they straddle an end,
//! [`plane_extrusion_at_rims`] splits it at its crossings with the rim
//! there, roots of a polynomial B-spline found to rounding, for a caller
//! that asks.
//!
//! The pair is left to the marcher, as before, wherever the closed form
//! does not hold cleanly: a profile this does not map (a helix, an offset,
//! a hyperbola, whose chart trace on the extrusion no planar curve states,
//! a reversed trim, a periodic B-spline), a B-spline whose heights straddle
//! an end of the sweep or that a plane holding `d` cuts, a profile lying in
//! a plane it is cut by, a crossing it only touches, a plane or a ruling
//! the extrusion meets within the sine the marcher takes as tangent, or a
//! section that does not land on both surfaces within the confusion
//! distance when checked.
//!
//! A cylinder whose axis runs along the sweep meets an extrusion of a
//! B-spline profile in the rulings through the profile's crossings with
//! it, found the same way ([`drum_extrusion`]).

use ogeom_core::Tolerances;
use ogeom_geom::{
    BSpline2d, BSplineCurve, Curve, Curve2d as _, Curve3d as _, LineCurve, PlanarCurve, Surface,
    SurfaceGeometry,
};
use ogeom_math::bspline::{self, Spline, Weighted};
use ogeom_math::{
    Axis, Axis2, Direction, Direction2, Frame, KnotVector, Point, Point2, Vector, Vector2,
};

use crate::march::TANGENTIAL_SINE;
use crate::section::{SectionCurve, exact_pcurve_over, inside_box, touches_box};

/// Where a plane meets an extrusion, exactly: `None` where the pair is not
/// a plane and an extrusion, or the closed form does not hold cleanly (see
/// the module's notes). An empty answer is a plane missing the extrusion's
/// window.
pub(crate) fn plane_extrusion(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<Vec<SectionCurve>> {
    plane_extrusion_with(a, b, false, tol)
}

/// Where a plane meets an extrusion of a B-spline profile whose heights to
/// the plane straddle an end of the sweep, the section split where it
/// crosses the rim there: `None` for any other pair, which
/// [`intersect_surfaces`](crate::intersect_surfaces) answers as before,
/// and where a crossing is not resolved cleanly.
///
/// The profile's crossings with the plane moved to the sweep's end are
/// the roots of a polynomial B-spline, isolated on its Bernstein
/// coefficients and bisected to rounding; each stretch of the profile whose height stays
/// within the sweep is carried onto the plane in closed form and ends
/// exactly on the rim, and a stretch past the end, where the section
/// leaves the extrusion, is dropped.
#[must_use]
pub fn plane_extrusion_at_rims(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<Vec<SectionCurve>> {
    if !crosses_extrusion_rim(a, b, tol) {
        return None;
    }
    plane_extrusion_with(a, b, true, tol)
}

/// Whether `a` and `b` are a plane and an extrusion of a B-spline profile
/// whose control points' heights to the plane straddle an end of the
/// sweep: the pair [`plane_extrusion_at_rims`] answers, and
/// [`intersect_surfaces`](crate::intersect_surfaces) leaves to the
/// marcher.
#[must_use]
pub fn crosses_extrusion_rim(a: &SurfaceGeometry, b: &SurfaceGeometry, tol: Tolerances) -> bool {
    let Some((carry, spline, sweep)) = spline_pair(a, b, tol) else {
        return false;
    };
    let heights = spline.heights(&carry);
    [sweep.0, sweep.1]
        .iter()
        .any(|&end| heights.iter().any(|&v| v < end) && heights.iter().any(|&v| v > end))
}

/// The carry of a plane across an extrusion of a B-spline profile, the
/// profile and the sweep, where the pair is one and the plane crosses the
/// sweep.
fn spline_pair(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<(Carry, SplineProfile, (f64, f64))> {
    let (p, extrusion) = match (a, b) {
        (SurfaceGeometry::Plane(p), SurfaceGeometry::Extrusion(e))
        | (SurfaceGeometry::Extrusion(e), SurfaceGeometry::Plane(p)) => (p, &**e),
        _ => return None,
    };
    let plane = p.plane();
    let (n, d) = (plane.normal().vector(), extrusion.direction().vector());
    let along = n.dot(d);
    if along.abs() < TANGENTIAL_SINE {
        return None;
    }
    let (window, sweep) = extrusion.domain();
    let spline = SplineProfile::of(extrusion.curve(), window, tol)?;
    let carry = Carry {
        n,
        o: plane.origin(),
        d,
        along,
    };
    Some((carry, spline, sweep))
}

/// Where a cylinder meets an extrusion of a B-spline profile swept along
/// its axis, exactly: the rulings through the profile's crossings with the
/// cylinder, each over the sweep and cut to the cylinder's window. `None`
/// where the pair is not one, the sweep leans off the axis by more than
/// the confusion distance over its reach, or a crossing is not resolved
/// cleanly; an empty answer is a profile missing the cylinder.
///
/// A point `P` of the profile lies on the cylinder where `X² + Y² = r²`,
/// `X` and `Y` its coordinates across the axis in the cylinder's frame.
/// With the profile's weights `w`, `w²(X² + Y² − r²)` is a polynomial
/// B-spline of twice the degree whose sign is the condition's, its
/// Bernstein coefficients over each Bézier piece the products of the
/// piece's `wX`, `wY` and `w` (see [`level_crossings`] for how its roots
/// are isolated). A root isolated that way is simple, so the profile
/// crosses the cylinder there however shallow the angle, as a line
/// crossing it does in [`surface_surface`](crate::surface_surface)'s
/// closed form.
pub(crate) fn drum_extrusion(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<Vec<SectionCurve>> {
    let (drum_surface, extrusion, drum_first) = match (a, b) {
        (SurfaceGeometry::Cylinder(_), SurfaceGeometry::Extrusion(e)) => (a, &**e, true),
        (SurfaceGeometry::Extrusion(e), SurfaceGeometry::Cylinder(_)) => (b, &**e, false),
        _ => return None,
    };
    let SurfaceGeometry::Cylinder(c) = drum_surface else {
        return None;
    };
    let drum = c.cylinder();
    let frame = drum.frame();
    let (o, r) = (frame.origin(), drum.radius());
    let (x, y) = (frame.x().vector(), frame.y().vector());
    let d = extrusion.direction();
    let (window, sweep) = extrusion.domain();
    let reach = sweep.0.abs().max(sweep.1.abs());
    if d.vector().cross(frame.z().vector()).magnitude() * reach > tol.confusion() {
        return None;
    }
    let spline = SplineProfile::of(extrusion.curve(), window, tol)?;
    // An end of the profile on the cylinder is the extrusion's border
    // lying in it, its sign there rounding's: left to the marcher.
    let ends = (spline.control.first()?, spline.control.last()?);
    if [ends.0, ends.1]
        .iter()
        .any(|w| drum.distance_to(w.point()) <= tol.confusion())
    {
        return None;
    }
    let across: Vec<Weighted<Point2>> = spline
        .control
        .iter()
        .map(|w| {
            let q = w.point() - o;
            Weighted {
                scaled: Point2::new(w.weight * q.dot(x), w.weight * q.dot(y)),
                weight: w.weight,
            }
        })
        .collect();
    // Each coefficient is a blend of products of the control points'
    // coordinates across the axis, each rounded to their size and to the
    // origin's.
    let noise = spline
        .control
        .iter()
        .map(|w| {
            let size = w.point().to_vector().magnitude() + o.to_vector().magnitude() + r;
            w.weight * w.weight * size * size
        })
        .fold(0.0, f64::max)
        * f64::EPSILON
        * 64.0;
    let pieces = bspline::to_bezier_segments(&spline.knots, &across, tol).ok()?;
    let mut roots = Vec::new();
    for ((lo, hi), piece) in pieces {
        let sx: Vec<f64> = piece.iter().map(|w| w.scaled.x).collect();
        let sy: Vec<f64> = piece.iter().map(|w| w.scaled.y).collect();
        let sw: Vec<f64> = piece.iter().map(|w| w.weight).collect();
        let (xx, yy, ww) = (
            bernstein_product(&sx, &sx),
            bernstein_product(&sy, &sy),
            bernstein_product(&sw, &sw),
        );
        let values: Vec<f64> = (0..xx.len())
            .map(|k| (-r * r).mul_add(ww[k], xx[k] + yy[k]))
            .collect();
        bezier_roots(&values, (lo, hi), noise, tol, &mut roots)?;
    }
    roots.sort_by(f64::total_cmp);
    roots.dedup_by(|p, q| (*p - *q).abs() <= tol.parametric());
    let mut out = Vec::with_capacity(roots.len());
    for u in roots {
        let at = bspline::evaluate_rational(&spline.knots, &spline.control, u, tol).ok()?;
        let axis = Axis {
            location: at,
            direction: d,
        };
        let line: Curve = LineCurve::over(axis, sweep.0, sweep.1).ok()?.into();
        let on_drum = exact_pcurve_over(&line, sweep, drum_surface, tol)?;
        let Some((lo, hi)) = inside_box(&on_drum, drum_surface) else {
            continue;
        };
        let (lo, hi) = (lo.max(sweep.0), hi.min(sweep.1));
        if hi - lo <= tol.parametric() {
            continue;
        }
        let line: Curve = LineCurve::over(axis, lo, hi).ok()?.into();
        let on_drum = exact_pcurve_over(&line, (lo, hi), drum_surface, tol)?;
        let on_extrusion: PlanarCurve =
            ogeom_geom::Line2d::over(Axis2::new(Point2::new(u, 0.0), Direction2::Y), lo, hi)
                .ok()?
                .into();
        lands_on_both(&line, &on_drum, &on_extrusion, drum_surface, extrusion, tol)?;
        let (on_a, on_b) = if drum_first {
            (on_drum, on_extrusion)
        } else {
            (on_extrusion, on_drum)
        };
        out.push(SectionCurve {
            curve: line,
            on_a: Some(on_a),
            on_b: Some(on_b),
            tolerance: 0.0,
            exact: true,
            closed: false,
            tangential: false,
        });
    }
    Some(out)
}

/// The Bernstein coefficients of the product of two polynomials given by
/// theirs over the same interval: `cₖ = Σ C(m, i) C(n, k − i) aᵢ bₖ₋ᵢ /
/// C(m + n, k)`, of degree `m + n`.
fn bernstein_product(a: &[f64], b: &[f64]) -> Vec<f64> {
    let (m, n) = (a.len() - 1, b.len() - 1);
    let choose = |top: usize, k: usize| -> f64 {
        (0..k).fold(1.0, |acc, j| {
            #[allow(clippy::cast_precision_loss, reason = "small degrees")]
            let (t, j) = ((top - j) as f64, (j + 1) as f64);
            acc * t / j
        })
    };
    (0..=m + n)
        .map(|k| {
            let terms: f64 = (k.saturating_sub(n)..=k.min(m))
                .map(|i| choose(m, i) * choose(n, k - i) * a[i] * b[k - i])
                .sum();
            terms / choose(m + n, k)
        })
        .collect()
}

/// [`plane_extrusion`], with a B-spline profile straddling an end of the
/// sweep split at the rim there where `at_rims` holds.
fn plane_extrusion_with(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    at_rims: bool,
    tol: Tolerances,
) -> Option<Vec<SectionCurve>> {
    let (plane_surface, extrusion, plane_first) = match (a, b) {
        (SurfaceGeometry::Plane(_), SurfaceGeometry::Extrusion(e)) => (a, &**e, true),
        (SurfaceGeometry::Extrusion(e), SurfaceGeometry::Plane(_)) => (b, &**e, false),
        _ => return None,
    };
    let SurfaceGeometry::Plane(p) = plane_surface else {
        return None;
    };
    let plane = p.plane();
    let n = plane.normal().vector();
    let d = extrusion.direction().vector();
    let o = plane.origin();
    let (window, sweep) = extrusion.domain();
    let along = n.dot(d);
    let sections = if let Some(spline) = SplineProfile::of(extrusion.curve(), window, tol) {
        // A spline profile's crossings with a plane holding `d` have no
        // closed form, so its rulings are not taken here.
        if along.abs() < TANGENTIAL_SINE {
            return None;
        }
        spline.carried(&Carry { n, o, d, along }, sweep, at_rims, tol)?
    } else {
        let (profile, window, whole_turn) = Profile::of(extrusion.curve(), window, tol)?;
        planar_sections(&profile, window, whole_turn, (n, o, d), extrusion, tol)?
    };
    finish(sections, plane_surface, extrusion, plane_first, tol)
}

/// The sections of an extrusion of a profile in closed form by the plane
/// through `o` with normal `n`, the extrusion swept along `d`.
fn planar_sections(
    profile: &Profile,
    window: (f64, f64),
    whole_turn: bool,
    (n, o, d): (Vector, Point, Vector),
    extrusion: &ogeom_geom::ExtrusionSurface,
    tol: Tolerances,
) -> Option<Vec<(Curve, PlanarCurve)>> {
    let (_, sweep) = extrusion.domain();
    let along = n.dot(d);
    let reach = sweep.0.abs().max(sweep.1.abs());
    Some(if along.abs() * reach <= tol.confusion() {
        rulings(profile, window, whole_turn, (n, o), extrusion, tol)?
    } else {
        // The extrusion's normal is square to `d`, so the sine between it
        // and `n` is at least `|n·d|`: past the tangency the marcher and
        // the contact handling answer, the two cross at an angle all
        // along every section.
        if along.abs() < TANGENTIAL_SINE {
            return None;
        }
        let carry = Carry { n, o, d, along };
        carried(profile, window, whole_turn, &carry, sweep, tol)?
    })
}

/// The pieces checked against both surfaces and stated as sections, each
/// trace on the surface asked first first.
fn finish(
    sections: Vec<(Curve, PlanarCurve)>,
    plane_surface: &SurfaceGeometry,
    extrusion: &ogeom_geom::ExtrusionSurface,
    plane_first: bool,
    tol: Tolerances,
) -> Option<Vec<SectionCurve>> {
    // Each piece checked against both surfaces, kept where it reaches the
    // plane's window, and cut to it where it is straight.
    let mut out = Vec::with_capacity(sections.len());
    for (curve, on_extrusion) in sections {
        let Some((curve, on_extrusion, on_plane)) =
            within_plane_window(curve, on_extrusion, plane_surface, tol)?
        else {
            continue;
        };
        lands_on_both(
            &curve,
            &on_plane,
            &on_extrusion,
            plane_surface,
            extrusion,
            tol,
        )?;
        let closed = curve.is_closed(tol);
        let (on_a, on_b) = if plane_first {
            (on_plane, on_extrusion)
        } else {
            (on_extrusion, on_plane)
        };
        out.push(SectionCurve {
            curve,
            on_a: Some(on_a),
            on_b: Some(on_b),
            tolerance: 0.0,
            exact: true,
            closed,
            tangential: false,
        });
    }
    Some(out)
}

/// The profile in the forms the map acts on, each at the profile's own
/// parameter `u`.
enum Profile {
    /// `p + u t`, `t` of unit length.
    Line { p: Point, t: Vector },
    /// `c + cos u · x + sin u · y`: a circle or an ellipse.
    Conic { c: Point, x: Vector, y: Vector },
    /// `apex + u² · x + u · y`.
    Parabola { apex: Point, x: Vector, y: Vector },
}

impl Profile {
    /// The profile of `curve` over `window`, the window, and whether that
    /// is a whole turn of a periodic profile. `None` for a profile this
    /// does not map.
    fn of(curve: &Curve, window: (f64, f64), tol: Tolerances) -> Option<(Self, (f64, f64), bool)> {
        let (basis, whole) = match curve {
            Curve::Trimmed(t) if !t.is_reversed() => (t.basis(), false),
            Curve::Trimmed(_) => return None,
            other => (other, other.is_periodic()),
        };
        let profile = match basis {
            Curve::Line(line) => {
                let axis = line.axis();
                Self::Line {
                    p: axis.location,
                    t: axis.direction.vector(),
                }
            }
            Curve::Circle(c) => {
                let circle = c.forward(tol).ok()?.circle();
                let frame = circle.frame();
                Self::Conic {
                    c: circle.centre(),
                    x: frame.x().vector() * circle.radius(),
                    y: frame.y().vector() * circle.radius(),
                }
            }
            Curve::Ellipse(e) => {
                let ellipse = e.forward(tol).ok()?.ellipse();
                let frame = ellipse.frame();
                Self::Conic {
                    c: ellipse.centre(),
                    x: frame.x().vector() * ellipse.major_radius(),
                    y: frame.y().vector() * ellipse.minor_radius(),
                }
            }
            Curve::Parabola(parabola) if !parabola.is_reversed() => {
                let parabola = parabola.parabola();
                let frame = parabola.frame();
                Self::Parabola {
                    apex: parabola.apex(),
                    x: frame.x().vector() / (4.0 * parabola.focal()),
                    y: frame.y().vector(),
                }
            }
            _ => return None,
        };
        Some((profile, window, whole))
    }

    /// The point at `u`.
    fn point(&self, u: f64) -> Point {
        match self {
            Self::Line { p, t } => *p + *t * u,
            Self::Conic { c, x, y } => *c + *x * u.cos() + *y * u.sin(),
            Self::Parabola { apex, x, y } => *apex + *x * (u * u) + *y * u,
        }
    }

    /// The first derivative at `u`.
    fn tangent(&self, u: f64) -> Vector {
        match self {
            Self::Line { t, .. } => *t,
            Self::Conic { x, y, .. } => *y * u.cos() - *x * u.sin(),
            Self::Parabola { x, y, .. } => *x * (2.0 * u) + *y,
        }
    }

    /// Where `n·C(u) = k` over `window`, in order. `None` where the
    /// profile lies within the confusion distance of that plane along a
    /// stretch, or touches it, or two crossings stand closer than the
    /// parametric confusion: there the plane is no clean crossing.
    fn crossings(
        &self,
        n: Vector,
        k: f64,
        (lo, hi): (f64, f64),
        tol: Tolerances,
    ) -> Option<Vec<f64>> {
        let near = tol.confusion();
        let mut found = Vec::new();
        match self {
            Self::Line { p, t } => {
                let (rate, at) = (n.dot(*t), n.dot(p.to_vector()) - k);
                if rate.abs() * (hi - lo) <= near {
                    if rate.mul_add(lo, at).abs() <= near {
                        return None;
                    }
                } else {
                    found.push(-at / rate);
                }
            }
            Self::Conic { c, x, y } => {
                // `A + B cos u + C sin u`: `R cos(u − φ) = −A`.
                let (a, b, s) = (n.dot(c.to_vector()) - k, n.dot(*x), n.dot(*y));
                let r = b.hypot(s);
                if (a.abs() - r).abs() <= near {
                    return None;
                }
                if a.abs() < r {
                    let phase = s.atan2(b);
                    let half = (-a / r).acos();
                    let tau = core::f64::consts::TAU;
                    for root in [phase - half, phase + half] {
                        let mut u = lo + (root - lo).rem_euclid(tau);
                        while u <= hi + tol.parametric() {
                            found.push(u);
                            u += tau;
                        }
                    }
                }
            }
            Self::Parabola { apex, x, y } => {
                let (a, b, c) = (n.dot(*x), n.dot(*y), n.dot(apex.to_vector()) - k);
                if a.abs() > 0.0 {
                    let turn = -b / (2.0 * a);
                    if (lo..=hi).contains(&turn)
                        && a.mul_add(turn * turn, b.mul_add(turn, c)).abs() <= near
                    {
                        return None;
                    }
                } else if b.abs() * (hi - lo) <= near && b.mul_add(lo, c).abs() <= near {
                    return None;
                }
                found.extend(ogeom_math::solve::quadratic_roots(a, b, c));
            }
        }
        let slack = tol.parametric();
        let mut kept: Vec<f64> = found
            .into_iter()
            .filter(|u| u.is_finite() && *u >= lo - slack && *u <= hi + slack)
            .map(|u| u.clamp(lo, hi))
            .collect();
        kept.sort_by(f64::total_cmp);
        kept.dedup_by(|x, y| (*x - *y).abs() <= slack);
        Some(kept)
    }
}

/// A clamped B-spline profile over the extrusion's window, polynomial or
/// rational: carried onto a plane, its control points carried with their
/// weights are the image's, since the carry is affine and the weights are
/// positive.
struct SplineProfile {
    knots: KnotVector,
    control: Vec<Weighted<Point>>,
}

impl SplineProfile {
    /// The profile of `curve` cut to `window`, where it is a B-spline, a
    /// forward trim of one, and not periodic.
    fn of(curve: &Curve, window: (f64, f64), tol: Tolerances) -> Option<Self> {
        let basis = match curve {
            Curve::Trimmed(t) if !t.is_reversed() => t.basis(),
            other => other,
        };
        let Curve::BSpline(spline) = basis else {
            return None;
        };
        if basis.is_periodic() || !spline.knots().is_clamped() {
            return None;
        }
        let whole: Spline<Weighted<Point>> =
            (spline.knots().clone(), spline.control_points().to_vec());
        let (knots, control) = cut_to(whole, window, tol)?;
        Some(Self { knots, control })
    }

    /// The profile carried onto the plane, with its chart trace on the
    /// extrusion: nothing where every control point's height, and so every
    /// point's, falls outside the sweep on one side. Where the heights
    /// straddle an end of the sweep, `None` unless `at_rims` holds, and
    /// then each stretch whose height stays within the sweep, split at its
    /// crossings with the rims (see [`level_crossings`]).
    fn carried(
        &self,
        carry: &Carry,
        sweep: (f64, f64),
        at_rims: bool,
        tol: Tolerances,
    ) -> Option<Vec<(Curve, PlanarCurve)>> {
        // A point of the profile is a convex blend of its control points,
        // and so is its height of theirs.
        let heights = self.heights(carry);
        if heights.iter().all(|&v| v < sweep.0) || heights.iter().all(|&v| v > sweep.1) {
            return Some(Vec::new());
        }
        if heights.iter().all(|&v| v >= sweep.0 && v <= sweep.1) {
            return Some(vec![self.carried_whole(carry, &heights, tol)?]);
        }
        if !at_rims {
            return None;
        }
        let window = self.knots.domain();
        let mut cuts = vec![window.0, window.1];
        for end in [sweep.0, sweep.1] {
            if heights.iter().any(|&v| v < end) && heights.iter().any(|&v| v > end) {
                cuts.extend(level_crossings(
                    &self.knots,
                    &self.control,
                    carry,
                    end,
                    tol,
                )?);
            }
        }
        cuts.sort_by(f64::total_cmp);
        cuts.dedup_by(|x, y| (*x - *y).abs() <= tol.parametric());
        let mut out = Vec::new();
        for pair in cuts.windows(2) {
            let (lo, hi) = (pair[0], pair[1]);
            if hi - lo <= tol.parametric() {
                continue;
            }
            let mid =
                bspline::evaluate_rational(&self.knots, &self.control, f64::midpoint(lo, hi), tol)
                    .ok()?;
            let v = carry.height(mid);
            if v < sweep.0 || v > sweep.1 {
                continue;
            }
            let (knots, control) =
                cut_to((self.knots.clone(), self.control.clone()), (lo, hi), tol)?;
            let piece = Self { knots, control };
            let heights = piece.heights(carry);
            out.push(piece.carried_whole(carry, &heights, tol)?);
        }
        Some(out)
    }

    /// Each control point's height to the plane.
    fn heights(&self, carry: &Carry) -> Vec<f64> {
        self.control
            .iter()
            .map(|w| carry.height(w.point()))
            .collect()
    }

    /// The whole profile carried onto the plane, with its chart trace.
    fn carried_whole(
        &self,
        carry: &Carry,
        heights: &[f64],
        tol: Tolerances,
    ) -> Option<(Curve, PlanarCurve)> {
        let image: Vec<Weighted<Point>> = self
            .control
            .iter()
            .map(|w| Weighted::new(carry.point(w.point()), w.weight, tol))
            .collect::<Result<_, _>>()
            .ok()?;
        let curve: Curve = BSplineCurve::rational(self.knots.clone(), image)
            .ok()?
            .into();
        let chart = self.chart(heights, tol)?;
        Some((curve, chart.into()))
    }

    /// The chart trace `(u, v(u))` at the profile's own parameter, where
    /// `v` is the height blended as the profile's points are. Over each
    /// Bézier piece on `[a, b]` with weights `wₖ`, `v` is the blend of the
    /// heights with those weights, and `u` is `u W(u) / W(u)`, `u W(u)` one
    /// degree up with coefficients `((p + 1 − k) wₖ a + k wₖ₋₁ b) / (p + 1)`.
    /// The pieces join on their shared ends.
    fn chart(&self, heights: &[f64], tol: Tolerances) -> Option<BSpline2d> {
        let p = self.knots.degree();
        let lifted: Vec<Weighted<Point2>> = self
            .control
            .iter()
            .zip(heights)
            .map(|(w, &v)| Weighted {
                scaled: Point2::new(0.0, v * w.weight),
                weight: w.weight,
            })
            .collect();
        let pieces = bspline::to_bezier_segments(&self.knots, &lifted, tol).ok()?;
        let up = p + 1;
        #[allow(clippy::cast_precision_loss)]
        let share = |k: usize| k as f64 / up as f64;
        let mut knots: Vec<f64> = Vec::with_capacity(pieces.len() * up + 2);
        let mut control: Vec<Weighted<Point2>> = Vec::with_capacity(pieces.len() * up + 1);
        for (index, ((a, b), piece)) in pieces.iter().enumerate() {
            for k in 0..=up {
                if index > 0 && k == 0 {
                    continue;
                }
                let (before, here) = ((k > 0).then(|| piece[k - 1]), (k <= p).then(|| piece[k]));
                let blend = |f: fn(&Weighted<Point2>) -> f64| {
                    before.map_or(0.0, |w| share(k) * f(&w))
                        + here.map_or(0.0, |w| (1.0 - share(k)) * f(&w))
                };
                let weight = blend(|w| w.weight);
                let rise = blend(|w| w.scaled.y);
                let run = (1.0 - share(k)) * here.map_or(0.0, |w| w.weight) * a
                    + share(k) * before.map_or(0.0, |w| w.weight) * b;
                control.push(Weighted {
                    scaled: Point2::new(run, rise),
                    weight,
                });
            }
            let ends = if index == 0 { up + 1 } else { up };
            knots.extend(core::iter::repeat_n(*a, ends));
            if index == pieces.len() - 1 {
                knots.extend(core::iter::repeat_n(*b, up + 1));
            }
        }
        BSpline2d::rational(KnotVector::new(knots, up).ok()?, control).ok()
    }
}

/// A clamped B-spline cut to `(lo, hi)` within its domain.
fn cut_to<P: bspline::Blend>(
    mut piece: Spline<P>,
    (lo, hi): (f64, f64),
    tol: Tolerances,
) -> Option<Spline<P>> {
    let (start, end) = piece.0.domain();
    if lo > start + tol.parametric() {
        piece = bspline::split(&piece.0, &piece.1, lo, tol).ok()?.1;
    }
    if hi < end - tol.parametric() {
        piece = bspline::split(&piece.0, &piece.1, hi, tol).ok()?.0;
    }
    Some(piece)
}

/// Where a B-spline profile's height to the plane is `level`, in order:
/// its crossings with the plane moved `level` along the sweep, which on
/// the extrusion is the rim the profile sweeps to there.
///
/// The height less `level` is `Σ Nₖ wₖ (vₖ − level) / Σ Nₖ wₖ`, so its sign
/// is that of a polynomial B-spline with coefficients `wₖ (vₖ − level)`, a
/// blend of the control points' heights with their weights. Over each
/// Bézier piece its values lie within its coefficients' hull, and the
/// number of the coefficients' sign changes bounds its roots from above
/// with the same parity (Bernstein's form of Descartes' rule). A piece
/// whose coefficients keep one sign holds no root, one with a single sign
/// change and its ends of opposite signs holds exactly one, bisected down
/// to rounding, and any other is halved by de Casteljau's construction.
/// `None` where a piece is halved down to the parametric confusion
/// without settling (a profile touching the level there, or meeting it at
/// a piece's joint), or where every coefficient of a piece is within the
/// rounding they were computed to (a profile lying in the level).
fn level_crossings(
    knots: &KnotVector,
    control: &[Weighted<Point>],
    carry: &Carry,
    level: f64,
    tol: Tolerances,
) -> Option<Vec<f64>> {
    let signed: Vec<Weighted<f64>> = control
        .iter()
        .map(|w| Weighted {
            scaled: w.weight * (carry.height(w.point()) - level),
            weight: w.weight,
        })
        .collect();
    // Each coefficient is a difference of two heights, rounded to the
    // larger's last bits, and so is every blend of them.
    let noise = control
        .iter()
        .map(|w| w.weight * (carry.height(w.point()).abs() + level.abs()))
        .fold(0.0, f64::max)
        * f64::EPSILON
        * 64.0;
    let pieces = bspline::to_bezier_segments(knots, &signed, tol).ok()?;
    let mut found = Vec::new();
    for ((a, b), piece) in pieces {
        let values: Vec<f64> = piece.iter().map(|w| w.scaled).collect();
        bezier_roots(&values, (a, b), noise, tol, &mut found)?;
    }
    found.sort_by(f64::total_cmp);
    found.dedup_by(|x, y| (*x - *y).abs() <= tol.parametric());
    Some(found)
}

/// The roots over `(a, b)` of the polynomial whose Bernstein coefficients
/// over it are `values`, pushed onto `found` (see [`level_crossings`]).
fn bezier_roots(
    values: &[f64],
    (a, b): (f64, f64),
    noise: f64,
    tol: Tolerances,
    found: &mut Vec<f64>,
) -> Option<()> {
    if values.iter().all(|&c| c > 0.0) || values.iter().all(|&c| c < 0.0) {
        return Some(());
    }
    if values.iter().all(|&c| c.abs() <= noise) {
        return None;
    }
    let (first, last) = (values[0], values[values.len() - 1]);
    let changes = values
        .windows(2)
        .filter(|pair| pair[0] * pair[1] < 0.0)
        .count();
    if changes == 1 && first * last < 0.0 && values.iter().all(|&c| c != 0.0) {
        let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
        loop {
            let mid = f64::midpoint(lo, hi);
            if mid <= lo || mid >= hi {
                break;
            }
            if de_casteljau(values, mid).0[values.len() - 1] * first > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        found.push((b - a).mul_add(f64::midpoint(lo, hi), a));
        return Some(());
    }
    if b - a <= tol.parametric() {
        return None;
    }
    let (left, right) = de_casteljau(values, 0.5);
    let mid = f64::midpoint(a, b);
    bezier_roots(&left, (a, mid), noise, tol, found)?;
    bezier_roots(&right, (mid, b), noise, tol, found)
}

/// Bernstein coefficients split at `s` by de Casteljau's construction:
/// those of the stretch before it and after it, the first's last value
/// the polynomial's at `s`.
fn de_casteljau(values: &[f64], s: f64) -> (Vec<f64>, Vec<f64>) {
    let mut row = values.to_vec();
    let mut left = Vec::with_capacity(values.len());
    let mut right = Vec::with_capacity(values.len());
    left.push(row[0]);
    right.push(row[row.len() - 1]);
    for len in (1..row.len()).rev() {
        for k in 0..len {
            row[k] = (row[k + 1] - row[k]).mul_add(s, row[k]);
        }
        left.push(row[0]);
        right.push(row[len - 1]);
    }
    right.reverse();
    (left, right)
}

/// The sections of an extrusion by a plane holding its direction: the
/// rulings through the profile's crossings with the plane, each over the
/// whole sweep and charted at its `u`. A crossing at an open end of the
/// profile is the extrusion's border lying in the plane, and is kept where
/// the extrusion leaves the plane there at an angle.
fn rulings(
    profile: &Profile,
    window: (f64, f64),
    whole_turn: bool,
    (n, o): (Vector, Point),
    extrusion: &ogeom_geom::ExtrusionSurface,
    tol: Tolerances,
) -> Option<Vec<(Curve, PlanarCurve)>> {
    let d = extrusion.direction();
    let (_, sweep) = extrusion.domain();
    let k = n.dot(o.to_vector());
    let mut roots = profile.crossings(n, k, window, tol)?;
    // An end of the profile on the plane to within the confusion distance
    // is a crossing there, whichever side rounding left it.
    for end in [window.0, window.1] {
        let gap = n.dot(profile.point(end).to_vector()) - k;
        if gap.abs() <= tol.confusion()
            && !roots.iter().any(|r| (r - end).abs() <= tol.parametric())
        {
            roots.push(end);
        }
    }
    roots.sort_by(f64::total_cmp);
    if whole_turn
        && roots.len() > 1
        && (roots[0] - window.0).abs() <= tol.parametric()
        && (roots[roots.len() - 1] - window.1).abs() <= tol.parametric()
    {
        roots.pop();
    }
    let mut out = Vec::with_capacity(roots.len());
    for u in roots {
        let tangent = profile.tangent(u);
        let across = tangent - d.vector() * tangent.dot(d.vector());
        let size = across.magnitude();
        if size <= tol.confusion() || n.dot(tangent).abs() < TANGENTIAL_SINE * size {
            return None;
        }
        let axis = Axis {
            location: profile.point(u),
            direction: d,
        };
        let line: Curve = LineCurve::over(axis, sweep.0, sweep.1).ok()?.into();
        let chart: PlanarCurve = ogeom_geom::Line2d::over(
            Axis2::new(Point2::new(u, 0.0), Direction2::Y),
            sweep.0,
            sweep.1,
        )
        .ok()?
        .into();
        out.push((line, chart));
    }
    Some(out)
}

/// The affine map carrying a point along `d` onto the plane through `o`
/// with normal `n`, where `along = n·d` is not zero.
struct Carry {
    n: Vector,
    o: Point,
    d: Vector,
    along: f64,
}

impl Carry {
    /// How far along `d` the point `p` is carried to reach the plane.
    fn height(&self, p: Point) -> f64 {
        self.n.dot(self.o - p) / self.along
    }

    /// The point `p` carried onto the plane.
    fn point(&self, p: Point) -> Point {
        p + self.d * self.height(p)
    }

    /// A vector carried as the difference of two points is.
    fn vector(&self, w: Vector) -> Vector {
        w - self.d * (self.n.dot(w) / self.along)
    }
}

/// The sections of an extrusion by a plane across its direction: the
/// profile carried onto the plane, over each stretch of the profile whose
/// height stays within the sweep, each with its chart trace on the
/// extrusion.
fn carried(
    profile: &Profile,
    window: (f64, f64),
    whole_turn: bool,
    carry: &Carry,
    sweep: (f64, f64),
    tol: Tolerances,
) -> Option<Vec<(Curve, PlanarCurve)>> {
    let level = |v: f64| carry.n.dot(carry.o.to_vector()) - v * carry.along;
    let mut cuts = vec![window.0, window.1];
    cuts.extend(profile.crossings(carry.n, level(sweep.0), window, tol)?);
    cuts.extend(profile.crossings(carry.n, level(sweep.1), window, tol)?);
    cuts.sort_by(f64::total_cmp);
    cuts.dedup_by(|x, y| (*x - *y).abs() <= tol.parametric());
    let mut stretches: Vec<(f64, f64)> = Vec::new();
    for pair in cuts.windows(2) {
        let (lo, hi) = (pair[0], pair[1]);
        let v = carry.height(profile.point(f64::midpoint(lo, hi)));
        if v < sweep.0 || v > sweep.1 {
            continue;
        }
        match stretches.last_mut() {
            Some(last) if (last.1 - lo).abs() <= tol.parametric() => last.1 = hi,
            _ => stretches.push((lo, hi)),
        }
    }
    // Round a whole turn, a stretch ending on its end runs on into one
    // starting on its start.
    let period = window.1 - window.0;
    if whole_turn && stretches.len() > 1 {
        let (first, last) = (stretches[0], stretches[stretches.len() - 1]);
        if (first.0 - window.0).abs() <= tol.parametric()
            && (last.1 - window.1).abs() <= tol.parametric()
        {
            stretches.pop();
            stretches[0] = (last.0, first.1 + period);
        }
    }
    let image = Image::of(profile, carry, tol)?;
    let mut out = Vec::with_capacity(stretches.len());
    for (lo, hi) in stretches {
        if hi - lo <= tol.parametric() {
            continue;
        }
        let whole = whole_turn && (hi - lo - period).abs() <= tol.parametric();
        out.push(image.over((lo, hi), whole, profile, carry, tol)?);
    }
    Some(out)
}

/// The profile carried onto the plane, in its own spelling.
enum Image {
    /// A line, at `scale` times the profile's parameter.
    Line { axis: Axis, scale: f64 },
    /// An ellipse or a circle, at the profile's parameter less `phase`.
    Conic { curve: Curve, phase: f64 },
    /// A parabola, at `scale` times the profile's parameter less `vertex`.
    Parabola {
        parabola: ogeom_math::Parabola,
        vertex: f64,
        scale: f64,
    },
}

impl Image {
    fn of(profile: &Profile, carry: &Carry, tol: Tolerances) -> Option<Self> {
        Some(match profile {
            Profile::Line { p, t } => {
                let run = carry.vector(*t);
                let scale = run.magnitude();
                Self::Line {
                    axis: Axis {
                        location: carry.point(*p),
                        direction: Direction::new(run, tol).ok()?,
                    },
                    scale,
                }
            }
            Profile::Conic { c, x, y } => {
                // Conjugate half-diameters `p`, `q`; turned by `phase` they
                // are the principal ones, `p cos φ + q sin φ` the major.
                let (p, q) = (carry.vector(*x), carry.vector(*y));
                let phase = 0.5 * (2.0 * p.dot(q)).atan2(p.dot(p) - q.dot(q));
                let (cos, sin) = (phase.cos(), phase.sin());
                let major = p * cos + q * sin;
                let minor = q * cos - p * sin;
                let (ra, rb) = (major.magnitude(), minor.magnitude());
                let xa = Direction::new(major, tol).ok()?;
                let ya = Direction::new(minor - xa.vector() * minor.dot(xa.vector()), tol).ok()?;
                let za = Direction::new(xa.vector().cross(ya.vector()), tol).ok()?;
                let frame = Frame::from_axes(carry.point(*c), xa, ya, za, tol).ok()?;
                let curve: Curve = if ra - rb <= tol.confusion() {
                    ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(frame, ra, tol).ok()?)
                        .into()
                } else {
                    ogeom_geom::EllipseCurve::new(
                        ogeom_math::Ellipse::new(frame, ra, rb, tol).ok()?,
                    )
                    .into()
                };
                Self::Conic { curve, phase }
            }
            Profile::Parabola { apex, x, y } => {
                // `V + s² p + s q⊥` from the vertex, where the tangent is
                // square to the axis `p`.
                let (p, q) = (carry.vector(*x), carry.vector(*y));
                let vertex = -p.dot(q) / (2.0 * p.dot(p));
                let across = q - p * (p.dot(q) / p.dot(p));
                let scale = across.magnitude();
                let xa = Direction::new(p, tol).ok()?;
                let ya = Direction::new(across, tol).ok()?;
                let za = Direction::new(xa.vector().cross(ya.vector()), tol).ok()?;
                let origin = carry.point(*apex) + p * (vertex * vertex) + q * vertex;
                let frame = Frame::from_axes(origin, xa, ya, za, tol).ok()?;
                let focal = scale * scale / (4.0 * p.magnitude());
                Self::Parabola {
                    parabola: ogeom_math::Parabola::new(frame, focal, tol).ok()?,
                    vertex,
                    scale,
                }
            }
        })
    }

    /// The image over the profile's stretch `(lo, hi)`, the whole curve
    /// where `whole`, and its chart trace on the extrusion: `(u, v)` at
    /// the image's own parameter.
    fn over(
        &self,
        (lo, hi): (f64, f64),
        whole: bool,
        profile: &Profile,
        carry: &Carry,
        tol: Tolerances,
    ) -> Option<(Curve, PlanarCurve)> {
        let height = |u: f64| carry.height(profile.point(u));
        match self {
            Self::Line { axis, scale } => {
                let (t0, t1) = (lo * scale, hi * scale);
                let curve: Curve = LineCurve::over(*axis, t0, t1).ok()?.into();
                let chart = BSpline2d::new(
                    KnotVector::new(vec![t0, t0, t1, t1], 1).ok()?,
                    vec![Point2::new(lo, height(lo)), Point2::new(hi, height(hi))],
                    tol,
                )
                .ok()?;
                Some((curve, chart.into()))
            }
            Self::Conic { curve, phase } => {
                // `v = α + β cos u + γ sin u` at `u = t + phase`.
                let Profile::Conic { c, x, y } = profile else {
                    return None;
                };
                let alpha = carry.height(*c);
                let (beta, gamma) = (
                    -carry.n.dot(*x) / carry.along,
                    -carry.n.dot(*y) / carry.along,
                );
                let (cos, sin) = (phase.cos(), phase.sin());
                let range = if whole {
                    curve.domain()
                } else {
                    (lo - phase, hi - phase)
                };
                let chart = ogeom_geom::Trig2d::new(
                    Point2::new(*phase, alpha),
                    Vector2::new(1.0, 0.0),
                    Vector2::new(0.0, beta.mul_add(cos, gamma * sin)),
                    Vector2::new(0.0, gamma.mul_add(cos, -beta * sin)),
                    range,
                )
                .ok()?;
                let piece = if whole {
                    curve.clone()
                } else {
                    ogeom_geom::TrimmedCurve::new(curve.clone(), range.0, range.1, tol)
                        .ok()?
                        .into()
                };
                Some((piece, PlanarCurve::Trig(chart)))
            }
            Self::Parabola {
                parabola,
                vertex,
                scale,
            } => {
                let (t0, t1) = ((lo - vertex) * scale, (hi - vertex) * scale);
                let curve: Curve = ogeom_geom::ParabolaCurve::over(*parabola, t0, t1)
                    .ok()?
                    .into();
                // Linear in `u`, quadratic in `v`: one quadratic Bézier, its
                // middle control point where the end tangents meet.
                let rise = -carry.n.dot(profile.tangent(lo)) / carry.along;
                let start = Point2::new(lo, height(lo));
                let middle = Point2::new(
                    f64::midpoint(lo, hi),
                    rise.mul_add((hi - lo) / 2.0, start.y),
                );
                let end = Point2::new(hi, height(hi));
                let chart = BSpline2d::new(
                    KnotVector::new(vec![t0, t0, t0, t1, t1, t1], 2).ok()?,
                    vec![start, middle, end],
                    tol,
                )
                .ok()?;
                Some((curve, chart.into()))
            }
        }
    }
}

/// The piece and its chart traces cut to the plane's window: a straight
/// one to the stretch the window holds, any other kept whole where it
/// reaches the window and dropped where it does not, as the closed forms'
/// sections are. `Some(None)` for a piece wholly off the window, `None`
/// where a trace has no closed form.
#[allow(clippy::type_complexity, reason = "a piece and its two traces")]
fn within_plane_window(
    curve: Curve,
    on_extrusion: PlanarCurve,
    plane: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<Option<(Curve, PlanarCurve, PlanarCurve)>> {
    let range = curve.domain();
    let on_plane = exact_pcurve_over(&curve, range, plane, tol)?;
    let Curve::Line(line) = &curve else {
        return Some(touches_box(&on_plane, plane, tol).then_some((curve, on_extrusion, on_plane)));
    };
    let Some((lo, hi)) = inside_box(&on_plane, plane) else {
        return Some(None);
    };
    let (lo, hi) = (lo.max(range.0), hi.min(range.1));
    if hi - lo <= tol.parametric() {
        return Some(None);
    }
    if lo <= range.0 && hi >= range.1 {
        return Some(Some((curve, on_extrusion, on_plane)));
    }
    let cut: Curve = LineCurve::over(line.axis(), lo, hi).ok()?.into();
    let on_extrusion: PlanarCurve = match on_extrusion {
        PlanarCurve::Line(l) => ogeom_geom::Line2d::over(l.axis(), lo, hi).ok()?.into(),
        other => {
            let (a, b) = (other.point_at(lo, tol).ok()?, other.point_at(hi, tol).ok()?);
            BSpline2d::new(
                KnotVector::new(vec![lo, lo, hi, hi], 1).ok()?,
                vec![a, b],
                tol,
            )
            .ok()?
            .into()
        }
    };
    let on_plane = exact_pcurve_over(&cut, (lo, hi), plane, tol)?;
    Some(Some((cut, on_extrusion, on_plane)))
}

/// Whether a section lands where it says: on the plane, and on each
/// surface at the chart point its trace gives, within the confusion
/// distance at stations along it. `None` where it does not.
fn lands_on_both(
    curve: &Curve,
    on_plane: &PlanarCurve,
    on_extrusion: &PlanarCurve,
    plane: &SurfaceGeometry,
    extrusion: &ogeom_geom::ExtrusionSurface,
    tol: Tolerances,
) -> Option<()> {
    const STATIONS: u32 = 16;
    let (lo, hi) = curve.domain();
    for i in 0..=STATIONS {
        let t = (hi - lo).mul_add(f64::from(i) / f64::from(STATIONS), lo);
        let point = curve.point_at(t, tol).ok()?;
        let flat = on_plane.point_at(t, tol).ok()?;
        let swept = on_extrusion.point_at(t, tol).ok()?;
        let on_plane_at = plane.point_at(flat.x, flat.y, tol).ok()?;
        let on_extrusion_at = extrusion.point_at(swept.x, swept.y, tol).ok()?;
        if point.distance(on_plane_at) > tol.confusion()
            || point.distance(on_extrusion_at) > tol.confusion()
        {
            return None;
        }
    }
    Some(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use ogeom_geom::{
        BSplineCurve, CircleCurve, Curve, Curve2d as _, Curve3d as _, EllipseCurve,
        ExtrusionSurface, LineCurve, ParabolaCurve, PlanarCurve, PlaneSurface, Surface as _,
        SurfaceGeometry, TrimmedCurve,
    };
    use ogeom_math::bspline::Weighted;
    use ogeom_math::{
        Axis, Circle, Direction, Ellipse, Frame, KnotVector, Parabola, Plane, Point, Vector,
    };
    use proptest::prelude::*;

    use crate::{IntersectOptions, SectionCurve, SurfaceIntersection, intersect_surfaces};

    const T: ogeom_core::Tolerances = ogeom_core::Tolerances::millimetres();

    fn direction(x: f64, y: f64, z: f64) -> Direction {
        Direction::new(Vector::new(x, y, z), T).unwrap()
    }

    /// The sections of `profile` swept along `along` over `sweep` by the
    /// plane through `at` square to `normal`, asked with the plane first;
    /// asked the other way round they must be the same curves, their
    /// traces swapped.
    fn cut(
        profile: Curve,
        along: Direction,
        sweep: (f64, f64),
        at: Point,
        normal: Direction,
    ) -> (Vec<SectionCurve>, SurfaceGeometry, SurfaceGeometry) {
        let extrusion: SurfaceGeometry = ExtrusionSurface::over(profile, along, sweep)
            .unwrap()
            .into();
        let plane: SurfaceGeometry = PlaneSurface::new(Plane::through(at, normal)).into();
        let ask = |a: &SurfaceGeometry, b: &SurfaceGeometry| -> Vec<SectionCurve> {
            match intersect_surfaces(a, b, IntersectOptions::default(), T).unwrap() {
                SurfaceIntersection::Along(found) => found,
                SurfaceIntersection::Apart => Vec::new(),
                other => panic!("not a section: {other:?}"),
            }
        };
        let first = ask(&plane, &extrusion);
        let second = ask(&extrusion, &plane);
        assert_eq!(first.len(), second.len());
        for (x, y) in first.iter().zip(&second) {
            assert_eq!(x.curve, y.curve);
            assert_eq!((&x.on_a, &x.on_b), (&y.on_b, &y.on_a));
        }
        (first, plane, extrusion)
    }

    /// The farthest a section stands off either surface, measured at its
    /// own points and at each trace's chart point, at stations along it.
    fn stray(section: &SectionCurve, plane: &SurfaceGeometry, extrusion: &SurfaceGeometry) -> f64 {
        let SurfaceGeometry::Plane(flat) = plane else {
            unreachable!()
        };
        let (lo, hi) = section.curve.domain();
        let (Some(on_plane), Some(on_extrusion)) = (&section.on_a, &section.on_b) else {
            panic!("a section without its traces");
        };
        (0..=64)
            .map(|i| {
                let t = (hi - lo).mul_add(f64::from(i) / 64.0, lo);
                let p = section.curve.point_at(t, T).unwrap();
                let a = on_plane.point_at(t, T).unwrap();
                let b = on_extrusion.point_at(t, T).unwrap();
                flat.plane()
                    .signed_distance_to(p)
                    .abs()
                    .max(plane.point_at(a.x, a.y, T).unwrap().distance(p))
                    .max(extrusion.point_at(b.x, b.y, T).unwrap().distance(p))
            })
            .fold(0.0, f64::max)
    }

    fn exact_on_both(found: &[SectionCurve], plane: &SurfaceGeometry, extrusion: &SurfaceGeometry) {
        for section in found {
            assert!(section.exact && section.tolerance == 0.0 && !section.tangential);
            let off = stray(section, plane, extrusion);
            assert!(off <= 1e-12, "a section stands {off} off");
        }
    }

    /// A line swept along a slant and cut level is the line where the
    /// sweep reaches the plane's height, over the profile's own stretch.
    #[test]
    fn a_swept_line_cut_level_is_the_line_at_that_height() {
        let profile: Curve = LineCurve::over(
            Axis {
                location: Point::ORIGIN,
                direction: Direction::X,
            },
            -2.0,
            3.0,
        )
        .unwrap()
        .into();
        let (found, plane, extrusion) = cut(
            profile,
            direction(0.0, 1.0, 1.0),
            (-5.0, 5.0),
            Point::new(0.0, 0.0, 1.0),
            Direction::Z,
        );
        exact_on_both(&found, &plane, &extrusion);
        let [section] = found.as_slice() else {
            panic!("{} sections", found.len());
        };
        assert!(matches!(section.curve, Curve::Line(_)));
        let (lo, hi) = section.curve.domain();
        let ends = [lo, hi].map(|t| section.curve.point_at(t, T).unwrap());
        assert!(ends[0].distance(Point::new(-2.0, 1.0, 1.0)) < 1e-12);
        assert!(ends[1].distance(Point::new(3.0, 1.0, 1.0)) < 1e-12);
    }

    /// A drum's circle swept along its axis and cut at a tilt `θ` is the
    /// closed ellipse of radii `R / cos θ` and `R` about the axis.
    #[test]
    fn a_swept_circle_cut_at_a_tilt_is_an_ellipse_of_the_tilt() {
        let (radius, tilt) = (2.0, 0.6_f64);
        let centre = Point::new(1.0, 0.0, 0.0);
        let circle = Circle::new(
            Frame::new(centre, Direction::Z, Direction::X, T).unwrap(),
            radius,
            T,
        )
        .unwrap();
        let (found, plane, extrusion) = cut(
            CircleCurve::new(circle).into(),
            Direction::Z,
            (-10.0, 10.0),
            centre,
            direction(tilt.sin(), 0.0, tilt.cos()),
        );
        exact_on_both(&found, &plane, &extrusion);
        let [section] = found.as_slice() else {
            panic!("{} sections", found.len());
        };
        let Curve::Ellipse(e) = &section.curve else {
            panic!("not an ellipse: {:?}", section.curve);
        };
        let e = e.ellipse();
        assert!(section.closed);
        assert!((e.major_radius() - radius / tilt.cos()).abs() < 1e-12);
        assert!((e.minor_radius() - radius).abs() < 1e-12);
        assert!(e.centre().distance(centre) < 1e-12);
    }

    /// An ellipse swept along a slant and cut by a plane parallel to its
    /// own is the same ellipse, moved along the slant to the plane.
    #[test]
    fn a_swept_ellipse_cut_parallel_to_it_is_the_ellipse_moved() {
        let ellipse = Ellipse::new(Frame::WORLD, 3.0, 1.0, T).unwrap();
        let (found, plane, extrusion) = cut(
            EllipseCurve::new(ellipse).into(),
            direction(1.0, 0.0, 1.0),
            (-1.0, 5.0),
            Point::new(0.0, 0.0, 2.0),
            Direction::Z,
        );
        exact_on_both(&found, &plane, &extrusion);
        let [section] = found.as_slice() else {
            panic!("{} sections", found.len());
        };
        let Curve::Ellipse(e) = &section.curve else {
            panic!("not an ellipse: {:?}", section.curve);
        };
        let e = e.ellipse();
        assert!((e.major_radius() - 3.0).abs() < 1e-12);
        assert!((e.minor_radius() - 1.0).abs() < 1e-12);
        assert!(e.centre().distance(Point::new(2.0, 0.0, 2.0)) < 1e-12);
    }

    /// A parabola of focal length `f` swept square to its plane and cut at
    /// a tilt `θ` about its axis's normal is the parabola of focal length
    /// `f cos θ` from the same apex: its axis stretched by `1 / cos θ`.
    #[test]
    fn a_swept_parabola_cut_at_a_tilt_has_its_focal_length_shortened() {
        let (focal, tilt) = (0.5, 0.4_f64);
        let parabola = Parabola::new(Frame::WORLD, focal, T).unwrap();
        let (found, plane, extrusion) = cut(
            ParabolaCurve::over(parabola, -2.0, 2.0).unwrap().into(),
            Direction::Z,
            (-10.0, 10.0),
            Point::ORIGIN,
            direction(-tilt.sin(), 0.0, tilt.cos()),
        );
        exact_on_both(&found, &plane, &extrusion);
        let [section] = found.as_slice() else {
            panic!("{} sections", found.len());
        };
        let Curve::Parabola(p) = &section.curve else {
            panic!("not a parabola: {:?}", section.curve);
        };
        assert!((p.parabola().focal() - focal * tilt.cos()).abs() < 1e-12);
        assert!(p.parabola().apex().distance(Point::ORIGIN) < 1e-12);
    }

    /// A plane holding the sweep's direction cuts it in the rulings
    /// through the profile's crossings with it, over the whole sweep, each
    /// charted at the angle it stands at.
    #[test]
    fn a_plane_along_the_sweep_cuts_it_in_rulings() {
        let circle = Circle::new(Frame::WORLD, 2.0, T).unwrap();
        let (found, plane, extrusion) = cut(
            CircleCurve::new(circle).into(),
            Direction::Z,
            (0.0, 3.0),
            Point::new(1.0, 0.0, 0.0),
            Direction::X,
        );
        exact_on_both(&found, &plane, &extrusion);
        assert_eq!(found.len(), 2);
        let third = core::f64::consts::FRAC_PI_3;
        for (section, (y, u)) in found
            .iter()
            .zip([(3.0_f64.sqrt(), third), (-(3.0_f64.sqrt()), 5.0 * third)])
        {
            let Curve::Line(line) = &section.curve else {
                panic!("not a ruling: {:?}", section.curve);
            };
            assert!(line.axis().direction.is_parallel(Direction::Z, T));
            let (lo, hi) = section.curve.domain();
            let ends = [lo, hi].map(|t| section.curve.point_at(t, T).unwrap());
            assert!(ends[0].distance(Point::new(1.0, y, 0.0)) < 1e-12);
            assert!(ends[1].distance(Point::new(1.0, y, 3.0)) < 1e-12);
            let Some(PlanarCurve::Line(chart)) = &section.on_b else {
                panic!("a ruling charted as {:?}", section.on_b);
            };
            assert!((chart.axis().location.x - u).abs() < 1e-12);
        }
    }

    /// A circle swept over `[0, 1]` and cut by `z = 1/2 + x`: the height
    /// `1/2 + cos u` leaves the sweep, and what stays is the two arcs over
    /// `cos u` in `[−1/2, 1/2]`, each ending on the sweep's two rims.
    #[test]
    fn a_section_leaving_the_sweep_is_cut_to_it() {
        let circle = Circle::new(Frame::WORLD, 1.0, T).unwrap();
        let (found, plane, extrusion) = cut(
            CircleCurve::new(circle).into(),
            Direction::Z,
            (0.0, 1.0),
            Point::new(0.0, 0.0, 0.5),
            direction(-1.0, 0.0, 1.0),
        );
        exact_on_both(&found, &plane, &extrusion);
        assert_eq!(found.len(), 2);
        for section in &found {
            assert!(!section.closed);
            let (lo, hi) = section.curve.domain();
            let mut heights = [lo, hi].map(|t| section.curve.point_at(t, T).unwrap().z);
            heights.sort_by(f64::total_cmp);
            assert!(heights[0].abs() < 1e-12 && (heights[1] - 1.0).abs() < 1e-12);
        }
    }

    /// Cut by `z = 1/2 + 4y/5`, the circle swept over `[0, 1]` keeps
    /// `sin u` in `[−5/8, 5/8]`: two arcs, one through `(1, 0, 1/2)` across
    /// the profile's own start, stated as one section rather than two
    /// halves.
    #[test]
    fn an_arc_across_the_profile_s_start_is_one_section() {
        let circle = Circle::new(Frame::WORLD, 1.0, T).unwrap();
        let (found, plane, extrusion) = cut(
            CircleCurve::new(circle).into(),
            Direction::Z,
            (0.0, 1.0),
            Point::new(0.0, 0.0, 0.5),
            direction(0.0, -0.8, 1.0),
        );
        exact_on_both(&found, &plane, &extrusion);
        assert_eq!(found.len(), 2);
        let across = found
            .iter()
            .filter(|section| {
                let (lo, hi) = section.curve.domain();
                let at = |t: f64| section.curve.point_at(t, T).unwrap();
                let mid = at(f64::midpoint(lo, hi));
                let ends = [at(lo).z, at(hi).z];
                assert!(ends.iter().any(|z| z.abs() < 1e-12));
                assert!(ends.iter().any(|z| (z - 1.0).abs() < 1e-12));
                mid.distance(Point::new(1.0, 0.0, 0.5)) < 1e-12
            })
            .count();
        assert_eq!(across, 1);
    }

    /// A circle's arc spelled as a rational quadratic B-spline, swept along
    /// its axis and cut at a tilt `θ`: every point of the section lies on
    /// the drum and in the plane, on the ellipse of radii `R / cos θ` and
    /// `R`, and the section at the profile's parameter `t` is the arc's
    /// point there lifted to the plane.
    #[test]
    fn a_swept_arc_spline_cut_at_a_tilt_is_the_arc_lifted_to_the_plane() {
        let (radius, tilt) = (2.0, 0.6_f64);
        let centre = Point::new(1.0, 0.0, 0.0);
        let circle = Circle::new(
            Frame::new(centre, Direction::Z, Direction::X, T).unwrap(),
            radius,
            T,
        )
        .unwrap();
        let arc: Curve = Curve::from(CircleCurve::new(circle))
            .to_bspline_over((0.3, 2.5), T)
            .unwrap()
            .into();
        let (found, plane, extrusion) = cut(
            arc.clone(),
            Direction::Z,
            (-10.0, 10.0),
            centre,
            direction(tilt.sin(), 0.0, tilt.cos()),
        );
        exact_on_both(&found, &plane, &extrusion);
        let [section] = found.as_slice() else {
            panic!("{} sections", found.len());
        };
        assert!(matches!(section.curve, Curve::BSpline(_)));
        assert_eq!(section.curve.domain(), arc.domain());
        let (lo, hi) = arc.domain();
        for i in 0..=32 {
            let t = (hi - lo).mul_add(f64::from(i) / 32.0, lo);
            let p = section.curve.point_at(t, T).unwrap();
            let q = arc.point_at(t, T).unwrap();
            let lifted = Point::new(q.x, q.y, -tilt.tan() * (q.x - 1.0));
            assert!(p.distance(lifted) < 1e-12);
            let (x, y) = ((p.x - 1.0).mul_add(tilt.cos(), -p.z * tilt.sin()), p.y);
            let on = (x * tilt.cos() / radius).powi(2) + (y / radius).powi(2);
            assert!((on - 1.0).abs() < 1e-12);
        }
    }

    /// A polynomial cubic winding in space, swept along a slant `d` and cut
    /// by `z = 1`: the section at `t` is `C(t) + d (1 − C_z(t)) / d_z`, and
    /// its trace on the extrusion is `(t, (1 − C_z(t)) / d_z)`. Swept short
    /// of the plane it is apart, and swept to straddle it the pair is left
    /// to the marcher.
    #[test]
    fn a_swept_cubic_cut_level_is_the_cubic_carried_along_the_sweep() {
        let control = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 2.0, 0.5),
            Point::new(2.0, -1.0, -0.4),
            Point::new(3.0, 1.0, 0.3),
            Point::new(4.0, 0.0, 0.0),
        ];
        let cubic: Curve =
            BSplineCurve::new(KnotVector::clamped_uniform(3, 5).unwrap(), control, T)
                .unwrap()
                .into();
        let along = direction(0.0, 1.0, 1.0);
        let (found, plane, extrusion) = cut(
            cubic.clone(),
            along,
            (-5.0, 5.0),
            Point::new(0.0, 0.0, 1.0),
            Direction::Z,
        );
        exact_on_both(&found, &plane, &extrusion);
        let [section] = found.as_slice() else {
            panic!("{} sections", found.len());
        };
        let Some(trace) = &section.on_b else {
            panic!("a section without its trace");
        };
        let d = along.vector();
        let (lo, hi) = cubic.domain();
        for i in 0..=32 {
            let t = (hi - lo).mul_add(f64::from(i) / 32.0, lo);
            let q = cubic.point_at(t, T).unwrap();
            let height = (1.0 - q.z) / d.z;
            assert!(
                section
                    .curve
                    .point_at(t, T)
                    .unwrap()
                    .distance(q + d * height)
                    < 1e-12
            );
            let uv = trace.point_at(t, T).unwrap();
            assert!((uv.x - t).abs() < 1e-12 && (uv.y - height).abs() < 1e-12);
        }
        let plane: SurfaceGeometry =
            PlaneSurface::new(Plane::through(Point::new(0.0, 0.0, 1.0), Direction::Z)).into();
        let short: SurfaceGeometry = ExtrusionSurface::over(cubic.clone(), along, (-5.0, 0.5))
            .unwrap()
            .into();
        assert_eq!(
            super::plane_extrusion(&plane, &short, T).map(|f| f.len()),
            Some(0)
        );
        let straddling: SurfaceGeometry = ExtrusionSurface::over(cubic, along, (-5.0, 1.0))
            .unwrap()
            .into();
        assert!(super::plane_extrusion(&plane, &straddling, T).is_none());
    }

    /// The cubic above swept over `[−5, 1.6]` and cut by `z = 1`, its
    /// heights straddling the sweep's top: split at its rim, each section
    /// lies on both, ends at the profile's own end or on the rim,
    /// and every point of the profile whose height falls inside the sweep
    /// is on one.
    #[test]
    fn a_swept_cubic_straddling_its_rim_is_split_there() {
        let control = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 2.0, 0.5),
            Point::new(2.0, -1.0, -0.4),
            Point::new(3.0, 1.0, 0.3),
            Point::new(4.0, 0.0, 0.0),
        ];
        let cubic: Curve =
            BSplineCurve::new(KnotVector::clamped_uniform(3, 5).unwrap(), control, T)
                .unwrap()
                .into();
        let along = direction(0.0, 1.0, 1.0);
        let plane: SurfaceGeometry =
            PlaneSurface::new(Plane::through(Point::new(0.0, 0.0, 1.0), Direction::Z)).into();
        let top = 1.6;
        let extrusion: SurfaceGeometry = ExtrusionSurface::over(cubic.clone(), along, (-5.0, top))
            .unwrap()
            .into();
        assert!(super::crosses_extrusion_rim(&plane, &extrusion, T));
        let found = super::plane_extrusion_at_rims(&plane, &extrusion, T).unwrap();
        assert!(!found.is_empty());
        exact_on_both(&found, &plane, &extrusion);
        let height = |t: f64| (1.0 - cubic.point_at(t, T).unwrap().z) / along.vector().z;
        let (lo, hi) = cubic.domain();
        let mut spans = Vec::new();
        for section in &found {
            let trace = section.on_b.as_ref().unwrap();
            let (a, b) = trace.domain();
            for end in [a, b] {
                let uv = trace.point_at(end, T).unwrap();
                assert!(
                    (uv.x - lo).abs() < 1e-12
                        || (uv.x - hi).abs() < 1e-12
                        || (uv.y - top).abs() < 1e-12,
                    "an end at {uv:?} neither on the profile's end nor on the rim"
                );
            }
            spans.push((
                trace.point_at(a, T).unwrap().x,
                trace.point_at(b, T).unwrap().x,
            ));
        }
        for i in 0..=400 {
            let t = (hi - lo).mul_add(f64::from(i) / 400.0, lo);
            if height(t) >= top - 1e-9 {
                continue;
            }
            assert!(
                spans
                    .iter()
                    .any(|&(a, b)| t >= a.min(b) - 1e-12 && t <= a.max(b) + 1e-12),
                "t {t} at height {} is on no section",
                height(t)
            );
        }
    }

    /// A half circle of radius five spelled as a rational quadratic, swept
    /// ten up, and cut by a plane leaning a hundred-thousandth of a radian
    /// that crosses the swept top a ten-thousandth inside the arc's crown,
    /// at the chord `d = r − 10⁻⁴` from its centre. Past the chord the plane
    /// stands above the top by at most `10⁻⁹`, within the confusion distance
    /// over the whole crown, and meets the top rim all but tangentially. It
    /// is split there, in two: each section lies on both, runs from the
    /// arc's end to the rim, and ends on the rim where the plane crosses
    /// it, at `x = 5 ± √(r² − d²)`. The section meets the rim within a
    /// rounding `δ` of height over a stretch `δ / (θ |dy/dx|)` long, at the
    /// lean `θ` and the arc's slope there `√(r² − d²) / d`, a few tenths of a
    /// micron for `δ` a few times `10⁻¹⁵`, and is held to a micron.
    #[test]
    fn a_rational_arc_grazing_its_rim_is_split_where_the_plane_crosses_it() {
        let (radius, lean, depth) = (5.0_f64, 1e-5_f64, 1e-4_f64);
        let centre = Point::new(5.0, 10.0, 0.0);
        let circle = Circle::new(
            Frame::new(centre, Direction::Z, Direction::X, T).unwrap(),
            radius,
            T,
        )
        .unwrap();
        let arc: Curve = Curve::from(CircleCurve::new(circle))
            .to_bspline_over((0.0, core::f64::consts::PI), T)
            .unwrap()
            .into();
        let chord = radius - depth;
        let plane: SurfaceGeometry = PlaneSurface::new(Plane::through(
            Point::new(5.0, 10.0 + chord, 10.0),
            direction(0.0, -lean.sin(), lean.cos()),
        ))
        .into();
        let extrusion: SurfaceGeometry = ExtrusionSurface::over(arc, Direction::Z, (0.0, 10.0))
            .unwrap()
            .into();
        let found = super::plane_extrusion_at_rims(&plane, &extrusion, T).unwrap();
        exact_on_both(&found, &plane, &extrusion);
        assert_eq!(found.len(), 2);
        let half = radius.mul_add(radius, -chord * chord).sqrt();
        let mut crossed = Vec::new();
        for section in &found {
            let (lo, hi) = section.curve.domain();
            let ends = [lo, hi].map(|t| section.curve.point_at(t, T).unwrap());
            let (base, rim): (Vec<Point>, Vec<Point>) =
                ends.iter().partition(|p| (p.y - 10.0).abs() < 1e-12);
            assert_eq!((base.len(), rim.len()), (1, 1), "ends {ends:?}");
            assert!((rim[0].z - 10.0).abs() < 1e-12);
            crossed.push(rim[0].x);
        }
        crossed.sort_by(f64::total_cmp);
        assert!((crossed[0] - (5.0 - half)).abs() < 1e-6, "{crossed:?}");
        assert!((crossed[1] - (5.0 + half)).abs() < 1e-6, "{crossed:?}");
        assert!(super::plane_extrusion(&plane, &extrusion, T).is_none());
    }

    /// The profiles the property test sweeps: a line, a circle, an
    /// ellipse, a parabola and a rational cubic, whole or trimmed.
    fn profile() -> impl Strategy<Value = Curve> {
        let point = || (-2.0..2.0f64, -2.0..2.0f64, -2.0..2.0f64);
        let frame = (point(), -3.0..3.0f64, -3.0..3.0f64).prop_map(|((x, y, z), a, b)| {
            let axis = direction(a.cos() * b.sin(), a.sin() * b.sin(), b.cos());
            let lean = Vector::new(1.0, 0.3, -0.2);
            let x_axis = Direction::new(lean - axis.vector() * lean.dot(axis.vector()), T).unwrap();
            Frame::new(Point::new(x, y, z), axis, x_axis, T).unwrap()
        });
        let kind = 0..5u32;
        let trim = proptest::option::of((0.05..0.45f64, 0.55..0.95f64));
        (kind, frame, 0.3..2.0f64, 0.2..1.0f64, trim).prop_map(
            |(kind, frame, size, ratio, trim)| {
                let curve: Curve = match kind {
                    0 => LineCurve::over(
                        Axis {
                            location: frame.origin(),
                            direction: frame.x(),
                        },
                        -size,
                        size * 2.0,
                    )
                    .unwrap()
                    .into(),
                    1 => CircleCurve::new(Circle::new(frame, size, T).unwrap()).into(),
                    2 => EllipseCurve::new(Ellipse::new(frame, size, size * ratio, T).unwrap())
                        .into(),
                    3 => ParabolaCurve::over(
                        Parabola::new(frame, size * ratio, T).unwrap(),
                        -size,
                        size,
                    )
                    .unwrap()
                    .into(),
                    // A rational cubic winding off its frame's plane.
                    _ => {
                        let control = (0..6u32)
                            .map(|i| {
                                let k = f64::from(i);
                                let p = frame.origin()
                                    + frame.x().vector() * (size * (k - 2.5))
                                    + frame.y().vector() * (size * ratio * (k * 1.7).sin())
                                    + frame.z().vector() * (size * 0.3 * (k * 0.9).cos());
                                Weighted::new(p, 1.0 + ratio * f64::from(i % 2), T).unwrap()
                            })
                            .collect();
                        BSplineCurve::rational(KnotVector::clamped_uniform(3, 6).unwrap(), control)
                            .unwrap()
                            .into()
                    }
                };
                match trim {
                    Some((a, b)) => {
                        let (lo, hi) = curve.domain();
                        let at = |f: f64| (hi - lo).mul_add(f, lo);
                        TrimmedCurve::new(curve, at(a), at(b), T).unwrap().into()
                    }
                    None => curve,
                }
            },
        )
    }

    proptest! {
        /// Whatever the profile, sweep and plane across it, every section
        /// lies on both surfaces to within 1e-12, at its own points and at
        /// its traces' chart points; and every point of the profile whose
        /// height to the plane falls inside the sweep is on a section.
        #[test]
        fn a_carried_section_lies_on_both_and_misses_nothing(
            profile in profile(),
            (a, b) in (-3.0..3.0f64, 0.0..3.0f64),
            (c, e) in (-3.0..3.0f64, 0.0..3.0f64),
            at in (-1.0..1.0f64, -1.0..1.0f64, -1.0..1.0f64),
            sweep in (-4.0..-0.5f64, 0.5..4.0f64),
        ) {
            let along = direction(a.cos() * b.sin(), a.sin() * b.sin(), b.cos());
            let normal = direction(c.cos() * e.sin(), c.sin() * e.sin(), e.cos());
            prop_assume!(normal.vector().dot(along.vector()).abs() >= 0.2);
            let at = Point::new(at.0, at.1, at.2);
            let extrusion = ExtrusionSurface::over(profile.clone(), along, sweep).unwrap();
            let plane_surface: SurfaceGeometry = PlaneSurface::new(Plane::through(at, normal)).into();
            let swept: SurfaceGeometry = extrusion.clone().into();
            let Some(found) = super::plane_extrusion(&plane_surface, &swept, T)
                .or_else(|| super::plane_extrusion_at_rims(&plane_surface, &swept, T))
            else {
                // Declined: a touch, or a section meeting the extrusion
                // all but tangentially. Neither is what this measures.
                return Ok(());
            };
            let extrusion_surface: SurfaceGeometry = extrusion.into();
            for section in &found {
                prop_assert!(stray(section, &plane_surface, &extrusion_surface) <= 1e-12);
            }
            // Every point whose height falls inside the sweep is charted
            // by a section: its `u` lies within one's trace.
            let spans: Vec<(f64, f64)> = found
                .iter()
                .map(|s| {
                    let trace = s.on_b.as_ref().unwrap();
                    let (lo, hi) = trace.domain();
                    let (p, q) = (trace.point_at(lo, T).unwrap().x, trace.point_at(hi, T).unwrap().x);
                    (p.min(q), p.max(q))
                })
                .collect();
            let (lo, hi) = profile.domain();
            let period = core::f64::consts::TAU;
            for i in 0..=200 {
                let u = (hi - lo).mul_add(f64::from(i) / 200.0, lo);
                let p = profile.point_at(u, T).unwrap();
                let v = normal.vector().dot(at - p) / normal.vector().dot(along.vector());
                if v <= sweep.0 + 1e-6 || v >= sweep.1 - 1e-6 {
                    continue;
                }
                let covered = spans.iter().any(|&(s0, s1)| {
                    [0.0, period, -period].iter().any(|k| u + k >= s0 - 1e-9 && u + k <= s1 + 1e-9)
                });
                prop_assert!(covered, "u {u} at height {v} is on no section");
            }
        }
    }

    /// A cylinder of radius 20 about `z`, 1 tall, and the extrusion along
    /// `z` of the degree-1 B-spline from `(-1, y)` to `(1, y)`.
    fn drum_and_side(y: f64) -> (SurfaceGeometry, SurfaceGeometry) {
        let drum: SurfaceGeometry = ogeom_geom::CylinderSurface::new(
            ogeom_math::Cylinder::new(Frame::WORLD, 20.0, T).unwrap(),
            (0.0, 1.0),
        )
        .unwrap()
        .into();
        let side = BSplineCurve::new(
            KnotVector::clamped_uniform(1, 2).unwrap(),
            vec![Point::new(-1.0, y, 0.0), Point::new(1.0, y, 0.0)],
            T,
        )
        .unwrap();
        let swept: SurfaceGeometry =
            ExtrusionSurface::over(side.into(), direction(0.0, 0.0, 1.0), (0.0, 1.0))
                .unwrap()
                .into();
        (drum, swept)
    }

    /// A straight B-spline side `d` inside the drum's wall crosses it on
    /// the rulings at `x = ±sqrt(2rd − d²)`, however shallow the chord,
    /// each exact, on both surfaces at its traces and over the drum's
    /// height; asked the other way round, the same rulings.
    #[test]
    fn a_spline_side_crosses_a_drum_on_rulings_at_its_crossings() {
        for d in [1e-2, 1e-3, 1e-5] {
            let (drum, swept) = drum_and_side(20.0 - d);
            let want = d.mul_add(-d, 40.0 * d).sqrt();
            for (a, b) in [(&drum, &swept), (&swept, &drum)] {
                let found = super::drum_extrusion(a, b, T).unwrap();
                assert_eq!(found.len(), 2, "dip {d}");
                for (section, sign) in found.iter().zip([-1.0, 1.0]) {
                    let Curve::Line(line) = &section.curve else {
                        panic!("dip {d}: a ruling is a line");
                    };
                    let at = line.axis().location;
                    assert!(
                        (at.x - sign * want).abs() <= 1e-10,
                        "dip {d}: x {} against {}",
                        at.x,
                        sign * want
                    );
                    assert!(section.exact && section.tolerance == 0.0);
                    let (lo, hi) = section.curve.domain();
                    assert!(
                        lo.abs() <= 1e-12 && (hi - 1.0).abs() <= 1e-12,
                        "dip {d}: {lo}..{hi}"
                    );
                    let (Some(on_a), Some(on_b)) = (&section.on_a, &section.on_b) else {
                        panic!("dip {d}: a section without its traces");
                    };
                    for t in [lo, f64::midpoint(lo, hi), hi] {
                        let p = section.curve.point_at(t, T).unwrap();
                        let (ua, ub) = (on_a.point_at(t, T).unwrap(), on_b.point_at(t, T).unwrap());
                        assert!(a.point_at(ua.x, ua.y, T).unwrap().distance(p) <= 1e-11);
                        assert!(b.point_at(ub.x, ub.y, T).unwrap().distance(p) <= 1e-11);
                    }
                }
            }
        }
    }

    /// A side touching the drum, or ending on it, is no crossing to take in
    /// closed form: both are left to the marcher and the contact handling,
    /// and the side clear of it misses.
    #[test]
    fn a_spline_side_touching_a_drum_is_left_to_the_marcher() {
        for y in [20.0, 399.0_f64.sqrt()] {
            let (drum, swept) = drum_and_side(y);
            assert!(super::drum_extrusion(&drum, &swept, T).is_none(), "y {y}");
        }
        let (drum, swept) = drum_and_side(20.5);
        assert!(super::drum_extrusion(&drum, &swept, T).unwrap().is_empty());
    }
}
