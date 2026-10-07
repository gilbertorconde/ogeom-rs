//! Mass properties of a trimmed face, integrated in its own chart.
//!
//! A face is a region of its surface's `(u, v)` chart bounded by its
//! pcurves. Green's theorem turns the integral over that region into one
//! round its boundary:
//!
//! ```text
//! ∬ f du dv = ∮ F dv,   F(u, v) = ∫ f(s, v) ds from u_ref to u
//! ```
//!
//! so each quadrature node on a boundary loop carries an inner integral
//! along `u` from a fixed `u_ref` to the node. Every sample of that inner
//! integral is a point of the surface with its `n dA` and a weight, which
//! is what the mass accumulators take: the region is never meshed and its
//! boundary is the exact pcurves, whatever surface and whatever trim.
//!
//! The surfaces taken are those whose integrands are trigonometric
//! polynomials or polynomials in the chart: the analytic ones, and lines
//! and conics swept or revolved. Panels break at a pcurve's knots and at
//! every quarter turn, where the ten-point Gauss rule integrates such an
//! integrand to rounding, and the whole is run again on panels twice as
//! fine until two runs agree to a part in ten billion on what the caller
//! sums: a volume's runs on the flux integrands of the volume and its
//! moments, an area's also on `|n dA|`. A short boundary panel on an
//! analytic surface takes as few points as its error bound allows, and the
//! inner integrals there, exact on quarter turns, are not refined between
//! runs. Near a fold of a spline surface `|n|` dips almost to nothing
//! between two nodes of an inner panel, and no uniform split takes
//! `|n dA|` there quickly: an area's inner panel across such a dip is
//! graded towards its bottom instead.
//!
//! A pcurve that lifts off its edge's own curve by more than a confusion
//! distance, while its neighbour across the edge runs along the curve,
//! would leave a slit in the boundary the divergence theorem closes over.
//! The strip between the lifted pcurve and the curve is integrated with
//! the face as a ruled surface, and a volume closes the strips' ends onto
//! the curves' ends at each corner, so the faces still close; an area
//! counts a strip as surface the face lacks, or takes it away where it
//! lies back over the face.

use core::f64::consts::FRAC_1_SQRT_2;
use std::sync::OnceLock;

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{
    Curve, Curve2d as _, Curve3d as _, PlanarCurve, Surface as _, SurfaceGeometry,
    Transformable as _,
};
use ogeom_math::{Point, Point2, Vector, Vector2, gauss_legendre_rule};
use ogeom_topo::{EdgeRepr, Model, NodeData, Orientation, Shape};

/// One pcurve piece of a boundary loop, walked from `t0` to `t1` and moved
/// by `shift` in the chart: a whole number of periods where the walk
/// crossed a periodic surface's join.
struct Segment {
    edge: Shape,
    curve: PlanarCurve,
    t0: f64,
    t1: f64,
    shift: Vector2,
    /// The edge's own curve where the pcurve lifts off it by more than the
    /// edge states, for the strip between the two.
    ribbon: Option<Ribbon>,
}

/// An edge's curve, placed, which a piece's lifted pcurve runs beside
/// instead of along. The face's neighbour across the edge is bounded by the
/// curve, so the strip between the lifted pcurve and the curve is
/// integrated with the face as a ruled surface, each ruling from a point of
/// the lifted pcurve to its nearest point on the curve. The face, its
/// strips and its neighbours' then close: what is left open is a gap no
/// wider than the strip at each end of it.
struct Ribbon {
    curve: Curve,
    range: (f64, f64),
    /// The curve's parameter at the piece's start and end, the guesses the
    /// nearest points are sought from.
    ends: (f64, f64),
    /// How far the piece was found to stand off the curve at most.
    reach: f64,
    /// How far a nearest point may lie before it is taken for a wrong one.
    limit: f64,
    /// Whether an area counts the strip: where the neighbour across the
    /// edge is bounded by another curve than the piece, a curve with a
    /// closed form it stands off by more than a hundred confusion
    /// distances or any curve past what the edge states. A fitted curve's
    /// own pieces within its tolerance bound the face as they lie, and the
    /// strip, standing up from the surface as much as lying in it, is
    /// counted by the volume alone, which needs it closed.
    area: bool,
    /// The curve's parameter nearest the lifted piece at [`FEET`] even
    /// shares of the piece and its ends, each sought from the last, so a
    /// nearest point is sought from where the strip's feet run and does
    /// not jump to another stretch of the curve passing as near.
    feet: Vec<f64>,
    /// The piece's parameters where the feet cross one of the curve's
    /// knots or come to rest at one of its ends: there the strip's
    /// integrand changes piece, and the panels break.
    breaks: Vec<f64>,
    /// The strip's lobes in walking order, each as the piece's parameter
    /// where it ends and whether it lies back over the face, which then
    /// counts its area already: an area takes such a lobe away, and adds
    /// one lying out past the face's edge or standing up from it. A lobe
    /// ends where the piece crosses the curve and the rulings turn round,
    /// or where their feet come to rest at an end of the curve; the last
    /// ends with the piece. Read for an area only, the first time one is
    /// asked; a volume takes the strip's `n dA` the same way either side.
    lobes: OnceLock<Vec<(f64, bool)>>,
}

/// How many even steps along a piece its strip's feet are traced at.
const FEET: u32 = 64;

impl Ribbon {
    /// Where to seek the nearest point from at `share` of the piece: read
    /// off the traced feet.
    fn guess(&self, share: f64) -> f64 {
        let last = self.feet.len().saturating_sub(1);
        if last == 0 {
            return self.ends.0 + (self.ends.1 - self.ends.0) * share;
        }
        let x = share.clamp(0.0, 1.0) * f64::from(FEET);
        // `x` lies in `[0, FEET]`, so its floor is a valid index.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let k = (x.floor() as usize).min(last - 1);
        #[allow(clippy::cast_precision_loss)]
        let f = x - k as f64;
        self.feet[k] + (self.feet[k + 1] - self.feet[k]) * f
    }

    /// The curve's parameter nearest `at`, by Newton's method from `guess`,
    /// held to the curve's range.
    fn foot(&self, at: Point, guess: f64, tol: Tolerances) -> OgeomResult<f64> {
        foot_on(&self.curve, self.range, at, guess, tol)
    }

    /// Traces [`Ribbon::feet`] along a piece, each foot sought from the
    /// last ones carried on, or from the curve's nearest sample where that
    /// lands further than twice `seen`, the widest gap sampled so far; and
    /// finds [`Ribbon::breaks`]. Returns how far the piece stands off the
    /// curve at most, or `None` where it stands further than `widest`.
    fn trace(
        &mut self,
        placed: &SurfaceGeometry,
        segment: &Segment,
        seen: f64,
        widest: f64,
        tol: Tolerances,
    ) -> OgeomResult<Option<f64>> {
        let lifted = |share: f64| -> OgeomResult<Point> {
            let t = segment.t0 + (segment.t1 - segment.t0) * share;
            let (at, _) = segment.at(t, tol)?;
            let at = into_domain(placed, at);
            placed.point_at(at.x, at.y, tol)
        };
        let off = |tau: f64, point: Point| -> OgeomResult<f64> {
            Ok(self.curve.point_at(tau, tol)?.distance(point))
        };
        let mut feet: Vec<f64> = Vec::with_capacity(FEET as usize + 1);
        let mut points = Vec::with_capacity(FEET as usize + 1);
        let mut reach: f64 = 0.0;
        for k in 0..=FEET {
            let point = lifted(f64::from(k) / f64::from(FEET))?;
            let guess = match feet.as_slice() {
                [] => self.ends.0,
                [one] => *one,
                [.., a, b] => 2.0 * b - a,
            };
            let mut tau = self.foot(point, guess, tol)?;
            let mut gap = off(tau, point)?;
            if gap > 2.0 * seen.max(reach) {
                let (from, _) = nearest_on(&self.curve, self.range, point, tol)?;
                tau = self.foot(point, from, tol)?;
                gap = off(tau, point)?;
            }
            if gap > widest {
                return Ok(None);
            }
            reach = reach.max(gap);
            feet.push(tau);
            points.push(point);
        }
        // Where the feet pass a knot `kappa` of the curve, or leave or come
        // to rest at an end of its range: the piece's point there lies on
        // the plane through the curve's point at `kappa` square to it.
        let (lo, hi) = (
            self.range.0.min(self.range.1),
            self.range.0.max(self.range.1),
        );
        let mut marks = spline_knots(&self.curve).unwrap_or_default();
        marks.retain(|k| *k > lo && *k < hi);
        marks.extend([lo, hi]);
        let mut breaks = Vec::new();
        for kappa in marks {
            let d = self.curve.derivatives_at(kappa, 1, tol)?;
            let (base, tangent) = (Point::ORIGIN + d[0], d[1]);
            let side = |tau: f64| {
                if tau == kappa {
                    0.0
                } else {
                    (tau - kappa).signum()
                }
            };
            for k in 0..FEET as usize {
                if side(feet[k]) == side(feet[k + 1]) {
                    continue;
                }
                let plane = |point: Point| (point - base).dot(tangent);
                #[allow(clippy::cast_precision_loss)]
                let (mut a, mut b) = (k as f64 / f64::from(FEET), (k + 1) as f64 / f64::from(FEET));
                let (mut fa, mut fb) = (plane(points[k]), plane(points[k + 1]));
                if fa * fb > 0.0 {
                    continue;
                }
                // Regula falsi, the end held twice in a row halved.
                let mut held = 0i8;
                for _ in 0..60 {
                    if fa == fb || b - a <= 1e-15 {
                        break;
                    }
                    let c = (a * fb - b * fa) / (fb - fa);
                    let fc = plane(lifted(c)?);
                    if fc == 0.0 {
                        (a, b) = (c, c);
                        break;
                    }
                    if fc * fa < 0.0 {
                        (b, fb) = (c, fc);
                        if held == -1 {
                            fa *= 0.5;
                        }
                        held = -1;
                    } else {
                        (a, fa) = (c, fc);
                        if held == 1 {
                            fb *= 0.5;
                        }
                        held = 1;
                    }
                }
                breaks.push(segment.t0 + (segment.t1 - segment.t0) * f64::midpoint(a, b));
            }
        }
        self.feet = feet;
        self.breaks = breaks;
        Ok(Some(reach))
    }

    /// The point of the curve nearest `at`, the curve's tangent there, how
    /// fast that point moves along the curve as `at` moves by `moving`, and
    /// whether it is held at an end of the edge's range, where it stands
    /// still. `guess` starts the search.
    fn nearest(
        &self,
        at: Point,
        moving: Vector,
        guess: f64,
        tol: Tolerances,
    ) -> OgeomResult<(Point, Vector, f64, bool)> {
        let (lo, hi) = (
            self.range.0.min(self.range.1),
            self.range.0.max(self.range.1),
        );
        let tau = self.foot(at, guess, tol)?;
        let d = self.curve.derivatives_at(tau, 2, tol)?;
        let point = Point::ORIGIN + d[0];
        if point.distance(at) > self.limit {
            ogeom_core::ogeom_bail!(
                NotDone,
                "a boundary point found no nearest point on its edge's curve"
            );
        }
        let slope = d[1].dot(d[1]) + (point - at).dot(d[2]);
        let held = tau <= lo || tau >= hi;
        // Differentiating `(C(tau) - at) . C'(tau) = 0` along the piece.
        let rate = if held || slope <= 0.0 {
            0.0
        } else {
            moving.dot(d[1]) / slope
        };
        Ok((point, d[1], rate, held))
    }
}

/// The parameter of `curve` over `range` nearest `at`, by Newton's method
/// from `guess`, held to the range.
fn foot_on(
    curve: &Curve,
    range: (f64, f64),
    at: Point,
    guess: f64,
    tol: Tolerances,
) -> OgeomResult<f64> {
    let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
    let mut tau = guess.clamp(lo, hi);
    for _ in 0..30 {
        let d = curve.derivatives_at(tau, 2, tol)?;
        let off = (Point::ORIGIN + d[0]) - at;
        let slope = d[1].dot(d[1]) + off.dot(d[2]);
        if slope <= 0.0 {
            break;
        }
        let next = (tau - off.dot(d[1]) / slope).clamp(lo, hi);
        let moved = (next - tau).abs();
        tau = next;
        if moved <= 1e-15 * (1.0 + tau.abs()) {
            break;
        }
    }
    Ok(tau)
}

impl Segment {
    fn at(&self, t: f64, tol: Tolerances) -> OgeomResult<(Point2, Vector2)> {
        let d = self.curve.derivatives_at(t, 1, tol)?;
        Ok((Point2::ORIGIN + d[0] + self.shift, d[1]))
    }
}

/// A face as its surface and the closed chart loops that bound it, each
/// with the sign that makes its enclosed region count positive for the
/// face's outer boundary and negative for a hole.
pub(crate) struct ChartFace {
    surface: SurfaceGeometry,
    loops: Vec<(Vec<Segment>, f64)>,
    /// The face's orientation: `-1` where it uses its surface reversed.
    sign: f64,
    /// Where the inner integrals start.
    u_ref: f64,
    /// The chart's size, for weighing a miss against it.
    scale: f64,
    /// The surface's own knot lines in `u` and in `v`, where panels break:
    /// read once, not once per inner integral.
    knot_lines: (Vec<f64>, Vec<f64>),
    /// The straight chart steps closing each junction where one piece ends
    /// short of where the next begins, as their start, their step and their
    /// loop's sign. A fitted pcurve may end anywhere within its edge's
    /// tolerance of the vertex, and an open loop's integral depends on
    /// where the inner integrals start: a step of `dv` costs the strip from
    /// `u_ref` across, the chart's whole width on a face wrapping round its
    /// axis. Closed by these steps, every loop measures the region it
    /// bounds whatever `u_ref` is.
    bridges: Vec<(Point2, Vector2, f64)>,
    /// Where a piece with a strip meets the next, the polygon closing the
    /// strips' ends onto the curves' ends (see [`corner`]), and the sign
    /// its loop counts with.
    corners: Vec<(Vec<Point>, f64)>,
}

/// A face's chart loops, or `None` where they cannot be had exactly: a
/// surface whose integrand is not a trigonometric polynomial, an edge
/// without a pcurve on the face or whose pcurve strays from it further than
/// [`BESIDE`] allows, a placement that scales, or a loop whose pieces do
/// not meet.
pub(crate) fn chart_face(model: &Model, face: &Shape, tol: Tolerances) -> Option<ChartFace> {
    loops_of(model, face, tol).ok().flatten()
}

/// Whether a surface's integrands are trigonometric polynomials in its
/// chart, or polynomials: the analytic surfaces, and a line or conic swept
/// or revolved. A spline's are neither, and a rule that integrates them to
/// rounding costs more than the mesh.
fn integrable(surface: &SurfaceGeometry) -> bool {
    fn conic(curve: &Curve) -> bool {
        match curve {
            Curve::Line(_)
            | Curve::Circle(_)
            | Curve::Ellipse(_)
            | Curve::Hyperbola(_)
            | Curve::Parabola(_) => true,
            Curve::Trimmed(trimmed) => conic(trimmed.basis()),
            _ => false,
        }
    }
    match surface {
        SurfaceGeometry::Plane(_)
        | SurfaceGeometry::Cylinder(_)
        | SurfaceGeometry::Cone(_)
        | SurfaceGeometry::Sphere(_)
        | SurfaceGeometry::Torus(_) => true,
        // A spline, and a spline curve swept or revolved, is a polynomial
        // (or rational) piece by knot span: taken with the panels broken
        // at its knots.
        SurfaceGeometry::BSpline(_) => true,
        SurfaceGeometry::Extrusion(e) => conic(e.curve()) || spline_knots(e.curve()).is_some(),
        SurfaceGeometry::Revolution(r) => conic(r.curve()) || spline_knots(r.curve()).is_some(),
        _ => false,
    }
}

/// A spline curve's distinct knots, through a trim.
fn spline_knots(curve: &Curve) -> Option<Vec<f64>> {
    match curve {
        Curve::BSpline(b) => Some(b.knots().distinct().into_iter().map(|(k, _)| k).collect()),
        Curve::Trimmed(t) => spline_knots(t.basis()),
        _ => None,
    }
}

/// Where a surface's integrand changes piece: its knot lines in `u` and
/// in `v`.
fn knot_lines(surface: &SurfaceGeometry) -> (Vec<f64>, Vec<f64>) {
    let distinct = |k: &ogeom_math::KnotVector| -> Vec<f64> {
        k.distinct().into_iter().map(|(x, _)| x).collect()
    };
    match surface {
        SurfaceGeometry::BSpline(b) => (distinct(b.u_knots()), distinct(b.v_knots())),
        SurfaceGeometry::Extrusion(e) => (spline_knots(e.curve()).unwrap_or_default(), Vec::new()),
        SurfaceGeometry::Revolution(r) => (Vec::new(), spline_knots(r.curve()).unwrap_or_default()),
        _ => (Vec::new(), Vec::new()),
    }
}

fn loops_of(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<ChartFace>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    if !integrable(surface) {
        return Ok(None);
    }
    let Some(walked) = walked(model, face, true, tol)? else {
        return Ok(None);
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    let u_ref = if surface.is_periodic_u() {
        walked.lo.x
    } else {
        walked.lo.x.clamp(u0, u1)
    };
    let knot_lines = knot_lines(&walked.placed);
    let period = Vector2::new(
        if surface.is_periodic_u() {
            u1 - u0
        } else {
            0.0
        },
        if surface.is_periodic_v() {
            v1 - v0
        } else {
            0.0
        },
    );
    let Walked {
        placed,
        handedness,
        loops,
        scale,
        ..
    } = walked;
    let mut bridges = Vec::new();
    for (segments, region) in &loops {
        for (k, segment) in segments.iter().enumerate() {
            let next = &segments[(k + 1) % segments.len()];
            let (end, _) = segment.at(segment.t1, tol)?;
            let (start, _) = next.at(next.t0, tol)?;
            // A loop wrapping round a periodic chart closes a whole period
            // from where it began; only what is left over is a step.
            let gap = start - end;
            let step = gap - fold(gap, period);
            if step.y != 0.0 {
                bridges.push((end, step, *region));
            }
        }
    }
    let sign = handedness
        * if face.orientation() == Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
    let mut corners = Vec::new();
    for (segments, region) in &loops {
        for (k, segment) in segments.iter().enumerate() {
            let next = &segments[(k + 1) % segments.len()];
            if let Some(polygon) = corner(&placed, segment, next, period, tol) {
                corners.push((polygon, sign * region));
            }
        }
    }
    Ok(Some(ChartFace {
        knot_lines,
        bridges,
        corners,
        surface: placed,
        loops,
        sign,
        u_ref,
        scale,
    }))
}

/// The polygon closing a face's boundary where `segment` meets `next`
/// and either has a strip, or `None` where neither has one.
///
/// Walked as the face's boundary runs, the strips end in rulings from the
/// lifted pieces' ends to their feet on the curves, and the face's own
/// boundary steps from one lifted end to the other; the neighbours across
/// the two edges end on the curves' own ends instead. The polygon runs
/// from the first strip's last foot along its curve to the curve's end, on
/// to the next curve's start and along it to the next strip's first foot,
/// then back up that ruling, across the step and down the first ruling.
/// Integrated with the face it leaves the face ending on the curves' ends,
/// as its neighbours do. A piece without a strip lies along its curve and
/// stands for it, its lifted end for the curve's end.
fn corner(
    placed: &SurfaceGeometry,
    segment: &Segment,
    next: &Segment,
    period: Vector2,
    tol: Tolerances,
) -> Option<Vec<Point>> {
    const PIECES: u32 = 4;
    if segment.ribbon.is_none() && next.ribbon.is_none() {
        return None;
    }
    let lift = |at: Point2| -> Option<Point> {
        let at = into_domain(placed, at);
        placed.point_at(at.x, at.y, tol).ok()
    };
    let (end, _) = segment.at(segment.t1, tol).ok()?;
    let (start, _) = next.at(next.t0, tol).ok()?;
    let (last, first) = (lift(end)?, lift(start)?);
    let along = |ribbon: &Ribbon, from: f64, to: f64| -> Option<Vec<Point>> {
        (0..=PIECES)
            .map(|j| {
                let tau = from + (to - from) * f64::from(j) / f64::from(PIECES);
                ribbon.curve.point_at(tau, tol).ok()
            })
            .collect()
    };
    let mut polygon = Vec::new();
    match &segment.ribbon {
        Some(ribbon) => {
            let foot = ribbon.foot(last, ribbon.guess(1.0), tol).ok()?;
            polygon.extend(along(ribbon, foot, ribbon.ends.1)?);
        }
        None => polygon.push(last),
    }
    if let Some(ribbon) = &next.ribbon {
        let foot = ribbon.foot(first, ribbon.guess(0.0), tol).ok()?;
        polygon.extend(along(ribbon, ribbon.ends.0, foot)?);
    }
    polygon.push(first);
    // Back across the step from the next piece's start to this one's end.
    let gap = start - end;
    let step = gap - fold(gap, period);
    for j in 1..PIECES {
        polygon.push(lift(end + step * (1.0 - f64::from(j) / f64::from(PIECES)))?);
    }
    if segment.ribbon.is_some() {
        polygon.push(last);
    }
    Some(polygon)
}

/// A piece's strip split into lobes where the piece crosses its edge's
/// curve (see [`Ribbon::lobes`]), and where its rulings' feet come to rest
/// at an end of the curve. A lobe lies back over the face where its
/// rulings, read at a few points, run on the whole within half a right
/// angle of straight into the face. Each lobe is one answer and its ends
/// break the panels, so the strip's area is smooth across every panel and
/// every run counts it the same way.
fn lobes(
    placed: &SurfaceGeometry,
    segment: &Segment,
    ribbon: &Ribbon,
    region: f64,
    tol: Tolerances,
) -> Vec<(f64, bool)> {
    // The face lies to the left of the walk in the chart where the region
    // counts positive.
    let turn = region * (segment.t1 - segment.t0).signum();
    // The ruling at `t`, the unit direction into the face there, and
    // whether the ruling's foot is held at an end of the curve.
    let ruling = |t: f64| -> Option<(Vector, Vector, bool)> {
        let (at, d) = segment.at(t, tol).ok()?;
        let at = into_domain(placed, at);
        let (point, du, dv) = placed.point_d1_at(at.x, at.y, tol).ok()?;
        let share = if segment.t1 == segment.t0 {
            0.0
        } else {
            (t - segment.t0) / (segment.t1 - segment.t0)
        };
        let guess = ribbon.guess(share);
        let (target, _, _, held) = ribbon
            .nearest(point, du * d.x + dv * d.y, guess, tol)
            .ok()?;
        let inward = (dv * d.x - du * d.y) * turn;
        let length = inward.magnitude();
        (length > 0.0).then(|| (target - point, inward / length, held))
    };
    const SAMPLES: u32 = 32;
    let at = |k: u32| segment.t0 + (segment.t1 - segment.t0) * f64::from(k) / f64::from(SAMPLES);
    let mut ends = Vec::new();
    // The last sample read: where, its ruling, and whether its foot was held.
    let mut last: Option<(f64, Vector, bool)> = None;
    for k in 0..=SAMPLES {
        let Some((across, _, held)) = ruling(at(k)) else {
            continue;
        };
        if let Some((from, before, was_held)) = last {
            let turned = across.dot(before) < 0.0;
            if turned || held != was_held {
                // Where the rulings shrink to nothing and turn round, or
                // where the foot comes to rest, by bisection.
                let changed = |x: &(Vector, Vector, bool)| {
                    if turned {
                        x.0.dot(before) < 0.0
                    } else {
                        x.2 != was_held
                    }
                };
                let (mut lo, mut hi) = (from, at(k));
                for _ in 0..50 {
                    let mid = f64::midpoint(lo, hi);
                    match ruling(mid) {
                        Some(x) if !changed(&x) => lo = mid,
                        Some(_) => hi = mid,
                        None => break,
                    }
                }
                let end = f64::midpoint(lo, hi);
                if end != segment.t0 && end != segment.t1 {
                    ends.push(end);
                }
            }
        }
        let before = match last {
            Some((_, before, _)) if across.magnitude() == 0.0 => before,
            _ => across,
        };
        last = Some((at(k), before, held));
    }
    ends.push(segment.t1);
    let mut start = segment.t0;
    let mut out = Vec::with_capacity(ends.len());
    for end in ends {
        let (mut into, mut width) = (0.0, 0.0);
        for k in 1..=5 {
            if let Some((across, inward, _)) = ruling(start + (end - start) * f64::from(k) / 6.0) {
                into += across.dot(inward);
                width += across.magnitude();
            }
        }
        out.push((end, into > FRAC_1_SQRT_2 * width));
        start = end;
    }
    out
}

/// `1` for a placement that keeps a chart's metric and handedness, `-1` for
/// one that keeps its metric and reflects it, `None` for one that scales:
/// the pcurves are the unplaced surface's, and a scale changes the metric
/// under them.
pub(crate) fn rigid_handedness(placement: &ogeom_math::Transform) -> Option<f64> {
    match placement.kind() {
        ogeom_math::TransformKind::Identity
        | ogeom_math::TransformKind::Translation
        | ogeom_math::TransformKind::Rotation => Some(1.0),
        ogeom_math::TransformKind::PlaneMirror | ogeom_math::TransformKind::PointMirror => {
            Some(-1.0)
        }
        _ => None,
    }
}

/// Whether a pcurve of the face runs beside its edge's curve (see
/// [`Fit::Beside`]) further than a confusion distance: the region its
/// pcurves bound then leaves a strip open against the neighbour across
/// that edge, which only the chart loops close. Only edges stating that
/// much are read.
pub(crate) fn runs_wide_of_an_edge(
    model: &Model,
    face: &Shape,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(false);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(false);
    };
    let mut placed = None;
    for wire in model.ordered_children_of(face)? {
        for edge in model.ordered_children_of(&wire)? {
            let Some(edge_data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                continue;
            };
            if edge_data.tolerance.get() <= tol.confusion() {
                continue;
            }
            let (ids, range) = match edge_data.pcurve_for(data.surface, edge.location()) {
                Some(EdgeRepr::PCurve { curve, range, .. }) => (vec![*curve], *range),
                Some(EdgeRepr::Seam {
                    forward,
                    reversed,
                    range,
                    ..
                }) => (vec![*forward, *reversed], *range),
                _ => continue,
            };
            if placed.is_none() {
                let placement = face.transform(model.datums())?;
                if rigid_handedness(&placement).is_none() {
                    return Ok(false);
                }
                placed = Some(surface.clone().transformed(&placement, tol)?);
            }
            let Some(placed) = placed.as_ref() else {
                continue;
            };
            for id in ids {
                let Some(curve) = model.geometry().pcurve(id) else {
                    continue;
                };
                let segment = Segment {
                    edge: edge.clone(),
                    curve: curve.clone(),
                    t0: range.0,
                    t1: range.1,
                    shift: Vector2::new(0.0, 0.0),
                    ribbon: None,
                };
                if let Fit::Beside(ribbon) = fit_to_edge(model, &edge, placed, &segment, tol)?
                    && ribbon.reach > tol.confusion()
                {
                    return Ok(true);
                }
            }
        }
    }
    Ok(false)
}

/// A face's boundary walked into closed chart loops.
struct Walked {
    /// The surface, placed as the face is.
    placed: SurfaceGeometry,
    /// `-1` where the placement reflects: the placed chart's own normal
    /// then points against the face's.
    handedness: f64,
    /// Each loop, with the sign that turns it round its region: `1` where
    /// the face lies to the left of the walk, `-1` where to the right.
    loops: Vec<(Vec<Segment>, f64)>,
    /// The chart's lower corner and its size.
    lo: Point2,
    scale: f64,
}

/// Points along a face's boundary, each with its edge and the chart
/// direction the face lies in from there, read off the walked loops'
/// windings. `None` where the boundary cannot be walked into closed loops.
pub(crate) fn material_sides(
    model: &Model,
    face: &Shape,
    tol: Tolerances,
) -> Option<Vec<(Shape, Point2, Vector2)>> {
    let walked = walked(model, face, false, tol).ok()??;
    let mut out = Vec::new();
    for (segments, region) in &walked.loops {
        for segment in segments {
            for k in 1..=4 {
                let t = segment.t0 + (segment.t1 - segment.t0) * f64::from(k) / 5.0;
                let (at, d) = segment.at(t, tol).ok()?;
                let heading = if segment.t1 < segment.t0 { -d } else { d };
                out.push((segment.edge.clone(), at, heading.perpendicular() * *region));
            }
        }
    }
    Some(out)
}

/// A face's loops in its chart, or `None` where they cannot be had: an
/// edge without a pcurve on the face, a placement that scales, or a loop
/// whose pieces do not meet; and where `strict`, a pcurve straying from
/// its edge.
fn walked(
    model: &Model,
    face: &Shape,
    strict: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Walked>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    // The pcurves are the unplaced surface's; a rigid placement keeps its
    // chart, a scaling one would change the metric under them.
    let placement = face.transform(model.datums())?;
    let Some(handedness) = rigid_handedness(&placement) else {
        return Ok(None);
    };
    let placed = surface.clone().transformed(&placement, tol)?;
    let ((u0, u1), (v0, v1)) = surface.domain();
    let period = Vector2::new(
        if surface.is_periodic_u() {
            u1 - u0
        } else {
            0.0
        },
        if surface.is_periodic_v() {
            v1 - v0
        } else {
            0.0
        },
    );

    let mut loops = Vec::new();
    for wire in model.ordered_children_of(face)? {
        let Some(segments) = walk(model, data.surface, &placed, &wire, period, strict, tol)? else {
            return Ok(None);
        };
        loops.push(segments);
    }
    if loops.is_empty() {
        return Ok(None);
    }

    // The chart's extent, from a few points of every piece.
    let (mut lo, mut hi) = (
        Point2::new(f64::INFINITY, f64::INFINITY),
        Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for segments in &loops {
        for segment in segments {
            for k in 0..=4 {
                let t = segment.t0 + (segment.t1 - segment.t0) * f64::from(k) / 4.0;
                let (p, _) = segment.at(t, tol)?;
                lo = Point2::new(lo.x.min(p.x), lo.y.min(p.y));
                hi = Point2::new(hi.x.max(p.x), hi.y.max(p.y));
            }
        }
    }
    let scale = (hi.x - lo.x).max(hi.y - lo.y);
    if !scale.is_finite() || scale <= 0.0 {
        return Ok(None);
    }

    // Each loop must close, piece to piece, to a millionth of the chart, or
    // to what the two edges meeting there own: fitted sections each miss
    // their junction by up to their stated tolerance, possibly opposite
    // ways, and the vertex there holds both ends within its tolerance, so
    // they may stand its diameter apart. The chart reads that slack
    // through the surface's stretch there.
    let reach = scale * 1e-6 + tol.parametric();
    let owned = |edge: &Shape| {
        let own = model
            .node(edge)
            .and_then(|n| n.data().as_edge())
            .map_or(0.0, |d| d.tolerance.get());
        let ends = crate::edge_vertices(model, edge)
            .ok()
            .flatten()
            .map_or(0.0, |(a, b)| {
                [a, b]
                    .iter()
                    .filter_map(|v| model.node(v).and_then(|n| n.data().as_vertex()))
                    .map(|d| d.tolerance.get())
                    .fold(0.0, f64::max)
            });
        (own, ends)
    };
    for segments in &loops {
        for (k, segment) in segments.iter().enumerate() {
            let next = &segments[(k + 1) % segments.len()];
            let (end, _) = segment.at(segment.t1, tol)?;
            let (start, _) = next.at(next.t0, tol)?;
            let (a, b) = (owned(&segment.edge), owned(&next.edge));
            let slack = (a.0 + b.0).max(2.0 * a.1).max(2.0 * b.1);
            let stretch = placed.d1_at(end.x, end.y, tol).map_or(1.0, |(du, dv)| {
                du.magnitude().min(dv.magnitude()).max(tol.confusion())
            });
            if end.distance(start) > reach + slack / stretch {
                return Ok(None);
            }
        }
    }

    // Which loop is the boundary and which the holes: the boundary encloses
    // the most chart, whichever way each was wound. A spline is integrated
    // knot span by knot span: one rule across a closed spline of many spans
    // can miss its area entirely and give the wrong winding.
    let mut areas = Vec::with_capacity(loops.len());
    for segments in &loops {
        let mut area = 0.0;
        for segment in segments {
            let (lo, hi) = (segment.t0.min(segment.t1), segment.t0.max(segment.t1));
            let mut breaks = vec![segment.t0, segment.t1];
            if let PlanarCurve::BSpline(spline) = &segment.curve {
                breaks.extend(
                    spline
                        .knots()
                        .distinct()
                        .into_iter()
                        .map(|(k, _)| k)
                        .filter(|k| *k > lo && *k < hi),
                );
            }
            breaks.sort_by(f64::total_cmp);
            if segment.t1 < segment.t0 {
                breaks.reverse();
            }
            for pair in breaks.windows(2) {
                for (t, w) in gauss_legendre_rule(pair[0], pair[1]) {
                    let (p, d) = segment.at(t, tol)?;
                    area += p.x * d.y * w;
                }
            }
        }
        areas.push(area);
    }
    let Some(outer) = (0..areas.len()).max_by(|a, b| areas[*a].abs().total_cmp(&areas[*b].abs()))
    else {
        return Ok(None);
    };
    let loops = loops
        .into_iter()
        .zip(&areas)
        .enumerate()
        .map(|(k, (segments, area))| {
            let region = if k == outer { 1.0 } else { -1.0 };
            (segments, region * area.signum())
        })
        .collect();
    Ok(Some(Walked {
        placed,
        handedness,
        loops,
        lo,
        scale,
    }))
}

/// A wire's pcurve pieces in walking order, or `None` where one is missing.
///
/// A seam bounds its face twice, once down each side of the chart; which
/// side an occurrence takes is the one continuing the point already walked
/// to, so the walk starts off a seam where it can. A piece that starts a
/// whole period from where the last one ended is moved by that period.
fn walk(
    model: &Model,
    surface: ogeom_topo::SurfaceId,
    placed: &SurfaceGeometry,
    wire: &Shape,
    period: Vector2,
    strict: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Segment>>> {
    let mut edges = model.ordered_children_of(wire)?;
    let is_seam = |e: &Shape| {
        model
            .node(e)
            .and_then(|n| n.data().as_edge())
            .and_then(|d| d.pcurve_for(surface, e.location()))
            .is_some_and(|r| matches!(r, EdgeRepr::Seam { .. }))
    };
    if let Some(start) = edges.iter().position(|e| !is_seam(e)) {
        edges.rotate_left(start);
    }
    let mut segments: Vec<Segment> = Vec::with_capacity(edges.len());
    let mut last: Option<Point2> = None;
    for edge in &edges {
        let Some(repr) = model
            .node(edge)
            .and_then(|n| n.data().as_edge())
            .and_then(|d| d.pcurve_for(surface, edge.location()))
        else {
            return Ok(None);
        };
        let reversed = edge.orientation() == Orientation::Reversed;
        let (ids, range) = match repr {
            EdgeRepr::PCurve { curve, range, .. } => (vec![*curve], *range),
            EdgeRepr::Seam {
                forward,
                reversed: back,
                range,
                ..
            } => {
                let preferred = if reversed {
                    [*back, *forward]
                } else {
                    [*forward, *back]
                };
                (preferred.to_vec(), *range)
            }
            _ => return Ok(None),
        };
        let (t0, t1) = if reversed { (range.1, range.0) } else { range };
        // Of the candidate sides, the one whose start lies nearest the walk,
        // after folding by whole periods.
        let mut best: Option<(f64, Segment)> = None;
        for id in ids {
            let Some(curve) = model.geometry().pcurve(id) else {
                return Ok(None);
            };
            let start = Point2::ORIGIN + curve.derivatives_at(t0, 0, tol)?[0];
            let shift = match last {
                Some(at) => fold(at - start, period),
                None => Vector2::new(0.0, 0.0),
            };
            let miss = last.map_or(0.0, |at| (start + shift).distance(at));
            if best.as_ref().is_none_or(|(held, _)| miss < *held) {
                best = Some((
                    miss,
                    Segment {
                        edge: edge.clone(),
                        curve: curve.clone(),
                        t0,
                        t1,
                        shift,
                        ribbon: None,
                    },
                ));
            }
        }
        let Some((_, mut segment)) = best else {
            return Ok(None);
        };
        if strict {
            match fit_to_edge(model, edge, placed, &segment, tol)? {
                Fit::Along => {}
                Fit::Beside(ribbon) => segment.ribbon = Some(*ribbon),
                Fit::Off => return Ok(None),
            }
        }
        last = Some(segment.at(segment.t1, tol)?.0);
        segments.push(segment);
    }
    if segments.is_empty() {
        return Ok(None);
    }
    Ok(Some(segments))
}

/// How far, in confusion distances, a piece's lifted pcurve on a curve
/// with a closed form may run beside the curve past what the edge states,
/// and still be integrated with the strip between them: a thousandth of a
/// millimetre at millimetre tolerances. The strip's rulings are straight,
/// and a ruled strip of width `w` departs from any smooth surface through
/// both its sides by about `w^2` times the curvature, so what it leaves
/// out of a volume goes as `w^3`.
const BESIDE: f64 = 1e4;

/// How a piece's lifted pcurve lies against its edge's own curve.
enum Fit {
    /// Along it, to a confusion distance; or the edge has no curve of its
    /// own (a pole) to stray from; or further than [`BESIDE`] confusion
    /// distances but within the edge's stated tolerance, which then owns
    /// the gap.
    Along,
    /// Beside it, further than a confusion distance but within the edge's
    /// stated tolerance, or on a curve with a closed form within [`BESIDE`]
    /// confusion distances: a fitted section, or a curve whose edge took a
    /// curve close by for its own, or a chord standing for an arc. The
    /// strip between is integrated with the face.
    Beside(Box<Ribbon>),
    /// Further: the region the pcurve bounds is not the face's, and the
    /// face is left to the mesh, which takes its boundary from the edge.
    Off,
}

/// How a piece's pcurve, lifted through the surface, lies against its
/// edge's own curve. The two curves need not share a parameter, so each
/// lifted point is measured against the nearest point of the edge's curve
/// over its range.
fn fit_to_edge(
    model: &Model,
    edge: &Shape,
    placed: &SurfaceGeometry,
    segment: &Segment,
    tol: Tolerances,
) -> OgeomResult<Fit> {
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .and_then(|d| d.curve3d())
    else {
        return Ok(Fit::Along);
    };
    let Some(curve) = model.geometry().curve(*curve) else {
        return Ok(Fit::Off);
    };
    let curve = curve
        .clone()
        .transformed(&edge.transform(model.datums())?, tol)?;
    let stated = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .map_or(0.0, |d| d.tolerance.get());
    let reach = (tol.confusion() * 100.0).max(stated);
    let lifted = |t: f64| -> OgeomResult<Point> {
        let (at, _) = segment.at(t, tol)?;
        let at = into_domain(placed, at);
        placed.point_at(at.x, at.y, tol)
    };
    let mut seen: f64 = 0.0;
    let mut within_stated = true;
    for k in 1..=5 {
        let t = segment.t0 + (segment.t1 - segment.t0) * f64::from(k) / 6.0;
        let off = nearest(&curve, *range, lifted(t)?, tol)?;
        seen = seen.max(off);
        if off > reach {
            within_stated = false;
            break;
        }
    }
    // A fitted pcurve can meet the curve at every inner sample and stray
    // towards its ends, where it answers to its vertices' tolerance rather
    // than the edge's: the ends ask for a strip but do not refuse one.
    if within_stated {
        for t in [segment.t0, segment.t1] {
            seen = seen.max(nearest(&curve, *range, lifted(t)?, tol)?);
        }
    }
    if seen <= tol.confusion() && within_stated {
        return Ok(Fit::Along);
    }
    // A curve with a closed form carries no fit error of its own, so a
    // pcurve standing off it past what the edge states is another curve's
    // description beside it: two edges a little apart taken for one. The
    // strip is integrated there all the same, to [`BESIDE`]. Within the
    // stated tolerance the strip is the face's to close however wide it
    // is: an edge taking a straight chord or an arc for a curved section
    // states the chord's sag, its neighbour across the edge is bounded by
    // the chord, and so is a fitted curve's neighbour by the fit.
    let exact = matches!(
        curve,
        Curve::Line(_)
            | Curve::Circle(_)
            | Curve::Ellipse(_)
            | Curve::Hyperbola(_)
            | Curve::Parabola(_)
    );
    if !within_stated && !exact {
        return Ok(Fit::Off);
    }
    let widest_strip = if within_stated {
        (tol.confusion() * BESIDE).max(stated)
    } else {
        tol.confusion() * BESIDE
    };
    // Which way the curve runs along the piece: the reading that keeps
    // the piece's points nearer the curve's at the same share of the way,
    // the quarter points included so a closed curve is told apart too.
    let mut ahead = 0.0;
    let mut behind = 0.0;
    for k in 0..=4 {
        let share = f64::from(k) / 4.0;
        let at = lifted(segment.t0 + (segment.t1 - segment.t0) * share)?;
        ahead += at.distance(curve.point_at(range.0 + (range.1 - range.0) * share, tol)?);
        behind += at.distance(curve.point_at(range.1 + (range.0 - range.1) * share, tol)?);
    }
    let ends = if ahead <= behind {
        *range
    } else {
        (range.1, range.0)
    };
    let mut ribbon = Ribbon {
        area: !within_stated || (exact && seen > tol.confusion() * 100.0),
        curve,
        range: *range,
        ends,
        reach: 0.0,
        limit: widest_strip,
        feet: Vec::new(),
        breaks: Vec::new(),
        lobes: OnceLock::new(),
    };
    // Beside the curve all the way, ends included.
    let Some(widest) = ribbon.trace(placed, segment, seen, widest_strip, tol)? else {
        return Ok(if within_stated { Fit::Along } else { Fit::Off });
    };
    ribbon.reach = widest;
    ribbon.limit = (4.0 * widest).max(widest_strip);
    Ok(Fit::Beside(Box::new(ribbon)))
}

/// A chart point brought into the surface's domain where a pcurve fitted
/// along its border strays past it by rounding, a millionth of the span or
/// less; any further is left for the evaluation to refuse.
pub(crate) fn into_domain(surface: &SurfaceGeometry, at: Point2) -> Point2 {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let onto = |x: f64, lo: f64, hi: f64, periodic: bool| {
        let slack = (hi - lo) * 1e-6;
        if periodic || !(lo..=hi).contains(&x) && (x < lo - slack || x > hi + slack) {
            x
        } else {
            x.clamp(lo, hi)
        }
    };
    Point2::new(
        onto(at.x, u0, u1, surface.is_periodic_u()),
        onto(at.y, v0, v1, surface.is_periodic_v()),
    )
}

/// The distance from `target` to `curve` over `range`: the best of a
/// sampling, narrowed by golden sections about it.
fn nearest(curve: &Curve, range: (f64, f64), target: Point, tol: Tolerances) -> OgeomResult<f64> {
    Ok(nearest_on(curve, range, target, tol)?.1)
}

/// The parameter of `curve` over `range` nearest `target`, and how far it
/// lies: the nearest of a sampling, settled by golden section.
fn nearest_on(
    curve: &Curve,
    range: (f64, f64),
    target: Point,
    tol: Tolerances,
) -> OgeomResult<(f64, f64)> {
    const SAMPLES: u32 = 32;
    let at = |t: f64| -> OgeomResult<f64> { Ok(curve.point_at(t, tol)?.distance(target)) };
    let step = (range.1 - range.0) / f64::from(SAMPLES);
    let mut best = (range.0, at(range.0)?);
    for k in 1..=SAMPLES {
        let t = range.0 + step * f64::from(k);
        let d = at(t)?;
        if d < best.1 {
            best = (t, d);
        }
    }
    let (mut a, mut b) = (
        (best.0 - step).max(range.0.min(range.1)),
        (best.0 + step).min(range.0.max(range.1)),
    );
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    for _ in 0..60 {
        let (c, d) = (b - (b - a) * ratio, a + (b - a) * ratio);
        if at(c)? < at(d)? {
            b = d;
        } else {
            a = c;
        }
    }
    let middle = f64::midpoint(a, b);
    let settled = at(middle)?;
    Ok(if settled < best.1 {
        (middle, settled)
    } else {
        best
    })
}

/// The whole number of periods nearest to `gap`, along each periodic
/// direction.
fn fold(gap: Vector2, period: Vector2) -> Vector2 {
    let along = |g: f64, p: f64| if p > 0.0 { (g / p).round() * p } else { 0.0 };
    Vector2::new(along(gap.x, period.x), along(gap.y, period.y))
}

/// The most times the panels are doubled before the face is left to the
/// mesh.
const DOUBLINGS: u32 = 5;

/// Agreement asked of two runs, the second on panels twice as fine,
/// relative to the size of what they integrate.
const AGREE: f64 = 1e-10;

/// A quarter turn: the widest panel on an angular parameter, over which a
/// trigonometric polynomial integrates to rounding under the ten-point
/// rule.
const QUARTER: f64 = core::f64::consts::FRAC_PI_2;

/// The boundary rules a short panel may take, by their number of points.
const ORDERS: [usize; 5] = [3, 4, 5, 6, 7];

/// The Gauss-Legendre rules of [`ORDERS`] on `[-1, 1]`: nodes and weights.
const RULES: [(&[f64], &[f64]); 5] = [
    (
        &[-0.7745966692414834, 0.0, 0.7745966692414834],
        &[0.5555555555555557, 0.8888888888888888, 0.5555555555555557],
    ),
    (
        &[
            -0.8611363115940526,
            -0.33998104358485626,
            0.33998104358485626,
            0.8611363115940526,
        ],
        &[
            0.34785484513745357,
            0.6521451548625464,
            0.6521451548625464,
            0.34785484513745357,
        ],
    ),
    (
        &[
            -0.906179845938664,
            -0.5384693101056831,
            0.0,
            0.5384693101056831,
            0.906179845938664,
        ],
        &[
            0.23692688505618928,
            0.4786286704993663,
            0.5688888888888887,
            0.4786286704993663,
            0.23692688505618928,
        ],
    ),
    (
        &[
            -0.9324695142031519,
            -0.6612093864662645,
            -0.2386191860831969,
            0.2386191860831969,
            0.6612093864662645,
            0.9324695142031519,
        ],
        &[
            0.17132449237917027,
            0.3607615730481387,
            0.46791393457269104,
            0.46791393457269104,
            0.3607615730481387,
            0.17132449237917027,
        ],
    ),
    (
        &[
            -0.9491079123427586,
            -0.7415311855993945,
            -0.4058451513773972,
            0.0,
            0.4058451513773972,
            0.7415311855993945,
            0.9491079123427586,
        ],
        &[
            0.12948496616886973,
            0.27970539148927687,
            0.3818300505051187,
            0.4179591836734693,
            0.3818300505051187,
            0.27970539148927687,
            0.12948496616886973,
        ],
    ),
];

/// The highest frequency taken for a boundary integrand in an angle: the
/// second moments' on a sphere or a torus, a cubic in the point times the
/// `n dA`.
const FREQUENCY: f64 = 6.0;

/// What a boundary panel's rule may miss by, against the panel's own size:
/// below what ten points over a quarter turn are bounded by.
const PANEL_MISS: f64 = 1e-16;

/// The `order`-point Gauss-Legendre rule on `[a, b]`: one of [`ORDERS`],
/// or ten points.
fn rule(order: usize, a: f64, b: f64) -> Vec<(f64, f64)> {
    let (half, middle) = ((b - a) * 0.5, f64::midpoint(a, b));
    match ORDERS.iter().position(|&n| n == order) {
        Some(i) => {
            let (nodes, weights) = RULES[i];
            nodes
                .iter()
                .zip(weights)
                .map(|(x, w)| (middle + half * x, w * half))
                .collect()
        }
        None => gauss_legendre_rule(a, b).to_vec(),
    }
}

/// One sample of a face's integral: the surface point, its outward
/// `Su x Sv`, and the quadrature weight.
type Sample = (Point, Vector, f64);

/// What a face is integrated for, which says which measures two runs must
/// agree on before the finer one is taken.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Measure {
    /// The area and its moments: the runs agree on `|n dA|`, `n dA` and
    /// the flux of the point.
    Area,
    /// The volume and its moments, which integrate `n dA` times a
    /// polynomial in the point: the runs agree on `n dA` and on the flux
    /// integrands of the volume, its first moments and its second, never
    /// on `|n dA|`. Near a fold of the surface `|n|` almost vanishes and
    /// its square root converges slowly, which says nothing about a volume.
    Volume,
}

/// How an inner integral treats a panel across which `|n|` dips.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Folds {
    /// Integrated as any other.
    Ignore,
    /// Integrated as any other, and reported.
    Watch,
    /// Graded towards the dip where there is one (see [`fold_in`]), and
    /// reported.
    Grade,
}

/// The measures two runs are compared by: `|n dA|` (an area's only),
/// `n dA`, the flux of the point, and a volume's first and second moment
/// integrands.
const PROXIES: usize = 17;
type Proxy = [f64; PROXIES];

impl ChartFace {
    /// A point of the face's surface, to take moments from.
    pub(crate) fn anchor(&self, tol: Tolerances) -> OgeomResult<Point> {
        let segment = &self.loops[0].0[0];
        let (at, _) = segment.at(segment.t0, tol)?;
        let at = into_domain(&self.surface, at);
        self.surface.point_at(at.x, at.y, tol)
    }

    /// The face's integral summed into an accumulator from `fresh`: every
    /// sample is handed to `contribute` as the surface point, its `n dA`
    /// with the weight's magnitude folded in, and the weight's sign (an area
    /// takes `|n dA|` times that sign, a volume the product). Each run sums
    /// into an accumulator of its own, and the one that settles is
    /// returned; `None` where none did.
    pub(crate) fn integrate<A>(
        &self,
        measure: Measure,
        reference: Point,
        tol: Tolerances,
        fresh: impl Fn() -> A,
        contribute: impl Fn(&mut A, Point, Vector, f64),
    ) -> Option<A> {
        // An area's first run watches for `|n|` dipping across an inner
        // panel; where it does, that run and every finer one grade such
        // panels towards the dip (see [`fold_in`]).
        let watch = if measure == Measure::Area {
            Folds::Watch
        } else {
            Folds::Ignore
        };
        let (mut held, mut reach, dipped) = self
            .run(1, measure, watch, reference, tol, &mut |_| {})
            .ok()?;
        let folds = if dipped { Folds::Grade } else { Folds::Ignore };
        if dipped {
            (held, reach, _) = self
                .run(1, measure, folds, reference, tol, &mut |_| {})
                .ok()?;
        }
        // The moments are weighed against the volume by the face's reach
        // from the reference, the same for every run: once for the first,
        // twice for the second.
        let reach = reach.sqrt();
        let weigh = |mut proxy: Proxy| {
            if reach > 0.0 {
                for x in &mut proxy[5..8] {
                    *x /= reach;
                }
                for x in &mut proxy[8..] {
                    *x /= reach * reach;
                }
            }
            proxy
        };
        held = weigh(held);
        for doubling in 1..=DOUBLINGS {
            // Each doubling costs twice the last: a cancelled watch is
            // honoured between them, and the caller's own checkpoint then
            // reports it.
            if ogeom_core::progress::checkpoint().is_err() {
                return None;
            }
            let mut sum = fresh();
            let (proxy, _, _) = self
                .run(
                    1 << doubling,
                    measure,
                    folds,
                    reference,
                    tol,
                    &mut |(p, n, w)| {
                        contribute(&mut sum, p, n * w.abs(), w.signum());
                    },
                )
                .ok()?;
            let proxy = weigh(proxy);
            if settled(held, proxy) {
                return Some(sum);
            }
            held = proxy;
        }
        None
    }

    /// The face's samples with every panel split `fine` ways, each handed to
    /// `sink` in a fixed order, the integrals of a few measures over them to
    /// compare runs by, the farthest sample's squared distance from
    /// `reference`, and whether `|n|` dipped across an inner panel where
    /// `folds` watches for it.
    fn run(
        &self,
        fine: u32,
        measure: Measure,
        folds: Folds,
        reference: Point,
        tol: Tolerances,
        sink: &mut dyn FnMut(Sample),
    ) -> OgeomResult<(Proxy, f64, bool)> {
        let mut proxy = [0.0; PROXIES];
        let mut reach = 0.0_f64;
        let mut dipped = false;
        let size = self.scale.max(1.0);
        let mut take = |(p, n, w): Sample| {
            let e = p - reference;
            let flux = e.dot(n) / size;
            match measure {
                Measure::Area => {
                    let row = [n.magnitude(), n.x, n.y, n.z, flux];
                    for (acc, x) in proxy.iter_mut().zip(row) {
                        *acc += x * w;
                    }
                }
                Measure::Volume => {
                    reach = reach.max(e.dot(e));
                    let (q, nq) = ([e.x, e.y, e.z], [n.x, n.y, n.z]);
                    for (acc, x) in proxy[1..5].iter_mut().zip([n.x, n.y, n.z, flux]) {
                        *acc += x * w;
                    }
                    for i in 0..3 {
                        let lift = q[i] * q[i] * nq[i] / size * w;
                        proxy[5 + i] += lift;
                        for j in 0..3 {
                            proxy[8 + 3 * i + j] += lift * q[j];
                        }
                    }
                }
            }
            sink((p, n, w));
        };
        for (segments, region) in &self.loops {
            for segment in segments {
                if let Some(ribbon) = segment
                    .ribbon
                    .as_ref()
                    .filter(|r| measure == Measure::Volume || r.area)
                {
                    let lobes: &[(f64, bool)] = match measure {
                        Measure::Area => ribbon
                            .lobes
                            .get_or_init(|| lobes(&self.surface, segment, ribbon, *region, tol)),
                        Measure::Volume => &[],
                    };
                    let mut breaks = self.outer_breaks(segment, tol)?;
                    breaks.extend(ribbon.breaks.iter().copied());
                    breaks.extend(lobes.iter().map(|(end, _)| *end));
                    breaks.sort_by(f64::total_cmp);
                    if segment.t1 < segment.t0 {
                        breaks.reverse();
                    }
                    breaks.dedup();
                    // A strip is as wide as the piece stands off its curve,
                    // a sliver of the face's own integral, and five points
                    // a panel take it far closer than the runs agree.
                    for pair in breaks.windows(2) {
                        for k in 0..fine {
                            let a = pair[0] + (pair[1] - pair[0]) * f64::from(k) / f64::from(fine);
                            let b =
                                pair[0] + (pair[1] - pair[0]) * f64::from(k + 1) / f64::from(fine);
                            for (t, wt) in rule(5, a, b) {
                                self.strip(segment, ribbon, lobes, t, region * wt, tol, &mut take)?;
                            }
                        }
                    }
                }
                // A straight piece along which `v` does not move adds
                // nothing.
                if let PlanarCurve::Line(_) = segment.curve {
                    let (_, d) = segment.at(segment.t0, tol)?;
                    if d.y.abs() <= 1e-14 * d.x.abs() {
                        continue;
                    }
                }
                let breaks = self.outer_breaks(segment, tol)?;
                for pair in breaks.windows(2) {
                    let order = self.outer_order(segment, pair[0], pair[1], tol)?;
                    for k in 0..fine {
                        let a = pair[0] + (pair[1] - pair[0]) * f64::from(k) / f64::from(fine);
                        let b = pair[0] + (pair[1] - pair[0]) * f64::from(k + 1) / f64::from(fine);
                        for (t, wt) in rule(order, a, b) {
                            let (at, d) = segment.at(t, tol)?;
                            dipped |=
                                self.inner(at, region * wt * d.y, fine, folds, tol, &mut take)?;
                        }
                    }
                }
            }
        }
        // The corners close a volume; an area has nothing to close.
        if measure == Measure::Volume {
            for (polygon, turn) in &self.corners {
                // A fan of triangles from the first corner, each
                // integrated on its unit square collapsed onto it.
                let hub = polygon[0];
                for k in 1..polygon.len() {
                    let (a, b) = (polygon[k], polygon[(k + 1) % polygon.len()]);
                    let n = (a - hub).cross(b - a) * *turn;
                    for (x, wx) in rule(3, 0.0, 1.0) {
                        for (y, wy) in rule(3, 0.0, 1.0) {
                            take((hub + (a - hub) * x + (b - a) * (x * y), n * x, wx * wy));
                        }
                    }
                }
            }
        }
        for &(start, step, region) in &self.bridges {
            for (t, wt) in gauss_legendre_rule(0.0, 1.0) {
                dipped |= self.inner(
                    start + step * t,
                    region * wt * step.y,
                    fine,
                    folds,
                    tol,
                    &mut take,
                )?;
            }
        }
        Ok((proxy, reach, dipped))
    }

    /// The rulings of a piece's strip at `t`, its samples weighted by
    /// `outer`: the boundary weight at `t`, signed as the walk runs. The
    /// face's boundary runs the lifted pcurve one way; the strip, to close
    /// on it, runs it the other, so its `n dA` is `-sign` times the
    /// ruling's sweep along the walk crossed with the ruling.
    #[allow(clippy::too_many_arguments)]
    fn strip(
        &self,
        segment: &Segment,
        ribbon: &Ribbon,
        lobes: &[(f64, bool)],
        t: f64,
        outer: f64,
        tol: Tolerances,
        sink: &mut dyn FnMut(Sample),
    ) -> OgeomResult<()> {
        let (at, d) = segment.at(t, tol)?;
        let at = into_domain(&self.surface, at);
        let (point, du, dv) = self.surface.point_d1_at(at.x, at.y, tol)?;
        let moving = du * d.x + dv * d.y;
        let share = if segment.t1 == segment.t0 {
            0.0
        } else {
            (t - segment.t0) / (segment.t1 - segment.t0)
        };
        let guess = ribbon.guess(share);
        let (target, tangent, rate, _) = ribbon.nearest(point, moving, guess, tol)?;
        let across = target - point;
        let ahead = (segment.t1 - segment.t0).signum();
        let over = lobes
            .iter()
            .find(|(end, _)| (end - t) * ahead >= 0.0)
            .is_some_and(|(_, over)| *over);
        for (s, ws) in rule(3, 0.0, 1.0) {
            let sweep = moving * (1.0 - s) + tangent * (rate * s);
            let (n, w) = (sweep.cross(across), -self.sign * outer * ws);
            // The weight's sign is the area's (see [`Ribbon::lobes`]); the
            // product `n w`, which a volume takes, is the strip's own
            // either way.
            if (w < 0.0) == over {
                sink((point + across * s, n, w));
            } else {
                sink((point + across * s, -n, -w));
            }
        }
        Ok(())
    }

    /// Where a boundary piece's panels break: its own ends, its knots, and
    /// every quarter turn it makes in an angular chart.
    fn outer_breaks(&self, segment: &Segment, tol: Tolerances) -> OgeomResult<Vec<f64>> {
        let (t0, t1) = (segment.t0, segment.t1);
        let mut breaks = vec![t0, t1];
        let (lo, hi) = (t0.min(t1), t0.max(t1));
        let knots = match &segment.curve {
            PlanarCurve::BSpline(spline) => spline
                .knots()
                .distinct()
                .into_iter()
                .map(|(k, _)| k)
                .collect(),
            _ => Vec::new(),
        };
        breaks.extend(knots.into_iter().filter(|k| *k > lo && *k < hi));
        // Where the piece crosses the surface's own knot lines, found on a
        // sampling and settled by bisection.
        let (u_knots, v_knots) = &self.knot_lines;
        if !u_knots.is_empty() || !v_knots.is_empty() {
            const N: usize = 32;
            let mut prev: Option<(f64, Point2)> = None;
            for k in 0..=N {
                #[allow(clippy::cast_precision_loss)]
                let t = t0 + (t1 - t0) * k as f64 / N as f64;
                let (p, _) = segment.at(t, tol)?;
                if let Some((tp, pp)) = prev {
                    for (lines, read) in [(&u_knots, 0usize), (&v_knots, 1usize)] {
                        let coord = |q: Point2| if read == 0 { q.x } else { q.y };
                        for line in lines.iter() {
                            let (a, b) = (coord(pp) - line, coord(p) - line);
                            if a * b >= 0.0 {
                                continue;
                            }
                            let (mut lo_t, mut hi_t, mut f_lo) = (tp, t, a);
                            for _ in 0..60 {
                                let mid = 0.5 * (lo_t + hi_t);
                                let f = coord(segment.at(mid, tol)?.0) - line;
                                if f.signum() == f_lo.signum() {
                                    (lo_t, f_lo) = (mid, f);
                                } else {
                                    hi_t = mid;
                                }
                            }
                            breaks.push(0.5 * (lo_t + hi_t));
                        }
                    }
                }
                prev = Some((t, p));
            }
        }
        // How far the piece reaches across the chart, in quarter turns of
        // whichever parameters are angles.
        let mut pieces = match segment.curve {
            PlanarCurve::Line(_) => 1.0,
            _ => ((t1 - t0).abs() / QUARTER).ceil(),
        };
        if !matches!(self.surface, SurfaceGeometry::Plane(_)) {
            let (mut u, mut v) = (
                (f64::INFINITY, f64::NEG_INFINITY),
                (f64::INFINITY, f64::NEG_INFINITY),
            );
            for k in 0..=8 {
                let t = t0 + (t1 - t0) * f64::from(k) / 8.0;
                let (p, _) = segment.at(t, tol)?;
                u = (u.0.min(p.x), u.1.max(p.x));
                v = (v.0.min(p.y), v.1.max(p.y));
            }
            pieces = pieces
                .max(((u.1 - u.0) / QUARTER).ceil())
                .max(((v.1 - v.0) / QUARTER).ceil());
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let pieces = pieces.clamp(1.0, 64.0) as u32;
        for k in 1..pieces {
            breaks.push(t0 + (t1 - t0) * f64::from(k) / f64::from(pieces));
        }
        breaks.sort_by(f64::total_cmp);
        if t1 < t0 {
            breaks.reverse();
        }
        breaks.dedup();
        Ok(breaks)
    }

    /// How many Gauss points a boundary panel from `a` to `b` takes: the
    /// fewest whose error bound is below [`PANEL_MISS`] of the panel's own
    /// size, or ten.
    ///
    /// The integrand along a panel is a trigonometric polynomial of the
    /// chart point, of frequency at most [`FREQUENCY`] in an angle (a
    /// polynomial in a length, whose chart size stands in for a radian),
    /// times the pcurve's `v` speed. On a polynomial piece of degree three
    /// or less, written `c(s) = c0 + A1 s + A2 s^2 + A3 s^3` over the panel
    /// as `s` runs over `[-1, 1]`, the integrand is analytic in the
    /// Bernstein ellipse of every `rho`, bounded there through the bounds
    /// of `|Im c|` and `|c'|`, and the `n`-point rule misses by at most
    /// `64/15 M rho^(2 - 2n) / (rho^2 - 1)`. Any other piece, or a surface
    /// whose integrand is of another kind, takes ten.
    fn outer_order(
        &self,
        segment: &Segment,
        a: f64,
        b: f64,
        tol: Tolerances,
    ) -> OgeomResult<usize> {
        let cubic = match &segment.curve {
            PlanarCurve::Line(_) => true,
            PlanarCurve::BSpline(spline) => !spline.is_rational() && spline.knots().degree() <= 3,
            _ => false,
        };
        let length = 1.0 / self.scale.max(f64::MIN_POSITIVE);
        let (wu, wv) = match self.surface {
            SurfaceGeometry::Sphere(_) | SurfaceGeometry::Torus(_) => (1.0, 1.0),
            SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Cone(_) => (1.0, length),
            SurfaceGeometry::Plane(_) => (length, length),
            _ => return Ok(10),
        };
        if !cubic {
            return Ok(10);
        }
        // The cubic's coefficients from its points at s = -1, -1/3, 1/3, 1.
        let (middle, half) = (f64::midpoint(a, b), (b - a) * 0.5);
        let mut y = [Point2::ORIGIN; 4];
        for (slot, s) in y.iter_mut().zip([-1.0, -1.0 / 3.0, 1.0 / 3.0, 1.0]) {
            *slot = segment.at(middle + half * s, tol)?.0;
        }
        let odd = (y[3] - y[0]) * 0.5;
        let odd_inner = (y[2] - y[1]) * 0.5;
        let a3 = (odd - odd_inner * 3.0) * (9.0 / 8.0);
        let a1 = odd - a3;
        let a2 = ((y[0].to_vector() + y[3].to_vector()) - (y[1].to_vector() + y[2].to_vector()))
            * (9.0 / 16.0);
        let size = a1.y.abs() + a2.y.abs() + a3.y.abs();
        if size == 0.0 {
            return Ok(ORDERS[0]);
        }
        let reach = |c: Vector2| c.x.abs() * wu + c.y.abs() * wv;
        let (r1, r2, r3) = (reach(a1), reach(a2), reach(a3));
        let mut best = [f64::INFINITY; ORDERS.len()];
        for k in 0..48 {
            let rho = 1.05 * 1.25_f64.powi(k);
            let (big, small) = (0.5 * (rho + 1.0 / rho), 0.5 * (rho - 1.0 / rho));
            // |Im c| on the ellipse, and |dv/ds| there against the panel's
            // own size.
            let im = r1 * small
                + r2 * 2.0 * big * small
                + r3 * (3.0 * big * big + small * small) * small;
            let exponent = FREQUENCY * im;
            if exponent > 600.0 {
                break;
            }
            let speed = (a1.y.abs() + 2.0 * big * a2.y.abs() + 3.0 * big * big * a3.y.abs()) / size;
            let common = 64.0 / 15.0 * exponent.exp() * speed * rho * rho / (rho * rho - 1.0);
            for (slot, &n) in best.iter_mut().zip(&ORDERS) {
                #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
                let miss = common * rho.powi(-2 * n as i32);
                *slot = slot.min(miss);
            }
        }
        Ok(ORDERS
            .iter()
            .zip(best)
            .find(|&(_, miss)| miss <= PANEL_MISS)
            .map_or(10, |(&n, _)| n))
    }

    /// The inner integral from `u_ref` to the boundary point `at`, along
    /// `u` at its `v`, its samples weighted by `outer`, and whether `|n|`
    /// dipped across one of its panels where `folds` watches for it.
    fn inner(
        &self,
        at: Point2,
        outer: f64,
        fine: u32,
        folds: Folds,
        tol: Tolerances,
        sink: &mut dyn FnMut(Sample),
    ) -> OgeomResult<bool> {
        let at = into_domain(&self.surface, at);
        let (ua, ub) = (self.u_ref, at.x);
        if ua == ub || outer == 0.0 {
            return Ok(false);
        }
        // Along `u` at a fixed `v`, an analytic surface's point, `n dA` and
        // its length are trigonometric polynomials in `u` (polynomials on a
        // plane), which a quarter-turn panel takes to rounding: only the
        // boundary's own panels are refined between runs there, and none
        // is graded.
        let analytic = matches!(
            self.surface,
            SurfaceGeometry::Plane(_)
                | SurfaceGeometry::Cylinder(_)
                | SurfaceGeometry::Cone(_)
                | SurfaceGeometry::Sphere(_)
                | SurfaceGeometry::Torus(_)
        );
        let refined = if analytic { 1 } else { fine };
        let folds = if analytic { Folds::Ignore } else { folds };
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let pieces = if matches!(self.surface, SurfaceGeometry::Plane(_)) {
            1
        } else {
            ((ub - ua).abs() / QUARTER).ceil().clamp(1.0, 64.0) as u32
        } * refined;
        let u_knots = &self.knot_lines.0;
        let (lo, hi) = (ua.min(ub), ua.max(ub));
        let mut cuts: Vec<f64> = (0..=pieces)
            .map(|k| ua + (ub - ua) * f64::from(k) / f64::from(pieces))
            .collect();
        cuts.extend(u_knots.iter().copied().filter(|k| *k > lo && *k < hi));
        cuts.sort_by(f64::total_cmp);
        if ub < ua {
            cuts.reverse();
        }
        cuts.dedup();
        let isoline = Isoline::of(&self.surface, at.y, tol);
        let mut dipped = false;
        let mut normals = [Vector::ZERO; 10];
        for pair in cuts.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            if folds == Folds::Ignore {
                for (u, wu) in gauss_legendre_rule(a, b) {
                    let (p, du, dv) = match isoline.as_ref().and_then(|line| line.at(u, tol)) {
                        Some(found) => found,
                        None => self.surface.point_d1_at(u, at.y, tol)?,
                    };
                    sink((p, du.cross(dv) * self.sign, outer * wu));
                }
                continue;
            }
            let (mut least, mut most) = (f64::INFINITY, 0.0_f64);
            for ((u, wu), slot) in gauss_legendre_rule(a, b).into_iter().zip(&mut normals) {
                let (p, du, dv) = match isoline.as_ref().and_then(|line| line.at(u, tol)) {
                    Some(found) => found,
                    None => self.surface.point_d1_at(u, at.y, tol)?,
                };
                let n = du.cross(dv) * self.sign;
                let q = n.dot(n);
                (least, most) = (least.min(q), most.max(q));
                *slot = n;
                sink((p, n, outer * wu));
            }
            // Near its bottom a dip of `|n|` runs up both sides in
            // proportion to the distance, so a dip anywhere on the panel,
            // or just past its end, brings some node below half the
            // largest; a panel whose nodes all stand above that has none.
            if least < 0.25 * most {
                dipped = true;
                if folds == Folds::Grade {
                    self.regrade(isoline.as_ref(), at.y, (a, b), &normals, outer, tol, sink)?;
                }
            }
        }
        Ok(dipped)
    }

    /// An inner panel from `a` to `b` along `u` at `v`, whose samples
    /// (with `n` at its nodes `normals`) were summed, graded towards where
    /// `|n|` dips if it does (see [`fold_in`]): its samples are taken back,
    /// each summed again with its weight negated, and the graded panels'
    /// summed instead.
    #[cold]
    #[allow(clippy::too_many_arguments)]
    fn regrade(
        &self,
        isoline: Option<&Isoline<'_>>,
        v: f64,
        (a, b): (f64, f64),
        normals: &[Vector; 10],
        outer: f64,
        tol: Tolerances,
        sink: &mut dyn FnMut(Sample),
    ) -> OgeomResult<()> {
        let Some((x, levels)) = fold_in(normals) else {
            return Ok(());
        };
        let point = |u: f64| -> OgeomResult<(Point, Vector)> {
            let (p, du, dv) = match isoline.and_then(|line| line.at(u, tol)) {
                Some(found) => found,
                None => self.surface.point_d1_at(u, v, tol)?,
            };
            Ok((p, du.cross(dv) * self.sign))
        };
        for (u, wu) in gauss_legendre_rule(a, b) {
            let (p, n) = point(u)?;
            sink((p, n, -outer * wu));
        }
        let fold = f64::midpoint(a, b) + (b - a) * 0.5 * x;
        for (from, to) in graded(a, fold, b, levels) {
            for (u, wu) in gauss_legendre_rule(from, to) {
                let (p, n) = point(u)?;
                sink((p, n, outer * wu));
            }
        }
        Ok(())
    }
}

/// How small `|n|` may grow across an inner panel, against its largest
/// at the panel's nodes, before the panel is graded towards the fold.
const DIP: f64 = 0.25;

/// The ten-point rule's Legendre transform on `[-1, 1]`: row `k` carries
/// each node's share of the `k`-th Legendre coefficient of the values
/// there, exact for a polynomial of degree nine or less.
fn legendre_transform() -> &'static [[f64; 10]; 10] {
    static TABLE: std::sync::OnceLock<[[f64; 10]; 10]> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| {
        let mut table = [[0.0; 10]; 10];
        for (i, (x, w)) in gauss_legendre_rule(-1.0, 1.0).into_iter().enumerate() {
            let (mut before, mut p) = (0.0, 1.0);
            for (k, row) in table.iter_mut().enumerate() {
                #[allow(clippy::cast_precision_loss)]
                let k = k as f64;
                row[i] = (2.0 * k + 1.0) * 0.5 * w * p;
                let next = ((2.0 * k + 1.0) * x * p - k * before) / (k + 1.0);
                (before, p) = (p, next);
            }
        }
        table
    })
}

/// The Legendre series `c` at `x` in `[-1, 1]`.
fn legendre_at(c: &[Vector; 10], x: f64) -> Vector {
    let (mut before, mut p) = (0.0, 1.0);
    let mut sum = Vector::ZERO;
    for (k, ck) in c.iter().enumerate() {
        sum += *ck * p;
        #[allow(clippy::cast_precision_loss)]
        let k = k as f64;
        let next = ((2.0 * k + 1.0) * x * p - k * before) / (k + 1.0);
        (before, p) = (p, next);
    }
    sum
}

/// Where along an inner panel `|n|` nearly vanishes, from `n` at the
/// panel's ten Gauss nodes: the place in `[-1, 1]` and how many times the
/// panel is halved towards it to resolve the dip. `None` where `|n|` stays
/// above [`DIP`] of its largest at the nodes, or vanishes only at the
/// panel's end.
///
/// Along `u` at a fixed `v` a polynomial patch's `n` is a polynomial of
/// degree `2p - 1` within a knot span, which the nodes' Legendre series
/// gives exactly up to degree five in `u`. Near a fold `|n|` dips to a sliver
/// the nodes can step over, and `sqrt(|n|^2)` there converges slowly under
/// any uniform split; the series finds the dip, and panels halving towards
/// it down to its width take it to rounding.
fn fold_in(normals: &[Vector; 10]) -> Option<(f64, u32)> {
    let largest = normals.iter().map(|n| n.magnitude()).fold(0.0, f64::max);
    if largest == 0.0 {
        return None;
    }
    let table = legendre_transform();
    let mut c = [Vector::ZERO; 10];
    for (ck, row) in c.iter_mut().zip(table) {
        for (n, share) in normals.iter().zip(row) {
            *ck += *n * *share;
        }
    }
    // `|n|` stays at least the mean's length less every other term's,
    // each Legendre polynomial at most one in size.
    let rest: f64 = c[1..].iter().map(|ck| ck.magnitude()).sum();
    if c[0].magnitude() - rest > DIP * largest {
        return None;
    }
    let squared = |x: f64| {
        let n = legendre_at(&c, x);
        n.dot(n)
    };
    const SAMPLES: u32 = 32;
    let step = 2.0 / f64::from(SAMPLES);
    let (mut best, mut low) = (-1.0, squared(-1.0));
    for k in 1..=SAMPLES {
        let x = -1.0 + step * f64::from(k);
        let q = squared(x);
        if q < low {
            (best, low) = (x, q);
        }
    }
    let (mut a, mut b) = ((best - step).max(-1.0), (best + step).min(1.0));
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    for _ in 0..50 {
        let (c1, c2) = (b - (b - a) * ratio, a + (b - a) * ratio);
        if squared(c1) < squared(c2) {
            b = c2;
        } else {
            a = c1;
        }
    }
    let x = f64::midpoint(a, b);
    let low = squared(x);
    if low.sqrt() > DIP * largest {
        return None;
    }
    // The dip's width from `|n|^2 = m^2 + q (x - x*)^2` about its bottom.
    let h = 1e-3;
    let bend: f64 = ((squared(x + h) + squared(x - h) - 2.0 * low) / (h * h) * 0.5).max(0.0);
    let width: f64 = if bend > 0.0 { (low / bend).sqrt() } else { 0.0 };
    // Each side of the bottom, `|n|` runs as `sqrt(t^2 + w^2)` in the
    // distance `t`, which the rule takes to within about `w^2` of the
    // panel's own integral: a dip narrower than [`NARROW`] needs only the
    // cut at its bottom, and none at the panel's end (a pole on its edge).
    if width < NARROW {
        return (x.abs() < 1.0 - 1e-6).then_some((x, 0));
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let levels = ((2.0 / width).log2().ceil() + 2.0) as u32;
    Some((x, levels))
}

/// How narrow a dip of `|n|`, against the panel's half-width, is taken
/// without grading towards it.
const NARROW: f64 = 1e-6;

/// The panels from `a` to `b` halving towards `fold` between them, `levels`
/// on each side, in order.
fn graded(a: f64, fold: f64, b: f64, levels: u32) -> Vec<(f64, f64)> {
    let mut panels = Vec::with_capacity(2 * levels as usize + 2);
    if fold != a {
        let mut from = a;
        for k in 1..=levels {
            let to = fold + (a - fold) * 0.5_f64.powi(k.cast_signed());
            panels.push((from, to));
            from = to;
        }
        panels.push((from, fold));
    }
    if fold != b {
        let mut from = fold;
        for k in (1..=levels).rev() {
            let to = fold + (b - fold) * 0.5_f64.powi(k.cast_signed());
            panels.push((from, to));
            from = to;
        }
        panels.push((from, b));
    }
    panels
}

/// A polynomial spline patch read along one `v`: its control net summed
/// across `v` once, into the control points of the row at `v` and of its
/// `v` derivative. Every sample of an inner integral shares that `v`, so a
/// sample costs one row's worth of basis functions instead of the grid's.
struct Isoline<'a> {
    knots: &'a ogeom_math::KnotVector,
    row: Vec<Point>,
    across: Vec<Vector>,
}

impl<'a> Isoline<'a> {
    fn of(surface: &'a SurfaceGeometry, v: f64, tol: Tolerances) -> Option<Self> {
        let SurfaceGeometry::BSpline(patch) = surface else {
            return None;
        };
        if patch.is_rational() {
            return None;
        }
        let v_knots = patch.v_knots();
        let span = v_knots.span(v, tol).ok()?;
        let q = v_knots.degree();
        let basis = v_knots.basis_derivatives(span, v, 1);
        let grid = patch.grid();
        let mut row = Vec::with_capacity(grid.u_count());
        let mut across = Vec::with_capacity(grid.u_count());
        for i in 0..grid.u_count() {
            let (mut p, mut d) = (Vector::ZERO, Vector::ZERO);
            for k in 0..=q {
                let c = grid.get(i, span - q + k)?.scaled.to_vector();
                p += c * basis[0][k];
                d += c * basis[1][k];
            }
            row.push(Point::ORIGIN + p);
            across.push(d);
        }
        Some(Self {
            knots: patch.u_knots(),
            row,
            across,
        })
    }

    /// The point, `du` and `dv` at `u`; `None` off the knots' domain.
    #[inline(always)]
    fn at(&self, u: f64, tol: Tolerances) -> Option<(Point, Vector, Vector)> {
        let span = self.knots.span(u, tol).ok()?;
        let p = self.knots.degree();
        let basis = self.knots.basis_derivatives(span, u, 1);
        let (mut point, mut du, mut dv) = (Vector::ZERO, Vector::ZERO, Vector::ZERO);
        for k in 0..=p {
            let c = self.row[span - p + k].to_vector();
            point += c * basis[0][k];
            du += c * basis[1][k];
            dv += self.across[span - p + k] * basis[0][k];
        }
        Some((Point::ORIGIN + point, du, dv))
    }
}

/// Whether two runs agree, against the size of what they integrate.
fn settled(a: Proxy, b: Proxy) -> bool {
    let size: f64 = b.iter().map(|x| x.abs()).sum();
    let miss = a
        .iter()
        .zip(&b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max);
    miss <= AGREE * size || size == 0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Each short rule integrates every polynomial of degree below twice
    /// its points exactly, on a panel away from the origin.
    #[test]
    fn the_short_rules_are_exact_to_their_degree() {
        let (a, b) = (0.75, 2.0);
        for n in ORDERS {
            let points = rule(n, a, b);
            assert_eq!(points.len(), n);
            for degree in 0..2 * n {
                #[allow(clippy::cast_possible_truncation, clippy::cast_possible_wrap)]
                let k = degree as i32;
                let exact = (b.powi(k + 1) - a.powi(k + 1)) / f64::from(k + 1);
                let sum: f64 = points.iter().map(|&(x, w)| w * x.powi(k)).sum();
                assert!(
                    (sum - exact).abs() <= 1e-14 * exact.abs(),
                    "{n} points, degree {degree}: {sum} against {exact}"
                );
            }
        }
    }

    /// A normal that nearly vanishes between two of a panel's nodes, as
    /// across a fold: `n = (x - b, e, 0)`, so `|n| = sqrt((x - b)^2 + e^2)`
    /// with a closed-form integral. The dip is found, and the panels graded
    /// towards it integrate `|n|` to rounding where the plain rule misses
    /// by orders of magnitude more. A normal vanishing at the panel's end
    /// (a pole on its edge) is left to the plain rule.
    #[test]
    fn a_dip_of_the_normal_between_nodes_is_graded_to_rounding() {
        let (bottom, e) = (0.3, 1e-4);
        let normal = |x: f64| Vector::new(x - bottom, e, 0.0);
        let primitive = |t: f64| 0.5 * (t * t.hypot(e) + e * e * (t / e).asinh());
        let exact = primitive(1.0 - bottom) - primitive(-1.0 - bottom);
        let rule = gauss_legendre_rule(-1.0, 1.0);
        let normals = rule.map(|(x, _)| normal(x));
        let plain: f64 = rule.iter().map(|&(x, w)| normal(x).magnitude() * w).sum();
        let Some((x, levels)) = fold_in(&normals) else {
            panic!("the dip is not found");
        };
        // The bottom is found on the nodes' Legendre series, which carries
        // this linear normal exactly.
        assert!((x - bottom).abs() < 1e-9, "dip found at {x}");
        let graded: f64 = graded(-1.0, x, 1.0, levels)
            .into_iter()
            .flat_map(|(a, b)| gauss_legendre_rule(a, b))
            .map(|(u, w)| normal(u).magnitude() * w)
            .sum();
        // Rounding over a few hundred samples of an integral of about one.
        assert!(
            (graded - exact).abs() < 1e-13,
            "graded {graded} against {exact}"
        );
        assert!(
            (plain - exact).abs() > 1e-6,
            "plain {plain} against {exact}"
        );

        let pole = rule.map(|(x, _)| Vector::new(x + 1.0, 0.0, 0.0));
        assert!(fold_in(&pole).is_none(), "a pole on the panel's end");
    }
}
