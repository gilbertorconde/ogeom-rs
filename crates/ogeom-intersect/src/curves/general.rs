//! The general path: curves sampled into segments, seeded, and polished onto crossings.

use super::*;

// --- the general path --------------------------------------------------------

/// Sampled segments of one curve, with the parameters they span.
pub(super) struct Sampled<P> {
    pub(super) points: Vec<P>,
    pub(super) parameters: Vec<f64>,
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

pub(super) fn sample_3d(curve: &Curve, n: usize, tol: Tolerances) -> Sampled<Point> {
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

pub(super) fn general_2d(
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

pub(super) fn general_3d(
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
    let segments_b = Segments::of(&sb);
    let feet: Vec<Option<(f64, f64)>> = sa
        .points
        .iter()
        .map(|p| foot_via_samples(b, &sb, &segments_b, *p, tol))
        .collect();
    let overlaps = if hugging_runs(&feet, options.gap).contains(&true) {
        shared_support_3d(a, b, &sa, (&sb, &segments_b), &feet, options, tol)
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
    // finds the same few crossings over and over.
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
    let (boxes_a, boxes_b) = (boxes(&sa), &segments_b.boxes);
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
                push_unique_3d(a, &mut crossings, found, tol);
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
        let ends = valley_ends_3d(a, c.on_a, gap, tol, &|q| {
            distance_to_curve_3d(b, &sb, &segments_b, q, tol)
        });
        let valley = valley_span_3d(a, c, &ends, tol);
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

/// A sampled curve's segments boxed, for the scans for a point's nearest
/// segment: a segment whose box stands further from the point than the
/// nearest found so far, past the rounding of either distance, is not
/// nearer, and is passed over without changing what the scan finds.
pub(super) struct Segments {
    boxes: Vec<(Point, Point)>,
    /// The largest coordinate among the samples.
    size: f64,
    /// The segment last found nearest, which the next scan starts from.
    last: core::cell::Cell<usize>,
}

impl Segments {
    pub(super) fn of(s: &Sampled<Point>) -> Self {
        let boxes = s
            .points
            .windows(2)
            .map(|w| {
                (
                    Point::new(w[0].x.min(w[1].x), w[0].y.min(w[1].y), w[0].z.min(w[1].z)),
                    Point::new(w[0].x.max(w[1].x), w[0].y.max(w[1].y), w[0].z.max(w[1].z)),
                )
            })
            .collect();
        let size = s.points.iter().fold(0.0_f64, |m, p| {
            m.max(p.x.abs()).max(p.y.abs()).max(p.z.abs())
        });
        Self {
            boxes,
            size,
            last: core::cell::Cell::new(0),
        }
    }

    /// Whether the box of the segment ending at sample `j` stands further
    /// from `p` than `nearest` by more than rounding. A segment's distance
    /// is formed from its ends and `p` in a few dozen rounded operations,
    /// each off by at most half an ulp of the largest coordinate, and the
    /// box's distance in a few more; no segment so boxed measures nearer.
    pub(super) fn beyond(&self, j: usize, p: Point, nearest: f64) -> bool {
        const ROUNDING: f64 = 64.0 * f64::EPSILON;
        let (lo, hi) = self.boxes[j - 1];
        let out = |x: f64, l: f64, h: f64| (l - x).max(x - h).max(0.0);
        let off = ogeom_math::Vector::new(
            out(p.x, lo.x, hi.x),
            out(p.y, lo.y, hi.y),
            out(p.z, lo.z, hi.z),
        );
        let size = self.size.max(p.x.abs()).max(p.y.abs()).max(p.z.abs());
        off.magnitude() > ROUNDING.mul_add(size, nearest)
    }
}

/// The foot of a point on the second curve, seeded from its sampled
/// polyline's nearest segment: the parameter and the distance there.
pub(super) fn foot_via_samples(
    b: &Curve,
    sb: &Sampled<Point>,
    segments: &Segments,
    p: Point,
    tol: Tolerances,
) -> Option<(f64, f64)> {
    let approach = |j: usize| segments_approach_3d((p, p), (sb.points[j - 1], sb.points[j]));
    // The segment nearest the last point asked about bounds the nearest
    // to this one, which the scan then finds in the same order as ever.
    let hint = segments.last.get();
    let mut bound = if hint >= 1 && hint < sb.points.len() {
        approach(hint).2
    } else {
        f64::INFINITY
    };
    let mut seed = (f64::INFINITY, 0.0);
    for j in 1..sb.points.len() {
        if segments.beyond(j, p, bound) {
            continue;
        }
        let (_, tb, gap) = approach(j);
        if gap < seed.0 {
            seed = (
                gap,
                sb.parameters[j - 1] + (sb.parameters[j] - sb.parameters[j - 1]) * tb,
            );
            bound = bound.min(gap);
            segments.last.set(j);
        }
    }
    if cfg!(debug_assertions) {
        let mut every = (f64::INFINITY, 0.0);
        for j in 1..sb.points.len() {
            let (_, tb, gap) = segments_approach_3d((p, p), (sb.points[j - 1], sb.points[j]));
            if gap < every.0 {
                every = (
                    gap,
                    sb.parameters[j - 1] + (sb.parameters[j] - sb.parameters[j - 1]) * tb,
                );
            }
        }
        assert!(
            every.0.to_bits() == seed.0.to_bits() && every.1.to_bits() == seed.1.to_bits(),
            "boxed scan found {seed:?}, every segment {every:?}"
        );
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

/// Runs of the first curve's samples whose feet on the second lie within
/// the gap, each bisected to its parametric ends.
fn shared_support_3d(
    a: &Curve,
    b: &Curve,
    sa: &Sampled<Point>,
    (sb, segments_b): (&Sampled<Point>, &Segments),
    feet: &[Option<(f64, f64)>],
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Vec<Overlap> {
    let hugs = |t: f64| -> Option<(f64, f64)> {
        let p = a.point_at(t, tol).ok()?;
        let (s, gap) = foot_via_samples(b, sb, segments_b, p, tol)?;
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
        let (hi_a, mut hi_b) = refine(end, (end + 1 < sa.points.len()).then_some(end + 1));
        // On a periodic second curve the run may pass its seam: its feet
        // are followed from the start, each taken a whole turn nearest the
        // last, and the far end read beside the last of them, so the run
        // keeps the way round it goes. A trim of one holds its feet in its
        // window, which no run passes the seam of.
        if b.is_periodic() && !matches!(b, Curve::Trimmed(_)) {
            let (from, to) = b.domain();
            let period = to - from;
            let mut last = lo_b;
            for foot in feet[start..=end].iter().flatten() {
                last = foot.0 + period * ((last - foot.0) / period).round();
            }
            hi_b += period * ((last - hi_b) / period).round();
        }
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
pub(super) fn foot_on_3d(
    curve: &Curve,
    p: Point,
    seed: f64,
    tol: Tolerances,
) -> Option<(f64, f64)> {
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

/// Gauss-Newton on the closest approach of two space curves.
///
/// Three equations would be overdetermined for two unknowns, so the system is
/// the two *stationarity* conditions (the gap vector perpendicular to both
/// tangents), whose solutions are the local closest approaches. The gap test
/// afterwards decides whether the approach found is a crossing.
///
/// Newton damped on the stationarity residual can settle where that residual
/// is least without being zero: a fitted curve that all but stops along its
/// parameter, as a fit through a sharp turn can, holds every seed nearby at
/// such a point, short of a crossing beside it. Where the solve ends there
/// without converging, the approach is walked on by alternating feet, each
/// curve's point dropped onto the other, which never lets the gap grow, and
/// polished again from where the walk stops.
fn polish_3d(
    a: &Curve,
    b: &Curve,
    seed_a: f64,
    seed_b: f64,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Option<Crossing<Point>> {
    let (t, s, converged) = stationary_3d(a, b, seed_a, seed_b, tol)?;
    let mut found = crossing_at_3d(a, b, t, s, tol)?;
    if found.gap > options.gap && !converged {
        let (t, s) = alternating_feet_3d(a, b, t, s, found.gap, tol)?;
        let walked = crossing_at_3d(a, b, t, s, tol)?;
        let polished = stationary_3d(a, b, t, s, tol)
            .and_then(|(t, s, _)| crossing_at_3d(a, b, t, s, tol))
            .filter(|p| p.gap < walked.gap);
        found = polished.unwrap_or(walked);
    }
    (found.gap <= options.gap).then_some(found)
}

/// The two curves' points at `t` and `s` as a crossing, with their gap.
fn crossing_at_3d(
    a: &Curve,
    b: &Curve,
    t: f64,
    s: f64,
    tol: Tolerances,
) -> Option<Crossing<Point>> {
    let pa = a.point_at(t, tol).ok()?;
    let pb = b.point_at(s, tol).ok()?;
    Some(Crossing {
        on_a: t,
        on_b: s,
        point: pa,
        gap: pa.distance(pb),
        reach: 0.0,
    })
}

/// From `t` on the first curve and `s` on the second, `gap` apart, each
/// point's foot on the other curve in turn while the gap falls by more than
/// a hundredth of the confusion distance a round, at most sixty-four rounds.
/// A foot that would widen the gap is not taken, so the pair returned is
/// never further apart than the pair given.
fn alternating_feet_3d(
    a: &Curve,
    b: &Curve,
    mut t: f64,
    mut s: f64,
    mut gap: f64,
    tol: Tolerances,
) -> Option<(f64, f64)> {
    for _ in 0..64 {
        let before = gap;
        let pa = a.point_at(t, tol).ok()?;
        if let Some((next, d)) = foot_on_3d(b, pa, s, tol)
            && d < gap
        {
            (s, gap) = (next, d);
        }
        let pb = b.point_at(s, tol).ok()?;
        if let Some((next, d)) = foot_on_3d(a, pb, t, tol)
            && d < gap
        {
            (t, gap) = (next, d);
        }
        if before - gap <= tol.confusion() * 0.01 {
            break;
        }
    }
    Some((t, s))
}

/// Newton on the stationarity system from a seed: the parameters it ends
/// at, clamped to each curve's domain, and whether the residual met its
/// tolerance there.
fn stationary_3d(
    a: &Curve,
    b: &Curve,
    seed_a: f64,
    seed_b: f64,
    tol: Tolerances,
) -> Option<(f64, f64, bool)> {
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
    let ([t, s], _, convergence, _) =
        solve::newton_system_2(system, [seed_a, seed_b], criteria).ok()?;
    Some((
        clamp_3d(a, t),
        clamp_3d(b, s),
        convergence == solve::Convergence::Residual,
    ))
}
