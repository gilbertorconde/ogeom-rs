//! Small helpers shared by the closed forms and the general path.

use super::*;

// --- small helpers -----------------------------------------------------------

pub(super) fn clamp_2d(curve: &PlanarCurve, t: f64) -> f64 {
    let (lo, hi) = curve.domain();
    if curve.is_periodic() {
        let span = hi - lo;
        if span > 0.0 {
            return lo + (t - lo).rem_euclid(span);
        }
    }
    t.clamp(lo, hi)
}

pub(super) fn clamp_3d(curve: &Curve, t: f64) -> f64 {
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
pub(super) fn segments_cross_2d(a: (Point2, Point2), b: (Point2, Point2)) -> Option<(f64, f64)> {
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
pub(super) fn distance_to_curve_3d(
    b: &Curve,
    sb: &Sampled<Point>,
    segments: &Segments,
    p: Point,
    tol: Tolerances,
) -> f64 {
    let mut best = (0_usize, f64::INFINITY);
    for i in 1..sb.points.len() {
        if segments.beyond(i, p, best.1) {
            continue;
        }
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

/// How far from `on_a`, along the first curve, a point stays within `gap`
/// of the second by `distance`: in each direction, the last parameter found
/// inside and the first found outside, if the walk left the valley before
/// the curve's end.
pub(super) fn valley_ends_3d(
    a: &Curve,
    on_a: f64,
    gap: f64,
    tol: Tolerances,
    distance: &dyn Fn(Point) -> f64,
) -> [(f64, Option<f64>); 2] {
    let (lo, hi) = a.domain();
    let span = hi - lo;
    let mut ends = [(on_a, None); 2];
    if span <= 0.0 {
        return ends;
    }
    let inside = |t: f64| -> bool {
        if t < lo || t > hi {
            return false;
        }
        a.point_at(t, tol).is_ok_and(|q| distance(q) <= gap)
    };
    for (end, direction) in ends.iter_mut().zip([-1.0, 1.0]) {
        let mut step = span * 1e-6;
        let mut last_in = on_a;
        let mut first_out: Option<f64> = None;
        while step <= span {
            let t = on_a + direction * step;
            if inside(t) {
                last_in = t;
                step *= 2.0;
            } else {
                first_out = Some(t);
                break;
            }
        }
        *end = match first_out {
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
                (r#in, Some(out))
            }
            None => (last_in, None),
        };
    }
    ends
}

/// How far a crossing's valley reaches in space: the larger of the distances
/// from the crossing to the valley's two ends on the first curve.
pub(super) fn valley_span_3d(
    a: &Curve,
    crossing: &Crossing<Point>,
    ends: &[(f64, Option<f64>); 2],
    tol: Tolerances,
) -> f64 {
    let (lo, hi) = a.domain();
    if hi - lo <= 0.0 {
        return 0.0;
    }
    ends.iter()
        .filter_map(|(edge, _)| a.point_at(*edge, tol).ok())
        .map(|q| q.distance(crossing.point))
        .fold(0.0, f64::max)
}

/// The closest approach of two spatial segments, as fractions and a distance.
pub(super) fn segments_approach_3d(a: (Point, Point), b: (Point, Point)) -> (f64, f64, f64) {
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

pub(super) fn order(a: f64, b: f64) -> (f64, f64) {
    if a <= b { (a, b) } else { (b, a) }
}

pub(super) fn sort_crossings<P>(crossings: &mut [Crossing<P>]) {
    crossings.sort_by(|x, y| {
        x.on_a
            .partial_cmp(&y.on_a)
            .unwrap_or(core::cmp::Ordering::Equal)
    });
}

pub(super) fn push_unique_2d(
    crossings: &mut Vec<Crossing<Point2>>,
    found: Crossing<Point2>,
    tol: Tolerances,
) {
    let reach = tol.confusion() * 100.0;
    if crossings
        .iter()
        .any(|c| c.point.distance(found.point) <= reach)
    {
        return;
    }
    crossings.push(found);
}

/// How close two polished crossings stand, with the first curve staying
/// there between them, to be one root found from two seeds.
pub(super) fn one_crossing_radius(tol: Tolerances) -> f64 {
    tol.confusion() * 100.0
}

/// Keep `found` unless it is a crossing already held: one standing at the
/// same point with the first curve staying there between the two, which is
/// one root polished from two seeds. A curve passing the same point twice
/// (a section shaped as a figure eight, crossing an edge at its double
/// point) crosses there twice, once per passage, and both are kept. A
/// closed first curve is measured the short way round its seam.
pub(super) fn push_unique_3d(
    a: &Curve,
    crossings: &mut Vec<Crossing<Point>>,
    found: Crossing<Point>,
    tol: Tolerances,
) {
    let reach = one_crossing_radius(tol);
    let (lo, hi) = a.domain();
    let period = hi - lo;
    let closed = a.is_periodic()
        || (period.is_finite()
            && a.point_at(lo, tol)
                .ok()
                .zip(a.point_at(hi, tol).ok())
                .is_some_and(|(p, q)| p.distance(q) <= reach));
    let same = |c: &Crossing<Point>| {
        if c.point.distance(found.point) > reach {
            return false;
        }
        let mut mid = f64::midpoint(c.on_a, found.on_a);
        if closed && (c.on_a - found.on_a).abs() > period * 0.5 {
            mid += period * 0.5;
            if mid > hi {
                mid -= period;
            }
        }
        a.point_at(mid, tol)
            .is_ok_and(|p| p.distance(found.point) <= reach)
    };
    if crossings.iter().any(same) {
        return;
    }
    crossings.push(found);
}
