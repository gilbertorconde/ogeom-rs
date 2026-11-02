//! Splitting a polyline into runs that a line or a circle each follow
//! within a tolerance: the segmentation the outline fit stands on.
//!
//! From the first point, each run is the longest stretch that one line or
//! one circle follows within the tolerance, found by doubling the stretch
//! until it no longer fits and then halving the step back; a circle is
//! taken only where it reaches farther than a line, and only where the
//! points bow from their chord by more than the tolerance. The next run
//! starts at the point the last one ended on, so consecutive runs share a
//! point, which is then moved to where the two fits leave the least
//! squared residual between them.
//!
//! A stretch fits when every point of the polyline along it, between the
//! points as well as at them, lies within the tolerance of the shape: for a
//! line the distance along a segment is largest at an end, and for a circle
//! a chord can also dip inside it, which is measured exactly.

use ogeom_math::{Point2, Vector2};

/// What follows a run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Shape {
    /// The line through `at` along the unit `along`.
    Line { at: Point2, along: Vector2 },
    /// The circle about `centre`.
    Circle { centre: Point2, radius: f64 },
}

/// A stretch of the polyline, `first..=last` by index, and its shape.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Run {
    pub(crate) first: usize,
    pub(crate) last: usize,
    pub(crate) shape: Shape,
}

/// Split an open polyline into runs, each followed by a line or a circle
/// within `tolerance`. The first run starts at the first point and the
/// last ends at the last; a closed loop is passed with its first point
/// repeated at the end.
pub(crate) fn runs(points: &[Point2], tolerance: f64) -> Vec<Run> {
    let mut out = Vec::new();
    if points.len() < 2 {
        return out;
    }
    let end = points.len() - 1;
    let mut first = 0;
    while first < end {
        let (line_to, line) = reach(first, end, 1, |last| {
            line_through(&points[first..=last])
                .filter(|s| deviation(&points[first..=last], *s) <= tolerance)
        });
        let circle = (first + 2 <= end)
            .then(|| {
                reach(first, end, 2, |last| {
                    let run = &points[first..=last];
                    circle_through(run)
                        .filter(|s| bow(run) > tolerance && deviation(run, *s) <= tolerance)
                })
            })
            .filter(|(to, shape)| *to > line_to && shape.is_some());
        let (last, shape) = match (circle, line) {
            (Some((to, Some(shape))), _) => (to, shape),
            (_, Some(shape)) => (line_to, shape),
            // Two points always have a line; reaching none means they
            // are one point, which a polyline of distinct points lacks.
            _ => (first + 1, chord(points[first], points[first + 1])),
        };
        out.push(Run { first, last, shape });
        first = last;
    }
    for k in 1..out.len() {
        balance(points, &mut out, k, tolerance);
    }
    out
}

/// The shape of the same kind fitted to `points`, if it follows them
/// within `tolerance`.
fn refit(kind: Shape, points: &[Point2], tolerance: f64) -> Option<Shape> {
    let shape = match kind {
        Shape::Line { .. } => line_through(points),
        Shape::Circle { .. } => circle_through(points),
    }?;
    (deviation(points, shape) <= tolerance).then_some(shape)
}

/// The sum of the points' squared distances from the shape.
fn squares(points: &[Point2], shape: Shape) -> f64 {
    match shape {
        Shape::Line { at, along } => points.iter().map(|p| along.cross(*p - at).powi(2)).sum(),
        Shape::Circle { centre, radius } => points
            .iter()
            .map(|p| ((*p - centre).magnitude() - radius).powi(2))
            .sum(),
    }
}

/// Move the point runs `k - 1` and `k` share to where the two fits leave
/// the least squared residual between them, each still fitting. Greedy
/// runs reach as far as they can, so a line eats into the arc after it;
/// the least residual is where one shape gives way to the other, at a
/// corner the corner itself.
fn balance(points: &[Point2], runs: &mut [Run], k: usize, tolerance: f64) {
    let (a, b) = (runs[k - 1], runs[k]);
    let least = |r: &Run| {
        if matches!(r.shape, Shape::Line { .. }) {
            1
        } else {
            2
        }
    };
    let first_fits = |m: usize| refit(a.shape, &points[a.first..=m], tolerance);
    let second_fits = |m: usize| refit(b.shape, &points[m..=b.last], tolerance);
    // The first fits up to some point and the second from some point:
    // find both by halving, the current shared point fitting both.
    let (mut lo, mut hi) = (a.first + least(&a), b.last);
    let mut from = b.first;
    while lo < hi {
        let m = lo + (hi - lo) / 2;
        if m <= b.last - least(&b) && second_fits(m).is_some() {
            hi = m;
        } else {
            lo = m + 1;
        }
    }
    if lo <= b.first {
        from = lo;
    }
    let (mut lo, mut hi) = (a.last, b.last - least(&b));
    let mut to = a.last;
    while lo < hi {
        let m = hi - (hi - lo) / 2;
        if m >= a.first + least(&a) && first_fits(m).is_some() {
            lo = m;
        } else {
            hi = m - 1;
        }
    }
    if lo >= a.last {
        to = lo;
    }
    let cost = |m: usize| -> Option<(f64, Shape, Shape)> {
        let (sa, sb) = (first_fits(m)?, second_fits(m)?);
        Some((
            squares(&points[a.first..=m], sa) + squares(&points[m..=b.last], sb),
            sa,
            sb,
        ))
    };
    // A ternary search for the least residual over the shared points
    // both fit from.
    let (mut l, mut h) = (from, to);
    while h - l > 2 {
        let (m1, m2) = (l + (h - l) / 3, h - (h - l) / 3);
        match (cost(m1), cost(m2)) {
            (Some(c1), Some(c2)) if c1.0 <= c2.0 => h = m2,
            (Some(_), Some(_)) => l = m1,
            _ => break,
        }
    }
    let best = (l..=h)
        .filter_map(|m| cost(m).map(|c| (m, c)))
        .min_by(|x, y| x.1.0.total_cmp(&y.1.0));
    if let Some((m, (_, sa, sb))) = best {
        runs[k - 1] = Run {
            last: m,
            shape: sa,
            ..a
        };
        runs[k] = Run {
            first: m,
            shape: sb,
            ..b
        };
    }
}

/// How far the points stand from the chord between the first and the
/// last: a circle bowing no more than the tolerance is taken for a line.
fn bow(points: &[Point2]) -> f64 {
    let (Some(&a), Some(&b)) = (points.first(), points.last()) else {
        return 0.0;
    };
    let Shape::Line { at, along } = chord(a, b) else {
        return 0.0;
    };
    if a == b {
        return f64::INFINITY;
    }
    deviation(points, Shape::Line { at, along })
}

/// The farthest `last` from `first + least` up to `end` for which `fit`
/// gives a shape, by doubling then halving the step, and that shape.
fn reach(
    first: usize,
    end: usize,
    least: usize,
    fit: impl Fn(usize) -> Option<Shape>,
) -> (usize, Option<Shape>) {
    let mut good = (first + least - 1, None);
    let mut step = least;
    let bad;
    loop {
        let at = (first + step).min(end);
        match fit(at) {
            Some(s) => {
                good = (at, Some(s));
                if at == end {
                    return good;
                }
                step *= 2;
            }
            None => {
                bad = at;
                break;
            }
        }
    }
    let (mut lo, mut hi) = (good.0, bad);
    while hi - lo > 1 {
        let mid = lo + (hi - lo) / 2;
        match fit(mid) {
            Some(s) => {
                lo = mid;
                good = (mid, Some(s));
            }
            None => hi = mid,
        }
    }
    good
}

fn chord(a: Point2, b: Point2) -> Shape {
    let d = b - a;
    let len = d.magnitude();
    Shape::Line {
        at: a,
        along: if len > 0.0 {
            d / len
        } else {
            Vector2::new(1.0, 0.0)
        },
    }
}

/// The least-squares line through points, along the direction from the
/// first to the last.
pub(crate) fn line_through(points: &[Point2]) -> Option<Shape> {
    let n = count(points.len());
    let (first, last) = (*points.first()?, *points.last()?);
    if points.len() == 2 {
        return (first != last).then(|| chord(first, last));
    }
    let mean = points
        .iter()
        .fold(Vector2::new(0.0, 0.0), |s, p| s + p.to_vector())
        / n;
    let (mut xx, mut xy, mut yy) = (0.0, 0.0, 0.0);
    for p in points {
        let d = p.to_vector() - mean;
        xx += d.x * d.x;
        xy += d.x * d.y;
        yy += d.y * d.y;
    }
    // The major axis of the scatter: angle half that of (xx - yy, 2xy).
    let angle = 0.5 * (2.0 * xy).atan2(xx - yy);
    let mut along = Vector2::new(angle.cos(), angle.sin());
    let span = last - first;
    if along.dot(span) < 0.0 {
        along = -along;
    }
    Some(Shape::Line {
        at: Point2::from_vector(mean),
        along,
    })
}

/// The circle nearest the points: the algebraic fit, then a few steps of
/// Gauss-Newton on the distances themselves. `None` where the points are
/// too near a line for a circle to be told from it.
pub(crate) fn circle_through(points: &[Point2]) -> Option<Shape> {
    if points.len() < 3 {
        return None;
    }
    let n = count(points.len());
    let mean = points
        .iter()
        .fold(Vector2::new(0.0, 0.0), |s, p| s + p.to_vector())
        / n;
    let scale = points
        .iter()
        .map(|p| (p.to_vector() - mean).magnitude())
        .fold(0.0, f64::max);
    if scale <= 0.0 {
        return None;
    }
    // x^2 + y^2 + d x + e y + f = 0 in coordinates about the mean, scaled.
    let mut m = [[0.0f64; 3]; 3];
    let mut r = [0.0f64; 3];
    for p in points {
        let v = (p.to_vector() - mean) / scale;
        let row = [v.x, v.y, 1.0];
        let rhs = -(v.x * v.x + v.y * v.y);
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += row[i] * row[j];
            }
            r[i] += row[i] * rhs;
        }
    }
    let [d, e, f] = solve3(m, r)?;
    let (cx, cy) = (-d / 2.0, -e / 2.0);
    let rr = cx.mul_add(cx, cy * cy) - f;
    if rr.is_nan() || rr <= 0.0 {
        return None;
    }
    let mut centre = Vector2::new(cx, cy) * scale + mean;
    let mut radius = rr.sqrt() * scale;
    // A radius thousands of times the points' spread is a line.
    if radius.is_nan() || radius >= 1e3 * scale {
        return None;
    }
    for _ in 0..4 {
        let (mut jtj, mut jtr) = ([[0.0f64; 3]; 3], [0.0f64; 3]);
        for p in points {
            let d = p.to_vector() - centre;
            let dist = d.magnitude();
            if dist == 0.0 {
                return None;
            }
            let row = [-d.x / dist, -d.y / dist, -1.0];
            let res = dist - radius;
            for i in 0..3 {
                for j in 0..3 {
                    jtj[i][j] += row[i] * row[j];
                }
                jtr[i] -= row[i] * res;
            }
        }
        let Some(step) = solve3(jtj, jtr) else { break };
        centre += Vector2::new(step[0], step[1]);
        radius += step[2];
        if radius.is_nan() || radius <= 0.0 {
            return None;
        }
    }
    Some(Shape::Circle {
        centre: Point2::from_vector(centre),
        radius,
    })
}

/// A count of points as a double.
#[allow(
    clippy::cast_precision_loss,
    reason = "point counts are far below 2^52"
)]
const fn count(n: usize) -> f64 {
    n as f64
}

fn solve3(m: [[f64; 3]; 3], r: [f64; 3]) -> Option<[f64; 3]> {
    let det = |a: [[f64; 3]; 3]| {
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    };
    let d = det(m);
    let size = m.iter().flatten().map(|v| v.abs()).fold(0.0, f64::max);
    if d.is_nan() || d.abs() <= 1e-12 * size * size * size {
        return None;
    }
    let mut out = [0.0; 3];
    for (k, o) in out.iter_mut().enumerate() {
        let mut a = m;
        for i in 0..3 {
            a[i][k] = r[i];
        }
        *o = det(a) / d;
    }
    out.iter().all(|v| v.is_finite()).then_some(out)
}

/// The largest distance from the shape of any point of the polyline.
pub(crate) fn deviation(points: &[Point2], shape: Shape) -> f64 {
    match shape {
        Shape::Line { at, along } => points
            .iter()
            .map(|p| along.cross(*p - at).abs())
            .fold(0.0, f64::max),
        Shape::Circle { centre, radius } => {
            let mut worst = 0.0f64;
            for (k, p) in points.iter().enumerate() {
                worst = worst.max(((*p - centre).magnitude() - radius).abs());
                if let Some(q) = points.get(k + 1) {
                    worst = worst.max(radius - nearest_on_segment(centre, *p, *q));
                }
            }
            worst
        }
    }
}

/// The distance from `c` to the segment from `a` to `b`.
fn nearest_on_segment(c: Point2, a: Point2, b: Point2) -> f64 {
    let d = b - a;
    let len2 = d.dot(d);
    let t = if len2 > 0.0 {
        ((c - a).dot(d) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (a + d * t - c).magnitude()
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;

    /// A slot outline, two straights and two half circles, splits into
    /// four runs alternating line and circle, the circles of the slot's
    /// radius.
    #[test]
    fn a_slot_splits_into_two_lines_and_two_arcs() {
        let mut points = Vec::new();
        for k in 0..=40 {
            points.push(Point2::new(f64::from(k) * 0.25, -1.0));
        }
        for k in 1..40 {
            let a = -core::f64::consts::FRAC_PI_2 + f64::from(k) / 40.0 * core::f64::consts::PI;
            points.push(Point2::new(10.0 + a.cos(), a.sin()));
        }
        for k in 0..=40 {
            points.push(Point2::new(10.0 - f64::from(k) * 0.25, 1.0));
        }
        for k in 1..=40 {
            let a = core::f64::consts::FRAC_PI_2 + f64::from(k) / 40.0 * core::f64::consts::PI;
            points.push(Point2::new(a.cos(), a.sin()));
        }
        let runs = runs(&points, 1e-3);
        assert_eq!(runs.len(), 4, "{runs:?}");
        for (k, run) in runs.iter().enumerate() {
            match run.shape {
                Shape::Line { .. } => assert_eq!(k % 2, 0),
                Shape::Circle { radius, .. } => {
                    assert_eq!(k % 2, 1);
                    assert!((radius - 1.0).abs() < 1e-3);
                }
            }
        }
    }
}
