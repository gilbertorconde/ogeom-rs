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
//!
//! [`fewest_runs`] splits the polyline into the fewest such runs instead,
//! for the segment fit: from every point it finds how far a line or a
//! circle reaches, and a greedy walk over those reaches is minimal. A line
//! reaches exactly as far as a line can, by the narrowest strip about the
//! points.

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
pub(crate) fn refit(kind: Shape, points: &[Point2], tolerance: f64) -> Option<Shape> {
    let shape = match kind {
        Shape::Line { .. } => line_through(points),
        Shape::Circle { .. } => circle_through(points),
    }?;
    (deviation(points, shape) <= tolerance).then_some(shape)
}

/// The sum of the points' squared distances from the shape.
pub(crate) fn squares(points: &[Point2], shape: Shape) -> f64 {
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
pub(crate) fn balance(points: &[Point2], runs: &mut [Run], k: usize, tolerance: f64) {
    let (a, b) = (runs[k - 1], runs[k]);
    if let Some((m, sa, sb, _)) = balanced(points, a, b, |kind, p| refit(kind, p, tolerance)) {
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

/// Where runs `a` and `b` would best share their point, as [`balance`]
/// finds it with `refit` giving each run's shape over a stretch, and the
/// two shapes and the squared residual there.
pub(crate) fn balanced(
    points: &[Point2],
    a: Run,
    b: Run,
    refit: impl Fn(Shape, &[Point2]) -> Option<Shape>,
) -> Option<(usize, Shape, Shape, f64)> {
    let least = |r: &Run| {
        if matches!(r.shape, Shape::Line { .. }) {
            1
        } else {
            2
        }
    };
    let first_fits = |m: usize| refit(a.shape, &points[a.first..=m]);
    let second_fits = |m: usize| refit(b.shape, &points[m..=b.last]);
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
    (l..=h)
        .filter_map(|m| cost(m).map(|(c, sa, sb)| (m, sa, sb, c)))
        .min_by(|x, y| x.3.total_cmp(&y.3))
}

/// How far the points stand from the chord between the first and the
/// last: a circle bowing no more than the tolerance is taken for a line.
pub(crate) fn bow(points: &[Point2]) -> f64 {
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
pub(crate) fn reach(
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

pub(crate) fn chord(a: Point2, b: Point2) -> Shape {
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
pub(crate) const fn count(n: usize) -> f64 {
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

/// How far one shape reaches from a start: it holds the points from the
/// start through `last`, and every stretch within them, within the
/// tolerance.
#[derive(Debug, Clone, Copy)]
struct Reach {
    last: usize,
    shape: Shape,
}

/// Split a polyline into the fewest runs a line or a circle each holds
/// within `tolerance`, among the runs the reaches below certify.
///
/// For an open polyline `ring` is its points and `period` is `None`; the
/// runs go from the first point to the last. For a closed one of `n`
/// points `ring` is them at least three times over and once more the
/// first, and `period` is `Some(n)`; the runs go once round, from a start
/// about the second turn to that start plus `n`, by index into `ring`.
///
/// From each start `i` the reach is the farthest point a line or a circle
/// fitted to the points from `i` holds within the tolerance, searched by
/// doubling then halving the step from the farthest reach of any earlier
/// start: a shape holding the points from an earlier start holds those
/// from `i` too. A line reaches exactly as far as any line can, since a
/// strip narrow enough for the points is narrow enough for any of them;
/// a circle reaches as far as the fit finds. The reaches never fall as
/// the start moves on, and any stretch within one is held by its shape,
/// so a walk taking the farthest reach from each run's end needs no more
/// runs than any split into stretches the reaches cover. A closed loop is
/// walked from every start within the first such run, one of which is a
/// break of the fewest.
pub(crate) fn fewest_runs(ring: &[Point2], period: Option<usize>, tolerance: f64) -> Vec<Run> {
    let n = period.unwrap_or_else(|| ring.len().saturating_sub(1));
    if ring.len() < 2 || n == 0 {
        return Vec::new();
    }
    let reaches = reaches(ring, n, period.is_some(), tolerance);
    let reach_at = |i: usize| -> Reach {
        let r = reaches[i % n];
        let shift = i / n * n;
        Reach {
            last: r.last + shift,
            ..r
        }
    };
    // The farthest reach of any start up to each point, never falling.
    let ends = if period.is_some() { 3 * n } else { n };
    let mut farthest: Vec<Reach> = Vec::with_capacity(ends);
    for i in 0..ends {
        let here = reach_at(i);
        let best = match farthest.last() {
            Some(b) if b.last >= here.last => *b,
            _ => here,
        };
        farthest.push(best);
    }
    let walk = |start: usize, stop: usize| -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        let mut at = start;
        while at < stop {
            let to = farthest[at].last.min(stop);
            out.push((at, to));
            at = to;
        }
        out
    };
    let spans = if period.is_some() {
        (n..=farthest[n].last.min(2 * n - 1))
            .map(|s| walk(s, s + n))
            .min_by_key(Vec::len)
            .unwrap_or_default()
    } else {
        walk(0, n)
    };
    let mut runs: Vec<Run> = spans
        .iter()
        .map(|&(first, last)| Run {
            first,
            last,
            shape: settle(&ring[first..=last], farthest[first].shape, tolerance),
        })
        .collect();
    relax(ring, &mut runs, period, tolerance);
    runs
}

/// Move each point consecutive runs share to where the two leave the
/// least squared residual, a line held by its least-squares fit or else by
/// its strip, until no move lowers the residual. A walk's runs reach as far
/// as they can, tight against the tolerance at their ends, so a run's fit
/// improves as the runs either side give way; each move lowers the total,
/// so the moves end.
fn relax(ring: &[Point2], runs: &mut [Run], period: Option<usize>, tolerance: f64) {
    let refit = |kind: Shape, points: &[Point2]| match kind {
        Shape::Line { .. } => line_through(points)
            .filter(|s| deviation(points, *s) <= tolerance)
            .or_else(|| line_fits(points, tolerance)),
        Shape::Circle { .. } => {
            circle_through(points).filter(|s| deviation(points, *s) <= tolerance)
        }
    };
    let m = runs.len();
    let joints = match period {
        Some(_) if m > 1 => m,
        _ => m.saturating_sub(1),
    };
    let mut moved = true;
    while moved {
        moved = false;
        for j in 0..joints {
            // Joint `j` is the end of run `j`; the last closes the loop,
            // the run after it numbered a turn on.
            let (k, shift) = if j + 1 < m {
                (j + 1, 0)
            } else {
                (0, period.unwrap_or(0))
            };
            let a = runs[j];
            let b = Run {
                first: runs[k].first + shift,
                last: runs[k].last + shift,
                ..runs[k]
            };
            let now = squares(&ring[a.first..=a.last], a.shape)
                + squares(&ring[b.first..=b.last], b.shape);
            let Some((at, sa, sb, cost)) = balanced(ring, a, b, refit) else {
                continue;
            };
            if cost < now {
                runs[j] = Run {
                    last: at,
                    shape: sa,
                    ..a
                };
                runs[k] = Run {
                    first: at - shift,
                    last: b.last - shift,
                    shape: sb,
                };
                moved = true;
            }
        }
    }
}

/// The reach from each of the first `count` starts. A closed ring's
/// stretches run at most once round.
fn reaches(ring: &[Point2], count: usize, closed: bool, tolerance: f64) -> Vec<Reach> {
    let mut out: Vec<Reach> = Vec::with_capacity(count);
    let mut held: Option<Reach> = None;
    for i in 0..count {
        let end = if closed { i + count } else { ring.len() - 1 };
        let base = match held {
            Some(r) if r.last > i => r,
            _ => Reach {
                last: i + 1,
                shape: chord(ring[i], ring[i + 1]),
            },
        };
        let mut best = base;
        let line = extend(base.last, end, |last| line_fits(&ring[i..=last], tolerance));
        if let Some(r) = line {
            best = r;
        }
        let circle = extend(best.last, end, |last| {
            let run = &ring[i..=last];
            circle_through(run).filter(|s| deviation(run, *s) <= tolerance)
        });
        if let Some(r) = circle {
            best = r;
        }
        held = Some(best);
        out.push(best);
    }
    out
}

/// The farthest `last` past `base`, up to `end`, for which `fit` gives a
/// shape, by doubling then halving the step from `base`.
fn extend(base: usize, end: usize, fit: impl Fn(usize) -> Option<Shape>) -> Option<Reach> {
    let mut good: Option<Reach> = None;
    let mut step = 1;
    let bad = loop {
        let at = (base + step).min(end);
        if at <= good.map_or(base, |g| g.last) {
            return good;
        }
        match fit(at) {
            Some(shape) => {
                good = Some(Reach { last: at, shape });
                if at == end {
                    return good;
                }
                step *= 2;
            }
            None => break at,
        }
    };
    let (mut lo, mut hi) = (good.map_or(base, |g| g.last), bad);
    while hi - lo > 1 {
        let mid = lo + (hi - lo) / 2;
        match fit(mid) {
            Some(shape) => {
                lo = mid;
                good = Some(Reach { last: mid, shape });
            }
            None => hi = mid,
        }
    }
    good
}

/// A line holding the points within `tolerance`, if any does: the
/// least-squares line where that holds them, and otherwise the middle of
/// the narrowest strip about them. Three of the points spanning a triangle
/// wider than the strip may be leave no line to look for.
pub(crate) fn line_fits(points: &[Point2], tolerance: f64) -> Option<Shape> {
    if let Some(s) = line_through(points).filter(|s| deviation(points, *s) <= tolerance) {
        return Some(s);
    }
    if wider_than(points, 2.0 * tolerance) {
        return None;
    }
    strip(points).filter(|s| deviation(points, *s) <= tolerance)
}

/// Whether the triangle of the first point, the last and the point
/// farthest from the line through them is wider than `width`: the width of
/// a set is at least that of any three of its points, and a triangle's
/// width is its least height.
fn wider_than(points: &[Point2], width: f64) -> bool {
    let (Some(&a), Some(&b)) = (points.first(), points.last()) else {
        return false;
    };
    let Some(&c) = points.iter().max_by(|p, q| {
        (b - a)
            .cross(**p - a)
            .abs()
            .total_cmp(&(b - a).cross(**q - a).abs())
    }) else {
        return false;
    };
    let height = |p: Point2, q: Point2, r: Point2| {
        let base = (r - q).magnitude();
        if base > 0.0 {
            (r - q).cross(p - q).abs() / base
        } else {
            0.0
        }
    };
    height(c, a, b).min(height(a, b, c)).min(height(b, c, a)) > width
}

/// The shape a run takes: a line where one holds it, the least-squares one
/// where that does; otherwise the least-squares circle where it holds, and
/// otherwise the shape whose reach covered the run.
fn settle(points: &[Point2], reached: Shape, tolerance: f64) -> Shape {
    if let Some(strip) = line_fits(points, tolerance) {
        return line_through(points)
            .filter(|s| deviation(points, *s) <= tolerance)
            .unwrap_or(strip);
    }
    circle_through(points)
        .filter(|s| deviation(points, *s) <= tolerance)
        .unwrap_or(reached)
}

/// The middle line of the narrowest strip holding the points, along the
/// direction from the first to the last. The narrowest strip has a side
/// along an edge of the points' convex hull, so each edge is tried
/// against the hull vertex farthest from it, found by walking round.
pub(crate) fn strip(points: &[Point2]) -> Option<Shape> {
    let (first, last) = (*points.first()?, *points.last()?);
    let hull = hull(points);
    let along = if hull.len() < 3 {
        let [a, b] = [*hull.first()?, *hull.last()?];
        if a == b {
            return None;
        }
        (b - a) / (b - a).magnitude()
    } else {
        let h = hull.len();
        let height = |i: usize, j: usize| (hull[(i + 1) % h] - hull[i]).cross(hull[j] - hull[i]);
        let mut j = 1;
        let mut best = (f64::INFINITY, Vector2::new(1.0, 0.0));
        for i in 0..h {
            let edge = hull[(i + 1) % h] - hull[i];
            let len = edge.magnitude();
            if len == 0.0 {
                continue;
            }
            let mut steps = 0;
            while steps < h && height(i, (j + 1) % h) > height(i, j) {
                j = (j + 1) % h;
                steps += 1;
            }
            let width = height(i, j) / len;
            if width < best.0 {
                best = (width, edge / len);
            }
        }
        best.1
    };
    let along = if along.dot(last - first) < 0.0 {
        -along
    } else {
        along
    };
    let normal = along.perpendicular();
    let (lo, hi) = points
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
            let d = normal.dot(*p - first);
            (lo.min(d), hi.max(d))
        });
    Some(Shape::Line {
        at: first + normal * (0.5 * (lo + hi)),
        along,
    })
}

/// The convex hull of the points, counter-clockwise, without points
/// along its edges: one point, or two, where the points are on a line.
fn hull(points: &[Point2]) -> Vec<Point2> {
    let mut sorted = points.to_vec();
    sorted.sort_by(|a, b| a.x.total_cmp(&b.x).then(a.y.total_cmp(&b.y)));
    sorted.dedup();
    if sorted.len() < 3 {
        return sorted;
    }
    let turn = |o: Point2, a: Point2, b: Point2| (a - o).cross(b - o);
    let mut out: Vec<Point2> = Vec::with_capacity(sorted.len() + 1);
    let reversed: Vec<Point2> = sorted.iter().rev().copied().collect();
    for order in [&sorted, &reversed] {
        let floor = out.len();
        for &p in order {
            while out.len() >= floor + 2 && turn(out[out.len() - 2], out[out.len() - 1], p) <= 0.0 {
                out.pop();
            }
            out.push(p);
        }
        out.pop();
    }
    if out.len() < 3 {
        let (a, b) = (sorted[0], sorted[sorted.len() - 1]);
        return vec![a, b];
    }
    out
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

    /// The fewest lines by brute force: over every split, each stretch
    /// held by the narrowest strip about it.
    fn fewest_lines(points: &[Point2], tolerance: f64) -> usize {
        let n = points.len();
        let mut best = vec![usize::MAX; n];
        best[0] = 0;
        for j in 1..n {
            for i in 0..j {
                if best[i] != usize::MAX && line_fits(&points[i..=j], tolerance).is_some() {
                    best[j] = best[j].min(best[i] + 1);
                }
            }
        }
        best[n - 1]
    }

    /// Sides of a few jittered points each, turning sharply between them,
    /// so that no circle holds more than a corner's worth.
    fn wander(seed: u64, n: usize) -> Vec<Point2> {
        let mut s = seed;
        let mut next = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            count(usize::try_from(s >> 11).unwrap()) / count(1 << 53)
        };
        let (mut at, mut heading) = (Point2::new(0.0, 0.0), 0.0f64);
        let mut points = Vec::with_capacity(n);
        while points.len() < n {
            heading += if next() < 0.5 { -1.0 } else { 1.0 } * (1.0 + next());
            let side = 2 + (0..5).filter(|_| next() < 0.5).count();
            for _ in 0..side {
                at += Vector2::new(heading.cos(), heading.sin()) * (0.5 + 0.5 * next());
                points.push(at + Vector2::new(next() - 0.5, next() - 0.5) * 0.2);
            }
        }
        points.truncate(n);
        points
    }

    /// On wandering polylines, open and closed, the walk over the reaches
    /// never needs more runs than the fewest lines found by trying every
    /// split, and every run holds its points.
    #[test]
    fn the_walk_needs_no_more_runs_than_the_fewest_lines() {
        let tolerance = 0.3;
        for seed in 1..40u64 {
            let points = wander(seed.wrapping_mul(0x9e37_79b9_7f4a_7c15), 30);
            let runs = fewest_runs(&points, None, tolerance);
            assert!(
                runs.len() <= fewest_lines(&points, tolerance),
                "seed {seed}"
            );
            let n = points.len();
            let ring: Vec<Point2> = (0..=3 * n).map(|k| points[k % n]).collect();
            let closed = fewest_runs(&ring, Some(n), tolerance);
            let fewest = (0..n)
                .map(|s| {
                    let turned: Vec<Point2> = (0..=n).map(|k| points[(s + k) % n]).collect();
                    fewest_lines(&turned, tolerance)
                })
                .min()
                .unwrap();
            assert!(closed.len() <= fewest, "seed {seed}");
            assert_eq!(closed[closed.len() - 1].last, closed[0].first + n);
            for run in runs.iter().chain(&closed) {
                let held = &ring[run.first..=run.last];
                assert!(deviation(held, run.shape) <= tolerance, "seed {seed}");
            }
        }
    }
}
