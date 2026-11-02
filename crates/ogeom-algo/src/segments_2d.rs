//! A chain of lines and arcs fitted to a 2D polyline within a tolerance.
//!
//! The fit has four stages, and the last decides what is returned.
//!
//! 1. **The fewest pieces.** From every point the fit finds how far one
//!    line or one circle holds the polyline within the tolerance: a line
//!    exactly as far as any line can, by the narrowest strip about the
//!    points, a circle as far as its least-squares fit holds them. A shape
//!    holding the points from one start holds those from any later start
//!    within its reach, so the reaches never fall as the start moves on,
//!    and a walk taking the farthest reach from each run's end needs no
//!    more runs than any split into stretches the reaches cover. A closed
//!    loop is walked from every start within one run, one of which is a
//!    break of the fewest. The point consecutive runs share is then moved
//!    to where their fits leave the least squared residual. A loop one
//!    circle holds is that circle.
//! 2. **The shapes.** A circle run bowing from its chord by no more than
//!    `min_bow` is split into the fewest lines instead. A line within
//!    `axis_snap` of an axis turns onto it about its points' centroid where
//!    it still holds them. A circle run next to lines is refitted tangent
//!    to both, or else to one, where the tangent circle still holds it.
//! 3. **The joints.** Consecutive pieces meet at one point: the tangent
//!    point of a tangent pair; with `corners`, where the two curves cross
//!    nearest the input point the runs share, so a corner the polyline
//!    rounds off comes back sharp; otherwise midway between the nearest
//!    points of the two curves to that shared point. A crossing moves no
//!    piece's end by more than half its run. An open chain starts and ends
//!    at the nearest points of its end pieces to the end points.
//! 4. **The check.** Every input point must lie within the tolerance of
//!    the piece whose range holds it, the piece being the line segment or
//!    the arc itself, not its whole curve. The point two ranges share may
//!    move to any point both pieces hold up to, since a crossing cuts one
//!    piece short of points the other then holds. A joint where there is
//!    none takes its next choice, the midway point and then the shared
//!    input point itself, and a joint out of choices turns the pieces
//!    either side into the polylines through their points, which hold them
//!    exactly.
//!
//! The tolerance is therefore a hard bound on what is returned. What is
//! minimal is the count from the first stage: no split into stretches a
//! line holds, or that lie within the reach of a circle fitted from a
//! point at or before them, has fewer pieces, and in particular no split
//! into lines alone. A circle the least-squares fit misses is not found,
//! and the arcs `min_bow` turns into lines and the polyline fallbacks add
//! pieces to that count.
//!
//! Each point's reach is a fit over the run it starts, so the work grows
//! as the number of points times the points in a run.

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::{Point2, Vector2};

use crate::fit::segments::{
    Run, Shape, bow, chord, circle_through, count, deviation, fewest_runs, line_fits, line_through,
    reach,
};
use crate::outline::{Path, Piece2, meet};

/// How a polyline is split into lines and arcs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segments {
    /// Every input point lies within this distance of its piece.
    pub tolerance: f64,
    /// Whether pieces turning against each other meet where their curves
    /// cross, bringing back the sharp corner the polyline rounds off,
    /// rather than midway between their nearest points to the input point
    /// they share.
    pub corners: bool,
    /// Lines within this angle of the `x` or the `y` axis, in radians up
    /// to a quarter turn, are turned onto it where they still hold their
    /// points.
    pub axis_snap: Option<f64>,
    /// An arc bowing from its chord by no more than this is taken as
    /// lines. Arcs are only fitted where a line does not hold the points,
    /// so this matters only above the tolerance.
    pub min_bow: f64,
    /// Whether the polyline closes back to its first point, which it then
    /// need not repeat.
    pub closed: bool,
}

/// One piece of a fitted chain.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segment2 {
    /// The line or the arc: a whole circle where one holds a closed loop.
    pub piece: Piece2,
    /// The input points it holds, first and last by index. Consecutive
    /// pieces share the point between them. On a closed polyline a range
    /// can run on past the last point to the first, so that `first >
    /// last`, and a whole circle has `first == last`.
    pub range: (usize, usize),
    /// The root mean square distance of those points from the piece.
    pub sigma: f64,
}

/// Fit a chain of lines and arcs to a polyline, every point within the
/// tolerance of its piece, consecutive pieces meeting end to end and, for a
/// closed polyline, the last meeting the first. See the module
/// documentation for how, and for what is minimal.
///
/// Repeated consecutive points count once, and a closed polyline's first
/// point repeated at its end is dropped.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// the tolerance is not positive, `min_bow` is negative, `axis_snap` is
/// not an angle from zero to a quarter turn, or there are fewer than two
/// distinct points, three for a closed polyline, and
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if a point is
/// not finite.
pub fn fit_segments_2d(points: &[Point2], options: &Segments) -> OgeomResult<Vec<Segment2>> {
    let Segments {
        tolerance,
        corners,
        axis_snap,
        min_bow,
        closed,
    } = *options;
    if !(tolerance > 0.0 && tolerance.is_finite()) {
        ogeom_bail!(
            Construction,
            "a segment tolerance of {tolerance} is no bound"
        );
    }
    if min_bow.is_nan() || min_bow < 0.0 {
        ogeom_bail!(Construction, "a least bow of {min_bow} is not a distance");
    }
    if let Some(a) = axis_snap
        && !(0.0..=core::f64::consts::FRAC_PI_4).contains(&a)
    {
        ogeom_bail!(
            Construction,
            "snapping lines within {a} rad of an axis is not an angle up to a quarter turn"
        );
    }
    if let Some(p) = points
        .iter()
        .find(|p| !(p.x.is_finite() && p.y.is_finite()))
    {
        ogeom_bail!(Numeric, "the polyline point {p:?} is not finite");
    }
    let (distinct, index) = distinct(points, closed);
    let n = distinct.len();
    if n < 2 || (closed && n < 3) {
        ogeom_bail!(
            Construction,
            "{n} distinct points are too few for a {} polyline",
            if closed { "closed" } else { "open" }
        );
    }
    if closed && let Some(circle) = whole_circle(&distinct, tolerance) {
        let sigma = rms(distinct.iter().map(|p| distance_to_piece(*p, &circle)));
        return Ok(vec![Segment2 {
            piece: circle,
            range: (index[0], index[0]),
            sigma,
        }]);
    }
    // A closed loop goes round four times: its runs take three, and the
    // joints look half a run past either end.
    let ring: Vec<Point2> = if closed {
        (0..=4 * n).map(|k| distinct[k % n]).collect()
    } else {
        distinct.clone()
    };
    let mut runs = fewest_runs(&ring, closed.then_some(n), tolerance);
    if closed && let Some(first) = runs.first().map(|r| r.first) {
        // Numbered from the second turn.
        let to = first % n + n;
        for run in &mut runs {
            (run.first, run.last) = (run.first + to - first, run.last + to - first);
        }
    }
    if closed && runs.len() == 1 {
        return Ok(whole_loop(&distinct, &index, runs[0].shape));
    }
    runs = straighten(&ring, runs, min_bow, tolerance);
    if let Some(angle) = axis_snap {
        for run in &mut runs {
            snap(&ring[run.first..=run.last], run, angle, tolerance);
        }
    }
    let tangent = tangents(&ring, &mut runs, closed, tolerance);
    let chain = Chain {
        ring: &ring,
        period: n,
        runs: &runs,
        closed,
        tolerance,
    };
    let pieces = chain.pieces(&tangent, corners);
    let original = |k: usize| index[if closed { k % n } else { k }];
    Ok(pieces
        .into_iter()
        .map(|(piece, first, last)| Segment2 {
            piece,
            range: (original(first), original(last)),
            sigma: rms(ring[first..=last]
                .iter()
                .map(|p| distance_to_piece(*p, &piece))),
        })
        .collect())
}

/// A closed loop one run holds whole: its circle, or for a loop folded
/// flat along a line, the polyline through its points.
fn whole_loop(points: &[Point2], index: &[usize], shape: Shape) -> Vec<Segment2> {
    if let Shape::Circle { centre, radius } = shape {
        let mut closed = points.to_vec();
        closed.push(points[0]);
        let circle = Piece2::Circle {
            centre,
            radius,
            counter_clockwise: turn_about(&closed, centre) > 0.0,
        };
        let sigma = rms(points.iter().map(|p| distance_to_piece(*p, &circle)));
        return vec![Segment2 {
            piece: circle,
            range: (index[0], index[0]),
            sigma,
        }];
    }
    let n = points.len();
    (0..n)
        .map(|k| Segment2 {
            piece: Piece2::Line {
                from: points[k],
                to: points[(k + 1) % n],
            },
            range: (index[k], index[(k + 1) % n]),
            sigma: 0.0,
        })
        .collect()
}

/// The points with consecutive repeats dropped, a closed polyline's
/// repeated first point too, and the input index of each.
fn distinct(points: &[Point2], closed: bool) -> (Vec<Point2>, Vec<usize>) {
    let (mut out, mut index) = (Vec::with_capacity(points.len()), Vec::new());
    for (k, p) in points.iter().enumerate() {
        if out.last() != Some(p) {
            out.push(*p);
            index.push(k);
        }
    }
    while closed && out.len() > 1 && out.first() == out.last() {
        out.pop();
        index.pop();
    }
    (out, index)
}

/// The circle holding a whole closed loop within the tolerance, if the
/// fit finds one.
fn whole_circle(points: &[Point2], tolerance: f64) -> Option<Piece2> {
    let Shape::Circle { centre, radius } = circle_through(points)? else {
        return None;
    };
    let mut closed = points.to_vec();
    closed.push(points[0]);
    (deviation(&closed, Shape::Circle { centre, radius }) <= tolerance).then(|| Piece2::Circle {
        centre,
        radius,
        counter_clockwise: turn_about(&closed, centre) > 0.0,
    })
}

/// The signed angle the points turn through about a centre.
fn turn_about(points: &[Point2], centre: Point2) -> f64 {
    points
        .windows(2)
        .map(|w| {
            let (u, v) = (w[0] - centre, w[1] - centre);
            u.cross(v).atan2(u.dot(v))
        })
        .sum()
}

/// The circle runs bowing no more than `min_bow`, each split into the
/// fewest lines: a line reaches exactly as far as any can, so taking the
/// farthest reach each time is the fewest.
fn straighten(ring: &[Point2], runs: Vec<Run>, min_bow: f64, tolerance: f64) -> Vec<Run> {
    let mut out = Vec::with_capacity(runs.len());
    for run in runs {
        let points = &ring[run.first..=run.last];
        if matches!(run.shape, Shape::Line { .. }) || bow(points) > min_bow {
            out.push(run);
            continue;
        }
        let mut first = run.first;
        while first < run.last {
            let (last, strip) = reach(first, run.last, 1, |last| {
                line_fits(&ring[first..=last], tolerance)
            });
            let points = &ring[first..=last];
            let shape = line_through(points)
                .filter(|s| deviation(points, *s) <= tolerance)
                .or(strip)
                .unwrap_or_else(|| chord(ring[first], ring[last]));
            out.push(Run { first, last, shape });
            first = last;
        }
    }
    out
}

/// Turn a line run onto the axis it lies within `angle` of, about its
/// points' centroid, where the turned line still holds them.
fn snap(points: &[Point2], run: &mut Run, angle: f64, tolerance: f64) {
    let Shape::Line { along, .. } = run.shape else {
        return;
    };
    let axis = if along.x.abs() >= along.y.abs() {
        Vector2::new(along.x.signum(), 0.0)
    } else {
        Vector2::new(0.0, along.y.signum())
    };
    if along.cross(axis).abs().asin() > angle {
        return;
    }
    let sum = points
        .iter()
        .fold(Vector2::new(0.0, 0.0), |s, p| s + p.to_vector());
    let centroid = Point2::from_vector(sum / count(points.len()));
    let snapped = Shape::Line {
        at: centroid,
        along: axis,
    };
    if deviation(points, snapped) <= tolerance {
        run.shape = snapped;
    }
}

/// Refit each circle run next to lines tangent to both, or else to the
/// one before, or else the one after, where the tangent circle holds the
/// run. Whether the pieces meeting at each joint are tangent: joint `k` is
/// the start of run `k`, and for an open chain joint `runs.len()` its end.
fn tangents(ring: &[Point2], runs: &mut [Run], closed: bool, tolerance: f64) -> Vec<bool> {
    let m = runs.len();
    let mut tangent = vec![false; m + 1];
    let line = |run: &Run| match run.shape {
        Shape::Line { at, along } => Some((at, along)),
        Shape::Circle { .. } => None,
    };
    for k in 0..m {
        let Shape::Circle { centre, radius } = runs[k].shape else {
            continue;
        };
        // A loop of two runs has one neighbour, on both sides.
        let before = (k > 0 || closed)
            .then(|| (k + m - 1) % m)
            .filter(|&b| b != k);
        let after = (k + 1 < m || closed)
            .then(|| (k + 1) % m)
            .filter(|&a| a != k && Some(a) != before);
        let before = before.and_then(|b| line(&runs[b]));
        let after = after.and_then(|a| line(&runs[a]));
        let points = &ring[runs[k].first..=runs[k].last];
        let holds = |s: Option<Shape>| s.filter(|s| deviation(points, *s) <= tolerance);
        let both = before
            .zip(after)
            .and_then(|(b, a)| holds(tangent_to_lines(points, centre, b, a)));
        let fitted = both
            .map(|s| (s, true, true))
            .or_else(|| {
                before.and_then(|b| {
                    holds(tangent_to_line(points, centre, radius, b)).map(|s| (s, true, false))
                })
            })
            .or_else(|| {
                after.and_then(|a| {
                    holds(tangent_to_line(points, centre, radius, a)).map(|s| (s, false, true))
                })
            });
        if let Some((shape, at_start, at_end)) = fitted {
            runs[k].shape = shape;
            tangent[k] |= at_start;
            if at_end {
                tangent[if closed { (k + 1) % m } else { k + 1 }] = true;
            }
        }
    }
    tangent
}

/// The unit normal of a line toward the side `p` lies on, `None` on it.
fn toward(at: Point2, along: Vector2, p: Point2) -> Option<Vector2> {
    let normal = along.perpendicular();
    let side = normal.dot(p - at);
    (side != 0.0).then(|| normal * side.signum())
}

/// The circle nearest the points tangent to a line, its centre on the side
/// of `centre`: Gauss-Newton on the centre's place along the line and the
/// radius, from the free fit.
fn tangent_to_line(
    points: &[Point2],
    centre: Point2,
    radius: f64,
    (at, along): (Point2, Vector2),
) -> Option<Shape> {
    let normal = toward(at, along, centre)?;
    let (mut s, mut r) = (along.dot(centre - at), radius);
    for _ in 0..8 {
        let c = at + along * s + normal * r;
        let (mut m, mut v) = ([[0.0f64; 2]; 2], [0.0f64; 2]);
        for p in points {
            let d = *p - c;
            let dist = d.magnitude();
            if dist == 0.0 {
                return None;
            }
            let e = d / dist;
            let row = [-e.dot(along), -e.dot(normal) - 1.0];
            let res = dist - r;
            for i in 0..2 {
                for j in 0..2 {
                    m[i][j] += row[i] * row[j];
                }
                v[i] -= row[i] * res;
            }
        }
        let det = m[0][0].mul_add(m[1][1], -m[0][1] * m[1][0]);
        if det.is_nan() || det == 0.0 {
            break;
        }
        s += (v[0] * m[1][1] - v[1] * m[0][1]) / det;
        r += (m[0][0] * v[1] - m[1][0] * v[0]) / det;
        if r.is_nan() || r <= 0.0 {
            return None;
        }
    }
    let c = at + along * s + normal * r;
    (s.is_finite() && r.is_finite()).then_some(Shape::Circle {
        centre: c,
        radius: r,
    })
}

/// The circle nearest the points tangent to two lines, its centre on the
/// side of `centre` of each. The centres equally far from both lines lie
/// on one line, a bisector or, for parallel lines, the line midway; the
/// centre's place along it is found by Gauss-Newton from the free fit's.
fn tangent_to_lines(
    points: &[Point2],
    centre: Point2,
    (a1, u1): (Point2, Vector2),
    (a2, u2): (Point2, Vector2),
) -> Option<Shape> {
    let (n1, n2) = (toward(a1, u1, centre)?, toward(a2, u2, centre)?);
    // n1 . (c - a1) = n2 . (c - a2) = r, about the free centre.
    let g = n1 - n2;
    let g2 = g.dot(g);
    if g2.is_nan() || g2 == 0.0 {
        return None;
    }
    let k = n1.dot(a1 - centre) - n2.dot(a2 - centre);
    let base = centre + g * (k / g2);
    let b = g.perpendicular() / g2.sqrt();
    let (r0, rate) = (n1.dot(base - a1), n1.dot(b));
    let mut t = 0.0;
    for _ in 0..8 {
        let (c, r) = (base + b * t, rate.mul_add(t, r0));
        if r.is_nan() || r <= 0.0 {
            return None;
        }
        let (mut jj, mut jr) = (0.0, 0.0);
        for p in points {
            let d = *p - c;
            let dist = d.magnitude();
            if dist == 0.0 {
                return None;
            }
            let row = -(d / dist).dot(b) - rate;
            jj += row * row;
            jr += row * (dist - r);
        }
        if jj.is_nan() || jj == 0.0 {
            break;
        }
        t -= jr / jj;
    }
    let r = rate.mul_add(t, r0);
    (r > 0.0 && t.is_finite()).then(|| Shape::Circle {
        centre: base + b * t,
        radius: r,
    })
}

fn path(shape: Shape) -> Path {
    match shape {
        Shape::Line { at, along } => Path::Line { at, along },
        Shape::Circle { centre, radius } => Path::Circle {
            centre,
            radius,
            ccw: true,
        },
    }
}

/// The runs and the points they hold, joined into pieces.
struct Chain<'a> {
    ring: &'a [Point2],
    /// The number of points of a closed loop, which the ring repeats.
    period: usize,
    runs: &'a [Run],
    closed: bool,
    tolerance: f64,
}

impl Chain<'_> {
    fn joints(&self) -> usize {
        self.runs.len() + usize::from(!self.closed)
    }

    /// The pieces, each with the range it holds, every input point within
    /// the tolerance of its own.
    ///
    /// Each joint's point is taken from its choices, best first, and the
    /// input point the two runs share moves to any point both pieces hold
    /// up to: a crossing cuts one piece short of the points the other then
    /// takes. A joint where no such point is found moves to its next
    /// choice, and one out of choices turns the pieces either side into the
    /// polylines through their runs.
    fn pieces(&self, tangent: &[bool], corners: bool) -> Vec<(Piece2, usize, usize)> {
        let m = self.runs.len();
        if m == 0 {
            return Vec::new();
        }
        let joints = self.joints();
        let choices: Vec<Vec<Point2>> = (0..joints)
            .map(|j| self.choices(j, tangent[j], corners))
            .collect();
        let mut level = vec![0usize; joints];
        let mut exact = vec![false; m];
        loop {
            for k in (0..m).filter(|&k| exact[k]) {
                for j in [k, self.end(k)] {
                    level[j] = choices[j].len() - 1;
                }
            }
            let built: Vec<Option<Piece2>> = (0..m)
                .map(|k| {
                    let (from, to) = (k, self.end(k));
                    if exact[k] {
                        None
                    } else {
                        self.piece(k, choices[from][level[from]], choices[to][level[to]])
                    }
                })
                .collect();
            let windows: Vec<Option<(usize, usize)>> = (0..m)
                .map(|k| {
                    let run = self.runs[k];
                    if exact[k] {
                        Some((run.first, run.last))
                    } else {
                        built[k].and_then(|p| self.window(k, &p))
                    }
                })
                .collect();
            let splits: Vec<Option<usize>> = (0..joints).map(|j| self.split(j, &windows)).collect();
            let mut failing = (0..joints).filter(|&j| splits[j].is_none()).peekable();
            if failing.peek().is_none() {
                return self.emit(&built, &splits);
            }
            for j in failing {
                if level[j] + 1 < choices[j].len() {
                    level[j] += 1;
                } else {
                    for k in self.beside(j) {
                        exact[k] = true;
                    }
                }
            }
        }
    }

    /// The joint at the end of run `k`.
    const fn end(&self, k: usize) -> usize {
        if self.closed {
            (k + 1) % self.runs.len()
        } else {
            k + 1
        }
    }

    /// The runs either side of joint `j`.
    fn beside(&self, j: usize) -> Vec<usize> {
        let m = self.runs.len();
        if self.closed {
            vec![(j + m - 1) % m, j]
        } else {
            (j.saturating_sub(1)..=j.min(m - 1)).collect()
        }
    }

    /// The middle point of run `k`, by index, in run `of`'s numbering.
    fn middle(&self, k: usize, of: usize) -> usize {
        let run = self.runs[k];
        let n = self.period;
        let mid = run.first + (run.last - run.first) / 2;
        let m = self.runs.len();
        if self.closed && k == m - 1 && of == 0 {
            mid - n
        } else if self.closed && k == 0 && of == m - 1 {
            mid + n
        } else {
            mid
        }
    }

    /// The stretch about the middle of run `k` that the piece holds
    /// within the tolerance, reaching no further than the middles of the
    /// runs either side: `None` if it does not hold the middle itself.
    fn window(&self, k: usize, piece: &Piece2) -> Option<(usize, usize)> {
        let m = self.runs.len();
        let holds = |i: usize| distance_to_piece(self.ring[i], piece) <= self.tolerance;
        let mid = self.middle(k, k);
        if !holds(mid) {
            return None;
        }
        let lower = if self.closed || k > 0 {
            self.middle((k + m - 1) % m, k)
        } else {
            self.runs[0].first
        };
        let upper = if self.closed || k + 1 < m {
            self.middle((k + 1) % m, k)
        } else {
            self.runs[m - 1].last
        };
        let (mut lo, mut hi) = (mid, mid);
        while lo > lower && holds(lo - 1) {
            lo -= 1;
        }
        while hi < upper && holds(hi + 1) {
            hi += 1;
        }
        Some((lo, hi))
    }

    /// Where the range of the run after joint `j` starts, in its own
    /// numbering: a point both pieces hold up to, nearest the one the
    /// runs share.
    fn split(&self, j: usize, windows: &[Option<(usize, usize)>]) -> Option<usize> {
        let m = self.runs.len();
        if !self.closed && j == 0 {
            let first = self.runs[0].first;
            return windows[0].filter(|w| w.0 <= first).map(|_| first);
        }
        if !self.closed && j == m {
            let last = self.runs[m - 1].last;
            return windows[m - 1].filter(|w| w.1 >= last).map(|_| last);
        }
        let (a, b) = ((j + m - 1) % m, j % m);
        let (wa, wb) = (windows[a]?, windows[b]?);
        // The run before the loop closes is numbered a turn on.
        let shift = if self.closed && j == 0 {
            self.period
        } else {
            0
        };
        let (lo, hi) = (
            wa.0.saturating_sub(shift).max(wb.0),
            (wa.1 - shift).min(wb.1),
        );
        let shared = self.runs[b].first;
        (lo <= hi).then(|| shared.clamp(lo, hi))
    }

    /// The pieces and their ranges, the runs left without a piece as the
    /// polylines through their points.
    fn emit(
        &self,
        built: &[Option<Piece2>],
        splits: &[Option<usize>],
    ) -> Vec<(Piece2, usize, usize)> {
        let m = self.runs.len();
        let n = self.period;
        let mut out = Vec::with_capacity(m);
        for (k, piece) in built.iter().enumerate() {
            let run = self.runs[k];
            let first = splits[k].unwrap_or(run.first);
            let last = splits[self.end(k)].map_or(run.last, |s| {
                if self.closed && k == m - 1 { s + n } else { s }
            });
            match piece {
                Some(p) => out.push((*p, first, last)),
                // Both ends are on input points: the polyline holds every
                // point of the range exactly.
                None => out.extend((run.first..run.last).map(|i| {
                    (
                        Piece2::Line {
                            from: self.ring[i],
                            to: self.ring[i + 1],
                        },
                        i,
                        i + 1,
                    )
                })),
            }
        }
        out
    }

    /// Where the pieces at joint `j` may meet, best first, the input point
    /// they share last.
    fn choices(&self, j: usize, tangent: bool, corners: bool) -> Vec<Point2> {
        let m = self.runs.len();
        if !self.closed && (j == 0 || j == m) {
            let (run, p) = if j == 0 {
                (self.runs[0], self.ring[self.runs[0].first])
            } else {
                (self.runs[m - 1], self.ring[self.runs[m - 1].last])
            };
            return path(run.shape).foot(p).into_iter().chain([p]).collect();
        }
        let (a, b) = (self.runs[(j + m - 1) % m], self.runs[j % m]);
        let q = self.ring[b.first];
        let (pa, pb) = (path(a.shape), path(b.shape));
        let near = pa.foot(q).zip(pb.foot(q)).map(|(x, y)| x.lerp(y, 0.5));
        let mut out = Vec::with_capacity(3);
        if tangent {
            // The circle touches the line where its centre's foot is.
            let touch = match (a.shape, b.shape) {
                (Shape::Circle { centre, .. }, Shape::Line { .. }) => pb.foot(centre),
                (Shape::Line { .. }, Shape::Circle { centre, .. }) => pa.foot(centre),
                _ => None,
            };
            out.extend(touch);
        } else if corners {
            let crossing = meet(pa, pb)
                .into_iter()
                .filter(|x| {
                    let past_a = self.past_end(&a, *x);
                    let before_b = self.past_end(&reversed(b), *x);
                    past_a.zip(before_b).is_some_and(|(s, t)| {
                        (-0.5..=0.5).contains(&s) && (-0.5..=0.5).contains(&t)
                    })
                })
                .min_by(|x, y| x.square_distance(q).total_cmp(&y.square_distance(q)));
            out.extend(crossing);
        }
        out.extend(near);
        out.push(q);
        out
    }

    /// How far past the end of a run a point of its curve lies, as a
    /// fraction of the run: negative before the end.
    fn past_end(&self, run: &Run, x: Point2) -> Option<f64> {
        let (from, to) = (self.ring[run.first], self.ring[run.last]);
        match run.shape {
            Shape::Line { along, .. } => {
                let span = (to - from).dot(along);
                (span != 0.0).then(|| (x - to).dot(along) / span)
            }
            Shape::Circle { centre, .. } => {
                let turn = self.turn(run, centre);
                let (u, v) = (to - centre, x - centre);
                (turn != 0.0).then(|| u.cross(v).atan2(u.dot(v)) / turn)
            }
        }
    }

    /// The signed angle a run's points turn through about a centre, from
    /// its first to its last.
    fn turn(&self, run: &Run, centre: Point2) -> f64 {
        let (lo, hi) = (run.first.min(run.last), run.first.max(run.last));
        let turn = turn_about(&self.ring[lo..=hi], centre);
        if run.first <= run.last { turn } else { -turn }
    }

    /// The piece of run `k` from one joint point to the next.
    fn piece(&self, k: usize, from: Point2, to: Point2) -> Option<Piece2> {
        let run = self.runs[k];
        match run.shape {
            Shape::Line { along, .. } => {
                ((to - from).dot(along) > 0.0).then_some(Piece2::Line { from, to })
            }
            Shape::Circle { centre, .. } => arc_through(centre, from, to, self.turn(&run, centre)),
        }
    }
}

/// The run with its ends swapped, so that its start is its end.
const fn reversed(run: Run) -> Run {
    Run {
        first: run.last,
        last: run.first,
        shape: run.shape,
    }
}

/// The arc from one point to another about the centre nearest `centre`
/// equally far from both, turning the way `turn` does.
fn arc_through(centre: Point2, from: Point2, to: Point2, turn: f64) -> Option<Piece2> {
    let d = to - from;
    let len = d.magnitude();
    if len == 0.0 || turn == 0.0 {
        return None;
    }
    let middle = from.lerp(to, 0.5);
    let w = d.perpendicular() / len;
    let c = middle + w * w.dot(centre - middle);
    let (u, v) = (from - c, to - c);
    let mut sweep = u.cross(v).atan2(u.dot(v));
    let tau = core::f64::consts::TAU;
    if turn > 0.0 && sweep <= 0.0 {
        sweep += tau;
    } else if turn < 0.0 && sweep >= 0.0 {
        sweep -= tau;
    }
    (sweep.abs() < tau).then(|| Piece2::Arc {
        centre: c,
        radius: u.magnitude(),
        start: u.y.atan2(u.x),
        sweep,
    })
}

/// The distance from a point to a piece: to the segment, or to the arc
/// within its sweep and to its nearer end outside it.
fn distance_to_piece(p: Point2, piece: &Piece2) -> f64 {
    match *piece {
        Piece2::Line { from, to } => {
            let d = to - from;
            let len2 = d.dot(d);
            let t = if len2 > 0.0 {
                ((p - from).dot(d) / len2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            from.lerp(to, t).distance(p)
        }
        Piece2::Arc {
            centre,
            radius,
            start,
            sweep,
        } => {
            let d = p - centre;
            let tau = core::f64::consts::TAU;
            let along = ((d.y.atan2(d.x) - start) * sweep.signum()).rem_euclid(tau);
            if along <= sweep.abs() {
                (d.magnitude() - radius).abs()
            } else {
                p.distance(piece.start()).min(p.distance(piece.end()))
            }
        }
        Piece2::Circle { centre, radius, .. } => ((p - centre).magnitude() - radius).abs(),
    }
}

/// The root mean square of distances.
fn rms(distances: impl Iterator<Item = f64>) -> f64 {
    let (sum, count) = distances.fold((0.0, 0.0), |(s, c), d| (d.mul_add(d, s), c + 1.0));
    if count > 0.0 {
        (sum / count).sqrt()
    } else {
        0.0
    }
}
