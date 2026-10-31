//! Tails peeled off a boundary ring that crosses itself.
//!
//! A face's boundary can run out along a path and come back along it: two
//! edges lying along each other within their tolerance, each fitted on its
//! own, so that in the chart the ring runs back over itself and crosses
//! its own way out. Over that stretch the face bounds nothing, and a
//! constrained triangulation, which takes the ring's segments as walls,
//! refuses the segments that cross and loses the boundary points the
//! neighbours across those edges keep.
//!
//! Each crossing of the ring cuts it into two loops. Where one of them is
//! a strip no wider than the widest tolerance of the face's edges, it is a
//! tail: it is taken out of the ring, which keeps the chord across its
//! mouth, and drawn on its own as a strip of triangles zipped across it,
//! at the points its edges put it, so the whole mesh still closes over it.
//! A ring left crossing itself, or one whose loops are both wide, is left
//! as it is.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Surface as _, SurfaceGeometry};
use ogeom_math::{Point, Point2, Transform};

/// A face's rings with their tails taken out, and the tails.
pub(crate) struct Peeled {
    /// The rings, each without its tails.
    pub(crate) rings: Vec<Vec<Point2>>,
    /// The tails taken out.
    pub(crate) strips: Vec<Strip>,
}

/// One tail taken out of a ring.
pub(crate) struct Strip {
    /// The chord across the tail's mouth, from and to as the ring kept
    /// walks it.
    pub(crate) mouth: (Point2, Point2),
    /// Triangles across the tail, each wound as the tail walks its own
    /// points: the mouth is walked the other way.
    pub(crate) triangles: Vec<[Point2; 3]>,
}

/// `rings` with every tail taken out (see the module documentation), or
/// `None` where no ring has one, or a ring crosses itself other than at
/// a tail. `reach` is the widest tolerance of the face's edges and
/// vertices. Widths are measured in space, on the surface placed as the
/// face places it.
///
/// # Errors
///
/// Where the surface cannot be evaluated at a ring point.
pub(crate) fn peel(
    rings: &[Vec<Point2>],
    surface: &SurfaceGeometry,
    placement: &Transform,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Peeled>> {
    let mut kept = Vec::with_capacity(rings.len());
    let mut strips = Vec::new();
    for ring in rings {
        let crossings = crossings(ring);
        if crossings.is_empty() {
            kept.push(ring.clone());
            continue;
        }
        let lifted = ring
            .iter()
            .map(|p| Ok(placement.apply(surface.point_at(p.x, p.y, tol)?)))
            .collect::<OgeomResult<Vec<Point>>>()?;
        let Some((ring_kept, ring_strips)) = peel_ring(ring, &lifted, &crossings, reach) else {
            return Ok(None);
        };
        kept.push(ring_kept);
        strips.extend(ring_strips);
    }
    if strips.is_empty() {
        return Ok(None);
    }
    Ok(Some(Peeled {
        rings: kept,
        strips,
    }))
}

/// One ring with its tails taken out, and the tails; `None` where a
/// crossing cuts off no tail, or what is left is no ring.
fn peel_ring(
    ring: &[Point2],
    lifted: &[Point],
    crossings: &[(usize, usize)],
    reach: f64,
) -> Option<(Vec<Point2>, Vec<Strip>)> {
    let n = ring.len();
    // Each crossing of segment `s` with segment `t` (from `s` to `s + 1`,
    // `t` to `t + 1`) leaves two loops: one from `s` on to `t + 1`, one
    // from `t` on to `s + 1`, each holding both crossing segments. The
    // tail is the one that is a strip, as `(first, last)` walked forward.
    let mut arcs: Vec<(usize, usize)> = Vec::with_capacity(crossings.len());
    for &(s, t) in crossings {
        let one = (s, (t + 1) % n);
        let other = (t, (s + 1) % n);
        let narrow =
            |(first, last): (usize, usize)| zip(&walk(lifted, first, last), reach).is_some();
        match (narrow(one), narrow(other)) {
            (true, false) => arcs.push(one),
            (false, true) => arcs.push(other),
            _ => return None,
        }
    }
    // Started at a point inside no tail, the walk meets every tail whole,
    // and tails that overlap are one.
    let inside = |k: usize, &(first, last): &(usize, usize)| {
        let along = (k + n - first) % n;
        along > 0 && along < (last + n - first) % n
    };
    let start = (0..n).find(|&k| !arcs.iter().any(|arc| inside(k, arc)))?;
    let mut spans: Vec<(usize, usize)> = arcs
        .iter()
        .map(|&(first, last)| {
            let from = (first + n - start) % n;
            let to = (last + n - start) % n;
            (from, if to < from { to + n } else { to })
        })
        .collect();
    spans.sort_unstable();
    let mut merged: Vec<(usize, usize)> = Vec::with_capacity(spans.len());
    for (from, to) in spans {
        match merged.last_mut() {
            Some(last) if from < last.1 => last.1 = last.1.max(to),
            _ => merged.push((from, to)),
        }
    }
    let at = |k: usize| (start + k) % n;
    let mut kept = Vec::with_capacity(n);
    let mut strips = Vec::with_capacity(merged.len());
    let mut next = 0;
    for &(from, to) in &merged {
        kept.extend((next..=from).map(|k| ring[at(k)]));
        let points: Vec<usize> = (from..=to).map(at).collect();
        let space: Vec<Point> = points.iter().map(|&k| lifted[k]).collect();
        let triangles = zip(&space, reach)?;
        strips.push(Strip {
            mouth: (ring[at(from)], ring[at(to)]),
            triangles: triangles
                .into_iter()
                .map(|t| t.map(|k| ring[points[k]]))
                .collect(),
        });
        next = to;
    }
    kept.extend((next..n).map(|k| ring[at(k)]));
    if merged.last().is_some_and(|&(_, to)| to >= n) {
        // The last tail ends past the walk's start, at a point already
        // kept as the first.
        kept.pop();
    }
    (kept.len() >= 3).then_some((kept, strips))
}

/// The ring's points from `first` on to `last`, wrapping.
fn walk(lifted: &[Point], first: usize, last: usize) -> Vec<Point> {
    let n = lifted.len();
    let count = (last + n - first) % n + 1;
    (0..count).map(|k| lifted[(first + k) % n]).collect()
}

/// Triangles zipped across a loop from its two ends (its mouth) in,
/// each the next point of one side or the other joined to the last of
/// the other, whichever joins the nearer two, wound as the loop walks its
/// points; `None` where a triangle is wider than `reach`, the loop no
/// strip.
fn zip(points: &[Point], reach: f64) -> Option<Vec<[usize; 3]>> {
    let mut triangles = Vec::with_capacity(points.len().saturating_sub(2));
    let (mut left, mut right) = (0, points.len().checked_sub(1)?);
    while right > left + 1 {
        let triangle = if points[left + 1].distance(points[right])
            <= points[left].distance(points[right - 1])
        {
            left += 1;
            [left - 1, left, right]
        } else {
            right -= 1;
            [left, right, right + 1]
        };
        let [a, b, c] = triangle.map(|k| points[k]);
        let longest = a.distance(b).max(b.distance(c)).max(c.distance(a));
        let twice_area = (b - a).cross(c - a).magnitude();
        if longest > 0.0 && twice_area / longest > reach {
            return None;
        }
        triangles.push(triangle);
    }
    Some(triangles)
}

/// Every pair of the ring's segments that cross properly, as `(s, t)`
/// with `s < t`, segment `k` running from point `k` to the next.
fn crossings(ring: &[Point2]) -> Vec<(usize, usize)> {
    let n = ring.len();
    if n < 4 {
        return Vec::new();
    }
    let segment = |k: usize| (ring[k], ring[(k + 1) % n]);
    let mut order: Vec<(f64, f64, usize)> = (0..n)
        .map(|k| {
            let (a, b) = segment(k);
            (a.x.min(b.x), a.x.max(b.x), k)
        })
        .collect();
    order.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    let mut found = Vec::new();
    for (i, &(_, high, s)) in order.iter().enumerate() {
        let (a, b) = segment(s);
        for &(low, _, t) in &order[i + 1..] {
            if low > high {
                break;
            }
            let (c, d) = segment(t);
            if a.y.max(b.y) < c.y.min(d.y) || c.y.max(d.y) < a.y.min(b.y) {
                continue;
            }
            if crate::triangulate::segments_cross(a, b, c, d) {
                found.push((s.min(t), s.max(t)));
            }
        }
    }
    found.sort_unstable();
    found
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    /// A ring in the plane z = 0, lifted as it is.
    fn flat(ring: &[(f64, f64)]) -> (Vec<Point2>, Vec<Point>) {
        let chart: Vec<Point2> = ring.iter().map(|&(x, y)| Point2::new(x, y)).collect();
        let space = chart.iter().map(|p| Point::new(p.x, p.y, 0.0)).collect();
        (chart, space)
    }

    /// A unit square, its top side run on along itself to a tip at 0.2
    /// and its way back a hair either side of it until it parts at 0.62
    /// and heads down to its start; with `parting` the last point before
    /// the way back leaves the top.
    fn tailed(parting: &[(f64, f64)]) -> (Vec<Point2>, Vec<Point>) {
        let mut ring = vec![
            (0.5, 0.0),
            (1.0, 0.0),
            (1.0, 1.0),
            (0.7, 1.0),
            (0.4, 1.0),
            (0.2, 1.0),
            (0.35, 1.0 - 1e-4),
            (0.55, 1.0 + 1e-4),
        ];
        ring.extend_from_slice(parting);
        ring.push((0.8, 0.5));
        flat(&ring)
    }

    /// The way back crosses the top twice and parts from it within the
    /// reach: the tail out to the tip comes off, and the square keeps the
    /// chord across its mouth.
    #[test]
    fn a_tail_crossing_its_way_out_comes_off() {
        let (ring, lifted) = tailed(&[(0.62, 1.0 - 1e-4)]);
        let crossings = crossings(&ring);
        assert_eq!(crossings, vec![(3, 6), (3, 7)]);
        let (kept, strips) = peel_ring(&ring, &lifted, &crossings, 1e-3).unwrap();
        assert_eq!(kept, [0, 1, 2, 3, 8, 9].map(|k| ring[k]));
        assert_eq!(strips.len(), 1);
        assert_eq!(strips[0].mouth, (ring[3], ring[8]));
        assert_eq!(strips[0].triangles.len(), 4);
    }

    /// The way back leaves the top across one long segment: cut at its
    /// far end, the tail would take a wide piece of the square with it,
    /// so nothing comes off.
    #[test]
    fn a_tail_leaving_across_a_long_segment_stays() {
        let (ring, lifted) = tailed(&[]);
        let crossings = crossings(&ring);
        assert_eq!(crossings, vec![(3, 6), (3, 7)]);
        assert!(peel_ring(&ring, &lifted, &crossings, 1e-3).is_none());
    }

    /// A ring crossing itself in a figure of eight, both loops wide: no
    /// tail, and nothing comes off.
    #[test]
    fn a_figure_of_eight_is_no_tail() {
        let (ring, lifted) = flat(&[(0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 1.0)]);
        let crossings = crossings(&ring);
        assert_eq!(crossings, vec![(0, 2)]);
        assert!(peel_ring(&ring, &lifted, &crossings, 1e-3).is_none());
    }
}
