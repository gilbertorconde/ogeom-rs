//! Where a seam runs in a face's chart, and what rings measure there: the
//! seam between two rims or from a cap's rim to its pole, rings moved onto
//! the branch the face's pcurves are read on, how many times a ring goes
//! round, the area it encloses, and distances to segments and lines.

use ogeom_core::Tolerances;
use ogeom_math::{Point, Point2};
use ogeom_topo::Shape;

use super::segment::{chart, periodic, unwrapped};
use super::weld::{Half, from_to};
use super::{Curved, Stop};
use crate::recognize::Canonical;

/// Where a rim's entry starts: the entry, its vertex, and where that
/// stands.
pub(super) type RimStart = ((usize, bool), Shape, Point);

/// A seam's two rim vertices, by their places in each rim, and the ends of
/// its line in the chart.
type SeamChoice = (usize, usize, (f64, f64), (f64, f64));

/// The tube angle of a torus's seam round its axis: half a turn from its
/// chart's centre.
pub(super) fn torus_seam_v(curved: &Curved) -> f64 {
    curved.centre.1 - core::f64::consts::PI
}

/// A face's holes in its chart: each a polygon of its mesh vertices,
/// carried round continuously.
pub(super) fn hole_polygons(
    shape: &Canonical,
    holes: &[&[Half]],
    triangles: &[[u32; 3]],
    points: &[Point],
    tol: Tolerances,
) -> Vec<Vec<(f64, f64)>> {
    holes
        .iter()
        .filter_map(|ring| {
            let mut out: Vec<(f64, f64)> = Vec::new();
            for &h in *ring {
                let (a, _) = from_to(triangles, h);
                let (u, v) = chart(shape, points[a as usize], tol)?;
                // Unwrapped along the ring both ways the chart closes on
                // itself: a hole across a torus's outer equator, where the
                // tube's angle starts, runs on past a whole turn rather
                // than jumping back to nought, so a seam there is seen to
                // cross it.
                let (u, v) = match out.last() {
                    Some(&(last_u, last_v)) => (
                        last_u + ogeom_math::elementary::wrap_signed_angle(u - last_u),
                        if matches!(shape, Canonical::Torus(_)) {
                            last_v + ogeom_math::elementary::wrap_signed_angle(v - last_v)
                        } else {
                            v
                        },
                    ),
                    None => (u, v),
                };
                out.push((u, v));
            }
            Some(out)
        })
        .collect()
}

/// The seam for a face round its axis (or round a torus's tube, with
/// `round_tube`): of the pairs of a vertex on one rim and a vertex on the
/// other, the one turning least between them whose straight chart line
/// crosses no hole, a whole turn either way included. Its indices and the
/// line's ends in the chart, the angle it goes round first (as `holes`
/// are given).
pub(super) fn choose_seam(
    curved: &Curved,
    round_tube: bool,
    from: &[Point],
    to: &[Point],
    holes: &[Vec<(f64, f64)>],
    tol: Tolerances,
) -> Option<SeamChoice> {
    let tau = core::f64::consts::TAU;
    // A line through a hole's vertex counts as crossing the hole, though
    // no side of it crosses the line strictly; so does one passing within
    // a millionth of the holes' mean step of a vertex, a hair from it.
    let steps: Vec<f64> = holes
        .iter()
        .flat_map(|ring| {
            (0..ring.len()).map(move |i| {
                let (p, q) = (ring[i], ring[(i + 1) % ring.len()]);
                (q.0 - p.0).hypot(q.1 - p.1)
            })
        })
        .collect();
    #[allow(clippy::cast_precision_loss, reason = "ring lengths are small")]
    let margin = 1e-6 * steps.iter().sum::<f64>() / steps.len().max(1) as f64;
    let crosses = |a: (f64, f64), b: (f64, f64)| {
        holes.iter().any(|ring| {
            [-tau, 0.0, tau].iter().any(|shift| {
                (0..ring.len()).any(|i| {
                    let (p, q) = (ring[i], ring[(i + 1) % ring.len()]);
                    segments_cross(a, b, (p.0 + shift, p.1), (q.0 + shift, q.1))
                }) || ring
                    .iter()
                    .any(|h| chart_distance_to_segment((h.0 + shift, h.1), a, b) < margin)
            })
        })
    };
    let mut best: Option<(f64, SeamChoice)> = None;
    let at = |p: Point| unwrapped(curved, p, tol).map(|c| swapped(c, round_tube));
    for (i, pa) in from.iter().enumerate() {
        let Some((ua, va)) = at(*pa) else {
            continue;
        };
        for (j, pb) in to.iter().enumerate() {
            let Some((ub, vb)) = at(*pb) else {
                continue;
            };
            let turn = ogeom_math::elementary::wrap_signed_angle(ub - ua);
            let (a, b) = ((ua, va), (ua + turn, vb));
            if crosses(a, b) {
                continue;
            }
            let score = turn.abs() * 1e3 + pa.distance(*pb);
            if best.is_none_or(|held| score < held.0) {
                best = Some((score, (i, j, a, b)));
            }
        }
    }
    best.map(|(_, choice)| choice)
}

/// A sphere's cap with holes: the one ring that goes round its axis (the
/// rim), which way it goes round, and the rings that do not (the holes).
pub(super) fn cap_rings(
    curved: &Curved,
    rings: &[Vec<Half>],
    triangles: &[[u32; 3]],
    points: &[Point],
    tol: Tolerances,
) -> Option<(usize, i32, Vec<usize>)> {
    let windings: Vec<i32> = rings
        .iter()
        .map(|ring| winding(&curved.shape, ring, triangles, points, tol))
        .collect::<Option<_>>()?;
    let rims: Vec<usize> = (0..rings.len())
        .filter(|&k| windings[k].abs() == 1)
        .collect();
    let [rim] = rims[..] else {
        return None;
    };
    let holes: Vec<usize> = (0..rings.len()).filter(|&k| windings[k] == 0).collect();
    (!holes.is_empty() && holes.len() + 1 == rings.len()).then_some((rim, windings[rim], holes))
}

/// The seam of a sphere's cap with holes, by the index of the rim vertex
/// it starts from and its line's ends in the chart, on the branch the
/// face's pcurves are read on.
type CapSeamChoice = (usize, (f64, f64), (f64, f64));

/// A straight segment in a chart, by its ends.
type ChartSegment = ((f64, f64), (f64, f64));

/// The seam of a sphere's cap with holes: from one of the rim's `starts`
/// to the pole, straight in the chart, crossing no hole and meeting the
/// rim only where it starts. A meridian where one is clear, the one
/// standing farthest round from the holes; otherwise the line that turns
/// least about the axis on its way up.
pub(super) fn cap_seam(
    curved: &Curved,
    rim: &[Half],
    holes: &[&[Half]],
    starts: &[Point],
    triangles: &[[u32; 3]],
    points: &[Point],
    tol: Tolerances,
) -> Option<CapSeamChoice> {
    use core::f64::consts::{FRAC_PI_2, TAU};
    const TURNS: i32 = 32;
    let holes = centred_rings(
        curved,
        hole_polygons(&curved.shape, holes, triangles, points, tol),
    );
    let rim = centred_rings(
        curved,
        hole_polygons(&curved.shape, &[rim], triangles, points, tol),
    )
    .pop()?;
    let closed = |ring: &[(f64, f64)]| -> Vec<ChartSegment> {
        (0..ring.len())
            .map(|k| {
                let p = ring[k];
                let q = if k + 1 < ring.len() {
                    ring[k + 1]
                } else {
                    // A ring closes where it began, the rim a whole turn on.
                    let first = ring[0];
                    (first.0 + ((p.0 - first.0) / TAU).round() * TAU, first.1)
                };
                (p, q)
            })
            .collect()
    };
    let rim_segments = closed(&rim);
    let hole_segments: Vec<_> = holes.iter().flat_map(|ring| closed(ring)).collect();
    let crosses = |segments: &[ChartSegment], a: (f64, f64), b: (f64, f64)| {
        segments.iter().any(|&(p, q)| {
            [-2.0 * TAU, -TAU, 0.0, TAU, 2.0 * TAU]
                .iter()
                .any(|s| segments_cross(a, b, (p.0 + s, p.1), (q.0 + s, q.1)))
        })
    };
    // A seam running along a hole's side, or within half its mean step
    // of a hole's vertex, counts as crossing it: it would leave the face
    // a sliver there.
    #[allow(clippy::cast_precision_loss, reason = "ring lengths are small")]
    let margin = 0.5
        * hole_segments
            .iter()
            .map(|(p, q)| (q.0 - p.0).hypot(q.1 - p.1))
            .sum::<f64>()
        / hole_segments.len().max(1) as f64;
    let grazes = |a: (f64, f64), b: (f64, f64)| {
        holes.iter().flatten().any(|h| {
            [-2.0 * TAU, -TAU, 0.0, TAU, 2.0 * TAU].iter().any(|s| {
                let (x, y) = (h.0 + s - a.0, h.1 - a.1);
                let (dx, dy) = (b.0 - a.0, b.1 - a.1);
                let f = ((x * dx + y * dy) / dx.mul_add(dx, dy * dy)).clamp(0.0, 1.0);
                (x - f * dx).hypot(y - f * dy) < margin
            })
        })
    };
    let at: Vec<Option<(f64, f64)>> = starts.iter().map(|p| unwrapped(curved, *p, tol)).collect();
    let clearance = |u: f64| {
        holes
            .iter()
            .flatten()
            .map(|h| ogeom_math::elementary::wrap_signed_angle(h.0 - u).abs())
            .fold(f64::INFINITY, f64::min)
    };
    let mut offsets = vec![0];
    for k in 1..=TURNS / 2 {
        offsets.extend([k, -k]);
    }
    for k in offsets {
        let turn = TAU * f64::from(k) / f64::from(TURNS);
        let mut best: Option<(f64, CapSeamChoice)> = None;
        for (i, a) in at.iter().enumerate() {
            let Some(a) = *a else {
                continue;
            };
            if a.1 >= FRAC_PI_2 - 1e-9 {
                continue;
            }
            let b = (a.0 + turn, FRAC_PI_2);
            // Just above the rim, so the rim's own sides at the start
            // vertex do not count.
            let lifted = (a.0 + (b.0 - a.0) * 1e-3, a.1 + (b.1 - a.1) * 1e-3);
            if crosses(&hole_segments, a, b) || grazes(a, b) || crosses(&rim_segments, lifted, b) {
                continue;
            }
            let score = clearance(a.0);
            if best.is_none_or(|held| score > held.0) {
                best = Some((score, (i, a, b)));
            }
        }
        if let Some((_, choice)) = best {
            return Some(choice);
        }
    }
    None
}

/// The straight pieces of a [`Thread`]'s chain in its working chart: from
/// the start to the first hole, from each hole to the next, and from the
/// last back to the start a whole turn on.
///
/// [`Thread`]: super::Thread
pub(super) fn thread_pieces(start: (f64, f64), stops: &[Stop]) -> Vec<((f64, f64), (f64, f64))> {
    let mut pieces = Vec::with_capacity(stops.len() + 1);
    let mut from = start;
    for stop in stops {
        pieces.push((from, stop.enter_at));
        from = stop.leave_at;
    }
    pieces.push((from, (start.0 + core::f64::consts::TAU, start.1)));
    pieces
}

/// A chart point with its coordinates exchanged where `yes`: a face round a
/// torus's tube is read with the tube's angle first, as one round the axis
/// reads the axis's.
pub(super) fn swapped(p: (f64, f64), yes: bool) -> (f64, f64) {
    if yes { (p.1, p.0) } else { p }
}

/// Rings in a face's chart, each moved by whole turns so that it starts
/// on the branch the face's pcurves are read on.
pub(super) fn centred_rings(curved: &Curved, rings: Vec<Vec<(f64, f64)>>) -> Vec<Vec<(f64, f64)>> {
    let (pu, pv) = periodic(&curved.shape);
    let tau = core::f64::consts::TAU;
    let shift = |x: f64, c: f64, wraps: bool| {
        if wraps {
            ((c - x) / tau).round() * tau
        } else {
            0.0
        }
    };
    rings
        .into_iter()
        .map(|ring| {
            let Some(&(u, v)) = ring.first() else {
                return ring;
            };
            let (du, dv) = (shift(u, curved.centre.0, pu), shift(v, curved.centre.1, pv));
            ring.into_iter().map(|(u, v)| (u + du, v + dv)).collect()
        })
        .collect()
}

/// [`swapped`] for every point of some rings.
pub(super) fn swapped_rings(rings: Vec<Vec<(f64, f64)>>, yes: bool) -> Vec<Vec<(f64, f64)>> {
    if !yes {
        return rings;
    }
    rings
        .into_iter()
        .map(|ring| ring.into_iter().map(|p| swapped(p, true)).collect())
        .collect()
}

/// Whether a face goes round a torus's tube but not round its axis: its
/// rims, and the seam between them, are read with the tube's angle first.
pub(super) fn round_the_tube(curved: &Curved) -> bool {
    matches!(curved.shape, Canonical::Torus(_)) && curved.wraps_v && !curved.wraps
}

/// How many times a ring goes round the way its face does: round the
/// tube for a face [`round_the_tube`], round the axis otherwise.
pub(super) fn turns_along(
    curved: &Curved,
    ring: &[Half],
    triangles: &[[u32; 3]],
    points: &[Point],
    tol: Tolerances,
) -> Option<i32> {
    if round_the_tube(curved) {
        windings(&curved.shape, ring, triangles, points, tol).map(|w| w.1)
    } else {
        winding(&curved.shape, ring, triangles, points, tol)
    }
}

/// The middle of the widest stretch of angle no hole covers.
pub(super) fn free_angle(holes: &[Vec<(f64, f64)>]) -> Option<f64> {
    let tau = core::f64::consts::TAU;
    let mut angles: Vec<f64> = holes
        .iter()
        .flatten()
        .map(|(u, _)| u.rem_euclid(tau))
        .collect();
    if angles.is_empty() {
        return None;
    }
    angles.sort_by(f64::total_cmp);
    let mut best = (
        angles[0] + tau - angles[angles.len() - 1],
        angles[angles.len() - 1],
    );
    for pair in angles.windows(2) {
        if pair[1] - pair[0] > best.0 {
            best = (pair[1] - pair[0], pair[0]);
        }
    }
    Some(best.1 + best.0 / 2.0)
}

/// The middle of the widest gap between `bounds` (angles) that holds at
/// least one of `inside`, in `[0, 2pi)`.
pub(super) fn widest_gap_holding(bounds: &mut [f64], inside: &[f64]) -> Option<f64> {
    let tau = core::f64::consts::TAU;
    for b in bounds.iter_mut() {
        *b = b.rem_euclid(tau);
    }
    bounds.sort_by(f64::total_cmp);
    let mut best: Option<(f64, f64)> = None;
    for k in 0..bounds.len() {
        let from = bounds[k];
        let to = if k + 1 < bounds.len() {
            bounds[k + 1]
        } else {
            bounds[0] + tau
        };
        let holds = inside
            .iter()
            .any(|&x| (x - from).rem_euclid(tau) < to - from && (x - from).rem_euclid(tau) > 0.0);
        if holds && best.is_none_or(|(width, _)| to - from > width) {
            best = Some((to - from, from));
        }
    }
    best.map(|(width, from)| ogeom_math::elementary::wrap_angle(from + width / 2.0))
}

/// Whether two chart segments cross, each at a point strictly inside both.
pub(super) fn segments_cross(a: (f64, f64), b: (f64, f64), p: (f64, f64), q: (f64, f64)) -> bool {
    let side = |o: (f64, f64), x: (f64, f64), y: (f64, f64)| {
        (x.0 - o.0).mul_add(y.1 - o.1, -((x.1 - o.1) * (y.0 - o.0)))
    };
    let (d1, d2) = (side(p, q, a), side(p, q, b));
    let (d3, d4) = (side(a, b, p), side(a, b, q));
    d1 * d2 < 0.0 && d3 * d4 < 0.0
}

/// How many times a ring of half-edges goes round a surface's angle, with
/// the sign of its sense; `None` where a vertex has no chart position.
fn winding(
    shape: &Canonical,
    ring: &[Half],
    triangles: &[[u32; 3]],
    points: &[Point],
    tol: Tolerances,
) -> Option<i32> {
    let mut turned = 0.0;
    for &h in ring {
        let (a, b) = from_to(triangles, h);
        let (ua, _) = chart(shape, points[a as usize], tol)?;
        let (ub, _) = chart(shape, points[b as usize], tol)?;
        turned += ogeom_math::elementary::wrap_signed_angle(ub - ua);
    }
    #[allow(clippy::cast_possible_truncation, reason = "a handful of turns")]
    Some((turned / core::f64::consts::TAU).round() as i32)
}

/// How many times a ring goes round each of a surface's chart directions,
/// as [`winding`] counts the first.
pub(super) fn windings(
    shape: &Canonical,
    ring: &[Half],
    triangles: &[[u32; 3]],
    points: &[Point],
    tol: Tolerances,
) -> Option<(i32, i32)> {
    let (_, wraps_v) = periodic(shape);
    let mut turned = (0.0, 0.0);
    for &h in ring {
        let (a, b) = from_to(triangles, h);
        let (ua, va) = chart(shape, points[a as usize], tol)?;
        let (ub, vb) = chart(shape, points[b as usize], tol)?;
        turned.0 += ogeom_math::elementary::wrap_signed_angle(ub - ua);
        if wraps_v {
            turned.1 += ogeom_math::elementary::wrap_signed_angle(vb - va);
        }
    }
    let tau = core::f64::consts::TAU;
    #[allow(clippy::cast_possible_truncation, reason = "a handful of turns")]
    Some((
        (turned.0 / tau).round() as i32,
        (turned.1 / tau).round() as i32,
    ))
}

/// How many times a loop of mesh vertices goes round each of a surface's
/// chart directions, as [`windings`] counts a ring of half-edges.
pub(super) fn loop_turns(
    shape: &Canonical,
    ring: &[u32],
    points: &[Point],
    tol: Tolerances,
) -> Option<(i32, i32)> {
    let (_, wraps_v) = periodic(shape);
    let mut turned = (0.0, 0.0);
    for (k, &a) in ring.iter().enumerate() {
        let b = ring[(k + 1) % ring.len()];
        let (ua, va) = chart(shape, points[a as usize], tol)?;
        let (ub, vb) = chart(shape, points[b as usize], tol)?;
        turned.0 += ogeom_math::elementary::wrap_signed_angle(ub - ua);
        if wraps_v {
            turned.1 += ogeom_math::elementary::wrap_signed_angle(vb - va);
        }
    }
    let tau = core::f64::consts::TAU;
    #[allow(clippy::cast_possible_truncation, reason = "a handful of turns")]
    Some((
        (turned.0 / tau).round() as i32,
        (turned.1 / tau).round() as i32,
    ))
}

/// Twice the area a ring of half-edges encloses in a chart, signed.
pub(super) fn ring_area(
    ring: &[Half],
    triangles: &[[u32; 3]],
    chart: impl Fn(Point) -> Option<Point2>,
    points: &[Point],
) -> f64 {
    let mut area = 0.0;
    for &h in ring {
        let (a, b) = from_to(triangles, h);
        if let (Some(a), Some(b)) = (chart(points[a as usize]), chart(points[b as usize])) {
            area += a.x * b.y - b.x * a.y;
        }
    }
    area
}

/// How far a chart point stands from a chart segment.
fn chart_distance_to_segment(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (dx, dy) = (b.0 - a.0, b.1 - a.1);
    let (x, y) = (p.0 - a.0, p.1 - a.1);
    let length = dx.mul_add(dx, dy * dy);
    let f = if length > 0.0 {
        (x.mul_add(dx, y * dy) / length).clamp(0.0, 1.0)
    } else {
        0.0
    };
    f.mul_add(-dx, x).hypot(f.mul_add(-dy, y))
}

/// The distance from `p` to the segment from `a` to `b`.
pub(super) fn distance_to_segment(p: Point, a: Point, b: Point) -> f64 {
    let along = b - a;
    let length = along.dot(along);
    if length <= 0.0 {
        return p.distance(a);
    }
    let t = ((p - a).dot(along) / length).clamp(0.0, 1.0);
    p.distance(a + along * t)
}

pub(super) fn distance_to_line(p: Point, a: Point, b: Point) -> f64 {
    let d = b - a;
    let m = d.magnitude();
    if m == 0.0 {
        return p.distance(a);
    }
    (p - a).cross(d).magnitude() / m
}
