//! Holes that enclose no area the face's tolerance can resolve.
//!
//! A face lies within its tolerance `t` of where it says, and so does every
//! point of its boundary. A hole every point of which stands within `t` of
//! the hole's own boundary is a strip that boundary collapses onto its
//! middle without moving any point farther than `t`: a hole of nothing, by
//! the face's own admission. Kept, it is a boundary every boolean through
//! the face cuts, and a slit the mesher has to draw.
//!
//! The area such a hole encloses is at most `t` times its boundary's
//! length `L`: every point of it lies within `t` of the boundary, and the
//! band reaching `t` in from a boundary of length `L` covers no more. So a
//! hole enclosing more than `t·L` is kept without further ado. That bound
//! is necessary, not sufficient: a round hole of radius well past `t` with
//! a long thin slit off it encloses less than `t·L`, and is a hole. What
//! decides is the hole's widest point: the region is covered by cells,
//! each found either outside the hole or within `t` of its boundary
//! throughout by the distance from its centre, which no point of the cell
//! exceeds by more than its half diagonal. A cell whose centre stands
//! inside farther than that keeps the hole.
//!
//! This is measured in the surface's chart drawn to scale, which is the
//! surface itself only where the chart has the same scale everywhere and
//! its two directions meet square: a plane, a cylinder (its girth scaled by
//! its radius). Distances along the surface are no shorter than through
//! space, so a point within `t` along the surface is within `t`. A hole on
//! any other surface is kept. The boundary is drawn as a polyline, each
//! pcurve divided until its middle points stand off their chords by no
//! more than a quarter of `t`, and the reach asked of the cells is `t` less
//! that.
//!
//! Only a hole no other face walks an edge of goes: a hole whose edges a
//! neighbour shares bounds that neighbour too. The face must check clean
//! without it.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve2d as _, Surface as _, SurfaceGeometry, SurfaceKind};
use ogeom_math::{Point2, Vector2};
use ogeom_topo::{EdgeRepr, Model, Orientation, Shape, ShapeType, explore_unique};

use crate::Reshape;

/// The most cells a hole is covered by before it is kept as unproven: a
/// cost, not a measure of the hole.
const CELLS: usize = 1 << 16;

/// The most pieces a pcurve is divided into before the hole is kept as
/// undrawn: a cost, not a measure of the hole.
const PIECES: usize = 1 << 12;

/// Stage the removal of every hole of a face of `shape` that encloses no
/// area its face's tolerance can resolve (see the module documentation).
/// Returns how many holes go.
pub(crate) fn stage_small_wires(
    model: &mut Model,
    shape: &Shape,
    reshape: &mut Reshape,
    tol: Tolerances,
) -> OgeomResult<usize> {
    let uses = ogeom_algo::edge_uses(model, shape)?;
    let mut removed = 0;
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()).cloned() else {
            continue;
        };
        if !face.location().is_identity() || !data.location.is_identity() {
            continue;
        }
        let Some(surface) = model.geometry().surface(data.surface) else {
            continue;
        };
        let Some(scale) = chart_scale(surface, tol)? else {
            continue;
        };
        let reach = data.tolerance.get();
        let stored = face.oriented(Orientation::Forward);
        let wires = model.ordered_children_of(&stored)?;
        let mut gone = Vec::new();
        for wire in wires.iter().skip(1) {
            let edges = model.ordered_children_of(&wire.oriented(Orientation::Forward))?;
            if edges.iter().any(|e| uses.get(&e.node()) != Some(&1)) {
                continue;
            }
            let Some((ring, sag)) = drawn(model, data.surface, surface, &edges, scale, reach, tol)?
            else {
                continue;
            };
            if within_reach_of_boundary(&ring, reach - sag) {
                gone.push(wire.clone());
            }
        }
        if gone.is_empty() {
            continue;
        }
        let mut local = Reshape::new();
        for wire in &gone {
            local.remove(&Shape::of(wire.node()));
        }
        let after = local.apply(model, &face)?.shape;
        if !ogeom_algo::check(model, &after, tol)?.is_valid() {
            continue;
        }
        for wire in &gone {
            reshape.remove(&Shape::of(wire.node()));
        }
        removed += gone.len();
    }
    Ok(removed)
}

/// How far a step of one in each chart parameter carries a point along
/// `surface`, where that is the same everywhere and the two directions
/// meet square: a plane, a cylinder. `None` on any other surface.
fn chart_scale(surface: &SurfaceGeometry, tol: Tolerances) -> OgeomResult<Option<Vector2>> {
    if !matches!(surface.kind(), SurfaceKind::Plane | SurfaceKind::Cylinder) {
        return Ok(None);
    }
    let ((u0, u1), (v0, v1)) = surface.domain();
    let (u, v) = (mid(u0, u1), mid(v0, v1));
    let (du, dv) = surface.d1_at(u, v, tol)?;
    let (su, sv) = (du.magnitude(), dv.magnitude());
    if su <= 0.0 || sv <= 0.0 || du.dot(dv).abs() > tol.angular() * su * sv {
        return Ok(None);
    }
    Ok(Some(Vector2::new(su, sv)))
}

/// The middle of a range, or 0 where it is unbounded.
fn mid(a: f64, b: f64) -> f64 {
    let m = f64::midpoint(a, b);
    if m.is_finite() { m } else { 0.0 }
}

/// The hole walked by `edges` drawn to scale as a closed polyline, with
/// the widest its pieces' middles stand off their chords; `None` where an
/// edge has no plain pcurve on the surface, the ring does not close in one
/// copy of the chart, it encloses more than `reach` times its length, or
/// it cannot be drawn within a quarter of `reach`.
fn drawn(
    model: &Model,
    surface_id: ogeom_topo::SurfaceId,
    surface: &SurfaceGeometry,
    edges: &[Shape],
    scale: Vector2,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Vec<Point2>, f64)>> {
    let ((u0, u1), _) = surface.domain();
    let period = (surface.is_periodic_u() || surface.is_closed_u(tol)).then_some(u1 - u0);
    let carried = |from: Point2, to: Point2| match period {
        Some(p) if p > 0.0 && p.is_finite() => ((to.x - from.x) / p).round() * p,
        _ => 0.0,
    };
    let to_scale = |p: Point2| Point2::new(p.x * scale.x, p.y * scale.y);
    let mut ring: Vec<Point2> = Vec::new();
    let mut sag: f64 = 0.0;
    for edge in edges {
        let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
            return Ok(None);
        };
        if data.degenerate {
            return Ok(None);
        }
        let Some(EdgeRepr::PCurve { curve, range, .. }) =
            data.pcurve_for(surface_id, &ogeom_topo::Location::identity())
        else {
            return Ok(None);
        };
        let Some(pcurve) = model.geometry().pcurve(*curve) else {
            return Ok(None);
        };
        let reversed = edge.orientation() == Orientation::Reversed;
        let at = |s: f64| -> OgeomResult<Point2> {
            let s = if reversed { 1.0 - s } else { s };
            pcurve.point_at(range.0 + (range.1 - range.0) * s, tol)
        };
        let mut pieces = 1;
        let (points, off) = loop {
            #[allow(clippy::cast_precision_loss, reason = "a count of pieces")]
            let step = 1.0 / pieces as f64;
            let mut points = Vec::with_capacity(pieces + 1);
            let mut off: f64 = 0.0;
            for k in 0..=pieces {
                #[allow(clippy::cast_precision_loss, reason = "a count of pieces")]
                let s = k as f64 * step;
                points.push(to_scale(at(s)?));
                if k > 0 {
                    let middle = to_scale(at(s - 0.5 * step)?);
                    let chord = points[k - 1].midpoint(points[k]);
                    off = off.max(middle.distance(chord));
                }
            }
            if off <= 0.25 * reach {
                break (points, off);
            }
            pieces *= 2;
            if pieces > PIECES {
                return Ok(None);
            }
        };
        sag = sag.max(off);
        let by = match ring.last() {
            Some(&last) => carried(points[0], last) * scale.x,
            None => 0.0,
        };
        ring.extend(points.into_iter().map(|p| p + Vector2::new(by, 0.0)));
    }
    let (Some(&first), Some(&last)) = (ring.first(), ring.last()) else {
        return Ok(None);
    };
    if carried(first, last) != 0.0 {
        return Ok(None);
    }
    let length: f64 = ring.windows(2).map(|w| w[0].distance(w[1])).sum();
    if crate::divide::signed_area(&ring).abs() > reach * length {
        return Ok(None);
    }
    Ok(Some((ring, sag)))
}

/// Whether every point inside the closed polyline `ring` stands within
/// `reach` of it, proved cell by cell (see the module documentation);
/// `false` where it is not, or where the proof runs past its cost.
fn within_reach_of_boundary(ring: &[Point2], reach: f64) -> bool {
    if reach <= 0.0 || ring.len() < 3 {
        return false;
    }
    let (mut low, mut high) = (ring[0], ring[0]);
    for p in ring {
        low = Point2::new(low.x.min(p.x), low.y.min(p.y));
        high = Point2::new(high.x.max(p.x), high.y.max(p.y));
    }
    let half = 0.5 * (high.x - low.x).max(high.y - low.y);
    let mut cells = vec![(low.midpoint(high), half)];
    let mut asked = 0;
    while let Some((centre, half)) = cells.pop() {
        asked += 1;
        if asked > CELLS {
            return false;
        }
        let near = nearest(ring, centre);
        let corner = half * std::f64::consts::SQRT_2;
        if near + corner <= reach {
            continue;
        }
        let inside = winds(ring, centre);
        if !inside && near >= corner {
            continue;
        }
        if inside && near > reach {
            return false;
        }
        let quarter = 0.5 * half;
        for (dx, dy) in [(-1.0, -1.0), (1.0, -1.0), (-1.0, 1.0), (1.0, 1.0)] {
            cells.push((centre + Vector2::new(dx * quarter, dy * quarter), quarter));
        }
    }
    true
}

/// The distance from `p` to the nearest point of the closed polyline.
fn nearest(ring: &[Point2], p: Point2) -> f64 {
    let mut best = f64::INFINITY;
    for k in 0..ring.len() {
        let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
        let ab = b - a;
        let length = ab.dot(ab);
        let s = if length > 0.0 {
            ((p - a).dot(ab) / length).clamp(0.0, 1.0)
        } else {
            0.0
        };
        best = best.min(p.distance(a + ab * s));
    }
    best
}

/// Whether the closed polyline winds round `p`.
fn winds(ring: &[Point2], p: Point2) -> bool {
    let mut inside = false;
    for k in 0..ring.len() {
        let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
        if (a.y > p.y) != (b.y > p.y) {
            let x = a.x + (p.y - a.y) / (b.y - a.y) * (b.x - a.x);
            if x > p.x {
                inside = !inside;
            }
        }
    }
    inside
}
