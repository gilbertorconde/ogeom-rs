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
//! This is measured in the surface's chart drawn to scale. Over a box of
//! the chart, a step of `du` carries a point along the surface no farther
//! than `|du|` times the most `∂S/∂u` reaches there, and likewise in `v`;
//! drawn with each direction stretched by that bound, the chart is a map
//! on which a straight run of length `l` within the box moves the point it
//! carries through space by at most `l` where the two directions meet
//! square everywhere (every elementary surface), and by at most `√2·l`
//! elsewhere (`|a| + |b| ≤ √2·|(a, b)|`). Each family's bound over the box
//! is its own:
//!
//! - a plane, a cylinder: the same everywhere, read off the surface;
//! - a cone: `|r(v)|` round, at the box's nearer or farther end, and the
//!   slant `1/cos α` along;
//! - a sphere: `R·|cos v|` round, at its widest over the box, and `R`
//!   along;
//! - a torus: `|R + r·cos v|` round, at its widest over the box, and `r`
//!   along;
//! - a B-spline patch, polynomial or rational: from its control net. On
//!   the stretch where the basis functions of control points `P_i` (weights
//!   `w_i`) are not zero, `S' = Σ N_i' w_i (P_i − S) / W`, so
//!   `|S'| ≤ max_i p·(w_{i+1}·|P_{i+1} − P_i| + |w_{i+1} − w_i|·D)/(t_{i+p+1} − t_{i+1}) / w_min`,
//!   with `D` the extent of those control points (each and `S` lie among
//!   them) and `w_min` their lightest weight, which `W` never falls below;
//! - a revolution of a B-spline: its points' distance from the axis round
//!   (the curve lies among its control points, and that distance is
//!   convex), the curve's own bound along;
//! - an extrusion of a B-spline: the curve's bound across, one along;
//! - a trimmed surface: its basis's.
//!
//! Any other surface keeps its holes. The box is where the hole may reach:
//! the drawn boundary's, widened by how far the boundary may stand off
//! its drawing. It depends on the stretch it is drawn with, so the hole is
//! drawn again with the bound over its box until the bound over the box
//! drawn is no wider than the one it was drawn with. Every cell is clipped
//! to the box, so a run from a cell's centre to any point of it or to the
//! boundary stays within it. A point within `t` through space is within
//! `t`, so a point of the hole within `t` of its boundary on the map drawn
//! to scale, divided by the stretch, is within `t`. The boundary is drawn
//! as a polyline, each pcurve divided until its middle points stand off
//! their chords by no more than a quarter of that reach, and the reach
//! asked of the cells is less that.
//!
//! Only a hole no other face walks an edge of goes: a hole whose edges a
//! neighbour shares bounds that neighbour too. The face must check clean
//! without it.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{BSplineSurface, Curve, Curve2d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{KnotVector, Point, Point2, Vector2, Weighted};
use ogeom_topo::{EdgeRepr, Model, Orientation, Shape, ShapeType, explore_unique};

use crate::Reshape;

/// The most cells a hole is covered by before it is kept as unproven: a
/// cost, not a measure of the hole.
const CELLS: usize = 1 << 16;

/// The most pieces a pcurve is divided into before the hole is kept as
/// undrawn: a cost, not a measure of the hole.
const PIECES: usize = 1 << 12;

/// The most times a hole is drawn again for the bound over the box it
/// reaches before it is kept as undrawn: a cost, not a measure of the hole.
const ROUNDS: usize = 8;

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
        let reach = data.tolerance.get();
        let stored = face.oriented(Orientation::Forward);
        let wires = model.ordered_children_of(&stored)?;
        let mut gone = Vec::new();
        for wire in wires.iter().skip(1) {
            let edges = model.ordered_children_of(&wire.oriented(Orientation::Forward))?;
            if edges.iter().any(|e| uses.get(&e.node()) != Some(&1)) {
                continue;
            }
            let Some(hole) = charted(model, data.surface, surface, &edges, reach, tol)? else {
                continue;
            };
            let length: f64 = hole.ring.windows(2).map(|w| w[0].distance(w[1])).sum();
            if crate::divide::signed_area(&hole.ring).abs() > hole.reach * length {
                continue;
            }
            if within_reach_of_boundary(&hole.ring, hole.reach - hole.sag, hole.within) {
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

/// How far a unit step in each chart parameter carries a point along a
/// surface at most, over a box of its chart.
#[derive(Debug, Clone, Copy)]
struct Metric {
    /// The bound in `u`, then in `v`.
    scale: Vector2,
    /// Whether the two directions meet square throughout.
    square: bool,
}

impl Metric {
    /// How much farther than the chart drawn to scale a run may carry a
    /// point through space.
    fn stretch(self) -> f64 {
        if self.square {
            1.0
        } else {
            std::f64::consts::SQRT_2
        }
    }

    /// Whether this bounds at least as widely as `other`.
    fn covers(self, other: Self) -> bool {
        self.scale.x >= other.scale.x
            && self.scale.y >= other.scale.y
            && (other.square || !self.square)
    }

    /// The wider of the two bounds.
    fn widest(self, other: Self) -> Self {
        Self {
            scale: Vector2::new(
                self.scale.x.max(other.scale.x),
                self.scale.y.max(other.scale.y),
            ),
            square: self.square && other.square,
        }
    }
}

/// A hole drawn to scale.
struct Hole {
    /// The boundary as a closed polyline on the chart drawn to scale.
    ring: Vec<Point2>,
    /// The widest its pieces' middles stand off their chords.
    sag: f64,
    /// The box of the chart drawn to scale the hole may reach.
    within: (Point2, Point2),
    /// The face's tolerance divided by the stretch: the reach on the map.
    reach: f64,
}

/// The hole walked by `edges` drawn to the scale of the bound over the box
/// it reaches (see the module documentation); `None` where the surface has
/// no bound here, an edge has no plain pcurve on the surface, the ring
/// does not close in one copy of the chart, or it cannot be drawn within a
/// quarter of the face's tolerance `reach`.
fn charted(
    model: &Model,
    surface_id: ogeom_topo::SurfaceId,
    surface: &SurfaceGeometry,
    edges: &[Shape],
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Hole>> {
    let Some(start) = drawn(model, surface_id, surface, edges, None, 0.0, tol)? else {
        return Ok(None);
    };
    let Some(mut metric) = metric_over(surface, start.low, start.high, tol)? else {
        return Ok(None);
    };
    for _ in 0..ROUNDS {
        let on_map = reach / metric.stretch();
        let Some(drawing) = drawn(
            model,
            surface_id,
            surface,
            edges,
            Some(metric.scale),
            on_map,
            tol,
        )?
        else {
            return Ok(None);
        };
        let Some(over) = metric_over(surface, drawing.low, drawing.high, tol)? else {
            return Ok(None);
        };
        if metric.covers(over) {
            let to_scale = |p: Point2| Point2::new(p.x * metric.scale.x, p.y * metric.scale.y);
            return Ok(Some(Hole {
                ring: drawing.ring,
                sag: drawing.sag,
                within: (to_scale(drawing.low), to_scale(drawing.high)),
                reach: on_map,
            }));
        }
        metric = metric.widest(over);
    }
    Ok(None)
}

/// A boundary drawn.
struct Drawing {
    /// The closed polyline, drawn to scale.
    ring: Vec<Point2>,
    /// The widest its pieces' middles stand off their chords, to scale.
    sag: f64,
    /// The corners of the chart box the boundary may reach: its drawing's,
    /// widened by the sag.
    low: Point2,
    high: Point2,
}

/// The hole walked by `edges` as a closed polyline on the chart drawn to
/// `scale`, each pcurve divided until its pieces' middles stand off their
/// chords by no more than a quarter of `reach`; the ends of the edges
/// alone, unscaled, where no scale is given. `None` where an edge has no
/// plain pcurve on the surface, the ring does not close in one copy of the
/// chart, or it cannot be drawn so.
fn drawn(
    model: &Model,
    surface_id: ogeom_topo::SurfaceId,
    surface: &SurfaceGeometry,
    edges: &[Shape],
    scale: Option<Vector2>,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Drawing>> {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let period = |periodic: bool, closed: bool, span: f64| {
        (periodic || closed)
            .then_some(span)
            .filter(|p| *p > 0.0 && p.is_finite())
    };
    let periods = (
        period(surface.is_periodic_u(), surface.is_closed_u(tol), u1 - u0),
        period(surface.is_periodic_v(), surface.is_closed_v(tol), v1 - v0),
    );
    let carry = |by: f64, period: Option<f64>| period.map_or(0.0, |p| (by / p).round() * p);
    let carried = |from: Point2, to: Point2| {
        Vector2::new(
            carry(to.x - from.x, periods.0),
            carry(to.y - from.y, periods.1),
        )
    };
    let unit = Vector2::new(1.0, 1.0);
    let stretch = scale.unwrap_or(unit);
    let to_scale = |p: Point2| Point2::new(p.x * stretch.x, p.y * stretch.y);
    let mut chart: Vec<Point2> = Vec::new();
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
                points.push(at(s)?);
                if k > 0 && scale.is_some() {
                    let middle = to_scale(at(s - 0.5 * step)?);
                    let chord = to_scale(points[k - 1]).midpoint(to_scale(points[k]));
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
        let by = match chart.last() {
            Some(&last) => carried(points[0], last),
            None => Vector2::new(0.0, 0.0),
        };
        chart.extend(points.into_iter().map(|p| p + by));
    }
    let (Some(&first), Some(&last)) = (chart.first(), chart.last()) else {
        return Ok(None);
    };
    if carried(first, last) != Vector2::new(0.0, 0.0) {
        return Ok(None);
    }
    let (mut low, mut high) = (first, first);
    for p in &chart {
        low = Point2::new(low.x.min(p.x), low.y.min(p.y));
        high = Point2::new(high.x.max(p.x), high.y.max(p.y));
    }
    let margin = Vector2::new(sag / stretch.x, sag / stretch.y);
    Ok(Some(Drawing {
        ring: chart.into_iter().map(to_scale).collect(),
        sag,
        low: low - margin,
        high: high + margin,
    }))
}

/// The bound on how far a unit step in each chart parameter carries a
/// point of `surface` over the chart box from `low` to `high` (see the
/// module documentation); `None` on a surface with no such bound here, or
/// where it is not positive both ways.
fn metric_over(
    surface: &SurfaceGeometry,
    low: Point2,
    high: Point2,
    tol: Tolerances,
) -> OgeomResult<Option<Metric>> {
    let metric = match surface {
        SurfaceGeometry::Plane(_) | SurfaceGeometry::Cylinder(_) => {
            let ((u0, u1), (v0, v1)) = surface.domain();
            let (du, dv) = surface.d1_at(mid(u0, u1), mid(v0, v1), tol)?;
            let (su, sv) = (du.magnitude(), dv.magnitude());
            if du.dot(dv).abs() > tol.angular() * su * sv {
                return Ok(None);
            }
            Metric {
                scale: Vector2::new(su, sv),
                square: true,
            }
        }
        SurfaceGeometry::Cone(c) => {
            let cone = c.cone();
            let round = cone
                .radius_at(low.y)
                .abs()
                .max(cone.radius_at(high.y).abs());
            Metric {
                scale: Vector2::new(round, cone.half_angle().tan().hypot(1.0)),
                square: true,
            }
        }
        SurfaceGeometry::Sphere(s) => {
            let radius = s.sphere().radius();
            let (least, most) = cosine_range(low.y, high.y);
            Metric {
                scale: Vector2::new(radius * least.abs().max(most.abs()), radius),
                square: true,
            }
        }
        SurfaceGeometry::Torus(t) => {
            let torus = t.torus();
            let (major, minor) = (torus.major_radius(), torus.minor_radius());
            let (least, most) = cosine_range(low.y, high.y);
            let round = minor
                .mul_add(least, major)
                .abs()
                .max(minor.mul_add(most, major).abs());
            Metric {
                scale: Vector2::new(round, minor),
                square: true,
            }
        }
        SurfaceGeometry::BSpline(patch) => match patch_metric(patch, low, high, tol) {
            Some(metric) => metric,
            None => return Ok(None),
        },
        SurfaceGeometry::Revolution(r) => {
            let Curve::BSpline(curve) = r.curve() else {
                return Ok(None);
            };
            let Some((first, last)) = active(curve.knots(), low.y, high.y, tol) else {
                return Ok(None);
            };
            let points = &curve.control_points()[first..=last];
            let axis = r.axis();
            let round = points
                .iter()
                .map(|w| {
                    let p = w.point();
                    p.distance(axis.project(p))
                })
                .fold(0.0, f64::max);
            let (spread, lightest) = extent(points.iter().copied());
            Metric {
                scale: Vector2::new(round, speed(curve.knots(), points, first, spread, lightest)),
                square: false,
            }
        }
        SurfaceGeometry::Extrusion(e) => {
            let Curve::BSpline(curve) = e.curve() else {
                return Ok(None);
            };
            let Some((first, last)) = active(curve.knots(), low.x, high.x, tol) else {
                return Ok(None);
            };
            let points = &curve.control_points()[first..=last];
            let (spread, lightest) = extent(points.iter().copied());
            Metric {
                scale: Vector2::new(speed(curve.knots(), points, first, spread, lightest), 1.0),
                square: false,
            }
        }
        SurfaceGeometry::Trimmed(t) => return metric_over(t.basis(), low, high, tol),
        SurfaceGeometry::Offset(_) => return Ok(None),
    };
    let positive = |s: f64| s > 0.0 && s.is_finite();
    Ok((positive(metric.scale.x) && positive(metric.scale.y)).then_some(metric))
}

/// The bound over the box from `low` to `high` of a B-spline patch, from
/// its control net (see the module documentation); `None` where the box
/// leaves the patch's domain.
fn patch_metric(
    patch: &BSplineSurface,
    low: Point2,
    high: Point2,
    tol: Tolerances,
) -> Option<Metric> {
    let (u_knots, v_knots) = (patch.u_knots(), patch.v_knots());
    let (i0, i1) = active(u_knots, low.x, high.x, tol)?;
    let (j0, j1) = active(v_knots, low.y, high.y, tol)?;
    let grid = patch.grid();
    let net: Vec<Vec<Weighted<Point>>> = (i0..=i1)
        .map(|i| {
            (j0..=j1)
                .map(|j| grid.get(i, j))
                .collect::<Option<Vec<_>>>()
        })
        .collect::<Option<Vec<_>>>()?;
    let (spread, lightest) = extent(net.iter().flatten().copied());
    let mut across: f64 = 0.0;
    for j in 0..net[0].len() {
        let column: Vec<Weighted<Point>> = net.iter().map(|row| row[j]).collect();
        across = across.max(speed(u_knots, &column, i0, spread, lightest));
    }
    let along = net
        .iter()
        .map(|row| speed(v_knots, row, j0, spread, lightest))
        .fold(0.0, f64::max);
    Some(Metric {
        scale: Vector2::new(across, along),
        square: false,
    })
}

/// The first and last control points whose basis functions are not zero
/// somewhere from `low` to `high`; `None` where that leaves the domain.
fn active(knots: &KnotVector, low: f64, high: f64, tol: Tolerances) -> Option<(usize, usize)> {
    let first = knots.span(low, tol).ok()?;
    let last = knots.span(high, tol).ok()?;
    Some((first.checked_sub(knots.degree())?, last))
}

/// The extent of the control points (the diagonal of their box, which no
/// two of them nor any point among them stand farther apart than) and
/// their lightest weight.
fn extent(points: impl Iterator<Item = Weighted<Point>>) -> (f64, f64) {
    let mut low = Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut high = Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut lightest = f64::INFINITY;
    for w in points {
        let p = w.point();
        low = Point::new(low.x.min(p.x), low.y.min(p.y), low.z.min(p.z));
        high = Point::new(high.x.max(p.x), high.y.max(p.y), high.z.max(p.z));
        lightest = lightest.min(w.weight);
    }
    (low.distance(high), lightest)
}

/// The most a unit step along `knots` carries a point of the B-spline
/// whose control points from index `first` on are `row`, over the stretch
/// where only those are not zero: `spread` bounds how far any of them
/// stands from any point of it, and `lightest` its weights (see the module
/// documentation).
fn speed(
    knots: &KnotVector,
    row: &[Weighted<Point>],
    first: usize,
    spread: f64,
    lightest: f64,
) -> f64 {
    let degree = knots.degree();
    let t = knots.knots();
    let mut most: f64 = 0.0;
    for (n, pair) in row.windows(2).enumerate() {
        let i = first + n;
        let span = t[i + degree + 1] - t[i + 1];
        if span <= 0.0 {
            continue;
        }
        let (a, b) = (pair[0], pair[1]);
        let step = b.weight.mul_add(
            a.point().distance(b.point()),
            (b.weight - a.weight).abs() * spread,
        );
        #[allow(clippy::cast_precision_loss, reason = "a degree")]
        let rate = degree as f64 * step / span;
        most = most.max(rate);
    }
    most / lightest
}

/// The least and the most `cos` reaches from `a` to `b`.
fn cosine_range(a: f64, b: f64) -> (f64, f64) {
    use std::f64::consts::{PI, TAU};
    let (ca, cb) = (a.cos(), b.cos());
    let crest = (a / TAU).ceil() * TAU <= b;
    let trough = ((a - PI) / TAU).ceil().mul_add(TAU, PI) <= b;
    (
        if trough { -1.0 } else { ca.min(cb) },
        if crest { 1.0 } else { ca.max(cb) },
    )
}

/// The middle of a range, or 0 where it is unbounded.
fn mid(a: f64, b: f64) -> f64 {
    let m = f64::midpoint(a, b);
    if m.is_finite() { m } else { 0.0 }
}

/// Whether every point inside the closed polyline `ring` stands within
/// `reach` of it, proved cell by cell over the cells' parts within the box
/// `within` (see the module documentation); `false` where it is not, or
/// where the proof runs past its cost.
fn within_reach_of_boundary(ring: &[Point2], reach: f64, within: (Point2, Point2)) -> bool {
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
        let from = Point2::new(
            (centre.x - half).max(within.0.x),
            (centre.y - half).max(within.0.y),
        );
        let to = Point2::new(
            (centre.x + half).min(within.1.x),
            (centre.y + half).min(within.1.y),
        );
        if from.x > to.x || from.y > to.y {
            continue;
        }
        let (middle, corner) = (from.midpoint(to), 0.5 * from.distance(to));
        let near = nearest(ring, middle);
        if near + corner <= reach {
            continue;
        }
        let inside = winds(ring, middle);
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
