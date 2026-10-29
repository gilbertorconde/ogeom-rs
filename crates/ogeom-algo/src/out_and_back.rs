//! Where a face's boundary runs out along a path and straight back along
//! it.
//!
//! Two consecutive edges of a wire, the first from one vertex to another
//! and the second from that one back along the first, bound nothing
//! between them: the region the wire encloses is the same with that stretch
//! of boundary as without. The tip is a vertex of the boundary that bounds
//! no corner, every boolean through the face meets a strand there, and a
//! mesher has to strip it before it can triangulate.
//!
//! - A **tail**: the second edge comes back to the vertex the first left
//!   from. Both edges run out and back whole. [`is_tail`] also takes one
//!   edge walked out and straight back (a slit), which the booleans cut and
//!   make with exact volumes; [`out_and_back`] does not report a slit,
//!   only two distinct edges over one path.
//! - A **notch**: one edge runs farther than the other. The shorter one's
//!   far end stands on the longer within the tolerance either edge or that
//!   vertex states, and the stretch of the longer from there to the tip
//!   runs along the shorter. Only that stretch and the shorter edge bound
//!   nothing.
//!
//! A pair is one only where all of these hold, so that nothing a neighbour
//! holds or the chart tells apart is taken for one:
//!
//! - neither edge is degenerate, and the pair's ends are distinct vertices;
//! - no other face of the shape walks either edge, nor this face again,
//!   except as the pair itself (one edge walked out and back, a tail);
//! - neither edge is a seam on the face, and in the face's chart the tip is
//!   one point, and a tail's walk back ends where its walk out began, with
//!   no whole period between them along a direction the surface closes in:
//!   two meridians a period apart on a sphere are one curve in space and
//!   bound the whole sphere between them in the chart;
//! - the two curves stand no farther apart than the wider of the two
//!   edges' tolerances anywhere along the stretch that runs out and back.

use ogeom_core::{FastMap, OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve2d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Point2, Transform};
use ogeom_topo::{
    EdgeRepr, Model, Orientation, Shape, ShapeType, SurfaceId, TShapeId, explore_unique,
};

use crate::build::edge_vertices;
use crate::pcurve_gap::{nearest_point_on_stretch, stretches_apart};

/// Two consecutive edges of a face's wire running out and back along one
/// path. See the module documentation.
#[derive(Debug, Clone, PartialEq)]
pub struct OutAndBack {
    /// The face, as the shape holds it.
    pub face: Shape,
    /// The wire's index among the face's wires as stored.
    pub wire: usize,
    /// The edge run out along and the edge run back along, by index in the
    /// wire's order: the second follows the first, the last wrapping to
    /// the first.
    pub edges: [usize; 2],
    /// Where one of them runs farther than the other: `None` for a tail.
    pub notch: Option<Notch>,
}

/// Where a notch's shorter edge turns back from the longer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Notch {
    /// Which of the pair runs farther: 0 the edge out, 1 the edge back.
    pub longer: usize,
    /// The parameter of the longer edge's curve nearest the shorter one's
    /// far end: from there to the tip, the longer runs along the shorter.
    pub parameter: f64,
}

/// How many times each edge node is walked by the faces of `shape`.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle
/// fails to resolve.
pub fn edge_uses(model: &Model, shape: &Shape) -> OgeomResult<FastMap<TShapeId, usize>> {
    let mut uses: FastMap<TShapeId, usize> = FastMap::default();
    for face in explore_unique(model, shape, ShapeType::Face)? {
        for wire in model.children_of(&face)? {
            for edge in model.children_of(&wire)? {
                *uses.entry(edge.node()).or_default() += 1;
            }
        }
    }
    Ok(uses)
}

/// Every pair of consecutive edges of a face of `shape` that runs out and
/// back along one path, tails of two distinct edges and notches alike, face
/// by face in the shape's order and wire by wire; a slit, one edge walked
/// out and back, is not one. See the module documentation.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle
/// fails to resolve; an evaluation error where a curve cannot be evaluated
/// over its range.
pub fn out_and_back(model: &Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Vec<OutAndBack>> {
    let uses = edge_uses(model, shape)?;
    let mut found = Vec::new();
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(surface_id) = model
            .node(&face)
            .and_then(|n| n.data().as_face())
            .map(|d| d.surface)
        else {
            continue;
        };
        let Some(surface) = model.geometry().surface(surface_id) else {
            continue;
        };
        let stored = face.oriented(Orientation::Forward);
        for (w, wire) in model.ordered_children_of(&stored)?.into_iter().enumerate() {
            let walked = model.ordered_children_of(&wire.oriented(Orientation::Forward))?;
            if walked.len() < 2 {
                continue;
            }
            for i in 0..walked.len() {
                let j = (i + 1) % walked.len();
                if walked.len() == 2 && j == 0 {
                    // Two edges: the pair round the wire's start is the
                    // same pair walked the other way.
                    continue;
                }
                let (out, back) = (&walked[i], &walked[j]);
                if out.node() == back.node() {
                    continue;
                }
                let notch = if is_tail(model, (surface_id, surface), out, back, &uses, tol)? {
                    None
                } else if let Some(notch) =
                    notch_of(model, (surface_id, surface), out, back, &uses, tol)?
                {
                    Some(notch)
                } else {
                    continue;
                };
                found.push(OutAndBack {
                    face: face.clone(),
                    wire: w,
                    edges: [i, j],
                    notch,
                });
            }
        }
    }
    Ok(found)
}

/// Whether `out`, then `back`, walked in that order by one wire of the
/// face on `surface`, are a tail. `uses` counts how many times the faces of
/// the shape walk each edge node. See the module documentation.
///
/// # Errors
///
/// As [`out_and_back`].
pub fn is_tail(
    model: &Model,
    (surface_id, surface): (SurfaceId, &SurfaceGeometry),
    out: &Shape,
    back: &Shape,
    uses: &FastMap<TShapeId, usize>,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let Some((out_data, back_data)) = alone(model, out, back, uses) else {
        return Ok(false);
    };
    let (Some((a, b)), Some((c, d))) = (edge_vertices(model, out)?, edge_vertices(model, back)?)
    else {
        return Ok(false);
    };
    if a.node() == b.node() || c.node() != b.node() || d.node() != a.node() {
        return Ok(false);
    }
    let (Some(out_chart), Some(back_chart)) = (
        chart_ends(model, surface_id, out, tol)?,
        chart_ends(model, surface_id, back, tol)?,
    ) else {
        return Ok(false);
    };
    if !within_a_period(surface, out_chart.0, back_chart.1, tol)
        || !within_a_period(surface, out_chart.1, back_chart.0, tol)
    {
        return Ok(false);
    }
    if out.node() == back.node() {
        return Ok(true);
    }
    let reach = out_data.tolerance.get().max(back_data.tolerance.get());
    Ok(crate::edges_apart(model, out, back, tol)?.is_some_and(|apart| apart <= reach))
}

/// The notch `out`, then `back`, walked in that order by one wire of the
/// face on `surface`, make, where they make one. See [`is_tail`] and the
/// module documentation.
///
/// # Errors
///
/// As [`out_and_back`].
pub fn notch_of(
    model: &Model,
    (surface_id, surface): (SurfaceId, &SurfaceGeometry),
    out: &Shape,
    back: &Shape,
    uses: &FastMap<TShapeId, usize>,
    tol: Tolerances,
) -> OgeomResult<Option<Notch>> {
    let Some((out_data, back_data)) = alone(model, out, back, uses) else {
        return Ok(None);
    };
    if out.node() == back.node() {
        return Ok(None);
    }
    let (Some((a, b)), Some((c, d))) = (edge_vertices(model, out)?, edge_vertices(model, back)?)
    else {
        return Ok(None);
    };
    if a.node() == b.node() || c.node() != b.node() || d.node() == c.node() || d.node() == a.node()
    {
        return Ok(None);
    }
    let (Some(out_chart), Some(back_chart)) = (
        chart_ends(model, surface_id, out, tol)?,
        chart_ends(model, surface_id, back, tol)?,
    ) else {
        return Ok(None);
    };
    if !within_a_period(surface, out_chart.1, back_chart.0, tol) {
        return Ok(None);
    }
    let (Some(out_curve), Some(back_curve)) =
        (placed_curve(model, out)?, placed_curve(model, back)?)
    else {
        return Ok(None);
    };
    let reach = out_data.tolerance.get().max(back_data.tolerance.get());
    // Each way round: the longer edge, the shorter, the shorter's far end,
    // and which end of the longer's range is the tip.
    let tip_at_end = |edge: &Shape, walked_into: bool| {
        (edge.orientation() == Orientation::Reversed) != walked_into
    };
    for (longer, (long_curve, short_curve, far, tip_high)) in [
        (out_curve, back_curve, &d, tip_at_end(out, true)),
        (back_curve, out_curve, &a, tip_at_end(back, false)),
    ]
    .into_iter()
    .enumerate()
    {
        let Some(vertex) = model.node(far).and_then(|n| n.data().as_vertex()) else {
            continue;
        };
        let placed = far.transform(model.datums())?.apply(vertex.point);
        let (curve, range, placement) = long_curve;
        let into = placement.inverse()?;
        let (t, off) = nearest_point_on_stretch(curve, range, into.apply(placed), tol)?;
        let off = off * placement.scale_factor().abs();
        if off > reach.max(vertex.tolerance.get()) {
            continue;
        }
        let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
        if !(t > lo && t < hi) {
            continue;
        }
        let stretch = if tip_high { (t, range.1) } else { (range.0, t) };
        if stretches_apart((curve, stretch, placement), short_curve, tol)? <= reach {
            return Ok(Some(Notch {
                longer,
                parameter: t,
            }));
        }
    }
    Ok(None)
}

/// The data of `out` and `back` where neither is degenerate and no face
/// walks either but the one walking them as a pair: each once, or one
/// edge out and back.
fn alone<'m>(
    model: &'m Model,
    out: &Shape,
    back: &Shape,
    uses: &FastMap<TShapeId, usize>,
) -> Option<(&'m ogeom_topo::EdgeData, &'m ogeom_topo::EdgeData)> {
    let data = |edge: &Shape| model.node(edge).and_then(|n| n.data().as_edge());
    let (out_data, back_data) = (data(out)?, data(back)?);
    if out_data.degenerate || back_data.degenerate {
        return None;
    }
    let used = |edge: &Shape| uses.get(&edge.node()).copied().unwrap_or(0);
    let alone = if out.node() == back.node() {
        used(out) == 2 && out.orientation() != back.orientation()
    } else {
        used(out) == 1 && used(back) == 1
    };
    alone.then_some((out_data, back_data))
}

/// A curve, the range of it an edge runs over, and where the edge's
/// occurrence places it.
type PlacedCurve<'m> = (&'m Curve, (f64, f64), Transform);

/// `edge`'s curve, its range, and where its occurrence places it.
fn placed_curve<'m>(model: &'m Model, edge: &Shape) -> OgeomResult<Option<PlacedCurve<'m>>> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        return Ok(None);
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        return Ok(None);
    };
    let Some(curve) = model.geometry().curve(*curve) else {
        return Ok(None);
    };
    Ok(Some((curve, *range, edge.transform(model.datums())?)))
}

/// Where `edge`'s one plain pcurve on `surface` begins and ends as its
/// occurrence walks it. `None` for a seam, or an edge with no pcurve there.
fn chart_ends(
    model: &Model,
    surface: SurfaceId,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(Point2, Point2)>> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        return Ok(None);
    };
    let Some(EdgeRepr::PCurve {
        curve,
        range,
        location,
        ..
    }) = data.pcurve_for(surface, edge.location())
    else {
        return Ok(None);
    };
    // Two plain pcurves on one surface at one placement are a seam stored
    // as two columns, walked once each way.
    let plain = data
        .representations
        .iter()
        .filter(|r| {
            matches!(r, EdgeRepr::PCurve { .. })
                && r.surface() == Some(surface)
                && r.location() == Some(location)
        })
        .count();
    if plain > 1 {
        return Ok(None);
    }
    let Some(pcurve) = model.geometry().pcurve(*curve) else {
        return Ok(None);
    };
    let (start, end) = (
        pcurve.point_at(range.0, tol)?,
        pcurve.point_at(range.1, tol)?,
    );
    Ok(Some(if edge.orientation() == Orientation::Reversed {
        (end, start)
    } else {
        (start, end)
    }))
}

/// Whether `p` and `q` stand less than half a period apart along every
/// direction `surface` closes in: the same copy of the chart.
fn within_a_period(surface: &SurfaceGeometry, p: Point2, q: Point2, tol: Tolerances) -> bool {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let closes = [
        surface.is_periodic_u() || surface.is_closed_u(tol),
        surface.is_periodic_v() || surface.is_closed_v(tol),
    ];
    [
        (closes[0], u1 - u0, (p.x - q.x).abs()),
        (closes[1], v1 - v0, (p.y - q.y).abs()),
    ]
    .into_iter()
    .all(|(closes, period, apart)| !closes || !period.is_finite() || 2.0 * apart < period)
}
