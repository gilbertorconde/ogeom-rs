//! Lacking edges made where a wire gaps past what merging its vertices may
//! close.
//!
//! A gap wider than a fitted trim may stand off its surface is a missing
//! piece of boundary, not slop. It is filled with a straight edge only
//! where that edge is the one the boundary lacks beyond doubt: the gap is
//! shorter than every edge of each face it is in, no other vertex of those
//! faces stands as near either of its ends as the other does, the straight
//! line lies on every such face's surface within the tolerance the edge is
//! given, read along the chart segment between where the edges either side
//! leave off, and each face, its gaps filled, checks clean.

use ogeom_algo::{edge_vertices, linear_properties, pcurve_gap};
use ogeom_core::{FastMap, OgeomResult, Tolerance, Tolerances};
use ogeom_geom::{Curve2d as _, Line2d, LineCurve, SurfaceGeometry};
use ogeom_math::{Point, Point2};
use ogeom_mesh::Deflection;
use ogeom_topo::{
    EdgeData, EdgeRepr, Location, Model, Shape, ShapeType, SurfaceId, TShapeId, explore_unique,
};

use crate::Reshape;

/// A gap in one face's boundary: after which edge of which wire, between
/// the vertex that edge ends at and the one the next begins at, and where
/// the two edges leave off in the face's chart.
struct Gap {
    face: usize,
    wire: usize,
    after: usize,
    from: Shape,
    to: Shape,
    chart: (Point2, Point2),
}

/// One face's boundary as read for its gaps.
struct Boundary {
    face: Shape,
    surface: SurfaceId,
    wires: Vec<(Shape, Vec<Shape>)>,
    shortest: f64,
    points: Vec<(TShapeId, Point)>,
}

/// Stage, for every wire with a gap past `cap` that a straight edge
/// provably fills (see the module documentation), a rebuilt wire with that
/// edge in it; returns how many edges were made. A gap two faces share
/// (they walk the missing edge's two ends each the other way) is filled by
/// one edge with a pcurve on each, or not at all.
pub(crate) fn stage_lacking_edges(
    model: &mut Model,
    shape: &Shape,
    reshape: &mut Reshape,
    cap: f64,
    tol: Tolerances,
) -> OgeomResult<usize> {
    let mut boundaries: Vec<Boundary> = Vec::new();
    let mut gaps: Vec<Gap> = Vec::new();
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        if !face.location().is_identity() || !data.location.is_identity() {
            continue;
        }
        let surface = data.surface;
        let stored = face.oriented(ogeom_topo::Orientation::Forward);
        let mut found = Vec::new();
        let mut wires = Vec::new();
        for (w, wire) in model.ordered_children_of(&stored)?.into_iter().enumerate() {
            let edges = model.ordered_children_of(&wire)?;
            for i in 0..edges.len() {
                let (Some((_, end)), Some((next, _))) = (
                    edge_vertices(model, &edges[i])?,
                    edge_vertices(model, &edges[(i + 1) % edges.len()])?,
                ) else {
                    continue;
                };
                if end.node() == next.node() {
                    continue;
                }
                let (Some(p), Some(q)) = (placed(model, &end)?, placed(model, &next)?) else {
                    continue;
                };
                if p.distance(q) <= cap || model.same_position(&end, &next, tol)? {
                    continue;
                }
                found.push((w, i, end, next));
            }
            wires.push((wire, edges));
        }
        if found.is_empty() {
            continue;
        }
        let index = boundaries.len();
        let mut shortest = f64::INFINITY;
        let mut points = Vec::new();
        for (_, edges) in &wires {
            for edge in edges {
                shortest =
                    shortest.min(linear_properties(model, edge, Deflection::default(), tol)?.mass);
                for v in model.children_of(&Shape::of(edge.node()))? {
                    if let Some(p) = placed(model, &v)? {
                        points.push((v.node(), p));
                    }
                }
            }
        }
        for (w, i, from, to) in found {
            let edges = &wires[w].1;
            let (Some(a), Some(b)) = (
                chart_end(model, surface, &edges[i], true, tol)?,
                chart_end(model, surface, &edges[(i + 1) % edges.len()], false, tol)?,
            ) else {
                continue;
            };
            gaps.push(Gap {
                face: index,
                wire: w,
                after: i,
                from,
                to,
                chart: (a, b),
            });
        }
        boundaries.push(Boundary {
            face,
            surface,
            wires,
            shortest,
            points,
        });
    }
    if gaps.is_empty() {
        return Ok(0);
    }

    // Each gap by the pair of vertices it lies between.
    let mut pairs: FastMap<(TShapeId, TShapeId), Vec<usize>> = FastMap::default();
    for (k, gap) in gaps.iter().enumerate() {
        let (a, b) = (gap.from.node(), gap.to.node());
        pairs
            .entry(if a < b { (a, b) } else { (b, a) })
            .or_default()
            .push(k);
    }
    // A pair some edge already joins is no lacking edge.
    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        let bounds = model.children_of(&Shape::of(edge.node()))?;
        if let [first, .., last] = &bounds[..] {
            let (a, b) = (first.node(), last.node());
            pairs.remove(&if a < b { (a, b) } else { (b, a) });
        }
    }

    let mut made: Vec<(Vec<usize>, Shape)> = Vec::new();
    let mut keys: Vec<_> = pairs.keys().copied().collect();
    keys.sort_unstable();
    for key in keys {
        let at = &pairs[&key];
        if let Some(edge) = lacking_edge(model, &boundaries, &gaps, at, tol)? {
            made.push((at.clone(), edge));
        }
    }

    // Each face's wires with its edges made, checked whole before any is
    // staged: a face that does not check clean keeps every gap, and so does
    // every face sharing an edge with it.
    let mut filled: FastMap<(usize, usize), Vec<(usize, Shape)>> = FastMap::default();
    for (at, edge) in &made {
        for &k in at {
            let gap = &gaps[k];
            let first = gap.from.node() == model.children_of(edge)?[0].node();
            let occurrence = if first { edge.clone() } else { edge.reversed() };
            filled
                .entry((gap.face, gap.wire))
                .or_default()
                .push((gap.after, occurrence));
        }
    }
    let mut rebuilt: FastMap<(usize, usize), Shape> = FastMap::default();
    let mut clean: FastMap<usize, bool> = FastMap::default();
    let mut faces: Vec<usize> = filled.keys().map(|(f, _)| *f).collect();
    faces.sort_unstable();
    faces.dedup();
    for f in faces {
        let boundary = &boundaries[f];
        let mut wires = Vec::with_capacity(boundary.wires.len());
        for (w, (wire, edges)) in boundary.wires.iter().enumerate() {
            let Some(inserts) = filled.get(&(f, w)) else {
                wires.push(wire.clone());
                continue;
            };
            let mut ring = Vec::with_capacity(edges.len() + inserts.len());
            for (i, edge) in edges.iter().enumerate() {
                ring.push(edge.clone());
                ring.extend(
                    inserts
                        .iter()
                        .filter(|(after, _)| *after == i)
                        .map(|(_, e)| e.clone()),
                );
            }
            let new = model.add_wire(&ring)?;
            rebuilt.insert((f, w), new.clone());
            wires.push(new);
        }
        let Some(data) = model
            .node(&boundary.face)
            .and_then(|n| n.data().as_face())
            .cloned()
        else {
            continue;
        };
        let candidate = model.add_face(data, &wires)?;
        clean.insert(f, ogeom_algo::check(model, &candidate, tol)?.is_valid());
    }
    // A face that failed spoils every edge made in it, and with each such
    // edge the other face it fills: repeated until nothing more spoils.
    let mut kept: Vec<bool> = vec![true; made.len()];
    loop {
        let mut changed = false;
        for (m, (at, _)) in made.iter().enumerate() {
            if kept[m]
                && at
                    .iter()
                    .any(|&k| !clean.get(&gaps[k].face).copied().unwrap_or(false))
            {
                kept[m] = false;
                for &k in at {
                    clean.insert(gaps[k].face, false);
                }
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    for ((f, w), new) in rebuilt {
        if clean.get(&f).copied().unwrap_or(false) {
            reshape.replace(&boundaries[f].wires[w].0, new);
        }
    }
    Ok(kept.into_iter().filter(|keep| *keep).count())
}

/// The straight edge between the ends of the gaps `at` (one vertex pair,
/// in one or two faces), with a pcurve on each face, where it provably
/// fills them; `None` where it does not.
fn lacking_edge(
    model: &mut Model,
    boundaries: &[Boundary],
    gaps: &[Gap],
    at: &[usize],
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    // Two faces at most, walking the gap each the other way.
    match at {
        [_] => {}
        [a, b] => {
            let (a, b) = (&gaps[*a], &gaps[*b]);
            if a.face == b.face || a.from.node() == b.from.node() {
                return Ok(None);
            }
        }
        _ => return Ok(None),
    }
    let first = &gaps[at[0]];
    let (Some(p), Some(q)) = (placed(model, &first.from)?, placed(model, &first.to)?) else {
        return Ok(None);
    };
    let length = p.distance(q);
    let mut stated = tol.confusion();
    for &k in at {
        let gap = &gaps[k];
        let boundary = &boundaries[gap.face];
        // Shorter than every edge of the face, and no other vertex of it as
        // near either end as the other end is.
        if length >= boundary.shortest {
            return Ok(None);
        }
        let crowded = boundary.points.iter().any(|(id, at)| {
            *id != gap.from.node()
                && *id != gap.to.node()
                && (at.distance(p) <= length || at.distance(q) <= length)
        });
        if crowded {
            return Ok(None);
        }
        if let Some(t) = model.tolerance_of(&boundary.face)? {
            stated = stated.max(t.get());
        }
    }
    let line = match LineCurve::segment(p, q, tol) {
        Ok(line) => ogeom_geom::Curve::from(line),
        Err(_) => return Ok(None),
    };
    let range = (0.0, length);
    // Each face's pcurve runs from where the edge before the gap leaves
    // off to where the edge after it takes up, as the edge runs: from its
    // first vertex to its second.
    let mut pcurves = Vec::with_capacity(at.len());
    for &k in at {
        let gap = &gaps[k];
        let (a, b) = if gap.from.node() == first.from.node() {
            gap.chart
        } else {
            (gap.chart.1, gap.chart.0)
        };
        let Ok(segment) = Line2d::segment(a, b, tol) else {
            return Ok(None);
        };
        let pcurve: ogeom_geom::PlanarCurve = segment.into();
        let domain = pcurve.domain();
        let surface_id = boundaries[gap.face].surface;
        let Some(surface): Option<&SurfaceGeometry> = model.geometry().surface(surface_id) else {
            return Ok(None);
        };
        let Ok(off) = pcurve_gap((&line, range), (&pcurve, domain), surface, stated, tol) else {
            return Ok(None);
        };
        if off > stated {
            return Ok(None);
        }
        pcurves.push((pcurve, domain, surface_id));
    }
    let id = model.geometry_mut().add_curve(line);
    let mut data = EdgeData::on_curve(id, Location::identity(), range);
    data.widen(Tolerance::new(stated)?);
    let edge = model.add_edge(data, &[first.from.clone(), first.to.clone()])?;
    for (pcurve, domain, surface) in pcurves {
        ogeom_algo::attach_pcurve(model, &edge, pcurve, surface, Location::identity(), domain)?;
    }
    Ok(Some(edge))
}

/// Where a vertex stands: its point, placed. `None` for a vertex under a
/// placement of its own, whose point is not in its edges' frame.
fn placed(model: &Model, v: &Shape) -> OgeomResult<Option<Point>> {
    if !v.location().is_identity() {
        return Ok(None);
    }
    Ok(model
        .node(v)
        .and_then(|n| n.data().as_vertex())
        .map(|d| d.point))
}

/// Where `edge`'s unplaced pcurve on `surface` leaves off as the wire walks
/// it (`last`) or takes up (not `last`): a seam's side by the occurrence's
/// orientation. `None` for a placed edge, or one with no unplaced pcurve
/// there.
fn chart_end(
    model: &Model,
    surface: SurfaceId,
    edge: &Shape,
    last: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Point2>> {
    if !edge.location().is_identity() {
        return Ok(None);
    }
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        return Ok(None);
    };
    let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
    let (id, range) = match data.pcurve_for(surface, edge.location()) {
        Some(EdgeRepr::PCurve {
            curve,
            range,
            location,
            ..
        }) if location.is_identity() => (*curve, *range),
        Some(EdgeRepr::Seam {
            forward,
            reversed: back,
            range,
            location,
            ..
        }) if location.is_identity() => (if reversed { *back } else { *forward }, *range),
        _ => return Ok(None),
    };
    let Some(pcurve) = model.geometry().pcurve(id) else {
        return Ok(None);
    };
    // The walk's last point is the range's end, or its start where the
    // edge is walked backwards.
    let at = if last != reversed { range.1 } else { range.0 };
    Ok(Some(pcurve.point_at(at, tol)?))
}
