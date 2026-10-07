//! Seams stored as plain pcurves, joined.
//!
//! An edge along a closed surface's seam that the surface's faces walk
//! once each way needs both columns of the seam, and a seam
//! representation holds them so each walk finds its own. Stored instead
//! as plain pcurves on the one surface (two a period apart, or one or two
//! on the same column), the lookup by surface finds the first for both
//! walks, and one of them reads its boundary a period away.
//! [`join_seam_columns`] puts such an edge right.

use std::collections::HashMap;

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve2d as _, Surface as _};
use ogeom_math::{Point2, Transform2, Vector2};
use ogeom_topo::{
    EdgeData, EdgeRepr, Location, Model, NodeData, Orientation, PCurveId, Shape, ShapeType,
    SurfaceId, TShapeId, explore_unique,
};

/// Make a seam of every edge that holds its column on a closed surface as
/// plain pcurves where the surface's faces walk it twice, once each way,
/// on two columns a period apart: two faces joined along the surface's
/// seam, or one face going round it. Returns how many edges were made
/// seams.
///
/// The edge's plain pcurves on the surface, one, or two on one column or a
/// whole period apart, state where the edge stands in the chart. Each walk
/// takes the copy of that column, a whole period over or none, whose ends
/// meet the walk's neighbouring edges. Where the two walks take different
/// copies, the edge becomes a seam carrying both over one range, each on
/// the side walked its way, and a lookup by surface finds the right one
/// for each walk. An edge whose copies the neighbours do not tell apart,
/// whose walks take one copy, or that the surface's faces walk other than
/// twice, is left as it is.
///
/// Nodes below `shape` that other shapes hold as well are copied first
/// ([`Model::unshare`]), so the repair reaches no other shape.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if the shape
/// or anything it names does not resolve.
pub fn join_seam_columns(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<usize> {
    model.unshare(shape)?;
    // Each edge's uses by faces on their surface: the face's boundary ring
    // as it stores it, and where in the ring the edge stands.
    let mut uses = HashMap::<(TShapeId, SurfaceId), Vec<(Vec<Shape>, usize)>>::new();
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(surface) = model
            .node(&face)
            .and_then(|n| n.data().as_face())
            .map(|d| d.surface)
        else {
            ogeom_bail!(Dangling, "face is not in this model");
        };
        let stored = face.oriented(Orientation::Forward);
        for wire in model.ordered_children_of(&stored)? {
            let ring = model.ordered_children_of(&wire)?;
            for (at, edge) in ring.iter().enumerate() {
                uses.entry((edge.node(), surface))
                    .or_default()
                    .push((ring.clone(), at));
            }
        }
    }
    let mut joined = 0;
    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            ogeom_bail!(Dangling, "edge is not in this model");
        };
        for stated in stated_columns(model, data, tol)? {
            let Some([(ring_x, at_x), (ring_y, at_y)]) = uses
                .get(&(edge.node(), stated.surface))
                .and_then(|faces| <&[_; 2]>::try_from(&faces[..]).ok())
            else {
                continue;
            };
            let (walk_x, walk_y) = (ring_x[*at_x].orientation(), ring_y[*at_y].orientation());
            if walk_x == walk_y {
                continue;
            }
            let (Some(sx), Some(sy)) = (
                copy_met(model, &stated, ring_x, *at_x, tol)?,
                copy_met(model, &stated, ring_y, *at_y, tol)?,
            ) else {
                continue;
            };
            if (sx - sy).magnitude() <= tol.parametric() {
                continue;
            }
            let Some(column) = model.geometry().pcurve(stated.pcurve).cloned() else {
                ogeom_bail!(Dangling, "an edge names geometry not in this model");
            };
            let mut copy_at = |shift: Vector2| -> OgeomResult<PCurveId> {
                if shift.magnitude() <= tol.parametric() {
                    return Ok(stated.pcurve);
                }
                let moved = column.transformed(&Transform2::translation(shift), tol)?;
                Ok(model.geometry_mut().add_pcurve(moved))
            };
            let (forward, reversed) = if walk_x == Orientation::Forward {
                (copy_at(sx)?, copy_at(sy)?)
            } else {
                (copy_at(sy)?, copy_at(sx)?)
            };
            let Some(node) = model.node_mut(&edge) else {
                ogeom_bail!(Dangling, "edge is not in this model");
            };
            let NodeData::Edge(data) = node.data_mut() else {
                ogeom_bail!(Dangling, "edge node holds no edge data");
            };
            // The columns stand where the edge's pcurves stood, or a whole
            // period over, so its same-parameter claim stands as it was.
            data.representations[stated.at[0]] = EdgeRepr::Seam {
                forward,
                reversed,
                surface: stated.surface,
                location: stated.location,
                range: stated.range,
            };
            if let Some(&second) = stated.at.get(1) {
                data.representations.remove(second);
            }
            joined += 1;
            break;
        }
    }
    Ok(joined)
}

/// An edge's column on a closed surface as its plain pcurves state it: the
/// pcurves' places among its representations, the first's pcurve, range
/// and ends in the chart, and the shifts by whole periods a walk may read
/// the column at.
struct StatedColumn {
    at: Vec<usize>,
    pcurve: PCurveId,
    surface: SurfaceId,
    location: Location,
    range: (f64, f64),
    ends: (Point2, Point2),
    shifts: Vec<Vector2>,
}

/// The edge's columns on closed surfaces: per surface and placement, one
/// plain pcurve, or two that stand on one column or a whole period apart
/// all along, each read over its own range.
fn stated_columns(
    model: &Model,
    data: &EdgeData,
    tol: Tolerances,
) -> OgeomResult<Vec<StatedColumn>> {
    const SAMPLES: u32 = 8;
    let plain: Vec<_> = data
        .representations
        .iter()
        .enumerate()
        .filter_map(|(k, r)| match r {
            EdgeRepr::PCurve {
                curve,
                surface,
                location,
                range,
            } => Some((k, *curve, *surface, location, *range)),
            _ => None,
        })
        .collect();
    let mut stated = Vec::new();
    for a in &plain {
        let group: Vec<_> = plain.iter().filter(|p| p.2 == a.2 && p.3 == a.3).collect();
        // Each group once, from its first member.
        if group[0].0 != a.0 || group.len() > 2 {
            continue;
        }
        let Some(surface) = model.geometry().surface(a.2) else {
            ogeom_bail!(Dangling, "an edge names geometry not in this model");
        };
        let ((u0, u1), (v0, v1)) = surface.domain();
        let periods = [
            (surface.is_periodic_u() || surface.is_closed_u(tol))
                .then_some(Vector2::new(u1 - u0, 0.0)),
            (surface.is_periodic_v() || surface.is_closed_v(tol))
                .then_some(Vector2::new(0.0, v1 - v0)),
        ];
        let periods: Vec<Vector2> = periods
            .into_iter()
            .flatten()
            .filter(|p| p.magnitude() > 0.0 && p.magnitude().is_finite())
            .collect();
        if periods.is_empty() {
            continue;
        }
        let Some(pa) = model.geometry().pcurve(a.1) else {
            ogeom_bail!(Dangling, "an edge names geometry not in this model");
        };
        if let Some(b) = group.get(1) {
            let Some(pb) = model.geometry().pcurve(b.1) else {
                ogeom_bail!(Dangling, "an edge names geometry not in this model");
            };
            if !a_whole_period_apart(pa, a.4, pb, b.4, &periods, SAMPLES, tol)? {
                continue;
            }
        }
        let mut shifts = vec![Vector2::new(0.0, 0.0)];
        for p in &periods {
            let more: Vec<Vector2> = shifts.iter().flat_map(|s| [*s + *p, *s - *p]).collect();
            shifts.extend(more);
        }
        stated.push(StatedColumn {
            at: group.iter().map(|p| p.0).collect(),
            pcurve: a.1,
            surface: a.2,
            location: a.3.clone(),
            range: a.4,
            ends: (pa.point_at(a.4.0, tol)?, pa.point_at(a.4.1, tol)?),
            shifts,
        });
    }
    Ok(stated)
}

/// Whether `b` stands on `a`'s column or a whole number of `periods` from
/// it, all along, each read over its own range.
fn a_whole_period_apart(
    a: &ogeom_geom::PlanarCurve,
    ra: (f64, f64),
    b: &ogeom_geom::PlanarCurve,
    rb: (f64, f64),
    periods: &[Vector2],
    samples: u32,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let mut offset: Option<Vector2> = None;
    for s in 0..=samples {
        let f = f64::from(s) / f64::from(samples);
        let d = b.point_at(rb.0 + (rb.1 - rb.0) * f, tol)?
            - a.point_at(ra.0 + (ra.1 - ra.0) * f, tol)?;
        match offset {
            Some(o) if (d - o).magnitude() > tol.parametric() => return Ok(false),
            Some(_) => {}
            None => offset = Some(d),
        }
    }
    let Some(mut rest) = offset else {
        return Ok(false);
    };
    // What is left once every whole period along each closed direction is
    // taken off.
    for p in periods {
        let along = rest.dot(*p) / p.dot(*p);
        rest -= *p * along.round();
    }
    Ok(rest.magnitude() <= tol.parametric())
}

/// The shift of the edge's column that one walk of it reads: the copy
/// whose ends meet the edges before and after it in the ring, walked as the
/// ring walks it; `None` where another copy meets them about as well.
fn copy_met(
    model: &Model,
    stated: &StatedColumn,
    ring: &[Shape],
    at: usize,
    tol: Tolerances,
) -> OgeomResult<Option<Vector2>> {
    let n = ring.len();
    let before = ends_on(model, stated.surface, &ring[(at + n - 1) % n], tol)?;
    let after = ends_on(model, stated.surface, &ring[(at + 1) % n], tol)?;
    let nearest = |to: Point2, ends: &[(Point2, Point2)], start: bool| {
        ends.iter()
            .map(|e| if start { e.0 } else { e.1 }.distance(to))
            .fold(f64::INFINITY, f64::min)
    };
    let mut misses: Vec<(f64, Vector2)> = stated
        .shifts
        .iter()
        .map(|&shift| {
            let (mut start, mut end) = (stated.ends.0 + shift, stated.ends.1 + shift);
            if ring[at].orientation() == Orientation::Reversed {
                core::mem::swap(&mut start, &mut end);
            }
            (
                nearest(start, &before, false) + nearest(end, &after, true),
                shift,
            )
        })
        .collect();
    misses.sort_by(|a, b| a.0.total_cmp(&b.0));
    // The copy a walk reads meets its neighbours within a fitted trim's
    // slack; every other stands a period off.
    match misses[..] {
        [(near, shift), (far, _), ..]
            if near.is_finite() && far > 10.0 * near.max(tol.parametric()) =>
        {
            Ok(Some(shift))
        }
        _ => Ok(None),
    }
}

/// Every reading of an edge's ends in a surface's chart, as its occurrence
/// walks it: each pcurve the edge holds on the surface, each side of a
/// seam.
fn ends_on(
    model: &Model,
    surface: SurfaceId,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Vec<(Point2, Point2)>> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        ogeom_bail!(Dangling, "edge is not in this model");
    };
    let mut ends = Vec::new();
    for repr in &data.representations {
        let (ids, range) = match repr {
            EdgeRepr::PCurve {
                curve,
                surface: s,
                range,
                ..
            } if *s == surface => (vec![*curve], *range),
            EdgeRepr::Seam {
                forward,
                reversed,
                surface: s,
                range,
                ..
            } if *s == surface => (vec![*forward, *reversed], *range),
            _ => continue,
        };
        for id in ids {
            let Some(pcurve) = model.geometry().pcurve(id) else {
                ogeom_bail!(Dangling, "an edge names geometry not in this model");
            };
            let (a, b) = (
                pcurve.point_at(range.0, tol)?,
                pcurve.point_at(range.1, tol)?,
            );
            ends.push(if edge.orientation() == Orientation::Reversed {
                (b, a)
            } else {
                (a, b)
            });
        }
    }
    Ok(ends)
}
