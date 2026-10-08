//! Pieces of one periodic face that meet across its seam, joined back into
//! one face.
//!
//! The arrangement splits a face in its own chart, where the seam is the
//! chart's edge. A piece of a cylinder reaching across the seam comes out
//! as two pieces, one at each end of the chart, sewn together along the
//! seam edge. They are one face. The joined face lies on a copy of the
//! surface turned about its axis so that the chart's join falls in a part
//! of the turn the face does not reach, and its boundary runs inside that
//! chart without a seam.

use ogeom_core::FastMap;

use ogeom_algo::{History, attach_pcurve, edge_vertices, make_wire};
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::Curve2d as _;
use ogeom_geom::Surface as _;
use ogeom_geom::{PlanarCurve, SurfaceGeometry, Transformable};
use ogeom_math::{Axis, Point2, Transform, Transform2, Vector2};
use ogeom_topo::{
    EdgeRepr, Location, Model, NodeData, Orientation, Shape, ShapeType, SurfaceId, TShapeId,
    explore_unique,
};

/// Samples taken along each boundary edge's pcurve to find the part of
/// the turn a joined face leaves free.
const SAMPLES: usize = 256;

/// How far apart in the chart one edge's end and the next edge's start
/// may lie and still be taken for one point of the boundary.
const CHART_GAP: f64 = 1e-6;

/// The narrowest free part of the turn the joined face's chart join is put
/// in. A face reaching nearer than this all the way round stays split.
const LEAST_GAP: f64 = 1e-2;

/// `v` held to the rows a surface has where they end.
///
/// A fitted image of an edge along a bounded surface's last row runs past
/// it by the fit's own slop, a fraction of a micron, and is asked there
/// at the row it ends on. The turn about the axis moves no row, so a
/// surface and its turned copy are held alike.
fn on_rows(surface: &SurfaceGeometry, v: f64) -> f64 {
    let (_, (low, high)) = surface.domain();
    if surface.is_periodic_v() || low > high {
        v
    } else {
        v.clamp(low, high)
    }
}

/// Join, in each shell, the faces that meet across their surface's seam.
///
/// Two faces are joined when they lie on one surface of revolution (a
/// cylinder, cone, sphere, torus or a full revolution) at no placement,
/// face the same way, and share an edge that is that surface's seam in
/// the angle about the axis. Faces chained that way are joined as a group.
/// A group whose union reaches all the way round, holds an edge any two
/// of its faces share other than such a seam, or whose boundary branches,
/// is left as it is: valid, only split. The returned history takes each
/// joined face to the face it became and deletes the seam edges that
/// dissolved; each shell holding a joined face is rebuilt in place.
///
/// A joined face's edges take pcurves on its turned surface, which no edge
/// of `held` may: a group holding one is refused.
pub(crate) fn join_across_seams(
    model: &mut Model,
    shells: &mut [Shape],
    held: &ogeom_core::FastSet<TShapeId>,
    tol: Tolerances,
) -> OgeomResult<History> {
    let mut history = History::new();
    for shell in shells.iter_mut() {
        let faces = model.children_of(shell)?;
        let mut users: FastMap<TShapeId, Vec<usize>> = FastMap::default();
        for (i, face) in faces.iter().enumerate() {
            for edge in explore_unique(model, face, ShapeType::Edge)? {
                users.entry(edge.node()).or_default().push(i);
            }
        }
        let mut group: Vec<usize> = (0..faces.len()).collect();
        fn root(group: &mut [usize], mut i: usize) -> usize {
            while group[i] != i {
                group[i] = group[group[i]];
                i = group[i];
            }
            i
        }
        let mut edges: Vec<(&TShapeId, &Vec<usize>)> = users.iter().collect();
        edges.sort_by_key(|(id, _)| **id);
        for (&node, faces_on) in edges {
            let [a, b] = faces_on.as_slice() else {
                continue;
            };
            let (Some(sa), Some(sb)) = (
                joinable_surface(model, &faces[*a])?,
                joinable_surface(model, &faces[*b])?,
            ) else {
                continue;
            };
            if sa != sb || faces[*a].orientation() != faces[*b].orientation() {
                continue;
            }
            if !is_angle_seam(model, node, sa, tol)? {
                continue;
            }
            let (ra, rb) = (root(&mut group, *a), root(&mut group, *b));
            group[ra] = rb;
        }
        let mut clusters: FastMap<usize, Vec<usize>> = FastMap::default();
        for i in 0..faces.len() {
            let r = root(&mut group, i);
            clusters.entry(r).or_default().push(i);
        }
        let mut clusters: Vec<Vec<usize>> =
            clusters.into_values().filter(|m| m.len() > 1).collect();
        clusters.sort();
        if clusters.is_empty() {
            continue;
        }
        let mut replaced: FastMap<usize, Option<Shape>> = FastMap::default();
        for members in &clusters {
            let chosen: Vec<Shape> = members.iter().map(|&i| faces[i].clone()).collect();
            if !held.is_empty() {
                for face in &chosen {
                    if explore_unique(model, face, ShapeType::Edge)?
                        .iter()
                        .any(|e| held.contains(&e.node()))
                    {
                        ogeom_core::ogeom_bail!(
                            NotDone,
                            "faces to join across a seam hold an edge shared with a face set aside"
                        );
                    }
                }
            }
            let Some((joined, dissolved)) = join_group(model, &chosen, tol)? else {
                continue;
            };
            for face in &chosen {
                history.modify(face, joined.clone());
            }
            for edge in &dissolved {
                history.delete(edge);
            }
            replaced.insert(members[0], Some(joined));
            for &other in &members[1..] {
                replaced.insert(other, None);
            }
        }
        if replaced.is_empty() {
            continue;
        }
        let kept: Vec<Shape> = faces
            .iter()
            .enumerate()
            .filter_map(|(i, f)| match replaced.get(&i) {
                Some(slot) => slot.clone(),
                None => Some(f.clone()),
            })
            .collect();
        *shell = model.add_shell(&kept)?;
    }
    Ok(history)
}

/// The surface of a face that can be joined across its seam: one turned
/// about an axis, periodic in the angle about it over a full turn, and
/// the face at no placement.
fn joinable_surface(model: &Model, face: &Shape) -> OgeomResult<Option<SurfaceId>> {
    if !face.location().is_identity() {
        return Ok(None);
    }
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    if !data.location.is_identity() {
        return Ok(None);
    }
    Ok(model
        .geometry()
        .surface(data.surface)
        .and_then(axis_of)
        .map(|_| data.surface))
}

/// The axis a surface turns about, when its first parameter is the angle
/// about that axis over one full period.
fn axis_of(surface: &SurfaceGeometry) -> Option<Axis> {
    let axis = match surface {
        SurfaceGeometry::Cylinder(s) => s.cylinder().frame().axis(),
        SurfaceGeometry::Cone(s) => s.cone().frame().axis(),
        SurfaceGeometry::Sphere(s) => s.sphere().frame().axis(),
        SurfaceGeometry::Torus(s) => s.torus().frame().axis(),
        SurfaceGeometry::Revolution(s) => s.axis(),
        _ => return None,
    };
    let (u, _) = surface.domain();
    (surface.is_periodic_u() && ((u.1 - u.0) - core::f64::consts::TAU).abs() <= 1e-9)
        .then_some(axis)
}

/// Whether edge `node` is the seam of `surface` in the angle: its two
/// pcurves one period apart in u and level in v.
fn is_angle_seam(
    model: &Model,
    node: TShapeId,
    surface: SurfaceId,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let Some(NodeData::Edge(data)) = model.node(&Shape::of(node)).map(|n| n.data()) else {
        return Ok(false);
    };
    let Some(EdgeRepr::Seam {
        forward,
        reversed,
        range,
        ..
    }) = data.pcurve_for(surface, &Location::identity())
    else {
        return Ok(false);
    };
    let (Some(f), Some(r)) = (
        model.geometry().pcurve(*forward),
        model.geometry().pcurve(*reversed),
    ) else {
        return Ok(false);
    };
    let middle = 0.5 * (range.0 + range.1);
    let shift = f.point_at(middle, tol)? - r.point_at(middle, tol)?;
    let level = shift.y.abs() <= tol.confusion();
    let turn = (shift.x.abs() - core::f64::consts::TAU).abs() <= 1e-6;
    Ok(level && turn)
}

/// One face for the faces of `members`, or `None` where they do not join
/// into one, with the seam edges that dissolved.
fn join_group(
    model: &mut Model,
    members: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Vec<Shape>)>> {
    let Some(surface_id) = joinable_surface(model, &members[0])? else {
        return Ok(None);
    };
    let Some(surface) = model.geometry().surface(surface_id).cloned() else {
        return Ok(None);
    };
    let Some(axis) = axis_of(&surface) else {
        return Ok(None);
    };
    let period = core::f64::consts::TAU;

    // Every edge occurrence of every member, as its wire runs it.
    let mut occurrences: Vec<Shape> = Vec::new();
    for face in members {
        for wire in model.ordered_children_of(face)? {
            for edge in model.ordered_children_of(&wire)? {
                occurrences.push(edge.composed(wire.orientation()));
            }
        }
    }
    let mut count: FastMap<TShapeId, usize> = FastMap::default();
    for edge in &occurrences {
        *count.entry(edge.node()).or_default() += 1;
    }
    // An edge the group holds twice is a seam that dissolves, or the group
    // is not one face.
    let mut dissolved: Vec<Shape> = Vec::new();
    let mut boundary: Vec<Shape> = Vec::new();
    for edge in &occurrences {
        match count[&edge.node()] {
            1 => boundary.push(edge.clone()),
            2 => {
                if !is_angle_seam(model, edge.node(), surface_id, tol)? {
                    return Ok(None);
                }
                if !dissolved.iter().any(|d| d.node() == edge.node()) {
                    dissolved.push(edge.clone());
                }
            }
            _ => return Ok(None),
        }
    }

    // Each boundary edge's pcurve on the surface, and the angles it covers.
    let mut charted: Vec<(Shape, PlanarCurve, (f64, f64))> = Vec::with_capacity(boundary.len());
    let mut angles: Vec<f64> = Vec::new();
    let mut widest_step = 0.0_f64;
    for edge in &boundary {
        let Some(NodeData::Edge(data)) = model.node(edge).map(|n| n.data()) else {
            return Ok(None);
        };
        let (id, range) = match data.pcurve_for(surface_id, edge.location()) {
            Some(EdgeRepr::PCurve { curve, range, .. }) => (*curve, *range),
            Some(EdgeRepr::Seam { forward, range, .. }) => (*forward, *range),
            _ => return Ok(None),
        };
        let Some(pcurve) = model.geometry().pcurve(id).cloned() else {
            return Ok(None);
        };
        let mut previous: Option<f64> = None;
        for i in 0..=SAMPLES {
            #[allow(clippy::cast_precision_loss)]
            let t = range.0 + (range.1 - range.0) * (i as f64 / SAMPLES as f64);
            let u = pcurve.point_at(t, tol)?.x;
            if let Some(p) = previous {
                widest_step = widest_step.max((u - p).abs());
            }
            previous = Some(u);
            angles.push(u.rem_euclid(period));
        }
        charted.push((edge.clone(), pcurve, range));
    }

    // The widest part of the turn no boundary reaches: the face's angles
    // are its boundary's, so the face does not reach it either.
    angles.sort_by(f64::total_cmp);
    let mut gap = (0.0_f64, 0.0_f64);
    for (i, &a) in angles.iter().enumerate() {
        let next = if i + 1 < angles.len() {
            angles[i + 1]
        } else {
            angles[0] + period
        };
        if next - a > gap.1 - gap.0 {
            gap = (a, next);
        }
    }
    if gap.1 - gap.0 <= LEAST_GAP.max(4.0 * widest_step) {
        return Ok(None);
    }
    let start = 0.5 * (gap.0 + gap.1);

    // The surface turned so its angle starts in the middle of that gap,
    // checked halfway up the face, away from any pole.
    let mut heights = (f64::INFINITY, f64::NEG_INFINITY);
    for (_, pcurve, range) in &charted {
        for t in [range.0, 0.5 * (range.0 + range.1), range.1] {
            let v = pcurve.point_at(t, tol)?.y;
            heights = (heights.0.min(v), heights.1.max(v));
        }
    }
    let height = 0.5 * (heights.0 + heights.1);
    let mut turned = None;
    for sense in [1.0, -1.0] {
        let candidate = surface.transformed(&Transform::rotation(axis, sense * start), tol)?;
        let here = candidate.point_at(0.5, height, tol)?;
        let there = surface.point_at(start + 0.5, height, tol)?;
        if here.distance(there) <= tol.confusion() {
            turned = Some(candidate);
            break;
        }
    }
    let Some(turned) = turned else {
        return Ok(None);
    };

    // Each boundary pcurve moved into the turned chart, where the face's
    // angles all lie strictly inside one period, measured there.
    let mut moved: Vec<(Shape, PlanarCurve, (f64, f64))> = Vec::with_capacity(charted.len());
    for (edge, pcurve, range) in charted {
        let middle = pcurve.point_at(0.5 * (range.0 + range.1), tol)?.x - start;
        let shift = -start - period * (middle / period).floor();
        let shifted =
            pcurve.transformed(&Transform2::translation(Vector2::new(shift, 0.0)), tol)?;
        for i in 0..=SAMPLES / 8 {
            #[allow(clippy::cast_precision_loss)]
            let t = range.0 + (range.1 - range.0) * (i as f64 / (SAMPLES / 8) as f64);
            let (was, now) = (pcurve.point_at(t, tol)?, shifted.point_at(t, tol)?);
            if !(-1e-9..=period + 1e-9).contains(&now.x)
                || surface
                    .point_at(was.x, on_rows(&surface, was.y), tol)?
                    .distance(turned.point_at(now.x, on_rows(&turned, now.y), tol)?)
                    > tol.confusion()
            {
                return Ok(None);
            }
        }
        moved.push((edge, shifted, range));
    }

    // The boundary chained into rings: each edge followed by the one that
    // leaves the vertex it reaches from where it reaches it in the chart. A
    // pole is one vertex with an edge along its row, so the vertex alone
    // does not say which edge comes next.
    let mut ends: Vec<((TShapeId, Point2), (TShapeId, Point2))> = Vec::with_capacity(moved.len());
    for (edge, pcurve, range) in &moved {
        let Some((from, to)) = edge_vertices(model, edge)? else {
            return Ok(None);
        };
        let (a, b) = (
            pcurve.point_at(range.0, tol)?,
            pcurve.point_at(range.1, tol)?,
        );
        let (a, b) = if edge.orientation() == Orientation::Reversed {
            (b, a)
        } else {
            (a, b)
        };
        ends.push(((from.node(), a), (to.node(), b)));
    }
    let after = |at: usize| -> Option<usize> {
        let (node, point) = ends[at].1;
        let mut next = (0..ends.len())
            .filter(|&j| ends[j].0.0 == node && ends[j].0.1.distance(point) <= CHART_GAP);
        let found = next.next()?;
        next.next().is_none().then_some(found)
    };
    let mut used = vec![false; moved.len()];
    let mut rings: Vec<Vec<Shape>> = Vec::new();
    for first in 0..moved.len() {
        if used[first] {
            continue;
        }
        let mut ring = Vec::new();
        let mut at = first;
        loop {
            used[at] = true;
            ring.push(moved[at].0.clone());
            let Some(next) = after(at) else {
                return Ok(None);
            };
            if next == first {
                break;
            }
            if used[next] {
                return Ok(None);
            }
            at = next;
        }
        rings.push(ring);
    }

    let mut wires = Vec::with_capacity(rings.len());
    for ring in &rings {
        wires.push(make_wire(model, ring, tol)?.shape);
    }
    let turned_id = model.geometry_mut().add_surface(turned);
    for (edge, pcurve, range) in moved {
        attach_pcurve(
            model,
            &edge,
            pcurve,
            turned_id,
            edge.location().clone(),
            range,
        )?;
    }
    let Some(NodeData::Face(data)) = model.node(&members[0]).map(|n| n.data().clone()) else {
        return Ok(None);
    };
    let mut data = *data;
    for face in &members[1..] {
        if let Some(NodeData::Face(other)) = model.node(face).map(|n| n.data()) {
            data.tolerance = data.tolerance.widen(other.tolerance);
        }
    }
    data.surface = turned_id;
    data.natural_restriction = false;
    data.triangulation = None;
    let widest = data.tolerance;
    let joined = model.add_face(data, &wires)?;
    model.widen(&joined, widest)?;
    let joined = joined.oriented(members[0].orientation());
    Ok(Some((joined, dissolved)))
}
