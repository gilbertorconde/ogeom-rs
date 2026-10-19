//! The upgrade family: same-domain face unification, collinear edge
//! merging, and tolerance reduction, undoing the splits an operation left
//! behind without changing the shape they describe.

use ogeom_core::{FastMap, FastSet};

use ogeom_algo::{Built, History};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve2d as _, Curve3d as _, PlanarCurve, SurfaceGeometry};
use ogeom_math::{Direction2, Plane, Point2, Transform};
use ogeom_topo::{
    EdgeRepr, Location, Model, NodeData, Orientation, PCurveId, Shape, ShapeType, SurfaceId,
    TShapeId, explore_unique,
};

use crate::reshape::Reshape;

mod curved;

/// Merge adjacent faces lying on one carrier into single faces.
///
/// Two faces qualify when they share an edge and lie on one carrier facing
/// one way, whatever frame each stored it in and whichever way round: the
/// split a boolean or an exchange leaves behind. A carrier is a plane; a
/// cylinder, cone, sphere or torus of one size about one axis; or an
/// extrusion of one basis curve along one line. The merged face keeps the
/// first face's surface, its window grown to hold the others'; the others'
/// pcurves are carried into its chart by the exact map between the two
/// charts (a turn and a shift on a plane; on a curved carrier each
/// parameter run either way and shifted, by a period where it repeats);
/// the shared edges dissolve; the remaining boundary re-chains into an
/// outer ring and one ring per hole. Faces whose boundary does not
/// re-chain so (two rings through one vertex, an edge a third face also
/// bounds, a band closing round its axis) are left split, the rest
/// merged. Faces on any other surface pass through untouched.
///
/// # Errors
///
/// As the model's own builders, where the shape does not resolve in
/// `model`.
pub fn unify_same_domain(
    model: &mut Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<(Built, usize)> {
    unify(model, shape, None, tol)
}

/// [`unify_same_domain`] for the faces `around` alone: each is merged with
/// the faces on its carrier it reaches across shared edges, and every other
/// face of `shape` is left as it is, split or not.
///
/// What an application refining after a feature wants: the faces the
/// feature made unified with their neighbours, and the splits the part
/// already held (kept on purpose, or not its business) untouched. A face
/// of `around` that is not in `shape` names nothing and is passed over.
///
/// # Errors
///
/// As [`unify_same_domain`].
pub fn unify_same_domain_around(
    model: &mut Model,
    shape: &Shape,
    around: &[Shape],
    tol: Tolerances,
) -> OgeomResult<(Built, usize)> {
    unify(model, shape, Some(around), tol)
}

/// The merged faces describe the edges they keep on their surfaces, so
/// what `shape` shares with other shapes is copied first and `around`
/// carried to the copies.
fn unify(
    model: &mut Model,
    shape: &Shape,
    around: Option<&[Shape]>,
    tol: Tolerances,
) -> OgeomResult<(Built, usize)> {
    ogeom_algo::on_own_nodes_with(
        model,
        shape,
        around.unwrap_or_default(),
        |model, shape, parts| unify_own(model, shape, around.map(|_| parts), tol),
        |(built, _)| &mut built.history,
    )
}

fn unify_own(
    model: &mut Model,
    shape: &Shape,
    around: Option<&[Shape]>,
    tol: Tolerances,
) -> OgeomResult<(Built, usize)> {
    let faces = explore_unique(model, shape, ShapeType::Face)?;

    // The faces in play: every face, or those `around` names and every
    // face their carriers reach across shared edges, found as they are
    // reached.
    let in_play: Vec<bool> = match around {
        None => vec![true; faces.len()],
        Some(seeds) => {
            let mut edge_users: FastMap<TShapeId, Vec<usize>> = FastMap::default();
            for (i, face) in faces.iter().enumerate() {
                for edge in explore_unique(model, face, ShapeType::Edge)? {
                    edge_users.entry(edge.node()).or_default().push(i);
                }
            }
            let mut reached = vec![false; faces.len()];
            let mut carriers: FastMap<usize, Option<Carrier>> = FastMap::default();
            let mut queue: Vec<usize> = Vec::new();
            for (i, face) in faces.iter().enumerate() {
                if seeds.iter().any(|s| s.node() == face.node()) {
                    reached[i] = true;
                    queue.push(i);
                }
            }
            while let Some(i) = queue.pop() {
                if let ogeom_core::collections::hash_map::Entry::Vacant(slot) = carriers.entry(i) {
                    slot.insert(carrier_of(model, &faces[i], tol)?);
                }
                for edge in explore_unique(model, &faces[i], ShapeType::Edge)? {
                    for &j in edge_users.get(&edge.node()).map_or(&[][..], Vec::as_slice) {
                        if reached[j] {
                            continue;
                        }
                        if let ogeom_core::collections::hash_map::Entry::Vacant(slot) =
                            carriers.entry(j)
                        {
                            slot.insert(carrier_of(model, &faces[j], tol)?);
                        }
                        if let (Some(Some(ci)), Some(Some(cj))) =
                            (carriers.get(&i), carriers.get(&j))
                            && same_carrier(ci, cj, tol)
                        {
                            reached[j] = true;
                            queue.push(j);
                        }
                    }
                }
            }
            reached
        }
    };

    // The carrier of each face in play: where it lies in the world, which
    // decides the grouping, and the chart it was stored in, which is where
    // any carried or synthesized pcurve has to land.
    let mut carriers: Vec<Option<Carrier>> = Vec::with_capacity(faces.len());
    for (face, &playing) in faces.iter().zip(&in_play) {
        carriers.push(if playing {
            carrier_of(model, face, tol)?
        } else {
            None
        });
    }

    // Union-find over faces sharing an edge on one carrier.
    let mut group: Vec<usize> = (0..faces.len()).collect();
    fn root(group: &mut [usize], mut i: usize) -> usize {
        while group[i] != i {
            group[i] = group[group[i]];
            i = group[i];
        }
        i
    }
    let mut edge_users: FastMap<TShapeId, Vec<usize>> = FastMap::default();
    for (i, face) in faces.iter().enumerate() {
        for edge in explore_unique(model, face, ShapeType::Edge)? {
            edge_users.entry(edge.node()).or_default().push(i);
        }
    }
    for users in edge_users.values() {
        for pair in users.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            // One carrier means one *surface*, not one frame: two halves of
            // a wall carry their own origins on the same flat.
            if let (Some(ca), Some(cb)) = (carriers[a].as_ref(), carriers[b].as_ref())
                && same_carrier(ca, cb, tol)
            {
                let (ra, rb) = (root(&mut group, a), root(&mut group, b));
                group[ra] = rb;
            }
        }
    }

    // Clusters in the order of their first face, so the merged faces enter
    // the arena in the same order on every run.
    let mut cluster_of_root: FastMap<usize, usize> = FastMap::default();
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    for i in 0..faces.len() {
        let r = root(&mut group, i);
        let c = *cluster_of_root.entry(r).or_insert_with(|| {
            clusters.push(Vec::new());
            clusters.len() - 1
        });
        clusters[c].push(i);
    }

    let mut reshape = Reshape::new();
    for members in &clusters {
        if members.len() < 2 {
            continue;
        }
        // A cluster whose boundary does not close into rings the merged
        // face can keep is left split, its faces untouched.
        let Ok(Some(merged)) = merge_cluster(model, &faces, members, &carriers, &edge_users, tol)
        else {
            continue;
        };
        reshape.replace(&faces[members[0]], merged);
        for &other in &members[1..] {
            reshape.remove(&faces[other]);
        }
    }
    if reshape.is_empty() {
        let mut built = Built::from_nothing(shape.clone());
        built.history.modify(shape, shape.clone());
        return Ok((built, 0));
    }
    let unified = reshape.len();
    Ok((reshape.apply(model, shape)?, unified))
}

/// The face `members` make together: one face on the first member's
/// surface, bounded by the rings their outer edges walk. `None` where they
/// make no one face bounding the same material.
///
/// Every edge is walked as the members walk it, seen from the side they
/// all face, so the rings come out wound as the material lies: the outer
/// one counter-clockwise about that side, the holes clockwise. An edge two
/// members walk once each, opposite ways, dissolves; one a face outside
/// the cluster also bounds, or one walked against the rule, keeps the
/// cluster as it is, as does a ring that does not close or a vertex two
/// rings pass through, where no one walk is the boundary.
fn merge_cluster(
    model: &mut Model,
    faces: &[Shape],
    members: &[usize],
    carriers: &[Option<Carrier>],
    edge_users: &FastMap<TShapeId, Vec<usize>>,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let keeper = members[0];
    let Some(carrier) = carriers[keeper].clone() else {
        return Ok(None);
    };
    let Some(NodeData::Face(keeper_data)) = model.node(&faces[keeper]).map(|n| n.data().clone())
    else {
        return Ok(None);
    };

    // Each member's edges as it walks them, seen from the side it faces.
    let mut walks: Vec<(usize, Shape)> = Vec::new();
    let mut uses: FastMap<TShapeId, Vec<usize>> = FastMap::default();
    for &i in members {
        let wires = model.children_of(&faces[i])?;
        if wires.is_empty() {
            return Ok(None);
        }
        for wire in wires {
            for edge in model.ordered_children_of(&wire)? {
                uses.entry(edge.node()).or_default().push(walks.len());
                walks.push((i, edge));
            }
        }
    }
    let mut kept = vec![false; walks.len()];
    for (node, at) in &uses {
        let across = match at.as_slice() {
            [_] => false,
            [a, b] if walks[*a].0 == walks[*b].0 => false,
            [a, b] => {
                if walks[*a].1.orientation() == walks[*b].1.orientation()
                    || edge_users.get(node).map_or(0, Vec::len) != 2
                {
                    return Ok(None);
                }
                true
            }
            _ => return Ok(None),
        };
        if !across {
            for &k in at {
                kept[k] = true;
            }
        }
    }
    let reversed = faces[keeper].orientation() == Orientation::Reversed;
    let rings = match &carrier {
        Carrier::Flat(flat) => {
            let Some(rings) =
                flat_rings(model, faces, members, carriers, &walks, &kept, flat, tol)?
            else {
                return Ok(None);
            };
            rings
        }
        Carrier::Curved(_) => {
            let mut on: Vec<&curved::Curved> = Vec::with_capacity(members.len());
            for &i in members {
                let Some(Carrier::Curved(c)) = &carriers[i] else {
                    return Ok(None);
                };
                on.push(c);
            }
            let cluster = curved::Members {
                members,
                carriers: &on,
            };
            let Some(rings) = curved::merge(model, &cluster, &walks, &kept, reversed, tol)? else {
                return Ok(None);
            };
            rings
        }
    };

    // The rings as the kept face's own node walks them: a member used
    // reversed in the solid walks its node's rings back.
    let mut wires = Vec::with_capacity(rings.len());
    for ring in &rings {
        let ring = ring.iter().map(|&w| walks[w].1.clone());
        let ring: Vec<Shape> = if reversed {
            ring.rev().map(|e| e.reversed()).collect()
        } else {
            ring.collect()
        };
        wires.push(ogeom_algo::make_wire(model, &ring, tol)?.shape);
    }
    // The kept face's mesh covers the kept face alone.
    let mut data = *keeper_data;
    data.triangulation = None;
    Ok(Some(model.add_face(data, &wires)?))
}

/// The rings a cluster on one plane makes, outer first, as indices into
/// `walks`: the kept walks chained head to tail through their vertices, the
/// outer ring the one turning about the side the faces face. Every outer
/// edge is given a pcurve on the keeper's plane, and its window grown to
/// hold the members'. `None` where a vertex starts two walks, a ring does
/// not close, or the rings do not turn as one outer ring round its holes.
#[allow(clippy::too_many_arguments, reason = "the cluster's parts, read once")]
fn flat_rings(
    model: &mut Model,
    faces: &[Shape],
    members: &[usize],
    carriers: &[Option<Carrier>],
    walks: &[(usize, Shape)],
    kept: &[bool],
    flat: &Flat,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<Vec<usize>>>> {
    let Some(NodeData::Face(keeper_data)) =
        model.node(&faces[members[0]]).map(|n| n.data().clone())
    else {
        return Ok(None);
    };
    let surface_id = keeper_data.surface;
    let boundary: Vec<usize> = (0..walks.len()).filter(|&w| kept[w]).collect();
    if boundary.is_empty() {
        return Ok(None);
    }
    let mut seen: FastSet<TShapeId> = FastSet::default();
    let outer_edges: Vec<Shape> = boundary
        .iter()
        .map(|&w| &walks[w].1)
        .filter(|e| seen.insert(e.node()))
        .cloned()
        .collect();

    // Chain the walks into rings, head to tail.
    let mut ends: Vec<(Shape, Shape)> = Vec::with_capacity(boundary.len());
    for &w in &boundary {
        let Some(pair) = ogeom_algo::edge_vertices(model, &walks[w].1)? else {
            return Ok(None);
        };
        ends.push(pair);
    }
    let mut leaving: FastMap<TShapeId, usize> = FastMap::default();
    for (k, (start, _)) in ends.iter().enumerate() {
        if leaving.insert(start.node(), k).is_some() {
            return Ok(None);
        }
    }
    let mut placed = vec![false; boundary.len()];
    let mut rings: Vec<Vec<usize>> = Vec::new();
    for first in 0..boundary.len() {
        if placed[first] {
            continue;
        }
        let mut ring = Vec::new();
        let mut k = first;
        loop {
            placed[k] = true;
            ring.push(boundary[k]);
            let Some(&next) = leaving.get(&ends[k].1.node()) else {
                return Ok(None);
            };
            if next == first {
                break;
            }
            if placed[next] {
                return Ok(None);
            }
            k = next;
        }
        rings.push(ring);
    }

    // The outer ring is the one turning about the side the faces face.
    let mut turns = Vec::with_capacity(rings.len());
    for ring in &rings {
        let edges: Vec<Shape> = ring.iter().map(|&w| walks[w].1.clone()).collect();
        turns.push(ogeom_algo::edges_turn(model, &edges, flat.outward, tol)?);
    }
    let area = tol.confusion() * tol.confusion();
    let mut outer = turns.iter().enumerate().filter(|(_, t)| **t > area);
    let (Some((outer, _)), None) = (outer.next(), outer.next()) else {
        return Ok(None);
    };
    if turns.iter().filter(|t| **t < -area).count() != rings.len() - 1 {
        return Ok(None);
    }
    rings.swap(0, outer);

    // Every outer edge needs a pcurve for the kept surface: the one it
    // carries on another member's chart, moved across, or else a
    // projection.
    let mut charts: Vec<(SurfaceId, Flat)> = Vec::with_capacity(members.len());
    for &i in &members[1..] {
        if let (Some(NodeData::Face(data)), Some(Carrier::Flat(c))) =
            (model.node(&faces[i]).map(|n| n.data()), &carriers[i])
        {
            charts.push((data.surface, c.clone()));
        }
    }
    for edge in &outer_edges {
        ensure_planar_pcurve(model, edge, surface_id, flat, &charts, tol)?;
    }
    widen_to_members(model, surface_id, flat, &charts, tol)?;
    Ok(Some(rings))
}

/// Grow the window of the plane `surface_id` stores, where it is bounded,
/// to hold every window of `charts` moved into its chart: the merged face
/// reaches over all of them, and a point of a face past its surface's
/// window reads as off the face.
fn widen_to_members(
    model: &mut Model,
    surface_id: SurfaceId,
    carrier: &Flat,
    charts: &[(SurfaceId, Flat)],
    tol: Tolerances,
) -> OgeomResult<()> {
    use ogeom_geom::Surface as _;
    let Some(SurfaceGeometry::Plane(kept)) = model.geometry().surface(surface_id).cloned() else {
        return Ok(());
    };
    let ((mut ua, mut ub), (mut va, mut vb)) = kept.domain();
    let was = ((ua, ub), (va, vb));
    for (id, chart) in charts {
        let Some(SurfaceGeometry::Plane(other)) = model.geometry().surface(*id) else {
            continue;
        };
        let Some(map) = chart_to_chart(chart, carrier, tol)? else {
            ogeom_bail!(
                Construction,
                "two faces on one plane do not share its chart"
            );
        };
        let ((a, b), (c, d)) = other.domain();
        for (u, v) in [(a, c), (a, d), (b, c), (b, d)] {
            let at = map.apply(Point2::new(u, v));
            (ua, ub) = (ua.min(at.x), ub.max(at.x));
            (va, vb) = (va.min(at.y), vb.max(at.y));
        }
    }
    if ((ua, ub), (va, vb)) != was
        && let Some(held) = model.geometry_mut().surface_mut(surface_id)
    {
        *held = ogeom_geom::PlaneSurface::over(kept.plane(), (ua, ub), (va, vb))?.into();
    }
    Ok(())
}

/// Where a face lies and which chart it was stored in, for a face on a
/// carrier unification knows.
#[derive(Debug, Clone)]
enum Carrier {
    Flat(Flat),
    Curved(curved::Curved),
}

/// Where a planar face lies and which chart it was stored in.
#[derive(Debug, Clone)]
struct Flat {
    /// The plane in world space: what "same carrier" is decided on.
    world: Plane,
    /// The side the face faces in the world: its plane's normal, turned
    /// back where the face is used reversed. Two faces on one plane are
    /// one carrier only when they face the same way, whichever way each
    /// stored its plane.
    outward: ogeom_math::Vector,
    /// The plane as the surface stores it: the chart pcurves speak in.
    stored: Plane,
    /// What takes the stored chart to the world.
    placement: Transform,
}

/// Whether two carriers are one surface facing one way, whatever their
/// frames and whichever way each stored its surface.
fn same_carrier(a: &Carrier, b: &Carrier, tol: Tolerances) -> bool {
    match (a, b) {
        (Carrier::Flat(a), Carrier::Flat(b)) => {
            a.outward.dot(b.outward) > 1.0 - tol.angular()
                && a.world.signed_distance_to(b.world.frame().origin()).abs() <= tol.confusion()
        }
        (Carrier::Curved(a), Carrier::Curved(b)) => curved::same_carrier(a, b, tol),
        _ => false,
    }
}

/// The carrier of a face, if it lies on one.
fn carrier_of(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<Carrier>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data().clone()) else {
        return Ok(None);
    };
    let Some(SurfaceGeometry::Plane(p)) = model.geometry().surface(data.surface) else {
        return Ok(curved::carrier_of(model, face, tol)?.map(Carrier::Curved));
    };
    let stored = p.plane();
    let placement = face.transform(model.datums())?;
    if (placement.scale_factor().abs() - 1.0).abs() > 1e-9 {
        return Ok(None);
    }
    let world = stored.transformed(&placement, tol)?;
    let normal = world.normal().vector();
    Ok(Some(Carrier::Flat(Flat {
        outward: if face.orientation() == Orientation::Reversed {
            -normal
        } else {
            normal
        },
        world,
        stored,
        placement,
    })))
}

/// Attach a pcurve for `surface_id` to `edge` if it does not carry one.
///
/// A pcurve the edge carries on one of `charts` (another face on the same
/// plane) is moved into the stored plane's frame by the rigid map between
/// the two charts, which keeps its kind and its parameter: an arc stays
/// the arc it was. With none, a line is projected directly, exact, and any
/// other curve is fitted by projection.
///
/// The chart is the *stored* plane's, so the world point is carried back
/// through the face's placement before it is flattened; a pcurve read under
/// that placement then lands where the edge is.
fn ensure_planar_pcurve(
    model: &mut Model,
    edge: &Shape,
    surface_id: SurfaceId,
    carrier: &Flat,
    charts: &[(SurfaceId, Flat)],
    tol: Tolerances,
) -> OgeomResult<()> {
    let (curve, range, has, sibling) = {
        let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
            ogeom_bail!(Construction, "edge holds no edge data");
        };
        let has = data.pcurve_for(surface_id, edge.location()).is_some();
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            ogeom_bail!(Construction, "an edge with no curve cannot be unified over");
        };
        let sibling = charts.iter().find_map(|(surface, chart)| {
            match data.pcurve_for(*surface, edge.location()) {
                Some(EdgeRepr::PCurve { curve, range, .. }) => Some((*curve, *range, chart)),
                _ => None,
            }
        });
        (*curve, *range, has, sibling)
    };
    if has {
        return Ok(());
    }
    if let Some((pcurve, prange, chart)) = sibling
        && let Some(moved) = model.geometry().pcurve(pcurve).cloned()
        && let Some(map) = chart_to_chart(chart, carrier, tol)?
    {
        return ogeom_algo::attach_pcurve(
            model,
            edge,
            moved.transformed(&map, tol)?,
            surface_id,
            edge.location().clone(),
            prange,
        );
    }
    let Some(geometry) = model.geometry().curve(curve).cloned() else {
        ogeom_bail!(Construction, "edge refers to a curve not in this model");
    };
    let back = carrier.placement.inverse()?;
    let pcurve = if matches!(geometry, ogeom_geom::Curve::Line(_)) {
        let flat = |p: ogeom_math::Point| {
            let local = carrier.stored.frame().to_local(back.apply(p));
            Point2::new(local.x, local.y)
        };
        let a = flat(geometry.point_at(range.0, tol)?);
        let b = flat(geometry.point_at(range.1, tol)?);
        chart_line(a, b, range, tol)?
    } else {
        let local = ogeom_geom::Transformable::transformed(&geometry, &back, tol)?;
        let plane = SurfaceGeometry::Plane(ogeom_geom::PlaneSurface::new(carrier.stored));
        ogeom_algo::pcurve_fit::fit_projected_pcurve(&local, range, &plane, tol)?.0
    };
    ogeom_algo::attach_pcurve(
        model,
        edge,
        pcurve,
        surface_id,
        edge.location().clone(),
        range,
    )
}

/// The rigid map taking `from`'s stored chart onto `to`'s, two charts of
/// one plane: a turn and a shift, with a mirror where their frames run
/// opposite ways round. `None` where the charts are not one plane.
fn chart_to_chart(
    from: &Flat,
    to: &Flat,
    tol: Tolerances,
) -> OgeomResult<Option<ogeom_math::Transform2>> {
    use ogeom_math::{Point, Transform2, Vector2};
    let into = to.placement.inverse()?;
    let lift = |u: f64, v: f64| {
        let world = from
            .placement
            .apply(from.stored.frame().to_world(Point::new(u, v, 0.0)));
        to.stored.frame().to_local(into.apply(world))
    };
    let (o, x, y) = (lift(0.0, 0.0), lift(1.0, 0.0), lift(0.0, 1.0));
    if [o, x, y].iter().any(|p| p.z.abs() > tol.confusion()) {
        return Ok(None);
    }
    let ex = Vector2::new(x.x - o.x, x.y - o.y);
    let ey = Vector2::new(y.x - o.x, y.y - o.y);
    let turn = Transform2::rotation(Point2::ORIGIN, ex.y.atan2(ex.x));
    let shift = Transform2::translation(Vector2::new(o.x, o.y));
    let map = if ex.x.mul_add(ey.y, -(ex.y * ey.x)) > 0.0 {
        shift * turn
    } else {
        let flip = Transform2::line_mirror(
            Point2::ORIGIN,
            Direction2::new(Vector2::new(0.0, 1.0), tol)?,
        );
        shift * turn * flip
    };
    // The charts are rigid images of each other, so the map built from the
    // first axis must carry the second one where it lands too.
    let lands = map.apply(Point2::new(0.0, 1.0));
    if lands.distance(Point2::new(y.x, y.y)) > tol.confusion() {
        return Ok(None);
    }
    Ok(Some(map))
}

/// The chart line running from `a` at `range.0` to `b` at `range.1`.
///
/// A [`Line2d`](ogeom_geom::Line2d) reads its parameter from the axis
/// origin outward, so the origin is where the parameter would have been
/// zero, not where the range starts.
fn chart_line(
    a: Point2,
    b: Point2,
    range: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<PlanarCurve> {
    let span = range.1 - range.0;
    let direction = Direction2::new((b - a) / span, tol)?;
    let origin = a - direction.vector() * range.0;
    Ok(
        ogeom_geom::Line2d::over(ogeom_math::Axis2::new(origin, direction), range.0, range.1)?
            .into(),
    )
}

/// Merge chains of edges lying on one curve into single edges.
///
/// Within each wire, consecutive edges continuing one curve (the same
/// stored curve over contiguous ranges, or two collinear lines meeting end
/// to end) become one edge over the joined range, and the vertex between
/// them dissolves. Every pcurve the pair carried is joined with it and then
/// *measured* against the joined curve: a pair whose parameterizations do
/// not join cleanly is left split rather than merged into a face that could
/// no longer be triangulated.
///
/// # Errors
///
/// As the model's own builders.
pub fn merge_edges(
    model: &mut Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<(Built, usize)> {
    // One pass merges disjoint pairs; a chain of three needs the next pass
    // to see the pair the first one made. Eight is well past any real chain.
    let mut current = shape.clone();
    let mut steps: Vec<History> = Vec::new();
    let mut merged = 0;
    for _ in 0..8 {
        let Some((step, joined)) = merge_pass(model, &current, tol)? else {
            break;
        };
        current = step.shape;
        steps.push(step.history);
        merged += joined;
    }
    if steps.is_empty() {
        let mut built = Built::from_nothing(shape.clone());
        built.history.modify(shape, shape.clone());
        return Ok((built, 0));
    }
    Ok((Built::new(current, History::chain(&steps)), merged))
}

/// One merge pass: every disjoint joinable pair, in one rebuild. `None`
/// when nothing joins.
fn merge_pass(
    model: &mut Model,
    shape: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(Built, usize)>> {
    // Which parameterizations each edge is actually *read* on: the faces
    // bounding it, not whatever reprs it accumulated along the way. An
    // earlier unification leaves pcurves on surfaces nothing carries any
    // more, and a join that tried to honour those would refuse work it can
    // do.
    let mut needed: FastMap<TShapeId, Vec<(SurfaceId, Location)>> = FastMap::default();
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data().clone()) else {
            continue;
        };
        for edge in explore_unique(model, &face, ShapeType::Edge)? {
            let slot = (data.surface, edge.location().clone());
            let slots = needed.entry(edge.node()).or_default();
            if !slots.contains(&slot) {
                slots.push(slot);
            }
        }
    }

    let wires = explore_unique(model, shape, ShapeType::Wire)?;
    let mut reshape = Reshape::new();
    let mut merged_nodes: FastSet<TShapeId> = FastSet::default();
    for wire in &wires {
        let edges = model.ordered_children_of(wire)?;
        if edges.len() < 2 {
            continue;
        }
        // A closed wire has no last edge: its end is the start again, and
        // the split sitting across that join is a split like any other.
        let closed = is_closed_wire(model, &edges, tol);
        let pairs = if closed { edges.len() } else { edges.len() - 1 };
        let mut i = 0;
        while i < pairs {
            let (a, b) = (edges[i].clone(), edges[(i + 1) % edges.len()].clone());
            if a.node() == b.node()
                || merged_nodes.contains(&a.node())
                || merged_nodes.contains(&b.node())
            {
                i += 1;
                continue;
            }
            let mut slots: Vec<(SurfaceId, Location)> =
                needed.get(&a.node()).cloned().unwrap_or_default();
            for slot in needed.get(&b.node()).into_iter().flatten() {
                if !slots.contains(slot) {
                    slots.push(slot.clone());
                }
            }
            let Some(join) = joinable(model, &a, &b, &slots, tol)? else {
                i += 1;
                continue;
            };
            let joined = ogeom_algo::make_edge(model, join.curve, join.range, tol)?.shape;
            for pcurve in join.pcurves {
                ogeom_algo::attach_pcurve(
                    model,
                    &joined,
                    pcurve.curve,
                    pcurve.surface,
                    pcurve.location,
                    pcurve.range,
                )?;
            }
            merged_nodes.insert(a.node());
            merged_nodes.insert(b.node());
            reshape.replace(&a, joined);
            reshape.remove(&b);
            i += 2;
        }
    }
    if reshape.is_empty() {
        return Ok(None);
    }
    let joined = reshape.len();
    Ok(Some((reshape.apply(model, shape)?, joined)))
}

/// Whether an ordered edge list closes back on its own start.
fn is_closed_wire(model: &Model, edges: &[Shape], tol: Tolerances) -> bool {
    let ends = |e: &Shape| -> Option<(ogeom_math::Point, ogeom_math::Point)> {
        let vertices = explore_unique(model, e, ShapeType::Vertex).ok()?;
        let first = model.node(vertices.first()?)?.data().as_vertex()?.point;
        let last = model.node(vertices.last()?)?.data().as_vertex()?.point;
        Some((first, last))
    };
    let (Some(first), Some(last)) = (ends(&edges[0]), ends(&edges[edges.len() - 1])) else {
        return false;
    };
    [first.0, first.1]
        .iter()
        .any(|p| p.distance(last.0) <= tol.confusion() || p.distance(last.1) <= tol.confusion())
}

/// One joined edge, ready to build.
struct Join {
    /// The curve the joined edge runs on.
    curve: ogeom_geom::Curve,
    /// The range it runs over.
    range: (f64, f64),
    /// The joined pcurve for every surface the pair was parameterized on.
    pcurves: Vec<PcurveJoin>,
}

/// A joined pcurve, in the slot it goes back into.
struct PcurveJoin {
    surface: SurfaceId,
    location: Location,
    curve: PlanarCurve,
    range: (f64, f64),
}

/// The three-dimensional half of a join: the curve, the range it runs
/// over, and whether `a` runs first along it.
type CurveJoin = (ogeom_geom::Curve, (f64, f64), bool);

/// The three-dimensional half of a join: the curve, its range, and whether
/// `a` runs first along it.
fn joinable_curve(
    model: &Model,
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<CurveJoin>> {
    let read = |e: &Shape| -> Option<(ogeom_topo::CurveId, (f64, f64))> {
        let data = model.node(e)?.data().as_edge()?;
        let EdgeRepr::Curve3d { curve, range, .. } = data.curve3d()? else {
            return None;
        };
        Some((*curve, *range))
    };
    let (Some((ca, ra)), Some((cb, rb))) = (read(a), read(b)) else {
        return Ok(None);
    };
    if ca == cb {
        let Some(geometry) = model.geometry().curve(ca).cloned() else {
            return Ok(None);
        };
        if (ra.1 - rb.0).abs() <= tol.parametric() {
            return Ok(Some((geometry, (ra.0, rb.1), true)));
        }
        if (rb.1 - ra.0).abs() <= tol.parametric() {
            return Ok(Some((geometry, (rb.0, ra.1), false)));
        }
        return Ok(None);
    }
    // Distinct curves: collinear lines meeting end to end still join.
    let (Some(ga), Some(gb)) = (
        model.geometry().curve(ca).cloned(),
        model.geometry().curve(cb).cloned(),
    ) else {
        return Ok(None);
    };
    let (ogeom_geom::Curve::Line(la), ogeom_geom::Curve::Line(lb)) = (&ga, &gb) else {
        return Ok(None);
    };
    if !la.axis().is_collinear(lb.axis(), tol) {
        return Ok(None);
    }
    let (a0, a1) = (ga.point_at(ra.0, tol)?, ga.point_at(ra.1, tol)?);
    let (b0, b1) = (gb.point_at(rb.0, tol)?, gb.point_at(rb.1, tol)?);
    let fresh = |from: ogeom_math::Point,
                 to: ogeom_math::Point,
                 a_first: bool|
     -> OgeomResult<Option<CurveJoin>> {
        let segment = ogeom_geom::LineCurve::segment(from, to, tol)?;
        let range = segment.domain();
        Ok(Some((segment.into(), range, a_first)))
    };
    if a1.distance(b0) <= tol.confusion() {
        return fresh(a0, b1, true);
    }
    if b1.distance(a0) <= tol.confusion() {
        return fresh(b0, a1, false);
    }
    Ok(None)
}

/// Whether two consecutive edges continue one curve, pcurves and all.
///
/// `slots` is what the joined edge will be read on: every one of them has
/// to come out of the join, or the merge would leave a face that cannot be
/// triangulated, worse than the split it came to fix.
fn joinable(
    model: &Model,
    a: &Shape,
    b: &Shape,
    slots: &[(SurfaceId, Location)],
    tol: Tolerances,
) -> OgeomResult<Option<Join>> {
    let Some((curve, range, a_first)) = joinable_curve(model, a, b, tol)? else {
        return Ok(None);
    };
    let (first, second) = if a_first { (a, b) } else { (b, a) };
    let reprs = |e: &Shape| -> Vec<EdgeRepr> {
        model
            .node(e)
            .and_then(|n| n.data().as_edge())
            .map(|d| d.representations.to_vec())
            .unwrap_or_default()
    };
    let (first_reprs, second_reprs) = (reprs(first), reprs(second));
    // A seam is a face's own doubling, not a split; leave those pairs be.
    if first_reprs
        .iter()
        .chain(&second_reprs)
        .any(|r| matches!(r, EdgeRepr::Seam { .. }))
    {
        return Ok(None);
    }
    let parametric = |reprs: &[EdgeRepr]| -> Vec<(SurfaceId, Location, PCurveId, (f64, f64))> {
        reprs
            .iter()
            .filter_map(|r| match r {
                EdgeRepr::PCurve {
                    curve,
                    range,
                    surface,
                    location,
                } => Some((*surface, location.clone(), *curve, *range)),
                _ => None,
            })
            .collect()
    };
    let (ones, twos) = (parametric(&first_reprs), parametric(&second_reprs));
    let mut pcurves = Vec::with_capacity(slots.len());
    for (surface, location) in slots {
        let (surface, location) = (*surface, location.clone());
        let find = |from: &[(SurfaceId, Location, PCurveId, (f64, f64))]| {
            from.iter()
                .find(|(s, l, _, _)| *s == surface && *l == location)
                .cloned()
        };
        let (Some((_, _, id_one, range_one)), Some((_, _, id_two, range_two))) =
            (find(&ones), find(&twos))
        else {
            return Ok(None);
        };
        let (Some(one), Some(two)) = (
            model.geometry().pcurve(id_one).cloned(),
            model.geometry().pcurve(id_two).cloned(),
        ) else {
            return Ok(None);
        };
        let joined = if id_one == id_two && (range_one.1 - range_two.0).abs() <= tol.parametric() {
            (one, (range_one.0, range_two.1))
        } else {
            // Two chart lines meeting end to end join over the 3D range the
            // joined edge will be read at.
            let (PlanarCurve::Line(_), PlanarCurve::Line(_)) = (&one, &two) else {
                return Ok(None);
            };
            let start = one.point_at(range_one.0, tol)?;
            let end = two.point_at(range_two.1, tol)?;
            if one
                .point_at(range_one.1, tol)?
                .distance(two.point_at(range_two.0, tol)?)
                > tol.confusion()
            {
                return Ok(None);
            }
            (chart_line(start, end, range, tol)?, range)
        };
        pcurves.push(PcurveJoin {
            surface,
            location,
            curve: joined.0,
            range: joined.1,
        });
    }
    // The plan is only a plan until it measures up against the joined curve.
    for pcurve in &pcurves {
        if !agrees(model, &curve, range, pcurve, tol)? {
            return Ok(None);
        }
    }
    Ok(Some(Join {
        curve,
        range,
        pcurves,
    }))
}

/// Whether a joined pcurve, read the way a face reads it (the same
/// fraction along both ranges) lands on the joined curve.
fn agrees(
    model: &Model,
    curve: &ogeom_geom::Curve,
    range: (f64, f64),
    pcurve: &PcurveJoin,
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_geom::Surface as _;
    let Some(surface) = model.geometry().surface(pcurve.surface) else {
        return Ok(false);
    };
    let placement = pcurve.location.composed(model.datums())?;
    for k in 0..=6 {
        let f = f64::from(k) / 6.0;
        let t = (range.1 - range.0).mul_add(f, range.0);
        let pt = (pcurve.range.1 - pcurve.range.0).mul_add(f, pcurve.range.0);
        let (Ok(on_curve), Ok(chart)) = (curve.point_at(t, tol), pcurve.curve.point_at(pt, tol))
        else {
            return Ok(false);
        };
        let Ok(lifted) = surface.point_at(chart.x, chart.y, tol) else {
            return Ok(false);
        };
        if placement.apply(lifted).distance(on_curve) > tol.confusion() {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Shrink every edge and vertex tolerance to what the geometry measures,
/// floored at the confusion tolerance: the inverse of the widenings repair
/// operations apply, safe because it is measured the same way.
///
/// Returns how many claims shrank.
///
/// Nodes below `shape` that other shapes hold as well are copied first
/// ([`Model::unshare`]), so no other shape's claims shrink.
///
/// # Errors
///
/// As evaluation.
pub fn reduce_tolerances(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<usize> {
    model.unshare(shape)?;
    let edges = explore_unique(model, shape, ShapeType::Edge)?;
    // Each edge is measured on its own and the model is only read, so the
    // measures run side by side.
    let measured = {
        let model = &*model;
        ogeom_core::parallel::map_ordered(&edges, |_, edge| measure_edge(model, edge, false, tol))
    };
    let mut shrunk = 0;
    for (edge, measured) in edges.iter().zip(measured) {
        let Some(worst) = measured else { continue };
        let target = ogeom_core::Tolerance::new(worst + tol.confusion())?;
        if let Some(node) = model.node_mut(edge)
            && let NodeData::Edge(data) = node.data_mut()
            && target.get() < data.tolerance.get()
        {
            data.tolerance = target;
            shrunk += 1;
        }
    }
    Ok(shrunk)
}

/// How far an edge's pcurves stand from its curve, or `None` where it has
/// no curve or no pcurve to measure. Each lifted point is held to the
/// curve's point at the same parameter where `paced` or where the edge
/// claims same parameter, and otherwise to the nearest point of the curve.
pub(crate) fn measure_edge(
    model: &Model,
    edge: &Shape,
    paced: bool,
    tol: Tolerances,
) -> Option<f64> {
    let data = model.node(edge).and_then(|n| n.data().as_edge())?;
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        return None;
    };
    let geometry = model.geometry().curve(*curve)?;
    // The claim an edge tolerance covers: every pcurve, both columns
    // of a seam, lifted against its curve, measured as the readers
    // and repairs that widen it measure.
    let mut worst = 0.0f64;
    let mut measurable = false;
    for representation in &data.representations {
        let (sides, prange, surface, location) = match representation {
            EdgeRepr::PCurve {
                curve,
                range,
                surface,
                location,
            } => ([Some(*curve), None], *range, *surface, location),
            EdgeRepr::Seam {
                forward,
                reversed,
                range,
                surface,
                location,
            } => (
                [Some(*forward), Some(*reversed)],
                *range,
                *surface,
                location,
            ),
            _ => continue,
        };
        let Some(surface_geometry) = model.geometry().surface(surface) else {
            continue;
        };
        let placed;
        let surface_geometry = if location.is_identity() {
            surface_geometry
        } else {
            use ogeom_geom::Transformable as _;
            let Ok(placement) = location.composed(model.datums()) else {
                continue;
            };
            let Ok(moved) = surface_geometry.transformed(&placement, tol) else {
                continue;
            };
            placed = moved;
            &placed
        };
        for pc in sides.into_iter().flatten() {
            let Some(pcurve) = model.geometry().pcurve(pc) else {
                continue;
            };
            let Ok(gap) = ogeom_algo::pcurve_fit::lifted_gap(
                (geometry, *range),
                (pcurve, prange),
                surface_geometry,
                paced || data.same_parameter(),
                tol,
            ) else {
                continue;
            };
            worst = worst.max(gap);
            measurable = true;
        }
    }
    measurable.then_some(worst)
}
