//! Defeaturing by face removal: delete faces, close the wound from the
//! neighbours' own surfaces.
//!
//! The input is a set of faces — what those faces *mean* is the caller's
//! business, and the operation works on a solid whose history is gone. Two
//! wounds exist, and they close differently.
//!
//! A feature whose rim is an **inner loop** of a surviving face — a bore in a
//! lid, a boss on a base, a pocket in the middle of a top — leaves survivors
//! whose boundary is already right except for that loop. The cure is wire
//! surgery: the surviving face is rebuilt without the rim wire, edges,
//! pcurves and all, and nothing is re-intersected because nothing new meets.
//!
//! A feature that **interrupts** its neighbours' outer boundaries — a fillet
//! band or a chamfer along an edge — leaves a gap no surviving boundary
//! closes. The cure is the neighbours themselves: the two side faces'
//! surfaces are re-intersected to recover the edge the blend replaced, the
//! end faces' edges are extended along their own curves to the recovered
//! corners, and the faces are rebuilt on the result. Extension here is the
//! surfaces' and curves' own unbounded carriers — no new geometry is
//! invented, only wider windows of what is already there.
//!
//! Several bands close together. Each removed band recovers its own
//! crease; where two creases meet — two blends that met at a corner, or
//! one blend's flush cap standing against another's band, the cap named
//! with its band — the corner is where one crease pierces the other's
//! side, and it is one vertex for both. What this does not yet close is
//! refused by name: a band whose sides do not meet in a curve, a removal
//! that would leave a face with no boundary, a band whose wound needs a
//! neighbour to meet itself.

use crate::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_algo::{Built, History, make_edge_between, make_solid, make_vertex, sew};
use ogeom_core::ogeom_err;
use ogeom_geom::Curve3d as _;
use ogeom_geom::Transformable as _;
use ogeom_geom::{Curve, SurfaceGeometry};
use ogeom_intersect::{
    CurveSurfaceOptions, IntersectOptions, SurfaceIntersection, intersect_curve_surface,
    intersect_surfaces,
};
use ogeom_math::Point;
use ogeom_topo::{Filter, Model, NodeData, Shape, ShapeType, TShapeId, explore};
use std::collections::{HashMap, HashSet};

/// Remove `faces` from `solid` and close the openings from the neighbours'
/// own geometry.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction), by
/// name, when the removal is not one this operation closes: no face named,
/// every face named, a named shape that is not a face of the solid, a wound
/// whose side surfaces do not meet in a single curve, more than one band, or
/// geometry whose pcurves have no closed form to rebuild with.
pub fn remove_faces(
    model: &mut Model,
    solid: &Shape,
    faces: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Built> {
    if faces.is_empty() {
        ogeom_bail!(Construction, "no faces named; there is nothing to remove");
    }
    // Separate features remove separately. Two bores named in one call are
    // two wounds; classifying their ring edges together declares the two
    // longest interrupted faces "the sides" across both and recovers a
    // nonsense edge. Faces group into features by shared edges, and each
    // feature runs the whole machinery on the previous feature's result —
    // sequential exactly as a caller would have called it, so one call
    // means what N calls mean, in the order given.
    let groups = feature_groups(model, faces)?;
    if groups.len() > 1 {
        let mut current = Built::from_nothing(solid.clone());
        for group in &groups {
            // A later group's faces survive the earlier surgeries untouched
            // — different regions — but the solid they belong to is new.
            let step = remove_faces(model, &current.shape, group, tol)?;
            current = Built {
                shape: step.shape,
                history: current.history.then(&step.history),
            };
        }
        return Ok(current);
    }
    let all_faces = explore(model, solid, Filter::OfType(ShapeType::Face))?;
    let removed: HashSet<TShapeId> = faces.iter().map(Shape::node).collect();
    for face in faces {
        if !all_faces.iter().any(|f| f.node() == face.node()) {
            ogeom_bail!(
                Construction,
                "a face named for removal is not a face of this solid"
            );
        }
    }
    let survivors: Vec<Shape> = all_faces
        .iter()
        .filter(|f| !removed.contains(&f.node()))
        .cloned()
        .collect();
    if survivors.is_empty() {
        ogeom_bail!(
            Construction,
            "every face was named for removal; nothing remains to close"
        );
    }

    // Which edges the removed set shares with the world: an edge is a ring
    // edge when a removed face and a surviving face both use it.
    let mut users: HashMap<TShapeId, Vec<Shape>> = HashMap::new();
    for face in &all_faces {
        for edge in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            users.entry(edge.node()).or_default().push(face.clone());
        }
    }
    let is_ring = |edge: &Shape| -> bool {
        users.get(&edge.node()).is_some_and(|fs| {
            fs.iter().any(|f| removed.contains(&f.node()))
                && fs.iter().any(|f| !removed.contains(&f.node()))
        })
    };

    // Sort survivors: untouched, rim-only (mode A), interrupted (mode B).
    let mut untouched: Vec<Shape> = Vec::new();
    let mut rim_surgery: Vec<(Shape, Vec<Shape>)> = Vec::new(); // face, kept wires
    let mut interrupted: Vec<Shape> = Vec::new();
    for face in &survivors {
        let wires = model.ordered_children_of(face)?;
        let mut kept = Vec::new();
        let mut touched = false;
        let mut partial = false;
        for (index, wire) in wires.iter().enumerate() {
            let edges = model.ordered_children_of(wire)?;
            let ring_count = edges.iter().filter(|e| is_ring(e)).count();
            if ring_count == 0 {
                kept.push(wire.clone());
            } else if ring_count == edges.len() {
                // The whole wire is the feature's rim. Dropping the outer
                // boundary would leave a face with no boundary at all.
                if index == 0 {
                    ogeom_bail!(
                        Construction,
                        "removing these faces erases a neighbour's whole outer \
                         boundary; that face has nothing left to stand on"
                    );
                }
                touched = true;
            } else {
                partial = true;
            }
        }
        if partial {
            interrupted.push(face.clone());
        } else if touched {
            rim_surgery.push((face.clone(), kept));
        } else {
            untouched.push(face.clone());
        }
    }

    let mut history = History::new();
    for face in faces {
        history.delete(face);
    }

    let mut rebuilt: Vec<Shape> = untouched;
    for (face, kept_wires) in rim_surgery {
        let new_face = {
            let Some(data) = model.node(&face).and_then(|n| match n.data() {
                NodeData::Face(d) => Some(d.clone()),
                _ => None,
            }) else {
                ogeom_bail!(Construction, "a surviving face holds no face data");
            };
            // The kept wires carry their edges, and the edges their pcurves
            // for this very surface: nothing to recompute.
            ogeom_algo::make_face_on(model, data.surface, &kept_wires, tol)?.shape
        };
        history.modify(&face, new_face.clone());
        rebuilt.push(new_face);
    }

    if !interrupted.is_empty() {
        let band = close_wound(model, faces, &interrupted, &removed, &users, &is_ring, tol)?;
        for (old, new) in band {
            history.modify(&old, new.clone());
            rebuilt.push(new);
        }
    }

    let sewn = sew(model, &rebuilt, tol)?;
    let [shell] = sewn.shells.as_slice() else {
        ogeom_bail!(
            Construction,
            "closing the wound left {} shells where one solid's worth was \
             expected; the removal disconnected the boundary",
            sewn.shells.len()
        );
    };
    if !ogeom_algo::is_shell_closed(model, shell)? {
        ogeom_bail!(
            Construction,
            "the boundary does not close after removal; the wound needs a \
             closure this operation does not construct yet"
        );
    }
    let built = make_solid(model, std::slice::from_ref(shell))?;
    let mut solid_history = history;
    solid_history.modify(solid, built.shape.clone());
    Ok(Built::new(built.shape, solid_history))
}

/// One crease the wound recovers: the edge a removed band replaced, from
/// its two side faces' own surfaces.
struct Crease {
    /// The side faces, by node, in a fixed order.
    sides: [Shape; 2],
    /// The recovered curve, the branch nearest the removed faces.
    curve: Curve,
    /// The removed faces' extent along the curve.
    extent: (f64, f64),
}

/// Close a wound: each removed band's two side faces re-intersected into
/// the crease it replaced, the creases' ends placed where they pierce the
/// other interrupted faces — or one another's sides, which is where two
/// bands meeting at a corner share their corner — every dangling edge
/// extended along its own curve to the corner standing on it, and every
/// interrupted face rebuilt with the creases it borders.
///
/// A band's sides are the two survivors it shares the most ring length
/// with; a wedge's cap named alongside its band shares the band's sides
/// and folds into the same crease. A crease's ends are the nearest
/// piercings just past the removed faces' own extent along it, so a
/// survivor the curve merely runs through far away is not mistaken for an
/// end. Corners are one vertex wherever two creases place them within
/// tolerance of each other.
#[allow(clippy::too_many_lines, reason = "one wound, one narrative")]
fn close_wound(
    model: &mut Model,
    removed_faces: &[Shape],
    interrupted: &[Shape],
    removed: &HashSet<TShapeId>,
    users: &HashMap<TShapeId, Vec<Shape>>,
    is_ring: &dyn Fn(&Shape) -> bool,
    tol: Tolerances,
) -> OgeomResult<Vec<(Shape, Shape)>> {
    let surface_of = |model: &Model, face: &Shape| -> OgeomResult<SurfaceGeometry> {
        let placement = face.transform(model.datums())?;
        let Some(data) = model.node(face).and_then(|n| n.data().as_face().cloned()) else {
            ogeom_bail!(Construction, "a band face holds no face data");
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            ogeom_bail!(Construction, "a band face's surface is not in this model");
        };
        surface.clone().transformed(&placement, tol)
    };
    let vertices_of = |model: &Model, face: &Shape| -> OgeomResult<Vec<Point>> {
        let mut out = Vec::new();
        for vertex in explore(model, face, Filter::OfType(ShapeType::Vertex))? {
            let placement = vertex.transform(model.datums())?;
            if let Some(d) = model.node(&vertex).and_then(|nd| nd.data().as_vertex()) {
                out.push(placement.apply(d.point));
            }
        }
        Ok(out)
    };
    let interrupted_by_node: HashMap<TShapeId, Shape> =
        interrupted.iter().map(|f| (f.node(), f.clone())).collect();

    // Each removed face's sides: the two interrupted survivors it shares
    // the most ring length with. Creases are keyed by the side pair. A
    // removed face with fewer than two such neighbours — a wedge's cap
    // standing against another blend's band, bordering one wall and two
    // removed faces — joins the crease of a removed neighbour it shares an
    // edge with, once that neighbour has one.
    let mut creases: Vec<(TShapeId, TShapeId, [Shape; 2], Vec<Point>)> = Vec::new();
    let mut crease_of: HashMap<TShapeId, usize> = HashMap::new();
    let mut leftovers: Vec<Shape> = Vec::new();
    for face in removed_faces {
        let mut shared: HashMap<TShapeId, f64> = HashMap::new();
        for edge in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            if !is_ring(&edge) {
                continue;
            }
            let length = edge_length(model, &edge, tol)?;
            for user in users.get(&edge.node()).into_iter().flatten() {
                if !removed.contains(&user.node()) && interrupted_by_node.contains_key(&user.node())
                {
                    *shared.entry(user.node()).or_default() += length;
                }
            }
        }
        let mut ranked: Vec<(TShapeId, f64)> = shared.into_iter().collect();
        ranked.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.index().cmp(&b.0.index())));
        let [(a, _), (b, _), ..] = ranked.as_slice() else {
            leftovers.push(face.clone());
            continue;
        };
        let (lo, hi) = if a.index() <= b.index() {
            (*a, *b)
        } else {
            (*b, *a)
        };
        let points = vertices_of(model, face)?;
        let index = match creases.iter().position(|c| c.0 == lo && c.1 == hi) {
            Some(i) => {
                creases[i].3.extend(points);
                i
            }
            None => {
                creases.push((
                    lo,
                    hi,
                    [
                        interrupted_by_node[&lo].clone(),
                        interrupted_by_node[&hi].clone(),
                    ],
                    points,
                ));
                creases.len() - 1
            }
        };
        crease_of.insert(face.node(), index);
    }
    for face in leftovers {
        let mut joined = None;
        for edge in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            for user in users.get(&edge.node()).into_iter().flatten() {
                if let Some(&index) = crease_of.get(&user.node()) {
                    joined = Some(index);
                }
            }
        }
        let Some(index) = joined else {
            ogeom_bail!(
                Construction,
                "a removed face shares ring edges with fewer than two \
                 interrupted neighbours and borders no removed face with a \
                 crease; closing it needs a neighbour to meet itself, which \
                 is not constructed yet"
            );
        };
        let points = vertices_of(model, &face)?;
        creases[index].3.extend(points);
    }

    // The recovered curve of each crease, and the removed faces' extent on it.
    let mut recovered: Vec<Crease> = Vec::new();
    for (_, _, sides, points) in creases {
        let sa = surface_of(model, &sides[0])?;
        let sb = surface_of(model, &sides[1])?;
        let meeting = intersect_surfaces(&sa, &sb, IntersectOptions::default(), tol)?;
        let SurfaceIntersection::Along(sections) = meeting else {
            ogeom_bail!(
                Construction,
                "the band's side surfaces do not meet along a curve; the edge \
                 the feature replaced cannot be recovered from them"
            );
        };
        let anchor = {
            let mut sum = ogeom_math::Vector::ZERO;
            for p in &points {
                sum += p.to_vector();
            }
            #[allow(clippy::cast_precision_loss)]
            let n = points.len().max(1) as f64;
            Point::ORIGIN + sum * (1.0 / n)
        };
        let section = sections
            .into_iter()
            .min_by(|p, q| {
                nearest_distance(&p.curve, anchor, tol)
                    .total_cmp(&nearest_distance(&q.curve, anchor, tol))
            })
            .ok_or_else(|| ogeom_err!(Construction, "the side surfaces meet along no branch"))?;
        let curve = section.curve;
        let mut extent = (f64::INFINITY, f64::NEG_INFINITY);
        for p in &points {
            let t = parameter_near(&curve, *p, tol)?;
            extent = (extent.0.min(t), extent.1.max(t));
        }
        recovered.push(Crease {
            sides,
            curve,
            extent,
        });
    }

    // Corners: where each crease pierces an interrupted face that is not
    // one of its sides, the nearest piercing past each end of its extent.
    // A shared corner is one vertex.
    let mut corner_vertices: Vec<(Point, Shape)> = Vec::new();
    let mut vertex_at = |model: &mut Model, p: Point| -> Shape {
        if let Some((_, v)) = corner_vertices
            .iter()
            .find(|(q, _)| q.distance(p) <= tol.confusion() * 1e3)
        {
            return v.clone();
        }
        let v = make_vertex(model, p).shape;
        corner_vertices.push((p, v.clone()));
        v
    };
    let mut new_edges: Vec<(Shape, [TShapeId; 2])> = Vec::new(); // edge, its sides
    let mut corners: Vec<Shape> = Vec::new();
    for crease in &recovered {
        let mut piercings: Vec<(f64, Point)> = Vec::new();
        for face in interrupted {
            if crease.sides.iter().any(|s| s.node() == face.node()) {
                continue;
            }
            let se = surface_of(model, face)?;
            let hit =
                intersect_curve_surface(&crease.curve, &se, CurveSurfaceOptions::default(), tol)?;
            for c in &hit.crossings {
                piercings.push((c.on_curve, c.point));
            }
        }
        let slack = tol.parametric().max(1e-6);
        let below = piercings
            .iter()
            .filter(|(t, _)| *t <= crease.extent.0 + slack)
            .max_by(|a, b| a.0.total_cmp(&b.0));
        let above = piercings
            .iter()
            .filter(|(t, _)| *t >= crease.extent.1 - slack)
            .min_by(|a, b| a.0.total_cmp(&b.0));
        match (below, above) {
            (None, None) if piercings.is_empty() => {
                // No ends: the band wraps, and the recovered edge closes.
                let (lo, hi) = crease.curve.domain();
                if !crease.curve.is_periodic() {
                    ogeom_bail!(
                        Construction,
                        "a wrapping band recovered an open curve; the closure is \
                         not constructible from it"
                    );
                }
                let v = make_vertex(model, crease.curve.point_at(lo, tol)?).shape;
                let edge =
                    make_edge_between(model, crease.curve.clone(), (lo, hi), &v, &v, tol)?.shape;
                new_edges.push((edge, [crease.sides[0].node(), crease.sides[1].node()]));
            }
            (Some(&(t0, p0)), Some(&(t1, p1))) => {
                let v0 = vertex_at(model, p0);
                let v1 = vertex_at(model, p1);
                let edge =
                    make_edge_between(model, crease.curve.clone(), (t0, t1), &v0, &v1, tol)?.shape;
                corners.push(v0);
                corners.push(v1);
                new_edges.push((edge, [crease.sides[0].node(), crease.sides[1].node()]));
            }
            _ => ogeom_bail!(
                Construction,
                "an end face's surface never meets the recovered edge; the \
                 corner cannot be placed"
            ),
        }
    }

    let mut out = Vec::new();
    let mut extended: HashMap<TShapeId, Shape> = HashMap::new();
    for face in interrupted {
        let rims: Vec<Shape> = explore(model, face, Filter::OfType(ShapeType::Edge))?
            .into_iter()
            .filter(|e| is_ring(e))
            .collect();
        let borders: Vec<Shape> = new_edges
            .iter()
            .filter(|(_, sides)| sides.contains(&face.node()))
            .map(|(e, _)| e.clone())
            .collect();
        let new_face =
            rebuild_interrupted(model, face, &rims, &borders, &corners, &mut extended, tol)?;
        out.push((face.clone(), new_face));
    }
    Ok(out)
}

/// The connected components of the removal set: faces joined by shared
/// edges belong to one feature and close as one wound.
fn feature_groups(model: &Model, faces: &[Shape]) -> OgeomResult<Vec<Vec<Shape>>> {
    let mut edge_sets: Vec<HashSet<TShapeId>> = Vec::with_capacity(faces.len());
    for face in faces {
        edge_sets.push(
            explore(model, face, Filter::OfType(ShapeType::Edge))?
                .iter()
                .map(Shape::node)
                .collect(),
        );
    }
    let mut group_of: Vec<usize> = (0..faces.len()).collect();
    // Union by scan: small sets, clarity over asymptotics.
    fn root(group_of: &mut [usize], mut i: usize) -> usize {
        while group_of[i] != i {
            group_of[i] = group_of[group_of[i]];
            i = group_of[i];
        }
        i
    }
    for i in 0..faces.len() {
        for j in i + 1..faces.len() {
            if edge_sets[i].intersection(&edge_sets[j]).next().is_some() {
                let (a, b) = (root(&mut group_of, i), root(&mut group_of, j));
                group_of[a.max(b)] = a.min(b);
            }
        }
    }
    let mut groups: HashMap<usize, Vec<Shape>> = HashMap::new();
    for (i, face) in faces.iter().enumerate() {
        groups
            .entry(root(&mut group_of, i))
            .or_default()
            .push(face.clone());
    }
    let mut out: Vec<Vec<Shape>> = groups.into_values().collect();
    // Deterministic order: by each group's smallest node index.
    out.sort_by_key(|g| g.iter().map(|f| f.node().index()).min());
    Ok(out)
}

/// Rebuild one interrupted face: drop its ring edges, extend the edges that
/// now dangle to the corner vertex standing on their own curve, add the
/// recovered edges this face borders, and rechain.
fn rebuild_interrupted(
    model: &mut Model,
    face: &Shape,
    rims: &[Shape],
    borders: &[Shape],
    corners: &[Shape],
    extended: &mut HashMap<TShapeId, Shape>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let placement = face.transform(model.datums())?;
    let Some(data) = model.node(face).and_then(|n| n.data().as_face().cloned()) else {
        ogeom_bail!(Construction, "an interrupted face holds no face data");
    };
    let Some(surface) = model.geometry().surface(data.surface).cloned() else {
        ogeom_bail!(
            Construction,
            "an interrupted face's surface is not in this model"
        );
    };
    let surface = surface.transformed(&placement, tol)?;
    let rim_nodes: HashSet<TShapeId> = rims.iter().map(Shape::node).collect();

    let mut wires: Vec<Vec<Shape>> = Vec::new();
    for wire in model.ordered_children_of(face)? {
        let edges = model.ordered_children_of(&wire)?;
        let touched = edges.iter().any(|e| rim_nodes.contains(&e.node()));
        if !touched {
            wires.push(edges);
            continue;
        }
        // Which vertices the dropped rim owned: an edge that shared one now
        // dangles there and must reach a corner instead.
        let mut rim_vertices: HashSet<TShapeId> = HashSet::new();
        for edge in &edges {
            if rim_nodes.contains(&edge.node()) {
                for v in model.ordered_children_of(edge)? {
                    rim_vertices.insert(v.node());
                }
            }
        }
        let mut kept: Vec<Shape> = Vec::new();
        for edge in &edges {
            if rim_nodes.contains(&edge.node()) {
                continue;
            }
            // A face that has already extended this edge decided for
            // everyone; sewing rejoins on the shared node.
            if let Some(found) = extended.get(&edge.node()) {
                kept.push(found.clone());
                continue;
            }
            let dangles = model
                .ordered_children_of(edge)?
                .iter()
                .any(|v| rim_vertices.contains(&v.node()));
            let e = match (dangles, corner_on_edge(model, edge, corners, tol)?) {
                (true, Some(corner)) => extend_to_corner(model, edge, &corner, extended, tol)?,
                _ => edge.clone(),
            };
            kept.push(e);
        }
        kept.extend(borders.iter().cloned());
        let chained = ogeom_algo::order_edges(model, &kept, tol)?;
        wires.push(chained);
    }
    Ok(ogeom_algo::make_face_with_pcurves(model, surface, &wires, tol)?.shape)
}

/// The corner vertex standing on an edge's own curve, nearest the edge,
/// when one does.
fn corner_on_edge(
    model: &Model,
    edge: &Shape,
    corners: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let placement = edge.transform(model.datums())?;
    let Some((curve, range)) = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .and_then(|d| match d.curve3d()? {
            ogeom_topo::EdgeRepr::Curve3d { curve, range, .. } => Some((*curve, *range)),
            _ => None,
        })
    else {
        return Ok(None);
    };
    let Some(geometry) = model.geometry().curve(curve).cloned() else {
        return Ok(None);
    };
    let geometry = geometry.transformed(&placement, tol)?;
    let head = geometry.point_at(range.0, tol)?;
    let tail = geometry.point_at(range.1, tol)?;
    let mut best: Option<(f64, Shape)> = None;
    for corner in corners {
        let Some(p) = model
            .node(corner)
            .and_then(|n| n.data().as_vertex())
            .map(|d| d.point)
        else {
            continue;
        };
        let t = parameter_near(&geometry, p, tol)?;
        if geometry.point_at(t, tol)?.distance(p) > tol.confusion() * 1e3 {
            continue;
        }
        let gap = head.distance(p).min(tail.distance(p));
        if best.as_ref().is_none_or(|(g, _)| gap < *g) {
            best = Some((gap, corner.clone()));
        }
    }
    Ok(best.map(|(_, c)| c))
}

/// The edge, extended along its own curve so its dangling end reaches the
/// corner vertex — shared across the faces that use it, so sewing rejoins
/// them on one node.
fn extend_to_corner(
    model: &mut Model,
    edge: &Shape,
    corner: &Shape,
    extended: &mut HashMap<TShapeId, Shape>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    if let Some(found) = extended.get(&edge.node()) {
        return Ok(found.clone());
    }
    let placement = edge.transform(model.datums())?;
    let Some((curve, range)) = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .and_then(|d| match d.curve3d()? {
            ogeom_topo::EdgeRepr::Curve3d { curve, range, .. } => Some((*curve, *range)),
            _ => None,
        })
    else {
        ogeom_bail!(Construction, "a dangling edge has no curve to extend");
    };
    let Some(geometry) = model.geometry().curve(curve).cloned() else {
        ogeom_bail!(Construction, "a dangling edge's curve is not in this model");
    };
    let geometry = geometry.transformed(&placement, tol)?;
    let corner_point = model
        .node(corner)
        .and_then(|n| n.data().as_vertex())
        .map(|d| d.point)
        .ok_or_else(|| ogeom_err!(Construction, "a corner vertex holds no point"))?;

    // Which end dangles: the one nearer the corner. The corner's parameter
    // on this curve comes from the geometry the curve already has.
    let head = geometry.point_at(range.0, tol)?;
    let tail = geometry.point_at(range.1, tol)?;
    let t_corner = parameter_near(&geometry, corner_point, tol)?;
    let (vertices, new_range, dangle_head) = {
        // Storage order, deliberately: the range is the stored curve's, and
        // `edge_vertices` would swap the pair for a reversed use.
        let bounds = model.children_of(edge)?;
        let (Some(va), Some(vb)) = (bounds.first().cloned(), bounds.last().cloned()) else {
            ogeom_bail!(Construction, "a dangling edge has no vertices");
        };
        if head.distance(corner_point) <= tail.distance(corner_point) {
            ((corner.clone(), vb), (t_corner, range.1), true)
        } else {
            ((va, corner.clone()), (range.0, t_corner), false)
        }
    };
    let _ = dangle_head;
    if new_range.1 <= new_range.0 {
        ogeom_bail!(
            Construction,
            "extending an edge to its corner inverted its range; the corner \
             sits on the wrong side of the edge"
        );
    }
    // A segment's stored domain ends at its own vertices; the extension is
    // the same line over a wider window.
    let geometry = match geometry {
        Curve::Line(line) => {
            let (lo, hi) = ogeom_geom::Curve3d::domain(&line);
            Curve::Line(ogeom_geom::LineCurve::over(
                line.axis(),
                lo.min(new_range.0),
                hi.max(new_range.1),
            )?)
        }
        other => other,
    };
    let built = make_edge_between(model, geometry, new_range, &vertices.0, &vertices.1, tol)?;
    extended.insert(edge.node(), built.shape.clone());
    Ok(built.shape)
}

/// The corner's parameter on a curve, by closed form where one exists and by
/// projection where not.
fn parameter_near(curve: &Curve, p: Point, tol: Tolerances) -> OgeomResult<f64> {
    match curve {
        Curve::Line(line) => {
            let axis = line.axis();
            Ok((p - axis.location).dot(axis.direction.vector()))
        }
        Curve::Circle(c) => {
            let local = c.circle().frame().to_local(p);
            Ok(local.y.atan2(local.x).rem_euclid(core::f64::consts::TAU))
        }
        _ => Ok(ogeom_algo::project_on_curve(curve, p, 64, tol)?.parameter),
    }
}

fn nearest_distance(curve: &Curve, p: Point, tol: Tolerances) -> f64 {
    ogeom_algo::project_on_curve(curve, p, 32, tol).map_or(f64::INFINITY, |pr| pr.distance)
}

fn edge_length(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<f64> {
    let Some((curve, range)) = model
        .node(edge)
        .and_then(|n| n.data().as_edge())
        .and_then(|d| match d.curve3d()? {
            ogeom_topo::EdgeRepr::Curve3d { curve, range, .. } => Some((*curve, *range)),
            _ => None,
        })
    else {
        return Ok(0.0);
    };
    let Some(geometry) = model.geometry().curve(curve) else {
        return Ok(0.0);
    };
    ogeom_algo::curve_length(geometry, range, tol)
}
