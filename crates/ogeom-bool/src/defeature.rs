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
//! What this does not yet close is refused by name: a band whose sides do
//! not meet in a curve, a removal that would leave a face with no boundary,
//! more than one band at once.

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
        let band = close_band(model, &interrupted, &removed, &users, &is_ring, tol)?;
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

/// Close a band wound: two side faces re-intersected, end faces' edges
/// extended to the recovered corners, all four-plus faces rebuilt.
#[allow(clippy::too_many_lines, reason = "one wound, one narrative")]
fn close_band(
    model: &mut Model,
    interrupted: &[Shape],
    removed: &HashSet<TShapeId>,
    users: &HashMap<TShapeId, Vec<Shape>>,
    is_ring: &dyn Fn(&Shape) -> bool,
    tol: Tolerances,
) -> OgeomResult<Vec<(Shape, Shape)>> {
    // Each interrupted survivor's ring edges, and their total length: the
    // two longest are the band's sides, the rest are its ends.
    let mut with_length: Vec<(Shape, Vec<Shape>, f64)> = Vec::new();
    for face in interrupted {
        let mut ring_edges = Vec::new();
        let mut length = 0.0;
        for edge in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            if is_ring(&edge) {
                length += edge_length(model, &edge, tol)?;
                ring_edges.push(edge);
            }
        }
        with_length.push((face.clone(), ring_edges, length));
    }
    with_length.sort_by(|a, b| b.2.total_cmp(&a.2));
    if with_length.len() < 2 {
        ogeom_bail!(
            Construction,
            "a band wound with a single interrupted neighbour; closing it \
             needs the neighbour to meet itself, which is not constructed yet"
        );
    }
    let (side_a, rims_a, _) = with_length[0].clone();
    let (side_b, rims_b, _) = with_length[1].clone();
    let ends: Vec<(Shape, Vec<Shape>)> = with_length[2..]
        .iter()
        .map(|(f, r, _)| (f.clone(), r.clone()))
        .collect();

    // Where the sides' own surfaces meet is the edge the feature replaced.
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
    let sa = surface_of(model, &side_a)?;
    let sb = surface_of(model, &side_b)?;
    let meeting = intersect_surfaces(&sa, &sb, IntersectOptions::default(), tol)?;
    let SurfaceIntersection::Along(sections) = meeting else {
        ogeom_bail!(
            Construction,
            "the band's side surfaces do not meet along a curve; the edge \
             the feature replaced cannot be recovered from them"
        );
    };
    // The branch that stands in for the band is the one nearest it.
    let anchor = band_anchor(model, removed, users, tol)?;
    let section = sections
        .into_iter()
        .min_by(|p, q| {
            nearest_distance(&p.curve, anchor, tol)
                .total_cmp(&nearest_distance(&q.curve, anchor, tol))
        })
        .ok_or_else(|| ogeom_err!(Construction, "the side surfaces meet along no branch"))?;
    let curve = section.curve;

    // Corners: the recovered curve against each end face's surface. No ends
    // means the band wraps — a rim fillet — and the recovered edge closes.
    let mut corners: Vec<(Shape, f64, Point)> = Vec::new(); // end face, param, point
    for (end, _) in &ends {
        let se = surface_of(model, end)?;
        let hit = intersect_curve_surface(&curve, &se, CurveSurfaceOptions::default(), tol)?;
        let pierce = hit
            .crossings
            .iter()
            .min_by(|p, q| {
                p.point
                    .distance(anchor)
                    .total_cmp(&q.point.distance(anchor))
            })
            .ok_or_else(|| {
                ogeom_err!(
                    Construction,
                    "an end face's surface never meets the recovered edge; \
                     the corner cannot be placed"
                )
            })?;
        corners.push((end.clone(), pierce.on_curve, pierce.point));
    }

    // The new edge: between the two corners for an interrupted band, closed
    // over the curve's own period for a wrapping one.
    let (new_edge, corner_vertices) = match corners.as_slice() {
        [] => {
            let (lo, hi) = curve.domain();
            if !curve.is_periodic() {
                ogeom_bail!(
                    Construction,
                    "a wrapping band recovered an open curve; the closure is \
                     not constructible from it"
                );
            }
            let v = make_vertex(model, curve.point_at(lo, tol)?).shape;
            let edge = make_edge_between(model, curve.clone(), (lo, hi), &v, &v, tol)?.shape;
            (edge, Vec::new())
        }
        [(fa, ta, pa), (fb, tb, pb)] => {
            let va = make_vertex(model, *pa).shape;
            let vb = make_vertex(model, *pb).shape;
            let (t0, t1, v0, v1) = if ta <= tb {
                (*ta, *tb, va.clone(), vb.clone())
            } else {
                (*tb, *ta, vb.clone(), va.clone())
            };
            let edge = make_edge_between(model, curve.clone(), (t0, t1), &v0, &v1, tol)?.shape;
            (edge, vec![(fa.clone(), va), (fb.clone(), vb)])
        }
        more => ogeom_bail!(
            Construction,
            "the band meets {} end faces; only a straight-through band with \
             two ends, or a wrapping band with none, closes today",
            more.len()
        ),
    };

    // Rebuild ends first: extending a cap's dangling edges to the corner
    // populates the shared-extension map, and the sides — which share those
    // very edges — then pick the extended versions up by node. Order is
    // load-bearing, not tidy.
    let mut out = Vec::new();
    let mut extended: HashMap<TShapeId, Shape> = HashMap::new();
    for (face, rims) in &ends {
        let new_face = rebuild_interrupted(
            model,
            face,
            rims,
            &new_edge,
            &corner_vertices,
            &mut extended,
            tol,
        )?;
        out.push((face.clone(), new_face));
    }
    for (face, rims) in [(side_a.clone(), rims_a), (side_b.clone(), rims_b)] {
        let new_face = rebuild_interrupted(
            model,
            &face,
            &rims,
            &new_edge,
            &corner_vertices,
            &mut extended,
            tol,
        )?;
        out.push((face, new_face));
    }
    Ok(out)
}

/// Rebuild one interrupted face: drop its ring edges, extend the edges that
/// now dangle to the band's corner vertices, add the recovered edge where
/// this face borders it, and rechain.
fn rebuild_interrupted(
    model: &mut Model,
    face: &Shape,
    rims: &[Shape],
    new_edge: &Shape,
    corners: &[(Shape, Shape)],
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

    // The corner vertex on this face, when this face is an end.
    let own_corner: Option<&Shape> = corners
        .iter()
        .find(|(f, _)| f.node() == face.node())
        .map(|(_, v)| v);

    let mut wires: Vec<Vec<Shape>> = Vec::new();
    for wire in model.ordered_children_of(face)? {
        let edges = model.ordered_children_of(&wire)?;
        let touched = edges.iter().any(|e| rim_nodes.contains(&e.node()));
        if !touched {
            wires.push(edges);
            continue;
        }
        // Which vertices the dropped rim owned: an edge that shared one now
        // dangles there and must reach the corner instead.
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
            // A face that has already extended this edge — the caps run
            // first — decided for everyone; sewing rejoins on the shared
            // node.
            if let Some(found) = extended.get(&edge.node()) {
                kept.push(found.clone());
                continue;
            }
            let dangles = model
                .ordered_children_of(edge)?
                .iter()
                .any(|v| rim_vertices.contains(&v.node()));
            let e = match (own_corner, dangles) {
                (Some(corner), true) => extend_to_corner(model, edge, corner, extended, tol)?,
                _ => edge.clone(),
            };
            kept.push(e);
        }
        // The sides gain the recovered edge; an end closes on its own
        // extended edges.
        if own_corner.is_none() {
            kept.push(new_edge.clone());
        }
        let chained = ogeom_algo::order_edges(model, &kept, tol)?;
        wires.push(chained);
    }
    Ok(ogeom_algo::make_face_with_pcurves(model, surface, &wires, tol)?.shape)
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

/// A point inside the removed region, to pick branches and corners by: the
/// mean of the removed faces' own vertices.
fn band_anchor(
    model: &Model,
    removed: &HashSet<TShapeId>,
    users: &HashMap<TShapeId, Vec<Shape>>,
    tol: Tolerances,
) -> OgeomResult<Point> {
    let _ = tol;
    let mut sum = ogeom_math::Vector::ZERO;
    let mut n = 0.0;
    for faces in users.values() {
        for face in faces {
            if !removed.contains(&face.node()) {
                continue;
            }
            for vertex in explore(model, face, Filter::OfType(ShapeType::Vertex))? {
                let placement = vertex.transform(model.datums())?;
                if let Some(d) = model.node(&vertex).and_then(|nd| nd.data().as_vertex()) {
                    sum += placement.apply(d.point).to_vector();
                    n += 1.0;
                }
            }
        }
    }
    if n == 0.0 {
        ogeom_bail!(
            Construction,
            "the removed faces carry no vertices to anchor by"
        );
    }
    Ok(Point::ORIGIN + sum * (1.0 / n))
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
