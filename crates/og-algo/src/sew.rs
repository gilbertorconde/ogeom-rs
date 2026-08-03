//! Joining loose topology: ordering a bag of edges into a wire, and sewing free
//! faces into a shell.
//!
//! Both exist because geometry arrives disconnected. An imported file gives a
//! pile of faces that *touch* but share nothing; a sketch gives edges in
//! whatever order they were drawn. Topologically these are unrelated pieces,
//! and every algorithm that walks a boundary treats them that way — a shell of
//! faces that merely abut has a free edge everywhere two of them meet, encloses
//! no volume, and cannot be classified against.
//!
//! # Sewing is a topological operation, not a geometric one
//!
//! It does not move anything. Two edges within tolerance of each other are
//! decided to be *one* edge, and every face that used either now uses that one
//! — so the shell closes because the topology says so, not because the geometry
//! was nudged until it did. A version that moved geometry to close gaps would
//! be a repair, would need to decide which of two positions is right, and would
//! quietly invalidate every tolerance in the neighbourhood.
//!
//! What it will not do is claim a closure it did not achieve. Faces that do not
//! meet within tolerance stay in separate shells, and the result says how many
//! there are.

use std::collections::HashMap;

use og_core::{OgResult, Tolerances, og_bail};
use og_geom::Curve3d;
use og_math::Point;
use og_topo::{EdgeRepr, Model, NodeData, Orientation, Shape, ShapeType, TShapeId, explore_unique};

use crate::build::{edge_vertices, make_face_on, make_shell, make_wire};
use crate::history::{Built, History};

/// Roles sewing assigns.
pub mod roles {
    use og_core::Role;

    /// An edge that two faces were found to share.
    pub const SEWN_EDGE: Role = Role::op_defined(40);
    /// A face rebuilt on shared edges.
    pub const SEWN_FACE: Role = Role::op_defined(41);
}

/// Put a bag of edges into an order that walks them end to end.
///
/// Reverses an edge where the chain reaches its far end first, so the result is
/// a path rather than a set. [`make_wire`] then accepts it — it checks that
/// consecutive edges meet, and a bag in the order it happened to be built in
/// almost never does.
///
/// Follows the chain from one end. Where an end meets more than two edges the
/// path is genuinely ambiguous — that is a branching network, not a wire — and
/// this refuses rather than picking one, because picking one silently discards
/// the branch nobody asked it to drop.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if the list is
/// empty, an edge is unbounded, the edges do not form a single connected path,
/// or a vertex joins three or more of them.
pub fn order_edges(model: &Model, edges: &[Shape], tol: Tolerances) -> OgResult<Vec<Shape>> {
    if edges.is_empty() {
        og_bail!(Construction, "there are no edges to order");
    }
    if edges.len() == 1 {
        return Ok(edges.to_vec());
    }

    let mut ends = Vec::with_capacity(edges.len());
    for edge in edges {
        let Some((start, finish)) = edge_vertices(model, edge)? else {
            og_bail!(
                Construction,
                "an unbounded edge cannot be shown to join anything"
            );
        };
        ends.push((placed(model, &start)?, placed(model, &finish)?));
    }

    // How many edge-ends meet at each position. Three is a branch, and a branch
    // has no single walk through it.
    for i in 0..edges.len() {
        for at in [ends[i].0, ends[i].1] {
            let meeting = ends
                .iter()
                .filter(|(a, b)| a.is_equal(at, tol) || b.is_equal(at, tol))
                .count();
            if meeting > 2 {
                og_bail!(
                    Construction,
                    "{meeting} edges meet at {at:?}; that is a branching \
                     network rather than a wire, and choosing a path through it \
                     would silently drop the branches not chosen"
                );
            }
        }
    }

    // Start from a free end if there is one, so an open chain comes out running
    // the way it reads. A closed loop has none, and any edge will do.
    let start = (0..edges.len())
        .find(|&i| {
            !ends.iter().enumerate().any(|(j, (a, b))| {
                j != i && (a.is_equal(ends[i].0, tol) || b.is_equal(ends[i].0, tol))
            })
        })
        .unwrap_or(0);

    let mut used = vec![false; edges.len()];
    let mut out = Vec::with_capacity(edges.len());
    used[start] = true;
    out.push(edges[start].clone());
    let mut reach = ends[start].1;

    while out.len() < edges.len() {
        let mut stepped = false;
        for i in 0..edges.len() {
            if used[i] {
                continue;
            }
            let (a, b) = ends[i];
            if a.is_equal(reach, tol) {
                out.push(edges[i].clone());
                reach = b;
            } else if b.is_equal(reach, tol) {
                // The chain arrived at this edge's far end, so it is walked
                // backwards. Reversing the occurrence is what keeps the wire a
                // path; leaving it would make `make_wire` report a gap that is
                // really a direction.
                out.push(edges[i].reversed());
                reach = a;
            } else {
                continue;
            }
            used[i] = true;
            stepped = true;
            break;
        }
        if !stepped {
            og_bail!(
                Construction,
                "the edges do not form one connected path: {} of {} could not \
                 be reached from the first",
                edges.len() - out.len(),
                edges.len()
            );
        }
    }
    Ok(out)
}

/// Build a wire from edges in any order.
///
/// [`order_edges`] then [`make_wire`].
///
/// # Errors
///
/// As [`order_edges`] and [`make_wire`].
pub fn make_wire_unordered(model: &mut Model, edges: &[Shape], tol: Tolerances) -> OgResult<Built> {
    let ordered = order_edges(model, edges, tol)?;
    make_wire(model, &ordered, tol)
}

/// What sewing produced.
#[derive(Debug, Clone)]
pub struct Sewn {
    /// One shell per connected group of faces.
    ///
    /// More than one means the faces did not all meet. That is reported rather
    /// than papered over: a single shell containing disconnected pieces would
    /// claim a closure that is not there.
    pub shells: Vec<Shape>,
    /// How many pairs of edges were found to be the same edge.
    pub joined: usize,
    /// Edges still used by exactly one face after sewing.
    ///
    /// Zero means every shell is closed. Anything else is the boundary that
    /// remains, and a caller that needs a solid needs this to be empty.
    pub free_edges: Vec<Shape>,
    /// History, as every operation reports.
    pub history: History,
}

/// Sew free faces into shells by finding the edges they share.
///
/// Two edges are the same edge when their ends coincide within tolerance —
/// either way round — *and* a point along them does too. The midpoint test is
/// what stops two different arcs between the same pair of vertices from being
/// merged into one, which is a real case: the two halves of a circle share both
/// ends.
///
/// Nothing is moved. See the module documentation.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if `faces` is
/// empty or holds something that is not a face;
/// [`OgError::Dangling`](og_core::OgError::Dangling) if a handle fails to
/// resolve.
pub fn sew(model: &mut Model, faces: &[Shape], tol: Tolerances) -> OgResult<Sewn> {
    if faces.is_empty() {
        og_bail!(Construction, "there are no faces to sew");
    }
    for face in faces {
        if model.kind_of(face)? != ShapeType::Face {
            og_bail!(Construction, "sewing joins faces");
        }
    }
    model.begin_operation();

    // Vertices first, and this is not an optimisation — it is what makes the
    // rest work. Deciding that two edges are one edge leaves the *neighbouring*
    // edges ending at the vertices they always had, which sit at the same
    // places as the survivor's but are different nodes. A wire built from that
    // mixture is reported to have a gap, because it has one: `is_same_position`
    // asks whether one node appears at two placements, which is the right
    // question and not this one.
    let vertices = merge_vertices(model, faces, tol)?;
    let rebuilt_edges = rebuild_edges(model, faces, &vertices)?;

    // Every distinct edge node used by the faces, with the geometry that
    // decides whether two of them are the same edge.
    let mut catalogue: Vec<(TShapeId, Fingerprint)> = Vec::new();
    for face in faces {
        for edge in explore_unique(model, face, ShapeType::Edge)? {
            let id = rebuilt_edges
                .get(&edge.node())
                .copied()
                .unwrap_or(edge.node());
            if catalogue.iter().any(|(seen, _)| *seen == id) {
                continue;
            }
            if let Some(print) = fingerprint(model, &Shape::of(id), tol)? {
                catalogue.push((id, print));
            }
        }
    }

    // Which node each edge is decided to *be*, and whether it runs the other
    // way from the one it replaced.
    let mut merged: HashMap<TShapeId, (TShapeId, bool)> = HashMap::new();
    let mut joined = 0;
    for i in 0..catalogue.len() {
        if merged.contains_key(&catalogue[i].0) {
            continue;
        }
        for j in (i + 1)..catalogue.len() {
            if merged.contains_key(&catalogue[j].0) {
                continue;
            }
            let Some(flipped) = catalogue[i].1.same_as(&catalogue[j].1, tol) else {
                continue;
            };
            merged.insert(catalogue[j].0, (catalogue[i].0, flipped));
            joined += 1;
        }
    }

    // The survivor has to carry the pcurves of the edge it replaced, or the
    // face that used the replaced one loses its description in parameter space
    // and stops being triangulable.
    for (dropped, (kept, _)) in merged.clone() {
        let carried: Vec<EdgeRepr> = model
            .node_by_id(dropped)
            .and_then(|n| n.data().as_edge())
            .map(|d| {
                d.representations
                    .iter()
                    .filter(|r| r.is_parametric())
                    .cloned()
                    .collect()
            })
            .unwrap_or_default();
        if carried.is_empty() {
            continue;
        }
        let survivor = Shape::of(kept);
        let Some(node) = model.node_mut(&survivor) else {
            og_bail!(Dangling, "an edge is not in this model");
        };
        let NodeData::Edge(data) = node.data_mut() else {
            og_bail!(Construction, "edge node holds no edge data");
        };
        for repr in carried {
            data.add(repr);
        }
    }

    // One map from every original edge to what it is now: rebuilt onto merged
    // vertices, then possibly merged with a coincident twin.
    let mut substitution: HashMap<TShapeId, (TShapeId, bool)> = HashMap::new();
    for (original, rebuilt) in &rebuilt_edges {
        let (final_id, flipped) = merged.get(rebuilt).copied().unwrap_or((*rebuilt, false));
        substitution.insert(*original, (final_id, flipped));
    }
    for (dropped, kept) in &merged {
        substitution.entry(*dropped).or_insert(*kept);
    }

    let mut history = History::new();
    let mut rebuilt = Vec::with_capacity(faces.len());
    for face in faces {
        let sewn = rebuild_face(model, face, &substitution, tol)?;
        model.set_derived(&sewn, std::slice::from_ref(face), roles::SEWN_FACE)?;
        history.modify(face, sewn.clone());
        rebuilt.push(sewn);
    }

    let groups = connected_groups(model, &rebuilt)?;
    let mut shells = Vec::with_capacity(groups.len());
    for group in groups {
        let shell = make_shell(model, &group)?.shape;
        for face in &group {
            history.generate(face, shell.clone());
        }
        shells.push(shell);
    }

    let free_edges = free_edges(model, &rebuilt)?;
    Ok(Sewn {
        shells,
        joined,
        free_edges,
        history,
    })
}

/// Decide which coincident vertices are the same vertex.
///
/// Returns only the ones that were replaced, mapping each to its survivor.
fn merge_vertices(
    model: &Model,
    faces: &[Shape],
    tol: Tolerances,
) -> OgResult<HashMap<TShapeId, TShapeId>> {
    let mut seen: Vec<(TShapeId, Point)> = Vec::new();
    let mut out = HashMap::new();
    for face in faces {
        for vertex in explore_unique(model, face, ShapeType::Vertex)? {
            if seen.iter().any(|(id, _)| *id == vertex.node()) {
                continue;
            }
            let at = placed(model, &vertex)?;
            match seen.iter().find(|(_, p)| p.is_equal(at, tol)) {
                Some((kept, _)) => {
                    out.insert(vertex.node(), *kept);
                }
                None => seen.push((vertex.node(), at)),
            }
        }
    }
    Ok(out)
}

/// Rebuild every edge whose bounding vertices were merged away.
///
/// An edge's bounds live in its node, so an edge cannot be pointed at a
/// different vertex — it has to be built again. Its data comes across whole,
/// representations included, so the new edge describes itself exactly as the
/// old one did and only its ends have changed.
fn rebuild_edges(
    model: &mut Model,
    faces: &[Shape],
    vertices: &HashMap<TShapeId, TShapeId>,
) -> OgResult<HashMap<TShapeId, TShapeId>> {
    let mut out = HashMap::new();
    if vertices.is_empty() {
        return Ok(out);
    }
    let mut done: Vec<TShapeId> = Vec::new();
    for face in faces {
        for edge in explore_unique(model, face, ShapeType::Edge)? {
            if done.contains(&edge.node()) {
                continue;
            }
            done.push(edge.node());

            let Some(node) = model.node(&edge) else {
                og_bail!(Dangling, "edge is not in this model");
            };
            let bounds: Vec<Shape> = node.children().to_vec();
            if !bounds.iter().any(|b| vertices.contains_key(&b.node())) {
                continue;
            }
            let NodeData::Edge(data) = node.data().clone() else {
                continue;
            };
            let moved: Vec<Shape> = bounds
                .iter()
                .map(|b| match vertices.get(&b.node()) {
                    Some(kept) => Shape::new(*kept, b.location().clone(), b.orientation()),
                    None => b.clone(),
                })
                .collect();
            let fresh = model.add_edge(*data, &moved)?;
            out.insert(edge.node(), fresh.node());
        }
    }
    Ok(out)
}

/// What decides whether two edges are the same edge.
#[derive(Debug, Clone, Copy)]
struct Fingerprint {
    start: Point,
    middle: Point,
    end: Point,
}

impl Fingerprint {
    /// Whether two edges coincide, and if so whether the second runs backwards.
    fn same_as(&self, other: &Self, tol: Tolerances) -> Option<bool> {
        // The midpoint is not a nicety. Two arcs between the same pair of
        // vertices — the two halves of a circle — agree at both ends and are
        // not the same edge, and merging them would fuse a shape to itself.
        if !self.middle.is_equal(other.middle, tol) {
            return None;
        }
        if self.start.is_equal(other.start, tol) && self.end.is_equal(other.end, tol) {
            return Some(false);
        }
        if self.start.is_equal(other.end, tol) && self.end.is_equal(other.start, tol) {
            return Some(true);
        }
        None
    }
}

/// An edge's ends and midpoint, in space.
fn fingerprint(model: &Model, edge: &Shape, tol: Tolerances) -> OgResult<Option<Fingerprint>> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        return Ok(None);
    };
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        // A degenerate edge has no curve and no length; there is nothing about
        // it that could match another edge's geometry.
        return Ok(None);
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        og_bail!(Dangling, "curve is not in this model");
    };
    let placement = edge.transform(model.datums())?;
    Ok(Some(Fingerprint {
        start: placement.apply(geometry.point_at(range.0, tol)?),
        middle: placement.apply(geometry.point_at(f64::midpoint(range.0, range.1), tol)?),
        end: placement.apply(geometry.point_at(range.1, tol)?),
    }))
}

/// A vertex's position in space.
fn placed(model: &Model, vertex: &Shape) -> OgResult<Point> {
    let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
        og_bail!(Construction, "expected a vertex");
    };
    Ok(vertex.transform(model.datums())?.apply(data.point))
}

/// Rebuild a face with merged edges in place of the ones they replaced.
fn rebuild_face(
    model: &mut Model,
    face: &Shape,
    merged: &HashMap<TShapeId, (TShapeId, bool)>,
    tol: Tolerances,
) -> OgResult<Shape> {
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()).cloned() else {
        og_bail!(Construction, "expected a face");
    };
    // Nothing to substitute: the face is already built on the shared edges, and
    // rebuilding it would only mint a node identical to the one there.
    let mut touched = false;

    let mut wires = Vec::new();
    for wire in model.ordered_children_of(face)? {
        let mut ring = Vec::new();
        for edge in model.ordered_children_of(&wire)? {
            match merged.get(&edge.node()) {
                Some((kept, flipped)) => {
                    touched = true;
                    let mut replacement =
                        Shape::new(*kept, edge.location().clone(), edge.orientation());
                    if *flipped {
                        replacement = replacement.reversed();
                    }
                    ring.push(replacement);
                }
                None => ring.push(edge),
            }
        }
        wires.push(make_wire(model, &ring, tol)?.shape);
    }
    if !touched {
        return Ok(face.clone());
    }
    let sewn = make_face_on(model, data.surface, &wires, tol)?.shape;
    Ok(if face.orientation() == Orientation::Reversed {
        sewn.reversed()
    } else {
        sewn
    })
}

/// Group faces by whether they share an edge, transitively.
fn connected_groups(model: &Model, faces: &[Shape]) -> OgResult<Vec<Vec<Shape>>> {
    let mut group_of: Vec<usize> = (0..faces.len()).collect();
    let mut edges_of = Vec::with_capacity(faces.len());
    for face in faces {
        edges_of.push(
            explore_unique(model, face, ShapeType::Edge)?
                .into_iter()
                .map(|e| e.node())
                .collect::<Vec<_>>(),
        );
    }

    // Union-find, flattened by hand: the counts here are face counts, so the
    // simple version is not the slow one.
    for i in 0..faces.len() {
        for j in (i + 1)..faces.len() {
            if edges_of[i].iter().any(|e| edges_of[j].contains(e)) {
                let (a, b) = (find(&group_of, i), find(&group_of, j));
                if a != b {
                    group_of[b] = a;
                }
            }
        }
    }

    let mut groups: HashMap<usize, Vec<Shape>> = HashMap::new();
    for (i, face) in faces.iter().enumerate() {
        groups
            .entry(find(&group_of, i))
            .or_default()
            .push(face.clone());
    }
    let mut out: Vec<Vec<Shape>> = groups.into_values().collect();
    // Deterministic: a result whose shells come back in a different order each
    // run is one nobody can compare against.
    out.sort_by_key(|group| group.first().map(Shape::node));
    Ok(out)
}

/// Follow a union-find chain to its root.
fn find(parent: &[usize], mut i: usize) -> usize {
    while parent[i] != i {
        i = parent[i];
    }
    i
}

/// Edges still used by exactly one face.
fn free_edges(model: &Model, faces: &[Shape]) -> OgResult<Vec<Shape>> {
    let mut uses: HashMap<TShapeId, (usize, Shape)> = HashMap::new();
    for face in faces {
        for wire in model.children_of(face)? {
            for edge in model.children_of(&wire)? {
                if model
                    .node(&edge)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| d.degenerate)
                {
                    continue;
                }
                let entry = uses.entry(edge.node()).or_insert((0, edge.clone()));
                entry.0 += 1;
            }
        }
    }
    let mut out: Vec<Shape> = uses
        .into_values()
        .filter(|(count, _)| count % 2 == 1)
        .map(|(_, edge)| edge)
        .collect();
    out.sort_by_key(Shape::node);
    Ok(out)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::{check_tessellation, is_shell_closed, make_box, make_polygon};
    use og_geom::PlaneSurface;
    use og_math::{Frame, Plane, Vector};
    use og_topo::Location;

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> og_mesh::Deflection {
        og_mesh::Deflection {
            chord: 0.02,
            ..og_mesh::Deflection::default()
        }
    }

    /// A square face in the z = `at` plane, built from its own fresh edges so
    /// it shares nothing with anything else.
    fn loose_square(model: &mut Model, corners: [Point; 4]) -> Shape {
        let wire = make_polygon(model, &corners, true, T).unwrap().shape;
        let normal =
            og_math::Direction::from_cross(corners[1] - corners[0], corners[2] - corners[1], T)
                .unwrap();
        let frame = og_math::Frame::new(
            corners[0],
            normal,
            og_math::Direction::new(corners[1] - corners[0], T).unwrap(),
            T,
        )
        .unwrap();
        let surface = model
            .geometry_mut()
            .add_surface(PlaneSurface::new(Plane::new(frame)).into());
        for edge in model.children_of(&wire).unwrap() {
            let (a, b) = crate::edge_vertices(model, &edge).unwrap().unwrap();
            let (pa, pb) = (placed(model, &a).unwrap(), placed(model, &b).unwrap());
            let flat = |p: Point| {
                let l = frame.to_local(p);
                og_math::Point2::new(l.x, l.y)
            };
            crate::attach_pcurve(
                model,
                &edge,
                og_geom::Line2d::segment(flat(pa), flat(pb), T)
                    .unwrap()
                    .into(),
                surface,
                Location::identity(),
                (0.0, pa.distance(pb)),
            )
            .unwrap();
        }
        crate::make_face_on(model, surface, std::slice::from_ref(&wire), T)
            .unwrap()
            .shape
    }

    #[test]
    fn edges_in_any_order_come_back_as_a_path() {
        let mut model = Model::new();
        let corners = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let wire = make_polygon(&mut model, &corners, true, T).unwrap().shape;
        let mut edges = model.children_of(&wire).unwrap();
        // Shuffled, and some of them turned round.
        edges.swap(0, 2);
        edges[1] = edges[1].reversed();
        edges[3] = edges[3].reversed();

        let ordered = order_edges(&model, &edges, T).unwrap();
        assert_eq!(ordered.len(), 4);
        // A wire only builds if consecutive edges actually meet, so this is the
        // property under test rather than a separate one.
        let rebuilt = make_wire(&mut model, &ordered, T).unwrap().shape;
        assert!(crate::is_wire_closed(&model, &rebuilt, T).unwrap());
    }

    #[test]
    fn edges_that_do_not_form_one_path_are_refused() {
        let mut model = Model::new();
        let a = make_polygon(
            &mut model,
            &[Point::ORIGIN, Point::new(1.0, 0.0, 0.0)],
            false,
            T,
        )
        .unwrap()
        .shape;
        let b = make_polygon(
            &mut model,
            &[Point::new(5.0, 0.0, 0.0), Point::new(6.0, 0.0, 0.0)],
            false,
            T,
        )
        .unwrap()
        .shape;
        let mut edges = model.children_of(&a).unwrap();
        edges.extend(model.children_of(&b).unwrap());

        let err = order_edges(&model, &edges, T).unwrap_err();
        assert!(
            err.to_string().contains("connected path"),
            "unexpected message: {err}"
        );
        assert!(order_edges(&model, &[], T).is_err());
    }

    #[test]
    fn a_branching_network_is_refused_rather_than_arbitrarily_walked() {
        // Three edges from one point. Any path through it drops a branch, and
        // dropping one silently is worse than saying there is no answer.
        let mut model = Model::new();
        let hub = Point::ORIGIN;
        let mut edges = Vec::new();
        for tip in [
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
        ] {
            let w = make_polygon(&mut model, &[hub, tip], false, T)
                .unwrap()
                .shape;
            edges.extend(model.children_of(&w).unwrap());
        }
        let err = order_edges(&model, &edges, T).unwrap_err();
        assert!(
            err.to_string().contains("branching"),
            "unexpected message: {err}"
        );
    }

    #[test]
    fn two_faces_that_touch_are_sewn_into_one_shell() {
        let mut model = Model::new();
        let left = loose_square(
            &mut model,
            [
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
                Point::new(1.0, 1.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
            ],
        );
        let right = loose_square(
            &mut model,
            [
                Point::new(1.0, 0.0, 0.0),
                Point::new(2.0, 0.0, 0.0),
                Point::new(2.0, 1.0, 0.0),
                Point::new(1.0, 1.0, 0.0),
            ],
        );

        // Before: eight edges, nothing shared.
        let before = explore_unique(&model, &left, ShapeType::Edge)
            .unwrap()
            .len()
            + explore_unique(&model, &right, ShapeType::Edge)
                .unwrap()
                .len();
        assert_eq!(before, 8);

        let sewn = sew(&mut model, &[left.clone(), right.clone()], T).unwrap();
        assert_eq!(sewn.shells.len(), 1, "they touch, so they are one shell");
        assert_eq!(sewn.joined, 1, "one shared edge");
        assert_eq!(
            explore_unique(&model, &sewn.shells[0], ShapeType::Edge)
                .unwrap()
                .len(),
            7,
            "the shared edge is one edge now, not two"
        );
        // A sheet, so it still has a boundary — six free edges round the
        // outside, and the shared one is not among them.
        assert_eq!(sewn.free_edges.len(), 6);
        assert!(!is_shell_closed(&model, &sewn.shells[0]).unwrap());
        assert!(sewn.history.is_affected(&left));
    }

    #[test]
    fn faces_that_do_not_meet_stay_in_separate_shells() {
        // Claiming one shell would claim a closure that is not there.
        let mut model = Model::new();
        let here = loose_square(
            &mut model,
            [
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
                Point::new(1.0, 1.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
            ],
        );
        let far = loose_square(
            &mut model,
            [
                Point::new(50.0, 0.0, 0.0),
                Point::new(51.0, 0.0, 0.0),
                Point::new(51.0, 1.0, 0.0),
                Point::new(50.0, 1.0, 0.0),
            ],
        );
        let sewn = sew(&mut model, &[here, far], T).unwrap();
        assert_eq!(sewn.shells.len(), 2);
        assert_eq!(sewn.joined, 0);
        assert_eq!(sewn.free_edges.len(), 8);
    }

    #[test]
    fn a_boxs_faces_taken_apart_and_sewn_back_close_again() {
        // The end-to-end case. The faces already share edges here, so what is
        // under test is that sewing does not *break* a shell that was closed —
        // and that the mesh still agrees with the topology afterwards, which is
        // the check a re-built face is most likely to fail.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T)
            .unwrap()
            .shape;
        let faces = explore_unique(&model, &solid, ShapeType::Face).unwrap();

        let sewn = sew(&mut model, &faces, T).unwrap();
        assert_eq!(sewn.shells.len(), 1);
        assert!(sewn.free_edges.is_empty(), "a box has no free edges");
        assert!(is_shell_closed(&model, &sewn.shells[0]).unwrap());
        assert!(
            check_tessellation(&model, &sewn.shells[0], fine(), T)
                .unwrap()
                .is_valid()
        );
    }

    #[test]
    fn two_arcs_between_the_same_vertices_are_not_the_same_edge() {
        // The reason the fingerprint samples the middle. Both halves of a
        // circle agree at both ends; merging them would fuse the shape to
        // itself and the mistake would look like a successful sew.
        let mut model = Model::new();
        let circle = og_math::Circle::new(Frame::WORLD, 1.0, T).unwrap();
        let upper = crate::make_edge(
            &mut model,
            og_geom::CircleCurve::new(circle).into(),
            (0.0, std::f64::consts::PI),
            T,
        )
        .unwrap()
        .shape;
        let lower = crate::make_edge(
            &mut model,
            og_geom::CircleCurve::new(circle).into(),
            (std::f64::consts::PI, std::f64::consts::TAU),
            T,
        )
        .unwrap()
        .shape;

        let a = fingerprint(&model, &upper, T).unwrap().unwrap();
        let b = fingerprint(&model, &lower, T).unwrap().unwrap();
        assert!(
            a.same_as(&b, T).is_none(),
            "two different arcs were called the same edge"
        );
        assert!(a.same_as(&a, T) == Some(false));
    }

    #[test]
    fn an_edge_found_the_other_way_round_is_reversed_rather_than_dropped() {
        let mut model = Model::new();
        let up = Point::new(0.0, 0.0, 1.0);
        let down = Point::new(0.0, 0.0, 0.0);
        let a = crate::make_edge(
            &mut model,
            og_geom::LineCurve::segment(down, up, T).unwrap().into(),
            (0.0, 1.0),
            T,
        )
        .unwrap()
        .shape;
        let b = crate::make_edge(
            &mut model,
            og_geom::LineCurve::segment(up, down, T).unwrap().into(),
            (0.0, 1.0),
            T,
        )
        .unwrap()
        .shape;

        let pa = fingerprint(&model, &a, T).unwrap().unwrap();
        let pb = fingerprint(&model, &b, T).unwrap().unwrap();
        assert_eq!(
            pa.same_as(&pb, T),
            Some(true),
            "the same edge, running the other way"
        );
    }

    #[test]
    fn sewing_nothing_and_sewing_the_wrong_kind_are_refused() {
        let mut model = Model::new();
        assert!(sew(&mut model, &[], T).is_err());
        let vertex = model.add_point(Point::ORIGIN);
        assert!(sew(&mut model, &[vertex], T).is_err());
    }

    #[test]
    fn sewing_does_not_move_anything() {
        // A repair that closes gaps by moving geometry has to decide which of
        // two positions is right, and would invalidate every tolerance nearby.
        // This decides that two edges *are* one edge and leaves the points
        // where they were.
        let mut model = Model::new();
        let solid = make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let faces = explore_unique(&model, &solid, ShapeType::Face).unwrap();
        let before: Vec<Point> = explore_unique(&model, &solid, ShapeType::Vertex)
            .unwrap()
            .iter()
            .map(|v| placed(&model, v).unwrap())
            .collect();

        let sewn = sew(&mut model, &faces, T).unwrap();
        let after: Vec<Point> = explore_unique(&model, &sewn.shells[0], ShapeType::Vertex)
            .unwrap()
            .iter()
            .map(|v| placed(&model, v).unwrap())
            .collect();
        assert_eq!(before.len(), after.len());
        for p in &after {
            assert!(
                before.iter().any(|q| q.is_equal(*p, T)),
                "a vertex moved: {p:?}"
            );
        }
        let _ = Vector::ZERO;
    }
}
