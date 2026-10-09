//! Joining loose topology: ordering a bag of edges into a wire, and sewing free
//! faces into a shell.
//!
//! Both exist because geometry arrives disconnected. An imported file gives a
//! pile of faces that *touch* but share nothing; a sketch gives edges in
//! whatever order they were drawn. Topologically these are unrelated pieces,
//! and every algorithm that walks a boundary treats them that way: a shell of
//! faces that merely abut has a free edge everywhere two of them meet, encloses
//! no volume, and cannot be classified against.
//!
//! # Sewing is a topological operation, not a geometric one
//!
//! It does not move anything. Two edges within tolerance of each other are
//! decided to be *one* edge, and every face that used either uses that one,
//! so the shell closes because the topology says so, not because the geometry
//! was nudged until it did. A version that moved geometry to close gaps would
//! be a repair, would need to decide which of two positions is right, and would
//! quietly invalidate every tolerance in the neighbourhood.
//!
//! What it will not do is claim a closure it did not achieve. Faces that do not
//! meet within tolerance stay in separate shells, and the result says how many
//! there are.

use ogeom_core::{FastMap, FastSet};

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve3d;
use ogeom_math::Point;
use ogeom_topo::{
    EdgeRepr, Filter, Model, NodeData, Orientation, Shape, ShapeType, TShapeId, explore,
    explore_unique,
};

use crate::bins::Bins;
use crate::build::{edge_vertices, make_face_on, make_shell, make_wire};
use crate::history::{Built, History};

/// Roles sewing assigns.
pub mod roles {
    use ogeom_core::Role;

    /// An edge that two faces were found to share.
    pub const SEWN_EDGE: Role = Role::op_defined(40);
    /// A face rebuilt on shared edges.
    pub const SEWN_FACE: Role = Role::op_defined(41);
}

/// Put a bag of edges into an order that walks them end to end.
///
/// Reverses an edge where the chain reaches its far end first, so the result is
/// a path rather than a set. [`make_wire`] then accepts it: it checks that
/// consecutive edges meet, and a bag in the order it happened to be built in
/// almost never does.
///
/// Follows the chain from one end. Where an end meets more than two edges the
/// path is genuinely ambiguous (that is a branching network, not a wire), and
/// this refuses rather than picking one, because picking one silently discards
/// the branch nobody asked it to drop.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the list is
/// empty, an edge is unbounded, the edges do not form a single connected path,
/// or a vertex joins three or more of them.
pub fn order_edges(model: &Model, edges: &[Shape], tol: Tolerances) -> OgeomResult<Vec<Shape>> {
    if edges.is_empty() {
        ogeom_bail!(Construction, "there are no edges to order");
    }
    if edges.len() == 1 {
        return Ok(edges.to_vec());
    }

    let mut ends = Vec::with_capacity(edges.len());
    for edge in edges {
        let Some((start, finish)) = edge_vertices(model, edge)? else {
            ogeom_bail!(
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
                ogeom_bail!(
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
            ogeom_bail!(
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
pub fn make_wire_unordered(
    model: &mut Model,
    edges: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Built> {
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

impl Sewn {
    /// The edges two faces of a shell both walk the same way, in shell
    /// order: none where every face keeps its material on the left of its
    /// rings.
    ///
    /// Sewing turns faces until every edge two of them share is walked once
    /// each way, so an edge named here belongs to a shell that admits no
    /// such orientation (a band with a half twist), whose faces sewing left
    /// as given.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a shell
    /// is not in `model`.
    pub fn edges_walked_one_way(&self, model: &Model) -> OgeomResult<Vec<Shape>> {
        let mut out = Vec::new();
        for shell in &self.shells {
            let mut walks: FastMap<(TShapeId, ogeom_topo::Location), (Shape, usize, usize)> =
                FastMap::default();
            let mut order = Vec::new();
            for face in explore(model, shell, Filter::OfType(ShapeType::Face))? {
                for edge in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
                    if model
                        .node(&edge)
                        .and_then(|n| n.data().as_edge())
                        .is_some_and(|d| d.degenerate)
                    {
                        continue;
                    }
                    let key = (edge.node(), edge.location().clone());
                    let walk = walks.entry(key.clone()).or_insert_with(|| {
                        order.push(key);
                        (edge.clone(), 0, 0)
                    });
                    walk.1 += 1;
                    if edge.orientation() == Orientation::Forward {
                        walk.2 += 1;
                    }
                }
            }
            out.extend(order.into_iter().filter_map(|key| {
                let (edge, count, forward) = walks.remove(&key)?;
                (count == 2 && forward != 1).then_some(edge)
            }));
        }
        Ok(out)
    }
}

/// Sew free faces into shells by finding the edges they share.
///
/// Two edges are the same edge when their ends coincide within tolerance
/// (either way round) *and* a point along them does too. The midpoint test is
/// what stops two different arcs between the same pair of vertices from being
/// merged into one, which is a real case: the two halves of a circle share both
/// ends.
///
/// Each face keeps its material on the left of its rings, so two faces
/// sharing an edge walk it opposite ways. Where two faces of a shell walk
/// a shared edge the same way, the shell's faces are turned over breadth
/// first across their shared edges until each is walked once each way,
/// keeping the orientation most of them already have (open shells too);
/// a closed shell whose faces were turned is then made to face out of the
/// volume it bounds. A closed shell whose faces already agree is kept as
/// given, facing in or out: it may bound a void, which only the solid it
/// joins can say, and [`make_solid`](crate::make_solid) turns it to face
/// out as an outer shell or in as a void. A shell that admits no such
/// orientation (a band with a half twist) keeps its faces as given, and
/// [`Sewn::edges_walked_one_way`] names the edges that show it.
///
/// Nothing is moved. See the module documentation. What the faces share
/// with other shapes is copied before the sewing edits it
/// ([`Model::unshare_each`]); the history runs from the faces given.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if `faces` is
/// empty or holds something that is not a face;
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle fails to
/// resolve.
pub fn sew(model: &mut Model, faces: &[Shape], tol: Tolerances) -> OgeomResult<Sewn> {
    sew_with(model, faces, 0.0, tol)
}

/// Sew faces whose edges meet across a gap of up to `gap`.
///
/// As [`sew`], with two edges taken for one where their ends and middle
/// lie within `gap` of each other rather than within the confusion
/// distance: surfaces built apart (an imported sheet, a fitted fill beside
/// an exact extrusion) meet a few micrometres apart. The edge kept, and the
/// vertices it ends on, widen their tolerances to reach the edge they
/// replace, as healing does, so the shell checks clean. An edge that runs
/// along two or more shorter edges of other faces is first split where
/// their vertices meet it within `gap`, so each piece has a twin.
///
/// # Errors
///
/// As [`sew`], and [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// if `gap` is not finite and non-negative.
pub fn sew_within(
    model: &mut Model,
    faces: &[Shape],
    gap: f64,
    tol: Tolerances,
) -> OgeomResult<Sewn> {
    if !(gap.is_finite() && gap >= 0.0) {
        ogeom_bail!(
            Construction,
            "the gap to sew across, {gap}, is not finite and non-negative"
        );
    }
    sew_with(model, faces, gap, tol)
}

/// Sew faces of which some already share their edges and vertices with
/// every neighbour and are to be passed through as they stand.
///
/// As [`sew`], with each face whose `settled` flag is set taken as already
/// sewn: it is neither compared nor rebuilt, and comes through as an exact
/// copy of itself into the shell it closes with the others. The question is
/// asked only of the other faces, so a large shell with a few faces to sew
/// costs what those few cost. Where sewing the others would move an edge or
/// a vertex a settled face holds, the faces are sewn as [`sew`] sews them,
/// all alike. Unlike [`sew`], nothing the faces share with other shapes is
/// copied first: the edges and vertices kept are edited where they stand.
///
/// # Errors
///
/// As [`sew`], and [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// if `settled` does not hold one flag per face.
pub fn sew_around(
    model: &mut Model,
    faces: &[Shape],
    settled: &[bool],
    tol: Tolerances,
) -> OgeomResult<Sewn> {
    if settled.len() != faces.len() {
        ogeom_bail!(
            Construction,
            "{} settled flags for {} faces; one flag per face",
            settled.len(),
            faces.len()
        );
    }
    check_faces(model, faces)?;
    model.begin_operation();
    if let Some(sewn) = sew_faces(model, faces, settled, 0.0, tol)? {
        return Ok(sewn);
    }
    let none = vec![false; faces.len()];
    match sew_faces(model, faces, &none, 0.0, tol)? {
        Some(sewn) => Ok(sewn),
        None => ogeom_bail!(
            Construction,
            "sewing with no face settled moved a settled face"
        ),
    }
}

fn check_faces(model: &Model, faces: &[Shape]) -> OgeomResult<()> {
    if faces.is_empty() {
        ogeom_bail!(Construction, "there are no faces to sew");
    }
    for face in faces {
        if model.kind_of(face)? != ShapeType::Face {
            ogeom_bail!(Construction, "sewing joins faces");
        }
    }
    Ok(())
}

fn sew_with(model: &mut Model, faces: &[Shape], gap: f64, tol: Tolerances) -> OgeomResult<Sewn> {
    check_faces(model, faces)?;
    model.begin_operation();
    // Sewing widens and describes the edges and vertices it keeps where
    // they stand, so what other shapes hold as well is copied first, a
    // face held whole included; the history runs from the faces given.
    let given = faces;
    let faces = model.unshare_each(given)?.shapes;
    let none = vec![false; faces.len()];
    let Some(mut sewn) = sew_faces(model, &faces, &none, gap, tol)? else {
        ogeom_bail!(
            Construction,
            "sewing with no face settled moved a settled face"
        );
    };
    let mut copied = History::new();
    for (given, face) in given.iter().zip(&faces) {
        if given.node() != face.node() {
            copied.copy(given, face.clone());
        }
    }
    if !copied.is_empty() {
        sewn.history = copied.then(&sewn.history);
    }
    Ok(sewn)
}

/// The sewing itself, the faces flagged `settled` passed through as they
/// stand; `None` where sewing the others would move an edge or a vertex a
/// settled face holds.
fn sew_faces(
    model: &mut Model,
    all: &[Shape],
    settled: &[bool],
    gap: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Sewn>> {
    let reach = gap.max(tol.confusion());
    let originals: Vec<Shape> = all
        .iter()
        .zip(settled)
        .filter(|(_, settled)| !**settled)
        .map(|(face, _)| face.clone())
        .collect();
    // What the settled faces hold, which the sewing of the others must
    // leave where it is.
    let mut held_edges: FastSet<TShapeId> = FastSet::default();
    let mut held_vertices: FastSet<TShapeId> = FastSet::default();
    for (face, _) in all.iter().zip(settled).filter(|(_, settled)| **settled) {
        for edge in explore_unique(model, face, ShapeType::Edge)? {
            held_edges.insert(edge.node());
        }
        for vertex in explore_unique(model, face, ShapeType::Vertex)? {
            held_vertices.insert(vertex.node());
        }
    }

    // Sewing tells edges apart by their nodes, so one node placed twice (a
    // prism's far end edges are its profile's, carried along it by a
    // location) would be read as one edge where it stands first. A face
    // holding such a node is baked: rebuilt with its placements in its
    // geometry, every edge a node of its own where it stands.
    let mut baked = unshared(model, &originals, tol)?;
    if gap > 0.0 {
        baked = split_at_vertices(model, &baked, reach, tol)?;
    }
    let faces = baked.as_slice();

    // Vertices first, and this is not an optimisation; it is what makes the
    // rest work. Deciding that two edges are one edge leaves the *neighbouring*
    // edges ending at the vertices they always had, which sit at the same
    // places as the survivor's but are different nodes. A wire built from that
    // mixture is reported to have a gap, because it has one: `is_same_position`
    // asks whether one node appears at two placements, which is the right
    // question and not this one.
    let mut vertices = merge_vertices(model, faces, reach, tol)?;
    // Two edges decided to be one must end on the same vertices, or the
    // faces that bounded the dropped one keep its neighbours ending where
    // the dropped one did, and their wires open. Twins whose ends are still
    // two vertices (each within its own span of the other edge's end, not
    // of the other vertex) have those vertices joined, and the edges are
    // rebuilt on them and matched again.
    let mut rounds = 0;
    let (rebuilt_edges, merged, joined) = loop {
        let rebuilt_edges = rebuild_edges(model, faces, &vertices)?;

        // Every distinct edge node used by the faces, with the geometry that
        // decides whether two of them are the same edge.
        let mut catalogue: Vec<(TShapeId, Fingerprint)> = Vec::new();
        let mut catalogued: FastSet<TShapeId> = FastSet::default();
        for face in faces {
            for edge in explore_unique(model, face, ShapeType::Edge)? {
                let id = rebuilt_edges
                    .get(&edge.node())
                    .copied()
                    .unwrap_or(edge.node());
                if !catalogued.insert(id) {
                    continue;
                }
                if let Some(mut print) = fingerprint(model, &Shape::of(id), tol)? {
                    // Every comparison honours an edge's width, so the gap
                    // sewn across is carried as one.
                    print.width = print.width.max(reach);
                    print.ends = print.ends.map(|w| w.max(reach));
                    catalogue.push((id, print));
                }
            }
        }

        // Which node each edge is decided to *be*, and whether it runs the other
        // way from the one it replaced.
        let mut merged: FastMap<TShapeId, (TShapeId, bool)> = FastMap::default();
        let mut joined = 0;
        let starts = StartBins::new(&catalogue, tol);
        for i in 0..catalogue.len() {
            if merged.contains_key(&catalogue[i].0) {
                continue;
            }
            for j in starts.twin_candidates(&catalogue, i) {
                if merged.contains_key(&catalogue[j].0) {
                    continue;
                }
                let Some(flipped) = catalogue[i].1.same_as(&catalogue[j].1, tol)? else {
                    continue;
                };
                // The survivor answers for both descriptions of the edge,
                // and its vertices must reach the twin's ends: two fingerprints
                // that matched within their stated widths may still disagree by
                // more than a fresh vertex's tolerance, and the disagreement is
                // recorded where the data model records it.
                let (kept_fp, dropped_fp) = (catalogue[i].1.clone(), catalogue[j].1.clone());
                let survivor = Shape::of(catalogue[i].0);
                let ends = if flipped {
                    [
                        (kept_fp.start, dropped_fp.end),
                        (kept_fp.end, dropped_fp.start),
                    ]
                } else {
                    [
                        (kept_fp.start, dropped_fp.start),
                        (kept_fp.end, dropped_fp.end),
                    ]
                };
                let bounds = model.children_of(&survivor)?;
                for vertex in &bounds {
                    if let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) {
                        let at = data.point;
                        let mut need = data.tolerance.get();
                        for (a, b) in &ends {
                            if at.distance(*a) <= need.max(tol.confusion() * 1e2) {
                                need = need.max(a.distance(*b) + tol.confusion());
                            }
                        }
                        if need > data.tolerance.get()
                            && let Some(node) = model.node_mut(vertex)
                            && let NodeData::Vertex(v) = node.data_mut()
                        {
                            v.tolerance = v.tolerance.widen_to(need);
                        }
                    }
                }
                // Across a gap, the kept edge answers for where the dropped
                // one ran: its tolerance reaches the farthest of the dropped
                // edge's ends and middle from it.
                if gap > 0.0 {
                    let mut need = 0.0_f64;
                    for p in [dropped_fp.start, dropped_fp.middle, dropped_fp.end] {
                        need = need.max(kept_fp.off(p, tol)?);
                    }
                    let need = need + tol.confusion();
                    if let Some(node) = model.node_mut(&survivor)
                        && let NodeData::Edge(e) = node.data_mut()
                        && need > e.tolerance.get()
                    {
                        e.tolerance = e.tolerance.widen_to(need);
                    }
                }
                merged.insert(catalogue[j].0, (catalogue[i].0, flipped));
                joined += 1;
            }
        }
        let apart = twin_ends_apart(model, &merged)?;
        if apart.is_empty() || rounds == 3 {
            break (rebuilt_edges, merged, joined);
        }
        rounds += 1;
        for (gone, keep) in apart {
            join_vertex(model, &mut vertices, gone, keep, tol)?;
        }
    };
    // A settled face keeps its edges and vertices, so none of them may be
    // merged away or rebuilt.
    if vertices.keys().any(|v| held_vertices.contains(v))
        || rebuilt_edges.keys().any(|e| held_edges.contains(e))
        || merged.keys().any(|e| held_edges.contains(e))
    {
        return Ok(None);
    }

    // The survivor has to carry the pcurves of the edge it replaced, or the
    // face that used the replaced one loses its description in parameter space
    // and stops being triangulable.
    //
    // Carrying is not copying. When the merge *flipped* (the two edges run
    // opposite ways), the dropped edge's pcurve traverses the shared points
    // backwards relative to the survivor's own curve, and copied unchanged it
    // makes the survivor's face walk one edge of its boundary the wrong way:
    // the parameter-space ring zigzags to zero area and the face stops being
    // triangulable. The pcurve is consumed by *proportional* same-parameter
    // mapping, so reversing its traversal exactly is swapping the stored
    // range's ends.
    //
    // Nor is it copying when the two edges describe one curve at different
    // paces: a fitted rim against the exact circle it traces, which the
    // match admits by asking each middle to lie on the other's stretch. The
    // dropped edge's pcurve is same-parameter with the *dropped* curve; on
    // the survivor's parameter it drifts along the edge, and a face walks
    // its boundary off the vertex it shares. Such a pcurve is refitted at
    // the survivor's own parameters: the survivor's point at each, found
    // on the dropped curve, read through the pcurve into the chart.
    // In a fixed order: this adds pcurves to the geometry tables and
    // appends to survivors' representations, and a map's order would make
    // both differ from run to run.
    let mut carries: Vec<(TShapeId, (TShapeId, bool))> =
        merged.iter().map(|(d, k)| (*d, *k)).collect();
    carries.sort_by_key(|(dropped, _)| dropped.index());
    for (dropped, (kept, flipped)) in carries {
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
        let paced = repaced_carry(
            model,
            &Shape::of(dropped),
            &survivor,
            flipped,
            &carried,
            reach,
            tol,
        )?;
        let off = carried_off(model, &survivor, &paced, tol)?;
        let Some(node) = model.node_mut(&survivor) else {
            ogeom_bail!(Dangling, "an edge is not in this model");
        };
        let NodeData::Edge(data) = node.data_mut() else {
            ogeom_bail!(Construction, "edge node holds no edge data");
        };
        for repr in paced {
            data.add(repr);
        }
        // Two edges a little apart (a fitted section beside the circle it
        // hugs) are one edge to their ends' tolerances, and the survivor's
        // curve then stands off the twin's pcurves by up to the gap between
        // them. The survivor states that gap, and its vertices hold it,
        // with the smallest tolerance as margin: a fuzzy sew's confusion
        // is its fuzz.
        if let Some(off) = off {
            model.widen(
                &survivor,
                ogeom_core::Tolerance::new(off + ogeom_core::Tolerance::MIN.get())?,
            )?;
        }
    }

    // One map from every original edge to what it becomes: rebuilt onto merged
    // vertices, then possibly merged with a coincident twin.
    let mut substitution: FastMap<TShapeId, (TShapeId, bool)> = FastMap::default();
    for (original, rebuilt) in &rebuilt_edges {
        let (final_id, flipped) = merged.get(rebuilt).copied().unwrap_or((*rebuilt, false));
        substitution.insert(*original, (final_id, flipped));
    }
    for (dropped, kept) in &merged {
        substitution.entry(*dropped).or_insert(*kept);
    }

    let mut history = History::new();
    let mut rebuilt = Vec::with_capacity(all.len());
    // Each rebuilt face's original, and whether it is that original's copy.
    let mut sources: Vec<(Shape, bool)> = Vec::with_capacity(all.len());
    let mut sewing = faces.iter().zip(&originals);
    for (face, settled) in all.iter().zip(settled) {
        if *settled {
            // Its edges and vertices stay as they were: the face is its own
            // copy.
            model.set_derived(face, std::slice::from_ref(face), roles::SEWN_FACE)?;
            rebuilt.push(face.clone());
            sources.push((face.clone(), true));
            continue;
        }
        let Some((face, original)) = sewing.next() else {
            ogeom_bail!(Construction, "a face to sew was lost");
        };
        let Some((sewn, whole)) = rebuild_face(model, face, &substitution, tol)? else {
            history.delete(original);
            continue;
        };
        model.set_derived(&sewn, std::slice::from_ref(original), roles::SEWN_FACE)?;
        // Moved onto the shared edges with every ring kept, the face is the
        // same face: a copy, on its own surface within its own boundary.
        rebuilt.push(sewn);
        sources.push((original.clone(), whole && face == original));
    }

    let groups = connected_groups(model, &rebuilt)?;
    let mut shells = Vec::with_capacity(groups.len());
    for group in groups {
        let members: Vec<Shape> = group.iter().map(|&i| rebuilt[i].clone()).collect();
        let mut shell = make_shell(model, &members)?.shape;
        if orient_group(model, &mut rebuilt, &group, tol)? {
            let members: Vec<Shape> = group.iter().map(|&i| rebuilt[i].clone()).collect();
            shell = make_shell(model, &members)?.shape;
        }
        for &i in &group {
            let (original, copy) = &sources[i];
            // A face turned over to face out of what it bounds is the same
            // face presented the other way: changed, not copied.
            if *copy && rebuilt[i].orientation() == original.orientation() {
                history.copy(original, rebuilt[i].clone());
            } else {
                history.modify(original, rebuilt[i].clone());
            }
            history.generate(&rebuilt[i], shell.clone());
        }
        shells.push(shell);
    }

    let free_edges = free_edges(model, &rebuilt)?;
    Ok(Some(Sewn {
        shells,
        joined,
        free_edges,
        history,
    }))
}

/// Turn the faces of one group (indices into `faces`) so that every edge
/// two of them share is walked once each way; whether any was turned.
///
/// Each face keeps its material on the left of its rings, so two faces
/// sharing an edge agree where they walk it opposite ways, and a face
/// that walks it the way its neighbour does is turned against it. The
/// faces are oriented breadth first over their shared edges from a seed,
/// and each connected run of them keeps the orientation most of its faces
/// already have (its first face's on a tie). A closed group whose faces
/// were turned is then made to face out of what it bounds: where its
/// meshed volume comes out negative, every face is turned. An edge with
/// one face or more than two, a degenerate edge and a seam walked both
/// ways by its one face constrain nothing. Where the faces admit no
/// consistent orientation (a band with a half twist), every face is left
/// as given, and [`Sewn::edges_walked_one_way`] names the edges that say
/// so.
fn orient_group(
    model: &mut Model,
    faces: &mut [Shape],
    group: &[usize],
    tol: Tolerances,
) -> OgeomResult<bool> {
    let uses = walks(model, group.iter().map(|&i| &faces[i]))?;
    if walked_each_way(&uses) {
        return Ok(false);
    }
    let closed = uses.values().all(|(count, _)| *count >= 2);
    // Each edge's walks, by position in the group: which faces walk it and
    // whether forward.
    let mut walkers: FastMap<(TShapeId, ogeom_topo::Location), Vec<(usize, bool)>> =
        FastMap::default();
    for (k, &i) in group.iter().enumerate() {
        for edge in explore(model, &faces[i], Filter::OfType(ShapeType::Edge))? {
            if model
                .node(&edge)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|d| d.degenerate)
            {
                continue;
            }
            walkers
                .entry((edge.node(), edge.location().clone()))
                .or_default()
                .push((k, edge.orientation() == Orientation::Forward));
        }
    }
    // Each face's neighbours across an edge with exactly two faces, and
    // whether the two must differ in sense to walk it opposite ways (they
    // must where they now walk it the same way).
    let mut across: Vec<Vec<(usize, bool)>> = vec![Vec::new(); group.len()];
    for walks in walkers.values() {
        if let [(a, a_forward), (b, b_forward)] = walks[..]
            && a != b
        {
            let differ = a_forward == b_forward;
            across[a].push((b, differ));
            across[b].push((a, differ));
        }
    }
    let mut turn: Vec<Option<bool>> = vec![None; group.len()];
    for seed in 0..group.len() {
        if turn[seed].is_some() {
            continue;
        }
        turn[seed] = Some(false);
        let mut run = vec![seed];
        let mut queue = std::collections::VecDeque::from([seed]);
        while let Some(k) = queue.pop_front() {
            let here = turn[k].unwrap_or(false);
            for &(n, differ) in &across[k] {
                let want = here != differ;
                match turn[n] {
                    None => {
                        turn[n] = Some(want);
                        run.push(n);
                        queue.push_back(n);
                    }
                    Some(has) if has != want => return Ok(false),
                    Some(_) => {}
                }
            }
        }
        // The run keeps the orientation most of its faces have.
        let against = run.iter().filter(|&&k| turn[k] == Some(true)).count();
        if 2 * against > run.len() {
            for &k in &run {
                turn[k] = turn[k].map(|t| !t);
            }
        }
    }
    if turn.iter().all(|t| *t != Some(true)) {
        return Ok(false);
    }
    let mut turned: Vec<Shape> = group
        .iter()
        .zip(&turn)
        .map(|(&i, t)| {
            if *t == Some(true) {
                faces[i].reversed()
            } else {
                faces[i].clone()
            }
        })
        .collect();
    if closed {
        let candidate = model.add_shell(&turned)?;
        if let Ok(mesh) =
            ogeom_mesh::triangulate(model, &candidate, ogeom_mesh::Deflection::default(), tol)
            && mesh.volume() < 0.0
        {
            for face in &mut turned {
                *face = face.reversed();
            }
        }
    }
    let mut changed = false;
    for (&i, face) in group.iter().zip(turned) {
        changed |= face.orientation() != faces[i].orientation();
        faces[i] = face;
    }
    Ok(changed)
}

/// The faces, each holding an edge node that some face places elsewhere
/// too baked into a copy with its placements in its geometry; the others
/// as they are. Edges are matched on their nodes' own geometry, which is
/// sound where every edge stands under one placement (an instance's faces
/// sewn among themselves); where placements differ, a face holding a
/// placed edge is baked too.
fn unshared(model: &mut Model, faces: &[Shape], tol: Tolerances) -> OgeomResult<Vec<Shape>> {
    let mut placements: FastMap<TShapeId, Vec<ogeom_topo::Location>> = FastMap::default();
    let mut held: Vec<Vec<TShapeId>> = Vec::with_capacity(faces.len());
    let mut placed: Vec<bool> = Vec::with_capacity(faces.len());
    let mut seen: Vec<ogeom_topo::Location> = Vec::new();
    for face in faces {
        let mut nodes = Vec::new();
        let mut moved = false;
        for edge in ogeom_topo::explore(model, face, ogeom_topo::Filter::OfType(ShapeType::Edge))? {
            let list = placements.entry(edge.node()).or_default();
            if !list.contains(edge.location()) {
                list.push(edge.location().clone());
            }
            if !seen.contains(edge.location()) {
                seen.push(edge.location().clone());
            }
            moved |= !edge.location().is_identity();
            nodes.push(edge.node());
        }
        held.push(nodes);
        placed.push(moved);
    }
    let mixed = seen.len() > 1;
    let placed_twice: Vec<bool> = held
        .iter()
        .zip(&placed)
        .map(|(nodes, moved)| {
            (mixed && *moved)
                || nodes
                    .iter()
                    .any(|n| placements.get(n).is_some_and(|l| l.len() > 1))
        })
        .collect();
    let to_bake: Vec<Shape> = faces
        .iter()
        .zip(&placed_twice)
        .filter(|(_, twice)| **twice)
        .map(|(face, _)| face.clone())
        .collect();
    if to_bake.is_empty() {
        return Ok(faces.to_vec());
    }
    // Baked together, in a solid of their own for the rebuild to walk, so
    // the edges these faces share stay shared.
    let shell = model.add_shell(&to_bake)?;
    let solid = model.add_solid(&[shell])?;
    let baked = crate::baked_shape(model, &solid, tol)?;
    let mut out = Vec::with_capacity(faces.len());
    for (face, twice) in faces.iter().zip(&placed_twice) {
        let now = if *twice {
            baked.history.modified(face).first().cloned()
        } else {
            None
        };
        out.push(now.unwrap_or_else(|| face.clone()));
    }
    Ok(out)
}

/// The faces with every edge split where a vertex of another edge meets
/// its interior within `reach`, the pieces ending on that vertex; faces
/// with no such edge as they are. Edges and vertices placed by a location
/// are left whole.
fn split_at_vertices(
    model: &mut Model,
    faces: &[Shape],
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let mut vertices: Vec<(TShapeId, Point)> = Vec::new();
    let mut seen: FastSet<TShapeId> = FastSet::default();
    for face in faces {
        for vertex in explore_unique(model, face, ShapeType::Vertex)? {
            if vertex.location().is_identity() && seen.insert(vertex.node()) {
                vertices.push((vertex.node(), placed(model, &vertex)?));
            }
        }
    }
    // Each edge's cuts: the parameter, and the vertex the pieces meet at.
    let mut cuts: Vec<(TShapeId, Vec<(f64, TShapeId)>)> = Vec::new();
    let mut walked: FastSet<TShapeId> = FastSet::default();
    for face in faces {
        for edge in explore_unique(model, face, ShapeType::Edge)? {
            if !edge.location().is_identity() || !walked.insert(edge.node()) {
                continue;
            }
            let Some(print) = fingerprint(model, &edge, tol)? else {
                continue;
            };
            let own: Vec<TShapeId> = model.children_of(&edge)?.iter().map(Shape::node).collect();
            let drawn =
                ogeom_mesh::polyline_of_edge(model, &edge, ogeom_mesh::Deflection::default(), tol)?;
            let (mut lo, mut hi) = (drawn[0], drawn[0]);
            for p in &drawn {
                lo = Point::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
                hi = Point::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
            }
            let (r0, r1) = (
                print.range.0.min(print.range.1),
                print.range.0.max(print.range.1),
            );
            let mut here: Vec<(f64, TShapeId)> = Vec::new();
            for &(v, p) in &vertices {
                let outside = p.x < lo.x - reach
                    || p.y < lo.y - reach
                    || p.z < lo.z - reach
                    || p.x > hi.x + reach
                    || p.y > hi.y + reach
                    || p.z > hi.z + reach;
                if own.contains(&v)
                    || outside
                    || p.distance(print.start) <= reach
                    || p.distance(print.end) <= reach
                {
                    continue;
                }
                let foot = crate::project_on_curve(&print.curve, p, 64, tol)?;
                if foot.distance > reach {
                    continue;
                }
                let mut t = foot.parameter;
                if print.curve.is_periodic() {
                    let (dlo, dhi) = print.curve.domain();
                    let period = dhi - dlo;
                    if period > 0.0 {
                        t = r0 + (t - r0).rem_euclid(period);
                    }
                }
                if t > r0 + tol.parametric() && t < r1 - tol.parametric() {
                    here.push((t, v));
                }
            }
            if !here.is_empty() {
                here.sort_by(|a, b| a.0.total_cmp(&b.0));
                here.dedup_by_key(|c| c.1);
                cuts.push((edge.node(), here));
            }
        }
    }
    if cuts.is_empty() {
        return Ok(faces.to_vec());
    }
    // Each cut edge's pieces, in its own direction.
    let mut pieces: FastMap<TShapeId, Vec<Shape>> = FastMap::default();
    for (node, at) in cuts {
        let edge = Shape::of(node);
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()).cloned() else {
            continue;
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d().cloned() else {
            continue;
        };
        let Some((first, last)) = edge_vertices(model, &edge)? else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(curve).cloned() else {
            continue;
        };
        // The vertex at the curve's start, and the one at its end.
        let start_at = geometry.point_at(range.0, tol)?;
        let (from, to) = if placed(model, &first)?.distance(start_at)
            <= placed(model, &last)?.distance(start_at)
        {
            (first, last)
        } else {
            (last, first)
        };
        let rising = range.1 >= range.0;
        let mut stops: Vec<f64> = at.iter().map(|c| c.0).collect();
        let mut through: Vec<Shape> = at.iter().map(|c| Shape::of(c.1)).collect();
        if !rising {
            stops.reverse();
            through.reverse();
        }
        let mut params = vec![range.0];
        params.extend(stops);
        params.push(range.1);
        let mut ends = vec![from];
        ends.extend(through);
        ends.push(to);
        let span = range.1 - range.0;
        let mut out = Vec::with_capacity(params.len() - 1);
        for k in 0..params.len() - 1 {
            let (a, b) = (params[k], params[k + 1]);
            // Each representation's own stretch for this piece, read
            // proportionally, as same-parameter pcurves are.
            let sub = |r: (f64, f64)| {
                let at = |t: f64| r.0 + (t - range.0) / span * (r.1 - r.0);
                (at(a), at(b))
            };
            let mut piece = data.clone();
            piece.representations.retain(|r| {
                matches!(
                    r,
                    EdgeRepr::Curve3d { .. } | EdgeRepr::PCurve { .. } | EdgeRepr::Seam { .. }
                )
            });
            for r in &mut piece.representations {
                match r {
                    EdgeRepr::Curve3d { range: own, .. } => *own = (a, b),
                    EdgeRepr::PCurve { range: own, .. } | EdgeRepr::Seam { range: own, .. } => {
                        *own = sub(*own);
                    }
                    _ => {}
                }
            }
            // The vertex the piece ends on reaches the curve where the
            // piece does.
            for (vertex, t) in [(&ends[k], a), (&ends[k + 1], b)] {
                let need =
                    placed(model, vertex)?.distance(geometry.point_at(t, tol)?) + tol.confusion();
                if let Some(node) = model.node_mut(vertex)
                    && let NodeData::Vertex(v) = node.data_mut()
                    && need > v.tolerance.get()
                {
                    v.tolerance = v.tolerance.widen_to(need);
                }
            }
            out.push(model.add_edge(piece, &[ends[k].clone(), ends[k + 1].clone()])?);
        }
        pieces.insert(node, out);
    }
    let mut out = Vec::with_capacity(faces.len());
    for face in faces {
        let touched = explore_unique(model, face, ShapeType::Edge)?
            .iter()
            .any(|e| pieces.contains_key(&e.node()));
        let Some(data) = model.node(face).and_then(|n| n.data().as_face()).cloned() else {
            continue;
        };
        if !touched || !face.location().is_identity() {
            out.push(face.clone());
            continue;
        }
        // Read as stored, the face's own sense put back below.
        let stored = face.oriented(Orientation::Forward);
        let mut wires = Vec::new();
        for wire in model.ordered_children_of(&stored)? {
            let mut ring = Vec::new();
            for edge in model.ordered_children_of(&wire)? {
                match pieces.get(&edge.node()) {
                    Some(list) if edge.orientation() == Orientation::Reversed => {
                        ring.extend(list.iter().rev().map(Shape::reversed));
                    }
                    Some(list) => ring.extend(list.iter().cloned()),
                    None => ring.push(edge),
                }
            }
            wires.push(make_wire(model, &ring, tol)?.shape);
        }
        let split = make_face_on(model, data.surface, &wires, tol)?.shape;
        out.push(if face.orientation() == Orientation::Reversed {
            split.reversed()
        } else {
            split
        });
    }
    Ok(out)
}

/// Decide which coincident vertices are the same vertex.
///
/// Returns only the ones that were replaced, mapping each to its survivor.
fn merge_vertices(
    model: &mut Model,
    faces: &[Shape],
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<FastMap<TShapeId, TShapeId>> {
    let tolerance_of = |model: &Model, vertex: &Shape| {
        model
            .node(vertex)
            .and_then(|n| n.data().as_vertex())
            .map_or(0.0, |d| d.tolerance.get())
    };
    // The survivors, binned by position on cells sized to the loosest
    // vertex's reach; a vertex is compared with the survivors within the widest
    // reach any comparison can have, in the order they were kept.
    let mut loosest = reach;
    for face in faces {
        for vertex in explore_unique(model, face, ShapeType::Vertex)? {
            loosest = loosest.max(tolerance_of(model, &vertex));
        }
    }
    let mut bins = Bins::new(loosest);
    let mut seen: Vec<(TShapeId, Point, f64)> = Vec::new();
    let mut index_of: FastMap<TShapeId, usize> = FastMap::default();
    let mut widest = reach;
    let mut out = FastMap::default();
    for face in faces {
        for vertex in explore_unique(model, face, ShapeType::Vertex)? {
            if index_of.contains_key(&vertex.node()) {
                continue;
            }
            let at = placed(model, &vertex)?;
            let own = tolerance_of(model, &vertex);
            // Two vertices are one junction within what their *stated*
            // tolerances allow, not within a fresh vertex's default: a
            // vertex that recorded a welded gap reaches that far, and
            // merging by raw confusion would leave its twin standing a
            // recorded-but-ignored distance away.
            let meets =
                |(_, p, w): &&(TShapeId, Point, f64)| p.distance(at) <= reach.max(*w).max(own);
            let hit = match bins.near(at, widest.max(own)) {
                Some(near) => near.into_iter().map(|i| &seen[i]).find(meets),
                None => seen.iter().find(meets),
            }
            .map(|(kept, p, w)| (*kept, *p, *w));
            match hit {
                Some((kept, p, w)) => {
                    // The survivor answers for the absorbed vertex: its
                    // tolerance widens to reach the absorbed position plus
                    // whatever that vertex itself was allowed to stray.
                    let need = p.distance(at) + own + tol.confusion();
                    if need > w
                        && let Some(node) = model.node_mut(&Shape::of(kept))
                        && let NodeData::Vertex(v) = node.data_mut()
                    {
                        v.tolerance = v.tolerance.widen_to(need);
                    }
                    if let Some(entry) = index_of.get(&kept).map(|&i| &mut seen[i]) {
                        entry.2 = entry.2.max(need);
                        widest = widest.max(entry.2);
                    }
                    out.insert(vertex.node(), kept);
                }
                None => {
                    bins.insert(at, seen.len());
                    index_of.insert(vertex.node(), seen.len());
                    seen.push((vertex.node(), at, own));
                    widest = widest.max(own);
                }
            }
        }
    }
    Ok(out)
}

/// The vertex pairs twin edges end on that are not yet one vertex: the
/// dropped edge's end first, the survivor's second.
fn twin_ends_apart(
    model: &Model,
    merged: &FastMap<TShapeId, (TShapeId, bool)>,
) -> OgeomResult<Vec<(TShapeId, TShapeId)>> {
    let ends = |id: TShapeId| -> Option<(TShapeId, TShapeId)> {
        let children = model.node_by_id(id)?.children();
        Some((children.first()?.node(), children.last()?.node()))
    };
    let mut out = Vec::new();
    let mut pairs: Vec<(&TShapeId, &(TShapeId, bool))> = merged.iter().collect();
    pairs.sort_by_key(|(dropped, _)| dropped.index());
    for (dropped, (kept, flipped)) in pairs {
        let (Some((d0, d1)), Some((k0, k1))) = (ends(*dropped), ends(*kept)) else {
            continue;
        };
        let matched = if *flipped {
            [(d0, k1), (d1, k0)]
        } else {
            [(d0, k0), (d1, k1)]
        };
        for (d, k) in matched {
            if d != k && !out.contains(&(d, k)) {
                out.push((d, k));
            }
        }
    }
    Ok(out)
}

/// Make `gone` one vertex with `keep`: every vertex mapped to either then
/// maps to `keep`'s survivor, whose tolerance reaches `gone`'s span.
fn join_vertex(
    model: &mut Model,
    vertices: &mut FastMap<TShapeId, TShapeId>,
    gone: TShapeId,
    keep: TShapeId,
    tol: Tolerances,
) -> OgeomResult<()> {
    let resolve = |vertices: &FastMap<TShapeId, TShapeId>, mut v: TShapeId| {
        while let Some(&next) = vertices.get(&v) {
            if next == v {
                break;
            }
            v = next;
        }
        v
    };
    let (gone, keep) = (resolve(vertices, gone), resolve(vertices, keep));
    if gone == keep {
        return Ok(());
    }
    let (Some(g), Some(k)) = (
        model
            .node_by_id(gone)
            .and_then(|n| n.data().as_vertex())
            .map(|d| (d.point, d.tolerance.get())),
        model
            .node_by_id(keep)
            .and_then(|n| n.data().as_vertex())
            .map(|d| d.point),
    ) else {
        return Ok(());
    };
    let need = g.0.distance(k) + g.1 + tol.confusion();
    if let Some(node) = model.node_mut(&Shape::of(keep))
        && let NodeData::Vertex(v) = node.data_mut()
    {
        v.tolerance = v.tolerance.widen_to(need);
    }
    for target in vertices.values_mut() {
        if *target == gone {
            *target = keep;
        }
    }
    vertices.insert(gone, keep);
    Ok(())
}

/// The dropped edge's parametric representations as the survivor carries
/// them: reversed when the merge flipped, and refitted at the survivor's
/// parameters when the two curves pace one stretch differently. Two edges
/// on one curve object, or on curves whose middle parameters land within
/// the pair's honesty of each other, are carried as they are.
fn repaced_carry(
    model: &mut Model,
    dropped: &Shape,
    survivor: &Shape,
    flipped: bool,
    carried: &[EdgeRepr],
    floor: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<EdgeRepr>> {
    let as_is = || -> Vec<EdgeRepr> {
        carried
            .iter()
            .cloned()
            .map(|r| if flipped { reversed_repr(r) } else { r })
            .collect()
    };
    let (Some(dropped_fp), Some(kept_fp)) = (
        fingerprint(&*model, dropped, tol)?,
        fingerprint(&*model, survivor, tol)?,
    ) else {
        return Ok(as_is());
    };
    let reach = floor.max(dropped_fp.width).max(kept_fp.width);
    if dropped_fp.middle.distance(kept_fp.middle) <= reach {
        if std::env::var_os("OGEOM_DEBUG_SEW").is_some() {
            eprintln!(
                "SEW carry as is: flipped {flipped}, middles {:.2e} apart, kept {:?} {:?} dropped {:?} {:?}",
                dropped_fp.middle.distance(kept_fp.middle),
                kept_fp.start,
                kept_fp.end,
                dropped_fp.start,
                dropped_fp.end
            );
        }
        return Ok(as_is());
    }
    const SAMPLES: usize = 24;
    let (klo, khi) = (kept_fp.range.0, kept_fp.range.1);
    let mut out = Vec::with_capacity(carried.len());
    let mut pending: Vec<(
        ogeom_topo::SurfaceId,
        ogeom_topo::Location,
        ogeom_geom::PlanarCurve,
        (f64, f64),
    )> = Vec::new();
    for repr in carried {
        let EdgeRepr::PCurve {
            curve: pc_id,
            surface,
            location,
            range: prange,
        } = repr
        else {
            out.push(if flipped {
                reversed_repr(repr.clone())
            } else {
                repr.clone()
            });
            continue;
        };
        let Some(pcurve) = model.geometry().pcurve(*pc_id) else {
            ogeom_bail!(Dangling, "pcurve is not in this model");
        };
        let periods = model.geometry().surface(*surface).map(|sg| {
            use ogeom_geom::Surface as _;
            let ((ua, ub), (va, vb)) = sg.domain();
            (
                if sg.is_periodic_u() { ub - ua } else { 0.0 },
                if sg.is_periodic_v() { vb - va } else { 0.0 },
            )
        });
        let (dlo, dhi) = (dropped_fp.range.0, dropped_fp.range.1);
        let mut params = Vec::with_capacity(SAMPLES + 1);
        let mut image: Vec<ogeom_math::Point2> = Vec::with_capacity(SAMPLES + 1);
        for k in 0..=SAMPLES {
            #[allow(clippy::cast_precision_loss)]
            let t = klo + (khi - klo) * (k as f64) / (SAMPLES as f64);
            let p = kept_fp.curve.point_at(t, tol)?;
            let foot = crate::project_on_curve(&dropped_fp.curve, p, 64, tol)?;
            // The foot's parameter on the dropped curve, then through the
            // proportional map onto the pcurve's own window. On a curve that
            // closes on itself the foot may come back a turn away from the
            // stretch (a rim's last piece, seen from its own points, sits
            // at the start of the loop as much as at its end) and is
            // carried across the turn before it is clamped.
            let (lo, hi) = (dlo.min(dhi), dlo.max(dhi));
            let (da, db) = dropped_fp.curve.domain();
            let turn = db - da;
            let mut s = foot.parameter;
            if turn > 0.0 {
                if s < lo - tol.parametric() && s + turn <= hi + tol.parametric() {
                    s += turn;
                } else if s > hi + tol.parametric() && s - turn >= lo - tol.parametric() {
                    s -= turn;
                }
            }
            let s = s.clamp(lo, hi);
            let pt = if (dhi - dlo).abs() <= f64::MIN_POSITIVE {
                prange.0
            } else {
                prange.0 + (prange.1 - prange.0) * (s - dlo) / (dhi - dlo)
            };
            let mut uv = ogeom_geom::Curve2d::point_at(pcurve, pt, tol)?;
            if let (Some((pu, pv)), Some(prev)) = (periods, image.last()) {
                for (coord, period, before) in [(&mut uv.x, pu, prev.x), (&mut uv.y, pv, prev.y)] {
                    if period > 0.0 {
                        while *coord - before > period / 2.0 {
                            *coord -= period;
                        }
                        while before - *coord > period / 2.0 {
                            *coord += period;
                        }
                    }
                }
            }
            params.push(t);
            image.push(uv);
        }
        let fitted =
            ogeom_geom::fit::fit_points_2d_at(&params, &image, 3, tol.confusion() * 10.0, tol)?;
        if std::env::var_os("OGEOM_DEBUG_SEW").is_some() {
            use ogeom_geom::Surface as _;
            let mut worst_sample = 0.0_f64;
            let mut worst_fit = 0.0_f64;
            if let Some(sg) = model.geometry().surface(*surface) {
                for (t, uv) in params.iter().zip(&image) {
                    let p = kept_fp.curve.point_at(*t, tol)?;
                    worst_sample = worst_sample.max(sg.point_at(uv.x, uv.y, tol)?.distance(p));
                    let f = ogeom_geom::Curve2d::point_at(&fitted.curve, *t, tol)?;
                    worst_fit = worst_fit.max(sg.point_at(f.x, f.y, tol)?.distance(p));
                }
            }
            eprintln!(
                "SEW refit: samples off surface by {worst_sample:.2e}, fit off by {worst_fit:.2e} (fit error {:.2e} met {}) domain {:?} over ({klo:.4}, {khi:.4}); dropped range {:?}",
                fitted.error,
                fitted.met,
                ogeom_geom::Curve2d::domain(&fitted.curve),
                dropped_fp.range
            );
        }
        if !fitted.met || fitted.error > reach.max(tol.confusion() * 1e3) {
            // A refit that misses its budget is no description of the edge;
            // the pcurve is carried as it came, its drift along the edge and
            // all, rather than replaced by a worse one.
            out.push(if flipped {
                reversed_repr(repr.clone())
            } else {
                repr.clone()
            });
            continue;
        }
        pending.push((*surface, location.clone(), fitted.curve.into(), (klo, khi)));
    }
    for (surface, location, planar, range) in pending {
        let curve = model.geometry_mut().add_pcurve(planar);
        out.push(EdgeRepr::PCurve {
            curve,
            surface,
            location,
            range,
        });
    }
    Ok(out)
}

/// How far the pcurves carried onto `survivor` leave its curve beyond its
/// stated tolerance, each lifted through its surface at 33 points along it
/// and measured against the nearest point of the curve's stretch; `None`
/// where they keep within it.
fn carried_off(
    model: &Model,
    survivor: &Shape,
    carried: &[EdgeRepr],
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let Some(data) = model.node(survivor).and_then(|n| n.data().as_edge()) else {
        ogeom_bail!(Dangling, "an edge is not in this model");
    };
    let Some(EdgeRepr::Curve3d {
        curve,
        range,
        location,
    }) = data.curve3d()
    else {
        return Ok(None);
    };
    let Some(curve) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    // The edge's own statement, not the confusion: a fuzzy sew's
    // confusion is its fuzz, and a pcurve standing off within it still
    // stands off the curve.
    let reach = data.tolerance.get();
    let mut widest: Option<f64> = None;
    for repr in carried {
        let (sides, prange, surface, at) = match repr {
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
        if at != location {
            continue;
        }
        let Some(surface) = model.geometry().surface(surface) else {
            ogeom_bail!(Dangling, "surface is not in this model");
        };
        for id in sides.into_iter().flatten() {
            let Some(pcurve) = model.geometry().pcurve(id) else {
                ogeom_bail!(Dangling, "pcurve is not in this model");
            };
            // A pcurve fitted through its edge's points strays most between
            // them, so the samples are many more than a fit takes per span.
            const SAMPLES: u32 = 32;
            for i in 0..=SAMPLES {
                let t = f64::from(i) / f64::from(SAMPLES);
                let on_curve = curve.point_at(range.0 + (range.1 - range.0) * t, tol)?;
                let uv = ogeom_geom::Curve2d::point_at(
                    pcurve,
                    prange.0 + (prange.1 - prange.0) * t,
                    tol,
                )?;
                let Ok(lifted) = ogeom_geom::Surface::point_at(surface, uv.x, uv.y, tol) else {
                    continue;
                };
                let mut gap = on_curve.distance(lifted);
                if gap > reach {
                    gap = gap.min(crate::pcurve_gap::nearest_on_stretch(
                        curve, *range, lifted, tol,
                    )?);
                }
                if gap > reach {
                    widest = Some(widest.map_or(gap, |w| w.max(gap)));
                }
            }
        }
    }
    Ok(widest)
}

/// A parametric representation running the other way.
///
/// The consumers map a 3D-curve parameter onto the stored range
/// proportionally, so swapping the range's ends reverses the traversal
/// exactly, with no new geometry. A seam also swaps which pcurve is the
/// forward one, since "forward" is defined by the traversal that just
/// reversed.
fn reversed_repr(repr: EdgeRepr) -> EdgeRepr {
    match repr {
        EdgeRepr::PCurve {
            curve,
            surface,
            location,
            range,
        } => EdgeRepr::PCurve {
            curve,
            surface,
            location,
            range: (range.1, range.0),
        },
        EdgeRepr::Seam {
            forward,
            reversed,
            surface,
            location,
            range,
        } => EdgeRepr::Seam {
            forward: reversed,
            reversed: forward,
            surface,
            location,
            range: (range.1, range.0),
        },
        other => other,
    }
}

/// Rebuild every edge whose bounding vertices were merged away.
///
/// An edge's bounds live in its node, so an edge cannot be pointed at a
/// different vertex; it has to be built again. Its data comes across whole,
/// representations included, so the new edge describes itself exactly as the
/// old one did and only its ends have changed.
fn rebuild_edges(
    model: &mut Model,
    faces: &[Shape],
    vertices: &FastMap<TShapeId, TShapeId>,
) -> OgeomResult<FastMap<TShapeId, TShapeId>> {
    let mut out = FastMap::default();
    if vertices.is_empty() {
        return Ok(out);
    }
    let mut done: FastSet<TShapeId> = FastSet::default();
    for face in faces {
        for edge in explore_unique(model, face, ShapeType::Edge)? {
            if !done.insert(edge.node()) {
                continue;
            }

            let Some(node) = model.node(&edge) else {
                ogeom_bail!(Dangling, "edge is not in this model");
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
#[derive(Debug, Clone)]
struct Fingerprint {
    start: Point,
    middle: Point,
    end: Point,
    /// The points a quarter and three quarters of the way along, which tell
    /// the sense of a closed edge whose ends are one point.
    quarters: [Point; 2],
    /// The edge's curve in space and the stretch it covers, for the middle
    /// of another edge to be asked whether it lies on this one.
    curve: ogeom_geom::Curve,
    range: (f64, f64),
    /// How far this edge's own stated tolerances let it stray: the widest of
    /// the edge's and its vertices'. An edge whose junction was welded across
    /// a recorded gap carries that gap here, and the comparison honours it:
    /// per-entity tolerances are the data model's, not a nicety of import.
    width: f64,
    /// How far each end may stray: the edge's own tolerance and that end's
    /// vertex's. A vertex widened where it was welded claims its span there
    /// and not at the edge's other end.
    ends: [f64; 2],
}

impl Fingerprint {
    /// How far `p` sits from this edge's own stretch of its curve.
    fn off(&self, p: Point, tol: Tolerances) -> OgeomResult<f64> {
        // A line's or a circle's foot is had in closed form; a point on a
        // circle's axis, every angle as near, is left to the search.
        let closed = match &self.curve {
            ogeom_geom::Curve::Line(line) => {
                let t = ogeom_math::elementary::line_parameter(line.axis(), p);
                let at = ogeom_math::elementary::line_at(line.axis(), t).point;
                Some((t, p.distance(at)))
            }
            ogeom_geom::Curve::Circle(circle) => {
                match ogeom_math::elementary::circle_parameter(&circle.circle(), p, tol) {
                    Ok(angle) => {
                        let t = if circle.is_reversed() {
                            (-angle).rem_euclid(core::f64::consts::TAU)
                        } else {
                            angle
                        };
                        let at = ogeom_math::elementary::circle_at(&circle.circle(), angle).point;
                        Some((t, p.distance(at)))
                    }
                    Err(_) => None,
                }
            }
            _ => None,
        };
        let (parameter, distance) = match closed {
            Some(foot) => foot,
            None => {
                let foot = crate::project_on_curve(&self.curve, p, 64, tol)?;
                (foot.parameter, foot.distance)
            }
        };
        let (lo, hi) = (
            self.range.0.min(self.range.1),
            self.range.0.max(self.range.1),
        );
        let mut t = parameter;
        if self.curve.is_periodic() {
            let (dlo, dhi) = self.curve.domain();
            let period = dhi - dlo;
            if period > 0.0 {
                t = lo + (t - lo).rem_euclid(period);
            }
        }
        if t >= lo - tol.parametric() && t <= hi + tol.parametric() {
            return Ok(distance);
        }
        Ok(p.distance(self.start).min(p.distance(self.end)))
    }

    /// Whether two edges coincide, and if so whether the second runs backwards.
    fn same_as(&self, other: &Self, tol: Tolerances) -> OgeomResult<Option<bool>> {
        let reach = tol.confusion().max(self.width).max(other.width);
        let near = |a: Point, b: Point| a.distance(b) <= reach;
        // Each pair of ends within what those two ends state: a vertex
        // widened where a section was welded into a junction claims that
        // span at the junction, and two edges of a sliver that meet there
        // stay apart at their other ends by however much they part.
        let meet =
            |a: Point, wa: f64, b: Point, wb: f64| a.distance(b) <= tol.confusion().max(wa).max(wb);
        let [self_start, self_end] = self.ends;
        let [other_start, other_end] = other.ends;
        // Ends first: they cost a distance each, and most candidates that
        // start near this edge end somewhere else. The middle is asked only
        // of a pair whose ends already agree.
        let along = meet(self.start, self_start, other.start, other_start)
            && meet(self.end, self_end, other.end, other_end);
        let against = meet(self.start, self_start, other.end, other_end)
            && meet(self.end, self_end, other.start, other_start);
        if !along && !against {
            return Ok(None);
        }
        // The midpoint is not a nicety. Two arcs between the same pair of
        // vertices (the two halves of a circle) agree at both ends and are
        // not the same edge, and merging them would fuse a shape to itself.
        // Two descriptions of one curve need not agree on where its middle
        // *parameter* falls (a fitted rim against the exact circle it
        // traces paces itself differently), so each middle is asked to lie
        // on the other's stretch instead, which the far half of a circle
        // still fails.
        if !near(self.middle, other.middle)
            && !(other.off(self.middle, tol)? <= reach && self.off(other.middle, tol)? <= reach)
        {
            if std::env::var_os("OGEOM_DEBUG_SEW").is_some() {
                let foot = crate::project_on_curve(&other.curve, self.middle, 64, tol)?;
                eprintln!(
                    "SEW near miss: ends agree within {reach:.2e}, middles off {:.2e} / {:.2e} (widths {:.2e}, {:.2e}) at {:?}; foot on other at {:.6} (range {:?}, domain {:?}, periodic {}) distance {:.2e}",
                    other.off(self.middle, tol)?,
                    self.off(other.middle, tol)?,
                    self.width,
                    other.width,
                    self.middle,
                    foot.parameter,
                    other.range,
                    other.curve.domain(),
                    other.curve.is_periodic(),
                    foot.distance
                );
            }
            return Ok(None);
        }
        // A closed edge meets the other's ends both ways round; its sense is
        // the one whose quarter points lie the nearer.
        if along && against {
            let [q1, q3] = self.quarters;
            let [p1, p3] = other.quarters;
            let forward = q1.distance(p1) + q3.distance(p3);
            let backward = q1.distance(p3) + q3.distance(p1);
            return Ok(Some(backward < forward));
        }
        Ok(Some(!along))
    }
}

/// Edge starts binned for the twin search. Two edges are one only when
/// each end of one meets an end of the other, so an edge's twin starts
/// within reach of one of its ends, and the edges binned near either end
/// are every edge that could match. A pair is compared within the wider
/// of its two edges' widths, so the few loose edges an import can carry
/// are binned apart: the tight edges' cells are sized by the tight edges,
/// and one loose edge does not coarsen the search for all of them.
struct StartBins {
    /// Edges no wider than `cell`.
    tight: Bins,
    cell: f64,
    /// Edges wider than `cell`, in cells sized to the widest.
    loose: Bins,
    widest: f64,
}

impl StartBins {
    fn new(catalogue: &[(TShapeId, Fingerprint)], tol: Tolerances) -> Self {
        let mut widths: Vec<f64> = catalogue
            .iter()
            .map(|(_, print)| print.width.max(tol.confusion()))
            .collect();
        widths.sort_by(f64::total_cmp);
        let widest = widths.last().copied().unwrap_or(tol.confusion());
        // The width nine in ten edges keep within.
        let cell = widths
            .get(widths.len().saturating_sub(1) * 9 / 10)
            .copied()
            .unwrap_or(tol.confusion());
        let (mut tight, mut loose) = (Bins::new(cell), Bins::new(widest));
        for (index, (_, print)) in catalogue.iter().enumerate() {
            if print.width <= cell {
                tight.insert(print.start, index);
            } else {
                loose.insert(print.start, index);
            }
        }
        Self {
            tight,
            cell,
            loose,
            widest,
        }
    }

    /// The edges after `i` whose start lies within reach of either of its
    /// ends: every edge that could be its twin, and perhaps some that are
    /// not, in order.
    fn twin_candidates(&self, catalogue: &[(TShapeId, Fingerprint)], i: usize) -> Vec<usize> {
        let print = &catalogue[i].1;
        // A tight edge is compared within this edge's width or its own,
        // which is at most `cell`; a loose one within at most the widest.
        let near_tight = self.cell.max(print.width);
        let found = [
            self.tight.near(print.start, near_tight),
            self.tight.near(print.end, near_tight),
            self.loose.near(print.start, self.widest),
            self.loose.near(print.end, self.widest),
        ];
        // A reach wider than the grid answers for is walked: every later
        // edge whose start lies within the widest reach of either end.
        if found.iter().any(Option::is_none) {
            let reach = self.widest.max(near_tight);
            return ((i + 1)..catalogue.len())
                .filter(|&j| {
                    let start = catalogue[j].1.start;
                    start.distance(print.start) <= reach || start.distance(print.end) <= reach
                })
                .collect();
        }
        let mut out: Vec<usize> = found
            .into_iter()
            .flatten()
            .flatten()
            .filter(|&j| j > i)
            .collect();
        out.sort_unstable();
        out.dedup();
        out
    }
}

/// An edge's ends and midpoint, in space.
fn fingerprint(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<Option<Fingerprint>> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        return Ok(None);
    };
    let Some(EdgeRepr::Curve3d {
        curve,
        range,
        location,
    }) = data.curve3d()
    else {
        // A degenerate edge has no curve and no length; there is nothing about
        // it that could match another edge's geometry.
        return Ok(None);
    };
    let Some(geometry) = model.geometry().curve(*curve) else {
        ogeom_bail!(Dangling, "curve is not in this model");
    };
    // The curve sits where its own location puts it, inside wherever the
    // edge is placed: a prism's far end edges are its profile's curves
    // carried along it by a location, not new curves.
    let own = location.composed(model.datums())?;
    let outer = edge.transform(model.datums())?;
    let own_width = data.tolerance.get();
    let bounds = model.children_of(edge)?;
    let reach_of = |vertex: Option<&Shape>| {
        vertex
            .and_then(|v| model.node(v))
            .and_then(|n| n.data().as_vertex())
            .map_or(own_width, |v| own_width.max(v.tolerance.get()))
    };
    let ends = [reach_of(bounds.first()), reach_of(bounds.last())];
    let mut width = own_width;
    for vertex in &bounds {
        if let Some(v) = model.node(vertex).and_then(|n| n.data().as_vertex()) {
            width = width.max(v.tolerance.get());
        }
    }
    use ogeom_geom::Transformable as _;
    let placement = |p: Point| outer.apply(own.apply(p));
    Ok(Some(Fingerprint {
        start: placement(geometry.point_at(range.0, tol)?),
        middle: placement(geometry.point_at(f64::midpoint(range.0, range.1), tol)?),
        end: placement(geometry.point_at(range.1, tol)?),
        quarters: [
            placement(geometry.point_at(range.0 + 0.25 * (range.1 - range.0), tol)?),
            placement(geometry.point_at(range.0 + 0.75 * (range.1 - range.0), tol)?),
        ],
        curve: geometry
            .clone()
            .transformed(&own, tol)?
            .transformed(&outer, tol)?,
        range: *range,
        width,
        ends,
    }))
}

/// A vertex's position in space.
fn placed(model: &Model, vertex: &Shape) -> OgeomResult<Point> {
    let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
        ogeom_bail!(Construction, "expected a vertex");
    };
    Ok(vertex.transform(model.datums())?.apply(data.point))
}

/// Rebuild a face with merged edges in place of the ones they replaced.
fn rebuild_face(
    model: &mut Model,
    face: &Shape,
    merged: &FastMap<TShapeId, (TShapeId, bool)>,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, bool)>> {
    let mut whole = true;
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()).cloned() else {
        ogeom_bail!(Construction, "expected a face");
    };
    // Nothing to substitute: the face is already built on the shared edges, and
    // rebuilding it would only mint a node identical to the one there.
    let mut touched = false;

    // The wires are read as the face stores them, its own sense left out:
    // the rebuilt face takes that sense back at the end, and reading them
    // under it as well would turn every ring of a reversed face around.
    let stored = face.oriented(Orientation::Forward);
    let mut wires = Vec::new();
    for wire in model.ordered_children_of(&stored)? {
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
        // A ring that merging collapsed onto one edge, walked out and back,
        // bounds no area: the sliver between two coincident strands that the
        // sew has just found to be one edge. Such a ring is dropped, and a
        // face whose outer ring it was with it.
        let one_edge = !ring.is_empty() && ring.iter().all(|e| e.node() == ring[0].node());
        if one_edge && ring.len() >= 2 {
            if wires.is_empty() {
                return Ok(None);
            }
            whole = false;
            continue;
        }
        let wire = match make_wire(model, &ring, tol) {
            Ok(w) => w.shape,
            Err(e) => {
                if std::env::var_os("OGEOM_DEBUG_SEW").is_some() {
                    eprintln!("SEW WIRE FAIL: {e}");
                    for edge in &ring {
                        if let Some((a, b)) = edge_vertices(model, edge)? {
                            let (pa, pb) = (placed(model, &a)?, placed(model, &b)?);
                            let fp = fingerprint(model, edge, tol)?;
                            eprintln!(
                                "   edge {:?}{} ({:.5},{:.5},{:.5}) v{} -> ({:.5},{:.5},{:.5}) v{} width {:.2e}",
                                edge.node(),
                                if edge.orientation() == ogeom_topo::Orientation::Reversed {
                                    " rev"
                                } else {
                                    ""
                                },
                                pa.x,
                                pa.y,
                                pa.z,
                                a.node().index(),
                                pb.x,
                                pb.y,
                                pb.z,
                                b.node().index(),
                                fp.map_or(0.0, |f| f.width)
                            );
                        }
                    }
                }
                return Err(e);
            }
        };
        wires.push(wire);
    }
    if !touched {
        return Ok(Some((face.clone(), whole)));
    }
    if wires.is_empty() {
        return Ok(None);
    }
    let sewn = match make_face_on(model, data.surface, &wires, tol) {
        Ok(built) => built.shape,
        Err(e) => {
            if std::env::var_os("OGEOM_DEBUG_SEW").is_some() {
                eprintln!("SEW FACE FAIL: {e}");
                for (wi, wire) in wires.iter().enumerate() {
                    for edge in model.ordered_children_of(wire)? {
                        if let Some((a, b)) = edge_vertices(model, &edge)? {
                            let (pa, pb) = (placed(model, &a)?, placed(model, &b)?);
                            eprintln!(
                                "   wire {wi} edge {:?}{} ({:.5},{:.5},{:.5}) v{} -> ({:.5},{:.5},{:.5}) v{}",
                                edge.node(),
                                if edge.orientation() == ogeom_topo::Orientation::Reversed {
                                    " rev"
                                } else {
                                    ""
                                },
                                pa.x,
                                pa.y,
                                pa.z,
                                a.node().index(),
                                pb.x,
                                pb.y,
                                pb.z,
                                b.node().index()
                            );
                        }
                    }
                }
            }
            return Err(e);
        }
    };
    Ok(Some((
        if face.orientation() == Orientation::Reversed {
            sewn.reversed()
        } else {
            sewn
        },
        whole,
    )))
}

/// How many times the faces use each edge, and how many of those uses walk
/// it forward.
fn walks<'a>(
    model: &Model,
    faces: impl Iterator<Item = &'a Shape>,
) -> OgeomResult<FastMap<(TShapeId, ogeom_topo::Location), (usize, usize)>> {
    let mut uses: FastMap<(TShapeId, ogeom_topo::Location), (usize, usize)> = FastMap::default();
    for face in faces {
        for edge in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            if model
                .node(&edge)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|d| d.degenerate)
            {
                continue;
            }
            let walk = uses
                .entry((edge.node(), edge.location().clone()))
                .or_default();
            walk.0 += 1;
            if edge.orientation() == Orientation::Forward {
                walk.1 += 1;
            }
        }
    }
    Ok(uses)
}

/// Whether every edge used by two faces is walked once each way.
fn walked_each_way(uses: &FastMap<(TShapeId, ogeom_topo::Location), (usize, usize)>) -> bool {
    uses.values()
        .all(|(count, forward)| *count != 2 || *forward == 1)
}

/// Group faces by whether they share an edge, transitively.
fn connected_groups(model: &Model, faces: &[Shape]) -> OgeomResult<Vec<Vec<usize>>> {
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

    // Union-find: each face joins the first face seen with each of its
    // edges, which joins it to every face sharing that edge transitively.
    let mut first_user: FastMap<TShapeId, usize> = FastMap::default();
    for (i, edges) in edges_of.iter().enumerate() {
        for edge in edges {
            let j = *first_user.entry(*edge).or_insert(i);
            let (a, b) = (find(&mut group_of, j), find(&mut group_of, i));
            if a != b {
                group_of[b] = a;
            }
        }
    }

    let mut groups: FastMap<usize, Vec<usize>> = FastMap::default();
    for i in 0..faces.len() {
        groups.entry(find(&mut group_of, i)).or_default().push(i);
    }
    let mut out: Vec<Vec<usize>> = groups.into_values().collect();
    // Deterministic: a result whose shells come back in a different order each
    // run is one nobody can compare against.
    out.sort_by_key(|group| group.first().map(|&i| faces[i].node()));
    Ok(out)
}

/// Follow a union-find chain to its root.
fn find(parent: &mut [usize], mut i: usize) -> usize {
    // Halving the path on the way up keeps every chain short.
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}

/// Edges still used by exactly one face.
fn free_edges(model: &Model, faces: &[Shape]) -> OgeomResult<Vec<Shape>> {
    let mut uses: FastMap<TShapeId, (usize, Shape)> = FastMap::default();
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
    use ogeom_geom::PlaneSurface;
    use ogeom_math::{Frame, Plane, Vector};
    use ogeom_topo::Location;

    const T: Tolerances = Tolerances::millimetres();

    fn fine() -> ogeom_mesh::Deflection {
        ogeom_mesh::Deflection {
            chord: 0.02,
            ..ogeom_mesh::Deflection::default()
        }
    }

    /// A flat face on `corners`, built from its own fresh edges so it shares
    /// nothing with anything else.
    fn loose_polygon<const N: usize>(model: &mut Model, corners: [Point; N]) -> Shape {
        let wire = make_polygon(model, &corners, true, T).unwrap().shape;
        let normal =
            ogeom_math::Direction::from_cross(corners[1] - corners[0], corners[2] - corners[1], T)
                .unwrap();
        let frame = ogeom_math::Frame::new(
            corners[0],
            normal,
            ogeom_math::Direction::new(corners[1] - corners[0], T).unwrap(),
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
                ogeom_math::Point2::new(l.x, l.y)
            };
            crate::attach_pcurve(
                model,
                &edge,
                ogeom_geom::Line2d::segment(flat(pa), flat(pb), T)
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

    /// Two unit squares side by side, the second a twentieth of a
    /// millimetre off: sewn across a tenth they are one shell of one joined
    /// pair with nothing broken; across a hundredth they stay apart.
    #[test]
    fn squares_a_little_apart_sew_across_the_gap_asked_for() {
        let at = |x: f64, y: f64| Point::new(x, y, 0.0);
        for (gap, shells, joined) in [(0.1, 1, 1), (0.01, 2, 0)] {
            let mut model = Model::new();
            let first = loose_polygon(
                &mut model,
                [at(0.0, 0.0), at(1.0, 0.0), at(1.0, 1.0), at(0.0, 1.0)],
            );
            let second = loose_polygon(
                &mut model,
                [at(1.05, 0.0), at(2.05, 0.0), at(2.05, 1.0), at(1.05, 1.0)],
            );
            let sewn = sew_within(&mut model, &[first, second], gap, T).unwrap();
            assert_eq!(sewn.shells.len(), shells, "gap {gap}");
            assert_eq!(sewn.joined, joined, "gap {gap}");
            if shells == 1 {
                let d = crate::check(&model, &sewn.shells[0], T).unwrap();
                assert!(d.is_usable(), "{d}");
            }
        }
    }

    /// A 10 mm edge against two 5 mm edges along it: the long edge is split
    /// where the short ones meet, and the three faces sew into one shell.
    #[test]
    fn a_long_edge_against_two_short_ones_is_split_and_sewn() {
        let at = |x: f64, y: f64| Point::new(x, y, 0.0);
        let mut model = Model::new();
        let long = loose_polygon(
            &mut model,
            [at(0.0, 0.0), at(10.0, 0.0), at(10.0, 10.0), at(0.0, 10.0)],
        );
        let left = loose_polygon(
            &mut model,
            [at(0.0, -5.0), at(5.0, -5.0), at(5.0, 0.0), at(0.0, 0.0)],
        );
        let right = loose_polygon(
            &mut model,
            [at(5.0, -5.0), at(10.0, -5.0), at(10.0, 0.0), at(5.0, 0.0)],
        );
        let sewn = sew_within(&mut model, &[long, left, right], 1e-3, T).unwrap();
        assert_eq!(sewn.shells.len(), 1);
        assert_eq!(sewn.joined, 3);
        let d = crate::check(&model, &sewn.shells[0], T).unwrap();
        assert!(d.is_usable(), "{d}");
    }

    /// A prism's walls with a lid built apart at each end: the walls' end
    /// edges are the profile's edges, the top ones carried up the prism,
    /// and both lids sew on, closing the box.
    #[test]
    fn a_prism_s_walls_sew_to_lids_at_both_ends() {
        let mut model = Model::new();
        let at = |x: f64, y: f64, z: f64| Point::new(x, y, z);
        let profile = make_polygon(
            &mut model,
            &[
                at(0.0, 0.0, 0.0),
                at(10.0, 0.0, 0.0),
                at(10.0, 10.0, 0.0),
                at(0.0, 10.0, 0.0),
            ],
            true,
            T,
        )
        .unwrap()
        .shape;
        let walls = crate::make_prism(&mut model, &profile, Vector::new(0.0, 0.0, 5.0), T)
            .unwrap()
            .shape;
        let mut faces = explore_unique(&model, &walls, ShapeType::Face).unwrap();
        faces.push(loose_polygon(
            &mut model,
            [
                at(0.0, 0.0, 0.0),
                at(0.0, 10.0, 0.0),
                at(10.0, 10.0, 0.0),
                at(10.0, 0.0, 0.0),
            ],
        ));
        faces.push(loose_polygon(
            &mut model,
            [
                at(0.0, 0.0, 5.0),
                at(10.0, 0.0, 5.0),
                at(10.0, 10.0, 5.0),
                at(0.0, 10.0, 5.0),
            ],
        ));
        let sewn = sew(&mut model, &faces, T).unwrap();
        assert_eq!(sewn.shells.len(), 1);
        assert_eq!(sewn.joined, 8);
        assert!(sewn.free_edges.is_empty());
        assert!(is_shell_closed(&model, &sewn.shells[0]).unwrap());
        let solid = crate::make_solid(&mut model, &sewn.shells).unwrap().shape;
        let volume = crate::volume_properties(&model, &solid, fine(), T)
            .unwrap()
            .mass;
        assert!((volume.abs() - 500.0).abs() < 1e-6, "{volume}");
    }

    /// A 10 x 10 x 5 box's six faces: a prism's four walls, turned to face
    /// in where `walls_in`, and two lids built apart, the bottom one wound
    /// to face up into the box where `bottom_in` and down out of it
    /// otherwise; the top one faces up.
    fn box_faces(model: &mut Model, walls_in: bool, bottom_in: bool) -> Vec<Shape> {
        box_faces_in(model, walls_in, bottom_in, Frame::WORLD)
    }

    /// [`box_faces`] laid out in `frame`.
    fn box_faces_in(
        model: &mut Model,
        walls_in: bool,
        bottom_in: bool,
        frame: Frame,
    ) -> Vec<Shape> {
        let at = |x: f64, y: f64, z: f64| {
            frame.origin()
                + frame.x().vector() * x
                + frame.y().vector() * y
                + frame.z().vector() * z
        };
        let square = |z: f64| {
            [
                at(0.0, 0.0, z),
                at(10.0, 0.0, z),
                at(10.0, 10.0, z),
                at(0.0, 10.0, z),
            ]
        };
        let profile = make_polygon(model, &square(0.0), true, T).unwrap().shape;
        let walls = crate::make_prism(model, &profile, frame.z().vector() * 5.0, T)
            .unwrap()
            .shape;
        let mut faces: Vec<Shape> = explore_unique(model, &walls, ShapeType::Face)
            .unwrap()
            .into_iter()
            .map(|f| if walls_in { f.reversed() } else { f })
            .collect();
        let mut bottom = square(0.0);
        if !bottom_in {
            bottom.reverse();
        }
        faces.push(loose_polygon(model, bottom));
        faces.push(loose_polygon(model, square(5.0)));
        faces
    }

    /// The solid a sewn closed shell bounds, checked valid with no face
    /// facing into it, and its volume.
    fn solid_volume(model: &mut Model, shell: &Shape) -> f64 {
        let solid = crate::make_solid(model, std::slice::from_ref(shell))
            .unwrap()
            .shape;
        let d = crate::check(model, &solid, T).unwrap();
        assert!(d.is_valid(), "{d}");
        assert!(
            crate::inside_out_faces(model, &solid, T)
                .unwrap()
                .is_empty()
        );
        crate::volume_properties(model, &solid, ogeom_mesh::Deflection::default(), T)
            .unwrap()
            .mass
    }

    /// Both lids of a box filled the same way round, so the bottom one
    /// faces into the box: sewn, that lid is turned over, every edge is
    /// walked once each way and the box measures its volume.
    #[test]
    fn a_lid_facing_into_the_box_is_turned_over_when_sewn() {
        let mut model = Model::new();
        let faces = box_faces(&mut model, false, true);
        let sewn = sew(&mut model, &faces, T).unwrap();
        assert!(sewn.free_edges.is_empty());
        assert_eq!(
            crate::check::edges_walked_one_way(&model, &sewn.shells[0]),
            0
        );
        let image = sewn.history.modified(&faces[4]);
        assert_eq!(image.len(), 1);
        assert_ne!(image[0].orientation(), faces[4].orientation());
        let volume = solid_volume(&mut model, &sewn.shells[0]);
        assert!((volume - 500.0).abs() < 1e-9, "{volume}");
    }

    /// A box whose walls face in and whose lids face out: the lids are
    /// turned to agree with the walls, and the shell, then facing in as a
    /// whole, is turned over to face out.
    #[test]
    fn a_closed_shell_whose_faces_disagree_is_turned_to_face_out() {
        let mut model = Model::new();
        let faces = box_faces(&mut model, true, false);
        let sewn = sew(&mut model, &faces, T).unwrap();
        for face in &faces[..4] {
            assert_ne!(
                sewn.history.modified(face)[0].orientation(),
                face.orientation()
            );
        }
        for face in &faces[4..] {
            assert!(sewn.history.copy_of(face).is_some());
        }
        let volume = solid_volume(&mut model, &sewn.shells[0]);
        assert!((volume - 500.0).abs() < 1e-9, "{volume}");
    }

    /// The box whose walls face in and whose lids face out, up to 2e7 from
    /// the origin: which way the sewn shell faces is read off the volume it
    /// encloses, a sum the box's own size, so the walls are turned there as
    /// at the origin and the shell faces out.
    #[test]
    fn a_closed_shell_far_from_the_origin_whose_faces_disagree_is_turned_to_face_out() {
        for o in [0.0, 1e6, 1e7, 2e7] {
            let mut model = Model::new();
            let n = ogeom_math::Direction::new(Vector::new(1.0, 1.0, 1.0), T).unwrap();
            let frame = Frame::new(Point::new(o, -o, o), n, ogeom_math::Direction::X, T).unwrap();
            let faces = box_faces_in(&mut model, true, false, frame);
            let sewn = sew(&mut model, &faces, T).unwrap();
            for face in &faces[..4] {
                assert_ne!(
                    sewn.history.modified(face)[0].orientation(),
                    face.orientation(),
                    "o {o:e}: a wall facing in is turned"
                );
            }
            let mesh = ogeom_mesh::triangulate(
                &model,
                &sewn.shells[0],
                ogeom_mesh::Deflection::default(),
                T,
            )
            .unwrap();
            let volume = mesh.volume();
            assert!((volume - 500.0).abs() < 1e-6, "o {o:e}: {volume}");
        }
    }

    /// A box whose six faces all face in, so every edge is already walked
    /// once each way: sewing turns nothing (the shell may bound a void),
    /// and the solid made of it faces out.
    #[test]
    fn a_closed_shell_facing_in_throughout_makes_a_solid_facing_out() {
        let mut model = Model::new();
        let mut faces = box_faces(&mut model, true, true);
        faces[5] = faces[5].reversed();
        let sewn = sew(&mut model, &faces, T).unwrap();
        assert!(sewn.free_edges.is_empty());
        for face in &faces {
            let image = sewn
                .history
                .copy_of(face)
                .cloned()
                .unwrap_or_else(|| sewn.history.modified(face)[0].clone());
            assert_eq!(image.orientation(), face.orientation());
        }
        let volume = solid_volume(&mut model, &sewn.shells[0]);
        assert!((volume - 500.0).abs() < 1e-9, "{volume}");
    }

    /// Three unit squares in a row, the middle one wound to face down:
    /// sewn, the middle one is turned to agree with the two beside it, the
    /// others are copies, and every shared edge is walked once each way.
    #[test]
    fn a_sheet_turns_the_one_face_walking_against_its_neighbours() {
        let mut model = Model::new();
        let at = |x: f64, y: f64| Point::new(x, y, 0.0);
        let square = |x: f64| [at(x, 0.0), at(x + 1.0, 0.0), at(x + 1.0, 1.0), at(x, 1.0)];
        let mut middle = square(1.0);
        middle.reverse();
        let faces = vec![
            loose_polygon(&mut model, square(0.0)),
            loose_polygon(&mut model, middle),
            loose_polygon(&mut model, square(2.0)),
        ];
        let sewn = sew(&mut model, &faces, T).unwrap();
        assert_eq!(sewn.shells.len(), 1);
        assert_eq!(sewn.joined, 2);
        assert!(sewn.edges_walked_one_way(&model).unwrap().is_empty());
        let image = sewn.history.modified(&faces[1]);
        assert_eq!(image.len(), 1);
        assert_ne!(image[0].orientation(), faces[1].orientation());
        for face in [&faces[0], &faces[2]] {
            assert!(sewn.history.copy_of(face).is_some());
        }
    }

    /// A band of twelve triangles around a circle, given a half twist so
    /// that its two sides are one: no orientation walks every shared edge
    /// once each way. Sewn, it is one open shell whose faces are all kept
    /// as given, and the edges walked the same way are named.
    #[test]
    fn a_half_twisted_band_is_kept_and_named() {
        let mut model = Model::new();
        let (radius, half_width, steps) = (5.0, 1.0, 6_u32);
        let at = |k: u32, side: f64| {
            let u = std::f64::consts::TAU * f64::from(k) / f64::from(steps);
            let (radial, up) = ((u / 2.0).cos(), (u / 2.0).sin());
            Point::new(
                (radius + side * half_width * radial) * u.cos(),
                (radius + side * half_width * radial) * u.sin(),
                side * half_width * up,
            )
        };
        let mut faces = Vec::new();
        for k in 0..steps {
            faces.push(loose_polygon(
                &mut model,
                [at(k, -1.0), at(k + 1, -1.0), at(k, 1.0)],
            ));
            faces.push(loose_polygon(
                &mut model,
                [at(k + 1, -1.0), at(k + 1, 1.0), at(k, 1.0)],
            ));
        }
        let sewn = sew(&mut model, &faces, T).unwrap();
        assert_eq!(sewn.shells.len(), 1);
        // Six rungs and six diagonals shared, the rim of twelve edges free.
        assert_eq!(sewn.joined, 12);
        assert_eq!(sewn.free_edges.len(), 12);
        for face in &faces {
            assert!(sewn.history.copy_of(face).is_some());
        }
        assert!(!sewn.edges_walked_one_way(&model).unwrap().is_empty());
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
        let left = loose_polygon(
            &mut model,
            [
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
                Point::new(1.0, 1.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
            ],
        );
        let right = loose_polygon(
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
        // A sheet, so it still has a boundary: six free edges round the
        // outside, and the shared one is not among them.
        assert_eq!(sewn.free_edges.len(), 6);
        assert!(!is_shell_closed(&model, &sewn.shells[0]).unwrap());
        assert!(sewn.history.is_affected(&left));
    }

    /// A face presented reversed, its rings wound against its surface's
    /// normal, beside one presented as built: both face up, so the edge
    /// they share is walked once each way. Rebuilt onto the shared edge,
    /// the reversed face still walks it against its neighbour.
    #[test]
    fn a_reversed_face_rebuilt_onto_a_shared_edge_keeps_its_walk() {
        let at = |x: f64, y: f64| Point::new(x, y, 0.0);
        let mut model = Model::new();
        let left = loose_polygon(
            &mut model,
            [at(0.0, 0.0), at(1.0, 0.0), at(1.0, 1.0), at(0.0, 1.0)],
        );
        // Wound clockwise seen from above: its plane faces down, and the
        // face is presented the other way round.
        let right = loose_polygon(
            &mut model,
            [at(1.0, 0.0), at(1.0, 1.0), at(2.0, 1.0), at(2.0, 0.0)],
        )
        .reversed();
        let sewn = sew(&mut model, &[left, right], T).unwrap();
        assert_eq!(sewn.joined, 1, "one shared edge");
        let mut walks: FastMap<TShapeId, (usize, usize)> = FastMap::default();
        for face in ogeom_topo::explore(
            &model,
            &sewn.shells[0],
            ogeom_topo::Filter::OfType(ShapeType::Face),
        )
        .unwrap()
        {
            for edge in
                ogeom_topo::explore(&model, &face, ogeom_topo::Filter::OfType(ShapeType::Edge))
                    .unwrap()
            {
                let entry = walks.entry(edge.node()).or_default();
                if edge.orientation() == Orientation::Forward {
                    entry.0 += 1;
                } else {
                    entry.1 += 1;
                }
            }
        }
        let shared: Vec<_> = walks.values().filter(|(f, r)| f + r == 2).collect();
        assert_eq!(shared, [&(1, 1)], "the shared edge, walked once each way");
    }

    /// Twin edges whose ends are two vertices apart by more than either
    /// vertex's own tolerance, though within the edges': the right square's
    /// corner sits a hundredth below the left's, and its shared edge is loose
    /// enough to be the left's. Merged, the edges must end on one vertex, or
    /// the right square's bottom edge still ends at its own corner and its
    /// wire opens.
    #[test]
    fn twin_edges_join_the_vertices_they_end_on() {
        let mut model = Model::new();
        let left = loose_polygon(
            &mut model,
            [
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
                Point::new(1.0, 1.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
            ],
        );
        let right = loose_polygon(
            &mut model,
            [
                Point::new(1.0, -0.01, 0.0),
                Point::new(2.0, 0.0, 0.0),
                Point::new(2.0, 1.0, 0.0),
                Point::new(1.0, 1.0, 0.0),
            ],
        );
        for edge in explore_unique(&model, &right, ShapeType::Edge).unwrap() {
            let (a, b) = crate::edge_vertices(&model, &edge).unwrap().unwrap();
            let (pa, pb) = (placed(&model, &a).unwrap(), placed(&model, &b).unwrap());
            if (pa.x - 1.0).abs() < 1e-9
                && (pb.x - 1.0).abs() < 1e-9
                && let Some(node) = model.node_mut(&edge)
                && let NodeData::Edge(data) = node.data_mut()
            {
                data.tolerance = data.tolerance.widen_to(0.012);
            }
        }
        let sewn = sew(&mut model, &[left, right], T).unwrap();
        assert_eq!(sewn.joined, 1, "the loose edge is the left square's");
        assert_eq!(sewn.shells.len(), 1);
        assert_eq!(
            explore_unique(&model, &sewn.shells[0], ShapeType::Vertex)
                .unwrap()
                .len(),
            6,
            "the two corners at the bottom of the shared edge are one"
        );
        assert_eq!(sewn.free_edges.len(), 6);
    }

    #[test]
    fn faces_that_do_not_meet_stay_in_separate_shells() {
        // Claiming one shell would claim a closure that is not there.
        let mut model = Model::new();
        let here = loose_polygon(
            &mut model,
            [
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
                Point::new(1.0, 1.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
            ],
        );
        let far = loose_polygon(
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
        // under test is that sewing does not *break* a shell that was closed,
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
        let circle = ogeom_math::Circle::new(Frame::WORLD, 1.0, T).unwrap();
        let upper = crate::make_edge(
            &mut model,
            ogeom_geom::CircleCurve::new(circle).into(),
            (0.0, std::f64::consts::PI),
            T,
        )
        .unwrap()
        .shape;
        let lower = crate::make_edge(
            &mut model,
            ogeom_geom::CircleCurve::new(circle).into(),
            (std::f64::consts::PI, std::f64::consts::TAU),
            T,
        )
        .unwrap()
        .shape;

        let a = fingerprint(&model, &upper, T).unwrap().unwrap();
        let b = fingerprint(&model, &lower, T).unwrap().unwrap();
        assert!(
            a.same_as(&b, T).unwrap().is_none(),
            "two different arcs were called the same edge"
        );
        assert!(a.same_as(&a, T).unwrap() == Some(false));
    }

    #[test]
    fn a_wide_vertex_at_one_end_does_not_join_edges_apart_at_the_other() {
        // Two sides of a sliver six microns long, a quarter of a micron
        // wide at its open end, meeting at a vertex widened to most of a
        // micron where a section was welded. The open end's vertices are
        // exact, so the sides are two edges and the sliver stays a face.
        let mut model = Model::new();
        let shared = Point::new(0.0, 0.0, 0.0);
        let side = |model: &mut Model, to: Point| {
            crate::make_edge(
                model,
                ogeom_geom::LineCurve::segment(shared, to, T)
                    .unwrap()
                    .into(),
                (0.0, to.distance(shared)),
                T,
            )
            .unwrap()
            .shape
        };
        let a = side(&mut model, Point::new(6e-3, 0.0, 0.0));
        let b = side(&mut model, Point::new(6e-3, 2.7e-4, 0.0));
        for edge in [&a, &b] {
            let start = model.children_of(edge).unwrap()[0].clone();
            model
                .widen(&start, ogeom_core::Tolerance::new(8.8e-4).unwrap())
                .unwrap();
        }
        let pa = fingerprint(&model, &a, T).unwrap().unwrap();
        let pb = fingerprint(&model, &b, T).unwrap().unwrap();
        assert!(pa.same_as(&pb, T).unwrap().is_none());
    }

    #[test]
    fn an_edge_found_the_other_way_round_is_reversed_rather_than_dropped() {
        let mut model = Model::new();
        let up = Point::new(0.0, 0.0, 1.0);
        let down = Point::new(0.0, 0.0, 0.0);
        let a = crate::make_edge(
            &mut model,
            ogeom_geom::LineCurve::segment(down, up, T).unwrap().into(),
            (0.0, 1.0),
            T,
        )
        .unwrap()
        .shape;
        let b = crate::make_edge(
            &mut model,
            ogeom_geom::LineCurve::segment(up, down, T).unwrap().into(),
            (0.0, 1.0),
            T,
        )
        .unwrap()
        .shape;

        let pa = fingerprint(&model, &a, T).unwrap().unwrap();
        let pb = fingerprint(&model, &b, T).unwrap().unwrap();
        assert_eq!(
            pa.same_as(&pb, T).unwrap(),
            Some(true),
            "the same edge, running the other way"
        );
    }

    #[test]
    fn flipped_merges_carry_their_pcurves_the_right_way_round() {
        // Six faces of a unit cube, each built loose with its own vertices
        // and pre-attached pcurves, wound counter-clockwise around the
        // outward normal as a shell is. Sewing merges all twelve edge pairs,
        // and every merge is *flipped*: the two faces walk their shared
        // edge opposite ways. A carried pcurve copied unchanged then makes
        // the losing face walk edges backwards in parameter space; with
        // several such edges in one ring the boundary zigzags, and faces
        // stop triangulating or triangulate degenerately. Two squares are
        // not enough to see it (a single backwards two-point edge self-heals
        // in ring assembly), which is why this test is a cube.
        //
        // The carry must reverse with the merge: consumers map 3D parameters
        // onto the pcurve range proportionally, so swapping the stored
        // range's ends reverses the traversal exactly.
        let mut model = Model::new();
        let c = Point::new;
        let faces = [
            [c(0., 0., 0.), c(0., 1., 0.), c(1., 1., 0.), c(1., 0., 0.)],
            [c(0., 0., 1.), c(1., 0., 1.), c(1., 1., 1.), c(0., 1., 1.)],
            [c(0., 0., 0.), c(1., 0., 0.), c(1., 0., 1.), c(0., 0., 1.)],
            [c(0., 1., 0.), c(0., 1., 1.), c(1., 1., 1.), c(1., 1., 0.)],
            [c(0., 0., 0.), c(0., 0., 1.), c(0., 1., 1.), c(0., 1., 0.)],
            [c(1., 0., 0.), c(1., 1., 0.), c(1., 1., 1.), c(1., 0., 1.)],
        ];
        let built: Vec<Shape> = faces
            .iter()
            .map(|f| loose_polygon(&mut model, *f))
            .collect();
        let sewn = sew(&mut model, &built, T).unwrap();
        assert_eq!(sewn.joined, 12, "every edge pair merged");
        assert!(sewn.free_edges.is_empty());
        for face in ogeom_topo::explore(
            &model,
            &sewn.shells[0],
            ogeom_topo::Filter::OfType(ShapeType::Face),
        )
        .unwrap()
        {
            let mesh = ogeom_mesh::triangulate(&model, &face, fine(), T).unwrap();
            assert_eq!(
                mesh.triangles.len(),
                2,
                "a unit square triangulates into two triangles, whichever twin \
                 survived and whichever way it runs"
            );
        }
    }

    /// One loose edge among tight ones: the tight edges' twin search stays
    /// local, and every pair within its wider width is still a candidate.
    #[test]
    fn one_loose_edge_does_not_widen_every_search() {
        let mut model = Model::new();
        let mut catalogue = Vec::new();
        for k in 0..200 {
            let x = f64::from(k);
            let wire = make_polygon(
                &mut model,
                &[Point::new(x, 0.0, 0.0), Point::new(x + 0.5, 0.0, 0.0)],
                false,
                T,
            )
            .unwrap()
            .shape;
            let edge = model.children_of(&wire).unwrap()[0].clone();
            let print = fingerprint(&model, &edge, T).unwrap().unwrap();
            catalogue.push((edge.node(), print));
        }
        let loose = 150;
        catalogue[loose].1.width = 20.0;
        let bins = StartBins::new(&catalogue, T);
        assert!(bins.twin_candidates(&catalogue, 0).len() <= 3);
        // Edge 140 starts 10 from the loose edge, within its width.
        assert!(bins.twin_candidates(&catalogue, 140).contains(&loose));
        let from_loose = bins.twin_candidates(&catalogue, loose);
        assert!((151..=170).all(|j| from_loose.contains(&j)));
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
