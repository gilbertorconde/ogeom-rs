//! One entry point over a shape nobody promised was well-formed.
//!
//! The exchange readers heal what their corpus exhibits, inline and in a
//! fixed order; a shape built by a caller, or read from a file the corpus
//! does not resemble, has no such pass. This is it: diagnose, mend what a
//! mend is known for, diagnose again, and say what changed and what did
//! not. Nothing is moved that the model's own tolerances do not already
//! call the same place, and nothing is dropped that the model says is
//! there: a face with area is a face, whatever its shape.
//!
//! What is mended, in order:
//!
//! - a wire whose edges are not end to end is put in the order that
//!   walks them, where one exists;
//! - an edge shorter than its own vertices' tolerances (the two ends
//!   the same point by the model's own admission) is collapsed, its two
//!   vertices made one;
//! - an edge holding its column on a closed surface as plain pcurves,
//!   walked once each way on two columns a period apart, is made a seam
//!   of the surface ([`join_seam_columns`](ogeom_algo::join_seam_columns));
//! - an edge with no pcurve on a face it bounds is given the trim
//!   projection can honestly fit, as the readers do;
//! - loose faces (a compound of them, or an open shell) are sewn where
//!   they share edges;
//! - a face facing into its solid's material is turned;
//! - tolerances are tightened to what the geometry needs.
//!
//! What is not, and where it lives: small faces and small solids stay
//! (removing a face opens the shell it is in; that is defeaturing, and
//! the boolean crate has it for the features it knows), and a face across
//! a grid of patches is not rebuilt as one.

use ogeom_core::FastMap;

use ogeom_algo::{
    Diagnosis, History, ProbeCache, check_kept, edge_vertices, linear_properties, make_wire,
    order_edges, sew,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_mesh::Deflection;
use ogeom_topo::{Model, Shape, ShapeType, TShapeId, explore_unique};

use crate::{Reshape, reduce_tolerances};

/// What [`fix_shape`] did, and what it found before and after.
#[derive(Debug, Clone)]
pub struct FixReport {
    /// The diagnosis the shape came in with.
    pub before: Diagnosis,
    /// The diagnosis it leaves with. Not necessarily valid: what has no
    /// mend here is still reported.
    pub after: Diagnosis,
    /// Wires whose edges were put end to end.
    pub wires_reordered: usize,
    /// Edges collapsed to a vertex.
    pub edges_collapsed: usize,
    /// Edges given a pcurve on a face they bound.
    pub edges_trimmed: usize,
    /// Edge pairs sewn, and edges still free after, when sewing ran.
    pub sewn: Option<(usize, usize)>,
    /// Tolerances tightened.
    pub tolerances_reduced: usize,
    /// Tolerances widened so that every vertex is at least as loose as the
    /// edges it bounds and every edge as the faces it bounds.
    pub tolerances_widened: usize,
    /// Faces turned the right way out: each faced into its solid's
    /// material.
    pub faces_turned: usize,
}

/// A fixed shape: the result, its history, and the report.
#[derive(Debug, Clone)]
pub struct Fixed {
    /// The shape as mended. The input where nothing was.
    pub shape: Shape,
    /// What became of every input node.
    pub history: History,
    /// What was done and what remains.
    pub report: FixReport,
}

/// The cap on how far a fitted trim may sit from its surface: the
/// readers' own, a millimetre at unit scale.
const TRIM_CAP: f64 = 1e7;

/// Mend `shape` where a mend is known for what is wrong with it.
///
/// See the module documentation for what is and is not done.
///
/// # Errors
///
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) if a handle
/// fails to resolve; [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
/// if a rebuilt container comes out empty.
pub fn fix_shape(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<Fixed> {
    // Every mend below edits in place what it reaches, so what other
    // shapes hold as well is copied first, and the copies recorded.
    let mut history = ogeom_algo::unshare(model, shape)?;
    // The orientation probe is asked three times, before, while mending
    // and after; each asks again only what the steps between changed.
    let kept = ProbeCache::new();
    let before = check_kept(model, shape, tol, &kept)?;
    let mut current = shape.clone();

    // Wires and small edges, in one rebuild.
    let mut reshape = Reshape::new();
    let wires_reordered = reorder_wires(model, &current, &mut reshape, tol)?;
    let edges_collapsed = collapse_small_edges(model, &current, &mut reshape, tol)?;
    if !reshape.is_empty() {
        let built = reshape.apply(model, &current)?;
        history = history.then(&built.history);
        current = built.shape;
    }

    // Seams stored as plain pcurves, joined.
    ogeom_algo::join_seam_columns(model, &current, tol)?;

    // Trims for edges that have none on a face they bound.
    // Every face planned, the fits made side by side, then attached face
    // by face in order.
    let plans = explore_unique(model, &current, ShapeType::Face)?
        .iter()
        .map(|face| crate::fix::plan_trims(model, face, tol))
        .collect::<OgeomResult<Vec<_>>>()?;
    let fits = crate::fix::fit_trims(&plans, tol.confusion() * TRIM_CAP, tol);
    let mut trims = crate::FixedTrims::default();
    for (plan, fits) in plans.into_iter().zip(fits) {
        crate::fix::attach_trims(model, plan, fits, tol, &mut trims)?;
    }
    let edges_trimmed = trims.fitted;

    // Loose faces sewn: a compound of faces, or an open shell.
    let mut sewn = None;
    let kind = model.kind_of(&current)?;
    if matches!(kind, ShapeType::Compound | ShapeType::Shell) {
        let faces = explore_unique(model, &current, ShapeType::Face)?;
        if !faces.is_empty() {
            let result = sew(model, &faces, tol)?;
            sewn = Some((result.joined, result.free_edges.len()));
            // Rebuilt when sewing joined anything, or when the faces
            // already closed among themselves and only the container said
            // otherwise: a compound of faces that is a shell becomes one.
            if result.joined > 0 || (kind == ShapeType::Compound && result.free_edges.is_empty()) {
                let rebuilt = match (kind, result.shells.len()) {
                    (ShapeType::Shell, 1) => result.shells[0].clone(),
                    _ => ogeom_algo::make_compound(model, &result.shells)?.shape,
                };
                let mut step = result.history;
                step.modify(&current, rebuilt.clone());
                history = history.then(&step);
                current = rebuilt;
            }
        }
    }

    // Faces facing into their solid's material, turned: their loops and
    // their neighbours already agree about which way is out, and the flag
    // is the one that says otherwise.
    let mut faces_turned = 0;
    let mut turn = Reshape::new();
    for solid in explore_unique(model, &current, ShapeType::Solid)? {
        for face in ogeom_algo::inside_out_faces_kept(model, &solid, tol, &kept)? {
            turn.replace(&face, turned_face(model, &face)?.reversed());
            faces_turned += 1;
        }
    }
    if !turn.is_empty() {
        let built = turn.apply(model, &current)?;
        history = history.then(&built.history);
        current = built.shape;
    }

    let tolerances_reduced = reduce_tolerances(model, &current, tol)?;
    // Last, because every step above may leave a vertex tighter than an
    // edge it bounds (a reduction tightens edges and faces, never below
    // what they bound, but a shape can arrive broken), and containment is
    // established only by widening what is bounded.
    let tolerances_widened = ogeom_algo::restore_containment(model, &current)?;
    let after = check_kept(model, &current, tol, &kept)?;
    Ok(Fixed {
        shape: current,
        history,
        report: FixReport {
            before,
            after,
            wires_reordered,
            edges_collapsed,
            edges_trimmed,
            sewn,
            tolerances_reduced,
            tolerances_widened,
            faces_turned,
        },
    })
}

/// Stage a rebuilt wire for every wire whose edges do not walk end to end
/// but can be put in an order that does.
fn reorder_wires(
    model: &mut Model,
    shape: &Shape,
    reshape: &mut Reshape,
    tol: Tolerances,
) -> OgeomResult<usize> {
    let mut count = 0;
    for wire in explore_unique(model, shape, ShapeType::Wire)? {
        let edges = model.ordered_children_of(&wire)?;
        if edges.len() < 2 || walks_end_to_end(model, &edges, tol)? {
            continue;
        }
        // A bag that is no path at all (a branch, a gap) is left as it
        // is and reported by the diagnosis; reordering cannot mend it.
        let Ok(ordered) = order_edges(model, &edges, tol) else {
            continue;
        };
        let rebuilt = make_wire(model, &ordered, tol)?.shape;
        reshape.replace(&wire, rebuilt);
        count += 1;
    }
    Ok(count)
}

/// Whether consecutive edges meet, the last back at the first.
fn walks_end_to_end(model: &Model, edges: &[Shape], tol: Tolerances) -> OgeomResult<bool> {
    for i in 0..edges.len() {
        let (Some((_, end)), Some((next, _))) = (
            edge_vertices(model, &edges[i])?,
            edge_vertices(model, &edges[(i + 1) % edges.len()])?,
        ) else {
            return Ok(false);
        };
        if !end.is_same(&next) && !model.same_position(&end, &next, tol)? {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Stage the collapse of every edge shorter than its vertices' tolerances:
/// the edge goes, and its far vertex becomes its near one.
///
/// A run of such edges collapses to one vertex, not a chain of
/// substitutions: the survivor of each merge is found through the merges
/// before it.
fn collapse_small_edges(
    model: &mut Model,
    shape: &Shape,
    reshape: &mut Reshape,
    tol: Tolerances,
) -> OgeomResult<usize> {
    let mut survivor: FastMap<TShapeId, Shape> = FastMap::default();
    fn root(survivor: &FastMap<TShapeId, Shape>, v: &Shape) -> Shape {
        let mut current = v.clone();
        while let Some(next) = survivor.get(&current.node()) {
            if next.node() == current.node() {
                break;
            }
            current = next.clone();
        }
        current
    }
    let mut count = 0;
    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        let Some((a, b)) = edge_vertices(model, &edge)? else {
            continue;
        };
        if a.is_same(&b) {
            // A closed edge is a loop, however short.
            continue;
        }
        let reach = [&a, &b]
            .iter()
            .filter_map(|v| model.tolerance_of(v).ok().flatten())
            .map(|t| t.get())
            .fold(tol.confusion(), f64::max);
        // An edge is at least as long as the gap between its ends: where
        // that alone passes the reach, nothing needs measuring.
        let ends = |v: &Shape| -> OgeomResult<Option<ogeom_math::Point>> {
            let Some(data) = model.node(v).and_then(|n| n.data().as_vertex()) else {
                return Ok(None);
            };
            Ok(Some(v.transform(model.datums())?.apply(data.point)))
        };
        if let (Some(pa), Some(pb)) = (ends(&a)?, ends(&b)?)
            && pa.distance(pb) > reach
        {
            continue;
        }
        let length = linear_properties(model, &edge, Deflection::default(), tol)?.mass;
        if length > reach {
            continue;
        }
        let (keep, drop) = (root(&survivor, &a), root(&survivor, &b));
        if keep.is_same(&drop) {
            // Already one vertex through earlier collapses; the edge is a
            // loop on it and goes.
            reshape.remove(&edge);
            count += 1;
            continue;
        }
        survivor.insert(drop.node(), keep.clone());
        reshape.remove(&edge);
        count += 1;
    }
    // The survivor stands where it stood, and the curves that ended at each
    // vertex it absorbs still end there: it widens to reach every one, as
    // far as the absorbed vertex stood plus that vertex's own tolerance.
    // Merged without it, the neighbours of a collapsed edge stop the
    // collapsed length short of their vertex and the wire gapes.
    let placed = |model: &Model, v: &Shape| -> OgeomResult<Option<(ogeom_math::Point, f64)>> {
        let Some(data) = model.node(v).and_then(|n| n.data().as_vertex()) else {
            return Ok(None);
        };
        let (point, own) = (data.point, data.tolerance.get());
        Ok(Some((v.transform(model.datums())?.apply(point), own)))
    };
    for vertex in explore_unique(model, shape, ShapeType::Vertex)? {
        let to = root(&survivor, &vertex);
        if !to.is_same(&vertex) {
            if let (Some((from, reach)), Some((at, _))) =
                (placed(model, &vertex)?, placed(model, &to)?)
            {
                let need = from.distance(at) + reach;
                model.widen(&to, ogeom_core::Tolerance::new(need.max(tol.confusion()))?)?;
            }
            reshape.replace(&vertex, to);
        }
    }
    if count > 0 && reshape.is_empty() {
        ogeom_bail!(Construction, "a collapse staged nothing");
    }
    Ok(count)
}

/// A new face node holding the same surface and wires as `face`'s, for
/// replacing it turned: a replacement by the node itself is no change. It
/// stands where `face` does: a face placed by a location (an assembly
/// part's, read from a file) keeps it, or the void it bounds would move out
/// of its solid.
fn turned_face(model: &mut Model, face: &Shape) -> OgeomResult<Shape> {
    let stored = if face.orientation() == ogeom_topo::Orientation::Reversed {
        face.reversed()
    } else {
        face.clone()
    };
    let Some(data) = model
        .node(&stored)
        .and_then(|n| n.data().as_face())
        .cloned()
    else {
        ogeom_bail!(Dangling, "face is not in this model");
    };
    let wires = model.children_of(&Shape::of(stored.node()))?;
    Ok(model
        .add_face(data, &wires)?
        .located(face.location().clone()))
}
