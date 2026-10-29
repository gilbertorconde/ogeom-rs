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
//! - where one edge of a wire ends at a vertex and the next begins at
//!   another, no farther apart than a fitted trim may stand off its
//!   surface, the two vertices are made one, wide enough to reach both:
//!   the same slop as a curve ending away from its vertex, and stated the
//!   same way. A gap past that reach is a missing piece of boundary;
//! - such a gap is filled with a straight edge where that edge is provably
//!   the one lacking: the gap is shorter than every edge of each face it
//!   is in and no other vertex of those faces stands as near its ends, the
//!   line lies on each face's surface within the face's tolerance along the
//!   chart segment between the edges either side, and each face, filled,
//!   checks clean. A gap two faces share gets one edge with a pcurve on
//!   each. Any other gap stays reported;
//! - a boundary that crosses a pole of its surface (a sphere's, a cone's
//!   apex) with no edge spanning the pole's row of the chart is given a
//!   degenerate edge there;
//! - an edge holding its column on a closed surface as plain pcurves,
//!   walked once each way on two columns a period apart, is made a seam
//!   of the surface ([`join_seam_columns`](ogeom_algo::join_seam_columns));
//! - an edge with no pcurve on a face it bounds is given the trim
//!   projection can honestly fit, as the readers do;
//! - a wire that runs out along an edge and straight back along the same
//!   curve, within the tolerance either edge states and on one copy of
//!   the chart, with no other face walking either edge, loses the pair (a
//!   tail); a hole made of such pairs alone goes. Where one of the two runs
//!   farther than the other (a notch), the longer is cut where the shorter
//!   turns back, and the stretch out to the tip goes with the shorter edge
//!   (see `tails`);
//! - a hole on a plane or a cylinder every point of which stands within the
//!   face's tolerance of the hole's boundary, which no neighbour walks,
//!   goes: by the face's own tolerance it encloses nothing (see
//!   `small_wires`);
//! - a boundary that crosses itself once in its surface's chart is mended
//!   where the mend is the one it says (see `crossing`): two consecutive
//!   edges crossing short of the vertex they meet at, by no more than a
//!   fitted trim may stand off its surface, are cut back to the crossing
//!   and the vertex moves there, widened to reach where it stood; any
//!   other crossing of one wire splits its face into the figure eight's
//!   two loops. Each face changed must check clean after and cross itself
//!   nowhere; the report counts the faces still crossing themselves;
//! - a face with no wires on a surface whose domain is its whole extent (a
//!   sphere, a torus, a patch) is bounded by that domain's sides, with a
//!   seam edge where the surface closes and a degenerate edge where it
//!   pinches to a pole, or else marked as covering it;
//! - loose faces (a compound of them, or an open shell) are sewn where
//!   they share edges. Two free edges of a shell's faces that are each
//!   other's twin across a gap no wider than a fitted trim may stand off
//!   its surface, the gap shorter than every edge of either face and no
//!   other vertex as near their ends, first state that gap, each edge and
//!   vertex as wide as its own, where sewing a copy of the shape shows the
//!   pair joins (see `free_gaps`);
//! - faces two of which walk an edge they share the same way are turned
//!   until each such edge is walked once each way, a shell's as sewing
//!   turns them, a solid's keeping the orientation most of each run of
//!   them has;
//! - a run of a solid's faces facing into its material is turned, where
//!   most of its faces read that way;
//! - tolerances are tightened to what the geometry needs;
//! - an edge whose pcurves stand off its curve by more than it states, and
//!   a vertex a curve it bounds ends away from by more than either states,
//!   widen to the gap measured, where that gap is within the reach a
//!   fitted trim is allowed: slop a reader would have recorded the same
//!   way. A gap past that reach is no slop and stays reported.
//!
//! The free edges of the shape as mended come back chained into contours
//! ([`free_bounds`](crate::free_bounds())), each of which
//! [`FreeBounds::wires`](crate::FreeBounds::wires) makes a wire.
//!
//! What is not, and where it lives: small faces and small solids stay
//! (removing a face opens the shell it is in; that is defeaturing, and
//! the boolean crate has it for the features it knows), and a face across
//! a grid of patches is not rebuilt as one.

use ogeom_core::FastMap;

use ogeom_algo::{
    Diagnosis, History, ProbeCache, check_kept, edge_pcurve_gap, edge_vertices, linear_properties,
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
    /// Vertices merged into the vertex beside them across a gap in a wire.
    pub gap_vertices_merged: usize,
    /// Straight edges made across gaps in wires too wide for their
    /// vertices to merge.
    pub lacking_edges_made: usize,
    /// Degenerate edges made where a boundary crosses a pole.
    pub pole_edges_made: usize,
    /// Edges given a pcurve on a face they bound.
    pub edges_trimmed: usize,
    /// Edge pairs sewn, and edges still free after, when sewing ran.
    pub sewn: Option<(usize, usize)>,
    /// Pairs of a shell's free edges that are each other's twin across a
    /// gap, widened to state it before sewing.
    pub free_gaps_stated: usize,
    /// Tolerances tightened.
    pub tolerances_reduced: usize,
    /// Tolerances widened so that every vertex is at least as loose as the
    /// edges it bounds and every edge as the faces it bounds.
    pub tolerances_widened: usize,
    /// Faces turned the right way out: each faced into its solid's
    /// material.
    pub faces_turned: usize,
    /// Edges widened to how far their pcurves stand off their curves.
    pub pcurve_gaps_stated: usize,
    /// Vertices widened to how far the curves they bound end from them.
    pub curve_ends_stated: usize,
    /// Pairs of edges stripped from wires that ran out along one and
    /// straight back along the other.
    pub tails_stripped: usize,
    /// Notches cut off: an edge run back along part of the edge before it,
    /// the shorter edge and the stretch of the longer it ran along gone.
    pub notches_cut: usize,
    /// Holes removed that enclosed no area the face's tolerance resolves:
    /// every point of each within that tolerance of its boundary.
    pub small_wires_removed: usize,
    /// Pairs of consecutive edges crossing short of the vertex they meet
    /// at, cut back to the crossing.
    pub crossings_trimmed: usize,
    /// Faces whose boundary crossed itself, split into its two loops.
    pub faces_split_at_crossings: usize,
    /// Faces of the shape as mended whose boundary still crosses itself
    /// in its surface's chart ([`chart_crossings`](ogeom_algo::chart_crossings)).
    pub crossings_left: usize,
    /// Faces with no wires bounded by the sides of their surface's domain.
    pub natural_bounds_made: usize,
    /// Faces with no wires marked as covering their whole surface, where
    /// the sides of its domain could not bound them.
    pub natural_bounds_marked: usize,
    /// The free edges of the shape as mended, chained into contours.
    pub free_bounds: crate::FreeBounds,
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

    // Gaps between consecutive edges, on the wires as reordered: a wire
    // out of order has no consecutive edges to speak of.
    let mut gaps = Reshape::new();
    let gap_vertices_merged = merge_wire_gaps(model, &current, &mut gaps, tol)?;
    if !gaps.is_empty() {
        let built = gaps.apply(model, &current)?;
        history = history.then(&built.history);
        current = built.shape;
    }

    // Lacking edges across the gaps merging leaves, where a straight edge
    // provably fills them.
    let mut lacking = Reshape::new();
    let lacking_edges_made = crate::lacking::stage_lacking_edges(
        model,
        &current,
        &mut lacking,
        tol.confusion() * TRIM_CAP,
        tol,
    )?;
    if !lacking.is_empty() {
        let built = lacking.apply(model, &current)?;
        history = history.then(&built.history);
        current = built.shape;
    }

    // Degenerate edges where a boundary crosses a pole, on the wires as
    // reordered, collapsed and closed: a collapsed edge at a pole leaves
    // the gap.
    let mut poles = Reshape::new();
    let pole_edges_made = crate::pole::stage_pole_edges(model, &current, &mut poles, tol)?;
    if !poles.is_empty() {
        let built = poles.apply(model, &current)?;
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

    // Tails, on the wires as trimmed: a tail is told by its pcurves, and
    // is stripped before sewing could join one of its edges to a
    // neighbour's.
    let mut tails = Reshape::new();
    let tails_stripped = crate::tails::stage_tails(model, &current, &mut tails, tol)?;
    if !tails.is_empty() {
        let built = tails.apply(model, &current)?;
        history = history.then(&built.history);
        current = built.shape;
    }
    // Notches, on the wires stripped of tails: a tail on a notch's tip
    // has gone, and what runs back is the notch's own edge.
    let mut notches = Reshape::new();
    let notches_cut = crate::tails::stage_notches(model, &current, &mut notches, tol)?;
    if !notches.is_empty() {
        let built = notches.apply(model, &current)?;
        history = history.then(&built.history);
        current = built.shape;
    }

    // Holes enclosing nothing the face's tolerance resolves, read on the
    // wires stripped of tails and notches.
    let mut small = Reshape::new();
    let small_wires_removed =
        crate::small_wires::stage_small_wires(model, &current, &mut small, tol)?;
    if !small.is_empty() {
        let built = small.apply(model, &current)?;
        history = history.then(&built.history);
        current = built.shape;
    }

    // Boundaries crossing themselves in their charts, read on the wires as
    // trimmed and stripped of tails.
    let mut crossings = Reshape::new();
    let crossed = crate::crossing::stage_crossings(
        model,
        &current,
        &mut crossings,
        tol.confusion() * TRIM_CAP,
        tol,
    )?;
    if !crossings.is_empty() {
        let built = crossings.apply(model, &current)?;
        history = history.then(&built.history);
        current = built.shape;
    }
    let mut natural = Reshape::new();
    let (natural_bounds_made, natural_bounds_marked) =
        bound_natural_faces(model, &current, &mut natural, tol)?;
    if !natural.is_empty() {
        let built = natural.apply(model, &current)?;
        history = history.then(&built.history);
        current = built.shape;
    }

    // Loose faces sewn: a compound of faces, or an open shell.
    let mut sewn = None;
    let kind = model.kind_of(&current)?;
    let mut free_gaps_stated = 0;
    if matches!(kind, ShapeType::Compound | ShapeType::Shell) {
        // Twins across a gap, stated as it so that sewing joins them.
        free_gaps_stated =
            crate::free_gaps::state_free_gaps(model, &current, tol.confusion() * TRIM_CAP, tol)?;
        let faces = explore_unique(model, &current, ShapeType::Face)?;
        if !faces.is_empty() {
            let result = sew(model, &faces, tol)?;
            sewn = Some((result.joined, result.free_edges.len()));
            // Rebuilt when sewing joined anything, when the faces already
            // closed among themselves and only the container said
            // otherwise (a compound of faces that is a shell becomes one),
            // or when sewing turned faces two of which walked an edge they
            // share the same way, where no solid holds them.
            let turned = explore_unique(model, &current, ShapeType::Solid)?.is_empty()
                && ogeom_algo::agreeing_orientation(model, &faces)?
                    .is_some_and(|agreement| agreement.turn.contains(&true));
            if result.joined > 0
                || (kind == ShapeType::Compound && result.free_edges.is_empty())
                || turned
            {
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

    // The faces of each solid turned, first by the topology alone, then
    // by the material. See `solid_turns`.
    let mut faces_turned = 0;
    for by_material in [false, true] {
        let mut turn = Reshape::new();
        for solid in explore_unique(model, &current, ShapeType::Solid)? {
            for face in solid_turns(model, &solid, by_material, tol, &kept)? {
                turn.replace(&face, turned_face(model, &face)?.reversed());
                faces_turned += 1;
            }
        }
        if !turn.is_empty() {
            let built = turn.apply(model, &current)?;
            history = history.then(&built.history);
            current = built.shape;
        }
    }

    let tolerances_reduced = reduce_tolerances(model, &current, tol)?;
    // After the reduction, which measures more coarsely than the checker
    // and may leave an edge tighter than its pcurves' gap.
    let pcurve_gaps_stated = state_pcurve_gaps(model, &current, tol)?;
    let curve_ends_stated = state_curve_ends(model, &current, tol)?;
    // Last, because every step above may leave a vertex tighter than an
    // edge it bounds (a reduction tightens edges and faces, never below
    // what they bound, but a shape can arrive broken), and containment is
    // established only by widening what is bounded.
    let tolerances_widened = ogeom_algo::restore_containment(model, &current)?;
    let after = check_kept(model, &current, tol, &kept)?;
    let mut crossings_left = 0;
    for face in explore_unique(model, &current, ShapeType::Face)? {
        if !ogeom_algo::chart_crossings(model, &face, tol)?.is_empty() {
            crossings_left += 1;
        }
    }
    let free_bounds = crate::free_bounds::free_bounds(model, &current)?;
    Ok(Fixed {
        shape: current,
        history,
        report: FixReport {
            before,
            after,
            wires_reordered,
            edges_collapsed,
            gap_vertices_merged,
            lacking_edges_made,
            pole_edges_made,
            edges_trimmed,
            sewn,
            free_gaps_stated,
            tolerances_reduced,
            tolerances_widened,
            faces_turned,
            pcurve_gaps_stated,
            curve_ends_stated,
            tails_stripped,
            notches_cut,
            small_wires_removed,
            crossings_trimmed: crossed.trimmed,
            faces_split_at_crossings: crossed.split,
            crossings_left,
            natural_bounds_made,
            natural_bounds_marked,
            free_bounds,
        },
    })
}

/// Bound every face of `shape` with no wires that is not marked as
/// covering its whole surface, where the surface's parameter domain is its
/// whole extent (a sphere, a torus, a Bézier or B-spline patch): what the
/// face is a face of, and what the mesher and every other reader take a
/// face with no wires to be. The face is rebuilt with the sides of that
/// domain as its wire, a periodic pair of sides one seam edge and a side
/// the surface pinches to a point a degenerate edge, as a whole surface is
/// bounded before it is divided; where those sides have no closed-form
/// curve, or the face so bounded does not check clean, it is marked
/// instead. On any other surface the domain is a window standing in for an
/// unbounded sheet (a plane's says where its face was built), and the face
/// stays reported. Returns how many faces are bounded and how many marked.
fn bound_natural_faces(
    model: &mut Model,
    shape: &Shape,
    reshape: &mut Reshape,
    tol: Tolerances,
) -> OgeomResult<(usize, usize)> {
    use ogeom_geom::{Surface as _, SurfaceKind};
    let (mut bounded, mut marked) = (0, 0);
    for face in explore_unique(model, shape, ShapeType::Face)? {
        if !model.children_of(&face)?.is_empty() {
            continue;
        }
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        if data.natural_restriction {
            continue;
        }
        let Some(surface) = model.geometry().surface(data.surface) else {
            continue;
        };
        if !matches!(
            surface.kind(),
            SurfaceKind::Sphere | SurfaceKind::Torus | SurfaceKind::Bezier | SurfaceKind::BSpline
        ) {
            continue;
        }
        let stored = face.oriented(ogeom_topo::Orientation::Forward);
        if let Some(whole) = crate::divide::bounded_natural(model, &stored, tol)?
            && ogeom_algo::check(model, &whole, tol)?.is_valid()
        {
            reshape.replace(&stored, whole);
            bounded += 1;
            continue;
        }
        if let Some(node) = model.node_mut(&face)
            && let ogeom_topo::NodeData::Face(data) = node.data_mut()
        {
            data.natural_restriction = true;
            marked += 1;
        }
    }
    Ok((bounded, marked))
}

/// Widen every edge whose pcurves stand off its curve farther than it
/// states, by no more than a fitted trim may stand off its surface, to the
/// gap, measured as densely as the checker samples and more. Returns how
/// many edges widened.
fn state_pcurve_gaps(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<usize> {
    let cap = tol.confusion() * TRIM_CAP;
    let edges = explore_unique(model, shape, ShapeType::Edge)?;
    let gaps = {
        let model = &*model;
        ogeom_core::parallel::map_ordered(&edges, |_, edge| edge_pcurve_gap(model, edge, tol))
    };
    let mut slack = Vec::new();
    for (edge, gap) in edges.into_iter().zip(gaps) {
        let Some(gap) = gap? else { continue };
        let stated = model.tolerance_of(&edge)?.map_or(0.0, |t| t.get());
        if gap > stated && gap <= cap {
            slack.push(edge);
        }
    }
    ogeom_algo::state_pcurve_gaps_of(model, &slack, tol)?;
    Ok(slack.len())
}

/// Widen every vertex that a curve it bounds ends away from, farther than
/// both the vertex and the edge state and by no more than a fitted trim
/// may stand off its surface, to the widest such gap. Measured as the
/// checker measures it: the curve's point at the end of its range, placed
/// as the edge is, against the vertex placed as it is. Returns how many
/// vertices widened.
fn state_curve_ends(model: &mut Model, shape: &Shape, tol: Tolerances) -> OgeomResult<usize> {
    let cap = tol.confusion() * TRIM_CAP;
    let mut widest: FastMap<TShapeId, usize> = FastMap::default();
    let mut widened: Vec<(Shape, f64)> = Vec::new();
    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            continue;
        };
        let reach = data.tolerance.get().max(tol.confusion());
        let placement = edge.transform(model.datums())?;
        let bounds = model.children_of(&edge)?;
        for (parameter, vertex) in [(range.0, bounds.first()), (range.1, bounds.last())] {
            let Some(vertex) = vertex else { continue };
            let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
                continue;
            };
            let placed = vertex.transform(model.datums())?.apply(data.point);
            let on_curve =
                placement.apply(ogeom_geom::Curve3d::point_at(geometry, parameter, tol)?);
            let gap = on_curve.distance(placed);
            if gap > reach.max(data.tolerance.get()) && gap <= cap {
                let at = *widest.entry(vertex.node()).or_insert_with(|| {
                    widened.push((vertex.clone(), gap));
                    widened.len() - 1
                });
                widened[at].1 = widened[at].1.max(gap);
            }
        }
    }
    for (vertex, gap) in &widened {
        model.widen(vertex, ogeom_core::Tolerance::new(*gap)?)?;
    }
    Ok(widened.len())
}

/// The faces of `solid` to turn: by the topology alone when not
/// `by_material`, else by the material.
///
/// By the topology: where two faces walk an edge they share the same way,
/// the faces of each run joined by such edges turn until every one is
/// walked once each way, the run keeping the orientation most of its faces
/// have ([`agreeing_orientation`](ogeom_algo::agreeing_orientation)).
///
/// By the material: a run whose faces agree faces one way throughout, so
/// it turns as a whole, and only where more than half its faces read
/// facing into the material. A minority reading in is a misreading, as on
/// a boundary that passes through itself, and turning those faces alone
/// would set them against every neighbour. Where the topology admits no
/// agreement (a band with a half twist), each face reading in is turned on
/// its own.
fn solid_turns(
    model: &Model,
    solid: &Shape,
    by_material: bool,
    tol: Tolerances,
    kept: &ProbeCache,
) -> OgeomResult<Vec<Shape>> {
    let faces = explore_unique(model, solid, ShapeType::Face)?;
    let agreement = ogeom_algo::agreeing_orientation(model, &faces)?;
    if !by_material {
        return Ok(agreement.map_or_else(Vec::new, |agreement| {
            faces
                .iter()
                .zip(&agreement.turn)
                .filter(|(_, turn)| **turn)
                .map(|(face, _)| face.clone())
                .collect()
        }));
    }
    let inward = ogeom_algo::inside_out_faces_kept(model, solid, tol, kept)?;
    let Some(agreement) = agreement else {
        return Ok(inward);
    };
    let reads_in = |k: &usize| inward.iter().any(|face| face.is_same(&faces[*k]));
    let mut out = Vec::new();
    for run in &agreement.runs {
        if 2 * run.iter().filter(|k| reads_in(k)).count() > run.len() {
            out.extend(run.iter().map(|&k| faces[k].clone()));
        }
    }
    Ok(out)
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
        // Built as the bag it is: consecutive ends that meet at two
        // vertices rather than one are no order problem, and are made one
        // with the other gaps.
        let rebuilt = model.add_wire(&ordered)?;
        reshape.replace(&wire, rebuilt);
        count += 1;
    }
    Ok(count)
}

/// Whether consecutive edges meet, the last back at the first: at one
/// vertex, or at two the model's tolerances call one place, which the gap
/// step makes one.
fn walks_end_to_end(model: &Model, edges: &[Shape], tol: Tolerances) -> OgeomResult<bool> {
    for i in 0..edges.len() {
        let (Some((_, end)), Some((next, _))) = (
            edge_vertices(model, &edges[i])?,
            edge_vertices(model, &edges[(i + 1) % edges.len()])?,
        ) else {
            return Ok(false);
        };
        if !end.is_same(&next)
            && !model.same_position(&end, &next, tol)?
            && !one_place(model, &end, &next, tol)?
        {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether two vertices stand within the wider of their tolerances.
fn one_place(model: &Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<bool> {
    let mut reach = tol.confusion();
    let mut points = Vec::with_capacity(2);
    for v in [a, b] {
        let Some(data) = model.node(v).and_then(|n| n.data().as_vertex()) else {
            return Ok(false);
        };
        reach = reach.max(data.tolerance.get());
        points.push(v.transform(model.datums())?.apply(data.point));
    }
    Ok(points[0].distance(points[1]) <= reach)
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

/// Stage the merge of the two vertices at every joint of a wire where one
/// edge ends at a vertex and the next begins at another no farther away
/// than a fitted trim may stand off its surface. Returns how many vertices
/// are merged into another.
///
/// Merges chain: vertices joined through joints are one class, kept as the
/// first met. A class is not merged where merging would change more than
/// the joints it closes: where two of its vertices bound one edge (the
/// edge would close on itself), where it would make two joints of one
/// face meet that share no vertex now (the boundary would pinch), where a
/// vertex of it is held by an edge under a placement of its own (its point
/// is not in the edge's frame), or where a vertex of it lies past the reach
/// from the one kept. What is not merged stays reported.
///
/// The vertex kept widens to reach every vertex it absorbs, as far as that
/// vertex stood plus that vertex's own tolerance, so the curves that ended
/// within each still end within it.
fn merge_wire_gaps(
    model: &mut Model,
    shape: &Shape,
    reshape: &mut Reshape,
    tol: Tolerances,
) -> OgeomResult<usize> {
    fn find(parent: &FastMap<TShapeId, TShapeId>, id: TShapeId) -> TShapeId {
        let mut current = id;
        while let Some(&up) = parent.get(&current) {
            if up == current {
                break;
            }
            current = up;
        }
        current
    }
    let cap = tol.confusion() * TRIM_CAP;
    // The joints of each face's wires together, and of each wire no face
    // holds on its own: the vertex the edge before ends at and the one the
    // edge after begins at.
    let mut boundaries: Vec<Vec<Shape>> = Vec::new();
    let mut held = ogeom_core::FastSet::default();
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let wires = model.children_of(&face)?;
        held.extend(wires.iter().map(Shape::node));
        boundaries.push(wires);
    }
    for wire in explore_unique(model, shape, ShapeType::Wire)? {
        if !held.contains(&wire.node()) {
            boundaries.push(vec![wire]);
        }
    }
    let point = |model: &Model, v: &Shape| -> OgeomResult<Option<ogeom_math::Point>> {
        let Some(data) = model.node(v).and_then(|n| n.data().as_vertex()) else {
            return Ok(None);
        };
        Ok(Some(v.transform(model.datums())?.apply(data.point)))
    };
    let mut joints: Vec<Vec<[TShapeId; 2]>> = Vec::with_capacity(boundaries.len());
    let mut parent: FastMap<TShapeId, TShapeId> = FastMap::default();
    let mut order: Vec<TShapeId> = Vec::new();
    for wires in &boundaries {
        let mut at = Vec::new();
        for wire in wires {
            let edges = model.ordered_children_of(wire)?;
            for i in 0..edges.len() {
                let (Some((_, end)), Some((next, _))) = (
                    edge_vertices(model, &edges[i])?,
                    edge_vertices(model, &edges[(i + 1) % edges.len()])?,
                ) else {
                    continue;
                };
                at.push([end.node(), next.node()]);
                if end.node() == next.node()
                    || !end
                        .location()
                        .is_same_placement(next.location(), model.datums(), tol)?
                {
                    continue;
                }
                let (Some(p), Some(q)) = (point(model, &end)?, point(model, &next)?) else {
                    continue;
                };
                if p.distance(q) > cap {
                    continue;
                }
                for id in [end.node(), next.node()] {
                    if !parent.contains_key(&id) {
                        parent.insert(id, id);
                        order.push(id);
                    }
                }
                let (keep, absorb) = (find(&parent, end.node()), find(&parent, next.node()));
                if keep != absorb {
                    parent.insert(absorb, keep);
                }
            }
        }
        joints.push(at);
    }
    if order.is_empty() {
        return Ok(0);
    }
    let mut refused = ogeom_core::FastSet::default();
    for edge in explore_unique(model, shape, ShapeType::Edge)? {
        let bounds = model.children_of(&Shape::of(edge.node()))?;
        for vertex in &bounds {
            if !vertex.location().is_identity() && parent.contains_key(&vertex.node()) {
                refused.insert(find(&parent, vertex.node()));
            }
        }
        if let [first, .., last] = &bounds[..]
            && first.node() != last.node()
            && parent.contains_key(&first.node())
            && find(&parent, first.node()) == find(&parent, last.node())
        {
            refused.insert(find(&parent, first.node()));
        }
    }
    let class = |id: TShapeId| parent.contains_key(&id).then(|| find(&parent, id));
    for at in &joints {
        let touched: Vec<&[TShapeId; 2]> = at
            .iter()
            .filter(|joint| joint.iter().any(|id| class(*id).is_some()))
            .collect();
        for (j, a) in touched.iter().enumerate() {
            for b in &touched[j + 1..] {
                if a.iter().any(|id| b.contains(id)) {
                    continue;
                }
                for x in a.iter().filter_map(|id| class(*id)) {
                    if b.iter().any(|id| class(*id) == Some(x)) {
                        refused.insert(x);
                    }
                }
            }
        }
    }
    let mut need: FastMap<TShapeId, f64> = FastMap::default();
    let mut absorbed = Vec::new();
    for &id in &order {
        let keep = find(&parent, id);
        if keep == id {
            continue;
        }
        let (Some(from), Some(to)) = (
            point(model, &Shape::of(id))?,
            point(model, &Shape::of(keep))?,
        ) else {
            refused.insert(keep);
            continue;
        };
        let reach = model.tolerance_of(&Shape::of(id))?.map_or(0.0, |t| t.get());
        if from.distance(to) > cap {
            refused.insert(keep);
        }
        let wide = need.entry(keep).or_insert(tol.confusion());
        *wide = wide.max(from.distance(to) + reach);
        absorbed.push((id, keep));
    }
    let mut count = 0;
    for (id, keep) in absorbed {
        if refused.contains(&keep) {
            continue;
        }
        reshape.replace(&Shape::of(id), Shape::of(keep));
        count += 1;
    }
    for (keep, wide) in need {
        if !refused.contains(&keep) {
            model.widen(&Shape::of(keep), ogeom_core::Tolerance::new(wide)?)?;
        }
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
