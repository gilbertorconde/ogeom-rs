//! A cut by blend wedges that sit flush on the solid, made by trimming the
//! faces they sit on rather than by the general fuse. This module trims
//! straight bands; [`ring`] trims rims, for a cut and for a fill.
//!
//! A straight blend's wedge is bounded by planar seats (its legs along the
//! two faces it rounds, its caps) and by the blend, a strip of a cylinder
//! between two of its rulings and two of its cross-sections. Where every
//! seat lies on a face of the solid, facing the same way, and no other face
//! of the solid comes into the wedge, the cut is known without intersecting
//! anything: each face a seat lies on loses the seat's region, a face a seat
//! covers whole goes, and the blend closes what is left. The curves the
//! general fuse would find, the contacts where the blend leaves each face,
//! are the wedge's own edges.
//!
//! Each condition is checked, and where one fails the cut is not made here:
//!
//! - Every vertex of the wedge is a vertex of the solid, or stands inside
//!   one of its edges, a line or an arc, which is split there. Every edge
//!   two seats share is an edge of the solid, and every other edge of a seat
//!   (one it shares with the blend) is an edge of the solid or crosses the
//!   face the seat lies on from one of its vertices to another.
//! - Each seat lies on one planar face of the solid with the same outward
//!   normal, bounded by a run of that face's boundary walked the same way and
//!   by the edges it shares with the blend, which meet no other edge of the
//!   face and leave none of its holes inside the seat. The seat is then the
//!   part of the face on its side of those edges.
//! - The wedge lies in the convex region bounded by its seats' planes and by
//!   the chord plane through the blend's two rulings, and the blend touches
//!   those planes only along its own edges. Every face of the solid but the
//!   ones the seats lie on stands on the far side of one of those planes,
//!   within the tolerance it states, so none reaches into the wedge.
//! - Each new edge leaves its ends clear of every other edge there, as the
//!   general fuse reads a meeting ([`ends_clear`]).
//! - Lumps of a compound tool stand apart, box from box.
//! - The faces the trim makes, with the solid's other faces, walk every edge
//!   they touch twice, once each way.
//!
//! The wedge's inside then lies inside the solid: no face of the solid
//! crosses it, and beside each seat the solid's material is on the wedge's
//! side. So the solid less the wedge is bounded by the solid's faces with
//! the seats taken out, and by the blend turned to face out of the material.
//! The result is put together as the general fuse puts its own together,
//! the faces it would set aside kept and the rest made again, so the two
//! record the same history.

use ogeom_algo::{Built, History, make_edge_between, make_face_with_pcurves, make_wire};
use ogeom_core::{FastMap, FastSet, OgeomError, OgeomResult, Tolerance, Tolerances};
use ogeom_geom::{Curve, Curve3d as _, SurfaceGeometry};
use ogeom_math::{Aabb, Point, Vector};
use ogeom_topo::{
    EdgeRepr, Location, Model, NodeData, Orientation, Shape, ShapeType, SurfaceId, TShapeId,
};
use smallvec::SmallVec;

mod ring;

/// Samples taken along an edge of the solid to say whether it lies along an
/// edge of the wedge.
const ALONG_SAMPLES: u32 = 16;

/// One edge of a face's boundary as the face walks it: the edge, oriented
/// as walked, and the vertices it runs from and to.
#[derive(Clone)]
struct Walked {
    edge: Shape,
    from: TShapeId,
    to: TShapeId,
}

/// What an edge of the wedge stands for in the result: an edge of the
/// solid it lies along, or a new edge across the face its seat lies on.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
enum Target {
    Held(TShapeId),
    New(usize),
}

/// A step of a loop the trim leaves: what the edge stands for, the
/// vertices of the solid it runs from and to, and for an edge of the solid
/// walked as one of its faces walks it, the sense that face walks it in
/// (which alone tells the way round a closed edge).
#[derive(Clone, Copy)]
struct Step {
    target: Target,
    from: TShapeId,
    to: TShapeId,
    sense: Option<Orientation>,
}

/// Per face of the solid the trim changes, the loops it is left with, or
/// none where a seat covers it whole.
type Trimmed = FastMap<usize, Option<Vec<Vec<Step>>>>;

/// An edge the trim adds: a copy of a wedge edge's curve between the
/// solid's vertices its ends stand on.
struct NewEdge {
    curve: Curve,
    range: (f64, f64),
    from: Shape,
    to: Shape,
    tolerance: Tolerance,
}

/// A seat and where it lies: the face of the solid it lies on, and either
/// the whole of that face or a run of one of its loops, with the edges the
/// seat shares with the blend in the order the seat walks them.
struct Seat {
    seat: Shape,
    host: usize,
    run: Option<(usize, usize, usize)>,
    free: Vec<Walked>,
}

/// The trim of one wedge: its seats, its blend, and what each of its edges
/// and vertices stands for.
struct Plan {
    seats: Vec<Seat>,
    blend: Shape,
    blend_loop: Vec<Walked>,
    targets: FastMap<TShapeId, Target>,
    vertices: FastMap<TShapeId, Shape>,
}

/// The solid as the trim reads it: its faces in shell order, each face's
/// loops walked, and which faces walk each edge.
struct Solid {
    faces: Vec<Shape>,
    loops: Vec<Vec<Vec<Walked>>>,
    walkers: FastMap<TShapeId, SmallVec<[usize; 2]>>,
    by_vertex: FastMap<TShapeId, SmallVec<[Shape; 4]>>,
    vertices: Vec<(Shape, Point, f64)>,
    vertex_shapes: FastMap<TShapeId, Shape>,
    edges: FastMap<TShapeId, Shape>,
    split: FastSet<usize>,
    boxes: Vec<Option<(Aabb, f64)>>,
}

impl Solid {
    /// A face's conservative box and the widest tolerance it or any of its
    /// edges and vertices states, read once.
    fn face_box(&mut self, model: &Model, i: usize, tol: Tolerances) -> OgeomResult<(Aabb, f64)> {
        if let Some(known) = self.boxes[i] {
            return Ok(known);
        }
        let face = &self.faces[i];
        let bounds = ogeom_algo::shape_bounds(model, face, tol)?;
        let mut stated = model.tolerance_of(face)?.map_or(0.0, |t| t.get());
        for kind in [ShapeType::Edge, ShapeType::Vertex] {
            for sub in ogeom_topo::explore_unique(model, face, kind)? {
                stated = stated.max(model.tolerance_of(&sub)?.map_or(0.0, |t| t.get()));
            }
        }
        self.boxes[i] = Some((bounds, stated));
        Ok((bounds, stated))
    }
}

/// `a` with the wedges of `b` cut from it, or fused on to it where `fill`
/// is set, by trimming the faces they sit on, with the faces' history as
/// the general fuse records it; `None` where a condition of the module's
/// fails, and the boolean is the general fuse's. A fill is made for rims
/// alone.
pub(crate) fn seated_cut(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    fill: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    #[cfg(feature = "test-hooks")]
    if crate::work::is_unseated() {
        return Ok(None);
    }
    let ringed = ring::is_ringed(model, b)?;
    if fill && !ringed {
        return Ok(None);
    }
    // What reading the wedges and making the result widen is put back
    // where the cut is not made here.
    let mark = model.note_widened();
    let cut = if ringed {
        ring::ring_cut(model, a, b, fill, b, tol)
    } else {
        trimmed_cut(model, a, b, tol)
    };
    match cut {
        Ok(Some(built)) => {
            model.keep_widened(mark);
            #[cfg(feature = "test-hooks")]
            crate::work::add(|w| {
                if fill {
                    w.seated_fills += 1;
                } else {
                    w.seated_cuts += 1;
                }
            });
            Ok(Some(built))
        }
        Err(e @ (OgeomError::Cancelled | OgeomError::Dangling(_))) => {
            model.keep_widened(mark);
            Err(e)
        }
        Ok(None) | Err(_) => {
            model.undo_widened(mark);
            Ok(None)
        }
    }
}

/// `a` less the lumps of `b`, or with them fused on where `fill` is set,
/// where `b` mixes rim wedges with lumps of another kind: the rims trimmed
/// first, as [`seated_cut`] trims them, and the other lumps applied to what
/// is left, by [`seated_cut`] where it makes the boolean and by `rest`
/// otherwise. Made only where every rim's box stands further from every
/// other lump's than the boolean welds past the two lumps' tolerances, so
/// that applying all of `b` at once relates nothing of a rim's to the other
/// lumps; `None` otherwise, and where the trim is not made or `rest` fails.
pub(crate) fn rings_first(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    fill: bool,
    tol: Tolerances,
    rest: impl FnOnce(&mut Model, &Shape, &Shape) -> OgeomResult<Built>,
) -> OgeomResult<Option<Built>> {
    #[cfg(feature = "test-hooks")]
    if crate::work::is_unseated() {
        return Ok(None);
    }
    let Some(all) = lumps(model, b)? else {
        return Ok(None);
    };
    let (mut rings, mut others) = (Vec::new(), Vec::new());
    for lump in all {
        let kind = if ring::is_ringed(model, &lump)? {
            &mut rings
        } else {
            &mut others
        };
        let bounds = ogeom_algo::shape_bounds(model, &lump, tol)?;
        let mut stated = 0.0_f64;
        for sub_kind in [ShapeType::Face, ShapeType::Edge, ShapeType::Vertex] {
            for sub in ogeom_topo::explore_unique(model, &lump, sub_kind)? {
                stated = stated.max(model.tolerance_of(&sub)?.map_or(0.0, |t| t.get()));
            }
        }
        kind.push((lump, bounds, stated));
    }
    if rings.is_empty() || others.is_empty() {
        return Ok(None);
    }
    let weld = crate::sub_shapes::reach(tol);
    for (_, ring_box, ring_stated) in &rings {
        for (_, other_box, other_stated) in &others {
            if ring_box.distance_to_box(other_box) <= weld + ring_stated + other_stated {
                return Ok(None);
            }
        }
    }
    let tool = |model: &mut Model, lumps: &[(Shape, Aabb, f64)]| -> OgeomResult<Shape> {
        match lumps {
            [(one, _, _)] => Ok(one.clone()),
            _ => {
                let shapes: Vec<Shape> = lumps.iter().map(|(lump, _, _)| lump.clone()).collect();
                model.add_compound(&shapes)
            }
        }
    };
    let ringed = tool(model, &rings)?;
    let other = tool(model, &others)?;
    // What the trim widens is put back where the rest is not applied.
    let mark = model.note_widened();
    let trimmed = match ring::ring_cut(model, a, &ringed, fill, b, tol) {
        Ok(Some(trimmed)) => trimmed,
        Err(e @ (OgeomError::Cancelled | OgeomError::Dangling(_))) => {
            model.keep_widened(mark);
            return Err(e);
        }
        Ok(None) | Err(_) => {
            model.undo_widened(mark);
            return Ok(None);
        }
    };
    #[cfg(feature = "test-hooks")]
    crate::work::add(|w| {
        if fill {
            w.seated_fills += 1;
        } else {
            w.seated_cuts += 1;
        }
    });
    let cut = match seated_cut(model, &trimmed.shape, &other, fill, tol) {
        Ok(Some(cut)) => Ok(cut),
        Ok(None) => rest(model, &trimmed.shape, &other),
        Err(e) => Err(e),
    };
    match cut {
        Ok(cut) => {
            model.keep_widened(mark);
            let mut history = trimmed.history.then(&cut.history);
            history.modify(a, cut.shape.clone());
            history.modify(b, cut.shape.clone());
            Ok(Some(Built::new(cut.shape, history)))
        }
        Err(e @ (OgeomError::Cancelled | OgeomError::Dangling(_))) => {
            model.keep_widened(mark);
            Err(e)
        }
        Err(_) => {
            model.undo_widened(mark);
            Ok(None)
        }
    }
}

/// [`seated_cut`], with what it widens left to its caller to keep or put
/// back.
fn trimmed_cut(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some(mut solid) = read_solid(model, a)? else {
        return Ok(None);
    };
    let Some(wedges) = lumps_apart(model, b, tol)? else {
        return Ok(None);
    };
    // Every wedge's vertices found on the solid first: a vertex made on an
    // edge of it changes the loops of the faces walking that edge.
    let mut reads = Vec::with_capacity(wedges.len());
    for wedge in &wedges {
        let Some(read) = read_wedge(model, &mut solid, wedge, tol)? else {
            return Ok(None);
        };
        reads.push(read);
    }
    let mut new_edges: Vec<NewEdge> = Vec::new();
    let mut plans = Vec::with_capacity(wedges.len());
    for read in reads {
        let Some(plan) = plan_wedge(model, &mut solid, read, &mut new_edges, tol)? else {
            return Ok(None);
        };
        plans.push(plan);
    }
    // Per face of the solid, the seats lying on it.
    let mut on_face: FastMap<usize, Vec<(usize, usize)>> = FastMap::default();
    for (w, plan) in plans.iter().enumerate() {
        for (s, seat) in plan.seats.iter().enumerate() {
            on_face.entry(seat.host).or_default().push((w, s));
        }
    }
    // The loops each face the seats lie on is left with, or none where a
    // seat covers it whole.
    let mut trimmed: Trimmed = FastMap::default();
    for (&host, seats) in &on_face {
        let whole = seats
            .iter()
            .filter(|&&(w, s)| plans[w].seats[s].run.is_none())
            .count();
        if whole > 0 {
            if seats.len() > 1 {
                return Ok(None);
            }
            trimmed.insert(host, None);
            continue;
        }
        let mut loops: Vec<Vec<Step>> = Vec::new();
        for (l, ring) in solid.loops[host].iter().enumerate() {
            // Each position of the ring the seats' runs take, and the run
            // starting there.
            let mut taken = vec![false; ring.len()];
            let mut starts: FastMap<usize, (usize, usize)> = FastMap::default();
            for &(w, s) in seats {
                let Some((at_loop, start, len)) = plans[w].seats[s].run else {
                    continue;
                };
                if at_loop != l {
                    continue;
                }
                for k in 0..len {
                    let p = (start + k) % ring.len();
                    if taken[p] {
                        return Ok(None);
                    }
                    taken[p] = true;
                }
                starts.insert(start, (w, s));
            }
            // Walked from a position no run takes, each run replaced by its
            // seat's free edges walked back.
            let Some(first) = (0..ring.len()).find(|&p| !taken[p]) else {
                return Ok(None);
            };
            let mut out = Vec::with_capacity(ring.len());
            for k in 0..ring.len() {
                let p = (first + k) % ring.len();
                if let Some(&(w, s)) = starts.get(&p) {
                    let plan = &plans[w];
                    for free in plan.seats[s].free.iter().rev() {
                        out.push(Step {
                            target: plan.targets[&free.edge.node()],
                            from: plan.vertices[&free.to].node(),
                            to: plan.vertices[&free.from].node(),
                            sense: None,
                        });
                    }
                }
                if !taken[p] {
                    out.push(held_step(&ring[p]));
                }
            }
            loops.push(out);
        }
        trimmed.insert(host, Some(loops));
    }
    // A face no seat lies on whose edge a wedge's vertex split walks the
    // pieces.
    for &f in &solid.split {
        if !trimmed.contains_key(&f) {
            let loops = solid.loops[f]
                .iter()
                .map(|ring| ring.iter().map(held_step).collect())
                .collect();
            trimmed.insert(f, Some(loops));
        }
    }
    // Each blend, turned to face out of the material, on what its edges
    // stand for.
    let blends: Vec<Vec<Step>> = plans
        .iter()
        .map(|plan| {
            plan.blend_loop
                .iter()
                .rev()
                .map(|walked| Step {
                    target: plan.targets[&walked.edge.node()],
                    from: plan.vertices[&walked.to].node(),
                    to: plan.vertices[&walked.from].node(),
                    sense: None,
                })
                .collect()
        })
        .collect();
    if !ends_clear(model, &solid, &new_edges, tol)? {
        return Ok(None);
    }
    if !closes(&solid, &trimmed, &blends) {
        return Ok(None);
    }
    // The faces the general fuse would set aside, which it passes through
    // as they stand: none is one a seat lies on.
    let mut aside = vec![false; solid.faces.len()];
    if let Ok(Some(sides)) = crate::aside::set_aside(model, a, b, tol) {
        for (face, _) in &sides[0].faces {
            if let Some(i) = solid.faces.iter().position(|f| f.is_same(face)) {
                aside[i] = true;
            }
        }
    }
    if trimmed.keys().any(|&host| aside[host]) {
        return Ok(None);
    }
    build(
        model, &solid, &plans, &trimmed, &blends, &new_edges, &aside, a, b, tol,
    )
}

/// Whether each new edge leaves its ends clear of every other edge there,
/// as the general fuse reads a meeting. A straight edge and a curved one
/// leaving a vertex the same way part as the square of the arc they run,
/// at the curve's curvature `k`, so they stand within the boolean's weld
/// reach of each other, and their tolerances (`w` in all), up to the turn
/// `sqrt(2 w k)`. Where that turn is past the one the boolean still reads
/// as tangential ([`TANGENT_TURN`](crate::TANGENT_TURN)), the boolean
/// splits the two where they part, and the cut is left to it; so it is
/// where two straight or two curved edges leave a vertex the same way,
/// which may run on together.
fn ends_clear(
    model: &Model,
    solid: &Solid,
    new_edges: &[NewEdge],
    tol: Tolerances,
) -> OgeomResult<bool> {
    use crate::TANGENT_TURN;
    type Leaves = SmallVec<[(Vector, f64); 2]>;
    // A curve leaving `at` over `range`: at each end of it standing there
    // (both, for a closed curve), the unit direction into it and its
    // curvature; `None` where it has no direction there.
    let leaving = |curve: &Curve, range: (f64, f64), at: Point| -> OgeomResult<Option<Leaves>> {
        let ends = [
            (range.0, 1.0, curve.point_at(range.0, tol)?.distance(at)),
            (range.1, -1.0, curve.point_at(range.1, tol)?.distance(at)),
        ];
        let nearest = ends[0].2.min(ends[1].2);
        let mut out = SmallVec::new();
        for (t, sign, off) in ends {
            if off > nearest + crate::sub_shapes::reach(tol) {
                continue;
            }
            let d = curve.d1_at(t, tol)? * sign;
            let m = d.magnitude();
            if m <= tol.confusion() {
                return Ok(None);
            }
            out.push((d / m, curve.curvature_at(t, tol)?.abs()));
        }
        Ok(Some(out))
    };
    for (i, new) in new_edges.iter().enumerate() {
        for end in [&new.from, &new.to] {
            let Some((at, _)) = vertex_at(model, end) else {
                return Ok(false);
            };
            let Some(mine) = leaving(&new.curve, new.range, at)? else {
                return Ok(false);
            };
            let [(own, own_bend)] = mine.as_slice() else {
                return Ok(false);
            };
            let (own, own_bend) = (*own, *own_bend);
            let mut others: Vec<(Curve, (f64, f64), f64)> = Vec::new();
            for held in solid
                .by_vertex
                .get(&end.node())
                .map_or(&[][..], |v| v.as_slice())
            {
                let degenerate = model
                    .node(held)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| d.degenerate);
                if degenerate {
                    continue;
                }
                let Some(read) = curve_of(model, held) else {
                    return Ok(false);
                };
                others.push(read);
            }
            for (j, other) in new_edges.iter().enumerate() {
                if j != i && (other.from.node() == end.node() || other.to.node() == end.node()) {
                    others.push((other.curve.clone(), other.range, other.tolerance.get()));
                }
            }
            for (curve, range, tolerance) in &others {
                let Some(leaves) = leaving(curve, *range, at)? else {
                    return Ok(false);
                };
                for (theirs, their_bend) in leaves {
                    if own.dot(theirs) <= 0.0 || own.cross(theirs).magnitude() >= TANGENT_TURN {
                        continue;
                    }
                    let straight = [&new.curve, curve].map(|c| matches!(c, Curve::Line(_)));
                    if straight[0] == straight[1] {
                        return Ok(false);
                    }
                    let weld = new.tolerance.get() + tolerance + crate::sub_shapes::reach(tol);
                    if (2.0 * weld * (own_bend + their_bend)).sqrt() >= TANGENT_TURN {
                        return Ok(false);
                    }
                }
            }
        }
    }
    Ok(true)
}

/// An edge of the solid as a face of it walks it, as a step.
fn held_step(walked: &Walked) -> Step {
    Step {
        target: Target::Held(walked.edge.node()),
        from: walked.from,
        to: walked.to,
        sense: Some(walked.edge.orientation()),
    }
}

/// Whether the faces the trim leaves walk every edge they touch twice, once
/// each way, or not at all: the faces of the solid it does not touch, the
/// faces it trims and the blends. A closed edge is walked each way by its
/// two faces' senses.
fn closes(solid: &Solid, trimmed: &Trimmed, blends: &[Vec<Step>]) -> bool {
    let mut walks: FastMap<Target, SmallVec<[Step; 2]>> = FastMap::default();
    let mut touched: Vec<TShapeId> = Vec::new();
    for (&host, loops) in trimmed {
        for ring in &solid.loops[host] {
            touched.extend(ring.iter().map(|w| w.edge.node()));
        }
        for &step in loops.iter().flatten().flatten() {
            walks.entry(step.target).or_default().push(step);
        }
    }
    for &step in blends.iter().flatten() {
        walks.entry(step.target).or_default().push(step);
    }
    for edge in touched {
        walks.entry(Target::Held(edge)).or_default();
    }
    for (target, uses) in &mut walks {
        if let Target::Held(edge) = target {
            for &face in solid.walkers.get(edge).map_or(&[][..], |v| v.as_slice()) {
                if trimmed.contains_key(&face) {
                    continue;
                }
                for ring in &solid.loops[face] {
                    for walked in ring.iter().filter(|w| w.edge.node() == *edge) {
                        uses.push(held_step(walked));
                    }
                }
            }
        }
        match uses.as_slice() {
            // An edge every face walking it lost is gone.
            [] => {}
            [one, two] if one.from != one.to => {
                if (one.from, one.to) != (two.to, two.from) {
                    return false;
                }
            }
            [one, two] => match (one.sense, two.sense) {
                (Some(a), Some(b)) if a != b && one.from == two.from && two.from == two.to => {}
                _ => return false,
            },
            _ => return false,
        }
    }
    true
}

/// The solid's edges and vertices as the result holds them: those of the
/// faces set aside as they stand, every other one copied once.
struct Copies {
    kept: FastSet<TShapeId>,
    vertices: FastMap<TShapeId, Shape>,
    edges: FastMap<TShapeId, Shape>,
}

impl Copies {
    /// The vertex as the result holds it, forward.
    fn vertex(&mut self, model: &mut Model, v: &Shape) -> OgeomResult<Shape> {
        if self.kept.contains(&v.node()) {
            return Ok(v.oriented(Orientation::Forward));
        }
        if let Some(copy) = self.vertices.get(&v.node()) {
            return Ok(copy.clone());
        }
        let Some(data) = model.node(v).and_then(|n| n.data().as_vertex()).cloned() else {
            ogeom_core::ogeom_bail!(Dangling, "a vertex of the solid is not in this model");
        };
        let copy = model.add_vertex(data);
        self.vertices.insert(v.node(), copy.clone());
        Ok(copy)
    }

    /// The edge as the result holds it, forward: its curve, its images on
    /// the faces' surfaces and its tolerance, between its vertices' copies.
    fn edge(&mut self, model: &mut Model, e: &Shape) -> OgeomResult<Shape> {
        if self.kept.contains(&e.node()) {
            return Ok(e.oriented(Orientation::Forward));
        }
        if let Some(copy) = self.edges.get(&e.node()) {
            return Ok(copy.clone());
        }
        let Some(node) = model.node(e) else {
            ogeom_core::ogeom_bail!(Dangling, "an edge of the solid is not in this model");
        };
        let Some(mut data) = node.data().as_edge().cloned() else {
            ogeom_core::ogeom_bail!(Construction, "an edge node holds no edge data");
        };
        let bounds = node.children().to_vec();
        // A path through a face's cached mesh speaks for that face alone.
        data.representations
            .retain(|r| !matches!(r, EdgeRepr::PolygonOnTriangulation { .. }));
        let mut ends = Vec::with_capacity(bounds.len());
        for bound in &bounds {
            ends.push(self.vertex(model, bound)?.oriented(bound.orientation()));
        }
        let copy = model.add_edge(data, &ends)?;
        self.edges.insert(e.node(), copy.clone());
        Ok(copy)
    }
}

/// Make the faces the plans call for and the solid they bound, with the
/// faces' history; `None` where a face made does not close in its chart.
///
/// The result is the one the general fuse assembles: the faces it sets
/// aside stand as they are, each its own copy; every other face of the
/// solid is made again, a face no seat lies on as its exact copy, on its
/// own edges and vertices copied, those it shares with a face set aside
/// kept; and every face made is sewn into the shell.
#[allow(clippy::too_many_arguments, reason = "one construction, all its data")]
fn build(
    model: &mut Model,
    solid: &Solid,
    plans: &[Plan],
    trimmed: &Trimmed,
    blends: &[Vec<Step>],
    new_edges: &[NewEdge],
    aside: &[bool],
    a: &Shape,
    b: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut copies = Copies {
        kept: FastSet::default(),
        vertices: FastMap::default(),
        edges: FastMap::default(),
    };
    for (i, rings) in solid.loops.iter().enumerate() {
        if !aside[i] {
            continue;
        }
        for walked in rings.iter().flatten() {
            copies.kept.insert(walked.edge.node());
            for end in model.node(&walked.edge).map_or(&[][..], |n| n.children()) {
                copies.kept.insert(end.node());
            }
        }
    }
    let mut made: Vec<Shape> = Vec::with_capacity(new_edges.len());
    for new in new_edges {
        let from = copies.vertex(model, &new.from)?;
        let to = copies.vertex(model, &new.to)?;
        let edge = make_edge_between(model, new.curve.clone(), new.range, &from, &to, tol)?.shape;
        model.widen(&edge, new.tolerance)?;
        made.push(edge);
    }
    let mut held: FastMap<TShapeId, Shape> = FastMap::default();
    for walked in solid.loops.iter().flatten().flatten() {
        held.entry(walked.edge.node())
            .or_insert_with(|| walked.edge.oriented(Orientation::Forward));
    }
    // An edge for a walk from one vertex of the solid to another: the
    // edge as the result holds it, turned so it runs that way.
    let walk = |model: &mut Model, copies: &mut Copies, step: Step| -> OgeomResult<Shape> {
        let edge = match step.target {
            Target::Held(id) => copies.edge(model, &held[&id])?,
            Target::New(i) => made[i].clone(),
        };
        if let Some(sense) = step.sense {
            return Ok(edge.oriented(sense));
        }
        if step.from == step.to {
            ogeom_core::ogeom_bail!(Construction, "a closed edge of the trim has no sense");
        }
        let start = copies.vertex(model, &solid.vertex_shapes[&step.from])?;
        let Some((first, _)) = ogeom_algo::edge_vertices(model, &edge)? else {
            ogeom_core::ogeom_bail!(Construction, "an edge of the trim has no vertices");
        };
        Ok(if first.node() == start.node() {
            edge
        } else {
            edge.reversed()
        })
    };
    let mut history = History::new();
    let mut faces: Vec<Shape> = Vec::with_capacity(solid.faces.len() + plans.len());
    let mut checked: Vec<Shape> = Vec::new();
    for (i, face) in solid.faces.iter().enumerate() {
        if aside[i] {
            history.copy(face, face.clone());
            faces.push(face.clone());
            continue;
        }
        let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
            return Ok(None);
        };
        let (surface_id, mut data) = (data.surface, (**data).clone());
        data.triangulation = None;
        let reversed = face.orientation() == Orientation::Reversed;
        let image = match trimmed.get(&i) {
            // Covered whole by a seat.
            Some(None) => {
                history.delete(face);
                continue;
            }
            // No seat on it: its exact copy, wire by wire as its node
            // holds them.
            None => {
                let mut wires = Vec::new();
                for wire in model.node(face).map_or(&[][..], |n| n.children()).to_vec() {
                    let mut edges = Vec::new();
                    for edge in model.node(&wire).map_or(&[][..], |n| n.children()).to_vec() {
                        edges.push(copies.edge(model, &edge)?.oriented(edge.orientation()));
                    }
                    wires.push(model.add_wire(&edges)?.oriented(wire.orientation()));
                }
                let image = model.add_face(data, &wires)?.oriented(face.orientation());
                history.copy(face, image.clone());
                image
            }
            Some(Some(loops)) => {
                let Some(surface) = model.geometry().surface(surface_id).cloned() else {
                    return Ok(None);
                };
                // The new edges the face walks take their image on its plane.
                for step in loops.iter().flatten() {
                    let Target::New(k) = step.target else {
                        continue;
                    };
                    let imaged = model
                        .node(&made[k])
                        .and_then(|n| n.data().as_edge())
                        .is_some_and(|d| d.pcurve_on(surface_id).is_some());
                    if imaged {
                        continue;
                    }
                    let Some(pcurve) = ogeom_intersect::exact_pcurve_over(
                        &new_edges[k].curve,
                        new_edges[k].range,
                        &surface,
                        tol,
                    ) else {
                        return Ok(None);
                    };
                    ogeom_algo::attach_pcurve(
                        model,
                        &made[k],
                        pcurve,
                        surface_id,
                        Location::identity(),
                        new_edges[k].range,
                    )?;
                }
                // The loops as the face's node holds them: walked back where
                // the face is used reversed.
                let mut wires = Vec::with_capacity(loops.len());
                for ring in loops {
                    let mut edges = Vec::with_capacity(ring.len());
                    for &step in ring {
                        edges.push(walk(model, &mut copies, step)?);
                    }
                    if reversed {
                        edges = edges.iter().rev().map(Shape::reversed).collect();
                    }
                    wires.push(make_wire(model, &edges, tol)?.shape);
                }
                data.natural_restriction = false;
                let image = model.add_face(data, &wires)?.oriented(face.orientation());
                history.modify(face, image.clone());
                checked.push(image.clone());
                image
            }
        };
        faces.push(image);
    }
    for (plan, ring) in plans.iter().zip(blends) {
        let Some(NodeData::Face(data)) = model.node(&plan.blend).map(|n| n.data()) else {
            return Ok(None);
        };
        let Some(surface) = model.geometry().surface(data.surface).cloned() else {
            return Ok(None);
        };
        // The blend faces out of the material, the other way from the
        // wedge's.
        let turned = match plan.blend.orientation() {
            Orientation::Forward => Orientation::Reversed,
            _ => Orientation::Forward,
        };
        let mut edges = Vec::with_capacity(ring.len());
        for &step in ring {
            edges.push(walk(model, &mut copies, step)?);
        }
        if turned == Orientation::Reversed {
            edges = edges.iter().rev().map(Shape::reversed).collect();
        }
        let face = make_face_with_pcurves(model, surface, &[edges], tol)?
            .shape
            .oriented(turned);
        history.modify(&plan.blend, face.clone());
        for seat in &plan.seats {
            history.delete(&seat.seat);
        }
        checked.push(face.clone());
        faces.push(face);
    }
    for face in &checked {
        if !ogeom_algo::chart_breaks(model, face, &|_| false, tol)?.is_empty() {
            return Ok(None);
        }
    }
    let shell = model.add_shell(&faces)?;
    let result = model.add_solid(&[shell])?;
    history.modify(a, result.clone());
    history.modify(b, result.clone());
    Ok(Some(Built::new(result, history)))
}

/// The solid as the trim reads it, where it is one shell at no placement
/// whose faces and edges stand at none either and walk their loops head to
/// tail.
fn read_solid(model: &Model, a: &Shape) -> OgeomResult<Option<Solid>> {
    if model.kind_of(a)? != ShapeType::Solid
        || !a.location().is_identity()
        || a.orientation() != Orientation::Forward
    {
        return Ok(None);
    }
    let shells = model.children_of(a)?;
    let [shell] = shells.as_slice() else {
        return Ok(None);
    };
    if !shell.location().is_identity() || shell.orientation() != Orientation::Forward {
        return Ok(None);
    }
    let faces = model.children_of(shell)?;
    let mut loops = Vec::with_capacity(faces.len());
    let mut walkers: FastMap<TShapeId, SmallVec<[usize; 2]>> = FastMap::default();
    let mut by_vertex: FastMap<TShapeId, SmallVec<[Shape; 4]>> = FastMap::default();
    let mut vertices: FastMap<TShapeId, (Shape, Point, f64)> = FastMap::default();
    for (i, face) in faces.iter().enumerate() {
        let Some(rings) = walk_face(model, face)? else {
            return Ok(None);
        };
        for ring in &rings {
            for walked in ring {
                let slot = walkers.entry(walked.edge.node()).or_default();
                if !slot.contains(&i) {
                    slot.push(i);
                }
                for v in [walked.from, walked.to] {
                    let list = by_vertex.entry(v).or_default();
                    if !list.iter().any(|e| e.node() == walked.edge.node()) {
                        list.push(walked.edge.oriented(Orientation::Forward));
                    }
                }
                if let Some((start, end)) = ogeom_algo::edge_vertices(model, &walked.edge)? {
                    for v in [start, end] {
                        if !vertices.contains_key(&v.node()) {
                            let Some((p, t)) = vertex_at(model, &v) else {
                                return Ok(None);
                            };
                            vertices.insert(v.node(), (v, p, t));
                        }
                    }
                }
            }
        }
        loops.push(rings);
    }
    let boxes = vec![None; faces.len()];
    let vertex_shapes = vertices
        .iter()
        .map(|(id, (v, _, _))| (*id, v.clone()))
        .collect();
    let mut edges: FastMap<TShapeId, Shape> = FastMap::default();
    for walked in loops.iter().flatten().flatten() {
        edges
            .entry(walked.edge.node())
            .or_insert_with(|| walked.edge.oriented(Orientation::Forward));
    }
    Ok(Some(Solid {
        faces,
        loops,
        walkers,
        by_vertex,
        vertices: vertices.into_values().collect(),
        vertex_shapes,
        edges,
        split: FastSet::default(),
        boxes,
    }))
}

/// The lumps of a tool: itself where it is a solid, or the solids of a
/// compound of solids, each of one shell, nothing at a placement.
fn lumps(model: &Model, b: &Shape) -> OgeomResult<Option<Vec<Shape>>> {
    if !b.location().is_identity() || b.orientation() != Orientation::Forward {
        return Ok(None);
    }
    let lumps = match model.kind_of(b)? {
        ShapeType::Solid => vec![b.clone()],
        ShapeType::Compound => model.children_of(b)?,
        _ => return Ok(None),
    };
    for lump in &lumps {
        if model.kind_of(lump)? != ShapeType::Solid
            || !lump.location().is_identity()
            || lump.orientation() != Orientation::Forward
            || model.children_of(lump)?.len() != 1
        {
            return Ok(None);
        }
    }
    Ok((!lumps.is_empty()).then_some(lumps))
}

/// [`lumps`] where they stand apart, box from box.
fn lumps_apart(model: &Model, b: &Shape, tol: Tolerances) -> OgeomResult<Option<Vec<Shape>>> {
    let Some(wedges) = lumps(model, b)? else {
        return Ok(None);
    };
    let mut boxes = Vec::with_capacity(wedges.len());
    for wedge in &wedges {
        boxes.push(ogeom_algo::shape_bounds(model, wedge, tol)?);
    }
    for i in 0..boxes.len() {
        for j in i + 1..boxes.len() {
            if boxes[i].intersects(&boxes[j]) {
                return Ok(None);
            }
        }
    }
    Ok(Some(wedges))
}

/// A face's loops as it walks them, each edge with the vertices it runs
/// between; `None` where the face or an edge stands at a placement, an edge
/// has no vertices, or a loop does not run head to tail.
fn walk_face(model: &Model, face: &Shape) -> OgeomResult<Option<Vec<Vec<Walked>>>> {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        return Ok(None);
    };
    if !face.location().is_identity() || !data.location.is_identity() {
        return Ok(None);
    }
    let mut rings = Vec::new();
    for wire in model.children_of(face)? {
        let mut ring: Vec<Walked> = Vec::new();
        for edge in model.ordered_children_of(&wire)? {
            if !edge.location().is_identity() {
                return Ok(None);
            }
            let Some((from, to)) = ogeom_algo::edge_vertices(model, &edge)? else {
                return Ok(None);
            };
            ring.push(Walked {
                edge,
                from: from.node(),
                to: to.node(),
            });
        }
        let closed = !ring.is_empty()
            && (0..ring.len()).all(|k| ring[k].to == ring[(k + 1) % ring.len()].from);
        if !closed {
            return Ok(None);
        }
        rings.push(ring);
    }
    Ok(Some(rings))
}

/// A vertex's point and tolerance, where it stands at no placement.
fn vertex_at(model: &Model, v: &Shape) -> Option<(Point, f64)> {
    if !v.location().is_identity() {
        return None;
    }
    let data = model.node(v)?.data().as_vertex()?;
    Some((data.point, data.tolerance.get()))
}

/// An edge's curve, the range it uses and its tolerance, where it stands at
/// no placement.
fn curve_of(model: &Model, edge: &Shape) -> Option<(Curve, (f64, f64), f64)> {
    let data = model.node(edge)?.data().as_edge()?;
    let Some(EdgeRepr::Curve3d {
        curve,
        location,
        range,
    }) = data.curve3d()
    else {
        return None;
    };
    if !location.is_identity() || !edge.location().is_identity() {
        return None;
    }
    let curve = model.geometry().curve(*curve)?.clone();
    Some((curve, *range, data.tolerance.get()))
}

/// A planar face's surface and its outward normal and a point of its plane,
/// where it is a plane at no placement.
fn plane_of(model: &Model, face: &Shape) -> Option<(SurfaceId, Point, Vector)> {
    let NodeData::Face(data) = model.node(face)?.data() else {
        return None;
    };
    if !face.location().is_identity() || !data.location.is_identity() {
        return None;
    }
    let SurfaceGeometry::Plane(plane) = model.geometry().surface(data.surface)? else {
        return None;
    };
    let normal = plane.plane().normal().vector();
    let outward = if face.orientation() == Orientation::Reversed {
        -normal
    } else {
        normal
    };
    Some((data.surface, plane.plane().origin(), outward))
}

/// A wedge as the trim reads it: its planar faces, each with its loop, a
/// point of its plane and its outward normal; its blend and the blend's
/// loop; the vertex of the solid each of its vertices stands on, and where
/// each stands with its tolerance; and the widest tolerance it states.
struct Read {
    planar: Vec<(Shape, Vec<Walked>, Point, Vector)>,
    blend: Shape,
    blend_loop: Vec<Walked>,
    vertices: FastMap<TShapeId, Shape>,
    points: FastMap<TShapeId, (Point, f64)>,
    tolerance: f64,
}

/// A wedge read, each of its vertices found among the solid's, or made on
/// an edge of the solid it stands on inside ([`split_at`]); `None` where
/// the wedge is not planar faces and one other, or a vertex stands on no
/// vertex or edge of the solid alone.
fn read_wedge(
    model: &mut Model,
    solid: &mut Solid,
    wedge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<Read>> {
    let shells = model.children_of(wedge)?;
    let [shell] = shells.as_slice() else {
        return Ok(None);
    };
    let faces = model.children_of(shell)?;
    let mut planar: Vec<(Shape, Vec<Walked>, Point, Vector)> = Vec::new();
    let mut blend: Option<(Shape, Vec<Walked>)> = None;
    for face in &faces {
        let Some(mut rings) = walk_face(model, face)? else {
            return Ok(None);
        };
        if rings.len() != 1 {
            return Ok(None);
        }
        let ring = rings.remove(0);
        if let Some((_, origin, outward)) = plane_of(model, face) {
            planar.push((face.clone(), ring, origin, outward));
        } else if blend.is_none() {
            blend = Some((face.clone(), ring));
        } else {
            return Ok(None);
        }
    }
    let Some((blend, blend_loop)) = blend else {
        return Ok(None);
    };
    if planar.is_empty() {
        return Ok(None);
    }
    // Every vertex of the wedge is one of the solid's: the one whose ball
    // meets its own, and no other.
    let mut vertices: FastMap<TShapeId, Shape> = FastMap::default();
    let mut points: FastMap<TShapeId, (Point, f64)> = FastMap::default();
    let mut wedge_tolerance = tol.confusion();
    for ring in planar
        .iter()
        .map(|p| &p.1)
        .chain(std::iter::once(&blend_loop))
    {
        for walked in ring {
            wedge_tolerance =
                wedge_tolerance.max(model.tolerance_of(&walked.edge)?.map_or(0.0, |t| t.get()));
            if vertices.contains_key(&walked.from) {
                continue;
            }
            let Some((start, _)) = ogeom_algo::edge_vertices(model, &walked.edge)? else {
                return Ok(None);
            };
            let Some((p, t)) = vertex_at(model, &start) else {
                return Ok(None);
            };
            wedge_tolerance = wedge_tolerance.max(t);
            let mut found = solid
                .vertices
                .iter()
                .filter(|(_, q, s)| p.distance(*q) <= t + s);
            let held = match (found.next(), found.next()) {
                (Some((held, _, _)), None) => held.clone(),
                (None, _) => match split_at(model, solid, p, t, tol)? {
                    Some(made) => made,
                    None => return Ok(None),
                },
                _ => return Ok(None),
            };
            vertices.insert(walked.from, held);
            points.insert(walked.from, (p, t));
        }
    }
    Ok(Some(Read {
        planar,
        blend,
        blend_loop,
        vertices,
        points,
        tolerance: wedge_tolerance,
    }))
}

/// The plan of one wedge's trim on the solid as every wedge's vertices
/// left it, `None` where a condition of the module's fails for it.
fn plan_wedge(
    model: &Model,
    solid: &mut Solid,
    read: Read,
    new_edges: &mut Vec<NewEdge>,
    tol: Tolerances,
) -> OgeomResult<Option<Plan>> {
    let Read {
        planar,
        blend,
        blend_loop,
        vertices,
        points,
        tolerance: wedge_tolerance,
    } = read;
    let Some(cylinder) = blend_cylinder(model, &blend, &blend_loop, tol)? else {
        return Ok(None);
    };
    // What each edge of the wedge stands for.
    let blend_edges: Vec<TShapeId> = blend_loop.iter().map(|w| w.edge.node()).collect();
    let mut targets: FastMap<TShapeId, Target> = FastMap::default();
    for ring in planar.iter().map(|p| &p.1) {
        for walked in ring {
            let id = walked.edge.node();
            if targets.contains_key(&id) {
                continue;
            }
            let (from, to) = (&vertices[&walked.from], &vertices[&walked.to]);
            if from.node() == to.node() {
                return Ok(None);
            }
            let Some((curve, range, tolerance)) = curve_of(model, &walked.edge) else {
                return Ok(None);
            };
            if !matches!(curve, Curve::Line(_) | Curve::Circle(_)) {
                return Ok(None);
            }
            let mut along = Vec::new();
            for held in solid
                .by_vertex
                .get(&from.node())
                .map_or(&[][..], |v| v.as_slice())
            {
                let Some((start, end)) = ogeom_algo::edge_vertices(model, held)? else {
                    continue;
                };
                let ends = [start.node(), end.node()];
                if !(ends.contains(&from.node()) && ends.contains(&to.node())) {
                    continue;
                }
                if lies_along(model, held, &curve, range, tolerance, tol)? {
                    along.push(held.node());
                }
            }
            let target = match along.as_slice() {
                [] => {
                    // A new edge, across the face the seat lies on: only one
                    // the seat shares with the blend.
                    if !blend_edges.contains(&id) {
                        return Ok(None);
                    }
                    let Some((first, last)) = ogeom_algo::edge_vertices(
                        model,
                        &walked.edge.oriented(Orientation::Forward),
                    )?
                    else {
                        return Ok(None);
                    };
                    new_edges.push(NewEdge {
                        curve,
                        range,
                        from: vertices[&first.node()].clone(),
                        to: vertices[&last.node()].clone(),
                        tolerance: Tolerance::new(tolerance.max(tol.confusion()))?,
                    });
                    Target::New(new_edges.len() - 1)
                }
                [one] => Target::Held(*one),
                _ => return Ok(None),
            };
            targets.insert(id, target);
        }
    }
    if blend_edges.iter().any(|id| !targets.contains_key(id)) {
        return Ok(None);
    }
    // Each seat on the face of the solid it lies on.
    let mut seats = Vec::with_capacity(planar.len());
    for (seat, ring, _, outward) in &planar {
        let Some(placed) = seat_on(model, solid, seat, ring, *outward, &targets, &vertices, tol)?
        else {
            return Ok(None);
        };
        seats.push(placed);
    }
    let hosts: Vec<usize> = seats.iter().map(|s| s.host).collect();
    // The convex region the wedge lies in: behind each seat's plane and the
    // chord plane through the blend's rulings.
    let mut walls: Vec<(Point, Vector)> = planar.iter().map(|p| (p.2, p.3)).collect();
    let chord = {
        let [l0, l1] = cylinder.rulings;
        let normal = cylinder.axis.cross(l1 - l0);
        let m = normal.magnitude();
        if m <= tol.confusion() {
            return Ok(None);
        }
        let normal = normal / m;
        // Turned away from the wedge's vertex furthest off the plane.
        let far = points
            .values()
            .map(|(p, _)| (*p - l0).dot(normal))
            .fold(0.0_f64, |m, d| if d.abs() > m.abs() { d } else { m });
        if far.abs() <= wedge_tolerance {
            return Ok(None);
        }
        (l0, if far > 0.0 { -normal } else { normal })
    };
    walls.push(chord);
    for &(origin, normal) in &walls {
        for (p, _) in points.values() {
            if (*p - origin).dot(normal) > wedge_tolerance {
                return Ok(None);
            }
        }
        for arc in &cylinder.arcs {
            let Some(top) = arc_reach(arc, origin, normal, tol)? else {
                return Ok(None);
            };
            if top > wedge_tolerance {
                return Ok(None);
            }
        }
    }
    // Every other face of the solid stands beyond one of those planes.
    for i in 0..solid.faces.len() {
        if hosts.contains(&i) {
            continue;
        }
        let (bounds, stated) = solid.face_box(model, i, tol)?;
        let margin = stated + tol.confusion();
        let corners = bounds.corners();
        let beyond = walls
            .iter()
            .any(|&(origin, normal)| corners.iter().all(|c| (*c - origin).dot(normal) >= -margin));
        if !beyond {
            return Ok(None);
        }
    }
    Ok(Some(Plan {
        seats,
        blend,
        blend_loop,
        targets,
        vertices,
    }))
}

/// Split the edge of the solid a point of a wedge stands on inside, where
/// it stands within both tolerances of one edge, a line or an arc, away
/// from its ends, and every image of the edge runs over the range its
/// curve does: the two pieces, each the edge over its part of that range,
/// take its place in every loop walking it. The vertex between them, or
/// `None` where the point stands on no edge so.
fn split_at(
    model: &mut Model,
    solid: &mut Solid,
    p: Point,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Shape>> {
    let mut found: Option<(Shape, f64)> = None;
    for edge in solid.edges.values() {
        let Some((curve, range, tolerance)) = curve_of(model, edge) else {
            continue;
        };
        let within = reach + tolerance + tol.confusion();
        let (s, e) = (curve.point_at(range.0, tol)?, curve.point_at(range.1, tol)?);
        let t = match &curve {
            Curve::Line(_) => {
                if segment_distance(p, s, e) > within {
                    continue;
                }
                let d = e - s;
                range.0 + (range.1 - range.0) * ((p - s).dot(d) / d.dot(d))
            }
            Curve::Circle(c) => {
                let circle = c.circle();
                let local = circle.frame().to_local(p);
                if (local.x.hypot(local.y) - circle.radius()).hypot(local.z) > within
                    || !in_arc(&curve, range, p, tol)?
                {
                    continue;
                }
                let centre = circle.centre();
                let x = curve.point_at(0.0, tol)? - centre;
                let y = curve.point_at(core::f64::consts::FRAC_PI_2, tol)? - centre;
                let d = p - centre;
                let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
                let t = d.dot(y).atan2(d.dot(x));
                let t = t + ((lo - t) / core::f64::consts::TAU).ceil() * core::f64::consts::TAU;
                if t > hi {
                    t - core::f64::consts::TAU
                } else {
                    t
                }
            }
            _ => continue,
        };
        let Some((start, end)) = ogeom_algo::edge_vertices(model, edge)? else {
            continue;
        };
        // Away from the ends: a point within reach of one is that vertex.
        let mut clear = true;
        for v in [&start, &end] {
            let Some((q, own)) = vertex_at(model, v) else {
                return Ok(None);
            };
            clear &= p.distance(q) > reach + own;
        }
        if !clear {
            continue;
        }
        if found.is_some() {
            return Ok(None);
        }
        found = Some((edge.clone(), t));
    }
    let Some((edge, t)) = found else {
        return Ok(None);
    };
    let Some(node) = model.node(&edge) else {
        return Ok(None);
    };
    let Some(data) = node.data().as_edge().cloned() else {
        return Ok(None);
    };
    let bounds = node.children().to_vec();
    let Some(EdgeRepr::Curve3d { range, .. }) = data.curve3d().cloned() else {
        return Ok(None);
    };
    if !(range.0.min(range.1) < t && t < range.0.max(range.1)) || bounds.len() != 2 {
        return Ok(None);
    }
    let Some((curve, _, tolerance)) = curve_of(model, &edge) else {
        return Ok(None);
    };
    // Each piece is the edge over its part of the range, every image with
    // it where each runs over the curve's own range.
    let piece = |over: (f64, f64)| -> Option<ogeom_topo::EdgeData> {
        let mut data = data.clone();
        data.representations.retain(|r| {
            matches!(
                r,
                EdgeRepr::Curve3d { .. } | EdgeRepr::PCurve { .. } | EdgeRepr::Seam { .. }
            )
        });
        for r in &mut data.representations {
            let (EdgeRepr::Curve3d { range: own, .. }
            | EdgeRepr::PCurve { range: own, .. }
            | EdgeRepr::Seam { range: own, .. }) = r
            else {
                continue;
            };
            if *own != range {
                return None;
            }
            *own = over;
        }
        Some(data)
    };
    let (Some(first), Some(second)) = (piece((range.0, t)), piece((t, range.1))) else {
        return Ok(None);
    };
    let at = curve.point_at(t, tol)?;
    let vertex = model.add_vertex(ogeom_topo::VertexData {
        point: at,
        tolerance: Tolerance::new(tolerance.max(reach + p.distance(at)).max(tol.confusion()))?,
    });
    let first = model.add_edge(first, &[bounds[0].clone(), vertex.clone()])?;
    let second = model.add_edge(second, &[vertex.clone(), bounds[1].clone()])?;
    let start = bounds[0].node();
    // In every loop walking the edge, the two pieces in its place, in the
    // sense the loop walks it.
    let walkers = solid.walkers.remove(&edge.node()).unwrap_or_default();
    for &f in &walkers {
        for ring in &mut solid.loops[f] {
            let mut out = Vec::with_capacity(ring.len() + 1);
            for walked in ring.drain(..) {
                if walked.edge.node() != edge.node() {
                    out.push(walked);
                } else if walked.edge.orientation() == Orientation::Forward {
                    out.push(Walked {
                        edge: first.clone(),
                        from: walked.from,
                        to: vertex.node(),
                    });
                    out.push(Walked {
                        edge: second.clone(),
                        from: vertex.node(),
                        to: walked.to,
                    });
                } else {
                    out.push(Walked {
                        edge: second.reversed(),
                        from: walked.from,
                        to: vertex.node(),
                    });
                    out.push(Walked {
                        edge: first.reversed(),
                        from: vertex.node(),
                        to: walked.to,
                    });
                }
            }
            *ring = out;
        }
        solid.split.insert(f);
    }
    solid.edges.remove(&edge.node());
    for (piece, ends) in [
        (&first, [start, vertex.node()]),
        (&second, [vertex.node(), bounds[1].node()]),
    ] {
        solid.walkers.insert(piece.node(), walkers.clone());
        solid.edges.insert(piece.node(), piece.clone());
        for v in ends {
            let list = solid.by_vertex.entry(v).or_default();
            list.retain(|e| e.node() != edge.node());
            list.push(piece.clone());
        }
    }
    solid
        .vertices
        .push((vertex.clone(), at, tolerance.max(reach + p.distance(at))));
    solid.vertex_shapes.insert(vertex.node(), vertex.clone());
    Ok(Some(vertex))
}

/// A blend that is a strip of a cylinder between two of its rulings and two
/// of its cross-sections: the axis direction, a point of each ruling, and
/// the two arcs.
struct Strip {
    axis: Vector,
    rulings: [Point; 2],
    arcs: Vec<(Curve, (f64, f64))>,
}

/// The blend as a [`Strip`], where its surface is a cylinder at no
/// placement and its loop two lines along the axis on the cylinder and two
/// arcs of its cross-sections, each under half a turn.
fn blend_cylinder(
    model: &Model,
    blend: &Shape,
    ring: &[Walked],
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
    let Some(NodeData::Face(data)) = model.node(blend).map(|n| n.data()) else {
        return Ok(None);
    };
    if !data.location.is_identity() || !blend.location().is_identity() {
        return Ok(None);
    }
    let Some(SurfaceGeometry::Cylinder(surface)) = model.geometry().surface(data.surface) else {
        return Ok(None);
    };
    let cylinder = surface.cylinder();
    let frame = cylinder.frame();
    let axis = frame.z().vector();
    let radius = cylinder.radius();
    if ring.len() != 4 {
        return Ok(None);
    }
    let mut rulings = Vec::new();
    let mut arcs = Vec::new();
    for walked in ring {
        let Some((curve, range, tolerance)) = curve_of(model, &walked.edge) else {
            return Ok(None);
        };
        let reach = tolerance.max(tol.confusion());
        match &curve {
            Curve::Line(_) => {
                let (p, q) = (curve.point_at(range.0, tol)?, curve.point_at(range.1, tol)?);
                let along = q - p;
                let length = along.magnitude();
                if length <= tol.confusion()
                    || along.cross(axis).magnitude() > length * tol.angular()
                {
                    return Ok(None);
                }
                if (cylinder.distance_to(p)).abs() > reach
                    || (cylinder.distance_to(q)).abs() > reach
                {
                    return Ok(None);
                }
                rulings.push(p);
            }
            Curve::Circle(c) => {
                let circle = c.circle();
                let normal = circle.frame().z().vector();
                let off = circle.centre() - frame.origin();
                if normal.cross(axis).magnitude() > tol.angular()
                    || off.cross(axis).magnitude() > reach
                    || (circle.radius() - radius).abs() > reach
                {
                    return Ok(None);
                }
                let sweep = (range.1 - range.0).abs();
                if sweep <= tol.angular() || sweep >= core::f64::consts::PI {
                    return Ok(None);
                }
                arcs.push((curve.clone(), range));
            }
            _ => return Ok(None),
        }
    }
    let ([l0, l1], 2) = (
        match rulings.as_slice() {
            [p, q] => [*p, *q],
            _ => return Ok(None),
        },
        arcs.len(),
    ) else {
        return Ok(None);
    };
    Ok(Some(Strip {
        axis,
        rulings: [l0, l1],
        arcs,
    }))
}

/// How far an arc reaches past a plane (a point of it and its outward unit
/// normal): the greatest signed distance of its points. `None` where the
/// plane runs along the axis of the arc's circle and the arc meets it at a
/// point inside its range: a blend touching a plane of its region inside
/// its own boundary.
fn arc_reach(
    (curve, range): &(Curve, (f64, f64)),
    origin: Point,
    normal: Vector,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let Curve::Circle(c) = curve else {
        return Ok(None);
    };
    let centre = c.circle().centre();
    // A point of a circle is its centre and a sinusoid in the angle: read
    // off the curve at two angles a quarter turn apart.
    let x = curve.point_at(0.0, tol)? - centre;
    let y = curve.point_at(core::f64::consts::FRAC_PI_2, tol)? - centre;
    let (base, bx, by) = ((centre - origin).dot(normal), x.dot(normal), y.dot(normal));
    let at = |t: f64| base + bx * t.cos() + by * t.sin();
    let (lo, hi) = if range.0 <= range.1 {
        (range.0, range.1)
    } else {
        (range.1, range.0)
    };
    let mut top = at(lo).max(at(hi));
    if bx.hypot(by) > tol.confusion() {
        let peak = by.atan2(bx);
        let turns = ((lo - peak) / core::f64::consts::TAU).ceil();
        let t = peak + turns * core::f64::consts::TAU;
        if t < hi {
            let along_axis = c.circle().frame().z().vector().dot(normal).abs() <= tol.angular();
            if along_axis && t > lo + tol.angular() && t < hi - tol.angular() {
                return Ok(None);
            }
            top = top.max(at(t));
        }
    }
    Ok(Some(top))
}

/// Whether an edge of the solid, from one end of a wedge edge to the other,
/// lies along it: its points within both edges' tolerances of the wedge
/// edge's curve over its range.
fn lies_along(
    model: &Model,
    held: &Shape,
    curve: &Curve,
    range: (f64, f64),
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let Some((own, own_range, own_tolerance)) = curve_of(model, held) else {
        return Ok(false);
    };
    let reach = tolerance + own_tolerance + tol.confusion();
    for k in 1..ALONG_SAMPLES {
        let t = own_range.0 + (own_range.1 - own_range.0) * f64::from(k) / f64::from(ALONG_SAMPLES);
        let p = own.point_at(t, tol)?;
        let off = match curve {
            Curve::Line(_) => {
                let (s, e) = (curve.point_at(range.0, tol)?, curve.point_at(range.1, tol)?);
                segment_distance(p, s, e)
            }
            Curve::Circle(c) => {
                let circle = c.circle();
                let local = circle.frame().to_local(p);
                let off = (local.x.hypot(local.y) - circle.radius()).hypot(local.z);
                if !in_arc(curve, range, p, tol)? {
                    return Ok(false);
                }
                off
            }
            _ => return Ok(false),
        };
        if off > reach {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether a point by a circle curve stands within the range the edge
/// uses: the angle the curve reads it at, a whole turn either way, lies in
/// the range.
fn in_arc(curve: &Curve, range: (f64, f64), p: Point, tol: Tolerances) -> OgeomResult<bool> {
    let Curve::Circle(c) = curve else {
        return Ok(false);
    };
    let (lo, hi) = if range.0 <= range.1 {
        (range.0, range.1)
    } else {
        (range.1, range.0)
    };
    // The curve is its centre and a sinusoid in its parameter, read off at
    // two parameters a quarter turn apart, whichever way it runs.
    let centre = c.circle().centre();
    let x = curve.point_at(0.0, tol)? - centre;
    let y = curve.point_at(core::f64::consts::FRAC_PI_2, tol)? - centre;
    let d = p - centre;
    let t = d.dot(y).atan2(d.dot(x));
    let t = t + ((lo - t) / core::f64::consts::TAU - tol.angular()).ceil() * core::f64::consts::TAU;
    Ok(t <= hi + tol.angular())
}

/// The distance from a point to a segment.
fn segment_distance(p: Point, s: Point, e: Point) -> f64 {
    let d = e - s;
    let length2 = d.dot(d);
    if length2 <= 0.0 {
        return p.distance(s);
    }
    let t = ((p - s).dot(d) / length2).clamp(0.0, 1.0);
    p.distance(s + d * t)
}

/// A seat on the face of the solid it lies on, `None` where it lies on none
/// as the module asks.
#[allow(clippy::too_many_arguments, reason = "one question, all its data")]
fn seat_on(
    model: &Model,
    solid: &Solid,
    seat: &Shape,
    ring: &[Walked],
    outward: Vector,
    targets: &FastMap<TShapeId, Target>,
    vertices: &FastMap<TShapeId, Shape>,
    tol: Tolerances,
) -> OgeomResult<Option<Seat>> {
    let n = ring.len();
    let held_at: Vec<Option<TShapeId>> = ring
        .iter()
        .map(|w| match targets[&w.edge.node()] {
            Target::Held(id) => Some(id),
            Target::New(_) => None,
        })
        .collect();
    let Some(any) = held_at.iter().flatten().next().copied() else {
        return Ok(None);
    };
    // The seat walked from its first held edge after a new one: the held
    // run, then the new edges, each one run.
    let start = (0..n)
        .find(|&k| held_at[k].is_some() && held_at[(k + n - 1) % n].is_none())
        .unwrap_or(0);
    let order: Vec<usize> = (0..n).map(|k| (start + k) % n).collect();
    let held_len = order.iter().take_while(|&&k| held_at[k].is_some()).count();
    if order[held_len..].iter().any(|&k| held_at[k].is_some()) {
        return Ok(None);
    }
    // The face of the solid walking the first held edge in the seat's plane,
    // facing the seat's way.
    let mut hosts = solid
        .walkers
        .get(&any)
        .map_or(&[][..], |v| v.as_slice())
        .iter()
        .copied()
        .filter(|&f| {
            plane_of(model, &solid.faces[f])
                .is_some_and(|(_, _, normal)| (normal - outward).magnitude() <= tol.angular())
        });
    let (Some(host), None) = (hosts.next(), hosts.next()) else {
        return Ok(None);
    };
    // The held run is a run of one of the host's loops, walked the same way.
    let first = &ring[order[0]];
    let mut found = None;
    for (l, face_ring) in solid.loops[host].iter().enumerate() {
        for (p, walked) in face_ring.iter().enumerate() {
            if Some(walked.edge.node()) == held_at[order[0]]
                && walked.from == vertices[&first.from].node()
            {
                found = Some((l, p));
            }
        }
    }
    let Some((l, p)) = found else {
        return Ok(None);
    };
    let face_ring = &solid.loops[host][l];
    if held_len > face_ring.len() {
        return Ok(None);
    }
    for (k, &at) in order.iter().take(held_len).enumerate() {
        let there = &face_ring[(p + k) % face_ring.len()];
        let mine = &ring[at];
        if Some(there.edge.node()) != held_at[at]
            || there.from != vertices[&mine.from].node()
            || there.to != vertices[&mine.to].node()
        {
            return Ok(None);
        }
    }
    let free: Vec<Walked> = order[held_len..].iter().map(|&k| ring[k].clone()).collect();
    if free.is_empty() {
        // The seat is the whole face, which has no hole.
        if held_len != face_ring.len() || solid.loops[host].len() != 1 {
            return Ok(None);
        }
        return Ok(Some(Seat {
            seat: seat.clone(),
            host,
            run: None,
            free,
        }));
    }
    // The new edges cross the face: they meet no edge of it off the run but
    // at their ends, and no hole of it stands inside the seat.
    let seat_box = ogeom_algo::shape_bounds(model, seat, tol)?;
    let run: Vec<TShapeId> = (0..held_len)
        .map(|k| face_ring[(p + k) % face_ring.len()].edge.node())
        .collect();
    for (other_l, other) in solid.loops[host].iter().enumerate() {
        if other_l != l {
            let mut ring_box = Aabb::EMPTY;
            for walked in other {
                ring_box = ring_box.union(&ogeom_algo::shape_bounds(model, &walked.edge, tol)?);
            }
            if ring_box.intersects(&seat_box) {
                return Ok(None);
            }
        }
        for walked in other {
            if run.contains(&walked.edge.node()) {
                continue;
            }
            for crossing in &free {
                if !apart(model, &walked.edge, crossing, vertices, tol)? {
                    return Ok(None);
                }
            }
        }
    }
    Ok(Some(Seat {
        seat: seat.clone(),
        host,
        run: Some((l, p, held_len)),
        free,
    }))
}

/// Whether an edge of a face and a new edge across it meet nowhere but at
/// the vertices they share: their boxes apart, or two segments, or a
/// segment and an arc, whose only common points are those vertices.
fn apart(
    model: &Model,
    edge: &Shape,
    crossing: &Walked,
    vertices: &FastMap<TShapeId, Shape>,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let (Some((own, own_range, own_tolerance)), Some((new, new_range, new_tolerance))) =
        (curve_of(model, edge), curve_of(model, &crossing.edge))
    else {
        return Ok(false);
    };
    let reach = own_tolerance + new_tolerance + crate::sub_shapes::reach(tol);
    let a = ogeom_algo::measure::curve_bounds_over(&own, own_range, tol)?.expanded(reach);
    let b = ogeom_algo::measure::curve_bounds_over(&new, new_range, tol)?.expanded(reach);
    if !a.intersects(&b) {
        return Ok(true);
    }
    let Some((start, end)) = ogeom_algo::edge_vertices(model, edge)? else {
        return Ok(false);
    };
    let new_ends = [
        vertices[&crossing.from].node(),
        vertices[&crossing.to].node(),
    ];
    let shared: Vec<Point> = [&start, &end]
        .into_iter()
        .filter(|v| new_ends.contains(&v.node()))
        .filter_map(|v| vertex_at(model, v).map(|(p, _)| p))
        .collect();
    match (&own, &new) {
        (Curve::Line(_), Curve::Line(_)) => {
            let p = [
                own.point_at(own_range.0, tol)?,
                own.point_at(own_range.1, tol)?,
            ];
            let q = [
                new.point_at(new_range.0, tol)?,
                new.point_at(new_range.1, tol)?,
            ];
            match shared.as_slice() {
                [] => Ok(segments_gap(p, q) > reach),
                // Two segments from one point meet only there unless they
                // leave it the same way.
                [at] => {
                    let away = |s: [Point; 2]| {
                        if s[0].distance(*at) <= s[1].distance(*at) {
                            s[1] - s[0]
                        } else {
                            s[0] - s[1]
                        }
                    };
                    let (u, v) = (away(p), away(q));
                    let (mu, mv) = (u.magnitude(), v.magnitude());
                    Ok(mu > tol.confusion()
                        && mv > tol.confusion()
                        && u.dot(v) / (mu * mv) < 1.0 - tol.angular())
                }
                _ => Ok(false),
            }
        }
        (Curve::Line(_), Curve::Circle(_)) => {
            segment_clear_of_arc(&own, own_range, (&new, new_range), &shared, reach, tol)
        }
        (Curve::Circle(_), Curve::Line(_)) => {
            segment_clear_of_arc(&new, new_range, (&own, own_range), &shared, reach, tol)
        }
        _ => Ok(false),
    }
}

/// Whether a segment and an arc lying in one plane have no common point
/// but the vertices `shared` they end on together: past those, the line
/// crosses the arc's circle nowhere on the arc, and where it misses the
/// circle it stands further than `reach` from the arc.
fn segment_clear_of_arc(
    line: &Curve,
    range: (f64, f64),
    (arc, arc_range): (&Curve, (f64, f64)),
    shared: &[Point],
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let Curve::Circle(c) = arc else {
        return Ok(false);
    };
    let circle = c.circle();
    let (centre, radius) = (circle.centre(), circle.radius());
    let (s, e) = (line.point_at(range.0, tol)?, line.point_at(range.1, tol)?);
    let along = e - s;
    let length = along.magnitude();
    if length <= tol.confusion()
        || along.dot(circle.frame().z().vector()).abs() > length * tol.angular()
    {
        return Ok(false);
    }
    let d = along / length;
    // Where the line meets the circle, as distances from `s` along it.
    let off = s - centre;
    let (half_b, c0) = (d.dot(off), off.dot(off) - radius * radius);
    let disc = half_b * half_b - c0;
    let mut meets: Vec<f64> = Vec::new();
    if disc >= 0.0 {
        let root = disc.sqrt();
        meets.extend([-half_b - root, -half_b + root]);
    }
    let at_shared = |u: f64| {
        let q = s + d * u;
        shared.iter().any(|v| v.distance(q) <= reach)
    };
    // An end of the arc off the shared vertices stands clear of the segment.
    for t in [arc_range.0, arc_range.1] {
        let q = arc.point_at(t, tol)?;
        if !shared.iter().any(|v| v.distance(q) <= reach) && segment_distance(q, s, e) <= reach {
            return Ok(false);
        }
    }
    for &u in &meets {
        if u < -reach || u > length + reach || at_shared(u) {
            continue;
        }
        if in_arc(arc, arc_range, s + d * u, tol)? {
            return Ok(false);
        }
    }
    if shared.is_empty() && meets.is_empty() {
        // The line misses the circle: the arc stands off it no nearer than
        // the circle does, at the foot of the centre on the line.
        let foot = (-half_b).clamp(0.0, length);
        let gap = (s + d * foot).distance(centre) - radius;
        return Ok(gap > reach);
    }
    Ok(true)
}

/// The least distance between two segments.
fn segments_gap(p: [Point; 2], q: [Point; 2]) -> f64 {
    let d1 = p[1] - p[0];
    let d2 = q[1] - q[0];
    let r = p[0] - q[0];
    let (a, e, f) = (d1.dot(d1), d2.dot(d2), d2.dot(r));
    let (s, t) = if a <= f64::EPSILON && e <= f64::EPSILON {
        (0.0, 0.0)
    } else if a <= f64::EPSILON {
        (0.0, (f / e).clamp(0.0, 1.0))
    } else {
        let c = d1.dot(r);
        if e <= f64::EPSILON {
            ((-c / a).clamp(0.0, 1.0), 0.0)
        } else {
            let b = d1.dot(d2);
            let denom = a * e - b * b;
            let mut s = if denom > 0.0 {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let mut t = (b * s + f) / e;
            if t < 0.0 {
                t = 0.0;
                s = (-c / a).clamp(0.0, 1.0);
            } else if t > 1.0 {
                t = 1.0;
                s = ((b - c) / a).clamp(0.0, 1.0);
            }
            (s, t)
        }
    };
    (p[0] + d1 * s).distance(q[0] + d2 * t)
}
