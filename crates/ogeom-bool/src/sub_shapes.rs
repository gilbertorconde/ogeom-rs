//! The history a boolean records for its arguments' edges and vertices.
//!
//! The boolean's own steps record faces: each face of either argument is
//! copied, modified into its pieces, or deleted. Edges and vertices are read
//! off the result against those records. One the result holds as it stands
//! (an edge of a face set aside) has no record: it traces to itself, which
//! the result holds. Any other edge's images are the result's edges lying
//! along it, among the edges of what its faces became; a vertex's images
//! are the vertices of its edges' images standing where it stands. An edge
//! rebuilt whole as one edge is its exact copy, and an edge or a vertex
//! with no image is deleted.

use ogeom_algo::{Built, History};
use ogeom_core::{FastMap, FastSet, OgeomResult, Tolerances};
use ogeom_geom::Curve3d as _;
use ogeom_geom::{Curve, Transformable as _, TrimmedCurve};
use ogeom_math::{Aabb, Point};
use ogeom_topo::{EdgeRepr, Filter, Model, SameKey, Shape, ShapeType, explore, explore_unique};
use smallvec::SmallVec;

/// Samples taken along an edge for its box and for a point's foot on it.
const SAMPLES: usize = 16;

/// `built` with every edge and vertex of `operands` the result does not
/// hold as it stands recorded against the result: copied where one exact
/// copy of it stands there, modified into its images where it was cut
/// short, split or rebuilt, and deleted where nothing of it is left.
///
/// The work is that of the faces the boolean changed: a face the result
/// holds as it stands holds its edges and vertices with it, and none of
/// them is read.
pub(crate) fn with_sub_shapes(
    model: &Model,
    operands: [&Shape; 2],
    built: Built,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let Built {
        shape: result,
        mut history,
    } = built;
    let held_faces: FastSet<SameKey> = explore(model, &result, Filter::OfType(ShapeType::Face))?
        .into_iter()
        .map(SameKey)
        .collect();
    // Every edge and vertex of the result, read only where an argument is a
    // sheet: its open boundary is walked once, where a solid's every edge
    // is walked twice.
    let mut held_parts: Option<FastSet<SameKey>> = None;
    let mut reads = Reads::default();
    for operand in operands {
        let sheet = !crate::is_solid_or_lumps(model, operand)?;
        if sheet && held_parts.is_none() {
            let mut parts = FastSet::default();
            for kind in [ShapeType::Edge, ShapeType::Vertex] {
                parts.extend(
                    explore(model, &result, Filter::OfType(kind))?
                        .into_iter()
                        .map(SameKey),
                );
            }
            held_parts = Some(parts);
        }
        let touched: Vec<Shape> = explore_unique(model, operand, ShapeType::Face)?
            .into_iter()
            .filter(|face| !held_faces.contains(&SameKey(face.clone())))
            .collect();
        if touched.is_empty() {
            continue;
        }
        // Each edge of the faces changed, with how often they walk it and
        // which of them do. A solid's boundary walks every edge an even
        // number of times, so an edge they walk an odd number of times is
        // walked by a face the result holds as well, and the result holds
        // it.
        let mut edges: Vec<(Shape, usize, SmallVec<[usize; 2]>)> = Vec::new();
        let mut slot: FastMap<SameKey, usize> = FastMap::default();
        for (index, face) in touched.iter().enumerate() {
            for edge in explore(model, face, Filter::OfType(ShapeType::Edge))? {
                let key = SameKey(edge);
                let at = match slot.get(&key) {
                    Some(&at) => at,
                    None => {
                        edges.push((key.0.clone(), 0, SmallVec::new()));
                        slot.insert(key, edges.len() - 1);
                        edges.len() - 1
                    }
                };
                edges[at].1 += 1;
                if !edges[at].2.contains(&index) {
                    edges[at].2.push(index);
                }
            }
        }
        // A degenerate edge (a pole, an apex) is walked once by the one
        // face it bounds.
        let held = |shape: &Shape, walks: usize| match &held_parts {
            Some(parts) if sheet => parts.contains(&SameKey(shape.clone())),
            _ => {
                walks % 2 == 1
                    && !model
                        .node(shape)
                        .and_then(|n| n.data().as_edge())
                        .is_some_and(|d| d.degenerate)
            }
        };
        // The edges of what each changed face became, read where asked.
        let mut near: Vec<Option<Vec<Shape>>> = vec![None; touched.len()];
        // Each vertex of the edges read, with whether the result holds it
        // as it stands, whether every one came through as a copy, and
        // their images.
        let mut ends: FastMap<SameKey, (Shape, bool, bool, Vec<Shape>)> = FastMap::default();
        let mut kept_edges: Vec<&Shape> = Vec::new();
        for (edge, walks, faces) in &edges {
            if held(edge, *walks) {
                kept_edges.push(edge);
                continue;
            }
            for &face in faces {
                if near[face].is_none() {
                    let mut found = Vec::new();
                    for image in history.trace(&touched[face]) {
                        if held_faces.contains(&SameKey(image.clone())) {
                            found.extend(explore_unique(model, image, ShapeType::Edge)?);
                        }
                    }
                    near[face] = Some(found);
                }
            }
            let candidates = faces.iter().filter_map(|&f| near[f].as_ref()).flatten();
            let (images, copied) = images_of(model, edge, candidates, &slot, &mut reads, tol)?;
            // The rebuild keeps an edge it shares with a face set aside as
            // it stands.
            let kept = images.iter().any(|image| image.is_same(edge));
            if kept {
                kept_edges.push(edge);
            } else {
                record(&mut history, edge, &images, copied);
            }
            for vertex in vertices_of(model, edge) {
                let entry = ends
                    .entry(SameKey(vertex.clone()))
                    .or_insert_with(|| (vertex, false, true, Vec::new()));
                entry.2 &= copied;
                entry.3.extend(images.iter().cloned());
            }
        }
        if ends.is_empty() {
            continue;
        }
        // A vertex of an edge the result holds is held with it.
        for edge in kept_edges {
            for vertex in vertices_of(model, edge) {
                if let Some(entry) = ends.get_mut(&SameKey(vertex)) {
                    entry.1 = true;
                }
            }
        }
        for (vertex, kept, copied, edges) in ends.into_values() {
            if kept || (sheet && held(&vertex, 0)) {
                continue;
            }
            let Some((at, own)) = vertex_at(model, &vertex)? else {
                continue;
            };
            let mut images: Vec<Shape> = Vec::new();
            for edge in &edges {
                for end in &reads.edge(model, edge, tol)?.ends {
                    if end.0.distance(at) <= own + end.1 + reach(tol)
                        && !images.iter().any(|v| v.is_same(&end.2))
                    {
                        images.push(end.2.clone());
                    }
                }
            }
            if !images.iter().any(|image| image.is_same(&vertex)) {
                record(&mut history, &vertex, &images, copied);
            }
        }
    }
    Ok(Built::new(result, history))
}

/// An edge's vertices, placed as it stands.
fn vertices_of<'m>(model: &'m Model, edge: &'m Shape) -> impl Iterator<Item = Shape> + 'm {
    model
        .node(edge)
        .map_or(&[][..], |n| n.children())
        .iter()
        .map(|child| child.beneath(edge.location(), edge.orientation()))
}

/// Record `part` against its images: deleted where it has none, an exact
/// copy where `copied` and it has one, modified into each otherwise. What
/// the history said of it before is replaced.
fn record(history: &mut History, part: &Shape, images: &[Shape], copied: bool) {
    history.delete(part);
    match images {
        [] => {}
        [image] if copied => history.copy(part, image.clone()),
        _ => {
            for image in images {
                history.modify(part, image.clone());
            }
        }
    }
}

/// How far two descriptions of one point may stand apart past their own
/// tolerances: the distance the boolean welds ends within.
fn reach(tol: Tolerances) -> f64 {
    tol.confusion() * 1e2
}

/// A vertex in space with its tolerance.
fn vertex_at(model: &Model, vertex: &Shape) -> OgeomResult<Option<(Point, f64)>> {
    let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
        return Ok(None);
    };
    Ok(Some((
        vertex.transform(model.datums())?.apply(data.point),
        data.tolerance.get(),
    )))
}

/// A result edge as the comparisons read it.
struct Read {
    /// Each end in space, with its vertex's tolerance and the vertex.
    ends: Vec<(Point, f64, Shape)>,
    /// The curve's middle in space, where the edge has a curve.
    middle: Option<Point>,
    /// The edge's own tolerance.
    tolerance: f64,
}

/// The result's edges, each read once.
#[derive(Default)]
struct Reads {
    read: FastMap<SameKey, Read>,
}

impl Reads {
    fn edge(&mut self, model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<&Read> {
        let key = SameKey(edge.clone());
        if !self.read.contains_key(&key) {
            let read = read_edge(model, edge, tol)?;
            self.read.insert(key.clone(), read);
        }
        Ok(&self.read[&key])
    }
}

fn read_edge(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<Read> {
    let node = model.node(edge);
    let mut ends: Vec<(Point, f64, Shape)> = Vec::with_capacity(2);
    for child in node.map_or(&[][..], |n| n.children()) {
        let vertex = child.beneath(edge.location(), edge.orientation());
        if !ends.iter().any(|end| end.2.is_same(&vertex))
            && let Some((at, own)) = vertex_at(model, &vertex)?
        {
            ends.push((at, own, vertex));
        }
    }
    let data = node.and_then(|n| n.data().as_edge());
    let middle = match data.and_then(|d| d.curve3d()) {
        Some(EdgeRepr::Curve3d { curve, range, .. }) if !data.is_some_and(|d| d.degenerate) => {
            match model.geometry().curve(*curve) {
                Some(g) => Some(
                    edge.transform(model.datums())?
                        .apply(g.point_at(f64::midpoint(range.0, range.1), tol)?),
                ),
                None => None,
            }
        }
        _ => None,
    };
    Ok(Read {
        ends,
        middle,
        tolerance: data.map_or(0.0, |d| d.tolerance.get()),
    })
}

/// An argument's edge in space, for asking whether a point lies along it.
struct Along {
    /// The curve restricted to the edge, where it is no straight line.
    curve: Option<Curve>,
    /// The edge's ends along its curve, in its parameter's order.
    span: (Point, Point),
    /// Whether the edge is a straight segment between `span`.
    straight: bool,
    /// A box round the edge, grown by how far the curve may sag between
    /// its samples.
    bound: Aabb,
    tolerance: f64,
}

impl Along {
    fn of(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<Option<Self>> {
        let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
            return Ok(None);
        };
        if data.degenerate {
            return Ok(None);
        }
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            return Ok(None);
        };
        let Some(stored) = model.geometry().curve(*curve) else {
            return Ok(None);
        };
        let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
        if hi - lo <= tol.parametric() {
            return Ok(None);
        }
        let placed = edge.transform(model.datums())?;
        let straight = matches!(stored, Curve::Line(_));
        let span = (
            placed.apply(stored.point_at(lo, tol)?),
            placed.apply(stored.point_at(hi, tol)?),
        );
        if straight {
            return Ok(Some(Self {
                curve: None,
                span,
                straight,
                bound: Aabb::of_corners(span.0, span.1),
                tolerance: data.tolerance.get(),
            }));
        }
        let world = stored.transformed(&placed, tol)?;
        let mut bound = Aabb::EMPTY;
        let mut sag = 0.0_f64;
        let mut previous: Option<(f64, Point)> = None;
        for i in 0..=SAMPLES {
            #[allow(clippy::cast_precision_loss)]
            let t = lo + (hi - lo) * (i as f64 / SAMPLES as f64);
            let p = world.point_at(t, tol)?;
            if let Some((s, q)) = previous {
                let mid = world.point_at(f64::midpoint(s, t), tol)?;
                sag = sag.max(mid.distance(Point::midpoint(q, p)) * 2.0);
            }
            bound = bound.with_point(p);
            previous = Some((t, p));
        }
        let trimmed = TrimmedCurve::new(world, lo, hi, tol)?;
        Ok(Some(Self {
            curve: Some(Curve::Trimmed(Box::new(trimmed))),
            span,
            straight,
            bound: bound.expanded(sag),
            tolerance: data.tolerance.get(),
        }))
    }

    /// Whether `p` lies within `within` of the edge.
    fn holds(&self, p: Point, within: f64, tol: Tolerances) -> OgeomResult<bool> {
        if !self.bound.expanded(within).contains(p) {
            return Ok(false);
        }
        if self.straight {
            let (a, b) = self.span;
            let ab = b - a;
            let t = ((p - a).dot(ab) / ab.dot(ab)).clamp(0.0, 1.0);
            return Ok((a + ab * t).distance(p) <= within);
        }
        let Some(curve) = &self.curve else {
            return Ok(false);
        };
        Ok(ogeom_algo::project_on_curve(curve, p, SAMPLES * 2, tol)?.distance <= within)
    }
}

/// The edges among `near` lying along `edge`, and whether the edge came
/// through whole as the one of them. `siblings` are the edges of the
/// argument's changed faces, by key.
fn images_of<'a>(
    model: &Model,
    edge: &Shape,
    near: impl Iterator<Item = &'a Shape>,
    siblings: &FastMap<SameKey, usize>,
    reads: &mut Reads,
    tol: Tolerances,
) -> OgeomResult<(Vec<Shape>, bool)> {
    let mut images: Vec<Shape> = Vec::new();
    let mut seen: FastSet<SameKey> = FastSet::default();
    let own = read_edge(model, edge, tol)?;
    let along = Along::of(model, edge, tol)?;
    for candidate in near {
        let key = SameKey(candidate.clone());
        // Another edge of the same argument, kept as it stands, lies along
        // no edge of that argument but itself.
        if (siblings.contains_key(&key) && !candidate.is_same(edge)) || !seen.insert(key) {
            continue;
        }
        let read = reads.edge(model, candidate, tol)?;
        let lies = match (&along, read.middle) {
            (Some(along), Some(middle)) => {
                let width = along.tolerance.max(read.tolerance) + reach(tol);
                let mut lies = along.holds(middle, width, tol)?;
                for end in &read.ends {
                    if !lies {
                        break;
                    }
                    lies = along.holds(end.0, width.max(end.1 + reach(tol)), tol)?;
                }
                lies
            }
            // A degenerate edge, or one with no curve: its image is an
            // edge of the same kind ending where it ends.
            (None, None) => {
                !read.ends.is_empty()
                    && read.ends.iter().all(|end| {
                        own.ends
                            .iter()
                            .any(|o| o.0.distance(end.0) <= o.1 + end.1 + reach(tol))
                    })
            }
            _ => false,
        };
        if lies {
            images.push(candidate.clone());
        }
    }
    // Whole: one image ending where the edge ends, at both ends.
    let copied = match images.as_slice() {
        [image] => {
            let read = reads.edge(model, image, tol)?;
            let meets = |a: &[(Point, f64, Shape)], b: &[(Point, f64, Shape)]| {
                a.iter().all(|p| {
                    b.iter()
                        .any(|q| p.0.distance(q.0) <= p.1 + q.1 + reach(tol))
                })
            };
            !own.ends.is_empty() && meets(&own.ends, &read.ends) && meets(&read.ends, &own.ends)
        }
        _ => false,
    };
    Ok((images, copied))
}
