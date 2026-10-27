//! A face's boundary crossing itself in its surface's chart, mended where
//! the mend is the one the boundary itself says.
//!
//! Two mends, each for a face whose boundary crosses itself exactly once,
//! within one wire, in one copy of the chart
//! ([`chart_crossings`](ogeom_algo::chart_crossings)):
//!
//! - Two consecutive edges that cross just short of the vertex they meet
//!   at enclose a loop from the crossing to the vertex and back, wound
//!   against the rest of the face. Where that loop is slop (the crossing
//!   stands from the vertex, and each edge runs past the crossing, no
//!   farther than a fitted trim may stand off its surface), both edges are
//!   cut back to the crossing and the vertex moves there. It widens to
//!   reach where it stood, plus what it stated, so every other curve that
//!   ended within it still does: the same slop, stated the same way, as
//!   two vertices merged across a gap.
//! - Any other crossing makes a figure eight of two loops, one wound each
//!   way. The face's region is both, facing the way the face does, so it
//!   becomes two faces on its surface, one per loop, each wound as an
//!   outer boundary is: the loop wound the other way is walked back. The
//!   two edges crossing are cut there, at one new vertex. A face whose
//!   neighbours walk the loop walked back as it was walked before now
//!   walks those edges the same way they do, and the orientation steps
//!   after this one settle which way that run faces, as for any face
//!   turned against its neighbours. Only a face with one wire (a hole's
//!   loop is not known) and held by a container (a face cannot stand for
//!   two) is split.
//!
//! Every face a mend changes, the face itself and those sharing an edge
//! or the vertex moved with it, must check clean after it and cross
//! itself nowhere in its chart, or the mend is not made. A face one mend
//! changes is not changed by another in the same pass. Anything else is
//! left untouched and reported.

use ogeom_algo::{ChartCrossing, chart_crossings, check, edge_vertices, linear_properties};
use ogeom_core::{FastMap, FastSet, OgeomResult, Tolerances};
use ogeom_geom::{Curve2d as _, Surface as _};
use ogeom_math::{Point, Point2, Vector2};
use ogeom_mesh::Deflection;
use ogeom_topo::{
    EdgeRepr, Model, Orientation, Shape, ShapeType, SurfaceId, TShapeId, VertexData, explore_unique,
};

use crate::Reshape;
use crate::divide::{edge_data, face_data, signed_area, split_edge, vertex_point};

/// How many of each kind of mend were staged.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct Mended {
    /// Pairs of consecutive edges cut back to where they cross.
    pub(crate) trimmed: usize,
    /// Faces split into two at a crossing.
    pub(crate) split: usize,
}

/// One mend, as the substitutions that make it.
struct Mend {
    replace: Vec<(Shape, Shape)>,
    splits: Vec<(Shape, Vec<Shape>)>,
    /// Faces made whole, to be checked as they are.
    made: Vec<Shape>,
    /// Edges and vertices it changes, for the faces they bound.
    touched: Vec<TShapeId>,
    /// Whether it splits the face.
    splits_face: bool,
}

/// Stage the mends of the module documentation over `shape`; `cap` is how
/// far a fitted trim may stand off its surface.
pub(crate) fn stage_crossings(
    model: &mut Model,
    shape: &Shape,
    reshape: &mut Reshape,
    cap: f64,
    tol: Tolerances,
) -> OgeomResult<Mended> {
    let alone = model.kind_of(shape)? == ShapeType::Face;
    let faces = explore_unique(model, shape, ShapeType::Face)?;
    let mut holders: FastMap<TShapeId, Vec<usize>> = FastMap::default();
    for (k, face) in faces.iter().enumerate() {
        for kind in [ShapeType::Edge, ShapeType::Vertex] {
            for sub in explore_unique(model, face, kind)? {
                holders.entry(sub.node()).or_default().push(k);
            }
        }
    }
    let mut changed: FastSet<usize> = FastSet::default();
    let mut mended = Mended::default();
    for (k, face) in faces.iter().enumerate() {
        if changed.contains(&k) {
            continue;
        }
        let crossings = chart_crossings(model, face, tol)?;
        let [crossing] = &crossings[..] else {
            continue;
        };
        let Some(mend) = plan(model, face, crossing, cap, alone, tol)? else {
            continue;
        };
        let mut hit: Vec<usize> = vec![k];
        for id in &mend.touched {
            hit.extend(holders.get(id).into_iter().flatten().copied());
        }
        hit.sort_unstable();
        hit.dedup();
        if hit.iter().any(|h| changed.contains(h)) {
            continue;
        }
        let mut local = Reshape::new();
        stage(&mend, &mut local);
        let mut clean = true;
        let mut changed_faces = mend.made.clone();
        for &h in &hit {
            if h != k || !mend.splits_face {
                changed_faces.push(local.apply(model, &faces[h])?.shape);
            }
        }
        for face in &changed_faces {
            clean &= check(model, face, tol)?.is_valid()
                && chart_crossings(model, face, tol)?.is_empty();
        }
        if !clean {
            continue;
        }
        stage(&mend, reshape);
        changed.extend(hit);
        if mend.splits_face {
            mended.split += 1;
        } else {
            mended.trimmed += 1;
        }
    }
    Ok(mended)
}

/// Stage `mend`'s substitutions in `reshape`.
fn stage(mend: &Mend, reshape: &mut Reshape) {
    for (old, new) in &mend.replace {
        reshape.replace(old, new.clone());
    }
    for (old, pieces) in &mend.splits {
        reshape.split(old, pieces.clone());
    }
}

/// The mend for `face`'s one crossing, where there is one; `None` where
/// the crossing is not one the module documentation mends. `alone` says
/// the face is the shape mended, which a face cannot be split in.
fn plan(
    model: &mut Model,
    face: &Shape,
    crossing: &ChartCrossing,
    cap: f64,
    alone: bool,
    tol: Tolerances,
) -> OgeomResult<Option<Mend>> {
    let data = face_data(model, face)?;
    if !face.location().is_identity() || !data.location.is_identity() || crossing.apart != [0, 0] {
        return Ok(None);
    }
    let [(w, i), (w2, j)] = crossing.edges;
    if w != w2 {
        return Ok(None);
    }
    let stored = face.oriented(Orientation::Forward);
    let wires = model.ordered_children_of(&stored)?;
    let edges = model.ordered_children_of(&wires[w])?;
    let n = edges.len();
    if n < 3 {
        return Ok(None);
    }
    // Each edge crossing walked once by the face, placed as stored, with a
    // curve and two vertices of its own.
    for at in [i, j] {
        let edge = &edges[at];
        let uses = wires
            .iter()
            .map(|wire| model.children_of(wire))
            .collect::<OgeomResult<Vec<_>>>()?
            .iter()
            .flatten()
            .filter(|e| e.node() == edge.node())
            .count();
        let ed = edge_data(model, edge)?;
        let Some((a, b)) = edge_vertices(model, edge)? else {
            return Ok(None);
        };
        if uses != 1
            || !edge.location().is_identity()
            || ed.degenerate
            || ed.curve3d().is_none()
            || a.node() == b.node()
            || !a.location().is_identity()
            || !b.location().is_identity()
        {
            return Ok(None);
        }
    }
    let surface = data.surface;
    let (a, b) = (crossing.parameters[0], crossing.parameters[1]);
    // The edge walked into the vertex they meet at first, where they are
    // consecutive.
    let consecutive = if j == i + 1 {
        Some(((i, a), (j, b)))
    } else if i == 0 && j == n - 1 {
        Some(((j, b), (i, a)))
    } else {
        None
    };
    if let Some(((into, on_into), (out, on_out))) = consecutive
        && let Some(mend) = trim(
            model,
            surface,
            (&edges[into], on_into),
            (&edges[out], on_out),
            crossing.at,
            cap,
            tol,
        )?
    {
        return Ok(Some(mend));
    }
    if alone || wires.len() != 1 {
        return Ok(None);
    }
    split(model, face, &edges, [(i, a), (j, b)], crossing.at, tol)
}

/// `edge` cut at `parameter` of its pcurve on `surface`: its two pieces as
/// stored, the vertex between them, and where that stands. `None` where
/// the parameter is not inside the pcurve's range.
fn cut(
    model: &mut Model,
    surface: SurfaceId,
    edge: &Shape,
    parameter: f64,
    tol: Tolerances,
) -> OgeomResult<Option<([Shape; 2], Shape, Point)>> {
    let data = edge_data(model, edge)?;
    let Some(EdgeRepr::PCurve { range, .. }) =
        data.pcurve_for(surface, &ogeom_topo::Location::identity())
    else {
        return Ok(None);
    };
    let fraction = (parameter - range.0) / (range.1 - range.0);
    if !(fraction > 0.0 && fraction < 1.0) {
        return Ok(None);
    }
    let stored = Shape::of(edge.node());
    let (pieces, vertices) = split_edge(model, &stored, &[fraction], Some(surface), tol)?;
    let ([first, second], [vertex]) = (&pieces[..], &vertices[..]) else {
        return Ok(None);
    };
    let point = vertex_point(model, vertex)?;
    Ok(Some((
        [first.clone(), second.clone()],
        vertex.clone(),
        point,
    )))
}

/// A new edge with `piece`'s geometry, bounded as it is but at `to` where
/// it is bounded at `from`.
fn rebound(model: &mut Model, piece: &Shape, from: &Shape, to: &Shape) -> OgeomResult<Shape> {
    let data = edge_data(model, piece)?;
    let bounds: Vec<Shape> = model
        .children_of(piece)?
        .iter()
        .map(|v| {
            if v.node() == from.node() {
                to.clone()
            } else {
                Shape::of(v.node())
            }
        })
        .collect();
    model.add_edge(data, &bounds)
}

/// The piece of a cut edge the walk meets first and the one it meets
/// second, as the wire walks them.
fn walked(edge: &Shape, pieces: &[Shape; 2]) -> (Shape, Shape) {
    if edge.orientation() == Orientation::Reversed {
        (pieces[1].reversed(), pieces[0].reversed())
    } else {
        (pieces[0].clone(), pieces[1].clone())
    }
}

/// Cut `before` and `after`, consecutive in a wire and crossing at `at`
/// (the point on the surface) at the given parameters of their pcurves,
/// back to the crossing, and move the vertex they share there; `None`
/// where they share no vertex or the loop cut off is not slop.
fn trim(
    model: &mut Model,
    surface: SurfaceId,
    (before, on_before): (&Shape, f64),
    (after, on_after): (&Shape, f64),
    at: Point,
    cap: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Mend>> {
    let (Some((_, joint)), Some((start, _))) =
        (edge_vertices(model, before)?, edge_vertices(model, after)?)
    else {
        return Ok(None);
    };
    if joint.node() != start.node() {
        return Ok(None);
    }
    let stood = vertex_point(model, &joint)?;
    if stood.distance(at) > cap {
        return Ok(None);
    }
    let (Some((cut_before, made_before, end_before)), Some((cut_after, made_after, end_after))) = (
        cut(model, surface, before, on_before, tol)?,
        cut(model, surface, after, on_after, tol)?,
    ) else {
        return Ok(None);
    };
    let (kept_before, lost_before) = walked(before, &cut_before);
    let (lost_after, kept_after) = walked(after, &cut_after);
    for lost in [&lost_before, &lost_after] {
        if linear_properties(model, lost, Deflection::default(), tol)?.mass > cap {
            return Ok(None);
        }
    }
    let stated = |model: &Model, s: &Shape| -> OgeomResult<f64> {
        Ok(model.tolerance_of(s)?.map_or(0.0, |t| t.get()))
    };
    let reach = [
        stated(model, &joint)? + stood.distance(at),
        end_before.distance(at),
        end_after.distance(at),
        stated(model, before)?,
        stated(model, after)?,
        tol.confusion(),
    ]
    .into_iter()
    .fold(0.0, f64::max);
    let moved = model.add_vertex(VertexData::with_tolerance(at, reach)?);
    let stored = |s: &Shape| {
        if s.orientation() == Orientation::Reversed {
            s.reversed()
        } else {
            s.clone()
        }
    };
    let new_before = rebound(model, &stored(&kept_before), &made_before, &moved)?;
    let new_after = rebound(model, &stored(&kept_after), &made_after, &moved)?;
    Ok(Some(Mend {
        replace: vec![
            (Shape::of(before.node()), new_before),
            (Shape::of(after.node()), new_after),
            (Shape::of(joint.node()), moved),
        ],
        splits: Vec::new(),
        made: Vec::new(),
        touched: vec![before.node(), after.node(), joint.node()],
        splits_face: false,
    }))
}

/// Split `face`, whose one wire walks `edges` and crosses itself where the
/// edges at the two indices cross, at the given parameters of their
/// pcurves, at `at`, into the figure eight's two loops, each wound as an
/// outer boundary; `None` where the loops do not close in the chart or are
/// not wound one each way.
fn split(
    model: &mut Model,
    face: &Shape,
    edges: &[Shape],
    [(i, on_first), (j, on_second)]: [(usize, f64); 2],
    at: Point,
    tol: Tolerances,
) -> OgeomResult<Option<Mend>> {
    let data = face_data(model, face)?;
    let surface = data.surface;
    let (first, second) = (&edges[i], &edges[j]);
    let (Some((cut_first, made_first, end_first)), Some((cut_second, made_second, end_second))) = (
        cut(model, surface, first, on_first, tol)?,
        cut(model, surface, second, on_second, tol)?,
    ) else {
        return Ok(None);
    };
    let stated = |model: &Model, s: &Shape| -> OgeomResult<f64> {
        Ok(model.tolerance_of(s)?.map_or(0.0, |t| t.get()))
    };
    let reach = [
        end_first.distance(at),
        end_second.distance(at),
        stated(model, first)?,
        stated(model, second)?,
        tol.confusion(),
    ]
    .into_iter()
    .fold(0.0, f64::max);
    let middle = model.add_vertex(VertexData::with_tolerance(at, reach)?);
    let pieces_first = [
        rebound(model, &cut_first[0], &made_first, &middle)?,
        rebound(model, &cut_first[1], &made_first, &middle)?,
    ];
    let pieces_second = [
        rebound(model, &cut_second[0], &made_second, &middle)?,
        rebound(model, &cut_second[1], &made_second, &middle)?,
    ];
    let (first_in, first_out) = walked(first, &pieces_first);
    let (second_in, second_out) = walked(second, &pieces_second);
    let mut inner = vec![first_out];
    inner.extend(edges[i + 1..j].iter().cloned());
    inner.push(second_in);
    let mut outer = vec![second_out];
    outer.extend(edges[j + 1..].iter().cloned());
    outer.extend(edges[..i].iter().cloned());
    outer.push(first_in);
    let (Some(a), Some(b)) = (
        chart_area(model, surface, &inner, tol)?,
        chart_area(model, surface, &outer, tol)?,
    ) else {
        return Ok(None);
    };
    let (kept, turned) = match (a > 0.0, b > 0.0) {
        (true, false) if b < 0.0 => (inner, outer),
        (false, true) if a < 0.0 => (outer, inner),
        _ => return Ok(None),
    };
    let back: Vec<Shape> = turned.iter().rev().map(Shape::reversed).collect();
    let kept_wire = model.add_wire(&kept)?;
    let back_wire = model.add_wire(&back)?;
    let one = model.add_face(data.clone(), &[kept_wire])?;
    let other = model.add_face(data, &[back_wire])?;
    Ok(Some(Mend {
        replace: Vec::new(),
        splits: vec![
            (Shape::of(first.node()), pieces_first.to_vec()),
            (Shape::of(second.node()), pieces_second.to_vec()),
            (Shape::of(face.node()), vec![one.clone(), other.clone()]),
        ],
        made: vec![one, other],
        touched: vec![first.node(), second.node()],
        splits_face: true,
    }))
}

/// The signed area a loop of edges encloses in `surface`'s chart, positive
/// where it winds as a face's outer boundary does, each pcurve carried
/// across the seam to where the last one ended; `None` where an edge has
/// no plain pcurve there, or the loop does not close in one copy of the
/// chart.
fn chart_area(
    model: &Model,
    surface: SurfaceId,
    edges: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    const SAMPLES: u32 = 32;
    let Some(geometry) = model.geometry().surface(surface) else {
        return Ok(None);
    };
    let ((u0, u1), (v0, v1)) = geometry.domain();
    let periods = [
        (geometry.is_periodic_u() || geometry.is_closed_u(tol)).then_some(u1 - u0),
        (geometry.is_periodic_v() || geometry.is_closed_v(tol)).then_some(v1 - v0),
    ];
    let carried = |from: Point2, to: Point2| {
        let mut by = Vector2::new(0.0, 0.0);
        for (axis, period) in periods.iter().enumerate() {
            let Some(period) = period.filter(|p| *p > 0.0 && p.is_finite()) else {
                continue;
            };
            if axis == 0 {
                by.x = ((to.x - from.x) / period).round() * period;
            } else {
                by.y = ((to.y - from.y) / period).round() * period;
            }
        }
        by
    };
    let mut ring: Vec<Point2> = Vec::new();
    for edge in edges {
        let data = edge_data(model, edge)?;
        let Some(EdgeRepr::PCurve { curve, range, .. }) =
            data.pcurve_for(surface, &ogeom_topo::Location::identity())
        else {
            return Ok(None);
        };
        let Some(pcurve) = model.geometry().pcurve(*curve) else {
            return Ok(None);
        };
        let reversed = edge.orientation() == Orientation::Reversed;
        let mut points = Vec::with_capacity(SAMPLES as usize + 1);
        for k in 0..=SAMPLES {
            let mut s = f64::from(k) / f64::from(SAMPLES);
            if reversed {
                s = 1.0 - s;
            }
            points.push(pcurve.point_at(range.0 + (range.1 - range.0) * s, tol)?);
        }
        let by = match ring.last() {
            Some(&last) => carried(points[0], last),
            None => Vector2::new(0.0, 0.0),
        };
        ring.extend(points.into_iter().map(|p| p + by));
    }
    let (Some(&first), Some(&last)) = (ring.first(), ring.last()) else {
        return Ok(None);
    };
    if carried(first, last) != Vector2::new(0.0, 0.0) {
        return Ok(None);
    }
    Ok(Some(signed_area(&ring)))
}
