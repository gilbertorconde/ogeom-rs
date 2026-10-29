//! Free edges of a shell that are each other's twin across a gap, stated
//! as the gap so that sewing joins them.
//!
//! Sewing takes two edges for one where each lies within the tolerance the
//! other states ([`sew`](ogeom_algo::sew())). Two faces of an open shell
//! built apart (fitted separately, read from separate entities) can meet a
//! little wider than that, each edge free, the shell open between them.
//! Sewing everything at a wider reach joins edges that are not twins at
//! all: on the community assemblies' open sheets, a millimetre joins
//! hundreds of pairs. So the reach is not widened; the pair is.
//!
//! Two free edges of two faces of a shell are twins across a gap only where
//! all of these hold, as for a lacking edge across a wire gap:
//!
//! - their ends stand pairwise, one way round or the other, and the two
//!   curves anywhere along them, no farther apart than a fitted trim may
//!   stand off its surface;
//! - each is the only free edge so near the other;
//! - the gap between the curves, and between each pair of ends, is shorter
//!   than every edge of either face, so no feature of the faces fits in
//!   it;
//! - no other vertex of the shape stands as near any of their ends as the
//!   gap at that end, so no end can be taken for another.
//!
//! Each edge then states the gap between the two curves, and each vertex
//! the gap at its end, where that is wider than it states: the slop a
//! reader would have recorded the same way. Sewing then joins them as it
//! joins any two edges within their tolerances.
//!
//! A pair is stated only where sewing is shown to join it: a copy of the
//! shape is sewn with the pair stated and with none, and the pair must add
//! exactly one join; then all the pairs kept together must add exactly
//! one each, or none is kept. The shape itself is widened only for the
//! pairs kept, so a pair sewing would not join leaves it exactly as it
//! was.

use ogeom_algo::{edge_vertices, edges_apart, linear_properties};
use ogeom_core::{FastMap, OgeomResult, Tolerance, Tolerances};
use ogeom_math::Point;
use ogeom_mesh::Deflection;
use ogeom_topo::{Model, Shape, ShapeType, TShapeId, explore_unique};

/// One free edge, with its face and its ends as that face walks it.
struct Free {
    edge: Shape,
    face: usize,
    ends: [(Shape, Point); 2],
}

/// Two free edges taken for twins across a gap: what each of them and of
/// their ends widens to, the two edges first.
#[derive(Debug, Clone)]
pub(crate) struct Twins {
    widen: Vec<(Shape, f64)>,
}

/// Widen the free edges of `shape`'s shells that are twins across a gap no
/// wider than `cap`, and their ends, to that gap, where sewing is shown to
/// join them (see the module documentation). Returns how many pairs are
/// stated.
pub(crate) fn state_free_gaps(
    model: &mut Model,
    shape: &Shape,
    cap: f64,
    tol: Tolerances,
) -> OgeomResult<usize> {
    let found = twins(model, shape, cap, tol)?;
    state_joining(model, shape, &found, tol)
}

/// Widen `twins` on `shape` where sewing is shown to join them: each pair
/// adds exactly one join to the sewing of a copy of the shape, and the
/// pairs kept together add one each. Returns how many pairs are widened.
pub(crate) fn state_joining(
    model: &mut Model,
    shape: &Shape,
    twins: &[Twins],
    tol: Tolerances,
) -> OgeomResult<usize> {
    if twins.is_empty() {
        return Ok(0);
    }
    let Some(alone) = trial(model, shape, twins, &[], tol)? else {
        return Ok(0);
    };
    let mut kept = Vec::new();
    for k in 0..twins.len() {
        if trial(model, shape, twins, &[k], tol)? == Some(alone + 1) {
            kept.push(k);
        }
    }
    if kept.len() > 1 && trial(model, shape, twins, &kept, tol)? != Some(alone + kept.len()) {
        kept.clear();
    }
    for &k in &kept {
        for (sub, reach) in &twins[k].widen {
            model.widen(sub, Tolerance::new(*reach)?)?;
        }
    }
    Ok(kept.len())
}

/// How many edge pairs sewing joins on a copy of `shape` with the twins
/// `which` widened; `None` where the copy cannot be made out.
fn trial(
    model: &mut Model,
    shape: &Shape,
    twins: &[Twins],
    which: &[usize],
    tol: Tolerances,
) -> OgeomResult<Option<usize>> {
    let copy = ogeom_algo::copied(model, shape)?;
    let mut image: FastMap<TShapeId, Shape> = FastMap::default();
    for given in copy.history.inputs() {
        if let [made] = copy.history.modified(&given) {
            image.insert(given.node(), Shape::of(made.node()));
        }
    }
    for &k in which {
        for (sub, reach) in &twins[k].widen {
            let Some(made) = image.get(&sub.node()) else {
                return Ok(None);
            };
            model.widen(made, Tolerance::new(*reach)?)?;
        }
    }
    let faces = explore_unique(model, &copy.shape, ShapeType::Face)?;
    if faces.is_empty() {
        return Ok(None);
    }
    Ok(Some(ogeom_algo::sew(model, &faces, tol)?.joined))
}

/// The free edges of `shape`'s shells that are twins across a gap no wider
/// than `cap`, as the module documentation says, with what each widens to.
fn twins(model: &Model, shape: &Shape, cap: f64, tol: Tolerances) -> OgeomResult<Vec<Twins>> {
    let mut faces: Vec<Shape> = Vec::new();
    for shell in explore_unique(model, shape, ShapeType::Shell)? {
        for face in explore_unique(model, &shell, ShapeType::Face)? {
            if !faces.iter().any(|f| f.is_same(&face)) {
                faces.push(face);
            }
        }
    }
    let uses = ogeom_algo::edge_uses(model, shape)?;
    let mut free: Vec<Free> = Vec::new();
    for (k, face) in faces.iter().enumerate() {
        for wire in model.children_of(face)? {
            for edge in model.ordered_children_of(&wire)? {
                if uses.get(&edge.node()) != Some(&1) {
                    continue;
                }
                let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
                    continue;
                };
                if data.degenerate {
                    continue;
                }
                let Some((start, end)) = edge_vertices(model, &edge)? else {
                    continue;
                };
                if start.node() == end.node() {
                    continue;
                }
                let ends = [
                    (start.clone(), placed(model, &start)?),
                    (end.clone(), placed(model, &end)?),
                ];
                free.push(Free {
                    edge,
                    face: k,
                    ends,
                });
            }
        }
    }
    if free.len() < 2 {
        return Ok(Vec::new());
    }

    // Each end of `a` with the end of `b` it meets, the nearer way round,
    // and the gap between them.
    let paired = |a: &Free, b: &Free| -> [(usize, usize, f64); 2] {
        let gap = |i: usize, j: usize| a.ends[i].1.distance(b.ends[j].1);
        let opposite = [(0, 1, gap(0, 1)), (1, 0, gap(1, 0))];
        let same = [(0, 0, gap(0, 0)), (1, 1, gap(1, 1))];
        let widest = |p: &[(usize, usize, f64); 2]| p[0].2.max(p[1].2);
        if widest(&opposite) <= widest(&same) {
            opposite
        } else {
            same
        }
    };
    let mut near: Vec<Vec<usize>> = vec![Vec::new(); free.len()];
    for i in 0..free.len() {
        for j in i + 1..free.len() {
            if free[i].face == free[j].face {
                continue;
            }
            let ends = paired(&free[i], &free[j]);
            if ends[0].2.max(ends[1].2) <= cap {
                near[i].push(j);
                near[j].push(i);
            }
        }
    }
    let mut shortest: FastMap<usize, f64> = FastMap::default();
    let mut vertices: Option<Vec<(TShapeId, Point)>> = None;
    let mut found = Vec::new();
    for i in 0..free.len() {
        let [j] = near[i][..] else { continue };
        if j < i || near[j][..] != [i] {
            continue;
        }
        let (a, b) = (&free[i], &free[j]);
        let Some(apart) = edges_apart(model, &a.edge, &b.edge, tol)? else {
            continue;
        };
        let ends = paired(a, b);
        let widest = apart.max(ends[0].2).max(ends[1].2);
        if apart > cap || widest <= tol.confusion() {
            continue;
        }
        let mut feature = f64::INFINITY;
        for face in [a.face, b.face] {
            let length = match shortest.get(&face) {
                Some(length) => *length,
                None => {
                    let length = shortest_edge(model, &faces[face], tol)?;
                    shortest.insert(face, length);
                    length
                }
            };
            feature = feature.min(length);
        }
        if widest >= feature {
            continue;
        }
        let vertices = match &mut vertices {
            Some(v) => v,
            None => vertices.insert(
                explore_unique(model, shape, ShapeType::Vertex)?
                    .iter()
                    .map(|v| Ok((v.node(), placed(model, v)?)))
                    .collect::<OgeomResult<Vec<_>>>()?,
            ),
        };
        let ours: Vec<TShapeId> = a
            .ends
            .iter()
            .chain(&b.ends)
            .map(|(v, _)| v.node())
            .collect();
        let crowded = a.ends.iter().chain(&b.ends).any(|(_, at)| {
            vertices
                .iter()
                .any(|(id, p)| !ours.contains(id) && p.distance(*at) <= widest)
        });
        if crowded {
            continue;
        }
        let mut widen: Vec<(Shape, f64)> = vec![
            (a.edge.clone(), apart.max(tol.confusion())),
            (b.edge.clone(), apart.max(tol.confusion())),
        ];
        for (i, j, gap) in ends {
            if gap > 0.0 {
                widen.push((a.ends[i].0.clone(), gap));
                widen.push((b.ends[j].0.clone(), gap));
            }
        }
        found.push(Twins { widen });
    }
    Ok(found)
}

/// Where `vertex` stands, placed.
fn placed(model: &Model, vertex: &Shape) -> OgeomResult<Point> {
    let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
        ogeom_core::ogeom_bail!(Dangling, "vertex is not in this model");
    };
    Ok(vertex.transform(model.datums())?.apply(data.point))
}

/// The length of `face`'s shortest edge with length.
fn shortest_edge(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<f64> {
    let mut shortest = f64::INFINITY;
    for edge in explore_unique(model, face, ShapeType::Edge)? {
        if model
            .node(&edge)
            .and_then(|n| n.data().as_edge())
            .is_some_and(|d| d.degenerate)
        {
            continue;
        }
        shortest = shortest.min(linear_properties(model, &edge, Deflection::default(), tol)?.mass);
    }
    Ok(shortest)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_algo::{make_face, make_polygon};
    use ogeom_geom::{PlaneSurface, SurfaceGeometry};
    use ogeom_math::{Frame, Plane};

    const T: Tolerances = Tolerances::millimetres();

    /// Every edge's and vertex's tolerance in `shape`, by node.
    fn tolerances(model: &Model, shape: &Shape) -> Vec<(TShapeId, f64)> {
        let mut out = Vec::new();
        for kind in [ShapeType::Edge, ShapeType::Vertex] {
            for sub in explore_unique(model, shape, kind).unwrap() {
                let t = model.tolerance_of(&sub).unwrap().map_or(0.0, |t| t.get());
                out.push((sub.node(), t));
            }
        }
        out
    }

    /// Two 10 mm squares side by side, 5 mm apart, in one shell: the
    /// first's right edge and the second's left edge, handed over as twins
    /// stated at a micron, do not sew. Nothing is widened, and every
    /// tolerance in the shell is as it was.
    #[test]
    fn twins_sewing_does_not_join_leave_the_shape_as_it_was() {
        let mut model = Model::new();
        let p = |x: f64, y: f64| Point::new(x, y, 0.0);
        let mut faces = Vec::new();
        for x in [0.0, 15.0] {
            let quad = [p(x, 0.0), p(x + 10.0, 0.0), p(x + 10.0, 10.0), p(x, 10.0)];
            let ring = make_polygon(&mut model, &quad, true, T).unwrap().shape;
            let plane = SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(Frame::WORLD)));
            faces.push(make_face(&mut model, plane, &[ring], T).unwrap().shape);
        }
        let shell = model.add_shell(&faces).unwrap();
        let side = |model: &Model, face: &Shape, x: f64| {
            explore_unique(model, face, ShapeType::Edge)
                .unwrap()
                .into_iter()
                .find(|e| {
                    edge_vertices(model, e).unwrap().is_some_and(|(a, b)| {
                        [a, b]
                            .iter()
                            .all(|v| (placed(model, v).unwrap().x - x).abs() < 1e-9)
                    })
                })
                .unwrap()
        };
        let (right, left) = (side(&model, &faces[0], 10.0), side(&model, &faces[1], 15.0));
        let mut widen = vec![(right.clone(), 1e-3), (left.clone(), 1e-3)];
        for edge in [&right, &left] {
            let (a, b) = edge_vertices(&model, edge).unwrap().unwrap();
            widen.extend([(a, 1e-3), (b, 1e-3)]);
        }
        let before = tolerances(&model, &shell);
        let stated = state_joining(&mut model, &shell, &[Twins { widen }], T).unwrap();
        assert_eq!(stated, 0);
        assert_eq!(tolerances(&model, &shell), before);
    }
}
