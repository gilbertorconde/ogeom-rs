//! Tails and notches stripped from wires: a boundary that runs out along
//! an edge and straight back along the same path
//! ([`out_and_back`](ogeom_algo::out_and_back) says which pairs of edges
//! do, and why each bounds nothing).
//!
//! A tail runs out and back whole: the pair goes. Pairs collapse as they
//! close, so a tail with a tail on its tip goes as a whole. A wire made of
//! such pairs alone (a slit) bounds nothing at all and, as a hole, goes;
//! the outer wire of a face is never stripped to nothing.
//!
//! A notch runs back short of where it ran out from, or past it: the
//! longer edge is cut where the shorter one's far end stands on it, the
//! piece out to the tip goes with the shorter edge, and the piece kept
//! ends at that far vertex. The vertex stays where it stands and widens to
//! reach where the piece kept ends, no farther than the tolerance the
//! notch was found within: the same slop as a curve ending away from its
//! vertex, stated the same way. Every face the cut changes (the face, and
//! those sharing that vertex) must check clean after it, or it is not made.

use ogeom_algo::{edge_vertices, linear_properties};
use ogeom_core::{FastMap, FastSet, OgeomResult, Tolerances};
use ogeom_mesh::Deflection;
use ogeom_topo::{
    EdgeRepr, Model, Orientation, SameKey, Shape, ShapeType, TShapeId, VertexData, explore_unique,
};

use crate::Reshape;
use crate::crossing::rebound;
use crate::divide::{split_edge, vertex_point};

/// Stage, for every wire of `shape` with a tail (see the module
/// documentation), the wire without it, and for every hole made of tails
/// alone, its removal. Returns how many pairs of edges are stripped.
///
/// Each placed occurrence of a wire is judged and stripped on its own: a
/// copy whose tail bounds a neighbour there keeps it, and copies stripped
/// alike share the stripped wire.
pub(crate) fn stage_tails(
    model: &mut Model,
    shape: &Shape,
    reshape: &mut Reshape,
    tol: Tolerances,
) -> OgeomResult<usize> {
    let faces = explore_unique(model, shape, ShapeType::Face)?;
    let uses = ogeom_algo::placed_edge_uses(model, shape)?;
    let mut seen: FastSet<SameKey> = FastSet::default();
    let mut staged: Vec<(Shape, Option<Vec<usize>>)> = Vec::new();
    let mut stripped = 0;
    for face in &faces {
        let Some(surface_id) = model
            .node(face)
            .and_then(|n| n.data().as_face())
            .map(|d| d.surface)
        else {
            continue;
        };
        let Some(surface) = model.geometry().surface(surface_id) else {
            continue;
        };
        let stored = face.oriented(Orientation::Forward);
        for (k, wire) in model.ordered_children_of(&stored)?.into_iter().enumerate() {
            if !seen.insert(SameKey(wire.clone())) {
                continue;
            }
            let walked = model.ordered_children_of(&wire.oriented(Orientation::Forward))?;
            let tail = |i: usize, j: usize| -> OgeomResult<bool> {
                ogeom_algo::is_tail(
                    model,
                    (surface_id, surface),
                    &walked[i],
                    &walked[j],
                    &uses,
                    tol,
                )
            };
            // Pairs struck as they close, so a tail with a tail on its tip
            // goes as a whole.
            let mut kept: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
            for i in 0..walked.len() {
                if let Some(&top) = kept.back()
                    && tail(top, i)?
                {
                    kept.pop_back();
                } else {
                    kept.push_back(i);
                }
            }
            // Round the wire's start: its last edge runs into its first.
            while kept.len() >= 2 {
                let (Some(&last), Some(&first)) = (kept.back(), kept.front()) else {
                    break;
                };
                if !tail(last, first)? {
                    break;
                }
                kept.pop_back();
                kept.pop_front();
            }
            if kept.len() == walked.len() || (kept.is_empty() && k == 0) {
                continue;
            }
            stripped += (walked.len() - kept.len()) / 2;
            staged.push((
                wire.clone(),
                (!kept.is_empty()).then(|| kept.into_iter().collect()),
            ));
        }
    }
    let mut made: FastMap<(TShapeId, Vec<usize>), Shape> = FastMap::default();
    for (wire, kept) in staged {
        match kept {
            None => reshape.remove_at(&wire),
            Some(kept) => {
                let key = (wire.node(), kept);
                let rebuilt = match made.get(&key) {
                    Some(rebuilt) => rebuilt.clone(),
                    None => {
                        let local = model.ordered_children_of(&Shape::of(wire.node()))?;
                        let edges: Vec<Shape> = key.1.iter().map(|&i| local[i].clone()).collect();
                        let rebuilt = model.add_wire(&edges)?;
                        made.insert(key, rebuilt.clone());
                        rebuilt
                    }
                };
                reshape.replace_at(&wire, placed_at(rebuilt, &wire));
            }
        }
    }
    Ok(stripped)
}

/// Stage, for every notch of `shape` (see the module documentation), the
/// cut that takes it off. A face is cut at one notch per pass, and a face
/// one cut changes is not changed by another. Returns how many notches are
/// cut.
pub(crate) fn stage_notches(
    model: &mut Model,
    shape: &Shape,
    reshape: &mut Reshape,
    tol: Tolerances,
) -> OgeomResult<usize> {
    let found = ogeom_algo::out_and_back(model, shape, tol)?;
    if found.iter().all(|pair| pair.notch.is_none()) {
        return Ok(0);
    }
    let faces = explore_unique(model, shape, ShapeType::Face)?;
    let mut holders: FastMap<SameKey, Vec<usize>> = FastMap::default();
    for (k, face) in faces.iter().enumerate() {
        for vertex in explore_unique(model, face, ShapeType::Vertex)? {
            holders.entry(SameKey(vertex)).or_default().push(k);
        }
    }
    // A cut is made on the nodes, in their own frame, once for every
    // placement it is staged at.
    let mut made: FastMap<(TShapeId, TShapeId), Option<(Shape, Shape)>> = FastMap::default();
    let mut changed: FastSet<usize> = FastSet::default();
    let mut cut = 0;
    for pair in found {
        let Some(notch) = pair.notch else { continue };
        let Some(k) = faces.iter().position(|f| f.is_same(&pair.face)) else {
            continue;
        };
        if changed.contains(&k) {
            continue;
        }
        let stored = pair.face.oriented(Orientation::Forward);
        let wire = model.ordered_children_of(&stored)?[pair.wire].clone();
        let walked = model.ordered_children_of(&wire.oriented(Orientation::Forward))?;
        let (out, back) = (&walked[pair.edges[0]], &walked[pair.edges[1]]);
        let (Some((a, _)), Some((_, d))) =
            (edge_vertices(model, out)?, edge_vertices(model, back)?)
        else {
            continue;
        };
        let (longer, shorter, far) = if notch.longer == 0 {
            (out, back, d)
        } else {
            (back, out, a)
        };
        let Some(EdgeRepr::Curve3d { range, .. }) = model
            .node(longer)
            .and_then(|n| n.data().as_edge())
            .and_then(|d| d.curve3d().cloned())
        else {
            continue;
        };
        // The edges and the vertex at the face's own placement, rigid, so
        // the cut made in their own frame is the cut where they stand.
        let at = pair.face.location();
        if [longer.location(), shorter.location(), far.location()]
            .iter()
            .any(|l| *l != at)
            || (!at.is_identity() && at.composed(model.datums())?.scale_factor().abs() != 1.0)
        {
            continue;
        }
        let key = (longer.node(), far.node());
        let cut_made = match made.get(&key) {
            Some(cut_made) => cut_made.clone(),
            None => {
                // The tip is where the edge out ends and the edge back begins.
                let tip_high =
                    (longer.orientation() == Orientation::Reversed) != (notch.longer == 0);
                let fraction = (notch.parameter - range.0) / (range.1 - range.0);
                let cut_made = notch_cut(
                    model,
                    &Shape::of(longer.node()),
                    &Shape::of(far.node()),
                    fraction,
                    tip_high,
                    tol,
                )?;
                made.insert(key, cut_made.clone());
                cut_made
            }
        };
        let Some((new_longer, widened)) = cut_made else {
            continue;
        };
        let stage = |into: &mut Reshape| {
            into.replace_at(longer, placed_at(new_longer.clone(), longer));
            into.remove_at(shorter);
            into.replace_at(&far, placed_at(widened.clone(), &far));
        };
        let mut local = Reshape::new();
        stage(&mut local);
        let mut hit: Vec<usize> = vec![k];
        hit.extend(
            holders
                .get(&SameKey(far.clone()))
                .into_iter()
                .flatten()
                .copied(),
        );
        hit.sort_unstable();
        hit.dedup();
        if hit.iter().any(|h| changed.contains(h)) {
            continue;
        }
        let mut clean = true;
        for &h in &hit {
            let after = local.apply(model, &faces[h])?.shape;
            clean &= ogeom_algo::check(model, &after, tol)?.is_valid();
        }
        if !clean {
            continue;
        }
        stage(reshape);
        changed.extend(hit);
        cut += 1;
    }
    Ok(cut)
}

/// The cut off a notch on `longer` at `fraction` of its range, the tip at
/// its high end where `tip_high`, its piece kept ending at `far`, all in
/// their own frame: the piece kept, rebounded on `far` widened to reach
/// where it ends, and that vertex. `None` where the piece kept is no
/// longer than that reach.
fn notch_cut(
    model: &mut Model,
    longer: &Shape,
    far: &Shape,
    fraction: f64,
    tip_high: bool,
    tol: Tolerances,
) -> OgeomResult<Option<(Shape, Shape)>> {
    let (pieces, made) = split_edge(model, longer, &[fraction], None, tol)?;
    let ([first, second], [middle]) = (&pieces[..], &made[..]) else {
        return Ok(None);
    };
    let kept = if tip_high { first } else { second };
    let stood = vertex_point(model, far)?;
    let ends = vertex_point(model, middle)?;
    let stated = model.tolerance_of(far)?.map_or(0.0, |t| t.get());
    let reach = stated.max(stood.distance(ends)).max(tol.confusion());
    let widened = model.add_vertex(VertexData::with_tolerance(stood, reach)?);
    let new_longer = rebound(model, kept, middle, &widened)?;
    if linear_properties(model, &new_longer, Deflection::default(), tol)?.mass <= reach {
        return Ok(None);
    }
    Ok(Some((new_longer, widened)))
}

/// `made`, a node in its own frame, placed where `occurrence` stands.
fn placed_at(made: Shape, occurrence: &Shape) -> Shape {
    if occurrence.location().is_identity() {
        made
    } else {
        made.located(occurrence.location().clone())
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use ogeom_algo::{
        attach_pcurve, make_edge_between, make_face_on, make_face_with_pcurves, make_vertex,
        make_wire,
    };
    use ogeom_geom::{CircleCurve, Line2d, LineCurve, PlaneSurface, SphereSurface};
    use ogeom_math::{Circle, Direction, Frame, Plane, Point, Sphere};
    use ogeom_topo::Location;

    const T: Tolerances = Tolerances::millimetres();

    fn line(model: &mut Model, from: &Shape, to: &Shape) -> Shape {
        let at = |v: &Shape| model.node(v).unwrap().data().as_vertex().unwrap().point;
        let (a, b) = (at(from), at(to));
        let curve = LineCurve::segment(a, b, T).unwrap();
        make_edge_between(model, curve.into(), (0.0, a.distance(b)), from, to, T)
            .unwrap()
            .shape
    }

    /// The corners of a 4 mm square on the XY plane, and a vertex 1 mm
    /// inside it from the second.
    fn square(model: &mut Model) -> ([Shape; 4], Shape) {
        let corners = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
            .map(|(x, y)| make_vertex(model, Point::new(x, y, 0.0)).shape);
        let tip = make_vertex(model, Point::new(3.0, 1.0, 0.0)).shape;
        (corners, tip)
    }

    fn plane_face(model: &mut Model, rings: &[Vec<Shape>]) -> Shape {
        let plane = PlaneSurface::new(Plane::new(Frame::WORLD)).into();
        make_face_with_pcurves(model, plane, rings, T)
            .unwrap()
            .shape
    }

    fn rim(model: &mut Model, corners: &[Shape; 4]) -> Vec<Shape> {
        (0..4)
            .map(|i| line(model, &corners[i], &corners[(i + 1) % 4]))
            .collect()
    }

    /// Edges per wire of `face` after the tails are stripped, and how many
    /// pairs went.
    fn stripped(model: &mut Model, face: &Shape) -> (Vec<usize>, usize) {
        let mut reshape = Reshape::new();
        let count = stage_tails(model, face, &mut reshape, T).unwrap();
        let face = if reshape.is_empty() {
            face.clone()
        } else {
            reshape.apply(model, face).unwrap().shape
        };
        let wires = model.ordered_children_of(&face).unwrap();
        let lengths = wires
            .iter()
            .map(|w| model.ordered_children_of(w).unwrap().len())
            .collect();
        (lengths, count)
    }

    /// A square whose boundary runs out from a corner along one edge and
    /// back along the same edge: the square, with its four sides.
    #[test]
    fn one_edge_walked_out_and_back_is_stripped() {
        let mut model = Model::new();
        let (corners, tip) = square(&mut model);
        let [a, b, c, d] = rim(&mut model, &corners).try_into().unwrap();
        let spur = line(&mut model, &corners[1], &tip);
        let face = plane_face(
            &mut model,
            &[vec![a, spur.clone(), spur.reversed(), b, c, d]],
        );
        assert_eq!(stripped(&mut model, &face), (vec![4], 1));
    }

    /// The square with a tail out to (3, 1) and back on a second edge, and
    /// a copy of it placed 10 above beside a triangle standing up from the
    /// tail's edge out, which the triangle walks too. The copy below goes
    /// without its tail; the copy above keeps it, still shared with the
    /// triangle, on the nodes it had.
    #[test]
    fn a_copy_whose_tail_another_face_walks_keeps_it() {
        use ogeom_math::{Transform, Vector};
        let mut model = Model::new();
        let (corners, tip) = square(&mut model);
        let [a, b, c, d] = rim(&mut model, &corners).try_into().unwrap();
        let out = line(&mut model, &corners[1], &tip);
        let back = line(&mut model, &tip, &corners[1]);
        let face = plane_face(&mut model, &[vec![a, out.clone(), back, b, c, d]]);
        let apex = make_vertex(&mut model, Point::new(3.5, 0.5, 2.0)).shape;
        let rise = line(&mut model, &tip, &apex);
        let fall = line(&mut model, &apex, &corners[1]);
        let upright = Frame::new(
            Point::new(4.0, 0.0, 0.0),
            Direction::new(Vector::new(1.0, 1.0, 0.0), T).unwrap(),
            Direction::new(Vector::new(-1.0, 1.0, 0.0), T).unwrap(),
            T,
        )
        .unwrap();
        let triangle = make_face_with_pcurves(
            &mut model,
            PlaneSurface::new(Plane::new(upright)).into(),
            &[vec![out.clone(), rise, fall]],
            T,
        )
        .unwrap()
        .shape;
        let beside = model
            .add_compound(&[face.clone(), triangle.clone()])
            .unwrap();
        let above = model.placed(&beside, Transform::translation(Vector::new(0.0, 0.0, 10.0)));
        let root = model.add_compound(&[face.clone(), above]).unwrap();

        let mut reshape = Reshape::new();
        assert_eq!(stage_tails(&mut model, &root, &mut reshape, T).unwrap(), 1);
        let built = reshape.apply(&mut model, &root).unwrap().shape;
        let faces = explore_unique(&model, &built, ShapeType::Face).unwrap();
        assert_eq!(faces.len(), 3);
        let walks = |model: &Model, face: &Shape| {
            let wire = model.ordered_children_of(face).unwrap().remove(0);
            model.ordered_children_of(&wire).unwrap()
        };
        let below = faces.iter().find(|f| f.location().is_identity()).unwrap();
        assert_eq!(walks(&model, below).len(), 4);
        assert_ne!(below.node(), face.node());
        let copy = faces
            .iter()
            .find(|f| !f.location().is_identity() && f.node() == face.node())
            .expect("the copy above keeps its nodes");
        let up = faces.iter().find(|f| f.node() == triangle.node()).unwrap();
        let kept = walks(&model, copy);
        assert_eq!(kept.len(), 6);
        assert!(
            walks(&model, up)
                .iter()
                .any(|e| e.node() == out.node() && kept.iter().any(|k| k.is_same(e))),
            "the tail's edge out is still the triangle's, above"
        );
    }

    /// Out along one line and back along a second over the same points:
    /// the two lie along each other, and go.
    #[test]
    fn two_edges_along_one_line_out_and_back_are_stripped() {
        let mut model = Model::new();
        let (corners, tip) = square(&mut model);
        let [a, b, c, d] = rim(&mut model, &corners).try_into().unwrap();
        let out = line(&mut model, &corners[1], &tip);
        let back = line(&mut model, &tip, &corners[1]);
        let face = plane_face(&mut model, &[vec![a, out, back, b, c, d]]);
        assert_eq!(stripped(&mut model, &face), (vec![4], 1));
    }

    /// Out along a line and back along an arc between the same two
    /// vertices: the two bound a half disc between them, and stay.
    #[test]
    fn a_line_out_and_an_arc_back_are_no_tail() {
        let mut model = Model::new();
        let (corners, _) = square(&mut model);
        let [a, b, c, d] = rim(&mut model, &corners).try_into().unwrap();
        let (from, to) = (Point::new(1.0, 1.0, 0.0), Point::new(3.0, 1.0, 0.0));
        let [p, q] = [from, to].map(|x| make_vertex(&mut model, x).shape);
        let out = line(&mut model, &p, &q);
        let circle =
            Circle::new(Frame::WORLD.with_origin(Point::new(2.0, 1.0, 0.0)), 1.0, T).unwrap();
        let arc = make_edge_between(
            &mut model,
            CircleCurve::new(circle).into(),
            (0.0, core::f64::consts::PI),
            &q,
            &p,
            T,
        )
        .unwrap()
        .shape;
        let rim = vec![a, b, c, d];
        let face = plane_face(&mut model, &[rim, vec![out, arc]]);
        assert_eq!(stripped(&mut model, &face), (vec![4, 2], 0));
    }

    /// A hole that is a slit, out along one line and back along another
    /// over the same points, encloses nothing and goes; the outer wire
    /// keeps its four sides.
    #[test]
    fn a_hole_of_tails_alone_goes() {
        let mut model = Model::new();
        let (corners, _) = square(&mut model);
        let rim = rim(&mut model, &corners);
        let [p, q] = [Point::new(1.0, 1.0, 0.0), Point::new(3.0, 2.0, 0.0)]
            .map(|x| make_vertex(&mut model, x).shape);
        let out = line(&mut model, &p, &q);
        let back = line(&mut model, &q, &p);
        let face = plane_face(&mut model, &[rim, vec![out, back]]);
        assert_eq!(stripped(&mut model, &face), (vec![4], 1));
    }

    /// A face whose only wire is a slit is left as it is: stripped, it
    /// would bound nothing.
    #[test]
    fn an_outer_wire_of_tails_alone_stays() {
        let mut model = Model::new();
        let [p, q] = [Point::new(1.0, 1.0, 0.0), Point::new(3.0, 2.0, 0.0)]
            .map(|x| make_vertex(&mut model, x).shape);
        let out = line(&mut model, &p, &q);
        let back = line(&mut model, &q, &p);
        let face = plane_face(&mut model, &[vec![out, back]]);
        assert_eq!(stripped(&mut model, &face), (vec![2], 0));
    }

    /// A sphere bounded by two edges on one meridian, from pole to pole
    /// and back, each with its plain pcurve a period from the other's:
    /// one curve in space, the whole sphere between them in the chart.
    #[test]
    fn two_edges_a_period_apart_in_the_chart_are_no_tail() {
        let mut model = Model::new();
        let pi = core::f64::consts::PI;
        let sphere = Sphere::new(Frame::WORLD, 1.0, T).unwrap();
        let surface = model
            .geometry_mut()
            .add_surface(SphereSurface::new(sphere).into());
        let [south, north] =
            [-1.0, 1.0].map(|z| make_vertex(&mut model, Point::new(0.0, 0.0, z)).shape);
        let frame = Frame::new(Point::ORIGIN, -Direction::Y, Direction::X, T).unwrap();
        let meridian = CircleCurve::new(Circle::new(frame, 1.0, T).unwrap());
        let mut edge = |u: f64| {
            let edge = make_edge_between(
                &mut model,
                meridian.into(),
                (-pi / 2.0, pi / 2.0),
                &south,
                &north,
                T,
            )
            .unwrap()
            .shape;
            let chart = Line2d::segment(
                ogeom_math::Point2::new(u, -pi / 2.0),
                ogeom_math::Point2::new(u, pi / 2.0),
                T,
            )
            .unwrap();
            attach_pcurve(
                &mut model,
                &edge,
                chart.into(),
                surface,
                Location::identity(),
                (0.0, pi),
            )
            .unwrap();
            edge
        };
        let (east, west) = (edge(0.0), edge(2.0 * pi));
        let wire = make_wire(&mut model, &[east, west.reversed()], T)
            .unwrap()
            .shape;
        let face = make_face_on(&mut model, surface, &[wire], T).unwrap().shape;
        // The pair is the face's whole outer wire, which is never stripped
        // to nothing: asked of the pair itself.
        let wire = model.ordered_children_of(&face).unwrap()[0].clone();
        let walked = model.ordered_children_of(&wire).unwrap();
        let uses = walked.iter().map(|e| (SameKey(e.clone()), 1)).collect();
        let on = (surface, model.geometry().surface(surface).unwrap());
        assert!(!ogeom_algo::is_tail(&model, on, &walked[0], &walked[1], &uses, T).unwrap());
    }
}
