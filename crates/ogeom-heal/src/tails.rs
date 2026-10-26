//! Tails stripped from wires: a boundary that runs out along an edge and
//! straight back along the same curve.
//!
//! Two consecutive edges of a wire, the first from one vertex to another
//! and the second from that one back, that lie along each other within the
//! tolerance either states, bound nothing between them: the region the
//! wire encloses is the same with them as without. Kept, the tip is a
//! vertex of the boundary that bounds no corner, every boolean through the
//! face meets a strand there, and a mesher has to strip it before it can
//! triangulate. A wire made of such pairs alone (a slit) bounds nothing at
//! all and, as a hole, goes; the outer wire of a face is never stripped to
//! nothing.
//!
//! A pair is a tail only where all of these hold, so that nothing a
//! neighbour holds or the chart tells apart is taken for one:
//!
//! - neither edge is degenerate, and the pair's ends are two vertices;
//! - no other face of the shape walks either edge, nor this face again,
//!   except as the pair itself (one edge walked out and back);
//! - neither edge is a seam on the face, and in the face's chart the walk
//!   back ends where the walk out began, and the tip is one point, with no
//!   whole period between them along a direction the surface closes in:
//!   two meridians a period apart on a sphere are one curve in space and
//!   bound the whole sphere between them in the chart;
//! - the two curves stand no farther apart than the wider of the two
//!   edges' tolerances anywhere along them.

use ogeom_algo::edge_vertices;
use ogeom_core::{FastMap, FastSet, OgeomResult, Tolerances};
use ogeom_geom::{Curve2d as _, Surface as _, SurfaceGeometry};
use ogeom_math::Point2;
use ogeom_topo::{
    EdgeRepr, Model, Orientation, Shape, ShapeType, SurfaceId, TShapeId, explore_unique,
};

use crate::Reshape;

/// Stage, for every wire of `shape` with a tail (see the module
/// documentation), the wire without it, and for every hole made of tails
/// alone, its removal. Returns how many pairs of edges are stripped.
pub(crate) fn stage_tails(
    model: &mut Model,
    shape: &Shape,
    reshape: &mut Reshape,
    tol: Tolerances,
) -> OgeomResult<usize> {
    let faces = explore_unique(model, shape, ShapeType::Face)?;
    let mut uses: FastMap<TShapeId, usize> = FastMap::default();
    for face in &faces {
        for wire in model.children_of(face)? {
            for edge in model.children_of(&wire)? {
                *uses.entry(edge.node()).or_default() += 1;
            }
        }
    }
    let mut seen: FastSet<TShapeId> = FastSet::default();
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
            if !seen.insert(wire.node()) {
                continue;
            }
            let walked = model.ordered_children_of(&wire.oriented(Orientation::Forward))?;
            let tail = |i: usize, j: usize| -> OgeomResult<bool> {
                is_tail(
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
    for (wire, kept) in staged {
        let node = Shape::of(wire.node());
        match kept {
            None => reshape.remove(&node),
            Some(kept) => {
                let local = model.ordered_children_of(&node)?;
                let edges: Vec<Shape> = kept.into_iter().map(|i| local[i].clone()).collect();
                let rebuilt = model.add_wire(&edges)?;
                reshape.replace(&node, rebuilt);
            }
        }
    }
    Ok(stripped)
}

/// Whether `out`, then `back`, walked in that order by one wire of the
/// face on `surface`, are a tail. See the module documentation.
fn is_tail(
    model: &Model,
    (surface_id, surface): (SurfaceId, &SurfaceGeometry),
    out: &Shape,
    back: &Shape,
    uses: &FastMap<TShapeId, usize>,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let data = |edge: &Shape| model.node(edge).and_then(|n| n.data().as_edge());
    let (Some(out_data), Some(back_data)) = (data(out), data(back)) else {
        return Ok(false);
    };
    if out_data.degenerate || back_data.degenerate {
        return Ok(false);
    }
    let used = |edge: &Shape| uses.get(&edge.node()).copied().unwrap_or(0);
    let alone = if out.node() == back.node() {
        used(out) == 2 && out.orientation() != back.orientation()
    } else {
        used(out) == 1 && used(back) == 1
    };
    if !alone {
        return Ok(false);
    }
    let (Some((a, b)), Some((c, d))) = (edge_vertices(model, out)?, edge_vertices(model, back)?)
    else {
        return Ok(false);
    };
    if a.node() == b.node() || c.node() != b.node() || d.node() != a.node() {
        return Ok(false);
    }
    let (Some(out_chart), Some(back_chart)) = (
        chart_ends(model, surface_id, out, tol)?,
        chart_ends(model, surface_id, back, tol)?,
    ) else {
        return Ok(false);
    };
    if !within_a_period(surface, out_chart.0, back_chart.1, tol)
        || !within_a_period(surface, out_chart.1, back_chart.0, tol)
    {
        return Ok(false);
    }
    if out.node() == back.node() {
        return Ok(true);
    }
    let reach = out_data.tolerance.get().max(back_data.tolerance.get());
    Ok(ogeom_algo::edges_apart(model, out, back, tol)?.is_some_and(|apart| apart <= reach))
}

/// Where `edge`'s one plain pcurve on `surface` begins and ends as its
/// occurrence walks it. `None` for a seam, or an edge with no pcurve there.
fn chart_ends(
    model: &Model,
    surface: SurfaceId,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(Point2, Point2)>> {
    let Some(data) = model.node(edge).and_then(|n| n.data().as_edge()) else {
        return Ok(None);
    };
    let Some(EdgeRepr::PCurve {
        curve,
        range,
        location,
        ..
    }) = data.pcurve_for(surface, edge.location())
    else {
        return Ok(None);
    };
    // Two plain pcurves on one surface at one placement are a seam stored
    // as two columns, walked once each way.
    let plain = data
        .representations
        .iter()
        .filter(|r| {
            matches!(r, EdgeRepr::PCurve { .. })
                && r.surface() == Some(surface)
                && r.location() == Some(location)
        })
        .count();
    if plain > 1 {
        return Ok(None);
    }
    let Some(pcurve) = model.geometry().pcurve(*curve) else {
        return Ok(None);
    };
    let (start, end) = (
        pcurve.point_at(range.0, tol)?,
        pcurve.point_at(range.1, tol)?,
    );
    Ok(Some(if edge.orientation() == Orientation::Reversed {
        (end, start)
    } else {
        (start, end)
    }))
}

/// Whether `p` and `q` stand less than half a period apart along every
/// direction `surface` closes in: the same copy of the chart.
fn within_a_period(surface: &SurfaceGeometry, p: Point2, q: Point2, tol: Tolerances) -> bool {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let closes = [
        surface.is_periodic_u() || surface.is_closed_u(tol),
        surface.is_periodic_v() || surface.is_closed_v(tol),
    ];
    [
        (closes[0], u1 - u0, (p.x - q.x).abs()),
        (closes[1], v1 - v0, (p.y - q.y).abs()),
    ]
    .into_iter()
    .all(|(closes, period, apart)| !closes || !period.is_finite() || 2.0 * apart < period)
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
        let uses = walked.iter().map(|e| (e.node(), 1)).collect();
        let on = (surface, model.geometry().surface(surface).unwrap());
        assert!(!is_tail(&model, on, &walked[0], &walked[1], &uses, T).unwrap());
    }
}
