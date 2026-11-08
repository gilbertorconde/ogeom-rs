//! Where three bevels meet at a vertex: a point or a facet.

use crate::Chamfer;
use crate::chamfer::{Legs, legs};
use crate::support::planar_face;
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape};

/// How a chamfer request closes a vertex where three of its bevels meet.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChamferCorner {
    /// The bevels mitre into one point, as
    /// [`chamfer_edges_with`](crate::chamfer_edges_with) leaves them.
    Point,
    /// Each bevel stops at the plane through the three points where two
    /// bevels' contact lines cross on the face they share, and the corner
    /// is the triangle there.
    Facet,
}

/// [`chamfer_edges_with`](crate::chamfer_edges_with), the vertices where
/// three of its bevels meet closed as `corner` says.
///
/// [`ChamferCorner::Point`] is exactly
/// [`chamfer_edges_with`](crate::chamfer_edges_with). For
/// [`ChamferCorner::Facet`] each such vertex is a convex corner of three
/// planes, each pair of its bevels sharing one face there: on that face the
/// two bevels' contact lines cross at a point, and the three points span
/// the facet's plane. The mitred bevels meet at the point their three planes
/// share, and the tetrahedron between it and the three points, whose sides
/// lie on the bevels, is cut away: each bevel ends on the facet, and a box's
/// corner of three bevels `d` gains the triangle cutting each edge at `2d`
/// from the corner, its point's tetrahedron `d³ / 12` taken off.
///
/// # Errors
///
/// As [`chamfer_edges_with`](crate::chamfer_edges_with), and for
/// [`ChamferCorner::Facet`]
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// a vertex holds more than three of the request's bevels, or three of
/// which one is not a straight convex edge between planes, two share no
/// face, or the facet's points do not stand on the edges' own lengths.
pub fn chamfer_edges_cornered(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    corner: ChamferCorner,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if corner == ChamferCorner::Point {
        return crate::chamfer_edges_with(model, solid, specs, tol);
    }
    let tetrahedra = facets(model, solid, specs, tol)?;
    let mut built = crate::chamfer_edges_with(model, solid, specs, tol)?;
    for corners in tetrahedra {
        let faces = tetrahedron(model, corners, tol)?;
        let step = crate::support::apply_wedge(model, &built.shape, None, &faces, false, tol)?;
        built = Built {
            shape: step.shape.clone(),
            history: built.history.then(&step.history),
        };
    }
    Ok(built)
}

/// The tetrahedron each vertex of three bevels sheds: the bevels' shared
/// point, then the three facet points.
fn facets(
    model: &Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Vec<[Point; 4]>> {
    let conf = tol.confusion();
    let mut ends: Vec<[ogeom_topo::TShapeId; 2]> = Vec::with_capacity(specs.len());
    let mut points: Vec<(ogeom_topo::TShapeId, Point)> = Vec::new();
    for (edge, _) in specs {
        let Some((a, b)) = ogeom_algo::edge_vertices(model, edge)? else {
            ogeom_bail!(Construction, "a chamfered edge has no vertices");
        };
        for v in [&a, &b] {
            let Some(data) = model.node(v).and_then(|n| n.data().as_vertex()) else {
                ogeom_bail!(Dangling, "a vertex is not in this model");
            };
            points.push((v.node(), v.transform(model.datums())?.apply(data.point)));
        }
        ends.push([a.node(), b.node()]);
    }
    let mut vertices: Vec<ogeom_topo::TShapeId> = Vec::new();
    for pair in &ends {
        for v in pair {
            if !vertices.contains(v) && ends.iter().filter(|e| e.contains(v)).count() >= 3 {
                vertices.push(*v);
            }
        }
    }
    let mut out = Vec::with_capacity(vertices.len());
    for v in vertices {
        let at: Vec<usize> = (0..specs.len()).filter(|&i| ends[i].contains(&v)).collect();
        let [i, j, k] = at.as_slice() else {
            ogeom_bail!(
                Construction,
                "a facet closes a corner of three bevels, not {}",
                at.len()
            );
        };
        // Each bevel's seat, its legs' distances, and its edge leaving the
        // vertex.
        let mut bevels = Vec::with_capacity(3);
        for &n in [i, j, k] {
            let (edge, spec) = &specs[n];
            let Legs::Straight(seat, distances) = legs(model, solid, edge, spec, tol)? else {
                ogeom_bail!(
                    Construction,
                    "a facet closes a corner of bevels on straight edges between planes"
                );
            };
            if !seat.convex {
                ogeom_bail!(Construction, "a facet closes a convex corner");
            }
            // The seat's end standing at the vertex.
            let Some(&(_, at)) = points.iter().find(|(id, _)| *id == v) else {
                ogeom_bail!(Dangling, "a vertex is not in this model");
            };
            let at_start = seat.start.distance(at) <= seat.end.distance(at);
            let (vertex, away) = if at_start {
                (seat.start, seat.along)
            } else {
                (seat.end, -seat.along)
            };
            bevels.push((seat, distances, vertex, away));
        }
        // A bevel's plane: through its contact lines, normal to them both.
        let plane = |n: usize| -> OgeomResult<(Point, Vector)> {
            let (seat, d, _, _) = &bevels[n];
            let a = seat.start + seat.leg(0, tol)? * d[0];
            let b = seat.start + seat.leg(1, tol)? * d[1];
            Ok((a, seat.along.cross(b - a)))
        };
        let planes = [plane(0)?, plane(1)?, plane(2)?];
        let Some(point) = meet(&planes) else {
            ogeom_bail!(Construction, "the corner's three bevels share no point");
        };
        // On the face two bevels share, where their contact lines cross.
        let mut feet = Vec::with_capacity(3);
        for (m, n) in [(0, 1), (1, 2), (2, 0)] {
            let (sm, dm, vm, um) = &bevels[m];
            let (sn, dn, vn, un) = &bevels[n];
            let shared: Vec<(usize, usize)> = (0..2)
                .flat_map(|a| (0..2).map(move |b| (a, b)))
                .filter(|&(a, b)| sm.faces[a].node() == sn.faces[b].node())
                .collect();
            let [(a, b)] = shared.as_slice() else {
                ogeom_bail!(
                    Construction,
                    "two bevels of a faceted corner share no single face"
                );
            };
            let pm = *vm + sm.leg(*a, tol)? * dm[*a];
            let pn = *vn + sn.leg(*b, tol)? * dn[*b];
            let Some((s, t, gap)) = crossing(pm, *um, pn, *un) else {
                ogeom_bail!(
                    Construction,
                    "two bevels of a faceted corner run parallel on their face"
                );
            };
            // The crossing stands on both edges' own lengths, and the lines
            // cross on the face.
            let reach = |seat: &crate::support::Seat| seat.start.distance(seat.end);
            if gap > conf * 10.0
                || s <= conf
                || t <= conf
                || s >= reach(sm) - conf
                || t >= reach(sn) - conf
            {
                ogeom_bail!(
                    Construction,
                    "the facet's point on a face stands off the edges' lengths"
                );
            }
            feet.push(pm + *um * s);
        }
        // The point is on the vertex's side of the facet.
        let normal = (feet[1] - feet[0]).cross(feet[2] - feet[0]);
        let side = |p: Point| normal.dot(p - feet[0]);
        let vertex = bevels[0].2;
        let depth = side(point) / normal.magnitude();
        if depth.abs() <= conf || depth.signum() != side(vertex).signum() {
            ogeom_bail!(
                Construction,
                "the bevels' shared point stands on the facet's far side"
            );
        }
        out.push([point, feet[0], feet[1], feet[2]]);
    }
    Ok(out)
}

/// The point three planes, each a point and a normal, share.
fn meet(planes: &[(Point, Vector); 3]) -> Option<Point> {
    let [(p0, n0), (p1, n1), (p2, n2)] = *planes;
    let det = n0.dot(n1.cross(n2));
    let scale = n0.magnitude() * n1.magnitude() * n2.magnitude();
    if det.abs() <= 1e-12 * scale {
        return None;
    }
    let (d0, d1, d2) = (
        n0.dot(p0 - Point::ORIGIN),
        n1.dot(p1 - Point::ORIGIN),
        n2.dot(p2 - Point::ORIGIN),
    );
    let v = (n1.cross(n2) * d0 + n2.cross(n0) * d1 + n0.cross(n1) * d2) / det;
    Some(Point::ORIGIN + v)
}

/// Where the lines `p + s u` and `q + t v` come closest: `s`, `t` and the
/// gap there; `None` for parallel lines.
fn crossing(p: Point, u: Vector, q: Point, v: Vector) -> Option<(f64, f64, f64)> {
    let w = p - q;
    let (a, b, c) = (u.dot(u), u.dot(v), v.dot(v));
    let (d, e) = (u.dot(w), v.dot(w));
    let den = a.mul_add(c, -b * b);
    if den.abs() <= 1e-12 * a * c {
        return None;
    }
    let s = b.mul_add(e, -c * d) / den;
    let t = a.mul_add(e, -b * d) / den;
    Some((s, t, (p + u * s).distance(q + v * t)))
}

/// The four faces of the tetrahedron on `corners`, each facing out.
fn tetrahedron(model: &mut Model, corners: [Point; 4], tol: Tolerances) -> OgeomResult<Vec<Shape>> {
    let centroid = Point::ORIGIN
        + (corners
            .iter()
            .fold(Vector::ZERO, |sum, p| sum + (*p - Point::ORIGIN))
            / 4.0);
    let mut faces = Vec::with_capacity(4);
    for skip in 0..4 {
        let tri: Vec<Point> = (0..4).filter(|&n| n != skip).map(|n| corners[n]).collect();
        let mut normal = (tri[1] - tri[0]).cross(tri[2] - tri[0]);
        if normal.dot(tri[0] - centroid) < 0.0 {
            normal = -normal;
        }
        let normal = normal / normal.magnitude();
        faces.push(planar_face(model, &tri, normal, tol)?);
    }
    Ok(faces)
}
