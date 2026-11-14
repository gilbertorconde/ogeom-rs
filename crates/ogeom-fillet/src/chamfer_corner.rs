//! Where three bevels meet at a vertex: a point or a facet.

use crate::Chamfer;
use crate::chamfer::{Legs, legs};
use crate::support::planar_face;
use ogeom_algo::Built;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape};

/// How a chamfer request closes a vertex where three or more of its bevels
/// meet.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChamferCorner {
    /// The bevels mitre into one point, by wedges built on the solid as given
    /// and applied together.
    Point,
    /// Each bevel stops at the plane through the points where two bevels'
    /// contact lines cross on the face they share, and the corner is the
    /// polygon there: a triangle at a corner of three bevels. What
    /// [`chamfer_edges_with`](crate::chamfer_edges_with) builds.
    Facet,
}

/// [`chamfer_edges_with`](crate::chamfer_edges_with), the vertices where
/// three or more of its bevels meet closed as `corner` says.
///
/// [`ChamferCorner::Facet`] is exactly
/// [`chamfer_edges_with`](crate::chamfer_edges_with): on a solid bounded by
/// planes the blend engine closes each such vertex, convex, concave or
/// mixed, with the polygon through its faces' setback corners. Elsewhere,
/// and past what the engine's setbacks hold, each such vertex is a convex corner of three
/// planes, each pair of its bevels sharing one face there: on that face the
/// two bevels' contact lines cross at a point, and the three points span
/// the facet's plane. The mitred bevels meet at the point their three planes
/// share, and the tetrahedron between it and the three points, whose sides
/// lie on the bevels, is cut away: each bevel ends on the facet, and a box's
/// corner of three bevels `d` gains the triangle cutting each edge at `2d`
/// from the corner, its point's tetrahedron `d³ / 12` taken off.
///
/// At a vertex of `k > 3` bevels, every face there between two of them in
/// turn round it, the `k` points where consecutive bevels' contact lines
/// cross must lie on one plane, within the confusion distance, for a facet
/// through them all to exist: as at a pyramid's apex beveled alike. The
/// solid cut away is then the convex piece on the vertex's side of that
/// plane and the solid's side of every bevel's plane, whose faces are the
/// facet and the bevels' parts past it. On a solid bounded by planes, a
/// vertex of four bevels or more whose points leave one plane closes
/// instead on the patch through them, its sides the bevels' ends: the
/// bilinear patch at four, the points blended by the regular polygon's
/// Wachspress coordinates at five or more, a face on each side's sector.
///
/// # Errors
///
/// As [`chamfer_edges_with`](crate::chamfer_edges_with), and for
/// [`ChamferCorner::Facet`]
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// a vertex of three or more of the request's bevels holds one that is not
/// a straight convex edge between planes, two in turn round it that share
/// no face, facet points that do not stand on the edges' own lengths, or,
/// past three, facet points off one plane that no patch closes (where the
/// engine hands the vertex back).
pub fn chamfer_edges_cornered(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    corner: ChamferCorner,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if specs.is_empty() {
        ogeom_bail!(Construction, "a chain of no edges bevels nothing");
    }
    let (specs, _) = crate::chamfer::without_seams(model, solid, specs, tol)?;
    if specs.is_empty() {
        return Ok(crate::chamfer::unchanged(solid));
    }
    let specs = &specs[..];
    if corner == ChamferCorner::Point {
        // On a solid bounded by planes the blend engine mitres in one
        // pass; elsewhere the wedges do.
        if let Some(built) = crate::chamfer::engine_chamfer_mitred(model, solid, specs, tol)? {
            return Ok(built);
        }
        return crate::chamfer::chamfer_edges_mitred(model, solid, specs, tol);
    }
    // A solid bounded by planes is chamfered by the blend engine in one
    // pass: every face set back, every vertex its facet.
    if let Some(built) = crate::chamfer::engine_chamfer(model, solid, specs, tol)? {
        return Ok(built);
    }
    let cuts = facets(model, solid, specs, tol)?;
    let mut built = crate::chamfer::chamfer_edges_mitred(model, solid, specs, tol)?;
    for cut in cuts {
        let faces = match cut {
            Cut::Tetrahedron(corners) => tetrahedron(model, corners, tol)?,
            Cut::Polytope(polygons) => {
                let mut faces = Vec::with_capacity(polygons.len());
                for (ring, outward) in &polygons {
                    faces.push(planar_face(model, ring, *outward, tol)?);
                }
                faces
            }
        };
        let step = crate::support::apply_wedge(model, &built.shape, None, &faces, false, tol)?;
        built = Built {
            shape: step.shape.clone(),
            history: built.history.then(&step.history),
        };
    }
    Ok(built)
}

/// What a faceted corner sheds.
enum Cut {
    /// At a corner of three bevels: the bevels' shared point, then the
    /// three facet points.
    Tetrahedron([Point; 4]),
    /// At a corner of more: the convex piece's faces, each a ring of
    /// points with its outward normal.
    Polytope(Vec<(Vec<Point>, Vector)>),
}

/// The piece each vertex of three or more bevels sheds.
fn facets(
    model: &Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Vec<Cut>> {
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
            out.push(Cut::Polytope(polygon_corner(
                model, solid, specs, &at, v, &points, tol,
            )?));
            continue;
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
        out.push(Cut::Tetrahedron([point, feet[0], feet[1], feet[2]]));
    }
    Ok(out)
}

/// A bevel's seat at vertex `v`, its legs' distances, its end at the
/// vertex and its direction away from it.
type Bevel = (crate::support::Seat, [f64; 2], Point, Vector);

/// The faces of the convex piece a vertex of the `at` bevels, more than
/// three, sheds: on the vertex's side of the plane through the points
/// where consecutive bevels' contact lines cross on the face they share,
/// and on the solid's side of every bevel's plane.
fn polygon_corner(
    model: &Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    at: &[usize],
    v: ogeom_topo::TShapeId,
    points: &[(ogeom_topo::TShapeId, Point)],
    tol: Tolerances,
) -> OgeomResult<Vec<(Vec<Point>, Vector)>> {
    let conf = tol.confusion();
    let Some(&(_, vertex)) = points.iter().find(|(id, _)| *id == v) else {
        ogeom_bail!(Dangling, "a vertex is not in this model");
    };
    let mut bevels: Vec<Bevel> = Vec::with_capacity(at.len());
    for &n in at {
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
        let at_start = seat.start.distance(vertex) <= seat.end.distance(vertex);
        let (end, away) = if at_start {
            (seat.start, seat.along)
        } else {
            (seat.end, -seat.along)
        };
        bevels.push((seat, distances, end, away));
    }
    // The face two bevels share, as each one's leg index.
    let shared = |m: usize, n: usize| -> Option<(usize, usize)> {
        let pairs: Vec<(usize, usize)> = (0..2)
            .flat_map(|a| (0..2).map(move |b| (a, b)))
            .filter(|&(a, b)| bevels[m].0.faces[a].node() == bevels[n].0.faces[b].node())
            .collect();
        match pairs.as_slice() {
            [one] => Some(*one),
            _ => None,
        }
    };
    // The bevels in turn round the vertex, each sharing a face with the
    // next.
    let mut ring = vec![0];
    while ring.len() < bevels.len() {
        let last = ring[ring.len() - 1];
        let Some(next) =
            (0..bevels.len()).find(|&n| !ring.contains(&n) && shared(last, n).is_some())
        else {
            ogeom_bail!(
                Construction,
                "two bevels in turn round a faceted corner share no single face"
            );
        };
        ring.push(next);
    }
    if shared(ring[ring.len() - 1], ring[0]).is_none() {
        ogeom_bail!(
            Construction,
            "two bevels in turn round a faceted corner share no single face"
        );
    }
    // Where consecutive bevels' contact lines cross on their face.
    let reach = |seat: &crate::support::Seat| seat.start.distance(seat.end);
    let mut feet = Vec::with_capacity(ring.len());
    for w in 0..ring.len() {
        let (m, n) = (ring[w], ring[(w + 1) % ring.len()]);
        let Some((a, b)) = shared(m, n) else {
            ogeom_bail!(
                Construction,
                "two bevels in turn round a faceted corner share no single face"
            );
        };
        let (sm, dm, vm, um) = &bevels[m];
        let (sn, dn, vn, un) = &bevels[n];
        let pm = *vm + sm.leg(a, tol)? * dm[a];
        let pn = *vn + sn.leg(b, tol)? * dn[b];
        let Some((s, t, gap)) = crossing(pm, *um, pn, *un) else {
            ogeom_bail!(
                Construction,
                "two bevels of a faceted corner run parallel on their face"
            );
        };
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
    // The facet's plane, by Newell's sum, and every point on it.
    let mut normal = Vector::ZERO;
    for w in 0..feet.len() {
        let (p, q) = (
            feet[w] - Point::ORIGIN,
            feet[(w + 1) % feet.len()] - Point::ORIGIN,
        );
        normal += p.cross(q);
    }
    let size = normal.magnitude();
    if size.is_nan() || size <= 0.0 {
        ogeom_bail!(Construction, "the facet's points span no plane");
    }
    let normal = normal / size;
    let centre = Point::ORIGIN
        + feet
            .iter()
            .fold(Vector::ZERO, |sum, p| sum + (*p - Point::ORIGIN))
            / f64::from(u32::try_from(feet.len()).unwrap_or(u32::MAX));
    let off = feet
        .iter()
        .fold(0.0_f64, |worst, p| worst.max(normal.dot(*p - centre).abs()));
    if off > conf {
        ogeom_bail!(
            Construction,
            "the {} facet points of a corner stand {off:.3e} off one plane",
            feet.len()
        );
    }
    // The bounding planes as `(point, normal)`, each normal facing out of
    // the piece: the facet's away from the vertex, each bevel's toward it.
    let toward = if normal.dot(vertex - centre) > 0.0 {
        1.0
    } else {
        -1.0
    };
    let mut planes: Vec<(Point, Vector)> = vec![(centre, normal * -toward)];
    for (seat, d, _, _) in &bevels {
        let a = seat.start + seat.leg(0, tol)? * d[0];
        let b = seat.start + seat.leg(1, tol)? * d[1];
        let mut n = seat.along.cross(b - a);
        if n.dot(vertex - a) < 0.0 {
            n = -n;
        }
        planes.push((a, n / n.magnitude()));
    }
    // The piece's corners: every three planes' shared point inside all.
    let inside = |p: Point| planes.iter().all(|(q, n)| n.dot(p - *q) <= conf * 10.0);
    let mut corners: Vec<Point> = Vec::new();
    for x in 0..planes.len() {
        for y in x + 1..planes.len() {
            for z in y + 1..planes.len() {
                let Some(p) = meet(&[planes[x], planes[y], planes[z]]) else {
                    continue;
                };
                if inside(p) && !corners.iter().any(|c| c.distance(p) <= conf * 10.0) {
                    corners.push(p);
                }
            }
        }
    }
    let middle = Point::ORIGIN
        + corners
            .iter()
            .fold(Vector::ZERO, |sum, p| sum + (*p - Point::ORIGIN))
            / f64::from(u32::try_from(corners.len().max(1)).unwrap_or(u32::MAX));
    // Each plane's face: the corners on it, in turn round its normal.
    let mut faces = Vec::with_capacity(planes.len());
    for (q, n) in &planes {
        let mut on: Vec<Point> = corners
            .iter()
            .copied()
            .filter(|p| n.dot(*p - *q).abs() <= conf * 10.0)
            .collect();
        if on.len() < 3 {
            ogeom_bail!(
                Construction,
                "a faceted corner's piece has a plane touching it at fewer than three corners"
            );
        }
        let mid = Point::ORIGIN
            + on.iter()
                .fold(Vector::ZERO, |sum, p| sum + (*p - Point::ORIGIN))
                / f64::from(u32::try_from(on.len()).unwrap_or(u32::MAX));
        let x = (on[0] - mid) / (on[0] - mid).magnitude();
        let y = n.cross(x);
        on.sort_by(|a, b| {
            let angle = |p: &Point| (*p - mid).dot(y).atan2((*p - mid).dot(x));
            angle(a).total_cmp(&angle(b))
        });
        let outward = if n.dot(mid - middle) >= 0.0 { *n } else { -*n };
        faces.push((on, outward));
    }
    Ok(faces)
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
