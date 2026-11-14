//! P4's fallback (ADR 0001, §4.5, M5a): where the one assembly cannot place
//! a request, the blend is applied by booleans instead. Every bevel takes the
//! material past its plane between its edge and its contact lines, as far as
//! its setback corners: the convex hull of the edge and its four corners.
//! Every faceted vertex takes the corner between the vertex and its facet:
//! the hull of the vertex and its faces' corners there. A convex edge's or
//! vertex's piece is cut from the solid, a concave one's fused into it.
//!
//! The pieces are what the local assembly bounds with its rings, so where
//! both build a request they agree; the fallback takes what the rings cannot
//! (a bevel reaching a face it does not touch, a step across a third face).

use ogeom_algo::{Built, History};
use ogeom_core::{FastMap, OgeomResult, Tolerances};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use super::Id;
use super::classify::{Planar, Seat, faces_round};

/// One convex piece: its corners, and the corners of the face it carries
/// (the bevel's or the facet's), by which its result faces are found.
struct Piece {
    points: Vec<Point>,
    carried: Vec<Point>,
    from: Shape,
    /// The vertices it ends at: an edge's two, a vertex itself.
    ends: Vec<Id>,
    additive: bool,
}

/// The convex hull of `points`, every face a ring of the hull's corners
/// wound counter-clockwise seen from outside, with its outward normal and
/// offset; `None` where the points span no volume.
#[allow(
    clippy::type_complexity,
    reason = "the hull's corners, rings and planes"
)]
fn hull(points: &[Point], eps: f64) -> Option<(Vec<Point>, Vec<Vec<usize>>, Vec<(Vector, f64)>)> {
    let mut at: Vec<Point> = Vec::with_capacity(points.len());
    for p in points {
        if !at.iter().any(|q| q.distance(*p) <= eps) {
            at.push(*p);
        }
    }
    let n = at.len();
    if n < 4 {
        return None;
    }
    let mut planes: Vec<(Vector, f64)> = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                let normal = (at[j] - at[i]).cross(at[k] - at[i]);
                let m = normal.magnitude();
                if m <= eps * eps {
                    continue;
                }
                let mut normal = normal / m;
                let mut off = normal.dot(at[i].to_vector());
                let side: Vec<f64> = at.iter().map(|p| normal.dot(p.to_vector()) - off).collect();
                if side.iter().all(|s| *s <= eps) {
                } else if side.iter().all(|s| *s >= -eps) {
                    normal = -normal;
                    off = -off;
                } else {
                    continue;
                }
                if !planes
                    .iter()
                    .any(|(m, o)| (*m - normal).magnitude() <= 1e-9 && (o - off).abs() <= eps)
                {
                    planes.push((normal, off));
                }
            }
        }
    }
    if planes.len() < 4 {
        return None;
    }
    let mut rings = Vec::with_capacity(planes.len());
    for (normal, off) in &planes {
        let on: Vec<usize> = (0..n)
            .filter(|i| (normal.dot(at[*i].to_vector()) - off).abs() <= eps)
            .collect();
        #[expect(clippy::cast_precision_loss, reason = "a face's few corners")]
        let count = on.len() as f64;
        let centre = on.iter().fold(Vector::ZERO, |s, i| s + at[*i].to_vector()) * (1.0 / count);
        let u = at[on[0]].to_vector() - centre;
        let u = u / u.magnitude().max(f64::MIN_POSITIVE);
        let w = normal.cross(u);
        let mut ring = on;
        ring.sort_by(|a, b| {
            let angle = |i: &usize| {
                let d = at[*i].to_vector() - centre;
                d.dot(w).atan2(d.dot(u))
            };
            angle(a).total_cmp(&angle(b))
        });
        rings.push(ring);
    }
    Some((at, rings, planes))
}

/// Each bevel's piece and each vertex's of two bevels or more, in key
/// order; `None` where a vertex mixes convex and concave bevels.
fn pieces(
    planar: &Planar,
    seats: &FastMap<(Id, Id), Seat>,
    blended: &ogeom_core::FastSet<Id>,
    corners: &FastMap<(usize, Id), Point>,
    steps: &FastMap<(usize, Id), [Point; 2]>,
) -> Option<Vec<Piece>> {
    // The corner of `face` at `v` on `seat`'s own line where the face steps.
    let corner = |seat: &Seat, k: usize, v: Id| -> Point {
        let face = seat.faces[k];
        match steps.get(&(face, v)) {
            Some([p, q]) => {
                let arriving = (k == 0) == (v == seat.ends.1);
                if arriving { *p } else { *q }
            }
            None => corners[&(face, v)],
        }
    };
    let mut pieces: Vec<Piece> = Vec::new();
    let mut keys: Vec<&(Id, Id)> = seats.keys().collect();
    keys.sort_unstable();
    for key in keys {
        let seat = &seats[key];
        if !blended.contains(&seat.edge) {
            continue;
        }
        let (a, b) = seat.ends;
        let carried = vec![
            corner(seat, 0, a),
            corner(seat, 0, b),
            corner(seat, 1, a),
            corner(seat, 1, b),
        ];
        let mut points = carried.clone();
        points.extend([planar.points[&a], planar.points[&b]]);
        pieces.push(Piece {
            points,
            carried,
            from: planar.shapes[&seat.edge].clone(),
            ends: vec![a, b],
            additive: !seat.convex,
        });
    }
    let mut vertices: Vec<Id> = planar.points.keys().copied().collect();
    vertices.sort_unstable();
    for v in vertices {
        let at: Vec<&Seat> = seats
            .values()
            .filter(|s| blended.contains(&s.edge) && (s.ends.0 == v || s.ends.1 == v))
            .collect();
        if at.len() < 2 {
            continue;
        }
        let convex = at.iter().filter(|s| s.convex).count();
        if convex != 0 && convex != at.len() {
            return None;
        }
        let order = faces_round(&planar.faces, v)?;
        let carried: Vec<Point> = order.iter().map(|fi| corners[&(*fi, v)]).collect();
        let mut points = carried.clone();
        points.push(planar.points[&v]);
        pieces.push(Piece {
            points,
            carried,
            from: planar.shapes[&v].clone(),
            ends: vec![v],
            additive: convex == 0,
        });
    }
    Some(pieces)
}

/// The non-local test (ADR 0001, §4.5, M5a): the first piece that reaches
/// or touches a face not round the vertices it ends at, with that face's
/// index, before any surface is built. Its blend would change a face the
/// request does not name, which the rings do not trim.
pub(crate) fn non_local(
    planar: &Planar,
    seats: &FastMap<(Id, Id), Seat>,
    blended: &ogeom_core::FastSet<Id>,
    corners: &FastMap<(usize, Id), Point>,
    steps: &FastMap<(usize, Id), [Point; 2]>,
    tol: Tolerances,
) -> Option<(Shape, usize)> {
    let eps = tol.confusion();
    let pieces = pieces(planar, seats, blended, corners, steps)?;
    for piece in &pieces {
        let Some((_, _, planes)) = hull(&piece.points, eps) else {
            continue;
        };
        let near: Vec<usize> = piece
            .ends
            .iter()
            .filter_map(|v| faces_round(&planar.faces, *v))
            .flatten()
            .collect();
        for (fi, face) in planar.faces.iter().enumerate() {
            if near.contains(&fi) {
                continue;
            }
            // The face's outline clipped to the piece, touching counted.
            let mut ring: Vec<Point> = face.corners.iter().map(|v| planar.points[v]).collect();
            for (n, off) in &planes {
                let side = |p: &Point| n.dot(p.to_vector()) - off;
                let mut out = Vec::with_capacity(ring.len() + 1);
                for k in 0..ring.len() {
                    let (p, q) = (ring[k], ring[(k + 1) % ring.len()]);
                    let (sp, sq) = (side(&p), side(&q));
                    if sp <= eps {
                        out.push(p);
                    }
                    if (sp <= eps) != (sq <= eps) {
                        let t = (sp - eps) / (sp - sq);
                        out.push(p + (q - p) * t);
                    }
                }
                ring = out;
                if ring.is_empty() {
                    break;
                }
            }
            if !ring.is_empty() {
                return Some((piece.from.clone(), fi));
            }
        }
    }
    None
}

/// The request applied by booleans: each bevel's piece and each faceted
/// vertex's, fused where concave, then cut where convex, the history traced
/// through the booleans ([`apply`]), each piece's bevel or facet face
/// named for its edge or vertex.
/// `None` where a vertex mixes convex and concave bevels, or a piece will
/// not build or apply.
///
/// # Errors
///
/// As the booleans, other than their refusals, which decline.
#[allow(
    clippy::too_many_arguments,
    reason = "the phases' products, as P4 takes them"
)]
pub(crate) fn one_boolean(
    model: &mut Model,
    solid: &Shape,
    planar: &Planar,
    seats: &FastMap<(Id, Id), Seat>,
    blended: &ogeom_core::FastSet<Id>,
    corners: &FastMap<(usize, Id), Point>,
    steps: &FastMap<(usize, Id), [Point; 2]>,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let eps = tol.confusion();
    // Each bevel must fit its two faces on its own: a contact line past a
    // face's far side leaves the face it is seated on, which is refused by
    // name, not built.
    for seat in seats.values().filter(|s| blended.contains(&s.edge)) {
        for (k, fi) in seat.faces.iter().enumerate() {
            let face = &planar.faces[*fi];
            let inward = face
                .normal
                .cross(if k == 0 { seat.along } else { -seat.along });
            let from = planar.points[&seat.ends.0];
            let width = face
                .loops()
                .flatten()
                .map(|c| (planar.points[c] - from).dot(inward))
                .fold(0.0, f64::max);
            let Some(corner) = corners.get(&(*fi, seat.ends.0)) else {
                return Ok(None);
            };
            if (*corner - from).dot(inward) > width + eps {
                return Ok(None);
            }
        }
    }
    // Every setback corner on its own face: one standing off the face (past
    // what the part holds) is refused by name, not built.
    for ((fi, _), corner) in corners {
        let face = &planar.faces[*fi];
        let n = face.normal;
        let u = if n.x.abs() < 0.9 {
            Vector::X
        } else {
            Vector::Y
        };
        let u = u - n * u.dot(n);
        let u = u / u.magnitude();
        let w = n.cross(u);
        let flat = |p: Point| (p.to_vector().dot(u), p.to_vector().dot(w));
        let at = flat(*corner);
        let ring =
            |l: &[Id]| -> Vec<(f64, f64)> { l.iter().map(|v| flat(planar.points[v])).collect() };
        let on = |l: &[(f64, f64)]| {
            super::guards::strictly_inside(l, at, eps)
                || (0..l.len()).any(|k| {
                    let (a, b) = (l[k], l[(k + 1) % l.len()]);
                    let ab = (b.0 - a.0, b.1 - a.1);
                    let t = (((at.0 - a.0) * ab.0 + (at.1 - a.1) * ab.1)
                        / ab.0.mul_add(ab.0, ab.1 * ab.1))
                    .clamp(0.0, 1.0);
                    (a.0 + ab.0 * t - at.0).hypot(a.1 + ab.1 * t - at.1) <= eps
                })
        };
        if !on(&ring(&face.corners))
            || face
                .holes
                .iter()
                .any(|h| super::guards::strictly_inside(&ring(h), at, eps))
        {
            return Ok(None);
        }
    }
    let Some(pieces) = pieces(planar, seats, blended, corners, steps) else {
        return Ok(None);
    };
    // Each piece built as a tool, the face carrying its bevel or facet
    // named for its edge or vertex.
    let mut tools: Vec<Tool> = Vec::with_capacity(pieces.len());
    for piece in &pieces {
        let Some((at, rings, planes)) = hull(&piece.points, eps) else {
            continue;
        };
        let Ok((made, faces)) =
            ogeom_algo::make_faceted_solid(model, &at, &rings, &vec![None; rings.len()], tol)
        else {
            return Ok(None);
        };
        let carried = planes
            .iter()
            .enumerate()
            .filter(|(_, (normal, off))| {
                piece
                    .carried
                    .iter()
                    .all(|p| (normal.dot(p.to_vector()) - off).abs() <= eps)
            })
            .map(|(fi, _)| (faces[fi].clone(), piece.from.clone()))
            .collect();
        tools.push(Tool {
            solid: made.shape,
            sense: if piece.additive {
                Sense::Fuse
            } else {
                Sense::Cut
            },
            faces: carried,
        });
    }
    let Some(built) = apply(model, solid, &tools, tol)? else {
        return Ok(None);
    };
    let mut history = built.history;
    for e in blended {
        history.delete(&planar.shapes[e]);
    }
    Ok(Some(Built::new(built.shape, history)))
}

/// Whether a tool is fused into the solid (a concave blend's fill) or cut
/// from it (a convex blend's).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Sense {
    Fuse,
    Cut,
}

/// One closed tool of the fallback: its solid, its sense, and each of its
/// faces that carries a blend's face, with the edge or vertex that face
/// comes from.
pub(crate) struct Tool {
    pub solid: Shape,
    pub sense: Sense,
    pub faces: Vec<(Shape, Shape)>,
}

/// `tools` applied to `solid` by booleans, every fuse first and then every
/// cut, one boolean each (pieces that share faces are not fused into one
/// tool: the boolean refuses same-domain contact between its arguments but
/// melts a tool's faces into the solid's own). The history is the
/// booleans', with each result face a tool's named face became generated
/// by that face's edge or vertex. `None` where a boolean refuses.
///
/// # Errors
///
/// As the booleans, other than their refusals, which decline.
pub(crate) fn apply(
    model: &mut Model,
    solid: &Shape,
    tools: &[Tool],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut built = Built::from_nothing(solid.clone());
    for sense in [Sense::Fuse, Sense::Cut] {
        for tool in tools.iter().filter(|t| t.sense == sense) {
            let step = match sense {
                Sense::Fuse => ogeom_bool::fuse(model, &built.shape, &tool.solid, tol),
                Sense::Cut => ogeom_bool::cut(model, &built.shape, &tool.solid, tol),
            };
            let Ok(step) = step else {
                return Ok(None);
            };
            built = Built {
                shape: step.shape.clone(),
                history: built.history.then(&step.history),
            };
        }
    }
    let kept: Vec<Shape> = explore_unique(model, &built.shape, ShapeType::Face)?;
    let mut history: History = built.history.clone();
    for tool in tools {
        for (face, from) in &tool.faces {
            for made in built
                .history
                .modified(face)
                .iter()
                .chain(built.history.generated(face))
            {
                if kept.iter().any(|k| k.is_same(made)) {
                    history.generate(from, made.clone());
                }
            }
        }
    }
    Ok(Some(Built::new(built.shape, history)))
}
