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
/// order; `None` where a vertex's faces do not close round it once.
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
        let order = faces_round(&planar.faces, v)?;
        let carried: Vec<Point> = order.iter().map(|fi| corners[&(*fi, v)]).collect();
        // A vertex of convex and concave bevels: the corner between its
        // facet and the vertex is taken where the facet faces the vertex
        // (the vertex stands outside the result), filled where it faces
        // away. The facet faces out of the solid, as P3 winds it.
        let additive = if convex == 0 {
            true
        } else if convex == at.len() {
            false
        } else {
            let mut n = Vector::ZERO;
            for k in 0..carried.len() {
                n += carried[k]
                    .to_vector()
                    .cross(carried[(k + 1) % carried.len()].to_vector());
            }
            let out = order
                .iter()
                .fold(Vector::ZERO, |s, fi| s + planar.faces[*fi].normal);
            if n.dot(out) < 0.0 {
                n = -n;
            }
            n.dot(planar.points[&v] - carried[0]) < 0.0
        };
        let mut points = carried.clone();
        points.push(planar.points[&v]);
        pieces.push(Piece {
            points,
            carried,
            from: planar.shapes[&v].clone(),
            ends: vec![v],
            additive,
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
        if let Some(fi) = reaches(planar, &planes, &near, eps) {
            return Some((piece.from.clone(), fi));
        }
    }
    None
}

/// The first face not in `near` whose outline reaches into or touches the
/// convex piece `planes` bound.
fn reaches(planar: &Planar, planes: &[(Vector, f64)], near: &[usize], eps: f64) -> Option<usize> {
    for (fi, face) in planar.faces.iter().enumerate() {
        if near.contains(&fi) {
            continue;
        }
        // The face's outline clipped to the piece, touching counted.
        let mut ring: Vec<Point> = face.corners.iter().map(|v| planar.points[v]).collect();
        for (n, off) in planes {
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
            return Some(fi);
        }
    }
    None
}

/// A mitred request (`ChamferCorner::Point`) reaching a face it does not end
/// at, built by the one boolean: each bevel runs its edge's whole length,
/// its piece the prism of the edge and its contact lines square at the
/// edge's ends, whose union at a convex vertex is the mitre to the point
/// the bevels share. `None` where no piece reaches such a face (the rings'
/// request), or a piece will not build or apply.
///
/// # Errors
///
/// As the booleans, other than their refusals.
pub(crate) fn mitred_non_local(
    model: &mut Model,
    solid: &Shape,
    planar: &Planar,
    seats: &FastMap<(Id, Id), Seat>,
    blended: &ogeom_core::FastSet<Id>,
    setback: &FastMap<(Id, Id), f64>,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let eps = tol.confusion();
    let mut keys: Vec<&(Id, Id)> = seats.keys().collect();
    keys.sort_unstable();
    let mut reached = false;
    let mut tools = Vec::new();
    for key in keys {
        let seat = &seats[key];
        if !blended.contains(&seat.edge) {
            continue;
        }
        let (a, b) = seat.ends;
        let (pa, pb) = (planar.points[&a], planar.points[&b]);
        let mut contacts = Vec::with_capacity(4);
        for (k, fi) in seat.faces.iter().enumerate() {
            let face = &planar.faces[*fi];
            let inward = face
                .normal
                .cross(if k == 0 { seat.along } else { -seat.along });
            let Some(d) = setback.get(&(super::id(&face.shape), seat.edge)) else {
                return Ok(None);
            };
            contacts.extend([pa + inward * *d, pb + inward * *d]);
        }
        let mut points = contacts.clone();
        points.extend([pa, pb]);
        let Some((at, rings, planes)) = hull(&points, eps) else {
            return Ok(None);
        };
        let near: Vec<usize> = [a, b]
            .iter()
            .filter_map(|v| faces_round(&planar.faces, *v))
            .flatten()
            .collect();
        reached |= reaches(planar, &planes, &near, eps).is_some();
        let Ok((made, faces)) =
            ogeom_algo::make_faceted_solid(model, &at, &rings, &vec![None; rings.len()], tol)
        else {
            return Ok(None);
        };
        let named = planes
            .iter()
            .enumerate()
            .filter(|(_, (n, off))| {
                contacts
                    .iter()
                    .all(|p| (n.dot(p.to_vector()) - off).abs() <= eps)
            })
            .map(|(fi, _)| (faces[fi].clone(), planar.shapes[&seat.edge].clone()))
            .collect();
        tools.push(Tool {
            solid: made.shape,
            sense: if seat.convex { Sense::Cut } else { Sense::Fuse },
            faces: named,
        });
    }
    if !reached {
        return Ok(None);
    }
    let mut edges: Vec<Id> = blended.iter().copied().collect();
    edges.sort_unstable();
    let edges: Vec<Shape> = edges.iter().map(|e| planar.shapes[e].clone()).collect();
    applied(model, solid, &edges, &tools, tol)
}

/// Whether every blend holds to what the rings hold it to: each contact
/// line within its face's far side, each setback corner on its own face. A
/// blend past either is past what the part holds, refused by name on the
/// current path, not built.
fn held(
    planar: &Planar,
    seats: &FastMap<(Id, Id), Seat>,
    blended: &ogeom_core::FastSet<Id>,
    corners: &FastMap<(usize, Id), Point>,
    eps: f64,
) -> bool {
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
                return false;
            };
            if (*corner - from).dot(inward) > width + eps {
                return false;
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
            return false;
        }
    }
    true
}

/// `tools` applied ([`apply`]), each blended edge deleted from the history.
/// `None` where a boolean refuses, or where a blend the others leave
/// nothing of is past what the part holds (refused by name, not built).
///
/// # Errors
///
/// As [`apply`].
fn finished(
    model: &mut Model,
    solid: &Shape,
    planar: &Planar,
    blended: &ogeom_core::FastSet<Id>,
    tools: &[Tool],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut ids: Vec<Id> = blended.iter().copied().collect();
    ids.sort_unstable();
    let edges: Vec<Shape> = ids.iter().map(|e| planar.shapes[e].clone()).collect();
    applied(model, solid, &edges, tools, tol)
}

/// `tools` applied ([`apply`]), each of `edges` deleted from the history.
/// `None` where a boolean refuses, or an edge's blend is left nothing of.
///
/// # Errors
///
/// As [`apply`].
pub(crate) fn applied(
    model: &mut Model,
    solid: &Shape,
    edges: &[Shape],
    tools: &[Tool],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some(built) = apply(model, solid, tools, tol)? else {
        return Ok(None);
    };
    if edges.iter().any(|e| built.history.generated(e).is_empty()) {
        return Ok(None);
    }
    let mut history = built.history;
    for e in edges {
        history.delete(e);
    }
    Ok(Some(Built::new(built.shape, history)))
}

/// Every face [`non_local`] finds a piece reaching or touching, in face
/// order.
pub(crate) fn non_local_faces(
    planar: &Planar,
    seats: &FastMap<(Id, Id), Seat>,
    blended: &ogeom_core::FastSet<Id>,
    corners: &FastMap<(usize, Id), Point>,
    steps: &FastMap<(usize, Id), [Point; 2]>,
    tol: Tolerances,
) -> Vec<usize> {
    let eps = tol.confusion();
    let Some(pieces) = pieces(planar, seats, blended, corners, steps) else {
        return Vec::new();
    };
    let mut found = Vec::new();
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
        let mut skip: Vec<usize> = near.into_iter().chain(found.iter().copied()).collect();
        while let Some(fi) = reaches(planar, &planes, &skip, eps) {
            found.push(fi);
            skip.push(fi);
        }
    }
    found.sort_unstable();
    found
}

/// The request applied by booleans: each bevel's piece and each faceted
/// vertex's, fused where concave, then cut where convex, the history traced
/// through the booleans ([`apply`]), each piece's bevel or facet face
/// named for its edge or vertex.
/// `None` where a bevel runs past its face, a setback corner stands off
/// its face, a bevel is left nothing of, or a piece will not build or
/// apply: each refused by name on the current path.
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
    if !held(planar, seats, blended, corners, eps) {
        return Ok(None);
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
    finished(model, solid, planar, blended, &tools, tol)
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

/// Planar rounds applied by booleans (ADR 0001, M5c): each round's tool is
/// the region between its edge and its strip, closed by its two flaps (its
/// faces' planes between the edge and its contact lines) and its end caps
/// (the end faces' planes): cut where the edge is convex, fused where it is
/// concave. The strip is the engine's own, the round of the same edge on
/// the wedge its faces and end faces bound alone (the hull of the edge and
/// its setback corners, pushed twice as far out), which no other face
/// reaches: at a concave edge that wedge is the room in front of it, whose
/// edge is convex, and its round's spandrel is the concave round's fill.
/// Rounds joined at vertices, or ending where other than three faces meet,
/// close together ([`group_tool`]). `None` where a round fails what the
/// rings hold it to, or a tool will not build.
///
/// # Errors
///
/// As the booleans, other than their refusals, which decline.
#[allow(
    clippy::too_many_arguments,
    reason = "the phases' products, as P4 takes them"
)]
pub(crate) fn rounds_one_boolean(
    model: &mut Model,
    solid: &Shape,
    planar: &Planar,
    seats: &FastMap<(Id, Id), Seat>,
    by_edge: &FastMap<Id, &super::Section>,
    blended: &ogeom_core::FastSet<Id>,
    corners: &FastMap<(usize, Id), Point>,
    steps: &FastMap<(usize, Id), [Point; 2]>,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let eps = tol.confusion();
    if !steps.is_empty() || !held(planar, seats, blended, corners, eps) {
        return Ok(None);
    }
    let mut keys: Vec<&(Id, Id)> = seats.keys().collect();
    keys.sort_unstable();
    let rounded: Vec<&Seat> = keys
        .iter()
        .map(|k| &seats[*k])
        .filter(|s| matches!(by_edge.get(&s.edge), Some(super::Section::Round { .. })))
        .collect();
    // Each contact on its own face, end to end: a contact running where its
    // face does not reach (a wall cut away under a boss) leaves the ball no
    // seat there, which is not a strip crossing a face but a ball leaving
    // its seat, refused by name, not built.
    for seat in &rounded {
        for fi in &seat.faces {
            let (p, q) = (corners[&(*fi, seat.ends.0)], corners[&(*fi, seat.ends.1)]);
            if !contact_on_face(planar, *fi, p, q, eps) {
                return Ok(None);
            }
        }
    }
    // The rounds joined at their vertices, in groups.
    let mut group: Vec<usize> = (0..rounded.len()).collect();
    fn top(group: &mut [usize], k: usize) -> usize {
        let up = group[k];
        if up == k {
            return k;
        }
        let t = top(group, up);
        group[k] = t;
        t
    }
    for i in 0..rounded.len() {
        for j in i + 1..rounded.len() {
            let (a, b) = (rounded[i].ends, rounded[j].ends);
            if [a.0, a.1].iter().any(|v| *v == b.0 || *v == b.1) {
                let (x, y) = (top(&mut group, i), top(&mut group, j));
                group[x] = y;
            }
        }
    }
    let mut tools = Vec::new();
    let mut tops: Vec<usize> = (0..rounded.len()).map(|k| top(&mut group, k)).collect();
    tops.sort_unstable();
    tops.dedup();
    for t in tops {
        let members: Vec<&Seat> = (0..rounded.len())
            .filter(|k| top(&mut group, *k) == t)
            .map(|k| rounded[k])
            .collect();
        let alone = members.len() == 1
            && [members[0].ends.0, members[0].ends.1]
                .iter()
                .all(|v| faces_round(&planar.faces, *v).is_some_and(|f| f.len() == 3));
        let tool = if alone {
            let seat = members[0];
            let Some(super::Section::Round { radius }) = by_edge.get(&seat.edge) else {
                return Ok(None);
            };
            round_tool(model, planar, seat, *radius, corners, tol)?
        } else {
            group_tool(model, planar, seats, by_edge, &members, corners, tol)?
        };
        let Some(tool) = tool else {
            return Ok(None);
        };
        tools.push(tool);
    }
    finished(model, solid, planar, blended, &tools, tol)
}

/// Whether the segment `p` to `q` lies on face `fi`: its middle within the
/// face, and the segment crossing none of the face's sides (touching one is
/// on it).
fn contact_on_face(planar: &Planar, fi: usize, p: Point, q: Point, eps: f64) -> bool {
    let face = &planar.faces[fi];
    if !super::rounds::within_face(planar, fi, p + (q - p) * 0.5, eps) {
        return false;
    }
    let n = face.normal;
    let x = q - p;
    let length = x.magnitude();
    if length <= eps {
        return true;
    }
    let x = x / length;
    let y = n.cross(x);
    let flat = |c: Point| ((c - p).dot(x), (c - p).dot(y));
    for walk in face.loops() {
        for k in 0..walk.len() {
            let (a, b) = (
                flat(planar.points[&walk[k]]),
                flat(planar.points[&walk[(k + 1) % walk.len()]]),
            );
            // The side crosses the contact's line strictly between its ends,
            // from one side of it to the other.
            if (a.1 > eps && b.1 < -eps) || (a.1 < -eps && b.1 > eps) {
                let at = a.0 + (b.0 - a.0) * (a.1 / (a.1 - b.1));
                if at > eps && at < length - eps {
                    return false;
                }
            }
        }
    }
    true
}

/// One round's tool ([`rounds_one_boolean`]), its strip named for its edge;
/// `None` where it will not build.
fn round_tool(
    model: &mut Model,
    planar: &Planar,
    seat: &Seat,
    radius: f64,
    corners: &FastMap<(usize, Id), Point>,
    tol: Tolerances,
) -> OgeomResult<Option<Tool>> {
    use ogeom_geom::{LineCurve, PlaneSurface};
    use ogeom_math::{Direction, Plane};
    let eps = tol.confusion();
    let (a, b) = seat.ends;
    let (pa, pb) = (planar.points[&a], planar.points[&b]);
    let [Some(c0a), Some(c0b), Some(c1a), Some(c1b)] = [
        corners.get(&(seat.faces[0], a)),
        corners.get(&(seat.faces[0], b)),
        corners.get(&(seat.faces[1], a)),
        corners.get(&(seat.faces[1], b)),
    ] else {
        return Ok(None);
    };
    // The wedge its faces and end faces bound alone, its edge `a` to `b`.
    let wide = [
        pa,
        pb,
        pa + (*c0a - pa) * 2.0,
        pb + (*c0b - pb) * 2.0,
        pa + (*c1a - pa) * 2.0,
        pb + (*c1b - pb) * 2.0,
    ];
    let Some((at, rings, _)) = hull(&wide, eps) else {
        return Ok(None);
    };
    let Ok((wedge, _)) =
        ogeom_algo::make_faceted_solid(model, &at, &rings, &vec![None; rings.len()], tol)
    else {
        return Ok(None);
    };
    let wedge = wedge.shape;
    let point = |model: &Model, v: &Shape| -> OgeomResult<Point> {
        Ok(model
            .node(v)
            .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a vertex not in the model"))?
            .data()
            .as_vertex()
            .map_or(Point::ORIGIN, |d| d.point))
    };
    let mut edge = None;
    for e in explore_unique(model, &wedge, ShapeType::Edge)? {
        if let Some((p, q)) = ogeom_algo::edge_vertices(model, &e)? {
            let (p, q) = (point(model, &p)?, point(model, &q)?);
            if (p.distance(pa) <= eps && q.distance(pb) <= eps)
                || (p.distance(pb) <= eps && q.distance(pa) <= eps)
            {
                edge = Some(e);
            }
        }
    }
    let Some(edge) = edge else {
        return Ok(None);
    };
    let request = super::Request {
        blends: vec![super::request::Blend {
            edge: edge.clone(),
            section: super::Section::Round { radius },
        }],
        mitre: false,
        splits: Vec::new(),
    };
    let Ok(Some(rounded)) = super::run(model, &wedge, &request, tol) else {
        return Ok(None);
    };
    let [strip] = rounded.history.generated(&edge) else {
        return Ok(None);
    };
    let strip = strip.clone();
    // The strip's sides: its contact on each face's plane, and its end on
    // each end face's plane, with the contacts' ends there.
    let on = |p: Point, at: Point, n: Vector| (p - at).dot(n).abs() <= 10.0 * eps;
    let unit = |n: Vector| n / n.magnitude().max(f64::MIN_POSITIVE);
    let caps = [
        (pa, unit((*c0a - pa).cross(*c1a - pa))),
        (pb, unit((*c0b - pb).cross(*c1b - pb))),
    ];
    let mut contacts: [Option<Shape>; 2] = [None, None];
    let mut ends: [Option<Shape>; 2] = [None, None];
    let mut corner_at: [[Option<Shape>; 2]; 2] = [[None, None], [None, None]];
    let sides = explore_unique(model, &strip, ShapeType::Edge)?;
    if sides.len() != 4 {
        return Ok(None);
    }
    for e in sides {
        let Some((p, q)) = ogeom_algo::edge_vertices(model, &e)? else {
            return Ok(None);
        };
        let (pp, pq) = (point(model, &p)?, point(model, &q)?);
        if let Some(k) =
            (0..2).find(|k| on(pp, pa, seat.normals[*k]) && on(pq, pa, seat.normals[*k]))
        {
            contacts[k] = Some(e);
        } else if let Some(j) =
            (0..2).find(|j| on(pp, caps[*j].0, caps[*j].1) && on(pq, caps[*j].0, caps[*j].1))
        {
            for (v, pv) in [(p, pp), (q, pq)] {
                if let Some(k) = (0..2).find(|k| on(pv, pa, seat.normals[*k])) {
                    corner_at[j][k] = Some(v);
                }
            }
            ends[j] = Some(e);
        } else {
            return Ok(None);
        }
    }
    let ([Some(l0), Some(l1)], [Some(ea), Some(eb)]) = (contacts, ends) else {
        return Ok(None);
    };
    let [[Some(v0a), Some(v1a)], [Some(v0b), Some(v1b)]] = corner_at else {
        return Ok(None);
    };
    // The flaps and caps, over the strip's own sides.
    let segment = |model: &mut Model, from: &Shape, to: &Shape| -> OgeomResult<Shape> {
        let (p, q) = (point(model, from)?, point(model, to)?);
        let line = LineCurve::segment(p, q, tol)?;
        let range = ogeom_geom::Curve3d::domain(&line);
        Ok(ogeom_algo::make_edge_between(model, line.into(), range, from, to, tol)?.shape)
    };
    let va = ogeom_algo::make_vertex(model, pa).shape;
    let vb = ogeom_algo::make_vertex(model, pb).shape;
    let ab = segment(model, &va, &vb)?;
    let a0 = segment(model, &va, &v0a)?;
    let a1 = segment(model, &va, &v1a)?;
    let b0 = segment(model, &vb, &v0b)?;
    let b1 = segment(model, &vb, &v1b)?;
    let middle = wide.iter().fold(Vector::ZERO, |s, p| s + p.to_vector()) * (1.0 / 6.0);
    let face = |model: &mut Model, at: Point, n: Vector, sides: &[Shape]| -> OgeomResult<Shape> {
        let n = if n.dot(middle - at.to_vector()) > 0.0 {
            -n
        } else {
            n
        };
        // Walked counter-clockwise about the face's outward normal.
        let mut ring = ogeom_algo::order_edges(model, sides, tol)?;
        let mut corners = Vec::with_capacity(ring.len());
        for k in 0..ring.len() {
            let (Some((p, q)), Some((r, t))) = (
                ogeom_algo::edge_vertices(model, &ring[k])?,
                ogeom_algo::edge_vertices(model, &ring[(k + 1) % ring.len()])?,
            ) else {
                ogeom_core::ogeom_bail!(Construction, "a round's tool side has no ends");
            };
            let shared = [&p, &q]
                .into_iter()
                .find(|v| v.is_same(&r) || v.is_same(&t));
            corners.push(point(model, shared.unwrap_or(&q))?);
        }
        let mut turn = Vector::ZERO;
        for k in 0..corners.len() {
            turn += (corners[k] - corners[0]).cross(corners[(k + 1) % corners.len()] - corners[0]);
        }
        if turn.dot(n) < 0.0 {
            ring = ring.iter().rev().map(Shape::reversed).collect();
        }
        let plane = PlaneSurface::new(Plane::through(at, Direction::new(n, tol)?));
        Ok(ogeom_algo::make_face_with_pcurves(model, plane.into(), &[ring], tol)?.shape)
    };
    let carried = strip.reversed();
    let mut faces = vec![carried.clone()];
    faces.push(face(
        model,
        pa,
        seat.normals[0],
        &[ab.clone(), b0.clone(), l0, a0.clone()],
    )?);
    faces.push(face(
        model,
        pa,
        seat.normals[1],
        &[ab, b1.clone(), l1, a1.clone()],
    )?);
    faces.push(face(model, pa, caps[0].1, &[a0, ea, a1])?);
    faces.push(face(model, pb, caps[1].1, &[b0, eb, b1])?);
    let shell = ogeom_algo::make_shell(model, &faces)?.shape;
    if !ogeom_algo::is_shell_closed(model, &shell)? {
        return Ok(None);
    }
    let made = ogeom_algo::make_solid(model, &[shell])?.shape;
    Ok(Some(Tool {
        solid: made,
        sense: if seat.convex { Sense::Cut } else { Sense::Fuse },
        faces: vec![(carried, planar.shapes[&seat.edge].clone())],
    }))
}

/// One tool for rounds joined at their vertices, or a round ending where
/// other than three faces meet, where the solid is convex round every
/// vertex they end at: the convex polytope the faces round those vertices
/// bound (within a box well past the rounds' setbacks), less that polytope
/// rounded by the engine (one boolean between the two, sharing their faces
/// but the rounds'). The rounds' strips and corners are the engine's own,
/// as on the solid where nothing else reaches them; each strip named for
/// its edge, each corner for its vertex. `None` where a vertex is not
/// convex all round, or the polytope will not build or round.
fn group_tool(
    model: &mut Model,
    planar: &Planar,
    seats: &FastMap<(Id, Id), Seat>,
    by_edge: &FastMap<Id, &super::Section>,
    members: &[&Seat],
    corners: &FastMap<(usize, Id), Point>,
    tol: Tolerances,
) -> OgeomResult<Option<Tool>> {
    let eps = tol.confusion();
    let mut vertices: Vec<Id> = members.iter().flat_map(|s| [s.ends.0, s.ends.1]).collect();
    vertices.sort_unstable();
    vertices.dedup();
    // Convex all round each vertex: the faces there bound the solid's own
    // corner, a convex cone.
    if seats
        .values()
        .any(|s| !s.convex && (vertices.contains(&s.ends.0) || vertices.contains(&s.ends.1)))
    {
        return Ok(None);
    }
    let mut faces: Vec<usize> = Vec::new();
    for v in &vertices {
        let Some(round) = faces_round(&planar.faces, *v) else {
            return Ok(None);
        };
        faces.extend(round);
    }
    faces.sort_unstable();
    faces.dedup();
    let mut planes: Vec<(Vector, f64)> = faces
        .iter()
        .map(|fi| {
            let face = &planar.faces[*fi];
            let n = face.normal;
            (n, n.dot(planar.points[&face.corners[0]].to_vector()))
        })
        .collect();
    // A box well past the rounds' setbacks.
    let mut reach: f64 = 0.0;
    let (mut lo, mut hi) = (
        Vector::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
        Vector::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for seat in members {
        for v in [seat.ends.0, seat.ends.1] {
            let p = planar.points[&v].to_vector();
            lo = Vector::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Vector::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
            for fi in &seat.faces {
                if let Some(c) = corners.get(&(*fi, v)) {
                    reach = reach.max(c.distance(planar.points[&v]));
                }
            }
        }
    }
    let margin = 3.0 * reach + 10.0 * eps;
    for (axis, k) in [(Vector::X, 0), (Vector::Y, 1), (Vector::Z, 2)] {
        let (l, h) = match k {
            0 => (lo.x, hi.x),
            1 => (lo.y, hi.y),
            _ => (lo.z, hi.z),
        };
        planes.push((axis, h + margin));
        planes.push((-axis, -(l - margin)));
    }
    // The polytope's corners: where three planes meet within all of them.
    let mut points: Vec<Point> = Vec::new();
    for i in 0..planes.len() {
        for j in i + 1..planes.len() {
            for k in j + 1..planes.len() {
                let ((a, da), (b, db), (c, dc)) = (planes[i], planes[j], planes[k]);
                let det = a.dot(b.cross(c));
                if det.abs() <= 1e-12 {
                    continue;
                }
                let p = (b.cross(c) * da + c.cross(a) * db + a.cross(b) * dc) / det;
                if planes.iter().all(|(n, d)| n.dot(p) <= d + eps) {
                    points.push(Point::ORIGIN + p);
                }
            }
        }
    }
    let Some((at, rings, _)) = hull(&points, eps) else {
        return Ok(None);
    };
    let Ok((made, _)) =
        ogeom_algo::make_faceted_solid(model, &at, &rings, &vec![None; rings.len()], tol)
    else {
        return Ok(None);
    };
    let helper = made.shape;
    // The helper's edges and vertices standing where the rounds' do.
    let point_of = |model: &Model, v: &Shape| -> Option<Point> {
        model
            .node(v)
            .and_then(|n| n.data().as_vertex().map(|d| d.point))
    };
    let mut asked: Vec<(Shape, f64)> = Vec::with_capacity(members.len());
    let mut names: Vec<(Shape, Shape)> = Vec::new();
    for seat in members {
        let Some(super::Section::Round { radius }) = by_edge.get(&seat.edge) else {
            return Ok(None);
        };
        let (pa, pb) = (planar.points[&seat.ends.0], planar.points[&seat.ends.1]);
        let mut found = None;
        for e in explore_unique(model, &helper, ShapeType::Edge)? {
            if let Some((p, q)) = ogeom_algo::edge_vertices(model, &e)? {
                let (Some(p), Some(q)) = (point_of(model, &p), point_of(model, &q)) else {
                    continue;
                };
                if (p.distance(pa) <= eps && q.distance(pb) <= eps)
                    || (p.distance(pb) <= eps && q.distance(pa) <= eps)
                {
                    found = Some(e);
                }
            }
        }
        let Some(e) = found else {
            return Ok(None);
        };
        names.push((e.clone(), planar.shapes[&seat.edge].clone()));
        asked.push((e, *radius));
    }
    for v in &vertices {
        let p = planar.points[v];
        for hv in explore_unique(model, &helper, ShapeType::Vertex)? {
            if point_of(model, &hv).is_some_and(|q| q.distance(p) <= eps) {
                names.push((hv, planar.shapes[v].clone()));
            }
        }
    }
    let request = super::Request::from_rounds(&asked, tol)?;
    let Ok(Some(rounded)) = super::run(model, &helper, &request, tol) else {
        return Ok(None);
    };
    let Ok(spandrel) = ogeom_bool::cut(model, &helper, &rounded.shape, tol) else {
        return Ok(None);
    };
    // Each face of the tool from a round's strip or corner, named for the
    // solid's edge or vertex it comes from.
    let kept = explore_unique(model, &spandrel.shape, ShapeType::Face)?;
    let mut faces_named: Vec<(Shape, Shape)> = Vec::new();
    for (from_helper, from_solid) in &names {
        for blend in rounded.history.generated(from_helper) {
            let images: Vec<Shape> = spandrel
                .history
                .modified(blend)
                .iter()
                .chain(spandrel.history.generated(blend))
                .cloned()
                .chain(std::iter::once(blend.clone()))
                .collect();
            for image in images {
                if let Some(k) = kept.iter().find(|k| k.is_same(&image)) {
                    faces_named.push((k.clone(), from_solid.clone()));
                }
            }
        }
    }
    if explore_unique(model, &spandrel.shape, ShapeType::Solid)?.is_empty() {
        return Ok(None);
    }
    Ok(Some(Tool {
        solid: spandrel.shape,
        sense: Sense::Cut,
        faces: faces_named,
    }))
}

/// Straight bevels between planes on a solid with curved faces, built by
/// the one boolean where one reaches a face it does not end at (ADR 0001,
/// amendment "M5a"): each bevel's piece is the convex hull of its edge and
/// its contact lines' ends on the planar faces at its two ends, and a face
/// not touching those ends that cannot be proved outside that hull makes
/// the request non-local. `None` (the curved path's) where an asked edge is
/// not straight between planes, two asked edges meet, an end's face is not
/// a plane, or no face is reached.
///
/// # Errors
///
/// As reading the solid and the booleans, other than their refusals.
pub(crate) fn straight_non_local(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, crate::Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    use ogeom_algo::HalfSpaces;
    let eps = tol.confusion();
    // Each bevel: its edge, its hull's points, its contact lines' ends, its
    // end vertices, and whether it fills.
    type Straight = (Shape, Vec<Point>, Vec<Point>, Vec<Shape>, bool);
    let mut pieces: Vec<Straight> = Vec::new();
    let mut ends_seen: Vec<Shape> = Vec::new();
    for (edge, spec) in specs {
        let Ok(crate::chamfer::Legs::Straight(seat, d)) =
            crate::chamfer::legs(model, solid, edge, spec, tol)
        else {
            return Ok(None);
        };
        let Some((va, vb)) = ogeom_algo::edge_vertices(model, edge)? else {
            return Ok(None);
        };
        if ends_seen.iter().any(|v| v.is_same(&va) || v.is_same(&vb)) {
            return Ok(None);
        }
        ends_seen.extend([va.clone(), vb.clone()]);
        let legs = [seat.leg(0, tol)?, seat.leg(1, tol)?];
        let mut carried = Vec::with_capacity(4);
        let mut points = Vec::with_capacity(6);
        for (v, at) in [(&va, seat.start), (&vb, seat.end)] {
            // The face at this end that is neither of the seat's: a plane
            // the contact lines end on.
            let around = ogeom_topo::ancestors_of(model, solid, v, ShapeType::Face)?;
            let end: Vec<&Shape> = around
                .iter()
                .filter(|f| !seat.faces.iter().any(|g| g.is_same(f)))
                .collect();
            let [end] = end[..] else {
                return Ok(None);
            };
            let Some((n, off)) = face_plane(model, end, tol)? else {
                return Ok(None);
            };
            for k in 0..2 {
                let from = at + legs[k] * d[k];
                let along = n.dot(seat.along);
                if along.abs() <= tol.angular() {
                    return Ok(None);
                }
                carried.push(from + seat.along * ((off - n.dot(from.to_vector())) / along));
            }
            points.push(at);
        }
        points.extend(carried.iter().copied());
        pieces.push((edge.clone(), points, carried, vec![va, vb], !seat.convex));
    }
    // A face not touching a piece's ends that is not proved outside it.
    let faces = explore_unique(model, solid, ShapeType::Face)?;
    let mut reached = false;
    let mut tools = Vec::with_capacity(pieces.len());
    for (edge, points, carried, ends, additive) in &pieces {
        let Some((at, rings, planes)) = hull(points, eps) else {
            return Ok(None);
        };
        let near: Vec<Shape> = ends
            .iter()
            .map(|v| ogeom_topo::ancestors_of(model, solid, v, ShapeType::Face))
            .collect::<OgeomResult<Vec<_>>>()?
            .into_iter()
            .flatten()
            .collect();
        let outside = HalfSpaces::Any(
            planes
                .iter()
                .map(|(n, off)| HalfSpaces::Below {
                    normal: -*n,
                    at: Point::ORIGIN + *n * *off,
                })
                .collect(),
        );
        for face in &faces {
            if near.iter().any(|f| f.is_same(face)) {
                continue;
            }
            if ogeom_algo::certify_face_within(model, face, &outside, eps, tol)?.is_err() {
                reached = true;
            }
        }
        let Ok((made, made_faces)) =
            ogeom_algo::make_faceted_solid(model, &at, &rings, &vec![None; rings.len()], tol)
        else {
            return Ok(None);
        };
        let named = planes
            .iter()
            .enumerate()
            .filter(|(_, (n, off))| {
                carried
                    .iter()
                    .all(|p| (n.dot(p.to_vector()) - off).abs() <= eps)
            })
            .map(|(fi, _)| (made_faces[fi].clone(), edge.clone()))
            .collect();
        tools.push(Tool {
            solid: made.shape,
            sense: if *additive { Sense::Fuse } else { Sense::Cut },
            faces: named,
        });
    }
    if !reached {
        return Ok(None);
    }
    let Some(built) = apply(model, solid, &tools, tol)? else {
        return Ok(None);
    };
    if specs
        .iter()
        .any(|(e, _)| built.history.generated(e).is_empty())
    {
        return Ok(None);
    }
    let mut history = built.history;
    for (e, _) in specs {
        history.delete(e);
    }
    Ok(Some(Built::new(built.shape, history)))
}

/// A planar face's outward unit normal and offset; `None` on another
/// surface.
fn face_plane(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<Option<(Vector, f64)>> {
    let planar = model
        .node(face)
        .and_then(|n| n.data().as_face())
        .and_then(|d| model.geometry().surface(d.surface))
        .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Plane(_)));
    if !planar {
        return Ok(None);
    }
    let (at, n) = ogeom_algo::face_normal(model, face, tol)?;
    Ok(Some((n, n.dot(at.to_vector()))))
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
    use super::super::{Request, run};
    use core::f64::consts::FRAC_PI_4;
    use ogeom_core::Tolerances;
    use ogeom_math::{Direction, Frame, Point};
    use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    fn block(m: &mut Model, at: (f64, f64, f64), size: (f64, f64, f64)) -> Shape {
        let frame =
            Frame::new(Point::new(at.0, at.1, at.2), Direction::Z, Direction::X, T).unwrap();
        ogeom_algo::make_box(m, frame, size, T).unwrap().shape
    }

    fn volume(m: &Model, s: &Shape) -> f64 {
        ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
            .unwrap()
            .mass
    }

    /// The edge of `s` along `x` at `y`, `z`, its whole length `length`.
    fn edge_along(m: &Model, s: &Shape, y: f64, z: f64, length: f64) -> Shape {
        explore_unique(m, s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let bb = ogeom_algo::shape_bounds(m, e, T).unwrap();
                let (l, h) = (bb.low().unwrap(), bb.high().unwrap());
                (l.y - y).abs() < 1e-6
                    && (h.y - y).abs() < 1e-6
                    && (l.z - z).abs() < 1e-6
                    && (h.z - z).abs() < 1e-6
                    && (h.x - l.x - length).abs() < 1e-5
            })
            .unwrap()
    }

    /// The area under `sqrt(r^2 - u^2)` from nought to `u`.
    fn under(r: f64, u: f64) -> f64 {
        0.5 * (u * (r * r - u * u).sqrt() + r * r * (u / r).asin())
    }

    /// The round built by the engine, through the one boolean: valid, the
    /// edge generating one strip, the volume `want`.
    fn rounds_to(m: &mut Model, s: &Shape, edge: &Shape, r: f64, want: f64) {
        let request = Request::from_rounds(&[(edge.clone(), r)], T).unwrap();
        let built = run(m, s, &request, T)
            .unwrap()
            .expect("built by the engine");
        let found = ogeom_algo::check(m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        assert_eq!(built.history.generated(edge).len(), 1, "the strip");
        let got = volume(m, &built.shape);
        assert!((got - want).abs() < 1e-6, "{got} against {want}");
    }

    /// A plate 4 thick, its plan the rectangle 20 by 10 or the trapezoid
    /// narrowing to 14 at its back (its sides leaning in 3 over 10, so the
    /// round's ends lie on them askew), with a square pocket 1 by 1/2
    /// through it, 1 in from its front; its top front edge rounded at 2.
    /// The round runs through the pocket, whose walls it crosses, faces not
    /// round its ends. Built by the engine through the one boolean: the
    /// plate less the pocket and the round's spandrel, less what the two
    /// share. The spandrel's section `A = (1 - pi/4) r^2` runs `20 - 2 l y`
    /// long at `y` behind the front (`l` the sides' lean), so it takes
    /// `20 A - 2 l (20/3 - 2 pi)`; across the pocket its height is
    /// `2 - sqrt(4 - u^2)` at `u = 2 - y`.
    #[test]
    fn a_round_crossing_a_pocket_goes_through_the_one_boolean() {
        for lean in [0.0, 0.3] {
            let mut m = Model::new();
            let corners = [
                Point::new(0.0, 0.0, 0.0),
                Point::new(20.0, 0.0, 0.0),
                Point::new(20.0 - 10.0 * lean, 10.0, 0.0),
                Point::new(10.0 * lean, 10.0, 0.0),
            ];
            let wire = ogeom_algo::make_polygon(&mut m, &corners, true, T)
                .unwrap()
                .shape;
            let plane = ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(Frame::WORLD));
            let base = ogeom_algo::make_face(&mut m, plane.into(), &[wire], T)
                .unwrap()
                .shape;
            let plate =
                ogeom_algo::make_prism(&mut m, &base, ogeom_math::Vector::new(0.0, 0.0, 4.0), T)
                    .unwrap()
                    .shape;
            let pocket = block(&mut m, (8.0, 1.0, -1.0), (1.0, 0.5, 6.0));
            let s = ogeom_bool::cut(&mut m, &plate, &pocket, T).unwrap().shape;
            let edge = edge_along(&m, &s, 0.0, 4.0, 20.0);
            let r = 2.0;
            let spandrel = (1.0 - FRAC_PI_4) * r * r * 20.0
                - 2.0 * lean * (20.0 / 3.0 - 2.0 * core::f64::consts::PI);
            let shared = 1.0 - (under(r, 1.0) - under(r, 0.5));
            let plate = (20.0 + 20.0 - 20.0 * lean) * 5.0 * 4.0;
            let want = plate - 2.0 - spandrel + shared;
            rounds_to(&mut m, &s, &edge, r, want);
        }
    }

    /// A 20 by 10 by 6 block with a square pocket 1 by 1/2 through it, 1 in
    /// from its front, the three edges at its corner (20, 0, 6) rounded at
    /// 2: the front top edge's round runs through the pocket, and the three
    /// close on the corner's sphere. One tool for the three, through the
    /// one boolean. The block less the pocket and what the rounds take
    /// (each edge's section `A = (1 - pi/4) r^2` along it short of the
    /// corner's cube, and the cube `r^3` less the sphere's octant), less
    /// what the front round and the pocket share.
    #[test]
    fn a_corner_s_rounds_crossing_a_pocket_go_through_the_one_boolean() {
        let mut m = Model::new();
        let solid = block(&mut m, (0.0, 0.0, 0.0), (20.0, 10.0, 6.0));
        let pocket = block(&mut m, (8.0, 1.0, -1.0), (1.0, 0.5, 8.0));
        let s = ogeom_bool::cut(&mut m, &solid, &pocket, T).unwrap().shape;
        let near = |p: Point| {
            explore_unique(&m, &s, ShapeType::Edge)
                .unwrap()
                .into_iter()
                .filter(|e| {
                    let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
                    let at = |v: &Shape| {
                        m.node(v)
                            .and_then(|n| n.data().as_vertex().map(|d| d.point))
                            .unwrap()
                    };
                    at(&a).distance(p) < 1e-9 || at(&b).distance(p) < 1e-9
                })
                .collect::<Vec<_>>()
        };
        let edges = near(Point::new(20.0, 0.0, 6.0));
        assert_eq!(edges.len(), 3);
        let r = 2.0;
        let asked: Vec<(Shape, f64)> = edges.iter().map(|e| (e.clone(), r)).collect();
        let request = Request::from_rounds(&asked, T).unwrap();
        let built = run(&mut m, &s, &request, T)
            .unwrap()
            .expect("built by the engine");
        let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        for e in &edges {
            assert_eq!(built.history.generated(e).len(), 1, "each strip");
        }
        let section = (1.0 - FRAC_PI_4) * r * r;
        let corner = r * r * r * (1.0 - core::f64::consts::PI / 6.0);
        let taken = section * (20.0 + 10.0 + 6.0 - 3.0 * r) + corner;
        let shared = 1.0 - (under(r, 1.0) - under(r, 0.5));
        let want = 1200.0 - 3.0 - taken + shared;
        let got = volume(&m, &built.shape);
        assert!((got - want).abs() < 1e-6, "{got} against {want}");
    }

    /// An L of a floor 20 by 10 by 2 and a wall 2 thick along its back, 10
    /// tall, with a pillar 4 by 1 by 3 standing on the floor 1/2 in front of
    /// the wall; the concave edge at the wall's foot rounded at 2. The
    /// round's fill reaches into the pillar, a face not round its ends.
    /// Built by the engine through the one boolean: the L and the pillar
    /// with the fill `(1 - pi/4) r^2` along the edge, less what the fill
    /// and the pillar share (the fill's height `2 - sqrt(4 - u^2)` at `u =
    /// y - 6` in front of the ball's centre, across the pillar).
    #[test]
    fn a_concave_round_reaching_a_pillar_goes_through_the_one_boolean() {
        let mut m = Model::new();
        let floor = block(&mut m, (0.0, 0.0, 0.0), (20.0, 10.0, 2.0));
        let wall = block(&mut m, (0.0, 8.0, 0.0), (20.0, 2.0, 10.0));
        let pillar = block(&mut m, (8.0, 6.5, 2.0), (4.0, 1.0, 3.0));
        let l = ogeom_bool::fuse(&mut m, &floor, &wall, T).unwrap().shape;
        let s = ogeom_bool::fuse(&mut m, &l, &pillar, T).unwrap().shape;
        let (s, _) = ogeom_heal::unify_same_domain(&mut m, &s, T).unwrap();
        let s = s.shape;
        let edge = edge_along(&m, &s, 8.0, 2.0, 20.0);
        let r = 2.0;
        let fill = (1.0 - FRAC_PI_4) * r * r;
        let shared = 4.0 * (2.0 - (under(r, 1.5) - under(r, 0.5)));
        let want = 400.0 + 20.0 * 2.0 * 8.0 + 12.0 + fill * 20.0 - shared;
        rounds_to(&mut m, &s, &edge, r, want);
    }
}
