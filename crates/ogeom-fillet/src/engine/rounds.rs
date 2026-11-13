//! P3 and P4 for rounds between planes: the strips, the vertex blends and
//! the faces set back to their corners, built as one sheet over shared
//! points and edges.
//!
//! P2 has set every face back to its setback corners. Each round's strip
//! is its section's cylinder between its contact lines, ended by its arc
//! at each end: where both its faces' setback corners there stand at one
//! station along the edge, its arc runs between them. Each vertex closes by
//! its [`RoundVertex`]: a sphere through the three strips' arcs, or a flush
//! end, where a strip's arc lies in the vertex's third face and becomes
//! that face's corner. A vertex no method closes yet leaves the request to
//! the current path.

use ogeom_algo::{Built, HalfSpaces, History};
use ogeom_core::{FastMap, OgeomResult, Tolerances};
use ogeom_geom::{Curve, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Orientation, Shape, ShapeType, explore_unique};

use super::Id;
use super::classify::{Planar, Seat, faces_round};
use super::curved::{SheetEdge, SheetFace, build_faces};
use super::request::Section;
use super::vertices::RoundVertex;
use crate::round::section::{PlanarSection, RoundSpec};

/// The rounded solid of `by_edge` (every section a round) on `planar`, its
/// faces set back to `corners`; `None` where a section is not a round, a
/// strip's ends do not stand square to its edge, or a vertex is not one
/// [`RoundVertex`] closes.
///
/// # Errors
///
/// As the strips', faces' and guards' own: a face that cannot be trimmed,
/// or a solid failing its check.
pub(crate) fn assemble(
    model: &mut Model,
    planar: &Planar,
    seats: &FastMap<(Id, Id), Seat>,
    by_edge: &FastMap<Id, &Section>,
    corners: &FastMap<(usize, Id), Point>,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut keys: Vec<(Id, Id)> = seats.keys().copied().collect();
    keys.sort_unstable();
    // Each round's section, keyed by its edge.
    let mut sections: FastMap<Id, (PlanarSection, f64)> = FastMap::default();
    for key in &keys {
        let seat = &seats[key];
        let Some(section) = by_edge.get(&seat.edge) else {
            continue;
        };
        let Section::Round { radius } = section else {
            return Ok(None);
        };
        let support = crate::support::Seat {
            start: planar.points[&seat.ends.0],
            end: planar.points[&seat.ends.1],
            along: seat.along,
            normals: seat.normals,
            faces: [
                planar.faces[seat.faces[0]].shape.clone(),
                planar.faces[seat.faces[1]].shape.clone(),
            ],
            convex: seat.convex,
        };
        let solved = PlanarSection::new(&support, &RoundSpec::Circular(*radius), tol)?;
        sections.insert(seat.edge, (solved, *radius));
    }
    // An edge asked that is not one of the solid's is the current path's to
    // refuse.
    if sections.len() != by_edge.len() {
        return Ok(None);
    }
    let edge_seat: FastMap<Id, &Seat> = keys.iter().map(|k| (seats[k].edge, &seats[k])).collect();

    // Each strip's end at each of its vertices: the station along its edge
    // of each face's setback corner there, and their middle. Where they
    // stand level the end is the section's arc there; else the strip ends
    // skew, through both corners.
    let reach = tol.confusion() * 100.0;
    let mut level: FastMap<(Id, Id), [f64; 2]> = FastMap::default();
    let mut station: FastMap<(Id, Id), f64> = FastMap::default();
    for (edge, (section, _)) in &sections {
        let seat = edge_seat[edge];
        for v in [seat.ends.0, seat.ends.1] {
            let s = [0, 1]
                .map(|k| (corners[&(seat.faces[k], v)] - section.crease(0.0)).dot(seat.along));
            level.insert((*edge, v), s);
            station.insert((*edge, v), 0.5 * (s[0] + s[1]));
        }
    }
    let stands_level = |edge: Id, v: Id| {
        let s = level[&(edge, v)];
        (s[0] - s[1]).abs() <= reach
    };

    // A face set back to a line or a point, or past one of its sides, is
    // not a class this builds.
    for (fi, face) in planar.faces.iter().enumerate() {
        let n = face.corners.len();
        for k in 0..n {
            let (v, w) = (face.corners[k], face.corners[(k + 1) % n]);
            let side = planar.points[&w] - planar.points[&v];
            let kept = corners[&(fi, w)] - corners[&(fi, v)];
            if kept.dot(side) <= reach * side.magnitude() {
                return Ok(None);
            }
        }
    }

    // Each vertex's method.
    let mut vertices: Vec<Id> = planar.points.keys().copied().collect();
    vertices.sort_unstable();
    let mut method: FastMap<Id, RoundVertex> = FastMap::default();
    // For a flush vertex, its third face and the round ending there.
    let mut flush: FastMap<Id, (usize, Id)> = FastMap::default();
    // For a mitre vertex, its two rounds and the mitre plane's normal.
    let mut mitres: FastMap<Id, ([Id; 2], Vector)> = FastMap::default();
    // For two rounds beside an edge left sharp that meet otherwise.
    let mut creases: FastMap<Id, Crease> = FastMap::default();
    // For a vertex the ball rounds rolling along a ridge.
    let mut ridges: FastMap<Id, Ridge> = FastMap::default();
    for &v in &vertices {
        let Some(round) = faces_round(&planar.faces, v) else {
            return Ok(None);
        };
        let here: Vec<&Seat> = keys
            .iter()
            .filter(|(a, b)| *a == v || *b == v)
            .map(|k| &seats[k])
            .collect();
        let blended: Vec<&Seat> = here
            .iter()
            .copied()
            .filter(|s| sections.contains_key(&s.edge))
            .collect();
        let chosen = match (here.len(), blended.len()) {
            (_, 0) => RoundVertex::Sharp,
            (3, 1) => {
                // The strip ends in the third face, cut by its plane: square
                // to the edge, its section's arc; leaning, the ellipse.
                let seat = blended[0];
                let Some(third) = round.iter().copied().find(|f| !seat.faces.contains(f)) else {
                    return Ok(None);
                };
                let normal = planar.faces[third].normal;
                let square = normal.cross(seat.along).magnitude() <= tol.angular();
                let ends = if square {
                    let s = station[&(seat.edge, v)];
                    let at_end =
                        s.abs() <= reach || (s - section_length(planar, seat)).abs() <= reach;
                    at_end && stands_level(seat.edge, v)
                } else {
                    // Leaning, the strip's contacts end at different
                    // stations, each corner in the third face's plane and
                    // within its own face, and the end within the third
                    // face; a face nearly holding the edge's direction cuts
                    // no end. A strip whose contact runs off a face narrower
                    // than its setback is left to the current path.
                    let on = |k: usize| {
                        let corner = corners[&(seat.faces[k], v)];
                        (corner - planar.points[&v]).dot(normal).abs() <= reach
                            && within_face(planar, seat.faces[k], corner, reach)
                    };
                    let cuts = normal.dot(seat.along).abs() > tol.angular() && on(0) && on(1);
                    cuts && {
                        let (section, _) = &sections[&seat.edge];
                        let order = section.order();
                        let s = level[&(seat.edge, v)];
                        let (end, (a, b)) =
                            section.plane_end(order, [s[order[0]], s[order[1]]], normal, tol)?;
                        (0..=16).all(|k| {
                            let t = a + (b - a) * f64::from(k) / 16.0;
                            ogeom_geom::Curve3d::point_at(&end, t, tol)
                                .is_ok_and(|p| within_face(planar, third, p, reach))
                        })
                    }
                };
                if !ends {
                    return Ok(None);
                }
                flush.insert(v, (third, seat.edge));
                RoundVertex::Flush
            }
            (3, 2) => {
                match mitre(planar, &sections, &here, &round, corners, v, reach, tol)? {
                    Some((pair, normal)) => {
                        mitres.insert(v, (pair, normal));
                    }
                    None => {
                        let Some(crease) =
                            crease_at(&sections, &here, &round, corners, v, planar, reach)
                        else {
                            return Ok(None);
                        };
                        creases.insert(v, crease);
                    }
                }
                RoundVertex::Mitre
            }
            // Rounds at every edge: one ball touching every face at its
            // setback corner, or the patch, at three edges or more of one
            // sense.
            (n, k) if n >= 3 && k == n => {
                let radius = sections[&blended[0].edge].1;
                let convex = blended[0].convex;
                let alike = blended.iter().all(|s| {
                    s.convex == convex && sections[&s.edge].1 == radius && stands_level(s.edge, v)
                });
                let ball = if alike {
                    one_ball(planar, &round, v, radius, convex, tol)
                } else {
                    None
                };
                let ridge = if n == 4
                    && blended
                        .iter()
                        .all(|s| s.convex && sections[&s.edge].1 == radius)
                {
                    ridge_at(planar, &round, v, radius, corners, reach, tol)
                } else {
                    None
                };
                let fits = ball.as_ref().is_some_and(|(_, feet)| {
                    round
                        .iter()
                        .zip(feet)
                        .all(|(f, foot)| foot.distance(corners[&(*f, v)]) <= reach)
                });
                if !fits && let Some(ridge) = ridge {
                    ridges.insert(v, ridge);
                    method.insert(v, RoundVertex::Ridge);
                    continue;
                }
                match ball {
                    Some((_, feet))
                        if round
                            .iter()
                            .zip(&feet)
                            .all(|(f, foot)| foot.distance(corners[&(*f, v)]) <= reach) =>
                    {
                        RoundVertex::Sphere
                    }
                    // Mixed convexity at four edges or more closes on the
                    // patch too (ADR 0001, §13, question 7), proved against
                    // the solid's own material there.
                    _ if n == 3
                        || blended.iter().all(|s| s.convex == convex)
                        || local_cone(planar, &here, &round, v).is_some() =>
                    {
                        RoundVertex::Patch
                    }
                    _ => return Ok(None),
                }
            }
            _ => return Ok(None),
        };
        method.insert(v, chosen);
    }

    // The sheet's points, shared by index.
    let mut points: Vec<Point> = Vec::new();
    let mut point = |p: Point| -> usize {
        if let Some(i) = points.iter().position(|q| q.distance(p) <= reach) {
            return i;
        }
        points.push(p);
        points.len() - 1
    };
    let at: FastMap<(usize, Id), usize> = {
        let mut sorted: Vec<((usize, Id), Point)> = corners.iter().map(|(k, p)| (*k, *p)).collect();
        sorted.sort_unstable_by_key(|(k, _)| *k);
        sorted.into_iter().map(|(k, p)| (k, point(p))).collect()
    };
    // Each ridge's balls' feet, by vertex, face and ball.
    let mut ridge_points: FastMap<(Id, usize, usize), usize> = FastMap::default();
    let mut ridge_list: Vec<(&Id, &Ridge)> = ridges.iter().collect();
    ridge_list.sort_unstable_by_key(|(v, _)| **v);
    for &(&v, ridge) in &ridge_list {
        for k in 0..2 {
            for f in [ridge.faces[0], ridge.faces[1 + 2 * k], ridge.faces[2]] {
                if let Some(foot) = ridge.foot(k, f) {
                    ridge_points.insert((v, f, k), point(foot));
                }
            }
        }
    }
    // Where face `fi`'s ring passes vertex `v`, coming from `prev`: its
    // setback corner, or at a ridge's faces either side the two balls'
    // feet, the one of the ball the arriving edge's strip ends on first.
    let passes = |fi: usize, prev: Id, v: Id| -> Option<Vec<usize>> {
        let Some(ridge) = ridges.get(&v) else {
            return Some(vec![at[&(fi, v)]]);
        };
        if fi != ridge.faces[0] && fi != ridge.faces[2] {
            return Some(vec![at[&(fi, v)]]);
        }
        let seat = seats.get(&(prev.min(v), prev.max(v)))?;
        let first = ridge.ball_of(seat)?;
        Some(vec![
            ridge_points[&(v, fi, first)],
            ridge_points[&(v, fi, 1 - first)],
        ])
    };
    // The faces' setback polygons: none crossing itself or another, none
    // turned over.
    let mut rings: Vec<super::assemble::Ring> = Vec::with_capacity(planar.faces.len());
    for (fi, face) in planar.faces.iter().enumerate() {
        let n = face.corners.len();
        let mut ring = Vec::with_capacity(n);
        for k in 0..n {
            let Some(pts) = passes(fi, face.corners[(k + n - 1) % n], face.corners[k]) else {
                return Ok(None);
            };
            ring.extend(pts);
        }
        rings.push((ring, vec![super::Origin::Face(fi)], Some(face.normal)));
    }
    super::guards::rings_fit(&points, &rings, tol)?;

    // The sheet's edges: straight segments by their ends, and each strip's
    // end arcs by edge and vertex.
    let mut edges: Vec<SheetEdge> = Vec::new();
    let mut segments: FastMap<(usize, usize), usize> = FastMap::default();
    let mut segment = |edges: &mut Vec<SheetEdge>, a: usize, b: usize| -> (usize, bool) {
        let key = (a.min(b), a.max(b));
        let index = *segments.entry(key).or_insert_with(|| {
            edges.push(SheetEdge {
                ends: key,
                curve: None,
            });
            edges.len() - 1
        });
        (index, a == key.0)
    };
    // Each strip's arc at each end: the edge's index, from the first face's
    // corner to the second's.
    let mut arcs: FastMap<(Id, Id), (usize, bool)> = FastMap::default();
    // For each crease vertex: the smaller round's corner on the far face,
    // the larger's corner on the sharp edge, and the larger's end there.
    let mut crease_corners: FastMap<Id, (usize, usize, usize)> = FastMap::default();
    let mut strip_list: Vec<Id> = sections.keys().copied().collect();
    strip_list.sort_unstable();
    for edge in &strip_list {
        let (section, _) = &sections[edge];
        let seat = edge_seat[edge];
        for v in [seat.ends.0, seat.ends.1] {
            let order = section.order();
            let [from_face, _] = order;
            let (p0, p1) = (at[&(seat.faces[0], v)], at[&(seat.faces[1], v)]);
            if creases.contains_key(&v) {
                continue;
            }
            // At a ridge the strip ends on its cap's ball, its section's arc
            // through the ball's feet.
            if let Some(ridge) = ridges.get(&v) {
                let Some(k) = ridge.ball_of(seat) else {
                    return Ok(None);
                };
                let (Some(&q0), Some(&q1)) = (
                    ridge_points.get(&(v, seat.faces[0], k)),
                    ridge_points.get(&(v, seat.faces[1], k)),
                ) else {
                    return Ok(None);
                };
                let s = (ridge.balls[k].0 - section.crease(0.0)).dot(seat.along);
                let (curve, range) = section.arc_at(s, tol)?;
                edges.push(SheetEdge {
                    ends: if from_face == 0 { (q0, q1) } else { (q1, q0) },
                    curve: Some((curve, range)),
                });
                arcs.insert((*edge, v), (edges.len() - 1, from_face == 0));
                continue;
            }
            // At a mitre the second round ends on the first's end: the two
            // equal cylinders meet in the mitre plane, on one ellipse.
            if let Some(([first, _], _)) = mitres.get(&v)
                && *first != *edge
            {
                let (shared, _) = arcs[&(*first, v)];
                let ends = edges[shared].ends;
                if !(ends == (p0, p1) || ends == (p1, p0)) {
                    return Ok(None);
                }
                arcs.insert((*edge, v), (shared, ends.0 == p0));
                continue;
            }
            // A flush end in a face leaning across the edge.
            let leaning = flush.get(&v).and_then(|(third, at_edge)| {
                let n = planar.faces[*third].normal;
                (*at_edge == *edge && n.cross(seat.along).magnitude() > tol.angular()).then_some(n)
            });
            let (curve, range) = if let Some((_, normal)) = mitres.get(&v) {
                let s = level[&(*edge, v)];
                section.plane_end(order, [s[order[0]], s[order[1]]], *normal, tol)?
            } else if let Some(normal) = leaning {
                let s = level[&(*edge, v)];
                section.plane_end(order, [s[order[0]], s[order[1]]], normal, tol)?
            } else if stands_level(*edge, v) {
                section.arc_at(station[&(*edge, v)], tol)?
            } else {
                let s = level[&(*edge, v)];
                section.skew_end(order, [s[order[0]], s[order[1]]], tol)?
            };
            // The arc runs from the first face of its order's contact.
            let ends = if from_face == 0 { (p0, p1) } else { (p1, p0) };
            edges.push(SheetEdge {
                ends,
                curve: Some((curve, range)),
            });
            arcs.insert((*edge, v), (edges.len() - 1, from_face == 0));
        }
    }

    // Where two rounds meet in a crease beside an edge left sharp: each
    // strip's end at the vertex as a path of edges from its first face's
    // corner to its second's, and the smaller's corner on the far face.
    let mut end_paths: FastMap<(Id, Id), Vec<usize>> = FastMap::default();
    let mut strip_corner: FastMap<(Id, usize, Id), usize> = FastMap::default();
    let mut crease_list: Vec<(&Id, &Crease)> = creases.iter().collect();
    crease_list.sort_unstable_by_key(|(v, _)| **v);
    for (&v, crease) in crease_list {
        let Some(built) = crease_edges(planar, &sections, &edge_seat, crease, v, &level, tol)?
        else {
            return Ok(None);
        };
        points.push(built.q);
        let q = points.len() - 1;
        let [small, large] = crease.pair;
        let c_f = at[&(crease.shared, v)];
        let p_b = at[&(crease.near, v)];
        edges.push(SheetEdge {
            ends: (c_f, q),
            curve: Some(built.crease),
        });
        let crease_edge = edges.len() - 1;
        let (t_q, t_b) = built.arc_ends;
        let (ends, range) = if t_b < t_q {
            ((p_b, q), (t_b, t_q))
        } else {
            ((q, p_b), (t_q, t_b))
        };
        edges.push(SheetEdge {
            ends,
            curve: Some((built.arc, range)),
        });
        let arc_edge = edges.len() - 1;
        end_paths.insert((small, v), vec![crease_edge]);
        end_paths.insert((large, v), vec![crease_edge, arc_edge]);
        strip_corner.insert((small, crease.far, v), q);
        crease_corners.insert(v, (q, p_b, arc_edge));
    }

    // At a ridge each strip's corners are its cap's ball's feet.
    for &(&v, ridge) in &ridge_list {
        for seat in seats_at(&keys, seats, v) {
            if !sections.contains_key(&seat.edge) {
                continue;
            }
            let Some(k) = ridge.ball_of(seat) else {
                return Ok(None);
            };
            for f in seat.faces {
                let Some(&q) = ridge_points.get(&(v, f, k)) else {
                    return Ok(None);
                };
                strip_corner.insert((seat.edge, f, v), q);
            }
        }
    }
    // The faces: each with whether its surface's own normal faces out of
    // the solid along its ring (wound counter-clockwise from outside).
    let mut faces: Vec<(SheetFace, bool)> = Vec::new();
    // What each face of the sheet is built from, in turn: an original
    // face, an edge's strip or a vertex's sphere.
    let mut origins: Vec<Result<usize, Id>> = Vec::new();
    for (fi, face) in planar.faces.iter().enumerate() {
        let n = face.corners.len();
        let mut ring: Vec<(usize, bool)> = Vec::new();
        let mut at_corner: Vec<Corner> = Vec::with_capacity(n);
        for k in 0..n {
            let (prev, v, next) = (
                face.corners[(k + n - 1) % n],
                face.corners[k],
                face.corners[(k + 1) % n],
            );
            if ridges
                .get(&v)
                .is_some_and(|r| fi == r.faces[0] || fi == r.faces[2])
            {
                // Along the ridge's line on this face, from the arriving
                // strip's ball's foot to the other's.
                let Some(pts) = passes(fi, prev, v) else {
                    return Ok(None);
                };
                let line = segment(&mut edges, pts[0], pts[1]);
                at_corner.push((pts, Some(line)));
                continue;
            }
            match flush.get(&v) {
                Some(&(third, edge)) if third == fi => {
                    // In from the face across the arriving edge, out to the
                    // face across the leaving one: each a strip's face.
                    let across = |a: Id, b: Id| -> Option<usize> {
                        let seat = seats.get(&(a.min(b), a.max(b)))?;
                        seat.faces.iter().copied().find(|f| *f != fi)
                    };
                    let (Some(g_in), Some(g_out)) = (across(prev, v), across(v, next)) else {
                        return Ok(None);
                    };
                    let (entry, exit) = (at[&(g_in, v)], at[&(g_out, v)]);
                    let (arc, _) = arcs[&(edge, v)];
                    let runs = if edges[arc].ends == (entry, exit) {
                        true
                    } else if edges[arc].ends == (exit, entry) {
                        false
                    } else {
                        return Ok(None);
                    };
                    at_corner.push((vec![entry, exit], Some((arc, runs))));
                }
                _ => match (creases.get(&v), crease_corners.get(&v)) {
                    (Some(crease), Some(&(q, p_b, arc))) if crease.far == fi => {
                        // In along the smaller round's contact to its corner,
                        // round the larger's end in this face to the sharp
                        // edge; or the other way.
                        let arriving = seats.get(&(prev.min(v), prev.max(v))).map(|s| s.edge);
                        let (entry, exit) = if arriving == Some(crease.pair[0]) {
                            (q, p_b)
                        } else {
                            (p_b, q)
                        };
                        at_corner
                            .push((vec![entry, exit], Some((arc, edges[arc].ends.0 == entry))));
                    }
                    _ => at_corner.push((vec![at[&(fi, v)]], None)),
                },
            }
        }
        for (k, (pts, arc)) in at_corner.iter().enumerate() {
            let next_first = at_corner[(k + 1) % n].0[0];
            if let Some(arc) = arc {
                ring.push(*arc);
            }
            let last = *pts.last().unwrap_or(&next_first);
            if last != next_first {
                ring.push(segment(&mut edges, last, next_first));
            }
        }
        let (surface, out) = match face.plane {
            Some(plane) => {
                let out = plane.normal().vector().dot(face.normal) > 0.0;
                (
                    Some(SurfaceGeometry::from(ogeom_geom::PlaneSurface::new(plane))),
                    out,
                )
            }
            None => (None, true),
        };
        faces.push((SheetFace { surface, ring }, out));
        origins.push(Ok(fi));
    }
    // Each strip's face in the sheet, by its edge.
    let mut strip_face: FastMap<Id, usize> = FastMap::default();
    // The strips: along the first face's contact backward, the arc at the
    // edge's start, along the second face's contact, the arc at its end.
    for edge in &strip_list {
        let (section, _) = &sections[edge];
        let seat = edge_seat[edge];
        let (a, b) = seat.ends;
        let corner = |k: usize, w: Id| {
            strip_corner
                .get(&(*edge, seat.faces[k], w))
                .copied()
                .unwrap_or(at[&(seat.faces[k], w)])
        };
        let path = |w: Id| -> Vec<usize> {
            end_paths
                .get(&(*edge, w))
                .cloned()
                .unwrap_or_else(|| vec![arcs[&(*edge, w)].0])
        };
        let mut ring = Vec::with_capacity(5);
        ring.push(segment(&mut edges, corner(0, b), corner(0, a)));
        let Some(end_a) = walk_path(&edges, &path(a), corner(0, a), corner(1, a)) else {
            return Ok(None);
        };
        ring.extend(end_a);
        ring.push(segment(&mut edges, corner(1, a), corner(1, b)));
        let Some(end_b) = walk_path(&edges, &path(b), corner(1, b), corner(0, b)) else {
            return Ok(None);
        };
        ring.extend(end_b);
        // Over the strip's own stations, from the least of its four corners'
        // and its ends' to the greatest, no farther: a surface reaching past
        // its face's trim is one more thing a later operation reads. An end
        // cut by a plane leaning across the crease reaches past its corners'
        // stations between them.
        // A ridge's strip stands level at its cap's ball.
        let level_at = |w: Id| -> [f64; 2] {
            ridges
                .get(&w)
                .and_then(|r| r.ball_of(seat))
                .map_or(level[&(*edge, w)], |k| {
                    let s = (ridges[&w].balls[k].0 - section.crease(0.0)).dot(seat.along);
                    [s, s]
                })
        };
        let (corners_lo, corners_hi) = (
            level_at(a)[0].min(level_at(a)[1]),
            level_at(b)[0].max(level_at(b)[1]),
        );
        let (mut lo, mut hi) = (corners_lo, corners_hi);
        for end in path(a).into_iter().chain(path(b)) {
            if let Some((curve, (t0, t1))) = &edges[end].curve {
                for k in 0..=END_SAMPLES {
                    let t = t0 + (t1 - t0) * f64::from(k) / f64::from(END_SAMPLES);
                    let p = ogeom_geom::Curve3d::point_at(curve, t, tol)?;
                    let station = (p - section.crease(0.0)).dot(seat.along);
                    lo = lo.min(station);
                    hi = hi.max(station);
                }
            }
        }
        // Between samples an end's station strays less than the sag of a
        // thirty-second of its arc, a thousandth of the radius: where an end
        // reaches past its corners, the span is widened by that much too.
        let margin = 1e-3 * sections[edge].1;
        let (lo, hi) = (
            if lo < corners_lo - reach {
                lo - margin
            } else {
                corners_lo
            },
            if hi > corners_hi + reach {
                hi + margin
            } else {
                corners_hi
            },
        );
        let surface = section.surface((lo, hi), tol)?;
        strip_face.insert(*edge, faces.len());
        faces.push((
            SheetFace {
                surface: Some(surface),
                ring,
            },
            seat.convex,
        ));
        origins.push(Err(*edge));
    }
    // The spheres: the three arcs round the vertex, each walked against
    // its strip, chained end to end.
    for &v in &vertices {
        if method[&v] != RoundVertex::Sphere {
            continue;
        }
        let here: Vec<&Seat> = keys
            .iter()
            .filter(|(a, b)| *a == v || *b == v)
            .map(|k| &seats[k])
            .collect();
        let Some(ring) = hole_at(&here, &arcs, &edges, &at, v) else {
            return Ok(None);
        };
        let seat = here[0];
        let radius = sections[&seat.edge].1;
        let round = faces_round(&planar.faces, v).unwrap_or_default();
        let Some((centre, feet)) = one_ball(planar, &round, v, radius, seat.convex, tol) else {
            return Ok(None);
        };
        let surface = crate::round::sphere::corner_sphere(centre, radius, &feet, tol)?;
        faces.push((
            SheetFace {
                surface: Some(surface),
                ring,
            },
            seat.convex,
        ));
        origins.push(Err(v));
    }

    // Each ridge: a sphere at each ball, through its strips' ends and the
    // ridge's arc there, and the cylinder the ball sweeps between them,
    // bounded by its two arcs and its lines on the faces either side.
    for &(&v, ridge) in &ridge_list {
        let radius = ridge
            .balls
            .iter()
            .map(|(c, feet)| c.distance(feet[0]))
            .fold(0.0, f64::max);
        let (axis_from, axis_to) = (ridge.balls[0].0, ridge.balls[1].0);
        let along = (axis_to - axis_from).normalized(tol)?;
        let [a, _, c, _] = ridge.faces;
        // The ridge's arc at each ball, from its foot on `A` to its foot
        // on `C`, the shorter way round.
        let mut ends = [0usize; 2];
        for k in 0..2 {
            let (centre, _) = ridge.balls[k];
            let (from, to) = (ridge_points[&(v, a, k)], ridge_points[&(v, c, k)]);
            let (x, w) = (points[from] - centre, points[to] - centre);
            let normal = x.cross(w);
            let frame = ogeom_math::Frame::new(
                centre,
                ogeom_math::Direction::new(normal, tol)?,
                ogeom_math::Direction::new(x, tol)?,
                tol,
            )?;
            let y = normal.cross(x);
            let angle = w.dot(y).atan2(w.dot(x) * normal.magnitude());
            let circle = ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(frame, radius, tol)?);
            edges.push(SheetEdge {
                ends: (from, to),
                curve: Some((circle.into(), (0.0, angle))),
            });
            ends[k] = edges.len() - 1;
        }
        let ends_of = |(e, walk): (usize, bool)| {
            let (p, q) = edges[e].ends;
            if walk { (p, q) } else { (q, p) }
        };
        for k in 0..2 {
            let here: Vec<&Seat> = seats_at(&keys, seats, v)
                .into_iter()
                .filter(|s| ridge.ball_of(s) == Some(k))
                .collect();
            // The strips' ends walked against them, as a hole's are, then
            // the ridge's arc back to the start.
            let mut ring: Vec<(usize, bool)> = Vec::with_capacity(3);
            for s in &here {
                let (arc, _) = arcs[&(s.edge, v)];
                let face = if v == s.ends.0 {
                    s.faces[0]
                } else {
                    s.faces[1]
                };
                ring.push((arc, edges[arc].ends.0 != ridge_points[&(v, face, k)]));
            }
            let [first, second] = ring[..] else {
                return Ok(None);
            };
            let ring = if ends_of(first).1 == ends_of(second).0 {
                vec![first, second]
            } else if ends_of(second).1 == ends_of(first).0 {
                vec![second, first]
            } else {
                return Ok(None);
            };
            let (start, tip) = (ends_of(ring[0]).0, ends_of(ring[1]).1);
            let arc = ends[k];
            let walk = if edges[arc].ends == (tip, start) {
                true
            } else if edges[arc].ends == (start, tip) {
                false
            } else {
                return Ok(None);
            };
            let mut ring = ring;
            ring.push((arc, walk));
            let (centre, _) = ridge.balls[k];
            let feet: Vec<Point> = ring.iter().map(|m| points[ends_of(*m).0]).collect();
            let surface = crate::round::sphere::corner_sphere(centre, radius, &feet, tol)?;
            faces.push((
                SheetFace {
                    surface: Some(surface),
                    ring,
                },
                true,
            ));
            origins.push(Err(v));
        }
        // The cylinder: each arc walked against its sphere, the lines on
        // `C` and `A` between.
        let walked: Vec<(usize, bool)> = faces[faces.len() - 2..]
            .iter()
            .map(|(f, _)| {
                f.ring
                    .iter()
                    .find(|(e, _)| ends.contains(e))
                    .map(|(e, w)| (*e, !*w))
            })
            .collect::<Option<Vec<_>>>()
            .unwrap_or_default();
        let [arc0, arc1] = walked[..] else {
            return Ok(None);
        };
        let (e0, e1) = (ends_of(arc0), ends_of(arc1));
        let line0 = segment(&mut edges, e0.1, e1.0);
        let line1 = segment(&mut edges, e1.1, e0.0);
        let frame = ogeom_math::Frame::new(
            axis_from,
            ogeom_math::Direction::new(along, tol)?,
            ogeom_math::Direction::new(points[ridge_points[&(v, a, 0)]] - axis_from, tol)?,
            tol,
        )?;
        let surface: SurfaceGeometry = ogeom_geom::CylinderSurface::new(
            ogeom_math::Cylinder::new(frame, radius, tol)?,
            (0.0, axis_from.distance(axis_to)),
        )?
        .into();
        faces.push((
            SheetFace {
                surface: Some(surface),
                ring: vec![arc0, line0, arc1, line1],
            },
            true,
        ));
        origins.push(Err(v));
    }
    // Built once, each face whose surface faces into the solid built on
    // its ring walked back and turned over.
    let sheet: Vec<SheetFace> = faces
        .iter()
        .map(|(f, out)| {
            if *out {
                f.clone()
            } else {
                SheetFace {
                    surface: f.surface.clone(),
                    ring: f.ring.iter().rev().map(|(e, w)| (*e, !w)).collect(),
                }
            }
        })
        .collect();
    let built = build_faces(model, &points, &edges, &sheet, tol)?;
    let mut built: Vec<_> = built
        .into_iter()
        .zip(&faces)
        .map(|(f, (_, out))| if *out { f } else { f.reversed() })
        .collect();
    // The patches, each through its hole's strip ends as built, tangent
    // to the strips, held to the bands' fit target.
    let fit = crate::marched::band_fit_target(tol);
    for &v in &vertices {
        if method[&v] != RoundVertex::Patch {
            continue;
        }
        let here = seats_at(&keys, seats, v);
        let Some(ring) = hole_at(&here, &arcs, &edges, &at, v) else {
            return Ok(None);
        };
        // Each end's strip, by its arc.
        let owner: FastMap<usize, Id> = here
            .iter()
            .map(|s| (arcs[&(s.edge, v)].0, s.edge))
            .collect();
        let mut sides = Vec::with_capacity(ring.len());
        for (e, walk) in &ring {
            let strip = strip_face[&owner[e]];
            let (a, b) = edges[*e].ends;
            let Some(edge) = edge_between(model, &built[strip], points[a], points[b], reach)?
            else {
                return Ok(None);
            };
            let (Some((curve, range)), Some(support)) =
                (edges[*e].curve.clone(), faces[strip].0.surface.clone())
            else {
                return Ok(None);
            };
            sides.push(crate::round::transfinite::HoleSide {
                edge,
                curve,
                range,
                forward: *walk,
                support,
                pieces: Vec::new(),
            });
        }
        let patch = crate::round::transfinite::corner_patch(model, &sides, tol)?;
        if !(patch.gap <= fit && patch.angle <= fit) {
            ogeom_core::ogeom_bail!(
                Construction,
                "a round vertex's patch misses the bands' fit target: {:.1e} off its strips' \
                 ends, {:.1e} rad off their tangent planes",
                patch.gap,
                patch.angle
            );
        }
        for face in patch.faces {
            built.push(face);
            origins.push(Err(v));
        }
    }
    let shell = ogeom_algo::make_shell(model, &built)?.shape;
    let solid = ogeom_algo::make_solid(model, &[shell])?.shape;
    super::guards::solid_fits(model, &solid, tol)?;
    // P5 on the curved faces: each strip proved within its edge's material,
    // each sphere within its vertex's, to the bands' fit target.
    let fit = crate::marched::band_fit_target(tol);
    for (face, from) in built.iter().zip(&origins) {
        let Err(element) = from else {
            continue;
        };
        let region = match edge_seat.get(element) {
            Some(seat) => edge_material(planar, seat, &sections),
            None => {
                let round = faces_round(&planar.faces, *element).unwrap_or_default();
                corner_material(
                    planar,
                    &sections,
                    &seats_at(&keys, seats, *element),
                    &round,
                    *element,
                )
            }
        };
        if let Err(out) = ogeom_algo::certify_face_within(model, face, &region, fit, tol)? {
            ogeom_core::ogeom_bail!(
                Construction,
                "a round's face stands {:.1e} past the material it blends, at {:?}{}",
                out.past,
                out.at,
                if out.found { "" } else { " (not cleared)" }
            );
        }
    }

    let mut history = History::new();
    for edge in &strip_list {
        history.delete(&planar.shapes[edge]);
    }
    for (face, from) in built.iter().zip(&origins) {
        match from {
            Ok(fi) => history.modify(&planar.faces[*fi].shape, face.clone()),
            Err(element) => history.generate(&planar.shapes[element], face.clone()),
        }
    }
    // Each edge left sharp, to what of it the setbacks leave: the edge of
    // its first face on its line. Each vertex left sharp, to the vertex at
    // its point. A later step finds them so.
    let at_point = |v: &Shape| {
        model
            .node(v)
            .and_then(|n| n.data().as_vertex().map(|d| d.point))
    };
    for key in &keys {
        let seat = &seats[key];
        if sections.contains_key(&seat.edge) {
            continue;
        }
        let (p, q) = (planar.points[&seat.ends.0], planar.points[&seat.ends.1]);
        let on_line = |x: Point| {
            let d = q - p;
            let t = ((x - p).dot(d) / d.dot(d)).clamp(0.0, 1.0);
            (p + d * t).distance(x) <= reach
        };
        for e in explore_unique(model, &built[seat.faces[0]], ShapeType::Edge)? {
            let Some((a, b)) = ogeom_algo::edge_vertices(model, &e)? else {
                continue;
            };
            if let (Some(x), Some(y)) = (at_point(&a), at_point(&b))
                && on_line(x)
                && on_line(y)
                && x.distance(y) > reach
            {
                history.modify(&planar.shapes[&seat.edge], e.oriented(Orientation::Forward));
                break;
            }
        }
    }
    for &v in &vertices {
        if method[&v] != RoundVertex::Sharp {
            continue;
        }
        let point = planar.points[&v];
        let Some(round) = faces_round(&planar.faces, v) else {
            continue;
        };
        if let Some(kept) = explore_unique(model, &built[round[0]], ShapeType::Vertex)?
            .into_iter()
            .find(|w| at_point(w).is_some_and(|x| x.distance(point) <= reach))
        {
            history.modify(&planar.shapes[&v], kept);
        }
    }
    Ok(Some(Built::new(solid, history)))
}

/// The strips' ends at `v`, each walked against its strip, chained end to
/// end into the ring of the hole they bound; `None` where they do not
/// close into one.
fn hole_at(
    here: &[&Seat],
    arcs: &FastMap<(Id, Id), (usize, bool)>,
    edges: &[SheetEdge],
    at: &FastMap<(usize, Id), usize>,
    v: Id,
) -> Option<Vec<(usize, bool)>> {
    // Each strip walks its end at `v` from its first face's contact at the
    // edge's start, from its second's at the edge's end; the hole the
    // other way.
    let mut pending: Vec<(usize, bool)> = here
        .iter()
        .map(|s| {
            let (arc, _) = arcs[&(s.edge, v)];
            let face = if v == s.ends.0 {
                s.faces[0]
            } else {
                s.faces[1]
            };
            (arc, edges[arc].ends.0 != at[&(face, v)])
        })
        .collect();
    let ends = |(e, walk): (usize, bool)| {
        let (a, b) = edges[e].ends;
        if walk { (a, b) } else { (b, a) }
    };
    let first = pending.remove(0);
    let start = ends(first).0;
    let mut tip = ends(first).1;
    let mut ring = vec![first];
    while !pending.is_empty() {
        let k = pending.iter().position(|m| ends(*m).0 == tip)?;
        let next = pending.remove(k);
        tip = ends(next).1;
        ring.push(next);
    }
    (tip == start).then_some(ring)
}

/// The edge of `face` running between `from` and `to`, as a node (its own
/// orientation forward).
fn edge_between(
    model: &Model,
    face: &Shape,
    from: Point,
    to: Point,
    reach: f64,
) -> OgeomResult<Option<Shape>> {
    for edge in explore_unique(model, face, ShapeType::Edge)? {
        let edge = edge.oriented(Orientation::Forward);
        let Some((a, b)) = ogeom_algo::edge_vertices(model, &edge)? else {
            continue;
        };
        let at = |v: &Shape| {
            model
                .node(v)
                .and_then(|n| n.data().as_vertex().map(|d| d.point))
        };
        let (Some(p), Some(q)) = (at(&a), at(&b)) else {
            continue;
        };
        if (p.distance(from) <= reach && q.distance(to) <= reach)
            || (p.distance(to) <= reach && q.distance(from) <= reach)
        {
            return Ok(Some(edge));
        }
    }
    Ok(None)
}

/// The seats of the edges meeting at `v`, in `keys`' order.
fn seats_at<'a>(keys: &[(Id, Id)], seats: &'a FastMap<(Id, Id), Seat>, v: Id) -> Vec<&'a Seat> {
    keys.iter()
        .filter(|(a, b)| *a == v || *b == v)
        .map(|k| &seats[k])
        .collect()
}

/// Face `fi`'s plane as a half-space, at `origin`, moved out by `rise`.
fn face_below(planar: &Planar, fi: usize, origin: Point, rise: f64) -> HalfSpaces {
    let normal = planar.faces[fi].normal;
    HalfSpaces::Below {
        normal,
        at: origin + normal * rise,
    }
}

/// A concave round's fill: the wedge between its two faces as far out as
/// its contact on each stands past the other's plane.
fn fill(planar: &Planar, seat: &Seat, section: &PlanarSection) -> HalfSpaces {
    let origin = planar.points[&seat.ends.0];
    let rise = |k: usize| {
        let other = section.contact(1 - k, 0.0) - section.crease(0.0);
        planar.faces[seat.faces[k]].normal.dot(other)
    };
    HalfSpaces::Every(vec![
        face_below(planar, seat.faces[0], origin, rise(0)),
        face_below(planar, seat.faces[1], origin, rise(1)),
    ])
}

/// The material along `seat`'s edge once rounded: inside both faces at a
/// convex edge; inside either, or the round's fill, at a concave one.
fn edge_material(
    planar: &Planar,
    seat: &Seat,
    sections: &FastMap<Id, (PlanarSection, f64)>,
) -> HalfSpaces {
    let origin = planar.points[&seat.ends.0];
    let sides = vec![
        face_below(planar, seat.faces[0], origin, 0.0),
        face_below(planar, seat.faces[1], origin, 0.0),
    ];
    match sections.get(&seat.edge) {
        Some((section, _)) if !seat.convex => {
            let mut any = sides;
            any.push(fill(planar, seat, section));
            HalfSpaces::Any(any)
        }
        _ if seat.convex => HalfSpaces::Every(sides),
        _ => HalfSpaces::Any(sides),
    }
}

/// The material round a vertex of three faces once rounded, by its edges'
/// convexity: every face's side at a convex corner; with one concave edge,
/// the opposite face's side and either of the edge's (or its fill); with
/// two, the opposite face's side, or both of the convex edge's, or either
/// fill; with three, any face's side or fill.
fn corner_material(
    planar: &Planar,
    sections: &FastMap<Id, (PlanarSection, f64)>,
    here: &[&Seat],
    round: &[usize],
    v: Id,
) -> HalfSpaces {
    let origin = planar.points[&v];
    let side = |fi: usize| face_below(planar, fi, origin, 0.0);
    let concave: Vec<&&Seat> = here.iter().filter(|s| !s.convex).collect();
    let fills: Vec<HalfSpaces> = concave
        .iter()
        .filter_map(|s| {
            sections
                .get(&s.edge)
                .map(|(section, _)| fill(planar, s, section))
        })
        .collect();
    // The face round `v` the edge `seat` does not bound.
    let opposite = |seat: &Seat| round.iter().copied().find(|f| !seat.faces.contains(f));
    if round.len() > 3 && !concave.is_empty() {
        // Four faces or more of mixed convexity: the solid's own material
        // there, or a fill.
        let mut any = vec![local_cone(planar, here, round, v).map_or_else(
            || HalfSpaces::Every(round.iter().map(|f| side(*f)).collect()),
            |cone| cone.half_spaces(&side),
        )];
        any.extend(fills);
        return HalfSpaces::Any(any);
    }
    match (concave.len(), here.iter().find(|s| s.convex)) {
        (0, _) => HalfSpaces::Every(round.iter().map(|f| side(*f)).collect()),
        (1, _) => {
            let seat = concave[0];
            let mut any = vec![side(seat.faces[0]), side(seat.faces[1])];
            any.extend(fills);
            let far = opposite(seat).map_or(HalfSpaces::Any(Vec::new()), side);
            HalfSpaces::Every(vec![far, HalfSpaces::Any(any)])
        }
        (2, Some(convex)) => {
            let mut any = vec![HalfSpaces::Every(vec![
                side(convex.faces[0]),
                side(convex.faces[1]),
            ])];
            if let Some(far) = opposite(convex) {
                any.push(side(far));
            }
            any.extend(fills);
            HalfSpaces::Any(any)
        }
        _ => {
            let mut any: Vec<HalfSpaces> = round.iter().map(|f| side(*f)).collect();
            any.extend(fills);
            HalfSpaces::Any(any)
        }
    }
}

/// The solid's material near a vertex as a formula of its faces' sides.
#[derive(Debug, Clone)]
enum Cone {
    Face(usize),
    Every(Vec<Cone>),
    Any(Vec<Cone>),
}

impl Cone {
    fn holds(&self, planar: &Planar, v: Id, p: Point) -> bool {
        match self {
            Self::Face(f) => (p - planar.points[&v]).dot(planar.faces[*f].normal) <= 0.0,
            Self::Every(all) => all.iter().all(|c| c.holds(planar, v, p)),
            Self::Any(any) => any.iter().any(|c| c.holds(planar, v, p)),
        }
    }

    fn half_spaces(&self, side: &dyn Fn(usize) -> HalfSpaces) -> HalfSpaces {
        match self {
            Self::Face(f) => side(*f),
            Self::Every(all) => {
                HalfSpaces::Every(all.iter().map(|c| c.half_spaces(side)).collect())
            }
            Self::Any(any) => HalfSpaces::Any(any.iter().map(|c| c.half_spaces(side)).collect()),
        }
    }
}

/// The directions the local solid at a vertex is read along: a Fibonacci
/// spiral's points over the sphere.
const CONE_SAMPLES: u32 = 256;

/// The solid's material at `v`, whose faces `round` (in turn round it) meet
/// along the edges `here`, as a formula of their sides (ADR 0001, §13,
/// question 7): around the cycle, faces joined by concave edges unite and
/// the runs intersect, or faces joined by convex edges intersect and the
/// runs unite, whichever agrees with the solid along every direction read
/// (each classified by a ray's crossings with the faces' sectors at `v`).
/// `None` where neither does, or the faces' edges cannot be read.
fn local_cone(planar: &Planar, here: &[&Seat], round: &[usize], v: Id) -> Option<Cone> {
    let n = round.len();
    let origin = planar.points[&v];
    // The edge between each face and the next, its convexity and its
    // direction from `v`.
    let mut joins = Vec::with_capacity(n);
    for k in 0..n {
        let (a, b) = (round[k], round[(k + 1) % n]);
        let seat = here
            .iter()
            .find(|s| s.faces.contains(&a) && s.faces.contains(&b))?;
        let far = if seat.ends.0 == v {
            seat.ends.1
        } else {
            seat.ends.0
        };
        let d = planar.points[&far] - origin;
        joins.push((seat.convex, d / d.magnitude()));
    }
    // Face `k`'s sector: between the edges it shares with the last face and
    // the next.
    let inside = |p: Point| -> bool {
        let rays = [
            Vector::new(0.371, 0.557, 0.743),
            Vector::new(-0.612, 0.288, 0.737),
            Vector::new(0.205, -0.859, 0.469),
        ];
        let mut votes = 0;
        for w in rays {
            let mut crossings = 0;
            for k in 0..n {
                let normal = planar.faces[round[k]].normal;
                let den = normal.dot(w);
                if den.abs() < 1e-12 {
                    continue;
                }
                let t = normal.dot(origin - p) / den;
                if t <= 0.0 {
                    continue;
                }
                let x = p + w * t - origin;
                let (e0, e1) = (joins[(k + n - 1) % n].1, joins[k].1);
                // `x = a e0 + b e1` in the face's plane.
                let (g00, g01, g11) = (e0.dot(e0), e0.dot(e1), e1.dot(e1));
                let (r0, r1) = (x.dot(e0), x.dot(e1));
                let det = g00 * g11 - g01 * g01;
                if det.abs() < 1e-12 {
                    continue;
                }
                let a = (r0 * g11 - r1 * g01) / det;
                let b = (r1 * g00 - r0 * g01) / det;
                if a >= 0.0 && b >= 0.0 {
                    crossings += 1;
                }
            }
            votes += crossings % 2;
        }
        votes >= 2
    };
    // The runs of faces joined by edges of convexity `joined`, in turn round
    // the cycle from just after an edge of the other.
    let runs = |joined: bool| -> Vec<Vec<usize>> {
        let Some(cut) = joins.iter().position(|(c, _)| *c != joined) else {
            return vec![round.to_vec()];
        };
        let mut out: Vec<Vec<usize>> = Vec::new();
        let mut run = Vec::new();
        for step in 1..=n {
            let k = (cut + step) % n;
            run.push(round[k]);
            if joins[k].0 != joined {
                out.push(core::mem::take(&mut run));
            }
        }
        out
    };
    let faces = |run: &[usize]| run.iter().map(|f| Cone::Face(*f)).collect::<Vec<_>>();
    let candidates = [
        Cone::Every(runs(false).iter().map(|r| Cone::Any(faces(r))).collect()),
        Cone::Any(runs(true).iter().map(|r| Cone::Every(faces(r))).collect()),
    ];
    let scale = joins.iter().map(|(_, d)| d.magnitude()).fold(0.0, f64::max) * 1e-3;
    let golden = core::f64::consts::PI * (3.0 - 5.0f64.sqrt());
    candidates.into_iter().find(|cone| {
        (0..CONE_SAMPLES).all(|i| {
            let z = 1.0 - 2.0 * (f64::from(i) + 0.5) / f64::from(CONE_SAMPLES);
            let rho = (1.0 - z * z).sqrt();
            let t = golden * f64::from(i);
            let p = origin + Vector::new(rho * t.cos(), rho * t.sin(), z) * scale;
            cone.holds(planar, v, p) == inside(p)
        })
    })
}

/// The points a face's boundary passes at one of its corners, in turn, and
/// the arc between them where a strip ends flush in the face.
type Corner = (Vec<usize>, Option<(usize, bool)>);

/// Two rounds of one sense meeting beside an edge left sharp other than in
/// a mitre (radii differing, or axes not meeting): the round whose corner on
/// the sharp edge stands nearer the vertex first, the face the two share,
/// the larger's other face (`near`, where its corner is) and the smaller's
/// (`far`, where the larger runs on to end flush).
#[derive(Debug, Clone, Copy)]
struct Crease {
    pair: [Id; 2],
    shared: usize,
    near: usize,
    far: usize,
}

/// The crease vertex's reading, or `None` where it is not one this engine
/// builds.
fn crease_at(
    sections: &FastMap<Id, (PlanarSection, f64)>,
    here: &[&Seat],
    round: &[usize],
    corners: &FastMap<(usize, Id), Point>,
    v: Id,
    planar: &Planar,
    reach: f64,
) -> Option<Crease> {
    let rounds: Vec<&&Seat> = here
        .iter()
        .filter(|s| sections.contains_key(&s.edge))
        .collect();
    let sharp = here.iter().find(|s| !sections.contains_key(&s.edge))?;
    let [a, b] = rounds.as_slice() else {
        return None;
    };
    if a.convex != b.convex || round.len() != 3 {
        return None;
    }
    let shared = a.faces.iter().copied().find(|f| b.faces.contains(f))?;
    let other = |s: &Seat| s.faces.iter().copied().find(|f| *f != shared);
    let (fa, fb) = (other(a)?, other(b)?);
    if !(sharp.faces.contains(&fa) && sharp.faces.contains(&fb)) {
        return None;
    }
    let at_v = planar.points[&v];
    let (da, db) = (
        corners[&(fa, v)].distance(at_v),
        corners[&(fb, v)].distance(at_v),
    );
    if (da - db).abs() <= reach {
        return None;
    }
    Some(if da < db {
        Crease {
            pair: [a.edge, b.edge],
            shared,
            near: fb,
            far: fa,
        }
    } else {
        Crease {
            pair: [b.edge, a.edge],
            shared,
            near: fa,
            far: fb,
        }
    })
}

/// A crease vertex's curves: where the smaller strip's contact on the far
/// face meets the larger strip (`q`), the crease where the two strips cross
/// from the shared face's corner to `q` (fitted through points solved on
/// both cylinders exactly, with its range), and the larger strip's end in
/// the far face's plane with the parameters of `q` and of its corner on
/// the sharp edge.
struct CreaseEdges {
    q: Point,
    crease: (Curve, (f64, f64)),
    arc: Curve,
    arc_ends: (f64, f64),
}

/// How many points the crease is solved at and fitted through.
const CREASE_POINTS: u32 = 64;

#[allow(
    clippy::too_many_arguments,
    reason = "one vertex's reading of the solid"
)]
fn crease_edges(
    planar: &Planar,
    sections: &FastMap<Id, (PlanarSection, f64)>,
    edge_seat: &FastMap<Id, &Seat>,
    crease: &Crease,
    v: Id,
    level: &FastMap<(Id, Id), [f64; 2]>,
    tol: Tolerances,
) -> OgeomResult<Option<CreaseEdges>> {
    let [small, large] = crease.pair;
    let (seat_s, seat_l) = (edge_seat[&small], edge_seat[&large]);
    let ((sec_s, r_s), (sec_l, r_l)) = (&sections[&small], &sections[&large]);
    // Each strip's axis: a point on it and its unit direction.
    // The ball stands inside the material at a convex edge and in the open
    // corner at a concave one.
    let axis = |section: &PlanarSection, seat: &Seat, r: f64| -> (Point, Vector) {
        let n = planar.faces[seat.faces[0]].normal;
        let side = if seat.convex { -r } else { r };
        (section.contact(0, 0.0) + n * side, seat.along)
    };
    let ((a_s, d_s), (a_l, d_l)) = (axis(sec_s, seat_s, *r_s), axis(sec_l, seat_l, *r_l));
    // The stations `t` along `p + d t` where it stands `r_l` off the larger
    // axis: a quadratic.
    let on_large = |p: Point, d: Vector| -> Vec<f64> {
        let w = p - a_l;
        let w = w - d_l * w.dot(d_l);
        let e = d - d_l * d.dot(d_l);
        let (qa, qb, qc) = (e.dot(e), 2.0 * w.dot(e), w.dot(w) - r_l * r_l);
        if qa <= f64::EPSILON {
            return Vec::new();
        }
        let disc = qb * qb - 4.0 * qa * qc;
        if disc < -1e-12 * qb.abs().max(1.0) {
            return Vec::new();
        }
        let root = disc.max(0.0).sqrt();
        vec![(-qb - root) / (2.0 * qa), (-qb + root) / (2.0 * qa)]
    };
    let at_v = planar.points[&v];
    // `q`: the smaller's contact on the far face, nearest the vertex.
    let k_far = usize::from(seat_s.faces[1] == crease.far);
    let c0 = sec_s.contact(k_far, 0.0);
    let Some(q) = on_large(c0, d_s)
        .into_iter()
        .map(|t| c0 + d_s * t)
        .min_by(|x, y| x.distance(at_v).total_cmp(&y.distance(at_v)))
    else {
        return Ok(None);
    };
    let k_shared = usize::from(seat_s.faces[1] == crease.shared);
    let c_f = {
        let s = level[&(small, v)][k_shared];
        sec_s.contact(k_shared, s)
    };
    // The smaller's section angle from the shared face's corner to `q`.
    let radial = |p: Point| {
        let w = p - a_s;
        w - d_s * w.dot(d_s)
    };
    let e1 = radial(c_f).normalized(tol)?;
    let rq = radial(q);
    let e2 = (rq - e1 * rq.dot(e1)).normalized(tol)?;
    let sweep = rq.dot(e2).atan2(rq.dot(e1));
    // Solved from `q` back to the corner, each root the one nearest the
    // last: the two strips touch at the corner, where both roots meet.
    let mut solved: Vec<Point> = Vec::with_capacity(CREASE_POINTS as usize + 1);
    let mut last = (q - a_s).dot(d_s);
    for k in (0..=CREASE_POINTS).rev() {
        let u =
            0.5 * (1.0 - (core::f64::consts::PI * f64::from(k) / f64::from(CREASE_POINTS)).cos());
        let theta = sweep * u;
        let base = a_s + (e1 * theta.cos() + e2 * theta.sin()) * *r_s;
        let Some(t) = on_large(base, d_s)
            .into_iter()
            .min_by(|x, y| (x - last).abs().total_cmp(&(y - last).abs()))
        else {
            return Ok(None);
        };
        last = t;
        solved.push(base + d_s * t);
    }
    solved.reverse();
    if solved[0].distance(c_f) > 1e3 * tol.confusion() {
        return Ok(None);
    }
    solved[0] = c_f;
    if let Some(end) = solved.last_mut() {
        *end = q;
    }
    let fitted = ogeom_algo::interpolate(&solved, 3, ogeom_algo::Spacing::default(), tol)?;
    let (t0, t1) = ogeom_geom::Curve3d::domain(&fitted);
    // Held to both cylinders between the points it was fitted through.
    let off = |p: Point, a: Point, d: Vector, r: f64| {
        let w = p - a;
        ((w - d * w.dot(d)).magnitude() - r).abs()
    };
    for k in 0..(4 * CREASE_POINTS) {
        let t = t0 + (t1 - t0) * (f64::from(k) + 0.5) / f64::from(4 * CREASE_POINTS);
        let p = ogeom_geom::Curve3d::point_at(&fitted, t, tol)?;
        if off(p, a_s, d_s, *r_s).max(off(p, a_l, d_l, *r_l)) > tol.confusion() {
            return Ok(None);
        }
    }
    // The larger's end in the far face's plane, from its corner on the sharp
    // edge (its contact on the near face) to its contact on the shared face.
    let normal = planar.faces[crease.far].normal;
    let k_near = usize::from(seat_l.faces[1] == crease.near);
    let k_shared_l = 1 - k_near;
    let mut stations = [0.0; 2];
    stations[k_near] = level[&(large, v)][k_near];
    let corner_b = sec_l.contact(k_near, stations[k_near]);
    let c_shared = sec_l.contact(k_shared_l, 0.0);
    let along_l = seat_l.along;
    if normal.dot(along_l).abs() <= tol.angular() {
        return Ok(None);
    }
    stations[k_shared_l] = normal.dot(corner_b - c_shared) / normal.dot(along_l);
    let order = sec_l.order();
    let (arc, _) = sec_l.plane_end(order, [stations[order[0]], stations[order[1]]], normal, tol)?;
    let t_b = if order[0] == k_near { 0.0 } else { 1.0 };
    let t_q = ogeom_algo::project_on_curve(&arc, q, 64, tol)?.parameter;
    if ogeom_geom::Curve3d::point_at(&arc, t_q, tol)?.distance(q) > 1e3 * tol.confusion()
        || !within_face(planar, crease.far, q, 1e3 * tol.confusion())
    {
        return Ok(None);
    }
    Ok(Some(CreaseEdges {
        q,
        crease: (Curve::BSpline(fitted), (t0, t1)),
        arc,
        arc_ends: (t_q, t_b),
    }))
}

/// `path`'s edges walked from `from` to `to`, in its order or against it.
fn walk_path(
    edges: &[SheetEdge],
    path: &[usize],
    from: usize,
    to: usize,
) -> Option<Vec<(usize, bool)>> {
    let walk = |order: &mut dyn Iterator<Item = &usize>| -> Option<Vec<(usize, bool)>> {
        let mut at = from;
        let mut out = Vec::new();
        for &e in order {
            let (a, b) = edges[e].ends;
            if a == at {
                out.push((e, true));
                at = b;
            } else if b == at {
                out.push((e, false));
                at = a;
            } else {
                return None;
            }
        }
        (at == to).then_some(out)
    };
    walk(&mut path.iter()).or_else(|| walk(&mut path.iter().rev()))
}

/// A vertex of three edges where two rounds of one radius and one sense
/// meet beside an edge left sharp, both rounds' setback corners on the
/// sharp edge one point and their axes meeting: the two rounds in their
/// edges' order, and the normal of the mitre plane they meet in (through
/// the shared face's setback corner, the sharp edge's and where the axes
/// meet, which holds the two cylinders' crossing). `None` for any other.
#[allow(
    clippy::too_many_arguments,
    reason = "one vertex's reading of the solid"
)]
fn mitre(
    planar: &Planar,
    sections: &FastMap<Id, (PlanarSection, f64)>,
    here: &[&Seat],
    round: &[usize],
    corners: &FastMap<(usize, Id), Point>,
    v: Id,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Option<([Id; 2], Vector)>> {
    let mut rounds: Vec<&&Seat> = here
        .iter()
        .filter(|s| sections.contains_key(&s.edge))
        .collect();
    rounds.sort_by_key(|s| s.edge);
    let (Some(sharp), [a, b]) = (
        here.iter().find(|s| !sections.contains_key(&s.edge)),
        rounds.as_slice(),
    ) else {
        return Ok(None);
    };
    let (ra, rb) = (sections[&a.edge].1, sections[&b.edge].1);
    if ra != rb || a.convex != b.convex {
        return Ok(None);
    }
    let Some(shared) = a.faces.iter().copied().find(|f| b.faces.contains(f)) else {
        return Ok(None);
    };
    let other = |s: &Seat| s.faces.iter().copied().find(|f| *f != shared);
    let (Some(fa), Some(fb)) = (other(a), other(b)) else {
        return Ok(None);
    };
    if !(sharp.faces.contains(&fa) && sharp.faces.contains(&fb)) || round.len() != 3 {
        return Ok(None);
    }
    let on_sharp = corners[&(fa, v)];
    if on_sharp.distance(corners[&(fb, v)]) > reach {
        return Ok(None);
    }
    // Each round's axis: its ball's centre along its crease.
    let axis = |seat: &Seat| -> (Point, Vector) {
        let (section, r) = &sections[&seat.edge];
        let n = planar.faces[seat.faces[0]].normal;
        let side = if seat.convex { -*r } else { *r };
        (section.contact(0, 0.0) + n * side, seat.along)
    };
    let ((pa, da), (pb, db)) = (axis(a), axis(b));
    // Where the two axes come closest, and how close.
    let w = pa - pb;
    let (aa, ab, bb, aw, bw) = (da.dot(da), da.dot(db), db.dot(db), da.dot(w), db.dot(w));
    let det = aa * bb - ab * ab;
    if det.abs() <= tol.angular() {
        return Ok(None);
    }
    let (t, u) = ((ab * bw - bb * aw) / det, (aa * bw - ab * aw) / det);
    let (qa, qb) = (pa + da * t, pb + db * u);
    if qa.distance(qb) > reach {
        return Ok(None);
    }
    let meet = Point::from_vector((qa.to_vector() + qb.to_vector()) * 0.5);
    let on_shared = corners[&(shared, v)];
    let normal = (on_sharp - on_shared).cross(meet - on_shared);
    if normal.magnitude() <= tol.confusion() {
        return Ok(None);
    }
    Ok(Some(([a.edge, b.edge], normal)))
}

/// The points read along a strip's end for the stations it reaches.
const END_SAMPLES: u32 = 32;

/// The length of a seat's edge.
/// A convex vertex of four planes rounded at one radius on every edge
/// where no one ball touches all four faces but one touches three at each
/// end of a ridge: the faces either side of the ridge and the two caps, in
/// the vertex's round order `A, B, C, D`, and each cap's ball, its centre
/// and its feet on `A`, the cap and `C`.
struct Ridge {
    faces: [usize; 4],
    balls: [(Point, [Point; 3]); 2],
}

impl Ridge {
    /// Which ball an edge's strip ends on: the one of the cap it bounds.
    fn ball_of(&self, seat: &Seat) -> Option<usize> {
        if seat.faces.contains(&self.faces[1]) {
            Some(0)
        } else if seat.faces.contains(&self.faces[3]) {
            Some(1)
        } else {
            None
        }
    }

    /// Ball `k`'s foot on face `f`, where it touches `f`.
    fn foot(&self, k: usize, f: usize) -> Option<Point> {
        let (_, feet) = self.balls[k];
        if f == self.faces[0] {
            Some(feet[0])
        } else if f == self.faces[1 + 2 * k] {
            Some(feet[1])
        } else if f == self.faces[2] {
            Some(feet[2])
        } else {
            None
        }
    }
}

/// The [`Ridge`] at `v`, whose faces `round` are four, rounded convex at
/// `radius`: of the two ways to pair the faces opposite, the one where the
/// ball touching a pair and the cap between them stands clear of the other
/// cap at each end, apart from the other ball, its feet on the caps their
/// setback corners and on the pair within their faces. `None` where
/// neither pairing does.
fn ridge_at(
    planar: &Planar,
    round: &[usize],
    v: Id,
    radius: f64,
    corners: &FastMap<(usize, Id), Point>,
    reach: f64,
    tol: Tolerances,
) -> Option<Ridge> {
    let [r0, r1, r2, r3] = round[..] else {
        return None;
    };
    let plane = |f: usize| (planar.faces[f].normal, planar.points[&v]);
    let inside = |centre: Point, f: usize| {
        (centre - planar.points[&v]).dot(planar.faces[f].normal) <= -radius + reach
    };
    for faces in [[r0, r1, r2, r3], [r1, r2, r3, r0]] {
        let [a, b, c, d] = faces;
        let ball = |cap: usize| {
            crate::round::sphere::corner_ball([plane(a), plane(cap), plane(c)], radius, true, tol)
        };
        let (Some((at_b, on_b)), Some((at_d, on_d))) = (ball(b), ball(d)) else {
            continue;
        };
        if at_b.distance(at_d) <= reach || !inside(at_b, d) || !inside(at_d, b) {
            continue;
        }
        let (Some(corner_b), Some(corner_d)) = (corners.get(&(b, v)), corners.get(&(d, v))) else {
            continue;
        };
        if on_b[1].distance(*corner_b) > reach || on_d[1].distance(*corner_d) > reach {
            continue;
        }
        let within = [on_b[0], on_d[0]]
            .iter()
            .all(|p| within_face(planar, a, *p, reach))
            && [on_b[2], on_d[2]]
                .iter()
                .all(|p| within_face(planar, c, *p, reach));
        if within {
            return Some(Ridge {
                faces,
                balls: [(at_b, on_b), (at_d, on_d)],
            });
        }
    }
    None
}

/// The ball of `radius` touching the planes of the faces `round` vertex
/// `v` (inside a convex vertex's material, in a concave one's open corner),
/// found from the first three: its centre and its foot on each face's
/// plane, in turn. `None` where the first three planes meet in no single
/// point; a foot off its plane where the ball misses that face.
fn one_ball(
    planar: &Planar,
    round: &[usize],
    v: Id,
    radius: f64,
    convex: bool,
    tol: Tolerances,
) -> Option<(Point, Vec<Point>)> {
    let plane = |f: usize| (planar.faces[f].normal, planar.points[&v]);
    let three = [
        plane(*round.first()?),
        plane(*round.get(1)?),
        plane(*round.get(2)?),
    ];
    let (centre, _) = crate::round::sphere::corner_ball(three, radius, convex, tol)?;
    let off = if convex { -radius } else { radius };
    let feet = round
        .iter()
        .map(|f| centre - planar.faces[*f].normal * off)
        .collect();
    Some((centre, feet))
}

fn section_length(planar: &Planar, seat: &Seat) -> f64 {
    planar.points[&seat.ends.0].distance(planar.points[&seat.ends.1])
}

/// Whether `p`, a point of face `fi`'s plane, lies within the face's
/// polygon or within `reach` of its boundary.
fn within_face(planar: &Planar, fi: usize, p: Point, reach: f64) -> bool {
    let face = &planar.faces[fi];
    let ring: Vec<Point> = face.corners.iter().map(|c| planar.points[c]).collect();
    let n = face.normal;
    let x = ring
        .iter()
        .map(|q| *q - ring[0])
        .find(|d| d.magnitude() > reach)
        .map_or(Vector::ZERO, |d| d / d.magnitude());
    let y = n.cross(x);
    let flat = |q: Point| ((q - p).dot(x), (q - p).dot(y));
    let mut winding = 0.0;
    for k in 0..ring.len() {
        let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
        // On the boundary, within reach.
        let side = b - a;
        let t = ((p - a).dot(side) / side.dot(side)).clamp(0.0, 1.0);
        if (a + side * t).distance(p) <= reach {
            return true;
        }
        let ((ax, ay), (bx, by)) = (flat(a), flat(b));
        winding += (ax * by - ay * bx).atan2(ax * bx + ay * by);
    }
    winding.abs() > core::f64::consts::PI
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::super::{Request, run};
    use ogeom_core::Tolerances;
    use ogeom_math::Frame;
    use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    fn volume(m: &Model, s: &Shape) -> f64 {
        let fine = ogeom_mesh::Deflection {
            chord: 1e-5,
            ..ogeom_mesh::Deflection::default()
        };
        ogeom_algo::volume_properties(m, s, fine, T).unwrap().mass
    }

    fn rounded(edges: impl Fn(&[Shape]) -> Vec<Shape>, r: f64) -> (Model, Option<Shape>) {
        let mut m = Model::new();
        let s = ogeom_algo::make_box(&mut m, Frame::WORLD, (10.0, 8.0, 6.0), T)
            .unwrap()
            .shape;
        let all = explore_unique(&m, &s, ShapeType::Edge).unwrap();
        let asked: Vec<(Shape, f64)> = edges(&all).into_iter().map(|e| (e, r)).collect();
        let request = Request::from_rounds(&asked, T).unwrap();
        let built = run(&mut m, &s, &request, T).unwrap().map(|b| b.shape);
        if let Some(b) = &built {
            let found = ogeom_algo::check(&m, b, T).unwrap();
            assert!(found.is_valid(), "{:?}", found.problems);
        }
        (m, built)
    }

    /// A box with one edge rounded: two flush ends, the box less the
    /// round's section `(1 - pi/4) r^2` along the edge.
    #[test]
    fn one_round_ends_flush_in_the_box() {
        let (m, built) = rounded(|all| vec![all[0].clone()], 1.5);
        let built = built.unwrap();
        assert_eq!(
            explore_unique(&m, &built, ShapeType::Edge).unwrap().len(),
            15
        );
        let v = volume(&m, &built);
        let lost = (1.0 - core::f64::consts::FRAC_PI_4) * 1.5 * 1.5;
        let edge_lengths = [10.0, 8.0, 6.0];
        let want: Vec<f64> = edge_lengths.iter().map(|l| 480.0 - lost * l).collect();
        assert!(
            want.iter().any(|w| (v - w).abs() < 1e-6),
            "{v} not one of {want:?}"
        );
    }

    /// An L-shaped prism, its profile in the plane x = 0 swept 10 along
    /// +x and stated in place, with its six sweep edges rounded: flush ends
    /// on both caps, the concave edge's strip filling what the convex five
    /// take off, the profile's area less `(1 - pi/4) r^2` per convex corner
    /// and more per concave one, along the 10.
    #[test]
    fn sweep_rounds_end_flush_on_the_caps_convex_or_concave() {
        use ogeom_math::{Direction, Plane, Point, Vector};
        let mut m = Model::new();
        let profile = [
            (0.0, 0.0),
            (6.0, 0.0),
            (6.0, 2.0),
            (2.0, 2.0),
            (2.0, 5.0),
            (0.0, 5.0),
        ];
        let points: Vec<Point> = profile
            .iter()
            .map(|&(y, z)| Point::new(0.0, y, z))
            .collect();
        let wire = ogeom_algo::make_polygon(&mut m, &points, true, T)
            .unwrap()
            .shape;
        let plane = ogeom_geom::PlaneSurface::over(
            Plane::through(Point::ORIGIN, Direction::X),
            (-100.0, 100.0),
            (-100.0, 100.0),
        )
        .unwrap();
        let face = ogeom_algo::make_face(&mut m, plane.into(), &[wire], T)
            .unwrap()
            .shape;
        let swept = ogeom_algo::make_prism(&mut m, &face, Vector::new(10.0, 0.0, 0.0), T)
            .unwrap()
            .shape;
        let s = ogeom_algo::baked_shape(&mut m, &swept, T).unwrap().shape;
        let r = 0.5;
        let along: Vec<(Shape, f64)> = explore_unique(&m, &s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
                let at = |v: &Shape| {
                    m.node(v)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .unwrap()
                };
                ((at(&a).x - at(&b).x).abs() - 10.0).abs() < 1e-9
            })
            .map(|e| (e, r))
            .collect();
        assert_eq!(along.len(), 6);
        let request = Request::from_rounds(&along, T).unwrap();
        let built = run(&mut m, &s, &request, T).unwrap().unwrap().shape;
        let corner = (1.0 - core::f64::consts::FRAC_PI_4) * r * r;
        let want = (18.0 - 5.0 * corner + corner) * 10.0;
        let v = volume(&m, &built);
        assert!((v - want).abs() < 1e-6, "{v} against {want}");
    }

    /// A box rounded at half its least side has its two smallest faces set
    /// back to lines: not a class the sheet builds, so the request is left
    /// to the current path; a hair less builds the exact rounded box.
    #[test]
    fn a_face_set_back_to_a_line_is_left_to_the_current_path() {
        let (_, built) = rounded(<[Shape]>::to_vec, 3.0);
        assert!(built.is_none());
        let r = 2.9;
        let (m, built) = rounded(<[Shape]>::to_vec, r);
        let (a, b, c) = (10.0 - 2.0 * r, 8.0 - 2.0 * r, 6.0 - 2.0 * r);
        let pi = core::f64::consts::PI;
        let want = a * b * c
            + 2.0 * r * (a * b + b * c + c * a)
            + pi * r * r * (a + b + c)
            + 4.0 / 3.0 * pi * r * r * r;
        let v = volume(&m, &built.unwrap());
        assert!((v - want).abs() < 1e-6, "{v} against {want}");
    }

    /// The request's order is not the solid's: every edge of the box asked
    /// backwards builds the same points, bit for bit, and the same volume.
    #[test]
    fn the_request_s_order_does_not_change_the_rounded_box() {
        let bits = |m: &Model, s: &Shape| {
            let mut out: Vec<[u64; 3]> = explore_unique(m, s, ShapeType::Vertex)
                .unwrap()
                .iter()
                .map(|v| {
                    let p = m
                        .node(v)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .unwrap();
                    [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
                })
                .collect();
            out.sort_unstable();
            (volume(m, s).to_bits(), out)
        };
        let (m, forward) = rounded(<[Shape]>::to_vec, 1.5);
        let (n, backward) = rounded(|all| all.iter().rev().cloned().collect(), 1.5);
        assert_eq!(bits(&m, &forward.unwrap()), bits(&n, &backward.unwrap()));
    }

    /// The faces of `s` on a surface of `kind`, counted.
    fn count_of(m: &Model, s: &Shape, kind: ogeom_geom::SurfaceKind) -> usize {
        use ogeom_geom::Surface as _;
        explore_unique(m, s, ShapeType::Face)
            .unwrap()
            .iter()
            .filter(|f| {
                let d = m.node(f).unwrap().data().as_face().unwrap();
                m.geometry().surface(d.surface).unwrap().kind() == kind
            })
            .count()
    }

    /// The vertices of `s` as bits, sorted, and its volume's bits.
    fn bits(m: &Model, s: &Shape) -> (u64, Vec<[u64; 3]>) {
        let mut out: Vec<[u64; 3]> = explore_unique(m, s, ShapeType::Vertex)
            .unwrap()
            .iter()
            .map(|v| {
                let p = m
                    .node(v)
                    .and_then(|n| n.data().as_vertex().map(|d| d.point))
                    .unwrap();
                [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
            })
            .collect();
        out.sort_unstable();
        (volume(m, s).to_bits(), out)
    }

    /// Every edge of `s` rounded at `radius(edge)`, in the order given and
    /// reversed: built by the engine both ways, valid, the same to the bit.
    fn rounded_both_ways(
        m: &mut Model,
        s: &Shape,
        radius: impl Fn(&Model, &Shape) -> f64,
    ) -> Shape {
        let all = explore_unique(m, s, ShapeType::Edge).unwrap();
        let asked: Vec<(Shape, f64)> = all.iter().map(|e| (e.clone(), radius(m, e))).collect();
        let mut built = Vec::new();
        for order in [asked.clone(), asked.iter().rev().cloned().collect()] {
            let request = Request::from_rounds(&order, T).unwrap();
            let b = run(m, s, &request, T).unwrap().unwrap().shape;
            let found = ogeom_algo::check(m, &b, T).unwrap();
            assert!(found.is_valid(), "{:?}", found.problems);
            built.push(b);
        }
        assert_eq!(bits(m, &built[0]), bits(m, &built[1]));
        built.remove(0)
    }

    /// A 20 by 14 by 9 box, its edges along `x`, `y` and `z` rounded at 1,
    /// 1.5 and 2: no ball fits any corner, so each is a patch through the
    /// three strips' ends, skew where two radii differ on a face. The
    /// solid is valid, its patches pass the engine's measures and P5, and
    /// it holds less than the box less the strips short of the corners and
    /// more than the box less the whole strips and a corner cell each.
    #[test]
    fn a_box_of_three_radii_closes_each_corner_with_a_patch() {
        let (a, b, c) = (20.0, 14.0, 9.0);
        let mut m = Model::new();
        let s = ogeom_algo::make_box(&mut m, Frame::WORLD, (a, b, c), T)
            .unwrap()
            .shape;
        let radius = |m: &Model, e: &Shape| {
            let (p, q) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
            let at = |v: &Shape| {
                m.node(v)
                    .and_then(|n| n.data().as_vertex().map(|d| d.point))
                    .unwrap()
            };
            let d = at(&q) - at(&p);
            if d.x.abs() > 1e-9 {
                1.0
            } else if d.y.abs() > 1e-9 {
                1.5
            } else {
                2.0
            }
        };
        let built = rounded_both_ways(&mut m, &s, radius);
        assert_eq!(count_of(&m, &built, ogeom_geom::SurfaceKind::BSpline), 8);
        assert_eq!(count_of(&m, &built, ogeom_geom::SurfaceKind::Cylinder), 12);
        let corner = |r: f64| (1.0 - core::f64::consts::FRAC_PI_4) * r * r;
        let edges = [(1.0, a), (1.5, b), (2.0, c)];
        let whole: f64 = edges.iter().map(|(r, l)| 4.0 * corner(*r) * l).sum();
        let short: f64 = edges
            .iter()
            .map(|(r, l)| 4.0 * corner(*r) * (l - 4.0))
            .sum();
        let v = volume(&m, &built);
        assert!(
            v < a * b * c - short && v > a * b * c - whole - 8.0 * 8.0,
            "{v}"
        );
    }

    /// An L-shaped prism with every edge rounded at 0.5: its two vertices
    /// at the concave edge's ends meet two convex rounds and one concave,
    /// where no ball fits, and close with patches across the caps' reflex
    /// corners; the rest are spheres. Valid, and the same to the bit in
    /// either order.
    #[test]
    fn an_l_prism_s_mixed_corners_close_with_patches() {
        use ogeom_math::{Direction, Plane, Point, Vector};
        let mut m = Model::new();
        let profile = [
            (0.0, 0.0),
            (6.0, 0.0),
            (6.0, 2.0),
            (2.0, 2.0),
            (2.0, 5.0),
            (0.0, 5.0),
        ];
        let points: Vec<Point> = profile
            .iter()
            .map(|&(y, z)| Point::new(0.0, y, z))
            .collect();
        let wire = ogeom_algo::make_polygon(&mut m, &points, true, T)
            .unwrap()
            .shape;
        let plane = ogeom_geom::PlaneSurface::over(
            Plane::through(Point::ORIGIN, Direction::X),
            (-100.0, 100.0),
            (-100.0, 100.0),
        )
        .unwrap();
        let face = ogeom_algo::make_face(&mut m, plane.into(), &[wire], T)
            .unwrap()
            .shape;
        let swept = ogeom_algo::make_prism(&mut m, &face, Vector::new(10.0, 0.0, 0.0), T)
            .unwrap()
            .shape;
        let s = ogeom_algo::baked_shape(&mut m, &swept, T).unwrap().shape;
        let built = rounded_both_ways(&mut m, &s, |_, _| 0.5);
        assert_eq!(count_of(&m, &built, ogeom_geom::SurfaceKind::BSpline), 2);
        assert_eq!(count_of(&m, &built, ogeom_geom::SurfaceKind::Sphere), 10);
    }

    /// The corpus's L channel (walls 5 thick, swept 30) with every edge
    /// rounded at 1: at each end of its concave edge the two convex strips'
    /// skew ends meet on the cap at its setback corner running straight on,
    /// a corner of the hole where the patch's chart has no normal; the
    /// patch is still tangent to every strip within the fit target, and
    /// the solid valid, the same in either order.
    #[test]
    fn an_l_channel_s_corners_patch_where_the_hole_runs_straight_on() {
        use ogeom_math::{Direction, Plane, Point, Vector};
        let mut m = Model::new();
        let profile = [
            (0.0, 0.0),
            (20.0, 0.0),
            (20.0, 5.0),
            (5.0, 5.0),
            (5.0, 20.0),
            (0.0, 20.0),
        ];
        let points: Vec<Point> = profile
            .iter()
            .map(|&(y, z)| Point::new(0.0, y, z))
            .collect();
        let wire = ogeom_algo::make_polygon(&mut m, &points, true, T)
            .unwrap()
            .shape;
        let plane = ogeom_geom::PlaneSurface::over(
            Plane::through(Point::ORIGIN, Direction::X),
            (-100.0, 100.0),
            (-100.0, 100.0),
        )
        .unwrap();
        let face = ogeom_algo::make_face(&mut m, plane.into(), &[wire], T)
            .unwrap()
            .shape;
        let swept = ogeom_algo::make_prism(&mut m, &face, Vector::new(30.0, 0.0, 0.0), T)
            .unwrap()
            .shape;
        let s = ogeom_algo::baked_shape(&mut m, &swept, T).unwrap().shape;
        let built = rounded_both_ways(&mut m, &s, |_, _| 1.0);
        assert_eq!(count_of(&m, &built, ogeom_geom::SurfaceKind::BSpline), 2);
    }

    /// A box with one edge rounded: every edge left sharp comes through the
    /// history as one edge of the result between two of its faces, on its
    /// own line; every vertex left sharp as a vertex at its own point.
    #[test]
    fn edges_and_vertices_left_sharp_come_through_the_history() {
        let mut m = Model::new();
        let s = ogeom_algo::make_box(&mut m, Frame::WORLD, (10.0, 8.0, 6.0), T)
            .unwrap()
            .shape;
        let all = explore_unique(&m, &s, ShapeType::Edge).unwrap();
        let request = Request::from_rounds(&[(all[0].clone(), 1.5)], T).unwrap();
        let built = run(&mut m, &s, &request, T).unwrap().unwrap();
        let point = |v: &Shape| {
            m.node(v)
                .and_then(|n| n.data().as_vertex().map(|d| d.point))
                .unwrap()
        };
        for e in &all[1..] {
            let images = built.history.trace(e);
            assert_eq!(images.len(), 1);
            let faces =
                ogeom_topo::ancestors_of(&m, &built.shape, &images[0], ShapeType::Face).unwrap();
            assert_eq!(faces.len(), 2);
            let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
            let (x, y) = ogeom_algo::edge_vertices(&m, &images[0]).unwrap().unwrap();
            let d = point(&b) - point(&a);
            for p in [point(&x), point(&y)] {
                assert!((p - point(&a)).cross(d).magnitude() < 1e-9);
            }
        }
        let (a, b) = ogeom_algo::edge_vertices(&m, &all[0]).unwrap().unwrap();
        for v in explore_unique(&m, &s, ShapeType::Vertex).unwrap() {
            if v.is_same(&a) || v.is_same(&b) {
                continue;
            }
            let images = built.history.trace(&v);
            assert_eq!(images.len(), 1);
            assert!(point(&images[0]).distance(point(&v)) < 1e-12);
        }
    }

    /// A 10 by 8 by 6 box, the two top edges at the corner `(0, 0, 6)`
    /// rounded at 1.5 and the upright edge there left sharp: the two
    /// strips run on to the mitre plane `x = y` and meet along its ellipse,
    /// the crease carried on from the sharp edge. The box less each strip's
    /// spandrel along its whole edge, less what the two would both take at
    /// the corner, `r^3 (5/3 - pi/2)` (the integral of the two spandrels'
    /// common square over the corner's height).
    #[test]
    fn two_rounds_beside_a_sharp_edge_meet_in_their_mitre() {
        let r = 1.5;
        let mut m = Model::new();
        let s = ogeom_algo::make_box(&mut m, Frame::WORLD, (10.0, 8.0, 6.0), T)
            .unwrap()
            .shape;
        let point = |m: &Model, v: &Shape| {
            m.node(v)
                .and_then(|n| n.data().as_vertex().map(|d| d.point))
                .unwrap()
        };
        let top_at_origin = |m: &Model, e: &Shape| {
            let (a, b) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
            let (p, q) = (point(m, &a), point(m, &b));
            (p.z - 6.0).abs() < 1e-9
                && (q.z - 6.0).abs() < 1e-9
                && (p.x.abs() + p.y.abs() < 1e-9 || q.x.abs() + q.y.abs() < 1e-9)
        };
        let pair: Vec<(Shape, f64)> = explore_unique(&m, &s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| top_at_origin(&m, e))
            .map(|e| (e, r))
            .collect();
        assert_eq!(pair.len(), 2);
        let mut built = Vec::new();
        for order in [pair.clone(), pair.iter().rev().cloned().collect()] {
            let request = Request::from_rounds(&order, T).unwrap();
            let b = run(&mut m, &s, &request, T).unwrap().unwrap().shape;
            let found = ogeom_algo::check(&m, &b, T).unwrap();
            assert!(found.is_valid(), "{:?}", found.problems);
            built.push(b);
        }
        assert_eq!(bits(&m, &built[0]), bits(&m, &built[1]));
        let pi = core::f64::consts::PI;
        let spandrel = (1.0 - pi / 4.0) * r * r;
        let want = 480.0 - spandrel * (10.0 + 8.0) + r * r * r * (5.0 / 3.0 - pi / 2.0);
        let v = volume(&m, &built[0]);
        assert!((v - want).abs() < 1e-6, "{v} against {want}");
    }

    /// The pyramid over the rectangle `|x| <= a`, `|y| <= b` at `z = 0`
    /// with its apex at `(0, 0, h)`.
    fn pyramid(m: &mut Model, a: f64, b: f64, h: f64) -> Shape {
        use ogeom_math::Point;
        let points = [
            Point::new(-a, -b, 0.0),
            Point::new(a, -b, 0.0),
            Point::new(a, b, 0.0),
            Point::new(-a, b, 0.0),
            Point::new(0.0, 0.0, h),
        ];
        let rings = vec![
            vec![0, 3, 2, 1],
            vec![0, 1, 4],
            vec![1, 2, 4],
            vec![2, 3, 4],
            vec![3, 0, 4],
        ];
        ogeom_algo::make_polyhedron(m, &points, &rings, T)
            .unwrap()
            .shape
    }

    /// The edges of `s` through the point `p`.
    fn edges_through(m: &Model, s: &Shape, p: ogeom_math::Point) -> Vec<Shape> {
        explore_unique(m, s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                let (a, b) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
                [a, b].iter().any(|v| {
                    let q = m.node(v).unwrap().data().as_vertex().unwrap().point;
                    q.distance(p) < 1e-9
                })
            })
            .collect()
    }

    /// The square pyramid of side 6 and height 4 has an insphere of radius
    /// 3/2, its four sides at the apex tangent to one ball as at every
    /// corner. Every edge rounded at 1/2, each vertex closes on a sphere,
    /// the apex's a spherical quadrilateral, and the solid is the inner
    /// parallel body at 1/2 (the pyramid scaled by 2/3 about the insphere's
    /// centre) swept by the ball, whose volume is Steiner's
    /// `V' + S' r + M' r^2 + 4/3 pi r^3`, `M'` the sum over the edges of
    /// their lengths times their outer angles over two: the base edges'
    /// `pi - atan(4/3)`, the slant edges' (of length `sqrt 34`)
    /// `acos(9/25)`, between their sides' normals `(0, -4, 3)/5` and
    /// `(4, 0, 3)/5`.
    #[test]
    fn a_square_pyramid_s_apex_rounds_on_the_ball_its_four_sides_touch() {
        use core::f64::consts::PI;
        let mut m = Model::new();
        let s = pyramid(&mut m, 3.0, 3.0, 4.0);
        let r = 0.5;
        let edges = explore_unique(&m, &s, ShapeType::Edge).unwrap();
        assert_eq!(edges.len(), 8);
        let request =
            Request::from_rounds(&edges.iter().map(|e| (e.clone(), r)).collect::<Vec<_>>(), T)
                .unwrap();
        let built = run(&mut m, &s, &request, T).unwrap().unwrap().shape;
        let found = ogeom_algo::check(&m, &built, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        assert_eq!(count_of(&m, &built, ogeom_geom::SurfaceKind::Sphere), 5);
        let k = (1.5 - r) / 1.5;
        let (v, area) = (48.0, 96.0);
        let mean = 4.0 * 6.0 * (PI - 4.0f64.atan2(3.0)) / 2.0
            + 4.0 * 34.0f64.sqrt() * (9.0f64 / 25.0).acos() / 2.0;
        let want = k.powi(3) * v + k * k * area * r + k * mean * r * r + 4.0 / 3.0 * PI * r.powi(3);
        let got = volume(&m, &built);
        assert!((got - want).abs() < 1e-6, "{got} against {want}");
    }

    /// The right pyramid over a regular hexagon of circumradius 3, height
    /// 4: its six sides at the apex touch one ball.
    /// Every edge rounded at 1/2, every vertex closes on a sphere, the
    /// apex's a spherical hexagon, and the volume is Steiner's.
    #[test]
    fn a_hexagonal_pyramid_s_apex_of_six_edges_rounds_on_one_sphere() {
        use ogeom_math::Point;
        let mut points: Vec<Point> = (0..6)
            .map(|k| {
                let t = f64::from(k) * core::f64::consts::PI / 3.0;
                Point::new(3.0 * t.cos(), 3.0 * t.sin(), 0.0)
            })
            .collect();
        points.push(Point::new(0.0, 0.0, 4.0));
        let mut rings = vec![vec![5, 4, 3, 2, 1, 0]];
        rings.extend((0..6).map(|k| vec![k, (k + 1) % 6, 6]));
        let mut m = Model::new();
        let s = ogeom_algo::make_polyhedron(&mut m, &points, &rings, T)
            .unwrap()
            .shape;
        let r = 0.5;
        let edges = explore_unique(&m, &s, ShapeType::Edge).unwrap();
        assert_eq!(edges.len(), 12);
        let request =
            Request::from_rounds(&edges.iter().map(|e| (e.clone(), r)).collect::<Vec<_>>(), T)
                .unwrap();
        let built = run(&mut m, &s, &request, T).unwrap().unwrap().shape;
        let found = ogeom_algo::check(&m, &built, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        assert_eq!(count_of(&m, &built, ogeom_geom::SurfaceKind::Sphere), 7);
        let want = steiner(&points, &rings, r);
        let got = volume(&m, &built);
        assert!((got - want).abs() < 1e-6, "{got} against {want}");
    }

    /// The rectangular pyramid's apex `edges` rounded at `r`, through the
    /// engine: valid, its two spheres and the cylinder between them the
    /// ball rolling along the ridge (every point of each the contact of a
    /// ball of radius `r` inside the four slant planes, touching two of them
    /// or more), no patch, and the pyramid's two mirrors kept (its centroid
    /// on the axis).
    fn rolls_along_its_ridge(a: f64, b: f64, h: f64, r: f64) {
        use ogeom_math::{Point, Vector};
        let mut m = Model::new();
        let s = pyramid(&mut m, a, b, h);
        let apex = Point::new(0.0, 0.0, h);
        let slant = edges_through(&m, &s, apex);
        assert_eq!(slant.len(), 4);
        let request =
            Request::from_rounds(&slant.iter().map(|e| (e.clone(), r)).collect::<Vec<_>>(), T)
                .unwrap();
        let built = run(&mut m, &s, &request, T).unwrap().unwrap().shape;
        let found = ogeom_algo::check(&m, &built, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        assert_eq!(count_of(&m, &built, ogeom_geom::SurfaceKind::BSpline), 0);
        assert_eq!(count_of(&m, &built, ogeom_geom::SurfaceKind::Sphere), 2);
        assert_eq!(count_of(&m, &built, ogeom_geom::SurfaceKind::Cylinder), 5);
        // The slant planes, each an outward unit normal through the apex.
        let base = [
            Point::new(-a, -b, 0.0),
            Point::new(a, -b, 0.0),
            Point::new(a, b, 0.0),
            Point::new(-a, b, 0.0),
        ];
        let planes: Vec<Vector> = (0..4)
            .map(|k| {
                let n = (base[(k + 1) % 4] - base[k]).cross(apex - base[k]);
                n / n.magnitude()
            })
            .collect();
        let fine = ogeom_mesh::Deflection {
            chord: 1e-3,
            ..ogeom_mesh::Deflection::default()
        };
        let mut checked = 0;
        for face in explore_unique(&m, &built, ShapeType::Face).unwrap() {
            let d = m.node(&face).unwrap().data().as_face().unwrap();
            let surface = m.geometry().surface(d.surface).unwrap();
            let apex_side = matches!(surface, ogeom_geom::SurfaceGeometry::Sphere(_))
                || matches!(surface, ogeom_geom::SurfaceGeometry::Cylinder(c)
                    if c.cylinder().axis().direction.vector().z.abs() < 1e-9);
            if !apex_side {
                continue;
            }
            let mesh = ogeom_mesh::triangulate_face(&m, &face, fine, T).unwrap();
            for (p, n) in mesh.positions.iter().zip(&mesh.normals) {
                // To the mesh's own normals' accuracy.
                let centre = *p - *n * r;
                let gaps: Vec<f64> = planes.iter().map(|q| -(centre - apex).dot(*q)).collect();
                assert!(
                    gaps.iter().all(|g| *g >= r - 1e-5),
                    "a ball past a side: {gaps:?}"
                );
                let touching = gaps.iter().filter(|g| (**g - r).abs() <= 1e-5).count();
                assert!(touching >= 2, "a ball touching {touching} sides at {p:?}");
                checked += 1;
            }
        }
        assert!(checked > 0);
        let centre = ogeom_algo::volume_properties(
            &m,
            &built,
            ogeom_mesh::Deflection {
                chord: 1e-5,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap()
        .centre;
        assert!(centre.x.abs() < 1e-6 && centre.y.abs() < 1e-6, "{centre:?}");
    }

    /// The pyramid over an 8 by 6 rectangle, height 4: its long sides lean
    /// at `atan(4/3)`, its short ones at 45 degrees, so no one ball touches
    /// the four at the apex; one touches the long sides and a short one at
    /// each end of the ridge along `x`. Its four slant edges rounded at 1/2,
    /// the apex rounds as the ball rolls along that ridge: a sphere at each
    /// end and the cylinder between ([`rolls_along_its_ridge`]).
    #[test]
    fn a_rectangular_pyramid_s_apex_rounds_as_the_ball_rolls_along_its_ridge() {
        rolls_along_its_ridge(4.0, 3.0, 4.0, 0.5);
    }

    /// The pyramid over a 20 by 10 rectangle, height 15, its slant edges
    /// rounded at 1.5: the hole the strips leave at the apex has corners
    /// turning back on themselves, which no patch over a square takes
    /// soundly; the ball rolling along the ridge rounds it exactly
    /// ([`rolls_along_its_ridge`]), as `edge_chains.rs` asks of the public
    /// call.
    #[test]
    fn an_elongated_pyramid_s_apex_rounds_as_the_ball_rolls_along_its_ridge() {
        rolls_along_its_ridge(10.0, 5.0, 15.0, 1.5);
    }

    /// A roof whose ridge rises along `x` through a saddle where a valley
    /// falls along `y`: four triangles round the vertex at the origin, its
    /// edges to `(+-2, 0, 1)` convex and to `(0, +-2, -1)` concave, walled
    /// down to `z = -3`. Every edge rounded at 0.25: the saddle, a vertex of
    /// four edges of mixed convexity, closes on one patch (ADR 0001, §13,
    /// question 7), held to the bands' fit target and proved within the
    /// solid's own material there (the ridge's sides meet across the
    /// convex edges, the runs between unite across the concave ones) by the
    /// engine. No closed form: valid, the same solid in either order, and
    /// its two mirrors kept (the centroid on the axis).
    #[test]
    fn a_saddle_vertex_of_four_edges_of_mixed_convexity_closes_on_one_patch() {
        use ogeom_math::Point;
        let points = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(2.0, 0.0, 1.0),
            Point::new(0.0, 2.0, -1.0),
            Point::new(-2.0, 0.0, 1.0),
            Point::new(0.0, -2.0, -1.0),
            Point::new(2.0, 0.0, -3.0),
            Point::new(0.0, 2.0, -3.0),
            Point::new(-2.0, 0.0, -3.0),
            Point::new(0.0, -2.0, -3.0),
        ];
        let rings = vec![
            vec![0, 1, 2],
            vec![0, 2, 3],
            vec![0, 3, 4],
            vec![0, 4, 1],
            vec![1, 5, 6, 2],
            vec![2, 6, 7, 3],
            vec![3, 7, 8, 4],
            vec![4, 8, 5, 1],
            vec![5, 8, 7, 6],
        ];
        let mut m = Model::new();
        let s = ogeom_algo::make_polyhedron(&mut m, &points, &rings, T)
            .unwrap()
            .shape;
        let edges = explore_unique(&m, &s, ShapeType::Edge).unwrap();
        let mut built = Vec::new();
        for order in [edges.clone(), edges.iter().rev().cloned().collect()] {
            let request = Request::from_rounds(
                &order.iter().map(|e| (e.clone(), 0.25)).collect::<Vec<_>>(),
                T,
            )
            .unwrap();
            let b = run(&mut m, &s, &request, T).unwrap().unwrap().shape;
            let found = ogeom_algo::check(&m, &b, T).unwrap();
            assert!(found.is_valid(), "{:?}", found.problems);
            built.push(b);
        }
        assert_eq!(bits(&m, &built[0]), bits(&m, &built[1]));
        assert!(count_of(&m, &built[0], ogeom_geom::SurfaceKind::BSpline) >= 1);
        let centre = ogeom_algo::volume_properties(
            &m,
            &built[0],
            ogeom_mesh::Deflection {
                chord: 1e-5,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap()
        .centre;
        assert!(centre.x.abs() < 1e-6 && centre.y.abs() < 1e-6, "{centre:?}");
    }

    /// The pyramid over the pentagon `(-3, -2)`, `(3, -2)`, `(3.5, 1.5)`,
    /// `(0, 3.5)`, `(-3.5, 1.5)`, height 4: no one ball touches its five
    /// sides at the apex. Its five slant edges rounded at 1/2, the apex
    /// closes on one five-sided patch, a face on each side's sector, held
    /// to the bands' fit target on its strips and between its faces and
    /// proved inside the five planes by the engine. The solid is valid and
    /// keeps the pyramid's mirror `x = 0`: its centroid stands on it.
    #[test]
    fn a_pentagonal_pyramid_s_apex_closes_on_one_five_sided_patch() {
        use ogeom_math::Point;
        let points = [
            Point::new(-3.0, -2.0, 0.0),
            Point::new(3.0, -2.0, 0.0),
            Point::new(3.5, 1.5, 0.0),
            Point::new(0.0, 3.5, 0.0),
            Point::new(-3.5, 1.5, 0.0),
            Point::new(0.0, 0.0, 4.0),
        ];
        let mut rings = vec![vec![4, 3, 2, 1, 0]];
        rings.extend((0..5).map(|k| vec![k, (k + 1) % 5, 5]));
        let mut m = Model::new();
        let s = ogeom_algo::make_polyhedron(&mut m, &points, &rings, T)
            .unwrap()
            .shape;
        let slant = edges_through(&m, &s, points[5]);
        assert_eq!(slant.len(), 5);
        let request = Request::from_rounds(
            &slant.iter().map(|e| (e.clone(), 0.5)).collect::<Vec<_>>(),
            T,
        )
        .unwrap();
        let built = run(&mut m, &s, &request, T).unwrap().unwrap().shape;
        let found = ogeom_algo::check(&m, &built, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        assert_eq!(count_of(&m, &built, ogeom_geom::SurfaceKind::BSpline), 5);
        let fine = ogeom_mesh::Deflection {
            chord: 1e-5,
            ..ogeom_mesh::Deflection::default()
        };
        let centre = ogeom_algo::volume_properties(&m, &built, fine, T)
            .unwrap()
            .centre;
        assert!(centre.x.abs() < 1e-6, "{centre:?}");
    }

    /// The edge of `s` whose ends' midpoint is `mid`.
    fn edge_at(m: &Model, s: &Shape, mid: ogeom_math::Point) -> Shape {
        explore_unique(m, s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let (a, b) = ogeom_algo::edge_vertices(m, e).unwrap().unwrap();
                let p = |v: &Shape| m.node(v).unwrap().data().as_vertex().unwrap().point;
                p(&a).lerp(p(&b), 0.5).distance(mid) < 1e-9
            })
            .unwrap()
    }

    /// The overlap of two rounds' corner regions at a box corner, the round
    /// along `x` of radius `r1` and the one along `y` of `r2`, with
    /// `r1 <= r2`: over the first's section angle `phi`, the first's region
    /// stands `r1 (1 - sin phi)` deep in `y` at height `z = r1 (1 - cos phi)`
    /// and the second's `r2 - sqrt(z (2 r2 - z))` deep in `x` there, so the
    /// overlap is the integral of their product over `dz`. Smooth in `phi`,
    /// so 64-point Gauss–Legendre takes it to rounding.
    fn corner_overlap(r1: f64, r2: f64) -> f64 {
        // Nodes and weights by Newton's method on the Legendre polynomial.
        let n = 64;
        let mut sum = 0.0;
        let half = core::f64::consts::FRAC_PI_4;
        for i in 0..n {
            let mut x =
                (core::f64::consts::PI * (f64::from(i) + 0.75) / (f64::from(n) + 0.5)).cos();
            let mut dp = 0.0;
            for _ in 0..100 {
                let (mut p0, mut p1) = (1.0, x);
                for k in 2..=n {
                    let p2 = ((2.0 * f64::from(k) - 1.0) * x * p1 - (f64::from(k) - 1.0) * p0)
                        / f64::from(k);
                    p0 = p1;
                    p1 = p2;
                }
                dp = f64::from(n) * (x * p1 - p0) / (x * x - 1.0);
                let step = p1 / dp;
                x -= step;
                if step.abs() < 1e-16 {
                    break;
                }
            }
            let w = 2.0 / ((1.0 - x * x) * dp * dp);
            let phi = half * (x + 1.0);
            let z = r1 * (1.0 - phi.cos());
            let y1 = r1 * (1.0 - phi.sin());
            let x2 = r2 - (z * (2.0 * r2 - z)).sqrt();
            sum += w * half * y1 * x2 * r1 * phi.sin();
        }
        sum
    }

    /// The overlap at equal radii is the mitre's, `r^3 (5/3 - pi/2)`: the
    /// quadrature checked against the closed form it generalises.
    #[test]
    fn the_corner_overlap_at_equal_radii_is_the_mitre_s() {
        let r: f64 = 1.5;
        let want = r.powi(3) * (5.0 / 3.0 - core::f64::consts::FRAC_PI_2);
        assert!(
            (corner_overlap(r, r) - want).abs() < 1e-13,
            "{}",
            corner_overlap(r, r)
        );
    }

    /// A 10 by 8 by 6 box, the corner at the origin's edge along `x` rounded
    /// at 1 and its edge along `y` at 2, the edge along `z` left sharp: the
    /// larger round runs on to the face `y = 0` and ends flush there on its
    /// section, and the smaller ends on the larger along the curve the two
    /// cylinders cross in. The solid is the box less each round's corner
    /// region along its edge, plus their overlap at the corner (the union
    /// of the two regions removed), in either order.
    #[test]
    fn two_rounds_of_different_radii_beside_a_sharp_edge_meet_in_their_crease() {
        let (r1, r2) = (1.0, 2.0);
        let corner = 1.0 - core::f64::consts::FRAC_PI_4;
        let want = 480.0 - corner * (r1 * r1 * 10.0 + r2 * r2 * 8.0) + corner_overlap(r1, r2);
        let mut volumes = Vec::new();
        for flip in [false, true] {
            let mut m = Model::new();
            let s = ogeom_algo::make_box(&mut m, Frame::WORLD, (10.0, 8.0, 6.0), T)
                .unwrap()
                .shape;
            let along_x = edge_at(&m, &s, ogeom_math::Point::new(5.0, 0.0, 0.0));
            let along_y = edge_at(&m, &s, ogeom_math::Point::new(0.0, 4.0, 0.0));
            let mut asked = vec![(along_x, r1), (along_y, r2)];
            if flip {
                asked.reverse();
            }
            let request = Request::from_rounds(&asked, T).unwrap();
            let built = run(&mut m, &s, &request, T).unwrap().unwrap().shape;
            let found = ogeom_algo::check(&m, &built, T).unwrap();
            assert!(found.is_valid(), "{:?}", found.problems);
            volumes.push(volume(&m, &built));
        }
        for v in &volumes {
            assert!((v - want).abs() < 1e-6, "{v} against {want}");
        }
    }

    /// A 10 by 8 by 6 block with the corner `x < 4`, `y < 3`, `z > 3` cut
    /// away: at the notch's inner floor corner the floor's edges along the
    /// walls `y = 3` (along `x`, 4 long) and `x = 4` (along `y`, 3 long)
    /// are concave, rounded at 1/2 and 1, and the walls' own edge is left
    /// sharp. The larger fill runs on to the wall `y = 3` and ends flush in
    /// it, the smaller ends on the larger along their cylinders' crease, and
    /// each fill's other end is flush in the block's side. The solid is the
    /// notched block plus each fill's corner region along its edge, less
    /// their overlap at the corner (the union of the two fills), in either
    /// order.
    #[test]
    fn two_concave_rounds_of_different_radii_beside_a_sharp_edge_meet_in_their_crease() {
        use ogeom_math::Point;
        let p = Point::new;
        let points = [
            p(0.0, 0.0, 0.0),
            p(10.0, 0.0, 0.0),
            p(10.0, 8.0, 0.0),
            p(0.0, 8.0, 0.0),
            p(10.0, 0.0, 6.0),
            p(10.0, 8.0, 6.0),
            p(0.0, 8.0, 6.0),
            p(0.0, 3.0, 6.0),
            p(4.0, 3.0, 6.0),
            p(4.0, 0.0, 6.0),
            p(0.0, 0.0, 3.0),
            p(4.0, 0.0, 3.0),
            p(4.0, 3.0, 3.0),
            p(0.0, 3.0, 3.0),
        ];
        let rings = vec![
            vec![0, 3, 2, 1],
            vec![9, 4, 5, 6, 7, 8],
            vec![10, 11, 12, 13],
            vec![0, 1, 4, 9, 11, 10],
            vec![3, 6, 5, 2],
            vec![0, 10, 13, 7, 6, 3],
            vec![1, 2, 5, 4],
            vec![11, 9, 8, 12],
            vec![13, 12, 8, 7],
        ];
        let (ra, rb) = (0.5, 1.0);
        let corner = 1.0 - core::f64::consts::FRAC_PI_4;
        let want = 444.0 + corner * (ra * ra * 4.0 + rb * rb * 3.0) - corner_overlap(ra, rb);
        for flip in [false, true] {
            let mut m = Model::new();
            let (s, _) = ogeom_algo::make_faceted_solid(&mut m, &points, &rings, &[], T).unwrap();
            let s = s.shape;
            let along_x = edge_at(&m, &s, p(2.0, 3.0, 3.0));
            let along_y = edge_at(&m, &s, p(4.0, 1.5, 3.0));
            let mut asked = vec![(along_x, ra), (along_y, rb)];
            if flip {
                asked.reverse();
            }
            let request = Request::from_rounds(&asked, T).unwrap();
            let built = run(&mut m, &s, &request, T).unwrap().unwrap().shape;
            let found = ogeom_algo::check(&m, &built, T).unwrap();
            assert!(found.is_valid(), "{:?}", found.problems);
            let v = volume(&m, &built);
            assert!((v - want).abs() < 1e-6, "{v} against {want}");
        }
    }

    /// A 4 by 3 rectangle at `x = 0` swept along `(10, 3, 0)`: its caps
    /// meet the sweep edges leaning, and each sweep edge's two side faces
    /// stand square to each other (both hold the sweep and the rectangle's
    /// sides are square). One sweep edge rounded at 1: its strip ends in
    /// each cap on the ellipse the cap's plane cuts, and the prism less the
    /// round's square corner along the edge's length between the caps,
    /// `sqrt 109`, is the volume (Cavalieri).
    #[test]
    fn a_round_ends_on_the_ellipse_a_leaning_cap_cuts() {
        use ogeom_math::{Direction, Plane, Point, Vector};
        let mut m = Model::new();
        let points = [
            Point::new(0.0, 0.0, 0.0),
            Point::new(0.0, 4.0, 0.0),
            Point::new(0.0, 4.0, 3.0),
            Point::new(0.0, 0.0, 3.0),
        ];
        let wire = ogeom_algo::make_polygon(&mut m, &points, true, T)
            .unwrap()
            .shape;
        let plane = ogeom_geom::PlaneSurface::over(
            Plane::through(Point::ORIGIN, Direction::X),
            (-100.0, 100.0),
            (-100.0, 100.0),
        )
        .unwrap();
        let face = ogeom_algo::make_face(&mut m, plane.into(), &[wire], T)
            .unwrap()
            .shape;
        let swept = ogeom_algo::make_prism(&mut m, &face, Vector::new(10.0, 3.0, 0.0), T)
            .unwrap()
            .shape;
        let s = ogeom_algo::baked_shape(&mut m, &swept, T).unwrap().shape;
        let point = |m: &Model, v: &Shape| {
            m.node(v)
                .and_then(|n| n.data().as_vertex().map(|d| d.point))
                .unwrap()
        };
        let edge = explore_unique(&m, &s, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let (a, b) = ogeom_algo::edge_vertices(&m, e).unwrap().unwrap();
                let (p, q) = (point(&m, &a), point(&m, &b));
                let start = if p.x < q.x { p } else { q };
                start.distance(Point::ORIGIN) < 1e-9 && (p.x - q.x).abs() > 1.0
            })
            .unwrap();
        let r = 1.0;
        let request = Request::from_rounds(&[(edge, r)], T).unwrap();
        let built = run(&mut m, &s, &request, T).unwrap().unwrap().shape;
        let found = ogeom_algo::check(&m, &built, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        let pi = core::f64::consts::PI;
        let want = 120.0 - (1.0 - pi / 4.0) * r * r * 109.0f64.sqrt();
        let v = volume(&m, &built);
        assert!((v - want).abs() < 1e-6, "{v} against {want}");
    }

    /// A box with every edge rounded: the exact rounded box,
    /// `(a-2r)(b-2r)(c-2r) + 2r(sum of pairs) + pi r^2 (sum) + 4/3 pi r^3`
    /// over the shrunk sides.
    #[test]
    fn every_round_meets_in_eight_spheres() {
        let r = 1.5;
        let (m, built) = rounded(<[Shape]>::to_vec, r);
        let built = built.unwrap();
        let (a, b, c) = (10.0 - 2.0 * r, 8.0 - 2.0 * r, 6.0 - 2.0 * r);
        let pi = core::f64::consts::PI;
        let want = a * b * c
            + 2.0 * r * (a * b + b * c + c * a)
            + pi * r * r * (a + b + c)
            + 4.0 / 3.0 * pi * r * r * r;
        let v = volume(&m, &built);
        assert!((v - want).abs() < 1e-6, "{v} against {want}");
        let kinds: Vec<String> = explore_unique(&m, &built, ShapeType::Face)
            .unwrap()
            .iter()
            .map(|f| {
                let Some(ogeom_topo::NodeData::Face(d)) = m.node(f).map(|n| n.data()) else {
                    return String::new();
                };
                use ogeom_geom::Surface as _;
                format!("{:?}", m.geometry().surface(d.surface).unwrap().kind())
            })
            .collect();
        let count = |k: &str| kinds.iter().filter(|s| *s == k).count();
        assert_eq!(
            [count("Plane"), count("Cylinder"), count("Sphere")],
            [6, 12, 8]
        );
    }

    /// The volume of the convex polytope of `points` and `rings` (wound
    /// outward) with every edge and corner rounded at `r`, where its faces
    /// moved in by `r` still meet as it does, each vertex's at one point:
    /// that inner parallel body swept by the ball, by Steiner,
    /// `V' + S' r + M' r^2 + 4/3 pi r^3`, `M'` the sum over its edges of
    /// their lengths times the angles between their faces' outward normals,
    /// over two.
    fn steiner(points: &[ogeom_math::Point], rings: &[Vec<usize>], r: f64) -> f64 {
        use ogeom_math::{Point, Vector};
        let normal = |ring: &Vec<usize>, at: &[Point]| {
            let mut n = Vector::ZERO;
            for k in 1..ring.len() - 1 {
                let (a, b, c) = (at[ring[0]], at[ring[k]], at[ring[k + 1]]);
                n += (b - a).cross(c - a);
            }
            n / n.magnitude()
        };
        let normals: Vec<Vector> = rings.iter().map(|ring| normal(ring, points)).collect();
        // Each vertex moved to where three of its faces' planes, moved in
        // by `r`, meet.
        let inner: Vec<Point> = (0..points.len())
            .map(|i| {
                let faces: Vec<usize> = (0..rings.len())
                    .filter(|f| rings[*f].contains(&i))
                    .collect();
                let n = [normals[faces[0]], normals[faces[1]], normals[faces[2]]];
                let d = n.map(|n| n.dot(points[i].to_vector()) - r);
                let det = n[0].dot(n[1].cross(n[2]));
                Point::from_vector(
                    (n[1].cross(n[2]) * d[0] + n[2].cross(n[0]) * d[1] + n[0].cross(n[1]) * d[2])
                        / det,
                )
            })
            .collect();
        let o = inner[0];
        let (mut volume, mut area, mut mean) = (0.0, 0.0, 0.0);
        for (i, ring) in rings.iter().enumerate() {
            let mut n = Vector::ZERO;
            for k in 1..ring.len() - 1 {
                let (a, b, c) = (inner[ring[0]], inner[ring[k]], inner[ring[k + 1]]);
                n += (b - a).cross(c - a) * 0.5;
                volume += (a - o).dot((b - o).cross(c - o)) / 6.0;
            }
            area += n.magnitude();
            for k in 0..ring.len() {
                let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
                // The other face holds the edge walked the other way.
                let j = rings
                    .iter()
                    .position(|q| (0..q.len()).any(|m| q[m] == b && q[(m + 1) % q.len()] == a))
                    .unwrap();
                if i < j {
                    let angle = normals[i].dot(normals[j]).clamp(-1.0, 1.0).acos();
                    mean += inner[a].distance(inner[b]) * angle / 2.0;
                }
            }
        }
        volume + area * r + mean * r * r + 4.0 / 3.0 * core::f64::consts::PI * r.powi(3)
    }

    /// A frustum of a hexagonal pyramid, a regular hexagon of circumradius
    /// 3 at `z = 0` narrowing to one of 1.5 at `z = 2`: every vertex joins
    /// three edges, the base's at an obtuse corner whose ball's foot on the
    /// base stands nearly a quarter turn from the middle of its triangle.
    /// Every edge rounded at 1/2, every vertex closes on a sphere and the
    /// volume is Steiner's; a sphere charted with a pole beside that foot
    /// leaves the solid unsound.
    #[test]
    fn a_hexagonal_frustum_s_obtuse_corners_round_on_their_spheres() {
        use ogeom_math::Point;
        let ring = |radius: f64, z: f64| {
            (0..6).map(move |k| {
                let t = f64::from(k) * core::f64::consts::PI / 3.0;
                Point::new(radius * t.cos(), radius * t.sin(), z)
            })
        };
        let points: Vec<Point> = ring(3.0, 0.0).chain(ring(1.5, 2.0)).collect();
        let mut rings = vec![vec![5, 4, 3, 2, 1, 0], vec![6, 7, 8, 9, 10, 11]];
        rings.extend((0..6).map(|k| vec![k, (k + 1) % 6, 6 + (k + 1) % 6, 6 + k]));
        let mut m = Model::new();
        let s = ogeom_algo::make_polyhedron(&mut m, &points, &rings, T)
            .unwrap()
            .shape;
        let r = 0.5;
        let edges = explore_unique(&m, &s, ShapeType::Edge).unwrap();
        assert_eq!(edges.len(), 18);
        let request =
            Request::from_rounds(&edges.iter().map(|e| (e.clone(), r)).collect::<Vec<_>>(), T)
                .unwrap();
        let built = run(&mut m, &s, &request, T).unwrap().unwrap().shape;
        let found = ogeom_algo::check(&m, &built, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        assert_eq!(count_of(&m, &built, ogeom_geom::SurfaceKind::Sphere), 12);
        let want = steiner(&points, &rings, r);
        let got = volume(&m, &built);
        assert!((got - want).abs() < 1e-6, "{got} against {want}");
    }
}
