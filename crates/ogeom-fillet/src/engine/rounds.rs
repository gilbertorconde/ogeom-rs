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

use ogeom_algo::{Built, History};
use ogeom_core::{FastMap, OgeomResult, Tolerances};
use ogeom_geom::SurfaceGeometry;
use ogeom_math::Point;
use ogeom_topo::Model;

use super::Id;
use super::assemble::{SheetEdge, SheetFace, build_faces};
use super::classify::{Planar, Seat, faces_round};
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

    // Each strip's station at each end: both its faces' setback corners
    // there stand at one station along the edge, or the strip's end is not
    // its section and this class does not hold it.
    let reach = tol.confusion() * 100.0;
    let mut station: FastMap<(Id, Id), f64> = FastMap::default();
    for (edge, (section, _)) in &sections {
        let seat = edge_seat[edge];
        for v in [seat.ends.0, seat.ends.1] {
            let s = [0, 1]
                .map(|k| (corners[&(seat.faces[k], v)] - section.crease(0.0)).dot(seat.along));
            if (s[0] - s[1]).abs() > reach {
                return Ok(None);
            }
            station.insert((*edge, v), 0.5 * (s[0] + s[1]));
        }
    }

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
                let seat = blended[0];
                let Some(third) = round.iter().copied().find(|f| !seat.faces.contains(f)) else {
                    return Ok(None);
                };
                let square =
                    planar.faces[third].normal.cross(seat.along).magnitude() <= tol.angular();
                let s = station[&(seat.edge, v)];
                let at_end = s.abs() <= reach || (s - section_length(planar, seat)).abs() <= reach;
                if !(square && at_end) {
                    return Ok(None);
                }
                flush.insert(v, (third, seat.edge));
                RoundVertex::Flush
            }
            (3, 3) => {
                let radius = sections[&blended[0].edge].1;
                let convex = blended[0].convex;
                if blended
                    .iter()
                    .any(|s| s.convex != convex || sections[&s.edge].1 != radius)
                {
                    return Ok(None);
                }
                let planes = [round[0], round[1], round[2]]
                    .map(|f| (planar.faces[f].normal, planar.points[&v]));
                let Some((_, feet)) =
                    crate::round::sphere::corner_ball(planes, radius, convex, tol)
                else {
                    return Ok(None);
                };
                for (k, f) in round.iter().enumerate() {
                    if feet[k].distance(corners[&(*f, v)]) > reach {
                        return Ok(None);
                    }
                }
                RoundVertex::Sphere
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
    // The faces' setback polygons: none crossing itself or another, none
    // turned over.
    let rings: Vec<super::assemble::Ring> = planar
        .faces
        .iter()
        .enumerate()
        .map(|(fi, face)| {
            let ring = face.corners.iter().map(|v| at[&(fi, *v)]).collect();
            (ring, vec![super::Origin::Face(fi)], Some(face.normal))
        })
        .collect();
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
    let mut strip_list: Vec<Id> = sections.keys().copied().collect();
    strip_list.sort_unstable();
    for edge in &strip_list {
        let (section, _) = &sections[edge];
        let seat = edge_seat[edge];
        for v in [seat.ends.0, seat.ends.1] {
            let (curve, range) = section.arc_at(station[&(*edge, v)], tol)?;
            let [from_face, _] = section.order();
            let (p0, p1) = (at[&(seat.faces[0], v)], at[&(seat.faces[1], v)]);
            // The arc runs from the first face of its order's contact.
            let ends = if from_face == 0 { (p0, p1) } else { (p1, p0) };
            edges.push(SheetEdge {
                ends,
                curve: Some((curve, range)),
            });
            arcs.insert((*edge, v), (edges.len() - 1, from_face == 0));
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
                _ => at_corner.push((vec![at[&(fi, v)]], None)),
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
    // The strips: along the first face's contact backward, the arc at the
    // edge's start, along the second face's contact, the arc at its end.
    for edge in &strip_list {
        let (section, _) = &sections[edge];
        let seat = edge_seat[edge];
        let (a, b) = seat.ends;
        let (s0, s1) = (station[&(*edge, a)], station[&(*edge, b)]);
        let mut ring = Vec::with_capacity(4);
        ring.push(segment(
            &mut edges,
            at[&(seat.faces[0], b)],
            at[&(seat.faces[0], a)],
        ));
        let (arc_a, _) = arcs[&(*edge, a)];
        ring.push((arc_a, edges[arc_a].ends.0 == at[&(seat.faces[0], a)]));
        ring.push(segment(
            &mut edges,
            at[&(seat.faces[1], a)],
            at[&(seat.faces[1], b)],
        ));
        let (arc_b, _) = arcs[&(*edge, b)];
        ring.push((arc_b, edges[arc_b].ends.0 == at[&(seat.faces[1], b)]));
        let length = (s1 - s0).abs().max(tol.confusion());
        let surface = section.surface((s0.min(s1) - length, s0.max(s1) + length), tol)?;
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
        let mut pending: Vec<(usize, usize, usize)> = here
            .iter()
            .map(|s| {
                let (arc, _) = arcs[&(s.edge, v)];
                (arc, edges[arc].ends.0, edges[arc].ends.1)
            })
            .collect();
        // Each strip walks its arc at `v` one way; the sphere the other.
        let strip_walk = |s: &Seat| -> bool {
            let (arc, _) = arcs[&(s.edge, v)];
            if v == s.ends.0 {
                edges[arc].ends.0 == at[&(s.faces[0], v)]
            } else {
                edges[arc].ends.0 == at[&(s.faces[1], v)]
            }
        };
        let walks: FastMap<usize, bool> = here
            .iter()
            .map(|s| (arcs[&(s.edge, v)].0, !strip_walk(s)))
            .collect();
        let mut ring: Vec<(usize, bool)> = Vec::new();
        let (first, ..) = pending.remove(0);
        ring.push((first, walks[&first]));
        let mut tip = if walks[&first] {
            edges[first].ends.1
        } else {
            edges[first].ends.0
        };
        while !pending.is_empty() {
            let Some(k) = pending
                .iter()
                .position(|(arc, a, b)| (walks[arc] && *a == tip) || (!walks[arc] && *b == tip))
            else {
                return Ok(None);
            };
            let (arc, a, b) = pending.remove(k);
            ring.push((arc, walks[&arc]));
            tip = if walks[&arc] { b } else { a };
        }
        let seat = here[0];
        let radius = sections[&seat.edge].1;
        let round = faces_round(&planar.faces, v).unwrap_or_default();
        let planes =
            [round[0], round[1], round[2]].map(|f| (planar.faces[f].normal, planar.points[&v]));
        let Some((centre, feet)) =
            crate::round::sphere::corner_ball(planes, radius, seat.convex, tol)
        else {
            return Ok(None);
        };
        let surface = crate::round::sphere::triangle_sphere(centre, radius, feet, tol)?;
        faces.push((
            SheetFace {
                surface: Some(surface),
                ring,
            },
            seat.convex,
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
    let built: Vec<_> = built
        .into_iter()
        .zip(&faces)
        .map(|(f, (_, out))| if *out { f } else { f.reversed() })
        .collect();
    let shell = ogeom_algo::make_shell(model, &built)?.shape;
    let solid = ogeom_algo::make_solid(model, &[shell])?.shape;
    super::guards::solid_fits(model, &solid, tol)?;

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
    Ok(Some(Built::new(solid, history)))
}

/// The points a face's boundary passes at one of its corners, in turn, and
/// the arc between them where a strip ends flush in the face.
type Corner = (Vec<usize>, Option<(usize, bool)>);

/// The length of a seat's edge.
fn section_length(planar: &Planar, seat: &Seat) -> f64 {
    planar.points[&seat.ends.0].distance(planar.points[&seat.ends.1])
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
}
