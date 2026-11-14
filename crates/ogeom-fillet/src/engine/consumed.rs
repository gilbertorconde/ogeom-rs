//! Corners the setback rings cannot close alone, cut back by the planes
//! the blend puts there: faces a request consumes, whole or in part (a
//! thin wall bevelled on both its edges, a plate bevelled all round past
//! half its thickness, a face whose bevels cross on it at an angle), and a
//! convex corner where two bevels of unequal setbacks meet beside a sharp
//! edge, which mitres, as does a convex corner of four bevels or more
//! asked to mitre whose planes share no point.
//!
//! A face is consumed where its setback corners no longer make a ring
//! turning its way. Round such faces, and at a convex corner, every bevel
//! is convex, so what the blend leaves there is the part on the inner side
//! of every plane round it: each face's, each bevel's, each corner facet's
//! where three bevels or more meet. Each bevel there runs on past the
//! corner, and each ring touching it (the faces round it, the bevels, the
//! facets, a consumed face's own outline) is clipped by every one of those
//! planes but its own. A consumed face keeps what of its outline survives,
//! or goes; two bevels meet along the ridge their planes share; a face
//! across a consumed one is notched by the bevels; two bevels at a corner
//! of two meet along the line their planes share, with no facet between.
//!
//! A consumed face may hold holes, each convex, and reflex corners, every
//! edge of a hole or at a reflex corner bevelled. There the material is
//! not the inner side of every plane but the union of the sides a group of
//! bevels leaves: round a hole, its bevels', and at a reflex corner, its
//! two bevels'. Every other ring keeps what of it lies under any bevel of
//! the group: it loses the funnel above all of them. The group's own
//! bevels and walls keep their setback corners, which mitre at the hole's
//! corners as on any face, and do not run on; the walls are no cut for the
//! rest, which lie outside the hole.
//!
//! The clip is local by construction: a ring may lose only the setback
//! corners at those vertices, a consumed face's own corners, a bevel's
//! run-on, or points the clip itself made. A clip that would take any
//! other corner reaches past them, and the engine declines. A consumed face
//! with a hole that is not convex or not bevelled all round, a reflex
//! corner beside an edge left sharp, or a concave edge at any of those
//! vertices but one left sharp between a group's walls, is outside the
//! class too.

use ogeom_core::{FastSet, Tolerances};
use ogeom_math::{Point, Vector};

use super::classify::{Face, Planar, Seat, VertexClass};
use super::{Id, Origin};
use ogeom_core::FastMap;

/// A plane the blend puts beside consumed faces, by what makes it: its
/// outward unit normal and its offset along it.
pub(crate) type Cut = (Origin, Vector, f64);

/// Bevels round a consumed face's hole or at one of its reflex corners,
/// whose material is the union of the sides they leave; or one bevel of a
/// hole that is not convex, with the band it takes.
pub(crate) struct Group {
    /// The group's bevels' edges.
    pub edges: FastSet<Id>,
    /// The faces beside the group's edges other than the consumed one.
    pub walls: FastSet<usize>,
    /// Each of the group's bevel planes: its outward unit normal and its
    /// offset along it. For a band, the bevel's plane, then its wall's and
    /// the two planes that split it from the bevels beside it, each
    /// turned so the band stands past them all.
    pub planes: Vec<(Vector, f64)>,
}

/// One bevel's band round a hole that is not convex: what its bevel takes
/// between its wall and the planes splitting it from the bevels at its two
/// ends, so the bands of the hole's bevels fill what they take together,
/// none overlapping another. Each end gives the bevel there and whether
/// the material is convex at that corner (the hole turning the other way).
struct Band {
    edge: Id,
    wall: usize,
    ends: [(Id, bool); 2],
}

impl Group {
    /// Whether a ring from `origin` is one of the group's own.
    pub fn owns(&self, origin: &Origin) -> bool {
        match origin {
            Origin::Edge(e) => self.edges.contains(e),
            Origin::Face(f) => self.walls.contains(f),
            Origin::Vertex(_) => false,
        }
    }
}

/// The corners cut back by planes, and the planes.
pub(crate) struct Cluster {
    /// Each consumed face, by index.
    pub faces: FastSet<usize>,
    /// Every vertex cut back: the consumed faces' and the mitred corners'.
    pub vertices: FastSet<Id>,
    /// The vertices whose facet ring the clip replaces: those of fewer than
    /// three bevels.
    pub unfaceted: FastSet<Id>,
    /// Each face, bevel and facet plane round those vertices but a
    /// group's.
    pub cuts: Vec<Cut>,
    /// The groups whose union the rest keep, and their vertices, which
    /// keep their setback corners.
    pub groups: Vec<Group>,
    pub grouped: FastSet<Id>,
    /// Each blended edge's two ends.
    ends: FastMap<Id, (Id, Id)>,
    /// Each face's vertices, outer ring and holes.
    face_vertices: Vec<FastSet<Id>>,
    /// Per consumed face, its edges, its vertices and the faces beside it:
    /// the bevels, facets and faces that meet each other across it.
    across: Vec<(FastSet<Id>, FastSet<Id>, FastSet<usize>)>,
}

impl Cluster {
    /// Whether a ring from `ring` may lose its corner at vertex `at` to the
    /// cut from `cut`: a cut takes only what stands at its own element (a
    /// bevel at its edge's ends, a facet at its vertex, a face at its
    /// vertices), except that round one consumed face its bevels, its
    /// facets and the faces beside it cut each other anywhere, meeting
    /// across it (each bevel there fits its faces on its own, so none runs
    /// past them). A consumed face's outline, and a point with no vertex,
    /// may go to any cut.
    pub fn may_take(&self, ring: &Origin, at: Option<Id>, cut: &Origin) -> bool {
        let Some(v) = at else {
            return true;
        };
        if let Origin::Face(fi) = ring
            && self.faces.contains(fi)
        {
            return true;
        }
        let touches = match cut {
            Origin::Edge(e) => self.ends.get(e).is_some_and(|(a, b)| *a == v || *b == v),
            Origin::Vertex(w) => *w == v,
            Origin::Face(f) => self.face_vertices[*f].contains(&v),
        };
        if touches {
            return true;
        }
        let member =
            |o: &Origin, (edges, vertices, beside): &(FastSet<Id>, FastSet<Id>, FastSet<usize>)| {
                match o {
                    Origin::Edge(e) => edges.contains(e),
                    Origin::Vertex(w) => vertices.contains(w),
                    Origin::Face(f) => beside.contains(f),
                }
            };
        self.across
            .iter()
            .any(|f| member(ring, f) && member(cut, f))
    }
}

/// What the setbacks make of the corners.
pub(crate) enum Consumed {
    /// Every corner closes by its rings alone.
    Nothing,
    /// Corners to cut back, inside the class.
    Cluster(Box<Cluster>),
    /// Corners past the class: the engine declines.
    Outside,
}

/// A ring's normal by Newell's sum, twice its area long.
fn newell(points: &[Point]) -> Vector {
    let mut n = Vector::ZERO;
    for k in 0..points.len() {
        let (p, q) = (
            points[k].to_vector(),
            points[(k + 1) % points.len()].to_vector(),
        );
        n += p.cross(q);
    }
    n
}

/// Whether the face's ring `ring` turns its way at its `k`th corner.
fn turns(planar: &Planar, face: &Face, ring: &[Id], k: usize) -> bool {
    let n = ring.len();
    let (a, b, c) = (
        planar.points[&ring[(k + n - 1) % n]],
        planar.points[&ring[k]],
        planar.points[&ring[(k + 1) % n]],
    );
    (b - a).cross(c - b).dot(face.normal) > 0.0
}

/// Whether the setback rings of a face keep apart: each turning its way
/// (a hole the other way) without crossing itself, no two meeting, and
/// each hole inside the outer one.
fn apart(
    planar: &Planar,
    fi: usize,
    corners: &FastMap<(usize, Id), Point>,
    tol: Tolerances,
) -> bool {
    let face = &planar.faces[fi];
    let rings: Vec<Vec<Point>> = face
        .loops()
        .map(|ring| ring.iter().map(|v| corners[&(fi, *v)]).collect())
        .collect();
    if !fair(&rings[0], face.normal, tol) {
        return false;
    }
    if rings.len() == 1 {
        return true;
    }
    if rings[1..].iter().any(|h| !fair(h, -face.normal, tol)) {
        return false;
    }
    let o = rings[0][0];
    let x = (rings[0][1] - o).cross(face.normal);
    let x = x / x.magnitude();
    let y = face.normal.cross(x) / face.normal.magnitude();
    let flat: Vec<Vec<(f64, f64)>> = rings
        .iter()
        .map(|r| {
            r.iter()
                .map(|p| ((*p - o).dot(x), (*p - o).dot(y)))
                .collect()
        })
        .collect();
    let eps = tol.confusion();
    (0..flat.len())
        .all(|i| (i + 1..flat.len()).all(|j| !super::guards::loops_meet(&flat[i], &flat[j], eps)))
        && flat[1..]
            .iter()
            .all(|h| super::guards::strictly_inside(&flat[0], h[0], eps))
}

/// Whether the ring `points` turns about `normal` without crossing itself.
fn fair(points: &[Point], normal: Vector, tol: Tolerances) -> bool {
    if newell(points).dot(normal) <= tol.confusion() * tol.confusion() {
        return false;
    }
    let x = (points[1] - points[0]).cross(normal);
    let x = x / x.magnitude();
    let y = normal.cross(x) / normal.magnitude();
    let flat: Vec<(f64, f64)> = points
        .iter()
        .map(|p| ((*p - points[0]).dot(x), (*p - points[0]).dot(y)))
        .collect();
    !super::guards::crosses_itself(&flat, tol.confusion())
}

/// The corners of a request whose setback corners are `corners` that the
/// clip cuts back, and the planes it cuts by.
pub(crate) fn consumed(
    planar: &Planar,
    seats: &FastMap<(Id, Id), Seat>,
    blended: &FastSet<Id>,
    corners: &FastMap<(usize, Id), Point>,
    classes: &FastMap<Id, VertexClass>,
    scattered: &FastSet<Id>,
    tol: Tolerances,
) -> Consumed {
    let mut faces: FastSet<usize> = FastSet::default();
    for fi in 0..planar.faces.len() {
        if !apart(planar, fi, corners, tol) {
            faces.insert(fi);
        }
    }
    let edge_of = |a: Id, b: Id| planar.edges.get(&(a.min(b), a.max(b))).copied();
    // Each consumed face's groups: a hole, every corner of it reflex in the
    // face and every edge bevelled; a reflex corner of the outline, both
    // its edges bevelled.
    let mut eaten: Vec<usize> = faces.iter().copied().collect();
    eaten.sort_unstable();
    let mut grouped_edges: Vec<(Vec<Id>, Vec<Id>)> = Vec::new();
    let mut bands: Vec<(usize, Band)> = Vec::new();
    // Reflex corners with one edge sharp: the group, the bevelled edge,
    // the sharp one, the consumed face and the corner.
    let mut halves: Vec<(usize, Id, Id, usize, Id)> = Vec::new();
    let mut vertices: FastSet<Id> = FastSet::default();
    for fi in &eaten {
        let face = &planar.faces[*fi];
        let n = face.corners.len();
        for k in 0..n {
            if turns(planar, face, &face.corners, k) {
                continue;
            }
            let (a, v, b) = (
                face.corners[(k + n - 1) % n],
                face.corners[k],
                face.corners[(k + 1) % n],
            );
            let (Some(e0), Some(e1)) = (edge_of(a, v), edge_of(v, b)) else {
                return Consumed::Outside;
            };
            // One edge left sharp: the other's bevel runs into the corner
            // and stops on the sharp edge's wall plane, past which its
            // group takes nothing; its end, a facet at the corner, lies in
            // that wall.
            match (blended.contains(&e0), blended.contains(&e1)) {
                (true, false) => halves.push((grouped_edges.len(), e0, e1, *fi, v)),
                (false, true) => halves.push((grouped_edges.len(), e1, e0, *fi, v)),
                _ => {}
            }
            grouped_edges.push((vec![e0, e1], vec![v]));
        }
        for hole in &face.holes {
            let m = hole.len();
            let mut edges = Vec::with_capacity(m);
            for k in 0..m {
                let Some(e) = edge_of(hole[k], hole[(k + 1) % m]) else {
                    return Consumed::Outside;
                };
                edges.push(e);
            }
            // A convex hole is one group: its bevels take the funnel past
            // them all. One that is not goes as its bevels' bands, split
            // at each corner by the plane through the two bevels' meeting
            // line (the hole's straight skeleton, swept up the bevels).
            if (0..m).all(|k| !turns(planar, face, hole, k)) {
                grouped_edges.push((edges, hole.clone()));
                continue;
            }
            for k in 0..m {
                let Some(wall) = seats
                    .values()
                    .find(|s| s.edge == edges[k])
                    .and_then(|s| s.faces.into_iter().find(|f| *f != *fi))
                else {
                    return Consumed::Outside;
                };
                bands.push((
                    grouped_edges.len(),
                    Band {
                        edge: edges[k],
                        wall,
                        ends: [
                            (edges[(k + m - 1) % m], turns(planar, face, hole, k)),
                            (edges[(k + 1) % m], turns(planar, face, hole, (k + 1) % m)),
                        ],
                    },
                ));
                grouped_edges.push((vec![edges[k]], vec![hole[k], hole[(k + 1) % m]]));
            }
        }
        vertices.extend(face.loops().flatten().copied());
    }
    let sharp_ok: FastSet<Id> = halves.iter().map(|h| h.2).collect();
    if grouped_edges
        .iter()
        .flat_map(|(edges, _)| edges)
        .any(|e| !blended.contains(e) && !sharp_ok.contains(e))
    {
        return Consumed::Outside;
    }
    let grouped: FastSet<Id> = grouped_edges
        .iter()
        .flat_map(|(_, vs)| vs.iter().copied())
        .collect();
    // Every edge round the cluster convex, but a sharp one at a group's
    // vertex (a hole's corner, the outline's reflex corner).
    let convex_at = |v: Id| {
        seats
            .iter()
            .filter(|((a, b), _)| *a == v || *b == v)
            .all(|(_, seat)| seat.convex || (grouped.contains(&v) && !blended.contains(&seat.edge)))
    };
    if vertices.iter().any(|v| !convex_at(*v)) {
        return Consumed::Outside;
    }
    // A convex corner of two bevels whose facet ring spans an area: the two
    // setbacks differ, and the bevels mitre instead.
    let mut all: Vec<Id> = planar.points.keys().copied().collect();
    all.sort_unstable();
    for v in &all {
        if vertices.contains(v) || classes.get(v).is_none_or(|c| c.blended != 2) || !convex_at(*v) {
            continue;
        }
        let Some(order) = super::classify::faces_round(&planar.faces, *v) else {
            continue;
        };
        let ring: Vec<Point> = order.iter().map(|fi| corners[&(*fi, *v)]).collect();
        if newell(&ring).magnitude() > tol.confusion() * tol.confusion() {
            vertices.insert(*v);
        }
    }
    // A mitre of four bevels or more sharing no point: cut back, convex.
    for v in scattered {
        if !convex_at(*v) {
            return Consumed::Outside;
        }
        vertices.insert(*v);
    }
    if vertices.is_empty() {
        return Consumed::Nothing;
    }
    // Each bevel at a consumed face must fit its two faces on its own: a
    // contact line past the face's far side leaves the face it is seated
    // on, which is refused, not cut back.
    for seat in seats.values() {
        let at_cluster = |v: Id| faces.iter().any(|f| planar.faces[*f].corners.contains(&v));
        if !blended.contains(&seat.edge) || !(at_cluster(seat.ends.0) || at_cluster(seat.ends.1)) {
            continue;
        }
        for (k, fi) in seat.faces.iter().enumerate() {
            let face = &planar.faces[*fi];
            let inward = face
                .normal
                .cross(if k == 0 { seat.along } else { -seat.along });
            let from = planar.points[&seat.ends.0];
            let width = face
                .corners
                .iter()
                .map(|c| (planar.points[c] - from).dot(inward))
                .fold(0.0, f64::max);
            let Some(d) = corners_setback(planar, corners, seat, k) else {
                return Consumed::Outside;
            };
            if d > width + tol.confusion() {
                return Consumed::Outside;
            }
        }
    }
    let ended: FastSet<Id> = halves.iter().map(|h| h.4).collect();
    let unfaceted: FastSet<Id> = vertices
        .iter()
        .copied()
        .filter(|v| {
            !ended.contains(v)
                && (scattered.contains(v) || classes.get(v).is_some_and(|c| c.blended < 3))
        })
        .collect();
    let mut groups: Vec<Group> = grouped_edges
        .iter()
        .map(|(edges, _)| {
            let edges: FastSet<Id> = edges.iter().copied().collect();
            let walls = seats
                .values()
                .filter(|s| edges.contains(&s.edge))
                .flat_map(|s| s.faces)
                .filter(|f| !faces.contains(f))
                .collect();
            Group {
                edges,
                walls,
                planes: Vec::new(),
            }
        })
        .collect();
    let walls: FastSet<usize> = groups
        .iter()
        .flat_map(|g| g.walls.iter().copied())
        .collect();
    let mut cuts: Vec<Cut> = Vec::new();
    let mut bevels: FastMap<Id, (Vector, f64)> = FastMap::default();
    let mut keys: Vec<&(Id, Id)> = seats.keys().collect();
    keys.sort_unstable();
    for key in keys {
        let seat = &seats[key];
        if !blended.contains(&seat.edge)
            || !(vertices.contains(&key.0) || vertices.contains(&key.1))
        {
            continue;
        }
        let (a, b) = seat.ends;
        let (p, q, r) = (
            corners[&(seat.faces[0], a)],
            corners[&(seat.faces[0], b)],
            corners[&(seat.faces[1], a)],
        );
        let n = (q - p).cross(r - p);
        let out = seat.normals[0] + seat.normals[1];
        let n = if n.dot(out) < 0.0 { -n } else { n };
        let n = n / n.magnitude();
        let off = n.dot(p.to_vector());
        bevels.insert(seat.edge, (n, off));
        let mut in_group = false;
        for g in &mut groups {
            if g.edges.contains(&seat.edge) {
                g.planes.push((n, off));
                in_group = true;
            }
        }
        if !in_group {
            cuts.push((Origin::Edge(seat.edge), n, off));
        }
    }
    // Each band past its bevel, behind its wall, and on its own side of
    // the plane splitting it from each bevel beside it: where the material
    // is convex the band holds what its bevel takes deeper than the other
    // does, else what it takes less deep.
    for (gi, band) in &bands {
        let Some(&(n, off)) = bevels.get(&band.edge) else {
            return Consumed::Outside;
        };
        let wall = &planar.faces[band.wall];
        let at = planar.points[&wall.corners[0]];
        let mut planes = vec![(n, off), (-wall.normal, -wall.normal.dot(at.to_vector()))];
        for (beside, convex) in band.ends {
            let Some(&(m, o)) = bevels.get(&beside) else {
                return Consumed::Outside;
            };
            let (d, e) = if convex {
                (n - m, off - o)
            } else {
                (m - n, o - off)
            };
            let len = d.magnitude();
            if len <= tol.confusion() {
                return Consumed::Outside;
            }
            planes.push((d / len, e / len));
        }
        groups[*gi].planes = planes;
    }
    for (gi, bevelled, sharp, fi, _) in &halves {
        let Some(&(n, off)) = bevels.get(bevelled) else {
            return Consumed::Outside;
        };
        let Some(wall) = seats
            .values()
            .find(|s| s.edge == *sharp)
            .and_then(|s| s.faces.into_iter().find(|f| *f != *fi))
        else {
            return Consumed::Outside;
        };
        let wall = &planar.faces[wall];
        let at = planar.points[&wall.corners[0]];
        groups[*gi].planes = vec![(n, off), (wall.normal, wall.normal.dot(at.to_vector()))];
    }
    let mut round: Vec<Id> = vertices.iter().copied().collect();
    round.sort_unstable();
    let mut faced: FastSet<usize> = FastSet::default();
    for v in round {
        let Some(order) = super::classify::faces_round(&planar.faces, v) else {
            return Consumed::Outside;
        };
        for fi in &order {
            if !walls.contains(fi) && faced.insert(*fi) {
                let face = &planar.faces[*fi];
                let at = planar.points[&face.corners[0]];
                cuts.push((
                    Origin::Face(*fi),
                    face.normal,
                    face.normal.dot(at.to_vector()),
                ));
            }
        }
        // A facet ending a bevel on a sharp wall lies in that wall: no cut.
        if unfaceted.contains(&v) || ended.contains(&v) {
            continue;
        }
        let ring: Vec<Point> = order.iter().map(|fi| corners[&(*fi, v)]).collect();
        let n = newell(&ring);
        let m = n.magnitude();
        if m <= tol.confusion() * tol.confusion() {
            continue;
        }
        let n = n / m;
        let at = ring[0];
        let out = order
            .iter()
            .fold(Vector::ZERO, |s, fi| s + planar.faces[*fi].normal);
        let n = if n.dot(out) < 0.0 { -n } else { n };
        cuts.push((Origin::Vertex(v), n, n.dot(at.to_vector())));
    }
    let ends: FastMap<Id, (Id, Id)> = seats.values().map(|s| (s.edge, s.ends)).collect();
    let face_vertices: Vec<FastSet<Id>> = planar
        .faces
        .iter()
        .map(|f| f.loops().flatten().copied().collect())
        .collect();
    let mut eaten: Vec<usize> = faces.iter().copied().collect();
    eaten.sort_unstable();
    let across = eaten
        .iter()
        .map(|fi| {
            let face = &planar.faces[*fi];
            let pairs: Vec<(Id, Id)> = face
                .loops()
                .flat_map(|ring| {
                    let n = ring.len();
                    (0..n).map(move |k| {
                        let (a, b) = (ring[k], ring[(k + 1) % n]);
                        (a.min(b), a.max(b))
                    })
                })
                .collect();
            let edges = pairs
                .iter()
                .filter_map(|p| planar.edges.get(p).copied())
                .collect();
            let beside = pairs
                .iter()
                .filter_map(|p| seats.get(p))
                .map(|s| {
                    if s.faces[0] == *fi {
                        s.faces[1]
                    } else {
                        s.faces[0]
                    }
                })
                .collect();
            (edges, face.loops().flatten().copied().collect(), beside)
        })
        .collect();
    Consumed::Cluster(Box::new(Cluster {
        faces,
        vertices,
        unfaceted,
        cuts,
        groups,
        grouped,
        ends,
        face_vertices,
        across,
    }))
}

/// The setback of `seat`'s bevel on its `k`th face, read back from its
/// setback corners: the distance of the corner at its first end from the
/// edge's line, in that face. `None` where the corner is not there.
fn corners_setback(
    planar: &Planar,
    corners: &FastMap<(usize, Id), Point>,
    seat: &Seat,
    k: usize,
) -> Option<f64> {
    let face = &planar.faces[seat.faces[k]];
    let inward = face
        .normal
        .cross(if k == 0 { seat.along } else { -seat.along });
    let from = planar.points[&seat.ends.0];
    Some((*corners.get(&(seat.faces[k], seat.ends.0))? - from).dot(inward))
}
