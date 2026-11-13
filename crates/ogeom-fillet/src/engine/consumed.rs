//! Corners the setback rings cannot close alone, cut back by the planes
//! the blend puts there: faces a request consumes, whole or in part (a
//! thin wall bevelled on both its edges, a plate bevelled all round past
//! half its thickness, a face whose bevels cross on it at an angle), and a
//! convex corner where two bevels of unequal setbacks meet beside a sharp
//! edge, which mitres.
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
//! The clip is local by construction: a ring may lose only the setback
//! corners at those vertices, a consumed face's own corners, a bevel's
//! run-on, or points the clip itself made. A clip that would take any
//! other corner reaches past them, and the engine declines. A consumed face
//! with holes or reflex corners, or a concave edge at any of those
//! vertices, is outside the class too.

use ogeom_core::{FastSet, Tolerances};
use ogeom_math::{Point, Vector};

use super::classify::{Face, Planar, Seat, VertexClass};
use super::{Id, Origin};
use ogeom_core::FastMap;

/// A plane the blend puts beside consumed faces, by what makes it: its
/// outward unit normal and its offset along it.
pub(crate) type Cut = (Origin, Vector, f64);

/// The corners cut back by planes, and the planes.
pub(crate) struct Cluster {
    /// Each consumed face, by index.
    pub faces: FastSet<usize>,
    /// Every vertex cut back: the consumed faces' and the mitred corners'.
    pub vertices: FastSet<Id>,
    /// The vertices whose facet ring the clip replaces: those of fewer than
    /// three bevels.
    pub unfaceted: FastSet<Id>,
    /// Each face, bevel and facet plane round those vertices.
    pub cuts: Vec<Cut>,
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

/// Whether the face's outline turns its way at every corner.
fn convex(planar: &Planar, face: &Face) -> bool {
    let n = face.corners.len();
    (0..n).all(|k| {
        let (a, b, c) = (
            planar.points[&face.corners[(k + n - 1) % n]],
            planar.points[&face.corners[k]],
            planar.points[&face.corners[(k + 1) % n]],
        );
        (b - a).cross(c - b).dot(face.normal) > 0.0
    })
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
    tol: Tolerances,
) -> Consumed {
    let mut faces: FastSet<usize> = FastSet::default();
    for (fi, face) in planar.faces.iter().enumerate() {
        let ring: Vec<Point> = face.corners.iter().map(|v| corners[&(fi, *v)]).collect();
        if !fair(&ring, face.normal, tol) {
            faces.insert(fi);
        }
    }
    let convex_at = |v: Id| {
        seats
            .iter()
            .filter(|((a, b), _)| *a == v || *b == v)
            .all(|(_, seat)| seat.convex)
    };
    let mut vertices: FastSet<Id> = FastSet::default();
    for fi in &faces {
        let face = &planar.faces[*fi];
        if !face.holes.is_empty() || !convex(planar, face) {
            return Consumed::Outside;
        }
        vertices.extend(face.corners.iter().copied());
    }
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
    let unfaceted: FastSet<Id> = vertices
        .iter()
        .copied()
        .filter(|v| classes.get(v).is_some_and(|c| c.blended < 3))
        .collect();
    let mut cuts: Vec<Cut> = Vec::new();
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
        cuts.push((Origin::Edge(seat.edge), n, n.dot(p.to_vector())));
    }
    let mut round: Vec<Id> = vertices.iter().copied().collect();
    round.sort_unstable();
    let mut faced: FastSet<usize> = FastSet::default();
    for v in round {
        let Some(order) = super::classify::faces_round(&planar.faces, v) else {
            return Consumed::Outside;
        };
        for fi in &order {
            if faced.insert(*fi) {
                let face = &planar.faces[*fi];
                let at = planar.points[&face.corners[0]];
                cuts.push((
                    Origin::Face(*fi),
                    face.normal,
                    face.normal.dot(at.to_vector()),
                ));
            }
        }
        if unfaceted.contains(&v) {
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
            let n = face.corners.len();
            let pairs: Vec<(Id, Id)> = (0..n)
                .map(|k| {
                    let (a, b) = (face.corners[k], face.corners[(k + 1) % n]);
                    (a.min(b), a.max(b))
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
            (edges, face.corners.iter().copied().collect(), beside)
        })
        .collect();
    Consumed::Cluster(Box::new(Cluster {
        faces,
        vertices,
        unfaceted,
        cuts,
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
