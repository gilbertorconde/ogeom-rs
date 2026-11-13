//! Faces a request consumes: a band between two parallel bevels set back
//! past each other, as a thin wall chamfered on both its edges.
//!
//! The band's setback lines cross, so it keeps nothing: it goes, and its
//! two bevels run on until they meet each other along the ridge their
//! planes share. Each face across one of the band's sharp ends loses the
//! band's edge: its outline there runs up one bevel's trace to the ridge
//! and down the other's. The band's four vertices must each join three
//! edges, so each corner is where three known planes meet:
//!
//! - a band corner is the ridge's end, where the two bevels meet the face
//!   across the band's end;
//! - that face's corner beside it is where it meets the bevel and the
//!   plane its other edge there keeps (that edge's bevel where it has one,
//!   else the face across it).
//!
//! Any other consumed face stays outside the class: its rings fail P5's
//! guards and the engine declines.

use ogeom_core::{FastMap, FastSet, Tolerances};
use ogeom_math::{Point, Vector};

use super::classify::{Planar, Seat};
use super::{Id, id};

/// The bands a request consumes, and what each changes.
#[derive(Default)]
pub(crate) struct Consumed {
    /// Each consumed face, by index.
    pub faces: FastSet<usize>,
    /// Each ridge end, put into the face across a band's sharp end between
    /// that face's corners at the end's two vertices, in its turn: keyed by
    /// the face and the two vertices.
    pub ridges: FastMap<(usize, Id, Id), Point>,
    /// The band's sharp edges and its vertices: nothing of them is left.
    pub gone: Vec<Id>,
}

/// A plane as its unit normal and its offset along it.
type Flat = (Vector, f64);

/// Where three planes meet, `None` where two of them run parallel.
fn meet(a: Flat, b: Flat, c: Flat) -> Option<Point> {
    let det = a.0.dot(b.0.cross(c.0));
    if det.abs() <= 1e-12 {
        return None;
    }
    let v = (b.0.cross(c.0) * a.1 + c.0.cross(a.0) * b.1 + a.0.cross(b.0) * c.1) / det;
    Some(Point::ORIGIN + v)
}

/// The plane through three points, `None` where they are on one line.
fn through(p: Point, q: Point, r: Point) -> Option<Flat> {
    let n = (q - p).cross(r - p);
    let length = n.magnitude();
    if length <= 1e-12 {
        return None;
    }
    let n = n / length;
    Some((n, n.dot(p.to_vector())))
}

/// The consumed bands of a request, the setback corners moved where they
/// change: a band's corners to its ridge ends, and the corners of each face
/// across its sharp ends beside them to where that face meets its bevels.
/// `None` where a band borders another or a face across one is set back
/// itself on that edge: past the class.
pub(crate) fn consumed(
    planar: &Planar,
    seats: &FastMap<(Id, Id), Seat>,
    setback: &FastMap<(Id, Id), f64>,
    corners: &mut FastMap<(usize, Id), Point>,
    tol: Tolerances,
) -> Option<Consumed> {
    let edge = |a: Id, b: Id| planar.edges.get(&(a.min(b), a.max(b))).copied();
    let seat = |a: Id, b: Id| seats.get(&(a.min(b), a.max(b)));
    let across = |fi: usize, a: Id, b: Id| {
        let s = seat(a, b)?;
        Some(if s.faces[0] == fi {
            s.faces[1]
        } else {
            s.faces[0]
        })
    };
    let degree = |v: Id| seats.keys().filter(|(a, b)| *a == v || *b == v).count();
    let plane_of = |fi: usize| {
        let face = &planar.faces[fi];
        let n = face.normal;
        (n, n.dot(planar.points[&face.corners[0]].to_vector()))
    };
    let mut out = Consumed::default();
    let mut moves: Vec<((usize, Id), Point)> = Vec::new();
    let mut touched: FastSet<usize> = FastSet::default();
    for (fi, face) in planar.faces.iter().enumerate() {
        let ([_, _, _, _], []) = (face.corners.as_slice(), face.holes.as_slice()) else {
            continue;
        };
        let c = &face.corners;
        let back = |e: Id| setback.get(&(id(&face.shape), e)).copied();
        for k in 0..2 {
            let (v1, v2, v3, v4) = (c[k], c[k + 1], c[k + 2], c[(k + 3) % 4]);
            let (Some(a), Some(b), Some(cc), Some(d)) =
                (edge(v1, v2), edge(v2, v3), edge(v3, v4), edge(v4, v1))
            else {
                return None;
            };
            let (Some(da), None, Some(dc), None) = (back(a), back(b), back(cc), back(d)) else {
                continue;
            };
            let p = |v: Id| planar.points[&v];
            let ua = (p(v2) - p(v1)) / (p(v2) - p(v1)).magnitude();
            let uc = (p(v4) - p(v3)) / (p(v4) - p(v3)).magnitude();
            if ua.cross(uc).magnitude() > tol.angular() || ua.dot(uc) > 0.0 {
                continue;
            }
            let width = (p(v3) - p(v1)).dot(face.normal.cross(ua));
            if da + dc < width - tol.confusion() {
                continue;
            }
            // A band: past the class where a vertex joins more than three
            // edges.
            if [v1, v2, v3, v4].iter().any(|v| degree(*v) != 3) {
                return None;
            }
            let (Some(fa), Some(fc), Some(gb), Some(gd)) = (
                across(fi, v1, v2),
                across(fi, v3, v4),
                across(fi, v2, v3),
                across(fi, v4, v1),
            ) else {
                return None;
            };
            let at = |f: usize, v: Id| corners.get(&(f, v)).copied();
            let (qa, qc) = (
                through(at(fi, v1)?, at(fi, v2)?, at(fa, v1)?)?,
                through(at(fi, v3)?, at(fi, v4)?, at(fc, v3)?)?,
            );
            let (rb, rd) = (meet(qa, qc, plane_of(gb))?, meet(qa, qc, plane_of(gd))?);
            // The plane the other edge at `v` of the face across an end
            // keeps: its bevel, or the face across it.
            let kept = |g: usize, v: Id, end: Id| -> Option<Flat> {
                let gface = &planar.faces[g];
                let (prev, next) = gface.beside(v)?;
                let other = if prev == end { next } else { prev };
                let e = edge(v, other)?;
                match setback.get(&(id(&gface.shape), e)) {
                    Some(_) => {
                        let far = across(g, v, other)?;
                        through(at(g, v)?, at(g, other)?, at(far, v)?)
                    }
                    None => Some(plane_of(across(g, v, other)?)),
                }
            };
            moves.extend([
                ((fi, v1), rd),
                ((fi, v2), rb),
                ((fi, v3), rb),
                ((fi, v4), rd),
                ((gb, v2), meet(plane_of(gb), qa, kept(gb, v2, v3)?)?),
                ((gb, v3), meet(plane_of(gb), qc, kept(gb, v3, v2)?)?),
                ((gd, v1), meet(plane_of(gd), qa, kept(gd, v1, v4)?)?),
                ((gd, v4), meet(plane_of(gd), qc, kept(gd, v4, v1)?)?),
            ]);
            // Each face across an end walks the end the other way round.
            out.ridges.insert((gb, v3, v2), rb);
            out.ridges.insert((gd, v1, v4), rd);
            out.faces.insert(fi);
            out.gone.extend([b, d, v1, v2, v3, v4]);
            touched.extend([fa, fc, gb, gd]);
            break;
        }
    }
    if out.faces.iter().any(|f| touched.contains(f)) {
        return None;
    }
    for (key, point) in moves {
        corners.insert(key, point);
    }
    out.gone.sort_unstable();
    out.gone.dedup();
    Some(out)
}
