//! Chamfers of solids with curved faces: P1's contacts on curved seats, P3's
//! vertex cases there, and P4's trimmed-face assembly, over the solid
//! [`curved`](super::curved) reads and the corners curved P2
//! ([`curved_corners`](super::setbacks::curved_corners)) finds.
//!
//! The seats:
//!
//! - a straight edge between two planes, each contact the edge moved in
//!   along its face;
//! - a straight edge between a plane and a cylinder it runs along, the
//!   plane's contact moved in along it and the cylinder's the ruling the
//!   distance round the cylinder from the edge, measured along it;
//! - a rim between a planar cap and a coaxial cylinder or cone, each
//!   contact the ring the distance in along the cap and down the wall.
//!
//! Each bevel is the plane through two straight contacts, or the cone
//! through two rings. Every face is rebuilt on its own surface, its
//! boundary each of its edges' pieces between the corners at their ends: a
//! blended edge's piece along the face's contact, an edge left sharp along
//! itself. A sharp edge's piece ends where a blend at its vertex moved the
//! corner on it.
//!
//! The vertex cases built:
//!
//! - one bevel, its end lying on the third face, a plane, there: that face
//!   runs round the bevel's end (a flush end);
//! - one bevel and no third face, as at a seam's foot: the bevel's end is
//!   its own seam;
//! - two bevels running on into each other along tangent edges: they share
//!   their end.
//!
//! Any other vertex, seat or section declines, and the caller keeps its own
//! path; so does a result failing its check, a bevel not proved within the
//! material its two faces bound (P5's certify, a plane's half-space and a
//! cylinder's or cone's wall), or a result not holding less than the solid
//! (every bevel here cuts a convex edge).

use ogeom_algo::{Built, History};
use ogeom_core::{FastMap, FastSet, OgeomResult, Tolerances};
use ogeom_geom::{CircleCurve, Curve, Curve3d as _, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, Shape};

use super::curved::{Curved, CurvedEdge, SheetCurve, SheetRings, build, read};
use super::setbacks::{Corner, curved_corners};
use super::{Id, id};
use crate::chamfer::{Chamfer, Legs, legs, rim_bevel_surface, wall_depth};

/// One blended edge on a curved solid, as a section law hands it to the
/// assembly: its two faces, each one's contact, the strip's surface, and
/// its section, which closes each end of the strip between the two faces'
/// corners there.
#[derive(Debug, Clone)]
pub(crate) struct Strip {
    pub edge: Id,
    pub faces: [usize; 2],
    pub contacts: [Curve; 2],
    pub surface: SurfaceGeometry,
    pub section: Section,
}

/// A strip's section, as it closes an end of the strip.
#[derive(Debug, Clone)]
pub(crate) enum Section {
    /// A chamfer's: the segment between the two corners.
    Line,
    /// A round's: the arc between the two corners about the spine's point
    /// at their station, the shorter way round. `spine` is a bounded curve
    /// along the strip.
    #[cfg_attr(
        not(test),
        allow(
            dead_code,
            reason = "the round sections on curved seats build it from M3 on"
        )
    )]
    Arc { spine: Curve },
}

/// The chamfer of `specs` on `solid` with curved faces, or `None` where a
/// seat, a vertex or the result lies outside the classes built here.
///
/// # Errors
///
/// As building the sheet and checking the solid.
pub(crate) fn chamfer(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    // P1: each blended edge's contacts and bevel, in key order.
    let mut asked: Vec<(Id, &Shape, &Chamfer)> = specs
        .iter()
        .map(|(edge, spec)| (id(edge), edge, spec))
        .collect();
    asked.sort_unstable_by_key(|(e, _, _)| *e);
    if asked.windows(2).any(|w| w[0].0 == w[1].0) {
        return Ok(None);
    }
    let mut bevels: Vec<Strip> = Vec::with_capacity(asked.len());
    for (edge, shape, spec) in &asked {
        let Some(bevel) = bevel_of(model, solid, &curved, *edge, shape, spec, tol)? else {
            return Ok(None);
        };
        bevels.push(bevel);
    }
    let mut contacts: FastMap<(usize, Id), Curve> = FastMap::default();
    for b in &bevels {
        for k in 0..2 {
            contacts.insert((b.faces[k], b.edge), b.contacts[k].clone());
        }
    }
    // P2.
    let Some(corners) = curved_corners(&curved, &contacts, tol)? else {
        return Ok(None);
    };
    let blended: FastSet<Id> = bevels.iter().map(|b| b.edge).collect();
    // P4 and P5. A face the bevels consume (two rims' bevels crossing on a
    // thin wall) is outside the cases built here: its assembly, or the
    // solid's check, fails, and the caller keeps its own path.
    let built = match assemble(model, &curved, &bevels, &blended, &corners, tol) {
        Ok(Some(built)) => built,
        Ok(None) => return Ok(None),
        Err(e @ ogeom_core::OgeomError::Cancelled) => return Err(e),
        Err(_) => return Ok(None),
    };
    match super::guards::solid_fits(model, &built.shape, tol) {
        Ok(()) => {}
        Err(e @ ogeom_core::OgeomError::Cancelled) => return Err(e),
        Err(_) => return Ok(None),
    }
    // Every bevel here cuts a convex edge: the result holds less than the
    // solid, and more than nothing.
    let volume = |model: &Model, s: &Shape| {
        ogeom_algo::volume_properties(model, s, ogeom_mesh::Deflection::default(), tol)
            .map(|v| v.mass)
    };
    let (before, after) = (volume(model, solid)?, volume(model, &built.shape)?);
    if !(after > 0.0 && after < before) {
        return Ok(None);
    }
    Ok(Some(built))
}

/// The faces of `curved` walking `edge`, each once.
fn faces_of(curved: &Curved, edge: Id) -> Vec<usize> {
    (0..curved.faces.len())
        .filter(|fi| {
            curved.faces[*fi]
                .loops
                .iter()
                .flatten()
                .any(|(e, _)| *e == edge)
        })
        .collect()
}

/// A face's outward unit normal at `p`, on a plane or a cylinder; `None` on
/// any other surface.
fn outward(curved: &Curved, face: usize, p: Point) -> Option<Vector> {
    let f = &curved.faces[face];
    let sign = if f.reversed { -1.0 } else { 1.0 };
    match &f.surface {
        SurfaceGeometry::Plane(plane) => Some(plane.plane().normal().vector() * sign),
        SurfaceGeometry::Cylinder(c) => {
            let axis = c.cylinder().axis();
            let (o, a) = (axis.location, axis.direction.vector());
            let r = (p - o) - a * (p - o).dot(a);
            let m = r.magnitude();
            (m > 0.0).then(|| r / m * sign)
        }
        _ => None,
    }
}

/// A point of face `face` away from edge `edge`: a vertex of its outline
/// not on the edge. `None` where it has none.
fn away(curved: &Curved, face: usize, edge: Id) -> Option<Point> {
    let ends = curved.edges[&edge].ends;
    curved.faces[face]
        .loops
        .iter()
        .flatten()
        .flat_map(|(e, _)| {
            let c = &curved.edges[e];
            [c.ends.0, c.ends.1]
        })
        .find(|v| *v != ends.0 && *v != ends.1)
        .map(|v| curved.points[&v])
}

/// The two distances of `spec` on `edge`'s faces in `faces`' order, `None`
/// for an angle, whose second distance [`angled`] reads from the section.
fn distances(curved: &Curved, faces: [usize; 2], spec: &Chamfer) -> Option<[f64; 2]> {
    match spec {
        Chamfer::Symmetric(d) => Some([*d, *d]),
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => {
            let k = faces
                .iter()
                .position(|f| curved.faces[*f].shape.is_same(face))?;
            let mut d = [*on_other; 2];
            d[k] = *on_face;
            Some(d)
        }
        Chamfer::Angle { .. } => None,
    }
}

/// The contacts and bevel of `spec` on `edge`, `None` past the seats built
/// here.
fn bevel_of(
    model: &Model,
    solid: &Shape,
    curved: &Curved,
    edge: Id,
    shape: &Shape,
    spec: &Chamfer,
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
    let Some(e) = curved.edges.get(&edge) else {
        return Ok(None);
    };
    let [f0, f1] = faces_of(curved, edge)[..] else {
        return Ok(None);
    };
    let faces = [f0, f1];
    match &e.curve {
        Curve::Line(_) => {
            let d = match spec {
                Chamfer::Angle {
                    face,
                    distance,
                    angle,
                } => {
                    let Some(k) = faces
                        .iter()
                        .position(|f| curved.faces[*f].shape.is_same(face))
                    else {
                        return Ok(None);
                    };
                    let Some(other) = angled(curved, edge, e, faces, k, *distance, *angle, tol)
                    else {
                        return Ok(None);
                    };
                    let mut d = [other; 2];
                    d[k] = *distance;
                    d
                }
                _ => {
                    let Some(d) = distances(curved, faces, spec) else {
                        return Ok(None);
                    };
                    d
                }
            };
            if d.iter().any(|x| !(x.is_finite() && *x > tol.confusion())) {
                return Ok(None);
            }
            straight(curved, edge, e, faces, d, tol)
        }
        Curve::Circle(_) => rim(model, solid, curved, edge, shape, faces, spec, tol),
        _ => Ok(None),
    }
}

/// A distance-angle chamfer's second distance on a straight seat: in the
/// edge's section, the bevel leaves face `k`'s contact, `distance` in from
/// the edge, at `angle` from that face into the material, and runs to the
/// other face's trace, a line or the cylinder's circle; the distance there
/// is measured along it from the edge. `None` where the bevel does not
/// reach the other face, or a face is not a plane or a cylinder along the
/// edge.
#[allow(clippy::too_many_arguments, reason = "the seat's whole section")]
fn angled(
    curved: &Curved,
    edge: Id,
    e: &CurvedEdge,
    faces: [usize; 2],
    k: usize,
    distance: f64,
    angle: f64,
    tol: Tolerances,
) -> Option<f64> {
    if !(angle.is_finite() && angle > tol.angular() && angle < std::f64::consts::PI)
        || !(distance.is_finite() && distance > tol.confusion())
    {
        return None;
    }
    let (a, b) = (curved.points[&e.ends.0], curved.points[&e.ends.1]);
    let u = (b - a) / a.distance(b);
    let at = a + (b - a) * 0.5;
    let normals = [
        outward(curved, faces[0], at)?,
        outward(curved, faces[1], at)?,
    ];
    // Each face's trace in the section: a line from the edge, or a circle.
    enum Trace {
        Line(Vector),
        Circle(Point, f64, Vector),
    }
    let trace = |j: usize| -> Option<Trace> {
        match &curved.faces[faces[j]].surface {
            SurfaceGeometry::Plane(_) => {
                let m = normals[j].cross(u);
                let inside = away(curved, faces[j], edge)?;
                let m = if (inside - at).dot(m) < 0.0 { -m } else { m };
                Some(Trace::Line(m / m.magnitude()))
            }
            SurfaceGeometry::Cylinder(c) => {
                let cyl = c.cylinder();
                let axis = cyl.axis();
                let (o, ax) = (axis.location, axis.direction.vector());
                if u.cross(ax).magnitude() > tol.angular() {
                    return None;
                }
                let centre = o + ax * (at - o).dot(ax);
                Some(Trace::Circle(centre, cyl.radius(), ax))
            }
            _ => None,
        }
    };
    let (named, other) = (trace(k)?, trace(1 - k)?);
    let into = |p: Point| -> Option<Vector> { Some(-outward(curved, faces[k], p)?) };
    // Where the bevel leaves the named face, and the way back along it.
    let (start, back) = match named {
        Trace::Line(m) => (at + m * distance, -m),
        Trace::Circle(centre, radius, ax) => {
            let r = at - centre;
            let turn = |sign: f64| {
                let w = ax.cross(r);
                centre + r * (distance / radius).cos() + w * (sign * (distance / radius).sin())
            };
            let n_other = normals[1 - k];
            let (plus, minus) = (turn(1.0), turn(-1.0));
            let p = if (plus - at).dot(n_other) < (minus - at).dot(n_other) {
                plus
            } else {
                minus
            };
            let w = ax.cross(p - centre);
            let w = w / w.magnitude();
            (p, if w.dot(at - p) > 0.0 { w } else { -w })
        }
    };
    let t = back * angle.cos() + into(start)? * angle.sin();
    match other {
        Trace::Line(m) => {
            let n = normals[1 - k];
            let s = n.dot(at - start) / n.dot(t);
            let hit = start + t * s;
            let lambda = (hit - at).dot(m);
            (s > 0.0 && lambda > tol.confusion()).then_some(lambda)
        }
        Trace::Circle(centre, radius, _) => {
            let w = start - centre;
            let (qa, qb, qc) = (t.dot(t), 2.0 * w.dot(t), w.dot(w) - radius * radius);
            let disc = qb * qb - 4.0 * qa * qc;
            if disc < 0.0 {
                return None;
            }
            let root = disc.sqrt();
            let s = [(-qb - root) / (2.0 * qa), (-qb + root) / (2.0 * qa)]
                .into_iter()
                .filter(|s| *s > tol.confusion())
                .fold(f64::INFINITY, f64::min);
            if !s.is_finite() {
                return None;
            }
            let hit = start + t * s;
            let (x, y) = (at - centre, hit - centre);
            let arc = radius
                * (x.dot(y) / (x.magnitude() * y.magnitude()))
                    .clamp(-1.0, 1.0)
                    .acos();
            (arc > tol.confusion()).then_some(arc)
        }
    }
}

/// A straight edge between two planes, or a plane and a cylinder it runs
/// along: each contact a segment the length of the edge and half again past
/// each end, so the corners at its ends lie on it.
fn straight(
    curved: &Curved,
    edge: Id,
    e: &CurvedEdge,
    faces: [usize; 2],
    d: [f64; 2],
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
    let (a, b) = (curved.points[&e.ends.0], curved.points[&e.ends.1]);
    let length = a.distance(b);
    if length <= tol.confusion() {
        return Ok(None);
    }
    let u = (b - a) / length;
    let run = u * (length / 2.0);
    let mid = a + (b - a) * 0.5;
    let (Some(n0), Some(n1)) = (
        outward(curved, faces[0], mid),
        outward(curved, faces[1], mid),
    ) else {
        return Ok(None);
    };
    let normals = [n0, n1];
    let mut lines: Vec<(Point, Point)> = Vec::with_capacity(2);
    for k in 0..2 {
        let (fi, other) = (faces[k], normals[1 - k]);
        let face = &curved.faces[fi];
        let shift: Vector = match &face.surface {
            SurfaceGeometry::Plane(_) => {
                let m = normals[k].cross(u);
                let Some(inside) = away(curved, fi, edge) else {
                    return Ok(None);
                };
                let m = if (inside - mid).dot(m) < 0.0 { -m } else { m };
                // A convex edge only: into the face goes under the other.
                if m.dot(other) >= 0.0 {
                    return Ok(None);
                }
                m / m.magnitude() * d[k]
            }
            SurfaceGeometry::Cylinder(c) => {
                let cyl = c.cylinder();
                let axis = cyl.axis();
                let (o, ax) = (axis.location, axis.direction.vector());
                if u.cross(ax).magnitude() > tol.angular() {
                    return Ok(None);
                }
                let radius = cyl.radius();
                let r = (mid - o) - ax * (mid - o).dot(ax);
                let angle = d[k] / radius;
                if angle >= std::f64::consts::FRAC_PI_2 {
                    return Ok(None);
                }
                // Round the axis, the way into the face: away from the other
                // face's outward side.
                let turn = |sign: f64| {
                    let w = ax.cross(r);
                    r * angle.cos() + w * (sign * angle.sin()) - r
                };
                let (plus, minus) = (turn(1.0), turn(-1.0));
                let shift = if plus.dot(other) < minus.dot(other) {
                    plus
                } else {
                    minus
                };
                if shift.dot(other) >= 0.0 {
                    return Ok(None);
                }
                shift
            }
            _ => return Ok(None),
        };
        lines.push((a + shift - run, b + shift + run));
    }
    let contacts = [
        Curve::from(LineCurve::segment(lines[0].0, lines[0].1, tol)?),
        Curve::from(LineCurve::segment(lines[1].0, lines[1].1, tol)?),
    ];
    // The bevel: the plane through both contacts.
    let across = lines[1].0 - lines[0].0;
    let n = u.cross(across);
    if n.magnitude() <= tol.confusion() {
        return Ok(None);
    }
    let frame = Frame::new(
        lines[0].0,
        Direction::new(n, tol)?,
        Direction::new(u, tol)?,
        tol,
    )?;
    Ok(Some(Strip {
        edge,
        faces,
        contacts,
        surface: PlaneSurface::new(Plane::new(frame)).into(),
        section: Section::Line,
    }))
}

/// A rim between a planar cap and a coaxial cylinder or cone, as the wedge
/// path reads it ([`legs`]): the ring along the cap and the ring down the
/// wall, and the cone through both.
#[allow(clippy::too_many_arguments, reason = "the seat's whole reading")]
fn rim(
    model: &Model,
    solid: &Shape,
    curved: &Curved,
    edge: Id,
    shape: &Shape,
    faces: [usize; 2],
    spec: &Chamfer,
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
    let Ok(Legs::Rim(seat, on_wall, on_cap)) = legs(model, solid, shape, spec, tol) else {
        return Ok(None);
    };
    if seat.additive() {
        return Ok(None);
    }
    let cap = faces
        .iter()
        .position(|f| curved.faces[*f].shape.is_same(&seat.cap_face));
    let wall = faces
        .iter()
        .position(|f| curved.faces[*f].shape.is_same(&seat.wall_face));
    let (Some(cap), Some(wall)) = (cap, wall) else {
        return Ok(None);
    };
    let cap_rho = seat.sigma.mul_add(-(seat.tau * on_cap), seat.radius);
    let depth = wall_depth(&seat, on_wall);
    let wall_rho = seat.radius_at(-seat.tau * depth);
    if cap_rho <= tol.confusion() || wall_rho <= tol.confusion() {
        return Ok(None);
    }
    let ring = |centre: Point, radius: f64| -> OgeomResult<Curve> {
        Ok(CircleCurve::new(Circle::new(seat.frame_at(centre, tol)?, radius, tol)?).into())
    };
    let mut contacts = [ring(seat.centre, cap_rho)?, ring(seat.centre, cap_rho)?];
    contacts[cap] = ring(seat.centre, cap_rho)?;
    contacts[wall] = ring(seat.centre + seat.up * (-seat.tau * depth), wall_rho)?;
    let surface = rim_bevel_surface(&seat, depth, cap_rho, wall_rho, tol)?;
    Ok(Some(Strip {
        edge,
        faces,
        contacts,
        surface,
        section: Section::Line,
    }))
}

/// The sheet's points, each held once within the confusion distance.
struct Points {
    points: Vec<Point>,
    tol: Tolerances,
}

impl Points {
    fn index(&mut self, p: Point) -> usize {
        if let Some(i) = self
            .points
            .iter()
            .position(|q| q.distance(p) <= self.tol.confusion())
        {
            return i;
        }
        self.points.push(p);
        self.points.len() - 1
    }
}

/// A piece of `curve` from parameter `from` to `to`, running the way the
/// curve's parameter does when `rising`: its stored range (always
/// increasing) and whether the walk runs against it. A closed curve's piece
/// from a parameter to itself is its whole turn.
fn piece(curve: &Curve, from: f64, to: f64, rising: bool, tol: Tolerances) -> ((f64, f64), bool) {
    let period = if curve.is_periodic() || curve.is_closed(tol) {
        let (lo, hi) = curve.domain();
        Some(hi - lo)
    } else {
        None
    };
    let same = (from - to).abs() <= tol.parametric();
    match (rising, period) {
        (true, Some(p)) => {
            let mut to = to;
            while to <= from + if same { 0.5 * p } else { 0.0 } {
                to += p;
            }
            ((from, to), false)
        }
        (false, Some(p)) => {
            let mut from = from;
            while from <= to + if same { 0.5 * p } else { 0.0 } {
                from += p;
            }
            ((to, from), true)
        }
        (true, None) => ((from, to), false),
        (false, None) => ((to, from), true),
    }
}

/// Whether `curve` runs at parameter `t` the way `along` points.
fn runs(curve: &Curve, t: f64, along: Vector, tol: Tolerances) -> OgeomResult<bool> {
    Ok(curve.d1_at(t, tol)?.dot(along) > 0.0)
}

/// The sheet's edges, each made once and keyed by what it is.
#[derive(PartialEq, Eq, Hash, Clone, Copy)]
enum Key {
    /// A blended edge's contact on a face.
    Contact(usize, Id),
    /// An edge left sharp, between the corners at its ends.
    Sharp(Id),
    /// A bevel's end, a segment between two of the sheet's points.
    Section(usize, usize),
}

/// P3, P4 and P5 for strips on a solid with curved faces, whatever section
/// law made them: the sheet of the trimmed faces and the strips, its shell
/// and solid with its seams joined, each strip proved within the material
/// its two faces bound, and the history; `None` at a vertex outside the
/// cases built, or a strip found past its material.
///
/// # Errors
///
/// As building the sheet.
pub(crate) fn assemble(
    model: &mut Model,
    curved: &Curved,
    bevels: &[Strip],
    blended: &FastSet<Id>,
    corners: &FastMap<(usize, Id), Corner>,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let mut points = Points {
        points: Vec::new(),
        tol,
    };
    let mut edges: Vec<SheetCurve> = Vec::new();
    let mut made: FastMap<Key, usize> = FastMap::default();
    let contact_of = |fi: usize, e: Id| -> Option<&Curve> {
        bevels.iter().find(|b| b.edge == e).and_then(|b| {
            b.faces
                .iter()
                .position(|f| *f == fi)
                .map(|k| &b.contacts[k])
        })
    };
    // Each sharp edge's end at each of its vertices: the corner a blend
    // moved there on it, or the vertex itself.
    let mut sharp_end: FastMap<(Id, Id), Point> = FastMap::default();
    for (e, edge) in &curved.edges {
        if blended.contains(e) {
            continue;
        }
        for v in [edge.ends.0, edge.ends.1] {
            let at = curved.points[&v];
            let mut moved: Vec<Point> = Vec::new();
            for fi in faces_of(curved, *e) {
                if let Some(c) = corners.get(&(fi, v))
                    && c.point.distance(at) > tol.confusion()
                    && c.on(*e).is_some()
                    && !moved.iter().any(|m| m.distance(c.point) <= tol.confusion())
                {
                    moved.push(c.point);
                }
            }
            match moved.as_slice() {
                [] => {
                    sharp_end.insert((*e, v), at);
                }
                [p] => {
                    sharp_end.insert((*e, v), *p);
                }
                _ => return Ok(None),
            }
        }
    }
    // A vertex's blends: one, or two running on into each other; each
    // bevel's end there, the segment between its two faces' corners.
    let mut sections: FastMap<Id, Vec<(Point, Point, usize)>> = FastMap::default();
    for (bi, b) in bevels.iter().enumerate() {
        let e = &curved.edges[&b.edge];
        for v in [e.ends.0, e.ends.1] {
            let (Some(c0), Some(c1)) =
                (corners.get(&(b.faces[0], v)), corners.get(&(b.faces[1], v)))
            else {
                return Ok(None);
            };
            if c0.point.distance(c1.point) > tol.confusion() {
                sections
                    .entry(v)
                    .or_default()
                    .push((c0.point, c1.point, bi));
            }
        }
    }
    let mut keyed: Vec<&Id> = curved.points.keys().collect();
    keyed.sort_unstable();
    for v in keyed {
        let mut at_v: Vec<Id> = curved
            .edges
            .iter()
            .filter(|(e, c)| blended.contains(*e) && (c.ends.0 == *v || c.ends.1 == *v))
            .map(|(e, _)| *e)
            .collect();
        at_v.sort_unstable();
        at_v.dedup();
        match at_v.len() {
            0 | 1 => {}
            2 => {
                // Two bevels running on into each other share their end.
                if let Some(ends) = sections.get(v) {
                    let first = ends[0];
                    let same = |(p, q, _): &(Point, Point, usize)| {
                        (p.distance(first.0) <= tol.confusion()
                            && q.distance(first.1) <= tol.confusion())
                            || (p.distance(first.1) <= tol.confusion()
                                && q.distance(first.0) <= tol.confusion())
                    };
                    if !ends.iter().all(same) {
                        return Ok(None);
                    }
                }
            }
            _ => return Ok(None),
        }
    }
    let section = |points: &mut Points,
                   edges: &mut Vec<SheetCurve>,
                   made: &mut FastMap<Key, usize>,
                   kind: &Section,
                   p: Point,
                   q: Point|
     -> OgeomResult<(usize, bool)> {
        let (i, j) = (points.index(p), points.index(q));
        let key = Key::Section(i.min(j), i.max(j));
        if let Some(k) = made.get(&key) {
            return Ok((*k, edges[*k].ends.0 == i));
        }
        let (curve, range, ends): (Curve, (f64, f64), (usize, usize)) = match kind {
            Section::Line => {
                let (a, b) = (points.points[i.min(j)], points.points[i.max(j)]);
                let line = LineCurve::segment(a, b, tol)?;
                let range = line.domain();
                (line.into(), range, (i.min(j), i.max(j)))
            }
            Section::Arc { spine } => {
                let (a, b) = (points.points[i], points.points[j]);
                let centre = spine.point_at(param_on(spine, spine.domain(), a, tol)?, tol)?;
                let (x, w) = (a - centre, b - centre);
                let normal = x.cross(w);
                let frame = Frame::new(
                    centre,
                    Direction::new(normal, tol)?,
                    Direction::new(x, tol)?,
                    tol,
                )?;
                let y = normal.cross(x);
                let angle = w.dot(y).atan2(w.dot(x) * normal.magnitude());
                let circle = CircleCurve::new(Circle::new(frame, x.magnitude(), tol)?);
                (circle.into(), (0.0, angle), (i, j))
            }
        };
        edges.push(SheetCurve { ends, curve, range });
        made.insert(key, edges.len() - 1);
        Ok((edges.len() - 1, ends.0 == i))
    };
    // P4: each face rebuilt, its loops in the stored order, each edge's
    // piece between the corners at its ends, walked the surface's way.
    let mut rings_of: Vec<SheetRings> = Vec::with_capacity(curved.faces.len() + bevels.len());
    for (fi, face) in curved.faces.iter().enumerate() {
        let mut rings = Vec::with_capacity(face.loops.len());
        for walk in &face.loops {
            let mut ring: Vec<(usize, bool)> = Vec::new();
            let mut last: Option<Point> = None;
            let mut first: Option<Point> = None;
            for (e, w) in walk {
                let forward = *w != face.reversed;
                let edge = &curved.edges[e];
                let (vs, ve) = if forward {
                    (edge.ends.0, edge.ends.1)
                } else {
                    (edge.ends.1, edge.ends.0)
                };
                // This edge's piece on this face: its contact, or itself.
                let (curve, start, end, key) = if blended.contains(e) {
                    let (Some(c), Some(a), Some(b)) = (
                        contact_of(fi, *e),
                        corners.get(&(fi, vs)),
                        corners.get(&(fi, ve)),
                    ) else {
                        return Ok(None);
                    };
                    let (Some(ta), Some(tb)) = (a.on(*e), b.on(*e)) else {
                        return Ok(None);
                    };
                    (
                        c.clone(),
                        (a.point, ta),
                        (b.point, tb),
                        Key::Contact(fi, *e),
                    )
                } else {
                    let (pa, pb) = (sharp_end[&(*e, vs)], sharp_end[&(*e, ve)]);
                    let ta = param_on(&edge.curve, edge.range, pa, tol)?;
                    let tb = param_on(&edge.curve, edge.range, pb, tol)?;
                    (edge.curve.clone(), (pa, ta), (pb, tb), Key::Sharp(*e))
                };
                // Where the last piece stopped short of this one, a
                // bevel's end at the vertex runs between them.
                if let Some(p) = last
                    && p.distance(start.0) > tol.confusion()
                {
                    let Some(ends) = sections.get(&vs) else {
                        return Ok(None);
                    };
                    let joins = ends.iter().find(|(x, y, _)| {
                        (x.distance(p) <= tol.confusion() && y.distance(start.0) <= tol.confusion())
                            || (y.distance(p) <= tol.confusion()
                                && x.distance(start.0) <= tol.confusion())
                    });
                    let Some((_, _, bi)) = joins else {
                        return Ok(None);
                    };
                    let kind = &bevels[*bi].section;
                    if !flush(face, kind, p, start.0, tol)? {
                        return Ok(None);
                    }
                    let (k, along) = section(&mut points, &mut edges, &mut made, kind, p, start.0)?;
                    ring.push((k, along));
                }
                let along_edge = {
                    let t = param_on(&edge.curve, edge.range, start.0, tol)?;
                    edge.curve.d1_at(t, tol)? * if forward { 1.0 } else { -1.0 }
                };
                let rising = runs(&curve, start.1, along_edge, tol)?;
                let (range, against) = piece(&curve, start.1, end.1, rising, tol);
                let (i, j) = (points.index(start.0), points.index(end.0));
                let k = if let Some(k) = made.get(&key) {
                    *k
                } else {
                    let (lo_end, hi_end) = if against { (j, i) } else { (i, j) };
                    edges.push(SheetCurve {
                        ends: (lo_end, hi_end),
                        curve: curve.clone(),
                        range,
                    });
                    made.insert(key, edges.len() - 1);
                    edges.len() - 1
                };
                ring.push((k, !against));
                if first.is_none() {
                    first = Some(start.0);
                }
                last = Some(end.0);
            }
            // The loop closes: a bevel's end where it stopped short.
            if let (Some(p), Some(q)) = (last, first)
                && p.distance(q) > tol.confusion()
            {
                let v = loop_start_vertex(curved, face, walk);
                let Some(ends) = sections.get(&v) else {
                    return Ok(None);
                };
                let joins = ends.iter().find(|(x, y, _)| {
                    (x.distance(p) <= tol.confusion() && y.distance(q) <= tol.confusion())
                        || (y.distance(p) <= tol.confusion() && x.distance(q) <= tol.confusion())
                });
                let Some((_, _, bi)) = joins else {
                    return Ok(None);
                };
                let kind = &bevels[*bi].section;
                if !flush(face, kind, p, q, tol)? {
                    return Ok(None);
                }
                let (k, along) = section(&mut points, &mut edges, &mut made, kind, p, q)?;
                ring.push((k, along));
            }
            rings.push(ring);
        }
        rings_of.push(SheetRings {
            surface: face.surface.clone(),
            rings,
            reversed: face.reversed,
        });
    }
    // Each bevel: its contact on its first face, walked against the way
    // that face walks it (faces beside each other walk a shared edge
    // opposite ways, each as it faces out), its end, its contact on its
    // second face the other way round, and its other end.
    for b in bevels {
        let e = &curved.edges[&b.edge];
        let k0 = made[&Key::Contact(b.faces[0], b.edge)];
        let k1 = made[&Key::Contact(b.faces[1], b.edge)];
        let Some(walked) = face_walks(&rings_of[b.faces[0]], k0) else {
            return Ok(None);
        };
        // As the first face walks it facing out, then against that.
        let out_first = walked != rings_of[b.faces[0]].reversed;
        let s0 = !out_first;
        let s1 = !s0;
        let ends = |k: usize, forward: bool| {
            let (a, z) = edges[k].ends;
            if forward { (a, z) } else { (z, a) }
        };
        let (p0, q0) = ends(k0, s0);
        let (p1, q1) = ends(k1, s1);
        let mut ring: Vec<(usize, bool)> = vec![(k0, s0)];
        if q0 != p1 {
            let (a, z) = (points.points[q0], points.points[p1]);
            ring.push(section(
                &mut points,
                &mut edges,
                &mut made,
                &b.section,
                a,
                z,
            )?);
        }
        ring.push((k1, s1));
        if q1 != p0 {
            let (a, z) = (points.points[q1], points.points[p0]);
            ring.push(section(
                &mut points,
                &mut edges,
                &mut made,
                &b.section,
                a,
                z,
            )?);
        }
        // The bevel faces out: away from both faces at the edge.
        let mid =
            curved.points[&e.ends.0] + (curved.points[&e.ends.1] - curved.points[&e.ends.0]) * 0.5;
        let out = outward_at(curved, b.faces[0], mid) + outward_at(curved, b.faces[1], mid);
        let natural = natural_normal(&b.surface, points.points[p0], tol)?;
        let reversed = natural.dot(out) < 0.0;
        // The ring above runs as the bevel faces out; the sheet wants it
        // the surface's way.
        let ring = if reversed {
            ring.iter().rev().map(|(k, f)| (*k, !f)).collect()
        } else {
            ring
        };
        rings_of.push(SheetRings {
            surface: b.surface.clone(),
            rings: vec![ring],
            reversed,
        });
    }
    let sheet = build(model, &points.points, &edges, &rings_of, tol)?;
    let shell = ogeom_algo::make_shell(model, &sheet.faces)?.shape;
    let solid = ogeom_algo::make_solid(model, &[shell])?.shape;
    // An edge a face walks twice, once each way (a wall's or a cone
    // bevel's seam), needs both columns of the seam: the builder fits it
    // one pcurve, and a later walk of the face in its chart would not
    // close. Nothing else holds these nodes, so nothing is copied.
    ogeom_algo::join_seam_columns(model, &solid, tol)?;
    // P5: each bevel proved within the material its two faces bound there,
    // to the bands' fit target: under a plane, inside a wall (outside a
    // bore's).
    let fit = crate::marched::band_fit_target(tol);
    for (bi, b) in bevels.iter().enumerate() {
        let (Some(first), Some(second)) =
            (material(curved, b.faces[0]), material(curved, b.faces[1]))
        else {
            return Ok(None);
        };
        let region = ogeom_algo::HalfSpaces::Every(vec![first, second]);
        let face = &sheet.faces[curved.faces.len() + bi];
        if ogeom_algo::certify_face_within(model, face, &region, fit, tol)?.is_err() {
            return Ok(None);
        }
    }
    // History: faces modified, bevelled edges deleted into their bevels,
    // sharp edges modified into their pieces, vertices kept where nothing
    // moved off them.
    let mut history = History::new();
    for (fi, face) in curved.faces.iter().enumerate() {
        history.modify(&face.shape, sheet.faces[fi].clone());
    }
    for (bi, b) in bevels.iter().enumerate() {
        let edge = &curved.shapes[&b.edge];
        history.delete(edge);
        history.generate(edge, sheet.faces[curved.faces.len() + bi].clone());
    }
    for (e, _) in &curved.edges {
        if let Some(k) = made.get(&Key::Sharp(*e)) {
            history.modify(&curved.shapes[e], sheet.edges[*k].clone());
        }
    }
    for (v, p) in &curved.points {
        let kept = curved.faces.iter().enumerate().all(|(fi, _)| {
            corners
                .get(&(fi, *v))
                .is_none_or(|c| c.point.distance(*p) <= tol.confusion())
        });
        let shape = &curved.shapes[v];
        if kept
            && let Some(i) = points
                .points
                .iter()
                .position(|q| q.distance(*p) <= tol.confusion())
        {
            history.modify(shape, sheet.vertices[i].clone());
        } else {
            history.delete(shape);
        }
    }
    Ok(Some(Built::new(solid, history)))
}

/// The material face `face` bounds, as a region the certify reads: the
/// side of its plane under its outward normal, or the side of its cylinder
/// or cone the material lies on (inside a wall, outside a bore); `None` on
/// any other surface.
fn material(curved: &Curved, face: usize) -> Option<ogeom_algo::HalfSpaces> {
    let f = &curved.faces[face];
    let inside = !f.reversed;
    match &f.surface {
        SurfaceGeometry::Plane(plane) => {
            let plane = plane.plane();
            let normal = plane.normal().vector() * if f.reversed { -1.0 } else { 1.0 };
            Some(ogeom_algo::HalfSpaces::Below {
                normal,
                at: plane.frame().origin(),
            })
        }
        SurfaceGeometry::Cylinder(c) => {
            let cyl = c.cylinder();
            let axis = cyl.axis();
            Some(ogeom_algo::HalfSpaces::Revolved {
                at: axis.location,
                axis: axis.direction.vector(),
                radius: cyl.radius(),
                slope: 0.0,
                inside,
            })
        }
        SurfaceGeometry::Cone(c) => {
            let cone = c.cone();
            let frame = cone.frame();
            Some(ogeom_algo::HalfSpaces::Revolved {
                at: frame.origin(),
                axis: frame.z().vector(),
                radius: cone.reference_radius(),
                slope: cone.half_angle().tan(),
                inside,
            })
        }
        _ => None,
    }
}

/// Whether a strip's end from `p` to `q`, by its section, lies on `face`:
/// a plane holding both ends and, for an arc, its centre.
fn flush(
    face: &super::curved::CurvedFace,
    kind: &Section,
    p: Point,
    q: Point,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let SurfaceGeometry::Plane(plane) = &face.surface else {
        return Ok(false);
    };
    let plane = plane.plane();
    let on = |x: Point| plane.signed_distance_to(x).abs() <= tol.confusion();
    Ok(on(p)
        && on(q)
        && match kind {
            Section::Line => true,
            Section::Arc { spine } => {
                on(spine.point_at(param_on(spine, spine.domain(), p, tol)?, tol)?)
            }
        })
}

/// The vertex a face's loop starts at, walked the surface's way.
fn loop_start_vertex(curved: &Curved, face: &super::curved::CurvedFace, walk: &[(Id, bool)]) -> Id {
    let (e, w) = walk[0];
    let edge = &curved.edges[&e];
    if w != face.reversed {
        edge.ends.0
    } else {
        edge.ends.1
    }
}

/// Whether a face's rings walk sheet edge `k` the way it is stored.
fn face_walks(rings: &SheetRings, k: usize) -> Option<bool> {
    rings
        .rings
        .iter()
        .flatten()
        .find(|(e, _)| *e == k)
        .map(|(_, f)| *f)
}

/// A face's outward normal at `p`, nought where it is not a plane or a
/// cylinder.
fn outward_at(curved: &Curved, face: usize, p: Point) -> Vector {
    outward(curved, face, p).unwrap_or(Vector::ZERO)
}

/// A bevel surface's own normal at `p`: a plane's, or a cone's or a
/// cylinder's off its axis there.
fn natural_normal(surface: &SurfaceGeometry, p: Point, tol: Tolerances) -> OgeomResult<Vector> {
    let _ = tol;
    Ok(match surface {
        SurfaceGeometry::Plane(plane) => plane.plane().normal().vector(),
        SurfaceGeometry::Cylinder(c) => {
            let axis = c.cylinder().axis();
            let (o, a) = (axis.location, axis.direction.vector());
            (p - o) - a * (p - o).dot(a)
        }
        SurfaceGeometry::Cone(c) => {
            let cone = c.cone();
            let frame = cone.frame();
            let (o, a) = (frame.origin(), frame.z().vector());
            let radial = (p - o) - a * (p - o).dot(a);
            let r = radial / radial.magnitude();
            // Off the axis, leaning back against the cone's opening.
            r - a * cone.half_angle().tan()
        }
        _ => Vector::ZERO,
    })
}

/// The parameter of `curve`, over `range`, at the point `p` on it: by a
/// fine sampling then Newton on the distance.
pub(crate) fn param_on(
    curve: &Curve,
    range: (f64, f64),
    p: Point,
    tol: Tolerances,
) -> OgeomResult<f64> {
    const SAMPLES: usize = 64;
    let (lo, hi) = range;
    let mut best = (lo, f64::INFINITY);
    for i in 0..=SAMPLES {
        #[expect(clippy::cast_precision_loss, reason = "a small sample count")]
        let t = lo + (hi - lo) * (i as f64) / (SAMPLES as f64);
        let d = curve.point_at(t, tol)?.distance(p);
        if d < best.1 {
            best = (t, d);
        }
    }
    let mut t = best.0;
    for _ in 0..50 {
        let at = curve.point_at(t, tol)?;
        let d1 = curve.d1_at(t, tol)?;
        let step = (p - at).dot(d1) / d1.dot(d1);
        t = (t + step).clamp(lo, hi);
        if step.abs() <= tol.parametric() {
            break;
        }
    }
    Ok(t)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
    use super::*;
    use ogeom_math::{Direction, Frame};
    use ogeom_topo::{ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    fn up_at(x: f64, y: f64, z: f64) -> Frame {
        Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
    }

    /// A face of a small box at `(x, y, z)`, 0.2 on a side: its top.
    fn probe(m: &mut Model, x: f64, y: f64, z: f64) -> Shape {
        let b = ogeom_algo::make_box(m, up_at(x, y, z), (0.2, 0.2, 0.2), T)
            .unwrap()
            .shape;
        explore_unique(m, &b, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| {
                let (p, n) = ogeom_algo::face_normal(m, f, T).unwrap();
                n.z > 0.5 && (p.z - (z + 0.2)).abs() < 1e-9
            })
            .unwrap()
    }

    fn within(m: &Model, face: &Shape, region: &ogeom_algo::HalfSpaces) -> bool {
        ogeom_algo::certify_face_within(m, face, region, 1e-6, T)
            .unwrap()
            .is_ok()
    }

    /// The material a face bounds, as the certify reads it: on the D-block
    /// (a drum of radius 3 halved along its axis) under the flat (`x >= 0`)
    /// and inside the wall; on a tube bored to 2, outside the bore. A face
    /// inside is proved within, one past the flat, outside the wall or in
    /// the bore is found past.
    #[test]
    fn a_face_s_material_is_its_side_of_its_plane_or_wall() {
        let mut m = Model::new();
        let full = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 10.0, T)
            .unwrap()
            .shape;
        let half = ogeom_algo::make_box(&mut m, up_at(-4.0, -4.0, -1.0), (4.0, 8.0, 12.0), T)
            .unwrap()
            .shape;
        let d = ogeom_bool::cut(&mut m, &full, &half, T).unwrap().shape;
        let curved = read(&m, &d, T).unwrap().unwrap();
        let of = |kind: fn(&SurfaceGeometry) -> bool| {
            let fi = (0..curved.faces.len())
                .find(|fi| {
                    kind(&curved.faces[*fi].surface) && {
                        // The flat at x = 0, not a cap.
                        match &curved.faces[*fi].surface {
                            SurfaceGeometry::Plane(p) => p.plane().normal().vector().x.abs() > 0.5,
                            _ => true,
                        }
                    }
                })
                .unwrap();
            material(&curved, fi).unwrap()
        };
        let flat = of(|s| matches!(s, SurfaceGeometry::Plane(_)));
        let wall = of(|s| matches!(s, SurfaceGeometry::Cylinder(_)));
        let inside = probe(&mut m, 1.0, 0.0, 5.0);
        let past_flat = probe(&mut m, -1.0, 0.0, 5.0);
        let past_wall = probe(&mut m, 3.5, 0.0, 5.0);
        assert!(within(&m, &inside, &flat) && within(&m, &inside, &wall));
        assert!(!within(&m, &past_flat, &flat));
        assert!(!within(&m, &past_wall, &wall));

        let outer = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 4.0, T)
            .unwrap()
            .shape;
        let inner = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, -1.0), 2.0, 6.0, T)
            .unwrap()
            .shape;
        let tube = ogeom_bool::cut(&mut m, &outer, &inner, T).unwrap().shape;
        let curved = read(&m, &tube, T).unwrap().unwrap();
        let bore = (0..curved.faces.len())
            .find(|fi| match &curved.faces[*fi].surface {
                SurfaceGeometry::Cylinder(c) => (c.cylinder().radius() - 2.0).abs() < 1e-9,
                _ => false,
            })
            .unwrap();
        let bore = material(&curved, bore).unwrap();
        let in_wall = probe(&mut m, 2.3, 0.0, 2.0);
        let in_bore = probe(&mut m, 0.5, 0.0, 2.0);
        assert!(within(&m, &in_wall, &bore));
        assert!(!within(&m, &in_bore, &bore));
    }

    /// The D-shaped block's two straight edges rounded at 0.5 by strips a
    /// section law would hand over, built here in closed form: the ball
    /// between the flat `x = 0` and the wall of radius 3 has its centre at
    /// `(r, ±sqrt((3 - r)^2 - r^2))`, touching the flat level with it and
    /// the wall on its radius through it; the strip is the cylinder of the
    /// ball's radius about that line, its ends arcs about it. The assembly
    /// builds them through the caps' flush ends to a valid solid, each strip
    /// proved within its faces' material, of the half disc less the
    /// line-and-circle corner twice along the 10: the corner the kite from
    /// the corner to the contacts and the centre, less the ball's sector,
    /// plus the wall's segment beyond the chord.
    #[test]
    fn round_strips_assemble_through_flush_ends() {
        use core::f64::consts::PI;
        use ogeom_geom::CylinderSurface;
        use ogeom_math::Cylinder;
        let mut m = Model::new();
        let full = ogeom_algo::make_cylinder(&mut m, up_at(0.0, 0.0, 0.0), 3.0, 10.0, T)
            .unwrap()
            .shape;
        let half = ogeom_algo::make_box(&mut m, up_at(-4.0, -4.0, -1.0), (4.0, 8.0, 12.0), T)
            .unwrap()
            .shape;
        let d = ogeom_bool::cut(&mut m, &full, &half, T).unwrap().shape;
        let curved = read(&m, &d, T).unwrap().unwrap();
        let (big, r) = (3.0_f64, 0.5_f64);
        let reach = ((big - r) * (big - r) - r * r).sqrt();
        let mut strips: Vec<Strip> = Vec::new();
        for side in [1.0_f64, -1.0] {
            let edge = *curved
                .edges
                .iter()
                .find(|(_, e)| {
                    let (a, b) = (curved.points[&e.ends.0], curved.points[&e.ends.1]);
                    a.x.abs() < 1e-9
                        && b.x.abs() < 1e-9
                        && (a.y - big * side).abs() < 1e-9
                        && (a.z - b.z).abs() > 1.0
                })
                .unwrap()
                .0;
            let faces = faces_of(&curved, edge);
            let (cx, cy) = (r, side * reach);
            let on_wall = (cx * big / (big - r), cy * big / (big - r));
            let line = |x: f64, y: f64| -> Curve {
                LineCurve::segment(Point::new(x, y, -5.0), Point::new(x, y, 15.0), T)
                    .unwrap()
                    .into()
            };
            let contact = |fi: usize| match &curved.faces[fi].surface {
                SurfaceGeometry::Plane(_) => line(0.0, cy),
                _ => line(on_wall.0, on_wall.1),
            };
            let surface: SurfaceGeometry = CylinderSurface::new(
                Cylinder::new(up_at(cx, cy, -5.0), r, T).unwrap(),
                (0.0, 20.0),
            )
            .unwrap()
            .into();
            strips.push(Strip {
                edge,
                faces: [faces[0], faces[1]],
                contacts: [contact(faces[0]), contact(faces[1])],
                surface,
                section: Section::Arc {
                    spine: line(cx, cy),
                },
            });
        }
        let mut contacts: FastMap<(usize, Id), Curve> = FastMap::default();
        for st in &strips {
            for k in 0..2 {
                contacts.insert((st.faces[k], st.edge), st.contacts[k].clone());
            }
        }
        let corners = curved_corners(&curved, &contacts, T).unwrap().unwrap();
        let blended: FastSet<Id> = strips.iter().map(|st| st.edge).collect();
        let built = assemble(&mut m, &curved, &strips, &blended, &corners, T)
            .unwrap()
            .expect("the round strips assemble");
        assert!(ogeom_algo::check(&m, &built.shape, T).unwrap().is_valid());
        // The corner: kite less sector plus the wall's segment.
        let (p0, a, c, b) = (
            (0.0, big),
            (0.0, reach),
            (r, reach),
            (r * big / (big - r), reach * big / (big - r)),
        );
        let quad = 0.5
            * ((p0.0 * a.1 - a.0 * p0.1)
                + (a.0 * c.1 - c.0 * a.1)
                + (c.0 * b.1 - b.0 * c.1)
                + (b.0 * p0.1 - p0.0 * b.1))
                .abs();
        let sector = {
            let (u, w) = ((a.0 - c.0, a.1 - c.1), (b.0 - c.0, b.1 - c.1));
            let angle = ((u.0 * w.0 + u.1 * w.1) / (r * r)).clamp(-1.0, 1.0).acos();
            angle * r * r / 2.0
        };
        let segment = {
            let phi = ((b.0 * p0.0 + b.1 * p0.1) / (big * big))
                .clamp(-1.0, 1.0)
                .acos();
            big * big / 2.0 * (phi - phi.sin())
        };
        let corner = quad - sector + segment;
        let want = (PI * big * big / 2.0 - 2.0 * corner) * 10.0;
        let got = ogeom_algo::volume_properties(
            &m,
            &built.shape,
            ogeom_mesh::Deflection {
                chord: 1e-5,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap()
        .mass;
        assert!((got - want).abs() < 1e-6, "volume {got} against {want}");
    }
}
