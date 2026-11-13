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
use crate::chamfer_sphere::{SphereSeat, sphere_bevel_rings};
use crate::support::RevolvedSeat;

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
/// seat, a vertex or the result lies outside the classes built here, or
/// anything in building it fails.
///
/// # Errors
///
/// [`OgeomError::Cancelled`](ogeom_core::OgeomError::Cancelled) only.
pub(crate) fn chamfer(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    declining(chamfer_with(model, solid, specs, false, tol))
}

/// What the curved path builds, any failure in it but a cancellation
/// handed back for the caller's own path.
fn declining(found: OgeomResult<Option<Built>>) -> OgeomResult<Option<Built>> {
    match found {
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(None),
        ok => ok,
    }
}

/// [`chamfer`] on a solid bounded by planes, which the planar engine
/// handed back: `None` unless a vertex of three bevels or more has its
/// faces' corners off one plane, where the facet is a patch.
///
/// # Errors
///
/// As [`chamfer`].
pub(crate) fn chamfer_patched(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    declining(chamfer_with(model, solid, specs, true, tol))
}

fn chamfer_with(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    patched_only: bool,
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
    // A wall the bevels of its two rims consume: its strips meet along
    // their ridge instead.
    let Some(bands) = consumed_bands(&curved, &mut bevels, tol)? else {
        return Ok(None);
    };
    let mut contacts: FastMap<(usize, Id), Curve> = FastMap::default();
    for b in &bevels {
        for k in 0..2 {
            if !bands.faces.contains(&b.faces[k]) {
                contacts.insert((b.faces[k], b.edge), b.contacts[k].clone());
            }
        }
    }
    // P2, a consumed wall's corners where the ridge stands at each vertex.
    let Some(mut corners) = curved_corners(&curved, &contacts, tol)? else {
        return Ok(None);
    };
    for ((fi, v), point) in &bands.corners {
        let on: Vec<(Id, f64)> = bevels
            .iter()
            .filter(|b| b.faces.contains(fi))
            .filter_map(|b| {
                let k = b.faces.iter().position(|f| f == fi)?;
                let c = &b.contacts[k];
                Some((b.edge, param_on(c, c.domain(), *point, tol).ok()?))
            })
            .collect();
        corners.insert((*fi, *v), Corner { point: *point, on });
    }
    let blended: FastSet<Id> = bevels.iter().map(|b| b.edge).collect();
    if patched_only && !off_plane_corner(&curved, &blended, &corners, tol) {
        return Ok(None);
    }
    // P4 and P5. Whatever fails here, its assembly or the solid's check,
    // is outside the cases built, and the caller keeps its own path.
    let built = match assemble_consuming(
        model,
        &curved,
        &bevels,
        &blended,
        &corners,
        &bands.faces,
        tol,
    ) {
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

/// The bilinear patch through the four corners `c`, wound about the side
/// it faces: its sides the segments between them in turn, its own normal
/// facing that side.
fn bilinear(c: &[Point], tol: Tolerances) -> OgeomResult<SurfaceGeometry> {
    let knots = ogeom_math::KnotVector::clamped_uniform(1, 2)?;
    let grid = ogeom_math::ControlGrid::new(vec![c[0], c[3], c[1], c[2]], 2, 2)?;
    Ok(ogeom_geom::BSplineSurface::new(knots.clone(), knots, &grid, tol)?.into())
}

/// Whether a vertex of three bevels or more has its faces' corners off
/// one plane.
fn off_plane_corner(
    curved: &Curved,
    blended: &FastSet<Id>,
    corners: &FastMap<(usize, Id), Corner>,
    tol: Tolerances,
) -> bool {
    curved.points.keys().any(|v| {
        let at_v = curved
            .edges
            .iter()
            .filter(|(e, c)| blended.contains(*e) && (c.ends.0 == *v || c.ends.1 == *v))
            .count();
        if at_v < 3 {
            return false;
        }
        let ring: Vec<Point> = (0..curved.faces.len())
            .filter_map(|fi| corners.get(&(fi, *v)).map(|c| c.point))
            .collect();
        if ring.len() < 4 {
            return false;
        }
        // The plane through the corner and the two farthest from it and
        // from each other; any corner off it.
        let a = ring[0];
        let Some(b) = ring
            .iter()
            .copied()
            .max_by(|x, y| x.distance(a).total_cmp(&y.distance(a)))
        else {
            return false;
        };
        let Some(c) = ring.iter().copied().max_by(|x, y| {
            (*x - a)
                .cross(b - a)
                .magnitude()
                .total_cmp(&(*y - a).cross(b - a).magnitude())
        }) else {
            return false;
        };
        let n = (b - a).cross(c - a);
        let m = n.magnitude();
        m > tol.confusion() * tol.confusion()
            && ring
                .iter()
                .any(|p| (*p - a).dot(n / m).abs() > tol.confusion())
    })
}

/// The walls a request's bevels consume, their ridges and corners.
#[derive(Default)]
struct Bands {
    /// Each consumed wall, by index.
    faces: FastSet<usize>,
    /// Each consumed wall's corner at each of its vertices: where the
    /// ridge stands at the vertex's station.
    corners: FastMap<(usize, Id), Point>,
}

/// The walls whose two rims' bevels cross on them (a disc too thin for
/// both rims' chamfers): a cylinder or a cone whose bevelled edges lie on
/// two coaxial circles, the contact of the higher at or below the lower's.
/// Each such wall's strips take as their contact on it the circle their
/// surfaces meet along between the rims, the ridge. `None` where a wall's
/// contacts cross and no such ridge is found.
fn consumed_bands(
    curved: &Curved,
    strips: &mut [Strip],
    tol: Tolerances,
) -> OgeomResult<Option<Bands>> {
    let mut bands = Bands::default();
    for (fi, face) in curved.faces.iter().enumerate() {
        let (o, a) = match &face.surface {
            SurfaceGeometry::Cylinder(c) => {
                let axis = c.cylinder().axis();
                (axis.location, axis.direction.vector())
            }
            SurfaceGeometry::Cone(c) => {
                let f = c.cone().frame();
                (f.origin(), f.z().vector())
            }
            _ => continue,
        };
        let height = |p: Point| (p - o).dot(a);
        // The strips on this wall, grouped by their rim's height: the
        // rim's height and the contact's.
        let mut groups: Vec<(f64, f64, Vec<usize>)> = Vec::new();
        for (si, st) in strips.iter().enumerate() {
            let Some(k) = st.faces.iter().position(|f| *f == fi) else {
                continue;
            };
            let (Curve::Circle(rim), Curve::Circle(contact)) =
                (&curved.edges[&st.edge].curve, &st.contacts[k])
            else {
                continue;
            };
            let (h_rim, h_contact) = (
                height(rim.circle().centre()),
                height(contact.circle().centre()),
            );
            match groups
                .iter_mut()
                .find(|g| (g.0 - h_rim).abs() <= tol.confusion())
            {
                Some(g) => g.2.push(si),
                None => groups.push((h_rim, h_contact, vec![si])),
            }
        }
        let [high, low] = match groups.as_slice() {
            [x, y] if x.0 > y.0 => [x.clone(), y.clone()],
            [x, y] => [y.clone(), x.clone()],
            _ => continue,
        };
        if high.1 > low.1 + tol.confusion() {
            continue;
        }
        // The ridge: where the high rim's bevel meets the low rim's,
        // between the two rims. Both turn about the wall's axis, each
        // radius running linearly with the height, so they meet in the
        // circle where the two radii agree.
        let (Some(ra), Some(rb)) = (
            radius_law(&strips[high.2[0]].surface, o, a, tol),
            radius_law(&strips[low.2[0]].surface, o, a, tol),
        ) else {
            return Ok(None);
        };
        if (ra.1 - rb.1).abs() <= tol.angular() {
            return Ok(None);
        }
        let h = (rb.0 - ra.0) / (ra.1 - rb.1);
        let r = ra.1.mul_add(h, ra.0);
        if !(h > low.0 + tol.confusion() && h < high.0 - tol.confusion() && r > tol.confusion()) {
            return Ok(None);
        }
        let x = if a.x.abs() < 0.9 {
            Vector::X
        } else {
            Vector::Y
        };
        let x = x - a * x.dot(a);
        let frame = Frame::new(
            o + a * h,
            Direction::new(a, tol)?,
            Direction::new(x, tol)?,
            tol,
        )?;
        let ridge: Curve = CircleCurve::new(Circle::new(frame, r, tol)?).into();
        for si in high.2.iter().chain(&low.2) {
            let st = &mut strips[*si];
            if let Some(k) = st.faces.iter().position(|f| *f == fi) {
                st.contacts[k] = ridge.clone();
            }
        }
        for (e, _) in face.loops.iter().flatten() {
            let c = &curved.edges[e];
            for v in [c.ends.0, c.ends.1] {
                let at = curved.points[&v];
                let t = param_on(&ridge, ridge.domain(), at, tol)?;
                bands.corners.insert((fi, v), ridge.point_at(t, tol)?);
            }
        }
        bands.faces.insert(fi);
    }
    Ok(Some(bands))
}

/// A cone's or cylinder's radius as a line in the height along the axis
/// through `o` along the unit `a`: its radius at height nought and how it
/// grows. `None` off a cone or a cylinder, or about another axis.
fn radius_law(
    surface: &SurfaceGeometry,
    o: Point,
    a: Vector,
    tol: Tolerances,
) -> Option<(f64, f64)> {
    let (origin, axis, r0, slope) = match surface {
        SurfaceGeometry::Cylinder(c) => {
            let cyl = c.cylinder();
            let axis = cyl.axis();
            (axis.location, axis.direction.vector(), cyl.radius(), 0.0)
        }
        SurfaceGeometry::Cone(c) => {
            let cone = c.cone();
            let f = cone.frame();
            (
                f.origin(),
                f.z().vector(),
                cone.reference_radius(),
                cone.half_angle().tan(),
            )
        }
        _ => return None,
    };
    let off = origin - o;
    if axis.cross(a).magnitude() > tol.angular()
        || (off - a * off.dot(a)).magnitude() > tol.confusion()
    {
        return None;
    }
    let sign = axis.dot(a).signum();
    let at = off.dot(a);
    // r(h) = r0 + slope * sign * (h - at).
    let grows = slope * sign;
    Some((grows.mul_add(-at, r0), grows))
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

/// A face's outward unit normal at `p`, on a plane, a cylinder, a sphere or
/// a torus; `None` on any other surface.
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
        SurfaceGeometry::Sphere(sphere) => {
            let r = p - sphere.sphere().centre();
            let m = r.magnitude();
            (m > 0.0).then(|| r / m * sign)
        }
        SurfaceGeometry::Torus(torus) => {
            let torus = torus.torus();
            let frame = torus.frame();
            let (o, a) = (frame.origin(), frame.z().vector());
            let q = p - o;
            let radial = q - a * q.dot(a);
            let m = radial.magnitude();
            if m <= 0.0 {
                return None;
            }
            let r = p - (o + radial / m * torus.major_radius());
            let n = r.magnitude();
            (n > 0.0).then(|| r / n * sign)
        }
        _ => None,
    }
}

/// A point of face `face` away from edge `edge`: a vertex of its outline
/// not on the edge, or where it has none, the middle of another of its
/// edges. `None` where it has neither.
fn away(curved: &Curved, face: usize, edge: Id) -> Option<Point> {
    let ends = curved.edges[&edge].ends;
    let walks = || curved.faces[face].loops.iter().flatten();
    walks()
        .flat_map(|(e, _)| {
            let c = &curved.edges[e];
            [c.ends.0, c.ends.1]
        })
        .find(|v| *v != ends.0 && *v != ends.1)
        .map(|v| curved.points[&v])
        .or_else(|| {
            // A face of two edges between the same two vertices: the
            // middle of its other edge.
            walks().filter(|(e, _)| *e != edge).find_map(|(e, _)| {
                let c = &curved.edges[e];
                c.curve
                    .point_at(f64::midpoint(c.range.0, c.range.1), Tolerances::default())
                    .ok()
            })
        })
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
        Curve::Circle(circle) => match legs(model, solid, shape, spec, tol) {
            Ok(Legs::Rim(seat, on_wall, on_cap)) => {
                rim(curved, edge, faces, &seat, on_wall, on_cap, tol)
            }
            Ok(Legs::Sphere(seat, on_wall, on_cap)) => {
                sphere_rim(curved, edge, faces, &seat, on_wall, on_cap, tol)
            }
            _ => torus_rim(curved, edge, circle, faces, spec, tol),
        },
        _ => extruded_crease(curved, edge, e, faces, spec, tol),
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
fn rim(
    curved: &Curved,
    edge: Id,
    faces: [usize; 2],
    seat: &RevolvedSeat,
    on_wall: f64,
    on_cap: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
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
    let depth = wall_depth(seat, on_wall);
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
    let surface = rim_bevel_surface(seat, depth, cap_rho, wall_rho, tol)?;
    Ok(Some(Strip {
        edge,
        faces,
        contacts,
        surface,
        section: Section::Line,
    }))
}

/// A rim between a planar cap and a sphere centred on its axis, as the
/// wedge path reads it ([`legs`]): the ring along the cap, the ring down
/// the sphere's great circle, and the cone through both. `None` where the
/// legs leave no such bevel.
fn sphere_rim(
    curved: &Curved,
    edge: Id,
    faces: [usize; 2],
    seat: &SphereSeat,
    on_wall: f64,
    on_cap: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
    let on = |face: &Shape| {
        faces
            .iter()
            .position(|f| curved.faces[*f].shape.is_same(face))
    };
    let (Some(cap), Some(wall)) = (on(&seat.cap_face), on(&seat.wall_face)) else {
        return Ok(None);
    };
    if cap == wall || on_wall >= seat.to_pole() {
        return Ok(None);
    }
    let Ok(bevel) = sphere_bevel_rings(seat, on_wall, on_cap, tol) else {
        return Ok(None);
    };
    let ring = |centre: Point, radius: f64| -> OgeomResult<Curve> {
        let frame = Frame::new(centre, Direction::new(seat.up, tol)?, seat.x_ref, tol)?;
        Ok(CircleCurve::new(Circle::new(frame, radius, tol)?).into())
    };
    let mut contacts = [
        ring(seat.centre, bevel.cap_rho)?,
        ring(seat.centre, bevel.cap_rho)?,
    ];
    contacts[wall] = ring(seat.centre + seat.up * bevel.wall_h, bevel.wall_rho)?;
    Ok(Some(Strip {
        edge,
        faces,
        contacts,
        surface: bevel.surface,
        section: Section::Line,
    }))
}

/// A rim between a planar cap square to a torus's axis and the torus, the
/// material inside its tube or outside it (a groove): in the meridian
/// half-plane through a point of the rim, `(rho, z)` from the torus's
/// centre along its axis, the leg on the cap runs straight along it and the
/// leg on the wall round the tube's circle, and the bevel is the cone
/// through the two contact rings. `None` where the faces are not such a
/// plane and torus, the rim is not coaxial, the edge is not convex, or the
/// legs leave no bevel in the material.
fn torus_rim(
    curved: &Curved,
    edge: Id,
    rim: &CircleCurve,
    faces: [usize; 2],
    spec: &Chamfer,
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
    let (eps, ang) = (tol.confusion(), tol.angular());
    let surface = |k: usize| &curved.faces[faces[k]].surface;
    let (cap, wall) = match (surface(0), surface(1)) {
        (SurfaceGeometry::Plane(_), SurfaceGeometry::Torus(_)) => (0, 1),
        (SurfaceGeometry::Torus(_), SurfaceGeometry::Plane(_)) => (1, 0),
        _ => return Ok(None),
    };
    let (SurfaceGeometry::Plane(plane), SurfaceGeometry::Torus(torus)) =
        (surface(cap), surface(wall))
    else {
        return Ok(None);
    };
    // The material inside the tube, or outside it (a groove).
    let within: f64 = if curved.faces[faces[wall]].reversed {
        -1.0
    } else {
        1.0
    };
    let torus = torus.torus();
    let frame = torus.frame();
    let (o, a) = (frame.origin(), frame.z().vector());
    let (big, small) = (torus.major_radius(), torus.minor_radius());
    let circle = rim.circle();
    let centre = circle.centre();
    let off = centre - o;
    let up = plane.plane().normal().vector()
        * if curved.faces[faces[cap]].reversed {
            -1.0
        } else {
            1.0
        };
    if up.cross(a).magnitude() > ang
        || circle.frame().z().vector().cross(a).magnitude() > ang
        || (off - a * off.dot(a)).magnitude() > eps
    {
        return Ok(None);
    }
    let rim_rho = circle.radius();
    let rim_z = off.dot(a);
    if ((rim_rho - big).hypot(rim_z) - small).abs() > eps {
        return Ok(None);
    }
    // The meridian: the cap's outward normal, the wall's at the rim, and
    // which way along the cap from the rim the cap lies (toward the axis
    // where it has no other edge).
    let cap_up = up.dot(a).signum();
    let wall_out = (within * (rim_rho - big) / small, within * rim_z / small);
    let edge_ids: Vec<Id> = curved.faces[faces[cap]]
        .loops
        .iter()
        .flatten()
        .map(|(e, _)| *e)
        .filter(|e| *e != edge)
        .collect();
    let mut side = 0.0_f64;
    for e in &edge_ids {
        let c = &curved.edges[e];
        let p = c.curve.point_at(f64::midpoint(c.range.0, c.range.1), tol)?;
        let q = p - o;
        let gap = (q - a * q.dot(a)).magnitude() - rim_rho;
        if gap.abs() <= eps {
            // The rim's other arcs.
            continue;
        }
        let s = gap.signum();
        if side == 0.0 {
            side = s;
        } else if side != s {
            return Ok(None);
        }
    }
    if side == 0.0 {
        side = -1.0;
    }
    // Convex: the cap leaves the rim behind the wall's outward side.
    if side * wall_out.0 >= -ang {
        return Ok(None);
    }
    // Round the tube from the rim, away from the cap's outward side.
    let alpha = rim_z.atan2(rim_rho - big);
    let turn: f64 = if alpha.cos() * cap_up < 0.0 {
        1.0
    } else {
        -1.0
    };
    let on_tube = |t: f64| {
        let b = turn.mul_add(t / small, alpha);
        (small.mul_add(b.cos(), big), small * b.sin())
    };
    let named = |face: &Shape| {
        if curved.faces[faces[cap]].shape.is_same(face) {
            Some(true)
        } else if curved.faces[faces[wall]].shape.is_same(face) {
            Some(false)
        } else {
            None
        }
    };
    // The legs: along the cap, and round the tube.
    let (on_cap, on_wall) = match spec {
        Chamfer::Symmetric(d) => (*d, *d),
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => match named(face) {
            Some(true) => (*on_face, *on_other),
            Some(false) => (*on_other, *on_face),
            None => return Ok(None),
        },
        Chamfer::Angle {
            face,
            distance,
            angle,
        } => {
            let (d, phi) = (*distance, *angle);
            match named(face) {
                Some(true) => {
                    // From the cap's contact back toward the rim and into
                    // the material, to where it meets the tube.
                    let c = (side.mul_add(d, rim_rho), rim_z);
                    let dir = (-side * phi.cos(), -cap_up * phi.sin());
                    let rel = (c.0 - big, c.1);
                    let b = rel.0.mul_add(dir.0, rel.1 * dir.1);
                    let k = rel.0.mul_add(rel.0, rel.1 * rel.1) - small * small;
                    let disc = b.mul_add(b, -k);
                    // Inside the tube the contact stands in it and the line
                    // leaves through the far root; outside, it first meets
                    // the tube at the near one.
                    if within * k >= 0.0 || disc < 0.0 {
                        return Ok(None);
                    }
                    let t = -b + within * disc.sqrt();
                    if t <= eps {
                        return Ok(None);
                    }
                    let w = (t.mul_add(dir.0, c.0), t.mul_add(dir.1, c.1));
                    let beta = w.1.atan2(w.0 - big);
                    let mut swept = turn * (beta - alpha);
                    if swept < 0.0 {
                        swept += core::f64::consts::TAU;
                    }
                    (d, swept * small)
                }
                Some(false) => {
                    // From the wall's contact back along the tube toward
                    // the rim and into the material, to the cap's line.
                    let w = on_tube(d);
                    let b = turn.mul_add(d / small, alpha);
                    let back = (turn * b.sin(), -turn * b.cos());
                    let into = (-within * b.cos(), -within * b.sin());
                    let dir = (
                        back.0.mul_add(phi.cos(), into.0 * phi.sin()),
                        back.1.mul_add(phi.cos(), into.1 * phi.sin()),
                    );
                    if dir.1.abs() <= ang {
                        return Ok(None);
                    }
                    let t = (rim_z - w.1) / dir.1;
                    if t <= eps {
                        return Ok(None);
                    }
                    let c_rho = t.mul_add(dir.0, w.0);
                    ((c_rho - rim_rho) * side, d)
                }
                None => return Ok(None),
            }
        }
    };
    if !(on_cap.is_finite() && on_cap > eps && on_wall.is_finite() && on_wall > eps)
        || on_wall >= core::f64::consts::PI * small
    {
        return Ok(None);
    }
    let c = (side.mul_add(on_cap, rim_rho), rim_z);
    let w = on_tube(on_wall);
    if c.0 <= eps || w.0 <= eps || (w.1 - c.1).abs() <= eps {
        return Ok(None);
    }
    // The chord from the wall's contact to the cap's meets the tube's
    // circle again `-2 (W - T).(C - W) / |C - W|^2` along it, which must
    // not fall between the two.
    let chord = (c.0 - w.0, c.1 - w.1);
    let spoke = (w.0 - big, w.1);
    let length2 = chord.0.mul_add(chord.0, chord.1 * chord.1);
    let again = -2.0 * spoke.0.mul_add(chord.0, spoke.1 * chord.1) / length2;
    if again * length2.sqrt() > eps && again < 1.0 {
        return Ok(None);
    }
    let x = circle.frame().x().vector();
    let x = Direction::new(x - a * x.dot(a), tol)?;
    let axis_dir = Direction::new(a, tol)?;
    let ring = |z: f64, radius: f64| -> OgeomResult<Curve> {
        let frame = Frame::new(o + a * z, axis_dir, x, tol)?;
        Ok(CircleCurve::new(Circle::new(frame, radius, tol)?).into())
    };
    let mut contacts = [ring(c.1, c.0)?, ring(c.1, c.0)?];
    contacts[wall] = ring(w.1, w.0)?;
    let surface = coaxial_cone(
        Frame::new(o + a * c.1, axis_dir, x, tol)?,
        c.0,
        (w.0, w.1 - c.1),
        tol,
    )?;
    Ok(Some(Strip {
        edge,
        faces,
        contacts,
        surface,
        section: Section::Line,
    }))
}

/// How many stations the cap's contact along an extruded crease is fitted
/// through.
const CREASE_STATIONS: u32 = 128;

/// How far an extruded crease's ruled bevel runs on past each contact, as
/// a share of its rulings' length.
const BEVEL_RUN_ON: f64 = 0.5;

/// A crease between a plane and the wall extruded from it square to it (a
/// prism's end cap and the side its profile's curve sweeps), the crease
/// any curve of the profile: in the section square to the crease the
/// corner is square, the leg on the wall runs straight down it and the leg
/// on the cap straight in along it. The wall's contact is the crease moved
/// down the wall, exactly on it; the cap's is the crease's offset in the
/// cap's plane, fitted through stations at the crease's own parameters, so
/// it lies in the plane; the bevel is the ruled surface between the two on
/// their shared knots, each ruling in its station's section, on the wall's
/// material side. `None` where the faces are not such a plane and wall, the
/// edge is not convex, or the offset's fit misses the bands' fit target.
fn extruded_crease(
    curved: &Curved,
    edge: Id,
    e: &CurvedEdge,
    faces: [usize; 2],
    spec: &Chamfer,
    tol: Tolerances,
) -> OgeomResult<Option<Strip>> {
    let (eps, ang) = (tol.confusion(), tol.angular());
    let surface = |k: usize| &curved.faces[faces[k]].surface;
    let (cap, wall) = match (surface(0), surface(1)) {
        (SurfaceGeometry::Plane(_), SurfaceGeometry::Extrusion(_)) => (0, 1),
        (SurfaceGeometry::Extrusion(_), SurfaceGeometry::Plane(_)) => (1, 0),
        _ => return Ok(None),
    };
    let (SurfaceGeometry::Plane(plane), SurfaceGeometry::Extrusion(extrusion)) =
        (surface(cap), surface(wall))
    else {
        return Ok(None);
    };
    let plane = plane.plane();
    let up = plane.normal().vector()
        * if curved.faces[faces[cap]].reversed {
            -1.0
        } else {
            1.0
        };
    if extrusion.direction().vector().cross(up).magnitude() > ang {
        return Ok(None);
    }
    let Curve::BSpline(crease) = &e.curve else {
        return Ok(None);
    };
    if crease.is_rational() {
        return Ok(None);
    }
    let (t0, t1) = e.range;
    let at = |k: u32| (t1 - t0).mul_add(f64::from(k) / f64::from(CREASE_STATIONS), t0);
    // The crease in the cap's plane, and which way into the cap and down
    // the wall: toward the faces' other points.
    let mid = e.curve.point_at(at(CREASE_STATIONS / 2), tol)?;
    for k in 0..=CREASE_STATIONS {
        if plane
            .signed_distance_to(e.curve.point_at(at(k), tol)?)
            .abs()
            > eps
        {
            return Ok(None);
        }
    }
    let (Some(in_cap), Some(in_wall)) = (
        away(curved, faces[cap], edge),
        away(curved, faces[wall], edge),
    ) else {
        return Ok(None);
    };
    // Convex: the wall runs down from the cap, under its outward side.
    if (in_wall - mid).dot(up) >= -eps {
        return Ok(None);
    }
    let across = |t: f64| -> OgeomResult<Vector> {
        let tangent = e.curve.d1_at(t, tol)?;
        let m = up.cross(tangent);
        Ok(m / m.magnitude())
    };
    let sign = if (in_cap - mid).dot(across(at(CREASE_STATIONS / 2))?) > 0.0 {
        1.0
    } else {
        -1.0
    };
    let named = |face: &Shape| {
        if curved.faces[faces[cap]].shape.is_same(face) {
            Some(true)
        } else if curved.faces[faces[wall]].shape.is_same(face) {
            Some(false)
        } else {
            None
        }
    };
    let (on_cap, on_wall) = match spec {
        Chamfer::Symmetric(d) => (*d, *d),
        Chamfer::Distances {
            face,
            on_face,
            on_other,
        } => match named(face) {
            Some(true) => (*on_face, *on_other),
            Some(false) => (*on_other, *on_face),
            None => return Ok(None),
        },
        // The corner is square: the other leg is the named one's along the
        // angle's tangent.
        Chamfer::Angle {
            face,
            distance,
            angle,
        } => match named(face) {
            Some(true) => (*distance, distance * angle.tan()),
            Some(false) => (distance * angle.tan(), *distance),
            None => return Ok(None),
        },
    };
    if !(on_cap.is_finite() && on_cap > eps && on_wall.is_finite() && on_wall > eps) {
        return Ok(None);
    }
    // The wall's contact: the crease moved down the wall.
    let mut on_side = ogeom_geom::BSplineCurve::new(
        crease.knots().clone(),
        crease
            .control_points()
            .iter()
            .map(|w| w.point() + up * -on_wall)
            .collect(),
        tol,
    )?;
    // The cap's: the offset through stations at the crease's parameters.
    let mut parameters = Vec::with_capacity(CREASE_STATIONS as usize + 1);
    let mut points = Vec::with_capacity(CREASE_STATIONS as usize + 1);
    for k in 0..=CREASE_STATIONS {
        let t = at(k);
        parameters.push(t);
        points.push(e.curve.point_at(t, tol)? + across(t)? * (sign * on_cap));
    }
    let fit = crate::marched::band_fit_target(tol);
    let fitted = ogeom_geom::fit::fit_points_at(&parameters, &points, 3, eps, tol)?;
    if !fitted.met && fitted.error > fit {
        return Ok(None);
    }
    let mut on_top = fitted.curve;
    // The fitted offset between its stations, measured at their middles.
    for k in 0..CREASE_STATIONS {
        let t = f64::midpoint(at(k), at(k + 1));
        let want = e.curve.point_at(t, tol)? + across(t)? * (sign * on_cap);
        if on_top.point_at(t, tol)?.distance(want) > fit {
            return Ok(None);
        }
    }
    if (on_top.domain().0 - on_side.domain().0).abs() > tol.parametric()
        || (on_top.domain().1 - on_side.domain().1).abs() > tol.parametric()
    {
        on_side = on_side.segment(on_top.domain(), tol)?;
    }
    crate::bridge::compatible(&mut on_top, &mut on_side, tol)?;
    let knots = on_top.knots().clone();
    let rows = on_top.control_points().len();
    // The rulings run on past both contacts by half their length, as a
    // swept bevel's band does, so a later blend reading the bevel's
    // surface past its edges finds it there; the contacts stay exact, the
    // surface linear across.
    let mut grid = Vec::with_capacity(2 * rows);
    for (a, b) in on_top.control_points().iter().zip(on_side.control_points()) {
        let (a, b) = (a.point(), b.point());
        grid.push(a + (a - b) * BEVEL_RUN_ON);
        grid.push(b + (b - a) * BEVEL_RUN_ON);
    }
    let bevel = ogeom_geom::BSplineSurface::new(
        knots,
        ogeom_math::KnotVector::clamped_uniform(1, 2)?,
        &ogeom_math::ControlGrid::new(grid, rows, 2)?,
        tol,
    )?;
    let mut contacts = [Curve::from(on_top.clone()), Curve::from(on_top)];
    contacts[wall] = Curve::from(on_side);
    Ok(Some(Strip {
        edge,
        faces,
        contacts,
        surface: bevel.into(),
        section: Section::Line,
    }))
}

/// The cone (a cylinder where the radii agree) about `frame`'s axis
/// through the ring of `radius` in its plane and the ring `(rho, dz)`
/// `dz` along its axis, its chart kept clear of its apex.
fn coaxial_cone(
    frame: Frame,
    radius: f64,
    (rho, dz): (f64, f64),
    tol: Tolerances,
) -> OgeomResult<SurfaceGeometry> {
    let slope = (rho - radius) / dz;
    let rows = (dz.min(0.0), dz.max(0.0));
    let depth = dz.abs();
    if slope.atan().abs() <= tol.angular() {
        let pad = 0.1 * depth;
        return Ok(ogeom_geom::CylinderSurface::new(
            ogeom_math::Cylinder::new(frame, radius, tol)?,
            (rows.0 - pad, rows.1 + pad),
        )?
        .into());
    }
    let apex = -radius / slope;
    let clear = (apex - rows.0).abs().min((apex - rows.1).abs());
    let pad = (0.1 * depth).min(0.5 * clear);
    let cone = ogeom_math::Cone::new(frame, radius, slope.atan(), tol)?;
    Ok(ogeom_geom::ConeSurface::new(cone, (rows.0 - pad, rows.1 + pad))?.into())
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
    /// A contact on a consumed face, between two of the sheet's points,
    /// shared by the strips meeting along it.
    Shared(usize, usize, usize),
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
#[cfg_attr(
    not(test),
    allow(
        dead_code,
        reason = "the round sections on curved seats assemble through it from M3 on"
    )
)]
pub(crate) fn assemble(
    model: &mut Model,
    curved: &Curved,
    bevels: &[Strip],
    blended: &FastSet<Id>,
    corners: &FastMap<(usize, Id), Corner>,
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    assemble_consuming(
        model,
        curved,
        bevels,
        blended,
        corners,
        &FastSet::default(),
        tol,
    )
}

/// [`assemble`] where the strips consume the faces of `consumed` whole:
/// each such face is not rebuilt, its strips meet each other along their
/// contacts on it (given as one curve, the ridge), its corners are where
/// that ridge stands at each vertex's station, and its sharp edges go.
///
/// # Errors
///
/// As [`assemble`].
pub(crate) fn assemble_consuming(
    model: &mut Model,
    curved: &Curved,
    bevels: &[Strip],
    blended: &FastSet<Id>,
    corners: &FastMap<(usize, Id), Corner>,
    consumed: &FastSet<usize>,
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
        let kept: Vec<usize> = faces_of(curved, *e)
            .into_iter()
            .filter(|fi| !consumed.contains(fi))
            .collect();
        if blended.contains(e) || kept.is_empty() {
            continue;
        }
        for v in [edge.ends.0, edge.ends.1] {
            let at = curved.points[&v];
            let mut moved: Vec<Point> = Vec::new();
            for fi in kept.iter().copied() {
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
    let mut facets: Vec<(Id, Vec<Point>)> = Vec::new();
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
            _ => {
                // Three bevels or more, every edge there bevelled by a
                // segment: the facet through the faces' corners, its sides
                // the bevels' ends, where those corners lie on one plane.
                let degree = curved
                    .edges
                    .values()
                    .filter(|c| c.ends.0 == *v || c.ends.1 == *v)
                    .count();
                let Some(ends) = sections.get(v) else {
                    return Ok(None);
                };
                if degree != at_v.len()
                    || ends
                        .iter()
                        .any(|(_, _, bi)| !matches!(bevels[*bi].section, Section::Line))
                {
                    return Ok(None);
                }
                let Some(cycle) = chained(ends, tol) else {
                    return Ok(None);
                };
                facets.push((*v, cycle));
            }
        }
    }
    // Each facet's plane through its corners, facing out, and its corners
    // in turn about that side.
    let mut facet_planes: FastMap<Id, SurfaceGeometry> = FastMap::default();
    let mut facet_rings: Vec<(Id, Vec<Point>)> = Vec::with_capacity(facets.len());
    let mut patches: FastMap<Id, SurfaceGeometry> = FastMap::default();
    let mut polygons: Vec<(Id, Vec<Point>)> = Vec::new();
    let mut middles: FastMap<Id, Point> = FastMap::default();
    for (v, cycle) in &facets {
        let mut out = Vector::ZERO;
        for (e, c) in &curved.edges {
            if c.ends.0 == *v || c.ends.1 == *v {
                for fi in faces_of(curved, *e) {
                    out += outward_at(curved, fi, curved.points[v]);
                }
            }
        }
        let mut newell = Vector::ZERO;
        for k in 0..cycle.len() {
            let (p, q) = (
                cycle[k].to_vector(),
                cycle[(k + 1) % cycle.len()].to_vector(),
            );
            newell += p.cross(q);
        }
        let m = newell.magnitude();
        if m <= tol.confusion() * tol.confusion() {
            return Ok(None);
        }
        let mut cycle = cycle.clone();
        if newell.dot(out) < 0.0 {
            cycle.reverse();
        }
        let normal = if newell.dot(out) < 0.0 {
            -newell
        } else {
            newell
        } / m;
        if cycle
            .iter()
            .any(|p| (*p - cycle[0]).dot(normal).abs() > tol.confusion())
        {
            // Off one plane: the patch through the corners, its sides the
            // segments between them, where each bevel there is a plane and
            // its end that segment; at four corners the bilinear patch.
            let curved_here = sections[v]
                .iter()
                .any(|(_, _, bi)| !matches!(bevels[*bi].surface, SurfaceGeometry::Plane(_)));
            if curved_here {
                // A curved bevel there: its end is where its surface meets
                // the plane through its two corners and the hole's middle,
                // and the patch takes that curve for its side.
                #[expect(clippy::cast_precision_loss, reason = "a handful of corners")]
                let middle = Point::from_vector(
                    cycle
                        .iter()
                        .fold(Vector::ZERO, |sum, p| sum + p.to_vector())
                        / cycle.len() as f64,
                );
                middles.insert(*v, middle);
                polygons.push((*v, cycle));
            } else if cycle.len() == 4 {
                patches.insert(*v, bilinear(&cycle, tol)?);
                facet_rings.push((*v, cycle));
            } else {
                // Five corners or more: their Wachspress blend, built on
                // the sheet's edges once it stands.
                polygons.push((*v, cycle));
            }
            continue;
        }
        let frame = Frame::new(
            cycle[0],
            Direction::new(normal, tol)?,
            Direction::new(cycle[1] - cycle[0], tol)?,
            tol,
        )?;
        facet_planes.insert(*v, PlaneSurface::new(Plane::new(frame)).into());
        facet_rings.push((*v, cycle));
    }
    // What a strip's end at `v`, from `p` to `q`, lies on: the vertex's
    // facet, or a plane face there holding both, where the strip ends flush
    // on it, or the other strip sharing the end where either is a plane,
    // or the plane through both and along an edge left sharp there, where
    // two strips meet across it mirrored. Where the
    // strips are not mirrored, the meeting found on that plane is off the
    // other strip and its face is not built.
    let cut_at = |strip: &Strip, v: Id, p: Point, q: Point| -> Option<SurfaceGeometry> {
        if let Some(plane) = facet_planes.get(&v) {
            return Some(plane.clone());
        }
        if let Some(middle) = middles.get(&v) {
            let normal = Direction::new((q - p).cross(*middle - p), tol).ok()?;
            let plane =
                Plane::new(Frame::new(p, normal, Direction::new(q - p, tol).ok()?, tol).ok()?);
            return Some(PlaneSurface::new(plane).into());
        }
        let at_v = || {
            curved
                .edges
                .iter()
                .filter(move |(_, c)| c.ends.0 == v || c.ends.1 == v)
        };
        let flush_plane = at_v()
            .flat_map(|(e, _)| faces_of(curved, *e))
            .filter(|fi| !strip.faces.contains(fi))
            .find_map(|fi| match &curved.faces[fi].surface {
                SurfaceGeometry::Plane(plane) => {
                    let pl = plane.plane();
                    (pl.signed_distance_to(p).abs() <= tol.confusion()
                        && pl.signed_distance_to(q).abs() <= tol.confusion())
                    .then(|| curved.faces[fi].surface.clone())
                }
                _ => None,
            });
        if flush_plane.is_some() {
            return flush_plane;
        }
        // Another strip sharing this end: where either of the two is a
        // plane, they meet along the other cut by it.
        let close = |a: Point, b: Point| a.distance(b) <= tol.confusion();
        let other = sections.get(&v).and_then(|ends| {
            ends.iter().find_map(|(x, y, bi)| {
                let same = (close(*x, p) && close(*y, q)) || (close(*x, q) && close(*y, p));
                (same && bevels[*bi].edge != strip.edge).then_some(&bevels[*bi])
            })
        });
        if let Some(other) = other
            && (matches!(other.surface, SurfaceGeometry::Plane(_))
                || matches!(strip.surface, SurfaceGeometry::Plane(_)))
        {
            return Some(other.surface.clone());
        }
        let at = curved.points[&v];
        at_v()
            .filter(|(e, _)| !blended.contains(*e))
            .find_map(|(_, c)| {
                let t = if c.ends.0 == v { c.range.0 } else { c.range.1 };
                let along = c.curve.d1_at(t, tol).ok()?;
                let normal = Direction::new(along.cross(q - p), tol).ok()?;
                let plane =
                    Plane::new(Frame::new(p, normal, Direction::new(q - p, tol).ok()?, tol).ok()?);
                (plane.signed_distance_to(at).abs() <= tol.confusion())
                    .then(|| PlaneSurface::new(plane).into())
            })
    };
    let section = |points: &mut Points,
                   edges: &mut Vec<SheetCurve>,
                   made: &mut FastMap<Key, usize>,
                   strip: &Strip,
                   cut: Option<&SurfaceGeometry>,
                   p: Point,
                   q: Point|
     -> OgeomResult<(usize, bool)> {
        let (i, j) = (points.index(p), points.index(q));
        let key = Key::Section(i.min(j), i.max(j));
        if let Some(k) = made.get(&key) {
            return Ok((*k, edges[*k].ends.0 == i));
        }
        let (curve, range, ends): (Curve, (f64, f64), (usize, usize)) = match &strip.section {
            Section::Line => {
                // Two planes meet in the segment; a plane and a curved
                // surface in the curved one cut by the plane.
                let met = match (&strip.surface, cut) {
                    (_, None) | (SurfaceGeometry::Plane(_), Some(SurfaceGeometry::Plane(_))) => {
                        None
                    }
                    (SurfaceGeometry::Plane(_), Some(cut)) => {
                        met_between(cut, &strip.surface, p, q, tol)?
                    }
                    (surface, Some(cut)) => met_between(surface, cut, p, q, tol)?,
                };
                match met {
                    Some((curve, range, p_low)) => {
                        (curve, range, if p_low { (i, j) } else { (j, i) })
                    }
                    None => {
                        let (a, b) = (points.points[i.min(j)], points.points[i.max(j)]);
                        let line = LineCurve::segment(a, b, tol)?;
                        let range = line.domain();
                        (line.into(), range, (i.min(j), i.max(j)))
                    }
                }
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
        if consumed.contains(&fi) {
            // Not rebuilt; a placeholder keeps the faces' indices.
            rings_of.push(SheetRings {
                surface: face.surface.clone(),
                rings: Vec::new(),
                reversed: face.reversed,
            });
            continue;
        }
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
                    let strip = &bevels[*bi];
                    if !flush(face, &strip.section, p, start.0, tol)? {
                        return Ok(None);
                    }
                    let (k, along) = section(
                        &mut points,
                        &mut edges,
                        &mut made,
                        strip,
                        Some(&face.surface),
                        p,
                        start.0,
                    )?;
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
                let strip = &bevels[*bi];
                if !flush(face, &strip.section, p, q, tol)? {
                    return Ok(None);
                }
                let (k, along) = section(
                    &mut points,
                    &mut edges,
                    &mut made,
                    strip,
                    Some(&face.surface),
                    p,
                    q,
                )?;
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
        // A contact's piece: the face's own, or, on a consumed face, made
        // here between the corners at the edge's ends, shared by every
        // strip meeting along it there.
        let piece_of = |k: usize,
                        points: &mut Points,
                        edges: &mut Vec<SheetCurve>,
                        made: &mut FastMap<Key, usize>|
         -> OgeomResult<Option<usize>> {
            let fi = b.faces[k];
            if let Some(found) = made.get(&Key::Contact(fi, b.edge)) {
                return Ok(Some(*found));
            }
            let (Some(a), Some(z)) = (corners.get(&(fi, e.ends.0)), corners.get(&(fi, e.ends.1)))
            else {
                return Ok(None);
            };
            let contact = &b.contacts[k];
            let (ta, tz) = (
                param_on(contact, contact.domain(), a.point, tol)?,
                param_on(contact, contact.domain(), z.point, tol)?,
            );
            let (i, j) = (points.index(a.point), points.index(z.point));
            let key = Key::Shared(fi, i.min(j), i.max(j));
            if let Some(found) = made.get(&key).copied() {
                made.insert(Key::Contact(fi, b.edge), found);
                return Ok(Some(found));
            }
            let (range, against) = piece(contact, ta, tz, true, tol);
            let ends = if against { (j, i) } else { (i, j) };
            edges.push(SheetCurve {
                ends,
                curve: contact.clone(),
                range,
            });
            made.insert(key, edges.len() - 1);
            made.insert(Key::Contact(fi, b.edge), edges.len() - 1);
            Ok(Some(edges.len() - 1))
        };
        let (Some(k0), Some(k1)) = (
            piece_of(0, &mut points, &mut edges, &mut made)?,
            piece_of(1, &mut points, &mut edges, &mut made)?,
        ) else {
            return Ok(None);
        };
        // Walk the contact on a face rebuilt against the way that face walks
        // it, as it faces out; then on to the other contact.
        let r = usize::from(consumed.contains(&b.faces[0]));
        if consumed.contains(&b.faces[r]) {
            return Ok(None);
        }
        let (kr, ko) = if r == 0 { (k0, k1) } else { (k1, k0) };
        let Some(walked) = face_walks(&rings_of[b.faces[r]], kr) else {
            return Ok(None);
        };
        let sr = !(walked != rings_of[b.faces[r]].reversed);
        let ends = |edges: &Vec<SheetCurve>, k: usize, forward: bool| {
            let (a, z) = edges[k].ends;
            if forward { (a, z) } else { (z, a) }
        };
        let (pr, qr) = ends(&edges, kr, sr);
        // The vertex an end of the strip stands at: the one whose corners
        // it joins.
        let vertex_of = |x: Point| -> Id {
            let near = |v: Id| {
                [b.faces[0], b.faces[1]].iter().any(|f| {
                    corners
                        .get(&(*f, v))
                        .is_some_and(|c| c.point.distance(x) <= tol.confusion())
                })
            };
            if near(e.ends.0) { e.ends.0 } else { e.ends.1 }
        };
        let other = b.faces[1 - r];
        // The other contact: from its corner where the first ended; a
        // closed one the other way round the strip's axis.
        let start_o = corners
            .get(&(other, vertex_of(points.points[qr])))
            .map(|c| points.index(c.point));
        let Some(start_o) = start_o else {
            return Ok(None);
        };
        let so = if edges[ko].ends.0 != edges[ko].ends.1 {
            edges[ko].ends.0 == start_o
        } else {
            let (Some(turn_r), Some(turn_o)) = (
                turning(&b.surface, &edges[kr], tol)?,
                turning(&b.surface, &edges[ko], tol)?,
            ) else {
                return Ok(None);
            };
            (turn_r * if sr { 1.0 } else { -1.0 }) * turn_o < 0.0
        };
        let (po, qo) = ends(&edges, ko, so);
        let mut ring: Vec<(usize, bool)> = vec![(kr, sr)];
        if qr != po {
            let (a, z) = (points.points[qr], points.points[po]);
            let cut = cut_at(b, vertex_of(a), a, z);
            ring.push(section(
                &mut points,
                &mut edges,
                &mut made,
                b,
                cut.as_ref(),
                a,
                z,
            )?);
        }
        ring.push((ko, so));
        if qo != pr {
            let (a, z) = (points.points[qo], points.points[pr]);
            let cut = cut_at(b, vertex_of(a), a, z);
            ring.push(section(
                &mut points,
                &mut edges,
                &mut made,
                b,
                cut.as_ref(),
                a,
                z,
            )?);
        }
        let p0 = pr;
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
    // Each facet: its ring the bevels' ends in turn about the side it
    // faces.
    for (v, cycle) in &facet_rings {
        let mut ring = Vec::with_capacity(cycle.len());
        for k in 0..cycle.len() {
            let (p, q) = (cycle[k], cycle[(k + 1) % cycle.len()]);
            let strip = sections.get(v).and_then(|ends| {
                ends.iter()
                    .find(|(x, y, _)| {
                        (x.distance(p) <= tol.confusion() && y.distance(q) <= tol.confusion())
                            || (y.distance(p) <= tol.confusion()
                                && x.distance(q) <= tol.confusion())
                    })
                    .map(|(_, _, bi)| &bevels[*bi])
            });
            let Some(strip) = strip else {
                return Ok(None);
            };
            ring.push(section(
                &mut points,
                &mut edges,
                &mut made,
                strip,
                facet_planes.get(v),
                p,
                q,
            )?);
        }
        let Some(surface) = facet_planes.get(v).or_else(|| patches.get(v)) else {
            return Ok(None);
        };
        rings_of.push(SheetRings {
            surface: surface.clone(),
            rings: vec![ring],
            reversed: false,
        });
    }
    // A consumed face is not built: the sheet is every other.
    let built_rings: Vec<SheetRings> = rings_of
        .iter()
        .enumerate()
        .filter(|(i, _)| !consumed.contains(i))
        .map(|(_, r)| r.clone())
        .collect();
    let mut sheet = build(model, &points.points, &edges, &built_rings, tol)?;
    // Each corner of five bevels or more off one plane: the patch through
    // its corners on the bevels' end segments, held to the bands' fit
    // target where it is interpolated.
    let fit = crate::marched::band_fit_target(tol);
    let mut polygon_faces: Vec<(Id, Vec<Shape>)> = Vec::with_capacity(polygons.len());
    for (v, cycle) in &polygons {
        let index = |p: Point| {
            points
                .points
                .iter()
                .position(|x| x.distance(p) <= tol.confusion())
        };
        let mut sides = Vec::with_capacity(cycle.len());
        for k in 0..cycle.len() {
            let (p, q) = (cycle[k], cycle[(k + 1) % cycle.len()]);
            let (Some(i), Some(j)) = (index(p), index(q)) else {
                return Ok(None);
            };
            let Some(&e) = made.get(&Key::Section(i.min(j), i.max(j))) else {
                return Ok(None);
            };
            let strip = sections[v].iter().find_map(|(x, y, bi)| {
                let at = |a: Point, b: Point| {
                    a.distance(p) <= tol.confusion() && b.distance(q) <= tol.confusion()
                };
                (at(*x, *y) || at(*y, *x)).then_some(&bevels[*bi])
            });
            let Some(strip) = strip else {
                return Ok(None);
            };
            sides.push(crate::round::transfinite::HoleSide {
                edge: sheet.edges[e].clone(),
                curve: edges[e].curve.clone(),
                range: edges[e].range,
                forward: edges[e].ends.0 == i,
                support: strip.surface.clone(),
            });
        }
        let patch =
            match crate::round::transfinite::corner_patch_through_corners(model, &sides, tol) {
                Ok(patch) => patch,
                Err(e @ ogeom_core::OgeomError::Cancelled) => return Err(e),
                Err(_) => return Ok(None),
            };
        if !(patch.gap <= fit && patch.off <= fit) {
            return Ok(None);
        }
        polygon_faces.push((*v, patch.faces));
    }
    let shell_faces: Vec<Shape> = sheet
        .faces
        .iter()
        .chain(polygon_faces.iter().flat_map(|(_, f)| f))
        .cloned()
        .collect();
    let shell = ogeom_algo::make_shell(model, &shell_faces)?.shape;
    // The sheet's faces again by the rings' indices, a consumed face's slot
    // holding the face's own shape, never used.
    let mut by_index: Vec<Shape> = Vec::with_capacity(rings_of.len());
    let mut next = sheet.faces.iter();
    for i in 0..rings_of.len() {
        match curved.faces.get(i) {
            Some(face) if consumed.contains(&i) => by_index.push(face.shape.clone()),
            _ => by_index.push(
                next.next()
                    .cloned()
                    .unwrap_or_else(|| curved.faces[0].shape.clone()),
            ),
        }
    }
    sheet.faces = by_index;
    let solid = ogeom_algo::make_solid(model, &[shell])?.shape;
    // An edge a face walks twice, once each way (a wall's or a cone
    // bevel's seam), needs both columns of the seam: the builder fits it
    // one pcurve, and a later walk of the face in its chart would not
    // close. Nothing else holds these nodes, so nothing is copied.
    ogeom_algo::join_seam_columns(model, &solid, tol)?;
    // P5: each bevel proved within the material its two faces bound there,
    // to the bands' fit target: under a plane, inside a wall (outside a
    // bore's).
    for (bi, b) in bevels.iter().enumerate() {
        // An extruded wall has no region the certify reads: its bevel's
        // rulings stand in their sections on its material side by
        // construction ([`extruded_crease`]), and only the cap is proved.
        let by_section = |k: usize| {
            matches!(
                curved.faces[b.faces[k]].surface,
                SurfaceGeometry::Extrusion(_)
            ) && matches!(b.surface, SurfaceGeometry::BSpline(_))
        };
        let pieces: Vec<ogeom_algo::HalfSpaces> =
            match (material(curved, b.faces[0]), material(curved, b.faces[1])) {
                (Some(first), Some(second)) => vec![first, second],
                (Some(first), None) if by_section(1) => vec![first],
                (None, Some(second)) if by_section(0) => vec![second],
                _ => return Ok(None),
            };
        let region = ogeom_algo::HalfSpaces::Every(pieces);
        let face = &sheet.faces[curved.faces.len() + bi];
        if ogeom_algo::certify_face_within(model, face, &region, fit, tol)?.is_err() {
            return Ok(None);
        }
    }
    // Each patch within the material every face at its vertex bounds.
    let patch_faces: Vec<(Id, Shape)> = facet_rings
        .iter()
        .enumerate()
        .filter(|(_, (v, _))| patches.contains_key(v))
        .map(|(fi, (v, _))| {
            (
                *v,
                sheet.faces[curved.faces.len() + bevels.len() + fi].clone(),
            )
        })
        .chain(
            polygon_faces
                .iter()
                .flat_map(|(v, faces)| faces.iter().map(move |f| (*v, f.clone()))),
        )
        .collect();
    for (v, face) in &patch_faces {
        let mut round: Vec<usize> = curved
            .edges
            .iter()
            .filter(|(_, c)| c.ends.0 == *v || c.ends.1 == *v)
            .flat_map(|(e, _)| faces_of(curved, *e))
            .collect();
        round.sort_unstable();
        round.dedup();
        let Some(pieces) = round
            .iter()
            .map(|fi| material(curved, *fi))
            .collect::<Option<Vec<_>>>()
        else {
            return Ok(None);
        };
        let region = ogeom_algo::HalfSpaces::Every(pieces);
        if ogeom_algo::certify_face_within(model, face, &region, fit, tol)?.is_err() {
            return Ok(None);
        }
    }
    // History: faces modified, bevelled edges deleted into their bevels,
    // sharp edges modified into their pieces, vertices kept where nothing
    // moved off them.
    let mut history = History::new();
    for (fi, face) in curved.faces.iter().enumerate() {
        if consumed.contains(&fi) {
            history.delete(&face.shape);
        } else {
            history.modify(&face.shape, sheet.faces[fi].clone());
        }
    }
    for (bi, b) in bevels.iter().enumerate() {
        let edge = &curved.shapes[&b.edge];
        history.delete(edge);
        history.generate(edge, sheet.faces[curved.faces.len() + bi].clone());
    }
    for (fi, (v, _)) in facet_rings.iter().enumerate() {
        history.generate(
            &curved.shapes[v],
            sheet.faces[curved.faces.len() + bevels.len() + fi].clone(),
        );
    }
    for (v, faces) in &polygon_faces {
        for face in faces {
            history.generate(&curved.shapes[v], face.clone());
        }
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

/// The bevels' ends at a vertex chained into one cycle of corners, each
/// end joining the corners of two faces; `None` where they do not close one
/// cycle through every corner once.
fn chained(ends: &[(Point, Point, usize)], tol: Tolerances) -> Option<Vec<Point>> {
    let same = |a: Point, b: Point| a.distance(b) <= tol.confusion();
    let mut cycle = vec![ends[0].0, ends[0].1];
    let mut used = vec![false; ends.len()];
    used[0] = true;
    loop {
        let last = *cycle.last()?;
        let next = ends
            .iter()
            .enumerate()
            .find(|(i, (p, q, _))| !used[*i] && (same(*p, last) || same(*q, last)))
            .map(|(i, (p, q, _))| (i, if same(*p, last) { *q } else { *p }));
        let Some((i, far)) = next else {
            break;
        };
        used[i] = true;
        if same(far, cycle[0]) {
            return used.iter().all(|u| *u).then_some(cycle);
        }
        cycle.push(far);
    }
    None
}

/// Which way a sheet edge, walked as stored, turns about a revolved strip's
/// axis: the sign of its tangent's turn at its start, `None` off a cone or
/// a cylinder.
fn turning(
    surface: &SurfaceGeometry,
    edge: &SheetCurve,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let (o, a) = match surface {
        SurfaceGeometry::Cone(c) => {
            let f = c.cone().frame();
            (f.origin(), f.z().vector())
        }
        SurfaceGeometry::Cylinder(c) => {
            let axis = c.cylinder().axis();
            (axis.location, axis.direction.vector())
        }
        _ => return Ok(None),
    };
    let t = edge.range.0;
    let (p, d) = (edge.curve.point_at(t, tol)?, edge.curve.d1_at(t, tol)?);
    Ok(Some((p - o).cross(d).dot(a).signum()))
}

/// A piece of a curve: the curve, its range, and whether the first point
/// asked for stands at the range's low end.
type Met = (Curve, (f64, f64), bool);

/// The piece of where `surface` meets the plane `cut` between `p` and `q`,
/// the shorter way along it: its curve, its range, and whether `p` stands
/// at the range's low end. `None` where the two have no closed-form meeting
/// through both points.
fn met_between(
    surface: &SurfaceGeometry,
    cut: &SurfaceGeometry,
    p: Point,
    q: Point,
    tol: Tolerances,
) -> OgeomResult<Option<Met>> {
    let Ok(ogeom_intersect::surface::Meeting::Along(curves)) =
        ogeom_intersect::surface::surface_surface(surface, cut, tol)
    else {
        return Ok(None);
    };
    for curve in curves {
        let domain = curve.domain();
        if !(domain.0.is_finite() && domain.1.is_finite()) {
            continue;
        }
        let (tp, tq) = (
            param_on(&curve, domain, p, tol)?,
            param_on(&curve, domain, q, tol)?,
        );
        if curve.point_at(tp, tol)?.distance(p) > tol.confusion()
            || curve.point_at(tq, tol)?.distance(q) > tol.confusion()
        {
            continue;
        }
        let period = curve.is_periodic().then_some(domain.1 - domain.0);
        let (lo, hi, p_low) = if tp <= tq {
            (tp, tq, true)
        } else {
            (tq, tp, false)
        };
        return Ok(Some(match period {
            Some(period) if hi - lo > period / 2.0 => (curve, (hi, lo + period), !p_low),
            _ => (curve, (lo, hi), p_low),
        }));
    }
    Ok(None)
}

/// The material face `face` bounds, as a region the certify reads: the
/// side of its plane under its outward normal, the side of its cylinder,
/// cone or sphere the material lies on (inside a wall, outside a bore), or
/// the inside of its torus's tube; `None` on any other surface or side.
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
        SurfaceGeometry::Sphere(sphere) => {
            let sphere = sphere.sphere();
            Some(ogeom_algo::HalfSpaces::Ball {
                centre: sphere.centre(),
                radius: sphere.radius(),
                inside,
            })
        }
        SurfaceGeometry::Torus(torus) => {
            let torus = torus.torus();
            let frame = torus.frame();
            Some(ogeom_algo::HalfSpaces::Tube {
                centre: frame.origin(),
                axis: frame.z().vector(),
                major: torus.major_radius(),
                minor: torus.minor_radius(),
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

/// A bevel surface's own normal at `p`: a plane's, a cone's or a
/// cylinder's off its axis there, or a B-spline's at the foot of `p`.
fn natural_normal(surface: &SurfaceGeometry, p: Point, tol: Tolerances) -> OgeomResult<Vector> {
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
        SurfaceGeometry::BSpline(_) => {
            let foot = ogeom_geom::project::project_on_surface(surface, p, 16, tol)?;
            let (u, v) = foot.parameters;
            let (du, dv) = ogeom_geom::Surface::d1_at(surface, u, v, tol)?;
            du.cross(dv)
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
