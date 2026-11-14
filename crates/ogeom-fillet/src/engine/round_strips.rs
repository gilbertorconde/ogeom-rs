//! Rounds on solids with curved faces through the trimmed assembly: each
//! edge's strip in closed form, handed to [`assemble`] with the corners
//! where the contacts cross ([`curved_corners`]).
//!
//! A straight edge between two planes, or between a plane and a cylinder
//! along one of its rulings (the plane holding the cylinder's direction),
//! rounds on the cylinder its ball sweeps: the ball's centre stands a radius
//! in from both faces, on the line where the plane moved in meets the other
//! plane moved in, or the cylinder moved in or out (by the side its
//! material is on), out of both at a concave edge. A rim between a planar
//! cap and a coaxial cylinder or cone rounds on the torus its ball sweeps.

use ogeom_algo::Built;
use ogeom_core::{FastMap, FastSet, OgeomResult, Tolerances};
use ogeom_geom::{
    CircleCurve, Curve, Curve3d as _, CylinderSurface, LineCurve, SurfaceGeometry, TorusSurface,
};
use ogeom_math::{Circle, Cylinder, Direction, Frame, Point, Torus, Vector};
use ogeom_topo::{Model, Shape};

use super::curved::{Curved, CurvedEdge, read};
use super::setbacks::{corner_for, curved_corners, roll_overs, rolls_over};
use super::trimmed::{RolledRun, Section, Strip, assemble_rolled, faces_of, outward};
use super::{Id, id};

/// The rounds of `rounds` (each an edge and its radius) on `solid` with
/// curved faces, or `None` where an edge, a vertex or the result lies
/// outside the classes built here.
///
/// # Errors
///
/// Where the operation is cancelled.
pub(crate) fn round(
    model: &mut Model,
    solid: &Shape,
    rounds: &[(Shape, f64)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let Some(curved) = read(model, solid, tol)? else {
        return Ok(None);
    };
    // A face's own seam (an edge one face walks both ways) names nothing.
    let mut asked: Vec<(Id, &Shape, f64)> = rounds
        .iter()
        .map(|(e, r)| (id(e), e, *r))
        .filter(|(e, _, _)| faces_of(&curved, *e).len() != 1)
        .collect();
    if asked.is_empty() {
        return Ok(None);
    }
    asked.sort_unstable_by_key(|(e, _, _)| *e);
    if asked.windows(2).any(|w| w[0].0 == w[1].0) {
        return Ok(None);
    }
    // In order of where each edge stands, its key only between two that
    // stand alike: a dissolved solid's keys change from one reading to the
    // next, and the strips' order is the result's face order.
    asked.sort_by_cached_key(|(e, _, _)| {
        let c = &curved.edges[e];
        let at = c
            .curve
            .point_at(0.5 * (c.range.0 + c.range.1), tol)
            .map_or([u64::MAX; 3], |p| {
                [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
            });
        (at, *e)
    });
    // P1.
    let mut strips: Vec<Strip> = Vec::with_capacity(asked.len());
    let mut boxes: Vec<(usize, ogeom_algo::HalfSpaces)> = Vec::new();
    for (edge, _, radius) in &asked {
        let Some((strip, clear)) = strip_of(&curved, *edge, *radius, tol)? else {
            return Ok(None);
        };
        if let Some(clear) = clear {
            boxes.push((strips.len(), clear));
        }
        strips.push(strip);
    }
    let mut contacts: FastMap<(usize, Id), Curve> = FastMap::default();
    for s in &strips {
        for k in 0..2 {
            contacts.insert((s.faces[k], s.edge), s.contacts[k].clone());
        }
    }
    // P2.
    let Some(corners) = curved_corners(&curved, &contacts, tol)? else {
        return Ok(None);
    };
    // A contact leaving its face rolls over onto what bounds it there: the
    // ball rests on the edge it crosses until it comes back (ADR 0001,
    // section 13, question 8), each run read from the crossing pairs.
    let crossings = roll_overs(&curved, &contacts, &corners, tol)?;
    let mut rolled: Vec<RolledRun> = Vec::new();
    if rolls_over(&crossings) {
        let mut k = 0;
        while k < crossings.len() {
            let (enter, Some(leave)) = (&crossings[k], crossings.get(k + 1)) else {
                return Ok(None);
            };
            if (enter.face, enter.edge) != (leave.face, leave.edge) {
                return Ok(None);
            }
            let Some(strip) = strips.iter().position(|s| s.edge == enter.edge) else {
                return Ok(None);
            };
            let s = &strips[strip];
            let Some(side) = s.faces.iter().position(|f| *f == enter.face) else {
                return Ok(None);
            };
            let Section::Arc { spine } = &s.section else {
                return Ok(None);
            };
            let radius = asked
                .iter()
                .find(|(e, _, _)| *e == s.edge)
                .map_or(0.0, |(_, _, r)| *r);
            let run = match super::seats::edge_run(
                model,
                &curved,
                spine,
                radius,
                enter.face,
                s.faces[1 - side],
                enter,
                leave,
                tol,
            )? {
                Ok(run) => run,
                // A ball that would pass a vertex of the edge it rests on, or
                // leave that edge at its end, takes a vertex for a seat, which
                // is refused by name (ADR 0001, section 13, question 8).
                Err(
                    why @ (super::seats::NotSeated::PastAVertex
                    | super::seats::NotSeated::OffTheEdge),
                ) => ogeom_core::ogeom_bail!(Construction, "{}", why.name()),
                // A seat change not built yet goes to the current path.
                Err(_) => return Ok(None),
            };
            // Nothing else of the solid in the box the run's piece can take
            // up: its arcs, the edge it rests on and the strip's edge there.
            if !run_clear(model, &curved, s, &run, tol)? {
                return Ok(None);
            }
            rolled.push(RolledRun {
                strip,
                side,
                enter: enter.clone(),
                leave: leave.clone(),
                run,
            });
            k += 2;
        }
    }
    // A rim's round reaching a face it does not blend (a cavity's wall
    // under a thin cap) is not built here: every face but the strip's two
    // proved outside the box the round can take up.
    let fit = crate::marched::band_fit_target(tol);
    for (k, clear) in &boxes {
        // The faces met at the strip's ends close its corners there, and
        // those of the edges its ball rests on hold its run.
        let e = &curved.edges[&strips[*k].edge];
        let local: Vec<usize> = (0..curved.faces.len())
            .filter(|fi| {
                curved.faces[*fi].loops.iter().flatten().any(|(x, _)| {
                    let c = &curved.edges[x];
                    [c.ends.0, c.ends.1]
                        .iter()
                        .any(|v| *v == e.ends.0 || *v == e.ends.1)
                })
            })
            .chain(
                rolled
                    .iter()
                    .filter(|r| r.strip == *k)
                    .flat_map(|r| r.run.over.iter().flat_map(|o| faces_of(&curved, *o))),
            )
            .collect();
        for (fi, face) in curved.faces.iter().enumerate() {
            if strips[*k].faces.contains(&fi) || local.contains(&fi) {
                continue;
            }
            if ogeom_algo::certify_face_within(model, &face.shape, clear, fit, tol)?.is_err() {
                return Ok(None);
            }
        }
    }
    // A straight strip's cylinder over its own stations, from the least of
    // its corners' to the greatest: a surface reaching past its face's trim
    // is one more thing a later operation reads. An end whose two corners
    // stand at different stations (a leaning cap) is widened by a radius,
    // for its arc between them.
    for strip in &mut strips {
        let SurfaceGeometry::Cylinder(c) = &strip.surface else {
            continue;
        };
        let cyl = c.cylinder();
        let frame = cyl.frame();
        let (o, u) = (frame.origin(), frame.z().vector());
        let e = &curved.edges[&strip.edge];
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for (end, v) in [e.ends.0, e.ends.1].into_iter().enumerate() {
            let (Some(c0), Some(c1)) = (
                corner_for(
                    &corners,
                    &curved,
                    strip.faces[0],
                    v,
                    strip.edge,
                    end == 0,
                    tol,
                )?,
                corner_for(
                    &corners,
                    &curved,
                    strip.faces[1],
                    v,
                    strip.edge,
                    end == 0,
                    tol,
                )?,
            ) else {
                return Ok(None);
            };
            let (s0, s1) = ((c0.point - o).dot(u), (c1.point - o).dot(u));
            let lean = if (s0 - s1).abs() > tol.confusion() {
                cyl.radius()
            } else {
                0.0
            };
            lo = lo.min(s0.min(s1) - lean);
            hi = hi.max(s0.max(s1) + lean);
        }
        if hi - lo <= tol.confusion() {
            return Ok(None);
        }
        let start = Frame::new(o + u * lo, frame.z(), frame.x(), tol)?;
        strip.surface =
            CylinderSurface::new(Cylinder::new(start, cyl.radius(), tol)?, (0.0, hi - lo))?.into();
        // Each contact the segment between its face's two corners, and the
        // spine over the surface's span: the lines run past the edge only to
        // find the corners, and a curve's domain past its edge is one more
        // thing a later operation reads.
        // Each keeps its own line and parameters (the corners hold theirs on
        // it), its domain cut to them.
        let trimmed = |curve: &Curve, from: Point, to: Point| -> OgeomResult<Option<Curve>> {
            let Curve::Line(line) = curve else {
                return Ok(None);
            };
            let axis = line.axis();
            let t = |p: Point| (p - axis.location).dot(axis.direction.vector());
            let (a, b) = (t(from), t(to));
            Ok(Some(LineCurve::over(axis, a.min(b), a.max(b))?.into()))
        };
        for k in 0..2 {
            let (Some(p), Some(q)) = (
                corner_for(
                    &corners,
                    &curved,
                    strip.faces[k],
                    e.ends.0,
                    strip.edge,
                    true,
                    tol,
                )?,
                corner_for(
                    &corners,
                    &curved,
                    strip.faces[k],
                    e.ends.1,
                    strip.edge,
                    false,
                    tol,
                )?,
            ) else {
                return Ok(None);
            };
            let (p, q) = (p.point, q.point);
            let Some(cut) = trimmed(&strip.contacts[k], p, q)? else {
                return Ok(None);
            };
            strip.contacts[k] = cut;
        }
        let Section::Arc { spine } = &strip.section else {
            return Ok(None);
        };
        let Some(spine) = trimmed(spine, o + u * lo, o + u * hi)? else {
            return Ok(None);
        };
        strip.section = Section::Arc { spine };
    }
    let blended: FastSet<Id> = strips.iter().map(|s| s.edge).collect();
    // P4 and P5; what the assembly or the solid's check refuses is handed
    // back.
    let built = match assemble_rolled(model, &curved, &strips, &blended, &corners, &rolled, tol) {
        Ok(Some(built)) => built,
        Ok(None) => return Ok(None),
        Err(e @ ogeom_core::OgeomError::Cancelled) => return Err(e),
        Err(_) => return Ok(None),
    };
    match super::guards::solid_fits(model, &built.shape, tol) {
        Ok(()) => Ok(Some(built)),
        Err(e @ ogeom_core::OgeomError::Cancelled) => Err(e),
        Err(_) => Ok(None),
    }
}

/// The strip of `edge` rounded at `radius`, with, for a rim, the outside of
/// the box its round can take up; `None` where its seat is not one built
/// here.
fn strip_of(
    curved: &Curved,
    edge: Id,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Strip, Option<ogeom_algo::HalfSpaces>)>> {
    let Some(e) = curved.edges.get(&edge) else {
        return Ok(None);
    };
    let [f0, f1] = faces_of(curved, edge)[..] else {
        return Ok(None);
    };
    if !(radius.is_finite() && radius > tol.confusion()) {
        return Ok(None);
    }
    match &e.curve {
        Curve::Line(_) => Ok(straight(curved, edge, e, [f0, f1], radius, tol)?),
        Curve::Circle(circle) => {
            Ok(rim(curved, edge, circle, [f0, f1], radius, tol)?.map(|(s, c)| (s, Some(c))))
        }
        _ => Ok(None),
    }
}

/// A straight edge's strip between two planes, a plane and a cylinder
/// along its ruling, or two cylinders along their rulings, their axes
/// parallel: the ball's centre a radius in from both faces at a
/// convex edge, out from both at a concave one. A concave strip's fill also
/// gives the region every other face must keep to: past the edge's ends or
/// farther from the spine than the edge stands, which holds the whole fill
/// (from the ball's centre the fill's boundary is farthest at the edge:
/// along the plane, and along the cylinder, whose nearest point is the
/// foot).
fn straight(
    curved: &Curved,
    edge: Id,
    e: &CurvedEdge,
    faces: [usize; 2],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Strip, Option<ogeom_algo::HalfSpaces>)>> {
    let (a, b) = (curved.points[&e.ends.0], curved.points[&e.ends.1]);
    let length = a.distance(b);
    if length <= tol.confusion() {
        return Ok(None);
    }
    let u = (b - a) / length;
    let mid = a + (b - a) * 0.5;
    let (Some(n0), Some(n1)) = (
        outward(curved, faces[0], mid),
        outward(curved, faces[1], mid),
    ) else {
        return Ok(None);
    };
    // Into each face from the edge, square to it: along the face's other
    // edge leaving one of the edge's ends.
    let (Some(into0), Some(into1)) = (
        into_face(curved, faces[0], edge, u, tol)?,
        into_face(curved, faces[1], edge, u, tol)?,
    ) else {
        return Ok(None);
    };
    let into = [into0, into1];
    // Convex where the second face leaves under the first's tangent plane;
    // the ball's centre a radius off each face, in at a convex edge and out
    // at a concave one.
    let convex = n0.dot(into1) < 0.0;
    let side = if convex { -1.0 } else { 1.0 };
    // Each face as the ball's centre sees it: a plane's offset, or a
    // cylinder's axis and the radius its centre stands at.
    enum Offset {
        Plane(Vector),
        Cylinder(Point, f64),
    }
    let mut offsets = Vec::with_capacity(2);
    for (k, fi) in faces.iter().enumerate() {
        let n = if k == 0 { n0 } else { n1 };
        match &curved.faces[*fi].surface {
            SurfaceGeometry::Plane(_) => offsets.push(Offset::Plane(n)),
            SurfaceGeometry::Cylinder(c) => {
                let cyl = c.cylinder();
                let axis = cyl.axis();
                let (o, ax) = (axis.location, axis.direction.vector());
                if u.cross(ax).magnitude() > tol.angular() {
                    return Ok(None);
                }
                // The axis where it crosses the edge's section through `mid`.
                let o = o + ax * (mid - o).dot(ax);
                // Material inside the cylinder where its outward normal
                // points off the axis.
                let inside = n.dot(mid - o) > 0.0;
                let rho = if inside == convex {
                    cyl.radius() - radius
                } else {
                    cyl.radius() + radius
                };
                if rho <= tol.confusion() {
                    return Ok(None);
                }
                offsets.push(Offset::Cylinder(o, rho));
            }
            _ => return Ok(None),
        }
    }
    // The ball's centre in the section through `mid`.
    let centre = match (&offsets[0], &offsets[1]) {
        (Offset::Plane(p), Offset::Plane(q)) => {
            let g = p.dot(*q);
            if (1.0 + g).abs() <= tol.angular() {
                return Ok(None);
            }
            let s = side * radius / (1.0 + g);
            mid + (*p + *q) * s
        }
        (Offset::Plane(n), Offset::Cylinder(o, rho))
        | (Offset::Cylinder(o, rho), Offset::Plane(n)) => {
            // Along the plane moved in: `mid - n r + t m`.
            let m = n.cross(u);
            let base = mid + *n * (side * radius);
            let w = base - *o;
            let (qb, qc) = (2.0 * w.dot(m), w.dot(w) - rho * rho);
            let disc = qb * qb - 4.0 * qc;
            if disc < 0.0 {
                return Ok(None);
            }
            let root = disc.sqrt();
            let pick = [(-qb - root) / 2.0, (-qb + root) / 2.0]
                .into_iter()
                .map(|t| base + m * t)
                .min_by(|x, y| x.distance(mid).total_cmp(&y.distance(mid)));
            let Some(c) = pick else {
                return Ok(None);
            };
            c
        }
        (Offset::Cylinder(o0, r0), Offset::Cylinder(o1, r1)) => {
            // Two walls along the edge (a boss's lobes fused side by side):
            // where the two circles their axes and offset radii draw in the
            // section cross, the crossing nearer the edge.
            let d = *o1 - *o0;
            let gap = d.magnitude();
            if gap <= tol.confusion() {
                return Ok(None);
            }
            let e = d / gap;
            let along = (r0 * r0 - r1 * r1 + gap * gap) / (2.0 * gap);
            let h2 = r0 * r0 - along * along;
            if h2 < 0.0 {
                return Ok(None);
            }
            let w = u.cross(e);
            let base = *o0 + e * along;
            let h = h2.sqrt();
            let Some(c) = [base + w * h, base - w * h]
                .into_iter()
                .min_by(|x, y| x.distance(mid).total_cmp(&y.distance(mid)))
            else {
                return Ok(None);
            };
            c
        }
    };
    // Each contact, and its standing on its own face's side of the edge.
    let mut feet = [Point::ORIGIN; 2];
    for (k, fi) in faces.iter().enumerate() {
        let foot = match &offsets[k] {
            Offset::Plane(n) => centre - *n * (side * radius),
            Offset::Cylinder(o, _) => {
                let cyl_r = match &curved.faces[*fi].surface {
                    SurfaceGeometry::Cylinder(c) => c.cylinder().radius(),
                    _ => return Ok(None),
                };
                let off = centre - *o;
                *o + off / off.magnitude() * cyl_r
            }
        };
        if (foot - mid).dot(into[k]) <= 0.0 {
            return Ok(None);
        }
        feet[k] = foot;
    }
    // Lines run on past both ends; the assembly trims them.
    let run = u * (2.0 * radius);
    // Each line's parameter nought at the edge's start, as the edge's own.
    let line = |p: Point| -> OgeomResult<Curve> {
        let at = p + (a - mid);
        let reach = run.magnitude();
        Ok(LineCurve::over(
            ogeom_math::Axis {
                location: at,
                direction: Direction::new(u, tol)?,
            },
            -reach,
            length + reach,
        )?
        .into())
    };
    let contacts = [line(feet[0])?, line(feet[1])?];
    // The surface over the edge and two radii past each end, no farther: a
    // surface reaching far past its face's trim is one more thing a later
    // operation reads. An end the assembly cuts beyond that refuses.
    let reach = 2.0 * radius;
    let start = centre + (a - mid) - u * reach;
    let x = Direction::new(feet[0] - centre, tol)?;
    let frame = Frame::new(start, Direction::new(u, tol)?, x, tol)?;
    let surface: SurfaceGeometry = CylinderSurface::new(
        Cylinder::new(frame, radius, tol)?,
        (0.0, length + 2.0 * reach),
    )?
    .into();
    let clear = (!convex).then(|| {
        let at = centre + (a - mid);
        ogeom_algo::HalfSpaces::Any(vec![
            ogeom_algo::HalfSpaces::Below { normal: u, at: a },
            ogeom_algo::HalfSpaces::Below { normal: -u, at: b },
            ogeom_algo::HalfSpaces::Revolved {
                at,
                axis: u,
                radius: (mid - centre).magnitude(),
                slope: 0.0,
                inside: false,
            },
        ])
    });
    Ok(Some((
        Strip {
            edge,
            faces,
            contacts,
            surface,
            section: Section::Arc {
                spine: line(centre)?,
            },
        },
        clear,
    )))
}

/// How far off a straight edge's line another edge must leave it, as a
/// sine, to tell which side of the edge its face lies.
const LEAVING: f64 = 1e-3;

/// Into `face` from the straight `edge` along `u`, square to it: the
/// face's other edge leaving one of the edge's ends, its tangent there with
/// its part along `u` taken off; where every one runs on along `u` (an arc
/// leaving tangent, a rim continuing a straight edge), on the left of the
/// edge as the face walks it, about the face's outward normal at the edge's
/// middle (every loop, a hole's too, walks with the face on its left: a
/// neighbour's chord would point into a hole beside it). `None` where
/// neither reads.
///
/// # Errors
///
/// As evaluating the curves.
fn into_face(
    curved: &Curved,
    face: usize,
    edge: Id,
    u: Vector,
    tol: Tolerances,
) -> OgeomResult<Option<Vector>> {
    use ogeom_geom::Curve3d as _;
    let ends = curved.edges[&edge].ends;
    for (e, _) in curved.faces[face].loops.iter().flatten() {
        if *e == edge {
            continue;
        }
        let c = &curved.edges[e];
        for (v, at, sign) in [(c.ends.0, c.range.0, 1.0), (c.ends.1, c.range.1, -1.0)] {
            if v != ends.0 && v != ends.1 {
                continue;
            }
            let t = c.curve.d1_at(at, tol)? * sign;
            let square = t - u * t.dot(u);
            if square.magnitude() > LEAVING * t.magnitude() {
                return Ok(Some(square));
            }
        }
    }
    let f = &curved.faces[face];
    let Some(&(_, w)) = f.loops.iter().flatten().find(|(e, _)| *e == edge) else {
        return Ok(None);
    };
    let c = &curved.edges[&edge];
    let mid = c.curve.point_at(0.5 * (c.range.0 + c.range.1), tol)?;
    let Some(n) = outward(curved, face, mid) else {
        return Ok(None);
    };
    let along = if w == f.reversed { -u } else { u };
    let into = n.cross(along);
    let square = into - u * into.dot(u);
    if square.magnitude() <= LEAVING * into.magnitude().max(f64::MIN_POSITIVE) {
        return Ok(None);
    }
    Ok(Some(square))
}

/// The points each other edge of a rim's face is read at, to keep the
/// contact ring clear of it.
const RING_SAMPLES: u32 = 32;

/// A rim's strip, between two faces of revolution about the rim's axis
/// (a plane square to it, a coaxial cylinder or cone): the torus its ball
/// sweeps. In the rim's meridian each face is a straight line leaving the
/// rim's point, read from a point of the face off the rim; the ball's centre
/// stands a radius in from both (out of the material at a concave rim), and
/// its feet on them are the strip's contact rings.
fn rim(
    curved: &Curved,
    edge: Id,
    circle: &CircleCurve,
    faces: [usize; 2],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<(Strip, ogeom_algo::HalfSpaces)>> {
    let c = circle.circle();
    let frame = c.frame();
    let (centre, up) = (frame.origin(), frame.z().vector());
    let r0 = c.radius();
    let x = frame.x();
    let meridian = |p: Point| {
        let off = p - centre;
        let h = off.dot(up);
        ((off - up * h).magnitude(), h)
    };
    let radial = x.vector();
    // Each face's way out of the rim and its outward normal, in (radial,
    // axial) coordinates.
    let mut ways = [(0.0, 0.0); 2];
    let mut outs = [(0.0, 0.0); 2];
    for (k, fi) in faces.iter().enumerate() {
        let face = &curved.faces[*fi];
        let sign = if face.reversed { -1.0 } else { 1.0 };
        let coaxial = |o: Point, a: Vector| {
            a.cross(up).magnitude() <= tol.angular()
                && ((o - centre) - up * (o - centre).dot(up)).magnitude() <= tol.confusion()
        };
        let natural: Vector = match &face.surface {
            SurfaceGeometry::Plane(p) => {
                let n = p.plane().normal().vector();
                if n.cross(up).magnitude() > tol.angular() {
                    return Ok(None);
                }
                n
            }
            SurfaceGeometry::Cylinder(cy) => {
                let axis = cy.cylinder().axis();
                if !coaxial(axis.location, axis.direction.vector())
                    || (cy.cylinder().radius() - r0).abs() > tol.confusion()
                {
                    return Ok(None);
                }
                radial
            }
            SurfaceGeometry::Cone(co) => {
                let cone = co.cone();
                let f = cone.frame();
                let a = f.z().vector();
                if !coaxial(f.origin(), a) {
                    return Ok(None);
                }
                radial - a * cone.half_angle().tan()
            }
            _ => return Ok(None),
        };
        let n = natural * sign;
        let out = (n.dot(radial), n.dot(up));
        let norm = out.0.hypot(out.1);
        if norm <= tol.angular() {
            return Ok(None);
        }
        outs[k] = (out.0 / norm, out.1 / norm);
        // A wall's way out of the rim, read at a point of it off the rim; a
        // plane's is read from the wall below, the plane being the one face
        // that can reach both sides of the circle.
        if !matches!(face.surface, SurfaceGeometry::Plane(_)) {
            // A vertex of the face off the rim's circle (another arc of the
            // rim stands on it).
            let off = face
                .loops
                .iter()
                .flatten()
                .flat_map(|(e, _)| {
                    let c = &curved.edges[e];
                    [c.ends.0, c.ends.1]
                })
                .map(|v| {
                    let (rho, h) = meridian(curved.points[&v]);
                    (rho - r0, h)
                })
                .find(|w| w.0.hypot(w.1) > tol.confusion());
            let Some(way) = off else {
                return Ok(None);
            };
            let length = way.0.hypot(way.1);
            ways[k] = (way.0 / length, way.1 / length);
        }
    }
    let dot = |a: (f64, f64), b: (f64, f64)| a.0 * b.0 + a.1 * b.1;
    let planes: Vec<usize> = (0..2)
        .filter(|k| matches!(curved.faces[faces[*k]].surface, SurfaceGeometry::Plane(_)))
        .collect();
    match planes.as_slice() {
        [] => {}
        [k] => {
            let (k, w) = (*k, 1 - *k);
            // Convex where the wall leaves on the plane's material side; the
            // plane then leaves on the wall's, and the other way round at a
            // concave rim.
            let convex = dot(ways[w], outs[k]) < 0.0;
            let square = (-outs[k].1, outs[k].0);
            let side = dot(square, outs[w]);
            if side.abs() <= tol.angular() {
                return Ok(None);
            }
            let toward = if convex {
                -side.signum()
            } else {
                side.signum()
            };
            ways[k] = (square.0 * toward, square.1 * toward);
        }
        _ => return Ok(None),
    }
    // Convex where the other face leaves on this one's material side.
    let convex = dot(ways[1], outs[0]) < 0.0;
    let s = if convex { -radius } else { radius };
    let g = dot(outs[0], outs[1]);
    if (1.0 + g).abs() <= tol.angular() {
        return Ok(None);
    }
    let t = s / (1.0 + g);
    let ball = (
        r0 + (outs[0].0 + outs[1].0) * t,
        (outs[0].1 + outs[1].1) * t,
    );
    let mut feet = [(0.0, 0.0); 2];
    for k in 0..2 {
        let foot = (ball.0 - outs[k].0 * s, ball.1 - outs[k].1 * s);
        if dot((foot.0 - r0, foot.1), ways[k]) <= tol.confusion() || foot.0 <= tol.confusion() {
            return Ok(None);
        }
        feet[k] = foot;
    }
    if ball.0 <= radius + tol.confusion() {
        return Ok(None);
    }
    // Each face holds its contact ring: nothing else of a plane face stands
    // between the rim and the ring, and nothing else of a wall between the
    // rim and the ring's height; a straight edge leaving the rim's circle is
    // read at its far end. A face narrower than the round (a cap the ring
    // runs off) is left to the current path.
    let on_rim =
        |p: (f64, f64)| (p.0 - r0).abs() <= tol.confusion() && p.1.abs() <= tol.confusion();
    for (k, fi) in faces.iter().enumerate() {
        let plane = matches!(curved.faces[*fi].surface, SurfaceGeometry::Plane(_));
        let foot = feet[k];
        let between = |p: (f64, f64)| {
            let (a, b, x) = if plane {
                (r0, foot.0, p.0)
            } else {
                (0.0, foot.1, p.1)
            };
            x > a.min(b) + tol.confusion() && x < a.max(b) - tol.confusion()
        };
        for (e, _) in curved.faces[*fi].loops.iter().flatten() {
            if *e == edge {
                continue;
            }
            let c = &curved.edges[e];
            let ends = (
                meridian(curved.points[&c.ends.0]),
                meridian(curved.points[&c.ends.1]),
            );
            let samples: Vec<(f64, f64)> =
                if matches!(c.curve, Curve::Line(_)) && (on_rim(ends.0) || on_rim(ends.1)) {
                    vec![if on_rim(ends.0) { ends.1 } else { ends.0 }]
                } else {
                    let (t0, t1) = c.range;
                    (0..=RING_SAMPLES)
                        .map(|j| {
                            let t = t0 + (t1 - t0) * f64::from(j) / f64::from(RING_SAMPLES);
                            ogeom_geom::Curve3d::point_at(&c.curve, t, tol).map(meridian)
                        })
                        .collect::<OgeomResult<_>>()?
                };
            // One crossing the ring too is where the contact leaves its
            // face, the roll-over P2 labels; one standing wholly inside it
            // is under the round, which is not built here.
            let beyond = |p: &(f64, f64)| {
                let (a, b, x) = if plane {
                    (r0, foot.0, p.0)
                } else {
                    (0.0, foot.1, p.1)
                };
                (x - b) * (b - a) > tol.confusion() * (b - a).abs()
            };
            if samples.iter().copied().any(between) && !samples.iter().any(beyond) {
                return Ok(None);
            }
        }
    }
    let at = |h: f64| Frame::new(centre + up * h, frame.z(), x, tol);
    let ring = |(rho, h): (f64, f64)| -> OgeomResult<Curve> {
        Ok(CircleCurve::new(Circle::new(at(h)?, rho, tol)?).into())
    };
    let contacts = [ring(feet[0])?, ring(feet[1])?];
    let surface: SurfaceGeometry =
        TorusSurface::new(Torus::new(at(ball.1)?, ball.0, radius, tol)?).into();
    // The box the round can take up: between the rim and both feet along
    // the axis and off it.
    let span = |v: [f64; 3]| {
        (
            v.iter().copied().fold(f64::INFINITY, f64::min),
            v.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        )
    };
    let (h_lo, h_hi) = span([0.0, feet[0].1, feet[1].1]);
    let (r_lo, r_hi) = span([r0, feet[0].0, feet[1].0]);
    let off_axis = |radius: f64, inside: bool| ogeom_algo::HalfSpaces::Revolved {
        at: centre,
        axis: up,
        radius,
        slope: 0.0,
        inside,
    };
    let clear = ogeom_algo::HalfSpaces::Any(vec![
        ogeom_algo::HalfSpaces::Below {
            normal: -up,
            at: centre + up * h_hi,
        },
        ogeom_algo::HalfSpaces::Below {
            normal: up,
            at: centre + up * h_lo,
        },
        off_axis(r_lo, true),
        off_axis(r_hi, false),
    ]);
    Ok(Some((
        Strip {
            edge,
            faces,
            contacts,
            surface,
            section: Section::Arc { spine: ring(ball)? },
        },
        clear,
    )))
}

/// Whether every face of `curved` but the strip's two and those of the edges
/// the run rests on stands outside the box the run's piece can take up:
/// round its balls' arcs, the edges it rests on between the crossings and
/// the strip's own edge over the run's stations.
///
/// # Errors
///
/// As certifying the faces.
fn run_clear(
    model: &Model,
    curved: &Curved,
    strip: &Strip,
    run: &super::seats::EdgeRun,
    tol: Tolerances,
) -> OgeomResult<bool> {
    use ogeom_geom::Curve3d as _;
    let mut lo = Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut hold = |p: Point| {
        lo = Point::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
        hi = Point::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
    };
    let SurfaceGeometry::BSpline(surface) = &run.surface else {
        return Ok(false);
    };
    // A B-spline surface lies in its control net's hull.
    for w in surface.grid().points() {
        hold(w.point());
    }
    let edge = &curved.edges[&strip.edge];
    for (c, f) in run.centres.iter().zip(&run.feet) {
        let t = super::trimmed::param_on(&edge.curve, edge.range, *c, tol)?;
        hold(edge.curve.point_at(t, tol)?);
        hold(*f);
    }
    for p in &run.on_edge {
        hold(*p);
    }
    let pad = 2.0 * crate::marched::band_fit_target(tol);
    let below = |normal: Vector, at: Point| ogeom_algo::HalfSpaces::Below { normal, at };
    let clear = ogeom_algo::HalfSpaces::Any(vec![
        below(-Vector::X, Point::new(lo.x - pad, 0.0, 0.0)),
        below(Vector::X, Point::new(hi.x + pad, 0.0, 0.0)),
        below(-Vector::Y, Point::new(0.0, lo.y - pad, 0.0)),
        below(Vector::Y, Point::new(0.0, hi.y + pad, 0.0)),
        below(-Vector::Z, Point::new(0.0, 0.0, lo.z - pad)),
        below(Vector::Z, Point::new(0.0, 0.0, hi.z + pad)),
    ]);
    let mut own: Vec<usize> = strip.faces.to_vec();
    for e in &run.over {
        own.extend(faces_of(curved, *e));
    }
    let fit = crate::marched::band_fit_target(tol);
    for (fi, face) in curved.faces.iter().enumerate() {
        if own.contains(&fi) {
            continue;
        }
        if ogeom_algo::certify_face_within(model, &face.shape, &clear, fit, tol)?.is_err() {
            return Ok(false);
        }
    }
    Ok(true)
}
