//! Chamfers of solids with curved faces: P1's contacts on curved seats, P3's
//! vertex cases there, and P4's trimmed-face assembly, over the solid
//! [`curved`](super::curved) reads and the corners curved P2
//! ([`curved_corners`]) finds.
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

use super::curved::{Curved, CurvedEdge, SheetCurve, SheetRings, read};
use super::setbacks::{Corner, corner_for, curved_corners, pinched};
use super::{Id, curved, fallback, id, seats, setbacks};
use crate::chamfer::{Chamfer, Legs, legs, rim_bevel_surface, wall_depth};
use crate::chamfer_sphere::{SphereSeat, sphere_bevel_rings};
use crate::support::RevolvedSeat;

mod assemble;
mod bevels;
mod cusps;
mod one_boolean;
mod regions;
mod rolling_tools;
#[cfg(test)]
mod tests;

use assemble::assemble_widened;
pub(crate) use assemble::{RolledRun, assemble_rolled};
use bevels::{bevel_of, distances, general_crease, measured_within};
use cusps::{
    HoleSpec, MERGE_CAP, Met, PatchRing, RIDGE_SAMPLES, arc_control, cusps_of, met_between,
    periodic_through, piece_of, rational_quadratics, segment_on, shared_point, skew_control,
    torus_skew,
};
use one_boolean::into_face;
pub(crate) use one_boolean::{convex_one_boolean, creases_one_boolean, rims_revolved};
use regions::{
    Fate, Key, Part, Points, ROLLED_EDGE_WIDTH, Span, chained_ring, grouped, period_of, piece,
    rejoined, resting_region, runs, span_between, wrapped,
};
pub(crate) use rolling_tools::{EndCap, close_tool, open_tool_over, rolling_tools};
use rolling_tools::{RimSection, asked_apart, pieces_close, revolved_tool};

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
    declining(chamfer_with(model, solid, specs, false, false, tol))
}

/// [`chamfer`] with each vertex of three bevels or more mitred: the bevels
/// run on to the point their surfaces share and meet each other along
/// their sections from there, with no facet.
///
/// # Errors
///
/// As [`chamfer`].
pub(crate) fn chamfer_mitred(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    declining(chamfer_with(model, solid, specs, false, true, tol))
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
    declining(chamfer_with(model, solid, specs, true, false, tol))
}

fn chamfer_with(
    model: &mut Model,
    solid: &Shape,
    specs: &[(Shape, Chamfer)],
    patched_only: bool,
    mitre: bool,
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
    // An edge not of the solid is no request of the curved path's: handed
    // back, refused by name.
    if asked.iter().any(|(e, _, _)| !curved.edges.contains_key(e)) {
        return Ok(None);
    }
    // In order of where each edge stands, its key only between two that
    // stand alike: a dissolved solid's keys change from one reading to the
    // next, and the bevels' order is the result's face order.
    let place = |e: Id| -> [u64; 3] {
        let c = &curved.edges[&e];
        c.curve
            .point_at(0.5 * (c.range.0 + c.range.1), tol)
            .map_or([u64::MAX; 3], |p| {
                [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]
            })
    };
    // The pieces of one curve the boolean left between the same two faces,
    // meeting at vertices of no other edge, are one edge of the form: the
    // chamfer asked on one runs on into the rest, one bevel round them. A
    // piece asked otherwise is handed back.
    let Some(chains) = chain_pieces(&curved, &mut asked) else {
        return Ok(None);
    };
    asked.sort_by_cached_key(|(e, _, _)| (place(*e), *e));
    let mut bevels: Vec<Strip> = Vec::with_capacity(asked.len());
    // Each fitted contact's distance off its face, which its edges are
    // widened to.
    let mut widths: FastMap<(usize, Id), f64> = FastMap::default();
    for (edge, shape, spec) in &asked {
        if let Some(bevel) = bevel_of(model, solid, &curved, *edge, shape, spec, false, tol)? {
            bevels.push(bevel);
            continue;
        }
        let [f0, f1] = faces_of(&curved, *edge)[..] else {
            return Ok(None);
        };
        // A piece of a closed run: the bevel once round the whole loop, from
        // the run's first vertex, its contacts and surface each piece's.
        let run = chains
            .iter()
            .find(|(r, closed)| *closed && r.contains(edge));
        // From the vertex where the crease's own curve starts (its fitted
        // seam's join), as a marched band runs.
        let whole = run.map(|(r, _)| {
            let start = r
                .iter()
                .find(|x| {
                    let c = &curved.edges[*x];
                    (c.range.0 - c.curve.domain().0).abs() <= tol.parametric()
                })
                .unwrap_or(&r[0]);
            curved.points[&curved.edges[start].ends.0]
        });
        let shared = run.and_then(|(r, _)| {
            bevels
                .iter()
                .find(|b| r.contains(&b.edge) && b.edge != *edge)
                .map(|b| {
                    (
                        b.clone(),
                        widths.get(&(b.faces[0], b.edge)).copied(),
                        widths.get(&(b.faces[1], b.edge)).copied(),
                    )
                })
        });
        let (bevel, off) = match shared {
            Some((mut b, Some(w0), Some(w1))) if b.faces == [f0, f1] => {
                b.edge = *edge;
                (b, [w0, w1])
            }
            _ => {
                let Some(found) = general_crease(
                    model,
                    solid,
                    &curved,
                    *edge,
                    shape,
                    [f0, f1],
                    spec,
                    whole,
                    tol,
                )?
                else {
                    return Ok(None);
                };
                found
            }
        };
        for (face, width) in bevel.faces.iter().zip(off) {
            widths.insert((*face, bevel.edge), width);
        }
        bevels.push(bevel);
    }
    // The pieces of one run whose bevels are one surface (each built on its
    // own, a cone's frame may differ) take the first piece's, so their
    // strips are one face.
    for (run, _) in &chains {
        // The piece starting lowest on the run's curve, not the first in
        // key order: the keys of a dissolved solid change from one reading
        // to the next.
        let Some(first) = (0..bevels.len())
            .filter(|bi| run.contains(&bevels[*bi].edge))
            .min_by(|x, y| {
                let start = |bi: &usize| curved.edges[&bevels[*bi].edge].range.0;
                start(x).total_cmp(&start(y))
            })
        else {
            continue;
        };
        let surface = bevels[first].surface.clone();
        for (bi, bevel) in bevels.iter_mut().enumerate() {
            if bi == first || !run.contains(&bevel.edge) || bevel.surface == surface {
                continue;
            }
            let e = &curved.edges[&bevel.edge];
            let mid = e.curve.point_at(0.5 * (e.range.0 + e.range.1), tol)?;
            let mut on = true;
            for c in &bevel.contacts {
                let (lo, hi) = c.domain();
                let p = c.point_at(param_on(c, (lo.max(-1e6), hi.min(1e6)), mid, tol)?, tol)?;
                let foot = ogeom_geom::project::project_on_surface(&surface, p, 16, tol)?;
                on &= foot.distance <= tol.confusion();
            }
            if on {
                bevel.surface = surface.clone();
            }
        }
    }
    // A wall the bevels of its two rims consume: its strips meet along
    // their ridge instead.
    let Some(bands) = consumed_bands(&curved, &mut bevels, tol)? else {
        return Ok(None);
    };
    widths.extend(bands.widths.iter().map(|(k, w)| (*k, *w)));
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
    // Each contact as far off its corners as it stands: two fitted pieces
    // of one run meet at their joint each within its own fit of the
    // corner, and the edges there widen to cover the gap.
    for ((fi, _), corner) in &corners {
        for (e, t) in &corner.on {
            let Some(b) = bevels.iter().find(|b| b.edge == *e) else {
                continue;
            };
            let Some(k) = b.faces.iter().position(|f| f == fi) else {
                continue;
            };
            let off = b.contacts[k].point_at(*t, tol)?.distance(corner.point);
            if off > tol.confusion() {
                let width = widths.entry((*fi, *e)).or_insert(0.0);
                *width = width.max(off + tol.confusion());
            }
        }
    }
    let blended: FastSet<Id> = bevels.iter().map(|b| b.edge).collect();
    if patched_only && !off_plane_corner(&curved, &blended, &corners, tol) {
        return Ok(None);
    }
    // P4 and P5. Whatever fails here, its assembly or the solid's check,
    // is outside the cases built, and the caller keeps its own path.
    let built = match assemble_widened(
        model,
        &curved,
        &bevels,
        &blended,
        &corners,
        &bands.faces,
        &widths,
        mitre,
        &[],
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

/// The runs of pieces of one curve between the same two faces, each
/// continuing the last through a vertex of no other edge, that hold an
/// asked edge: every piece of such a run is asked as that edge is. Back
/// come the runs, each in order along the curve with whether it closes;
/// `None` where a piece of a run was asked with another chamfer.
#[allow(
    clippy::type_complexity,
    reason = "the asked list as the caller holds it"
)]
fn chain_pieces<'a>(
    curved: &'a Curved,
    asked: &mut Vec<(Id, &'a Shape, &'a Chamfer)>,
) -> Option<Vec<(Vec<Id>, bool)>> {
    let mut runs: Vec<(Vec<Id>, bool)> = Vec::new();
    let mut k = 0;
    while k < asked.len() {
        let (e, _, spec) = asked[k];
        k += 1;
        if runs.iter().any(|(r, _)| r.contains(&e)) {
            continue;
        }
        let c = &curved.edges[&e];
        if c.ends.0 == c.ends.1 {
            continue;
        }
        let mut faces = faces_of(curved, e);
        faces.sort_unstable();
        // The piece continuing `x` through its end (or start).
        let next = |x: Id, forward: bool| -> Option<Id> {
            let c = &curved.edges[&x];
            let v = if forward { c.ends.1 } else { c.ends.0 };
            let at: Vec<Id> = curved
                .edges
                .iter()
                .filter(|(y, o)| **y != x && (o.ends.0 == v || o.ends.1 == v))
                .map(|(y, _)| *y)
                .collect();
            let [y] = at[..] else {
                return None;
            };
            let o = &curved.edges[&y];
            let mut beside = faces_of(curved, y);
            beside.sort_unstable();
            let joins = if forward {
                o.ends.0 == v
            } else {
                o.ends.1 == v
            };
            (o.curve == c.curve && beside == faces && joins).then_some(y)
        };
        let mut run = vec![e];
        let mut closed = false;
        while let Some(y) = next(*run.last()?, true) {
            if y == run[0] {
                closed = true;
                break;
            }
            if run.contains(&y) {
                break;
            }
            run.push(y);
        }
        if !closed {
            while let Some(y) = next(run[0], false) {
                if run.contains(&y) {
                    break;
                }
                run.insert(0, y);
            }
        }
        if run.len() < 2 {
            continue;
        }
        // It closes where its last piece ends at its first's start, other
        // edges meeting there or not.
        let closed =
            closed || curved.edges[&run[0]].ends.0 == curved.edges[&run[run.len() - 1]].ends.1;
        for piece in &run {
            match asked.iter().find(|(x, _, _)| x == piece) {
                Some((_, _, theirs)) if !super::seams::alike(theirs, spec) => return None,
                Some(_) => {}
                None => asked.push((*piece, &curved.edges[piece].shape, spec)),
            }
        }
        runs.push((run, closed));
    }
    Some(runs)
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
    /// The width a fitted ridge's edges are built to, by the consumed face
    /// and each strip's edge.
    widths: FastMap<(usize, Id), f64>,
}

/// The walls whose two rims' bevels cross on them (a disc too thin for
/// both rims' chamfers): a cylinder, a cone, a sphere's zone or a torus's
/// band whose bevelled edges lie on two coaxial circles, the contact of the
/// higher at or below the lower's.
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
            SurfaceGeometry::Torus(t) => {
                let axis = t.torus().axis();
                (axis.location, axis.direction.vector())
            }
            // A sphere's axis is its rims': through its centre, square to
            // the first rim on it.
            SurfaceGeometry::Sphere(sphere) => {
                let rim = strips.iter().find_map(|st| {
                    if !st.faces.contains(&fi) {
                        return None;
                    }
                    match &curved.edges[&st.edge].curve {
                        Curve::Circle(c) => Some(c.circle().frame().z().vector()),
                        _ => None,
                    }
                });
                let Some(a) = rim else {
                    continue;
                };
                (sphere.sphere().centre(), a)
            }
            _ => continue,
        };
        let height = |p: Point| (p - o).dot(a);
        // A circle about the wall's axis: centred on it, square to it.
        let coaxial = |c: &Circle| {
            let off = c.centre() - o;
            (off - a * off.dot(a)).magnitude() <= tol.confusion()
                && c.frame().z().vector().cross(a).magnitude() <= tol.angular()
        };
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
            if !coaxial(&rim.circle()) || !coaxial(&contact.circle()) {
                continue;
            }
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
    // A face two bevels' contacts cross on, any surface, the two edges
    // across it from each other (a web too narrow for its two uprights'
    // chamfers): the bevels meet in a ridge, each face's vertex standing on
    // it at its foot.
    for (fi, face) in curved.faces.iter().enumerate() {
        if bands.faces.contains(&fi) {
            continue;
        }
        let on: Vec<usize> = strips
            .iter()
            .enumerate()
            .filter(|(_, st)| st.faces.contains(&fi))
            .map(|(si, _)| si)
            .collect();
        // The strips in two groups, each one bevel's surface (the pieces of
        // one rim the boolean left share their cone).
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for si in on {
            match groups
                .iter_mut()
                .find(|g| strips[g[0]].surface == strips[si].surface)
            {
                Some(g) => g.push(si),
                None => groups.push(vec![si]),
            }
        }
        let [group_a, group_b] = &groups[..] else {
            continue;
        };
        let (sa, sb) = (group_a[0], group_b[0]);
        // Two edges across the face from each other, not meeting at one
        // of its corners.
        let ends = |g: &[usize]| -> Vec<Id> {
            g.iter()
                .flat_map(|si| {
                    let c = &curved.edges[&strips[*si].edge];
                    [c.ends.0, c.ends.1]
                })
                .collect()
        };
        let (ends_a, ends_b) = (ends(group_a), ends(group_b));
        if ends_a.iter().any(|v| ends_b.contains(v)) {
            continue;
        }
        // How far a point stands past a bevel, out of the material its two
        // faces bound: along its surface's normal at the foot of the point.
        let past = |st: &Strip, p: Point| -> OgeomResult<f64> {
            let e = &curved.edges[&st.edge];
            let mid = e.curve.point_at(0.5 * (e.range.0 + e.range.1), tol)?;
            let out = outward_at(curved, st.faces[0], mid) + outward_at(curved, st.faces[1], mid);
            let foot = ogeom_geom::project::project_on_surface(&st.surface, p, 16, tol)?;
            let n = ogeom_geom::Surface::normal_at(
                &st.surface,
                foot.parameters.0,
                foot.parameters.1,
                tol,
            )?
            .vector();
            let n = if n.dot(out) < 0.0 { -n } else { n };
            Ok((p - foot.point).dot(n))
        };
        // Crossed: each contact's middle (beside its edge's own, a rim's
        // piece reading its whole circle) stands past the other's bevel.
        let middle = |st: &Strip| -> OgeomResult<Point> {
            let k = usize::from(st.faces[1] == fi);
            let c = &st.contacts[k];
            let (lo, hi) = c.domain();
            let e = &curved.edges[&st.edge];
            let mid = e.curve.point_at(0.5 * (e.range.0 + e.range.1), tol)?;
            c.point_at(param_on(c, (lo.max(-1e6), hi.min(1e6)), mid, tol)?, tol)
        };
        let (ma, mb) = (middle(&strips[sa])?, middle(&strips[sb])?);
        if past(&strips[sb], ma)? <= tol.confusion() || past(&strips[sa], mb)? <= tol.confusion() {
            continue;
        }
        // And near: each middle within two of the other bevel's widths of
        // one of its edges, not past its surface where that runs on far
        // from them.
        let near = |st: &Strip, p: Point| -> OgeomResult<bool> {
            let e = &curved.edges[&st.edge];
            let t = param_on(&e.curve, e.range, p, tol)?;
            let reach = st
                .contacts
                .iter()
                .map(|c| -> OgeomResult<f64> {
                    let (lo, hi) = c.domain();
                    let (lo, hi) = (lo.max(-1e6), hi.min(1e6));
                    let q = e.curve.point_at(t, tol)?;
                    Ok(c.point_at(param_on(c, (lo, hi), q, tol)?, tol)?.distance(q))
                })
                .collect::<OgeomResult<Vec<f64>>>()?
                .into_iter()
                .fold(0.0, f64::max);
            Ok(e.curve.point_at(t, tol)?.distance(p) <= 2.0 * reach)
        };
        let near_any = |group: &[usize], p: Point| -> OgeomResult<bool> {
            for si in group {
                if near(&strips[*si], p)? {
                    return Ok(true);
                }
            }
            Ok(false)
        };
        if !near_any(group_b, ma)? || !near_any(group_a, mb)? {
            continue;
        }
        let vertices: Vec<Id> = face
            .loops
            .iter()
            .flatten()
            .flat_map(|(e, _)| {
                let c = &curved.edges[e];
                [c.ends.0, c.ends.1]
            })
            .collect();
        // The ridge: two planes meet in a line, over the face's vertices'
        // feet on it; any other two bevels in their section, marched, the
        // curve passing nearest both contacts, its edges widened to its
        // tolerance.
        let (ridge, width): (Curve, f64) = match (&strips[sa].surface, &strips[sb].surface) {
            (SurfaceGeometry::Plane(a), SurfaceGeometry::Plane(b)) => {
                let (na, pa) = (a.plane().normal().vector(), a.plane().frame().origin());
                let (nb, pb) = (b.plane().normal().vector(), b.plane().frame().origin());
                let along = na.cross(nb);
                if along.magnitude() <= tol.angular() {
                    return Ok(None);
                }
                let along = along / along.magnitude();
                // A point on both planes: from `pa`, along the first plane
                // square to the line, to the second.
                let across = along.cross(na);
                let step = (pb - pa).dot(nb) / across.dot(nb);
                let through = pa + across * step;
                let at = |v: &Id| (curved.points[v] - through).dot(along);
                let (lo, hi) = vertices
                    .iter()
                    .map(at)
                    .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), t| {
                        (l.min(t), h.max(t))
                    });
                if hi - lo <= tol.confusion() {
                    return Ok(None);
                }
                (
                    LineCurve::segment(through + along * lo, through + along * hi, tol)?.into(),
                    0.0,
                )
            }
            (a, b) => {
                use ogeom_intersect::section::{SurfaceIntersection, intersect_surfaces};
                let Ok(SurfaceIntersection::Along(curves)) =
                    intersect_surfaces(a, b, section_options(), tol)
                else {
                    return Ok(None);
                };
                let between = ma.midpoint(mb);
                let mut best: Option<(f64, ogeom_intersect::section::SectionCurve)> = None;
                for section in curves {
                    let dom = section.curve.domain();
                    let t = param_on(&section.curve, dom, between, tol)?;
                    let d = section.curve.point_at(t, tol)?.distance(between);
                    if best.as_ref().is_none_or(|b| d < b.0) {
                        best = Some((d, section));
                    }
                }
                let Some((_, section)) = best else {
                    return Ok(None);
                };
                // A closed section fitted clamped is made periodic, so its
                // pieces may run round its join.
                if section.closed && !section.curve.is_periodic() {
                    let (curve, off) = periodic_through(&section.curve, RIDGE_SAMPLES, tol)?;
                    (curve.into(), section.tolerance + off + tol.confusion())
                } else {
                    (section.curve, section.tolerance + tol.confusion())
                }
            }
        };
        for si in group_a.iter().chain(group_b) {
            let st = &mut strips[*si];
            let k = usize::from(st.faces[1] == fi);
            st.contacts[k] = ridge.clone();
            if width > tol.confusion() {
                bands.widths.insert((fi, st.edge), width);
            }
        }
        let dom = ridge.domain();
        for v in &vertices {
            let at = curved.points[v];
            let mut t = param_on(&ridge, dom, at, tol)?;
            // On a cone or cylinder bevel the strip's end at the vertex is
            // its generator there: the corner is where the ridge crosses the
            // plane through its axis and the vertex.
            let axis = group_a
                .iter()
                .chain(group_b)
                .filter(|si| {
                    let c = &curved.edges[&strips[**si].edge];
                    c.ends.0 == *v || c.ends.1 == *v
                })
                .find_map(|si| match &strips[*si].surface {
                    SurfaceGeometry::Cone(c) => {
                        let f = c.cone().frame();
                        Some((f.origin(), f.z().vector()))
                    }
                    SurfaceGeometry::Cylinder(c) => {
                        let a = c.cylinder().axis();
                        Some((a.location, a.direction.vector()))
                    }
                    _ => None,
                });
            if let Some((o, a)) = axis {
                let n = a.cross(at - o);
                if n.magnitude() > tol.confusion() {
                    let n = n / n.magnitude();
                    for _ in 0..50 {
                        let g = (ridge.point_at(t, tol)? - at).dot(n);
                        let dg = ridge.d1_at(t, tol)?.dot(n);
                        if g.abs() <= 1e-3 * tol.confusion() || dg.abs() <= f64::EPSILON {
                            break;
                        }
                        t -= g / dg;
                    }
                }
            }
            bands.corners.insert((fi, *v), ridge.point_at(t, tol)?);
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
pub(crate) fn faces_of(curved: &Curved, edge: Id) -> Vec<usize> {
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

/// A face's outward unit normal at `p`: in closed form on a plane, a
/// cylinder, a sphere or a torus, else the surface's own normal at the
/// foot of `p`; `None` where it has none there.
pub(crate) fn outward(curved: &Curved, face: usize, p: Point) -> Option<Vector> {
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
        // Any other surface: its own normal at the foot of `p`.
        surface => {
            let foot =
                ogeom_geom::project::project_on_surface(surface, p, 16, Tolerances::default())
                    .ok()?;
            let (u, v) = foot.parameters;
            let (du, dv) = ogeom_geom::Surface::d1_at(surface, u, v, Tolerances::default()).ok()?;
            let n = du.cross(dv);
            let m = n.magnitude();
            (m > 0.0).then(|| n / m * sign)
        }
    }
}

/// A point of face `face` away from edge `edge`: a vertex of its outline
/// not on the edge, or where it has none, the middle of another of its
/// edges. `None` where it has neither.
pub(crate) fn away(curved: &Curved, face: usize, edge: Id) -> Option<Point> {
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

/// How a section with no closed form is marched: traced to ten confusion
/// distances and fitted to one. The edge is built to what that leaves, and
/// the faces either side must meet it.
fn section_options() -> ogeom_intersect::section::IntersectOptions {
    ogeom_intersect::section::IntersectOptions {
        tolerance: 1e-7,
        marching: ogeom_intersect::march::Marching {
            chord: 1e-6,
            ..Default::default()
        },
    }
}

/// Whether two curved strips ending on `p` and `q` mirror each other
/// through `plane`: the curve `plane` cuts from `surface` between the two
/// lies on `other` too, to the confusion distance at every sample.
fn mirrored(
    surface: &SurfaceGeometry,
    other: &SurfaceGeometry,
    plane: &SurfaceGeometry,
    p: Point,
    q: Point,
    tol: Tolerances,
) -> OgeomResult<bool> {
    const SAMPLES: u32 = 8;
    let Some((curve, (lo, hi), _)) = met_between(surface, plane, p, q, tol)? else {
        return Ok(false);
    };
    for i in 1..SAMPLES {
        let t = lo + (hi - lo) * f64::from(i) / f64::from(SAMPLES);
        let at = curve.point_at(t, tol)?;
        let foot = ogeom_geom::project::project_on_surface(other, at, 16, tol)?;
        if foot.point.distance(at) > tol.confusion() {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Where two curved surfaces with no closed-form meeting meet between `p`
/// and `q`, both on each: the marched section's curve through the two,
/// its range between them (the low end first, and whether that is `p`),
/// and the width its edge is built to: the section's own tolerance plus
/// how far the curve passes from either end. `None` where no section
/// curve passes through both, or the piece between them is not one run of
/// its parameter.
fn marched_between(
    surface: &SurfaceGeometry,
    cut: &SurfaceGeometry,
    p: Point,
    q: Point,
    tol: Tolerances,
) -> OgeomResult<Option<(Met, f64)>> {
    use ogeom_intersect::section::{SurfaceIntersection, intersect_surfaces};
    let Ok(SurfaceIntersection::Along(curves)) =
        intersect_surfaces(surface, cut, section_options(), tol)
    else {
        return Ok(None);
    };
    let middle = p.midpoint(q);
    for section in curves {
        let curve = section.curve;
        let domain = curve.domain();
        let (tp, tq) = (
            param_on(&curve, domain, p, tol)?,
            param_on(&curve, domain, q, tol)?,
        );
        let off = curve
            .point_at(tp, tol)?
            .distance(p)
            .max(curve.point_at(tq, tol)?.distance(q));
        if off > section.tolerance + tol.confusion() {
            continue;
        }
        let (lo, hi, p_low) = if tp <= tq {
            (tp, tq, true)
        } else {
            (tq, tp, false)
        };
        if section.closed {
            // A loop holds two pieces between the ends: the one taken must
            // be the nearer to the ends' midpoint of the two.
            let inner = curve.point_at(0.5 * (lo + hi), tol)?.distance(middle);
            let (a, b) = (
                curve.point_at(0.5 * (domain.0 + lo), tol)?.distance(middle),
                curve.point_at(0.5 * (hi + domain.1), tol)?.distance(middle),
            );
            if a.min(b) < inner {
                continue;
            }
        }
        return Ok(Some(((curve, (lo, hi), p_low), section.tolerance + off)));
    }
    Ok(None)
}

/// The points outside `region`: each plane and wall turned to its other
/// side, the greatest and the least of pieces exchanged.
fn outside(region: &ogeom_algo::HalfSpaces) -> ogeom_algo::HalfSpaces {
    use ogeom_algo::HalfSpaces;
    match region {
        HalfSpaces::Below { normal, at } => HalfSpaces::Below {
            normal: -*normal,
            at: *at,
        },
        HalfSpaces::Revolved {
            at,
            axis,
            radius,
            slope,
            inside,
        } => HalfSpaces::Revolved {
            at: *at,
            axis: *axis,
            radius: *radius,
            slope: *slope,
            inside: !*inside,
        },
        HalfSpaces::Ball {
            centre,
            radius,
            inside,
        } => HalfSpaces::Ball {
            centre: *centre,
            radius: *radius,
            inside: !*inside,
        },
        HalfSpaces::Tube {
            centre,
            axis,
            major,
            minor,
            inside,
        } => HalfSpaces::Tube {
            centre: *centre,
            axis: *axis,
            major: *major,
            minor: *minor,
            inside: !*inside,
        },
        HalfSpaces::Every(pieces) => HalfSpaces::Any(pieces.iter().map(outside).collect()),
        HalfSpaces::Any(pieces) => HalfSpaces::Every(pieces.iter().map(outside).collect()),
    }
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
            // A fitted strip (no closed form) ends on the plane wherever
            // its corners stand on it: the strip cut by it.
            Section::Arc {
                spine: Curve::BSpline(_),
            } => true,
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

/// A strip surface's own normal at `p`: a plane's, a cone's or a
/// cylinder's off its axis there, a torus's off its tube's core, or a
/// B-spline's at the foot of `p`.
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
        SurfaceGeometry::Torus(t) => {
            let torus = t.torus();
            let frame = torus.frame();
            let (o, a) = (frame.origin(), frame.z().vector());
            let radial = (p - o) - a * (p - o).dot(a);
            let core = o + radial / radial.magnitude() * torus.major_radius();
            p - core
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
