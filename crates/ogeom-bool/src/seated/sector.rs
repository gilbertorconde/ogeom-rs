//! A piece of a rim's blend wedge, the turn of it a tangent chain's arc
//! rounds, read for the straight trim ([`super::trimmed_cut`]).
//!
//! A rim's wedge turned through less than a whole turn is five faces about
//! one axis: a band of a cylinder (the wall seat), a sector of an annulus of
//! a plane square to the axis (the cap seat), the quarter of a torus's tube
//! between them (the blend), and two planar ends, each through the axis.
//! Through the axis, each end is the square between the tube's centre, the
//! wall and the cap, less the quarter disc of the tube. The trim reads it as
//! it reads a straight band's wedge: each seat, the ends with them, lies on
//! a face of the solid facing the same way, bounded by a run of that face's
//! boundary and by the edges the seat shares with the blend; an end where a
//! band rounding the chain on past the arc has already been cut is the face
//! that cut left, covered whole.
//!
//! The conditions read here, where one fails the cut is not made by the
//! trim:
//!
//! - The blend is a torus at no placement; the wall a cylinder about its
//!   axis whose radius is the tube's distance from the axis and the tube's
//!   radius apart, facing away from the tube; the cap a plane square to the
//!   axis standing the tube's radius from the tube's level, facing away from
//!   the tube; and each end a plane through the axis.
//! - The edges each two faces share are one edge each: the rim (cap and
//!   wall), the cap's contact at the tube's distance (cap and blend), the
//!   wall's contact at the tube's level (wall and blend), and at each end a
//!   meridian of the tube (end and blend), a line across the cap and a line
//!   along the wall. Each face's loop is those edges and no other. The three
//!   arcs about the axis are where the square's corners put them and turn
//!   alike, at most half a turn, the middle of the rim standing inside both
//!   ends; each meridian is a quarter turn about the tube's centre whose
//!   middle stands inside the square.
//!
//! The wedge is then the sector of the ring its square sweeps between the
//! two ends: the part of the slab between the tube's level and the cap's
//! plane, between the contact's and the wall's radius, behind both ends'
//! planes, less the tube. Behind both planes is the turn between them, as it
//! is no more than half a turn. Every face of the solid other than the ones
//! the seats lie on standing out of that sector ([`Sector::stands_clear`]),
//! the sector's inside meets no face of the solid: the faces the seats lie
//! on lie on the planes and the cylinder that bound it. The sector is one
//! connected region, and beside each seat the solid's material is on its
//! side, so it lies inside the solid and the wedge with it. The solid less
//! the wedge is bounded by the solid's faces with the seats taken out, and
//! by the blend turned to face out of the material, as for a straight band.

use super::ring::{self, Axis};
use super::{Lies, Solid, Walked, curve_of};
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, Curve3d as _, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Model, Orientation, Shape, TShapeId};

/// A piece of a rim's wedge as the trim reads it: the region it lies in.
pub(super) struct Sector {
    /// The axis, from the tube's level toward the cap.
    axis: Axis,
    /// The tube's distance from the axis, where the cap's contact runs.
    tube: f64,
    /// The wall's radius.
    wall_radius: f64,
    /// The tube's radius, the cap's height over the tube's level.
    height: f64,
    /// Each end's plane: a point of it and its outward unit normal.
    ends: [(Point, Vector); 2],
    /// The widest tolerance the wedge states, with the confusion distance.
    reach: f64,
}

/// An arc about the axis: its radius, how far along the axis it stands and
/// the turn it makes, where its curve is a circle whose centre lies on the
/// axis and whose plane is square to it, both within `reach`.
fn arc_on(
    model: &Model,
    axis: &Axis,
    edge: &Shape,
    reach: f64,
    tol: Tolerances,
) -> Option<(f64, f64, f64)> {
    let Some((Curve::Circle(c), range, _)) = curve_of(model, edge) else {
        return None;
    };
    let circle = c.circle();
    if axis.off(circle.centre()) > reach
        || !ring::parallel(circle.frame().z().vector(), axis.up, tol)
    {
        return None;
    }
    Some((
        circle.radius(),
        axis.along(circle.centre()),
        (range.1 - range.0).abs(),
    ))
}

/// The one edge two loops share, `None` where they share none or more.
fn shared(x: &[Walked], y: &[Walked]) -> Option<Shape> {
    let mut out: Option<Shape> = None;
    for w in x {
        if y.iter().any(|v| v.edge.node() == w.edge.node()) {
            if out.as_ref().is_some_and(|e| e.node() != w.edge.node()) {
                return None;
            }
            out = Some(w.edge.oriented(Orientation::Forward));
        }
    }
    out
}

/// Whether a loop walks each of `edges` once and nothing else.
fn walks_only(ring: &[Walked], edges: &[&Shape]) -> bool {
    let ids: Vec<TShapeId> = edges.iter().map(|e| e.node()).collect();
    ring.len() == ids.len()
        && ids
            .iter()
            .all(|id| ring.iter().filter(|w| w.edge.node() == *id).count() == 1)
}

/// The wedge whose seats (cap, wall and ends) are `faces` and whose blend
/// is `blend` read as a [`Sector`]; `None` where it is not one as the
/// module asks.
pub(super) fn read_sector(
    model: &Model,
    faces: &[(Shape, Vec<Walked>, Lies)],
    blend: &Shape,
    blend_loop: &[Walked],
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Sector>> {
    let Some((_, SurfaceGeometry::Torus(t))) = ring::surface_of(model, blend) else {
        return Ok(None);
    };
    let torus = t.torus();
    let frame = torus.frame();
    let (tube, radius) = (torus.major_radius(), torus.minor_radius());
    let reach = tolerance + tol.confusion();
    let mut axis = Axis {
        origin: frame.origin(),
        up: frame.z().vector(),
    };
    // The seats by the surfaces they lie on.
    let mut wall: Option<(&[Walked], Point, Vector, f64, f64)> = None;
    let mut cap: Option<(&[Walked], Point, Vector)> = None;
    let mut ends: Vec<(&[Walked], Point, Vector)> = Vec::new();
    for (_, ring, lies) in faces {
        match *lies {
            Lies::Cylinder {
                origin,
                axis: along,
                radius,
                sense,
            } if wall.is_none() => wall = Some((ring, origin, along, radius, sense)),
            Lies::Plane { origin, outward } if ring::parallel(outward, axis.up, tol) => {
                if cap.is_some() {
                    return Ok(None);
                }
                cap = Some((ring, origin, outward));
            }
            Lies::Plane { origin, outward } if outward.dot(axis.up).abs() <= tol.angular() => {
                ends.push((ring, origin, outward));
            }
            _ => return Ok(None),
        }
    }
    let (Some(wall), Some(cap), [end0, end1]) = (wall, cap, ends.as_slice()) else {
        return Ok(None);
    };
    let (wall_loop, wall_origin, wall_axis, wall_radius, wall_sense) = wall;
    // The cap's side of the tube's level.
    let height = axis.along(cap.1);
    if height < 0.0 {
        axis.up = -axis.up;
    }
    let height = height.abs();
    let through = |&(_, origin, outward): &(&[Walked], Point, Vector)| {
        (axis.origin - origin).dot(outward).abs() <= reach
    };
    if tube <= reach
        || (height - radius).abs() > reach
        || ((wall_radius - tube).abs() - radius).abs() > reach
        || cap.2.dot(axis.up) <= 0.0
        || !ring::parallel(wall_axis, axis.up, tol)
        || axis.off(wall_origin) > reach
        || wall_sense != if wall_radius > tube { 1.0 } else { -1.0 }
        || !through(end0)
        || !through(end1)
    {
        return Ok(None);
    }
    // The edges each two faces share, and every loop those alone.
    let (Some(rim), Some(cap_ring), Some(wall_ring)) = (
        shared(cap.0, wall_loop),
        shared(cap.0, blend_loop),
        shared(wall_loop, blend_loop),
    ) else {
        return Ok(None);
    };
    let mut meridians = Vec::with_capacity(2);
    let mut across = Vec::with_capacity(4);
    for end in [end0, end1] {
        let (Some(meridian), Some(radial), Some(upright)) = (
            shared(end.0, blend_loop),
            shared(end.0, cap.0),
            shared(end.0, wall_loop),
        ) else {
            return Ok(None);
        };
        if !walks_only(end.0, &[&meridian, &radial, &upright]) {
            return Ok(None);
        }
        meridians.push(meridian);
        across.extend([radial, upright]);
    }
    let [m0, m1] = meridians.as_slice() else {
        return Ok(None);
    };
    let [r0, u0, r1, u1] = across.as_slice() else {
        return Ok(None);
    };
    if !walks_only(cap.0, &[&rim, &cap_ring, r0, r1])
        || !walks_only(wall_loop, &[&rim, &wall_ring, u0, u1])
        || !walks_only(blend_loop, &[&cap_ring, &wall_ring, m0, m1])
    {
        return Ok(None);
    }
    for line in &across {
        let Some((Curve::Line(_), _, _)) = curve_of(model, line) else {
            return Ok(None);
        };
    }
    // The arcs about the axis where the square's corners put them, each
    // turning alike, at most half a turn.
    let mut turn: Option<f64> = None;
    for (edge, at_radius, at_along) in [
        (&rim, wall_radius, height),
        (&cap_ring, tube, height),
        (&wall_ring, wall_radius, 0.0),
    ] {
        let Some((own, along, sweep)) = arc_on(model, &axis, edge, reach, tol) else {
            return Ok(None);
        };
        if (own - at_radius).abs() > reach
            || (along - at_along).abs() > reach
            || sweep <= tol.angular()
            || sweep > core::f64::consts::PI + tol.angular()
            || turn.is_some_and(|t| (t - sweep).abs() > tol.angular())
        {
            return Ok(None);
        }
        turn = Some(sweep);
    }
    // Each meridian a quarter turn about the tube's centre, its middle
    // standing inside the square.
    let (lo, hi) = (tube.min(wall_radius), tube.max(wall_radius));
    for meridian in [m0, m1] {
        let Some((curve @ Curve::Circle(c), range, _)) = curve_of(model, meridian) else {
            return Ok(None);
        };
        let circle = c.circle();
        let sweep = (range.1 - range.0).abs();
        let middle = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
        let (off, along) = (axis.off(middle), axis.along(middle));
        if (circle.radius() - radius).abs() > reach
            || axis.along(circle.centre()).abs() > reach
            || (axis.off(circle.centre()) - tube).abs() > reach
            || circle.frame().z().vector().dot(axis.up).abs() > tol.angular()
            || (sweep - core::f64::consts::FRAC_PI_2).abs() > tol.angular()
            || off <= lo + reach
            || off >= hi - reach
            || along <= reach
            || along >= height - reach
        {
            return Ok(None);
        }
    }
    // Every end of an arc behind both ends' planes, and the rim's middle
    // inside them: the turn between them is the rim's.
    let ends = [(end0.1, end0.2), (end1.1, end1.2)];
    for edge in [&rim, &cap_ring, &wall_ring] {
        let Some((curve, range, _)) = curve_of(model, edge) else {
            return Ok(None);
        };
        for t in [range.0, range.1] {
            let p = curve.point_at(t, tol)?;
            if ends.iter().any(|&(o, n)| (p - o).dot(n) > reach) {
                return Ok(None);
            }
        }
        if edge.node() == rim.node() {
            let p = curve.point_at(f64::midpoint(range.0, range.1), tol)?;
            if ends.iter().any(|&(o, n)| (p - o).dot(n) >= -reach) {
                return Ok(None);
            }
        }
    }
    Ok(Some(Sector {
        axis,
        tube,
        wall_radius,
        height,
        ends,
        reach,
    }))
}

impl Sector {
    /// Whether every face of the solid but `hosts` stands, by its
    /// conservative box and within the tolerance it states, out of the
    /// sector: past the cap's plane or the tube's level, past the wall's or
    /// the contact's radius, or past one end's plane.
    pub(super) fn stands_clear(
        &self,
        model: &Model,
        solid: &mut Solid,
        hosts: &[usize],
        tol: Tolerances,
    ) -> OgeomResult<bool> {
        let (lo, hi) = (
            self.tube.min(self.wall_radius),
            self.tube.max(self.wall_radius),
        );
        for i in 0..solid.faces.len() {
            if hosts.contains(&i) {
                continue;
            }
            let (bounds, stated) = solid.face_box(model, i, tol)?;
            let margin = stated + tol.confusion();
            let ((low, high), (least, most)) =
                ring::face_reach(model, &self.axis, &solid.faces[i], &bounds, self.reach, tol);
            if low >= self.height - margin
                || high <= margin
                || most <= lo + margin
                || least >= hi - margin
            {
                continue;
            }
            let corners = bounds.corners();
            let behind = self
                .ends
                .iter()
                .any(|&(o, n)| corners.iter().all(|c| (*c - o).dot(n) >= -margin));
            if !behind {
                return Ok(false);
            }
        }
        Ok(true)
    }
}
