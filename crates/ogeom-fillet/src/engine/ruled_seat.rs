//! The ruled fillet: a straight crease whose faces are planes or drums
//! running parallel to it, blended exactly.
//!
//! A wall meeting a drum along one of its rulings (the vertical edge where
//! an extruded profile's line meets its arc at an angle) or two drums
//! meeting along a shared ruling: every section square to the edge is the
//! same, so the blend is a section swept straight. The section is a plane
//! problem: each host cuts the section plane in a line or a circle, the
//! ball in a circle tangent to both, and the wedge is the region between
//! the crease point, the two touch points and the ball's arc. Extruded
//! along the edge, its faces are planes and drums: the legs lie exactly on
//! the hosts and the band is a drum of the fillet's radius, so the melt is
//! exact and a band's cap meets a corner ball's rim on the same circle.

use crate::march::{Sides, seat_section};
use crate::support::{edge_curve, same_occurrence};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _, SurfaceGeometry};
use ogeom_math::{Point, Vector};
use ogeom_topo::{Filter, Model, Orientation, Shape, ShapeType, explore};

/// A host's trace in the section plane: a line or a circle about a point.
pub(crate) enum Trace {
    Line,
    Circle { centre_on_axis: Point, radius: f64 },
}

/// A straight crease between two planes or drums parallel to it, at least
/// one a drum: the seat the ruled fillet blends.
pub(crate) struct RuledSeat {
    pub(crate) guide: Curve,
    pub(crate) range: (f64, f64),
    pub(crate) start: Point,
    pub(crate) end: Point,
    pub(crate) length: f64,
    pub(crate) along: Vector,
    pub(crate) faces: Vec<Shape>,
    pub(crate) hosts: Vec<(SurfaceGeometry, f64, Trace)>,
}

/// The ruled seat of a straight `edge` of `solid`; `None` where the seat is
/// not of that kind.
///
/// # Errors
///
/// A construction error where the edge has no length.
pub(crate) fn ruled_seat(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<RuledSeat>> {
    let (guide, range) = edge_curve(model, edge, tol)?;
    let Curve::Line(_) = &guide else {
        return Ok(None);
    };
    let start = guide.point_at(range.0, tol)?;
    let end = guide.point_at(range.1, tol)?;
    let length = start.distance(end);
    if length <= tol.confusion() {
        ogeom_bail!(Construction, "the edge has no length to blend along");
    }
    let along = (end - start) / length;

    let mut faces: Vec<Shape> = Vec::new();
    let mut hosts: Vec<(SurfaceGeometry, f64, Trace)> = Vec::new();
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let touches = explore(model, &face, Filter::OfType(ShapeType::Edge))?
            .iter()
            .any(|e| same_occurrence(model, e, edge, tol));
        if !touches {
            continue;
        }
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face().cloned()) else {
            ogeom_bail!(Construction, "a face node holds no face data");
        };
        let Some(stored) = model.geometry().surface(data.surface).cloned() else {
            ogeom_bail!(Dangling, "a face's surface is not in this model");
        };
        use ogeom_geom::Transformable as _;
        let surface = stored.transformed(&face.transform(model.datums())?, tol)?;
        let trace = match &surface {
            SurfaceGeometry::Plane(_) => Trace::Line,
            SurfaceGeometry::Cylinder(c) => {
                let cylinder = c.cylinder();
                let axis = cylinder.frame().z().vector();
                if axis.cross(along).magnitude() > tol.angular() {
                    return Ok(None);
                }
                // The axis's crossing with the section plane at the start.
                let o = cylinder.frame().origin();
                Trace::Circle {
                    centre_on_axis: o + axis * (start - o).dot(axis),
                    radius: cylinder.radius(),
                }
            }
            _ => return Ok(None),
        };
        let sign = if face.orientation() == Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        faces.push(face);
        hosts.push((surface, sign, trace));
    }
    if hosts.len() != 2 {
        return Ok(None);
    }
    if hosts.iter().all(|(_, _, t)| matches!(t, Trace::Line)) {
        return Ok(None);
    }
    Ok(Some(RuledSeat {
        guide,
        range,
        start,
        end,
        length,
        along,
        faces,
        hosts,
    }))
}

/// Exact ball contacts for the same ruled seat used by the builder. The
/// tangent-plane setback is not a width on a cylindrical wall: solve the
/// two offset supports, then measure the actual contact's reach instead.
pub(crate) fn contacts(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Vec<crate::support::Contact>>> {
    use ogeom_geom::Surface as _;
    let Some(seat) = ruled_seat(model, solid, edge, tol)? else {
        return Ok(None);
    };
    let convex = crate::marched::crease_convexity(
        model,
        &seat.faces[0],
        [
            (&seat.hosts[0].0, seat.hosts[0].1),
            (&seat.hosts[1].0, seat.hosts[1].1),
        ],
        &seat.guide,
        seat.range,
        radius,
        tol,
    )?;
    let sign = if convex { -1.0 } else { 1.0 };
    let sides = Sides {
        first: if sign * seat.hosts[0].1 > 0.0 { 1 } else { -1 },
        second: if sign * seat.hosts[1].1 > 0.0 { 1 } else { -1 },
    };
    let mut contacts = Vec::new();
    for fraction in [0.25, 0.5, 0.75] {
        let t = seat.range.0 + (seat.range.1 - seat.range.0) * fraction;
        let p = seat.guide.point_at(t, tol)?;
        let a = ogeom_algo::project_on_surface(&seat.hosts[0].0, p, 32, tol)?.parameters;
        let b = ogeom_algo::project_on_surface(&seat.hosts[1].0, p, 32, tol)?.parameters;
        let x = seat_section(
            &seat.hosts[0].0,
            &seat.hosts[1].0,
            radius,
            &seat.guide,
            sides,
            t,
            [a.0, a.1, b.0, b.1],
            tol,
        )?;
        for (i, (u, v)) in [(x[0], x[1]), (x[2], x[3])].into_iter().enumerate() {
            let at = seat.hosts[i].0.point_at(u, v, tol)?;
            if !crate::support::chart_holds(
                model,
                &seat.faces[i],
                ogeom_math::Point2::new(u, v),
                tol,
            )? {
                ogeom_bail!(
                    Construction,
                    "a ruled blend of radius {radius} touches outside its face at {at:?}"
                );
            }
            contacts.push(crate::support::Contact {
                face: seat.faces[i].clone(),
                at,
                setback: crate::support::distance_to_edge(&seat.guide, seat.range, at, tol)?,
            });
        }
    }
    Ok(Some(contacts))
}
