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
use ogeom_topo::{Filter, Model, Orientation, Shape, ShapeType, explore};

/// A host's trace in the section plane: a line or a circle about a point.
pub(crate) enum Trace {
    Line,
    Circle,
}

/// A straight crease between two planes or drums parallel to it, at least
/// one a drum: the seat the ruled fillet blends.
pub(crate) struct RuledSeat {
    pub(crate) guide: Curve,
    pub(crate) range: (f64, f64),
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
                Trace::Circle
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
            // A drum's angle is read a turn either way: its face's chart
            // may run across the drum's seam.
            let turns: &[f64] = match seat.hosts[i].2 {
                Trace::Circle => &[0.0, core::f64::consts::TAU, -core::f64::consts::TAU],
                Trace::Line => &[0.0],
            };
            let mut held = false;
            for turn in turns {
                held = crate::support::chart_holds(
                    model,
                    &seat.faces[i],
                    ogeom_math::Point2::new(u + turn, v),
                    tol,
                )?;
                if held {
                    break;
                }
            }
            if !held {
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

/// The exact tool of a round of `radius` on the ruled seat of `edge`: the
/// ball's section is the same all along, so the band is the cylinder of
/// the radius about the line its centre runs on, the legs lie on the two
/// hosts and the tool ends in the sections square to the edge at its
/// vertices, where a chain's junctions join it. `None` where the seat is
/// not ruled, or the tool will not close.
///
/// # Errors
///
/// As reading the seat and building the tool.
#[allow(clippy::too_many_lines, reason = "one tool, its pieces in turn")]
pub(crate) fn exact_wedge(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<crate::support::Wedge>> {
    use super::overflow::{ToolFace, closed_within, segment};
    use super::switch_over::outward_at;
    use ogeom_geom::Surface as _;
    use ogeom_math::{Circle, Direction, Frame, Plane, Point};
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
    let (t0, t1) = seat.range;
    let tm = 0.5 * (t0 + t1);
    let (a, b, m) = (
        seat.guide.point_at(t0, tol)?,
        seat.guide.point_at(t1, tol)?,
        seat.guide.point_at(tm, tol)?,
    );
    let length = a.distance(b);
    let along = (b - a) / length;
    let ua = ogeom_algo::project_on_surface(&seat.hosts[0].0, m, 32, tol)?.parameters;
    let ub = ogeom_algo::project_on_surface(&seat.hosts[1].0, m, 32, tol)?.parameters;
    let x = seat_section(
        &seat.hosts[0].0,
        &seat.hosts[1].0,
        radius,
        &seat.guide,
        sides,
        tm,
        [ua.0, ua.1, ub.0, ub.1],
        tol,
    )?;
    let hosts = [
        (seat.hosts[0].0.clone(), seat.hosts[0].1),
        (seat.hosts[1].0.clone(), seat.hosts[1].1),
    ];
    let p0 = hosts[0].0.point_at(x[0], x[1], tol)?;
    let p1 = hosts[1].0.point_at(x[2], x[3], tol)?;
    // The centre a radius off each host, outward of a concave crease and
    // within a convex one's material.
    let (_, n0) = outward_at(&hosts[0], p0, tol)?;
    let centre = p0 + n0 * (sign * radius);
    if (centre.distance(p1) - radius).abs() > crate::marched::band_fit_target(tol) {
        return Ok(None);
    }
    // The section is the same at every point of the crease.
    let (o0, o1, oc) = (p0 - m, p1 - m, centre - m);
    let points = vec![a, b, a + o0, b + o0, a + o1, b + o1];
    let (pa, pb, a0, b0, a1, b1) = (0, 1, 2, 3, 4, 5);
    let mut curves: Vec<super::curved::SheetCurve> = Vec::new();
    let push = |curves: &mut Vec<super::curved::SheetCurve>, ends, (curve, range)| {
        curves.push(super::curved::SheetCurve { ends, curve, range });
        curves.len() - 1
    };
    let crease = push(&mut curves, (pa, pb), (seat.guide.clone(), seat.range));
    let rail0 = push(&mut curves, (a0, b0), segment(points[a0], points[b0], tol)?);
    let rail1 = push(&mut curves, (a1, b1), segment(points[a1], points[b1], tol)?);
    // The ball's arc at each end, from the first host's contact to the
    // second's, in the section's plane.
    let arc = |at: Point| -> OgeomResult<(Curve, (f64, f64))> {
        let c = at + oc;
        let x = Direction::new(o0 - oc, tol)?;
        let mut z = along;
        let to = o1 - oc;
        let y = z.cross(x.vector());
        if y.dot(to) < 0.0 {
            z = -z;
        }
        let y = z.cross(x.vector());
        let sweep = y.dot(to).atan2(x.vector().dot(to));
        let frame = Frame::new(c, Direction::new(z, tol)?, x, tol)?;
        let circle = Circle::new(frame, radius, tol)?;
        Ok((
            Curve::Circle(ogeom_geom::CircleCurve::new(circle)),
            (0.0, sweep),
        ))
    };
    let arc_a = push(&mut curves, (a0, a1), arc(a)?);
    let arc_b = push(&mut curves, (b0, b1), arc(b)?);
    // Each end's legs where its section cuts the hosts.
    let mut legs = [[0; 2]; 2];
    for (k, (at, crease_point, contacts)) in [(a, pa, [a0, a1]), (b, pb, [b0, b1])]
        .into_iter()
        .enumerate()
    {
        for i in 0..2 {
            let Some(leg) = super::pipe_corner::leg(
                &hosts[i].0,
                (at, along),
                points[crease_point],
                points[contacts[i]],
                tol,
            )?
            else {
                return Ok(None);
            };
            legs[k][i] = push(&mut curves, (crease_point, contacts[i]), leg);
        }
    }
    let mid_dir = (o0 + o1) * 0.5 - oc;
    let mid_dir = mid_dir / mid_dir.magnitude();
    let band_frame = Frame::new(
        a + oc,
        Direction::new(along, tol)?,
        Direction::new(-mid_dir, tol)?,
        tol,
    )?;
    let band: SurfaceGeometry = ogeom_geom::CylinderSurface::new(
        ogeom_math::Cylinder::new(band_frame, radius, tol)?,
        (-radius, length + radius),
    )?
    .into();
    let on_band = m + oc + mid_dir * radius;
    let mut faces = vec![ToolFace {
        surface: band,
        curves: vec![rail0, arc_b, rail1, arc_a],
        at: on_band,
        out: -mid_dir,
    }];
    for (i, rail) in [rail0, rail1].into_iter().enumerate() {
        let o = if i == 0 { o0 } else { o1 };
        let (at, n) = outward_at(&hosts[i], m + o * 0.5, tol)?;
        faces.push(ToolFace {
            surface: crate::marched::ruled_on(&hosts[i].0, radius)?,
            curves: vec![crease, legs[1][i], rail, legs[0][i]],
            at,
            out: n * -sign,
        });
    }
    for (k, (at, outward, arc)) in [(a, -along, arc_a), (b, along, arc_b)]
        .into_iter()
        .enumerate()
    {
        let reach = (radius * 4.0).max(1.0);
        let plane = Plane::through(at, Direction::new(outward, tol)?);
        faces.push(ToolFace {
            surface: ogeom_geom::PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?
                .into(),
            curves: vec![legs[k][0], arc, legs[k][1]],
            at: at + (o0 + o1 + oc) * (1.0 / 3.0),
            out: outward,
        });
    }
    let widths = vec![crate::marched::band_fit_target(tol); curves.len()];
    let Some((tool, _)) = closed_within(model, &points, &curves, &faces, &widths, tol)? else {
        return Ok(None);
    };
    Ok(Some(crate::support::Wedge {
        solid: tool,
        additive: !convex,
        edges: vec![edge.clone()],
        melt_refusal: None,
        credit: Vec::new(),
    }))
}
