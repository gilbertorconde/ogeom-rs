//! The exact round of an arc of a rim: the torus band of a revolved seat
//! over the arc's angular window, closed by planar end caps in the meridian
//! half-planes at its ends. A strip builder shared by the engine's marched
//! seats (a rim meeting a marched crease), the round distances on rims and
//! the one-edge entries.

use crate::fillet::Mate;
use crate::support::{edge_curve, face_from_edges, revolved_seat, segment_between};
use ogeom_algo::{Built, make_edge_between, make_vertex};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{CircleCurve, Curve, PlaneSurface, SurfaceGeometry, TorusSurface};
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Torus, Vector};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

/// The rim fillet over an arc of the rim rather than its whole turn.
///
/// The same four-seat torus geometry as the whole rim's `revolved_fillet`, restricted to
/// the arc's angular window and closed by two planar end caps in the
/// meridian half-planes at its ends. Every face is a chart rectangle or a
/// planar triangle, every edge an arc or a segment, every pcurve exact,
/// which is what lets a chain's neighbouring wedges melt at their shared
/// caps.
#[allow(clippy::too_many_arguments, reason = "the seat's own, and the section")]
pub(crate) fn revolved_arc_fillet(
    model: &mut Model,
    solid: &Shape,
    edge: &Shape,
    rim: &CircleCurve,
    radius: f64,
    conic: Option<&crate::asym::RimConic>,
    mates: Option<(usize, &[Mate])>,
    aside: crate::support::Aside<'_>,
    tol: Tolerances,
) -> OgeomResult<Built> {
    use ogeom_geom::Curve3d as _;
    let seat = revolved_seat(model, solid, edge, rim, tol)?;

    // The ball's place in the meridian, or the conic section's contacts;
    // and how far the section reaches, which scales its run-out.
    let ball = conic.map_or_else(|| seat.ball(radius), |c| c.contacts(&seat));
    let radius = conic.map_or(radius, crate::asym::RimConic::size);
    let tube_rho = ball.rho;
    if tube_rho <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "a fillet of radius {radius} swallows the axis of a rim of \
             radius {}",
            seat.radius
        );
    }
    let tube_level = seat.centre + seat.up * ball.height;
    let wall_level = seat.centre + seat.up * ball.wall_height;

    // The arc's window in the seat's own parameterization: the edge's stored
    // circle may be framed anyhow, so the ends are read back through the
    // seat frame and the window is the one that holds the arc's midpoint.
    let seat_frame = seat.frame_at(seat.centre, tol)?;
    let (curve, crange) = edge_curve(model, edge, tol)?;
    let angle_of = |t: f64| -> OgeomResult<f64> {
        let local = seat_frame.to_local(curve.point_at(t, tol)?);
        Ok(local.y.atan2(local.x))
    };
    let tau_turn = core::f64::consts::TAU;
    let raw0 = angle_of(crange.0)?;
    let raw1 = angle_of(crange.1)?;
    let mid = angle_of(f64::midpoint(crange.0, crange.1))?;
    let (theta0, theta1) = {
        let span = (raw1 - raw0).rem_euclid(tau_turn);
        let inside = (mid - raw0).rem_euclid(tau_turn) <= span;
        if inside {
            (raw0, raw0 + span)
        } else {
            (raw1, raw1 + (tau_turn - span))
        }
    };
    if theta1 - theta0 <= tol.angular() {
        ogeom_bail!(Construction, "the rim arc has no length to blend along");
    }
    let dir_of = |theta: f64| -> Vector {
        seat_frame.x().vector() * theta.cos() + seat_frame.y().vector() * theta.sin()
    };
    // At an end where the meridian cap would stop short of an oblique face
    // the arc ends on, the band runs on past the end, its window widened
    // until the ball's section stands clear of the material (a plane leaves
    // the section's triangle last at a corner: the rim, the wall's contact
    // and the cap's), and a quarter radius further; the cut trims it
    // against the face, as a straight band's run-out is.
    let mut trims: Vec<(Point, Vector)> = Vec::new();
    let (theta0, theta1) = {
        let mut window = [theta0, theta1];
        let step = (radius / seat.radius) / 8.0;
        for (slot, sign) in [(0_usize, -1.0), (1, 1.0)] {
            let end = window[slot];
            let at = seat.centre + dir_of(end) * seat.radius;
            let ahead = seat.up.cross(dir_of(end)) * sign;
            let Some(face) = rim_end_stops_short(model, solid, edge, &seat, mates, at, ahead, tol)?
            else {
                continue;
            };
            let probe = crate::support::material_probe(model, solid, tol)?;
            let (on, back) = face;
            let mut past = None;
            let mut k = 1.0_f64;
            // Over four radii along the rim, as a straight band's run-out
            // is read: a section that has not left the material by then, or
            // whose corner leaves it other than across the face the arc ends
            // on, meets a step or a ledge past the end rather than a ramp
            // across it, and the band stops at the end.
            while k * step * seat.radius <= radius * 4.0 {
                let theta = end + sign * k * step;
                // The corners stand on the hosts; each is read a hundredth
                // of the way in towards the others, inside the corner.
                let (apex, wall, cap) = (
                    seat.centre + dir_of(theta) * seat.radius,
                    wall_level + dir_of(theta) * ball.wall_rho,
                    seat.centre + dir_of(theta) * tube_rho,
                );
                let section = [
                    apex + ((wall - apex) + (cap - apex)) * 0.01,
                    wall + (apex - wall) * 0.01 + (cap - wall) * 0.01,
                    cap + (apex - cap) * 0.01 + (wall - cap) * 0.01,
                ];
                let mut clear = true;
                let mut ledge = false;
                for p in section {
                    if probe.holds(model, p, tol)? == ogeom_algo::Containment::In {
                        clear = false;
                    } else if (p - on).dot(back) > tol.confusion() {
                        ledge = true;
                    }
                }
                if ledge {
                    break;
                }
                if clear {
                    past = Some(k * step + 2.0 * step);
                    break;
                }
                k += 1.0;
            }
            if let Some(past) = past {
                window[slot] = end + sign * past;
                trims.push(face);
            }
        }
        if window[1] - window[0] >= tau_turn - tol.angular() {
            ogeom_bail!(
                Construction,
                "a rim arc's band run on past both its ends would close the turn"
            );
        }
        (window[0], window[1])
    };
    let sweep = theta1 - theta0;

    // Corner points at each end: A on the rim, W down the wall, C in on the
    // cap; and the frames the arcs run on.
    let dir_at = |theta: f64| -> Vector {
        seat_frame.x().vector() * theta.cos() + seat_frame.y().vector() * theta.sin()
    };
    // Every wedge chart is anchored at the arc's start: the window then
    // lives in (0, sweep), clear of the chart's own seam whatever part of
    // the rim the arc covers.
    let wedge_frame_at = |origin: Point| -> OgeomResult<Frame> {
        Frame::new(
            origin,
            Direction::new(seat.up, tol)?,
            Direction::new(dir_at(theta0), tol)?,
            tol,
        )
    };
    let apex_at = |theta: f64| seat.centre + dir_at(theta) * seat.radius;
    let wall_at = |theta: f64| wall_level + dir_at(theta) * ball.wall_rho;
    let cap_at = |theta: f64| seat.centre + dir_at(theta) * tube_rho;

    // Six shared corner vertices, nine shared edges: a wire chains through
    // vertex objects, and neighbouring faces sew through shared edges.
    let va0 = make_vertex(model, apex_at(theta0)).shape;
    let va1 = make_vertex(model, apex_at(theta1)).shape;
    let vw0 = make_vertex(model, wall_at(theta0)).shape;
    let vw1 = make_vertex(model, wall_at(theta1)).shape;
    let vc0 = make_vertex(model, cap_at(theta0)).shape;
    let vc1 = make_vertex(model, cap_at(theta1)).shape;

    let ring_arc = |model: &mut Model,
                    origin: Point,
                    r: f64,
                    from: &Shape,
                    to: &Shape|
     -> OgeomResult<Shape> {
        let circle = Circle::new(wedge_frame_at(origin)?, r, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        Ok(make_edge_between(model, curve, (0.0, sweep), from, to, tol)?.shape)
    };
    // The ball's arc in the meridian at `theta`, from the wall contact to
    // the cap contact: the tube's own circle, in the meridian plane, a
    // quarter of it beside a cylinder.
    let quarter = |model: &mut Model, theta: f64, from: &Shape, to: &Shape| -> OgeomResult<Shape> {
        if let Some(c) = conic {
            return Ok(make_edge_between(
                model,
                c.meridian(&seat, theta, tol)?,
                (0.0, 1.0),
                from,
                to,
                tol,
            )?
            .shape);
        }
        let dir = dir_at(theta);
        let centre_m = tube_level + dir * tube_rho;
        let x_dir = Direction::new((dir - seat.up * seat.slope) * (seat.sigma * seat.tau), tol)?;
        let z_dir = Direction::new(dir.cross(seat.up) * seat.sigma, tol)?;
        let frame = Frame::new(centre_m, z_dir, x_dir, tol)?;
        let circle = Circle::new(frame, radius, tol)?;
        let curve = Curve::Circle(CircleCurve::new(circle));
        Ok(make_edge_between(model, curve, (0.0, ball.sweep), from, to, tol)?.shape)
    };

    let arc_apex = ring_arc(model, seat.centre, seat.radius, &va0, &va1)?;
    let arc_wall = ring_arc(model, wall_level, ball.wall_rho, &vw0, &vw1)?;
    let arc_cap = ring_arc(model, seat.centre, tube_rho, &vc0, &vc1)?;
    let ruling0 = segment_between(model, (&va0, apex_at(theta0)), (&vw0, wall_at(theta0)), tol)?;
    let ruling1 = segment_between(model, (&vw1, wall_at(theta1)), (&va1, apex_at(theta1)), tol)?;
    let radial0 = segment_between(model, (&vc0, cap_at(theta0)), (&va0, apex_at(theta0)), tol)?;
    let radial1 = segment_between(model, (&va1, apex_at(theta1)), (&vc1, cap_at(theta1)), tol)?;
    let quarter0 = quarter(model, theta0, &vw0, &vc0)?;
    let quarter1 = quarter(model, theta1, &vw1, &vc1)?;
    // Each face keeps itself on the left of its ring in its own chart: a
    // ring given turning the other way is walked back.
    let walked = |edges: Vec<Shape>, back: bool| -> Vec<Shape> {
        if back {
            edges.iter().rev().map(Shape::reversed).collect()
        } else {
            edges
        }
    };

    // The wall patch: the solid wall's own chart between the contact ring
    // and the rim, over the window, its chart reaching past the rim on
    // both sides whichever side of the cap the wall stands on. The ring
    // runs the contact ring forward and the rim back, anticlockwise in the
    // chart where the contact ring stands below the rim.
    let wall_patch = {
        let below = radius + 1.0;
        let surface = seat.wall_surface(
            (-below, below),
            (ball.wall_height.min(0.0), ball.wall_height.max(0.0)),
            Direction::new(dir_at(theta0), tol)?,
            tol,
        )?;
        let edges = walked(
            vec![
                arc_wall.clone(),
                ruling1.clone(),
                arc_apex.reversed(),
                ruling0.clone(),
            ],
            seat.tau < 0.0,
        );
        let face = face_from_edges(model, surface, &edges, tol)?;
        if seat.sigma * seat.tau < 0.0 {
            face.reversed()
        } else {
            face
        }
    };

    // The cap sector: the cap plane between the rim and the contact ring.
    // The ring runs the rim forward about `up`, which turns it about the
    // plane's normal `up * tau` where the rim is the wider ring and the cap
    // faces up, or the narrower and it faces down.
    let cap_sector = {
        let plane = Plane::through(seat.centre, Direction::new(seat.up * seat.tau, tol)?);
        let reach = (seat.radius + radius + 1.0) * 2.0;
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let rim_wider = seat.radius > tube_rho;
        let edges = walked(
            vec![
                arc_apex.clone(),
                radial1.clone(),
                arc_cap.reversed(),
                radial0.clone(),
            ],
            rim_wider != (seat.tau > 0.0),
        );
        face_from_edges(model, surface, &edges, tol)?
    };

    // The blend: the quarter-tube patch over the window, reversed as the
    // full rim's band is: its natural normal points into the wedge. The
    // ring runs the wall's contact forward and the cap's back, anticlockwise
    // in the torus chart where the tube's angle grows from the wall's contact
    // to the cap's, which is where `sigma` is positive.
    let blend_patch = if let Some(c) = conic {
        c.band_sector(
            model,
            &seat,
            (theta0, sweep),
            [&arc_wall, &arc_cap, &quarter0, &quarter1],
            tol,
        )?
    } else {
        let surface: SurfaceGeometry = TorusSurface::new(Torus::new(
            wedge_frame_at(tube_level)?,
            tube_rho,
            radius,
            tol,
        )?)
        .into();
        let edges = walked(
            vec![
                arc_wall.clone(),
                quarter1.clone(),
                arc_cap.reversed(),
                quarter0.reversed(),
            ],
            seat.sigma < 0.0,
        );
        face_from_edges(model, surface, &edges, tol)?.reversed()
    };

    // The end caps: planar triangles in the meridian half-planes, outward
    // along the rim's travel: behind it at the start, ahead at the end.
    // Each ring is walked back where the corners it passes in turn
    // against `outward`, so the cap keeps itself on the left: the quarter
    // arc bulges toward the rim corner and never past it, so the corners'
    // turn is the ring's.
    let tangent_at = |theta: f64| seat.up.cross(dir_at(theta));
    let end_cap = |model: &mut Model,
                   theta: f64,
                   outward: Vector,
                   edges: Vec<Shape>,
                   [p, q, r]: [Point; 3]|
     -> OgeomResult<Shape> {
        let plane = Plane::through(apex_at(theta), Direction::new(outward, tol)?);
        let reach = (seat.radius + radius + 1.0) * 2.0;
        let surface: SurfaceGeometry =
            PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
        let edges = walked(edges, (q - p).cross(r - p).dot(outward) < 0.0);
        face_from_edges(model, surface, &edges, tol)
    };
    let cap_start = end_cap(
        model,
        theta0,
        -tangent_at(theta0),
        vec![ruling0.clone(), quarter0.clone(), radial0.clone()],
        [apex_at(theta0), wall_at(theta0), cap_at(theta0)],
    )?;
    let cap_end = end_cap(
        model,
        theta1,
        tangent_at(theta1),
        vec![ruling1.clone(), radial1.clone(), quarter1.reversed()],
        [wall_at(theta1), apex_at(theta1), cap_at(theta1)],
    )?;

    let faces = [wall_patch, cap_sector, blend_patch, cap_start, cap_end];
    if !seat.additive() || trims.is_empty() {
        return crate::support::finish_wedge(
            model,
            solid,
            Some(edge),
            &faces,
            seat.additive(),
            aside,
            tol,
        );
    }
    let mut wedge = crate::support::wedge_solid(model, &faces, false, tol)?;
    let span = (seat.radius + radius) * 4.0;
    for plane in trims {
        wedge = crate::support::trimmed_wedge(model, &wedge, plane, span, tol)?;
    }
    crate::support::finish_wedge_solid(model, solid, Some(edge), &wedge, true, aside, tol)
}

/// Whether an open rim arc's torus band, capped in the meridian half-plane
/// through its end `at` (the crease leaving along `ahead` past it), would
/// stop short of a face the arc ends on: the crease terminates there, no
/// mate of the request ends there, and a planar face there stands oblique
/// to the crease, its normal neither square to the crease's direction (a
/// wall continuing the crease tangentially, which meets the cap square) nor
/// along it (a wall square to the crease, which is the cap). Past such an
/// end the ball still rolls through material between the meridian cap and
/// the face while its contact on the cap stands on it. The face's plane,
/// facing back over the arc, where it does.
#[allow(clippy::too_many_arguments, reason = "one question, all its data")]
fn rim_end_stops_short(
    model: &Model,
    solid: &Shape,
    edge: &Shape,
    seat: &crate::support::RevolvedSeat,
    mates: Option<(usize, &[Mate])>,
    at: Point,
    ahead: Vector,
    tol: Tolerances,
) -> OgeomResult<Option<(Point, Vector)>> {
    let hosts = [&seat.cap_face, &seat.wall_face];
    let met = mates.is_some_and(|(index, mates)| {
        mates.iter().enumerate().any(|(i, mate)| {
            i != index
                && mate
                    .ends
                    .iter()
                    .any(|(p, _)| p.distance(at) <= tol.confusion() * 1e3)
        })
    });
    if met || !crate::marched::crease_terminates_at(model, solid, edge, hosts, at, tol)? {
        return Ok(None);
    }
    for face in explore_unique(model, solid, ShapeType::Face)? {
        if hosts.iter().any(|h| h.is_same(&face)) {
            continue;
        }
        let Some(data) = model.node(&face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        let Some(SurfaceGeometry::Plane(plane)) = model.geometry().surface(data.surface) else {
            continue;
        };
        use ogeom_geom::Transformable as _;
        let SurfaceGeometry::Plane(placed) =
            SurfaceGeometry::Plane(*plane).transformed(&face.transform(model.datums())?, tol)?
        else {
            continue;
        };
        let plane = placed.plane();
        if plane.signed_distance_to(at).abs() > tol.confusion() * 10.0 {
            continue;
        }
        let meets = explore_unique(model, &face, ShapeType::Vertex)?
            .iter()
            .any(|v| {
                model
                    .node(v)
                    .and_then(|n| n.data().as_vertex())
                    .zip(v.transform(model.datums()).ok())
                    .is_some_and(|(d, place)| {
                        place.apply(d.point).distance(at) <= tol.confusion() * 1e3
                    })
            });
        let normal = plane.normal().vector();
        let lean = normal.dot(ahead).abs();
        if meets && lean > 1e-6 && lean < 1.0 - 1e-6 {
            // Facing back over the arc: the side of the face the arc lies on.
            let back = if normal.dot(ahead) > 0.0 {
                -normal
            } else {
                normal
            };
            return Ok(Some((plane.origin(), back)));
        }
    }
    Ok(None)
}
