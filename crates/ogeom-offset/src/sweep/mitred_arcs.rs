//! A disc swept down an open spine of lines and arcs whose sharp corners
//! turn in the planes of the arcs beside them, exactly.
//!
//! Under the rotation-minimizing frame a corner turning in its arcs' plane
//! keeps the section's frame in that plane, and the turn at the corner is
//! the reflection in the plane halving it followed by the reflection in the
//! plane square to the next leg. A point of the disc at offset `u` outside
//! the turn and height `w` above its plane runs at that height along the
//! line offset by `u` from a line leg and round the circle of radius `R ± u`
//! about an arc leg's axis, and its paths down the two legs cross where the
//! pipe shell ends its walls:
//!
//! - outside the turn, on each leg's run straight on past the corner. The
//!   two run-ons are mirror images in the plane halving the corner, so
//!   they cross on it: each is the prism of its leg's end section, cut on
//!   that plane, and the two cuts are one half ellipse;
//! - inside the turn, on the legs themselves (see [`InsideTurn`]). There
//!   the slice of either leg at height `w` is the band of half-width
//!   `√(r² - w²)` about its track, whose edge is the path of the disc's
//!   rim point at that height, so a rim point's path enters the other leg
//!   exactly where it crosses that leg's own path for it: the legs' union
//!   ends each wall on that crossing.
//!
//! The pipe is the union of its legs (each the cylinder or torus segment
//! its section sweeps) and, at each corner, the two run-ons outside the
//! turn. The union is that pipe where every rim point inside a turn finds
//! its crossing within the first half of either leg, as the general sweep
//! asks; elsewhere `None` leaves the corner to it.

use ogeom_algo::History;
use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, PlaneSurface, SurfaceGeometry};
use ogeom_math::{Axis, Direction, Plane, Point, Transform, Vector};
use ogeom_topo::{Filter, Model, Shape, ShapeType, explore};

use super::inside_turn::InsideTurn;
use super::{
    Built, DOUBLED_BACK, Leg, LegKind, fuse_in_order, realized_profile_wound, spine_curve_of,
};

/// The pipe of the disc `profile` down the open spine of `legs`, its sharp
/// corners turning in the planes of the arcs beside them. `None` where the
/// profile is not a disc centred on the spine's start, the spine closes, a
/// corner turns an arc out of its plane or all but doubles back, the disc
/// is as wide as an arc bends, or a rim point inside a turn finds no
/// crossing, or finds it past the middle of a leg.
pub(super) fn mitred_arcs(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    legs: &[Leg],
    tol: Tolerances,
) -> OgeomResult<Option<Built>> {
    let (Some(first), Some(last)) = (legs.first(), legs.last()) else {
        return Ok(None);
    };
    let Some(r) = disc_radius(model, profile, first, tol)? else {
        return Ok(None);
    };
    let end_of = |leg: &Leg| leg.motion.apply(leg.start);
    if first.start.distance(end_of(last)) <= tol.confusion() * 100.0 {
        return Ok(None);
    }
    let unit = |v: Vector| v / v.magnitude();
    // Each arc's centre, radius and angle.
    let arc_of = |leg: &Leg| match leg.along {
        LegKind::Arc(axis, angle) => {
            let radial = leg.start - axis.location;
            let up = axis.direction.vector();
            Some((axis, (radial - up * radial.dot(up)).magnitude(), angle))
        }
        LegKind::Line(_) => None,
    };
    for leg in legs {
        if let Some((_, radius, _)) = arc_of(leg)
            && r >= radius - tol.confusion()
        {
            return Ok(None);
        }
    }

    // Each sharp corner: its point, the headings either side, its axis and
    // angle.
    struct Corner {
        at: Point,
        before: Vector,
        after: Vector,
        axis: Vector,
        angle: f64,
    }
    let mut corners: Vec<Option<Corner>> = Vec::with_capacity(legs.len() - 1);
    for pair in legs.windows(2) {
        let (d, n) = (unit(pair[0].heading.1), unit(pair[1].heading.0));
        let turn = d.cross(n);
        if turn.magnitude() <= tol.angular() && d.dot(n) > 0.0 {
            corners.push(None);
            continue;
        }
        if d.dot(n) < DOUBLED_BACK || turn.magnitude() <= tol.angular() {
            return Ok(None);
        }
        let up = unit(turn);
        let arcs = [&pair[0], &pair[1]].map(arc_of);
        if arcs
            .iter()
            .flatten()
            .any(|(axis, _, _)| axis.direction.vector().cross(up).magnitude() > tol.angular())
        {
            return Ok(None);
        }
        let at = pair[1].start;
        let angle = d.dot(n).clamp(-1.0, 1.0).acos();
        // The rim point deepest inside the turn reaches farthest along both
        // legs; between two lines its crossing is on the plane halving the
        // corner, `r·tan(θ/2)` along either.
        let crossing = if arcs.iter().all(Option::is_none) {
            at + unit(n - d * d.dot(n)) * r - d * (r * (angle * 0.5).tan())
        } else {
            let inside = InsideTurn::new(
                at,
                up,
                (d, n),
                arcs.map(|a| a.map(|(axis, ..)| axis.location)),
            );
            let outward = inside.offset_plane().1;
            match inside.meeting(at - outward * r) {
                Some(crossing) => crossing,
                None => return Ok(None),
            }
        };
        // Each leg's share the crossing takes, at most half of it, as the
        // general sweep asks.
        for (leg, heading) in [(&pair[0], -d), (&pair[1], n)] {
            let depth = (crossing - at).dot(heading);
            if depth < 0.0 {
                return Ok(None);
            }
            let share = match arc_of(leg) {
                // The crossing's turn about the arc's axis, against the arc's.
                Some((axis, _, arc_angle)) => {
                    let centre = axis.location + up * (at - axis.location).dot(up);
                    let (a, b) = (crossing - centre, at - centre);
                    a.cross(b).magnitude().atan2(a.dot(b)) / arc_angle
                }
                None => depth / leg_length(leg),
            };
            if share >= 0.5 {
                return Ok(None);
            }
        }
        corners.push(Some(Corner {
            at,
            before: d,
            after: n,
            axis: up,
            angle,
        }));
    }
    // The legs, and the run-ons outside each corner: the prism of the
    // section where the corner stands, run past it straight along each leg,
    // cut on the plane halving the corner.
    let mut carried = Transform::IDENTITY;
    let mut pieces: Vec<Shape> = Vec::with_capacity(legs.len() * 3);
    for (i, leg) in legs.iter().enumerate() {
        let section = realized_profile_wound(model, profile, &carried, Some(leg.heading.0), tol)?;
        pieces.push(match leg.along {
            LegKind::Line(v) => ogeom_algo::make_prism(model, &section, v, tol)?.shape,
            LegKind::Arc(axis, angle) => {
                ogeom_algo::make_revolution(model, &section, axis, angle, tol)?.shape
            }
        });
        carried = leg.motion * carried;
        let Some(Some(corner)) = corners.get(i) else {
            continue;
        };
        let run = r * ((corner.angle * 0.5).tan() + 1.0);
        let halving = unit(corner.before + corner.after);
        // The leg's run-on, kept behind the plane.
        let end = realized_profile_wound(model, profile, &carried, Some(corner.before), tol)?;
        let ahead = ogeom_algo::make_prism(model, &end, corner.before * run, tol)?.shape;
        pieces.push(cut_on(model, &ahead, corner.at, -halving, r + run, tol)?);
        // Turned onto the next leg by the least rotation between headings.
        carried = Transform::rotation(
            Axis {
                location: corner.at,
                direction: Direction::new(corner.axis, tol)?,
            },
            corner.angle,
        ) * carried;
        // The next leg's run-back, kept ahead of the plane.
        let back = realized_profile_wound(
            model,
            profile,
            &(Transform::translation(-corner.after * run) * carried),
            Some(corner.after),
            tol,
        )?;
        let behind = ogeom_algo::make_prism(model, &back, corner.after * run, tol)?.shape;
        pieces.push(cut_on(model, &behind, corner.at, halving, r + run, tol)?);
    }
    let Some(shape) = fuse_in_order(model, pieces, tol)? else {
        return Ok(None);
    };
    let mut history = History::new();
    history.generate(spine, shape.clone());
    history.generate(profile, shape.clone());
    Ok(Some(Built::new(shape, history)))
}

/// The radius of `profile` where it is a disc: one ring of one circle,
/// centred on the start of `first` and square to its heading there.
fn disc_radius(
    model: &Model,
    profile: &Shape,
    first: &Leg,
    tol: Tolerances,
) -> OgeomResult<Option<f64>> {
    let wires = explore(model, profile, Filter::OfType(ShapeType::Wire))?;
    let [wire] = wires.as_slice() else {
        return Ok(None);
    };
    let edges = model.ordered_children_of(wire)?;
    let [edge] = edges.as_slice() else {
        return Ok(None);
    };
    if !edge.location().is_identity() || !profile.location().is_identity() {
        return Ok(None);
    }
    let (curve, range) = spine_curve_of(model, edge)?;
    let Curve::Circle(c) = curve else {
        return Ok(None);
    };
    let circle = c.circle();
    let heading = first.heading.0 / first.heading.0.magnitude();
    let whole = (range.1 - range.0 - core::f64::consts::TAU).abs() <= tol.angular();
    let centred = circle.frame().origin().distance(first.start) <= tol.confusion();
    let square = circle.frame().z().vector().cross(heading).magnitude() <= tol.angular();
    Ok((whole && centred && square).then(|| circle.radius()))
}

/// A leg's length along its spine.
fn leg_length(leg: &Leg) -> f64 {
    match leg.along {
        LegKind::Line(v) => v.magnitude(),
        LegKind::Arc(axis, angle) => {
            let radial = leg.start - axis.location;
            let up = axis.direction.vector();
            (radial - up * radial.dot(up)).magnitude() * angle
        }
    }
}

/// What of `solid` lies on the side `keep` points to of the plane through
/// `at` square to it, `reach` bounding the solid's extent from `at`.
fn cut_on(
    model: &mut Model,
    solid: &Shape,
    at: Point,
    keep: Vector,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let plane = Plane::through(at, Direction::new(keep, tol)?);
    let span = (-reach * 4.0, reach * 4.0);
    let face = ogeom_algo::make_natural_face(
        model,
        SurfaceGeometry::Plane(PlaneSurface::over(plane, span, span)?),
    )?
    .shape;
    let half = ogeom_algo::make_half_space(model, &face, at + keep * reach, tol)?.shape;
    Ok(ogeom_bool::common(model, solid, &half, tol)?.shape)
}
