//! An open wire swept down an open spine of lines and arcs, as one sheet,
//! or round a closed one, as one ring.
//!
//! An open wire bounds no face, so there is no solid to fuse round a
//! corner and take the walls of. The sheet is built from its walls. The
//! rotation-minimizing frame carries the section to each leg: along a line
//! by its displacement, round an arc by the turn about the arc's axis, and
//! at each sharp corner by the least rotation between the headings. Each
//! line leg is the prism of its section, each arc leg the section revolved
//! about the arc's axis, and a corner is mitred or rounded on the plane
//! halving it:
//!
//! - mitred, each leg runs on straight past the corner, a line leg by its
//!   own prism and an arc leg by the prism of its end section, and is cut
//!   there. The turn is the reflection in that plane followed by the
//!   reflection in the plane square to the next leg, so the two run-ons
//!   are mirror images in the plane and a profile point's paths down them
//!   meet it at one point: the two legs' cuts are one curve;
//! - round, the leg ends square at the corner, and the plane cuts only the
//!   part of the section inside the turn. The part outside is turned about
//!   the corner's axis through the corner's angle, onto the next leg's
//!   section there; a profile edge along the axis turns on itself and
//!   sweeps nothing.
//!
//! A point of the section inside a turn reaches the other leg's wall on
//! the legs themselves, not on their run-ons, mitred or round alike: each
//! leg runs on until it meets the other's wall, which is where the two
//! paths of each point cross. Two lines are mirror images there too, so
//! their paths meet on the plane. A line and an arc, or two arcs, are not:
//!
//! - where the corner turns in the arc's plane the frame keeps the
//!   section's frame in that plane, and a point's two paths cross off the
//!   plane, on a conic for each straight profile edge (see
//!   [`InsideTurn`]). The profile's edges are cut where they cross the
//!   corner's axis and swept one by one, and each edge inside the corner
//!   has both its walls split on that one curve, so they share it;
//!   A curved profile edge's walls there meet on a curve with no closed
//!   form: each of its points' crossing is found exactly, and the edge
//!   both walls are split on is fitted through them, its error stated as
//!   its tolerance;
//! - where the corner turns an arc out of its plane, a point's two paths
//!   pass each other without meeting, and the walls leave a gap only a
//!   face of the plane halving the corner could close, as a face's sweep
//!   closes it. A point there is refused.
//!
//! Round a closed spine the sheet is swept from the middle of the first
//! leg round to there, the profile carried there as the leg carries it,
//! so the corner where the spine closes is turned as any other and the
//! sheet closes on itself where it runs on smoothly. It closes only where
//! the frame comes home unturned: a spine turning in one plane, or out of
//! it by turns that cancel. Elsewhere the ring is refused.
//!
//! Each leg lies between its corners' planes and each round between the
//! planes square to its legs, so neighbours meet only on the edges they
//! share. The pieces are sewn and the sheet must come out a disc: one
//! shell facing one way, every edge of one face or two, one boundary loop
//! running only along the two end sections and the paths of the profile's
//! two ends, and no two faces meeting other than along what they share.
//!
//! The Frenet, binormal and auxiliary laws place the section themselves.
//! Where every frame they set is the rotation-minimizing one turned by one
//! fixed angle, they place it as that frame carries it, and sweep the same
//! sheet; elsewhere their sections at a corner differ, and the two legs'
//! cuts with them.

use ogeom_algo::{
    Built, History, check, check_self_intersection, edge_vertices, make_half_space,
    make_natural_face, make_prism, make_revolution, make_wire, sew_within, shape_bounds,
};
use ogeom_core::{FastMap, OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _, PlaneSurface, SurfaceGeometry, Transformable as _};
use ogeom_math::{Axis, Direction, Plane, Point, Transform, Vector};
use ogeom_topo::{
    EdgeRepr, Filter, Location, Model, Shape, ShapeType, TShapeId, explore, explore_unique,
};

use super::inside_turn::InsideTurn;
use super::{
    DOUBLED_BACK, PipeCorners, PipeLaw, SpineStation, frenet_normals, law_normals, realized_rings,
    rmf_normals, spine_curve_of, wire_ends,
};

/// How a leg carries the section: down a line by its displacement, or
/// round an arc about its axis through its angle.
#[derive(Debug, Clone, Copy)]
enum Along {
    Line(Vector),
    Arc(Axis, f64),
}

/// A leg of the spine where it stands: its ends, its unit headings there
/// in the direction of travel, and how it carries the section.
#[derive(Debug, Clone, Copy)]
struct Leg {
    start: Point,
    end: Point,
    heading_in: Vector,
    heading_out: Vector,
    along: Along,
}

impl Leg {
    /// The motion carrying the section from the leg's start to its end.
    fn motion(&self) -> Transform {
        match self.along {
            Along::Line(v) => Transform::translation(v),
            Along::Arc(axis, angle) => Transform::rotation(axis, angle),
        }
    }
}

/// A leg's turn at its far corner: the corner's axis and the angle, where
/// the leg turns.
type Turn = Option<(Axis, f64)>;

/// Whether every edge of `spine` is a line or an arc of a circle.
pub(super) fn lines_and_arcs(model: &Model, spine: &Shape) -> OgeomResult<bool> {
    for edge in spine_edges(model, spine)? {
        let (curve, _) = spine_curve_of(model, &edge)?;
        let basis = match &curve {
            Curve::Trimmed(t) => t.basis().clone(),
            other => other.clone(),
        };
        if !matches!(basis, Curve::Line(_) | Curve::Circle(_)) {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The spine's edges in the order it runs them.
fn spine_edges(model: &Model, spine: &Shape) -> OgeomResult<Vec<Shape>> {
    Ok(match model.kind_of(spine)? {
        ShapeType::Edge => vec![spine.clone()],
        _ => model.ordered_children_of(spine)?,
    })
}

/// Sweep an open wire down an open spine of lines and arcs round sharp
/// corners, or round a closed one, under a frame law other than the
/// rotation-minimizing one: as [`open_wire_legs`] where the law's frames
/// are the rotation-minimizing ones turned by one fixed angle, so that it
/// places every section as that frame carries it. The Frenet frame is read run by run between the sharp
/// corners, as a face's runs are swept, and is refused along a straight
/// run, which has none.
#[allow(clippy::too_many_arguments)]
pub(super) fn open_wire_law(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    corners: PipeCorners,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let legs = placed_legs(model, spine, tol)?;
    let sharp = |k: usize| {
        let (a, b) = (legs[k].heading_out, legs[k + 1].heading_in);
        a.cross(b).magnitude() > tol.angular() || a.dot(b) < 0.0
    };
    let stations = placed_stations(model, spine, &legs)?;
    let normals = match law {
        PipeLaw::Frenet => {
            // Each run between sharp corners has its own Frenet frame.
            let mut normals = Vec::with_capacity(stations.len());
            let mut from = 0;
            for k in 0..legs.len() {
                if k + 1 < legs.len() && !sharp(k) {
                    continue;
                }
                let run: Vec<SpineStation> = stations
                    .iter()
                    .filter(|s| (from..=k).contains(&s.edge))
                    .copied()
                    .collect();
                if legs[from..=k]
                    .iter()
                    .all(|l| matches!(l.along, Along::Line(_)))
                {
                    ogeom_bail!(
                        Construction,
                        "the Frenet frame is not defined along a straight run \
                         between sharp corners"
                    );
                }
                normals.extend(frenet_normals(&run, tol)?);
                from = k + 1;
            }
            normals
        }
        _ => law_normals(model, &stations, law, tolerance, tol)?,
    };
    if !turned_rotation_minimizing(&stations, &normals, tol) {
        ogeom_bail!(
            Construction,
            "the law turns an open wire's section other than the \
             rotation-minimizing frame does; its legs' sections differ at a \
             corner, and their cuts with them, so they do not join into one \
             sheet"
        );
    }
    open_wire_legs(model, profile, spine, corners, tolerance, tol)
}

/// Whether the frames `normals` at `stations` are the rotation-minimizing
/// frames turned by one fixed angle: the angle between the law's normal
/// and the carried one the same at every station, to an angle's
/// tolerance.
pub(super) fn turned_rotation_minimizing(
    stations: &[SpineStation],
    normals: &[Vector],
    tol: Tolerances,
) -> bool {
    let carried = rmf_normals(stations);
    let angle = |i: usize| -> f64 {
        let (r, n) = (carried[i], normals[i]);
        r.cross(n).dot(stations[i].tangent).atan2(r.dot(n))
    };
    let first = angle(0);
    (0..stations.len()).all(|i| {
        let d = (angle(i) - first).rem_euclid(2.0 * core::f64::consts::PI);
        d.min(2.0 * core::f64::consts::PI - d) <= tol.angular()
    })
}

/// Stations along the spine where it stands: each leg's ends, and points
/// between them at most a sixty-fourth of a turn apart round an arc; a corner's
/// two stations stand on one point with the two legs' headings.
fn placed_stations(model: &Model, spine: &Shape, legs: &[Leg]) -> OgeomResult<Vec<SpineStation>> {
    let edges = spine_edges(model, spine)?;
    let mut stations = Vec::new();
    for (k, (edge, leg)) in edges.iter().zip(legs).enumerate() {
        let count: u32 = match leg.along {
            Along::Line(_) => 1,
            #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
            Along::Arc(_, angle) => {
                (angle / (core::f64::consts::TAU / 64.0)).ceil().max(1.0) as u32
            }
        };
        let (_, range) = spine_curve_of(model, edge)?;
        for i in 0..=count {
            let f = f64::from(i) / f64::from(count);
            let (at, tangent) = match leg.along {
                Along::Line(v) => (leg.start + v * f, leg.heading_in),
                Along::Arc(axis, angle) => {
                    let turn = Transform::rotation(axis, angle * f);
                    (turn.apply(leg.start), turn.apply_vector(leg.heading_in))
                }
            };
            stations.push(SpineStation {
                at,
                tangent,
                edge: k,
                t: range.0 + (range.1 - range.0) * f,
            });
        }
    }
    Ok(stations)
}

/// Sweep an open wire down an open spine of lines and arcs, or round a
/// closed one, under the rotation-minimizing frame, each corner mitred or
/// rounded; see the module documentation. The sheet faces the way the
/// prism of the profile along the first leg's heading does.
pub(super) fn open_wire_legs(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    corners: PipeCorners,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let given = profile;
    let legs = placed_legs(model, spine, tol)?;
    let unit = |v: Vector| v / v.magnitude();
    let sharp = |a: Vector, b: Vector| a.cross(b).magnitude() > tol.angular() || a.dot(b) < 0.0;
    let Some(first) = legs.first().copied() else {
        ogeom_bail!(Construction, "the spine has no leg");
    };
    match ogeom_algo::find_plane(model, profile, tol)? {
        Some(plane) => plane_square_at_start(&plane, first.start, first.heading_in, tol)?,
        None => lines_square_at_start(model, profile, first.start, first.heading_in, tol)?,
    }
    // Round a closed spine the sheet is swept from the middle of the first
    // leg, where the spine runs on smoothly, round to there again: the
    // profile is carried there as the leg carries it, and the corner where
    // the spine closes, if it turns there, is turned as any other.
    let closed = legs[legs.len() - 1].end.distance(first.start) <= tol.confusion();
    let (legs, section) = if closed {
        let (before, after) = halved(&first);
        let section = ogeom_algo::transformed(model, profile, before.motion())?.shape;
        let mut ring = vec![after];
        ring.extend_from_slice(&legs[1..]);
        ring.push(before);
        (ring, section)
    } else {
        (legs, profile.clone())
    };
    let profile = &section;
    if corners == PipeCorners::Extended
        && legs
            .windows(2)
            .any(|pair| sharp(pair[0].heading_out, pair[1].heading_in))
    {
        ogeom_bail!(
            Construction,
            "an extended corner squares the outside of a solid's section; an open \
             wire sweeps a sheet, whose corners are mitred or round"
        );
    }
    let first = legs[0];
    let (carried, turns) = carried_down(&legs, tol)?;
    // The motion carrying the profile to each leg's end.
    let carried_end: Vec<Transform> = legs
        .iter()
        .zip(&carried)
        .map(|(leg, motion)| leg.motion() * *motion)
        .collect();
    let hull = hull_points(model, profile, tol)?;
    // The frame carried once round a closed spine comes home turned by the
    // loop's holonomy, but for a spine turning in one plane, or as such
    // turns cancel; the sheet closes onto itself only where it comes home
    // unturned, every point of the section where it set out. How far a
    // point is moved is convex in the point, so the hull's points bound
    // it.
    if closed {
        let home = carried_end[carried_end.len() - 1];
        if hull
            .iter()
            .any(|p| home.apply(*p).distance(*p) > tol.confusion())
        {
            ogeom_bail!(
                Construction,
                "the rotation-minimizing frame comes home turned round this \
                 closed spine, so an open wire's sheet round it would not \
                 close onto itself"
            );
        }
    }
    let width = {
        let bounds = shape_bounds(model, profile, tol)?;
        match (bounds.low(), bounds.high()) {
            (Some(a), Some(b)) => a.distance(b),
            _ => ogeom_bail!(Construction, "the profile has no extent"),
        }
    };
    let straight = |i: usize| matches!(legs[i].along, Along::Line(_));
    let last = legs.len() - 1;
    // Whether each corner is rounded: a corner that turns nothing has no
    // round, and is cut on the plane square to both legs.
    let rounded: Vec<bool> = (0..last)
        .map(|k| corners == PipeCorners::Round && turns[k].is_some())
        .collect();
    // How far the section reaches outside and inside each turning corner:
    // its farthest points from the corner's axis on either side, in the
    // plane square to the leg there.
    let mut inside_turns: Vec<Option<InsideTurn>> = vec![None; last];
    for k in 0..last {
        let Some((axis, angle)) = turns[k] else {
            continue;
        };
        let (d, next) = (legs[k].heading_out, legs[k + 1].heading_in);
        let outward = unit(d * next.dot(d) - next);
        let (mut outside, mut inside) = (0.0_f64, 0.0_f64);
        for p in &hull {
            let alpha = (carried_end[k].apply(*p) - axis.location).dot(outward);
            outside = outside.max(alpha);
            inside = inside.max(-alpha);
        }
        // A point inside the turn crosses into the other leg's wall on the
        // legs themselves. Two lines are mirror images in the plane halving
        // the corner and meet on it; against an arc the walls meet where
        // the two paths of each point cross, which they do only where the
        // corner turns in the arc's plane.
        if inside > tol.confusion() && !(straight(k) && straight(k + 1)) {
            let up = axis.direction.vector();
            let arcs = [k, k + 1].map(|j| match legs[j].along {
                Along::Arc(arc, _) => Some(arc),
                Along::Line(_) => None,
            });
            if arcs
                .iter()
                .flatten()
                .any(|arc| arc.direction.vector().cross(up).magnitude() > tol.angular())
            {
                ogeom_bail!(
                    Construction,
                    "an open wire reaches inside a corner that turns an arc out of \
                     its plane: a point's paths down the two legs pass each other \
                     without meeting, so the legs' walls leave a gap no sheet \
                     closes"
                );
            }
            inside_turns[k] = Some(InsideTurn::new(
                legs[k].end,
                up,
                (d, next),
                arcs.map(|arc| arc.map(|arc| arc.location)),
            ));
        }
        // A round corner is told from its square end by the arc its
        // section's farthest point outside the turn sweeps, and by the
        // depth the plane halving the corner cuts its farthest point inside
        // the turn short. A section within the confusion distance of the
        // axis on a side has nothing on that side; past it, an arc or a
        // depth within that distance is no curve or cut the kernel holds.
        let held = |reach: f64, run: f64| reach <= tol.confusion() || run > tol.confusion();
        if rounded[k]
            && (!held(outside, outside * angle) || !held(inside, inside * (angle / 2.0).tan()))
        {
            ogeom_bail!(
                Construction,
                "a round corner turns so little that its round or its cut lies \
                 within the confusion distance of its square end; a mitred \
                 corner sweeps it"
            );
        }
    }

    // Each corner's plane halving it, its normal heading on down the spine:
    // where each leg starts and where it ends.
    let cuts: Vec<(Cut, Cut)> = (0..legs.len())
        .map(|i| {
            (
                (i > 0).then(|| {
                    (
                        legs[i].start,
                        unit(legs[i - 1].heading_out + legs[i].heading_in),
                    )
                }),
                (i < last).then(|| {
                    (
                        legs[i].end,
                        unit(legs[i].heading_out + legs[i + 1].heading_in),
                    )
                }),
            )
        })
        .collect();
    let sweep = Sweep {
        legs: &legs,
        carried: &carried,
        carried_end: &carried_end,
        turns: &turns,
        cuts: &cuts,
        rounded: &rounded,
        inside: &inside_turns,
        width,
        tolerance,
        tol,
    };

    let mut faces: Vec<Shape> = Vec::new();
    if inside_turns.iter().all(Option::is_none) {
        for i in 0..legs.len() {
            faces.extend(sweep.leg_faces(model, profile, i, [None, None])?);
        }
    } else {
        // Each edge of the profile, cut where it crosses a turn's axis,
        // sweeps on its own: inside a corner against an arc its two legs'
        // walls are cut on the one curve where they meet.
        for piece in sweep.profile_pieces(model, profile)? {
            let mut meetings: Vec<Option<(Shape, Point)>> = Vec::with_capacity(last);
            for k in 0..last {
                meetings.push(sweep.meeting_edge(model, &piece, k)?);
            }
            for i in 0..legs.len() {
                let trims = [
                    i.checked_sub(1).and_then(|k| meetings[k].clone()),
                    meetings.get(i).cloned().flatten(),
                ];
                faces.extend(sweep.leg_faces(model, &piece.wire, i, trims)?);
            }
        }
    }

    // The sheet's boundary: the profile where the spine starts and its image
    // where it ends, and the paths of the profile's two ends.
    let ends = wire_ends(model, profile, tol)?.map(|(at, _)| at);
    let mut boundary = Boundary {
        sections: [
            (ends[0], first.heading_in),
            (carried_end[last].apply(ends[0]), legs[last].heading_out),
        ],
        rails: Vec::new(),
    };
    for end in ends {
        for (i, leg) in legs.iter().enumerate() {
            let (at, there) = (carried[i].apply(end), carried_end[i].apply(end));
            let outside = match leg.along {
                Along::Line(v) => {
                    let (d, length) = (unit(v), v.magnitude());
                    let (enters, leaves) = sweep.path(i, v, end)?;
                    boundary
                        .rails
                        .push(Rail::Line(at + d * enters, at + d * leaves));
                    leaves >= length - tol.confusion()
                }
                Along::Arc(axis, _) => {
                    let (enter, leave) = sweep.arc_span(i, end)?;
                    boundary.rails.push(Rail::Arc(
                        axis,
                        Transform::rotation(axis, enter).apply(at),
                        leave - enter,
                    ));
                    let (back, on) = sweep.run_on(i, end);
                    if back > 0.0 {
                        boundary
                            .rails
                            .push(Rail::Line(at - leg.heading_in * back, at));
                    }
                    if on > 0.0 {
                        boundary
                            .rails
                            .push(Rail::Line(there, there + leg.heading_out * on));
                    }
                    i == last || sweep.meeting(i, end)?.is_none()
                }
            };
            if i < last
                && rounded[i]
                && outside
                && let Some((axis, angle)) = turns[i]
            {
                boundary.rails.push(Rail::Arc(axis, there, angle));
            }
        }
    }
    let shape = if closed {
        one_ring(model, &faces, &boundary, tol)?
    } else {
        one_disc(model, &faces, &boundary, tol)?
    };
    let mut history = History::new();
    history.generate(spine, shape.clone());
    history.generate(given, shape.clone());
    Ok(Built::new(shape, history))
}

/// A leg cut in two at its middle: the half from its start, and the half
/// to its end.
fn halved(leg: &Leg) -> (Leg, Leg) {
    let (along, half) = match leg.along {
        Along::Line(v) => (Along::Line(v * 0.5), Transform::translation(v * 0.5)),
        Along::Arc(axis, angle) => (
            Along::Arc(axis, angle * 0.5),
            Transform::rotation(axis, angle * 0.5),
        ),
    };
    let (middle, heading) = (half.apply(leg.start), half.apply_vector(leg.heading_in));
    (
        Leg {
            end: middle,
            heading_out: heading,
            along,
            ..*leg
        },
        Leg {
            start: middle,
            heading_in: heading,
            along,
            ..*leg
        },
    )
}

/// An open wire's sweep down its legs, the corners read: what each leg
/// carries the profile by, how each corner turns and is cut, and the
/// corners against an arc the section reaches inside.
struct Sweep<'a> {
    legs: &'a [Leg],
    /// The motion carrying the profile to each leg's start, and to its
    /// end.
    carried: &'a [Transform],
    carried_end: &'a [Transform],
    turns: &'a [Turn],
    cuts: &'a [(Cut, Cut)],
    rounded: &'a [bool],
    inside: &'a [Option<InsideTurn>],
    /// The profile's extent, the run-on past a mitred end's farthest
    /// crossing.
    width: f64,
    /// What a fitted meeting of a curved piece's walls is held to.
    tolerance: f64,
    tol: Tolerances,
}

/// A piece of the profile swept on its own: a wire of one edge, and the
/// ends of a straight one or a curved one's curve and the range it runs,
/// where the profile stands.
struct ProfilePiece {
    wire: Shape,
    line: Option<[Point; 2]>,
    curve: Option<(Curve, (f64, f64))>,
}

/// A cut on the curve where two legs' walls meet inside a corner: the
/// edge, and the piece's point farthest inside the corner, whose section
/// at the corner the cut takes off.
type Trim = Option<(Shape, Point)>;

impl Sweep<'_> {
    /// Where the profile point `p`'s paths down the legs either side of
    /// corner `k` meet, where it lies inside a corner against an arc.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// where the two paths miss each other.
    fn meeting(&self, k: usize, p: Point) -> OgeomResult<Option<Point>> {
        let Some(Some(turn)) = self.inside.get(k) else {
            return Ok(None);
        };
        let at = self.carried_end[k].apply(p);
        if turn.offset_and_height(at).0 >= -self.tol.confusion() {
            return Ok(None);
        }
        match turn.meeting(at) {
            Some(x) => Ok(Some(x)),
            None => ogeom_bail!(
                Construction,
                "inside a corner against an arc, a point's paths down the two \
                 legs miss each other"
            ),
        }
    }

    /// Where a profile point's path down line leg `i`, from the leg's start,
    /// leaves the start's cut and meets the end's: on the plane halving a
    /// corner, linear in the point, or where it meets the other leg's path
    /// inside a corner against an arc. A rounded corner stops the path of a
    /// point outside the turn square at the corner.
    fn path(&self, i: usize, v: Vector, p: Point) -> OgeomResult<(f64, f64)> {
        let last = self.legs.len() - 1;
        let (d, length) = (v / v.magnitude(), v.magnitude());
        let from = self.carried[i].apply(p);
        let r = from - self.legs[i].start;
        let (start_cut, end_cut) = self.cuts[i];
        let past = |at: f64, lo: f64, hi: f64| -> OgeomResult<f64> {
            if at < lo - self.tol.confusion() || at > hi + self.tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "inside a corner against an arc, a point's paths down the \
                     two legs meet past the line leg's end"
                );
            }
            Ok(at)
        };
        let enters = match i.checked_sub(1).map(|k| self.meeting(k, p)).transpose()? {
            Some(Some(x)) => past((x - from).dot(d), 0.0, length)?,
            _ => {
                let enters = start_cut.map_or(0.0, |(_, n)| -r.dot(n) / d.dot(n));
                if i > 0 && self.rounded[i - 1] {
                    enters.max(0.0)
                } else {
                    enters
                }
            }
        };
        let leaves = match (i < last).then(|| self.meeting(i, p)).transpose()? {
            Some(Some(x)) => past((x - from).dot(d), 0.0, length)?,
            _ => {
                let leaves = end_cut.map_or(length, |(_, n)| length - r.dot(n) / d.dot(n));
                if i < last && self.rounded[i] {
                    leaves.min(length)
                } else {
                    leaves
                }
            }
        };
        Ok((enters, leaves))
    }

    /// The angles about arc leg `i`'s axis between which a profile point's
    /// path runs: the whole arc, or from and to where it meets the other
    /// leg's path inside a corner against an arc.
    fn arc_span(&self, i: usize, p: Point) -> OgeomResult<(f64, f64)> {
        let Along::Arc(axis, angle) = self.legs[i].along else {
            ogeom_bail!(Construction, "the leg is not an arc");
        };
        let (c, a) = (axis.location, axis.direction.vector());
        // The angle about the axis from `from` to `x`, in `(-π, π]`.
        let turned = |from: Point, x: Point| {
            let (f, r) = (from - c, x - c);
            let (pf, pr) = (f - a * f.dot(a), r - a * r.dot(a));
            pf.cross(pr).dot(a).atan2(pf.dot(pr))
        };
        let last = self.legs.len() - 1;
        let enter = match i.checked_sub(1).map(|k| self.meeting(k, p)).transpose()? {
            Some(Some(x)) => turned(self.carried[i].apply(p), x),
            _ => 0.0,
        };
        let leave = match (i < last).then(|| self.meeting(i, p)).transpose()? {
            Some(Some(x)) => angle + turned(self.carried_end[i].apply(p), x),
            _ => angle,
        };
        let reach = {
            let f = self.carried[i].apply(p) - c;
            (f - a * f.dot(a)).magnitude()
        };
        let slack = self.tol.confusion() / reach.max(self.tol.confusion());
        let on_arc = |at: f64| at >= -slack && at <= angle + slack;
        if !on_arc(enter) || !on_arc(leave) {
            ogeom_bail!(
                Construction,
                "inside a corner against an arc, a point's paths down the two \
                 legs meet past the arc's end"
            );
        }
        Ok((enter, leave))
    }

    /// How far a profile point's path runs on straight past arc leg `i`'s
    /// ends to its mitred corners' planes: back from its start, and on from
    /// its end. Negative for a point inside the corner, which runs on
    /// nowhere.
    fn run_on(&self, i: usize, p: Point) -> (f64, f64) {
        let (start_cut, end_cut) = self.cuts[i];
        let back = match start_cut {
            Some((at, n)) if self.turns[i - 1].is_some() && !self.rounded[i - 1] => {
                (self.carried[i].apply(p) - at).dot(n) / self.legs[i].heading_in.dot(n)
            }
            _ => 0.0,
        };
        let on = match end_cut {
            Some((at, n)) if self.turns[i].is_some() && !self.rounded[i] => {
                -(self.carried_end[i].apply(p) - at).dot(n) / self.legs[i].heading_out.dot(n)
            }
            _ => 0.0,
        };
        (back, on)
    }

    /// The profile's edges, each a wire of its own, cut where it crosses
    /// the axis of a corner against an arc it reaches inside, so each
    /// piece lies on one side of every such corner. Each piece runs the way
    /// the profile runs it.
    fn profile_pieces(&self, model: &mut Model, profile: &Shape) -> OgeomResult<Vec<ProfilePiece>> {
        let tol = self.tol;
        let mut pieces = Vec::new();
        for edge in model.ordered_children_of(profile)? {
            let Some((p, q)) = placed_line(model, &edge, tol)? else {
                let (curve, range) = spine_curve_of(model, &edge)?;
                let curve = curve.transformed(&edge.transform(model.datums())?, tol)?;
                let (from, to) = if edge.orientation() == ogeom_topo::Orientation::Reversed {
                    (range.1, range.0)
                } else {
                    range
                };
                let mut splits: Vec<f64> = Vec::new();
                for (k, turn) in self.inside.iter().enumerate() {
                    let Some(turn) = turn else {
                        continue;
                    };
                    // The axis where the section stands at the corner, read
                    // back onto the profile.
                    let back = self.carried_end[k].inverse()?;
                    let (at, outward) = turn.offset_plane();
                    splits.extend(crossings(
                        &curve,
                        (from, to),
                        back.apply(at),
                        back.apply_vector(outward),
                        tol,
                    )?);
                }
                if splits.is_empty() && !self.reaches_inside(&curve, (from, to))? {
                    let wire = make_wire(model, std::slice::from_ref(&edge), tol)?.shape;
                    pieces.push(ProfilePiece {
                        wire,
                        line: None,
                        curve: Some((curve, (from, to))),
                    });
                    continue;
                }
                splits.sort_by(|a, b| {
                    if from <= to {
                        a.total_cmp(b)
                    } else {
                        b.total_cmp(a)
                    }
                });
                let mut ends = vec![from];
                for t in splits {
                    let x = curve.point_at(t, tol)?;
                    let near = |s: f64| {
                        curve
                            .point_at(s, tol)
                            .map(|y| y.distance(x) <= tol.confusion())
                    };
                    if !near(ends[ends.len() - 1])? && !near(to)? {
                        ends.push(t);
                    }
                }
                ends.push(to);
                for pair in ends.windows(2) {
                    let (curve, (from, to)) = own_piece(&curve, (pair[0], pair[1]), tol)?
                        .unwrap_or_else(|| (curve.clone(), (pair[0], pair[1])));
                    let (lo, hi) = (from.min(to), from.max(to));
                    let piece = ogeom_algo::make_edge(model, curve.clone(), (lo, hi), tol)?.shape;
                    let piece = if from > to { piece.reversed() } else { piece };
                    let wire = make_wire(model, &[piece], tol)?.shape;
                    pieces.push(ProfilePiece {
                        wire,
                        line: None,
                        curve: Some((curve, (from, to))),
                    });
                }
                continue;
            };
            let (a, b) = if edge.orientation() == ogeom_topo::Orientation::Reversed {
                (q, p)
            } else {
                (p, q)
            };
            let mut splits: Vec<f64> = Vec::new();
            for (k, turn) in self.inside.iter().enumerate() {
                let Some(turn) = turn else {
                    continue;
                };
                let u = |x: Point| turn.offset_and_height(self.carried_end[k].apply(x)).0;
                let (ua, ub) = (u(a), u(b));
                if (ua < -tol.confusion() && ub > tol.confusion())
                    || (ua > tol.confusion() && ub < -tol.confusion())
                {
                    splits.push(ua / (ua - ub));
                }
            }
            splits.sort_by(f64::total_cmp);
            let mut points = vec![a];
            for s in splits {
                let x = a + (b - a) * s;
                if points
                    .last()
                    .is_some_and(|y| y.distance(x) > tol.confusion())
                    && x.distance(b) > tol.confusion()
                {
                    points.push(x);
                }
            }
            points.push(b);
            for pair in points.windows(2) {
                let line = ogeom_geom::LineCurve::segment(pair[0], pair[1], tol)?;
                let range = ogeom_geom::Curve3d::domain(&line);
                let edge = ogeom_algo::make_edge(model, line.into(), range, tol)?.shape;
                let wire = make_wire(model, &[edge], tol)?.shape;
                pieces.push(ProfilePiece {
                    wire,
                    line: Some([pair[0], pair[1]]),
                    curve: None,
                });
            }
        }
        Ok(pieces)
    }

    /// Whether any point of the curve between `from` and `to` lies inside a
    /// corner against an arc.
    fn reaches_inside(&self, curve: &Curve, (from, to): (f64, f64)) -> OgeomResult<bool> {
        for (k, turn) in self.inside.iter().enumerate() {
            let Some(turn) = turn else {
                continue;
            };
            for f in fractions(32) {
                let p = curve.point_at(from + (to - from) * f, self.tol)?;
                if turn.offset_and_height(self.carried_end[k].apply(p)).0 < -self.tol.confusion() {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }

    /// The edge where the walls `piece` sweeps down the legs either side of
    /// corner `k` meet, where it lies inside a corner against an arc, with
    /// the piece's point farthest inside it. A straight piece's walls meet
    /// on a conic held exactly (see [`InsideTurn::meeting_curve`]). A
    /// curved piece's meet on a curve with no closed form: each of its
    /// points' crossing is found exactly, the curve fitted through them at
    /// the piece's own parameters and held to the confusion distance where
    /// it can be, and its error, measured between the samples, stated as
    /// the edge's tolerance.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// where a point's paths miss each other, or the fitted meeting strays
    /// past the tolerance.
    fn meeting_edge(&self, model: &mut Model, piece: &ProfilePiece, k: usize) -> OgeomResult<Trim> {
        let tol = self.tol;
        let Some(turn) = &self.inside[k] else {
            return Ok(None);
        };
        let u = |x: Point| turn.offset_and_height(self.carried_end[k].apply(x)).0;
        let Some([a, b]) = piece.line else {
            let Some((curve, (from, to))) = &piece.curve else {
                ogeom_bail!(Construction, "a profile piece has no curve");
            };
            let at = |f: f64| curve.point_at(from + (to - from) * f, tol);
            // The piece lies on one side of the axis: inside it where any
            // point is, the point farthest in the one the cut takes off.
            let mut deepest = (f64::INFINITY, at(0.0)?);
            for f in fractions(32) {
                let p = at(f)?;
                let offset = u(p);
                if offset < deepest.0 {
                    deepest = (offset, p);
                }
            }
            if deepest.0 >= -tol.confusion() {
                return Ok(None);
            }
            let meeting = |f: f64| -> OgeomResult<Point> {
                match turn.meeting(self.carried_end[k].apply(at(f)?)) {
                    Some(x) => Ok(x),
                    None => ogeom_bail!(
                        Construction,
                        "inside a corner against an arc, a point's paths down the \
                         two legs miss each other"
                    ),
                }
            };
            let fitted = ogeom_geom::fit::fit_curve_sampled(
                meeting,
                &fractions(16),
                false,
                3,
                tol.confusion(),
                tol,
            )?;
            let fit: Curve = fitted.curve.into();
            let mut error = fitted.error;
            for f in fractions(256) {
                error = error.max(fit.point_at(f, tol)?.distance(meeting(f)?));
            }
            if error > self.tolerance.max(tol.confusion()) {
                ogeom_bail!(
                    Construction,
                    "the walls a curved profile edge sweeps inside a corner against \
                     an arc meet on a curve no fit holds within {}: it strays {error}",
                    self.tolerance
                );
            }
            let edge = ogeom_algo::make_edge(model, fit, (0.0, 1.0), tol)?.shape;
            if error > tol.confusion() {
                let held = ogeom_core::Tolerance::new(error)?;
                model.widen(&edge, held)?;
                if let Some((v0, v1)) = edge_vertices(model, &edge)? {
                    for v in [&v0, &v1] {
                        model.widen(v, held)?;
                    }
                }
            }
            return Ok(Some((edge, deepest.1)));
        };
        let (ua, ub) = (u(a), u(b));
        if ua.min(ub) >= -tol.confusion() {
            return Ok(None);
        }
        let curve = turn.meeting_curve(
            self.carried_end[k].apply(a),
            self.carried_end[k].apply(b),
            tol,
        )?;
        let range = ogeom_geom::Curve3d::domain(&curve);
        let edge = ogeom_algo::make_edge(model, curve, range, tol)?.shape;
        Ok(Some((edge, if ua < ub { a } else { b })))
    }

    /// The faces `section` sweeps down leg `i`: a line leg the prism of its
    /// section and an arc leg its revolution, run on and cut on the planes
    /// halving mitred corners, with the round of a rounded corner after
    /// it; each end given a trim cut on its curve instead.
    fn leg_faces(
        &self,
        model: &mut Model,
        section: &Shape,
        i: usize,
        trims: [Trim; 2],
    ) -> OgeomResult<Vec<Shape>> {
        let tol = self.tol;
        let leg = self.legs[i];
        let last = self.legs.len() - 1;
        let (start_cut, end_cut) = self.cuts[i];
        let rounded = |k: usize| self.rounded[k];
        // The extremes over the section are taken at its hull points; a
        // curved piece's meeting inside a corner is read along it.
        let hull = if trims.iter().any(Option::is_some) {
            along_points(model, section, tol)?
        } else {
            hull_points(model, section, tol)?
        };
        let mut faces: Vec<Shape> = Vec::new();
        // The piece whose end section a round corner after the leg turns.
        let piece = match leg.along {
            Along::Line(v) => {
                let (d, length) = (v / v.magnitude(), v.magnitude());
                let mut shortest = f64::INFINITY;
                let (mut before, mut after) = (0.0_f64, 0.0_f64);
                for p in &hull {
                    let (enters, leaves) = self.path(i, v, *p)?;
                    shortest = shortest.min(leaves - enters);
                    before = before.max(-enters);
                    after = after.max(leaves - length);
                }
                if shortest <= tol.confusion() {
                    ogeom_bail!(
                        Construction,
                        "a leg is shorter than its corners' reach: the cuts at its ends \
                         cross within the section"
                    );
                }
                // A mitred end runs on a section's width past its farthest
                // crossing, so the cut falls inside the prism. A rounded end
                // stops square at the corner, as does an end trimmed where
                // the legs' walls meet before it.
                let before = if trims[0].is_none() && start_cut.is_some() && !rounded(i - 1) {
                    before + self.width
                } else {
                    0.0
                };
                let after = if trims[1].is_none() && end_cut.is_some() && !rounded(i) {
                    after + self.width
                } else {
                    0.0
                };
                let motion = Transform::translation(-d * before) * self.carried[i];
                let section = section_wire(model, section, &motion, tol)?;
                let mut piece =
                    make_prism(model, &section, d * (length + before + after), tol)?.shape;
                if trims[0].is_none()
                    && let Some((at, n)) = start_cut
                {
                    piece = kept_side(model, &piece, at, n, tol)?;
                }
                if trims[1].is_none()
                    && let Some((at, n)) = end_cut
                {
                    piece = kept_side(model, &piece, at, -n, tol)?;
                }
                piece
            }
            Along::Arc(axis, angle) => {
                let mut narrowest = f64::INFINITY;
                for p in &hull {
                    let (enter, leave) = self.arc_span(i, *p)?;
                    let f = self.carried[i].apply(*p) - axis.location;
                    let a = axis.direction.vector();
                    narrowest = narrowest.min((leave - enter) * (f - a * f.dot(a)).magnitude());
                }
                if trims.iter().any(Option::is_some) && narrowest <= tol.confusion() {
                    ogeom_bail!(
                        Construction,
                        "a leg is shorter than its corners' reach: the cuts at its ends \
                         cross within the section"
                    );
                }
                let section_at = section_wire(model, section, &self.carried[i], tol)?;
                let piece = make_revolution(model, &section_at, axis, angle, tol)?.shape;
                // The run-ons past mitred corners, each the prism of the
                // end section a section's width past its farthest crossing,
                // cut on the corner's plane.
                let (mut back, mut on) = (0.0_f64, 0.0_f64);
                for p in &hull {
                    let (b, o) = self.run_on(i, *p);
                    back = back.max(b);
                    on = on.max(o);
                }
                if back > tol.confusion()
                    && let Some((at, n)) = start_cut
                {
                    let d = leg.heading_in;
                    let motion = Transform::translation(-d * (back + self.width)) * self.carried[i];
                    let section = section_wire(model, section, &motion, tol)?;
                    let run = make_prism(model, &section, d * (back + self.width), tol)?.shape;
                    let run = kept_side(model, &run, at, n, tol)?;
                    faces.extend(explore(model, &run, Filter::OfType(ShapeType::Face))?);
                }
                if on > tol.confusion()
                    && let Some((at, n)) = end_cut
                {
                    let d = leg.heading_out;
                    let section = section_wire(model, section, &self.carried_end[i], tol)?;
                    let run = make_prism(model, &section, d * (on + self.width), tol)?.shape;
                    let run = kept_side(model, &run, at, -n, tol)?;
                    faces.extend(explore(model, &run, Filter::OfType(ShapeType::Face))?);
                }
                piece
            }
        };
        let mut piece = piece;
        for (j, trim) in trims.iter().enumerate() {
            if let Some((edge, deepest)) = trim {
                let corner = if j == 0 {
                    self.carried[i].apply(*deepest)
                } else {
                    self.carried_end[i].apply(*deepest)
                };
                piece = trimmed_short_of(model, &piece, edge, corner, tol)?;
            }
        }
        faces.extend(explore(model, &piece, Filter::OfType(ShapeType::Face))?);
        if i < last
            && rounded(i)
            && let Some((axis, angle)) = self.turns[i]
        {
            // A round turns the piece's own edges, which then stand in the
            // round at a second placement; the sewing bakes every face
            // holding such an edge, and a trimmed wall's fitted trace with
            // it. Turned from a copy, the round alone holds them so.
            let turned = if trims.iter().any(Option::is_some) {
                ogeom_algo::copied(model, &piece)?.shape
            } else {
                piece.clone()
            };
            faces.extend(round_corner(
                model,
                &turned,
                axis,
                angle,
                leg.heading_out,
                tol,
            )?);
        }
        Ok(faces)
    }
}

/// The one face of `piece` split on `edge`, a curve running across it
/// from side to side, and the part of it away from `corner`, a point on
/// the boundary of the part the split takes off and of no other.
fn trimmed_short_of(
    model: &mut Model,
    piece: &Shape,
    edge: &Shape,
    corner: Point,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let [face] = explore(model, piece, Filter::OfType(ShapeType::Face))?
        .try_into()
        .map_err(|_| ogeom_core::ogeom_err!(Construction, "a wall to trim is not one face"))?;
    // Read in place, so the trace given the curve is on the surface the
    // split reads.
    let face = ogeom_algo::baked_where_placed(model, &face, tol)?.shape;
    ring_unbroken(model, &face, tol)?;
    charted_on(model, edge, &face, corner, tol)?;
    let split = ogeom_heal::split_face(
        model,
        &face,
        &face,
        std::slice::from_ref(edge),
        ogeom_heal::split::Projection::OnFace,
        tol,
    )?;
    let parts = explore(model, &split.shape, Filter::OfType(ShapeType::Face))?;
    // The part off is told by `corner` as a vertex of it; where it is none
    // (a curved piece's point farthest in lies inside its section), by
    // `corner` on one of its edges.
    let mut off = Vec::with_capacity(parts.len());
    for part in &parts {
        let mut at_vertex = false;
        for v in explore_unique(model, part, ShapeType::Vertex)? {
            let Some(data) = model.node(&v).and_then(|n| n.data().as_vertex()) else {
                continue;
            };
            let at = v.transform(model.datums())?.apply(data.point);
            at_vertex |= at.distance(corner) <= tol.confusion();
        }
        off.push(at_vertex);
    }
    if !off.iter().any(|o| *o) {
        for (part, off) in parts.iter().zip(off.iter_mut()) {
            for e in explore_unique(model, part, ShapeType::Edge)? {
                if model
                    .node(&e)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| d.degenerate)
                {
                    continue;
                }
                let (curve, (t0, t1)) = spine_curve_of(model, &e)?;
                let curve: Curve = ogeom_geom::TrimmedCurve::new(curve, t0, t1, tol)?.into();
                let curve = curve.transformed(&e.transform(model.datums())?, tol)?;
                *off |= ogeom_algo::project_on_curve(&curve, corner, 64, tol)?.distance
                    <= tol.confusion();
            }
        }
    }
    let kept: Vec<&Shape> = parts
        .iter()
        .zip(&off)
        .filter(|(_, off)| !**off)
        .map(|(part, _)| part)
        .collect();
    let ([part], 2) = (kept.as_slice(), parts.len()) else {
        ogeom_bail!(
            Construction,
            "the curve two legs' walls meet on does not cut a wall in two"
        );
    };
    Ok((*part).clone())
}

/// Move the traces of `face`'s outer ring on a periodic surface by whole
/// periods, each to start where the one before it ends: the same traces,
/// the ring run on unbroken in the chart where it was closed only a period
/// apart, as a split reads its outline.
fn ring_unbroken(model: &mut Model, face: &Shape, tol: Tolerances) -> OgeomResult<()> {
    use ogeom_geom::{Curve2d as _, Surface as _};
    let Some(surface_id) = model
        .node(face)
        .and_then(|n| n.data().as_face())
        .map(|d| d.surface)
    else {
        ogeom_bail!(Dangling, "a wall to trim holds no face");
    };
    let Some(surface) = model.geometry().surface(surface_id).cloned() else {
        ogeom_bail!(Dangling, "a wall's surface is not in this model");
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    let period = (
        if surface.is_periodic_u() {
            u1 - u0
        } else {
            0.0
        },
        if surface.is_periodic_v() {
            v1 - v0
        } else {
            0.0
        },
    );
    if period == (0.0, 0.0) {
        return Ok(());
    }
    let Some(outer) = model.children_of(face)?.first().cloned() else {
        return Ok(());
    };
    let mut previous: Option<ogeom_math::Point2> = None;
    for edge in model.ordered_children_of(&outer)? {
        let Some(EdgeRepr::PCurve {
            curve: id, range, ..
        }) = model
            .node(&edge)
            .and_then(|n| n.data().as_edge())
            .and_then(|d| d.pcurve_for(surface_id, edge.location()))
            .cloned()
        else {
            // A seam's two sides stand a period apart by design.
            return Ok(());
        };
        let Some(trace) = model.geometry().pcurve(id).cloned() else {
            continue;
        };
        let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
        let (a, b) = (trace.point_at(range.0, tol)?, trace.point_at(range.1, tol)?);
        let (start, end) = if reversed { (b, a) } else { (a, b) };
        let mut shift = ogeom_math::Vector2::new(0.0, 0.0);
        if let Some(p) = previous {
            if period.0 > 0.0 {
                shift.x = ((start.x - p.x) / period.0).round() * period.0;
            }
            if period.1 > 0.0 {
                shift.y = ((start.y - p.y) / period.1).round() * period.1;
            }
        }
        if shift.x != 0.0 || shift.y != 0.0 {
            let moved = trace.transformed(&ogeom_math::Transform2::translation(-shift), tol)?;
            if let Some(stored) = model.geometry_mut().pcurve_mut(id) {
                *stored = moved;
            }
        }
        previous = Some(ogeom_math::Point2::new(end.x - shift.x, end.y - shift.y));
    }
    Ok(())
}

/// Give `edge`, a curve lying on `face`'s surface, its trace in the
/// surface's chart at its own parameters, held to the confusion distance,
/// where the surface has no closed-form one: on a periodic chart, read in
/// the period `corner`, a point of the face beside it, stands in. Left to
/// the split's own fit where the face is placed or the trace misses.
fn charted_on(
    model: &mut Model,
    edge: &Shape,
    face: &Shape,
    corner: Point,
    tol: Tolerances,
) -> OgeomResult<()> {
    use ogeom_geom::{Curve2d as _, Surface as _};
    if !face.location().is_identity() || !edge.location().is_identity() {
        return Ok(());
    }
    let Some(surface_id) = model
        .node(face)
        .and_then(|n| n.data().as_face())
        .map(|d| d.surface)
    else {
        ogeom_bail!(Dangling, "a wall to trim holds no face");
    };
    let Some(surface) = model.geometry().surface(surface_id).cloned() else {
        ogeom_bail!(Dangling, "a wall's surface is not in this model");
    };
    let (curve, range) = spine_curve_of(model, edge)?;
    if ogeom_intersect::exact_pcurve_over(&curve, range, &surface, tol).is_some() {
        return Ok(());
    }
    let ((u0, u1), (v0, v1)) = surface.domain();
    let period = (
        if surface.is_periodic_u() {
            u1 - u0
        } else {
            0.0
        },
        if surface.is_periodic_v() {
            v1 - v0
        } else {
            0.0
        },
    );
    let periodic = period.0 > 0.0 || period.1 > 0.0;
    // A chart parameter in the period of `near`.
    let wrapped = |x: f64, near: f64, period: f64| {
        if period > 0.0 {
            near + (x - near + period / 2.0).rem_euclid(period) - period / 2.0
        } else {
            x
        }
    };
    let chart = |p: Point, near: ogeom_math::Point2| -> OgeomResult<ogeom_math::Point2> {
        let uv = ogeom_algo::project_on_surface(&surface, p, 16, tol)?.parameters;
        Ok(ogeom_math::Point2::new(
            wrapped(uv.0, near.x, period.0),
            wrapped(uv.1, near.y, period.1),
        ))
    };
    // The period the face's chart holds the corner in, read off its edges'
    // own traces; or, where the corner is no end of theirs, the period the
    // edge it starts on holds its start in.
    let mut near = None;
    let start = curve.point_at(range.0, tol)?;
    if periodic {
        let mut along = None;
        for e in explore(model, face, Filter::OfType(ShapeType::Edge))? {
            let Some(EdgeRepr::PCurve {
                curve: id, range, ..
            }) = model
                .node(&e)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.pcurve_for(surface_id, e.location()))
                .cloned()
            else {
                continue;
            };
            let Some(trace) = model.geometry().pcurve(id).cloned() else {
                continue;
            };
            for t in [range.0, range.1] {
                let uv = trace.point_at(t, tol)?;
                if surface.point_at(uv.x, uv.y, tol)?.distance(corner) <= tol.confusion() {
                    near = Some(uv);
                }
            }
            if along.is_none()
                && !model
                    .node(&e)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| d.degenerate)
            {
                let (own, _) = spine_curve_of(model, &e)?;
                let foot = ogeom_algo::project_on_curve(&own, start, 64, tol)?;
                if foot.distance <= tol.confusion()
                    && foot.parameter >= range.0.min(range.1) - tol.parametric()
                    && foot.parameter <= range.0.max(range.1) + tol.parametric()
                {
                    along = Some(trace.point_at(foot.parameter, tol)?);
                }
            }
        }
        near = near.or(along);
        if near.is_none() {
            return Ok(());
        }
    }
    // Each point charted in the period of the one before, so the trace runs
    // on across the chart's seam.
    const STEPS: u32 = 256;
    let mut unwrapped: Vec<ogeom_math::Point2> = Vec::with_capacity(STEPS as usize + 1);
    let mut previous = near.unwrap_or(ogeom_math::Point2::new(0.0, 0.0));
    for k in 0..=STEPS {
        let t = range.0 + (range.1 - range.0) * f64::from(k) / f64::from(STEPS);
        previous = chart(curve.point_at(t, tol)?, previous)?;
        unwrapped.push(previous);
    }
    let trace = |t: f64| -> OgeomResult<ogeom_math::Point2> {
        let f = ((t - range.0) / (range.1 - range.0)).clamp(0.0, 1.0) * f64::from(STEPS);
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let k = (f.floor() as usize).min(STEPS as usize - 1);
        #[allow(clippy::cast_precision_loss)]
        let g = f - k as f64;
        let (a, b) = (unwrapped[k], unwrapped[k + 1]);
        chart(
            curve.point_at(t, tol)?,
            ogeom_math::Point2::new(a.x + (b.x - a.x) * g, a.y + (b.y - a.y) * g),
        )
    };
    // Held in the chart to the confusion distance over the surface's
    // fastest speed along the trace.
    let mut speed = 0.0_f64;
    let ts: Vec<f64> = (0..=32)
        .map(|k| range.0 + (range.1 - range.0) * f64::from(k) / 32.0)
        .collect();
    for &t in &ts {
        let uv = trace(t)?;
        let (su, sv) = surface.d1_at(uv.x, uv.y, tol)?;
        speed = speed.max(su.magnitude()).max(sv.magnitude());
    }
    let fitted = ogeom_geom::fit::fit_curve_2d_sampled(
        trace,
        &ts,
        false,
        3,
        tol.confusion() / (4.0 * speed.max(1.0)),
        tol,
    )?;
    let pcurve = ogeom_geom::PlanarCurve::BSpline(fitted.curve);
    let mut worst = 0.0_f64;
    for k in 0..=256 {
        let t = range.0 + (range.1 - range.0) * f64::from(k) / 256.0;
        let uv = pcurve.point_at(t, tol)?;
        worst = worst.max(
            surface
                .point_at(uv.x, uv.y, tol)?
                .distance(curve.point_at(t, tol)?),
        );
    }
    // A fitted edge stands off the surface by up to its own tolerance.
    let held = model
        .tolerance_of(edge)?
        .map_or(tol.confusion(), |t| t.get().max(tol.confusion()));
    if worst > held {
        return Ok(());
    }
    ogeom_algo::attach_pcurve(model, edge, pcurve, surface_id, Location::identity(), range)?;
    if let Some(node) = model.node_mut(edge)
        && let ogeom_topo::NodeData::Edge(data) = node.data_mut()
    {
        data.assert_same_parameter(true);
    }
    Ok(())
}

/// The profile's wire rebuilt under `motion`.
fn section_wire(
    model: &mut Model,
    profile: &Shape,
    motion: &Transform,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let Some(ring) = realized_rings(model, profile, motion, tol)?
        .into_iter()
        .next()
    else {
        ogeom_bail!(Construction, "the profile has no wire");
    };
    Ok(make_wire(model, &ring, tol)?.shape)
}

/// The motion carrying the profile to each leg's start, down the legs
/// before it and turned at each corner by the least rotation between the
/// headings, as the rotation-minimizing frame turns; and each leg's turn
/// at its far corner, about the corner's axis, where it turns.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// the spine all but doubles back at a corner.
fn carried_down(legs: &[Leg], tol: Tolerances) -> OgeomResult<(Vec<Transform>, Vec<Turn>)> {
    let mut carried: Vec<Transform> = Vec::with_capacity(legs.len());
    let mut turns: Vec<Turn> = Vec::with_capacity(legs.len());
    let mut motion = Transform::IDENTITY;
    for (i, leg) in legs.iter().enumerate() {
        carried.push(motion);
        motion = leg.motion() * motion;
        let turn = match legs.get(i + 1) {
            Some(next) => {
                let (d, n) = (leg.heading_out, next.heading_in);
                if d.dot(n) < DOUBLED_BACK {
                    ogeom_bail!(
                        Construction,
                        "the spine all but doubles back at a corner; no corner \
                         of that kind turns it"
                    );
                }
                // The turn is read off the cross product's length, not
                // the dot product's, which rounds a turn of an angle's
                // tolerance to none.
                let axis = d.cross(n);
                let sine = axis.magnitude();
                if sine > tol.angular() {
                    let about = Axis {
                        location: next.start,
                        direction: Direction::unit(axis / sine, tol)?,
                    };
                    let angle = sine.atan2(d.dot(n));
                    motion = Transform::rotation(about, angle) * motion;
                    Some((about, angle))
                } else {
                    None
                }
            }
            None => None,
        };
        turns.push(turn);
    }
    Ok((carried, turns))
}

/// The legs of a spine of lines and arcs where they stand, each edge's
/// placement applied, in the order the spine runs them.
fn placed_legs(model: &Model, spine: &Shape, tol: Tolerances) -> OgeomResult<Vec<Leg>> {
    let edges = spine_edges(model, spine)?;
    let mut legs = Vec::with_capacity(edges.len());
    for edge in &edges {
        let (curve, range) = spine_curve_of(model, edge)?;
        let placed = curve.transformed(&edge.transform(model.datums())?, tol)?;
        let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
        let (t0, t1) = if reversed { (range.1, range.0) } else { range };
        let (a, b) = (placed.point_at(t0, tol)?, placed.point_at(t1, tol)?);
        let sense = if reversed { -1.0 } else { 1.0 };
        let heading = |t: f64| -> OgeomResult<Vector> {
            let d = placed.d1_at(t, tol)? * sense;
            let m = d.magnitude();
            if m <= tol.confusion() {
                ogeom_bail!(Construction, "the spine is degenerate at {t}");
            }
            Ok(d / m)
        };
        let (heading_in, heading_out) = (heading(t0)?, heading(t1)?);
        let basis = match &placed {
            Curve::Trimmed(t) => t.basis().clone(),
            other => other.clone(),
        };
        let along = match basis {
            Curve::Line(_) => Along::Line(b - a),
            Curve::Circle(c) => {
                let frame = c.circle().frame();
                // About the axis the way the walk turns at its start.
                let turn = (a - frame.origin()).cross(heading_in);
                let direction = if turn.dot(frame.z().vector()) > 0.0 {
                    frame.z()
                } else {
                    -frame.z()
                };
                let axis = Axis {
                    location: frame.origin(),
                    direction,
                };
                Along::Arc(axis, (t1 - t0).abs())
            }
            _ => ogeom_bail!(
                Construction,
                "a leg of the spine is neither a line nor an arc"
            ),
        };
        legs.push(Leg {
            start: a,
            end: b,
            heading_in,
            heading_out,
            along,
        });
    }
    Ok(legs)
}

/// Refuse a profile on `plane` that does not stand square to the spine at
/// its start `at`, heading `heading`: square to the heading to an angle's
/// tolerance, and through the start, as [`profile_square_at_start`] asks.
///
/// [`profile_square_at_start`]: super::profile_square_at_start
pub(super) fn plane_square_at_start(
    plane: &Plane,
    at: Point,
    heading: Vector,
    tol: Tolerances,
) -> OgeomResult<()> {
    if plane.normal().vector().cross(heading).magnitude() > tol.angular() {
        ogeom_bail!(
            Construction,
            "the profile leans along its spine; a pipe shell runs square to \
             the start"
        );
    }
    if plane.distance_to(at) > tol.confusion() * 100.0 {
        ogeom_bail!(
            Construction,
            "the profile does not sit at the spine's start"
        );
    }
    Ok(())
}

/// Refuse a profile of lines spanning no plane that does not lie in the
/// plane square to the spine at its start `at`, heading `heading`: each
/// line square to the heading to an angle's tolerance, and through that
/// plane, as [`plane_square_at_start`] asks of a planar profile.
pub(super) fn lines_square_at_start(
    model: &Model,
    profile: &Shape,
    at: Point,
    heading: Vector,
    tol: Tolerances,
) -> OgeomResult<()> {
    for edge in model.ordered_children_of(profile)? {
        let (a, b) = match placed_line(model, &edge, tol)? {
            Some(ends) => ends,
            None => ogeom_bail!(Construction, "a pipe shell sweeps a planar profile"),
        };
        let along = b - a;
        if along.dot(heading).abs() > tol.angular() * along.magnitude() {
            ogeom_bail!(
                Construction,
                "the profile leans along its spine; a pipe shell runs square to \
                 the start"
            );
        }
        if (a - at).dot(heading).abs() > tol.confusion() * 100.0 {
            ogeom_bail!(
                Construction,
                "the profile does not sit at the spine's start"
            );
        }
    }
    Ok(())
}

/// A line edge's ends where it stands; `None` for any other curve.
fn placed_line(
    model: &Model,
    edge: &Shape,
    tol: Tolerances,
) -> OgeomResult<Option<(Point, Point)>> {
    let (curve, range) = spine_curve_of(model, edge)?;
    let basis = match &curve {
        Curve::Trimmed(t) => t.basis().clone(),
        other => other.clone(),
    };
    if !matches!(basis, Curve::Line(_)) {
        return Ok(None);
    }
    let placed = curve.transformed(&edge.transform(model.datums())?, tol)?;
    Ok(Some((
        placed.point_at(range.0, tol)?,
        placed.point_at(range.1, tol)?,
    )))
}

/// Points whose convex hull holds the profile: a line's ends, and the
/// corners of any other edge's bound.
fn hull_points(model: &Model, profile: &Shape, tol: Tolerances) -> OgeomResult<Vec<Point>> {
    let mut out = Vec::new();
    for edge in model.ordered_children_of(profile)? {
        match placed_line(model, &edge, tol)? {
            Some((a, b)) => out.extend([a, b]),
            None => out.extend(shape_bounds(model, &edge, tol)?.corners()),
        }
    }
    Ok(out)
}

/// The piece of a curve from `from` to `to` as a curve of its own: a
/// circle's arc framed anew to start at naught, the same points run the
/// same way over `(0, |to - from|)`, so the walls it sweeps hold it in one
/// period of their chart; a spline's piece cut out exactly, its parameters
/// kept, so a surface it sweeps ends where the piece does. `None` for any
/// other curve.
fn own_piece(
    curve: &Curve,
    (from, to): (f64, f64),
    tol: Tolerances,
) -> OgeomResult<Option<(Curve, (f64, f64))>> {
    let basis = match curve {
        Curve::Trimmed(t) => t.basis().clone(),
        other => other.clone(),
    };
    if let Curve::BSpline(spline) = &basis {
        let piece = spline.segment((from.min(to), from.max(to)), tol)?;
        return Ok(Some((piece.into(), (from, to))));
    }
    let Curve::Circle(circle) = basis else {
        return Ok(None);
    };
    let circle = circle.forward(tol)?.circle();
    let frame = circle.frame();
    let start = curve.point_at(from.min(to), tol)?;
    let framed = ogeom_math::Frame::new(
        frame.origin(),
        frame.z(),
        Direction::new(start - frame.origin(), tol)?,
        tol,
    )?;
    let arc: Curve =
        ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(framed, circle.radius(), tol)?).into();
    let span = (to - from).abs();
    Ok(Some((
        arc,
        if from <= to { (0.0, span) } else { (span, 0.0) },
    )))
}

/// Points along every edge of `profile` where it stands: a line's ends, and
/// a curved edge at even steps of its range.
fn along_points(model: &Model, profile: &Shape, tol: Tolerances) -> OgeomResult<Vec<Point>> {
    let mut out = Vec::new();
    for edge in model.ordered_children_of(profile)? {
        match placed_line(model, &edge, tol)? {
            Some((a, b)) => out.extend([a, b]),
            None => {
                let (curve, (t0, t1)) = spine_curve_of(model, &edge)?;
                let placed = curve.transformed(&edge.transform(model.datums())?, tol)?;
                for f in fractions(32) {
                    out.push(placed.point_at(t0 + (t1 - t0) * f, tol)?);
                }
            }
        }
    }
    Ok(out)
}

/// The cut points of a curve, between `from` and `to`, on the plane through
/// `origin` square to `across`: a line's in closed form, any other curve's
/// where it pierces the plane.
pub(super) fn crossings(
    curve: &Curve,
    (from, to): (f64, f64),
    origin: Point,
    across: Vector,
    tol: Tolerances,
) -> OgeomResult<Vec<f64>> {
    let (lo, hi) = (from.min(to), from.max(to));
    let inside = |t: f64| t > lo + tol.parametric() && t < hi - tol.parametric();
    let basis = match curve {
        Curve::Trimmed(t) => t.basis().clone(),
        other => other.clone(),
    };
    if matches!(basis, Curve::Line(_)) {
        let (a, b) = (curve.point_at(from, tol)?, curve.point_at(to, tol)?);
        let (ua, ub) = ((a - origin).dot(across), (b - origin).dot(across));
        if (ua < -tol.confusion() && ub > tol.confusion())
            || (ua > tol.confusion() && ub < -tol.confusion())
        {
            let t = from + (to - from) * (ua / (ua - ub));
            return Ok(if inside(t) { vec![t] } else { Vec::new() });
        }
        return Ok(Vec::new());
    }
    let plane = Plane::through(origin, Direction::new(across, tol)?);
    // The plane's patch reaches past the curve's bound, so its rim stays
    // clear of every crossing.
    let mut reach = 1.0_f64;
    for f in fractions(16) {
        let p = curve.point_at(from + (to - from) * f, tol)?;
        reach = reach.max(p.distance(origin) * 4.0);
    }
    let surface: SurfaceGeometry =
        PlaneSurface::over(plane, (-reach, reach), (-reach, reach))?.into();
    let met = ogeom_intersect::intersect_curve_surface(
        curve,
        &surface,
        ogeom_intersect::CurveSurfaceOptions {
            gap: tol.confusion(),
            ..ogeom_intersect::CurveSurfaceOptions::default()
        },
        tol,
    )?;
    Ok(met
        .crossings
        .iter()
        .map(|p| p.on_curve)
        .filter(|t| inside(*t))
        .collect())
}

/// `n + 1` even fractions of one, from naught to one.
fn fractions(n: u32) -> Vec<f64> {
    (0..=n).map(|i| f64::from(i) / f64::from(n)).collect()
}

/// The part of `piece` on the side of the plane through `at` square to
/// `toward` that `toward` points into.
fn kept_side(
    model: &mut Model,
    piece: &Shape,
    at: Point,
    toward: Vector,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    // The plane's patch reaches twice as far as anything of the piece, so
    // its rim stays clear of the piece.
    let reach = shape_bounds(model, piece, tol)?
        .corners()
        .into_iter()
        .map(|p| p.distance(at))
        .fold(0.0_f64, f64::max);
    let plane = Plane::through(at, Direction::new(toward, tol)?);
    let face = make_natural_face(
        model,
        SurfaceGeometry::Plane(PlaneSurface::over(
            plane,
            (-reach * 2.0, reach * 2.0),
            (-reach * 2.0, reach * 2.0),
        )?),
    )?
    .shape;
    let side = make_half_space(model, &face, at + toward / toward.magnitude() * reach, tol)?.shape;
    Ok(ogeom_bool::common(model, piece, &side, tol)?.shape)
}

/// The round at a corner: each edge the leg `piece` ends on square to its
/// heading `d` at the corner, turned about `axis` through `angle`. An edge
/// along the axis turns on itself and is left out.
fn round_corner(
    model: &mut Model,
    piece: &Shape,
    axis: Axis,
    angle: f64,
    d: Vector,
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    let c = axis.location;
    let a = axis.direction.vector();
    let off_axis = |p: Point| {
        let r = p - c;
        (r - a * r.dot(a)).magnitude()
    };
    let mut out = Vec::new();
    for edge in explore_unique(model, piece, ShapeType::Edge)? {
        let points = edge_points(model, &edge, tol)?;
        if points
            .iter()
            .any(|p| (*p - c).dot(d).abs() > tol.confusion())
            || points.iter().all(|p| off_axis(*p) <= tol.confusion())
        {
            continue;
        }
        // An edge with both ends on the axis would end where its turned
        // copy does, on the same vertices placed again, and the sewing
        // would take the two for one: it is turned in two halves, each
        // with one end off the axis, which the sewing joins to the whole.
        let halves =
            if off_axis(points[0]) <= tol.confusion() && off_axis(points[2]) <= tol.confusion() {
                let (curve, (t0, t1)) = spine_curve_of(model, &edge)?;
                let placed = curve.transformed(&edge.transform(model.datums())?, tol)?;
                let middle = f64::midpoint(t0, t1);
                vec![
                    ogeom_algo::make_edge(model, placed.clone(), (t0, middle), tol)?.shape,
                    ogeom_algo::make_edge(model, placed, (middle, t1), tol)?.shape,
                ]
            } else {
                vec![edge]
            };
        for half in halves {
            let swept = make_revolution(model, &half, axis, angle, tol)?.shape;
            out.extend(explore(model, &swept, Filter::OfType(ShapeType::Face))?);
        }
    }
    Ok(out)
}

/// An edge's ends and the middle of its range, where it stands.
fn edge_points(model: &Model, edge: &Shape, tol: Tolerances) -> OgeomResult<[Point; 3]> {
    let (curve, (t0, t1)) = spine_curve_of(model, edge)?;
    let placed = curve.transformed(&edge.transform(model.datums())?, tol)?;
    Ok([
        placed.point_at(t0, tol)?,
        placed.point_at(0.5 * (t0 + t1), tol)?,
        placed.point_at(t1, tol)?,
    ])
}

/// A corner's plane halving it where a leg ends: a point on it and its
/// normal, heading on down the spine.
type Cut = Option<(Point, Vector)>;

/// Where a swept sheet may be open: its two end sections, each on the
/// plane through a point square to a heading, and the paths of the
/// profile's two ends.
pub(super) struct Boundary {
    pub(super) sections: [(Point, Vector); 2],
    pub(super) rails: Vec<Rail>,
}

/// A piece of the path of one of the profile's ends: a segment down a leg,
/// an arc turned about a round corner's axis from a point through an
/// angle, or a curve fitted along a skinned leg, with the error it is held
/// to.
pub(super) enum Rail {
    Line(Point, Point),
    Arc(Axis, Point, f64),
    Fitted(Curve, f64),
}

impl Rail {
    /// Whether `x` lies on the rail, to the confusion distance, or a fitted
    /// rail's own error.
    fn holds(&self, x: Point, tol: Tolerances) -> OgeomResult<bool> {
        Ok(match *self {
            Rail::Fitted(ref curve, error) => {
                ogeom_algo::project_on_curve(curve, x, 64, tol)?.distance
                    <= error.max(tol.confusion())
            }
            Rail::Line(a, b) => {
                let ab = b - a;
                let t = ((x - a).dot(ab) / ab.dot(ab)).clamp(0.0, 1.0);
                x.distance(a + ab * t) <= tol.confusion()
            }
            Rail::Arc(axis, from, angle) => {
                let (c, a) = (axis.location, axis.direction.vector());
                let (f, r) = (from - c, x - c);
                let (pf, pr) = (f - a * f.dot(a), r - a * r.dot(a));
                let radius = pf.magnitude();
                if (f.dot(a) - r.dot(a)).abs() > tol.confusion()
                    || (radius - pr.magnitude()).abs() > tol.confusion()
                {
                    return Ok(false);
                }
                if radius <= tol.confusion() {
                    return Ok(true);
                }
                let turned = pf.cross(pr).dot(a).atan2(pf.dot(pr));
                let slack = tol.confusion() / radius;
                turned >= -slack && turned <= angle + slack
            }
        })
    }
}

/// The faces sewn into one sheet that is a disc: one shell facing one way,
/// every edge of one face or two, the edges of one face forming one loop
/// that runs only along `boundary`, Euler characteristic one, no finding
/// but its being open, and no two faces meeting other than along what
/// they share.
pub(super) fn one_disc(
    model: &mut Model,
    faces: &[Shape],
    boundary: &Boundary,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    one_sheet(model, faces, boundary, 1, tol)
}

/// The faces sewn into one sheet that is a ring, as [`one_disc`] asks of a
/// disc but bounded by two loops, Euler characteristic naught: a sheet
/// swept round a closed spine, open only along the paths of the
/// profile's two ends.
fn one_ring(
    model: &mut Model,
    faces: &[Shape],
    boundary: &Boundary,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    one_sheet(model, faces, boundary, 2, tol)
}

/// The faces sewn into one sheet of genus naught bounded by `loops` loops
/// running only along `boundary`; see [`one_disc`].
fn one_sheet(
    model: &mut Model,
    faces: &[Shape],
    boundary: &Boundary,
    loops: usize,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    // Sewn within the confusion distance, as any sewing joins, and with an
    // edge split where a neighbour's vertices meet it: the trim may split
    // one leg's cut where the next leg's runs whole.
    let sewn = sew_within(model, faces, tol.confusion(), tol)?;
    let [shell] = sewn.shells.as_slice() else {
        ogeom_bail!(
            Construction,
            "the swept sheet's pieces do not sew into one sheet"
        );
    };
    if !sewn.edges_walked_one_way(model)?.is_empty() {
        ogeom_bail!(Construction, "the swept sheet's pieces do not face one way");
    }
    type Key = (TShapeId, Location);
    let key = |s: &Shape| (s.node(), s.location().clone());
    let mut uses: FastMap<Key, (Shape, usize)> = FastMap::default();
    let mut face_count = 0_i64;
    for face in explore(model, shell, Filter::OfType(ShapeType::Face))? {
        face_count += 1;
        for edge in explore(model, &face, Filter::OfType(ShapeType::Edge))? {
            if model
                .node(&edge)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|data| data.degenerate)
            {
                continue;
            }
            uses.entry(key(&edge))
                .or_insert_with(|| (edge.clone(), 0))
                .1 += 1;
        }
    }
    if uses.values().any(|(_, n)| *n > 2) {
        ogeom_bail!(
            Construction,
            "the swept sheet has an edge of more than two faces"
        );
    }
    // The boundary: each of its edges along `boundary`.
    let mut around: FastMap<Key, Vec<Key>> = FastMap::default();
    for (edge, _) in uses.values().filter(|(_, n)| *n == 1) {
        let points = edge_points(model, edge, tol)?;
        // A fitted end section stands within its own error of the plane.
        let slack = model
            .tolerance_of(edge)?
            .map_or(tol.confusion(), |t| t.get().max(tol.confusion()));
        // A ring has no end sections: the sheet closes on itself there.
        let sections: &[(Point, Vector)] = if loops == 1 { &boundary.sections } else { &[] };
        let on_section = sections
            .iter()
            .any(|&(at, n)| points.iter().all(|p| (*p - at).dot(n).abs() <= slack));
        let mut on_rails = true;
        for p in &points {
            let mut held = false;
            for rail in &boundary.rails {
                if rail.holds(*p, tol)? {
                    held = true;
                    break;
                }
            }
            on_rails &= held;
        }
        if !on_section && !on_rails {
            ogeom_bail!(
                Construction,
                "the swept sheet is open where its pieces should join"
            );
        }
        let Some((a, b)) = edge_vertices(model, edge)? else {
            ogeom_bail!(Construction, "a swept sheet's edge has no vertices");
        };
        around.entry(key(&a)).or_default().push(key(&b));
        around.entry(key(&b)).or_default().push(key(&a));
    }
    // The boundary's loops, each vertex of it on two of its edges.
    let mut seen = ogeom_core::FastSet::default();
    let mut found = 0_usize;
    for start in around.keys() {
        if seen.contains(start) {
            continue;
        }
        found += 1;
        let mut pending: Vec<&Key> = vec![start];
        while let Some(at) = pending.pop() {
            if seen.insert(at.clone()) {
                pending.extend(around[at].iter());
            }
        }
    }
    let bounded = around.values().all(|next| next.len() == 2) && found == loops;
    let vertex_count = explore_unique(model, shell, ShapeType::Vertex)?.len();
    #[allow(clippy::cast_possible_wrap, reason = "counts of a model's shapes")]
    let euler = vertex_count as i64 - uses.len() as i64 + face_count;
    #[allow(clippy::cast_possible_wrap, reason = "a count of loops")]
    if !bounded || euler != 2 - loops as i64 {
        if loops == 1 {
            ogeom_bail!(
                Construction,
                "the swept sheet's pieces do not join into one sheet with one boundary"
            );
        }
        ogeom_bail!(
            Construction,
            "the swept sheet's pieces do not join into one ring bounded by the paths \
             of the profile's two ends"
        );
    }
    if check(model, shell, tol)?
        .problems
        .iter()
        .any(|problem| problem.kind != ShapeType::Shell)
    {
        ogeom_bail!(Construction, "the swept sheet does not check clean");
    }
    if !check_self_intersection(model, shell, tol)?.is_empty() {
        ogeom_bail!(Construction, "the swept sheet runs into itself");
    }
    Ok(shell.clone())
}
