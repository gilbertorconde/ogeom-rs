//! An open wire swept down an open spine of lines and arcs, as one sheet.
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
//! - where the corner turns an arc out of its plane, a point's two paths
//!   pass each other without meeting, and the walls leave a gap only a
//!   face of the plane halving the corner could close, as a face's sweep
//!   closes it. A point there is refused, as is a curved profile edge
//!   inside a corner against an arc, whose walls meet on a curve with no
//!   closed form.
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

/// Sweep an open wire down an open spine of lines and arcs under a frame
/// law other than the rotation-minimizing one, round sharp corners: as
/// [`open_wire_legs`] where the law's frames are the rotation-minimizing
/// ones turned by one fixed angle, so that it places every section as that
/// frame carries it. The Frenet frame is read run by run between the sharp
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
    if let PipeLaw::Auxiliary { guide } = law {
        let placed = !guide.location().is_identity()
            || spine_edges(model, guide)?
                .iter()
                .any(|edge| !edge.location().is_identity());
        if placed {
            ogeom_bail!(
                Construction,
                "an auxiliary spine placed by a location is not read where it \
                 stands; an open wire's sheet under it is refused"
            );
        }
    }
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
    open_wire_legs(model, profile, spine, corners, tol)
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

/// Sweep an open wire down an open spine of lines and arcs under the
/// rotation-minimizing frame, each corner mitred or rounded; see the
/// module documentation. The sheet faces the way the prism of the profile
/// along the first leg's heading does.
pub(super) fn open_wire_legs(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    corners: PipeCorners,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let legs = placed_legs(model, spine, tol)?;
    let unit = |v: Vector| v / v.magnitude();
    let sharp = |a: Vector, b: Vector| a.cross(b).magnitude() > tol.angular() || a.dot(b) < 0.0;
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
    let Some(first) = legs.first().copied() else {
        ogeom_bail!(Construction, "the spine has no leg");
    };
    match ogeom_algo::find_plane(model, profile, tol)? {
        Some(plane) => plane_square_at_start(&plane, first.start, first.heading_in, tol)?,
        None => lines_square_at_start(model, profile, first.start, first.heading_in, tol)?,
    }
    let (carried, turns) = carried_down(&legs, tol)?;
    // The motion carrying the profile to each leg's end.
    let carried_end: Vec<Transform> = legs
        .iter()
        .zip(&carried)
        .map(|(leg, motion)| leg.motion() * *motion)
        .collect();
    let hull = hull_points(model, profile, tol)?;
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
    let shape = one_disc(model, &faces, &boundary, tol)?;
    let mut history = History::new();
    history.generate(spine, shape.clone());
    history.generate(profile, shape.clone());
    Ok(Built::new(shape, history))
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
    tol: Tolerances,
}

/// A piece of the profile swept on its own: a wire of one edge, and the
/// ends of a straight one, where the profile stands.
struct ProfilePiece {
    wire: Shape,
    line: Option<[Point; 2]>,
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

    /// The profile's edges, each a wire of its own, a straight edge cut
    /// where it crosses the axis of a corner against an arc it reaches
    /// inside, so each piece lies on one side of every such corner. Each
    /// piece runs the way the profile runs it.
    fn profile_pieces(&self, model: &mut Model, profile: &Shape) -> OgeomResult<Vec<ProfilePiece>> {
        let tol = self.tol;
        let mut pieces = Vec::new();
        for edge in model.ordered_children_of(profile)? {
            let Some((p, q)) = placed_line(model, &edge, tol)? else {
                let wire = make_wire(model, std::slice::from_ref(&edge), tol)?.shape;
                pieces.push(ProfilePiece { wire, line: None });
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
                });
            }
        }
        Ok(pieces)
    }

    /// The edge where the walls `piece` sweeps down the legs either side of
    /// corner `k` meet, where it lies inside a corner against an arc, with
    /// the piece's end farthest inside it.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// where a curved piece reaches inside such a corner: its walls meet
    /// on a curve with no closed form.
    fn meeting_edge(&self, model: &mut Model, piece: &ProfilePiece, k: usize) -> OgeomResult<Trim> {
        let tol = self.tol;
        let Some(turn) = &self.inside[k] else {
            return Ok(None);
        };
        let u = |x: Point| turn.offset_and_height(self.carried_end[k].apply(x)).0;
        let Some([a, b]) = piece.line else {
            if hull_points(model, &piece.wire, tol)?
                .iter()
                .any(|p| u(*p) < -tol.confusion())
            {
                ogeom_bail!(
                    Construction,
                    "a curved edge of an open wire reaches inside a corner against \
                     an arc, where its two legs' walls meet on a curve with no \
                     closed form"
                );
            }
            return Ok(None);
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
        // The extremes over the section are taken at its hull points.
        let hull = hull_points(model, section, tol)?;
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
/// from side to side, and the part of it away from `corner`, a vertex of
/// the part the split takes off.
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
    charted_on(model, edge, &face, corner, tol)?;
    let split = ogeom_heal::split_face(
        model,
        &face,
        &face,
        std::slice::from_ref(edge),
        ogeom_heal::split::Projection::OnFace,
        tol,
    )?;
    let mut kept = Vec::new();
    let mut parts = 0;
    for part in explore(model, &split.shape, Filter::OfType(ShapeType::Face))? {
        parts += 1;
        let mut off = false;
        for v in explore_unique(model, &part, ShapeType::Vertex)? {
            let Some(data) = model.node(&v).and_then(|n| n.data().as_vertex()) else {
                continue;
            };
            let at = v.transform(model.datums())?.apply(data.point);
            off |= at.distance(corner) <= tol.confusion();
        }
        if !off {
            kept.push(part);
        }
    }
    let ([part], 2) = (kept.as_slice(), parts) else {
        ogeom_bail!(
            Construction,
            "the curve two legs' walls meet on does not cut a wall in two"
        );
    };
    Ok(part.clone())
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
    let period = if surface.is_periodic_u() {
        let ((u0, u1), _) = surface.domain();
        u1 - u0
    } else {
        0.0
    };
    let chart = |p: Point, near: f64| -> OgeomResult<ogeom_math::Point2> {
        let mut uv = ogeom_algo::project_on_surface(&surface, p, 16, tol)?.parameters;
        if period > 0.0 {
            uv.0 = near + (uv.0 - near + period / 2.0).rem_euclid(period) - period / 2.0;
        }
        Ok(ogeom_math::Point2::new(uv.0, uv.1))
    };
    // The period the face's chart holds the corner in, read off its edges'
    // own traces.
    let mut near = None;
    if period > 0.0 {
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
            let Some(trace) = model.geometry().pcurve(id) else {
                continue;
            };
            for t in [range.0, range.1] {
                let uv = trace.point_at(t, tol)?;
                if surface.point_at(uv.x, uv.y, tol)?.distance(corner) <= tol.confusion() {
                    near = Some(uv.x);
                }
            }
        }
        if near.is_none() {
            return Ok(());
        }
    }
    let near = near.unwrap_or(0.0);
    let first = chart(curve.point_at(range.0, tol)?, near)?;
    let last = chart(curve.point_at(range.1, tol)?, first.x)?;
    let trace = |t: f64| -> OgeomResult<ogeom_math::Point2> {
        let f = (t - range.0) / (range.1 - range.0);
        chart(curve.point_at(t, tol)?, first.x + (last.x - first.x) * f)
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
    if worst > tol.confusion() {
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
fn plane_square_at_start(
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
fn lines_square_at_start(
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
        let swept = make_revolution(model, &edge, axis, angle, tol)?.shape;
        out.extend(explore(model, &swept, Filter::OfType(ShapeType::Face))?);
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
struct Boundary {
    sections: [(Point, Vector); 2],
    rails: Vec<Rail>,
}

/// A piece of the path of one of the profile's ends: a segment down a leg,
/// or an arc turned about a round corner's axis from a point through an
/// angle.
enum Rail {
    Line(Point, Point),
    Arc(Axis, Point, f64),
}

impl Rail {
    /// Whether `x` lies on the rail, to the confusion distance.
    fn holds(&self, x: Point, tol: Tolerances) -> bool {
        match *self {
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
                    return false;
                }
                if radius <= tol.confusion() {
                    return true;
                }
                let turned = pf.cross(pr).dot(a).atan2(pf.dot(pr));
                let slack = tol.confusion() / radius;
                turned >= -slack && turned <= angle + slack
            }
        }
    }
}

/// The faces sewn into one sheet that is a disc: one shell facing one way,
/// every edge of one face or two, the edges of one face forming one loop
/// that runs only along `boundary`, Euler characteristic one, no finding
/// but its being open, and no two faces meeting other than along what
/// they share.
fn one_disc(
    model: &mut Model,
    faces: &[Shape],
    boundary: &Boundary,
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
    // The boundary: each of its vertices on two of its edges, and all of
    // them one loop.
    let mut around: FastMap<Key, Vec<Key>> = FastMap::default();
    for (edge, _) in uses.values().filter(|(_, n)| *n == 1) {
        let points = edge_points(model, edge, tol)?;
        let on_section = boundary.sections.iter().any(|&(at, n)| {
            points
                .iter()
                .all(|p| (*p - at).dot(n).abs() <= tol.confusion())
        });
        let on_rails = points
            .iter()
            .all(|p| boundary.rails.iter().any(|rail| rail.holds(*p, tol)));
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
    let one_loop = !around.is_empty() && around.values().all(|next| next.len() == 2) && {
        let mut seen = ogeom_core::FastSet::default();
        let mut pending: Vec<&Key> = around.keys().take(1).collect();
        while let Some(at) = pending.pop() {
            if seen.insert(at.clone()) {
                pending.extend(around[at].iter());
            }
        }
        seen.len() == around.len()
    };
    let vertex_count = explore_unique(model, shell, ShapeType::Vertex)?.len();
    #[allow(clippy::cast_possible_wrap, reason = "counts of a model's shapes")]
    let euler = vertex_count as i64 - uses.len() as i64 + face_count;
    if !one_loop || euler != 1 {
        ogeom_bail!(
            Construction,
            "the swept sheet's pieces do not join into one sheet with one boundary"
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
