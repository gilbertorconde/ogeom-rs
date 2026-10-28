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
//! A point of the section inside a turn crosses the plane halving it on
//! the legs themselves, not on their run-ons. Two lines are mirror images
//! there too, so their paths still meet on the plane; a line and an arc,
//! or two arcs, are not, and a point inside a corner against an arc is
//! refused.
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
use ogeom_topo::{Filter, Location, Model, Shape, ShapeType, TShapeId, explore, explore_unique};

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
    // Whether the corner after leg `k` is rounded: a corner that turns
    // nothing has no round, and is cut on the plane square to both legs.
    let rounded = |k: usize| corners == PipeCorners::Round && turns[k].is_some();
    let last = legs.len() - 1;
    // How far the section reaches outside and inside each turning corner:
    // its farthest points from the corner's axis on either side, in the
    // plane square to the leg there.
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
        // A point inside the turn crosses the plane halving it on the legs
        // themselves; only two lines meet it there at one point.
        if inside > tol.confusion() && !(straight(k) && straight(k + 1)) {
            ogeom_bail!(
                Construction,
                "an open wire reaches inside a corner against an arc: its \
                 points' paths down the two legs cross the plane halving the \
                 corner at different points, so the legs' cuts do not join"
            );
        }
        // A round corner is told from its square end by the arc its
        // section's farthest point outside the turn sweeps, and by the
        // depth the plane halving the corner cuts its farthest point inside
        // the turn short. A section within the confusion distance of the
        // axis on a side has nothing on that side; past it, an arc or a
        // depth within that distance is no curve or cut the kernel holds.
        let held = |reach: f64, run: f64| reach <= tol.confusion() || run > tol.confusion();
        if rounded(k)
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
    // Where a profile point's path down line leg `i`, from the leg's start,
    // leaves the start's cut and meets the end's, both linear in the point.
    // A rounded corner stops the path of a point outside the turn square at
    // the corner.
    let path = |i: usize, v: Vector, p: Point| -> (f64, f64) {
        let (d, length) = (unit(v), v.magnitude());
        let r = carried[i].apply(p) - legs[i].start;
        let (start_cut, end_cut) = cuts[i];
        let enters = start_cut.map_or(0.0, |(_, n)| -r.dot(n) / d.dot(n));
        let leaves = end_cut.map_or(length, |(_, n)| length - r.dot(n) / d.dot(n));
        let enters = if i > 0 && rounded(i - 1) {
            enters.max(0.0)
        } else {
            enters
        };
        let leaves = if i < last && rounded(i) {
            leaves.min(length)
        } else {
            leaves
        };
        (enters, leaves)
    };
    // How far a profile point's path runs on straight past arc leg `i`'s
    // ends to its mitred corners' planes: back from its start, and on from
    // its end. No point lies inside such a corner, so neither is negative.
    let run_on = |i: usize, p: Point| -> (f64, f64) {
        let (start_cut, end_cut) = cuts[i];
        let back = match start_cut {
            Some((at, n)) if turns[i - 1].is_some() && !rounded(i - 1) => {
                (carried[i].apply(p) - at).dot(n) / legs[i].heading_in.dot(n)
            }
            _ => 0.0,
        };
        let on = match end_cut {
            Some((at, n)) if turns[i].is_some() && !rounded(i) => {
                -(carried_end[i].apply(p) - at).dot(n) / legs[i].heading_out.dot(n)
            }
            _ => 0.0,
        };
        (back, on)
    };

    let mut faces: Vec<Shape> = Vec::new();
    for (i, leg) in legs.iter().enumerate() {
        let (start_cut, end_cut) = cuts[i];
        // The piece whose end section a round corner after the leg turns.
        let piece = match leg.along {
            Along::Line(v) => {
                let (d, length) = (unit(v), v.magnitude());
                // The extremes of the crossings over the profile are taken
                // at the hull points.
                let mut shortest = f64::INFINITY;
                let (mut before, mut after) = (0.0_f64, 0.0_f64);
                for p in &hull {
                    let (enters, leaves) = path(i, v, *p);
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
                // stops square at the corner.
                let before = if start_cut.is_some() && !rounded(i - 1) {
                    before + width
                } else {
                    0.0
                };
                let after = if end_cut.is_some() && !rounded(i) {
                    after + width
                } else {
                    0.0
                };
                let motion = Transform::translation(-d * before) * carried[i];
                let section = section_wire(model, profile, &motion, tol)?;
                let mut piece =
                    make_prism(model, &section, d * (length + before + after), tol)?.shape;
                if let Some((at, n)) = start_cut {
                    piece = kept_side(model, &piece, at, n, tol)?;
                }
                if let Some((at, n)) = end_cut {
                    piece = kept_side(model, &piece, at, -n, tol)?;
                }
                piece
            }
            Along::Arc(axis, angle) => {
                let section = section_wire(model, profile, &carried[i], tol)?;
                let piece = make_revolution(model, &section, axis, angle, tol)?.shape;
                // The run-ons past mitred corners, each the prism of the
                // end section a section's width past its farthest crossing,
                // cut on the corner's plane.
                let (mut back, mut on) = (0.0_f64, 0.0_f64);
                for p in &hull {
                    let (b, o) = run_on(i, *p);
                    back = back.max(b);
                    on = on.max(o);
                }
                if back > tol.confusion()
                    && let Some((at, n)) = start_cut
                {
                    let d = leg.heading_in;
                    let motion = Transform::translation(-d * (back + width)) * carried[i];
                    let section = section_wire(model, profile, &motion, tol)?;
                    let run = make_prism(model, &section, d * (back + width), tol)?.shape;
                    let run = kept_side(model, &run, at, n, tol)?;
                    faces.extend(explore(model, &run, Filter::OfType(ShapeType::Face))?);
                }
                if on > tol.confusion()
                    && let Some((at, n)) = end_cut
                {
                    let d = leg.heading_out;
                    let section = section_wire(model, profile, &carried_end[i], tol)?;
                    let run = make_prism(model, &section, d * (on + width), tol)?.shape;
                    let run = kept_side(model, &run, at, -n, tol)?;
                    faces.extend(explore(model, &run, Filter::OfType(ShapeType::Face))?);
                }
                piece
            }
        };
        faces.extend(explore(model, &piece, Filter::OfType(ShapeType::Face))?);
        if i < last
            && rounded(i)
            && let Some((axis, angle)) = turns[i]
        {
            faces.extend(round_corner(
                model,
                &piece,
                axis,
                angle,
                leg.heading_out,
                tol,
            )?);
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
                    let (enters, leaves) = path(i, v, end);
                    boundary
                        .rails
                        .push(Rail::Line(at + d * enters, at + d * leaves));
                    leaves >= length - tol.confusion()
                }
                Along::Arc(axis, angle) => {
                    boundary.rails.push(Rail::Arc(axis, at, angle));
                    let (back, on) = run_on(i, end);
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
                    true
                }
            };
            if i < last
                && rounded(i)
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
