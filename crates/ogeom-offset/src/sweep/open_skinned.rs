//! An open wire swept as one sheet of fitted skins: down a spine with a
//! leg that is neither a line nor an arc, or under a law whose frames are
//! not the rotation-minimizing ones turned by one fixed angle.
//!
//! Each edge of the profile sweeps a strip along each run of the spine
//! between its twin stations, fitted to the profile point's own path and
//! held to the caller's tolerance, as a face's walls are. Neighbouring
//! strips share the rails between them, and neighbouring runs the row they
//! meet on: the later strip adopts the earlier one's border, so no two fits
//! of one curve stand side by side. Nothing caps the ends.
//!
//! At a sharp corner the sheet follows the open-sheet rules of the exact
//! sweep down lines and arcs (see [`open_lines`](super::open_lines)), each
//! profile point on its own:
//!
//! - mitred, a point's path down each leg runs on until it meets its path
//!   down the other: between two straight legs on the plane halving the
//!   corner, the frame mirroring the two legs' sections in it, and against
//!   a curved leg where the two paths cross, the leg's own curve before
//!   the corner and its tangent run on straight past it;
//! - round, a point outside the turn stops square at the corner and turns
//!   about the corner's axis through its angle onto the next leg, the
//!   section's end row revolved; a point inside meets as a mitred one does.
//!
//! The profile is cut where it crosses a rounded corner's axis, or a
//! corner's against a curved leg, so each piece lies on one side of it and
//! each strip's border is one smooth curve. A point inside a corner whose
//! two paths pass each other without meeting (a corner turning a curved leg
//! out of its plane) has no sheet joining the legs, and is refused.
//!
//! A law other than the rotation-minimizing or Frenet frame is asked for
//! its normal wherever the skin is read. Round a sharp corner its sections
//! on the two legs differ unless its frames are the rotation-minimizing
//! ones turned by one fixed angle, and the two legs' cuts with them; only
//! a face of the corner's plane could close the gap, which a sheet has
//! not, so it is refused there.

use ogeom_algo::{Built, History, make_edge_between, make_vertex};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve, Curve3d as _, Surface as _, Transformable as _};
use ogeom_math::{Axis, Direction, Point, Transform, Vector};
use ogeom_topo::{Filter, Model, Shape, ShapeType, explore};

use super::open_lines::{
    Boundary, Rail, crossings, lines_and_arcs, lines_square_at_start, one_disc, open_wire_legs,
    plane_square_at_start, turned_rotation_minimizing,
};
use super::{
    DOUBLED_BACK, PipeCorners, PipeLaw, Skin, SpineStation, SpineWalk, fractions, frenet_normals,
    guide_curves, law_normals, law_normals_on, law_stations, placed_path, rmf_normals,
    shell_stations, skinned_strip, spine_curve_of, walk_curves,
};

/// The samples each strip starts from along its profile piece.
const ALONG_PIECE: usize = 8;

/// Sweep an open wire down an open spine as one sheet of fitted skins; see
/// the module documentation. The sheet faces the way the prism of the
/// profile along the spine's start heading does.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// the profile is not planar or does not stand square at the spine's
/// start, a corner is extended, the spine all but doubles back, the Frenet
/// frame is asked to mitre a corner, a law sets a corner's two sections
/// other than the rotation-minimizing frame turned by one fixed angle
/// does, a point's two paths inside a corner pass each other, a leg is
/// shorter than its corners' reach, or the pieces do not join into one
/// disc; [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where a
/// skin cannot reach `tolerance`.
#[allow(clippy::too_many_lines)]
pub(super) fn open_wire_skinned(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    corners: PipeCorners,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let given_spine = spine;
    let spine = &placed_path(model, spine, tol)?;
    let stations = shell_stations(model, spine, tol)?;
    let turns = |i: usize, stations: &[SpineStation]| {
        let (a, b) = (stations[i].tangent, stations[i + 1].tangent);
        stations[i].at.distance(stations[i + 1].at) <= tol.confusion()
            && (a.cross(b).magnitude() > tol.angular() || a.dot(b) < 0.0)
    };
    let sharp = (0..stations.len() - 1).any(|i| turns(i, &stations));
    if sharp && corners == PipeCorners::Extended {
        ogeom_bail!(
            Construction,
            "an extended corner squares the outside of a solid's section; an open \
             wire sweeps a sheet, whose corners are mitred or round"
        );
    }
    // The frames the sheet is swept in: the rotation-minimizing or Frenet
    // frame carried between stations, or the law's own, asked wherever the
    // skin is read.
    let (law, stations, normals) = match law {
        PipeLaw::RotationMinimizing => {
            let normals = rmf_normals(&stations);
            (None, stations, normals)
        }
        PipeLaw::Frenet if sharp && corners == PipeCorners::Mitre => ogeom_bail!(
            Construction,
            "a Frenet frame has no direction at a corner; sweep a cornered \
             spine with the rotation-minimizing frame"
        ),
        _ if sharp => {
            // Round a sharp corner the law must carry the section as the
            // rotation-minimizing frame does, turned by one fixed angle;
            // the sheet is then that frame's.
            let normals = match law {
                PipeLaw::Frenet => {
                    let mut normals = Vec::with_capacity(stations.len());
                    let mut from = 0;
                    for k in 0..stations.len() {
                        if k + 1 < stations.len() && !turns(k, &stations) {
                            continue;
                        }
                        normals.extend(frenet_normals(&stations[from..=k], tol)?);
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
            let normals = rmf_normals(&stations);
            (None, stations, normals)
        }
        _ => {
            // Down lines and arcs a law carrying the section as the
            // rotation-minimizing frame turned by one fixed angle sweeps
            // that frame's exact sheet.
            let normals = law_normals(model, &stations, law, tolerance, tol)?;
            if lines_and_arcs(model, spine)? && turned_rotation_minimizing(&stations, &normals, tol)
            {
                return open_wire_legs(model, profile, given_spine, corners, tolerance, tol);
            }
            if matches!(law, PipeLaw::Frenet) {
                (None, stations, normals)
            } else {
                let stations = law_stations(model, spine, stations, law, tolerance, tol)?;
                let normals = law_normals(model, &stations, law, tolerance, tol)?;
                (Some(*law), stations, normals)
            }
        }
    };
    let last = stations.len() - 1;
    let start = stations[0];
    match ogeom_algo::find_plane(model, profile, tol)? {
        Some(plane) => plane_square_at_start(&plane, start.at, start.tangent, tol)?,
        None => lines_square_at_start(model, profile, start.at, start.tangent, tol)?,
    }

    let frames = Frames {
        walk: SpineWalk {
            curves: walk_curves(model, spine)?,
            stations: &stations,
            normals: &normals,
        },
        law: match law {
            Some(law) => Some((law, guide_curves(model, &law)?)),
            None => None,
        },
        reach: tolerance,
        tol,
    };
    // The runs between twin stations, each skinned on its own.
    let twins: Vec<usize> = (0..last)
        .filter(|&i| stations[i].at.distance(stations[i + 1].at) <= tol.confusion())
        .collect();
    let runs: Vec<(usize, usize)> = {
        let mut out = Vec::with_capacity(twins.len() + 1);
        let mut from = 0;
        for &k in &twins {
            out.push((from, k));
            from = k + 1;
        }
        out.push((from, last));
        out
    };
    if runs.iter().any(|&(rs, re)| re <= rs) {
        ogeom_bail!(Construction, "a run of the spine has a single station");
    }
    let straight = |(rs, re): (usize, usize)| -> bool {
        law.is_none()
            && (rs..=re).all(|i| {
                stations[i].tangent.cross(stations[rs].tangent).magnitude() <= tol.angular()
            })
    };
    // The profile in the start frame's coordinates.
    let x0 = normals[0];
    let y0 = start.tangent.cross(x0);
    let flat = |p: Point| -> (f64, f64) { ((p - start.at).dot(x0), (p - start.at).dot(y0)) };

    // Each corner between two runs, as the profile meets it.
    let mut joints: Vec<Joint> = Vec::with_capacity(twins.len());
    for (r, &k) in twins.iter().enumerate() {
        let (d, n) = (stations[k].tangent, stations[k + 1].tangent);
        let turn = if turns(k, &stations) {
            if d.dot(n) < DOUBLED_BACK {
                ogeom_bail!(
                    Construction,
                    "the spine all but doubles back at a corner; no corner \
                     of that kind turns it"
                );
            }
            let axis = d.cross(n);
            let sine = axis.magnitude();
            let outward = d * n.dot(d) - n;
            Some(Turn {
                axis: Axis {
                    location: stations[k].at,
                    direction: Direction::unit(axis / sine, tol)?,
                },
                angle: sine.atan2(d.dot(n)),
                outward: outward / outward.magnitude(),
                straight: straight(runs[r]) && straight(runs[r + 1]),
                rounded: corners == PipeCorners::Round,
            })
        } else {
            None
        };
        joints.push(Joint {
            station: k,
            before: runs[r],
            after: runs[r + 1],
            turn,
        });
    }
    let sheet = Sheet {
        frames: &frames,
        reach: tolerance.max(tol.confusion() * 100.0),
        tol,
    };

    // The profile's pieces: its edges where they stand, each cut where it
    // crosses the axis of a corner it is to lie on one side of.
    let mut pieces: Vec<Piece> = Vec::new();
    for edge in model.ordered_children_of(profile)? {
        let (curve, range) = spine_curve_of(model, &edge)?;
        let curve = curve.transformed(&edge.transform(model.datums())?, tol)?;
        let (from, to) = if edge.orientation() == ogeom_topo::Orientation::Reversed {
            (range.1, range.0)
        } else {
            range
        };
        let mut cuts: Vec<f64> = Vec::new();
        for joint in &joints {
            let Some(turn) = &joint.turn else {
                continue;
            };
            if !turn.squares() {
                continue;
            }
            let across = sheet.offset_gradient(joint, x0, y0);
            cuts.extend(crossings(&curve, (from, to), start.at, across, tol)?);
        }
        cuts.sort_by(|a, b| {
            if from <= to {
                a.total_cmp(b)
            } else {
                b.total_cmp(a)
            }
        });
        let mut ends = vec![from];
        for t in cuts {
            let p = curve.point_at(t, tol)?;
            let near = |s: f64| {
                curve
                    .point_at(s, tol)
                    .map(|q| q.distance(p) <= tol.confusion())
            };
            if !near(ends[ends.len() - 1])? && !near(to)? {
                ends.push(t);
            }
        }
        ends.push(to);
        for pair in ends.windows(2) {
            pieces.push(Piece {
                curve: curve.clone(),
                from: pair[0],
                to: pair[1],
            });
        }
    }
    let count = pieces.len();
    // Where each piece lies against each corner: outside its turn, or not.
    let mut outside: Vec<Vec<bool>> = vec![vec![false; joints.len()]; count];
    for (j, piece) in pieces.iter().enumerate() {
        for (c, joint) in joints.iter().enumerate() {
            let Some(turn) = &joint.turn else {
                continue;
            };
            let (mut most, mut least) = (f64::NEG_INFINITY, f64::INFINITY);
            for f in fractions(32) {
                let alpha = sheet.offset(joint, flat(piece.at(f, tol)?))?;
                most = most.max(alpha);
                least = least.min(alpha);
            }
            if turn.squares() && most > tol.confusion() && least < -tol.confusion() {
                ogeom_bail!(
                    Construction,
                    "a profile piece lies across a corner's axis after the profile \
                     was cut there"
                );
            }
            outside[j][c] = most > tol.confusion();
            if turn.squares() && outside[j][c] {
                // A round's arc, or a run-on to the plane halving the
                // corner, within the confusion distance of the square end
                // is no curve or face the kernel holds.
                let reach = if turn.rounded {
                    most * turn.angle
                } else {
                    most * (turn.angle / 2.0).tan()
                };
                if reach <= tol.confusion() {
                    ogeom_bail!(
                        Construction,
                        "a corner against a curved leg turns so little that the round \
                         or the run-on past it lies within the confusion distance of \
                         its square end"
                    );
                }
            }
        }
    }
    // The chain's vertices: each piece's start, and the last piece's end.
    let mut chain: Vec<(f64, f64)> = Vec::with_capacity(count + 1);
    for piece in &pieces {
        chain.push(flat(piece.at(0.0, tol)?));
    }
    chain.push(flat(pieces[count - 1].at(1.0, tol)?));
    // Where two edges of the spine run on into each other the two runs
    // share the section's row, which the frame must place on one curve
    // from both: placing is affine in the profile's coordinates, so the
    // vertices of a polyline and points along each piece settle it.
    for joint in joints.iter().filter(|joint| joint.turn.is_none()) {
        for piece in &pieces {
            for f in fractions(4) {
                let ab = flat(piece.at(f, tol)?);
                let apart = sheet
                    .end(joint, ab, Side::Before, false)?
                    .at
                    .distance(sheet.end(joint, ab, Side::After, false)?.at);
                if apart > tol.confusion() {
                    ogeom_bail!(
                        Construction,
                        "where two edges of the spine run on into each other the \
                         frame places the section {apart} apart on the two, so their \
                         walls do not meet"
                    );
                }
            }
        }
    }
    // Whether a chain vertex's path stops square at a corner: outside a
    // rounded turn, or a mitred one against a curved leg.
    let squared = |v: usize, c: usize| -> bool {
        let neighbours = [v.checked_sub(1), (v < count).then_some(v)];
        neighbours.iter().flatten().any(|&j| outside[j][c])
            && joints[c].turn.is_some_and(|t| t.squares())
            && sheet
                .offset(&joints[c], chain[v])
                .is_ok_and(|alpha| alpha > tol.confusion())
    };

    // The vertices each run starts and ends on.
    let at_station = |model: &mut Model, i: usize| -> OgeomResult<Vec<Shape>> {
        let mut out = Vec::with_capacity(chain.len());
        #[allow(clippy::cast_precision_loss, reason = "a station index")]
        let s = i as f64;
        for ab in &chain {
            out.push(make_vertex(model, frames.generator(s, run_of(&runs, i), *ab)?).shape);
        }
        Ok(out)
    };
    let mut starts: Vec<Vec<Shape>> = Vec::with_capacity(runs.len());
    let mut ends: Vec<Vec<Shape>> = Vec::with_capacity(runs.len());
    // Where each vertex's run-ons past a mitred corner meet its plane.
    let mut cuts: Vec<Vec<Shape>> = Vec::with_capacity(joints.len());
    starts.push(at_station(model, 0)?);
    for (c, joint) in joints.iter().enumerate() {
        let mut end = Vec::with_capacity(chain.len());
        let mut cut = Vec::with_capacity(chain.len());
        let mut next = Vec::with_capacity(chain.len());
        for (v, ab) in chain.iter().enumerate() {
            let square = squared(v, c);
            let here = sheet.end(joint, *ab, Side::Before, square)?.at;
            let vertex = make_vertex(model, here).shape;
            match joint.turn {
                Some(turn) if square && !turn.rounded => {
                    cut.push(make_vertex(model, sheet.run_on(joint, here)).shape);
                }
                _ => cut.push(vertex.clone()),
            }
            next.push(if square {
                make_vertex(model, sheet.end(joint, *ab, Side::After, true)?.at).shape
            } else {
                vertex.clone()
            });
            end.push(vertex);
        }
        ends.push(end);
        cuts.push(cut);
        starts.push(next);
    }
    ends.push(at_station(model, last)?);

    // The strips, run by run along each piece, each run's strip adopting
    // the row the run before ended it on (or that row turned about a
    // rounded corner), and each piece's strip the rail the piece before
    // ran along.
    let mut faces: Vec<Shape> = Vec::new();
    let mut rails: Vec<Rail> = Vec::new();
    let mut tops: Vec<Option<Shape>> = vec![None; count];
    // Each piece's end row turned onto the next run, where it stops square.
    let mut carried: Vec<Option<Shape>> = vec![None; count];
    for (r, &run) in runs.iter().enumerate() {
        let before = r.checked_sub(1).map(|c| &joints[c]);
        let after = joints.get(r);
        let mut previous_rail: Option<Shape> = None;
        // The rulings of the run-ons past the corner after the run, each
        // shared by the two pieces either side of its vertex.
        let mut rulings: Vec<[Option<Shape>; 2]> = vec![[None, None]; count + 1];
        for (j, piece) in pieces.iter().enumerate() {
            let along = |f: f64| -> OgeomResult<(f64, f64)> { Ok(flat(piece.at(f, tol)?)) };
            let square_at = |joint: Option<&Joint>, c: usize| {
                joint.is_some_and(|joint| joint.turn.is_some_and(|t| t.squares()) && outside[j][c])
            };
            let (square_start, square_end) =
                (square_at(before, r.wrapping_sub(1)), square_at(after, r));
            let column = Column {
                sheet: &sheet,
                run,
                before: before.map(|joint| (joint, square_start)),
                after: after.map(|joint| (joint, square_end)),
            };
            let held = std::cell::RefCell::default();
            let skin = if straight(run) {
                Skin::columns(
                    |f| column.straight(along(f)?),
                    fractions(ALONG_PIECE),
                    false,
                    false,
                    false,
                )?
            } else if column.plain() {
                Skin::swept(
                    |f, s| frames.generator(s, run, along(f)?),
                    fractions(ALONG_PIECE),
                    run,
                    false,
                )?
            } else {
                Skin::swept(
                    |f, s| column.spread(&held, f, along(f)?, s),
                    fractions(ALONG_PIECE),
                    run,
                    false,
                )?
            };
            let bottom = if square_start {
                carried[j].take()
            } else {
                tops[j].clone()
            };
            let strip = skinned_strip(
                model,
                &skin,
                (
                    &starts[r][j],
                    &starts[r][j + 1],
                    &ends[r][j],
                    &ends[r][j + 1],
                ),
                [bottom.as_ref(), None, previous_rail.as_ref(), None],
                stations[run.0].at,
                false,
                tolerance,
                tol,
            )?;
            faces.push(facing_as_prism(model, &strip.face, &skin, tol)?);
            for (v, rail) in [(0, &strip.rail0), (count, &strip.rail1)] {
                if (v == 0 && j == 0) || (v == count && j + 1 == count) {
                    let (curve, _) = spine_curve_of(model, rail)?;
                    let error = model.tolerance_of(rail)?.map_or(0.0, |t| t.get());
                    rails.push(Rail::Fitted(curve, error));
                }
            }
            previous_rail = Some(strip.rail1.clone());
            if square_end
                && let Some(joint) = after
                && let Some(turn) = joint.turn
            {
                let turned = turned_row(
                    model,
                    &strip.top,
                    &turn,
                    [&starts[r + 1][j], &starts[r + 1][j + 1]],
                    tol,
                )?;
                if turn.rounded {
                    // The round: the strip's end row turned about the
                    // corner's axis, from a copy so the sewing alone joins
                    // the two.
                    let copy = ogeom_algo::copied(model, &strip.top)?.shape;
                    let swept =
                        ogeom_algo::make_revolution(model, &copy, turn.axis, turn.angle, tol)?
                            .shape;
                    for face in explore(model, &swept, Filter::OfType(ShapeType::Face))? {
                        faces.push(face);
                    }
                } else {
                    let run_on = RunOn {
                        joint,
                        top: &strip.top,
                        turned: &turned,
                        ends: [&ends[r][j], &ends[r][j + 1]],
                        cuts: [&cuts[r][j], &cuts[r][j + 1]],
                        starts: [&starts[r + 1][j], &starts[r + 1][j + 1]],
                    };
                    let (pair, held) =
                        run_on.faces(model, &sheet, [&rulings[j], &rulings[j + 1]])?;
                    faces.extend(pair);
                    rulings[j] = held[0].clone();
                    rulings[j + 1] = held[1].clone();
                    for (v, ruling) in [(0, &held[0]), (count, &held[1])] {
                        if (v == 0 && j == 0) || (v == count && j + 1 == count) {
                            for edge in ruling.iter().flatten() {
                                let (curve, _) = spine_curve_of(model, edge)?;
                                let error = model.tolerance_of(edge)?.map_or(0.0, |t| t.get());
                                rails.push(Rail::Fitted(curve, error));
                            }
                        }
                    }
                }
                carried[j] = Some(turned);
            }
            tops[j] = Some(strip.top);
        }
    }
    // The rounds' arcs the profile's two ends turn along.
    for (c, joint) in joints.iter().enumerate() {
        if let Some(turn) = joint.turn.filter(|t| t.rounded) {
            for v in [0, count] {
                if squared(v, c) {
                    let from = sheet.end(joint, chain[v], Side::Before, true)?.at;
                    rails.push(Rail::Arc(turn.axis, from, turn.angle));
                }
            }
        }
    }
    let boundary = Boundary {
        sections: [
            (start.at, start.tangent),
            (stations[last].at, stations[last].tangent),
        ],
        rails,
    };
    let shape = one_disc(model, &faces, &boundary, tol)?;
    let mut history = History::new();
    history.generate(given_spine, shape.clone());
    history.generate(profile, shape.clone());
    Ok(Built::new(shape, history))
}

/// The run holding station `i`.
fn run_of(runs: &[(usize, usize)], i: usize) -> (usize, usize) {
    runs.iter()
        .copied()
        .find(|&(rs, re)| (rs..=re).contains(&i))
        .unwrap_or(runs[0])
}

/// `top`, a piece's end row square at a corner, turned about the corner's
/// axis through its angle onto the next run, between the next run's start
/// vertices `between`: the row the next run's strip adopts.
fn turned_row(
    model: &mut Model,
    top: &Shape,
    turn: &Turn,
    between: [&Shape; 2],
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let (curve, range) = spine_curve_of(model, top)?;
    let turned = curve.transformed(&Transform::rotation(turn.axis, turn.angle), tol)?;
    let edge = make_edge_between(model, turned, range, between[0], between[1], tol)?.shape;
    if let Some(held) = model.tolerance_of(top)? {
        model.widen(&edge, held)?;
    }
    Ok(edge)
}

/// A piece of the profile: a curve where it stands, run from `from` to
/// `to`.
struct Piece {
    curve: Curve,
    from: f64,
    to: f64,
}

impl Piece {
    /// The point a fraction `f` of the way along the piece.
    fn at(&self, f: f64, tol: Tolerances) -> OgeomResult<Point> {
        self.curve
            .point_at(self.from + (self.to - self.from) * f, tol)
    }
}

/// A curve of an auxiliary law's guide, and its range.
type GuideCurve = (Curve, (f64, f64));

/// The frames the sheet is swept in.
struct Frames<'a> {
    walk: SpineWalk<'a>,
    /// A law asked for its normal wherever the skin is read, with its
    /// guide's curves; `None` for a frame carried between stations.
    law: Option<(PipeLaw<'a>, Vec<GuideCurve>)>,
    /// How far an auxiliary guide may end short of a station's plane.
    reach: f64,
    tol: Tolerances,
}

impl Frames<'_> {
    /// The point `(a, b)` of the profile's start frame where the section
    /// stands at `s` within `run`, `s` a station index.
    fn generator(&self, s: f64, run: (usize, usize), (a, b): (f64, f64)) -> OgeomResult<Point> {
        let (p, t, x) = self.walk.frame_at(s, run, self.tol)?;
        let x = match &self.law {
            Some((law, guide)) if s.fract() != 0.0 => {
                let here = SpineStation {
                    at: p,
                    tangent: t,
                    ..self.walk.stations[run.0]
                };
                law_normals_on(guide, &[here], law, self.reach, self.tol)?[0]
            }
            _ => x,
        };
        let y = t.cross(x);
        Ok(p + x * a + y * b)
    }
}

/// Where the spine splits between two runs: at a corner it turns, or at a
/// smooth junction of its edges.
struct Joint {
    /// The first of the twin stations.
    station: usize,
    before: (usize, usize),
    after: (usize, usize),
    turn: Option<Turn>,
}

/// How a corner turns.
#[derive(Clone, Copy)]
struct Turn {
    axis: Axis,
    angle: f64,
    /// The unit normal of the leg into the corner, in its plane, away from
    /// the turn.
    outward: Vector,
    /// Whether both runs either side are straight.
    straight: bool,
    rounded: bool,
}

impl Turn {
    /// Whether a point outside the turn stops square at the corner: where
    /// it is rounded, or mitred against a curved leg, whose run-on past
    /// the corner is a face of its own.
    fn squares(&self) -> bool {
        self.rounded || !self.straight
    }
}

/// Which side of a corner a path is read on.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Side {
    Before,
    After,
}

/// Where a profile point's path ends at a corner: the point, and the
/// station index along its run, past the corner's station where it runs
/// on straight.
struct End {
    at: Point,
    s: f64,
}

/// The sheet's corners as the profile meets them.
struct Sheet<'a> {
    frames: &'a Frames<'a>,
    /// How far two paths may miss each other and still be taken to meet.
    reach: f64,
    tol: Tolerances,
}

impl Sheet<'_> {
    /// The profile point `ab`'s offset outside the corner's turn, where
    /// the section stands on the leg into it.
    fn offset(&self, joint: &Joint, ab: (f64, f64)) -> OgeomResult<f64> {
        let Some(turn) = &joint.turn else {
            return Ok(0.0);
        };
        #[allow(clippy::cast_precision_loss, reason = "a station index")]
        let s = joint.station as f64;
        let p = self.frames.generator(s, joint.before, ab)?;
        Ok((p - turn.axis.location).dot(turn.outward))
    }

    /// Where a point `p` of the section square at `joint`'s corner, on the
    /// leg into it, runs on straight to the plane halving the corner.
    fn run_on(&self, joint: &Joint, p: Point) -> Point {
        let stations = self.frames.walk.stations;
        let (d, n) = (
            stations[joint.station].tangent,
            stations[joint.station + 1].tangent,
        );
        let m = d + n;
        p + d * ((stations[joint.station].at - p).dot(m) / d.dot(m))
    }

    /// The gradient of [`Self::offset`] over the profile where it stands
    /// at the spine's start, whose start frame's axes are `x0` and `y0`.
    fn offset_gradient(&self, joint: &Joint, x0: Vector, y0: Vector) -> Vector {
        let Some(turn) = &joint.turn else {
            return Vector::ZERO;
        };
        let station = &self.frames.walk.stations[joint.station];
        let x = self.frames.walk.normals[joint.station];
        let y = station.tangent.cross(x);
        x0 * x.dot(turn.outward) + y0 * y.dot(turn.outward)
    }

    /// Where profile point `ab`'s path down the run on `side` of `joint`
    /// ends: square at the corner where `square`, on the plane halving a
    /// corner between straight runs, and where its paths down the two runs
    /// cross against a curved one.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction)
    /// where the two paths miss each other, or meet past the middle of a
    /// run.
    fn end(&self, joint: &Joint, ab: (f64, f64), side: Side, square: bool) -> OgeomResult<End> {
        let (station, run) = match side {
            Side::Before => (joint.station, joint.before),
            Side::After => (joint.station + 1, joint.after),
        };
        #[allow(clippy::cast_precision_loss, reason = "a station index")]
        let s = station as f64;
        let raw = || -> OgeomResult<End> {
            Ok(End {
                at: self.frames.generator(s, run, ab)?,
                s,
            })
        };
        let Some(turn) = &joint.turn else {
            return raw();
        };
        if square {
            return raw();
        }
        let stations = self.frames.walk.stations;
        if turn.straight {
            // A straight run's path is ruled between its ends; its station
            // index is not read.
            let p = self.frames.generator(s, run, ab)?;
            let t = stations[station].tangent;
            let m = stations[joint.station].tangent + stations[joint.station + 1].tangent;
            return Ok(End {
                at: p + t * ((turn.axis.location - p).dot(m) / t.dot(m)),
                s,
            });
        }
        let join = self
            .frames
            .walk
            .join(joint.before, joint.after, ab, self.tol)?;
        if join.gap > self.reach {
            ogeom_bail!(
                Construction,
                "inside a corner against a curved leg, a point's paths down the \
                 two legs pass each other {} apart without meeting, so the legs' \
                 walls leave a gap no sheet closes",
                join.gap
            );
        }
        #[allow(clippy::cast_precision_loss, reason = "station indices")]
        if join.s1 <= joint.before.0 as f64 + 0.5 || join.s2 >= joint.after.1 as f64 - 0.5 {
            ogeom_bail!(
                Construction,
                "a leg is shorter than its corner's reach; the mitre would \
                 run off its far end"
            );
        }
        Ok(End {
            at: join.at,
            s: match side {
                Side::Before => join.s1,
                Side::After => join.s2,
            },
        })
    }
}

/// One piece's paths down one run, between the corners either side: each
/// given with whether the piece stops square there.
struct Column<'a> {
    sheet: &'a Sheet<'a>,
    run: (usize, usize),
    before: Option<(&'a Joint, bool)>,
    after: Option<(&'a Joint, bool)>,
}

impl Column<'_> {
    /// Where the path of `ab` starts and ends on the run.
    fn ends(&self, ab: (f64, f64)) -> OgeomResult<(Option<End>, Option<End>)> {
        let start = match self.before {
            Some((joint, square)) => Some(self.sheet.end(joint, ab, Side::After, square)?),
            None => None,
        };
        let end = match self.after {
            Some((joint, square)) => Some(self.sheet.end(joint, ab, Side::Before, square)?),
            None => None,
        };
        Ok((start, end))
    }

    /// Whether every path down the run starts and ends on the run's own
    /// stations: no corner either side cuts it short or runs it on.
    fn plain(&self) -> bool {
        let at_station = |side: Option<(&Joint, bool)>| {
            side.is_none_or(|(joint, square)| joint.turn.is_none() || square)
        };
        at_station(self.before) && at_station(self.after)
    }

    /// A straight run's path of `ab`: its two ends, ruled.
    fn straight(&self, ab: (f64, f64)) -> OgeomResult<Vec<Point>> {
        let frames = self.sheet.frames;
        #[allow(clippy::cast_precision_loss, reason = "station indices")]
        let (s0, s1) = (self.run.0 as f64, self.run.1 as f64);
        let (start, end) = self.ends(ab)?;
        let a = match start {
            Some(e) => e.at,
            None => frames.generator(s0, self.run, ab)?,
        };
        let b = match end {
            Some(e) => e.at,
            None => frames.generator(s1, self.run, ab)?,
        };
        let t = frames.walk.stations[self.run.0].tangent;
        if (b - a).dot(t) <= self.sheet.tol.confusion() {
            ogeom_bail!(
                Construction,
                "a leg is shorter than its corners' reach: the cuts at its ends \
                 cross within the section"
            );
        }
        Ok(vec![a, b])
    }

    /// The point of `ab`'s path at `s` stations' worth of the way along
    /// the run: the same fraction of the way from where the path starts to
    /// where it ends, in station index. A path inside a corner meets the
    /// other leg's on the leg itself, so it never leaves the leg's curve,
    /// along which the frame moves smoothly with the station index. Each
    /// path's ends are kept in `held` by `u`.
    fn spread(
        &self,
        held: &std::cell::RefCell<ogeom_core::FastMap<u64, (f64, f64)>>,
        u: f64,
        ab: (f64, f64),
        s: f64,
    ) -> OgeomResult<Point> {
        let known = held.borrow().get(&u.to_bits()).copied();
        let (lo, hi) = match known {
            Some(ends) => ends,
            None => {
                let (start, end) = self.ends(ab)?;
                #[allow(clippy::cast_precision_loss, reason = "station indices")]
                let lo = start.map_or(self.run.0 as f64, |e| e.s);
                #[allow(clippy::cast_precision_loss, reason = "station indices")]
                let hi = end.map_or(self.run.1 as f64, |e| e.s);
                if hi <= lo {
                    ogeom_bail!(
                        Construction,
                        "a leg is shorter than its corners' reach: the cuts at its ends \
                         cross within the section"
                    );
                }
                held.borrow_mut().insert(u.to_bits(), (lo, hi));
                (lo, hi)
            }
        };
        #[allow(clippy::cast_precision_loss, reason = "station indices")]
        let f = (s - self.run.0 as f64) / (self.run.1 - self.run.0) as f64;
        self.sheet
            .frames
            .generator(lo + (hi - lo) * f, self.run, ab)
    }
}

/// `face`, a strip skinned over `skin`'s rows, turned to face as the prism
/// of its profile piece along the spine does: its normal along the
/// piece's heading crossed with the spine's, read at the middle of the
/// rows.
fn facing_as_prism(
    model: &Model,
    face: &Shape,
    skin: &Skin<'_>,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let rows = &skin.rows;
    let (j, i) = (rows.len() / 2, rows[0].len() / 2);
    let across = |a: usize, b: usize| rows[b][i] - rows[a][i];
    let along = rows[j][(i + 1).min(rows[j].len() - 1)] - rows[j][i.saturating_sub(1)];
    let down = across(j.saturating_sub(1), (j + 1).min(rows.len() - 1));
    let want = along.cross(down);
    let at = rows[j][i];
    let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
        ogeom_bail!(Construction, "a strip holds no face");
    };
    let Some(surface) = model.geometry().surface(data.surface) else {
        ogeom_bail!(Dangling, "a strip's surface is not in this model");
    };
    let (u, v) = ogeom_algo::project_on_surface(surface, at, 16, tol)?.parameters;
    let (du, dv) = surface.d1_at(u, v, tol)?;
    let mut normal = du.cross(dv);
    if face.orientation() == ogeom_topo::Orientation::Reversed {
        normal = -normal;
    }
    Ok(if normal.dot(want) * PRISM_FACING >= 0.0 {
        face.clone()
    } else {
        face.reversed()
    })
}

/// The sign of a prism's face normal against the profile's heading crossed
/// with the travel.
const PRISM_FACING: f64 = 1.0;

/// The rulings at a piece's two ends, before the plane halving a corner
/// and after it.
type Rulings = [[Option<Shape>; 2]; 2];

/// The run-ons of a piece outside a mitred corner against a curved leg:
/// its end row square at the corner run on straight down the leg's
/// heading to the plane halving the corner, and back from there down the
/// next leg's heading to the row turned onto it. Each is the ruled surface
/// between the rows: the oblique projection onto a plane is affine, so the
/// cut is the row's own spline with its control points projected, exact,
/// and the two run-ons share it.
struct RunOn<'a> {
    joint: &'a Joint,
    /// The end row on the leg into the corner.
    top: &'a Shape,
    /// That row turned onto the leg out of it.
    turned: &'a Shape,
    /// The vertices the piece's two ends stop square on, run on to and
    /// start the next run from.
    ends: [&'a Shape; 2],
    cuts: [&'a Shape; 2],
    starts: [&'a Shape; 2],
}

impl RunOn<'_> {
    /// The two run-on faces, each facing as the prism of the piece along
    /// its leg, and the rulings at the piece's two ends (before the plane,
    /// then after it), `shared` where a neighbouring piece built them.
    fn faces(
        &self,
        model: &mut Model,
        sheet: &Sheet<'_>,
        shared: [&[Option<Shape>; 2]; 2],
    ) -> OgeomResult<(Vec<Shape>, Rulings)> {
        let tol = sheet.tol;
        let stations = sheet.frames.walk.stations;
        let headings = [
            stations[self.joint.station].tangent,
            stations[self.joint.station + 1].tangent,
        ];
        let spline =
            |model: &Model, edge: &Shape| -> OgeomResult<(ogeom_geom::BSplineCurve, (f64, f64))> {
                match spine_curve_of(model, edge)? {
                    (Curve::BSpline(curve), range) => Ok((curve, range)),
                    _ => ogeom_bail!(Construction, "a skinned row is not a spline"),
                }
            };
        let (row, range) = spline(model, self.top)?;
        let (turned, _) = spline(model, self.turned)?;
        let rows: [Vec<ogeom_math::Weighted<Point>>; 3] = [
            row.control_points().to_vec(),
            row.control_points()
                .iter()
                .map(|w| {
                    ogeom_math::Weighted::new(sheet.run_on(self.joint, w.point()), w.weight, tol)
                })
                .collect::<OgeomResult<_>>()?,
            turned.control_points().to_vec(),
        ];
        let knots = row.knots().clone();
        let cut_curve: Curve =
            ogeom_geom::BSplineCurve::rational(knots.clone(), rows[1].clone())?.into();
        // The cut's ends stand where the row's do, run on.
        for (vertex, t) in self.cuts.iter().zip([range.0, range.1]) {
            let at = cut_curve.point_at(t, tol)?;
            let miss = vertex_point(model, vertex)?.distance(at);
            if miss > tol.confusion() {
                model.widen(vertex, ogeom_core::Tolerance::new(miss + tol.confusion())?)?;
            }
        }
        let cut =
            make_edge_between(model, cut_curve, range, self.cuts[0], self.cuts[1], tol)?.shape;
        let mut faces = Vec::with_capacity(2);
        let mut rulings: Rulings = [[None, None], [None, None]];
        for side in 0..2 {
            let (near, far) = (&rows[side], &rows[side + 1]);
            let mut net = Vec::with_capacity(near.len() * 2);
            for (a, b) in near.iter().zip(far) {
                net.push(*a);
                net.push(*b);
            }
            let surface: ogeom_geom::SurfaceGeometry = ogeom_geom::BSplineSurface::rational(
                knots.clone(),
                ogeom_math::KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1)?,
                ogeom_math::ControlGrid::new(net, near.len(), 2)?,
            )?
            .into();
            let id = model.geometry_mut().add_surface(surface.clone());
            let (from, to) = match side {
                0 => (self.ends, self.cuts),
                _ => (self.cuts, self.starts),
            };
            let (bottom, top) = match side {
                0 => (self.top, &cut),
                _ => (&cut, self.turned),
            };
            let row_line = |v: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
                Ok(ogeom_geom::Line2d::over(
                    ogeom_math::Axis2::new(
                        ogeom_math::Point2::new(0.0, v),
                        ogeom_math::Direction2::X,
                    ),
                    range.0 - 1.0,
                    range.1 + 1.0,
                )?
                .into())
            };
            let column_line = |u: f64| -> OgeomResult<ogeom_geom::PlanarCurve> {
                Ok(ogeom_geom::Line2d::over(
                    ogeom_math::Axis2::new(
                        ogeom_math::Point2::new(u, 0.0),
                        ogeom_math::Direction2::Y,
                    ),
                    -1.0,
                    2.0,
                )?
                .into())
            };
            let identity = ogeom_topo::Location::identity();
            ogeom_algo::attach_pcurve(model, bottom, row_line(0.0)?, id, identity.clone(), range)?;
            ogeom_algo::attach_pcurve(model, top, row_line(1.0)?, id, identity.clone(), range)?;
            let mut sides: [Option<Shape>; 2] = [None, None];
            for (k, u) in [range.0, range.1].into_iter().enumerate() {
                let corner = [
                    near[k * (near.len() - 1)].point(),
                    far[k * (far.len() - 1)].point(),
                ];
                let ruling = if from[k].is_same(to[k]) {
                    // On the corner's axis the row runs on nowhere: the
                    // surface's column collapses to the point.
                    let mut data = ogeom_topo::EdgeData::new();
                    data.degenerate = true;
                    model.add_edge(data, &[from[k].clone(), from[k].clone()])?
                } else if let Some(edge) = shared[k][side].clone() {
                    edge
                } else {
                    let (a, b) = (vertex_point(model, from[k])?, vertex_point(model, to[k])?);
                    let line = ogeom_geom::BSplineCurve::new(
                        ogeom_math::KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1)?,
                        vec![a, b],
                        tol,
                    )?;
                    let edge =
                        make_edge_between(model, line.into(), (0.0, 1.0), from[k], to[k], tol)?
                            .shape;
                    let miss = a.distance(corner[0]).max(b.distance(corner[1]));
                    if miss > tol.confusion() {
                        model.widen(&edge, ogeom_core::Tolerance::new(miss + tol.confusion())?)?;
                    }
                    edge
                };
                ogeom_algo::attach_pcurve(
                    model,
                    &ruling,
                    column_line(u)?,
                    id,
                    identity.clone(),
                    (0.0, 1.0),
                )?;
                let degenerate = model
                    .node(&ruling)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| d.degenerate);
                if !degenerate {
                    sides[k] = Some(ruling.clone());
                }
                rulings[k][side] = Some(ruling);
            }
            let [Some(r0), Some(r1)] = rulings.clone().map(|pair| pair[side].clone()) else {
                ogeom_bail!(Construction, "a run-on has no ruling");
            };
            let wire = ogeom_algo::make_wire(
                model,
                &[bottom.clone(), r1, top.reversed(), r0.reversed()],
                tol,
            )?
            .shape;
            let face = ogeom_algo::make_face_on(model, id, std::slice::from_ref(&wire), tol)?.shape;
            // Facing as the prism of the piece along the leg's heading.
            let mid = f64::midpoint(range.0, range.1);
            let (du, dv) = surface.d1_at(mid, 0.5, tol)?;
            let want = du.cross(headings[side]);
            faces.push(if du.cross(dv).dot(want) * PRISM_FACING >= 0.0 {
                face
            } else {
                face.reversed()
            });
            for k in 0..2 {
                if sides[k].is_none() {
                    rulings[k][side] = None;
                }
            }
        }
        Ok((faces, rulings))
    }
}

/// A vertex's point where it stands.
fn vertex_point(model: &Model, vertex: &Shape) -> OgeomResult<Point> {
    let Some(data) = model.node(vertex).and_then(|n| n.data().as_vertex()) else {
        ogeom_bail!(Construction, "a vertex holds no point");
    };
    Ok(vertex.transform(model.datums())?.apply(data.point))
}
