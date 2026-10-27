//! Pipe shells round sharp corners, run by run.
//!
//! The spine is cut at its sharp corners into runs, each smooth. Every run
//! is swept on its own with the section it starts with: a straight run
//! whose section keeps its frame along it is the prism of that section,
//! any other run the pipe shell along it. The runs then meet at each corner
//! the way the corner asks:
//!
//! - mitred, each run is run on straight past the corner by the prism of
//!   its section there and trimmed by the plane halving the corner, so the
//!   two meet on that plane whatever their sections are; a straight run
//!   meeting a skinned one there takes the skin's own fit of each ring the
//!   two share, so their cuts of the plane are one;
//! - extended, each run is run on straight past the corner until its
//!   section clears the other run's far side;
//! - round, the half of the section outside the corner is turned about the
//!   corner's axis through the corner's angle, which joins the two runs
//!   only where the next run's section is this one's turned so.
//!
//! The rotation-minimizing law carries the section across a corner by
//! that turn: the next run starts with the section the last one ended
//! with, turned about the corner. The Frenet, binormal and auxiliary laws
//! place the section at every point of the spine themselves; the Frenet
//! frame of a run is its own curvature's, so a straight run has none. The pieces are prisms,
//! revolutions, the pipe shells of smooth runs and the planes trimming
//! them, so the solid is exact wherever the runs are swept exactly. A
//! closed wire sweeps the walls of the face it bounds.

use ogeom_algo::{
    Built, History, curve_length, make_box, make_prism, make_revolution, parameter_at_length,
    project_on_curve, shape_bounds,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve3d as _, SurfaceGeometry};
use ogeom_math::{Axis, Direction, Frame, Point, Transform, Vector};
use ogeom_topo::{Filter, Model, NodeData, Shape, ShapeType, explore};

use super::evolving::{Evolution, evolved_run};
use super::{
    DOUBLED_BACK, LEAST_TURN_SINE, PipeCorners, PipeLaw, bounded_face, fuse_in_order, law_normals,
    make_pipe_shell_with, outside_half, profile_square_at_start, realized_profile_wound,
    rmf_normals, sample_wire, shell_stations, spine_curve_of, station_frame, walls_of,
};

/// A spine edge's ends and its unit tangents there, in the direction of
/// travel.
struct Leg {
    edge: Shape,
    straight: bool,
    start: Point,
    end: Point,
    heading_in: Vector,
    heading_out: Vector,
}

/// A run of spine edges between sharp corners.
struct Run {
    legs: std::ops::Range<usize>,
}

/// Samples round a section's outer ring for the reaches past a corner, as
/// many as the straight-legged corners take.
const AROUND: usize = 256;

/// Whether two unit headings meet at a sharp corner: turned by more than
/// an angle's tolerance, or doubled back.
fn sharp(a: Vector, b: Vector, tol: Tolerances) -> bool {
    a.cross(b).magnitude() > tol.angular() || a.dot(b) < 0.0
}

/// Sweep a face along a spine with sharp corners (or round a closed one)
/// under any frame law but the fixed one, each corner turned as `corners`
/// asks; see the module documentation. Where `evolution` is given the
/// section also scales and twists by it: each run by the laws over its
/// share of the spine's length, from the section the run starts with.
#[allow(clippy::too_many_arguments)]
pub(super) fn law_cornered(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    corners: PipeCorners,
    evolution: Option<&Evolution>,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if model.kind_of(profile)? == ShapeType::Wire {
        // A closed planar wire sweeps the walls of the face it bounds: that
        // face is swept and its end caps, on the planes square to the
        // spine's ends, taken off.
        if !ogeom_algo::is_wire_closed(model, profile, tol)? {
            ogeom_bail!(
                Construction,
                "a pipe turned round a sharp corner run by run sweeps a face or \
                 a closed wire; an open wire bounds no face to join the runs with"
            );
        }
        let face = bounded_face(model, profile, tol)?;
        let solid =
            law_cornered(model, &face, spine, law, corners, evolution, tolerance, tol)?.shape;
        let legs = legs_of(model, spine, tol)?;
        let (first, last) = (&legs[0], &legs[legs.len() - 1]);
        let ends = if first.start.distance(last.end) <= tol.confusion() * 10.0 {
            Vec::new()
        } else {
            vec![
                (first.start, first.heading_in),
                (last.end, last.heading_out),
            ]
        };
        let shape = walls_of(model, &solid, &ends, tol)?;
        let mut history = History::new();
        history.generate(spine, shape.clone());
        history.generate(profile, shape.clone());
        return Ok(Built::new(shape, history));
    }
    if model.kind_of(profile)? != ShapeType::Face {
        ogeom_bail!(
            Construction,
            "a pipe turned round a sharp corner run by run sweeps a face or a closed wire"
        );
    }
    let Some(plane) = ogeom_algo::find_plane(model, profile, tol)? else {
        ogeom_bail!(Construction, "a pipe shell sweeps a planar profile");
    };
    let stations = shell_stations(model, spine, tol)?;
    profile_square_at_start(model, &plane, spine, stations[0].at, tol)?;

    let legs = legs_of(model, spine, tol)?;
    let n = legs.len();
    let closed = legs[0].start.distance(legs[n - 1].end) <= tol.confusion() * 10.0;
    // A run starts at the spine's start and after every sharp corner.
    let mut starts: Vec<usize> = vec![0];
    for k in 1..n {
        if sharp(legs[k - 1].heading_out, legs[k].heading_in, tol) {
            starts.push(k);
        }
    }
    let runs: Vec<Run> = (0..starts.len())
        .map(|r| Run {
            legs: starts[r]..starts.get(r + 1).copied().unwrap_or(n),
        })
        .collect();
    if matches!(law, PipeLaw::Frenet)
        && runs
            .iter()
            .any(|run| legs[run.legs.clone()].iter().all(|l| l.straight))
    {
        ogeom_bail!(
            Construction,
            "the Frenet frame is not defined along a straight run between sharp corners"
        );
    }
    let wrap_sharp = closed && sharp(legs[n - 1].heading_out, legs[0].heading_in, tol);
    if runs.len() == 1 && !wrap_sharp {
        ogeom_bail!(
            Construction,
            "the spine turns no sharp corner; a pipe shell sweeps it whole"
        );
    }

    // The share of the spine's length run to each leg's start, and to its
    // end.
    let mut run_to = Vec::with_capacity(n + 1);
    run_to.push(0.0);
    for leg in &legs {
        let (curve, range) = spine_curve_of(model, &leg.edge)?;
        run_to.push(run_to[run_to.len() - 1] + curve_length(&curve, range, tol)?);
    }
    let total = run_to[n];
    let share = |k: usize| if k == n { 1.0 } else { run_to[k] / total };
    // The law's frame at the spine's start, which the profile stands in.
    let start_frame = law_frame_at(model, spine, law, tolerance, tol)?;
    // Whether a join between runs `a` and `b` turns a sharp corner.
    let turns = |a: usize, b: usize| {
        sharp(
            legs[runs[a].legs.end - 1].heading_out,
            legs[runs[b].legs.start].heading_in,
            tol,
        )
    };
    let after_run = |r: usize| {
        if r + 1 < runs.len() {
            Some(r + 1)
        } else if closed {
            Some(0)
        } else {
            None
        }
    };
    let before_run = |r: usize| {
        if r > 0 {
            Some(r - 1)
        } else if closed {
            Some(runs.len() - 1)
        } else {
            None
        }
    };
    let mut wires: Vec<Shape> = Vec::with_capacity(runs.len());
    let mut windows: Vec<Option<Evolution>> = Vec::with_capacity(runs.len());
    let mut prisms: Vec<bool> = Vec::with_capacity(runs.len());
    // How far each run's solid may stand off the law's placement of its
    // sections: a prism places them exactly, a pipe shell holds its skin to
    // `tolerance`.
    let mut held: Vec<f64> = Vec::with_capacity(runs.len());
    for run in &runs {
        let first = &legs[run.legs.start];
        let wire = run_wire(model, &legs[run.legs.clone()], tol)?;
        let window = evolution
            .map(|e| e.window(share(run.legs.start), share(run.legs.end)))
            .filter(|w| !w.is_still());
        let prism = window.is_none()
            && run.legs.len() == 1
            && first.straight
            && keeps_frame(model, &wire, law, tolerance, tol)?;
        held.push(if prism { 0.0 } else { tolerance });
        wires.push(wire);
        windows.push(window);
        prisms.push(prism);
    }
    // Under a law that places each run's section itself, a prism meeting a
    // skinned run at a mitre is swept once the skinned run is, from the
    // rings that run meets the mitre with (`mitre_shared`).
    let deferred = |r: usize| {
        corners == PipeCorners::Mitre
            && !matches!(law, PipeLaw::RotationMinimizing)
            && prisms[r]
            && (before_run(r).is_some_and(|p| !prisms[p] && turns(p, r))
                || after_run(r).is_some_and(|q| !prisms[q] && turns(r, q)))
    };
    let mut sections: Vec<Option<Shape>> = vec![None; runs.len()];
    let mut solids: Vec<Option<Shape>> = vec![None; runs.len()];
    let mut ends: Vec<Option<Shape>> = vec![None; runs.len()];
    // The face each run starts with: its section, or, after a prism swept
    // once this run is, the skin's own fit of it where that stands off the
    // section.
    let mut starts: Vec<Option<Shape>> = vec![None; runs.len()];
    // Whether the prism before a run took rings from the face it starts
    // with: then that face, not the section, is what the run is run back
    // on from at the mitre.
    let mut shared_start: Vec<bool> = vec![false; runs.len()];
    let order: Vec<usize> = (0..runs.len())
        .filter(|&r| !deferred(r))
        .chain((0..runs.len()).filter(|&r| deferred(r)))
        .collect();
    for r in order {
        let run = &runs[r];
        let first = &legs[run.legs.start];
        let last = &legs[run.legs.end - 1];
        let wire = &wires[r];
        let section = match law {
            PipeLaw::RotationMinimizing if r > 0 => {
                let before = &legs[run.legs.start - 1];
                let turn = least_turn(before.end, before.heading_out, first.heading_in, tol)?;
                let Some(end) = &ends[r - 1] else {
                    ogeom_bail!(Construction, "a run is swept before the one it follows");
                };
                realized_profile_wound(model, end, &turn, Some(first.heading_in), tol)?
            }
            PipeLaw::RotationMinimizing => realized_profile_wound(
                model,
                profile,
                &Transform::IDENTITY,
                Some(first.heading_in),
                tol,
            )?,
            _ => {
                let here = law_frame_at(model, wire, law, tolerance, tol)?;
                let evolved = match evolution {
                    Some(evolution) => evolution.local(share(run.legs.start), tol)?,
                    None => Transform::IDENTITY,
                };
                let motion =
                    Transform::from_frame(&here) * evolved * Transform::to_frame(&start_frame);
                realized_profile_wound(model, profile, &motion, Some(first.heading_in), tol)?
            }
        };
        let section = if deferred(r) {
            let mut section = section;
            let mut kept: Vec<bool> = Vec::new();
            if let Some(p) = before_run(r).filter(|&p| !prisms[p] && turns(p, r))
                && let Some(end) = &ends[p]
            {
                let out = legs[runs[p].legs.end - 1].heading_out;
                let mirror = mirror_at(first.start, out, first.heading_in, tol)?;
                (section, kept) = mitre_shared(
                    model,
                    &section,
                    end,
                    &mirror,
                    first.heading_in,
                    held[p],
                    &kept,
                    tol,
                )?;
            }
            if let Some(q) = after_run(r).filter(|&q| !prisms[q] && turns(r, q))
                && let Some(start) = &starts[q]
            {
                let into = legs[runs[q].legs.start].heading_in;
                let back = Transform::translation(first.start - first.end);
                let mirror = back * mirror_at(last.end, last.heading_out, into, tol)?;
                let taken;
                (section, taken) = mitre_shared(
                    model,
                    &section,
                    start,
                    &mirror,
                    first.heading_in,
                    held[q],
                    &kept,
                    tol,
                )?;
                shared_start[q] = taken
                    .iter()
                    .enumerate()
                    .any(|(k, t)| *t && !kept.get(k).copied().unwrap_or(false));
            }
            section
        } else {
            section
        };
        let solid = if let Some(window) = &windows[r] {
            let lengths = &run_to[run.legs.start..=run.legs.end];
            evolved_stretches(
                model,
                &section,
                &legs[run.legs.clone()],
                lengths,
                law,
                window,
                tolerance,
                tol,
            )?
        } else if prisms[r] {
            make_prism(model, &section, first.end - first.start, tol)?.shape
        } else {
            make_pipe_shell_with(
                model,
                &section,
                wire,
                law,
                PipeCorners::Mitre,
                tolerance,
                tol,
            )?
            .shape
        };
        ends[r] = Some(end_section(model, &solid, last.end, last.heading_out, tol)?);
        let after_deferred = before_run(r).is_some_and(|p| deferred(p) && turns(p, r));
        starts[r] = Some(if !prisms[r] && after_deferred {
            let own = end_section(model, &solid, first.start, first.heading_in, tol)?;
            if same_section(model, &own, &section, tol)? {
                section.clone()
            } else {
                own
            }
        } else {
            section.clone()
        });
        sections[r] = Some(section);
        solids[r] = Some(solid);
    }
    let unwrap_all = |v: Vec<Option<Shape>>| -> OgeomResult<Vec<Shape>> {
        v.into_iter()
            .map(|s| s.ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a run was not swept")))
            .collect()
    };
    let sections = unwrap_all(sections)?;
    let mut solids = unwrap_all(solids)?;
    let ends = unwrap_all(ends)?;
    let starts = unwrap_all(starts)?;

    // The joins: after every run but the last, and the wrap of a closed
    // spine, smooth or sharp.
    let joins: Vec<(usize, usize)> = (0..runs.len())
        .filter_map(|r| {
            if r + 1 < runs.len() {
                Some((r, r + 1))
            } else if closed {
                Some((r, 0))
            } else {
                None
            }
        })
        .collect();
    let mut after: Vec<Vec<Shape>> = vec![Vec::new(); runs.len()];
    for &(a, b) in &joins {
        let out = &legs[runs[a].legs.end - 1];
        let into = &legs[runs[b].legs.start];
        let corner = out.end;
        let (d1, d2) = (out.heading_out, into.heading_in);
        if !sharp(d1, d2, tol) {
            // The wrap of a closed spine through a smooth point: the runs
            // butt on the section there, which must be one section.
            if !same_section(model, &ends[a], &sections[b], tol)? {
                ogeom_bail!(
                    Construction,
                    "the section comes home turned round the closed spine; the runs do \
                     not meet on one section"
                );
            }
            continue;
        }
        if d1.dot(d2) < DOUBLED_BACK {
            ogeom_bail!(
                Construction,
                "the spine all but doubles back at a corner; no corner of that kind turns it"
            );
        }
        let turn = d1.dot(d2).clamp(-1.0, 1.0).acos();
        match corners {
            PipeCorners::Mitre => {
                let halving = (d1 + d2) / (d1 + d2).magnitude();
                let reach_a = mitre_reach(model, &ends[a], corner, d1, halving, tol)?;
                let from_b = if shared_start[b] {
                    &starts[b]
                } else {
                    &sections[b]
                };
                let reach_b = mitre_reach(model, from_b, corner, d2, halving, tol)?;
                let on_a = make_prism(model, &ends[a], d1 * reach_a, tol)?.shape;
                let on_b = make_prism(model, from_b, -d2 * reach_b, tol)?.shape;
                let reach = reach_a.max(reach_b)
                    + bound_reach(model, &ends[a], corner, tol)?
                    + bound_reach(model, from_b, corner, tol)?;
                let a_kept = ogeom_bool::fuse(model, &solids[a], &on_a, tol)?.shape;
                let ahead = half_box(model, corner, halving, reach, tol)?;
                solids[a] = ogeom_bool::cut(model, &a_kept, &ahead, tol)?.shape;
                let b_kept = ogeom_bool::fuse(model, &solids[b], &on_b, tol)?.shape;
                let behind = half_box(model, corner, -halving, reach, tol)?;
                solids[b] = ogeom_bool::cut(model, &b_kept, &behind, tol)?.shape;
            }
            PipeCorners::Extended => {
                let past_a = reach_past(model, &sections[b], corner, d1, turn, tol)?;
                let past_b = reach_past(model, &ends[a], corner, -d2, turn, tol)?;
                let on_a = make_prism(model, &ends[a], d1 * past_a, tol)?.shape;
                let on_b = make_prism(model, &sections[b], -d2 * past_b, tol)?.shape;
                after[a].push(on_a);
                after[a].push(on_b);
            }
            PipeCorners::Round => {
                let about = Axis {
                    location: corner,
                    direction: Direction::new(d1.cross(d2), tol)?,
                };
                let rotation = Transform::rotation(about, turn);
                if !matches!(law, PipeLaw::RotationMinimizing) || b == 0 {
                    let turned = realized_profile_wound(model, &ends[a], &rotation, Some(d2), tol)?;
                    // The two sections are one within what either run is
                    // held to: a run's end stands off the law's placement
                    // by as much as its skin may.
                    let near = (tol.confusion() * 10.0).max(held[a]).max(held[b]);
                    if !lies_on(model, &turned, &sections[b], near, tol)?
                        || !lies_on(model, &sections[b], &turned, near, tol)?
                    {
                        ogeom_bail!(
                            Construction,
                            "a round corner turns the section about the corner, and the law \
                             places the next run's section otherwise"
                        );
                    }
                }
                let inward = d2 - d1 * d2.dot(d1);
                let outer = outside_half(model, &ends[a], corner, d1, -inward, tol)?;
                after[a].push(make_revolution(model, &outer, about, turn, tol)?.shape);
            }
        }
    }
    let mut pieces = Vec::with_capacity(runs.len() * 3);
    for (solid, joined) in solids.into_iter().zip(after) {
        pieces.push(solid);
        pieces.extend(joined);
    }
    let Some(shape) = fuse_in_order(model, pieces, tol)? else {
        ogeom_bail!(Construction, "the spine has no run");
    };
    let mut history = History::new();
    history.generate(spine, shape.clone());
    history.generate(profile, shape.clone());
    Ok(Built::new(shape, history))
}

/// The spine's edges in order, each with its ends and headings.
fn legs_of(model: &Model, spine: &Shape, tol: Tolerances) -> OgeomResult<Vec<Leg>> {
    let edges: Vec<Shape> = match model.kind_of(spine)? {
        ShapeType::Edge => vec![spine.clone()],
        _ => model.ordered_children_of(spine)?,
    };
    let mut legs = Vec::with_capacity(edges.len());
    for edge in edges {
        let (curve, range) = spine_curve_of(model, &edge)?;
        let reversed = edge.orientation() == ogeom_topo::Orientation::Reversed;
        let (t0, t1) = if reversed { (range.1, range.0) } else { range };
        let sense = if reversed { -1.0 } else { 1.0 };
        let heading = |t: f64| -> OgeomResult<Vector> {
            let d = curve.d1_at(t, tol)? * sense;
            let m = d.magnitude();
            if m <= tol.confusion() {
                ogeom_bail!(Construction, "the spine is degenerate at {t}");
            }
            Ok(d / m)
        };
        let straight = match &curve {
            ogeom_geom::Curve::Trimmed(t) => matches!(t.basis(), ogeom_geom::Curve::Line(_)),
            other => matches!(other, ogeom_geom::Curve::Line(_)),
        };
        legs.push(Leg {
            straight,
            start: curve.point_at(t0, tol)?,
            end: curve.point_at(t1, tol)?,
            heading_in: heading(t0)?,
            heading_out: heading(t1)?,
            edge,
        });
    }
    Ok(legs)
}

/// The run's edges as a wire of their own.
fn run_wire(model: &mut Model, legs: &[Leg], tol: Tolerances) -> OgeomResult<Shape> {
    let edges: Vec<Shape> = legs.iter().map(|l| l.edge.clone()).collect();
    Ok(ogeom_algo::make_wire(model, &edges, tol)?.shape)
}

/// The law's frame at the start of `spine`.
fn law_frame_at(
    model: &Model,
    spine: &Shape,
    law: &PipeLaw<'_>,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Frame> {
    let stations = shell_stations(model, spine, tol)?;
    let normals = law_normals(model, &stations, law, tolerance, tol)?;
    station_frame(&stations[0], normals[0], tol)
}

/// Whether the law keeps the section's frame all along a straight run:
/// the rotation-minimizing and binormal laws always do; an auxiliary
/// spine does where its normals agree at every station.
fn keeps_frame(
    model: &Model,
    run: &Shape,
    law: &PipeLaw<'_>,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    match law {
        PipeLaw::RotationMinimizing | PipeLaw::Binormal(_) => Ok(true),
        _ => {
            let stations = shell_stations(model, run, tol)?;
            let normals = law_normals(model, &stations, law, tolerance, tol)?;
            Ok(normals.iter().all(|n| {
                n.cross(normals[0]).magnitude() <= tol.angular() && n.dot(normals[0]) > 0.0
            }))
        }
    }
}

/// The least rotation about the corner at `at` carrying heading `from` to
/// heading `to`.
fn least_turn(at: Point, from: Vector, to: Vector, tol: Tolerances) -> OgeomResult<Transform> {
    let axis = from.cross(to);
    if axis.magnitude() <= tol.angular() {
        return Ok(Transform::IDENTITY);
    }
    let about = Axis {
        location: at,
        direction: Direction::new(axis, tol)?,
    };
    Ok(Transform::rotation(
        about,
        from.dot(to).clamp(-1.0, 1.0).acos(),
    ))
}

/// The face of a run's solid on the plane square to `heading` through the
/// run's end `at`: the section it ends with, as a face of its own.
fn end_section(
    model: &mut Model,
    solid: &Shape,
    at: Point,
    heading: Vector,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    for face in explore(model, solid, Filter::OfType(ShapeType::Face))? {
        let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data()) else {
            continue;
        };
        if !matches!(
            model.geometry().surface(data.surface),
            Some(SurfaceGeometry::Plane(_))
        ) {
            continue;
        }
        let (on, normal) = ogeom_algo::face_normal(model, &face, tol)?;
        if normal.cross(heading).magnitude() <= tol.angular()
            && (on - at).dot(heading).abs() <= tol.confusion() * 100.0
        {
            return realized_profile_wound(model, &face, &Transform::IDENTITY, Some(heading), tol);
        }
    }
    ogeom_bail!(
        Construction,
        "a run's solid has no end section square to it"
    )
}

/// Whether two sections are one: every sample round each ring of either
/// lies within ten confusions of the other's boundary.
fn same_section(model: &Model, a: &Shape, b: &Shape, tol: Tolerances) -> OgeomResult<bool> {
    let near = tol.confusion() * 10.0;
    Ok(lies_on(model, a, b, near, tol)? && lies_on(model, b, a, near, tol)?)
}

/// Whether every sample round each ring of `a` lies within `near` of an
/// edge of `b`.
fn lies_on(model: &Model, a: &Shape, b: &Shape, near: f64, tol: Tolerances) -> OgeomResult<bool> {
    for wire in model.ordered_children_of(a)? {
        if !ring_lies_on(model, &wire, b, near, tol)? {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Whether every sample round the ring `wire` lies within `near` of an
/// edge of `b`.
fn ring_lies_on(
    model: &Model,
    wire: &Shape,
    b: &Shape,
    near: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    const SAMPLES: usize = 64;
    let mut edges = Vec::new();
    for edge in explore(model, b, Filter::OfType(ShapeType::Edge))? {
        edges.push(spine_curve_of(model, &edge)?);
    }
    for p in sample_wire(model, wire, SAMPLES, tol)? {
        let mut best = f64::INFINITY;
        for (curve, range) in &edges {
            // Nearest on the whole curve where that falls in the edge's
            // range, else at one of its ends.
            let found = project_on_curve(curve, p, SAMPLES, tol)?;
            let inside = found.parameter >= range.0.min(range.1) - tol.parametric()
                && found.parameter <= range.0.max(range.1) + tol.parametric();
            let d = if inside {
                found.distance
            } else {
                p.distance(curve.point_at(range.0, tol)?)
                    .min(p.distance(curve.point_at(range.1, tol)?))
            };
            best = best.min(d);
        }
        if best > near {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The mirror in the plane through `corner` halving the turn from heading
/// `d1` to heading `d2`: it carries the plane square to `d1` there onto
/// the plane square to `d2`.
fn mirror_at(corner: Point, d1: Vector, d2: Vector, tol: Tolerances) -> OgeomResult<Transform> {
    let halving = Direction::new((d1 + d2) / (d1 + d2).magnitude(), tol)?;
    Ok(Transform::plane_mirror(corner, halving))
}

/// The section a prism sweeps where it meets a skinned run at a mitre,
/// with which of its rings were taken from that run.
///
/// A skinned run starts and ends on its skin's fit of the section, which
/// wobbles across the exact ring the law places. Run on straight past the
/// corner, the fit and the exact ring each cut the mitre plane, and the two
/// cuts cross there again and again a skin's width apart: slivers no
/// boolean can tell from a touch. So each ring of `placed` that a ring of
/// `skinned` (the skinned run's face at the corner), carried by `motion`
/// (the mitre's mirror, and back along the prism for its far end), matches
/// within `held` (what that run is held to) but not within ten confusions
/// is taken as that carried ring, unless `kept` marks it taken already;
/// every other ring stays as placed. The prism then cuts the mitre plane
/// where the skinned run does, and stands off the law's placement by at
/// most what that run already does.
#[allow(clippy::too_many_arguments)]
fn mitre_shared(
    model: &mut Model,
    placed: &Shape,
    skinned: &Shape,
    motion: &Transform,
    heading: Vector,
    held: f64,
    kept: &[bool],
    tol: Tolerances,
) -> OgeomResult<(Shape, Vec<bool>)> {
    let same = tol.confusion() * 10.0;
    let carried = realized_profile_wound(model, skinned, motion, Some(heading), tol)?;
    let theirs = model.ordered_children_of(&carried)?;
    let mut rings: Vec<Vec<Shape>> = Vec::new();
    let mut taken: Vec<bool> = Vec::new();
    for (k, ring) in model.ordered_children_of(placed)?.into_iter().enumerate() {
        let mut chosen = (ring.clone(), kept.get(k).copied().unwrap_or(false));
        if !chosen.1 {
            for other in &theirs {
                if ring_lies_on(model, other, &ring, held, tol)?
                    && ring_lies_on(model, &ring, other, held, tol)?
                    && !(ring_lies_on(model, other, &ring, same, tol)?
                        && ring_lies_on(model, &ring, other, same, tol)?)
                {
                    chosen = (other.clone(), true);
                    break;
                }
            }
        }
        rings.push(model.ordered_children_of(&chosen.0)?);
        taken.push(chosen.1);
    }
    if taken
        .iter()
        .zip(kept.iter().chain(core::iter::repeat(&false)))
        .all(|(t, k)| t == k)
    {
        return Ok((placed.clone(), taken));
    }
    let Some(NodeData::Face(data)) = model.node(placed).map(|n| n.data()) else {
        ogeom_bail!(Construction, "a run's section is not a face");
    };
    let Some(surface) = model.geometry().surface(data.surface).cloned() else {
        ogeom_bail!(Dangling, "surface is not in this model");
    };
    let face = ogeom_algo::make_face_with_pcurves(model, surface, &rings, tol)?.shape;
    Ok((face, taken))
}

/// The largest distance from `corner` to the section's bounding box: a
/// bound on how far any of its points stands from the corner.
fn bound_reach(model: &Model, section: &Shape, corner: Point, tol: Tolerances) -> OgeomResult<f64> {
    let bounds = shape_bounds(model, section, tol)?;
    let (Some(lo), Some(hi)) = (bounds.low(), bounds.high()) else {
        ogeom_bail!(Construction, "the section has no extent");
    };
    let mut far = 0.0_f64;
    for k in 0..8 {
        let p = Point::new(
            if k & 1 == 0 { lo.x } else { hi.x },
            if k & 2 == 0 { lo.y } else { hi.y },
            if k & 4 == 0 { lo.z } else { hi.z },
        );
        far = far.max(p.distance(corner));
    }
    Ok(far)
}

/// How far a section at `corner` is run on along `heading` for every point
/// of it to pass the plane through the corner square to `halving`: no
/// point stands farther from the corner than its bounding box's corners,
/// and a point that far needs that over `heading . halving` to cross.
fn mitre_reach(
    model: &Model,
    section: &Shape,
    corner: Point,
    heading: Vector,
    halving: Vector,
    tol: Tolerances,
) -> OgeomResult<f64> {
    let along = heading.dot(halving).abs();
    Ok(bound_reach(model, section, corner, tol)? / along.max(tol.angular()))
}

/// How far a run goes on past a corner for an extended corner: until its
/// section passes the far side of the other run's section there.
fn reach_past(
    model: &Model,
    other: &Shape,
    corner: Point,
    along: Vector,
    turn: f64,
    tol: Tolerances,
) -> OgeomResult<f64> {
    let wires = model.ordered_children_of(other)?;
    let far = sample_wire(model, &wires[0], AROUND, tol)?
        .iter()
        .map(|p| (*p - corner).dot(along))
        .fold(0.0_f64, f64::max);
    Ok(far / turn.sin().max(LEAST_TURN_SINE))
}

/// The box `reach` across about `corner` on the side `toward` of the plane
/// through it square to `toward`: what is cut from a run past its mitre.
fn half_box(
    model: &mut Model,
    corner: Point,
    toward: Vector,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let z = Direction::new(toward, tol)?;
    let seed = if toward.x.abs() < 0.9 {
        Vector::new(1.0, 0.0, 0.0)
    } else {
        Vector::new(0.0, 1.0, 0.0)
    };
    let x = toward.cross(seed);
    let x = Direction::new(x, tol)?;
    let y = z.vector().cross(x.vector());
    let origin = corner - x.vector() * reach - y * reach;
    let frame = Frame::new(origin, z, x, tol)?;
    Ok(make_box(model, frame, (2.0 * reach, 2.0 * reach, reach), tol)?.shape)
}

/// Whether a law's frames along a closed spine with no sharp corner are
/// the rotation-minimizing frames turned by one fixed angle: the angle
/// between the law's normal and the carried one the same at every
/// station, to an angle's tolerance. A frame carried round a ring that
/// comes home turned never is.
pub(super) fn law_is_rotation_minimizing(
    model: &Model,
    spine: &Shape,
    law: &PipeLaw<'_>,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<bool> {
    let stations = shell_stations(model, spine, tol)?;
    let normals = law_normals(model, &stations, law, tolerance, tol)?;
    let carried = rmf_normals(&stations);
    let angle = |i: usize| -> f64 {
        let (r, n) = (carried[i], normals[i]);
        r.cross(n).dot(stations[i].tangent).atan2(r.dot(n))
    };
    let first = angle(0);
    Ok((0..stations.len()).all(|i| {
        let d = (angle(i) - first).rem_euclid(2.0 * core::f64::consts::PI);
        d.min(2.0 * core::f64::consts::PI - d) <= tol.angular()
    }))
}

/// A run swept under the evolving laws `window`, read over the run: the
/// stretches at either end where the laws stand still swept exactly (the
/// prism of the section down a straight stretch the law keeps it level
/// along, else the pipe shell), the stretch between under the laws, each
/// from the section the one before ends with, the stretches fused in
/// order. `run_to` holds the length run along the spine to each leg's
/// start and to the run's end.
#[allow(clippy::too_many_arguments)]
fn evolved_stretches(
    model: &mut Model,
    section: &Shape,
    legs: &[Leg],
    run_to: &[f64],
    law: &PipeLaw<'_>,
    window: &Evolution,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Shape> {
    let (lead, trail) = window.still_ends();
    let (from, to) = (run_to[0], run_to[run_to.len() - 1]);
    let at = |share: f64| {
        if share == 1.0 {
            to
        } else {
            from + (to - from) * share
        }
    };
    let mut cuts = vec![0.0];
    if lead > 0.0 && lead < trail {
        cuts.push(lead);
    }
    if trail < 1.0 && lead < trail {
        cuts.push(trail);
    }
    cuts.push(1.0);
    let mut pieces = Vec::with_capacity(cuts.len() - 1);
    let mut start = section.clone();
    for w in cuts.windows(2) {
        let (a, b) = (w[0], w[1]);
        let (wire, last) = stretch_wire(model, legs, run_to, (at(a), at(b)), tol)?;
        let moving = a < trail && b > lead;
        let solid = if moving {
            let laws = if a == 0.0 && b == 1.0 {
                window.clone()
            } else {
                window.window(a, b)
            };
            evolved_run(model, &start, &wire, law, &laws, tolerance, tol)?
        } else if legs.len() == 1
            && legs[0].straight
            && keeps_frame(model, &wire, law, tolerance, tol)?
        {
            let (p, q) = (
                point_along(model, legs, run_to, at(a), tol)?,
                point_along(model, legs, run_to, at(b), tol)?,
            );
            make_prism(model, &start, q - p, tol)?.shape
        } else {
            make_pipe_shell_with(
                model,
                &start,
                &wire,
                law,
                PipeCorners::Mitre,
                tolerance,
                tol,
            )?
            .shape
        };
        let (end, heading) = last;
        start = end_section(model, &solid, end, heading, tol)?;
        pieces.push(solid);
    }
    match fuse_in_order(model, pieces, tol)? {
        Some(shape) => Ok(shape),
        None => ogeom_bail!(Construction, "a run has no stretch"),
    }
}

/// The point the spine reaches at length `along`, among the run's legs.
fn point_along(
    model: &Model,
    legs: &[Leg],
    run_to: &[f64],
    along: f64,
    tol: Tolerances,
) -> OgeomResult<Point> {
    let k = run_to[1..legs.len()].partition_point(|r| *r <= along);
    let (curve, t) = leg_parameter(model, &legs[k], along - run_to[k], tol)?;
    curve.point_at(t, tol)
}

/// A leg's curve, and the parameter a length `run` along it from its
/// start in the direction of travel.
fn leg_parameter(
    model: &Model,
    leg: &Leg,
    run: f64,
    tol: Tolerances,
) -> OgeomResult<(ogeom_geom::Curve, f64)> {
    let (curve, range) = spine_curve_of(model, &leg.edge)?;
    let length = curve_length(&curve, range, tol)?;
    let run = run.clamp(0.0, length);
    let t = if leg.edge.orientation() == ogeom_topo::Orientation::Reversed {
        parameter_at_length(&curve, range, length - run, tol)?
    } else {
        parameter_at_length(&curve, range, run, tol)?
    };
    Ok((curve, t))
}

/// The wire of the run's spine from length `span.0` to length `span.1`
/// along it, each leg it crosses cut where the stretch ends, with the
/// stretch's end point and its heading there.
fn stretch_wire(
    model: &mut Model,
    legs: &[Leg],
    run_to: &[f64],
    span: (f64, f64),
    tol: Tolerances,
) -> OgeomResult<(Shape, (Point, Vector))> {
    let mut edges = Vec::new();
    let mut joint: Option<Shape> = None;
    let mut last = (legs[0].start, legs[0].heading_in);
    for (k, leg) in legs.iter().enumerate() {
        let (lo, hi) = (span.0.max(run_to[k]), span.1.min(run_to[k + 1]));
        if hi <= lo {
            continue;
        }
        let whole = lo == run_to[k] && hi == run_to[k + 1];
        let (curve, t0) = leg_parameter(model, leg, lo - run_to[k], tol)?;
        let (_, t1) = leg_parameter(model, leg, hi - run_to[k], tol)?;
        let (a, b) = (curve.point_at(t0, tol)?, curve.point_at(t1, tol)?);
        let edge = if whole {
            leg.edge.clone()
        } else {
            let va = match &joint {
                Some(v) => v.clone(),
                None => ogeom_algo::make_vertex(model, a).shape,
            };
            let vb = ogeom_algo::make_vertex(model, b).shape;
            if leg.edge.orientation() == ogeom_topo::Orientation::Reversed {
                ogeom_algo::make_edge_between(model, curve.clone(), (t1, t0), &vb, &va, tol)?
                    .shape
                    .reversed()
            } else {
                ogeom_algo::make_edge_between(model, curve.clone(), (t0, t1), &va, &vb, tol)?.shape
            }
        };
        joint = ogeom_algo::edge_vertices(model, &edge)?.map(|(_, end)| end);
        let d = curve.d1_at(t1, tol)?;
        let d = if leg.edge.orientation() == ogeom_topo::Orientation::Reversed {
            -d
        } else {
            d
        };
        last = (b, d / d.magnitude());
        edges.push(edge);
    }
    Ok((ogeom_algo::make_wire(model, &edges, tol)?.shape, last))
}
