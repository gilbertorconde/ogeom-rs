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
//!   two meet on that plane whatever their sections are;
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
//! them, so the solid is exact wherever the runs are swept exactly.

use ogeom_algo::{
    Built, History, make_box, make_prism, make_revolution, project_on_curve, shape_bounds,
};
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::{Curve3d as _, SurfaceGeometry};
use ogeom_math::{Axis, Direction, Frame, Point, Transform, Vector};
use ogeom_topo::{Filter, Model, NodeData, Shape, ShapeType, explore};

use super::{
    DOUBLED_BACK, LEAST_TURN_SINE, PipeCorners, PipeLaw, fuse_in_order, law_normals,
    make_pipe_shell_with, outside_half, profile_square_at_start, realized_profile_wound,
    rmf_normals, sample_wire, shell_stations, spine_curve_of, station_frame,
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
/// asks; see the module documentation.
pub(super) fn law_cornered(
    model: &mut Model,
    profile: &Shape,
    spine: &Shape,
    law: &PipeLaw<'_>,
    corners: PipeCorners,
    tolerance: f64,
    tol: Tolerances,
) -> OgeomResult<Built> {
    if model.kind_of(profile)? != ShapeType::Face {
        ogeom_bail!(
            Construction,
            "a pipe turned round a sharp corner run by run sweeps a face"
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

    // The law's frame at the spine's start, which the profile stands in.
    let start_frame = law_frame_at(model, spine, law, tolerance, tol)?;
    let mut sections: Vec<Shape> = Vec::with_capacity(runs.len());
    let mut solids: Vec<Shape> = Vec::with_capacity(runs.len());
    let mut ends: Vec<Shape> = Vec::with_capacity(runs.len());
    for (r, run) in runs.iter().enumerate() {
        let first = &legs[run.legs.start];
        let last = &legs[run.legs.end - 1];
        let wire = run_wire(model, &legs[run.legs.clone()], tol)?;
        let section = match law {
            PipeLaw::RotationMinimizing if r > 0 => {
                let before = &legs[run.legs.start - 1];
                let turn = least_turn(before.end, before.heading_out, first.heading_in, tol)?;
                realized_profile_wound(model, &ends[r - 1], &turn, Some(first.heading_in), tol)?
            }
            PipeLaw::RotationMinimizing => realized_profile_wound(
                model,
                profile,
                &Transform::IDENTITY,
                Some(first.heading_in),
                tol,
            )?,
            _ => {
                let here = law_frame_at(model, &wire, law, tolerance, tol)?;
                let motion = Transform::from_frame(&here) * Transform::to_frame(&start_frame);
                realized_profile_wound(model, profile, &motion, Some(first.heading_in), tol)?
            }
        };
        let solid = if run.legs.len() == 1
            && first.straight
            && keeps_frame(model, &wire, law, tolerance, tol)?
        {
            make_prism(model, &section, first.end - first.start, tol)?.shape
        } else {
            make_pipe_shell_with(
                model,
                &section,
                &wire,
                law,
                PipeCorners::Mitre,
                tolerance,
                tol,
            )?
            .shape
        };
        ends.push(end_section(model, &solid, last.end, last.heading_out, tol)?);
        sections.push(section);
        solids.push(solid);
    }

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
                let reach_b = mitre_reach(model, &sections[b], corner, d2, halving, tol)?;
                let on_a = make_prism(model, &ends[a], d1 * reach_a, tol)?.shape;
                let on_b = make_prism(model, &sections[b], -d2 * reach_b, tol)?.shape;
                let reach = reach_a.max(reach_b)
                    + bound_reach(model, &ends[a], corner, tol)?
                    + bound_reach(model, &sections[b], corner, tol)?;
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
                    if !same_section(model, &turned, &sections[b], tol)? {
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
    const SAMPLES: usize = 64;
    let mut edges = Vec::new();
    for edge in explore(model, b, Filter::OfType(ShapeType::Edge))? {
        edges.push(spine_curve_of(model, &edge)?);
    }
    for wire in model.ordered_children_of(a)? {
        for p in sample_wire(model, &wire, SAMPLES, tol)? {
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
    }
    Ok(true)
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
