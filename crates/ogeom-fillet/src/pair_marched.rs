//! Rounds between two faces that share no edge where the ball's seat has no
//! closed form: surfaces at a slant to each other, or a B-spline face.
//!
//! With no edge of the solid to run along, the ball is guided by where the
//! two surfaces cross (the corner it rounds) or, between faces whose
//! surfaces do not cross, by where their offsets by the radius cross (the
//! path of the ball's centre). It is marched along that guide on the two
//! surfaces, carried on past the faces, and the round spans the stations
//! where it touches both faces: the band fitted through the ball's arcs,
//! as the marched edge blend fits it.
//!
//! Where the ball touches each face is a line on it: the touch points, at
//! the guide's parameters, fitted together with their places in the face's
//! chart. A face the ball leaves is cut along that line, which settles
//! exactly where the line enters and leaves the face, and the round spans
//! the overlap of the two faces' stretches.
//!
//! Between faces of a solid the corner between the band and the crease is
//! built as the marched edge blend builds it and cut off or fused on: a
//! ring where the ball touches both faces all the way round a closed seat,
//! and otherwise capped at each end of the run in the ball's section there,
//! as the open marched edge blend caps its run. Between separate faces each
//! face is cut back to its line of contact and the round sewn between them.

use ogeom_algo::{
    Built, History, attach_pcurve, make_edge, make_edge_between, make_vertex, project_on_curve,
    project_on_surface, shape_bounds,
};
use ogeom_core::{OgeomError, OgeomResult, Tolerance, Tolerances, ogeom_bail};
use ogeom_geom::{
    Curve, Curve3d as _, CylinderSurface, OffsetSurface, PlanarCurve, PlaneSurface, Surface as _,
    SurfaceGeometry,
};
use ogeom_intersect::{IntersectOptions, Marching, SurfaceIntersection};
use ogeom_math::{Point, Point2, Vector};
use ogeom_topo::{Location, Model, Orientation, Shape, ShapeType, SurfaceId, explore_unique};

use crate::march::{BlendStop, MarchedBlend, Sides, march_blend_seeded, seat_section};
use crate::marched::{band_fit_target, fit_open_band};
use crate::sheet_curved::{Marched, cut_to_run, marched_round, unplaced_face};
use crate::support::edge_curve;

/// The two faces' surfaces as the ball rolls on them, carried on past the
/// faces: a plane's or a cylinder's window widened, a B-spline patch
/// continued as itself on every open side. Each is written back as the
/// face's own surface, the same parameters over a wider window, so the
/// lines of contact cut the faces on the surfaces they were marched on.
#[derive(Clone)]
pub(crate) struct Hosts {
    ids: [SurfaceId; 2],
    surfaces: [SurfaceGeometry; 2],
    /// One where the face runs with its surface's normal, minus one where
    /// it runs against it.
    signs: [f64; 2],
}

/// The two faces' surfaces, carried on past the faces once for every seat
/// read between them.
pub(crate) fn hosts_of(
    model: &mut Model,
    faces: [&Shape; 2],
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Hosts> {
    let mut reach = radius * 6.0;
    for face in faces {
        reach += shape_bounds(model, face, tol)?.diagonal();
    }
    let mut ids = Vec::with_capacity(2);
    let mut surfaces = Vec::with_capacity(2);
    let mut signs = [1.0; 2];
    for (f, face) in faces.iter().enumerate() {
        let (id, surface, _) = unplaced_face(model, face)?;
        let surface = match surface {
            SurfaceGeometry::Plane(plane) => {
                let ((u0, u1), (v0, v1)) = plane.domain();
                PlaneSurface::over(
                    plane.plane(),
                    (u0 - reach, u1 + reach),
                    (v0 - reach, v1 + reach),
                )?
                .into()
            }
            SurfaceGeometry::Cylinder(cylinder) => {
                let (_, (h0, h1)) = cylinder.domain();
                CylinderSurface::new(cylinder.cylinder(), (h0 - reach, h1 + reach))?.into()
            }
            SurfaceGeometry::BSpline(patch) => {
                let mut longer = patch;
                for (along_u, closed) in [
                    (true, longer.is_closed_u(tol) || longer.is_periodic_u()),
                    (false, longer.is_closed_v(tol) || longer.is_periodic_v()),
                ] {
                    if closed {
                        continue;
                    }
                    for at_end in [false, true] {
                        if let Ok(grown) = longer.extended(along_u, at_end, radius * 6.0, 2, tol) {
                            longer = grown;
                        }
                    }
                }
                SurfaceGeometry::BSpline(longer)
            }
            SurfaceGeometry::Trimmed(_) | SurfaceGeometry::Offset(_) => ogeom_bail!(
                Construction,
                "a round between faces that share no edge rolls on surfaces with their own \
                 chart; a trimmed or offset face is refused"
            ),
            other => other,
        };
        if let Some(held) = model.geometry_mut().surface_mut(id) {
            *held = surface.clone();
        }
        if face.orientation() == Orientation::Reversed {
            signs[f] = -1.0;
        }
        ids.push(id);
        surfaces.push(surface);
    }
    let [first, second] = [surfaces[0].clone(), surfaces[1].clone()];
    Ok(Hosts {
        ids: [ids[0], ids[1]],
        surfaces: [first, second],
        signs,
    })
}

/// Where a face's stretch of the ball's line of contact lies on it.
pub(crate) struct FaceStretch {
    /// The piece of the face on the far side of the line from the round.
    piece: Shape,
    /// The edge the line cut it along.
    edge: Shape,
    /// The line's curve, parameterized as the guide.
    line: Curve,
    /// The guide's parameters where the line enters and leaves the face.
    span: (f64, f64),
}

/// How much of the marched seat the round takes.
pub(crate) enum Run {
    /// The whole loop: the ball touches both faces all the way round.
    Closed,
    /// The overlap of the faces' two stretches.
    Open(Box<[FaceStretch; 2]>),
    /// Between faces of a solid, the stretch of the seat where the ball
    /// touches both, ended where it leaves one: the corner is capped in
    /// the ball's section at each end.
    Capped,
}

/// A marched seat between two faces: the stations the round spans, the
/// guide they were marched along, and the surfaces they stand on.
pub(crate) struct PairSeat {
    hosts: Hosts,
    guide: Curve,
    /// The stations: the loop for a closed run, the run with its exact end
    /// sections otherwise, in order along the guide.
    blend: MarchedBlend,
    run: Run,
    /// For a seat between faces of a solid, a point just off the middle of
    /// the ball's arc toward the crease, where the solid says whether the
    /// corner is material.
    probe: Option<Point>,
}

impl PairSeat {
    pub(crate) const fn probe(&self) -> Option<Point> {
        self.probe
    }
}

/// The marched seat of a ball of `radius` touching both faces, on the side
/// each face's normal points to, or behind both with `behind`. `hosts` are
/// the faces' surfaces as [`hosts_of`] carries them on.
///
/// Between faces of a solid (`solid`) the corner is where the surfaces
/// cross: only that crossing guides the ball, and the round runs all the
/// way round a closed seat or ends where the ball leaves a face. Between
/// separate faces, where the surfaces do not cross, their offsets'
/// crossing (the ball's own centre line) guides it.
pub(crate) fn pair_seat(
    model: &mut Model,
    faces: [&Shape; 2],
    hosts: &Hosts,
    radius: f64,
    behind: bool,
    solid: bool,
    tol: Tolerances,
) -> OgeomResult<PairSeat> {
    let hosts = hosts.clone();
    let seat = if behind { -1.0 } else { 1.0 };
    #[allow(clippy::cast_possible_truncation, reason = "a sign")]
    let sides = Sides {
        first: (seat * hosts.signs[0]) as i8,
        second: (seat * hosts.signs[1]) as i8,
    };
    let [first, second] = &hosts.surfaces;

    // The guides: where the surfaces cross, or where their offsets do.
    let options = IntersectOptions::default();
    let crossing = |a: &SurfaceGeometry, b: &SurfaceGeometry| -> OgeomResult<Vec<Curve>> {
        Ok(
            match ogeom_intersect::intersect_surfaces(a, b, options, tol)? {
                SurfaceIntersection::Along(curves) => curves
                    .into_iter()
                    .filter(|c| !c.tangential)
                    .map(|c| c.curve)
                    .collect(),
                _ => Vec::new(),
            },
        )
    };
    let mut guides = crossing(first, second)?;
    if guides.is_empty() {
        if solid {
            ogeom_bail!(
                Construction,
                "the two faces' surfaces do not meet, so there is no corner between them for \
                 the round to take off or fill"
            );
        }
        let offset = |surface: &SurfaceGeometry, side: i8| -> OgeomResult<SurfaceGeometry> {
            Ok(SurfaceGeometry::Offset(Box::new(OffsetSurface::new(
                surface.clone(),
                f64::from(side) * radius,
            )?)))
        };
        guides = crossing(&offset(first, sides.first)?, &offset(second, sides.second)?)?;
        if guides.is_empty() {
            ogeom_bail!(
                Construction,
                "no ball of radius {radius} touches both faces' surfaces: neither the surfaces \
                 nor their offsets by the radius meet"
            );
        }
    }

    let mut found = Vec::new();
    let mut why: Option<OgeomError> = None;
    for guide in guides {
        match seat_along(model, faces, &hosts, &guide, sides, radius, solid, tol) {
            Ok(Some((blend, run, probe))) => found.push((guide, blend, run, probe)),
            Ok(None) => {}
            Err(e @ (OgeomError::Construction(_) | OgeomError::NotDone(_))) => why = Some(e),
            Err(e) => return Err(e),
        }
    }
    match found.len() {
        1 => {
            let (guide, blend, run, probe) = found.remove(0);
            Ok(PairSeat {
                hosts,
                guide,
                blend,
                run,
                probe,
            })
        }
        0 => Err(why.unwrap_or_else(|| {
            ogeom_core::ogeom_err!(
                Construction,
                "a ball of radius {radius} rolling along where the surfaces meet does not touch \
                 both faces; it misses at least one of them"
            )
        })),
        _ => ogeom_bail!(
            Construction,
            "a ball of radius {radius} touches both faces along more than one corner; which \
             corner to round is ambiguous"
        ),
    }
}

/// The stations' fields, reordered by `order`.
fn reordered(blend: &MarchedBlend, order: &[usize]) -> MarchedBlend {
    MarchedBlend {
        spine: order.iter().map(|&i| blend.spine[i]).collect(),
        on_first: order.iter().map(|&i| blend.on_first[i]).collect(),
        on_second: order.iter().map(|&i| blend.on_second[i]).collect(),
        touch_first: order.iter().map(|&i| blend.touch_first[i]).collect(),
        touch_second: order.iter().map(|&i| blend.touch_second[i]).collect(),
        along: order.iter().map(|&i| blend.along[i]).collect(),
        sides: blend.sides,
        stopped: blend.stopped,
    }
}

/// The deflection faces are read at for whether the ball touches them.
fn reading(radius: f64, tol: Tolerances) -> ogeom_mesh::Deflection {
    ogeom_mesh::Deflection {
        chord: (radius * 1e-3).max(tol.confusion() * 1e2),
        ..ogeom_mesh::Deflection::default()
    }
}

/// Of `count` stations, at most `most` spread evenly.
fn spread(count: usize, most: usize) -> Vec<usize> {
    if count <= most {
        return (0..count).collect();
    }
    (0..most).map(|k| k * (count - 1) / (most - 1)).collect()
}

/// The march along one guide, and the run it gives: `None` where the ball
/// does not touch both faces along it. With `whole`, for faces of a solid,
/// the run is the closed seat or the stretch of it the faces give ends to,
/// and the probe beside it is read.
#[allow(clippy::too_many_arguments, reason = "one seat, all its data")]
#[allow(clippy::type_complexity, reason = "the seat's three parts")]
fn seat_along(
    model: &mut Model,
    faces: [&Shape; 2],
    hosts: &Hosts,
    guide: &Curve,
    sides: Sides,
    radius: f64,
    whole: bool,
    tol: Tolerances,
) -> OgeomResult<Option<(MarchedBlend, Run, Option<Point>)>> {
    let [first, second] = &hosts.surfaces;
    // Seeded where the guide passes nearest each face, then at its middle.
    let mut seeds = Vec::with_capacity(3);
    for face in faces {
        let (p, _) = ogeom_algo::face_normal(model, face, tol)?;
        seeds.push(project_on_curve(guide, p, 256, tol)?.parameter);
    }
    let (lo, hi) = guide.domain();
    seeds.push(f64::midpoint(lo, hi));
    let options = Marching {
        // The blend's own scale, as the marched edge blend steps it.
        chord: (radius * 3e-6).max(tol.confusion() * 0.1),
        ..Marching::default()
    };
    let mut marched = None;
    for seed in seeds {
        match march_blend_seeded(first, second, radius, guide, sides, seed, options, tol) {
            Ok(blend) if blend.len() >= 8 => {
                marched = Some(blend);
                break;
            }
            Ok(_) | Err(OgeomError::NotDone(_) | OgeomError::Construction(_)) => {}
            Err(e) => return Err(e),
        }
    }
    let Some(blend) = marched else {
        ogeom_bail!(
            Construction,
            "no ball of radius {radius} seats on that side of both faces' surfaces along where \
             they meet"
        );
    };

    // Whether the ball touches each face at a spread of stations, and a
    // station where it misses one. Between faces of a solid every station
    // is read: where each face lets the ball go settles where the round
    // ends.
    let deflection = reading(radius, tol);
    let mut touching = [0usize; 2];
    let mut missing = None;
    let mut missed = [vec![false; blend.len()], vec![false; blend.len()]];
    let sampled = if whole {
        (0..blend.len()).collect()
    } else {
        spread(blend.len(), 128)
    };
    for &i in &sampled {
        for (f, face) in faces.iter().enumerate() {
            let p = if f == 0 {
                blend.touch_first[i]
            } else {
                blend.touch_second[i]
            };
            if ogeom_algo::classify_on_face(model, face, p, deflection, tol)?
                == ogeom_algo::Containment::Out
            {
                missing.get_or_insert(i);
                missed[f][i] = true;
            } else {
                touching[f] += 1;
            }
        }
    }
    if touching.contains(&0) {
        return Ok(None);
    }

    let closed = blend.stopped == BlendStop::Closed;
    let (blend, missed) = match (closed, missing) {
        (true, None) => {
            let probe = whole
                .then(|| probe_at(&blend, guide, blend.len() / 2, radius, tol))
                .transpose()?;
            return Ok(Some((blend, Run::Closed, probe)));
        }
        // A loop the faces hold part of: opened where the ball misses a
        // face, and read from there round as one run. Between faces of a
        // solid it is opened where the ball misses every face it misses
        // anywhere, so each such face's line of contact enters and leaves
        // it once along the run.
        (true, Some(at)) => {
            let n = blend.len();
            let at = if whole {
                // The middle of the longest stretch where the ball is off
                // every face it leaves, so the run starts and ends inside it.
                let left: Vec<usize> = (0..2).filter(|&f| missed[f].contains(&true)).collect();
                let gone = |i: usize| left.iter().all(|&f| missed[f][i % n]);
                let Some(start) = (0..n).find(|&i| gone(i) && !gone(i + n - 1)) else {
                    ogeom_bail!(
                        Construction,
                        "the ball leaves the two faces at different places round its seat; a \
                         round between faces of a solid is ended where the ball leaves both, or \
                         one face it leaves while the other holds it throughout"
                    );
                };
                let mut best = (0, start);
                for i in (0..n).filter(|&i| gone(i) && !gone(i + n - 1)) {
                    let length = (0..n).take_while(|&k| gone(i + k)).count();
                    if length > best.0 {
                        best = (length, i);
                    }
                }
                (best.1 + best.0 / 2) % n
            } else {
                at
            };
            let order: Vec<usize> = (0..n).map(|k| (at + k) % n).collect();
            let missed = [
                order.iter().map(|&i| missed[0][i]).collect(),
                order.iter().map(|&i| missed[1][i]).collect(),
            ];
            (reordered(&blend, &order), missed)
        }
        (false, _) => (blend, missed),
    };

    // In order along the guide: a looping guide's parameters unwrapped
    // into one sweep, then sorted, every station that fails to advance
    // dropped.
    let mut along = blend.along.clone();
    let loops = guide.is_periodic()
        || guide
            .point_at(lo, tol)
            .and_then(|p| guide.point_at(hi, tol).map(|q| p.distance(q)))
            .is_ok_and(|d| d <= tol.confusion() * 10.0);
    if loops {
        let period = hi - lo;
        for i in 1..along.len() {
            let mut t = along[i];
            while t - along[i - 1] > period / 2.0 {
                t -= period;
            }
            while along[i - 1] - t > period / 2.0 {
                t += period;
            }
            along[i] = t;
        }
    }
    let mut order: Vec<usize> = (0..along.len()).collect();
    order.sort_by(|&i, &j| along[i].total_cmp(&along[j]));
    order.dedup_by(|j, i| along[*j] <= along[*i] + tol.parametric());
    let mut blend = reordered(&blend, &order);
    blend.along = order.iter().map(|&i| along[i]).collect();
    let missed: [Vec<bool>; 2] = [
        order.iter().map(|&i| missed[0][i]).collect(),
        order.iter().map(|&i| missed[1][i]).collect(),
    ];

    if whole {
        // A face the ball never leaves holds the whole run; one it leaves
        // is cut along its line of contact, which settles exactly where
        // the line enters and leaves it. The round spans the overlap, and
        // each of its ends must be where the ball leaves a face: where the
        // march stopped on both faces is no end the faces give.
        let (first_at, last_at) = (blend.along[0], blend.along[blend.len() - 1]);
        let (mut w0, mut w1) = (first_at, last_at);
        for (f, face) in faces.iter().enumerate() {
            if !missed[f].contains(&true) {
                continue;
            }
            let Some(cut) = contact_cut(model, face, f, &blend, hosts, tol)? else {
                return Ok(None);
            };
            w0 = w0.max(cut.span.0);
            w1 = w1.min(cut.span.1);
        }
        let margin = (last_at - first_at) * 1e-6;
        if w0 <= first_at + margin || w1 >= last_at - margin {
            ogeom_bail!(
                Construction,
                "the march along where the surfaces meet stops while the ball still touches \
                 both faces; the round has no end the faces give it"
            );
        }
        if w1 - w0 <= tol.parametric() {
            ogeom_bail!(
                Construction,
                "the two faces reach the round over stretches of the corner that do not overlap"
            );
        }
        let run = trimmed(&blend, hosts, guide, sides, radius, (w0, w1), tol)?;
        let probe = probe_at(&run, guide, run.len() / 2, radius, tol)?;
        return Ok(Some((run, Run::Capped, Some(probe))));
    }

    let mut stretches = Vec::with_capacity(2);
    for (f, face) in faces.iter().enumerate() {
        stretches.push(stretch_on(model, face, f, &blend, hosts, tol)?);
    }
    let [Some(a), Some(b)] = [stretches.remove(0), stretches.remove(0)] else {
        return Ok(None);
    };
    let (w0, w1) = (a.span.0.max(b.span.0), a.span.1.min(b.span.1));
    if w1 - w0 <= tol.parametric() {
        ogeom_bail!(
            Construction,
            "the two faces reach the round over stretches of the corner that do not overlap"
        );
    }
    let run = trimmed(&blend, hosts, guide, sides, radius, (w0, w1), tol)?;
    Ok(Some((run, Run::Open(Box::new([a, b])), None)))
}

/// Just off the middle of the ball's arc at station `i`, toward the crease.
fn probe_at(
    blend: &MarchedBlend,
    guide: &Curve,
    i: usize,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Point> {
    let centre = blend.spine[i];
    let bisector = (blend.touch_first[i] - centre) + (blend.touch_second[i] - centre);
    let on_arc = centre + bisector.normalized(tol)? * radius;
    let crease = guide.point_at(on_guide(guide, blend.along[i]), tol)?;
    let toward = crease - on_arc;
    let gap = toward.magnitude();
    if gap <= tol.confusion() {
        ogeom_bail!(
            Construction,
            "the ball's arc runs through the crease; there is no corner between the round and \
             the crease"
        );
    }
    Ok(on_arc + toward / gap * (gap.min(radius) * 0.1))
}

/// A face cut along the ball's line of contact.
struct ContactCut {
    /// The face split along the line.
    split: Shape,
    /// The edge the line cut the face along.
    edge: Shape,
    /// The line's curve, parameterized as the guide.
    line: Curve,
    /// The guide's parameters where the line enters and leaves the face.
    span: (f64, f64),
}

/// The face's stretch of the ball's line of contact: the face cut along
/// the line, and the piece away from the round kept. `None` where the
/// line misses the face.
fn stretch_on(
    model: &mut Model,
    face: &Shape,
    f: usize,
    blend: &MarchedBlend,
    hosts: &Hosts,
    tol: Tolerances,
) -> OgeomResult<Option<FaceStretch>> {
    let surface = &hosts.surfaces[f];
    let Some(ContactCut {
        split,
        edge: cut,
        line,
        span,
    }) = contact_cut(model, face, f, blend, hosts, tol)?
    else {
        return Ok(None);
    };
    let (touch, other) = if f == 0 {
        (&blend.touch_first, &blend.touch_second)
    } else {
        (&blend.touch_second, &blend.touch_first)
    };
    let (cut_curve, cut_range) = edge_curve(model, &cut, tol)?;

    // The piece away from the round: probed just off the line, half way
    // along the stretch, on the far side from the round.
    let middle = f64::midpoint(span.0, span.1);
    let i = blend
        .along
        .partition_point(|&t| t < middle)
        .min(blend.len() - 1);
    let (centre, at, to) = (blend.spine[i], touch[i], other[i]);
    let (a, b) = (at - centre, to - centre);
    let into_arc = b - a * (a.dot(b) / a.dot(a));
    let away: Vector = -into_arc.normalized(tol)?;
    let reach = cut_curve
        .point_at(cut_range.0, tol)?
        .distance(cut_curve.point_at(cut_range.1, tol)?)
        .max(tol.confusion() * 1e3)
        .min(a.magnitude());
    let pieces = explore_unique(model, &split, ShapeType::Face)?;
    for scale in [1e-3, 1e-2, 5e-2] {
        let eps = reach * scale;
        let deflection = ogeom_mesh::Deflection {
            chord: eps * 0.1,
            ..ogeom_mesh::Deflection::default()
        };
        let probe = project_on_surface(surface, at + away * eps, 32, tol)?.point;
        let mut holding = Vec::new();
        for piece in &pieces {
            if ogeom_algo::classify_on_face(model, piece, probe, deflection, tol)?
                == ogeom_algo::Containment::In
            {
                holding.push(piece.clone());
            }
        }
        let [piece] = holding.as_slice() else {
            continue;
        };
        if !explore_unique(model, piece, ShapeType::Edge)?
            .iter()
            .any(|e| e.is_same(&cut))
        {
            return Ok(None);
        }
        return Ok(Some(FaceStretch {
            piece: piece.clone(),
            edge: cut,
            line,
            span,
        }));
    }
    Ok(None)
}

/// The face cut along the ball's line of contact: the touch points fitted
/// at the guide's parameters with their places in the face's chart.
/// `None` where the line misses the face.
fn contact_cut(
    model: &mut Model,
    face: &Shape,
    f: usize,
    blend: &MarchedBlend,
    hosts: &Hosts,
    tol: Tolerances,
) -> OgeomResult<Option<ContactCut>> {
    let surface = &hosts.surfaces[f];
    let (touch, on) = if f == 0 {
        (&blend.touch_first, &blend.on_first)
    } else {
        (&blend.touch_second, &blend.on_second)
    };
    let target = tol.confusion() * 10.0;
    let line = ogeom_geom::fit::fit_points_at(&blend.along, touch, 3, target, tol)?;
    if !line.met {
        ogeom_bail!(
            NotDone,
            "the ball's line of contact on a face reached {} against a target of {target}",
            line.error
        );
    }
    // The chart image at the same parameters, unwrapped across a periodic
    // chart's seam.
    let ((u0, u1), (v0, v1)) = surface.domain();
    let mut chart: Vec<Point2> = on.iter().map(|&(u, v)| Point2::new(u, v)).collect();
    for (periodic, span, pick) in [
        (surface.is_periodic_u(), u1 - u0, 0),
        (surface.is_periodic_v(), v1 - v0, 1),
    ] {
        if !periodic {
            continue;
        }
        for i in 1..chart.len() {
            let get = |p: Point2| if pick == 0 { p.x } else { p.y };
            let mut x = get(chart[i]);
            let last = get(chart[i - 1]);
            while x - last > span / 2.0 {
                x -= span;
            }
            while last - x > span / 2.0 {
                x += span;
            }
            if pick == 0 {
                chart[i].x = x;
            } else {
                chart[i].y = x;
            }
        }
    }
    let image = ogeom_geom::fit::fit_points_2d_at(&blend.along, &chart, 3, target * 1e-3, tol)?;
    if !image.met {
        ogeom_bail!(
            NotDone,
            "the chart image of the ball's line of contact on a face reached {} against a \
             target of {}",
            image.error,
            target * 1e-3
        );
    }
    let line = Curve::BSpline(line.curve);
    let domain = line.domain();
    let edge = make_edge(model, line.clone(), domain, tol)?.shape;
    attach_pcurve(
        model,
        &edge,
        PlanarCurve::from(image.curve),
        hosts.ids[f],
        Location::identity(),
        domain,
    )?;
    crate::sheet_curved::same_parameter(model, &edge);
    model.widen(&edge, Tolerance::new(target.max(image.error * 10.0))?)?;
    let split = ogeom_heal::split_face(
        model,
        face,
        face,
        std::slice::from_ref(&edge),
        ogeom_heal::Projection::OnFace,
        tol,
    )
    .map_err(|e| match e {
        OgeomError::Construction(why) => ogeom_core::ogeom_err!(
            Construction,
            "the line along which the ball touches a face does not cut it from boundary to \
             boundary ({why}); the march ends on the face or the line runs along its boundary"
        ),
        other => other,
    })?;
    let cut: Vec<Shape> = split
        .history
        .generated(&edge)
        .iter()
        .filter(|s| model.kind_of(s).is_ok_and(|k| k == ShapeType::Edge))
        .cloned()
        .collect();
    let cut = match cut.as_slice() {
        [] => return Ok(None),
        [one] => one.clone(),
        _ => ogeom_bail!(
            Construction,
            "the line the ball touches a face along crosses it more than once; the stretch the \
             round spans is ambiguous"
        ),
    };
    let (cut_curve, cut_range) = edge_curve(model, &cut, tol)?;
    let mut span = [0.0; 2];
    for (k, t) in [cut_range.0, cut_range.1].into_iter().enumerate() {
        let p = cut_curve.point_at(t, tol)?;
        span[k] = project_on_curve(&line, p, 512, tol)?.parameter;
    }
    let span = (span[0].min(span[1]), span[0].max(span[1]));
    Ok(Some(ContactCut {
        split: split.shape,
        edge: cut,
        line,
        span,
    }))
}

/// The stations over the run `w0..w1` of the guide, with the ball's exact
/// section at each end.
fn trimmed(
    blend: &MarchedBlend,
    hosts: &Hosts,
    guide: &Curve,
    sides: Sides,
    radius: f64,
    (w0, w1): (f64, f64),
    tol: Tolerances,
) -> OgeomResult<MarchedBlend> {
    let [first, second] = &hosts.surfaces;
    let margin = (w1 - w0) * 1e-6;
    let inner: Vec<usize> = (0..blend.len())
        .filter(|&i| blend.along[i] > w0 + margin && blend.along[i] < w1 - margin)
        .collect();
    let mut run = reordered(blend, &inner);
    let near = |i: usize| -> [f64; 4] {
        [
            blend.on_first[i].0,
            blend.on_first[i].1,
            blend.on_second[i].0,
            blend.on_second[i].1,
        ]
    };
    let nearest = |w: f64| -> usize {
        (0..blend.len())
            .min_by(|&i, &j| {
                (blend.along[i] - w)
                    .abs()
                    .total_cmp(&(blend.along[j] - w).abs())
            })
            .unwrap_or(0)
    };
    for (w, front) in [(w0, true), (w1, false)] {
        let x = seat_section(
            first,
            second,
            radius,
            guide,
            sides,
            on_guide(guide, w),
            near(nearest(w)),
            tol,
        )
        .map_err(|_| {
            ogeom_core::ogeom_err!(
                Construction,
                "the ball does not seat between the faces where the round ends"
            )
        })?;
        let p1 = first.point_at(x[0], x[1], tol)?;
        let p2 = second.point_at(x[2], x[3], tol)?;
        let (du, dv) = first.d1_at(x[0], x[1], tol)?;
        let centre = p1 + du.cross(dv).normalized(tol)? * (f64::from(sides.first) * radius);
        let at = if front { 0 } else { run.len() };
        run.spine.insert(at, centre);
        run.touch_first.insert(at, p1);
        run.touch_second.insert(at, p2);
        run.on_first.insert(at, (x[0], x[1]));
        run.on_second.insert(at, (x[2], x[3]));
        run.along.insert(at, w);
    }
    // A station crowding its neighbour puts two of the band's columns in
    // one place; the ends stay, their crowding neighbours go.
    let mut steps: Vec<f64> = run.spine.windows(2).map(|w| w[0].distance(w[1])).collect();
    steps.sort_by(f64::total_cmp);
    let cramped = steps.get(steps.len() / 2).copied().unwrap_or(0.0) * 0.25;
    let mut i = 1;
    while i < run.len() {
        let victim = if i == run.len() - 1 { i - 1 } else { i };
        if victim > 0 && run.spine[i].distance(run.spine[i - 1]) <= cramped && run.len() > 4 {
            run.spine.remove(victim);
            run.touch_first.remove(victim);
            run.touch_second.remove(victim);
            run.on_first.remove(victim);
            run.on_second.remove(victim);
            run.along.remove(victim);
        } else {
            i += 1;
        }
    }
    if run.len() < 4 {
        ogeom_bail!(
            NotDone,
            "the rolling ball's march gave too few stations over the run to fit a round"
        );
    }
    Ok(run)
}

/// Blend two faces of `solid` over a marched seat: the corner between the
/// round and the crease cut off (`behind`, the ball rolling in the
/// material) or filled. Round a closed seat the corner is a ring; over a
/// stretch of the seat it is capped at each end in the ball's section
/// there, as the open marched edge blend caps its run.
pub(crate) fn apply_to_solid(
    model: &mut Model,
    solid: &Shape,
    seat: PairSeat,
    radius: f64,
    behind: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let PairSeat {
        hosts,
        guide,
        blend,
        run,
        ..
    } = seat;
    let [first, second] = &hosts.surfaces;
    let surfaces = [(first, hosts.signs[0]), (second, hosts.signs[1])];
    match run {
        Run::Closed => {
            let range = guide.domain();
            crate::marched::closed_band_wedge(
                model, solid, None, blend, &guide, range, true, surfaces, radius, behind, tol,
            )
        }
        Run::Capped => {
            let (w0, w1) = (blend.along[0], blend.along[blend.len() - 1]);
            let crease = crease_over(&guide, (w0, w1), tol)?;
            let mut blend = blend;
            for t in &mut blend.along {
                *t =
                    project_on_curve(&crease, guide.point_at(on_guide(&guide, *t), tol)?, 64, tol)?
                        .parameter;
            }
            crate::marched::build_open_band(
                model,
                solid,
                None,
                &blend,
                &crease,
                surfaces,
                radius,
                behind,
                [false, false],
                tol,
            )
        }
        Run::Open(_) => ogeom_bail!(
            Invariant,
            "a seat between faces of a solid is taken round a closed seat or capped"
        ),
    }
}

/// A parameter unwrapped round a closed guide, brought back into the
/// guide's domain where the guide is not periodic and so does not wrap
/// it itself.
fn on_guide(guide: &Curve, t: f64) -> f64 {
    let (lo, hi) = guide.domain();
    if guide.is_periodic() || (lo..=hi).contains(&t) {
        return t;
    }
    lo + (t - lo).rem_euclid(hi - lo)
}

/// The crease under a capped round as one spline over the run `w0..w1`
/// of the guide and a little past it, exact for a conic. The run's
/// parameters are unwrapped round a closed guide; a run across the start
/// of a closed guide that is not periodic takes the guide's two stretches
/// either side of the start joined end to end.
fn crease_over(guide: &Curve, (w0, w1): (f64, f64), tol: Tolerances) -> OgeomResult<Curve> {
    let (lo, hi) = guide.domain();
    let period = hi - lo;
    let margin = ((w1 - w0) * 0.02).min((period - (w1 - w0)) * 0.25).max(0.0);
    let (mut from, mut to) = (w0 - margin, w1 + margin);
    if guide.is_periodic() {
        return Ok(Curve::BSpline(guide.to_bspline_over((from, to), tol)?));
    }
    while from >= hi {
        (from, to) = (from - period, to - period);
    }
    while from < lo {
        (from, to) = (from + period, to + period);
    }
    if to <= hi {
        return Ok(Curve::BSpline(guide.to_bspline_over((from, to), tol)?));
    }
    let joined = guide
        .point_at(lo, tol)
        .and_then(|p| guide.point_at(hi, tol).map(|q| p.distance(q)))
        .is_ok_and(|d| d <= tol.confusion() * 10.0);
    if !joined || to - period >= from {
        ogeom_bail!(
            Invariant,
            "a capped run reaches past the ends of a guide it did not go round"
        );
    }
    let before = guide.to_bspline_over((from, hi), tol)?;
    let after = guide.to_bspline_over((lo, to - period), tol)?;
    let whole = ogeom_math::bspline::join(
        &(before.knots().clone(), before.control_points().to_vec()),
        &(after.knots().clone(), after.control_points().to_vec()),
    )?;
    Ok(Curve::BSpline(ogeom_geom::BSplineCurve::rational(
        whole.0, whole.1,
    )?))
}

/// Round the corner between two faces of separate shapes over a marched
/// seat, the ball in front of both; with `trim`, each face cut back to its
/// line of contact and the three sewn into one shell.
pub(crate) fn fillet_faces(
    model: &mut Model,
    a: &Shape,
    b: &Shape,
    radius: f64,
    trim: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let hosts = hosts_of(model, [a, b], radius, tol)?;
    let seat = pair_seat(model, [a, b], &hosts, radius, false, false, tol)?;
    let PairSeat {
        hosts, blend, run, ..
    } = seat;
    let stretches = match run {
        Run::Open(stretches) => stretches,
        Run::Closed => ogeom_bail!(
            Construction,
            "the ball rolls round a closed seat between the two faces; a marched round closing \
             on itself between separate faces is not built"
        ),
        Run::Capped => ogeom_bail!(
            Invariant,
            "a seat between separate faces is taken over their stretches"
        ),
    };
    let fit_target = band_fit_target(tol);
    let band = fit_open_band(&blend, radius, [false, false], tol)?;
    let (w0, w1) = (blend.along[0], blend.along[blend.len() - 1]);
    let marched = Marched {
        blend,
        band,
        fit_target,
    };
    let widened = Tolerance::new(fit_target)?;
    let mut history = History::new();

    if !trim {
        let mut corners = Vec::with_capacity(2);
        let mut rails = Vec::with_capacity(2);
        for f in 0..2 {
            let v0 = make_vertex(model, marched.corner(f, 0)).shape;
            let v1 = make_vertex(model, marched.corner(f, 1)).shape;
            model.widen(&v0, widened)?;
            model.widen(&v1, widened)?;
            let border = marched.border(f, tol)?;
            let domain = border.domain();
            let rail = make_edge_between(model, border, domain, &v0, &v1, tol)?.shape;
            model.widen(&rail, widened)?;
            corners.push([v0, v1]);
            rails.push(rail);
        }
        let round = marched_round(
            model,
            &marched,
            &[rails[0].clone(), rails[1].clone()],
            &[corners[0].clone(), corners[1].clone()],
            true,
            false,
            tol,
        )?;
        history.generate(a, round.clone());
        history.generate(b, round.clone());
        return Ok(Built::new(round, history));
    }

    let mut kept = Vec::with_capacity(2);
    let mut rails = Vec::with_capacity(2);
    let mut corners = Vec::with_capacity(2);
    for (f, stretch) in stretches.iter().enumerate() {
        let (_, _, tolerance) = unplaced_face(model, &stretch.piece)?;
        let ends = [
            stretch.line.point_at(w0, tol)?,
            stretch.line.point_at(w1, tol)?,
        ];
        let (piece, rail, ends) = cut_to_run(
            model,
            &stretch.piece,
            &stretch.edge,
            hosts.ids[f],
            tolerance,
            ends,
            tol,
        )?;
        // The rail's image on the band is fitted; it and its ends hold the
        // band's own fit.
        model.widen(&rail, widened)?;
        for v in &ends {
            model.widen(v, widened)?;
        }
        kept.push(piece);
        rails.push(rail);
        corners.push(ends);
    }
    let round = marched_round(
        model,
        &marched,
        &[rails[0].clone(), rails[1].clone()],
        &[corners[0].clone(), corners[1].clone()],
        false,
        false,
        tol,
    )?;
    let shell = model.add_shell(&[kept[0].clone(), round.clone(), kept[1].clone()])?;
    history.modify(a, kept[0].clone());
    history.modify(b, kept[1].clone());
    history.generate(a, round.clone());
    history.generate(b, round);
    Ok(Built::new(shell, history))
}
