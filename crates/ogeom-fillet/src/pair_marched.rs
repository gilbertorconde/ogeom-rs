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
//! face is cut back to its line of contact and the round sewn between them:
//! over the overlap of the faces' stretches, or, where the ball touches
//! both all the way round a closed seat, as a band closed round the loop
//! between two closed lines of contact.

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
use ogeom_math::{Point, Point2, Vector, solve};
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
    let mut reach = radius * crate::marched::HOST_REACH;
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
            SurfaceGeometry::BSpline(patch) => SurfaceGeometry::BSpline(
                crate::marched::continued_patch(patch, radius * crate::marched::HOST_REACH, tol),
            ),
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
    /// For a seat between faces of a solid, two points beside the middle of
    /// the ball's arc where the solid says which side of the faces the ball
    /// rolls on: one just off the arc toward the crease, where it says
    /// whether the corner is material, and one as far inside the ball,
    /// where it says whether the ball's lune is.
    probes: Option<[Point; 2]>,
}

impl PairSeat {
    pub(crate) const fn probes(&self) -> Option<[Point; 2]> {
        self.probes
    }

    /// Whether the round runs all the way round a closed seat.
    pub(crate) const fn closed(&self) -> bool {
        matches!(self.run, Run::Closed)
    }
}

/// The marched seat of a ball of `radius` touching both faces, on the side
/// each face's normal points to, or behind both with `behind`. `hosts` are
/// the faces' surfaces as [`hosts_of`] carries them on.
///
/// Between faces of `solid` the corner is where the surfaces cross: only
/// that crossing guides the ball, and the round runs all the way round a
/// closed seat or ends where the ball leaves a face or where the solid's
/// edges along the crease end. Between separate faces, where the surfaces
/// do not cross, their offsets' crossing (the ball's own centre line)
/// guides it.
pub(crate) fn pair_seat(
    model: &mut Model,
    faces: [&Shape; 2],
    hosts: &Hosts,
    radius: f64,
    behind: bool,
    solid: Option<&Shape>,
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
        if solid.is_some() {
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
            Ok(Some((blend, run, probes))) => found.push((guide, blend, run, probes)),
            Ok(None) => {}
            Err(e @ (OgeomError::Construction(_) | OgeomError::NotDone(_))) => why = Some(e),
            Err(e) => return Err(e),
        }
    }
    match found.len() {
        1 => {
            let (guide, blend, run, probes) = found.remove(0);
            Ok(PairSeat {
                hosts,
                guide,
                blend,
                run,
                probes,
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
pub(crate) fn reading(radius: f64, tol: Tolerances) -> ogeom_mesh::Deflection {
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
/// does not touch both faces along it. With `whole`, the solid the faces
/// bound, the run is the closed seat or the stretch of it the faces and the
/// solid's edges along the crease give ends to, and the probes beside it
/// are read.
#[allow(clippy::too_many_arguments, reason = "one seat, all its data")]
#[allow(clippy::type_complexity, reason = "the seat's three parts")]
fn seat_along(
    model: &mut Model,
    faces: [&Shape; 2],
    hosts: &Hosts,
    guide: &Curve,
    sides: Sides,
    radius: f64,
    whole: Option<&Shape>,
    tol: Tolerances,
) -> OgeomResult<Option<(MarchedBlend, Run, Option<[Point; 2]>)>> {
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
    let sampled = if whole.is_some() {
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
            let probes = whole
                .is_some()
                .then(|| probes_at(&blend, guide, blend.len() / 2, radius, tol))
                .transpose()?;
            return Ok(Some((blend, Run::Closed, probes)));
        }
        // A loop the faces hold part of: opened where the ball misses a
        // face, and read from there round as one run. Between faces of a
        // solid it is opened where the ball misses every face it misses
        // anywhere, so each such face's line of contact enters and leaves
        // it once along the run.
        (true, Some(at)) => {
            let n = blend.len();
            let at = if whole.is_some() {
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

    if let Some(solid) = whole {
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
            let Some(cut) = contact_cut(model, face, f, &blend, hosts, guide, radius, tol)? else {
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
        // Where the solid's edges along the crease end inside the run, the
        // corner ends there: past it the crease is no corner of the solid
        // (a crease running off the side of a block the ball still rolls
        // over). That end caps in the section plane through the edges' end.
        let ends = crease_edge_ends(model, solid, guide, loops, (w0, w1), tol)?;
        if let Some((w, _)) = ends[0] {
            w0 = w;
        }
        if let Some((w, _)) = ends[1] {
            w1 = w;
        }
        let through = [ends[0].map(|(_, p)| p), ends[1].map(|(_, p)| p)];
        let run = trimmed(&blend, hosts, guide, sides, radius, (w0, w1), through, tol)?;
        let probes = probes_at(&run, guide, run.len() / 2, radius, tol)?;
        return Ok(Some((run, Run::Capped, Some(probes))));
    }

    let mut stretches = Vec::with_capacity(2);
    for (f, face) in faces.iter().enumerate() {
        stretches.push(stretch_on(
            model, face, f, &blend, hosts, guide, radius, tol,
        )?);
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
    let run = trimmed(
        &blend,
        hosts,
        guide,
        sides,
        radius,
        (w0, w1),
        [None; 2],
        tol,
    )?;
    Ok(Some((run, Run::Open(Box::new([a, b])), None)))
}

/// Two points beside the middle of the ball's arc at station `i`: one just
/// off it toward the crease, and one as far inside the ball, toward its
/// centre.
fn probes_at(
    blend: &MarchedBlend,
    guide: &Curve,
    i: usize,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<[Point; 2]> {
    let centre = blend.spine[i];
    let bisector =
        ((blend.touch_first[i] - centre) + (blend.touch_second[i] - centre)).normalized(tol)?;
    let on_arc = centre + bisector * radius;
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
    let off = gap.min(radius) * 0.1;
    Ok([on_arc + toward / gap * off, on_arc - bisector * off])
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
#[allow(clippy::too_many_arguments, reason = "one stretch, all its data")]
fn stretch_on(
    model: &mut Model,
    face: &Shape,
    f: usize,
    blend: &MarchedBlend,
    hosts: &Hosts,
    guide: &Curve,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<FaceStretch>> {
    let surface = &hosts.surfaces[f];
    let Some(ContactCut {
        split,
        edge: cut,
        line,
        span,
    }) = contact_cut(model, face, f, blend, hosts, guide, radius, tol)?
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

/// The ball's line of contact on one face, fitted at the guide's
/// parameters and read against the ball's true contact between the
/// stations it was fitted through.
struct ContactLine {
    /// The guide's parameter at each station the line was fitted through:
    /// the march's, and the sections solved between them.
    along: Vec<f64>,
    /// The touch point's place in the face's chart at each station.
    on: Vec<(f64, f64)>,
    /// The line, parameterized as the guide.
    line: ogeom_geom::BSplineCurve,
    /// The farthest the line stands from the ball's true contact: at its
    /// stations, and at the quarters of the way between them.
    error: f64,
}

/// The ball's line of contact on face `f` (0 the first, 1 the second),
/// fitted through the march's touch points at their guide parameters and
/// held to `target` between them as well as at them.
///
/// A fit through the stations states its error at the stations only, and
/// where the touch point runs fast or unevenly in the guide's parameter the
/// fit can leave the true contact between them by many times its target.
/// So the ball is seated exactly at the quarters of the way between every
/// two neighbouring stations, the same solve the march made at each
/// station, and wherever the line passes farther than `target` from those
/// touch points, or reaches them a guide parameter late or early by more
/// than the round's own fit, the one it misses most joins the stations
/// and the line is fitted again. A line that cannot be held to `target`
/// so is stated as far off as it was read, the closest of the fits tried.
#[allow(clippy::too_many_arguments, reason = "one line, all its data")]
fn contact_line(
    blend: &MarchedBlend,
    f: usize,
    hosts: &Hosts,
    guide: &Curve,
    radius: f64,
    target: f64,
    tol: Tolerances,
) -> OgeomResult<ContactLine> {
    /// The most rounds of sections added between stations.
    const ROUNDS: usize = 8;
    let [first, second] = &hosts.surfaces;
    let surface = &hosts.surfaces[f];
    // The round's band is fitted to this, and its corners on the line are
    // widened to it: where the line stands at a guide parameter may stray
    // from the ball's section there by no more.
    let usable = band_fit_target(tol);
    let mut along = blend.along.clone();
    let mut touch = if f == 0 {
        blend.touch_first.clone()
    } else {
        blend.touch_second.clone()
    };
    let mut seats: Vec<[f64; 4]> = blend
        .on_first
        .iter()
        .zip(&blend.on_second)
        .map(|(a, b)| [a.0, a.1, b.0, b.1])
        .collect();
    let mut best: Option<(ContactLine, f64)> = None;
    for round in 0..=ROUNDS {
        let fitted = ogeom_geom::fit::fit_points_at(&along, &touch, 3, target, tol)?;
        if !fitted.met && round == 0 {
            ogeom_bail!(
                NotDone,
                "the ball's line of contact on a face reached {} against a target of {target}",
                fitted.error
            );
        }
        let line = fitted.curve;
        let mut error = fitted.error;
        let mut slip = fitted.error;
        let mut added: Vec<(usize, f64, [f64; 4], Point)> = Vec::new();
        for i in 0..along.len() - 1 {
            // Seated from either neighbour; a section that lands off the
            // stretch between them is another seating, not this one.
            let reach = touch[i].distance(touch[i + 1]).mul_add(3.0, target);
            // Where the line passes nearest a touch point, sought over the
            // stations either side.
            let (from, to) = (i.saturating_sub(1), (i + 2).min(along.len() - 1));
            let around = (along[from], along[to]);
            // The ball's contact anywhere over those stations, seated from
            // the nearest of them and kept to the stretch they span.
            let span = touch[from].distance(touch[to]);
            let contact_at = |t: f64| -> Option<Point> {
                let k = (from..=to)
                    .min_by(|&j, &k| (along[j] - t).abs().total_cmp(&(along[k] - t).abs()))?;
                let x = seat_section(
                    first,
                    second,
                    radius,
                    guide,
                    blend.sides,
                    on_guide(guide, t),
                    seats[k],
                    tol,
                )
                .ok()?;
                let p = surface.point_at(x[2 * f], x[2 * f + 1], tol).ok()?;
                (p.distance(touch[from]) + p.distance(touch[to]) <= span.mul_add(3.0, target))
                    .then_some(p)
            };
            let mut worst: Option<(f64, f64, [f64; 4], Point)> = None;
            for share in [0.25, 0.5, 0.75] {
                let w = (along[i + 1] - along[i]).mul_add(share, along[i]);
                if w <= along[i] || w >= along[i + 1] {
                    continue;
                }
                let mut found = None;
                for near in [seats[i], seats[i + 1]] {
                    let Ok(x) = seat_section(
                        first,
                        second,
                        radius,
                        guide,
                        blend.sides,
                        on_guide(guide, w),
                        near,
                        tol,
                    ) else {
                        continue;
                    };
                    let truth = surface.point_at(x[2 * f], x[2 * f + 1], tol)?;
                    if truth.distance(touch[i]) + truth.distance(touch[i + 1]) <= reach {
                        found = Some(([x[0], x[1], x[2], x[3]], truth));
                        break;
                    }
                }
                let Some((seat, truth)) = found else {
                    ogeom_bail!(
                        NotDone,
                        "the ball does not seat between two stations of its march, so its \
                         line of contact there cannot be read"
                    );
                };
                let late = line.point_at(w, tol)?.distance(truth);
                // Off by more than the target at the same parameter, the
                // line may still run along the contact, only a parameter
                // ahead or behind: how far it stands off is then read both
                // ways, the touch point from the line and the line's point
                // from the ball's contact anywhere round the stations.
                let off = if late <= target {
                    late
                } else {
                    let criteria = solve::Criteria {
                        residual: 0.0,
                        step: tol.parametric(),
                        max_iterations: 100,
                    };
                    let to_line = |t: f64| {
                        line.point_at(t, tol)
                            .map_or(f64::INFINITY, |p| p.square_distance(truth))
                    };
                    let nearest = solve::minimize(to_line, around.0, around.1, criteria)?;
                    let point = line.point_at(w, tol)?;
                    let to_contact =
                        |t: f64| contact_at(t).map_or(f64::INFINITY, |p| p.square_distance(point));
                    let touching = solve::minimize(to_contact, around.0, around.1, criteria)?;
                    to_line(nearest.value)
                        .max(to_contact(touching.value))
                        .sqrt()
                        .min(late)
                };
                error = error.max(off);
                slip = slip.max(late);
                let miss = (off / target).max(late / usable);
                if miss > 1.0 && worst.as_ref().is_none_or(|held| miss > held.0) {
                    worst = Some((miss, w, seat, truth));
                }
            }
            if let Some((_, w, seat, truth)) = worst {
                added.push((i + 1, w, seat, truth));
            }
        }
        let settled = fitted.met && added.is_empty();
        let closer = best.as_ref().is_none_or(|(held, held_slip)| {
            (slip > usable, error) < (*held_slip > usable, held.error)
        });
        if closer {
            let on = seats.iter().map(|s| (s[2 * f], s[2 * f + 1])).collect();
            best = Some((
                ContactLine {
                    along: along.clone(),
                    on,
                    line,
                    error,
                },
                slip,
            ));
        }
        // A fit that misses its target through the stations it was given
        // gains nothing from more of them.
        if settled || !fitted.met {
            break;
        }
        for (at, w, seat, truth) in added.into_iter().rev() {
            along.insert(at, w);
            touch.insert(at, truth);
            seats.insert(at, seat);
        }
    }
    let Some((best, slip)) = best else {
        ogeom_bail!(Invariant, "a line of contact is fitted at least once");
    };
    if best.error.max(slip) > usable {
        ogeom_bail!(
            NotDone,
            "the ball's line of contact on a face strays {} from the ball's true contact \
             between the march's stations, past the round's own fit of {usable}",
            best.error.max(slip)
        );
    }
    Ok(best)
}

/// The face cut along the ball's line of contact: the touch points fitted
/// at the guide's parameters with their places in the face's chart.
/// `None` where the line misses the face.
#[allow(clippy::too_many_arguments, reason = "one cut, all its data")]
fn contact_cut(
    model: &mut Model,
    face: &Shape,
    f: usize,
    blend: &MarchedBlend,
    hosts: &Hosts,
    guide: &Curve,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<ContactCut>> {
    let surface = &hosts.surfaces[f];
    let target = tol.confusion() * 10.0;
    let ContactLine {
        along,
        on,
        line,
        error,
    } = contact_line(blend, f, hosts, guide, radius, target, tol)?;
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
    let image = ogeom_geom::fit::fit_points_2d_at(&along, &chart, 3, target * 1e-3, tol)?;
    if !image.met {
        ogeom_bail!(
            NotDone,
            "the chart image of the ball's line of contact on a face reached {} against a \
             target of {}",
            image.error,
            target * 1e-3
        );
    }
    let line = Curve::BSpline(line);
    let domain = line.domain();
    let image = PlanarCurve::from(image.curve);
    // The line and its image were fitted apart through the same stations;
    // the edge states how far the image, lifted, stands from the line
    // anywhere along it, how far the line stands from the ball's true
    // contact between the stations, and never less than the fit's target.
    let gap =
        ogeom_algo::pcurve_fit::lifted_gap((&line, domain), (&image, domain), surface, false, tol)?;
    let edge = make_edge(model, line.clone(), domain, tol)?.shape;
    attach_pcurve(
        model,
        &edge,
        image,
        hosts.ids[f],
        Location::identity(),
        domain,
    )?;
    crate::sheet_curved::same_parameter(model, &edge);
    model.widen(
        &edge,
        Tolerance::new(target.max(gap.max(error) + tol.confusion()))?,
    )?;
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
/// section at each end. An end given a point in `through` takes the
/// section whose arc plane holds that point, a whisker off its parameter.
#[allow(clippy::too_many_arguments, reason = "one run, all its data")]
fn trimmed(
    blend: &MarchedBlend,
    hosts: &Hosts,
    guide: &Curve,
    sides: Sides,
    radius: f64,
    (w0, w1): (f64, f64),
    through: [Option<Point>; 2],
    tol: Tolerances,
) -> OgeomResult<MarchedBlend> {
    let [first, second] = &hosts.surfaces;
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
    let mut ends = Vec::with_capacity(2);
    for (w, point) in [(w0, through[0]), (w1, through[1])] {
        let seed = near(nearest(w));
        let section = |w: f64| -> OgeomResult<Section> {
            let x = seat_section(
                first,
                second,
                radius,
                guide,
                sides,
                on_guide(guide, w),
                seed,
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
            Ok(Section {
                w,
                x,
                p1,
                p2,
                centre,
            })
        };
        let found = match point {
            Some(point) => section_through(&section, w, (w1 - w0) * 1e-4, point, tol)?,
            None => section(w)?,
        };
        ends.push(found);
    }
    let (start, end) = (ends.remove(0), ends.remove(0));
    if end.w <= start.w {
        ogeom_bail!(
            Construction,
            "the round's end sections cross; nothing to blend"
        );
    }
    let margin = (end.w - start.w) * 1e-6;
    let inner: Vec<usize> = (0..blend.len())
        .filter(|&i| blend.along[i] > start.w + margin && blend.along[i] < end.w - margin)
        .collect();
    let mut run = reordered(blend, &inner);
    for (found, front) in [(start, true), (end, false)] {
        let at = if front { 0 } else { run.len() };
        run.spine.insert(at, found.centre);
        run.touch_first.insert(at, found.p1);
        run.touch_second.insert(at, found.p2);
        run.on_first.insert(at, (found.x[0], found.x[1]));
        run.on_second.insert(at, (found.x[2], found.x[3]));
        run.along.insert(at, found.w);
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

/// The ball's exact section at one guide parameter: the seat solve, the
/// touch points and the ball's centre.
struct Section {
    w: f64,
    x: [f64; 5],
    p1: Point,
    p2: Point,
    centre: Point,
}

/// The section near `w` whose arc plane holds `point`, by secant steps on
/// the guide parameter from `w` and `w + step`.
fn section_through(
    section: &impl Fn(f64) -> OgeomResult<Section>,
    w: f64,
    step: f64,
    point: Point,
    tol: Tolerances,
) -> OgeomResult<Section> {
    // The point's signed distance from a section's arc plane.
    let off = |s: &Section| -> f64 {
        let n = (s.p1 - s.centre).cross(s.p2 - s.centre);
        let m = n.magnitude();
        if m <= f64::MIN_POSITIVE {
            0.0
        } else {
            (point - s.centre).dot(n / m)
        }
    };
    let mut found = section(w)?;
    let mut f = off(&found);
    let other = section(w + step)?;
    let mut previous = (other.w, off(&other));
    for _ in 0..30 {
        if f.abs() <= tol.confusion() * 1e-2 {
            return Ok(found);
        }
        let slope = (previous.1 - f) / (previous.0 - found.w);
        if !slope.is_finite() || slope.abs() <= f64::MIN_POSITIVE {
            break;
        }
        previous = (found.w, f);
        found = section(found.w - f / slope)?;
        f = off(&found);
    }
    if f.abs() <= tol.confusion() {
        return Ok(found);
    }
    ogeom_bail!(
        Construction,
        "no section of the ball holds the end of the solid's edges along the crease"
    )
}

/// Where the solid's edges along the crease end inside the run `w0..w1` of
/// the guide (unwrapped round a guide that `loops`): for each end of the
/// run, the guide parameter and point where the stretch those edges cover
/// round the run's middle ends before the run does. `None` at an end the
/// edges reach past, and at both where no edge runs along the crease at
/// the run's middle.
pub(crate) fn crease_edge_ends(
    model: &Model,
    solid: &Shape,
    guide: &Curve,
    loops: bool,
    (w0, w1): (f64, f64),
    tol: Tolerances,
) -> OgeomResult<[Option<(f64, Point)>; 2]> {
    let (lo, hi) = guide.domain();
    let period = hi - lo;
    let middle = f64::midpoint(w0, w1);
    let near_to = |t: f64, near: f64| -> f64 {
        if loops {
            t + ((near - t) / period).round() * period
        } else {
            t
        }
    };
    let mut covered: Vec<(f64, f64, Point, Point)> = Vec::new();
    for edge in explore_unique(model, solid, ShapeType::Edge)? {
        let (curve, (a, b)) = match edge_curve(model, &edge, tol) {
            Ok(found) => found,
            Err(OgeomError::Construction(_)) => continue,
            Err(e) => return Err(e),
        };
        if b - a <= tol.parametric() {
            continue;
        }
        let reach = model
            .node(&edge)
            .and_then(|n| n.data().as_edge().map(|d| d.tolerance.get()))
            .unwrap_or(0.0_f64)
            .mul_add(2.0, tol.confusion() * 100.0);
        // Five samples on the guide, unwrapped outward from the middle one.
        let mut params = [0.0; 5];
        let mut along = true;
        for (k, slot) in params.iter_mut().enumerate() {
            #[allow(clippy::cast_precision_loss, reason = "a small index")]
            let t = a + (b - a) * (k as f64) / 4.0;
            let found = project_on_curve(guide, curve.point_at(t, tol)?, 64, tol)?;
            if found.distance > reach {
                along = false;
                break;
            }
            *slot = found.parameter;
        }
        if !along {
            continue;
        }
        params[2] = near_to(params[2], middle);
        for k in (0..2).rev() {
            params[k] = near_to(params[k], params[k + 1]);
        }
        for k in 3..5 {
            params[k] = near_to(params[k], params[k - 1]);
        }
        let (start, end) = (curve.point_at(a, tol)?, curve.point_at(b, tol)?);
        let (from, to, p, q) = if params[0] <= params[4] {
            (params[0], params[4], start, end)
        } else {
            (params[4], params[0], end, start)
        };
        // Round a loop the same edge stands a period either way too.
        let shifts: &[f64] = if loops { &[-1.0, 0.0, 1.0] } else { &[0.0] };
        for &k in shifts {
            covered.push((from + k * period, to + k * period, p, q));
        }
    }
    covered.sort_by(|x, y| x.0.total_cmp(&y.0));
    let gap = (w1 - w0) * 1e-6 + tol.parametric();
    let mut merged: Vec<(f64, f64, Point, Point)> = Vec::new();
    for c in covered {
        match merged.last_mut() {
            Some(last) if c.0 <= last.1 + gap => {
                if c.1 > last.1 {
                    last.1 = c.1;
                    last.3 = c.3;
                }
            }
            _ => merged.push(c),
        }
    }
    let Some(&(from, to, p, q)) = merged.iter().find(|c| c.0 <= middle && middle <= c.1) else {
        return Ok([None, None]);
    };
    let margin = (w1 - w0) * 1e-6;
    Ok([
        (from > w0 + margin).then_some((from, p)),
        (to < w1 - margin).then_some((to, q)),
    ])
}

/// A point of a face of `solid` other than the two `faces` inside the ball
/// of `radius` at some station of a marched seat's run, as
/// [`face_inside_balls`](crate::between::face_inside_balls) finds it; `None`
/// where none is.
///
/// The ball is seated exactly between the march's stations, from the
/// nearest of them. A closed run goes round its loop and has no ends. A
/// capped run ends at its first and last stations, and a face holding a
/// contact or the crease there is a face it runs off.
pub(crate) fn face_inside(
    model: &mut Model,
    solid: &Shape,
    faces: [&Shape; 2],
    seat: &PairSeat,
    radius: f64,
    tol: Tolerances,
) -> OgeomResult<Option<Point>> {
    let PairSeat {
        hosts,
        guide,
        blend,
        run,
        ..
    } = seat;
    let [first, second] = &hosts.surfaces;
    let n = blend.len();
    let (lo, hi) = guide.domain();
    let period = hi - lo;
    // The stations' guide parameters in order along the run, unwrapped
    // round a loop, a closed run's last station followed by its first.
    let mut stations = blend.along.clone();
    let closed = matches!(run, Run::Closed);
    if closed {
        for i in 1..n {
            while stations[i] - stations[i - 1] > period / 2.0 {
                stations[i] -= period;
            }
            while stations[i - 1] - stations[i] > period / 2.0 {
                stations[i] += period;
            }
        }
        let back = if stations[n - 1] >= stations[0] {
            period
        } else {
            -period
        };
        stations.push(stations[0] + back);
    }
    let index = |k: usize| k % n;
    // How far the centre moves between neighbouring stations at most: a
    // seating farther off than three times that from its seed is another
    // seating.
    let mut step = 0.0_f64;
    for k in 0..stations.len() - 1 {
        step = step.max(blend.spine[index(k)].distance(blend.spine[index(k + 1)]));
    }
    let reach = step.mul_add(3.0, tol.confusion());
    let params = stations.clone();
    let centre = |t: f64| -> OgeomResult<Point> {
        let k = (0..params.len())
            .min_by(|&j, &k| (params[j] - t).abs().total_cmp(&(params[k] - t).abs()))
            .unwrap_or(0);
        let i = index(k);
        if (params[k] - t).abs() <= tol.parametric() {
            return Ok(blend.spine[i]);
        }
        let near = [
            blend.on_first[i].0,
            blend.on_first[i].1,
            blend.on_second[i].0,
            blend.on_second[i].1,
        ];
        let x = seat_section(
            first,
            second,
            radius,
            guide,
            blend.sides,
            on_guide(guide, t),
            near,
            tol,
        )?;
        let p = first.point_at(x[0], x[1], tol)?;
        let Some(normal) = crate::march::unit_normal(first, x[0], x[1], tol) else {
            ogeom_bail!(
                NotDone,
                "the ball's seat between two stations has no normal"
            );
        };
        let c = p + normal * (f64::from(blend.sides.first) * radius);
        if c.distance(blend.spine[i]) > reach {
            ogeom_bail!(
                NotDone,
                "the ball does not seat between two stations of its march, so whether it \
                 clears the solid's other faces there cannot be read"
            );
        }
        Ok(c)
    };
    let (lines, ends) = if closed {
        (vec![(guide.clone(), (lo, hi))], Vec::new())
    } else {
        let (w0, w1) = (stations[0], stations[n - 1]);
        let crease = crease_over(guide, (w0, w1), tol)?;
        let range = crease.domain();
        let mut ends = Vec::with_capacity(6);
        for i in [0, n - 1] {
            ends.push(blend.touch_first[i]);
            ends.push(blend.touch_second[i]);
            ends.push(guide.point_at(on_guide(guide, blend.along[i]), tol)?);
        }
        (vec![(crease, range)], ends)
    };
    let rolling = crate::between::Rolling {
        radius,
        stations,
        centre: &centre,
        lines,
        ends,
        lunes: None,
    };
    crate::between::face_inside_balls(model, solid, faces, &rolling, tol)
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
    let seat = pair_seat(model, [a, b], &hosts, radius, false, None, tol)?;
    let PairSeat {
        hosts, blend, run, ..
    } = seat;
    let stretches = match run {
        Run::Open(stretches) => stretches,
        Run::Closed => return closed_round(model, [a, b], &hosts, blend, radius, trim, tol),
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

/// Round the corner between two faces of separate shapes where the ball
/// rolls round a closed seat touching both all the way: the band fitted
/// closed round the loop of the ball's arcs, bounded by its two borders,
/// one along each line of contact, and a seam across it at the loop's
/// first station. With `trim`, each face is cut along its border and the
/// piece away from the round kept, and the three are one shell.
#[allow(clippy::too_many_arguments, reason = "one round, all its data")]
fn closed_round(
    model: &mut Model,
    faces: [&Shape; 2],
    hosts: &Hosts,
    mut blend: MarchedBlend,
    radius: f64,
    trim: bool,
    tol: Tolerances,
) -> OgeomResult<Built> {
    crate::marched::owned_once(&mut blend, tol);
    let fit_target = band_fit_target(tol);
    // Aimed at half the target, as the open band is: the borders are the
    // rails the faces are cut along, and the band's normal there, which
    // the faces' must match, strays by the fit over the arc's sampling.
    let (fitted, _) = crate::marched::fit_closed_band(&blend, radius, fit_target * 0.5, tol)?;
    if fitted.error > fit_target {
        ogeom_bail!(
            NotDone,
            "the blend surface reached {} against a target of {fit_target}",
            fitted.error
        );
    }
    let band = fitted.curve;
    let marched = Marched {
        blend,
        band,
        fit_target,
    };
    let band: SurfaceGeometry = marched.band.clone().into();
    let (u_dom, v_dom) = band.domain();
    let band_id = model.geometry_mut().add_surface(band.clone());
    let widened = Tolerance::new(fit_target)?;
    let mut rails = Vec::with_capacity(2);
    for f in 0..2 {
        let rail = make_edge(model, marched.border(f, tol)?, v_dom, tol)?.shape;
        model.widen(&rail, widened)?;
        rails.push(rail);
    }

    let mut history = History::new();
    let mut kept = Vec::with_capacity(2);
    let mut walked = Vec::with_capacity(2);
    for (f, face) in faces.iter().enumerate() {
        if !trim {
            walked.push(vec![rails[f].clone()]);
            continue;
        }
        let (piece, pieces) = cut_along_rail(model, face, f, &rails[f], hosts, &marched, tol)?;
        kept.push(piece);
        walked.push(pieces);
    }
    for (f, pieces) in walked.iter().enumerate() {
        let u = if f == 0 { u_dom.0 } else { u_dom.1 };
        for piece in pieces {
            let (_, range) = edge_curve(model, piece, tol)?;
            attach_pcurve(
                model,
                &piece.oriented(Orientation::Forward),
                column(u, v_dom)?,
                band_id,
                Location::identity(),
                range,
            )?;
        }
    }
    let anchors = [
        ogeom_algo::edge_vertices(model, &walked[0][0])?,
        ogeom_algo::edge_vertices(model, &walked[1][0])?,
    ];
    let [Some((anchor0, _)), Some((anchor1, _))] = anchors else {
        ogeom_bail!(Construction, "a rail of the round has no vertex");
    };
    let seam_curve = {
        let grid = marched.band.grid();
        let nv = grid.v_count();
        let control: Vec<Point> = (0..grid.u_count())
            .map(|i| grid.points()[i * nv].point())
            .collect();
        Curve::BSpline(ogeom_geom::BSplineCurve::new(
            marched.band.u_knots().clone(),
            control,
            tol,
        )?)
    };
    let seam = make_edge_between(model, seam_curve, u_dom, &anchor0, &anchor1, tol)?.shape;
    model.widen(&seam, widened)?;
    let row = |v: f64| -> OgeomResult<PlanarCurve> {
        Ok(ogeom_geom::Line2d::over(
            ogeom_math::Axis2::new(Point2::new(0.0, v), ogeom_math::Direction2::X),
            u_dom.0 - 1.0,
            u_dom.1 + 1.0,
        )?
        .into())
    };
    ogeom_algo::attach_seam(
        model,
        &seam,
        row(v_dom.0)?,
        row(v_dom.1)?,
        band_id,
        Location::identity(),
        u_dom,
    )?;
    let ring: Vec<Shape> = [
        vec![seam.clone()],
        walked[1].clone(),
        vec![seam.reversed()],
        walked[0].iter().rev().map(Shape::reversed).collect(),
    ]
    .concat();
    let wire = ogeom_algo::make_wire(model, &ring, tol)?.shape;
    let face = ogeom_algo::make_face_on(model, band_id, std::slice::from_ref(&wire), tol)?.shape;
    let spine = marched.blend.spine.clone();
    let round = crate::sheet_curved::oriented_round(
        model,
        face,
        false,
        |q| crate::marched::touching_ball(&spine, q),
        tol,
    )?;
    for face in faces {
        history.generate(face, round.clone());
    }
    if !trim {
        return Ok(Built::new(round, history));
    }
    let shell = model.add_shell(&[kept[0].clone(), round.clone(), kept[1].clone()])?;
    history.modify(faces[0], kept[0].clone());
    history.modify(faces[1], kept[1].clone());
    Ok(Built::new(shell, history))
}

/// The line `v` of a band's chart at its `u`, over the band's `v` window
/// and a little past it.
fn column(u: f64, (v0, v1): (f64, f64)) -> OgeomResult<PlanarCurve> {
    Ok(ogeom_geom::Line2d::over(
        ogeom_math::Axis2::new(Point2::new(u, 0.0), ogeom_math::Direction2::Y),
        v0 - 1.0,
        v1 + 1.0,
    )?
    .into())
}

/// Face `f` cut along the closed `rail` the round's band borders it with,
/// and the piece away from the round: back come the piece and the edges
/// the cut left along the rail, in order round it from the rail's own
/// vertex.
fn cut_along_rail(
    model: &mut Model,
    face: &Shape,
    f: usize,
    rail: &Shape,
    hosts: &Hosts,
    marched: &Marched,
    tol: Tolerances,
) -> OgeomResult<(Shape, Vec<Shape>)> {
    let split = ogeom_heal::split_face(
        model,
        face,
        face,
        std::slice::from_ref(rail),
        ogeom_heal::Projection::OnFace,
        tol,
    )
    .map_err(|e| match e {
        OgeomError::Construction(why) => ogeom_core::ogeom_err!(
            Construction,
            "cutting a face along the line the ball touches it on fails ({why}); the line \
             crosses a seam of the face its chart does not wrap round, or leaves the face"
        ),
        other => other,
    })?;
    let cut: Vec<Shape> = split
        .history
        .generated(rail)
        .iter()
        .filter(|s| model.kind_of(s).is_ok_and(|k| k == ShapeType::Edge))
        .cloned()
        .collect();
    // The cut's edges stand on the rail's own curve: put in order along it
    // from its start, where the round's seam meets it.
    let (_, whole) = edge_curve(model, rail, tol)?;
    let mut cut = cut
        .into_iter()
        .map(|e| {
            let e = e.oriented(Orientation::Forward);
            edge_curve(model, &e, tol).map(|(_, range)| (range, e))
        })
        .collect::<OgeomResult<Vec<_>>>()?;
    cut.sort_by(|x, y| x.0.0.total_cmp(&y.0.0));
    let gap = tol.parametric() * 1e3;
    let chained = !cut.is_empty()
        && (cut[0].0.0 - whole.0).abs() <= gap
        && (cut[cut.len() - 1].0.1 - whole.1).abs() <= gap
        && cut.windows(2).all(|w| (w[1].0.0 - w[0].0.1).abs() <= gap);
    if !chained {
        ogeom_bail!(
            Construction,
            "the line the ball touches a face along does not stay on the face all the way round; \
             the cut along it leaves no closed rail"
        );
    }
    let cut: Vec<Shape> = cut.into_iter().map(|(_, e)| e).collect();
    let blend = &marched.blend;
    let surface = &hosts.surfaces[f];
    let i = blend.len() / 2;
    let (touch, other) = if f == 0 {
        (&blend.touch_first, &blend.touch_second)
    } else {
        (&blend.touch_second, &blend.touch_first)
    };
    let (centre, at, to) = (blend.spine[i], touch[i], other[i]);
    let (a, b) = (at - centre, to - centre);
    let away: Vector = -(b - a * (a.dot(b) / a.dot(a))).normalized(tol)?;
    let pieces = explore_unique(model, &split.shape, ShapeType::Face)?;
    // The pieces cover the face once: a cut the face's chart cannot hold
    // (a line round a closed spline face, across a seam its chart does not
    // wrap) parts it into pieces that do not. Their areas, integrated on
    // the exact surface, are held to the face's within the most a cut
    // standing within its tolerance of the line could move across it.
    let (rail_curve, rail_range) = edge_curve(model, rail, tol)?;
    let length = ogeom_algo::curve_length(&rail_curve, rail_range, tol)?;
    let slack = model
        .node(rail)
        .and_then(|n| n.data().as_edge().map(|d| d.tolerance.get()))
        .unwrap_or(marched.fit_target)
        * length;
    let area = |model: &Model, face: &Shape| -> Option<f64> {
        ogeom_algo::surface_properties(model, face, ogeom_mesh::Deflection::default(), tol)
            .ok()
            .filter(|found| found.deflection == 0.0)
            .map(|found| found.mass)
    };
    let whole = area(model, face);
    let parts: Option<f64> = pieces.iter().map(|p| area(model, p)).sum();
    if !matches!((whole, parts), (Some(w), Some(p)) if (w - p).abs() <= slack) {
        ogeom_bail!(
            Construction,
            "cutting a face along the line the ball touches it on does not part it into \
             pieces covering it once; the line crosses a seam of the face its chart does not \
             wrap round"
        );
    }
    for scale in [1e-3, 1e-2, 5e-2] {
        let eps = a.magnitude() * scale;
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
        return closed_up(model, piece, cut, tol);
    }
    ogeom_bail!(
        Construction,
        "no piece of a face cut along its line of contact stands away from the round"
    )
}

/// The kept `piece` with the edges `cut` along a closed line, in order
/// round it, joined end to end. A cut across a periodic face's seam ends
/// on a second vertex beside the one it started from, with a sliver of the
/// seam between the two: the second vertex is merged into the first, the
/// sliver dropped, and every edge of the piece at the second vertex moved
/// onto the first. Back come the piece and the cut's edges as it holds
/// them.
fn closed_up(
    model: &mut Model,
    piece: &Shape,
    cut: Vec<Shape>,
    tol: Tolerances,
) -> OgeomResult<(Shape, Vec<Shape>)> {
    let ends = |model: &Model, e: &Shape| -> OgeomResult<(Shape, Shape)> {
        ogeom_algo::edge_vertices(model, e)?
            .ok_or_else(|| ogeom_core::ogeom_err!(Construction, "a cut edge has no vertices"))
    };
    let tolerance_of = |model: &Model, v: &Shape| -> f64 {
        model
            .node(v)
            .and_then(|n| n.data().as_vertex().map(|d| d.tolerance.get()))
            .unwrap_or(0.0)
    };
    // Each vertex a cut edge ends on that the next does not start from,
    // paired with the one it does start from.
    let mut moves: Vec<(Shape, Shape)> = Vec::new();
    for k in 0..cut.len() {
        let (_, end) = ends(model, &cut[k])?;
        let (start, _) = ends(model, &cut[(k + 1) % cut.len()])?;
        if end.is_same(&start) {
            continue;
        }
        let (p, q) = (
            crate::sheet_curved::vertex_point(model, &end)?,
            crate::sheet_curved::vertex_point(model, &start)?,
        );
        if p.distance(q) > tolerance_of(model, &end).max(tolerance_of(model, &start)) {
            ogeom_bail!(
                Construction,
                "the cut along the line the ball touches a face on does not close; its ends \
                 stand {} apart",
                p.distance(q)
            );
        }
        moves.push((end, start));
    }
    if moves.is_empty() {
        return Ok((piece.clone(), cut));
    }
    let (surface_id, _, tolerance) = unplaced_face(model, piece)?;
    let moved = |v: &Shape| -> Shape {
        moves
            .iter()
            .find(|(gone, _)| gone.is_same(v))
            .map_or_else(|| v.clone(), |(_, kept)| kept.clone())
    };
    let paired = |v: &Shape| moves.iter().any(|(g, k)| g.is_same(v) || k.is_same(v));
    let mut swap: Vec<(Shape, Vec<Shape>)> = Vec::new();
    let mut restated = cut.clone();
    for edge in explore_unique(model, piece, ShapeType::Edge)? {
        let forward = edge.oriented(Orientation::Forward);
        let (a, b) = ends(model, &forward)?;
        let on_cut = cut.iter().position(|c| c.is_same(&forward));
        if on_cut.is_none()
            && !a.is_same(&b)
            && paired(&a)
            && paired(&b)
            && moved(&a).is_same(&moved(&b))
        {
            // The sliver between a pair.
            swap.push((forward, Vec::new()));
            continue;
        }
        let (from, to) = (moved(&a), moved(&b));
        if from.is_same(&a) && to.is_same(&b) {
            continue;
        }
        let (_, range) = edge_curve(model, &forward, tol)?;
        let new =
            crate::sheet_curved::piece_of(model, &forward, surface_id, range, &from, &to, tol)?;
        if let Some(k) = on_cut {
            restated[k] = new.clone();
        }
        swap.push((forward, vec![new]));
    }
    for (gone, kept) in &moves {
        let p = crate::sheet_curved::vertex_point(model, gone)?;
        let reach =
            p.distance(crate::sheet_curved::vertex_point(model, kept)?) + tolerance_of(model, gone);
        model.widen(kept, Tolerance::new(reach.max(tolerance_of(model, kept)))?)?;
    }
    let rebuilt =
        crate::sheet_curved::rebuilt_face(model, piece, surface_id, tolerance, &swap, tol)?;
    Ok((rebuilt, restated))
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
mod tests {
    use super::*;
    use ogeom_math::{ControlGrid, Direction, KnotVector, Plane};

    const T: Tolerances = Tolerances::millimetres();

    /// A floor z = 0 and a cubic patch rising through it over y from -8 to
    /// 8, leaning `lean` degrees off the floor across y, rippled by
    /// `ripple` and bowed by `bow`, each with the sign that sets the ball
    /// in the acute corner between them.
    fn floor_and_patch(lean: &dyn Fn(f64) -> f64, ripple: f64, rows: usize, bow: f64) -> Hosts {
        let mut points = Vec::with_capacity(4 * rows);
        #[allow(clippy::cast_precision_loss, reason = "a small count")]
        let last = (rows - 1) as f64;
        for i in 0..4 {
            let t = 11.0_f64.mul_add(f64::from(i) / 3.0, -3.0);
            for j in 0..rows {
                #[allow(clippy::cast_precision_loss, reason = "a small count")]
                let y = 16.0_f64.mul_add(j as f64 / last, -8.0);
                let (s, c) = lean(y).to_radians().sin_cos();
                let r = if j % 2 == 0 { ripple } else { -ripple };
                let b = bow * (y / 8.0).powi(2) + 0.3 * (t / 8.0).powi(2) + r;
                points.push(Point::new(t.mul_add(c, b), y, t * s));
            }
        }
        let cubic = KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3).unwrap();
        let mut knots = vec![0.0; 4];
        for k in 1..rows - 3 {
            #[allow(clippy::cast_precision_loss, reason = "a small count")]
            knots.push(k as f64 / (rows - 3) as f64);
        }
        knots.extend([1.0; 4]);
        let along = KnotVector::new(knots, 3).unwrap();
        let patch = ogeom_geom::BSplineSurface::new(
            cubic,
            along,
            &ControlGrid::new(points, 4, rows).unwrap(),
            T,
        )
        .unwrap();
        let floor = PlaneSurface::over(
            Plane::through(Point::ORIGIN, Direction::Z),
            (-40.0, 40.0),
            (-40.0, 40.0),
        )
        .unwrap();
        let (s, c) = lean(0.0).to_radians().sin_cos();
        let (du, dv) = patch.d1_at(0.5, 0.5, T).unwrap();
        let sign = if du.cross(dv).dot(Vector::new(s, 0.0, -c)) > 0.0 {
            1.0
        } else {
            -1.0
        };
        let mut model = Model::new();
        let surfaces: [SurfaceGeometry; 2] = [floor.into(), patch.into()];
        let ids = [
            model.geometry_mut().add_surface(surfaces[0].clone()),
            model.geometry_mut().add_surface(surfaces[1].clone()),
        ];
        Hosts {
            ids,
            surfaces,
            signs: [1.0, sign],
        }
    }

    /// The marches of a ball of `radius` along where the supports' offsets
    /// cross, each in order along its guide.
    fn marches(hosts: &Hosts, radius: f64) -> Vec<(Curve, MarchedBlend)> {
        #[allow(clippy::cast_possible_truncation, reason = "a sign")]
        let sides = Sides {
            first: 1,
            second: hosts.signs[1] as i8,
        };
        let [first, second] = &hosts.surfaces;
        let offset = |surface: &SurfaceGeometry, side: i8| -> SurfaceGeometry {
            SurfaceGeometry::Offset(Box::new(
                OffsetSurface::new(surface.clone(), f64::from(side) * radius).unwrap(),
            ))
        };
        let SurfaceIntersection::Along(guides) = ogeom_intersect::intersect_surfaces(
            &offset(first, sides.first),
            &offset(second, sides.second),
            IntersectOptions::default(),
            T,
        )
        .unwrap() else {
            return Vec::new();
        };
        let options = Marching {
            chord: (radius * 3e-6).max(T.confusion() * 0.1),
            ..Marching::default()
        };
        let mut found = Vec::new();
        for guide in guides.into_iter().map(|c| c.curve) {
            let (lo, hi) = guide.domain();
            let Ok(blend) = march_blend_seeded(
                first,
                second,
                radius,
                &guide,
                sides,
                f64::midpoint(lo, hi),
                options,
                T,
            ) else {
                continue;
            };
            if blend.len() < 8 {
                continue;
            }
            let mut order: Vec<usize> = (0..blend.len()).collect();
            order.sort_by(|&i, &j| blend.along[i].total_cmp(&blend.along[j]));
            order.dedup_by(|j, i| blend.along[*j] <= blend.along[*i] + T.parametric());
            let blend = reordered(&blend, &order);
            found.push((guide, blend));
        }
        found
    }

    /// How far `line` stands from the ball's true contact on the first
    /// support at guide parameter `w`, between stations `i` and `i + 1`:
    /// at the same parameter, or, where that is past `within`, both ways
    /// over the stations either side, the contact from the line and the
    /// line's point from the contact anywhere there.
    #[allow(clippy::too_many_arguments, reason = "one reading, all its data")]
    fn off_contact(
        blend: &MarchedBlend,
        hosts: &Hosts,
        guide: &Curve,
        radius: f64,
        line: &ogeom_geom::BSplineCurve,
        i: usize,
        w: f64,
        within: f64,
    ) -> f64 {
        let [first, second] = &hosts.surfaces;
        let contact = |t: f64| -> Option<Point> {
            let k = if (t - blend.along[i]).abs() < (t - blend.along[i + 1]).abs() {
                i
            } else {
                i + 1
            };
            let near = [
                blend.on_first[k].0,
                blend.on_first[k].1,
                blend.on_second[k].0,
                blend.on_second[k].1,
            ];
            let x = seat_section(first, second, radius, guide, blend.sides, t, near, T).ok()?;
            first.point_at(x[0], x[1], T).ok()
        };
        let truth = contact(w).expect("the ball seats between two stations");
        let point = line.point_at(w, T).unwrap();
        let late = point.distance(truth);
        if late <= within {
            return late;
        }
        let (lo, hi) = (
            blend.along[i.saturating_sub(1)],
            blend.along[(i + 2).min(blend.len() - 1)],
        );
        // Scanned finely, then the best sample refined.
        let nearest = |g: &dyn Fn(f64) -> f64| -> f64 {
            let at = |k: u32| (hi - lo).mul_add(f64::from(k) / 64.0, lo);
            let best = (0..=64)
                .min_by(|&a, &b| g(at(a)).total_cmp(&g(at(b))))
                .unwrap();
            let step = (hi - lo) / 64.0;
            let refined = solve::minimize(
                g,
                (at(best) - step).max(lo),
                (at(best) + step).min(hi),
                solve::Criteria {
                    residual: 0.0,
                    step: 1e-15,
                    max_iterations: 200,
                },
            )
            .unwrap();
            g(at(best)).min(g(refined.value))
        };
        let to_line = nearest(&|t| line.point_at(t, T).unwrap().distance(truth));
        let to_contact = nearest(&|t| contact(t).map_or(f64::INFINITY, |p| p.distance(point)));
        to_line.max(to_contact).min(late)
    }

    /// The farthest the bare fit through the stations and the held line
    /// stand from the ball's true contact, read every eighth of the way
    /// between stations, against the error the held line states.
    fn read(hosts: &Hosts, radius: f64) -> Vec<(f64, f64, f64)> {
        let target = T.confusion() * 10.0;
        let mut out = Vec::new();
        for (guide, blend) in marches(hosts, radius) {
            let Ok(held) = contact_line(&blend, 0, hosts, &guide, radius, target, T) else {
                continue;
            };
            let stated = target.max(held.error + T.confusion());
            let bare =
                ogeom_geom::fit::fit_points_at(&blend.along, &blend.touch_first, 3, target, T)
                    .unwrap()
                    .curve;
            let (mut held_off, mut bare_off) = (0.0_f64, 0.0_f64);
            for i in 0..blend.len() - 1 {
                for k in 1..8 {
                    let w = (blend.along[i + 1] - blend.along[i])
                        .mul_add(f64::from(k) / 8.0, blend.along[i]);
                    let off = off_contact(&blend, hosts, &guide, radius, &held.line, i, w, stated);
                    held_off = held_off.max(off);
                    let off = off_contact(&blend, hosts, &guide, radius, &bare, i, w, target);
                    bare_off = bare_off.max(off);
                }
            }
            out.push((stated, held_off, bare_off));
        }
        out
    }

    /// Where the patch twists and ripples, the ball's touch point on the
    /// floor runs unevenly in the guide's parameter, and a fit through
    /// the march's stations alone strays from the true contact between
    /// them by ten times its target. The line held to the contact between
    /// the stations stands within the error it states, read every eighth
    /// of the way between them both ways: the contact from the line and
    /// the line from the contact.
    #[test]
    fn a_line_of_contact_holds_its_stated_error_between_the_stations() {
        let hosts = floor_and_patch(&|y| 25.0_f64.mul_add(y / 8.0, 55.0), 0.1, 9, 0.4);
        let target = T.confusion() * 10.0;
        let lines = read(&hosts, 2.0);
        assert!(!lines.is_empty(), "a line of contact is read on the floor");
        for (stated, held, _) in &lines {
            assert!(
                held <= stated,
                "the line stands {held:e} off, stating {stated:e}"
            );
        }
        let strayed = lines.iter().map(|l| l.2).fold(0.0, f64::max);
        assert!(
            strayed > target * 5.0,
            "the stations alone leave the contact between them: {strayed:e}"
        );
    }
}
