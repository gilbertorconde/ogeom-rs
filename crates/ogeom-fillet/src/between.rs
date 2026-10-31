//! Whether a face of the solid stands inside the ball rolling between two
//! faces it blends: a strip between them (a chamfer, a smaller round, a
//! wedge fused into the corner over part of its run) that a round of the
//! two faces alone would leave standing through the round.
//!
//! The ball's centre is known exactly at any station of its run, so the
//! question is asked over the whole run, not at samples of it. A face's
//! distance from the centre falls between two stations by at most the
//! centre's travel from the nearer of them, half the path's length there.
//! A convex arc is at most its chord and twice its height long, so that
//! is half the chord and the path's sagitta, read at the exactly seated
//! middle station. An interval is halved until that travel fits under the
//! face's clearance at its ends and what is read as touching, or a point
//! of the face is found inside a ball. A face's distance is the kernel's
//! point-to-face distance. A refusal rests on the point it names.

use ogeom_algo::PointDistance;
use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_geom::Curve;
use ogeom_math::Point;
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

/// How many times an interval of the run is halved, at most, to show that
/// the ball stays clear of a face between two stations.
const HALVINGS: u32 = 24;

/// The balls a round between two faces rolls along its run.
pub(crate) struct Rolling<'a> {
    pub(crate) radius: f64,
    /// Stations along the run, in order, from its start to its end.
    pub(crate) stations: Vec<f64>,
    /// The ball's centre at a station, seated exactly.
    pub(crate) centre: &'a dyn Fn(f64) -> OgeomResult<Point>,
    /// The crease and the lines of contact, each a curve over the run. A
    /// face one of them meets is where the run ends (a crease running off
    /// a block's side, a face a line of contact leaves by) or a face along
    /// the corner itself (a wedge's end standing over the crease), and is
    /// not asked: the round runs off it there, and a ball reaching through
    /// it says nothing of a face between the two.
    pub(crate) lines: Vec<(Curve, (f64, f64))>,
    /// Where the run ends: the contacts and the crease in its end
    /// sections. A face holding one ends the run and is not asked.
    pub(crate) ends: Vec<Point>,
}

/// A point of a face of `solid` other than the two `faces`, and other than
/// the faces the run ends on or the corner runs along, lying more than the
/// confusion inside the ball at some station of the run; `None` where every
/// such face is shown clear of every ball, or within what is read as the
/// ball touching it.
///
/// The ball touching both faces at a station is the one the round's arc
/// there lies on, so such a point is a point of the solid's boundary inside
/// it: no ball of the radius rolls there touching the two faces alone, and
/// a round of the two would leave that face standing through it.
///
/// # Errors
///
/// [`OgeomError::NotDone`](ogeom_core::OgeomError::NotDone) where a face
/// comes within the centre's travel of a ball and halving the interval
/// as far as [`HALVINGS`] neither clears it nor finds it inside, or where
/// the ball cannot be seated between two stations.
pub(crate) fn face_inside_balls(
    model: &mut Model,
    solid: &Shape,
    faces: [&Shape; 2],
    rolling: &Rolling<'_>,
    tol: Tolerances,
) -> OgeomResult<Option<Point>> {
    let conf = tol.confusion();
    let radius = rolling.radius;
    // A face reaching into the ball by less than the faces are read at for
    // whether the ball touches them is read as touching it, as a contact
    // on a face is. One reaching in by more than the confusion at a station
    // is inside it.
    let touching = crate::pair_marched::reading(radius, tol).chord;
    let stations = &rolling.stations;
    if stations.len() < 2 {
        ogeom_bail!(Invariant, "a run has a start and an end");
    }
    let centres = stations
        .iter()
        .map(|&s| (rolling.centre)(s))
        .collect::<OgeomResult<Vec<Point>>>()?;
    let mut lines = Vec::with_capacity(rolling.lines.len());
    for (curve, range) in &rolling.lines {
        lines.push(ogeom_algo::make_edge(model, curve.clone(), *range, tol)?.shape);
    }
    // Each interval's middle centre, and the centre's travel over it.
    let mut middles = Vec::with_capacity(stations.len() - 1);
    for i in 0..stations.len() - 1 {
        let pm = (rolling.centre)(f64::midpoint(stations[i], stations[i + 1]))?;
        middles.push(pm);
    }
    let travel = |pa: Point, pm: Point, pb: Point| {
        pa.distance(pb)
            .mul_add(0.5, pm.distance(Point::midpoint(pa, pb)))
    };
    // No interval's travel reaches past this, so a face whose bound stands
    // farther off than a radius and this is clear of every ball.
    let reach_past = (0..middles.len())
        .map(|i| travel(centres[i], middles[i], centres[i + 1]))
        .fold(0.0_f64, f64::max)
        + conf;

    for face in explore_unique(model, solid, ShapeType::Face)? {
        if faces.iter().any(|f| f.is_same(&face)) {
            continue;
        }
        let bounds = ogeom_algo::shape_bounds(model, &face, tol)?;
        if centres
            .iter()
            .all(|c| bounds.distance_to(*c) - radius > reach_past)
        {
            continue;
        }
        let measure = PointDistance::of(model, &face, tol)?;
        // A point the face holds, within its own tolerance.
        let slack = face_slack(model, &face, tol)?;
        let mut ends_here = false;
        for &p in &rolling.ends {
            if bounds.distance_to(p) <= slack && measure.nearest(p, tol)?.0 <= slack {
                ends_here = true;
                break;
            }
        }
        if !ends_here {
            for line in &lines {
                let gap = ogeom_algo::distance_between_shapes(
                    model,
                    line,
                    &face,
                    ogeom_intersect::ExtremaOptions::default(),
                    tol,
                )?;
                if gap.distance <= slack {
                    ends_here = true;
                    break;
                }
            }
        }
        if ends_here {
            continue;
        }

        // The clearance from the ball: the face's bound where that clears
        // any interval, else the face's own distance and its nearest point.
        let clearance = |c: Point| -> OgeomResult<(f64, Option<Point>)> {
            let far = bounds.distance_to(c) - radius;
            if far > reach_past {
                return Ok((far, None));
            }
            let (d, q) = measure.nearest(c, tol)?;
            Ok((d - radius, Some(q)))
        };
        let inside = |(gap, at): (f64, Option<Point>)| (gap < -conf).then_some(at).flatten();
        let at_stations: Vec<(f64, Option<Point>)> = centres
            .iter()
            .map(|c| clearance(*c))
            .collect::<OgeomResult<_>>()?;
        if let Some(at) = at_stations.iter().find_map(|e| inside(*e)) {
            return Ok(Some(at));
        }
        let mut intervals: Vec<Interval> = (0..middles.len())
            .map(|i| Interval {
                at: (stations[i], stations[i + 1]),
                ends: (centres[i], centres[i + 1]),
                middle: Some(middles[i]),
                gaps: (at_stations[i].0, at_stations[i + 1].0),
                depth: 0,
            })
            .collect();
        while let Some(Interval {
            at: (a, b),
            ends: (pa, pb),
            middle,
            gaps: (ga, gb),
            depth,
        }) = intervals.pop()
        {
            let m = f64::midpoint(a, b);
            let pm = match middle {
                Some(pm) => pm,
                None => (rolling.centre)(m)?,
            };
            if ga.min(gb) + touching >= travel(pa, pm, pb) {
                continue;
            }
            let gm = clearance(pm)?;
            if let Some(at) = inside(gm) {
                return Ok(Some(at));
            }
            if depth >= HALVINGS || m <= a || m >= b {
                ogeom_bail!(
                    NotDone,
                    "whether another face of the solid comes inside the ball of radius \
                     {radius} rolling between the two faces cannot be told: it stands within \
                     {} of the ball between two stations {} apart",
                    ga.min(gb).max(0.0),
                    pa.distance(pb)
                );
            }
            for (at, ends, gaps) in [
                ((a, m), (pa, pm), (ga, gm.0)),
                ((m, b), (pm, pb), (gm.0, gb)),
            ] {
                intervals.push(Interval {
                    at,
                    ends,
                    middle: None,
                    gaps,
                    depth: depth + 1,
                });
            }
        }
    }
    Ok(None)
}

/// A stretch of the run still to be shown clear of a face.
struct Interval {
    /// Its stations.
    at: (f64, f64),
    /// The ball's centre at each.
    ends: (Point, Point),
    /// The centre at its middle station, where already seated.
    middle: Option<Point>,
    /// The face's clearance from the ball at each end, or a lower bound.
    gaps: (f64, f64),
    /// How many times the run was halved to reach it.
    depth: u32,
}

/// How far a point may stand from a face and still be held by it: twice
/// the largest tolerance of its edges and vertices, and a hundred times the
/// confusion.
fn face_slack(model: &Model, face: &Shape, tol: Tolerances) -> OgeomResult<f64> {
    let mut most = 0.0_f64;
    for edge in explore_unique(model, face, ShapeType::Edge)? {
        if let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) {
            most = most.max(data.tolerance.get());
        }
    }
    for vertex in explore_unique(model, face, ShapeType::Vertex)? {
        if let Some(data) = model.node(&vertex).and_then(|n| n.data().as_vertex()) {
            most = most.max(data.tolerance.get());
        }
    }
    Ok(most.mul_add(2.0, tol.confusion() * 100.0))
}
