//! Every branch of an intersection: seed, trace each, and keep the distinct ones.

use super::*;

/// Every branch of the intersection: seed, trace each, and keep the distinct
/// ones.
///
/// A branch crossing many sampling cells produces many seeds, and tracing from
/// any of them gives the same curve. So a seed already lying on something
/// traced is dropped rather than followed again, which is what makes the
/// *number* of branches returned meaningful, and it is the number a boolean
/// will act on.
///
/// Where the surfaces stand within the confusion distance `ε` of each other
/// across a band rather than along a curve, the section is that band to
/// tolerance. Two curves of it in one band, with nowhere between them that
/// the surfaces part by more than `ε`, may come back as one: each sample
/// lies on both surfaces within `ε`, and each point of the curve left out
/// joins one returned by a straight way that does too. A branch is placed
/// across the section only to within `ε` over the sine of the angle the
/// surfaces cross at, so where they cross shallowly its samples may zigzag
/// within that band and its length exceed the curve's.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if the settings are
/// unusable.
pub fn branches(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> OgeomResult<Vec<Traced>> {
    branches_cut_at_poles(a, b, options, 0.0, tol)
}

/// Every branch of the intersection, as [`branches`], each cut at a pole
/// of either chart lying on the other surface within `reach` of it as well
/// as within the chord (see [`cut_at_poles`]).
///
/// A section fitted whole past a pole states a tolerance about its own
/// distance from the pole, since no image follows the chart's half turn
/// there; where it passes the pole within that tolerance
/// ([`passes_a_pole`]), the pole is a point of the section as far as the
/// fit can say, and the pair is traced again cut there.
///
/// # Errors
///
/// As [`branches`].
pub fn branches_cut_at_poles(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    reach: f64,
    tol: Tolerances,
) -> OgeomResult<Vec<Traced>> {
    let (traced, touched) = branches_from_touches(a, b, options, reach, tol)?;
    Ok(traced
        .into_iter()
        .map(|(branch, _)| branch)
        .chain(touched)
        .collect())
}

/// Each open branch's end that stopped short of a pole of either chart the
/// other surface passes through, carried onto the pole.
///
/// Where the section runs into a pole the two surfaces nearly touch at (a
/// transfinite patch's point side lying on a plane it is all but tangent
/// to there), the surfaces stay within the march's chord of each other for
/// a stretch, and the walk (or a seed standing in that stretch) ends short
/// of the pole: the section's end then stands off the vertex the other
/// curves end at. An end within [`joining_reach`] of such a pole is the
/// section running into it, its last state carried onto the pole: on the
/// collapsed chart at the end's own free parameter, on the other surface
/// at the pole's foot found from the end's state. Every surface pair whose
/// section ends there is carried alike, so the curves meet at the pole.
fn join_ends_to_poles(
    branches: &mut [Traced],
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) {
    let reach = joining_reach(options, tol);
    let mut poles: Vec<(usize, usize, f64, Point)> = Vec::new();
    for (surface, other, first) in [(a, b, 0_usize), (b, a, 2)] {
        for (k, bound, point) in chart_poles(surface, tol) {
            if ogeom_geom::project_on_surface(other, point, 16, tol)
                .is_ok_and(|foot| foot.distance <= tol.confusion() * 10.0)
            {
                poles.push((first, k, bound, point));
            }
        }
    }
    if poles.is_empty() {
        return;
    }
    for branch in branches.iter_mut() {
        let n = branch.points.len();
        if n == 0 || branch.stopped == Stopped::Closed {
            continue;
        }
        for at_start in [true, false] {
            let i = if at_start { 0 } else { branch.points.len() - 1 };
            let end = branch.points[i];
            let Some(&(first, k, bound, pole)) = poles
                .iter()
                .min_by(|x, y| x.3.distance(end).total_cmp(&y.3.distance(end)))
            else {
                continue;
            };
            let gap = pole.distance(end);
            if gap <= tol.confusion() || gap > reach {
                continue;
            }
            let mut state = [
                branch.on_a[i].0,
                branch.on_a[i].1,
                branch.on_b[i].0,
                branch.on_b[i].1,
            ];
            state[first + k] = bound;
            let other = if first == 0 { 2 } else { 0 };
            let other_surface = if first == 0 { b } else { a };
            let Ok(foot) = ogeom_geom::project::project_on_surface_from(
                other_surface,
                pole,
                (state[other], state[other + 1]),
                tol,
            ) else {
                continue;
            };
            if foot.distance > tol.confusion() * 10.0 {
                continue;
            }
            state[other] = foot.parameters.0;
            state[other + 1] = foot.parameters.1;
            let (on_a, on_b) = ((state[0], state[1]), (state[2], state[3]));
            if at_start {
                branch.points.insert(0, pole);
                branch.on_a.insert(0, on_a);
                branch.on_b.insert(0, on_b);
            } else {
                branch.points.push(pole);
                branch.on_a.push(on_a);
                branch.on_b.push(on_b);
            }
        }
    }
}

/// Every branch of the intersection, as [`branches`], each with whether it
/// leaves a touch where the section crosses itself along one of the
/// touch's arms: see [`ArmsOf::leaves`]; then the tangential contacts
/// walked through the isolated touches the seeds stood beside.
pub(crate) type Branches = (Vec<(Traced, bool)>, Vec<Traced>);

/// See [`Branches`].
pub(crate) fn branches_from_touches(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    pole_reach: f64,
    tol: Tolerances,
) -> OgeomResult<Branches> {
    let (mut traced, mut touched) = branches_traced(a, b, options, pole_reach, tol)?;
    let mut ends: Vec<Traced> = traced.iter().map(|(t, _)| t.clone()).collect();
    join_ends_to_poles(&mut ends, a, b, options, tol);
    for ((t, _), joined) in traced.iter_mut().zip(ends) {
        *t = joined;
    }
    join_ends_to_poles(&mut touched, a, b, options, tol);
    Ok((traced, touched))
}

/// [`branches_from_touches`] before its ends are carried onto the poles.
fn branches_traced(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    pole_reach: f64,
    tol: Tolerances,
) -> OgeomResult<Branches> {
    let Seeded {
        found,
        beside_tangency,
        lying,
    } = seeded(a, b, options, tol)?;
    // A stretch of a border lying in the other surface is a branch already
    // traced, and a seed on it is not followed again.
    let mut out: Vec<Traced> = lying;
    let mut again: Vec<Traced> = Vec::new();
    let reach = options.chord.max(tol.confusion()) * 8.0;
    // The isolated touches decided, the contacts walked through them, and
    // those whose seeds are traced as any other's.
    let isolated: core::cell::RefCell<Vec<crate::contact::Touch>> =
        core::cell::RefCell::new(Vec::new());
    let touched: core::cell::RefCell<Vec<Traced>> = core::cell::RefCell::new(Vec::new());
    let as_any: core::cell::RefCell<Vec<crate::contact::Touch>> =
        core::cell::RefCell::new(Vec::new());
    let follow = |seed: Contact, out: &mut Vec<Traced>, again: &mut Vec<Traced>| {
        // A seed within an isolated touch's reach (see
        // [`isolated_touch`](crate::contact::isolated_touch)) lies on no
        // crossing: a trace from it keeps within the confusion distance of
        // both surfaces, so within the reach, until it runs out of points,
        // and is read as a tangential fragment whose only use is to start
        // the contact walk from its middle. No other branch comes near it.
        // Whether that walk starts is decided once, at the touch: where the
        // walk names a direction there, the contact is walked from the
        // touch, and where it names none, the reach holds nothing. Main
        // reads it at the fragment's middle instead, which agrees with the
        // touch except where the curvatures' ratio lies within the form's
        // drift over the reach of the walk's bound, and there turns on
        // rounding. A touch whose walk fails is traced as any other.
        let beside = |touches: &[crate::contact::Touch]| {
            touches
                .iter()
                .any(|t| t.at.point.distance(seed.point) <= t.reach)
        };
        if beside(&isolated.borrow()) {
            return;
        }
        if !beside(&as_any.borrow())
            && let Some(touch) = crate::contact::isolated_touch(a, b, seed, options, tol)
        {
            if !touch.walks {
                isolated.borrow_mut().push(touch);
                return;
            }
            match crate::contact::trace_tangential(a, b, touch.at, options, tol) {
                Ok(walked) if walked.points.len() >= 4 => {
                    touched.borrow_mut().push(walked);
                    isolated.borrow_mut().push(touch);
                    return;
                }
                _ => as_any.borrow_mut().push(touch),
            }
        }
        // A seed that will not trace (a tangency) is reported by being
        // absent rather than by an error, since the other branches are still
        // real answers.
        let Ok(branch) = trace(a, b, seed, options, tol) else {
            return;
        };
        if branch.points.len() < 2 || is_fragment(&branch, options) {
            return;
        }
        for piece in cut_at_edge_touches(branch, a, b, tol) {
            // A branch whose middle lies on one already traced is that
            // branch again, reached from a seed its trace stopped short of,
            // or one that turned a corner where the section crosses itself
            // and ran on along another: held to find that crossing from.
            // A branch traced before that lies along this piece all its
            // length is, the other way round, this piece again: its walk
            // stopped at a touch this one ran over and was cut at, or at the
            // closure of a chart this one crossed. So a piece whose middle
            // lies only on such branches, and which runs on past them, is
            // the longer trace of their stretch and takes their place.
            let along: Vec<bool> = out
                .iter()
                .map(|other| other.points.iter().all(|&p| passes_near(&piece, p, reach)))
                .collect();
            let n = piece.points.len();
            let middle = if n > 2 {
                piece.points[n / 2]
            } else {
                piece.points[0].midpoint(piece.points[n - 1])
            };
            let on_middle: Vec<bool> = out
                .iter()
                .zip(&along)
                .filter(|(other, _)| passes_near(other, middle, reach))
                .map(|(_, &along)| along)
                .collect();
            let runs_past = || {
                piece.points.iter().any(|&p| {
                    !out.iter()
                        .zip(&along)
                        .any(|(other, &along)| along && passes_near(other, p, reach))
                })
            };
            if !on_middle.is_empty() && (on_middle.contains(&false) || !runs_past()) {
                again.push(piece);
                continue;
            }
            for (other, along) in core::mem::take(out).into_iter().zip(along) {
                if along {
                    again.push(other);
                } else {
                    out.push(other);
                }
            }
            out.push(piece);
        }
    };
    let mut followed: Vec<Point> = Vec::new();
    let follow_all = |seeds: &[Contact],
                      out: &mut Vec<Traced>,
                      again: &mut Vec<Traced>,
                      followed: &mut Vec<Point>| {
        for &seed in seeds {
            // Already on something followed.
            if out
                .iter()
                .any(|branch| passes_near(branch, seed.point, reach))
            {
                continue;
            }
            follow(seed, out, again);
        }
        // A touch of an edge is where two arcs of the section end, one
        // either side. A seed there lies on the arc already followed to it
        // and is skipped, and the arc on the other side is found only if
        // some seed lies on it: so each touch ending one arc alone is
        // followed from, each way, and the arc it runs back along is that
        // arc again.
        loop {
            let alone: Vec<Contact> = lone_touches(out, a, b, tol)
                .into_iter()
                .filter(|touch| {
                    followed
                        .iter()
                        .all(|p| p.distance(touch.point) > tol.confusion())
                })
                .collect();
            if alone.is_empty() {
                break;
            }
            for touch in alone {
                followed.push(touch.point);
                follow(touch, out, again);
            }
        }
    };
    follow_all(&found, &mut out, &mut again, &mut followed);
    let joined = |out: Vec<Traced>, again: Vec<Traced>| {
        let stitched: Vec<Traced> = stitch_stalled(out, a, b, options, tol)
            .into_iter()
            .map(|branch| first_lap(branch, options.chord.max(tol.confusion())))
            .collect();
        split_at_touches(stitched, again, a, b, options, tol)
            .into_iter()
            .flat_map(|(branch, crossing)| {
                cut_at_poles(branch, a, b, options.chord.max(pole_reach), tol)
                    .into_iter()
                    .map(move |piece| (piece, crossing))
            })
            .collect::<Vec<_>>()
    };
    if beside_tangency.is_empty() {
        return Ok((joined(out, again), touched.take()));
    }
    // A contact held back beside a tangency lies on the section, so some
    // branch must pass it: where none does, it is followed in its own
    // right, and the section joined again with what it finds.
    let first = joined(out.clone(), again.clone());
    let missed: Vec<Contact> = beside_tangency
        .into_iter()
        .filter(|c| {
            !first
                .iter()
                .any(|(branch, _)| passes_near(branch, c.point, reach))
        })
        .collect();
    // Each piece followed is kept in one list or the other.
    let pieces = out.len() + again.len();
    follow_all(&missed, &mut out, &mut again, &mut followed);
    if out.len() + again.len() == pieces {
        return Ok((first, touched.take()));
    }
    Ok((joined(out, again), touched.take()))
}

/// The touches of an edge that end one branch alone, as contacts to follow
/// the section from.
fn lone_touches(
    branches: &[Traced],
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Vec<Contact> {
    let pair = SurfacePair { a, b };
    let ends: Vec<(Point, [f64; 4])> = branches
        .iter()
        .filter(|branch| !branch.closed() && !branch.points.is_empty())
        .flat_map(|branch| {
            let last = branch.points.len() - 1;
            [0, last].map(|i| {
                let (on_a, on_b) = (branch.on_a[i], branch.on_b[i]);
                (branch.points[i], [on_a.0, on_a.1, on_b.0, on_b.1])
            })
        })
        .collect();
    let mut alone: Vec<Contact> = Vec::new();
    for &(point, state) in &ends {
        let touch = [(a, 0_usize), (b, 2)]
            .into_iter()
            .find_map(|(surface, first)| {
                let ((ua, ub), (va, vb)) = surface.domain();
                [
                    (first, (ua, ub), surface.is_periodic_u()),
                    (first + 1, (va, vb), surface.is_periodic_v()),
                ]
                .into_iter()
                .filter(|&(_, _, periodic)| !periodic)
                .flat_map(|(k, (lo, hi), _)| [(k, lo), (k, hi)])
                .filter(|&(k, bound)| (state[k] - bound).abs() <= tol.parametric())
                .find_map(|(k, bound)| touch_on_edge(&pair, &state, k, bound, tol))
            });
        let Some(at) = touch else {
            continue;
        };
        if a.point_at(at[0], at[1], tol)
            .is_ok_and(|p| p.distance(point) > tol.confusion())
        {
            continue;
        }
        let meeting = ends
            .iter()
            .filter(|(other, _)| other.distance(point) <= tol.confusion())
            .count();
        if meeting == 1 {
            alone.push(Contact {
                on_a: (at[0], at[1]),
                on_b: (at[2], at[3]),
                point,
            });
        }
    }
    alone
}

/// A branch cut where it runs over a touch of either surface's edge: the
/// pieces between the touches, each ending on the touch.
///
/// A walk ends at a touch of an edge where its step lands on the edge
/// ([`end_at_touch`]). A step can as well land beyond the touch, back inside
/// the domain, and the walk runs on over it. A section is then split at a
/// touch or not by where the steps fell, and the same stretch comes back
/// from two seeds as two different curves: one ending at the touch, one
/// running over it. So a branch is cut at every touch it runs over.
///
/// A touch run over is beside a sample where one parameter turns back
/// toward the edge it bounds. It is solved there by [`touch_on_edge`], and
/// kept where it lies within the steps either side of the sample.
fn cut_at_edge_touches(
    branch: Traced,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Vec<Traced> {
    let pair = SurfacePair { a, b };
    let closed = branch.closed();
    let n = branch.points.len();
    if n < 3 {
        return vec![branch];
    }
    let state = |i: usize| {
        [
            branch.on_a[i].0,
            branch.on_a[i].1,
            branch.on_b[i].0,
            branch.on_b[i].1,
        ]
    };
    // A loop's last sample is its first.
    let count = if closed { n - 1 } else { n };
    let period = |surface: &SurfaceGeometry, along_u: bool| {
        let ((ua, ub), (va, vb)) = surface.domain();
        if along_u {
            surface.is_periodic_u().then_some(ub - ua)
        } else {
            surface.is_periodic_v().then_some(vb - va)
        }
    };
    let periods = [
        period(a, true),
        period(a, false),
        period(b, true),
        period(b, false),
    ];
    let mut touches: Vec<(usize, [f64; 4], Point)> = Vec::new();
    for i in 0..count {
        let (before, after) = if closed {
            ((i + count - 1) % count, (i + 1) % count)
        } else if i == 0 || i + 1 == n {
            continue;
        } else {
            (i - 1, i + 1)
        };
        let here = state(i);
        let (back, ahead) = (
            branch.points[i].distance(branch.points[before]),
            branch.points[i].distance(branch.points[after]),
        );
        if back <= 0.0 || ahead <= 0.0 {
            continue;
        }
        let stride = back.max(ahead);
        // A periodic parameter is read on the turn of the middle sample's.
        let unwrap = |x: [f64; 4]| -> [f64; 4] {
            core::array::from_fn(|j| match periods[j] {
                Some(period) => x[j] - ((x[j] - here[j]) / period).round() * period,
                None => x[j],
            })
        };
        let around = [unwrap(state(before)), here, unwrap(state(after))];
        'bounds: for (surface, first) in [(a, 0_usize), (b, 2)] {
            let ((ua, ub), (va, vb)) = surface.domain();
            for (k, (lo, hi), periodic) in [
                (first, (ua, ub), surface.is_periodic_u()),
                (first + 1, (va, vb), surface.is_periodic_v()),
            ] {
                if periodic || !(hi - lo).is_finite() {
                    continue;
                }
                for bound in [lo, hi] {
                    let off = |x: &[f64; 4]| (x[k] - bound).abs();
                    if off(&here) >= off(&around[0]) || off(&here) > off(&around[2]) {
                        continue;
                    }
                    // A touch within the steps either side lies on an edge
                    // passing within them.
                    let mut on_edge = here;
                    on_edge[k] = bound;
                    let edge = surface.point_at(on_edge[first], on_edge[first + 1], tol);
                    if edge.is_ok_and(|p| p.distance(branch.points[i]) > stride) {
                        continue;
                    }
                    let start = turning(around, (back, ahead), k);
                    let Some(at) = touch_on_edge(&pair, &start, k, bound, tol) else {
                        continue;
                    };
                    let Ok(point) = a.point_at(at[0], at[1], tol) else {
                        continue;
                    };
                    if point.distance(branch.points[i]) <= stride {
                        touches.push((i, at, point));
                        break 'bounds;
                    }
                }
            }
        }
    }
    if touches.is_empty() {
        return vec![branch];
    }
    // The samples with each touch in its place, and where the touches are.
    let mut samples: Vec<(Point, [f64; 4])> = Vec::with_capacity(n + touches.len());
    let mut cuts: Vec<usize> = Vec::new();
    let mut next = touches.iter().peekable();
    for i in 0..count {
        let Some(&(_, at, point)) = next.next_if(|t| t.0 == i) else {
            samples.push((branch.points[i], state(i)));
            continue;
        };
        let (before, after) = (
            branch.points[if i == 0 { count - 1 } else { i - 1 }],
            branch.points[(i + 1) % n],
        );
        if point.distance(branch.points[i]) <= tol.confusion() {
            cuts.push(samples.len());
            samples.push((point, at));
        } else if (point - branch.points[i]).dot(after - before) > 0.0 {
            samples.push((branch.points[i], state(i)));
            cuts.push(samples.len());
            samples.push((point, at));
        } else {
            cuts.push(samples.len());
            samples.push((point, at));
            samples.push((branch.points[i], state(i)));
        }
    }
    let piece = |run: &[(Point, [f64; 4])], stopped: Stopped| Traced {
        points: run.iter().map(|s| s.0).collect(),
        on_a: run.iter().map(|s| (s.1[0], s.1[1])).collect(),
        on_b: run.iter().map(|s| (s.1[2], s.1[3])).collect(),
        stopped,
    };
    let mut pieces = Vec::new();
    if closed {
        // Round the loop from its first touch back to it.
        let first = cuts[0];
        let mut round: Vec<(Point, [f64; 4])> = samples[first..].to_vec();
        round.extend_from_slice(&samples[..=first]);
        let mut from = 0;
        for cut in cuts
            .iter()
            .skip(1)
            .map(|c| c - first)
            .chain([round.len() - 1])
        {
            if cut > from {
                pieces.push(piece(&round[from..=cut], Stopped::LeftTheDomain));
            }
            from = cut;
        }
    } else {
        let mut from = 0;
        for cut in cuts.iter().copied().chain([samples.len() - 1]) {
            if cut > from {
                pieces.push(piece(&samples[from..=cut], Stopped::LeftTheDomain));
            }
            from = cut;
        }
        // The ends the trace stopped at are its own: a branch truncated or
        // stalled at either end is reported so in its first and last piece.
        if let Some(first) = pieces.first_mut() {
            first.stopped = branch.stopped;
        }
        if let Some(last) = pieces.last_mut() {
            last.stopped = branch.stopped;
        }
    }
    pieces
}
