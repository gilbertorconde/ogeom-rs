//! Branches cut where the surfaces touch and the section crosses itself.

use super::*;

/// Branches cut where the surfaces touch and the section crosses itself.
///
/// Where two surfaces touch at a point and bend apart differently across
/// it (a drill lying inside a ring's outer wall and touching it), the
/// section is a figure eight: two loops meeting at the touch, four arms
/// leaving it. The walk has no single direction there. It runs straight
/// through on one arm, so one branch goes round both loops, or it turns
/// onto the next arm and closes one loop, or it stalls a short way off
/// the touch. As one branch the figure eight is no ring a face can be
/// split by.
///
/// So each touch is found exactly from the low-angle stretches of the
/// branches, those traced `again` (a walk that turned a corner and ran on
/// over a branch already held) among them, and a branch crossing itself
/// there, turning a corner on it, ending on it or stopping short of it is
/// cut there: the samples within reach of it go and the arcs either side
/// end on it, a stalled end heading for it carried on to it. A branch
/// running once straight through each touch it passes is a smooth curve
/// another crosses there, and stays whole; where the touch's two arms are
/// known, straight through means in along one arm and on along it. An arc
/// starting and ending on one
/// touch is a whole loop, cut again so each piece is open with two
/// distinct ends: off its middle sample, which on a loop symmetric about
/// the touch lands on whatever else lies on the plane of symmetry, a seam
/// crossing another seam among them. An arm of a touch that nothing runs
/// along is then followed out from the touch, and a branch or arc lying
/// wholly on a longer one is dropped, or cut back where its ends do.
/// Branches nowhere near a touch, and those running along a tangency, are
/// returned as they came.
///
/// An arc shorter than ten chords, or of fewer than four samples, is
/// dropped as debris, unless it joins two touches or runs from a touch
/// out along its arm to the border of either chart: a branch running off
/// the edge of a surface within a few chords of a touch keeps the stub
/// between the two, which is the section there however short.
///
/// Each arc comes with whether it leaves a touch along one of its arms:
/// see [`ArmsOf::leaves`].
pub(super) fn split_at_touches(
    found: Vec<Traced>,
    again: Vec<Traced>,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> Vec<(Traced, bool)> {
    let near = options.chord.max(tol.confusion());
    let mut touches: Vec<Touch> = Vec::new();
    let solved_ends: Vec<[Option<usize>; 2]> = found
        .iter()
        .map(|branch| find_touches(branch, &mut touches, (a, b), options, tol))
        .collect();
    for branch in &again {
        find_touches(branch, &mut touches, (a, b), options, tol);
    }
    let uncertified = |arcs: Vec<Traced>| arcs.into_iter().map(|arc| (arc, false)).collect();
    if touches.is_empty() {
        return uncertified(found);
    }
    settle_touches(&mut touches, (a, b), options, tol);
    // The march strides up to its longest step, and two touches within one
    // stride of each other may stand in one step, the section's direction
    // turning about twice in it and reading the same at both ends: a touch
    // between them goes unfound, and cutting at those found loses the
    // stretch the step covered.
    let stride = found
        .iter()
        .chain(&again)
        .flat_map(|branch| branch.points.windows(2).map(|w| w[0].distance(w[1])))
        .fold(0.0, f64::max);
    let crowded = touches.iter().enumerate().any(|(i, t)| {
        t.arms.is_some()
            && touches
                .iter()
                .enumerate()
                .any(|(j, u)| j != i && t.point.distance(u.point) <= stride)
    });
    // Every crossing of the pair is then found by cutting both charts into
    // patches, whatever the march's stride: a crossing found alone may
    // have another within a stride of it just the same. Where crossings
    // crowd and that search does not complete, the branches are cut as
    // they were traced.
    let crossing = touches.iter().any(|t| t.arms.is_some());
    let complete = crossing
        && subdivided_touches(a, b, options, tol).is_some_and(|crossings| {
            for touch in crossings {
                if !touches.iter().any(|t| same_touch(t, &touch, near)) {
                    touches.push(touch);
                }
            }
            true
        });
    if crowded && !complete {
        return split_at_touches_as_traced(found, a, b, options, tol);
    }
    let certified =
        |arcs: Vec<Traced>, touches: &[Touch]| certified(arcs, touches, (a, b), near, tol);
    if complete {
        settle_touches(&mut touches, (a, b), options, tol);
    }

    // Branches no touch cuts stay as they are, the others are cut into
    // arcs; each notes the arms of the touches it runs along.
    let mut kept: Vec<Traced> = Vec::with_capacity(found.len());
    let mut arcs: Vec<Traced> = Vec::new();
    let mut covered: Vec<(usize, usize)> = Vec::new();
    for (branch, solved) in found.into_iter().zip(solved_ends) {
        let cut = cut_at_touches(branch, solved, &touches, (a, b), options, tol);
        covered.extend(cut.covered);
        if cut.cut {
            arcs.extend(cut.pieces);
        } else {
            kept.extend(cut.pieces);
        }
    }
    let (joins, walked_out) = join_arms(&touches, &mut covered, &[], (a, b), options, tol);
    arcs.extend(joins);
    if walked_out.is_empty() && arcs.is_empty() {
        return certified(kept, &touches);
    }
    let solved_out: Vec<[Option<usize>; 2]> = walked_out
        .iter()
        .map(|arm| find_touches(arm, &mut touches, (a, b), options, tol))
        .collect();
    settle_touches(&mut touches, (a, b), options, tol);
    let walked_any = !walked_out.is_empty();
    for (arm, solved) in walked_out.into_iter().zip(solved_out) {
        arcs.extend(cut_at_touches(arm, solved, &touches, (a, b), options, tol).pieces);
    }

    // An arc lying wholly on another, every sample clear of the touches
    // within the chord of it, is that stretch twice: longest first, each
    // is kept only if some sample of it lies off everything kept so far,
    // and cut back where its ends run along what is kept. A branch no
    // touch cut is weighed too where arms were followed out, as one may
    // run along it. Two arcs between the same two touches, leaving each
    // along the same arm, are one stretch: near the touches neither may
    // have a sample clear of them to weigh.
    let mut weighed = arcs;
    if walked_any {
        weighed.append(&mut kept);
    }
    weighed.sort_by(|x, y| arc_length(y).total_cmp(&arc_length(x)));
    let weigh = |arc: Traced, kept: &mut Vec<Traced>, touches: &[Touch]| {
        for arc in straightened(arc, touches, (a, b), near, tol) {
            let across = across_section(&arc, touches, near, (a, b), tol);
            let arms = ArmsOf::new(&arc, touches, (a, b), near, tol);
            if kept.iter().any(|other| {
                lies_on(&arc, other, near, (&across, &arms), touches)
                    || same_stretch(&arc, other, near, touches)
                    || stub_lies_on(&arc, other, near, touches)
            }) {
                continue;
            }
            if let Some(arc) = trimmed(arc, kept, near, (&across, &arms), touches) {
                kept.push(arc);
            }
        }
    };
    for arc in weighed {
        weigh(arc, &mut kept, &touches);
    }
    // With every touch found, an arm the arcs still do not run along is
    // section left out: joined or walked out once more from what is kept.
    if complete {
        let mut covered = arms_covered(&kept, &touches);
        let ends = loose_ends(&kept, &touches, near);
        let (joins, walks) = join_arms(&touches, &mut covered, &ends, (a, b), options, tol);
        let mut more = joins;
        for arm in walks {
            more.extend(cut_at_touches(arm, [None, None], &touches, (a, b), options, tol).pieces);
        }
        more.sort_by(|x, y| arc_length(y).total_cmp(&arc_length(x)));
        for arc in more {
            weigh(arc, &mut kept, &touches);
        }
    }
    // An arc is cut back only against those kept before it, and an end on
    // a touch is never cut back: where a loose end of one runs along an arc
    // kept after it that ends on a touch, the stretch between is traced
    // twice. Each is cut back once more against all the others.
    let mut k = 0;
    while k < kept.len() {
        let arc = kept.remove(k);
        let across = across_section(&arc, &touches, near, (a, b), tol);
        let arms = ArmsOf::new(&arc, &touches, (a, b), near, tol);
        if let Some(arc) = trimmed(arc, &kept, near, (&across, &arms), &touches) {
            kept.insert(k, arc);
            k += 1;
        }
    }
    certified(kept, &touches)
}

/// Arcs cut at the touches, each straightened ([`straightened`]) and with
/// whether it leaves a touch along one of its arms ([`ArmsOf::leaves`]).
fn certified(
    arcs: Vec<Traced>,
    touches: &[Touch],
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    near: f64,
    tol: Tolerances,
) -> Vec<(Traced, bool)> {
    arcs.into_iter()
        .flat_map(|arc| straightened(arc, touches, (a, b), near, tol))
        .map(|arc| {
            let leaves = ArmsOf::new(&arc, touches, (a, b), near, tol).leaves(&arc);
            (arc, leaves)
        })
        .collect()
}

/// The arms of the touches `covered` does not list, joined to a facing
/// touch along the section, or else walked out from the touch: the joins,
/// each noted as covering both arms it runs along, and the walks.
fn join_arms(
    touches: &[Touch],
    covered: &mut Vec<(usize, usize)>,
    ends: &[(Contact, Vector)],
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    options: Marching,
    tol: Tolerances,
) -> (Vec<Traced>, Vec<Traced>) {
    let near = options.chord.max(tol.confusion());
    let carry = near * 60.0 * 40.0;
    let mut arcs: Vec<Traced> = Vec::new();
    // Two touches facing each other along arms nothing runs along are
    // joined along the section between them: a stretch between two
    // crossings may be too shallow for a walk to set out on, and is
    // followed by samples placed where the chord calls for them, each
    // settling on the section within half its chord.
    let touch_contact = |t: &Touch| Contact {
        on_a: t.on_a,
        on_b: t.on_b,
        point: t.point,
    };
    for t in 0..touches.len() {
        if touches[t].arms.is_none() {
            continue;
        }
        for h in 0..4 {
            if covered.contains(&(t, h)) {
                continue;
            }
            let facing = (0..touches.len())
                .filter(|&u| u != t && touches[u].arms.is_some())
                .filter(|&u| half_arm(&touches[t], touches[u].point) == Some(h))
                .filter_map(|u| {
                    let back = half_arm(&touches[u], touches[t].point)?;
                    (!covered.contains(&(u, back))).then_some((u, back))
                })
                .min_by(|&(x, _), &(y, _)| {
                    let d = |u: usize| touches[u].point.distance(touches[t].point);
                    d(x).total_cmp(&d(y))
                });
            let from = touches[t].point;
            let ahead = |q: Point| {
                touches.iter().enumerate().any(|(v, w)| {
                    v != t
                        && half_arm(&touches[t], w.point) == Some(h)
                        && w.point.distance(from) < q.distance(from)
                })
            };
            // An arc that stopped short of the touch, heading at it with no
            // other touch on the way, is joined to it the same way.
            let loose = ends
                .iter()
                .filter(|(e, heading)| {
                    half_arm(&touches[t], e.point) == Some(h)
                        && heading.dot(from - e.point) > 0.0
                        && !ahead(e.point)
                })
                .min_by(|x, y| {
                    x.0.point
                        .distance(from)
                        .total_cmp(&y.0.point.distance(from))
                });
            let target = match (facing, loose) {
                (Some((u, _)), Some((e, _)))
                    if e.point.distance(from) < touches[u].point.distance(from) =>
                {
                    (*e, None)
                }
                (Some((u, back)), _) => (touch_contact(&touches[u]), Some((u, back))),
                (None, Some((e, _))) => (*e, None),
                (None, None) => continue,
            };
            let (p, q) = (touch_contact(&touches[t]), target.0);
            let q = Contact {
                on_a: beside(a, q.on_a, p.on_a),
                on_b: beside(b, q.on_b, p.on_b),
                point: q.point,
            };
            let mut budget = options.max_points;
            let Some(between) = fill((a, b), p, q, near, &mut budget, tol) else {
                continue;
            };
            let leaves = between
                .first()
                .is_none_or(|c| half_arm(&touches[t], c.point) == Some(h))
                && target.1.is_none_or(|(u, back)| {
                    between
                        .last()
                        .is_none_or(|c| half_arm(&touches[u], c.point) == Some(back))
                });
            if !leaves {
                continue;
            }
            let mut stretch = Traced {
                points: vec![p.point],
                on_a: vec![p.on_a],
                on_b: vec![p.on_b],
                stopped: Stopped::Stalled,
            };
            for c in between.iter().chain([&q]) {
                stretch.points.push(c.point);
                stretch.on_a.push(c.on_a);
                stretch.on_b.push(c.on_b);
            }
            covered.push((t, h));
            covered.extend(target.1);
            arcs.push(stretch);
        }
    }
    // An arm of a touch nothing runs along is followed out from the touch:
    // a walk turning a corner there, or stalling, leaves the arm it did not
    // take untraced, and every seed on that arm may have walked onto
    // another one first.
    let mut walked_out: Vec<Traced> = Vec::new();
    for (t, touch) in touches.iter().enumerate() {
        let Some(arms) = touch.arms else { continue };
        for (h, half) in halves(arms).into_iter().enumerate() {
            if covered.contains(&(t, h)) {
                continue;
            }
            // The walk starts where the section has a direction it can
            // follow: out from the touch's radius, twice as far each time,
            // as far as a stall is carried onto a touch. An arm that leaves
            // the chart inside the radius is started inside it, half as far
            // each time, no nearer than the touch's snap, inside which the
            // samples of a walk stand anywhere in the band the confusion
            // distance leaves round the arms. An arm leaving the chart
            // short of every start is solved out to where it does.
            let floor = touch.snap.max(near);
            let outwards = core::iter::successors(Some(touch.radius), |r| {
                Some(r * 2.0).filter(|&r| r <= carry)
            });
            let inwards = core::iter::successors(Some(touch.radius * 0.5), |r| Some(r * 0.5))
                .take_while(|&r| r >= floor);
            // The nearest start tried past the snap.
            let shortest = outwards
                .clone()
                .chain(inwards.clone())
                .filter(|&r| r >= floor)
                .fold(f64::INFINITY, f64::min);
            let leaves = |out: f64| {
                matches!(
                    arm_point((a, b), touch, half, (out, arms), near, tol),
                    Landing::Leaves(..)
                )
            };
            let inside_snap = leaves(floor);
            let start = if inside_snap {
                None
            } else {
                outwards.chain(inwards).find_map(|out| {
                    let start = along_arm(a, b, touch, half, (out, arms), options.chord, tol)?;
                    Some((start, tangent_at(a, b, start, tol)?))
                })
            };
            let Some((start, tangent)) = start else {
                let beyond = if inside_snap { floor } else { shortest };
                if (inside_snap || (beyond.is_finite() && leaves(beyond)))
                    && let Some(arm) = arm_to_border(
                        (a, b),
                        touch,
                        half,
                        (beyond, arms),
                        (near, options.max_points),
                        tol,
                    )
                {
                    covered.push((t, h));
                    arcs.push(arm);
                }
                continue;
            };
            let sense = if tangent.dot(half) >= 0.0 { 1.0 } else { -1.0 };
            let Ok(walked) = walk(a, b, start, sense, options, tol) else {
                continue;
            };
            let mut arm = Traced {
                points: vec![touch.point],
                on_a: vec![beside(a, touch.on_a, start.on_a)],
                on_b: vec![beside(b, touch.on_b, start.on_b)],
                stopped: walked.stopped,
            };
            arm.points.extend(walked.points);
            arm.on_a.extend(walked.on_a);
            arm.on_b.extend(walked.on_b);
            walked_out.push(arm);
        }
    }
    (arcs, walked_out)
}

/// The arms of the touches the arcs run along: each run of samples within
/// a touch's radius, or each step passing within it, and the samples
/// either side clear of it.
fn arms_covered(arcs: &[Traced], touches: &[Touch]) -> Vec<(usize, usize)> {
    let mut covered = Vec::new();
    for arc in arcs {
        let m = arc.points.len();
        for (t, touch) in touches.iter().enumerate() {
            if touch.arms.is_none() {
                continue;
            }
            let inside = |k: usize| arc.points[k].distance(touch.point) <= touch.radius;
            for k in 0..m {
                let passes = inside(k)
                    || (k + 1 < m
                        && distance_to_segment(touch.point, arc.points[k], arc.points[k + 1])
                            <= touch.radius);
                if !passes {
                    continue;
                }
                let before = (0..=k).rev().find(|&j| !inside(j));
                let after = (k + 1..m).find(|&j| !inside(j));
                for j in [before, after].into_iter().flatten() {
                    if let Some(h) = half_arm(touch, arc.points[j]) {
                        covered.push((t, h));
                    }
                }
            }
        }
    }
    covered
}

/// The ends of arcs that stop on no touch and meet no other end, each with
/// the way the arc runs out there.
fn loose_ends(arcs: &[Traced], touches: &[Touch], near: f64) -> Vec<(Contact, Vector)> {
    let mut ends = Vec::new();
    let all: Vec<(usize, Point)> = arcs
        .iter()
        .enumerate()
        .filter(|(_, arc)| !arc.points.is_empty())
        .flat_map(|(i, arc)| {
            [
                (2 * i, arc.points[0]),
                (2 * i + 1, arc.points[arc.points.len() - 1]),
            ]
        })
        .collect();
    for (i, arc) in arcs.iter().enumerate() {
        let m = arc.points.len();
        if m < 2 {
            continue;
        }
        for (end, (k, before)) in [(0, 1), (m - 1, m - 2)].into_iter().enumerate() {
            let p = arc.points[k];
            if touches.iter().any(|t| p.distance(t.point) <= t.radius)
                || all
                    .iter()
                    .any(|&(j, q)| j != 2 * i + end && q.distance(p) <= near)
            {
                continue;
            }
            let contact = Contact {
                on_a: arc.on_a[k],
                on_b: arc.on_b[k],
                point: p,
            };
            ends.push((contact, p - arc.points[before]));
        }
    }
    ends
}

/// Whether a branch runs straight through a touch between the samples
/// `before` and `after` it, taken from the samples either side clear of the
/// touch's radius: a stretch between two touches starts on one. Where the
/// touch's two arms are known, the branch runs straight through when it
/// comes in along one arm and goes on along it, the same way; it turns a
/// corner when it goes on along the other. Which arm a sample is on is
/// read from the section's own direction there, the cross of the normals,
/// not the chord to the touch, which a long stride tilts.
pub(super) fn straight_through(
    branch: &Traced,
    (before, after): (&[usize], &[usize]),
    touch: &Touch,
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    tol: Tolerances,
) -> bool {
    let unit = |v: Vector| {
        let m = v.magnitude();
        (m > tol.confusion()).then(|| v / m)
    };
    let clear = |&&k: &&usize| branch.points[k].distance(touch.point) > touch.radius;
    let (Some(&i), Some(&o)) = (before.iter().rev().find(clear), after.iter().find(clear)) else {
        return false;
    };
    let into = unit(touch.point - branch.points[i]);
    let onward = unit(branch.points[o] - touch.point);
    let (Some(into), Some(onward)) = (into, onward) else {
        return false;
    };
    match touch.arms {
        Some(arms) => {
            let heading = |k: usize, chord: Vector| {
                let at = Contact {
                    on_a: branch.on_a[k],
                    on_b: branch.on_b[k],
                    point: branch.points[k],
                };
                tangent_at(a, b, at, tol).unwrap_or(chord)
            };
            let line = |w: Vector| w.dot(arms[1]).abs() > w.dot(arms[0]).abs();
            line(heading(i, into)) == line(heading(o, onward)) && into.dot(onward) > 0.0
        }
        None => into.dot(onward) > 0.5,
    }
}

/// The samples of an open branch freed from the touches they visit, back
/// from an end on the border of either chart: see [`cut_at_touches`]. For
/// each end, head then tail, whether any were, or whether it lies on a
/// border visiting no touch at all. A branch walked out both ways from its
/// seed stops for the reason its one end did, and an end on a border clear
/// of every touch is where the walk left the chart, whatever stopped the
/// other: it heads for no touch either.
pub(super) fn freed_from_borders(
    branch: &Traced,
    visits: &mut [Option<usize>],
    touches: &[Touch],
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    near: f64,
    tol: Tolerances,
) -> [bool; 2] {
    let n = branch.points.len();
    let mut freed = [false; 2];
    if n == 0 || branch.closed() {
        return freed;
    }
    for (side, at) in [0, n - 1].into_iter().enumerate() {
        if !outside(a, branch.on_a[at], tol) && !outside(b, branch.on_b[at], tol) {
            continue;
        }
        let Some(t) = visits[at] else {
            freed[side] = true;
            continue;
        };
        let touch = &touches[t];
        let Some(h) = half_arm(touch, branch.points[at]) else {
            continue;
        };
        if branch.points[at].distance(touch.point) <= near {
            continue;
        }
        let inward: Vec<usize> = if at == 0 {
            (0..n).collect()
        } else {
            (0..n).rev().collect()
        };
        for i in inward {
            let p = branch.points[i];
            if visits[i] != Some(t)
                || p.distance(touch.point) <= near
                || half_arm(touch, p) != Some(h)
                || !on_arm((a, b), touch, (h, false), p, near, tol)
            {
                break;
            }
            visits[i] = None;
            freed[side] = true;
        }
    }
    freed
}

/// Branches cut at the touches as they were traced, where the touches stand
/// within a stride of the march of each other (see [`split_at_touches`]):
/// each touch found from the low-angle samples of transversal branches,
/// within `reach` of a sample, a branch cut there unless it runs straight
/// through, judged by the directions in and out, and the arcs either side
/// carried onto it along the straight way. Debris goes as
/// [`split_at_touches`] drops it: an arc joining two touches, or running
/// from a touch out to the border of a chart, is kept however short. Each
/// comes certified as [`split_at_touches`] certifies its arcs.
pub(super) fn split_at_touches_as_traced(
    found: Vec<Traced>,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> Vec<(Traced, bool)> {
    let near = options.chord.max(tol.confusion());
    let reach = near * 60.0;
    // How far a touch solved from a low-angle sample may stand from it and
    // still be the one it lies beside, and how far short of a touch any
    // stalled end is carried onto it: the walk gives out where the arms
    // close in on each other. A touch solved from a stalled end itself is
    // that end's at any distance: a stall cut back to the last sample
    // clear of it stands a whole stride short, and `carry_onto` follows
    // the section out to the touch or leaves the end where it is. Both ends
    // of a branch marked stalled count as stalled, though one may have
    // left the domain: `carry_onto` leaves such an end where it is unless
    // the section runs on from it, ahead, to the touch.
    let carry = reach * 40.0;

    let mut touches: Vec<Touch> = Vec::new();
    // The touch solved from each stalled branch's ends, head then tail.
    let mut solved_ends: Vec<[Option<usize>; 2]> = vec![[None; 2]; found.len()];
    for (branch, solved) in found.iter().zip(&mut solved_ends) {
        let n = branch.points.len();
        if n < 5 {
            continue;
        }
        // A branch with no transversal stretch to certify it may be
        // tangential contact rather than a curve, and only crossings are
        // taken from it, as [`find_touches`] takes them.
        let transversal = interior_is_transversal(branch, a, b, tol);
        let sines: Vec<f64> = (0..n)
            .map(|i| crossing_sine(a, b, branch.on_a[i], branch.on_b[i], tol))
            .collect();
        for i in 0..n {
            let low = sines[i] < BRANCH_POINT_SINE
                && (i == 0 || sines[i] <= sines[i - 1])
                && (i + 1 == n || sines[i] <= sines[i + 1]);
            if !low
                || touches
                    .iter()
                    .any(|t| t.point.distance(branch.points[i]) <= reach)
            {
                continue;
            }
            let Some(mut touch) = touching_point(a, b, branch.on_a[i], branch.on_b[i], tol) else {
                continue;
            };
            let form = touch_arms(a, b, touch.on_a, touch.on_b, tol);
            if form.is_some_and(|form| !crosses(&touch, form, near, tol)) {
                continue;
            }
            touch.blur = form.map_or(0.0, |form| blur_of(form.3, tol));
            touch.arms = form
                .map(|(arms, _, _, _)| arms)
                .filter(|&arms| apart(arms, reach, near));
            if !transversal && touch.arms.is_none() {
                continue;
            }
            let stall = branch.stopped == Stopped::Stalled && (i == 0 || i + 1 == n);
            if !stall && touch.point.distance(branch.points[i]) > carry {
                continue;
            }
            let at = touches
                .iter()
                .position(|t| t.point.distance(touch.point) <= reach)
                .unwrap_or_else(|| {
                    touches.push(touch);
                    touches.len() - 1
                });
            if stall {
                solved[usize::from(i != 0)] = Some(at);
            }
        }
    }
    if touches.is_empty() {
        return found.into_iter().map(|arc| (arc, false)).collect();
    }
    settle_touches(&mut touches, (a, b), options, tol);
    let mut out = Vec::with_capacity(found.len());
    for (branch, solved) in found.into_iter().zip(solved_ends) {
        let n = branch.points.len();
        // A branch with no transversal stretch is cut only at crossings.
        let transversal = interior_is_transversal(&branch, a, b, tol);
        let counts = |t: &Touch| transversal || t.arms.is_some();
        let touch_near = |p: Point, within: f64| -> Option<usize> {
            touches
                .iter()
                .position(|t| counts(t) && t.point.distance(p) <= within)
        };
        // A sample on a touch, or either end of a step passing within
        // reach of one: the walk strides straight arms in long steps.
        let mut visits: Vec<Option<usize>> = branch
            .points
            .iter()
            .map(|p| touch_near(*p, reach))
            .collect();
        for i in 1..n {
            let (p, q) = (branch.points[i - 1], branch.points[i]);
            if let Some(t) = touches
                .iter()
                .position(|t| counts(t) && distance_to_segment(t.point, p, q) <= reach)
            {
                visits[i - 1].get_or_insert(t);
                visits[i].get_or_insert(t);
            }
        }
        // The touch a stalled end heads for: one within `carry` of it, or
        // failing that the one solved from the end itself. Where the end
        // will not carry onto the first, the one solved from it is tried
        // next: two crossings a little apart on one branch both lie within
        // `carry` of its end, and the one listed first may be the other.
        let freed = freed_from_borders(&branch, &mut visits, &touches, (a, b), near, tol);
        let solved_at = |at: usize| if at == 0 { solved[0] } else { solved[1] };
        let end_touch = |at: usize| -> Option<usize> {
            if branch.stopped != Stopped::Stalled || freed[usize::from(at != 0)] {
                return None;
            }
            touch_near(branch.points[at], carry).or_else(|| solved_at(at))
        };
        let end_fallback = |at: usize| -> Option<usize> {
            let first = end_touch(at)?;
            solved_at(at).filter(|&t| t != first)
        };
        let stall_end = |at: usize| end_touch(at).is_some();
        let touched =
            visits.iter().any(Option::is_some) || (n > 0 && (stall_end(0) || stall_end(n - 1)));
        if !touched {
            out.push(branch);
            continue;
        }
        let closed = branch.closed();
        // A loop is read from the start of a pass over a touch, so each arc
        // of it is one run of samples clear of every touch.
        let order: Vec<usize> = if closed {
            let last = if branch.points[0].distance(branch.points[n - 1]) <= tol.confusion() {
                n - 1
            } else {
                n
            };
            let Some(first) =
                (0..last).find(|&i| visits[i].is_some() && visits[(i + last - 1) % last].is_none())
            else {
                out.push(branch);
                continue;
            };
            (0..last).map(|k| (first + k) % last).collect()
        } else {
            (0..n).collect()
        };
        let mut runs: Vec<(Vec<usize>, Option<usize>, Option<usize>)> = Vec::new();
        let mut run: Vec<usize> = Vec::new();
        let mut head: Option<usize> = None;
        for (k, &i) in order.iter().enumerate() {
            if let Some(t) = visits[i] {
                if !run.is_empty() {
                    runs.push((core::mem::take(&mut run), head, Some(t)));
                }
                head = Some(t);
                continue;
            }
            if run.is_empty() && k == 0 && !closed && stall_end(i) {
                head = end_touch(i);
            }
            run.push(i);
        }
        if !run.is_empty() {
            let last = *run.last().unwrap_or(&0);
            let tail = if closed {
                visits[order[0]]
            } else if last == n - 1 {
                end_touch(last)
            } else {
                None
            };
            runs.push((run, head, tail));
        }
        // A branch passing each touch once, running straight through,
        // stays whole: the crossing is the caller's to split like any other.
        let pass = |before: &[usize], after: &[usize], t: usize| -> (usize, bool) {
            (
                t,
                straight_through(&branch, (before, after), &touches[t], (a, b), tol),
            )
        };
        let mut passes: Vec<(usize, bool)> = runs
            .windows(2)
            .filter_map(|w| Some(pass(&w[0].0, &w[1].0, w[0].2?)))
            .collect();
        if closed
            && let (Some(last), Some(first)) = (runs.last(), runs.first())
            && let Some(t) = last.2
        {
            passes.push(pass(&last.0, &first.0, t));
        }
        let mut met: Vec<usize> = passes.iter().map(|(t, _)| *t).collect();
        met.sort_unstable();
        let repeated = met.windows(2).any(|w| w[0] == w[1]);
        let straight = passes.iter().all(|(_, s)| *s);
        let stopped_short = !closed && (stall_end(0) || stall_end(n - 1));
        if !repeated && straight && !stopped_short {
            out.push(branch);
            continue;
        }

        let mut arcs: Vec<Traced> = Vec::new();
        for (run, head, tail) in runs {
            let mut arc = Traced {
                points: run.iter().map(|&i| branch.points[i]).collect(),
                on_a: run.iter().map(|&i| branch.on_a[i]).collect(),
                on_b: run.iter().map(|&i| branch.on_b[i]).collect(),
                stopped: Stopped::Stalled,
            };
            if arc.points.len() < 2 {
                continue;
            }
            let fallback = |at_end: bool, index: usize| {
                (at_end && !closed && visits[index].is_none())
                    .then(|| end_fallback(index))
                    .flatten()
            };
            let ends = [
                (tail, false, fallback(run.last() == Some(&(n - 1)), n - 1)),
                (head, true, fallback(run.first() == Some(&0), 0)),
            ];
            for (end, at_head, next) in ends {
                if let Some(t) = end
                    && !carry_onto_straight(&mut arc, &touches[t], at_head, (a, b), options, tol)
                    && let Some(t) = next
                {
                    carry_onto_straight(&mut arc, &touches[t], at_head, (a, b), options, tol);
                }
            }
            let length: f64 = arc.points.windows(2).map(|w| w[0].distance(w[1])).sum();
            let m = arc.points.len();
            let on = |end: Option<usize>, at: usize| {
                end.is_some_and(|t| arc.points[at].distance(touches[t].point) <= tol.confusion())
            };
            let joins = head != tail && on(head, 0) && on(tail, m - 1);
            // A run freed back from the border of a chart, carried onto its
            // touch, is the section's arm from the touch out to the border.
            let to_border = (freed[0] && run.first() == Some(&0) && on(tail, m - 1))
                || (freed[1] && run.last() == Some(&(n - 1)) && on(head, 0));
            if !joins
                && !to_border
                && (length <= options.chord * 10.0 || (m < 4 && length <= near * 60.0))
            {
                continue;
            }
            if head.is_some() && head == tail && arc.points[0].distance(arc.points[m - 1]) <= reach
            {
                // The golden section of the samples, clear of the middle
                // and quarters where a symmetric loop's features fall.
                let middle = m * 382 / 1000;
                let part = |r: core::ops::RangeInclusive<usize>| Traced {
                    points: arc.points[r.clone()].to_vec(),
                    on_a: arc.on_a[r.clone()].to_vec(),
                    on_b: arc.on_b[r].to_vec(),
                    stopped: Stopped::Stalled,
                };
                arcs.push(part(0..=middle));
                arcs.push(part(middle..=m - 1));
            } else {
                arcs.push(arc);
            }
        }
        out.extend(arcs);
    }
    certified(out, &touches, (a, b), near, tol)
}

/// An arc's end carried onto a touch along the straight way there, as
/// [`split_at_touches_as_traced`] carries it: samples `reach` apart down to
/// the touch, then halving to an eighth of it, each corrected onto both
/// surfaces. Nothing is added where the end does not head for the touch,
/// a sample further than `reach` from it will not settle onto the section,
/// or the samples would take the arc past the points a branch may hold.
/// Whether the end was carried onto the touch.
fn carry_onto_straight(
    arc: &mut Traced,
    touch: &Touch,
    at_head: bool,
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    options: Marching,
    tol: Tolerances,
) -> bool {
    let reach = options.chord.max(tol.confusion()) * 60.0;
    let n = arc.points.len();
    let (end, inner) = if at_head { (0, 1) } else { (n - 1, n - 2) };
    let from = arc.points[end];
    let gap = touch.point - from;
    let distance = gap.magnitude();
    let mut added: Vec<Contact> = Vec::new();
    if distance > tol.confusion() {
        let along = gap * (1.0 / distance);
        let heading = from - arc.points[inner];
        if distance > reach && heading.dot(along) < 0.5 * heading.magnitude() {
            return false;
        }
        // How far short of the touch each sample stands: steps of at most
        // `reach` down to it, then halving to an eighth of it, where the
        // two arms still stand far enough apart to correct onto the right
        // one and the last straight stretch bows from the arm well inside
        // the chord.
        let pieces = (distance / reach).ceil().max(1.0);
        #[allow(clippy::cast_precision_loss)]
        if pieces + 3.0 > options.max_points.saturating_sub(n) as f64 {
            return false;
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let count = pieces as usize;
        #[allow(clippy::cast_precision_loss)]
        let mut short: Vec<f64> = (1..count)
            .map(|k| distance - distance * k as f64 / pieces)
            .collect();
        let mut close = short.last().copied().unwrap_or(distance).min(reach);
        for _ in 0..3 {
            close *= 0.5;
            short.push(close);
        }
        let (mut on_a, mut on_b) = (arc.on_a[end], arc.on_b[end]);
        for remaining in short {
            let s = distance - remaining;
            let guess = from + along * s;
            let start = [on_a.0, on_a.1, on_b.0, on_b.1];
            let settled = correct(
                a,
                b,
                start,
                guess,
                Some((from, along, s)),
                options.chord,
                tol,
            )
            .filter(|c| c.point.distance(guess) <= remaining * 0.5);
            let Some(c) = settled else {
                // Beside the touch a sample that will not settle is left
                // out; further off, the arm is not there to follow.
                if remaining > reach {
                    return false;
                }
                continue;
            };
            on_a = beside(a, c.on_a, on_a);
            on_b = beside(b, c.on_b, on_b);
            added.push(Contact {
                on_a,
                on_b,
                point: c.point,
            });
        }
        added.push(Contact {
            on_a: beside(a, touch.on_a, on_a),
            on_b: beside(b, touch.on_b, on_b),
            point: touch.point,
        });
    } else {
        let (on_a, on_b) = (arc.on_a[end], arc.on_b[end]);
        arc.points[end] = touch.point;
        arc.on_a[end] = beside(a, touch.on_a, on_a);
        arc.on_b[end] = beside(b, touch.on_b, on_b);
        return true;
    }
    if at_head {
        added.reverse();
        arc.points.splice(0..0, added.iter().map(|c| c.point));
        arc.on_a.splice(0..0, added.iter().map(|c| c.on_a));
        arc.on_b.splice(0..0, added.iter().map(|c| c.on_b));
    } else {
        arc.points.extend(added.iter().map(|c| c.point));
        arc.on_a.extend(added.iter().map(|c| c.on_a));
        arc.on_b.extend(added.iter().map(|c| c.on_b));
    }
    true
}

/// Whether every sample of `arc` clear of the touches lies within `near` of
/// `other`, the middles of its chords within twice that: the samples of
/// both lie on the section, and the chords of a march stand within its
/// chord of it, so a stretch traced twice holds each of its samples within
/// the chord of the other tracing. Across the section that is widened to
/// where the confusion distance leaves it: see [`across_section`]. Round a
/// touch the samples are cut and carried, not marched, and are not
/// weighed. A sample on an arm of a touch is held by a tracing that runs
/// along that arm past it: see [`ArmsOf`].
fn lies_on(
    arc: &Traced,
    other: &Traced,
    near: f64,
    (across, arms): (&[Across], &ArmsOf<'_>),
    touches: &[Touch],
) -> bool {
    let mut last = 0;
    let mut weighed = false;
    let mut reaches = Reaches::default();
    let mut on_arm = |k: usize| arms.held(k, other, &mut reaches);
    // The samples, and the middle of each chord: an arc joining two touches
    // may have no sample clear of them. A chord between two samples on one
    // arm of a touch lies on that arm.
    let m = arc.points.len();
    for k in 0..2 * m.saturating_sub(1) + 1 {
        let (p, within, wide, held) = if k < m {
            (arc.points[k], near, across[k], on_arm(k))
        } else {
            let i = k - m;
            let step = arc.points[i + 1] - arc.points[i];
            let length = step.magnitude();
            let way = (length > f64::MIN_POSITIVE).then(|| step * (1.0 / length));
            (
                arc.points[i] + step * 0.5,
                2.0 * near,
                (2.0 * across[i].0.min(across[i + 1].0), way),
                arms.same_arm(i, i + 1) && on_arm(i) && on_arm(i + 1),
            )
        };
        if touches.iter().any(|t| t.point.distance(p) <= t.radius) {
            continue;
        }
        if !held && !holds(other, p, within, wide, &mut last) {
            return false;
        }
        weighed = true;
    }
    weighed
}

/// Where each sample of an arc lies on a half arm of a touch: the touch,
/// the half, and how far out along it, for a sample past the touch's
/// radius and within its `around` that [`on_arm`] places on that half
/// within the band the confusion distance leaves round it.
///
/// Out along an arm `ρ` from a touch the surfaces part across it at a sine
/// of about `ρ · spread`, so a sample within the confusion distance `ε` of
/// both stands up to `ε / (ρ · spread)` off it, and the way the section
/// runs through a sample that far off turns by as much over `ρ`. Two
/// tracings of the arm near the touch stand apart across it by more than
/// the chord, and the way at either sample does not run along the other:
/// neither is held by the other, sample by sample. The section leaving a
/// touch along one half arm is one curve, so a sample on that half is a
/// point of any tracing that leaves the touch along it and runs out past
/// the sample: see [`ArmsOf::held`].
struct ArmsOf<'t> {
    samples: Vec<Option<(usize, usize, f64)>>,
    touches: &'t [Touch],
    surfaces: (&'t SurfaceGeometry, &'t SurfaceGeometry),
    near: f64,
    tol: Tolerances,
}

/// How far out along each half arm of a touch one tracing runs from the
/// touch: see [`ArmsOf::held`].
type Reaches = Vec<((usize, usize), f64)>;

impl<'t> ArmsOf<'t> {
    fn new(
        arc: &Traced,
        touches: &'t [Touch],
        surfaces: (&'t SurfaceGeometry, &'t SurfaceGeometry),
        near: f64,
        tol: Tolerances,
    ) -> Self {
        let mut out = Self {
            samples: Vec::new(),
            touches,
            surfaces,
            near,
            tol,
        };
        out.samples = arc.points.iter().map(|&p| out.on_arm(p)).collect();
        out
    }

    /// The touch, half arm and distance out along it a point lies on.
    fn on_arm(&self, p: Point) -> Option<(usize, usize, f64)> {
        self.touches.iter().enumerate().find_map(|(t, touch)| {
            let arms = touch.arms?;
            let r = p.distance(touch.point);
            if r <= touch.radius || r > touch.around {
                return None;
            }
            let h = half_arm(touch, p)?;
            on_arm(self.surfaces, touch, (h, true), p, self.near, self.tol)
                .then(|| (t, h, (p - touch.point).dot(halves(arms)[h])))
        })
    }

    /// Whether samples `i` and `j` lie on one half arm of one touch.
    fn same_arm(&self, i: usize, j: usize) -> bool {
        match (self.samples[i], self.samples[j]) {
            (Some((t, h, _)), Some((u, k, _))) => t == u && h == k,
            _ => false,
        }
    }

    /// Whether sample `k`, on a half arm of a touch, is a point of `other`:
    /// `other` ends on the touch, its samples from there on lie on that
    /// half, and it runs out along it at least as far as the sample, or on
    /// past the touch's `around` along it. `reaches` keeps how far `other`
    /// runs along each half it has been asked about.
    fn held(&self, k: usize, other: &Traced, reaches: &mut Reaches) -> bool {
        let Some((t, h, out)) = self.samples[k] else {
            return false;
        };
        if let Some(&(_, reach)) = reaches.iter().find(|(key, _)| *key == (t, h)) {
            return reach >= out;
        }
        let reach = self.reach(other, t, h);
        reaches.push(((t, h), reach));
        reach >= out
    }

    /// Whether `arc` leaves a touch where the section crosses itself along
    /// one of its arms: an end on the touch, and the first sample past its
    /// radius on a half arm of it, or where the arc ends within the radius,
    /// its sample furthest from the touch, clear of where an end counts as
    /// on it. There the second forms' difference `Q` is indefinite and the
    /// surfaces stand off each other by `xᵀ Q x / 2` round the touch, which
    /// changes sign across each arm: the arc is a crossing, however shallow
    /// the surfaces cross along it, and no contact runs along it.
    fn leaves(&self, arc: &Traced) -> bool {
        let m = arc.points.len();
        if m < 2 {
            return false;
        }
        [(0, true), (m - 1, false)]
            .into_iter()
            .any(|(end, forward)| {
                self.touches.iter().enumerate().any(|(t, touch)| {
                    if touch.arms.is_none()
                        || arc.points[end].distance(touch.point) > self.near.max(touch.blur)
                    {
                        return false;
                    }
                    let inward: Vec<usize> = if forward {
                        (1..m).collect()
                    } else {
                        (0..m - 1).rev().collect()
                    };
                    let from = |i: &usize| arc.points[*i].distance(touch.point);
                    let past = inward.iter().copied().find(|i| from(i) > touch.radius);
                    let Some(i) = past.or_else(|| {
                        inward
                            .iter()
                            .copied()
                            .max_by(|x, y| from(x).total_cmp(&from(y)))
                            .filter(|i| from(i) > self.near.max(touch.blur))
                    }) else {
                        return false;
                    };
                    let p = arc.points[i];
                    match self.samples[i] {
                        Some((u, _, _)) => u == t,
                        None => {
                            (from(&i) <= touch.radius || from(&i) > touch.around)
                                && half_arm(touch, p).is_some_and(|h| {
                                    on_arm(self.surfaces, touch, (h, true), p, self.near, self.tol)
                                })
                        }
                    }
                })
            })
    }

    /// How far `other` runs out from touch `t` along its half arm `h`.
    fn reach(&self, other: &Traced, t: usize, h: usize) -> f64 {
        let touch = &self.touches[t];
        let m = other.points.len();
        let mut reach: f64 = 0.0;
        if m < 2 {
            return reach;
        }
        let on_touch = |p: Point| p.distance(touch.point) <= self.near.max(touch.blur);
        for end in [0, m - 1] {
            if !on_touch(other.points[end]) {
                continue;
            }
            let inward: Vec<usize> = if end == 0 {
                (1..m).collect()
            } else {
                (0..m - 1).rev().collect()
            };
            for i in inward {
                let p = other.points[i];
                let r = p.distance(touch.point);
                if r <= touch.radius {
                    continue;
                }
                if r > touch.around {
                    if half_arm(touch, p) == Some(h) {
                        reach = reach.max(touch.around);
                    }
                    break;
                }
                match self.on_arm(p) {
                    Some((u, k, out)) if u == t && k == h => reach = reach.max(out),
                    _ => break,
                }
            }
        }
        reach
    }
}

/// Whether `arc`, lying wholly within one touch's radius from the touch
/// out to an end of its own, is a stretch of `other`: every sample of it
/// off the touch lies within `near` of `other`. [`lies_on`] weighs no
/// sample so near a touch; the arm's stretch out to a chart's border is
/// still traced twice where two walks reach it.
fn stub_lies_on(arc: &Traced, other: &Traced, near: f64, touches: &[Touch]) -> bool {
    let m = arc.points.len();
    if m < 2 {
        return false;
    }
    let Some(touch) = touches
        .iter()
        .find(|t| arc.points.iter().all(|p| t.point.distance(*p) <= t.radius))
    else {
        return false;
    };
    let on_touch = |p: Point| p.distance(touch.point) <= near.max(touch.blur);
    if on_touch(arc.points[0]) == on_touch(arc.points[m - 1]) {
        return false;
    }
    let mut last = 0;
    let mut weighed = false;
    for &p in &arc.points {
        if on_touch(p) {
            continue;
        }
        if !holds(other, p, near, (near, None), &mut last) {
            return false;
        }
        weighed = true;
    }
    weighed
}

/// How far across the section two tracings of it may stand apart at each
/// sample of `arc`, and the way the section runs there: a sample is held
/// to both surfaces within the confusion distance `ε`, and a step `δ`
/// across the section opens a gap of `δ · sin θ`, `θ` the angle the
/// surfaces cross at, so each tracing lies within `ε / sin θ` of the
/// section across it, and two within twice that of each other. Round a
/// touch, where the surfaces cross at no angle, and where the angle is not
/// known, no further than `near`.
fn across_section(
    arc: &Traced,
    touches: &[Touch],
    near: f64,
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    tol: Tolerances,
) -> Vec<Across> {
    (0..arc.points.len())
        .map(|k| {
            let p = arc.points[k];
            let none = (near, None);
            if touches.iter().any(|t| t.point.distance(p) <= t.radius) {
                return none;
            }
            let (Some(na), Some(nb)) = (
                normal_at(a, arc.on_a[k], tol),
                normal_at(b, arc.on_b[k], tol),
            ) else {
                return none;
            };
            let cross = na.cross(nb);
            let sine = cross.magnitude();
            if sine <= f64::MIN_POSITIVE {
                return none;
            }
            (
                near.max(2.0 * tol.confusion() / sine),
                Some(cross * (1.0 / sine)),
            )
        })
        .collect()
}

/// How far across the section a point may stand off another tracing of
/// it, and the way the section runs there: see [`across_section`].
type Across = (f64, Option<Vector>);

/// Whether two arcs are one stretch between the same two touches: each
/// runs from one touch to the other, leaving each along the same half arm.
fn same_stretch(arc: &Traced, other: &Traced, near: f64, touches: &[Touch]) -> bool {
    let ends = |t: &Traced| -> Option<[(usize, usize); 2]> {
        let m = t.points.len();
        if m < 2 {
            return None;
        }
        let at = |end: usize, beside: usize| {
            let k = touches
                .iter()
                .position(|u| u.point.distance(t.points[end]) <= near.max(u.blur))?;
            Some((k, half_arm(&touches[k], t.points[beside])?))
        };
        Some([at(0, 1)?, at(m - 1, m - 2)?])
    };
    match (ends(arc), ends(other)) {
        (Some(x), Some(y)) => x == y || x == [y[1], y[0]],
        _ => false,
    }
}

/// Whether a point lies within `near` of a polyline, or beside a segment
/// of it within `across` (see [`beside_segment`]), searched outward from
/// the segment `last` a point beside it was found on, which it is left at:
/// points taken in order along one curve are found at once.
fn holds(other: &Traced, p: Point, near: f64, across: Across, last: &mut usize) -> bool {
    let segments = other.points.len().saturating_sub(1);
    let on = |k: usize| beside_segment(p, other.points[k], other.points[k + 1], near, across);
    let from = (*last).min(segments.saturating_sub(1));
    let found = (0..segments).find_map(|step| {
        let ahead = from + step;
        if ahead < segments && on(ahead) {
            return Some(ahead);
        }
        let behind = from.checked_sub(step)?;
        on(behind).then_some(behind)
    });
    found.is_some_and(|k| {
        *last = k;
        true
    })
}

/// Whether a point lies within `near` of the segment `from`–`to`, or, with
/// the way `along` the section runs at the point, stands off a point of the
/// segment within `near` along the section and within `across` across it.
fn beside_segment(p: Point, from: Point, to: Point, near: f64, (across, along): Across) -> bool {
    if distance_to_segment(p, from, to) <= near {
        return true;
    }
    let Some(along) = along else {
        return false;
    };
    let step = to - from;
    let length = step.square_magnitude();
    if across <= near || length <= f64::MIN_POSITIVE {
        return false;
    }
    let t = (p - from).dot(step) / length;
    if !(0.0..=1.0).contains(&t) {
        return false;
    }
    let off = p - (from + step * t);
    let ahead = off.dot(along);
    ahead.abs() <= near && (off - along * ahead).magnitude() <= across
}

/// An arc with the samples at either end that lie on arcs already kept cut
/// back, so a stretch traced twice where two tracings overlap at their
/// ends is kept once. Each cut end keeps the last sample on the other
/// tracing, so the two still meet; an end on a touch is kept. A sample on
/// an arm of a touch lies on a tracing running along that arm past it: see
/// [`ArmsOf`].
fn trimmed(
    arc: Traced,
    kept: &[Traced],
    near: f64,
    (across, arms): (&[Across], &ArmsOf<'_>),
    touches: &[Touch],
) -> Option<Traced> {
    let n = arc.points.len();
    // How many samples from one end run along arcs already kept, and the
    // arc the last of them lies on: each sample is looked for first on the
    // arc the one before was found on.
    let run = |order: &mut dyn Iterator<Item = usize>| -> (usize, Option<usize>) {
        let mut holder: Option<(usize, usize)> = None;
        let mut count = 0;
        for k in order {
            let p = arc.points[k];
            if touches.iter().any(|t| t.point.distance(p) <= t.radius) {
                break;
            }
            let on = |other: &Traced, last: &mut usize| {
                holds(other, p, near, across[k], last)
                    || arms.held(k, other, &mut Reaches::default())
            };
            let on_holder =
                holder.and_then(|(o, mut last)| on(&kept[o], &mut last).then_some((o, last)));
            let found = on_holder.or_else(|| {
                kept.iter().enumerate().find_map(|(o, other)| {
                    let mut last = 0;
                    on(other, &mut last).then_some((o, last))
                })
            });
            if found.is_none() {
                break;
            }
            holder = found;
            count += 1;
        }
        (count, holder.map(|(o, _)| o))
    };
    let (head, head_on) = run(&mut (0..n));
    if head == n {
        return None;
    }
    let (back, tail_on) = run(&mut (0..n).rev());
    let tail = n - 1 - back;
    if tail < head {
        return None;
    }
    let mut points = arc.points[head..=tail].to_vec();
    let mut on_a = arc.on_a[head..=tail].to_vec();
    let mut on_b = arc.on_b[head..=tail].to_vec();
    // Each cut end meets the arc it ran along: at that arc's own end where
    // it stops between the last sample cut and the first kept, so the two
    // do not run side by side, and otherwise at the last sample cut.
    let meet = |cut: usize, kept_at: usize, on: Option<usize>| -> (Point, (f64, f64), (f64, f64)) {
        let (u, h) = (arc.points[kept_at], arc.points[cut]);
        on.map(|o| &kept[o])
            .and_then(|other| {
                let m = other.points.len();
                [0, m - 1]
                    .into_iter()
                    .filter(|&e| {
                        let wide = (across[cut].0.min(across[kept_at].0), across[cut].1);
                        beside_segment(other.points[e], u, h, near, wide)
                    })
                    .min_by(|&x, &y| {
                        other.points[x]
                            .distance(u)
                            .total_cmp(&other.points[y].distance(u))
                    })
                    .map(|e| (other.points[e], other.on_a[e], other.on_b[e]))
            })
            .unwrap_or((h, arc.on_a[cut], arc.on_b[cut]))
    };
    if head > 0 {
        let (p, a, b) = meet(head - 1, head, head_on);
        points.insert(0, p);
        on_a.insert(0, a);
        on_b.insert(0, b);
    }
    if tail + 1 < n {
        let (p, a, b) = meet(tail + 1, tail, tail_on);
        points.push(p);
        on_a.push(a);
        on_b.push(b);
    }
    (points.len() >= 2).then_some(Traced {
        points,
        on_a,
        on_b,
        stopped: arc.stopped,
    })
}

/// Whether the chord of the step into sample `i` stands off the section by
/// more than `near`, the march's chord: the section solved on the plane
/// square to the step through its middle lands further than that from it,
/// so the step is not one along a single curve.
pub(super) fn off_the_section(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    branch: &Traced,
    i: usize,
    near: f64,
    tol: Tolerances,
) -> bool {
    let (p, q) = (branch.points[i - 1], branch.points[i]);
    let step = q - p;
    let length = step.magnitude();
    if length <= near {
        return false;
    }
    let along = step / length;
    let middle = p + step * 0.5;
    let start = [
        branch.on_a[i - 1].0,
        branch.on_a[i - 1].1,
        branch.on_b[i - 1].0,
        branch.on_b[i - 1].1,
    ];
    correct(
        a,
        b,
        start,
        middle,
        Some((p, along, length * 0.5)),
        near,
        tol,
    )
    .is_some_and(|landed| landed.point.distance(middle) > near)
}

/// The four halves of a touch's two arms, each pointing away from it.
pub(super) fn halves(arms: [Vector; 2]) -> [Vector; 4] {
    [arms[0], -arms[0], arms[1], -arms[1]]
}

/// The half arm of a touch a point lies nearest in direction.
pub(super) fn half_arm(touch: &Touch, p: Point) -> Option<usize> {
    let arms = touch.arms?;
    let off = p - touch.point;
    let m = off.magnitude();
    if m <= f64::MIN_POSITIVE {
        return None;
    }
    halves(arms)
        .iter()
        .map(|h| h.dot(off) / m)
        .enumerate()
        .max_by(|x, y| x.1.total_cmp(&y.1))
        .map(|(k, _)| k)
}

/// Whether a sample of the section lies on one half arm of a touch: the
/// section solved on the plane square to that half through the sample
/// lands within `near` of it, or, with `blurred`, within the band the
/// confusion distance leaves round the arm that far out.
pub(super) fn on_arm(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    touch: &Touch,
    (half, blurred): (usize, bool),
    p: Point,
    near: f64,
    tol: Tolerances,
) -> bool {
    let Some(arms) = touch.arms else {
        return false;
    };
    let direction = halves(arms)[half];
    let along = (p - touch.point).dot(direction);
    let band = if blurred && touch.spread > 0.0 {
        tol.confusion() / (along * touch.spread)
    } else {
        0.0
    };
    along > 0.0
        && along_arm(a, b, touch, direction, (along, arms), near, tol)
            .is_some_and(|landed| landed.point.distance(p) <= near.max(band))
}
