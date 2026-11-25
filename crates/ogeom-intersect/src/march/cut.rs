//! One branch cut at the touches it visits.

use super::*;

/// A branch as [`cut_at_touches`] leaves it.
pub(super) struct Cut {
    /// The branch whole, or the arcs cut from it.
    pub(super) pieces: Vec<Traced>,
    /// Whether it was cut.
    pub(super) cut: bool,
    /// The touches it runs into or through and the half arm it runs along
    /// at each, by index into the touches and into [`halves`].
    pub(super) covered: Vec<(usize, usize)>,
}

/// A branch with each touch a step of it passes over, both ends of the
/// step clear of the touch's radius, put in as a sample where the step
/// passes it; with each sample, the touch it was put in for.
///
/// A step passes over a touch it passes within the touch's radius, or
/// within its `snap`, held to its `around`: the straight way from the
/// step to a touch that near lies on both surfaces within the confusion
/// distance, so the walk, which places a sample across the section only
/// to within that distance over the sine of the crossing angle, may stride
/// past the touch that far off it.
fn with_touches_on_steps(
    branch: Traced,
    touches: &[Touch],
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    counts: &dyn Fn(&Touch) -> bool,
) -> (Traced, Vec<Option<usize>>) {
    let n = branch.points.len();
    let mut out = Traced {
        points: Vec::with_capacity(n),
        on_a: Vec::with_capacity(n),
        on_b: Vec::with_capacity(n),
        stopped: branch.stopped,
    };
    let mut on_touch = Vec::with_capacity(n);
    for i in 0..n {
        if i > 0 {
            let (p, q) = (branch.points[i - 1], branch.points[i]);
            let step = q - p;
            let length2 = step.square_magnitude();
            let mut passed: Vec<(f64, usize)> = touches
                .iter()
                .enumerate()
                .filter(|(_, t)| {
                    let over = (t.point - p).dot(step);
                    counts(t)
                        && t.point.distance(p) > t.radius
                        && t.point.distance(q) > t.radius
                        && over > 0.0
                        && over < length2
                        && distance_to_segment(t.point, p, q) <= t.radius.max(t.snap.min(t.around))
                })
                .map(|(k, t)| ((t.point - p).dot(step) / length2.max(f64::MIN_POSITIVE), k))
                .collect();
            passed.sort_by(|x, y| x.0.total_cmp(&y.0));
            for (_, k) in passed {
                let t = &touches[k];
                out.points.push(t.point);
                out.on_a.push(beside(a, t.on_a, branch.on_a[i - 1]));
                out.on_b.push(beside(b, t.on_b, branch.on_b[i - 1]));
                on_touch.push(Some(k));
            }
        }
        out.points.push(branch.points[i]);
        out.on_a.push(branch.on_a[i]);
        out.on_b.push(branch.on_b[i]);
        on_touch.push(None);
    }
    (out, on_touch)
}

/// One branch cut at the touches it visits, as [`split_at_touches`] cuts
/// it. `solved` holds the touches solved from its ends, head then tail.
pub(super) fn cut_at_touches(
    branch: Traced,
    solved: [Option<usize>; 2],
    touches: &[Touch],
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    options: Marching,
    tol: Tolerances,
) -> Cut {
    let near = options.chord.max(tol.confusion());
    // A branch with no transversal stretch to certify it is cut only at
    // crossings, as only crossings are found from it.
    let transversal = interior_is_transversal(&branch, a, b, tol);
    let counts = |t: &Touch| transversal || t.arms.is_some();
    let touch_near = |p: Point| -> Option<usize> {
        (0..touches.len())
            .filter(|&t| counts(&touches[t]) && touches[t].point.distance(p) <= touches[t].around)
            .min_by(|&x, &y| {
                touches[x]
                    .point
                    .distance(p)
                    .total_cmp(&touches[y].point.distance(p))
            })
    };
    // A step passing over touches, its ends clear of them, gets each as a
    // sample where it passes: a straight arm walked in long strides passes
    // several touches a step apart, and each is cut at on its own.
    let (branch, on_touch) = with_touches_on_steps(branch, touches, (a, b), &counts);
    let n = branch.points.len();
    // A sample on a touch, or either end of a step passing within its
    // radius: the walk strides straight arms in long steps.
    let mut visits: Vec<Option<usize>> = branch
        .points
        .iter()
        .map(|p| {
            touches
                .iter()
                .position(|t| counts(t) && t.point.distance(*p) <= t.radius)
        })
        .collect();
    for i in 1..n {
        let (p, q) = (branch.points[i - 1], branch.points[i]);
        if let Some(t) = touches
            .iter()
            .position(|t| counts(t) && distance_to_segment(t.point, p, q) <= t.radius)
        {
            visits[i - 1].get_or_insert(t);
            visits[i].get_or_insert(t);
        }
    }
    // A step from one arm of a touch onto another passes over the touch
    // however far off it the step runs: a walk turning a corner there cuts
    // it by as much as its step. One stepping from an arm of one touch onto
    // an arm of another passes over both where its chord leaves the
    // section, as a step along one arm through the other touch does not.
    // Both ends must lie on the arms, so a stretch of the section near a
    // touch but on neither is not taken for one.
    let arm_of = |t: usize, p: Point| -> Option<usize> {
        let touch = &touches[t];
        if !counts(touch) || touch.point.distance(p) > touch.around {
            return None;
        }
        let h = half_arm(touch, p)?;
        on_arm((a, b), touch, (h, true), p, near, tol).then_some(h)
    };
    for i in 1..n {
        if visits[i - 1].is_some() || visits[i].is_some() {
            continue;
        }
        let (p, q) = (branch.points[i - 1], branch.points[i]);
        let mut ends: [Option<(usize, usize)>; 2] = [None; 2];
        let mut over = None;
        for t in 0..touches.len() {
            let (hp, hq) = (arm_of(t, p), arm_of(t, q));
            if let (Some(hp), Some(hq)) = (hp, hq) {
                if hp != hq {
                    over = Some((t, t));
                }
                ends = [Some((t, hp)), Some((t, hq))];
                break;
            }
            ends[0] = ends[0].or(hp.map(|h| (t, h)));
            ends[1] = ends[1].or(hq.map(|h| (t, h)));
        }
        if over.is_none()
            && let [Some((tp, _)), Some((tq, _))] = ends
            && tp != tq
            && off_the_section((a, b), &branch, i, near, tol)
        {
            over = Some((tp, tq));
        }
        if let Some((tp, tq)) = over {
            visits[i - 1] = Some(tp);
            visits[i] = Some(tq);
        }
    }
    // An end of an open branch on the border of either chart is an end of
    // the section, however near a touch it stands: the samples from it
    // back along the half arm it lies on, each solved onto that arm, are
    // the arm's and no visit, and the stretch from the touch out to the
    // border is cut as an arc of its own, however short: an arm leaving
    // the chart within the radius has no sample clear of the touch. Such
    // an end heads for no touch.
    let freed = freed_from_borders(&branch, &mut visits, touches, (a, b), near, tol);
    // The touch a stalled end heads for: one within `carry` of it, or
    // failing that the one solved from the end itself. Where the end
    // will not carry onto the first, the one solved from it is tried
    // next: two crossings a little apart on one branch both lie within
    // `carry` of its end, and the one listed first may be the other.
    let solved_at = |at: usize| if at == 0 { solved[0] } else { solved[1] };
    let end_touch = |at: usize| -> Option<usize> {
        if branch.stopped != Stopped::Stalled || freed[usize::from(at != 0)] {
            return None;
        }
        touch_near(branch.points[at]).or_else(|| solved_at(at))
    };
    let end_fallback = |at: usize| -> Option<usize> {
        let first = end_touch(at)?;
        solved_at(at).filter(|&t| t != first)
    };
    let stall_end = |at: usize| end_touch(at).is_some();
    let touched =
        visits.iter().any(Option::is_some) || (n > 0 && (stall_end(0) || stall_end(n - 1)));
    if !touched {
        return Cut {
            pieces: vec![branch],
            cut: false,
            covered: Vec::new(),
        };
    }
    let closed = branch.closed();
    // The half arm a branch leaves each touch along: the sample beside each
    // run of samples on it, and an end carried onto it.
    let mut covered: Vec<(usize, usize)> = Vec::new();
    for i in 0..n {
        let Some(t) = visits[i] else { continue };
        let beside = [
            i.checked_sub(1).or_else(|| closed.then_some(n - 1)),
            (i + 1 < n).then_some(i + 1).or_else(|| closed.then_some(0)),
        ];
        for j in beside.into_iter().flatten() {
            // The sample beside, clear of the touch or on another one a
            // step put beside it.
            if visits[j].is_none_or(|u| u != t && on_touch[j].is_some())
                && let Some(h) = half_arm(&touches[t], branch.points[j])
            {
                covered.push((t, h));
            }
        }
    }
    for at in [0, n - 1] {
        if !closed
            && visits[at].is_none()
            && let Some(t) = end_touch(at)
            && let Some(h) = half_arm(&touches[t], branch.points[at])
        {
            covered.push((t, h));
        }
    }
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
            return Cut {
                pieces: vec![branch],
                cut: false,
                covered,
            };
        };
        (0..last).map(|k| (first + k) % last).collect()
    } else {
        (0..n).collect()
    };
    let mut runs: Vec<(Vec<usize>, Option<usize>, Option<usize>)> = Vec::new();
    let mut run: Vec<usize> = Vec::new();
    let mut head: Option<usize> = None;
    let mut before: Option<usize> = None;
    for (k, &i) in order.iter().enumerate() {
        if let Some(t) = visits[i] {
            if !run.is_empty() {
                runs.push((core::mem::take(&mut run), head, Some(t)));
            } else if let Some(j) = before
                && let (Some(u), Some(v)) = (on_touch[j], on_touch[i])
                && u != v
            {
                // Two touches a step put side by side: the stretch between
                // them is that much of the step.
                runs.push((vec![j, i], Some(u), Some(v)));
            }
            head = Some(t);
            before = Some(i);
            continue;
        }
        before = Some(i);
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
    // Runs either side of a step over two touches meet neither straight.
    let mut passes: Vec<(usize, bool)> = runs
        .windows(2)
        .filter_map(|w| {
            let t = w[0].2?;
            Some(if w[1].1 == Some(t) {
                pass(&w[0].0, &w[1].0, t)
            } else {
                (t, false)
            })
        })
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
    // An open branch ending on a touch or short of one is cut there too:
    // the samples it wanders round the touch on go, and its end is carried
    // onto the touch.
    let ends_at_touch = !closed
        && n > 0
        && (stall_end(0) || stall_end(n - 1) || visits[0].is_some() || visits[n - 1].is_some());
    if !repeated && straight && !ends_at_touch {
        return Cut {
            pieces: vec![branch],
            cut: false,
            covered,
        };
    }

    let mut arcs: Vec<Traced> = Vec::new();
    for (run, head, tail) in runs {
        let mut arc = Traced {
            points: run.iter().map(|&i| branch.points[i]).collect(),
            on_a: run.iter().map(|&i| branch.on_a[i]).collect(),
            on_b: run.iter().map(|&i| branch.on_b[i]).collect(),
            stopped: Stopped::Stalled,
        };
        // A run between two touches is the stretch between two crossings,
        // whatever few samples a long stride left on it: one sample is
        // joined to both touches along the section.
        let between = match (head, tail) {
            (Some(h), Some(t)) if h != t => Some((h, t)),
            _ => None,
        };
        if arc.points.len() == 1
            && let Some((h, t)) = between
        {
            let touch_contact = |touch: &Touch, near_to: Contact| Contact {
                on_a: beside(a, touch.on_a, near_to.on_a),
                on_b: beside(b, touch.on_b, near_to.on_b),
                point: touch.point,
            };
            let sample = Contact {
                on_a: arc.on_a[0],
                on_b: arc.on_b[0],
                point: arc.points[0],
            };
            let (from, to) = (
                touch_contact(&touches[h], sample),
                touch_contact(&touches[t], sample),
            );
            let mut budget = options.max_points;
            let joined = fill((a, b), from, sample, near, &mut budget, tol).zip(fill(
                (a, b),
                sample,
                to,
                near,
                &mut budget,
                tol,
            ));
            if let Some((before, after)) = joined {
                let all: Vec<Contact> = core::iter::once(from)
                    .chain(before)
                    .chain([sample])
                    .chain(after)
                    .chain([to])
                    .collect();
                arcs.push(Traced {
                    points: all.iter().map(|c| c.point).collect(),
                    on_a: all.iter().map(|c| c.on_a).collect(),
                    on_b: all.iter().map(|c| c.on_b).collect(),
                    stopped: Stopped::Stalled,
                });
            }
            continue;
        }
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
            let Some(t) = end else { continue };
            // Samples by the touch that stand off its arms by more than the
            // chord are the walk's wandering round it, where the surfaces
            // stand within the confusion distance of each other across the
            // corner, and go before the end is carried on.
            let mut trial = arc.clone();
            if touches[t].arms.is_some() {
                while trial.points.len() > 2 {
                    let at = if at_head { 0 } else { trial.points.len() - 1 };
                    let p = trial.points[at];
                    let wanders = p.distance(touches[t].point) <= touches[t].around
                        && half_arm(&touches[t], p).is_some_and(|h| {
                            !on_arm((a, b), &touches[t], (h, false), p, near, tol)
                        });
                    if !wanders {
                        break;
                    }
                    trial.points.remove(at);
                    trial.on_a.remove(at);
                    trial.on_b.remove(at);
                }
            }
            let carried = carry_onto(&mut trial, &touches[t], at_head, (a, b), options, tol)
                || next.is_some_and(|t| {
                    carry_onto(&mut trial, &touches[t], at_head, (a, b), options, tol)
                });
            if carried {
                arc = trial;
                continue;
            }
            // An end that will not carry onto the touch keeps the samples
            // the walk took past it, as the walk left them.
            let at = if at_head { run.first() } else { run.last() };
            let Some(&at) = at else { continue };
            let Some(k) = order.iter().position(|&i| i == at) else {
                continue;
            };
            let m = order.len();
            let step = |k: usize| -> Option<usize> {
                if at_head {
                    k.checked_sub(1).or_else(|| closed.then_some(m - 1))
                } else {
                    (k + 1 < m).then_some(k + 1).or_else(|| closed.then_some(0))
                }
            };
            let mut j = step(k);
            let mut seen = 0;
            while let Some(jj) = j
                && visits[order[jj]].is_some()
                && seen < m
            {
                let i = order[jj];
                if at_head {
                    arc.points.insert(0, branch.points[i]);
                    arc.on_a.insert(0, branch.on_a[i]);
                    arc.on_b.insert(0, branch.on_b[i]);
                } else {
                    arc.points.push(branch.points[i]);
                    arc.on_a.push(branch.on_a[i]);
                    arc.on_b.push(branch.on_b[i]);
                }
                j = step(jj);
                seen += 1;
            }
        }
        let length: f64 = arc.points.windows(2).map(|w| w[0].distance(w[1])).sum();
        let joins = between.is_some_and(|(h, t)| {
            let m = arc.points.len();
            arc.points[0].distance(touches[h].point) <= tol.confusion()
                && arc.points[m - 1].distance(touches[t].point) <= tol.confusion()
        });
        // A run freed back from the border of a chart, carried onto its
        // touch, is the section's arm from the touch out to the border.
        let to_border = {
            let m = arc.points.len();
            let carried = |end: Option<usize>, at: usize| {
                end.is_some_and(|t| arc.points[at].distance(touches[t].point) <= tol.confusion())
            };
            (freed[0] && run.first() == Some(&0) && carried(tail, m - 1))
                || (freed[1] && run.last() == Some(&(n - 1)) && carried(head, 0))
        };
        if !joins
            && !to_border
            && (length <= options.chord * 10.0 || (arc.points.len() < 4 && length <= near * 60.0))
        {
            continue;
        }
        let m = arc.points.len();
        if let (Some(h), Some(t)) = (head, tail)
            && h == t
            && arc.points[0].distance(arc.points[m - 1]) <= touches[h].radius
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
    // The arcs cut from it run along the arms they leave the touches by,
    // where they reach the touches.
    let mut covered = Vec::new();
    for arc in &arcs {
        let m = arc.points.len();
        for (end, beside) in [(0, 1), (m - 1, m.saturating_sub(2))] {
            if m < 2 {
                continue;
            }
            if let Some(t) = touches
                .iter()
                .position(|t| t.point.distance(arc.points[end]) <= tol.confusion())
                && let Some(h) = half_arm(&touches[t], arc.points[beside])
            {
                covered.push((t, h));
            }
        }
    }
    Cut {
        pieces: arcs,
        cut: true,
        covered,
    }
}
