//! Each touch settled: its radius and arms, the section along them, samples beside it.

use super::*;

/// Each touch's radius, arms and snap, once all are known: a sample counts
/// as on the touch within the reach of [`split_at_touches`], and a stalled
/// end or a step between its arms as heading for it within the `carry` of
/// [`find_touches`], each held to half the way to the nearest other touch
/// where that is less, so nothing counts as on or heading for two; the
/// arms are kept where they stand apart at the radius.
pub(super) fn settle_touches(
    touches: &mut [Touch],
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    options: Marching,
    tol: Tolerances,
) {
    let near = options.chord.max(tol.confusion());
    let reach = near * 60.0;
    let points: Vec<Point> = touches.iter().map(|t| t.point).collect();
    for (i, touch) in touches.iter_mut().enumerate() {
        let nearest = points
            .iter()
            .enumerate()
            .filter(|&(j, _)| j != i)
            .map(|(_, p)| p.distance(touch.point))
            .fold(f64::INFINITY, f64::min);
        let form = touch_arms(a, b, touch.on_a, touch.on_b, tol);
        let spread = form.map_or(0.0, |(_, spread, _, _)| spread);
        touch.snap = form.map_or(0.0, |(_, _, bend, _)| {
            if bend > 0.0 {
                (2.0 * tol.confusion() / bend).sqrt()
            } else {
                0.0
            }
        });
        touch.radius = reach.min(nearest * 0.5);
        touch.around = (reach * 40.0).min(nearest * 0.5);
        touch.spread = spread;
        touch.arms = form
            .map(|(arms, _, _, _)| arms)
            .filter(|&arms| apart(arms, touch.radius, near));
    }
}

/// The length of a traced polyline.
pub(super) fn arc_length(arc: &Traced) -> f64 {
    arc.points.windows(2).map(|w| w[0].distance(w[1])).sum()
}

/// The parameters of both surfaces moved from a touch by `step` in the
/// shared tangent plane, each through its own chart's first derivatives.
fn out_along(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    touch: &Touch,
    step: Vector,
    tol: Tolerances,
) -> Option<[f64; 4]> {
    let moved = |surface: &SurfaceGeometry, (u, v): (f64, f64)| -> Option<(f64, f64)> {
        let (su, sv) = surface.d1_at(u, v, tol).ok()?;
        let (e, f, g) = (su.dot(su), su.dot(sv), sv.dot(sv));
        let det = e.mul_add(g, -(f * f));
        if det <= f64::MIN_POSITIVE {
            return None;
        }
        let (p, q) = (su.dot(step), sv.dot(step));
        Some((
            u + g.mul_add(p, -(f * q)) / det,
            v + e.mul_add(q, -(f * p)) / det,
        ))
    };
    let (ua, va) = moved(a, touch.on_a)?;
    let (ub, vb) = moved(b, touch.on_b)?;
    Some([ua, va, ub, vb])
}

/// The section a distance `reach` out from a touch along one half of an
/// arm: corrected onto both surfaces on the plane square to the arm there,
/// from parameters moved out along it on each chart. `None` where the
/// correction does not settle nearer that arm than the other.
pub(super) fn along_arm(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    touch: &Touch,
    half: Vector,
    (reach, arms): (f64, [Vector; 2]),
    within: f64,
    tol: Tolerances,
) -> Option<Contact> {
    let [ua, va, ub, vb] = out_along((a, b), touch, half * reach, tol)?;
    let guess = touch.point + half * reach;
    let landed = correct(
        a,
        b,
        [ua, va, ub, vb],
        guess,
        Some((touch.point, half, reach)),
        within,
        tol,
    )?;
    // The other arm crosses that plane `reach · tan α` off this one.
    let apart = arms[0].cross(arms[1]).magnitude();
    (landed.point.distance(guess) * 2.0 < reach * apart).then_some(landed)
}

/// Where a solve for the section lands: see [`settled_from`].
#[derive(Debug, Clone, Copy)]
pub(super) enum Landing {
    /// On the section, within the distance asked.
    At(Contact),
    /// Short of the section, a chart's border in the way: the parameter
    /// held to the chart, and the bound it was held to.
    Leaves(usize, f64),
    /// The solve did not settle.
    Unsettled,
}

/// One more equation for the section's three, in the four parameters of
/// both surfaces: its value and gradient, given the parameters, the first
/// surface's point and its derivatives there.
type Constraint<'c> = &'c dyn Fn(&[f64; 4], Point, [Vector; 2]) -> (f64, [f64; 4]);

/// The section's three equations and `constraint`, solved from `start` to
/// the residual that places the point within `within` of where the section
/// meets the constraint ([`crate::walk::settling_residual`]), stopping no
/// later than [`correct`] does. The landing is kept where the Jacobian
/// changes over the whole move, seen through the start's inverse, by no
/// more than a half, the condition [`crate::walk::settle`] keeps a further
/// solve to: the iteration ran to the one root beside its start. A solve
/// that cannot reach the section without leaving a chart says which bound
/// held it.
pub(super) fn settled_from(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    start: [f64; 4],
    constraint: Constraint<'_>,
    within: f64,
    tol: Tolerances,
) -> Landing {
    let measure = |x: &[f64; 4]| -> Option<Measured> {
        let (ua, va) = clamp(a, x[0], x[1]);
        let (ub, vb) = clamp(b, x[2], x[3]);
        let (pa, au, av) = a.point_d1_at(ua, va, tol).ok()?;
        let (pb, bu, bv) = b.point_d1_at(ub, vb, tol).ok()?;
        let gap = pa - pb;
        let (value, row) = constraint(x, pa, [au, av]);
        let residual = [gap.x, gap.y, gap.z, value];
        let jacobian = [
            [au.x, av.x, -bu.x, -bv.x],
            [au.y, av.y, -bu.y, -bv.y],
            [au.z, av.z, -bu.z, -bv.z],
            row,
        ];
        Some((*x, jacobian, residual, [au, av]))
    };
    let system = |x: &[f64; 4]| {
        measure(x).map_or(([f64::INFINITY; 4], [[0.0; 4]; 4]), |(_, j, r, _)| (r, j))
    };
    let Some((_, first, _, [au, av])) = measure(&start) else {
        return Landing::Unsettled;
    };
    let gradient = [au, av, Vector::ZERO, Vector::ZERO];
    let residual =
        crate::walk::settling_residual(&first, &gradient, within, tol.confusion() * 0.01);
    let criteria = solve::Criteria {
        residual,
        step: 0.0,
        max_iterations: 50,
    };
    let Ok((x, norm, _, _)) = solve::newton_system_fixed(system, start, criteria) else {
        return Landing::Unsettled;
    };
    let held = [clamp(a, x[0], x[1]), clamp(b, x[2], x[3])];
    let held = [held[0].0, held[0].1, held[1].0, held[1].1];
    if norm > residual {
        // A periodic parameter wraps rather than being held.
        let bounded = |k: usize| match k {
            0 => !a.is_periodic_u(),
            1 => !a.is_periodic_v(),
            2 => !b.is_periodic_u(),
            _ => !b.is_periodic_v(),
        };
        return (0..4)
            .find(|&k| bounded(k) && held[k] != x[k])
            .map_or(Landing::Unsettled, |k| Landing::Leaves(k, held[k]));
    }
    let Some((_, last, _, _)) = measure(&x) else {
        return Landing::Unsettled;
    };
    let mut drift = 0.0;
    for column in 0..4 {
        let change: [f64; 4] = core::array::from_fn(|row| last[row][column] - first[row][column]);
        let Some(moved) = solve::linear_fixed(first, change) else {
            return Landing::Unsettled;
        };
        drift += moved.iter().map(|v| v * v).sum::<f64>();
    }
    if drift.sqrt() > 0.5 {
        return Landing::Unsettled;
    }
    a.point_at(held[0], held[1], tol)
        .map_or(Landing::Unsettled, |point| {
            Landing::At(Contact {
                on_a: (held[0], held[1]),
                on_b: (held[2], held[3]),
                point,
            })
        })
}

/// The section `along` out from a touch on one half of an arm, solved on
/// the plane square to that half there to within `within` of the section.
///
/// Out along an arm `ρ` from a touch the surfaces cross at a sine of about
/// `ρ · spread`, so a point within the confusion distance of both may
/// stand `ε / (ρ · spread)` off the arm, and a solve stopped at that
/// residual places its point anywhere across that band. This one starts
/// on the arm's straight line, the parameters moved out along it on each
/// chart, off the arm by its bow, `O(ρ²)`, against a band of `O(1/ρ)`, and
/// is solved down as [`settled_from`] solves. The point is kept where it
/// stands nearer this arm than the other, which crosses the plane
/// `ρ · tan α` off it.
pub(super) fn arm_point(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    touch: &Touch,
    half: Vector,
    (along, arms): (f64, [Vector; 2]),
    within: f64,
    tol: Tolerances,
) -> Landing {
    let Some(start) = out_along((a, b), touch, half * along, tol) else {
        return Landing::Unsettled;
    };
    let plane = |_: &[f64; 4], pa: Point, [au, av]: [Vector; 2]| {
        (
            (pa - touch.point).dot(half) - along,
            [au.dot(half), av.dot(half), 0.0, 0.0],
        )
    };
    let landed = settled_from((a, b), start, &plane, within, tol);
    let guess = touch.point + half * along;
    let apart = arms[0].cross(arms[1]).magnitude();
    match landed {
        Landing::At(c) if c.point.distance(guess) * 2.0 >= along * apart => Landing::Unsettled,
        other => other,
    }
}

/// Where an arm leaving a touch along `half` leaves either chart, short of
/// `beyond` out along it, where [`arm_point`] has it leave: the point of
/// the section on the chart's border there, solved from the last point of
/// the arm found, within the confusion distance of the way out. `None`
/// where a solve on the way does not settle, or the arm leaves within the
/// confusion distance of the touch.
fn arm_exit(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    touch: &Touch,
    half: Vector,
    (beyond, arms): (f64, [Vector; 2]),
    near: f64,
    tol: Tolerances,
) -> Option<Contact> {
    let Landing::Leaves(k, at) = arm_point((a, b), touch, half, (beyond, arms), near, tol) else {
        return None;
    };
    let mut inside: Option<(Contact, f64)> = None;
    let (mut out, mut bound) = (beyond, (k, at));
    loop {
        let from = inside.map_or(0.0, |(_, along)| along);
        if out - from <= tol.confusion() {
            break;
        }
        let middle = f64::midpoint(from, out);
        match arm_point((a, b), touch, half, (middle, arms), near, tol) {
            Landing::At(c) => inside = Some((c, middle)),
            Landing::Leaves(k, at) => (out, bound) = (middle, (k, at)),
            Landing::Unsettled => return None,
        }
    }
    let ((c, _), (k, at)) = (inside?, bound);
    let mut start = [c.on_a.0, c.on_a.1, c.on_b.0, c.on_b.1];
    start[k] = at;
    let mut row = [0.0; 4];
    row[k] = 1.0;
    let border = |x: &[f64; 4], _: Point, _: [Vector; 2]| (x[k] - at, row);
    match settled_from((a, b), start, &border, near, tol) {
        Landing::At(end) => Some(end),
        _ => None,
    }
}

/// An arm leaving a chart within `beyond` of its touch, where no walk can
/// start: the touch, the points of the arm the chord between neighbours
/// calls for, each solved on it by [`arm_point`], and where it leaves
/// ([`arm_exit`]). `None` where a solve does not settle, the arm leaves
/// within the chord of the touch, or it would take more than `budget`
/// points.
pub(super) fn arm_to_border(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    touch: &Touch,
    half: Vector,
    (beyond, arms): (f64, [Vector; 2]),
    (near, budget): (f64, usize),
    tol: Tolerances,
) -> Option<Traced> {
    let exit = arm_exit((a, b), touch, half, (beyond, arms), near, tol)?;
    let length = (exit.point - touch.point).dot(half);
    if length <= near {
        return None;
    }
    let start = Contact {
        on_a: touch.on_a,
        on_b: touch.on_b,
        point: touch.point,
    };
    // Each chord between neighbours is halved where the arm solved at its
    // middle stands further than `near` off it.
    let mut points = vec![(start, 0.0), (exit, length)];
    let mut k = 0;
    while k + 1 < points.len() {
        let ((p, from), (q, to)) = (points[k], points[k + 1]);
        let middle = f64::midpoint(from, to);
        let Landing::At(m) = arm_point((a, b), touch, half, (middle, arms), near, tol) else {
            return None;
        };
        if distance_to_segment(m.point, p.point, q.point) <= near {
            k += 1;
            continue;
        }
        if points.len() >= budget {
            return None;
        }
        points.insert(k + 1, (m, middle));
    }
    let mut arm = Traced {
        points: Vec::with_capacity(points.len()),
        on_a: Vec::with_capacity(points.len()),
        on_b: Vec::with_capacity(points.len()),
        stopped: Stopped::LeftTheDomain,
    };
    for (c, _) in points {
        arm.points.push(c.point);
        arm.on_a.push(beside(a, c.on_a, exit.on_a));
        arm.on_b.push(beside(b, c.on_b, exit.on_b));
    }
    Some(arm)
}

/// What [`straightened`] makes of one sample.
#[derive(Debug, Clone, Copy)]
enum Placed {
    /// Left as it is.
    Kept,
    /// Moved onto the arm it lies out on.
    Moved(Contact),
    /// No point of that arm within the chord: past where the arm leaves
    /// the chart, with where it does where that is found.
    Off(Option<Contact>),
}

/// Where a sample stands about the touches, for [`straightened`]: on a
/// touch, on a half arm of one, or clear of both.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Beside {
    Touch(usize),
    Arm(usize, usize),
    Clear,
}

/// A sample beside a touch where the section crosses itself, placed on the
/// arm it lies out on: see [`straightened`].
fn placed(
    touches: &[Touch],
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    sample: Contact,
    near: f64,
    tol: Tolerances,
) -> (Placed, Beside) {
    let p = sample.point;
    let Some((t, touch)) = touches
        .iter()
        .enumerate()
        .filter(|(_, t)| t.arms.is_some())
        .min_by(|x, y| x.1.point.distance(p).total_cmp(&y.1.point.distance(p)))
    else {
        return (Placed::Kept, Beside::Clear);
    };
    let r = p.distance(touch.point);
    if r <= near {
        return (Placed::Kept, Beside::Touch(t));
    }
    let (Some(arms), Some(h)) = (touch.arms, half_arm(touch, p)) else {
        return (Placed::Kept, Beside::Clear);
    };
    let half = halves(arms)[h];
    let along = (p - touch.point).dot(half);
    if r > touch.around || along <= 0.0 || touch.spread <= 0.0 {
        return (Placed::Kept, Beside::Clear);
    }
    let band = tol.confusion() / (along * touch.spread);
    let inside = r <= touch.snap;
    if !inside && band <= near {
        return (Placed::Kept, Beside::Clear);
    }
    let on_border = outside(a, sample.on_a, tol) || outside(b, sample.on_b, tol);
    let placed = match arm_point((a, b), touch, half, (along, arms), near, tol) {
        Landing::At(c) if on_border => {
            if c.point.distance(p) <= near || !inside {
                Placed::Kept
            } else {
                Placed::Off(None)
            }
        }
        Landing::At(c) if inside || c.point.distance(p) <= band => Placed::Moved(c),
        Landing::Leaves(..) if inside => {
            match arm_exit((a, b), touch, half, (along, arms), near, tol) {
                Some(exit) if exit.point.distance(p) <= near => Placed::Moved(exit),
                exit => Placed::Off(exit),
            }
        }
        _ => return (Placed::Kept, Beside::Clear),
    };
    (placed, Beside::Arm(t, h))
}

/// An arc with its samples beside the touches where the section crosses
/// itself placed on the section, in the pieces left where some are not.
///
/// The march settles a sample within the chord of the section, but beside a
/// touch the solve stops in the band the confusion distance leaves round
/// an arm, `ε / (ρ · spread)` across at `ρ` out, and the samples zigzag
/// across it. Each sample between the touch and its `around`, where that
/// band is wider than the chord, is solved again on the arm it lies out
/// on, from the arm's straight line ([`arm_point`]), and moved there where
/// it lies within the band of it.
///
/// Within the touch's snap every point stands within `ε` of both surfaces
/// and the band holds the whole of it; the section there is the two arms,
/// so a sample there is moved onto the arm it lies out on. One out past
/// where that arm leaves the chart, further than the chord from where it
/// does, is a point of the band: next to the touch, or to a sample of that
/// same arm, it is the arc's end on the chart's border, misplaced along it,
/// and is moved to where the arm leaves; otherwise it lies on no arm the
/// arc runs along, and goes, the arc cut there. A piece left lying wholly
/// on one touch is that touch.
pub(super) fn straightened(
    arc: Traced,
    touches: &[Touch],
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    near: f64,
    tol: Tolerances,
) -> Vec<Traced> {
    let m = arc.points.len();
    let mut placements: Vec<(Placed, Beside)> = (0..m)
        .map(|k| {
            let sample = Contact {
                on_a: arc.on_a[k],
                on_b: arc.on_b[k],
                point: arc.points[k],
            };
            placed(touches, (a, b), sample, near, tol)
        })
        .collect();
    for k in 0..m {
        let (Placed::Off(Some(exit)), Beside::Arm(t, h)) = placements[k] else {
            continue;
        };
        let continues = [k.checked_sub(1), (k + 1 < m).then_some(k + 1)]
            .into_iter()
            .flatten()
            .any(|j| match placements[j] {
                (Placed::Off(_), _) => false,
                (_, Beside::Touch(u)) => u == t,
                (_, side) => side == Beside::Arm(t, h),
            });
        if continues {
            placements[k].0 = Placed::Moved(exit);
        }
    }
    let empty = || Traced {
        points: Vec::new(),
        on_a: Vec::new(),
        on_b: Vec::new(),
        stopped: arc.stopped,
    };
    let mut pieces = Vec::new();
    let mut piece = empty();
    for (k, (placed, _)) in placements.into_iter().enumerate() {
        let (p, on_a, on_b) = match placed {
            Placed::Kept => (arc.points[k], arc.on_a[k], arc.on_b[k]),
            Placed::Moved(c) => (
                c.point,
                beside(a, c.on_a, arc.on_a[k]),
                beside(b, c.on_b, arc.on_b[k]),
            ),
            Placed::Off(_) => {
                pieces.push(core::mem::replace(&mut piece, empty()));
                continue;
            }
        };
        piece.points.push(p);
        piece.on_a.push(on_a);
        piece.on_b.push(on_b);
    }
    pieces.push(piece);
    if pieces.len() > 1 {
        for cut in &mut pieces {
            cut.stopped = Stopped::Stalled;
        }
    }
    let on_one_touch = |piece: &Traced| {
        touches
            .iter()
            .any(|t| piece.points.iter().all(|p| p.distance(t.point) <= near))
    };
    pieces.retain(|piece| piece.points.len() >= 2 && !on_one_touch(piece));
    pieces
}

/// The touches a branch's low-angle samples are solved onto, added to
/// `touches` unless one is already there, and the touch solved from each of
/// its ends where it stalled, head then tail.
pub(super) fn find_touches(
    branch: &Traced,
    touches: &mut Vec<Touch>,
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    options: Marching,
    tol: Tolerances,
) -> [Option<usize>; 2] {
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
    let mut solved = [None; 2];
    let n = branch.points.len();
    if n < 5 {
        return solved;
    }
    // A branch with no transversal stretch to certify it may be tangential
    // contact rather than a curve, and only crossings are taken from it.
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
                .any(|t| t.point.distance(branch.points[i]) <= near)
        {
            continue;
        }
        let Some(mut touch) = touching_point(a, b, branch.on_a[i], branch.on_b[i], tol) else {
            continue;
        };
        // A touch whose section passes it without crossing there is none.
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
        // Two solves land on one touch to within the chord the march
        // resolves, or where either may land; touches further apart are
        // distinct, however close.
        let at = touches
            .iter()
            .position(|t| same_touch(t, &touch, near))
            .unwrap_or_else(|| {
                touches.push(touch);
                touches.len() - 1
            });
        if stall {
            solved[usize::from(i != 0)] = Some(at);
        }
    }
    solved
}
