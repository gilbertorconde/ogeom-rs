//! Chart poles, branches cut at them, and the measures the branch tests share.

use super::*;

/// The poles of a surface's chart: each bound of a direction that does not
/// repeat along which the surface collapses to a point (a sphere's pole, a
/// patch's collapsed side), as the unknown that holds it (0 for `u`, 1 for
/// `v`), the bound, and the point.
#[must_use]
pub fn chart_poles(surface: &SurfaceGeometry, tol: Tolerances) -> Vec<(usize, f64, Point)> {
    const ALONG: u32 = 8;
    let ((ua, ub), (va, vb)) = surface.domain();
    let mut out = Vec::new();
    for (k, (lo, hi), periodic, (flo, fhi)) in [
        (0_usize, (ua, ub), surface.is_periodic_u(), (va, vb)),
        (1, (va, vb), surface.is_periodic_v(), (ua, ub)),
    ] {
        if periodic || !(hi - lo).is_finite() || !(fhi - flo).is_finite() || fhi <= flo {
            continue;
        }
        for bound in [lo, hi] {
            let at = |f: f64| {
                if k == 0 {
                    surface.point_at(bound, f, tol)
                } else {
                    surface.point_at(f, bound, tol)
                }
            };
            let Ok(first) = at(flo) else {
                continue;
            };
            let collapsed = (1..=ALONG).all(|j| {
                at((fhi - flo).mul_add(f64::from(j) / f64::from(ALONG), flo))
                    .is_ok_and(|p| p.distance(first) <= tol.confusion())
            });
            if collapsed {
                out.push((k, bound, first));
            }
        }
    }
    out
}

/// A surface's state at one of its chart's poles, on the curve leaving the
/// pole along `leaving`: the pole's bound, and the free parameter whose
/// line leaves the pole that way.
///
/// At the pole every value of the free parameter names the one point, and
/// a curve through it is read at one value arriving and another leaving
/// (half a turn round a sphere's). The value a curve's image ends on is
/// the limit along it: the line of the free parameter whose derivative
/// across the bound, turned into the domain, points along the curve. It is
/// scanned over the free parameter's span and refined by golden section.
fn pole_parameter(
    surface: &SurfaceGeometry,
    (k, bound): (usize, f64),
    leaving: Vector,
    tol: Tolerances,
) -> Option<(f64, f64)> {
    const SCAN: u32 = 72;
    let ((ua, ub), (va, vb)) = surface.domain();
    let ((lo, hi), (flo, fhi)) = if k == 0 {
        ((ua, ub), (va, vb))
    } else {
        ((va, vb), (ua, ub))
    };
    let inward = if (bound - hi).abs() <= (bound - lo).abs() {
        -1.0
    } else {
        1.0
    };
    let place = |f: f64| if k == 0 { (bound, f) } else { (f, bound) };
    let alignment = |f: f64| -> f64 {
        let (u, v) = place(f);
        surface
            .d1_at(u, v, tol)
            .map_or(f64::NEG_INFINITY, |(du, dv)| {
                let across = if k == 0 { du } else { dv } * inward;
                let length = across.magnitude();
                if length <= f64::MIN_POSITIVE {
                    f64::NEG_INFINITY
                } else {
                    across.dot(leaving) / length
                }
            })
    };
    let at = |j: u32| (fhi - flo).mul_add(f64::from(j) / f64::from(SCAN), flo);
    let best = (0..=SCAN).max_by(|&x, &y| alignment(at(x)).total_cmp(&alignment(at(y))))?;
    if !alignment(at(best)).is_finite() {
        return None;
    }
    let (mut left, mut right) = (at(best.saturating_sub(1)), at((best + 1).min(SCAN)));
    let golden = (5.0_f64.sqrt() - 1.0) / 2.0;
    while right - left > tol.parametric() * 1e-3 {
        let (x, y) = (
            right - golden * (right - left),
            left + golden * (right - left),
        );
        if alignment(x) < alignment(y) {
            left = x;
        } else {
            right = y;
        }
    }
    Some(place(f64::midpoint(left, right)))
}

/// A closed branch held to its first lap: where an interior segment passes
/// over the branch's own start within `near`, heading the way the branch
/// left it, the loop was complete there and what follows goes round again.
/// A section traced twice round states a curve of twice its length whose
/// image crosses the chart's seams and poles twice, and no face is bounded
/// by it. A single lap passes its start only at its ends, or across it.
pub(super) fn first_lap(branch: Traced, near: f64) -> Traced {
    let n = branch.points.len();
    if !branch.closed() || n < 5 {
        return branch;
    }
    let start = branch.points[0];
    let leaving = branch.points[1] - start;
    for j in 1..n - 2 {
        let (p, q) = (branch.points[j], branch.points[j + 1]);
        let along = q - p;
        let length = along.square_magnitude();
        if length <= f64::MIN_POSITIVE || along.dot(leaving) <= 0.0 {
            continue;
        }
        let t = (start - p).dot(along) / length;
        if t > 0.0 && t < 1.0 && start.distance(p + along * t) <= near {
            let mut points = branch.points[..=j].to_vec();
            let mut on_a = branch.on_a[..=j].to_vec();
            let mut on_b = branch.on_b[..=j].to_vec();
            points.push(start);
            on_a.push(branch.on_a[0]);
            on_b.push(branch.on_b[0]);
            return Traced {
                points,
                on_a,
                on_b,
                stopped: branch.stopped,
            };
        }
    }
    branch
}

/// A branch cut where it runs through a pole of either surface's chart:
/// the pieces between the poles, each ending on the pole, its image there
/// on the line of the chart it arrives or leaves along.
///
/// At a pole the chart's free parameter turns (half a turn round a
/// sphere's pole) where the curve does not, so no one image follows a
/// curve through it: fitted whole, the image swings across the chart
/// between two samples, or a loop closing on the pole joins two ends a
/// half turn apart. A pole lying on the other surface, within the chord
/// of the branch's polyline, is a point of the section the branch passes
/// through; it is put among the samples there (or a sample on it taken
/// for it), and the branch cut. Each piece's end on the pole is read on
/// the line of the free parameter its curve leaves along, the section's
/// tangent there oriented into the piece.
///
/// `reach` is how far the pole may stand off the other surface and off the
/// polyline and still be a point of the section: the march's chord as it
/// traces, or a fitted section's stated tolerance where the fit passes the
/// pole within it (see [`passes_a_pole`]).
#[must_use]
pub fn cut_at_poles(
    branch: Traced,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    reach: f64,
    tol: Tolerances,
) -> Vec<Traced> {
    let n = branch.points.len();
    if n < 3 {
        return vec![branch];
    }
    let poles: Vec<(usize, (usize, f64), Point)> = [(a, 0_usize), (b, 2)]
        .into_iter()
        .flat_map(|(surface, first)| {
            chart_poles(surface, tol)
                .into_iter()
                .map(move |(k, bound, point)| (first, (k, bound), point))
        })
        .collect();
    if poles.is_empty() {
        return vec![branch];
    }
    let closed = branch.closed();
    let reach = reach.max(tol.confusion());
    let state = |i: usize| {
        [
            branch.on_a[i].0,
            branch.on_a[i].1,
            branch.on_b[i].0,
            branch.on_b[i].1,
        ]
    };
    let segments = n - 1;
    let gap = |p: Point, j: usize| distance_to_segment(p, branch.points[j], branch.points[j + 1]);
    // Each visit to a pole: the segment it lies beside, and the pole.
    let mut visits: Vec<(usize, usize)> = Vec::new();
    for (index, &(first, _, point)) in poles.iter().enumerate() {
        // The pole on the other surface, read from the nearest sample.
        let nearest = (0..n)
            .min_by(|&x, &y| {
                branch.points[x]
                    .distance(point)
                    .total_cmp(&branch.points[y].distance(point))
            })
            .unwrap_or(0);
        let other = if first == 0 { b } else { a };
        let guess = if first == 0 {
            branch.on_b[nearest]
        } else {
            branch.on_a[nearest]
        };
        if ogeom_geom::project_on_surface_from(other, point, guess, tol)
            .map_or(true, |foot| foot.distance > reach)
        {
            continue;
        }
        for j in 0..segments {
            let here = gap(point, j);
            let (before, after) = if closed {
                (
                    gap(point, (j + segments - 1) % segments),
                    gap(point, (j + 1) % segments),
                )
            } else {
                (
                    if j > 0 {
                        gap(point, j - 1)
                    } else {
                        f64::INFINITY
                    },
                    if j + 1 < segments {
                        gap(point, j + 1)
                    } else {
                        f64::INFINITY
                    },
                )
            };
            if here <= reach && here <= before && here < after {
                visits.push((j, index));
            }
        }
    }
    if visits.is_empty() {
        return vec![branch];
    }
    // The state on each surface at a pole visited, the curve leaving it
    // along `leaving`: on a surface whose pole it is, the pole's line that
    // way; on the other, the foot of the pole.
    let pole_state = |point: Point, guess: [f64; 4], leaving: Vector| -> Option<[f64; 4]> {
        let mut out = [0.0; 4];
        for (surface, first) in [(a, 0_usize), (b, 2)] {
            let own = poles
                .iter()
                .find(|(f, _, p)| *f == first && p.distance(point) <= tol.confusion());
            let (u, v) = match own {
                Some(&(_, bound, _)) => pole_parameter(surface, bound, leaving, tol)?,
                None => {
                    ogeom_geom::project_on_surface_from(
                        surface,
                        point,
                        (guess[first], guess[first + 1]),
                        tol,
                    )
                    .ok()?
                    .parameters
                }
            };
            out[first] = u;
            out[first + 1] = v;
        }
        Some(out)
    };
    // The section's direction at a pole: the line both tangent planes hold.
    // Across a collapsed bound the derivative sweeps the tangent plane, so
    // two of its lines a quarter of the span apart give the pole's normal.
    let normal_at_pole = |surface: &SurfaceGeometry, (k, bound): (usize, f64)| -> Option<Vector> {
        let ((ua, ub), (va, vb)) = surface.domain();
        let (flo, fhi) = if k == 0 { (va, vb) } else { (ua, ub) };
        let across = |f: f64| {
            let (u, v) = if k == 0 { (bound, f) } else { (f, bound) };
            surface
                .d1_at(u, v, tol)
                .ok()
                .map(|(du, dv)| if k == 0 { du } else { dv })
        };
        let quarter = (fhi - flo) / 4.0;
        let cross = across(flo)?.cross(across(flo + quarter)?);
        let length = cross.magnitude();
        (length > f64::MIN_POSITIVE).then(|| cross * (1.0 / length))
    };
    let direction_at = |point: Point, guess: [f64; 4]| -> Option<Vector> {
        let mut normals = [Vector::ZERO; 2];
        for (slot, (surface, first)) in [(a, 0_usize), (b, 2)].into_iter().enumerate() {
            let own = poles
                .iter()
                .find(|(f, _, p)| *f == first && p.distance(point) <= tol.confusion());
            normals[slot] = match own {
                Some(&(_, bound, _)) => normal_at_pole(surface, bound)?,
                None => {
                    let foot = ogeom_geom::project_on_surface_from(
                        surface,
                        point,
                        (guess[first], guess[first + 1]),
                        tol,
                    )
                    .ok()?;
                    normal_at(surface, foot.parameters, tol)?
                }
            };
        }
        tangent_of_normals(normals[0], normals[1], tol)
    };
    // The samples with each pole visited in its place, marked: a sample on
    // the pole is taken for it, and a pole between two samples is put
    // between them, on the side of the nearer it lies.
    let mut snaps = vec![false; n];
    let mut inserted: Vec<(usize, Point)> = Vec::new();
    for &(j, pole) in &visits {
        let point = poles[pole].2;
        // The sample of the segment nearer the pole, and which side of it
        // the pole lies on along the branch there.
        let i = if branch.points[j].distance(point) <= branch.points[j + 1].distance(point) {
            j
        } else {
            j + 1
        };
        if branch.points[i].distance(point) <= tol.confusion() {
            snaps[i] = true;
            continue;
        }
        let (before, after) = if closed {
            let count = n - 1;
            ((i + count - 1) % count, (i + 1) % count)
        } else {
            (i.saturating_sub(1), (i + 1).min(n - 1))
        };
        let ahead =
            (point - branch.points[i]).dot(branch.points[after] - branch.points[before]) > 0.0;
        inserted.push(match (ahead, i) {
            (true, _) => (i, point),
            (false, 0) if closed => (n - 2, point),
            (false, 0) => (usize::MAX, point),
            (false, i) => (i - 1, point),
        });
    }
    if closed && snaps[n - 1] {
        snaps[0] = true;
    }
    let mut samples: Vec<(Point, [f64; 4], bool)> = Vec::with_capacity(n + inserted.len());
    for &(_, point) in inserted.iter().filter(|(after, _)| *after == usize::MAX) {
        samples.push((point, state(0), true));
    }
    for (i, &snap) in snaps.iter().enumerate() {
        let mut at = branch.points[i];
        if snap {
            let pole = visits
                .iter()
                .map(|&(_, pole)| poles[pole].2)
                .find(|p| p.distance(at) <= tol.confusion());
            if let Some(pole) = pole {
                at = pole;
            }
        }
        samples.push((at, state(i), snap));
        for &(_, point) in inserted.iter().filter(|(after, _)| *after == i) {
            samples.push((point, state(i), true));
        }
    }
    // A loop's last sample is its first.
    if closed {
        let last = samples.len() - 1;
        samples[last] = samples[0];
    }
    let cuts: Vec<usize> = (0..samples.len()).filter(|&i| samples[i].2).collect();
    if cuts.is_empty() {
        return vec![branch];
    }
    let mut runs: Vec<Vec<(Point, [f64; 4], bool)>> = Vec::new();
    if closed {
        let first = cuts[0];
        let mut round: Vec<(Point, [f64; 4], bool)> = samples[first..samples.len() - 1].to_vec();
        round.extend_from_slice(&samples[..=first]);
        let mut from = 0;
        for i in 1..round.len() {
            if round[i].2 {
                if i > from {
                    runs.push(round[from..=i].to_vec());
                }
                from = i;
            }
        }
    } else {
        let mut from = 0;
        for i in 1..samples.len() {
            if samples[i].2 || i + 1 == samples.len() {
                runs.push(samples[from..=i].to_vec());
                from = i;
            }
        }
    }
    let mut pieces = Vec::with_capacity(runs.len());
    // The samples a piece's end on a pole may stand for: those its chord to
    // a sample further along passes within the trace's chord of. Beside the
    // pole a sample's image is its point's noise over its distance from the
    // pole, and the chord is all the trace claims of them.
    let beside_pole = |run: &[(Point, [f64; 4], bool)]| -> usize {
        let pole = run[0].0;
        let mut dropped = 0;
        for far in 2..run.len() {
            let chord = (1..far).all(|k| distance_to_segment(run[k].0, pole, run[far].0) <= reach);
            if !chord {
                break;
            }
            dropped = far - 1;
        }
        dropped
    };
    for mut run in runs {
        if run[0].2 {
            let dropped = beside_pole(&run);
            run.drain(1..=dropped);
        }
        if run.last().is_some_and(|s| s.2) {
            run.reverse();
            let dropped = beside_pole(&run);
            run.drain(1..=dropped);
            run.reverse();
        }
        let m = run.len();
        if m < 2 {
            continue;
        }
        for (end, inner) in [(0, 1), (m - 1, m - 2)] {
            if !run[end].2 {
                continue;
            }
            let point = run[end].0;
            let guess = run[inner].1;
            let Some(tangent) = direction_at(point, guess) else {
                return vec![branch];
            };
            let leaving = if tangent.dot(run[inner].0 - point) >= 0.0 {
                tangent
            } else {
                -tangent
            };
            let Some(at) = pole_state(point, guess, leaving) else {
                return vec![branch];
            };
            run[end].1 = at;
        }
        pieces.push(Traced {
            points: run.iter().map(|s| s.0).collect(),
            on_a: run.iter().map(|s| (s.1[0], s.1[1])).collect(),
            on_b: run.iter().map(|s| (s.1[2], s.1[3])).collect(),
            stopped: Stopped::LeftTheDomain,
        });
    }
    if !closed {
        if let Some(first) = pieces.first_mut()
            && !samples[0].2
        {
            first.stopped = branch.stopped;
        }
        if let Some(last) = pieces.last_mut()
            && !samples[samples.len() - 1].2
        {
            last.stopped = branch.stopped;
        }
    }
    pieces
}

/// Whether `curve` passes a pole of either surface's chart within `reach`
/// between its ends: nearer the pole than `reach` somewhere, its ends both
/// farther.
///
/// A section fitted whole past a pole states a tolerance about its own
/// distance from the pole, since no image follows the chart's half turn
/// there. Within that tolerance the pole is a point of the section as the
/// fit states it, and the branch is cut there as [`cut_at_poles`] cuts it
/// within the march's chord.
///
/// The nearest approach is read at its worst among evenly spaced samples,
/// refined by golden section between their neighbours.
#[must_use]
pub fn passes_a_pole(
    curve: &ogeom_geom::Curve,
    surfaces: [&SurfaceGeometry; 2],
    reach: f64,
    tol: Tolerances,
) -> bool {
    const SAMPLES: u32 = 256;
    use ogeom_geom::Curve3d;
    let (lo, hi) = curve.domain();
    let (Ok(start), Ok(end)) = (curve.point_at(lo, tol), curve.point_at(hi, tol)) else {
        return false;
    };
    let at = |k: u32| (hi - lo).mul_add(f64::from(k) / f64::from(SAMPLES), lo);
    for surface in surfaces {
        for (_, _, pole) in chart_poles(surface, tol) {
            if start.distance(pole) <= reach || end.distance(pole) <= reach {
                continue;
            }
            let off = |t: f64| {
                curve
                    .point_at(t, tol)
                    .map_or(f64::INFINITY, |p| p.distance(pole))
            };
            let nearest = (0..=SAMPLES)
                .min_by(|&x, &y| off(at(x)).total_cmp(&off(at(y))))
                .unwrap_or(0);
            let (mut l, mut h) = (
                at(nearest.saturating_sub(1)),
                at((nearest + 1).min(SAMPLES)),
            );
            let golden = (5.0_f64.sqrt() - 1.0) / 2.0;
            for _ in 0..100 {
                let (x, y) = (h - golden * (h - l), l + golden * (h - l));
                if off(x) < off(y) {
                    h = y;
                } else {
                    l = x;
                }
            }
            if off(f64::midpoint(l, h)).min(off(at(nearest))) <= reach {
                return true;
            }
        }
    }
    false
}

/// Where unknown `k` turns along three samples of a branch, to second
/// order: the states of the quadratics through them, as functions of the
/// length along their chords, `back` behind the middle one and `ahead`
/// past it, at the quadratic's vertex in `k`, held between the outer two.
pub(super) fn turning(around: [[f64; 4]; 3], (back, ahead): (f64, f64), k: usize) -> [f64; 4] {
    let [prior, here, next] = around;
    let (rise, fall) = ((here[k] - prior[k]) / back, (next[k] - here[k]) / ahead);
    let t = if fall == rise {
        0.0
    } else {
        ((rise * (ahead + back) / (rise - fall) - back) / 2.0).clamp(-back, ahead)
    };
    let weights = [
        t * (t - ahead) / (back * (back + ahead)),
        (t + back) * (ahead - t) / (back * ahead),
        (t + back) * t / (ahead * (ahead + back)),
    ];
    core::array::from_fn(|j| prior[j] * weights[0] + here[j] * weights[1] + next[j] * weights[2])
}

/// Below this sine the surfaces count as tangent at a point: the
/// branch-point certificate a stall end must carry to participate in
/// stitching.
pub(crate) const BRANCH_POINT_SINE: f64 = 0.05;

/// The sine of the normal angle at a contact: the transversality measure.
pub(super) fn crossing_sine(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    on_a: (f64, f64),
    on_b: (f64, f64),
    tol: Tolerances,
) -> f64 {
    let Ok(na) = a.normal_at(on_a.0, on_a.1, tol) else {
        return 0.0;
    };
    let Ok(nb) = b.normal_at(on_b.0, on_b.1, tol) else {
        return 0.0;
    };
    na.vector().cross(nb.vector()).magnitude()
}

/// Whether a stalled trace is a fragment rather than a curve.
///
/// Coincident or near-coincident surfaces defeat the tangency check at a seed:
/// rounding in the corrected parameters leaves the two normals a whisker apart,
/// the walk takes a couple of steps, and then stalls where the arithmetic gives
/// out. What comes back lies on both surfaces perfectly and describes nothing:
/// identical spheres yield such fragments, each a few points long.
///
/// A stalled branch shorter than a handful of chords carries no information the
/// seed did not, so it is noise from a degenerate configuration and dropped. A
/// *real* stalled branch (one that ran into a genuine tangency) has length
/// behind it and is kept, because a truncated real answer is still an answer.
///
/// The marcher is deliberately not a coincidence detector: for the pairs with
/// closed forms, [`surface_surface`](crate::surface_surface) answers
/// [`Same`](crate::Meeting::Same), and that check belongs before this one.
pub(super) fn is_fragment(branch: &Traced, options: Marching) -> bool {
    if branch.stopped != Stopped::Stalled {
        return false;
    }
    let length: f64 = branch
        .points
        .windows(2)
        .map(|pair| pair[0].distance(pair[1]))
        .sum();
    length < options.chord * 10.0
}

/// Whether a traced branch passes within a distance of a point.
///
/// Measured against the polyline's *segments*, not its vertices. The vertices
/// are a marching step apart (far more than the chord tolerance), so a seed
/// sitting neatly between two of them looks distant from both, and comparing to
/// vertices alone would report one circle many times.
pub(super) fn passes_near(branch: &Traced, p: Point, reach: f64) -> bool {
    branch
        .points
        .windows(2)
        .any(|pair| distance_to_segment(p, pair[0], pair[1]) <= reach)
}

/// Distance from a point to a segment.
pub(crate) fn distance_to_segment(p: Point, a: Point, b: Point) -> f64 {
    let along = b - a;
    let length = along.square_magnitude();
    if length <= f64::MIN_POSITIVE {
        return p.distance(a);
    }
    let t = ((p - a).dot(along) / length).clamp(0.0, 1.0);
    p.distance(a + along * t)
}

/// The farthest the march carries a stalled end onto a touch after the
/// tracing (the `carry` of [`split_at_touches_as_traced`] and its kin):
/// no join of branch ends spans more.
pub(crate) fn joining_reach(options: Marching, tol: Tolerances) -> f64 {
    options.chord.max(tol.confusion()) * 60.0 * 40.0
}
