//! The section solved between samples, and arcs carried onto touches.

use super::*;

/// Samples of the section between two of its points, each chord between
/// neighbours standing within `near` of it: the section solved on the plane
/// square to a chord through its middle ([`middle_of`]), the chord halved
/// where that lands further off. `None` where a solve does not settle
/// within half the chord of its middle, or more samples than `budget` would
/// be needed.
pub(super) fn fill(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    p: Contact,
    q: Contact,
    near: f64,
    budget: &mut usize,
    tol: Tolerances,
) -> Option<Vec<Contact>> {
    if p.point.distance(q.point) <= near {
        return Some(Vec::new());
    }
    let Middle {
        landed: m,
        off,
        on_both,
    } = middle_of((a, b), p, q, near, tol)?;
    // A chord whose middle lies within the confusion distance of both
    // surfaces is on the section to tolerance, though the section solved
    // through it lands further off: an end settled onto a touch stands off
    // the exact crossing by up to the touch's own tolerance.
    if off <= near || on_both {
        return Some(Vec::new());
    }
    *budget = budget.checked_sub(1)?;
    let mut out = fill((a, b), p, m, near, budget, tol)?;
    out.push(m);
    out.extend(fill((a, b), m, q, near, budget, tol)?);
    Some(out)
}

/// The section beside the middle of a chord between two of its points:
/// see [`middle_of`].
struct Middle {
    /// The section there, its parameters beside the chord's start.
    landed: Contact,
    /// How far it stands off the chord's middle.
    off: f64,
    /// Whether the chord's middle itself lies within the confusion
    /// distance of both surfaces.
    on_both: bool,
}

/// The section solved on the plane square to the chord between two of its
/// points through its middle, to within `near` of it. `None` where the
/// solve does not settle within half the chord of the middle.
fn middle_of(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    p: Contact,
    q: Contact,
    near: f64,
    tol: Tolerances,
) -> Option<Middle> {
    let chord = q.point - p.point;
    let length = chord.magnitude();
    let middle = p.point + chord * 0.5;
    // The solve starts from the points of both surfaces nearest the middle
    // of the chord, so it lands on the stretch of section beside it: from
    // an end standing on a touch, where the surfaces are tangent and the
    // solve is singular, it may land on another stretch crossing there.
    let feet = nearest_on(a, p.on_a, middle, tol).zip(nearest_on(b, p.on_b, middle, tol));
    let start = feet.map_or(
        [p.on_a.0, p.on_a.1, p.on_b.0, p.on_b.1],
        |((ua, _), (ub, _))| [ua.0, ua.1, ub.0, ub.1],
    );
    // Where the solve gives out, as it does where the surfaces part too
    // slowly to tell apart, the middle of the chord lying within the
    // confusion distance of both surfaces is on the section as far as
    // anything can be.
    let on_both = || -> Option<Contact> {
        let ((on_a, pa), (on_b, pb)) = feet?;
        (pa.distance(middle) <= tol.confusion() && pb.distance(middle) <= tol.confusion())
            .then_some(Contact {
                on_a,
                on_b,
                point: middle,
            })
    };
    let landed = correct(
        a,
        b,
        start,
        middle,
        Some((p.point, chord / length, length * 0.5)),
        near,
        tol,
    )
    .filter(|c| c.point.distance(middle) <= length * 0.5)
    .or_else(on_both)?;
    let off = landed.point.distance(middle);
    if off > length * 0.5 {
        return None;
    }
    Some(Middle {
        landed: Contact {
            on_a: beside(a, landed.on_a, p.on_a),
            on_b: beside(b, landed.on_b, p.on_b),
            point: landed.point,
        },
        off,
        on_both: on_both().is_some(),
    })
}

/// A branch with the section solved between each two neighbouring samples,
/// on the plane square to their chord through its middle, to within `near`
/// of the section ([`middle_of`]). The trace stands within `chord` of the
/// section, so a point solved further than that off the chord it stands by
/// lies on some other stretch of the section, and the branch is not
/// filled: `None` then, or where the branch would pass `budget` samples. A
/// chord whose solve does not settle is left as it was.
pub(crate) fn with_middles(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    branch: &Traced,
    (near, chord): (f64, f64),
    budget: usize,
    tol: Tolerances,
) -> Option<Traced> {
    let n = branch.points.len();
    let contact = |i: usize| Contact {
        on_a: branch.on_a[i],
        on_b: branch.on_b[i],
        point: branch.points[i],
    };
    let mut out = Traced {
        points: Vec::with_capacity(2 * n),
        on_a: Vec::with_capacity(2 * n),
        on_b: Vec::with_capacity(2 * n),
        stopped: branch.stopped,
    };
    for i in 0..n {
        if i > 0 {
            let (p, q) = (contact(i - 1), contact(i));
            if let Some(Middle { landed: m, .. }) = middle_of((a, b), p, q, near, tol) {
                if distance_to_segment(m.point, p.point, q.point) > chord {
                    return None;
                }
                out.points.push(m.point);
                out.on_a.push(m.on_a);
                out.on_b.push(m.on_b);
            }
        }
        out.points.push(branch.points[i]);
        out.on_a.push(branch.on_a[i]);
        out.on_b.push(branch.on_b[i]);
    }
    (out.points.len() <= budget).then_some(out)
}

/// Whether a sample beside a touch runs along the arm it lies out on: the
/// section's direction there, the cross of the normals, nearer that arm's
/// line than the other's. A sample of a walk wandering across the corner
/// does not.
fn along_its_arm(
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    touch: &Touch,
    at: Contact,
    tol: Tolerances,
) -> bool {
    let (Some(arms), Some(h)) = (touch.arms, half_arm(touch, at.point)) else {
        return false;
    };
    let Some(tangent) = tangent_at(a, b, at, tol) else {
        return false;
    };
    let (own, other) = (arms[h / 2], arms[1 - h / 2]);
    tangent.dot(own).abs() >= tangent.dot(other).abs()
}

/// An arc's end carried onto a touch: samples along the straight way there
/// corrected onto both surfaces, then the touch itself. The arc stays as it
/// was where the end does not head for the touch, a sample further than
/// `reach` from it will not settle onto the section, or the samples would
/// take the arc past the points a branch may hold. Whether the end was
/// carried onto the touch.
pub(super) fn carry_onto(
    arc: &mut Traced,
    touch: &Touch,
    at_head: bool,
    (a, b): (&SurfaceGeometry, &SurfaceGeometry),
    options: Marching,
    tol: Tolerances,
) -> bool {
    let near = options.chord.max(tol.confusion());
    let reach = near * 60.0;
    let carry = reach * 40.0;
    let n = arc.points.len();
    let (end, inner) = if at_head { (0, 1) } else { (n - 1, n - 2) };
    let from = arc.points[end];
    let gap = touch.point - from;
    let distance = gap.magnitude();
    let mut added: Vec<Contact> = Vec::new();
    let end_contact = Contact {
        on_a: arc.on_a[end],
        on_b: arc.on_b[end],
        point: from,
    };
    if distance > tol.confusion()
        && distance <= touch.snap
        && along_its_arm((a, b), touch, end_contact, tol)
    {
        // The straight way to the touch lies on both surfaces to within
        // the confusion distance.
        let (on_a, on_b) = (arc.on_a[end], arc.on_b[end]);
        added.push(Contact {
            on_a: beside(a, touch.on_a, on_a),
            on_b: beside(b, touch.on_b, on_b),
            point: touch.point,
        });
    } else if distance > tol.confusion() {
        let along = gap * (1.0 / distance);
        let heading = from - arc.points[inner];
        if distance > reach && heading.dot(along) < 0.5 * heading.magnitude() {
            return false;
        }
        // How far short of the touch each sample stands: steps of at most
        // `reach` down to it, then halving to an eighth of it, where the
        // two arms still stand far enough apart to correct onto the right
        // one and the last straight stretch bows from the arm well inside
        // the chord. An end further off than `carry` is followed in to
        // `carry` short of the touch by samples placed where the chord
        // calls for them.
        let span = distance.min(carry);
        let pieces = (span / reach).ceil().max(1.0);
        #[allow(clippy::cast_precision_loss)]
        if pieces + 3.0 > options.max_points.saturating_sub(n) as f64 {
            return false;
        }
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let count = pieces as usize;
        #[allow(clippy::cast_precision_loss)]
        let mut short: Vec<f64> = (1..count)
            .map(|k| span - span * k as f64 / pieces)
            .collect();
        if distance > carry {
            short.insert(0, carry);
        }
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
            let c = Contact {
                on_a,
                on_b,
                point: c.point,
            };
            if added.is_empty() && remaining >= carry && distance > carry {
                let mut budget = options.max_points.saturating_sub(n);
                let Some(between) = fill((a, b), end_contact, c, near, &mut budget, tol) else {
                    return false;
                };
                added.extend(between);
            }
            added.push(c);
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

/// Newton projection of a point onto a surface, warm-started.
pub(crate) fn nearest_on(
    surface: &SurfaceGeometry,
    seed: (f64, f64),
    target: Point,
    tol: Tolerances,
) -> Option<((f64, f64), Point)> {
    let (mut u, mut v) = seed;
    for _ in 0..16 {
        let (u_ok, v_ok) = surface.normalize_parameters(u, v, tol).ok()?;
        u = u_ok;
        v = v_ok;
        let p = surface.point_at(u, v, tol).ok()?;
        let (su, sv) = surface.d1_at(u, v, tol).ok()?;
        let r = p - target;
        let (a11, a12, a22) = (su.dot(su), su.dot(sv), sv.dot(sv));
        let det = a11.mul_add(a22, -(a12 * a12));
        if det.abs() <= f64::MIN_POSITIVE {
            break;
        }
        let (b1, b2) = (-su.dot(r), -sv.dot(r));
        let du = b1.mul_add(a22, -(b2 * a12)) / det;
        let dv = a11.mul_add(b2, -(a12 * b1)) / det;
        u += du;
        v += dv;
        if du.hypot(dv) < 1e-14 {
            break;
        }
    }
    let (u, v) = surface.normalize_parameters(u, v, tol).ok()?;
    Some(((u, v), surface.point_at(u, v, tol).ok()?))
}
