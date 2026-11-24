//! Stalled fragments joined at branch points into the curves they belong to.

use super::*;

/// One arc between branch points, cut from a stalled fragment.
struct Arc {
    points: Vec<Point>,
    on_a: Vec<(f64, f64)>,
    on_b: Vec<(f64, f64)>,
    /// The branch-point cluster each end attaches to, if any.
    head_bp: Option<usize>,
    tail_bp: Option<usize>,
}

impl Arc {
    fn length(&self) -> f64 {
        self.points
            .windows(2)
            .map(|pair| pair[0].distance(pair[1]))
            .sum()
    }

    /// The direction the arc leaves its end, read across a window deep
    /// enough to stand clear of any residual wander.
    fn outgoing(&self, tail: bool) -> Option<Vector> {
        let n = self.points.len();
        if n < 2 {
            return None;
        }
        let window = (n - 1).min(24);
        let (at, back) = if tail {
            (n - 1, n - 1 - window)
        } else {
            (0, window)
        };
        let out = self.points[at] - self.points[back];
        let m = out.magnitude();
        (m > f64::MIN_POSITIVE).then(|| out / m)
    }
}

/// Whether a stalled branch is transversal *somewhere* in its interior:
/// the certificate that it is a curve passing branch points rather than
/// tangential-contact debris. A plane resting on a torus produces
/// fragments tangent along their whole length; a real curve through a
/// branch point is tangent only in passing.
pub(super) fn interior_is_transversal(
    branch: &Traced,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> bool {
    let n = branch.points.len();
    if n < 5 {
        return false;
    }
    [n / 4, n / 2, 3 * n / 4]
        .into_iter()
        .any(|i| crossing_sine(a, b, branch.on_a[i], branch.on_b[i], tol) > BRANCH_POINT_SINE)
}

/// Join stalled fragments that meet at branch points into the curves they
/// belong to.
///
/// Where the two normals become parallel the intersection has no single
/// direction: the walk stalls there, wanders in place while the correction
/// gives out, and a curve that passes *through* the singularity comes back
/// as fragments with parked ends. The reassembly is geometric, not
/// bookkeeping: cluster the near-tangent stall ends into branch points;
/// cut every fragment at its branch-point visits, which trims the wander
/// and separates the arcs a walk-through glued together; drop debris and
/// duplicate coverage; then, at each branch point, pair arc ends whose
/// tangents continue one another (a smooth curve crosses the singularity
/// collinearly, and the crossing curve turns through the crossing angle)
/// and chain the pairs into whole curves, closing the loops that close.
///
/// Only fragments transversal somewhere in their interior participate:
/// tangential contact along a whole curve is a different phenomenon and
/// keeps its honest fragments.
pub(super) fn stitch_stalled(
    found: Vec<Traced>,
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> Vec<Traced> {
    let reach = options.chord.max(tol.confusion()) * 60.0;
    const CONTINUES: f64 = 0.5;

    let (candidates, mut out): (Vec<Traced>, Vec<Traced>) = found.into_iter().partition(|branch| {
        branch.stopped == Stopped::Stalled && interior_is_transversal(branch, a, b, tol)
    });
    if candidates.is_empty() {
        return out;
    }

    // Branch points: the near-tangent stall ends, clustered, each the
    // first such end met, with the parameters it was walked at.
    let mut bps: Vec<Point> = Vec::new();
    let mut bp_params: Vec<((f64, f64), (f64, f64))> = Vec::new();
    for branch in &candidates {
        let n = branch.points.len();
        for at in [0, n - 1] {
            if crossing_sine(a, b, branch.on_a[at], branch.on_b[at], tol) < BRANCH_POINT_SINE {
                let p = branch.points[at];
                if !bps.iter().any(|held| held.distance(p) <= reach) {
                    bps.push(p);
                    bp_params.push((branch.on_a[at], branch.on_b[at]));
                }
            }
        }
    }
    if bps.is_empty() {
        out.extend(candidates);
        return out;
    }
    let bp_of =
        |p: Point| -> Option<usize> { bps.iter().position(|held| held.distance(p) <= reach) };

    // Cut each fragment at its branch-point visits: maximal runs of points
    // clear of every branch point become arcs, attached to the branch
    // points beside them.
    let mut arcs: Vec<Arc> = Vec::new();
    for branch in &candidates {
        let n = branch.points.len();
        let mut run_start: Option<usize> = None;
        for i in 0..=n {
            let near = i < n && bp_of(branch.points[i]).is_some();
            match (run_start, near, i == n) {
                (None, false, false) => run_start = Some(i),
                (Some(s), true, _) | (Some(s), _, true) => {
                    let e = i;
                    if e > s + 1 {
                        let head_bp = if s > 0 {
                            bp_of(branch.points[s - 1])
                        } else {
                            None
                        };
                        let tail_bp = if e < n { bp_of(branch.points[e]) } else { None };
                        arcs.push(Arc {
                            points: branch.points[s..e].to_vec(),
                            on_a: branch.on_a[s..e].to_vec(),
                            on_b: branch.on_b[s..e].to_vec(),
                            head_bp,
                            tail_bp,
                        });
                    }
                    run_start = None;
                }
                _ => {}
            }
        }
    }

    // Debris and duplicate coverage out; longest first so the fuller
    // tracing of a doubly-walked arc is the one kept.
    arcs.retain(|arc| arc.length() > options.chord * 10.0 && arc.points.len() >= 4);
    arcs.sort_by(|x, y| {
        y.length()
            .partial_cmp(&x.length())
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    let mut kept: Vec<Arc> = Vec::new();
    'candidate: for arc in arcs {
        let n = arc.points.len();
        for probe in [n / 4, n / 2, 3 * n / 4] {
            let p = arc.points[probe];
            if kept.iter().any(|held| {
                held.points
                    .windows(2)
                    .any(|pair| distance_to_segment(p, pair[0], pair[1]) <= reach)
            }) {
                continue 'candidate;
            }
        }
        kept.push(arc);
    }

    // Pair arc ends at each branch point by tangent continuation: the
    // smooth curve runs straight through, the crossing one turns.
    let ends: Vec<(usize, bool, usize, Vector)> = kept
        .iter()
        .enumerate()
        .flat_map(|(i, arc)| {
            [(false, arc.head_bp), (true, arc.tail_bp)]
                .into_iter()
                .filter_map(move |(tail, bp)| Some((i, tail, bp?, arc.outgoing(tail)?)))
        })
        .collect();
    let mut partner: Vec<Option<usize>> = vec![None; ends.len()];
    for bp in 0..bps.len() {
        loop {
            let mut best: Option<(usize, usize, f64)> = None;
            for x in 0..ends.len() {
                if partner[x].is_some() || ends[x].2 != bp {
                    continue;
                }
                for y in (x + 1)..ends.len() {
                    if partner[y].is_some() || ends[y].2 != bp {
                        continue;
                    }
                    let score = -ends[x].3.dot(ends[y].3);
                    if score > CONTINUES && best.is_none_or(|(_, _, held)| score > held) {
                        best = Some((x, y, score));
                    }
                }
            }
            let Some((x, y, _)) = best else { break };
            partner[x] = Some(y);
            partner[y] = Some(x);
        }
    }

    // Chain the arcs through the pairings into whole curves.
    let end_index = |arc: usize, tail: bool| -> Option<usize> {
        ends.iter().position(|e| e.0 == arc && e.1 == tail)
    };
    let mut used = vec![false; kept.len()];
    for start in 0..kept.len() {
        if used[start] {
            continue;
        }
        // Walk backwards first to a free entry, unless the chain loops.
        let mut first = start;
        let mut first_reversed = false;
        let mut seen_back = vec![false; kept.len()];
        loop {
            seen_back[first] = true;
            // Forward traversal enters an arc at its head, reversed at its
            // tail, so the entry end is named by the orientation flag.
            let Some(entry) = end_index(first, first_reversed) else {
                break;
            };
            let Some(p) = partner[entry] else { break };
            let (prev, prev_tail, _, _) = ends[p];
            if seen_back[prev] {
                break; // the chain is a loop; any start serves
            }
            first = prev;
            // The previous arc *leaves* through the paired end: leaving at
            // its tail means it runs forward.
            first_reversed = !prev_tail;
        }

        // Now walk forwards from `first`, consuming arcs.
        let mut points: Vec<Point> = Vec::new();
        let mut on_a: Vec<(f64, f64)> = Vec::new();
        let mut on_b: Vec<(f64, f64)> = Vec::new();
        let mut current = first;
        let mut reversed = first_reversed;
        let mut closed = false;
        loop {
            used[current] = true;
            let arc = &kept[current];
            type Run = (Vec<Point>, Vec<(f64, f64)>, Vec<(f64, f64)>);
            let (pts, pa, pb): Run = if reversed {
                (
                    arc.points.iter().rev().copied().collect(),
                    arc.on_a.iter().rev().copied().collect(),
                    arc.on_b.iter().rev().copied().collect(),
                )
            } else {
                (arc.points.clone(), arc.on_a.clone(), arc.on_b.clone())
            };
            // Insert the branch point itself at the junction, at the
            // parameters it was walked at: the arc's first sample can stand
            // up to the cluster's reach from it, and across a chart's
            // closure.
            if !points.is_empty() {
                let joint_bp = if reversed { arc.tail_bp } else { arc.head_bp };
                if let Some(bp) = joint_bp {
                    points.push(bps[bp]);
                    on_a.push(bp_params[bp].0);
                    on_b.push(bp_params[bp].1);
                }
            }
            points.extend(pts);
            on_a.extend(pa);
            on_b.extend(pb);

            let leaving = end_index(current, !reversed);
            let Some(l) = leaving else { break };
            let Some(p) = partner[l] else { break };
            let (next, next_tail, _, _) = ends[p];
            if used[next] {
                closed = next == first;
                break;
            }
            current = next;
            reversed = next_tail;
        }
        if closed && points.len() > 3 {
            let bridge = points[0];
            let ba = on_a[0];
            let bb = on_b[0];
            points.push(bridge);
            on_a.push(ba);
            on_b.push(bb);
        }
        out.push(Traced {
            points,
            on_a,
            on_b,
            stopped: if closed {
                Stopped::Closed
            } else {
                Stopped::Stalled
            },
        });
    }
    out
}
