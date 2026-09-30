//! Splitting a face by section curves, in its own parameter space.
//!
//! This is the builder's 2D half, and it deliberately knows nothing about
//! geometry beyond polylines. The caller (the pave filler) has already done
//! the exact work: every boundary edge and every section curve is split at
//! every mutual crossing by the intersectors, so what arrives here is
//! *strands*: polyline scaffolding in parameter space, each carrying a tag
//! naming the exact sub-curve it stands for, meeting other strands only at
//! endpoints. This module decides the combinatorics (which strands bound
//! which region) and hands back pieces as sequences of directed tags, so the
//! rebuild works from the exact geometry and the polylines are never promoted
//! to an answer.

use ogeom_core::predicates::{Exact, Predicates, Sign};
use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::Point2;

/// One pre-split piece of boundary or section, as scaffolding plus a name.
#[derive(Debug, Clone)]
pub(crate) struct Strand<T> {
    /// The curve's course through parameter space, finely enough sampled
    /// that its first few snaps from each end give the way it leaves.
    pub polyline: Vec<Point2>,
    /// Which exact sub-curve this stands for.
    pub tag: T,
    /// Whether this is a piece of the face's own boundary, as opposed to a
    /// section. Boundary strands define the material; a dangling boundary is
    /// an error where a dangling section is a pruning.
    pub boundary: bool,
}

/// One directed traversal of a strand inside a piece's ring.
#[derive(Debug, Clone)]
pub(crate) struct Traversal<T> {
    /// The strand's tag.
    pub tag: T,
    /// Whether the ring runs the strand backwards.
    pub reversed: bool,
}

/// One piece of a split face.
#[derive(Debug, Clone)]
pub(crate) struct Piece<T> {
    /// The boundary as directed strand traversals: ring `[0]` is the outer
    /// contour, counter-clockwise in parameter space; further rings are
    /// holes, clockwise.
    pub rings: Vec<Vec<Traversal<T>>>,
    /// The same rings as parameter-space polylines, in the same order.
    ///
    /// The caller needs the region itself, not only its name: deciding
    /// whether two coincident faces' pieces stand for the *same* patch of one
    /// surface is a containment question, and containment is asked of an
    /// outline.
    pub outlines: Vec<Vec<Point2>>,
    /// Points strictly inside the piece, best first.
    ///
    /// More than one, because a single probe can be unlucky: a piece that
    /// merely *touches* the other solid has a probe on that contact reading
    /// neither in nor out, and the way past it is to ask somewhere else in
    /// the same piece.
    pub interiors: Vec<Point2>,
}

/// Assemble pre-split strands into the pieces they bound.
///
/// Dangling sections (chains that separate no material) are pruned;
/// regions outside the boundary strands' material (a hole's inside, say) are
/// dropped by an even-odd test against the boundary polylines.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if a boundary
/// strand dangles, or the graph yields no piece at all.
pub(crate) fn assemble<T: Clone>(strands: &[Strand<T>], snap: f64) -> OgeomResult<Vec<Piece<T>>> {
    let mut live: Vec<&Strand<T>> = strands
        .iter()
        .filter(|s| s.polyline.len() >= 2 && polyline_length(&s.polyline) > snap)
        .collect();
    if live.is_empty() {
        ogeom_bail!(Construction, "a face with no boundary bounds nothing");
    }

    // Endpoints snapped to canonical nodes. Only endpoints: the pre-split
    // contract says strands meet nowhere else.
    let mut nodes: Vec<Point2> = Vec::new();
    // Nodes binned on a grid of the snap's own size, so an endpoint looks
    // only at the nine bins about it; of the nodes in reach it takes the
    // first made, the one a scan in order would have met first.
    let cell = |p: Point2| -> (i64, i64) {
        let at = |x: f64| -> i64 {
            #[allow(clippy::cast_possible_truncation)]
            let k = (x / snap.max(f64::MIN_POSITIVE)).floor() as i64;
            k
        };
        (at(p.x), at(p.y))
    };
    let mut bins: std::collections::HashMap<(i64, i64), Vec<usize>> =
        std::collections::HashMap::new();
    let mut canon = |p: Point2, nodes: &mut Vec<Point2>| -> usize {
        let (cx, cy) = cell(p);
        let mut first: Option<usize> = None;
        for dx in -1..=1_i64 {
            for dy in -1..=1_i64 {
                let key = (cx.saturating_add(dx), cy.saturating_add(dy));
                for &i in bins.get(&key).map_or(&[][..], Vec::as_slice) {
                    if nodes[i].distance(p) <= snap && first.is_none_or(|f| i < f) {
                        first = Some(i);
                    }
                }
            }
        }
        if let Some(i) = first {
            return i;
        }
        nodes.push(p);
        bins.entry((cx, cy)).or_default().push(nodes.len() - 1);
        nodes.len() - 1
    };
    let mut ends: Vec<(usize, usize)> = Vec::new();
    for strand in &live {
        let from = canon(strand.polyline[0], &mut nodes);
        let to = canon(
            *strand.polyline.last().unwrap_or(&strand.polyline[0]),
            &mut nodes,
        );
        ends.push((from, to));
    }

    // Dust welds through. A strand shorter than the snap is dropped from
    // the graph, but its two ends were still one junction of the boundary:
    // a run of consecutive dust pieces (the paving's crossing clusters on
    // a tolerant rail leave them) can jointly span more than the snap, and
    // dropping the pieces one by one would tear a gap no positional weld
    // reaches. Each dropped strand therefore aliases its endpoints, and a
    // chain of dust aliases end to end, so the surviving neighbours meet at
    // one node however long the run.
    {
        let mut parent: Vec<usize> = (0..nodes.len()).collect();
        let find = |parent: &mut Vec<usize>, mut i: usize| -> usize {
            while parent[i] != i {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            i
        };
        let mut any = false;
        for strand in strands {
            if strand.polyline.len() < 2 || polyline_length(&strand.polyline) > snap {
                continue;
            }
            let a = canon(strand.polyline[0], &mut nodes);
            let b = canon(
                *strand.polyline.last().unwrap_or(&strand.polyline[0]),
                &mut nodes,
            );
            while parent.len() < nodes.len() {
                parent.push(parent.len());
            }
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            if ra != rb {
                parent[rb] = ra;
                any = true;
            }
        }
        if any {
            for (from, to) in &mut ends {
                *from = find(&mut parent, *from);
                *to = find(&mut parent, *to);
            }
        }
    }

    // Prune dangling chains. A section that fails to separate material hangs
    // by an end; a *boundary* strand doing so means the face's own boundary
    // does not close, which no amount of pruning repairs.
    //
    // Peeled by a queue: a strand goes when either end is left on fewer
    // than two strands, and its going may leave its neighbours so, which
    // the queue takes up in turn. What goes is what peeling round after
    // round would take, and the survivors keep their order.
    {
        let mut degree = vec![0_usize; nodes.len()];
        let mut at_node: Vec<Vec<usize>> = vec![Vec::new(); nodes.len()];
        for (k, (u, v)) in ends.iter().enumerate() {
            degree[*u] += 1;
            degree[*v] += 1;
            at_node[*u].push(k);
            at_node[*v].push(k);
        }
        let mut gone = vec![false; ends.len()];
        let mut queue: std::collections::VecDeque<usize> =
            (0..nodes.len()).filter(|n| degree[*n] < 2).collect();
        while let Some(n) = queue.pop_front() {
            for &k in &at_node[n] {
                if gone[k] {
                    continue;
                }
                let (u, v) = ends[k];
                if degree[u] >= 2 && degree[v] >= 2 {
                    continue;
                }
                let strand = live[k];
                if strand.boundary {
                    if std::env::var("OGEOM_ARRANGE_DEBUG").is_ok() {
                        eprintln!(
                            "DANGLE: boundary strand {:?} .. {:?} (deg {} / {})",
                            strand.polyline[0],
                            strand.polyline[strand.polyline.len() - 1],
                            degree[u],
                            degree[v]
                        );
                    }
                    ogeom_bail!(
                        Construction,
                        "a face boundary strand dangles; the boundary does not \
                         close in parameter space"
                    );
                }
                gone[k] = true;
                for end in [u, v] {
                    degree[end] -= 1;
                    if degree[end] < 2 {
                        queue.push_back(end);
                    }
                }
            }
        }
        let mut keep_live = Vec::with_capacity(live.len());
        let mut keep_ends = Vec::with_capacity(ends.len());
        for (k, (strand, end)) in live.iter().zip(&ends).enumerate() {
            if !gone[k] {
                keep_live.push(*strand);
                keep_ends.push(*end);
            }
        }
        live = keep_live;
        ends = keep_ends;
    }
    if live.is_empty() {
        ogeom_bail!(Construction, "the face's boundary vanished in arrangement");
    }

    // Darts: twins adjacent by construction: dart 2k runs strand k forward,
    // dart 2k + 1 backward.
    let dart_count = live.len() * 2;
    let head = |d: usize| -> usize {
        let (u, v) = ends[d / 2];
        if d.is_multiple_of(2) { v } else { u }
    };
    let tail = |d: usize| -> usize {
        let (u, v) = ends[d / 2];
        if d.is_multiple_of(2) { u } else { v }
    };
    // A dart's polyline point `i` steps from its tail, that way round.
    let point = |d: usize, i: usize| -> Point2 {
        let line = &live[d / 2].polyline;
        if d.is_multiple_of(2) {
            line[i]
        } else {
            line[line.len() - 1 - i]
        }
    };
    let mut around: Vec<Vec<usize>> = vec![Vec::new(); nodes.len()];
    for d in 0..dart_count {
        around[tail(d)].push(d);
    }
    // The order of the darts round a node is read where each first leaves
    // a small circle about it, one circle for all of them. Strands that
    // meet only at their ends leave any such circle in the order they
    // stand round the node, where their first steps need not: a section
    // leaving a boundary strand all but tangentially (a side swept from a
    // face's outline crossing a round that meets that outline at a
    // tangent) can set out a hair to the wrong side of it before turning
    // in, and read by its first step it would cross the boundary. The
    // circle is a few snaps wide, and within half of every dart's reach
    // from the node, so each dart leaves it.
    let reach = |d: usize| -> f64 {
        let at = nodes[tail(d)];
        live[d / 2]
            .polyline
            .iter()
            .map(|p| p.distance(at))
            .fold(0.0_f64, f64::max)
    };
    let mut radius = vec![snap * 4.0; nodes.len()];
    for d in 0..dart_count {
        radius[tail(d)] = radius[tail(d)].min(reach(d) * 0.5);
    }
    let leaving = |d: usize| -> Point2 {
        let at = nodes[tail(d)];
        let r = radius[tail(d)];
        let count = live[d / 2].polyline.len();
        let mut previous = point(d, 0);
        for i in 1..count {
            let p = point(d, i);
            let (dp, dq) = (previous.distance(at), p.distance(at));
            if dq >= r {
                // Where the step crosses the circle, along the step.
                let t = if dq - dp > f64::MIN_POSITIVE {
                    ((r - dp) / (dq - dp)).clamp(0.0, 1.0)
                } else {
                    1.0
                };
                return Point2::new(
                    previous.x + (p.x - previous.x) * t,
                    previous.y + (p.y - previous.y) * t,
                );
            }
            previous = p;
        }
        point(d, count - 1)
    };
    // Each dart's angle once, not once per comparison.
    let heading: Vec<f64> = (0..dart_count)
        .map(|d| angle(nodes[tail(d)], leaving(d)))
        .collect();
    for ring in &mut around {
        ring.sort_by(|&x, &y| {
            heading[x]
                .partial_cmp(&heading[y])
                .unwrap_or(core::cmp::Ordering::Equal)
        });
    }
    let mut position = vec![0_usize; dart_count];
    for ring in &around {
        for (at, &d) in ring.iter().enumerate() {
            position[d] = at;
        }
    }
    let next = |d: usize| -> usize {
        let ring = &around[head(d)];
        let at = position[d ^ 1];
        // The dart before the twin in counter-clockwise order: the face on
        // the left continues there.
        ring[(at + ring.len() - 1) % ring.len()]
    };

    let mut seen = vec![false; dart_count];
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    for start in 0..dart_count {
        if seen[start] {
            continue;
        }
        let mut cycle = Vec::new();
        let mut d = start;
        loop {
            seen[d] = true;
            cycle.push(d);
            d = next(d);
            if d == start {
                break;
            }
        }
        cycles.push(cycle);
    }

    // The polyline a cycle traces, for area, containment and interior points.
    let outline = |cycle: &[usize]| -> Vec<Point2> {
        let mut out: Vec<Point2> = Vec::new();
        for &d in cycle {
            let line = &live[d / 2].polyline;
            let mut points: Vec<Point2> = if d.is_multiple_of(2) {
                line.clone()
            } else {
                line.iter().rev().copied().collect()
            };
            points.pop();
            out.append(&mut points);
        }
        out
    };

    let mut positives: Vec<(&Vec<usize>, Vec<Point2>)> = Vec::new();
    let mut negatives: Vec<(&Vec<usize>, Vec<Point2>)> = Vec::new();
    for cycle in &cycles {
        let line = outline(cycle);
        let a = area(&line);
        if a > snap * snap {
            positives.push((cycle, line));
        } else if a < -(snap * snap) {
            negatives.push((cycle, line));
        }
    }

    // The boundary strands' polylines, for the material test.
    let material: Vec<&[Point2]> = live
        .iter()
        .filter(|s| s.boundary)
        .map(|s| s.polyline.as_slice())
        .collect();

    // The nodes each cycle passes through: a hole that shares one with a
    // positive cycle is the same component, not a hole in it. Asked of the
    // nodes, not of the polylines' nearness: a hole can pass within the weld
    // of a loose boundary without meeting it, and read by distance it was
    // taken for part of the boundary and dropped.
    let nodes_of = |cycle: &[usize]| -> std::collections::BTreeSet<usize> {
        cycle.iter().map(|&d| tail(d)).collect()
    };
    // Each cycle's nodes, area and box, found once: the nesting below asks
    // them of every hole against every positive cycle, twice over.
    let boxed = |line: &[Point2]| -> (Point2, Point2) {
        line.iter().fold(
            (
                Point2::new(f64::INFINITY, f64::INFINITY),
                Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
            ),
            |(lo, hi), p| {
                (
                    Point2::new(lo.x.min(p.x), lo.y.min(p.y)),
                    Point2::new(hi.x.max(p.x), hi.y.max(p.y)),
                )
            },
        )
    };
    let within = |b: &(Point2, Point2), p: Point2| {
        p.x >= b.0.x && p.x <= b.1.x && p.y >= b.0.y && p.y <= b.1.y
    };
    let positive_nodes: Vec<_> = positives.iter().map(|(c, _)| nodes_of(c)).collect();
    let positive_area: Vec<f64> = positives.iter().map(|(_, l)| area(l).abs()).collect();
    let positive_box: Vec<_> = positives.iter().map(|(_, l)| boxed(l)).collect();
    let negative_nodes: Vec<_> = negatives.iter().map(|(c, _)| nodes_of(c)).collect();
    let meet = |x: &std::collections::BTreeSet<usize>, y: &std::collections::BTreeSet<usize>| {
        x.iter().any(|n| y.contains(n))
    };
    // A point outside a polygon's box is outside the polygon.
    let contains = |k: usize, p: Point2| within(&positive_box[k], p) && inside(&positives[k].1, p);
    let mut pieces = Vec::new();
    for (pi, (cycle, line)) in positives.iter().enumerate() {
        let mut rings = vec![traversals(cycle, &live)];
        let mut rings_outline = vec![line.clone()];
        for (hi, (hole_cycle, hole)) in negatives.iter().enumerate() {
            // A hole belongs to the smallest positive cycle strictly
            // containing it; sharing a node means same component, not a hole.
            if meet(&negative_nodes[hi], &positive_nodes[pi]) {
                continue;
            }
            if !contains(pi, hole[0]) {
                continue;
            }
            let direct = !(0..positives.len()).any(|oi| {
                oi != pi
                    && contains(pi, positives[oi].1[0])
                    && positive_area[oi] < positive_area[pi]
                    && contains(oi, hole[0])
                    && !meet(&negative_nodes[hi], &positive_nodes[oi])
            });
            if direct {
                rings.push(traversals(hole_cycle, &live));
                rings_outline.push(hole.clone());
            }
        }
        let interiors = interior_points(&rings_outline, snap);
        let Some(interior) = interiors.first().copied() else {
            continue;
        };
        if !inside_many(&material, interior) {
            continue;
        }
        // Only probes inside the material this arrangement bounds are of
        // any use to a caller asking "where does this piece stand".
        let interiors: Vec<Point2> = interiors
            .into_iter()
            .filter(|p| inside_many(&material, *p))
            .collect();
        pieces.push(Piece {
            rings,
            outlines: rings_outline,
            interiors,
        });
    }
    if pieces.is_empty() {
        ogeom_bail!(Construction, "arrangement left no piece of the face");
    }
    Ok(pieces)
}

fn traversals<T: Clone>(cycle: &[usize], live: &[&Strand<T>]) -> Vec<Traversal<T>> {
    cycle
        .iter()
        .map(|&d| Traversal {
            tag: live[d / 2].tag.clone(),
            reversed: d % 2 == 1,
        })
        .collect()
}

fn polyline_length(line: &[Point2]) -> f64 {
    line.windows(2).map(|w| w[0].distance(w[1])).sum()
}

/// The angle of the direction from `from` towards `to`.
fn angle(from: Point2, to: Point2) -> f64 {
    let v = to - from;
    v.y.atan2(v.x)
}

/// The signed area: positive for counter-clockwise.
fn area(ring: &[Point2]) -> f64 {
    let mut doubled = 0.0;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        doubled += a.x.mul_add(b.y, -(b.x * a.y));
    }
    doubled * 0.5
}

/// Whether the ray from `p` toward `+x` crosses the segment `a`→`b`,
/// decided by the predicate rather than by a computed intersection.
///
/// The straddle test picks the segments the ray could cross; for those, the
/// crossing question is exactly "which side of the segment's line does `p`
/// lie on", which is `orient2d`'s question. A point exactly on the segment's
/// line reads as no crossing, matching the strict comparison this replaces;
/// the callers probe interior sample points, never boundary ones.
fn ray_crosses_segment<P: Predicates>(a: Point2, b: Point2, p: Point2) -> bool {
    if (a.y > p.y) == (b.y > p.y) {
        return false;
    }
    let side = P::orient2d([a.x, a.y], [b.x, b.y], [p.x, p.y]);
    if b.y > a.y {
        side == Sign::Positive
    } else {
        side == Sign::Negative
    }
}

/// As [`ray_crosses_segment`], with the ray leaning off the axes.
fn slanted_ray_crosses_segment<P: Predicates>(a: Point2, b: Point2, p: Point2) -> bool {
    // A tangent chain puts whole strands exactly along an axis-aligned
    // junction line, where a horizontal ray grazes corner after corner
    // inside rounding noise and counts them at random. No real boundary
    // runs along this slope.
    const SLANT: f64 = 0.618_033_988_749_894_9;
    let shifted = [p.x + 1.0, p.y + SLANT];
    let above = |q: Point2| P::orient2d([p.x, p.y], shifted, [q.x, q.y]) == Sign::Positive;
    if above(a) == above(b) {
        return false;
    }
    let side = P::orient2d([a.x, a.y], [b.x, b.y], [p.x, p.y]);
    if above(b) {
        side == Sign::Positive
    } else {
        side == Sign::Negative
    }
}

/// Even-odd containment with the leaning ray: the entry for probes that may
/// legitimately sit along an axis-aligned line of the boundary (a contact
/// strand down a tangent junction), where the horizontal ray is degenerate.
pub(crate) fn inside_many_slanted(lines: &[&[Point2]], p: Point2) -> bool {
    let mut inside = false;
    for line in lines {
        for w in line.windows(2) {
            if slanted_ray_crosses_segment::<Exact>(w[0], w[1], p) {
                inside = !inside;
            }
        }
    }
    inside
}

/// Even-odd containment of a point in one closed polyline.
fn inside(ring: &[Point2], p: Point2) -> bool {
    let mut inside = false;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        if ray_crosses_segment::<Exact>(a, b, p) {
            inside = !inside;
        }
    }
    inside
}

/// Even-odd containment of a point in a piece's rings.
///
/// A ring is a closed loop stored without repeating its first point, so each
/// one is closed as it is walked. Ring `[0]` is the outer contour and the
/// rest are holes, and the even-odd count over all of them answers for the
/// region they jointly bound.
///
/// This is the test for a *piece*. [`inside_many`] is the test for the
/// strands a piece is assembled from, which close only jointly and must not
/// be closed one by one.
pub(crate) fn inside_rings(rings: &[Vec<Point2>], p: Point2) -> bool {
    rings.iter().fold(false, |acc, ring| acc != inside(ring, p))
}

/// Even-odd containment against open polylines that jointly close.
///
/// The strands are pieces of closed rings, so counting crossings segment by
/// segment over all of them gives the same even-odd answer the assembled
/// rings would.
pub(crate) fn inside_many(lines: &[&[Point2]], p: Point2) -> bool {
    let mut inside = false;
    for line in lines {
        for w in line.windows(2) {
            if ray_crosses_segment::<Exact>(w[0], w[1], p) {
                inside = !inside;
            }
        }
    }
    inside
}

/// Points strictly inside the region the rings bound, by scanline.
///
/// A horizontal line through the widest gap between distinct vertex heights
/// cannot pass through a vertex or run along a horizontal segment, so its
/// crossings are transversal and the midpoint of the first inside interval
/// is interior with room to spare.
///
/// The rest are for the caller who has to ask again. They vary in *both*
/// chart directions (other heights, and other positions along each) because
/// a piece that merely touches the other solid touches it somewhere, and a
/// second opinion taken from the same place is not one.
fn interior_points(rings: &[Vec<Point2>], snap: f64) -> Vec<Point2> {
    let mut heights: Vec<f64> = rings.iter().flatten().map(|p| p.y).collect();
    heights.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    heights.dedup_by(|a, b| (*a - *b).abs() <= snap);
    // Widest gap first (that is the scanline with the most room), then the
    // rest, so a caller that needs a second opinion has one.
    let mut levels: Vec<(f64, f64)> = heights
        .windows(2)
        .map(|pair| (pair[1] - pair[0], f64::midpoint(pair[0], pair[1])))
        .collect();
    levels.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(core::cmp::Ordering::Equal));
    // An unsplit chart (a whole sphere, pole to pole) has exactly one gap
    // and therefore one scanline, straight across its middle. That is
    // precisely where a solid seated on its equator touches it, so the widest
    // gap also offers its quarter heights: same piece, different latitude.
    if let Some(&(gap, level)) = levels.first() {
        levels.insert(1, (gap, gap.mul_add(-0.25, level)));
        levels.insert(2, (gap, gap.mul_add(0.25, level)));
    }

    // Every inside interval of every scanline is a candidate, and they are
    // ranked by width. The first interval of the roomiest scanline is not
    // good enough: a piece with a cusp (the sliver beside a line tangent
    // to a circle, which is what a ball inscribed in a cylinder leaves on
    // any plane through both) puts that interval inside the cusp, where
    // the "interior" point is within rounding of the boundary and reads as
    // lying on it.
    let mut candidates: Vec<(f64, Point2)> = Vec::new();
    for (gap, level) in levels {
        if gap <= snap {
            continue;
        }
        let mut crossings: Vec<f64> = Vec::new();
        for ring in rings {
            for i in 0..ring.len() {
                let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
                if (a.y > level) != (b.y > level) {
                    crossings.push((b.x - a.x).mul_add((level - a.y) / (b.y - a.y), a.x));
                }
            }
        }
        crossings.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
        for pair in crossings.as_chunks::<2>().0 {
            let width = pair[1] - pair[0];
            if width > snap {
                // Room is the lesser of the interval's width and its
                // scanline's gap: a level in a thin gap may pass just beyond
                // an arc's chords and inside the arc itself, where the
                // chords report an interval as wide as the whole piece.
                let room = width.min(gap);
                candidates.push((room, Point2::new(f64::midpoint(pair[0], pair[1]), level)));
                // And its quarter positions, ranked below the midpoint. Moving
                // the scanline is not enough on its own: a piece symmetric
                // about a chart-vertical line (a cylinder band, a revolved
                // wall, a chart rectangle) has the same midpoint at every
                // height, so a solid touching it along that line is met by
                // every one of these "different" probes at once. The quarter
                // heights above exist for the same reason in the other
                // direction; this is that rule, applied to the width.
                candidates.push((room * 0.5, Point2::new(width.mul_add(0.25, pair[0]), level)));
                candidates.push((room * 0.5, Point2::new(width.mul_add(0.75, pair[0]), level)));
            }
        }
    }
    candidates.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(core::cmp::Ordering::Equal));

    // Room first, but never room alone. Ranked purely by width, every scanline
    // of a piece symmetric about a chart-vertical line offers its own midpoint
    // before any of them offers a different position, so a caller asking nine
    // times asks the same column nine times. Take the roomiest candidate at
    // each distinct column first, then the rest in width order, which leaves
    // the other columns at the front for a touch running down one, and the
    // other heights right behind them for a touch running across.
    //
    // Nine are returned, so nine columns end the search, and nine of the
    // rest are all that can follow them: a face with hundreds of holes
    // offers hundreds of thousands of candidates, and comparing each with
    // every column already chosen cost seconds for nine points.
    const PROBES: usize = 9;
    let mut chosen: Vec<(f64, Point2)> = Vec::new();
    let mut rest: Vec<(f64, Point2)> = Vec::new();
    for candidate in candidates {
        if chosen.len() == PROBES {
            break;
        }
        if chosen
            .iter()
            .any(|(_, p)| (p.x - candidate.1.x).abs() <= snap)
        {
            if rest.len() < PROBES {
                rest.push(candidate);
            }
        } else {
            chosen.push(candidate);
        }
    }
    chosen.extend(rest);
    chosen.truncate(PROBES);
    chosen.into_iter().map(|(_, p)| p).collect()
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    /// The even-odd containment is a sign decision, and it goes through the
    /// predicate: the crossing question near a slanted edge is answered by
    /// `orient2d` on both sides of the edge, at separations where a computed
    /// intersection abscissa would be pure rounding.
    #[test]
    fn containment_is_decided_by_the_predicate() {
        let ring = vec![
            Point2::new(0.0, -1.0),
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
            Point2::new(-1.0, 0.0),
        ];
        assert!(inside(&ring, Point2::new(0.0, 0.0)));
        assert!(!inside(&ring, Point2::new(0.75, 0.5)));
        // Astride the slanted edge x + y = 1: the side flips with the sign.
        let eps = 1e-14;
        assert!(inside(&ring, Point2::new(0.5 - eps, 0.5)));
        assert!(!inside(&ring, Point2::new(0.5 + eps, 0.5)));
        // A point exactly on the edge's line, straddle satisfied, reads as
        // no crossing from either side: the convention the strict
        // comparison had, now stated by `Sign::Zero`.
        assert!(!ray_crosses_segment::<Exact>(
            Point2::new(1.0, 0.0),
            Point2::new(0.0, 1.0),
            Point2::new(0.5, 0.5),
        ));
    }

    /// A square's four sides as boundary strands, pre-split at the corners.
    fn square(side: f64) -> Vec<Strand<usize>> {
        let c = [
            Point2::new(0.0, 0.0),
            Point2::new(side, 0.0),
            Point2::new(side, side),
            Point2::new(0.0, side),
        ];
        (0..4)
            .map(|i| Strand {
                polyline: vec![c[i], c[(i + 1) % 4]],
                tag: i,
                boundary: true,
            })
            .collect()
    }

    fn strand(from: Point2, to: Point2, tag: usize) -> Strand<usize> {
        Strand {
            polyline: vec![from, to],
            tag,
            boundary: false,
        }
    }

    #[test]
    fn an_uncut_face_is_one_piece_with_its_tags_in_order() {
        let pieces = assemble(&square(2.0), 1e-7).unwrap();
        assert_eq!(pieces.len(), 1);
        assert_eq!(pieces[0].rings.len(), 1);
        assert_eq!(pieces[0].rings[0].len(), 4);
    }

    #[test]
    fn a_rectangle_offers_probes_away_from_its_own_centre_line() {
        // A piece symmetric about a chart-vertical line has the same midpoint
        // at every height, so ranking probes by room alone offers one column
        // over and over. A solid touching this piece down that column would
        // meet every probe at once, and the classifier (which asks again
        // precisely because the first answer was "on the boundary") would
        // have nowhere else to ask.
        let pieces = assemble(&square(2.0), 1e-7).unwrap();
        assert_eq!(pieces.len(), 1);

        let columns = pieces[0].interiors.iter().map(|p| p.x);
        let spread = columns.clone().fold(f64::NEG_INFINITY, f64::max)
            - columns.fold(f64::INFINITY, f64::min);
        assert!(
            spread > 1e-6,
            "every probe stands in the same column: {:?}",
            pieces[0].interiors
        );
    }

    /// A section setting out all but along the boundary, its first step a
    /// hair to the outside before it turns in (a section starting at a
    /// tangent, sampled): read where the strands leave a circle about their
    /// node, it stands inside, and the face splits in two.
    #[test]
    fn a_section_leaving_along_the_boundary_still_splits_the_face() {
        let p = Point2::new;
        let boundary = |a: Point2, b: Point2, tag: usize| Strand {
            polyline: vec![a, b],
            tag,
            boundary: true,
        };
        let strands = vec![
            boundary(p(0.0, 0.0), p(0.5, 0.0), 0),
            boundary(p(0.5, 0.0), p(1.0, 0.0), 1),
            boundary(p(1.0, 0.0), p(1.0, 1.0), 2),
            boundary(p(1.0, 1.0), p(0.5, 1.0), 3),
            boundary(p(0.5, 1.0), p(0.0, 1.0), 4),
            boundary(p(0.0, 1.0), p(0.0, 0.0), 5),
            Strand {
                polyline: vec![
                    p(0.5, 0.0),
                    p(0.501, -1e-7),
                    p(0.51, 0.001),
                    p(0.55, 0.05),
                    p(0.6, 0.5),
                    p(0.5, 1.0),
                ],
                tag: 6,
                boundary: false,
            },
        ];
        let pieces = assemble(&strands, 1e-3).unwrap();
        assert_eq!(pieces.len(), 2);
    }

    #[test]
    fn a_chord_split_at_the_boundary_makes_two_pieces() {
        // The caller's pre-split contract: the chord arrives as one strand
        // spanning wall to wall, and the walls arrive split at its feet.
        let mut strands = vec![
            Strand {
                polyline: vec![Point2::new(0.0, 0.0), Point2::new(2.0, 0.0)],
                tag: 0,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(2.0, 0.0), Point2::new(2.0, 1.0)],
                tag: 1,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(2.0, 1.0), Point2::new(2.0, 2.0)],
                tag: 2,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(2.0, 2.0), Point2::new(0.0, 2.0)],
                tag: 3,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(0.0, 2.0), Point2::new(0.0, 1.0)],
                tag: 4,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(0.0, 1.0), Point2::new(0.0, 0.0)],
                tag: 5,
                boundary: true,
            },
        ];
        strands.push(strand(Point2::new(0.0, 1.0), Point2::new(2.0, 1.0), 6));
        let pieces = assemble(&strands, 1e-7).unwrap();
        assert_eq!(pieces.len(), 2);
        // The chord appears in both pieces, once each way round.
        let uses: Vec<bool> = pieces
            .iter()
            .flat_map(|p| p.rings[0].iter())
            .filter(|t| t.tag == 6)
            .map(|t| t.reversed)
            .collect();
        assert_eq!(uses.len(), 2);
        assert_ne!(uses[0], uses[1]);
    }

    #[test]
    fn a_dangling_section_is_pruned_not_walked() {
        let mut strands = square(2.0);
        strands.push(strand(Point2::new(1.0, 0.0), Point2::new(1.0, 1.0), 9));
        let pieces = assemble(&strands, 1e-7).unwrap();
        assert_eq!(pieces.len(), 1);
        assert!(pieces[0].rings[0].iter().all(|t| t.tag != 9));
    }

    #[test]
    fn a_dangling_boundary_is_an_error() {
        let mut strands = square(2.0);
        strands.pop();
        assert!(assemble(&strands, 1e-7).is_err());
    }

    #[test]
    fn a_closed_section_loop_gives_a_hole_with_its_tags() {
        // The loop arrives as two arcs (the closed-curve pre-split), wholly
        // inside the face: the inner region is a piece, and the rest is a
        // piece with a hole whose ring carries the arcs' tags.
        let mut strands = square(4.0);
        let top = Point2::new(2.0, 3.0);
        let bottom = Point2::new(2.0, 1.0);
        strands.push(Strand {
            polyline: vec![bottom, Point2::new(3.0, 2.0), top],
            tag: 10,
            boundary: false,
        });
        strands.push(Strand {
            polyline: vec![top, Point2::new(1.0, 2.0), bottom],
            tag: 11,
            boundary: false,
        });
        let pieces = assemble(&strands, 1e-7).unwrap();
        assert_eq!(pieces.len(), 2);
        let with_hole = pieces.iter().find(|p| p.rings.len() == 2).unwrap();
        let tags: Vec<usize> = with_hole.rings[1].iter().map(|t| t.tag).collect();
        assert!(tags.contains(&10) && tags.contains(&11));
        let inner = pieces
            .iter()
            .find(|p| p.rings.len() == 1 && p.rings[0].len() == 2);
        assert!(inner.is_some(), "the loop's inside is its own piece");
    }

    #[test]
    fn curved_strands_walk_like_straight_ones() {
        // A wavy section spanning the face: the angular sort works from the
        // first polyline step, so curvature is invisible to the walk.
        let strands = vec![
            Strand {
                polyline: vec![Point2::new(0.0, 0.0), Point2::new(4.0, 0.0)],
                tag: 0,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(4.0, 0.0), Point2::new(4.0, 2.0)],
                tag: 1,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(4.0, 2.0), Point2::new(4.0, 4.0)],
                tag: 2,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(4.0, 4.0), Point2::new(0.0, 4.0)],
                tag: 3,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(0.0, 4.0), Point2::new(0.0, 2.0)],
                tag: 4,
                boundary: true,
            },
            Strand {
                polyline: vec![Point2::new(0.0, 2.0), Point2::new(0.0, 0.0)],
                tag: 5,
                boundary: true,
            },
            Strand {
                polyline: vec![
                    Point2::new(0.0, 2.0),
                    Point2::new(1.0, 2.4),
                    Point2::new(2.0, 2.0),
                    Point2::new(3.0, 1.6),
                    Point2::new(4.0, 2.0),
                ],
                tag: 6,
                boundary: false,
            },
        ];
        let pieces = assemble(&strands, 1e-7).unwrap();
        assert_eq!(pieces.len(), 2);
    }
}
