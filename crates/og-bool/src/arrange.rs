//! Splitting a face by section curves, in its own parameter space.
//!
//! This is the builder's 2D half, and it deliberately knows nothing about
//! geometry beyond polylines. The caller — the pave filler — has already done
//! the exact work: every boundary edge and every section curve is split at
//! every mutual crossing by the intersectors, so what arrives here is
//! *strands*: polyline scaffolding in parameter space, each carrying a tag
//! naming the exact sub-curve it stands for, meeting other strands only at
//! endpoints. This module decides the combinatorics — which strands bound
//! which region — and hands back pieces as sequences of directed tags, so the
//! rebuild works from the exact geometry and the polylines are never promoted
//! to an answer.

use og_core::{OgResult, og_bail};
use og_math::Point2;

/// One pre-split piece of boundary or section, as scaffolding plus a name.
#[derive(Debug, Clone)]
pub(crate) struct Strand<T> {
    /// The curve's course through parameter space, finely enough sampled
    /// that the first segment approximates the tangent at the start.
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
    /// A point strictly inside the piece.
    pub interior: Point2,
}

/// Assemble pre-split strands into the pieces they bound.
///
/// Dangling sections — chains that separate no material — are pruned;
/// regions outside the boundary strands' material (a hole's inside, say) are
/// dropped by an even-odd test against the boundary polylines.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if a boundary
/// strand dangles, or the graph yields no piece at all.
pub(crate) fn assemble<T: Clone>(strands: &[Strand<T>], snap: f64) -> OgResult<Vec<Piece<T>>> {
    let mut live: Vec<&Strand<T>> = strands
        .iter()
        .filter(|s| s.polyline.len() >= 2 && polyline_length(&s.polyline) > snap)
        .collect();
    if live.is_empty() {
        og_bail!(Construction, "a face with no boundary bounds nothing");
    }

    // Endpoints snapped to canonical nodes. Only endpoints: the pre-split
    // contract says strands meet nowhere else.
    let mut nodes: Vec<Point2> = Vec::new();
    let canon = |p: Point2, nodes: &mut Vec<Point2>| -> usize {
        if let Some(i) = nodes.iter().position(|n| n.distance(p) <= snap) {
            return i;
        }
        nodes.push(p);
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

    // Prune dangling chains. A section that fails to separate material hangs
    // by an end; a *boundary* strand doing so means the face's own boundary
    // does not close, which no amount of pruning repairs.
    loop {
        let mut degree = vec![0_usize; nodes.len()];
        for (u, v) in &ends {
            degree[*u] += 1;
            degree[*v] += 1;
        }
        let mut dropped = false;
        let mut keep_ends = Vec::with_capacity(ends.len());
        let mut keep_live = Vec::with_capacity(live.len());
        for (strand, (u, v)) in live.iter().zip(&ends) {
            if degree[*u] < 2 || degree[*v] < 2 {
                if strand.boundary {
                    og_bail!(
                        Construction,
                        "a face boundary strand dangles; the boundary does not \
                         close in parameter space"
                    );
                }
                dropped = true;
                continue;
            }
            keep_live.push(*strand);
            keep_ends.push((*u, *v));
        }
        live = keep_live;
        ends = keep_ends;
        if !dropped {
            break;
        }
    }
    if live.is_empty() {
        og_bail!(Construction, "the face's boundary vanished in arrangement");
    }

    // Darts: twins adjacent by construction — dart 2k runs strand k forward,
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
    // The direction a dart leaves its tail node: the polyline's first step
    // that way round.
    let leaving = |d: usize| -> Point2 {
        let line = &live[d / 2].polyline;
        if d.is_multiple_of(2) {
            line[1]
        } else {
            line[line.len() - 2]
        }
    };
    let mut around: Vec<Vec<usize>> = vec![Vec::new(); nodes.len()];
    for d in 0..dart_count {
        around[tail(d)].push(d);
    }
    for (node, ring) in around.iter_mut().enumerate() {
        ring.sort_by(|&x, &y| {
            let a = angle(nodes[node], leaving(x));
            let b = angle(nodes[node], leaving(y));
            a.partial_cmp(&b).unwrap_or(core::cmp::Ordering::Equal)
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

    let mut pieces = Vec::new();
    for (cycle, line) in &positives {
        let mut rings = vec![traversals(cycle, &live)];
        let mut rings_outline = vec![line.clone()];
        for (hole_cycle, hole) in &negatives {
            // A hole belongs to the smallest positive cycle strictly
            // containing it; sharing a node means same component, not a hole.
            if hole
                .iter()
                .any(|p| line.iter().any(|q| q.distance(*p) <= snap))
            {
                continue;
            }
            if !inside(line, hole[0]) {
                continue;
            }
            let direct = !positives.iter().any(|(_, other)| {
                !core::ptr::eq(other, line)
                    && inside(line, other[0])
                    && area(other).abs() < area(line).abs()
                    && inside(other, hole[0])
                    && !hole
                        .iter()
                        .any(|p| other.iter().any(|q| q.distance(*p) <= snap))
            });
            if direct {
                rings.push(traversals(hole_cycle, &live));
                rings_outline.push(hole.clone());
            }
        }
        let Some(interior) = interior_point(&rings_outline, snap) else {
            continue;
        };
        if !inside_many(&material, interior) {
            continue;
        }
        pieces.push(Piece { rings, interior });
    }
    if pieces.is_empty() {
        og_bail!(Construction, "arrangement left no piece of the face");
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

/// The signed area — positive for counter-clockwise.
fn area(ring: &[Point2]) -> f64 {
    let mut doubled = 0.0;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        doubled += a.x.mul_add(b.y, -(b.x * a.y));
    }
    doubled * 0.5
}

/// Even-odd containment of a point in one closed polyline.
fn inside(ring: &[Point2], p: Point2) -> bool {
    let mut inside = false;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        if (a.y > p.y) != (b.y > p.y) {
            let x = (b.x - a.x).mul_add((p.y - a.y) / (b.y - a.y), a.x);
            if p.x < x {
                inside = !inside;
            }
        }
    }
    inside
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
            let (a, b) = (w[0], w[1]);
            if (a.y > p.y) != (b.y > p.y) {
                let x = (b.x - a.x).mul_add((p.y - a.y) / (b.y - a.y), a.x);
                if p.x < x {
                    inside = !inside;
                }
            }
        }
    }
    inside
}

/// A point strictly inside the region the rings bound, by scanline.
///
/// A horizontal line through the widest gap between distinct vertex heights
/// cannot pass through a vertex or run along a horizontal segment, so its
/// crossings are transversal and the midpoint of the first inside interval
/// is interior with room to spare.
fn interior_point(rings: &[Vec<Point2>], snap: f64) -> Option<Point2> {
    let mut heights: Vec<f64> = rings.iter().flatten().map(|p| p.y).collect();
    heights.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    heights.dedup_by(|a, b| (*a - *b).abs() <= snap);
    let mut best: Option<(f64, f64)> = None;
    for pair in heights.windows(2) {
        let gap = pair[1] - pair[0];
        if best.is_none_or(|(g, _)| gap > g) {
            best = Some((gap, f64::midpoint(pair[0], pair[1])));
        }
    }
    let (_, level) = best?;

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
    if crossings.len() < 2 {
        return None;
    }
    Some(Point2::new(
        f64::midpoint(crossings[0], crossings[1]),
        level,
    ))
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

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
