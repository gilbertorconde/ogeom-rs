//! Splitting a planar face by section segments, in its own parameter plane.
//!
//! This is the builder's 2D half: a face's boundary rings and the section
//! segments crossing it become a planar graph, the graph's bounded regions
//! become the face's pieces, and each piece knows a point strictly inside
//! itself for the classification that decides its fate.
//!
//! Everything here is segments — the vertical slice covers planar faces with
//! straight boundaries, so the arrangement is exact up to rounding and the
//! node snap absorbs the rounding.

use og_core::{OgResult, og_bail};
use og_math::{Point2, Vector2};

/// One piece of a split face.
#[derive(Debug, Clone)]
pub(crate) struct Piece {
    /// The boundary: ring `[0]` is the outer contour, counter-clockwise;
    /// any further rings are holes, clockwise.
    pub rings: Vec<Vec<Point2>>,
    /// A point strictly inside the piece.
    pub interior: Point2,
}

/// Split the region bounded by `boundary` along `sections`.
///
/// Sections that dangle — fail to separate material — are pruned rather than
/// carried into piece boundaries as zero-width spikes. Regions of the graph
/// outside the boundary's material (a hole's inside, say) are dropped by an
/// even-odd test against the original rings.
///
/// # Errors
///
/// [`OgError::Construction`](og_core::OgError::Construction) if a section
/// runs along a boundary segment — tangential, same-domain contact the
/// vertical slice refuses rather than resolves — or the graph produces no
/// usable region.
pub(crate) fn split(
    boundary: &[Vec<Point2>],
    sections: &[(Point2, Point2)],
    snap: f64,
) -> OgResult<Vec<Piece>> {
    let mut segments: Vec<(Point2, Point2, bool)> = Vec::new();
    for ring in boundary {
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            if a.distance(b) > snap {
                segments.push((a, b, true));
            }
        }
    }
    for (a, b) in sections {
        if a.distance(*b) > snap {
            segments.push((*a, *b, false));
        }
    }
    if segments.is_empty() {
        og_bail!(Construction, "a face with no boundary bounds nothing");
    }

    // Split every segment wherever another touches or crosses it, so the
    // graph's edges meet only at nodes.
    let mut cuts: Vec<Vec<f64>> = segments.iter().map(|_| vec![0.0, 1.0]).collect();
    for i in 0..segments.len() {
        for j in (i + 1)..segments.len() {
            let (a0, a1, a_boundary) = segments[i];
            let (b0, b1, b_boundary) = segments[j];
            match meet(a0, a1, b0, b1, snap) {
                Met::Apart => {}
                Met::At { on_a, on_b } => {
                    cuts[i].push(on_a);
                    cuts[j].push(on_b);
                }
                Met::Along => {
                    if a_boundary != b_boundary {
                        og_bail!(
                            Construction,
                            "a section runs along a face boundary; tangential \
                             same-domain contact is refused by the vertical \
                             slice rather than resolved"
                        );
                    }
                    // Two collinear overlapping segments of the same kind:
                    // boundary rings of a valid face cannot do this, and
                    // duplicate sections were removed before arrival.
                    og_bail!(Construction, "overlapping segments in one face");
                }
            }
        }
    }

    // Nodes snapped to canonical points, edges deduplicated.
    let mut nodes: Vec<Point2> = Vec::new();
    let canon = |p: Point2, nodes: &mut Vec<Point2>| -> usize {
        if let Some(i) = nodes.iter().position(|n| n.distance(p) <= snap) {
            return i;
        }
        nodes.push(p);
        nodes.len() - 1
    };
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for (i, (a, b, _)) in segments.iter().enumerate() {
        let direction = *b - *a;
        let mut ts = cuts[i].clone();
        ts.sort_by(|x, y| x.partial_cmp(y).unwrap_or(core::cmp::Ordering::Equal));
        let length = a.distance(*b);
        let mut previous: Option<f64> = None;
        let mut stops: Vec<f64> = Vec::new();
        for t in ts {
            if previous.is_none_or(|p| (t - p) * length > snap) {
                stops.push(t);
                previous = Some(t);
            }
        }
        for pair in stops.windows(2) {
            let from = canon(*a + direction * pair[0], &mut nodes);
            let to = canon(*a + direction * pair[1], &mut nodes);
            if from == to {
                continue;
            }
            let key = (from.min(to), from.max(to));
            if !edges.contains(&key) {
                edges.push(key);
            }
        }
    }

    // Dangling edges separate nothing; prune them until every node has at
    // least two edges, so no piece boundary carries a zero-width spike.
    loop {
        let mut degree = vec![0_usize; nodes.len()];
        for (u, v) in &edges {
            degree[*u] += 1;
            degree[*v] += 1;
        }
        let before = edges.len();
        edges.retain(|(u, v)| degree[*u] >= 2 && degree[*v] >= 2);
        if edges.len() == before {
            break;
        }
    }
    if edges.is_empty() {
        og_bail!(Construction, "the face's boundary vanished in arrangement");
    }

    // Darts, sorted around each node by angle; walking each dart's face with
    // the interior on the left extracts every region exactly once.
    let mut darts: Vec<(usize, usize)> = Vec::new();
    for (u, v) in &edges {
        darts.push((*u, *v));
        darts.push((*v, *u));
    }
    let mut around: Vec<Vec<usize>> = vec![Vec::new(); nodes.len()];
    for (d, (u, _)) in darts.iter().enumerate() {
        around[*u].push(d);
    }
    for (node, list) in around.iter_mut().enumerate() {
        list.sort_by(|&x, &y| {
            let a = angle(nodes[node], nodes[darts[x].1]);
            let b = angle(nodes[node], nodes[darts[y].1]);
            a.partial_cmp(&b).unwrap_or(core::cmp::Ordering::Equal)
        });
    }
    // Twins are adjacent by construction — dart `2k` and `2k + 1` — and each
    // dart's position in its node's angular ring is precomputed, so the walk
    // below never searches.
    let mut position = vec![0_usize; darts.len()];
    for ring in &around {
        for (at, &d) in ring.iter().enumerate() {
            position[d] = at;
        }
    }
    let next = |d: usize| -> usize {
        let (_, v) = darts[d];
        let twin = d ^ 1;
        let ring = &around[v];
        let at = position[twin];
        // The dart before the twin in counter-clockwise order: the face on
        // the left continues there.
        ring[(at + ring.len() - 1) % ring.len()]
    };

    let mut seen = vec![false; darts.len()];
    let mut cycles: Vec<Vec<usize>> = Vec::new();
    for start in 0..darts.len() {
        if seen[start] {
            continue;
        }
        let mut cycle = Vec::new();
        let mut d = start;
        loop {
            seen[d] = true;
            cycle.push(darts[d].0);
            d = next(d);
            if d == start {
                break;
            }
        }
        cycles.push(cycle);
    }

    // Positive-area cycles are candidate pieces; negative ones are either a
    // component's unbounded contour or a piece's hole, told apart by nesting.
    let polygon = |cycle: &[usize]| -> Vec<Point2> { cycle.iter().map(|&n| nodes[n]).collect() };
    let mut positives: Vec<Vec<Point2>> = Vec::new();
    let mut negatives: Vec<Vec<Point2>> = Vec::new();
    for cycle in &cycles {
        let ring = polygon(cycle);
        let a = area(&ring);
        if a > snap * snap {
            positives.push(ring);
        } else if a < -(snap * snap) {
            negatives.push(ring);
        }
    }

    let mut pieces = Vec::new();
    for outer in &positives {
        let mut rings = vec![outer.clone()];
        for hole in &negatives {
            // A hole belongs to the smallest positive cycle strictly
            // containing it — for this arrangement's one level of nesting per
            // component, "contained and from another component" suffices,
            // and sharing a node means same component.
            if hole
                .iter()
                .any(|p| outer.iter().any(|q| q.distance(*p) <= snap))
            {
                continue;
            }
            if inside(outer, hole[0]) {
                // Contained in this outer — but only a *direct* hole if no
                // smaller positive sits between.
                let direct = !positives.iter().any(|other| {
                    !core::ptr::eq(other, outer)
                        && inside(outer, other[0])
                        && area(other).abs() < area(outer).abs()
                        && inside(other, hole[0])
                        && !hole
                            .iter()
                            .any(|p| other.iter().any(|q| q.distance(*p) <= snap))
                });
                if direct {
                    rings.push(hole.clone());
                }
            }
        }
        let Some(interior) = interior_point(&rings, snap) else {
            continue;
        };
        // Regions outside the face's material — a source hole's inside —
        // are graph regions but not face pieces.
        if !inside_rings(boundary, interior) {
            continue;
        }
        pieces.push(Piece { rings, interior });
    }
    if pieces.is_empty() {
        og_bail!(Construction, "arrangement left no piece of the face");
    }
    Ok(pieces)
}

/// How two segments meet.
enum Met {
    Apart,
    /// At one point, with the parameter on each.
    At {
        on_a: f64,
        on_b: f64,
    },
    /// Collinear with overlap.
    Along,
}

/// Where segment `a0..a1` meets `b0..b1`, tolerantly.
fn meet(a0: Point2, a1: Point2, b0: Point2, b1: Point2, snap: f64) -> Met {
    let d = a1 - a0;
    let e = b1 - b0;
    let la = a0.distance(a1);
    let lb = b0.distance(b1);
    let denominator = d.cross(e);
    if denominator.abs() <= snap * (la + lb) {
        // Parallel. Overlapping only if collinear and the projections meet.
        let offset = (b0 - a0).cross(d) / la;
        if offset.abs() > snap {
            return Met::Apart;
        }
        let project = |p: Point2| (p - a0).dot(d) / (la * la);
        let (s0, s1) = (project(b0), project(b1));
        let (lo, hi) = (s0.min(s1), s0.max(s1));
        if hi * la < snap || lo * la > la - snap {
            return Met::Apart;
        }
        // Overlap in more than a point.
        if (hi.min(1.0) - lo.max(0.0)) * la > snap * 2.0 {
            return Met::Along;
        }
        return Met::Apart;
    }
    let t = (b0 - a0).cross(e) / denominator;
    let s = (b0 - a0).cross(d) / denominator;
    let slack_a = snap / la;
    let slack_b = snap / lb;
    if t < -slack_a || t > 1.0 + slack_a || s < -slack_b || s > 1.0 + slack_b {
        return Met::Apart;
    }
    Met::At {
        on_a: t.clamp(0.0, 1.0),
        on_b: s.clamp(0.0, 1.0),
    }
}

/// The angle of the direction from `from` to `to`.
fn angle(from: Point2, to: Point2) -> f64 {
    let v = to - from;
    v.y.atan2(v.x)
}

/// Twice the signed area — positive for counter-clockwise — halved.
fn area(ring: &[Point2]) -> f64 {
    let mut doubled = 0.0;
    for i in 0..ring.len() {
        let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
        doubled += a.x.mul_add(b.y, -(b.x * a.y));
    }
    doubled * 0.5
}

/// Even-odd containment of a point in a set of rings.
pub(crate) fn inside_rings(rings: &[Vec<Point2>], p: Point2) -> bool {
    let mut inside = false;
    for ring in rings {
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
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

fn inside(ring: &[Point2], p: Point2) -> bool {
    inside_rings(core::slice::from_ref(&ring.to_vec()), p)
}

/// A point strictly inside the region the rings bound, by scanline.
///
/// A horizontal line through the widest gap between distinct node heights
/// cannot pass through a node or run along a horizontal edge, so its
/// crossings with the rings are all transversal and the midpoint of the
/// first inside interval is interior with room to spare.
pub(crate) fn interior_point(rings: &[Vec<Point2>], snap: f64) -> Option<Point2> {
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
    // Even-odd: the span between the first pair of crossings is inside.
    Some(Point2::new(
        f64::midpoint(crossings[0], crossings[1]),
        level,
    ))
}

/// Unused import guard: `Vector2` appears only through operators.
#[allow(dead_code)]
fn keep_vector2(v: Vector2) -> Vector2 {
    v
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    fn square(side: f64) -> Vec<Point2> {
        vec![
            Point2::new(0.0, 0.0),
            Point2::new(side, 0.0),
            Point2::new(side, side),
            Point2::new(0.0, side),
        ]
    }

    #[test]
    fn an_uncut_face_is_one_piece() {
        let pieces = split(&[square(2.0)], &[], 1e-7).unwrap();
        assert_eq!(pieces.len(), 1);
        assert!((area(&pieces[0].rings[0]) - 4.0).abs() < 1e-9);
        assert!(inside_rings(&pieces[0].rings, pieces[0].interior));
    }

    #[test]
    fn one_chord_makes_two_pieces() {
        let chord = (Point2::new(0.0, 1.0), Point2::new(2.0, 1.0));
        let pieces = split(&[square(2.0)], &[chord], 1e-7).unwrap();
        assert_eq!(pieces.len(), 2);
        let mut areas: Vec<f64> = pieces.iter().map(|p| area(&p.rings[0])).collect();
        areas.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((areas[0] - 2.0).abs() < 1e-9 && (areas[1] - 2.0).abs() < 1e-9);
    }

    #[test]
    fn two_crossing_chords_make_four_pieces() {
        let chords = [
            (Point2::new(0.0, 1.0), Point2::new(2.0, 1.0)),
            (Point2::new(1.0, 0.0), Point2::new(1.0, 2.0)),
        ];
        let pieces = split(&[square(2.0)], &chords, 1e-7).unwrap();
        assert_eq!(pieces.len(), 4);
        for p in &pieces {
            assert!((area(&p.rings[0]) - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn a_closed_section_loop_makes_a_piece_and_a_pierced_piece() {
        // A square section wholly inside the face: the inner square is one
        // piece, and the rest of the face is a piece with a hole.
        let loop_sections = [
            (Point2::new(1.0, 1.0), Point2::new(3.0, 1.0)),
            (Point2::new(3.0, 1.0), Point2::new(3.0, 3.0)),
            (Point2::new(3.0, 3.0), Point2::new(1.0, 3.0)),
            (Point2::new(1.0, 3.0), Point2::new(1.0, 1.0)),
        ];
        let pieces = split(&[square(4.0)], &loop_sections, 1e-7).unwrap();
        assert_eq!(pieces.len(), 2);
        let with_hole = pieces.iter().find(|p| p.rings.len() == 2).unwrap();
        let solid = pieces.iter().find(|p| p.rings.len() == 1).unwrap();
        assert!((area(&solid.rings[0]) - 4.0).abs() < 1e-9);
        let outer = area(&with_hole.rings[0]);
        let hole = area(&with_hole.rings[1]);
        assert!((outer - 16.0).abs() < 1e-9);
        assert!((hole + 4.0).abs() < 1e-9, "holes run clockwise: {hole}");
        assert!(inside_rings(&with_hole.rings, with_hole.interior));
    }

    #[test]
    fn a_dangling_section_splits_nothing() {
        // A segment poking into the face but separating no material: pruned,
        // one piece, no spike in its boundary.
        let dangle = (Point2::new(1.0, 0.0), Point2::new(1.0, 1.0));
        let pieces = split(&[square(2.0)], &[dangle], 1e-7).unwrap();
        assert_eq!(pieces.len(), 1);
        assert!((area(&pieces[0].rings[0]) - 4.0).abs() < 1e-9);
    }

    #[test]
    fn a_corner_crossing_chord_makes_a_triangle_and_the_rest() {
        let chord = (Point2::new(1.0, 0.0), Point2::new(0.0, 1.0));
        let pieces = split(&[square(2.0)], &[chord], 1e-7).unwrap();
        assert_eq!(pieces.len(), 2);
        let mut areas: Vec<f64> = pieces.iter().map(|p| area(&p.rings[0])).collect();
        areas.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((areas[0] - 0.5).abs() < 1e-9);
        assert!((areas[1] - 3.5).abs() < 1e-9);
    }

    #[test]
    fn a_section_along_the_boundary_is_refused() {
        let along = (Point2::new(0.5, 0.0), Point2::new(1.5, 0.0));
        assert!(split(&[square(2.0)], &[along], 1e-7).is_err());
    }

    #[test]
    fn a_face_with_a_hole_keeps_the_hole_out_of_its_pieces() {
        // Original face: big square with a square hole. No sections. One
        // piece, with the hole attached and the hole's inside excluded.
        let hole: Vec<Point2> = vec![
            Point2::new(1.0, 1.0),
            Point2::new(1.0, 3.0),
            Point2::new(3.0, 3.0),
            Point2::new(3.0, 1.0),
        ];
        let pieces = split(&[square(4.0), hole], &[], 1e-7).unwrap();
        assert_eq!(pieces.len(), 1);
        assert_eq!(pieces[0].rings.len(), 2);
        assert!(inside_rings(&pieces[0].rings, pieces[0].interior));
    }
}
