//! Rounds and blends put on the surfaces their neighbours fix: a round
//! between two planes on the cylinder or cone tangent to both, whether
//! recognition grew it or left it as facets, and fillets and corner balls
//! among rounds on the tori and spheres tangent to their supports.

use ogeom_core::{FastMap, FastSet, Tolerances};
use ogeom_math::{Cylinder, Direction, Frame, Plane, Point, Sphere, Torus, Vector};

use super::planner::{ENCLOSED_CLUSTER, is_sliver, sags_as_the_surface};
use super::seams::distance_to_line;
use super::segment::{axis_frame, leans_as_the_surface, unit_normal};
use super::weld::{Adjacency, from_to};
use super::{Carrier, Curved, Groups};
use crate::recognize::{Canonical, worst_deviation};

/// The normals of the planes the triangles left over from recognition
/// gather into (`planes`, the groups with those planes grown), the largest
/// plane first: a small plane's fit leans with its few vertices' slop.
pub(super) fn plane_normals(planes: &Groups) -> Vec<Direction> {
    let mut size = vec![0_usize; planes.carriers.len()];
    for &g in &planes.of {
        if let Some(n) = size.get_mut(g) {
            *n += 1;
        }
    }
    let mut normals: Vec<(usize, Direction)> = planes
        .carriers
        .iter()
        .zip(&size)
        .filter_map(|(c, &n)| match c {
            Carrier::Plane(plane) => Some((n, plane.frame().z())),
            _ => None,
        })
        .collect();
    normals.sort_by_key(|&(n, _)| core::cmp::Reverse(n));
    normals.into_iter().map(|(_, d)| d).collect()
}

/// Put a round between two flat faces on the cylinder tangent to both.
///
/// A round only a row or two of facets across has its vertices on a few
/// lines, which a cone, or a cylinder leaning from the faces, fits as well
/// as the round's own cylinder; the one fitted then meets the faces along
/// lines that are not where the facets end. The two faces fix the round's
/// axis (along their line of meeting) and leave only the radius, which each
/// vertex gives: the circle through it tangent to both faces. A curved
/// region meeting non-parallel planes (of `planes`, the groups with the
/// planes grown) across smooth edges takes the cylinder tangent to the two
/// largest at the vertices' median radius where it holds every vertex
/// within the distance. With `only`, that region alone is put on its
/// round.
///
/// A long round meshed a few rows across is fitted leaning, and the fit
/// misses pieces of a row by more than the distance: each row is flat, so
/// the pieces are left to planes of a few triangles beside the round, or
/// between it and its flank. Those planes are no flanks; their triangles,
/// where they are free in `groups` and lie on the round put in place as
/// its own facets would, join it. Returns whether any did.
#[allow(clippy::too_many_arguments, reason = "the segmentation's inputs")]
pub(super) fn tangent_rounds(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    planes: &Groups,
    cos_crease: f64,
    flat: f64,
    only: Option<usize>,
    tol: Tolerances,
) -> bool {
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); groups.carriers.len()];
    for (t, &g) in groups.of.iter().enumerate() {
        if let Some(list) = members.get_mut(g) {
            list.push(t);
        }
    }
    let mut plane_members: Vec<Vec<usize>> = vec![Vec::new(); planes.carriers.len()];
    for (t, &g) in planes.of.iter().enumerate() {
        if let Some(list) = plane_members.get_mut(g) {
            list.push(t);
        }
    }
    let area = |tris: &[usize]| -> f64 {
        tris.iter()
            .map(|&t| {
                let [a, b, c] = triangles[t].map(|v| points[v as usize]);
                (b - a).cross(c - a).magnitude() / 2.0
            })
            .sum()
    };
    let mut joined = false;
    for (i, region) in members.iter().enumerate() {
        if only.is_some_and(|o| o != i) {
            continue;
        }
        let Carrier::Curved(curved) = &groups.carriers[i] else {
            continue;
        };
        if matches!(curved.shape, Canonical::Sphere(_) | Canonical::Torus(_)) {
            continue;
        }
        let pts: Vec<Point> = curved
            .vertices
            .iter()
            .map(|&v| points[v as usize])
            .collect();
        // The planes beside the region, the largest in area first: the two
        // largest are its flanks, and any other smaller than both, of a few
        // triangles, is a row of its own facets a leaning fit left out.
        // Where the round on the two largest does not hold the region, the
        // planes of a few triangles past the largest are taken as such rows,
        // and the planes beyond them looked at: a row can stand between the
        // region and its flank. A round found so must hold the region closer
        // than its own fit does.
        let mut reached: FastSet<usize> = region.iter().copied().collect();
        let mut rows: Vec<usize> = Vec::new();
        let mut found = None;
        for _ in 0..3 {
            let mut beside: Vec<usize> = Vec::new();
            for &t in &reached {
                for h in 3 * t..3 * t + 3 {
                    let Some(g) = adjacency.twin[h] else {
                        continue;
                    };
                    let other = g / 3;
                    if reached.contains(&other)
                        || unit_normal(points, triangles[t])
                            .dot(unit_normal(points, triangles[other]))
                            < cos_crease
                    {
                        continue;
                    }
                    beside.push(planes.of[other]);
                }
            }
            beside.sort_unstable();
            beside.dedup();
            let mut flanking: Vec<(usize, &Plane, f64)> = beside
                .iter()
                .filter_map(|&j| match planes.carriers.get(j) {
                    Some(Carrier::Plane(plane)) => Some((j, plane, area(&plane_members[j]))),
                    _ => None,
                })
                .collect();
            flanking.sort_by(|a, b| b.2.total_cmp(&a.2).then(a.0.cmp(&b.0)));
            if flanking.len() < 2 {
                break;
            }
            let few = |j: usize| plane_members[j].len() <= ENCLOSED_CLUSTER;
            let others = &flanking[2..];
            if others.iter().all(|&(j, _, a)| few(j) && a < flanking[1].2)
                && let Some(shape) = tangent_cylinder(flanking[0].1, flanking[1].1, &pts, None, tol)
                && worst_deviation(&shape, &pts)
                    <= if rows.is_empty() {
                        flat
                    } else {
                        flat.min(curved.deviation)
                    }
            {
                rows.extend(others.iter().map(|o| o.0));
                found = Some(shape);
                break;
            }
            let more: Vec<usize> = flanking[1..]
                .iter()
                .map(|f| f.0)
                .filter(|&j| few(j))
                .collect();
            if more.is_empty() {
                break;
            }
            for j in more {
                rows.push(j);
                reached.extend(plane_members[j].iter().copied());
            }
        }
        let Some(shape) = found else {
            continue;
        };
        let deviation = worst_deviation(&shape, &pts);
        let mut taken: Vec<usize> = Vec::new();
        for j in rows {
            let tris = &plane_members[j];
            let on = tris.iter().all(|&t| {
                let corners = triangles[t].map(|v| points[v as usize]);
                groups.of[t] == usize::MAX
                    && corners.iter().all(|p| shape.distance_to(*p) <= flat)
                    && sags_as_the_surface(&shape, corners, flat)
                    && (is_sliver(corners)
                        || leans_as_the_surface(&shape, corners, unit_normal(points, triangles[t])))
            });
            if on {
                taken.extend(tris);
            }
        }
        let Carrier::Curved(curved) = &mut groups.carriers[i] else {
            continue;
        };
        curved.shape = shape;
        curved.deviation = deviation;
        if !taken.is_empty() {
            joined = true;
            for &t in &taken {
                curved.vertices.extend(triangles[t]);
                groups.of[t] = i;
            }
            curved.vertices.sort_unstable();
            curved.vertices.dedup();
            let pts: Vec<Point> = curved
                .vertices
                .iter()
                .map(|&v| points[v as usize])
                .collect();
            curved.deviation = worst_deviation(&curved.shape, &pts);
        }
    }
    joined
}

/// The cylinder tangent to two planes, on the side of them the points are
/// on, at `radius` or else through the points at their median radius.
///
/// With the planes' unit normals `a` and `b` and `w = (a + b) / (1 + a·b)`,
/// the axis of a circle of radius `r` tangent to both runs through
/// `c0 + s·r·w`, `c0` on both planes and `s` the side (-1 within both, +1
/// beyond both). A point `q` from `c0` (square to the axis) is on that
/// circle where `r²(|w|² - 1) - 2s(q·w)r + |q|² = 0`, the larger root.
fn tangent_cylinder(
    a: &Plane,
    b: &Plane,
    pts: &[Point],
    radius: Option<f64>,
    tol: Tolerances,
) -> Option<Canonical> {
    let (na, nb) = (a.frame().z().vector(), b.frame().z().vector());
    let g = na.dot(nb);
    // Nearly parallel planes leave the axis to their slop; nearly opposite
    // ones a radius the points barely fix.
    if g.abs() > 5.0_f64.to_radians().cos() {
        return None;
    }
    let axis = Direction::new(na.cross(nb), tol).ok()?;
    let (pa, pb) = (
        a.frame().origin().to_vector(),
        b.frame().origin().to_vector(),
    );
    // The point on both planes nearest the origin of the axis's normal
    // plane: solve c·na = pa·na, c·nb = pb·nb, c·axis = 0.
    let (ha, hb) = (pa.dot(na), pb.dot(nb));
    let c0 = (nb.cross(axis.vector()) * ha + axis.vector().cross(na) * hb)
        / na.cross(nb).dot(axis.vector());
    let w = (na + nb) / (1.0 + g);
    let side = pts
        .iter()
        .map(|p| (p.to_vector() - c0).dot(na) + (p.to_vector() - c0).dot(nb))
        .sum::<f64>()
        .signum();
    let k = w.dot(w) - 1.0;
    let mut radii: Vec<f64> = pts
        .iter()
        .filter_map(|p| {
            let d = p.to_vector() - c0;
            let q = d - axis.vector() * d.dot(axis.vector());
            let qw = side * q.dot(w);
            let disc = qw * qw - k * q.dot(q);
            (qw > 0.0).then(|| (qw + disc.max(0.0).sqrt()) / k)
        })
        .collect();
    if radii.len() < pts.len() / 2 + 1 {
        return None;
    }
    let radius = match radius {
        Some(radius) => radius,
        None => {
            let at = radii.len() / 2;
            *radii.select_nth_unstable_by(at, f64::total_cmp).1
        }
    };
    let through = Point::from_vector(c0 + w * (side * radius));
    Some(Canonical::Cylinder(
        Cylinder::new(Frame::about(through, axis), radius, tol).ok()?,
    ))
}

/// The most facets a round left faceted between two planes may span and
/// still be put on its surface by [`faceted_rounds`].
const FACETED_ROUND_FACETS: usize = 12;

/// Put a round recognition left as facets between two flat faces tangent
/// to it on its cylinder or cone.
///
/// A round a few facets across and short along its rulings (a drafted
/// wall's corner, met by a fillet at its foot) has too few vertices for
/// any fit to say what it is, and stays a run of planar facets, each a
/// chord between two rulings. The faces either side of the run meet it
/// across smooth edges along rulings, and are tangent to the round there:
/// each such ruling and its face's normal span a plane holding the axis,
/// so the two fix the axis (a cone's through the rulings' meeting point,
/// a cylinder's where they are parallel, square to both normals), and
/// the run's vertices fix the radius. The rulings between the facets
/// verify it: every vertex of the run within `flat` of the surface, and
/// each facet sagging from it as a chord does. A polygon whose sides are
/// all chords, the flanking ones too, is no such round: a surface tangent
/// to its flanking sides along their edges misses its corners.
///
/// Runs are the planar groups of `planes` (the groups with the planes
/// grown) joined across smooth edges, each joined so to two others; runs
/// of two to [`FACETED_ROUND_FACETS`] facets between two planes are tried
/// longest first, none overlapping or flanked by one already taken. Each
/// verified run becomes a curved region of `groups`. Returns whether any
/// did.
#[allow(clippy::too_many_arguments, reason = "the segmentation's inputs")]
pub(super) fn faceted_rounds(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    planes: &Groups,
    cos_crease: f64,
    flat: f64,
    tol: Tolerances,
) -> bool {
    let count = planes.carriers.len();
    let plane = |j: usize| match planes.carriers.get(j) {
        Some(Carrier::Plane(plane)) => Some(*plane),
        _ => None,
    };
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); count];
    for (t, &j) in planes.of.iter().enumerate() {
        if groups.of[t] == usize::MAX && plane(j).is_some() {
            members[j].push(t);
        }
    }
    let mut beside: Vec<Vec<usize>> = vec![Vec::new(); count];
    for (h, twin) in adjacency.twin.iter().enumerate() {
        let Some(g) = *twin else {
            continue;
        };
        let (t, u) = (h / 3, g / 3);
        let (a, b) = (planes.of[t], planes.of[u]);
        if a == b
            || members[a].is_empty()
            || members[b].is_empty()
            || unit_normal(points, triangles[t]).dot(unit_normal(points, triangles[u])) < cos_crease
        {
            continue;
        }
        beside[a].push(b);
    }
    for list in &mut beside {
        list.sort_unstable();
        list.dedup();
    }
    // Chains of facets each joined smoothly to two others, in order, with
    // whether they close on themselves and the groups past their two ends.
    let mut chained = vec![false; count];
    let mut chains: Vec<(Vec<usize>, bool, [usize; 2])> = Vec::new();
    for start in 0..count {
        if chained[start] || beside[start].len() != 2 {
            continue;
        }
        let walk = |from: usize, towards: usize, chained: &mut [bool]| {
            let (mut previous, mut at) = (from, towards);
            let mut run = Vec::new();
            while beside[at].len() == 2 && !chained[at] {
                chained[at] = true;
                run.push(at);
                let next = if beside[at][0] == previous {
                    beside[at][1]
                } else {
                    beside[at][0]
                };
                (previous, at) = (at, next);
            }
            (run, at)
        };
        chained[start] = true;
        let (left, left_end) = walk(start, beside[start][0], &mut chained);
        let (right, right_end) = walk(start, beside[start][1], &mut chained);
        let ring = left_end == start || right_end == start;
        let mut chain: Vec<usize> = left.into_iter().rev().collect();
        chain.push(start);
        chain.extend(right);
        chains.push((chain, ring, [left_end, right_end]));
    }
    // Every window of each chain with what flanks it, longest first.
    let mut windows: Vec<(Vec<usize>, usize, usize)> = Vec::new();
    for (chain, ring, ends) in &chains {
        let m = chain.len();
        let longest = if *ring { m.saturating_sub(2) } else { m };
        for length in 2..=longest.min(FACETED_ROUND_FACETS) {
            let starts = if *ring { m } else { m + 1 - length };
            for i in 0..starts {
                let run: Vec<usize> = (i..i + length).map(|k| chain[k % m]).collect();
                let (before, after) = if *ring {
                    (chain[(i + m - 1) % m], chain[(i + length) % m])
                } else {
                    (
                        if i == 0 { ends[0] } else { chain[i - 1] },
                        if i + length == m {
                            ends[1]
                        } else {
                            chain[i + length]
                        },
                    )
                };
                windows.push((run, before, after));
            }
        }
    }
    windows.sort_by_key(|(run, _, _)| core::cmp::Reverse(run.len()));
    let mut claimed = vec![false; count];
    let mut changed = false;
    let mut sizes = vec![0_usize; count];
    for &j in &planes.of {
        if let Some(n) = sizes.get_mut(j) {
            *n += 1;
        }
    }
    let planes_size = |j: usize| sizes[j];
    for (run, before, after) in windows {
        if [before, after].iter().chain(&run).any(|&j| claimed[j])
            || before == after
            || (planes_size(before) < 2 && planes_size(after) < 2)
        {
            continue;
        }
        let (Some(a), Some(b)) = (plane(before), plane(after)) else {
            continue;
        };
        let tris: Vec<usize> = run.iter().flat_map(|&j| members[j].clone()).collect();
        let Some(shape) = round_between(
            points,
            triangles,
            adjacency,
            planes,
            &tris,
            [(before, a), (after, b)],
            flat,
            tol,
        ) else {
            continue;
        };
        let mut vertices: Vec<u32> = tris.iter().flat_map(|&t| triangles[t]).collect();
        vertices.sort_unstable();
        vertices.dedup();
        let pts: Vec<Point> = vertices.iter().map(|&v| points[v as usize]).collect();
        let deviation = worst_deviation(&shape, &pts);
        let g = groups.carriers.len();
        for &t in &tris {
            groups.of[t] = g;
        }
        groups.carriers.push(Carrier::Curved(Curved {
            shape,
            deviation,
            fitted: deviation,
            centre: (0.0, 0.0),
            wraps: false,
            wraps_v: false,
            fixed: false,
            vertices,
            patch: None,
        }));
        for &j in &run {
            claimed[j] = true;
        }
        changed = true;
    }
    changed
}

/// The cylinder or cone through the triangles `tris` tangent to the two
/// flanking planes along the rulings the triangles share with them, where
/// it holds every vertex within `flat` and every triangle sags from it as
/// a chord: see [`faceted_rounds`]. Each flank is its group in `planes`
/// and its plane.
#[allow(clippy::too_many_arguments, reason = "the run and its flanks")]
fn round_between(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    planes: &Groups,
    tris: &[usize],
    flanks: [(usize, Plane); 2],
    flat: f64,
    tol: Tolerances,
) -> Option<Canonical> {
    // The ruling the run shares with a flank: the line through the
    // farthest two of the shared edges' ends, holding the others.
    let ruling = |flank: usize| -> Option<(Point, Vector)> {
        let mut ends: Vec<Point> = Vec::new();
        for &t in tris {
            for h in 3 * t..3 * t + 3 {
                if adjacency.twin[h].is_some_and(|g| planes.of[g / 3] == flank) {
                    let (p, q) = from_to(triangles, h);
                    ends.extend([points[p as usize], points[q as usize]]);
                }
            }
        }
        let (mut far, mut length) = ((*ends.first()?, *ends.first()?), 0.0);
        for (i, &p) in ends.iter().enumerate() {
            for &q in &ends[i + 1..] {
                if p.distance(q) > length {
                    (far, length) = ((p, q), p.distance(q));
                }
            }
        }
        if length <= flat
            || ends
                .iter()
                .any(|&p| distance_to_line(p, far.0, far.1) > flat)
        {
            return None;
        }
        Some((far.0, (far.1 - far.0) / length))
    };
    let (p1, d1) = ruling(flanks[0].0)?;
    let (p2, d2) = ruling(flanks[1].0)?;
    let (n1, n2) = (
        flanks[0].1.frame().z().vector(),
        flanks[1].1.frame().z().vector(),
    );
    let mut vertices: Vec<u32> = tris.iter().flat_map(|&t| triangles[t]).collect();
    vertices.sort_unstable();
    vertices.dedup();
    let pts: Vec<Point> = vertices.iter().map(|&v| points[v as usize]).collect();
    let holds = |shape: &Canonical| {
        worst_deviation(shape, &pts) <= flat
            && tris.iter().all(|&t| {
                sags_as_the_surface(shape, triangles[t].map(|v| points[v as usize]), flat)
            })
    };
    // A cone's axis lies in the plane of each tangent ruling and its face's
    // normal, and passes through where the rulings meet.
    let axis = n1.cross(d1).cross(n2.cross(d2));
    let across = d1.cross(d2);
    if across.magnitude() > tol.angular() && axis.magnitude() > tol.angular() {
        // The rulings' closest points: their meeting point, the apex.
        let w = p1 - p2;
        let (b, d, e) = (d1.dot(d2), d1.dot(w), d2.dot(w));
        let s = (b * e - d) / (1.0 - b * b);
        let apex = p1 + d1 * s;
        if let Ok(axis) = Direction::new(axis, tol)
            && let Some(shape) = crate::recognize::ruled_about(&pts, apex, axis, flat, tol)
            && holds(&shape)
        {
            return Some(shape);
        }
    }
    tangent_cylinder(&flanks[0].1, &flanks[1].1, &pts, None, tol).filter(holds)
}

/// Put fillets and corner balls on the surfaces their neighbours fix.
///
/// A fillet is fitted as freely as any region, so it meets the faces it
/// blends into at a slight angle or a slight gap, and the seam solved
/// between two nearly tangent surfaces wanders along them. A fillet's
/// supports fix it but for its radius:
///
/// - rounds between two planes that meet at corner balls are one rolling
///   ball's: they take their median radius together where their own
///   radii agree to within `flat`, and each ball is centred a radius off
///   the planes its rounds run between, the point all their axes pass
///   through;
/// - a torus between a plane and a cylinder or a cone whose axis is square
///   to the plane sits on that axis, its tube's centre a radius off both;
///   so does a torus between any other two supports on its axis: two of
///   cylinders, cones and planes square to it, or one of them and a
///   sphere centred on it ([`torus_between`]);
/// - a sphere where cylinders of its own radius meet otherwise is centred
///   nearest their axes; a sphere beside one plane and a cylinder or cone
///   whose axis is square to it is a piece of the torus between them, met
///   a few facets round (each piece's vertices lie on a sphere as exactly
///   as on the torus), where that torus holds it.
///
/// The surface derived so is tangent to its supports by construction, and
/// replaces the fitted one where it holds every vertex of the region within
/// `flat`, on the fitted one's frame so the chart branch fixed for it still
/// holds. Other regions keep their fits. `planes` are the groups with the
/// planes grown, as for [`tangent_rounds`]. With `only`, the fillets are
/// derived as for all of them and that region alone is put on its own.
/// Returns whether a sphere was put on a torus.
#[allow(clippy::too_many_arguments, reason = "the segmentation's inputs")]
pub(super) fn tangent_blends(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    planes: &Groups,
    cos_crease: f64,
    flat: f64,
    only: Option<usize>,
    tol: Tolerances,
) -> bool {
    let count = groups.carriers.len();
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); count];
    for (t, &g) in groups.of.iter().enumerate() {
        if let Some(list) = members.get_mut(g) {
            list.push(t);
        }
    }
    // Each curved region's planes and curved regions across its smooth
    // edges, and its vertices.
    let mut flanks: Vec<Vec<Plane>> = vec![Vec::new(); count];
    let mut supports: Vec<Vec<usize>> = vec![Vec::new(); count];
    let mut samples: Vec<Vec<Point>> = vec![Vec::new(); count];
    for (i, region) in members.iter().enumerate() {
        let Carrier::Curved(curved) = &groups.carriers[i] else {
            continue;
        };
        let mut beside: Vec<usize> = Vec::new();
        for &t in region {
            for h in 3 * t..3 * t + 3 {
                let Some(g) = adjacency.twin[h] else {
                    continue;
                };
                let other = g / 3;
                if groups.of[other] == i
                    || unit_normal(points, triangles[t]).dot(unit_normal(points, triangles[other]))
                        < cos_crease
                {
                    continue;
                }
                beside.push(planes.of[other]);
            }
        }
        beside.sort_unstable();
        beside.dedup();
        for &j in &beside {
            match (planes.carriers.get(j), groups.carriers.get(j)) {
                (Some(Carrier::Plane(plane)), _) => flanks[i].push(*plane),
                (_, Some(Carrier::Curved(_))) if j != i => supports[i].push(j),
                _ => {}
            }
        }
        samples[i] = curved
            .vertices
            .iter()
            .map(|&v| points[v as usize])
            .collect();
    }
    let taken = |i: usize| only.is_none_or(|o| o == i);
    let shape_of = |groups: &Groups, i: usize| match &groups.carriers[i] {
        Carrier::Curved(c) => Some(c.shape.clone()),
        _ => None,
    };
    let is_round = |groups: &Groups, i: usize| {
        flanks[i].len() == 2 && matches!(shape_of(groups, i), Some(Canonical::Cylinder(_)))
    };
    let is_ball =
        |groups: &Groups, i: usize| matches!(shape_of(groups, i), Some(Canonical::Sphere(_)));
    // Rounds joined through the balls they meet at.
    let mut chain: Vec<usize> = (0..count).collect();
    fn root(chain: &mut [usize], mut i: usize) -> usize {
        while chain[i] != i {
            chain[i] = chain[chain[i]];
            i = chain[i];
        }
        i
    }
    for (i, beside) in supports.iter().enumerate() {
        if !is_ball(groups, i) {
            continue;
        }
        for &j in beside {
            if is_round(groups, j) {
                let (a, b) = (root(&mut chain, i), root(&mut chain, j));
                chain[a] = b;
            }
        }
    }
    let mut parts: FastMap<usize, Vec<usize>> = FastMap::default();
    for i in 0..count {
        if is_round(groups, i) || is_ball(groups, i) {
            parts.entry(root(&mut chain, i)).or_default().push(i);
        }
    }
    let mut settled = vec![false; count];
    for part in parts.values() {
        let rounds: Vec<usize> = part
            .iter()
            .copied()
            .filter(|&j| is_round(groups, j))
            .collect();
        let balls: Vec<usize> = part
            .iter()
            .copied()
            .filter(|&j| is_ball(groups, j))
            .collect();
        if balls.is_empty() || rounds.is_empty() {
            continue;
        }
        let mut radii: Vec<f64> = rounds
            .iter()
            .filter_map(|&j| match shape_of(groups, j) {
                Some(Canonical::Cylinder(c)) => Some(c.radius()),
                _ => None,
            })
            .collect();
        let (lo, hi) = radii
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &r| {
                (lo.min(r), hi.max(r))
            });
        if hi - lo > flat {
            continue;
        }
        let at = radii.len() / 2;
        let radius = *radii.select_nth_unstable_by(at, f64::total_cmp).1;
        for &j in &rounds {
            let Some(old) = shape_of(groups, j) else {
                continue;
            };
            let derived =
                tangent_cylinder(&flanks[j][0], &flanks[j][1], &samples[j], Some(radius), tol)
                    .and_then(|d| on_frame_of(&d, &old, tol));
            if let Some(shape) = derived
                && taken(j)
            {
                settled[j] = put(groups, j, shape, &samples[j], flat);
            }
        }
        for &j in &balls {
            let Some(Canonical::Sphere(old)) = shape_of(groups, j) else {
                continue;
            };
            let mut around: Vec<Plane> = Vec::new();
            for &k in &supports[j] {
                if rounds.contains(&k) {
                    around.extend(flanks[k].iter().copied());
                }
            }
            if let Some(shape) = ball_off_planes(&old, &around, &samples[j], radius, tol)
                && taken(j)
            {
                settled[j] = put(groups, j, shape, &samples[j], flat);
            }
        }
    }
    let mut tori = false;
    for i in 0..count {
        if settled[i] || !taken(i) {
            continue;
        }
        let derived = match shape_of(groups, i) {
            Some(Canonical::Torus(torus)) => {
                // The support is the band on the torus's own axis; other
                // fillets running into it smoothly are its neighbours, not
                // its supports.
                let (o, z) = (torus.frame().origin(), torus.frame().z().vector());
                let coaxial: Vec<Canonical> = supports[i]
                    .iter()
                    .filter_map(|&j| shape_of(groups, j))
                    .filter(|s| {
                        axis_frame(s).is_some_and(|f| {
                            let w = f.origin() - o;
                            f.z().vector().cross(z).magnitude() <= 1e-3
                                && (w - z * w.dot(z)).magnitude() <= flat * 10.0
                        })
                    })
                    .collect();
                let around: Vec<Canonical> = supports[i]
                    .iter()
                    .filter_map(|&j| shape_of(groups, j))
                    .collect();
                match (&flanks[i][..], &coaxial[..]) {
                    ([plane], [support]) => tangent_torus(plane, support, &samples[i], tol),
                    _ => None,
                }
                .or_else(|| {
                    torus_between(&torus, &flanks[i], &around, &samples[i], flat * 10.0, tol)
                })
            }
            Some(Canonical::Sphere(sphere)) => {
                let around: Vec<Canonical> = supports[i]
                    .iter()
                    .filter_map(|&j| shape_of(groups, j))
                    .collect();
                let holds = |torus: Canonical| {
                    (worst_deviation(&torus, &samples[i]) <= flat).then_some(torus)
                };
                corner_ball(&sphere, &around, flat, tol)
                    .or_else(|| {
                        let [plane] = &flanks[i][..] else {
                            return None;
                        };
                        let torus = around.iter().find_map(|support| {
                            holds(tangent_torus(plane, support, &samples[i], tol)?)
                        })?;
                        tori = true;
                        Some(torus)
                    })
                    .or_else(|| {
                        let torus = flanks[i].iter().find_map(|plane| {
                            supports[i].iter().find_map(|&j| {
                                let ends: Vec<Point> = samples[j]
                                    .iter()
                                    .copied()
                                    .filter(|p| samples[i].contains(p))
                                    .collect();
                                let round = shape_of(groups, j)?;
                                holds(turning_round(plane, &round, &samples[i], &ends, flat, tol)?)
                            })
                        })?;
                        tori = true;
                        Some(torus)
                    })
            }
            _ => None,
        };
        if let Some(shape) = derived {
            put(groups, i, shape, &samples[i], flat);
        }
    }
    tori
}

/// `shape` as group `i`'s surface where it holds every sample within
/// `flat`; whether it does.
fn put(groups: &mut Groups, i: usize, shape: Canonical, samples: &[Point], flat: f64) -> bool {
    let deviation = worst_deviation(&shape, samples);
    if deviation > flat {
        return false;
    }
    if let Carrier::Curved(curved) = &mut groups.carriers[i] {
        curved.shape = shape;
        curved.deviation = deviation;
    }
    true
}

/// A cylinder put on the frame of the one it replaces: its axis turned to
/// run the same way, its angle measured from the same direction and its
/// height from the same level, so the chart branch fixed for the old one
/// reads the same on the new.
fn on_frame_of(shape: &Canonical, old: &Canonical, tol: Tolerances) -> Option<Canonical> {
    let (Canonical::Cylinder(new), Some(was)) = (shape, axis_frame(old)) else {
        return None;
    };
    let z = new.frame().z();
    let z = if z.vector().dot(was.z().vector()) >= 0.0 {
        z
    } else {
        -z
    };
    let x = was.x().vector() - z.vector() * was.x().vector().dot(z.vector());
    let origin =
        new.frame().origin() + z.vector() * (was.origin() - new.frame().origin()).dot(z.vector());
    let frame = Frame::new(origin, z, Direction::new(x, tol).ok()?, tol).ok()?;
    Some(Canonical::Cylinder(
        Cylinder::new(frame, new.radius(), tol).ok()?,
    ))
}

/// The ball a radius off the planes its rounds run between, on the side of
/// each the samples are on, on the sphere's own frame. `None` unless the
/// planes are three that meet in a point.
fn ball_off_planes(
    sphere: &Sphere,
    around: &[Plane],
    samples: &[Point],
    radius: f64,
    tol: Tolerances,
) -> Option<Canonical> {
    let mut distinct: Vec<Plane> = Vec::new();
    for plane in around {
        let n = plane.frame().z().vector();
        let same = distinct.iter().any(|q| {
            let m = q.frame().z().vector();
            n.cross(m).magnitude() <= 1e-9
                && (plane.frame().origin() - q.frame().origin()).dot(m).abs() <= tol.confusion()
        });
        if !same {
            distinct.push(*plane);
        }
    }
    let [a, b, c] = distinct[..] else {
        return None;
    };
    #[allow(clippy::cast_precision_loss, reason = "vertex counts are small")]
    let count = samples.len().max(1) as f64;
    let mut m = nalgebra::Matrix3::<f64>::zeros();
    let mut rhs = nalgebra::Vector3::<f64>::zeros();
    for (row, plane) in [a, b, c].iter().enumerate() {
        let (o, n) = (plane.frame().origin(), plane.frame().z().vector());
        let side = (samples.iter().map(|p| (*p - o).dot(n)).sum::<f64>() / count).signum();
        m.set_row(row, &nalgebra::RowVector3::new(n.x, n.y, n.z));
        rhs[row] = o.to_vector().dot(n) + side * radius;
    }
    let solved = m.lu().solve(&rhs)?;
    let centre = Point::new(solved[0], solved[1], solved[2]);
    let frame = sphere.frame();
    let on = Frame::new(centre, frame.z(), frame.x(), tol).ok()?;
    Some(Canonical::Sphere(Sphere::new(on, radius, tol).ok()?))
}

/// The torus tangent to a plane and to a cylinder or cone whose axis is
/// square to it, on the side of each the points are on, through the
/// points at their median tube radius.
///
/// In the support's axial half-plane, with `rho` the distance from the axis
/// and `z` the height along it, the support's line is `rho = a + b z` and
/// the plane is `z = h`. A tube circle of radius `r` tangent to both has its
/// centre at `A + B r`, with `A = (a + b h, h)` and
/// `B = (s2 √(1 + b²) + s1 b, s1)`, `s1` and `s2` the sides of the plane and
/// the support the points are on. A point `q` lies on that circle where
/// `(|B|² - 1) r² - 2 (q - A)·B r + |q - A|² = 0`, the larger root.
fn tangent_torus(
    plane: &Plane,
    support: &Canonical,
    pts: &[Point],
    tol: Tolerances,
) -> Option<Canonical> {
    let (frame, a, b) = match support {
        Canonical::Cylinder(c) => (c.frame(), c.radius(), 0.0),
        Canonical::Cone(c) => (
            c.frame(),
            c.radius_at(0.0),
            c.radius_at(1.0) - c.radius_at(0.0),
        ),
        _ => return None,
    };
    let (o, z) = (frame.origin(), frame.z().vector());
    let n = plane.frame().z().vector();
    if n.cross(z).magnitude() > 1e-9 {
        return None;
    }
    let h = (plane.frame().origin() - o).dot(z);
    let local: Vec<(f64, f64)> = pts
        .iter()
        .map(|p| {
            let w = *p - o;
            let along = w.dot(z);
            ((w - z * along).magnitude(), along)
        })
        .collect();
    #[allow(clippy::cast_precision_loss, reason = "vertex counts are small")]
    let count = local.len().max(1) as f64;
    let s1 = (local.iter().map(|q| q.1 - h).sum::<f64>() / count).signum();
    let s2 = (local.iter().map(|q| q.0 - a - b * q.1).sum::<f64>() / count).signum();
    let base = (a + b * h, h);
    let step = (s2 * b.hypot(1.0) + s1 * b, s1);
    let (centre, radius) = rolled_between_lines(&local, base, step)?;
    if centre.0 <= radius {
        return None;
    }
    let on = Frame::new(o + z * centre.1, frame.z(), frame.x(), tol).ok()?;
    Some(Canonical::Torus(
        Torus::new(on, centre.0, radius, tol).ok()?,
    ))
}

/// The tube circle tangent to two crossing lines of the axial half-plane,
/// its centre at `base + step r`, through `local` at their median radius:
/// a point `q` lies on it where
/// `(|step|² - 1) r² - 2 (q - base)·step r + |q - base|² = 0`, the larger
/// root (the arc between the tangent points faces the lines' crossing).
/// The centre and the radius; `None` where fewer than half the points
/// have a root.
fn rolled_between_lines(
    local: &[(f64, f64)],
    base: (f64, f64),
    step: (f64, f64),
) -> Option<((f64, f64), f64)> {
    let k = step.0.mul_add(step.0, step.1 * step.1) - 1.0;
    if k <= 0.0 {
        return None;
    }
    let mut radii: Vec<f64> = local
        .iter()
        .filter_map(|q| {
            let d = (q.0 - base.0, q.1 - base.1);
            let along = d.0.mul_add(step.0, d.1 * step.1);
            let disc = along * along - k * d.0.mul_add(d.0, d.1 * d.1);
            (along > 0.0).then(|| (along + disc.max(0.0).sqrt()) / k)
        })
        .collect();
    if radii.len() < local.len() / 2 + 1 {
        return None;
    }
    let at = radii.len() / 2;
    let (_, radius, _) = radii.select_nth_unstable_by(at, f64::total_cmp);
    let radius = *radius;
    let centre = (
        step.0.mul_add(radius, base.0),
        step.1.mul_add(radius, base.1),
    );
    Some((centre, radius))
}

/// A support's trace in the half-plane through an axis, a point there
/// being `(rho, h)`, its distance from the axis and its height along it:
/// the line `n · p = d` (`n` a unit vector), or a circle.
#[derive(Clone, Copy)]
enum Trace {
    Line((f64, f64), f64),
    Circle((f64, f64), f64),
}

impl Trace {
    /// How far `q` stands off the trace, signed: positive off the side
    /// `n` points to, or outside the circle.
    fn side(self, q: (f64, f64)) -> f64 {
        match self {
            Self::Line(n, d) => n.0.mul_add(q.0, n.1 * q.1) - d,
            Self::Circle(c, radius) => (q.0 - c.0).hypot(q.1 - c.1) - radius,
        }
    }
}

/// `support`'s trace about the axis through `o` along `z`: a plane square
/// to the axis, or a cylinder, cone or sphere on it. `None` for any other
/// surface or placing.
fn trace_about(support: &Canonical, o: Point, z: Direction, tol: Tolerances) -> Option<Trace> {
    let z = z.vector();
    let square = |d: Direction| d.vector().cross(z).magnitude() <= 1e-9;
    let on_axis = |p: Point| {
        let w = p - o;
        (w - z * w.dot(z)).magnitude() <= tol.confusion()
    };
    match support {
        Canonical::Plane(plane) if square(plane.frame().z()) => {
            Some(Trace::Line((0.0, 1.0), (plane.frame().origin() - o).dot(z)))
        }
        Canonical::Cylinder(c) if square(c.frame().z()) && on_axis(c.frame().origin()) => {
            Some(Trace::Line((1.0, 0.0), c.radius()))
        }
        Canonical::Cone(c) if square(c.frame().z()) && on_axis(c.frame().origin()) => {
            // `rho = a + slope h` along this axis.
            let slope =
                (c.radius_at(1.0) - c.radius_at(0.0)) * c.frame().z().vector().dot(z).signum();
            let a = c.radius_at(0.0) - slope * (c.frame().origin() - o).dot(z);
            let k = slope.hypot(1.0);
            Some(Trace::Line((1.0 / k, -slope / k), a / k))
        }
        Canonical::Sphere(s) if on_axis(s.centre()) => {
            Some(Trace::Circle((0.0, (s.centre() - o).dot(z)), s.radius()))
        }
        _ => None,
    }
}

/// The torus tangent to the two supports on its axis where they are other
/// than one plane and one cylinder or cone ([`tangent_torus`]): two of
/// cylinders, cones and planes square to the axis, or one of them and a
/// sphere centred on it. The axis is a cylinder's or cone's among them, or
/// a sphere's centre, or the fitted torus's, turned square to a plane
/// among them. `near` is how far off the fitted torus's axis a support's
/// may stand and still be on it; supports off it are other fillets
/// running into this one. `None` where there are not exactly two
/// supports on the axis, or no such torus.
///
/// A tube circle of radius `r` tangent to two traces, on the side of each
/// the points are on, has its centre `c` at `n · c = d + s r` for a line
/// and `|c - centre| = radius + s r` for a circle, `s` the side. Across
/// two crossing lines `c = A + B r` and the radius is found as by
/// [`rolled_between_lines`]. Between two parallel lines the radius is half
/// their distance and the centre slides along them, placed where the
/// points put it on one side of the tube (the arc is a half circle). For a
/// line and a circle the centre at each `r` is where the offset line and
/// circle meet, the meeting nearer the fitted tube's centre taken, and the
/// radius the one the points stand closest to in the least squares.
fn torus_between(
    torus: &Torus,
    planes: &[Plane],
    around: &[Canonical],
    pts: &[Point],
    near: f64,
    tol: Tolerances,
) -> Option<Canonical> {
    let fitted = torus.frame();
    let tz = fitted.z().vector();
    let off_axis = |p: Point| {
        let w = p - fitted.origin();
        (w - tz * w.dot(tz)).magnitude()
    };
    let parallel = |d: Direction| d.vector().cross(tz).magnitude() <= 1e-3;
    let mut supports: Vec<Canonical> = planes
        .iter()
        .filter(|p| parallel(p.frame().z()))
        .map(|p| Canonical::Plane(*p))
        .collect();
    for shape in around {
        let on = match shape {
            Canonical::Sphere(s) => off_axis(s.centre()) <= near,
            _ => axis_frame(shape).is_some_and(|f| parallel(f.z()) && off_axis(f.origin()) <= near),
        };
        if on {
            supports.push(shape.clone());
        }
    }
    let [first, second] = &supports[..] else {
        return None;
    };
    let normal = supports.iter().find_map(|s| match s {
        Canonical::Plane(p) => Some(p.frame().z()),
        _ => None,
    });
    let (o, z) = if let Some(f) = supports.iter().find_map(|s| match s {
        Canonical::Cylinder(_) | Canonical::Cone(_) => axis_frame(s),
        _ => None,
    }) {
        (f.origin(), f.z())
    } else {
        let o = supports
            .iter()
            .find_map(|s| match s {
                Canonical::Sphere(s) => Some(s.centre()),
                _ => None,
            })
            .unwrap_or(fitted.origin());
        (o, normal.unwrap_or(fitted.z()))
    };
    let traces = [
        trace_about(first, o, z, tol)?,
        trace_about(second, o, z, tol)?,
    ];
    let zv = z.vector();
    let local: Vec<(f64, f64)> = pts
        .iter()
        .map(|p| {
            let w = *p - o;
            let along = w.dot(zv);
            ((w - zv * along).magnitude(), along)
        })
        .collect();
    #[allow(clippy::cast_precision_loss, reason = "vertex counts are small")]
    let count = local.len().max(1) as f64;
    let sides = traces.map(|t| (local.iter().map(|&q| t.side(q)).sum::<f64>() / count).signum());
    let (centre, radius) = match traces {
        [Trace::Line(n1, d1), Trace::Line(n2, d2)] => {
            let det = n1.0.mul_add(n2.1, -(n1.1 * n2.0));
            if det.abs() > 1e-9 {
                let solve = |e1: f64, e2: f64| {
                    (
                        e1.mul_add(n2.1, -(n1.1 * e2)) / det,
                        n1.0.mul_add(e2, -(e1 * n2.0)) / det,
                    )
                };
                rolled_between_lines(&local, solve(d1, d2), solve(sides[0], sides[1]))?
            } else {
                rolled_between_parallels((n1, d1), (n2, d2), sides, &local)?
            }
        }
        [Trace::Line(n, d), Trace::Circle(c, r)] => {
            rolled_by_circle((n, d, sides[0]), (c, r, sides[1]), torus, o, zv, &local)?
        }
        [Trace::Circle(c, r), Trace::Line(n, d)] => {
            rolled_by_circle((n, d, sides[1]), (c, r, sides[0]), torus, o, zv, &local)?
        }
        [Trace::Circle(..), Trace::Circle(..)] => return None,
    };
    if !radius.is_finite() || radius <= 0.0 || centre.0 <= radius {
        return None;
    }
    let sense = if fitted.z().vector().dot(zv) >= 0.0 {
        z
    } else {
        -z
    };
    let x = fitted.x().vector() - zv * fitted.x().vector().dot(zv);
    let mut on = Frame::new(o + zv * centre.1, sense, Direction::new(x, tol).ok()?, tol).ok()?;
    if on.handedness() != fitted.handedness() {
        on = on.mirrored();
    }
    Some(Canonical::Torus(
        Torus::new(on, centre.0, radius, tol).ok()?,
    ))
}

/// The tube circle between two parallel lines of the axial half-plane,
/// tangent to both: its radius half their distance, its centre on the line
/// midway and placed along it by the points, which lie on one side of the
/// tube. Each point puts the centre at two places, one per side; the side
/// whose places agree the more closely is taken, at their median.
fn rolled_between_parallels(
    (n1, d1): ((f64, f64), f64),
    (n2, d2): ((f64, f64), f64),
    sides: [f64; 2],
    local: &[(f64, f64)],
) -> Option<((f64, f64), f64)> {
    let flip = n1.0.mul_add(n2.0, n1.1 * n2.1).signum();
    let across = sides[0] - flip * sides[1];
    if across == 0.0 {
        return None;
    }
    let radius = flip.mul_add(d2, -d1) / across;
    if radius <= 0.0 {
        return None;
    }
    let level = sides[0].mul_add(radius, d1);
    let m = (-n1.1, n1.0);
    let mut places = [Vec::new(), Vec::new()];
    for q in local {
        let off = n1.0.mul_add(q.0, n1.1 * q.1) - level;
        let along = m.0.mul_add(q.0, m.1 * q.1);
        let half = radius.mul_add(radius, -(off * off)).max(0.0).sqrt();
        places[0].push(along - half);
        places[1].push(along + half);
    }
    let settled = places.map(|mut t: Vec<f64>| {
        if t.is_empty() {
            return (f64::INFINITY, 0.0);
        }
        let at = t.len() / 2;
        let median = *t.select_nth_unstable_by(at, f64::total_cmp).1;
        let spread = t.iter().map(|v| (v - median).abs()).fold(0.0, f64::max);
        (spread, median)
    });
    let (_, t) = if settled[0].0 <= settled[1].0 {
        settled[0]
    } else {
        settled[1]
    };
    let centre = (n1.0.mul_add(level, m.0 * t), n1.1.mul_add(level, m.1 * t));
    Some((centre, radius))
}

/// The tube circle tangent to a line `n · p = d` and a circle of the axial
/// half-plane, on sides `s` of each, through `local` in the least squares.
/// Its centre at radius `r` is where the line offset by `s r` meets the
/// circle offset by `s r`, the meeting nearer the fitted torus's tube
/// centre; `r` is searched over half to twice the fitted tube's.
fn rolled_by_circle(
    (n, d, s1): ((f64, f64), f64, f64),
    (c, big, s2): ((f64, f64), f64, f64),
    torus: &Torus,
    o: Point,
    z: Vector,
    local: &[(f64, f64)],
) -> Option<((f64, f64), f64)> {
    let w = torus.frame().origin() - o;
    let fitted = (torus.major_radius(), w.dot(z));
    let m = (-n.1, n.0);
    let (cn, cm) = (n.0.mul_add(c.0, n.1 * c.1), m.0.mul_add(c.0, m.1 * c.1));
    let centre_at = |r: f64| -> Option<(f64, f64)> {
        let level = s1.mul_add(r, d);
        let reach = s2.mul_add(r, big);
        let disc = reach.mul_add(reach, -((level - cn) * (level - cn)));
        if reach <= 0.0 || disc < 0.0 {
            return None;
        }
        let at = |t: f64| (n.0.mul_add(level, m.0 * t), n.1.mul_add(level, m.1 * t));
        let (a, b) = (at(cm - disc.sqrt()), at(cm + disc.sqrt()));
        let gap = |p: (f64, f64)| (p.0 - fitted.0).hypot(p.1 - fitted.1);
        Some(if gap(a) <= gap(b) { a } else { b })
    };
    let cost = |r: f64| {
        centre_at(r).map_or(f64::INFINITY, |p| {
            local
                .iter()
                .map(|q| {
                    let e = (q.0 - p.0).hypot(q.1 - p.1) - r;
                    e * e
                })
                .sum::<f64>()
        })
    };
    let r0 = torus.minor_radius();
    let (lo, hi) = (r0 * 0.5, r0 * 2.0);
    let steps = 64;
    let width = (hi - lo) / f64::from(steps);
    let mut best = (f64::INFINITY, r0);
    for i in 0..=steps {
        let r = width.mul_add(f64::from(i), lo);
        let e = cost(r);
        if e < best.0 {
            best = (e, r);
        }
    }
    if !best.0.is_finite() {
        return None;
    }
    let (mut a, mut b) = ((best.1 - width).max(lo), (best.1 + width).min(hi));
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    for _ in 0..80 {
        let (p, q) = (b - ratio * (b - a), a + ratio * (b - a));
        if cost(p) <= cost(q) {
            b = q;
        } else {
            a = p;
        }
    }
    let radius = (a + b) / 2.0;
    Some((centre_at(radius)?, radius))
}

/// The torus a round along a plane turns on where the plane's edge it
/// follows turns: its axis square to the plane, its tube the round's, and
/// the circle of its tube's centres tangent to the round's axis (which
/// lies at the round's radius from the plane, within `flat`) where the
/// round ends. `ends` are the points the round and `pts` share, its end
/// section: they fix where along the axis the round ends. `None` where the
/// round does not lie so, or its end is no section square to its axis.
///
/// With `a` a point of the round's axis, `d` its direction, `n` the plane's
/// normal and `e = n × d`, a point `q` is at `(u, v, z)` along `d`, `e` and
/// `n` from `a`. The torus centred at `a + t d + R e`, of major radius `|R|`
/// and tube `r`, holds `q` where `(u - t)² + v² - s² = 2 R (v ± s)`,
/// `s = √(r² - z²)`, the sign that of `R` times the side of the tube the
/// point is on: `R` is solved by least squares for each sign, the one
/// holding the points closer taken. Without the end fixing `t` the points
/// do not tell the torus from a sphere through them (`R = 0`), which a
/// few facets of it fit as closely.
fn turning_round(
    plane: &Plane,
    support: &Canonical,
    pts: &[Point],
    ends: &[Point],
    flat: f64,
    tol: Tolerances,
) -> Option<Canonical> {
    let Canonical::Cylinder(round) = support else {
        return None;
    };
    let (a, d, r) = (
        round.frame().origin(),
        round.frame().z().vector(),
        round.radius(),
    );
    let n = plane.frame().z().vector();
    if d.dot(n).abs() > 1e-9
        || ((a - plane.frame().origin()).dot(n).abs() - r).abs() > flat
        || ends.len() < 2
    {
        return None;
    }
    let along: Vec<f64> = ends.iter().map(|p| (*p - a).dot(d)).collect();
    #[allow(clippy::cast_precision_loss, reason = "vertex counts are small")]
    let t = along.iter().sum::<f64>() / along.len() as f64;
    if along.iter().any(|u| (u - t).abs() > flat) {
        return None;
    }
    let e = n.cross(d);
    let local: Vec<(f64, f64, f64)> = pts
        .iter()
        .map(|p| {
            let w = *p - a;
            let z = w.dot(n).clamp(-r, r);
            (w.dot(d) - t, w.dot(e), (r * r - z * z).sqrt())
        })
        .collect();
    let mut best: Option<(f64, Canonical)> = None;
    for side in [1.0_f64, -1.0] {
        let (mut cc, mut cy) = (0.0, 0.0);
        for &(u, v, s) in &local {
            let c = 2.0 * side.mul_add(s, v);
            cc += c * c;
            cy += c * (u.mul_add(u, v * v) - s * s);
        }
        if cc <= flat * flat {
            continue;
        }
        let major = cy / cc;
        if major.abs() <= r {
            continue;
        }
        let centre = a + d * t + e * major;
        let Ok(frame) = Frame::new(centre, plane.frame().z(), round.frame().z(), tol) else {
            continue;
        };
        let Ok(torus) = Torus::new(frame, major.abs(), r, tol) else {
            continue;
        };
        let torus = Canonical::Torus(torus);
        let off = worst_deviation(&torus, pts);
        if best.as_ref().is_none_or(|(o, _)| off < *o) {
            best = Some((off, torus));
        }
    }
    best.map(|(_, torus)| torus)
}

/// The ball where cylinders of its own radius meet: centred at the point
/// nearest their axes, its radius theirs, on the sphere's own frame. The
/// cylinders taken are those within a hundredth of the fitted sphere's
/// radius; `None` unless there are at least two, their radii agree to
/// within `flat`, their axes are not parallel, and every axis passes within
/// `flat` of the centre.
fn corner_ball(
    sphere: &Sphere,
    supports: &[Canonical],
    flat: f64,
    tol: Tolerances,
) -> Option<Canonical> {
    let rounds: Vec<&Cylinder> = supports
        .iter()
        .filter_map(|s| match s {
            Canonical::Cylinder(c)
                if (c.radius() - sphere.radius()).abs() <= sphere.radius() * 0.01 =>
            {
                Some(c)
            }
            _ => None,
        })
        .collect();
    if rounds.len() < 2 {
        return None;
    }
    let (lo, hi) = rounds
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), c| {
            (lo.min(c.radius()), hi.max(c.radius()))
        });
    if hi - lo > flat {
        return None;
    }
    // Least squares: the sum over the axes of the projections square to
    // each, applied to the centre, equals the same applied to their origins.
    let mut m = nalgebra::Matrix3::<f64>::zeros();
    let mut rhs = nalgebra::Vector3::<f64>::zeros();
    for c in &rounds {
        let z = c.frame().z().vector();
        let zv = nalgebra::Vector3::new(z.x, z.y, z.z);
        let across = nalgebra::Matrix3::identity() - zv * zv.transpose();
        let o = c.frame().origin();
        m += across;
        rhs += across * nalgebra::Vector3::new(o.x, o.y, o.z);
    }
    let eigen = m.symmetric_eigen();
    if eigen.eigenvalues.min() < 1e-3 {
        return None;
    }
    let solved = m.lu().solve(&rhs)?;
    let centre = Point::new(solved[0], solved[1], solved[2]);
    let met = rounds.iter().all(|c| {
        let z = c.frame().z().vector();
        let w = centre - c.frame().origin();
        (w - z * w.dot(z)).magnitude() <= flat
    });
    if !met {
        return None;
    }
    #[allow(clippy::cast_precision_loss, reason = "a few cylinders")]
    let radius = rounds.iter().map(|c| c.radius()).sum::<f64>() / rounds.len() as f64;
    let frame = sphere.frame();
    let on = Frame::new(centre, frame.z(), frame.x(), tol).ok()?;
    Some(Canonical::Sphere(Sphere::new(on, radius, tol).ok()?))
}
