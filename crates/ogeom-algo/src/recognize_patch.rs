//! A fitted B-spline patch for a smooth region of a mesh that no plane,
//! canonical surface or sweep fits.
//!
//! Only a region that is one disk bounded by one loop is tried, and only
//! one with vertices enough inside its boundary to hold a surface across
//! it. Its chart (a parameter pair per vertex, in the unit square) comes
//! from a canonical surface that nearly fits it, where one does and its
//! chart does not fold the region over; otherwise from the mean-value map
//! of the region onto the unit square (Floater 2003): the boundary laid on
//! the square's sides by arc length, the square's corners at four of its
//! vertices, and each interior vertex the weighted mean of its neighbours,
//! one sparse linear system solved in a band. The map of a disk onto a
//! convex boundary by positive weights folds nowhere. Where the patch over
//! that map does not verify, the chart is the region's projection onto the
//! plane across its mean normal, if that folds nowhere: a region with a
//! rounded boundary, running out tangentially into its neighbours, is
//! pinched into the square's corners by the map but not by the projection.
//!
//! A cubic patch over the unit square is fitted to the vertices at those
//! parameters, faired by a small thin-plate term. After each fit every
//! vertex's parameters move to its foot on the patch; knots are inserted
//! in the spans holding vertices whose feet are farther than half the
//! distance, no finer than the vertices' own spacing in the chart, until
//! none is, and three more rounds of fitting and moving follow at those
//! knots.
//!
//! It is then verified both ways: every vertex lies within the distance of
//! the patch; at points inside every triangle the triangle lies within the
//! distance plus its own sag of the patch, and the patch within that of
//! the triangle's plane; the patch's Jacobian does not vanish over the
//! region, and its normal agrees with every triangle's. A patch that fails
//! any of it is refused, and the region stays faceted. A verified patch is
//! continued a little past the square on every side, for the seams solved
//! onto it.

use std::collections::{HashMap, HashSet};

use ogeom_core::Tolerances;
use ogeom_geom::{Surface as _, SurfaceGeometry};
use ogeom_math::{KnotVector, Point, Vector};

use crate::recognize::Canonical;

/// A verified patch.
#[derive(Debug, Clone)]
pub(crate) struct Patch {
    /// The patch, a [`SurfaceGeometry::BSpline`].
    pub(crate) surface: SurfaceGeometry,
    /// The worst distance from a vertex of the region to it.
    pub(crate) deviation: f64,
    /// Whether its chart is the mean-value map, rather than a nearly
    /// fitting canonical surface's.
    pub(crate) mapped: bool,
}

/// Why a region got no patch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Refused {
    /// The region is not one disk bounded by one loop.
    NotDisk,
    /// The region is a disk too narrow to have a surface of its own: so few
    /// of its vertices lie inside its boundary that nothing between its rims
    /// would hold the patch (a row or two of facets left along a seam).
    Narrow,
    /// No patch fitted to it passed the verification.
    Unverified,
}

/// A region of a mesh to fit.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Region<'a> {
    pub(crate) points: &'a [Point],
    pub(crate) triangles: &'a [[u32; 3]],
    /// The region's triangles, by index.
    pub(crate) members: &'a [usize],
    /// The region's boundary vertices at which the face across its boundary
    /// changes: where the square's corners go first.
    pub(crate) corners: &'a HashSet<u32>,
}

/// The fewest vertices a region must hold inside its boundary: a few dozen
/// triangles left between recognized faces (a fillet's run-out, a corner)
/// say too little about a surface across them, and the seams round them
/// cost more than their facets.
const MIN_INSIDE: usize = 32;

/// A region must hold at least one vertex inside its boundary for this many
/// on it.
const INSIDE_SHARE: usize = 4;

/// How much farther than the distance a canonical surface may stand off
/// the region and still lend it its chart.
const NEAR: f64 = 10.0;

/// The share of the distance the fit is held to at every vertex, leaving
/// the rest to the triangles' interiors and the foot's own search.
const FIT_SHARE: f64 = 0.5;

/// How far past the unit square the verified patch is continued on every
/// side, against the region's size: the seams solved onto it lie a little
/// off the region's boundary, and a patch that stopped at that boundary
/// would stop short of them. It is fitted on the square alone, where the
/// vertices hold every control; the continuation is each side's own Taylor
/// polynomial to second order.
const EXTENSION: f64 = 0.02;

/// The thin-plate term's weight against the mean squared distance, in the
/// unit square's parameters, over the square of the fit's target against the
/// region's size: small enough that where vertices hold the patch it moves
/// it by less than the target, and where none do (a span they leave empty,
/// the square's corners beyond a canonical chart's region) it alone
/// settles the patch, at any weight.
const FAIRING: f64 = 1.0;

/// How far below its first weight the fairing is lowered, a hundredfold at a
/// time, where the knots run out before the vertices are held.
const FAIRING_FLOOR: f64 = 1e-4;

/// The most knot spans the patch takes along either direction.
const MAX_SPANS: usize = 128;

/// How many of the vertices' steps in the chart a knot span must hold to be
/// split: the halves keep one each.
const SPAN_STEPS: f64 = 2.0;

/// How many times the parameters are moved to their feet on the patch.
const CORRECTIONS: usize = 3;

/// How far a patch's normal may lean from a triangle's past the turn of the
/// patch's normals across the triangle: about eleven degrees, as for any
/// recognized surface.
const LEAN: f64 = 0.2;

/// How far a boundary vertex must turn, past straight on, to take a corner
/// of the square where the neighbouring face does not change there: thirty
/// degrees, a crease's turn.
const SHARP: f64 = core::f64::consts::FRAC_PI_6;

/// The largest banded solve the mean-value map is allowed, in multiplies.
const MAP_WORK: f64 = 2e9;

#[allow(clippy::cast_precision_loss, reason = "counts are far below 2^52")]
fn count(n: usize) -> f64 {
    n as f64
}

/// The region with local vertex numbers: each local vertex's mesh index,
/// and the triangles in local numbers.
struct Local {
    vertices: Vec<u32>,
    triangles: Vec<[usize; 3]>,
    points: Vec<Point>,
}

/// The region's local mesh, and its one boundary loop with the region on
/// its left; `None` where the region is not one disk bounded by one loop.
fn disk(region: &Region) -> Option<(Local, Vec<usize>)> {
    let mut index: HashMap<u32, usize> = HashMap::new();
    let mut vertices = Vec::new();
    let mut triangles = Vec::with_capacity(region.members.len());
    for &t in region.members {
        let tri = region.triangles[t].map(|v| {
            *index.entry(v).or_insert_with(|| {
                vertices.push(v);
                vertices.len() - 1
            })
        });
        triangles.push(tri);
    }
    // Each directed edge at most once: a second use is a fold or a seam
    // three triangles share.
    let mut directed: HashSet<(usize, usize)> = HashSet::new();
    for t in &triangles {
        for k in 0..3 {
            if !directed.insert((t[k], t[(k + 1) % 3])) {
                return None;
            }
        }
    }
    let edges = directed
        .iter()
        .filter(|&&(a, b)| a < b || !directed.contains(&(b, a)))
        .count();
    // The boundary runs along the directed edges with no twin, the region
    // on their left; at each vertex one leaves, or the boundary touches
    // itself there.
    let mut next: HashMap<usize, usize> = HashMap::new();
    for &(a, b) in &directed {
        if !directed.contains(&(b, a)) && next.insert(a, b).is_some() {
            return None;
        }
    }
    let start = *next.keys().min()?;
    let mut ring = vec![start];
    let mut at = next[&start];
    while at != start {
        ring.push(at);
        at = *next.get(&at)?;
        if ring.len() > next.len() {
            return None;
        }
    }
    // One loop, and Euler's count of a disk.
    let euler = isize::try_from(vertices.len()).ok()? - isize::try_from(edges).ok()?
        + isize::try_from(triangles.len()).ok()?;
    if ring.len() != next.len() || euler != 1 || ring.len() < 4 {
        return None;
    }
    let points = vertices
        .iter()
        .map(|&v| region.points[v as usize])
        .collect();
    Some((
        Local {
            vertices,
            triangles,
            points,
        },
        ring,
    ))
}

/// Fit and verify a patch over the region, at the coplanar distance `flat`.
pub(crate) fn fit_patch(region: &Region, flat: f64, tol: Tolerances) -> Result<Patch, Refused> {
    let Some((local, ring)) = disk(region) else {
        return Err(Refused::NotDisk);
    };
    let inside = local.points.len() - ring.len();
    if inside < MIN_INSIDE || inside * INSIDE_SHARE < ring.len() {
        return Err(Refused::Narrow);
    }
    let normals: Vec<Vector> = local
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| local.points[i]);
            let n = (b - a).cross(c - a);
            let m = n.magnitude();
            if m > 0.0 { n / m } else { Vector::ZERO }
        })
        .collect();
    let mut vertex_normals = vec![Vector::ZERO; local.points.len()];
    for (t, n) in local.triangles.iter().zip(&normals) {
        for &i in t {
            vertex_normals[i] += *n;
        }
    }
    for n in &mut vertex_normals {
        let m = n.magnitude();
        *n = if m > 0.0 { *n / m } else { Vector::Z };
    }
    let near = crate::recognize::recognize_points(&local.points, &vertex_normals, flat * NEAR, tol)
        .ok()
        .flatten();
    if let Some(found) = near
        && let Some(chart) = canonical_chart(&found.surface, &local, tol)
        && let Some((surface, deviation)) = fitted(&local, &normals, chart, flat, tol)
    {
        return Ok(Patch {
            surface,
            deviation,
            mapped: false,
        });
    }
    let corners: HashSet<usize> = (0..local.vertices.len())
        .filter(|&i| region.corners.contains(&local.vertices[i]))
        .collect();
    if let Some(chart) = mean_value_map(&local, &ring, &corners)
        && let Some((surface, deviation)) = fitted(&local, &normals, chart, flat, tol)
    {
        return Ok(Patch {
            surface,
            deviation,
            mapped: true,
        });
    }
    // A region that is a height over the plane across its mean normal
    // projects onto that plane without folding (the chart's winding says
    // whether it does) and keeps the surface's own spacing there, where the
    // mean-value map pinches a rounded boundary (a region running out
    // tangentially into its neighbours) into the square's corners.
    let plane = mean_plane(&local, &normals, tol).ok_or(Refused::Unverified)?;
    let chart =
        canonical_chart(&Canonical::Plane(plane), &local, tol).ok_or(Refused::Unverified)?;
    let (surface, deviation) =
        fitted(&local, &normals, chart, flat, tol).ok_or(Refused::Unverified)?;
    Ok(Patch {
        surface,
        deviation,
        mapped: false,
    })
}

/// The plane through the region's centroid across its triangles' mean
/// normal, weighted by area; `None` where the region turns so far that the
/// mean normal vanishes.
fn mean_plane(local: &Local, normals: &[Vector], tol: Tolerances) -> Option<ogeom_math::Plane> {
    let mut mean = Vector::ZERO;
    for (t, n) in local.triangles.iter().zip(normals) {
        let [a, b, c] = t.map(|i| local.points[i]);
        mean += *n * (b - a).cross(c - a).magnitude();
    }
    let z = ogeom_math::Direction::new(mean, tol).ok()?;
    let mut centre = Vector::ZERO;
    for p in &local.points {
        centre += p.to_vector();
    }
    let centre = Point::from_vector(centre / count(local.points.len()));
    Some(ogeom_math::Plane::new(ogeom_math::Frame::about(centre, z)))
}

/// The region's chart on a canonical surface, normalized to the unit
/// square and turned so the parameter triangles wind as the mesh's do.
/// `None` where a vertex has no chart position (a pole), the region goes
/// all the way round the surface, or the chart folds it over.
fn canonical_chart(shape: &Canonical, local: &Local, tol: Tolerances) -> Option<Vec<(f64, f64)>> {
    use ogeom_math::elementary as e;
    let raw: Option<Vec<(f64, f64)>> = local
        .points
        .iter()
        .map(|&p| match shape {
            Canonical::Plane(plane) => Some(e::plane_parameters(plane, p)),
            Canonical::Cylinder(c) => e::cylinder_parameters(c, p, tol).ok(),
            Canonical::Cone(c) => e::cone_parameters(c, p, tol).ok(),
            Canonical::Sphere(s) => e::sphere_parameters(s, p, tol).ok(),
            Canonical::Torus(t) => e::torus_parameters(t, p, tol).ok(),
            Canonical::Swept(_) => None,
        })
        .collect();
    let mut chart = raw?;
    let (wraps_u, wraps_v) = match shape {
        Canonical::Cylinder(_) | Canonical::Cone(_) | Canonical::Sphere(_) => (true, false),
        Canonical::Torus(_) => (true, true),
        Canonical::Plane(_) | Canonical::Swept(_) => (false, false),
    };
    // An angle is read on the branch opposite the widest gap between the
    // vertices' angles; a gap of under a tenth of a turn says the region
    // goes round.
    let unwrap = |values: &mut Vec<f64>| -> Option<()> {
        let tau = core::f64::consts::TAU;
        let mut sorted: Vec<f64> = values.iter().map(|x| x.rem_euclid(tau)).collect();
        sorted.sort_by(f64::total_cmp);
        let (mut gap, mut after) = (sorted[0] + tau - sorted[sorted.len() - 1], sorted[0]);
        for w in sorted.windows(2) {
            if w[1] - w[0] > gap {
                gap = w[1] - w[0];
                after = w[1];
            }
        }
        if gap < tau * 0.1 {
            return None;
        }
        for x in values.iter_mut() {
            *x = after + (*x - after).rem_euclid(tau);
        }
        Some(())
    };
    for (wraps, pick) in [(wraps_u, 0), (wraps_v, 1)] {
        if !wraps {
            continue;
        }
        let mut values: Vec<f64> = chart
            .iter()
            .map(|c| if pick == 0 { c.0 } else { c.1 })
            .collect();
        unwrap(&mut values)?;
        for (c, x) in chart.iter_mut().zip(values) {
            if pick == 0 {
                c.0 = x;
            } else {
                c.1 = x;
            }
        }
    }
    let span = |pick: fn(&(f64, f64)) -> f64| {
        chart
            .iter()
            .map(pick)
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
                (lo.min(x), hi.max(x))
            })
    };
    let (u0, u1) = span(|c| c.0);
    let (v0, v1) = span(|c| c.1);
    if u1 - u0 <= 0.0 || v1 - v0 <= 0.0 {
        return None;
    }
    let mut chart: Vec<(f64, f64)> = chart
        .iter()
        .map(|&(u, v)| ((u - u0) / (u1 - u0), (v - v0) / (v1 - v0)))
        .collect();
    match winding(&local.triangles, &chart) {
        Some(true) => {}
        Some(false) => {
            for c in &mut chart {
                c.0 = 1.0 - c.0;
            }
        }
        None => return None,
    }
    Some(chart)
}

/// Whether every parameter triangle winds the mesh's way (`true`) or every
/// one the other way (`false`); `None` where they disagree or one has no
/// area, the chart folding the region over.
fn winding(triangles: &[[usize; 3]], chart: &[(f64, f64)]) -> Option<bool> {
    let area = |t: &[usize; 3]| {
        let [a, b, c] = t.map(|i| chart[i]);
        (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
    };
    let areas: Vec<f64> = triangles.iter().map(area).collect();
    let typical = areas.iter().map(|a| a.abs()).sum::<f64>() / count(areas.len().max(1));
    let floor = typical * 1e-9;
    if areas.iter().all(|&a| a > floor) {
        Some(true)
    } else if areas.iter().all(|&a| a < -floor) {
        Some(false)
    } else {
        None
    }
}

/// The mean-value map of the region onto the unit square: its boundary on
/// the square's sides by arc length, the corners at the boundary vertices
/// where the face across the boundary changes (the four sharpest of them)
/// or else at the sharpest of the rest; each interior vertex the mean of
/// its neighbours under Floater's weights, found by one banded solve in an
/// ordering that keeps the band narrow. `None` where the solve is too wide
/// to make or the map folds.
fn mean_value_map(
    local: &Local,
    ring: &[usize],
    corners: &HashSet<usize>,
) -> Option<Vec<(f64, f64)>> {
    let n = local.points.len();
    let angle_at = |t: &[usize; 3], k: usize| {
        let [a, b, c] = [t[k], t[(k + 1) % 3], t[(k + 2) % 3]].map(|i| local.points[i]);
        let (x, y) = (b - a, c - a);
        x.cross(y).magnitude().atan2(x.dot(y))
    };
    // The region's angle at each vertex, and Floater's weights: across
    // each triangle at vertex i, tan(angle / 2) over the length of each of
    // its two sides from i.
    let mut inside = vec![0.0_f64; n];
    let mut weights: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for t in &local.triangles {
        for k in 0..3 {
            let i = t[k];
            let angle = angle_at(t, k);
            inside[i] += angle;
            let half = (angle / 2.0).tan();
            for j in [t[(k + 1) % 3], t[(k + 2) % 3]] {
                let length = local.points[i].distance(local.points[j]);
                if length > 0.0 {
                    weights[i].push((j, half / length));
                }
            }
        }
    }
    // The four corners, in the ring's order: the sharpest of the vertices
    // where the neighbouring face changes; then the sharpest of the rest
    // that turn by a crease's worth, a ring vertex at least between any two
    // corners; then, while fewer than four, the vertex halfway along the
    // longest stretch of the ring between corners.
    let mut by_angle: Vec<usize> = (0..ring.len()).collect();
    by_angle.sort_by(|&a, &b| inside[ring[a]].total_cmp(&inside[ring[b]]));
    let mut chosen: Vec<usize> = by_angle
        .iter()
        .copied()
        .filter(|&k| corners.contains(&ring[k]))
        .take(4)
        .collect();
    let apart = |k: usize, chosen: &[usize]| {
        chosen.iter().all(|&c| {
            let d = k.abs_diff(c);
            d.min(ring.len() - d) >= 2
        })
    };
    for &k in &by_angle {
        if chosen.len() == 4 || core::f64::consts::PI - inside[ring[k]] < SHARP {
            break;
        }
        if !chosen.contains(&k) && apart(k, &chosen) {
            chosen.push(k);
        }
    }
    let mut along = vec![0.0];
    for k in 1..=ring.len() {
        let last = along[k - 1];
        along.push(last + local.points[ring[k - 1]].distance(local.points[ring[k % ring.len()]]));
    }
    let perimeter = along[ring.len()];
    if chosen.is_empty() {
        chosen.push(by_angle[0]);
    }
    while chosen.len() < 4 {
        chosen.sort_unstable();
        // The longest stretch, from a corner to the next round the ring.
        let (start, length) = (0..chosen.len())
            .map(|c| {
                let (a, b) = (chosen[c], chosen[(c + 1) % chosen.len()]);
                let mut stretch = along[b] - along[a];
                if stretch <= 0.0 {
                    stretch += perimeter;
                }
                (a, stretch)
            })
            .max_by(|x, y| x.1.total_cmp(&y.1))?;
        let middle = (along[start] + length / 2.0) % perimeter;
        let k = (0..ring.len())
            .filter(|k| !chosen.contains(k))
            .min_by(|&x, &y| {
                let off = |k: usize| {
                    let d = (along[k] - middle).abs();
                    d.min(perimeter - d)
                };
                off(x).total_cmp(&off(y))
            })?;
        chosen.push(k);
    }
    chosen.sort_unstable();
    // The boundary laid on the square counter-clockwise, side by side, each
    // by its own arc length.
    let squares = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let mut chart = vec![(f64::NAN, f64::NAN); n];
    for side in 0..4 {
        let (from, to) = (chosen[side], chosen[(side + 1) % 4]);
        let steps = (to + ring.len() - from) % ring.len();
        let run: Vec<usize> = (0..=steps).map(|s| ring[(from + s) % ring.len()]).collect();
        let mut lengths = vec![0.0];
        for w in run.windows(2) {
            let last = lengths[lengths.len() - 1];
            lengths.push(last + local.points[w[0]].distance(local.points[w[1]]));
        }
        let total = lengths[lengths.len() - 1];
        if total <= 0.0 {
            return None;
        }
        let (a, b) = (squares[side], squares[(side + 1) % 4]);
        for (i, s) in run.iter().zip(&lengths) {
            let f = s / total;
            chart[*i] = (a.0 + (b.0 - a.0) * f, a.1 + (b.1 - a.1) * f);
        }
    }
    // The interior: Σ w_ij (x_i - x_j) = 0 at each interior vertex.
    let on_ring: HashSet<usize> = ring.iter().copied().collect();
    let interior: Vec<usize> = (0..n).filter(|i| !on_ring.contains(i)).collect();
    if !interior.is_empty() {
        let order = narrow_order(&interior, &weights);
        let mut slot = vec![usize::MAX; n];
        for (s, &i) in order.iter().enumerate() {
            slot[i] = s;
        }
        let slot = slot;
        let mut band = 0;
        for &i in &order {
            for &(j, _) in &weights[i] {
                if slot[j] != usize::MAX {
                    band = band.max(slot[i].abs_diff(slot[j]));
                }
            }
        }
        let size = order.len();
        if count(size) * count(band + 1).powi(2) > MAP_WORK {
            return None;
        }
        let width = 2 * band + 1;
        let mut matrix = vec![0.0; size * width];
        let mut rhs = vec![[0.0; 2]; size];
        for (s, &i) in order.iter().enumerate() {
            for &(j, w) in &weights[i] {
                matrix[s * width + band] += w;
                if slot[j] == usize::MAX {
                    rhs[s][0] += w * chart[j].0;
                    rhs[s][1] += w * chart[j].1;
                } else {
                    matrix[s * width + band + slot[j] - s] -= w;
                }
            }
        }
        let solved = banded_solve(&mut matrix, &mut rhs, band)?;
        for (s, &i) in order.iter().enumerate() {
            chart[i] = (solved[s][0], solved[s][1]);
        }
    }
    (winding(&local.triangles, &chart) == Some(true)).then_some(chart)
}

/// The interior vertices in reverse Cuthill-McKee order: breadth first
/// from a vertex of least degree in each connected piece, neighbours by
/// rising degree, the whole reversed. Neighbours then sit close in the
/// order, and the system's band is narrow.
fn narrow_order(interior: &[usize], weights: &[Vec<(usize, f64)>]) -> Vec<usize> {
    let member: HashSet<usize> = interior.iter().copied().collect();
    let neighbours = |i: usize| -> Vec<usize> {
        let mut out: Vec<usize> = weights[i]
            .iter()
            .map(|&(j, _)| j)
            .filter(|j| member.contains(j))
            .collect();
        out.sort_unstable();
        out.dedup();
        out
    };
    let degree: HashMap<usize, usize> =
        interior.iter().map(|&i| (i, neighbours(i).len())).collect();
    let mut by_degree: Vec<usize> = interior.to_vec();
    by_degree.sort_by_key(|i| (degree[i], *i));
    let mut seen: HashSet<usize> = HashSet::new();
    let mut order = Vec::with_capacity(interior.len());
    for &start in &by_degree {
        if !seen.insert(start) {
            continue;
        }
        let mut queue = std::collections::VecDeque::from([start]);
        while let Some(i) = queue.pop_front() {
            order.push(i);
            let mut next = neighbours(i);
            next.sort_by_key(|j| (degree[j], *j));
            for j in next {
                if seen.insert(j) {
                    queue.push_back(j);
                }
            }
        }
    }
    order.reverse();
    order
}

/// Solve a banded system by elimination without pivoting, which the
/// mean-value system (diagonally dominant, irreducibly so through the rows
/// beside the boundary) does not need. The matrix holds each row's band,
/// `2 band + 1` wide, the diagonal in the middle. `None` at a vanishing
/// pivot.
fn banded_solve(matrix: &mut [f64], rhs: &mut [[f64; 2]], band: usize) -> Option<Vec<[f64; 2]>> {
    let n = rhs.len();
    let width = 2 * band + 1;
    // Entry (r, c) of the full matrix.
    let at = |r: usize, c: usize| r * width + band + c - r;
    for c in 0..n {
        let pivot = matrix[at(c, c)];
        if pivot.abs() <= f64::MIN_POSITIVE || !pivot.is_finite() {
            return None;
        }
        for r in c + 1..n.min(c + band + 1) {
            let factor = matrix[at(r, c)] / pivot;
            if factor == 0.0 {
                continue;
            }
            for k in c..n.min(c + band + 1) {
                matrix[at(r, k)] -= factor * matrix[at(c, k)];
            }
            rhs[r][0] -= factor * rhs[c][0];
            rhs[r][1] -= factor * rhs[c][1];
        }
    }
    let mut x = vec![[0.0; 2]; n];
    for r in (0..n).rev() {
        let mut sum = rhs[r];
        for k in r + 1..n.min(r + band + 1) {
            sum[0] -= matrix[at(r, k)] * x[k][0];
            sum[1] -= matrix[at(r, k)] * x[k][1];
        }
        let pivot = matrix[at(r, r)];
        x[r] = [sum[0] / pivot, sum[1] / pivot];
    }
    x.iter()
        .all(|v| v[0].is_finite() && v[1].is_finite())
        .then_some(x)
}

/// A clamped cubic knot vector over the domain with `spans` equal spans.
fn uniform_knots(spans: usize) -> Option<KnotVector> {
    let (lo, hi) = (0.0, 1.0);
    let mut knots = vec![lo; 4];
    for k in 1..spans {
        knots.push(lo + (hi - lo) * count(k) / count(spans));
    }
    knots.extend([hi; 4]);
    KnotVector::new(knots, 3).ok()
}

/// The knot vector with a knot inserted in the middle of each span listed.
fn split(knots: &KnotVector, spans: &[usize]) -> Option<KnotVector> {
    let distinct = knots.distinct();
    let mut out = knots.clone();
    for &s in spans {
        let (a, b) = (distinct[s].0, distinct[s + 1].0);
        out = out.with_knot_inserted(f64::midpoint(a, b), 1).ok()?;
    }
    Some(out)
}

/// Which of a knot vector's distinct spans holds `t`.
fn span_of(knots: &KnotVector, t: f64) -> usize {
    let distinct = knots.distinct();
    distinct
        .windows(2)
        .position(|w| t < w[1].0)
        .unwrap_or(distinct.len().saturating_sub(2))
}

/// The patch fitted at the chart's parameters, and its worst vertex
/// distance, once it passes the verification, continued past the square on
/// every side; `None` where it does not pass, or
/// the knots run out before the vertices are held.
fn fitted(
    local: &Local,
    normals: &[Vector],
    chart: Vec<(f64, f64)>,
    flat: f64,
    tol: Tolerances,
) -> Option<(SurfaceGeometry, f64)> {
    let points = &local.points;
    let target = flat * FIT_SHARE;
    let start = (count(points.len()).sqrt() / 8.0).ceil().clamp(1.0, 8.0);
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "a small positive span count"
    )]
    let start = start as usize;
    let (mut u_knots, mut v_knots) = (uniform_knots(start)?, uniform_knots(start)?);
    // The vertices' spacing in the chart along each direction: the median
    // of the mesh's edges that run more that way than the other. A span
    // is split only while it holds two such steps, so the knots never
    // outrun the vertices that settle them.
    let (mut along_u, mut along_v) = (Vec::new(), Vec::new());
    for t in &local.triangles {
        for k in 0..3 {
            let (a, b) = (chart[t[k]], chart[t[(k + 1) % 3]]);
            let (du, dv) = ((a.0 - b.0).abs(), (a.1 - b.1).abs());
            if du >= dv {
                along_u.push(du);
            } else {
                along_v.push(dv);
            }
        }
    }
    let median = |mut steps: Vec<f64>| -> f64 {
        if steps.is_empty() {
            return 0.0;
        }
        steps.sort_by(f64::total_cmp);
        steps[steps.len() / 2]
    };
    let (step_u, step_v) = (median(along_u), median(along_v));
    let mut parameters = chart;
    let size = {
        let (mut lo, mut hi) = (points[0], points[0]);
        for p in points {
            lo = Point::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Point::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        lo.distance(hi)
    };
    let first = FAIRING * (target / size).powi(2);
    let mut fairing = first;
    // Each round fits at the parameters, then moves every vertex's
    // parameters to its foot on the fit: that takes out the error along the
    // surface the chart's own spacing leaves, and the distance across it
    // (the foot's) is what the knots answer. The rounds end once every foot
    // has stayed within the target through the corrections at fixed knots.
    let mut settled = 0;
    let surface = loop {
        let fit = ogeom_geom::fit::fit_surface_scattered_at(
            points,
            &parameters,
            &u_knots,
            &v_knots,
            fairing,
            tol,
        )
        .ok()?;
        let geometry = SurfaceGeometry::BSpline(fit.curve);
        let mut distances = Vec::with_capacity(points.len());
        for (p, at) in points.iter().zip(parameters.iter_mut()) {
            match crate::measure::project_on_surface_from(&geometry, *p, *at, tol) {
                Ok(foot) => {
                    *at = foot.parameters;
                    distances.push(foot.distance);
                }
                Err(_) => distances.push(f64::INFINITY),
            }
        }
        if distances.iter().all(|&d| d <= target) {
            settled += 1;
            if settled > CORRECTIONS {
                break geometry;
            }
            continue;
        }
        settled = 0;
        // The spans holding a vertex past the target, each way.
        let mut u_bad: Vec<usize> = Vec::new();
        let mut v_bad: Vec<usize> = Vec::new();
        for (&(u, v), d) in parameters.iter().zip(&distances) {
            if *d > target {
                u_bad.push(span_of(&u_knots, u));
                v_bad.push(span_of(&v_knots, v));
            }
        }
        for list in [&mut u_bad, &mut v_bad] {
            list.sort_unstable();
            list.dedup();
            list.reverse();
        }
        let wide = |knots: &KnotVector, s: usize, step: f64| {
            let distinct = knots.distinct();
            distinct[s + 1].0 - distinct[s].0 >= step * SPAN_STEPS
        };
        u_bad.retain(|&s| wide(&u_knots, s, step_u));
        v_bad.retain(|&s| wide(&v_knots, s, step_v));
        let room = |knots: &KnotVector, more: usize| knots.distinct().len() - 1 + more <= MAX_SPANS;
        if (u_bad.is_empty() && v_bad.is_empty())
            || !room(&u_knots, u_bad.len())
            || !room(&v_knots, v_bad.len())
        {
            if fairing <= first * FAIRING_FLOOR {
                return None;
            }
            fairing /= 100.0;
            continue;
        }
        u_knots = split(&u_knots, &u_bad)?;
        v_knots = split(&v_knots, &v_bad)?;
    };
    let deviation = verified(local, normals, &surface, &parameters, flat, tol)?;
    let SurfaceGeometry::BSpline(mut patch) = surface else {
        return None;
    };
    for (along_u, at_end) in [(true, false), (true, true), (false, false), (false, true)] {
        patch = patch
            .extended(along_u, at_end, size * EXTENSION, 2, tol)
            .ok()?;
    }
    Some((SurfaceGeometry::BSpline(patch), deviation))
}

/// The patch's worst distance from a vertex, where it passes the checks
/// both ways; `None` where it does not.
fn verified(
    local: &Local,
    normals: &[Vector],
    surface: &SurfaceGeometry,
    parameters: &[(f64, f64)],
    flat: f64,
    tol: Tolerances,
) -> Option<f64> {
    let foot = |p: Point, at: (f64, f64)| {
        crate::measure::project_on_surface_from(surface, p, at, tol).ok()
    };
    // Every vertex, at its foot.
    let mut deviation = 0.0_f64;
    let mut feet = Vec::with_capacity(parameters.len());
    for (p, &at) in local.points.iter().zip(parameters) {
        let found = foot(*p, at)?;
        if found.distance > flat {
            return None;
        }
        deviation = deviation.max(found.distance);
        feet.push(found.parameters);
    }
    // The feet must keep the chart unfolded.
    if winding(&local.triangles, &feet) != Some(true) {
        return None;
    }
    let unit = |(u, v): (f64, f64)| -> Option<(Vector, f64)> {
        let (du, dv) = surface.d1_at(u, v, tol).ok()?;
        let n = du.cross(dv);
        let m = n.magnitude();
        (m > 0.0 && m.is_finite()).then(|| (n / m, m))
    };
    let mut jacobians: Vec<f64> = Vec::with_capacity(feet.len());
    for &at in &feet {
        jacobians.push(unit(at)?.1);
    }
    let typical = {
        let mut sorted = jacobians.clone();
        sorted.sort_by(f64::total_cmp);
        sorted[sorted.len() / 2]
    };
    // A Jacobian a millionth of the typical one is a pinch of the patch.
    let pinch = typical * 1e-6;
    if jacobians.iter().any(|&j| j <= pinch) {
        return None;
    }
    let angle = |a: Vector, b: Vector| a.dot(b).clamp(-1.0, 1.0).acos();
    const INSIDE: [[f64; 3]; 7] = [
        [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        [0.5, 0.5, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
    ];
    for (t, normal) in local.triangles.iter().zip(normals) {
        let corners = t.map(|i| local.points[i]);
        let at = t.map(|i| feet[i]);
        let mut turned = [Vector::ZERO; 3];
        for (n, &uv) in turned.iter_mut().zip(&at) {
            *n = unit(uv)?.0;
        }
        let turn = angle(turned[0], turned[1])
            .max(angle(turned[1], turned[2]))
            .max(angle(turned[0], turned[2]));
        let longest = corners[0]
            .distance(corners[1])
            .max(corners[1].distance(corners[2]))
            .max(corners[0].distance(corners[2]));
        // The triangle's own sag: a chord of that length over that much
        // turn of the surface stands off it by about this much.
        let allowed = flat + longest * turn / 6.0;
        for b in INSIDE {
            let q = Point::from_vector(
                corners[0].to_vector() * b[0]
                    + corners[1].to_vector() * b[1]
                    + corners[2].to_vector() * b[2],
            );
            let uv = (
                at[0].0 * b[0] + at[1].0 * b[1] + at[2].0 * b[2],
                at[0].1 * b[0] + at[1].1 * b[1] + at[2].1 * b[2],
            );
            // The triangle onto the patch, and the patch onto the
            // triangle: the patch's point at the same blend of the corners'
            // parameters stands off the triangle's plane by no more, and
            // lies over the triangle or beside it. The chart need not be
            // affine across a triangle, so that point may stand a little
            // along the plane from the blend of the corners, which is no
            // distance of the surface from the mesh.
            if foot(q, uv)?.distance > allowed {
                return None;
            }
            let over = surface.point_at(uv.0, uv.1, tol).ok()?;
            if (over - q).dot(*normal).abs() > allowed
                || crate::mesh_solid::distance_to_triangle(over, corners[0], corners[1], corners[2])
                    > longest
            {
                return None;
            }
            let (n, j) = unit(uv)?;
            if j <= pinch {
                return None;
            }
            // The patch's normal leans from the triangle's no more than
            // the patch turns under it, and a little more.
            if *normal != Vector::ZERO && angle(n, *normal) > turn + LEAN {
                return None;
            }
        }
    }
    Some(deviation)
}

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code")]
mod tests {
    use super::*;

    const T: Tolerances = Tolerances::millimetres();

    /// A grid of `n` by `n` cells over the unit square lifted by `lift`,
    /// each cell two triangles wound counter-clockwise seen from `+z`.
    fn sheet(n: usize, lift: impl Fn(f64, f64) -> Point) -> (Vec<Point>, Vec<[u32; 3]>) {
        let mut points = Vec::new();
        for j in 0..=n {
            for i in 0..=n {
                points.push(lift(count(i) / count(n), count(j) / count(n)));
            }
        }
        let at = |i: usize, j: usize| u32::try_from(j * (n + 1) + i).unwrap();
        let mut triangles = Vec::new();
        for j in 0..n {
            for i in 0..n {
                triangles.push([at(i, j), at(i + 1, j), at(i + 1, j + 1)]);
                triangles.push([at(i, j), at(i + 1, j + 1), at(i, j + 1)]);
            }
        }
        (points, triangles)
    }

    /// The patch's worst distance from the points.
    fn worst(patch: &Patch, points: &[Point]) -> f64 {
        let shape = crate::recognize::SweptShape::new(patch.surface.clone(), T);
        points
            .iter()
            .map(|p| shape.foot(*p).map_or(f64::INFINITY, |f| f.distance))
            .fold(0.0, f64::max)
    }

    /// A lopsided hill over a square sheet, which no canonical surface comes
    /// near, comes back as one patch charted by the mean-value map, every
    /// vertex within the distance of it.
    #[test]
    fn a_hill_is_charted_by_the_mean_value_map() {
        let (points, triangles) = sheet(40, |x, y| {
            let hill = (core::f64::consts::PI * x).sin() * (core::f64::consts::PI * y).sin();
            Point::new(x * 10.0, y * 10.0, 2.0 * hill * (1.0 + 0.3 * x))
        });
        let members: Vec<usize> = (0..triangles.len()).collect();
        let corners: HashSet<u32> = [0, 40, 1680, 1640].into_iter().collect();
        let region = Region {
            points: &points,
            triangles: &triangles,
            members: &members,
            corners: &corners,
        };
        let patch = fit_patch(&region, 1e-4, T).unwrap();
        assert!(patch.mapped);
        assert!(patch.deviation <= 1e-4, "{}", patch.deviation);
        assert!(worst(&patch, &points) <= 1e-4);
    }

    /// A quarter of a cylinder rippled by three times the distance is no
    /// cylinder at the distance, but one within ten times it, and lends the
    /// patch its chart: angle round and height along.
    #[test]
    fn a_rippled_cylinder_is_charted_by_the_cylinder() {
        let (points, triangles) = sheet(24, |x, y| {
            let (angle, z) = (x * core::f64::consts::FRAC_PI_2, y * 10.0);
            let r = 5.0 + 3e-4 * (3.0 * angle).sin() * (z * 0.5).cos();
            Point::new(r * angle.cos(), r * angle.sin(), z)
        });
        let members: Vec<usize> = (0..triangles.len()).collect();
        let corners = HashSet::new();
        let region = Region {
            points: &points,
            triangles: &triangles,
            members: &members,
            corners: &corners,
        };
        let patch = fit_patch(&region, 1e-4, T).unwrap();
        assert!(!patch.mapped);
        assert!(patch.deviation <= 1e-4, "{}", patch.deviation);
        assert!(worst(&patch, &points) <= 1e-4);
    }

    /// A ring of triangles round a hole is no disk, and gets no patch.
    #[test]
    fn a_ring_is_refused_as_no_disk() {
        let (points, mut triangles) = sheet(6, |x, y| Point::new(x, y, (x * y).sin()));
        // The middle cells taken out leave a hole.
        let hole: HashSet<usize> = [28, 29, 30, 31, 40, 41, 42, 43].into_iter().collect();
        triangles = triangles
            .into_iter()
            .enumerate()
            .filter(|(i, _)| !hole.contains(i))
            .map(|(_, t)| t)
            .collect();
        let members: Vec<usize> = (0..triangles.len()).collect();
        let corners = HashSet::new();
        let region = Region {
            points: &points,
            triangles: &triangles,
            members: &members,
            corners: &corners,
        };
        assert_eq!(fit_patch(&region, 1e-4, T).unwrap_err(), Refused::NotDisk);
    }
}
