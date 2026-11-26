//! The zones with exact solutions: the narrowest slab, the narrowest
//! ring, the smallest enclosing circle and the largest empty circle.
//!
//! The slab and the ring are found over a working subset of the points,
//! seeded from the least-squares fit's extremes. A subset's narrowest zone
//! is found exactly by trying every configuration the optimum can take,
//! and is no wider than the whole set's, since fewer points are easier to
//! hold. If it holds every point, to the rounding of their offsets, it is
//! the whole set's optimum; if not, the points furthest outside it join
//! the subset, and the subset's zone is found again. Each round adds a
//! point the subset lacked, so the rounds end, and the subset stays small:
//! about as many points as the zone has contacts, and those it took to
//! find them.

use super::hull;
use super::lp::ROUNDING;
use ogeom_core::predicates::{Exact, Predicates, Sign};
use ogeom_math::{Point, Point2, Vector, Vector2};
use spade::handles::FixedVertexHandle;
use spade::{DelaunayTriangulation, Point2 as SpadePoint, Triangulation as _};

/// A slab between two parallel planes.
#[derive(Debug, Clone, Copy)]
pub(super) struct Slab {
    /// The planes' unit normal.
    pub normal: Vector,
    /// The planes' offsets along it, from the origin the points were
    /// measured from.
    pub low: f64,
    pub high: f64,
}

/// The narrowest slab holding every point, its offsets measured from
/// `origin`; `seeds` index points that span a plane.
pub(super) fn narrowest_slab(points: &[Point], origin: Point, seeds: &[usize]) -> Option<Slab> {
    let local: Vec<Vector> = points.iter().map(|p| *p - origin).collect();
    let mut subset = distinct(seeds);
    loop {
        let few: Vec<Vector> = subset.iter().map(|&i| local[i]).collect();
        let mut slab = match hull::hull(&few) {
            Some(faces) => slab_of(&few, &faces)?,
            None => {
                // The subset lies in one plane: the whole set does too, and
                // is flat, or a point off that plane joins it.
                let [a, b, c] = spread(&few)?;
                let off = local.iter().position(|v| {
                    Exact::orient3d(
                        few[a].to_array(),
                        few[b].to_array(),
                        few[c].to_array(),
                        v.to_array(),
                    ) != Sign::Zero
                });
                if let Some(i) = off {
                    subset.push(i);
                    continue;
                }
                let n = unit((few[b] - few[a]).cross(few[c] - few[a]))?;
                Slab {
                    normal: n,
                    low: 0.0,
                    high: 0.0,
                }
            }
        };
        (slab.low, slab.high) = extent(&few, slab.normal);
        // A point outside the subset's slab by more than the rounding of
        // its offset is one the subset lacks, the subset's own offsets
        // being measured alike.
        let past = |side: f64| {
            let mut worst: Option<(f64, usize)> = None;
            for (i, v) in local.iter().enumerate() {
                let bound = if side > 0.0 { slab.high } else { slab.low };
                let excess = side * (slab.normal.dot(*v) - bound) - ROUNDING * v.magnitude();
                if excess > 0.0 && worst.is_none_or(|(w, _)| excess > w) {
                    worst = Some((excess, i));
                }
            }
            worst.map(|(_, i)| i)
        };
        let outside: Vec<usize> = [past(1.0), past(-1.0)].into_iter().flatten().collect();
        if outside.is_empty() {
            (slab.low, slab.high) = extent(&local, slab.normal);
            return Some(slab);
        }
        subset.extend(outside);
    }
}

/// The least and greatest offset of `vectors` along `n`.
fn extent(vectors: &[Vector], n: Vector) -> (f64, f64) {
    vectors
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| {
            let h = n.dot(*v);
            (lo.min(h), hi.max(h))
        })
}

/// The narrowest slab of a few points. Its planes touch the points' hull
/// face to vertex, one plane on a face and the other through the vertex
/// farthest from it, or edge to edge, along a normal square to both edges
/// that each edge's two faces bracket, one edge's on the normal's side and
/// the other's on the opposite (an antipodal pair: the edges' arcs on the
/// sphere of normals cross). Every face's normal is tried, and every such
/// pair's.
fn slab_of(few: &[Vector], faces: &[[usize; 3]]) -> Option<Slab> {
    let mut corners: Vec<usize> = faces.iter().flatten().copied().collect();
    corners.sort_unstable();
    corners.dedup();

    let mut best: Option<Slab> = None;
    let mut consider = |n: Vector| {
        let bound = best.map_or(f64::INFINITY, |b| b.high - b.low);
        let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
        for &i in &corners {
            let h = n.dot(few[i]);
            low = low.min(h);
            high = high.max(h);
            if high - low >= bound {
                return;
            }
        }
        best = Some(Slab {
            normal: n,
            low,
            high,
        });
    };
    let normals: Vec<Option<Vector>> = faces
        .iter()
        .map(|&[a, b, c]| unit((few[b] - few[a]).cross(few[c] - few[a])))
        .collect();
    for n in normals.iter().flatten() {
        consider(*n);
    }

    // Each edge once, with the normals of the faces either side.
    let mut edges: Vec<(usize, usize, Vector, Vector)> = Vec::new();
    for (f, face) in faces.iter().enumerate() {
        for k in 0..3 {
            let (u, v) = (face[k], face[(k + 1) % 3]);
            if u > v {
                continue;
            }
            let twin = faces
                .iter()
                .position(|g| (0..3).any(|j| g[j] == v && g[(j + 1) % 3] == u));
            if let (Some(n), Some(m)) = (normals[f], twin.and_then(|g| normals[g])) {
                edges.push((u, v, n, m));
            }
        }
    }
    for (a, &(u, v, n1, n2)) in edges.iter().enumerate() {
        let e = few[v] - few[u];
        for &(p, q, m1, m2) in &edges[a + 1..] {
            let Some(d) = unit(e.cross(few[q] - few[p])) else {
                continue;
            };
            if (between(d, n1, n2) && between(-d, m1, m2))
                || (between(-d, n1, n2) && between(d, m1, m2))
            {
                consider(d);
            }
        }
    }
    best
}

/// Whether unit `d`, square to the edge between faces of unit normals `n`
/// and `m`, lies on the arc between them: `d = αn + βm` with `α, β ≥ 0`.
/// Faces in one plane bound no arc; their normal is a face's, tried as
/// such.
fn between(d: Vector, n: Vector, m: Vector) -> bool {
    let c = n.dot(m);
    let gap = 1.0 - c * c;
    if gap <= 0.0 {
        return false;
    }
    let (a, b) = (d.dot(n), d.dot(m));
    a - c * b >= 0.0 && b - c * a >= 0.0
}

fn unit(v: Vector) -> Option<Vector> {
    let m = v.magnitude();
    (m > 0.0 && m.is_finite()).then(|| v / m)
}

/// Three points spanning a triangle, the first among them; `None` where
/// they all lie on one line.
fn spread(points: &[Vector]) -> Option<[usize; 3]> {
    let b = (1..points.len()).find(|&i| points[i] != points[0])?;
    let e = points[b] - points[0];
    let c = (0..points.len()).max_by(|&i, &j| {
        let s = |k: usize| e.cross(points[k] - points[0]).magnitude();
        s(i).total_cmp(&s(j))
    })?;
    (e.cross(points[c] - points[0]).magnitude() > 0.0).then_some([0, b, c])
}

/// What the narrowest ring about a few points came to.
#[derive(Debug, Clone)]
pub(super) enum Ring {
    /// Two concentric circles about `centre`: each point's distance from
    /// it less the distance of a reference point, `radius`, and the least
    /// and greatest of them.
    Circles {
        centre: Point2,
        radius: f64,
        deviations: Vec<f64>,
        inner: f64,
        outer: f64,
    },
    /// No ring is narrower than the strip between two parallel lines that
    /// holds every point: a ring is a strip's limit as its centre goes off
    /// to infinity, and the narrowest ring is not reached.
    Straight { width: f64 },
}

/// The narrowest ring holding every point; `seeds` index points that span
/// the plane, and `origin` is a point near the centre the coordinates are
/// taken from.
pub(super) fn narrowest_ring(points: &[Point2], origin: Point2, seeds: &[usize]) -> Ring {
    let local: Vec<Vector2> = points.iter().map(|p| *p - origin).collect();
    let mut subset = distinct(seeds);
    loop {
        let (ring, strip) = ring_of(&local, &subset);
        let ring_width = ring.map_or(f64::INFINITY, |(_, w)| w);
        // Each point's offset across the zone, and the rounding it carries;
        // the subset's extent is measured alike, so a point outside it by
        // more than that rounding is one the subset lacks.
        let (offsets, rounding, inner, outer): (Vec<f64>, Vec<f64>, f64, f64) = match strip {
            Some(Strip { normal, low, high }) if high - low < ring_width => (
                local.iter().map(|v| normal.dot(*v)).collect(),
                local.iter().map(|v| ROUNDING * v.magnitude()).collect(),
                low,
                high,
            ),
            _ => {
                let Some((centre, _)) = ring else {
                    // No three points span a circle and none a strip: the
                    // subset is one point, which the seeds never are.
                    return Ring::Straight { width: 0.0 };
                };
                let reference = subset[0];
                let offsets: Vec<f64> = (0..local.len())
                    .map(|i| beyond(&local, centre, reference, i))
                    .collect();
                let (inner, outer) = subset
                    .iter()
                    .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &i| {
                        (lo.min(offsets[i]), hi.max(offsets[i]))
                    });
                let size = local[reference].magnitude();
                let rounding = local
                    .iter()
                    .map(|v| ROUNDING * (v.magnitude() + size))
                    .collect();
                (offsets, rounding, inner, outer)
            }
        };
        let past = |side: f64| {
            let bound = if side > 0.0 { outer } else { inner };
            let mut worst: Option<(f64, usize)> = None;
            for (i, h) in offsets.iter().enumerate() {
                let excess = side * (h - bound) - rounding[i];
                if excess > 0.0 && worst.is_none_or(|(w, _)| excess > w) {
                    worst = Some((excess, i));
                }
            }
            worst.map(|(_, i)| i)
        };
        let outside: Vec<usize> = [past(1.0), past(-1.0)].into_iter().flatten().collect();
        if !outside.is_empty() {
            subset.extend(outside);
            continue;
        }
        let (least, most) = offsets
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), h| {
                (lo.min(*h), hi.max(*h))
            });
        return match ring {
            Some((centre, _)) if strip.is_none_or(|s| s.high - s.low >= ring_width) => {
                Ring::Circles {
                    centre: Point2::from_vector(origin.to_vector() + centre),
                    radius: (local[subset[0]] - centre).magnitude(),
                    deviations: offsets,
                    inner: least,
                    outer: most,
                }
            }
            _ => Ring::Straight {
                width: most - least,
            },
        };
    }
}

/// How much farther point `i` is from `c` than point `reference`, as
/// `(|a|² − |b|²) / (|a| + |b|)`: exact to rounding of the difference
/// itself, where subtracting the two distances loses it to their size
/// when `c` is far off.
fn beyond(local: &[Vector2], c: Vector2, reference: usize, i: usize) -> f64 {
    let (a, b) = (local[i] - c, local[reference] - c);
    let sum = a.magnitude() + b.magnitude();
    if sum == 0.0 {
        return 0.0;
    }
    (local[i] - local[reference]).dot(a + b) / sum
}

/// A strip between two parallel lines, by their unit normal and offsets.
#[derive(Debug, Clone, Copy)]
struct Strip {
    normal: Vector2,
    low: f64,
    high: f64,
}

/// The narrowest ring of a few points at a finite centre, by its centre
/// and width, and their narrowest strip, by its unit normal and offsets.
///
/// The ring's centre lies where the farthest-point and nearest-point
/// Voronoi diagrams cross: at a vertex of either, the centre of a circle
/// through three points, or where an edge of one crosses an edge of the
/// other, equidistant from two points and from two others. Within a cell
/// of the overlay, and along an edge, the width has no minimum, since a
/// unit vector cannot lie strictly between two others; every vertex is
/// tried. The strip's sides touch along an edge of the points' hull: every
/// pair's direction is tried.
fn ring_of(local: &[Vector2], subset: &[usize]) -> (Option<(Vector2, f64)>, Option<Strip>) {
    let mut best: Option<(Vector2, f64)> = None;
    let mut consider = |c: Vector2| {
        if !(c.x.is_finite() && c.y.is_finite()) {
            return;
        }
        let bound = best.map_or(f64::INFINITY, |(_, w)| w);
        let (mut inner, mut outer) = (f64::INFINITY, f64::NEG_INFINITY);
        for &i in subset {
            let r = beyond(local, c, subset[0], i);
            inner = inner.min(r);
            outer = outer.max(r);
            if outer - inner >= bound {
                return;
            }
        }
        best = Some((c, outer - inner));
    };
    let edges = pairs(subset);
    for (a, &(i, j)) in edges.iter().enumerate() {
        for &k in subset.iter().filter(|&&k| k > j) {
            if let Some(c) = circumcentre(local[i], local[j], local[k]) {
                consider(c);
            }
        }
        for &(k, l) in &edges[a + 1..] {
            if let Some(c) = bisectors_meet(local[i], local[j], local[k], local[l]) {
                consider(c);
            }
        }
    }

    let mut strip: Option<Strip> = None;
    for &(i, j) in &edges {
        let e = local[j] - local[i];
        let m = e.magnitude();
        if m == 0.0 {
            continue;
        }
        let normal = e.perpendicular() / m;
        let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
        for &k in subset {
            let h = normal.dot(local[k]);
            low = low.min(h);
            high = high.max(h);
        }
        if strip.is_none_or(|s| high - low < s.high - s.low) {
            strip = Some(Strip { normal, low, high });
        }
    }
    (best, strip)
}

/// The centre of the circle through three points; `None` where they are
/// collinear.
pub(super) fn circumcentre(a: Vector2, b: Vector2, c: Vector2) -> Option<Vector2> {
    let (u, v) = (b - a, c - a);
    let d = 2.0 * u.cross(v);
    if d == 0.0 {
        return None;
    }
    let (uu, vv) = (u.square_magnitude(), v.square_magnitude());
    Some(a + Vector2::new(v.y * uu - u.y * vv, u.x * vv - v.x * uu) / d)
}

/// The point equidistant from `a` and `b` and from `c` and `d`; `None`
/// where the two bisectors are parallel.
fn bisectors_meet(a: Vector2, b: Vector2, c: Vector2, d: Vector2) -> Option<Vector2> {
    let (u, v) = (b - a, d - c);
    let det = u.cross(v);
    if det == 0.0 {
        return None;
    }
    // u·x = u·(a + b)/2 and v·x = v·(c + d)/2.
    let (p, q) = (u.dot(a + b) / 2.0, v.dot(c + d) / 2.0);
    Some(Vector2::new(p * v.y - q * u.y, u.x * q - v.x * p) / det)
}

fn pairs(subset: &[usize]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for (a, &i) in subset.iter().enumerate() {
        for &j in &subset[a + 1..] {
            out.push((i.min(j), i.max(j)));
        }
    }
    out
}

fn distinct(seeds: &[usize]) -> Vec<usize> {
    let mut out = seeds.to_vec();
    out.sort_unstable();
    out.dedup();
    out
}

/// The smallest circle enclosing the points, by Welzl's algorithm over a
/// fixed shuffle: its centre and radius, the radius the greatest distance
/// from that centre.
pub(super) fn enclosing_circle(points: &[Point2]) -> (Point2, f64) {
    let origin = points[0];
    let mut order: Vec<Vector2> = points.iter().map(|p| *p - origin).collect();
    let mut draws = super::super::primitive::Draws::new();
    for i in (1..order.len()).rev() {
        order.swap(i, draws.below(i + 1));
    }
    let inside = |c: Vector2, r2: f64, p: Vector2| (p - c).square_magnitude() <= r2;
    let mut c = order[0];
    let mut r2 = 0.0;
    for i in 1..order.len() {
        if inside(c, r2, order[i]) {
            continue;
        }
        c = order[i];
        r2 = 0.0;
        for j in 0..i {
            if inside(c, r2, order[j]) {
                continue;
            }
            c = (order[i] + order[j]) / 2.0;
            r2 = (order[i] - c).square_magnitude();
            for k in 0..j {
                if inside(c, r2, order[k]) {
                    continue;
                }
                (c, r2) = through_three(order[i], order[j], order[k]);
            }
        }
    }
    let radius = order
        .iter()
        .fold(0.0_f64, |m, p| m.max((*p - c).magnitude()));
    (Point2::from_vector(origin.to_vector() + c), radius)
}

/// The circle through three points, or across the farthest two where
/// rounding has left them collinear.
fn through_three(a: Vector2, b: Vector2, c: Vector2) -> (Vector2, f64) {
    if let Some(centre) = circumcentre(a, b, c) {
        let r2 = [a, b, c]
            .iter()
            .fold(0.0_f64, |m, p| m.max((*p - centre).square_magnitude()));
        return (centre, r2);
    }
    let [(p, q), ..] = {
        let mut spans = [(a, b), (a, c), (b, c)];
        spans.sort_by(|x, y| {
            (y.0 - y.1)
                .square_magnitude()
                .total_cmp(&(x.0 - x.1).square_magnitude())
        });
        spans
    };
    let centre = (p + q) / 2.0;
    (centre, (p - centre).square_magnitude())
}

/// The largest circle whose centre lies in the points' convex hull and
/// which holds none of them inside: its centre, or `None` where the points
/// lie on one line and enclose nothing.
///
/// Its centre is a vertex of the points' Voronoi diagram inside the hull,
/// or where an edge of the diagram crosses the hull's boundary: within a
/// Voronoi cell the distance to its point has no maximum, and along an
/// edge the distance to its two points is convex.
pub(super) fn largest_empty_circle(points: &[Point2]) -> Option<Point2> {
    let origin = points[0];
    let local: Vec<SpadePoint<f64>> = points
        .iter()
        .map(|p| {
            let v = *p - origin;
            SpadePoint::new(v.x, v.y)
        })
        .collect();
    let dt: DelaunayTriangulation<SpadePoint<f64>> =
        DelaunayTriangulation::bulk_load(local).ok()?;
    if dt.num_inner_faces() == 0 {
        return None;
    }
    let at = |p: SpadePoint<f64>| Vector2::new(p.x, p.y);
    // The hull, counter-clockwise.
    let mut hull: Vec<Vector2> = dt.convex_hull().map(|e| at(e.from().position())).collect();
    let twice_area: f64 = (0..hull.len())
        .map(|k| hull[k].cross(hull[(k + 1) % hull.len()]))
        .sum();
    if twice_area < 0.0 {
        hull.reverse();
    }

    // The distance from `c` to its nearest point, by a greedy walk over
    // the triangulation from `from`: a point that is not the nearest has a
    // Delaunay neighbour nearer.
    let nearest = |from: FixedVertexHandle, c: Vector2| {
        let mut v = dt.vertex(from);
        let mut d = (at(v.position()) - c).magnitude();
        loop {
            let closer = v
                .out_edges()
                .map(|e| e.to())
                .map(|w| ((at(w.position()) - c).magnitude(), w))
                .filter(|(dw, _)| *dw < d)
                .min_by(|a, b| a.0.total_cmp(&b.0));
            match closer {
                Some((dw, w)) => {
                    d = dw;
                    v = w;
                }
                None => return d,
            }
        }
    };
    // Each candidate is weighed by its true distance to the nearest point,
    // so a centre that rounding has moved off its Voronoi edge is still
    // weighed truly. A vertex counts inside the hull only; a crossing of the
    // hull's boundary lies on it to rounding, and counts as on it.
    let mut best: Option<(Vector2, f64)> = None;
    let mut consider = |c: Vector2, crossing: bool, from: FixedVertexHandle| {
        if !(c.x.is_finite() && c.y.is_finite()) || !(crossing || in_hull(&hull, c)) {
            return;
        }
        let r = nearest(from, c);
        if best.is_none_or(|(_, b)| r > b) {
            best = Some((c, r));
        }
    };
    let centre_of = |face: spade::handles::FaceHandle<'_, spade::handles::InnerTag, _, _, _, _>| {
        let [a, b, c] = face.vertices().map(|v| at(v.position()));
        circumcentre(a, b, c)
    };
    for face in dt.inner_faces() {
        if let Some(centre) = centre_of(face) {
            consider(centre, false, face.vertices()[0].fix());
        }
    }
    // A Voronoi edge lies on its two points' bisector, between the
    // circumcentres of the triangles either side of their Delaunay edge,
    // or from the one triangle's outward where the edge is on the hull. It
    // is measured along the bisector from the two points' midpoint, so a
    // circumcentre that a near-collinear triangle throws far off moves only
    // the edge's end, not where it crosses the hull.
    for edge in dt.undirected_edges() {
        let directed = edge.as_directed();
        let site = at(directed.from().position());
        let other = at(directed.to().position());
        let middle = (site + other) / 2.0;
        let e = other - site;
        let across = e.perpendicular() / e.magnitude();
        let faces = [directed.face(), directed.rev().face()].map(|f| f.as_inner());
        let ends = faces.map(|f| f.and_then(centre_of));
        let (low, high) = match (faces, ends) {
            (_, [Some(a), Some(b)]) => {
                if in_hull(&hull, a) && in_hull(&hull, b) {
                    continue;
                }
                let (s, t) = (across.dot(a - middle), across.dot(b - middle));
                (s.min(t), s.max(t))
            }
            ([Some(face), _] | [_, Some(face)], [Some(a), _] | [_, Some(a)]) => {
                // A ray from inside the hull leaves it through the hull
                // edge it bisects, at that edge's midpoint, nearer both
                // points than its start, a vertex already weighed.
                if in_hull(&hull, a) {
                    continue;
                }
                let s = across.dot(a - middle);
                let [p, q, r] = face.vertices().map(|v| at(v.position()));
                let inward = across.dot((p + q + r) / 3.0 - middle) > 0.0;
                if inward {
                    (f64::NEG_INFINITY, s)
                } else {
                    (s, f64::INFINITY)
                }
            }
            _ => continue,
        };
        if let Some((s, t)) = clip(&hull, middle, across, low, high) {
            let from = directed.from().fix();
            consider(middle + across * s, s > low, from);
            consider(middle + across * t, t < high, from);
        }
    }
    let (centre, _) = best?;
    Some(Point2::from_vector(origin.to_vector() + centre))
}

/// Whether `q` lies in the convex polygon `hull`, counter-clockwise, its
/// boundary included.
fn in_hull(hull: &[Vector2], q: Vector2) -> bool {
    let o2 =
        |a: Vector2, b: Vector2, c: Vector2| Exact::orient2d([a.x, a.y], [b.x, b.y], [c.x, c.y]);
    let h = hull.len();
    if h < 3 {
        return false;
    }
    let o = hull[0];
    if o2(o, hull[1], q) == Sign::Negative || o2(o, hull[h - 1], q) == Sign::Positive {
        return false;
    }
    let (mut lo, mut hi) = (1, h - 1);
    while hi - lo > 1 {
        let mid = usize::midpoint(lo, hi);
        if o2(o, hull[mid], q) == Sign::Negative {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    o2(hull[lo], hull[hi], q) != Sign::Negative
}

/// The part of the line `from + t·direction`, `low ≤ t ≤ high`, inside
/// the convex polygon `hull`, counter-clockwise, by its parameter range;
/// `None` where it misses the polygon.
fn clip(
    hull: &[Vector2],
    from: Vector2,
    direction: Vector2,
    low: f64,
    high: f64,
) -> Option<(f64, f64)> {
    let (mut t0, mut t1) = (low, high);
    for (k, &a) in hull.iter().enumerate() {
        let b = hull[(k + 1) % hull.len()];
        // Inside is to the left of a → b: (b − a) × (p − a) ≥ 0.
        let e = b - a;
        let at_start = e.cross(from - a);
        let rate = e.cross(direction);
        if rate == 0.0 {
            if at_start < 0.0 {
                return None;
            }
            continue;
        }
        let t = -at_start / rate;
        if rate > 0.0 {
            t0 = t0.max(t);
        } else {
            t1 = t1.min(t);
        }
        if t0 > t1 {
            return None;
        }
    }
    (t0.is_finite() && t1.is_finite()).then_some((t0, t1))
}
