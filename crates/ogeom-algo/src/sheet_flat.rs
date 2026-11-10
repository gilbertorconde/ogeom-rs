//! The least-distortion sheet: the conformal flattening of the geodesic
//! disc the sheet covers whose scale factor is one along its rim, scaled
//! to keep the disc's area.
//!
//! The disc is meshed in geodesic polar coordinates, rings of vertices about
//! the origin, each vertex landed on the surface exactly by the geodesic
//! sheet. The flat metric conformal to the surface's with no stretch along
//! the rim is unique (its scale's logarithm solves a Poisson problem held
//! to zero on the rim), and so is its map to the plane once the origin and
//! the sheet's `+x` there are fixed: on a sphere's geodesic disc it is the
//! stereographic projection from the antipode.
//!
//! The flattening is solved on `n` rings and on `2n`, and the most a
//! vertex the two share moves between them is taken as the estimated
//! departure of the finer one from the smooth map, which it bounds only
//! where the doubling at least halves that departure (measured near four
//! times on a sphere and a torus, not proved); the rings double until it
//! is within the resolution asked. The map from the flattened sheet back to geodesic polar
//! coordinates is then a polynomial fitted to the finer mesh's vertices,
//! its miss measured at every vertex and every triangle's centroid against
//! the mesh's own linear map there, raised in degree until within the
//! resolution as well. The sheet point is carried by the polynomial, then
//! landed by the geodesic sheet: smooth, and on the surface exactly.

use core::f64::consts::TAU;

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::{Point, Point2, SparseMatrix};

/// The most rings the finer flattening is meshed with.
const MOST_RINGS: usize = 64;

/// The highest total degree the polynomial back to geodesic coordinates
/// is raised to.
const MOST_DEGREE: usize = 16;

/// The flattened sheet's map back to geodesic polar coordinates.
#[derive(Debug, Clone)]
pub(crate) struct Flattening {
    degree: usize,
    /// The polynomial's coefficients for each coordinate, by monomial, in
    /// the flattened coordinates divided by `scale`.
    x: Vec<f64>,
    y: Vec<f64>,
    scale: f64,
    /// The estimated departure from the smooth conformal map, on the
    /// sheet: the mesh's and the polynomial's together. For the rigid
    /// flattening, the polynomial's miss of the mesh's own map alone.
    pub(crate) error: f64,
    /// The flattened mesh's vertices and triangle centroids within the
    /// reach, where the rigid flattening's stretch is measured.
    pub(crate) samples: Vec<Point2>,
}

/// How many rings the rigid flattening is meshed with.
const RIGID_RINGS: usize = 32;

/// The most local-global rounds the rigid flattening takes.
const RIGID_ROUNDS: usize = 400;

/// A polar mesh of the geodesic disc of radius `radius`, `rings` rings.
struct Mesh {
    /// Each vertex's geodesic polar position on the sheet.
    sheet: Vec<Point2>,
    /// Each vertex landed on the surface.
    space: Vec<Point>,
    triangles: Vec<[usize; 3]>,
    /// The vertices of the outer ring, in turn.
    outer: Vec<usize>,
}

impl Flattening {
    /// The flattening of the geodesic disc whose flattened reach covers
    /// `reach`, resolved to `resolution`. `land` puts a geodesic polar
    /// position on the surface; `most` is the largest geodesic radius the
    /// disc may take.
    pub(crate) fn build(
        land: &dyn Fn(Point2) -> OgeomResult<Point>,
        reach: f64,
        most: f64,
        resolution: f64,
        rigid: bool,
    ) -> OgeomResult<Self> {
        if !(resolution > 0.0 && resolution.is_finite()) {
            ogeom_bail!(
                Construction,
                "a flattening resolution of {resolution} is not a length"
            );
        }
        let mut radius = reach;
        // The flattened disc's rim stands nearer the origin than the
        // geodesic one wherever the surface curves away; the disc is
        // widened by the shortfall until its flattening covers the reach.
        loop {
            let (mesh, flat, departure) = if rigid {
                let mesh = mesh(land, radius, RIGID_RINGS)?;
                let flat = rigid_flattening(&mesh, flattened(&mesh)?)?;
                (mesh, flat, 0.0)
            } else {
                resolved(land, radius, resolution)?
            };
            let rim = rim_distance(&mesh, &flat);
            if rim >= reach {
                let method = if rigid {
                    "as-rigid-as-possible"
                } else {
                    "least-distortion"
                };
                let (degree, x, y, scale, miss) = fitted(&mesh, &flat, rim, resolution, method)?;
                let mut samples: Vec<Point2> = flat.clone();
                samples.extend(mesh.triangles.iter().map(|tri| {
                    let [a, b, c] = tri.map(|v| flat[v]);
                    Point2::new((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0)
                }));
                samples.retain(|p| p.x.hypot(p.y) <= reach);
                return Ok(Self {
                    degree,
                    x,
                    y,
                    scale,
                    error: departure + miss,
                    samples,
                });
            }
            let wider = radius * reach / rim;
            if wider >= most {
                ogeom_bail!(
                    Construction,
                    "the least-distortion sheet of reach {reach} needs a geodesic disc of \
                     radius {wider}, past the {most} the surface allows"
                );
            }
            // A ring's width more than the shortfall, so the next disc
            // passes the reach rather than closing on it.
            radius = wider + radius / f64::from(u32::try_from(MOST_RINGS).unwrap_or(1));
        }
    }

    /// The geodesic polar position a flattened sheet point is carried to,
    /// and the map's derivative there, rows by output coordinate.
    pub(crate) fn at(&self, p: Point2) -> (Point2, [[f64; 2]; 2]) {
        let (s, t) = (p.x / self.scale, p.y / self.scale);
        let mut q = [0.0_f64; 2];
        let mut d = [[0.0_f64; 2]; 2];
        for (k, (a, b)) in monomials(self.degree).enumerate() {
            let value = powi(s, a) * powi(t, b);
            let ds = if a == 0 {
                0.0
            } else {
                f64::from(a) * powi(s, a - 1) * powi(t, b)
            };
            let dt = if b == 0 {
                0.0
            } else {
                f64::from(b) * powi(s, a) * powi(t, b - 1)
            };
            for (row, c) in [&self.x, &self.y].into_iter().enumerate() {
                q[row] += c[k] * value;
                d[row][0] += c[k] * ds / self.scale;
                d[row][1] += c[k] * dt / self.scale;
            }
        }
        (Point2::new(q[0], q[1]), d)
    }
}

/// `x` to the non-negative power `n`.
fn powi(x: f64, n: u32) -> f64 {
    x.powi(i32::try_from(n).unwrap_or(i32::MAX))
}

/// The exponents `(a, b)` of every monomial `s^a t^b` of total degree up
/// to `degree`, in a fixed order.
fn monomials(degree: usize) -> impl Iterator<Item = (u32, u32)> {
    let top = u32::try_from(degree).unwrap_or(0);
    (0..=top).flat_map(move |total| (0..=total).map(move |b| (total - b, b)))
}

/// The flattening on rings doubled until the finer and the coarser agree
/// within `resolution` at every vertex they share: the finer mesh, its
/// flattening, and the most a shared vertex moved.
fn resolved(
    land: &dyn Fn(Point2) -> OgeomResult<Point>,
    radius: f64,
    resolution: f64,
) -> OgeomResult<(Mesh, Vec<Point2>, f64)> {
    let mut rings = 4;
    let mut coarse = {
        let mesh = mesh(land, radius, rings)?;
        let flat = flattened(&mesh)?;
        (mesh, flat)
    };
    loop {
        let mesh = mesh(land, radius, 2 * rings)?;
        let flat = flattened(&mesh)?;
        // Ring k of the coarser mesh is ring 2k of the finer, and its
        // j-th vertex the finer ring's 2j-th.
        let mut departure = flat[0].distance(coarse.1[0]);
        for k in 1..=rings {
            for j in 0..6 * k {
                let a = vertex(k, j);
                let b = vertex(2 * k, 2 * j);
                departure = departure.max(flat[b].distance(coarse.1[a]));
            }
        }
        if departure <= resolution / 2.0 {
            return Ok((mesh, flat, departure));
        }
        if 4 * rings > MOST_RINGS {
            ogeom_bail!(
                Construction,
                "the least-distortion sheet cannot be resolved to {resolution} on {} \
                 rings: its flattening still moves by {departure:.3e} between meshes",
                2 * rings
            );
        }
        rings *= 2;
        coarse = (mesh, flat);
    }
}

/// The index of vertex `j` of ring `k`, the origin vertex 0.
const fn vertex(k: usize, j: usize) -> usize {
    if k == 0 { 0 } else { 1 + 3 * k * (k - 1) + j }
}

/// The polar mesh of the geodesic disc of radius `radius` on `rings`
/// rings, ring `k` holding `6k` vertices from the sheet's `+x` round.
fn mesh(
    land: &dyn Fn(Point2) -> OgeomResult<Point>,
    radius: f64,
    rings: usize,
) -> OgeomResult<Mesh> {
    let ring_count = f64::from(u32::try_from(rings).unwrap_or(1));
    let mut sheet = vec![Point2::ORIGIN];
    for k in 1..=rings {
        let r = radius * f64::from(u32::try_from(k).unwrap_or(0)) / ring_count;
        let count = 6 * k;
        for j in 0..count {
            let angle = TAU * f64::from(u32::try_from(j).unwrap_or(0))
                / f64::from(u32::try_from(count).unwrap_or(1));
            sheet.push(Point2::new(r * angle.cos(), r * angle.sin()));
        }
    }
    let space = sheet
        .iter()
        .map(|q| land(*q))
        .collect::<OgeomResult<Vec<_>>>()?;
    let mut triangles = Vec::with_capacity(6 * rings * rings);
    for k in 1..=rings {
        // Zip ring k - 1 and ring k together by angle: each step takes the
        // ring whose next vertex turns least.
        let (inner, outer) = (6 * (k - 1), 6 * k);
        let turn = |count: usize, j: usize| {
            if count == 0 {
                0.0
            } else {
                f64::from(u32::try_from(j).unwrap_or(0))
                    / f64::from(u32::try_from(count).unwrap_or(1))
            }
        };
        let at_inner = |j: usize| {
            if inner == 0 {
                0
            } else {
                vertex(k - 1, j % inner)
            }
        };
        let (mut i, mut o) = (0, 0);
        while i < inner || o < outer {
            let take_outer = i >= inner || (o < outer && turn(outer, o + 1) <= turn(inner, i + 1));
            if take_outer {
                triangles.push([at_inner(i), vertex(k, o), vertex(k, (o + 1) % outer)]);
                o += 1;
            } else {
                triangles.push([at_inner(i), vertex(k, o % outer), at_inner(i + 1)]);
                i += 1;
            }
        }
    }
    let outer = (0..6 * rings).map(|j| vertex(rings, j)).collect();
    Ok(Mesh {
        sheet,
        space,
        triangles,
        outer,
    })
}

/// The conformal flattening of `mesh` whose scale factor is one on its
/// rim, the origin at the origin and the sheet's `+x` kept there, scaled
/// to keep the mesh's area.
///
/// The flat metric `e^{2u} g` of the surface's metric `g` has `Δu = K`,
/// `K` the Gaussian curvature, and here `u = 0` on the rim: on the mesh,
/// the cotangent Laplacian against each inner vertex's angle defect. The
/// rim keeps its lengths, and turns at each vertex by its own turn less
/// the Laplacian of `u` there, which by Gauss-Bonnet sum to a whole turn;
/// it is closed by the least change of its lengths, each weighed by
/// itself. The inner vertices are the harmonic extension of the rim, both
/// coordinates: a conformal map's coordinates are harmonic in the
/// surface's own conformal structure, which the cotangent Laplacian reads.
fn flattened(mesh: &Mesh) -> OgeomResult<Vec<Point2>> {
    let n = mesh.sheet.len();
    let mut on_rim = vec![false; n];
    for &v in &mesh.outer {
        on_rim[v] = true;
    }
    // The angle at each corner of each triangle, and each vertex's sum.
    let mut weight: Vec<ogeom_core::FastMap<usize, f64>> = vec![Default::default(); n];
    let mut angles = vec![0.0_f64; n];
    let mut space_area = 0.0;
    for tri in &mesh.triangles {
        let p = tri.map(|v| mesh.space[v]);
        let double = (p[1] - p[0]).cross(p[2] - p[0]).magnitude();
        if double.is_nan() || double <= 0.0 {
            ogeom_bail!(
                Construction,
                "the least-distortion sheet's mesh has a triangle of no area on the surface"
            );
        }
        space_area += double / 2.0;
        for k in 0..3 {
            let (i, j) = ((k + 1) % 3, (k + 2) % 3);
            let (a, b) = (p[i] - p[k], p[j] - p[k]);
            let angle = a.cross(b).magnitude().atan2(a.dot(b));
            angles[tri[k]] += angle;
            // The cotangent of the corner opposite edge (i, j), halved.
            let half = a.dot(b) / a.cross(b).magnitude() / 2.0;
            *weight[tri[i]].entry(tri[j]).or_insert(0.0) += half;
            *weight[tri[j]].entry(tri[i]).or_insert(0.0) += half;
        }
    }
    // `-L` over the inner vertices, which is positive definite: row i holds
    // the sum of i's weights on the diagonal and minus each weight to an
    // inner neighbour.
    let inner: Vec<usize> = (0..n).filter(|v| !on_rim[*v]).collect();
    let mut index = vec![usize::MAX; n];
    for (k, &v) in inner.iter().enumerate() {
        index[v] = k;
    }
    let mut triplets = Vec::new();
    for (k, &v) in inner.iter().enumerate() {
        let mut diagonal = 0.0;
        for (&w, &c) in &weight[v] {
            diagonal += c;
            if !on_rim[w] {
                triplets.push((k, index[w], -c));
            }
        }
        triplets.push((k, k, diagonal));
    }
    let system = SparseMatrix::from_triplets(inner.len(), inner.len(), &triplets);
    // The scale's logarithm: `-L u = -K` inside, zero on the rim.
    let defect: Vec<f64> = inner.iter().map(|&v| -(TAU - angles[v])).collect();
    let u_inner = conjugate_gradient(&system, &defect)?;
    let mut u = vec![0.0_f64; n];
    for (k, &v) in inner.iter().enumerate() {
        u[v] = u_inner[k];
    }
    // The rim laid flat: each vertex turns by its own turn on the surface
    // less the Laplacian of `u` there.
    let rim = &mesh.outer;
    let m = rim.len();
    let mut turns = Vec::with_capacity(m);
    for &v in rim {
        let laplacian: f64 = weight[v].iter().map(|(&w, &c)| c * (u[w] - u[v])).sum();
        turns.push(core::f64::consts::PI - angles[v] - laplacian);
    }
    let lengths: Vec<f64> = (0..m)
        .map(|k| mesh.space[rim[(k + 1) % m]].distance(mesh.space[rim[k]]))
        .collect();
    // Edge k runs from rim vertex k to k + 1, its heading the sum of the
    // turns at rim vertices 1 to k.
    let mut heading = 0.0;
    let mut tangents = Vec::with_capacity(m);
    for (k, turn) in turns.iter().enumerate() {
        if k > 0 {
            heading += turn;
        }
        tangents.push((heading.cos(), heading.sin()));
    }
    // The least change of lengths, each weighed by itself, that closes the
    // rim: `delta = -diag(l) T^t (T diag(l) T^t)^-1 sum(l t)`.
    let gap = (0..m).fold((0.0, 0.0), |g, k| {
        (
            lengths[k].mul_add(tangents[k].0, g.0),
            lengths[k].mul_add(tangents[k].1, g.1),
        )
    });
    let mut gram = [[0.0_f64; 2]; 2];
    for k in 0..m {
        let t = tangents[k];
        gram[0][0] += lengths[k] * t.0 * t.0;
        gram[0][1] += lengths[k] * t.0 * t.1;
        gram[1][1] += lengths[k] * t.1 * t.1;
    }
    let det = gram[0][0].mul_add(gram[1][1], -gram[0][1] * gram[0][1]);
    if det.is_nan() || det <= 0.0 {
        ogeom_bail!(
            Construction,
            "the least-distortion sheet's rim cannot be closed"
        );
    }
    let lambda = (
        gram[1][1].mul_add(gap.0, -gram[0][1] * gap.1) / det,
        gram[0][0].mul_add(gap.1, -gram[0][1] * gap.0) / det,
    );
    let mut flat = vec![Point2::ORIGIN; n];
    let mut at = Point2::ORIGIN;
    for k in 0..m {
        flat[rim[k]] = at;
        let t = tangents[k];
        let length = lengths[k] - lengths[k] * (t.0 * lambda.0 + t.1 * lambda.1);
        at = Point2::new(t.0.mul_add(length, at.x), t.1.mul_add(length, at.y));
    }
    // The inside: harmonic in both coordinates, the rim held.
    for axis in 0..2 {
        let mut rhs = vec![0.0_f64; inner.len()];
        for (k, &v) in inner.iter().enumerate() {
            for (&w, &c) in &weight[v] {
                if on_rim[w] {
                    rhs[k] += c * if axis == 0 { flat[w].x } else { flat[w].y };
                }
            }
        }
        let values = conjugate_gradient(&system, &rhs)?;
        for (k, &v) in inner.iter().enumerate() {
            if axis == 0 {
                flat[v].x = values[k];
            } else {
                flat[v].y = values[k];
            }
        }
    }
    let mut flat_area = 0.0;
    for tri in &mesh.triangles {
        let [a, b, c] = tri.map(|v| flat[v]);
        let signed = ((b - a).x * (c - a).y - (b - a).y * (c - a).x) / 2.0;
        if signed.is_nan() || signed <= 0.0 {
            ogeom_bail!(
                Construction,
                "the least-distortion sheet's flattening folds a triangle over"
            );
        }
        flat_area += signed;
    }
    // The origin to the origin, the first ring turned to its geodesic
    // directions on the mean (the sheet's `+x` at the origin: each ring
    // vertex's turn from its own direction is the map's turn there plus a
    // term in the vertex's position, which the ring's six vertices sum to
    // nothing), and scaled to keep the mesh's area.
    let centre = flat[0];
    let turn = (1..=6)
        .map(|v| {
            let (f, q) = (flat[v] - centre, mesh.sheet[v]);
            let d = f.y.atan2(f.x) - q.y.atan2(q.x);
            (d + core::f64::consts::PI).rem_euclid(TAU) - core::f64::consts::PI
        })
        .sum::<f64>()
        / 6.0;
    let scale = (space_area / flat_area).sqrt();
    let (sin, cos) = (-turn).sin_cos();
    for f in &mut flat {
        let d = *f - centre;
        *f = Point2::new(
            scale * cos.mul_add(d.x, -sin * d.y),
            scale * sin.mul_add(d.x, cos * d.y),
        );
    }
    Ok(flat)
}

/// The as-rigid-as-possible flattening of `mesh`, started from `start`:
/// the flat positions whose triangles are each nearest a rotation of the
/// triangle on the surface, every edge weighed by its cotangent weight.
///
/// Local and global steps alternate: each triangle's best rotation onto
/// its flat image (in the plane, the angle of the weighted sum of its
/// edges' cross and dot products), then the flat positions that best take
/// every rotated edge, a cotangent Laplacian solve with the origin held.
/// Neither step raises the energy, and the rounds stop where a round
/// lowers it by less than `EPSILON^(1/2)` of itself. The energy is not
/// convex and the positions found are the minimum the rounds settle in
/// from the conformal start, not a unique one: the sheet is judged by the
/// stretch it is measured to have, not by its distance from any smooth
/// map. Pinned as the conformal flattening is, the origin at the origin
/// and the first ring turned to its geodesic directions, and not scaled.
///
/// Refused where a triangle folds over, or the rim crosses itself, so the
/// mesh's map from the sheet is one to one where it is accepted.
#[allow(
    clippy::needless_range_loop,
    reason = "a triangle's corners and the edges opposite them are indexed together"
)]
fn rigid_flattening(mesh: &Mesh, start: Vec<Point2>) -> OgeomResult<Vec<Point2>> {
    let n = mesh.sheet.len();
    // Each triangle laid flat by itself, and its edges' half cotangent
    // weights: `own[t][k]` its corner `k`, `half[t][k]` the edge opposite.
    let mut own: Vec<[Point2; 3]> = Vec::with_capacity(mesh.triangles.len());
    let mut half: Vec<[f64; 3]> = Vec::with_capacity(mesh.triangles.len());
    for tri in &mesh.triangles {
        let p = tri.map(|v| mesh.space[v]);
        let (ab, ac) = (p[1] - p[0], p[2] - p[0]);
        let along = ab.magnitude();
        let across = ab.cross(ac).magnitude();
        if !(along > 0.0 && across > 0.0) || across.is_nan() {
            ogeom_bail!(
                Construction,
                "the as-rigid-as-possible sheet's mesh has a triangle of no area on the surface"
            );
        }
        own.push([
            Point2::ORIGIN,
            Point2::new(along, 0.0),
            Point2::new(ab.dot(ac) / along, across / along),
        ]);
        let mut h = [0.0_f64; 3];
        for (k, w) in h.iter_mut().enumerate() {
            let (a, b) = (p[(k + 1) % 3] - p[k], p[(k + 2) % 3] - p[k]);
            *w = a.dot(b) / a.cross(b).magnitude() / 2.0;
        }
        half.push(h);
    }
    // The cotangent Laplacian over every vertex but the origin, held at
    // the origin: positive definite, the Laplacian being semi-definite
    // with the constants its only null space.
    let mut weight: Vec<ogeom_core::FastMap<usize, f64>> = vec![Default::default(); n];
    for (tri, h) in mesh.triangles.iter().zip(&half) {
        for k in 0..3 {
            let (i, j) = (tri[(k + 1) % 3], tri[(k + 2) % 3]);
            *weight[i].entry(j).or_insert(0.0) += h[k];
            *weight[j].entry(i).or_insert(0.0) += h[k];
        }
    }
    let mut triplets = Vec::new();
    for v in 1..n {
        let mut diagonal = 0.0;
        for (&w, &c) in &weight[v] {
            diagonal += c;
            if w != 0 {
                triplets.push((v - 1, w - 1, -c));
            }
        }
        triplets.push((v - 1, v - 1, diagonal));
    }
    let system = SparseMatrix::from_triplets(n - 1, n - 1, &triplets);
    let mut flat = start;
    let origin = flat[0];
    for f in &mut flat {
        *f = Point2::new(f.x - origin.x, f.y - origin.y);
    }
    let energy = |flat: &[Point2], turns: &[(f64, f64)]| -> f64 {
        let mut e = 0.0;
        for (t, tri) in mesh.triangles.iter().enumerate() {
            let (cos, sin) = turns[t];
            for k in 0..3 {
                let (i, j) = ((k + 1) % 3, (k + 2) % 3);
                let l = own[t][i] - own[t][j];
                let f = flat[tri[i]] - flat[tri[j]];
                let r = (cos.mul_add(l.x, -sin * l.y), sin.mul_add(l.x, cos * l.y));
                e += half[t][k] * ((f.x - r.0).powi(2) + (f.y - r.1).powi(2));
            }
        }
        e
    };
    let best_turns = |flat: &[Point2]| -> Vec<(f64, f64)> {
        mesh.triangles
            .iter()
            .enumerate()
            .map(|(t, tri)| {
                let (mut cross, mut dot) = (0.0, 0.0);
                for k in 0..3 {
                    let (i, j) = ((k + 1) % 3, (k + 2) % 3);
                    let l = own[t][i] - own[t][j];
                    let f = flat[tri[i]] - flat[tri[j]];
                    cross += half[t][k] * l.x.mul_add(f.y, -l.y * f.x);
                    dot += half[t][k] * l.x.mul_add(f.x, l.y * f.y);
                }
                let angle = cross.atan2(dot);
                (angle.cos(), angle.sin())
            })
            .collect()
    };
    let mut turns = best_turns(&flat);
    let mut previous = energy(&flat, &turns);
    for _ in 0..RIGID_ROUNDS {
        // The global step: each vertex takes its edges' rotated lengths.
        let mut rhs = [vec![0.0_f64; n - 1], vec![0.0_f64; n - 1]];
        for (t, tri) in mesh.triangles.iter().enumerate() {
            let (cos, sin) = turns[t];
            for k in 0..3 {
                let (i, j) = ((k + 1) % 3, (k + 2) % 3);
                let l = own[t][i] - own[t][j];
                let r = (cos.mul_add(l.x, -sin * l.y), sin.mul_add(l.x, cos * l.y));
                let c = half[t][k];
                if tri[i] != 0 {
                    rhs[0][tri[i] - 1] += c * r.0;
                    rhs[1][tri[i] - 1] += c * r.1;
                }
                if tri[j] != 0 {
                    rhs[0][tri[j] - 1] -= c * r.0;
                    rhs[1][tri[j] - 1] -= c * r.1;
                }
            }
        }
        let xs = conjugate_gradient(&system, &rhs[0])?;
        let ys = conjugate_gradient(&system, &rhs[1])?;
        for v in 1..n {
            flat[v] = Point2::new(xs[v - 1], ys[v - 1]);
        }
        turns = best_turns(&flat);
        let now = energy(&flat, &turns);
        let settled = previous - now <= f64::EPSILON.sqrt() * previous;
        previous = now;
        if settled {
            break;
        }
    }
    for tri in &mesh.triangles {
        let [a, b, c] = tri.map(|v| flat[v]);
        let signed = ((b - a).x * (c - a).y - (b - a).y * (c - a).x) / 2.0;
        if signed.is_nan() || signed <= 0.0 {
            ogeom_bail!(
                Construction,
                "the as-rigid-as-possible sheet's flattening folds a triangle over"
            );
        }
    }
    if rim_crosses(&mesh.outer, &flat) {
        ogeom_bail!(
            Construction,
            "the as-rigid-as-possible sheet's flattened rim crosses itself"
        );
    }
    // The first ring turned to its geodesic directions on the mean, as the
    // conformal flattening is.
    let turn = (1..=6)
        .map(|v| {
            let (f, q) = (flat[v], mesh.sheet[v]);
            let d = f.y.atan2(f.x) - q.y.atan2(q.x);
            (d + core::f64::consts::PI).rem_euclid(TAU) - core::f64::consts::PI
        })
        .sum::<f64>()
        / 6.0;
    let (sin, cos) = (-turn).sin_cos();
    for f in &mut flat {
        *f = Point2::new(cos.mul_add(f.x, -sin * f.y), sin.mul_add(f.x, cos * f.y));
    }
    Ok(flat)
}

/// Whether the flattened rim, the closed path through `rim` in turn,
/// crosses itself: two edges that share no vertex meet.
fn rim_crosses(rim: &[usize], flat: &[Point2]) -> bool {
    let m = rim.len();
    let orient = |a: Point2, b: Point2, c: Point2| {
        let (u, v) = (b - a, c - a);
        u.x.mul_add(v.y, -u.y * v.x)
    };
    for i in 0..m {
        let (a, b) = (flat[rim[i]], flat[rim[(i + 1) % m]]);
        for j in i + 2..m {
            if (j + 1) % m == i {
                continue;
            }
            let (c, d) = (flat[rim[j]], flat[rim[(j + 1) % m]]);
            let (o1, o2) = (orient(a, b, c), orient(a, b, d));
            let (o3, o4) = (orient(c, d, a), orient(c, d, b));
            if o1 * o2 <= 0.0 && o3 * o4 <= 0.0 {
                return true;
            }
        }
    }
    false
}

/// The solution of `a x = b`, `a` symmetric positive definite, by
/// conjugate gradients from zero, until the residual stands within
/// `EPSILON^(3/4)` of `b`, which leaves the solution's own error far under
/// the arithmetic's square root and its accumulated rounding over the
/// steps clear of it; in exact arithmetic it is reached in as many steps
/// as there are unknowns, and twice that is allowed.
fn conjugate_gradient(a: &SparseMatrix, b: &[f64]) -> OgeomResult<Vec<f64>> {
    let dot = |x: &[f64], y: &[f64]| x.iter().zip(y).fold(0.0, |s, (p, q)| p.mul_add(*q, s));
    let n = b.len();
    let target = f64::EPSILON.powf(0.75) * dot(b, b).sqrt();
    let mut x = vec![0.0; n];
    let mut r = b.to_vec();
    let mut p = r.clone();
    let mut rr = dot(&r, &r);
    if rr.sqrt() <= target {
        return Ok(x);
    }
    for _ in 0..2 * n.max(1) {
        let ap = a.multiply(&p);
        let pap = dot(&p, &ap);
        if pap.is_nan() || pap <= 0.0 {
            ogeom_bail!(
                Construction,
                "the least-distortion sheet's Laplacian is not positive definite"
            );
        }
        let alpha = rr / pap;
        for k in 0..n {
            x[k] = alpha.mul_add(p[k], x[k]);
            r[k] = alpha.mul_add(-ap[k], r[k]);
        }
        let next = dot(&r, &r);
        if next.sqrt() <= target {
            return Ok(x);
        }
        let beta = next / rr;
        rr = next;
        for k in 0..n {
            p[k] = beta.mul_add(p[k], r[k]);
        }
    }
    ogeom_bail!(
        Construction,
        "the least-distortion sheet's solve did not converge"
    )
}

/// The least distance from the origin to the flattened mesh's rim.
fn rim_distance(mesh: &Mesh, flat: &[Point2]) -> f64 {
    let mut least = f64::INFINITY;
    for (k, &v) in mesh.outer.iter().enumerate() {
        let (a, b) = (flat[v], flat[mesh.outer[(k + 1) % mesh.outer.len()]]);
        let along = b - a;
        let f = (-(a.x * along.x + a.y * along.y) / along.dot(along)).clamp(0.0, 1.0);
        let foot = Point2::new(along.x.mul_add(f, a.x), along.y.mul_add(f, a.y));
        least = least.min(foot.x.hypot(foot.y));
    }
    least
}

/// The polynomial from the flattened sheet back to geodesic polar
/// coordinates, its degree raised until its miss at every vertex and
/// every triangle's centroid is within half `resolution`: the degree, the
/// coefficients of each coordinate, the scale they are read in, and the
/// miss. A refusal names the sheet by `method`.
type Fitted = (usize, Vec<f64>, Vec<f64>, f64, f64);

fn fitted(
    mesh: &Mesh,
    flat: &[Point2],
    rim: f64,
    resolution: f64,
    method: &str,
) -> OgeomResult<Fitted> {
    // The samples: every vertex, and every triangle's centroid with the
    // mesh's own linear map there.
    let mut from: Vec<Point2> = flat.to_vec();
    let mut to: Vec<Point2> = mesh.sheet.clone();
    for tri in &mesh.triangles {
        let mean = |points: &[Point2]| {
            let [a, b, c] = tri.map(|v| points[v]);
            Point2::new((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0)
        };
        from.push(mean(flat));
        to.push(mean(&mesh.sheet));
    }
    let scale = rim;
    let mut best: Option<Fitted> = None;
    for degree in 1..=MOST_DEGREE {
        let terms: Vec<(u32, u32)> = monomials(degree).collect();
        let basis = nalgebra::DMatrix::from_fn(from.len(), terms.len(), |i, k| {
            let (a, b) = terms[k];
            powi(from[i].x / scale, a) * powi(from[i].y / scale, b)
        });
        let svd = basis.svd(true, true);
        let solve = |values: Vec<f64>| -> OgeomResult<Vec<f64>> {
            let column = nalgebra::DVector::from_vec(values);
            match svd.solve(&column, f64::EPSILON) {
                Ok(c) => Ok(c.iter().copied().collect()),
                Err(e) => ogeom_bail!(Construction, "the flattening's polynomial fit failed: {e}"),
            }
        };
        let x = solve(to.iter().map(|q| q.x).collect())?;
        let y = solve(to.iter().map(|q| q.y).collect())?;
        let candidate = Flattening {
            degree,
            x,
            y,
            scale,
            error: 0.0,
            samples: Vec::new(),
        };
        let mut miss: f64 = 0.0;
        for (f, q) in from.iter().zip(&to) {
            miss = miss.max(candidate.at(*f).0.distance(*q));
        }
        let folds = from.iter().any(|f| {
            let d = candidate.at(*f).1;
            let det = d[0][0].mul_add(d[1][1], -d[0][1] * d[1][0]);
            det.is_nan() || det <= 0.0
        });
        if !folds && best.as_ref().is_none_or(|b| miss < b.4) {
            best = Some((degree, candidate.x, candidate.y, scale, miss));
        }
        if !folds && miss <= resolution / 2.0 {
            break;
        }
    }
    match best {
        Some(b) if b.4 <= resolution / 2.0 => Ok(b),
        Some(b) => ogeom_bail!(
            Construction,
            "the {method} sheet's map cannot be held within {:.3e}, half its resolution of \
             {resolution}: its best polynomial misses the flattening by {:.3e}",
            resolution / 2.0,
            b.4
        ),
        None => ogeom_bail!(
            Construction,
            "the {method} sheet's map folds over at every degree tried"
        ),
    }
}
