//! The least-distortion sheet: a least-squares conformal flattening of the
//! geodesic disc the sheet covers, scaled to keep its area.
//!
//! The disc is meshed in geodesic polar coordinates, rings of vertices about
//! the origin, each vertex landed on the surface exactly by the geodesic
//! sheet. The mesh is flattened by least-squares conformal maps: each
//! triangle's map to the plane is asked to be conformal, `sum f_j d_j = 0`
//! for its flattened corners `f_j` and opposite edges `d_j` read as complex
//! numbers in the triangle's own plane, weighted by the root of its area,
//! the origin pinned and its first neighbour on the sheet's `+x`; the
//! sparse system is solved by conjugate gradients on its normal equations.
//! The flattening is scaled so the mesh keeps its area on the plane.
//!
//! The flattening is solved on `n` rings and on `2n`, and the most a
//! vertex the two share moves between them is the stated departure from
//! the smooth map; the rings double until it is within the resolution
//! asked. The map from the flattened sheet back to geodesic polar
//! coordinates is then a polynomial fitted to the finer mesh's vertices,
//! its miss measured at every vertex and every triangle's centroid against
//! the mesh's own linear map there, raised in degree until within the
//! resolution as well. The sheet point is carried by the polynomial, then
//! landed by the geodesic sheet: smooth, and on the surface exactly.

use core::f64::consts::TAU;

use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::{Point, Point2, SparseMatrix, least_squares_cgnr};

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
    /// The flattened sheet's reach: the least distance from the origin to
    /// the flattened mesh's boundary.
    pub(crate) reach: f64,
    /// The stated departure from the smooth least-squares conformal map,
    /// on the sheet: the mesh's and the polynomial's together.
    pub(crate) error: f64,
}

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
            let (mesh, flat, departure) = resolved(land, radius, resolution)?;
            let rim = rim_distance(&mesh, &flat);
            if rim >= reach {
                let (degree, x, y, scale, miss) = fitted(&mesh, &flat, rim, resolution)?;
                return Ok(Self {
                    degree,
                    x,
                    y,
                    scale,
                    reach: rim,
                    error: departure + miss,
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
                f64::from(u32::try_from(j).unwrap_or(0)) / f64::from(u32::try_from(count).unwrap_or(1))
            }
        };
        let at_inner = |j: usize| if inner == 0 { 0 } else { vertex(k - 1, j % inner) };
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

/// The least-squares conformal flattening of `mesh`, the origin pinned at
/// the origin and its first neighbour on `+x`, scaled to keep the mesh's
/// area.
fn flattened(mesh: &Mesh) -> OgeomResult<Vec<Point2>> {
    let n = mesh.sheet.len();
    // The origin and the rim's vertex on the sheet's `+x` are pinned: two
    // points a whole radius apart fix the flattening's similarity with no
    // more than the solve's own error, where two a ring apart would turn it
    // by the error of one short edge.
    let rim = mesh.outer[0];
    let pinned = |v: usize| v == 0 || v == rim;
    let pin = |v: usize| {
        if v == 0 {
            Point2::ORIGIN
        } else {
            Point2::new(mesh.sheet[rim].x, 0.0)
        }
    };
    // The free vertices' columns, two each.
    let column = |v: usize| if v < rim { 2 * (v - 1) } else { 2 * (v - 2) };
    let mut triplets = Vec::with_capacity(mesh.triangles.len() * 12);
    let mut rhs = Vec::with_capacity(mesh.triangles.len() * 2);
    let mut space_area = 0.0;
    for (row, tri) in mesh.triangles.iter().enumerate() {
        let [a, b, c] = tri.map(|v| mesh.space[v]);
        let (ab, ac) = (b - a, c - a);
        let normal = ab.cross(ac);
        let double = normal.magnitude();
        if double.is_nan() || double <= 0.0 || ab.magnitude() <= 0.0 {
            ogeom_bail!(
                Construction,
                "the least-distortion sheet's mesh has a triangle of no area on the surface"
            );
        }
        space_area += double / 2.0;
        // The triangle in its own plane: `a` at the origin, `b` on the
        // first axis, `c` on the side the normal turns `b` towards.
        let e1 = ab * (1.0 / ab.magnitude());
        let e2 = normal.cross(ab) * (1.0 / (double * ab.magnitude()));
        let local = [
            (0.0, 0.0),
            (ab.magnitude(), 0.0),
            (ac.dot(e1), ac.dot(e2)),
        ];
        let weight = 1.0 / (double / 2.0).sqrt();
        let mut fixed = (0.0, 0.0);
        for j in 0..3 {
            // The edge opposite corner j, from j + 1 to j + 2.
            let (p1, p2) = (local[(j + 1) % 3], local[(j + 2) % 3]);
            let (da, db) = ((p2.0 - p1.0) * weight, (p2.1 - p1.1) * weight);
            let v = tri[j];
            // f d with f = U + iV, d = da + i db: real aU - bV, imaginary
            // bU + aV.
            if pinned(v) {
                let f = pin(v);
                fixed.0 += da * f.x - db * f.y;
                fixed.1 += db * f.x + da * f.y;
            } else {
                triplets.push((2 * row, column(v), da));
                triplets.push((2 * row, column(v) + 1, -db));
                triplets.push((2 * row + 1, column(v), db));
                triplets.push((2 * row + 1, column(v) + 1, da));
            }
        }
        rhs.push(-fixed.0);
        rhs.push(-fixed.1);
    }
    let unknowns = 2 * (n - 2);
    let system = SparseMatrix::from_triplets(rhs.len(), unknowns, &triplets);
    // The relative gradient the solve is held to, a hundred times the
    // arithmetic's resolution over the system's own entries; the
    // iterations bounded by twice the unknowns, past which conjugate
    // gradients has had every direction twice.
    let Some(solution) = least_squares_cgnr(&system, &rhs, 1e2 * f64::EPSILON.sqrt(), 2 * unknowns)
    else {
        ogeom_bail!(
            Construction,
            "the least-distortion sheet's flattening did not converge"
        );
    };
    let mut flat: Vec<Point2> = (0..n)
        .map(|v| {
            if pinned(v) {
                pin(v)
            } else {
                Point2::new(solution[column(v)], solution[column(v) + 1])
            }
        })
        .collect();
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
    // Scaled to keep the mesh's area, and turned so the first ring stands
    // at its geodesic directions on the mean: the sheet's `+x` at the
    // origin, the turn of each of the ring's vertices from its own
    // direction being the map's turn there plus a term in the vertex's
    // position, which the ring's six vertices sum to nothing.
    let scale = (space_area / flat_area).sqrt();
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
        *f = Point2::new(
            scale * cos.mul_add(f.x, -sin * f.y),
            scale * sin.mul_add(f.x, cos * f.y),
        );
    }
    Ok(flat)
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
/// miss.
type Fitted = (usize, Vec<f64>, Vec<f64>, f64, f64);

fn fitted(mesh: &Mesh, flat: &[Point2], rim: f64, resolution: f64) -> OgeomResult<Fitted> {
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
            reach: rim,
            error: 0.0,
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
            "the least-distortion sheet's map cannot be held to {resolution}: its best \
             polynomial misses the flattening by {:.3e}",
            b.4
        ),
        None => ogeom_bail!(
            Construction,
            "the least-distortion sheet's map folds over at every degree tried"
        ),
    }
}
