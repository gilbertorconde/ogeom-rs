//! Hollowing a closed mesh: a sealed cavity inside it, its wall a set
//! distance from the part's surface.
//!
//! # The distance field
//!
//! The part's distance field is sampled on a grid of cubic cells, a third
//! of the wall wide unless the caller says otherwise, over the part's box.
//! Each node is inside, outside or on the part by the crossings of the line
//! along `x` through it: the triangles that line crosses are found exactly,
//! a line through an edge or a vertex resolved by simulation of
//! simplicity, as [`boolean`](crate::boolean()) resolves its rays, and a
//! node is inside where an odd number of them lie ahead of it, each
//! decided by the exact side of the node on the triangle's plane. A node
//! inside takes its exact distance to the nearest triangle, through the
//! part's [`Bvh`], negated.
//!
//! Only distances near the wall decide anything. The field is 1-Lipschitz,
//! so a node deeper than the wall by more than the longest edge of the
//! tetrahedra below, a cell's diagonal `√3 h`, has every neighbour deeper
//! than the wall: a node with no triangle within `wall + 2h` takes that
//! distance and loses nothing. For the same reason a node outside needs no
//! distance unless `√3 h` reaches across the wall.
//!
//! # The wall
//!
//! Each cell is split into six tetrahedra about its diagonal along `x + y +
//! z`, the same in every cell, so neighbours share their faces' splits. The
//! cavity is where the field is below `-wall`; a node exactly at it counts
//! outside the cavity, so no node lies on the surface. In each tetrahedron
//! the surface crosses the edges joining a node in the cavity to one out of
//! it, at the point the field interpolates linearly to `-wall`, held at
//! least [`CLAMP`] of the edge from either end. Points anywhere inside the
//! edges give triangles that meet only along their shared edges, so the
//! cavity's surface is closed, 2-manifold and crosses itself nowhere; each
//! triangle is wound by the exact side of a node in the cavity, so that it
//! faces into the cavity, out of the part's material.
//!
//! The surface is then decimated with [`Target::Verified`] within a
//! sixteenth of a cell, every point of each side measured.
//!
//! # Checked, not assumed
//!
//! The cavity is added to the part as a second shell, with no boolean, so
//! the result is checked exactly as [`boolean`](crate::boolean()) checks
//! an operand: closed, 2-manifold, no degenerate triangle, no triangle of
//! the cavity meeting another triangle but at shared corners, and every
//! shell wound outward of the material, which puts each cavity strictly
//! inside the part. A decimated cavity that fails, or that misses its
//! tolerance, is replaced by the one marched; one that fails too is an
//! error.
//!
//! # Error
//!
//! A point of the cavity's wall, before decimation, is where the linear
//! interpolation of the field is `-wall`; the field differs from its
//! interpolation within a tetrahedron by at most the tetrahedron's
//! diameter, `√3 h`, and where the field is smooth with curvature `κ`, by
//! about `3 h² κ / 2`. Holding the point off the ends moves it by at most
//! `√3 h` [`CLAMP`], and the decimation by `h / 16`.

use ogeom_core::parallel;
use ogeom_core::{Exact, FastMap, OgeomResult, P3, Predicates, Sign, Tolerances, ogeom_bail};
use ogeom_math::Point;
use ogeom_topo::Triangulation;

use crate::Bvh;
use crate::boolean::check_solid;
use crate::simplify::{Deviation, Target, simplify};

/// The share of a grid edge's length a crossing is held from either end.
///
/// It moves a crossing by at most `√3 h / 1024`, under a five-hundredth of
/// the decimation's tolerance, and keeps the corners of every triangle at
/// least `h / 1024` apart along their edges, far above the rounding of any
/// coordinate a grid of under `2^40` cells along an axis holds, so no
/// triangle's corners round onto one line.
pub const CLAMP: f64 = 1.0 / 1024.0;

/// The grid a hollow's distance field is sampled on.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Grid {
    /// The spacing of the nodes; `None` for a third of the wall.
    pub cell: Option<f64>,
}

/// A hollowed mesh, and what made it.
#[derive(Debug, Clone, PartialEq)]
pub struct Hollowed {
    /// The part and its cavities: the input's vertices and triangles first,
    /// as they were, then each cavity's, wound into the cavity. It carries
    /// no normals and no parameters.
    pub mesh: Triangulation,
    /// The first triangle of the cavities in [`Hollowed::mesh`].
    pub first_cavity_triangle: usize,
    /// How many closed shells the cavities add: one for a cavity like a
    /// ball's, two for one about a hole of the part's own. Zero where no
    /// point of the part is deeper than the wall.
    pub shells: usize,
    /// The grid's spacing.
    pub cell: f64,
    /// The decimation's measurement against the marched cavity, or `None`
    /// where the marched cavity was kept as it was.
    pub deviation: Option<Deviation>,
}

/// Hollow a closed mesh, leaving a wall `wall` thick.
///
/// The part must bound a solid as [`boolean`](crate::boolean()) asks of an
/// operand: closed, 2-manifold, no degenerate triangle, no triangle meeting
/// another but at shared corners, wound outward; it may have cavities of
/// its own, and the wall keeps clear of them. The cavity is sampled on
/// `grid` and decimated within a sixteenth of a cell, with `tol`'s
/// confusion as the decimation's; see the [module documentation](self).
///
/// The field costs a nearest-triangle query for each node inside the part,
/// about `V / h³` of them, and the checks one exact pass over the part;
/// the decimation, which measures and repeats, costs most on a curved
/// cavity.
///
/// # Errors
///
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) if
/// `wall` or the cell is not a positive distance, and as
/// [`boolean`](crate::boolean()) for an operand if the part does not bound
/// a solid, named "part". [`OgeomError::Range`](ogeom_core::OgeomError::Range)
/// if the grid would have more than `u32::MAX` nodes along an axis or the
/// cavity more than `u32::MAX` vertices.
/// [`OgeomError::Invariant`](ogeom_core::OgeomError::Invariant) if the
/// marched cavity fails the check of the result, named "hollowed".
pub fn hollow(
    mesh: &Triangulation,
    wall: f64,
    grid: &Grid,
    tol: Tolerances,
) -> OgeomResult<Hollowed> {
    if !(wall.is_finite() && wall > 0.0) {
        ogeom_bail!(Construction, "a wall of {wall} is not a distance");
    }
    let cell = grid.cell.unwrap_or(wall / 3.0);
    if !(cell.is_finite() && cell > 0.0) {
        ogeom_bail!(Construction, "a cell of {cell} is not a distance");
    }
    let bvh = check_solid(mesh, "part", 0)?;
    let lattice = Lattice::over(&bvh, cell, wall)?;
    let marched = lattice.march(&bvh)?;
    let first = mesh.triangles.len();
    if marched.triangles.is_empty() {
        return Ok(Hollowed {
            mesh: join(mesh, &marched)?,
            first_cavity_triangle: first,
            shells: 0,
            cell,
            deviation: None,
        });
    }

    let decimated = simplify(
        &marched,
        Target::Verified {
            tolerance: cell / 16.0,
            share: 1.0,
            max_triangles: None,
        },
        tol,
    )?;
    if decimated.target_met {
        let joined = join(mesh, &decimated.mesh)?;
        if check_solid(&joined, "hollowed", first).is_ok() {
            return Ok(Hollowed {
                shells: shells(&decimated.mesh),
                mesh: joined,
                first_cavity_triangle: first,
                cell,
                deviation: decimated.deviation,
            });
        }
    }
    let joined = join(mesh, &marched)?;
    if let Err(e) = check_solid(&joined, "hollowed", first) {
        ogeom_bail!(Invariant, "the marched cavity does not seal the part: {e}");
    }
    Ok(Hollowed {
        shells: shells(&marched),
        mesh: joined,
        first_cavity_triangle: first,
        cell,
        deviation: None,
    })
}

/// The part's vertices and triangles, then the cavity's.
fn join(part: &Triangulation, cavity: &Triangulation) -> OgeomResult<Triangulation> {
    let Ok(base) = u32::try_from(part.positions.len()) else {
        ogeom_bail!(Range, "a part of more than four billion vertices");
    };
    if u32::try_from(part.positions.len() + cavity.positions.len()).is_err() {
        ogeom_bail!(Range, "a hollowed part of more than four billion vertices");
    }
    let mut out = Triangulation::new();
    out.positions = part.positions.clone();
    out.positions.extend_from_slice(&cavity.positions);
    out.triangles = part.triangles.clone();
    out.triangles
        .extend(cavity.triangles.iter().map(|t| t.map(|v| v + base)));
    out.deflection_met = part.deflection_met;
    Ok(out)
}

/// The root of `x`'s set, halving the path to it.
fn find(parent: &mut [u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        parent[x as usize] = parent[parent[x as usize] as usize];
        x = parent[x as usize];
    }
    x
}

/// How many connected pieces a mesh's triangles make, joined at vertices.
fn shells(mesh: &Triangulation) -> usize {
    let n = u32::try_from(mesh.positions.len()).unwrap_or(u32::MAX);
    let mut parent: Vec<u32> = (0..n).collect();
    let mut used = vec![false; parent.len()];
    for t in &mesh.triangles {
        let r = find(&mut parent, t[0]);
        for &v in t {
            used[v as usize] = true;
            let s = find(&mut parent, v);
            if s != r {
                parent[s as usize] = r;
            }
        }
    }
    (0..n)
        .filter(|&v| used[v as usize] && find(&mut parent, v) == v)
        .count()
}

/// Where a node lies against the part.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Lies {
    Inside,
    Outside,
    On,
}

/// A triangle the line through a row of nodes crosses.
struct Crossing {
    corners: [P3; 3],
    /// The triangle's extent along the line.
    low: f64,
    high: f64,
    /// The sign of its normal along the line.
    facing: Sign,
}

/// The grid: its first node, spacing and node counts.
struct Lattice {
    origin: P3,
    cell: f64,
    count: [u32; 3],
    wall: f64,
    /// The distance beyond which a node's own distance does not matter.
    reach: f64,
    /// Whether nodes outside the part need their distances.
    outside_counts: bool,
}

fn coords(p: Point) -> P3 {
    [p.x, p.y, p.z]
}

/// The sign of `a - b`, exactly.
fn order(a: f64, b: f64) -> Sign {
    if a > b {
        Sign::Positive
    } else if a < b {
        Sign::Negative
    } else {
        Sign::Zero
    }
}

/// Whether the line along `x` through `(y, z)` crosses the triangle, moved
/// by `(ε, ε²)` in `(y, z)`, `ε` infinitesimal; if so, the sign of the
/// triangle's normal along `x`.
fn crosses(p: &[P3; 3], y: f64, z: f64) -> Option<Sign> {
    let flat = p.map(|c| [c[1], c[2]]);
    let facing = Exact::orient2d(flat[0], flat[1], flat[2]);
    if facing.is_zero() {
        return None;
    }
    for k in 0..3 {
        let (a, b) = (p[k], p[(k + 1) % 3]);
        let mut o = Exact::orient2d(flat[k], flat[(k + 1) % 3], [y, z]);
        if o.is_zero() {
            o = order(a[2], b[2]);
        }
        if o.is_zero() {
            o = order(b[1], a[1]);
        }
        if o != facing {
            return None;
        }
    }
    Some(facing)
}

/// How many nodes of a grid fit `extent` at `cell`, the last at or past it.
fn nodes_along(extent: f64, cell: f64) -> OgeomResult<u32> {
    let count = (extent / cell).ceil() + 1.0;
    if !(count.is_finite() && count <= f64::from(u32::MAX)) {
        ogeom_bail!(
            Range,
            "a grid of {count} nodes along an axis, a cell of {cell} across {extent}"
        );
    }
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "a whole number in 1..=u32::MAX, checked above"
    )]
    Ok(count as u32)
}

/// The six tetrahedra of a cell, as corners numbered by their offsets'
/// bits, `x` the lowest: each runs from the low corner to the high one,
/// one axis at a time.
const TETRAHEDRA: [[u8; 4]; 6] = [
    [0, 1, 3, 7],
    [0, 1, 5, 7],
    [0, 2, 3, 7],
    [0, 2, 6, 7],
    [0, 4, 5, 7],
    [0, 4, 6, 7],
];

/// The cavity's surface while it is being marched.
struct Marching {
    out: Triangulation,
    /// The vertex on each grid edge crossed, keyed by the edge's low node
    /// times eight plus the bits of its offset.
    on_edge: FastMap<u64, u32>,
}

impl Lattice {
    fn over(bvh: &Bvh, cell: f64, wall: f64) -> OgeomResult<Self> {
        let bounds = bvh.bounds();
        let (Some(low), Some(high)) = (bounds.low(), bounds.high()) else {
            ogeom_bail!(NullObject, "the part has no triangles");
        };
        let size = high - low;
        Ok(Self {
            origin: coords(low),
            cell,
            count: [
                nodes_along(size.x, cell)?,
                nodes_along(size.y, cell)?,
                nodes_along(size.z, cell)?,
            ],
            wall,
            reach: 2.0f64.mul_add(cell, wall),
            outside_counts: 3.0f64.sqrt() * cell >= wall,
        })
    }

    fn at(&self, i: u32, j: u32, k: u32) -> P3 {
        [
            f64::from(i).mul_add(self.cell, self.origin[0]),
            f64::from(j).mul_add(self.cell, self.origin[1]),
            f64::from(k).mul_add(self.cell, self.origin[2]),
        ]
    }

    /// Where each node of row `j` of layer `k` lies against the part.
    fn row(&self, bvh: &Bvh, j: u32, k: u32) -> Vec<Lies> {
        let [_, y, z] = self.at(0, j, k);
        let mut crossings = Vec::new();
        for t in bvh.triangles_where(|b| match (b.low(), b.high()) {
            (Some(l), Some(h)) => l.y <= y && y <= h.y && l.z <= z && z <= h.z,
            _ => false,
        }) {
            let Some(corners) = bvh.triangle(t).map(|c| c.map(coords)) else {
                continue;
            };
            if let Some(facing) = crosses(&corners, y, z) {
                let xs = corners.map(|c| c[0]);
                crossings.push(Crossing {
                    corners,
                    low: xs[0].min(xs[1]).min(xs[2]),
                    high: xs[0].max(xs[1]).max(xs[2]),
                    facing,
                });
            }
        }
        (0..self.count[0])
            .map(|i| {
                let node = self.at(i, j, k);
                let x = node[0];
                let mut ahead = 0_u32;
                for c in &crossings {
                    if x < c.low {
                        ahead += 1;
                    } else if x <= c.high {
                        // The node's side of the plane is the crossing's
                        // side of the node along the line, times its facing.
                        let s = Exact::orient3d(c.corners[0], c.corners[1], c.corners[2], node);
                        if s.is_zero() {
                            return Lies::On;
                        }
                        if s == c.facing {
                            ahead += 1;
                        }
                    }
                }
                if ahead % 2 == 1 {
                    Lies::Inside
                } else {
                    Lies::Outside
                }
            })
            .collect()
    }

    /// The field at each node of layer `k`, row by row: the wall plus the
    /// node's signed distance, so the cavity is where it is negative.
    fn layer(&self, bvh: &Bvh, k: u32) -> Vec<f64> {
        let rows: Vec<u32> = (0..self.count[1]).collect();
        let lies: Vec<Lies> = parallel::map_ordered(&rows, |_, &j| self.row(bvh, j, k))
            .into_iter()
            .flatten()
            .collect();
        let nx = self.count[0];
        let mut asked = Vec::new();
        let mut points = Vec::new();
        let nx = nx as usize;
        for (at, &l) in lies.iter().enumerate() {
            if l == Lies::Inside || (l == Lies::Outside && self.outside_counts) {
                let (i, j) = (at % nx, at / nx);
                let [x, y, z] = self.at(
                    u32::try_from(i).unwrap_or(u32::MAX),
                    u32::try_from(j).unwrap_or(u32::MAX),
                    k,
                );
                asked.push(at);
                points.push(Point::new(x, y, z));
            }
        }
        let mut field = vec![self.wall; lies.len()];
        for (at, hit) in asked.into_iter().zip(bvh.nearest_many(&points, self.reach)) {
            let d = hit.map_or(self.reach, |h| h.distance);
            field[at] = if lies[at] == Lies::Inside {
                self.wall - d
            } else {
                self.wall + d
            };
        }
        field
    }

    /// The cavity's surface, wound into the cavity.
    fn march(&self, bvh: &Bvh) -> OgeomResult<Triangulation> {
        let mut m = Marching {
            out: Triangulation::new(),
            on_edge: FastMap::default(),
        };
        let [nx, ny, nz] = self.count;
        if nx < 2 || ny < 2 || nz < 2 {
            return Ok(m.out);
        }
        let mut below = self.layer(bvh, 0);
        for k in 1..nz {
            let above = self.layer(bvh, k);
            for j in 0..ny - 1 {
                for i in 0..nx - 1 {
                    let mut f = [0.0; 8];
                    for (bits, value) in f.iter_mut().enumerate() {
                        let layer = if bits & 4 == 0 { &below } else { &above };
                        let (di, dj) = (usize::from(bits & 1 == 1), usize::from(bits & 2 == 2));
                        *value = layer[(j as usize + dj) * nx as usize + i as usize + di];
                    }
                    let cavity = f.iter().filter(|&&v| v < 0.0).count();
                    if cavity == 0 || cavity == 8 {
                        continue;
                    }
                    for tet in TETRAHEDRA {
                        self.tetrahedron(&mut m, [i, j, k - 1], tet, &f)?;
                    }
                }
            }
            below = above;
        }
        m.out.deflection_met = false;
        Ok(m.out)
    }

    fn corner(&self, cube: [u32; 3], bits: u8) -> [u32; 3] {
        [
            cube[0] + u32::from(bits & 1),
            cube[1] + u32::from((bits >> 1) & 1),
            cube[2] + u32::from((bits >> 2) & 1),
        ]
    }

    fn point(&self, n: [u32; 3]) -> P3 {
        self.at(n[0], n[1], n[2])
    }

    /// The vertex where the surface crosses the edge from corner `inside`,
    /// in the cavity, to corner `outside`.
    fn vertex(
        &self,
        m: &mut Marching,
        cube: [u32; 3],
        (inside, outside): (u8, u8),
        f: &[f64; 8],
    ) -> OgeomResult<u32> {
        let low = self.corner(cube, inside & outside);
        let node = (u64::from(low[2]) * u64::from(self.count[1]) + u64::from(low[1]))
            * u64::from(self.count[0])
            + u64::from(low[0]);
        let key = node * 8 + u64::from(inside ^ outside);
        if let Some(&v) = m.on_edge.get(&key) {
            return Ok(v);
        }
        let (a, b) = (f[inside as usize], f[outside as usize]);
        let t = (a / (a - b)).clamp(CLAMP, 1.0 - CLAMP);
        let p = self.point(self.corner(cube, inside));
        let q = self.point(self.corner(cube, outside));
        let at = Point::new(
            t.mul_add(q[0] - p[0], p[0]),
            t.mul_add(q[1] - p[1], p[1]),
            t.mul_add(q[2] - p[2], p[2]),
        );
        let Ok(v) = u32::try_from(m.out.positions.len()) else {
            ogeom_bail!(Range, "a cavity of more than four billion vertices");
        };
        m.out.positions.push(at);
        m.on_edge.insert(key, v);
        Ok(v)
    }

    /// Add the triangle on three vertices, wound so that the corner
    /// `inside` of the cell, in the cavity, is in front of it.
    fn triangle(
        &self,
        m: &mut Marching,
        cube: [u32; 3],
        mut t: [u32; 3],
        inside: u8,
    ) -> OgeomResult<()> {
        let [a, b, c] = t.map(|v| coords(m.out.positions[v as usize]));
        match Exact::orient3d(a, b, c, self.point(self.corner(cube, inside))) {
            Sign::Negative => {}
            Sign::Positive => t.swap(1, 2),
            Sign::Zero => ogeom_bail!(
                Invariant,
                "a cavity triangle's plane passes through a node it separates"
            ),
        }
        m.out.triangles.push(t);
        Ok(())
    }

    fn tetrahedron(
        &self,
        m: &mut Marching,
        cube: [u32; 3],
        tet: [u8; 4],
        f: &[f64; 8],
    ) -> OgeomResult<()> {
        let (inside, outside): (Vec<u8>, Vec<u8>) = tet.iter().partition(|&&c| f[c as usize] < 0.0);
        match (inside.as_slice(), outside.as_slice()) {
            (&[a], &[b, c, d]) => {
                let t = [
                    self.vertex(m, cube, (a, b), f)?,
                    self.vertex(m, cube, (a, c), f)?,
                    self.vertex(m, cube, (a, d), f)?,
                ];
                self.triangle(m, cube, t, a)
            }
            (&[a, b, c], &[d]) => {
                let t = [
                    self.vertex(m, cube, (a, d), f)?,
                    self.vertex(m, cube, (b, d), f)?,
                    self.vertex(m, cube, (c, d), f)?,
                ];
                self.triangle(m, cube, t, a)
            }
            (&[a, b], &[c, d]) => {
                let ac = self.vertex(m, cube, (a, c), f)?;
                let ad = self.vertex(m, cube, (a, d), f)?;
                let bd = self.vertex(m, cube, (b, d), f)?;
                let bc = self.vertex(m, cube, (b, c), f)?;
                // Each triangle holds two edges from one corner in the
                // cavity, so its plane parts that corner from the corners
                // out of it.
                self.triangle(m, cube, [ac, ad, bd], a)?;
                self.triangle(m, cube, [ac, bd, bc], b)
            }
            _ => Ok(()),
        }
    }
}

#[cfg(test)]
#[path = "hollow/tests.rs"]
mod tests;
