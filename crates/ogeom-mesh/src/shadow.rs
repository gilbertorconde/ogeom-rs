//! The shadow a triangle mesh casts on a plane: the outline of everything
//! it covers seen along the plane's normal, with its holes.
//!
//! The shadow is the union of the triangles projected onto the plane. A
//! pocket or a dent in a face seen from above is covered by the face
//! around it and leaves no mark; an indentation in the outline itself
//! stays; a hole through the part, which no triangle covers, is a hole.
//! A triangle seen edge on covers no area and adds nothing, so a sheet
//! seen edge on casts no shadow.
//!
//! The vertices are projected and snapped to the grid of
//! [`planar`](crate::planar) once each, and every triangle's orientation is
//! read from those grid points exactly. Each triangle, turned to run
//! counter-clockwise, adds its three edges to a chain, and the edges two
//! triangles share in opposite directions cancel: those are the edges with
//! a triangle either side of them in the projection, which bound nothing.
//! What is left is the outline the mesh shows seen along the normal, its
//! folds and its open edges, and only that goes on to the region tracing:
//! the points covered by at least one triangle.
//!
//! The tracing runs in two rounds. The triangles are dealt to square tiles
//! of some thousands each by their first corner, and each tile's union is
//! traced on its own, in parallel: the folds a noisy scan has all over its
//! faces are covered within their tile and go no further. The tiles'
//! outlines are then traced together, the edges two tiles share
//! cancelling as before.

use ogeom_core::collections::FastMap;
use ogeom_core::parallel::map_ordered;
use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::{Plane, Point2};
use ogeom_topo::Triangulation;

use crate::planar::{Grid, Loop2, Q, Seg, estimate, orient, slot, trace, whole};

/// The shadow of a mesh on a plane, seen along its normal.
///
/// The loops are in the plane's own coordinates, along its frame's `x` and
/// `y`, with the shadow on their left: outer boundaries counter-clockwise,
/// holes clockwise. Each round of tracing moves the boundary by at most a
/// cell of the snapping grid, `2^-38` of the projected mesh's extent
/// rounded up to a power of two: within two cells of the exact union of the
/// projected triangles, as the projection rounds them.
///
/// # Errors
///
/// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if a triangle names
/// a vertex the mesh does not have, and
/// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if a position or
/// the plane's origin is not finite.
pub fn shadow(mesh: &Triangulation, plane: &Plane) -> OgeomResult<Vec<Loop2>> {
    if !plane.origin().is_finite() {
        ogeom_bail!(Numeric, "the shadow plane's origin is not finite");
    }
    if let Some(k) = mesh.positions.iter().position(|p| !p.is_finite()) {
        ogeom_bail!(Numeric, "mesh vertex {k} is not finite");
    }
    let count = mesh.positions.len();
    if let Some(&k) = mesh
        .triangles
        .iter()
        .flatten()
        .find(|&&k| k as usize >= count)
    {
        ogeom_bail!(Range, "a triangle names vertex {k} of {count}");
    }
    let frame = plane.frame();
    let flat: Vec<Point2> = mesh
        .positions
        .iter()
        .map(|&p| frame.to_local(p).xy())
        .collect();
    let Some(grid) = Grid::around(flat.iter().copied()) else {
        return Ok(Vec::new());
    };

    // One index per grid point, so vertices projecting together are one.
    let mut index: FastMap<Q, u32> = FastMap::default();
    let mut points: Vec<Q> = Vec::new();
    let id: Vec<u32> = flat
        .iter()
        .map(|&p| {
            let q = grid.snap(p);
            *index.entry(q).or_insert_with(|| {
                points.push(q);
                slot((points.len() - 1) as i64)
            })
        })
        .collect();

    // Every edge of every triangle seen counter-clockwise, keyed by the
    // tile the triangle falls in and its ends smaller first, with +1
    // running that way and -1 the other.
    let tiles = Tiles::over(&points, mesh.triangles.len());
    let mut edges: Vec<(u32, u64, i32)> = Vec::with_capacity(mesh.triangles.len() * 3);
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|k| id[k as usize]);
        let turn = orient(points[a as usize], points[b as usize], points[c as usize]);
        if turn == 0 {
            continue;
        }
        let ring = if turn > 0 { [a, b, c] } else { [a, c, b] };
        let tile = tiles.of(points[a as usize]);
        for k in 0..3 {
            let (u, v) = (ring[k], ring[(k + 1) % 3]);
            let (lo, hi, w) = if u < v { (u, v, 1) } else { (v, u, -1) };
            edges.push((tile, (u64::from(lo) << 32) | u64::from(hi), w));
        }
    }
    edges.sort_unstable_by_key(|e| (e.0, e.1));
    let mut ranges = Vec::new();
    let mut k = 0;
    while k < edges.len() {
        let from = k;
        while k < edges.len() && edges[k].0 == edges[from].0 {
            k += 1;
        }
        ranges.push(from..k);
    }
    let covered = |w: [i64; 2]| w[0] > 0;
    let per_tile = map_ordered(&ranges, |_, range| {
        let mut segments = Vec::new();
        let here = &edges[range.clone()];
        let mut k = 0;
        while k < here.len() {
            let key = here[k].1;
            let mut w = 0;
            while k < here.len() && here[k].1 == key {
                w += here[k].2;
                k += 1;
            }
            if w != 0 {
                segments.push(Seg {
                    a: points[(key >> 32) as usize],
                    b: points[(key & 0xffff_ffff) as usize],
                    w: [w, 0],
                });
            }
        }
        if ranges.len() == 1 {
            return Ok(segments);
        }
        Ok(rings(trace(segments, &covered)?))
    });
    let mut segments = Vec::new();
    for tile in per_tile {
        segments.extend(tile?);
    }
    Ok(trace(segments, &covered)?
        .into_iter()
        .map(|ring| Loop2 {
            points: ring.into_iter().map(|q| grid.point(q)).collect(),
        })
        .collect())
}

/// Loops as a chain, each edge counted once.
fn rings(loops: Vec<Vec<Q>>) -> Vec<Seg> {
    let mut out = Vec::new();
    for ring in loops {
        let n = ring.len();
        for k in 0..n {
            out.push(Seg {
                a: ring[k],
                b: ring[(k + 1) % n],
                w: [1, 0],
            });
        }
    }
    out
}

/// Square tiles over the grid, holding some thousands of triangles each.
struct Tiles {
    lo: Q,
    size: i64,
    columns: i64,
}

/// About how many triangles a tile holds.
const PER_TILE: usize = 4096;

impl Tiles {
    fn over(points: &[Q], triangles: usize) -> Self {
        let mut lo = [i64::MAX; 2];
        let mut hi = [i64::MIN; 2];
        for p in points {
            lo = [lo[0].min(p[0]), lo[1].min(p[1])];
            hi = [hi[0].max(p[0]), hi[1].max(p[1])];
        }
        let per_side = whole(estimate(triangles / PER_TILE).sqrt().ceil()).max(1);
        let span = (hi[0] - lo[0]).max(hi[1] - lo[1]) + 1;
        let size = (span + per_side - 1) / per_side;
        Self {
            lo,
            size,
            columns: (hi[0] - lo[0]) / size + 1,
        }
    }

    fn of(&self, p: Q) -> u32 {
        let (c, r) = (
            (p[0] - self.lo[0]) / self.size,
            (p[1] - self.lo[1]) / self.size,
        );
        slot(r * self.columns + c)
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, reason = "test code")]
    use super::*;
    use ogeom_math::Point;

    /// A unit cube seen from above casts a unit square, whatever the
    /// winding of its triangles.
    #[test]
    fn a_cube_casts_its_square() {
        let mut mesh = Triangulation::new();
        for z in [0.0, 1.0] {
            for (x, y) in [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)] {
                mesh.positions.push(Point::new(x, y, z));
            }
        }
        mesh.triangles = vec![
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 0, 4],
            [3, 4, 7],
        ];
        let loops = shadow(&mesh, &Plane::XY).unwrap();
        assert_eq!(loops.len(), 1);
        assert_eq!(loops[0].points.len(), 4);
        assert!((loops[0].area() - 1.0).abs() < 1e-12);
    }
}
