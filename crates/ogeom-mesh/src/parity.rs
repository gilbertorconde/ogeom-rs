//! Inside or outside a closed mesh, by the parity of its crossings along a
//! line parallel to `x`.
//!
//! The triangles the line crosses are found exactly. A line through an edge
//! or a vertex is moved by `(ε, ε²)` in `(y, z)`, `ε` infinitesimal, a
//! simulation of simplicity, so that two triangles sharing an edge never
//! both count it and never both miss it. A point on the line is inside
//! where an odd number of the crossings lie ahead of it, each decided by
//! the exact side of the point on the triangle's plane, and on the mesh
//! where it lies exactly in a crossed triangle.
//!
//! The answer is the mesh's own wherever every edge is used an even number
//! of times, as [`Triangulation::is_closed`](ogeom_topo::Triangulation::is_closed)
//! asks of a closed mesh: such a mesh bounds the points it winds about an
//! odd number of times, and only those.

use ogeom_core::{Exact, P3, Predicates, Sign};
use ogeom_math::Point;

use crate::Bvh;

/// Where a point lies against a closed mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Lies {
    /// Strictly inside.
    Inside,
    /// Strictly outside.
    Outside,
    /// Exactly in a triangle the line through it crosses.
    On,
}

impl Lies {
    /// Whether the point is inside or on the mesh.
    #[must_use]
    pub const fn is_inside_or_on(self) -> bool {
        matches!(self, Self::Inside | Self::On)
    }
}

/// A triangle the line crosses.
#[derive(Debug, Clone, Copy)]
struct Crossing {
    corners: [P3; 3],
    /// The triangle's extent along the line.
    low: f64,
    high: f64,
    /// The sign of its normal along the line.
    facing: Sign,
}

/// The line along `x` through `(y, z)`, with the triangles of a mesh it
/// crosses.
#[derive(Debug, Clone)]
pub struct Line {
    y: f64,
    z: f64,
    crossings: Vec<Crossing>,
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
/// by `(ε, ε²)` in `(y, z)`; if so, the sign of the triangle's normal along
/// `x`.
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

impl Line {
    /// The line along `x` through `(y, z)`, and the triangles of `bvh` it
    /// crosses.
    #[must_use]
    pub fn through(bvh: &Bvh, y: f64, z: f64) -> Self {
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
        Self { y, z, crossings }
    }

    /// Where the point of the line at `x` lies against the mesh.
    #[must_use]
    pub fn lies(&self, x: f64) -> Lies {
        let node = [x, self.y, self.z];
        let mut ahead = 0_u32;
        for c in &self.crossings {
            if x < c.low {
                ahead += 1;
            } else if x <= c.high {
                // The node's side of the plane is the crossing's side of
                // the node along the line, times its facing.
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
    }
}
