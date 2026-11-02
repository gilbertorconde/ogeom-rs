//! Smoothing an editable mesh by projecting each vertex onto a local
//! quadric, keeping creases sharp where asked.
//!
//! # A pass
//!
//! Each vertex is moved along its neighbourhood's normal onto the local
//! quadric fitted about it, the one [`curvature`](crate::curvature()) reads:
//! a weighted least-squares height field over the vertex's tangent plane,
//! fitted to the [`NEIGHBOURS`](crate::curvature::NEIGHBOURS) vertices
//! nearest it within the smoothing's `size`. Every vertex of a pass is
//! moved from the positions the pass started from, so the result does not
//! depend on the order the vertices are taken in, and is the same at any
//! thread count. Smoothing makes [`PASSES`] passes.
//!
//! # Creases
//!
//! With [`Smooth::keep_edges`] an angle, a neighbour counts towards a
//! vertex's quadric only if its normal turns from the vertex's own by no
//! more than that angle. A vertex beside a crease then fits the faces on
//! its own side, and stays on them. A vertex on the crease, whose normal
//! lies between the faces', is fitted to the vertices along the crease, or
//! where they fix no quadric, as at a corner, stays where it is.
//!
//! # How far a vertex moves
//!
//! No vertex ends farther than [`Smooth::max_move`] from where it started:
//! a move past it is shortened, along its line, until the distance as
//! computed is no more than `max_move`.
//!
//! With [`Smooth::region`], only the vertices in it move, each by its
//! pass's move scaled by its distance to the nearest vertex outside the
//! region over `size`, up to the whole move from one `size` in. The region
//! thus blends into the rest of the mesh rather than stepping off it.
//!
//! # What stays
//!
//! Only positions change. A vertex with no face, and every vertex
//! [`Editable::non_manifold`] reports, the copies of a point several fans
//! meet and the ends of edges not cleanly joined, stays where it is, so a
//! point's copies stay together. How faces are joined across an edge of
//! three or more of them depends on where their far corners stand, so
//! those edges are joined again afterwards, by the rule in the [module
//! documentation](super).

use super::{Editable, NONE, Scope, index};
use crate::curvature::{Cloud, check_size, in_parallel, vertex_normals};
use crate::point_index::PointIndex;
use ogeom_core::{OgeomResult, ogeom_bail};
use ogeom_math::Point;

/// How many passes a smoothing makes.
pub const PASSES: usize = 2;

/// How a mesh is smoothed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Smooth<'a> {
    /// How far about each vertex its neighbours are taken from.
    pub size: f64,
    /// The farthest a vertex may end from where it started.
    pub max_move: f64,
    /// The angle, in radians, beyond which a neighbour's normal turns too
    /// far from a vertex's for it to count; `None` to count every
    /// neighbour.
    pub keep_edges: Option<f64>,
    /// Which vertices may move, one entry per vertex number; `None` for
    /// all.
    pub region: Option<&'a [bool]>,
}

/// What a smoothing moved.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SmoothReport {
    /// The farthest a vertex ended from where it started.
    pub max: f64,
    /// The root mean square of how far each vertex with a face ended from
    /// where it started.
    pub rms: f64,
    /// How many vertices `max_move` held back.
    pub limited: usize,
}

impl Editable {
    /// Smooth the mesh, moving its vertices onto local quadrics fitted to
    /// their neighbours, as the [module documentation](crate::editable::smooth)
    /// describes. Vertex and face numbers stay; only an edge of three or
    /// more faces may be joined anew, as its faces now stand.
    ///
    /// # Errors
    ///
    /// [`OgeomError::Numeric`](ogeom_core::OgeomError::Numeric) if `size` is
    /// not positive and finite, `max_move` is negative or not a number, or
    /// `keep_edges` is not an angle from zero to `π`;
    /// [`OgeomError::Range`](ogeom_core::OgeomError::Range) if `region` does
    /// not have one entry per vertex number.
    pub fn smooth(&mut self, smooth: &Smooth<'_>) -> OgeomResult<SmoothReport> {
        let Smooth {
            size,
            max_move,
            keep_edges,
            region,
        } = *smooth;
        check_size(size)?;
        if max_move.is_nan() || max_move < 0.0 {
            ogeom_bail!(Numeric, "a largest move of {max_move} is not a distance");
        }
        if let Some(angle) = keep_edges
            && !(0.0..=core::f64::consts::PI).contains(&angle)
        {
            ogeom_bail!(Numeric, "a crease angle of {angle} is not from 0 to π");
        }
        let n = self.positions.len();
        if let Some(region) = region
            && region.len() != n
        {
            ogeom_bail!(
                Range,
                "a region of {} entries for {n} vertices",
                region.len()
            );
        }

        let with_face: Vec<u32> = (0..n).filter(|&v| self.out[v] != NONE).map(index).collect();
        let blend = self.blend(&with_face, size, region)?;
        let movers: Vec<u32> = with_face
            .iter()
            .copied()
            .filter(|&v| blend[v as usize] > 0.0)
            .collect();
        let faces: Vec<[u32; 3]> = self
            .corners
            .iter()
            .zip(&self.live)
            .filter(|&(_, &live)| live)
            .map(|(&t, _)| t)
            .collect();
        let cos = keep_edges.map(f64::cos);
        let start = self.positions.clone();
        let mut limited = vec![false; n];
        for _ in 0..PASSES {
            let normals = vertex_normals(&self.positions, faces.iter().copied());
            let members: Vec<u32> = with_face
                .iter()
                .copied()
                .filter(|&v| normals[v as usize].is_some())
                .collect();
            let cloud = Cloud::build(&self.positions, &normals, members)?;
            let targets = in_parallel(&movers, |v| {
                let p = self.positions[v as usize];
                let near = cloud.near(p, size);
                let local = match cos {
                    Some(cos) => {
                        let own = normals[v as usize]?;
                        cloud.fit(p, &near, |w| {
                            normals[w as usize].is_some_and(|m| m.dot(own) >= cos)
                        })
                    }
                    None => cloud.fit(p, &near, |_| true),
                }?;
                Some(local.projected())
            });
            drop(cloud);
            for (&v, target) in movers.iter().zip(targets) {
                let Some(target) = target else { continue };
                let k = v as usize;
                let p = self.positions[k];
                let moved = p + (target - p) * blend[k];
                let (q, held) = capped(start[k], moved, max_move);
                self.positions[k] = q;
                limited[k] = held;
            }
        }
        self.rejoin_crowded_edges();

        let (mut max, mut sum) = (0.0_f64, 0.0);
        for &v in &with_face {
            let d = (self.positions[v as usize] - start[v as usize]).magnitude();
            max = max.max(d);
            sum += d * d;
        }
        let rms = if with_face.is_empty() {
            0.0
        } else {
            (sum / f64::from(index(with_face.len()))).sqrt()
        };
        Ok(SmoothReport {
            max,
            rms,
            limited: limited.iter().filter(|&&l| l).count(),
        })
    }

    /// The share of its move each vertex makes: zero for a vertex with no
    /// face, one [`Editable::non_manifold`] reports, or one outside the
    /// region; otherwise its distance to the nearest vertex with a face
    /// outside the region over `size`, up to one.
    fn blend(
        &self,
        with_face: &[u32],
        size: f64,
        region: Option<&[bool]>,
    ) -> OgeomResult<Vec<f64>> {
        let mut blend = vec![0.0; self.positions.len()];
        let outside: Vec<Point> = with_face
            .iter()
            .filter(|&&v| region.is_some_and(|r| !r[v as usize]))
            .map(|&v| self.positions[v as usize])
            .collect();
        let outside = PointIndex::build(&outside)?;
        for &v in with_face {
            let k = v as usize;
            if region.is_none_or(|r| r[k]) {
                blend[k] = outside
                    .nearest(self.positions[k], size)
                    .map_or(1.0, |near| near.distance / size);
            }
        }
        for v in self.non_manifold() {
            blend[v as usize] = 0.0;
        }
        Ok(blend)
    }

    /// Join anew the faces on every edge that three or more of them use.
    pub(super) fn rejoin_crowded_edges(&mut self) {
        let mut edges: Vec<(u64, u32)> = Vec::new();
        for (f, &live) in self.live.iter().enumerate() {
            if live {
                for h in index(3 * f)..index(3 * f + 3) {
                    edges.push((self.edge_key(h), h));
                }
            }
        }
        edges.sort_unstable();
        let mut points: Vec<u32> = Vec::new();
        for run in edges.chunk_by(|a, b| a.0 == b.0) {
            if run.len() > 2 {
                let h = run[0].1;
                points.push(self.point[self.tail(h) as usize]);
                points.push(self.point[self.head(h) as usize]);
            }
        }
        if points.is_empty() {
            return;
        }
        points.sort_unstable();
        points.dedup();
        let around = self.about(&points);
        self.relink(&Scope::Points(points), &around);
    }
}

/// `to`, or where a move from `from` towards it is shortened to end no
/// farther than `most` from `from` as computed, and whether it was.
pub(super) fn capped(from: Point, to: Point, most: f64) -> (Point, bool) {
    let d = to - from;
    let length = d.magnitude();
    if length <= most {
        return (to, false);
    }
    let mut s = most / length;
    loop {
        let q = from + d * s;
        if (q - from).magnitude() <= most {
            return (q, true);
        }
        s = s.next_down();
    }
}
