//! Rigid registration: the rotation and translation that best carry one set
//! of points onto another, onto a mesh, or a part onto its own mirror image.
//!
//! - [`absolute_orientation`] is the closed form for paired points.
//! - [`icp`] registers points onto a mesh through its [`Bvh`], pairing each
//!   point with its nearest point on the mesh and solving again until the
//!   pairs stop moving.
//! - [`symmetry_plane`] finds the plane a part is the mirror image of
//!   itself in.
//!
//! The linear algebra is `nalgebra`'s: a symmetric eigen-decomposition for
//! the quaternion of the closed form and for the normal matrix of each
//! point-to-plane step. Every sum runs over the pairs in their order and the
//! nearest-point queries are [`to_mesh`](crate::distance::to_mesh)'s, so
//! every answer is the same bit for bit on any thread count.
//!
//! # Iterative closest points
//!
//! Each round moves the samples by the transform so far, finds each one's
//! nearest point on the mesh (within [`Icp::max_distance`], and with
//! [`Icp::facing_limit`] among the triangles facing within the limit of the
//! sample's own normal), keeps the pairs within [`Icp::cutoff`]'s threshold
//! over the round's distances, and solves for an increment from the kept
//! pairs:
//!
//! - **point to plane** (the default): the small rotation and translation
//!   that minimise the squared distances of the samples to the planes
//!   through their nearest points, square to the mesh's normal there (the
//!   face's, or the pseudonormal at an edge or a corner). The increment is
//!   solved about the kept samples' centroid with the rotation measured in
//!   lengths at their RMS radius, so its six unknowns are commensurate, and
//!   its rotation is applied exactly;
//! - **point to point**: the [`absolute_orientation`] of the kept pairs.
//!
//! The rounds stop once an increment moves no kept sample further than the
//! confusion distance, which is convergence, or after [`Icp::max_iter`]
//! rounds. Each start of [`Icp::starts`] is run so, and the one with the
//! least [`Registration::score`] wins, ties going to the earlier start. The
//! score is the mean distance of *every* sample, each capped (at
//! [`Icp::max_distance`], or the mesh's diagonal where there is none) and a
//! sample with no pair counting as the cap, so a start that pairs a few
//! samples closely and leaves the rest does not beat one that pairs them
//! all.
//!
//! # Free motions
//!
//! Where the pairs cannot hold a motion (the samples of one flat face can
//! slide in its plane and spin about its normal), the point-to-plane normal
//! matrix is singular along it. A step moves only along the directions whose
//! eigenvalues stand above what rounding of the matrix's sums could make of
//! zero, and the registration reports the rest as
//! [`Registration::free_directions`]: the eigenvectors whose eigenvalues are
//! no larger than the largest times `(σ/r)²`, `σ` the RMS point-to-plane
//! residual and `r` the samples' RMS radius. Moving along such a direction
//! by the whole radius changes the residuals no more than moving along the
//! best-held direction by the noise does.

mod icp;
mod symmetry;

use ogeom_core::{OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Plane, Point, Quaternion, Transform, Vector};
#[cfg(doc)]
use ogeom_mesh::Bvh;

use crate::Cutoff;

pub use icp::icp;
pub use symmetry::symmetry_plane;

/// Where [`icp`] starts from.
#[derive(Debug, Clone, PartialEq, Default)]
pub enum Starts {
    /// The points as they are.
    #[default]
    Identity,
    /// The 24 rotations that carry the samples' principal axes onto the
    /// mesh's (every way of matching the three axes with their senses that
    /// keeps a right hand), each moving the samples' centroid onto the
    /// mesh's area centroid.
    ///
    /// Each start runs its own rounds over every sample, the starts in
    /// parallel. Starts far from the answer run all [`Icp::max_iter`]
    /// rounds, so a coarse search over a small [`Icp::sample`] refined
    /// from its answer by [`Starts::Given`] costs a fraction of one over
    /// many samples.
    PrincipalAxes24,
    /// The rigid transforms given, in their order.
    Given(Vec<Transform>),
}

/// How [`icp`] and [`symmetry_plane`] register.
#[derive(Debug, Clone, PartialEq)]
pub struct Icp {
    /// Whether a step minimises the distances to the planes through the
    /// nearest points rather than to the points themselves.
    pub point_to_plane: bool,
    /// Which pairs a round keeps, by their distances.
    pub cutoff: Cutoff,
    /// The farthest a sample's nearest point may be for the two to pair;
    /// `None` pairs at any distance.
    pub max_distance: Option<f64>,
    /// The greatest angle, in radians, between a sample's own normal and
    /// the triangles it may pair on, in `[0, π]`; `None` pairs on all of
    /// them.
    pub facing_limit: Option<f64>,
    /// Where to start from.
    pub starts: Starts,
    /// The most rounds per start.
    pub max_iter: u32,
    /// How many of the points to register by, chosen by a generator seeded
    /// with [`Icp::seed`]; zero, or at least the points' count, takes every
    /// point.
    pub sample: usize,
    /// The seed the sample is chosen with.
    pub seed: u64,
    /// The tolerances: an increment moving no sample further than their
    /// confusion distance has converged.
    pub tol: Tolerances,
}

impl Default for Icp {
    /// Point to plane, keeping the pairs within three times the median
    /// distance, at any distance and facing, from the identity, for up to
    /// a hundred rounds, every point, in millimetres.
    fn default() -> Self {
        Self {
            point_to_plane: true,
            cutoff: Cutoff::Median(3.0),
            max_distance: None,
            facing_limit: None,
            starts: Starts::Identity,
            max_iter: 100,
            sample: 0,
            seed: 0x5eed_1c9e,
            tol: Tolerances::millimetres(),
        }
    }
}

impl Icp {
    /// The options checked: the cutoff's own number, a positive farthest
    /// distance, a facing limit in `[0, π]`, and given starts that are
    /// there and rigid.
    ///
    /// # Errors
    ///
    /// As [`Cutoff::check`];
    /// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for any of the
    /// others.
    pub fn check(&self) -> OgeomResult<()> {
        self.cutoff.check()?;
        if let Some(d) = self.max_distance
            && (d.is_nan() || d <= 0.0)
        {
            ogeom_bail!(Domain, "the farthest pairing distance {d} must be positive");
        }
        if let Some(limit) = self.facing_limit
            && !(0.0..=core::f64::consts::PI).contains(&limit)
        {
            ogeom_bail!(Domain, "the facing limit {limit} must be in [0, π]");
        }
        if let Starts::Given(starts) = &self.starts {
            if starts.is_empty() {
                ogeom_bail!(Domain, "no starts are given");
            }
            for start in starts {
                if !is_rigid(start) {
                    ogeom_bail!(Domain, "a start is not a rigid motion: {start:?}");
                }
            }
        }
        Ok(())
    }
}

/// Whether a transform is a rotation and a translation: unit scale, and
/// a linear part that keeps a right hand. A transform composed of rigid
/// ones has a scale of exactly one.
fn is_rigid(t: &Transform) -> bool {
    t.scale_factor() == 1.0 && t.preserves_handedness()
}

/// What [`icp`] found.
#[derive(Debug, Clone, PartialEq)]
pub struct Registration {
    /// The rigid transform carrying the points onto the mesh.
    pub transform: Transform,
    /// The RMS distance of the kept samples to their nearest points, once
    /// moved.
    pub rms: f64,
    /// How many samples the last round kept.
    pub inliers: usize,
    /// How many samples the registration was made by.
    pub sampled: usize,
    /// Whether the last increment moved no kept sample further than the
    /// confusion distance.
    pub converged: bool,
    /// The capped mean distance of every sample, once moved, the starts
    /// were ranked by.
    pub score: f64,
    /// The motions the kept pairs do not hold, the most translational
    /// first. Empty where they hold all six.
    pub free_directions: Vec<FreeMotion>,
}

/// A rigid motion pairs do not hold: a turn at [`FreeMotion::rotation`]
/// about an axis through [`FreeMotion::about`], with a slide of `about` at
/// [`FreeMotion::translation`].
///
/// Scaled so that `|rotation| · radius` and `|translation|` have a unit
/// root sum of squares. Where the motions free are a subspace, they are
/// the basis of it that splits rotation from translation as far as it
/// allows: a selection free to slide and to spin reports pure slides and
/// a pure spin.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FreeMotion {
    /// The kept samples' centroid, once moved.
    pub about: Point,
    /// The kept samples' RMS distance from [`FreeMotion::about`].
    pub radius: f64,
    /// The angular velocity: the axis's direction, times the rate.
    pub rotation: Vector,
    /// The velocity of [`FreeMotion::about`].
    pub translation: Vector,
}

/// What [`symmetry_plane`] found.
#[derive(Debug, Clone, PartialEq)]
pub struct SymmetryFit {
    /// The plane the part is its own mirror image in.
    pub plane: Plane,
    /// The mirror residual: the RMS distance to the mesh of the kept
    /// samples reflected in the plane.
    pub rms: f64,
    /// How many reflected samples were kept.
    pub matched: usize,
    /// How many samples were reflected.
    pub sampled: usize,
    /// The capped mean distance of every reflected sample to the mesh, the
    /// seeds were ranked by.
    pub score: f64,
    /// Whether the last bisection moved no sample's reflection further
    /// than the confusion distance.
    pub converged: bool,
}

/// The rigid transform carrying the first point of each pair onto the
/// second with the least weighted sum of squared distances, in closed form.
///
/// The rotation is the unit quaternion of the largest eigenvalue of the
/// symmetric 4×4 matrix of the pairs' weighted cross-covariance about their
/// weighted centroids, and the translation then carries the first centroid
/// onto the second. A quaternion is a rotation by construction, so a set a
/// mirror would fit better (coplanar points, or noisy ones) still gets the
/// best proper rotation, never a reflection. Pairs that hold no rotation
/// (all one point) get none; pairs on one line get one about it that
/// carries them exactly. Equal weights where `weights` is `None`.
///
/// The cross-covariance holds the turn about a set's long axis in its
/// spread off that axis squared, so pairs spread by `δ` off a line of
/// length `L` are carried to within about `ε · L² / δ`: within `1e-12 · L`
/// for `δ` above some `2e-4 · L`, and exactly on the line itself.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for no pairs, a
/// point or weight not finite, a negative weight, or weights summing to
/// zero; [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) for
/// weights not one per pair.
pub fn absolute_orientation(
    pairs: &[(Point, Point)],
    weights: Option<&[f64]>,
) -> OgeomResult<Transform> {
    if pairs.is_empty() {
        ogeom_bail!(Domain, "no pairs to orient");
    }
    if let Some(w) = weights
        && w.len() != pairs.len()
    {
        ogeom_bail!(Dimension, "{} weights for {} pairs", w.len(), pairs.len());
    }
    if let Some(w) = weights
        && w.iter().any(|w| !(w.is_finite() && *w >= 0.0))
    {
        ogeom_bail!(Domain, "a weight is negative or not finite");
    }
    if pairs.iter().any(|(a, b)| !(is_finite(*a) && is_finite(*b))) {
        ogeom_bail!(Domain, "a paired point is not finite");
    }
    let weight = |i: usize| weights.map_or(1.0, |w| w[i]);
    let Some((q, t)) = orientation(
        pairs
            .iter()
            .copied()
            .enumerate()
            .map(|(i, (a, b))| (a, b, weight(i))),
    ) else {
        ogeom_bail!(Domain, "the weights sum to zero");
    };
    Ok(rigid(q, t))
}

fn is_finite(p: Point) -> bool {
    p.x.is_finite() && p.y.is_finite() && p.z.is_finite()
}

/// The rigid transform rotating by the unit quaternion `q` about the
/// origin, then translating by `t`.
fn rigid(q: Quaternion, t: Vector) -> Transform {
    Transform::translation(t) * Transform::from_quaternion(q)
}

/// The closed form over weighted pairs: the unit quaternion and the
/// translation, `p ↦ R p + t`. `None` where the weights sum to zero.
fn orientation(
    pairs: impl Iterator<Item = (Point, Point, f64)> + Clone,
) -> Option<(Quaternion, Vector)> {
    // The centroids are measured from the first pair's points, so far from
    // the origin the sums hold the spread rather than the offset.
    let (a0, b0, _) = pairs.clone().next()?;
    let mut total = 0.0;
    let (mut sa, mut sb) = (Vector::ZERO, Vector::ZERO);
    for (a, b, w) in pairs.clone() {
        total += w;
        sa += (a - a0) * w;
        sb += (b - b0) * w;
    }
    if total.is_nan() || total <= 0.0 {
        return None;
    }
    let (ca, cb) = (a0 + sa / total, b0 + sb / total);
    let mut s = [[0.0; 3]; 3];
    for (a, b, w) in pairs {
        let (u, v) = (a - ca, b - cb);
        let (u, v) = ([u.x, u.y, u.z], [v.x, v.y, v.z]);
        for (row, ui) in s.iter_mut().zip(u) {
            for (cell, vj) in row.iter_mut().zip(v) {
                *cell += w * ui * vj;
            }
        }
    }
    let q = if s.iter().flatten().all(|&x| x == 0.0) {
        Quaternion::IDENTITY
    } else {
        let [[sxx, sxy, sxz], [syx, syy, syz], [szx, szy, szz]] = s;
        let n = nalgebra::Matrix4::new(
            sxx + syy + szz,
            syz - szy,
            szx - sxz,
            sxy - syx,
            syz - szy,
            sxx - syy - szz,
            sxy + syx,
            szx + sxz,
            szx - sxz,
            sxy + syx,
            -sxx + syy - szz,
            syz + szy,
            sxy - syx,
            szx + sxz,
            syz + szy,
            -sxx - syy + szz,
        );
        let eigen = nalgebra::SymmetricEigen::new(n);
        let top = (0..4)
            .max_by(|&i, &j| {
                eigen.eigenvalues[i]
                    .total_cmp(&eigen.eigenvalues[j])
                    .then(j.cmp(&i))
            })
            .unwrap_or(0);
        let v = eigen.eigenvectors.column(top);
        let q = Quaternion::new(v[0], v[1], v[2], v[3]);
        let q = q * (1.0 / q.norm());
        // One sign of the two that are the same rotation, so the answer
        // does not depend on which the solver returned.
        if q.w < 0.0 { q * -1.0 } else { q }
    };
    let rotated = q.to_matrix() * ca.to_vector();
    Some((q, cb.to_vector() - rotated))
}
