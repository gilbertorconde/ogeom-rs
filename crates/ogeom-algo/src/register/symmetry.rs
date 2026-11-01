//! The plane a part is its own mirror image in.

use nalgebra::SymmetricEigen;
use ogeom_core::{OgeomResult, ogeom_bail, parallel};
use ogeom_math::{Direction, Plane, Point, Vector};
use ogeom_mesh::Bvh;

use super::icp::{Mesh, Pair, Pose, Samples, rms};
use super::{Icp, SymmetryFit};

/// The plane the points are the mirror image of the mesh `target` was
/// built over in, from each of `seeds`, the best scored.
///
/// From a seed, each round reflects the samples in the plane, registers
/// the reflection back onto the mesh from where it lies (as [`icp`] does
/// from the identity, whatever [`Icp::starts`] says), and bisects each kept
/// sample and its registered image: the new plane's normal is the
/// direction the samples' displacements to their images share (the
/// eigenvector of the largest eigenvalue of the sum of their outer
/// products, in the sense of the old normal), and it passes through the
/// mean of their midpoints. Where the reflection in a plane registers
/// with no motion the plane bisects every sample and its image already,
/// so the rounds stop once a round moves no sample's reflection further
/// than the confusion distance, or after [`Icp::max_iter`] rounds.
///
/// The seeds are ranked by the capped mean distance to the mesh of every
/// sample reflected in the plane each ends at, ties going to the earlier
/// seed; [`SymmetryFit::rms`] is the RMS of the kept ones'. `normals`, one
/// per point, are reflected with them and used only for a facing limit.
///
/// [`icp`]: fn@super::icp
///
/// # Errors
///
/// As [`Icp::check`];
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for no points, a
/// point not finite, an empty mesh or no seeds;
/// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) for normals
/// not one per point.
pub fn symmetry_plane(
    points: &[Point],
    normals: &[Vector],
    target: &Bvh,
    seeds: &[Plane],
    opts: &Icp,
) -> OgeomResult<SymmetryFit> {
    if seeds.is_empty() {
        ogeom_bail!(Domain, "no seed planes to start from");
    }
    let mesh = Mesh::of(target)?;
    let samples = Samples::of(points, Some(normals), opts)?;
    let fits = parallel::map_ordered(seeds, |_, seed| from_seed(&mesh, &samples, *seed, opts));
    let mut best: Option<SymmetryFit> = None;
    for found in fits {
        let found = found?;
        if best.as_ref().is_none_or(|b| found.score < b.score) {
            best = Some(found);
        }
    }
    best.map_or_else(|| ogeom_bail!(Domain, "no seed planes to start from"), Ok)
}

/// The rounds from one seed, and the fit of the plane they end at.
fn from_seed(mesh: &Mesh, samples: &Samples, seed: Plane, opts: &Icp) -> OgeomResult<SymmetryFit> {
    let reach = opts.tol.confusion();
    let mut plane = Mirror::of(&seed);
    let mut converged = false;
    for _ in 0..opts.max_iter {
        let mirrored = plane.reflected(samples);
        let run = mesh.register(&mirrored, opts, Pose::IDENTITY)?;
        let Some(next) = plane.bisecting(samples, &run.pose, &run.paired.kept, opts) else {
            break;
        };
        let moves = samples
            .points
            .iter()
            .map(|&p| (next.reflect(p) - plane.reflect(p)).magnitude())
            .fold(0.0, f64::max);
        plane = next;
        if moves <= reach {
            converged = true;
            break;
        }
    }
    let paired = mesh.pair(&plane.reflected(samples), opts, &Pose::IDENTITY)?;
    let normal = Direction::new(plane.normal, opts.tol)?;
    Ok(SymmetryFit {
        plane: Plane::through(plane.origin, normal),
        rms: rms(&paired.kept),
        matched: paired.kept.len(),
        sampled: samples.points.len(),
        score: paired.score,
        converged,
    })
}

/// A plane as a point on it and its unit normal.
#[derive(Debug, Clone, Copy)]
struct Mirror {
    origin: Point,
    normal: Vector,
}

impl Mirror {
    fn of(plane: &Plane) -> Self {
        let frame = plane.frame();
        Self {
            origin: frame.origin(),
            normal: frame.z().vector(),
        }
    }

    fn reflect(&self, p: Point) -> Point {
        p - self.normal * (2.0 * (p - self.origin).dot(self.normal))
    }

    fn reflect_vector(&self, v: Vector) -> Vector {
        v - self.normal * (2.0 * v.dot(self.normal))
    }

    /// The samples and their normals reflected in the plane.
    fn reflected(&self, samples: &Samples) -> Samples {
        samples.mapped(|p| self.reflect(p), |v| self.reflect_vector(v))
    }

    /// The plane bisecting each kept sample and its reflection registered
    /// by `pose`; `None` where the displacements share no direction.
    fn bisecting(&self, samples: &Samples, pose: &Pose, kept: &[Pair], opts: &Icp) -> Option<Self> {
        let first = samples.points[kept.first()?.sample];
        let m = pose.q.to_matrix();
        let mut spread = nalgebra::Matrix3::<f64>::zeros();
        let mut middle = Vector::ZERO;
        for pair in kept {
            let p = samples.points[pair.sample];
            let image = pose.apply(&m, self.reflect(p));
            let d = p - image;
            let d = nalgebra::Vector3::new(d.x, d.y, d.z);
            spread += d * d.transpose();
            middle += (p - first + (image - first)) * 0.5;
        }
        let eigen = SymmetricEigen::new(spread);
        let top = eigen.eigenvalues.imax();
        if eigen.eigenvalues[top].is_nan() || eigen.eigenvalues[top] <= 0.0 {
            return None;
        }
        let v = eigen.eigenvectors.column(top);
        let normal = Direction::new(Vector::new(v[0], v[1], v[2]), opts.tol)
            .ok()?
            .vector();
        let normal = if normal.dot(self.normal) < 0.0 {
            -normal
        } else {
            normal
        };
        #[allow(clippy::cast_precision_loss, reason = "pair counts are far below 2^52")]
        let count = kept.len() as f64;
        Some(Self {
            origin: first + middle / count,
            normal,
        })
    }
}
