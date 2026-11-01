//! Iterative closest points onto a mesh, from one start or many.

use nalgebra::{Matrix6, SymmetricEigen, Vector6};
use ogeom_core::{OgeomResult, ogeom_bail, parallel};
use ogeom_math::{Matrix3, Point, Quaternion, Transform, Vector};
use ogeom_mesh::Bvh;
use ogeom_topo::Triangulation;

use super::{FreeMotion, Icp, Registration, Starts, is_finite, orientation, rigid};
use crate::distance::{DistanceOptions, Pseudonormals, to_mesh_with};

/// Register `source` onto the mesh `target` was built over: the rigid
/// transform carrying the points onto it, from each of [`Icp::starts`],
/// the best scored. See the [module](super) for the rounds, the score and
/// the free motions.
///
/// `source_normals`, one per point, are needed with a facing limit and
/// used for nothing else.
///
/// # Errors
///
/// As [`Icp::check`];
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) for no points, a
/// point not finite, an empty mesh, or a facing limit without normals;
/// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) for normals
/// not one per point; and as
/// [`principal_frame`](crate::orient::principal_frame) for the principal
/// axes of a mesh with no area.
pub fn icp(
    source: &[Point],
    source_normals: Option<&[Vector]>,
    target: &Bvh,
    opts: &Icp,
) -> OgeomResult<Registration> {
    let mesh = Mesh::of(target)?;
    let samples = Samples::of(source, source_normals, opts)?;
    let starts = starts(&samples, target, opts)?;
    // The starts run in parallel, each on one worker, and are ranked in
    // their order.
    let runs = parallel::map_ordered(&starts, |_, &start| {
        mesh.register(&samples, opts, start)
            .map(|run| run.registration(&samples))
    });
    let mut best: Option<Registration> = None;
    for found in runs {
        let found = found?;
        if best.as_ref().is_none_or(|b| found.score < b.score) {
            best = Some(found);
        }
    }
    best.map_or_else(|| ogeom_bail!(Domain, "no start to register from"), Ok)
}

/// The mesh registered onto: its hierarchy, its normals and the cap on a
/// sample's distance in a score.
pub(super) struct Mesh<'a> {
    index: &'a Bvh,
    normals: Pseudonormals,
    cap: f64,
}

impl<'a> Mesh<'a> {
    pub(super) fn of(index: &'a Bvh) -> OgeomResult<Self> {
        if index.is_empty() {
            ogeom_bail!(Domain, "an empty mesh registers nothing");
        }
        Ok(Self {
            index,
            normals: Pseudonormals::of(index),
            cap: index.bounds().diagonal(),
        })
    }

    /// The score's cap: the farthest pairing distance, or the diagonal of
    /// the mesh's box, past which no point of the mesh lies from a sample
    /// inside it.
    fn cap(&self, opts: &Icp) -> f64 {
        opts.max_distance.unwrap_or(self.cap)
    }

    /// The samples moved by `pose`, paired with the mesh, and the pairs
    /// within the cutoff's threshold over them kept.
    pub(super) fn pair(&self, samples: &Samples, opts: &Icp, pose: &Pose) -> OgeomResult<Paired> {
        let m = pose.q.to_matrix();
        let moved: Vec<Point> = samples.points.iter().map(|&p| pose.apply(&m, p)).collect();
        let normals: Option<Vec<Vector>> = samples
            .normals
            .as_ref()
            .map(|n| n.iter().map(|&v| m * v).collect());
        let options = DistanceOptions {
            max: opts.max_distance.unwrap_or(f64::INFINITY),
            facing_limit: opts.facing_limit,
            signed: false,
        };
        let found = to_mesh_with(
            self.index,
            &self.normals,
            &moved,
            normals.as_deref(),
            &options,
        )?;
        let cap = self.cap(opts);
        let mut capped = 0.0;
        let mut pairs = Vec::with_capacity(found.len());
        for (i, closest) in found.into_iter().enumerate() {
            match closest {
                Some(c) => {
                    capped += c.distance.min(cap);
                    pairs.push(Pair {
                        sample: i,
                        at: moved[i],
                        foot: c.point,
                        normal: c.normal,
                        distance: c.distance,
                    });
                }
                None => capped += cap,
            }
        }
        let distances: Vec<f64> = pairs.iter().map(|p| p.distance).collect();
        let threshold = opts.cutoff.threshold(&distances);
        pairs.retain(|p| p.distance <= threshold);
        #[allow(
            clippy::cast_precision_loss,
            reason = "sample counts are far below 2^52"
        )]
        let score = capped / samples.points.len() as f64;
        Ok(Paired { kept: pairs, score })
    }

    /// The rounds from `start`, and the pairs at the pose they end at.
    pub(super) fn register(&self, samples: &Samples, opts: &Icp, start: Pose) -> OgeomResult<Run> {
        let reach = opts.tol.confusion();
        let mut pose = start;
        let mut converged = false;
        for _ in 0..opts.max_iter {
            let paired = self.pair(samples, opts, &pose)?;
            if paired.kept.is_empty() {
                break;
            }
            let step = if opts.point_to_plane {
                NormalMatrix::of(&paired.kept).step()
            } else {
                let pairs = paired.kept.iter().map(|p| (p.at, p.foot, 1.0));
                orientation(pairs).map_or(Pose::IDENTITY, |(q, t)| Pose { q, t })
            };
            let m = step.q.to_matrix();
            let moves = paired
                .kept
                .iter()
                .map(|p| (step.apply(&m, p.at) - p.at).magnitude())
                .fold(0.0, f64::max);
            pose = step.then(&pose);
            if moves <= reach {
                converged = true;
                break;
            }
        }
        let paired = self.pair(samples, opts, &pose)?;
        Ok(Run {
            pose,
            converged,
            paired,
        })
    }
}

/// The points registered by, and their normals where there are some.
pub(super) struct Samples {
    pub(super) points: Vec<Point>,
    pub(super) normals: Option<Vec<Vector>>,
}

impl Samples {
    /// The points chosen by `opts`' sample count and seed, in their order,
    /// with their normals.
    pub(super) fn of(
        points: &[Point],
        normals: Option<&[Vector]>,
        opts: &Icp,
    ) -> OgeomResult<Self> {
        opts.check()?;
        if points.is_empty() {
            ogeom_bail!(Domain, "no points to register");
        }
        if points.iter().any(|&p| !is_finite(p)) {
            ogeom_bail!(Domain, "a point to register is not finite");
        }
        if let Some(n) = normals
            && n.len() != points.len()
        {
            ogeom_bail!(Dimension, "{} normals for {} points", n.len(), points.len());
        }
        if opts.facing_limit.is_some() && normals.is_none() {
            ogeom_bail!(Domain, "a facing limit needs the points' normals");
        }
        let chosen = chosen(points.len(), opts.sample, opts.seed);
        Ok(Self {
            points: chosen.iter().map(|&i| points[i]).collect(),
            normals: normals.map(|n| chosen.iter().map(|&i| n[i]).collect()),
        })
    }

    /// The same samples moved by `f`, their normals by `g`.
    pub(super) fn mapped(&self, f: impl Fn(Point) -> Point, g: impl Fn(Vector) -> Vector) -> Self {
        Self {
            points: self.points.iter().map(|&p| f(p)).collect(),
            normals: self
                .normals
                .as_ref()
                .map(|n| n.iter().map(|&v| g(v)).collect()),
        }
    }
}

/// `want` of the indices below `count`, in increasing order, chosen
/// uniformly by Floyd's method from a generator seeded with `seed`; every
/// index where `want` is zero or not below `count`.
fn chosen(count: usize, want: usize, seed: u64) -> Vec<usize> {
    if want == 0 || want >= count {
        return (0..count).collect();
    }
    let mut state = seed;
    let mut next = move || {
        // SplitMix64.
        state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    };
    let mut taken = vec![false; count];
    for j in count - want..count {
        let bound = j as u64 + 1;
        let t = usize::try_from(next() % bound).unwrap_or(j);
        let pick = if taken[t] { j } else { t };
        taken[pick] = true;
    }
    (0..count).filter(|&i| taken[i]).collect()
}

/// A rigid motion, `p ↦ R p + t` with `R` the unit quaternion's rotation.
#[derive(Debug, Clone, Copy)]
pub(super) struct Pose {
    pub(super) q: Quaternion,
    pub(super) t: Vector,
}

impl Pose {
    pub(super) const IDENTITY: Self = Self {
        q: Quaternion::IDENTITY,
        t: Vector::ZERO,
    };

    /// The motion a rigid transform makes.
    fn of(t: &Transform, opts: &Icp) -> OgeomResult<Self> {
        Ok(Self {
            q: Quaternion::from_matrix(&t.linear(), opts.tol)?,
            t: t.translation_vector(),
        })
    }

    /// `p` moved, `m` being the rotation's matrix.
    pub(super) fn apply(&self, m: &Matrix3, p: Point) -> Point {
        Point::from_vector(*m * p.to_vector() + self.t)
    }

    /// This motion after `first`.
    fn then(&self, first: &Self) -> Self {
        let q = self.q * first.q;
        Self {
            q: q * (1.0 / q.norm()),
            t: self.q.to_matrix() * first.t + self.t,
        }
    }

    pub(super) fn transform(&self) -> Transform {
        rigid(self.q, self.t)
    }
}

/// A moved sample and its nearest point on the mesh.
#[derive(Debug, Clone, Copy)]
pub(super) struct Pair {
    pub(super) sample: usize,
    pub(super) at: Point,
    pub(super) foot: Point,
    /// The mesh's unit normal at the foot; zero where it has none.
    normal: Vector,
    pub(super) distance: f64,
}

/// One round's pairs.
pub(super) struct Paired {
    /// The pairs the cutoff kept, in the samples' order.
    pub(super) kept: Vec<Pair>,
    /// The capped mean distance of every sample.
    pub(super) score: f64,
}

/// The rounds from one start, ended.
pub(super) struct Run {
    pub(super) pose: Pose,
    pub(super) converged: bool,
    pub(super) paired: Paired,
}

impl Run {
    fn registration(self, samples: &Samples) -> Registration {
        let kept = &self.paired.kept;
        let free_directions = if kept.is_empty() {
            Vec::new()
        } else {
            NormalMatrix::of(kept).free()
        };
        Registration {
            transform: self.pose.transform(),
            rms: rms(kept),
            inliers: kept.len(),
            sampled: samples.points.len(),
            converged: self.converged,
            score: self.paired.score,
            free_directions,
        }
    }
}

/// The RMS distance of the pairs; zero for none.
pub(super) fn rms(pairs: &[Pair]) -> f64 {
    if pairs.is_empty() {
        return 0.0;
    }
    #[allow(clippy::cast_precision_loss, reason = "pair counts are far below 2^52")]
    let count = pairs.len() as f64;
    (pairs.iter().map(|p| p.distance * p.distance).sum::<f64>() / count).sqrt()
}

/// The point-to-plane least squares of a round's kept pairs, linearised
/// about their centroid: the unknowns are the rotation scaled by the RMS
/// radius and the translation, so each moves the samples by lengths.
struct NormalMatrix {
    a: Matrix6<f64>,
    g: Vector6<f64>,
    centre: Point,
    radius: f64,
    /// The RMS point-to-plane residual.
    noise: f64,
    count: usize,
}

impl NormalMatrix {
    fn of(pairs: &[Pair]) -> Self {
        let first = pairs[0].at;
        #[allow(clippy::cast_precision_loss, reason = "pair counts are far below 2^52")]
        let count = pairs.len() as f64;
        let shift = pairs.iter().fold(Vector::ZERO, |s, p| s + (p.at - first)) / count;
        let centre = first + shift;
        let spread = pairs
            .iter()
            .map(|p| (p.at - centre).square_magnitude())
            .sum::<f64>()
            / count;
        let radius = spread.sqrt();
        // A set all at one point turns no differently at any scale.
        let scale = if radius > 0.0 { radius } else { 1.0 };
        let mut a = Matrix6::zeros();
        let mut g = Vector6::zeros();
        let mut squares = 0.0;
        for p in pairs {
            let turn = (p.at - centre).cross(p.normal) / scale;
            let row = Vector6::new(turn.x, turn.y, turn.z, p.normal.x, p.normal.y, p.normal.z);
            let b = (p.foot - p.at).dot(p.normal);
            a += row * row.transpose();
            g += row * b;
            squares += b * b;
        }
        Self {
            a,
            g,
            centre,
            radius: scale,
            noise: (squares / count).sqrt(),
            count: pairs.len(),
        }
    }

    /// The eigenvalue below which a direction is zero to the rounding of
    /// the matrix's sums: each entry sums `count` terms, so its rounding
    /// is within `count · ε` of the largest eigenvalue.
    fn rounding(&self, largest: f64) -> f64 {
        #[allow(clippy::cast_precision_loss, reason = "pair counts are far below 2^52")]
        let count = self.count as f64;
        largest * count * f64::EPSILON
    }

    /// The increment minimising the linearised residuals, along the
    /// directions the pairs hold above rounding; the rotation applied
    /// exactly, about the centroid.
    fn step(&self) -> Pose {
        let eigen = SymmetricEigen::new(self.a);
        let largest = eigen.eigenvalues.max();
        let floor = self.rounding(largest);
        let mut x = Vector6::zeros();
        for (j, &lambda) in eigen.eigenvalues.iter().enumerate() {
            if lambda > floor && lambda > 0.0 {
                let v = eigen.eigenvectors.column(j);
                x += v * (v.dot(&self.g) / lambda);
            }
        }
        let omega = Vector::new(x[0], x[1], x[2]) / self.radius;
        let t = Vector::new(x[3], x[4], x[5]);
        let angle = omega.magnitude();
        let q = if angle > 0.0 {
            let axis = omega / angle;
            let (s, c) = (angle * 0.5).sin_cos();
            Quaternion::new(c, axis.x * s, axis.y * s, axis.z * s)
        } else {
            Quaternion::IDENTITY
        };
        // About the centroid: p ↦ c + R (p − c) + t.
        let c = self.centre.to_vector();
        Pose {
            q,
            t: c - q.to_matrix() * c + t,
        }
    }

    /// The motions the pairs do not hold: the eigenvectors whose
    /// eigenvalues are no larger than the largest times `(σ/r)²` or its
    /// rounding, in the basis that diagonalises their rotation parts, the
    /// smallest rotation first.
    fn free(&self) -> Vec<FreeMotion> {
        let eigen = SymmetricEigen::new(self.a);
        let largest = eigen.eigenvalues.max();
        let relative = self.noise / self.radius;
        let threshold = (largest * relative * relative).max(self.rounding(largest));
        let free: Vec<Vector6<f64>> = (0..6)
            .filter(|&j| eigen.eigenvalues[j] <= threshold)
            .map(|j| eigen.eigenvectors.column(j).into_owned())
            .collect();
        if free.is_empty() {
            return Vec::new();
        }
        let k = free.len();
        let turn = |v: &Vector6<f64>| nalgebra::Vector3::new(v[0], v[1], v[2]);
        let gram = nalgebra::DMatrix::from_fn(k, k, |i, j| turn(&free[i]).dot(&turn(&free[j])));
        let split = SymmetricEigen::new(gram);
        let mut order: Vec<usize> = (0..k).collect();
        order.sort_by(|&i, &j| {
            split.eigenvalues[i]
                .total_cmp(&split.eigenvalues[j])
                .then(i.cmp(&j))
        });
        order
            .into_iter()
            .map(|c| {
                let w = split.eigenvectors.column(c);
                let mut v = free
                    .iter()
                    .zip(w.iter())
                    .fold(Vector6::zeros(), |s, (f, &wi)| s + f * wi);
                // The sense whose largest component is positive.
                let big = (0..6)
                    .max_by(|&i, &j| v[i].abs().total_cmp(&v[j].abs()).then(j.cmp(&i)))
                    .unwrap_or(0);
                if v[big] < 0.0 {
                    v = -v;
                }
                FreeMotion {
                    about: self.centre,
                    radius: self.radius,
                    rotation: Vector::new(v[0], v[1], v[2]) / self.radius,
                    translation: Vector::new(v[3], v[4], v[5]),
                }
            })
            .collect()
    }
}

/// The poses `opts.starts` names.
fn starts(samples: &Samples, target: &Bvh, opts: &Icp) -> OgeomResult<Vec<Pose>> {
    match &opts.starts {
        Starts::Identity => Ok(vec![Pose::IDENTITY]),
        Starts::Given(given) => given.iter().map(|t| Pose::of(t, opts)).collect(),
        Starts::PrincipalAxes24 => principal_starts(samples, target, opts),
    }
}

/// The 24 rotations carrying the samples' principal axes onto the mesh's,
/// each with the samples' centroid moved onto the mesh's area centroid.
fn principal_starts(samples: &Samples, target: &Bvh, opts: &Icp) -> OgeomResult<Vec<Pose>> {
    let mut soup = Triangulation::new();
    for t in 0..u32::try_from(target.triangle_count()).unwrap_or(u32::MAX) {
        if let Some(corners) = target.triangle(t) {
            let first = u32::try_from(soup.positions.len()).unwrap_or(u32::MAX);
            soup.positions.extend(corners);
            soup.triangles.push([first, first + 1, first + 2]);
        }
    }
    let frame = crate::orient::principal_frame(&soup, opts.tol)?;
    let to = frame.to_matrix();
    let (centre, from) = point_axes(&samples.points);
    let mut poses = Vec::with_capacity(24);
    for perm in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        for signs in 0..8_u8 {
            let mut q = [[0.0; 3]; 3];
            for (row, &col) in perm.iter().enumerate() {
                q[row][col] = if signs >> row & 1 == 0 { 1.0 } else { -1.0 };
            }
            let q = Matrix3::new(q);
            if q.determinant() < 0.0 {
                continue;
            }
            let r = to * q * from.transposed();
            let rotation = Quaternion::from_matrix(&r, opts.tol)?;
            poses.push(Pose {
                q: rotation,
                t: frame.origin().to_vector() - rotation.to_matrix() * centre.to_vector(),
            });
        }
    }
    Ok(poses)
}

/// The points' centroid and their principal axes, the columns of a
/// rotation: the eigenvectors of their covariance, the third the cross of
/// the first two.
fn point_axes(points: &[Point]) -> (Point, Matrix3) {
    let first = points[0];
    #[allow(
        clippy::cast_precision_loss,
        reason = "point counts are far below 2^52"
    )]
    let count = points.len() as f64;
    let centre = first + points.iter().fold(Vector::ZERO, |s, &p| s + (p - first)) / count;
    let mut c = nalgebra::Matrix3::<f64>::zeros();
    for &p in points {
        let d = p - centre;
        let d = nalgebra::Vector3::new(d.x, d.y, d.z);
        c += d * d.transpose();
    }
    let eigen = SymmetricEigen::new(c);
    let column = |i: usize| {
        let v = eigen.eigenvectors.column(i);
        Vector::new(v[0], v[1], v[2])
    };
    let (x, y) = (column(0), column(1));
    (centre, Matrix3::from_columns(x, y, x.cross(y)))
}
