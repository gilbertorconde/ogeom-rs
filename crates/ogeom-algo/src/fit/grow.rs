//! The region of a mesh one surface covers, or the crease one curve runs
//! along, grown from points clicked on it: [`grow`].

use std::collections::VecDeque;

use ogeom_core::{OgeomError, OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Axis, Circle, Cone, Cylinder, Plane, Point, Sphere, Torus, Vector};
use ogeom_mesh::Bvh;
use ogeom_topo::Triangulation;

use super::primitive::{self, Fit, FitFixed, FitOptions, Kind, Model, Side};
use super::tails;
use crate::Cutoff;
use crate::orient::{check_indices, edge_uses, key};

mod rim;

/// The surface kinds a region is grown for, and the curve kinds a crease is
/// traced for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum FitKind {
    /// A plane, [`plane`](super::plane).
    Plane,
    /// A sphere, [`sphere`](super::sphere).
    Sphere,
    /// A cylinder, [`cylinder`](super::cylinder).
    Cylinder,
    /// A cone, [`cone`](super::cone).
    Cone,
    /// A torus, [`torus`](super::torus).
    Torus,
    /// A circle in space, [`circle3`](super::circle3), traced along a
    /// crease: a bore's or a boss's rim.
    Circle3,
    /// A line in space, [`line3`](super::line3), traced along a crease: a
    /// box's edge.
    Line3,
}

impl FitKind {
    const fn kind(self) -> Kind {
        match self {
            Self::Plane => Kind::Plane,
            Self::Sphere => Kind::Sphere,
            Self::Cylinder => Kind::Cylinder,
            Self::Cone => Kind::Cone,
            Self::Torus => Kind::Torus,
            Self::Circle3 => Kind::Circle3,
            Self::Line3 => Kind::Line3,
        }
    }

    /// A curve, traced along a crease rather than grown over triangles.
    const fn curve(self) -> bool {
        matches!(self, Self::Circle3 | Self::Line3)
    }
}

/// How to grow a region.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GrowOptions {
    /// Which points each fit keeps, and how far from the model a vertex may
    /// lie and its triangle still join: for [`Cutoff::Sigma`] and
    /// [`Cutoff::Median`] the gate [`grow`] describes, `k σ` over
    /// many points; for [`Cutoff::Distance`] the distance, widened only by
    /// the fit's own uncertainty. [`Cutoff::All`] puts no limit on the
    /// distance, and growth is held by the normals alone.
    pub cutoff: Cutoff,
    /// The largest angle, in radians, between a triangle's normal and the
    /// model's for it to join; above zero and at most π/2. A scan's noise
    /// turns a triangle's normal by about the noise over the triangle's
    /// size, which the limit must clear.
    pub normal_limit: f64,
    /// Take the region's rim off, ring by ring, while it turns away from
    /// the model by clearly more than the noise; needs a cutoff with a
    /// threshold.
    pub peel: bool,
    /// The most triangles growth reaches; a seed's first rings are kept
    /// whole even past it.
    pub max_triangles: Option<usize>,
    /// What every fit holds fixed.
    pub fixed: FitFixed,
    /// Which side the material is on; where not given, the mesh's normals
    /// decide.
    pub side: Option<Side>,
    /// The angle, in radians, the surface must turn through across an
    /// edge, between its two triangles' normals, for a curve to be traced
    /// along it; above zero and below π. As with the normal limit, the
    /// noise turns an edge by about the noise over the triangles' size,
    /// which the limit must clear.
    pub crease_limit: f64,
}

impl Default for GrowOptions {
    /// A 3σ cutoff, which keeps 99.7 % of Gaussian noise, and a normal
    /// limit and a crease limit of π/6, which the noise of a mesh whose
    /// triangles are ten times the noise across turns past rarely.
    fn default() -> Self {
        Self {
            cutoff: Cutoff::Sigma(3.0),
            normal_limit: core::f64::consts::FRAC_PI_6,
            peel: false,
            max_triangles: None,
            fixed: FitFixed::default(),
            side: None,
            crease_limit: core::f64::consts::FRAC_PI_6,
        }
    }
}

/// A grown region's or a traced crease's fit, of the kind asked for.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum AnyFit {
    /// A plane's fit.
    Plane(Fit<Plane>),
    /// A sphere's fit.
    Sphere(Fit<Sphere>),
    /// A cylinder's fit.
    Cylinder(Fit<Cylinder>),
    /// A cone's fit.
    Cone(Fit<Cone>),
    /// A torus's fit.
    Torus(Fit<Torus>),
    /// A circle's fit, [`FitKind::Circle3`].
    Circle(Fit<Circle>),
    /// A line's fit, [`FitKind::Line3`].
    Line(Fit<Axis>),
}

impl AnyFit {
    /// The kind fitted.
    #[must_use]
    pub const fn kind(&self) -> FitKind {
        match self {
            Self::Plane(_) => FitKind::Plane,
            Self::Sphere(_) => FitKind::Sphere,
            Self::Cylinder(_) => FitKind::Cylinder,
            Self::Cone(_) => FitKind::Cone,
            Self::Torus(_) => FitKind::Torus,
            Self::Circle(_) => FitKind::Circle3,
            Self::Line(_) => FitKind::Line3,
        }
    }

    /// The RMS of the residuals of the points kept.
    #[must_use]
    pub fn sigma(&self) -> f64 {
        self.measured().0
    }

    /// Each point's signed residual.
    #[must_use]
    pub fn residuals(&self) -> &[f64] {
        self.measured().1
    }

    /// Whether each point was kept by the last round.
    #[must_use]
    pub fn used(&self) -> &[bool] {
        self.measured().2
    }

    fn measured(&self) -> (f64, &[f64], &[bool]) {
        match self {
            Self::Plane(f) => (f.sigma, &f.residuals, &f.used),
            Self::Sphere(f) => (f.sigma, &f.residuals, &f.used),
            Self::Cylinder(f) => (f.sigma, &f.residuals, &f.used),
            Self::Cone(f) => (f.sigma, &f.residuals, &f.used),
            Self::Torus(f) => (f.sigma, &f.residuals, &f.used),
            Self::Circle(f) => (f.sigma, &f.residuals, &f.used),
            Self::Line(f) => (f.sigma, &f.residuals, &f.used),
        }
    }

    fn model(&self) -> Model {
        match self {
            Self::Plane(f) => Model::plane(f),
            Self::Sphere(f) => Model::sphere(f),
            Self::Cylinder(f) => Model::cylinder(f),
            Self::Cone(f) => Model::cone(f),
            Self::Torus(f) => Model::torus(f),
            Self::Circle(f) => Model::circle3(f),
            Self::Line(f) => Model::line3(f),
        }
    }

    /// The cutoff's threshold over the residuals of the points kept.
    fn threshold(&self, cutoff: Cutoff) -> f64 {
        let (_, residuals, used) = self.measured();
        let kept: Vec<f64> = residuals
            .iter()
            .zip(used)
            .filter_map(|(r, keep)| keep.then_some(*r))
            .collect();
        cutoff.threshold(&kept)
    }
}

/// How a seed stands against the region's fit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeedReport {
    /// The triangle nearest the seed.
    pub triangle: u32,
    /// Whether that triangle is in the region; for a curve, whether a
    /// crease edge about the seed has both ends on the chain.
    pub in_region: bool,
    /// The seed's residual from the fit.
    pub offset: f64,
    /// The angle, in radians, between the triangle's normal and the fit's
    /// normal at the point of the triangle nearest the seed; for a curve,
    /// the angle the surface turns through across the nearest crease edge
    /// about the seed on the chain, or the nearest about it where none is.
    pub angle: f64,
}

/// A region grown over a mesh, or a crease traced along it, and its fit.
#[derive(Debug, Clone, PartialEq)]
pub struct Grown {
    /// The region's triangles, in increasing order; none for a curve.
    pub triangles: Vec<u32>,
    /// The fit over the region's vertices.
    pub fit: AnyFit,
    /// The mesh vertices the fit measured, in the order of its residuals:
    /// those of the region's triangles, in increasing order; for a curve,
    /// the chain's, in their order along it.
    pub vertices: Vec<u32>,
    /// Each seed's report, in the seeds' order.
    pub seeds: Vec<SeedReport>,
    /// The rings [`GrowOptions::peel`] took off.
    pub peeled_rings: u32,
}

/// No neighbour across an edge.
const NONE: u32 = u32::MAX;

/// Grow the region of `mesh` a surface of `kind` covers around `seeds`, and
/// fit it.
///
/// `index` is the hierarchy built over `mesh`. The mesh's triangles are
/// read as wound counter-clockwise about the material's outward normal,
/// which tells the fit the material's side where
/// [`GrowOptions::side`] does not.
///
/// # Seeds
///
/// Each seed starts at the triangle nearest it. For a plane, its
/// neighbourhood is that triangle and the rings of triangles round it whose
/// normals lie within [`GrowOptions::normal_limit`] of its own, until they
/// hold twice the points the fit needs. For a curved kind it is the region
/// a plane grows from there: the part of the surface round the seed that
/// stays within the normal limit of one plane, over which the kind is
/// determined far better than over a few triangles, fitted from a consensus
/// start so a neighbouring surface's points that joined it are cut off.
///
/// The first seed's neighbourhood is fitted, keeps the triangles its fit
/// admits, and grows. Each later seed whose neighbourhood the region has not
/// reached joins where its neighbourhood agrees with the grown fit, and the
/// region grows again; one that does not agree is left out, and reported
/// with its offset and angle from the fit.
///
/// # Growth
///
/// A neighbour across an edge joins where its three vertices lie within the
/// gate of the model and its normal within the normal limit of the model's
/// at its centroid. The gate is a prediction interval. A cutoff of `k`
/// sigmas over a fit of `p` parameters to `n` points measures the noise with
/// `n - p` degrees of freedom, so the gate is that estimate of the noise
/// times the Student t as rare as `k` sigmas of a known noise: as wide as a
/// few points' noise could be, and `k σ` over many. It widens further at a
/// point by how far the fit's own uncertainty moves the residual there, so a
/// model a small region determines poorly turns no triangle away for its
/// own error.
///
/// Where every coordinate of the mesh is a single-precision number, no gate
/// is below what rounding to it could do: a vertex is near within its own
/// coordinates' rounding along the model's normal, plus the cutoff's sigmas
/// of the fit's error from the region's rounding, so a noise-free `f32`
/// mesh grows whole however its rounding fell.
///
/// The model is refitted from the last fit whenever the region has doubled
/// since, and once more whenever no neighbour joins, every neighbour refused
/// before asked again; growth ends when a whole pass against the latest fit
/// admits none. So a bore stops at its end faces and a plane at its edges,
/// however the mesh runs on across them.
///
/// # Peeling
///
/// A scanner rounds a sharp edge, and the round's first rings lie within
/// the noise of the face they leave: they join, and pull the fit towards
/// them. [`GrowOptions::peel`] takes the region's rim off, ring by ring,
/// while the rim's normals turn away from the model's by clearly more than
/// the noise turns them.
///
/// How far a rim triangle turns, across its rim edges, is to first order a
/// weighted sum of its vertices' residuals: the gradient of the residuals
/// interpolated over the triangle, along the in-plane direction out of the
/// region. Summed over the rim, it is one weighted sum `S = Σ cᵥ rᵥ` of the
/// rim's vertices' residuals, a vertex shared by several rim triangles
/// weighing their weights together. Where the residuals are independent
/// noise of sigma `σ`, `S` is normal with standard deviation `σ √Σ cᵥ²`. A
/// ring is peeled where `|S|` exceeds the gate times `√Σ cᵥ²`, the gate
/// being `k σ` for [`Cutoff::Sigma`]`(k)` over many points: a ring of plain
/// noise is peeled as rarely as a point of plain noise is cut off. The
/// model is refitted after each ring.
///
/// # Curves
///
/// A [`FitKind::Circle3`] or a [`FitKind::Line3`] is traced along the
/// mesh's creases rather than grown over its triangles: a crease edge is
/// one two triangles share across which the surface turns, between their
/// normals, by more than [`GrowOptions::crease_limit`]. A seed reaches the
/// crease edges about the smooth part of the surface it lies on, the
/// triangles reached from the one nearest it across edges that are not
/// creases. The chain starts at the nearest of them whose chain holds
/// twice the vertices the kind's fit needs, so an edge the noise alone
/// turned past the limit is passed over; a seed with no such chain about
/// it is refused, by the nearest chain's fit's own refusal where one
/// refused and as lying on a smooth part otherwise.
///
/// The chain runs from the first seed's edge both ways, a vertex at each
/// end in turn. Before a fit has a finite gate (with [`Cutoff::All`],
/// never), it runs on through a vertex only where one crease edge leaves
/// it. Once one has, it runs on
/// along the crease edge whose far vertex lies within the gate of the fit,
/// the gate growth uses, its rounding floor included; where none does, or
/// where several do and the fit cannot tell which the crease follows, that
/// end stops. So a rim stops where its vertices leave the circle, and a
/// box's edge at its corners, where every edge on turns past the line. The
/// fit is refitted whenever the chain has doubled since the last, and once
/// more whenever both ends stop, each asked again; tracing ends when a
/// whole pass against the latest fit adds no vertex, or the chain closes.
/// A crease that turns at a vertex no other crease meets, before a fit
/// over few vertices has a gate narrow enough to see it, is followed round
/// the turn.
///
/// The answer's triangles are empty, and its vertices are the chain's.
/// Later seeds are reported against the fit, and do not extend the chain.
/// [`GrowOptions::normal_limit`] and [`GrowOptions::max_triangles`] apply
/// to regions only.
///
/// # Determinism
///
/// Everything runs on the calling thread in an order fixed by the seeds and
/// the triangle indices, and the consensus start draws a fixed sequence:
/// the same mesh and seeds give the same region and fit bit for bit on any
/// thread count.
///
/// # Errors
///
/// [`OgeomError::Domain`](ogeom_core::OgeomError::Domain) with no seeds, a
/// seed that is not finite, a normal limit or a crease limit out of range,
/// a cutoff whose number is not positive, peeling with [`Cutoff::All`] or
/// for a curve, a triangle cap of zero, or a triangle naming a vertex the
/// mesh lacks;
/// [`OgeomError::Dimension`](ogeom_core::OgeomError::Dimension) where the
/// hierarchy holds a different number of triangles from the mesh;
/// [`OgeomError::Construction`](ogeom_core::OgeomError::Construction) where
/// the mesh has no triangle to be near a seed, or for a curve where no
/// crease about a seed holds twice the vertices the kind needs; and the
/// kind's own refusals, as its fit gives them.
pub fn grow(
    mesh: &Triangulation,
    index: &Bvh,
    seeds: &[Point],
    kind: FitKind,
    options: &GrowOptions,
    tol: Tolerances,
) -> OgeomResult<Grown> {
    checked(mesh, index, seeds, kind, options)?;
    if kind.curve() {
        return rim::traced(mesh, index, seeds, kind, options, tol);
    }
    let surface = Surface::of(mesh);
    let least = 2 * primitive::needed(kind.kind(), &options.fixed);

    let mut hits = Vec::with_capacity(seeds.len());
    for (i, seed) in seeds.iter().enumerate() {
        let Some(hit) = index.nearest(*seed, f64::INFINITY) else {
            ogeom_bail!(Construction, "no triangle of the mesh lies near seed {i}");
        };
        hits.push(hit);
    }
    let mut grower = Grower {
        surface: &surface,
        kind,
        options,
        tol,
        in_region: vec![false; mesh.triangles.len()],
        region: Vec::new(),
        fit: None,
        gate: None,
        vertices: Vec::new(),
    };
    grower.seeded(&hits)?;
    let peeled_rings = if options.peel { grower.peel(least)? } else { 0 };

    let Grower {
        in_region,
        mut region,
        fit,
        gate,
        vertices,
        ..
    } = grower;
    let (Some(fit), Some(Gate { model, .. })) = (fit, gate) else {
        ogeom_bail!(Construction, "the seeds' neighbourhoods gave no fit");
    };
    let reports = seeds
        .iter()
        .zip(&hits)
        .map(|(seed, hit)| {
            let t = hit.triangle as usize;
            SeedReport {
                triangle: hit.triangle,
                in_region: in_region[t],
                offset: model.residual(*seed),
                angle: angle(&model, surface.normal[t], hit.point),
            }
        })
        .collect();
    region.sort_unstable();
    Ok(Grown {
        triangles: region,
        fit,
        vertices,
        seeds: reports,
        peeled_rings,
    })
}

fn checked(
    mesh: &Triangulation,
    index: &Bvh,
    seeds: &[Point],
    kind: FitKind,
    options: &GrowOptions,
) -> OgeomResult<()> {
    if seeds.is_empty() {
        ogeom_bail!(Domain, "growing a region needs at least one seed");
    }
    if let Some(i) = seeds.iter().position(|p| !p.is_finite()) {
        ogeom_bail!(Domain, "seed {i} is not finite");
    }
    let limit = options.normal_limit;
    if !(limit > 0.0 && limit <= core::f64::consts::FRAC_PI_2) {
        ogeom_bail!(
            Domain,
            "the normal limit {limit} must lie above 0 and at most pi/2"
        );
    }
    let crease = options.crease_limit;
    if !(crease > 0.0 && crease < core::f64::consts::PI) {
        ogeom_bail!(
            Domain,
            "the crease limit {crease} must lie above 0 and below pi"
        );
    }
    options.cutoff.check()?;
    if options.peel && kind.curve() {
        ogeom_bail!(
            Domain,
            "peeling takes off a region's rim, and a curve traced along a crease has none"
        );
    }
    if options.peel && options.cutoff == Cutoff::All {
        ogeom_bail!(
            Domain,
            "peeling measures a rim against the cutoff's threshold, and Cutoff::All has none"
        );
    }
    if options.max_triangles == Some(0) {
        ogeom_bail!(Domain, "a region capped at no triangles");
    }
    check_indices(mesh)?;
    if index.triangle_count() != mesh.triangles.len() {
        ogeom_bail!(
            Dimension,
            "the hierarchy holds {} triangles and the mesh {}",
            index.triangle_count(),
            mesh.triangles.len()
        );
    }
    Ok(())
}

/// The angle between a triangle's normal and the model's at a point; for
/// a model without the material's sense, the lesser of the two senses'.
fn angle(model: &Model, normal: Vector, at: Point) -> f64 {
    let cos = normal.dot(model.normal(at));
    let cos = if model.oriented() { cos } else { cos.abs() };
    cos.clamp(-1.0, 1.0).acos()
}

/// A mesh read for growing: each triangle's normal, doubled vector area,
/// centroid and neighbours across its edges, and the triangles about each
/// vertex.
struct Surface<'a> {
    mesh: &'a Triangulation,
    /// Unit, or zero for a triangle without area.
    normal: Vec<Vector>,
    /// Twice the area, along the normal.
    cross: Vec<Vector>,
    centroid: Vec<Point>,
    /// The triangle across each edge `(k, k + 1)` where exactly two use it.
    across: Vec<[u32; 3]>,
    fan_start: Vec<usize>,
    fan: Vec<u32>,
    /// Whether every coordinate is a single-precision number, as a mesh
    /// read from an `f32` file's are.
    single: bool,
}

impl<'a> Surface<'a> {
    fn of(mesh: &'a Triangulation) -> Self {
        let count = mesh.triangles.len();
        let mut normal = Vec::with_capacity(count);
        let mut cross = Vec::with_capacity(count);
        let mut centroid = Vec::with_capacity(count);
        for t in &mesh.triangles {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            let doubled = (b - a).cross(c - a);
            let length = doubled.magnitude();
            if length.is_finite() && length > 0.0 {
                normal.push(doubled / length);
                cross.push(doubled);
            } else {
                normal.push(Vector::ZERO);
                cross.push(Vector::ZERO);
            }
            centroid.push(a + ((b - a) + (c - a)) / 3.0);
        }
        let uses = edge_uses(&mesh.triangles);
        let across = mesh
            .triangles
            .iter()
            .enumerate()
            .map(|(t, tri)| {
                let mut out = [NONE; 3];
                for (k, slot) in out.iter_mut().enumerate() {
                    if let Some(u) = uses.get(&key(tri[k], tri[(k + 1) % 3]))
                        && u.count == 2
                    {
                        *slot = if u.first as usize == t {
                            u.second
                        } else {
                            u.first
                        };
                    }
                }
                out
            })
            .collect();
        let mut fan_start = vec![0_usize; mesh.positions.len() + 1];
        for tri in &mesh.triangles {
            for &v in tri {
                fan_start[v as usize + 1] += 1;
            }
        }
        for v in 0..mesh.positions.len() {
            fan_start[v + 1] += fan_start[v];
        }
        let mut fill = fan_start.clone();
        let mut fan = vec![0_u32; fan_start[mesh.positions.len()]];
        for (t, tri) in mesh.triangles.iter().enumerate() {
            for &v in tri {
                let slot = &mut fill[v as usize];
                fan[*slot] = as_index(t);
                *slot += 1;
            }
        }
        let single = mesh
            .positions
            .iter()
            .all(|p| [p.x, p.y, p.z].iter().all(|&x| is_single(x)));
        Self {
            mesh,
            normal,
            cross,
            centroid,
            across,
            fan_start,
            fan,
            single,
        }
    }

    fn corners(&self, t: u32) -> [u32; 3] {
        self.mesh.triangles[t as usize]
    }

    fn position(&self, v: u32) -> Point {
        self.mesh.positions[v as usize]
    }

    fn fan(&self, v: u32) -> &[u32] {
        let v = v as usize;
        &self.fan[self.fan_start[v]..self.fan_start[v + 1]]
    }

    /// The triangle and the rings of triangles about it, each ring those
    /// sharing a vertex with the last, at least one ring and until they
    /// hold `least` vertices or there are no more. Only triangles whose
    /// normals lie within `cos_limit` of the first's are taken, so a
    /// triangle beside an edge takes its own side's.
    fn patch(&self, start: u32, least: usize, cos_limit: f64) -> Vec<u32> {
        let facing = self.normal[start as usize];
        let mut triangles = vec![start];
        let mut vertices: Vec<u32> = self.corners(start).to_vec();
        vertices.sort_unstable();
        vertices.dedup();
        loop {
            let mut ring: Vec<u32> = vertices
                .iter()
                .flat_map(|&v| self.fan(v).iter().copied())
                .filter(|&t| {
                    facing == Vector::ZERO || self.normal[t as usize].dot(facing) >= cos_limit
                })
                .filter(|t| triangles.binary_search(t).is_err())
                .collect();
            ring.sort_unstable();
            ring.dedup();
            if ring.is_empty() {
                break;
            }
            for &t in &ring {
                vertices.extend(self.corners(t));
            }
            vertices.sort_unstable();
            vertices.dedup();
            triangles.extend(ring);
            triangles.sort_unstable();
            if vertices.len() >= least {
                break;
            }
        }
        triangles
    }

    /// The vertices of the triangles, in increasing order, each with the
    /// sum of its triangles' vector areas; a vertex of triangles without
    /// area alone is left out.
    fn vertices(&self, triangles: &[u32]) -> (Vec<u32>, Vec<Point>, Vec<Vector>) {
        let mut sums: Vec<(u32, Vector)> = triangles
            .iter()
            .flat_map(|&t| self.corners(t).map(|v| (v, self.cross[t as usize])))
            .collect();
        sums.sort_by_key(|(v, _)| *v);
        let mut vertices = Vec::new();
        let mut points = Vec::new();
        let mut normals = Vec::new();
        let mut i = 0;
        while i < sums.len() {
            let v = sums[i].0;
            let mut sum = Vector::ZERO;
            while i < sums.len() && sums[i].0 == v {
                sum += sums[i].1;
                i += 1;
            }
            if sum.magnitude() > 0.0 {
                vertices.push(v);
                points.push(self.position(v));
                normals.push(sum);
            }
        }
        (vertices, points, normals)
    }
}

/// Whether `x` is a single-precision number.
fn is_single(x: f64) -> bool {
    #[allow(
        clippy::cast_possible_truncation,
        reason = "the round trip is the test"
    )]
    let narrow = x as f32;
    f64::from(narrow) == x
}

/// The most rounding to single precision moves a number of `x`'s size:
/// half the gap from there to the next single-precision number up, which
/// is the wider gap about it.
fn half_ulp(x: f64) -> f64 {
    #[allow(
        clippy::cast_possible_truncation,
        reason = "the size of a single-precision coordinate, or past the largest"
    )]
    let size = x.abs() as f32;
    let up = f64::from(size.next_up());
    let gap = if up.is_finite() {
        up - f64::from(size)
    } else {
        f64::from(size) - f64::from(size.next_down())
    };
    0.5 * gap
}

fn as_index(t: usize) -> u32 {
    #[allow(
        clippy::cast_possible_truncation,
        reason = "the hierarchy holds the mesh's triangles, which it refuses past u32"
    )]
    let t = t as u32;
    t
}

/// What decides whether a triangle joins: the model and the distance and
/// angle it allows.
#[derive(Debug, Clone, Copy)]
struct Gate {
    model: Model,
    threshold: f64,
    cos_limit: f64,
    /// For a single-precision mesh, the most rounding moved a coordinate
    /// of the region's size; zero otherwise.
    quantum: f64,
    /// How many sigmas of the noise the cutoff keeps.
    sigmas: f64,
}

impl Gate {
    /// The gate of a fit to `points`.
    ///
    /// A cutoff of `k` sigmas (or a median's worth of them) measured over
    /// a fit of `p` parameters to `n` points, each missing the model in `c`
    /// directions (two for a curve in space, one for a surface), is a noise
    /// estimated with `c n - p` degrees of freedom: the threshold is the
    /// sigma that estimate gives, `σ √(c n / (c n - p))`, times the Student
    /// t at which a residual is as rare as `k` sigmas of a known noise. A few points then admit
    /// as much as their noise could, and many admit `k σ`. No threshold is
    /// below the confusion distance, within which points are one.
    ///
    /// On a single-precision mesh, no threshold is below the rounding's
    /// either. Rounding moves each coordinate by at most `q`, half the gap
    /// between single-precision numbers at the region's size, so it moves a
    /// residual along a unit normal by a sigma of at most `q`: the
    /// threshold is at least as many of those sigmas as the cutoff keeps of
    /// the fit's own.
    fn of(
        fit: &AnyFit,
        points: &[Point],
        single: bool,
        options: &GrowOptions,
        tol: Tolerances,
    ) -> Self {
        let model = fit
            .model()
            .measured_over(points, fit.used(), &options.fixed);
        let (mut threshold, mut quantum, mut sigmas) = (f64::INFINITY, 0.0, 0.0);
        match options.cutoff {
            Cutoff::All => {}
            Cutoff::Distance(d) => threshold = d,
            cutoff @ (Cutoff::Sigma(k) | Cutoff::Median(k)) => {
                let kept = model.components() * fit.used().iter().filter(|u| **u).count();
                let dof = kept.saturating_sub(model.parameters());
                let base = fit.threshold(cutoff);
                let rms = fit.sigma();
                sigmas = if rms > 0.0 { base / rms } else { k };
                threshold = if dof == 0 {
                    f64::INFINITY
                } else if rms > 0.0 {
                    #[allow(
                        clippy::cast_precision_loss,
                        reason = "point counts are far below 2^52"
                    )]
                    let unbiased = rms * (kept as f64 / dof as f64).sqrt();
                    tails::t_at_normal(sigmas, dof) * unbiased
                } else {
                    0.0
                };
                if single {
                    let size = points
                        .iter()
                        .flat_map(|p| [p.x, p.y, p.z])
                        .fold(0.0, |m: f64, x| m.max(x.abs()));
                    quantum = half_ulp(size);
                    threshold = threshold.max(sigmas * quantum);
                }
            }
        }
        Self {
            model,
            threshold: threshold.max(tol.confusion()),
            cos_limit: options.normal_limit.cos(),
            quantum,
            sigmas,
        }
    }

    /// Whether a point lies within the threshold, widened by how far the
    /// fit's own uncertainty moves its residual there: the prediction
    /// interval of a point the fit did not see, so a model a small region
    /// determines poorly turns no triangle away for its own error.
    ///
    /// On a single-precision mesh, a point is near at least where its
    /// residual is within what rounding alone could have moved it: its own
    /// coordinates' rounding, along the normal, at most (for a curve, the
    /// rounding's half diagonal), and the fit's
    /// error from the region's rounding, as many sigmas of it as the
    /// cutoff keeps.
    fn near(&self, p: Point) -> bool {
        let r = self.model.residual(p);
        let leverage = self.model.leverage(p);
        let mut bound = self.threshold * (1.0 + leverage).sqrt();
        if self.quantum > 0.0 {
            let half = Vector::new(half_ulp(p.x), half_ulp(p.y), half_ulp(p.z));
            let own = self.model.moved_by(p, half);
            bound = bound.max(own + self.sigmas * self.quantum * leverage.sqrt());
        }
        r.abs() <= bound
    }

    fn admits(&self, surface: &Surface<'_>, t: u32) -> bool {
        let i = t as usize;
        if surface.cross[i] == Vector::ZERO {
            return false;
        }
        let near = surface
            .corners(t)
            .iter()
            .all(|&v| self.near(surface.position(v)));
        near && {
            let cos = surface.normal[i].dot(self.model.normal(surface.centroid[i]));
            let cos = if self.model.oriented() {
                cos
            } else {
                cos.abs()
            };
            cos >= self.cos_limit
        }
    }
}

/// A region as it grows, and its latest fit.
struct Grower<'s, 'a> {
    surface: &'s Surface<'a>,
    kind: FitKind,
    options: &'s GrowOptions,
    tol: Tolerances,
    in_region: Vec<bool>,
    /// The region's triangles, in the order they joined.
    region: Vec<u32>,
    fit: Option<AnyFit>,
    gate: Option<Gate>,
    vertices: Vec<u32>,
}

impl Grower<'_, '_> {
    /// Fit the region, from the last fit where there is one.
    fn refit(&mut self) -> OgeomResult<Gate> {
        let (vertices, points, normals) = self.surface.vertices(&self.region);
        let fit = fitted(
            self.kind,
            &points,
            Some(&normals),
            self.options,
            self.fit.as_ref(),
            self.tol,
        )?;
        self.vertices = vertices;
        let gate = Gate::of(&fit, &points, self.surface.single, self.options, self.tol);
        self.gate = Some(gate);
        self.fit = Some(fit);
        Ok(gate)
    }

    fn set_region(&mut self, triangles: Vec<u32>) {
        for &t in &self.region {
            self.in_region[t as usize] = false;
        }
        for &t in &triangles {
            self.in_region[t as usize] = true;
        }
        self.region = triangles;
    }

    /// A seed's neighbourhood: for a plane, the rings about the hit
    /// triangle; for a curved kind, the region a plane grows from them
    /// with these options.
    fn neighbourhood(&self, hit: &ogeom_mesh::Hit) -> OgeomResult<Vec<u32>> {
        let least = 2 * primitive::needed(Kind::Plane, &FitFixed::default());
        let start = self
            .surface
            .patch(hit.triangle, least, self.options.normal_limit.cos());
        if self.kind == FitKind::Plane {
            return Ok(start);
        }
        let options = GrowOptions {
            peel: false,
            fixed: FitFixed::default(),
            side: None,
            ..*self.options
        };
        let mut plane = Grower {
            surface: self.surface,
            kind: FitKind::Plane,
            options: &options,
            tol: self.tol,
            in_region: vec![false; self.in_region.len()],
            region: Vec::new(),
            fit: None,
            gate: None,
            vertices: Vec::new(),
        };
        plane.set_region(start);
        plane.refit()?;
        plane.spread()?;
        plane.region.sort_unstable();
        Ok(plane.region)
    }

    /// Grow from the seeds in their order. The first neighbourhoods are
    /// fitted together until they determine a fit, which keeps the
    /// triangles it admits and grows. Each later neighbourhood the region
    /// has not reached joins where it agrees with the grown fit, with the
    /// triangles that fit admits, and the region grows again; one that does
    /// not agree is left out.
    fn seeded(&mut self, hits: &[ogeom_mesh::Hit]) -> OgeomResult<()> {
        let patches = hits
            .iter()
            .map(|h| self.neighbourhood(h))
            .collect::<OgeomResult<Vec<_>>>()?;
        let mut first: Vec<u32> = Vec::new();
        for (i, patch) in patches.iter().enumerate() {
            let Some(gate) = self.gate else {
                first.extend(patch);
                first.sort_unstable();
                first.dedup();
                self.set_region(first.clone());
                match self.refit() {
                    Ok(gate) => {
                        let kept: Vec<u32> = first
                            .iter()
                            .copied()
                            .filter(|&t| gate.admits(self.surface, t))
                            .collect();
                        if kept.len() != first.len() {
                            self.set_region(kept);
                            self.refit()?;
                        }
                        self.spread()?;
                    }
                    Err(OgeomError::Construction(_)) if i + 1 < patches.len() => {}
                    Err(e) => return Err(e),
                }
                continue;
            };
            if patch.iter().all(|&t| self.in_region[t as usize])
                || self.disagreement(&gate, patch) > 1.0
            {
                continue;
            }
            let joining: Vec<u32> = patch
                .iter()
                .copied()
                .filter(|&t| !self.in_region[t as usize] && gate.admits(self.surface, t))
                .collect();
            if joining.is_empty() {
                continue;
            }
            for t in joining {
                self.in_region[t as usize] = true;
                self.region.push(t);
            }
            self.refit()?;
            self.spread()?;
        }
        Ok(())
    }

    /// How far a neighbourhood is from agreeing with the fit, as the larger
    /// of its vertices' RMS residual over the gate, each residual shrunk by
    /// the widening the gate gives it, and its mean normal's angle over the
    /// limit: it agrees at one or less.
    fn disagreement(&self, gate: &Gate, patch: &[u32]) -> f64 {
        let surface = self.surface;
        let (_, points, _) = surface.vertices(patch);
        let squares: f64 = points
            .iter()
            .map(|p| gate.model.residual(*p).powi(2) / (1.0 + gate.model.leverage(*p)))
            .sum();
        #[allow(
            clippy::cast_precision_loss,
            reason = "point counts are far below 2^52"
        )]
        let rms = (squares / points.len().max(1) as f64).sqrt();
        let (mut cos, mut area) = (0.0, 0.0);
        for &t in patch {
            let i = t as usize;
            let a = surface.cross[i].magnitude();
            let c = surface.normal[i].dot(gate.model.normal(surface.centroid[i]));
            cos += a * if gate.model.oriented() { c } else { c.abs() };
            area += a;
        }
        let mean_angle = if area > 0.0 {
            (cos / area).clamp(-1.0, 1.0).acos()
        } else {
            core::f64::consts::PI
        };
        ratio(rms, gate.threshold).max(ratio(mean_angle, self.options.normal_limit))
    }

    fn capped(&self) -> bool {
        self.options
            .max_triangles
            .is_some_and(|most| self.region.len() >= most)
    }

    /// Grow across edges until a whole pass against the latest fit admits
    /// no triangle.
    fn spread(&mut self) -> OgeomResult<()> {
        let surface = self.surface;
        let count = self.in_region.len();
        let mut queued = vec![false; count];
        let mut refused = vec![false; count];
        let mut waiting: Vec<u32> = Vec::new();
        let mut queue: VecDeque<u32> = VecDeque::new();
        for &t in &self.region {
            for &n in &surface.across[t as usize] {
                if n != NONE && !self.in_region[n as usize] && !queued[n as usize] {
                    queued[n as usize] = true;
                    queue.push_back(n);
                }
            }
        }
        let Some(mut gate) = self.gate else {
            ogeom_bail!(Construction, "the seeds' neighbourhoods gave no fit");
        };
        let mut fitted_at = self.region.len();
        loop {
            while let Some(t) = queue.pop_front() {
                queued[t as usize] = false;
                if self.capped() {
                    queue.clear();
                    break;
                }
                if self.in_region[t as usize] {
                    continue;
                }
                if !gate.admits(surface, t) {
                    if !refused[t as usize] {
                        refused[t as usize] = true;
                        waiting.push(t);
                    }
                    continue;
                }
                self.in_region[t as usize] = true;
                self.region.push(t);
                for &n in &surface.across[t as usize] {
                    if n != NONE && !self.in_region[n as usize] && !queued[n as usize] {
                        queued[n as usize] = true;
                        queue.push_back(n);
                    }
                }
                if self.region.len() >= 2 * fitted_at {
                    gate = self.refit()?;
                    fitted_at = self.region.len();
                    requeue(&mut waiting, &mut refused, &mut queued, &mut queue);
                }
            }
            if self.region.len() == fitted_at {
                break;
            }
            gate = self.refit()?;
            fitted_at = self.region.len();
            if self.capped() {
                break;
            }
            requeue(&mut waiting, &mut refused, &mut queued, &mut queue);
            if queue.is_empty() {
                break;
            }
        }
        Ok(())
    }

    /// Take the rim off, ring by ring, while it turns away from the model
    /// by clearly more than the noise, keeping at least `least` vertices.
    fn peel(&mut self, least: usize) -> OgeomResult<u32> {
        let surface = self.surface;
        let mut rings = 0;
        let mut weight = vec![0.0_f64; surface.mesh.positions.len()];
        let mut marked = vec![false; surface.mesh.positions.len()];
        loop {
            let Some(Gate {
                model, threshold, ..
            }) = self.gate
            else {
                ogeom_bail!(Construction, "the seeds' neighbourhoods gave no fit");
            };
            self.region.sort_unstable();
            let mut inner = Vec::new();
            let mut touched: Vec<u32> = Vec::new();
            for &t in &self.region {
                let i = t as usize;
                let corners = surface.corners(t);
                let p = corners.map(|v| surface.position(v));
                let n = surface.normal[i];
                let mut out = Vector::ZERO;
                let mut on_rim = false;
                for k in 0..3 {
                    let other = surface.across[i][k];
                    if other != NONE && self.in_region[other as usize] {
                        continue;
                    }
                    on_rim = true;
                    let (a, b, c) = (p[k], p[(k + 1) % 3], p[(k + 2) % 3]);
                    let side = (b - a).cross(n);
                    out += if side.dot(c - a) > 0.0 { -side } else { side };
                }
                if !on_rim {
                    inner.push(t);
                    continue;
                }
                let length = out.magnitude();
                if length.is_nan() || length <= 0.0 {
                    continue;
                }
                let out = out / length;
                let doubled = surface.cross[i].magnitude();
                for k in 0..3 {
                    // The gradient of the barycentric weight of corner k.
                    let opposite = p[(k + 2) % 3] - p[(k + 1) % 3];
                    let gradient = n.cross(opposite) / doubled;
                    let v = corners[k];
                    if !marked[v as usize] {
                        marked[v as usize] = true;
                        touched.push(v);
                    }
                    weight[v as usize] += gradient.dot(out);
                }
            }
            let (mut sum, mut squares) = (0.0, 0.0);
            for &v in &touched {
                let c = weight[v as usize];
                sum += c * model.residual(surface.position(v));
                squares += c * c;
                weight[v as usize] = 0.0;
                marked[v as usize] = false;
            }
            let clear = sum.abs() > threshold * squares.sqrt();
            if !clear || inner.is_empty() || surface.vertices(&inner).0.len() < least {
                break;
            }
            self.set_region(inner);
            self.refit()?;
            rings += 1;
        }
        Ok(rings)
    }
}

fn requeue(
    waiting: &mut Vec<u32>,
    refused: &mut [bool],
    queued: &mut [bool],
    queue: &mut VecDeque<u32>,
) {
    for t in waiting.drain(..) {
        refused[t as usize] = false;
        if !queued[t as usize] {
            queued[t as usize] = true;
            queue.push_back(t);
        }
    }
}

/// `a` over `b`: zero for an `a` of zero or an infinite `b`, and
/// infinite for a positive `a` over a zero `b`.
fn ratio(a: f64, b: f64) -> f64 {
    if a <= 0.0 || b.is_infinite() {
        0.0
    } else if b > 0.0 {
        a / b
    } else {
        f64::INFINITY
    }
}

/// The kind's fit to the points, with their normals, from the last fit's
/// shape where there is one.
fn fitted(
    kind: FitKind,
    points: &[Point],
    normals: Option<&[Vector]>,
    options: &GrowOptions,
    last: Option<&AnyFit>,
    tol: Tolerances,
) -> OgeomResult<AnyFit> {
    fn with<'n, S>(
        options: &GrowOptions,
        normals: Option<&'n [Vector]>,
        guess: Option<S>,
        robust: bool,
    ) -> FitOptions<'n, S> {
        FitOptions {
            normals,
            cutoff: options.cutoff,
            robust_start: robust && guess.is_none() && options.cutoff != Cutoff::All,
            fixed: options.fixed,
            guess,
            side: options.side,
            radius_range: None,
        }
    }
    // A consensus start for a curved kind, whose start from a
    // neighbourhood a plane grew may hold another surface's points; a
    // plane's and a line's closed forms start from nothing.
    let robust = !matches!(kind, FitKind::Plane | FitKind::Line3);
    Ok(match kind {
        FitKind::Plane => {
            let guess = match last {
                Some(AnyFit::Plane(f)) => Some(f.shape),
                _ => None,
            };
            AnyFit::Plane(primitive::plane(
                points,
                &with(options, normals, guess, robust),
                tol,
            )?)
        }
        FitKind::Sphere => {
            let guess = match last {
                Some(AnyFit::Sphere(f)) => Some(f.shape),
                _ => None,
            };
            AnyFit::Sphere(primitive::sphere(
                points,
                &with(options, normals, guess, robust),
                tol,
            )?)
        }
        FitKind::Cylinder => {
            let guess = match last {
                Some(AnyFit::Cylinder(f)) => Some(f.shape),
                _ => None,
            };
            AnyFit::Cylinder(primitive::cylinder(
                points,
                &with(options, normals, guess, robust),
                tol,
            )?)
        }
        FitKind::Cone => {
            let guess = match last {
                Some(AnyFit::Cone(f)) => Some(f.shape),
                _ => None,
            };
            AnyFit::Cone(primitive::cone(
                points,
                &with(options, normals, guess, robust),
                tol,
            )?)
        }
        FitKind::Torus => {
            let guess = match last {
                Some(AnyFit::Torus(f)) => Some(f.shape),
                _ => None,
            };
            AnyFit::Torus(primitive::torus(
                points,
                &with(options, normals, guess, robust),
                tol,
            )?)
        }
        FitKind::Circle3 => {
            let guess = match last {
                Some(AnyFit::Circle(f)) => Some(f.shape),
                _ => None,
            };
            AnyFit::Circle(primitive::circle3(
                points,
                &with(options, normals, guess, robust),
                tol,
            )?)
        }
        FitKind::Line3 => AnyFit::Line(primitive::line3(
            points,
            &with(options, normals, None, robust),
            tol,
        )?),
    })
}
