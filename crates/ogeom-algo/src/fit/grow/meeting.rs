//! A crease line where two planes meet either side of a rounded edge, for
//! [`grow`](super::grow): the sharp edge the round stands in for.

use ogeom_core::{OgeomError, OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Axis, Direction, Plane, Point, Vector};
use ogeom_mesh::Bvh;

use super::{
    AnyFit, Crease, FitKind, Gate, GrowOptions, Grown, Meeting, NONE, SeedReport, Surface, region,
};
use crate::Cutoff;
use crate::fit::primitive::{self, Fit, FitFixed, Kind};

/// The line where two planes meet across the round the first seed lies
/// on, as [`grow`](super::grow) describes under its rounded edges.
pub(super) fn met(
    surface: &Surface<'_>,
    index: &Bvh,
    seeds: &[Point],
    options: &GrowOptions,
    tol: Tolerances,
) -> OgeomResult<Grown> {
    if options.cutoff == Cutoff::All {
        ogeom_bail!(
            Construction,
            "the planes beside a round are peeled, and Cutoff::All has no threshold to peel by"
        );
    }
    if options.fixed != FitFixed::default() {
        ogeom_bail!(
            Construction,
            "two planes' meeting cannot hold a line's fixed parts"
        );
    }
    let mut hits = Vec::with_capacity(seeds.len());
    for (i, seed) in seeds.iter().enumerate() {
        let Some(hit) = index.nearest(*seed, f64::INFINITY) else {
            ogeom_bail!(Construction, "no triangle of the mesh lies near seed {i}");
        };
        hits.push(hit);
    }
    let start = hits[0].triangle;
    let Some(axis) = turning_axis(surface, start, options.crease_limit) else {
        ogeom_bail!(
            Construction,
            "the surface about seed 0 turns by no more than the crease limit {} one way more \
             than another",
            options.crease_limit
        );
    };
    let section = Section {
        surface,
        origin: hits[0].point,
        normal: axis,
    };
    let walks = Walks::of(&section, start);
    let peeled = GrowOptions {
        peel: true,
        fixed: FitFixed::default(),
        side: None,
        ..*options
    };

    // Each side's planes, grown from the triangles the section reaches 1,
    // 2, 4, ... triangles out; a triangle a plane already holds grows no
    // other, and a plane holding the seed's triangle is no side of a round
    // about it.
    let mut sides: [Vec<Grown>; 2] = [Vec::new(), Vec::new()];
    let mut refusal = None;
    for (j, side) in sides.iter_mut().enumerate() {
        let mut step = 1;
        while step <= walks.steps[j].len() {
            let t = walks.steps[j][step - 1].0;
            step *= 2;
            if side
                .iter()
                .any(|plane| plane.triangles.binary_search(&t).is_ok())
            {
                continue;
            }
            let from = surface.centroid[t as usize];
            match region(surface, index, &[from], FitKind::Plane, &peeled, tol) {
                Ok(plane) if plane.triangles.binary_search(&start).is_err() => side.push(plane),
                Ok(_) => {}
                Err(e @ OgeomError::Construction(_)) => refusal = Some(e),
                Err(e) => return Err(e),
            }
        }
    }

    // Of the pairs that meet across a round, the one whose planes are
    // widest for the way across it: the faces either side of the round
    // rather than a facet of the round beside one of them, or faces
    // further off across others.
    let mut best: Option<Found> = None;
    for first in &sides[0] {
        for second in &sides[1] {
            let planes = [first, second];
            match across(surface, &section, &walks, start, planes, &peeled, tol) {
                Ok(found) => {
                    if best.as_ref().is_none_or(|b| found.clear > b.clear) {
                        best = Some(found);
                    }
                }
                Err(e @ OgeomError::Construction(_)) => refusal = Some(e),
                Err(e) => return Err(e),
            }
        }
    }
    if let Some(found) = best {
        return Ok(found.answer(&hits, seeds));
    }
    Err(refusal.unwrap_or_else(|| {
        ogeom_core::ogeom_err!(
            Construction,
            "the section across the surface through seed 0 reaches no plane on either side"
        )
    }))
}

/// The section's walks out of the seed's triangle both ways, each with the
/// way along it from the seed's point to where it enters each triangle.
struct Walks {
    steps: [Vec<(u32, Point)>; 2],
    length: [Vec<f64>; 2],
}

impl Walks {
    fn of(section: &Section<'_, '_>, start: u32) -> Self {
        let steps = [0, 1].map(|way| section.walk(start, way));
        let length = [0, 1].map(|j| {
            let mut from = section.origin;
            let mut sum = 0.0;
            steps[j]
                .iter()
                .map(|(_, p)| {
                    sum += (*p - from).magnitude();
                    from = *p;
                    sum
                })
                .collect()
        });
        Self { steps, length }
    }

    /// The first step of walk `j` into a triangle `held` holds.
    fn entry(&self, j: usize, held: impl Fn(u32) -> bool) -> Option<usize> {
        self.steps[j].iter().position(|(t, _)| held(*t))
    }
}

/// The direction the surface about `start` turns least, the axis of the
/// round it lies on: the rings about the triangle grow until their
/// normals, weighted by area, spread one way across their mean by more
/// than they do the other by at least as much as normals spread evenly
/// over the crease limit do; the axis is the way they spread least.
/// `None` where they never do.
fn turning_axis(surface: &Surface<'_>, start: u32, limit: f64) -> Option<Vector> {
    // Normals spread evenly over an arc of angle `L` have a second moment
    // across their mean of `(1 - sin L / L) / 2`.
    let wanted = 0.5 * (1.0 - limit.sin() / limit);
    let count = surface.mesh.triangles.len();
    let mut taken = vec![false; count];
    let mut ring = vec![start];
    taken[start as usize] = true;
    let mut moment = nalgebra::Matrix3::<f64>::zeros();
    let mut area = 0.0;
    while !ring.is_empty() {
        for &t in &ring {
            let a = surface.cross[t as usize].magnitude();
            let n = surface.normal[t as usize];
            let n = nalgebra::Vector3::new(n.x, n.y, n.z);
            moment += n * n.transpose() * a;
            area += a;
        }
        if area > 0.0 {
            let eigen = nalgebra::SymmetricEigen::new(moment / area);
            let mut order = [0, 1, 2];
            order.sort_by(|&i, &j| eigen.eigenvalues[j].total_cmp(&eigen.eigenvalues[i]));
            let values = order.map(|k| eigen.eigenvalues[k]);
            if values[1] - values[2] >= wanted {
                let v = eigen.eigenvectors.column(order[2]);
                return Some(Vector::new(v[0], v[1], v[2]));
            }
        }
        let mut next: Vec<u32> = ring
            .iter()
            .flat_map(|&t| surface.corners(t))
            .flat_map(|v| surface.fan(v).iter().copied())
            .filter(|&t| !taken[t as usize])
            .collect();
        next.sort_unstable();
        next.dedup();
        for &t in &next {
            taken[t as usize] = true;
        }
        ring = next;
    }
    None
}

/// Whether a plane's points lie on it to within the noise the mesh shows
/// there, and not on a curved patch the plane only approximates.
///
/// The noise is measured at the plane's vertices whose every triangle is
/// the plane's, each by its offset along the plane's normal from the mean
/// of its `m` neighbours: over a flat with noise of sigma `σ` that varies
/// as `σ² (1 + 1/m)`; over a patch of constant curvature it moves by a
/// constant as well, which taking the offsets about their mean leaves out,
/// while the plane's own residuals keep the patch's bend. A plane is flat
/// where its sigma squared over the noise's is no larger than the
/// cutoff's `k` sigmas of chance allow: each estimate's logarithm varies as
/// `2/ν` over its `ν` degrees of freedom, so the ratio is at most
/// `exp(k √(2/ν₁ + 2/ν₂))`. With [`Cutoff::Distance`] the sigma is at most
/// the distance. Neither is held below what rounding to single precision
/// or the confusion distance leaves. A plane with fewer than two vertices
/// inside it shows no noise, and is not flat.
fn flat(
    surface: &Surface<'_>,
    plane: &Grown,
    fit: &Fit<Plane>,
    gate: &Gate,
    inside: impl Fn(u32) -> bool,
    cutoff: Cutoff,
    tol: Tolerances,
) -> OgeomResult<()> {
    let floor = (gate.sigmas * gate.quantum).max(tol.confusion());
    let k = match cutoff {
        Cutoff::Distance(d) => {
            if fit.sigma > d.max(floor) {
                ogeom_bail!(
                    Construction,
                    "a plane beside the round misses its points by a sigma of {}, past the \
                     cutoff's distance {d}",
                    fit.sigma
                );
            }
            return Ok(());
        }
        Cutoff::Sigma(k) | Cutoff::Median(k) => k,
        Cutoff::All => return Ok(()),
    };
    let normal = fit.shape.normal().vector();
    let mut offsets: Vec<(f64, f64)> = Vec::new();
    let mut around: Vec<u32> = Vec::new();
    for &v in &plane.vertices {
        let fan = surface.fan(v);
        if fan.is_empty() || !fan.iter().all(|&t| inside(t)) {
            continue;
        }
        around.clear();
        around.extend(
            fan.iter()
                .flat_map(|&t| surface.corners(t))
                .filter(|&w| w != v),
        );
        around.sort_unstable();
        around.dedup();
        #[allow(
            clippy::cast_precision_loss,
            reason = "a vertex has far fewer than 2^52 neighbours"
        )]
        let m = around.len() as f64;
        let mean = around.iter().fold(Vector::ZERO, |sum, &w| {
            sum + (surface.position(w) - Point::ORIGIN)
        }) / m;
        let off = (surface.position(v) - Point::ORIGIN - mean).dot(normal);
        offsets.push((off, 1.0 + 1.0 / m));
    }
    let count = offsets.len();
    if count < 2 {
        ogeom_bail!(
            Construction,
            "a plane beside the round has {count} vertices inside it, too few to show its noise"
        );
    }
    #[allow(
        clippy::cast_precision_loss,
        reason = "point counts are far below 2^52"
    )]
    let (n, kept) = (count as f64, fit.used.iter().filter(|u| **u).count() as f64);
    let middle = offsets.iter().map(|(o, _)| o).sum::<f64>() / n;
    let spread: f64 = offsets.iter().map(|(o, _)| (o - middle).powi(2)).sum();
    let weight: f64 = offsets.iter().map(|(_, w)| w).sum::<f64>() * (n - 1.0) / n;
    let noise = (spread / weight).sqrt().max(floor);
    let fitted = (kept - 3.0).max(1.0);
    let allowed = (k * (2.0 / fitted + 2.0 / (n - 1.0)).sqrt()).exp();
    let sigma = fit.sigma.max(floor);
    if sigma * sigma > noise * noise * allowed {
        ogeom_bail!(
            Construction,
            "a plane beside the round misses its points by a sigma of {}, past the noise {} \
             its mesh shows inside it by more than chance allows",
            fit.sigma,
            noise
        );
    }
    Ok(())
}

/// Two planes' meeting across a round, as measured.
struct Found {
    planes: [Grown; 2],
    fits: [Fit<Plane>; 2],
    angle: f64,
    apex: Point,
    along: Direction,
    span: (f64, f64),
    width: f64,
    /// How many times the way across the round the narrower plane is wide.
    clear: f64,
}

/// Whether the planes on either side of `start` meet across a round, as
/// [`grow`](super::grow) gives the tests under its rounded edges.
fn across(
    surface: &Surface<'_>,
    section: &Section<'_, '_>,
    walks: &Walks,
    start: u32,
    planes: [&Grown; 2],
    peeled: &GrowOptions,
    tol: Tolerances,
) -> OgeomResult<Found> {
    let count = surface.mesh.triangles.len();
    let mut in_plane = vec![[false; 2]; count];
    let least = 2 * primitive::needed(Kind::Plane, &FitFixed::default());
    let mut gates = Vec::with_capacity(2);
    let mut fits = Vec::with_capacity(2);
    for (j, plane) in planes.iter().enumerate() {
        let AnyFit::Plane(fit) = &plane.fit else {
            ogeom_bail!(Construction, "a side gave no plane");
        };
        for &t in &plane.triangles {
            in_plane[t as usize][j] = true;
        }
        if fit.used.iter().filter(|u| **u).count() < least {
            ogeom_bail!(
                Construction,
                "a plane beside the round keeps fewer than {least} points"
            );
        }
        let points: Vec<Point> = plane
            .vertices
            .iter()
            .map(|&v| surface.position(v))
            .collect();
        gates.push(Gate::of(&plane.fit, &points, surface.single, peeled, tol));
        fits.push(fit.clone());
    }
    let Ok(fits) = <[Fit<Plane>; 2]>::try_from(fits) else {
        ogeom_bail!(Construction, "two planes were not fitted");
    };
    if in_plane[start as usize].iter().any(|on| *on) {
        ogeom_bail!(
            Construction,
            "seed 0 lies on a plane beside it, not on a round"
        );
    }
    if in_plane.iter().any(|on| on[0] && on[1]) {
        ogeom_bail!(
            Construction,
            "the two sides of the round grew into one plane"
        );
    }
    for (j, plane) in planes.iter().enumerate() {
        let inside = |t: u32| in_plane[t as usize][j];
        flat(
            surface,
            plane,
            &fits[j],
            &gates[j],
            inside,
            peeled.cutoff,
            tol,
        )?;
    }

    // Where they meet: the line's direction, and its point nearest where
    // the section starts.
    let n = [0, 1].map(|j| fits[j].shape.normal().vector());
    let cos = n[0].dot(n[1]).clamp(-1.0, 1.0);
    let angle = cos.acos();
    let Ok(along) = Direction::new(n[0].cross(n[1]), tol) else {
        ogeom_bail!(Construction, "the planes beside the round are parallel");
    };
    let d = along.vector();

    // They turn past the crease limit beyond the doubt their own
    // uncertainty puts on the angle: a plane's residual is quadratic in the
    // point, so its leverage's second difference across the line, a unit
    // either side of its centroid, is its tilt there over the noise.
    let tilts = [0, 1].map(|j| {
        let kept: Vec<Point> = planes[j]
            .vertices
            .iter()
            .zip(&fits[j].used)
            .filter(|(_, keep)| **keep)
            .map(|(&v, _)| surface.position(v))
            .collect();
        #[allow(
            clippy::cast_precision_loss,
            reason = "point counts are far below 2^52"
        )]
        let centre = Point::ORIGIN
            + kept
                .iter()
                .fold(Vector::ZERO, |sum, p| sum + (*p - Point::ORIGIN))
                / kept.len() as f64;
        let model = &gates[j].model;
        let across = d.cross(n[j]);
        let bend = 0.5 * (model.leverage(centre + across) + model.leverage(centre - across))
            - model.leverage(centre);
        gates[j].threshold * bend.max(0.0).sqrt()
    });
    let unsure = tilts[0].hypot(tilts[1]);
    if angle - unsure <= peeled.crease_limit {
        ogeom_bail!(
            Construction,
            "the planes beside the round turn by {angle}, give or take {unsure}, not past the \
             crease limit {}",
            peeled.crease_limit
        );
    }
    let at = section.origin;
    let short = [0, 1].map(|j| -fits[j].shape.signed_distance_to(at));
    let det = 1.0 - cos * cos;
    let apex = at
        + n[0] * ((short[0] - cos * short[1]) / det)
        + n[1] * ((short[1] - cos * short[0]) / det);

    // Where the section through the seed enters each plane: the round's
    // ends there.
    let mut ends = [Point::ORIGIN; 2];
    let mut way = 0.0;
    for (j, end) in ends.iter_mut().enumerate() {
        let Some(at) = walks.entry(j, |t| in_plane[t as usize][0] || in_plane[t as usize][1])
        else {
            ogeom_bail!(
                Construction,
                "the section through seed 0 does not reach the plane on its side"
            );
        };
        if !in_plane[walks.steps[j][at].0 as usize][j] {
            ogeom_bail!(
                Construction,
                "the section through seed 0 reaches the other side's plane first"
            );
        }
        *end = walks.steps[j][at].1;
        way += walks.length[j][at];
    }
    let chord = ends[1] - ends[0];
    let width = chord.magnitude();
    let Ok(chord_way) = Direction::new(chord, tol) else {
        ogeom_bail!(Construction, "the round has no width");
    };
    let chord_way = chord_way.vector();

    // The apex lies over the round, beyond the doubt the planes' own
    // uncertainty puts on it: moving plane `j` by `δ` along its normal
    // moves the apex by `δ m_j`, in the other plane and across the line.
    let doubt = [0, 1]
        .map(|j| {
            let m = n[1 - j].cross(d);
            let moved = m / m.dot(n[j]);
            let sway = gates[j].threshold * gates[j].model.leverage(apex).sqrt();
            (sway * moved.dot(chord_way)).powi(2)
        })
        .iter()
        .sum::<f64>()
        .sqrt();
    let over = (apex - ends[0]).dot(chord_way);
    if !(over >= doubt && width - over >= doubt) {
        ogeom_bail!(
            Construction,
            "the planes beside the round meet {over} along its width {width}, not over it \
             beyond their doubt {doubt}"
        );
    }

    // The way across the round along the section is shorter than each
    // plane is wide across the line, which a facet of the round, the round
    // turning on past it, is not; and the round runs where both planes run
    // beside it.
    let (mut first, mut last) = (f64::NEG_INFINITY, f64::INFINITY);
    let mut narrowest = f64::INFINITY;
    for j in 0..2 {
        let kept: Vec<Point> = planes[j]
            .vertices
            .iter()
            .zip(&fits[j].used)
            .filter(|(_, keep)| **keep)
            .map(|(&v, _)| surface.position(v))
            .collect();
        let spread = |way: Vector| {
            kept.iter()
                .map(|p| (*p - apex).dot(way))
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), x| {
                    (lo.min(x), hi.max(x))
                })
        };
        let (lo, hi) = spread(d.cross(n[j]));
        if way >= hi - lo {
            ogeom_bail!(
                Construction,
                "the way across the round, {way}, is not shorter than a plane beside it is \
                 wide, {}",
                hi - lo
            );
        }
        narrowest = narrowest.min(hi - lo);
        let (lo, hi) = spread(d);
        first = first.max(lo);
        last = last.min(hi);
    }
    if first >= last {
        ogeom_bail!(
            Construction,
            "the planes beside the round do not run beside each other along it"
        );
    }
    Ok(Found {
        planes: planes.map(Clone::clone),
        fits,
        angle,
        apex,
        along,
        span: (first, last),
        width,
        clear: narrowest / way,
    })
}

impl Found {
    /// The answer: the line, its residuals the planes', and each seed's
    /// report against it.
    fn answer(self, hits: &[ogeom_mesh::Hit], seeds: &[Point]) -> Grown {
        let middle = 0.5 * (self.span.0 + self.span.1);
        let d = self.along.vector();
        let line = Axis::new(self.apex + d * middle, self.along);
        let mut vertices = Vec::new();
        let mut residuals = Vec::new();
        let mut used = Vec::new();
        for (plane, fit) in self.planes.iter().zip(&self.fits) {
            vertices.extend_from_slice(&plane.vertices);
            residuals.extend_from_slice(&fit.residuals);
            used.extend_from_slice(&fit.used);
        }
        let kept: Vec<f64> = residuals
            .iter()
            .zip(&used)
            .filter_map(|(r, keep)| keep.then_some(*r))
            .collect();
        #[allow(
            clippy::cast_precision_loss,
            reason = "point counts are far below 2^52"
        )]
        let sigma = (kept.iter().map(|r| r * r).sum::<f64>() / kept.len() as f64).sqrt();
        let form = kept
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), r| {
                (lo.min(*r), hi.max(*r))
            });
        let [a, b] = &self.fits;
        let fit = Fit {
            shape: line,
            sigma,
            form,
            used,
            residuals,
            iterations: a.iterations + b.iterations,
            rounds: a.rounds.max(b.rounds),
            converged: a.converged && b.converged,
            side: None,
        };
        let model = AnyFit::Line(fit.clone()).model();
        let surface_of = |t: u32| {
            (
                self.planes[0].triangles.binary_search(&t).is_ok(),
                self.planes[1].triangles.binary_search(&t).is_ok(),
            )
        };
        let reports = seeds
            .iter()
            .zip(hits)
            .map(|(seed, hit)| {
                let (on_first, on_second) = surface_of(hit.triangle);
                let along = (hit.point - line.location).dot(d);
                SeedReport {
                    triangle: hit.triangle,
                    in_region: !on_first
                        && !on_second
                        && along >= self.span.0 - middle
                        && along <= self.span.1 - middle
                        && model.residual(hit.point) <= self.width,
                    offset: model.residual(*seed),
                    angle: self.angle,
                }
            })
            .collect();
        let [first, second] = self.planes.map(|p| p.triangles);
        Grown {
            triangles: Vec::new(),
            fit: AnyFit::Line(fit),
            vertices,
            seeds: reports,
            peeled_rings: 0,
            crease: Some(Crease::Meeting(Box::new(Meeting {
                planes: self.fits,
                triangles: [first, second],
                angle: self.angle,
                span: (self.span.0 - middle, self.span.1 - middle),
            }))),
        }
    }
}

/// The mesh cut by the plane through `origin` across `normal`: a walk
/// from triangle to triangle along the cut.
struct Section<'s, 'a> {
    surface: &'s Surface<'a>,
    origin: Point,
    normal: Vector,
}

impl Section<'_, '_> {
    /// Whether the vertex lies on the plane's positive side, zero counted
    /// as positive, so each edge is cut or not without doubt.
    fn above(&self, v: u32) -> bool {
        (self.surface.position(v) - self.origin).dot(self.normal) >= 0.0
    }

    /// The edges `k` of the triangle, `(k, k + 1)`, the plane cuts: none or
    /// two.
    fn cut(&self, t: u32) -> Vec<usize> {
        let c = self.surface.corners(t);
        (0..3)
            .filter(|&k| self.above(c[k]) != self.above(c[(k + 1) % 3]))
            .collect()
    }

    /// Where the plane cuts the edge `a b`.
    fn crossing(&self, a: u32, b: u32) -> Point {
        let (p, q) = (self.surface.position(a), self.surface.position(b));
        let (s, e) = (
            (p - self.origin).dot(self.normal),
            (q - self.origin).dot(self.normal),
        );
        p + (q - p) * (s / (s - e))
    }

    /// Walk from `start` out through the `way`-th edge the plane cuts in
    /// it: each triangle reached in turn, with where the walk entered it,
    /// until the walk meets the mesh's border, an edge three triangles
    /// share or a triangle the plane does not cross, or comes round to
    /// `start`.
    fn walk(&self, start: u32, way: usize) -> Vec<(u32, Point)> {
        let surface = self.surface;
        let mut reached = Vec::new();
        let cut = self.cut(start);
        if cut.len() != 2 {
            return reached;
        }
        let (mut t, mut k) = (start, cut[way]);
        for _ in 0..surface.mesh.triangles.len() {
            let c = surface.corners(t);
            let (a, b) = (c[k], c[(k + 1) % 3]);
            let next = surface.across[t as usize][k];
            if next == NONE || next == start {
                break;
            }
            reached.push((next, self.crossing(a, b)));
            let d = surface.corners(next);
            let Some(leaving) = self.cut(next).into_iter().find(|&j| {
                let (x, y) = (d[j], d[(j + 1) % 3]);
                !((x == a && y == b) || (x == b && y == a))
            }) else {
                break;
            };
            t = next;
            k = leaving;
        }
        reached
    }
}
