//! A circle or a line traced along a mesh's creases from a seed, for
//! [`grow`](super::grow).

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use ogeom_core::{OgeomError, OgeomResult, Tolerances, ogeom_bail};
use ogeom_math::{Point, Vector};
use ogeom_mesh::Bvh;

use super::{
    AnyFit, Crease, FitKind, Gate, GrowOptions, Grown, NONE, SeedReport, Surface, fitted, meeting,
};
use crate::Cutoff;
use crate::fit::{primitive, tails};

/// Trace the crease the first seed is nearest and fit `kind` along it.
///
/// The first seed's crease edges within reach of it are tried nearest
/// first, the reach being as many of its triangle's longest edge as the
/// shortest chain accepted has edges. The first whose chain holds twice
/// the vertices the fit needs, and turns past the limit by more than the
/// noise its fit measures could turn it, is traced: an edge the noise
/// turned past the limit on its own, or with a few neighbours, makes a
/// chain too short to measure the kind over or too little past the limit,
/// and is passed over. Where none does, a line is sought where two planes
/// meet across a round, and failing that the nearest edge's fit's own
/// refusal is given, or where no fit refused, the seed is refused as lying
/// on a smooth part.
pub(super) fn traced(
    surface: &Surface<'_>,
    index: &Bvh,
    seeds: &[Point],
    kind: FitKind,
    options: &GrowOptions,
    tol: Tolerances,
) -> OgeomResult<Grown> {
    let creases = Creases::of(surface, options.crease_limit);
    let least = 2 * primitive::needed(kind.kind(), &options.fixed);
    let smooth = |i: usize| -> OgeomError {
        ogeom_core::ogeom_err!(
            Construction,
            "seed {i} lies on a smooth part of the mesh: no chain of {least} vertices near \
             it, enough to measure a {} over, runs along edges turning by more than the \
             crease limit {}, past it by more than the noise",
            kind.kind().name(),
            options.crease_limit
        )
    };
    let mut reached = Vec::with_capacity(seeds.len());
    for (i, seed) in seeds.iter().enumerate() {
        let Some(hit) = index.nearest(*seed, f64::INFINITY) else {
            ogeom_bail!(Construction, "no triangle of the mesh lies near seed {i}");
        };
        let edges = creases.about(surface, hit.triangle, *seed);
        if edges.is_empty() {
            return if i == 0 {
                meeting_or(surface, index, seeds, kind, options, tol, smooth(0))
            } else {
                Err(smooth(i))
            };
        }
        reached.push((hit.triangle, edges));
    }
    let corners = surface.corners(reached[0].0).map(|v| surface.position(v));
    let longest = (0..3)
        .map(|k| (corners[(k + 1) % 3] - corners[k]).magnitude())
        .fold(0.0, f64::max);
    #[allow(
        clippy::cast_precision_loss,
        reason = "a chain needs a handful of vertices"
    )]
    let reach = (least - 1) as f64 * longest;
    let within = reached[0].1.partition_point(|edge| edge.3 <= reach);

    let mut chain = Chain {
        surface,
        creases: &creases,
        kind,
        options,
        tol,
        vertices: VecDeque::new(),
        on_chain: vec![false; surface.mesh.positions.len()],
        closed: false,
        fit: None,
        gate: None,
    };
    let mut tried: BTreeSet<(u32, u32)> = BTreeSet::new();
    let mut found = false;
    let mut refusal = None;
    for &(a, b, ..) in &reached[0].1[..within] {
        if tried.contains(&(a, b)) {
            continue;
        }
        chain.restart(a, b);
        let outcome = chain.trace();
        if chain.vertices.len() >= least {
            outcome?;
            if chain.beyond_noise(within) {
                found = true;
                break;
            }
        } else {
            match outcome {
                Ok(()) => {}
                Err(e @ OgeomError::Construction(_)) => {
                    refusal.get_or_insert(e);
                }
                Err(e) => return Err(e),
            }
        }
        let mut walked: Vec<u32> = chain.vertices.iter().copied().collect();
        if chain.closed {
            walked.push(walked[0]);
        }
        for pair in walked.windows(2) {
            tried.insert((pair[0].min(pair[1]), pair[0].max(pair[1])));
        }
    }
    if !found {
        let refused = refusal.unwrap_or_else(|| smooth(0));
        return meeting_or(surface, index, seeds, kind, options, tol, refused);
    }

    let Chain {
        vertices,
        on_chain,
        fit,
        gate,
        ..
    } = chain;
    let (Some(fit), Some(Gate { model, .. })) = (fit, gate) else {
        ogeom_bail!(Construction, "the crease gave no fit");
    };
    let reports = seeds
        .iter()
        .zip(&reached)
        .map(|(seed, (triangle, edges))| {
            let on = edges
                .iter()
                .find(|(a, b, ..)| on_chain[*a as usize] && on_chain[*b as usize]);
            SeedReport {
                triangle: *triangle,
                in_region: on.is_some(),
                offset: model.residual(*seed),
                angle: on.unwrap_or(&edges[0]).2,
            }
        })
        .collect();
    Ok(Grown {
        triangles: Vec::new(),
        fit,
        vertices: vertices.into(),
        seeds: reports,
        peeled_rings: 0,
        crease: Some(Crease::Traced),
    })
}

/// For a line, where two planes meet across a round about the first seed;
/// for a circle, or where they do not, `refused`.
fn meeting_or(
    surface: &Surface<'_>,
    index: &Bvh,
    seeds: &[Point],
    kind: FitKind,
    options: &GrowOptions,
    tol: Tolerances,
    refused: OgeomError,
) -> OgeomResult<Grown> {
    if kind != FitKind::Line3 {
        return Err(refused);
    }
    match meeting::met(surface, index, seeds, options, tol) {
        Err(OgeomError::Construction(_)) => Err(refused),
        met => met,
    }
}

/// A mesh's crease edges: which edges of each triangle are creases, and
/// the vertices each vertex meets by a crease, in increasing order.
struct Creases {
    /// Whether the edge `(k, k + 1)` of each triangle is a crease.
    edge: Vec<[bool; 3]>,
    /// The angle the surface turns through across each such edge.
    turn: Vec<[f64; 3]>,
    start: Vec<usize>,
    next: Vec<u32>,
}

impl Creases {
    /// The edges two triangles share across which their normals turn by
    /// more than `limit`, each triangle's normal read with its neighbour's
    /// winding where the two disagree about the edge's direction.
    fn of(surface: &Surface<'_>, limit: f64) -> Self {
        let triangles = &surface.mesh.triangles;
        let cos_limit = limit.cos();
        let mut edge = vec![[false; 3]; triangles.len()];
        let mut turn = vec![[0.0; 3]; triangles.len()];
        let mut pairs: Vec<(u32, u32)> = Vec::new();
        for (t, tri) in triangles.iter().enumerate() {
            let n = surface.normal[t];
            for k in 0..3 {
                let other = surface.across[t][k];
                if other == NONE || n == Vector::ZERO {
                    continue;
                }
                let m = surface.normal[other as usize];
                if m == Vector::ZERO {
                    continue;
                }
                let (a, b) = (tri[k], tri[(k + 1) % 3]);
                let o = triangles[other as usize];
                let same = (0..3).any(|j| o[j] == a && o[(j + 1) % 3] == b);
                let cos = if same { -n.dot(m) } else { n.dot(m) };
                turn[t][k] = cos.clamp(-1.0, 1.0).acos();
                if cos < cos_limit {
                    edge[t][k] = true;
                    pairs.push((a, b));
                    pairs.push((b, a));
                }
            }
        }
        pairs.sort_unstable();
        pairs.dedup();
        let count = surface.mesh.positions.len();
        let mut start = vec![0_usize; count + 1];
        for &(a, _) in &pairs {
            start[a as usize + 1] += 1;
        }
        for v in 0..count {
            start[v + 1] += start[v];
        }
        let next = pairs.into_iter().map(|(_, b)| b).collect();
        Self {
            edge,
            turn,
            start,
            next,
        }
    }

    /// The vertices `v` meets by a crease.
    fn from(&self, v: u32) -> &[u32] {
        let v = v as usize;
        &self.next[self.start[v]..self.start[v + 1]]
    }

    /// The crease edges bounding the smooth part of the surface `triangle`
    /// lies on, the triangles reached from it across edges that are not
    /// creases, each with the angle the surface turns through there and
    /// its distance from `seed`: nearest first, the lesser vertex pair
    /// first on a tie.
    fn about(
        &self,
        surface: &Surface<'_>,
        triangle: u32,
        seed: Point,
    ) -> Vec<(u32, u32, f64, f64)> {
        let mut seen = vec![false; self.edge.len()];
        let mut queue = VecDeque::from([triangle]);
        seen[triangle as usize] = true;
        let mut found: Vec<(f64, (u32, u32), f64)> = Vec::new();
        while let Some(t) = queue.pop_front() {
            let i = t as usize;
            let tri = surface.corners(t);
            for k in 0..3 {
                if self.edge[i][k] {
                    let (a, b) = (tri[k].min(tri[(k + 1) % 3]), tri[k].max(tri[(k + 1) % 3]));
                    let d = to_segment(seed, surface.position(a), surface.position(b));
                    found.push((d, (a, b), self.turn[i][k]));
                    continue;
                }
                let other = surface.across[i][k];
                if other != NONE && !seen[other as usize] {
                    seen[other as usize] = true;
                    queue.push_back(other);
                }
            }
        }
        found.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        found.dedup_by_key(|x| x.1);
        found
            .into_iter()
            .map(|(d, (a, b), turn)| (a, b, turn, d))
            .collect()
    }
}

/// The angle the surface turns through across the edge `a b`, between the
/// triangles `a b c` and `a b d`, and its gradient in `a`, `b`, `c` and
/// `d`; `None` where a triangle has no area.
///
/// The angle is π less the one between the half-planes, each from the
/// edge out to `c` or `d`. Moving `c` across its half-plane by `δ` turns
/// that half-plane about the edge by `δ / h_c`, `h_c` the distance of `c`
/// from the edge; moving `a` or `b` turns it by the share of `δ` the
/// triangle's plane takes at the edge's other end, `1 - s` or `s`, `s` how
/// far along the edge `c`'s foot lies.
fn fold([a, b, c, d]: [Point; 4]) -> Option<(f64, [Vector; 4])> {
    let edge = b - a;
    let length = edge.dot(edge);
    if length.is_nan() || length <= 0.0 {
        return None;
    }
    let along = edge / length.sqrt();
    let foot = |q: Point| {
        let s = (q - a).dot(edge) / length;
        let out = (q - a) - edge * s;
        (s, out, out.magnitude())
    };
    let ((sc, uc, hc), (sd, ud, hd)) = (foot(c), foot(d));
    if !(hc > 0.0 && hd > 0.0) {
        return None;
    }
    let (u1, u2) = (uc / hc, ud / hd);
    let sine = u1.cross(u2).dot(along);
    let turn = core::f64::consts::PI - sine.abs().atan2(u1.dot(u2));
    let way = if sine < 0.0 { -1.0 } else { 1.0 };
    let (r1, r2) = (along.cross(u1) * (way / hc), along.cross(u2) * (way / hd));
    Some((
        turn,
        [
            r2 * (1.0 - sd) - r1 * (1.0 - sc),
            r2 * sd - r1 * sc,
            r1,
            -r2,
        ],
    ))
}

/// The distance from `p` to the segment `a b`.
fn to_segment(p: Point, a: Point, b: Point) -> f64 {
    let ab = b - a;
    let length = ab.dot(ab);
    let s = if length > 0.0 {
        ((p - a).dot(ab) / length).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (p - (a + ab * s)).magnitude()
}

/// One end of the chain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum End {
    Front,
    Back,
}

/// A chain of crease vertices as it is traced, and its latest fit.
struct Chain<'s, 'a> {
    surface: &'s Surface<'a>,
    creases: &'s Creases,
    kind: FitKind,
    options: &'s GrowOptions,
    tol: Tolerances,
    /// In their order along the crease.
    vertices: VecDeque<u32>,
    on_chain: Vec<bool>,
    closed: bool,
    fit: Option<AnyFit>,
    gate: Option<Gate>,
}

impl Chain<'_, '_> {
    /// Start again from the edge `a b`, with no fit.
    fn restart(&mut self, a: u32, b: u32) {
        for &v in &self.vertices {
            self.on_chain[v as usize] = false;
        }
        self.vertices.clear();
        self.closed = false;
        self.fit = None;
        self.gate = None;
        self.push(a, End::Back);
        self.push(b, End::Back);
    }

    fn push(&mut self, v: u32, end: End) {
        self.on_chain[v as usize] = true;
        match end {
            End::Front => self.vertices.push_front(v),
            End::Back => self.vertices.push_back(v),
        }
    }

    /// Fit the chain's vertices, from the last fit where there is one, with
    /// each vertex's normal (the sum of its triangles' vector areas) where
    /// every vertex has one.
    fn refit(&mut self) -> OgeomResult<()> {
        let surface = self.surface;
        let points: Vec<Point> = self.vertices.iter().map(|&v| surface.position(v)).collect();
        let normals: Vec<Vector> = self
            .vertices
            .iter()
            .map(|&v| {
                surface
                    .fan(v)
                    .iter()
                    .fold(Vector::ZERO, |sum, &t| sum + surface.cross[t as usize])
            })
            .collect();
        let facing = normals.iter().all(|n| n.magnitude() > 0.0);
        let fit = fitted(
            self.kind,
            &points,
            facing.then_some(normals.as_slice()),
            self.options,
            self.fit.as_ref(),
            self.tol,
        )?;
        self.gate = Some(Gate::of(
            &fit,
            &points,
            surface.single,
            self.options,
            self.tol,
        ));
        self.fit = Some(fit);
        Ok(())
    }

    /// Whether the chain turns past the crease limit by more than the noise
    /// could turn it, one of `tries` chains tried.
    ///
    /// Across each edge `a b` of the chain, between the triangles `a b c`
    /// and `a b d`, the surface turns through an angle whose excess over
    /// the limit, summed along the chain, is `S`. Each point moved by
    /// independent noise of sigma `σ` along each axis moves `S`, to first
    /// order, by its gradient there, summed over the edges it bounds: `S`
    /// varies as `σ² W`, `W` the sum of those gradients' squares. With the
    /// chain's fit measuring `σ` over `ν` degrees of freedom, the turn is
    /// past the limit where `S` exceeds the Student t at `ν` for `k'` times
    /// `σ √W`; `k'` is the normal deviate whose tail is the cutoff's `k`
    /// sigmas' over `tries`, so of all the chains tried, one whose edges
    /// turn by the limit at most passes as rarely as a point of plain noise
    /// lies `k` sigmas out. A short chain's fit measures its noise poorly
    /// and needs a larger excess. A cutoff that counts no sigmas,
    /// [`Cutoff::Distance`] or [`Cutoff::All`], measures no noise, and
    /// every chain passes.
    fn beyond_noise(&self, tries: usize) -> bool {
        let surface = self.surface;
        if !matches!(self.options.cutoff, Cutoff::Sigma(_) | Cutoff::Median(_)) {
            return true;
        }
        let Some(gate) = self.gate else {
            return false;
        };
        let mut edges: Vec<(u32, u32)> = self
            .vertices
            .iter()
            .zip(self.vertices.iter().skip(1))
            .map(|(&a, &b)| (a, b))
            .collect();
        if self.closed
            && let (Some(&last), Some(&first)) = (self.vertices.back(), self.vertices.front())
        {
            edges.push((last, first));
        }
        let mut excess = 0.0;
        let mut gradient: BTreeMap<u32, Vector> = BTreeMap::new();
        for (a, b) in edges {
            let beside: Vec<u32> = surface
                .fan(a)
                .iter()
                .map(|&t| surface.corners(t))
                .filter(|c| c.contains(&b))
                .filter_map(|c| c.into_iter().find(|&v| v != a && v != b))
                .collect();
            let &[c, d] = beside.as_slice() else {
                return false;
            };
            let points = [a, b, c, d];
            let Some((turn, slopes)) = fold(points.map(|v| surface.position(v))) else {
                return false;
            };
            excess += turn - self.options.crease_limit;
            for (v, slope) in points.into_iter().zip(slopes) {
                *gradient.entry(v).or_insert(Vector::ZERO) += slope;
            }
        }
        let spread = gradient.values().map(|g| g.dot(*g)).sum::<f64>().sqrt();
        #[allow(clippy::cast_precision_loss, reason = "far below 2^52 tries")]
        let tail = tails::normal_tail(gate.sigmas) / tries.max(1) as f64;
        let deviate = tails::t_at_normal(tails::normal_at_tail(tail), gate.dof);
        excess > deviate * gate.noise * spread
    }

    /// Run on from one end by one vertex, and say whether it did: along
    /// the one crease edge leaving it, before a fit with a finite gate;
    /// along the one whose far vertex the gate admits, after. Reaching the
    /// other end closes the chain.
    fn step(&mut self, end: End) -> bool {
        let count = self.vertices.len();
        let (at, came, other) = match end {
            End::Front => (self.vertices[0], self.vertices[1], self.vertices[count - 1]),
            End::Back => (
                self.vertices[count - 1],
                self.vertices[count - 2],
                self.vertices[0],
            ),
        };
        let open = |w: u32| w != came && (!self.on_chain[w as usize] || (w == other && count > 2));
        let gate = self.gate.filter(|g| g.threshold.is_finite());
        let mut onward = self.creases.from(at).iter().copied().filter(|&w| open(w));
        let next = onward.next();
        let Some(mut next) = next else {
            return false;
        };
        match gate {
            None => {
                if onward.next().is_some() {
                    return false;
                }
            }
            Some(gate) => {
                let mut admitted = core::iter::once(next)
                    .chain(onward)
                    .filter(|&w| gate.near(self.surface.position(w)));
                let Some(only) = admitted.next() else {
                    return false;
                };
                if admitted.next().is_some() {
                    return false;
                }
                next = only;
            }
        }
        if next == other {
            self.closed = true;
        } else {
            self.push(next, end);
        }
        true
    }

    /// Trace both ways until a whole pass against the latest fit adds no
    /// vertex or the chain closes, refitting whenever the chain has
    /// doubled since the last fit and once more whenever both ends stop;
    /// with no fit where the chain ends before it holds the vertices the
    /// fit needs.
    fn trace(&mut self) -> OgeomResult<()> {
        let needed = primitive::needed(self.kind.kind(), &self.options.fixed);
        let mut fitted_at = 0;
        loop {
            let mut stopped = [false; 2];
            while !self.closed && !(stopped[0] && stopped[1]) {
                for (k, end) in [End::Back, End::Front].into_iter().enumerate() {
                    if stopped[k] || self.closed {
                        continue;
                    }
                    if !self.step(end) {
                        stopped[k] = true;
                    }
                }
                let count = self.vertices.len();
                let due = if self.fit.is_none() {
                    count >= needed
                } else {
                    count >= 2 * fitted_at
                };
                if due && count > fitted_at {
                    self.refit()?;
                    fitted_at = count;
                    stopped = [false; 2];
                }
            }
            let count = self.vertices.len();
            if count == fitted_at {
                break;
            }
            if count < needed {
                break;
            }
            self.refit()?;
            fitted_at = count;
            if self.closed {
                break;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, reason = "test code")]
mod tests {
    use super::*;

    /// The fold's gradient matches its central differences, either way
    /// round the edge, across a flat, a convex and a concave edge, and a
    /// near fold back.
    #[test]
    fn the_fold_gradient_is_its_finite_difference() {
        let a = Point::new(0.1, -0.2, 0.05);
        let b = Point::new(0.9, 0.15, -0.1);
        for (c, d) in [
            (Point::new(0.4, 0.6, 0.0), Point::new(0.5, -0.7, 0.1)),
            (Point::new(0.3, 0.5, -0.3), Point::new(0.6, -0.4, -0.5)),
            (Point::new(0.3, 0.5, 0.4), Point::new(0.6, -0.4, 0.7)),
            (Point::new(0.2, 0.7, 0.1), Point::new(0.7, 0.6, 0.3)),
        ] {
            for points in [[a, b, c, d], [b, a, c, d], [a, b, d, c]] {
                let (_, slopes) = fold(points).unwrap();
                let h = 1e-6;
                for (i, slope) in slopes.iter().enumerate() {
                    for (k, unit) in [Vector::X, Vector::Y, Vector::Z].into_iter().enumerate() {
                        let moved = |by: f64| {
                            let mut p = points;
                            p[i] += unit * by;
                            fold(p).unwrap().0
                        };
                        let difference = (moved(h) - moved(-h)) / (2.0 * h);
                        let exact = slope.dot(unit);
                        assert!(
                            (difference - exact).abs() < 1e-6,
                            "point {i}, axis {k}: {difference} against {exact}"
                        );
                    }
                }
            }
        }
    }
}
