//! Curved regions grown by recognition: the mesh cut into flat patches and
//! the noise of its flat faces measured, regions grown seed by seed on one
//! verified canonical surface, the smooth regions left over rebuilt as
//! sweeps or B-spline patches, regions on one surface merged and
//! disconnected ones split, and spheres turned so their boundary circles
//! are latitudes.

use ogeom_core::{FastMap, OgeomResult, Tolerances};
use ogeom_math::{Direction, Frame, Point, Sphere, Vector};

use super::frames::border_loops;
use super::planner::{ENCLOSED_CLUSTER, REACH, is_sliver, sags_as_the_surface};
use super::segment::{
    SMALL_STAGE_TURN, angular_spread, chart, coplanar_groups, gradient, leans_as_the_surface,
    periodic, turn_of, unit_normal,
};
use super::snap::plane_through;
use super::weld::{Adjacency, Half, from_to};
use super::{Carrier, Curved, Groups, MeshSolidOptions};
use crate::recognize::{Canonical, recognize_curved, worst_deviation};

/// How much larger than a curved region's own facets a flat patch must be
/// to be a face of its own rather than facets of the curve.
const PATCH_SCALE: f64 = 20.0;

/// How much larger than a seed's flat patch another must be for a first
/// sample to take its facets without growing on from them.
const LEAF_SCALE: f64 = 16.0;

/// The mesh cut into flat patches: each grown from its largest facet across
/// edges to neighbours whose normals stay within a couple of degrees of
/// that facet's and whose corners stay on its plane. A curved region takes
/// a large patch whole or not at all.
struct FlatPatches {
    /// The patch each triangle is in.
    of: Vec<usize>,
    /// Each patch's triangles.
    members: Vec<Vec<usize>>,
    /// Each patch's area.
    area: Vec<f64>,
    /// Each triangle's area.
    facet: Vec<f64>,
}

impl FlatPatches {
    fn of(
        points: &[Point],
        triangles: &[[u32; 3]],
        adjacency: &Adjacency,
        normals: &[Vector],
        flat: f64,
    ) -> Self {
        // Two degrees: the folds a curve's facets make are larger except on
        // a very finely drawn one, whose facets' corners then leave the seed's
        // plane within a few rows.
        let cos_fold = 2.0_f64.to_radians().cos();
        let n = triangles.len();
        let facet: Vec<f64> = triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|v| points[v as usize]);
                (b - a).cross(c - a).magnitude() / 2.0
            })
            .collect();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&x, &y| facet[y].total_cmp(&facet[x]).then(x.cmp(&y)));
        let mut of = vec![usize::MAX; n];
        let (mut members, mut area) = (Vec::new(), Vec::new());
        for seed in order {
            if of[seed] != usize::MAX {
                continue;
            }
            let id = members.len();
            let normal = normals[seed];
            let origin = points[triangles[seed][0] as usize];
            let mut patch = vec![seed];
            of[seed] = id;
            let mut i = 0;
            while i < patch.len() {
                let t = patch[i];
                i += 1;
                for h in 3 * t..3 * t + 3 {
                    let Some(g) = adjacency.twin[h] else {
                        continue;
                    };
                    let u = g / 3;
                    if of[u] != usize::MAX || normals[u].dot(normal) < cos_fold {
                        continue;
                    }
                    if triangles[u]
                        .iter()
                        .all(|&v| (points[v as usize] - origin).dot(normal).abs() <= flat)
                    {
                        of[u] = id;
                        patch.push(u);
                    }
                }
            }
            area.push(patch.iter().map(|&t| facet[t]).sum());
            members.push(patch);
        }
        Self {
            of,
            members,
            area,
            facet,
        }
    }
}

/// How far the mesh's vertices stand off the flat faces they lie on: across
/// every edge between two triangles as good as coplanar (a turn of under a
/// twentieth of a degree), the height of the one's far corner over the
/// other's plane, at the ninetieth percentile. Zero where the mesh has no
/// such edges. A curve drawn finely turns that little too, steadily from
/// one edge to the next, where scatter turns at random; each height is
/// scaled by [`STEADY_GAIN`] times how far its edge's turn strays from its
/// neighbours' ([`steady_miss`]), up to the whole height, and an edge
/// whose neighbours stray from it by more than [`STRAY_LIMIT`] times its
/// turn is a lull in a surface that bends, and counts for nothing.
pub(super) fn flat_noise(points: &[Point], triangles: &[[u32; 3]], adjacency: &Adjacency) -> f64 {
    let cos_level = 0.05_f64.to_radians().cos();
    let normals: Vec<Vector> = triangles.iter().map(|&t| unit_normal(points, t)).collect();
    let fans = Fans::new(points.len(), triangles);
    let mut edges: Vec<(f64, Half, Half)> = Vec::new();
    for (h, twin) in adjacency.twin.iter().enumerate() {
        let Some(g) = *twin else {
            continue;
        };
        if g < h {
            continue;
        }
        let (t, u) = (h / 3, g / 3);
        let (a, b) = (normals[t], normals[u]);
        if a.dot(b) < cos_level {
            continue;
        }
        let far = triangles[u][(g % 3 + 2) % 3];
        let base = points[triangles[t][0] as usize];
        edges.push(((points[far as usize] - base).dot(a).abs(), h, g));
    }
    if edges.is_empty() {
        return 0.0;
    }
    // The percentile is the `kept`-th largest scaled height. Scaling only
    // lowers a height, so the edges are scaled tallest first, and once
    // `kept` scaled heights reach the next edge's unscaled one, no edge
    // left can change the answer.
    let kept = edges.len() - (edges.len() * 9 / 10).min(edges.len() - 1);
    edges.sort_unstable_by(|x, y| y.0.total_cmp(&x.0));
    // The scaled heights are not negative, so their bits order as they do.
    let mut tallest: std::collections::BinaryHeap<std::cmp::Reverse<u64>> =
        std::collections::BinaryHeap::with_capacity(kept + 1);
    for &(height, h, g) in &edges {
        if tallest.len() == kept
            && tallest
                .peek()
                .is_some_and(|least| f64::from_bits(least.0) >= height)
        {
            break;
        }
        let value = if height > 0.0 {
            match steady_miss(points, triangles, &normals, adjacency, &fans, h, g) {
                Some(miss) if miss > STRAY_LIMIT => 0.0,
                Some(miss) => height * (STEADY_GAIN * miss).min(1.0),
                None => height,
            }
        } else {
            0.0
        };
        tallest.push(std::cmp::Reverse(value.to_bits()));
        if tallest.len() > kept {
            tallest.pop();
        }
    }
    tallest.peek().map_or(0.0, |least| f64::from_bits(least.0))
}

/// How far the turn across a mesh edge strays from the turns of the edges
/// around it, against the turn itself. Two readings are taken, and the
/// nearer kept. Across: the edges beside it, about as long, within
/// [`PARALLEL_TURN`] of its direction, overlapping it along its length and
/// sharing no vertex with it; between the nearest on either side the turn
/// should lie, and with one side only, on the line through the nearest two
/// there, or where those two turn the same way, between the nearest's
/// turn and twice or half it as they grow or shrink toward the edge; with
/// one edge only, at its turn. Along: the edges continuing it past
/// either end, the straightest within [`ALONG_TURN`] of its direction;
/// between their turns, or at the one's. Every edge compared lies between
/// triangles within [`NEAR_TURN`] of the edge's own. A surface drawn finely
/// turns steadily from one edge to the next and strays little; scatter
/// turns either way at random and strays by as much as it turns. `None`
/// where no edge is found to compare with.
fn steady_miss(
    points: &[Point],
    triangles: &[[u32; 3]],
    normals: &[Vector],
    adjacency: &Adjacency,
    fans: &Fans,
    h: Half,
    g: Half,
) -> Option<f64> {
    let cos_near = NEAR_TURN.to_radians().cos();
    let cos_parallel = PARALLEL_TURN.to_radians().cos();
    let cos_along = ALONG_TURN.to_radians().cos();
    let (t, u) = (h / 3, g / 3);
    let normal = normals[t];
    let turn = signed_turn(points, triangles, h, g);
    if turn == 0.0 {
        return None;
    }
    let (i, j) = from_to(triangles, h);
    let (p, q) = (points[i as usize], points[j as usize]);
    let length = p.distance(q);
    let along = (q - p) / length;
    let across = normal.cross(along);
    let middle = p + (q - p) * 0.5;
    let near = |s: usize| normals[s].dot(normal) >= cos_near;
    // The triangles within three steps of the edge's two.
    let mut ring = vec![t, u];
    let mut start = 0;
    for _ in 0..3 {
        let end = ring.len();
        for r in start..end {
            for k in 3 * ring[r]..3 * ring[r] + 3 {
                if let Some(o) = adjacency.twin[k]
                    && near(o / 3)
                    && !ring.contains(&(o / 3))
                {
                    ring.push(o / 3);
                }
            }
        }
        start = end;
    }
    let mut below: Vec<(f64, f64)> = Vec::new();
    let mut above: Vec<(f64, f64)> = Vec::new();
    for &s in &ring {
        for k in 3 * s..3 * s + 3 {
            let Some(o) = adjacency.twin[k] else {
                continue;
            };
            if (o < k && ring.contains(&(o / 3))) || !near(o / 3) {
                continue;
            }
            let (c, d) = from_to(triangles, k);
            if c == i || c == j || d == i || d == j {
                continue;
            }
            let (c, d) = (points[c as usize], points[d as usize]);
            let other = c.distance(d);
            if other < 0.5 * length
                || other > 2.0 * length
                || ((d - c) / other).dot(along).abs() < cos_parallel
            {
                continue;
            }
            let offset = (c + (d - c) * 0.5) - middle;
            let off = offset.dot(across);
            if off.abs() <= 1e-3 * length
                || offset.dot(along).abs() > 0.5 * (length + other)
                || (d - c).dot(across).abs() > 0.25 * off.abs()
            {
                continue;
            }
            let side = if off < 0.0 { &mut below } else { &mut above };
            side.push((off.abs(), signed_turn(points, triangles, k, o)));
        }
    }
    for side in [&mut below, &mut above] {
        side.sort_by(|x, y| x.0.total_cmp(&y.0));
    }
    let mut ends: [Option<(f64, f64)>; 2] = [None, None];
    for (end, at, out) in [(0, i, -along), (1, j, along)] {
        for &s in fans.around(at) {
            for k in 3 * s as usize..3 * s as usize + 3 {
                let Some(o) = adjacency.twin[k] else {
                    continue;
                };
                let (c, d) = from_to(triangles, k);
                if o < k || (c != at && d != at) || !(near(s as usize) && near(o / 3)) {
                    continue;
                }
                let far = if c == at { d } else { c };
                if far == i || far == j {
                    continue;
                }
                let step = points[far as usize] - points[at as usize];
                let straight = step.dot(out) / step.magnitude();
                if straight >= cos_along && ends[end].is_none_or(|(best, _)| straight > best) {
                    ends[end] = Some((straight, signed_turn(points, triangles, k, o)));
                }
            }
        }
    }
    // Each reading is the span of the turns either side and the turn the
    // line between them gives at the edge's place.
    let beside = match (below.first(), above.first()) {
        (Some(&(o0, t0)), Some(&(o1, t1))) => {
            Some((t0.min(t1), t0.max(t1), (t0 * o1 + t1 * o0) / (o0 + o1)))
        }
        _ => {
            let side = if below.is_empty() { &above } else { &below };
            side.first().map(|&(o0, t0)| {
                match side.iter().find(|&&(o, _)| o > 1.5 * o0) {
                    // Two turns the same way that grow toward the edge, as
                    // where a bend tightens toward a surface's last row,
                    // put the edge's turn between the nearest and twice
                    // it; two that shrink, between the nearest and half
                    // of it.
                    Some(&(_, t1)) if t0 * t1 > 0.0 && t0.abs() != t1.abs() => {
                        let reach = if t0.abs() > t1.abs() { 2.0 } else { 0.5 };
                        let line = t0 * (1.0 + reach) / 2.0;
                        (t0.min(reach * t0), t0.max(reach * t0), line)
                    }
                    Some(&(o1, t1)) => {
                        let line = (t0 * o1 - t1 * o0) / (o1 - o0);
                        (line, line, line)
                    }
                    None => (t0, t0, t0),
                }
            })
        }
    };
    let onward = match ends {
        [Some((_, t0)), Some((_, t1))] => Some((t0.min(t1), t0.max(t1), 0.5 * (t0 + t1))),
        [Some((_, t0)), None] | [None, Some((_, t0))] => Some((t0, t0, t0)),
        [None, None] => None,
    };
    // Turning the same way as its neighbours, the edge need only lie
    // between them; where the turns change sign, a curve passing through
    // flat changes linearly, and the edge must lie on the line.
    [beside, onward]
        .into_iter()
        .flatten()
        .map(|(low, high, line)| {
            let miss = if low * turn > 0.0 && high * turn > 0.0 {
                (low - turn).max(turn - high).max(0.0)
            } else {
                (turn - line).abs()
            };
            miss / turn.abs()
        })
        .min_by(f64::total_cmp)
}

/// The triangles round each vertex.
struct Fans {
    start: Vec<usize>,
    triangles: Vec<u32>,
}

impl Fans {
    fn new(vertices: usize, triangles: &[[u32; 3]]) -> Self {
        let mut start = vec![0; vertices + 1];
        for t in triangles {
            for &v in t {
                start[v as usize + 1] += 1;
            }
        }
        for v in 0..vertices {
            start[v + 1] += start[v];
        }
        let mut fill = start.clone();
        let mut list = vec![0; start[vertices]];
        for (t, tri) in triangles.iter().enumerate() {
            for &v in tri {
                list[fill[v as usize]] = u32::try_from(t).unwrap_or(u32::MAX);
                fill[v as usize] += 1;
            }
        }
        Self {
            start,
            triangles: list,
        }
    }

    fn around(&self, v: u32) -> &[u32] {
        &self.triangles[self.start[v as usize]..self.start[v as usize + 1]]
    }
}

/// How many times its stray from its neighbours' turns an edge's height
/// counts as scatter: scatter strays by about as much as it turns, and a
/// height is kept whole unless the stray is under half the turn.
const STEADY_GAIN: f64 = 2.0;

/// How many times its own turn the edges an edge is compared with may
/// stray from it before the edge is taken for a lull in a surface that
/// bends, not part of a flat face: on a flat face they turn by about as
/// much as it does, either way, while in a bend they turn far more.
const STRAY_LIMIT: f64 = 3.0;

/// The widest turn, in degrees, between the triangles either side of an
/// edge whose scatter is measured and the triangles beside it whose edges
/// it is compared with.
const NEAR_TURN: f64 = 5.0;

/// The widest angle, in degrees, between an edge and one it is compared
/// with as running beside it.
const PARALLEL_TURN: f64 = 10.0;

/// The widest angle, in degrees, between an edge and one it is compared
/// with as continuing it past an end: a quarter circle drawn in three
/// chords turns by thirty.
const ALONG_TURN: f64 = 35.0;

/// The turn across a mesh edge from the triangle of half-edge `h` to the
/// triangle of its twin `g`, signed about the half-edge's direction, so
/// that on a consistently wound mesh every fold one way has one sign.
fn signed_turn(points: &[Point], triangles: &[[u32; 3]], h: Half, g: Half) -> f64 {
    let a = unit_normal(points, triangles[h / 3]);
    let b = unit_normal(points, triangles[g / 3]);
    let (p, q) = from_to(triangles, h);
    let along = points[q as usize] - points[p as usize];
    a.cross(b).dot(along / along.magnitude()).atan2(a.dot(b))
}

/// Keep only the largest edge-connected piece of a set of triangles.
fn keep_largest_piece(region: &mut Vec<usize>, adjacency: &Adjacency) {
    let members: ogeom_core::FastSet<usize> = region.iter().copied().collect();
    let mut piece: FastMap<usize, usize> =
        FastMap::with_capacity_and_hasher(region.len(), Default::default());
    let mut sizes: Vec<usize> = Vec::new();
    for &start in region.iter() {
        if piece.contains_key(&start) {
            continue;
        }
        let id = sizes.len();
        piece.insert(start, id);
        let mut stack = vec![start];
        let mut size = 0;
        while let Some(t) = stack.pop() {
            size += 1;
            for h in 3 * t..3 * t + 3 {
                if let Some(twin) = adjacency.twin[h] {
                    let other = twin / 3;
                    if members.contains(&other) && !piece.contains_key(&other) {
                        piece.insert(other, id);
                        stack.push(other);
                    }
                }
            }
        }
        sizes.push(size);
    }
    if sizes.len() < 2 {
        return;
    }
    let largest = (0..sizes.len()).max_by_key(|&i| sizes[i]).unwrap_or(0);
    region.retain(|t| piece.get(t) == Some(&largest));
}

/// Recognized regions sharing a mesh edge and lying on one surface are one
/// region: of one kind, each within twice the distance of the other's (a
/// band peeled of a flat face's facets can leave its two ends to be fitted
/// apart); of any canonical kind, the smaller within the distance of the
/// larger's (a few facets of a torus can lie on a sphere as exactly as on
/// the torus). The merged region keeps the larger's surface.
pub(super) fn merge_same_surface(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    flat: f64,
) {
    loop {
        let mut pair: Option<(usize, usize)> = None;
        'find: for (t, tri) in triangles.iter().enumerate() {
            let g = groups.of[t];
            let Some(Carrier::Curved(a)) = groups.carriers.get(g) else {
                continue;
            };
            for h in 3 * t..3 * t + 3 {
                let Some(twin) = adjacency.twin[h] else {
                    continue;
                };
                let o = groups.of[twin / 3];
                if o == g {
                    continue;
                }
                let Some(Carrier::Curved(b)) = groups.carriers.get(o) else {
                    continue;
                };
                let canonical = |c: &Curved| {
                    c.patch.is_none()
                        && !matches!(c.shape, Canonical::Swept(_) | Canonical::Plane(_))
                };
                let same = core::mem::discriminant(&a.shape) == core::mem::discriminant(&b.shape);
                if !same && !(canonical(a) && canonical(b)) {
                    continue;
                }
                // Each was fitted to its own vertices within the distance,
                // so the other's lie within twice it where the two are one.
                let on = |shape: &Canonical, vertices: &[u32], within: f64| {
                    vertices
                        .iter()
                        .all(|&v| shape.distance_to(points[v as usize]) <= within)
                };
                if same
                    && on(&a.shape, &b.vertices, flat * 2.0)
                    && on(&b.shape, &a.vertices, flat * 2.0)
                {
                    pair = Some((g.min(o), g.max(o)));
                    break 'find;
                }
                // A small piece fitted on its own leans with its few
                // vertices' slop, and the larger's fit may miss it by more;
                // where its vertices lie on the larger's surface, it is part
                // of the larger.
                let (large, small) = if a.vertices.len() >= b.vertices.len() {
                    ((g, a), (o, b))
                } else {
                    ((o, b), (g, a))
                };
                if on(&large.1.shape, &small.1.vertices, flat) {
                    pair = Some((large.0, small.0));
                    break 'find;
                }
            }
            let _ = tri;
        }
        let Some((keep, gone)) = pair else {
            return;
        };
        let Carrier::Curved(absorbed) =
            core::mem::replace(&mut groups.carriers[gone], Carrier::Gone)
        else {
            return;
        };
        for of in &mut groups.of {
            if *of == gone {
                *of = keep;
            }
        }
        if let Carrier::Curved(kept) = &mut groups.carriers[keep] {
            let pts: Vec<Point> = absorbed
                .vertices
                .iter()
                .map(|&v| points[v as usize])
                .collect();
            kept.vertices.extend(absorbed.vertices);
            kept.vertices.sort_unstable();
            kept.vertices.dedup();
            kept.deviation = kept
                .deviation
                .max(absorbed.deviation)
                .max(worst_deviation(&kept.shape, &pts));
            kept.fitted = kept.fitted.max(absorbed.fitted);
        }
    }
}

/// One recognized region per connected patch. Triangles dropped from a
/// region for touching a vertex its fit refused can take with them the only
/// triangles joining the rest, and a face is built from one region's
/// boundary: the patches past the first would be lost from the shell. Each
/// patch past the first becomes a region of its own on the same surface.
pub(super) fn split_disconnected(
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
) {
    let count = groups.carriers.len();
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); count];
    for (t, &g) in groups.of.iter().enumerate() {
        if g < count && matches!(groups.carriers[g], Carrier::Curved(_)) {
            members[g].push(t);
        }
    }
    for (g, tris) in members.into_iter().enumerate() {
        if tris.len() < 2 {
            continue;
        }
        let mut patch = vec![usize::MAX; triangles.len()];
        let mut patches = 0_usize;
        for &start in &tris {
            if patch[start] != usize::MAX {
                continue;
            }
            patch[start] = patches;
            let mut stack = vec![start];
            while let Some(t) = stack.pop() {
                for h in 3 * t..3 * t + 3 {
                    if let Some(twin) = adjacency.twin[h] {
                        let other = twin / 3;
                        if groups.of[other] == g && patch[other] == usize::MAX {
                            patch[other] = patches;
                            stack.push(other);
                        }
                    }
                }
            }
            patches += 1;
        }
        if patches < 2 {
            continue;
        }
        let Carrier::Curved(template) = &groups.carriers[g] else {
            continue;
        };
        let template = template.clone();
        let base = groups.carriers.len();
        for k in 1..patches {
            let mut vertices: Vec<u32> = tris
                .iter()
                .filter(|&&t| patch[t] == k)
                .flat_map(|&t| triangles[t])
                .collect();
            vertices.sort_unstable();
            vertices.dedup();
            groups.carriers.push(Carrier::Curved(Curved {
                vertices,
                ..template.clone()
            }));
        }
        let mut first: Vec<u32> = tris
            .iter()
            .filter(|&&t| patch[t] == 0)
            .flat_map(|&t| triangles[t])
            .collect();
        first.sort_unstable();
        first.dedup();
        if let Carrier::Curved(curved) = &mut groups.carriers[g] {
            curved.vertices = first;
        }
        for &t in &tris {
            if patch[t] > 0 {
                groups.of[t] = base + patch[t] - 1;
            }
        }
    }
}

/// What a seed's first samples came to: the triangles and vertices
/// gathered, those fitted, how many triangles the smallest sample held, the
/// fit if one held, with which of the fitted vertices it kept.
struct FirstFit {
    region: Vec<usize>,
    vertices: Vec<u32>,
    shared: Vec<u32>,
    first_sample: usize,
    /// Whether the first sample already turned through a tight bend: on a
    /// coarse mesh it can span several surfaces, and a failed fit says
    /// nothing about the triangles in it but the seed.
    wide: bool,
    found: Option<(crate::recognize::Recognized, Vec<bool>)>,
}

/// The mesh as recognition reads it.
struct Surfaces<'a> {
    points: &'a [Point],
    triangles: &'a [[u32; 3]],
    adjacency: &'a Adjacency,
    normals: Vec<Vector>,
    patches: &'a FlatPatches,
    cos_crease: f64,
    cos_flat: f64,
    flat: f64,
    tol: Tolerances,
}

impl Surfaces<'_> {
    fn turn(&self, h: Half) -> Option<f64> {
        self.adjacency.twin[h].map(|g| self.normals[h / 3].dot(self.normals[g / 3]))
    }

    fn curved(&self, h: Half) -> bool {
        self.turn(h)
            .is_some_and(|c| c >= self.cos_crease && c < self.cos_flat)
    }

    fn smooth(&self, h: Half) -> bool {
        self.turn(h).is_some_and(|c| c >= self.cos_crease)
    }

    fn bends(&self, t: usize) -> bool {
        (3 * t..3 * t + 3).any(|h| self.curved(h))
    }

    /// Sample points with normals averaged over the region's triangles.
    fn samples(&self, vertices: &[u32], region: &[usize]) -> (Vec<Point>, Vec<Vector>) {
        let mut sum: FastMap<u32, Vector> =
            FastMap::with_capacity_and_hasher(vertices.len(), Default::default());
        for &t in region {
            for &v in &self.triangles[t] {
                *sum.entry(v).or_insert(Vector::ZERO) += self.normals[t];
            }
        }
        let pts = vertices.iter().map(|&v| self.points[v as usize]).collect();
        let nrm = vertices
            .iter()
            .map(|v| {
                let s = sum.get(v).copied().unwrap_or(Vector::Z);
                let m = s.magnitude();
                if m > 0.0 { s / m } else { Vector::Z }
            })
            .collect();
        (pts, nrm)
    }

    /// The region's edges as segments, a few hundred at most, evenly
    /// through the region.
    fn chords(&self, region: &[usize]) -> Vec<(Point, Point)> {
        let stride = region.len().div_ceil(100).max(1);
        region
            .iter()
            .step_by(stride)
            .flat_map(|&t| {
                let [a, b, c] = self.triangles[t].map(|v| self.points[v as usize]);
                [(a, b), (b, c), (c, a)]
            })
            .collect()
    }

    /// A seed's first samples and their fit, read from where the faces and
    /// the retired seeds stand; it changes nothing.
    ///
    /// Samples grow across smooth edges into triangles where the surface
    /// is seen to bend, each with an edge across which it turns. Within a
    /// strip of a cylinder or a torus the two triangles of a cell meet
    /// flat, and each still bends across its other side; a flat face met
    /// tangentially (a fillet's run-out) lends only its triangles along
    /// the tangent line, whose far corners the trimmed fit drops. The fit
    /// is tried as the sample grows: a narrow fillet fits from a few dozen
    /// vertices and would take in its neighbours' by a hundred; a patch of a
    /// thick torus needs the hundred to show its tube.
    fn first_fit(&self, seed: usize, of: &[usize], tried: &[bool]) -> FirstFit {
        const STAGES: [usize; 4] = [12, 24, 60, 150];
        let mut held: ogeom_core::FastSet<usize> = ogeom_core::FastSet::from_iter([seed]);
        let mut seen: ogeom_core::FastSet<u32> = ogeom_core::FastSet::default();
        let mut region = vec![seed];
        let mut vertices: Vec<u32> = Vec::new();
        let take = |t: usize, vertices: &mut Vec<u32>, seen: &mut ogeom_core::FastSet<u32>| {
            for &v in &self.triangles[t] {
                if seen.insert(v) {
                    vertices.push(v);
                }
            }
        };
        take(seed, &mut vertices, &mut seen);
        let mut queue: std::collections::VecDeque<usize> = std::collections::VecDeque::from([seed]);
        let mut found = None;
        let mut shared = Vec::new();
        let mut first_sample = usize::MAX;
        let mut wide = false;
        for target in STAGES {
            while vertices.len() < target {
                let Some(next) = queue.pop_front() else {
                    break;
                };
                for h in 3 * next..3 * next + 3 {
                    let Some(g) = self.adjacency.twin[h] else {
                        continue;
                    };
                    let other = g / 3;
                    if !self.smooth(h) || held.contains(&other) {
                        continue;
                    }
                    if of[other] != usize::MAX
                        || tried[other]
                        || !(self.curved(h) || self.bends(other))
                    {
                        continue;
                    }
                    held.insert(other);
                    region.push(other);
                    take(other, &mut vertices, &mut seen);
                    // A facet of a flat patch much larger than the seed's
                    // own (a flat face beside a fillet, where the seed's
                    // patch is a row of it) lends its corners but leads
                    // nowhere: on a coarse mesh it borders other fillets
                    // than the seed's.
                    let (patch, own) = (self.patches.of[other], self.patches.of[seed]);
                    if self.patches.area[patch] <= self.patches.area[own] * LEAF_SCALE {
                        queue.push_back(other);
                    }
                }
            }
            // Fitted on the vertices two or more of the sampled triangles
            // share: a corner of a flat face across a tangent line is
            // touched by the one triangle that reached it, and a single
            // such corner far off the surface decides a least-squares axis.
            shared = {
                let mut count: FastMap<u32, u32> =
                    FastMap::with_capacity_and_hasher(vertices.len(), Default::default());
                for &t in &region {
                    for &v in &self.triangles[t] {
                        *count.entry(v).or_insert(0) += 1;
                    }
                }
                vertices
                    .iter()
                    .copied()
                    .filter(|v| count[v] >= 2)
                    .collect::<Vec<u32>>()
            };
            first_sample = first_sample.min(region.len());
            if target == STAGES[0] {
                let normals: Vec<Vector> = region.iter().map(|&t| self.normals[t]).collect();
                wide = turn_of(&normals) >= SMALL_STAGE_TURN;
            }
            if shared.len() < 8 {
                if queue.is_empty() {
                    break;
                }
                continue;
            }
            let (pts, nrm) = self.samples(&shared, &region);
            // The first, smallest stage is for a coarse mesh, where a dozen
            // vertices already span a tight bend; on a fine one they lie
            // nearly flat and say little about the surface they are on.
            if target == STAGES[0] && !wide {
                continue;
            }
            let chords = self.chords(&region);
            match crate::recognize::recognize_trimmed(&pts, &nrm, &chords, self.flat, self.tol) {
                Ok(fit) => {
                    found = Some(fit);
                    break;
                }
                // A larger sample is worth fitting only where this one came
                // near for its size (within a hundredth of its own span),
                // as a small patch of a thick torus does, its tube not yet
                // seen; a free-form patch misses by more, and is left.
                // A sample too small to be fitted as a torus has not been
                // asked whether it lies on one, and its miss says nothing of
                // what a larger sample fits: the next stage is taken.
                Err(closest)
                    if queue.is_empty()
                        || (closest > span(&pts) * 1e-2
                            && pts.len() >= crate::recognize::TORUS_SAMPLES) =>
                {
                    break;
                }
                Err(_) => {}
            }
        }
        FirstFit {
            region,
            vertices,
            shared,
            first_sample,
            wide,
            found,
        }
    }
}

/// Grow the recognized regions, seed by seed in triangle order.
///
/// A seed's first fit is the costly part, and it only reads: so seeds are
/// fitted a batch at a time in parallel against where things stand, then
/// taken in order, and a seed whose gathering took a triangle an earlier
/// seed of the same batch has since claimed or retired is fitted again.
/// The answer is the one seed-by-seed order gives, at any thread count.
#[allow(clippy::too_many_lines, reason = "one growth, read in one place")]
pub(super) fn recognized_regions(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    options: &MeshSolidOptions,
    flat: f64,
    groups: &mut Groups,
    tol: Tolerances,
) {
    let n = triangles.len();
    let normals: Vec<Vector> = triangles.iter().map(|t| unit_normal(points, *t)).collect();
    let patches = FlatPatches::of(points, triangles, adjacency, &normals, flat);
    let mesh = Surfaces {
        points,
        triangles,
        adjacency,
        normals,
        patches: &patches,
        cos_crease: options.crease.cos(),
        cos_flat: options.coplanar_angle.cos(),
        flat,
        tol,
    };
    // A facet of a coarse mesh leans from the surface at its centre by up
    // to half the turn between facets, which the crease bounds.
    let agree = options.crease.cos();
    let mut tried = vec![false; n];
    // The batch in which each triangle's standing last changed.
    let mut changed = vec![0_u32; n];
    let batch_size = ogeom_core::parallel::threads().max(1) * 4;
    let mut batch = 0_u32;
    let eligible =
        |t: usize, of: &[usize], tried: &[bool]| of[t] == usize::MAX && !tried[t] && mesh.bends(t);
    // Seeds are taken in a fixed stride through the triangles rather than
    // one after another: a mesh lists neighbours together, and a batch of
    // neighbouring seeds would mostly gather what the first of them took.
    // The stride is a prime not dividing the count, so every triangle comes
    // up once.
    let stride = [7919_usize, 7907, 7901]
        .into_iter()
        .find(|p| !n.is_multiple_of(*p))
        .unwrap_or(1);
    let mut order: Vec<usize> = (0..n).map(|i| (i * stride) % n).collect();
    // Seeds whose wide first sample fitted nothing, with that sample. Once
    // every seed has been taken, those whose sample has since lost
    // triangles to other regions are taken once more, the rest of their
    // samples free again: a fillet sampled together with its unclaimed
    // neighbours fits nothing, and its triangles, retired one by one, leave
    // ever smaller samples; on its own, its neighbours claimed, it fits.
    let mut straddled: Vec<(usize, Vec<usize>)> = Vec::new();
    let mut again = true;
    let mut next = 0;
    loop {
        if next == order.len() {
            let samples: Vec<(usize, Vec<usize>)> = std::mem::take(&mut straddled)
                .into_iter()
                .filter(|(seed, sample)| {
                    groups.of[*seed] == usize::MAX
                        && sample.iter().any(|&t| groups.of[t] != usize::MAX)
                })
                .collect();
            if !std::mem::take(&mut again) || samples.is_empty() {
                break;
            }
            for (seed, sample) in samples {
                for &t in &sample {
                    if groups.of[t] == usize::MAX {
                        tried[t] = false;
                    }
                }
                order.push(seed);
            }
        }
        batch += 1;
        let mut seeds = Vec::with_capacity(batch_size);
        while next < order.len() && seeds.len() < batch_size {
            if eligible(order[next], &groups.of, &tried) {
                seeds.push(order[next]);
            }
            next += 1;
        }
        let fits = {
            let (of, tried) = (&groups.of, &tried);
            ogeom_core::parallel::map_ordered(&seeds, |_, &seed| mesh.first_fit(seed, of, tried))
        };
        for (seed, fit) in seeds.into_iter().zip(fits) {
            if !eligible(seed, &groups.of, &tried) {
                continue;
            }
            // A triangle's standing only moves one way (free to taken), so
            // a triangle the gathering passed over stays passed over, and
            // only one it took can have changed what it gathers.
            let fit = if fit.region.iter().any(|&t| changed[t] == batch) {
                mesh.first_fit(seed, &groups.of, &tried)
            } else {
                fit
            };
            let FirstFit {
                mut region,
                mut vertices,
                shared,
                first_sample,
                wide,
                found,
                ..
            } = fit;
            let Some((found, keep)) = found else {
                // A fine sample that fits nothing lies on nothing canonical,
                // and its triangles are not seeded again; a wide one may
                // have straddled a fillet and its neighbours, and only the
                // seed is retired.
                let retired = if wide {
                    straddled.push((seed, region.clone()));
                    1
                } else {
                    first_sample
                };
                for &t in &region[..retired.min(region.len())] {
                    tried[t] = true;
                    changed[t] = batch;
                }
                continue;
            };
            // What the fit dropped, and what it never saw but misses the
            // fit: those vertices go, and the triangles that brought them.
            let kept: FastMap<u32, bool> = shared.iter().copied().zip(keep).collect();
            let dropped: ogeom_core::FastSet<u32> = vertices
                .iter()
                .copied()
                .filter(|v| {
                    !kept
                        .get(v)
                        .copied()
                        .unwrap_or_else(|| found.surface.distance_to(points[*v as usize]) <= flat)
                })
                .collect();
            if !dropped.is_empty() {
                let gathered = region.clone();
                region.retain(|&t| triangles[t].iter().all(|v| !dropped.contains(v)));
                if region.is_empty() {
                    for &t in &gathered[..first_sample.min(gathered.len())] {
                        tried[t] = true;
                        changed[t] = batch;
                    }
                    continue;
                }
                // What is left may have come apart where the dropped
                // triangles joined it. The largest piece is the claim; the
                // others stay free, to seed regions of their own or fall to
                // the planes.
                keep_largest_piece(&mut region, adjacency);
                let mut seen = ogeom_core::FastSet::default();
                vertices.clear();
                for &t in &region {
                    for &v in &triangles[t] {
                        if seen.insert(v) {
                            vertices.push(v);
                        }
                    }
                }
            }
            let mut shape = found.surface;
            // A facet of the region as it was first fitted, for telling a
            // flat face's facets from its own.
            let typical = {
                let mut areas: Vec<f64> = region.iter().map(|&t| patches.facet[t]).collect();
                areas.sort_by(f64::total_cmp);
                areas.get(areas.len() / 2).copied().unwrap_or(0.0)
            };
            let mut mine: ogeom_core::FastSet<usize> = region.iter().copied().collect();
            let mut seen: ogeom_core::FastSet<u32> = vertices.iter().copied().collect();

            // Then across any smooth edge, while the surface holds. A fit
            // from a small patch extrapolates only so far; when the growth
            // stalls with vertices gained since the last fit, the surface is
            // refitted to everything held and the rim tried again.
            let mut fitted_at = vertices.len();
            loop {
                let mut i = 0;
                while i < region.len() {
                    let t = region[i];
                    i += 1;
                    for h in 3 * t..3 * t + 3 {
                        let Some(g) = adjacency.twin[h] else {
                            continue;
                        };
                        let other = g / 3;
                        if mine.contains(&other) || groups.of[other] != usize::MAX {
                            continue;
                        }
                        let corners = triangles[other].map(|v| points[v as usize]);
                        if corners.iter().any(|p| shape.distance_to(*p) > flat) {
                            continue;
                        }
                        // Across a smooth edge, the facet's normal agrees
                        // with the surface's. Across a sharper one (a coarse
                        // mesh spanning two rows of a small fillet in one
                        // triangle), the surface must account for the whole
                        // lean.
                        if mesh.smooth(h) {
                            let centroid = Point::from_vector(
                                (corners[0].to_vector()
                                    + corners[1].to_vector()
                                    + corners[2].to_vector())
                                    / 3.0,
                            );
                            let direction = gradient(&shape, centroid);
                            let m = direction.magnitude();
                            if m == 0.0 || (direction.dot(mesh.normals[other]) / m).abs() < agree {
                                continue;
                            }
                        } else if !leans_as_the_surface(&shape, corners, mesh.normals[other])
                            || !sags_as_the_surface(&shape, corners, flat)
                        {
                            continue;
                        }
                        // A facet of a flat face much larger than the
                        // region's own (a plane tangent to it) comes with its
                        // whole face or not at all: its corners by the
                        // tangent line lie on the surface, its far ones not.
                        // Every facet of it must lie on the surface as one
                        // of the surface's own would: a disc capping a cone
                        // has every corner on the cone's rim, and the facet
                        // across its middle leans from the cone and sags
                        // far below it.
                        let patch = patches.of[other];
                        let taken: Vec<usize> = if patches.area[patch] > typical * PATCH_SCALE {
                            let whole = &patches.members[patch];
                            let on = whole.iter().all(|&t| {
                                let corners = triangles[t].map(|v| points[v as usize]);
                                mine.contains(&t)
                                    || (groups.of[t] == usize::MAX
                                        && corners.iter().all(|p| shape.distance_to(*p) <= flat)
                                        && sags_as_the_surface(&shape, corners, flat)
                                        && (is_sliver(corners)
                                            || leans_as_the_surface(
                                                &shape,
                                                corners,
                                                mesh.normals[t],
                                            )))
                            });
                            if !on {
                                continue;
                            }
                            whole
                                .iter()
                                .copied()
                                .filter(|t| !mine.contains(t))
                                .collect()
                        } else {
                            vec![other]
                        };
                        for other in taken {
                            mine.insert(other);
                            region.push(other);
                            for &v in &triangles[other] {
                                if seen.insert(v) {
                                    vertices.push(v);
                                }
                            }
                        }
                    }
                }
                if vertices.len() <= fitted_at {
                    break;
                }
                fitted_at = vertices.len();
                let (pts, nrm) = mesh.samples(&vertices, &region);
                match recognize_curved(&pts, &nrm, &mesh.chords(&region), flat, tol) {
                    Some(better) => shape = better.surface,
                    None => break,
                }
            }
            // A few triangles the region surrounds on every side, with their
            // corners on the surface, belong to it whatever their normals:
            // a sliver's plane through three nearly collinear points on the
            // surface tilts far from the surface's normal, as in the fans
            // round a sphere's pole, where several such slivers touch.
            let mut claimed: ogeom_core::FastSet<usize> = ogeom_core::FastSet::default();
            let rim: Vec<usize> = region
                .iter()
                .flat_map(|&t| (3 * t..3 * t + 3).filter_map(|h| adjacency.twin[h]))
                .map(|g| g / 3)
                .filter(|&other| !mine.contains(&other) && groups.of[other] == usize::MAX)
                .collect();
            for start in rim {
                if claimed.contains(&start) {
                    continue;
                }
                let mut cluster = vec![start];
                let mut inside: ogeom_core::FastSet<usize> =
                    ogeom_core::FastSet::from_iter([start]);
                let mut surrounded = true;
                let mut i = 0;
                while i < cluster.len() && surrounded {
                    let t = cluster[i];
                    i += 1;
                    for h in 3 * t..3 * t + 3 {
                        let Some(g) = adjacency.twin[h] else {
                            surrounded = false;
                            break;
                        };
                        let next = g / 3;
                        if mine.contains(&next) || inside.contains(&next) {
                            continue;
                        }
                        if groups.of[next] != usize::MAX || cluster.len() >= ENCLOSED_CLUSTER {
                            surrounded = false;
                            break;
                        }
                        inside.insert(next);
                        cluster.push(next);
                    }
                }
                let on_surface = cluster.iter().all(|&t| {
                    let corners = triangles[t].map(|v| points[v as usize]);
                    corners.iter().all(|p| shape.distance_to(*p) <= flat)
                        && sags_as_the_surface(&shape, corners, flat)
                        && (is_sliver(corners)
                            || leans_as_the_surface(&shape, corners, mesh.normals[t]))
                });
                if surrounded && on_surface {
                    for t in cluster {
                        claimed.insert(t);
                        if mine.insert(t) {
                            region.push(t);
                        }
                    }
                }
            }
            // The vertices on the surface do not put the triangles on it: a
            // long facet across a flat stretch has its corners where a
            // surface through both ends of the stretch passes, and its
            // middle far from it. Where that facet is a large flat face's (a
            // plane meeting the surface along a tangent circle, triangulated
            // across it), it is the face's and is peeled off, the largest
            // piece left being the region; any other such facet says the
            // surface is wrong.
            let peel: Vec<usize> = region
                .iter()
                .copied()
                .filter(|&t| {
                    patches.area[patches.of[t]] > typical * PATCH_SCALE
                        && !sags_as_the_surface(
                            &shape,
                            triangles[t].map(|v| points[v as usize]),
                            flat,
                        )
                })
                .collect();
            if !peel.is_empty() {
                region.retain(|t| !peel.contains(t));
                keep_largest_piece(&mut region, adjacency);
                let mut seen = ogeom_core::FastSet::default();
                vertices.clear();
                for &t in &region {
                    for &v in &triangles[t] {
                        if seen.insert(v) {
                            vertices.push(v);
                        }
                    }
                }
                // What is left is asked again what it is: a band a row or
                // two high lies on a whole family of surfaces, and the one
                // chosen with the flat facets in is not the one without.
                // Where nothing is found, the surface held so far stands if
                // it still holds what is left: the flat facets gone, the
                // normals at the rim they shared lean to one side, and can
                // mislead the fit that the surface already answers.
                let (pts, nrm) = mesh.samples(&vertices, &region);
                match recognize_curved(&pts, &nrm, &mesh.chords(&region), flat, tol) {
                    Some(better) => shape = better.surface,
                    None if worst_deviation(&shape, &pts) <= flat => {}
                    None => {
                        for &t in &region {
                            tried[t] = true;
                            changed[t] = batch;
                        }
                        continue;
                    }
                }
            }
            let spans = |t: usize| {
                !sags_as_the_surface(&shape, triangles[t].map(|v| points[v as usize]), flat)
            };
            let spans_off = region.iter().any(|&t| spans(t));
            let pts: Vec<Point> = vertices.iter().map(|&v| points[v as usize]).collect();
            let deviation = worst_deviation(&shape, &pts);
            let flat_too = crate::recognize::is_flat(&pts, flat, tol);
            // A sphere or a torus through a handful of noisy vertices is
            // as likely the noise's as the part's: twice its unknowns are
            // asked for.
            let few = match shape {
                Canonical::Sphere(_) => vertices.len() < 8,
                Canonical::Torus(_) => vertices.len() < 14,
                _ => false,
            };
            if deviation > flat || flat_too || spans_off || region.len() < 2 || few {
                for &t in &region {
                    tried[t] = true;
                    changed[t] = batch;
                }
                continue;
            }
            let g = groups.carriers.len();
            for &t in &region {
                groups.of[t] = g;
                changed[t] = batch;
            }
            groups.carriers.push(Carrier::Curved(Curved {
                shape,
                deviation,
                fitted: deviation,
                centre: (0.0, 0.0),
                wraps: false,
                wraps_v: false,
                fixed: false,
                vertices,
                patch: None,
            }));
        }
    }
}

/// The smooth regions recognition left free, rebuilt as extrusions or
/// surfaces of revolution where one fits.
///
/// A region is the free triangles joined across edges that turn but do not
/// crease. Its vertices, with the normals averaged over its triangles, are
/// fitted as an extrusion and then as a surface of revolution, each held
/// to the coplanar distance at every vertex; a region that goes all the way
/// round (a closed profile, a whole turn) is left to its facets, as is one
/// whose triangles do not sag as the surface does.
pub(super) fn swept_regions(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    options: &MeshSolidOptions,
    flat: f64,
    groups: &mut Groups,
    tol: Tolerances,
) {
    let normals: Vec<Vector> = triangles.iter().map(|t| unit_normal(points, *t)).collect();
    let cos_crease = options.crease.cos();
    let turn = |h: Half| adjacency.twin[h].map(|g| normals[h / 3].dot(normals[g / 3]));
    for region in smooth_regions(triangles, adjacency, &normals, options, groups) {
        // The recognized bands the region runs into smoothly, offered with
        // it: a stretch of a wavy profile straight enough to pass for a
        // cone is part of the one sweep, and the sweep takes it where the
        // whole fits.
        let mut bands: Vec<usize> = Vec::new();
        for &t in &region {
            for h in 3 * t..3 * t + 3 {
                let (Some(c), Some(g)) = (turn(h), adjacency.twin[h]) else {
                    continue;
                };
                let other = groups.of[g / 3];
                if c >= cos_crease
                    && other != usize::MAX
                    && matches!(groups.carriers.get(other), Some(Carrier::Curved(_)))
                    && !bands.contains(&other)
                {
                    bands.push(other);
                }
            }
        }
        if !bands.is_empty() {
            let mut widened = region.clone();
            widened.extend((0..triangles.len()).filter(|&t| bands.contains(&groups.of[t])));
            if let Some(claim) = swept_claim(points, triangles, &normals, &widened, flat, tol) {
                for &b in &bands {
                    groups.carriers[b] = Carrier::Gone;
                }
                claim_swept(groups, &widened, claim);
                continue;
            }
        }
        if let Some(claim) = swept_claim(points, triangles, &normals, &region, flat, tol) {
            claim_swept(groups, &region, claim);
        }
    }
}

/// The free triangles joined across edges that turn but do not crease,
/// each such region holding at least [`SWEPT_TRIANGLES`] and an edge that
/// turns: the smooth regions no canonical surface took.
fn smooth_regions(
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    normals: &[Vector],
    options: &MeshSolidOptions,
    groups: &Groups,
) -> Vec<Vec<usize>> {
    let (cos_crease, cos_flat) = (options.crease.cos(), options.coplanar_angle.cos());
    let turn = |h: Half| adjacency.twin[h].map(|g| normals[h / 3].dot(normals[g / 3]));
    let mut seen = vec![false; triangles.len()];
    let mut regions = Vec::new();
    for seed in 0..triangles.len() {
        if seen[seed] || groups.of[seed] != usize::MAX {
            continue;
        }
        let mut region = vec![seed];
        seen[seed] = true;
        let mut bends = false;
        let mut i = 0;
        while i < region.len() {
            let t = region[i];
            i += 1;
            for h in 3 * t..3 * t + 3 {
                let (Some(c), Some(g)) = (turn(h), adjacency.twin[h]) else {
                    continue;
                };
                if c < cos_crease {
                    continue;
                }
                let other = g / 3;
                // Only a turn between its own triangles: a flat face meeting
                // a recognized one across a smooth edge does not bend.
                if groups.of[other] != usize::MAX {
                    continue;
                }
                if c < cos_flat {
                    bends = true;
                }
                if !seen[other] {
                    seen[other] = true;
                    region.push(other);
                }
            }
        }
        if bends && region.len() >= SWEPT_TRIANGLES {
            regions.push(region);
        }
    }
    regions
}

/// A sweep fitted to a region's vertices, with the normals averaged over
/// its triangles: held at every vertex, its triangles sagging as it does,
/// and its chart centre the vertices' mean. `None` where none fits or the
/// region goes all the way round.
pub(super) fn swept_claim(
    points: &[Point],
    triangles: &[[u32; 3]],
    normals: &[Vector],
    region: &[usize],
    flat: f64,
    tol: Tolerances,
) -> Option<Curved> {
    let mut sums: FastMap<u32, Vector> = FastMap::default();
    for &t in region {
        for &v in &triangles[t] {
            *sums.entry(v).or_insert(Vector::ZERO) += normals[t];
        }
    }
    let mut vertices: Vec<u32> = sums.keys().copied().collect();
    vertices.sort_unstable();
    let pts: Vec<Point> = vertices.iter().map(|&v| points[v as usize]).collect();
    let nrm: Vec<Vector> = vertices
        .iter()
        .map(|v| {
            let s = sums[v];
            let m = s.magnitude();
            if m > 0.0 { s / m } else { Vector::Z }
        })
        .collect();
    let found = crate::recognize_swept::fit_extrusion(&pts, &nrm, flat, tol)
        .or_else(|| crate::recognize_swept::fit_revolution(&pts, &nrm, flat, tol))?;
    let shape = Canonical::Swept(Box::new(crate::recognize::SweptShape::new(
        found.surface,
        tol,
    )));
    if !region
        .iter()
        .all(|&t| sags_as_the_surface(&shape, triangles[t].map(|v| points[v as usize]), flat))
    {
        return None;
    }
    let charts: Vec<(f64, f64)> = pts.iter().filter_map(|p| chart(&shape, *p, tol)).collect();
    if charts.len() != pts.len() {
        return None;
    }
    let (pu, pv) = periodic(&shape);
    // Each chart direction's centre, and whether the region goes all the
    // way round it: then its centre is half a turn from the seam, as a
    // canonical band's is.
    let centre_of = |mut values: Vec<f64>, periodic: bool| -> (f64, bool) {
        if periodic {
            let (mean, gap) = angular_spread(&mut values);
            if gap < core::f64::consts::FRAC_PI_2 {
                (core::f64::consts::PI, true)
            } else {
                (mean, false)
            }
        } else {
            #[allow(clippy::cast_precision_loss, reason = "vertex counts are small")]
            (
                values.iter().sum::<f64>() / values.len().max(1) as f64,
                false,
            )
        }
    };
    let (cu, wraps) = centre_of(charts.iter().map(|c| c.0).collect(), pu);
    let (cv, wraps_v) = centre_of(charts.iter().map(|c| c.1).collect(), pv);
    Some(Curved {
        shape,
        deviation: found.deviation,
        fitted: found.deviation,
        centre: (cu, cv),
        wraps,
        wraps_v,
        fixed: true,
        vertices,
        patch: None,
    })
}

/// The smooth regions recognition and the sweeps left free, each rebuilt
/// on one fitted B-spline patch where the region is one disk and the patch
/// verifies.
///
/// A region is the free triangles joined across edges that turn but do not
/// crease, as for the sweeps. The square its chart maps onto takes its
/// corners first where the face across the region's boundary changes; the
/// faces there are read as the planar pass would make them, every region a
/// face of its own and the other free triangles grouped into planes.
/// Regions joined by the pieces they enclose ([`enclosed_unions`]) are
/// tried as one region first. A region refused a patch on its own keeps its
/// triangles for the planes, and is counted by why.
pub(super) fn patch_regions(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    options: &MeshSolidOptions,
    flat: f64,
    groups: &mut Groups,
    tol: Tolerances,
) -> OgeomResult<()> {
    let normals: Vec<Vector> = triangles.iter().map(|t| unit_normal(points, *t)).collect();
    // A free-form surface running on tangentially into a plane crosses no
    // crease, and the smooth region would take the plane with it. The
    // planes the free triangles gather into, where they hold a smooth
    // region's worth of triangles, bound the smooth regions instead.
    let mut held = groups.clone();
    {
        let mut planes = groups.clone();
        coplanar_groups(
            points,
            triangles,
            adjacency,
            options,
            flat,
            &mut planes,
            tol,
        )?;
        let mut sizes = vec![0_usize; planes.carriers.len()];
        for &g in &planes.of {
            if g != usize::MAX {
                sizes[g] += 1;
            }
        }
        for (t, &g) in planes.of.iter().enumerate() {
            if held.of[t] == usize::MAX && sizes[g] >= SWEPT_TRIANGLES {
                held.of[t] = g;
            }
        }
        held.carriers = planes.carriers;
    }
    let regions = smooth_regions(triangles, adjacency, &normals, options, &held);
    if regions.is_empty() {
        return Ok(());
    }
    let mut across = groups.clone();
    for region in &regions {
        let g = across.carriers.len();
        across.carriers.push(Carrier::Gone);
        for &t in region {
            across.of[t] = g;
        }
    }
    coplanar_groups(
        points,
        triangles,
        adjacency,
        options,
        flat,
        &mut across,
        tol,
    )?;
    // On a fine mesh a free-form surface is cut into small canonical pieces
    // before this pass, and the regions left round them are no disks. The
    // pieces the smooth area encloses are offered with the regions they
    // join, as one region, kept only where its patch verifies.
    let mut done = vec![false; regions.len()];
    for union in enclosed_unions(adjacency, &normals, options, groups, &regions) {
        let mut members: Vec<usize> = union
            .regions
            .iter()
            .flat_map(|&r| regions[r].iter().copied())
            .collect();
        members.extend_from_slice(&union.pieces);
        if let Ok(claim) = region_patch(points, triangles, adjacency, &across, &members, flat, tol)
        {
            for &b in &union.carriers {
                groups.carriers[b] = Carrier::Gone;
            }
            claim_swept(groups, &members, claim);
            for &r in &union.regions {
                done[r] = true;
            }
        }
    }
    for (region, _) in regions.iter().zip(&done).filter(|(_, d)| !**d) {
        match region_patch(points, triangles, adjacency, &across, region, flat, tol) {
            Ok(claim) => claim_swept(groups, region, claim),
            Err(crate::recognize_patch::Refused::NotDisk) => groups.refused.not_disk += 1,
            Err(crate::recognize_patch::Refused::Narrow) => groups.refused.narrow += 1,
            Err(crate::recognize_patch::Refused::Unverified) => groups.refused.unverified += 1,
        }
    }
    Ok(())
}

/// Smooth regions joined by the pieces between them into one region to
/// try a patch on.
struct EnclosedUnion {
    /// The smooth regions, by index.
    regions: Vec<usize>,
    /// The pieces' triangles.
    pieces: Vec<usize>,
    /// The carriers of the curved regions among the pieces.
    carriers: Vec<usize>,
}

/// The pieces the smooth regions enclose, gathered with the smooth regions
/// they join; only unions holding a curved region are returned.
///
/// A piece is a curved region not yet a patch, or a pocket of free
/// triangles too small to be a smooth region of its own. It is enclosed
/// when it meets the smooth regions or other enclosed pieces across an
/// edge that does not crease, and meets nothing else except across creases
/// and free edges. A curved region tangent to a face outside the smooth
/// area (a round beside a plane) is never enclosed, nor one holding more
/// triangles than the smooth regions it meets (a cylinder a hill rises
/// from), which bounds them.
fn enclosed_unions(
    adjacency: &Adjacency,
    normals: &[Vector],
    options: &MeshSolidOptions,
    groups: &Groups,
    regions: &[Vec<usize>],
) -> Vec<EnclosedUnion> {
    let cos_crease = options.crease.cos();
    let n = groups.of.len();
    let mut region_of = vec![usize::MAX; n];
    for (r, region) in regions.iter().enumerate() {
        for &t in region {
            region_of[t] = r;
        }
    }
    // The triangles across a triangle's edges that do not crease.
    let smooth = |t: usize| {
        (3 * t..3 * t + 3).filter_map(move |h| {
            adjacency.twin[h]
                .filter(|g| normals[h / 3].dot(normals[g / 3]) >= cos_crease)
                .map(|g| g / 3)
        })
    };
    // The pieces: each its carrier where it is a curved region, and its
    // triangles.
    let mut piece_of = vec![usize::MAX; n];
    let mut pieces: Vec<(Option<usize>, Vec<usize>)> = Vec::new();
    let mut by_carrier: FastMap<usize, usize> = FastMap::default();
    for (t, &g) in groups.of.iter().enumerate() {
        if matches!(groups.carriers.get(g), Some(Carrier::Curved(c)) if c.patch.is_none()) {
            let p = *by_carrier.entry(g).or_insert_with(|| {
                pieces.push((Some(g), Vec::new()));
                pieces.len() - 1
            });
            pieces[p].1.push(t);
            piece_of[t] = p;
        }
    }
    let free = |t: usize| groups.of[t] == usize::MAX && region_of[t] == usize::MAX;
    let mut seen = vec![false; n];
    for seed in 0..n {
        if !free(seed) || seen[seed] {
            continue;
        }
        let mut pocket = vec![seed];
        seen[seed] = true;
        let mut i = 0;
        while i < pocket.len() {
            for o in smooth(pocket[i]) {
                if free(o) && !seen[o] {
                    seen[o] = true;
                    pocket.push(o);
                }
            }
            i += 1;
        }
        if pocket.len() < SWEPT_TRIANGLES {
            for &t in &pocket {
                piece_of[t] = pieces.len();
            }
            pieces.push((None, pocket));
        }
    }
    // A curved region larger than the smooth regions it meets is a face of
    // its own that they run into (a hill on a cylinder's side), not a piece
    // cut from them: it bounds them, and is never enclosed.
    let bounding: Vec<bool> = pieces
        .iter()
        .map(|(carrier, triangles)| {
            let mut met: Vec<usize> = triangles
                .iter()
                .flat_map(|&t| smooth(t))
                .map(|o| region_of[o])
                .filter(|&r| r != usize::MAX)
                .collect();
            met.sort_unstable();
            met.dedup();
            let beside: usize = met.iter().map(|&r| regions[r].len()).sum();
            carrier.is_some() && !met.is_empty() && triangles.len() > beside
        })
        .collect();
    // Every piece reached from the smooth regions through smooth edges,
    // then those meeting anything else dropped until none does.
    let mut enclosed = vec![false; pieces.len()];
    let mut queue: Vec<usize> = Vec::new();
    for region in regions {
        for &t in region {
            for o in smooth(t) {
                let p = piece_of[o];
                if p != usize::MAX && !enclosed[p] && !bounding[p] {
                    enclosed[p] = true;
                    queue.push(p);
                }
            }
        }
    }
    while let Some(p) = queue.pop() {
        for &t in &pieces[p].1 {
            for o in smooth(t) {
                let q = piece_of[o];
                if q != usize::MAX && !enclosed[q] && !bounding[q] {
                    enclosed[q] = true;
                    queue.push(q);
                }
            }
        }
    }
    loop {
        let open: Vec<usize> = (0..pieces.len())
            .filter(|&p| {
                enclosed[p]
                    && pieces[p].1.iter().any(|&t| {
                        smooth(t).any(|o| {
                            region_of[o] == usize::MAX
                                && (piece_of[o] == usize::MAX || !enclosed[piece_of[o]])
                        })
                    })
            })
            .collect();
        if open.is_empty() {
            break;
        }
        for p in open {
            enclosed[p] = false;
        }
    }
    // The smooth regions and enclosed pieces joined through smooth edges.
    let mut seen_piece = vec![false; pieces.len()];
    let mut seen_region = vec![false; regions.len()];
    let mut unions = Vec::new();
    for start in 0..pieces.len() {
        if !enclosed[start] || seen_piece[start] {
            continue;
        }
        seen_piece[start] = true;
        let (mut joined, mut held) = (Vec::new(), vec![start]);
        // Each entry a piece (`Ok`) or a smooth region (`Err`) still to
        // search from.
        let mut queue: Vec<Result<usize, usize>> = vec![Ok(start)];
        while let Some(next) = queue.pop() {
            let from: &[usize] = match next {
                Ok(p) => &pieces[p].1,
                Err(r) => &regions[r],
            };
            for &t in from {
                for o in smooth(t) {
                    let (r, p) = (region_of[o], piece_of[o]);
                    if r != usize::MAX {
                        if !seen_region[r] {
                            seen_region[r] = true;
                            joined.push(r);
                            queue.push(Err(r));
                        }
                    } else if p != usize::MAX && enclosed[p] && !seen_piece[p] {
                        seen_piece[p] = true;
                        held.push(p);
                        queue.push(Ok(p));
                    }
                }
            }
        }
        let carriers: Vec<usize> = held.iter().filter_map(|&p| pieces[p].0).collect();
        if !joined.is_empty() && !carriers.is_empty() {
            unions.push(EnclosedUnion {
                regions: joined,
                pieces: held
                    .iter()
                    .flat_map(|&p| pieces[p].1.iter().copied())
                    .collect(),
                carriers,
            });
        }
    }
    unions
}

/// A patch fitted to a set of triangles, the square's corners first where
/// the face across its boundary (as read in `across`) changes.
pub(super) fn region_patch(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    across: &Groups,
    region: &[usize],
    flat: f64,
    tol: Tolerances,
) -> Result<Curved, crate::recognize_patch::Refused> {
    let inside: ogeom_core::FastSet<usize> = region.iter().copied().collect();
    // The faces across the boundary on either side of each boundary
    // vertex: a corner where they differ.
    let mut sides: FastMap<u32, [Vec<usize>; 2]> = FastMap::default();
    for &t in region {
        for h in 3 * t..3 * t + 3 {
            let other = match adjacency.twin[h] {
                Some(g) if inside.contains(&(g / 3)) => continue,
                Some(g) => across.of[g / 3],
                None => usize::MAX,
            };
            let (a, b) = from_to(triangles, h);
            sides.entry(a).or_default()[0].push(other);
            sides.entry(b).or_default()[1].push(other);
        }
    }
    let corners: ogeom_core::FastSet<u32> = sides
        .into_iter()
        .filter(|(_, [leaving, arriving])| {
            let mut faces = leaving.iter().chain(arriving);
            faces.next().is_some_and(|first| faces.any(|f| f != first))
        })
        .map(|(v, _)| v)
        .collect();
    let patch = crate::recognize_patch::fit_patch(
        &crate::recognize_patch::Region {
            points,
            triangles,
            members: region,
            corners: &corners,
        },
        flat,
        tol,
    )?;
    let mut vertices: Vec<u32> = region.iter().flat_map(|&t| triangles[t]).collect();
    vertices.sort_unstable();
    vertices.dedup();
    Ok(Curved {
        shape: Canonical::Swept(Box::new(crate::recognize::SweptShape::new(
            patch.surface,
            tol,
        ))),
        deviation: patch.deviation,
        fitted: patch.deviation,
        centre: (0.5, 0.5),
        wraps: false,
        wraps_v: false,
        fixed: true,
        vertices,
        patch: Some(patch.mapped),
    })
}

/// A region's triangles given to a new swept face.
fn claim_swept(groups: &mut Groups, region: &[usize], claim: Curved) {
    let g = groups.carriers.len();
    for &t in region {
        groups.of[t] = g;
    }
    groups.carriers.push(Carrier::Curved(claim));
}

/// The fewest triangles a smooth region needs to be tried as a sweep.
const SWEPT_TRIANGLES: usize = 32;

/// The largest distance between two of the points, from the first.
fn span(points: &[Point]) -> f64 {
    points.first().map_or(0.0, |a| {
        points.iter().map(|p| p.distance(*a)).fold(0.0, f64::max)
    })
}

/// Turn each recognized sphere whose boundary is circles in parallel planes
/// about their common normal, so the circles are latitudes and the face a
/// cap or a zone about a pole; the axis points into a cap, so its pole is
/// the frame's north one. A sphere with no boundary is whole and keeps its
/// frame, and one bounded otherwise keeps it for the re-framing that moves
/// its poles away.
pub(super) fn sphere_axes(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    flat: f64,
    tol: Tolerances,
) {
    for g in 0..groups.carriers.len() {
        sphere_axis(points, triangles, adjacency, groups, g, flat, tol);
    }
}

/// [`sphere_axes`] for the region `g`.
pub(super) fn sphere_axis(
    points: &[Point],
    triangles: &[[u32; 3]],
    adjacency: &Adjacency,
    groups: &mut Groups,
    g: usize,
    flat: f64,
    tol: Tolerances,
) {
    let reach = flat * REACH;
    let Carrier::Curved(curved) = &groups.carriers[g] else {
        return;
    };
    let Canonical::Sphere(sphere) = curved.shape else {
        return;
    };
    let Some(loops) = border_loops(triangles, adjacency, &groups.of, g) else {
        return;
    };
    let mut axis: Option<Vector> = None;
    let mut planar = true;
    for ring in &loops {
        let pts: Vec<Point> = ring.iter().map(|&v| points[v as usize]).collect();
        let Some((through, normal)) = (pts.len() >= 3).then(|| plane_through(&pts, tol)).flatten()
        else {
            planar = false;
            break;
        };
        let n = normal.vector();
        if pts.iter().any(|p| (*p - through).dot(n).abs() > reach) {
            planar = false;
            break;
        }
        match axis {
            None => axis = Some(n),
            Some(a) if a.cross(n).magnitude() <= 1e-3 => {}
            Some(_) => {
                planar = false;
                break;
            }
        }
    }
    let (true, Some(mut z)) = (planar, axis) else {
        return;
    };
    // Into the region: the pole a cap covers is the north one.
    let side: f64 = curved
        .vertices
        .iter()
        .map(|&v| (points[v as usize] - sphere.centre()).dot(z))
        .sum();
    if side < 0.0 {
        z = -z;
    }
    let Ok(z) = Direction::new(z, tol) else {
        return;
    };
    let Ok(frame) = Frame::new(sphere.centre(), z, z.any_perpendicular(), tol) else {
        return;
    };
    let Ok(turned) = Sphere::new(frame, sphere.radius(), tol) else {
        return;
    };
    if let Carrier::Curved(curved) = &mut groups.carriers[g] {
        curved.shape = Canonical::Sphere(turned);
        curved.fixed = true;
    }
}
