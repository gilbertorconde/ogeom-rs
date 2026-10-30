//! Where a solid's boundary passes through itself.
//!
//! A boundary that crosses itself bounds no single region: a point between
//! the two sheets is inside by one and outside by the other, and which a
//! probe meets decides what every question about the solid answers. The
//! faces are meshed, the triangles that cross each other are found by their
//! boxes, and each crossing is only a guess until it is settled on the
//! surfaces themselves: a point on both, where they meet at an angle,
//! inside both faces' trims and clear of their boundaries. A mesh can cross
//! where the surfaces do not (two walls closer than the chord), and two
//! faces always meet along the edges they share; neither is reported.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Surface as _, SurfaceGeometry};
use ogeom_math::{Aabb, Point, Transform, Vector};
use ogeom_mesh::Deflection;
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

use crate::classify::{Containment, FaceTrim};

/// The smallest sine of the angle two sheets cross at to count: below it
/// they touch, which is where faces meet along an edge they share or bend
/// smoothly into one another.
pub(crate) const CROSSING_SINE: f64 = 1e-2;

/// How many crossings of one pair of faces' meshes are asked of the
/// surfaces before the pair is let go: a pair whose meshes cross
/// everywhere without the surfaces doing so meets its walls closer than
/// the chord, and asking every triangle would cost without telling more.
const ASKED_PER_PAIR: usize = 16;

/// Below this many triangles a node of the box tree holds them as they are.
const LEAF: usize = 8;

/// One face of the solid the crossing test reads.
struct Sheet<'m> {
    face: &'m Shape,
    surface: &'m SurfaceGeometry,
    placement: Transform,
    /// How far from its trim a point counts as on it: the widest of the
    /// face's tolerance, its edges' and its vertices'. Edges meet at a
    /// vertex only within the vertex's tolerance, so near a loose vertex
    /// the drawn trim can stand that far from where the face ends.
    reach: f64,
    trim: std::sync::OnceLock<Option<FaceTrim<'m>>>,
}

/// One triangle of a face's mesh.
struct Facet {
    sheet: usize,
    corners: [Point; 3],
    chart: [(f64, f64); 3],
    bound: Aabb,
}

/// Two crossing triangles, by index, and a point where their meshes cross.
type Guess = (usize, usize, Point);

/// A point where a solid's boundary passes through itself.
#[derive(Debug, Clone)]
pub struct SelfCrossing {
    /// The face crossing.
    pub face: Shape,
    /// The face it crosses: the same one where a face folds through itself.
    pub other: Shape,
    /// Where, on both surfaces.
    pub at: Point,
}

/// The pairs of faces of `solid` that pass through each other, a face
/// folding through itself among them, one point for each: what
/// [`check`](crate::check()) reports as a boundary crossing itself.
///
/// # Errors
///
/// [`OgeomError::Cancelled`](ogeom_core::OgeomError::Cancelled) and
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) only.
pub fn self_crossings(
    model: &Model,
    solid: &Shape,
    tol: Tolerances,
) -> OgeomResult<Vec<SelfCrossing>> {
    crossings(model, solid, &ogeom_mesh::FaceMeshCache::default(), tol)
}

/// The pairs of faces of `solid` (a face with itself among them) that
/// cross each other, one point for each, in the order of the faces.
///
/// A face that will not mesh is not asked: what the test cannot draw it
/// cannot settle, and names nothing on a guess.
///
/// # Errors
///
/// [`OgeomError::Cancelled`](ogeom_core::OgeomError::Cancelled) and
/// [`OgeomError::Dangling`](ogeom_core::OgeomError::Dangling) only.
pub(crate) fn crossings(
    model: &Model,
    solid: &Shape,
    kept: &ogeom_mesh::FaceMeshCache,
    tol: Tolerances,
) -> OgeomResult<Vec<SelfCrossing>> {
    let faces = explore_unique(model, solid, ShapeType::Face)?;
    let meshed = ogeom_core::parallel::map_ordered(&faces, |_, face| {
        ogeom_mesh::triangulate_face_kept(model, face, Deflection::default(), kept, tol).ok()
    });
    let mut sheets = Vec::with_capacity(faces.len());
    let mut facets = Vec::new();
    for (face, mesh) in faces.iter().zip(meshed) {
        let Some(data) = model.node(face).and_then(|n| n.data().as_face()) else {
            continue;
        };
        let Some(surface) = model.geometry().surface(data.surface) else {
            continue;
        };
        let Some(mesh) = mesh.filter(|m| m.parameters.len() == m.positions.len()) else {
            continue;
        };
        let mut reach = data.tolerance.get().max(tol.confusion());
        for kind in [ShapeType::Edge, ShapeType::Vertex] {
            for bound in explore_unique(model, face, kind)? {
                if let Some(t) = model.node(&bound).and_then(|n| n.data().tolerance()) {
                    reach = reach.max(t.get());
                }
            }
        }
        let sheet = sheets.len();
        sheets.push(Sheet {
            face,
            surface,
            placement: face.transform(model.datums())?,
            reach,
            trim: std::sync::OnceLock::new(),
        });
        for triangle in &mesh.triangles {
            let corners = triangle.map(|i| mesh.positions[i as usize]);
            facets.push(Facet {
                sheet,
                corners,
                chart: triangle.map(|i| mesh.parameters[i as usize]),
                bound: Aabb::of_points(&corners).expanded(tol.confusion()),
            });
        }
    }
    if facets.len() < 2 {
        return Ok(Vec::new());
    }

    // The crossing triangles, a few per pair of faces, found in runs of
    // triangles so the search shares out.
    let tree = BoxTree::new(&facets);
    let runs: Vec<(usize, usize)> = (0..facets.len())
        .step_by(512)
        .map(|start| (start, (start + 512).min(facets.len())))
        .collect();
    let found = ogeom_core::parallel::map_ordered(&runs, |_, &(start, end)| {
        let mut out: Vec<(usize, usize, usize, usize, Point)> = Vec::new();
        let mut near = Vec::new();
        for i in start..end {
            if i % 64 == 0 {
                ogeom_core::progress::checkpoint()?;
            }
            tree.meeting(&facets, &facets[i].bound, &mut near);
            for &j in &near {
                if j <= i || share_a_corner(&facets[i], &facets[j]) {
                    continue;
                }
                if let Some(at) = triangles_cross(&facets[i].corners, &facets[j].corners) {
                    let (a, b) = (facets[i].sheet, facets[j].sheet);
                    if a <= b {
                        out.push((a, b, i, j, at));
                    } else {
                        out.push((b, a, j, i, at));
                    }
                }
            }
        }
        Ok(out)
    });
    let mut guesses = Vec::new();
    for run in found {
        guesses.extend(run?);
    }
    guesses.sort_by_key(|x| (x.0, x.1, x.2, x.3));
    let mut pairs: Vec<((usize, usize), Vec<Guess>)> = Vec::new();
    for (a, b, i, j, at) in guesses {
        match pairs.last_mut() {
            Some((pair, all)) if *pair == (a, b) => all.push((i, j, at)),
            _ => pairs.push(((a, b), vec![(i, j, at)])),
        }
    }
    for ((a, b), all) in &mut pairs {
        if a == b {
            // A face's mesh crosses itself where the face folds back on
            // itself sharply as well as where it passes through itself,
            // and there the two triangles stand side by side in the chart
            // and the surfaces meet only where they are one point. The
            // triangles farthest apart in the chart are asked.
            let surface = sheets[*a].surface;
            let mut apart: Vec<(f64, Guess)> = all
                .iter()
                .map(|&(i, j, at)| (chart_apart(surface, &facets[i], &facets[j]), (i, j, at)))
                .collect();
            apart.sort_by(|x, y| y.0.total_cmp(&x.0));
            *all = apart
                .into_iter()
                .take(ASKED_PER_PAIR)
                .map(|(_, guess)| guess)
                .collect();
        } else if all.len() > ASKED_PER_PAIR {
            // Spread along the crossing, rather than crowded at one end.
            *all = (0..ASKED_PER_PAIR)
                .map(|k| all[k * all.len() / ASKED_PER_PAIR])
                .collect();
        }
    }

    let settled = ogeom_core::parallel::map_ordered(&pairs, |_, ((a, b), asked)| {
        for &(i, j, guess) in asked {
            ogeom_core::progress::checkpoint()?;
            if let Some(at) = settle(model, &sheets, &facets[i], &facets[j], guess, tol)? {
                return Ok(Some(SelfCrossing {
                    face: sheets[*a].face.clone(),
                    other: sheets[*b].face.clone(),
                    at,
                }));
            }
        }
        Ok(None)
    });
    let mut out = Vec::new();
    for one in settled {
        out.extend(one?);
    }
    Ok(out)
}

/// Whether two triangles have a corner at the same place: neighbours in one
/// face's mesh, or across an edge two faces share, whose boundary points
/// both draw from the edge's own curve.
fn share_a_corner(p: &Facet, q: &Facet) -> bool {
    p.corners.iter().any(|c| q.corners.contains(c))
}

/// A point on the segment along which two triangles cross, where each
/// passes through the other's plane; `None` where they only touch, lie in
/// one plane, or miss.
fn triangles_cross(p: &[Point; 3], q: &[Point; 3]) -> Option<Point> {
    let plane = |t: &[Point; 3]| (t[1] - t[0]).cross(t[2] - t[0]);
    let (np, nq) = (plane(p), plane(q));
    let side = |n: Vector, origin: Point, t: &[Point; 3]| t.map(|c| n.dot(c - origin));
    let dp = side(nq, q[0], p);
    let dq = side(np, p[0], q);
    let straddles = |d: [f64; 3]| d.iter().any(|&x| x > 0.0) && d.iter().any(|&x| x < 0.0);
    if !straddles(dp) || !straddles(dq) {
        return None;
    }
    let line = np.cross(nq);
    if line.dot(line) <= f64::MIN_POSITIVE {
        return None;
    }
    // Where each triangle's edges pass through the other's plane, as an
    // interval along the line the two planes share.
    let interval = |t: &[Point; 3], d: [f64; 3]| -> Option<(f64, f64)> {
        let mut along: Vec<f64> = Vec::with_capacity(3);
        for k in 0..3 {
            let (a, b) = (k, (k + 1) % 3);
            if (d[a] > 0.0) != (d[b] > 0.0) && d[a] != d[b] {
                let s = d[a] / (d[a] - d[b]);
                let at = t[a] + (t[b] - t[a]) * s;
                along.push(line.dot(at - Point::ORIGIN));
            }
        }
        let lo = along.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = along.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        (lo <= hi).then_some((lo, hi))
    };
    let (pl, ph) = interval(p, dp)?;
    let (ql, qh) = interval(q, dq)?;
    let (lo, hi) = (pl.max(ql), ph.min(qh));
    if lo >= hi {
        return None;
    }
    // The middle of the shared stretch, on the line: a point of the first
    // triangle's crossing segment.
    let middle = 0.5 * (lo + hi);
    let mut ends: Vec<Point> = Vec::with_capacity(2);
    for k in 0..3 {
        let (a, b) = (k, (k + 1) % 3);
        if (dp[a] > 0.0) != (dp[b] > 0.0) && dp[a] != dp[b] {
            let s = dp[a] / (dp[a] - dp[b]);
            ends.push(p[a] + (p[b] - p[a]) * s);
        }
    }
    let [from, to] = [*ends.first()?, *ends.last()?];
    let (sf, st) = (line.dot(from - Point::ORIGIN), line.dot(to - Point::ORIGIN));
    if (st - sf).abs() <= f64::MIN_POSITIVE {
        return Some(from);
    }
    Some(from + (to - from) * ((middle - sf) / (st - sf)))
}

/// How far apart two triangles of one face stand in its chart, as a
/// fraction of the chart's extent along the direction they are farthest
/// apart in, the short way round a direction that repeats.
fn chart_apart(surface: &SurfaceGeometry, p: &Facet, q: &Facet) -> f64 {
    let middle = |f: &Facet| {
        let [a, b, c] = f.chart;
        ((a.0 + b.0 + c.0) / 3.0, (a.1 + b.1 + c.1) / 3.0)
    };
    let ((u0, u1), (v0, v1)) = surface.domain();
    let along = |d: f64, span: f64, periodic: bool| {
        if !(span.is_finite() && span > 0.0) {
            return d.abs();
        }
        let d = if periodic {
            let wrapped = d.rem_euclid(span);
            wrapped.min(span - wrapped)
        } else {
            d.abs()
        };
        d / span
    };
    let (pm, qm) = (middle(p), middle(q));
    along(pm.0 - qm.0, u1 - u0, surface.is_periodic_u()).max(along(
        pm.1 - qm.1,
        v1 - v0,
        surface.is_periodic_v(),
    ))
}

/// The chart point of `facet` under `at`, by the barycentric weights of
/// its foot in the triangle's plane.
fn chart_under(facet: &Facet, at: Point) -> Option<(f64, f64)> {
    let [a, b, c] = facet.corners;
    let (e1, e2, w) = (b - a, c - a, at - a);
    let (d11, d12, d22) = (e1.dot(e1), e1.dot(e2), e2.dot(e2));
    let (d1w, d2w) = (e1.dot(w), e2.dot(w));
    let det = d11.mul_add(d22, -d12 * d12);
    if det.abs() <= f64::MIN_POSITIVE {
        return None;
    }
    let s = d22.mul_add(d1w, -d12 * d2w) / det;
    let t = d11.mul_add(d2w, -d12 * d1w) / det;
    let r = 1.0 - s - t;
    let [pa, pb, pc] = facet.chart;
    Some((
        r.mul_add(pa.0, s.mul_add(pb.0, t * pc.0)),
        r.mul_add(pa.1, s.mul_add(pb.1, t * pc.1)),
    ))
}

/// A sheet's point, its two chart tangents and its unit normal, placed.
fn placed_at(
    sheet: &Sheet<'_>,
    at: (f64, f64),
    tol: Tolerances,
) -> Option<(Point, Vector, Vector, Vector)> {
    let (p, du, dv) = sheet.surface.point_d1_at(at.0, at.1, tol).ok()?;
    let (du, dv) = (
        sheet.placement.apply_vector(du),
        sheet.placement.apply_vector(dv),
    );
    let n = du.cross(dv);
    let length = n.magnitude();
    if length <= f64::MIN_POSITIVE || !length.is_finite() {
        return None;
    }
    Some((sheet.placement.apply(p), du, dv, n / length))
}

/// Hold a chart point inside the surface's domain along its directions
/// that do not repeat.
fn held(surface: &SurfaceGeometry, u: f64, v: f64) -> (f64, f64) {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let u = if surface.is_periodic_u() || u0 > u1 {
        u
    } else {
        u.clamp(u0, u1)
    };
    let v = if surface.is_periodic_v() || v0 > v1 {
        v
    } else {
        v.clamp(v0, v1)
    };
    (u, v)
}

/// The point where two crossing triangles' surfaces really cross, near
/// `guess`; `None` where they do not cross there.
///
/// Newton's method from the triangles' chart points onto a point of both
/// surfaces, held from sliding along the curve they meet in. The point
/// counts where the surfaces cross at an angle there, it lies near the
/// triangles, and it is inside both faces, clear of their boundaries by
/// their tolerances.
fn settle<'m>(
    model: &'m Model,
    sheets: &[Sheet<'m>],
    p: &Facet,
    q: &Facet,
    guess: Point,
    tol: Tolerances,
) -> OgeomResult<Option<Point>> {
    let (a, b) = (&sheets[p.sheet], &sheets[q.sheet]);
    let (Some(start_a), Some(start_b)) = (chart_under(p, guess), chart_under(q, guess)) else {
        return Ok(None);
    };
    let (Some((_, _, _, na)), Some((_, _, _, nb))) =
        (placed_at(a, start_a, tol), placed_at(b, start_b, tol))
    else {
        return Ok(None);
    };
    let along = na.cross(nb);
    if along.magnitude() < CROSSING_SINE {
        return Ok(None);
    }
    let along = along / along.magnitude();
    let system = |x: &[f64; 4]| {
        let (ua, va) = held(a.surface, x[0], x[1]);
        let (ub, vb) = held(b.surface, x[2], x[3]);
        let (Some((pa, au, av, _)), Some((pb, bu, bv, _))) =
            (placed_at(a, (ua, va), tol), placed_at(b, (ub, vb), tol))
        else {
            return ([f64::INFINITY; 4], [[0.0; 4]; 4]);
        };
        let gap = pa - pb;
        (
            [gap.x, gap.y, gap.z, (pa - guess).dot(along)],
            [
                [au.x, av.x, -bu.x, -bv.x],
                [au.y, av.y, -bu.y, -bv.y],
                [au.z, av.z, -bu.z, -bv.z],
                [au.dot(along), av.dot(along), 0.0, 0.0],
            ],
        )
    };
    let criteria = ogeom_math::solve::Criteria {
        residual: tol.confusion() * 0.01,
        step: tol.parametric(),
        max_iterations: 40,
    };
    let start = [start_a.0, start_a.1, start_b.0, start_b.1];
    let Ok((x, residual, ..)) = ogeom_math::solve::newton_system_fixed(system, start, criteria)
    else {
        return Ok(None);
    };
    if residual > tol.confusion() {
        return Ok(None);
    }
    let (on_a, on_b) = (held(a.surface, x[0], x[1]), held(b.surface, x[2], x[3]));
    let (Some((at, _, _, na)), Some((_, _, _, nb))) =
        (placed_at(a, on_a, tol), placed_at(b, on_b, tol))
    else {
        return Ok(None);
    };
    if na.cross(nb).magnitude() < CROSSING_SINE {
        return Ok(None);
    }
    let size = |f: &Facet| {
        let [x, y, z] = f.corners;
        x.distance(y).max(y.distance(z)).max(z.distance(x))
    };
    if at.distance(guess) > 2.0 * size(p).max(size(q)) {
        return Ok(None);
    }
    for (sheet, chart) in [(a, on_a), (b, on_b)] {
        let Some(trim) = trim_of(model, sheet, tol) else {
            return Ok(None);
        };
        if trim.holds(model, chart, sheet.reach, tol)? != Containment::In {
            return Ok(None);
        }
    }
    Ok(Some(at))
}

/// A sheet's trim, drawn the first time a crossing asks about it; `None`
/// for a face whose rings will not draw.
fn trim_of<'s, 'm>(
    model: &'m Model,
    sheet: &'s Sheet<'m>,
    tol: Tolerances,
) -> Option<&'s FaceTrim<'m>> {
    sheet
        .trim
        .get_or_init(|| {
            // Drawn at the chord the inside-out probe holds faces to.
            let chord = Deflection {
                chord: tol.confusion() * 1e4,
                ..Deflection::default()
            };
            FaceTrim::of(model, sheet.face, chord, tol).ok()
        })
        .as_ref()
}

/// The facets' boxes held in a tree, so the boxes meeting one are found by
/// descending where it reaches.
struct BoxTree {
    nodes: Vec<Node>,
    order: Vec<usize>,
}

struct Node {
    bound: Aabb,
    run: (usize, usize),
    children: Option<(usize, usize)>,
}

impl BoxTree {
    fn new(facets: &[Facet]) -> Self {
        let mut order: Vec<usize> = (0..facets.len()).collect();
        let mut tree = Self {
            nodes: Vec::new(),
            order: Vec::new(),
        };
        let len = order.len();
        if len > 0 {
            tree.build(facets, &mut order, 0, len);
        }
        tree.order = order;
        tree
    }

    fn build(&mut self, facets: &[Facet], order: &mut [usize], start: usize, end: usize) -> usize {
        let bound = order[start..end]
            .iter()
            .fold(Aabb::EMPTY, |acc, &i| acc.union(&facets[i].bound));
        let at = self.nodes.len();
        self.nodes.push(Node {
            bound,
            run: (start, end),
            children: None,
        });
        if end - start <= LEAF {
            return at;
        }
        let size = bound.size();
        let axis = if size.x >= size.y && size.x >= size.z {
            0
        } else if size.y >= size.z {
            1
        } else {
            2
        };
        let centre = |i: usize| -> f64 {
            let c = facets[i].bound.centre().unwrap_or(Point::ORIGIN);
            [c.x, c.y, c.z][axis]
        };
        let middle = start + (end - start) / 2;
        order[start..end]
            .select_nth_unstable_by(middle - start, |&p, &q| centre(p).total_cmp(&centre(q)));
        let left = self.build(facets, order, start, middle);
        let right = self.build(facets, order, middle, end);
        self.nodes[at].children = Some((left, right));
        at
    }

    /// The facets whose boxes meet `probe`, ascending, into `out`.
    fn meeting(&self, facets: &[Facet], probe: &Aabb, out: &mut Vec<usize>) {
        out.clear();
        if self.nodes.is_empty() {
            return;
        }
        let mut stack = vec![0_usize];
        while let Some(at) = stack.pop() {
            let node = &self.nodes[at];
            if !node.bound.intersects(probe) {
                continue;
            }
            match node.children {
                Some((left, right)) => {
                    stack.push(left);
                    stack.push(right);
                }
                None => {
                    for &i in &self.order[node.run.0..node.run.1] {
                        if facets[i].bound.intersects(probe) {
                            out.push(i);
                        }
                    }
                }
            }
        }
        out.sort_unstable();
    }
}
