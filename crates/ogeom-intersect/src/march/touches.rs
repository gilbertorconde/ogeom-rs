//! Touches, where the surfaces meet with parallel normals, and the search for them.

use super::*;

/// A point where the two surfaces touch, their normals parallel, and the
/// section crosses itself.
#[derive(Debug, Clone, Copy)]
pub(super) struct Touch {
    pub(super) point: Point,
    pub(super) on_a: (f64, f64),
    pub(super) on_b: (f64, f64),
    /// How far the second surface stands off the first along its normal
    /// there, within the confusion distance.
    apart: f64,
    /// How far the solve for the touch may land from it: where a gap the
    /// size of its residual opens, `√(2r / |λ|)` for the least eigenvalue
    /// `λ` of `Q`.
    pub(super) blur: f64,
    /// The two lines the section crosses itself along, where they stand
    /// far enough apart to be followed one by one: see [`touch_arms`].
    pub(super) arms: Option<[Vector; 2]>,
    /// How far round the touch a sample counts as on it: see
    /// [`settle_touches`].
    pub(super) radius: f64,
    /// How far round the touch a stalled end or a step between its arms is
    /// taken as heading for it: see [`settle_touches`].
    pub(super) around: f64,
    /// How fast the surfaces part across an arm, `√(−det Q)` of
    /// [`touch_arms`]: a sample `ρ` out along an arm with the surfaces
    /// within `ε` of each other stands up to `ε / (ρ · spread)` off it.
    pub(super) spread: f64,
    /// How far from the touch a straight chord to it stays on both
    /// surfaces within the confusion distance `ε`: `√(2ε / k)`, `k` the
    /// largest bend of either surface or of their difference there.
    pub(super) snap: f64,
}

/// The two lines along which the section leaves a touch, as unit vectors
/// in the shared tangent plane.
///
/// Near the touch each surface stands off the tangent plane by half its
/// second fundamental form, so the section is where the two forms agree:
/// the null lines of their difference `Q`. Where `Q` is indefinite there
/// are two, and the section crosses itself along them; where it is
/// definite the touch is an isolated point, and where it is singular the
/// surfaces touch along a curve. `None` unless there are two; with them,
/// how fast the surfaces part across either, `√(−det Q)`, the largest bend
/// of either surface or of `Q`, the spectral radius of the three forms,
/// and the eigenvalues of `Q`, least first.
pub(super) fn touch_arms(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    on_a: (f64, f64),
    on_b: (f64, f64),
    tol: Tolerances,
) -> Option<([Vector; 2], f64, f64, [f64; 2])> {
    let n = normal_at(a, on_a, tol)?;
    let (au, _) = a.d1_at(on_a.0, on_a.1, tol).ok()?;
    let e1 = au - n * au.dot(n);
    let e1 = e1 / e1.magnitude();
    let e2 = n.cross(e1);
    // The second fundamental form of one surface along `n`, in the basis
    // (e1, e2): each basis vector written in the chart through the first
    // fundamental form, then the second derivatives read off along `n`.
    let form = |surface: &SurfaceGeometry, (u, v): (f64, f64)| -> Option<[f64; 3]> {
        let (su, sv) = surface.d1_at(u, v, tol).ok()?;
        let (suu, suv, svv) = surface.d2_at(u, v, tol).ok()?;
        let (e, f, g) = (su.dot(su), su.dot(sv), sv.dot(sv));
        let det = e.mul_add(g, -(f * f));
        if det <= f64::MIN_POSITIVE {
            return None;
        }
        let chart = |w: Vector| {
            let (p, q) = (su.dot(w), sv.dot(w));
            (g.mul_add(p, -(f * q)) / det, e.mul_add(q, -(f * p)) / det)
        };
        let (x, y) = (chart(e1), chart(e2));
        let second = |s: (f64, f64), t: (f64, f64)| {
            n.dot(suu) * s.0 * t.0
                + n.dot(suv) * s.0.mul_add(t.1, s.1 * t.0)
                + n.dot(svv) * s.1 * t.1
        };
        Some([second(x, x), second(x, y), second(y, y)])
    };
    let (fa, fb) = (form(a, on_a)?, form(b, on_b)?);
    let q = [fa[0] - fb[0], fa[1] - fb[1], fa[2] - fb[2]];
    // Along the line at angle θ to e1 the difference is A + B cos 2θ +
    // C sin 2θ = A + R cos(2θ − φ), zero at 2θ = φ ± acos(−A / R): two
    // lines where |A| < R, the form indefinite.
    let mean = f64::midpoint(q[0], q[2]);
    let (half_gap, skew) = ((q[0] - q[2]) * 0.5, q[1]);
    let swing = half_gap.hypot(skew);
    if mean.abs().partial_cmp(&swing) != Some(core::cmp::Ordering::Less) {
        return None;
    }
    let (phi, open) = (skew.atan2(half_gap), (-mean / swing).acos());
    let line = |theta: f64| e1 * theta.cos() + e2 * theta.sin();
    let arms = [line((phi + open) * 0.5), line((phi - open) * 0.5)];
    let bend = |m: [f64; 3]| f64::midpoint(m[0], m[2]).abs() + ((m[0] - m[2]) * 0.5).hypot(m[1]);
    Some((
        arms,
        swing.mul_add(swing, -(mean * mean)).sqrt(),
        bend(fa).max(bend(fb)).max(bend(q)),
        [mean - swing, mean + swing],
    ))
}

/// Whether the section crosses itself at a touch, not just passes it: the
/// surfaces stand `apart` along the normal there, so the section is where
/// half of `Q` makes that up, a hyperbola whose branches pass
/// `√(2 |apart| / |λ|)` from the touch, `λ` the eigenvalue of `Q` of the
/// sign of `apart`. Two surfaces within the confusion distance of each
/// other at a near tangency leave such a point where nothing crosses; the
/// section crosses at the touch where the branches pass within the chord
/// `near` of it, or where `apart` is within the residual the solve for the
/// touch settles to, so no different from none.
pub(super) fn crosses(
    touch: &Touch,
    (_, _, _, eigen): ([Vector; 2], f64, f64, [f64; 2]),
    near: f64,
    tol: Tolerances,
) -> bool {
    let lambda = if touch.apart >= 0.0 {
        eigen[1]
    } else {
        eigen[0]
    };
    touch.apart.abs() <= TOUCH_RESIDUAL * tol.confusion()
        || (lambda.abs() > 0.0 && 2.0 * touch.apart.abs() / lambda.abs() <= near * near)
}

/// How far a solve for a touch may land from it, from the eigenvalues of
/// `Q` there: see [`Touch::blur`].
pub(super) fn blur_of(eigen: [f64; 2], tol: Tolerances) -> f64 {
    let least = eigen[0].abs().min(eigen[1].abs());
    if least > 0.0 {
        (2.0 * TOUCH_RESIDUAL * tol.confusion() / least).sqrt()
    } else {
        0.0
    }
}

/// Whether two solved touches are one: they stand within the chord `near`,
/// or within where either solve may land.
pub(super) fn same_touch(t: &Touch, u: &Touch, near: f64) -> bool {
    t.point.distance(u.point) <= near.max(t.blur).max(u.blur)
}

/// Whether two arms stand far enough apart that a sample `radius` out
/// along one lies more than twice `near`, the march's chord, off the
/// other, so each can be told from the other and followed on its own.
pub(super) fn apart(arms: [Vector; 2], radius: f64, near: f64) -> bool {
    radius * arms[0].cross(arms[1]).magnitude() > 2.0 * near
}

/// The residual a solve for a touch settles to, as a fraction of the
/// confusion distance.
const TOUCH_RESIDUAL: f64 = 1e-3;

/// Where the two surfaces touch near a contact: on both, normals parallel.
///
/// Posed as a pair of points standing apart by `t` along the first
/// surface's normal, with the second surface's tangent plane square to that
/// normal: five equations in the four parameters and `t`, regular where the
/// surfaces bend apart differently across the touch, as they do where a
/// section crosses itself. A touch is a solution with `t` inside the
/// confusion distance; a pair standing further apart there is a near miss,
/// which the march follows as it is. `None` where the solve does not settle.
pub(super) fn touching_point(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    on_a: (f64, f64),
    on_b: (f64, f64),
    tol: Tolerances,
) -> Option<Touch> {
    let eval = |x: &[f64]| -> Option<[f64; 5]> {
        let (ua, va) = clamp(a, x[0], x[1]);
        let (ub, vb) = clamp(b, x[2], x[3]);
        let pa = a.point_at(ua, va, tol).ok()?;
        let pb = b.point_at(ub, vb, tol).ok()?;
        let na = normal_at(a, (ua, va), tol)?;
        let (bu, bv) = b.d1_at(ub, vb, tol).ok()?;
        let (lu, lv) = (bu.magnitude(), bv.magnitude());
        if lu <= tol.confusion() || lv <= tol.confusion() {
            return None;
        }
        let gap = pa + na * x[4] - pb;
        Some([gap.x, gap.y, gap.z, na.dot(bu) / lu, na.dot(bv) / lv])
    };
    // Forward differences: the Jacobian's own error slows the iteration and
    // moves no root. They are taken at the point the residual measures,
    // held to the charts, and step back into a chart where a step forward
    // would leave it: a touch on a bound that does not repeat (a profile
    // ending on its apex) draws the iteration onto the bound and past it,
    // and a step clamped back onto the bound would read that parameter as
    // moving nothing, the Jacobian as singular.
    const STEP: f64 = 1e-7;
    let step_from = |k: usize, at: f64| -> f64 {
        let surface = if k < 2 { a } else { b };
        let ((_, ub), (_, vb)) = surface.domain();
        let (periodic, high) = if matches!(k, 0 | 2) {
            (surface.is_periodic_u(), ub)
        } else {
            (surface.is_periodic_v(), vb)
        };
        if periodic || at + STEP <= high {
            STEP
        } else {
            -STEP
        }
    };
    let system = |x: &[f64; 5]| {
        let failed = || ([f64::INFINITY; 5], [[0.0; 5]; 5]);
        let mut at = *x;
        (at[0], at[1]) = clamp(a, x[0], x[1]);
        (at[2], at[3]) = clamp(b, x[2], x[3]);
        let Some(here) = eval(&at) else {
            return failed();
        };
        let mut jacobian = [[0.0; 5]; 5];
        for k in 0..5 {
            let step = if k < 4 { step_from(k, at[k]) } else { STEP };
            let mut moved = at;
            moved[k] += step;
            let Some(there) = eval(&moved) else {
                return failed();
            };
            for (row, (t, h)) in jacobian.iter_mut().zip(there.iter().zip(here.iter())) {
                row[k] = (t - h) / step;
            }
        }
        (here, jacobian)
    };
    let criteria = solve::Criteria {
        residual: tol.confusion() * TOUCH_RESIDUAL,
        step: tol.parametric() * 1e-3,
        max_iterations: 50,
    };
    let start = [on_a.0, on_a.1, on_b.0, on_b.1, 0.0];
    let (x, norm, _, _) = solve::newton_system_fixed(system, start, criteria).ok()?;
    if norm > tol.confusion() || x[4].abs() > tol.confusion() {
        return None;
    }
    let (on_a, on_b) = (clamp(a, x[0], x[1]), clamp(b, x[2], x[3]));
    let point = a.point_at(on_a.0, on_a.1, tol).ok()?;
    if b.point_at(on_b.0, on_b.1, tol).ok()?.distance(point) > tol.confusion() {
        return None;
    }
    Some(Touch {
        point,
        on_a,
        on_b,
        apart: x[4],
        blur: 0.0,
        arms: None,
        radius: 0.0,
        around: 0.0,
        spread: 0.0,
        snap: 0.0,
    })
}

/// A rectangle of a surface's chart, sampled on a 3 × 3 grid: where the
/// surface lies, within how far it bows from the samples, and which way it
/// faces.
#[derive(Debug, Clone, Copy)]
struct Patch {
    u: (f64, f64),
    v: (f64, f64),
    /// The box of the samples, grown by the bow.
    low: Point,
    high: Point,
    /// The normal at the middle, and the angle every normal of the patch
    /// stays within of it.
    axis: Vector,
    spread: f64,
    /// The diagonal of the box of the samples.
    size: f64,
    /// The normals at the samples.
    normals: [Vector; 9],
    /// Whether the patch is a piece of a plane: its samples in one plane
    /// and its normals one direction, so no quarter of it says more.
    flat: bool,
    /// Across `u` and across `v`: how far the samples run, and how far the
    /// normals turn.
    across: [(f64, f64); 2],
}

impl Patch {
    /// The patch over `u` × `v`, or `None` where the surface cannot be
    /// evaluated or has no normal there.
    fn over(
        surface: &SurfaceGeometry,
        u: (f64, f64),
        v: (f64, f64),
        tol: Tolerances,
    ) -> Option<Self> {
        let mut points = [[Point::ORIGIN; 3]; 3];
        let mut normals = [[Vector::new(0.0, 0.0, 0.0); 3]; 3];
        for (i, row) in points.iter_mut().enumerate() {
            for (j, point) in row.iter_mut().enumerate() {
                #[allow(clippy::cast_precision_loss)]
                let (s, t) = (i as f64 * 0.5, j as f64 * 0.5);
                let (pu, pv) = (u.0 + (u.1 - u.0) * s, v.0 + (v.1 - v.0) * t);
                *point = surface.point_at(pu, pv, tol).ok()?;
                normals[i][j] = normal_at(surface, (pu, pv), tol)?;
            }
        }
        // How far the surface bows from the samples: the middle samples'
        // distance from the straight blend of the ones either side, the
        // sag of a quadratic through them, taken twice over as the
        // seeding's cells take it.
        let mid = |p: Point, q: Point| p + (q - p) * 0.5;
        let mut bow: f64 = 0.0;
        for (k, row) in points.iter().enumerate() {
            bow = bow
                .max(points[1][k].distance(mid(points[0][k], points[2][k])))
                .max(row[1].distance(mid(row[0], row[2])));
        }
        let grow = BOW * bow;
        let mut low = Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
        let mut high = Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
        for p in points.iter().flatten() {
            low = Point::new(low.x.min(p.x), low.y.min(p.y), low.z.min(p.z));
            high = Point::new(high.x.max(p.x), high.y.max(p.y), high.z.max(p.z));
        }
        let pad = Vector::new(grow, grow, grow);
        let axis = normals[1][1];
        // The normals turn between samples by as much again as they turn
        // across them, so the spread is taken twice over too.
        let mut turn: f64 = 0.0;
        for n in normals.iter().flatten() {
            turn = turn.max(axis.dot(*n).clamp(-1.0, 1.0).acos());
        }
        let mut flat = [axis; 9];
        for (k, n) in normals.iter().flatten().enumerate() {
            flat[k] = *n;
        }
        let angle = |m: Vector, n: Vector| m.dot(n).clamp(-1.0, 1.0).acos();
        // Taken through the middle sample, as a patch over a whole closed
        // chart ends where it starts.
        let mut across = [(0.0_f64, 0.0_f64); 2];
        for k in 0..3 {
            let (p, n) = (&points, &normals);
            across[0].0 = across[0]
                .0
                .max(p[0][k].distance(p[1][k]) + p[1][k].distance(p[2][k]));
            across[0].1 = across[0]
                .1
                .max(angle(n[0][k], n[1][k]) + angle(n[1][k], n[2][k]));
            across[1].0 = across[1]
                .0
                .max(p[k][0].distance(p[k][1]) + p[k][1].distance(p[k][2]));
            across[1].1 = across[1]
                .1
                .max(angle(n[k][0], n[k][1]) + angle(n[k][1], n[k][2]));
        }
        Some(Self {
            u,
            v,
            low: low - pad,
            high: high + pad,
            axis,
            spread: BOW * turn,
            size: (high - low).magnitude(),
            normals: flat,
            flat: grow == 0.0 && turn == 0.0,
            across,
        })
    }

    /// The way the patch is halved next, if it still runs further than
    /// `near` any way, and how much that says: whether the normals turn
    /// that way, and by how much, or else how far it runs. Where two
    /// patches could touch, the normals that turn across a patch are what
    /// must part for the pair to be set aside, so the patch is halved
    /// across the way they turn most, a thin strip along a line of
    /// near-tangency parting from its neighbours without being cut along
    /// its length; the ways the normals do not turn are halved once those
    /// that do run within `near`.
    fn cut(&self, near: f64) -> Option<(usize, (bool, f64))> {
        if self.flat {
            return None;
        }
        let open = |k: usize| self.across[k].0 > near;
        // Halved across the way that says most, or quartered where the
        // other says at least half as much.
        let pick = |says: &dyn Fn(usize) -> f64| {
            let k = (0..2)
                .filter(|&k| open(k) && says(k) > 0.0)
                .max_by(|&x, &y| says(x).total_cmp(&says(y)))?;
            let other = 1 - k;
            let way = if open(other) && says(other) * 2.0 >= says(k) {
                2
            } else {
                k
            };
            Some((way, says(k)))
        };
        if let Some((way, turn)) = pick(&|k| self.across[k].1) {
            return Some((way, (true, turn)));
        }
        pick(&|k| self.across[k].0).map(|(way, run)| (way, (false, run)))
    }

    /// The two halves of the patch across `u` (`way` 0) or `v` (1), or its
    /// four quarters (2).
    fn halves(&self, way: usize, surface: &SurfaceGeometry, tol: Tolerances) -> Vec<Self> {
        let (um, vm) = self.middle();
        let (u0, u1) = ((self.u.0, um), (um, self.u.1));
        let (v0, v1) = ((self.v.0, vm), (vm, self.v.1));
        let parts = match way {
            0 => vec![(u0, self.v), (u1, self.v)],
            1 => vec![(self.u, v0), (self.u, v1)],
            _ => vec![(u0, v0), (u1, v0), (u0, v1), (u1, v1)],
        };
        parts
            .into_iter()
            .filter_map(|(u, v)| Self::over(surface, u, v, tol))
            .collect()
    }

    /// The middle of the patch in its chart.
    fn middle(&self) -> (f64, f64) {
        (
            f64::midpoint(self.u.0, self.u.1),
            f64::midpoint(self.v.0, self.v.1),
        )
    }
}

/// Whether no normal of one patch can be parallel to a normal of the
/// other although their cones meet: some component of the cross product
/// of the normals keeps one sign over every pair of samples, its value
/// at the middles further from zero than it strays across the samples,
/// taken twice over as the bow is.
fn normals_part(pa: &Patch, pb: &Patch) -> bool {
    let middle = pa.normals[4].cross(pb.normals[4]);
    let mut stray = [0.0_f64; 3];
    for na in &pa.normals {
        for nb in &pb.normals {
            let c = na.cross(*nb);
            stray[0] = stray[0].max((c.x - middle.x).abs());
            stray[1] = stray[1].max((c.y - middle.y).abs());
            stray[2] = stray[2].max((c.z - middle.z).abs());
        }
    }
    [middle.x, middle.y, middle.z]
        .iter()
        .zip(stray)
        .any(|(m, d)| m.abs() > BOW * d)
}

/// A surface's patches, and for each the indices of its halves once cut.
type PatchTree = (Vec<Patch>, Vec<Option<Vec<usize>>>);

/// The crossings where the two surfaces touch, found whatever the march's
/// stride: both charts are cut into patches, and each pair is set aside
/// where their boxes do not meet, or where no normal of one can be
/// parallel to a normal of the other, their normal cones apart; every
/// other pair is halved, across the way the normals turn most (see
/// [`Patch::cut`]), until both stand within the chord, and the touch is
/// solved from the middles of what is left. Only touches the section crosses itself at, with two arms, are
/// kept. `None` where the charts are not bounded or the pairs to weigh
/// would outnumber the points a branch may hold, so the search is not
/// complete.
pub(super) fn subdivided_touches(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
) -> Option<Vec<Touch>> {
    search_touches(a, b, options, tol, true).or_else(|| search_touches(a, b, options, tol, false))
}

/// The search of [`subdivided_touches`], each patch halved across the way
/// its normals turn most where `across` holds (see [`Patch::cut`]), and
/// otherwise quartered, the larger first, until it stands within the
/// chord. Halving parts a strip along a line of near-tangency from its
/// neighbours without cutting it along its length; where the normals turn
/// alike every way, as about many crossings close together, quartering
/// weighs fewer pairs, and the search is run so where halving does not
/// complete.
fn search_touches(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    options: Marching,
    tol: Tolerances,
    across: bool,
) -> Option<Vec<Touch>> {
    let near = options.chord.max(tol.confusion());
    let reach = near * 60.0;
    let whole = |surface: &SurfaceGeometry| {
        let (u, v) = surface.domain();
        [u.0, u.1, v.0, v.1]
            .iter()
            .all(|x| x.is_finite())
            .then(|| Patch::over(surface, u, v, tol))
            .flatten()
    };
    let cut_of = |p: &Patch| {
        if across {
            p.cut(near)
        } else {
            (!p.flat && p.size > near).then_some((2, (false, p.size)))
        }
    };
    // Each surface's patches as a tree, each halved once however many
    // patches of the other it is weighed against.
    let mut trees: [PatchTree; 2] = [(vec![whole(a)?], vec![None]), (vec![whole(b)?], vec![None])];
    let halve = |tree: &mut PatchTree, surface: &SurfaceGeometry, k: usize| -> Vec<usize> {
        if let Some(children) = &tree.1[k] {
            return children.clone();
        }
        let start = tree.0.len();
        let Some((way, _)) = cut_of(&tree.0[k]) else {
            return Vec::new();
        };
        let halves = tree.0[k].halves(way, surface, tol);
        let count = halves.len();
        tree.0.extend(halves);
        tree.1.extend(core::iter::repeat_with(|| None).take(count));
        let children: Vec<usize> = (start..start + count).collect();
        tree.1[k] = Some(children.clone());
        children
    };
    let mut pairs = vec![(0_usize, 0_usize)];
    let mut touches: Vec<Touch> = Vec::new();
    let mut weighed = 0_usize;
    while let Some((ia, ib)) = pairs.pop() {
        weighed += 1;
        if weighed > options.max_points {
            return None;
        }
        let (pa, pb) = (trees[0].0[ia], trees[1].0[ib]);
        let meet = pa.low.x <= pb.high.x
            && pb.low.x <= pa.high.x
            && pa.low.y <= pb.high.y
            && pb.low.y <= pa.high.y
            && pa.low.z <= pb.high.z
            && pb.low.z <= pa.high.z;
        if !meet {
            continue;
        }
        let facing = pa.axis.dot(pb.axis).abs().clamp(-1.0, 1.0).acos();
        if facing > pa.spread + pb.spread || normals_part(&pa, &pb) {
            continue;
        }
        let (cut_a, cut_b) = (cut_of(&pa), cut_of(&pb));
        if cut_a.is_none() && cut_b.is_none() {
            // A touch already found within the pair is the one there:
            // touches stand further apart than the chord.
            let held = touches.iter().any(|t| {
                t.point.x >= pa.low.x - near
                    && t.point.x <= pa.high.x + near
                    && t.point.y >= pa.low.y - near
                    && t.point.y <= pa.high.y + near
                    && t.point.z >= pa.low.z - near
                    && t.point.z <= pa.high.z + near
            });
            if held {
                continue;
            }
            // A flat patch is started from beside the other's middle.
            let (start_a, start_b) = match (pa.flat, pb.flat) {
                (true, false) => {
                    let at = b.point_at(pb.middle().0, pb.middle().1, tol).ok();
                    let on = at.and_then(|p| nearest_on(a, pa.middle(), p, tol));
                    (on.map_or(pa.middle(), |(uv, _)| uv), pb.middle())
                }
                (false, true) => {
                    let at = a.point_at(pa.middle().0, pa.middle().1, tol).ok();
                    let on = at.and_then(|p| nearest_on(b, pb.middle(), p, tol));
                    (pa.middle(), on.map_or(pb.middle(), |(uv, _)| uv))
                }
                _ => (pa.middle(), pb.middle()),
            };
            let Some(mut touch) = touching_point(a, b, start_a, start_b, tol) else {
                continue;
            };
            let Some(form) = touch_arms(a, b, touch.on_a, touch.on_b, tol) else {
                continue;
            };
            if !crosses(&touch, form, near, tol) || !apart(form.0, reach, near) {
                continue;
            }
            touch.arms = Some(form.0);
            touch.blur = blur_of(form.3, tol);
            if !touches.iter().any(|t| same_touch(t, &touch, near)) {
                touches.push(touch);
            }
            continue;
        }
        // The patch whose cut says more is halved: one across which the
        // normals turn most, or else the longest.
        let says = |c: Option<(usize, (bool, f64))>| c.map(|(_, w)| w);
        if says(cut_a) >= says(cut_b) {
            for q in halve(&mut trees[0], a, ia) {
                pairs.push((q, ib));
            }
        } else {
            for q in halve(&mut trees[1], b, ib) {
                pairs.push((ia, q));
            }
        }
    }
    Some(touches)
}

/// A parameter pair moved by whole periods to stand nearest `near`, so a
/// point joined onto a run of samples continues it in the chart.
pub(super) fn beside(surface: &SurfaceGeometry, at: (f64, f64), near: (f64, f64)) -> (f64, f64) {
    let ((ua, ub), (va, vb)) = surface.domain();
    let shift = |x: f64, to: f64, span: f64, periodic: bool| {
        if !periodic || !span.is_finite() || span <= 0.0 {
            return x;
        }
        x + ((to - x) / span).round() * span
    };
    (
        shift(at.0, near.0, ub - ua, surface.is_periodic_u()),
        shift(at.1, near.1, vb - va, surface.is_periodic_v()),
    )
}
