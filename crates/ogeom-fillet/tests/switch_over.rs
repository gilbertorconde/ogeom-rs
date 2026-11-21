//! Rounds whose ball changes walls off the chain's vertex (ADR 0001,
//! amendment "a ball changing walls off the vertex").
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

/// The free-ball reference for a concave foot loop round a convex wall,
/// independent of the kernel.
///
/// The solid is a plate (z ≤ `base`) under the opening `S = (E ⊕ B_r)` of
/// a convex body by balls of radius `r`, where `E` is the set of the
/// opening balls' centres. Rounding the foot loop at `r` fills what no
/// ball of radius `r` clear of the solid reaches: the free centres are
/// `F = {z ≥ base + r} ∖ int D`, `D = E ⊕ B_2r`. For a point `p` outside
/// `S` with `base ≤ z < base + r`, `dist(p, F) ≤ r` exactly where one of
/// three holds:
/// - its lift to `z = base + r` lies outside `D`;
/// - `p` moved `2r − d` along the normal from its nearest point of `E`
///   (`d` its distance to `E`, the nearest point of `∂D`) stands at or
///   above `z = base + r`;
/// - `p` lies within `r` of the curve `Γ = ∂D ∩ {z = base + r}`.
///
/// Every other candidate nearest point of `F` is either on the plane `z =
/// base + r` (so at the lift or on `Γ`) or a critical point on `∂D`
/// above it, which for `E` convex is the second case or lies farther than
/// `2r`.
mod reference {
    /// A convex body `E`: its nearest point to `p` (and the distance),
    /// `p` itself when inside.
    pub trait Body {
        fn nearest(&self, p: [f64; 3]) -> (f64, [f64; 3]);
    }

    pub struct Fill<'a, B: Body> {
        pub body: &'a B,
        pub base: f64,
        pub r: f64,
        /// A point of `E` at every height the loop's round spans, about
        /// which each section of `S` and `D` is star-shaped.
        pub centre: [f64; 2],
        gamma: Vec<[f64; 2]>,
    }

    fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
        ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
    }

    fn seg_dist(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
        let (ux, uy) = (b[0] - a[0], b[1] - a[1]);
        let (vx, vy) = (p[0] - a[0], p[1] - a[1]);
        let l2 = ux * ux + uy * uy;
        let t = if l2 > 0.0 {
            ((ux * vx + uy * vy) / l2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        (vx - t * ux).hypot(vy - t * uy)
    }

    /// The first `t` in `[lo, hi]` where `inside` turns false, `inside(lo)`
    /// true and `inside(hi)` false.
    fn bisect(mut lo: f64, mut hi: f64, inside: impl Fn(f64) -> bool) -> f64 {
        for _ in 0..48 {
            let m = 0.5 * (lo + hi);
            if inside(m) {
                lo = m;
            } else {
                hi = m;
            }
        }
        0.5 * (lo + hi)
    }

    impl<'a, B: Body> Fill<'a, B> {
        pub fn new(body: &'a B, base: f64, r: f64, centre: [f64; 2], samples: u32) -> Self {
            let z = base + r;
            let gamma = (0..samples)
                .map(|i| {
                    let phi = core::f64::consts::TAU * f64::from(i) / f64::from(samples);
                    let (c, s) = (phi.cos(), phi.sin());
                    let at = |t: f64| [centre[0] + c * t, centre[1] + s * t, z];
                    let mut hi = 1.0;
                    while body.nearest(at(hi)).0 < 2.0 * r {
                        hi *= 2.0;
                    }
                    let t = bisect(0.0, hi, |t| body.nearest(at(t)).0 < 2.0 * r);
                    [centre[0] + c * t, centre[1] + s * t]
                })
                .collect();
            Self {
                body,
                base,
                r,
                centre,
                gamma,
            }
        }

        fn to_gamma(&self, p: [f64; 2]) -> f64 {
            let n = self.gamma.len();
            let stride = 16;
            let mut best = (f64::INFINITY, 0);
            for i in (0..n).step_by(stride) {
                let g = self.gamma[i];
                let d = (g[0] - p[0]).hypot(g[1] - p[1]);
                if d < best.0 {
                    best = (d, i);
                }
            }
            let mut least = f64::INFINITY;
            for k in 0..4 * stride {
                let i = (best.1 + n + k - 2 * stride) % n;
                least = least.min(seg_dist(p, self.gamma[i], self.gamma[(i + 1) % n]));
            }
            least
        }

        fn in_solid(&self, p: [f64; 3]) -> bool {
            p[2] <= self.base || self.body.nearest(p).0 <= self.r
        }

        /// Whether `p`, outside the solid with `base ≤ z < base + r`, is
        /// filled.
        fn filled(&self, p: [f64; 3]) -> bool {
            let (r, top) = (self.r, self.base + self.r);
            if self.body.nearest([p[0], p[1], top]).0 >= 2.0 * r {
                return false;
            }
            let (d, e) = self.body.nearest(p);
            if d > 0.0 && p[2] + (2.0 * r - d) * (p[2] - e[2]) / d >= top {
                return false;
            }
            let h = top - p[2];
            self.to_gamma([p[0], p[1]]) > (r * r - h * h).max(0.0).sqrt()
        }

        /// The fill's area in the section at height `z` between `phi`
        /// and `phi + dphi`, read on the ray at `phi` (per unit angle).
        fn ray(&self, phi: f64, z: f64) -> f64 {
            let (c, s) = (phi.cos(), phi.sin());
            let at = |t: f64| [self.centre[0] + c * t, self.centre[1] + s * t, z];
            let mut t = 0.0;
            while self.in_solid(at(t)) {
                t += 0.05;
            }
            let wall = bisect(t - 0.05, t, |t| self.in_solid(at(t)));
            let eps = 1e-9;
            if !self.filled(at(wall + eps)) {
                return 0.0;
            }
            let end = bisect(wall + eps, wall + 3.0 * self.r, |t| self.filled(at(t)));
            0.5 * (end * end - wall * wall)
        }

        /// The volume the loop's round adds: `rays` rays by `panels`
        /// six-point Gauss panels in `s`, `z = base + r s²`, on each
        /// stretch between the heights in `breaks` (where the wall's
        /// contact stands, and the integrand kinks).
        pub fn volume(&self, rays: u32, panels: u32, breaks: &[f64]) -> f64 {
            const GL: [(f64, f64); 6] = [
                (-0.932_469_514_203_152, 0.171_324_492_379_170),
                (-0.661_209_386_466_265, 0.360_761_573_048_139),
                (-0.238_619_186_083_197, 0.467_913_934_572_691),
                (0.238_619_186_083_197, 0.467_913_934_572_691),
                (0.661_209_386_466_265, 0.360_761_573_048_139),
                (0.932_469_514_203_152, 0.171_324_492_379_170),
            ];
            let mut cuts = vec![0.0];
            cuts.extend(breaks.iter().map(|z| ((z - self.base) / self.r).sqrt()));
            cuts.push(1.0);
            let dphi = core::f64::consts::TAU / f64::from(rays);
            let mut sum = 0.0;
            for w in cuts.windows(2) {
                let h = (w[1] - w[0]) / f64::from(panels);
                for k in 0..panels {
                    let mid = w[0] + h * (f64::from(k) + 0.5);
                    for (x, g) in GL {
                        let s = mid + 0.5 * h * x;
                        let z = self.base + self.r * s * s;
                        let jac = 2.0 * self.r * s * g * 0.5 * h;
                        let mut ring = 0.0;
                        for i in 0..rays {
                            ring += self.ray(f64::from(i) * dphi, z);
                        }
                        sum += jac * ring * dphi;
                    }
                }
            }
            sum
        }
    }

    /// The half frustum's eroded body: `{y ≥ 10 + r, ρ ≤ R(z) − r / cos α}`
    /// about the axis through (10, 10), `R(z) = r0 − (z − 4) slope`.
    pub struct Frustum {
        pub r0: f64,
        pub slope: f64,
        pub r: f64,
    }

    impl Frustum {
        fn cos_a(&self) -> f64 {
            1.0 / self.slope.hypot(1.0)
        }

        fn eroded(&self, z: f64) -> f64 {
            self.r0 - (z - 4.0) * self.slope - self.r / self.cos_a()
        }

        fn inside(&self, p: [f64; 3]) -> bool {
            p[1] >= 10.0 + self.r - 1e-12
                && (p[0] - 10.0).hypot(p[1] - 10.0) <= self.eroded(p[2]) + 1e-12
        }
    }

    impl Body for Frustum {
        fn nearest(&self, p: [f64; 3]) -> (f64, [f64; 3]) {
            if self.inside(p) {
                return (0.0, p);
            }
            let (ca, sa) = (self.cos_a(), self.slope * self.cos_a());
            let y0 = 10.0 + self.r;
            let mut best = (f64::INFINITY, p);
            if p[1] < y0 && self.inside([p[0], y0, p[2]]) {
                best = (y0 - p[1], [p[0], y0, p[2]]);
            }
            let (dx, dy) = (p[0] - 10.0, p[1] - 10.0);
            let rho = dx.hypot(dy);
            let d = (rho - self.eroded(p[2])) * ca;
            if d > 0.0 && rho > 0.0 {
                let (rf, zf) = (rho - d * ca, p[2] - d * sa);
                if rf >= 0.0 && 10.0 + dy / rho * rf >= y0 - 1e-12 && d < best.0 {
                    best = (d, [10.0 + dx / rho * rf, 10.0 + dy / rho * rf, zf]);
                }
            }
            for side in [-1.0, 1.0] {
                let on = |z: f64| -> Option<[f64; 3]> {
                    let e = self.eroded(z);
                    (e >= self.r).then(|| [10.0 + side * (e * e - self.r * self.r).sqrt(), y0, z])
                };
                let far = |z: f64| on(z).map_or(f64::INFINITY, |q| dist(p, q));
                let (mut a, mut b) = (p[2] - 4.0, p[2] + 4.0);
                let (mut at, mut least) = (a, f64::INFINITY);
                for i in 0..=64 {
                    let z = a + (b - a) * f64::from(i) / 64.0;
                    if far(z) < least {
                        (at, least) = (z, far(z));
                    }
                }
                let step = (b - a) / 64.0;
                (a, b) = (at - step, at + step);
                let g = (5.0_f64.sqrt() - 1.0) / 2.0;
                for _ in 0..80 {
                    let (c, d) = (b - g * (b - a), a + g * (b - a));
                    if far(c) < far(d) {
                        b = d;
                    } else {
                        a = c;
                    }
                }
                let z = 0.5 * (a + b);
                if let Some(q) = on(z)
                    && dist(p, q) < best.0
                {
                    best = (dist(p, q), q);
                }
            }
            best
        }
    }
}

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_fillet::{FilletCorner, fillet_edges, fillet_edges_corner};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry, Transformable as _};
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The D of radius `big` on its chord with its two corners rounded `r`:
/// its perimeter.
fn rounded_d_perimeter(big: f64, r: f64) -> f64 {
    let x = ((big - r).powi(2) - r * r).sqrt();
    let phi = (r / x).atan();
    2.0 * x + 2.0 * r * (PI / 2.0 + phi) + big * (PI - 2.0 * phi)
}

/// The reference on an upright D (no slope) gives the closed form: the
/// round of radius `r` round a convex outline of perimeter `p` adds
/// `p r²(1 − π/4) + π r³(5/3 − π/2)`.
#[test]
fn the_reference_rounds_an_upright_d_as_the_closed_form() {
    let body = reference::Frustum {
        r0: 5.0,
        slope: 0.0,
        r: 1.0,
    };
    let fill = reference::Fill::new(&body, 4.0, 1.0, [10.0, 12.0], 16384);
    let got = fill.volume(512, 2, &[]);
    let p = rounded_d_perimeter(5.0, 1.0);
    let want = p * (1.0 - PI / 4.0) + PI * (5.0 / 3.0 - PI / 2.0);
    eprintln!(
        "got {got}, want {want}, error {:.2e}",
        (got - want).abs() / want
    );
    assert!((got - want).abs() <= 1e-5 * want, "{got} against {want}");
}

/// The half frustum's foot loop rounded 1 once its uprights are: the
/// reference's fill at 2048 rays by 4 panels a stretch; 1024 rays agree
/// within 3e-8, and 8 panels within 3e-9.
const FRUSTUM_FILL: f64 = 3.438_297_00;

/// Where the ball on the cone's foot touches the cone: the integrand's
/// kink in height.
fn frustum_contact() -> f64 {
    5.0 - (1.0 / 3.0) / (1.0_f64 / 9.0 + 1.0).sqrt()
}

fn frustum_body() -> reference::Frustum {
    reference::Frustum {
        r0: 5.0,
        slope: 1.0 / 3.0,
        r: 1.0,
    }
}

/// [`FRUSTUM_FILL`] is the reference's, converged.
#[test]
#[ignore = "heavy: the reference at full resolution takes over two minutes"]
fn the_reference_gives_the_frustum_fill_the_tests_read() {
    let body = frustum_body();
    let fill = reference::Fill::new(&body, 4.0, 1.0, [10.0, 12.0], 16384);
    let fine = fill.volume(2048, 4, &[frustum_contact()]);
    assert!(
        (fine - FRUSTUM_FILL).abs() <= 1e-8,
        "{fine} against {FRUSTUM_FILL}"
    );
    let coarse = fill.volume(1024, 4, &[frustum_contact()]);
    assert!(
        (coarse - FRUSTUM_FILL).abs() <= 1e-7,
        "{coarse} against {FRUSTUM_FILL}"
    );
}

fn frame_at(p: Point, z: Direction) -> Frame {
    let x = if z == Direction::Z {
        Direction::X
    } else {
        Direction::Z
    };
    Frame::new(p, z, x, T).unwrap()
}

fn only_solid(m: &Model, s: &Shape) -> Shape {
    let solids = explore_unique(m, s, ShapeType::Solid).unwrap();
    assert_eq!(solids.len(), 1);
    solids[0].clone()
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// The edge of `s` passing nearest `q`, 32 samples along each.
fn edge_near(m: &Model, s: &Shape, q: Point) -> Shape {
    let d = |e: &Shape| {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
            return f64::INFINITY;
        };
        (0..=32)
            .map(|i| {
                c.point_at(a + (b - a) * f64::from(i) / 32.0, T)
                    .unwrap()
                    .distance(q)
            })
            .fold(f64::INFINITY, f64::min)
    };
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| d(a).total_cmp(&d(b)))
        .unwrap()
}

/// A plate 20 × 20 × 4 with a frustum (radii 5 to 3, 6 high) on it, cut
/// through its axis and kept on the side y ≥ 10.
fn half_frustum(m: &mut Model) -> Shape {
    let plate = ogeom_algo::make_box(
        m,
        frame_at(Point::new(0.0, 0.0, 0.0), Direction::Z),
        (20.0, 20.0, 4.0),
        T,
    )
    .unwrap()
    .shape;
    let cone = ogeom_algo::make_cone(
        m,
        frame_at(Point::new(10.0, 10.0, 4.0), Direction::Z),
        5.0,
        3.0,
        6.0,
        T,
    )
    .unwrap()
    .shape;
    let cutter = ogeom_algo::make_box(
        m,
        frame_at(Point::new(4.0, 4.0, 3.0), Direction::Z),
        (12.0, 6.0, 8.0),
        T,
    )
    .unwrap()
    .shape;
    let half = ogeom_bool::cut(m, &cone, &cutter, T).unwrap().shape;
    let half = only_solid(m, &half);
    let fused = ogeom_bool::fuse(m, &plate, &half, T).unwrap().shape;
    only_solid(m, &fused)
}

fn frustum_uprights(m: &Model, s: &Shape) -> Vec<Shape> {
    [Point::new(6.0, 10.0, 7.0), Point::new(14.0, 10.0, 7.0)]
        .iter()
        .map(|&p| edge_near(m, s, p))
        .collect()
}

/// The edges of `s` lying on the plane `z = top` away from its outline
/// `[0, size]²`: a boss's foot loop.
fn foot_loop(m: &Model, s: &Shape, top: f64, size: f64) -> Vec<Shape> {
    let mut foot = Vec::new();
    for e in explore_unique(m, s, ShapeType::Edge).unwrap() {
        let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap() else {
            continue;
        };
        let on = (0..=8).all(|i| {
            let p = c.point_at(a + (b - a) * f64::from(i) / 8.0, T).unwrap();
            (p.z - top).abs() < 1e-6
                && p.x > 1.0
                && p.x < size - 1.0
                && p.y > 1.0
                && p.y < size - 1.0
        });
        if on {
            foot.push(e);
        }
    }
    foot
}

fn face_surface(m: &Model, f: &Shape) -> SurfaceGeometry {
    m.node(f)
        .and_then(|n| n.data().as_face())
        .and_then(|d| m.geometry().surface(d.surface))
        .unwrap()
        .clone()
        .transformed(&f.transform(m.datums()).unwrap(), T)
        .unwrap()
}

/// The widest turn between the normals of `faces` and those of the faces
/// beside them across each of their edges, read at 15 points along each.
fn widest_turn(m: &Model, s: &Shape, faces: &[Shape]) -> f64 {
    let all = explore_unique(m, s, ShapeType::Face).unwrap();
    let normal = |f: &Shape, p: Point| {
        let surface = face_surface(m, f);
        let (u, v) = ogeom_algo::project_on_surface(&surface, p, 16, T)
            .unwrap()
            .parameters;
        let n = surface.normal_at(u, v, T).unwrap().vector();
        if f.orientation() == ogeom_topo::Orientation::Reversed {
            -n
        } else {
            n
        }
    };
    let mut worst = 0.0_f64;
    for face in faces {
        for edge in explore_unique(m, face, ShapeType::Edge).unwrap() {
            let Some(other) = all.iter().find(|f| {
                !f.is_same(face)
                    && explore_unique(m, f, ShapeType::Edge)
                        .unwrap()
                        .iter()
                        .any(|e| e.is_same(&edge))
            }) else {
                continue;
            };
            let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, &edge, T).unwrap() else {
                continue;
            };
            for i in 1..16 {
                let p = c.point_at(a + (b - a) * f64::from(i) / 16.0, T).unwrap();
                let turn = normal(face, p)
                    .dot(normal(other, p))
                    .clamp(-1.0, 1.0)
                    .acos();
                worst = worst.max(turn);
            }
        }
    }
    worst
}

/// The half frustum with its uprights rounded 1 (each band leaning into
/// the plate, carried past it and clipped there), then its foot loop
/// rounded 1: the plate's line, the two bands' feet and the cone's arc.
/// Where the ball along the line or the arc reaches a band, its contact
/// crosses the band's rail at its own height, off the vertex where the
/// creases meet (0.35 from it on the line, 0.09 on the arc): between
/// those crossings it touches the plate and the band alone, and the round
/// there is the pipe about its centre's curve. The fill is the reference's
/// within 1e-5, sound, the same in either order, G1 across every join.
#[test]
fn a_half_frustum_s_foot_rounds_through_its_pipe_corners() {
    let mut m = Model::new();
    let s = half_frustum(&mut m);
    let uprights = frustum_uprights(&m, &s);
    let up = fillet_edges(&mut m, &s, &uprights, 1.0, T).unwrap();
    let before = volume(&m, &up.shape);
    let foot = foot_loop(&m, &up.shape, 4.0, 20.0);
    assert_eq!(foot.len(), 4);
    let mut first: Option<f64> = None;
    for order in [foot.clone(), foot.iter().rev().cloned().collect()] {
        let mut rounded = m.clone();
        let built = fillet_edges(&mut rounded, &up.shape, &order, 1.0, T).unwrap();
        let found = ogeom_algo::check(&rounded, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        let fill = volume(&rounded, &built.shape) - before;
        let error = (fill - FRUSTUM_FILL).abs() / FRUSTUM_FILL;
        eprintln!("fill {fill}, reference {FRUSTUM_FILL}, relative error {error:.2e}");
        assert!(error <= 1e-5, "fill {fill} against {FRUSTUM_FILL}");
        match first {
            None => first = Some(fill),
            Some(v) => assert!((fill - v).abs() <= 1e-9 * v, "{fill} against {v}"),
        }
        let made: Vec<Shape> = foot
            .iter()
            .flat_map(|e| built.history.generated(e).to_vec())
            .collect();
        assert!(!made.is_empty());
        let turn = widest_turn(&rounded, &built.shape, &made);
        eprintln!("widest turn across a join {turn:.2e} rad");
        assert!(turn <= 1e-3, "a crease of {turn} rad");
    }
    // In one request under `Rolling` the uprights go first and the foot
    // rolls round them: the same solid in either order.
    let after = volume(&m, &up.shape) + first.unwrap();
    let mut asked = vec![
        edge_near(&m, &s, Point::new(10.0, 10.0, 4.0)),
        edge_near(&m, &s, Point::new(10.0, 15.0, 4.0)),
    ];
    asked.extend(uprights);
    for order in [asked.clone(), asked.iter().rev().cloned().collect()] {
        let mut rolled = m.clone();
        let built =
            fillet_edges_corner(&mut rolled, &s, &order, 1.0, FilletCorner::Rolling, T).unwrap();
        let found = ogeom_algo::check(&rolled, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        let got = volume(&rolled, &built.shape);
        assert!((got - after).abs() <= 1e-9 * after, "{got} against {after}");
    }
}

/// `edges` of `s` rounded `r`: refused, naming a vertex where the ball
/// changes walls off its section by a distance reading `by`.
fn refused_off_vertex(m: &Model, s: &Shape, edges: &[Shape], r: f64, by: &str) {
    for order in [edges.to_vec(), edges.iter().rev().cloned().collect()] {
        let got = fillet_edges(&mut m.clone(), s, &order, r, T)
            .unwrap_err()
            .to_string();
        assert!(
            got.contains("the ball changes walls off the vertex at (")
                && got.contains(&format!("seam between the walls {by}")),
            "{got}"
        );
    }
}

/// The half frustum with its uprights rounded 1, then its top loop
/// rounded 0.5: convex, the ball within the material changes walls where
/// its contact crosses each band's rail, 0.19 past the vertex on the line
/// (0.08 on the arc). No pipe is built for a convex chain yet: refused by
/// name.
#[test]
fn a_half_frustum_s_top_loop_is_refused_by_name() {
    let mut m = Model::new();
    let s = half_frustum(&mut m);
    let uprights = frustum_uprights(&m, &s);
    let up = fillet_edges(&mut m, &s, &uprights, 1.0, T).unwrap();
    let top = foot_loop(&m, &up.shape, 10.0, 20.0);
    assert_eq!(top.len(), 4);
    refused_off_vertex(&m, &up.shape, &top, 0.5, "0.19");
}

/// A half cylinder of radius 4 lying on a plate, its axis on the plate's
/// top, its end arcs rounded 1, then its foot loop rounded 1. Each arc's
/// torus band meets the end disc along a circle upright at the plate, but
/// at the ball's height the circle stands 0.17 along the disc's foot from
/// the vertex: the ball changes walls there. Its other end hands over in
/// the vertex's section on the cylinder, and a pipe ending there is not
/// built yet: refused by name.
#[test]
fn a_half_cylinder_s_foot_is_refused_by_name() {
    let mut m = Model::new();
    let plate = ogeom_algo::make_box(
        &mut m,
        frame_at(Point::new(0.0, 0.0, 0.0), Direction::Z),
        (20.0, 20.0, 6.0),
        T,
    )
    .unwrap()
    .shape;
    let drum = ogeom_algo::make_cylinder(
        &mut m,
        frame_at(Point::new(5.0, 10.0, 6.0), Direction::X),
        4.0,
        10.0,
        T,
    )
    .unwrap()
    .shape;
    let fused = ogeom_bool::fuse(&mut m, &plate, &drum, T).unwrap().shape;
    let s = only_solid(&m, &fused);
    let arcs: Vec<Shape> = explore_unique(&m, &s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(&m, e, T).unwrap().unwrap();
            let p = c.point_at(0.5 * (a + b), T).unwrap();
            ((p.x - 5.0).abs() < 1e-9 || (p.x - 15.0).abs() < 1e-9) && p.z > 6.5
        })
        .collect();
    assert_eq!(arcs.len(), 4);
    let up = fillet_edges(&mut m, &s, &arcs, 1.0, T).unwrap();
    let foot = foot_loop(&m, &up.shape, 6.0, 20.0);
    assert_eq!(foot.len(), 8);
    refused_off_vertex(&m, &up.shape, &foot, 1.0, "0.17");
}
