//! A marched round on a straight crease that leans into the face it ends
//! on (ADR 0001, amendment "a marched band leaning into its end face"):
//! carried past the vertex and clipped by that face, the band ends on it,
//! the face keeps all it holds, and the volume is the rolling ball's.
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{FilletCorner, fillet_edges, fillet_edges_corner};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry, Transformable as _};
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The half frustum's cone: axis through (10, 10), radius 5 at z = 4
/// falling `SLOPE` per unit of height to 3 at z = 10.
const R0: f64 = 5.0;
const SLOPE: f64 = 1.0 / 3.0;
/// The round's radius.
const R: f64 = 1.0;

fn frame_at(p: Point, z: Direction) -> Frame {
    Frame::new(p, z, Direction::X, T).unwrap()
}

fn only_solid(m: &Model, s: &Shape) -> Shape {
    let solids = explore_unique(m, s, ShapeType::Solid).unwrap();
    assert_eq!(solids.len(), 1);
    solids[0].clone()
}

/// A plate 20 × 20 × 4 with a frustum (radii 5 to 3, 6 high) on it, cut
/// through its axis and kept on the side y ≥ 10: its uprights lie on the
/// cone's rulings and lean into the plate.
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
        R0,
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

fn uprights(m: &Model, s: &Shape) -> Vec<Shape> {
    [Point::new(6.0, 10.0, 7.0), Point::new(14.0, 10.0, 7.0)]
        .iter()
        .map(|&p| edge_near(m, s, p))
        .collect()
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// The rolling ball's reference, independent of the kernel: the half
/// frustum run on below the plate is convex, so rounding its two uprights
/// is its opening, the balls of radius `R` inside it. Those balls' centres
/// fill `E` = {y ≥ 10 + R, ρ ≤ R(z) − R / cos α}, and a point is kept
/// where it lies within `R` of `E`. The plate and the frustum's top trim
/// what is kept to 4 ≤ z ≤ 10.
mod reference {
    use super::{R, R0, SLOPE};

    fn cos_a() -> f64 {
        1.0 / SLOPE.hypot(1.0)
    }

    fn sin_a() -> f64 {
        SLOPE * cos_a()
    }

    /// `E`'s radius at height `z`.
    fn eroded(z: f64) -> f64 {
        R0 - (z - 4.0) * SLOPE - R / cos_a()
    }

    fn in_e(p: [f64; 3]) -> bool {
        p[1] >= 10.0 + R && (p[0] - 10.0).hypot(p[1] - 10.0) <= eroded(p[2])
    }

    /// The distance from `p` to `E`: to its plane face, its cone face, or
    /// the two curves where they meet, whichever holds the nearest point.
    fn to_e(p: [f64; 3]) -> f64 {
        if in_e(p) {
            return 0.0;
        }
        let mut best = f64::INFINITY;
        let y0 = 10.0 + R;
        if p[1] < y0 && in_e([p[0], y0 + 1e-12, p[2]]) {
            best = best.min(y0 - p[1]);
        }
        let (dx, dy) = (p[0] - 10.0, p[1] - 10.0);
        let rho = dx.hypot(dy);
        let d = (rho - eroded(p[2])) * cos_a();
        if d > 0.0 && rho > 0.0 {
            let (r, z) = (rho - d * cos_a(), p[2] - d * sin_a());
            if 10.0 + dy / rho * r >= y0 - 1e-12 && z.is_finite() {
                best = best.min(d);
            }
        }
        for side in [-1.0, 1.0] {
            let far = |z: f64| -> f64 {
                let r = eroded(z);
                if r < R {
                    return f64::INFINITY;
                }
                let x = 10.0 + side * (r * r - R * R).sqrt();
                ((p[0] - x).powi(2) + (p[1] - y0).powi(2) + (p[2] - z).powi(2)).sqrt()
            };
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
            for _ in 0..60 {
                let (c, d) = (b - g * (b - a), a + g * (b - a));
                if far(c) < far(d) {
                    b = d;
                } else {
                    a = c;
                }
            }
            best = best.min(far(0.5 * (a + b)));
        }
        best
    }

    /// Whether `p`, at 4 ≤ z ≤ 10, is kept of the half frustum.
    fn kept(p: [f64; 3]) -> bool {
        to_e(p) <= R
    }

    /// The area kept of the section at height `z`, read in polar form
    /// about (10, 12), which every section holds, over `rays` rays.
    fn section(z: f64, rays: u32) -> f64 {
        let mut sum = 0.0;
        for i in 0..rays {
            let phi = core::f64::consts::TAU * f64::from(i) / f64::from(rays);
            let (c, s) = (phi.cos(), phi.sin());
            let at = |t: f64| kept([10.0 + c * t, 12.0 + s * t, z]);
            let mut t = 8.0;
            while !at(t) {
                t -= 0.04;
            }
            let (mut lo, mut hi) = (t, t + 0.04);
            for _ in 0..36 {
                let m = 0.5 * (lo + hi);
                if at(m) {
                    lo = m;
                } else {
                    hi = m;
                }
            }
            sum += 0.5 * lo.midpoint(hi).powi(2);
        }
        sum * core::f64::consts::TAU / f64::from(rays)
    }

    /// What rounding the uprights takes off: the half frustum's section
    /// less the kept section, integrated over 4 ≤ z ≤ 10 by six-point
    /// Gauss-Legendre on `panels` panels.
    pub fn taken(rays: u32, panels: u32) -> f64 {
        const GL: [(f64, f64); 6] = [
            (-0.932_469_514_203_152, 0.171_324_492_379_170),
            (-0.661_209_386_466_265, 0.360_761_573_048_139),
            (-0.238_619_186_083_197, 0.467_913_934_572_691),
            (0.238_619_186_083_197, 0.467_913_934_572_691),
            (0.661_209_386_466_265, 0.360_761_573_048_139),
            (0.932_469_514_203_152, 0.171_324_492_379_170),
        ];
        let h = 6.0 / f64::from(panels);
        let mut sum = 0.0;
        for k in 0..panels {
            let mid = 4.0 + h * (f64::from(k) + 0.5);
            for (x, w) in GL {
                let z = mid + 0.5 * h * x;
                let r = R0 - (z - 4.0) * SLOPE;
                let half = 0.5 * core::f64::consts::PI * r * r;
                sum += w * 0.5 * h * (half - section(z, rays));
            }
        }
        sum
    }
}

/// What rounding the uprights takes off, by [`reference::taken`] at 4096
/// rays and 8 panels; 2048 rays give the same within 3e-8.
const TAKEN: f64 = 4.393_196_94;

/// [`TAKEN`] is the reference's, converged: at 2048 rays and 8 panels it
/// comes within 1e-7 of it, and at 1024 rays and 4 within 1e-6.
#[test]
#[ignore = "heavy: the reference at full resolution takes most of a minute"]
fn the_reference_gives_what_the_tests_read() {
    let fine = reference::taken(2048, 8);
    assert!((fine - TAKEN).abs() <= 1e-7, "{fine} against {TAKEN}");
    let coarse = reference::taken(1024, 4);
    assert!((coarse - TAKEN).abs() <= 1e-6, "{coarse} against {TAKEN}");
}

/// The uprights rounded 1: each band runs past the plate's top and is
/// clipped by it, past the top face and out of the solid. Sound, the
/// volume taken off the rolling ball's within 1e-5 of it (measured 2e-6), the plate whole,
/// each band G1 to its faces, the same in either order and under
/// `Rolling`, and only the bands named as the rounds' faces.
#[test]
fn uprights_leaning_into_the_plate_round_to_the_rolling_ball() {
    let mut m = Model::new();
    let s = half_frustum(&mut m);
    let before = volume(&m, &s);
    let edges = uprights(&m, &s);
    let want = TAKEN;
    let mut first: Option<f64> = None;
    for order in [edges.clone(), edges.iter().rev().cloned().collect()] {
        let mut rounded = m.clone();
        let built = fillet_edges(&mut rounded, &s, &order, R, T).unwrap();
        let found = ogeom_algo::check(&rounded, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        let taken = before - volume(&rounded, &built.shape);
        let error = (taken - want).abs() / want;
        eprintln!("taken {taken}, reference {want}, relative error {error:.2e}");
        assert!(error <= 1e-5, "taken {taken} against {want}");
        match first {
            None => first = Some(taken),
            Some(v) => assert!((taken - v).abs() <= 1e-9 * v, "{taken} against {v}"),
        }
        // The plate keeps all it holds below its top, beside each foot.
        for q in [
            Point::new(5.2, 10.5, 3.99),
            Point::new(14.8, 10.5, 3.99),
            Point::new(5.5, 10.2, 3.999),
        ] {
            assert_eq!(
                ogeom_algo::classify_in_solid_exact(&rounded, &built.shape, q, T).unwrap(),
                ogeom_algo::Containment::In,
                "{q:?}"
            );
        }
        for edge in &edges {
            let made = built.history.generated(edge);
            assert!(!made.is_empty());
            for face in made {
                let surface = face_surface(&rounded, face);
                assert!(
                    matches!(surface, SurfaceGeometry::BSpline(_)),
                    "{surface:?}"
                );
                // The march holds its fit in distance (the bands' fit
                // target, 2e-4), its rails exactly on the faces; their
                // normals agree to about 5e-4 rad.
                let turn = widest_turn_at_rails(&rounded, &built.shape, face);
                eprintln!("widest turn at a rail {turn:.2e} rad");
                assert!(turn <= 1e-3, "a crease of {turn} rad");
            }
        }
    }
    let mut rolled = m.clone();
    let built = fillet_edges_corner(&mut rolled, &s, &edges, R, FilletCorner::Rolling, T).unwrap();
    let taken = before - volume(&rolled, &built.shape);
    assert!((taken - first.unwrap()).abs() <= 1e-9 * taken);
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

/// The widest turn between the band's normal and the normal of the face
/// beside it, read at 39 points along each of its two rails.
fn widest_turn_at_rails(m: &Model, s: &Shape, band: &Shape) -> f64 {
    let faces = explore_unique(m, s, ShapeType::Face).unwrap();
    let normal = |f: &Shape, p: Point| {
        let surface = face_surface(m, f);
        let (u, v) = ogeom_algo::project_on_surface(&surface, p, 16, T)
            .unwrap()
            .parameters;
        surface.normal_at(u, v, T).unwrap().vector()
    };
    let mut rails = 0;
    let mut worst = 0.0_f64;
    for edge in explore_unique(m, band, ShapeType::Edge).unwrap() {
        let Some(other) = faces.iter().find(|f| {
            !f.is_same(band)
                && explore_unique(m, f, ShapeType::Edge)
                    .unwrap()
                    .iter()
                    .any(|e| e.is_same(&edge))
        }) else {
            continue;
        };
        let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, &edge, T).unwrap().unwrap();
        let length = c
            .point_at(a, T)
            .unwrap()
            .distance(c.point_at(b, T).unwrap());
        // The rails run the band's length; its ends are short arcs on the
        // plate and the top, where it meets them at an angle.
        if length < 3.0 {
            continue;
        }
        rails += 1;
        for i in 1..40 {
            let p = c.point_at(a + (b - a) * f64::from(i) / 40.0, T).unwrap();
            let turn = normal(band, p)
                .dot(normal(other, p))
                .abs()
                .clamp(0.0, 1.0)
                .acos();
            worst = worst.max(turn);
        }
    }
    assert_eq!(rails, 2);
    worst
}
