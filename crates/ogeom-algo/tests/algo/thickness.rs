//! Wall thickness by rays and by medial balls, against closed forms: a
//! plate reads its thickness, a wedge reads square across it with balls
//! and along the normal with rays, a ribbed plate reads its rib and the
//! rolling ball under the rib's root, a tube reads its wall, and a noisy
//! scan of a plate reads within the noise; as meshes and as B-reps.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_algo::thickness::{Method, Thickness, thickness, thickness_of_shape};
use ogeom_core::Tolerances;
use ogeom_math::{Axis, Frame, Point, Vector};
use ogeom_mesh::Bvh;
use ogeom_topo::{Model, Triangulation};

const T: Tolerances = Tolerances::millimetres();

const RAY: Method = Method::Ray { cone: 0.0, rays: 1 };

fn options(method: Method) -> Thickness {
    Thickness {
        method,
        ..Thickness::default()
    }
}

/// A deterministic generator: xorshift.
struct Seeded(u64);

impl Seeded {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    fn unit(&mut self) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let unit = (self.next() >> 11) as f64 / (1_u64 << 53) as f64;
        unit
    }

    /// A normal deviate, by Box and Muller.
    fn gaussian(&mut self) -> f64 {
        let (a, b) = (1.0 - self.unit(), self.unit());
        (-2.0 * a.ln()).sqrt() * (std::f64::consts::TAU * b).cos()
    }
}

/// A grid of quads over `corner + i·du + j·dv`, `n × m` of them, wound so
/// the triangles' normals run along `du × dv`; the grid's positions
/// first, row by row in `i`.
fn grid(mesh: &mut Triangulation, corner: Point, du: Vector, dv: Vector, n: u32, m: u32) {
    let first = u32::try_from(mesh.positions.len()).unwrap();
    for i in 0..=n {
        for j in 0..=m {
            mesh.positions.push(
                corner + du * (f64::from(i) / f64::from(n)) + dv * (f64::from(j) / f64::from(m)),
            );
        }
    }
    let index = |i: u32, j: u32| first + i * (m + 1) + j;
    for i in 0..n {
        for j in 0..m {
            let (a, b, c, d) = (
                index(i, j),
                index(i + 1, j),
                index(i + 1, j + 1),
                index(i, j + 1),
            );
            mesh.triangles.extend([[a, b, c], [a, c, d]]);
        }
    }
}

/// A closed box from `low` spanning `size`, each face a grid of
/// `cells[axis]` along each axis it spans, wound outward.
fn box_mesh(low: Point, size: [f64; 3], cells: [u32; 3]) -> Triangulation {
    let mut mesh = Triangulation::new();
    let axis = |k: usize| {
        let mut v = [0.0; 3];
        v[k] = size[k];
        Vector::new(v[0], v[1], v[2])
    };
    for k in 0..3 {
        let (u, v) = ((k + 1) % 3, (k + 2) % 3);
        let (du, dv) = (axis(u), axis(v));
        // (u, v, k) is right-handed: du × dv runs along +k.
        grid(&mut mesh, low + axis(k), du, dv, cells[u], cells[v]);
        grid(&mut mesh, low, dv, du, cells[v], cells[u]);
    }
    mesh
}

/// The points of a mesh's positions from `first` on that pass `keep`.
fn vertices(mesh: &Triangulation, first: usize, keep: impl Fn(Point) -> bool) -> Vec<Point> {
    mesh.positions[first..]
        .iter()
        .copied()
        .filter(|&p| keep(p))
        .collect()
}

fn read(bvh: &Bvh, points: &[Point], normal: impl Fn(Point) -> Vector, method: Method) -> Vec<f64> {
    let normals: Vec<Vector> = points.iter().map(|&p| normal(p)).collect();
    thickness(bvh, points, &normals, &options(method))
        .unwrap()
        .into_iter()
        .map(|w| w.expect("a wall"))
        .collect()
}

/// A plate 3 mm thick reads 3 mm at every vertex inside its two faces with
/// both methods, to within rounding, however near the rim: the side walls
/// are edge-on to the normal and passed over. A cone of rays reads it too
/// where none of its rays reaches a side wall.
#[test]
fn a_plate_reads_its_thickness_at_every_interior_vertex() {
    let mesh = box_mesh(Point::ORIGIN, [40.0, 30.0, 3.0], [40, 30, 3]);
    let bvh = Bvh::build(&mesh).unwrap();
    let inside = |p: Point| p.x > 0.0 && p.x < 40.0 && p.y > 0.0 && p.y < 30.0;
    let top = vertices(&mesh, 0, |p| inside(p) && p.z == 3.0);
    let bottom = vertices(&mesh, 0, |p| inside(p) && p.z == 0.0);
    assert_eq!(
        top.len(),
        39 * 29,
        "the face's grid holds each interior vertex once"
    );
    // A cone's rays reach the side walls from within `3 tan(cone)` of them.
    let cone = 0.5_f64;
    let margin = 3.0 * cone.tan();
    let clear =
        |p: &Point| p.x > margin && p.x < 40.0 - margin && p.y > margin && p.y < 30.0 - margin;
    let cone = Method::Ray { cone, rays: 16 };
    for method in [RAY, Method::Sphere, cone] {
        for (points, up) in [(&top, 1.0), (&bottom, -1.0)] {
            let points: Vec<Point> = points
                .iter()
                .copied()
                .filter(|p| method != cone || clear(p))
                .collect();
            for w in read(&bvh, &points, |_| Vector::Z * up, method) {
                assert!((w - 3.0).abs() <= 1e-9, "{method:?} read {w}");
            }
        }
    }
}

/// The wedge between `z = 0` and `z = x tan α` for `x` in `[0, length]`,
/// `width` deep along `y`, closed by its back at `x = length` and two
/// ends; the floor and the slope gridded `cells` along `x`.
fn wedge_mesh(alpha: f64, length: f64, width: f64, cells: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    let rise = length * alpha.tan();
    let (x, y) = (Vector::X * length, Vector::Y * width);
    let slope = Vector::new(length, 0.0, rise);
    grid(&mut mesh, Point::ORIGIN, y, x, cells, cells);
    grid(&mut mesh, Point::ORIGIN, slope, y, cells, cells);
    grid(
        &mut mesh,
        Point::new(length, 0.0, 0.0),
        y,
        Vector::Z * rise,
        1,
        1,
    );
    let first = u32::try_from(mesh.positions.len()).unwrap();
    for y in [0.0, width] {
        mesh.positions.extend([
            Point::new(0.0, y, 0.0),
            Point::new(length, y, 0.0),
            Point::new(length, y, rise),
        ]);
    }
    mesh.triangles.extend([
        [first, first + 1, first + 2],
        [first + 3, first + 5, first + 4],
    ]);
    mesh
}

/// A wedge of angle α read from its floor at `s` from its edge: rays along
/// the normal read `s tan α`, balls read `2 s tan(α/2)`, the diameter of
/// the ball tangent to both faces, square across the wedge; and a cone of
/// rays wider than α reads the perpendicular to the slope, `s sin α`, as
/// nearly as its rays come to it.
#[test]
fn a_wedge_reads_square_across_with_balls() {
    for alpha in [0.2_f64, 30_f64.to_radians(), 60_f64.to_radians()] {
        let mesh = wedge_mesh(alpha, 10.0, 10.0, 40);
        let bvh = Bvh::build(&mesh).unwrap();
        let floor = vertices(&mesh, 0, |p| {
            p.z == 0.0 && p.x > 0.0 && p.x <= 7.0 && p.y > 2.0 && p.y < 8.0
        });
        assert!(floor.len() > 100);
        let down = |_| -Vector::Z;
        let rays = read(&bvh, &floor, down, RAY);
        let balls = read(&bvh, &floor, down, Method::Sphere);
        for ((p, ray), ball) in floor.iter().zip(&rays).zip(&balls) {
            let s = p.x;
            let along = s * alpha.tan();
            let across = 2.0 * s * (alpha / 2.0).tan();
            assert!(
                (ray - along).abs() <= 1e-9 * along,
                "ray {ray} for {along} at {s}"
            );
            assert!(
                (ball - across).abs() <= 1e-9 * across,
                "ball {ball} for {across} at {s}"
            );
        }
        let cone = Method::Ray {
            cone: alpha + 0.1,
            rays: 2000,
        };
        // Near enough the edge that no ray of the cone reaches the back or
        // the ends before the slope's perpendicular length.
        let middle: Vec<Point> = floor
            .iter()
            .copied()
            .filter(|p| p.y == 5.0 && p.x <= 3.0)
            .collect();
        for (p, w) in middle.iter().zip(read(&bvh, &middle, down, cone)) {
            let square = p.x * alpha.sin();
            assert!(
                w >= square * (1.0 - 1e-12) && w <= square * 1.001,
                "cone {w} for {square}"
            );
        }
    }
}

/// A plate 3 mm thick standing a rib 0.8 mm wide and 10 mm tall, extruded
/// 20 mm along `y`: its section's outline is gridded at 0.1 mm across and
/// 0.5 mm up, and the two ends close it.
fn ribbed_plate() -> Triangulation {
    let section = [
        (-10.0, 0.0),
        (10.0, 0.0),
        (10.0, 3.0),
        (0.4, 3.0),
        (0.4, 13.0),
        (-0.4, 13.0),
        (-0.4, 3.0),
        (-10.0, 3.0),
    ];
    let mut mesh = Triangulation::new();
    let along = Vector::Y * 20.0;
    for k in 0..section.len() {
        let (a, b) = (section[k], section[(k + 1) % section.len()]);
        let (from, to) = (Point::new(a.0, 0.0, a.1), Point::new(b.0, 0.0, b.1));
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let cells = ((to - from).magnitude() / if a.1 == b.1 { 0.1 } else { 0.5 }).round() as u32;
        grid(&mut mesh, from, to - from, along, cells.max(1), 40);
    }
    for y in [0.0, 20.0] {
        let first = u32::try_from(mesh.positions.len()).unwrap();
        for (x, z) in [
            (-10.0, 0.0),
            (10.0, 0.0),
            (10.0, 3.0),
            (-10.0, 3.0),
            (-0.4, 3.0),
            (0.4, 3.0),
            (0.4, 13.0),
            (-0.4, 13.0),
        ] {
            mesh.positions.push(Point::new(x, y, z));
        }
        let quads = [[0, 1, 2, 3], [4, 5, 6, 7]];
        for q in quads {
            let [a, b, c, d] = q.map(|k| first + k);
            if y == 0.0 {
                mesh.triangles.extend([[a, b, c], [a, c, d]]);
            } else {
                mesh.triangles.extend([[a, c, b], [a, d, c]]);
            }
        }
    }
    mesh
}

/// A rib reads its width from either wall with both methods, its root's
/// floor edge-on and passed over; the plate under the rib reads through
/// the rib with a ray, and with a ball the ball tangent at the point that
/// touches the rib's root corner nearer it, `((0.4 − |x|)² + 9) / 3`.
#[test]
fn a_thin_rib_reads_its_width_and_the_ball_under_it_rolls_on_its_root() {
    let mesh = ribbed_plate();
    let bvh = Bvh::build(&mesh).unwrap();
    let mid = |p: Point| p.y > 1.0 && p.y < 19.0;
    for side in [-0.4, 0.4] {
        let wall = vertices(&mesh, 0, |p| {
            mid(p) && p.x == side && p.z > 3.0 && p.z < 13.0
        });
        assert!(wall.len() > 500);
        for method in [RAY, Method::Sphere] {
            for w in read(&bvh, &wall, |_| Vector::X * side.signum(), method) {
                assert!((w - 0.8).abs() <= 1e-12, "{method:?} read {w} on the rib");
            }
        }
    }
    let under = vertices(&mesh, 0, |p| mid(p) && p.z == 0.0 && p.x.abs() < 0.35);
    assert!(under.len() > 100);
    let down = |_| -Vector::Z;
    for w in read(&bvh, &under, down, RAY) {
        assert!((w - 13.0).abs() <= 1e-12, "a ray under the rib read {w}");
    }
    for (p, w) in under.iter().zip(read(&bvh, &under, down, Method::Sphere)) {
        let rolled = ((0.4 - p.x.abs()).powi(2) + 9.0) / 3.0;
        assert!(
            (w - rolled).abs() <= 1e-9,
            "a ball under the rib read {w} for {rolled}"
        );
    }
}

/// A tube from radius 8 to 10 and 20 long, `around` facets about it, its
/// two walls' vertices at the same angles, with its two rings of ends.
fn tube(around: u32, rows: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    let at = |r: f64, i: u32, z: f64| {
        let a = std::f64::consts::TAU * f64::from(i % around) / f64::from(around);
        Point::new(r * a.cos(), r * a.sin(), z)
    };
    let wall = |mesh: &mut Triangulation, r: f64, outward: bool| {
        let first = u32::try_from(mesh.positions.len()).unwrap();
        for i in 0..around {
            for j in 0..=rows {
                mesh.positions
                    .push(at(r, i, 20.0 * f64::from(j) / f64::from(rows)));
            }
        }
        let index = |i: u32, j: u32| first + (i % around) * (rows + 1) + j;
        for i in 0..around {
            for j in 0..rows {
                let (a, b, c, d) = (
                    index(i, j),
                    index(i + 1, j),
                    index(i + 1, j + 1),
                    index(i, j + 1),
                );
                if outward {
                    mesh.triangles.extend([[a, b, c], [a, c, d]]);
                } else {
                    mesh.triangles.extend([[a, c, b], [a, d, c]]);
                }
            }
        }
        first
    };
    let outer = wall(&mut mesh, 10.0, true);
    let inner = wall(&mut mesh, 8.0, false);
    for (j, up) in [(0, false), (rows, true)] {
        for i in 0..around {
            let o = |i: u32| outer + (i % around) * (rows + 1) + j;
            let n = |i: u32| inner + (i % around) * (rows + 1) + j;
            let (a, b, c, d) = (n(i), o(i), o(i + 1), n(i + 1));
            if up {
                mesh.triangles.extend([[a, b, c], [a, c, d]]);
            } else {
                mesh.triangles.extend([[a, c, b], [a, d, c]]);
            }
        }
    }
    mesh
}

/// A faceted tube's wall reads 2 mm from outside and from inside with rays:
/// a ray along the radius meets the other wall's vertex at the same angle.
/// From outside the ball across touches the inner wall's vertex there,
/// clear of its facets, and reads 2 mm too. From inside it meets the outer
/// wall's vertex where the outer facets turn in round it, `δ = 2π/96`
/// apart, and reads the ball tangent to them, `4 cos(δ/2) / (1 + cos(δ/2))`,
/// as a ball fitting the faceted wall must.
#[test]
fn a_tube_reads_its_wall_from_both_sides() {
    let (around, rows) = (96, 20);
    let mesh = tube(around, rows);
    let bvh = Bvh::build(&mesh).unwrap();
    let rows_inside = |p: Point| p.z > 0.0 && p.z < 20.0;
    let count = (around * (rows + 1)) as usize;
    let outer = vertices(&mesh, 0, rows_inside);
    let outer = &outer[..outer.len() / 2];
    let inner: Vec<Point> = mesh.positions[count..2 * count]
        .iter()
        .copied()
        .filter(|&p| rows_inside(p))
        .collect();
    let radial = |p: Point| Vector::new(p.x, p.y, 0.0);
    let half = (std::f64::consts::PI / f64::from(around)).cos();
    let faceted = 4.0 * half / (1.0 + half);
    for (method, inside) in [(RAY, 2.0), (Method::Sphere, faceted)] {
        for w in read(&bvh, outer, radial, method) {
            assert!((w - 2.0).abs() <= 1e-9, "{method:?} read {w} from outside");
        }
        for w in read(&bvh, &inner, |p| -radial(p), method) {
            assert!(
                (w - inside).abs() <= 1e-9,
                "{method:?} read {w} from inside"
            );
        }
    }
}

/// A plate 3 mm thick scanned with a normal noise of 0.02 mm on each
/// vertex's height, read at each vertex along the scan's own normal there
/// (its triangles' area-weighted normals). A ray drops between two noisy
/// sheets, each within `a` of its plane, `a` the largest noise drawn, along
/// a normal turned `θ` from the plane's, so it reads within
/// `[3 − 2a, (3 + 2a) / cos θ]`; across all vertices its median is within
/// the noise of 3 mm. A ball is never longer than the ray, and its median
/// is within the noise too.
#[test]
fn a_noisy_scan_of_a_plate_reads_within_its_noise() {
    let sigma = 0.02;
    let n = 60;
    let mut mesh = Triangulation::new();
    grid(
        &mut mesh,
        Point::new(0.0, 0.0, 3.0),
        Vector::X * 30.0,
        Vector::Y * 30.0,
        n,
        n,
    );
    grid(
        &mut mesh,
        Point::ORIGIN,
        Vector::Y * 30.0,
        Vector::X * 30.0,
        n,
        n,
    );
    let mut rng = Seeded(0x5ca9_7e57);
    let mut largest: f64 = 0.0;
    for p in &mut mesh.positions {
        let e = sigma * rng.gaussian();
        largest = largest.max(e.abs());
        p.z += e;
    }
    let mut normals = vec![Vector::ZERO; mesh.positions.len()];
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|k| mesh.positions[k as usize]);
        let area = (b - a).cross(c - a);
        for &k in t {
            normals[k as usize] += area;
        }
    }
    let bvh = Bvh::build(&mesh).unwrap();
    let inside = |p: Point| p.x > 3.0 && p.x < 27.0 && p.y > 3.0 && p.y < 27.0;
    let chosen: Vec<usize> = (0..mesh.positions.len())
        .filter(|&k| inside(mesh.positions[k]))
        .collect();
    assert_eq!(chosen.len(), 2 * 47 * 47);
    let points: Vec<Point> = chosen.iter().map(|&k| mesh.positions[k]).collect();
    let normals: Vec<Vector> = chosen.iter().map(|&k| normals[k]).collect();
    let rays = thickness(&bvh, &points, &normals, &options(RAY)).unwrap();
    let balls = thickness(&bvh, &points, &normals, &options(Method::Sphere)).unwrap();
    let median = |mut v: Vec<f64>| {
        v.sort_by(f64::total_cmp);
        v[v.len() / 2]
    };
    let mut ray_walls = Vec::new();
    let mut ball_walls = Vec::new();
    for ((ray, ball), m) in rays.iter().zip(&balls).zip(&normals) {
        let (ray, ball) = (ray.unwrap(), ball.unwrap());
        let cos = m.z.abs() / m.magnitude();
        assert!(
            ray >= 3.0 - 2.0 * largest && ray <= (3.0 + 2.0 * largest) / cos,
            "ray {ray}"
        );
        assert!(
            ball <= ray * (1.0 + 1e-12),
            "ball {ball} longer than its ray {ray}"
        );
        ray_walls.push(ray);
        ball_walls.push(ball);
    }
    let (ray, ball) = (median(ray_walls), median(ball_walls));
    assert!((ray - 3.0).abs() <= sigma, "median ray {ray}");
    assert!((ball - 3.0).abs() <= sigma, "median ball {ball}");
}

/// The answers are in the points' order and the same on one thread as on
/// many; a point with no normal, or a ray leaving the mesh, answers
/// nothing, and so does a wall thicker than the thickest asked for.
#[test]
fn answers_keep_their_order_on_any_thread_count_and_say_where_there_is_no_wall() {
    let mesh = wedge_mesh(0.5, 10.0, 10.0, 30);
    let bvh = Bvh::build(&mesh).unwrap();
    let floor = vertices(&mesh, 0, |p| p.z == 0.0 && p.x > 0.0 && p.x < 10.0);
    let normals = vec![-Vector::Z; floor.len()];
    for method in [RAY, Method::Sphere, Method::Ray { cone: 0.3, rays: 7 }] {
        let o = options(method);
        ogeom_core::parallel::set_threads(1);
        let one = thickness(&bvh, &floor, &normals, &o).unwrap();
        ogeom_core::parallel::set_threads(8);
        let many = thickness(&bvh, &floor, &normals, &o).unwrap();
        ogeom_core::parallel::set_threads(0);
        assert_eq!(one, many);
    }
    let points = [
        Point::new(5.0, 5.0, 0.0),
        Point::new(5.0, 5.0, 0.0),
        Point::new(5.0, 5.0, 0.0),
    ];
    let normals = [-Vector::Z, Vector::ZERO, Vector::Z];
    let found = thickness(&bvh, &points, &normals, &options(RAY)).unwrap();
    assert!(found[0].is_some() && found[1].is_none() && found[2].is_none());
    let thin = Thickness {
        max: 1.0,
        ..options(Method::Sphere)
    };
    // At 5 from the edge the ray reads 5 tan 0.5 and the ball
    // 10 tan 0.25, both over 1.
    assert_eq!(
        thickness(&bvh, &points[..1], &normals[..1], &thin).unwrap(),
        vec![None]
    );
    assert!(thickness(&bvh, &points, &normals[..2], &options(RAY)).is_err());
    let wide = options(Method::Ray {
        cone: std::f64::consts::FRAC_PI_2,
        rays: 3,
    });
    assert!(thickness(&bvh, &points, &normals, &wide).is_err());
}

/// The same walls on B-reps, measured on their own surfaces: a plate, a
/// wedge along the normal and square across, and a tube revolved from a
/// rectangle, from both sides.
#[test]
fn b_reps_read_their_walls_on_their_own_surfaces() {
    let mut model = Model::new();
    let plate = ogeom_algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 3.0), T)
        .unwrap()
        .shape;
    let top: Vec<Point> = (0..50)
        .map(|k| {
            let k = f64::from(k);
            Point::new(0.5 + 0.78 * k, 29.5 - 0.58 * k, 3.0)
        })
        .collect();
    let up = vec![Vector::Z; top.len()];
    for method in [RAY, Method::Sphere] {
        for w in thickness_of_shape(&model, &plate, &top, &up, &options(method)).unwrap() {
            let w = w.unwrap();
            assert!((w - 3.0).abs() <= 1e-9, "{method:?} read {w} on the plate");
        }
    }

    // A ridge at x = 0 over the floor's edge at x = 10: the slope rises
    // from the floor's far edge at the angle atan(4 / 10).
    let wedge = ogeom_algo::make_wedge(&mut model, Frame::WORLD, (10.0, 10.0, 4.0), (0.0, 10.0), T)
        .unwrap()
        .shape;
    let alpha = (4.0_f64 / 10.0).atan();
    let floor: Vec<Point> = (1..=14)
        .map(|k| Point::new(10.0 - 0.5 * f64::from(k), 5.0, 0.0))
        .collect();
    let down = vec![-Vector::Z; floor.len()];
    let rays = thickness_of_shape(&model, &wedge, &floor, &down, &options(RAY)).unwrap();
    let balls =
        thickness_of_shape(&model, &wedge, &floor, &down, &options(Method::Sphere)).unwrap();
    for ((p, ray), ball) in floor.iter().zip(rays).zip(balls) {
        let s = 10.0 - p.x;
        let (ray, ball) = (ray.unwrap(), ball.unwrap());
        assert!((ray - s * alpha.tan()).abs() <= 1e-9, "ray {ray} at {s}");
        let across = 2.0 * s * (alpha / 2.0).tan();
        assert!(
            (ball - across).abs() <= 4.0 * T.confusion(),
            "ball {ball} for {across} at {s}"
        );
    }

    let corners =
        [(8.0, 0.0), (10.0, 0.0), (10.0, 20.0), (8.0, 20.0)].map(|(x, z)| Point::new(x, 0.0, z));
    let ring = ogeom_algo::make_polygon(&mut model, &corners, true, T)
        .unwrap()
        .shape;
    let tube = ogeom_algo::make_revolution(&mut model, &ring, Axis::Z, std::f64::consts::TAU, T)
        .unwrap()
        .shape;
    let around = |r: f64| -> Vec<Point> {
        (0..24)
            .map(|k| {
                let a = 0.3 + 0.25 * f64::from(k);
                Point::new(r * a.cos(), r * a.sin(), 1.0 + 0.75 * f64::from(k))
            })
            .collect()
    };
    let radial = |p: &Point| Vector::new(p.x, p.y, 0.0);
    let (outer, inner) = (around(10.0), around(8.0));
    let outward: Vec<Vector> = outer.iter().map(radial).collect();
    let inward: Vec<Vector> = inner.iter().map(|p| -radial(p)).collect();
    for method in [RAY, Method::Sphere] {
        for (points, normals) in [(&outer, &outward), (&inner, &inward)] {
            for w in thickness_of_shape(&model, &tube, points, normals, &options(method)).unwrap() {
                let w = w.unwrap();
                assert!((w - 2.0).abs() <= 1e-9, "{method:?} read {w} on the tube");
            }
        }
    }
}
