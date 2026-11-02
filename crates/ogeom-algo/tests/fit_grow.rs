//! Regions grown from clicks over meshes built in code, against the faces
//! the meshes were built from: a bore through a plate, a box with rounded
//! edges, a ball, and a cylinder in two patches beside a cap, every vertex
//! moved by seeded Gaussian noise of 0.02 mm along each axis.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    reason = "test code"
)]

use core::f64::consts::TAU;
use std::collections::BTreeMap;

use ogeom_algo::Cutoff;
use ogeom_algo::fit::{AnyFit, FitKind, GrowOptions, Grown, Side, grow};
use ogeom_core::{OgeomError, Tolerances};
use ogeom_math::{Axis, Direction, Point, Transform, Vector};
use ogeom_mesh::Bvh;
use ogeom_topo::Triangulation;

const T: Tolerances = Tolerances::millimetres();
const SIGMA: f64 = 0.02;

/// A fixed sequence of draws (xorshift), so every run sees the same noise.
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
        (self.next() >> 11) as f64 / (1_u64 << 53) as f64
    }

    /// A standard normal draw, by Box-Muller.
    fn gaussian(&mut self) -> f64 {
        let u = 1.0 - self.unit();
        let v = self.unit();
        (-2.0 * u.ln()).sqrt() * (TAU * v).cos()
    }
}

/// A rigid placement far from the origin about an oblique axis, so nothing
/// lines up with the world's axes by accident.
fn placement() -> Transform {
    let axis = Axis::new(
        Point::ORIGIN,
        Direction::from_coords(0.3, -0.5, 0.8, T).unwrap(),
    );
    Transform::translation(Vector::new(120.0, -75.0, 40.0)) * Transform::rotation(axis, 0.7)
}

/// A mesh with each triangle's face named, built vertex by vertex.
#[derive(Default)]
struct Builder {
    positions: Vec<Point>,
    triangles: Vec<[u32; 3]>,
    face: Vec<u8>,
}

impl Builder {
    fn vertex(&mut self, p: Point) -> u32 {
        self.positions.push(p);
        (self.positions.len() - 1) as u32
    }

    /// A triangle of `face`, wound counter-clockwise about `outward`.
    fn triangle(&mut self, [a, b, c]: [u32; 3], outward: Vector, face: u8) {
        let [pa, pb, pc] = [a, b, c].map(|i| self.positions[i as usize]);
        let tri = if (pb - pa).cross(pc - pa).dot(outward) >= 0.0 {
            [a, b, c]
        } else {
            [a, c, b]
        };
        self.triangles.push(tri);
        self.face.push(face);
    }

    /// The quad `a b c d` as two triangles.
    fn quad(&mut self, [a, b, c, d]: [u32; 4], outward: Vector, face: u8) {
        self.triangle([a, b, c], outward, face);
        self.triangle([a, c, d], outward, face);
    }

    /// The mesh placed and moved by the noise.
    fn mesh(&self, placed: &Transform, seed: u64) -> Triangulation {
        let mut draw = Seeded(seed);
        Triangulation {
            positions: self
                .positions
                .iter()
                .map(|p| {
                    placed.apply(*p)
                        + Vector::new(draw.gaussian(), draw.gaussian(), draw.gaussian()) * SIGMA
                })
                .collect(),
            triangles: self.triangles.clone(),
            ..Triangulation::new()
        }
    }

    fn count(&self, face: u8) -> usize {
        self.face.iter().filter(|f| **f == face).count()
    }

    fn in_region(&self, grown: &Grown, face: u8) -> usize {
        grown
            .triangles
            .iter()
            .filter(|t| self.face[**t as usize] == face)
            .count()
    }
}

/// Rings of `around` vertices on circles about the z axis, one per radius
/// and height.
fn rings(b: &mut Builder, around: usize, at: &[(f64, f64)]) -> Vec<Vec<u32>> {
    at.iter()
        .map(|&(r, z)| {
            (0..around)
                .map(|i| {
                    let a = TAU * i as f64 / around as f64;
                    b.vertex(Point::new(r * a.cos(), r * a.sin(), z))
                })
                .collect()
        })
        .collect()
}

/// Quads between consecutive rings, closing round, facing `outward` at a
/// point.
fn band(b: &mut Builder, rings: &[Vec<u32>], outward: impl Fn(Point) -> Vector, face: u8) {
    for pair in rings.windows(2) {
        let n = pair[0].len();
        for i in 0..n {
            let quad = [
                pair[0][i],
                pair[0][(i + 1) % n],
                pair[1][(i + 1) % n],
                pair[1][i],
            ];
            let p = b.positions[quad[0] as usize];
            b.quad(quad, outward(p), face);
        }
    }
}

const PLATE: u8 = 0;
const WALL: u8 = 1;
const BORE_RADIUS: f64 = 5.0;

/// A plate 8 mm thick, 24 mm across, with a 10 mm bore through it: the top
/// and the bottom (`PLATE`) meet the bore's wall (`WALL`) at sharp edges,
/// sharing its rim vertices. The plate is the material, so the wall faces
/// the bore's axis.
fn bored_plate() -> Builder {
    let (outer, depth, step, around) = (12.0, 8.0, 0.25, 128);
    let mut b = Builder::default();
    let plate_rings = ((outer - BORE_RADIUS) / step) as usize;
    let wall_rings = (depth / step) as usize;
    let radii = |z: f64| -> Vec<(f64, f64)> {
        (0..=plate_rings)
            .map(|j| (BORE_RADIUS + step * j as f64, z))
            .collect()
    };
    let top = rings(&mut b, around, &radii(0.0));
    let bottom = rings(&mut b, around, &radii(-depth));
    let mut wall = vec![top[0].clone()];
    let between: Vec<(f64, f64)> = (1..wall_rings)
        .map(|k| (BORE_RADIUS, -step * k as f64))
        .collect();
    wall.extend(rings(&mut b, around, &between));
    wall.push(bottom[0].clone());
    band(&mut b, &top, |_| Vector::Z, PLATE);
    band(&mut b, &bottom, |_| -Vector::Z, PLATE);
    band(&mut b, &wall, |p| Vector::new(-p.x, -p.y, 0.0), WALL);
    b
}

const FLAT: u8 = 0;
const ROUND: u8 = 1;

/// A cube of half side `a` with every edge rounded to radius `r`, meshed as
/// the cube's faces in `n × n` cells mapped onto it: a point of the cube
/// moves to the point at `r` from the inner cube of half side `a - r`,
/// along the line from its nearest point there. A triangle of the top
/// face's flat (`FLAT`) has every corner within the inner cube's top; the
/// rest are `ROUND`, the side faces included. With `r = a` it is a ball,
/// and with `r = 0` a sharp cube.
fn rounded_cube(a: f64, r: f64, n: usize) -> Builder {
    let mut b = Builder::default();
    let h = 2.0 * a / n as f64;
    let half = (n / 2) as i64;
    let inner = a - r;
    let mut ids: BTreeMap<[i64; 3], u32> = BTreeMap::new();
    let mut at = |b: &mut Builder, g: [i64; 3]| -> u32 {
        *ids.entry(g).or_insert_with(|| {
            let p = g.map(|c| c as f64 * h);
            let q = p.map(|c| c.clamp(-inner, inner));
            let d = Vector::new(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
            let on = if r > 0.0 {
                Point::new(q[0], q[1], q[2]) + d * (r / d.magnitude())
            } else {
                Point::new(p[0], p[1], p[2])
            };
            b.vertex(on)
        })
    };
    let flat_limit = (inner / h).round() as i64;
    for axis in 0..3 {
        for sign in [-1_i64, 1] {
            let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
            let mut outward = [0.0; 3];
            outward[axis] = sign as f64;
            let outward = Vector::new(outward[0], outward[1], outward[2]);
            for i in -half..half {
                for j in -half..half {
                    let corners = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)].map(|(s, t)| {
                        let mut g = [0; 3];
                        g[axis] = sign * half;
                        g[u] = s;
                        g[v] = t;
                        g
                    });
                    let flat = axis == 2
                        && sign == 1
                        && corners
                            .iter()
                            .all(|g| g[0].abs() <= flat_limit && g[1].abs() <= flat_limit);
                    let quad = corners.map(|g| at(&mut b, g));
                    b.quad(quad, outward, if flat { FLAT } else { ROUND });
                }
            }
        }
    }
    b
}

const SHELL: u8 = 0;
const CAP: u8 = 2;
const BOSS_RADIUS: f64 = 4.0;

/// A boss of radius 4 mm in two bands, z in [0, 4] (`SHELL`) and [6, 10]
/// (`WALL`), that share no vertex, the upper one closed by a flat cap
/// (`CAP`). The boss is the material, so its wall faces away from the axis.
fn split_boss() -> Builder {
    let (around, step) = (96, 0.25);
    let mut b = Builder::default();
    let heights = |from: f64| -> Vec<(f64, f64)> {
        (0..=16)
            .map(|k| (BOSS_RADIUS, from + step * k as f64))
            .collect()
    };
    let lower = rings(&mut b, around, &heights(0.0));
    let upper = rings(&mut b, around, &heights(6.0));
    let outward = |p: Point| Vector::new(p.x, p.y, 0.0);
    band(&mut b, &lower, outward, SHELL);
    band(&mut b, &upper, outward, WALL);
    let cap_rings = (BOSS_RADIUS / step) as usize;
    let mut cap = vec![upper[16].clone()];
    let radii: Vec<(f64, f64)> = (1..cap_rings)
        .map(|j| (BOSS_RADIUS - step * j as f64, 10.0))
        .collect();
    cap.extend(rings(&mut b, around, &radii));
    band(&mut b, &cap, |_| Vector::Z, CAP);
    let centre = b.vertex(Point::new(0.0, 0.0, 10.0));
    let last = cap.last().unwrap().clone();
    for i in 0..around {
        b.triangle([centre, last[i], last[(i + 1) % around]], Vector::Z, CAP);
    }
    b
}

fn grown(mesh: &Triangulation, seeds: &[Point], kind: FitKind, options: &GrowOptions) -> Grown {
    let index = Bvh::build(mesh).unwrap();
    grow(mesh, &index, seeds, kind, options, T).unwrap()
}

fn kept(fit: &AnyFit) -> usize {
    fit.used().iter().filter(|u| **u).count()
}

/// Four standard errors of a parameter that moves each point's distance by
/// as much as it moves itself.
fn bound(n: usize) -> f64 {
    4.0 * SIGMA / (n as f64).sqrt()
}

#[test]
fn a_click_in_a_bore_grows_over_its_wall_and_stops_at_the_plate() {
    let placed = placement();
    let b = bored_plate();
    let mesh = b.mesh(&placed, 7);
    let seed = placed.apply(Point::new(-BORE_RADIUS, 0.0, -3.0));
    let grown = grown(&mesh, &[seed], FitKind::Cylinder, &GrowOptions::default());

    let wall = b.in_region(&grown, WALL);
    assert_eq!(
        b.in_region(&grown, PLATE),
        0,
        "the region reached the plate"
    );
    assert!(
        wall as f64 >= 0.95 * b.count(WALL) as f64,
        "{wall} of the wall's {} triangles",
        b.count(WALL)
    );
    let AnyFit::Cylinder(fit) = &grown.fit else {
        panic!("not a cylinder: {:?}", grown.fit.kind());
    };
    let error = fit.shape.radius() - BORE_RADIUS;
    let n = kept(&grown.fit);
    assert!(
        error.abs() <= bound(n),
        "radius off by {error:e} over {n} points, more than {:e}",
        bound(n)
    );
    assert_eq!(fit.side, Some(Side::Outside));
    assert!(
        (fit.sigma / SIGMA - 1.0).abs() <= 0.1,
        "sigma {}",
        fit.sigma
    );
    let report = grown.seeds[0];
    assert!(report.in_region);
    assert!(report.offset.abs() <= 3.0 * SIGMA, "{report:?}");
}

/// The fitted plane's distance from the flat's corners, the largest.
fn plane_error(fit: &AnyFit, placed: &Transform, inner: f64, a: f64) -> f64 {
    let AnyFit::Plane(plane) = fit else {
        panic!("not a plane");
    };
    [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
        .iter()
        .map(|(x, y)| {
            let corner = placed.apply(Point::new(x * inner, y * inner, a));
            plane.shape.signed_distance_to(corner).abs()
        })
        .fold(0.0, f64::max)
}

#[test]
fn peeling_a_box_face_takes_off_its_rounded_edges() {
    let placed = placement();
    let (a, r) = (5.0, 0.5);
    let b = rounded_cube(a, r, 40);
    let mesh = b.mesh(&placed, 11);
    let seed = placed.apply(Point::new(0.3, -0.2, a));

    let options = GrowOptions {
        peel: true,
        ..GrowOptions::default()
    };
    let peeled = grown(&mesh, &[seed], FitKind::Plane, &options);
    assert_eq!(b.in_region(&peeled, ROUND), 0, "a round is left in");
    assert!(peeled.peeled_rings >= 1);
    let flat = b.in_region(&peeled, FLAT);
    assert!(
        flat as f64 >= 0.9 * b.count(FLAT) as f64,
        "{flat} of the flat's {} triangles",
        b.count(FLAT)
    );
    let error = plane_error(&peeled.fit, &placed, a - r, a);
    assert!(
        error <= 0.005,
        "the plane is {error} from the flat's corners"
    );

    // Unpeeled, the rounds' first triangles lie within the noise and join.
    let unpeeled = grown(&mesh, &[seed], FitKind::Plane, &GrowOptions::default());
    assert!(b.in_region(&unpeeled, ROUND) > 0);
    assert_eq!(unpeeled.peeled_rings, 0);
}

#[test]
fn a_plane_clicked_on_a_ball_stays_small_and_fits_badly() {
    let placed = placement();
    let radius = 5.0;
    let b = rounded_cube(radius, radius, 40);
    let mesh = b.mesh(&placed, 13);
    let centre = placed.apply(Point::ORIGIN);
    let top = placed.apply_vector(Vector::Z);
    let seed = centre + top * radius;
    let options = GrowOptions::default();
    let grown = grown(&mesh, &[seed], FitKind::Plane, &options);

    // Every triangle joined within the normal limit of the plane, so of a
    // cap that reaches that far round the ball and a ring of noise more.
    let reach = 2.0 * options.normal_limit;
    for &t in &grown.triangles {
        let [p, q, s] = mesh.triangles[t as usize].map(|v| mesh.positions[v as usize]);
        let towards = (p + ((q - p) + (s - p)) / 3.0) - centre;
        let off = (towards.dot(top) / towards.magnitude()).acos();
        assert!(off <= reach, "triangle {t} lies {off} rad round the ball");
    }
    let share = grown.triangles.len() as f64 / mesh.triangles.len() as f64;
    let cap = (1.0 - reach.cos()) / 2.0;
    assert!(share < cap, "the region holds {share} of the ball");
    assert!(
        grown.fit.sigma() > 5.0 * SIGMA,
        "sigma {} is the noise's, not a ball's",
        grown.fit.sigma()
    );
}

#[test]
fn two_patches_of_one_cylinder_feed_one_fit_and_a_click_off_it_is_left_out() {
    let placed = placement();
    let b = split_boss();
    let mesh = b.mesh(&placed, 17);
    let on = |x: f64, y: f64, z: f64| placed.apply(Point::new(x, y, z));
    let seeds = [
        on(BOSS_RADIUS, 0.0, 2.0),
        on(0.0, -BOSS_RADIUS, 8.0),
        on(1.0, 1.0, 10.0),
    ];
    let both = grown(&mesh, &seeds, FitKind::Cylinder, &GrowOptions::default());
    for face in [SHELL, WALL] {
        let got = b.in_region(&both, face);
        assert!(
            got as f64 >= 0.95 * b.count(face) as f64,
            "face {face}: {got} of {}",
            b.count(face)
        );
    }
    assert_eq!(b.in_region(&both, CAP), 0);
    assert!(both.seeds[0].in_region && both.seeds[1].in_region);
    let cap = both.seeds[2];
    assert!(!cap.in_region, "{cap:?}");
    assert!(
        (cap.angle - core::f64::consts::FRAC_PI_2).abs() < 0.2,
        "{cap:?}"
    );
    let AnyFit::Cylinder(fit) = &both.fit else {
        panic!("not a cylinder");
    };
    let n = kept(&both.fit);
    let error = fit.shape.radius() - BOSS_RADIUS;
    assert!(error.abs() <= bound(n), "radius off by {error:e} over {n}");
    assert_eq!(fit.side, Some(Side::Inside));

    // One click alone grows over its own band only.
    let one = grown(
        &mesh,
        &seeds[..1],
        FitKind::Cylinder,
        &GrowOptions::default(),
    );
    assert_eq!(b.in_region(&one, WALL), 0);
    assert!(b.in_region(&one, SHELL) > 0);
}

#[test]
fn the_region_does_not_depend_on_threads() {
    let placed = placement();
    let mesh = bored_plate().mesh(&placed, 7);
    let seed = placed.apply(Point::new(-BORE_RADIUS, 0.0, -3.0));
    let options = GrowOptions {
        peel: true,
        ..GrowOptions::default()
    };
    ogeom_core::parallel::set_threads(1);
    let one = grown(&mesh, &[seed], FitKind::Cylinder, &options);
    ogeom_core::parallel::set_threads(8);
    let many = grown(&mesh, &[seed], FitKind::Cylinder, &options);
    ogeom_core::parallel::set_threads(0);
    assert_eq!(one, many);
}

#[test]
fn a_cap_on_the_triangles_holds_growth_there() {
    let placed = placement();
    let mesh = bored_plate().mesh(&placed, 7);
    let seed = placed.apply(Point::new(-BORE_RADIUS, 0.0, -3.0));
    let options = GrowOptions {
        max_triangles: Some(500),
        ..GrowOptions::default()
    };
    let capped = grown(&mesh, &[seed], FitKind::Cylinder, &options);
    assert_eq!(capped.triangles.len(), 500);
}

#[test]
fn what_cannot_be_grown_is_refused_by_name() {
    let placed = placement();
    let mesh = bored_plate().mesh(&placed, 7);
    let index = Bvh::build(&mesh).unwrap();
    let seed = placed.apply(Point::new(-BORE_RADIUS, 0.0, -3.0));
    let fine = GrowOptions::default();
    let refused = |seeds: &[Point], kind: FitKind, options: &GrowOptions| {
        grow(&mesh, &index, seeds, kind, options, T).unwrap_err()
    };
    let domain = |e: OgeomError| assert!(matches!(e, OgeomError::Domain(_)), "{e:?}");

    domain(refused(&[], FitKind::Cylinder, &fine));
    domain(refused(
        &[Point::new(f64::NAN, 0.0, 0.0)],
        FitKind::Cylinder,
        &fine,
    ));
    for limit in [0.0, -0.1, 2.0, f64::NAN] {
        let options = GrowOptions {
            normal_limit: limit,
            ..fine
        };
        domain(refused(&[seed], FitKind::Cylinder, &options));
    }
    let all = GrowOptions {
        cutoff: Cutoff::All,
        peel: true,
        ..fine
    };
    domain(refused(&[seed], FitKind::Cylinder, &all));
    let none = GrowOptions {
        max_triangles: Some(0),
        ..fine
    };
    domain(refused(&[seed], FitKind::Cylinder, &none));
    let sided = GrowOptions {
        side: Some(Side::Inside),
        ..fine
    };
    domain(refused(&[seed], FitKind::Plane, &sided));

    let other = Bvh::build(&Triangulation::new()).unwrap();
    let e = grow(&mesh, &other, &[seed], FitKind::Cylinder, &fine, T).unwrap_err();
    assert!(matches!(e, OgeomError::Dimension(_)), "{e:?}");
    let empty = Triangulation::new();
    let e = grow(&empty, &other, &[seed], FitKind::Cylinder, &fine, T).unwrap_err();
    assert!(matches!(e, OgeomError::Construction(_)), "{e:?}");
}

#[test]
fn a_click_beside_an_edge_grows_on_its_own_side() {
    // Half a tenth of a millimetre below the plate's top, in the wall's
    // first row of triangles, whose corners the plate's triangles share.
    let placed = placement();
    let b = bored_plate();
    let mesh = b.mesh(&placed, 9);
    let seed = placed.apply(Point::new(-BORE_RADIUS, 0.0, -0.05));
    let grown = grown(&mesh, &[seed], FitKind::Cylinder, &GrowOptions::default());
    assert_eq!(b.in_region(&grown, PLATE), 0);
    let wall = b.in_region(&grown, WALL);
    assert!(wall as f64 >= 0.95 * b.count(WALL) as f64, "{wall}");
}

const SIDE: u8 = 0;
const END: u8 = 1;

/// A closed cylinder of radius 10 mm and height 30 mm, 64 facets round and
/// 24 rows up, its side (`SIDE`) closed by two fans (`END`), with no noise:
/// every vertex is where `place` puts the point at an angle and a height,
/// then stored as single precision, as a mesh read from an `f32` file is.
fn single_cylinder(place: impl Fn(usize, f64) -> Point) -> Builder {
    let (around, rows, height) = (64, 24, 30.0);
    let mut b = Builder::default();
    let rings: Vec<Vec<u32>> = (0..=rows)
        .map(|k| {
            let z = height * k as f64 / rows as f64;
            (0..around).map(|i| b.vertex(place(i, z))).collect()
        })
        .collect();
    band(&mut b, &rings, |p| Vector::new(p.x, p.y, 0.0), SIDE);
    for (ring, z, up) in [
        (&rings[0], 0.0, -Vector::Z),
        (&rings[rows], height, Vector::Z),
    ] {
        let centre = b.vertex(Point::new(0.0, 0.0, z));
        for i in 0..around {
            b.triangle([centre, ring[i], ring[(i + 1) % around]], up, END);
        }
    }
    b.positions = b.positions.iter().map(|p| single(*p)).collect();
    b
}

/// `p` with each coordinate rounded to single precision.
fn single(p: Point) -> Point {
    let round = |x: f64| f64::from(x as f32);
    Point::new(round(p.x), round(p.y), round(p.z))
}

/// The mesh as built, with no placement and no noise.
fn exact(b: &Builder) -> Triangulation {
    Triangulation {
        positions: b.positions.clone(),
        triangles: b.triangles.clone(),
        ..Triangulation::new()
    }
}

#[test]
fn a_noise_free_single_precision_cylinder_grows_over_its_whole_side() {
    let r = 10.0;
    let wide = single_cylinder(|i, z| {
        let a = TAU * i as f64 / 64.0;
        Point::new(r * a.cos(), r * a.sin(), z)
    });
    let narrow = single_cylinder(|i, z| {
        let a = core::f32::consts::TAU * i as f32 / 64.0;
        Point::new(
            f64::from(10.0_f32 * a.cos()),
            f64::from(10.0_f32 * a.sin()),
            z,
        )
    });
    // Moved off the origin before rounding, the same cylinder rounds
    // another way.
    for (name, b) in [("f64 trig", wide), ("f32 trig", narrow)] {
        for shift in [0.0, 1.9, 517.1] {
            let by = Vector::new(shift, -0.5 * shift, 0.25 * shift);
            let mesh = Triangulation {
                positions: b.positions.iter().map(|p| single(*p + by)).collect(),
                ..exact(&b)
            };
            let seed = Point::new(r, 0.0, 15.0) + by;
            let grown = grown(&mesh, &[seed], FitKind::Cylinder, &GrowOptions::default());
            let side = b.in_region(&grown, SIDE);
            assert!(
                side as f64 >= 0.95 * b.count(SIDE) as f64,
                "{name}, shifted {shift}: {side} of the side's {} triangles",
                b.count(SIDE)
            );
            assert_eq!(
                b.in_region(&grown, END),
                0,
                "{name}, shifted {shift}: an end joined"
            );
            let AnyFit::Cylinder(fit) = &grown.fit else {
                panic!("{name}, shifted {shift}: not a cylinder");
            };
            let error = fit.shape.radius() - r;
            assert!(
                error.abs() <= 0.05,
                "{name}, shifted {shift}: radius off by {error}"
            );
        }
    }
}

#[test]
fn a_noise_free_single_precision_box_face_and_ball_grow_whole_and_peel_nothing() {
    let a = 5.0;
    let cube = rounded_cube(a, 0.5, 40);
    let ball = rounded_cube(a, a, 40);
    for shift in [0.0, 1.9, 13.3, 101.7] {
        let placed =
            Transform::translation(Vector::new(shift, -0.5 * shift, 0.25 * shift)) * placement();
        let rounded = |b: &Builder| Triangulation {
            positions: b
                .positions
                .iter()
                .map(|p| single(placed.apply(*p)))
                .collect(),
            ..exact(b)
        };

        let mesh = rounded(&cube);
        let seed = placed.apply(Point::new(0.3, -0.2, a));
        for peel in [false, true] {
            let options = GrowOptions {
                peel,
                ..GrowOptions::default()
            };
            let face = grown(&mesh, &[seed], FitKind::Plane, &options);
            assert_eq!(
                cube.in_region(&face, FLAT),
                cube.count(FLAT),
                "shifted {shift}, peel {peel}: the flat is not whole"
            );
            assert_eq!(
                cube.in_region(&face, ROUND),
                0,
                "shifted {shift}, peel {peel}"
            );
            assert_eq!(face.peeled_rings, 0, "shifted {shift}, peel {peel}");
        }

        let mesh = rounded(&ball);
        let seed = placed.apply(Point::new(0.0, 0.0, a));
        let grown = grown(&mesh, &[seed], FitKind::Sphere, &GrowOptions::default());
        assert_eq!(
            grown.triangles.len(),
            mesh.triangles.len(),
            "shifted {shift}: the ball is not whole"
        );
        let AnyFit::Sphere(fit) = &grown.fit else {
            panic!("not a sphere");
        };
        let error = fit.shape.radius() - a;
        assert!(
            error.abs() <= 1e-4,
            "shifted {shift}: radius off by {error}"
        );
    }
}

/// The top rim of [`bored_plate`]: the first ring built.
fn top_rim() -> Vec<u32> {
    (0..128).collect()
}

/// The chain's vertices, sorted.
fn sorted(grown: &Grown) -> Vec<u32> {
    let mut v = grown.vertices.clone();
    v.sort_unstable();
    v
}

#[test]
fn a_click_by_a_bore_rim_traces_the_rim_circle() {
    let placed = placement();
    let b = bored_plate();
    let mesh = b.mesh(&placed, 7);
    let centre = placed.apply(Point::ORIGIN);
    let axis = placed.apply_vector(Vector::Z);
    // On the plate a tenth of a millimetre out from the rim, and on the
    // wall six tenths down from it.
    for at in [
        Point::new(-BORE_RADIUS - 0.1, 0.0, 0.0),
        Point::new(0.0, BORE_RADIUS, -0.6),
    ] {
        let seed = placed.apply(at);
        let grown = grown(&mesh, &[seed], FitKind::Circle3, &GrowOptions::default());
        assert!(grown.triangles.is_empty());
        assert_eq!(sorted(&grown), top_rim(), "seed {at:?}");
        let AnyFit::Circle(fit) = &grown.fit else {
            panic!("not a circle: {:?}", grown.fit.kind());
        };
        let n = kept(&grown.fit);
        assert!(n >= 125, "{n} of the rim's 128 vertices kept");
        let error = fit.shape.radius() - BORE_RADIUS;
        assert!(
            error.abs() <= bound(n),
            "radius off by {error:e} over {n} points, more than {:e}",
            bound(n)
        );
        let off = (fit.shape.centre() - centre).magnitude();
        assert!(
            off <= bound(n),
            "centre off by {off:e} over {n} points, more than {:e}",
            bound(n)
        );
        let tilt = fit.shape.frame().z().vector().cross(axis).magnitude();
        assert!(tilt <= 0.01, "axis tilted by {tilt}");
        assert_eq!(fit.side, Some(Side::Outside));
        let report = grown.seeds[0];
        assert!(report.in_region, "{report:?}");
        assert!(
            report.angle > GrowOptions::default().crease_limit,
            "{report:?}"
        );
    }
}

/// The fitted line's distance from a point and its angle from a direction.
fn line_error(fit: &AnyFit, at: Point, along: Vector) -> (f64, f64) {
    let AnyFit::Line(line) = fit else {
        panic!("not a line: {:?}", fit.kind());
    };
    let d = line.shape.direction.vector();
    let w = at - line.shape.location;
    let off = (w - d * w.dot(d)).magnitude();
    (off, d.cross(along).magnitude().asin())
}

#[test]
fn a_click_by_a_box_edge_traces_its_line_from_corner_to_corner() {
    let placed = placement();
    let (a, n) = (5.0, 20);
    let b = rounded_cube(a, 0.0, n);
    let mesh = b.mesh(&placed, 19);
    let seed = placed.apply(Point::new(0.3, a - 0.2, a));
    let grown = grown(&mesh, &[seed], FitKind::Line3, &GrowOptions::default());
    assert!(grown.triangles.is_empty());

    // Every vertex of the edge y = z = a, the corners at the chain's ends.
    let on_edge = |v: u32| {
        let p = b.positions[v as usize];
        p.y == a && p.z == a
    };
    assert_eq!(grown.vertices.len(), n + 1, "{:?}", grown.vertices);
    assert!(grown.vertices.iter().all(|&v| on_edge(v)));
    let ends = [grown.vertices[0], grown.vertices[n]].map(|v| b.positions[v as usize].x);
    assert_eq!(ends[0].abs(), a);
    assert_eq!(ends[0], -ends[1]);

    // Four standard errors: across the line at its middle, the noise over
    // the points, each way; its tilt, the noise over the points' spread
    // along it, `L / √12` for points spread evenly over a length `L`.
    let m = kept(&grown.fit) as f64;
    let (off, angle) = line_error(
        &grown.fit,
        placed.apply(Point::new(0.0, a, a)),
        placed.apply_vector(Vector::X),
    );
    let across = 4.0 * SIGMA * (2.0 / m).sqrt();
    assert!(
        off <= across,
        "the line is {off:e} off the edge, past {across:e}"
    );
    let tilt = 4.0 * SIGMA * (24.0 / m).sqrt() / (2.0 * a);
    assert!(
        angle <= tilt,
        "the line is {angle:e} off the edge's direction, past {tilt:e}"
    );
    assert!(grown.seeds[0].in_region);
}

/// A pipe's wall of radius 10 mm and length 20 mm, open at both ends, 64
/// facets round and 20 rows along: triangles about a millimetre across,
/// fifty times the noise, so no edge turns by more than a few degrees.
fn pipe() -> Builder {
    let mut b = Builder::default();
    let at: Vec<(f64, f64)> = (0..=20).map(|k| (10.0, f64::from(k))).collect();
    let rings = rings(&mut b, 64, &at);
    band(&mut b, &rings, |p| Vector::new(p.x, p.y, 0.0), SIDE);
    b
}

#[test]
fn a_click_on_a_smooth_surface_finds_no_crease_and_is_refused() {
    let placed = placement();
    let b = pipe();
    let ball = rounded_cube(5.0, 5.0, 40);
    let noisy = b.mesh(&placed, 13);
    let exact_ball = Triangulation {
        positions: ball
            .positions
            .iter()
            .map(|p| single(placed.apply(*p)))
            .collect(),
        ..exact(&ball)
    };
    for (mesh, at) in [
        (&noisy, Point::new(10.0, 0.0, 10.0)),
        (&exact_ball, Point::new(0.0, 0.0, 5.0)),
    ] {
        let index = Bvh::build(mesh).unwrap();
        let seed = placed.apply(at);
        for kind in [FitKind::Circle3, FitKind::Line3] {
            let e = grow(mesh, &index, &[seed], kind, &GrowOptions::default(), T).unwrap_err();
            assert!(
                matches!(e, OgeomError::Construction(_)) && e.to_string().contains("smooth"),
                "{e:?}"
            );
        }
    }
}

#[test]
fn a_noise_free_single_precision_rim_and_edge_trace_whole() {
    let b = bored_plate();
    let cube = rounded_cube(5.0, 0.0, 20);
    // Moved by 221.1 before rounding, the edge's rounding stops a gate
    // without the floor four vertices in.
    for shift in [0.0, 1.9, 221.1, 517.1] {
        let placed =
            Transform::translation(Vector::new(shift, -0.5 * shift, 0.25 * shift)) * placement();
        let rounded = |b: &Builder| Triangulation {
            positions: b
                .positions
                .iter()
                .map(|p| single(placed.apply(*p)))
                .collect(),
            ..exact(b)
        };

        let mesh = rounded(&b);
        let seed = placed.apply(Point::new(-BORE_RADIUS - 0.1, 0.0, 0.0));
        let rim = grown(&mesh, &[seed], FitKind::Circle3, &GrowOptions::default());
        assert_eq!(sorted(&rim), top_rim(), "shifted {shift}");
        assert_eq!(kept(&rim.fit), 128, "shifted {shift}");
        let AnyFit::Circle(fit) = &rim.fit else {
            panic!("not a circle");
        };
        let error = fit.shape.radius() - BORE_RADIUS;
        assert!(
            error.abs() <= 1e-4,
            "shifted {shift}: radius off by {error}"
        );

        let mesh = rounded(&cube);
        let seed = placed.apply(Point::new(0.3, 4.8, 5.0));
        let edge = grown(&mesh, &[seed], FitKind::Line3, &GrowOptions::default());
        assert_eq!(
            edge.vertices.len(),
            21,
            "shifted {shift}: {:?}",
            edge.vertices
        );
        assert_eq!(kept(&edge.fit), 21, "shifted {shift}");
        let ends = [edge.vertices[0], edge.vertices[20]].map(|v| cube.positions[v as usize]);
        for p in ends {
            assert_eq!((p.x.abs(), p.y, p.z), (5.0, 5.0, 5.0), "shifted {shift}");
        }
    }
}

#[test]
fn a_traced_crease_does_not_depend_on_threads() {
    let placed = placement();
    let mesh = bored_plate().mesh(&placed, 7);
    let seed = placed.apply(Point::new(-BORE_RADIUS - 0.1, 0.0, 0.0));
    let cube = rounded_cube(5.0, 0.0, 20).mesh(&placed, 19);
    let edge = placed.apply(Point::new(0.3, 4.8, 5.0));
    let both = || {
        (
            grown(&mesh, &[seed], FitKind::Circle3, &GrowOptions::default()),
            grown(&cube, &[edge], FitKind::Line3, &GrowOptions::default()),
        )
    };
    ogeom_core::parallel::set_threads(1);
    let one = both();
    ogeom_core::parallel::set_threads(8);
    let many = both();
    ogeom_core::parallel::set_threads(0);
    assert_eq!(one, many);
}

#[test]
fn what_cannot_be_traced_is_refused_by_name() {
    let placed = placement();
    let mesh = bored_plate().mesh(&placed, 7);
    let index = Bvh::build(&mesh).unwrap();
    let seed = placed.apply(Point::new(-BORE_RADIUS - 0.1, 0.0, 0.0));
    let fine = GrowOptions::default();
    let domain = |kind: FitKind, options: &GrowOptions| {
        let e = grow(&mesh, &index, &[seed], kind, options, T).unwrap_err();
        assert!(matches!(e, OgeomError::Domain(_)), "{e:?}");
    };
    for limit in [0.0, -0.1, core::f64::consts::PI, f64::NAN] {
        let options = GrowOptions {
            crease_limit: limit,
            ..fine
        };
        domain(FitKind::Circle3, &options);
    }
    let peel = GrowOptions { peel: true, ..fine };
    domain(FitKind::Circle3, &peel);
    let sided = GrowOptions {
        side: Some(Side::Inside),
        ..fine
    };
    domain(FitKind::Line3, &sided);

    // A circle along a box's exactly straight edges: the circle fit's own
    // refusal, not a smooth part's.
    let cube = exact(&rounded_cube(5.0, 0.0, 20));
    let index = Bvh::build(&cube).unwrap();
    let edge = Point::new(0.3, 4.8, 5.0);
    let e = grow(&cube, &index, &[edge], FitKind::Circle3, &fine, T).unwrap_err();
    assert!(
        matches!(e, OgeomError::Construction(_)) && e.to_string().contains("one line"),
        "{e:?}"
    );
}
