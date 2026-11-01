//! Rigid registration against known motions: the closed form carries
//! noise-free pairs exactly, a scan displaced by a random rigid motion
//! registers back onto its mesh, a flat selection reports what it leaves
//! free, and a part's symmetry plane is found from a seed off it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_algo::Cutoff;
use ogeom_algo::register::{Icp, Starts, absolute_orientation, icp, symmetry_plane};
use ogeom_core::{FastMap, Tolerances};
use ogeom_math::{Axis, Direction, Plane, Point, Transform, Vector};
use ogeom_mesh::Bvh;
use ogeom_topo::Triangulation;

const T: Tolerances = Tolerances::millimetres();

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

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let unit = (self.next() >> 11) as f64 / (1_u64 << 53) as f64;
        lo + (hi - lo) * unit
    }

    fn direction(&mut self) -> Vector {
        loop {
            let v = Vector::new(
                self.range(-1.0, 1.0),
                self.range(-1.0, 1.0),
                self.range(-1.0, 1.0),
            );
            let m = v.magnitude();
            if m > 0.1 && m <= 1.0 {
                return v / m;
            }
        }
    }

    /// A standard normal deviate, by Box and Muller.
    fn gaussian(&mut self) -> f64 {
        let u = self.range(f64::MIN_POSITIVE, 1.0);
        let v = self.range(0.0, std::f64::consts::TAU);
        (-2.0 * u.ln()).sqrt() * v.cos()
    }

    /// A rigid motion turning by up to `angle` about a random axis through
    /// `centre` and moving by up to `shift`.
    fn motion(&mut self, centre: Point, angle: f64, shift: f64) -> Transform {
        let axis = Axis::new(centre, Direction::new(self.direction(), T).unwrap());
        let turn = Transform::rotation(axis, self.range(-angle, angle));
        Transform::translation(self.direction() * self.range(0.0, shift)) * turn
    }
}

fn unit(v: Vector) -> Vector {
    v / v.magnitude()
}

/// An icosahedron subdivided `levels` times onto the unit sphere.
fn icosphere(levels: u32) -> (Vec<Vector>, Vec<[u32; 3]>) {
    let g = f64::midpoint(1.0, 5.0_f64.sqrt());
    let mut positions: Vec<Vector> = [
        (-1.0, g, 0.0),
        (1.0, g, 0.0),
        (-1.0, -g, 0.0),
        (1.0, -g, 0.0),
        (0.0, -1.0, g),
        (0.0, 1.0, g),
        (0.0, -1.0, -g),
        (0.0, 1.0, -g),
        (g, 0.0, -1.0),
        (g, 0.0, 1.0),
        (-g, 0.0, -1.0),
        (-g, 0.0, 1.0),
    ]
    .iter()
    .map(|&(x, y, z)| unit(Vector::new(x, y, z)))
    .collect();
    let mut triangles: Vec<[u32; 3]> = vec![
        [0, 11, 5],
        [0, 5, 1],
        [0, 1, 7],
        [0, 7, 10],
        [0, 10, 11],
        [1, 5, 9],
        [5, 11, 4],
        [11, 10, 2],
        [10, 7, 6],
        [7, 1, 8],
        [3, 9, 4],
        [3, 4, 2],
        [3, 2, 6],
        [3, 6, 8],
        [3, 8, 9],
        [4, 9, 5],
        [2, 4, 11],
        [6, 2, 10],
        [8, 6, 7],
        [9, 8, 1],
    ];
    for _ in 0..levels {
        let mut middles: FastMap<(u32, u32), u32> = FastMap::default();
        let mut middle = |a: u32, b: u32, positions: &mut Vec<Vector>| {
            *middles.entry((a.min(b), a.max(b))).or_insert_with(|| {
                positions.push(unit(positions[a as usize] + positions[b as usize]));
                u32::try_from(positions.len() - 1).unwrap()
            })
        };
        let mut next = Vec::with_capacity(triangles.len() * 4);
        for [a, b, c] in triangles {
            let ab = middle(a, b, &mut positions);
            let bc = middle(b, c, &mut positions);
            let ca = middle(c, a, &mut positions);
            next.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
        }
        triangles = next;
    }
    (positions, triangles)
}

/// A closed blob some 100 mm long: an icosphere's directions `u` stretched
/// to semi-axes 50, 32 and 22 mm and scaled by `1 + bump(u)`.
fn blob(levels: u32, bump: impl Fn(Vector) -> f64) -> Triangulation {
    let (directions, triangles) = icosphere(levels);
    let mut mesh = Triangulation::new();
    mesh.positions = directions
        .iter()
        .map(|&u| {
            let s = 1.0 + bump(u);
            Point::new(50.0 * u.x * s, 32.0 * u.y * s, 22.0 * u.z * s)
        })
        .collect();
    mesh.triangles = triangles;
    mesh
}

/// A blob with no symmetry.
fn lopsided(levels: u32) -> Triangulation {
    blob(levels, |u| {
        0.25 * u.x * u.y + 0.15 * u.z + 0.1 * u.x * u.z * u.z
    })
}

/// A rounded block some 180 by 160 by 140 mm with no symmetry: an
/// icosphere's directions pushed towards a box's faces (each coordinate's
/// size raised to the power one half) and scaled by a lopsided bump.
fn block(levels: u32) -> Triangulation {
    let (directions, triangles) = icosphere(levels);
    let squared = |c: f64| c.signum() * c.abs().sqrt();
    let mut mesh = Triangulation::new();
    mesh.positions = directions
        .iter()
        .map(|&u| {
            let s = 1.0 + 0.15 * u.x * u.y + 0.1 * u.z + 0.08 * u.x * u.z * u.z;
            Point::new(
                90.0 * squared(u.x) * s,
                80.0 * squared(u.y) * s,
                70.0 * squared(u.z) * s,
            )
        })
        .collect();
    mesh.triangles = triangles;
    mesh
}

/// A blob whose one symmetry is the mirror in `x = 0`: its bump is even in
/// `x`, and the icosphere's vertices are symmetric in each coordinate.
fn mirrored(levels: u32) -> Triangulation {
    blob(levels, |u| {
        0.2 * u.y + 0.1 * u.z + 0.15 * u.y * u.z + 0.1 * u.x * u.x
    })
}

fn corners(mesh: &Triangulation, t: usize) -> [Point; 3] {
    mesh.triangles[t].map(|k| mesh.positions[k as usize])
}

fn unit_normal(t: [Point; 3]) -> Vector {
    unit((t[1] - t[0]).cross(t[2] - t[0]))
}

/// `count` points on the mesh's triangles, chosen by area, each moved off
/// along its triangle's normal by Gaussian noise of deviation `sigma`,
/// with that normal.
fn scan(
    mesh: &Triangulation,
    count: usize,
    sigma: f64,
    rng: &mut Seeded,
) -> (Vec<Point>, Vec<Vector>) {
    let mut cumulative = Vec::with_capacity(mesh.triangles.len());
    let mut total = 0.0;
    for t in 0..mesh.triangles.len() {
        let [a, b, c] = corners(mesh, t);
        total += (b - a).cross(c - a).magnitude();
        cumulative.push(total);
    }
    (0..count)
        .map(|_| {
            let at = rng.range(0.0, total);
            let t = cumulative
                .partition_point(|&c| c < at)
                .min(cumulative.len() - 1);
            let [a, b, c] = corners(mesh, t);
            let (mut u, mut v) = (rng.range(0.0, 1.0), rng.range(0.0, 1.0));
            if u + v > 1.0 {
                (u, v) = (1.0 - u, 1.0 - v);
            }
            let n = unit_normal([a, b, c]);
            (
                a + (b - a) * u + (c - a) * v + n * (sigma * rng.gaussian()),
                n,
            )
        })
        .unzip()
}

/// The angle `r`'s rotation turns by, from its skew part.
fn turn_of(r: &Transform) -> f64 {
    let m = r.linear().rows;
    let skew = Vector::new(m[2][1] - m[1][2], m[0][2] - m[2][0], m[1][0] - m[0][1]);
    let cos = (m[0][0] + m[1][1] + m[2][2] - 1.0) * 0.5;
    (skew.magnitude() * 0.5).atan2(cos)
}

/// The farthest `t` moves any of the points.
fn farthest_move(t: &Transform, points: &[Point]) -> f64 {
    points
        .iter()
        .map(|&p| (t.apply(p) - p).magnitude())
        .fold(0.0, f64::max)
}

fn diagonal(mesh: &Triangulation) -> f64 {
    ogeom_math::Aabb::of_points(&mesh.positions).diagonal()
}

/// A set of points drawn from a generator.
type Shape = fn(&mut Seeded) -> Point;

/// Noise-free pairs are carried exactly: generic sets, sets far from the
/// origin, coplanar and collinear sets (where an unfixed singular-value
/// solution may reflect), sets a thousandth as thick as they are long
/// off a line or a plane, and weighted sets whose
/// zero-weight pairs are wild.
#[test]
fn absolute_orientation_carries_noise_free_pairs_exactly() {
    let mut rng = Seeded(0x0ab5_017e);
    let shapes: [(&str, Shape); 6] = [
        ("generic", |r| {
            Point::new(
                r.range(-50.0, 50.0),
                r.range(-30.0, 30.0),
                r.range(-20.0, 20.0),
            )
        }),
        ("far", |r| {
            Point::new(
                1e4 + r.range(-5.0, 5.0),
                -2e4 + r.range(-5.0, 5.0),
                3e3 + r.range(-5.0, 5.0),
            )
        }),
        ("coplanar", |r| {
            Point::new(r.range(-50.0, 50.0), r.range(-30.0, 30.0), 7.0)
        }),
        ("collinear", |r| {
            let s = r.range(-50.0, 50.0);
            Point::new(1.0 + s, 2.0 - 0.5 * s, 3.0 + 0.25 * s)
        }),
        ("nearly collinear", |r| {
            let s = r.range(-50.0, 50.0);
            Point::new(s, 0.05 * r.range(-1.0, 1.0), 0.05 * r.range(-1.0, 1.0))
        }),
        ("nearly coplanar", |r| {
            Point::new(
                r.range(-50.0, 50.0),
                r.range(-30.0, 30.0),
                0.05 * r.range(-1.0, 1.0),
            )
        }),
    ];
    for (name, shape) in shapes {
        for count in [3, 4, 10, 200] {
            for _ in 0..20 {
                let motion = rng.motion(Point::new(5.0, -3.0, 2.0), std::f64::consts::PI, 100.0);
                let from: Vec<Point> = (0..count).map(|_| shape(&mut rng)).collect();
                let pairs: Vec<(Point, Point)> =
                    from.iter().map(|&p| (p, motion.apply(p))).collect();
                let found = absolute_orientation(&pairs, None).unwrap();
                assert!(found.preserves_handedness(), "{name}: a reflection");
                let size = from
                    .iter()
                    .map(|p| p.to_vector().magnitude())
                    .fold(1.0, f64::max);
                let worst = pairs
                    .iter()
                    .map(|&(a, b)| (found.apply(a) - b).magnitude())
                    .fold(0.0, f64::max);
                assert!(
                    worst <= 1e-12 * size,
                    "{name}, {count} pairs: off by {worst:e} at size {size}"
                );
            }
        }
    }
    // Weighted: the pairs of weight zero are wild and change nothing, and
    // integer weights equal repeating the pairs.
    let motion = rng.motion(Point::ORIGIN, 2.0, 50.0);
    let from: Vec<Point> = (0..12).map(|_| shapes[0].1(&mut rng)).collect();
    let mut pairs: Vec<(Point, Point)> = from.iter().map(|&p| (p, motion.apply(p))).collect();
    let mut weights = vec![1.0; pairs.len()];
    for k in 0..4 {
        pairs[k].1 += rng.direction() * 30.0;
        weights[k] = 0.0;
    }
    let found = absolute_orientation(&pairs, Some(&weights)).unwrap();
    for (&(a, _), &p) in pairs.iter().zip(&from).skip(4) {
        assert!((found.apply(a) - motion.apply(p)).magnitude() <= 1e-12 * 60.0);
    }
    let noisy: Vec<(Point, Point)> = from
        .iter()
        .map(|&p| (p, motion.apply(p) + rng.direction() * 0.5))
        .collect();
    let twice: Vec<(Point, Point)> = noisy.iter().chain(&noisy[..5]).copied().collect();
    let mut doubled = vec![1.0; noisy.len()];
    doubled[..5].fill(2.0);
    let (a, b) = (
        absolute_orientation(&noisy, Some(&doubled)).unwrap(),
        absolute_orientation(&twice, None).unwrap(),
    );
    assert!(farthest_move(&(a * b.inverse().unwrap()), &from) <= 1e-11);
    assert!(absolute_orientation(&[], None).is_err());
    assert!(absolute_orientation(&noisy, Some(&[1.0])).is_err());
    assert!(absolute_orientation(&noisy[..1], Some(&[-1.0])).is_err());
    assert!(absolute_orientation(&noisy[..1], Some(&[0.0])).is_err());
    let one = absolute_orientation(&noisy[..1], None).unwrap();
    assert_eq!(one.linear(), ogeom_math::Matrix3::IDENTITY);
}

fn farthest_from_centroid(points: &[Point]) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let c = Point::ORIGIN
        + points.iter().fold(Vector::ZERO, |s, p| s + p.to_vector()) / points.len() as f64;
    points
        .iter()
        .map(|&p| (p - c).magnitude())
        .fold(0.0, f64::max)
}

/// The deviation of the rotation and of the translation that `σ` of noise
/// along the normals leaves a point-to-plane registration of exact pairs
/// with: the largest over directions of `σ` times the square root of the
/// inverse normal matrix's diagonal blocks, about the points' centroid.
fn noise_floor(points: &[Point], normals: &[Vector], sigma: f64) -> (f64, f64) {
    #[allow(clippy::cast_precision_loss)]
    let c = points.iter().fold(Vector::ZERO, |s, p| s + p.to_vector()) / points.len() as f64;
    let mut a = nalgebra::Matrix6::<f64>::zeros();
    for (p, n) in points.iter().zip(normals) {
        let r = (p.to_vector() - c).cross(*n);
        let row = nalgebra::Vector6::new(r.x, r.y, r.z, n.x, n.y, n.z);
        a += row * row.transpose();
    }
    let covariance = a.try_inverse().unwrap() * (sigma * sigma);
    let worst = |k: usize| {
        let block = covariance.fixed_view::<3, 3>(k, k).into_owned();
        nalgebra::SymmetricEigen::new(block)
            .eigenvalues
            .max()
            .sqrt()
    };
    (worst(0), worst(3))
}

/// A scan of a rounded block displaced by random rigid motions of up to
/// 30° and a tenth of its diagonal registers back from the identity, the
/// same at one thread and at eight: noise-free, to within a nanometre;
/// with 0.02 mm of noise, to within a micron and 1e-5 rad, and within four
/// deviations of what the noise leaves any registration of these samples.
#[test]
fn a_displaced_scan_registers_back_onto_its_mesh() {
    let mesh = block(6);
    let index = Bvh::build(&mesh).unwrap();
    let mut rng = Seeded(0x1c9_0b1b);
    let reach = 0.1 * diagonal(&mesh);
    for sigma in [0.0, 0.02] {
        let (points, normals) = scan(&mesh, 100_000, sigma, &mut rng);
        let (turn_floor, shift_floor) = noise_floor(&points, &normals, sigma);
        let radius = farthest_from_centroid(&points);
        for _ in 0..4 {
            let displaced = rng.motion(Point::new(3.0, -2.0, 1.0), 30_f64.to_radians(), reach);
            let source: Vec<Point> = points.iter().map(|&p| displaced.apply(p)).collect();
            let turned: Vec<Vector> = normals.iter().map(|&n| displaced.apply_vector(n)).collect();
            let opts = Icp::default();
            let mut answers = Vec::new();
            for threads in [1, 8] {
                ogeom_core::parallel::set_threads(threads);
                answers.push(icp(&source, Some(&turned), &index, &opts).unwrap());
            }
            ogeom_core::parallel::set_threads(0);
            assert_eq!(
                answers[0], answers[1],
                "the thread count changed the answer"
            );
            let found = &answers[0];
            assert!(found.converged, "σ {sigma}: not converged: {found:?}");
            let error = found.transform * displaced;
            let (turn, moved) = (turn_of(&error), farthest_move(&error, &points));
            if sigma == 0.0 {
                assert!(
                    turn <= 1e-11 && moved <= 1e-9,
                    "back to within {turn:e} rad and {moved:e} mm"
                );
            } else {
                assert!(
                    turn <= 1e-5 && moved <= 1e-3,
                    "back to within {turn:e} rad and {moved:e} mm"
                );
                assert!(
                    turn <= 4.0 * turn_floor,
                    "{turn:e} rad against a floor of {turn_floor:e}"
                );
                let floor = turn_floor * radius + shift_floor;
                assert!(
                    moved <= 4.0 * floor,
                    "{moved:e} mm against a floor of {floor:e}"
                );
            }
            assert!(
                found.free_directions.is_empty(),
                "{:?}",
                found.free_directions
            );
            assert!(found.rms <= 2.0 * sigma + 1e-9, "rms {}", found.rms);
            assert_eq!(
                (found.sampled, found.inliers <= found.sampled),
                (100_000, true)
            );
        }
    }
}

/// Point to point registers back too, and the principal axes' starts
/// recover a half turn that the identity alone does not: the score picks
/// the start that pairs every sample closely.
#[test]
fn point_to_point_and_the_principal_starts_register_back() {
    let mesh = lopsided(5);
    let index = Bvh::build(&mesh).unwrap();
    let mut rng = Seeded(0x5_7a27);
    let (points, _) = scan(&mesh, 5_000, 0.0, &mut rng);
    let displaced = rng.motion(Point::ORIGIN, 20_f64.to_radians(), 5.0);
    let source: Vec<Point> = points.iter().map(|&p| displaced.apply(p)).collect();
    let opts = Icp {
        point_to_plane: false,
        max_iter: 400,
        ..Icp::default()
    };
    let found = icp(&source, None, &index, &opts).unwrap();
    let error = found.transform * displaced;
    assert!(
        turn_of(&error) <= 1e-4 && farthest_move(&error, &points) <= 1e-2,
        "{found:?}"
    );

    let axis = Axis::new(
        Point::new(1.0, 2.0, 0.0),
        Direction::new(Vector::new(0.2, 0.3, 1.0), T).unwrap(),
    );
    let flipped =
        Transform::translation(Vector::new(4.0, -7.0, 2.0)) * Transform::rotation(axis, 3.0);
    let source: Vec<Point> = points.iter().map(|&p| flipped.apply(p)).collect();
    let alone = icp(&source, None, &index, &Icp::default()).unwrap();
    assert!(
        turn_of(&(alone.transform * flipped)) > 0.1,
        "the identity found it alone"
    );
    let opts = Icp {
        starts: Starts::PrincipalAxes24,
        ..Icp::default()
    };
    let found = icp(&source, None, &index, &opts).unwrap();
    let error = found.transform * flipped;
    assert!(
        turn_of(&error) <= 1e-5 && farthest_move(&error, &points) <= 1e-3,
        "{found:?}"
    );
    assert!(found.score < alone.score);
    // Given starts: the truth's inverse among wrong ones wins.
    let truth = flipped.inverse().unwrap();
    let opts = Icp {
        starts: Starts::Given(vec![Transform::IDENTITY, truth, flipped]),
        ..Icp::default()
    };
    let found = icp(&source, None, &index, &opts).unwrap();
    assert!(farthest_move(&(found.transform * flipped), &points) <= 1e-3);
}

/// A box's faces, each an `n` by `n` grid of quads, wound outward.
fn box_mesh(size: [f64; 3], n: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    let at = |axis: usize, i: u32| size[axis] * f64::from(i) / f64::from(n);
    for axis in 0..3 {
        let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
        for high in [false, true] {
            let first = u32::try_from(mesh.positions.len()).unwrap();
            for i in 0..=n {
                for j in 0..=n {
                    let mut c = [0.0; 3];
                    c[axis] = if high { size[axis] } else { 0.0 };
                    c[u] = at(u, i);
                    c[v] = at(v, j);
                    mesh.positions.push(Point::new(c[0], c[1], c[2]));
                }
            }
            let index = |i: u32, j: u32| first + i * (n + 1) + j;
            for i in 0..n {
                for j in 0..n {
                    let (a, b, c, d) = (
                        index(i, j),
                        index(i + 1, j),
                        index(i + 1, j + 1),
                        index(i, j + 1),
                    );
                    if high {
                        mesh.triangles.extend([[a, b, c], [a, c, d]]);
                    } else {
                        mesh.triangles.extend([[a, c, b], [a, d, c]]);
                    }
                }
            }
        }
    }
    mesh
}

/// Samples of a box's top face alone, noisy, leave two slides in its plane
/// and a spin about its normal free, and the registration moves them only
/// back onto the plane; the whole box leaves nothing free.
#[test]
fn a_flat_selection_reports_two_slides_and_a_spin_free() {
    let size = [40.0, 30.0, 20.0];
    let mesh = box_mesh(size, 8);
    let index = Bvh::build(&mesh).unwrap();
    let mut rng = Seeded(0xf1a7);
    let top = |sigma: f64, rng: &mut Seeded| -> Vec<Point> {
        (0..4_000)
            .map(|_| {
                Point::new(
                    rng.range(2.0, 38.0),
                    rng.range(2.0, 28.0),
                    20.0 + sigma * rng.gaussian(),
                )
            })
            .collect()
    };
    let lifted = Transform::translation(Vector::new(0.3, -0.2, 0.5));
    let z = Vector::new(0.0, 0.0, 1.0);
    for sigma in [0.0, 0.02] {
        let source: Vec<Point> = top(sigma, &mut rng)
            .iter()
            .map(|&p| lifted.apply(p))
            .collect();
        let found = icp(&source, None, &index, &Icp::default()).unwrap();
        assert!(found.converged);
        let free = &found.free_directions;
        assert_eq!(free.len(), 3, "σ {sigma}: {free:?}");
        for slide in &free[..2] {
            assert!(
                slide.rotation.magnitude() * slide.radius <= 1e-9,
                "{slide:?}"
            );
            assert!((slide.translation.magnitude() - 1.0).abs() <= 1e-9);
            assert!(slide.translation.dot(z).abs() <= 1e-9, "{slide:?}");
        }
        let spin = &free[2];
        assert!(spin.translation.magnitude() <= 1e-9, "{spin:?}");
        assert!((spin.rotation.magnitude() * spin.radius - 1.0).abs() <= 1e-9);
        assert!(
            (unit(spin.rotation).dot(z).abs() - 1.0).abs() <= 1e-9,
            "{spin:?}"
        );
        // Back down onto the plane, neither slid nor spun: the samples'
        // centroid moves only along the normal, the turn is about an axis
        // in the plane, and noise-free, the lift's slide is left as it was.
        let moved = found.transform * lifted;
        let m = moved.linear().rows;
        // Turns about two axes in the plane compose to a spin of the order
        // of their product, and no more.
        let tilt = turn_of(&moved);
        assert!(
            (m[1][0] - m[0][1]).abs() * 0.5 <= tilt * tilt + 1e-15,
            "{moved:?}"
        );
        #[allow(clippy::cast_precision_loss)]
        let centre = Point::ORIGIN
            + source.iter().fold(Vector::ZERO, |s, p| s + p.to_vector()) / source.len() as f64;
        let back = found.transform.apply(centre) - centre;
        assert!(back.x.abs() <= 1e-9 && back.y.abs() <= 1e-9, "{back:?}");
        if sigma == 0.0 {
            let left = moved.translation_vector() - Vector::new(0.3, -0.2, 0.0);
            assert!(
                turn_of(&moved) <= 1e-12 && left.magnitude() <= 1e-9,
                "{moved:?}"
            );
        }
    }

    let (all, _) = scan(&mesh, 6_000, 0.02, &mut rng);
    let found = icp(&all, None, &index, &Icp::default()).unwrap();
    assert!(
        found.free_directions.is_empty(),
        "{:?}",
        found.free_directions
    );
}

/// A sheet `z = κ x²` over 40 by 30 mm, its slope a ten-thousandth at its
/// sides, as a grid of quads.
fn sagging_sheet() -> Triangulation {
    let kappa = 1e-4 / 40.0;
    let mut mesh = Triangulation::new();
    let (nx, ny) = (80_u32, 60_u32);
    for i in 0..=nx {
        for j in 0..=ny {
            let x = -20.0 + 40.0 * f64::from(i) / f64::from(nx);
            let y = -15.0 + 30.0 * f64::from(j) / f64::from(ny);
            mesh.positions.push(Point::new(x, y, kappa * x * x));
        }
    }
    let at = |i: u32, j: u32| i * (ny + 1) + j;
    for i in 0..nx {
        for j in 0..ny {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            mesh.triangles.extend([[a, b, c], [a, c, d]]);
        }
    }
    mesh
}

/// A selection on a sheet that sags by a ten-thousandth of a radian at its
/// sides. Noise-free, two motions are free: the slide along its straight
/// lines, and the roll about its centre of curvature, which over the sheet
/// a circle matches to far below rounding; the spin is held, the sag
/// moving the residuals under it by some 1e-4 of the motion. With 0.02 mm
/// of noise the spin is free too: over the selection's radius it moves the
/// residuals by less than the noise does.
#[test]
fn what_the_noise_hides_is_free() {
    let mesh = sagging_sheet();
    let index = Bvh::build(&mesh).unwrap();
    let mut rng = Seeded(0x5a6);
    let lifted = Transform::translation(Vector::new(0.0, 0.0, 0.3));
    for (sigma, free) in [(0.0, 2), (0.02, 3)] {
        let (points, _) = scan(&mesh, 4_000, sigma, &mut rng);
        let source: Vec<Point> = points.iter().map(|&p| lifted.apply(p)).collect();
        let found = icp(&source, None, &index, &Icp::default()).unwrap();
        assert_eq!(
            found.free_directions.len(),
            free,
            "σ {sigma}: {:?}",
            found.free_directions
        );
        let slide = found.free_directions[0];
        assert!(
            slide.rotation.magnitude() * slide.radius <= 1e-6,
            "{slide:?}"
        );
        assert!(slide.translation.y.abs() >= 1.0 - 1e-6, "{slide:?}");
        let roll = found.free_directions[1];
        assert!(roll.translation.x.abs() >= 1.0 - 1e-6, "{roll:?}");
    }
}

/// The plane of a blob that is its own mirror image, placed at a tilt, is
/// found from seeds 10° and 3 mm off it to within the noise, and a seed
/// at right angles to it scores worse.
#[test]
fn a_symmetric_part_s_plane_is_found_from_a_seed_off_it() {
    let mesh = mirrored(5);
    let axis = Axis::new(
        Point::new(4.0, 1.0, -2.0),
        Direction::new(Vector::new(0.3, -0.5, 1.0), T).unwrap(),
    );
    let placed =
        Transform::translation(Vector::new(10.0, -5.0, 3.0)) * Transform::rotation(axis, 0.7);
    let mut moved = mesh.clone();
    for p in &mut moved.positions {
        *p = placed.apply(*p);
    }
    let index = Bvh::build(&moved).unwrap();
    let truth_normal = placed.apply_vector(Vector::new(1.0, 0.0, 0.0));
    let truth_origin = placed.apply(Point::ORIGIN);
    let mut rng = Seeded(0x5e_ed);
    for sigma in [0.0, 0.02] {
        let (points, normals) = scan(&moved, 20_000, sigma, &mut rng);
        let off = Transform::rotation(
            Axis::new(
                truth_origin,
                Direction::new(truth_normal.cross(rng.direction()), T).unwrap(),
            ),
            10_f64.to_radians(),
        );
        let seed = Plane::through(
            truth_origin + rng.direction() * 3.0,
            Direction::new(off.apply_vector(truth_normal), T).unwrap(),
        );
        let across = Plane::through(
            truth_origin,
            Direction::new(placed.apply_vector(Vector::new(0.0, 1.0, 0.0)), T).unwrap(),
        );
        let fit =
            symmetry_plane(&points, &normals, &index, &[across, seed], &Icp::default()).unwrap();
        assert!(fit.converged, "{fit:?}");
        let n = fit.plane.frame().z().vector();
        let angle = n
            .cross(truth_normal)
            .magnitude()
            .atan2(n.dot(truth_normal).abs());
        let offset = (fit.plane.frame().origin() - truth_origin)
            .dot(truth_normal)
            .abs();
        // Within the noise: the plane tilted across the part's half length
        // and moved along its normal by less than the scan's deviation.
        let tolerance = sigma.max(1e-6);
        println!(
            "σ {sigma}: {angle:.2e} rad, {offset:.2e} mm, rms {:.4}",
            fit.rms
        );
        assert!(
            angle * 50.0 <= tolerance && offset <= tolerance,
            "σ {sigma}: {angle:e} rad, {offset:e} mm, {fit:?}"
        );
        assert!(fit.rms <= 2.0 * sigma + 1e-6, "{fit:?}");
        assert_eq!(fit.sampled, 20_000);
        let alone = symmetry_plane(&points, &normals, &index, &[across], &Icp::default()).unwrap();
        assert!(alone.score > fit.score, "{alone:?}");
    }
}

/// The options are checked, and a facing limit needs normals.
#[test]
fn bad_options_are_refused() {
    let mesh = lopsided(2);
    let index = Bvh::build(&mesh).unwrap();
    let points = vec![Point::new(50.0, 0.0, 0.0); 3];
    for opts in [
        Icp {
            cutoff: Cutoff::Median(0.0),
            ..Icp::default()
        },
        Icp {
            max_distance: Some(0.0),
            ..Icp::default()
        },
        Icp {
            facing_limit: Some(4.0),
            ..Icp::default()
        },
        Icp {
            starts: Starts::Given(Vec::new()),
            ..Icp::default()
        },
        Icp {
            starts: Starts::Given(vec![Transform::plane_mirror(
                Point::ORIGIN,
                Direction::new(Vector::new(1.0, 0.0, 0.0), T).unwrap(),
            )]),
            ..Icp::default()
        },
        Icp {
            facing_limit: Some(1.0),
            ..Icp::default()
        },
    ] {
        assert!(icp(&points, None, &index, &opts).is_err(), "{opts:?}");
    }
    assert!(icp(&[], None, &index, &Icp::default()).is_err());
    assert!(
        icp(
            &points,
            Some(&[Vector::new(1.0, 0.0, 0.0)]),
            &index,
            &Icp::default()
        )
        .is_err()
    );
    let normals = vec![Vector::new(1.0, 0.0, 0.0); 3];
    assert!(symmetry_plane(&points, &normals, &index, &[], &Icp::default()).is_err());
    let empty = Bvh::build(&Triangulation::new()).unwrap();
    assert!(icp(&points, None, &empty, &Icp::default()).is_err());
}

/// Registration onto a rounded block's mesh of 1.3 M triangles, from the
/// identity with 10k and 100k samples and from the 24 principal starts
/// with 10k, at one thread and at eight. Run in release with `--ignored
/// --nocapture`.
#[test]
#[ignore = "a timing, for release builds"]
fn throughput_on_a_million_triangles() {
    let mesh = block(8);
    let started = ogeom_core::clock::Instant::now();
    let index = Bvh::build(&mesh).unwrap();
    println!(
        "{} triangles, hierarchy built in {:.3} s",
        mesh.triangles.len(),
        started.elapsed().as_secs_f64()
    );
    let mut rng = Seeded(0x7e57);
    let reach = 0.1 * diagonal(&mesh);
    for count in [10_000, 100_000] {
        let (points, normals) = scan(&mesh, count, 0.02, &mut rng);
        let displaced = rng.motion(Point::ORIGIN, 30_f64.to_radians(), reach);
        let source: Vec<Point> = points.iter().map(|&p| displaced.apply(p)).collect();
        let turned: Vec<Vector> = normals.iter().map(|&n| displaced.apply_vector(n)).collect();
        for (name, opts) in [
            ("identity", Icp::default()),
            (
                "identity, facing 60°",
                Icp {
                    facing_limit: Some(60_f64.to_radians()),
                    ..Icp::default()
                },
            ),
            (
                "24 principal starts",
                Icp {
                    starts: Starts::PrincipalAxes24,
                    ..Icp::default()
                },
            ),
        ] {
            if name.starts_with("24") && count > 10_000 {
                continue;
            }
            for threads in [1, 8] {
                ogeom_core::parallel::set_threads(threads);
                let started = ogeom_core::clock::Instant::now();
                let found = icp(&source, Some(&turned), &index, &opts).unwrap();
                let seconds = started.elapsed().as_secs_f64();
                let error = found.transform * displaced;
                println!(
                    "{count} samples, {name}, {threads} threads: {seconds:.3} s, back within {:.2e} rad, {:.2e} mm",
                    turn_of(&error),
                    farthest_move(&error, &points)
                );
            }
        }
    }
    ogeom_core::parallel::set_threads(0);
}
