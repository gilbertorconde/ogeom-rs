//! Smoothing an editable mesh onto local quadrics, and the curvature those
//! quadrics read, against shapes with closed forms: a sphere, a cylinder,
//! a torus and a box.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{FRAC_PI_2, FRAC_PI_6, TAU};
use std::time::Duration;

#[path = "support/cpu_time.rs"]
mod cpu_time;
use cpu_time::cpu_time;

use ogeom::core::FastMap;
use ogeom::math::{Point, Vector};
use ogeom::mesh::{Curvature, Editable, Smooth, curvature};
use ogeom::topo::Triangulation;

/// Gaussian samples from a fixed seed: splitmix64 through Box–Muller.
struct Noise(u64);

impl Noise {
    #[allow(clippy::cast_precision_loss, reason = "53 bits convert exactly")]
    fn uniform(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        ((z >> 11) as f64 + 0.5) / (1_u64 << 53) as f64
    }

    fn gaussian(&mut self, sigma: f64) -> f64 {
        let (a, b) = (self.uniform(), self.uniform());
        sigma * (-2.0 * a.ln()).sqrt() * (TAU * b).cos()
    }

    /// Every vertex moved by an independent Gaussian of `sigma` along
    /// each axis.
    fn shake(&mut self, mesh: &mut Triangulation, sigma: f64) {
        for p in &mut mesh.positions {
            *p += Vector::new(
                self.gaussian(sigma),
                self.gaussian(sigma),
                self.gaussian(sigma),
            );
        }
    }
}

/// A sphere about the origin: an icosahedron split `level` times, every
/// vertex on the sphere, wound outward.
fn sphere(radius: f64, level: u32) -> Triangulation {
    let t = f64::midpoint(1.0, 5.0_f64.sqrt());
    let mut points: Vec<Vector> = [
        (-1.0, t, 0.0),
        (1.0, t, 0.0),
        (-1.0, -t, 0.0),
        (1.0, -t, 0.0),
        (0.0, -1.0, t),
        (0.0, 1.0, t),
        (0.0, -1.0, -t),
        (0.0, 1.0, -t),
        (t, 0.0, -1.0),
        (t, 0.0, 1.0),
        (-t, 0.0, -1.0),
        (-t, 0.0, 1.0),
    ]
    .iter()
    .map(|&(x, y, z)| {
        let v = Vector::new(x, y, z);
        v / v.magnitude()
    })
    .collect();
    let mut faces: Vec<[u32; 3]> = vec![
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
    for _ in 0..level {
        let mut middle: FastMap<(u32, u32), u32> = FastMap::default();
        let mut split = |a: u32, b: u32, points: &mut Vec<Vector>| {
            *middle.entry((a.min(b), a.max(b))).or_insert_with(|| {
                let m = points[a as usize] + points[b as usize];
                points.push(m / m.magnitude());
                u32::try_from(points.len() - 1).unwrap()
            })
        };
        let mut finer = Vec::with_capacity(faces.len() * 4);
        for [a, b, c] in faces {
            let ab = split(a, b, &mut points);
            let bc = split(b, c, &mut points);
            let ca = split(c, a, &mut points);
            finer.extend([[a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]]);
        }
        faces = finer;
    }
    Triangulation {
        positions: points
            .iter()
            .map(|v| Point::new(v.x * radius, v.y * radius, v.z * radius))
            .collect(),
        triangles: faces,
        ..Triangulation::new()
    }
}

/// An open cylinder about the `z` axis from `0` to `height`, `around` by
/// `along` quads, two triangles each, wound outward.
fn cylinder(radius: f64, height: f64, around: u32, along: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    for j in 0..=along {
        let z = height * f64::from(j) / f64::from(along);
        for i in 0..around {
            let u = TAU * f64::from(i) / f64::from(around);
            mesh.positions
                .push(Point::new(radius * u.cos(), radius * u.sin(), z));
        }
    }
    let at = |i: u32, j: u32| j * around + i % around;
    for j in 0..along {
        for i in 0..around {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            mesh.triangles.push([a, b, c]);
            mesh.triangles.push([a, c, d]);
        }
    }
    mesh
}

/// A torus about the `z` axis, `rings` by `sides` quads, two triangles
/// each, wound outward; the vertex `(i, j)` at angle `TAU i / rings` about
/// the axis and `TAU j / sides` about the tube from its outer equator.
fn torus(major: f64, minor: f64, rings: u32, sides: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    for i in 0..rings {
        let u = f64::from(i) / f64::from(rings) * TAU;
        for j in 0..sides {
            let v = f64::from(j) / f64::from(sides) * TAU;
            let r = major + minor * v.cos();
            mesh.positions
                .push(Point::new(r * u.cos(), r * u.sin(), minor * v.sin()));
        }
    }
    let at = |i: u32, j: u32| (i % rings) * sides + (j % sides);
    for i in 0..rings {
        for j in 0..sides {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            mesh.triangles.push([a, b, c]);
            mesh.triangles.push([a, c, d]);
        }
    }
    mesh
}

/// The cube `[-half, half]³`, each face a grid of `cells` by `cells`
/// quads, welded into one closed mesh, wound outward.
fn cube(half: f64, cells: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    for axis in 0..3 {
        for sign in [-1.0, 1.0] {
            let first = u32::try_from(mesh.positions.len()).unwrap();
            for j in 0..=cells {
                for i in 0..=cells {
                    let s = -half + 2.0 * half * f64::from(i) / f64::from(cells);
                    let t = -half + 2.0 * half * f64::from(j) / f64::from(cells);
                    let mut c = [0.0; 3];
                    c[axis] = sign * half;
                    c[(axis + 1) % 3] = s;
                    c[(axis + 2) % 3] = t;
                    mesh.positions.push(Point::new(c[0], c[1], c[2]));
                }
            }
            let at = |i: u32, j: u32| first + j * (cells + 1) + i;
            for j in 0..cells {
                for i in 0..cells {
                    let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
                    // The grid's axes, crossed, point along `+axis`.
                    if sign > 0.0 {
                        mesh.triangles.extend([[a, b, c], [a, c, d]]);
                    } else {
                        mesh.triangles.extend([[a, c, b], [a, d, c]]);
                    }
                }
            }
        }
    }
    Editable::from(&mesh, half * 1e-9)
        .unwrap()
        .to_triangulation()
}

/// The root mean square of `f` over `items`.
fn rms<T>(items: &[T], f: impl Fn(&T) -> f64) -> f64 {
    (items.iter().map(|x| f(x).powi(2)).sum::<f64>()
        / f64::from(u32::try_from(items.len()).unwrap()))
    .sqrt()
}

/// A ball of radius 10 under Gaussian noise of 0.02 per axis, smoothed
/// with neighbourhoods of 1: its radius, fitted about the known centre,
/// stays within a micrometre, and the spread of the vertices' radii falls
/// from the noise's to below a quarter of it.
#[test]
fn a_noisy_sphere_smooths_onto_its_radius() {
    let (radius, sigma) = (10.0, 0.02);
    let mut mesh = sphere(radius, 6);
    Noise(1).shake(&mut mesh, sigma);
    let error = |m: &Triangulation| {
        let r: Vec<f64> = m
            .positions
            .iter()
            .map(|p| (*p - Point::new(0.0, 0.0, 0.0)).magnitude())
            .collect();
        let fitted = r.iter().sum::<f64>() / f64::from(u32::try_from(r.len()).unwrap());
        (fitted - radius, rms(&r, |x| x - fitted))
    };
    let (bias, before) = error(&mesh);
    let mut e = Editable::from(&mesh, 0.0).unwrap();
    let report = e
        .smooth(&Smooth {
            size: 1.0,
            max_move: 10.0 * sigma,
            keep_edges: None,
            region: None,
        })
        .unwrap();
    let smoothed = e.to_triangulation();
    let (after_bias, after) = error(&smoothed);
    eprintln!(
        "sphere: radius off by {bias:.2e} then {after_bias:.2e}; spread {before:.5} then \
         {after:.5}; {report:?}"
    );
    assert!(before > 0.9 * sigma, "the noise is there: {before}");
    assert!(after_bias.abs() < 1e-3, "radius off by {after_bias}");
    assert!(after < sigma / 4.0, "spread {after} after smoothing");
    assert_eq!(report.limited, 0);
    assert_eq!(smoothed.triangles, mesh.triangles);
}

/// No vertex ends farther than `max_move` from where it started, the
/// report's largest move is the largest there is, and every vertex held
/// back ends exactly that far, as computed.
#[test]
fn max_move_holds() {
    let mut mesh = sphere(10.0, 4);
    Noise(2).shake(&mut mesh, 0.05);
    let mut e = Editable::from(&mesh, 0.0).unwrap();
    let max_move = 0.01;
    let report = e
        .smooth(&Smooth {
            size: 2.0,
            max_move,
            keep_edges: None,
            region: None,
        })
        .unwrap();
    let moved: Vec<f64> = mesh
        .positions
        .iter()
        .zip(&e.to_triangulation().positions)
        .map(|(a, b)| (*b - *a).magnitude())
        .collect();
    let largest = moved.iter().copied().fold(0.0, f64::max);
    eprintln!("max_move {max_move}: {report:?}");
    assert!(moved.iter().all(|&d| d <= max_move));
    assert_eq!(report.max, largest);
    assert!(report.limited > moved.len() / 2, "{report:?}");
    let at_cap = moved
        .iter()
        .filter(|&&d| d > max_move * (1.0 - 1e-12))
        .count();
    assert!(at_cap >= report.limited, "{at_cap} at the cap");

    let mut still = Editable::from(&mesh, 0.0).unwrap();
    let report = still
        .smooth(&Smooth {
            size: 2.0,
            max_move: 0.0,
            keep_edges: None,
            region: None,
        })
        .unwrap();
    assert_eq!(report.max, 0.0);
    assert_eq!(still.to_triangulation().positions, mesh.positions);
}

/// The dihedral angle, for each of a cube's twelve edges, between the
/// planes fitted to its two faces over the band within `band` of it; and the root mean square distance from the edge's
/// line of the vertices that started on it, corners left out.
fn cube_edges(start: &Triangulation, now: &Triangulation, half: f64, band: f64) -> (Vec<f64>, f64) {
    let on = |p: Point, axis: usize, sign: f64| (p.coord(axis).unwrap() - sign * half).abs() < 1e-9;
    let mut dihedrals = Vec::new();
    let mut sag = Vec::new();
    for a in 0..3 {
        for b in a + 1..3 {
            for sa in [-1.0, 1.0] {
                for sb in [-1.0, 1.0] {
                    // The plane fitted to the vertices of face `axis = sign`
                    // within `band` of the edge, as the axis's coordinate
                    // over the other two; its outward normal.
                    let plane = |axis: usize, sign: f64, other: usize, other_sign: f64| {
                        let (s, t) = ((axis + 1) % 3, (axis + 2) % 3);
                        let mut m = [[0.0; 3]; 3];
                        let mut r = [0.0; 3];
                        for (p0, p) in start.positions.iter().zip(&now.positions) {
                            let near = (p0.coord(other).unwrap() - other_sign * half).abs();
                            if on(*p0, axis, sign) && !on(*p0, other, other_sign) && near <= band {
                                let row = [1.0, p.coord(s).unwrap(), p.coord(t).unwrap()];
                                for i in 0..3 {
                                    r[i] += row[i] * p.coord(axis).unwrap();
                                    for j in 0..3 {
                                        m[i][j] += row[i] * row[j];
                                    }
                                }
                            }
                        }
                        let c = solve3(m, r);
                        let mut n = [0.0; 3];
                        n[axis] = sign;
                        n[s] = -sign * c[1];
                        n[t] = -sign * c[2];
                        let n = Vector::new(n[0], n[1], n[2]);
                        n / n.magnitude()
                    };
                    let na = plane(a, sa, b, sb);
                    let nb = plane(b, sb, a, sa);
                    dihedrals.push(std::f64::consts::PI - na.dot(nb).clamp(-1.0, 1.0).acos());
                    let c = 3 - a - b;
                    for (p0, p) in start.positions.iter().zip(&now.positions) {
                        if on(*p0, a, sa)
                            && on(*p0, b, sb)
                            && p0.coord(c).unwrap().abs() < half - 1e-9
                        {
                            let da = p.coord(a).unwrap() - sa * half;
                            let db = p.coord(b).unwrap() - sb * half;
                            sag.push(da.hypot(db));
                        }
                    }
                }
            }
        }
    }
    (dihedrals, rms(&sag, |&d| d))
}

/// Solve a 3 by 3 system by Cramer's rule.
fn solve3(m: [[f64; 3]; 3], r: [f64; 3]) -> [f64; 3] {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(m);
    let mut x = [0.0; 3];
    for (k, x) in x.iter_mut().enumerate() {
        let mut mk = m;
        for i in 0..3 {
            mk[i][k] = r[i];
        }
        *x = det(mk) / d;
    }
    x
}

/// A cube under noise, smoothed with `keep_edges`, keeps its edges: the
/// faces beside each stay at a right angle and the vertices on it stay
/// within the noise of it, while the faces flatten. Smoothed without, the
/// same cube's edges round over.
#[test]
fn a_noisy_cube_keeps_its_edges() {
    let (half, sigma, size) = (5.0, 0.02, 1.5);
    let clean = cube(half, 40);
    let mut noisy = clean.clone();
    Noise(3).shake(&mut noisy, sigma);
    let flatness = |m: &Triangulation| {
        let off: Vec<f64> = clean
            .positions
            .iter()
            .zip(&m.positions)
            .filter(|(p0, _)| {
                [p0.x, p0.y, p0.z]
                    .iter()
                    .filter(|c| (c.abs() - half).abs() < 1e-9)
                    .count()
                    == 1
            })
            .map(|(p0, p)| {
                let axis = (0..3)
                    .find(|&k| (p0.coord(k).unwrap().abs() - half).abs() < 1e-9)
                    .unwrap();
                p.coord(axis).unwrap() - p0.coord(axis).unwrap()
            })
            .collect();
        rms(&off, |&d| d)
    };
    let smoothed = |keep_edges| {
        let mut e = Editable::from(&noisy, 0.0).unwrap();
        let report = e
            .smooth(&Smooth {
                size,
                max_move: 1.0,
                keep_edges,
                region: None,
            })
            .unwrap();
        (e.to_triangulation(), report)
    };
    let (kept, kept_report) = smoothed(Some(FRAC_PI_6));
    let (rounded, rounded_report) = smoothed(None);
    let band = size / 2.0;
    let (d0, sag0) = cube_edges(&clean, &noisy, half, band);
    let (dk, sagk) = cube_edges(&clean, &kept, half, band);
    let (dr, sagr) = cube_edges(&clean, &rounded, half, band);
    let worst = |d: &[f64]| {
        d.iter()
            .map(|a| (a - FRAC_PI_2).abs().to_degrees())
            .fold(0.0, f64::max)
    };
    eprintln!(
        "cube: dihedral off by at most {:.3}° noisy, {:.3}° kept, {:.3}° rounded; edge sag {sag0:.4} \
         noisy, {sagk:.4} kept, {sagr:.4} rounded; faces off their planes {:.4} noisy, {:.4} kept; \
         {kept_report:?} {rounded_report:?}",
        worst(&d0),
        worst(&dk),
        worst(&dr),
        flatness(&noisy),
        flatness(&kept),
    );
    assert!(worst(&dk) < 1.0, "kept edges bend by {}°", worst(&dk));
    assert!(sagk < 2.0 * sigma, "kept edges sag by {sagk}");
    assert!(flatness(&kept) < sigma / 2.0);
    assert!(worst(&dr) > 5.0, "rounded edges bend by {}°", worst(&dr));
    assert!(sagr > 5.0 * sigma, "rounded edges sag by {sagr}");
}

/// The curvature of a sphere of radius `R` is `1/R` both ways at every
/// vertex, to within the quadric's truncation: over a neighbourhood of
/// radius `r` the sphere's height departs from its osculating paraboloid
/// by `ρ⁴/8R³`, which a least-squares quadric over `ρ ≤ r` absorbs as a
/// curvature off by no more than `r²/2R²` of itself.
#[test]
fn curvature_of_a_sphere() {
    let radius = 10.0;
    let size = 1.0;
    let mesh = sphere(radius, 6);
    let found = curvature(&mesh, size).unwrap();
    let bound = (size / radius).powi(2) / 2.0;
    let mut worst = 0.0_f64;
    for c in &found {
        let c = c.unwrap();
        for k in [c.k1, c.k2, c.mean] {
            worst = worst.max((k * radius - 1.0).abs());
        }
        worst = worst.max(((c.gaussian * radius * radius).sqrt() - 1.0).abs());
    }
    eprintln!("sphere curvature off by at most {worst:.2e} of itself, bound {bound:.2e}");
    assert!(worst < bound);
}

/// A cylinder of radius `r` bends by `1/r` about its axis and not at all
/// along it, the directions about and along it, at every vertex at least
/// `size` from its ends, to the truncation bound of the sphere's.
#[test]
fn curvature_of_a_cylinder() {
    let (radius, height, size) = (5.0, 10.0, 1.0);
    let mesh = cylinder(radius, height, 180, 60);
    let found = curvature(&mesh, size).unwrap();
    let bound = (size / radius).powi(2) / 2.0;
    let (mut k1, mut k2, mut dirs) = (0.0_f64, 0.0_f64, 0.0_f64);
    let axis = Vector::new(0.0, 0.0, 1.0);
    for (p, c) in mesh.positions.iter().zip(&found) {
        if p.z < size || p.z > height - size {
            continue;
        }
        let Curvature {
            k1: a,
            k2: b,
            dir1,
            dir2,
            mean,
            gaussian,
        } = c.unwrap();
        k1 = k1.max((a * radius - 1.0).abs());
        k2 = k2.max((b * radius).abs());
        dirs = dirs
            .max(dir1.dot(axis).abs())
            .max(1.0 - dir2.dot(axis).abs());
        assert!((mean - (a + b) / 2.0).abs() < 1e-15 && (gaussian - a * b).abs() < 1e-15);
    }
    eprintln!(
        "cylinder: k1 off by {k1:.2e}, k2 by {k2:.2e}, directions by {dirs:.2e}; bound {bound:.2e}"
    );
    assert!(k1 < bound && k2 < bound && dirs < bound);
}

/// On a torus of radii `R` and `r`, at angle `v` about the tube from its
/// outer equator, the tube bends by `1/r` and the ring by
/// `cos v / (R + r cos v)`, the tube's direction across the ring's. Each
/// is read within the truncation bound for the largest curvature.
#[test]
fn curvature_of_a_torus() {
    let (major, minor, size) = (10.0, 3.0, 0.5);
    let (rings, sides) = (240, 72);
    let mesh = torus(major, minor, rings, sides);
    let found = curvature(&mesh, size).unwrap();
    let largest = 1.0 / minor;
    let bound = largest * (size * largest).powi(2) / 2.0;
    let mut worst = 0.0_f64;
    for j in 0..sides {
        let v = TAU * f64::from(j) / f64::from(sides);
        let tube = 1.0 / minor;
        let ring = v.cos() / (major + minor * v.cos());
        let c = found[j as usize].unwrap();
        // At `u = 0` the ring runs along `y`.
        let along_ring = Vector::new(0.0, 1.0, 0.0);
        worst = worst
            .max((c.k1 - tube).abs())
            .max((c.k2 - ring).abs())
            .max(largest * c.dir1.dot(along_ring).abs())
            .max(largest * (1.0 - c.dir2.dot(along_ring).abs()));
    }
    eprintln!("torus curvature off by at most {worst:.2e}, bound {bound:.2e}");
    assert!(worst < bound);
}

/// Only the region moves, blending in: the far side of a sphere outside it
/// is left to the bit, its border moves by about its share of a full move,
/// and its middle, far enough in that every neighbour it reads was moved in
/// full, is smoothed in full.
#[test]
fn region_blends_in() {
    let (radius, sigma, size) = (10.0, 0.02, 1.0);
    let mut mesh = sphere(radius, 6);
    Noise(4).shake(&mut mesh, sigma);
    let region: Vec<bool> = mesh.positions.iter().map(|p| p.z > 0.0).collect();
    let mut e = Editable::from(&mesh, 0.0).unwrap();
    e.smooth(&Smooth {
        size,
        max_move: 1.0,
        keep_edges: None,
        region: Some(&region),
    })
    .unwrap();
    let mut whole = Editable::from(&mesh, 0.0).unwrap();
    whole
        .smooth(&Smooth {
            size,
            max_move: 1.0,
            keep_edges: None,
            region: None,
        })
        .unwrap();
    let (part, full) = (e.to_triangulation(), whole.to_triangulation());
    let outside: Vec<Point> = mesh
        .positions
        .iter()
        .copied()
        .filter(|p| p.z <= 0.0)
        .collect();
    let (mut ours, mut theirs, mut share) = (Vec::new(), Vec::new(), 0.0_f64);
    for (k, p) in mesh.positions.iter().enumerate() {
        if p.z <= 0.0 {
            assert_eq!(part.positions[k], *p, "vertex {k} outside the region moved");
        } else if p.z > 3.0 * size {
            assert_eq!(
                part.positions[k], full.positions[k],
                "vertex {k} not smoothed in full"
            );
        } else if p.z < 0.05 * size {
            let gap = outside
                .iter()
                .map(|q| (*q - *p).magnitude())
                .fold(f64::MAX, f64::min);
            share = share.max(gap / size);
            ours.push((part.positions[k] - *p).magnitude());
            theirs.push((full.positions[k] - *p).magnitude());
        }
    }
    let ratio = rms(&ours, |&d| d) / rms(&theirs, |&d| d);
    eprintln!(
        "region: {} vertices at its border moved {ratio:.3} as far as smoothed in full, each \
         pass's move at most {share:.3} of theirs",
        ours.len()
    );
    assert!(!ours.is_empty() && ratio < 2.0 * share && share < 0.5);
}

/// Smoothing and curvature give the same answers to the bit on one thread
/// as on several.
#[test]
fn the_same_at_any_thread_count() {
    let mut mesh = sphere(10.0, 4);
    Noise(5).shake(&mut mesh, 0.05);
    let run = |threads: usize| {
        ogeom::core::parallel::set_threads(threads);
        let mut e = Editable::from(&mesh, 0.0).unwrap();
        let report = e
            .smooth(&Smooth {
                size: 2.0,
                max_move: 0.1,
                keep_edges: Some(FRAC_PI_6),
                region: None,
            })
            .unwrap();
        let read = curvature(&mesh, 2.0).unwrap();
        ogeom::core::parallel::set_threads(0);
        (e.to_triangulation().positions, report, read)
    };
    assert_eq!(run(1), run(4));
}

/// Two cubes sharing an edge, four faces on it, smoothed under noise:
/// the vertices on the shared edge stay where they are, each cube still
/// closes on its own, and the joins are the ones a mesh built afresh from
/// the result makes.
#[test]
fn an_edge_of_four_faces_is_joined_as_it_stands() {
    let half = 5.0;
    let a = cube(half, 10);
    let mut both = a.clone();
    let shift = Vector::new(2.0 * half, 2.0 * half, 0.0);
    let count = u32::try_from(a.positions.len()).unwrap();
    both.positions
        .extend(a.positions.iter().map(|&p| p + shift));
    both.triangles
        .extend(a.triangles.iter().map(|t| t.map(|v| v + count)));
    // Noise drawn from each position, so the two cubes' copies of a point
    // on the shared edge move alike and weld again.
    for p in &mut both.positions {
        let mut noise =
            Noise(p.x.to_bits() ^ p.y.to_bits().rotate_left(21) ^ p.z.to_bits().rotate_left(42));
        *p += Vector::new(
            noise.gaussian(0.02),
            noise.gaussian(0.02),
            noise.gaussian(0.02),
        );
    }
    let before = Editable::from(&both, half * 1e-9).unwrap();
    let mut e = Editable::from(&both, half * 1e-9).unwrap();
    let pinned = e.non_manifold();
    assert!(!pinned.is_empty());
    e.smooth(&Smooth {
        size: 3.0,
        max_move: 1.0,
        keep_edges: Some(FRAC_PI_6),
        region: None,
    })
    .unwrap();
    for &v in &pinned {
        assert_eq!(e.position(v), before.position(v));
    }
    let afresh = Editable::from(&e.to_triangulation(), half * 1e-9).unwrap();
    assert!(e.is_closed());
    assert_eq!(e.components().len(), 2);
    assert_eq!(e.non_manifold(), afresh.non_manifold());
    assert_eq!(e.to_triangulation(), afresh.to_triangulation());
}

/// A size that is not a length, a negative `max_move`, a crease angle past
/// `π` and a region of the wrong length are refused, the mesh untouched.
#[test]
fn bad_settings_are_refused() {
    let mesh = sphere(1.0, 2);
    let mut e = Editable::from(&mesh, 0.0).unwrap();
    let fine = Smooth {
        size: 0.5,
        max_move: 0.1,
        keep_edges: None,
        region: None,
    };
    let short = vec![true; 3];
    for bad in [
        Smooth { size: 0.0, ..fine },
        Smooth {
            size: f64::NAN,
            ..fine
        },
        Smooth {
            max_move: -1.0,
            ..fine
        },
        Smooth {
            keep_edges: Some(4.0),
            ..fine
        },
        Smooth {
            region: Some(&short),
            ..fine
        },
    ] {
        assert!(e.smooth(&bad).is_err(), "{bad:?}");
    }
    assert_eq!(e.to_triangulation().positions, mesh.positions);
    assert!(curvature(&mesh, -1.0).is_err());
}

/// Five million triangles of a noisy torus smoothed, and their curvature
/// read, each within a bound on processor time; the smoothing takes the
/// noise down to below a quarter.
#[test]
#[ignore = "heavy"]
fn five_million_triangles_smoothed() {
    let (major, minor, sigma) = (10.0, 3.0, 0.005);
    let mut mesh = torus(major, minor, 2500, 1000);
    Noise(7).shake(&mut mesh, sigma);
    let off = |m: &Triangulation| {
        rms(&m.positions, |p| {
            (Vector::new(p.x, p.y, 0.0).magnitude() - major).hypot(p.z) - minor
        })
    };
    let mut e = Editable::from(&mesh, 0.0).unwrap();
    let (report, smoothing) = cpu_time(|| {
        e.smooth(&Smooth {
            size: 0.15,
            max_move: 0.05,
            keep_edges: None,
            region: None,
        })
        .unwrap()
    });
    let smoothed = e.to_triangulation();
    let (read, reading) = cpu_time(|| curvature(&smoothed, 0.15).unwrap());
    let (before, after) = (off(&mesh), off(&smoothed));
    eprintln!(
        "5M triangles: smoothed in {smoothing:?}, off the torus by {before:.5} then {after:.5}, \
         {report:?}; curvature read in {reading:?}"
    );
    assert!(after < before / 4.0);
    assert!(read.iter().all(Option::is_some));
    assert!(
        smoothing < Duration::from_secs(600),
        "{smoothing:?} to smooth"
    );
    assert!(
        reading < Duration::from_secs(300),
        "{reading:?} to read curvature"
    );
}
