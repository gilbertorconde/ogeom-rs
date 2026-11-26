//! Smooth free-form regions grown from clicks over meshes built in code,
//! against the faces the meshes were built from: a doubly curved lid with
//! sharp walls round it, a floor that runs up a round into a wall in three
//! patches meshed apart, and a strip folded once gently and once sharply,
//! every vertex moved by seeded Gaussian noise along each axis.
#![allow(
    clippy::unwrap_used,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    reason = "test code"
)]

use core::f64::consts::{FRAC_PI_2, TAU};
use std::collections::{BTreeMap, BTreeSet};

use ogeom_core::{OgeomError, Tolerances};
use ogeom_fit::Cutoff;
use ogeom_fit::fit::{GrowOptions, GrownRegion, grow_smooth};
use ogeom_math::{Axis, Direction, Point, Transform, Vector};
use ogeom_mesh::Bvh;
use ogeom_topo::Triangulation;

const T: Tolerances = Tolerances::millimetres();

/// The bound every test grows with: 20 degrees between smoothed normals.
fn bound() -> f64 {
    20.0_f64.to_radians()
}

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

    /// A grid of `nu + 1` by `nv + 1` vertices at `at(i, j)`, its quads
    /// split in two and wound about `outward` at their first corner.
    fn grid(
        &mut self,
        nu: usize,
        nv: usize,
        at: impl Fn(usize, usize) -> Point,
        outward: impl Fn(Point) -> Vector,
        face: u8,
    ) -> Vec<Vec<u32>> {
        let ids: Vec<Vec<u32>> = (0..=nu)
            .map(|i| (0..=nv).map(|j| self.vertex(at(i, j))).collect())
            .collect();
        self.quads(&ids, outward, face);
        ids
    }

    /// The quads between the rows of `ids`, split in two.
    fn quads(&mut self, ids: &[Vec<u32>], outward: impl Fn(Point) -> Vector, face: u8) {
        for i in 0..ids.len() - 1 {
            for j in 0..ids[i].len() - 1 {
                let [a, b, c, d] = [ids[i][j], ids[i + 1][j], ids[i + 1][j + 1], ids[i][j + 1]];
                let out = outward(self.positions[a as usize]);
                self.triangle([a, b, c], out, face);
                self.triangle([a, c, d], out, face);
            }
        }
    }

    /// The mesh placed and moved by noise of `sigma` along each axis,
    /// vertices at one point moved alike, as a scan's would be.
    fn mesh(&self, sigma: f64, seed: u64) -> Triangulation {
        let placed = placement();
        let mut draw = Seeded(seed);
        let mut moved: Vec<Point> = Vec::with_capacity(self.positions.len());
        for (i, p) in self.positions.iter().enumerate() {
            let same = (0..i).find(|&k| self.positions[k].distance(*p) <= T.confusion());
            moved.push(match same {
                Some(k) => moved[k],
                None => {
                    placed.apply(*p)
                        + Vector::new(draw.gaussian(), draw.gaussian(), draw.gaussian()) * sigma
                }
            });
        }
        Triangulation {
            positions: moved,
            triangles: self.triangles.clone(),
            ..Triangulation::new()
        }
    }

    /// The triangles of the faces named, in increasing order.
    fn of(&self, faces: &[u8]) -> Vec<u32> {
        (0..self.face.len() as u32)
            .filter(|&t| faces.contains(&self.face[t as usize]))
            .collect()
    }

    /// The largest angle, in the mesh as built, between the normals of two
    /// triangles across an edge they share, where one is of face `a` and
    /// the other of face `b`; two triangles share an edge where two of
    /// their corners lie at the same point.
    fn turn_between(&self, a: u8, b: u8) -> f64 {
        let normal = |t: usize| {
            let [p, q, r] = self.triangles[t].map(|v| self.positions[v as usize]);
            (q - p).cross(r - p).normalized(T).unwrap()
        };
        let near = |u: u32, v: u32| {
            self.positions[u as usize].distance(self.positions[v as usize]) <= T.confusion()
        };
        let mut most = 0.0_f64;
        for s in 0..self.triangles.len() {
            for t in 0..self.triangles.len() {
                if self.face[s] != a || self.face[t] != b || s == t {
                    continue;
                }
                let shared = self.triangles[s]
                    .iter()
                    .filter(|&&u| self.triangles[t].iter().any(|&v| near(u, v)))
                    .count();
                if shared == 2 {
                    let cos = normal(s).dot(normal(t)).clamp(-1.0, 1.0);
                    most = most.max(cos.acos());
                }
            }
        }
        most
    }
}

fn grown(mesh: &Triangulation, seeds: &[Point], options: &GrowOptions) -> GrownRegion {
    let index = Bvh::build(mesh).unwrap();
    let placed = placement();
    let seeds: Vec<Point> = seeds.iter().map(|p| placed.apply(*p)).collect();
    grow_smooth(mesh, &index, &seeds, bound(), options, T).unwrap()
}

/// The lid: `z = 10 - x²/80 - y²/160` over a 24 mm square in 2 mm quads,
/// face 0, with a wall 6 mm deep straight down from each border, face 1.
fn lid() -> Builder {
    let n = 12;
    let h = |x: f64, y: f64| 10.0 - x * x / 80.0 - y * y / 160.0;
    let at = |i: usize| -12.0 + 2.0 * i as f64;
    let mut b = Builder::default();
    let top = b.grid(
        n,
        n,
        |i, j| Point::new(at(i), at(j), h(at(i), at(j))),
        |_| Vector::new(0.0, 0.0, 1.0),
        0,
    );
    // The border, counter-clockwise from the corner at (-12, -12).
    let mut border: Vec<u32> = Vec::new();
    border.extend((0..n).map(|i| top[i][0]));
    border.extend((0..n).map(|j| top[n][j]));
    border.extend((0..n).map(|i| top[n - i][n]));
    border.extend((0..n).map(|j| top[0][n - j]));
    border.push(border[0]);
    let mut rows = vec![border.clone()];
    for k in 1..=3 {
        let row = border
            .iter()
            .map(|&v| {
                let p = b.positions[v as usize];
                b.vertex(Point::new(p.x, p.y, p.z - 2.0 * k as f64))
            })
            .collect::<Vec<_>>();
        rows.push(row);
    }
    // Close the loop on the first column's vertices.
    for row in &mut rows {
        let first = row[0];
        *row.last_mut().unwrap() = first;
    }
    let columns: Vec<Vec<u32>> = (0..border.len())
        .map(|c| rows.iter().map(|r| r[c]).collect())
        .collect();
    for c in 0..columns.len() - 1 {
        let p = b.positions[columns[c][0] as usize];
        let q = b.positions[columns[c + 1][0] as usize];
        let along = q - p;
        let out = Vector::new(along.y, -along.x, 0.0);
        b.quads(&[columns[c].clone(), columns[c + 1].clone()], |_| out, 1);
    }
    b
}

#[test]
fn a_doubly_curved_lid_stops_at_its_sharp_walls() {
    let b = lid();
    let expected = b.of(&[0]);
    // The lid turns at most 2.9 degrees across an edge of its own and at
    // least 70 degrees into its walls, either side of the 20 degree bound.
    let inside = b.turn_between(0, 0);
    let walls = b.turn_between(0, 1);
    assert!(inside < 3.0_f64.to_radians(), "lid turns {inside}");
    let least_wall = {
        // The smallest turn into a wall: the lid's border slopes at most
        // atan(12/40), so it meets the vertical walls at 73 degrees or more.
        FRAC_PI_2 - (12.0_f64 / 40.0).atan()
    };
    assert!(walls >= least_wall - 1e-9);
    for (sigma, seed) in [(0.0, 1), (0.1, 2), (0.1, 3), (0.1, 4)] {
        let mesh = b.mesh(sigma, seed);
        let region = grown(
            &mesh,
            &[Point::new(0.5, 0.3, 10.0)],
            &GrowOptions::default(),
        );
        assert_eq!(region.triangles, expected, "noise {sigma}, seed {seed}");
        assert_eq!(region.vertices.len(), 13 * 13);
        assert_eq!(region.points.len(), 13 * 13);
        assert_eq!(region.peeled_rings, 0);
    }
}

#[test]
fn raw_normals_alone_stop_short_on_the_noisy_lid() {
    // Without smoothing, a bound tight enough for a scan's shallow creases
    // is past the noise's turn: the lid comes back cut short.
    let b = lid();
    let mesh = b.mesh(0.1, 2);
    let raw = GrowOptions {
        smoothing_rounds: 0,
        ..GrowOptions::default()
    };
    let index = Bvh::build(&mesh).unwrap();
    let seed = placement().apply(Point::new(0.5, 0.3, 10.0));
    let tight = 12.0_f64.to_radians();
    let cut = grow_smooth(&mesh, &index, &[seed], tight, &raw, T).unwrap();
    let smoothed = grow_smooth(&mesh, &index, &[seed], tight, &GrowOptions::default(), T).unwrap();
    assert!(cut.triangles.len() < b.of(&[0]).len());
    assert_eq!(smoothed.triangles, b.of(&[0]));
}

/// A floor, face 0, `z = 0` for `x` in [-20, 0]; a round of radius 8 up
/// from it, face 1, in 12 facets; and a wall, face 2, `x = 8` for `z` in
/// [8, 20]; each 16 mm wide and meshed apart, with its own vertices along
/// the seams. A step down at the floor's far end, face 3, and a flange out
/// from the wall's top, face 4, meet them at right angles.
fn cove() -> Builder {
    let w = |j: usize| 2.0 * j as f64;
    let mut b = Builder::default();
    let floor = |i: usize, j: usize| Point::new(-20.0 + 2.0 * i as f64, w(j), 0.0);
    b.grid(10, 8, floor, |_| Vector::new(0.0, 0.0, 1.0), 0);
    let round = |i: usize, j: usize| {
        let a = FRAC_PI_2 * i as f64 / 12.0;
        Point::new(8.0 * a.sin(), w(j), 8.0 - 8.0 * a.cos())
    };
    b.grid(12, 8, round, |p| Vector::new(-p.x, 0.0, 8.0 - p.z), 1);
    let wall = |i: usize, j: usize| Point::new(8.0, w(j), 8.0 + 2.0 * i as f64);
    b.grid(6, 8, wall, |_| Vector::new(-1.0, 0.0, 0.0), 2);
    let step = |i: usize, j: usize| Point::new(-20.0, w(j), -2.0 * i as f64);
    b.grid(2, 8, step, |_| Vector::new(-1.0, 0.0, 0.0), 3);
    let flange = |i: usize, j: usize| Point::new(8.0 + 2.0 * i as f64, w(j), 20.0);
    b.grid(2, 8, flange, |_| Vector::new(0.0, 0.0, 1.0), 4);
    b
}

#[test]
fn a_region_runs_on_across_tangent_seams_between_patches_meshed_apart() {
    let b = cove();
    // The seams are tangent: the floor and the wall meet the round's first
    // and last facets at half a facet's 7.5 degrees, and the facets meet
    // each other at 7.5; the step and the flange turn 90.
    let facet = (90.0_f64 / 12.0).to_radians();
    for (a, c, expected) in [
        (0, 1, 0.5 * facet),
        (1, 1, facet),
        (1, 2, 0.5 * facet),
        (0, 3, FRAC_PI_2),
        (2, 4, FRAC_PI_2),
    ] {
        let turn = b.turn_between(a, c);
        assert!(
            (turn - expected).abs() < 1e-9,
            "faces {a} and {c} turn {turn}, not {expected}"
        );
    }
    let expected = b.of(&[0, 1, 2]);
    for (sigma, seed) in [(0.0, 5), (0.05, 6), (0.05, 7)] {
        let mesh = b.mesh(sigma, seed);
        // Seeds on the floor alone, and on the wall alone: the same region.
        for click in [Point::new(-10.1, 7.3, 0.0), Point::new(8.0, 9.1, 14.2)] {
            let region = grown(&mesh, &[click], &GrowOptions::default());
            assert_eq!(region.triangles, expected, "noise {sigma}, seed {seed}");
            // Each seam's points are named once: 11 + 13 + 7 rows of 9,
            // less the two seams' rows.
            assert_eq!(region.vertices.len(), (11 + 13 + 7 - 2) * 9);
        }
    }
}

#[test]
fn patches_whose_seam_vertices_differ_stay_apart() {
    // The round moved off the floor by ten times the confusion distance:
    // the seam's vertices are not one, so the floor grows alone.
    let mut b = cove();
    // The round's vertices follow the floor's 11 x 9.
    for p in &mut b.positions[99..99 + 13 * 9] {
        p.z += 10.0 * T.confusion();
    }
    let mesh = b.mesh(0.0, 1);
    let region = grown(
        &mesh,
        &[Point::new(-10.1, 7.3, 0.0)],
        &GrowOptions::default(),
    );
    assert_eq!(region.triangles, b.of(&[0]));
}

/// A strip 20 mm wide in 2.5 mm quads: a flat middle, face 1, for `x` in
/// [-10, 10]; a side rising 10 degrees off it, face 0, for `x` in
/// [-30, -10]; and a side rising 35 degrees off it, face 2, for `x` in
/// [10, 30]; one mesh with the folds along its grid lines.
fn folds() -> Builder {
    let (gentle, sharp) = (10.0_f64.to_radians(), 35.0_f64.to_radians());
    let mut b = Builder::default();
    let x = |i: usize| -30.0 + 2.5 * i as f64;
    b.grid(
        24,
        8,
        |i, j| {
            let x = x(i);
            let z = if x < -10.0 {
                (-10.0 - x) * gentle.tan()
            } else if x > 10.0 {
                (x - 10.0) * sharp.tan()
            } else {
                0.0
            };
            Point::new(x, 2.5 * j as f64, z)
        },
        |_| Vector::new(0.0, 0.0, 1.0),
        0,
    );
    // Name the faces by where each triangle's centroid lies.
    for t in 0..b.triangles.len() {
        let c = b.triangles[t]
            .iter()
            .map(|&v| b.positions[v as usize].x)
            .sum::<f64>()
            / 3.0;
        b.face[t] = if c < -10.0 {
            0
        } else if c > 10.0 {
            2
        } else {
            1
        };
    }
    b
}

#[test]
fn a_crease_past_the_bound_stops_the_region_and_a_gentler_fold_does_not() {
    let b = folds();
    let gentle = b.turn_between(0, 1);
    let sharp = b.turn_between(1, 2);
    assert!((gentle - 10.0_f64.to_radians()).abs() < 1e-9);
    assert!((sharp - 35.0_f64.to_radians()).abs() < 1e-9);
    assert!(b.turn_between(0, 0) < 1e-6);
    for (sigma, seed) in [(0.0, 8), (0.1, 9), (0.1, 10), (0.1, 11)] {
        let mesh = b.mesh(sigma, seed);
        let on_gentle_side = grown(
            &mesh,
            &[Point::new(-20.0, 10.0, 1.8)],
            &GrowOptions::default(),
        );
        assert_eq!(
            on_gentle_side.triangles,
            b.of(&[0, 1]),
            "noise {sigma}, seed {seed}"
        );
        let on_sharp_side = grown(
            &mesh,
            &[Point::new(20.0, 10.0, 7.0)],
            &GrowOptions::default(),
        );
        assert_eq!(
            on_sharp_side.triangles,
            b.of(&[2]),
            "noise {sigma}, seed {seed}"
        );
        // Seeds either side of the crease give the union.
        let both = grown(
            &mesh,
            &[Point::new(-20.0, 10.0, 1.8), Point::new(20.0, 10.0, 7.0)],
            &GrowOptions::default(),
        );
        assert_eq!(both.triangles, b.of(&[0, 1, 2]));
        assert_eq!(both.seeds.len(), 2);
    }
}

/// Whether the triangles in one of `a` and `b` only all lie within a ring
/// of `a`'s border: each has a corner on an edge of `a` that only one
/// triangle of `a` uses.
fn within_a_ring(mesh: &Triangulation, a: &[u32], b: &[u32]) -> bool {
    let a: BTreeSet<u32> = a.iter().copied().collect();
    let b: BTreeSet<u32> = b.iter().copied().collect();
    let corners = |t: u32| mesh.triangles[t as usize];
    let mut uses: BTreeMap<(u32, u32), usize> = BTreeMap::new();
    for &t in &a {
        let c = corners(t);
        for k in 0..3 {
            let (p, q) = (c[k], c[(k + 1) % 3]);
            *uses.entry((p.min(q), p.max(q))).or_default() += 1;
        }
    }
    let border: BTreeSet<u32> = uses
        .iter()
        .filter(|(_, n)| **n == 1)
        .flat_map(|((p, q), _)| [*p, *q])
        .collect();
    a.symmetric_difference(&b)
        .all(|&t| corners(t).iter().any(|v| border.contains(v)))
}

#[test]
fn noisy_and_clean_copies_grow_the_same_region_within_a_ring() {
    // The folds under noise of 0.1 mm, which turns single edges by up to
    // 20 degrees, against the clean copy.
    let b = folds();
    let clean = b.mesh(0.0, 1);
    let seeds = [Point::new(-20.0, 10.0, 1.8)];
    let reference = grown(&clean, &seeds, &GrowOptions::default());
    assert_eq!(reference.triangles, b.of(&[0, 1]));
    for seed in 20..28 {
        let noisy = b.mesh(0.1, seed);
        let region = grown(&noisy, &seeds, &GrowOptions::default());
        assert!(
            within_a_ring(&clean, &reference.triangles, &region.triangles),
            "seed {seed}"
        );
    }
}

/// A floor, face 0, whose edge is rounded through 90 degrees over four
/// rows, turning 15, 30, 30 and 15 degrees, faces 1, 2, 2 and 3, into a
/// wall, face 4: a scanned sharp edge. Rows of 2 mm, 20 mm wide.
fn rounded_edge() -> Builder {
    let turns: [f64; 5] = [0.0, 15.0, 45.0, 75.0, 90.0];
    let mut b = Builder::default();
    let mut profile = vec![];
    for i in 0..=10 {
        profile.push(Point::new(-20.0 + 2.0 * i as f64, 0.0, 0.0));
    }
    let mut p = *profile.last().unwrap();
    for a in &turns[1..] {
        let a = a.to_radians();
        p = Point::new(p.x + 2.0 * a.cos(), 0.0, p.z + 2.0 * a.sin());
        profile.push(p);
    }
    for _ in 0..4 {
        p = Point::new(p.x, 0.0, p.z + 2.0);
        profile.push(p);
    }
    let rows = profile.len();
    b.grid(
        rows - 1,
        10,
        |i, j| Point::new(profile[i].x, 2.0 * j as f64, profile[i].z),
        |p| {
            if p.z > 0.5 {
                Vector::new(-1.0, 0.0, 1.0)
            } else {
                Vector::new(0.0, 0.0, 1.0)
            }
        },
        0,
    );
    for t in 0..b.triangles.len() {
        // Each quad's two triangles, row by row.
        let row = t / 20;
        b.face[t] = match row {
            0..=9 => 0,
            10 => 1,
            11 | 12 => 2,
            13 => 3,
            _ => 4,
        };
    }
    b
}

#[test]
fn peeling_takes_back_the_round_a_region_climbed() {
    let b = rounded_edge();
    for (sigma, seed) in [(0.0, 30), (0.05, 31), (0.05, 32)] {
        let mesh = b.mesh(sigma, seed);
        let seeds = [Point::new(-10.1, 9.3, 0.0)];
        let plain = grown(&mesh, &seeds, &GrowOptions::default());
        // The round's first row turns 15 degrees, within the bound, and
        // its second 30, past it: the region climbs one row.
        assert_eq!(plain.triangles, b.of(&[0, 1]), "noise {sigma}, seed {seed}");
        let peel = GrowOptions {
            peel: true,
            ..GrowOptions::default()
        };
        let peeled = grown(&mesh, &seeds, &peel);
        if sigma == 0.0 {
            assert_eq!(peeled.triangles, b.of(&[0]));
            assert_eq!(peeled.peeled_rings, 1);
        } else {
            // Under noise a rim triangle the noise turned as far as the
            // round goes too.
            let band = peeled
                .triangles
                .iter()
                .filter(|&&t| b.face[t as usize] != 0)
                .count();
            assert_eq!(band, 0, "seed {seed}, {} triangles", peeled.triangles.len());
            let missing: Vec<u32> = b
                .of(&[0])
                .into_iter()
                .filter(|t| peeled.triangles.binary_search(t).is_err())
                .collect();
            assert!(
                within_a_ring(&mesh, &b.of(&[0]), &peeled.triangles),
                "seed {seed}, passes {}, missing {missing:?}",
                peeled.peeled_rings
            );
        }
    }
}

#[test]
fn peeling_leaves_a_lid_with_sharp_walls_whole_within_the_noise() {
    let b = lid();
    let peel = GrowOptions {
        peel: true,
        ..GrowOptions::default()
    };
    for (sigma, seed) in [(0.0, 1), (0.1, 2), (0.1, 3)] {
        let mesh = b.mesh(sigma, seed);
        let region = grown(&mesh, &[Point::new(0.5, 0.3, 10.0)], &peel);
        if sigma == 0.0 {
            assert_eq!(region.triangles, b.of(&[0]));
            assert_eq!(region.peeled_rings, 0);
        } else {
            assert!(
                within_a_ring(&mesh, &b.of(&[0]), &region.triangles),
                "seed {seed}"
            );
        }
    }
}

#[test]
fn the_cap_holds_growth_to_its_count() {
    let b = lid();
    let mesh = b.mesh(0.0, 1);
    let capped = GrowOptions {
        max_triangles: Some(40),
        ..GrowOptions::default()
    };
    let region = grown(&mesh, &[Point::new(0.5, 0.3, 10.0)], &capped);
    assert_eq!(region.triangles.len(), 40);
    assert!(region.triangles.iter().all(|&t| b.face[t as usize] == 0));
}

#[test]
fn options_out_of_range_are_refused() {
    let b = lid();
    let mesh = b.mesh(0.0, 1);
    let index = Bvh::build(&mesh).unwrap();
    let seed = [placement().apply(Point::new(0.0, 0.0, 10.0))];
    let o = GrowOptions::default();
    let domain = |r: Result<GrownRegion, OgeomError>| matches!(r, Err(OgeomError::Domain(_)));
    assert!(domain(grow_smooth(&mesh, &index, &[], bound(), &o, T)));
    assert!(domain(grow_smooth(
        &mesh,
        &index,
        &[Point::new(f64::NAN, 0.0, 0.0)],
        bound(),
        &o,
        T
    )));
    assert!(domain(grow_smooth(&mesh, &index, &seed, 0.0, &o, T)));
    assert!(domain(grow_smooth(
        &mesh,
        &index,
        &seed,
        core::f64::consts::PI,
        &o,
        T
    )));
    for cutoff in [Cutoff::All, Cutoff::Distance(0.1)] {
        let peel = GrowOptions {
            peel: true,
            cutoff,
            ..o
        };
        assert!(domain(grow_smooth(&mesh, &index, &seed, bound(), &peel, T)));
    }
    let none = GrowOptions {
        max_triangles: Some(0),
        ..o
    };
    assert!(domain(grow_smooth(&mesh, &index, &seed, bound(), &none, T)));
    let mut fewer = mesh.clone();
    fewer.triangles.pop();
    let short = Bvh::build(&fewer).unwrap();
    assert!(matches!(
        grow_smooth(&mesh, &short, &seed, bound(), &o, T),
        Err(OgeomError::Dimension(_))
    ));
}
