//! Offsetting an editable mesh against shapes with closed forms: boxes and
//! a notched block, whose offsets with sharp edges are boxes and blocks
//! again, a sphere, whose offset is a sphere, and a plate with a region
//! stepped up or down.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::TAU;
use std::time::Duration;

#[path = "support/cpu_time.rs"]
mod cpu_time;
use cpu_time::cpu_time;

use ogeom::algo::distance::{DistanceOptions, to_mesh};
use ogeom::core::FastMap;
use ogeom::math::{Point, Vector};
use ogeom::mesh::{Bvh, Editable, OFFSET_FARTHEST, OffsetReport};
use ogeom::topo::Triangulation;

/// The closed surface of a union of cubes `h` wide, the cube at `(i, j, k)`
/// spanning `[i h, (i + 1) h]` along `x` and so on: one square of two
/// triangles for each side a cube shares with no other, wound outward,
/// corners shared by lattice point.
fn voxels(cells: &[[i32; 3]], h: f64) -> Triangulation {
    let mut filled: FastMap<[i32; 3], ()> = FastMap::default();
    for &c in cells {
        filled.insert(c, ());
    }
    let mut at: FastMap<[i32; 3], u32> = FastMap::default();
    let mut mesh = Triangulation::new();
    let mut vertex = |l: [i32; 3], mesh: &mut Triangulation| {
        *at.entry(l).or_insert_with(|| {
            let [x, y, z] = l.map(|c| f64::from(c) * h);
            mesh.positions.push(Point::new(x, y, z));
            u32::try_from(mesh.positions.len() - 1).unwrap()
        })
    };
    for &c in cells {
        for axis in 0..3 {
            for side in [0, 1] {
                let mut beyond = c;
                beyond[axis] += if side == 1 { 1 } else { -1 };
                if filled.contains_key(&beyond) {
                    continue;
                }
                // The other two axes in cyclic order, `u × v` along `axis`.
                let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
                let mut corner = |du: i32, dv: i32| {
                    let mut l = c;
                    l[axis] += side;
                    l[u] += du;
                    l[v] += dv;
                    vertex(l, &mut mesh)
                };
                let mut quad = [corner(0, 0), corner(1, 0), corner(1, 1), corner(0, 1)];
                if side == 0 {
                    quad.reverse();
                }
                mesh.triangles.push([quad[0], quad[1], quad[2]]);
                mesh.triangles.push([quad[0], quad[2], quad[3]]);
            }
        }
    }
    mesh
}

/// The cells of an `a` by `b` by `c` box from the origin, less those in
/// `cut` where given, as half-open ranges per axis.
fn block(a: i32, b: i32, c: i32, cut: Option<[(i32, i32); 3]>) -> Vec<[i32; 3]> {
    let mut cells = Vec::new();
    for i in 0..a {
        for j in 0..b {
            for k in 0..c {
                let inside = |r: [(i32, i32); 3]| {
                    [i, j, k]
                        .iter()
                        .zip(r)
                        .all(|(&x, (lo, hi))| lo <= x && x < hi)
                };
                if !cut.is_some_and(inside) {
                    cells.push([i, j, k]);
                }
            }
        }
    }
    cells
}

/// The volume a closed mesh encloses, by the divergence theorem.
fn volume(m: &Triangulation) -> f64 {
    m.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| m.positions[v as usize] - Point::new(0.0, 0.0, 0.0));
            a.dot(b.cross(c)) / 6.0
        })
        .sum()
}

/// The mesh is closed, a manifold, and no face is folded.
fn sound(e: &Editable, report: &OffsetReport) {
    assert!(e.is_closed(), "the offset is open");
    assert!(e.non_manifold().is_empty(), "the offset is not a manifold");
    assert_eq!(report.folded, 0, "faces left folded");
}

/// Each vertex's unsigned distance to `original`, measured through its
/// hierarchy.
fn distances(original: &Triangulation, e: &Editable) -> Vec<f64> {
    let bvh = Bvh::build(original).unwrap();
    let points = e.to_triangulation().positions;
    let options = DistanceOptions {
        signed: false,
        ..DistanceOptions::default()
    };
    to_mesh(&bvh, &points, None, &options)
        .unwrap()
        .into_iter()
        .map(|c| c.unwrap().distance)
        .collect()
}

/// Rounding of a sum of a few lengths up to `scale`.
fn rounding(scale: f64) -> f64 {
    64.0 * f64::EPSILON * scale
}

/// How many coordinates of `p` stand at `lo` or `hi` on their axis, each
/// within rounding; `None` if one stands outside them.
fn on_box(p: Point, lo: [f64; 3], hi: [f64; 3]) -> Option<usize> {
    let mut at = 0;
    for (x, (l, h)) in [p.x, p.y, p.z].into_iter().zip(lo.into_iter().zip(hi)) {
        let tol = rounding(l.abs().max(h.abs()));
        if x < l - tol || x > h + tol {
            return None;
        }
        if (x - l).abs() <= tol || (x - h).abs() <= tol {
            at += 1;
        }
    }
    Some(at)
}

/// A cube of `n` cells `h` wide offset by `d` is the cube `d` larger on
/// every side, its edges and corners sharp: every vertex on its surface,
/// none outside it, at `d`, `√2 d` or `√3 d` from the cube as it stands on
/// one, two or three of its sides.
fn cube_offset_is_the_cube_with_sharp_edges(n: i32, h: f64, d: f64) -> (Editable, OffsetReport) {
    let original = voxels(&block(n, n, n, None), h);
    let mut e = Editable::from(&original, 0.0).unwrap();
    let report = e.offset(d, None).unwrap();
    eprintln!("cube of {n} cells {h} wide offset by {d}: {report:?}");
    sound(&e, &report);
    let side = f64::from(n) * h;
    let (lo, hi) = ([-d; 3], [side + d; 3]);
    let out = e.to_triangulation();
    let near = distances(&original, &e);
    let mut corners = 0;
    for (p, r) in out.positions.iter().zip(near) {
        let at = on_box(*p, lo, hi).unwrap_or_else(|| panic!("{p:?} is off the box"));
        assert!(at >= 1, "{p:?} is inside the box");
        corners += usize::from(at == 3);
        // Inward, every vertex is |d| from the faces it stands over.
        let want = if d > 0.0 {
            d * f64::from(u32::try_from(at).unwrap()).sqrt()
        } else {
            -d
        };
        assert!(
            (r - want).abs() <= rounding(side),
            "{p:?} is {r} from the cube, not {want}"
        );
    }
    assert_eq!(corners, 8);
    let s = side + 2.0 * d;
    assert!((volume(&out) - s * s * s).abs() <= rounding(s * s * s));
    (e, report)
}

#[test]
fn a_cube_offset_outward_is_the_larger_cube() {
    let (_, report) = cube_offset_is_the_cube_with_sharp_edges(6, 0.5, 0.2);
    assert_eq!(
        report,
        OffsetReport {
            merged: 0,
            eased: 0,
            capped: 0,
            folded: 0
        }
    );
}

#[test]
fn a_cube_offset_inward_past_its_cells_folds_and_is_repaired_to_the_smaller_cube() {
    // A step of 0.25 overshoots the edges by two rows of 0.1 cells.
    let (e, report) = cube_offset_is_the_cube_with_sharp_edges(10, 0.1, -0.25);
    assert!(report.merged > 0);
    assert!(e.vertex_count() > e.to_triangulation().positions.len());
}

#[test]
fn a_cube_offset_inward_within_its_cells_needs_no_repair() {
    let (_, report) = cube_offset_is_the_cube_with_sharp_edges(4, 0.5, -0.2);
    assert_eq!(report.merged, 0);
}

/// A block 2 by 1 by 1 of cells 0.1 wide with a slot 0.4 wide and 0.3 deep
/// across its top. Offset by `d`, with sharp edges, it is the block `d`
/// larger and the slot `d` narrower on each side, its floor `d` higher.
fn notched(d: f64) -> (Triangulation, Editable, OffsetReport) {
    let h = 0.1;
    let original = voxels(&block(20, 10, 10, Some([(8, 12), (-1, 11), (7, 11)])), h);
    let mut e = Editable::from(&original, 0.0).unwrap();
    let report = e.offset(d, None).unwrap();
    sound(&e, &report);
    let out = e.to_triangulation();
    let (l, w, t) = (2.0 + 2.0 * d, 1.0 + 2.0 * d, 1.0 + 2.0 * d);
    let (slot, depth) = (0.4 - 2.0 * d, 0.3);
    let want = l * w * t - slot * w * depth;
    assert!(
        (volume(&out) - want).abs() <= rounding(10.0),
        "the volume is {}, not {want}",
        volume(&out)
    );
    // Every vertex stands in the offset block, out of the offset slot.
    for p in &out.positions {
        let tol = rounding(2.0);
        let inside = |x: f64, lo: f64, hi: f64| x >= lo - tol && x <= hi + tol;
        assert!(inside(p.x, -d, 2.0 + d) && inside(p.y, -d, 1.0 + d) && inside(p.z, -d, 1.0 + d));
        let in_slot = p.x > 0.8 + d + tol && p.x < 1.2 - d - tol && p.z > 0.7 + d + tol;
        assert!(!in_slot, "{p:?} is in the slot");
    }
    (original, e, report)
}

#[test]
fn a_notched_block_offset_inward_is_repaired_to_the_offset_block() {
    let d = -0.25;
    let (original, e, report) = notched(d);
    assert!(report.merged > 0);
    // No vertex comes nearer the block than |d|, and none stands farther
    // than a corner's √3 |d|: those beyond |d| lie on the planes extended
    // past the slot's concave edges, which the offset keeps sharp.
    let reach = d.abs();
    let (mut nearest, mut farthest, mut at_reach) = (f64::INFINITY, 0.0_f64, 0);
    for r in distances(&original, &e) {
        nearest = nearest.min(r);
        farthest = farthest.max(r);
        at_reach += usize::from((r - reach).abs() <= rounding(2.0));
    }
    eprintln!(
        "notched block offset by {d}: {at_reach} vertices at |d|, distances from {nearest} to {farthest}, {report:?}"
    );
    assert!(nearest >= reach - rounding(2.0));
    assert!(farthest <= 3.0_f64.sqrt() * reach + rounding(2.0));
}

#[test]
fn a_notched_block_offset_outward_closes_its_slot_by_twice_the_step() {
    let (_, _, report) = notched(0.15);
    assert!(report.merged > 0);
}

/// The largest angle between a vertex's direction from the origin and the
/// normal of a face at it.
fn largest_tilt(m: &Triangulation) -> f64 {
    let mut most = 0.0_f64;
    for t in &m.triangles {
        let [a, b, c] = t.map(|v| m.positions[v as usize]);
        let n = (b - a).cross(c - a);
        for p in [a, b, c] {
            let r = p - Point::new(0.0, 0.0, 0.0);
            most = most.max((n.dot(r) / (n.magnitude() * r.magnitude())).acos());
        }
    }
    most
}

#[test]
fn a_sphere_offset_either_way_is_the_sphere_of_radius_r_plus_d() {
    let radius = 10.0;
    let original = sphere(radius, 4);
    // Each face tilts from its corners' radii by at most `φ`, so least
    // squares over the faces puts a vertex between `|d|` and `|d| / cos φ`
    // from where it stood.
    let phi = largest_tilt(&original);
    for d in [0.5, -0.5] {
        let mut e = Editable::from(&original, 0.0).unwrap();
        let report = e.offset(d, None).unwrap();
        sound(&e, &report);
        assert_eq!(report.merged + report.eased + report.capped, 0);
        let bound = d.abs() * (1.0 / phi.cos() - 1.0) + rounding(radius);
        let out = e.to_triangulation();
        let mut worst = 0.0_f64;
        for p in &out.positions {
            let r = (*p - Point::new(0.0, 0.0, 0.0)).magnitude();
            worst = worst.max((r - (radius + d)).abs());
        }
        let deviation = distances(&original, &e)
            .into_iter()
            .map(|r| (r - d.abs()).abs())
            .fold(0.0, f64::max);
        eprintln!(
            "sphere offset by {d}: radius off by {worst:e}, distance off by {deviation:e}, bound {bound:e}"
        );
        assert!(worst <= bound);
        assert!(deviation <= bound);
    }
}

/// A plate 1 by 1 by 0.1 of cells 0.1 wide, its top's vertices within the
/// square `[0.3, 0.7]²` the region.
fn plate_with_region() -> (Triangulation, Vec<bool>) {
    let plate = voxels(&block(10, 10, 1, None), 0.1);
    let tol = rounding(1.0);
    let region = plate
        .positions
        .iter()
        .map(|p| {
            (p.z - 0.1).abs() <= tol
                && (0.3 - tol..=0.7 + tol).contains(&p.x)
                && (0.3 - tol..=0.7 + tol).contains(&p.y)
        })
        .collect();
    (plate, region)
}

#[test]
fn a_region_of_a_plate_steps_up_or_down_by_d_behind_a_closed_wall() {
    for d in [0.05, -0.05] {
        let (plate, region) = plate_with_region();
        let mut e = Editable::from(&plate, 0.0).unwrap();
        let before = (e.vertex_count(), e.face_count());
        let report = e.offset(d, Some(&region)).unwrap();
        sound(&e, &report);
        // The square's rim has 16 vertices and 16 edges: each vertex gains
        // a foot and each edge a wall of two faces.
        assert_eq!(e.vertex_count(), before.0 + 16);
        assert_eq!(e.face_count(), before.1 + 32);
        for (v, &inside) in region.iter().enumerate() {
            let v = u32::try_from(v).unwrap();
            let (p, q) = (plate.positions[v as usize], e.position(v).unwrap());
            let lift = if inside { d } else { 0.0 };
            assert_eq!(q.x, p.x);
            assert_eq!(q.y, p.y);
            assert!((q.z - (p.z + lift)).abs() <= rounding(1.0));
        }
        let want = 0.1 + 0.4 * 0.4 * d;
        let got = volume(&e.to_triangulation());
        assert!(
            (got - want).abs() <= rounding(1.0),
            "volume {got}, not {want}"
        );
    }
}

#[test]
fn a_region_left_empty_moves_nothing() {
    let (plate, region) = plate_with_region();
    let mut e = Editable::from(&plate, 0.0).unwrap();
    let report = e.offset(0.05, Some(&vec![false; region.len()])).unwrap();
    assert_eq!(report.merged, 0);
    assert_eq!(e.to_triangulation().positions, plate.positions);
}

#[test]
fn a_knife_edge_moves_no_farther_than_the_cap() {
    // A flat tetrahedron: along its base, faces whose normals turn by
    // about 170° meet, and their planes offset meet far beyond the cap.
    let positions = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
        Point::new(0.3, 0.3, 0.05),
    ];
    let triangles = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]];
    let mesh = Triangulation {
        positions,
        triangles,
        ..Triangulation::new()
    };
    let mut e = Editable::from(&mesh, 0.0).unwrap();
    let d = 0.01;
    let report = e.offset(d, None).unwrap();
    assert_eq!(report.capped, 3);
    for v in 0..4 {
        let moved = (e.position(v).unwrap() - mesh.positions[v as usize]).magnitude();
        assert!(moved <= OFFSET_FARTHEST * d);
    }
}

#[test]
fn the_offset_is_deterministic() {
    let run = || {
        let original = voxels(&block(20, 10, 10, Some([(8, 12), (-1, 11), (7, 11)])), 0.1);
        let mut e = Editable::from(&original, 0.0).unwrap();
        let report = e.offset(-0.25, None).unwrap();
        (e.to_triangulation(), report)
    };
    let (a, ra) = run();
    let (b, rb) = run();
    assert_eq!(ra, rb);
    assert_eq!(a.triangles, b.triangles);
    assert_eq!(a.positions, b.positions);
}

#[test]
fn a_bad_distance_or_region_is_refused() {
    let mesh = voxels(&block(2, 2, 2, None), 1.0);
    let mut e = Editable::from(&mesh, 0.0).unwrap();
    assert!(e.offset(f64::NAN, None).is_err());
    assert!(e.offset(f64::INFINITY, None).is_err());
    assert!(e.offset(1.0, Some(&[true])).is_err());
    let report = e.offset(0.0, None).unwrap();
    assert_eq!(report.merged, 0);
    assert_eq!(e.to_triangulation().positions, mesh.positions);
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
}

/// A scan-sized ball, five million triangles about 0.1 apart under noise of
/// σ = 0.01 along its radius, offset inward by 0.2: closed, a manifold, and
/// its mean radius 0.2 smaller, to within what least squares over faces
/// tilted by up to `φ` adds, `|d| (1 / cos φ − 1)`.
#[test]
#[ignore = "heavy"]
fn five_million_triangles_offset() {
    let radius = 50.0;
    let mut mesh = sphere(radius, 9);
    let mut noise = Noise(11);
    for p in &mut mesh.positions {
        let r = *p - Point::new(0.0, 0.0, 0.0);
        *p += r * (noise.gaussian(0.01) / r.magnitude());
    }
    let mean = |m: &Triangulation| {
        let sum: f64 = m
            .positions
            .iter()
            .map(|p| (*p - Point::new(0.0, 0.0, 0.0)).magnitude())
            .sum();
        sum / f64::from(u32::try_from(m.positions.len()).unwrap())
    };
    let d = -0.2;
    let mut e = Editable::from(&mesh, 0.0).unwrap();
    let (report, took) = cpu_time(|| e.offset(d, None).unwrap());
    let out = e.to_triangulation();
    let (was, now) = (mean(&mesh), mean(&out));
    eprintln!("5M triangles: offset in {took:?}, mean radius {was:.5} to {now:.5}, {report:?}");
    assert!(e.is_closed());
    assert!(e.non_manifold().is_empty());
    let phi = largest_tilt(&mesh);
    assert!((now - was - d).abs() <= d.abs() * (1.0 / phi.cos() - 1.0));
    assert!(took < Duration::from_secs(120), "{took:?} to offset");
}
