#![allow(clippy::unwrap_used, clippy::cast_precision_loss)]

use super::*;

/// An axis-aligned box from `lo` to `hi`, wound outward, with its square
/// faces split along the diagonal `split` picks.
fn cuboid(lo: P3, hi: P3) -> Triangulation {
    let mut mesh = Triangulation::new();
    for k in 0..8_u32 {
        let pick = |b: u32, axis: usize| {
            if (k >> b) & 1 == 1 {
                hi[axis]
            } else {
                lo[axis]
            }
        };
        mesh.positions
            .push(Point::new(pick(0, 0), pick(1, 1), pick(2, 2)));
    }
    let quads = [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ];
    for [a, b, c, d] in quads {
        mesh.triangles.push([a, b, c]);
        mesh.triangles.push([a, c, d]);
    }
    mesh
}

/// A sphere of latitude rings and longitude slices, wound outward.
fn sphere(centre: P3, radius: f64, rings: u32, slices: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    let at = |theta: f64, phi: f64| {
        Point::new(
            radius.mul_add(theta.sin() * phi.cos(), centre[0]),
            radius.mul_add(theta.sin() * phi.sin(), centre[1]),
            radius.mul_add(theta.cos(), centre[2]),
        )
    };
    mesh.positions.push(at(0.0, 0.0));
    for r in 1..rings {
        let theta = std::f64::consts::PI * f64::from(r) / f64::from(rings);
        for s in 0..slices {
            let phi = std::f64::consts::TAU * f64::from(s) / f64::from(slices);
            mesh.positions.push(at(theta, phi));
        }
    }
    mesh.positions.push(at(std::f64::consts::PI, 0.0));
    let south = rings - 1;
    let ring = |r: u32, s: u32| 1 + (r - 1) * slices + s % slices;
    let last = id(mesh.positions.len() - 1);
    for s in 0..slices {
        mesh.triangles.push([0, ring(1, s), ring(1, s + 1)]);
        for r in 1..south {
            let (a, b) = (ring(r, s), ring(r, s + 1));
            let (c, d) = (ring(r + 1, s), ring(r + 1, s + 1));
            mesh.triangles.push([a, c, d]);
            mesh.triangles.push([a, d, b]);
        }
        mesh.triangles
            .push([ring(south, s), last, ring(south, s + 1)]);
    }
    mesh
}

/// A prism about the z axis through `(cx, cy)`: a regular `n`-gon of
/// circumradius `radius`, from `z0` to `z1`.
fn prism(cx: f64, cy: f64, radius: f64, n: u32, z0: f64, z1: f64) -> Triangulation {
    let mut mesh = Triangulation::new();
    for z in [z0, z1] {
        mesh.positions.push(Point::new(cx, cy, z));
        for s in 0..n {
            let phi = std::f64::consts::TAU * f64::from(s) / f64::from(n);
            mesh.positions.push(Point::new(
                radius.mul_add(phi.cos(), cx),
                radius.mul_add(phi.sin(), cy),
                z,
            ));
        }
    }
    let (bottom, top) = (0, n + 1);
    let rim = |base: u32, s: u32| base + 1 + s % n;
    for s in 0..n {
        mesh.triangles
            .push([bottom, rim(bottom, s + 1), rim(bottom, s)]);
        mesh.triangles.push([top, rim(top, s), rim(top, s + 1)]);
        let (a, b) = (rim(bottom, s), rim(bottom, s + 1));
        let (c, d) = (rim(top, s), rim(top, s + 1));
        mesh.triangles.push([a, b, d]);
        mesh.triangles.push([a, d, c]);
    }
    mesh
}

fn moved(mesh: &Triangulation, f: impl Fn(Point) -> Point) -> Triangulation {
    let mut out = mesh.clone();
    for p in &mut out.positions {
        *p = f(*p);
    }
    out
}

fn closed(mesh: &Triangulation) {
    check_result(mesh, &vec![true; mesh.positions.len()]).unwrap();
}

fn run(a: &Triangulation, b: &Triangulation, op: BoolKind) -> MeshBoolean {
    let out = boolean(a, b, op).unwrap();
    closed(&out.mesh);
    assert_eq!(out.cut.len(), out.mesh.triangles.len());
    out
}

fn message(e: &ogeom_core::OgeomError) -> String {
    e.to_string()
}

#[test]
fn overlapping_cubes_give_the_closed_form_volumes() {
    let a = cuboid([0.0; 3], [2.0; 3]);
    let b = cuboid([1.0; 3], [3.0; 3]);
    for (op, want) in [
        (BoolKind::Union, 15.0),
        (BoolKind::Intersection, 1.0),
        (BoolKind::Difference, 7.0),
    ] {
        let out = run(&a, &b, op);
        assert!((out.mesh.volume() - want).abs() < 1e-12, "{op:?}");
        assert_eq!(out.rounding, 0.0, "every crossing here is a dyadic point");
    }
}

#[test]
fn a_cube_less_a_sphere_inside_it_is_the_cube_with_the_sphere_as_a_cavity() {
    let a = cuboid([0.0; 3], [4.0; 3]);
    let b = sphere([2.0, 2.1, 1.9], 1.5, 24, 32);
    let vb = b.volume();
    let diff = run(&a, &b, BoolKind::Difference);
    // Nothing is cut, so every vertex is an input vertex.
    assert_eq!(diff.first_new_vertex, diff.mesh.positions.len());
    assert!(diff.cut.iter().all(|c| !c));
    assert!((diff.mesh.volume() - (64.0 - vb)).abs() < 1e-12);
    let common = run(&a, &b, BoolKind::Intersection);
    assert!((common.mesh.volume() - vb).abs() < 1e-12);
    let union = run(&a, &b, BoolKind::Union);
    assert!((union.mesh.volume() - 64.0).abs() < 1e-12);
}

/// The volumes of union, intersection and difference agree with the
/// operands' own to rounding: `|A ∪ B| + |A ∩ B| = |A| + |B|` and
/// `|A - B| + |A ∩ B| = |A|`.
fn inclusion_exclusion(a: &Triangulation, b: &Triangulation) -> [MeshBoolean; 3] {
    let union = run(a, b, BoolKind::Union);
    let common = run(a, b, BoolKind::Intersection);
    let diff = run(a, b, BoolKind::Difference);
    let (va, vb) = (a.volume(), b.volume());
    let scale = va.abs() + vb.abs();
    let (u, i, d) = (
        union.mesh.volume(),
        common.mesh.volume(),
        diff.mesh.volume(),
    );
    assert!(
        (u + i - va - vb).abs() <= 1e-9 * scale,
        "{u} + {i} vs {va} + {vb}"
    );
    assert!((d + i - va).abs() <= 1e-9 * scale, "{d} + {i} vs {va}");
    assert!(i >= -1e-12 && i <= va.min(vb) + 1e-9 * scale);
    for out in [&union, &common, &diff] {
        let n = a.positions.len();
        // The input vertices kept are where they were.
        for p in &out.mesh.positions[..out.first_new_vertex] {
            assert!(a.positions.contains(p) || b.positions.contains(p));
        }
        let _ = n;
    }
    [union, common, diff]
}

#[test]
fn a_cube_less_a_sphere_through_its_face_balances_its_volumes() {
    let a = cuboid([0.0; 3], [4.0; 3]);
    let b = sphere([2.0, 2.1, 4.0], 1.5, 24, 32);
    let [_, common, diff] = inclusion_exclusion(&a, &b);
    // The sphere's cap below its equator, and a cut face flagged.
    assert!(common.mesh.volume() > 0.0 && common.mesh.volume() < b.volume());
    assert!(diff.cut.iter().any(|&c| c));
    assert!(diff.first_new_vertex < diff.mesh.positions.len());
    assert!(diff.rounding > 0.0 && diff.rounding < 1e-14);
}

#[test]
fn a_box_less_a_prism_through_it_leaves_the_prism_s_hole() {
    // The issue's box and cylinder mesh: the hole is a regular 64-gon
    // prism, whose volume is a closed form.
    let a = cuboid([0.0, 0.0, 0.0], [4.0, 4.0, 2.0]);
    let n = 64;
    let b = prism(2.0, 2.0, 1.0, n, -1.0, 3.0);
    let out = run(&a, &b, BoolKind::Difference);
    let polygon = f64::from(n) / 2.0 * (std::f64::consts::TAU / f64::from(n)).sin();
    let hole = 2.0 * polygon;
    assert!((out.mesh.volume() - (32.0 - hole)).abs() < 1e-12);
    // A torus: Euler characteristic 0.
    let edges = out.mesh.triangles.len() * 3 / 2;
    assert_eq!(out.mesh.positions.len() + out.mesh.triangles.len(), edges);
}

#[test]
fn stacked_cubes_sharing_a_face() {
    let a = cuboid([0.0; 3], [1.0; 3]);
    let b = cuboid([0.0, 0.0, 1.0], [1.0, 1.0, 2.0]);
    let union = run(&a, &b, BoolKind::Union);
    assert!((union.mesh.volume() - 2.0).abs() < 1e-12);
    // The shared face is gone: a box of eight corners.
    assert_eq!(union.mesh.positions.len(), 12);
    let common = run(&a, &b, BoolKind::Intersection);
    assert!(common.mesh.is_empty());
    let diff = run(&a, &b, BoolKind::Difference);
    assert!((diff.mesh.volume() - 1.0).abs() < 1e-12);
    // Shifted across the shared plane: the faces overlap in part.
    let c = cuboid([0.5, 0.25, 1.0], [1.5, 1.25, 2.0]);
    let [union, common, diff] = inclusion_exclusion(&a, &c);
    assert!((union.mesh.volume() - 2.0).abs() < 1e-12);
    assert!(common.mesh.is_empty());
    assert!((diff.mesh.volume() - 1.0).abs() < 1e-12);
}

#[test]
fn a_cube_flush_cut_by_a_box_sharing_four_of_its_planes() {
    let a = cuboid([0.0; 3], [2.0; 3]);
    let b = cuboid([1.0, 0.0, 0.0], [3.0, 2.0, 2.0]);
    let [union, common, diff] = inclusion_exclusion(&a, &b);
    assert!((union.mesh.volume() - 12.0).abs() < 1e-12);
    assert!((common.mesh.volume() - 4.0).abs() < 1e-12);
    assert!((diff.mesh.volume() - 4.0).abs() < 1e-12);
    // Flush on one face only, the other operand standing proud of the rest.
    let c = cuboid([1.0, -1.0, -1.0], [3.0, 2.0, 3.0]);
    let [union, common, diff] = inclusion_exclusion(&a, &c);
    assert!((diff.mesh.volume() - 4.0).abs() < 1e-12);
    assert!((common.mesh.volume() - 4.0).abs() < 1e-12);
    assert!((union.mesh.volume() - 28.0).abs() < 1e-12);
}

#[test]
fn refusals_name_their_cause() {
    let a = cuboid([0.0; 3], [2.0; 3]);
    let b = cuboid([1.0; 3], [3.0; 3]);
    let mut open = b.clone();
    open.triangles.pop();
    let e = boolean(&a, &open, BoolKind::Union).unwrap_err();
    assert!(message(&e).contains("second mesh is open"), "{e}");
    let far = cuboid([5.0; 3], [6.0; 3]);
    let e = boolean(&a, &far, BoolKind::Union).unwrap_err();
    assert!(message(&e).contains("do not meet"), "{e}");
    let mut inside_out = b.clone();
    for t in &mut inside_out.triangles {
        t.swap(1, 2);
    }
    let e = boolean(&a, &inside_out, BoolKind::Union).unwrap_err();
    assert!(message(&e).contains("not wound outward"), "{e}");
    // Two cubes in one mesh, overlapping.
    let mut crossing = a.clone();
    let shifted = cuboid([1.0; 3], [3.0; 3]);
    let base = id(crossing.positions.len());
    crossing.positions.extend(&shifted.positions);
    crossing
        .triangles
        .extend(shifted.triangles.iter().map(|t| t.map(|v| v + base)));
    let e = boolean(&crossing, &b, BoolKind::Union).unwrap_err();
    assert!(message(&e).contains("first mesh crosses itself"), "{e}");
    let mut flat = a.clone();
    flat.positions[7] = flat.positions[6];
    assert!(boolean(&flat, &b, BoolKind::Union).is_err());
    // Touching along an edge only: the union would not be 2-manifold.
    let edge = cuboid([2.0, 2.0, 0.0], [3.0, 3.0, 2.0]);
    let e = boolean(&a, &edge, BoolKind::Union).unwrap_err();
    assert!(message(&e).contains("not be 2-manifold"), "{e}");
}

#[test]
fn a_held_operand_meets_without_touching() {
    let a = cuboid([0.0; 3], [4.0; 3]);
    let b = cuboid([1.0; 3], [2.0; 3]);
    let diff = run(&a, &b, BoolKind::Difference);
    assert!((diff.mesh.volume() - 63.0).abs() < 1e-12);
    let union = run(&b, &a, BoolKind::Union);
    assert!((union.mesh.volume() - 64.0).abs() < 1e-12);
    let none = run(&b, &a, BoolKind::Difference);
    assert!(none.mesh.is_empty());
}

#[test]
fn the_answer_is_the_same_every_time() {
    let a = sphere([0.0; 3], 1.0, 16, 24);
    let b = moved(&sphere([0.0; 3], 0.8, 12, 20), |p| {
        Point::new(p.x + 0.7, p.y + 0.31, p.z - 0.2)
    });
    let first = boolean(&a, &b, BoolKind::Difference).unwrap();
    for _ in 0..3 {
        assert_eq!(boolean(&a, &b, BoolKind::Difference).unwrap(), first);
    }
}

/// A small deterministic generator for the fuzz tests.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }

    /// A rotation, as a unit quaternion's matrix, and a translation.
    fn motion(&mut self, reach: f64) -> impl Fn(Point) -> Point + use<> {
        let mut q = [0.0_f64; 4];
        for v in &mut q {
            *v = self.next() - 0.5;
        }
        let n = q.iter().map(|v| v * v).sum::<f64>().sqrt();
        let [w, x, y, z] = q.map(|v| v / n);
        let m = [
            [
                1.0 - 2.0 * (y * y + z * z),
                2.0 * (x * y - w * z),
                2.0 * (x * z + w * y),
            ],
            [
                2.0 * (x * y + w * z),
                1.0 - 2.0 * (x * x + z * z),
                2.0 * (y * z - w * x),
            ],
            [
                2.0 * (x * z - w * y),
                2.0 * (y * z + w * x),
                1.0 - 2.0 * (x * x + y * y),
            ],
        ];
        let t = [0, 1, 2].map(|_| (self.next() - 0.5) * reach);
        move |p: Point| {
            let r = |k: usize| m[k][0] * p.x + m[k][1] * p.y + m[k][2] * p.z + t[k];
            Point::new(r(0), r(1), r(2))
        }
    }
}

#[test]
fn random_rigid_motions_balance_their_volumes() {
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    let shapes = [
        cuboid([-0.5; 3], [0.5; 3]),
        sphere([0.0; 3], 0.6, 10, 14),
        prism(0.0, 0.0, 0.5, 9, -0.7, 0.7),
    ];
    for case in 0..60 {
        let a = moved(&shapes[case % 3], rng.motion(0.6));
        let b = moved(&shapes[(case / 3) % 3], rng.motion(0.6));
        inclusion_exclusion(&a, &b);
    }
}

#[test]
fn grid_aligned_cubes_meet_in_every_degenerate_way() {
    // Cubes on a half-unit grid share faces, edges and corners in every
    // combination. Each answer balances its volumes; the only refusal
    // allowed is the one for a result that touches itself.
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let unit = cuboid([0.0; 3], [1.0; 3]);
    let mut answered = 0;
    for _ in 0..200 {
        let step = |rng: &mut Rng| (rng.next() * 5.0).floor() / 2.0 - 1.0;
        let (dx, dy, dz) = (step(&mut rng), step(&mut rng), step(&mut rng));
        let b = moved(&unit, |p| Point::new(p.x + dx, p.y + dy, p.z + dz));
        let ops = [
            BoolKind::Union,
            BoolKind::Intersection,
            BoolKind::Difference,
        ];
        let outs: Vec<_> = ops.iter().map(|&op| boolean(&unit, &b, op)).collect();
        if outs.iter().any(Result::is_err) {
            for e in outs.iter().filter_map(|r| r.as_ref().err()) {
                let m = message(e);
                assert!(
                    m.contains("not be 2-manifold") || m.contains("do not meet"),
                    "{m}"
                );
            }
            continue;
        }
        let v: Vec<f64> = outs
            .iter()
            .map(|r| {
                let out = r.as_ref().unwrap();
                closed(&out.mesh);
                out.mesh.volume()
            })
            .collect();
        assert!((v[0] + v[1] - 2.0).abs() < 1e-12, "{dx} {dy} {dz}: {v:?}");
        assert!((v[2] + v[1] - 1.0).abs() < 1e-12, "{dx} {dy} {dz}: {v:?}");
        let exact =
            (1.0 - dx.abs()).max(0.0) * (1.0 - dy.abs()).max(0.0) * (1.0 - dz.abs()).max(0.0);
        assert!((v[1] - exact).abs() < 1e-12, "{dx} {dy} {dz}: {v:?}");
        answered += 1;
    }
    assert!(answered > 50);
}

#[test]
fn an_operand_against_itself_is_all_shared_faces() {
    for a in [
        cuboid([0.0; 3], [1.0; 3]),
        sphere([0.1, 0.2, 0.3], 0.7, 9, 13),
    ] {
        let va = a.volume();
        let union = run(&a, &a, BoolKind::Union);
        assert!((union.mesh.volume() - va).abs() < 1e-12);
        assert_eq!(union.mesh.triangles.len(), a.triangles.len());
        let common = run(&a, &a, BoolKind::Intersection);
        assert!((common.mesh.volume() - va).abs() < 1e-12);
        let diff = run(&a, &a, BoolKind::Difference);
        assert!(diff.mesh.is_empty());
    }
}

#[test]
fn a_prism_against_its_moved_copy_shares_faces_or_nearly() {
    // Raised by half its height, the copy shares the upper half of every
    // side face exactly.
    let a = prism(0.0, 0.0, 1.0, 12, 0.0, 1.0);
    let up = moved(&a, |p| Point::new(p.x, p.y, p.z + 0.5));
    let [union, common, diff] = inclusion_exclusion(&a, &up);
    let va = a.volume();
    assert!((union.mesh.volume() - 1.5 * va).abs() < 1e-12);
    assert!((common.mesh.volume() - 0.5 * va).abs() < 1e-12);
    assert!((diff.mesh.volume() - 0.5 * va).abs() < 1e-12);
    // Turned a half turn, each vertex lands within rounding of another, so
    // the side faces cross at angles of the order of 1e-16.
    let turned = moved(&a, |p| Point::new(-p.x, -p.y, p.z + 0.25));
    let [union, _, _] = inclusion_exclusion(&a, &turned);
    assert!((union.mesh.volume() - 1.25 * va).abs() < 1e-9);
    assert!(union.rounding < 1e-14);
}

#[test]
fn nearly_coincident_copies_balance_their_volumes() {
    // A mesh against itself moved by a rotation and a translation of the
    // order of 1e-15 to 1e-9: faces cross at vanishing angles everywhere.
    let mut rng = Rng(0xd1b5_4a32_d192_ed03);
    let shapes = [
        sphere([0.0; 3], 0.6, 8, 11),
        prism(0.0, 0.0, 0.5, 7, -0.4, 0.4),
    ];
    for case in 0..24 {
        let a = &shapes[case % 2];
        let size = 10_f64.powi(-15 + i32::try_from(case % 7).unwrap());
        let (ax, ay) = ((rng.next() - 0.5) * size, (rng.next() - 0.5) * size);
        let t = [0, 1, 2].map(|_| (rng.next() - 0.5) * size);
        let b = moved(a, |p| {
            // A small rotation about x then y, to first order.
            let (y, z) = (p.y - ax * p.z, p.z + ax * p.y);
            let (x, z) = (p.x + ay * z, z - ay * p.x);
            Point::new(x + t[0], y + t[1], z + t[2])
        });
        inclusion_exclusion(a, &b);
    }
}

#[test]
#[ignore = "heavy"]
fn two_meshes_of_a_hundred_thousand_triangles() {
    let a = sphere([0.0; 3], 1.0, 224, 224);
    let b = moved(&sphere([0.0; 3], 0.9, 224, 224), |p| {
        Point::new(p.x + 0.55, p.y + 0.21, p.z - 0.13)
    });
    assert!(a.triangles.len() >= 99_000 && b.triangles.len() >= 99_000);
    let start = ogeom_core::clock::Instant::now();
    let diff = run(&a, &b, BoolKind::Difference);
    let took = start.elapsed();
    eprintln!(
        "difference of {} and {} triangles: {} triangles in {took:?}",
        a.triangles.len(),
        b.triangles.len(),
        diff.mesh.triangles.len()
    );
    let common = run(&a, &b, BoolKind::Intersection);
    assert!((diff.mesh.volume() + common.mesh.volume() - a.volume()).abs() < 1e-9);
    assert!(took < std::time::Duration::from_secs(3), "{took:?}");
}

#[test]
fn a_box_mesh_less_a_cylinder_solid_holes_it_within_the_chord() {
    use ogeom_algo::make_cylinder;
    use ogeom_math::{Direction, Frame};
    const T: Tolerances = Tolerances::millimetres();
    let mut model = Model::new();
    let frame = Frame::new(Point::new(0.0, 0.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let cylinder = make_cylinder(&mut model, frame, 1.0, 4.0, T).unwrap();
    let block = cuboid([-2.0, -2.0, 0.0], [2.0, 2.0, 2.0]);
    let chord = 0.01;
    let out = boolean_with_shape(
        &block,
        &model,
        &cylinder.shape,
        BoolKind::Difference,
        chord,
        T,
    )
    .unwrap();
    closed(&out.mesh);
    // The hole is a polygon inscribed in the circle, its vertices on it:
    // short of the disc by at most the chord across its perimeter.
    let disc = std::f64::consts::PI;
    let hole = 32.0 - out.mesh.volume();
    assert!(hole <= 2.0 * disc + 1e-9, "{hole}");
    assert!(
        hole >= 2.0 * (disc - std::f64::consts::TAU * chord),
        "{hole}"
    );
    assert!(out.cut.iter().any(|&c| c));
}
