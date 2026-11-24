#![allow(
    clippy::unwrap_used,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]

use super::super::*;
use super::*;

/// An axis-aligned box from `lo` to `hi`, wound outward.
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

/// A sphere of latitude rings and longitude slices, wound outward, turned
/// about two axes so that no face lines up with another sphere's.
fn sphere(centre: P3, radius: f64, rings: u32, slices: u32, turn: f64) -> Triangulation {
    let mut mesh = Triangulation::new();
    let (s1, c1) = turn.sin_cos();
    let (s2, c2) = (1.7 * turn).sin_cos();
    let at = |theta: f64, phi: f64| {
        let (x, y, z) = (
            theta.sin() * phi.cos(),
            theta.sin() * phi.sin(),
            theta.cos(),
        );
        let (y, z) = (c1 * y - s1 * z, s1 * y + c1 * z);
        let (x, y) = (c2 * x - s2 * y, s2 * x + c2 * y);
        Point::new(
            radius.mul_add(x, centre[0]),
            radius.mul_add(y, centre[1]),
            radius.mul_add(z, centre[2]),
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

/// The meshes one after another as one mesh, sharing no vertex.
fn joined(parts: &[&Triangulation]) -> Triangulation {
    let mut out = Triangulation::new();
    for part in parts {
        let base = id(out.positions.len());
        out.positions.extend(&part.positions);
        out.triangles
            .extend(part.triangles.iter().map(|t| t.map(|v| v + base)));
    }
    out
}

fn inside_out(mesh: &Triangulation) -> Triangulation {
    let mut out = mesh.clone();
    for t in &mut out.triangles {
        t.swap(1, 2);
    }
    out
}

/// Every pair that meets away from its shared corners, by testing all of
/// them with the boolean's own predicate.
fn brute_force(mesh: &Triangulation) -> Vec<[u32; 2]> {
    let tri = |t: usize, side: Side| {
        let v = mesh.triangles[t];
        Tri {
            side,
            id: id(t),
            v,
            p: v.map(|k| coords(mesh.positions[k as usize])),
        }
    };
    let mut out = Vec::new();
    for t in 0..mesh.triangles.len() {
        for u in t + 1..mesh.triangles.len() {
            if !meet_at_shared_only(&tri(t, 0), &tri(u, 1)) {
                out.push([id(t), id(u)]);
            }
        }
    }
    out
}

/// The pairs, checked against [`brute_force`].
fn pairs_of(mesh: &Triangulation) -> Vec<[u32; 2]> {
    let pairs = self_intersections(mesh).unwrap();
    assert_eq!(pairs, brute_force(mesh));
    pairs
}

/// The answer, checked to be closed, 2-manifold, outward and embedded,
/// with a flag for each triangle and flags exactly on pieces that are not
/// input triangles.
fn resolved(mesh: &Triangulation) -> Resolved {
    let out = resolve_self_intersections(mesh).unwrap();
    assert_eq!(out.cut.len(), out.mesh.triangles.len());
    assert_eq!(out.crossings, self_intersections(mesh).unwrap().len());
    if out.mesh.triangles.is_empty() {
        return out;
    }
    let solid = Operand::new(&out.mesh, 0, "resolved").unwrap();
    solid.check_wound("resolved", None, 0).unwrap();
    assert_eq!(
        self_intersections(&out.mesh).unwrap(),
        Vec::<[u32; 2]>::new()
    );
    let key = |m: &Triangulation, t: [u32; 3]| {
        let mut c = t.map(|v| coords(m.positions[v as usize]).map(f64::to_bits));
        c.sort_unstable();
        c
    };
    let input: FastSet<[[u64; 3]; 3]> = mesh.triangles.iter().map(|&t| key(mesh, t)).collect();
    for (t, &cut) in out.mesh.triangles.iter().zip(&out.cut) {
        assert_eq!(cut, !input.contains(&key(&out.mesh, *t)), "{t:?}");
    }
    out
}

fn message(e: &ogeom_core::OgeomError) -> String {
    e.to_string()
}

#[test]
fn two_overlapping_cubes_cross_one_from_each_and_resolve_to_their_union() {
    let a = cuboid([0.0; 3], [10.0; 3]);
    let b = cuboid([5.0, 0.0, 0.0], [15.0, 10.0, 10.0]);
    let two = joined(&[&a, &b]);
    let pairs = pairs_of(&two);
    assert!(!pairs.is_empty());
    assert!(pairs.iter().all(|&[t, u]| t < 12 && u >= 12), "{pairs:?}");
    let out = resolved(&two);
    assert!((out.mesh.volume() - 1500.0).abs() < 1e-9);
    assert!(out.cut.iter().any(|&c| c));
    assert_eq!(out.rounding, 0.0, "every crossing here is a dyadic point");
}

#[test]
fn cubes_overlapping_in_a_corner_resolve_to_the_closed_form_union() {
    let a = cuboid([0.0; 3], [2.0; 3]);
    let b = cuboid([1.0; 3], [3.0; 3]);
    let out = resolved(&joined(&[&a, &b]));
    assert!((out.mesh.volume() - 15.0).abs() < 1e-12);
}

#[test]
fn embedded_meshes_have_no_pairs_and_come_back_as_they_were() {
    let one = cuboid([0.0; 3], [10.0; 3]);
    let two = joined(&[&one, &cuboid([20.0, 0.0, 0.0], [30.0, 10.0, 10.0])]);
    let ball = sphere([0.0; 3], 3.0, 9, 14, 0.3);
    for mesh in [&one, &two, &ball] {
        assert!(pairs_of(mesh).is_empty());
        let out = resolved(mesh);
        assert_eq!(out.crossings, 0);
        assert_eq!(out.mesh.positions, mesh.positions);
        assert_eq!(out.mesh.triangles, mesh.triangles);
        assert!(out.cut.iter().all(|&c| !c));
        assert_eq!(out.rounding, 0.0);
    }
}

#[test]
fn a_vertex_pulled_through_the_opposite_face_folds_its_triangles_through_it() {
    let mut cube = cuboid([0.0; 3], [10.0; 3]);
    // Vertex 7 is the corner (10, 10, 10); the face z = 0 is triangles 0, 1.
    cube.positions[7] = Point::new(6.0, 7.0, -5.0);
    let pairs = pairs_of(&cube);
    assert!(!pairs.is_empty());
    let at_seven = |t: u32| cube.triangles[t as usize].contains(&7);
    for &[t, u] in &pairs {
        assert!(at_seven(t) != at_seven(u), "{t} {u}");
    }
    let passed: FastSet<u32> = pairs
        .iter()
        .map(|&[t, u]| if at_seven(t) { u } else { t })
        .collect();
    assert!(passed.contains(&0) || passed.contains(&1), "{pairs:?}");
    let out = resolved(&cube);
    assert!(out.mesh.volume() > 0.0 && out.mesh.volume() < 1000.0);
}

#[test]
fn a_fold_keeps_where_the_winding_is_positive_and_drops_where_negative() {
    // The top edge at y = 10 pushed down to z = -5: every section across x
    // is the quadrilateral (0,0) (10,0) (10,-5) (0,10) in (y, z), which
    // crosses itself at y = 20/3. The upper triangle winds once, the lower
    // one minus once, so only the prism over the upper one is kept.
    let mut cube = cuboid([0.0; 3], [10.0; 3]);
    for v in [6, 7] {
        cube.positions[v].z = -5.0;
    }
    let pairs = pairs_of(&cube);
    assert!(!pairs.is_empty());
    let out = resolved(&cube);
    assert!((cube.volume() - 250.0).abs() < 1e-12);
    assert!(
        (out.mesh.volume() - 1000.0 / 3.0).abs() < 1e-9,
        "{}",
        out.mesh.volume()
    );
    assert!(out.rounding > 0.0 && out.rounding < 1e-13);
}

#[test]
fn a_nested_shell_is_dropped_and_an_inside_out_one_kept_as_a_cavity() {
    let outer = cuboid([0.0; 3], [10.0; 3]);
    let inner = cuboid([2.0; 3], [8.0; 3]);
    let nested = joined(&[&outer, &inner]);
    assert!(pairs_of(&nested).is_empty());
    let out = resolved(&nested);
    assert_eq!(out.mesh.triangles, outer.triangles);
    assert!((out.mesh.volume() - 1000.0).abs() < 1e-12);

    let hollow = joined(&[&outer, &inside_out(&inner)]);
    let out = resolved(&hollow);
    assert_eq!(out.mesh.triangles.len(), 24);
    assert!((out.mesh.volume() - 784.0).abs() < 1e-12);

    let out = resolved(&inside_out(&outer));
    assert!(out.mesh.triangles.is_empty(), "nowhere winds positively");
}

#[test]
fn an_inside_out_shell_crossing_another_is_cut_out_of_it() {
    let a = cuboid([0.0; 3], [10.0; 3]);
    let b = cuboid([5.0, 0.0, 0.0], [15.0, 10.0, 10.0]);
    let out = resolved(&joined(&[&a, &inside_out(&b)]));
    assert!((out.mesh.volume() - 500.0).abs() < 1e-9);
}

#[test]
fn cubes_sharing_a_face_touch_there_only_and_resolve_to_one_box() {
    let a = cuboid([0.0; 3], [10.0; 3]);
    let b = cuboid([10.0, 0.0, 0.0], [20.0, 10.0, 10.0]);
    let two = joined(&[&a, &b]);
    let pairs = pairs_of(&two);
    assert!(!pairs.is_empty());
    let on_contact = |t: u32| {
        two.triangles[t as usize]
            .iter()
            .any(|&v| two.positions[v as usize].x == 10.0)
    };
    for &[t, u] in &pairs {
        assert!(t < 12 && u >= 12 && on_contact(t) && on_contact(u));
    }
    let out = resolved(&two);
    assert!((out.mesh.volume() - 2000.0).abs() < 1e-12);
    assert_eq!(out.mesh.triangles.len(), 20);
    assert!(out.cut.iter().all(|&c| !c));
}

#[test]
fn an_open_sheet_folded_through_itself_is_found_and_not_resolved() {
    let mut sheet = Triangulation::new();
    sheet.positions = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(10.0, 0.0, 0.0),
        Point::new(0.0, 10.0, 0.0),
        Point::new(5.0, 5.0, 5.0),
        Point::new(2.0, 2.0, -5.0),
    ];
    sheet.triangles = vec![[0, 1, 2], [1, 3, 2], [2, 3, 4]];
    assert_eq!(pairs_of(&sheet), vec![[0, 2]]);
    let e = resolve_self_intersections(&sheet).unwrap_err();
    assert!(message(&e).contains("is open"), "{e}");
}

#[test]
fn neighbours_meeting_at_their_corners_only_are_never_reported() {
    // Two tetrahedra sharing one vertex, and a fan about a vertex whose
    // last blade passes through the first beyond the vertex they share.
    let mut tets = Triangulation::new();
    tets.positions = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(0.0, 1.0, 0.0),
        Point::new(0.0, 0.0, 1.0),
        Point::new(-1.0, 0.0, 0.0),
        Point::new(0.0, -1.0, 0.0),
        Point::new(0.0, 0.0, -1.0),
    ];
    tets.triangles = vec![
        [0, 2, 1],
        [0, 1, 3],
        [0, 3, 2],
        [1, 2, 3],
        [0, 4, 5],
        [0, 6, 4],
        [0, 5, 6],
        [4, 6, 5],
    ];
    assert!(pairs_of(&tets).is_empty());
    let mut blade = Triangulation::new();
    blade.positions = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(10.0, 0.0, 0.0),
        Point::new(0.0, 10.0, 0.0),
        Point::new(2.0, 2.0, -5.0),
        Point::new(2.0, 2.0, 5.0),
    ];
    blade.triangles = vec![[0, 1, 2], [0, 3, 4]];
    assert_eq!(pairs_of(&blade), vec![[0, 1]]);
}

#[test]
fn degenerate_and_empty_meshes_are_refused_by_name() {
    let mut flat = cuboid([0.0; 3], [1.0; 3]);
    flat.positions.push(Point::new(3.0, 3.0, 3.0));
    flat.triangles.push([0, 7, 8]);
    for e in [
        self_intersections(&flat).unwrap_err(),
        resolve_self_intersections(&flat).unwrap_err(),
    ] {
        assert!(
            message(&e).contains("triangle 12 of the input mesh is degenerate"),
            "{e}"
        );
    }
    let empty = Triangulation::new();
    assert!(self_intersections(&empty).unwrap().is_empty());
    let e = resolve_self_intersections(&empty).unwrap_err();
    assert!(message(&e).contains("no triangles"), "{e}");
}

#[test]
fn spheres_resolve_to_the_volume_the_two_operand_union_gives() {
    // Three balls crossing pairwise and all together: triple points where
    // three surfaces meet away from any edge.
    let a = sphere([0.0, 0.0, 0.0], 1.0, 10, 16, 0.31);
    let b = sphere([0.9, 0.13, 0.07], 0.8, 9, 15, 0.77);
    let c = sphere([0.4, 0.71, -0.23], 0.7, 8, 13, 1.13);
    let ab = boolean(&a, &b, BoolKind::Union).unwrap();
    let want2 = ab.mesh.volume();
    let out = resolved(&joined(&[&a, &b]));
    assert!(
        (out.mesh.volume() - want2).abs() < 1e-12,
        "{} {want2}",
        out.mesh.volume()
    );
    let want3 = boolean(&ab.mesh, &c, BoolKind::Union)
        .unwrap()
        .mesh
        .volume();
    let three = joined(&[&a, &b, &c]);
    let pairs = pairs_of(&three);
    assert!(!pairs.is_empty());
    let out = resolved(&three);
    assert!(
        (out.mesh.volume() - want3).abs() < 1e-12,
        "{} {want3}",
        out.mesh.volume()
    );
}

#[test]
fn a_shell_doubled_on_itself_resolves_to_one_copy() {
    let a = cuboid([0.0; 3], [1.0; 3]);
    let out = resolved(&joined(&[&a, &a]));
    assert!((out.mesh.volume() - 1.0).abs() < 1e-12);
    let ball = sphere([0.0; 3], 1.0, 7, 11, 0.4);
    let out = resolved(&joined(&[&ball, &ball]));
    assert!((out.mesh.volume() - ball.volume()).abs() < 1e-12);
    let out = resolved(&joined(&[&ball, &inside_out(&ball)]));
    assert!(out.mesh.triangles.is_empty());
}

#[test]
#[ignore = "heavy: a hundred thousand triangles"]
fn a_large_mesh_with_a_few_crossings_is_scanned_quickly() {
    let ball = sphere([0.0; 3], 10.0, 225, 224, 0.2);
    let probe = cuboid([9.5, -0.3, -0.2], [10.6, 0.4, 0.3]);
    let mesh = joined(&[&ball, &probe]);
    assert!(mesh.triangles.len() > 100_000);
    let start = ogeom_core::clock::Instant::now();
    let pairs = self_intersections(&mesh).unwrap();
    let took = start.elapsed();
    let first = id(ball.triangles.len());
    assert!(!pairs.is_empty() && pairs.len() < 100, "{}", pairs.len());
    assert!(pairs.iter().all(|&[t, u]| t < first && u >= first));
    // The probe's triangles against every triangle of the ball, by brute
    // force: the ball alone is embedded.
    let mut want = Vec::new();
    let tri = |t: usize, side: Side| {
        let v = mesh.triangles[t];
        Tri {
            side,
            id: id(t),
            v,
            p: v.map(|k| coords(mesh.positions[k as usize])),
        }
    };
    for t in 0..ball.triangles.len() {
        for u in ball.triangles.len()..mesh.triangles.len() {
            if !meet_at_shared_only(&tri(t, 0), &tri(u, 1)) {
                want.push([id(t), id(u)]);
            }
        }
    }
    assert_eq!(pairs, want);
    assert!(took.as_secs_f64() < 30.0, "{took:?}");
    let out = resolve_self_intersections(&mesh).unwrap();
    assert_eq!(out.crossings, pairs.len());
}

struct Rng(u64);

impl Rng {
    /// A number in `[0, 1)`.
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }
}

#[test]
fn grid_aligned_cubes_resolve_to_the_cells_they_wind_about() {
    // Cubes of side 2 at the corners of the unit cube share faces, edges
    // and corners in every combination, some wound inside out. The answer
    // is the unit cells the mesh winds about at least once; the only
    // refusal allowed is for an answer that touches itself.
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    let mut answered = 0;
    for _ in 0..60 {
        let mut parts = Vec::new();
        let mut wind = [[[0_i32; 3]; 3]; 3];
        for _ in 0..2 + (rng.next() * 4.0) as usize {
            let at = [0, 1, 2].map(|_| f64::from(u8::from(rng.next() < 0.5)));
            let outward = rng.next() < 0.75;
            let cube = cuboid(at, at.map(|x| x + 2.0));
            parts.push(if outward { cube } else { inside_out(&cube) });
            for x in 0..2 {
                for y in 0..2 {
                    for z in 0..2 {
                        let c = [x, y, z].map(|k| k as usize);
                        let cell = |k: usize| c[k] + at[k] as usize;
                        wind[cell(0)][cell(1)][cell(2)] += if outward { 1 } else { -1 };
                    }
                }
            }
        }
        let mesh = joined(&parts.iter().collect::<Vec<_>>());
        pairs_of(&mesh);
        let want = wind.iter().flatten().flatten().filter(|&&w| w >= 1).count();
        match resolve_self_intersections(&mesh) {
            Ok(_) => {
                let out = resolved(&mesh);
                assert!((out.mesh.volume() - want as f64).abs() < 1e-12, "{wind:?}");
                answered += 1;
            }
            Err(e) => assert!(message(&e).contains("not be 2-manifold"), "{e}"),
        }
    }
    assert!(answered > 40, "{answered}");
}

#[test]
fn turned_cubes_resolve_to_the_volume_of_their_union_taken_in_turn() {
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let unit = cuboid([-0.5; 3], [0.5; 3]);
    let mut answered = 0;
    for _ in 0..12 {
        let mut parts = Vec::new();
        for _ in 0..3 + (rng.next() * 2.0) as usize {
            let axis = [0, 1, 2].map(|_| rng.next() - 0.5);
            let norm = axis.iter().map(|x| x * x).sum::<f64>().sqrt();
            let [ux, uy, uz] = axis.map(|x| x / norm);
            let (s, c) = (rng.next() * std::f64::consts::TAU).sin_cos();
            let shift = [0, 1, 2].map(|_| rng.next() * 0.6 - 0.3);
            let mut cube = unit.clone();
            for p in &mut cube.positions {
                // Rodrigues' rotation about the unit axis, then the shift.
                let (x, y, z) = (p.x, p.y, p.z);
                let dot = ux * x + uy * y + uz * z;
                let cr = [uy * z - uz * y, uz * x - ux * z, ux * y - uy * x];
                let r = [x, y, z]
                    .iter()
                    .zip([ux, uy, uz])
                    .zip(cr)
                    .map(|((&v, u), w)| v * c + w * s + u * dot * (1.0 - c))
                    .collect::<Vec<_>>();
                *p = Point::new(r[0] + shift[0], r[1] + shift[1], r[2] + shift[2]);
            }
            parts.push(cube);
        }
        let mesh = joined(&parts.iter().collect::<Vec<_>>());
        assert!(!pairs_of(&mesh).is_empty());
        let out = resolved(&mesh);
        // Unions taken in turn round each result before the next, which
        // can leave a sliver that the next union refuses as crossing
        // itself; resolving does all of them at once.
        let union = parts[1..].iter().try_fold(parts[0].clone(), |acc, p| {
            boolean(&acc, p, BoolKind::Union).map(|r| r.mesh)
        });
        if let Ok(union) = union {
            let want = union.volume();
            assert!(
                (out.mesh.volume() - want).abs() < 1e-12,
                "{} {want}",
                out.mesh.volume()
            );
            answered += 1;
        }
    }
    assert!(answered >= 8, "{answered}");
}

#[test]
fn a_flat_triangle_is_flipped_away_without_moving_a_vertex() {
    // The cube's bottom triangle 0 2 3 split at the middle m of its side
    // 0-3, with the flat triangle m 3 0 closing the gap against 0 3 1.
    let mut mesh = cuboid([0.0; 3], [10.0; 3]);
    mesh.positions.push(Point::new(5.0, 5.0, 0.0));
    let m = 8;
    mesh.triangles[0] = [0, 2, m];
    mesh.triangles.push([m, 2, 3]);
    mesh.triangles.push([m, 3, 0]);
    let mut cut = vec![false; mesh.triangles.len()];
    let positions = mesh.positions.clone();
    flip_needles(&mut mesh, &mut cut).unwrap();
    assert_eq!(mesh.positions, positions);
    assert!((mesh.volume() - 1000.0).abs() < 1e-12);
    let solid = Operand::new(&mesh, 0, "flipped").unwrap();
    solid.check_wound("flipped", None, 0).unwrap();
    assert!(self_intersections(&mesh).unwrap().is_empty());
    assert_eq!(cut.iter().filter(|&&c| c).count(), 2);
    assert!(mesh.triangles.contains(&[0, m, 1]) && mesh.triangles.contains(&[m, 3, 1]));
}

/// The point of triangle `abc` nearest `p`, by the regions of its plane
/// (C. Ericson, *Real-Time Collision Detection*, 5.1.5).
fn nearest_on_triangle(p: P3, [a, b, c]: [P3; 3]) -> P3 {
    let sub = |x: P3, y: P3| [x[0] - y[0], x[1] - y[1], x[2] - y[2]];
    let dot = |x: P3, y: P3| x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    let along = |o: P3, d: P3, t: f64| [o[0] + t * d[0], o[1] + t * d[1], o[2] + t * d[2]];
    let (ab, ac, ap) = (sub(b, a), sub(c, a), sub(p, a));
    let (d1, d2) = (dot(ab, ap), dot(ac, ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = sub(p, b);
    let (d3, d4) = (dot(ab, bp), dot(ac, bp));
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return along(a, ab, d1 / (d1 - d3));
    }
    let cp = sub(p, c);
    let (d5, d6) = (dot(ab, cp), dot(ac, cp));
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return along(a, ac, d2 / (d2 - d6));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        return along(b, sub(c, b), (d4 - d3) / ((d4 - d3) + (d5 - d6)));
    }
    let denom = 1.0 / (va + vb + vc);
    let (v, w) = (vb * denom, vc * denom);
    along(along(a, ab, v), ac, w)
}

/// How far `p` lies from the surface of `mesh`, by every triangle.
fn distance_to(mesh: &Triangulation, p: P3) -> f64 {
    mesh.triangles
        .iter()
        .map(|t| {
            let q = nearest_on_triangle(p, t.map(|v| coords(mesh.positions[v as usize])));
            ((p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2) + (p[2] - q[2]).powi(2)).sqrt()
        })
        .fold(f64::INFINITY, f64::min)
}

fn area(mesh: &Triangulation) -> f64 {
    mesh.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| mesh.positions[v as usize]);
            (b - a).cross(c - a).magnitude() / 2.0
        })
        .sum()
}

/// What rounding to `quantum` did: the farthest a vertex of the result
/// lies from the input surface, the bound on it, and the change in volume
/// from the `f64` result with the first-order bound on that.
#[derive(Debug, Default, Clone, Copy)]
struct Measured {
    moved: f64,
    bound: f64,
    volume: f64,
    volume_bound: f64,
}

/// The answer under `quantum`, checked as [`resolved`] checks and then on
/// the rounded coordinates: every coordinate on the grid, no crossing, no
/// flat triangle, every shell outward; every vertex within
/// [`Resolved::rounding`] of the input surface, and that within the
/// documented bound.
fn snapped(mesh: &Triangulation, quantum: Quantum) -> (Resolved, Measured) {
    let out = resolve_self_intersections_with(mesh, quantum).unwrap();
    let exact = resolve_self_intersections(mesh).unwrap();
    assert_eq!(out.cut.len(), out.mesh.triangles.len());
    for p in &out.mesh.positions {
        for x in coords(*p) {
            assert_eq!(quantum.round(x), Some(x), "{x} is off the grid");
        }
    }
    let largest = mesh
        .positions
        .iter()
        .chain(&exact.mesh.positions)
        .flat_map(|p| coords(*p))
        .fold(0.0_f64, |m, x| m.max(x.abs()));
    // Five roundings onto the grid at most, plus the f64 rounding of each
    // pass's new vertices, a few units in the last place of `largest`.
    let bound = 5.0 * (3.0_f64.sqrt() * quantum.half_step(largest) + 64.0 * f64::EPSILON * largest);
    assert!(out.rounding <= bound, "{} > {bound}", out.rounding);
    let moved = out
        .mesh
        .positions
        .iter()
        .map(|p| distance_to(mesh, coords(*p)))
        .fold(0.0_f64, f64::max);
    assert!(
        moved <= out.rounding * (1.0 + 1e-9) + 4.0 * f64::EPSILON * largest,
        "{moved} > {}",
        out.rounding
    );
    let volume = (out.mesh.volume() - exact.mesh.volume()).abs();
    let volume_bound = area(&exact.mesh).max(area(&out.mesh)) * out.rounding;
    assert!(volume <= volume_bound, "{volume} > {volume_bound}");
    if !out.mesh.triangles.is_empty() {
        let solid = Operand::new(&out.mesh, 0, "rounded").unwrap();
        solid.check_wound("rounded", None, 0).unwrap();
        assert_eq!(pairs_of(&out.mesh), Vec::<[u32; 2]>::new());
    }
    (
        out,
        Measured {
            moved,
            bound,
            volume,
            volume_bound,
        },
    )
}

/// `mesh` turned by `angle` about the axis `axis`, through `centre`.
fn turned(mesh: &Triangulation, axis: P3, angle: f64, centre: P3) -> Triangulation {
    let norm = axis.iter().map(|x| x * x).sum::<f64>().sqrt();
    let [ux, uy, uz] = axis.map(|x| x / norm);
    let (s, c) = angle.sin_cos();
    let mut out = mesh.clone();
    for p in &mut out.positions {
        let (x, y, z) = (p.x - centre[0], p.y - centre[1], p.z - centre[2]);
        let dot = ux * x + uy * y + uz * z;
        let cr = [uy * z - uz * y, uz * x - ux * z, ux * y - uy * x];
        let r: Vec<f64> = [x, y, z]
            .iter()
            .zip([ux, uy, uz])
            .zip(cr)
            .map(|((&v, u), w)| v * c + w * s + u * dot * (1.0 - c))
            .collect();
        *p = Point::new(r[0] + centre[0], r[1] + centre[1], r[2] + centre[2]);
    }
    out
}

/// `mesh` moved by `by`.
fn shifted(mesh: &Triangulation, by: P3) -> Triangulation {
    let mut out = mesh.clone();
    for p in &mut out.positions {
        *p = Point::new(p.x + by[0], p.y + by[1], p.z + by[2]);
    }
    out
}

#[test]
fn double_precision_is_the_plain_resolution_bit_for_bit() {
    let a = cuboid([0.0; 3], [10.0; 3]);
    let b = cuboid([5.0, 3.0, 4.0], [15.0, 12.0, 13.0]);
    let ball = sphere([0.0; 3], 1.0, 10, 16, 0.31);
    let other = sphere([0.9, 0.13, 0.07], 0.8, 9, 15, 0.77);
    for mesh in [joined(&[&a, &b]), joined(&[&ball, &other])] {
        assert_eq!(
            resolve_self_intersections_with(&mesh, Quantum::F64).unwrap(),
            resolve_self_intersections(&mesh).unwrap()
        );
    }
}

#[test]
fn the_kernel_fixtures_rounded_to_single_precision_never_cross() {
    let a = cuboid([0.0; 3], [10.0; 3]);
    let b = cuboid([5.0, 3.0, 4.0], [15.0, 12.0, 13.0]);
    let ball = sphere([0.0; 3], 1.0, 10, 16, 0.31);
    let b2 = sphere([0.9, 0.13, 0.07], 0.8, 9, 15, 0.77);
    let c2 = sphere([0.4, 0.71, -0.23], 0.7, 8, 13, 1.13);
    let mut meshes = vec![
        joined(&[&a, &b]),
        joined(&[&a, &a]),
        joined(&[&ball, &b2]),
        joined(&[&ball, &b2, &c2]),
        joined(&[&a, &cuboid([2.0; 3], [8.0; 3])]),
        joined(&[&a, &inside_out(&cuboid([2.0; 3], [8.0; 3]))]),
    ];
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    let unit = cuboid([-0.5; 3], [0.5; 3]);
    for _ in 0..12 {
        let parts: Vec<Triangulation> = (0..3 + (rng.next() * 2.0) as usize)
            .map(|_| {
                let axis = [0, 1, 2].map(|_| rng.next() - 0.5);
                let shift = [0, 1, 2].map(|_| rng.next() * 0.6 - 0.3);
                shifted(
                    &turned(&unit, axis, rng.next() * std::f64::consts::TAU, [0.0; 3]),
                    shift,
                )
            })
            .collect();
        meshes.push(joined(&parts.iter().collect::<Vec<_>>()));
    }
    let (naive, refused) = report("kernel fixtures", &meshes, Quantum::F32);
    assert_eq!(refused, 0);
    assert!(naive > 0, "no fixture needed more than rounding");
}

/// Each mesh resolved under `quantum` and checked by [`snapped`]; prints
/// the largest of each measure over them. Answers how many of the `f64`
/// results rounded naively, every coordinate to the grid, would have
/// crossed themselves or held a flat triangle, and how many were refused,
/// each by name.
fn report(what: &str, meshes: &[Triangulation], quantum: Quantum) -> (usize, usize) {
    let mut worst = Measured::default();
    let (mut naive, mut refused) = (0, 0);
    for mesh in meshes {
        if let Err(e) = resolve_self_intersections_with(mesh, quantum) {
            assert!(
                message(&e).contains(&format!("rounding the result to {quantum}")),
                "{e}"
            );
            refused += 1;
            continue;
        }
        let (_, m) = snapped(mesh, quantum);
        worst.moved = worst.moved.max(m.moved);
        worst.bound = worst.bound.max(m.bound);
        worst.volume = worst.volume.max(m.volume);
        worst.volume_bound = worst.volume_bound.max(m.volume_bound);
        let mut plain = resolve_self_intersections(mesh).unwrap().mesh;
        for p in &mut plain.positions {
            let [x, y, z] = coords(*p).map(|x| quantum.round(x).unwrap());
            *p = Point::new(x, y, z);
        }
        if !self_intersections(&plain).is_ok_and(|pairs| pairs.is_empty()) {
            naive += 1;
        }
    }
    eprintln!(
        "{what}: {} meshes, {refused} refused, {naive} cross when rounded naively; vertex moved {:.3e} (bound {:.3e}), volume changed {:.3e} (bound {:.3e})",
        meshes.len(),
        worst.moved,
        worst.bound,
        worst.volume,
        worst.volume_bound
    );
    (naive, refused)
}

#[test]
fn boxes_whose_faces_lie_closer_than_a_single_step_round_to_one_solid() {
    // A box against a face of another, overlapping it, apart from it, or
    // turned on it by less than a step of `f32` at that size: rounded,
    // the near faces fall on one plane.
    let a = cuboid([0.0; 3], [1.0; 3]);
    let mut meshes = Vec::new();
    for gap in [-1e-9, -3e-8, 1e-9, 3e-8, 1e-12, -1e-12] {
        let b = cuboid([1.0 + gap, 0.25, 0.3], [2.0, 0.75, 0.8]);
        meshes.push(joined(&[&a, &b]));
        let flush = cuboid([1.0 + gap, 0.0, 0.0], [2.0, 1.0, 1.0]);
        meshes.push(joined(&[&a, &flush]));
        let over = cuboid([0.5, 0.2, 1.0 + gap], [1.5, 0.9, 1.7]);
        meshes.push(joined(&[&a, &over]));
    }
    for angle in [1e-9, 3e-8, -2e-8] {
        let b = cuboid([0.9, 0.25, 0.3], [2.0, 0.75, 0.8]);
        meshes.push(joined(&[
            &a,
            &turned(&b, [0.0, 0.0, 1.0], angle, [1.0, 0.5, 0.5]),
        ]));
    }
    let (naive, refused) = report("near-coincident boxes", &meshes, Quantum::F32);
    assert_eq!(refused, 0);
    assert!(naive > meshes.len() / 2, "{naive}");
}

#[test]
fn shells_turned_on_themselves_by_a_hair_leave_slivers_that_round_clean() {
    // A shell crossed with itself turned by a few steps of `f32`: every
    // pair of faces meets at a vanishing angle, so the pieces are slivers
    // a step or two across.
    let ball = sphere([0.0; 3], 1.0, 8, 12, 0.4);
    let cube = cuboid([-1.0; 3], [1.0; 3]);
    let mut meshes = Vec::new();
    for angle in [2e-7, 1e-6, 5e-6] {
        for axis in [[1.0, 0.3, 0.2], [0.1, 1.0, -0.4]] {
            meshes.push(joined(&[
                &ball,
                &turned(&ball, axis, angle, [0.01, 0.0, 0.0]),
            ]));
            meshes.push(joined(&[&cube, &turned(&cube, axis, angle, [0.0; 3])]));
        }
    }
    // Where the rounded pieces touch along an edge or at a point, the
    // answer would not be a 2-manifold and is refused by name.
    let (naive, refused) = report("hair-turned shells", &meshes, Quantum::F32);
    assert!(naive > 0, "no sliver here needed more than rounding");
    assert!(refused <= 4, "{refused}");
}

#[test]
fn randomly_jittered_overlaps_round_clean_to_single_precision() {
    // Boxes and balls jittered by a few steps of `f32`, vertex by vertex,
    // crossing a copy of themselves jittered differently.
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    let mut meshes = Vec::new();
    for k in 0..16 {
        let base = if k % 2 == 0 {
            cuboid([-0.7, -0.4, -0.5], [0.7, 0.4, 0.5])
        } else {
            sphere([0.0; 3], 0.8, 7, 10, 0.2 + f64::from(k))
        };
        let jitter = |m: &Triangulation, rng: &mut Rng| {
            let mut out = m.clone();
            for p in &mut out.positions {
                let d = [0, 1, 2].map(|_| (rng.next() - 0.5) * 4e-7);
                *p = Point::new(p.x + d[0], p.y + d[1], p.z + d[2]);
            }
            out
        };
        let (a, b) = (jitter(&base, &mut rng), jitter(&base, &mut rng));
        meshes.push(joined(&[&a, &b]));
    }
    let (naive, refused) = report("jittered copies", &meshes, Quantum::F32);
    assert!(naive > meshes.len() / 2, "{naive}");
    assert!(refused <= 1, "{refused}");
}

#[test]
fn a_coarse_grid_holds_the_result_and_a_bad_spacing_is_refused() {
    let a = cuboid([0.0; 3], [10.0; 3]);
    let b = turned(
        &cuboid([5.0, 3.0, 4.0], [15.0, 12.0, 13.0]),
        [1.0, 2.0, 3.0],
        0.3,
        [0.0; 3],
    );
    report("grid 1/64", &[joined(&[&a, &b])], Quantum::Grid(1.0 / 64.0));
    for h in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let e = resolve_self_intersections_with(&a, Quantum::Grid(h)).unwrap_err();
        assert!(message(&e).contains("positive and finite"), "{e}");
    }
    let far = cuboid([1e39; 3], [2e39; 3]);
    let e = resolve_self_intersections_with(&far, Quantum::F32).unwrap_err();
    assert!(message(&e).contains("single precision"), "{e}");
}
