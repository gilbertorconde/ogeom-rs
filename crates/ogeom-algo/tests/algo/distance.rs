//! Signed distances from points to meshes and shapes, against closed forms:
//! offsets of convex meshes and of boxes, spheres and cylinders read their
//! offset, points round edges and corners read their side, a thin plate
//! scanned from both sides reads each side against its own face, and the
//! mesh answers equal brute force with and without a facing filter.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::FastMap;

use ogeom_algo::distance::{
    Closest, DistanceOptions, MeshFeature, ShapeFeature, to_mesh, to_shape,
};
use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_mesh::{Bvh, Deflection};
use ogeom_topo::{Model, Triangulation};

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

    fn below(&mut self, n: usize) -> usize {
        usize::try_from(self.next() % n as u64).unwrap()
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

    /// A distance spread over many scales, from a micron to a centimetre.
    fn scale(&mut self) -> f64 {
        10_f64.powf(self.range(-3.0, 1.0))
    }
}

fn unit(v: Vector) -> Vector {
    v / v.magnitude()
}

fn tilted() -> Frame {
    Frame::new(
        Point::new(3.0, -2.0, 5.0),
        Direction::new(Vector::new(0.3, -0.4, 1.0), T).unwrap(),
        Direction::new(Vector::new(1.0, 0.5, 0.2), T).unwrap(),
        T,
    )
    .unwrap()
}

/// The exact signed distance to a box spanning `size` from the origin.
fn box_distance(p: Point, size: [f64; 3]) -> f64 {
    let q: Vec<f64> = (0..3)
        .map(|k| {
            let c = [p.x, p.y, p.z][k];
            (c - size[k] / 2.0).abs() - size[k] / 2.0
        })
        .collect();
    let outside = q.iter().map(|c| c.max(0.0).powi(2)).sum::<f64>().sqrt();
    let inside = q[0].max(q[1]).max(q[2]).min(0.0);
    outside + inside
}

/// A point near a box spanning `size`: near a face, an edge or a corner,
/// at a distance over many scales, inside or out.
fn near_box(rng: &mut Seeded, size: [f64; 3], kind: usize) -> Point {
    // How many of the three coordinates sit on a side: one for a face, two
    // for an edge, three for a corner.
    let pinned = 1 + kind % 3;
    let mut c = [0.0; 3];
    let first = rng.below(3);
    for k in 0..3 {
        let axis = (first + k) % 3;
        c[axis] = if k < pinned {
            let side = if rng.below(2) == 0 { 0.0 } else { size[axis] };
            side + rng.scale() * rng.range(-1.0, 1.0) * 0.1
        } else {
            rng.range(0.0, size[axis])
        };
    }
    Point::new(c[0], c[1], c[2])
}

/// An icosahedron subdivided `levels` times onto a sphere of `radius`,
/// wound outward.
fn icosphere(radius: f64, levels: u32) -> Triangulation {
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
    let mut mesh = Triangulation::new();
    mesh.positions = positions
        .iter()
        .map(|v| Point::ORIGIN + *v * radius)
        .collect();
    mesh.triangles = triangles;
    mesh
}

/// A box spanning `size` from the origin, each side an `n` by `n` grid of
/// quads, wound outward; the sides do not share their corners' indices,
/// only their positions.
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
                    // (u, v, axis) is right-handed, so a, b, c turns about
                    // +axis: outward on the high side.
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

/// A closed, non-convex surface: a torus with its vertices moved off it
/// radially by up to `lump` of their radius, wound outward.
fn torus(around: u32, across: u32, lump: f64, rng: &mut Seeded) -> Triangulation {
    let mut mesh = Triangulation::new();
    for i in 0..around {
        for j in 0..across {
            let u = std::f64::consts::TAU * f64::from(i) / f64::from(around);
            let v = std::f64::consts::TAU * f64::from(j) / f64::from(across);
            let lump = 1.0 + lump * rng.range(-1.0, 1.0);
            let r = (2.0 + 0.7 * v.cos()) * lump;
            mesh.positions
                .push(Point::new(r * u.cos(), r * u.sin(), 0.7 * v.sin() * lump));
        }
    }
    let at = |i: u32, j: u32| (i % around) * across + (j % across);
    for i in 0..around {
        for j in 0..across {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            mesh.triangles.push([a, b, c]);
            mesh.triangles.push([a, c, d]);
        }
    }
    mesh
}

/// The winding number of a closed mesh about `p`: the solid angles its
/// triangles subtend there, over 4π.
fn winding(mesh: &Triangulation, p: Point) -> f64 {
    let mut sum = 0.0;
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|k| mesh.positions[k as usize] - p);
        let (la, lb, lc) = (a.magnitude(), b.magnitude(), c.magnitude());
        let numerator = a.dot(b.cross(c));
        let denominator = la * lb * lc + a.dot(b) * lc + b.dot(c) * la + c.dot(a) * lb;
        sum += 2.0 * numerator.atan2(denominator);
    }
    sum / (4.0 * std::f64::consts::PI)
}

fn unit_normal(t: [Point; 3]) -> Vector {
    let n = (t[1] - t[0]).cross(t[2] - t[0]);
    let m = n.magnitude();
    if m > 0.0 && m.is_finite() {
        n / m
    } else {
        Vector::ZERO
    }
}

/// Every vertex of a convex mesh moved out by `d` along its vertex normal,
/// every face's centroid along its normal and every edge's middle along
/// its two faces' bisector: each is `d` from the mesh, nearest the feature
/// it left.
fn offsets_of_convex(mesh: &Triangulation, d: f64) -> Vec<(Point, f64)> {
    let mut out = Vec::new();
    let mut at_vertex = vec![Vector::ZERO; mesh.positions.len()];
    let mut at_edge: FastMap<(u32, u32), Vector> = FastMap::default();
    for t in &mesh.triangles {
        let corners = t.map(|k| mesh.positions[k as usize]);
        let n = unit_normal(corners);
        let centre = Point::ORIGIN
            + ((corners[0] - Point::ORIGIN)
                + (corners[1] - Point::ORIGIN)
                + (corners[2] - Point::ORIGIN))
                / 3.0;
        out.push((centre + n * d, d));
        for k in 0..3 {
            at_vertex[t[k] as usize] += n;
            let (a, b) = (t[k], t[(k + 1) % 3]);
            *at_edge.entry((a.min(b), a.max(b))).or_insert(Vector::ZERO) += n;
        }
    }
    for (p, n) in mesh.positions.iter().zip(&at_vertex) {
        out.push((*p + unit(*n) * d, d));
    }
    for (&(a, b), n) in &at_edge {
        let middle = mesh.positions[a as usize]
            + (mesh.positions[b as usize] - mesh.positions[a as usize]) * 0.5;
        out.push((middle + unit(*n) * d, d));
    }
    out
}

fn distances<F>(found: &[Option<Closest<F>>]) -> Vec<f64> {
    found.iter().map(|c| c.as_ref().unwrap().distance).collect()
}

/// A convex mesh offset outward by `d`, at its vertices, edges' middles
/// and faces' centroids, reads `+d` against the original, within 1e-6: a
/// sphere's and a box's, the box written with its sides' corners repeated.
/// Moved inward, a face's centroid reads `-d`.
#[test]
fn a_convex_mesh_offset_outward_reads_the_offset() {
    for (mesh, ds) in [
        (icosphere(10.0, 4), [0.01, 0.5, 3.0]),
        (box_mesh([4.0, 3.0, 2.0], 6), [1e-4, 0.25, 2.0]),
    ] {
        let index = Bvh::build(&mesh).unwrap();
        for d in ds {
            let (points, wanted): (Vec<Point>, Vec<f64>) =
                offsets_of_convex(&mesh, d).into_iter().unzip();
            let found = to_mesh(&index, &points, None, &DistanceOptions::default()).unwrap();
            for ((got, want), p) in distances(&found).iter().zip(&wanted).zip(&points) {
                assert!((got - want).abs() <= 1e-6, "{p:?}: {got} for {want}");
            }
            // Inward from the centroids, by less than any face is from
            // its neighbours' planes there.
            let inward: Vec<Point> = points[..mesh.triangles.len()]
                .iter()
                .zip(&mesh.triangles)
                .map(|(p, t)| *p - unit_normal(t.map(|k| mesh.positions[k as usize])) * (d + 0.01))
                .collect();
            let found = to_mesh(&index, &inward, None, &DistanceOptions::default()).unwrap();
            for got in distances(&found) {
                assert!((got + 0.01).abs() <= 1e-6, "{got}");
            }
        }
    }
}

/// Points all round a box's faces, edges and corners, from a micron to a
/// centimetre off them, inside and out, read the box's exact signed
/// distance against an exact box mesh: no sign flips where the nearest
/// point is on an edge or a corner.
#[test]
fn points_round_a_box_mesh_s_edges_and_corners_read_their_side() {
    let size = [4.0, 3.0, 2.0];
    let mesh = box_mesh(size, 3);
    let index = Bvh::build(&mesh).unwrap();
    let mut rng = Seeded(0x0b0c_5eed);
    let points: Vec<Point> = (0..60_000).map(|i| near_box(&mut rng, size, i)).collect();
    let found = to_mesh(&index, &points, None, &DistanceOptions::default()).unwrap();
    let mut on = [0_usize; 3];
    for (p, c) in points.iter().zip(&found) {
        let c = c.as_ref().unwrap();
        let want = box_distance(*p, size);
        assert!(
            (c.distance - want).abs() <= 1e-12,
            "{p:?}: {} for {want}",
            c.distance
        );
        on[match c.feature {
            MeshFeature::Face(_) => 0,
            MeshFeature::Edge(..) => 1,
            MeshFeature::Vertex(..) => 2,
        }] += 1;
    }
    // Each feature was nearest often enough to have been tested.
    assert!(on.iter().all(|&n| n > 2000), "{on:?}");
}

/// Points round the vertices and edges of a lumpy torus, convex and
/// saddle-shaped alike, take the side its winding number puts them on:
/// written with shared corners, and with every triangle's corners its own.
#[test]
fn points_round_a_non_convex_mesh_s_vertices_and_edges_read_their_side() {
    let mut rng = Seeded(0x7041_0005);
    let mesh = torus(40, 16, 0.08, &mut rng);
    let mut apart = Triangulation::new();
    for t in &mesh.triangles {
        let first = u32::try_from(apart.positions.len()).unwrap();
        apart
            .positions
            .extend(t.map(|k| mesh.positions[k as usize]));
        apart.triangles.push([first, first + 1, first + 2]);
    }
    let mut points = Vec::new();
    for _ in 0..20_000 {
        let t = mesh.triangles[rng.below(mesh.triangles.len())];
        let (a, b) = (mesh.positions[t[0] as usize], mesh.positions[t[1] as usize]);
        let at = if rng.below(2) == 0 {
            a
        } else {
            a + (b - a) * rng.range(0.0, 1.0)
        };
        points.push(at + rng.direction() * rng.scale() * 0.01);
    }
    for written in [&mesh, &apart] {
        let index = Bvh::build(written).unwrap();
        let found = to_mesh(&index, &points, None, &DistanceOptions::default()).unwrap();
        let mut sides = [0_usize; 2];
        for (p, c) in points.iter().zip(&found) {
            let c = c.as_ref().unwrap();
            let inside = winding(&mesh, *p) > 0.5;
            assert_eq!(c.distance < 0.0, inside, "{p:?}: {c:?}");
            sides[usize::from(inside)] += 1;
        }
        assert!(sides.iter().all(|&n| n > 5000), "{sides:?}");
    }
}

/// The answer on a random mesh, with and without a facing filter, is brute
/// force's: the nearest of the triangles kept, ties to the lowest, each
/// triangle's own distance as the hierarchy measures it alone. The same at
/// any thread count.
#[test]
fn mesh_answers_equal_brute_force_with_and_without_facing() {
    let mut rng = Seeded(0x00b2_07e0);
    for round in 0..4 {
        let mesh = if round % 2 == 0 {
            torus(16 + u32::try_from(rng.below(8)).unwrap(), 8, 0.08, &mut rng)
        } else {
            let mut soup = Triangulation::new();
            for _ in 0..300 {
                let a = Point::ORIGIN + rng.direction() * rng.range(0.0, 2.0);
                let first = u32::try_from(soup.positions.len()).unwrap();
                soup.positions
                    .extend([a, a + rng.direction() * 0.4, a + rng.direction() * 0.4]);
                soup.triangles.push([first, first + 1, first + 2]);
            }
            soup
        };
        let index = Bvh::build(&mesh).unwrap();
        let corners: Vec<[Point; 3]> = (0..mesh.triangles.len())
            .map(|t| index.triangle(u32::try_from(t).unwrap()).unwrap())
            .collect();
        let points: Vec<Point> = (0..1500)
            .map(|_| Point::ORIGIN + rng.direction() * rng.range(0.0, 3.5))
            .collect();
        let normals: Vec<Vector> = (0..points.len()).map(|_| rng.direction() * 2.0).collect();
        for limit in [
            None,
            Some(0.4),
            Some(std::f64::consts::FRAC_PI_2),
            Some(2.5),
        ] {
            let options = DistanceOptions {
                max: if round < 2 { f64::INFINITY } else { 1.2 },
                facing_limit: limit,
                ..DistanceOptions::default()
            };
            let found = to_mesh(&index, &points, Some(&normals), &options).unwrap();
            for (i, p) in points.iter().enumerate() {
                let m = unit(normals[i]);
                let mut best: Option<(f64, f64, u32)> = None;
                for (t, triangle) in corners.iter().enumerate() {
                    let t = u32::try_from(t).unwrap();
                    if let Some(limit) = limit {
                        let n = unit_normal(*triangle);
                        if n == Vector::ZERO || n.dot(m) < limit.cos() {
                            continue;
                        }
                    }
                    let Some(hit) = index.nearest_where(*p, options.max, |x| x == t) else {
                        continue;
                    };
                    // Compared squared, summed by axis, as the hierarchy
                    // compares them.
                    let v = hit.point - *p;
                    let d2 = v.x * v.x + v.y * v.y + v.z * v.z;
                    if best.is_none_or(|(_, held, _)| d2 < held) {
                        best = Some((hit.distance, d2, t));
                    }
                }
                let got = found[i].as_ref();
                assert_eq!(
                    got.map(|c| c.distance.abs()),
                    best.map(|b| b.0),
                    "{p:?} {limit:?}"
                );
                if let (Some(c), Some((_, _, t))) = (got, best) {
                    let held = match c.feature {
                        MeshFeature::Face(x)
                        | MeshFeature::Edge(x, _)
                        | MeshFeature::Vertex(x, _) => x,
                    };
                    assert_eq!(held, t);
                }
            }
            for threads in [1, 8] {
                ogeom_core::parallel::set_threads(threads);
                assert_eq!(
                    to_mesh(&index, &points, Some(&normals), &options).unwrap(),
                    found
                );
            }
            ogeom_core::parallel::set_threads(0);
        }
    }
}

/// Points about a box in a tilted frame, on its exact offsets out and in,
/// round its faces, edges and corners, from a micron to a centimetre off
/// them, read the box's exact signed distance within 1e-6 against the
/// B-rep's own faces: no sign flips at the edges or the corners.
#[test]
fn points_round_a_box_read_its_exact_signed_distance() {
    let size = [4.0, 3.0, 2.0];
    let frame = tilted();
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, frame, (size[0], size[1], size[2]), T)
        .unwrap()
        .shape;
    let mut rng = Seeded(0x0b0c_b4e9);
    let local: Vec<Point> = (0..6000).map(|i| near_box(&mut rng, size, i)).collect();
    let points: Vec<Point> = local.iter().map(|&p| frame.to_world(p)).collect();
    let found = to_shape(
        &model,
        &solid,
        &points,
        None,
        &DistanceOptions::default(),
        T,
    )
    .unwrap();
    let mut on = [0_usize; 3];
    for (p, c) in local.iter().zip(&found) {
        let c = c.as_ref().unwrap();
        let want = box_distance(*p, size);
        assert!(
            (c.distance - want).abs() <= 1e-6,
            "{p:?}: {} for {want}",
            c.distance
        );
        on[match c.feature {
            ShapeFeature::Face(_) => 0,
            ShapeFeature::Edge(_) => 1,
            ShapeFeature::Vertex(_) => 2,
        }] += 1;
    }
    assert!(on.iter().all(|&n| n > 100), "{on:?}");
}

/// Points on a sphere's offsets out and in, all round it and densely at
/// its seam and its poles, read the offset within 1e-6.
#[test]
fn points_on_a_sphere_s_offsets_read_the_offset() {
    let radius = 10.0;
    let frame = tilted();
    let mut model = Model::new();
    let solid = ogeom_algo::make_sphere(&mut model, frame, radius, T)
        .unwrap()
        .shape;
    let mut rng = Seeded(0x5_e4e7e);
    let mut points = Vec::new();
    let mut wanted = Vec::new();
    for i in 0..6000 {
        // Anywhere, beside the seam (the frame's x half-plane), or by a
        // pole.
        let direction = match i % 3 {
            0 => rng.direction(),
            1 => unit(Vector::new(
                1.0,
                rng.scale() * rng.range(-1e-3, 1e-3),
                rng.range(-1.0, 1.0),
            )),
            _ => unit(Vector::new(
                rng.scale() * rng.range(-1e-3, 1e-3),
                rng.scale() * rng.range(-1e-3, 1e-3),
                if rng.below(2) == 0 { 1.0 } else { -1.0 },
            )),
        };
        let d = rng.scale() * rng.range(-1.0, 1.0);
        points.push(frame.to_world(Point::ORIGIN + direction * (radius + d)));
        wanted.push(d);
    }
    let found = to_shape(
        &model,
        &solid,
        &points,
        None,
        &DistanceOptions::default(),
        T,
    )
    .unwrap();
    for ((p, c), want) in points.iter().zip(&found).zip(&wanted) {
        let got = c.as_ref().unwrap().distance;
        assert!((got - want).abs() <= 1e-6, "{p:?}: {got} for {want}");
    }
}

/// Points round a cylinder, near its circular edges where the trim's rings
/// are polylines, read its exact signed distance within 1e-6.
#[test]
fn points_round_a_cylinder_s_circular_edges_read_its_exact_signed_distance() {
    let (radius, height) = (5.0, 4.0);
    let frame = tilted();
    let mut model = Model::new();
    let solid = ogeom_algo::make_cylinder(&mut model, frame, radius, height, T)
        .unwrap()
        .shape;
    let mut rng = Seeded(0x00c1_11d5);
    let mut local = Vec::new();
    for _ in 0..4000 {
        let angle = rng.range(0.0, std::f64::consts::TAU);
        let r = radius + rng.scale() * rng.range(-0.1, 0.1);
        let z = if rng.below(2) == 0 { 0.0 } else { height } + rng.scale() * rng.range(-0.1, 0.1);
        local.push(Point::new(r * angle.cos(), r * angle.sin(), z));
    }
    let points: Vec<Point> = local.iter().map(|&p| frame.to_world(p)).collect();
    let found = to_shape(
        &model,
        &solid,
        &points,
        None,
        &DistanceOptions::default(),
        T,
    )
    .unwrap();
    for (p, c) in local.iter().zip(&found) {
        let across = (
            p.x.hypot(p.y) - radius,
            (p.z - height / 2.0).abs() - height / 2.0,
        );
        let want = across.0.max(0.0).hypot(across.1.max(0.0)) + across.0.max(across.1).min(0.0);
        let got = c.as_ref().unwrap().distance;
        assert!((got - want).abs() <= 1e-6, "{p:?}: {got} for {want}");
    }
}

/// A plate a millimetre thick, scanned from both sides with a noise that
/// puts some of each side's points nearer the other face: with a facing
/// limit of 90°, each point reads against the face its own side's scan
/// faces, as a mesh and as a B-rep; without it, the far ones read the
/// other face.
#[test]
fn a_thin_plate_scanned_on_both_sides_reads_each_side_against_its_own_face() {
    let mut model = Model::new();
    let plate = ogeom_algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 1.0), T)
        .unwrap()
        .shape;
    let mut rng = Seeded(0x91a7e);
    let mut points = Vec::new();
    let mut normals = Vec::new();
    let mut wanted = Vec::new();
    let mut owns = Vec::new();
    for i in 0..4000 {
        let (x, y) = (rng.range(5.0, 35.0), rng.range(5.0, 25.0));
        let off = rng.range(-0.75, 0.3);
        let tilt = Vector::new(rng.range(-0.3, 0.3), rng.range(-0.3, 0.3), 0.0);
        if i % 2 == 0 {
            points.push(Point::new(x, y, 1.0 + off));
            normals.push(Vector::Z + tilt);
            owns.push(1.0);
        } else {
            points.push(Point::new(x, y, -off));
            normals.push(-Vector::Z + tilt);
            owns.push(0.0);
        }
        wanted.push(off);
    }
    let facing = DistanceOptions {
        facing_limit: Some(90_f64.to_radians()),
        ..DistanceOptions::default()
    };
    let index = Bvh::of_shape(&model, &plate, Deflection::default(), T).unwrap();
    let mesh = to_mesh(&index, &points, Some(&normals), &facing).unwrap();
    let brep = to_shape(&model, &plate, &points, Some(&normals), &facing, T).unwrap();
    for ((((p, want), own), m), b) in points.iter().zip(&wanted).zip(&owns).zip(&mesh).zip(&brep) {
        for c in [m.as_ref().unwrap().point, b.as_ref().unwrap().point] {
            assert!((c.z - own).abs() <= 1e-12, "{p:?} read {c:?}");
        }
        assert!((m.as_ref().unwrap().distance - want).abs() <= 1e-9);
        assert!((b.as_ref().unwrap().distance - want).abs() <= 1e-6);
    }
    let unfiltered = to_shape(
        &model,
        &plate,
        &points,
        Some(&normals),
        &DistanceOptions::default(),
        T,
    )
    .unwrap();
    let crossed = owns
        .iter()
        .zip(&unfiltered)
        .filter(|(own, c)| (c.as_ref().unwrap().point.z - **own).abs() > 0.5)
        .count();
    assert!(crossed > 200, "{crossed}");
}

/// Beyond `max` there is no answer; a facing limit with no normals, a
/// normal per point missing, or a limit past π is refused; a point with
/// no normal of its own under a facing limit answers nothing.
#[test]
fn reach_and_arguments_are_honoured() {
    let mesh = box_mesh([1.0, 1.0, 1.0], 1);
    let index = Bvh::build(&mesh).unwrap();
    let points = [Point::new(0.5, 0.5, 3.0), Point::new(0.5, 0.5, 1.5)];
    let near = DistanceOptions {
        max: 1.0,
        ..DistanceOptions::default()
    };
    let found = to_mesh(&index, &points, None, &near).unwrap();
    assert!(found[0].is_none());
    assert!((found[1].as_ref().unwrap().distance - 0.5).abs() <= 1e-15);
    let facing = DistanceOptions {
        facing_limit: Some(1.0),
        ..DistanceOptions::default()
    };
    assert!(to_mesh(&index, &points, None, &facing).is_err());
    assert!(to_mesh(&index, &points, Some(&[Vector::Z]), &facing).is_err());
    let past = DistanceOptions {
        facing_limit: Some(4.0),
        ..DistanceOptions::default()
    };
    assert!(to_mesh(&index, &points, Some(&[Vector::Z; 2]), &past).is_err());
    let found = to_mesh(&index, &points, Some(&[Vector::ZERO, Vector::Z]), &facing).unwrap();
    assert!(found[0].is_none() && found[1].is_some());
    let unsigned = DistanceOptions {
        signed: false,
        ..DistanceOptions::default()
    };
    let inside = to_mesh(&index, &[Point::new(0.5, 0.5, 0.4)], None, &unsigned).unwrap();
    assert!((inside[0].as_ref().unwrap().distance - 0.4).abs() <= 1e-15);
}

/// Throughput on a million-triangle mesh, a million points about it with
/// a scan's noise, at one thread and at eight, with and without a facing
/// limit. Run in release with `--ignored --nocapture`.
#[test]
#[ignore = "a timing, for release builds"]
fn throughput_on_a_million_triangles() {
    let mut rng = Seeded(0x1_000_000);
    // Scaled to a scan's: a torus 80 mm about its axis and 28 mm about its
    // tube, its triangles some 0.7 by 0.35 mm, its vertices off it by up
    // to 0.05 mm.
    let mut mesh = torus(1000, 500, 5e-4, &mut rng);
    for p in &mut mesh.positions {
        *p = Point::ORIGIN + (*p - Point::ORIGIN) * 40.0;
    }
    assert_eq!(mesh.triangles.len(), 1_000_000);
    let started = ogeom_core::clock::Instant::now();
    let index = Bvh::build(&mesh).unwrap();
    println!(
        "hierarchy built in {:.3} s",
        started.elapsed().as_secs_f64()
    );
    let mut points = Vec::with_capacity(1_000_000);
    let mut normals = Vec::with_capacity(1_000_000);
    for _ in 0..1_000_000 {
        let t = mesh.triangles[rng.below(mesh.triangles.len())];
        let corners = t.map(|k| mesh.positions[k as usize]);
        let n = unit_normal(corners);
        let (u, v) = (rng.range(0.0, 0.5), rng.range(0.0, 0.5));
        let at = corners[0] + (corners[1] - corners[0]) * u + (corners[2] - corners[0]) * v;
        points.push(at + n * rng.range(-0.1, 0.1) + rng.direction() * 0.02);
        normals.push(n + rng.direction() * 0.2);
    }
    for limit in [None, Some(90_f64.to_radians())] {
        let options = DistanceOptions {
            facing_limit: limit,
            ..DistanceOptions::default()
        };
        for threads in [1, 8] {
            ogeom_core::parallel::set_threads(threads);
            let started = ogeom_core::clock::Instant::now();
            let found = to_mesh(&index, &points, Some(&normals), &options).unwrap();
            let seconds = started.elapsed().as_secs_f64();
            assert!(found.iter().all(Option::is_some));
            println!(
                "facing {limit:?}, {threads} threads: {seconds:.3} s, {:.2} M points/s",
                1.0 / seconds
            );
        }
    }
    ogeom_core::parallel::set_threads(0);
}

/// Throughput against a box's and a cylinder's own faces, a hundred
/// thousand points about each, at one thread and at eight. Run in release
/// with `--ignored --nocapture`.
#[test]
#[ignore = "a timing, for release builds"]
fn throughput_on_shapes() {
    let mut rng = Seeded(0x5_4a9e);
    let mut model = Model::new();
    let size = [40.0, 30.0, 20.0];
    let cube = ogeom_algo::make_box(&mut model, tilted(), (size[0], size[1], size[2]), T)
        .unwrap()
        .shape;
    let can = ogeom_algo::make_cylinder(&mut model, tilted(), 15.0, 20.0, T)
        .unwrap()
        .shape;
    let points: Vec<Point> = (0..100_000)
        .map(|i| tilted().to_world(near_box(&mut rng, size, i)))
        .collect();
    for (name, shape) in [("box", &cube), ("cylinder", &can)] {
        for threads in [1, 8] {
            ogeom_core::parallel::set_threads(threads);
            let started = ogeom_core::clock::Instant::now();
            let found =
                to_shape(&model, shape, &points, None, &DistanceOptions::default(), T).unwrap();
            let seconds = started.elapsed().as_secs_f64();
            assert!(found.iter().all(Option::is_some));
            println!(
                "{name}, {threads} threads: {seconds:.3} s, {:.0} points/s",
                1e5 / seconds
            );
        }
    }
    ogeom_core::parallel::set_threads(0);
}
