//! Face directions, axes of revolution, principal frames and open rims read
//! off meshes built in code, against the closed forms the meshes were built
//! from: a box drafted a degree on every wall, a stepped shaft turned about
//! a tilted axis, an open-topped housing with a hole through its floor,
//! and the same with seeded Gaussian vertex noise of 0.02 mm.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    reason = "test code"
)]

use ogeom_algo::orient::{
    Directions, FaceDirection, Rim, axis_of_revolution, face_directions, open_rims, principal_frame,
};
use ogeom_core::Tolerances;
use ogeom_math::{Axis, Direction, Point, Transform, Vector};
use ogeom_topo::Triangulation;

const T: Tolerances = Tolerances::millimetres();
const SIGMA: f64 = 0.02;
const DRAFT_DEGREES: f64 = 1.0;

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
        (-2.0 * u.ln()).sqrt() * (core::f64::consts::TAU * v).cos()
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

/// A planar face from its four corners, wound counter-clockwise about its
/// outward normal: the outward unit normal and the area, in closed form.
#[derive(Debug, Clone, Copy)]
struct Quad {
    corners: [Point; 4],
}

impl Quad {
    fn vector_area(&self) -> Vector {
        let [a, b, c, d] = self.corners;
        (c - a).cross(d - b) * 0.5
    }

    fn area(&self) -> f64 {
        self.vector_area().magnitude()
    }

    fn normal(&self) -> Vector {
        self.vector_area() / self.area()
    }
}

#[derive(Default)]
struct Builder {
    positions: Vec<Point>,
    triangles: Vec<[u32; 3]>,
}

impl Builder {
    /// A planar quad meshed `ns × nt` cells, `s` along its first edge and
    /// `t` along its last, leaving out the cells `skip` names.
    fn quad(&mut self, quad: Quad, ns: usize, nt: usize, skip: impl Fn(usize, usize) -> bool) {
        let [p00, p10, p11, p01] = quad.corners.map(Point::to_vector);
        let base = self.positions.len() as u32;
        for j in 0..=nt {
            for i in 0..=ns {
                let (s, t) = (i as f64 / ns as f64, j as f64 / nt as f64);
                let p = p00 * ((1.0 - s) * (1.0 - t))
                    + p10 * (s * (1.0 - t))
                    + p11 * (s * t)
                    + p01 * ((1.0 - s) * t);
                self.positions.push(Point::from_vector(p));
            }
        }
        let at = |i: usize, j: usize| base + (j * (ns + 1) + i) as u32;
        for j in 0..nt {
            for i in 0..ns {
                if skip(i, j) {
                    continue;
                }
                self.triangles
                    .push([at(i, j), at(i + 1, j), at(i + 1, j + 1)]);
                self.triangles
                    .push([at(i, j), at(i + 1, j + 1), at(i, j + 1)]);
            }
        }
    }

    fn mesh(self, placed: &Transform) -> Triangulation {
        Triangulation {
            positions: self.positions.iter().map(|p| placed.apply(*p)).collect(),
            triangles: self.triangles,
            ..Triangulation::new()
        }
        .welded(T)
    }
}

/// A box `2a × 2b × h` standing on the plane z = 0, every wall drafted
/// inwards by `DRAFT_DEGREES` towards its top, with its faces: bottom, top,
/// +X, −X, +Y, −Y.
fn drafted_faces(a: f64, b: f64, h: f64) -> [Quad; 6] {
    let inset = h * DRAFT_DEGREES.to_radians().tan();
    let bottom = |x: f64, y: f64| Point::new(x * a, y * b, 0.0);
    let top = |x: f64, y: f64| Point::new(x * (a - inset), y * (b - inset), h);
    let q = |corners| Quad { corners };
    [
        q([
            bottom(-1., -1.),
            bottom(-1., 1.),
            bottom(1., 1.),
            bottom(1., -1.),
        ]),
        q([top(-1., -1.), top(1., -1.), top(1., 1.), top(-1., 1.)]),
        q([bottom(1., -1.), bottom(1., 1.), top(1., 1.), top(1., -1.)]),
        q([
            bottom(-1., 1.),
            bottom(-1., -1.),
            top(-1., -1.),
            top(-1., 1.),
        ]),
        q([bottom(1., 1.), bottom(-1., 1.), top(-1., 1.), top(1., 1.)]),
        q([
            bottom(-1., -1.),
            bottom(1., -1.),
            top(1., -1.),
            top(-1., -1.),
        ]),
    ]
}

/// The drafted box meshed `nx × ny × nz` cells, its faces' edges meshed
/// alike so they weld; `open_top` leaves the top out, and `hole` the floor's
/// cells in that range of each axis.
fn drafted_box(
    (a, b, h): (f64, f64, f64),
    (nx, ny, nz): (usize, usize, usize),
    open_top: bool,
    hole: Option<((usize, usize), (usize, usize))>,
) -> Builder {
    let faces = drafted_faces(a, b, h);
    let mut built = Builder::default();
    built.quad(faces[0], ny, nx, |i, j| {
        hole.is_some_and(|((x0, x1), (y0, y1))| (y0..y1).contains(&i) && (x0..x1).contains(&j))
    });
    if !open_top {
        built.quad(faces[1], nx, ny, |_, _| false);
    }
    built.quad(faces[2], ny, nz, |_, _| false);
    built.quad(faces[3], ny, nz, |_, _| false);
    built.quad(faces[4], nx, nz, |_, _| false);
    built.quad(faces[5], nx, nz, |_, _| false);
    built
}

/// Every vertex moved by Gaussian noise of `SIGMA` along each axis.
fn noisy(mut mesh: Triangulation, seed: u64) -> Triangulation {
    let mut draw = Seeded(seed);
    for p in &mut mesh.positions {
        *p += Vector::new(draw.gaussian(), draw.gaussian(), draw.gaussian()) * SIGMA;
    }
    mesh
}

/// The axial angle between two directions, in degrees.
fn degrees_apart(a: Vector, b: Vector) -> f64 {
    a.cross(b).magnitude().atan2(a.dot(b).abs()).to_degrees()
}

/// The area-weighted mean of the normals of `faces`, each turned to agree
/// with the first: the direction the drafted planes average to.
fn average(faces: &[Quad], placed: &Transform) -> Vector {
    let first = faces[0].normal();
    let sum = faces.iter().fold(Vector::ZERO, |s, f| {
        let n = f.normal();
        let n = if n.dot(first) < 0.0 { -n } else { n };
        s + n * f.area()
    });
    placed.apply_vector(sum / sum.magnitude())
}

/// The face direction nearest `want`, and how far from it in degrees.
fn nearest(found: &[FaceDirection], want: Vector) -> (FaceDirection, f64) {
    found
        .iter()
        .map(|f| (*f, degrees_apart(f.direction.vector(), want)))
        .min_by(|p, q| p.1.total_cmp(&q.1))
        .unwrap()
}

/// The box's three face directions, each the average of the two drafted
/// planes it merges, with the area of the two.
fn box_expectations(faces: &[Quad; 6], placed: &Transform) -> [(Vector, f64); 3] {
    [(0, 1), (2, 3), (4, 5)].map(|(i, j)| {
        (
            average(&[faces[i], faces[j]], placed),
            faces[i].area() + faces[j].area(),
        )
    })
}

#[test]
fn a_drafted_box_faces_the_average_of_each_pair_of_drafted_planes() {
    let size = (20.0, 15.0, 12.0);
    let placed = placement();
    let mesh = drafted_box(size, (10, 8, 6), false, None).mesh(&placed);
    assert!(mesh.is_closed());
    let found = face_directions(&mesh, &Directions::default(), T).unwrap();
    assert_eq!(found.len(), 3, "{found:?}");
    let faces = drafted_faces(size.0, size.1, size.2);
    for (want, area) in box_expectations(&faces, &placed) {
        let (got, apart) = nearest(&found, want);
        assert!(apart <= 0.05, "{apart}° from the drafted planes' average");
        assert!(
            apart <= 1e-9,
            "{apart}°: the faces are exact, so is the mean"
        );
        assert!(
            (got.area - area).abs() <= 1e-9 * area,
            "{} for {area}",
            got.area
        );
    }
    // The walls are a degree either side of their mean, the top and bottom
    // on it; the bottom is the larger, so the vertical direction faces down.
    let up = placed.apply_vector(Vector::Z);
    let (vertical, _) = nearest(&found, up);
    assert!(vertical.spread <= 1e-9);
    assert!(vertical.direction.vector().dot(up) < 0.0);
    for wall in [Vector::X, Vector::Y] {
        let (w, _) = nearest(&found, placed.apply_vector(wall));
        assert!(
            (w.spread.to_degrees() - DRAFT_DEGREES).abs() <= 1e-9,
            "{}",
            w.spread
        );
    }
}

#[test]
fn a_noisy_drafted_box_still_faces_the_planes_average() {
    // Cells of half a millimetre, so each triangle's normal is a few
    // degrees off its face; the face's mean normal is its vector area over
    // its area, which the noise moves only through the face's border.
    let size = (20.0, 15.0, 12.0);
    let placed = placement();
    let mesh = noisy(
        drafted_box(size, (80, 60, 24), false, None).mesh(&placed),
        7,
    );
    let found = face_directions(&mesh, &Directions::default(), T).unwrap();
    let faces = drafted_faces(size.0, size.1, size.2);
    let total: f64 = faces.iter().map(Quad::area).sum();
    for (want, area) in box_expectations(&faces, &placed) {
        let (got, apart) = nearest(&found, want);
        assert!(apart <= 0.05, "{apart}° from the drafted planes' average");
        assert!(
            (got.area - area).abs() <= 0.02 * area,
            "{} for {area}",
            got.area
        );
    }
    let rest: f64 = found.iter().skip(3).map(|f| f.area).sum();
    assert!(
        rest <= 0.01 * total,
        "{rest} of {total} off the faces: {found:?}"
    );
}

/// A stepped shaft about the z axis, `segments` round: its profile as
/// (radius, height) from the bottom centre to the top centre.
const PROFILE: [(f64, f64); 8] = [
    (0.0, 0.0),
    (6.0, 0.0),
    (6.0, 20.0),
    (10.0, 20.0),
    (10.0, 35.0),
    (8.0, 35.0),
    (8.0, 60.0),
    (0.0, 60.0),
];

/// The shaft meshed `segments` round, each profile segment cut into pieces
/// a millimetre or two long; with its wall area and its whole area.
fn shaft(segments: usize, placed: &Transform) -> (Triangulation, f64, f64) {
    let mut profile: Vec<(f64, f64)> = vec![PROFILE[0]];
    for pair in PROFILE.windows(2) {
        let ((r0, z0), (r1, z1)) = (pair[0], pair[1]);
        let pieces = ((r1 - r0).hypot(z1 - z0) / 2.0).ceil() as usize;
        for k in 1..=pieces {
            let s = k as f64 / pieces as f64;
            profile.push((r0 + (r1 - r0) * s, z0 + (z1 - z0) * s));
        }
    }
    let mut positions = Vec::new();
    let mut rings: Vec<Vec<u32>> = Vec::new();
    for &(r, z) in &profile {
        let mut ring = Vec::new();
        if r == 0.0 {
            ring.push(positions.len() as u32);
            positions.push(Point::new(0.0, 0.0, z));
        } else {
            for k in 0..segments {
                let angle = core::f64::consts::TAU * k as f64 / segments as f64;
                ring.push(positions.len() as u32);
                positions.push(Point::new(r * angle.cos(), r * angle.sin(), z));
            }
        }
        rings.push(ring);
    }
    let mut triangles = Vec::new();
    for pair in rings.windows(2) {
        let (lower, upper) = (&pair[0], &pair[1]);
        for k in 0..segments {
            let at = |ring: &Vec<u32>, k: usize| ring[k % ring.len()];
            let (a, b) = (at(lower, k), at(lower, k + 1));
            let (c, d) = (at(upper, k + 1), at(upper, k));
            if lower.len() > 1 {
                triangles.push([a, b, c]);
            }
            if upper.len() > 1 {
                triangles.push([a, c, d]);
            }
        }
    }
    let mut mesh = Triangulation {
        positions,
        triangles,
        ..Triangulation::new()
    };
    if mesh.volume() < 0.0 {
        for t in &mut mesh.triangles {
            t.swap(1, 2);
        }
    }
    let chord = 2.0 * (core::f64::consts::PI / segments as f64).sin();
    let mut walls = 0.0;
    let mut flats = 0.0;
    for pair in PROFILE.windows(2) {
        let ((r0, z0), (r1, z1)) = (pair[0], pair[1]);
        if r0 == r1 {
            walls += segments as f64 * chord * r0 * (z1 - z0).abs();
        } else {
            let polygon = |r: f64| {
                segments as f64
                    * 0.5
                    * r
                    * r
                    * chord
                    * (core::f64::consts::PI / segments as f64).cos()
            };
            flats += (polygon(r1) - polygon(r0)).abs();
        }
    }
    mesh.positions = mesh.positions.iter().map(|p| placed.apply(*p)).collect();
    (mesh, walls, walls + flats)
}

#[test]
fn a_turned_shaft_turns_about_its_axis() {
    let placed = placement();
    let (mesh, walls, total) = shaft(96, &placed);
    assert!(mesh.is_closed());
    let found = axis_of_revolution(&mesh, &Directions::default(), T)
        .unwrap()
        .unwrap();
    let axis = placed.apply_vector(Vector::Z);
    let through = placed.apply(Point::ORIGIN);
    let apart = degrees_apart(found.axis.direction.vector(), axis);
    let off = (found.axis.location - through).cross(axis).magnitude();
    assert!(apart <= 1e-9, "{apart}° off the axis");
    assert!(off <= 1e-9, "{off} from the axis");
    assert!(
        (found.share - walls / total).abs() <= 1e-12,
        "{} for {}",
        found.share,
        walls / total
    );

    // Its ends and shoulders face along the axis; the walls sweep a circle
    // of normals and face no direction.
    let directions = face_directions(&mesh, &Directions::default(), T).unwrap();
    assert_eq!(directions.len(), 1, "{directions:?}");
    assert!(degrees_apart(directions[0].direction.vector(), axis) <= 1e-9);
    assert!((directions[0].area - (total - walls)).abs() <= 1e-9 * total);
}

#[test]
fn a_noisy_turned_shaft_turns_about_its_axis() {
    let placed = placement();
    let (mesh, walls, total) = shaft(180, &placed);
    let mesh = noisy(mesh, 11);
    let found = axis_of_revolution(&mesh, &Directions::default(), T)
        .unwrap()
        .unwrap();
    let axis = placed.apply_vector(Vector::Z);
    let through = placed.apply(Point::ORIGIN);
    let apart = degrees_apart(found.axis.direction.vector(), axis);
    let off = (found.axis.location - through).cross(axis).magnitude();
    assert!(apart <= 0.01, "{apart}° off the axis");
    assert!(off <= 2.0 * SIGMA, "{off} from the axis");
    assert!(
        (found.share - walls / total).abs() <= 0.02,
        "{} for {}",
        found.share,
        walls / total
    );
}

#[test]
fn a_box_is_no_turned_part() {
    let placed = placement();
    let mesh = drafted_box((20.0, 20.0, 20.0), (10, 10, 10), false, None).mesh(&placed);
    let found = axis_of_revolution(&mesh, &Directions::default(), T)
        .unwrap()
        .unwrap();
    // Its walls point at an axis only along a strip the resolution wide.
    assert!(
        found.share <= Directions::default().resolution.sin(),
        "{}",
        found.share
    );
}

#[test]
fn an_open_topped_housing_reports_its_opening_and_the_hole_through_its_floor() {
    let size = (20.0, 15.0, 12.0);
    let placed = placement();
    // A hole 6 × 3 mm through the floor: three cells of 2 mm by two of 1.5.
    let mesh = drafted_box(size, (20, 20, 6), true, Some(((10, 13), (8, 10)))).mesh(&placed);
    let rims = open_rims(&mesh).unwrap();
    assert_eq!(rims.len(), 2, "{rims:?}");
    let faces = drafted_faces(size.0, size.1, size.2);
    let opening = placed.apply_vector(faces[1].vector_area());
    let top_centre = placed.apply(Point::new(0.0, 0.0, size.2));
    let hole = placed.apply_vector(Vector::new(0.0, 0.0, -6.0 * 3.0));
    let hole_centre = placed.apply(Point::new(3.0, -1.5, 0.0));
    let check = |rim: &Rim, area: Vector, centre: Point| {
        assert!(rim.closed);
        assert!(
            (rim.vector_area - area).magnitude() <= 1e-9 * area.magnitude(),
            "{:?} for {area:?}",
            rim.vector_area
        );
        assert!(
            rim.centre.distance(centre) <= 1e-9,
            "{:?} for {centre:?}",
            rim.centre
        );
        assert!(rim.flatness <= 1e-9);
    };
    check(&rims[0], opening, top_centre);
    check(&rims[1], hole, hole_centre);
    assert_eq!(rims[0].vertices.len(), 2 * (20 + 20));
    assert_eq!(rims[1].vertices.len(), 2 * (3 + 2));

    // The floor faces down, and with no top to balance it so does the
    // vertical direction.
    let found = face_directions(&mesh, &Directions::default(), T).unwrap();
    let (vertical, _) = nearest(&found, opening);
    assert!(vertical.direction.vector().dot(opening) < 0.0);
}

#[test]
fn a_closed_mesh_has_no_rims() {
    let mesh = drafted_box((5.0, 4.0, 3.0), (3, 3, 3), false, None).mesh(&placement());
    assert!(open_rims(&mesh).unwrap().is_empty());
}

/// A box `2a × 2b × 2c` about the origin, each face meshed `n × n`.
fn cuboid((a, b, c): (f64, f64, f64), n: usize) -> Builder {
    let faces = [
        [
            (-1., -1., -1.),
            (-1., 1., -1.),
            (1., 1., -1.),
            (1., -1., -1.),
        ],
        [(-1., -1., 1.), (1., -1., 1.), (1., 1., 1.), (-1., 1., 1.)],
        [(1., -1., -1.), (1., 1., -1.), (1., 1., 1.), (1., -1., 1.)],
        [
            (-1., 1., -1.),
            (-1., -1., -1.),
            (-1., -1., 1.),
            (-1., 1., 1.),
        ],
        [(1., 1., -1.), (-1., 1., -1.), (-1., 1., 1.), (1., 1., 1.)],
        [
            (-1., -1., -1.),
            (1., -1., -1.),
            (1., -1., 1.),
            (-1., -1., 1.),
        ],
    ];
    let mut built = Builder::default();
    for f in faces {
        let corners = f.map(|(x, y, z)| Point::new(a * x, b * y, c * z));
        built.quad(Quad { corners }, n, n, |_, _| false);
    }
    built
}

#[test]
fn a_box_spreads_along_its_edges_longest_first_towards_where_it_reaches_further() {
    let placed = placement();
    let mut built = cuboid((15.0, 10.0, 5.0), 4);
    // A tab of 2 × 2 off the +x end, so the area reaches further that way.
    built.quad(
        Quad {
            corners: [
                Point::new(20.0, -1.0, 0.0),
                Point::new(22.0, -1.0, 0.0),
                Point::new(22.0, 1.0, 0.0),
                Point::new(20.0, 1.0, 0.0),
            ],
        },
        1,
        1,
        |_, _| false,
    );
    let mesh = built.mesh(&placed);
    let frame = principal_frame(&mesh, T).unwrap();
    let box_area = 2.0 * (30.0 * 20.0 + 20.0 * 10.0 + 30.0 * 10.0);
    let centroid = placed.apply(Point::new(4.0 * 21.0 / (box_area + 4.0), 0.0, 0.0));
    assert!(
        frame.origin().distance(centroid) <= 1e-9,
        "{:?}",
        frame.origin()
    );
    let x = placed.apply_vector(Vector::X);
    assert!(frame.x().vector().dot(x) > 0.0);
    assert!(degrees_apart(frame.x().vector(), x) <= 1e-9);
    assert!(degrees_apart(frame.y().vector(), placed.apply_vector(Vector::Y)) <= 1e-9);
    assert!(degrees_apart(frame.z().vector(), placed.apply_vector(Vector::Z)) <= 1e-9);
}

#[test]
fn a_cubes_equal_spreads_take_the_world_axes() {
    let shift = Transform::translation(Vector::new(1.0, 2.0, 3.0));
    let frame = principal_frame(&cuboid((5.0, 5.0, 5.0), 3).mesh(&shift), T).unwrap();
    assert!(frame.origin().distance(Point::new(1.0, 2.0, 3.0)) <= 1e-12);
    assert_eq!(frame.x().vector(), Vector::X);
    assert_eq!(frame.y().vector(), Vector::Y);
    assert_eq!(frame.z().vector(), Vector::Z);
}

#[test]
fn the_answers_do_not_depend_on_threads_or_on_the_order_of_the_triangles() {
    let placed = placement();
    let mesh = noisy(
        drafted_box((20.0, 15.0, 12.0), (40, 30, 12), false, None).mesh(&placed),
        3,
    );
    let options = Directions::default();
    ogeom_core::parallel::set_threads(1);
    let one = face_directions(&mesh, &options, T).unwrap();
    let axis_one = axis_of_revolution(&mesh, &options, T).unwrap();
    ogeom_core::parallel::set_threads(8);
    let many = face_directions(&mesh, &options, T).unwrap();
    let axis_many = axis_of_revolution(&mesh, &options, T).unwrap();
    ogeom_core::parallel::set_threads(0);
    assert_eq!(one, many);
    assert_eq!(axis_one, axis_many);

    let mut shuffled = mesh.clone();
    shuffled.triangles.reverse();
    for t in &mut shuffled.triangles {
        t.rotate_left(1);
    }
    let again = face_directions(&shuffled, &options, T).unwrap();
    assert_eq!(one.len(), again.len());
    for (p, q) in one.iter().zip(&again) {
        assert!(
            p.direction.vector().dot(q.direction.vector()) > 1.0 - 1e-12,
            "{p:?} {q:?}"
        );
        assert!((p.area - q.area).abs() <= 1e-9 * p.area);
    }
}
