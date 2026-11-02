//! Wall thickness on a B-rep at points on its edges and vertices: the
//! faces a point lies on are not walls for its own rays or balls. A tube's
//! rims read as single-precision points read the wall, rotated boxes read
//! their thickness at an edge and at a corner, and a cylinder reads its
//! diameter and its height on its seam; a ball past the thickest wall
//! asked for reads none; points inside a corpus part's faces read what
//! they read before.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    reason = "points rounded to single precision, as a mesh stores them"
)]

use ogeom::algo::thickness::{Method, Thickness, thickness_of_shape};
use ogeom::algo::{make_box, make_cylinder};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::topo::Model;

const T: Tolerances = Tolerances::millimetres();

const RAY: Method = Method::Ray { cone: 0.0, rays: 1 };

fn options(method: Method, facing_limit: f64, max: f64) -> Thickness {
    Thickness {
        method,
        facing_limit,
        max,
        tol: T,
    }
}

/// A point or a vector rounded to single precision, as a mesh stores it.
fn single(x: f64, y: f64, z: f64) -> [f64; 3] {
    [
        f64::from(x as f32),
        f64::from(y as f32),
        f64::from(z as f32),
    ]
}

/// A frame turned off every axis, so no face's normal is exactly square to
/// a direction along another.
fn turned() -> Frame {
    Frame::new(
        Point::new(1.0, 2.0, 3.0),
        Direction::new(Vector::new(0.2, -0.3, 0.9), T).unwrap(),
        Direction::new(Vector::new(1.0, 0.4, 0.1), T).unwrap(),
        T,
    )
    .unwrap()
}

/// A tube of radii 10 and 8 and height 30, cut from one cylinder by
/// another, with its rims sampled as a mesh's single-precision vertices
/// are: up to 6.7e-7 off the walls, past their tolerance, on the end faces.
/// Each rim point with its wall's normal reads the wall, 2, by rays and by
/// balls (the point is off the wall by less than a millionth, which the
/// reading takes in); with the axial normal it reads none within 7, the
/// tube being 30 long.
#[test]
fn a_tube_read_at_single_precision_rim_points_reads_its_wall() {
    let mut model = Model::new();
    let outer = make_cylinder(&mut model, Frame::WORLD, 10.0, 30.0, T)
        .unwrap()
        .shape;
    let inner = make_cylinder(&mut model, Frame::WORLD, 8.0, 30.0, T)
        .unwrap()
        .shape;
    let tube = ogeom::boolean::cut(&mut model, &outer, &inner, T)
        .unwrap()
        .shape;

    let (mut points, mut walls, mut axial) = (Vec::new(), Vec::new(), Vec::new());
    for z in [0.0, 30.0] {
        for (radius, side) in [(10.0, 1.0), (8.0, -1.0)] {
            for k in 0..64 {
                let a = std::f64::consts::TAU * f64::from(k) / 64.0;
                let [x, y, z] = single(radius * a.cos(), radius * a.sin(), z);
                points.push(Point::new(x, y, z));
                let [nx, ny, _] = single(side * a.cos(), side * a.sin(), 0.0);
                walls.push(Vector::new(nx, ny, 0.0));
                axial.push(Vector::new(0.0, 0.0, if z == 0.0 { -1.0 } else { 1.0 }));
            }
        }
    }
    // The issue's own point, a mesh vertex 6.2e-7 outside the outer wall.
    points.push(Point::new(9.569_404, 2.902_846_8, 0.0));
    walls.push(Vector::new(0.956_940_35, 0.290_284_66, 0.0));
    axial.push(-Vector::Z);

    let limit = 60f64.to_radians();
    for method in [RAY, Method::Sphere] {
        let read = thickness_of_shape(&model, &tube, &points, &walls, &options(method, limit, 7.0))
            .unwrap();
        for (p, w) in points.iter().zip(read) {
            let w = w.unwrap_or_else(|| panic!("{method:?} read nothing at {p:?}"));
            assert!((w - 2.0).abs() <= 1e-6, "{method:?} read {w} at {p:?}");
        }
        let read = thickness_of_shape(&model, &tube, &points, &axial, &options(method, limit, 7.0))
            .unwrap();
        for (p, w) in points.iter().zip(read) {
            assert!(w.is_none(), "{method:?} read {w:?} along the axis at {p:?}");
        }
    }
}

/// A 20 x 15 x 3 box turned off the axes, read at a point of the edge
/// between its top and its side x = 0, at its corner on those and y = 0,
/// and at a point of its bottom edge on y = 0: along the normal of the top
/// or the bottom both methods read the 3 between them, the side faces
/// through the point running edge-on along the ray and the ball's radius.
/// Along the side faces' normals, past the facing limit of the top and
/// bottom, the corner reads the box's length and width.
#[test]
fn a_turned_box_reads_its_thickness_at_an_edge_and_a_corner() {
    let mut model = Model::new();
    let frame = turned();
    let (a, b, c) = (20.0, 15.0, 3.0);
    let block = make_box(&mut model, frame, (a, b, c), T).unwrap().shape;
    let at = |x: f64, y: f64, z: f64| frame.to_world(Point::new(x, y, z));
    let along = |v: Vector| frame.vector_to_world(v);
    let square = [
        (at(0.0, 6.0, c), along(Vector::Z), c),
        (at(0.0, 0.0, c), along(Vector::Z), c),
        (at(7.0, 0.0, 0.0), along(-Vector::Z), c),
    ];
    let limit = 60f64.to_radians();
    let lengthwise = [
        (at(0.0, 0.0, c), along(-Vector::X), a),
        (at(0.0, 0.0, c), along(-Vector::Y), b),
    ];
    for method in [RAY, Method::Sphere] {
        for (cases, facing) in [(&square[..], 0.0), (&lengthwise[..], limit)] {
            for &(p, n, wall) in cases {
                let w =
                    thickness_of_shape(&model, &block, &[p], &[n], &options(method, facing, 100.0))
                        .unwrap()[0]
                        .unwrap();
                assert!(
                    (w - wall).abs() <= 1e-9,
                    "{method:?} read {w} for {wall} at {p:?}"
                );
            }
        }
    }
}

/// A cylinder of radius 4 and height 30 turned off the axes, read on its
/// seam: across, at mid-height and on its bottom rim, it reads its
/// diameter, the far side of the face the point lies on met as a wall;
/// along the axis from the rim it reads its height, the side face and the
/// seam running along the ray.
#[test]
fn a_turned_cylinder_reads_its_diameter_and_height_on_its_seam() {
    let mut model = Model::new();
    let frame = turned();
    let (radius, height) = (4.0, 30.0);
    let rod = make_cylinder(&mut model, frame, radius, height, T)
        .unwrap()
        .shape;
    let at = |z: f64| frame.to_world(Point::new(radius, 0.0, z));
    let across = frame.vector_to_world(Vector::X);
    let down = frame.vector_to_world(-Vector::Z);
    let cases = [
        (at(10.0), across, 2.0 * radius),
        (at(0.0), across, 2.0 * radius),
        (at(0.0), down, height),
    ];
    for method in [RAY, Method::Sphere] {
        for &(p, n, wall) in &cases {
            let w = thickness_of_shape(&model, &rod, &[p], &[n], &options(method, 0.0, 100.0))
                .unwrap()[0]
                .unwrap();
            assert!(
                (w - wall).abs() <= 1e-9,
                "{method:?} read {w} for {wall} at {p:?}"
            );
        }
    }
}

/// A ball whose wall is thicker than the thickest asked for reads none, as
/// a ray does: the face the point lies on, which every ball tangent there
/// touches at the point, does not bound it. A 40 x 30 x 3 plate read at
/// its top with walls to 2, and a rod of radius 4 read on its side with
/// walls to 7, read none; with walls to 9 the rod reads its diameter.
#[test]
fn a_ball_past_the_thickest_wall_reads_none() {
    let mut model = Model::new();
    let plate = make_box(&mut model, Frame::WORLD, (40.0, 30.0, 3.0), T)
        .unwrap()
        .shape;
    let rod = make_cylinder(&mut model, Frame::WORLD, 4.0, 30.0, T)
        .unwrap()
        .shape;
    let top = Point::new(10.0, 10.0, 3.0);
    let side = Point::new(0.0, 4.0, 10.0);
    let limit = 30f64.to_radians();
    for method in [RAY, Method::Sphere] {
        let read = |shape, p, n, max| {
            thickness_of_shape(&model, shape, &[p], &[n], &options(method, limit, max)).unwrap()[0]
        };
        assert_eq!(
            read(&plate, top, Vector::Z, 2.0),
            None,
            "{method:?} on the plate"
        );
        assert_eq!(
            read(&rod, side, Vector::Y, 7.0),
            None,
            "{method:?} on the rod"
        );
        let w = read(&rod, side, Vector::Y, 9.0).unwrap();
        assert!((w - 8.0).abs() <= 1e-9, "{method:?} read {w} on the rod");
    }
}

/// A point, its outward normal, and the bits of what it read by rays and
/// by a ball before points on edges took on the faces at their nearest
/// point.
type Reading = ([f64; 3], [f64; 3], Option<u64>, Option<u64>);

/// Points inside the faces of the first body of `m5x16_bhcs.step`, up to
/// three to a face, at triangle centres of each face's mesh carried to the
/// face and held inside its trim, clear of the band; read by five rays in
/// a cone of 0.2 and by a ball, facing limit 30 degrees, walls to 50.
const INSIDE: &[Reading] = &[
    (
        [7.935, -0.6468212050930554, -2.414875220097711],
        [-0.0, -0.25872848203722215, -0.9659500880390844],
        Some(0x4013cdbd2ee25fc7),
        Some(0x4014000000000000),
    ),
    (
        [0.0, 1.5242721801906214, 0.9398264344738518],
        [-1.0, 0.0, 0.0],
        Some(0x4032c00000000000),
        Some(0x3ff9d288821ad499),
    ),
    (
        [0.0, -1.0314592863576935, -1.8896514063109506],
        [-1.0, 0.0, 0.0],
        Some(0x4006728c6248f03d),
        Some(0x3ffacc4412df2666),
    ),
    (
        [0.0, -1.5067147794403448, 1.4996193029132014],
        [-1.0, 0.0, 0.0],
        Some(0x4006557a06300053),
        Some(0x3ffc1d8a7e32427d),
    ),
    (
        [1.5622991014825192, 0.009453062392351072, -4.249384785126207],
        [
            -0.6006825390847741,
            0.001778513415261079,
            -0.7994856622408563,
        ],
        Some(0x3ffc20effb2f4409),
        Some(0x3ff7be6f8530dfe8),
    ),
    (
        [1.9730482210835105, -0.38688404524936243, -4.512402627597604],
        [
            -0.5234035802461962,
            -0.07278894775762723,
            -0.848970235797346,
        ],
        Some(0x3ff4a5f4ea7176d4),
        Some(0x3ff0520184e048a3),
    ),
    (
        [2.4333333333333336, -2.1009210851280185, 4.260120960027369],
        [0.0, -0.4422991758164249, 0.8968675705320778],
        Some(0x400ced7bf97f50cd),
        Some(0x3fe0e0ec4138c409),
    ),
    (
        [2.75, 2.859166289229755, 1.5282562211284068],
        [1.0, -0.0, -0.0],
        Some(0x3fff515d99e4506c),
        Some(0x4001a7c58ed89bfa),
    ),
    (
        [2.75, -2.932774605691824, -2.713297421126595],
        [1.0, -0.0, -0.0],
        Some(0x3ff4704a5c1cafb4),
        Some(0x3ff2121df22bec6e),
    ),
    (
        [2.75, 2.932774605691824, 2.7132974211265943],
        [1.0, -0.0, -0.0],
        Some(0x3ff41c4f930f9d88),
        Some(0x3ff801782569e396),
    ),
    (
        [18.458513836438264, -0.6164760368860648, -2.2640571131780023],
        [
            0.7071067811865424,
            -0.18577326084014625,
            -0.6822670265789095,
        ],
        Some(0x40183b35515812e6),
        Some(0x3ff0c5aa3a444be6),
    ),
    (
        [18.75, 1.8874141641974749, 0.6427037607885143],
        [1.0, -0.0, -0.0],
        Some(0x4030148398fc1334),
        Some(0x3ff054e30dd05663),
    ),
    (
        [18.75, 0.5695566844272428, 0.38056560961842756],
        [1.0, -0.0, -0.0],
        Some(0x40303e2a330ba325),
        Some(0x401f64194ff6cbab),
    ),
    (
        [18.75, -1.8742520812736876, -0.7763406313695627],
        [1.0, -0.0, -0.0],
        Some(0x4030293c10326e4e),
        Some(0x3fc044f4493faedc),
    ),
    (
        [1.04, -1.0103629710818467, 1.250000000000001],
        [0.0, 0.8660254037844385, -0.5000000000000003],
        Some(0x3ffe5e50e571f5e8),
        Some(0x3ffbdbb53a0c193f),
    ),
    (
        [1.04, -1.5877132402714735, 0.25000000000000006],
        [0.0, 0.8660254037844385, -0.5000000000000003],
        Some(0x3ffe12a7260931c2),
        Some(0x3ffbdbb53a0c1940),
    ),
    (
        [0.52, -1.2990381056766602, 0.7500000000000004],
        [0.0, 0.8660254037844385, -0.5000000000000003],
        Some(0x3ff5ae30543ac5c8),
        Some(0x3ffb6364d88e8da8),
    ),
    (
        [1.04, -1.5877132402714715, -0.24999999999999867],
        [-0.0, 0.8660254037844385, 0.5000000000000003],
        Some(0x3ffe5e50e571f5ee),
        Some(0x3ffbdbb53a0c1949),
    ),
    (
        [1.04, -1.0103629710818447, -1.2499999999999996],
        [-0.0, 0.8660254037844385, 0.5000000000000003],
        Some(0x3ffe12a7260931cc),
        Some(0x3ffbdbb53a0c1949),
    ),
    (
        [0.52, -1.2990381056766582, -0.7499999999999992],
        [-0.0, 0.8660254037844385, 0.5000000000000003],
        Some(0x3ff5ae30543ac5ce),
        Some(0x3ffb6364d88e8db2),
    ),
    (
        [1.04, -0.5773502691896253, -1.5],
        [0.0, 0.0, 1.0],
        Some(0x3ffe5e50e571f5f0),
        Some(0x3ffbdbb53a0c1948),
    ),
    (
        [1.04, 0.5773502691896261, -1.5],
        [0.0, 0.0, 1.0],
        Some(0x3ffe12a7260931c8),
        Some(0x40001b15ba94252f),
    ),
    (
        [0.52, 3.3306690738754696e-16, -1.5],
        [0.0, 0.0, 1.0],
        Some(0x3ff5ae30543ac5cc),
        Some(0x3ffb6364d88e8db0),
    ),
    (
        [1.04, 1.0103629710818474, -1.2500000000000013],
        [0.0, -0.8660254037844388, 0.4999999999999999],
        Some(0x3ffe5e50e571f5e4),
        Some(0x40020ea11016087e),
    ),
    (
        [1.04, 1.5877132402714735, -0.2500000000000002],
        [0.0, -0.8660254037844388, 0.4999999999999999],
        Some(0x3ffe12a7260931c0),
        Some(0x40020ea11016087e),
    ),
    (
        [0.52, 1.2990381056766607, -0.7500000000000007],
        [0.0, -0.8660254037844388, 0.4999999999999999],
        Some(0x3ff5ae30543ac5c2),
        Some(0x3ffb6364d88e8da1),
    ),
    (
        [1.04, 1.5877132402714718, 0.24999999999999867],
        [0.0, -0.8660254037844388, -0.4999999999999999],
        Some(0x3ffe5e50e571f5ea),
        Some(0x40020ea110160882),
    ),
    (
        [1.04, 1.0103629710818456, 1.2499999999999996],
        [0.0, -0.8660254037844388, -0.4999999999999999],
        Some(0x3ffe12a7260931c6),
        Some(0x40020ea110160884),
    ),
    (
        [0.52, 1.2990381056766587, 0.7499999999999991],
        [0.0, -0.8660254037844388, -0.4999999999999999],
        Some(0x3ff5ae30543ac5cb),
        Some(0x3ffb6364d88e8dac),
    ),
    (
        [1.04, 0.5773502691896261, 1.5],
        [0.0, 0.0, -1.0],
        Some(0x3ffe5e50e571f5ec),
        Some(0x40001b15ba94252f),
    ),
    (
        [1.04, -0.5773502691896254, 1.5],
        [0.0, 0.0, -1.0],
        Some(0x3ffe12a7260931cc),
        Some(0x3ffbdbb53a0c1948),
    ),
    (
        [0.52, 3.3306690738754696e-16, 1.5],
        [0.0, 0.0, -1.0],
        Some(0x3ff5ae30543ac5cc),
        Some(0x3ffb6364d88e8db0),
    ),
    (
        [1.56, 0.6680846571460737, -1.4159756150367533],
        [-1.0, -0.0, -0.0],
        Some(0x40310bb0885e503e),
        Some(0x3ffec705e9cc01b0),
    ),
    (
        [1.56, 0.8922285251880873, -1.286566092485493],
        [-1.0, -0.0, -0.0],
        Some(0x403130a3d70a3d71),
        Some(0x3ffec705e9cc01ae),
    ),
    (
        [1.56, 1.075241227080307, -1.1035533905932735],
        [-1.0, -0.0, -0.0],
        Some(0x4031191bccd747ce),
        Some(0x3fff6950099686d8),
    ),
    (
        [1.56, 1.5603131823341612, -0.12940952255126018],
        [-1.0, -0.0, -0.0],
        Some(0x403130a3d70a3d71),
        Some(0x3ffec705e9cc01a9),
    ),
    (
        [1.56, 1.5603131823341612, 0.12940952255126026],
        [-1.0, -0.0, -0.0],
        Some(0x403130a3d70a3d71),
        Some(0x3ffec705e9cc01ab),
    ),
    (
        [1.56, 1.4933258842263808, 0.37940952255126026],
        [-1.0, -0.0, -0.0],
        Some(0x403130a3d70a3d71),
        Some(0x3fff6950099686c9),
    ),
    (
        [1.56, 0.8922285251880873, 1.286566092485493],
        [-1.0, -0.0, -0.0],
        Some(0x403130a3d70a3d71),
        Some(0x3ffec705e9cc01b0),
    ),
    (
        [1.56, 0.6680846571460737, 1.4159756150367533],
        [-1.0, -0.0, -0.0],
        Some(0x403130a3d70a3d71),
        Some(0x3ffec705e9cc01ae),
    ),
    (
        [1.56, 0.4180846571460737, 1.482962913144534],
        [-1.0, -0.0, -0.0],
        Some(0x4031108aad20f0c1),
        Some(0x3fff6950099686d1),
    ),
    (
        [1.56, -0.6680846571460732, 1.4159756150367533],
        [-1.0, -0.0, -0.0],
        Some(0x403130a3d70a3d71),
        Some(0x3fdf5ef6c18e8b0f),
    ),
    (
        [1.56, -0.892228525188087, 1.2865660924854925],
        [-1.0, -0.0, -0.0],
        Some(0x403130a3d70a3d71),
        Some(0x3fdf5ef6c18e8aee),
    ),
    (
        [1.56, -1.075241227080307, 1.1035533905932733],
        [-1.0, -0.0, -0.0],
        Some(0x4030fa14dd47d471),
        Some(0x3fd379e6c414cd2e),
    ),
    (
        [1.56, -1.560313182334161, -0.1294095225512604],
        [-1.0, -0.0, -0.0],
        Some(0x403130a3d70a3d71),
        Some(0x3fdf5ef6c18e8bb0),
    ),
    (
        [1.56, -1.4933258842263808, -0.3794095225512603],
        [-1.0, -0.0, -0.0],
        Some(0x403129894d2fea85),
        Some(0x3fd379e6c414cdcc),
    ),
    (
        [1.56, -1.560313182334161, 0.12940952255126006],
        [-1.0, -0.0, -0.0],
        Some(0x403130a3d70a3d71),
        Some(0x3fdf5ef6c18e8ba3),
    ),
    (
        [1.56, -0.6680846571460732, -1.4159756150367535],
        [-1.0, -0.0, -0.0],
        Some(0x403130a3d70a3d71),
        Some(0x3fdf5ef6c18e8b14),
    ),
    (
        [1.56, -0.418084657146073, -1.482962913144534],
        [-1.0, -0.0, -0.0],
        Some(0x403130a3d70a3d71),
        Some(0x3fd379e6c414cd00),
    ),
    (
        [1.56, -0.8922285251880866, -1.286566092485493],
        [-1.0, -0.0, -0.0],
        Some(0x403130a3d70a3d71),
        Some(0x3fdf5ef6c18e8b09),
    ),
    (
        [1.8579964905283077, -0.2447198561902864, -0.9529337494572974],
        [-0.8660254037844398, 0.12436785483223059, 0.4842856973775163],
        Some(0x4002839ddd293062),
        Some(0x40002c75411f792c),
    ),
    (
        [2.0383802741745285, -0.5519513173267857, 0.38230352185267896],
        [
            -0.8660254037844398,
            0.4110321751788743,
            -0.28469729708538033,
        ],
        Some(0x4006767fbdbe193a),
        Some(0x40041faec3025c8c),
    ),
];

/// Points inside a corpus part's faces read the same, bit for bit, as
/// before: only a point in the band of a face's trim takes on the faces at
/// its nearest point, and a ball passes over only what its own faces hold
/// on its radius.
#[test]
fn points_inside_a_corpus_part_s_faces_read_as_before() {
    let path = format!(
        "{}/../../tests/corpus/m5x16_bhcs.step",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(path).expect("the corpus file is committed");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let screw = &import.solids[0];
    let points: Vec<Point> = INSIDE
        .iter()
        .map(|&([x, y, z], ..)| Point::new(x, y, z))
        .collect();
    let normals: Vec<Vector> = INSIDE
        .iter()
        .map(|&(_, [x, y, z], ..)| Vector::new(x, y, z))
        .collect();
    let limit = 30f64.to_radians();
    let cone = Method::Ray { cone: 0.2, rays: 5 };
    let rays =
        thickness_of_shape(model, screw, &points, &normals, &options(cone, limit, 50.0)).unwrap();
    let balls = thickness_of_shape(
        model,
        screw,
        &points,
        &normals,
        &options(Method::Sphere, limit, 50.0),
    )
    .unwrap();
    for (((p, ray), ball), &(_, _, was_ray, was_ball)) in
        points.iter().zip(rays).zip(balls).zip(INSIDE)
    {
        assert_eq!(ray.map(f64::to_bits), was_ray, "rays at {p:?}");
        assert_eq!(ball.map(f64::to_bits), was_ball, "ball at {p:?}");
    }
}
