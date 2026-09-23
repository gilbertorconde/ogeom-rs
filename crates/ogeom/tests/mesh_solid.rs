//! A triangle mesh becomes a B-rep solid: planar faces from its coplanar
//! regions, topology from its own connectivity, windings made to agree.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use std::time::{Duration, Instant};

use ogeom::algo::{MeshSolidOptions, check, solid_from_mesh, volume_properties};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, Triangulation, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A cube's twelve triangles as an STL holds them: every vertex repeated
/// once per triangle that uses it.
fn cube_soup(size: f64) -> Triangulation {
    let corner = |i: u32| {
        Point::new(
            f64::from(i & 1) * size,
            f64::from((i >> 1) & 1) * size,
            f64::from((i >> 2) & 1) * size,
        )
    };
    let triangles: [[u32; 3]; 12] = [
        [0, 2, 1],
        [1, 2, 3],
        [4, 5, 6],
        [5, 7, 6],
        [0, 1, 4],
        [1, 5, 4],
        [2, 6, 3],
        [3, 6, 7],
        [0, 4, 2],
        [2, 4, 6],
        [1, 3, 5],
        [3, 7, 5],
    ];
    soup(triangles.iter().map(|t| t.map(corner)))
}

fn soup(triangles: impl Iterator<Item = [Point; 3]>) -> Triangulation {
    let mut mesh = Triangulation::new();
    for t in triangles {
        let base = u32::try_from(mesh.positions.len()).unwrap();
        mesh.positions.extend(t);
        mesh.triangles.push([base, base + 1, base + 2]);
    }
    mesh
}

fn count(model: &Model, shape: &Shape, kind: ShapeType) -> usize {
    explore_unique(model, shape, kind).unwrap().len()
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

#[test]
fn an_stl_cube_becomes_a_six_faced_solid() {
    let mut model = Model::new();
    let out = solid_from_mesh(
        &mut model,
        &cube_soup(10.0),
        &MeshSolidOptions::default(),
        T,
    )
    .unwrap();
    assert!(out.closed);
    assert_eq!(model.kind_of(&out.shape).unwrap(), ShapeType::Solid);
    assert_eq!(count(&model, &out.shape, ShapeType::Face), 6);
    assert_eq!(count(&model, &out.shape, ShapeType::Edge), 12);
    assert_eq!(count(&model, &out.shape, ShapeType::Vertex), 8);
    assert_eq!(out.report.vertices_welded, 36 - 8);
    assert!((volume(&model, &out.shape) - 1000.0).abs() < 1e-6);
    let diagnosis = check(&model, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");

    // Unmerged, every triangle is a face and the solid is as valid.
    let mut model = Model::new();
    let options = MeshSolidOptions {
        merge_coplanar: false,
        ..MeshSolidOptions::default()
    };
    let out = solid_from_mesh(&mut model, &cube_soup(10.0), &options, T).unwrap();
    assert_eq!(count(&model, &out.shape, ShapeType::Face), 12);
    assert_eq!(count(&model, &out.shape, ShapeType::Edge), 18);
    assert!((volume(&model, &out.shape) - 1000.0).abs() < 1e-6);
    assert!(check(&model, &out.shape, T).unwrap().is_valid());
}

#[test]
fn an_open_mesh_comes_back_as_a_shell_with_its_holes_counted() {
    let mut mesh = cube_soup(10.0);
    mesh.triangles.pop();
    let mut model = Model::new();
    let out = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(!out.closed);
    assert_eq!(out.report.edges_used_once, 3);
    assert_eq!(model.kind_of(&out.shape).unwrap(), ShapeType::Shell);
}

/// Triangles wound every which way (some reversed against their
/// neighbours, then the whole inside out) are wound to agree and to face
/// outward, and the report counts what was turned.
#[test]
fn inconsistent_windings_are_turned_outward() {
    let mut mesh = cube_soup(10.0);
    for t in &mut mesh.triangles {
        t.swap(1, 2);
    }
    for t in mesh.triangles.iter_mut().step_by(3) {
        t.swap(1, 2);
    }
    let mut model = Model::new();
    let out = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(out.closed);
    assert_eq!(out.report.windings_flipped, 8);
    assert!((volume(&model, &out.shape) - 1000.0).abs() < 1e-6);
    assert!(check(&model, &out.shape, T).unwrap().is_valid());
}

/// A closed piece inside another is a void of the solid around it.
#[test]
fn a_piece_inside_another_is_a_void() {
    let inner = cube_soup(5.0);
    let shift = ogeom::math::Vector::new(2.5, 2.5, 2.5);
    let mut mesh = cube_soup(10.0);
    let base = u32::try_from(mesh.positions.len()).unwrap();
    mesh.positions
        .extend(inner.positions.iter().map(|p| *p + shift));
    mesh.triangles
        .extend(inner.triangles.iter().map(|t| t.map(|i| i + base)));
    let mut model = Model::new();
    let out = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(out.closed);
    assert_eq!(model.kind_of(&out.shape).unwrap(), ShapeType::Solid);
    assert_eq!(count(&model, &out.shape, ShapeType::Shell), 2);
    assert!((volume(&model, &out.shape) - 875.0).abs() < 1e-6);
    assert!(check(&model, &out.shape, T).unwrap().is_valid());
}

#[test]
fn a_converted_solid_takes_a_boolean() {
    let mut model = Model::new();
    let solid = solid_from_mesh(
        &mut model,
        &cube_soup(10.0),
        &MeshSolidOptions::default(),
        T,
    )
    .unwrap()
    .shape;
    let frame = Frame::new(Point::new(5.0, 5.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, frame, 2.0, 12.0, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(&mut model, &solid, &drill, T)
        .unwrap()
        .shape;
    let diagnosis = check(&model, &cut, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let v = volume_properties(&model, &cut, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass;
    let expected = 1000.0 - core::f64::consts::PI * 4.0 * 10.0;
    assert!(
        (v - expected).abs() / expected < 1e-3,
        "{v} against {expected}"
    );
}

fn drilled_block(model: &mut Model) -> Shape {
    let block = ogeom::algo::make_box(model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(10.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(model, frame, 4.0, 12.0, T)
        .unwrap()
        .shape;
    ogeom::boolean::cut(model, &block, &drill, T).unwrap().shape
}

fn meshed(model: &Model, shape: &Shape) -> Triangulation {
    ogeom::mesh::triangulate(model, shape, Deflection::with_chord(0.01).unwrap(), T).unwrap()
}

/// The kinds of surface a shape's faces are built on, counted.
fn kinds(model: &Model, shape: &Shape) -> [usize; 5] {
    use ogeom::geom::SurfaceGeometry as S;
    let mut out = [0; 5];
    for face in explore_unique(model, shape, ShapeType::Face).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        out[match model.geometry().surface(data.surface).unwrap() {
            S::Plane(_) => 0,
            S::Cylinder(_) => 1,
            S::Cone(_) => 2,
            S::Sphere(_) => 3,
            S::Torus(_) => 4,
            _ => panic!("a surface recognition does not build"),
        }] += 1;
    }
    out
}

/// A converted shape is valid and holds the volume the original does, to
/// within what the measurement's own meshing leaves.
fn holds(original: (&Model, &Shape), converted: (&Model, &Shape)) {
    let diagnosis = check(converted.0, converted.1, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let fine = Deflection::with_chord(1e-3).unwrap();
    let (a, b) = (
        volume_properties(original.0, original.1, fine, T)
            .unwrap()
            .mass,
        volume_properties(converted.0, converted.1, fine, T)
            .unwrap()
            .mass,
    );
    assert!((a - b).abs() / a < 2e-4, "{a} went in, {b} came out");
}

/// Without recognition, a drilled block's flat faces come back whole (the
/// drilled ones with the bore's polygon as a hole) and the bore faceted.
#[test]
fn planar_faces_with_holes_come_back_whole() {
    let mut model = Model::new();
    let drilled = drilled_block(&mut model);
    let mesh = ogeom::mesh::triangulate(&model, &drilled, Deflection::with_chord(0.05).unwrap(), T)
        .unwrap();
    let mut back = Model::new();
    let options = MeshSolidOptions {
        recognize: false,
        ..MeshSolidOptions::default()
    };
    let out = solid_from_mesh(&mut back, &mesh, &options, T).unwrap();
    assert!(out.closed, "{:?}", out.report);
    let faces = explore_unique(&back, &out.shape, ShapeType::Face).unwrap();
    let with_holes = faces
        .iter()
        .filter(|f| back.children_of(f).unwrap().len() == 2)
        .count();
    assert_eq!(with_holes, 2, "top and bottom keep the bore as a hole");
    assert!(faces.len() - 6 >= 8, "the bore's facets");
    let diagnosis = check(&back, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let (meshed, rebuilt) = (volume(&model, &drilled), volume(&back, &out.shape));
    assert!(
        (meshed - rebuilt).abs() / meshed < 1e-2,
        "{meshed} against {rebuilt}"
    );
}

/// With it, the bore is a cylinder again: seven faces, the bore one band
/// between the two circles it cuts from the top and bottom, and a seam.
#[test]
fn a_meshed_bore_comes_back_a_cylinder() {
    let mut model = Model::new();
    let drilled = drilled_block(&mut model);
    let mut back = Model::new();
    let out = solid_from_mesh(
        &mut back,
        &meshed(&model, &drilled),
        &MeshSolidOptions::default(),
        T,
    )
    .unwrap();
    assert!(out.closed);
    assert_eq!(kinds(&back, &out.shape), [6, 1, 0, 0, 0]);
    assert_eq!(out.report.curved_faces, 1);
    holds((&model, &drilled), (&back, &out.shape));
    let bore = explore_unique(&back, &out.shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let data = back.node(f).unwrap().data().as_face().unwrap();
            matches!(
                back.geometry().surface(data.surface).unwrap(),
                ogeom::geom::SurfaceGeometry::Cylinder(_)
            )
        })
        .unwrap();
    let data = back.node(&bore).unwrap().data().as_face().unwrap();
    let ogeom::geom::SurfaceGeometry::Cylinder(c) = back.geometry().surface(data.surface).unwrap()
    else {
        unreachable!()
    };
    assert!((c.cylinder().radius() - 4.0).abs() < 1e-9);

    // And takes a boolean as the original does.
    let frame = Frame::new(Point::new(0.0, 10.0, 5.0), Direction::X, Direction::Y, T).unwrap();
    let cross = ogeom::algo::make_cylinder(&mut back, frame, 2.0, 20.0, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(&mut back, &out.shape, &cross, T)
        .unwrap()
        .shape;
    assert!(check(&back, &cut, T).unwrap().is_valid());
}

/// Cylinders, cones, and the fillets and corner blends of a rounded box (
/// cylinders along the edges, spheres at the corners) come back on the
/// surfaces they were meshed from.
#[test]
fn primitives_and_fillets_come_back_on_their_surfaces() {
    let mut model = Model::new();
    let cylinder = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 12.0, T)
        .unwrap()
        .shape;
    let cone = ogeom::algo::make_cone(&mut model, Frame::WORLD, 6.0, 3.0, 10.0, T)
        .unwrap()
        .shape;
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let one = ogeom::fillet::fillet_edges(&mut model, &block, &edges[..1], 3.0, T)
        .unwrap()
        .shape;
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let rounded = ogeom::fillet::fillet_edges(&mut model, &block, &edges, 2.0, T)
        .unwrap()
        .shape;
    for (shape, expected) in [
        (&cylinder, [2, 1, 0, 0, 0]),
        (&cone, [2, 0, 1, 0, 0]),
        (&one, [6, 1, 0, 0, 0]),
        (&rounded, [6, 12, 0, 8, 0]),
    ] {
        let mut back = Model::new();
        let out = solid_from_mesh(
            &mut back,
            &meshed(&model, shape),
            &MeshSolidOptions::default(),
            T,
        )
        .unwrap();
        assert!(out.closed);
        assert_eq!(kinds(&back, &out.shape), expected);
        holds((&model, shape), (&back, &out.shape));
    }
}

/// Converts a shape's mesh and checks the surfaces it comes back on, its
/// validity and its volume. The volume is held to a thousandth: a face
/// built on a sphere or torus from boundary edges alone meshes a little
/// coarser at the measuring chord than the face it was meshed from.
fn comes_back_as(model: &Model, shape: &Shape, expected: [usize; 5]) {
    let mut back = Model::new();
    let out = solid_from_mesh(
        &mut back,
        &meshed(model, shape),
        &MeshSolidOptions::default(),
        T,
    )
    .unwrap();
    assert!(out.closed, "{:?}", out.report);
    assert_eq!(kinds(&back, &out.shape), expected);
    assert_eq!(out.report.curved_faceted, 0);
    let diagnosis = check(&back, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let fine = Deflection::with_chord(1e-3).unwrap();
    let (a, b) = (
        volume_properties(model, shape, fine, T).unwrap().mass,
        volume_properties(&back, &out.shape, fine, T).unwrap().mass,
    );
    assert!((a - b).abs() / a < 1e-3, "{a} went in, {b} came out");
}

/// A whole sphere and a whole torus have no boundary: each comes back one
/// face, closed on itself.
#[test]
fn whole_spheres_and_tori_come_back_one_face() {
    let mut model = Model::new();
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 7.0, T)
        .unwrap()
        .shape;
    let ring = ogeom::algo::make_torus(&mut model, Frame::WORLD, 10.0, 3.0, T)
        .unwrap()
        .shape;
    comes_back_as(&model, &ball, [0, 0, 0, 1, 0]);
    comes_back_as(&model, &ring, [0, 0, 0, 0, 1]);
}

/// A sphere cut by one plane is a cap, by two parallel planes a zone, and
/// a torus cut across its tube twice is a bent tube: each comes back on its
/// sphere or torus, bounded by full circles.
#[test]
fn caps_zones_and_bent_tubes_come_back_exact() {
    let mut model = Model::new();
    let rod = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 4.0, 10.0, T)
        .unwrap()
        .shape;
    let top = Frame::new(Point::new(0.0, 0.0, 10.0), Direction::Z, Direction::X, T).unwrap();
    let end = ogeom::algo::make_sphere(&mut model, top, 4.0, T)
        .unwrap()
        .shape;
    let pin = ogeom::boolean::fuse(&mut model, &rod, &end, T)
        .unwrap()
        .shape;

    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 7.0, T)
        .unwrap()
        .shape;
    let under = Frame::new(
        Point::new(-10.0, -10.0, -3.0),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let slab = ogeom::algo::make_box(&mut model, under, (20.0, 20.0, 5.0), T)
        .unwrap()
        .shape;
    let zone = ogeom::boolean::common(&mut model, &ball, &slab, T)
        .unwrap()
        .shape;

    let ring = ogeom::algo::make_torus(&mut model, Frame::WORLD, 10.0, 3.0, T)
        .unwrap()
        .shape;
    let turned = Direction::new(Vector::new(3.0, 1.0, 0.0), T).unwrap();
    let below = Frame::new(Point::new(0.0, 0.0, -5.0), Direction::Z, turned, T).unwrap();
    let quarter = ogeom::algo::make_box(&mut model, below, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let elbow = ogeom::boolean::common(&mut model, &ring, &quarter, T)
        .unwrap()
        .shape;

    comes_back_as(&model, &pin, [1, 1, 0, 1, 0]);
    comes_back_as(&model, &zone, [2, 0, 0, 1, 0]);
    comes_back_as(&model, &elbow, [2, 0, 0, 0, 1]);
}

/// A sphere bored through twice, across, meets its bores in four circles
/// that no single frame makes parallels: it is recognized, said to be
/// faceted, and still converts to a valid solid.
#[test]
fn what_cannot_be_built_exactly_stays_faceted() {
    let mut model = Model::new();
    let mut shape = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 7.0, T)
        .unwrap()
        .shape;
    for (at, along, across, radius) in [
        (Point::new(0.0, 0.0, -10.0), Direction::Z, Direction::X, 2.0),
        (Point::new(-10.0, 0.0, 0.0), Direction::X, Direction::Y, 1.5),
    ] {
        let frame = Frame::new(at, along, across, T).unwrap();
        let bore = ogeom::algo::make_cylinder(&mut model, frame, radius, 20.0, T)
            .unwrap()
            .shape;
        shape = ogeom::boolean::cut(&mut model, &shape, &bore, T)
            .unwrap()
            .shape;
    }
    let mut back = Model::new();
    let out = solid_from_mesh(
        &mut back,
        &meshed(&model, &shape),
        &MeshSolidOptions::default(),
        T,
    )
    .unwrap();
    assert!(out.closed);
    assert!(out.report.curved_faceted >= 1);
    assert_eq!(kinds(&back, &out.shape)[3], 0, "the sphere is not built");
    assert!(check(&back, &out.shape, T).unwrap().is_valid());
}

/// A torus tessellated into 200 000 triangles converts in seconds, to the
/// one face of the torus it is.
#[test]
fn a_large_mesh_converts_in_seconds() {
    let (rings, sides) = (500_u32, 200_u32);
    let (major, minor) = (40.0, 10.0);
    let mut mesh = Triangulation::new();
    for i in 0..rings {
        let u = f64::from(i) / f64::from(rings) * core::f64::consts::TAU;
        for j in 0..sides {
            let v = f64::from(j) / f64::from(sides) * core::f64::consts::TAU;
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
    assert_eq!(mesh.triangles.len(), 200_000);
    let mut model = Model::new();
    let started = Instant::now();
    let out = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    let took = started.elapsed();
    assert!(out.closed);
    assert_eq!(out.report.triangles, 200_000);
    assert_eq!(out.report.faces, 1);
    assert!(took < Duration::from_secs(10), "{took:?}");
}
