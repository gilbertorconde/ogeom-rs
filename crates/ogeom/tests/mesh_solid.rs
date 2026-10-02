//! A triangle mesh becomes a B-rep solid: planar faces from its coplanar
//! regions, topology from its own connectivity, windings made to agree.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use std::time::{Duration, Instant};

use ogeom::algo::{
    MeshSolidOptions, Severity, check, solid_from_mesh, tight_bounds, volume_properties,
};
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

/// A cube whose bottom square is split two ways across its diagonal: one
/// triangle along the whole diagonal, two meeting at its middle, and a flat
/// sliver along the diagonal sealing the difference, as exporters leave at a
/// T-junction. The sliver has no area and is dropped. The triangle across
/// the diagonal is split at its middle instead, so the solid stays closed.
#[test]
fn a_flat_sliver_sealing_a_t_junction_leaves_the_cube_closed() {
    let size = 10.0;
    let corner = |i: u32| {
        Point::new(
            f64::from(i & 1) * size,
            f64::from((i >> 1) & 1) * size,
            f64::from((i >> 2) & 1) * size,
        )
    };
    let middle = Point::new(size / 2.0, size / 2.0, 0.0);
    let mut triangles: Vec<[Point; 3]> = [
        [0, 2, 1],
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
    ]
    .iter()
    .map(|t: &[u32; 3]| t.map(corner))
    .collect();
    triangles.push([corner(1), middle, corner(3)]);
    triangles.push([middle, corner(2), corner(3)]);
    triangles.push([corner(1), corner(2), middle]);
    let mut model = Model::new();
    let out = solid_from_mesh(
        &mut model,
        &soup(triangles.into_iter()),
        &MeshSolidOptions::default(),
        T,
    )
    .unwrap();
    assert_eq!(out.report.degenerate_dropped, 1);
    assert!(out.closed, "{:?}", out.report);
    assert!(check(&model, &out.shape, T).unwrap().is_valid());
    assert!((volume(&model, &out.shape) - 1000.0).abs() < 1e-6);
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
    comes_back_from(model, shape, &meshed(model, shape), expected);
}

fn comes_back_from(model: &Model, shape: &Shape, mesh: &Triangulation, expected: [usize; 5]) {
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, mesh, &MeshSolidOptions::default(), T).unwrap();
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

/// A mesh as coarse as a printer's STL export (the default deflection, a
/// tenth of a millimetre and half a radian, and a chord five times that)
/// spans a two millimetre fillet in three or four facets, some of them
/// across two of the tessellator's rows at once. The fillets, the corner
/// balls and a torus round a hole's edge still come back on their
/// surfaces.
#[test]
fn a_coarse_mesh_keeps_its_fillets() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let rounded = ogeom::fillet::fillet_edges(&mut model, &block, &edges, 2.0, T)
        .unwrap()
        .shape;
    let drilled = drilled_block(&mut model);
    let rim = explore_unique(&model, &drilled, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|edge| {
            let data = model.node(edge).unwrap().data().as_edge().unwrap();
            let Some(ogeom::topo::EdgeRepr::Curve3d { curve, .. }) = data.curve3d() else {
                return false;
            };
            matches!(
                model.geometry().curve(*curve),
                Some(ogeom::geom::Curve::Circle(_))
            )
        })
        .unwrap();
    let eased = ogeom::fillet::fillet_edges(&mut model, &drilled, &[rim], 1.0, T)
        .unwrap()
        .shape;
    for deflection in [
        Deflection::default(),
        Deflection {
            chord: 0.5,
            ..Deflection::default()
        },
    ] {
        for (shape, expected) in [(&rounded, [6, 12, 0, 8, 0]), (&eased, [6, 1, 0, 0, 1])] {
            let mesh = ogeom::mesh::triangulate(&model, shape, deflection, T).unwrap();
            comes_back_from(&model, shape, &mesh, expected);
        }
    }
}

/// What a coarse mesh comes back as does not hang on the order its
/// triangles are listed in: a fillet two facets wide, sampled before the
/// corners beside it are claimed, fits nothing then and is taken again.
#[test]
fn a_coarse_rounded_block_comes_back_whatever_its_triangle_order() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let rounded = ogeom::fillet::fillet_edges(&mut model, &block, &edges, 2.0, T)
        .unwrap()
        .shape;
    let coarse = Deflection {
        chord: 0.5,
        ..Deflection::default()
    };
    let mesh = ogeom::mesh::triangulate(&model, &rounded, coarse, T).unwrap();
    let n = mesh.triangles.len();
    for turn in [0, n / 7, n / 3, n / 2, 2 * n / 3] {
        let mut turned = mesh.clone();
        turned.triangles.rotate_left(turn);
        comes_back_from(&model, &rounded, &turned, [6, 12, 0, 8, 0]);
    }
    let mut reversed = mesh;
    reversed.triangles.reverse();
    comes_back_from(&model, &rounded, &reversed, [6, 12, 0, 8, 0]);
}

/// Two cylinders meeting along a curve that is no circle and no line: a
/// bar drilled across, and a pipe with a branch. Each meeting is fitted
/// once on both surfaces, and each wall is built round its axis with a seam
/// of its own clear of the holes in it, at a fine mesh and a printer's.
#[test]
fn crossing_cylinders_meet_along_a_fitted_curve() {
    let mut model = Model::new();
    let across = |p: (f64, f64, f64), z: Direction, x: Direction| {
        Frame::new(Point::new(p.0, p.1, p.2), z, x, T).unwrap()
    };
    let bar = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T)
        .unwrap()
        .shape;
    let drill = ogeom::algo::make_cylinder(
        &mut model,
        across((0.0, -10.0, 10.0), Direction::Y, Direction::X),
        1.5,
        20.0,
        T,
    )
    .unwrap()
    .shape;
    let drilled = ogeom::boolean::cut(&mut model, &bar, &drill, T)
        .unwrap()
        .shape;
    let main = ogeom::algo::make_cylinder(
        &mut model,
        across((-15.0, 0.0, 0.0), Direction::X, Direction::Y),
        4.0,
        30.0,
        T,
    )
    .unwrap()
    .shape;
    let branch = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 2.5, 12.0, T)
        .unwrap()
        .shape;
    let tee = ogeom::boolean::fuse(&mut model, &main, &branch, T)
        .unwrap()
        .shape;
    for deflection in [Deflection::with_chord(0.01).unwrap(), Deflection::default()] {
        for (shape, expected) in [(&drilled, [2, 2, 0, 0, 0]), (&tee, [3, 2, 0, 0, 0])] {
            let mesh = ogeom::mesh::triangulate(&model, shape, deflection, T).unwrap();
            comes_back_from(&model, shape, &mesh, expected);
        }
    }
}

/// A ball or a ring whose boundary only makes holes in it: a boss fused
/// on a ball, a ball bored twice across, a ring pierced through its tube
/// and through its outer equator. Each comes back as the whole surface
/// with the holes as inner wires, its seams and poles turned clear of
/// them.
#[test]
fn balls_and_rings_with_holes_come_back_whole_but_for_them() {
    let mut model = Model::new();
    let at = |p: (f64, f64, f64), z: Direction, x: Direction| {
        Frame::new(Point::new(p.0, p.1, p.2), z, x, T).unwrap()
    };
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 6.0, T)
        .unwrap()
        .shape;
    let boss = ogeom::algo::make_cylinder(
        &mut model,
        at((2.0, 1.0, 0.0), Direction::Z, Direction::X),
        1.5,
        9.0,
        T,
    )
    .unwrap()
    .shape;
    let knob = ogeom::boolean::fuse(&mut model, &ball, &boss, T)
        .unwrap()
        .shape;
    let mut bored = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 7.0, T)
        .unwrap()
        .shape;
    for (from, along, across, radius) in [
        ((0.0, 0.0, -10.0), Direction::Z, Direction::X, 2.0),
        ((-10.0, 0.0, 0.0), Direction::X, Direction::Y, 1.5),
    ] {
        let bore = ogeom::algo::make_cylinder(&mut model, at(from, along, across), radius, 20.0, T)
            .unwrap()
            .shape;
        bored = ogeom::boolean::cut(&mut model, &bored, &bore, T)
            .unwrap()
            .shape;
    }
    let ring = ogeom::algo::make_torus(&mut model, Frame::WORLD, 10.0, 3.0, T)
        .unwrap()
        .shape;
    let pin = ogeom::algo::make_cylinder(
        &mut model,
        at((10.0, 0.0, -10.0), Direction::Z, Direction::X),
        1.0,
        20.0,
        T,
    )
    .unwrap()
    .shape;
    let pierced = ogeom::boolean::cut(&mut model, &ring, &pin, T)
        .unwrap()
        .shape;
    let radial = ogeom::algo::make_cylinder(
        &mut model,
        at((5.0, 0.0, 0.0), Direction::X, Direction::Y),
        1.0,
        10.0,
        T,
    )
    .unwrap()
    .shape;
    let through = ogeom::boolean::cut(&mut model, &ring, &radial, T)
        .unwrap()
        .shape;
    // A ball less the corner two square planes cut off: one loop of two
    // arcs, which leaves the ball's larger side.
    let cornered = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 7.0, T)
        .unwrap()
        .shape;
    let above = ogeom::algo::make_box(
        &mut model,
        at((-20.0, -20.0, 2.0), Direction::Z, Direction::X),
        (40.0, 40.0, 40.0),
        T,
    )
    .unwrap()
    .shape;
    let beside = ogeom::algo::make_box(
        &mut model,
        at((3.0, -20.0, -20.0), Direction::X, Direction::Y),
        (40.0, 40.0, 40.0),
        T,
    )
    .unwrap()
    .shape;
    let cornered = ogeom::boolean::cut(&mut model, &cornered, &above, T)
        .unwrap()
        .shape;
    let cornered = ogeom::boolean::cut(&mut model, &cornered, &beside, T)
        .unwrap()
        .shape;
    for (shape, expected) in [
        (&knob, [1, 1, 0, 1, 0]),
        (&bored, [0, 3, 0, 1, 0]),
        (&pierced, [0, 1, 0, 0, 1]),
        (&through, [0, 1, 0, 0, 1]),
        (&cornered, [2, 0, 0, 1, 0]),
    ] {
        comes_back_as(&model, shape, expected);
        let mesh = ogeom::mesh::triangulate(&model, shape, Deflection::default(), T).unwrap();
        comes_back_from(&model, shape, &mesh, expected);
    }
}

/// A chamfered hole whose mesh has one vertex of the chamfer's rim a
/// little off the cone (pushed outward along the plate, which still holds
/// it): the triangles on it cannot join the cone, and the chamfer's band is
/// cut open along a slit. It comes back a cone all the same, a patch
/// running round from one side of the slit to the other.
#[test]
fn a_chamfer_with_a_stray_vertex_is_still_a_cone() {
    let mut model = Model::new();
    let at = |z: f64| Frame::new(Point::new(10.0, 10.0, z), Direction::Z, Direction::X, T).unwrap();
    let plate = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 5.0), T)
        .unwrap()
        .shape;
    let sink = ogeom::algo::make_cone(&mut model, at(4.0), 2.0, 3.0, 1.0, T)
        .unwrap()
        .shape;
    let bore = ogeom::algo::make_cylinder(&mut model, at(-1.0), 2.0, 7.0, T)
        .unwrap()
        .shape;
    let part = ogeom::boolean::cut(&mut model, &plate, &sink, T)
        .unwrap()
        .shape;
    let part = ogeom::boolean::cut(&mut model, &part, &bore, T)
        .unwrap()
        .shape;
    let mut mesh = meshed(&model, &part);
    let k = mesh
        .positions
        .iter()
        .position(|p| {
            (p.z - 5.0).abs() < 1e-9 && ((p.x - 10.0).hypot(p.y - 10.0) - 3.0).abs() < 1e-9
        })
        .unwrap();
    let p = mesh.positions[k];
    mesh.positions[k] = p + Vector::new(p.x - 10.0, p.y - 10.0, 0.0) * (1e-4 / 3.0);
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(out.closed);
    assert_eq!(kinds(&back, &out.shape)[2], 1, "the chamfer is a cone");
    assert_eq!(out.report.curved_faceted, 0);
    assert!(check(&back, &out.shape, T).unwrap().is_valid());
}

/// A plate bored through whose mesh has one vertex of the bore's rim a
/// little off the cylinder: the bore comes back as one patch running round
/// from one side of a slit to the other, its chart window starting at the
/// slit rather than at the cylinder's zero. Holes drilled across the plate
/// through the bore cut and fill valid and share the plate's volume, their
/// sections on the bore placed on the patch's side of the slit. The last
/// passes within the slit's loose tolerance of it without meeting it, and
/// its rim is still a hole in the bore.
#[test]
fn holes_across_a_bore_opened_along_a_slit_cut_and_fill_valid() {
    let mut model = Model::new();
    let at = |z: f64| Frame::new(Point::new(10.0, 10.0, z), Direction::Z, Direction::X, T).unwrap();
    let plate = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 5.0), T)
        .unwrap()
        .shape;
    let bore = ogeom::algo::make_cylinder(&mut model, at(-1.0), 2.0, 7.0, T)
        .unwrap()
        .shape;
    let part = ogeom::boolean::cut(&mut model, &plate, &bore, T)
        .unwrap()
        .shape;
    let mut mesh = meshed(&model, &part);
    let k = mesh
        .positions
        .iter()
        .position(|p| {
            (p.z - 5.0).abs() < 1e-9
                && ((p.x - 10.0).hypot(p.y - 10.0) - 2.0).abs() < 1e-9
                && p.x > 11.0
                && p.y > 10.5
        })
        .unwrap();
    let p = mesh.positions[k];
    mesh.positions[k] = p + Vector::new(p.x - 10.0, p.y - 10.0, 0.0) * (1e-4 / 2.0);
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
    let bores: Vec<Shape> = explore_unique(&back, &out.shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            let data = back.node(f).unwrap().data().as_face().unwrap();
            matches!(
                back.geometry().surface(data.surface).unwrap(),
                ogeom::geom::SurfaceGeometry::Cylinder(_)
            )
        })
        .collect();
    assert_eq!(bores.len(), 1);
    assert_eq!(
        back.children_of(&bores[0]).unwrap().len(),
        1,
        "a patch, not a band"
    );
    let fine = Deflection::with_chord(0.01).unwrap();
    let whole = volume_properties(&back, &out.shape, fine, T).unwrap().mass;
    for (y, z, r) in [
        (8.7, 2.5, 0.8),
        (10.0, 1.2, 0.4),
        (10.6, 3.7, 0.8),
        (11.2, 2.5, 0.4),
    ] {
        let across = Frame::new(Point::new(-5.0, y, z), Direction::X, Direction::Y, T).unwrap();
        let drill = ogeom::algo::make_cylinder(&mut back, across, r, 30.0, T)
            .unwrap()
            .shape;
        let mut shares = 0.0;
        for (name, made) in [
            ("cut", ogeom::boolean::cut(&mut back, &out.shape, &drill, T)),
            (
                "common",
                ogeom::boolean::common(&mut back, &out.shape, &drill, T),
            ),
        ] {
            let made = made.unwrap_or_else(|e| panic!("{name} at ({y}, {z}): {e}"));
            let diagnosis = check(&back, &made.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{name} at ({y}, {z}): {diagnosis}");
            shares += volume_properties(&back, &made.shape, fine, T).unwrap().mass;
        }
        assert!(
            (shares - whole).abs() < whole * 5e-4,
            "at ({y}, {z}): {shares} against {whole}"
        );
    }
}

/// Each face of a converted shape with its surface, and the pairs of
/// faces sharing an edge.
fn surfaces_and_neighbours(
    model: &Model,
    shape: &Shape,
) -> (Vec<ogeom::geom::SurfaceGeometry>, Vec<(usize, usize)>) {
    let faces = explore_unique(model, shape, ShapeType::Face).unwrap();
    let surfaces = faces
        .iter()
        .map(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            model.geometry().surface(data.surface).unwrap().clone()
        })
        .collect();
    let edges: Vec<Vec<Shape>> = faces
        .iter()
        .map(|f| explore_unique(model, f, ShapeType::Edge).unwrap())
        .collect();
    let mut pairs = Vec::new();
    for i in 0..faces.len() {
        for j in i + 1..faces.len() {
            if edges[i]
                .iter()
                .any(|e| edges[j].iter().any(|f| e.is_same(f)))
            {
                pairs.push((i, j));
            }
        }
    }
    (surfaces, pairs)
}

/// Fillets and corner balls come back on the surfaces their neighbours
/// fix, not on free fits, so they meet them tangentially to rounding: from
/// a single-precision mesh, each corner ball of a rounded block is centred
/// on the axes of the fillets it meets with their radius, and the torus
/// easing a bore's rim stands a tube radius off the top face and off the
/// bore.
#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
#[test]
fn fillets_and_corner_balls_meet_their_neighbours_tangentially() {
    use ogeom::geom::SurfaceGeometry as S;
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let rounded = ogeom::fillet::fillet_edges(&mut model, &block, &edges, 2.0, T)
        .unwrap()
        .shape;
    let drilled = drilled_block(&mut model);
    let rim = explore_unique(&model, &drilled, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|edge| {
            let data = model.node(edge).unwrap().data().as_edge().unwrap();
            let Some(ogeom::topo::EdgeRepr::Curve3d { curve, .. }) = data.curve3d() else {
                return false;
            };
            matches!(
                model.geometry().curve(*curve),
                Some(ogeom::geom::Curve::Circle(_))
            )
        })
        .unwrap();
    let eased = ogeom::fillet::fillet_edges(&mut model, &drilled, &[rim], 1.0, T)
        .unwrap()
        .shape;
    let (mut balls, mut tubes) = (0, 0);
    for shape in [&rounded, &eased] {
        let mut mesh = meshed(&model, shape);
        for p in &mut mesh.positions {
            *p = Point::new(
                f64::from(p.x as f32),
                f64::from(p.y as f32),
                f64::from(p.z as f32),
            );
        }
        let mut back = Model::new();
        let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
        assert!(check(&back, &out.shape, T).unwrap().is_valid());
        let (surfaces, pairs) = surfaces_and_neighbours(&back, &out.shape);
        for (i, j) in pairs {
            for (a, b) in [(&surfaces[i], &surfaces[j]), (&surfaces[j], &surfaces[i])] {
                match (a, b) {
                    (S::Sphere(s), S::Cylinder(c)) => {
                        let (s, c) = (s.sphere(), c.cylinder());
                        let z = c.frame().z().vector();
                        let w = s.centre() - c.frame().origin();
                        assert!((w - z * w.dot(z)).magnitude() < 1e-9);
                        assert!((s.radius() - c.radius()).abs() < 1e-9);
                        balls += 1;
                    }
                    (S::Torus(t), S::Plane(p)) => {
                        let (t, p) = (t.torus(), p.plane());
                        let n = p.normal().vector();
                        let height = (t.centre() - p.frame().origin()).dot(n).abs();
                        assert!(t.frame().z().vector().cross(n).magnitude() < 1e-9);
                        assert!((height - t.minor_radius()).abs() < 1e-9);
                        tubes += 1;
                    }
                    (S::Torus(t), S::Cylinder(c)) => {
                        let (t, c) = (t.torus(), c.cylinder());
                        let off = (t.major_radius() - c.radius()).abs();
                        assert!((off - t.minor_radius()).abs() < 1e-9);
                        tubes += 1;
                    }
                    _ => {}
                }
            }
        }
    }
    assert_eq!(balls, 24, "eight corner balls, each meeting three fillets");
    assert_eq!(tubes, 2, "the rim's torus, on the top face and the bore");
}

/// A rounded box whose corner balls are roughened into free-form facets,
/// as a mesh.
fn rough_rounded_box() -> Triangulation {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let rounded = ogeom::fillet::fillet_edges(&mut model, &block, &edges, 2.0, T)
        .unwrap()
        .shape;
    let mut mesh = meshed(&model, &rounded);
    for (i, p) in mesh.positions.iter_mut().enumerate() {
        let near = |v: f64, low: f64, high: f64| if v < (low + high) / 2.0 { low } else { high };
        let centre = Point::new(
            near(p.x, 2.0, 18.0),
            near(p.y, 2.0, 18.0),
            near(p.z, 2.0, 8.0),
        );
        let inside = |v: f64, c: f64, mid: f64| {
            if c < mid { v < c - 1e-3 } else { v > c + 1e-3 }
        };
        if inside(p.x, centre.x, 10.0) && inside(p.y, centre.y, 10.0) && inside(p.z, centre.z, 5.0)
        {
            let d = *p - centre;
            #[allow(clippy::cast_precision_loss, reason = "a vertex index")]
            let bump = 0.02 * (i as f64 * 1.7).sin();
            *p += d / d.magnitude() * bump;
        }
    }
    mesh
}

/// Each fillet of the rough box meets its corner all but tangentially,
/// along a chain no curve both surfaces share can be solved for. The
/// fillets are still cylinders, bounded there by curves threaded through
/// the chain's own vertices, good to how far they stray.
#[test]
fn fillets_ending_on_rough_corners_are_still_cylinders() {
    let mut back = Model::new();
    let mesh = rough_rounded_box();
    let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(out.closed);
    assert_eq!(kinds(&back, &out.shape)[1], 12, "every fillet a cylinder");
    assert_eq!(out.report.curved_faceted, 0);
    assert!(check(&back, &out.shape, T).unwrap().is_valid());
    // A single facet left beside a fillet is thinner than the curve the
    // two surfaces meet along bulges; their seam is threaded straight
    // instead, so the facet's trim does not fold over itself and the
    // whole shape tessellates closed.
    let drawn = ogeom::mesh::triangulate(&back, &out.shape, Deflection::default(), T).unwrap();
    assert!(drawn.is_closed());
    let got = volume(&back, &out.shape);
    assert!(
        (got - mesh.volume()).abs() < mesh.volume() * 2e-3,
        "{got} against {}",
        mesh.volume()
    );
}

/// Holes drilled through the rough box beside a corner cross planes,
/// fillets and facets bounded by threaded curves a few hundredths loose.
/// That looseness is a length in space: it welds a section's ends in each
/// chart only as far as it reaches there, and it collapses no section or
/// edge that is longer than its own doubt. Each hole cuts and fills valid.
#[test]
fn holes_beside_a_rough_corner_cut_and_fill_valid() {
    let mut model = Model::new();
    let out = solid_from_mesh(
        &mut model,
        &rough_rounded_box(),
        &MeshSolidOptions::default(),
        T,
    )
    .unwrap();
    let whole = volume(&model, &out.shape);
    for (x, y, r) in [(1.1, 0.6, 1.4), (2.05, 0.35, 1.0), (1.35, 1.9, 1.4)] {
        let at = Frame::new(Point::new(x, y, -5.0), Direction::Z, Direction::X, T).unwrap();
        let drill = ogeom::algo::make_cylinder(&mut model, at, r, 20.0, T)
            .unwrap()
            .shape;
        let mut shares = 0.0;
        for (name, made) in [
            (
                "cut",
                ogeom::boolean::cut(&mut model, &out.shape, &drill, T),
            ),
            (
                "common",
                ogeom::boolean::common(&mut model, &out.shape, &drill, T),
            ),
        ] {
            let made = made.unwrap_or_else(|e| panic!("{name} at ({x}, {y}): {e}"));
            let diagnosis = check(&model, &made.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{name} at ({x}, {y}): {diagnosis}");
            shares += volume(&model, &made.shape);
        }
        assert!(
            (shares - whole).abs() < whole * 2e-3,
            "at ({x}, {y}): {shares} against {whole}"
        );
    }
}

/// Holes drilled sideways past the rough box's corner, across faces whose
/// threaded edges set a weld of a few hundredths. A hole's rim is cut where
/// it crosses the hole's own seam, a few microns from where the rim starts,
/// and the sliver between collapses in both faces' charts. It collapses in
/// space as well, or the rest of the rim ends on two vertices and its ring
/// stays open.
#[test]
fn holes_sideways_past_a_rough_corner_cut_and_fill_valid() {
    let mut model = Model::new();
    let out = solid_from_mesh(
        &mut model,
        &rough_rounded_box(),
        &MeshSolidOptions::default(),
        T,
    )
    .unwrap();
    let fine = Deflection::with_chord(0.01).unwrap();
    let whole = volume_properties(&model, &out.shape, fine, T).unwrap().mass;
    let along_x =
        |y: f64, z: f64| Frame::new(Point::new(-5.0, y, z), Direction::X, Direction::Y, T).unwrap();
    let along_y =
        |x: f64, z: f64| Frame::new(Point::new(x, -5.0, z), Direction::Y, Direction::Z, T).unwrap();
    for (at, r) in [
        (along_x(0.65, 2.66), 0.64),
        (along_y(3.59, 2.32), 0.61),
        (along_y(2.19, 2.82), 0.27),
    ] {
        let drill = ogeom::algo::make_cylinder(&mut model, at, r, 30.0, T)
            .unwrap()
            .shape;
        let mut shares = 0.0;
        for (name, made) in [
            (
                "cut",
                ogeom::boolean::cut(&mut model, &out.shape, &drill, T),
            ),
            (
                "common",
                ogeom::boolean::common(&mut model, &out.shape, &drill, T),
            ),
        ] {
            let made = made.unwrap_or_else(|e| panic!("{name} at {:?}: {e}", at.origin()));
            let diagnosis = check(&model, &made.shape, T).unwrap();
            assert!(
                diagnosis.is_valid(),
                "{name} at {:?}: {diagnosis}",
                at.origin()
            );
            shares += volume_properties(&model, &made.shape, fine, T)
                .unwrap()
                .mass;
        }
        assert!(
            (shares - whole).abs() < whole * 5e-4,
            "at {:?}: {shares} against {whole}",
            at.origin()
        );
    }
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

/// A rounded rectangle `w` by `h`, corner radius `r` in `n` steps per
/// corner, drawn in by `inset` with its corners' centres kept.
fn rounded_outline(w: f64, h: f64, r: f64, n: u32, inset: f64) -> Vec<(f64, f64)> {
    let r = r - inset;
    let mut out = Vec::new();
    for (cx, cy, start) in [
        (w - r - inset, r + inset, -0.5),
        (w - r - inset, h - r - inset, 0.0),
        (r + inset, h - r - inset, 0.5),
        (r + inset, r + inset, 1.0),
    ] {
        for k in 0..=n {
            let a = core::f64::consts::PI * (start + 0.5 * f64::from(k) / f64::from(n));
            out.push((cx + r * a.cos(), cy + r * a.sin()));
        }
    }
    out
}

/// A rounded-rectangle slab from z = 0 to `top`: its walls stand between
/// z = `fillet` and `top - fillet`, each turned into its cap by a single
/// row of triangles, and its long sides are one quad each.
fn rounded_slab(r: f64, n: u32, fillet: f64, top: f64) -> Triangulation {
    let (w, h) = (47.0, 21.5);
    let rings = [
        (rounded_outline(w, h, r, n, fillet), 0.0),
        (rounded_outline(w, h, r, n, 0.0), fillet),
        (rounded_outline(w, h, r, n, 0.0), top - fillet),
        (rounded_outline(w, h, r, n, fillet), top),
    ];
    let m = u32::try_from(rings[0].0.len()).unwrap();
    let mut t = Triangulation::new();
    for (ring, z) in &rings {
        for &(x, y) in ring {
            t.positions.push(Point::new(x, y, *z));
        }
    }
    for k in 0..3 {
        for i in 0..m {
            let j = (i + 1) % m;
            let (a, b, c, d) = (k * m + i, k * m + j, (k + 1) * m + j, (k + 1) * m + i);
            t.triangles.push([a, b, c]);
            t.triangles.push([a, c, d]);
        }
    }
    let (bottom, cap) = (4 * m, 4 * m + 1);
    t.positions.push(Point::new(w / 2.0, h / 2.0, 0.0));
    t.positions.push(Point::new(w / 2.0, h / 2.0, top));
    for i in 0..m {
        let j = (i + 1) % m;
        t.triangles.push([bottom, j, i]);
        t.triangles.push([cap, 3 * m + i, 3 * m + j]);
    }
    t
}

/// A slab whose long walls are single facets between its rounded corners:
/// a surface through both corners' vertices passes through the long
/// facets' corners too, and must still not be taken for them. Every face
/// recognized stays on the triangles it replaced, the solid facing out and
/// bounded as the mesh is.
#[test]
fn recognized_faces_stay_on_the_triangles_they_replace() {
    for (r, n, fillet) in [
        (8.0, 19, 0.5),
        (8.0, 26, 0.3),
        (8.0, 32, 0.5),
        (8.0, 8, 0.5),
        (2.0, 3, 0.5),
    ] {
        let top = 5.0;
        let mesh = rounded_slab(r, n, fillet, top);
        let mut model = Model::new();
        let built = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
        let at = format!("corner {r} in {n} steps, fillet {fillet}");
        let diagnosis = check(&model, &built.shape, T).unwrap();
        assert!(
            diagnosis.of(Severity::Broken).is_empty(),
            "{at}: {diagnosis}"
        );
        let fine = ogeom::mesh::triangulate(
            &model,
            &built.shape,
            Deflection::with_chord(0.01).unwrap(),
            T,
        )
        .unwrap();
        let bounds = tight_bounds(&model, &built.shape, T).unwrap();
        let (lo, hi) = (bounds.low().unwrap(), bounds.high().unwrap());
        for p in fine.positions.iter().chain([&lo, &hi]) {
            assert!(
                p.z > -1e-4 && p.z < top + 1e-4 && p.y > -1e-4 && p.y < 21.5 + 1e-4,
                "{at}: {p:?} outside the slab"
            );
        }
    }
}

/// A slab over a rounded outline whose points are rounded to single
/// precision, as an STL stores them, from z = 0 to `top`.
#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
fn stl_slab(r: f64, n: u32, top: f64) -> (Vec<(f64, f64)>, Triangulation) {
    let ring: Vec<(f64, f64)> = rounded_outline(47.0, 21.5, r, n, 0.0)
        .into_iter()
        .map(|(x, y)| (f64::from(x as f32), f64::from(y as f32)))
        .collect();
    let count = u32::try_from(ring.len()).unwrap();
    let mut t = Triangulation::new();
    for z in [0.0, top] {
        for &(x, y) in &ring {
            t.positions.push(Point::new(x, y, z));
        }
    }
    t.positions.push(Point::new(23.5, 10.75, 0.0));
    t.positions.push(Point::new(23.5, 10.75, top));
    let (bottom, cap) = (2 * count, 2 * count + 1);
    for i in 0..count {
        let j = (i + 1) % count;
        t.triangles.push([i, j, count + j]);
        t.triangles.push([i, count + j, count + i]);
        t.triangles.push([bottom, j, i]);
        t.triangles.push([cap, count + i, count + j]);
    }
    (ring, t)
}

/// A pad sketched on a converted slab's own outline and sunk into it: its
/// sides run along the slab's flat walls and are chords of its recognized
/// corners, touching them at every vertex. The union is the slab and the
/// intersection the pad. A pad standing out through the top is never taken
/// for either.
#[test]
fn a_pad_within_a_converted_solid_fuses_to_the_solid() {
    for (r, n) in [(8.0, 20), (2.0, 12), (3.0, 7)] {
        let at = format!("corner {r} in {n} steps");
        let mut model = Model::new();
        let (ring, mesh) = stl_slab(r, n, 4.0);
        let base = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T)
            .unwrap()
            .shape;
        let pad = |model: &mut Model, height: f64| {
            let points: Vec<Point> = ring.iter().map(|&(x, y)| Point::new(x, y, 1.0)).collect();
            let wire = ogeom::algo::make_polygon(model, &points, true, T)
                .unwrap()
                .shape;
            let plane = ogeom::math::Plane::new(
                Frame::new(Point::new(0.0, 0.0, 1.0), Direction::Z, Direction::X, T).unwrap(),
            );
            let face = ogeom::algo::make_face(
                model,
                ogeom::geom::PlaneSurface::new(plane).into(),
                &[wire],
                T,
            )
            .unwrap()
            .shape;
            ogeom::algo::make_prism(model, &face, Vector::new(0.0, 0.0, height), T)
                .unwrap()
                .shape
        };
        let inside = pad(&mut model, 3.0);
        let (whole, sunk) = (volume(&model, &base), volume(&model, &inside));
        let fused = ogeom::boolean::fuse(&mut model, &base, &inside, T)
            .unwrap_or_else(|e| panic!("{at}: {e}"))
            .shape;
        assert!(check(&model, &fused, T).unwrap().is_valid(), "{at}");
        let v = volume(&model, &fused);
        assert!(
            (v - whole).abs() < whole * 1e-9,
            "{at}: {v} against {whole}"
        );
        let common = ogeom::boolean::common(&mut model, &base, &inside, T)
            .unwrap_or_else(|e| panic!("{at}: {e}"))
            .shape;
        let v = volume(&model, &common);
        assert!((v - sunk).abs() < sunk * 1e-9, "{at}: {v} against {sunk}");

        let proud = pad(&mut model, 3.5);
        if let Ok(fused) = ogeom::boolean::fuse(&mut model, &base, &proud, T) {
            let v = volume(&model, &fused.shape);
            assert!(v > whole * (1.0 + 1e-3), "{at}: a proud pad fused to {v}");
        }
    }
}

/// A slab's points rounded to single precision, and moved by up to `noise`
/// along each axis in a fixed pseudo-random pattern, as an exporter that
/// rounds its own vertices leaves them.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    reason = "the rounding to single precision is the point"
)]
fn exported(mut mesh: Triangulation, noise: f64) -> Triangulation {
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut next = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((state >> 11) as f64 / (1u64 << 53) as f64).mul_add(2.0, -1.0) * noise
    };
    for p in &mut mesh.positions {
        let moved = Point::new(p.x + next(), p.y + next(), p.z + next());
        *p = Point::new(
            f64::from(moved.x as f32),
            f64::from(moved.y as f32),
            f64::from(moved.z as f32),
        );
    }
    mesh
}

/// A rounded slab whose rims are chamfered by a single row of triangles:
/// each row's corners lie on two circles, which a sphere and a torus fit as
/// well as the cone. The cone is the one taken, and meets the corner's
/// cylinder on one circle. Exported points, up to a little over the default
/// distance off their surfaces, still come back the same faces, on axes
/// square to a cap.
#[test]
fn one_row_chamfers_come_back_cones() {
    for (r, n, noise) in [
        (8.0, 20, 0.0),
        (3.0, 12, 0.0),
        (8.0, 40, 0.0),
        (8.0, 20, 7.5e-5),
    ] {
        let at = format!("corner {r} in {n} steps, noise {noise}");
        let mesh = exported(rounded_slab(r, n, 0.5, 5.0), noise);
        let mut model = Model::new();
        let built = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T)
            .unwrap_or_else(|e| panic!("{at}: {e}"));
        let diagnosis = check(&model, &built.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{at}: {diagnosis}");
        assert_eq!(kinds(&model, &built.shape), [14, 4, 8, 0, 0], "{at}");
        // Every axis square to a cap: a lean of the fit's slop would meet
        // the cap in an ellipse where the part has a circle. The caps are
        // fitted too, so square to one of them, not to the world's `z`.
        use ogeom::geom::SurfaceGeometry as S;
        let faces = explore_unique(&model, &built.shape, ShapeType::Face).unwrap();
        let surfaces: Vec<&S> = faces
            .iter()
            .map(|f| {
                let data = model.node(f).unwrap().data().as_face().unwrap();
                model.geometry().surface(data.surface).unwrap()
            })
            .collect();
        let caps: Vec<Vector> = surfaces
            .iter()
            .filter_map(|s| match s {
                S::Plane(p) if p.plane().frame().z().vector().z.abs() > 0.99 => {
                    Some(p.plane().frame().z().vector())
                }
                _ => None,
            })
            .collect();
        for s in &surfaces {
            let axis = match s {
                S::Cylinder(c) => c.cylinder().frame().z().vector(),
                S::Cone(c) => c.cone().frame().z().vector(),
                _ => continue,
            };
            let lean = caps
                .iter()
                .map(|n| axis.cross(*n).magnitude())
                .fold(f64::INFINITY, f64::min);
            assert!(lean < 1e-12, "{at}: an axis leans {lean:e} off the caps");
        }
    }
}

/// The edges of `shape` whose bounds pass `pick`.
fn edges_where(model: &Model, shape: &Shape, pick: impl Fn(Point, Point) -> bool) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let b = tight_bounds(model, e, T).unwrap();
            pick(b.low().unwrap(), b.high().unwrap())
        })
        .collect()
}

/// The upward planar face at height `z` lying left of `x`.
fn cap_at(model: &Model, shape: &Shape, z: f64, x: f64) -> Shape {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let (p, n) = ogeom::algo::face_normal(model, f, T).unwrap();
            let b = tight_bounds(model, f, T).unwrap();
            (p.z - z).abs() < 1e-4 && n.z > 0.5 && b.high().unwrap().x < x
        })
        .unwrap()
}

/// A block with corners rounded to 4 and its bottom rim chamfered, carrying
/// a raised end 2 wide whose inner corners are rounded to 0.5 against the
/// big corners' walls.
fn stepped_part(model: &mut Model) -> Shape {
    let block = ogeom::algo::make_box(model, Frame::WORLD, (20.0, 10.0, 5.0), T)
        .unwrap()
        .shape;
    let upright = edges_where(model, &block, |lo, hi| hi.z - lo.z > 4.0);
    let block = ogeom::fillet::fillet_edges(model, &block, &upright, 4.0, T)
        .unwrap()
        .shape;
    let rim = edges_where(model, &block, |lo, hi| hi.z < 1e-9 && lo.z > -1e-9);
    let block = ogeom::fillet::chamfer_edges(model, &block, &rim, 0.5, T)
        .unwrap()
        .shape;
    let top = cap_at(model, &block, 5.0, 100.0);
    let column = ogeom::algo::make_prism(model, &top, Vector::new(0.0, 0.0, 3.0), T)
        .unwrap()
        .shape;
    let end = ogeom::algo::make_box(
        model,
        Frame::new(Point::new(-1.0, -1.0, 4.0), Direction::Z, Direction::X, T).unwrap(),
        (3.0, 12.0, 5.0),
        T,
    )
    .unwrap()
    .shape;
    let raised = ogeom::boolean::common(model, &column, &end, T)
        .unwrap()
        .shape;
    let inner = edges_where(model, &raised, |lo, hi| {
        hi.z - lo.z > 2.0 && (lo.x - 2.0).abs() < 1e-6 && (hi.x - 2.0).abs() < 1e-6
    });
    let raised = ogeom::fillet::fillet_edges(model, &raised, &inner, 0.5, T)
        .unwrap()
        .shape;
    ogeom::boolean::fuse(model, &block, &raised, T)
        .unwrap()
        .shape
}

/// A pad sketched on a converted part's raised end and pushed down through
/// it: its walls run along the part's own walls, round its corners, and out
/// through the chamfered bottom, where the chamfer's rims lie on them. The
/// converted part's surfaces sit a few single-precision roundings off the
/// exact part's, and every boolean with the pad matches the exact part's
/// own, including a pad a micron below the top and one out through the
/// bottom.
/// A rounded block in single precision, converted at a coplanar distance
/// a few times its rounding: facets left between the corner balls lie in
/// slivers whose two seams cross in their plane, and threaded straight they
/// let the whole shape tessellate closed.
#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
#[test]
fn a_fine_conversion_of_a_rounded_block_tessellates_closed() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let rounded = ogeom::fillet::fillet_edges(&mut model, &block, &edges, 2.0, T)
        .unwrap()
        .shape;
    let mut mesh = meshed(&model, &rounded);
    for p in &mut mesh.positions {
        *p = Point::new(
            f64::from(p.x as f32),
            f64::from(p.y as f32),
            f64::from(p.z as f32),
        );
    }
    let options = MeshSolidOptions {
        coplanar_distance: Some(8e-7),
        ..MeshSolidOptions::default()
    };
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, &mesh, &options, T).unwrap();
    assert!(check(&back, &out.shape, T).unwrap().is_valid());
    let drawn = ogeom::mesh::triangulate(&back, &out.shape, Deflection::default(), T).unwrap();
    assert!(drawn.is_closed());
}

#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
#[test]
fn a_pad_through_a_converted_part_matches_the_exact_part() {
    let mut model = Model::new();
    let exact = stepped_part(&mut model);
    let mut mesh =
        ogeom::mesh::triangulate(&model, &exact, Deflection::with_chord(0.05).unwrap(), T).unwrap();
    for p in &mut mesh.positions {
        *p = Point::new(
            f64::from(p.x as f32),
            f64::from(p.y as f32),
            f64::from(p.z as f32),
        );
    }
    let converted = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T)
        .unwrap()
        .shape;
    let fine = |model: &Model, shape: &Shape| {
        volume_properties(model, shape, Deflection::with_chord(0.001).unwrap(), T)
            .unwrap()
            .mass
    };
    for (lift, depth) in [(0.0, 3.0), (0.0, 8.0), (0.0, 9.5), (0.001, 9.5)] {
        let mut volumes = Vec::new();
        for base in [&exact, &converted] {
            let cap = cap_at(&model, base, 8.0, 2.5);
            let lowered = ogeom::algo::transformed(
                &mut model,
                &cap,
                ogeom::math::Transform::translation(Vector::new(0.0, 0.0, -lift)),
            )
            .unwrap()
            .shape;
            let pad =
                ogeom::algo::make_prism(&mut model, &lowered, Vector::new(0.0, 0.0, -depth), T)
                    .unwrap()
                    .shape;
            let mut three = Vec::new();
            for op in 0..3 {
                let at = format!("lift {lift}, depth {depth}, op {op}");
                let built = match op {
                    0 => ogeom::boolean::fuse(&mut model, base, &pad, T),
                    1 => ogeom::boolean::cut(&mut model, base, &pad, T),
                    _ => ogeom::boolean::common(&mut model, base, &pad, T),
                }
                .unwrap_or_else(|e| panic!("{at}: {e}"));
                let diagnosis = check(&model, &built.shape, T).unwrap();
                assert!(diagnosis.is_valid(), "{at}: {diagnosis}");
                three.push(fine(&model, &built.shape));
            }
            volumes.push(three);
        }
        // The conversion's own slop: surfaces a few roundings of a
        // single-precision coordinate off the exact ones, over walls of a
        // few hundred square millimetres.
        for (op, (e, c)) in volumes[0].iter().zip(&volumes[1]).enumerate() {
            assert!(
                (e - c).abs() < 1e-3,
                "lift {lift}, depth {depth}, op {op}: converted {c} against exact {e}"
            );
        }
    }
}

/// A sliver folded back under its three neighbours, as an exporter leaves
/// where it moved a vertex across a thin triangle, closed below by a
/// pyramid. The sliver faces against every neighbour while its winding
/// agrees with theirs. The diagonal it shares with its largest neighbour
/// is swapped, and the solid comes back valid and bounded by the unfolded
/// surface.
#[test]
fn a_folded_sliver_is_unfolded() {
    let at = |x: f64, y: f64, z: f64| {
        Point::new((x - 133.3) * 10.0, (y - 131.0) * 10.0, (z - 5.0) * 10.0)
    };
    let [a, b, c, d, e, f] = [
        at(133.3513, 130.9648, 4.9962),
        at(133.9104, 130.7647, 4.9981),
        at(133.0158, 131.0467, 4.9987),
        at(133.1917, 131.7724, 4.9975),
        at(132.8242, 130.9111, 4.9741),
        at(133.2858, 130.8868, 4.9726),
    ];
    let apex = Point::new(0.0, 0.0, -5.0);
    let mut triangles = vec![[a, b, c], [b, d, c], [e, a, c], [f, b, a]];
    for (u, v) in [(b, d), (d, c), (c, e), (e, a), (a, f), (f, b)] {
        triangles.push([v, u, apex]);
    }
    let mesh = soup(triangles.iter().copied());
    let mut model = Model::new();
    let options = MeshSolidOptions {
        recognize: false,
        ..MeshSolidOptions::default()
    };
    let built = solid_from_mesh(&mut model, &mesh, &options, T).unwrap();
    let diagnosis = check(&model, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    // The unfolded surface: the sliver and its neighbour across B C become
    // A B D and D C A, and the divergence theorem gives the volume.
    triangles[0] = [a, b, d];
    triangles[1] = [d, c, a];
    let want: f64 = triangles
        .iter()
        .map(|[p, q, r]| {
            let (p, q, r) = (*p - Point::ORIGIN, *q - Point::ORIGIN, *r - Point::ORIGIN);
            p.dot(q.cross(r)) / 6.0
        })
        .sum();
    let v = volume(&model, &built.shape);
    assert!((v - want).abs() < want * 1e-9, "{v} against {want}");
}

/// A tube with windows cut through its wall, meshed and moved by an
/// exporter's noise: the fits refuse a few vertices, and dropping the
/// triangles that bring them can cut one wall's region into patches joined
/// by nothing. Each patch is a face of its own, and the solid holds the
/// exact part's volume.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    reason = "the rounding to single precision is the point"
)]
#[test]
fn a_windowed_tube_s_wall_patches_each_become_a_face() {
    let mut model = Model::new();
    let outer = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 16.0, 30.0, T)
        .unwrap()
        .shape;
    let bore = Frame::new(Point::new(0.0, 0.0, 3.0), Direction::Z, Direction::X, T).unwrap();
    let inner = ogeom::algo::make_cylinder(&mut model, bore, 14.0, 30.0, T)
        .unwrap()
        .shape;
    let mut part = ogeom::boolean::cut(&mut model, &outer, &inner, T)
        .unwrap()
        .shape;
    for (angle, z0, z1, width) in [
        (0.3_f64, 6.0, 20.0, 6.0),
        (1.5, 6.0, 26.0, 5.0),
        (2.9, 22.0, 26.0, 8.0),
        (4.2, 8.0, 14.0, 4.0),
    ] {
        let x = Direction::new(Vector::new(angle.cos(), angle.sin(), 0.0), T).unwrap();
        let at = Frame::new(Point::new(0.0, 0.0, z0), Direction::Z, x, T).unwrap();
        let local = Frame::new(
            Point::new(10.0, -width / 2.0, 0.0),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        let window = ogeom::algo::make_box(&mut model, local, (10.0, width, z1 - z0), T)
            .unwrap()
            .shape;
        let window =
            ogeom::algo::transformed(&mut model, &window, ogeom::math::Transform::from_frame(&at))
                .unwrap()
                .shape;
        part = ogeom::boolean::cut(&mut model, &part, &window, T)
            .unwrap()
            .shape;
    }
    let exact = volume_properties(&model, &part, Deflection::with_chord(0.001).unwrap(), T)
        .unwrap()
        .mass;
    let mut mesh =
        ogeom::mesh::triangulate(&model, &part, Deflection::with_chord(0.05).unwrap(), T).unwrap();
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut noise = || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((state >> 11) as f64 / (1u64 << 53) as f64).mul_add(2.0, -1.0) * 1e-4
    };
    for p in &mut mesh.positions {
        let moved = Point::new(p.x + 128.0 + noise(), p.y + 128.0 + noise(), p.z + noise());
        *p = Point::new(
            f64::from(moved.x as f32) - 128.0,
            f64::from(moved.y as f32) - 128.0,
            f64::from(moved.z as f32),
        );
    }
    let built = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(!built.report.recognition_withdrawn);
    assert!(built.report.curved_faces >= 2);
    let diagnosis = check(&model, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let v = volume_properties(
        &model,
        &built.shape,
        Deflection::with_chord(0.001).unwrap(),
        T,
    )
    .unwrap()
    .mass;
    assert!((v - exact).abs() < exact * 1e-3, "{v} against {exact}");
}

/// A rounded block a hundred millimetres from the origin, written as ASCII
/// STL to six significant digits: every coordinate is rounded to the
/// thousandth, a hundred times the millionth of the diagonal the converter
/// allows by default. Read with its quantum, the rounding is the file's and
/// not the surfaces', and the fillets come back cylinders.
#[test]
fn an_ascii_file_s_rounding_is_allowed_its_vertices() {
    use std::fmt::Write as _;
    let mut model = Model::new();
    let at = Frame::new(Point::new(100.0, 100.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let block = ogeom::algo::make_box(&mut model, at, (40.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let upright = edges_where(&model, &block, |lo, hi| hi.z - lo.z > 5.0);
    let rounded = ogeom::fillet::fillet_edges(&mut model, &block, &upright, 3.0, T)
        .unwrap()
        .shape;
    let mesh = ogeom::mesh::triangulate(&model, &rounded, Deflection::with_chord(0.01).unwrap(), T)
        .unwrap();
    let mut text = String::from("solid block\n");
    for t in &mesh.triangles {
        text += "facet normal 0 0 0\nouter loop\n";
        for &i in t {
            let p = mesh.positions[i as usize];
            let _ = writeln!(text, "vertex {:.5e} {:.5e} {:.5e}", p.x, p.y, p.z);
        }
        text += "endloop\nendfacet\n";
    }
    text += "endsolid block\n";
    let read = ogeom::io::stl::read_with_quantum(text.as_bytes(), T).unwrap();
    let options = MeshSolidOptions {
        quantum: Some(read.quantum),
        ..MeshSolidOptions::default()
    };
    let mut back = Model::new();
    let built = solid_from_mesh(&mut back, &read.mesh, &options, T).unwrap();
    assert_eq!(kinds(&back, &built.shape), [6, 4, 0, 0, 0]);
    assert!(check(&back, &built.shape, T).unwrap().is_valid());
}

/// Two blocks touching along one edge, written as one mesh: four triangles
/// use that edge, two from each block. The blocks are what the edges used
/// exactly twice join, and at the shared edge each block's own two
/// triangles pair up: both close, and come back two solids.
#[test]
fn blocks_glued_along_an_edge_come_back_two_solids() {
    let corner = |base: Point, i: u32, size: f64| {
        Point::new(
            base.x + f64::from(i & 1) * size,
            base.y + f64::from((i >> 1) & 1) * size,
            base.z + f64::from((i >> 2) & 1) * size,
        )
    };
    let faces: [[u32; 3]; 12] = [
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
    // The second block starts where the first ends in x and y, so the
    // first's edge x = 10, y = 10 is the second's edge x = 10, y = 10.
    let triangles = [Point::ORIGIN, Point::new(10.0, 10.0, 0.0)]
        .into_iter()
        .flat_map(|base| faces.iter().map(move |f| f.map(|i| corner(base, i, 10.0))));
    let mesh = soup(triangles);
    let mut model = Model::new();
    let built = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(built.closed);
    assert_eq!(built.report.edges_used_more, 1);
    assert_eq!(count(&model, &built.shape, ShapeType::Solid), 2);
    assert!(check(&model, &built.shape, T).unwrap().is_valid());
    let v = volume(&model, &built.shape);
    assert!((v - 2000.0).abs() < 1e-9, "{v}");
}

/// A thin disc whose top rim is rounded, meshed coarsely: the rim's round
/// meets the disc's flat top along a tangent circle, and the top's facets
/// across that circle have their corners on the round. They are the top's,
/// peeled from the round's region, and the disc comes back a sound solid on
/// its exact surfaces.
#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
#[test]
fn a_rounded_rim_leaves_the_flat_top_its_facets() {
    let mut model = Model::new();
    let disc = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 12.7, 0.4, T)
        .unwrap()
        .shape;
    let rim = edges_where(&model, &disc, |lo, _| lo.z > 0.4 - 1e-6);
    let rounded = ogeom::fillet::fillet_edges(&mut model, &disc, &rim, 0.2, T)
        .unwrap()
        .shape;
    let mut mesh =
        ogeom::mesh::triangulate(&model, &rounded, Deflection::with_chord(0.1).unwrap(), T)
            .unwrap();
    for p in &mut mesh.positions {
        *p = Point::new(
            f64::from(p.x as f32),
            f64::from(p.y as f32),
            f64::from(p.z as f32),
        );
    }
    let exact = volume(&model, &rounded);
    let built = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(built.closed);
    let diagnosis = check(&model, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    assert_eq!(
        kinds(&model, &built.shape)[0],
        2,
        "the top and bottom stay planes"
    );
    let v = volume(&model, &built.shape);
    assert!((v - exact).abs() < exact * 1e-5, "{v} against {exact}");
}

/// A plate with rounded corners and a rounded bottom rim, meshed finely: the
/// rim's rounds meet the plate's flat walls and bottom along tangents, and
/// the walls' long facets lie on a round's surface at their near corners.
/// A flat facet much larger than a round's own comes with its whole flat
/// face or not at all, and a flat face's facets across a tangent circle are
/// peeled from a round's region: the rounds stop where the walls begin, and
/// the plate comes back on its eighteen exact surfaces where it took a
/// thousand faces.
#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
#[test]
fn a_round_does_not_take_a_flat_face_s_facets() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (47.0, 21.5, 5.0), T)
        .unwrap()
        .shape;
    let upright = edges_where(&model, &block, |lo, hi| hi.z - lo.z > 4.0);
    let plate = ogeom::fillet::fillet_edges(&mut model, &block, &upright, 8.0, T)
        .unwrap()
        .shape;
    let rim = edges_where(&model, &plate, |lo, hi| hi.z < 1e-9 && lo.z > -1e-9);
    let rounded = ogeom::fillet::fillet_edges(&mut model, &plate, &rim, 0.5, T)
        .unwrap()
        .shape;
    let exact = volume_properties(&model, &rounded, Deflection::with_chord(0.001).unwrap(), T)
        .unwrap()
        .mass;
    let mut mesh =
        ogeom::mesh::triangulate(&model, &rounded, Deflection::with_chord(0.005).unwrap(), T)
            .unwrap();
    for p in &mut mesh.positions {
        *p = Point::new(
            f64::from(p.x as f32),
            f64::from(p.y as f32),
            f64::from(p.z as f32),
        );
    }
    let built = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(check(&model, &built.shape, T).unwrap().is_valid());
    // Top, bottom, four walls, four corners, four rounds along the walls
    // and four at the corners.
    assert_eq!(kinds(&model, &built.shape), [6, 8, 0, 0, 4]);
    let v = volume_properties(
        &model,
        &built.shape,
        Deflection::with_chord(0.001).unwrap(),
        T,
    )
    .unwrap()
    .mass;
    assert!((v - exact).abs() < exact * 1e-4, "{v} against {exact}");
}

/// A hemisphere on a cylinder of its radius, meshed finely and coarsely:
/// the cylinder's seam and the cap's run from the circle they share, and
/// the sphere takes the cylinder's frame so both start at one angle. Three
/// faces, sound, holding the exact volume.
#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
#[test]
fn a_dome_on_a_cylinder_shares_its_seam() {
    for chord in [0.01, 0.05] {
        let mut model = Model::new();
        let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 10.0, 10.0, T)
            .unwrap()
            .shape;
        let top = Frame::new(Point::new(0.0, 0.0, 10.0), Direction::Z, Direction::X, T).unwrap();
        let ball = ogeom::algo::make_sphere(&mut model, top, 10.0, T)
            .unwrap()
            .shape;
        let dome = ogeom::boolean::fuse(&mut model, &drum, &ball, T)
            .unwrap()
            .shape;
        let pi = core::f64::consts::PI;
        let exact = pi * 100.0 * 10.0 + 2.0 / 3.0 * pi * 1000.0;
        let mut mesh =
            ogeom::mesh::triangulate(&model, &dome, Deflection::with_chord(chord).unwrap(), T)
                .unwrap();
        for p in &mut mesh.positions {
            *p = Point::new(
                f64::from(p.x as f32),
                f64::from(p.y as f32),
                f64::from(p.z as f32),
            );
        }
        let built = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
        let at = format!("chord {chord}");
        let diagnosis = check(&model, &built.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{at}: {diagnosis}");
        assert_eq!(kinds(&model, &built.shape), [1, 1, 0, 1, 0], "{at}");
        let v = volume(&model, &built.shape);
        assert!(
            (v - exact).abs() < exact * 1e-5,
            "{at}: {v} against {exact}"
        );
    }
}

/// A drilled block exported with noise a few times the default distance:
/// the coplanar distance follows how far the flat faces' own vertices stand
/// off them, and the block comes back its seven faces.
#[test]
fn a_noisy_drilled_block_keeps_its_faces_whole() {
    let mut model = Model::new();
    let drilled = drilled_block(&mut model);
    let mesh = ogeom::mesh::triangulate(&model, &drilled, Deflection::with_chord(0.05).unwrap(), T)
        .unwrap();
    let mesh = exported(mesh, 5e-5);
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert_eq!(
        kinds(&back, &out.shape),
        [6, 1, 0, 0, 0],
        "{:?}",
        out.report
    );
    holds((&model, &drilled), (&back, &out.shape));
}

/// A round two facets across fits many surfaces through its three rows of
/// vertices. The one it comes back on is the cylinder tangent to the two
/// faces it joins, on their line of meeting's direction, at the round's
/// radius.
#[test]
fn a_narrow_round_comes_back_tangent_to_its_faces() {
    use ogeom::geom::SurfaceGeometry as S;
    let mut model = Model::new();
    let rounded = narrow_round(&mut model);
    let mesh = ogeom::mesh::triangulate(&model, &rounded, Deflection::with_chord(0.05).unwrap(), T)
        .unwrap();
    let mesh = exported(mesh, 0.0);
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert_eq!(
        kinds(&back, &out.shape),
        [6, 1, 0, 0, 0],
        "{:?}",
        out.report
    );
    let cylinder = explore_unique(&back, &out.shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find_map(|f| {
            let data = back.node(&f).unwrap().data().as_face().unwrap();
            match back.geometry().surface(data.surface).unwrap() {
                S::Cylinder(c) => Some(c.cylinder()),
                _ => None,
            }
        })
        .unwrap();
    let (o, z) = (cylinder.frame().origin(), cylinder.frame().z().vector());
    assert!(z.cross(Vector::Y).magnitude() < 1e-12, "axis {z:?}");
    let r = cylinder.radius();
    assert!((r - 0.5).abs() < 1e-5, "radius {r}");
    // Tangent: the axis stands the radius in from each face.
    assert!((20.0 - o.x - r).abs() < 1e-12, "{o:?} against radius {r}");
    assert!((5.0 - o.z - r).abs() < 1e-12, "{o:?} against radius {r}");
    holds((&model, &rounded), (&back, &out.shape));
}

/// A disc 0.8 thick with its top rim rounded 0.3.
fn rounded_disc(model: &mut Model) -> Shape {
    let disc = ogeom::algo::make_cylinder(model, Frame::WORLD, 12.7, 0.8, T)
        .unwrap()
        .shape;
    let top = edges_where(model, &disc, |lo, _| lo.z > 0.8 - 1e-6);
    ogeom::fillet::fillet_edges(model, &disc, &top, 0.3, T)
        .unwrap()
        .shape
}

/// A thin disc's rounded rim, four rows of facets round a tube much
/// narrower than the rim: the round takes the flat top's facets in with it
/// (their corners are all on the rim circle), and without them its top
/// row's normals lean to one side and the fit finds nothing. The torus
/// already found still holds, and the rim comes back one face.
#[test]
fn a_thin_disc_s_rounded_rim_comes_back_a_torus() {
    let mut model = Model::new();
    let rounded = rounded_disc(&mut model);
    let mesh = meshed(&model, &rounded);
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert_eq!(
        kinds(&back, &out.shape),
        [2, 1, 0, 0, 1],
        "{:?}",
        out.report
    );
    holds((&model, &rounded), (&back, &out.shape));
}

/// The same disc in single precision: the torus fitted to the rounded
/// points has its tube's crest a few slops off the flat top, and the top's
/// rim, level at the crest, is still the crest's circle, which bounds the
/// torus as a band between two circles.
#[test]
fn a_rounded_rim_s_crest_is_a_circle() {
    let mut model = Model::new();
    let rounded = rounded_disc(&mut model);
    let mesh = exported(meshed(&model, &rounded), 0.0);
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert_eq!(
        kinds(&back, &out.shape),
        [2, 1, 0, 0, 1],
        "{:?}",
        out.report
    );
    holds((&model, &rounded), (&back, &out.shape));
}

/// A block with one edge rounded 0.5, meshed two facets across the round.
fn narrow_round(model: &mut Model) -> Shape {
    let block = ogeom::algo::make_box(model, Frame::WORLD, (20.0, 10.0, 5.0), T)
        .unwrap()
        .shape;
    let edge = edges_where(model, &block, |lo, hi| {
        lo.x > 20.0 - 1e-6 && lo.z > 5.0 - 1e-6 && hi.y - lo.y > 9.0
    });
    ogeom::fillet::fillet_edges(model, &block, &edge, 0.5, T)
        .unwrap()
        .shape
}

/// Converted without recognition, keeping every vertex, a block rounded on
/// every edge keeps its rounds' facets. Refined, the facets are the
/// cylinders and the corners' spheres again, as converting with
/// recognition gives in one step.
#[test]
fn a_verbatim_conversion_refines_to_its_surfaces() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (30.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let rounded = ogeom::fillet::fillet_edges(&mut model, &block, &edges, 2.0, T)
        .unwrap()
        .shape;
    let mesh = ogeom::mesh::triangulate(&model, &rounded, Deflection::with_chord(0.05).unwrap(), T)
        .unwrap();
    let mut back = Model::new();
    let verbatim = MeshSolidOptions {
        recognize: false,
        keep_vertices: true,
        ..MeshSolidOptions::default()
    };
    let first = solid_from_mesh(&mut back, &mesh, &verbatim, T).unwrap();
    assert_eq!(kinds(&back, &first.shape)[1], 0, "{:?}", first.report);
    let refined =
        ogeom::algo::refine_solid(&mut back, &first.shape, &MeshSolidOptions::default(), T)
            .unwrap();
    assert_eq!(
        kinds(&back, &refined.shape),
        [6, 12, 0, 8, 0],
        "{:?}",
        refined.report
    );
    holds((&model, &rounded), (&back, &refined.shape));
}

/// A part read from the STL file named by `OGEOM_TEST_77777` (a file that is
/// not ours to bundle), converted with the default options; `None`, and the
/// test skipped, where the variable is unset.
fn converted_part(model: &mut Model) -> Option<Shape> {
    let path = std::env::var_os("OGEOM_TEST_77777")?;
    let bytes = std::fs::read(path).expect("the file the variable names reads");
    let mesh = ogeom::io::stl::read(&bytes, T).unwrap();
    let out = solid_from_mesh(model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(out.closed, "{:?}", out.report);
    Some(out.shape)
}

/// The upward flat face at the top of `shape` whose bounds hold `(x, y)`.
fn top_face_over(model: &Model, shape: &Shape, x: f64, y: f64) -> Shape {
    let top = tight_bounds(model, shape, T).unwrap().high().unwrap().z;
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let (p, n) = ogeom::algo::face_normal(model, f, T).unwrap();
            let b = tight_bounds(model, f, T).unwrap();
            let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
            n.z > 1.0 - 1e-9
                && (p.z - top).abs() < 1e-6
                && (lo.x..=hi.x).contains(&x)
                && (lo.y..=hi.y).contains(&y)
        })
        .expect("a top face over the point")
}

/// A flat face's outline as a sketch holding its points in single precision
/// draws it: each line between its ends, each arc through its ends and its
/// middle, every point rounded to `f32`. The arcs come back on other
/// centres and radii, the sides swept from them a few microns off the
/// walls they were drawn from.
#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
fn outline_in_single_precision(model: &mut Model, face: &Shape) -> Shape {
    use ogeom::geom::{Curve, Curve3d as _, LineCurve};
    let round = |p: Point| {
        Point::new(
            f64::from(p.x as f32),
            f64::from(p.y as f32),
            f64::from(p.z as f32),
        )
    };
    // Each edge's ends (its vertices) and middle, and whether it is an arc.
    let vertex = |model: &Model, v: &Shape| {
        (
            v.node().index(),
            model.node(v).unwrap().data().as_vertex().unwrap().point,
        )
    };
    let mut spans: Vec<([u32; 2], [Point; 3], bool)> = Vec::new();
    for edge in explore_unique(model, face, ShapeType::Edge).unwrap() {
        let (first, last) = ogeom::algo::edge_vertices(model, &edge).unwrap().unwrap();
        let ((ia, a), (ib, b)) = (vertex(model, &first), vertex(model, &last));
        let data = model.node(&edge).unwrap().data().as_edge().unwrap();
        let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            panic!("an edge with no curve");
        };
        let curve = model.geometry().curve(*curve).unwrap();
        let middle = curve.point_at(f64::midpoint(range.0, range.1), T).unwrap();
        spans.push((
            [ia, ib],
            [a, middle, b].map(round),
            matches!(curve, Curve::Circle(_)),
        ));
    }
    // Chained end to start.
    let first = spans.remove(0);
    let mut end = first.0[1];
    let mut chain = vec![(first.1, first.2)];
    while !spans.is_empty() {
        let next = spans
            .iter()
            .position(|(ends, _, _)| ends.contains(&end))
            .expect("the outline closes");
        let (ends, mut points, arc) = spans.remove(next);
        if ends[1] == end {
            points.reverse();
            end = ends[0];
        } else {
            end = ends[1];
        }
        chain.push((points, arc));
    }
    // Anticlockwise seen from above, as the face's own outline runs.
    let area: f64 = chain
        .iter()
        .map(|(p, _)| p[0].x * p[2].y - p[2].x * p[0].y)
        .sum();
    if area < 0.0 {
        chain.reverse();
        for (points, _) in &mut chain {
            points.reverse();
        }
    }
    let corners: Vec<Shape> = chain
        .iter()
        .map(|(p, _)| model.add_vertex(ogeom::topo::VertexData::new(p[0])))
        .collect();
    let mut edges = Vec::new();
    for (i, (p, arc)) in chain.iter().enumerate() {
        let (from, to) = (&corners[i], &corners[(i + 1) % corners.len()]);
        let edge = if *arc {
            let circle = ogeom::math::Circle::through(p[0], p[1], p[2], T).unwrap();
            let param =
                |q: Point| ogeom::math::elementary::circle_parameter(&circle, q, T).unwrap();
            let tau = core::f64::consts::TAU;
            let (a, b, c) = (param(p[0]), param(p[1]), param(p[2]));
            let ahead = |t: f64| (t - a).rem_euclid(tau);
            let curve: Curve = ogeom::geom::CircleCurve::new(circle).into();
            if ahead(b) < ahead(c) {
                ogeom::algo::make_edge_between(model, curve, (a, a + ahead(c)), from, to, T)
                    .unwrap()
                    .shape
            } else {
                let back = (c - a).rem_euclid(tau) - tau;
                ogeom::algo::make_edge_between(model, curve, (a + back, a), to, from, T)
                    .unwrap()
                    .shape
                    .reversed()
            }
        } else {
            let line: Curve = LineCurve::segment(p[0], p[2], T).unwrap().into();
            ogeom::algo::make_edge_between(model, line, (0.0, p[0].distance(p[2])), from, to, T)
                .unwrap()
                .shape
        };
        edges.push(edge);
    }
    let wire = ogeom::algo::make_wire(model, &edges, T).unwrap().shape;
    let (at, _) = ogeom::algo::face_normal(model, face, T).unwrap();
    let plane = ogeom::math::Plane::new(
        Frame::new(Point::new(0.0, 0.0, at.z), Direction::Z, Direction::X, T).unwrap(),
    );
    let surface = ogeom::geom::PlaneSurface::new(plane).into();
    ogeom::algo::make_face(model, surface, &[wire], T)
        .unwrap()
        .shape
}

/// A converted part and a prism of one of its own top faces pushed `depth`
/// straight down: the prism's sides lie on the part's walls (a flat one and
/// recognized rounds) over their height, exactly or, where `rounded`, as a
/// single-precision sketch of the face's outline leaves them. Fuse, cut and
/// common each give a sound solid, and the fuse and the common hold what
/// the two inputs do.
fn pad_along_a_converted_part_s_walls(depth: f64, rounded: bool) {
    let mut model = Model::new();
    let Some(part) = converted_part(&mut model) else {
        return;
    };
    let face = top_face_over(&model, &part, 105.157, 89.76);
    let face = if rounded {
        outline_in_single_precision(&mut model, &face)
    } else {
        face
    };
    let pad = ogeom::algo::make_prism(&mut model, &face, Vector::new(0.0, 0.0, -depth), T)
        .unwrap()
        .shape;
    let fine = Deflection::with_chord(1e-3).unwrap();
    let volume =
        |model: &Model, shape: &Shape| volume_properties(model, shape, fine, T).unwrap().mass;
    let (v_part, v_pad) = (volume(&model, &part), volume(&model, &pad));
    let mut results = Vec::new();
    for (name, op) in [
        (
            "fuse",
            ogeom::boolean::fuse as fn(&mut Model, &Shape, &Shape, Tolerances) -> _,
        ),
        ("cut", ogeom::boolean::cut),
        ("common", ogeom::boolean::common),
    ] {
        let out = op(&mut model, &part, &pad, T).unwrap_or_else(|e| panic!("{name}: {e}"));
        let diagnosis = check(&model, &out.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
        results.push(volume(&model, &out.shape));
    }
    let [fused, cut, common] = results[..] else {
        unreachable!()
    };
    let scale = v_part + v_pad;
    assert!(
        (fused + common - scale).abs() < 1e-6 * scale,
        "fuse {fused} + common {common} against part {v_part} + pad {v_pad}"
    );
    assert!(
        (cut + common - v_part).abs() < 1e-6 * scale,
        "cut {cut} + common {common} against part {v_part}"
    );
    // Wholly inside, the pad adds nothing and the common is the pad.
    if depth < 5.0 {
        assert!(
            (fused - v_part).abs() < 1e-6 * scale,
            "{fused} against {v_part}"
        );
        assert!(
            (common - v_pad).abs() < 1e-6 * scale,
            "{common} against {v_pad}"
        );
    } else {
        assert!(fused > v_part + 1.0, "{fused} against {v_part}");
    }
}

/// The pad 3 mm deep, within the part.
#[test]
fn a_pad_along_a_converted_part_s_walls_within_it() {
    pad_along_a_converted_part_s_walls(3.0, false);
}

/// The pad 10 mm deep, out through the part's bottom.
#[test]
fn a_pad_along_a_converted_part_s_walls_out_through_its_bottom() {
    pad_along_a_converted_part_s_walls(10.0, false);
}

/// The pad 3 mm deep, from a single-precision sketch of the face.
#[test]
fn a_sketched_pad_along_a_converted_part_s_walls_within_it() {
    pad_along_a_converted_part_s_walls(3.0, true);
}

/// The pad 10 mm deep, from a single-precision sketch of the face, out
/// through the part's bottom.
#[test]
fn a_sketched_pad_along_a_converted_part_s_walls_out_through_its_bottom() {
    pad_along_a_converted_part_s_walls(10.0, true);
}

/// A whole sphere or torus meshed in single precision far from the origin
/// comes back as one face whose looser tolerance its seam and poles share,
/// as the checker requires of everything a face bounds.
#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
#[test]
fn a_whole_sphere_or_torus_from_single_precision_checks_valid() {
    let mut source = Model::new();
    let ball = ogeom::algo::make_sphere(&mut source, Frame::WORLD, 70.0, T)
        .unwrap()
        .shape;
    let ring = ogeom::algo::make_torus(&mut source, Frame::WORLD, 100.0, 30.0, T)
        .unwrap()
        .shape;
    for whole in [ball, ring] {
        let mut mesh =
            ogeom::mesh::triangulate(&source, &whole, Deflection::with_chord(0.1).unwrap(), T)
                .unwrap();
        for p in &mut mesh.positions {
            *p = Point::new(
                f64::from((p.x + 100.0) as f32),
                f64::from((p.y + 200.0) as f32),
                f64::from((p.z + 300.0) as f32),
            );
        }
        let options = MeshSolidOptions {
            quantum: Some(ogeom::algo::single_precision_quantum(&mesh)),
            ..MeshSolidOptions::default()
        };
        let mut model = Model::new();
        let built = solid_from_mesh(&mut model, &mesh, &options, T).unwrap();
        assert_eq!(count(&model, &built.shape, ShapeType::Face), 1);
        let diagnosis = check(&model, &built.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{diagnosis}");
    }
}

/// Options that say nothing a distance or a turn can mean are refused by
/// name rather than read as something else.
#[test]
fn unusable_converter_options_are_refused() {
    let cube = cube_soup(2.0);
    let bad = [
        MeshSolidOptions {
            coplanar_distance: Some(f64::NAN),
            ..MeshSolidOptions::default()
        },
        MeshSolidOptions {
            coplanar_distance: Some(0.0),
            ..MeshSolidOptions::default()
        },
        MeshSolidOptions {
            coplanar_distance: Some(-1e-3),
            ..MeshSolidOptions::default()
        },
        MeshSolidOptions {
            crease: f64::NAN,
            ..MeshSolidOptions::default()
        },
        MeshSolidOptions {
            crease: 0.0,
            ..MeshSolidOptions::default()
        },
        MeshSolidOptions {
            quantum: Some(-1.0),
            ..MeshSolidOptions::default()
        },
    ];
    for options in bad {
        let mut model = Model::new();
        assert!(
            solid_from_mesh(&mut model, &cube, &options, T).is_err(),
            "{options:?}"
        );
    }
    let mut model = Model::new();
    assert!(solid_from_mesh(&mut model, &cube, &MeshSolidOptions::default(), T).is_ok());
}

/// A small cube inside a large one is the large one's void, wherever a ray
/// from it happens to meet the large one: here along the diagonal two of
/// its triangles share, where a single ray counts the crossing twice.
#[test]
fn a_void_whose_ray_meets_a_shared_edge_stays_a_void() {
    let shifted = |mesh: Triangulation, by: Vector| {
        let mut mesh = mesh;
        for p in &mut mesh.positions {
            *p += by;
        }
        mesh
    };
    // A ray from the small cube's first corner along the probe direction
    // meets the large cube's +x face on its diagonal.
    let probe = Vector::new(0.577_215_664_9, 0.618_033_988_7, 0.533_751_168_7);
    let (reach, hit) = (21.1, Point::new(30.0, 15.31, 14.69));
    let corner = hit - probe * reach;
    let outer = cube_soup(30.0);
    let inner = shifted(cube_soup(5.0), corner - Point::ORIGIN);
    let mut mesh = outer;
    let base = u32::try_from(mesh.positions.len()).unwrap();
    mesh.positions.extend(inner.positions);
    mesh.triangles
        .extend(inner.triangles.iter().map(|t| t.map(|i| i + base)));
    let mut model = Model::new();
    let built = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert_eq!(count(&model, &built.shape, ShapeType::Solid), 1);
    let v = volume(&model, &built.shape);
    assert!((v - (27000.0 - 125.0)).abs() < 1e-6, "{v}");
}

/// A drum meshed finely and stored in single precision a thousand
/// millimetres out keeps its flat caps whole: their small triangles lean by
/// the rounding of their corners, which is no turn of the cap.
#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
#[test]
fn caps_of_a_fine_mesh_far_out_stay_whole() {
    let mut source = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut source, Frame::WORLD, 5.0, 12.0, T)
        .unwrap()
        .shape;
    let mut mesh =
        ogeom::mesh::triangulate(&source, &drum, Deflection::with_chord(0.01).unwrap(), T).unwrap();
    let (sin, cos) = 0.3_f64.sin_cos();
    for p in &mut mesh.positions {
        let (y, z) = (p.y * cos - p.z * sin, p.y * sin + p.z * cos);
        let (x, y) = (p.x * cos - y * sin, p.x * sin + y * cos);
        *p = Point::new(
            f64::from((x + 1000.0) as f32),
            f64::from((y + 2000.0) as f32),
            f64::from((z + 3000.0) as f32),
        );
    }
    let options = MeshSolidOptions {
        quantum: Some(ogeom::algo::single_precision_quantum(&mesh)),
        ..MeshSolidOptions::default()
    };
    let mut model = Model::new();
    let built = solid_from_mesh(&mut model, &mesh, &options, T).unwrap();
    assert_eq!(kinds(&model, &built.shape), [2, 1, 0, 0, 0]);
    assert!(check(&model, &built.shape, T).unwrap().is_valid());
}

/// A slab drafted inward from its top by `draft` times the depth to the
/// power one and a half, one corner rounded in `counts` segments row by
/// row, tessellated and rounded to `f32` a hundred millimetres from the
/// origin, converted face for facet, and a pad of its top face pushed
/// `depth` down into it: fused back and in common with it, each valid, the
/// two volumes adding up to the slab's and the pad's.
fn pad_on_a_drafted_slab(depth: f64, draft: f64, counts: [u32; 7]) {
    let rows = [0.0_f64, 1.5, 3.0, 4.5, 5.8333, 7.1667, 8.5];
    let top = 8.5;
    #[allow(
        clippy::cast_possible_truncation,
        reason = "the rounding to single precision is the point"
    )]
    let single = |v: f64| f64::from(v as f32);
    let mut mesh = Triangulation::default();
    let mut rings: Vec<Vec<u32>> = Vec::new();
    for (&z, &segments) in rows.iter().zip(&counts) {
        // The wall leans in as it falls, faster further down.
        let inset = draft * (top - z).powf(1.5);
        let mut outline = vec![(inset, inset), (10.0 - inset, inset)];
        for k in 0..=segments {
            let a = core::f64::consts::FRAC_PI_2 * f64::from(k) / f64::from(segments);
            outline.push((8.0 + (2.0 - inset) * a.cos(), 8.0 + (2.0 - inset) * a.sin()));
        }
        outline.push((inset, 10.0 - inset));
        let mut ring = Vec::new();
        for (x, y) in outline {
            mesh.positions
                .push(Point::new(single(x + 100.0), single(y + 80.0), single(z)));
            ring.push(u32::try_from(mesh.positions.len() - 1).unwrap());
        }
        rings.push(ring);
    }
    // Rows of different counts are stitched by their place round the
    // outline.
    let place = |ring: &[u32], mesh: &Triangulation| -> Vec<f64> {
        let points: Vec<Point> = ring.iter().map(|&i| mesh.positions[i as usize]).collect();
        let mut along = vec![0.0];
        for k in 0..points.len() {
            let step = points[k].distance(points[(k + 1) % points.len()]);
            along.push(along[k] + step);
        }
        let total = along[points.len()];
        along.iter().map(|a| a / total).collect()
    };
    for j in 0..rows.len() - 1 {
        let (lo, hi) = (&rings[j], &rings[j + 1]);
        let (at_lo, at_hi) = (place(lo, &mesh), place(hi, &mesh));
        let (nl, nh) = (lo.len(), hi.len());
        let (mut a, mut b) = (0_usize, 0_usize);
        while a < nl || b < nh {
            if b >= nh || (a < nl && at_lo[a + 1] <= at_hi[b + 1]) {
                mesh.triangles
                    .push([lo[a % nl], lo[(a + 1) % nl], hi[b % nh]]);
                a += 1;
            } else {
                mesh.triangles
                    .push([lo[a % nl], hi[(b + 1) % nh], hi[b % nh]]);
                b += 1;
            }
        }
    }
    let (bottom, upper) = (&rings[0], &rings[rows.len() - 1]);
    for i in 1..bottom.len() - 1 {
        mesh.triangles.push([bottom[0], bottom[i + 1], bottom[i]]);
    }
    for i in 1..upper.len() - 1 {
        mesh.triangles.push([upper[0], upper[i], upper[i + 1]]);
    }
    let options = MeshSolidOptions {
        recognize: false,
        keep_vertices: true,
        quantum: Some(ogeom::algo::single_precision_quantum(&mesh)),
        ..MeshSolidOptions::default()
    };
    let mut model = Model::new();
    let slab = solid_from_mesh(&mut model, &mesh, &options, T)
        .unwrap()
        .shape;
    let lid = explore_unique(&model, &slab, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let (p, n) = ogeom::algo::face_normal(&model, f, T).unwrap();
            n.z > 0.999 && (p.z - top).abs() < 1e-3
        })
        .unwrap();
    let pad = ogeom::algo::make_prism(&mut model, &lid, Vector::new(0.0, 0.0, -depth), T)
        .unwrap()
        .shape;
    let fuse = ogeom::boolean::fuse(&mut model, &slab, &pad, T)
        .unwrap_or_else(|e| panic!("depth {depth}: fuse: {e}"))
        .shape;
    let common = ogeom::boolean::common(&mut model, &slab, &pad, T)
        .unwrap_or_else(|e| panic!("depth {depth}: common: {e}"))
        .shape;
    for (name, made) in [("fuse", &fuse), ("common", &common)] {
        let diagnosis = check(&model, made, T).unwrap();
        assert!(diagnosis.is_valid(), "depth {depth}: {name}: {diagnosis}");
    }
    let fine = Deflection::with_chord(1e-3).unwrap();
    let v = |s: &Shape| volume_properties(&model, s, fine, T).unwrap().mass;
    let gap = v(&fuse) + v(&common) - v(&slab) - v(&pad);
    assert!(gap.abs() < 1e-4, "depth {depth}: off by {gap}");
}

/// Each wall of the pad stands on the facet just under the top edge they
/// share, the two at a hundredth of a radian, the facet's plane missing
/// the edge's ends by their rounding: the walls cross on that edge, not on
/// a line their planes' solve leaves a sliver under it.
#[test]
fn a_pad_into_a_drafted_single_precision_slab_fuses_back() {
    for depth in [3.0, 10.0] {
        pad_on_a_drafted_slab(depth, 0.02, [6; 7]);
    }
}

/// The same slab drafted a quarter as much, its corner's rows holding
/// eight facets down to five: the rows' diagonal edges run all but in the
/// pad's walls, and the walls' sections on the facets either side of one
/// reach it tens of microns apart, each to within its own doubt. Two such
/// crossings are one point where they stand no further apart than both
/// doubts together.
#[test]
fn a_pad_into_a_slab_whose_rows_differ_fuses_back() {
    for depth in [3.0, 10.0] {
        pad_on_a_drafted_slab(depth, 0.005, [8, 8, 8, 7, 7, 6, 5]);
    }
}

/// Pads from the large top faces of the part `OGEOM_TEST_77777` names,
/// converted face for facet in single precision, fused back into it. Where
/// a pad's wall crosses a facet row's edge that runs all but in the wall's
/// plane, the sections on the facets either side reach the edge microns
/// apart, and they are one junction.
#[test]
fn pads_on_a_part_converted_face_for_facet_fuse_back() {
    let Some(path) = std::env::var_os("OGEOM_TEST_77777") else {
        return;
    };
    let bytes = std::fs::read(path).expect("the file the variable names reads");
    let mesh = ogeom::io::stl::read(&bytes, T).unwrap();
    let options = MeshSolidOptions {
        recognize: false,
        keep_vertices: true,
        quantum: Some(ogeom::algo::single_precision_quantum(&mesh)),
        ..MeshSolidOptions::default()
    };
    let mut model = Model::new();
    let part = solid_from_mesh(&mut model, &mesh, &options, T)
        .unwrap()
        .shape;
    let fine = Deflection::with_chord(1e-3).unwrap();
    let whole = volume_properties(&model, &part, fine, T).unwrap().mass;
    let tops: Vec<Shape> = explore_unique(&model, &part, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            let (_, n) = ogeom::algo::face_normal(&model, f, T).unwrap();
            n.z > 0.999_999
                && ogeom::algo::surface_properties(&model, f, fine, T).is_ok_and(|p| p.mass > 1.0)
        })
        .collect();
    assert!(!tops.is_empty());
    for top in &tops {
        for depth in [3.0, 10.0] {
            let pad = ogeom::algo::make_prism(&mut model, top, Vector::new(0.0, 0.0, -depth), T)
                .unwrap()
                .shape;
            let fused = ogeom::boolean::fuse(&mut model, &part, &pad, T)
                .unwrap_or_else(|e| panic!("depth {depth}: {e}"))
                .shape;
            let diagnosis = check(&model, &fused, T).unwrap();
            assert!(diagnosis.is_valid(), "depth {depth}: {diagnosis}");
            let v = volume_properties(&model, &fused, fine, T).unwrap().mass;
            assert!(
                v > whole - 1e-3 * whole,
                "depth {depth}: {v} against {whole}"
            );
        }
    }
}

/// A mesh converted with the default options, in a model of its own.
fn converted(mesh: &Triangulation) -> (Model, Shape) {
    let mut model = Model::new();
    let out = solid_from_mesh(&mut model, mesh, &MeshSolidOptions::default(), T).unwrap();
    (model, out.shape)
}

/// The rounded block, the rough box and a stepped part in single precision,
/// each converted.
#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
fn converted_parts() -> Vec<(&'static str, Model, Shape)> {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let rounded = ogeom::fillet::fillet_edges(&mut model, &block, &edges, 2.0, T)
        .unwrap()
        .shape;
    let stepped = stepped_part(&mut model);
    let mut single =
        ogeom::mesh::triangulate(&model, &stepped, Deflection::with_chord(0.05).unwrap(), T)
            .unwrap();
    for p in &mut single.positions {
        *p = Point::new(
            f64::from(p.x as f32),
            f64::from(p.y as f32),
            f64::from(p.z as f32),
        );
    }
    let mut out = Vec::new();
    for (name, mesh) in [
        ("rounded block", meshed(&model, &rounded)),
        ("rough box", rough_rounded_box()),
        ("stepped part", single),
    ] {
        let (model, shape) = converted(&mesh);
        out.push((name, model, shape));
    }
    out
}

/// A converted solid written to STEP reads back valid, on the same kinds of
/// surface and with the same volume.
#[test]
fn converted_solids_come_back_through_step() {
    for (name, model, shape) in converted_parts() {
        let before = (kinds(&model, &shape), volume(&model, &shape));
        let mut document = ogeom::doc::Document::over(model);
        document.add_part(name, shape);
        let text = ogeom::io::write_step(&document, T).unwrap();
        let import = ogeom::io::read_step(&text, T).unwrap();
        let back = import.document.model();
        let [solid] = import.solids.as_slice() else {
            panic!("{name}: {} solids came back", import.solids.len());
        };
        let diagnosis = check(back, solid, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
        assert_eq!(kinds(back, solid), before.0, "{name}");
        let after = volume(back, solid);
        assert!(
            (after - before.1).abs() <= before.1 * 1e-6,
            "{name}: {} went out, {after} came back",
            before.1
        );
    }
}

/// The bore of a converted block offset and a side face moved come out as
/// the same edits on the exact block do.
#[test]
fn offsets_and_moves_on_a_converted_block_match_the_exact_block() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 12.0, 8.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(10.0, 6.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drum = ogeom::algo::make_cylinder(&mut model, at, 3.0, 10.0, T)
        .unwrap()
        .shape;
    let exact = ogeom::boolean::cut(&mut model, &block, &drum, T)
        .unwrap()
        .shape;
    let (mut back, conv) = converted(&meshed(&model, &exact));
    assert_eq!(kinds(&back, &conv), [6, 1, 0, 0, 0]);
    let mut changes = Vec::new();
    for (m, shape) in [(&mut model, exact), (&mut back, conv)] {
        let faces = explore_unique(m, &shape, ShapeType::Face).unwrap();
        let bore = faces
            .iter()
            .find(|f| {
                let data = m.node(f).unwrap().data().as_face().unwrap();
                matches!(
                    m.geometry().surface(data.surface),
                    Some(ogeom::geom::SurfaceGeometry::Cylinder(_))
                )
            })
            .cloned()
            .unwrap();
        let side = faces
            .iter()
            .find(|f| ogeom::algo::face_normal(m, f, T).is_ok_and(|(_, n)| n.x < -0.999))
            .cloned()
            .unwrap();
        let before = volume(m, &shape);
        let offset = ogeom::offset::offset_faces(m, &shape, std::slice::from_ref(&bore), 0.5, T)
            .unwrap()
            .shape;
        let step = ogeom::math::Transform::translation(Vector::new(-1.0, 0.0, 0.0));
        let moved = ogeom::offset::move_faces(m, &shape, std::slice::from_ref(&side), &step, T)
            .unwrap()
            .shape;
        for edited in [&offset, &moved] {
            assert!(check(m, edited, T).unwrap().is_valid());
        }
        changes.push([volume(m, &offset) - before, volume(m, &moved) - before]);
    }
    // The bore's radius 3 to 2.5 through 8, and 12 by 8 moved out by 1.
    let want = [core::f64::consts::PI * (9.0 - 6.25) * 8.0, 96.0];
    for change in changes {
        for (got, want) in change.iter().zip(want) {
            assert!((got - want).abs() < want * 1e-6, "{got} against {want}");
        }
    }
}

/// A side face of a block rounded all round, moved out past the fillets it
/// was tangent to, meets them nowhere: the move is refused, on the exact
/// block and on its conversion alike, rather than leaving the face's edges
/// where they were and its surface away from them.
#[test]
fn a_face_moved_clear_of_its_tangent_fillets_is_refused() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let rounded = ogeom::fillet::fillet_edges(&mut model, &block, &edges, 2.0, T)
        .unwrap()
        .shape;
    let (mut back, conv) = converted(&meshed(&model, &rounded));
    let step = ogeom::math::Transform::translation(Vector::new(-0.5, 0.0, 0.0));
    for (m, shape) in [(&mut model, rounded), (&mut back, conv)] {
        let side = explore_unique(m, &shape, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| ogeom::algo::face_normal(m, f, T).is_ok_and(|(_, n)| n.x < -0.999))
            .unwrap();
        let moved = ogeom::offset::move_faces(m, &shape, std::slice::from_ref(&side), &step, T);
        assert!(moved.is_err(), "{:?}", moved.map(|b| volume(m, &b.shape)));
    }
}

/// A spline through a wave: 3 + sin(0.4 t) over t in [0, 10], placed by
/// `at`.
fn wave(at: impl Fn(f64, f64) -> Point) -> ogeom::geom::Curve {
    let points: Vec<Point> = (0..=20)
        .map(|k| {
            let t = f64::from(k) * 0.5;
            at(t, 3.0 + (t * 0.4).sin())
        })
        .collect();
    ogeom::algo::fit::interpolate(&points, 3, ogeom::algo::fit::Spacing::Centripetal, T)
        .unwrap()
        .into()
}

/// A closed wire through `corners`, the edge from the first corner to the
/// second on `curve` and the rest straight.
fn wire_with(model: &mut Model, corners: &[Point], curve: ogeom::geom::Curve) -> Shape {
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|p| ogeom::algo::build::make_vertex(model, *p).shape)
        .collect();
    let mut edges = Vec::new();
    for i in 0..corners.len() {
        let (a, b) = (&vertices[i], &vertices[(i + 1) % corners.len()]);
        let c = if i == 0 {
            curve.clone()
        } else {
            ogeom::geom::LineCurve::segment(corners[i], corners[(i + 1) % corners.len()], T)
                .unwrap()
                .into()
        };
        let range = ogeom::geom::Curve3d::domain(&c);
        edges.push(
            ogeom::algo::build::make_edge_between(model, c, range, a, b, T)
                .unwrap()
                .shape,
        );
    }
    ogeom::algo::make_wire(model, &edges, T).unwrap().shape
}

/// The surface kinds of a shape's faces, by name.
fn surface_kinds(model: &Model, shape: &Shape) -> Vec<ogeom::geom::SurfaceKind> {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .iter()
        .map(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            ogeom::geom::Surface::kind(model.geometry().surface(data.surface).unwrap())
        })
        .collect()
}

/// A mesh of `exact` converted with sweeps comes back valid, closed, of its
/// volume, with one face on `kind`.
fn comes_back_swept(model: &Model, exact: &Shape, kind: ogeom::geom::SurfaceKind, faces: usize) {
    let mesh = meshed(model, exact);
    let options = MeshSolidOptions {
        sweeps: true,
        ..MeshSolidOptions::default()
    };
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, &mesh, &options, T).unwrap();
    let diagnosis = check(&back, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let kinds = surface_kinds(&back, &out.shape);
    assert_eq!(kinds.len(), faces, "{kinds:?}");
    assert_eq!(kinds.iter().filter(|k| **k == kind).count(), 1, "{kinds:?}");
    let fine = Deflection::with_chord(1e-4).unwrap();
    let (want, got) = (
        volume_properties(model, exact, fine, T).unwrap().mass,
        volume_properties(&back, &out.shape, fine, T).unwrap().mass,
    );
    assert!((got - want).abs() < want * 1e-4, "{want} drawn, {got} back");
    let drawn = ogeom::mesh::triangulate(&back, &out.shape, Deflection::default(), T).unwrap();
    assert!(drawn.is_closed());
}

/// A plate whose one side is a wave, pushed up: its wavy wall comes back
/// one extrusion of a fitted profile, not a row of facets.
#[test]
fn a_wavy_wall_comes_back_an_extrusion() {
    let mut model = Model::new();
    let profile = wave(|t, y| Point::new(10.0 - t, y, 0.0));
    let start =
        ogeom::geom::Curve3d::point_at(&profile, ogeom::geom::Curve3d::domain(&profile).0, T)
            .unwrap();
    let end = ogeom::geom::Curve3d::point_at(&profile, ogeom::geom::Curve3d::domain(&profile).1, T)
        .unwrap();
    let wire = wire_with(
        &mut model,
        &[
            start,
            end,
            Point::new(0.0, 0.0, 0.0),
            Point::new(10.0, 0.0, 0.0),
        ],
        profile,
    );
    let plate = ogeom::algo::make_face(
        &mut model,
        ogeom::geom::PlaneSurface::new(ogeom::math::Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let wall = ogeom::algo::make_prism(&mut model, &plate, Vector::new(0.0, 0.0, 5.0), T)
        .unwrap()
        .shape;
    comes_back_swept(&model, &wall, ogeom::geom::SurfaceKind::Extrusion, 6);
}

/// A wavy profile closed down to its axis, turned by `angle` about it.
fn turned_wave(angle: f64) -> (Model, Shape) {
    let mut model = Model::new();
    let profile = wave(|t, r| Point::new(r, 0.0, t));
    let start =
        ogeom::geom::Curve3d::point_at(&profile, ogeom::geom::Curve3d::domain(&profile).0, T)
            .unwrap();
    let end = ogeom::geom::Curve3d::point_at(&profile, ogeom::geom::Curve3d::domain(&profile).1, T)
        .unwrap();
    let wire = wire_with(
        &mut model,
        &[start, end, Point::new(0.0, 0.0, end.z), Point::ORIGIN],
        profile,
    );
    let plane =
        ogeom::math::Plane::new(Frame::new(Point::ORIGIN, Direction::Y, Direction::X, T).unwrap());
    let face = ogeom::algo::make_face(
        &mut model,
        ogeom::geom::PlaneSurface::new(plane).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let turned = ogeom::algo::make_revolution(
        &mut model,
        &face,
        ogeom::math::Axis::new(Point::ORIGIN, Direction::Z),
        angle,
        T,
    )
    .unwrap()
    .shape;
    (model, turned)
}

/// A wavy profile turned half round its axis: its wavy face comes back one
/// surface of revolution, the two profile faces one plane.
#[test]
fn a_turned_wave_comes_back_a_surface_of_revolution() {
    let (model, turned) = turned_wave(core::f64::consts::PI);
    comes_back_swept(&model, &turned, ogeom::geom::SurfaceKind::Revolution, 4);
}

/// A wavy profile turned all the way round: its wavy face comes back one
/// surface of revolution about the axis, closed on itself, between the two
/// flat ends. Over a whole turn a tilt of the axis read from the mesh's
/// normals leaves each parallel's radius changing round it, and the
/// samples fall off any one profile.
#[test]
fn a_wave_turned_all_the_way_round_comes_back_a_surface_of_revolution() {
    let (model, turned) = turned_wave(core::f64::consts::TAU);
    comes_back_swept(&model, &turned, ogeom::geom::SurfaceKind::Revolution, 3);
}
