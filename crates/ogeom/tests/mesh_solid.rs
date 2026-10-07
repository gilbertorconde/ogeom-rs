//! A triangle mesh becomes a B-rep solid: planar faces from its coplanar
//! regions, topology from its own connectivity, windings made to agree.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/pcurves.rs"]
mod pcurves;
#[path = "support/walks.rs"]
mod walks;

use std::time::{Duration, Instant};

use ogeom::algo::{
    MeshSolidOptions, Severity, check, solid_from_mesh, tight_bounds, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, Triangulation, explore_unique};
use walks::edges_walked_one_way;

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

/// A cube whose top is drawn with a needle along one of its edges: a
/// triangle whose third corner stands a micron off the edge, and a fraction
/// of a micron below the top, so its normal leans by a quarter turn and
/// more. It has no plane of its own at the coplanar distance, and comes
/// back part of the top, not a face of almost no area turned whichever way
/// its rounding points.
#[test]
fn a_needle_on_a_flat_face_is_part_of_it() {
    let corner = |i: u32| {
        Point::new(
            f64::from(i & 1) * 10.0,
            f64::from((i >> 1) & 1) * 10.0,
            f64::from((i >> 2) & 1) * 10.0,
        )
    };
    let needle = Point::new(5.0, 1e-6, 10.0 - 5e-7);
    let mut triangles: Vec<[Point; 3]> = [
        [0, 2, 1],
        [1, 2, 3],
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
    let [c4, c5, c6, c7] = [4, 5, 6, 7].map(corner);
    triangles.extend([
        [c4, c5, needle],
        [needle, c5, c7],
        [needle, c7, c6],
        [needle, c6, c4],
    ]);
    let mesh = soup(triangles.into_iter());
    let mut model = Model::new();
    let out = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(check(&model, &out.shape, T).unwrap().is_valid());
    assert_eq!(out.report.faces, 6);
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

/// The kinds of surface a shape's faces are built on, counted; it has no
/// fitted patch.
fn kinds(model: &Model, shape: &Shape) -> [usize; 5] {
    let (out, patches) = kinds_and_patches(model, shape);
    assert_eq!(patches, 0, "a fitted patch among {out:?}");
    out
}

/// The kinds of canonical surface a shape's faces are built on, counted,
/// and its fitted B-spline patches.
fn kinds_and_patches(model: &Model, shape: &Shape) -> ([usize; 5], usize) {
    use ogeom::geom::SurfaceGeometry as S;
    let (mut out, mut patches) = ([0; 5], 0);
    for face in explore_unique(model, shape, ShapeType::Face).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        out[match model.geometry().surface(data.surface).unwrap() {
            S::Plane(_) => 0,
            S::Cylinder(_) => 1,
            S::Cone(_) => 2,
            S::Sphere(_) => 3,
            S::Torus(_) => 4,
            S::BSpline(_) => {
                patches += 1;
                continue;
            }
            _ => panic!("a surface recognition does not build"),
        }] += 1;
    }
    (out, patches)
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
#[ignore = "heavy"]
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
#[ignore = "heavy"]
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
#[ignore = "heavy"]
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
#[ignore = "heavy"]
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
/// on a ball, a ball bored twice across, a ring pierced through its tube.
/// Each comes back as the whole surface with the holes as inner wires, its
/// seams and poles turned clear of them.
#[test]
#[ignore = "heavy"]
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
        (&cornered, [2, 0, 0, 1, 0]),
    ] {
        comes_back_as(&model, shape, expected);
        let mesh = ogeom::mesh::triangulate(&model, shape, Deflection::default(), T).unwrap();
        comes_back_from(&model, shape, &mesh, expected);
    }
}

/// A ring drilled across its tube: blind from outside through the outer
/// equator, where the tube's angle starts; blind from the hole through the
/// inner equator; straight through both; and blind from above through the
/// top parallel. Each comes back as the whole torus with the holes as
/// inner wires, its seam round the axis turned to a parallel no hole
/// crosses, the drilled walls cylinders and the blind ends planes. The
/// volume, measured on the exact surfaces both sides, agrees to a
/// millionth.
#[test]
#[ignore = "heavy"]
fn rings_drilled_across_their_equators_come_back_whole_but_for_the_holes() {
    let mut model = Model::new();
    let along_x =
        |x: f64| Frame::new(Point::new(x, 0.0, 0.0), Direction::X, Direction::Y, T).unwrap();
    let down_from_top =
        Frame::new(Point::new(10.0, 0.0, 1.0), Direction::Z, Direction::X, T).unwrap();
    for (name, frame, length, expected) in [
        ("outer equator", along_x(11.0), 5.0, [1, 1, 0, 0, 1]),
        ("inner equator", along_x(4.0), 5.0, [1, 1, 0, 0, 1]),
        ("both equators", along_x(5.0), 10.0, [0, 1, 0, 0, 1]),
        ("top parallel", down_from_top, 5.0, [1, 1, 0, 0, 1]),
    ] {
        let ring = ogeom::algo::make_torus(&mut model, Frame::WORLD, 10.0, 3.0, T)
            .unwrap()
            .shape;
        let drill = ogeom::algo::make_cylinder(&mut model, frame, 1.0, length, T)
            .unwrap()
            .shape;
        let drilled = ogeom::boolean::cut(&mut model, &ring, &drill, T)
            .unwrap()
            .shape;
        comes_back_as(&model, &drilled, expected);
        let mut back = Model::new();
        let out = solid_from_mesh(
            &mut back,
            &meshed(&model, &drilled),
            &MeshSolidOptions::default(),
            T,
        )
        .unwrap();
        let drawn = ogeom::mesh::triangulate(&back, &out.shape, Deflection::default(), T).unwrap();
        assert!(drawn.is_closed(), "{name}");
        let (a, b) = (volume(&model, &drilled), volume(&back, &out.shape));
        assert!(
            (a - b).abs() / a < 1e-6,
            "{name}: {a} went in, {b} came out"
        );
    }
}

/// Converts a ring's mesh and checks it comes back on `expected` surfaces
/// with no torus left as facets, valid, closed and meshing closed, and with
/// the volume, measured on the exact surfaces both sides, within a
/// millionth of the original's.
fn ring_comes_back_exact(model: &Model, shape: &Shape, expected: [usize; 5], name: &str) {
    let mut back = Model::new();
    let out = solid_from_mesh(
        &mut back,
        &meshed(model, shape),
        &MeshSolidOptions::default(),
        T,
    )
    .unwrap();
    assert!(out.closed, "{name}: {:?}", out.report);
    assert_eq!(kinds(&back, &out.shape), expected, "{name}");
    assert_eq!(out.report.curved_faceted, 0, "{name}");
    let diagnosis = check(&back, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
    let drawn = ogeom::mesh::triangulate(&back, &out.shape, Deflection::default(), T).unwrap();
    assert!(drawn.is_closed(), "{name}");
    let (a, b) = (volume(model, shape), volume(&back, &out.shape));
    assert!(
        (a - b).abs() / a < 1e-6,
        "{name}: {a} went in, {b} came out"
    );
}

fn ring(model: &mut Model) -> Shape {
    ogeom::algo::make_torus(model, Frame::WORLD, 10.0, 3.0, T)
        .unwrap()
        .shape
}

/// A ring grooved all the way round its axis by a second torus over its
/// outer equator: a band round the axis between two parallels, however
/// narrow the groove, and the groove's floor a band of the second torus.
/// And a ring grooved all the way round its tube by a flat-sided ring cut
/// about the tube: a band round the tube between two rims that are no
/// circles, a seam of its own joining them.
#[test]
fn rings_grooved_all_the_way_round_come_back_bands() {
    let mut model = Model::new();
    let rim = ring(&mut model);
    let cutter = ogeom::algo::make_torus(&mut model, Frame::WORLD, 13.0, 0.5, T)
        .unwrap()
        .shape;
    let round_axis = ogeom::boolean::cut(&mut model, &rim, &cutter, T)
        .unwrap()
        .shape;
    ring_comes_back_exact(&model, &round_axis, [0, 0, 0, 0, 2], "round the axis");

    let rim = ring(&mut model);
    let about_tube =
        Frame::new(Point::new(10.0, -0.25, 0.0), Direction::Y, Direction::X, T).unwrap();
    let outer = ogeom::algo::make_cylinder(&mut model, about_tube, 4.0, 0.5, T)
        .unwrap()
        .shape;
    let inner = ogeom::algo::make_cylinder(&mut model, about_tube, 2.5, 0.5, T)
        .unwrap()
        .shape;
    let cutter = ogeom::boolean::cut(&mut model, &outer, &inner, T)
        .unwrap()
        .shape;
    let round_tube = ogeom::boolean::cut(&mut model, &rim, &cutter, T)
        .unwrap()
        .shape;
    ring_comes_back_exact(&model, &round_tube, [2, 1, 0, 0, 1], "round the tube");
}

/// A ring with a long slot over its top on one side and under its bottom
/// on the other, the two overlapping round the axis: no meridian is free
/// of them. The torus comes back one face, its seam round the tube running
/// through a slot.
#[test]
fn a_ring_slotted_across_every_meridian_comes_back_one_torus() {
    let mut model = Model::new();
    let mut slotted = ring(&mut model);
    for (corner, size) in [
        (Point::new(0.0, -15.0, 1.5), (15.0, 30.0, 5.0)),
        (Point::new(-15.0, -15.0, -6.5), (16.0, 30.0, 5.0)),
    ] {
        let frame = Frame::new(corner, Direction::Z, Direction::X, T).unwrap();
        let cutter = ogeom::algo::make_box(&mut model, frame, size, T)
            .unwrap()
            .shape;
        slotted = ogeom::boolean::cut(&mut model, &slotted, &cutter, T)
            .unwrap()
            .shape;
    }
    ring_comes_back_exact(&model, &slotted, [6, 0, 0, 0, 1], "slots");
}

/// A ring drilled straight through its tube three times, a third of a
/// turn apart round the axis and a sixth of a turn apart round the tube:
/// the six holes between them cross every parallel. The torus comes back
/// one face, its seam round the axis running through a hole.
#[test]
#[ignore = "heavy"]
fn a_ring_drilled_across_every_parallel_comes_back_one_torus() {
    let mut model = Model::new();
    let mut drilled = ring(&mut model);
    for k in 0..3 {
        let round_axis = f64::from(k) * core::f64::consts::TAU / 3.0;
        let round_tube = f64::from(k) * core::f64::consts::PI / 3.0;
        let out = Vector::new(round_axis.cos(), round_axis.sin(), 0.0);
        let along = out * round_tube.cos() + Vector::Z * round_tube.sin();
        let centre = Point::new(10.0 * round_axis.cos(), 10.0 * round_axis.sin(), 0.0);
        let axis = Direction::new(along, T).unwrap();
        let frame = Frame::new(centre - along * 5.0, axis, axis.any_perpendicular(), T).unwrap();
        let drill = ogeom::algo::make_cylinder(&mut model, frame, 2.0, 10.0, T)
            .unwrap()
            .shape;
        drilled = ogeom::boolean::cut(&mut model, &drilled, &drill, T)
            .unwrap()
            .shape;
    }
    ring_comes_back_exact(&model, &drilled, [0, 3, 0, 0, 1], "drilled");
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
#[ignore = "heavy"]
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

/// A leg of a profile in the half-plane `(rho, h)`: a line to a point, or
/// an arc about a centre to a point, turning left or not.
#[derive(Clone, Copy)]
enum Leg {
    To(f64, f64),
    Arc((f64, f64), bool, (f64, f64)),
}

/// A closed profile from `start` along `legs`, drawn in the `xz` plane
/// (`x` the distance from the axis) and turned all the way round `z`.
fn turned_profile(model: &mut Model, start: (f64, f64), legs: &[Leg]) -> Shape {
    use ogeom::geom::{CircleCurve, Curve, Curve3d as _, LineCurve};
    let at = |p: (f64, f64)| Point::new(p.0, 0.0, p.1);
    let first = ogeom::algo::build::make_vertex(model, at(start)).shape;
    let (mut from, mut vertex) = (start, first.clone());
    let mut edges = Vec::new();
    for (k, leg) in legs.iter().enumerate() {
        let (Leg::To(x, y) | Leg::Arc(_, _, (x, y))) = *leg;
        let next = if k + 1 == legs.len() {
            first.clone()
        } else {
            ogeom::algo::build::make_vertex(model, at((x, y))).shape
        };
        let (curve, range): (Curve, (f64, f64)) = match *leg {
            Leg::To(..) => {
                let line: Curve = LineCurve::segment(at(from), at((x, y)), T).unwrap().into();
                let range = line.domain();
                (line, range)
            }
            Leg::Arc(c, left, _) => {
                // About -y the angle runs from x towards z.
                let normal = if left { -Direction::Y } else { Direction::Y };
                let frame = Frame::new(at(c), normal, Direction::X, T).unwrap();
                let angle = |p: (f64, f64)| {
                    let a = (p.1 - c.1).atan2(p.0 - c.0);
                    (if left { a } else { -a }).rem_euclid(core::f64::consts::TAU)
                };
                let (a0, mut a1) = (angle(from), angle((x, y)));
                if a1 <= a0 {
                    a1 += core::f64::consts::TAU;
                }
                let radius = (from.0 - c.0).hypot(from.1 - c.1);
                let circle = ogeom::math::Circle::new(frame, radius, T).unwrap();
                (CircleCurve::new(circle).into(), (a0, a1))
            }
        };
        edges.push(
            ogeom::algo::build::make_edge_between(model, curve, range, &vertex, &next, T)
                .unwrap()
                .shape,
        );
        (from, vertex) = ((x, y), next);
    }
    let wire = ogeom::algo::make_wire(model, &edges, T).unwrap().shape;
    let plane =
        ogeom::math::Plane::new(Frame::new(Point::ORIGIN, Direction::Y, Direction::X, T).unwrap());
    let face = ogeom::algo::make_face(
        model,
        ogeom::geom::PlaneSurface::new(plane).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    ogeom::algo::make_revolution(
        model,
        &face,
        ogeom::math::Axis::new(Point::ORIGIN, Direction::Z),
        core::f64::consts::TAU,
        T,
    )
    .unwrap()
    .shape
}

/// How far a torus's tube stands from tangency with a surface on its axis
/// (a plane square to it, or a cylinder, cone or sphere about it): the
/// distance from the tube's centre circle to the surface less the tube's
/// radius, in the half-plane through the axis. `None` off the axis.
fn tube_gap(torus: ogeom::math::Torus, other: &ogeom::geom::SurfaceGeometry) -> Option<f64> {
    use ogeom::geom::SurfaceGeometry as S;
    let (o, z) = (torus.frame().origin(), torus.frame().z().vector());
    let (major, minor) = (torus.major_radius(), torus.minor_radius());
    let half = |p: Point| {
        let w = p - o;
        let h = w.dot(z);
        ((w - z * h).magnitude(), h)
    };
    let along = |d: Direction| d.vector().cross(z).magnitude() <= 1e-9;
    let distance = match other {
        S::Plane(p) if along(p.plane().normal()) => half(p.plane().frame().origin()).1.abs(),
        S::Cylinder(c) if along(c.cylinder().frame().z()) => {
            let c = c.cylinder();
            (half(c.frame().origin()).0 <= 1e-9).then_some((major - c.radius()).abs())?
        }
        S::Cone(c) if along(c.cone().frame().z()) => {
            let c = c.cone();
            let (rho, h) = half(c.frame().origin());
            if rho > 1e-9 {
                return None;
            }
            // The cone's trace through (radius_at(0), h) and one unit on.
            let s = c.frame().z().vector().dot(z);
            let (dx, dy) = (c.radius_at(1.0) - c.radius_at(0.0), s);
            ((major - c.radius_at(0.0)) * dy + h * dx).abs() / dx.hypot(dy)
        }
        S::Sphere(s) => {
            let (rho, h) = half(s.sphere().centre());
            (rho <= 1e-9).then_some((major.hypot(h) - s.sphere().radius()).abs())?
        }
        _ => return None,
    };
    Some((distance - minor).abs())
}

/// Rounds between two surfaces of revolution about one axis come back on
/// the torus those surfaces fix, not on a free fit, so they meet them
/// tangentially to rounding: a round between a shaft and a taper (two
/// lines crossing in the half-plane through the axis), a tube's lip
/// rounded right across (two parallel lines) and a round under a dome
/// wider than its shaft (a line and a circle), each from a
/// single-precision mesh.
#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
#[test]
fn rounds_between_surfaces_of_revolution_meet_them_tangentially() {
    use ogeom::geom::SurfaceGeometry as S;
    let mut model = Model::new();
    // A shaft of radius 6 into a taper to radius 3, rounded 2 between them
    // and 1 at the foot: the tangent points worked out on the profile.
    let (cos, sin) = (8.0 / 73.0_f64.sqrt(), 3.0 / 73.0_f64.sqrt());
    let corner = 2.0 * (1.0 - cos) / sin;
    let taper = turned_profile(
        &mut model,
        (0.0, 0.0),
        &[
            Leg::To(5.0, 0.0),
            Leg::Arc((5.0, 1.0), true, (6.0, 1.0)),
            Leg::To(6.0, 10.0 - corner),
            Leg::Arc(
                (4.0, 10.0 - corner),
                true,
                (6.0 - 2.0 * (1.0 - cos), 10.0 - corner + 2.0 * sin),
            ),
            Leg::To(3.0, 18.0),
            Leg::To(0.0, 18.0),
            Leg::To(0.0, 0.0),
        ],
    );
    let lip = turned_profile(
        &mut model,
        (4.0, 0.0),
        &[
            Leg::To(6.0, 0.0),
            Leg::To(6.0, 10.0),
            Leg::Arc((5.0, 10.0), true, (4.0, 10.0)),
            Leg::To(4.0, 0.0),
        ],
    );
    // A dome of radius 8 on a shaft of radius 6, rounded 1.5 between them.
    let (dome, shaft, round) = (8.0_f64, 6.0_f64, 1.5_f64);
    let low = 10.0 - (dome * dome - shaft * shaft).sqrt();
    let centre = (
        shaft - round,
        low + ((dome - round).powi(2) - (shaft - round).powi(2)).sqrt(),
    );
    let out = dome / (dome - round);
    let domed = turned_profile(
        &mut model,
        (0.0, 0.0),
        &[
            Leg::To(shaft, 0.0),
            Leg::To(shaft, centre.1),
            Leg::Arc(centre, true, (centre.0 * out, low + (centre.1 - low) * out)),
            Leg::Arc((0.0, low), true, (0.0, low + dome)),
            Leg::To(0.0, 0.0),
        ],
    );
    for (name, part, tubes) in [("taper", &taper, 4), ("lip", &lip, 2), ("dome", &domed, 2)] {
        assert!(check(&model, part, T).unwrap().is_valid(), "{name}");
        let mut mesh = meshed(&model, part);
        for p in &mut mesh.positions {
            *p = Point::new(
                f64::from(p.x as f32),
                f64::from(p.y as f32),
                f64::from(p.z as f32),
            );
        }
        let mut back = Model::new();
        let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
        assert!(check(&back, &out.shape, T).unwrap().is_valid(), "{name}");
        let (surfaces, pairs) = surfaces_and_neighbours(&back, &out.shape);
        let mut met = 0;
        for (i, j) in pairs {
            for (a, b) in [(&surfaces[i], &surfaces[j]), (&surfaces[j], &surfaces[i])] {
                let S::Torus(t) = a else {
                    continue;
                };
                let gap = tube_gap(t.torus(), b)
                    .unwrap_or_else(|| panic!("{name}: a neighbour off the torus's axis"));
                assert!(gap < 1e-9, "{name}: a torus {gap:e} off tangency");
                met += 1;
            }
        }
        assert_eq!(met, tubes, "{name}: each round meets both its supports");
    }
}

/// A plate with rounded corners whose bottom rim is filleted, meshed
/// coarsely: where the rim's fillets turn the corners the mesh leaves
/// facets of a triangle or two that no surface claims, each bounded by
/// seams as wide as itself. They join the fillets beside them, and the
/// plate comes back its own eighteen faces, not fifty.
#[test]
fn facets_left_between_fillets_join_them() {
    let mut model = Model::new();
    let edge_ends = |model: &Model, e: &Shape| {
        let data = model.node(e).unwrap().data().as_edge().unwrap();
        let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            return None;
        };
        let c = model.geometry().curve(*curve).unwrap();
        Some((
            ogeom::geom::Curve3d::point_at(c, range.0, T).unwrap(),
            ogeom::geom::Curve3d::point_at(c, range.1, T).unwrap(),
        ))
    };
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (47.0, 21.5, 5.0), T)
        .unwrap()
        .shape;
    let upright: Vec<Shape> = explore_unique(&model, &block, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| edge_ends(&model, e).is_some_and(|(a, b)| (a.z - b.z).abs() > 4.0))
        .collect();
    let plate = ogeom::fillet::fillet_edges(&mut model, &block, &upright, 8.0, T)
        .unwrap()
        .shape;
    let rim: Vec<Shape> = explore_unique(&model, &plate, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| edge_ends(&model, e).is_some_and(|(a, b)| a.z.abs() < 1e-9 && b.z.abs() < 1e-9))
        .collect();
    let filleted = ogeom::fillet::fillet_edges(&mut model, &plate, &rim, 0.5, T)
        .unwrap()
        .shape;
    let mesh =
        ogeom::mesh::triangulate(&model, &filleted, Deflection::with_chord(0.05).unwrap(), T)
            .unwrap();
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(check(&back, &out.shape, T).unwrap().is_valid());
    assert_eq!(out.report.faces, 18, "{:?}", out.report);
    let drawn = ogeom::mesh::triangulate(&back, &out.shape, Deflection::default(), T).unwrap();
    assert!(drawn.is_closed());
}

/// A plate with rounded corners whose bottom rim is filleted, meshed in
/// single precision: each fillet along a side meets the bottom tangentially,
/// where the mesh's boundary between them wanders and can step past its
/// corner and back. The seam threaded along it still runs from corner to
/// corner without hooking past either, so no face's trim folds, and the
/// plate comes back its own eighteen faces.
#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
#[test]
#[ignore = "heavy"]
fn a_filleted_plate_s_tangent_seams_end_at_their_corners() {
    let mut model = Model::new();
    let edge_ends = |model: &Model, e: &Shape| {
        let data = model.node(e).unwrap().data().as_edge().unwrap();
        let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            return None;
        };
        let c = model.geometry().curve(*curve).unwrap();
        Some((
            ogeom::geom::Curve3d::point_at(c, range.0, T).unwrap(),
            ogeom::geom::Curve3d::point_at(c, range.1, T).unwrap(),
        ))
    };
    let block = ogeom::algo::make_box(
        &mut model,
        Frame::new(Point::new(104.5, 76.96, 0.0), Direction::Z, Direction::X, T).unwrap(),
        (47.0, 21.5, 5.0),
        T,
    )
    .unwrap()
    .shape;
    let upright: Vec<Shape> = explore_unique(&model, &block, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| edge_ends(&model, e).is_some_and(|(a, b)| (a.z - b.z).abs() > 1.0))
        .collect();
    let plate = ogeom::fillet::fillet_edges(&mut model, &block, &upright, 8.0, T)
        .unwrap()
        .shape;
    let rim: Vec<Shape> = explore_unique(&model, &plate, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| edge_ends(&model, e).is_some_and(|(a, b)| a.z.abs() < 1e-9 && b.z.abs() < 1e-9))
        .collect();
    let filleted = ogeom::fillet::fillet_edges(&mut model, &plate, &rim, 0.5, T)
        .unwrap()
        .shape;
    let deflection = Deflection {
        chord: 0.005,
        angular: 0.35,
        ..Deflection::default()
    };
    let mut mesh = ogeom::mesh::triangulate(&model, &filleted, deflection, T).unwrap();
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
    assert_eq!(out.report.faces, 18, "{:?}", out.report);
    // Every fitted seam stays between its ends: no point of it lies farther
    // from either end than the ends lie from each other, past its tolerance.
    for e in explore_unique(&back, &out.shape, ShapeType::Edge).unwrap() {
        let data = back.node(&e).unwrap().data().as_edge().unwrap();
        let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let c = back.geometry().curve(*curve).unwrap();
        if !matches!(c, ogeom::geom::Curve::BSpline(_)) {
            continue;
        }
        let at = |t: f64| ogeom::geom::Curve3d::point_at(c, t, T).unwrap();
        let (a, b) = (at(range.0), at(range.1));
        let span = a.distance(b) + data.tolerance.get();
        for k in 0..=200 {
            let p = at(range.0 + (range.1 - range.0) * f64::from(k) / 200.0);
            assert!(
                p.distance(a) <= span && p.distance(b) <= span,
                "{p:?} beyond {a:?} .. {b:?}"
            );
        }
    }
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
#[ignore = "heavy"]
fn fillets_ending_on_rough_corners_are_still_cylinders() {
    let mut back = Model::new();
    let mesh = rough_rounded_box();
    let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
    eprintln!("rough box: {:?}", out.report);
    assert!(out.closed);
    assert_eq!(
        kinds_and_patches(&back, &out.shape).0[1],
        12,
        "every fillet a cylinder"
    );
    assert_eq!(out.report.curved_faceted, 0);
    assert!(check(&back, &out.shape, T).unwrap().is_valid());
    // A single facet left beside a fillet is thinner than the curve the
    // two surfaces meet along bulges; their seam is threaded straight
    // instead, so the facet's trim does not fold over itself and the
    // whole shape tessellates closed.
    let drawn = ogeom::mesh::triangulate(&back, &out.shape, Deflection::default(), T).unwrap();
    assert!(drawn.is_closed());
    // Measured finely: tessellated at the default deflection the twelve
    // fitted cylinders lose about as much volume as the tolerance allows.
    let got = volume_properties(&back, &out.shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass;
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
#[ignore = "heavy"]
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
#[ignore = "heavy"]
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
/// one face of the torus it is. The bound allows for the rest of the suite
/// sharing the machine; a conversion that grows faster than the mesh takes
/// minutes.
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
    assert!(took < Duration::from_secs(30), "{took:?}");
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
#[ignore = "heavy"]
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
#[ignore = "heavy"]
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

/// A coplanar distance below the scatter of a single-precision mesh asks
/// for more than the data holds: no surface fitted to the vertices holds
/// them that close, and the regions break into fragments whose seams open.
/// It is raised to the scatter, the report says so, and the result
/// tessellates closed.
#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
#[test]
#[ignore = "heavy"]
fn a_distance_below_the_mesh_s_scatter_is_raised_and_said_so() {
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
    for asked in [1e-7, 5e-7] {
        let options = MeshSolidOptions {
            coplanar_distance: Some(asked),
            ..MeshSolidOptions::default()
        };
        let mut back = Model::new();
        let out = solid_from_mesh(&mut back, &mesh, &options, T).unwrap();
        assert!(out.report.coplanar_distance_raised);
        assert!(out.coplanar_distance > asked);
        assert!(check(&back, &out.shape, T).unwrap().is_valid());
        let drawn = ogeom::mesh::triangulate(&back, &out.shape, Deflection::default(), T).unwrap();
        assert!(drawn.is_closed(), "asked {asked}");
    }
    // Asked for at or above the scatter, the distance is kept.
    let options = MeshSolidOptions {
        coplanar_distance: Some(1e-5),
        ..MeshSolidOptions::default()
    };
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, &mesh, &options, T).unwrap();
    assert!(!out.report.coplanar_distance_raised);
    assert!((out.coplanar_distance - 1e-5).abs() < 1e-12);
}

#[allow(
    clippy::cast_possible_truncation,
    reason = "the rounding to single precision is the point"
)]
#[test]
#[ignore = "heavy"]
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
#[ignore = "heavy"]
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
#[ignore = "heavy"]
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

/// The area of a polygon in plan, positive when it runs anticlockwise.
fn plan_area(polygon: &[(f64, f64)]) -> f64 {
    let n = polygon.len();
    (0..n)
        .map(|i| {
            let (a, b) = (polygon[i], polygon[(i + 1) % n]);
            a.0 * b.1 - b.0 * a.1
        })
        .sum::<f64>()
        / 2.0
}

/// `subject` clipped to the convex anticlockwise polygon `window`.
fn clip_to_convex(subject: &[(f64, f64)], window: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut out = subject.to_vec();
    for i in 0..window.len() {
        let (a, b) = (window[i], window[(i + 1) % window.len()]);
        let side = |p: (f64, f64)| (b.0 - a.0) * (p.1 - a.1) - (b.1 - a.1) * (p.0 - a.0);
        let input = std::mem::take(&mut out);
        for k in 0..input.len() {
            let (p, q) = (input[k], input[(k + 1) % input.len()]);
            let (sp, sq) = (side(p), side(q));
            if sp >= 0.0 {
                out.push(p);
            }
            if (sp >= 0.0) != (sq >= 0.0) {
                let t = sp / (sp - sq);
                out.push((p.0 + t * (q.0 - p.0), p.1 + t * (q.1 - p.1)));
            }
        }
        if out.is_empty() {
            break;
        }
    }
    out
}

/// A slab drafted inward from its top by `draft` times the depth to the
/// power one and a half, one corner rounded in `counts` segments row by
/// row, tessellated and rounded to `f32` a hundred millimetres from the
/// origin, converted face for facet, and a pad of its top face pushed
/// `depth` down into it. With `rows_kept` the conversion's coplanar
/// distance is twice what single precision resolves, so each row stays a
/// face of its own however little it leans; without, it is the
/// converter's own.
///
/// The pad fused back, in common with the slab and cut from it must each
/// be valid, with a volume within `2e-5` cubic millimetres of one measured
/// from the mesh alone, plus the doubt the converted slab carries where the
/// pad's floor crosses it. Where rows that lean alike merge into one face
/// on a fitted plane, its vertices stand off the plane by up to their
/// tolerance; a face the floor splits is tessellated on those vertices in
/// pieces, not whole, and its volume moves by up to how far they stand off
/// the plane times its area. That sum over the faces the floor crosses is
/// the doubt, returned; faces holding their vertices to rounding carry
/// none. The slab's volume is its triangles' signed
/// tetrahedra, the pad's its lid's area times its depth, and their common
/// part the integral down the pad of the lid clipped by the slab's
/// section, which between two rows is the polygon through the edges
/// joining them (Simpson's rule, two thousand steps a band); the fuse and
/// the cut are what is left of the three. A result's volume is its own
/// tessellation's, which must close: on planes the tessellation stands on
/// the result's vertices, where the exact integral runs round each face's
/// own chart and so differs by what the vertices' tolerances allow.
fn pad_on_a_drafted_slab(
    depth: f64,
    draft: f64,
    counts: [u32; 7],
    rows_kept: bool,
) -> Result<f64, String> {
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
    // outline; each band keeps the edges joining its two rows in order.
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
    let mut rungs: Vec<Vec<(u32, u32)>> = Vec::new();
    for j in 0..rows.len() - 1 {
        let (lo, hi) = (&rings[j], &rings[j + 1]);
        let (at_lo, at_hi) = (place(lo, &mesh), place(hi, &mesh));
        let (nl, nh) = (lo.len(), hi.len());
        let (mut a, mut b) = (0_usize, 0_usize);
        let mut band = Vec::new();
        while a < nl || b < nh {
            band.push((lo[a % nl], hi[b % nh]));
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
        rungs.push(band);
    }
    let (bottom, upper) = (&rings[0], &rings[rows.len() - 1]);
    for i in 1..bottom.len() - 1 {
        mesh.triangles.push([bottom[0], bottom[i + 1], bottom[i]]);
    }
    for i in 1..upper.len() - 1 {
        mesh.triangles.push([upper[0], upper[i], upper[i + 1]]);
    }

    // The measurement from the mesh, about the slab's middle in plan.
    let local = |p: Point| Vector::new(p.x - 105.0, p.y - 85.0, p.z);
    let signed = |mesh: &Triangulation| -> f64 {
        mesh.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| local(mesh.positions[i as usize]));
                a.dot(b.cross(c)) / 6.0
            })
            .sum()
    };
    let slab_volume = signed(&mesh);
    let window: Vec<(f64, f64)> = upper
        .iter()
        .map(|&i| {
            let p = local(mesh.positions[i as usize]);
            (p.x, p.y)
        })
        .collect();
    let pad_volume = plan_area(&window) * depth;
    let section = |band: &[(u32, u32)], z: f64| -> f64 {
        let polygon: Vec<(f64, f64)> = band
            .iter()
            .map(|&(l, h)| {
                let (p, q) = (
                    local(mesh.positions[l as usize]),
                    local(mesh.positions[h as usize]),
                );
                let t = (z - p.z) / (q.z - p.z);
                (p.x + t * (q.x - p.x), p.y + t * (q.y - p.y))
            })
            .collect();
        plan_area(&clip_to_convex(&polygon, &window))
    };
    let floor = mesh.positions[upper[0] as usize].z - depth;
    let mut common_volume = 0.0;
    for (j, band) in rungs.iter().enumerate() {
        let from = mesh.positions[rings[j][0] as usize].z.max(floor);
        let to = mesh.positions[rings[j + 1][0] as usize].z;
        if to <= from {
            continue;
        }
        let steps = 2000_u32;
        let h = (to - from) / f64::from(steps);
        let mut sum = section(band, from) + section(band, to);
        for k in 1..steps {
            let weight = if k % 2 == 1 { 4.0 } else { 2.0 };
            sum += weight * section(band, from + h * f64::from(k));
        }
        common_volume += sum * h / 3.0;
    }

    let quantum = ogeom::algo::single_precision_quantum(&mesh);
    let options = MeshSolidOptions {
        recognize: false,
        keep_vertices: true,
        quantum: Some(quantum),
        coplanar_distance: rows_kept.then_some(2.0 * quantum),
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
    let fine = Deflection::with_chord(1e-3).unwrap();
    let mut doubt = 0.0;
    for face in explore_unique(&model, &slab, ShapeType::Face).unwrap() {
        let (on, normal) = ogeom::algo::face_normal(&model, &face, T).unwrap();
        let (mut off, mut low, mut high) = (0.0_f64, f64::MAX, f64::MIN);
        for vertex in explore_unique(&model, &face, ShapeType::Vertex).unwrap() {
            let p = model
                .node(&vertex)
                .unwrap()
                .data()
                .as_vertex()
                .unwrap()
                .point;
            off = off.max((p - on).dot(normal).abs());
            (low, high) = (low.min(p.z), high.max(p.z));
        }
        if low < floor && high > floor {
            let area = ogeom::algo::surface_properties(&model, &face, fine, T)
                .unwrap()
                .mass;
            doubt += off * area;
        }
    }
    for (name, expected) in [
        ("fuse", slab_volume + pad_volume - common_volume),
        ("common", common_volume),
        ("cut", slab_volume - common_volume),
    ] {
        let made = match name {
            "fuse" => ogeom::boolean::fuse(&mut model, &slab, &pad, T),
            "common" => ogeom::boolean::common(&mut model, &slab, &pad, T),
            _ => ogeom::boolean::cut(&mut model, &slab, &pad, T),
        }
        .map_err(|e| format!("depth {depth}: {name}: {e}"))?
        .shape;
        let diagnosis = check(&model, &made, T).unwrap();
        if !diagnosis.is_valid() {
            return Err(format!("depth {depth}: {name}: {diagnosis}"));
        }
        // A cut of a slab the pad holds whole leaves nothing, and nothing
        // is closed with no volume.
        let tessellated = ogeom::mesh::triangulate(&model, &made, fine, T).unwrap();
        if !tessellated.triangles.is_empty() && !tessellated.is_closed() {
            return Err(format!("depth {depth}: {name}: its tessellation is open"));
        }
        let v = signed(&tessellated);
        if (v - expected).abs() > 2e-5 + doubt {
            return Err(format!(
                "depth {depth}: {name}: volume {v}, measured {expected}, off by {:.2e}, \
                 the slab's doubt {doubt:.2e}",
                v - expected
            ));
        }
    }
    Ok(doubt)
}

/// Each wall of the pad stands on the facet just under the top edge they
/// share, the two at a hundredth of a radian, the facet's plane missing
/// the edge's ends by their rounding: the walls cross on that edge, not on
/// a line their planes' solve leaves a sliver under it.
#[test]
fn a_pad_into_a_drafted_single_precision_slab_fuses_back() {
    for depth in [3.0, 10.0] {
        pad_on_a_drafted_slab(depth, 0.02, [6; 7], false).unwrap();
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
        pad_on_a_drafted_slab(depth, 0.005, [8, 8, 8, 7, 7, 6, 5], false).unwrap();
    }
}

/// The slab drafted a thousandth by the rows' depth to the power one and a
/// half, less than the corner's chords sag, so the rows' edges cross the
/// pad's walls in plan and every corner facet meets a wall at a small
/// angle, row counts uniform, rising, falling, alternating and mixed. A rung
/// between two rows lies all but in a wall, and the sections on the facets
/// either side reach it a few tenths of a micron apart along it; a row's
/// edge crosses a wall a few hundredths of a radian off it, and the
/// sections either side reach it a tenth of a micron apart. Each pair is
/// one junction, and a stub of section between two such crossings is that
/// junction too.
#[test]
#[ignore = "heavy"]
fn pads_into_slabs_drafted_less_than_their_chords_sag_cut_and_fuse() {
    let families: [[u32; 7]; 15] = [
        [6; 7],
        [8; 7],
        [4; 7],
        [8, 8, 8, 7, 7, 6, 5],
        [5, 6, 7, 7, 8, 8, 8],
        [4, 5, 6, 7, 8, 9, 10],
        [10, 9, 8, 7, 6, 5, 4],
        [8, 7, 8, 7, 8, 7, 8],
        [6, 5, 6, 5, 6, 5, 6],
        [3, 4, 5, 6, 7, 8, 9],
        [12, 11, 10, 9, 8, 7, 6],
        [7, 7, 7, 6, 6, 6, 5],
        [8, 9, 9, 5, 4, 5, 5],
        [8, 5, 9, 9, 5, 8, 7],
        [6, 6, 6, 4, 8, 7, 5],
    ];
    let mut failed = Vec::new();
    for counts in families {
        for depth in [3.0, 10.0] {
            if let Err(e) = pad_on_a_drafted_slab(depth, 0.001, counts, true) {
                failed.push(format!("{counts:?}: {e}"));
            }
        }
    }
    assert!(failed.is_empty(), "{failed:#?}");
}

/// The same slabs converted at the converter's own coplanar distance. The
/// rows lean by under a twentieth of a degree more at each row down, as
/// little as the scatter of a mesh's flat faces, but steadily: the scatter
/// estimate reads them as a curve drawn finely, the distance stays at what
/// single precision resolves, and each row stays a face of its own.
#[test]
#[ignore = "heavy"]
fn pads_into_slabs_drafted_less_than_their_chords_sag_at_the_default_distance() {
    let families: [[u32; 7]; 15] = [
        [6; 7],
        [8; 7],
        [4; 7],
        [8, 8, 8, 7, 7, 6, 5],
        [5, 6, 7, 7, 8, 8, 8],
        [4, 5, 6, 7, 8, 9, 10],
        [10, 9, 8, 7, 6, 5, 4],
        [8, 7, 8, 7, 8, 7, 8],
        [6, 5, 6, 5, 6, 5, 6],
        [3, 4, 5, 6, 7, 8, 9],
        [12, 11, 10, 9, 8, 7, 6],
        [7, 7, 7, 6, 6, 6, 5],
        [8, 9, 9, 5, 4, 5, 5],
        [8, 5, 9, 9, 5, 8, 7],
        [6, 6, 6, 4, 8, 7, 5],
    ];
    let mut failed = Vec::new();
    for counts in families {
        for depth in [3.0, 10.0] {
            if let Err(e) = pad_on_a_drafted_slab(depth, 0.001, counts, false) {
                failed.push(format!("{counts:?}: {e}"));
            }
        }
    }
    assert!(failed.is_empty(), "{failed:#?}");
}

/// The slab drafted 0.00002 to 0.0007, every row holding the same count, at
/// both coplanar distances: each top facet stands parallel to its wall in
/// plan, leaning off it by under a thousandth of a radian, down to about
/// two hundred-thousandths, and their planes' solve crosses microns to a
/// tenth of a millimetre off the top edge they share to their points'
/// rounding. The edge is where they cross: the facet's far side stands off
/// the wall by more than the weld distance, though the band within it runs
/// a third of the way down the facet.
#[test]
fn pads_into_slabs_drafted_under_a_thousandth_with_uniform_rows_cut_and_fuse() {
    let mut failed = Vec::new();
    for draft in [0.0007, 0.0005, 0.0003, 0.0001, 0.00005, 0.00002] {
        for counts in [[8; 7], [6; 7], [4; 7]] {
            for rows_kept in [true, false] {
                for depth in [3.0, 10.0] {
                    if let Err(e) = pad_on_a_drafted_slab(depth, draft, counts, rows_kept) {
                        failed.push(format!("{draft} {counts:?} {rows_kept}: {e}"));
                    }
                }
            }
        }
    }
    assert!(failed.is_empty(), "{failed:#?}");
}

/// The slab drafted 0.0001 to 0.0007 with rows of differing counts,
/// converted at the converter's own coplanar distance, and the pad's floor
/// crossing it below the top band or standing below the slab. The rows'
/// turns grow toward the top as the draft steepens, and in the corner the
/// nearly flat edges between rows lie among edges that turn far more:
/// neither is scatter, the distance stays at what single precision
/// resolves, and each row stays a face of its own. The faces the floor
/// splits hold their vertices, so the results' volumes are the mesh's to
/// within a few hundred-thousandths of a cubic millimetre, the doubt among
/// them.
///
/// Drafted 0.0005, a rung between two rows of equal count leans out of the
/// pad's wall by a couple of thousandths of a radian, and the sections on
/// the facets either side reach it a micron apart along it: one
/// junction. Drafted 0.0003, a row's vertex stands a sixth of a micron
/// outside a wall, and the facets round it poke through the wall in a ring
/// a micron and a third wide on average, which the wall's mesh keeps as a
/// hole.
///
/// Drafted 0.00005, rows merge into faces whose vertices stand off them by
/// up to a hundredth of a micron, the doubt a thousandth of a cubic
/// millimetre. Where the pad's floor crosses a corner row, the floor keeps
/// a sliver a quarter of a micron wide between its own edge and the
/// section, the section's end vertex widened to most of a micron where it
/// was welded higher up: the sliver's two sides are two edges, not one,
/// and the wall above keeps its corner.
#[test]
fn pads_into_slabs_drafted_under_a_thousandth_with_differing_rows_cut_and_fuse() {
    let mut failed = Vec::new();
    for (draft, counts, depth) in [
        (0.0007, [6, 5, 6, 5, 6, 5, 6], 3.0),
        (0.0007, [8, 9, 9, 5, 4, 5, 5], 3.0),
        (0.0007, [8, 5, 9, 9, 5, 8, 7], 3.0),
        (0.0005, [8, 7, 8, 7, 8, 7, 8], 3.0),
        (0.0005, [6, 5, 6, 5, 6, 5, 6], 3.0),
        (0.0005, [8, 9, 9, 5, 4, 5, 5], 3.0),
        (0.0005, [7, 7, 7, 6, 6, 6, 5], 10.0),
        (0.0005, [8, 8, 8, 7, 7, 6, 5], 10.0),
        (0.0003, [3, 4, 5, 6, 7, 8, 9], 10.0),
        (0.0003, [6, 6, 6, 4, 8, 7, 5], 3.0),
        (0.0003, [6, 6, 6, 4, 8, 7, 5], 10.0),
        (0.0001, [8, 9, 9, 5, 4, 5, 5], 3.0),
        (0.0001, [6, 6, 6, 4, 8, 7, 5], 3.0),
        (0.00005, [8, 8, 8, 7, 7, 6, 5], 3.0),
    ] {
        match pad_on_a_drafted_slab(depth, draft, counts, false) {
            Err(e) => failed.push(format!("{draft} {counts:?}: {e}")),
            Ok(doubt) if draft >= 1e-4 && doubt > 5e-5 => {
                failed.push(format!("{draft} {counts:?}: the slab's doubt {doubt:.2e}"));
            }
            Ok(_) => {}
        }
    }
    assert!(failed.is_empty(), "{failed:#?}");
}

/// The slab drafted a hundred-thousandth, each row a face of its own: a top
/// facet leans off the pad's wall by about a hundred-thousandth of a radian.
/// With eight facets to a corner row, the lower triangle of each top quad
/// meets its wall only at a vertex, and its plane, fitted with the rows
/// below, crosses the wall's on a line the rounding places a tenth of a
/// millimetre down the facet; rising counts leave a facet's side edge beside
/// the next wall to within the rounding. Each face lies on one side of the
/// other's plane, and neither line splits anything.
#[test]
fn pads_into_a_slab_drafted_a_hundred_thousandth_with_rows_kept_cut_and_fuse() {
    let mut failed = Vec::new();
    for counts in [[8; 7], [5, 6, 7, 7, 8, 8, 8]] {
        for depth in [3.0, 10.0] {
            if let Err(e) = pad_on_a_drafted_slab(depth, 0.00001, counts, true) {
                failed.push(format!("{counts:?}: {e}"));
            }
        }
    }
    assert!(failed.is_empty(), "{failed:#?}");
}

/// The same slab at the converter's own coplanar distance, where rows merge
/// into faces on fitted planes: the top facets touch their walls only at
/// corners shared to rounding, and are cut and fused alike.
#[test]
fn pads_into_a_slab_drafted_a_hundred_thousandth_cut_and_fuse() {
    let mut failed = Vec::new();
    for (counts, depth) in [
        ([6; 7], 3.0),
        ([6; 7], 10.0),
        ([4; 7], 3.0),
        ([8; 7], 3.0),
        ([5, 6, 7, 7, 8, 8, 8], 10.0),
    ] {
        if let Err(e) = pad_on_a_drafted_slab(depth, 0.00001, counts, false) {
            failed.push(format!("{counts:?}: {e}"));
        }
    }
    assert!(failed.is_empty(), "{failed:#?}");
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
/// surface, with the pcurves it went out with and the same volume.
#[test]
#[ignore = "heavy"]
fn converted_solids_come_back_through_step() {
    for (name, model, shape) in converted_parts() {
        let before = (kinds_and_patches(&model, &shape), volume(&model, &shape));
        let mut document = ogeom::doc::Document::over(model.clone());
        document.add_part(name, shape.clone());
        let text = ogeom::io::write_step(&document, T).unwrap();
        let import = ogeom::io::read_step(&text, T).unwrap();
        let back = import.document.model();
        let [solid] = import.solids.as_slice() else {
            panic!("{name}: {} solids came back", import.solids.len());
        };
        let diagnosis = check(back, solid, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
        assert_eq!(kinds_and_patches(back, solid), before.0, "{name}");
        let (same, total) = pcurves::kept((&model, &shape), (back, solid), T);
        assert_eq!(
            same, total,
            "{name}: the pcurves come back as they went out"
        );
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
            assert_eq!(edges_walked_one_way(m, edited), 0);
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
#[ignore = "heavy"]
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

/// A solid over an `n` by `n` grid of cells on a square of side `size`:
/// its top lifted by `height`, its bottom flat at nought, a vertical wall
/// along every side between a cell kept and one not (or the square's
/// edge). Only the cells `keep` names are solid.
fn grid_solid(
    n: u32,
    size: f64,
    keep: impl Fn(u32, u32) -> bool,
    height: impl Fn(f64, f64) -> f64,
) -> Triangulation {
    let mut mesh = Triangulation::new();
    let at = |i: u32| size * f64::from(i) / f64::from(n);
    for j in 0..=n {
        for i in 0..=n {
            mesh.positions
                .push(Point::new(at(i), at(j), height(at(i), at(j))));
        }
    }
    for j in 0..=n {
        for i in 0..=n {
            mesh.positions.push(Point::new(at(i), at(j), 0.0));
        }
    }
    let top = |i: u32, j: u32| j * (n + 1) + i;
    let bottom = |i: u32, j: u32| (n + 1) * (n + 1) + j * (n + 1) + i;
    let kept = |i: i64, j: i64| {
        i >= 0
            && j >= 0
            && i < i64::from(n)
            && j < i64::from(n)
            && keep(u32::try_from(i).unwrap(), u32::try_from(j).unwrap())
    };
    for j in 0..n {
        for i in 0..n {
            if !keep(i, j) {
                continue;
            }
            let (a, b, c, d) = ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1));
            let t = |p: (u32, u32)| top(p.0, p.1);
            let w = |p: (u32, u32)| bottom(p.0, p.1);
            mesh.triangles.push([t(a), t(b), t(c)]);
            mesh.triangles.push([t(a), t(c), t(d)]);
            mesh.triangles.push([w(a), w(c), w(b)]);
            mesh.triangles.push([w(a), w(d), w(c)]);
            // Each side of the cell, as the top runs it, with the cell
            // across it.
            let (x, y) = (i64::from(i), i64::from(j));
            for (p, q, across) in [
                (a, b, (x, y - 1)),
                (b, c, (x + 1, y)),
                (c, d, (x, y + 1)),
                (d, a, (x - 1, y)),
            ] {
                if !kept(across.0, across.1) {
                    mesh.triangles.push([t(q), t(p), w(p)]);
                    mesh.triangles.push([t(q), w(p), w(q)]);
                }
            }
        }
    }
    mesh
}

/// A bump on a square plate twenty across: the height of a bicubic
/// B-spline over four spans each way, its control heights a lopsided hill
/// with a twist, so it is no surface of revolution, extrusion or canonical
/// surface.
fn bump(x: f64, y: f64) -> f64 {
    let (s, t) = (x / 20.0, y / 20.0);
    5.0 + 2.5
        * (core::f64::consts::PI * s).sin()
        * (core::f64::consts::PI * t).sin()
        * (1.0 + 0.4 * s)
        + 0.8 * (s - 0.5) * (t - 0.5)
}

/// The B-spline faces of a shape.
fn spline_faces(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            matches!(
                model.geometry().surface(data.surface).unwrap(),
                ogeom::geom::SurfaceGeometry::BSpline(_)
            )
        })
        .collect()
}

/// A bumped plate's top, meshed from its exact height and converted, comes
/// back one fitted B-spline face: within the coplanar distance of the exact
/// surface at points between the mesh's vertices both ways (the exact
/// surface's points to the patch, and the patch's to the exact surface),
/// the solid valid, a hundred times nearer the exact volume than the mesh
/// is, and tessellating closed.
#[test]
#[ignore = "heavy"]
fn a_free_form_bump_comes_back_one_fitted_patch() {
    let mesh = grid_solid(40, 20.0, |_, _| true, bump);
    let mut model = Model::new();
    let started = Instant::now();
    let out = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    let took = started.elapsed();
    let diagnosis = check(&model, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let patches = spline_faces(&model, &out.shape);
    eprintln!(
        "bump: {} faces, {} patches in {took:?}, distance {:e}, report {:?}",
        out.report.faces,
        patches.len(),
        out.coplanar_distance,
        out.report
    );
    assert_eq!(patches.len(), 1);
    assert_eq!(out.report.patch_faces, 1);
    assert_eq!(out.report.faces, 6);
    let face = &patches[0];
    let data = model.node(face).unwrap().data().as_face().unwrap();
    let surface = model.geometry().surface(data.surface).unwrap().clone();
    let flat = out.coplanar_distance;
    // The exact surface's points, a third of a cell off the mesh's grid, to
    // the patch.
    let mut worst: f64 = 0.0;
    for j in 0..60 {
        for i in 0..60 {
            let (x, y) = (
                20.0 * (f64::from(i) + 0.37) / 60.0,
                20.0 * (f64::from(j) + 0.61) / 60.0,
            );
            let exact = Point::new(x, y, bump(x, y));
            let foot = ogeom::algo::project_on_surface(&surface, exact, 16, T).unwrap();
            worst = worst.max(foot.distance);
        }
    }
    // The patch's points over the face to the exact surface: the height
    // between them over the slope's secant, the distance to first order.
    let drawn =
        ogeom::mesh::triangulate(&model, face, Deflection::with_chord(0.01).unwrap(), T).unwrap();
    let mut back: f64 = 0.0;
    for p in &drawn.positions {
        let h = 1e-6;
        let gx = (bump(p.x + h, p.y) - bump(p.x - h, p.y)) / (2.0 * h);
        let gy = (bump(p.x, p.y + h) - bump(p.x, p.y - h)) / (2.0 * h);
        back = back.max((p.z - bump(p.x, p.y)).abs() / (1.0 + gx * gx + gy * gy).sqrt());
    }
    eprintln!("bump: exact to patch {worst:e}, patch to exact {back:e}, distance {flat:e}");
    assert!(worst <= flat, "{worst} past {flat}");
    assert!(back <= flat, "{back} past {flat}");
    let volume = volume_properties(&model, &out.shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass;
    // The exact volume: the plate's, and the hill's, whose sines integrate
    // to 2.4 / pi and 2 / pi over the square; the twist integrates to none.
    let exact = 2000.0 + 1000.0 * 4.8 / core::f64::consts::PI.powi(2);
    let mesh_volume = mesh.volume();
    eprintln!("bump: volume {volume}, the mesh's {mesh_volume}, exact {exact}");
    assert!((volume - exact).abs() * 100.0 < (mesh_volume - exact).abs());
    let closed = ogeom::mesh::triangulate(&model, &out.shape, Deflection::default(), T).unwrap();
    assert!(closed.is_closed());
}

/// The heading of an S at arc length `s`: its curvature runs smoothly from
/// a quarter one way to a fifth the other, through a tanh two long about
/// fourteen along, so it turns back on itself by about 180 degrees and then
/// by about 210 the other way.
fn s_heading(s: f64) -> f64 {
    let lc = |x: f64| x.cosh().ln();
    -0.025 * s + 0.45 * (lc((s - 14.0) / 2.0) - lc(-7.0))
}

/// The plan of the S at arc length `s` from its start (at the origin,
/// heading along `x`): the point, by Simpson's rule on its heading, and
/// the unit normal on its left, the side the wall's thickness lies on.
fn s_plan(s: f64) -> (Point, Vector) {
    let steps = 2
        * (1..)
            .find(|k| f64::from(*k) * 0.02 >= s)
            .unwrap_or(1)
            .max(1);
    let h = s / f64::from(steps);
    let (mut x, mut y) = (0.0, 0.0);
    for k in 0..=steps {
        let w = if k == 0 || k == steps {
            1.0
        } else if k % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let a = s_heading(h * f64::from(k));
        x += w * a.cos();
        y += w * a.sin();
    }
    let a = s_heading(s);
    (
        Point::new(x * h / 3.0, y * h / 3.0, 0.0),
        Vector::new(-a.sin(), a.cos(), 0.0),
    )
}

/// The S's length.
fn s_length() -> f64 {
    34.0
}

/// A wall one thick along the S, eight high, its plan scaled about the
/// point where the S's two end normals meet by a factor growing with the
/// square of the height: so its two sides are free-form, and its ends
/// (each in a vertical plane through that point) and its top and bottom
/// are flat. `nt` cells along the S, `nz` up it.
fn s_wall(nt: u32, nz: u32) -> (Triangulation, impl Fn(f64, f64) -> Point) {
    let height = 8.0;
    // Where the end normals meet.
    let (a, na) = s_plan(0.0);
    let (b, nb) = s_plan(s_length());
    let det = na.x * (-nb.y) - na.y * (-nb.x);
    let k = ((b.x - a.x) * (-nb.y) - (b.y - a.y) * (-nb.x)) / det;
    let centre = a + na * k;
    let place = move |p: Point, z: f64| {
        let scale = 1.0 + 0.15 * (z / height).powi(2);
        let q = centre + (p - centre) * scale;
        Point::new(q.x, q.y, z)
    };
    let front = move |s: f64, z: f64| place(s_plan(s).0, z);
    let mut mesh = Triangulation::new();
    for side in 0..2 {
        for k in 0..=nz {
            for i in 0..=nt {
                let s = s_length() * f64::from(i) / f64::from(nt);
                let z = height * f64::from(k) / f64::from(nz);
                let (p, n) = s_plan(s);
                let p = if side == 0 { p } else { p + n };
                mesh.positions.push(place(p, z));
            }
        }
    }
    let at = |side: u32, i: u32, k: u32| side * (nt + 1) * (nz + 1) + k * (nt + 1) + i;
    let (f, w) = (|i, k| at(0, i, k), |i, k| at(1, i, k));
    for k in 0..nz {
        for i in 0..nt {
            mesh.triangles.push([f(i, k), f(i + 1, k), f(i + 1, k + 1)]);
            mesh.triangles.push([f(i, k), f(i + 1, k + 1), f(i, k + 1)]);
            mesh.triangles.push([w(i, k), w(i + 1, k + 1), w(i + 1, k)]);
            mesh.triangles.push([w(i, k), w(i, k + 1), w(i + 1, k + 1)]);
        }
    }
    for i in 0..nt {
        mesh.triangles.push([f(i + 1, 0), f(i, 0), w(i, 0)]);
        mesh.triangles.push([f(i + 1, 0), w(i, 0), w(i + 1, 0)]);
        mesh.triangles.push([f(i, nz), f(i + 1, nz), w(i + 1, nz)]);
        mesh.triangles.push([f(i, nz), w(i + 1, nz), w(i, nz)]);
    }
    for k in 0..nz {
        mesh.triangles.push([f(0, k), f(0, k + 1), w(0, k + 1)]);
        mesh.triangles.push([f(0, k), w(0, k + 1), w(0, k)]);
        mesh.triangles.push([f(nt, k), w(nt, k + 1), f(nt, k + 1)]);
        mesh.triangles.push([f(nt, k), w(nt, k), w(nt, k + 1)]);
    }
    (mesh, front)
}

/// Whether the triangles fold over the plane the points lie nearest: their
/// shadows on it wind both ways.
fn folds_over_its_plane(points: &[Point], triangles: &[[Point; 3]]) -> bool {
    let n = f64::from(u32::try_from(points.len()).unwrap());
    let mean = points.iter().fold(Vector::ZERO, |s, p| s + p.to_vector()) / n;
    let mut c = [[0.0_f64; 3]; 3];
    for p in points {
        let d = p.to_vector() - mean;
        let d = [d.x, d.y, d.z];
        for r in 0..3 {
            for k in 0..3 {
                c[r][k] += d[r] * d[k];
            }
        }
    }
    let apply = |c: &[[f64; 3]; 3], v: Vector| {
        Vector::new(
            c[0][0] * v.x + c[0][1] * v.y + c[0][2] * v.z,
            c[1][0] * v.x + c[1][1] * v.y + c[1][2] * v.z,
            c[2][0] * v.x + c[2][1] * v.y + c[2][2] * v.z,
        )
    };
    // The two widest directions by power iteration, the second with the
    // first taken out.
    let widest = |c: &[[f64; 3]; 3], start: Vector| {
        let mut v = start;
        for _ in 0..500 {
            let w = apply(c, v);
            v = w / w.magnitude();
        }
        v
    };
    let e1 = widest(&c, Vector::new(0.3, 0.9, 0.2));
    let l1 = apply(&c, e1).dot(e1);
    let mut d = c;
    let e = [e1.x, e1.y, e1.z];
    for r in 0..3 {
        for k in 0..3 {
            d[r][k] -= l1 * e[r] * e[k];
        }
    }
    let e2 = widest(&d, Vector::new(0.2, -0.1, 0.9));
    let e2 = (e2 - e1 * e2.dot(e1)) / (e2 - e1 * e2.dot(e1)).magnitude();
    let mut signs = [false; 2];
    for [a, b, c] in triangles {
        let (x, y) = (*b - *a, *c - *a);
        let area = x.dot(e1) * y.dot(e2) - x.dot(e2) * y.dot(e1);
        signs[usize::from(area > 0.0)] = true;
    }
    signs[0] && signs[1]
}

/// A thick S-shaped wall, its plan widening with height: each of its two
/// sides is one smooth region that folds over the plane it lies nearest
/// (no projection onto a plane charts it), so its chart comes from the
/// mean-value map onto a square, and the patch fitted on it still verifies.
/// Both sides come back fitted patches, the solid valid and tessellating
/// closed. Between the mesh's vertices the front patch keeps far closer to
/// the exact wall than the facets it was fitted to.
#[test]
#[ignore = "heavy"]
fn a_folding_s_wall_is_charted_by_the_mean_value_map() {
    let (nt, nz) = (160, 24);
    let (mesh, front) = s_wall(nt, nz);
    // The front side's vertices, first in the mesh, and its triangles, the
    // first two of the four each cell lists.
    let side = 4 * (nt * nz) as usize;
    let corners: Vec<[Point; 3]> = mesh.triangles[..side]
        .iter()
        .enumerate()
        .filter(|(k, _)| k % 4 < 2)
        .map(|(_, t)| t.map(|v| mesh.positions[v as usize]))
        .collect();
    let vertices = &mesh.positions[..((nt + 1) * (nz + 1)) as usize];
    assert!(folds_over_its_plane(vertices, &corners));

    let mut model = Model::new();
    let started = Instant::now();
    let out = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    let took = started.elapsed();
    eprintln!(
        "wall: {} faces, {} triangles in {took:?}, distance {:e}, report {:?}",
        out.report.faces,
        mesh.triangles.len(),
        out.coplanar_distance,
        out.report
    );
    let diagnosis = check(&model, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let patches = spline_faces(&model, &out.shape);
    assert_eq!(patches.len(), 2);
    assert_eq!(out.report.patch_faces, 2);
    assert_eq!(out.report.patch_charts_mapped, 2);
    assert_eq!(out.report.faces, 6);
    // The exact front side between the mesh's vertices, to the nearer
    // patch.
    let surfaces: Vec<ogeom::geom::SurfaceGeometry> = patches
        .iter()
        .map(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            model.geometry().surface(data.surface).unwrap().clone()
        })
        .collect();
    let flat = out.coplanar_distance;
    let (mut worst, mut facets): (f64, f64) = (0.0, 0.0);
    for k in 0..nz {
        for i in 0..nt {
            let s = s_length() * (f64::from(i) + 0.37) / f64::from(nt);
            let z = 8.0 * (f64::from(k) + 0.61) / f64::from(nz);
            let exact = front(s, z);
            let near = surfaces
                .iter()
                .map(|g| {
                    ogeom::algo::project_on_surface(g, exact, 16, T)
                        .unwrap()
                        .distance
                })
                .fold(f64::INFINITY, f64::min);
            worst = worst.max(near);
            // The cell's two facets on the front, the first two of the
            // four triangles the cell lists.
            let cell = 4 * (k * nt + i) as usize;
            let off = mesh.triangles[cell..cell + 2]
                .iter()
                .map(|t| {
                    let [a, b, c] = t.map(|v| mesh.positions[v as usize]);
                    let n = (b - a).cross(c - a);
                    ((exact - a).dot(n) / n.magnitude()).abs()
                })
                .fold(f64::INFINITY, f64::min);
            facets = facets.max(off);
        }
    }
    eprintln!("wall: exact to patch {worst:e}, to the facets {facets:e}, distance {flat:e}");
    assert!(
        worst * 10.0 < facets,
        "{worst} against the facets' {facets}"
    );
    let closed = ogeom::mesh::triangulate(&model, &out.shape, Deflection::default(), T).unwrap();
    assert!(closed.is_closed());
}

/// A bumped plate with a square hole through it: its top is one smooth
/// region no surface fits, but a ring round the hole, not a disk, so no
/// patch is tried. It stays faceted, counted, and the solid is valid.
#[test]
fn a_smooth_region_round_a_hole_stays_faceted() {
    let mesh = grid_solid(
        40,
        20.0,
        |i, j| !(16..24).contains(&i) || !(16..24).contains(&j),
        bump,
    );
    let mut model = Model::new();
    let out = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    eprintln!("ring: {} faces, report {:?}", out.report.faces, out.report);
    let diagnosis = check(&model, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    assert!(spline_faces(&model, &out.shape).is_empty());
    assert_eq!(out.report.patch_faces, 0);
    assert_eq!(out.report.patches_not_disk, 1);
    assert_eq!(out.report.patches_unverified, 0);
}

/// A plate twenty across, flat at five, with a lopsided twisted hill on the
/// middle twelve of each side that runs out into the flat tangentially:
/// height and slope are continuous where they meet, the curvature is not,
/// as along a round.
fn tangent_hill(x: f64, y: f64) -> f64 {
    if !(4.0..=16.0).contains(&x) || !(4.0..=16.0).contains(&y) {
        return 5.0;
    }
    let (s, t) = ((x - 4.0) / 12.0, (y - 4.0) / 12.0);
    let w = |a: f64| (core::f64::consts::PI * a).sin().powi(2);
    5.0 + w(s) * w(t) * (0.5 * (1.0 + 0.4 * s) + 0.5 * (s - 0.5) * (t - 0.5))
}

/// A hill running out tangentially into a flat top, meshed with vertices
/// along the line where they meet: no crease bounds the hill, and its smooth
/// region stops at the plane round it instead. It comes back one fitted
/// patch inside one plane, the seam between them a closed curve on the
/// plane within the coplanar distance of the patch. The fit holds the
/// vertices to half that distance, and between them, where its knots are
/// sparser than the vertices, the exact hill keeps within twice it. The
/// volume is within the distance over the hill's footprint, and the solid
/// valid and tessellating closed.
#[test]
fn a_hill_tangent_to_a_flat_comes_back_a_patch_inside_the_plane() {
    use ogeom::geom::Curve3d as _;
    let mesh = grid_solid(40, 20.0, |_, _| true, tangent_hill);
    let mut model = Model::new();
    let out = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    eprintln!(
        "tangent hill: {} faces, distance {:e}, report {:?}",
        out.report.faces, out.coplanar_distance, out.report
    );
    let diagnosis = check(&model, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let patches = spline_faces(&model, &out.shape);
    assert_eq!(patches.len(), 1);
    assert_eq!(out.report.faces, 7);
    let flat = out.coplanar_distance;
    let face = &patches[0];
    let edges = explore_unique(&model, face, ShapeType::Edge).unwrap();
    assert_eq!(edges.len(), 1);
    let data = model.node(&edges[0]).unwrap().data().as_edge().unwrap();
    let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        panic!("a seam without a curve");
    };
    let seam = model.geometry().curve(*curve).unwrap().clone();
    let tolerance = data.tolerance.get();
    let mut height: f64 = 0.0;
    for k in 0..=400 {
        let t = range.0 + (range.1 - range.0) * f64::from(k) / 400.0;
        height = height.max((seam.point_at(t, T).unwrap().z - 5.0).abs());
    }
    eprintln!("tangent hill: seam tolerance {tolerance:e}, off the plane {height:e}");
    assert!(tolerance <= flat, "{tolerance} past {flat}");
    assert!(height <= 1e-9, "{height}");
    let data = model.node(face).unwrap().data().as_face().unwrap();
    let surface = model.geometry().surface(data.surface).unwrap().clone();
    let mut worst: f64 = 0.0;
    for j in 0..36 {
        for i in 0..36 {
            let (x, y) = (
                4.0 + 12.0 * (f64::from(i) + 0.37) / 36.0,
                4.0 + 12.0 * (f64::from(j) + 0.61) / 36.0,
            );
            let exact = Point::new(x, y, tangent_hill(x, y));
            let foot = ogeom::algo::project_on_surface(&surface, exact, 16, T).unwrap();
            worst = worst.max(foot.distance);
        }
    }
    let volume = volume_properties(&model, &out.shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass;
    // The plate, and the hill: its sines squared integrate to a half each
    // way, the lean to 0.6 along x, the twist to nothing.
    let exact = 2000.0 + 0.5 * 0.6 * 0.5 * 144.0;
    eprintln!("tangent hill: exact to patch {worst:e}; volume {volume}, exact {exact}");
    assert!(worst <= flat * 2.0, "{worst} past {flat}");
    assert!((volume - exact).abs() <= flat * 144.0);
    let closed = ogeom::mesh::triangulate(&model, &out.shape, Deflection::default(), T).unwrap();
    assert!(closed.is_closed());
}

/// The lopsided twisted hill of [`tangent_hill`] over the unit square, its
/// height and slope nought along the square's sides.
fn unit_hill(s: f64, t: f64) -> f64 {
    if !(0.0..=1.0).contains(&s) || !(0.0..=1.0).contains(&t) {
        return 0.0;
    }
    let w = |a: f64| (core::f64::consts::PI * a).sin().powi(2);
    w(s) * w(t) * (0.5 * (1.0 + 0.4 * s) + 0.5 * (s - 0.5) * (t - 0.5))
}

/// The edges of a shape's fitted patches: each one's curve, range and
/// tolerance.
fn patch_edges(model: &Model, shape: &Shape) -> Vec<(ogeom::geom::Curve, (f64, f64), f64)> {
    let mut out = Vec::new();
    for face in spline_faces(model, shape) {
        for edge in explore_unique(model, &face, ShapeType::Edge).unwrap() {
            let data = model.node(&edge).unwrap().data().as_edge().unwrap();
            let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                panic!("an edge without a curve");
            };
            let curve = model.geometry().curve(*curve).unwrap().clone();
            out.push((curve, *range, data.tolerance.get()));
        }
    }
    out
}

/// The farthest a curve strays over its range, by a measure of its points.
fn worst_along(curve: &ogeom::geom::Curve, range: (f64, f64), off: impl Fn(Point) -> f64) -> f64 {
    use ogeom::geom::Curve3d as _;
    (0..=2000)
        .map(|k| {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 2000.0;
            off(curve.point_at(t, T).unwrap())
        })
        .fold(0.0, f64::max)
}

/// A round bar of radius five along x, twenty long, with a hill on its
/// side over a quarter turn and the middle twelve of its length: the radius
/// at angle `a` and abscissa `x`. The hill runs out into the cylinder with
/// its height and slope.
fn drum_radius(a: f64, x: f64) -> f64 {
    let (lo, hi) = (0.75 * core::f64::consts::PI, 1.25 * core::f64::consts::PI);
    5.0 + 0.5 * unit_hill((x - 4.0) / 12.0, (a - lo) / (hi - lo))
}

/// The bar of [`drum_radius`] meshed on `turn` cells round and `along`
/// cells along it, lines of the grid on the hill's foot; each end a fan
/// about its centre.
fn drum_mesh(turn: u32, along: u32) -> Triangulation {
    let mut mesh = Triangulation::new();
    for j in 0..=along {
        let x = 20.0 * f64::from(j) / f64::from(along);
        for k in 0..turn {
            let a = core::f64::consts::TAU * f64::from(k) / f64::from(turn);
            let r = drum_radius(a, x);
            mesh.positions.push(Point::new(x, r * a.cos(), r * a.sin()));
        }
    }
    let ring = |j: u32, k: u32| j * turn + k % turn;
    for j in 0..along {
        for k in 0..turn {
            let (a, b, c, d) = (
                ring(j, k),
                ring(j + 1, k),
                ring(j + 1, k + 1),
                ring(j, k + 1),
            );
            mesh.triangles.push([a, c, b]);
            mesh.triangles.push([a, d, c]);
        }
    }
    let base = u32::try_from(mesh.positions.len()).unwrap();
    mesh.positions.push(Point::ORIGIN);
    mesh.positions.push(Point::new(20.0, 0.0, 0.0));
    for k in 0..turn {
        mesh.triangles.push([base, ring(0, k + 1), ring(0, k)]);
        mesh.triangles
            .push([base + 1, ring(along, k), ring(along, k + 1)]);
    }
    mesh
}

/// A hill on a round bar's side, running out into the cylinder
/// tangentially, meshed with vertices along the line where they meet: the
/// cylinder bounds the hill's smooth region, which comes back one patch.
/// The section solve finds nothing at a tangency, and the seam is threaded
/// through the chain in the cylinder's chart, cornered where the chain
/// turns, and lifted onto it: it lies on the cylinder, within twice the
/// coplanar distance of the patch. The volume is within the distance over
/// the hill's area, and the solid valid and tessellating closed.
#[test]
fn a_hill_tangent_to_a_cylinder_meets_it_on_the_cylinder() {
    let mesh = drum_mesh(96, 40);
    let mut model = Model::new();
    let out = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    eprintln!("drum: {} faces, report {:?}", out.report.faces, out.report);
    assert_eq!(out.report.windings_flipped, 0);
    let diagnosis = check(&model, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    assert_eq!(kinds_and_patches(&model, &out.shape), ([2, 1, 0, 0, 0], 1));
    let flat = out.coplanar_distance;
    let [(seam, range, tolerance)] = &patch_edges(&model, &out.shape)[..] else {
        panic!("the patch is bounded by more than its seam");
    };
    let off = worst_along(seam, *range, |p| (p.y.hypot(p.z) - 5.0).abs());
    let volume = volume_properties(&model, &out.shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass;
    // The bar, and the hill: (r^2 - 25) / 2 over its angle and length, by
    // the midpoint rule.
    let (n, pi) = (400, core::f64::consts::PI);
    let mut hill = 0.0;
    for i in 0..n {
        for k in 0..n {
            let x = 4.0 + 12.0 * (f64::from(i) + 0.5) / f64::from(n);
            let a = 0.75 * pi + 0.5 * pi * (f64::from(k) + 0.5) / f64::from(n);
            hill += 0.5 * (drum_radius(a, x).powi(2) - 25.0);
        }
    }
    hill *= 12.0 * 0.5 * pi / f64::from(n * n);
    let exact = pi * 25.0 * 20.0 + hill;
    eprintln!(
        "drum: seam tolerance {tolerance:e} against {flat:e}, off the cylinder {off:e}; volume {volume}, exact {exact}"
    );
    assert!(off <= 1e-8, "{off}");
    assert!(*tolerance <= flat * 2.0, "{tolerance} past {flat}");
    assert!((volume - exact).abs() <= flat * 12.0 * 2.5 * pi);
    let closed = ogeom::mesh::triangulate(&model, &out.shape, Deflection::default(), T).unwrap();
    assert!(closed.is_closed());
}

/// The top of a block twenty long along x, fifteen deep along y and eight
/// high, its front top edge rounded at radius three: past the round it
/// rises from the round's tangent line as a free-form face, its height and
/// slope across the line those of the round's top.
fn rounded_top(x: f64, y: f64) -> f64 {
    let s = (y - 3.0) / 12.0;
    let along = x / 20.0;
    8.0 + 2.0 * s * s * (1.0 + 0.4 * (core::f64::consts::PI * along).sin() + 0.3 * s * along)
}

/// The block of [`rounded_top`] meshed on `along` cells in x, `arc` round
/// the round and `top` across the free-form face: every cross-section's
/// vertices at the same depths, so the walls and the bottom are grids.
fn rounded_block_mesh(along: u32, arc: u32, top: u32) -> Triangulation {
    let quarter = |k: u32| core::f64::consts::FRAC_PI_2 * f64::from(k) / f64::from(arc);
    let mut depths: Vec<f64> = (0..=arc).map(|k| 3.0 - 3.0 * quarter(k).cos()).collect();
    depths.extend((1..=top).map(|k| 3.0 + 12.0 * f64::from(k) / f64::from(top)));
    let height = |x: f64, k: usize| match u32::try_from(k).unwrap() {
        k if k <= arc => 5.0 + 3.0 * quarter(k).sin(),
        _ => rounded_top(x, depths[k]),
    };
    let (m, n) = (depths.len(), along as usize + 1);
    let at = |i: usize| 20.0 * f64::from(u32::try_from(i).unwrap()) / f64::from(along);
    let mut mesh = Triangulation::new();
    for i in 0..n {
        for (k, &y) in depths.iter().enumerate() {
            mesh.positions.push(Point::new(at(i), y, height(at(i), k)));
        }
    }
    for i in 0..n {
        for &y in &depths {
            mesh.positions.push(Point::new(at(i), y, 0.0));
        }
    }
    let up = |i: usize, k: usize| u32::try_from(i * m + k).unwrap();
    let down = |i: usize, k: usize| u32::try_from(n * m + i * m + k).unwrap();
    for i in 0..n - 1 {
        for k in 0..m - 1 {
            mesh.triangles
                .push([up(i, k), up(i + 1, k), up(i + 1, k + 1)]);
            mesh.triangles
                .push([up(i, k), up(i + 1, k + 1), up(i, k + 1)]);
            mesh.triangles
                .push([down(i, k), down(i + 1, k + 1), down(i + 1, k)]);
            mesh.triangles
                .push([down(i, k), down(i, k + 1), down(i + 1, k + 1)]);
        }
        // The front wall and the back.
        let (f, b) = (0, m - 1);
        mesh.triangles
            .push([down(i, f), down(i + 1, f), up(i + 1, f)]);
        mesh.triangles.push([down(i, f), up(i + 1, f), up(i, f)]);
        mesh.triangles
            .push([down(i, b), up(i + 1, b), down(i + 1, b)]);
        mesh.triangles.push([down(i, b), up(i, b), up(i + 1, b)]);
    }
    // The two ends, a column of the section between each two depths.
    let (e, l) = (0, n - 1);
    for k in 0..m - 1 {
        mesh.triangles
            .push([down(e, k), up(e, k + 1), down(e, k + 1)]);
        mesh.triangles.push([down(e, k), up(e, k), up(e, k + 1)]);
        mesh.triangles
            .push([down(l, k), down(l, k + 1), up(l, k + 1)]);
        mesh.triangles.push([down(l, k), up(l, k + 1), up(l, k)]);
    }
    mesh
}

/// A round along a block's edge running tangentially into a free-form top,
/// meshed with vertices along the line where they meet: the round bounds
/// the top's smooth region, which comes back one patch, and the seam
/// between them is the round's ruling through the chain, exactly on it and
/// within the coplanar distance of the patch. The volume is within the
/// distance over the top's area, and the solid valid and tessellating
/// closed.
#[test]
fn a_round_running_into_a_free_form_top_meets_it_on_a_ruling() {
    let mesh = rounded_block_mesh(40, 12, 24);
    let mut model = Model::new();
    let out = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    eprintln!("round: {} faces, report {:?}", out.report.faces, out.report);
    assert_eq!(out.report.windings_flipped, 0);
    let diagnosis = check(&model, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    assert_eq!(kinds_and_patches(&model, &out.shape), ([5, 1, 0, 0, 0], 1));
    let flat = out.coplanar_distance;
    let round = |p: Point| ((p.y - 3.0).hypot(p.z - 5.0) - 3.0).abs();
    let on_round: Vec<_> = patch_edges(&model, &out.shape)
        .into_iter()
        .filter(|(curve, range, _)| worst_along(curve, *range, round) <= 1e-9)
        .collect();
    let [(seam, _, tolerance)] = &on_round[..] else {
        panic!("{} seams on the round", on_round.len());
    };
    assert!(matches!(seam, ogeom::geom::Curve::Line(_)));
    assert!(*tolerance <= flat, "{tolerance} past {flat}");
    let volume = volume_properties(&model, &out.shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass;
    // The section's area along x (the wall under the round, the round's
    // quarter disc, the block under the top, and the top's rise over
    // eight) integrated in closed form.
    let pi = core::f64::consts::PI;
    let exact = 20.0 * (5.0 * 3.0 + pi * 9.0 / 4.0 + 8.0 * 12.0)
        + 2.0 * 12.0 * (20.0 / 3.0 + 0.4 * 40.0 / (3.0 * pi) + 0.3 * 20.0 / 8.0);
    eprintln!(
        "round: seam tolerance {tolerance:e} against {flat:e}; volume {volume}, exact {exact}"
    );
    assert!((volume - exact).abs() <= flat * 20.0 * 12.0);
    let closed = ogeom::mesh::triangulate(&model, &out.shape, Deflection::default(), T).unwrap();
    assert!(closed.is_closed());
}

/// A ball of radius ten with a hill on its side between polar angles of an
/// eighth and three eighths of a half turn, over a quarter turn about its
/// axis: the radius at polar angle `p` and azimuth `a`.
fn ball_radius(p: f64, a: f64) -> f64 {
    let pi = core::f64::consts::PI;
    10.0 + 0.5 * unit_hill((p - pi / 8.0) / (pi / 4.0), a / (pi / 2.0))
}

/// The ball of [`ball_radius`] meshed on `rings` cells from pole to pole
/// and `turn` round, a fan at either pole.
fn ball_mesh(rings: u32, turn: u32) -> Triangulation {
    let pi = core::f64::consts::PI;
    let mut mesh = Triangulation::new();
    mesh.positions.push(Point::new(0.0, 0.0, 10.0));
    for j in 1..rings {
        let p = pi * f64::from(j) / f64::from(rings);
        for k in 0..turn {
            let a = 2.0 * pi * f64::from(k) / f64::from(turn);
            let r = ball_radius(p, a);
            mesh.positions.push(Point::new(
                r * p.sin() * a.cos(),
                r * p.sin() * a.sin(),
                r * p.cos(),
            ));
        }
    }
    let at = |j: u32, k: u32| 1 + (j - 1) * turn + k % turn;
    let south = u32::try_from(mesh.positions.len()).unwrap();
    mesh.positions.push(Point::new(0.0, 0.0, -10.0));
    for k in 0..turn {
        mesh.triangles.push([0, at(1, k), at(1, k + 1)]);
        mesh.triangles
            .push([south, at(rings - 1, k + 1), at(rings - 1, k)]);
    }
    for j in 1..rings - 1 {
        for k in 0..turn {
            let (a, b, c, d) = (at(j, k), at(j + 1, k), at(j + 1, k + 1), at(j, k + 1));
            mesh.triangles.push([a, b, c]);
            mesh.triangles.push([a, c, d]);
        }
    }
    mesh
}

/// A hill on a ball's side, running out into the sphere tangentially: the
/// sphere, far larger than the hill, bounds it rather than joining it, and
/// the seam threaded in the sphere's chart lies on the sphere. The cylinder
/// case above is the quick one; a ball coarse enough to be quick leaves
/// the hill too few vertices for a patch.
#[test]
#[ignore = "heavy"]
fn a_hill_tangent_to_a_sphere_meets_it_on_the_sphere() {
    let mesh = ball_mesh(96, 128);
    let mut model = Model::new();
    let out = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    eprintln!("ball: {} faces, report {:?}", out.report.faces, out.report);
    assert_eq!(out.report.windings_flipped, 0);
    let diagnosis = check(&model, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    assert_eq!(kinds_and_patches(&model, &out.shape), ([0, 0, 0, 1, 0], 1));
    let flat = out.coplanar_distance;
    let [(seam, range, tolerance)] = &patch_edges(&model, &out.shape)[..] else {
        panic!("the patch is bounded by more than its seam");
    };
    let off = worst_along(seam, *range, |p| (p.to_vector().magnitude() - 10.0).abs());
    let volume = volume_properties(&model, &out.shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass;
    // The ball, and the hill: (r^3 - 1000) / 3 over its solid angle, by the
    // midpoint rule.
    let (n, pi) = (400, core::f64::consts::PI);
    let mut hill = 0.0;
    for i in 0..n {
        for k in 0..n {
            let p = pi / 8.0 + pi / 4.0 * (f64::from(i) + 0.5) / f64::from(n);
            let a = pi / 2.0 * (f64::from(k) + 0.5) / f64::from(n);
            hill += (ball_radius(p, a).powi(3) - 1000.0) / 3.0 * p.sin();
        }
    }
    hill *= pi / 4.0 * pi / 2.0 / f64::from(n * n);
    let exact = 4.0 / 3.0 * pi * 1000.0 + hill;
    eprintln!(
        "ball: seam tolerance {tolerance:e} against {flat:e}, off the sphere {off:e}; volume {volume}, exact {exact}"
    );
    assert!(off <= 1e-6, "{off}");
    assert!(*tolerance <= flat * 3.0, "{tolerance} past {flat}");
    assert!((volume - exact).abs() <= flat * 100.0 * pi * pi / 8.0);
    let closed = ogeom::mesh::triangulate(&model, &out.shape, Deflection::default(), T).unwrap();
    assert!(closed.is_closed());
}

/// The mesh of a shape's faces on one kind of surface alone: the shape
/// with every other face taken away, open where they were.
fn faces_meshed(
    model: &Model,
    shape: &Shape,
    keep: fn(&ogeom::geom::SurfaceGeometry) -> bool,
) -> Triangulation {
    let mut out = Triangulation::new();
    for face in explore_unique(model, shape, ShapeType::Face).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        if !keep(model.geometry().surface(data.surface).unwrap()) {
            continue;
        }
        let mesh = meshed(model, &face);
        let base = u32::try_from(out.positions.len()).unwrap();
        out.positions.extend(&mesh.positions);
        out.triangles
            .extend(mesh.triangles.iter().map(|t| t.map(|v| v + base)));
    }
    out
}

/// The segments of a mesh that one triangle uses, its vertices welded at a
/// nanometre.
fn mesh_boundary(mesh: &Triangulation) -> Vec<[Point; 2]> {
    #[allow(clippy::cast_possible_truncation, reason = "a grid cell")]
    let key = |p: Point| {
        (
            (p.x * 1e6).round() as i64,
            (p.y * 1e6).round() as i64,
            (p.z * 1e6).round() as i64,
        )
    };
    let mut uses: std::collections::HashMap<_, (usize, [Point; 2])> =
        std::collections::HashMap::new();
    for t in &mesh.triangles {
        for k in 0..3 {
            let (p, q) = (
                mesh.positions[t[k] as usize],
                mesh.positions[t[(k + 1) % 3] as usize],
            );
            let (a, b) = (key(p), key(q));
            if a == b {
                continue;
            }
            uses.entry((a.min(b), a.max(b))).or_insert((0, [p, q])).0 += 1;
        }
    }
    uses.into_values()
        .filter(|(n, _)| *n == 1)
        .map(|(_, s)| s)
        .collect()
}

/// The free edges of a converted shell (bounding one face, and no seam
/// of it): each one's curve trimmed to its range, and its tolerance.
fn free_edges(model: &Model, shape: &Shape) -> Vec<(ogeom::geom::Curve, f64)> {
    let mut count: std::collections::HashMap<ogeom::topo::SameKey, usize> =
        std::collections::HashMap::new();
    for face in explore_unique(model, shape, ShapeType::Face).unwrap() {
        for edge in ogeom::topo::explore(model, &face, ogeom::topo::Filter::OfType(ShapeType::Edge))
            .unwrap()
        {
            *count.entry(ogeom::topo::SameKey(edge)).or_default() += 1;
        }
    }
    let mut out = Vec::new();
    for (key, n) in count {
        let data = model.node(&key.0).unwrap().data().as_edge().unwrap();
        if n != 1 || data.degenerate {
            continue;
        }
        let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            panic!("an edge without a curve");
        };
        let curve = model.geometry().curve(*curve).unwrap().clone();
        let trimmed = ogeom::geom::TrimmedCurve::new(curve, range.0, range.1, T).unwrap();
        out.push((
            ogeom::geom::Curve::Trimmed(Box::new(trimmed)),
            data.tolerance.get(),
        ));
    }
    out
}

/// The curve under a trimmed free edge.
fn basis(curve: &ogeom::geom::Curve) -> &ogeom::geom::Curve {
    match curve {
        ogeom::geom::Curve::Trimmed(t) => t.basis(),
        other => other,
    }
}

/// A converted open mesh's free boundary against the mesh's: every vertex
/// of the mesh's boundary within its edge's tolerance of a free edge, and
/// every point of the converted shell's own boundary, as tessellated,
/// within `sag` (what the mesh's chords cut off its boundary) and the
/// edges' tolerance of the mesh's boundary. The worst of each, over the
/// allowance it was held to.
fn holds_the_mesh_boundary(
    mesh: &Triangulation,
    model: &Model,
    shape: &Shape,
    sag: f64,
) -> (f64, f64) {
    let boundary = mesh_boundary(mesh);
    let edges = free_edges(model, shape);
    let widest = edges.iter().map(|e| e.1).fold(0.0, f64::max);
    let mut on_edges: f64 = 0.0;
    for p in boundary.iter().flatten() {
        let (off, tolerance) = edges
            .iter()
            .map(|(curve, tolerance)| {
                let near = ogeom::algo::project_on_curve(curve, *p, 256, T).unwrap();
                (near.distance, *tolerance)
            })
            .min_by(|a, b| a.0.total_cmp(&b.0))
            .unwrap();
        assert!(off <= tolerance, "{p:?} {off:e} off, past {tolerance:e}");
        on_edges = on_edges.max(off / tolerance);
    }
    let drawn = ogeom::mesh::triangulate(model, shape, Deflection::default(), T).unwrap();
    let mut on_mesh: f64 = 0.0;
    for p in mesh_boundary(&drawn).iter().flatten() {
        let off = boundary
            .iter()
            .map(|s| {
                let (d, along) = (s[1] - s[0], *p - s[0]);
                let f = (along.dot(d) / d.dot(d)).clamp(0.0, 1.0);
                p.distance(s[0] + d * f)
            })
            .fold(f64::INFINITY, f64::min);
        assert!(off <= sag + widest, "{p:?} {off:e} off the mesh's boundary");
        on_mesh = on_mesh.max(off / (sag + widest));
    }
    (on_edges, on_mesh)
}

/// Converts an open mesh, checks the shell usable and open, and gives the
/// model, the conversion and the shell's free edges.
fn converted_open(
    mesh: &Triangulation,
) -> (
    Model,
    ogeom::algo::MeshSolid,
    Vec<(ogeom::geom::Curve, f64)>,
) {
    let mut model = Model::new();
    let out = solid_from_mesh(&mut model, mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(!out.closed);
    assert!(out.report.edges_used_once > 0);
    assert_eq!(model.kind_of(&out.shape).unwrap(), ShapeType::Shell);
    let diagnosis = check(&model, &out.shape, T).unwrap();
    assert!(diagnosis.is_usable(), "{diagnosis}");
    let edges = free_edges(&model, &out.shape);
    (model, out, edges)
}

/// A cylinder's side meshed alone, both caps taken away: a tube open at
/// both ends. It comes back one cylinder face bounded by two exact circles
/// of the cylinder's radius at its two ends, and the converted boundary
/// keeps to the mesh's.
#[test]
fn an_open_tube_comes_back_one_cylinder_between_two_circles() {
    let mut model = Model::new();
    let rod = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 4.0, 10.0, T)
        .unwrap()
        .shape;
    let mesh = faces_meshed(&model, &rod, |s| {
        matches!(s, ogeom::geom::SurfaceGeometry::Cylinder(_))
    });
    let (back, out, edges) = converted_open(&mesh);
    assert_eq!(
        kinds(&back, &out.shape),
        [0, 1, 0, 0, 0],
        "{:?}",
        out.report
    );
    assert_eq!(out.report.free_edges_fitted, 0);
    assert_eq!(edges.len(), 2);
    let mut heights = Vec::new();
    for (curve, tolerance) in &edges {
        let ogeom::geom::Curve::Circle(c) = basis(curve) else {
            panic!("a rim that is no circle: {curve:?}");
        };
        let circle = c.circle();
        eprintln!(
            "tube: rim radius {} at {:?}, tolerance {tolerance:e}",
            circle.radius(),
            circle.centre()
        );
        assert!((circle.radius() - 4.0).abs() <= *tolerance);
        let centre = circle.centre();
        assert!(Vector::new(centre.x, centre.y, 0.0).magnitude() <= *tolerance);
        heights.push(centre.z);
    }
    heights.sort_by(f64::total_cmp);
    assert!(
        heights[0].abs() <= 1e-6 && (heights[1] - 10.0).abs() <= 1e-6,
        "{heights:?}"
    );
    let held = holds_the_mesh_boundary(&mesh, &back, &out.shape, 0.01);
    eprintln!("tube: boundary held to {held:?} of its allowances");
}

/// Half a cylinder's side, meshed alone: a sheet open all round. It comes
/// back one cylinder face, its boundary cut where it turns square into
/// four edges: two rulings along the axis at the cut, and two half circles
/// at the ends.
#[test]
fn a_half_cylinder_sheet_comes_back_bounded_by_rulings_and_arcs() {
    use ogeom::geom::Curve3d as _;
    let mut model = Model::new();
    let rod = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 4.0, 10.0, T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(-5.0, 0.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let half_box = ogeom::algo::make_box(&mut model, at, (10.0, 5.0, 12.0), T)
        .unwrap()
        .shape;
    let half = ogeom::boolean::common(&mut model, &rod, &half_box, T)
        .unwrap()
        .shape;
    let mesh = faces_meshed(&model, &half, |s| {
        matches!(s, ogeom::geom::SurfaceGeometry::Cylinder(_))
    });
    let (back, out, edges) = converted_open(&mesh);
    assert_eq!(
        kinds(&back, &out.shape),
        [0, 1, 0, 0, 0],
        "{:?}",
        out.report
    );
    assert_eq!(out.report.free_edges_fitted, 0);
    assert_eq!(edges.len(), 4);
    let (mut rulings, mut arcs) = (0, 0);
    for (curve, tolerance) in &edges {
        let (a, b) = curve.domain();
        let (p, q) = (curve.point_at(a, T).unwrap(), curve.point_at(b, T).unwrap());
        match basis(curve) {
            ogeom::geom::Curve::Line(line) => {
                let along = line.axis().direction.vector();
                eprintln!("half: ruling from {p:?} to {q:?}, tolerance {tolerance:e}");
                assert!(along.cross(Vector::Z).magnitude() <= 1e-9);
                assert!(p.y.abs() <= *tolerance && (p.x.abs() - 4.0).abs() <= *tolerance);
                assert!((p.distance(q) - 10.0).abs() <= 1e-6);
                rulings += 1;
            }
            ogeom::geom::Curve::Circle(c) => {
                let circle = c.circle();
                eprintln!(
                    "half: arc of radius {} at {:?} over {}, tolerance {tolerance:e}",
                    circle.radius(),
                    circle.centre(),
                    b - a
                );
                assert!((circle.radius() - 4.0).abs() <= *tolerance);
                assert!(((b - a) - core::f64::consts::PI).abs() <= 1e-6);
                assert!(p.y.abs() <= *tolerance && q.y.abs() <= *tolerance);
                arcs += 1;
            }
            other => panic!("a free edge that is no ruling or arc: {other:?}"),
        }
    }
    assert_eq!((rulings, arcs), (2, 2));
    let held = holds_the_mesh_boundary(&mesh, &back, &out.shape, 0.01);
    eprintln!("half: boundary held to {held:?} of its allowances");
}

/// A ball's top above a plane, meshed alone: a dome open along its rim. It
/// comes back one sphere face whose rim is the latitude circle the plane
/// cut.
#[test]
fn a_dome_s_rim_comes_back_a_circle_on_its_sphere() {
    let mut model = Model::new();
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 7.0, T)
        .unwrap()
        .shape;
    let above = Frame::new(Point::new(-8.0, -8.0, 3.0), Direction::Z, Direction::X, T).unwrap();
    let top_box = ogeom::algo::make_box(&mut model, above, (16.0, 16.0, 8.0), T)
        .unwrap()
        .shape;
    let dome = ogeom::boolean::common(&mut model, &ball, &top_box, T)
        .unwrap()
        .shape;
    let mesh = faces_meshed(&model, &dome, |s| {
        matches!(s, ogeom::geom::SurfaceGeometry::Sphere(_))
    });
    let (back, out, edges) = converted_open(&mesh);
    assert_eq!(
        kinds(&back, &out.shape),
        [0, 0, 0, 1, 0],
        "{:?}",
        out.report
    );
    assert_eq!(out.report.free_edges_fitted, 0);
    assert_eq!(edges.len(), 1);
    let (curve, tolerance) = &edges[0];
    let ogeom::geom::Curve::Circle(c) = basis(curve) else {
        panic!("a rim that is no circle: {curve:?}");
    };
    let circle = c.circle();
    eprintln!(
        "dome: rim radius {} at {:?}, tolerance {tolerance:e}",
        circle.radius(),
        circle.centre()
    );
    assert!((circle.radius() - 40.0_f64.sqrt()).abs() <= *tolerance);
    assert!(circle.centre().distance(Point::new(0.0, 0.0, 3.0)) <= *tolerance);
    assert!(circle.frame().z().vector().cross(Vector::Z).magnitude() <= 1e-9);
    let held = holds_the_mesh_boundary(&mesh, &back, &out.shape, 0.01);
    eprintln!("dome: boundary held to {held:?} of its allowances");
}

/// The bump's height with a bowl added, so the sheet curves everywhere and
/// no two of its triangles read as one flat face.
fn bowl(x: f64, y: f64) -> f64 {
    bump(x, y) + 0.05 * ((x - 10.0).powi(2) + (y - 10.0).powi(2))
}

/// A free-form sheet meshed alone, open all round its square: it comes back
/// one fitted patch whose border, cut at its four corners, is four fitted
/// curves. Each holds every vertex of the mesh's border within the
/// tolerance it states, and the shell's own border, as tessellated, lies
/// within that tolerance of the exact surface's border.
#[test]
#[ignore = "heavy"]
fn a_free_form_sheet_s_border_comes_back_four_fitted_curves() {
    let full = grid_solid(40, 20.0, |_, _| true, bowl);
    let mut mesh = Triangulation::new();
    mesh.positions.clone_from(&full.positions);
    mesh.triangles = full
        .triangles
        .iter()
        .copied()
        .filter(|t| t.iter().all(|&v| v < 41 * 41))
        .collect();
    let (back, out, edges) = converted_open(&mesh);
    eprintln!(
        "sheet: distance {:e}, {:?}",
        out.coplanar_distance, out.report
    );
    assert_eq!(out.report.faces, 1);
    assert_eq!(out.report.patch_faces, 1);
    assert_eq!(out.report.free_edges_fitted, 4);
    assert_eq!(edges.len(), 4);
    let widest = edges.iter().map(|e| e.1).fold(0.0, f64::max);
    for (curve, tolerance) in &edges {
        assert!(
            matches!(basis(curve), ogeom::geom::Curve::BSpline(_)),
            "{curve:?}"
        );
        eprintln!("sheet: border edge tolerance {tolerance:e}");
        assert!(*tolerance <= out.coplanar_distance * 20.0);
    }
    let held = holds_the_mesh_boundary(&mesh, &back, &out.shape, 0.01);
    eprintln!("sheet: boundary held to {held:?} of its allowances");
    // The shell's border, as drawn, against the exact border: the side of
    // the square nearest each point, and the height over it.
    let drawn = ogeom::mesh::triangulate(&back, &out.shape, Deflection::default(), T).unwrap();
    let mut worst: f64 = 0.0;
    for p in mesh_boundary(&drawn).iter().flatten() {
        let (x, y) = (p.x.clamp(0.0, 20.0), p.y.clamp(0.0, 20.0));
        let off = [
            ((0.0, y), p.x.abs()),
            ((20.0, y), (p.x - 20.0).abs()),
            ((x, 0.0), p.y.abs()),
            ((x, 20.0), (p.y - 20.0).abs()),
        ]
        .into_iter()
        .map(|(on, across)| across.hypot(p.z - bowl(on.0, on.1)))
        .fold(f64::INFINITY, f64::min);
        worst = worst.max(off);
    }
    eprintln!("sheet: drawn border {worst:e} off the exact border, tolerance {widest:e}");
    assert!(worst <= widest, "{worst} past {widest}");
}

/// The bump alone, meshed open all round its square, with no flat stretch
/// anywhere: its few nearly flat triangle pairs turn steadily, row after
/// row, and are not read as scatter. The coplanar distance stays near the
/// closed plate's, and the sheet comes back one fitted patch.
#[test]
#[ignore = "heavy"]
fn a_bumped_sheet_alone_comes_back_one_fitted_patch() {
    let full = grid_solid(40, 20.0, |_, _| true, bump);
    let mut mesh = Triangulation::new();
    mesh.positions.clone_from(&full.positions);
    mesh.triangles = full
        .triangles
        .iter()
        .copied()
        .filter(|t| t.iter().all(|&v| v < 41 * 41))
        .collect();
    let (_, out, edges) = converted_open(&mesh);
    eprintln!(
        "bumped sheet: distance {:e}, {:?}",
        out.coplanar_distance, out.report
    );
    assert!(out.coplanar_distance <= 2e-4, "{}", out.coplanar_distance);
    assert_eq!(out.report.faces, 1);
    assert_eq!(out.report.patch_faces, 1);
    assert_eq!(edges.len(), 4);
}

/// The bump's plate over the cells `cells` each way of an `n` by `n` grid
/// on its square, meshed finely enough that small spheres, cones and
/// cylinders fit stretches of its top within the coplanar distance. Its
/// top comes back one fitted B-spline face all the same: the solid valid
/// with six faces, the exact surface within the coplanar distance of the
/// patch at points off the mesh's grid, the volume a hundred times nearer
/// the exact one than the mesh's, and tessellating closed.
fn fine_bump_comes_back_one_patch(n: u32, cells: core::ops::Range<u32>, samples: u32) {
    let mesh = grid_solid(
        n,
        20.0,
        |i, j| cells.contains(&i) && cells.contains(&j),
        bump,
    );
    let mut model = Model::new();
    let started = Instant::now();
    let out = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    let took = started.elapsed();
    let diagnosis = check(&model, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let patches = spline_faces(&model, &out.shape);
    eprintln!(
        "fine bump {n}: {} faces, {} patches in {took:?}, distance {:e}, report {:?}",
        out.report.faces,
        patches.len(),
        out.coplanar_distance,
        out.report
    );
    assert_eq!(patches.len(), 1);
    assert_eq!(out.report.patch_faces, 1);
    assert_eq!(out.report.faces, 6);
    let data = model.node(&patches[0]).unwrap().data().as_face().unwrap();
    let surface = model.geometry().surface(data.surface).unwrap().clone();
    let flat = out.coplanar_distance;
    let at = |k: u32| 20.0 * f64::from(k) / f64::from(n);
    let (low, high) = (at(cells.start), at(cells.end));
    let mut worst: f64 = 0.0;
    for j in 0..samples {
        for i in 0..samples {
            let x = low + (high - low) * (f64::from(i) + 0.37) / f64::from(samples);
            let y = low + (high - low) * (f64::from(j) + 0.61) / f64::from(samples);
            let exact = Point::new(x, y, bump(x, y));
            let foot = ogeom::algo::project_on_surface(&surface, exact, 16, T).unwrap();
            worst = worst.max(foot.distance);
        }
    }
    // The exact volume over the cells by Simpson's rule each way.
    let steps = 200_u32;
    let weight = |k: u32| {
        if k == 0 || k == steps {
            1.0
        } else if k % 2 == 1 {
            4.0
        } else {
            2.0
        }
    };
    let h = (high - low) / f64::from(steps);
    let mut exact = 0.0;
    for j in 0..=steps {
        for i in 0..=steps {
            let (x, y) = (low + h * f64::from(i), low + h * f64::from(j));
            exact += weight(i) * weight(j) * bump(x, y);
        }
    }
    exact *= h * h / 9.0;
    let volume = volume_properties(&model, &out.shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass;
    let mesh_volume = mesh.volume();
    eprintln!(
        "fine bump {n}: exact to patch {worst:e}, distance {flat:e}; volume {volume}, the mesh's {mesh_volume}, exact {exact}"
    );
    assert!(worst <= flat, "{worst} past {flat}");
    assert!((volume - exact).abs() * 100.0 < (mesh_volume - exact).abs());
    let closed = ogeom::mesh::triangulate(&model, &out.shape, Deflection::default(), T).unwrap();
    assert!(closed.is_closed());
}

/// A quarter of the bump's width meshed at a hundred and twenty cells
/// across the whole square: recognition takes small spheres out of its top
/// at every distance it tries, and the free triangles left round them are
/// no disk.
#[test]
fn a_finely_meshed_bump_patch_takes_the_small_spheres_cut_from_it() {
    fine_bump_comes_back_one_patch(120, 45..75, 12);
}

/// The whole bumped plate at eighty cells across: its top is cut into a
/// cylinder, two spheres and a torus among the free triangles.
#[test]
#[ignore = "heavy"]
fn a_finely_meshed_bump_comes_back_one_fitted_patch() {
    fine_bump_comes_back_one_patch(80, 0..80, 40);
}

/// A dome of radius 10 on its base: on a drum of the same radius 10
/// high, or flat on a plate (a hemisphere). Drilled straight down by a
/// cylinder of `radius` whose axis stands at `(x, y)`, where given.
fn dome(model: &mut Model, on_drum: bool, drill: Option<((f64, f64), f64)>) -> Shape {
    let base = if on_drum { 10.0 } else { 0.0 };
    let top = Frame::new(Point::new(0.0, 0.0, base), Direction::Z, Direction::X, T).unwrap();
    let ball = ogeom::algo::make_sphere(model, top, 10.0, T).unwrap().shape;
    let dome = if on_drum {
        let drum = ogeom::algo::make_cylinder(model, Frame::WORLD, 10.0, 10.0, T)
            .unwrap()
            .shape;
        ogeom::boolean::fuse(model, &drum, &ball, T).unwrap().shape
    } else {
        let above =
            Frame::new(Point::new(-20.0, -20.0, 0.0), Direction::Z, Direction::X, T).unwrap();
        let half = ogeom::algo::make_box(model, above, (40.0, 40.0, 20.0), T)
            .unwrap()
            .shape;
        ogeom::boolean::common(model, &ball, &half, T)
            .unwrap()
            .shape
    };
    let Some(((x, y), radius)) = drill else {
        return dome;
    };
    let below = Frame::new(Point::new(x, y, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(model, below, radius, 30.0, T)
        .unwrap()
        .shape;
    ogeom::boolean::cut(model, &dome, &drill, T).unwrap().shape
}

/// Converts a drilled [`dome`]'s mesh and checks it comes back on
/// `expected` surfaces with no curved region faceted, valid and meshing
/// closed, and with the volume, measured on the exact surfaces both sides,
/// within a millionth of the original's.
fn drilled_dome_comes_back(name: &str, on_drum: bool, at: (f64, f64), expected: [usize; 5]) {
    let mut model = Model::new();
    let drilled = dome(&mut model, on_drum, Some((at, 1.5)));
    let mesh = ogeom::mesh::triangulate(&model, &drilled, Deflection::with_chord(0.05).unwrap(), T)
        .unwrap();
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(out.closed, "{name}: {:?}", out.report);
    assert_eq!(out.report.curved_faceted, 0, "{name}: {:?}", out.report);
    assert_eq!(kinds(&back, &out.shape), expected, "{name}");
    let diagnosis = check(&back, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
    assert_eq!(edges_walked_one_way(&back, &out.shape), 0, "{name}");
    let drawn = ogeom::mesh::triangulate(&back, &out.shape, Deflection::default(), T).unwrap();
    assert!(drawn.is_closed(), "{name}");
    let (a, b) = (volume(&model, &drilled), volume(&back, &out.shape));
    eprintln!("{name}: {a} went in, {b} came out");
    assert!(
        (a - b).abs() / a < 1e-6,
        "{name}: {a} went in, {b} came out"
    );
}

/// A hemisphere drilled off its pole is a sphere's cap with a hole that
/// is no latitude: the sphere comes back one face about a pole clear of
/// the hole, its seam from the rim to the pole clear of it too, and the
/// hole an inner wire. Drilled through its pole off its axis, the sphere
/// is a zone between the rim and the hole, each going round an axis
/// through the hole, and comes back one face with a seam between them.
#[test]
fn hemispheres_drilled_off_and_through_the_pole_keep_their_sphere_one_face() {
    drilled_dome_comes_back("off the pole", false, (4.0, 1.0), [1, 1, 0, 1, 0]);
    drilled_dome_comes_back("through the pole", false, (0.8, 0.3), [1, 1, 0, 1, 0]);
}

/// The same on a dome standing on a drum, whose axis the sphere shares.
#[test]
fn domes_on_a_drum_drilled_off_and_through_the_pole_keep_their_sphere_one_face() {
    drilled_dome_comes_back("off the pole", true, (4.0, 1.0), [1, 2, 0, 1, 0]);
    drilled_dome_comes_back("through the pole", true, (0.8, 0.3), [1, 2, 0, 1, 0]);
}

/// A hemisphere of radius 10 on its flat base, with the hill of
/// [`ball_radius`] on its side, meshed on `rings` cells from the pole to
/// the base's rim and `turn` round, the base a fan about its centre.
fn hill_dome_mesh(rings: u32, turn: u32) -> Triangulation {
    let pi = core::f64::consts::PI;
    let mut mesh = Triangulation::new();
    mesh.positions.push(Point::new(0.0, 0.0, 10.0));
    for j in 1..=rings {
        let p = pi / 2.0 * f64::from(j) / f64::from(rings);
        for k in 0..turn {
            let a = 2.0 * pi * f64::from(k) / f64::from(turn);
            let r = ball_radius(p, a);
            mesh.positions.push(Point::new(
                r * p.sin() * a.cos(),
                r * p.sin() * a.sin(),
                if j == rings { 0.0 } else { r * p.cos() },
            ));
        }
    }
    let at = |j: u32, k: u32| 1 + (j - 1) * turn + k % turn;
    let centre = u32::try_from(mesh.positions.len()).unwrap();
    mesh.positions.push(Point::ORIGIN);
    for k in 0..turn {
        mesh.triangles.push([0, at(1, k), at(1, k + 1)]);
        mesh.triangles
            .push([centre, at(rings, k + 1), at(rings, k)]);
    }
    for j in 1..rings {
        for k in 0..turn {
            let (a, b, c, d) = (at(j, k), at(j + 1, k), at(j + 1, k + 1), at(j, k + 1));
            mesh.triangles.push([a, b, c]);
            mesh.triangles.push([a, c, d]);
        }
    }
    mesh
}

/// A hill on a hemisphere's side, running out into the sphere
/// tangentially: the sphere is a cap whose hole is the hill's free-form
/// border, one of whose sides runs along a meridian. It comes back one
/// face, its seam clear of that side, the hill one fitted patch and the
/// base a plane, valid and meshing closed, the volume within the mesh's
/// own sag over the hill.
#[test]
fn a_hill_on_a_hemisphere_keeps_the_sphere_one_face() {
    hill_dome_comes_back(40, 96);
}

/// The same hill meshed finer.
#[test]
#[ignore = "heavy"]
fn a_finely_meshed_hill_on_a_hemisphere_keeps_the_sphere_one_face() {
    hill_dome_comes_back(48, 128);
}

/// [`a_hill_on_a_hemisphere_keeps_the_sphere_one_face`] on `rings` cells
/// from the pole to the base and `turn` round.
fn hill_dome_comes_back(rings: u32, turn: u32) {
    let mesh = hill_dome_mesh(rings, turn);
    let mut model = Model::new();
    let out = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    eprintln!(
        "hill dome: {} faces, report {:?}",
        out.report.faces, out.report
    );
    assert_eq!(out.report.curved_faceted, 0, "{:?}", out.report);
    assert_eq!(kinds_and_patches(&model, &out.shape), ([1, 0, 0, 1, 0], 1));
    let diagnosis = check(&model, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let closed = ogeom::mesh::triangulate(&model, &out.shape, Deflection::default(), T).unwrap();
    assert!(closed.is_closed());
    let flat = out.coplanar_distance;
    let volume = volume_properties(&model, &out.shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass;
    // The hemisphere, and the hill: (r^3 - 1000) / 3 over its solid
    // angle, by the midpoint rule.
    let (n, pi) = (400, core::f64::consts::PI);
    let mut hill = 0.0;
    for i in 0..n {
        for k in 0..n {
            let p = pi / 8.0 + pi / 4.0 * (f64::from(i) + 0.5) / f64::from(n);
            let a = pi / 2.0 * (f64::from(k) + 0.5) / f64::from(n);
            hill += (ball_radius(p, a).powi(3) - 1000.0) / 3.0 * p.sin();
        }
    }
    hill *= pi / 4.0 * pi / 2.0 / f64::from(n * n);
    let exact = 2.0 / 3.0 * pi * 1000.0 + hill;
    let mesh_volume = mesh.volume();
    eprintln!(
        "hill dome: volume {volume}, the mesh's {mesh_volume}, exact {exact}, distance {flat:e}"
    );
    assert!((volume - exact).abs() <= flat * 100.0 * pi * pi / 8.0);
    assert!((volume - exact).abs() < (mesh_volume - exact).abs());
}

/// A boss on a plate: a square with round corners, its walls drafted 2
/// degrees, filleted to the plate at its foot. Its profile (top edge, the
/// wall, the fillet in `fillet_steps` facets, the plate's top and outer
/// wall) is swept along the outline, each corner in `corner_steps`
/// facets: along the sides the wall is a plane and the fillet a cylinder,
/// round the corners a cone and a torus. The top and bottom are fans.
fn drafted_boss(corner_steps: u32, fillet_steps: u32) -> Triangulation {
    let (half, rho, height, fillet) = (10.0_f64, 7.0_f64, 2.2_f64, 0.5_f64);
    let (plate, thick) = (5.0, 3.0);
    let draft = 2.0_f64.to_radians();
    // The wall leans in as it rises: r = rho - z tan(draft). The fillet's
    // centre stands a fillet radius off both the plate and the wall.
    let centre = (rho - fillet * draft.tan() + fillet / draft.cos(), fillet);
    let mut profile = vec![(rho - height * draft.tan(), height)];
    for k in 0..=fillet_steps {
        let a = core::f64::consts::PI
            + draft
            + (core::f64::consts::FRAC_PI_2 - draft) * f64::from(k) / f64::from(fillet_steps);
        profile.push((centre.0 + fillet * a.cos(), centre.1 + fillet * a.sin()));
    }
    profile.push((centre.0 + plate, 0.0));
    profile.push((centre.0 + plate, -thick));
    let mut ring = Vec::new();
    for (q, (cx, cy)) in [
        (0.0, (half, half)),
        (1.0, (-half, half)),
        (2.0, (-half, -half)),
        (3.0, (half, -half)),
    ] {
        for k in 0..=corner_steps {
            let a = (q + f64::from(k) / f64::from(corner_steps)) * core::f64::consts::FRAC_PI_2;
            ring.push(((cx, cy), (a.cos(), a.sin())));
        }
    }
    let (m, rows) = (
        u32::try_from(ring.len()).unwrap(),
        u32::try_from(profile.len()).unwrap(),
    );
    let mut t = Triangulation::new();
    for &(r, z) in &profile {
        for &((cx, cy), (nx, ny)) in &ring {
            t.positions.push(Point::new(cx + nx * r, cy + ny * r, z));
        }
    }
    for k in 0..rows - 1 {
        for i in 0..m {
            let j = (i + 1) % m;
            let (a, b, c, d) = (k * m + i, k * m + j, (k + 1) * m + j, (k + 1) * m + i);
            t.triangles.push([a, c, b]);
            t.triangles.push([a, d, c]);
        }
    }
    let (cap, bottom) = (rows * m, rows * m + 1);
    t.positions.push(Point::new(0.0, 0.0, height));
    t.positions.push(Point::new(0.0, 0.0, -thick));
    let last = (rows - 1) * m;
    for i in 0..m {
        let j = (i + 1) % m;
        t.triangles.push([cap, i, j]);
        t.triangles.push([bottom, last + j, last + i]);
    }
    t
}

/// A drafted boss's corners meshed four facets round and one high, each
/// facet a chord between the walls tangent to the corner: too few
/// vertices for a cone to be fitted to a corner alone, but the walls fix
/// its axis and the rulings between the facets verify it. The corners
/// come back cones, and the fillet at their foot tori seamed to them;
/// left as facets, the corners fold the tori under them.
#[test]
fn a_drafted_boss_s_faceted_corners_come_back_cones() {
    let mesh = drafted_boss(4, 4);
    let mut model = Model::new();
    let built = solid_from_mesh(&mut model, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(built.closed, "{:?}", built.report);
    let diagnosis = check(&model, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    assert_eq!(
        built.report.curved_faceted, 0,
        "{:?}",
        built.report.fallbacks
    );
    let (mut cones, mut tori) = (0, 0);
    for face in explore_unique(&model, &built.shape, ShapeType::Face).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        match model.geometry().surface(data.surface).unwrap() {
            ogeom::geom::SurfaceGeometry::Cone(c) => {
                assert!(
                    (c.cone().half_angle().abs() - 2.0_f64.to_radians()).abs() < 1e-9,
                    "{c:?}"
                );
                cones += 1;
            }
            ogeom::geom::SurfaceGeometry::Torus(_) => tori += 1,
            _ => {}
        }
    }
    assert_eq!((cones, tori), (4, 4));
    let drawn = ogeom::mesh::triangulate(&model, &built.shape, Deflection::default(), T).unwrap();
    assert!(drawn.is_closed());
}

/// A blind hole ending in a drill's point, a cone of 59 degrees half angle
/// closing at its apex inside its one rim. The cone comes back one face
/// running to its apex, with a ruling for its seam and the apex an edge of
/// no length.
/// Faces facing their surfaces away (a bore crossed by a hole, a drill
/// point, a pocket rounded into a ball) keep their material on the left of
/// their rings: each edge between two converted faces is walked once each
/// way, and the part measures what it was drawn from.
#[test]
fn faces_facing_their_surfaces_away_walk_each_edge_once_each_way() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 12.0), T)
        .unwrap()
        .shape;
    let at = |z: f64| Frame::new(Point::new(10.0, 10.0, z), Direction::Z, Direction::X, T).unwrap();
    let (radius, half) = (3.175_f64, 59.0_f64.to_radians());
    let depth = radius / half.tan();
    let shank = ogeom::algo::make_cylinder(&mut model, at(5.0), radius, 8.0, T)
        .unwrap()
        .shape;
    let point = ogeom::algo::make_cone(&mut model, at(5.0 - depth), 0.0, radius, depth, T)
        .unwrap()
        .shape;
    let drill = ogeom::boolean::fuse(&mut model, &shank, &point, T)
        .unwrap()
        .shape;
    let ball = ogeom::algo::make_sphere(&mut model, at(12.0), 6.0, T)
        .unwrap()
        .shape;
    let bore = ogeom::algo::make_cylinder(&mut model, at(-1.0), 3.0, 14.0, T)
        .unwrap()
        .shape;
    let across = Frame::new(Point::new(-1.0, 10.0, 4.0), Direction::X, Direction::Y, T).unwrap();
    let cross = ogeom::algo::make_cylinder(&mut model, across, 1.5, 22.0, T)
        .unwrap()
        .shape;
    let bored = ogeom::boolean::cut(&mut model, &block, &bore, T)
        .unwrap()
        .shape;
    let parts = [
        (
            "a drill point",
            ogeom::boolean::cut(&mut model, &block, &drill, T),
        ),
        (
            "a ball pocket",
            ogeom::boolean::cut(&mut model, &block, &ball, T),
        ),
        (
            "a crossed bore",
            ogeom::boolean::cut(&mut model, &bored, &cross, T),
        ),
    ];
    for (name, part) in parts {
        let part = part.unwrap().shape;
        let mesh =
            ogeom::mesh::triangulate(&model, &part, Deflection::with_chord(0.02).unwrap(), T)
                .unwrap();
        let mut back = Model::new();
        let built = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
        assert!(built.closed, "{name}: {:?}", built.report);
        let turned = explore_unique(&back, &built.shape, ShapeType::Face)
            .unwrap()
            .iter()
            .filter(|f| f.orientation() == ogeom::topo::Orientation::Reversed)
            .count();
        assert!(turned > 0, "{name}: no face faces its surface away");
        let diagnosis = check(&back, &built.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
        assert_eq!(edges_walked_one_way(&back, &built.shape), 0, "{name}");
        let (want, got) = (volume(&model, &part), volume(&back, &built.shape));
        assert!(
            (got - want).abs() < want * 1e-4,
            "{name}: {got} against {want}"
        );
    }
}

#[test]
fn a_drill_point_comes_back_a_cone_to_its_apex() {
    let mut model = Model::new();
    let (radius, half) = (3.175_f64, 59.0_f64.to_radians());
    let depth = radius / half.tan();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 12.0), T)
        .unwrap()
        .shape;
    let at = |z: f64| Frame::new(Point::new(10.0, 10.0, z), Direction::Z, Direction::X, T).unwrap();
    let bore = ogeom::algo::make_cylinder(&mut model, at(5.0), radius, 8.0, T)
        .unwrap()
        .shape;
    let point = ogeom::algo::make_cone(&mut model, at(5.0 - depth), 0.0, radius, depth, T)
        .unwrap()
        .shape;
    let drill = ogeom::boolean::fuse(&mut model, &bore, &point, T)
        .unwrap()
        .shape;
    let part = ogeom::boolean::cut(&mut model, &block, &drill, T)
        .unwrap()
        .shape;
    let mesh =
        ogeom::mesh::triangulate(&model, &part, Deflection::with_chord(0.02).unwrap(), T).unwrap();
    let mut back = Model::new();
    let built = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(built.closed, "{:?}", built.report);
    let diagnosis = check(&back, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    assert_eq!(
        built.report.curved_faceted, 0,
        "{:?}",
        built.report.fallbacks
    );
    let cones: Vec<f64> = explore_unique(&back, &built.shape, ShapeType::Face)
        .unwrap()
        .iter()
        .filter_map(|face| {
            let data = back.node(face).unwrap().data().as_face().unwrap();
            match back.geometry().surface(data.surface).unwrap() {
                ogeom::geom::SurfaceGeometry::Cone(c) => Some(c.cone().half_angle()),
                _ => None,
            }
        })
        .collect();
    assert_eq!(cones.len(), 1, "{:?}", built.report);
    assert!((cones[0].abs() - half).abs() < 1e-6, "{cones:?}");
    let want = volume(&model, &part);
    let got = volume(&back, &built.shape);
    assert!((got - want).abs() < want * 1e-5, "{got} against {want}");
    let drawn = ogeom::mesh::triangulate(&back, &built.shape, Deflection::default(), T).unwrap();
    assert!(drawn.is_closed());
}
