//! Non-manifold topology, pinned: the model permits an edge bounding three
//! faces and a compound mixing dimensions, and every query over them answers
//! honestly — traversal finds what is there, closure says not-closed, the
//! native format gives the same bytes back.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::{LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Plane, Point};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// Three rectangular sheets meeting along one shared spine edge — a T-joint
/// plus one more leaf, the textbook non-manifold edge.
fn t_joint(model: &mut Model) -> (Shape, Vec<Shape>) {
    let spine = ogeom::algo::make_edge(
        model,
        LineCurve::segment(Point::ORIGIN, Point::new(10.0, 0.0, 0.0), T)
            .unwrap()
            .into(),
        (0.0, 10.0),
        T,
    )
    .unwrap()
    .shape;
    let spine_ends = explore_unique(model, &spine, ShapeType::Vertex).unwrap();
    // The spine runs origin → (10,0,0); its vertices come back in that order.
    let (v0, v1) = (spine_ends[0].clone(), spine_ends[1].clone());

    // Each leaf rises from the spine along its own direction.
    let leaves = [
        (Direction::Y, Direction::Z),
        (
            Direction::Z,
            Direction::from_coords(0.0, -1.0, 0.0, T).unwrap(),
        ),
        (
            Direction::from_coords(0.0, -0.6, 0.8, T).unwrap(),
            Direction::from_coords(0.0, -0.8, -0.6, T).unwrap(),
        ),
    ];
    let mut faces = Vec::new();
    for (rise, normal) in leaves {
        let d = rise.vector() * 10.0;
        let c1 = ogeom::algo::make_vertex(model, Point::new(10.0, 0.0, 0.0) + d).shape;
        let c0 = ogeom::algo::make_vertex(model, Point::ORIGIN + d).shape;
        let segment = |model: &mut Model, a: &Shape, b: &Shape, pa: Point, pb: Point| {
            ogeom::algo::make_edge_between(
                model,
                LineCurve::segment(pa, pb, T).unwrap().into(),
                (0.0, pa.distance(pb)),
                a,
                b,
                T,
            )
            .unwrap()
            .shape
        };
        let up = segment(
            model,
            &v1,
            &c1,
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0, 0.0, 0.0) + d,
        );
        let across = segment(
            model,
            &c1,
            &c0,
            Point::new(10.0, 0.0, 0.0) + d,
            Point::ORIGIN + d,
        );
        let down = segment(model, &c0, &v0, Point::ORIGIN + d, Point::ORIGIN);
        let frame = Frame::new(Point::ORIGIN, normal, Direction::X, T).unwrap();
        let face = ogeom::algo::make_face_with_pcurves(
            model,
            SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(frame))),
            &[vec![spine.clone(), up, across, down]],
            T,
        )
        .unwrap()
        .shape;
        faces.push(face);
    }
    (spine, faces)
}

#[test]
fn an_edge_bounding_three_faces_is_permitted_and_read_honestly() {
    let mut model = Model::new();
    let (spine, faces) = t_joint(&mut model);

    // The builder accepts the shell; nothing silently rejects legitimate
    // non-manifold topology.
    let shell = ogeom::algo::make_shell(&mut model, &faces).unwrap().shape;

    // Traversal sees one spine edge shared by all three faces.
    let spine_users = faces
        .iter()
        .filter(|face| {
            explore_unique(&model, face, ShapeType::Edge)
                .unwrap()
                .iter()
                .any(|edge| edge.node() == spine.node())
        })
        .count();
    assert_eq!(spine_users, 3, "the spine bounds all three leaves");
    assert_eq!(
        explore_unique(&model, &shell, ShapeType::Edge)
            .unwrap()
            .len(),
        10,
        "one spine and three private edges per leaf"
    );

    // Closure answers what the word means: an edge used three times cannot
    // bound a volume, and the shell says so instead of guessing.
    assert!(!ogeom::algo::is_shell_closed(&model, &shell).unwrap());

    // The mesh layer meshes what is there — each leaf triangulates.
    let done =
        ogeom::mesh::tessellate(&mut model, &shell, ogeom::mesh::Deflection::default(), T).unwrap();
    assert_eq!(done.faces, 3);
    assert!(done.triangles >= 6);
}

#[test]
fn a_mixed_dimension_compound_round_trips_byte_stable() {
    let mut model = Model::new();
    let (spine, faces) = t_joint(&mut model);
    let solid = ogeom::algo::make_box(&mut model, Frame::WORLD, (4.0, 5.0, 6.0), T)
        .unwrap()
        .shape;
    let loose_vertex = ogeom::algo::make_vertex(&mut model, Point::new(-3.0, -3.0, -3.0)).shape;
    let compound = ogeom::algo::make_compound(
        &mut model,
        &[
            solid,
            faces[0].clone(),
            faces[1].clone(),
            faces[2].clone(),
            spine.clone(),
            loose_vertex,
        ],
    )
    .unwrap()
    .shape;

    // Traversal by type over the mixed bag counts what was put in.
    assert_eq!(
        explore_unique(&model, &compound, ShapeType::Solid)
            .unwrap()
            .len(),
        1
    );
    assert_eq!(
        explore_unique(&model, &compound, ShapeType::Face)
            .unwrap()
            .len(),
        6 + 3
    );

    // The native format carries it whole: write, read, write the same
    // bytes, the shared spine still one node used by three faces.
    let text = ogeom::io::native::write(
        &model,
        &[compound],
        ogeom::io::native::WriteOptions::default(),
    )
    .unwrap();
    let (back, roots) = ogeom::io::native::read(&text).unwrap();
    let again = ogeom::io::native::write(&back, &roots, ogeom::io::native::WriteOptions::default())
        .unwrap();
    assert_eq!(text, again);
    // Sharing survives structurally: some single edge node in the
    // read-back model is used by three distinct faces.
    let read_faces = explore_unique(&back, &roots[0], ShapeType::Face).unwrap();
    let mut users: std::collections::BTreeMap<u64, usize> = std::collections::BTreeMap::new();
    for face in &read_faces {
        for edge in explore_unique(&back, face, ShapeType::Edge).unwrap() {
            let mut hasher = std::hash::DefaultHasher::new();
            std::hash::Hash::hash(&edge.node(), &mut hasher);
            *users.entry(std::hash::Hasher::finish(&hasher)).or_default() += 1;
        }
    }
    assert_eq!(
        users.values().max().copied(),
        Some(3),
        "one edge node is used by three faces after the round trip"
    );
}
