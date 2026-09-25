//! `fix_shape`: one pass over a shape nobody promised was well-formed.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, make_box, make_compound, make_face, make_polygon, make_vertex, surface_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{PlaneSurface, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Plane, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{FaceData, Location, Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A wire whose edges were handed over in the wrong order is put in the
/// order that walks them.
#[test]
fn a_wire_out_of_order_is_reordered() {
    let mut model = Model::new();
    let square = make_polygon(
        &mut model,
        &[
            Point::new(0.0, 0.0, 0.0),
            Point::new(4.0, 0.0, 0.0),
            Point::new(4.0, 4.0, 0.0),
            Point::new(0.0, 4.0, 0.0),
        ],
        true,
        T,
    )
    .unwrap()
    .shape;
    let mut edges = model.ordered_children_of(&square).unwrap();
    edges.swap(1, 3);
    // The model takes a bag; only the checker knows it is no path. The
    // face is built straight on it, since the builder would refuse.
    let shuffled = model.add_wire(&edges).unwrap();
    let plane = SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(Frame::WORLD)));
    let surface = model.geometry_mut().add_surface(plane);
    let face = model
        .add_face(FaceData::new(surface, Location::default()), &[shuffled])
        .unwrap();
    let before = check(&model, &face, T).unwrap();
    assert!(!before.is_usable(), "the checker sees the gap: {before}");

    let fixed = ogeom::heal::fix_shape(&mut model, &face, T).unwrap();
    assert_eq!(fixed.report.wires_reordered, 1);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let area = surface_properties(&model, &fixed.shape, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((area - 16.0).abs() < 1e-9, "the same square: {area}");
}

/// An edge shorter than its vertices' tolerances collapses to one vertex.
#[test]
fn a_tiny_edge_is_collapsed() {
    let mut model = Model::new();
    // A square with a fifth corner a micron from the fourth.
    let ring = make_polygon(
        &mut model,
        &[
            Point::new(0.0, 0.0, 0.0),
            Point::new(4.0, 0.0, 0.0),
            Point::new(4.0, 4.0, 0.0),
            Point::new(0.0, 4.0, 0.0),
            Point::new(0.0, 4.0 - 1e-6, 0.0),
        ],
        true,
        T,
    )
    .unwrap()
    .shape;
    // The model admits the two corners are one place: a vertex tolerance
    // of ten microns, as an imported vertex might carry.
    let vertices = explore_unique(&model, &ring, ShapeType::Vertex).unwrap();
    let near = vertices
        .iter()
        .find(|v| {
            model
                .node(v)
                .and_then(|n| n.data().as_vertex())
                .is_some_and(|d| d.point.distance(Point::new(0.0, 4.0, 0.0)) < 1e-9)
        })
        .expect("the fourth corner");
    if let Some(node) = model.node_mut(near)
        && let NodeData::Vertex(data) = node.data_mut()
    {
        data.tolerance = data.tolerance.widen_to(1e-5);
    }
    let plane = SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(Frame::WORLD)));
    let face = make_face(&mut model, plane, &[ring], T).unwrap().shape;
    assert_eq!(
        explore_unique(&model, &face, ShapeType::Edge)
            .unwrap()
            .len(),
        5
    );

    let fixed = ogeom::heal::fix_shape(&mut model, &face, T).unwrap();
    assert_eq!(fixed.report.edges_collapsed, 1);
    assert_eq!(
        explore_unique(&model, &fixed.shape, ShapeType::Edge)
            .unwrap()
            .len(),
        4,
        "a square again"
    );
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let area = surface_properties(&model, &fixed.shape, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((area - 16.0).abs() < 1e-4, "the same square: {area}");
}

/// A compound of loose faces is sewn into the shell they make.
#[test]
fn loose_faces_are_sewn() {
    let mut model = Model::new();
    // Six faces of a unit cube, each with vertices and edges of its own:
    // nothing shared but position.
    let p = |x: f64, y: f64, z: f64| Point::new(x, y, z);
    let quads: [[Point; 4]; 6] = [
        [p(0., 0., 0.), p(0., 1., 0.), p(1., 1., 0.), p(1., 0., 0.)],
        [p(0., 0., 1.), p(1., 0., 1.), p(1., 1., 1.), p(0., 1., 1.)],
        [p(0., 0., 0.), p(1., 0., 0.), p(1., 0., 1.), p(0., 0., 1.)],
        [p(0., 1., 0.), p(0., 1., 1.), p(1., 1., 1.), p(1., 1., 0.)],
        [p(0., 0., 0.), p(0., 0., 1.), p(0., 1., 1.), p(0., 1., 0.)],
        [p(1., 0., 0.), p(1., 1., 0.), p(1., 1., 1.), p(1., 0., 1.)],
    ];
    let mut faces = Vec::new();
    for quad in &quads {
        let ring = make_polygon(&mut model, quad, true, T).unwrap().shape;
        let normal = Direction::new((quad[1] - quad[0]).cross(quad[2] - quad[0]), T).unwrap();
        let reference = Direction::new(quad[1] - quad[0], T).unwrap();
        let frame = Frame::new(quad[0], normal, reference, T).unwrap();
        let plane = SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(frame)));
        faces.push(make_face(&mut model, plane, &[ring], T).unwrap().shape);
    }
    let bag = make_compound(&mut model, &faces).unwrap().shape;

    let fixed = ogeom::heal::fix_shape(&mut model, &bag, T).unwrap();
    let (joined, free) = fixed
        .report
        .sewn
        .expect("sewing ran over a compound of faces");
    assert_eq!(free, 0, "the box closes");
    assert!(
        joined >= 12,
        "every edge of the box was found shared: {joined}"
    );
    let shells = explore_unique(&model, &fixed.shape, ShapeType::Shell).unwrap();
    assert_eq!(shells.len(), 1);
    assert!(fixed.report.after.is_usable(), "{}", fixed.report.after);
}

/// A shape with nothing wrong is handed back as it is.
#[test]
fn a_sound_shape_is_left_alone() {
    let mut model = Model::new();
    let solid = make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T)
        .unwrap()
        .shape;
    let fixed = ogeom::heal::fix_shape(&mut model, &solid, T).unwrap();
    assert!(fixed.shape.is_same(&solid), "the same node");
    assert_eq!(fixed.report.wires_reordered, 0);
    assert_eq!(fixed.report.edges_collapsed, 0);
    assert_eq!(fixed.report.edges_trimmed, 0);
    assert!(fixed.report.sewn.is_none(), "a solid is not sewn");
    assert!(fixed.report.after.is_valid());
    let _ = make_vertex(&mut model, Point::ORIGIN);
}

/// A substitution rebuilds every edge, wire and face above the substituted
/// vertex, and a rebuilt node is shared by every occurrence of the old one.
/// Rebuilt through an occurrence that reverses it, a wire came back with
/// its walk reversed once too often, and the collapse of a degenerate
/// edge, which substitutes a vertex, left the neighbouring wire gaping.
/// Every vertex of a prism, whose near cap is its profile reversed, and of
/// a fused pair, whose faces come back in either sense, is substituted in
/// turn, and the result must stay valid.
#[test]
fn a_substituted_vertex_rebuilds_reversed_occurrences_head_to_tail() {
    use ogeom::math::Vector;
    let solids: [fn(&mut Model) -> Shape; 2] = [
        |model| {
            let square = [
                Point::new(0.0, 0.0, 0.0),
                Point::new(2.0, 0.0, 0.0),
                Point::new(2.0, 1.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
            ];
            let wire = ogeom::algo::make_polygon(model, &square, true, T)
                .unwrap()
                .shape;
            let surface = ogeom::geom::PlaneSurface::over(
                Plane::through(Point::ORIGIN, Direction::Z),
                (-1.0, 3.0),
                (-1.0, 2.0),
            )
            .unwrap();
            let face = ogeom::algo::make_face(model, surface.into(), &[wire], T)
                .unwrap()
                .shape;
            let prism = ogeom::algo::make_prism(model, &face, Vector::new(0.0, 0.0, 1.0), T)
                .unwrap()
                .shape;
            // Baked, so the far cap's vertices are nodes of their own rather
            // than the near cap's placed: a substitution replaces a node
            // everywhere it occurs.
            ogeom::algo::baked_shape(model, &prism, T).unwrap().shape
        },
        |model| {
            let a = ogeom::algo::make_box(model, Frame::WORLD, (2.0, 2.0, 2.0), T)
                .unwrap()
                .shape;
            let at = Frame::new(Point::new(1.0, 1.0, 1.0), Direction::Z, Direction::X, T).unwrap();
            let b = ogeom::algo::make_box(model, at, (2.0, 2.0, 2.0), T)
                .unwrap()
                .shape;
            ogeom::boolean::fuse(model, &a, &b, T).unwrap().shape
        },
    ];
    for (which, build) in solids.iter().enumerate() {
        let mut probe = Model::new();
        let built = build(&mut probe);
        let count = explore_unique(&probe, &built, ShapeType::Vertex)
            .unwrap()
            .len();
        for index in 0..count {
            let mut model = Model::new();
            let solid = build(&mut model);
            let vertex = explore_unique(&model, &solid, ShapeType::Vertex).unwrap()[index].clone();
            let point = model
                .node(&vertex)
                .unwrap()
                .data()
                .as_vertex()
                .unwrap()
                .point;
            let point = vertex.transform(model.datums()).unwrap().apply(point);
            let twin = ogeom::algo::make_vertex(&mut model, point).shape;
            let mut reshape = ogeom::heal::Reshape::new();
            reshape.replace(&vertex, twin);
            let rebuilt = reshape.apply(&mut model, &solid).unwrap().shape;
            let diagnosis = ogeom::algo::check(&model, &rebuilt, T).unwrap();
            assert!(
                diagnosis.is_valid(),
                "solid {which} vertex {index}: {diagnosis}"
            );
        }
    }
}

fn volume_of(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// A cube whose corner is bevelled a micron along its whole length has a
/// strip face there: it collapses to one of its long sides, the two walls
/// meet on it, and the cube is six faces again.
#[test]
fn a_strip_face_collapses_to_an_edge() {
    let mut model = Model::new();
    let bevel = 1e-3;
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(10.0, 0.0, 0.0),
        Point::new(10.0, 10.0 - bevel, 0.0),
        Point::new(10.0 - bevel, 10.0, 0.0),
        Point::new(0.0, 10.0, 0.0),
    ];
    let wire = make_polygon(&mut model, &corners, true, T).unwrap().shape;
    let edges =
        ogeom::topo::explore(&model, &wire, ogeom::topo::Filter::OfType(ShapeType::Edge)).unwrap();
    let plane = PlaneSurface::new(Plane::through(Point::ORIGIN, Direction::Z));
    let profile = ogeom::algo::make_face_with_pcurves(&mut model, plane.into(), &[edges], T)
        .unwrap()
        .shape;
    let chamfered = ogeom::algo::make_prism(&mut model, &profile, ogeom::math::Vector::Z * 10.0, T)
        .unwrap()
        .shape;
    assert_eq!(
        explore_unique(&model, &chamfered, ShapeType::Face)
            .unwrap()
            .len(),
        7
    );
    let fixed = ogeom::heal::fix_small_faces(&mut model, &chamfered, 1e-2, T).unwrap();
    assert_eq!((fixed.spots, fixed.strips), (0, 1));
    let shape = &fixed.built.shape;
    assert_eq!(
        explore_unique(&model, shape, ShapeType::Face)
            .unwrap()
            .len(),
        6
    );
    let diagnosis = check(&model, shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    // The walls keep their geometry and meet across the bevel's width,
    // which the merged edge owns as tolerance: the volume is the cube's to
    // within that width over a wall's area.
    let v = volume_of(&model, shape);
    assert!((v - 1000.0).abs() < bevel * 100.0, "{v}");
}

/// A cube's corner cut a micron in is a spot: the little triangle
/// collapses to a point, where the three walls meet again.
#[test]
fn a_spot_face_collapses_to_a_point() {
    let mut model = Model::new();
    let cube = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let normal = Direction::new(ogeom::math::Vector::new(1.0, 1.0, 1.0), T).unwrap();
    let cutting = Plane::through(Point::new(10.0 - 1e-3, 10.0, 10.0), normal);
    let face = ogeom::algo::make_natural_face(
        &mut model,
        PlaneSurface::over(cutting, (-50.0, 50.0), (-50.0, 50.0))
            .unwrap()
            .into(),
    )
    .unwrap()
    .shape;
    let keep = ogeom::algo::make_half_space(&mut model, &face, Point::ORIGIN, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::common(&mut model, &cube, &keep, T)
        .unwrap()
        .shape;
    assert_eq!(
        explore_unique(&model, &cut, ShapeType::Face).unwrap().len(),
        7
    );
    let fixed = ogeom::heal::fix_small_faces(&mut model, &cut, 1e-2, T).unwrap();
    assert_eq!((fixed.spots, fixed.strips), (1, 0));
    let shape = &fixed.built.shape;
    assert_eq!(
        explore_unique(&model, shape, ShapeType::Face)
            .unwrap()
            .len(),
        6
    );
    let diagnosis = check(&model, shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
}

/// A speck beside a part is debris: the solid below the volume goes, the
/// part stays.
#[test]
fn a_small_solid_beside_a_part_is_removed() {
    let mut model = Model::new();
    let part = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(20.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let speck = make_box(&mut model, at, (1e-2, 1e-2, 1e-2), T)
        .unwrap()
        .shape;
    let both = make_compound(&mut model, &[part.clone(), speck.clone()])
        .unwrap()
        .shape;
    let (built, removed) = ogeom::heal::remove_small_solids(&mut model, &both, 1e-3, T).unwrap();
    assert_eq!(removed, 1);
    let solids = explore_unique(&model, &built.shape, ShapeType::Solid).unwrap();
    assert_eq!(solids.len(), 1);
    assert_eq!(solids[0].node(), part.node());
    assert!(built.history.is_deleted(&speck));
}
