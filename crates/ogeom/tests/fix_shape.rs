//! `fix_shape`: one pass over a shape nobody promised was well-formed.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, make_box, make_compound, make_face, make_polygon, make_vertex, surface_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{PlaneSurface, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Plane, Point, Point2};
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
    // The model takes a bag. Only the checker knows it is no path. The
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
    assert_eq!(fixed.report.pcurve_gaps_stated, 0);
    assert_eq!(fixed.report.curve_ends_stated, 0);
    assert!(fixed.report.after.is_valid());
    let _ = make_vertex(&mut model, Point::ORIGIN);
}

/// A substitution rebuilds every edge, wire and face above the substituted
/// vertex, and a rebuilt node is shared by every occurrence of the old one.
/// Rebuilt through an occurrence that reverses it, a wire must keep its
/// walk reversed exactly once, or the collapse of a degenerate edge, which
/// substitutes a vertex, leaves the neighbouring wire gaping.
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

/// A box with one face flagged the wrong way: its loops and its
/// neighbours agree about which way is out, the flag does not, and `check`
/// says so. The fix turns that face and no other, and the box measures
/// right again.
#[test]
fn a_face_facing_into_the_material_is_turned() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let face = explore_unique(&model, &block, ShapeType::Face).unwrap()[2].clone();
    // A fresh node with the face's data, flagged the other way.
    let data = model.node(&face).unwrap().data().as_face().unwrap().clone();
    let wires = model.children_of(&Shape::of(face.node())).unwrap();
    let copy = model.add_face(data, &wires).unwrap();
    // The replacement takes the occurrence's own direction on top.
    let mut reshape = ogeom::heal::Reshape::new();
    reshape.replace(&face, copy.reversed());
    let turned = reshape.apply(&mut model, &block).unwrap().shape;
    let broken = ogeom::algo::check(&model, &turned, T).unwrap();
    assert!(!broken.is_valid());
    // The flag turned, the ring did not: each of the face's four edges is
    // walked the same way by the face beside it.
    let same_way = broken
        .problems
        .iter()
        .filter(|p| p.what.contains("walked the same way"))
        .count();
    assert_eq!(same_way, 4, "{broken}");

    let fixed = ogeom::heal::fix_shape(&mut model, &turned, T).unwrap();
    assert_eq!(fixed.report.faces_turned, 1);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let v = ogeom::algo::volume_properties(&model, &fixed.shape, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((v - 1000.0).abs() < 1e-9, "{v}");
}

/// A box with every face flagged the wrong way: the faces agree with each
/// other along every edge, so only the material says it is inside out.
/// `check` names every face, and the fix turns them all back.
#[test]
fn a_solid_inside_out_as_a_whole_is_turned() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let faces: Vec<Shape> = explore_unique(&model, &block, ShapeType::Face)
        .unwrap()
        .into_iter()
        .map(|f| f.reversed())
        .collect();
    let shell = ogeom::algo::make_shell(&mut model, &faces).unwrap().shape;
    // Put together raw: `make_solid` would turn the shell to face out.
    let inside_out = model.add_solid(&[shell]).unwrap();
    let broken = ogeom::algo::check(&model, &inside_out, T).unwrap();
    assert!(!broken.is_valid());
    assert_eq!(broken.problems.len(), 6, "{broken}");

    let fixed = ogeom::heal::fix_shape(&mut model, &inside_out, T).unwrap();
    assert_eq!(fixed.report.faces_turned, 6);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let v = ogeom::algo::volume_properties(&model, &fixed.shape, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((v - 1000.0).abs() < 1e-9, "{v}");
}

/// A void whose faces are placed by a location, and read inside out (its
/// faces turned into the solid and its shell reversed as well): the faces
/// are turned right way out where they are, the void stays inside its
/// solid, and the repaired solid checks clean.
#[test]
fn a_located_void_turned_right_way_out_stays_in_place() {
    use ogeom::topo::{Filter, explore};
    let mut model = Model::new();
    let outer = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let outer_shell = explore(&model, &outer, Filter::OfType(ShapeType::Shell)).unwrap()[0].clone();
    let local = make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T)
        .unwrap()
        .shape;
    let local_shell = explore(&model, &local, Filter::OfType(ShapeType::Shell)).unwrap()[0].clone();
    let place = Location::of(model.add_datum(ogeom::math::Transform::translation(
        ogeom::math::Vector::new(4.0, 4.0, 4.0),
    )));
    let faces: Vec<Shape> = explore(&model, &local_shell, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .map(|f| f.reversed().moved(&place))
        .collect();
    let void = model.add_shell(&faces).unwrap().reversed();
    // Put together raw: `make_solid` would turn the void to face in.
    let solid = model.add_solid(&[outer_shell, void]).unwrap();

    let fixed = ogeom::heal::fix_shape(&mut model, &solid, T).unwrap();
    assert_eq!(fixed.report.faces_turned, 6, "{:?}", fixed.report);
    let shells = explore(&model, &fixed.shape, Filter::OfType(ShapeType::Shell)).unwrap();
    assert_eq!(shells.len(), 2);
    let bounds = ogeom::algo::tight_bounds(&model, &shells[1], T).unwrap();
    let (lo, hi) = (bounds.low().unwrap(), bounds.high().unwrap());
    for (got, want) in [
        (lo, Point::new(4.0, 4.0, 4.0)),
        (hi, Point::new(6.0, 6.0, 6.0)),
    ] {
        assert!(got.distance(want) < 1e-6, "void bounds {lo:?}..{hi:?}");
    }
    let volume = ogeom::algo::volume_properties(&model, &fixed.shape, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((volume - 992.0).abs() < 1e-6, "{volume}");
    let diagnosis = check(&model, &fixed.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
}

const DRUM_RADIUS: f64 = 5.0;
const DRUM_HEIGHT: f64 = 6.0;

/// A drum of radius 5 and height 6 on the XY plane whose round side is two
/// faces on one cylinder: `u` from 0 to pi and from pi to 2 pi, joined at
/// `u = pi` by a ruling with one pcurve and at the cylinder's seam by a
/// ruling `seam` holding two plain pcurves, the face on `u` from 0 to pi's
/// column at `u = 0` first and the other's at `u = 2 pi` after it.
fn drum_with_a_seam_as_two_pcurves(model: &mut Model) -> (Shape, Shape) {
    use ogeom::geom::{CircleCurve, CylinderSurface, Line2d, LineCurve};
    use ogeom::math::{Circle, Cylinder, Point2};
    let pi = core::f64::consts::PI;
    let (r, h) = (DRUM_RADIUS, DRUM_HEIGHT);
    let cylinder = model.geometry_mut().add_surface(
        CylinderSurface::new(Cylinder::new(Frame::WORLD, r, T).unwrap(), (0.0, h))
            .unwrap()
            .into(),
    );
    let at = |x: f64, z: f64| Point::new(x, 0.0, z);
    let [p0, p1, q0, q1] = [at(r, 0.0), at(-r, 0.0), at(r, h), at(-r, h)];
    let [v_p0, v_p1, v_q0, v_q1] = [p0, p1, q0, q1].map(|p| make_vertex(model, p).shape);
    let level = |z: f64| {
        let frame = Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X, T).unwrap();
        CircleCurve::new(Circle::new(frame, r, T).unwrap())
    };
    let arc = |model: &mut Model, z: f64, range: (f64, f64), from: &Shape, to: &Shape| {
        ogeom::algo::make_edge_between(model, level(z).into(), range, from, to, T)
            .unwrap()
            .shape
    };
    let b1 = arc(model, 0.0, (0.0, pi), &v_p0, &v_p1);
    let b2 = arc(model, 0.0, (pi, 2.0 * pi), &v_p1, &v_p0);
    let t1 = arc(model, h, (0.0, pi), &v_q0, &v_q1);
    let t2 = arc(model, h, (pi, 2.0 * pi), &v_q1, &v_q0);
    let rule = |model: &mut Model, a: Point, b: Point, from: &Shape, to: &Shape| {
        let line = LineCurve::segment(a, b, T).unwrap();
        ogeom::algo::make_edge_between(model, line.into(), (0.0, h), from, to, T)
            .unwrap()
            .shape
    };
    let seam = rule(model, p0, q0, &v_p0, &v_q0);
    let middle = rule(model, p1, q1, &v_p1, &v_q1);
    let chart = |model: &mut Model, edge: &Shape, from: (f64, f64), to: (f64, f64)| {
        let (a, b) = (Point2::new(from.0, from.1), Point2::new(to.0, to.1));
        let line = Line2d::segment(a, b, T).unwrap();
        ogeom::algo::attach_pcurve(
            model,
            edge,
            line.into(),
            cylinder,
            Location::identity(),
            (0.0, a.distance(b)),
        )
        .unwrap();
    };
    chart(model, &b1, (0.0, 0.0), (pi, 0.0));
    chart(model, &b2, (pi, 0.0), (2.0 * pi, 0.0));
    chart(model, &t1, (0.0, h), (pi, h));
    chart(model, &t2, (pi, h), (2.0 * pi, h));
    chart(model, &middle, (pi, 0.0), (pi, h));
    chart(model, &seam, (0.0, 0.0), (0.0, h));
    chart(model, &seam, (2.0 * pi, 0.0), (2.0 * pi, h));
    let side = |model: &mut Model, ring: &[Shape]| {
        let wire = ogeom::algo::make_wire(model, ring, T).unwrap().shape;
        ogeom::algo::make_face_on(model, cylinder, &[wire], T)
            .unwrap()
            .shape
    };
    let front = side(
        model,
        &[b1.clone(), middle.clone(), t1.reversed(), seam.reversed()],
    );
    let back = side(
        model,
        &[b2.clone(), seam.clone(), t2.reversed(), middle.reversed()],
    );
    let cap = |model: &mut Model, z: f64, normal: Direction, ring: Vec<Shape>| {
        let frame = Frame::new(Point::new(0.0, 0.0, z), normal, Direction::X, T).unwrap();
        ogeom::algo::make_face_with_pcurves(
            model,
            SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(frame))),
            &[ring],
            T,
        )
        .unwrap()
        .shape
    };
    let bottom = cap(
        model,
        0.0,
        -Direction::Z,
        vec![b2.reversed(), b1.reversed()],
    );
    let top = cap(model, h, Direction::Z, vec![t1, t2]);
    let shell = ogeom::algo::make_shell(model, &[front, back, bottom, top])
        .unwrap()
        .shape;
    let solid = ogeom::algo::make_solid(model, &[shell]).unwrap().shape;
    (solid, seam)
}

/// A seam between two faces on one surface stored as two plain pcurves a
/// period apart is found by `check`: the lookup by surface hands both faces
/// the first column, and the face on the other reads its boundary a period
/// off. Made a seam, each face walks its own column, `check` finds nothing,
/// and a boolean through the seam gives half the drum.
#[test]
fn a_seam_stored_as_two_plain_pcurves_is_made_a_seam() {
    use ogeom::topo::EdgeRepr;
    let mut model = Model::new();
    let (drum, seam) = drum_with_a_seam_as_two_pcurves(&mut model);
    let found = check(&model, &drum, T).unwrap();
    let broken = found.of(ogeom::algo::Severity::Broken);
    assert!(
        broken
            .iter()
            .any(|p| p.kind == ShapeType::Edge && p.at.is_same(&seam)),
        "the edge holding both columns: {found}"
    );
    assert!(
        broken.iter().any(|p| p.kind == ShapeType::Face),
        "the face reading the other's column: {found}"
    );

    let fixed = ogeom::heal::fix_shape(&mut model, &drum, T).unwrap();
    assert!(!fixed.report.before.is_valid());
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let mut seams = Vec::new();
    for edge in explore_unique(&model, &fixed.shape, ShapeType::Edge).unwrap() {
        let data = model.node(&edge).unwrap().data().as_edge().unwrap();
        for r in &data.representations {
            if let EdgeRepr::Seam {
                forward, reversed, ..
            } = r
            {
                assert_eq!(data.representations.len(), 2, "the curve and the seam");
                seams.push((*forward, *reversed));
            }
        }
    }
    assert_eq!(seams.len(), 1, "one edge made a seam: {seams:?}");
    // The face walking the edge up its length stands on `u = 2 pi`.
    let column = |id| {
        ogeom::geom::Curve2d::point_at(model.geometry().pcurve(id).unwrap(), 0.0, T)
            .unwrap()
            .x
    };
    let tau = core::f64::consts::TAU;
    assert!((column(seams[0].0) - tau).abs() < 1e-9);
    assert!(column(seams[0].1).abs() < 1e-9);

    let whole = core::f64::consts::PI * DRUM_RADIUS * DRUM_RADIUS * DRUM_HEIGHT;
    assert!((volume_of(&model, &fixed.shape) - whole).abs() < whole * 1e-3);
    let half = make_box(
        &mut model,
        Frame::new(Point::new(0.0, -10.0, -1.0), Direction::Z, Direction::X, T).unwrap(),
        (10.0, 20.0, DRUM_HEIGHT + 2.0),
        T,
    )
    .unwrap()
    .shape;
    let common = ogeom::boolean::common(&mut model, &fixed.shape, &half, T)
        .unwrap()
        .shape;
    let found = check(&model, &common, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let v = volume_of(&model, &common);
    assert!(
        (v - whole / 2.0).abs() < whole * 1e-3,
        "{v} against {whole}"
    );
}

/// The repair reaches only the shape it is given: in a compound noted as
/// holding the drum's nodes, as an operation's result holds what it passed
/// through from its operands, joining the seam's columns leaves the drum,
/// which still holds them, as it was.
#[test]
fn joining_seam_columns_leaves_an_operand_as_it_was() {
    let mut model = Model::new();
    let (drum, seam) = drum_with_a_seam_as_two_pcurves(&mut model);
    let since = model.node_count();
    let holder = make_compound(&mut model, std::slice::from_ref(&drum))
        .unwrap()
        .shape;
    model.note_held(&holder, since).unwrap();
    assert_eq!(
        ogeom::heal::join_seam_columns(&mut model, &holder, T).unwrap(),
        1
    );
    let found = check(&model, &holder, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let plain = model
        .node(&seam)
        .unwrap()
        .data()
        .as_edge()
        .unwrap()
        .representations
        .iter()
        .filter(|r| matches!(r, ogeom::topo::EdgeRepr::PCurve { .. }))
        .count();
    assert_eq!(plain, 2, "the drum's own edge keeps both pcurves");
    assert!(!check(&model, &drum, T).unwrap().is_valid());
}

/// A box with one straight pcurve moved `off` across its plane, parallel
/// to where it was: lifted, it runs `off` beside its edge's curve.
fn box_with_a_pcurve_off(model: &mut Model, off: f64) -> Shape {
    use ogeom::geom::{Line2d, PlanarCurve};
    use ogeom::topo::EdgeRepr;
    let block = make_box(model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let edge = explore_unique(model, &block, ShapeType::Edge).unwrap()[0].clone();
    let NodeData::Edge(data) = model.node(&edge).unwrap().data().clone() else {
        panic!("an edge");
    };
    let (slot, id) = data
        .representations
        .iter()
        .enumerate()
        .find_map(|(i, repr)| match repr {
            EdgeRepr::PCurve { curve, .. } => Some((i, *curve)),
            _ => None,
        })
        .expect("a box edge has a pcurve on each plane");
    let Some(PlanarCurve::Line(line)) = model.geometry().pcurve(id).cloned() else {
        panic!("a box's pcurves are lines");
    };
    let mut axis = line.axis();
    let d = axis.direction;
    axis.location = Point2::new(axis.location.x - d.y() * off, axis.location.y + d.x() * off);
    let (start, end) = ogeom::geom::Curve2d::domain(&line);
    let moved = model
        .geometry_mut()
        .add_pcurve(Line2d::over(axis, start, end).unwrap().into());
    if let NodeData::Edge(data) = model.node_mut(&edge).unwrap().data_mut()
        && let EdgeRepr::PCurve { curve, .. } = &mut data.representations[slot]
    {
        *curve = moved;
    }
    block
}

/// A pcurve a hundredth of a millimetre beside its curve, farther than the
/// edge states: the edge widens to the gap, as a reader records the same
/// slop, and the box checks clean and measures what it did.
#[test]
fn a_pcurve_off_its_curve_within_reach_widens_its_edge() {
    let mut model = Model::new();
    let off = 0.01;
    let block = box_with_a_pcurve_off(&mut model, off);
    let before = check(&model, &block, T).unwrap();
    assert!(
        before
            .problems
            .iter()
            .any(|p| p.what.contains("leaves its curve")),
        "{before}"
    );

    let fixed = ogeom::heal::fix_shape(&mut model, &block, T).unwrap();
    assert_eq!(fixed.report.pcurve_gaps_stated, 1);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let widened = explore_unique(&model, &fixed.shape, ShapeType::Edge)
        .unwrap()
        .iter()
        .map(|e| model.tolerance_of(e).unwrap().unwrap().get())
        .fold(0.0, f64::max);
    assert!(
        widened >= off && widened < off * (1.0 + 1e-3),
        "the widest edge states the gap: {widened}"
    );
    let v = ogeom::algo::volume_properties(&model, &fixed.shape, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((v - 1000.0).abs() < 1e-6, "{v}");
}

/// A pcurve five millimetres off its curve is no slop: widening the edge
/// that far would call two places five millimetres apart one. It is left
/// as it is and still reported.
#[test]
fn a_pcurve_off_its_curve_past_reach_is_left_reported() {
    let mut model = Model::new();
    let block = box_with_a_pcurve_off(&mut model, 5.0);
    let fixed = ogeom::heal::fix_shape(&mut model, &block, T).unwrap();
    assert_eq!(fixed.report.pcurve_gaps_stated, 0);
    assert!(
        fixed
            .report
            .after
            .problems
            .iter()
            .any(|p| p.what.contains("leaves its curve")),
        "{}",
        fixed.report.after
    );
}

/// A box corner moved two hundredths of a millimetre off the three lines
/// that end at it, farther than any of them or the corner states: the
/// corner widens to the gap and the box checks clean.
#[test]
fn a_vertex_its_curves_miss_within_reach_widens() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let corner = explore_unique(&model, &block, ShapeType::Vertex).unwrap()[0].clone();
    let miss = 0.02;
    if let NodeData::Vertex(data) = model.node_mut(&corner).unwrap().data_mut() {
        data.point = Point::new(data.point.x, data.point.y, data.point.z + miss);
    }
    let before = check(&model, &block, T).unwrap();
    let stops = before
        .problems
        .iter()
        .filter(|p| p.what.starts_with("curve stops"))
        .count();
    assert_eq!(stops, 3, "{before}");

    let fixed = ogeom::heal::fix_shape(&mut model, &block, T).unwrap();
    assert_eq!(fixed.report.curve_ends_stated, 1);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let reach = explore_unique(&model, &fixed.shape, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| model.tolerance_of(v).unwrap().unwrap().get())
        .fold(0.0, f64::max);
    assert!(
        reach >= miss && reach < miss * (1.0 + 1e-9),
        "the corner states the gap: {reach}"
    );
}

/// A corner five millimetres off its lines is a missing piece, not slop:
/// it stays as it is and reported.
#[test]
fn a_vertex_its_curves_miss_past_reach_is_left_reported() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let corner = explore_unique(&model, &block, ShapeType::Vertex).unwrap()[0].clone();
    if let NodeData::Vertex(data) = model.node_mut(&corner).unwrap().data_mut() {
        data.point = Point::new(data.point.x, data.point.y, data.point.z + 5.0);
    }
    let fixed = ogeom::heal::fix_shape(&mut model, &block, T).unwrap();
    assert_eq!(fixed.report.curve_ends_stated, 0);
    assert!(
        fixed
            .report
            .after
            .problems
            .iter()
            .any(|p| p.what.starts_with("curve stops")),
        "{}",
        fixed.report.after
    );
}

/// A box's shell, not a solid, with one face flagged the wrong way: its
/// four edges are walked the same way by it and the face beside each.
/// Nothing is left to sew, and the face is still turned back, as sewing
/// turns it.
#[test]
fn a_shell_s_face_walking_its_edges_against_its_neighbours_is_turned() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let shell = explore_unique(&model, &block, ShapeType::Shell).unwrap()[0].clone();
    let face = explore_unique(&model, &shell, ShapeType::Face).unwrap()[2].clone();
    let data = model.node(&face).unwrap().data().as_face().unwrap().clone();
    let wires = model.children_of(&Shape::of(face.node())).unwrap();
    let copy = model.add_face(data, &wires).unwrap();
    let mut reshape = ogeom::heal::Reshape::new();
    reshape.replace(&face, copy.reversed());
    let turned = reshape.apply(&mut model, &shell).unwrap().shape;
    let same_way = |d: &ogeom::algo::Diagnosis| {
        d.problems
            .iter()
            .filter(|p| p.what.starts_with("walked the same way"))
            .count()
    };
    assert_eq!(same_way(&check(&model, &turned, T).unwrap()), 4);

    let fixed = ogeom::heal::fix_shape(&mut model, &turned, T).unwrap();
    assert_eq!(fixed.report.sewn, Some((0, 0)), "nothing to join");
    assert_eq!(model.kind_of(&fixed.shape).unwrap(), ShapeType::Shell);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let solid = model.add_solid(std::slice::from_ref(&fixed.shape)).unwrap();
    let v = ogeom::algo::volume_properties(&model, &solid, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((v - 1000.0).abs() < 1e-9, "{v}");
}

/// Two boxes overlapping, as two shells of one solid: the boundary passes
/// through itself, and a face of each box inside the other reads facing
/// into material. Each box's faces agree along every edge, and most of
/// them read facing out, so no face is turned: turning the faces that
/// misread alone would set them against every neighbour.
#[test]
fn faces_misread_on_a_boundary_crossing_itself_are_not_turned() {
    let mut model = Model::new();
    let a = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let b = make_box(
        &mut model,
        Frame::WORLD.with_origin(Point::new(5.0, -1.0, -1.0)),
        (10.0, 12.0, 12.0),
        T,
    )
    .unwrap()
    .shape;
    let shells: Vec<Shape> = [a, b]
        .iter()
        .map(|s| explore_unique(&model, s, ShapeType::Shell).unwrap()[0].clone())
        .collect();
    let solid = model.add_solid(&shells).unwrap();
    let inward = ogeom::algo::inside_out_faces(&model, &solid, T).unwrap();
    assert!(!inward.is_empty(), "a face inside the other box reads in");

    let fixed = ogeom::heal::fix_shape(&mut model, &solid, T).unwrap();
    assert_eq!(fixed.report.faces_turned, 0);
    assert!(
        !fixed
            .report
            .after
            .problems
            .iter()
            .any(|p| p.what.starts_with("walked the same way")),
        "{}",
        fixed.report.after
    );
    assert_eq!(fixed.report.after, fixed.report.before);
}

/// The shell of an axis-aligned box from `low` to `high`.
fn box_shell(model: &mut Model, low: Point, high: Point) -> Shape {
    let solid = make_box(
        model,
        Frame::WORLD.with_origin(low),
        (high.x - low.x, high.y - low.y, high.z - low.z),
        T,
    )
    .unwrap()
    .shape;
    explore_unique(model, &solid, ShapeType::Shell).unwrap()[0].clone()
}

/// Four closed shells handed over in no order: a block, a cavity inside
/// it, an island inside the cavity, and a block apart. Three solids come
/// back: the block with its cavity a void, the island, the block apart.
/// Every shell faces out as given; the cavity's is turned to face into it.
#[test]
fn solids_are_made_of_closed_shells_by_which_holds_which() {
    let mut model = Model::new();
    let p = Point::new;
    let block = box_shell(&mut model, p(0.0, 0.0, 0.0), p(10.0, 10.0, 10.0));
    let cavity = box_shell(&mut model, p(2.0, 2.0, 2.0), p(8.0, 8.0, 8.0));
    let island = box_shell(&mut model, p(4.0, 4.0, 4.0), p(6.0, 6.0, 6.0));
    let apart = box_shell(&mut model, p(20.0, 0.0, 0.0), p(25.0, 5.0, 5.0));
    let built = ogeom::heal::solids_from_shells(
        &mut model,
        &[island, cavity.clone(), apart, block.clone()],
        T,
    )
    .unwrap();
    let solids = explore_unique(&model, &built.shape, ShapeType::Solid).unwrap();
    let mut volumes: Vec<f64> = solids
        .iter()
        .map(|s| {
            assert!(check(&model, s, T).unwrap().is_valid());
            ogeom::algo::volume_properties(&model, s, Deflection::default(), T)
                .unwrap()
                .mass
        })
        .collect();
    volumes.sort_by(f64::total_cmp);
    let expected = [8.0, 125.0, 1000.0 - 216.0];
    assert_eq!(volumes.len(), 3, "{volumes:?}");
    for (v, e) in volumes.iter().zip(expected) {
        assert!((v - e).abs() < 1e-9, "{volumes:?}");
    }
    let held = built.history.generated(&cavity);
    assert_eq!(
        held,
        built.history.generated(&block),
        "the cavity is the block's void"
    );
}

/// Two blocks that overlap: a corner of each inside the other, the rest
/// outside. Neither holds the other, nor do they stand apart, and no
/// solid is made of them.
#[test]
fn crossing_shells_make_no_solid() {
    let mut model = Model::new();
    let p = Point::new;
    let a = box_shell(&mut model, p(0.0, 0.0, 0.0), p(10.0, 10.0, 10.0));
    let b = box_shell(&mut model, p(5.0, 5.0, 5.0), p(15.0, 15.0, 15.0));
    let refused = ogeom::heal::solids_from_shells(&mut model, &[a, b], T).unwrap_err();
    assert!(refused.to_string().contains("cross"), "{refused}");
}

/// `solid` with every degenerate edge taken out of its wires: the
/// boundary still meets itself in space at each pole, and its chart is
/// open along the pole's row.
fn without_degenerate_edges(model: &mut Model, solid: &Shape) -> Shape {
    let degenerate: Vec<Shape> = explore_unique(model, solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| model.node(e).unwrap().data().as_edge().unwrap().degenerate)
        .collect();
    assert!(!degenerate.is_empty());
    let mut reshape = ogeom::heal::Reshape::new();
    for edge in &degenerate {
        reshape.remove(edge);
    }
    reshape.apply(model, solid).unwrap().shape
}

fn volume_of_solid(model: &Model, solid: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, solid, Deflection::default(), T)
        .unwrap()
        .mass
}

/// A sphere whose boundary walks its seam up and down with no edge at
/// either pole: `check` passes it, since both joints meet at a pole, but
/// it meshes to nothing and measures no volume. Each pole gets its
/// degenerate edge back, and the sphere measures its closed form.
#[test]
fn a_sphere_missing_its_pole_edges_gets_them_back() {
    let mut model = Model::new();
    let radius = 2.0;
    let sphere = ogeom::algo::make_sphere(&mut model, Frame::WORLD, radius, T)
        .unwrap()
        .shape;
    let stripped = without_degenerate_edges(&mut model, &sphere);
    assert!(volume_of_solid(&model, &stripped).abs() < 1e-9);

    let fixed = ogeom::heal::fix_shape(&mut model, &stripped, T).unwrap();
    assert_eq!(fixed.report.pole_edges_made, 2);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let exact = 4.0 / 3.0 * core::f64::consts::PI * radius.powi(3);
    let v = volume_of_solid(&model, &fixed.shape);
    assert!((v - exact).abs() < 1e-9 * exact, "{v} against {exact}");

    let intact = ogeom::heal::fix_shape(&mut model, &sphere, T).unwrap();
    assert_eq!(intact.report.pole_edges_made, 0, "a sphere with its poles");
}

/// A cone with no edge at its apex: it still meshes, but measures only
/// approximately. The apex gets its degenerate edge back, and the cone
/// measures a third of its base times its height.
#[test]
fn a_cone_missing_its_apex_edge_gets_it_back() {
    let mut model = Model::new();
    let (radius, height) = (2.0, 3.0);
    let cone = ogeom::algo::make_cone(&mut model, Frame::WORLD, radius, 0.0, height, T)
        .unwrap()
        .shape;
    let exact = core::f64::consts::PI * radius * radius * height / 3.0;
    let stripped = without_degenerate_edges(&mut model, &cone);
    assert!((volume_of_solid(&model, &stripped) - exact).abs() > 1e-6);

    let fixed = ogeom::heal::fix_shape(&mut model, &stripped, T).unwrap();
    assert_eq!(fixed.report.pole_edges_made, 1);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let v = volume_of_solid(&model, &fixed.shape);
    assert!((v - exact).abs() < 1e-9 * exact, "{v} against {exact}");
}

/// A 10 mm box with one edge cut short by `gap` at its start and given a
/// vertex of its own there: the two faces that walk the edge each have a
/// gap of `gap` in their boundary, along the edge's own line, so closing
/// it leaves the box as it was.
fn box_with_a_wire_gap(model: &mut Model, gap: f64) -> Shape {
    use ogeom::topo::EdgeRepr;
    let block = make_box(model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let edge = explore_unique(model, &block, ShapeType::Edge).unwrap()[0].clone();
    let stored = Shape::of(edge.node());
    let bounds = model.children_of(&stored).unwrap();
    let point = |model: &Model, v: &Shape| model.node(v).unwrap().data().as_vertex().unwrap().point;
    let (from, to) = (point(model, &bounds[0]), point(model, &bounds[1]));
    let NodeData::Edge(mut data) = model.node(&stored).unwrap().data().clone() else {
        panic!("an edge");
    };
    // A box's curves and pcurves are lines, parameterized alike.
    for repr in &mut data.representations {
        if let EdgeRepr::Curve3d { range, .. } | EdgeRepr::PCurve { range, .. } = repr {
            range.0 += gap / from.distance(to) * (range.1 - range.0);
        }
    }
    let along = from + (to - from) * (gap / from.distance(to));
    let start = make_vertex(model, along).shape;
    let cut = model.add_edge(*data, &[start, bounds[1].clone()]).unwrap();
    let mut reshape = ogeom::heal::Reshape::new();
    reshape.replace(&stored, cut);
    reshape.apply(model, &block).unwrap().shape
}

fn wire_gaps(diagnosis: &ogeom::algo::Diagnosis) -> usize {
    diagnosis
        .problems
        .iter()
        .filter(|p| p.what.contains("does not begin"))
        .count()
}

fn widest_vertex(model: &Model, shape: &Shape) -> f64 {
    explore_unique(model, shape, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| model.tolerance_of(v).unwrap().unwrap().get())
        .fold(0.0, f64::max)
}

/// A box whose boundary gapes a hundredth of a millimetre where one edge
/// begins, at a vertex of its own: the two vertices become one, wide
/// enough to reach both, and the box checks clean and measures 1000.
#[test]
fn a_wire_gap_within_reach_is_closed() {
    let mut model = Model::new();
    let gap = 0.01;
    let block = box_with_a_wire_gap(&mut model, gap);
    let before = check(&model, &block, T).unwrap();
    assert_eq!(wire_gaps(&before), 2, "{before}");
    assert_eq!(
        explore_unique(&model, &block, ShapeType::Vertex)
            .unwrap()
            .len(),
        9
    );

    let fixed = ogeom::heal::fix_shape(&mut model, &block, T).unwrap();
    assert_eq!(fixed.report.gap_vertices_merged, 1);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    assert_eq!(
        explore_unique(&model, &fixed.shape, ShapeType::Vertex)
            .unwrap()
            .len(),
        8
    );
    let reach = widest_vertex(&model, &fixed.shape);
    // As far as the absorbed vertex stood, plus its own tolerance.
    assert!((reach - (gap + T.confusion())).abs() < 1e-12, "{reach}");
    let v = volume_of_solid(&model, &fixed.shape);
    assert!((v - 1000.0).abs() < 1e-9, "{v}");
}

/// Two vertices at one place, one per edge, as edges built apart carry
/// them: they become one, and nothing widens past what they stated.
#[test]
fn a_wire_meeting_at_two_coincident_vertices_shares_one() {
    let mut model = Model::new();
    let block = box_with_a_wire_gap(&mut model, 0.0);
    assert_eq!(wire_gaps(&check(&model, &block, T).unwrap()), 2);

    let fixed = ogeom::heal::fix_shape(&mut model, &block, T).unwrap();
    assert_eq!(fixed.report.gap_vertices_merged, 1);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    assert!(widest_vertex(&model, &fixed.shape) <= T.confusion() * 2.0);
    let v = volume_of_solid(&model, &fixed.shape);
    assert!((v - 1000.0).abs() < 1e-9, "{v}");
}

/// A gap of five millimetres is a missing piece of boundary, not slop: no
/// vertices merge and the gap stays reported.
#[test]
fn a_wire_gap_past_reach_is_left_reported() {
    let mut model = Model::new();
    let block = box_with_a_wire_gap(&mut model, 5.0);
    let fixed = ogeom::heal::fix_shape(&mut model, &block, T).unwrap();
    assert_eq!(fixed.report.gap_vertices_merged, 0);
    assert_eq!(wire_gaps(&fixed.report.after), 2, "{}", fixed.report.after);
}

/// A square face gaping a hundredth of a millimetre at one corner, beside
/// a triangle with an edge across that same gap: merging the two corners
/// would close the triangle's edge on itself and pinch its boundary to a
/// point, so neither merges.
#[test]
fn a_wire_gap_spanned_by_an_edge_elsewhere_is_left() {
    let mut model = Model::new();
    let gap = 0.01;
    let corner =
        |model: &mut Model, x: f64, y: f64| make_vertex(model, Point::new(x, y, 0.0)).shape;
    let a = corner(&mut model, 0.0, 0.0);
    let b = corner(&mut model, 4.0, 0.0);
    let c = corner(&mut model, 4.0, 4.0);
    let d = corner(&mut model, 0.0, 4.0);
    let a2 = corner(&mut model, 0.0, -gap);
    let e = corner(&mut model, -3.0, 0.0);
    let line = |model: &mut Model, p: &Shape, q: &Shape| {
        let at = |v: &Shape| model.node(v).unwrap().data().as_vertex().unwrap().point;
        let (from, to) = (at(p), at(q));
        let curve = ogeom::geom::LineCurve::segment(from, to, T).unwrap().into();
        ogeom::algo::make_edge_between(model, curve, (0.0, from.distance(to)), p, q, T)
            .unwrap()
            .shape
    };
    let square = [
        line(&mut model, &a, &b),
        line(&mut model, &b, &c),
        line(&mut model, &c, &d),
        line(&mut model, &d, &a2),
    ];
    let triangle = [
        line(&mut model, &a2, &a),
        line(&mut model, &a, &e),
        line(&mut model, &e, &a2),
    ];
    let plane = || SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(Frame::WORLD)));
    let mut faces = Vec::new();
    for edges in [&square[..], &triangle[..]] {
        let wire = model.add_wire(edges).unwrap();
        let surface = model.geometry_mut().add_surface(plane());
        faces.push(
            model
                .add_face(FaceData::new(surface, Location::default()), &[wire])
                .unwrap(),
        );
    }
    let both = make_compound(&mut model, &faces).unwrap().shape;
    assert_eq!(wire_gaps(&check(&model, &both, T).unwrap()), 1);

    let fixed = ogeom::heal::fix_shape(&mut model, &both, T).unwrap();
    assert_eq!(fixed.report.gap_vertices_merged, 0);
    assert_eq!(wire_gaps(&fixed.report.after), 1, "{}", fixed.report.after);
    assert_eq!(
        explore_unique(&model, &fixed.shape, ShapeType::Vertex)
            .unwrap()
            .len(),
        6
    );
}
