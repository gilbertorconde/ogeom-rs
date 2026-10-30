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

/// A 10 mm box whose boundary gapes two millimetres where one edge
/// begins, past the reach within which vertices merge: the gap lies along
/// the cut edge's own line on both faces that walk it, so one straight
/// edge with a pcurve on each fills it, and the box checks clean and
/// measures 1000 again.
#[test]
fn a_wire_gap_past_reach_along_a_line_gets_its_edge() {
    let mut model = Model::new();
    let block = box_with_a_wire_gap(&mut model, 2.0);
    let edges = |model: &Model, s: &Shape| explore_unique(model, s, ShapeType::Edge).unwrap().len();
    assert_eq!(edges(&model, &block), 12);

    let fixed = ogeom::heal::fix_shape(&mut model, &block, T).unwrap();
    assert_eq!(fixed.report.gap_vertices_merged, 0);
    assert_eq!(fixed.report.lacking_edges_made, 1);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    assert_eq!(edges(&model, &fixed.shape), 13);
    assert_eq!(
        explore_unique(&model, &fixed.shape, ShapeType::Vertex)
            .unwrap()
            .len(),
        9
    );
    let v = volume_of_solid(&model, &fixed.shape);
    assert!((v - 1000.0).abs() < 1e-9, "{v}");
}

/// A 10 mm square face whose corner at the origin is missing: one edge
/// stops at (0, 1), the next starts at (1, 0). The gap, √2 across, is
/// shorter than every edge and no other corner is as near, so the
/// straight edge across it is made, and the face, now chamfered, measures
/// 100 − 1/2.
#[test]
fn a_gap_across_a_plane_face_s_corner_is_filled_straight() {
    let mut model = Model::new();
    let at = |model: &mut Model, x: f64, y: f64| make_vertex(model, Point::new(x, y, 0.0)).shape;
    let corners = [
        at(&mut model, 1.0, 0.0),
        at(&mut model, 10.0, 0.0),
        at(&mut model, 10.0, 10.0),
        at(&mut model, 0.0, 10.0),
        at(&mut model, 0.0, 1.0),
    ];
    let mut edges = Vec::new();
    for pair in corners.windows(2) {
        let point = |v: &Shape| model.node(v).unwrap().data().as_vertex().unwrap().point;
        let (from, to) = (point(&pair[0]), point(&pair[1]));
        let curve = ogeom::geom::LineCurve::segment(from, to, T).unwrap().into();
        let edge = ogeom::algo::make_edge_between(
            &mut model,
            curve,
            (0.0, from.distance(to)),
            &pair[0],
            &pair[1],
            T,
        )
        .unwrap()
        .shape;
        let line =
            ogeom::geom::Line2d::segment(Point2::new(from.x, from.y), Point2::new(to.x, to.y), T)
                .unwrap();
        edges.push((edge, line));
    }
    let surface = model
        .geometry_mut()
        .add_surface(SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(
            Frame::WORLD,
        ))));
    for (edge, line) in &edges {
        ogeom::algo::attach_pcurve(
            &mut model,
            edge,
            (*line).into(),
            surface,
            Location::default(),
            ogeom::geom::Curve2d::domain(line),
        )
        .unwrap();
    }
    let ring: Vec<Shape> = edges.into_iter().map(|(e, _)| e).collect();
    let wire = model.add_wire(&ring).unwrap();
    let face = model
        .add_face(FaceData::new(surface, Location::default()), &[wire])
        .unwrap();
    assert_eq!(wire_gaps(&check(&model, &face, T).unwrap()), 1);

    let fixed = ogeom::heal::fix_shape(&mut model, &face, T).unwrap();
    assert_eq!(fixed.report.lacking_edges_made, 1);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let area = surface_properties(&model, &fixed.shape, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((area - 99.5).abs() < 1e-9, "{area}");
}

/// A 10 mm cylinder whose top rim is cut two millimetres of arc short at
/// each end, so it stops either side of the seam's top vertex: each chord
/// from that vertex to a cut end stands off the cylinder by its sagitta,
/// 10 (1 − cos 0.1), about 0.05 mm, and the top's gap lies between the two
/// ends of one edge, which no lacking edge joins. No edge is made and
/// every gap stays reported.
#[test]
fn a_gap_whose_chord_leaves_the_surface_stays_reported() {
    use ogeom::topo::EdgeRepr;
    let mut model = Model::new();
    let radius = 10.0;
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, radius, 10.0, T)
        .unwrap()
        .shape;
    // The top rim: the closed edge whose vertex stands at the top.
    let rim = explore_unique(&model, &drum, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let bounds = model.children_of(&Shape::of(e.node())).unwrap();
            let p = model
                .node(&bounds[0])
                .unwrap()
                .data()
                .as_vertex()
                .unwrap()
                .point;
            bounds.len() == 1 || (bounds[0].node() == bounds[1].node() && p.z > 5.0)
        })
        .unwrap();
    let stored = Shape::of(rim.node());
    let bounds = model.children_of(&stored).unwrap();
    let NodeData::Edge(mut data) = model.node(&stored).unwrap().data().clone() else {
        panic!("an edge");
    };
    let cut = 2.0 / radius;
    let mut ends = None;
    for repr in &mut data.representations {
        if let EdgeRepr::Curve3d { range, curve, .. } = repr {
            *range = (range.0 + cut, range.1 - cut);
            let c = model.geometry().curve(*curve).unwrap();
            use ogeom::geom::Curve3d as _;
            ends = Some((
                c.point_at(range.0, T).unwrap(),
                c.point_at(range.1, T).unwrap(),
            ));
        }
        if let EdgeRepr::PCurve { range, .. } | EdgeRepr::Seam { range, .. } = repr {
            *range = (range.0 + cut, range.1 - cut);
        }
    }
    let (from, to) = ends.unwrap();
    let ends = [
        make_vertex(&mut model, from).shape,
        make_vertex(&mut model, to).shape,
    ];
    assert!(bounds.iter().all(|b| b.node() == bounds[0].node()));
    let cut_rim = model.add_edge(*data, &ends).unwrap();
    let mut reshape = ogeom::heal::Reshape::new();
    reshape.replace(&stored, cut_rim);
    let drum = reshape.apply(&mut model, &drum).unwrap().shape;
    let before = wire_gaps(&check(&model, &drum, T).unwrap());
    assert_eq!(before, 3);

    let fixed = ogeom::heal::fix_shape(&mut model, &drum, T).unwrap();
    assert_eq!(fixed.report.lacking_edges_made, 0);
    assert_eq!(
        wire_gaps(&fixed.report.after),
        before,
        "{}",
        fixed.report.after
    );
}

/// A face with no wires is a face of its whole surface; one that lost the
/// mark saying so is bounded by its domain's sides where the domain is the
/// surface's whole extent (a sphere: a seam and a degenerate edge at each
/// pole; a torus: two seams), and left reported where the domain is only a
/// window on an unbounded sheet (a plane). Bounded, a solid of the sphere's
/// one face is closed by its edges and measures the ball.
#[test]
fn a_face_with_no_wires_is_bounded_where_its_domain_is_its_extent() {
    use ogeom::geom::{SphereSurface, TorusSurface};
    use ogeom::math::{Sphere, Torus};
    use std::f64::consts::PI;
    let mut model = Model::new();
    let unmark = |model: &mut Model, face: &Shape| {
        if let Some(node) = model.node_mut(face)
            && let NodeData::Face(data) = node.data_mut()
        {
            data.natural_restriction = false;
        }
    };
    let no_wires = |d: &ogeom::algo::Diagnosis| {
        d.problems
            .iter()
            .filter(|p| p.what.starts_with("has no wires"))
            .count()
    };
    let area = |model: &Model, shape: &Shape| {
        surface_properties(model, shape, Deflection::default(), T)
            .unwrap()
            .mass
    };
    let sphere = SphereSurface::new(Sphere::new(Frame::WORLD, 2.0, T).unwrap());
    let ball = ogeom::algo::make_natural_face(&mut model, sphere.into())
        .unwrap()
        .shape;
    unmark(&mut model, &ball);
    assert_eq!(no_wires(&check(&model, &ball, T).unwrap()), 1);
    let fixed = ogeom::heal::fix_shape(&mut model, &ball, T).unwrap();
    assert_eq!(fixed.report.natural_bounds_made, 1);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let edges = explore_unique(&model, &fixed.shape, ShapeType::Edge).unwrap();
    let degenerate = edges
        .iter()
        .filter(|e| model.node(e).unwrap().data().as_edge().unwrap().degenerate)
        .count();
    assert_eq!((edges.len(), degenerate), (3, 2), "a seam and two poles");
    let exact = 4.0 * PI * 4.0;
    assert!((area(&model, &fixed.shape) - exact).abs() < 1e-9 * exact);

    let shell = model.add_shell(std::slice::from_ref(&ball)).unwrap();
    let solid = model.add_solid(&[shell]).unwrap();
    let fixed = ogeom::heal::fix_shape(&mut model, &solid, T).unwrap();
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let volume = volume_of(&model, &fixed.shape);
    let exact = 4.0 / 3.0 * PI * 8.0;
    assert!(
        (volume - exact).abs() < 1e-6 * exact,
        "{volume} against {exact}"
    );

    let torus = TorusSurface::new(Torus::new(Frame::WORLD, 5.0, 1.5, T).unwrap());
    let ring = ogeom::algo::make_natural_face(&mut model, torus.into())
        .unwrap()
        .shape;
    unmark(&mut model, &ring);
    let fixed = ogeom::heal::fix_shape(&mut model, &ring, T).unwrap();
    assert_eq!(fixed.report.natural_bounds_made, 1);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let exact = 4.0 * PI * PI * 5.0 * 1.5;
    assert!((area(&model, &fixed.shape) - exact).abs() < 1e-9 * exact);

    let plane = SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(Frame::WORLD)));
    let sheet = ogeom::algo::make_natural_face(&mut model, plane)
        .unwrap()
        .shape;
    unmark(&mut model, &sheet);
    let fixed = ogeom::heal::fix_shape(&mut model, &sheet, T).unwrap();
    assert_eq!(fixed.report.natural_bounds_made, 0);
    assert_eq!(fixed.report.natural_bounds_marked, 0);
    assert_eq!(no_wires(&fixed.report.after), 1);
}

/// A box's shell with its lid taken off is open round the lid's rim, and
/// with its floor off too, round two rims: each a closed contour of four
/// edges, walked as the faces beside it walk them. A closed shell has
/// none.
#[test]
fn free_edges_chain_into_the_rims_of_each_opening() {
    let mut model = Model::new();
    let shell = box_shell(&mut model, Point::ORIGIN, Point::new(10.0, 10.0, 10.0));
    assert!(ogeom::heal::free_bounds(&model, &shell).unwrap().is_empty());
    let faces = explore_unique(&model, &shell, ShapeType::Face).unwrap();
    let level = |model: &Model, face: &Shape| {
        let props = surface_properties(model, face, Deflection::default(), T).unwrap();
        props.centre.z
    };
    let lid = faces
        .iter()
        .find(|f| (level(&model, f) - 10.0).abs() < 1e-9)
        .unwrap();
    let floor = faces
        .iter()
        .find(|f| level(&model, f).abs() < 1e-9)
        .unwrap();
    for (open, rims) in [
        (vec![lid.clone()], 1),
        (vec![lid.clone(), floor.clone()], 2),
    ] {
        let kept: Vec<Shape> = faces
            .iter()
            .filter(|f| !open.iter().any(|o| o.is_same(f)))
            .cloned()
            .collect();
        let tube = model.add_shell(&kept).unwrap();
        let bounds = ogeom::heal::free_bounds(&model, &tube).unwrap();
        assert!(bounds.open.is_empty());
        assert_eq!(bounds.closed.len(), rims);
        for rim in &bounds.closed {
            assert_eq!(rim.len(), 4);
            for pair in 0..4 {
                let (_, end) = ogeom::algo::edge_vertices(&model, &rim[pair])
                    .unwrap()
                    .unwrap();
                let (start, _) = ogeom::algo::edge_vertices(&model, &rim[(pair + 1) % 4])
                    .unwrap()
                    .unwrap();
                assert!(end.is_same(&start), "head to tail");
            }
        }
    }
}

/// A face on the XY plane walking `ring`: each edge a straight line from
/// one point to another (its curve and its pcurve alike), bounded by the
/// vertices named by index, each vertex at its point with its tolerance.
/// The curves need not end at their vertices: the tolerance says how far
/// they may stand off.
fn plane_face(
    model: &mut Model,
    vertices: &[(Point, f64)],
    ring: &[(usize, usize, Point, Point)],
) -> Shape {
    let vertices: Vec<Shape> = vertices
        .iter()
        .map(|&(at, reach)| {
            let v = make_vertex(model, at).shape;
            model
                .widen(&v, ogeom::core::Tolerance::new(reach).unwrap())
                .unwrap();
            v
        })
        .collect();
    let surface = model
        .geometry_mut()
        .add_surface(SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(
            Frame::WORLD,
        ))));
    let mut edges = Vec::new();
    for &(a, b, from, to) in ring {
        let line = ogeom::geom::LineCurve::segment(from, to, T).unwrap().into();
        let curve = model.geometry_mut().add_curve(line);
        let range = (0.0, from.distance(to));
        let data = ogeom::topo::EdgeData::on_curve(curve, Location::identity(), range);
        let edge = model
            .add_edge(data, &[vertices[a].clone(), vertices[b].clone()])
            .unwrap();
        let chart =
            ogeom::geom::Line2d::segment(Point2::new(from.x, from.y), Point2::new(to.x, to.y), T)
                .unwrap();
        ogeom::algo::attach_pcurve(
            &mut *model,
            &edge,
            chart.into(),
            surface,
            Location::identity(),
            ogeom::geom::Curve2d::domain(&chart),
        )
        .unwrap();
        edges.push(edge);
    }
    let wire = model.add_wire(&edges).unwrap();
    model
        .add_face(FaceData::new(surface, Location::default()), &[wire])
        .unwrap()
}

fn crossings_in(model: &Model, face: &Shape) -> usize {
    ogeom::algo::chart_crossings(model, face, T).unwrap().len()
}

/// The 10 mm bow tie: (0, 0) to (10, 10) to (10, 0) to (0, 10) and back,
/// its first and third edges crossing at (5, 5).
fn bow_tie(model: &mut Model) -> Shape {
    let p = |x: f64, y: f64| Point::new(x, y, 0.0);
    let corners = [p(0.0, 0.0), p(10.0, 10.0), p(10.0, 0.0), p(0.0, 10.0)];
    let vertices: Vec<(Point, f64)> = corners.iter().map(|&c| (c, 1e-7)).collect();
    let ring: Vec<_> = (0..4)
        .map(|i| (i, (i + 1) % 4, corners[i], corners[(i + 1) % 4]))
        .collect();
    plane_face(model, &vertices, &ring)
}

/// A bow tie's boundary crosses itself at its middle, and `fix_shape`
/// splits the face, held in a compound, into its two
/// triangles, each wound as an outer boundary on the plane and facing up
/// as the face did: 25 mm² each, checking clean.
#[test]
fn a_bow_tie_face_is_split_into_its_two_loops() {
    let mut model = Model::new();
    let face = bow_tie(&mut model);
    assert_eq!(crossings_in(&model, &face), 1);
    let held = make_compound(&mut model, &[face]).unwrap().shape;

    let fixed = ogeom::heal::fix_shape(&mut model, &held, T).unwrap();
    assert_eq!(fixed.report.faces_split_at_crossings, 1);
    assert_eq!(fixed.report.crossings_left, 0);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let faces = explore_unique(&model, &fixed.shape, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 2);
    for face in &faces {
        let area = surface_properties(&model, face, Deflection::default(), T)
            .unwrap()
            .mass;
        assert!((area - 25.0).abs() < 1e-9, "{area}");
        let (_, normal) = ogeom::algo::face_normal(&model, face, T).unwrap();
        assert!(normal.z > 0.0, "faces up as the bow tie did: {normal:?}");
    }
}

/// A bow tie standing alone has no container to hold two faces: it is
/// left as it was, and counted as still crossing itself.
#[test]
fn a_bow_tie_face_alone_stays_reported() {
    let mut model = Model::new();
    let face = bow_tie(&mut model);
    let fixed = ogeom::heal::fix_shape(&mut model, &face, T).unwrap();
    assert_eq!(fixed.report.faces_split_at_crossings, 0);
    assert_eq!(fixed.report.crossings_left, 1);
}

/// `check` reports a face whose boundary crosses itself in its chart as
/// broken, once, where its edges cross: the bow tie at its middle. A
/// sliver within its edges' tolerances of their corner is no crossing,
/// and checks valid.
#[test]
fn check_reports_a_boundary_crossing_itself_as_broken() {
    let mut model = Model::new();
    let tie = bow_tie(&mut model);
    let diagnosis = check(&model, &tie, T).unwrap();
    let broken = diagnosis.of(ogeom::algo::Severity::Broken);
    assert_eq!(broken.len(), 1, "{diagnosis}");
    assert!(broken[0].at.is_same(&tie), "{diagnosis}");
    assert!(
        broken[0].what.contains("(5.000000, 5.000000, 0.000000)"),
        "{diagnosis}"
    );
    let sliver = sharp_corner(&mut model, 1e-3, 3e-4);
    let diagnosis = check(&model, &sliver, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
}

/// A thin triangle with its sharp corner at the origin, between (50, ±0.5),
/// whose edge out of the corner starts `lift` above it, inside the
/// corner's tolerance of 2 µm, and so crosses the edge into the corner
/// about 50 `lift` along it. Each edge states `stated`.
fn sharp_corner(model: &mut Model, lift: f64, stated: f64) -> Shape {
    let p = |x: f64, y: f64| Point::new(x, y, 0.0);
    let (top, corner, bottom) = (p(50.0, 0.5), p(0.0, 0.0), p(50.0, -0.5));
    let face = plane_face(
        model,
        &[(top, 1e-7), (corner, 2e-3), (bottom, 1e-7)],
        &[
            (0, 1, top, corner),
            (1, 2, p(0.0, lift), bottom),
            (2, 0, bottom, top),
        ],
    );
    for edge in explore_unique(model, &face, ShapeType::Edge).unwrap() {
        model
            .widen(&edge, ogeom::core::Tolerance::new(stated).unwrap())
            .unwrap();
    }
    face
}

/// The sharp corner's two edges cross 0.05 mm from it, past ten times its
/// tolerance, and the loop they close is 1 µm wide, far wider than twice
/// the edges' tolerances together:
/// both are cut back to the crossing and the corner moves there. What is
/// left is the triangle from the crossing, by the closed form of its area.
#[test]
fn edges_crossing_short_of_their_corner_are_cut_back_to_it() {
    let mut model = Model::new();
    let lift = 1e-3;
    let face = sharp_corner(&mut model, lift, 1e-7);
    assert_eq!(crossings_in(&model, &face), 1);

    let fixed = ogeom::heal::fix_shape(&mut model, &face, T).unwrap();
    assert_eq!(fixed.report.crossings_trimmed, 1);
    assert_eq!(fixed.report.crossings_left, 0);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    // y = x / 100 meets y = lift - (0.5 + lift) x / 50 at x = lift /
    // (0.02 + lift / 50).
    let x = lift / (0.02 + lift / 50.0);
    let crossing = Point2::new(x, x / 100.0);
    let (top, bottom) = (Point2::new(50.0, 0.5), Point2::new(50.0, -0.5));
    let expected = 0.5 * (bottom - crossing).cross(top - crossing);
    let area = surface_properties(&model, &fixed.shape, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((area - expected).abs() < 1e-9, "{area} against {expected}");
    let moved = explore_unique(&model, &fixed.shape, ShapeType::Vertex)
        .unwrap()
        .into_iter()
        .map(|v| model.node(&v).unwrap().data().as_vertex().unwrap().point)
        .any(|p| p.distance(Point::new(crossing.x, crossing.y, 0.0)) < 1e-9);
    assert!(moved, "the corner stands at the crossing");
}

/// The same corner, its edges each stating 0.3 µm: the 1 µm loop is
/// within twice their tolerances together, a sliver the images of two
/// curves meeting at the corner may close, and is no crossing.
#[test]
fn edges_crossing_within_their_tolerances_of_their_corner_are_no_crossing() {
    let mut model = Model::new();
    let face = sharp_corner(&mut model, 1e-3, 3e-4);
    assert_eq!(crossings_in(&model, &face), 0);
}

/// How many pairs of edges `check` names as running out and back.
fn out_and_back_found(diagnosis: &ogeom::algo::Diagnosis) -> usize {
    diagnosis
        .problems
        .iter()
        .filter(|p| p.what.contains("bounds nothing"))
        .count()
}

fn area_of(model: &Model, shape: &Shape) -> f64 {
    surface_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// A 4 mm square with an edge run from its corner (4, 0) to (3, 1) and
/// another straight back: `check` names the pair, and `fix_shape` strips
/// it, leaving the square.
#[test]
fn check_names_a_tail_and_fix_shape_strips_it() {
    let mut model = Model::new();
    let p = |x: f64, y: f64| (Point::new(x, y, 0.0), 1e-7);
    let at = |x: f64, y: f64| Point::new(x, y, 0.0);
    let face = plane_face(
        &mut model,
        &[
            p(0.0, 0.0),
            p(4.0, 0.0),
            p(3.0, 1.0),
            p(4.0, 4.0),
            p(0.0, 4.0),
        ],
        &[
            (0, 1, at(0.0, 0.0), at(4.0, 0.0)),
            (1, 2, at(4.0, 0.0), at(3.0, 1.0)),
            (2, 1, at(3.0, 1.0), at(4.0, 0.0)),
            (1, 3, at(4.0, 0.0), at(4.0, 4.0)),
            (3, 4, at(4.0, 4.0), at(0.0, 4.0)),
            (4, 0, at(0.0, 4.0), at(0.0, 0.0)),
        ],
    );
    let compound = make_compound(&mut model, &[face]).unwrap().shape;
    assert_eq!(out_and_back_found(&check(&model, &compound, T).unwrap()), 1);
    let fixed = ogeom::heal::fix_shape(&mut model, &compound, T).unwrap();
    assert_eq!(fixed.report.tails_stripped, 1);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    assert!((area_of(&model, &fixed.shape) - 16.0).abs() < 1e-9);
}

/// A 4 mm square whose edge up its right side, (4, 0) to (4, 4), is
/// followed by one back down it to (4, 2) before the boundary turns to
/// (0, 4): a notch. `check` names it; `fix_shape` cuts the side at
/// (4, 2) and drops the stretch above with the edge back, leaving the
/// quadrilateral (0, 0), (4, 0), (4, 2), (0, 4) of 12 mm². Walked the
/// other way round (mirrored, so it still faces up), the edge back is the
/// longer one, and it is cut instead.
#[test]
fn check_names_a_notch_and_fix_shape_cuts_it_off() {
    for mirrored in [false, true] {
        let mut model = Model::new();
        let x = |x: f64| if mirrored { -x } else { x };
        let p = |px: f64, py: f64| (Point::new(x(px), py, 0.0), 1e-7);
        let at = |px: f64, py: f64| Point::new(x(px), py, 0.0);
        let vertices = [
            p(0.0, 0.0),
            p(4.0, 0.0),
            p(4.0, 4.0),
            p(4.0, 2.0),
            p(0.0, 4.0),
        ];
        let ring = if mirrored {
            vec![
                (4, 3, at(0.0, 4.0), at(4.0, 2.0)),
                (3, 2, at(4.0, 2.0), at(4.0, 4.0)),
                (2, 1, at(4.0, 4.0), at(4.0, 0.0)),
                (1, 0, at(4.0, 0.0), at(0.0, 0.0)),
                (0, 4, at(0.0, 0.0), at(0.0, 4.0)),
            ]
        } else {
            vec![
                (0, 1, at(0.0, 0.0), at(4.0, 0.0)),
                (1, 2, at(4.0, 0.0), at(4.0, 4.0)),
                (2, 3, at(4.0, 4.0), at(4.0, 2.0)),
                (3, 4, at(4.0, 2.0), at(0.0, 4.0)),
                (4, 0, at(0.0, 4.0), at(0.0, 0.0)),
            ]
        };
        let face = plane_face(&mut model, &vertices, &ring);
        let compound = make_compound(&mut model, &[face]).unwrap().shape;
        let before = check(&model, &compound, T).unwrap();
        assert_eq!(out_and_back_found(&before), 1, "{before}");
        let fixed = ogeom::heal::fix_shape(&mut model, &compound, T).unwrap();
        assert_eq!(fixed.report.notches_cut, 1, "mirrored {mirrored}");
        assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
        assert_eq!(
            explore_unique(&model, &fixed.shape, ShapeType::Edge)
                .unwrap()
                .len(),
            4
        );
        assert!((area_of(&model, &fixed.shape) - 12.0).abs() < 1e-9);
    }
}

/// A plane face of tolerance `reach`, its outer ring the square from
/// (0, 0) to (10, 10), with a hole through each of `holes`' points in
/// turn, wound against the outer ring.
fn plane_face_with_holes(model: &mut Model, reach: f64, holes: &[Vec<Point2>]) -> Shape {
    let surface = model
        .geometry_mut()
        .add_surface(SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(
            Frame::WORLD,
        ))));
    let ring = |model: &mut Model, points: &[Point2]| -> Shape {
        let vertices: Vec<Shape> = points
            .iter()
            .map(|p| make_vertex(model, Point::new(p.x, p.y, 0.0)).shape)
            .collect();
        let mut edges = Vec::new();
        for k in 0..points.len() {
            let (from, to) = (points[k], points[(k + 1) % points.len()]);
            let (a, b) = (Point::new(from.x, from.y, 0.0), Point::new(to.x, to.y, 0.0));
            let line = ogeom::geom::LineCurve::segment(a, b, T).unwrap().into();
            let curve = model.geometry_mut().add_curve(line);
            let data =
                ogeom::topo::EdgeData::on_curve(curve, Location::identity(), (0.0, a.distance(b)));
            let edge = model
                .add_edge(
                    data,
                    &[
                        vertices[k].clone(),
                        vertices[(k + 1) % points.len()].clone(),
                    ],
                )
                .unwrap();
            let chart = ogeom::geom::Line2d::segment(from, to, T).unwrap();
            ogeom::algo::attach_pcurve(
                &mut *model,
                &edge,
                chart.into(),
                surface,
                Location::identity(),
                ogeom::geom::Curve2d::domain(&chart),
            )
            .unwrap();
            edges.push(edge);
        }
        model.add_wire(&edges).unwrap()
    };
    let square =
        [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)].map(|(x, y)| Point2::new(x, y));
    let mut wires = vec![ring(model, &square)];
    for hole in holes {
        wires.push(ring(model, hole));
    }
    let mut data = FaceData::new(surface, Location::identity());
    data.widen(ogeom::core::Tolerance::new(reach).unwrap());
    model.add_face(data, &wires).unwrap()
}

/// On a 10 mm square of tolerance 1e-3, a triangular hole 1 mm long and
/// 1e-4 mm high holds no point farther than 1e-3 from its boundary and
/// goes. A hole of 1 mm² stays, and so does a 0.01 mm square hole with a
/// slit 1 mm long and 1.5e-3 wide off it, which encloses less than the
/// tolerance times its boundary's length but whose middle stands 5e-3 from
/// its boundary. The face measures the square less the two holes kept.
#[test]
fn a_hole_within_the_face_s_tolerance_of_its_boundary_goes() {
    let pt = |x: f64, y: f64| Point2::new(x, y);
    let sliver = vec![pt(4.0, 5.0), pt(4.5, 5.0 + 1e-4), pt(5.0, 5.0)];
    let square = vec![pt(2.0, 2.0), pt(2.0, 3.0), pt(3.0, 3.0), pt(3.0, 2.0)];
    let flagged = vec![
        pt(7.0, 7.0),
        pt(7.0, 7.01),
        pt(7.01, 7.01),
        pt(7.01, 7.005 + 7.5e-4),
        pt(8.01, 7.005 + 7.5e-4),
        pt(8.01, 7.005 - 7.5e-4),
        pt(7.01, 7.005 - 7.5e-4),
        pt(7.01, 7.0),
    ];
    let mut model = Model::new();
    let face = plane_face_with_holes(&mut model, 1e-3, &[sliver, square, flagged]);
    let compound = make_compound(&mut model, &[face]).unwrap().shape;
    let fixed = ogeom::heal::fix_shape(&mut model, &compound, T).unwrap();
    assert_eq!(fixed.report.small_wires_removed, 1);
    assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    let flag_area = 1e-4 + 1.0 * 1.5e-3;
    let exact = 100.0 - 1.0 - flag_area;
    assert!(
        (area_of(&model, &fixed.shape) - exact).abs() < 1e-9,
        "{}",
        area_of(&model, &fixed.shape)
    );
    let faces = explore_unique(&model, &fixed.shape, ShapeType::Face).unwrap();
    assert_eq!(model.children_of(&faces[0]).unwrap().len(), 3);
}

/// A face of tolerance `reach` over the chart rectangle `outer` of
/// `surface`, every edge a chart segment carried onto the surface, with
/// `holes` (walked clockwise in the chart) cut from it.
fn chart_face_with_holes(
    model: &mut Model,
    surface: &SurfaceGeometry,
    outer: (Point2, Point2),
    reach: f64,
    holes: &[Vec<Point2>],
) -> Shape {
    use ogeom::geom::Surface as _;
    let id = model.geometry_mut().add_surface(surface.clone());
    let ring = |model: &mut Model, points: &[Point2]| -> Shape {
        let vertices: Vec<Shape> = points
            .iter()
            .map(|p| make_vertex(model, surface.point_at(p.x, p.y, T).unwrap()).shape)
            .collect();
        let mut edges = Vec::new();
        for k in 0..points.len() {
            let (from, to) = (points[k], points[(k + 1) % points.len()]);
            let chart = ogeom::geom::Line2d::segment(from, to, T).unwrap();
            let range = ogeom::geom::Curve2d::domain(&chart);
            let carried = ogeom::geom::Curve::OnSurface(Box::new(
                ogeom::geom::CurveOnSurface::new(chart.into(), surface.clone()),
            ));
            let curve = model.geometry_mut().add_curve(carried);
            let data = ogeom::topo::EdgeData::on_curve(curve, Location::identity(), range);
            let edge = model
                .add_edge(
                    data,
                    &[
                        vertices[k].clone(),
                        vertices[(k + 1) % points.len()].clone(),
                    ],
                )
                .unwrap();
            ogeom::algo::attach_pcurve(
                &mut *model,
                &edge,
                chart.into(),
                id,
                Location::identity(),
                range,
            )
            .unwrap();
            edges.push(edge);
        }
        model.add_wire(&edges).unwrap()
    };
    let (low, high) = outer;
    let corners = [
        low,
        Point2::new(high.x, low.y),
        high,
        Point2::new(low.x, high.y),
    ];
    let mut wires = vec![ring(model, &corners)];
    for hole in holes {
        wires.push(ring(model, hole));
    }
    let mut data = FaceData::new(id, Location::identity());
    data.widen(ogeom::core::Tolerance::new(reach).unwrap());
    model.add_face(data, &wires).unwrap()
}

/// The chart rectangle from `(u0, v0)` to `(u1, v1)`, walked clockwise.
fn chart_square(u0: f64, v0: f64, u1: f64, v1: f64) -> Vec<Point2> {
    vec![
        Point2::new(u0, v0),
        Point2::new(u0, v1),
        Point2::new(u1, v1),
        Point2::new(u1, v0),
    ]
}

/// On a cone, a sphere, a torus and a rational B-spline patch, a face of
/// tolerance 1e-3 holds a thin triangular hole about a millimetre long
/// and 1e-4 high, every point of it within 1e-3 of its boundary through
/// space, and a hole about a millimetre square. The thin hole goes and the
/// square stays: the face measures the closed form of its chart rectangle
/// less the square's. On the sphere a 0.01 mm square hole with a slit
/// 1 mm long and 1.5e-3 wide off it, enclosing less than the tolerance
/// times its boundary's length but whose middle stands 5e-3 from its
/// boundary, stays too.
#[test]
fn a_hole_within_the_face_s_tolerance_goes_on_curved_surfaces() {
    use ogeom::geom::{BSplineSurface, ConeSurface, SphereSurface, TorusSurface};
    use ogeom::math::{Cone, ControlGrid, KnotVector, Sphere, Torus, Weighted};
    let pt = |u: f64, v: f64| Point2::new(u, v);
    let sliver =
        |u: f64, v: f64, du: f64, dv: f64| vec![pt(u, v), pt(u + 0.5 * du, v + dv), pt(u + du, v)];

    // A cone of radius 5 at height 0 and half angle 0.3.
    let (r0, alpha) = (5.0, 0.3_f64);
    let cone: SurfaceGeometry =
        ConeSurface::new(Cone::new(Frame::WORLD, r0, alpha, T).unwrap(), (0.0, 10.0))
            .unwrap()
            .into();
    let cone_area = move |u0: f64, v0: f64, u1: f64, v1: f64| {
        let primitive = |v: f64| r0 * v + 0.5 * alpha.tan() * v * v;
        (u1 - u0) / alpha.cos() * (primitive(v1) - primitive(v0))
    };
    // A sphere of radius 10.
    let radius = 10.0;
    let sphere: SurfaceGeometry =
        SphereSurface::new(Sphere::new(Frame::WORLD, radius, T).unwrap()).into();
    let sphere_area = move |u0: f64, v0: f64, u1: f64, v1: f64| {
        radius * radius * (u1 - u0) * (v1.sin() - v0.sin())
    };
    // A torus of radii 10 and 3.
    let (major, minor) = (10.0, 3.0);
    let torus: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, major, minor, T).unwrap()).into();
    let torus_area = move |u0: f64, v0: f64, u1: f64, v1: f64| {
        minor * (u1 - u0) * (major * (v1 - v0) + minor * (v1.sin() - v0.sin()))
    };
    // A flat rational patch: x a rational quadratic of 0, 5, 10 weighted
    // 1, 2, 1 along u, y = 10 v.
    let mut net = Vec::new();
    for (x, w) in [(0.0, 1.0), (5.0, 2.0), (10.0, 1.0)] {
        for y in [0.0, 10.0] {
            net.push(Weighted::new(Point::new(x, y, 0.0), w, T).unwrap());
        }
    }
    let patch: SurfaceGeometry = BSplineSurface::rational(
        KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap(),
        KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap(),
        ControlGrid::new(net, 3, 2).unwrap(),
    )
    .unwrap()
    .into();
    let x = |u: f64| (20.0 * u * (1.0 - u) + 10.0 * u * u) / (1.0 + 2.0 * u * (1.0 - u));
    let patch_area = move |u0: f64, v0: f64, u1: f64, v1: f64| (x(u1) - x(u0)) * 10.0 * (v1 - v0);

    type Area = Box<dyn Fn(f64, f64, f64, f64) -> f64>;
    let flag = (2.0, 0.2, 1e-3, 0.1, 1.5e-4);
    let (fu, fv, side, slit, width) = flag;
    let flagged = vec![
        pt(fu, fv),
        pt(fu, fv + side),
        pt(fu + side, fv + side),
        pt(fu + side, fv + 0.5 * (side + width)),
        pt(fu + side + slit, fv + 0.5 * (side + width)),
        pt(fu + side + slit, fv + 0.5 * (side - width)),
        pt(fu + side, fv + 0.5 * (side - width)),
        pt(fu + side, fv),
    ];
    // A name, the surface, its chart rectangles' area, the face's
    // rectangle, the thin hole, the square hole, and whether the flagged
    // hole is cut too.
    type Case = (
        &'static str,
        SurfaceGeometry,
        Area,
        [f64; 4],
        Vec<Point2>,
        [f64; 4],
        bool,
    );
    let cases: [Case; 4] = [
        (
            "cone",
            cone,
            Box::new(cone_area),
            [0.5, 1.0, 2.0, 5.0],
            sliver(1.0, 3.0, 0.15, 1e-4),
            [1.5, 2.0, 1.6, 3.0],
            false,
        ),
        (
            "sphere",
            sphere,
            Box::new(sphere_area),
            [0.5, -0.5, 2.5, 0.5],
            sliver(1.0, 0.1, 0.1, 1e-5),
            [1.5, -0.3, 1.6, -0.2],
            true,
        ),
        (
            "torus",
            torus,
            Box::new(torus_area),
            [0.5, 0.5, 1.5, 1.5],
            sliver(0.8, 0.8, 0.05, 3e-5),
            [1.0, 1.0, 1.1, 1.1],
            false,
        ),
        (
            "patch",
            patch,
            Box::new(patch_area),
            [0.1, 0.1, 0.9, 0.9],
            sliver(0.5, 0.5, 0.02, 1e-5),
            [0.2, 0.2, 0.3, 0.3],
            false,
        ),
    ];
    for (name, surface, area, [u0, v0, u1, v1], thin, [a0, b0, a1, b1], with_flag) in cases {
        let mut model = Model::new();
        let mut holes = vec![thin, chart_square(a0, b0, a1, b1)];
        let mut exact = area(u0, v0, u1, v1) - area(a0, b0, a1, b1);
        if with_flag {
            holes.push(flagged.clone());
            exact -= area(fu, fv, fu + side, fv + side);
            exact -= area(
                fu + side,
                fv + 0.5 * (side - width),
                fu + side + slit,
                fv + 0.5 * (side + width),
            );
        }
        let kept = holes.len() - 1;
        let face =
            chart_face_with_holes(&mut model, &surface, (pt(u0, v0), pt(u1, v1)), 1e-3, &holes);
        let compound = make_compound(&mut model, &[face]).unwrap().shape;
        let fixed = ogeom::heal::fix_shape(&mut model, &compound, T).unwrap();
        assert_eq!(fixed.report.small_wires_removed, 1, "{name}");
        assert!(
            fixed.report.after.is_valid(),
            "{name}: {}",
            fixed.report.after
        );
        let measured = area_of(&model, &fixed.shape);
        assert!(
            (measured - exact).abs() < 1e-9 * exact,
            "{name}: {measured} against {exact}"
        );
        let faces = explore_unique(&model, &fixed.shape, ShapeType::Face).unwrap();
        assert_eq!(
            model.children_of(&faces[0]).unwrap().len(),
            1 + kept,
            "{name}"
        );
    }
}

/// A 10 mm box's faces with its lid lifted 1e-3 mm off the walls, held in
/// a shell and again loose in a compound: each of the lid's edges is free
/// beside a wall's, a gap the edges do not state. `fix_shape` states the
/// gap on each pair and sews them: the faces close and check clean, and
/// measure the box's 600 mm² to within the lift along the rim. A lid
/// lifted 2 mm, past the millimetre a fitted trim may stand off its
/// surface, stays apart, its faces as they were.
#[test]
fn free_edges_twinned_across_a_gap_are_sewn() {
    for (lift, sewn, loose) in [
        (1e-3, true, false),
        (2.0, false, false),
        (1e-3, true, true),
        (2.0, false, true),
    ] {
        let mut model = Model::new();
        let p = |x: f64, y: f64, z: f64| Point::new(x, y, z);
        let top = 10.0 + lift;
        let quads: [[Point; 4]; 6] = [
            [
                p(0., 0., 0.),
                p(0., 10., 0.),
                p(10., 10., 0.),
                p(10., 0., 0.),
            ],
            [
                p(0., 0., 0.),
                p(10., 0., 0.),
                p(10., 0., 10.),
                p(0., 0., 10.),
            ],
            [
                p(0., 10., 0.),
                p(0., 10., 10.),
                p(10., 10., 10.),
                p(10., 10., 0.),
            ],
            [
                p(0., 0., 0.),
                p(0., 0., 10.),
                p(0., 10., 10.),
                p(0., 10., 0.),
            ],
            [
                p(10., 0., 0.),
                p(10., 10., 0.),
                p(10., 10., 10.),
                p(10., 0., 10.),
            ],
            [
                p(0., 0., top),
                p(10., 0., top),
                p(10., 10., top),
                p(0., 10., top),
            ],
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
        let lid = faces.pop().unwrap();
        let walls = ogeom::algo::sew(&mut model, &faces, T).unwrap();
        assert_eq!(walls.shells.len(), 1);
        let mut all = explore_unique(&model, &walls.shells[0], ShapeType::Face).unwrap();
        all.push(lid);
        let open = if loose {
            make_compound(&mut model, &all).unwrap().shape
        } else {
            model.add_shell(&all).unwrap()
        };
        assert!(!check(&model, &open, T).unwrap().is_valid());
        let fixed = ogeom::heal::fix_shape(&mut model, &open, T).unwrap();
        if sewn {
            assert_eq!(fixed.report.free_gaps_stated, 4);
            assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
            assert!(fixed.report.free_bounds.is_empty());
        } else {
            assert_eq!(fixed.report.free_gaps_stated, 0);
            assert_eq!(fixed.report.free_bounds.closed.len(), 2);
        }
        // Sewn, the walls reach to whichever edge sewing keeps: at most
        // the lift times the rim's 40 mm more.
        let area = area_of(&model, &fixed.shape);
        let slack = if sewn { 40.0 * lift } else { 0.0 };
        assert!(
            (area - 600.0).abs() <= slack + 1e-9,
            "lift {lift}, loose {loose}: {area}"
        );
    }
}

/// `fix_shape` hands every sound primitive back as it is.
#[test]
fn sound_primitives_are_left_alone() {
    let mut model = Model::new();
    let solids = [
        make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T)
            .unwrap()
            .shape,
        ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 5.0, T)
            .unwrap()
            .shape,
        ogeom::algo::make_cone(&mut model, Frame::WORLD, 2.0, 0.0, 5.0, T)
            .unwrap()
            .shape,
        ogeom::algo::make_sphere(&mut model, Frame::WORLD, 2.0, T)
            .unwrap()
            .shape,
        ogeom::algo::make_torus(&mut model, Frame::WORLD, 5.0, 1.5, T)
            .unwrap()
            .shape,
    ];
    for solid in solids {
        assert!(check(&model, &solid, T).unwrap().is_valid());
        let fixed = ogeom::heal::fix_shape(&mut model, &solid, T).unwrap();
        assert!(fixed.shape.is_same(&solid), "the same node");
        assert!(fixed.report.after.is_valid(), "{}", fixed.report.after);
    }
}
