//! §10's tail: what a blend achieved, measured; blends between faces that
//! share no edge; edges whose envelope has no closed form; and the corner
//! where three of them meet.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The edge of `shape` whose midpoint is nearest `near`.
fn edge_near(model: &Model, shape: &Shape, near: Point) -> Shape {
    use ogeom::geom::Curve3d as _;
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| {
            let mid = |e: &Shape| {
                let data = model.node(e).unwrap().data().as_edge().unwrap();
                let ogeom::topo::EdgeRepr::Curve3d { curve, range, .. } = data.curve3d().unwrap()
                else {
                    unreachable!()
                };
                model
                    .geometry()
                    .curve(*curve)
                    .unwrap()
                    .point_at(f64::midpoint(range.0, range.1), T)
                    .unwrap()
                    .distance(near)
            };
            mid(a)
                .partial_cmp(&mid(b))
                .unwrap_or(core::cmp::Ordering::Equal)
        })
        .expect("some edge")
}

/// The planar face of `shape` whose plane passes through `on` and whose own
/// vertices bracket it.
fn planar_face_at(model: &Model, shape: &Shape, on: Point) -> Shape {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let ogeom::topo::NodeData::Face(data) = model.node(f).unwrap().data() else {
                return false;
            };
            let Some(ogeom::geom::SurfaceGeometry::Plane(plane)) =
                model.geometry().surface(data.surface)
            else {
                return false;
            };
            if plane.plane().distance_to(on).abs() > 1e-9 {
                return false;
            }
            let mut bound = ogeom::math::Aabb::EMPTY;
            for v in explore_unique(model, f, ShapeType::Vertex).unwrap() {
                bound = bound.with_point(model.node(&v).unwrap().data().as_vertex().unwrap().point);
            }
            bound.expanded(1e-6).contains(on)
        })
        .expect("a planar face there")
}

#[test]
fn a_fillet_reports_its_own_tangency_instead_of_claiming_it() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 12.0), T)
        .unwrap()
        .shape;
    let edge = edge_near(&model, &block, Point::new(20.0, 0.0, 12.0));
    let blended = ogeom::fillet::fillet_edge(&mut model, &block, &edge, 2.0, T)
        .unwrap()
        .shape;

    // The blend is the one cylindrical face on the result.
    let blend = explore_unique(&model, &blended, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let ogeom::topo::NodeData::Face(data) = model.node(f).unwrap().data() else {
                return false;
            };
            matches!(
                model.geometry().surface(data.surface),
                Some(ogeom::geom::SurfaceGeometry::Cylinder(_))
            )
        })
        .expect("the rolling ball left a cylinder");

    let contacts = ogeom::fillet::analyse_blend(&model, &blended, &blend, 9, T).unwrap();
    assert_eq!(contacts.len(), 4, "two tangency edges and two end caps");
    // The two long edges are the tangency lines: smooth to rounding. The
    // two ends are the cap arcs, where the blend meets a face it is *not*
    // tangent to — a right angle, and it should say so.
    let mut smooth = 0;
    let mut square = 0;
    for contact in &contacts {
        assert!(
            contact.gap < 1e-9,
            "the shared edge lies on both surfaces: {}",
            contact.gap
        );
        if contact.tangency_error < 1e-9 {
            smooth += 1;
        } else if (contact.tangency_error - core::f64::consts::FRAC_PI_2).abs() < 1e-9 {
            square += 1;
        }
    }
    assert_eq!(
        (smooth, square),
        (2, 2),
        "two tangent joins, two square ones: {contacts:?}"
    );
}

#[test]
fn a_blend_bridges_two_faces_that_share_no_edge() {
    // A step: a tall block and a low one side by side, their vertical wall
    // and horizontal lid meeting at no edge at all. The rolling ball still
    // has a seat — it touches both — and the blend is the fillet that seat
    // implies.
    let mut model = Model::new();
    let tall = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    let low = ogeom::algo::make_box(
        &mut model,
        Frame::new(Point::new(10.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap(),
        (20.0, 20.0, 10.0),
        T,
    )
    .unwrap()
    .shape;
    let step = ogeom::boolean::fuse(&mut model, &tall, &low, T)
        .unwrap()
        .shape;

    let wall = planar_face_at(&model, &step, Point::new(10.0, 10.0, 15.0));
    let lid = planar_face_at(&model, &step, Point::new(20.0, 10.0, 10.0));

    let blended = ogeom::fillet::blend_faces(&mut model, &step, &wall, &lid, 4.0, T)
        .unwrap()
        .shape;
    let volume =
        ogeom::algo::volume_properties(&model, &blended, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    // The step is 10*20*20 + 20*20*10 = 8000, and its inner corner is
    // concave: the ball rolls in the notch, so the blend *fills* it with
    // what a square corner would have held minus the quarter disc,
    // (r^2 - pi r^2 / 4), along the 20 of run.
    let r: f64 = 4.0;
    let filled = r.mul_add(r, -(core::f64::consts::PI * r * r / 4.0)) * 20.0;
    assert!(
        (volume - (8000.0 + filled)).abs() < 8000.0 * 2e-3,
        "the notch is filled, not cut: {volume} against {}",
        8000.0 + filled
    );
}
