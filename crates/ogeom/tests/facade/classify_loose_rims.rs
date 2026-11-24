//! A solid whose faces meet along edges they lie only within the edges'
//! tolerance of, as a solid converted back from a mesh has them: a box
//! roofed by a cylinder whose two ends are closed by planes cut off along
//! straight chords of the roof. Between each chord and the roof's arc above
//! it the boundary leaves a gap as wide as the chord's tolerance, and a ray
//! out through that gap crosses no face. Such a ray cannot count crossings.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    Containment, attach_pcurve, make_edge_between, make_face, make_face_on, make_shell, make_solid,
    make_vertex, make_wire,
};
use ogeom::core::Tolerances;
use ogeom::geom::{
    Curve2d as _, Curve3d as _, CylinderSurface, Line2d, LineCurve, PlaneSurface, Surface as _,
};
use ogeom::math::{Cylinder, Direction, Frame, Plane, Point, Point2};
use ogeom::topo::{Location, Model, NodeData, Orientation, Shape};

const T: Tolerances = Tolerances::millimetres();

/// The roofed box: two by two by one, roofed by a cylinder of radius root
/// two along x through the box's top long edges, its axis on the floor.
/// The roof's ends are straight chords across the box's ends, which stand
/// up to 0.414 below the roof's arc.
fn roofed_box(model: &mut Model) -> Shape {
    let pi = core::f64::consts::PI;
    let roof = CylinderSurface::new(
        Cylinder::new(
            Frame::new(Point::new(0.0, 1.0, 0.0), Direction::X, Direction::Y, T).unwrap(),
            2f64.sqrt(),
            T,
        )
        .unwrap(),
        (-1.0, 3.0),
    )
    .unwrap();
    let lift = |u: f64, v: f64| roof.point_at(u, v, T).unwrap();
    let cylinder = model.geometry_mut().add_surface(roof.into());
    let at = |x: f64, y: f64, z: f64| Point::new(x, y, z);
    let corner = |model: &mut Model, x: f64, y: f64, z: f64| make_vertex(model, at(x, y, z)).shape;
    let [b00, b20, b22, b02] =
        [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)].map(|(x, y)| corner(model, x, y, 0.0));
    let [t00, t20, t22, t02] =
        [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)].map(|(x, y)| corner(model, x, y, 1.0));
    let point = |model: &Model, v: &Shape| model.node(v).unwrap().data().as_vertex().unwrap().point;
    let line = |model: &mut Model, from: &Shape, to: &Shape| {
        let (a, b) = (point(model, from), point(model, to));
        let curve = LineCurve::segment(a, b, T).unwrap();
        make_edge_between(model, curve.into(), (0.0, a.distance(b)), from, to, T)
            .unwrap()
            .shape
    };
    // An edge of the roof, its pcurve the chart segment between `from` and
    // `to`, its tolerance how far its line stands from the roof along it.
    let on_roof = |model: &mut Model, edge: &Shape, from: (f64, f64), to: (f64, f64)| {
        let (a, b) = (Point2::new(from.0, from.1), Point2::new(to.0, to.1));
        let chart = Line2d::segment(a, b, T).unwrap();
        let data = model.node(edge).unwrap().data().as_edge().unwrap().clone();
        let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            unreachable!("built on a line")
        };
        let curve = model.geometry().curve(*curve).unwrap().clone();
        let mut apart = 0.0_f64;
        for i in 0..=256 {
            let s = f64::from(i) / 256.0;
            let p = curve
                .point_at(range.0 + (range.1 - range.0) * s, T)
                .unwrap();
            let q = chart.point_at(a.distance(b) * s, T).unwrap();
            apart = apart.max(p.distance(lift(q.x, q.y)));
        }
        attach_pcurve(
            model,
            edge,
            chart.into(),
            cylinder,
            Location::identity(),
            (0.0, a.distance(b)),
        )
        .unwrap();
        if let Some(NodeData::Edge(data)) = model.node_mut(edge).map(|n| n.data_mut()) {
            data.tolerance = data.tolerance.widen_to(apart);
        }
    };
    let [floor_front, floor_right, floor_back, floor_left] =
        [(&b00, &b20), (&b20, &b22), (&b22, &b02), (&b02, &b00)].map(|(a, b)| line(model, a, b));
    let [up00, up20, up22, up02] =
        [(&b00, &t00), (&b20, &t20), (&b22, &t22), (&b02, &t02)].map(|(a, b)| line(model, a, b));
    let (near, far) = (0.75 * pi, 0.25 * pi);
    let eave_front = line(model, &t00, &t20);
    on_roof(model, &eave_front, (near, 0.0), (near, 2.0));
    let eave_back = line(model, &t02, &t22);
    on_roof(model, &eave_back, (far, 0.0), (far, 2.0));
    let chord_left = line(model, &t00, &t02);
    on_roof(model, &chord_left, (near, 0.0), (far, 0.0));
    let chord_right = line(model, &t20, &t22);
    on_roof(model, &chord_right, (near, 2.0), (far, 2.0));

    let back = |e: &Shape| e.oriented(Orientation::Reversed);
    let ring = |model: &mut Model, edges: &[Shape]| make_wire(model, edges, T).unwrap().shape;
    let plane = |origin: Point, normal: Direction| -> ogeom::geom::SurfaceGeometry {
        PlaneSurface::new(Plane::through(origin, normal)).into()
    };
    let mut faces = Vec::new();
    let walls = [
        (
            vec![
                floor_front.clone(),
                floor_right.clone(),
                floor_back.clone(),
                floor_left.clone(),
            ],
            plane(at(0.0, 0.0, 0.0), -Direction::Z),
        ),
        (
            vec![
                floor_front.clone(),
                up20.clone(),
                back(&eave_front),
                back(&up00),
            ],
            plane(at(0.0, 0.0, 0.0), -Direction::Y),
        ),
        (
            vec![
                back(&floor_back),
                up22.clone(),
                back(&eave_back),
                back(&up02),
            ],
            plane(at(0.0, 2.0, 0.0), Direction::Y),
        ),
        (
            vec![
                back(&floor_left),
                up02.clone(),
                back(&chord_left),
                back(&up00),
            ],
            plane(at(0.0, 0.0, 0.0), -Direction::X),
        ),
        (
            vec![
                floor_right.clone(),
                up22.clone(),
                back(&chord_right),
                back(&up20),
            ],
            plane(at(2.0, 0.0, 0.0), Direction::X),
        ),
    ];
    for (edges, surface) in walls {
        let wire = ring(model, &edges);
        let face = make_face(model, surface, &[wire], T).unwrap().shape;
        ogeom::algo::build::trimmed_where_bare(model, &face, T).unwrap();
        faces.push(face);
    }
    // Counter-clockwise in the roof's chart, whose normal points away from
    // its axis: out of the box.
    let wire = ring(
        model,
        &[
            back(&chord_left),
            eave_front.clone(),
            chord_right.clone(),
            back(&eave_back),
        ],
    );
    faces.push(make_face_on(model, cylinder, &[wire], T).unwrap().shape);
    let shell = make_shell(model, &faces).unwrap().shape;
    make_solid(model, &[shell]).unwrap().shape
}

/// A point inside the box whose first ray leaves through the gap under
/// the roof's far end: it crosses no face there, and counted it reads the
/// point outside. Crossing a face within its edges' tolerance of its trim,
/// where the gap is, the ray cannot count, and a ray clear of every gap
/// reads the point inside.
#[test]
fn a_ray_out_through_a_gap_its_edges_tolerance_covers_is_not_counted() {
    let mut model = Model::new();
    let solid = roofed_box(&mut model);
    let inside = Point::new(1.5, 0.5, 0.7);
    assert_eq!(
        ogeom::algo::classify_in_solid_exact(&model, &solid, inside, T).unwrap(),
        Containment::In,
        "a point inside the box, under its roof"
    );
}
