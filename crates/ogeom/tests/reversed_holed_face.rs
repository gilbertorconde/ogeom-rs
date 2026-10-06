//! A face used reversed and carrying a hole: the base of a plate with a
//! bore, extruded upward, faces down and is its profile turned over.
//!
//! Walking a reversed face lists its wires backward, the hole before the
//! outer wire; the face stores the outer wire first. Every operation here
//! is measured against the plate's own numbers, so one that took the walk's
//! first wire for the outer one answers with the hole's region instead.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{Containment, make_edge, make_face_with_pcurves, make_polygon, make_prism};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, PlaneSurface};
use ogeom::math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Filter, Model, Orientation, Shape, ShapeType, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;
// The plate: 30 by 20, 5 thick, a bore of radius 4 through (15, 10).
const W: f64 = 30.0;
const H: f64 = 20.0;
const THICK: f64 = 5.0;
const R: f64 = 4.0;
const CENTRE: (f64, f64) = (15.0, 10.0);

fn fine() -> Deflection {
    Deflection::with_chord(1e-3).unwrap()
}

fn net_area() -> f64 {
    W * H - PI * R * R
}

/// The profile at z = 0: the outer square anticlockwise about +Z, the
/// bore's circle clockwise. With `turned`, the square's edges run
/// clockwise and the wire uses each of them reversed.
fn profile(model: &mut Model, turned: bool) -> Shape {
    let mut corners = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(W, 0.0, 0.0),
        Point::new(W, H, 0.0),
        Point::new(0.0, H, 0.0),
    ];
    if turned {
        corners.reverse();
    }
    let square = make_polygon(model, &corners, true, T).unwrap().shape;
    let mut outer = explore(model, &square, Filter::OfType(ShapeType::Edge)).unwrap();
    if turned {
        outer = outer.iter().rev().map(Shape::reversed).collect();
    }
    let frame = Frame::new(
        Point::new(CENTRE.0, CENTRE.1, 0.0),
        -Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let circle = Circle::new(frame, R, T).unwrap();
    let bore = make_edge(model, CircleCurve::new(circle).into(), (0.0, 2.0 * PI), T)
        .unwrap()
        .shape;
    let plane = PlaneSurface::new(Plane::through(Point::ORIGIN, Direction::Z));
    make_face_with_pcurves(model, plane.into(), &[outer, vec![bore]], T)
        .unwrap()
        .shape
}

/// The plate, and its base: the face at z = 0, used reversed.
fn plate(model: &mut Model) -> (Shape, Shape) {
    plate_of(model, false)
}

fn plate_of(model: &mut Model, turned: bool) -> (Shape, Shape) {
    let face = profile(model, turned);
    let solid = make_prism(model, &face, Vector::new(0.0, 0.0, THICK), T)
        .unwrap()
        .shape;
    let base = explore_unique(model, &solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| f.orientation() == Orientation::Reversed && wire_count(model, f) == 2)
        .expect("the base is used reversed and keeps the bore");
    (solid, base)
}

fn wire_count(model: &Model, face: &Shape) -> usize {
    model.children_of(face).unwrap().len()
}

/// How wide a wire's vertices spread: the square's corners span the plate,
/// the bore's single vertex nothing.
fn spread(model: &Model, wire: &Shape) -> f64 {
    let points: Vec<Point> = explore_unique(model, wire, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| {
            let data = model.node(v).unwrap().data().as_vertex().unwrap();
            v.transform(model.datums()).unwrap().apply(data.point)
        })
        .collect();
    let mut widest = 0.0_f64;
    for a in &points {
        for b in &points {
            widest = widest.max(a.distance(*b));
        }
    }
    widest
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, fine(), T)
        .unwrap()
        .mass
}

fn area(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::surface_properties(model, shape, fine(), T)
        .unwrap()
        .mass
}

/// The premise: the face stores the outer wire first, the walk lists the
/// bore first, and the outer-wire accessor answers with the square.
#[test]
fn the_reversed_base_walks_its_bore_first_and_stores_its_outer_wire_first() {
    let mut model = Model::new();
    let (_, base) = plate(&mut model);
    let walked = model.ordered_children_of(&base).unwrap();
    let stored = model.children_of(&base).unwrap();
    assert!(
        spread(&model, &walked[0]) < 1e-9,
        "the walk starts at the bore"
    );
    assert!(
        spread(&model, &stored[0]) > W,
        "the store starts at the square"
    );
    let outer = model.outer_wire(&base).unwrap().expect("a bounded face");
    assert_eq!(outer.node(), stored[0].node());
    assert_eq!(outer.orientation(), Orientation::Reversed);
}

/// Inside the bore is off the face and out of the solid; beside it is on
/// the face and in the solid.
#[test]
fn a_point_in_the_bore_is_out_and_one_beside_it_in() {
    let mut model = Model::new();
    let (solid, base) = plate(&mut model);
    let face = |p| ogeom::algo::classify_on_face(&model, &base, p, fine(), T).unwrap();
    assert_eq!(face(Point::new(CENTRE.0, CENTRE.1, 0.0)), Containment::Out);
    assert_eq!(
        face(Point::new(CENTRE.0 + 3.0, CENTRE.1, 0.0)),
        Containment::Out
    );
    assert_eq!(face(Point::new(5.0, 5.0, 0.0)), Containment::In);
    assert_eq!(
        face(Point::new(CENTRE.0 + 4.5, CENTRE.1, 0.0)),
        Containment::In
    );
    let body = |p| ogeom::algo::classify_in_solid(&model, &solid, p, fine(), T).unwrap();
    assert_eq!(body(Point::new(CENTRE.0, CENTRE.1, 2.5)), Containment::Out);
    assert_eq!(body(Point::new(5.0, 5.0, 2.5)), Containment::In);
    // On the base's own plane, where only its rings decide.
    let exact = |p| ogeom::algo::classify_in_solid_exact(&model, &solid, p, T).unwrap();
    assert_eq!(exact(Point::new(CENTRE.0, CENTRE.1, 0.0)), Containment::Out);
    assert_eq!(
        exact(Point::new(CENTRE.0 + 3.0, CENTRE.1, 0.0)),
        Containment::Out
    );
    assert_eq!(exact(Point::new(5.0, 5.0, 0.0)), Containment::On);
    assert_eq!(exact(Point::new(CENTRE.0, CENTRE.1, 2.5)), Containment::Out);
}

/// The base's mesh covers the plate less the bore, and no triangle stands
/// in the bore.
#[test]
fn the_reversed_base_meshes_around_its_bore() {
    let mut model = Model::new();
    let (_, base) = plate(&mut model);
    let mesh = ogeom::mesh::triangulate_face(&model, &base, fine(), T).unwrap();
    let mut total = 0.0;
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        total += 0.5 * (b - a).cross(c - a).magnitude();
        let mid = Point::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
            (a.z + b.z + c.z) / 3.0,
        );
        let off = (mid.x - CENTRE.0).hypot(mid.y - CENTRE.1);
        assert!(
            off > R - 0.01,
            "a triangle in the bore, {off} from its axis"
        );
    }
    assert!(
        (total - net_area()).abs() < 1e-2,
        "{total} against {}",
        net_area()
    );
}

/// Area and volume: the plate less the bore.
#[test]
fn the_reversed_base_measures_the_plate_less_the_bore() {
    let mut model = Model::new();
    let (solid, base) = plate(&mut model);
    let a = area(&model, &base);
    assert!((a - net_area()).abs() < 1e-3, "{a} against {}", net_area());
    let v = volume(&model, &solid);
    let want = net_area() * THICK;
    assert!((v - want).abs() < 1e-3, "{v} against {want}");
}

/// Grown by a half: every plane moves out a half and the bore closes in
/// by one, and the bore stays empty.
#[test]
fn the_plate_offsets_with_its_bore_narrowed() {
    let mut model = Model::new();
    let (solid, _) = plate(&mut model);
    let d = 0.5;
    let grown = ogeom::offset::shape::offset_shape(&mut model, &solid, d, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &grown, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let v = volume(&model, &grown);
    let want = (W + 2.0 * d) * (H + 2.0 * d) * (THICK + 2.0 * d)
        - PI * (R - d).powi(2) * (THICK + 2.0 * d);
    assert!((v - want).abs() < 1e-2, "{v} against {want}");
    let at = Point::new(CENTRE.0, CENTRE.1, THICK / 2.0);
    assert_eq!(
        ogeom::algo::classify_in_solid(&model, &grown, at, fine(), T).unwrap(),
        Containment::Out
    );
}

/// A drill narrower than the bore, down its axis, meets no material.
#[test]
fn a_drill_down_the_bore_meets_no_material() {
    let mut model = Model::new();
    let (solid, _) = plate(&mut model);
    let frame = Frame::new(
        Point::new(CENTRE.0, CENTRE.1, -1.0),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, frame, R - 1.0, THICK + 2.0, T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(&mut model, &solid, &drill, T)
        .unwrap()
        .shape;
    let want = net_area() * THICK;
    let v = volume(&model, &cut);
    assert!((v - want).abs() < 1e-3, "{v} against {want}");
    let common = ogeom::boolean::common(&mut model, &solid, &drill, T)
        .unwrap()
        .shape;
    let solids = explore_unique(&model, &common, ShapeType::Solid).unwrap();
    let shared: f64 = solids.iter().map(|s| volume(&model, s)).sum();
    assert!(
        shared.abs() < 1e-6,
        "the drill shares {shared} with the plate"
    );
}
