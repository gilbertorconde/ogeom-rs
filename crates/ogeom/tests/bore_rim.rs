//! The rim of a bore cut into a planar face, and the top edge of a padded
//! drum: a circle where a plane meets a cylinder, rounded or bevelled. A
//! circular profile swept against its own axis, as a tool pushed down
//! through a block is, makes a cylinder, so the rim is the circle it is and
//! the circular blend takes it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, make_edge, make_edge_between, make_face_with_pcurves, make_prism, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve, Curve3d, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{EdgeRepr, Model, Shape, ShapeType, VertexData, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn volume(model: &Model, shape: &Shape) -> f64 {
    let v = volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T).unwrap();
    assert_eq!(v.deflection, 0.0, "integrated on the exact surfaces");
    v.mass
}

/// How many faces of each kind: planes, cylinders, cones, tori.
fn kinds(model: &Model, shape: &Shape) -> [usize; 4] {
    let mut out = [0; 4];
    for face in explore_unique(model, shape, ShapeType::Face).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        match model.geometry().surface(data.surface).unwrap() {
            SurfaceGeometry::Plane(_) => out[0] += 1,
            SurfaceGeometry::Cylinder(_) => out[1] += 1,
            SurfaceGeometry::Cone(_) => out[2] += 1,
            SurfaceGeometry::Torus(_) => out[3] += 1,
            other => panic!("an unexpected surface: {other:?}"),
        }
    }
    out
}

/// The circular edge of `shape` passing nearest `near`, placement applied.
fn circle_near(model: &Model, shape: &Shape, near: Point) -> Shape {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter_map(|edge| {
            let data = model.node(&edge).unwrap().data().as_edge().unwrap();
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                return None;
            };
            let curve = model.geometry().curve(*curve).unwrap();
            if !matches!(curve, Curve::Circle(_)) {
                return None;
            }
            let placed = edge.transform(model.datums()).unwrap();
            let gap = (0..=32)
                .map(|k| {
                    let t = range.0 + (range.1 - range.0) * f64::from(k) / 32.0;
                    placed.apply(curve.point_at(t, T).unwrap()).distance(near)
                })
                .fold(f64::INFINITY, f64::min);
            Some((gap, edge))
        })
        .min_by(|a, b| a.0.total_cmp(&b.0))
        .expect("a circular edge")
        .1
}

/// A 40 by 30 block 12 high, a disc of radius 6 on its top face at
/// (20, 15), and the disc's face pushed down through the block or up off
/// it.
fn parts(model: &mut Model) -> (Shape, Shape) {
    let corners = [(0.0, 0.0), (40.0, 0.0), (40.0, 30.0), (0.0, 30.0)];
    let points: Vec<Point> = corners
        .iter()
        .map(|(x, y)| Point::new(*x, *y, 0.0))
        .collect();
    let vertices: Vec<Shape> = points
        .iter()
        .map(|p| model.add_vertex(VertexData::new(*p)))
        .collect();
    let edges: Vec<Shape> = (0..4)
        .map(|i| {
            let j = (i + 1) % 4;
            let line = LineCurve::segment(points[i], points[j], T).unwrap();
            let range = line.domain();
            make_edge_between(
                model,
                Curve::Line(line),
                range,
                &vertices[i],
                &vertices[j],
                T,
            )
            .unwrap()
            .shape
        })
        .collect();
    let floor = PlaneSurface::over(Plane::XY, (-100.0, 100.0), (-100.0, 100.0)).unwrap();
    let base = make_face_with_pcurves(model, floor.into(), &[edges], T)
        .unwrap()
        .shape;
    let block = make_prism(model, &base, Vector::new(0.0, 0.0, 12.0), T)
        .unwrap()
        .shape;
    let top = Frame::new(Point::new(20.0, 15.0, 12.0), Direction::Z, Direction::X, T).unwrap();
    let circle = Curve::Circle(CircleCurve::new(Circle::new(top, 6.0, T).unwrap()));
    let range = circle.domain();
    let rim = make_edge(model, circle, range, T).unwrap().shape;
    let lid = PlaneSurface::over(Plane::new(top), (-10.0, 10.0), (-10.0, 10.0)).unwrap();
    let disc = make_face_with_pcurves(model, lid.into(), &[vec![rim]], T)
        .unwrap()
        .shape;
    (block, disc)
}

/// The rim of a bore pushed down through the block rounds with one torus
/// tangent to the top face and the bore; by Pappus the fillet removes a
/// ring of area rho^2 (1 - pi/4) whose centroid stands out from the rim.
#[test]
fn a_bore_rim_rounds_and_bevels() {
    let mut model = Model::new();
    let (block, disc) = parts(&mut model);
    let tool = make_prism(&mut model, &disc, Vector::new(0.0, 0.0, -14.0), T)
        .unwrap()
        .shape;
    assert_eq!(
        kinds(&model, &tool),
        [2, 1, 0, 0],
        "the tool's wall is a cylinder"
    );
    assert!((volume(&model, &tool) - PI * 36.0 * 14.0).abs() < 1e-9);
    let bored = ogeom::boolean::cut(&mut model, &block, &tool, T)
        .unwrap()
        .shape;
    let rim = circle_near(&model, &bored, Point::new(26.0, 15.0, 12.0));

    let (radius, rho) = (6.0, 2.0);
    let ring = rho * rho * (1.0 - PI / 4.0);
    let out = rho * (10.0 - 3.0 * PI) / (3.0 * (4.0 - PI));
    let bore = 40.0 * 30.0 * 12.0 - PI * radius * radius * 12.0;
    let want = bore - 2.0 * PI * (radius + out) * ring;
    let rounded =
        ogeom::fillet::fillet_edges(&mut model, &bored, std::slice::from_ref(&rim), rho, T)
            .unwrap()
            .shape;
    assert!(check(&model, &rounded, T).unwrap().is_valid());
    assert_eq!(kinds(&model, &rounded), [6, 1, 0, 1]);
    let got = volume(&model, &rounded);
    assert!((got - want).abs() < want * 1e-9, "{got} against {want}");

    // A symmetric bevel of 1: a right triangle of area a half, its
    // centroid a third of the way out from the rim.
    let bevel = ogeom::fillet::Chamfer::Symmetric(1.0);
    let bevelled = ogeom::fillet::chamfer_edges_with(&mut model, &bored, &[(rim, bevel)], T)
        .unwrap()
        .shape;
    assert!(check(&model, &bevelled, T).unwrap().is_valid());
    assert_eq!(kinds(&model, &bevelled), [6, 1, 1, 0]);
    let want = bore - 2.0 * PI * (radius + 1.0 / 3.0) * 0.5;
    let got = volume(&model, &bevelled);
    assert!((got - want).abs() < want * 1e-9, "{got} against {want}");
}

/// The top edge of a drum padded up off the block's face rounds on the
/// outside, the ring's centroid standing in from the edge.
#[test]
fn a_padded_drum_top_rounds() {
    let mut model = Model::new();
    let (_, disc) = parts(&mut model);
    let pad = make_prism(&mut model, &disc, Vector::new(0.0, 0.0, 8.0), T)
        .unwrap()
        .shape;
    let edge = circle_near(&model, &pad, Point::new(26.0, 15.0, 20.0));
    let (radius, rho) = (6.0, 2.0);
    let rounded = ogeom::fillet::fillet_edges(&mut model, &pad, &[edge], rho, T)
        .unwrap()
        .shape;
    assert!(check(&model, &rounded, T).unwrap().is_valid());
    assert_eq!(kinds(&model, &rounded), [2, 1, 0, 1]);
    let ring = rho * rho * (1.0 - PI / 4.0);
    let inward = rho * (10.0 - 3.0 * PI) / (3.0 * (4.0 - PI));
    let want = PI * radius * radius * 8.0 - 2.0 * PI * (radius - inward) * ring;
    let got = volume(&model, &rounded);
    assert!((got - want).abs() < want * 1e-9, "{got} against {want}");
}
