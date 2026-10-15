//! Lofts and pipes closing to a point through several sections.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/walks.rs"]
mod walks;

use ogeom::algo::{make_edge, make_polygon, make_wire, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::CircleCurve;
use ogeom::math::{Circle, Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, VertexData, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

fn rectangle(model: &mut Model, (x0, y0): (f64, f64), (x1, y1): (f64, f64), z: f64) -> Shape {
    let corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)].map(|(x, y)| Point::new(x, y, z));
    make_polygon(model, &corners, true, T).unwrap().shape
}

fn circle(model: &mut Model, radius: f64, z: f64) -> Shape {
    let frame = Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X, T).unwrap();
    let curve = CircleCurve::new(Circle::new(frame, radius, T).unwrap());
    let edge = make_edge(model, curve.into(), (0.0, 2.0 * PI), T)
        .unwrap()
        .shape;
    make_wire(model, &[edge], T).unwrap().shape
}

fn has_vertex_at(model: &Model, shape: &Shape, at: Point) -> bool {
    explore_unique(model, shape, ShapeType::Vertex)
        .unwrap()
        .iter()
        .any(|v| {
            let p = model.node(v).unwrap().data().as_vertex().unwrap().point;
            p.distance(at) <= 1e-9
        })
}

/// Two rectangles and a point: the skin passes through the middle section
/// and closes on the point. The middle section is one of the skin's own
/// rows, so it is probed by classifying points either side of it.
#[test]
fn a_skinned_loft_of_rectangles_closes_to_a_point() {
    let mut model = Model::new();
    let base = rectangle(&mut model, (0.0, 0.0), (10.0, 20.0), 0.0);
    let middle = rectangle(&mut model, (1.0, 2.0), (9.0, 18.0), 5.0);
    let apex = model.add_vertex(VertexData::new(Point::new(5.0, 10.0, 10.0)));
    let loft = ogeom::offset::make_loft_skinned(&mut model, &[base, middle, apex], 1e-3, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &loft, T).unwrap().is_valid());
    assert_eq!(walks::edges_walked_one_way(&model, &loft), 0);
    assert!(has_vertex_at(&model, &loft, Point::new(5.0, 10.0, 10.0)));
    // The skin passes through the middle rectangle: just inside its corners
    // and sides is material, just outside is not.
    let fine = Deflection::with_chord(1e-3).unwrap();
    let side = |at: Point| ogeom::algo::classify_in_solid(&model, &loft, at, fine, T).unwrap();
    for (x, y) in [(1.0_f64, 2.0_f64), (9.0, 18.0), (5.0, 2.0), (9.0, 10.0)] {
        let inward = Point::new(
            x + (5.0 - x).signum() * 0.02,
            y + (10.0 - y).signum() * 0.02,
            5.0,
        );
        let outward = Point::new(
            x - (5.0 - x).signum() * 0.02,
            y - (10.0 - y).signum() * 0.02,
            5.0,
        );
        assert_eq!(side(inward), ogeom::algo::Containment::In, "{inward:?}");
        assert_eq!(side(outward), ogeom::algo::Containment::Out, "{outward:?}");
    }
}

/// Two coaxial circles and a point on their axis: a solid of revolution
/// about the axis, closing exactly on the point. Its volume is the one its
/// own meridian (read off the seam, which runs along y = 0) sweeps.
#[test]
fn a_skinned_loft_of_circles_closes_to_a_point_on_the_axis() {
    let mut model = Model::new();
    let base = circle(&mut model, 5.0, 0.0);
    let middle = circle(&mut model, 3.0, 5.0);
    let apex = model.add_vertex(VertexData::new(Point::new(0.0, 0.0, 10.0)));
    let loft = ogeom::offset::make_loft_skinned(&mut model, &[base, middle, apex], 1e-3, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &loft, T).unwrap().is_valid());
    assert_eq!(walks::edges_walked_one_way(&model, &loft), 0);
    assert!(has_vertex_at(&model, &loft, Point::new(0.0, 0.0, 10.0)));
    let mesh =
        ogeom::mesh::triangulate(&model, &loft, Deflection::with_chord(1e-3).unwrap(), T).unwrap();
    let mut meridian: Vec<(f64, f64)> = mesh
        .positions
        .iter()
        .filter(|p| p.y.abs() < 1e-9 && p.x > 0.0 && p.z > 1e-9)
        .map(|p| (p.z, p.x))
        .collect();
    meridian.sort_by(|a, b| a.0.total_cmp(&b.0));
    meridian.insert(0, (0.0, 5.0));
    meridian.push((10.0, 0.0));
    let swept: f64 = meridian
        .windows(2)
        .map(|w| PI * (w[1].0 - w[0].0) * (w[0].1 * w[0].1 + w[1].1 * w[1].1) / 2.0)
        .sum();
    let v = volume(&model, &loft);
    assert!((v - swept).abs() < swept * 1e-3, "{v} against {swept}");
}

/// A square or a rectangle turned by `turn` about Z and lifted to `z`.
fn turned(model: &mut Model, z: f64, turn: f64, half: (f64, f64)) -> Shape {
    let (c, s) = (turn.cos(), turn.sin());
    let corners = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)].map(|(x, y)| {
        Point::new(
            c * x * half.0 - s * y * half.1,
            s * x * half.0 + c * y * half.1,
            z,
        )
    });
    make_polygon(model, &corners, true, T).unwrap().shape
}

/// Twisted squares and rectangles closing to a point on their axis, off
/// it, and far off it. A skin carried on through the point bends its rows
/// to meet it: past the last rectangle they shrink through a figure eight
/// and come out wound the other way, the wall near the point folded through
/// itself and facing into the material. So the last rectangle closes on
/// the point by the straight cone over it, a face of its own, and every
/// solid is the same skin through the sections plus that cone: the
/// rectangle's area times the point's height over it, a third of it,
/// wherever the point stands over the plane.
#[test]
fn a_twisted_skinned_loft_closes_outward_to_any_point() {
    let sections = |model: &mut Model| {
        vec![
            turned(model, 0.0, 0.0, (3.0, 3.0)),
            turned(model, 3.0, 0.7, (6.0, 1.0)),
            turned(model, 6.0, 1.5, (3.0, 3.0)),
            turned(model, 9.0, 2.3, (6.0, 1.0)),
        ]
    };
    let mut skin: Option<f64> = None;
    for apex_at in [
        Point::new(1.0, 2.0, 14.0),
        Point::new(0.0, 0.0, 14.0),
        Point::new(-2.0, 1.0, 12.0),
        Point::new(8.0, -6.0, 13.0),
    ] {
        let mut model = Model::new();
        let mut through = sections(&mut model);
        through.push(model.add_vertex(VertexData::new(apex_at)));
        let loft = ogeom::offset::make_loft_skinned(&mut model, &through, 1e-3, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom::algo::check(&model, &loft, T).unwrap();
        assert!(diagnosis.is_valid(), "{apex_at:?}: {diagnosis:?}");
        assert_eq!(walks::edges_walked_one_way(&model, &loft), 0, "{apex_at:?}");
        assert!(has_vertex_at(&model, &loft, apex_at));
        let faces = explore_unique(&model, &loft, ShapeType::Face).unwrap();
        assert_eq!(
            faces.len(),
            3,
            "{apex_at:?}: the wall, the cone and the cap"
        );
        let cone = 24.0 * (apex_at.z - 9.0) / 3.0;
        let v = volume(&model, &loft);
        let skin = *skin.get_or_insert(v - cone);
        assert!(skin > 0.0);
        assert!(
            (v - (skin + cone)).abs() < v * 1e-4,
            "{apex_at:?}: {v} against {}",
            skin + cone
        );
        // Short of the point, on the way to it from the last
        // rectangle's middle, is material; past it is not.
        let coarse = Deflection::default();
        let side =
            |at: Point| ogeom::algo::classify_in_solid(&model, &loft, at, coarse, T).unwrap();
        let towards = apex_at - Point::new(0.0, 0.0, 9.0);
        let near = Point::new(0.0, 0.0, 9.0) + towards * 0.8;
        let past = Point::new(0.0, 0.0, 9.0) + towards * 1.2;
        assert_eq!(side(near), ogeom::algo::Containment::In, "{apex_at:?}");
        assert_eq!(side(past), ogeom::algo::Containment::Out, "{apex_at:?}");
    }
}

/// A 2 x 2 square face centred on the origin in the XY plane.
fn square_face(model: &mut Model) -> Shape {
    let wire = rectangle(model, (-1.0, -1.0), (1.0, 1.0), 0.0);
    let plane = ogeom::math::Plane::new(Frame::WORLD);
    ogeom::algo::make_face(
        model,
        ogeom::geom::PlaneSurface::new(plane).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape
}

/// A pipe through a square and the point where its straight spine ends:
/// the pyramid, from either end.
#[test]
fn a_pipe_closes_to_a_point_at_either_end() {
    let want = 4.0 * 10.0 / 3.0;
    for point_first in [false, true] {
        let mut model = Model::new();
        let square = square_face(&mut model);
        let apex = model.add_vertex(VertexData::new(Point::new(0.0, 0.0, 10.0)));
        let spine = make_polygon(
            &mut model,
            &[Point::ORIGIN, Point::new(0.0, 0.0, 10.0)],
            false,
            T,
        )
        .unwrap()
        .shape;
        // The point first: the spine runs from it to the square.
        let (sections, spine) = if point_first {
            let reversed = make_polygon(
                &mut model,
                &[Point::new(0.0, 0.0, 10.0), Point::ORIGIN],
                false,
                T,
            )
            .unwrap()
            .shape;
            (vec![apex, square], reversed)
        } else {
            (vec![square, apex], spine)
        };
        let pipe = ogeom::offset::make_pipe_sections(&mut model, &sections, &spine, false, 1e-3, T)
            .unwrap()
            .shape;
        assert!(ogeom::algo::check(&model, &pipe, T).unwrap().is_valid());
        assert!(has_vertex_at(&model, &pipe, Point::new(0.0, 0.0, 10.0)));
        let v = volume(&model, &pipe);
        assert!((v - want).abs() < want * 1e-3, "{v} against {want}");
    }
}

/// Along a quarter circle the pipe closes exactly where the spine ends. Its
/// squares shrink about the spine, so it holds each square's area along the
/// arc: 4 L / 3 by Pappus.
#[test]
fn a_bent_pipe_closes_to_the_spine_s_end() {
    let mut model = Model::new();
    let square = square_face(&mut model);
    // About (10, 0, 0) in the XZ plane, from the origin heading up Z round
    // to (10, 0, 10).
    let frame = Frame::new(Point::new(10.0, 0.0, 0.0), Direction::Y, -Direction::X, T).unwrap();
    let curve = CircleCurve::new(Circle::new(frame, 10.0, T).unwrap());
    let edge = make_edge(&mut model, curve.into(), (0.0, PI / 2.0), T)
        .unwrap()
        .shape;
    let spine = make_wire(&mut model, &[edge], T).unwrap().shape;
    let end = Point::new(10.0, 0.0, 10.0);
    let apex = model.add_vertex(VertexData::new(end));
    let pipe =
        ogeom::offset::make_pipe_sections(&mut model, &[square, apex], &spine, false, 1e-3, T)
            .unwrap()
            .shape;
    assert!(ogeom::algo::check(&model, &pipe, T).unwrap().is_valid());
    assert!(has_vertex_at(&model, &pipe, end));
    let want = 4.0 * (10.0 * PI / 2.0) / 3.0;
    let v = volume(&model, &pipe);
    assert!((v - want).abs() < want * 5e-3, "{v} against {want}");
}

/// Squares shrinking without a turn close smoothly to a point over them,
/// on their axis or off it: one wall through the sections and on to the
/// point, its normal turning gently across the last section (no crease
/// there), and the solid valid.
#[test]
fn a_gentle_skinned_loft_stays_smooth_to_the_point() {
    use ogeom::geom::Surface as _;
    for apex_at in [Point::new(0.0, 0.0, 10.0), Point::new(1.0, 1.0, 10.0)] {
        let mut model = Model::new();
        let sections = vec![
            turned(&mut model, 0.0, 0.0, (4.0, 4.0)),
            turned(&mut model, 3.0, 0.0, (3.0, 3.0)),
            turned(&mut model, 6.0, 0.0, (2.0, 2.0)),
            model.add_vertex(VertexData::new(apex_at)),
        ];
        let loft = ogeom::offset::make_loft_skinned(&mut model, &sections, 1e-3, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom::algo::check(&model, &loft, T).unwrap();
        assert!(diagnosis.is_valid(), "{apex_at:?}: {diagnosis:?}");
        assert_eq!(walks::edges_walked_one_way(&model, &loft), 0, "{apex_at:?}");
        let faces = explore_unique(&model, &loft, ShapeType::Face).unwrap();
        assert_eq!(faces.len(), 2, "{apex_at:?}: the wall and the cap");
        assert!(volume(&model, &loft) > 0.0);
        // The wall is the face reaching the point; across its chart, short
        // of the point, neighbouring normals up the wall stay close, the
        // last section's row among them.
        let wall = faces
            .iter()
            .map(|f| model.node(f).unwrap().data().as_face().unwrap().surface)
            .map(|id| model.geometry().surface(id).unwrap().clone())
            .find(|s| {
                let ((u0, _), (_, v1)) = s.domain();
                s.point_at(u0, v1, T).unwrap().distance(apex_at) < 1e-6
            })
            .unwrap();
        let ((u0, u1), (v0, v1)) = wall.domain();
        for i in 0..64 {
            let u = u0 + (u1 - u0) * f64::from(i) / 64.0;
            for j in 0..200 {
                let v = v0 + (v1 - v0) * 0.9 * f64::from(j) / 200.0;
                let dv = (v1 - v0) * 0.9 / 200.0;
                let a = wall.normal_at(u, v, T).unwrap().vector();
                let b = wall.normal_at(u, v + dv, T).unwrap().vector();
                assert!(a.dot(b) > 0.9, "{apex_at:?}: a crease at ({u}, {v})");
            }
        }
    }
}
