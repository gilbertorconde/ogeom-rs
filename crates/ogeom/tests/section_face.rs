//! Sectioning a face or an open shell by a plane: the curves where its
//! faces cross the plane, trimmed to the faces, with no solid behind them.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::{Curve3d as _, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Direction, Frame, Plane, Point};
use ogeom::topo::{EdgeRepr, Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn edges(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Edge).unwrap()
}

/// An edge's two ends.
fn ends(model: &Model, edge: &Shape) -> (Point, Point) {
    let data = model.node(edge).unwrap().data().as_edge().unwrap().clone();
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        panic!("an edge with no curve")
    };
    let c = model.geometry().curve(*curve).unwrap();
    (
        c.point_at(range.0, T).unwrap(),
        c.point_at(range.1, T).unwrap(),
    )
}

/// The cube's face whose points all have the given x.
fn face_at_x(model: &Model, solid: &Shape, x: f64) -> Shape {
    explore_unique(model, solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let b = ogeom::algo::shape_bounds(model, f, T).unwrap();
            (b.low().unwrap().x - x).abs() < 1e-6 && (b.high().unwrap().x - x).abs() < 1e-6
        })
        .unwrap()
}

fn half_space_below(model: &mut Model, height: f64) -> Shape {
    let plane = Plane::through(Point::new(0.0, 0.0, height), Direction::Z);
    let surface = PlaneSurface::over(plane, (-1e6, 1e6), (-1e6, 1e6)).unwrap();
    let face = ogeom::algo::make_natural_face(model, SurfaceGeometry::Plane(surface))
        .unwrap()
        .shape;
    ogeom::algo::make_half_space(model, &face, Point::new(0.0, 0.0, height - 1.0), T)
        .unwrap()
        .shape
}

#[test]
fn a_face_on_its_own_is_sectioned_by_a_half_space() {
    let mut model = Model::new();
    let cube = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let side = face_at_x(&model, &cube, 0.0);
    let half = half_space_below(&mut model, 4.0);
    let section = ogeom::boolean::section(&mut model, &side, &half, T)
        .unwrap()
        .shape;
    let found = edges(&model, &section);
    assert_eq!(found.len(), 1);
    let (a, b) = ends(&model, &found[0]);
    let (a, b) = if a.y <= b.y { (a, b) } else { (b, a) };
    assert!(a.distance(Point::new(0.0, 0.0, 4.0)) < 1e-9, "{a:?}");
    assert!(b.distance(Point::new(0.0, 10.0, 4.0)) < 1e-9, "{b:?}");
}

#[test]
fn a_plane_that_misses_a_face_sections_nothing_and_one_holding_it_its_boundary() {
    let mut model = Model::new();
    let cube = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let side = face_at_x(&model, &cube, 0.0);
    let clear = Plane::through(Point::new(0.0, 0.0, 12.0), Direction::Z);
    let none = ogeom::boolean::section_face(&mut model, &side, &clear, T)
        .unwrap()
        .shape;
    assert!(edges(&model, &none).is_empty());
    let holding = Plane::through(Point::ORIGIN, Direction::X);
    let boundary = ogeom::boolean::section_face(&mut model, &side, &holding, T)
        .unwrap()
        .shape;
    assert_eq!(edges(&model, &boundary).len(), 4);
}

#[test]
fn a_drum_s_side_sections_to_its_circle() {
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 3.0, 10.0, T)
        .unwrap()
        .shape;
    let wall = explore_unique(&model, &drum, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            matches!(
                model.geometry().surface(data.surface),
                Some(SurfaceGeometry::Cylinder(_))
            )
        })
        .unwrap();
    let plane = Plane::through(Point::new(0.0, 0.0, 6.0), Direction::Z);
    let section = ogeom::boolean::section_face(&mut model, &wall, &plane, T)
        .unwrap()
        .shape;
    let found = edges(&model, &section);
    assert!(!found.is_empty());
    let length: f64 = found
        .iter()
        .map(|e| {
            let data = model.node(e).unwrap().data().as_edge().unwrap().clone();
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                panic!()
            };
            let c = model.geometry().curve(*curve).unwrap();
            (0..400)
                .map(|k| {
                    let t = |k: i32| range.0 + (range.1 - range.0) * f64::from(k) / 400.0;
                    c.point_at(t(k), T)
                        .unwrap()
                        .distance(c.point_at(t(k + 1), T).unwrap())
                })
                .sum::<f64>()
        })
        .sum();
    let want = core::f64::consts::TAU * 3.0;
    assert!(
        (length - want).abs() < want * 1e-4,
        "{length} against {want}"
    );
}

#[test]
fn an_open_shell_s_sections_meet_at_the_edge_its_faces_share() {
    let mut model = Model::new();
    let cube = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let two = [
        face_at_x(&model, &cube, 0.0),
        face_at_x(&model, &cube, 10.0),
    ];
    let front = explore_unique(&model, &cube, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let b = ogeom::algo::shape_bounds(&model, f, T).unwrap();
            b.high().unwrap().y.abs() < 1e-6
        })
        .unwrap();
    let shell = ogeom::algo::make_shell(&mut model, &[two[0].clone(), front, two[1].clone()])
        .unwrap()
        .shape;
    let plane = Plane::through(Point::new(0.0, 0.0, 4.0), Direction::Z);
    let section = ogeom::boolean::section_face(&mut model, &shell, &plane, T)
        .unwrap()
        .shape;
    assert_eq!(edges(&model, &section).len(), 3);
    // Three edges end to end: four distinct vertices.
    assert_eq!(
        explore_unique(&model, &section, ShapeType::Vertex)
            .unwrap()
            .len(),
        4
    );
}

/// The face of `solid` whose surface is not a plane.
fn curved_face(model: &Model, solid: &Shape) -> Shape {
    explore_unique(model, solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            !matches!(
                model.geometry().surface(data.surface),
                Some(SurfaceGeometry::Plane(_))
            )
        })
        .unwrap()
}

/// Points along each edge of `section`, `n` per edge, ends included.
fn samples(model: &Model, section: &Shape, n: u32) -> Vec<Vec<Point>> {
    edges(model, section)
        .iter()
        .map(|e| {
            let data = model.node(e).unwrap().data().as_edge().unwrap().clone();
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                panic!("an edge with no curve")
            };
            let c = model.geometry().curve(*curve).unwrap();
            (0..=n)
                .map(|k| {
                    let t = range.0 + (range.1 - range.0) * f64::from(k) / f64::from(n);
                    c.point_at(t, T).unwrap()
                })
                .collect()
        })
        .collect()
}

fn polyline_length(points: &[Point]) -> f64 {
    points.windows(2).map(|w| w[0].distance(w[1])).sum()
}

#[test]
fn a_cylinder_cut_through_its_seam_keeps_the_line_on_the_seam() {
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 12.0, T)
        .unwrap()
        .shape;
    let wall = curved_face(&model, &drum);
    let plane = Plane::through(Point::ORIGIN, Direction::Y);
    let section = ogeom::boolean::section_face(&mut model, &wall, &plane, T)
        .unwrap()
        .shape;
    let lines = samples(&model, &section, 16);
    assert_eq!(
        lines.len(),
        2,
        "one line at x = -5, one on the seam at x = +5"
    );
    let mut xs: Vec<f64> = lines.iter().map(|l| l[0].x).collect();
    xs.sort_by(f64::total_cmp);
    let mut worst = 0.0_f64;
    for (line, x) in lines.iter().zip([xs[0], xs[1]]) {
        // On the plane y = 0 and the cylinder x^2 + y^2 = 25, from z = 0
        // to z = 12.
        for p in line {
            worst = worst.max(p.y.abs()).max((p.x.hypot(p.y) - 5.0).abs());
        }
        assert!((polyline_length(line) - 12.0).abs() < 1e-9, "{x}");
    }
    assert!(
        (xs[0] + 5.0).abs() < 1e-9 && (xs[1] - 5.0).abs() < 1e-9,
        "{xs:?}"
    );
    assert!(worst < 1e-9, "{worst}");
}

#[test]
fn a_torus_cut_through_its_seam_keeps_both_equators() {
    let mut model = Model::new();
    let ring = ogeom::algo::make_torus(&mut model, Frame::WORLD, 8.0, 2.0, T)
        .unwrap()
        .shape;
    let skin = curved_face(&model, &ring);
    let plane = Plane::through(Point::ORIGIN, Direction::Z);
    let section = ogeom::boolean::section_face(&mut model, &skin, &plane, T)
        .unwrap()
        .shape;
    let circles = samples(&model, &section, 720);
    let mut radii: Vec<f64> = Vec::new();
    let mut worst = 0.0_f64;
    for circle in &circles {
        let r = circle[0].x.hypot(circle[0].y);
        // On z = 0 and the torus (sqrt(x^2 + y^2) - 8)^2 + z^2 = 4, at one
        // radius all the way round.
        for p in circle {
            let rho = p.x.hypot(p.y);
            worst = worst
                .max(p.z.abs())
                .max(((rho - 8.0).hypot(p.z) - 2.0).abs())
                .max((rho - r).abs());
        }
        let length = polyline_length(circle);
        let want = core::f64::consts::TAU * r;
        assert!(
            (length - want).abs() < want * 1e-4,
            "{length} against {want}"
        );
        radii.push(r);
    }
    radii.sort_by(f64::total_cmp);
    assert_eq!(radii.len(), 2, "{radii:?}");
    assert!(
        (radii[0] - 6.0).abs() < 1e-9 && (radii[1] - 10.0).abs() < 1e-9,
        "{radii:?}"
    );
    assert!(worst < 1e-9, "{worst}");
}

#[test]
fn a_sphere_cut_through_its_seam_keeps_the_seam_half_of_the_circle() {
    let mut model = Model::new();
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 4.0, T)
        .unwrap()
        .shape;
    let skin = curved_face(&model, &ball);
    let plane = Plane::through(Point::ORIGIN, Direction::Y);
    let section = ogeom::boolean::section_face(&mut model, &skin, &plane, T)
        .unwrap()
        .shape;
    let arcs = samples(&model, &section, 360);
    let mut worst = 0.0_f64;
    for arc in &arcs {
        for p in arc {
            worst = worst
                .max(p.y.abs())
                .max((p.distance(Point::ORIGIN) - 4.0).abs());
        }
    }
    let length: f64 = arcs.iter().map(|a| polyline_length(a)).sum();
    let want = core::f64::consts::TAU * 4.0;
    assert!(
        (length - want).abs() < want * 1e-4,
        "{length} against {want}"
    );
    assert!(worst < 1e-9, "{worst}");
}

#[test]
fn a_box_cut_through_one_of_its_edges_keeps_that_edge_once() {
    let mut model = Model::new();
    let cube = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    // x + y = 0 holds the edge x = y = 0 and leaves the rest of the box on
    // one side.
    let normal = Direction::new(ogeom::math::Vector::new(1.0, 1.0, 0.0), T).unwrap();
    let plane = Plane::through(Point::ORIGIN, normal);
    for shape in [cube.clone(), face_at_x(&model, &cube, 0.0)] {
        let section = ogeom::boolean::section_face(&mut model, &shape, &plane, T)
            .unwrap()
            .shape;
        let lines = samples(&model, &section, 4);
        assert_eq!(lines.len(), 1);
        let line = &lines[0];
        assert!(line.iter().all(|p| p.x.abs() < 1e-9 && p.y.abs() < 1e-9));
        assert!((polyline_length(line) - 10.0).abs() < 1e-9);
    }
}
