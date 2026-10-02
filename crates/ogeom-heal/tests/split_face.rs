//! Faces split along curves on them or projected onto them, measured
//! against closed-form areas and against quadrature of the surface.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{PI, TAU};

use ogeom_algo::{make_box, make_edge, make_face, make_polygon, surface_properties};
use ogeom_core::Tolerances;
use ogeom_geom::{
    BSplineSurface, CircleCurve, Curve, LineCurve, PlaneSurface, Surface as _, SurfaceGeometry,
};
use ogeom_geom::{Curve2d as _, Curve3d as _};
use ogeom_heal::{Projection, split_face};
use ogeom_math::{Axis, Circle, ControlGrid, Direction, Frame, KnotVector, Plane, Point, Vector};
use ogeom_mesh::Deflection;
use ogeom_topo::{EdgeRepr, Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn fine() -> Deflection {
    Deflection {
        chord: 1e-4,
        angular: 0.01,
        ..Deflection::default()
    }
}

fn faces(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Face).unwrap()
}

fn area(model: &Model, face: &Shape) -> f64 {
    surface_properties(model, face, fine(), T).unwrap().mass
}

fn edge_nodes(model: &Model, face: &Shape) -> Vec<ogeom_topo::TShapeId> {
    explore_unique(model, face, ShapeType::Edge)
        .unwrap()
        .iter()
        .map(Shape::node)
        .collect()
}

fn wires(model: &Model, face: &Shape) -> usize {
    explore_unique(model, face, ShapeType::Wire).unwrap().len()
}

/// The square `[0, 10] x [0, 10]` in the plane `z = 0`.
fn square(model: &mut Model) -> Shape {
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(10.0, 0.0, 0.0),
        Point::new(10.0, 10.0, 0.0),
        Point::new(0.0, 10.0, 0.0),
    ];
    let wire = make_polygon(model, &corners, true, T).unwrap().shape;
    let plane: SurfaceGeometry = PlaneSurface::new(Plane::new(Frame::WORLD)).into();
    make_face(model, plane, &[wire], T).unwrap().shape
}

fn segment(model: &mut Model, from: Point, to: Point) -> Shape {
    let along = to - from;
    let line = LineCurve::new(Axis {
        location: from,
        direction: Direction::new(along, T).unwrap(),
    });
    make_edge(model, Curve::from(line), (0.0, along.magnitude()), T)
        .unwrap()
        .shape
}

/// A full circle in a plane parallel to `z = 0`, its parameter starting on
/// the `x` side the frame names.
fn circle(model: &mut Model, centre: Point, x: Vector, radius: f64) -> Shape {
    let frame = Frame::new(
        centre,
        Direction::new(Vector::new(0.0, 0.0, 1.0), T).unwrap(),
        Direction::new(x, T).unwrap(),
        T,
    )
    .unwrap();
    let curve = CircleCurve::new(Circle::new(frame, radius, T).unwrap());
    make_edge(model, Curve::from(curve), (0.0, TAU), T)
        .unwrap()
        .shape
}

fn refusal(result: ogeom_core::OgeomResult<ogeom_algo::Built>) -> String {
    match result {
        Ok(_) => panic!("the split should have been refused"),
        Err(e) => e.to_string(),
    }
}

#[test]
fn a_square_split_across_its_middle_gives_two_halves_sharing_one_edge() {
    let mut model = Model::new();
    let face = square(&mut model);
    let line = segment(
        &mut model,
        Point::new(0.0, 5.0, 0.0),
        Point::new(10.0, 5.0, 0.0),
    );
    let built = split_face(
        &mut model,
        &face,
        &face,
        std::slice::from_ref(&line),
        Projection::OnFace,
        T,
    )
    .unwrap();
    let pieces = faces(&model, &built.shape);
    assert_eq!(pieces.len(), 2);
    for piece in &pieces {
        assert!(
            (area(&model, piece) - 50.0).abs() <= 1e-9,
            "{}",
            area(&model, piece)
        );
    }
    let (a, b) = (
        edge_nodes(&model, &pieces[0]),
        edge_nodes(&model, &pieces[1]),
    );
    let shared: Vec<_> = a.iter().filter(|e| b.contains(e)).collect();
    assert_eq!(shared.len(), 1, "one edge between the halves");
    // The history names the halves and the edge cut along the line.
    let mut images: Vec<_> = built
        .history
        .modified(&face)
        .iter()
        .map(Shape::node)
        .collect();
    images.sort();
    let mut found: Vec<_> = pieces.iter().map(Shape::node).collect();
    found.sort();
    assert_eq!(images, found);
    assert_eq!(built.history.generated(&line).len(), 1);
    assert_eq!(built.history.generated(&line)[0].node(), *shared[0]);
}

#[test]
fn a_circle_inside_the_square_cuts_out_a_disk_and_leaves_a_hole() {
    let mut model = Model::new();
    let face = square(&mut model);
    let ring = circle(
        &mut model,
        Point::new(5.0, 5.0, 0.0),
        Vector::new(1.0, 0.0, 0.0),
        2.0,
    );
    let built = split_face(&mut model, &face, &face, &[ring], Projection::OnFace, T).unwrap();
    let pieces = faces(&model, &built.shape);
    assert_eq!(pieces.len(), 2);
    let disk = pieces.iter().find(|f| wires(&model, f) == 1).unwrap();
    let holed = pieces.iter().find(|f| wires(&model, f) == 2).unwrap();
    let band = 1e-5;
    assert!((area(&model, disk) - 4.0 * PI).abs() <= 4.0 * PI * band);
    assert!((area(&model, holed) - (100.0 - 4.0 * PI)).abs() <= 100.0 * band);
    let shared: Vec<_> = edge_nodes(&model, disk)
        .into_iter()
        .filter(|e| edge_nodes(&model, holed).contains(e))
        .collect();
    assert_eq!(shared.len(), 1, "the circle bounds both");
}

#[test]
fn a_closed_curve_starting_outside_is_cut_where_it_leaves_the_face() {
    // A circle about the middle of the left side, starting outside: the
    // half inside cuts a half disk off the square.
    let mut model = Model::new();
    let face = square(&mut model);
    let ring = circle(
        &mut model,
        Point::new(0.0, 5.0, 0.0),
        Vector::new(-1.0, 0.0, 0.0),
        2.0,
    );
    let built = split_face(&mut model, &face, &face, &[ring], Projection::OnFace, T).unwrap();
    let mut areas: Vec<f64> = faces(&model, &built.shape)
        .iter()
        .map(|f| area(&model, f))
        .collect();
    areas.sort_by(f64::total_cmp);
    assert_eq!(areas.len(), 2);
    assert!((areas[0] - 2.0 * PI).abs() <= 2.0 * PI * 1e-5, "{areas:?}");
    assert!(
        (areas[1] - (100.0 - 2.0 * PI)).abs() <= 100.0 * 1e-5,
        "{areas:?}"
    );
}

#[test]
fn two_crossing_lines_quarter_the_square() {
    let mut model = Model::new();
    let face = square(&mut model);
    let across = segment(
        &mut model,
        Point::new(0.0, 5.0, 0.0),
        Point::new(10.0, 5.0, 0.0),
    );
    let up = segment(
        &mut model,
        Point::new(5.0, 0.0, 0.0),
        Point::new(5.0, 10.0, 0.0),
    );
    let built = split_face(
        &mut model,
        &face,
        &face,
        &[across.clone(), up.clone()],
        Projection::OnFace,
        T,
    )
    .unwrap();
    let pieces = faces(&model, &built.shape);
    assert_eq!(pieces.len(), 4);
    for piece in &pieces {
        assert!((area(&model, piece) - 25.0).abs() <= 1e-9);
    }
    assert_eq!(built.history.modified(&face).len(), 4);
    // Each line ends as two edges, halved where the other crosses it, and
    // every one of them bounds two quarters.
    for line in [&across, &up] {
        let made = built.history.generated(line);
        assert_eq!(made.len(), 2, "{line:?}");
        for edge in made {
            let holders = pieces
                .iter()
                .filter(|p| edge_nodes(&model, p).contains(&edge.node()))
                .count();
            assert_eq!(holders, 2);
        }
    }
}

#[test]
fn a_line_ending_on_an_earlier_cut_divides_a_half() {
    // The second line starts on the first: a T, three pieces.
    let mut model = Model::new();
    let face = square(&mut model);
    let across = segment(
        &mut model,
        Point::new(0.0, 5.0, 0.0),
        Point::new(10.0, 5.0, 0.0),
    );
    let up = segment(
        &mut model,
        Point::new(4.0, 5.0, 0.0),
        Point::new(4.0, 10.0, 0.0),
    );
    // Given in the order that makes the second wait for the first.
    let built = split_face(
        &mut model,
        &face,
        &face,
        &[up, across],
        Projection::OnFace,
        T,
    )
    .unwrap();
    let mut areas: Vec<f64> = faces(&model, &built.shape)
        .iter()
        .map(|f| area(&model, f))
        .collect();
    areas.sort_by(f64::total_cmp);
    let expected = [20.0, 30.0, 50.0];
    assert_eq!(areas.len(), 3);
    for (got, want) in areas.iter().zip(expected) {
        assert!((got - want).abs() <= 1e-9, "{areas:?}");
    }
}

#[test]
fn a_box_face_split_across_keeps_the_box_closed() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let top = faces(&model, &block)
        .into_iter()
        .find(|f| {
            let Some(NodeData::Face(data)) = model.node(f).map(|n| n.data()) else {
                return false;
            };
            let s = model.geometry().surface(data.surface).unwrap();
            s.point_at(0.0, 0.0, T).unwrap().z > 9.0 && s.point_at(1.0, 1.0, T).unwrap().z > 9.0
        })
        .unwrap();
    let line = segment(
        &mut model,
        Point::new(-1.0, 3.0, 10.0),
        Point::new(11.0, 3.0, 10.0),
    );
    let built = split_face(&mut model, &block, &top, &[line], Projection::OnFace, T).unwrap();
    assert_eq!(faces(&model, &built.shape).len(), 7);
    let diagnosis = ogeom_algo::check(&model, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let volume = ogeom_algo::volume_properties(&model, &built.shape, fine(), T)
        .unwrap()
        .mass;
    assert!((volume - 1000.0).abs() <= 1e-6, "{volume}");
    let mut areas: Vec<f64> = built
        .history
        .modified(&top)
        .iter()
        .map(|f| area(&model, f))
        .collect();
    areas.sort_by(f64::total_cmp);
    assert_eq!(areas.len(), 2);
    assert!((areas[0] - 30.0).abs() <= 1e-9 && (areas[1] - 70.0).abs() <= 1e-9);
}

/// A bicubic patch over `[0, 10] x [0, 10]` whose `x` and `y` run with its
/// parameters, `x = 10 u` and `y = 10 v`, rising and falling in `z`.
fn bump() -> BSplineSurface {
    let heights = [
        [0.0, 1.0, 0.5, 0.0],
        [0.5, 3.0, 2.0, 1.0],
        [1.0, 2.0, 4.0, 0.5],
        [0.0, 0.5, 1.0, 0.0],
    ];
    let mut points = Vec::new();
    for (i, row) in heights.iter().enumerate() {
        for (j, z) in row.iter().enumerate() {
            points.push(Point::new(
                10.0 * f64::from(u8::try_from(i).unwrap()) / 3.0,
                10.0 * f64::from(u8::try_from(j).unwrap()) / 3.0,
                *z,
            ));
        }
    }
    let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3).unwrap();
    let grid = ControlGrid::new(points, 4, 4).unwrap();
    BSplineSurface::new(knots.clone(), knots, &grid, T).unwrap()
}

/// The area of the patch over `v` in `[v0, v1]`, by Gauss-Legendre
/// quadrature of `|Su x Sv|` on a grid of cells.
fn quadrature(surface: &SurfaceGeometry, v0: f64, v1: f64) -> f64 {
    let nodes = [
        (-0.906_179_845_938_664, 0.236_926_885_056_189),
        (-0.538_469_310_105_683, 0.478_628_670_499_366),
        (0.0, 0.568_888_888_888_889),
        (0.538_469_310_105_683, 0.478_628_670_499_366),
        (0.906_179_845_938_664, 0.236_926_885_056_189),
    ];
    let cells = 16;
    let mut total = 0.0;
    for a in 0..cells {
        for b in 0..cells {
            let (ua, ub) = (f64::from(a) / 16.0, f64::from(a + 1) / 16.0);
            let (va, vb) = (
                v0 + (v1 - v0) * f64::from(b) / 16.0,
                v0 + (v1 - v0) * f64::from(b + 1) / 16.0,
            );
            for (x, wx) in nodes {
                for (y, wy) in nodes {
                    let u = 0.5 * (ub - ua).mul_add(x, ua + ub);
                    let v = 0.5 * (vb - va).mul_add(y, va + vb);
                    let (su, sv) = surface.d1_at(u, v, T).unwrap();
                    total += wx * wy * su.cross(sv).magnitude() * 0.25 * (ub - ua) * (vb - va);
                }
            }
        }
    }
    total
}

#[test]
fn a_spline_face_split_by_a_line_dropped_along_z_halves_its_area() {
    let mut model = Model::new();
    let surface: SurfaceGeometry = bump().into();
    let face = ogeom_algo::make_natural_face(&mut model, surface.clone())
        .unwrap()
        .shape;
    let whole = area(&model, &face);
    let full = quadrature(&surface, 0.0, 1.0);
    assert!(
        (whole - full).abs() <= full * 1e-5,
        "{whole} against {full}"
    );
    // Overshooting both sides; dropped straight down it lands on y = 5,
    // the patch's line v = 1/2.
    let line = segment(
        &mut model,
        Point::new(-2.0, 5.0, 20.0),
        Point::new(12.0, 5.0, 20.0),
    );
    let down = Direction::new(Vector::new(0.0, 0.0, -1.0), T).unwrap();
    let built = split_face(
        &mut model,
        &face,
        &face,
        &[line],
        Projection::Along(down),
        T,
    )
    .unwrap();
    let pieces = faces(&model, &built.shape);
    assert_eq!(pieces.len(), 2);
    let areas: Vec<f64> = pieces.iter().map(|f| area(&model, f)).collect();
    assert!(
        (areas[0] + areas[1] - whole).abs() <= whole * 1e-5,
        "{areas:?} against {whole}"
    );
    let low = quadrature(&surface, 0.0, 0.5);
    let high = quadrature(&surface, 0.5, 1.0);
    let (small, large) = (areas[0].min(areas[1]), areas[0].max(areas[1]));
    assert!(
        (small - low.min(high)).abs() <= low * 1e-5,
        "{small} against {low}, {high}"
    );
    assert!(
        (large - low.max(high)).abs() <= high * 1e-5,
        "{large} against {low}, {high}"
    );
}

#[test]
fn a_spline_face_split_by_a_line_dropped_along_its_normals_keeps_its_area() {
    let mut model = Model::new();
    let surface: SurfaceGeometry = bump().into();
    let face = ogeom_algo::make_natural_face(&mut model, surface.clone())
        .unwrap()
        .shape;
    let whole = area(&model, &face);
    // Far past the patch both ways, so the ends land on its edges.
    let line = segment(
        &mut model,
        Point::new(-30.0, 4.0, 5.0),
        Point::new(40.0, 6.0, 5.0),
    );
    let built = split_face(
        &mut model,
        &face,
        &face,
        std::slice::from_ref(&line),
        Projection::AlongNormals,
        T,
    )
    .unwrap();
    let pieces = faces(&model, &built.shape);
    assert_eq!(pieces.len(), 2);
    let sum: f64 = pieces.iter().map(|f| area(&model, f)).sum();
    assert!((sum - whole).abs() <= whole * 1e-5, "{sum} against {whole}");
    // The edge cut along the projection: its curve and the surface along
    // its trace agree within the edge's tolerance, and each point of it is
    // where the normal through the line meets the surface.
    let cut = built.history.generated(&line)[0].clone();
    let data = model.node(&cut).unwrap().data().as_edge().unwrap().clone();
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d().cloned() else {
        panic!("the cut edge has a curve");
    };
    let curve = model.geometry().curve(curve).unwrap().clone();
    let Some(EdgeRepr::PCurve { curve: trace, .. }) = data
        .representations
        .iter()
        .find(|r| matches!(r, EdgeRepr::PCurve { .. }))
        .cloned()
    else {
        panic!("the cut edge has a trace");
    };
    let trace = model.geometry().pcurve(trace).unwrap().clone();
    let (a, b) = (Point::new(-30.0, 4.0, 5.0), Point::new(40.0, 6.0, 5.0));
    let along = (b - a) * (1.0 / (b - a).magnitude());
    for k in 0..=50 {
        let t = (range.1 - range.0).mul_add(f64::from(k) / 50.0, range.0);
        let p = curve.point_at(t, T).unwrap();
        let uv = trace.point_at(t, T).unwrap();
        let on = surface.point_at(uv.x, uv.y, T).unwrap();
        assert!(p.distance(on) <= data.tolerance.get(), "{}", p.distance(on));
        // The normal at `p` meets the line: their distance apart is zero.
        let (su, sv) = surface.d1_at(uv.x, uv.y, T).unwrap();
        let across = su.cross(sv).cross(along);
        let skew = across.dot(p - a).abs() / across.magnitude();
        assert!(
            skew <= 1e-5,
            "the normal at {p:?} misses the line by {skew}"
        );
    }
}

#[test]
fn a_line_ending_inside_the_face_is_refused() {
    let mut model = Model::new();
    let face = square(&mut model);
    let line = segment(
        &mut model,
        Point::new(0.0, 5.0, 0.0),
        Point::new(5.0, 5.0, 0.0),
    );
    let why = refusal(split_face(
        &mut model,
        &face,
        &face,
        &[line],
        Projection::OnFace,
        T,
    ));
    assert!(why.contains("ends inside the face"), "{why}");
}

#[test]
fn a_line_off_the_face_is_refused_unless_projected() {
    let mut model = Model::new();
    let face = square(&mut model);
    let line = segment(
        &mut model,
        Point::new(0.0, 5.0, 1.0),
        Point::new(10.0, 5.0, 1.0),
    );
    let why = refusal(split_face(
        &mut model,
        &face,
        &face,
        std::slice::from_ref(&line),
        Projection::OnFace,
        T,
    ));
    assert!(why.contains("does not lie on the face"), "{why}");
    let built = split_face(
        &mut model,
        &face,
        &face,
        &[line],
        Projection::AlongNormals,
        T,
    )
    .unwrap();
    for piece in faces(&model, &built.shape) {
        assert!((area(&model, &piece) - 50.0).abs() <= 1e-6);
    }
}

#[test]
fn a_line_missing_the_face_is_refused() {
    let mut model = Model::new();
    let face = square(&mut model);
    let line = segment(
        &mut model,
        Point::new(20.0, 0.0, 0.0),
        Point::new(20.0, 10.0, 0.0),
    );
    let why = refusal(split_face(
        &mut model,
        &face,
        &face,
        &[line],
        Projection::OnFace,
        T,
    ));
    assert!(why.contains("does not cross the face"), "{why}");
}

#[test]
fn a_projection_that_never_lands_is_refused() {
    let mut model = Model::new();
    let face = square(&mut model);
    let line = segment(
        &mut model,
        Point::new(0.0, 5.0, 1.0),
        Point::new(10.0, 5.0, 1.0),
    );
    let sideways = Direction::new(Vector::new(0.0, 1.0, 0.0), T).unwrap();
    let why = refusal(split_face(
        &mut model,
        &face,
        &face,
        &[line],
        Projection::Along(sideways),
        T,
    ));
    assert!(why.contains("misses the face"), "{why}");
}

#[test]
fn a_line_along_the_boundary_is_refused() {
    let mut model = Model::new();
    let face = square(&mut model);
    let line = segment(
        &mut model,
        Point::new(-1.0, 0.0, 0.0),
        Point::new(11.0, 0.0, 0.0),
    );
    let why = refusal(split_face(
        &mut model,
        &face,
        &face,
        &[line],
        Projection::OnFace,
        T,
    ));
    assert!(why.contains("runs along the face's boundary"), "{why}");
}

#[test]
fn a_line_from_the_outside_to_a_hole_is_refused() {
    let mut model = Model::new();
    let (shell, holed) = holed_square(&mut model);
    let line = segment(
        &mut model,
        Point::new(0.0, 5.0, 0.0),
        Point::new(3.0, 5.0, 0.0),
    );
    let why = refusal(split_face(
        &mut model,
        &shell,
        &holed,
        &[line],
        Projection::OnFace,
        T,
    ));
    assert!(why.contains("from one boundary loop"), "{why}");
}

#[test]
fn a_closed_curve_starting_inside_is_walked_from_where_it_meets_the_boundary() {
    // The circle about the middle of the left side again, its start now
    // inside the square: the same half disk comes off.
    let mut model = Model::new();
    let face = square(&mut model);
    let ring = circle(
        &mut model,
        Point::new(0.0, 5.0, 0.0),
        Vector::new(1.0, 0.0, 0.0),
        2.0,
    );
    let built = split_face(&mut model, &face, &face, &[ring], Projection::OnFace, T).unwrap();
    let mut areas: Vec<f64> = faces(&model, &built.shape)
        .iter()
        .map(|f| area(&model, f))
        .collect();
    areas.sort_by(f64::total_cmp);
    assert_eq!(areas.len(), 2);
    assert!((areas[0] - 2.0 * PI).abs() <= 2.0 * PI * 1e-5, "{areas:?}");
    assert!(
        (areas[1] - (100.0 - 2.0 * PI)).abs() <= 100.0 * 1e-5,
        "{areas:?}"
    );
}

#[test]
fn a_split_without_edges_or_of_a_foreign_face_is_refused() {
    let mut model = Model::new();
    let face = square(&mut model);
    let other = square(&mut model);
    let why = refusal(split_face(
        &mut model,
        &face,
        &face,
        &[],
        Projection::OnFace,
        T,
    ));
    assert!(why.contains("at least one edge"), "{why}");
    let line = segment(
        &mut model,
        Point::new(0.0, 5.0, 0.0),
        Point::new(10.0, 5.0, 0.0),
    );
    let why = refusal(split_face(
        &mut model,
        &face,
        &other,
        &[line],
        Projection::OnFace,
        T,
    ));
    assert!(why.contains("not a face of the shape"), "{why}");
}

#[test]
fn a_polyline_wire_is_followed_corner_to_corner() {
    // Two segments joined inside the square, given as one open wire.
    let mut model = Model::new();
    let face = square(&mut model);
    let corners = [
        Point::new(0.0, 2.0, 0.0),
        Point::new(4.0, 8.0, 0.0),
        Point::new(10.0, 3.0, 0.0),
    ];
    let wire = make_polygon(&mut model, &corners, false, T).unwrap().shape;
    let built = split_face(&mut model, &face, &face, &[wire], Projection::OnFace, T).unwrap();
    let mut areas: Vec<f64> = faces(&model, &built.shape)
        .iter()
        .map(|f| area(&model, f))
        .collect();
    areas.sort_by(f64::total_cmp);
    // Below the polyline: the pentagon (0,0) (10,0) (10,3) (4,8) (0,2),
    // 53 by the shoelace formula.
    assert_eq!(areas.len(), 2);
    assert!((areas[0] - 47.0).abs() <= 1e-9, "{areas:?}");
    assert!((areas[1] - 53.0).abs() <= 1e-9, "{areas:?}");
}

/// The drum of radius 2 and height 5 on the world frame, and its wall.
fn drum(model: &mut Model) -> (Shape, Shape) {
    let drum = ogeom_algo::make_cylinder(model, Frame::WORLD, 2.0, 5.0, T)
        .unwrap()
        .shape;
    let wall = faces(model, &drum)
        .into_iter()
        .find(|f| {
            let Some(NodeData::Face(data)) = model.node(f).map(|n| n.data()) else {
                return false;
            };
            matches!(
                model.geometry().surface(data.surface).unwrap(),
                SurfaceGeometry::Cylinder(_)
            )
        })
        .unwrap();
    (drum, wall)
}

/// The solid is valid and holds the drum's volume.
fn holds_the_drum(model: &Model, shape: &Shape) {
    let diagnosis = ogeom_algo::check(model, shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let volume = ogeom_algo::volume_properties(model, shape, fine(), T)
        .unwrap()
        .mass;
    assert!((volume - 20.0 * PI).abs() <= 20.0 * PI * 1e-5, "{volume}");
}

/// The areas of the pieces `face` became, smallest first.
fn piece_areas(model: &Model, built: &ogeom_algo::Built, face: &Shape) -> Vec<f64> {
    let mut areas: Vec<f64> = built
        .history
        .modified(face)
        .iter()
        .map(|f| area(model, f))
        .collect();
    areas.sort_by(f64::total_cmp);
    areas
}

#[test]
fn a_circle_dropped_onto_the_drum_wall_cuts_it_in_two() {
    let mut model = Model::new();
    let (drum, wall) = drum(&mut model);
    let whole = area(&model, &wall);
    // A circle of radius 3 about the axis, dropped onto the wall.
    let ring = circle(
        &mut model,
        Point::new(0.0, 0.0, 1.5),
        Vector::new(1.0, 0.0, 0.0),
        3.0,
    );
    let built = split_face(
        &mut model,
        &drum,
        &wall,
        &[ring],
        Projection::AlongNormals,
        T,
    )
    .unwrap();
    assert_eq!(faces(&model, &built.shape).len(), 4);
    holds_the_drum(&model, &built.shape);
    let areas = piece_areas(&model, &built, &wall);
    assert_eq!(areas.len(), 2);
    let (low, high) = (2.0 * PI * 2.0 * 1.5, 2.0 * PI * 2.0 * 3.5);
    assert!((whole - (low + high)).abs() <= whole * 1e-5);
    assert!((areas[0] - low).abs() <= low * 1e-5, "{areas:?}");
    assert!((areas[1] - high).abs() <= high * 1e-5, "{areas:?}");
}

#[test]
fn a_circle_round_the_drum_across_its_seam_cuts_the_wall_in_two() {
    // The circle starts a quarter turn from the wall's seam, so its trace
    // leaves the chart's window there and comes back a period on.
    let mut model = Model::new();
    let (drum, wall) = drum(&mut model);
    let ring = circle(
        &mut model,
        Point::new(0.0, 0.0, 2.0),
        Vector::new(0.0, 1.0, 0.0),
        2.0,
    );
    let built = split_face(&mut model, &drum, &wall, &[ring], Projection::OnFace, T).unwrap();
    assert_eq!(faces(&model, &built.shape).len(), 4);
    holds_the_drum(&model, &built.shape);
    let areas = piece_areas(&model, &built, &wall);
    let (low, high) = (2.0 * PI * 2.0 * 2.0, 2.0 * PI * 2.0 * 3.0);
    assert_eq!(areas.len(), 2);
    assert!((areas[0] - low).abs() <= low * 1e-5, "{areas:?}");
    assert!((areas[1] - high).abs() <= high * 1e-5, "{areas:?}");
}

#[test]
fn a_ruling_halves_the_drum_wall_with_its_seam() {
    // Bottom to top opposite the seam: with the seam, the wall's halves.
    let mut model = Model::new();
    let (drum, wall) = drum(&mut model);
    let line = segment(
        &mut model,
        Point::new(-2.0, 0.0, 0.0),
        Point::new(-2.0, 0.0, 5.0),
    );
    let built = split_face(&mut model, &drum, &wall, &[line], Projection::OnFace, T).unwrap();
    holds_the_drum(&model, &built.shape);
    let areas = piece_areas(&model, &built, &wall);
    assert_eq!(areas.len(), 2);
    for a in areas {
        assert!((a - 10.0 * PI).abs() <= 10.0 * PI * 1e-5, "{a}");
    }
}

#[test]
fn a_slanted_line_dropped_across_the_drum_seam_cuts_two_corners() {
    // The line x = 3 from (y, z) = (-3, 0) to (3, 5), dropped along the
    // normals, lands at angle a = atan(y / 3) and height z, so on
    // v(a) = 2.5 (tan a + 1), running from the bottom at a = -pi/4 to the
    // top at a = pi/4 through the seam at a = 0. With the seam it cuts off
    // the region under it on the seam's one side and over it on the other,
    // each 5 (pi/4 - ln sqrt 2) by integrating 2 v(a) and 2 (5 - v(a)).
    let mut model = Model::new();
    let (drum, wall) = drum(&mut model);
    let line = segment(
        &mut model,
        Point::new(3.0, -3.0, 0.0),
        Point::new(3.0, 3.0, 5.0),
    );
    let built = split_face(
        &mut model,
        &drum,
        &wall,
        &[line],
        Projection::AlongNormals,
        T,
    )
    .unwrap();
    holds_the_drum(&model, &built.shape);
    let areas = piece_areas(&model, &built, &wall);
    let corner = 5.0 * (PI / 4.0 - 2.0_f64.sqrt().ln());
    assert_eq!(areas.len(), 3);
    assert!((areas[0] - corner).abs() <= corner * 1e-5, "{areas:?}");
    assert!((areas[1] - corner).abs() <= corner * 1e-5, "{areas:?}");
    let rest = 20.0 * PI - 2.0 * corner;
    assert!((areas[2] - rest).abs() <= rest * 1e-5, "{areas:?}");
}

/// The square with the disk of radius 2 about its middle cut out: the
/// shell of the two, and the holed piece.
fn holed_square(model: &mut Model) -> (Shape, Shape) {
    let face = square(model);
    let ring = circle(
        model,
        Point::new(5.0, 5.0, 0.0),
        Vector::new(1.0, 0.0, 0.0),
        2.0,
    );
    let built = split_face(model, &face, &face, &[ring], Projection::OnFace, T).unwrap();
    let holed = faces(model, &built.shape)
        .into_iter()
        .find(|f| wires(model, f) == 2)
        .unwrap();
    (built.shape, holed)
}

#[test]
fn a_circle_round_a_hole_cuts_out_a_ring() {
    let mut model = Model::new();
    let (shell, holed) = holed_square(&mut model);
    let ring = circle(
        &mut model,
        Point::new(5.0, 5.0, 0.0),
        Vector::new(1.0, 0.0, 0.0),
        3.0,
    );
    let built = split_face(&mut model, &shell, &holed, &[ring], Projection::OnFace, T).unwrap();
    let areas = piece_areas(&model, &built, &holed);
    assert_eq!(areas.len(), 2);
    let annulus = 9.0 * PI - 4.0 * PI;
    assert!((areas[0] - annulus).abs() <= annulus * 1e-5, "{areas:?}");
    let rest = 100.0 - 9.0 * PI;
    assert!((areas[1] - rest).abs() <= rest * 1e-5, "{areas:?}");
    for piece in built.history.modified(&holed) {
        assert_eq!(wires(&model, piece), 2, "both pieces hold a hole");
    }
}

#[test]
fn a_circle_crossing_a_hole_cuts_off_the_part_outside_it() {
    // A circle of radius 2 two above the hole's middle crosses the hole's
    // edge twice; its part outside the hole is a new piece, the disk less
    // the lens the two disks share, 8 pi / 3 - sqrt 12.
    let mut model = Model::new();
    let (shell, holed) = holed_square(&mut model);
    let ring = circle(
        &mut model,
        Point::new(5.0, 7.0, 0.0),
        Vector::new(1.0, 0.0, 0.0),
        2.0,
    );
    let built = split_face(&mut model, &shell, &holed, &[ring], Projection::OnFace, T).unwrap();
    let areas = piece_areas(&model, &built, &holed);
    let lens = 8.0 * PI / 3.0 - 12.0_f64.sqrt();
    let crescent = 4.0 * PI - lens;
    let rest = 100.0 - 4.0 * PI - crescent;
    assert_eq!(areas.len(), 2);
    assert!((areas[0] - crescent).abs() <= crescent * 1e-5, "{areas:?}");
    assert!((areas[1] - rest).abs() <= rest * 1e-5, "{areas:?}");
}
