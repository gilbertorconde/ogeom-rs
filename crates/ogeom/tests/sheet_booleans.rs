//! Booleans on sheets: a face or an open shell trimmed by a solid or a half
//! space, and split where another sheet crosses it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, make_box, make_edge, make_face, make_half_space, make_natural_face, make_polygon,
    make_prism, make_wire, shape_bounds, surface_properties,
};
use ogeom::boolean::{common, cut, split_sheet};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Aabb, Circle, Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// A planar face bounded by the closed polygon through `corners`.
fn polygon_face(model: &mut Model, corners: &[Point]) -> Shape {
    let wire = make_polygon(model, corners, true, T).unwrap().shape;
    let plane = Plane::through_points(corners[0], corners[1], corners[2], T).unwrap();
    make_face(model, PlaneSurface::new(plane).into(), &[wire], T)
        .unwrap()
        .shape
}

/// The 10 x 10 square on z = 0 with a corner at the origin.
fn square(model: &mut Model) -> Shape {
    polygon_face(
        model,
        &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)].map(|(x, y)| Point::new(x, y, 0.0)),
    )
}

/// The half space x <= `at`.
fn below_x(model: &mut Model, at: f64) -> Shape {
    let plane = Plane::through(Point::new(at, 0.0, 0.0), Direction::X);
    let surface = PlaneSurface::over(plane, (-1e3, 1e3), (-1e3, 1e3)).unwrap();
    let face = make_natural_face(model, SurfaceGeometry::Plane(surface))
        .unwrap()
        .shape;
    make_half_space(model, &face, Point::new(at - 1.0, 0.0, 0.0), T)
        .unwrap()
        .shape
}

/// The segment from `from` to `to` swept up by `height`: a planar sheet.
fn wall(model: &mut Model, from: Point, to: Point, height: f64) -> Shape {
    let wire = make_polygon(model, &[from, to], false, T).unwrap().shape;
    make_prism(model, &wire, Vector::new(0.0, 0.0, height), T)
        .unwrap()
        .shape
}

/// The circle of `radius` round the z axis swept up by `height`: an open
/// tube, one cylindrical face with a seam.
fn tube(model: &mut Model, radius: f64, height: f64) -> Shape {
    let circle = Curve::Circle(CircleCurve::new(
        Circle::new(Frame::WORLD, radius, T).unwrap(),
    ));
    let edge = make_edge(model, circle, (0.0, 2.0 * PI), T).unwrap().shape;
    let wire = make_wire(model, &[edge], T).unwrap().shape;
    make_prism(model, &wire, Vector::new(0.0, 0.0, height), T)
        .unwrap()
        .shape
}

fn faces(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Face).unwrap()
}

/// Area integrated on the exact surfaces.
fn area(model: &Model, shape: &Shape) -> f64 {
    let found = surface_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T).unwrap();
    assert_eq!(found.deflection, 0.0, "integrated on the exact surfaces");
    found.mass
}

fn bounds(model: &Model, shape: &Shape) -> (Point, Point) {
    let b: Aabb = shape_bounds(model, shape, T).unwrap();
    (b.low().unwrap(), b.high().unwrap())
}

/// Whether the checker finds nothing wrong past the openness a sheet has.
fn usable(model: &Model, shape: &Shape) {
    let found = check(model, shape, T).unwrap();
    assert!(found.is_usable(), "{found}");
}

/// The surface a face lies on.
fn surface_of(model: &Model, face: &Shape) -> SurfaceGeometry {
    let data = model.node(face).unwrap().data().as_face().unwrap();
    model.geometry().surface(data.surface).unwrap().clone()
}

/// Every edge of every face carries a pcurve on that face's surface, and
/// the edges two faces share are counted.
fn shared_edges_with_pcurves(model: &Model, shape: &Shape) -> usize {
    let mut users: std::collections::HashMap<ogeom::topo::TShapeId, usize> =
        std::collections::HashMap::new();
    for face in faces(model, shape) {
        let surface = model.node(&face).unwrap().data().as_face().unwrap().surface;
        for edge in explore_unique(model, &face, ShapeType::Edge).unwrap() {
            let data = model.node(&edge).unwrap().data().as_edge().unwrap();
            if data.curve3d().is_some() {
                assert!(
                    data.pcurve_for(surface, edge.location()).is_some(),
                    "an edge with no pcurve on a face it bounds"
                );
            }
            *users.entry(edge.node()).or_default() += 1;
        }
    }
    users.values().filter(|&&n| n == 2).count()
}

#[test]
fn a_square_kept_on_one_side_of_a_half_space_is_one_face_of_half_its_area() {
    let mut model = Model::new();
    let sheet = square(&mut model);
    let half = below_x(&mut model, 5.0);
    let kept = common(&mut model, &sheet, &half, T).unwrap();
    let found = faces(&model, &kept.shape);
    assert_eq!(found.len(), 1);
    assert!((area(&model, &kept.shape) - 50.0).abs() < 1e-9);
    let (lo, hi) = bounds(&model, &found[0]);
    assert!(
        lo.x.abs() < 1e-6 && (hi.x - 5.0).abs() < 1e-6,
        "{lo:?} {hi:?}"
    );
    assert!(lo.z.abs() < 1e-6 && hi.z.abs() < 1e-6);
    assert!(matches!(
        surface_of(&model, &found[0]),
        SurfaceGeometry::Plane(_)
    ));
    // The piece is recorded against the face it came from.
    assert_eq!(kept.history.trace(&sheet), &[found[0].clone()]);
    usable(&model, &kept.shape);

    let away = cut(&mut model, &sheet, &half, T).unwrap();
    let found = faces(&model, &away.shape);
    assert_eq!(found.len(), 1);
    assert!((area(&model, &away.shape) - 50.0).abs() < 1e-9);
    let (lo, hi) = bounds(&model, &found[0]);
    assert!(
        (lo.x - 5.0).abs() < 1e-6 && (hi.x - 10.0).abs() < 1e-6,
        "{lo:?} {hi:?}"
    );
}

#[test]
fn a_sheet_trimmed_by_a_box_keeps_the_patch_inside_or_the_frame_round_it() {
    let mut model = Model::new();
    let sheet = square(&mut model);
    let frame = Frame::new(Point::new(2.0, 3.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let block = make_box(&mut model, frame, (4.0, 5.0, 2.0), T)
        .unwrap()
        .shape;

    let inside = common(&mut model, &sheet, &block, T).unwrap();
    assert_eq!(faces(&model, &inside.shape).len(), 1);
    assert!((area(&model, &inside.shape) - 20.0).abs() < 1e-9);
    let (lo, hi) = bounds(&model, &inside.shape);
    assert!(lo.distance(Point::new(2.0, 3.0, 0.0)) < 1e-6, "{lo:?}");
    assert!(hi.distance(Point::new(6.0, 8.0, 0.0)) < 1e-6, "{hi:?}");

    // Common takes the sheet on either side.
    let swapped = common(&mut model, &block, &sheet, T).unwrap();
    assert!((area(&model, &swapped.shape) - 20.0).abs() < 1e-9);

    // Outside, the square with a 4 x 5 hole: one face, two wires.
    let outside = cut(&mut model, &sheet, &block, T).unwrap();
    let found = faces(&model, &outside.shape);
    assert_eq!(found.len(), 1);
    assert_eq!(model.children_of(&found[0]).unwrap().len(), 2);
    assert!((area(&model, &outside.shape) - 80.0).abs() < 1e-9);
    usable(&model, &outside.shape);
}

#[test]
fn a_sheet_cut_in_two_apart_is_a_compound_of_two_shells() {
    let mut model = Model::new();
    let sheet = square(&mut model);
    // A bar across the square from x = 4 to 6, past it in y and z.
    let frame = Frame::new(Point::new(4.0, -1.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let bar = make_box(&mut model, frame, (2.0, 12.0, 2.0), T)
        .unwrap()
        .shape;
    let apart = cut(&mut model, &sheet, &bar, T).unwrap();
    assert_eq!(model.kind_of(&apart.shape).unwrap(), ShapeType::Compound);
    let shells = model.children_of(&apart.shape).unwrap();
    assert_eq!(shells.len(), 2);
    let mut areas: Vec<f64> = shells.iter().map(|s| area(&model, s)).collect();
    areas.sort_by(f64::total_cmp);
    assert!((areas[0] - 40.0).abs() < 1e-9 && (areas[1] - 40.0).abs() < 1e-9);
    let mut spans: Vec<(f64, f64)> = shells
        .iter()
        .map(|s| {
            let (lo, hi) = bounds(&model, s);
            (lo.x, hi.x)
        })
        .collect();
    spans.sort_by(|a, b| a.0.total_cmp(&b.0));
    assert!(spans[0].0.abs() < 1e-6 && (spans[0].1 - 4.0).abs() < 1e-6);
    assert!((spans[1].0 - 6.0).abs() < 1e-6 && (spans[1].1 - 10.0).abs() < 1e-6);
    // The one face is split into both: it traces to two faces.
    assert_eq!(apart.history.trace(&sheet).len(), 2);
}

#[test]
fn a_sheet_wholly_inside_the_tool_is_kept_whole_by_common_and_removed_by_cut() {
    let mut model = Model::new();
    let sheet = square(&mut model);
    let half = below_x(&mut model, 20.0);
    let kept = common(&mut model, &sheet, &half, T).unwrap();
    assert_eq!(faces(&model, &kept.shape).len(), 1);
    assert!((area(&model, &kept.shape) - 100.0).abs() < 1e-9);
    let gone = cut(&mut model, &sheet, &half, T).unwrap();
    assert!(faces(&model, &gone.shape).is_empty());
    assert!(gone.history.is_deleted(&sheet));
}

/// The plane through the z axis at `angle` from the x axis, as a face past
/// a tube of `height` at both ends.
fn axial_plane(model: &mut Model, angle: f64, height: f64) -> Shape {
    let (c, s) = (5.0 * angle.cos(), 5.0 * angle.sin());
    polygon_face(
        model,
        &[
            Point::new(-c, -s, -1.0),
            Point::new(c, s, -1.0),
            Point::new(c, s, height + 1.0),
            Point::new(-c, -s, height + 1.0),
        ],
    )
}

/// The tube's pieces split by the plane at `angle`, each checked to lie on
/// the cylinder on one side of the plane, with the area on each side.
fn tube_halves(angle: f64) -> (usize, [f64; 2], usize) {
    let (radius, height) = (2.0, 7.0);
    let mut model = Model::new();
    let sheet = tube(&mut model, radius, height);
    assert_eq!(model.kind_of(&sheet).unwrap(), ShapeType::Shell);
    let source = faces(&model, &sheet)[0].clone();
    let by = axial_plane(&mut model, angle, height);
    let split = split_sheet(&mut model, &sheet, &by, T).unwrap();
    let found = faces(&model, &split.shape);
    // The plane's normal: each piece stands on one side of it.
    let normal = Vector::new(-angle.sin(), angle.cos(), 0.0);
    let mut sides = [0.0, 0.0];
    for face in &found {
        assert!(matches!(
            surface_of(&model, face),
            SurfaceGeometry::Cylinder(_)
        ));
        let (lo, hi) = bounds(&model, face);
        assert!(lo.z.abs() < 1e-6 && (hi.z - height).abs() < 1e-6);
        let mesh = ogeom::mesh::triangulate_face(&model, face, Deflection::default(), T).unwrap();
        let reach: Vec<f64> = mesh
            .positions
            .iter()
            .map(|p| (*p - Point::ORIGIN).dot(normal))
            .collect();
        let above = reach.iter().all(|&d| d > -1e-6);
        let below = reach.iter().all(|&d| d < 1e-6);
        assert!(above != below, "a piece on both sides of the plane");
        sides[usize::from(above)] += area(&model, face);
    }
    let traced: std::collections::HashSet<Shape> =
        split.history.trace(&source).iter().cloned().collect();
    assert_eq!(traced, found.iter().cloned().collect());
    usable(&model, &split.shape);
    (
        found.len(),
        sides,
        shared_edges_with_pcurves(&model, &split.shape),
    )
}

#[test]
fn a_tube_split_by_a_plane_through_its_axis_is_two_half_tubes() {
    let half = PI * 2.0 * 7.0;
    // Through the seam, square to it, and at angles between: two faces,
    // half each, the two lines the plane cuts along shared by both with a
    // pcurve on the cylinder for each. Where the seam falls inside a half,
    // that half is one face across it.
    for angle in [0.0, PI / 2.0, PI / 3.0, 0.75 * PI, 1.0, 2.5] {
        let (count, sides, shared) = tube_halves(angle);
        assert_eq!(count, 2, "at {angle}");
        assert_eq!(shared, 2, "at {angle}");
        for side in sides {
            assert!((side - half).abs() < 1e-9 * half, "at {angle}: {sides:?}");
        }
    }
}

#[test]
fn two_crossing_walls_split_each_other_into_four_faces() {
    let mut model = Model::new();
    let first = wall(
        &mut model,
        Point::new(-5.0, 0.0, 0.0),
        Point::new(5.0, 0.0, 0.0),
        10.0,
    );
    let second = wall(
        &mut model,
        Point::new(1.0, -4.0, 0.0),
        Point::new(1.0, 6.0, 0.0),
        10.0,
    );
    // Each is split by the other, in a call of its own.
    let a = split_sheet(&mut model, &first, &second, T).unwrap();
    let b = split_sheet(&mut model, &second, &first, T).unwrap();
    let mut areas: Vec<f64> = faces(&model, &a.shape)
        .iter()
        .chain(faces(&model, &b.shape).iter())
        .map(|f| area(&model, f))
        .collect();
    areas.sort_by(f64::total_cmp);
    let expected = [40.0, 40.0, 60.0, 60.0];
    assert_eq!(areas.len(), 4);
    for (got, want) in areas.iter().zip(expected) {
        assert!((got - want).abs() < 1e-9, "{areas:?}");
    }
    // Each split is one shell of two faces sharing the line x = 1, y = 0.
    for split in [&a, &b] {
        assert_eq!(model.kind_of(&split.shape).unwrap(), ShapeType::Shell);
        assert_eq!(shared_edges_with_pcurves(&model, &split.shape), 1);
        usable(&model, &split.shape);
    }
}

#[test]
fn an_open_shell_is_split_only_on_the_faces_the_other_sheet_crosses() {
    let mut model = Model::new();
    // An L of two walls: y = 0 from x = 0 to 10, then x = 10 from y = 0 to 6.
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(10.0, 0.0, 0.0),
        Point::new(10.0, 6.0, 0.0),
    ];
    let wire = make_polygon(&mut model, &corners, false, T).unwrap().shape;
    let sheet = make_prism(&mut model, &wire, Vector::new(0.0, 0.0, 4.0), T)
        .unwrap()
        .shape;
    let untouched = faces(&model, &sheet)
        .into_iter()
        .find(|f| bounds(&model, f).0.x > 9.0)
        .unwrap();
    let by = wall(
        &mut model,
        Point::new(3.0, -2.0, -1.0),
        Point::new(3.0, 2.0, -1.0),
        6.0,
    );
    let split = split_sheet(&mut model, &sheet, &by, T).unwrap();
    let mut areas: Vec<f64> = faces(&model, &split.shape)
        .iter()
        .map(|f| area(&model, f))
        .collect();
    areas.sort_by(f64::total_cmp);
    let expected = [12.0, 24.0, 28.0];
    assert_eq!(areas.len(), 3, "{areas:?}");
    for (got, want) in areas.iter().zip(expected) {
        assert!((got - want).abs() < 1e-9, "{areas:?}");
    }
    // The corner edge and the new line are each shared by two faces.
    assert_eq!(shared_edges_with_pcurves(&model, &split.shape), 2);
    assert_eq!(split.history.trace(&untouched).len(), 1);
    usable(&model, &split.shape);
}

#[test]
fn a_sheet_crossing_partway_separates_nothing() {
    let mut model = Model::new();
    let sheet = square(&mut model);
    // The wall x = 5 reaches from below the square to y = 6 inside it: a
    // slit, which leaves the square one piece.
    let by = wall(
        &mut model,
        Point::new(5.0, -2.0, -1.0),
        Point::new(5.0, 6.0, -1.0),
        2.0,
    );
    let split = split_sheet(&mut model, &sheet, &by, T).unwrap();
    assert_eq!(faces(&model, &split.shape).len(), 1);
    assert!((area(&model, &split.shape) - 100.0).abs() < 1e-9);
}

#[test]
fn sheet_booleans_refuse_what_they_do_not_resolve() {
    let mut model = Model::new();
    let sheet = square(&mut model);
    let other = wall(
        &mut model,
        Point::new(5.0, -1.0, -1.0),
        Point::new(5.0, 11.0, -1.0),
        2.0,
    );
    let block = make_box(&mut model, Frame::WORLD, (4.0, 4.0, 4.0), T)
        .unwrap()
        .shape;
    let said = |r: ogeom::core::OgeomResult<ogeom::algo::Built>| r.unwrap_err().to_string();

    // A sheet bounds no volume to cut away, and two share none.
    assert!(said(cut(&mut model, &sheet, &other, T)).contains("a sheet bounds no volume"));
    assert!(said(cut(&mut model, &block, &sheet, T)).contains("a sheet bounds no volume"));
    assert!(said(common(&mut model, &sheet, &other, T)).contains("two sheets share no volume"));
    // A union with a sheet is not a sheet operation.
    assert!(
        said(ogeom::boolean::fuse(&mut model, &sheet, &block, T)).contains("a union takes solids")
    );
    // A solid splits nothing through split_sheet.
    assert!(said(split_sheet(&mut model, &sheet, &block, T)).contains("splits a sheet by a sheet"));
    assert!(said(split_sheet(&mut model, &block, &sheet, T)).contains("splits a sheet by a sheet"));
    // The square lies on the box's bottom face over a 4 x 4 patch.
    assert!(said(common(&mut model, &sheet, &block, T)).contains("lies on the tool's boundary"));
    // Two sheets overlapping on one plane.
    let overlapping = polygon_face(
        &mut model,
        &[(5.0, 5.0), (15.0, 5.0), (15.0, 15.0), (5.0, 15.0)].map(|(x, y)| Point::new(x, y, 0.0)),
    );
    assert!(said(split_sheet(&mut model, &sheet, &overlapping, T)).contains("lie on one surface"));
}
