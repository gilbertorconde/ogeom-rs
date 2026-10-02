//! Rounds on sheets: the edge where two faces of an open shell meet, and
//! the corner between two faces of separate shapes, measured against the
//! rolling ball's closed form.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, edge_vertices, face_normal, make_edge_between, make_face_with_pcurves, make_polygon,
    make_prism, make_vertex, make_wire, surface_properties,
};
use ogeom::core::Tolerances;
use ogeom::fillet::{analyse_blend, fillet_faces, fillet_sheet_edges};
use ogeom::geom::{CircleCurve, Curve, Curve3d as _, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// The open polyline through `corners` (in the xy plane) swept up by
/// `height`: a sheet of planar faces.
fn extruded(model: &mut Model, corners: &[(f64, f64)], travel: Vector) -> Shape {
    let points: Vec<Point> = corners
        .iter()
        .map(|&(x, y)| Point::new(x, y, 0.0))
        .collect();
    let wire = make_polygon(model, &points, false, T).unwrap().shape;
    let sheet = make_prism(model, &wire, travel, T).unwrap().shape;
    assert_eq!(model.kind_of(&sheet).unwrap(), ShapeType::Shell);
    sheet
}

/// The open L: a floor along x and a wall along y, 10 wide and `height`
/// tall, meeting at the z axis.
fn open_l(model: &mut Model, height: f64) -> Shape {
    extruded(
        model,
        &[(10.0, 0.0), (0.0, 0.0), (0.0, 10.0)],
        Vector::new(0.0, 0.0, height),
    )
}

/// A vertex's point, placed.
fn point_of(model: &Model, vertex: &Shape) -> Point {
    let point = model
        .node(vertex)
        .unwrap()
        .data()
        .as_vertex()
        .unwrap()
        .point;
    vertex.transform(model.datums()).unwrap().apply(point)
}

/// The edges of `shape` both of whose ends satisfy `keep`.
fn edges_where(model: &Model, shape: &Shape, keep: impl Fn(Point) -> bool) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = edge_vertices(model, e).unwrap().unwrap();
            keep(point_of(model, &a)) && keep(point_of(model, &b))
        })
        .collect()
}

/// The sheet's edge standing on the vertical line through `(x, y)`.
fn vertical_edge_at(model: &Model, shape: &Shape, x: f64, y: f64) -> Shape {
    let found = edges_where(model, shape, |p| {
        (p.x - x).abs() < 1e-9 && (p.y - y).abs() < 1e-9
    });
    assert_eq!(found.len(), 1, "one edge stands at ({x}, {y})");
    found[0].clone()
}

fn faces(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Face).unwrap()
}

fn surface_of(model: &Model, face: &Shape) -> SurfaceGeometry {
    let data = model.node(face).unwrap().data().as_face().unwrap();
    model.geometry().surface(data.surface).unwrap().clone()
}

/// The one cylindrical face of `shape`.
fn the_round(model: &Model, shape: &Shape) -> (Shape, ogeom::math::Cylinder) {
    let mut rounds = faces(model, shape)
        .into_iter()
        .filter_map(|f| match surface_of(model, &f) {
            SurfaceGeometry::Cylinder(c) => Some((f, c.cylinder())),
            _ => None,
        });
    let found = rounds.next().expect("a cylindrical round");
    assert!(rounds.next().is_none(), "one round");
    found
}

/// Area integrated on the exact surfaces.
fn area(model: &Model, shape: &Shape) -> f64 {
    let found = surface_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T).unwrap();
    assert_eq!(found.deflection, 0.0, "integrated on the exact surfaces");
    found.mass
}

/// How many faces use each edge: (shared by two, used by one).
fn edge_use(model: &Model, shape: &Shape) -> (usize, usize) {
    let mut users: std::collections::HashMap<ogeom::topo::TShapeId, usize> =
        std::collections::HashMap::new();
    for face in faces(model, shape) {
        for edge in explore_unique(model, &face, ShapeType::Edge).unwrap() {
            *users.entry(edge.node()).or_default() += 1;
        }
    }
    let shared = users.values().filter(|&&n| n == 2).count();
    let free = users.values().filter(|&&n| n == 1).count();
    assert_eq!(shared + free, users.len(), "no edge is used three times");
    (shared, free)
}

/// Whether the checker finds nothing wrong past the openness a sheet has.
fn usable(model: &Model, shape: &Shape) {
    let found = check(model, shape, T).unwrap();
    assert!(found.is_usable(), "{found}");
}

/// The distance from `p` to the cylinder's axis line.
fn off_axis(cylinder: &ogeom::math::Cylinder, p: Point) -> f64 {
    let axis = cylinder.axis();
    let d = p - axis.location;
    let along = axis.direction.vector();
    (d - along * d.dot(along)).magnitude()
}

/// The round's tangency and curvature against each face it shares an
/// edge with: the tangency error is zero and the curvature jump across
/// the edge is the cylinder's own `1 / r`, the plane's being nothing.
fn tangent_to_both(model: &Model, shape: &Shape, round: &Shape, radius: f64) {
    let contacts = analyse_blend(model, shape, round, 15, T).unwrap();
    assert_eq!(contacts.len(), 2, "the round meets two faces: {contacts:?}");
    for contact in &contacts {
        assert!(
            contact.tangency_error < 0.1_f64.to_radians(),
            "tangent within a tenth of a degree: {}",
            contact.tangency_error
        );
        assert!(
            contact.tangency_error < 1e-9,
            "tangent to rounding: {}",
            contact.tangency_error
        );
        assert!(
            contact.gap < 1e-9,
            "the shared edge lies on both: {}",
            contact.gap
        );
        assert!(
            (contact.curvature_error - 1.0 / radius).abs() < 1e-6,
            "the curvature jumps by the cylinder's 1/r: {}",
            contact.curvature_error
        );
    }
}

/// The acceptance case: an extruded open L with its corner rounded at
/// radius 2 is three faces, the round a quarter cylinder of radius 2 on
/// the axis through (2, 2), tangent to both flats.
#[test]
fn an_open_l_rounds_to_a_quarter_cylinder() {
    let mut model = Model::new();
    let (r, h) = (2.0, 5.0);
    let sheet = open_l(&mut model, h);
    let corner = vertical_edge_at(&model, &sheet, 0.0, 0.0);
    let built =
        fillet_sheet_edges(&mut model, &sheet, std::slice::from_ref(&corner), r, T).unwrap();
    let rounded = built.shape.clone();

    assert_eq!(model.kind_of(&rounded).unwrap(), ShapeType::Shell);
    assert_eq!(faces(&model, &rounded).len(), 3, "two flats and the round");
    usable(&model, &rounded);
    // Two shared edges, the lines of contact; the rest is free boundary:
    // three edges round each flat's other sides and two arcs at the ends.
    assert_eq!(edge_use(&model, &rounded), (2, 8));

    let (round, cylinder) = the_round(&model, &rounded);
    assert!(
        (cylinder.radius() - r).abs() < 1e-12,
        "radius {}",
        cylinder.radius()
    );
    let axis = cylinder.axis();
    assert!(
        axis.direction
            .vector()
            .cross(Vector::new(0.0, 0.0, 1.0))
            .magnitude()
            < 1e-12,
        "the axis runs along the edge"
    );
    assert!(
        off_axis(&cylinder, Point::new(r, r, 0.0)) < 1e-12,
        "the axis passes through (r, r)"
    );
    // A quarter turn: the area is a quarter of the cylinder's over the height.
    let round_area = area(&model, &round);
    assert!(
        (round_area - PI / 2.0 * r * h).abs() < 1e-9,
        "a quarter cylinder: {round_area} against {}",
        PI / 2.0 * r * h
    );
    // Each flat lost the strip r wide beside the corner.
    for face in faces(&model, &rounded) {
        if face.is_same(&round) {
            continue;
        }
        let a = area(&model, &face);
        assert!((a - (10.0 - r) * h).abs() < 1e-9, "a trimmed flat: {a}");
    }
    // Every vertex of the round is a line of contact's end, on a flat.
    for vertex in explore_unique(&model, &round, ShapeType::Vertex).unwrap() {
        let p = point_of(&model, &vertex);
        let on_floor = p.y.abs() < 1e-12 && (p.x - r).abs() < 1e-12;
        let on_wall = p.x.abs() < 1e-12 && (p.y - r).abs() < 1e-12;
        assert!(on_floor || on_wall, "a contact at {p:?}");
    }
    tangent_to_both(&model, &rounded, &round, r);

    // The round faces the way the sheet does: the prism's flats face into
    // the corner, so the round faces its own axis.
    let (p, n) = face_normal(&model, &round, T).unwrap();
    let axis_point =
        axis.location + axis.direction.vector() * (p - axis.location).dot(axis.direction.vector());
    let to_axis = (axis_point - p) / r;
    assert!((n - to_axis).magnitude() < 1e-9, "faces its axis: {n:?}");
    for face in faces(&model, &rounded) {
        if face.is_same(&round) {
            continue;
        }
        let (q, m) = face_normal(&model, &face, T).unwrap();
        let into_corner = if q.y.abs() < 1e-12 {
            Vector::new(0.0, 1.0, 0.0)
        } else {
            Vector::new(1.0, 0.0, 0.0)
        };
        assert!(
            (m - into_corner).magnitude() < 1e-12,
            "a flat faces the corner"
        );
    }

    // History: the edge gave the round and is gone; both flats were
    // rebuilt; the sheet became the result.
    assert!(built.history.is_deleted(&corner));
    assert!(
        built
            .history
            .generated(&corner)
            .iter()
            .any(|f| f.node() == round.node())
    );
    assert!(
        built
            .history
            .modified(&sheet)
            .iter()
            .any(|s| s.is_same(&rounded))
    );
}

/// A corner that is not square: the wall leaves the floor at 120 degrees,
/// so the ball sets back `r / tan 60°` along each face, sits `r / sin 60°`
/// from the corner, and its arc sweeps a sixth of a turn.
#[test]
fn an_obtuse_corner_sets_back_by_the_half_angle() {
    let mut model = Model::new();
    let (r, h) = (1.5, 4.0);
    let (c, s) = (
        (120.0_f64).to_radians().cos(),
        (120.0_f64).to_radians().sin(),
    );
    let sheet = extruded(
        &mut model,
        &[(10.0, 0.0), (0.0, 0.0), (10.0 * c, 10.0 * s)],
        Vector::new(0.0, 0.0, h),
    );
    let corner = vertical_edge_at(&model, &sheet, 0.0, 0.0);
    let rounded = fillet_sheet_edges(&mut model, &sheet, &[corner], r, T)
        .unwrap()
        .shape;
    usable(&model, &rounded);
    let (round, cylinder) = the_round(&model, &rounded);
    let half = 60.0_f64.to_radians();
    let setback = r / half.tan();
    let depth = r / half.sin();
    let bisector = Vector::new(
        60.0_f64.to_radians().cos(),
        60.0_f64.to_radians().sin(),
        0.0,
    );
    let centre = Point::new(0.0, 0.0, 0.0) + bisector * depth;
    assert!(
        off_axis(&cylinder, centre) < 1e-12,
        "the axis at the closed form"
    );
    let round_area = area(&model, &round);
    assert!(
        (round_area - r * (PI / 3.0) * h).abs() < 1e-9,
        "a sixth of a turn: {round_area}"
    );
    for face in faces(&model, &rounded) {
        if face.is_same(&round) {
            continue;
        }
        let a = area(&model, &face);
        assert!(
            (a - (10.0 - setback) * h).abs() < 1e-9,
            "each flat loses the setback: {a} against {}",
            (10.0 - setback) * h
        );
    }
    tangent_to_both(&model, &rounded, &round, r);
}

/// Two corners of a channel rounded in one call: the floor between them
/// is trimmed from both sides, and the history follows both edges.
#[test]
fn a_channel_rounds_both_corners_in_one_call() {
    let mut model = Model::new();
    let (r, h) = (1.0, 3.0);
    let sheet = extruded(
        &mut model,
        &[(0.0, 10.0), (0.0, 0.0), (6.0, 0.0), (6.0, 10.0)],
        Vector::new(0.0, 0.0, h),
    );
    let left = vertical_edge_at(&model, &sheet, 0.0, 0.0);
    let right = vertical_edge_at(&model, &sheet, 6.0, 0.0);
    let built =
        fillet_sheet_edges(&mut model, &sheet, &[left.clone(), right.clone()], r, T).unwrap();
    let rounded = built.shape.clone();
    usable(&model, &rounded);
    assert_eq!(faces(&model, &rounded).len(), 5);
    assert_eq!(edge_use(&model, &rounded).0, 4, "four lines of contact");
    let expected = 2.0 * (10.0 - r) * h + (6.0 - 2.0 * r) * h + 2.0 * (PI / 2.0 * r * h);
    let total = area(&model, &rounded);
    assert!(
        (total - expected).abs() < 1e-9,
        "{total} against {expected}"
    );
    for edge in [&left, &right] {
        assert!(built.history.is_deleted(edge));
        assert_eq!(built.history.generated(edge).len(), 1, "one round each");
    }
    for edge in [&left, &right] {
        let round = &built.history.generated(edge)[0];
        let contacts = analyse_blend(&model, &rounded, round, 9, T).unwrap();
        assert_eq!(contacts.len(), 2);
        assert!(contacts.iter().all(|c| c.tangency_error < 1e-9));
    }
}

/// Two separate planes tilted 30 degrees up from either side of a valley,
/// their normals facing into it, 8 long along y.
fn valley(model: &mut Model, near: f64, a_span: (f64, f64), b_span: (f64, f64)) -> (Shape, Shape) {
    let k = 30.0_f64.to_radians().tan();
    let face = |model: &mut Model, corners: [Point; 4]| {
        let wire = make_polygon(model, &corners, true, T).unwrap().shape;
        let edges = model.children_of(&wire).unwrap();
        let plane = Plane::through_points(corners[0], corners[1], corners[2], T).unwrap();
        make_face_with_pcurves(model, PlaneSurface::new(plane).into(), &[edges], T)
            .unwrap()
            .shape
    };
    let a = face(
        model,
        [
            Point::new(near, a_span.0, near * k),
            Point::new(10.0, a_span.0, 10.0 * k),
            Point::new(10.0, a_span.1, 10.0 * k),
            Point::new(near, a_span.1, near * k),
        ],
    );
    let b = face(
        model,
        [
            Point::new(-10.0, b_span.0, 10.0 * k),
            Point::new(-near, b_span.0, near * k),
            Point::new(-near, b_span.1, near * k),
            Point::new(-10.0, b_span.1, 10.0 * k),
        ],
    );
    for f in [&a, &b] {
        assert!(
            face_normal(model, f, T).unwrap().1.z > 0.0,
            "faces into the valley"
        );
    }
    (a, b)
}

/// The acceptance case for separate faces: two tilted planes rounded with
/// `trim` sew into one shell of three faces, the round tangent to both.
///
/// The valley opens 120 degrees; a ball of radius 2 sits `2 / sin 60°`
/// above its floor line and touches each plane `2 / tan 60°` up it, at
/// `x = ±1`.
#[test]
fn two_tilted_planes_trim_into_one_shell() {
    let mut model = Model::new();
    let r = 2.0;
    let (a, b) = valley(&mut model, 0.5, (0.0, 8.0), (0.0, 8.0));
    let built = fillet_faces(&mut model, &a, &b, r, true, T).unwrap();
    let shell = built.shape.clone();
    assert_eq!(model.kind_of(&shell).unwrap(), ShapeType::Shell);
    assert_eq!(faces(&model, &shell).len(), 3);
    usable(&model, &shell);
    // One shell: the round shares a whole edge with each trimmed face.
    let three = faces(&model, &shell);
    let sewn = ogeom::algo::sew(&mut model, &three, T).unwrap();
    assert_eq!(sewn.shells.len(), 1, "the three faces sew into one shell");
    assert_eq!(edge_use(&model, &shell), (2, 8));

    let (round, cylinder) = the_round(&model, &shell);
    let depth = r / 60.0_f64.to_radians().sin();
    assert!((cylinder.radius() - r).abs() < 1e-12);
    assert!(off_axis(&cylinder, Point::new(0.0, 0.0, depth)) < 1e-9);
    assert!(off_axis(&cylinder, Point::new(0.0, 5.0, depth)) < 1e-9);
    let round_area = area(&model, &round);
    assert!(
        (round_area - r * (PI / 3.0) * 8.0).abs() < 1e-9,
        "{round_area}"
    );
    let cos30 = 30.0_f64.to_radians().cos();
    let setback = r / 60.0_f64.to_radians().tan();
    for face in faces(&model, &shell) {
        if face.is_same(&round) {
            continue;
        }
        let a = area(&model, &face);
        let expected = (10.0 / cos30 - setback) * 8.0;
        assert!(
            (a - expected).abs() < 1e-9,
            "a trimmed plane: {a} against {expected}"
        );
    }
    for vertex in explore_unique(&model, &round, ShapeType::Vertex).unwrap() {
        let p = point_of(&model, &vertex);
        assert!(
            (p.x.abs() - setback * cos30).abs() < 1e-9,
            "a contact at {p:?}"
        );
    }
    tangent_to_both(&model, &shell, &round, r);
    assert!(
        built
            .history
            .generated(&a)
            .iter()
            .any(|f| f.node() == round.node())
    );
    assert_eq!(built.history.modified(&a).len(), 1);
}

/// Without `trim` the faces are left alone and the round comes back by
/// itself, its lines of contact on the two planes.
#[test]
fn untrimmed_faces_are_left_as_they_are() {
    let mut model = Model::new();
    let r = 2.0;
    let (a, b) = valley(&mut model, 0.5, (0.0, 8.0), (0.0, 8.0));
    let before = (area(&model, &a), area(&model, &b));
    let built = fillet_faces(&mut model, &a, &b, r, false, T).unwrap();
    assert_eq!(model.kind_of(&built.shape).unwrap(), ShapeType::Face);
    assert_eq!((area(&model, &a), area(&model, &b)), before);
    let round_area = area(&model, &built.shape);
    assert!((round_area - r * (PI / 3.0) * 8.0).abs() < 1e-9);
    let cos30 = 30.0_f64.to_radians().cos();
    let setback = r / 60.0_f64.to_radians().tan();
    let k = 30.0_f64.to_radians().tan();
    for vertex in explore_unique(&model, &built.shape, ShapeType::Vertex).unwrap() {
        let p = point_of(&model, &vertex);
        let x = setback * cos30;
        assert!(
            (p.x.abs() - x).abs() < 1e-9 && (p.z - x * k).abs() < 1e-9,
            "{p:?}"
        );
    }
    // Its normal faces the ball's centre, as the faces' do.
    let (p, n) = face_normal(&model, &built.shape, T).unwrap();
    let centre = Point::new(0.0, p.y, r / 60.0_f64.to_radians().sin());
    assert!((n - (centre - p) / r).magnitude() < 1e-9);
}

/// Faces that reach the corner over different stretches: the round spans
/// where both do, and the longer face's line of contact is split there.
#[test]
fn the_round_spans_where_both_faces_reach() {
    let mut model = Model::new();
    let r = 2.0;
    let (a, b) = valley(&mut model, 0.5, (0.0, 8.0), (2.0, 6.0));
    let shell = fillet_faces(&mut model, &a, &b, r, true, T).unwrap().shape;
    usable(&model, &shell);
    let (round, _) = the_round(&model, &shell);
    let round_area = area(&model, &round);
    assert!(
        (round_area - r * (PI / 3.0) * 4.0).abs() < 1e-9,
        "{round_area}"
    );
    for vertex in explore_unique(&model, &round, ShapeType::Vertex).unwrap() {
        let y = point_of(&model, &vertex).y;
        assert!((y - 2.0).abs() < 1e-9 || (y - 6.0).abs() < 1e-9, "{y}");
    }
    assert_eq!(
        edge_use(&model, &shell).0,
        2,
        "both lines of contact shared"
    );
    tangent_to_both(&model, &shell, &round, r);
}

fn refusal(result: ogeom::core::OgeomResult<ogeom::algo::Built>) -> String {
    match result {
        Ok(_) => panic!("expected a refusal"),
        Err(e) => e.to_string(),
    }
}

#[test]
fn sheet_rounds_refuse_by_name() {
    let mut model = Model::new();
    let sheet = open_l(&mut model, 5.0);
    let corner = vertical_edge_at(&model, &sheet, 0.0, 0.0);
    let free = vertical_edge_at(&model, &sheet, 10.0, 0.0);
    let face = faces(&model, &sheet)[0].clone();

    let said = refusal(fillet_sheet_edges(&mut model, &sheet, &[free], 1.0, T));
    assert!(said.contains("bounds one face"), "{said}");
    let said = refusal(fillet_sheet_edges(
        &mut model,
        &sheet,
        std::slice::from_ref(&corner),
        12.0,
        T,
    ));
    assert!(said.contains("sets back"), "{said}");
    let said = refusal(fillet_sheet_edges(
        &mut model,
        &sheet,
        std::slice::from_ref(&corner),
        0.0,
        T,
    ));
    assert!(said.contains("rounds nothing"), "{said}");
    let said = refusal(fillet_sheet_edges(
        &mut model,
        &face,
        std::slice::from_ref(&corner),
        1.0,
        T,
    ));
    assert!(said.contains("rounded on a shell"), "{said}");
    let said = refusal(fillet_sheet_edges(
        &mut model,
        &sheet,
        &[corner.clone(), corner.clone()],
        1.0,
        T,
    ));
    assert!(said.contains("not an edge of the sheet"), "{said}");

    // One flat turned over: the two faces disagree about the corner.
    let flats = faces(&model, &sheet);
    let turned = model
        .add_shell(&[flats[0].clone(), flats[1].reversed()])
        .unwrap();
    let said = refusal(fillet_sheet_edges(
        &mut model,
        &turned,
        std::slice::from_ref(&corner),
        1.0,
        T,
    ));
    assert!(said.contains("oriented inconsistently"), "{said}");

    // Two flats in one plane: no corner.
    let flat = extruded(
        &mut model,
        &[(10.0, 0.0), (0.0, 0.0), (-10.0, 0.0)],
        Vector::new(0.0, 0.0, 5.0),
    );
    let seam = vertical_edge_at(&model, &flat, 0.0, 0.0);
    let said = refusal(fillet_sheet_edges(&mut model, &flat, &[seam], 1.0, T));
    assert!(said.contains("continue each other"), "{said}");

    // Swept on a slant: the ends leave the edge at an angle.
    let leaning = extruded(
        &mut model,
        &[(10.0, 0.0), (0.0, 0.0), (0.0, 10.0)],
        Vector::new(1.0, 1.0, 5.0),
    );
    let slanted = edges_where(&model, &leaning, |p| (p.x - p.y).abs() < 1e-9 && p.x < 1.5)
        .into_iter()
        .next()
        .unwrap();
    let said = refusal(fillet_sheet_edges(&mut model, &leaning, &[slanted], 1.0, T));
    assert!(said.contains("at an angle"), "{said}");
}

/// A third face across the L's foot meets the corner edge's lower end:
/// the round would run into it.
#[test]
fn a_round_into_another_face_is_refused() {
    let mut model = Model::new();
    let sheet = open_l(&mut model, 5.0);
    let corner = vertical_edge_at(&model, &sheet, 0.0, 0.0);
    let feet = edges_where(&model, &sheet, |p| p.z.abs() < 1e-12);
    assert_eq!(feet.len(), 2);
    let ends: Vec<(Shape, Shape)> = feet
        .iter()
        .map(|e| edge_vertices(&model, e).unwrap().unwrap())
        .collect();
    let outer: Vec<Shape> = ends
        .iter()
        .flat_map(|(a, b)| [a.clone(), b.clone()])
        .filter(|v| point_of(&model, v).distance(Point::ORIGIN) > 1.0)
        .collect();
    let (p, q) = (point_of(&model, &outer[0]), point_of(&model, &outer[1]));
    let line = Curve::Line(LineCurve::segment(p, q, T).unwrap());
    let hypotenuse = make_edge_between(
        &mut model,
        line.clone(),
        line.domain(),
        &outer[0],
        &outer[1],
        T,
    )
    .unwrap()
    .shape;
    let mut ring = feet.clone();
    ring.push(hypotenuse);
    let ordered = ogeom::algo::order_edges(&model, &ring, T).unwrap();
    let floor = Plane::through(Point::ORIGIN, Direction::Z);
    let cap = make_face_with_pcurves(
        &mut model,
        PlaneSurface::over(floor, (-50.0, 50.0), (-50.0, 50.0))
            .unwrap()
            .into(),
        &[ordered],
        T,
    )
    .unwrap()
    .shape;
    let mut all = faces(&model, &sheet);
    all.push(cap);
    let capped = model.add_shell(&all).unwrap();
    let said = refusal(fillet_sheet_edges(&mut model, &capped, &[corner], 1.0, T));
    assert!(said.contains("run into another face"), "{said}");
}

/// A floor meeting a cylindrical wall: refused by name, not rounded wrong.
#[test]
fn a_curved_face_is_refused_by_name() {
    let mut model = Model::new();
    let far = make_vertex(&mut model, Point::new(10.0, 0.0, 0.0)).shape;
    let corner = make_vertex(&mut model, Point::ORIGIN).shape;
    let top = make_vertex(&mut model, Point::new(-5.0, 5.0, 0.0)).shape;
    let line =
        Curve::Line(LineCurve::segment(Point::new(10.0, 0.0, 0.0), Point::ORIGIN, T).unwrap());
    let floor = make_edge_between(&mut model, line.clone(), line.domain(), &far, &corner, T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(-5.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let circle = Curve::Circle(CircleCurve::new(Circle::new(frame, 5.0, T).unwrap()));
    let wall = make_edge_between(&mut model, circle, (0.0, PI / 2.0), &corner, &top, T)
        .unwrap()
        .shape;
    let wire = make_wire(&mut model, &[floor, wall], T).unwrap().shape;
    let sheet = make_prism(&mut model, &wire, Vector::new(0.0, 0.0, 5.0), T)
        .unwrap()
        .shape;
    let edge = vertical_edge_at(&model, &sheet, 0.0, 0.0);
    let said = refusal(fillet_sheet_edges(&mut model, &sheet, &[edge], 1.0, T));
    assert!(said.contains("only planar faces"), "{said}");
}

#[test]
fn face_rounds_refuse_by_name() {
    let mut model = Model::new();
    let (a, b) = valley(&mut model, 0.5, (0.0, 8.0), (0.0, 8.0));
    let said = refusal(fillet_faces(&mut model, &a, &a, 1.0, true, T));
    assert!(said.contains("against itself"), "{said}");
    let said = refusal(fillet_faces(&mut model, &a, &b, -1.0, true, T));
    assert!(said.contains("rounds nothing"), "{said}");
    let shell = model.add_shell(std::slice::from_ref(&a)).unwrap();
    let said = refusal(fillet_faces(&mut model, &shell, &b, 1.0, true, T));
    assert!(said.contains("between two faces"), "{said}");

    // A small ball touches the planes near the corner, short of the faces.
    let (far_a, far_b) = valley(&mut model, 6.0, (0.0, 8.0), (0.0, 8.0));
    let said = refusal(fillet_faces(&mut model, &far_a, &far_b, 0.5, true, T));
    assert!(said.contains("wholly beyond"), "{said}");
    // A huge one touches them past their far ends.
    let said = refusal(fillet_faces(&mut model, &a, &b, 100.0, false, T));
    assert!(said.contains("wholly between"), "{said}");
    // Faces that never share a stretch of the corner.
    let (c, d) = valley(&mut model, 0.5, (0.0, 3.0), (5.0, 8.0));
    let said = refusal(fillet_faces(&mut model, &c, &d, 1.0, true, T));
    assert!(said.contains("do not overlap"), "{said}");

    // Parallel planes.
    let copy = ogeom::algo::transformed(
        &mut model,
        &a,
        ogeom::math::Transform::translation(Vector::new(0.0, 0.0, 3.0)),
    )
    .unwrap()
    .shape;
    let said = refusal(fillet_faces(&mut model, &a, &copy, 1.0, true, T));
    assert!(said.contains("parallel"), "{said}");

    // A curved face.
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 3.0, 5.0, T)
        .unwrap()
        .shape;
    let side = faces(&model, &drum)
        .into_iter()
        .find(|f| matches!(surface_of(&model, f), SurfaceGeometry::Cylinder(_)))
        .unwrap();
    let said = refusal(fillet_faces(&mut model, &side, &a, 1.0, true, T));
    assert!(said.contains("only planar faces"), "{said}");
}

/// A floor whose far side notches back toward the corner, to within 1.5
/// of it, beside a plain wall: a round of radius 2 would take a strip 2
/// wide and cross the notch, so it is refused.
#[test]
fn a_notch_in_the_strip_is_refused() {
    let mut model = Model::new();
    let at = |x: f64, y: f64, z: f64| Point::new(x, y, z);
    let floor_points = [
        at(0.0, 0.0, 0.0),
        at(10.0, 0.0, 0.0),
        at(10.0, 0.0, 5.0),
        at(6.0, 0.0, 5.0),
        at(1.5, 0.0, 2.5),
        at(3.0, 0.0, 5.0),
        at(0.0, 0.0, 5.0),
    ];
    let wall_points = [
        at(0.0, 0.0, 5.0),
        at(0.0, 10.0, 5.0),
        at(0.0, 10.0, 0.0),
        at(0.0, 0.0, 0.0),
    ];
    let floor_vertices: Vec<Shape> = floor_points
        .iter()
        .map(|p| make_vertex(&mut model, *p).shape)
        .collect();
    let mut wall_vertices = vec![floor_vertices[6].clone()];
    for p in &wall_points[1..3] {
        wall_vertices.push(make_vertex(&mut model, *p).shape);
    }
    wall_vertices.push(floor_vertices[0].clone());
    let segment = |model: &mut Model, (a, p): (&Shape, Point), (b, q): (&Shape, Point)| {
        let line = Curve::Line(LineCurve::segment(p, q, T).unwrap());
        let range = line.domain();
        make_edge_between(model, line, range, a, b, T)
            .unwrap()
            .shape
    };
    let mut floor_edges = Vec::new();
    for i in 0..floor_points.len() {
        let j = (i + 1) % floor_points.len();
        floor_edges.push(segment(
            &mut model,
            (&floor_vertices[i], floor_points[i]),
            (&floor_vertices[j], floor_points[j]),
        ));
    }
    let corner = floor_edges[6].clone();
    let mut wall_edges = Vec::new();
    for i in 0..3 {
        wall_edges.push(segment(
            &mut model,
            (&wall_vertices[i], wall_points[i]),
            (&wall_vertices[i + 1], wall_points[i + 1]),
        ));
    }
    wall_edges.push(corner.reversed());
    let plane = |normal: Direction| {
        PlaneSurface::over(
            Plane::through(Point::ORIGIN, normal),
            (-50.0, 50.0),
            (-50.0, 50.0),
        )
        .unwrap()
        .into()
    };
    let floor = make_face_with_pcurves(&mut model, plane(Direction::Y), &[floor_edges], T)
        .unwrap()
        .shape;
    let wall = make_face_with_pcurves(&mut model, plane(Direction::X), &[wall_edges], T)
        .unwrap()
        .shape;
    let sheet = model.add_shell(&[floor, wall]).unwrap();
    // A smaller ball clears the notch.
    assert!(fillet_sheet_edges(&mut model, &sheet, std::slice::from_ref(&corner), 1.0, T).is_ok());
    let said = refusal(fillet_sheet_edges(&mut model, &sheet, &[corner], 2.0, T));
    assert!(said.contains("within the strip"), "{said}");
}
