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
fn a_placed_open_l_rounds_where_it_stands() {
    let mut model = Model::new();
    let (r, h) = (2.0, 5.0);
    let sheet = open_l(&mut model, h);
    let moved = model.placed(
        &sheet,
        ogeom::math::Transform::translation(Vector::new(5.0, 0.0, 0.0)),
    );
    let corner = vertical_edge_at(&model, &moved, 5.0, 0.0);
    let built =
        fillet_sheet_edges(&mut model, &moved, std::slice::from_ref(&corner), r, T).unwrap();
    let rounded = built.shape.clone();
    assert_eq!(faces(&model, &rounded).len(), 3, "two flats and the round");
    usable(&model, &rounded);
    assert_eq!(edge_use(&model, &rounded), (2, 8));
    let (round, cylinder) = the_round(&model, &rounded);
    assert!(
        off_axis(&cylinder, Point::new(5.0 + r, r, 0.0)) < 1e-12,
        "the axis passes through (5 + r, r)"
    );
    let round_area = area(&model, &round);
    assert!(
        (round_area - PI / 2.0 * r * h).abs() < 1e-9,
        "a quarter cylinder: {round_area}"
    );
    tangent_to_both(&model, &rounded, &round, r);
    assert!(built.history.is_deleted(&corner));
    assert!(
        built
            .history
            .modified(&moved)
            .iter()
            .any(|s| s.is_same(&rounded))
    );
}

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

/// A floor along x from the z axis to x = 10 and a quarter cylinder wall of
/// radius 5 about the vertical line through (-5, 0), rising from the z axis
/// and curving back over -x, both 5 tall: they meet along the z axis at a
/// right angle, the floor facing +y and the wall facing out of its axis.
fn floor_and_quarter_wall(model: &mut Model) -> Shape {
    floor_and_wall(model, false)
}

/// As [`floor_and_quarter_wall`], the wall's profile given as a B-spline
/// where `spline` holds, so its surface has no closed form to read.
fn floor_and_wall(model: &mut Model, spline: bool) -> Shape {
    let far = make_vertex(model, Point::new(10.0, 0.0, 0.0)).shape;
    let corner = make_vertex(model, Point::ORIGIN).shape;
    let top = make_vertex(model, Point::new(-5.0, 5.0, 0.0)).shape;
    let line =
        Curve::Line(LineCurve::segment(Point::new(10.0, 0.0, 0.0), Point::ORIGIN, T).unwrap());
    let floor = make_edge_between(model, line.clone(), line.domain(), &far, &corner, T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(-5.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let circle = Curve::Circle(CircleCurve::new(Circle::new(frame, 5.0, T).unwrap()));
    let (profile, range) = if spline {
        let curve = Curve::BSpline(circle.to_bspline_over((0.0, PI / 2.0), T).unwrap());
        let domain = curve.domain();
        (curve, domain)
    } else {
        (circle, (0.0, PI / 2.0))
    };
    let wall = make_edge_between(model, profile, range, &corner, &top, T)
        .unwrap()
        .shape;
    let wire = make_wire(model, &[floor, wall], T).unwrap().shape;
    make_prism(model, &wire, Vector::new(0.0, 0.0, 5.0), T)
        .unwrap()
        .shape
}

/// The round's faces against its neighbours: tangent within 1e-5 rad
/// along each shared edge, the edge on both surfaces.
fn tangent_within(model: &Model, shape: &Shape, round: &Shape, angle: f64, gap: f64) {
    let contacts = analyse_blend(model, shape, round, 15, T).unwrap();
    assert_eq!(contacts.len(), 2, "the round meets two faces: {contacts:?}");
    for contact in &contacts {
        assert!(
            contact.tangency_error < angle,
            "tangent within {angle}: {}",
            contact.tangency_error
        );
        assert!(
            contact.gap < gap,
            "the shared edge lies on both: {}",
            contact.gap
        );
    }
}

/// The face of `shape` whose surface `pick` accepts.
fn face_where(model: &Model, shape: &Shape, pick: impl Fn(&SurfaceGeometry) -> bool) -> Shape {
    let found: Vec<Shape> = faces(model, shape)
        .into_iter()
        .filter(|f| pick(&surface_of(model, f)))
        .collect();
    assert_eq!(found.len(), 1, "one such face");
    found[0].clone()
}

/// A plane and a cylinder meeting along a line: the floor and the quarter
/// wall rounded at radius 1 along the z axis. The ball sits 1 above the
/// floor and 5 + 1 from the wall's axis, so its centre is at
/// x = sqrt(6² - 1) - 5; the round is the cylinder of radius 1 about the
/// vertical line through it, touching the floor below it and the wall where
/// the line from the wall's axis to the centre crosses the wall.
#[test]
fn a_floor_and_a_cylindrical_wall_round_to_a_cylinder() {
    let mut model = Model::new();
    let (r, h) = (1.0, 5.0);
    let sheet = floor_and_quarter_wall(&mut model);
    let corner = vertical_edge_at(&model, &sheet, 0.0, 0.0);
    let built =
        fillet_sheet_edges(&mut model, &sheet, std::slice::from_ref(&corner), r, T).unwrap();
    let rounded = built.shape.clone();
    assert_eq!(model.kind_of(&rounded).unwrap(), ShapeType::Shell);
    assert_eq!(faces(&model, &rounded).len(), 3, "floor, wall and round");
    usable(&model, &rounded);
    assert_eq!(edge_use(&model, &rounded), (2, 8));

    let xc = 35.0_f64.sqrt() - 5.0;
    let round = face_where(
        &model,
        &rounded,
        |s| matches!(s, SurfaceGeometry::Cylinder(c) if (c.cylinder().radius() - r).abs() < 1e-12),
    );
    let SurfaceGeometry::Cylinder(c) = surface_of(&model, &round) else {
        unreachable!()
    };
    let cylinder = c.cylinder();
    assert!((cylinder.radius() - r).abs() < 1e-12, "the radius asked");
    assert!(
        cylinder
            .axis()
            .direction
            .vector()
            .cross(Vector::new(0.0, 0.0, 1.0))
            .magnitude()
            < 1e-12
    );
    assert!(
        off_axis(&cylinder, Point::new(xc, r, 0.0)) < 1e-12,
        "the axis at the closed form"
    );
    tangent_within(&model, &rounded, &round, 1e-5, 1e-9);
    // The arc sweeps from straight down to the direction toward the wall's
    // axis, which stands `theta` below -x.
    let theta = 1.0_f64.atan2(35.0_f64.sqrt());
    let sweep = PI / 2.0 - theta;
    let round_area = area(&model, &round);
    assert!(
        (round_area - r * sweep * h).abs() < 1e-9,
        "{round_area} against {}",
        r * sweep * h
    );
    let floor = face_where(&model, &rounded, |s| matches!(s, SurfaceGeometry::Plane(_)));
    assert!((area(&model, &floor) - (10.0 - xc) * h).abs() < 1e-9);
    let wall = face_where(
        &model,
        &rounded,
        |s| matches!(s, SurfaceGeometry::Cylinder(c) if (c.cylinder().radius() - 5.0).abs() < 1e-12),
    );
    assert!((area(&model, &wall) - 5.0 * (PI / 2.0 - theta) * h).abs() < 1e-9);
    // The round faces its axis, as the floor and wall face the corner.
    let (p, n) = face_normal(&model, &round, T).unwrap();
    let to_axis = Vector::new(xc - p.x, r - p.y, 0.0) / r;
    assert!((n - to_axis).magnitude() < 1e-9, "{n:?}");
    assert!(built.history.is_deleted(&corner));
    assert!(
        built
            .history
            .generated(&corner)
            .iter()
            .any(|f| f.node() == round.node())
    );
}

/// The same floor and wall with the wall's profile a B-spline: no closed
/// form is read from the wall, so the ball is marched, and the round lands
/// on the closed form's within the fit's reach.
#[test]
fn a_spline_wall_rounds_by_marching_onto_the_closed_form() {
    let mut model = Model::new();
    let r = 1.0;
    let sheet = floor_and_wall(&mut model, true);
    let wall = faces(&model, &sheet)
        .into_iter()
        .find(|f| !matches!(surface_of(&model, f), SurfaceGeometry::Plane(_)))
        .unwrap();
    assert!(
        !matches!(surface_of(&model, &wall), SurfaceGeometry::Cylinder(_)),
        "the wall is not read as a cylinder"
    );
    let corner = vertical_edge_at(&model, &sheet, 0.0, 0.0);
    let rounded = fillet_sheet_edges(&mut model, &sheet, &[corner], r, T)
        .unwrap()
        .shape;
    assert_eq!(faces(&model, &rounded).len(), 3);
    usable(&model, &rounded);
    assert_eq!(
        edge_use(&model, &rounded).0,
        2,
        "two lines of contact shared"
    );
    let before: Vec<Shape> = faces(&model, &sheet);
    let round = faces(&model, &rounded)
        .into_iter()
        .find(|f| {
            !before.iter().any(|b| b.node() == f.node())
                && explore_unique(&model, f, ShapeType::Edge).unwrap().len() == 4
                && analyse_blend(&model, &rounded, f, 3, T).unwrap().len() == 2
        })
        .expect("the round");
    // The band is fitted through the ball's arcs to two tenths of a
    // micron, so its normal holds to the faces' to the fit's angle: inside
    // a tenth of a degree.
    tangent_within(&model, &rounded, &round, 0.1_f64.to_radians(), 1e-3);
    // Every corner of the round sits where the closed form's ball touches.
    let xc = 35.0_f64.sqrt() - 5.0;
    let wall_contact =
        Point::new(-5.0, 0.0, 0.0) + Vector::new(xc + 5.0, r, 0.0) * (5.0 / (5.0 + r));
    for vertex in explore_unique(&model, &round, ShapeType::Vertex).unwrap() {
        let p = point_of(&model, &vertex);
        let on_floor = Point::new(xc, 0.0, p.z).distance(p);
        let on_wall = Point::new(wall_contact.x, wall_contact.y, p.z).distance(p);
        assert!(on_floor.min(on_wall) < 1e-3, "a contact at {p:?}");
        assert!(
            p.z.abs() < 1e-9 || (p.z - 5.0).abs() < 1e-9,
            "at an end: {p:?}"
        );
    }
    // The area of the closed form's round, within the fit's reach.
    let sweep = PI / 2.0 - 1.0_f64.atan2(35.0_f64.sqrt());
    let round_area = area(&model, &round);
    assert!(
        (round_area - r * sweep * 5.0).abs() < 1e-3,
        "{round_area} against {}",
        r * sweep * 5.0
    );
}

/// A tube of radius 3 and height 6 standing in a round hole in a square
/// floor, as a sheet of two faces sharing the hole's circle: the floor
/// faces up and the tube out of its axis.
fn tube_on_floor(model: &mut Model) -> (Shape, Shape) {
    tube_on_floor_as(model, false)
}

/// As [`tube_on_floor`], the hole's circle given as a B-spline where
/// `spline` holds.
fn tube_on_floor_as(model: &mut Model, spline: bool) -> (Shape, Shape) {
    tube_on_floor_of(model, spline, 10.0)
}

/// As [`tube_on_floor_as`], the floor's square reaching `half` from the
/// axis.
fn tube_on_floor_of(model: &mut Model, spline: bool, half: f64) -> (Shape, Shape) {
    let circle = Curve::Circle(CircleCurve::new(Circle::new(Frame::WORLD, 3.0, T).unwrap()));
    let (profile, range) = if spline {
        let curve = Curve::BSpline(circle.to_bspline_over((0.0, 2.0 * PI), T).unwrap());
        let domain = curve.domain();
        (curve, domain)
    } else {
        (circle, (0.0, 2.0 * PI))
    };
    let ring = ogeom::algo::make_edge(model, profile, range, T)
        .unwrap()
        .shape;
    let tube = make_prism(model, &ring, Vector::new(0.0, 0.0, 6.0), T)
        .unwrap()
        .shape;
    let outer = make_polygon(
        model,
        &[
            Point::new(-half, -half, 0.0),
            Point::new(half, -half, 0.0),
            Point::new(half, half, 0.0),
            Point::new(-half, half, 0.0),
        ],
        true,
        T,
    )
    .unwrap()
    .shape;
    let outer_edges = model.children_of(&outer).unwrap();
    let floor = make_face_with_pcurves(
        model,
        PlaneSurface::new(Plane::through(Point::ORIGIN, Direction::Z)).into(),
        &[outer_edges, vec![ring.reversed()]],
        T,
    )
    .unwrap()
    .shape;
    assert!(face_normal(model, &floor, T).unwrap().1.z > 0.5);
    let sheet = model.add_shell(&[floor, tube]).unwrap();
    (sheet, ring)
}

/// A plane meeting a cylinder square to it along a circle: the tube's
/// foot rounded at radius 1. The ball rolls round the outside of the tube
/// on the floor, its centre on the circle of radius 3 + 1 at height 1; the
/// round is the torus about the tube's axis with that circle for its core
/// and the radius asked for its tube, touching the floor on the circle of
/// radius 4 and the tube on its circle at height 1.
#[test]
fn a_tube_on_a_floor_rounds_to_a_torus() {
    let mut model = Model::new();
    let r = 1.0;
    let (sheet, ring) = tube_on_floor(&mut model);
    usable(&model, &sheet);
    let built = fillet_sheet_edges(&mut model, &sheet, std::slice::from_ref(&ring), r, T).unwrap();
    let rounded = built.shape.clone();
    assert_eq!(faces(&model, &rounded).len(), 3);
    usable(&model, &rounded);
    let round = face_where(&model, &rounded, |s| matches!(s, SurfaceGeometry::Torus(_)));
    let SurfaceGeometry::Torus(t) = surface_of(&model, &round) else {
        unreachable!()
    };
    let torus = t.torus();
    assert!((torus.minor_radius() - r).abs() < 1e-12, "the radius asked");
    assert!((torus.major_radius() - 4.0).abs() < 1e-12);
    assert!(torus.centre().distance(Point::new(0.0, 0.0, r)) < 1e-12);
    assert!(
        torus
            .axis()
            .direction
            .vector()
            .cross(Vector::new(0.0, 0.0, 1.0))
            .magnitude()
            < 1e-12
    );
    tangent_within(&model, &rounded, &round, 1e-5, 1e-9);
    // A quarter of the tube's turn about its core, from below the core to
    // its inside: 2π r ∫ (4 + r cos v) dv over that quarter.
    let expected = 2.0 * PI * r * (4.0 * PI / 2.0 - r);
    let round_area = area(&model, &round);
    assert!(
        (round_area - expected).abs() < 1e-9,
        "{round_area} against {expected}"
    );
    let floor = face_where(&model, &rounded, |s| matches!(s, SurfaceGeometry::Plane(_)));
    assert!((area(&model, &floor) - (400.0 - PI * 16.0)).abs() < 1e-9);
    let tube = face_where(&model, &rounded, |s| {
        matches!(s, SurfaceGeometry::Cylinder(_))
    });
    assert!((area(&model, &tube) - 2.0 * PI * 3.0 * 5.0).abs() < 1e-9);
    // The round faces its core, as the floor faces up toward it and the
    // tube out toward it.
    let (p, n) = face_normal(&model, &round, T).unwrap();
    let radial = Vector::new(p.x, p.y, 0.0).normalized(T).unwrap();
    let core = Point::new(0.0, 0.0, r) + radial * 4.0;
    assert!((n - (core - p) / r).magnitude() < 1e-9, "{n:?}");
    assert!(built.history.is_deleted(&ring));
}

/// The floor z = 0 over the rectangle `x` by `y`, facing up.
fn floor_rectangle(model: &mut Model, x: (f64, f64), y: (f64, f64)) -> Shape {
    let corners = [
        Point::new(x.0, y.0, 0.0),
        Point::new(x.1, y.0, 0.0),
        Point::new(x.1, y.1, 0.0),
        Point::new(x.0, y.1, 0.0),
    ];
    let wire = make_polygon(model, &corners, true, T).unwrap().shape;
    let edges = model.children_of(&wire).unwrap();
    let floor = make_face_with_pcurves(
        model,
        PlaneSurface::new(Plane::through(Point::ORIGIN, Direction::Z)).into(),
        &[edges],
        T,
    )
    .unwrap()
    .shape;
    assert!(face_normal(model, &floor, T).unwrap().1.z > 0.5);
    floor
}

/// The quarter of the cylinder of radius 3 about the line x = 0, z = 5
/// running from +x round to straight down, over `y`, facing out of its
/// axis.
fn drum_quarter(model: &mut Model, y: (f64, f64)) -> Shape {
    let frame = Frame::new(Point::new(0.0, y.0, 5.0), Direction::Y, Direction::X, T).unwrap();
    let circle = Curve::Circle(CircleCurve::new(Circle::new(frame, 3.0, T).unwrap()));
    // The frame's y is z × x = -z here, so angles run from +x down.
    let arc = ogeom::algo::make_edge(model, circle, (0.0, PI / 2.0), T)
        .unwrap()
        .shape;
    let face = make_prism(model, &arc, Vector::new(0.0, y.1 - y.0, 0.0), T)
        .unwrap()
        .shape;
    let (p, n) = face_normal(model, &face, T).unwrap();
    let out = Vector::new(p.x, 0.0, p.z - 5.0);
    if n.dot(out) > 0.0 {
        face
    } else {
        face.reversed()
    }
}

/// Half a tube on a floor along an arc: the floor z = 0 over y ≥ 0 out to
/// 10, notched by the half disc of radius 3 about the z axis, and the half
/// tube of radius 3 standing 6 tall on the notch's arc. Rounded at radius
/// 1 the round is half the torus of the whole tube's foot, ending in the
/// half-plane y = 0 at both ends, where the floor's and the tube's
/// boundaries leave the arc.
#[test]
fn a_half_tube_on_a_floor_rounds_to_half_a_torus() {
    let mut model = Model::new();
    let r = 1.0;
    let at = |x: f64, y: f64| Point::new(x, y, 0.0);
    let points = [
        at(10.0, 0.0),
        at(10.0, 10.0),
        at(-10.0, 10.0),
        at(-10.0, 0.0),
        at(-3.0, 0.0),
        at(3.0, 0.0),
    ];
    let v: Vec<Shape> = points
        .iter()
        .map(|p| make_vertex(&mut model, *p).shape)
        .collect();
    let mut edges = Vec::new();
    for i in 0..4 {
        let line = Curve::Line(LineCurve::segment(points[i], points[i + 1], T).unwrap());
        let range = line.domain();
        edges.push(
            make_edge_between(&mut model, line, range, &v[i], &v[i + 1], T)
                .unwrap()
                .shape,
        );
    }
    let circle = Curve::Circle(CircleCurve::new(Circle::new(Frame::WORLD, 3.0, T).unwrap()));
    let arc = make_edge_between(&mut model, circle, (0.0, PI), &v[5], &v[4], T)
        .unwrap()
        .shape;
    edges.push(arc.reversed());
    let line = Curve::Line(LineCurve::segment(points[5], points[0], T).unwrap());
    let range = line.domain();
    edges.push(
        make_edge_between(&mut model, line, range, &v[5], &v[0], T)
            .unwrap()
            .shape,
    );
    let floor = make_face_with_pcurves(
        &mut model,
        PlaneSurface::new(Plane::through(Point::ORIGIN, Direction::Z)).into(),
        &[edges],
        T,
    )
    .unwrap()
    .shape;
    let floor = if face_normal(&model, &floor, T).unwrap().1.z > 0.0 {
        floor
    } else {
        floor.reversed()
    };
    let tube = make_prism(&mut model, &arc, Vector::new(0.0, 0.0, 6.0), T)
        .unwrap()
        .shape;
    let (p, n) = face_normal(&model, &tube, T).unwrap();
    let tube = if n.dot(Vector::new(p.x, p.y, 0.0)) > 0.0 {
        tube
    } else {
        tube.reversed()
    };
    let sheet = model.add_shell(&[floor, tube]).unwrap();
    usable(&model, &sheet);
    let rounded = fillet_sheet_edges(&mut model, &sheet, &[arc], r, T)
        .unwrap()
        .shape;
    assert_eq!(faces(&model, &rounded).len(), 3);
    usable(&model, &rounded);
    let round = face_where(&model, &rounded, |s| matches!(s, SurfaceGeometry::Torus(_)));
    let SurfaceGeometry::Torus(t) = surface_of(&model, &round) else {
        unreachable!()
    };
    assert!(
        (t.torus().minor_radius() - r).abs() < 1e-12,
        "the radius asked"
    );
    assert!((t.torus().major_radius() - 4.0).abs() < 1e-12);
    tangent_within(&model, &rounded, &round, 1e-5, 1e-9);
    let expected = PI * r * (4.0 * PI / 2.0 - r);
    let round_area = area(&model, &round);
    assert!(
        (round_area - expected).abs() < 1e-9,
        "{round_area} against {expected}"
    );
    // The floor loses the half annulus between radius 3 and 4.
    let floor = face_where(&model, &rounded, |s| matches!(s, SurfaceGeometry::Plane(_)));
    let floor_area = 200.0 - PI * 16.0 / 2.0;
    assert!((area(&model, &floor) - floor_area).abs() < 1e-9);
    for vertex in explore_unique(&model, &round, ShapeType::Vertex).unwrap() {
        assert!(point_of(&model, &vertex).y.abs() < 1e-12, "an end on y = 0");
    }
}

/// A plane and a separate cylinder parallel to it: the floor z = 0 and a
/// quarter of the cylinder of radius 3 about the line x = 0, z = 5 along
/// y, both 8 long, facing the gap between them. A ball of radius 2 sits 2
/// above the floor and 5 from the cylinder's axis, at x = 4; it touches
/// the floor at x = 4 and the cylinder at (2.4, 3.2). Trimmed, the three
/// faces sew into one shell.
#[test]
fn a_plane_and_a_separate_cylinder_trim_into_one_shell() {
    let mut model = Model::new();
    let r = 2.0;
    let floor = floor_rectangle(&mut model, (-10.0, 10.0), (0.0, 8.0));
    let drum = drum_quarter(&mut model, (0.0, 8.0));
    let built = fillet_faces(&mut model, &floor, &drum, r, true, T).unwrap();
    let shell = built.shape.clone();
    assert_eq!(model.kind_of(&shell).unwrap(), ShapeType::Shell);
    assert_eq!(faces(&model, &shell).len(), 3);
    usable(&model, &shell);
    let three = faces(&model, &shell);
    let sewn = ogeom::algo::sew(&mut model, &three, T).unwrap();
    assert_eq!(sewn.shells.len(), 1, "the three faces sew into one shell");
    assert_eq!(
        edge_use(&model, &shell).0,
        2,
        "both lines of contact shared"
    );

    let round = face_where(
        &model,
        &shell,
        |s| matches!(s, SurfaceGeometry::Cylinder(c) if (c.cylinder().radius() - r).abs() < 1e-12),
    );
    let SurfaceGeometry::Cylinder(c) = surface_of(&model, &round) else {
        unreachable!()
    };
    assert!(off_axis(&c.cylinder(), Point::new(4.0, 0.0, 2.0)) < 1e-12);
    assert!(off_axis(&c.cylinder(), Point::new(4.0, 5.0, 2.0)) < 1e-12);
    tangent_within(&model, &shell, &round, 1e-5, 1e-9);
    // From straight down to toward the drum's axis: the angle whose
    // cosine is the two directions' dot, (0, -1) · (-0.8, 0.6) = -0.6.
    let sweep = (-0.6_f64).acos();
    assert!((area(&model, &round) - r * sweep * 8.0).abs() < 1e-9);
    let kept_floor = face_where(&model, &shell, |s| matches!(s, SurfaceGeometry::Plane(_)));
    assert!(
        (area(&model, &kept_floor) - 6.0 * 8.0).abs() < 1e-9,
        "the floor keeps x ≥ 4"
    );
    let kept_drum = face_where(
        &model,
        &shell,
        |s| matches!(s, SurfaceGeometry::Cylinder(c) if (c.cylinder().radius() - 3.0).abs() < 1e-12),
    );
    let kept = 0.6_f64.atan2(0.8);
    assert!(
        (area(&model, &kept_drum) - 3.0 * kept * 8.0).abs() < 1e-9,
        "the drum keeps the part above its contact"
    );
    // The round faces the ball's centre, as the floor and drum face it.
    let (p, n) = face_normal(&model, &round, T).unwrap();
    assert!((n - (Point::new(4.0, p.y, 2.0) - p) / r).magnitude() < 1e-9);
    assert_eq!(built.history.modified(&floor).len(), 1);
}

/// A plane and a separate cylinder square to it: the floor z = 0 and a
/// tube of radius 3 standing clear of it from z = 1 to z = 6, both facing
/// the gap. A ball of radius 2 rolls round the tube on the floor, its
/// centre on the circle of radius 5 at height 2: the round is the torus
/// with that core, touching the floor on the circle of radius 5 and the
/// tube on its circle at height 2, closed all the way round.
#[test]
fn a_plane_and_a_separate_tube_trim_to_a_torus() {
    let mut model = Model::new();
    let r = 2.0;
    let floor = {
        let corners = [
            Point::new(-10.0, -10.0, 0.0),
            Point::new(10.0, -10.0, 0.0),
            Point::new(10.0, 10.0, 0.0),
            Point::new(-10.0, 10.0, 0.0),
        ];
        let wire = make_polygon(&mut model, &corners, true, T).unwrap().shape;
        let edges = model.children_of(&wire).unwrap();
        make_face_with_pcurves(
            &mut model,
            PlaneSurface::new(Plane::through(Point::ORIGIN, Direction::Z)).into(),
            &[edges],
            T,
        )
        .unwrap()
        .shape
    };
    let tube = {
        let frame = Frame::new(Point::new(0.0, 0.0, 1.0), Direction::Z, Direction::X, T).unwrap();
        let ring = ogeom::algo::make_edge(
            &mut model,
            Curve::Circle(CircleCurve::new(Circle::new(frame, 3.0, T).unwrap())),
            (0.0, 2.0 * PI),
            T,
        )
        .unwrap()
        .shape;
        make_prism(&mut model, &ring, Vector::new(0.0, 0.0, 5.0), T)
            .unwrap()
            .shape
    };
    let (p, n) = face_normal(&model, &tube, T).unwrap();
    assert!(
        n.dot(Vector::new(p.x, p.y, 0.0)) > 0.0,
        "the tube faces out"
    );
    let built = fillet_faces(&mut model, &floor, &tube, r, true, T).unwrap();
    let shell = built.shape.clone();
    assert_eq!(faces(&model, &shell).len(), 3);
    usable(&model, &shell);
    let three = faces(&model, &shell);
    let sewn = ogeom::algo::sew(&mut model, &three, T).unwrap();
    assert_eq!(sewn.shells.len(), 1, "the three faces sew into one shell");
    let round = face_where(&model, &shell, |s| matches!(s, SurfaceGeometry::Torus(_)));
    let SurfaceGeometry::Torus(t) = surface_of(&model, &round) else {
        unreachable!()
    };
    let torus = t.torus();
    assert!((torus.minor_radius() - r).abs() < 1e-12, "the radius asked");
    assert!((torus.major_radius() - 5.0).abs() < 1e-12);
    assert!(torus.centre().distance(Point::new(0.0, 0.0, r)) < 1e-12);
    tangent_within(&model, &shell, &round, 1e-5, 1e-9);
    let expected = 2.0 * PI * r * (5.0 * PI / 2.0 - r);
    let round_area = area(&model, &round);
    assert!(
        (round_area - expected).abs() < 1e-9,
        "{round_area} against {expected}"
    );
    let kept_floor = face_where(&model, &shell, |s| matches!(s, SurfaceGeometry::Plane(_)));
    assert!((area(&model, &kept_floor) - (400.0 - 25.0 * PI)).abs() < 1e-9);
    let kept_tube = face_where(&model, &shell, |s| {
        matches!(s, SurfaceGeometry::Cylinder(_))
    });
    assert!((area(&model, &kept_tube) - 2.0 * PI * 3.0 * 4.0).abs() < 1e-9);
}

/// A quarter of the same tube over the same floor: the round is the
/// quarter of the torus over the quarter turn the tube spans, ending in
/// the half-planes through the axis at its ends.
#[test]
fn a_plane_and_a_quarter_tube_round_over_its_quarter_turn() {
    let mut model = Model::new();
    let r = 2.0;
    let floor = floor_rectangle(&mut model, (-10.0, 10.0), (-10.0, 10.0));
    let tube = {
        let frame = Frame::new(Point::new(0.0, 0.0, 1.0), Direction::Z, Direction::X, T).unwrap();
        let arc = ogeom::algo::make_edge(
            &mut model,
            Curve::Circle(CircleCurve::new(Circle::new(frame, 3.0, T).unwrap())),
            (0.0, PI / 2.0),
            T,
        )
        .unwrap()
        .shape;
        let face = make_prism(&mut model, &arc, Vector::new(0.0, 0.0, 5.0), T)
            .unwrap()
            .shape;
        let (p, n) = face_normal(&model, &face, T).unwrap();
        if n.dot(Vector::new(p.x, p.y, 0.0)) > 0.0 {
            face
        } else {
            face.reversed()
        }
    };
    let built = fillet_faces(&mut model, &floor, &tube, r, true, T).unwrap();
    let shell = built.shape.clone();
    assert_eq!(faces(&model, &shell).len(), 3);
    usable(&model, &shell);
    let three = faces(&model, &shell);
    let sewn = ogeom::algo::sew(&mut model, &three, T).unwrap();
    assert_eq!(sewn.shells.len(), 1, "the three faces sew into one shell");
    let round = face_where(&model, &shell, |s| matches!(s, SurfaceGeometry::Torus(_)));
    tangent_within(&model, &shell, &round, 1e-5, 1e-9);
    let expected = 2.0 * PI * r * (5.0 * PI / 2.0 - r) / 4.0;
    let round_area = area(&model, &round);
    assert!(
        (round_area - expected).abs() < 1e-9,
        "{round_area} against {expected}"
    );
    // The floor is cut back to its whole line of contact, the circle of
    // radius 5, and shares the quarter of it the round spans; the rest of
    // the circle is free boundary.
    let kept_floor = face_where(&model, &shell, |s| matches!(s, SurfaceGeometry::Plane(_)));
    assert!(
        (area(&model, &kept_floor) - (400.0 - 25.0 * PI)).abs() < 1e-9,
        "{}",
        area(&model, &kept_floor)
    );
    assert_eq!(
        edge_use(&model, &shell).0,
        2,
        "both lines of contact shared"
    );
}

/// What curved faces are still refused, each by name.
#[test]
fn curved_rounds_refuse_by_name() {
    // A closed edge whose seat has no closed form: the tube swept from a
    // B-spline ring would have to be marched round a loop.
    let mut model = Model::new();
    let (splined, ring) = tube_on_floor_as(&mut model, true);
    let said = refusal(fillet_sheet_edges(&mut model, &splined, &[ring], 1.0, T));
    assert!(said.contains("would have to be marched"), "{said}");

    // A floor reaching 4.5 from the tube's axis: a ball of radius 2 would
    // touch it on the circle of radius 5, across its edges.
    let mut model = Model::new();
    let (narrow, ring) = tube_on_floor_of(&mut model, false, 4.5);
    assert!(fillet_sheet_edges(&mut model, &narrow, std::slice::from_ref(&ring), 1.0, T).is_ok());
    let said = refusal(fillet_sheet_edges(&mut model, &narrow, &[ring], 2.0, T));
    assert!(said.contains("within the strip"), "{said}");

    // A cup: a disc floor inside the tube, both facing in. A ball of radius
    // 2 inside the radius 3 has its centre 1 from the axis, so its round
    // would be a torus crossing its own axis.
    let mut model = Model::new();
    let (cup, rim) = {
        let (sheet, ring) = tube_on_floor(&mut model);
        let tube = face_where(&model, &sheet, |s| {
            matches!(s, SurfaceGeometry::Cylinder(_))
        });
        let disc = make_face_with_pcurves(
            &mut model,
            PlaneSurface::new(Plane::through(Point::ORIGIN, Direction::Z)).into(),
            &[vec![ring.clone()]],
            T,
        )
        .unwrap()
        .shape;
        assert!(face_normal(&model, &disc, T).unwrap().1.z > 0.5);
        (model.add_shell(&[disc, tube.reversed()]).unwrap(), ring)
    };
    assert!(fillet_sheet_edges(&mut model, &cup, std::slice::from_ref(&rim), 1.0, T).is_ok());
    let said = refusal(fillet_sheet_edges(&mut model, &cup, &[rim], 2.0, T));
    assert!(said.contains("crossing its own axis"), "{said}");

    // A ball too big for the inside of a cylinder it would roll in: a floor
    // meeting a quarter wall of radius 5 that curves over it.
    let mut model = Model::new();
    let sheet = {
        let far = make_vertex(&mut model, Point::new(10.0, 0.0, 0.0)).shape;
        let corner = make_vertex(&mut model, Point::ORIGIN).shape;
        let top = make_vertex(&mut model, Point::new(5.0, 5.0, 0.0)).shape;
        let line =
            Curve::Line(LineCurve::segment(Point::new(10.0, 0.0, 0.0), Point::ORIGIN, T).unwrap());
        let floor = make_edge_between(&mut model, line.clone(), line.domain(), &far, &corner, T)
            .unwrap()
            .shape;
        let frame = Frame::new(Point::new(5.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
        let circle = Curve::Circle(CircleCurve::new(Circle::new(frame, 5.0, T).unwrap()));
        let wall = make_edge_between(&mut model, circle, (PI / 2.0, PI), &top, &corner, T)
            .unwrap()
            .shape;
        let wire = make_wire(&mut model, &[floor, wall.reversed()], T)
            .unwrap()
            .shape;
        make_prism(&mut model, &wire, Vector::new(0.0, 0.0, 5.0), T)
            .unwrap()
            .shape
    };
    let corner = vertical_edge_at(&model, &sheet, 0.0, 0.0);
    let said = refusal(fillet_sheet_edges(&mut model, &sheet, &[corner], 6.0, T));
    assert!(said.contains("does not fit inside"), "{said}");

    // A tube and a separate plane tilted across it: the ball rolls round
    // the tube, and its line of contact crosses one of the faces more than
    // once, so which stretch of it to round is ambiguous.
    let mut model = Model::new();
    let (a, _) = valley(&mut model, 0.5, (0.0, 8.0), (0.0, 8.0));
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 3.0, 5.0, T)
        .unwrap()
        .shape;
    let side = faces(&model, &drum)
        .into_iter()
        .find(|f| matches!(surface_of(&model, f), SurfaceGeometry::Cylinder(_)))
        .unwrap();
    let said = refusal(fillet_faces(&mut model, &side, &a, 1.0, true, T));
    assert!(said.contains("crosses it more than once"), "{said}");

    // A floor short of where the ball beside the drum touches it (x = 4),
    // on either side.
    let mut model = Model::new();
    let floor = floor_rectangle(&mut model, (5.0, 10.0), (0.0, 8.0));
    let drum = drum_quarter(&mut model, (0.0, 8.0));
    let said = refusal(fillet_faces(&mut model, &floor, &drum, 2.0, true, T));
    assert!(said.contains("does not touch both faces"), "{said}");
    // A floor and a drum side by side along y, never beside each other.
    let floor = floor_rectangle(&mut model, (-10.0, 10.0), (0.0, 8.0));
    let drum = drum_quarter(&mut model, (10.0, 18.0));
    let said = refusal(fillet_faces(&mut model, &floor, &drum, 2.0, false, T));
    assert!(said.contains("do not overlap"), "{said}");

    // A ball so large it touches the floor past the floor's far edge.
    let mut model = Model::new();
    let sheet = floor_and_quarter_wall(&mut model);
    let corner = vertical_edge_at(&model, &sheet, 0.0, 0.0);
    let said = refusal(fillet_sheet_edges(&mut model, &sheet, &[corner], 25.0, T));
    assert!(said.contains("sets back past"), "{said}");

    // A floor whose foot rises from the corner, (0, 0, 0) to (10, 0, 2),
    // beside the quarter wall standing on z = 0: the round ends in the
    // plane z = 0, which the floor's foot leaves.
    let mut model = Model::new();
    let at = |x: f64, y: f64, z: f64| Point::new(x, y, z);
    let points = [
        at(0.0, 0.0, 0.0),
        at(10.0, 0.0, 2.0),
        at(10.0, 0.0, 5.0),
        at(0.0, 0.0, 5.0),
        at(-5.0, 5.0, 0.0),
        at(-5.0, 5.0, 5.0),
    ];
    let v: Vec<Shape> = points
        .iter()
        .map(|p| make_vertex(&mut model, *p).shape)
        .collect();
    let segment = |model: &mut Model, i: usize, j: usize| {
        let line = Curve::Line(LineCurve::segment(points[i], points[j], T).unwrap());
        let range = line.domain();
        make_edge_between(model, line, range, &v[i], &v[j], T)
            .unwrap()
            .shape
    };
    let arc = |model: &mut Model, z: f64, i: usize, j: usize| {
        let frame = Frame::new(Point::new(-5.0, 0.0, z), Direction::Z, Direction::X, T).unwrap();
        let circle = Curve::Circle(CircleCurve::new(Circle::new(frame, 5.0, T).unwrap()));
        make_edge_between(model, circle, (0.0, PI / 2.0), &v[i], &v[j], T)
            .unwrap()
            .shape
    };
    let corner = segment(&mut model, 0, 3);
    let floor_edges = vec![
        segment(&mut model, 0, 1),
        segment(&mut model, 1, 2),
        segment(&mut model, 2, 3),
        corner.reversed(),
    ];
    let wall_edges = vec![
        arc(&mut model, 0.0, 0, 4),
        segment(&mut model, 4, 5),
        arc(&mut model, 5.0, 3, 5).reversed(),
        corner.reversed(),
    ];
    let floor = make_face_with_pcurves(
        &mut model,
        PlaneSurface::new(Plane::through(Point::ORIGIN, Direction::Y)).into(),
        &[floor_edges],
        T,
    )
    .unwrap()
    .shape;
    let wall_frame = Frame::new(Point::new(-5.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let wall = make_face_with_pcurves(
        &mut model,
        ogeom::geom::CylinderSurface::new(
            ogeom::math::Cylinder::new(wall_frame, 5.0, T).unwrap(),
            (-1.0, 6.0),
        )
        .unwrap()
        .into(),
        &[wall_edges],
        T,
    )
    .unwrap()
    .shape;
    let facing = |model: &Model, face: Shape, into: Vector| {
        if face_normal(model, &face, T).unwrap().1.dot(into) > 0.0 {
            face
        } else {
            face.reversed()
        }
    };
    let floor = facing(&model, floor, Vector::new(0.0, 1.0, 0.0));
    let wall = facing(&model, wall, Vector::new(1.0, 0.0, 0.0));
    let sheet = model.add_shell(&[floor, wall]).unwrap();
    let said = refusal(fillet_sheet_edges(&mut model, &sheet, &[corner], 1.0, T));
    assert!(said.contains("outside the round's end section"), "{said}");

    // The spline wall swept on a slant: at the edge's ends the ball's
    // section stands square to the edge, and there it would touch the wall
    // below its foot.
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
    let spline = Curve::BSpline(circle.to_bspline_over((0.0, PI / 2.0), T).unwrap());
    let domain = spline.domain();
    let wall = make_edge_between(&mut model, spline, domain, &corner, &top, T)
        .unwrap()
        .shape;
    let wire = make_wire(&mut model, &[floor, wall], T).unwrap().shape;
    let leaning = make_prism(&mut model, &wire, Vector::new(2.0, 0.0, 5.0), T)
        .unwrap()
        .shape;
    let slanted = edges_where(&model, &leaning, |p| {
        p.y.abs() < 1e-9 && (p.x - 0.4 * p.z).abs() < 1e-9
    })
    .into_iter()
    .next()
    .unwrap();
    let said = refusal(fillet_sheet_edges(&mut model, &leaning, &[slanted], 1.0, T));
    assert!(said.contains("does not seat"), "{said}");
}

/// A face over the whole of a B-spline patch, bounded by its four border
/// iso-curves with their lines in the chart.
fn patch_face(model: &mut Model, patch: ogeom::geom::BSplineSurface) -> Shape {
    use ogeom::geom::{Line2d, Surface as _};
    use ogeom::math::{Axis2, Direction2, Point2};
    let ((u0, u1), (v0, v1)) = patch.domain();
    let corners: Vec<Shape> = [(u0, v0), (u1, v0), (u1, v1), (u0, v1)]
        .iter()
        .map(|&(u, v)| make_vertex(model, patch.point_at(u, v, T).unwrap()).shape)
        .collect();
    let id = model.geometry_mut().add_surface(patch.clone().into());
    let mut edges = Vec::with_capacity(4);
    // Each side: whether it is a u iso-line, where, its corners, and
    // whether the loop runs it backwards.
    for (iso_u, at, from, to, backwards) in [
        (false, v0, 0, 1, false),
        (true, u1, 1, 2, false),
        (false, v1, 3, 2, true),
        (true, u0, 0, 3, true),
    ] {
        let (curve, range, image) = if iso_u {
            let line = Axis2::new(Point2::new(at, 0.0), Direction2::Y);
            (
                patch.iso_u_curve(at, T).unwrap(),
                (v0, v1),
                Line2d::over(line, v0 - 1.0, v1 + 1.0).unwrap(),
            )
        } else {
            let line = Axis2::new(Point2::new(0.0, at), Direction2::X);
            (
                patch.iso_v_curve(at, T).unwrap(),
                (u0, u1),
                Line2d::over(line, u0 - 1.0, u1 + 1.0).unwrap(),
            )
        };
        let edge = make_edge_between(
            model,
            Curve::BSpline(curve),
            range,
            &corners[from],
            &corners[to],
            T,
        )
        .unwrap()
        .shape;
        ogeom::algo::attach_pcurve(
            model,
            &edge,
            image.into(),
            id,
            ogeom::topo::Location::identity(),
            range,
        )
        .unwrap();
        edges.push(if backwards { edge.reversed() } else { edge });
    }
    let wire = make_wire(model, &edges, T).unwrap().shape;
    ogeom::algo::make_face_on(model, id, &[wire], T)
        .unwrap()
        .shape
}

/// A cubic patch rising from below the floor z = 0 to above it, leaning
/// 60 degrees toward +x and bowed both ways, over y from -8 to 8: it
/// crosses the floor along a curve near the y axis. It faces the acute
/// corner it makes with the floor's +x side.
fn leaning_patch(model: &mut Model) -> Shape {
    use ogeom::math::{ControlGrid, KnotVector};
    let (s, c) = 60.0_f64.to_radians().sin_cos();
    let mut points = Vec::with_capacity(16);
    for i in 0..4 {
        let t = 11.0_f64.mul_add(f64::from(i) / 3.0, -3.0);
        for j in 0..4 {
            let y = 16.0_f64.mul_add(f64::from(j) / 3.0, -8.0);
            let bow = 0.4 * (y / 8.0).powi(2) + 0.3 * (t / 8.0).powi(2);
            points.push(Point::new(t.mul_add(c, bow), y, t * s));
        }
    }
    let cubic = || KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3).unwrap();
    let patch = ogeom::geom::BSplineSurface::new(
        cubic(),
        cubic(),
        &ControlGrid::new(points, 4, 4).unwrap(),
        T,
    )
    .unwrap();
    let face = patch_face(model, patch);
    if face_normal(model, &face, T)
        .unwrap()
        .1
        .dot(Vector::new(s, 0.0, -c))
        > 0.0
    {
        face
    } else {
        face.reversed()
    }
}

/// The round's tangency and fit against each face it meets: tangent
/// within a tenth of a degree, its shared edge on both surfaces within the
/// band's fit.
fn marched_tangent(model: &Model, shape: &Shape, round: &Shape) {
    let contacts = analyse_blend(model, shape, round, 15, T).unwrap();
    assert_eq!(contacts.len(), 2, "the round meets two faces: {contacts:?}");
    for contact in &contacts {
        assert!(
            contact.tangency_error < 0.1_f64.to_radians(),
            "tangent within a tenth of a degree: {contact:?}"
        );
        assert!(
            contact.gap < 2e-4,
            "the shared edge lies on both: {contact:?}"
        );
    }
}

/// The face of `shape` holding `p`.
fn face_holding(model: &Model, shape: &Shape, p: Point) -> Option<Shape> {
    faces(model, shape).into_iter().find(|f| {
        ogeom::algo::classify_on_face(model, f, p, Deflection::with_chord(1e-3).unwrap(), T)
            .unwrap()
            == ogeom::algo::Containment::In
    })
}

/// The round's corners standing on the floor z = 0, within the band's fit.
fn ends_on_floor(model: &Model, round: &Shape) -> Vec<Point> {
    explore_unique(model, round, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| point_of(model, v))
        .filter(|p| p.z.abs() < 2e-4)
        .collect()
}

/// A floor and a separate B-spline patch crossing it at a slant: they
/// share no direction or axis, so the ball is marched along where they
/// cross and the round is a band fitted through its arcs. Trimmed, the
/// floor keeps its side away from the corner and the patch its top, and
/// the three faces sew into one shell; the round runs the floor's width,
/// ending at its edges y = ±5 where the patch is 16 wide, and rides both
/// faces within a tenth of a degree.
#[test]
fn a_plane_and_a_leaning_spline_patch_trim_into_one_shell() {
    let mut model = Model::new();
    let r = 1.5;
    let floor = floor_rectangle(&mut model, (-10.0, 10.0), (-5.0, 5.0));
    let patch = leaning_patch(&mut model);
    let built = fillet_faces(&mut model, &floor, &patch, r, true, T).unwrap();
    let shell = built.shape.clone();
    assert_eq!(model.kind_of(&shell).unwrap(), ShapeType::Shell);
    assert_eq!(faces(&model, &shell).len(), 3);
    usable(&model, &shell);
    let three = faces(&model, &shell);
    let sewn = ogeom::algo::sew(&mut model, &three, T).unwrap();
    assert_eq!(sewn.shells.len(), 1, "the three faces sew into one shell");
    assert_eq!(
        edge_use(&model, &shell).0,
        2,
        "both lines of contact shared"
    );

    let kept_floor = face_holding(&model, &shell, Point::new(9.0, 0.0, 0.0)).unwrap();
    assert!(
        face_holding(&model, &shell, Point::new(0.0, 0.0, 0.0)).is_none(),
        "the floor loses the corner"
    );
    let round = faces(&model, &shell)
        .into_iter()
        .find(|f| {
            !f.is_same(&kept_floor) && built.history.modified(&patch).iter().all(|k| !k.is_same(f))
        })
        .unwrap();
    marched_tangent(&model, &shell, &round);
    // The run ends where the floor does: the round's corners on the floor
    // stand on its edges y = ±5, those on the patch in the same sections.
    let on_floor = ends_on_floor(&model, &round);
    assert_eq!(on_floor.len(), 2, "two corners on the floor");
    for p in on_floor {
        assert!(
            (p.y.abs() - 5.0).abs() < 1e-9,
            "a corner off the floor's edge: {p:?}"
        );
    }
    // The patch keeps what stands above its line of contact, which runs
    // about 2.6 up it from the floor (the ball's radius over tan 30°): the
    // patch over y = 0 stands 0.74 over the floor at u = 0.35 and about 6
    // at u = 0.9.
    let surface = surface_of(&model, &patch);
    let at = |u: f64| {
        use ogeom::geom::Surface as _;
        surface.point_at(u, 0.5, T).unwrap()
    };
    assert!(face_holding(&model, &shell, at(0.9)).is_some());
    assert!(face_holding(&model, &shell, at(0.35)).is_none());

    // Untrimmed, the round alone, over the same run.
    let mut model = Model::new();
    let floor = floor_rectangle(&mut model, (-10.0, 10.0), (-5.0, 5.0));
    let patch = leaning_patch(&mut model);
    let round = fillet_faces(&mut model, &floor, &patch, r, false, T)
        .unwrap()
        .shape;
    assert_eq!(model.kind_of(&round).unwrap(), ShapeType::Face);
    usable(&model, &round);
    let on_floor = ends_on_floor(&model, &round);
    assert_eq!(on_floor.len(), 2, "two corners on the floor");
    for p in on_floor {
        assert!(
            (p.y.abs() - 5.0).abs() < 2e-4,
            "a corner off the floor's edge: {p:?}"
        );
    }
}

/// The quarter drum of [`drum_quarter`] over y from -8 to 8, its surface
/// the rational B-spline that restates the cylinder exactly, facing out
/// of its axis.
fn spline_drum_quarter(model: &mut Model) -> Shape {
    use ogeom::math::{ControlGrid, KnotVector, Weighted};
    let frame = Frame::new(Point::new(0.0, -8.0, 5.0), Direction::Y, Direction::X, T).unwrap();
    let circle = Curve::Circle(CircleCurve::new(Circle::new(frame, 3.0, T).unwrap()));
    let arc = circle.to_bspline_over((0.0, PI / 2.0), T).unwrap();
    let mut net = Vec::with_capacity(arc.control_points().len() * 2);
    for w in arc.control_points() {
        for along in [0.0, 16.0] {
            let p = w.point() + Vector::new(0.0, along, 0.0);
            net.push(Weighted::new(p, w.weight, T).unwrap());
        }
    }
    let count = arc.control_points().len();
    let patch = ogeom::geom::BSplineSurface::rational(
        arc.knots().clone(),
        KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap(),
        ControlGrid::new(net, count, 2).unwrap(),
    )
    .unwrap();
    let face = patch_face(model, patch);
    let (p, n) = face_normal(model, &face, T).unwrap();
    if n.dot(Vector::new(p.x, 0.0, p.z - 5.0)) > 0.0 {
        face
    } else {
        face.reversed()
    }
}

/// The floor and the drum of [`a_plane_and_a_separate_cylinder_trim_into_one_shell`]
/// with the drum's surface a B-spline: the surfaces do not cross (the
/// round bridges the gap), so the ball is guided by where the two surfaces
/// offset by the radius cross, the line its centre runs along. The marched
/// round lands on the closed form: it touches the floor along x = 4 and
/// the drum along (2.4, 3.2), and stays within the band's fit of the
/// cylinder of radius 2 about the line through (4, 0, 2).
#[test]
fn a_plane_and_a_spline_drum_across_a_gap_round_as_the_cylinder() {
    use ogeom::geom::Surface as _;
    let mut model = Model::new();
    let r = 2.0;
    let floor = floor_rectangle(&mut model, (-10.0, 10.0), (-5.0, 5.0));
    let drum = spline_drum_quarter(&mut model);
    let built = fillet_faces(&mut model, &floor, &drum, r, true, T).unwrap();
    let shell = built.shape.clone();
    assert_eq!(faces(&model, &shell).len(), 3);
    usable(&model, &shell);
    let three = faces(&model, &shell);
    assert_eq!(
        ogeom::algo::sew(&mut model, &three, T)
            .unwrap()
            .shells
            .len(),
        1
    );
    let kept_floor = face_holding(&model, &shell, Point::new(9.0, 0.0, 0.0)).unwrap();
    assert!(
        (area(&model, &kept_floor) - 6.0 * 10.0).abs() < 1e-6,
        "the floor keeps x ≥ 4"
    );
    let round = faces(&model, &shell)
        .into_iter()
        .find(|f| {
            !f.is_same(&kept_floor) && built.history.modified(&drum).iter().all(|k| !k.is_same(f))
        })
        .unwrap();
    marched_tangent(&model, &shell, &round);
    for vertex in explore_unique(&model, &round, ShapeType::Vertex).unwrap() {
        let p = point_of(&model, &vertex);
        let on_floor = (p.x - 4.0).hypot(p.z);
        let on_drum = (p.x - 2.4).hypot(p.z - 3.2);
        assert!(
            on_floor.min(on_drum) < 1e-9,
            "a corner off the closed form: {p:?}"
        );
        assert!((p.y.abs() - 5.0).abs() < 1e-9);
    }
    let data = model.node(&round).unwrap().data().as_face().unwrap();
    let band = model.geometry().surface(data.surface).unwrap();
    let ((u0, u1), (v0, v1)) = band.domain();
    for i in 0..=8 {
        for j in 0..=8 {
            let u = (u1 - u0).mul_add(f64::from(i) / 8.0, u0);
            let v = (v1 - v0).mul_add(f64::from(j) / 8.0, v0);
            let p = band.point_at(u, v, T).unwrap();
            let off = ((p.x - 4.0).hypot(p.z - 2.0) - r).abs();
            assert!(off < 2e-4, "the band leaves the cylinder by {off} at {p:?}");
        }
    }
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
