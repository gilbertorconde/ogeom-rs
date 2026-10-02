//! Sheet bodies through STEP and IGES: open shells and lone faces written
//! out and read back, alone and beside a solid in one part, measured by
//! face count, area, volume and the side each face presents.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Point, Vector};
use ogeom_mesh::Deflection;
use ogeom_topo::{Filter, Model, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

/// An L of two segments, 10 then 6 long, extruded 5: an open shell of two
/// planar faces of 50 and 30.
fn open_l(model: &mut Model) -> Shape {
    let profile = ogeom_algo::make_polygon(
        model,
        &[
            Point::new(0.0, 0.0, 0.0),
            Point::new(10.0, 0.0, 0.0),
            Point::new(10.0, 6.0, 0.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    ogeom_algo::make_prism(model, &profile, Vector::new(0.0, 0.0, 5.0), T)
        .unwrap()
        .shape
}

/// A lone 10 by 10 face, standing clear of a box at the origin.
fn square_face(model: &mut Model) -> Shape {
    let profile = ogeom_algo::make_polygon(
        model,
        &[Point::new(30.0, 0.0, 0.0), Point::new(40.0, 0.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let edge = faces_or_edges(model, &profile, ShapeType::Edge).remove(0);
    ogeom_algo::make_prism(model, &edge, Vector::new(0.0, 0.0, 10.0), T)
        .unwrap()
        .shape
}

fn faces_or_edges(model: &Model, shape: &Shape, kind: ShapeType) -> Vec<Shape> {
    explore(model, shape, Filter::OfType(kind)).unwrap()
}

fn area(model: &Model, shape: &Shape) -> f64 {
    ogeom_algo::surface_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom_algo::volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// Each face's area and outward normal, smallest area first.
fn face_records(model: &Model, shape: &Shape) -> Vec<(f64, Vector)> {
    let mut out: Vec<(f64, Vector)> = faces_or_edges(model, shape, ShapeType::Face)
        .iter()
        .map(|face| {
            let (_, normal) = ogeom_algo::face_normal(model, face, T).unwrap();
            (area(model, face), normal.normalized(T).unwrap())
        })
        .collect();
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out
}

fn assert_same_faces(original: &[(f64, Vector)], read: &[(f64, Vector)]) {
    assert_eq!(original.len(), read.len(), "{original:?} against {read:?}");
    for ((a, n), (b, m)) in original.iter().zip(read) {
        assert!((a - b).abs() < 1e-6, "area {b} against {a}");
        assert!(n.dot(*m) > 1.0 - 1e-9, "normal {m:?} against {n:?}");
    }
}

fn document_of(model: Model, shape: Shape) -> ogeom_doc::Document {
    let mut document = ogeom_doc::Document::over(model);
    document.add_part("part", shape);
    document
}

/// The shapes of a document's parts.
fn part_shapes(document: &ogeom_doc::Document) -> Vec<Shape> {
    document
        .products()
        .filter_map(|(_, product)| match &product.kind {
            ogeom_doc::ProductKind::Part { shape } => Some(shape.clone()),
            ogeom_doc::ProductKind::Assembly { .. } => None,
        })
        .collect()
}

#[test]
fn an_open_l_sheet_round_trips_through_step_as_two_faces_and_no_solid() {
    let mut model = Model::new();
    let sheet = open_l(&mut model);
    let original = face_records(&model, &sheet);
    assert_eq!(original.len(), 2);
    assert!((original[0].0 - 30.0).abs() < 1e-9 && (original[1].0 - 50.0).abs() < 1e-9);

    let text = ogeom_io::write_step(&document_of(model, sheet), T).unwrap();
    assert!(text.contains("MANIFOLD_SURFACE_SHAPE_REPRESENTATION("));
    assert!(text.contains("SHELL_BASED_SURFACE_MODEL("));
    assert!(text.contains("OPEN_SHELL("));
    assert!(!text.contains("MANIFOLD_SOLID_BREP("));

    let import = ogeom_io::read_step(&text, T).unwrap();
    assert!(import.solids.is_empty());
    assert_eq!(import.shells.len(), 1);
    let model = import.document.model();
    let parts = part_shapes(&import.document);
    assert_eq!(parts.len(), 1, "one part");
    assert_same_faces(&original, &face_records(model, &parts[0]));
    // The two faces share their corner edge again: one shell, not two.
    assert_eq!(model.kind_of(&import.shells[0]).unwrap(), ShapeType::Shell);
}

#[test]
fn a_box_and_a_square_sheet_in_one_part_round_trip_through_step() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let face = square_face(&mut model);
    let compound = ogeom_algo::make_compound(&mut model, &[solid, face.clone()])
        .unwrap()
        .shape;
    let sheet_faces = face_records(&model, &face);

    let text = ogeom_io::write_step(&document_of(model, compound), T).unwrap();
    assert!(text.contains("ADVANCED_BREP_SHAPE_REPRESENTATION("));
    assert!(text.contains("MANIFOLD_SURFACE_SHAPE_REPRESENTATION("));

    let import = ogeom_io::read_step(&text, T).unwrap();
    let model = import.document.model();
    assert_eq!(part_shapes(&import.document).len(), 1, "one part");
    assert_eq!(import.solids.len(), 1);
    assert_eq!(import.shells.len(), 1);
    assert!((volume(model, &import.solids[0]) - 1000.0).abs() < 1e-6);
    assert_same_faces(&sheet_faces, &face_records(model, &import.shells[0]));
    // The part holds both: six faces of the box and the sheet's one.
    let part = &part_shapes(&import.document)[0];
    assert_eq!(faces_or_edges(model, part, ShapeType::Face).len(), 7);
}

#[test]
fn an_open_l_sheet_round_trips_through_iges_as_two_faces_and_no_solid() {
    let mut model = Model::new();
    let sheet = open_l(&mut model);
    let original = face_records(&model, &sheet);

    let text = ogeom_io::write_iges(&document_of(model, sheet), T).unwrap();
    let import = ogeom_io::read_iges(&text, T).unwrap();
    assert!(import.solids.is_empty());
    assert_eq!(
        import.sheets.len(),
        1,
        "the two faces sew back into one sheet"
    );
    assert_same_faces(
        &original,
        &face_records(import.document.model(), &import.sheets[0]),
    );
    assert!(
        import.report.skipped.is_empty(),
        "{:?}",
        import.report.skipped
    );
}

#[test]
fn a_box_and_a_square_sheet_round_trip_through_iges() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let face = square_face(&mut model);
    let compound = ogeom_algo::make_compound(&mut model, &[solid, face.clone()])
        .unwrap()
        .shape;
    let sheet_faces = face_records(&model, &face);

    let text = ogeom_io::write_iges(&document_of(model, compound), T).unwrap();
    let import = ogeom_io::read_iges(&text, T).unwrap();
    let model = import.document.model();
    assert_eq!(import.solids.len(), 1);
    assert_eq!(import.sheets.len(), 1);
    assert!((volume(model, &import.solids[0]) - 1000.0).abs() < 1e-6);
    assert_same_faces(&sheet_faces, &face_records(model, &import.sheets[0]));
}

/// A reversed face presents its other side, and the file says so: STEP by
/// the face's sense, IGES by turning the surface over. A planar face and a
/// cylinder wall, the second going to IGES as an exchanged B-spline patch.
#[test]
fn a_reversed_sheet_face_keeps_its_side_through_step_and_iges() {
    for format in ["step", "iges"] {
        // The plane.
        let mut model = Model::new();
        let face = square_face(&mut model).reversed();
        let original = face_records(&model, &face);
        let document = document_of(model, face);
        let read = round_trip(&document, format);
        let model = read.model();
        let parts = part_shapes(&read);
        assert_same_faces(&original, &face_records(model, &parts[0]));

        // The cylinder wall, which faces the axis once reversed.
        let mut model = Model::new();
        let cylinder = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 4.0, 7.0, T)
            .unwrap()
            .shape;
        let wall = faces_or_edges(&model, &cylinder, ShapeType::Face)
            .into_iter()
            .find(|f| (area(&model, f) - std::f64::consts::TAU * 4.0 * 7.0).abs() < 1e-6)
            .expect("the cylinder has a wall")
            .reversed();
        let wall_area = area(&model, &wall);
        assert!(
            radial_side(&model, &wall) < 0.0,
            "the reversed wall faces in"
        );
        let read = round_trip(&document_of(model, wall), format);
        let model = read.model();
        let part = &part_shapes(&read)[0];
        let faces = faces_or_edges(model, part, ShapeType::Face);
        assert_eq!(faces.len(), 1, "{format}");
        assert!(
            (area(model, &faces[0]) - wall_area).abs() < 1e-6 * wall_area,
            "{format}: {} against {wall_area}",
            area(model, &faces[0])
        );
        assert!(radial_side(model, &faces[0]) < 0.0, "{format}: faces in");
    }
}

/// Which way a face round the world z axis presents: positive away from
/// the axis, negative towards it.
fn radial_side(model: &Model, face: &Shape) -> f64 {
    let (at, normal) = ogeom_algo::face_normal(model, face, T).unwrap();
    let radial = Vector::new(at.x, at.y, 0.0);
    radial.dot(normal)
}

fn round_trip(document: &ogeom_doc::Document, format: &str) -> ogeom_doc::Document {
    if format == "step" {
        let text = ogeom_io::write_step(document, T).unwrap();
        ogeom_io::read_step(&text, T).unwrap().document
    } else {
        let text = ogeom_io::write_iges(document, T).unwrap();
        ogeom_io::read_iges(&text, T).unwrap().document
    }
}

#[test]
fn wireframe_and_empty_parts_are_refused_by_name() {
    for format in ["step", "iges"] {
        let write = |document: &ogeom_doc::Document| {
            if format == "step" {
                ogeom_io::write_step(document, T).map(|_| ())
            } else {
                ogeom_io::write_iges(document, T).map(|_| ())
            }
        };

        // A box with a free wire beside it: the wire is not dropped.
        let mut model = Model::new();
        let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T)
            .unwrap()
            .shape;
        let wire = ogeom_algo::make_polygon(
            &mut model,
            &[Point::new(5.0, 0.0, 0.0), Point::new(6.0, 0.0, 0.0)],
            false,
            T,
        )
        .unwrap()
        .shape;
        let compound = ogeom_algo::make_compound(&mut model, &[solid, wire])
            .unwrap()
            .shape;
        let error = write(&document_of(model, compound))
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("free wire") && error.contains("wireframe"),
            "{error}"
        );

        // A compound of nothing.
        let mut model = Model::new();
        let empty = ogeom_algo::make_compound(&mut model, &[]).unwrap().shape;
        let error = write(&document_of(model, empty)).unwrap_err().to_string();
        assert!(error.contains("no solid, shell or face"), "{error}");
    }
}

#[test]
fn a_coloured_sheet_keeps_its_colour_through_step() {
    let mut model = Model::new();
    let sheet = open_l(&mut model);
    let teal = ogeom_doc::Colour::rgb(0.0, 0.5, 0.5);
    let mut document = document_of(model, sheet.clone());
    document.set_colour(&sheet, teal);

    let text = ogeom_io::write_step(&document, T).unwrap();
    let import = ogeom_io::read_step(&text, T).unwrap();
    let colour = import.document.colour_of(&import.shells[0]).unwrap();
    assert!(
        (colour.r - teal.r).abs() < 1e-6
            && (colour.g - teal.g).abs() < 1e-6
            && (colour.b - teal.b).abs() < 1e-6,
        "{colour:?}"
    );
}
