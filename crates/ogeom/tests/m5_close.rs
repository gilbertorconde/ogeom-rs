//! M5's closing bar, both criteria end to end.
//!
//! A multi-body STEP assembly imported with colours, names and PMI,
//! modified, exported, reimported, and compared against what went in; and a
//! 2D drawing generated from a 3D model — hidden lines removed, visible and
//! hidden edges classified, sections taken.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::doc::{Colour, MeasureKind, ProductKind};
use ogeom::hlr::{Source, View, Visibility};
use ogeom::math::{Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

fn fine() -> Deflection {
    Deflection {
        chord: 1e-3,
        ..Deflection::default()
    }
}

#[test]
fn an_assembly_imports_modifies_exports_and_returns_the_same() {
    // Import: products, names, colours and PMI all arrive.
    let text = corpus("ogeom_asm_bolted_plate.stp");
    let mut import = ogeom::io::read_step(&text, T).unwrap();
    let names: Vec<String> = import
        .document
        .products()
        .map(|(_, p)| p.name.clone())
        .collect();
    assert_eq!(names, ["plate", "bolt", "bolted-plate"]);
    assert_eq!(import.document.pmi().dimensions.len(), 1);
    assert_eq!(import.document.pmi().tolerances.len(), 1);
    assert_eq!(import.document.pmi().datums.len(), 1);
    let plate = import.document.products().next().unwrap().0;
    let bolt = import.document.products().nth(1).unwrap().0;
    let assembly = import.document.products().nth(2).unwrap().0;

    // Modify, three ways: drill the plate, add a third bolt, annotate.
    let plate_shape = match &import.document.get(plate).unwrap().kind {
        ProductKind::Part { shape } => shape.clone(),
        ProductKind::Assembly { .. } => unreachable!(),
    };
    let drill_frame =
        Frame::new(Point::new(20.0, 20.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill =
        ogeom::algo::make_cylinder(import.document.model_mut(), drill_frame, 3.0, 7.0, T).unwrap();
    let drilled =
        ogeom::boolean::cut(import.document.model_mut(), &plate_shape, &drill.shape, T).unwrap();
    import
        .document
        .replace_part_shape(plate, drilled.shape.clone())
        .unwrap();
    import
        .document
        .add_instance(
            assembly,
            bolt,
            ogeom::math::Transform::translation(Vector::new(10.0, 30.0, 5.0)),
            Some("bolt-3".into()),
        )
        .unwrap();
    import
        .document
        .set_colour(&drilled.shape, Colour::rgb(0.5, 0.5, 0.1));
    import
        .document
        .pmi_mut()
        .dimensions
        .push(ogeom::doc::Dimension {
            name: "diameter".into(),
            values: vec![6.0],
            kind: MeasureKind::Length,
            plus: Some(0.05),
            minus: Some(-0.05),
            features: Vec::new(),
            location: false,
        });

    // Export, reimport, compare against what went in.
    let written = ogeom::io::write_step(&import.document, T).unwrap();
    let back = ogeom::io::read_step(&written, T).unwrap();

    let names: Vec<String> = back
        .document
        .products()
        .map(|(_, p)| p.name.clone())
        .collect();
    assert_eq!(names, ["plate", "bolt", "bolted-plate"]);
    let root = back.document.roots()[0];
    let mut paths: Vec<String> = back
        .document
        .occurrences_of(root)
        .unwrap()
        .iter()
        .map(|o| o.path.clone())
        .collect();
    paths.sort();
    assert_eq!(
        paths,
        [
            "bolted-plate/bolt-1",
            "bolted-plate/bolt-2",
            "bolted-plate/bolt-3",
            "bolted-plate/plate-1"
        ]
    );

    // The drilled volume survives: plate minus the bore, plus three bolts.
    let total: f64 = back
        .document
        .occurrences_of(root)
        .unwrap()
        .iter()
        .map(|o| {
            ogeom::algo::volume_properties(back.document.model(), &o.shape, fine(), T)
                .unwrap()
                .mass
        })
        .sum();
    let plate_drilled = 40.0f64 * 40.0 * 5.0 - core::f64::consts::PI * 9.0 * 5.0;
    let expected = 3.0f64.mul_add(8.0 * 8.0 * 20.0, plate_drilled);
    assert!(
        (total - expected).abs() < 0.05,
        "assembly volume {total} against {expected}"
    );

    // The annotations went with it: the original dimension plus the added
    // one, the flatness, the datum, and the new plate colour.
    assert_eq!(back.document.pmi().dimensions.len(), 2);
    assert_eq!(back.document.pmi().tolerances.len(), 1);
    assert_eq!(back.document.pmi().tolerances[0].kind, "flatness");
    assert_eq!(back.document.pmi().datums.len(), 1);
    assert_eq!(back.document.pmi().datums[0].label, "A");
    let plate_occ = back
        .document
        .occurrences_of(root)
        .unwrap()
        .into_iter()
        .find(|o| o.path.ends_with("plate-1"))
        .unwrap();
    assert_eq!(
        back.document
            .resolved_colour(plate_occ.part, &plate_occ.shape),
        Some(Colour::rgb(0.5, 0.5, 0.1))
    );
}

#[test]
fn a_drawing_comes_off_the_model_with_hidden_lines_and_a_section() {
    // The model: a plate with a bore through it.
    let mut model = ogeom::topo::Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 6.0, 4.0), T).unwrap();
    let bore_frame = Frame::new(Point::new(5.0, 3.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let bore = ogeom::algo::make_cylinder(&mut model, bore_frame, 1.0, 6.0, T).unwrap();
    let part = ogeom::boolean::cut(&mut model, &block.shape, &bore.shape, T).unwrap();

    // The three-quarter view: edges classified, both classes present, and
    // the bore visible only as far as the eye can see into it.
    let view = View::looking(Vector::new(-1.0, -1.2, -0.9), Vector::new(0.0, 0.0, 1.0), T).unwrap();
    let drawing = ogeom::hlr::project(&model, &part.shape, &view, fine(), T).unwrap();
    assert!(!drawing.visible.is_empty(), "visible edges classified");
    assert!(!drawing.hidden.is_empty(), "hidden edges classified");
    let silhouettes = drawing
        .curves()
        .filter(|c| matches!(c.source, Source::Silhouette))
        .count();
    assert!(silhouettes > 0, "the bore's wall draws by silhouette");
    let hidden_edges = drawing
        .hidden
        .iter()
        .filter(|c| matches!(c.source, Source::Edge(_)))
        .count();
    assert!(hidden_edges > 0, "the far side is dashed, not gone");
    assert!(
        drawing
            .curves()
            .all(|c| matches!(c.visibility, Visibility::Visible | Visibility::Hidden)),
        "every curve is classified"
    );

    // The section, half a radius off the bore's axis: revealed material is
    // the rectangle minus the chord-width slot — arithmetic, not opinion.
    let plane =
        Plane::new(Frame::new(Point::new(0.0, 2.5, 0.0), Direction::Y, Direction::Z, T).unwrap());
    let section = ogeom::hlr::section(&mut model, &part.shape, &plane, fine(), T).unwrap();
    let area: f64 = section
        .outline
        .iter()
        .map(|l| {
            let mut sum = 0.0;
            for pair in l.windows(2) {
                sum += pair[0].x.mul_add(pair[1].y, -(pair[1].x * pair[0].y));
            }
            if let (Some(first), Some(last)) = (l.first(), l.last()) {
                sum += last.x.mul_add(first.y, -(first.x * last.y));
            }
            sum / 2.0
        })
        .sum::<f64>()
        .abs();
    let chord = 2.0 * (1.0_f64 - 0.25).sqrt();
    let expected = 10.0f64.mul_add(4.0, -(chord * 4.0));
    assert!(
        (area - expected).abs() < 1e-3,
        "section reveals {area} against {expected}"
    );
    assert!(!section.drawing.visible.is_empty());
}
