//! Callouts written to STEP and read back: the words a note carries, and
//! the part and faces a note or a dimension is about.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::doc::{Annotated, Callout, Dimension, Document, MeasureKind};
use ogeom::math::{Direction, Frame, Point};
use ogeom::topo::{ShapeType, TShapeId, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn note(name: &str, text: &str, items: Vec<TShapeId>, annotates: Option<Annotated>) -> Callout {
    Callout {
        name: name.into(),
        plane: Some(Frame::WORLD),
        polylines: vec![vec![Point::new(0.0, 0.0, 12.0), Point::new(8.0, 0.0, 12.0)]],
        annotates,
        items,
        text: Some(text.into()),
    }
}

/// A note's words are written as a text literal beside its strokes and
/// come back as the note's text, the strokes as they were.
#[test]
fn a_note_reads_back_with_its_words() {
    let mut document = Document::new();
    let block = ogeom::algo::make_box(document.model_mut(), Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    document.add_part("block", block);
    document
        .pmi_mut()
        .callouts
        .push(note("note", "Deburr this edge", Vec::new(), None));
    let text = ogeom::io::write_step(&document, T).unwrap();
    let back = ogeom::io::read_step(&text, T).unwrap();
    let callouts = &back.document.pmi().callouts;
    assert_eq!(callouts.len(), 1);
    assert_eq!(callouts[0].text.as_deref(), Some("Deburr this edge"));
    assert_eq!(callouts[0].polylines.len(), 1);
    assert_eq!(callouts[0].polylines[0].len(), 2);
    assert!(callouts[0].polylines[0][1].distance(Point::new(8.0, 0.0, 12.0)) < 1e-9);
}

/// Two parts, and a dimension and a note about a face of the second: each
/// is anchored to the part that owns its face, and both read back on that
/// part's face, whichever part the file wrote first.
#[test]
fn a_callout_reads_back_on_the_part_it_is_about() {
    let mut document = Document::new();
    let first = ogeom::algo::make_box(document.model_mut(), Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let beside = Frame::new(Point::new(50.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let second = ogeom::algo::make_box(document.model_mut(), beside, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let bore = explore_unique(document.model(), &second, ShapeType::Face).unwrap()[0].node();
    document.add_part("A", first);
    document.add_part("B", second);
    document.pmi_mut().dimensions.push(Dimension {
        name: "size".into(),
        values: vec![10.0],
        kind: MeasureKind::Length,
        plus: None,
        minus: None,
        features: vec![vec![bore]],
        location: false,
    });
    document.pmi_mut().callouts.push(note(
        "size",
        "10",
        Vec::new(),
        Some(Annotated::Dimension(0)),
    ));
    document
        .pmi_mut()
        .callouts
        .push(note("thread", "M6 thread", vec![bore], None));

    let text = ogeom::io::write_step(&document, T).unwrap();
    let back = ogeom::io::read_step(&text, T).unwrap();
    let model = back.document.model();
    let faces_of = |at: usize| -> Vec<TShapeId> {
        explore_unique(model, &back.solids[at], ShapeType::Face)
            .unwrap()
            .iter()
            .map(|f| f.node())
            .collect()
    };
    // The second part is the one beside the first.
    let second_at = (0..back.solids.len())
        .find(|&at| {
            ogeom::algo::tight_bounds(model, &back.solids[at], T)
                .unwrap()
                .low()
                .unwrap()
                .x
                > 40.0
        })
        .unwrap();
    let on_second = faces_of(second_at);
    let pmi = back.document.pmi();
    assert!(
        pmi.dimensions[0]
            .items()
            .all(|item| on_second.contains(&item)),
        "the dimension reads back on the second part"
    );
    let thread = pmi
        .callouts
        .iter()
        .find(|c| c.text.as_deref() == Some("M6 thread"))
        .unwrap();
    assert!(!thread.items.is_empty(), "the note names its face");
    assert!(thread.items.iter().all(|item| on_second.contains(item)));
    // Every aspect hangs from the product definition its faces are in: the
    // file names the second part's definition shape for the second part's
    // faces, and the first part's nowhere.
    let lines: Vec<&str> = text.lines().collect();
    let part_a = lines
        .iter()
        .position(|l| l.contains("PRODUCT('A','A'"))
        .unwrap();
    let pds_a = lines[part_a..]
        .iter()
        .find(|l| l.contains("PRODUCT_DEFINITION_SHAPE("))
        .and_then(|l| l.split('=').next())
        .unwrap();
    let on_a = format!("SHAPE_ASPECT('','',{pds_a},");
    let aspects = text.matches("SHAPE_ASPECT('','',").count();
    assert!(aspects >= 2, "{aspects} aspects");
    assert!(
        !text.contains(&on_a),
        "an aspect of the second part's face hangs from the first part"
    );
}

/// The native format keeps a note's words and the faces it is about.
#[test]
fn a_note_survives_the_native_format() {
    let mut document = Document::new();
    let block = ogeom::algo::make_box(document.model_mut(), Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let face = explore_unique(document.model(), &block, ShapeType::Face).unwrap()[2].clone();
    document.add_part("block", block);
    document
        .pmi_mut()
        .callouts
        .push(note("thread", "M6 thread", vec![face.node()], None));
    let text =
        ogeom::io::native::write_document(&document, ogeom::io::native::WriteOptions::default())
            .unwrap();
    let back = ogeom::io::native::read_document(&text).unwrap();
    let callout = &back.pmi().callouts[0];
    assert_eq!(callout.text.as_deref(), Some("M6 thread"));
    assert_eq!(callout.items.len(), 1);
    let model = back.model();
    let faces: Vec<TShapeId> = back
        .products()
        .flat_map(|(_, product)| match &product.kind {
            ogeom::doc::ProductKind::Part { shape } => {
                explore_unique(model, shape, ShapeType::Face)
                    .unwrap()
                    .iter()
                    .map(|f| f.node())
                    .collect::<Vec<_>>()
            }
            ogeom::doc::ProductKind::Assembly { .. } => Vec::new(),
        })
        .collect();
    assert!(faces.contains(&callout.items[0]));
}
