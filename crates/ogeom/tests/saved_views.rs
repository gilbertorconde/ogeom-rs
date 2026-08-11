//! §E3 of `docs/PLAN.md`: saved views and notes. An annotated part organises
//! its PMI presentation into named views; the views survive the native
//! document format and STEP, and notes survive the document.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::topo::Model;

const T: Tolerances = Tolerances::millimetres();

fn annotated_document() -> ogeom::doc::Document {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 10.0, 5.0), T)
        .unwrap()
        .shape;
    let mut document = ogeom::doc::Document::over(model);
    let part = document.add_part("block", block);

    // Two drawn callouts, no semantics behind them — a note-shaped drawing
    // and a frame-shaped one.
    for (name, y) in [("width", 0.0), ("height", 4.0)] {
        document.pmi_mut().callouts.push(ogeom::doc::Callout {
            name: name.into(),
            plane: Some(Frame::WORLD),
            polylines: vec![vec![
                Point::new(0.0, y, 6.0),
                Point::new(20.0, y, 6.0),
            ]],
            annotates: None,
        });
    }
    let front = Frame::new(
        Point::new(10.0, -30.0, 2.5),
        Direction::Y,
        Direction::X,
        T,
    )
    .unwrap();
    let top = Frame::new(Point::new(10.0, 5.0, 40.0), Direction::Z, Direction::X, T).unwrap();
    document.add_view(ogeom::doc::View {
        name: "front".into(),
        frame: front,
        clipping: None,
        callouts: vec![0],
    });
    document.add_view(ogeom::doc::View {
        name: "top".into(),
        frame: top,
        clipping: None,
        callouts: vec![0, 1],
    });
    document.add_note(ogeom::doc::Note {
        author: "gil".into(),
        text: "check the width against the fixture".into(),
        product: Some(part),
    });
    document
}

/// The native document format carries views, notes and the callouts the
/// views index — written, read, and equal field for field.
#[test]
fn views_and_notes_survive_the_native_document() {
    let document = annotated_document();
    let text = ogeom::io::native::write_document(&document, ogeom::io::native::WriteOptions::default())
        .unwrap();
    let back = ogeom::io::native::read_document(&text).unwrap();

    assert_eq!(back.pmi().callouts.len(), 2);
    assert_eq!(back.views().len(), 2);
    assert_eq!(back.views()[0].name, "front");
    assert_eq!(back.views()[0].callouts, vec![0]);
    assert_eq!(back.views()[1].callouts, vec![0, 1]);
    assert!(
        back.views()[1]
            .frame
            .origin()
            .distance(Point::new(10.0, 5.0, 40.0))
            < 1e-9
    );
    assert_eq!(back.notes().len(), 1);
    assert_eq!(back.notes()[0].author, "gil");
    assert!(back.notes()[0].product.is_some());
}

/// STEP carries the views as the named draughting models they are there:
/// name, camera placement, and which callouts each presents.
#[test]
fn views_survive_step() {
    let document = annotated_document();
    let text = ogeom::io::write_step(&document, T).unwrap();
    let import = ogeom::io::read_step(&text, T).unwrap();
    let back = &import.document;

    assert_eq!(back.pmi().callouts.len(), 2, "{:?}", import.report.warnings);
    assert_eq!(back.views().len(), 2);
    let top = back
        .views()
        .iter()
        .find(|v| v.name == "top")
        .expect("the top view survives by name");
    assert!(top.frame.origin().distance(Point::new(10.0, 5.0, 40.0)) < 1e-9);
    assert_eq!(top.callouts.len(), 2);
    let front = back.views().iter().find(|v| v.name == "front").unwrap();
    assert_eq!(front.callouts.len(), 1);
    assert_eq!(back.pmi().callouts[front.callouts[0]].name, "width");
}
