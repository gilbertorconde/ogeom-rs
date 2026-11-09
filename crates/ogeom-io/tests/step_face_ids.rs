//! Each face of a read solid answers with the file's id for it: the face
//! entity it was built from, found by partner among the import's faces.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_doc::Document;
use ogeom_io::step::parse::parse;
use ogeom_math::Frame;
use ogeom_topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// Every face of every solid has exactly one id among the import's faces,
/// each id is a face entity of the file, and no two faces share one.
fn every_face_has_its_entity(text: &str) -> usize {
    let import = ogeom_io::read_step(text, T).unwrap();
    let exchange = parse(text).unwrap();
    let model = import.document.model();
    let mut seen = Vec::new();
    for solid in &import.solids {
        for face in explore_unique(model, solid, ShapeType::Face).unwrap() {
            let ids: Vec<u64> = import
                .faces
                .iter()
                .filter(|(_, f)| f.is_partner(&face))
                .map(|(id, _)| *id)
                .collect();
            assert_eq!(ids.len(), 1, "a face answers with one id, not {ids:?}");
            let entity = &exchange.data[&ids[0]];
            assert!(
                entity.part("ADVANCED_FACE").is_some() || entity.part("FACE_SURFACE").is_some(),
                "#{} is {}",
                ids[0],
                entity.keyword()
            );
            seen.push(ids[0]);
        }
    }
    let count = seen.len();
    seen.sort_unstable();
    seen.dedup();
    assert_eq!(seen.len(), count, "no two faces share an id");
    assert!(
        import.faces.windows(2).all(|w| w[0].0 < w[1].0),
        "faces come in id order"
    );
    count
}

#[test]
fn the_faces_of_a_written_box_carry_their_entity_ids() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (3.0, 2.0, 1.0), T)
        .unwrap()
        .shape;
    let mut document = Document::over(model);
    document.add_part("box", solid);
    let text = ogeom_io::write_step(&document, T).unwrap();
    assert_eq!(every_face_has_its_entity(&text), 6);
}

#[test]
fn the_faces_of_a_corpus_part_carry_their_entity_ids() {
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../tests/corpus/nist_ctc_01_asme1_rd.stp"
    );
    let text = std::fs::read_to_string(path).unwrap();
    let faces = every_face_has_its_entity(&text);
    let written = text.matches("ADVANCED_FACE(").count();
    assert_eq!(
        faces, written,
        "every face entity of the file is a face read"
    );
}
