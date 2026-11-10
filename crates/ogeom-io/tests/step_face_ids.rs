//! The file's own ids for the faces a STEP read builds.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::Frame;
use ogeom_topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

/// The ids of the file's face entities, `ADVANCED_FACE` and `FACE_SURFACE`.
fn face_entities(text: &str) -> Vec<u64> {
    let mut ids = Vec::new();
    for statement in text.split(';') {
        let statement = statement.trim();
        let Some((name, body)) = statement.split_once('=') else {
            continue;
        };
        let body = body.trim_start();
        if (body.starts_with("ADVANCED_FACE(") || body.starts_with("FACE_SURFACE("))
            && let Ok(id) = name.trim().trim_start_matches('#').parse()
        {
            ids.push(id);
        }
    }
    ids
}

/// Every face of every solid read has exactly one id among the import's
/// faces, each id is a face entity of the file, no two faces share one,
/// and there are as many as the file has face entities, in id order.
fn faces_named_once(text: &str) {
    let import = ogeom_io::read_step(text, T).unwrap();
    let model = import.document.model();
    let entities = face_entities(text);
    assert_eq!(import.faces.len(), entities.len());
    assert!(
        import.faces.windows(2).all(|w| w[0].0 < w[1].0),
        "in id order"
    );
    for (id, _) in &import.faces {
        assert!(entities.contains(id), "#{id} is no face entity of the file");
    }
    for solid in &import.solids {
        for face in explore_unique(model, solid, ShapeType::Face).unwrap() {
            let named = import
                .faces
                .iter()
                .filter(|(_, f)| f.is_partner(&face))
                .count();
            assert_eq!(named, 1, "a face of the solid is named {named} times");
        }
    }
}

/// A box written and read back: six faces, six ids.
#[test]
fn a_written_box_reads_with_an_id_for_each_face() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T)
        .unwrap()
        .shape;
    let mut document = ogeom_doc::Document::over(model);
    document.add_part("part", solid);
    let text = ogeom_io::write_step(&document, T).unwrap();
    assert_eq!(face_entities(&text).len(), 6);
    faces_named_once(&text);
}

/// The NIST CTC 01 part: its 139 faces, each named by its `ADVANCED_FACE`.
#[test]
fn a_nist_part_reads_with_an_id_for_each_face() {
    let text = corpus("nist_ctc_01_asme1_rd.stp");
    assert_eq!(face_entities(&text).len(), 139);
    faces_named_once(&text);
}

/// A box written, its `ADVANCED_FACE`s then made `FACE_SURFACE`s, the
/// supertype with the same arguments: each face is named by its entity
/// all the same.
#[test]
fn face_surfaces_read_with_an_id_for_each_face() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T)
        .unwrap()
        .shape;
    let mut document = ogeom_doc::Document::over(model);
    document.add_part("part", solid);
    let text = ogeom_io::write_step(&document, T)
        .unwrap()
        .replace("ADVANCED_FACE(", "FACE_SURFACE(");
    assert_eq!(face_entities(&text).len(), 6);
    faces_named_once(&text);
}

/// Two solids of one shell: each face is held by both, and named once.
#[test]
fn faces_two_solids_share_are_named_once() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T)
        .unwrap()
        .shape;
    let mut document = ogeom_doc::Document::over(model);
    document.add_part("part", solid);
    let text = ogeom_io::write_step(&document, T).unwrap();
    let start = text.find("MANIFOLD_SOLID_BREP(").unwrap();
    let end = start + text[start..].find(';').unwrap();
    let brep = &text[start..end];
    let data_end = text.rfind("ENDSEC;").unwrap();
    let text = format!(
        "{}#999999={brep};\n{}",
        &text[..data_end],
        &text[data_end..]
    );
    let import = ogeom_io::read_step(&text, T).unwrap();
    assert_eq!(import.solids.len(), 2);
    assert_eq!(import.faces.len(), 6);
    faces_named_once(&text);
}
