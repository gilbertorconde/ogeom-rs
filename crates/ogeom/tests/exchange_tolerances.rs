//! What a reader builds keeps the tolerance containment rule, what the
//! healer is handed it restores, and what the IGES writer writes its reader
//! reads back.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{Severity, check, restore_containment};
use ogeom::core::{Tolerance, Tolerances};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

fn containment_findings(model: &Model, shape: &Shape) -> usize {
    check(model, shape, T)
        .unwrap()
        .of(Severity::Broken)
        .iter()
        .filter(|p| p.what.contains("tighter than the"))
        .count()
}

fn volume(model: &Model, shape: &Shape) -> Option<f64> {
    ogeom::algo::volume_properties(model, shape, Deflection::default(), T)
        .ok()
        .map(|p| p.mass)
}

/// The STEP reader widens an edge to how far its pcurves sit from its
/// curve; the edge's vertices widen with it, so a freshly read solid holds
/// every vertex at least as loose as the edges it bounds. These files read
/// with hundreds of vertices tighter than their edges otherwise.
#[test]
fn a_step_read_keeps_tolerance_containment() {
    for name in [
        "nist_ctc_02_asme1_rc.stp",
        "nist_ctc_05_asme1_rd.stp",
        "long_bore_between_cross_holes.step",
    ] {
        let import = ogeom::io::read_step(&corpus(name), T).unwrap();
        assert!(!import.solids.is_empty(), "{name}");
        for solid in &import.solids {
            assert_eq!(
                containment_findings(import.document.model(), solid),
                0,
                "{name}: a vertex tighter than an edge it bounds"
            );
        }
    }
}

/// The healer restores containment. A real part whose edges genuinely
/// need their width — pcurves sitting microns off their curves — has its
/// vertices reset to the confusion tolerance, which is the state a reader
/// that widened only the edges left behind. The reduction cannot tighten
/// those edges, so the vertices are what must grow, and the report says
/// how many did.
#[test]
fn fix_shape_restores_tolerance_containment() {
    let mut import = ogeom::io::read_step(&corpus("nist_ctc_02_asme1_rc.stp"), T).unwrap();
    let solid = import.solids[0].clone();
    let model = import.document.model_mut();
    let tight = Tolerance::new(T.confusion()).unwrap();
    for vertex in explore_unique(model, &solid, ShapeType::Vertex).unwrap() {
        if let Some(node) = model.node_mut(&vertex)
            && let NodeData::Vertex(data) = node.data_mut()
        {
            data.tolerance = tight;
        }
    }
    assert!(containment_findings(model, &solid) > 0);

    let fixed = ogeom::heal::fix_shape(model, &solid, T).unwrap();
    assert_eq!(containment_findings(model, &fixed.shape), 0);
    assert!(fixed.report.tolerances_widened > 0);
    assert!(
        fixed.report.after.of(Severity::Broken).is_empty(),
        "{}",
        fixed.report.after
    );

    // The pass alone, run again, has nothing left to grow.
    assert_eq!(restore_containment(model, &fixed.shape).unwrap(), 0);
}

/// Whatever the IGES writer writes, its reader reads: every record exactly
/// eighty columns, and every solid back, valid, at the volume it went out
/// with. One of these parts carries a coefficient near 1e-51, which
/// positional notation spells longer than a record; the others have
/// vertices their curves miss by a few nanometres, or project a hair past
/// a bounded curve's end.
#[test]
fn a_solid_written_as_iges_reads_back() {
    for name in [
        "nist_ctc_03_asme1_rc.stp",
        "nist_ctc_01_asme1_rd.stp",
        "nist_ftc_07_asme1_rd.stp",
        "closed_tube_face_crosses_its_join.step",
        "fillet_strip_with_a_narrow_chart.step",
    ] {
        let read = ogeom::io::read_step(&corpus(name), T).unwrap();
        let iges = ogeom::io::write_iges(&read.document, T).unwrap();
        assert!(
            iges.lines().all(|line| line.len() == 80),
            "{name}: a record not eighty columns"
        );
        let back = ogeom::io::read_iges(&iges, T).unwrap();
        assert_eq!(back.solids.len(), read.solids.len(), "{name}");
        for (sent, came) in read.solids.iter().zip(&back.solids) {
            let model = back.document.model();
            assert!(check(model, came, T).unwrap().is_usable(), "{name}");
            assert_eq!(containment_findings(model, came), 0, "{name}");
            let faces = |m: &Model, s: &Shape| explore_unique(m, s, ShapeType::Face).unwrap().len();
            assert_eq!(
                faces(read.document.model(), sent),
                faces(model, came),
                "{name}"
            );
            // A part that encloses a volume encloses the same one.
            if let (Some(a), Some(b)) = (volume(read.document.model(), sent), volume(model, came)) {
                assert!(
                    (a - b).abs() <= a.abs() * 1e-3,
                    "{name}: volume {a} went out, {b} came back"
                );
            }
        }
    }
}
