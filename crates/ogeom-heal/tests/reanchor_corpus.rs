//! Healing's first customer: the torus fillets of the smallest NIST part.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;

const T: Tolerances = Tolerances::millimetres();

#[test]
fn reanchoring_makes_the_smallest_nist_part_whole() {
    let path = format!(
        "{}/../../tests/corpus/nist_ftc_11_asme1_rb.stp",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(path).unwrap();
    let mut import = ogeom_io::read_step(&text, T).unwrap();
    assert_eq!(
        import.report.warnings.len(),
        2,
        "the reader names the two misaligned fillets"
    );

    let healed =
        ogeom_heal::reanchor_periodic_rings(import.document.model_mut(), &import.solids[0], T)
            .unwrap();
    assert!(!healed.history.is_empty(), "something was healed");

    // Every face of the healed solid now triangulates, and the part
    // measures: the reader's warning became the healer's work order, and
    // the work is done.
    let fine = ogeom_mesh::Deflection {
        chord: 1e-2,
        ..ogeom_mesh::Deflection::default()
    };
    let faces = ogeom_topo::explore(
        import.document.model(),
        &healed.shape,
        ogeom_topo::Filter::OfType(ogeom_topo::ShapeType::Face),
    )
    .unwrap();
    assert_eq!(faces.len(), 6);
    for face in &faces {
        ogeom_mesh::triangulate(import.document.model(), face, fine, T)
            .unwrap_or_else(|e| panic!("a healed face fails to mesh: {e}"));
    }
    let shell = ogeom_topo::explore_unique(
        import.document.model(),
        &healed.shape,
        ogeom_topo::ShapeType::Shell,
    )
    .unwrap()
    .remove(0);
    assert!(
        ogeom_algo::is_shell_closed(import.document.model(), &shell).unwrap(),
        "the healed topology closes"
    );
    // The part measures: every revolution face was re-annotated with
    // window-coherent pcurves, boundary vertices anchor to their edges'
    // curves, and the volume integral has a watertight boundary.
    let props =
        ogeom_algo::volume_properties(import.document.model(), &healed.shape, fine, T).unwrap();
    eprintln!("REPORT healed volume {:.3} mm^3", props.mass);
    assert!(props.mass > 0.0, "the healed part encloses volume");
}
