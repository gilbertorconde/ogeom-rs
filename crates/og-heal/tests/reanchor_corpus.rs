//! Healing's first customer: the torus fillets of the smallest NIST part.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use og_core::Tolerances;

const T: Tolerances = Tolerances::millimetres();

#[test]
fn reanchoring_makes_the_smallest_nist_part_whole() {
    let path = format!(
        "{}/../../tests/corpus/nist_ftc_11_asme1_rb.stp",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(path).unwrap();
    let mut import = og_io::read_step(&text, T).unwrap();
    assert_eq!(
        import.report.warnings.len(),
        2,
        "the reader names the two misaligned fillets"
    );

    let healed = og_heal::reanchor_periodic_rings(&mut import.model, &import.solids[0], T).unwrap();
    assert!(!healed.history.is_empty(), "something was healed");

    // Every face of the healed solid now triangulates, and the part
    // measures: the reader's warning became the healer's work order, and
    // the work is done.
    let fine = og_mesh::Deflection {
        chord: 1e-2,
        ..og_mesh::Deflection::default()
    };
    let faces = og_topo::explore(
        &import.model,
        &healed.shape,
        og_topo::Filter::OfType(og_topo::ShapeType::Face),
    )
    .unwrap();
    assert_eq!(faces.len(), 6);
    for face in &faces {
        og_mesh::triangulate(&import.model, face, fine, T)
            .unwrap_or_else(|e| panic!("a healed face fails to mesh: {e}"));
    }
    let shell = og_topo::explore_unique(&import.model, &healed.shape, og_topo::ShapeType::Shell)
        .unwrap()
        .remove(0);
    assert!(
        og_algo::is_shell_closed(&import.model, &shell).unwrap(),
        "the healed topology closes"
    );
    // The part measures: every revolution face was re-annotated with
    // window-coherent pcurves, boundary vertices anchor to their edges'
    // curves, and the volume integral has a watertight boundary.
    let props = og_algo::volume_properties(&import.model, &healed.shape, fine, T).unwrap();
    eprintln!("REPORT healed volume {:.3} mm^3", props.mass);
    assert!(props.mass > 0.0, "the healed part encloses volume");
}
