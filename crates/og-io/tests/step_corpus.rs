//! The corpus, finally consumed: reading the NIST test parts.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use og_core::Tolerances;

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

#[test]
fn the_smallest_nist_part_reads_into_a_closed_solid() {
    let text = corpus("nist_ftc_11_asme1_rb.stp");
    let import = og_io::read_step(&text, T).unwrap();
    eprintln!("REPORT scale={}mm", import.report.scale_mm);
    for w in &import.report.warnings {
        eprintln!("REPORT warn: {w}");
    }
    for (k, n) in &import.report.skipped {
        eprintln!("REPORT skipped {n:4}  {k}");
    }
    assert_eq!(import.solids.len(), 1);
    assert!((import.report.scale_mm - 1.0).abs() < 1e-12, "millimetres");

    // The topology closes even where meshing cannot yet follow.
    let solid = &import.solids[0];
    let shell = og_topo::explore_unique(&import.model, solid, og_topo::ShapeType::Shell)
        .unwrap()
        .remove(0);
    assert!(og_algo::is_shell_closed(&import.model, &shell).unwrap());

    // Four of the six faces triangulate — planes and full cylinder bands,
    // the bands through synthesised seams. The two torus fillets are the
    // honest remainder: their two ring vertices sit at different angles, so
    // no seam can join them without re-anchoring a shared edge, which is
    // healing's first named import case. The warnings say exactly that.
    let fine = og_mesh::Deflection {
        chord: 1e-2,
        ..og_mesh::Deflection::default()
    };
    let mut meshed = 0;
    for face in og_topo::explore(
        &import.model,
        solid,
        og_topo::Filter::OfType(og_topo::ShapeType::Face),
    )
    .unwrap()
    {
        if og_mesh::triangulate(&import.model, &face, fine, T).is_ok() {
            meshed += 1;
        }
    }
    assert_eq!(meshed, 4, "planes and seamed cylinder bands mesh");
    assert_eq!(
        import.report.warnings.len(),
        2,
        "two torus fillets await re-anchoring: {:?}",
        import.report.warnings
    );
}

#[test]
fn every_nist_part_reads_and_reports_honestly() {
    let files = [
        "nist_ctc_01_asme1_rd.stp",
        "nist_ctc_02_asme1_rc.stp",
        "nist_ctc_03_asme1_rc.stp",
        "nist_ctc_04_asme1_rd.stp",
        "nist_ctc_05_asme1_rd.stp",
        "nist_ftc_06_asme1_rd.stp",
        "nist_ftc_07_asme1_rd.stp",
        "nist_ftc_08_asme1_rc.stp",
        "nist_ftc_09_asme1_rd.stp",
        "nist_ftc_10_asme1_rb.stp",
        "nist_ftc_11_asme1_rb.stp",
    ];
    for name in files {
        let text = corpus(name);
        let import =
            og_io::read_step(&text, T).unwrap_or_else(|e| panic!("{name} failed to read: {e}"));
        let mut closed = 0;
        for solid in &import.solids {
            let shell = og_topo::explore_unique(&import.model, solid, og_topo::ShapeType::Shell)
                .unwrap()
                .remove(0);
            if og_algo::is_shell_closed(&import.model, &shell).unwrap() {
                closed += 1;
            }
        }
        eprintln!(
            "REPORT {name}: solids={} closed={} warnings={} skipped_kinds={}",
            import.solids.len(),
            closed,
            import.report.warnings.len(),
            import.report.skipped.len()
        );
        assert!(!import.solids.is_empty(), "{name}: no solid read");
    }
}
