//! The corpus, finally consumed: reading the NIST test parts.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

#[test]
fn the_smallest_nist_part_reads_into_a_closed_solid() {
    let text = corpus("nist_ftc_11_asme1_rb.stp");
    let import = ogeom_io::read_step(&text, T).unwrap();
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
    let shell =
        ogeom_topo::explore_unique(import.document.model(), solid, ogeom_topo::ShapeType::Shell)
            .unwrap()
            .remove(0);
    assert!(ogeom_algo::is_shell_closed(import.document.model(), &shell).unwrap());

    // Four of the six faces triangulate — planes and full cylinder bands,
    // the bands through synthesised seams. The two torus fillets are the
    // honest remainder: their two ring vertices sit at different angles, so
    // no seam can join them without re-anchoring a shared edge, which is
    // healing's first named import case. The warnings say exactly that.
    let fine = ogeom_mesh::Deflection {
        chord: 1e-2,
        ..ogeom_mesh::Deflection::default()
    };
    let mut meshed = 0;
    for face in ogeom_topo::explore(
        import.document.model(),
        solid,
        ogeom_topo::Filter::OfType(ogeom_topo::ShapeType::Face),
    )
    .unwrap()
    {
        if ogeom_mesh::triangulate(import.document.model(), &face, fine, T).is_ok() {
            meshed += 1;
        }
    }
    // Four planes and seamed cylinder bands, plus the two torus fillets
    // whose wound rings now close against their own translates — every face
    // of the raw import meshes.
    assert_eq!(meshed, 6, "every face meshes, wound torus rings included");
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
            ogeom_io::read_step(&text, T).unwrap_or_else(|e| panic!("{name} failed to read: {e}"));
        let mut closed = 0;
        for solid in &import.solids {
            let shell = ogeom_topo::explore_unique(
                import.document.model(),
                solid,
                ogeom_topo::ShapeType::Shell,
            )
            .unwrap()
            .remove(0);
            if ogeom_algo::is_shell_closed(import.document.model(), &shell).unwrap() {
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
        assert_eq!(
            closed,
            import.solids.len(),
            "{name}: every shell closes as read"
        );
    }
}

/// A cone trimmed to its own apex triangulates.
///
/// The face is bounded the way ST-Developer writes a countersink drilled to
/// a point: the rim, and one slant line used twice — down to the apex and
/// back. No vertex loop, no degenerate edge; the apex exists only as the
/// vertex the slant line ends at. In the chart that line is a seam, and its
/// second traversal must take the other side of the parameter rectangle —
/// continuity cannot say so, because at the apex both sides start at the
/// same 3D point, and choosing by nearness closes the ring over nothing.
/// Found as two invisible countersinks in a real frame assembly (issue #14).
#[test]
fn a_cone_walked_to_its_apex_triangulates() {
    let text = r#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''),'2;1');
FILE_NAME('apex','2026-08-25',(''),(''),'','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=CONICAL_SURFACE('',#4,1.,0.7853981633974483);
#6=PLANE('',#4);
#7=CARTESIAN_POINT('',(0.,0.,-1.));
#8=CARTESIAN_POINT('',(1.,0.,0.));
#9=VERTEX_POINT('',#7);
#10=VERTEX_POINT('',#8);
#11=DIRECTION('',(0.7071067811865476,0.,0.7071067811865476));
#12=VECTOR('',#11,1.);
#13=LINE('',#7,#12);
#14=CIRCLE('',#4,1.);
#15=EDGE_CURVE('',#9,#10,#13,.T.);
#16=EDGE_CURVE('',#10,#10,#14,.T.);
#17=ORIENTED_EDGE('',*,*,#15,.T.);
#18=ORIENTED_EDGE('',*,*,#16,.T.);
#19=ORIENTED_EDGE('',*,*,#15,.F.);
#20=EDGE_LOOP('',(#17,#18,#19));
#21=FACE_OUTER_BOUND('',#20,.T.);
#22=ADVANCED_FACE('',(#21),#5,.F.);
#23=ORIENTED_EDGE('',*,*,#16,.F.);
#24=EDGE_LOOP('',(#23));
#25=FACE_OUTER_BOUND('',#24,.T.);
#26=ADVANCED_FACE('',(#25),#6,.T.);
#27=CLOSED_SHELL('',(#22,#26));
#28=MANIFOLD_SOLID_BREP('',#27);
ENDSEC;
END-ISO-10303-21;
"#;
    let import = ogeom_io::read_step(text, T).unwrap();
    assert_eq!(import.solids.len(), 1);
    let model = import.document.model();
    let solid = &import.solids[0];
    for face in ogeom_topo::explore(
        model,
        solid,
        ogeom_topo::Filter::OfType(ogeom_topo::ShapeType::Face),
    )
    .unwrap()
    {
        let mesh = ogeom_mesh::triangulate_face(model, &face, ogeom_mesh::Deflection::default(), T)
            .expect("every face of the cone-to-apex solid meshes");
        assert!(!mesh.triangles.is_empty());
    }
}

/// A determinate progress bar's contract: reading a file with N solids
/// announces `step: solid` exactly N times, as `(1, N) … (N, N)`.
///
/// The denominator arrives with the *first* event — the host never counts
/// events or guesses the total (issue #9).
#[test]
fn a_step_read_announces_each_solid_with_its_total() {
    use std::sync::{Arc, Mutex};
    let text = corpus("nist_ftc_11_asme1_rb.stp");
    let heard: Arc<Mutex<Vec<(u64, u64)>>> = Arc::new(Mutex::new(Vec::new()));
    let record = Arc::clone(&heard);
    let watch = ogeom_core::progress::Watch::with_stage_sink(move |stage| {
        if stage.name == "step: solid"
            && let Some(at) = stage.progress
        {
            record.lock().unwrap().push(at);
        }
    });
    let import = ogeom_core::progress::watched(&watch, || ogeom_io::read_step(&text, T)).unwrap();
    let n = import.solids.len() as u64;
    let expected: Vec<(u64, u64)> = (1..=n).map(|i| (i, n)).collect();
    assert_eq!(*heard.lock().unwrap(), expected);
}
