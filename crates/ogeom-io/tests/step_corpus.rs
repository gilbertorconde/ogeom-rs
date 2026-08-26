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

/// A boundary too far from its surface to trim is *named*, not just lamented.
///
/// The face's edges sit 3 mm above the patch — past the one-millimetre
/// healing cap, the wrong-pairing regime — so the fit refuses, the face
/// reads without a trim, and `report.untrimmed_faces` carries its STEP id
/// for the consumer to mark. The warnings still tell the story in prose;
/// this is the form a UI can act on (issue #15).
#[test]
fn a_face_the_fit_refuses_is_named_in_the_report() {
    let text = r#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''),'2;1');
FILE_NAME('hover','2026-08-26',(''),(''),'','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=CARTESIAN_POINT('',(10.,0.,0.));
#3=CARTESIAN_POINT('',(0.,10.,0.));
#4=CARTESIAN_POINT('',(10.,10.,0.));
#5=B_SPLINE_SURFACE_WITH_KNOTS('',1,1,((#1,#3),(#2,#4)),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),(2,2),(0.,10.),(0.,10.),.UNSPECIFIED.);
#10=CARTESIAN_POINT('',(0.,0.,3.));
#11=CARTESIAN_POINT('',(10.,0.,3.));
#12=CARTESIAN_POINT('',(10.,10.,3.));
#13=CARTESIAN_POINT('',(0.,10.,3.));
#14=VERTEX_POINT('',#10);
#15=VERTEX_POINT('',#11);
#16=VERTEX_POINT('',#12);
#17=VERTEX_POINT('',#13);
#20=DIRECTION('',(1.,0.,0.));
#21=DIRECTION('',(0.,1.,0.));
#22=DIRECTION('',(-1.,0.,0.));
#23=DIRECTION('',(0.,-1.,0.));
#24=VECTOR('',#20,1.);
#25=VECTOR('',#21,1.);
#26=VECTOR('',#22,1.);
#27=VECTOR('',#23,1.);
#30=LINE('',#10,#24);
#31=LINE('',#11,#25);
#32=LINE('',#12,#26);
#33=LINE('',#13,#27);
#40=EDGE_CURVE('',#14,#15,#30,.T.);
#41=EDGE_CURVE('',#15,#16,#31,.T.);
#42=EDGE_CURVE('',#16,#17,#32,.T.);
#43=EDGE_CURVE('',#17,#14,#33,.T.);
#50=ORIENTED_EDGE('',*,*,#40,.T.);
#51=ORIENTED_EDGE('',*,*,#41,.T.);
#52=ORIENTED_EDGE('',*,*,#42,.T.);
#53=ORIENTED_EDGE('',*,*,#43,.T.);
#54=EDGE_LOOP('',(#50,#51,#52,#53));
#55=FACE_OUTER_BOUND('',#54,.T.);
#56=ADVANCED_FACE('',(#55),#5,.T.);
#57=CLOSED_SHELL('',(#56));
#58=MANIFOLD_SOLID_BREP('',#57);
ENDSEC;
END-ISO-10303-21;
"#;
    let import = ogeom_io::read_step(text, T).unwrap();
    assert_eq!(
        import.report.untrimmed_faces.len(),
        1,
        "the hovering face is named exactly once"
    );
    let refused = &import.report.untrimmed_faces[0];
    assert_eq!(refused.entity, 56, "by its file id");
    let of_solid = ogeom_topo::explore_unique(
        import.document.model(),
        &import.solids[0],
        ogeom_topo::ShapeType::Face,
    )
    .unwrap();
    assert!(
        of_solid.iter().any(|f| f.node() == refused.face.node()),
        "and the carried shape is a face of the solid"
    );
    assert!(
        import
            .report
            .warnings
            .iter()
            .any(|w| w.contains("no pcurve")),
        "and the prose still tells the story"
    );
}

/// The warning flood, counted: the hovering face's read carries its kinds
/// as summary entries — count, worst measured value, an exemplar id — so a
/// consumer shows four lines where the prose runs to hundreds (issue #24).
#[test]
fn warnings_summarise_by_kind_with_counts_and_worsts() {
    let text = corpus("nist_ftc_11_asme1_rb.stp");
    let import = ogeom_io::read_step(&text, T).unwrap();
    // The corpus file is imprecise the way real files are: something
    // tallies. Every entry is coherent — counted, and its worst finite.
    for entry in &import.report.summary {
        assert!(entry.count > 0);
        assert!(entry.worst.is_finite());
    }
    // The summary is a digest, not a second flood.
    assert!(
        import.report.summary.len() <= 8,
        "kinds, not occurrences: {}",
        import.report.summary.len()
    );
    let prose = import.report.warnings.len();
    let counted: usize = import.report.summary.iter().map(|e| e.count).sum();
    assert!(
        counted <= prose + import.report.untrimmed_faces.len(),
        "the summary counts what the prose says: {counted} vs {prose}"
    );
}
