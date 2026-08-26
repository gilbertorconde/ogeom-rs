//! Healing a face the reader refused to trim (issue #21).
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;

const T: Tolerances = Tolerances::millimetres();

/// The reader's refusal is the healer's instruction: a boundary 3 mm off
/// its surface reads untrimmed and refuses to mesh; `fix_face_pcurves` at
/// a caller's cap of 5 mm fits the trims the reader would not, widens the
/// edges to the measured offset, and the face triangulates.
#[test]
fn a_face_the_reader_refused_heals_at_the_callers_cap() {
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
    let import = ogeom::io::read_step(text, T).unwrap();
    // The report hands over the face itself — the instructed follow-up
    // needs no search (issue #34).
    assert_eq!(import.report.untrimmed_faces.len(), 1);
    let refused = &import.report.untrimmed_faces[0];
    assert_eq!(refused.entity, 56, "named by the id the warnings use");
    let face = refused.face.clone();
    let mut document = import.document;

    // Before: the face cannot draw.
    assert!(
        ogeom::mesh::triangulate_face(
            document.model(),
            &face,
            ogeom::mesh::Deflection::default(),
            T
        )
        .is_err(),
        "the untrimmed face refuses to mesh"
    );

    // The instructed heal, above the offset, and the face draws.
    let report = ogeom::heal::fix_face_pcurves(document.model_mut(), &face, 5.0, T).unwrap();
    assert_eq!(report.fitted, 4, "all four hovering edges gained trims");
    assert!(report.refused.is_empty());
    assert!(
        (report.worst - 3.0).abs() < 1e-6,
        "the offset is measured: {}",
        report.worst
    );
    let mesh = ogeom::mesh::triangulate_face(
        document.model(),
        &face,
        ogeom::mesh::Deflection::default(),
        T,
    )
    .unwrap();
    assert!(!mesh.triangles.is_empty());

    // A cap *below* the offset still refuses, and says how far.
    let import2 = ogeom::io::read_step(text, T).unwrap();
    let face2 = import2.report.untrimmed_faces[0].face.clone();
    let mut document2 = import2.document;

    let report2 = ogeom::heal::fix_face_pcurves(document2.model_mut(), &face2, 1.0, T).unwrap();
    assert_eq!(report2.fitted, 0);
    assert_eq!(report2.refused.len(), 4);
    assert!(
        report2
            .refused
            .iter()
            .all(|(_, off)| (*off - 3.0).abs() < 0.1)
    );
}
