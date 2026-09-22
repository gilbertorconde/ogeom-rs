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

#[test]
fn a_slit_sphere_zone_meshes_its_own_region_not_the_complement() {
    // Issue #37: the button head's sphere zone is slit along a meridian
    // that sits at the chart's own seam. The doubly-used slit edge was
    // bracketed like a period-wrapping seam, which wound the ring, pulled
    // in a pole row the face never touches, and meshed the complement of
    // the head — a fan the size of the whole sphere. The honest measure is
    // chart-consistency: every triangle's centre must sit at sag distance
    // from the surface evaluated at its own chart centre.
    use ogeom_geom::{Surface as _, Transformable as _};
    use ogeom_topo::{Filter, NodeData, ShapeType, explore};
    let text = corpus("m5x16_bhcs.step");
    let import = ogeom_io::step::read_step(&text, T).unwrap();
    let model = import.document.model();
    let solid = &import.solids[0];
    let mut checked = 0;
    for face in explore(model, solid, Filter::OfType(ShapeType::Face)).unwrap() {
        let NodeData::Face(d) = model.node(&face).unwrap().data() else {
            continue;
        };
        let surface = model
            .geometry()
            .surface(d.surface)
            .unwrap()
            .clone()
            .transformed(&face.transform(model.datums()).unwrap(), T)
            .unwrap();
        let Ok(mesh) =
            ogeom_mesh::triangulate_face(model, &face, ogeom_mesh::Deflection::default(), T)
        else {
            continue;
        };
        let mut worst = 0.0_f64;
        for t in &mesh.triangles {
            let [a, b, c] = [
                mesh.positions[t[0] as usize],
                mesh.positions[t[1] as usize],
                mesh.positions[t[2] as usize],
            ];
            let mid = ogeom_math::Point::from_vector(
                (a.to_vector() + b.to_vector() + c.to_vector()) / 3.0,
            );
            let mu = (mesh.parameters[t[0] as usize].0
                + mesh.parameters[t[1] as usize].0
                + mesh.parameters[t[2] as usize].0)
                / 3.0;
            let mv = (mesh.parameters[t[0] as usize].1
                + mesh.parameters[t[1] as usize].1
                + mesh.parameters[t[2] as usize].1)
                / 3.0;
            if let Ok(on) = surface.point_at(mu, mv, T) {
                worst = worst.max(on.distance(mid));
            }
        }
        assert!(
            worst < 0.2,
            "a face's triangles wander {worst} mm from their own chart"
        );
        checked += 1;
    }
    assert!(checked > 10, "the screw has faces to check: {checked}");
}

#[test]
fn two_oblique_rims_on_a_sphere_bound_a_face_not_a_band() {
    // Issue #37, the assembly-side half: the button head's sphere zone is
    // bounded by two closed circles cut square to the screw, while the
    // sphere's chart runs along z. They lie on the sphere but are not its
    // parallels, and the reader used to synthesise a band between them —
    // a phantom meridian slit and a latitude-line pcurve per rim that the
    // rims never follow, meshing the sphere's complement. Now they bound
    // the face on their own: two wires, no seam, and every triangle at sag
    // distance from the surface at its own chart centre.
    use ogeom_geom::{Surface as _, Transformable as _};
    use ogeom_topo::{EdgeRepr, Filter, NodeData, ShapeType, explore};
    let text = corpus("m5x16_bhcs_loops.step");
    let import = ogeom_io::step::read_step(&text, T).unwrap();
    assert!(
        !import
            .report
            .warnings
            .iter()
            .any(|w| w.contains("no seam could be synthesised")),
        "the reader must not warn about a face that is whole on its own bounds"
    );
    let model = import.document.model();
    let mut heads = 0;
    for face in explore(model, &import.solids[0], Filter::OfType(ShapeType::Face)).unwrap() {
        let NodeData::Face(d) = model.node(&face).unwrap().data() else {
            continue;
        };
        let Some(stored) = model.geometry().surface(d.surface) else {
            continue;
        };
        if !matches!(stored, ogeom_geom::SurfaceGeometry::Sphere(_)) {
            continue;
        }
        heads += 1;
        assert_eq!(
            model.ordered_children_of(&face).unwrap().len(),
            2,
            "the head keeps its two rims as two wires"
        );
        for e in explore(model, &face, Filter::OfType(ShapeType::Edge)).unwrap() {
            let NodeData::Edge(ed) = model.node(&e).unwrap().data() else {
                continue;
            };
            assert!(
                !matches!(
                    ed.pcurve_for(d.surface, e.location()),
                    Some(EdgeRepr::Seam { .. })
                ),
                "no slit was manufactured"
            );
        }
        let surface = stored
            .clone()
            .transformed(&face.transform(model.datums()).unwrap(), T)
            .unwrap();
        let mesh = ogeom_mesh::triangulate_face(model, &face, ogeom_mesh::Deflection::default(), T)
            .unwrap();
        let mut worst = 0.0_f64;
        for t in &mesh.triangles {
            let [a, b, c] = [
                mesh.positions[t[0] as usize],
                mesh.positions[t[1] as usize],
                mesh.positions[t[2] as usize],
            ];
            let mid = ogeom_math::Point::from_vector(
                (a.to_vector() + b.to_vector() + c.to_vector()) / 3.0,
            );
            let mu = (mesh.parameters[t[0] as usize].0
                + mesh.parameters[t[1] as usize].0
                + mesh.parameters[t[2] as usize].0)
                / 3.0;
            let mv = (mesh.parameters[t[0] as usize].1
                + mesh.parameters[t[1] as usize].1
                + mesh.parameters[t[2] as usize].1)
                / 3.0;
            if let Ok(on) = surface.point_at(mu, mv, T) {
                worst = worst.max(on.distance(mid));
            }
        }
        assert!(
            worst < 0.2,
            "the head's triangles sit on its own zone: {worst} mm"
        );
    }
    assert_eq!(heads, 1, "one button head on the screw");
}

/// A part exported as faces — `SHELL_BASED_SURFACE_MODEL` bodies, one open
/// shell of one face each — reads as shells under its product, every face
/// meshing, nothing of it left in the skipped table.
#[test]
fn a_surface_model_reads_as_shells_under_its_product() {
    let text = corpus("nema17_coupler_faces.step");
    let import = ogeom_io::read_step(&text, T).unwrap();
    let model = import.document.model();
    assert!(import.solids.is_empty(), "the file names no solid");
    assert_eq!(import.shells.len(), 6, "six surface bodies");
    for shell in &import.shells {
        assert_eq!(model.kind_of(shell).unwrap(), ogeom_topo::ShapeType::Shell);
        assert!(!ogeom_algo::is_shell_closed(model, shell).unwrap());
    }
    for keyword in ["SHELL_BASED_SURFACE_MODEL", "OPEN_SHELL", "ADVANCED_FACE"] {
        assert!(
            !import.report.skipped.contains_key(keyword),
            "{keyword} was read"
        );
    }

    let (_, coupler) = import
        .document
        .products()
        .find(|(_, p)| p.name == "NEMA17_Coupler")
        .expect("the product the file names");
    let ogeom_doc::ProductKind::Part { shape } = &coupler.kind else {
        panic!("a part");
    };
    let faces = ogeom_topo::explore(
        model,
        shape,
        ogeom_topo::Filter::OfType(ogeom_topo::ShapeType::Face),
    )
    .unwrap();
    assert_eq!(faces.len(), 6, "every body's face is under the product");
    for face in &faces {
        ogeom_mesh::triangulate(model, face, ogeom_mesh::Deflection::default(), T)
            .expect("a surface body's face meshes like any other");
    }
}

/// A hole loop straddling a drum's seam: its two edges' pcurves arrive on
/// different branches of the chart, and read as they come the loop never
/// closes and the hole is drawn but not cut. Chained onto one branch it
/// cuts: no mesh vertex lands inside the hole.
#[test]
fn a_hole_across_the_chart_seam_is_cut_from_the_face() {
    use ogeom_geom::Curve3d as _;
    let text = corpus("nema17_coupler_hole.step");
    let import = ogeom_io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let shell = &import.shells[0];
    let face = ogeom_topo::explore(
        model,
        shell,
        ogeom_topo::Filter::OfType(ogeom_topo::ShapeType::Face),
    )
    .unwrap()
    .remove(0);
    let wires = model.ordered_children_of(&face).unwrap();
    assert_eq!(wires.len(), 2, "the wall and its hole");
    // The hole's centre and reach, from its own edges.
    let mut points = Vec::new();
    for edge in model.ordered_children_of(&wires[1]).unwrap() {
        let data = model.node(&edge).unwrap().data().as_edge().unwrap().clone();
        let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            panic!("a hole edge has a curve");
        };
        let curve = model.geometry().curve(*curve).unwrap();
        for k in 0..=16 {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 16.0;
            points.push(curve.point_at(t, T).unwrap());
        }
    }
    let mut sum = ogeom_math::Vector::ZERO;
    for p in &points {
        sum += p.to_vector();
    }
    #[allow(clippy::cast_precision_loss)]
    let centre = ogeom_math::Point::ORIGIN + sum * (1.0 / points.len() as f64);
    let reach = points
        .iter()
        .map(|p| p.distance(centre))
        .fold(f64::INFINITY, f64::min);
    assert!(reach > 1.0, "a real hole: {reach}");
    let mesh = ogeom_mesh::triangulate(model, &face, ogeom_mesh::Deflection::default(), T).unwrap();
    let inside = mesh
        .positions
        .iter()
        .filter(|v| v.distance(centre) < reach * 0.7)
        .count();
    assert_eq!(inside, 0, "the hole is cut: no mesh vertex inside it");
}

/// A surface whose own placement sits half a kilometre from the face it
/// carries still gets a window wide enough to be asked about.
///
/// A plane, a cylinder and a cone are unbounded, so the window a reader
/// gives them is a convention — generous, and a guess. A real assembly
/// falsifies the guess: a community printer assembly places a cylinder's
/// origin at `z = 500000` and trims the face it carries near the world
/// origin, so the trim's height parameter runs to −5e5 where the window
/// stopped at −1e5. The surface then refused to be evaluated where its own
/// face lies, and fifteen faces of that assembly drew as holes.
///
/// The window is measured from the edges that bound the surface now, with
/// the convention kept as a floor.
#[test]
fn a_surface_placed_far_from_its_face_still_meshes() {
    let text = r#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''),'2;1');
FILE_NAME('far','2026-09-17',(''),(''),'','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,500000.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=CYLINDRICAL_SURFACE('',#4,3.5);
#6=CARTESIAN_POINT('',(0.,0.,0.));
#7=AXIS2_PLACEMENT_3D('',#6,#2,#3);
#8=PLANE('',#7);
#9=CARTESIAN_POINT('',(0.,0.,10.));
#10=AXIS2_PLACEMENT_3D('',#9,#2,#3);
#11=PLANE('',#10);
#12=CIRCLE('',#7,3.5);
#13=CIRCLE('',#10,3.5);
#14=CARTESIAN_POINT('',(3.5,0.,0.));
#15=VERTEX_POINT('',#14);
#16=CARTESIAN_POINT('',(3.5,0.,10.));
#17=VERTEX_POINT('',#16);
#18=EDGE_CURVE('',#15,#15,#12,.T.);
#19=EDGE_CURVE('',#17,#17,#13,.T.);
#20=ORIENTED_EDGE('',*,*,#18,.T.);
#21=EDGE_LOOP('',(#20));
#22=FACE_OUTER_BOUND('',#21,.T.);
#23=ORIENTED_EDGE('',*,*,#19,.T.);
#24=EDGE_LOOP('',(#23));
#25=FACE_BOUND('',#24,.T.);
#26=ADVANCED_FACE('',(#22,#25),#5,.T.);
#27=ORIENTED_EDGE('',*,*,#18,.F.);
#28=EDGE_LOOP('',(#27));
#29=FACE_OUTER_BOUND('',#28,.T.);
#30=ADVANCED_FACE('',(#29),#8,.F.);
#31=ORIENTED_EDGE('',*,*,#19,.F.);
#32=EDGE_LOOP('',(#31));
#33=FACE_OUTER_BOUND('',#32,.T.);
#34=ADVANCED_FACE('',(#33),#11,.T.);
#35=CLOSED_SHELL('',(#26,#30,#34));
#36=MANIFOLD_SOLID_BREP('',#35);
ENDSEC;
END-ISO-10303-21;
"#;
    let import = ogeom_io::read_step(text, T).unwrap();
    assert_eq!(import.solids.len(), 1);
    let model = import.document.model();
    let faces = ogeom_topo::explore(
        model,
        &import.solids[0],
        ogeom_topo::Filter::OfType(ogeom_topo::ShapeType::Face),
    )
    .unwrap();
    assert_eq!(faces.len(), 3, "a wall and two caps");
    for face in &faces {
        let mesh = ogeom_mesh::triangulate_face(model, face, ogeom_mesh::Deflection::default(), T)
            .expect("every face of a far-placed surface meshes");
        assert!(!mesh.triangles.is_empty());
    }
}

/// A trim on a sliver patch stays in its own chart.
///
/// The fixture is a degree 3×3 patch whose `u` direction is degenerate the
/// whole way across: `du` is exactly zero along `v = 0` and
/// four ten-thousandths at the far edge — a sliver four microns wide and a
/// tenth of a millimetre long. The projector answered `u = 0` at the first
/// sample and `u = 1` at every other, the fit swung across the chart to
/// join them, and its control points were dragged back into the window by
/// seven hundred and sixty-seven chart units, which the reader reported as
/// millimetres of mesh error on a face a tenth of a millimetre across.
///
/// Every `u` on such a patch describes the same points to within the
/// sliver's own width, so the trace takes one of them and the fit stays put.
#[test]
fn a_sliver_patch_s_trim_stays_in_its_chart() {
    let text = corpus("spline_face_fit_runs_away.step");
    let import = ogeom_io::read_step(&text, T).unwrap();
    for w in &import.report.warnings {
        eprintln!("REPORT warn: {w}");
    }
    let short: Vec<&String> = import
        .report
        .warnings
        .iter()
        .filter(|w| w.contains("pcurve fit stopped"))
        .collect();
    assert!(short.is_empty(), "the trims fit: {short:?}");

    // And the face is still there, drawn where it belongs.
    let model = import.document.model();
    let faces =
        ogeom_topo::explore_unique(model, &import.solids[0], ogeom_topo::ShapeType::Face).unwrap();
    let mesh =
        ogeom_mesh::triangulate(model, &faces[0], ogeom_mesh::Deflection::default(), T).unwrap();
    assert!(!mesh.triangles.is_empty(), "the sliver meshes");
    let bounds = ogeom_algo::shape_bounds(model, &import.solids[0], T).unwrap();
    assert!(
        bounds.diagonal() < 1.0,
        "and stays a tenth of a millimetre across: {bounds:?}"
    );
}

/// A solid with a cavity is read, cavity and all.
///
/// `BREP_WITH_VOIDS` is a subtype of `MANIFOLD_SOLID_BREP` — same name and
/// outer shell in the same two places, plus the shells that bound its
/// voids — so a file writes it under its own keyword and a reader matching
/// on the leading keyword alone never sees it. Three bodies of a community
/// printer assembly are written that way — a printed housing with six
/// cavities among them — and all three were simply absent from the import.
///
/// The cavity arrives as an `ORIENTED_CLOSED_SHELL` pointing the other way
/// round, which is what makes the volume come out as the difference rather
/// than the sum.
#[test]
fn a_solid_with_a_cavity_keeps_its_cavity() {
    let text = corpus("box_with_a_cavity.step");
    let import = ogeom_io::read_step(&text, T).unwrap();
    assert_eq!(import.solids.len(), 1, "the voided solid is a solid");

    let model = import.document.model();
    let solid = &import.solids[0];
    let shells = ogeom_topo::explore_unique(model, solid, ogeom_topo::ShapeType::Shell).unwrap();
    assert_eq!(shells.len(), 2, "the hull and the cavity");
    let faces = ogeom_topo::explore_unique(model, solid, ogeom_topo::ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 12, "six walls each");

    // A 10 mm cube with a 4 mm cube taken out of its middle: the cavity is
    // subtracted, not added, and the centroid stays at the centre.
    let props =
        ogeom_algo::volume_properties(model, solid, ogeom_mesh::Deflection::default(), T).unwrap();
    assert!(
        (props.mass - 936.0).abs() < 1e-6,
        "1000 less 64: {}",
        props.mass
    );
    for axis in [props.centre.x, props.centre.y, props.centre.z] {
        assert!((axis - 5.0).abs() < 1e-6, "centred: {axis}");
    }
}

/// A closed edge's seam is moved to its vertex, not the vertex's tolerance
/// to the seam.
///
/// A fitted loop written with its start wherever the fit began, the edge's
/// one vertex 2.18 mm along it. Held to the curve's own seam, the vertex
/// missed it by that much and its tolerance was widened to 2.18 mm — a
/// reach a solid's border weld then used. The seam is moved to the vertex:
/// the same curve, begun where the edge does, and the vertex met exactly.
#[test]
fn a_closed_edge_s_seam_is_moved_to_its_vertex() {
    let text = corpus("closed_edge_vertex_off_the_seam.step");
    let import = ogeom_io::read_step(&text, T).unwrap();
    assert!(
        import
            .report
            .warnings
            .iter()
            .any(|w| w.contains("the seam was moved to the vertex")),
        "the reseam is reported: {:?}",
        import.report.warnings
    );
    let worst_miss = import
        .report
        .warnings
        .iter()
        .filter(|w| w.contains("misses its vertex by"))
        .filter_map(|w| {
            w.split("misses its vertex by ")
                .nth(1)?
                .split(';')
                .next()?
                .parse::<f64>()
                .ok()
        })
        .fold(0.0_f64, f64::max);
    assert!(
        worst_miss < 1e-3,
        "no vertex is missed by millimetres any more: worst {worst_miss:.2e}"
    );
}

/// A product named the plainer way still reads under its name.
///
/// A modeller writing an assembly, and a mesh converter writing a shell,
/// spell the formation with its source — `..._WITH_SPECIFIED_SOURCE`, the
/// product in the same third slot — and the converter leaves the product's
/// name blank and fills its id instead. The reader took the first slot of
/// the source-spelt formation, which is its blank id, and named every such
/// product after its definition's entity number.
#[test]
fn a_product_spelt_with_its_source_and_a_blank_name_reads_by_name() {
    let text = corpus("ogeom_asm_bolted_plate.stp")
        .replace(
            "PRODUCT_DEFINITION_FORMATION('','',#14);",
            "PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE('','',#14,.NOT_KNOWN.);",
        )
        .replace(
            "PRODUCT_DEFINITION_FORMATION('','',#158);",
            "PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE('','',#158,.NOT_KNOWN.);",
        )
        .replace(
            "PRODUCT('bolt','bolt','',(#8));",
            "PRODUCT('bolt','','',(#8));",
        );
    let import = ogeom_io::read_step(&text, T).unwrap();
    let names: Vec<String> = import
        .document
        .products()
        .map(|(_, p)| p.name.clone())
        .collect();
    assert_eq!(names, ["plate", "bolt", "bolted-plate"]);
}
