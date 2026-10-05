//! Healing a face the reader refused to trim.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;

const T: Tolerances = Tolerances::millimetres();

const HOVER: &str = r#"ISO-10303-21;
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

/// The reader's refusal is the healer's instruction: a boundary 3 mm off
/// its surface reads untrimmed and refuses to mesh. `fix_face_pcurves` at
/// a caller's cap of 5 mm fits the trims the reader would not, widens the
/// edges to the measured offset, and the face triangulates.
#[test]
fn a_face_the_reader_refused_heals_at_the_callers_cap() {
    let text = HOVER;
    let import = ogeom::io::read_step(text, T).unwrap();
    // The report hands over the face itself. The instructed follow-up
    // needs no search.
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

/// The stronger fix: the hovering boundary *moves* onto its surface, the
/// displacement recorded in widened tolerances, and the trims then fit as
/// on-surface trims: the composition `reanchor_boundaries` then
/// `fix_face_pcurves` that a 3 mm-off face deserves.
#[test]
fn a_hovering_boundary_reanchors_then_trims_then_meshes() {
    let import = ogeom::io::read_step(HOVER, T).unwrap();
    let refused = import.report.untrimmed_faces[0].face.clone();
    let mut document = import.document;
    let solid = import.solids[0].clone();

    let (built, report) =
        ogeom::heal::reanchor_boundaries(document.model_mut(), &solid, 5.0, T).unwrap();
    assert_eq!(report.moved, 4, "all four hovering edges moved");
    assert!(
        (report.worst_before - 3.0).abs() < 1e-6,
        "{}",
        report.worst_before
    );
    assert!(report.worst_after < 1e-3, "{}", report.worst_after);
    assert!(report.refused.is_empty());

    // The face in the rebuilt solid, its boundary now on the surface: the
    // trims fit with next to no offset, and it draws.
    let face = ogeom::topo::explore(
        document.model(),
        &built.shape,
        ogeom::topo::Filter::OfType(ogeom::topo::ShapeType::Face),
    )
    .unwrap()
    .remove(0);
    let trims = ogeom::heal::fix_face_pcurves(document.model_mut(), &face, 1.0, T).unwrap();
    assert_eq!(trims.fitted, 4);
    assert!(trims.worst < 1e-3, "on-surface now: {}", trims.worst);
    let mesh = ogeom::mesh::triangulate_face(
        document.model(),
        &face,
        ogeom::mesh::Deflection::default(),
        T,
    )
    .unwrap();
    assert!(!mesh.triangles.is_empty());
    let _ = refused;
}

/// A boundary with more wiggle than any fixed sample count catches moves
/// onto its plane whole: a 64-span spline weaving `+/-0.05` across its
/// line, hovering 0.01 above the plane. The moved curve follows the
/// wiggle's shadow between the samples it was fitted through, the report
/// states the real residual, and the edge's tolerance covers its pcurves.
#[test]
fn a_wiggling_boundary_reanchors_with_its_wiggle() {
    use ogeom::geom::{Curve, Curve3d as _, PlaneSurface, SurfaceGeometry};
    use ogeom::math::{Frame, KnotVector, Plane, Point};
    use ogeom::topo::{EdgeRepr, Filter, Model, ShapeType, explore};

    const SPANS: usize = 64;
    const HOVER: f64 = 0.01;
    let mut model = Model::new();
    let control: Vec<Point> = (0..SPANS + 3)
        .map(|i| {
            #[allow(clippy::cast_precision_loss)]
            let x = 10.0 * i as f64 / (SPANS + 2) as f64;
            let y = if i == 0 || i == SPANS + 2 {
                0.0
            } else if i.is_multiple_of(2) {
                0.15
            } else {
                -0.15
            };
            Point::new(x, y, HOVER)
        })
        .collect();
    let knots = KnotVector::clamped_uniform(3, control.len()).unwrap();
    let wiggle = Curve::BSpline(ogeom::geom::BSplineCurve::new(knots, control, T).unwrap());
    let corners = [
        Point::new(0.0, 0.0, HOVER),
        Point::new(10.0, 0.0, HOVER),
        Point::new(10.0, 10.0, HOVER),
        Point::new(0.0, 10.0, HOVER),
    ];
    let v: Vec<_> = corners
        .iter()
        .map(|p| ogeom::algo::make_vertex(&mut model, *p).shape)
        .collect();
    let mut edges = vec![
        ogeom::algo::make_edge_between(&mut model, wiggle.clone(), (0.0, 1.0), &v[0], &v[1], T)
            .unwrap()
            .shape,
    ];
    for i in 1..4 {
        let (a, b) = (corners[i], corners[(i + 1) % 4]);
        let line = ogeom::geom::LineCurve::new(ogeom::math::Axis {
            location: a,
            direction: ogeom::math::Direction::new(b - a, T).unwrap(),
        });
        edges.push(
            ogeom::algo::make_edge_between(
                &mut model,
                Curve::from(line),
                (0.0, a.distance(b)),
                &v[i],
                &v[(i + 1) % 4],
                T,
            )
            .unwrap()
            .shape,
        );
    }
    let wire = ogeom::algo::make_wire(&mut model, &edges, T).unwrap().shape;
    let plane: SurfaceGeometry = PlaneSurface::new(Plane::new(Frame::WORLD)).into();
    let face = ogeom::algo::make_face(&mut model, plane, &[wire], T)
        .unwrap()
        .shape;
    let trims = ogeom::heal::fix_face_pcurves(&mut model, &face, 1.0, T).unwrap();
    assert_eq!(trims.fitted, 4);
    let trimmed = ogeom::algo::check(&model, &face, T).unwrap();
    assert!(trimmed.is_valid(), "{:?}", trimmed.problems);

    let (built, report) = ogeom::heal::reanchor_boundaries(&mut model, &face, 1.0, T).unwrap();
    assert_eq!(report.moved, 4);
    let moved = explore(&model, &built.shape, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            let data = model.node(e).unwrap().data().as_edge().unwrap();
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                return false;
            };
            let c = model.geometry().curve(*curve).unwrap();
            let mid = c.point_at(f64::midpoint(range.0, range.1), T).unwrap();
            mid.y.abs() < 1.0 && mid.x > 1.0 && mid.x < 9.0
        })
        .unwrap();
    let data = model.node(&moved).unwrap().data().as_edge().unwrap();
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        unreachable!()
    };
    let curve = model.geometry().curve(*curve).unwrap().clone();
    let (mut off_plane, mut off_shadow) = (0.0_f64, 0.0_f64);
    for i in 0..=4000 {
        let t = range.0 + (range.1 - range.0) * f64::from(i) / 4000.0;
        let p = curve.point_at(t, T).unwrap();
        let w = wiggle.point_at(t, T).unwrap();
        off_plane = off_plane.max(p.z.abs());
        off_shadow = off_shadow.max(p.distance(Point::new(w.x, w.y, 0.0)));
    }
    let stated = data.tolerance.get();
    let gap = ogeom::algo::edge_pcurve_gap(&model, &moved, T)
        .unwrap()
        .unwrap();
    assert!(
        off_shadow <= 1e-4,
        "the moved curve misses the shadow by {off_shadow}"
    );
    assert!(
        report.worst_after >= off_shadow * 0.9,
        "the report says {} for a miss of {off_shadow}",
        report.worst_after
    );
    assert!(gap <= stated, "pcurve gap {gap} over the stated {stated}");
    let diagnosis = ogeom::algo::check(&model, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
}
