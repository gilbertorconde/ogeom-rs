//! Healing an imported part whose torus fillets carry misaligned pcurves.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;

const T: Tolerances = Tolerances::millimetres();

/// The part's two torus fillets are bounded by rims whose vertices stand
/// on different columns. The reader splits one rim of each on the other's
/// column and seams them, so it warns of nothing, and healing finds nothing
/// left to move: the part comes back as itself, every face meshes, the
/// shell closes and the part measures.
#[test]
fn the_smallest_nist_part_reads_seamed_and_heals_to_itself() {
    let path = format!(
        "{}/../../tests/corpus/nist_ftc_11_asme1_rb.stp",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(path).unwrap();
    let mut import = ogeom_io::read_step(&text, T).unwrap();
    assert!(
        import.report.warnings.is_empty(),
        "{:?}",
        import.report.warnings
    );

    let healed =
        ogeom_heal::reanchor_periodic_rings(import.document.model_mut(), &import.solids[0], T)
            .unwrap();
    assert_eq!(healed.1, 0, "nothing is left to heal");
    assert!(healed.0.shape.is_same(&import.solids[0]));

    let fine = ogeom_mesh::Deflection {
        chord: 1e-2,
        ..ogeom_mesh::Deflection::default()
    };
    let model = import.document.model();
    let faces = ogeom_topo::explore(
        model,
        &healed.0.shape,
        ogeom_topo::Filter::OfType(ogeom_topo::ShapeType::Face),
    )
    .unwrap();
    assert_eq!(faces.len(), 6);
    for face in &faces {
        ogeom_mesh::triangulate(model, face, fine, T)
            .unwrap_or_else(|e| panic!("a face fails to mesh: {e}"));
    }
    let shell = ogeom_topo::explore_unique(model, &healed.0.shape, ogeom_topo::ShapeType::Shell)
        .unwrap()
        .remove(0);
    assert!(ogeom_algo::is_shell_closed(model, &shell).unwrap());
    let props = ogeom_algo::volume_properties(model, &healed.0.shape, fine, T).unwrap();
    assert!(props.mass > 0.0, "the part encloses volume");
}

/// A drum whose wall is bounded by two whole circles with no seam, their
/// vertices half a turn apart: healing moves one rim's vertex onto the
/// other's column, rebuilds the cap that shares the rim, and seams the
/// wall. The healed drum closes, checks valid and holds pi r^2 h.
#[test]
fn a_drum_whose_rims_start_half_a_turn_apart_is_reanchored() {
    use ogeom_geom::{CircleCurve, CylinderSurface, PlaneSurface, SurfaceGeometry};
    use ogeom_math::{Circle, Cylinder, Direction, Frame, Plane, Point};

    let mut model = ogeom_topo::Model::new();
    let (r, h) = (2.0, 3.0);
    let at =
        |z: f64, x: Direction| Frame::new(Point::new(0.0, 0.0, z), Direction::Z, x, T).unwrap();
    let mut ring = |frame: Frame| {
        let vertex = ogeom_algo::make_vertex(&mut model, frame.origin() + frame.x() * r).shape;
        ogeom_algo::make_edge_between(
            &mut model,
            CircleCurve::new(Circle::new(frame, r, T).unwrap()).into(),
            (0.0, core::f64::consts::TAU),
            &vertex,
            &vertex,
            T,
        )
        .unwrap()
        .shape
    };
    let low = ring(at(0.0, Direction::X));
    let high = ring(at(h, -Direction::X));
    let wall_surface: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, r, T).unwrap(), (-1.0, h + 1.0))
            .unwrap()
            .into();
    let wall = ogeom_algo::make_face_with_pcurves(
        &mut model,
        wall_surface,
        &[vec![low.clone()], vec![high.reversed()]],
        T,
    )
    .unwrap()
    .shape;
    let down = Frame::new(Point::ORIGIN, -Direction::Z, Direction::X, T).unwrap();
    let bottom = ogeom_algo::make_face_with_pcurves(
        &mut model,
        PlaneSurface::new(Plane::new(down)).into(),
        &[vec![low.reversed()]],
        T,
    )
    .unwrap()
    .shape;
    let top = ogeom_algo::make_face_with_pcurves(
        &mut model,
        PlaneSurface::new(Plane::new(at(h, Direction::X))).into(),
        &[vec![high]],
        T,
    )
    .unwrap()
    .shape;
    let shell = ogeom_algo::make_shell(&mut model, &[bottom, wall, top])
        .unwrap()
        .shape;
    let solid = ogeom_algo::make_solid(&mut model, &[shell]).unwrap().shape;

    let (healed, count) = ogeom_heal::reanchor_periodic_rings(&mut model, &solid, T).unwrap();
    assert!(count > 0, "a rim was moved");
    let diagnosis = ogeom_algo::check(&model, &healed.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let fine = ogeom_mesh::Deflection {
        chord: 1e-3,
        ..ogeom_mesh::Deflection::default()
    };
    let volume = ogeom_algo::volume_properties(&model, &healed.shape, fine, T)
        .unwrap()
        .mass;
    approx::assert_relative_eq!(
        volume,
        core::f64::consts::PI * r * r * h,
        max_relative = 1e-3
    );
}

#[test]
fn the_imported_part_repairs_its_same_parameter_claims() {
    let path = format!(
        "{}/../../tests/corpus/nist_ftc_11_asme1_rb.stp",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(path).unwrap();
    let mut import = ogeom_io::read_step(&text, T).unwrap();
    let solid = import.solids[0].clone();
    let report = ogeom_heal::repair_same_parameter(import.document.model_mut(), &solid, T).unwrap();
    assert!(report.checked > 0);
    // Imported pcurves are fitted against the file's own slop. Some edges
    // widen, and afterwards every claim is true.
    let all_true = ogeom_topo::explore(
        import.document.model(),
        &solid,
        ogeom_topo::Filter::OfType(ogeom_topo::ShapeType::Edge),
    )
    .unwrap()
    .iter()
    .all(|e| {
        import
            .document
            .model()
            .node(e)
            .and_then(|n| n.data().as_edge())
            .is_some_and(ogeom_topo::EdgeData::same_parameter)
    });
    assert!(all_true, "every edge's claim holds after the repair");
}
