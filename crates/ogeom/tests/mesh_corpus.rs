//! The STEP corpus meshed and converted back: meshes whose answer is known,
//! the part they were drawn from. Each part is meshed at a thousandth of its
//! diagonal and rebuilt with the default options.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{MeshSolidOptions, check, shape_bounds, solid_from_mesh, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::SurfaceGeometry;
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// How many of a shape's faces are on a surface other than a plane.
fn curved(model: &Model, shape: &Shape) -> usize {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .iter()
        .filter(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            !matches!(
                model.geometry().surface(data.surface),
                Some(SurfaceGeometry::Plane(_))
            )
        })
        .count()
}

/// The corpus part `name` meshed and converted back comes out a valid solid
/// of its volume, tessellating closed where `closed` says it does, and with
/// curved faces where the part has them.
fn comes_back(name: &str, closed: bool) {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(path).expect("the corpus file is committed");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let part = &import.solids[0];
    let diagonal = shape_bounds(model, part, T).unwrap().diagonal();
    let mesh = ogeom::mesh::triangulate(
        model,
        part,
        Deflection::with_chord(diagonal * 1e-3).unwrap(),
        T,
    )
    .unwrap();
    assert!(mesh.is_closed(), "{name}: the source mesh is open");
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(out.closed, "{name}: {:?}", out.report);
    let diagnosis = check(&back, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
    let fine = Deflection::with_chord(diagonal * 1e-5).unwrap();
    let (want, got) = (
        volume_properties(model, part, fine, T).unwrap().mass,
        volume_properties(&back, &out.shape, fine, T).unwrap().mass,
    );
    assert!(
        (got - want).abs() <= want * 3e-3,
        "{name}: {want} drawn, {got} came back"
    );
    if closed {
        let drawn = ogeom::mesh::triangulate(&back, &out.shape, Deflection::default(), T).unwrap();
        assert!(
            drawn.is_closed(),
            "{name}: the converted part tessellates open"
        );
    }
    if curved(model, part) > 0 {
        assert!(
            curved(&back, &out.shape) > 0,
            "{name}: nothing curved came back"
        );
    }
}

#[test]
fn nist_ctc_01() {
    comes_back("nist_ctc_01_asme1_rd.stp", true);
}

#[test]
fn nist_ctc_03() {
    comes_back("nist_ctc_03_asme1_rc.stp", true);
}

/// A fitted face here ran past a facet it should have ended on, and the
/// solid came back with the facet facing into material.
#[test]
fn nist_ctc_04() {
    comes_back("nist_ctc_04_asme1_rd.stp", true);
}

#[test]
fn nist_ftc_06() {
    comes_back("nist_ftc_06_asme1_rd.stp", true);
}

/// As `nist_ctc_04`, on tori.
#[test]
fn nist_ftc_07() {
    comes_back("nist_ftc_07_asme1_rd.stp", false);
}

#[test]
fn nist_ftc_08() {
    comes_back("nist_ftc_08_asme1_rc.stp", true);
}

#[test]
fn nist_ftc_09() {
    comes_back("nist_ftc_09_asme1_rd.stp", true);
}

#[test]
fn nist_ftc_11() {
    comes_back("nist_ftc_11_asme1_rb.stp", true);
}

#[test]
fn a_socket_head_screw() {
    comes_back("m5x16_bhcs_loops.step", true);
}
