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

/// How many edges of a shape's faces, each meshed on its own with the
/// edges drawn alike for every face, are used other than twice: where a
/// face's boundary leaves a gap its neighbours do not fill, with no weld to
/// close it.
fn unmatched_face_edges(model: &Model, shape: &Shape) -> usize {
    let deflection = Deflection::default();
    let chords = ogeom::mesh::edge_chords_for(model, shape, deflection, T).unwrap();
    let key = |p: &ogeom::math::Point| (p.x.to_bits(), p.y.to_bits(), p.z.to_bits());
    let mut uses: std::collections::HashMap<_, usize> = std::collections::HashMap::new();
    for face in explore_unique(model, shape, ShapeType::Face).unwrap() {
        let mesh =
            ogeom::mesh::triangulate_face_with(model, &face, deflection, &chords, T).unwrap();
        for t in &mesh.triangles {
            for k in 0..3 {
                let (a, b) = (
                    key(&mesh.positions[t[k] as usize]),
                    key(&mesh.positions[t[(k + 1) % 3] as usize]),
                );
                *uses.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
    }
    uses.values().filter(|&&n| n != 2).count()
}

/// The corpus part `name` meshed and converted back comes out a valid solid
/// of its volume, tessellating closed where `closed` says it does, its faces
/// meshed one by one meeting edge to edge where `meet` says they do, and
/// with curved faces where the part has them.
fn comes_back(name: &str, closed: bool, meet: bool) {
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
    if meet {
        assert_eq!(
            unmatched_face_edges(&back, &out.shape),
            0,
            "{name}: its faces meshed one by one do not meet"
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
    comes_back("nist_ctc_01_asme1_rd.stp", true, true);
}

#[test]
fn nist_ctc_03() {
    comes_back("nist_ctc_03_asme1_rc.stp", true, true);
}

/// A fitted face here ran past a facet it should have ended on, and the
/// solid came back with the facet facing into material.
#[test]
fn nist_ctc_04() {
    comes_back("nist_ctc_04_asme1_rd.stp", true, true);
}

/// Its whole mesh closes, welded; its faces meshed one by one still leave
/// six edges unmatched.
#[test]
fn nist_ftc_06() {
    comes_back("nist_ftc_06_asme1_rd.stp", true, false);
}

/// As `nist_ctc_04`, on tori.
#[test]
fn nist_ftc_07() {
    comes_back("nist_ftc_07_asme1_rd.stp", true, false);
}

#[test]
fn nist_ftc_08() {
    comes_back("nist_ftc_08_asme1_rc.stp", true, true);
}

#[test]
fn nist_ftc_09() {
    comes_back("nist_ftc_09_asme1_rd.stp", true, true);
}

/// Its seams end a little off the corners they meet at, each on its own
/// side, within the corners' tolerance; drawn from the corners themselves,
/// the faces round each one meet at one point and the mesh closes.
#[test]
fn nist_ftc_10() {
    comes_back("nist_ftc_10_asme1_rb.stp", true, true);
}

#[test]
fn nist_ftc_11() {
    comes_back("nist_ftc_11_asme1_rb.stp", true, true);
}

#[test]
fn a_socket_head_screw() {
    comes_back("m5x16_bhcs_loops.step", true, true);
}
