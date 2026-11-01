//! Face directions and axes read off the tessellations of corpus parts,
//! against the parts' own surfaces: a planar face tessellates exactly, so
//! each face direction is the area-weighted mean of the normals of the
//! planar faces within the resolution of it, and the area is theirs; a
//! screw turns about its cylinders' axis.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::orient::{Directions, axis_of_revolution, face_directions};
use ogeom::algo::{face_normal, shape_bounds, surface_properties};
use ogeom::core::Tolerances;
use ogeom::geom::SurfaceGeometry;
use ogeom::math::Vector;
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, Triangulation, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The corpus part `name`, read, and tessellated at a thousandth of its
/// diagonal.
fn tessellated(name: &str) -> (Model, Shape, Triangulation) {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(path).expect("the corpus file is committed");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let part = import.solids[0].clone();
    let model = import.document.model().clone();
    let diagonal = shape_bounds(&model, &part, T).unwrap().diagonal();
    let mesh = ogeom::mesh::triangulate(
        &model,
        &part,
        Deflection::with_chord(diagonal * 1e-3).unwrap(),
        T,
    )
    .unwrap();
    assert!(mesh.is_closed(), "{name}: the tessellation is open");
    (model, part, mesh)
}

/// Each planar face's outward normal and area.
fn planes(model: &Model, part: &Shape) -> Vec<(Vector, f64)> {
    explore_unique(model, part, ShapeType::Face)
        .unwrap()
        .iter()
        .filter(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            matches!(
                model.geometry().surface(data.surface),
                Some(SurfaceGeometry::Plane(_))
            )
        })
        .map(|f| {
            let (_, n) = face_normal(model, f, T).unwrap();
            let area = surface_properties(model, f, Deflection::default(), T)
                .unwrap()
                .mass;
            (n, area)
        })
        .collect()
}

fn degrees_apart(a: Vector, b: Vector) -> f64 {
    a.cross(b).magnitude().atan2(a.dot(b).abs()).to_degrees()
}

/// The part's three largest face directions are each the area-weighted
/// mean of the planar faces within the resolution of it, within 0.05°, and
/// hold their area to within half a percent of the part's: a triangle of a
/// curved face whose normal passes within the resolution counts too, and a
/// planar face between two directions counts with the nearer.
fn faces_the_mean_of_its_planes(name: &str) {
    let (model, part, mesh) = tessellated(name);
    let options = Directions::default();
    let total = mesh.area();
    let planes = planes(&model, &part);
    let found = face_directions(&mesh, &options, T).unwrap();
    assert!(found.len() >= 3, "{name}: {found:?}");
    for f in &found[..3] {
        let d = f.direction.vector();
        let (sum, area) = planes
            .iter()
            .filter(|(n, _)| n.dot(d).abs() >= options.resolution.cos())
            .fold((Vector::ZERO, 0.0), |(s, a), (n, w)| {
                let n = if n.dot(d) < 0.0 { -*n } else { *n };
                (s + n * *w, a + w)
            });
        let apart = degrees_apart(d, sum / sum.magnitude());
        assert!(
            apart <= 0.05,
            "{name}: {apart}° from its planes' mean: {f:?}"
        );
        assert!(
            (f.area - area).abs() <= 5e-3 * total,
            "{name}: {} of planes holding {area}, of {total}",
            f.area
        );
    }
}

#[test]
fn a_prismatic_part_faces_its_planes() {
    faces_the_mean_of_its_planes("nist_ctc_01_asme1_rd.stp");
}

#[test]
fn a_drafted_part_faces_the_mean_of_its_drafted_planes() {
    // Walls drafted 2° and 3° either side of the axes they face.
    faces_the_mean_of_its_planes("nist_ctc_02_asme1_rc.stp");
}

#[test]
fn a_part_with_a_face_at_45_degrees_faces_its_planes() {
    faces_the_mean_of_its_planes("nist_ctc_03_asme1_rc.stp");
}

#[test]
fn a_screw_turns_about_its_cylinders_axis() {
    let name = "m5x16_bhcs_loops.step";
    let (model, part, mesh) = tessellated(name);
    let cylinders: Vec<_> = explore_unique(&model, &part, ShapeType::Face)
        .unwrap()
        .iter()
        .filter_map(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            match model.geometry().surface(data.surface) {
                Some(SurfaceGeometry::Cylinder(c)) => {
                    let placed = f.transform(model.datums()).unwrap();
                    let frame = c.cylinder().frame();
                    Some((
                        placed.apply(frame.origin()),
                        placed.apply_vector(frame.z().vector()),
                    ))
                }
                _ => None,
            }
        })
        .collect();
    assert!(!cylinders.is_empty());
    let found = axis_of_revolution(&mesh, &Directions::default(), T)
        .unwrap()
        .unwrap();
    for (through, along) in cylinders {
        let apart = degrees_apart(found.axis.direction.vector(), along);
        let off = (found.axis.location - through).cross(along).magnitude();
        assert!(apart <= 0.05, "{name}: {apart}° off the cylinder's axis");
        assert!(off <= 0.01, "{name}: {off} from the cylinder's axis");
    }
    assert!(found.share >= 0.5, "{name}: {found:?}");
}
