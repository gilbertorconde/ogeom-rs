//! §17's mesh and drawing half: what this crate writes, it reads back.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Frame, Point2};
use ogeom::mesh::Deflection;
use ogeom::topo::Model;

const T: Tolerances = Tolerances::millimetres();

fn block_mesh() -> ogeom::topo::Triangulation {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 20.0, 30.0), T)
        .unwrap()
        .shape;
    ogeom::mesh::triangulate(&model, &block, Deflection::default(), T).unwrap()
}

#[test]
fn obj_and_ply_come_back_as_the_mesh_they_were() {
    let mesh = block_mesh();
    let written =
        ogeom::io::mesh_formats::write_obj(&[ogeom::io::mesh_formats::ExportMesh::plain(&mesh)]);
    let read = ogeom::io::mesh_formats::read_obj(&written).unwrap();
    assert_eq!(read.positions.len(), mesh.positions.len());
    assert_eq!(read.triangles.len(), mesh.triangles.len());
    for (a, b) in read.positions.iter().zip(&mesh.positions) {
        assert!(a.distance(*b) < 1e-6, "{a:?} against {b:?}");
    }
    // The box's own volume survives the trip, which is the only claim a
    // mesh format really makes.
    let volume = |m: &ogeom::topo::Triangulation| -> f64 {
        m.triangles
            .iter()
            .map(|[a, b, c]| {
                let (a, b, c) = (
                    m.positions[*a as usize],
                    m.positions[*b as usize],
                    m.positions[*c as usize],
                );
                a.to_vector().dot(b.to_vector().cross(c.to_vector())) / 6.0
            })
            .sum::<f64>()
            .abs()
    };
    assert!((volume(&read) - 6000.0).abs() < 1e-6, "{}", volume(&read));

    let written =
        ogeom::io::mesh_formats::write_ply(&ogeom::io::mesh_formats::ExportMesh::plain(&mesh));
    let read = ogeom::io::mesh_formats::read_ply(&written).unwrap();
    assert_eq!(read.positions.len(), mesh.positions.len());
    assert!((volume(&read) - 6000.0).abs() < 1e-6);
    // A binary PLY is refused by name rather than mis-parsed.
    let binary = written.replace("format ascii 1.0", "format binary_little_endian 1.0");
    assert!(ogeom::io::mesh_formats::read_ply(&binary).is_err());
}

#[test]
fn a_drawing_written_as_dxf_reads_back_layer_by_layer() {
    let visible = vec![
        vec![
            Point2::new(0.0, 0.0),
            Point2::new(10.0, 0.0),
            Point2::new(10.0, 5.0),
        ],
        vec![Point2::new(0.0, 0.0), Point2::new(0.0, 5.0)],
    ];
    let hidden = vec![vec![Point2::new(2.0, 2.0), Point2::new(8.0, 2.0)]];

    let text = ogeom::io::dxf::write_dxf(&visible, &hidden);
    let read = ogeom::io::dxf::read_dxf(&text).unwrap();
    assert_eq!(read.visible.len(), 2, "{:?}", read.visible);
    assert_eq!(read.hidden.len(), 1, "{:?}", read.hidden);
    for (a, b) in read.visible.iter().zip(&visible) {
        assert_eq!(a.len(), b.len());
        for (p, q) in a.iter().zip(b) {
            assert!(p.distance(*q) < 1e-9, "{p:?} against {q:?}");
        }
    }
    assert!(read.hidden[0][1].distance(Point2::new(8.0, 2.0)) < 1e-9);
}

#[test]
fn a_3mf_is_an_archive_this_crate_can_open_again() {
    let mesh = block_mesh();
    let bytes = ogeom::io::threemf::write_3mf(&[ogeom::io::threemf::Object {
        mesh: &mesh,
        name: Some("block".into()),
    }]);
    // It is a ZIP: the local header's own signature opens it.
    assert_eq!(&bytes[0..4], &[0x50, 0x4b, 0x03, 0x04]);

    let parts = ogeom::io::threemf::read_package(&bytes).unwrap();
    let names: Vec<&str> = parts.iter().map(|(n, _)| n.as_str()).collect();
    assert!(names.contains(&"[Content_Types].xml"));
    assert!(names.contains(&"_rels/.rels"));
    assert!(names.contains(&"3D/3dmodel.model"));

    let (_, model) = parts
        .iter()
        .find(|(n, _)| n == "3D/3dmodel.model")
        .expect("the model part");
    let text = String::from_utf8(model.clone()).unwrap();
    assert!(text.contains("unit=\"millimeter\""));
    assert!(text.contains("name=\"block\""));
    assert_eq!(
        text.matches("<vertex ").count(),
        mesh.positions.len(),
        "every vertex is in the package"
    );
    assert_eq!(text.matches("<triangle ").count(), mesh.triangles.len());
}

#[test]
fn vrml_says_what_it_holds() {
    let mesh = block_mesh();
    let text = ogeom::io::mesh_formats::write_vrml(&[ogeom::io::mesh_formats::ExportMesh {
        mesh: &mesh,
        colour: Some([1.0, 0.5, 0.0, 1.0]),
        name: Some("block".into()),
    }]);
    assert!(text.starts_with("#VRML V2.0 utf8"));
    assert!(text.contains("diffuseColor 1 0.5 0"));
    assert_eq!(text.matches("-1,").count(), mesh.triangles.len());
}
