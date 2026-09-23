//! §17's mesh and drawing half: what this crate writes, it reads back.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Frame, Point, Point2};
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

fn corpus_bytes(name: &str) -> Vec<u8> {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read(path).expect("the corpus file is committed")
}

/// Six times the volume a closed mesh encloses, signed: positive when
/// every triangle winds counter-clockwise about its outward normal.
fn signed_volume(mesh: &ogeom::topo::Triangulation) -> f64 {
    mesh.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize] - Point::ORIGIN);
            a.dot(b.cross(c)) / 6.0
        })
        .sum()
}

/// A package as a slicer writes it — deflated entries streamed with their
/// sizes after the data, the start part named by the relationships — reads
/// to one cube, placed by its build item and coloured by its material.
#[test]
fn a_deflated_3mf_package_reads() {
    let bytes = corpus_bytes("threemf_cube_deflated.3mf");
    let import = ogeom::io::read_3mf(&bytes, T).unwrap();
    assert!(import.warnings.is_empty(), "{:?}", import.warnings);
    assert_eq!(import.objects.len(), 1);
    let cube = &import.objects[0];
    assert_eq!(cube.name.as_deref(), Some("cube"));
    assert_eq!(cube.object_type, ogeom::io::ObjectType::Model);
    assert_eq!(cube.colour, Some([1.0, 0.0, 0.0, 1.0]));
    assert_eq!(cube.mesh.triangles.len(), 12);
    assert_eq!(cube.mesh.positions.len(), 8);
    let min_x = cube
        .mesh
        .positions
        .iter()
        .map(|p| p.x)
        .fold(f64::MAX, f64::min);
    assert!(
        (min_x - 10.0).abs() < 1e-9,
        "the build transform is applied"
    );
    assert!((signed_volume(&cube.mesh) - 1000.0).abs() < 1e-9);

    // A damaged byte in the model part fails its checksum, or its
    // inflation, and is refused rather than read wrong.
    let header = bytes
        .windows(13)
        .position(|w| w == b"3D/cube.model")
        .unwrap();
    let mut damaged = bytes.clone();
    damaged[header + 60] ^= 0x55;
    assert!(ogeom::io::read_3mf(&damaged, T).is_err());
    assert!(ogeom::io::read_3mf(&bytes[..bytes.len() / 2], T).is_err());
}

#[test]
fn write_3mf_round_trips_through_read_3mf() {
    let mesh = block_mesh();
    let bytes = ogeom::io::write_3mf(&[ogeom::io::threemf::Object {
        mesh: &mesh,
        name: Some("block".into()),
    }]);
    let back = ogeom::io::read_3mf(&bytes, T).unwrap();
    assert_eq!(back.objects.len(), 1);
    assert_eq!(back.objects[0].name.as_deref(), Some("block"));
    assert_eq!(back.objects[0].mesh.triangles.len(), mesh.triangles.len());
    assert!((signed_volume(&back.objects[0].mesh) - signed_volume(&mesh)).abs() < 1e-9);
}

/// Components flatten into their object, from a model part of their own
/// reached through the production extension's path; the model's inch
/// scales to millimetres; and the component placed by a mirroring
/// transform is rewound, so the pair still encloses its volume outward.
#[test]
fn components_and_units_flatten_to_millimetres() {
    let import = ogeom::io::read_3mf(&corpus_bytes("threemf_components_inch.3mf"), T).unwrap();
    assert!(import.warnings.is_empty(), "{:?}", import.warnings);
    let pair = &import.objects[0];
    assert_eq!(pair.name.as_deref(), Some("pair"));
    assert_eq!(pair.mesh.triangles.len(), 24);
    let (lo, hi) =
        pair.mesh
            .positions
            .iter()
            .fold(([f64::MAX; 3], [f64::MIN; 3]), |(lo, hi), p| {
                (
                    [lo[0].min(p.x), lo[1].min(p.y), lo[2].min(p.z)],
                    [hi[0].max(p.x), hi[1].max(p.y), hi[2].max(p.z)],
                )
            });
    assert!(lo.iter().all(|c| c.abs() < 1e-9));
    assert!((hi[0] - 3.0 * 25.4).abs() < 1e-9 && (hi[2] - 25.4).abs() < 1e-9);
    let inch3 = 25.4_f64.powi(3);
    assert!((signed_volume(&pair.mesh) - 2.0 * inch3).abs() < 1e-6 * inch3);
}

/// A mesh painted in more than one colour keeps its geometry and drops the
/// colours with a warning; a support structure is read, and said to be one.
#[test]
fn per_triangle_colours_and_supports_are_read_with_a_warning() {
    let import = ogeom::io::read_3mf(&corpus_bytes("threemf_painted_sphere.3mf"), T).unwrap();
    assert_eq!(import.objects.len(), 2);
    let (sphere, prop) = (&import.objects[0], &import.objects[1]);
    assert_eq!(sphere.mesh.triangles.len(), 4680);
    assert_eq!(
        sphere.colour,
        Some([224.0 / 255.0, 224.0 / 255.0, 224.0 / 255.0, 1.0])
    );
    let full = 4.0 / 3.0 * std::f64::consts::PI * 20.0_f64.powi(3);
    let v = signed_volume(&sphere.mesh);
    assert!(v > 0.98 * full && v < full, "{v} against {full}");
    assert_eq!(prop.object_type, ogeom::io::ObjectType::Support);
    assert!(
        import
            .warnings
            .iter()
            .any(|w| w.contains("several colours"))
    );
    assert!(import.warnings.iter().any(|w| w.contains("support")));
}

/// A package as a streaming writer lays it out whatever its size — ZIP64
/// end record and locator, every size and offset saturated in the headers
/// and carried in ZIP64 extra fields — reads to the same cube as the
/// classic form of the same package, whose entries stream their sizes in
/// data descriptors instead.
#[test]
fn zip64_and_classic_forms_of_one_package_read_the_same() {
    let classic = ogeom::io::read_3mf(&corpus_bytes("threemf_cube_deflated.3mf"), T).unwrap();
    let zip64 = ogeom::io::read_3mf(&corpus_bytes("threemf_cube_zip64.3mf"), T).unwrap();
    assert_eq!(zip64.objects.len(), 1);
    assert_eq!(zip64.objects[0].mesh.triangles.len(), 12);
    assert_eq!(
        classic.objects[0].mesh.triangles,
        zip64.objects[0].mesh.triangles
    );
    assert_eq!(
        classic.objects[0].mesh.positions,
        zip64.objects[0].mesh.positions
    );
    assert_eq!(classic.objects[0].colour, zip64.objects[0].colour);
}
