//! Reading VRML: what the writer wrote comes back, and hand-written scenes
//! in both dialects draw where their transforms put them.
#![allow(clippy::unwrap_used, reason = "test code")]

use ogeom::io::{ImportedMesh, read_vrml};
use ogeom::math::Point;

/// The volume a closed mesh encloses, by the divergence theorem.
fn volume(mesh: &ImportedMesh) -> f64 {
    mesh.mesh
        .triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.mesh.positions[i as usize].to_vector());
            a.dot(b.cross(c)) / 6.0
        })
        .sum()
}

fn bounds(mesh: &ImportedMesh) -> (Point, Point) {
    let mut lo = Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
    for p in &mesh.mesh.positions {
        lo = Point::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
        hi = Point::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
    }
    (lo, hi)
}

#[test]
fn what_the_writer_writes_reads_back() {
    let mut model = ogeom::topo::Model::new();
    let tol = ogeom::core::Tolerances::millimetres();
    let block = ogeom::algo::make_box(&mut model, ogeom::math::Frame::WORLD, (2.0, 3.0, 4.0), tol)
        .unwrap()
        .shape;
    let mesh =
        ogeom::mesh::triangulate(&model, &block, ogeom::mesh::Deflection::default(), tol).unwrap();
    let text = ogeom::io::mesh_formats::write_vrml(&[ogeom::io::mesh_formats::ExportMesh {
        mesh: &mesh,
        colour: Some([0.2, 0.4, 0.6, 1.0]),
        name: None,
    }]);
    let read = read_vrml(&text).unwrap();
    assert_eq!(read.len(), 1);
    assert_eq!(read[0].mesh.triangles.len(), mesh.triangles.len());
    assert!((volume(&read[0]) - 24.0).abs() < 1e-9);
    assert_eq!(read[0].colour, Some([0.2, 0.4, 0.6, 1.0]));
}

const SCENE: &str = r#"#VRML V2.0 utf8
# A unit cube spelt as quads, placed twice, and a switched-off ball.
PROTO Unused [ field SFFloat size 1 ] { Group {} }
DEF Cube Shape {
  appearance Appearance { material Material { diffuseColor 1 0 0 transparency 0.25 } }
  geometry IndexedFaceSet {
    coord Coordinate { point [ 0 0 0, 1 0 0, 1 1 0, 0 1 0, 0 0 1, 1 0 1, 1 1 1, 0 1 1 ] }
    coordIndex [ 0 3 2 1 -1, 4 5 6 7 -1, 0 1 5 4 -1, 2 3 7 6 -1, 0 4 7 3 -1, 1 2 6 5 ]
  }
}
Transform {
  translation 10 0 0
  rotation 0 0 1 1.5707963267948966
  scale 2 2 2
  children [ USE Cube ]
}
Switch { whichChoice -1 choice [ Shape { geometry Sphere { radius 5 } } ] }
Transform { translation 0 -5 0 children Shape { geometry Box { size 1 2 3 } } }
ROUTE A.b TO C.d
"#;

#[test]
fn a_scene_places_its_shapes_and_skips_what_it_switches_off() {
    let read = read_vrml(SCENE).unwrap();
    assert_eq!(read.len(), 3, "the cube twice and the box; the ball is off");
    assert!((volume(&read[0]) - 1.0).abs() < 1e-9);
    assert_eq!(read[0].colour, Some([1.0, 0.0, 0.0, 0.75]));
    // Scaled by two, turned a quarter about z, moved ten along x.
    assert!((volume(&read[1]) - 8.0).abs() < 1e-9);
    let (lo, hi) = bounds(&read[1]);
    assert!(lo.distance(Point::new(8.0, 0.0, 0.0)) < 1e-9, "{lo:?}");
    assert!(hi.distance(Point::new(10.0, 2.0, 2.0)) < 1e-9, "{hi:?}");
    assert!((volume(&read[2]) - 6.0).abs() < 1e-9);
    let (lo, hi) = bounds(&read[2]);
    assert!(lo.distance(Point::new(-0.5, -6.0, -1.5)) < 1e-9);
    assert!(hi.distance(Point::new(0.5, -4.0, 1.5)) < 1e-9);
}

const OLD: &str = r"#VRML V1.0 ascii
Separator {
  Material { diffuseColor 0 1 0 }
  Translation { translation 0 0 5 }
  Coordinate3 { point [ 0 0 0, 1 0 0, 0 1 0, 0 0 1 ] }
  ShapeHints { vertexOrdering COUNTERCLOCKWISE shapeType SOLID }
  IndexedFaceSet { coordIndex [ 0, 2, 1, -1, 0, 1, 3, -1, 0, 3, 2, -1, 1, 2, 3, -1 ] }
  Cylinder { parts (SIDES | TOP | BOTTOM) radius 1 height 2 }
}
Cube { width 1 height 1 depth 1 }
";

#[test]
fn a_first_version_file_draws_with_its_state() {
    let read = read_vrml(OLD).unwrap();
    assert_eq!(read.len(), 3);
    assert!(
        (volume(&read[0]) - 1.0 / 6.0).abs() < 1e-9,
        "{}",
        volume(&read[0])
    );
    assert_eq!(read[0].colour, Some([0.0, 1.0, 0.0, 1.0]));
    assert!(bounds(&read[0]).0.z > 4.9, "moved up by five");
    // A 32-sided drum, a little under the round one.
    let drum = volume(&read[1]);
    assert!(drum > 0.99 * 2.0 * core::f64::consts::PI * 0.99 && drum < 2.0 * core::f64::consts::PI);
    // The cube stands outside the separator: no colour, not moved.
    assert_eq!(read[2].colour, None);
    assert!((volume(&read[2]) - 1.0).abs() < 1e-9);
    assert!(bounds(&read[2]).0.z < 0.0);
}

#[test]
fn what_is_not_vrml_is_refused() {
    assert!(read_vrml("solid nothing\nendsolid").is_err());
    assert!(read_vrml("#VRML V2.0 utf8\nShape { geometry USE Missing }").is_err());
}
