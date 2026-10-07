//! Each face's mesh as the shape-wide pass draws it, and as a cache of face
//! meshes kept across models hands it back: the same as drawing that face
//! on its own with the shape's agreed chords.
#![allow(clippy::unwrap_used, reason = "test code")]

use ogeom::algo::{make_box, make_cylinder};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::{Deflection, FaceMeshCache};
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

/// A block drilled through by a hole of `radius`, in a model of its own.
fn drilled(radius: f64) -> (Model, Shape) {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (30.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(12.0, 9.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let hole = make_cylinder(&mut model, at, radius, 12.0, T)
        .unwrap()
        .shape;
    let solid = ogeom::boolean::cut(&mut model, &block, &hole, T)
        .unwrap()
        .shape;
    (model, solid)
}

/// Every face's mesh beside the whole is what drawing that face alone to
/// the returned chords gives, and the whole is the shape's own mesh;
/// through `kept`, the same again.
fn assert_drawn_alike(model: &Model, solid: &Shape, kept: Option<&FaceMeshCache>) {
    let deflection = Deflection::default();
    let (whole, chords, faces) =
        ogeom::mesh::triangulate_with_face_meshes(model, solid, deflection, kept, T).unwrap();
    assert_eq!(
        whole,
        ogeom::mesh::triangulate(model, solid, deflection, T).unwrap()
    );
    assert!(!faces.is_empty());
    for (face, mesh) in faces {
        let alone = ogeom::mesh::triangulate_face_with(model, &face, deflection, &chords, T);
        assert_eq!(mesh.unwrap(), alone.unwrap());
    }
}

#[test]
fn face_meshes_beside_the_whole_are_each_face_drawn_alone() {
    let (model, solid) = drilled(4.0);
    assert_drawn_alike(&model, &solid, None);
}

/// A cache filled from one model answers for another only where a face is
/// the same to the bit: the same block built again takes its meshes from
/// it, and a block whose hole is a little wider (the bore, and the two
/// faces it pierces, differ; the other four do not) gets its own.
#[test]
fn kept_face_meshes_answer_as_each_face_drawn_anew() {
    let kept = FaceMeshCache::new();
    for radius in [4.0, 4.0, 4.001, 4.0] {
        let (model, solid) = drilled(radius);
        assert_drawn_alike(&model, &solid, Some(&kept));
        for face in
            ogeom::topo::explore_unique(&model, &solid, ogeom::topo::ShapeType::Face).unwrap()
        {
            let deflection = Deflection::with_chord(0.05).unwrap();
            let fresh = ogeom::mesh::triangulate_face(&model, &face, deflection, T).unwrap();
            let held =
                ogeom::mesh::triangulate_face_kept(&model, &face, deflection, &kept, T).unwrap();
            assert_eq!(held, fresh);
        }
    }
}

/// The same face built in two models reads the same; a face of another
/// size reads otherwise.
#[test]
fn face_content_tells_faces_apart_by_geometry_not_handles() {
    let faces = |radius: f64| {
        let (model, solid) = drilled(radius);
        ogeom::topo::explore_unique(&model, &solid, ogeom::topo::ShapeType::Face)
            .unwrap()
            .iter()
            .map(|f| ogeom::topo::face_content(&model, f).unwrap())
            .collect::<Vec<_>>()
    };
    let (a, again, wider) = (faces(4.0), faces(4.0), faces(4.001));
    assert_eq!(a, again);
    assert_eq!(a.len(), wider.len());
    let differing = a.iter().zip(&wider).filter(|(x, y)| x != y).count();
    // The bore and the top and bottom faces it pierces.
    assert_eq!(differing, 3);
}
