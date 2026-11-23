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

/// The surfaces of [`drilled`]'s block, in closed form: the six planes of
/// the box, then the bore about (12, 9) of `radius`.
fn surface_distances(p: Point, radius: f64) -> [f64; 7] {
    [
        p.x.abs(),
        (p.x - 30.0).abs(),
        p.y.abs(),
        (p.y - 20.0).abs(),
        p.z.abs(),
        (p.z - 10.0).abs(),
        ((p.x - 12.0).hypot(p.y - 9.0) - radius).abs(),
    ]
}

/// The closed-form areas of those surfaces' faces: the two faces the bore
/// pierces lose its disc, the bore is its full wall.
fn surface_areas(radius: f64) -> [f64; 7] {
    let disc = std::f64::consts::PI * radius * radius;
    let wall = 2.0 * std::f64::consts::PI * radius * 10.0;
    [200.0, 200.0, 300.0, 300.0, 600.0 - disc, 600.0 - disc, wall]
}

/// Each welded triangle names the face it lies on: its centroid is on that
/// face's surface (exactly on a plane, within the chord on the bore), each
/// face's named triangles cover its closed-form area, and the weld joins
/// exactly the face meshes' coincident vertices.
#[test]
fn every_welded_triangle_names_the_face_it_lies_on() {
    let radius = 4.0;
    let (model, solid) = drilled(radius);
    let deflection = Deflection::default();
    let welded = ogeom::mesh::triangulate_with_faces(&model, &solid, deflection, None, T).unwrap();
    assert_eq!(
        welded.mesh,
        ogeom::mesh::triangulate(&model, &solid, deflection, T).unwrap()
    );
    assert!(welded.mesh.is_closed());
    let mesh = &welded.mesh;
    assert_eq!(welded.triangle_faces.len(), mesh.triangles.len());

    // Each face's surface, read from its own mesh: the one every vertex of
    // it lies on.
    let chord = deflection.chord;
    let surface_of: Vec<usize> = welded
        .faces
        .iter()
        .map(|(_, piece)| {
            let piece = piece.as_ref().unwrap();
            let on: Vec<usize> = (0..7)
                .filter(|&s| {
                    piece
                        .positions
                        .iter()
                        .all(|&p| surface_distances(p, radius)[s] <= 1e-9)
                })
                .collect();
            assert_eq!(on.len(), 1, "a face lies on one surface");
            on[0]
        })
        .collect();

    let centroid = |t: [u32; 3]| {
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        Point::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
            (a.z + b.z + c.z) / 3.0,
        )
    };
    let mut areas = [0.0_f64; 7];
    let (mut plane_off, mut bore_off) = (0.0_f64, 0.0_f64);
    for (&t, &face) in mesh.triangles.iter().zip(&welded.triangle_faces) {
        let surface = surface_of[face as usize];
        let off = surface_distances(centroid(t), radius)[surface];
        if surface == 6 {
            bore_off = bore_off.max(off);
        } else {
            plane_off = plane_off.max(off);
        }
        let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
        areas[surface] += (b - a).cross(c - a).magnitude() / 2.0;
    }
    eprintln!("centroid off its surface: planes {plane_off:.3e}, bore {bore_off:.3e}");
    assert!(
        plane_off <= 1e-9,
        "a plane's triangle off it by {plane_off}"
    );
    assert!(bore_off <= chord, "a bore triangle off it by {bore_off}");

    // An inscribed polygon falls short of a circle's disc by at most its
    // sag times the perimeter, and a faceted wall short of the cylinder by
    // its area times sag over radius; flat rectangles are exact.
    let closed = surface_areas(radius);
    let rim = 2.0 * std::f64::consts::PI * radius;
    let bounds = [
        1e-9,
        1e-9,
        1e-9,
        1e-9,
        chord * rim,
        chord * rim,
        closed[6] * chord / radius,
    ];
    for s in 0..7 {
        let error = (areas[s] - closed[s]).abs();
        eprintln!(
            "surface {s}: area {:.6} closed {:.6} error {error:.3e} bound {:.3e}",
            areas[s], closed[s], bounds[s]
        );
        assert!(error <= bounds[s], "surface {s}: {error} > {}", bounds[s]);
    }

    // The weld, done here by brute force: the face meshes' vertices merged
    // where they coincide within the confusion tolerance.
    let mut kept: Vec<Point> = Vec::new();
    for (_, piece) in &welded.faces {
        for &p in &piece.as_ref().unwrap().positions {
            if !kept.iter().any(|&q| q.distance(p) <= T.confusion()) {
                kept.push(p);
            }
        }
    }
    let mut used = vec![false; mesh.positions.len()];
    for t in &mesh.triangles {
        for &i in t {
            used[i as usize] = true;
        }
    }
    let used = used.into_iter().filter(|&u| u).count();
    eprintln!("welded vertices {used}, test-side weld {}", kept.len());
    assert_eq!(used, kept.len());
}
