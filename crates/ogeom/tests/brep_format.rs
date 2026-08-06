//! The `.brep` interchange format: what the writer emits, the reader
//! returns — and a file written by hand from the specification reads as the
//! shape it describes.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

fn counts(model: &Model, shape: &Shape) -> (usize, usize, usize, usize) {
    (
        explore_unique(model, shape, ShapeType::Vertex)
            .unwrap()
            .len(),
        explore_unique(model, shape, ShapeType::Edge).unwrap().len(),
        explore_unique(model, shape, ShapeType::Face).unwrap().len(),
        explore_unique(model, shape, ShapeType::Solid)
            .unwrap()
            .len(),
    )
}

#[test]
fn a_drilled_block_survives_the_round_trip() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 12.0), T)
        .unwrap()
        .shape;
    let bore = ogeom::algo::make_cylinder(
        &mut model,
        Frame::new(Point::new(30.0, 15.0, -1.0), Direction::Z, Direction::X, T).unwrap(),
        4.0,
        14.0,
        T,
    )
    .unwrap()
    .shape;
    let part = ogeom::boolean::cut(&mut model, &block, &bore, T)
        .unwrap()
        .shape;

    let text = ogeom::io::brep::write(&model, &part, T).unwrap();
    assert!(
        text.starts_with("DBRep_DrawableShape"),
        "the file says what it is"
    );

    let (back, shape) = ogeom::io::brep::read(&text, T).unwrap();
    assert_eq!(
        counts(&back, &shape),
        counts(&model, &part),
        "every vertex, edge, face and solid comes back"
    );
    let (was, now) = (volume(&model, &part), volume(&back, &shape));
    assert!(
        (now - was).abs() < was * 1e-9,
        "and it is the same solid: {now} against {was}"
    );

    // And the text is a fixed point from there: write, read, write again,
    // and the bytes are identical, so the format carries everything the
    // reader put back and the reader takes everything the format carries.
    //
    // The comparison starts at the *second* text rather than the first,
    // because one flag legitimately changes on the way in: the boolean's
    // own edges never had their same-parameter claim re-established, and
    // the reader measures the representations and establishes it.
    let again = ogeom::io::brep::write(&back, &shape, T).unwrap();
    let (twice, shape) = ogeom::io::brep::read(&again, T).unwrap();
    let third = ogeom::io::brep::write(&twice, &shape, T).unwrap();
    assert_eq!(third, again, "the text is a fixed point");
}

#[test]
fn seams_and_poles_come_back_as_themselves() {
    // A drum and a ball: between them a cylinder with a seam, a sphere with
    // a seam and two degenerate poles, two planes, and circular edges. The
    // format has to carry all of it, including the edges that are points.
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T)
        .unwrap()
        .shape;
    let ball = ogeom::algo::make_sphere(
        &mut model,
        Frame::new(Point::new(30.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap(),
        5.0,
        T,
    )
    .unwrap()
    .shape;
    let both = ogeom::algo::make_compound(&mut model, &[drum, ball])
        .unwrap()
        .shape;

    let text = ogeom::io::brep::write(&model, &both, T).unwrap();
    let (back, shape) = ogeom::io::brep::read(&text, T).unwrap();
    assert_eq!(counts(&back, &shape), counts(&model, &both));

    // The sphere is still a sphere and the cylinder still a cylinder — read
    // back as the surfaces they are, not fitted into splines.
    let kinds = |model: &Model, shape: &Shape| -> (usize, usize) {
        let faces = explore_unique(model, shape, ShapeType::Face).unwrap();
        let mut spheres = 0;
        let mut cylinders = 0;
        for face in faces {
            let ogeom::topo::NodeData::Face(data) = model.node(&face).unwrap().data() else {
                continue;
            };
            match model.geometry().surface(data.surface) {
                Some(ogeom::geom::SurfaceGeometry::Sphere(_)) => spheres += 1,
                Some(ogeom::geom::SurfaceGeometry::Cylinder(_)) => cylinders += 1,
                _ => {}
            }
        }
        (spheres, cylinders)
    };
    assert_eq!(kinds(&back, &shape), kinds(&model, &both));
    assert_eq!(kinds(&back, &shape), (1, 1));

    // And the degenerate edges — the sphere's poles — survived as degenerate
    // edges rather than as nothing.
    let poles = |model: &Model, shape: &Shape| -> usize {
        explore_unique(model, shape, ShapeType::Edge)
            .unwrap()
            .iter()
            .filter(|e| {
                model
                    .node(e)
                    .and_then(|n| n.data().as_edge())
                    .is_some_and(|d| d.degenerate)
            })
            .count()
    };
    assert_eq!(poles(&back, &shape), poles(&model, &both));
    assert_eq!(poles(&back, &shape), 2, "the ball has two poles");
}

#[test]
fn a_file_written_by_hand_reads_as_the_square_it_describes() {
    // Four vertices, four edges on four lines, one wire, one planar face —
    // spelled out against the specification rather than produced by the
    // writer, so the reader is tested against the format and not against
    // its own dialect.
    let text = "\
DBRep_DrawableShape

CASCADE Topology V1, (c) Matra-Datavision
Locations 0
Curve2ds 4
1 0 0 1 0
1 10 0 0 1
1 10 10 -1 0
1 0 10 0 -1
Curves 4
1 0 0 0 1 0 0
1 10 0 0 0 1 0
1 10 10 0 -1 0 0
1 0 10 0 0 -1 0
Polygon3D 0
PolygonOnTriangulations 0
Surfaces 1
1 0 0 0 0 0 1 1 0 0 0 1 0
Triangulations 0

TShapes 10
Ve
1e-07
0 0 0
0 0

0001000
*
Ve
1e-07
10 0 0
0 0

0001000
*
Ed
 1e-07 1 1 0
1 1 0 0 10
2 1 1 0 0 10
0

0001000
+10 0 -9 0 *
Ve
1e-07
10 10 0
0 0

0001000
*
Ed
 1e-07 1 1 0
1 2 0 0 10
2 2 1 0 0 10
0

0001000
+9 0 -7 0 *
Ve
1e-07
0 10 0
0 0

0001000
*
Ed
 1e-07 1 1 0
1 3 0 0 10
2 3 1 0 0 10
0

0001000
+7 0 -5 0 *
Ed
 1e-07 1 1 0
1 4 0 0 10
2 4 1 0 0 10
0

0001000
+5 0 -10 0 *
Wi

0001000
+8 0 +6 0 +4 0 +3 0 *
Fa
0 1e-07 1 0

0001000
+2 0 *

+1 0
";
    let (model, face) = ogeom::io::brep::read(text, T).unwrap();
    assert_eq!(model.kind_of(&face).unwrap(), ShapeType::Face);
    let (vertices, edges, faces, _) = counts(&model, &face);
    assert_eq!((vertices, edges, faces), (4, 4, 1));

    // It is the ten-by-ten square in the z = 0 plane: its area says so.
    let area = ogeom::algo::surface_properties(&model, &face, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!((area - 100.0).abs() < 1e-9, "a ten by ten square: {area}");
}
