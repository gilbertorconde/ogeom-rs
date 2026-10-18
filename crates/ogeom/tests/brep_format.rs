//! The `.brep` interchange format: what the writer emits, the reader
//! returns, and a file written by hand from the specification reads as the
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

    // The sphere is still a sphere and the cylinder still a cylinder, read
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

    // And the degenerate edges (the sphere's poles) survived as degenerate
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
    // Four vertices, four edges on four lines, one wire, one planar face,
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

#[test]
fn placed_shapes_and_cone_charts_survive_the_round_trip() {
    // A cone and an L-shaped prism, the pair placed by a turn and a shift:
    // the placement has to be written as a placement record the reader
    // takes back, the cone's charts in the format's own terms (along the
    // ruling), and the prism's edges charted where each face uses them.
    let mut model = Model::new();
    let cone = ogeom::algo::make_cone(&mut model, Frame::WORLD, 5.0, 2.0, 10.0, T)
        .unwrap()
        .shape;
    let corners = [
        Point::new(20.0, 0.0, 0.0),
        Point::new(40.0, 0.0, 0.0),
        Point::new(40.0, 5.0, 0.0),
        Point::new(25.0, 5.0, 0.0),
        Point::new(25.0, 15.0, 0.0),
        Point::new(20.0, 15.0, 0.0),
    ];
    let outline = ogeom::algo::make_polygon(&mut model, &corners, true, T)
        .unwrap()
        .shape;
    let plane = ogeom::geom::PlaneSurface::new(ogeom::math::Plane::new(Frame::WORLD)).into();
    let ell = ogeom::algo::make_face(&mut model, plane, &[outline], T)
        .unwrap()
        .shape;
    let prism =
        ogeom::algo::make_prism(&mut model, &ell, ogeom::math::Vector::new(0.0, 0.0, 7.0), T)
            .unwrap()
            .shape;
    let both = ogeom::algo::make_compound(&mut model, &[cone, prism])
        .unwrap()
        .shape;
    let turn = ogeom::math::Transform::rotation(
        ogeom::math::Axis::new(Point::new(1.0, 2.0, 3.0), Direction::X),
        0.7,
    ) * ogeom::math::Transform::translation(ogeom::math::Vector::new(5.0, -4.0, 2.0));
    let datum = model.add_datum(turn);
    let placed = both.moved(&ogeom::topo::Location::of(datum));

    let text = ogeom::io::brep::write(&model, &placed, T).unwrap();
    let (back, shape) = ogeom::io::brep::read(&text, T).unwrap();
    assert_eq!(counts(&back, &shape), counts(&model, &placed));
    // A frustum of radii 5 and 2 over 10, and an L of area 150 over 7: both
    // integrate in closed form, so a billionth is rounding's share.
    let frustum = core::f64::consts::PI * 10.0 / 3.0 * (25.0 + 10.0 + 4.0);
    let exact = frustum + 150.0 * 7.0;
    let now = ogeom::algo::volume_properties(&back, &shape, Deflection::default(), T).unwrap();
    assert!(
        (now.mass - exact).abs() < exact * 1e-9,
        "{} against {exact}",
        now.mass
    );
    let was = ogeom::algo::volume_properties(&model, &placed, Deflection::default(), T).unwrap();
    assert!(
        now.centre.distance(was.centre) < 1e-9,
        "and where it was: {:?} against {:?}",
        now.centre,
        was.centre
    );
}

#[test]
fn a_read_extrusion_survives_the_native_format() {
    // An elliptic prism's wall is an extrusion; read back from `.brep` it
    // reaches both ways from its curve, and the native format has to carry
    // that window through a write and a read.
    let mut model = Model::new();
    let ellipse = ogeom::geom::Curve::Ellipse(ogeom::geom::EllipseCurve::new(
        ogeom::math::Ellipse::new(Frame::WORLD, 10.0, 5.0, T).unwrap(),
    ));
    let range = ogeom::geom::Curve3d::domain(&ellipse);
    let edge = ogeom::algo::make_edge(&mut model, ellipse, range, T)
        .unwrap()
        .shape;
    let plane =
        ogeom::geom::PlaneSurface::over(ogeom::math::Plane::XY, (-20.0, 20.0), (-20.0, 20.0))
            .unwrap();
    let face = ogeom::algo::make_face_with_pcurves(&mut model, plane.into(), &[vec![edge]], T)
        .unwrap()
        .shape;
    let prism = ogeom::algo::make_prism(
        &mut model,
        &face,
        ogeom::math::Vector::new(0.0, 0.0, 4.0),
        T,
    )
    .unwrap()
    .shape;

    let text = ogeom::io::brep::write(&model, &prism, T).unwrap();
    let (read, solid) = ogeom::io::brep::read(&text, T).unwrap();
    assert!(
        read.geometry()
            .surfaces()
            .any(|(_, s)| matches!(s, ogeom::geom::SurfaceGeometry::Extrusion(_))),
        "the wall is an extrusion"
    );
    let options = ogeom::io::native::WriteOptions::default();
    let native = ogeom::io::native::write(&read, std::slice::from_ref(&solid), options).unwrap();
    let (back, roots) = ogeom::io::native::read(&native).unwrap();
    assert_eq!(
        ogeom::io::native::write(&back, &roots, options).unwrap(),
        native
    );
    assert_eq!(counts(&back, &roots[0]), counts(&read, &solid));
    assert!(ogeom::algo::check(&back, &roots[0], T).unwrap().is_valid());
    let (now, was) = (volume(&back, &roots[0]), volume(&model, &prism));
    assert!((now - was).abs() < 1e-9 * was, "{now} against {was}");
}
