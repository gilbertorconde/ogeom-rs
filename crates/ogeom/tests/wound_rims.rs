//! A face bounded by two rims that each wind a periodic direction of its
//! chart — a torus band between two parallels, a ball's belt between two
//! latitudes — is the strip between them, and draws as that strip.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::SurfaceGeometry;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, Triangulation, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

fn area_of(mesh: &Triangulation) -> f64 {
    mesh.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            (b - a).cross(c - a).magnitude() * 0.5
        })
        .sum()
}

fn surface_of(model: &Model, face: &Shape) -> SurfaceGeometry {
    let data = model.node(face).unwrap().data().as_face().unwrap().clone();
    model.geometry().surface(data.surface).unwrap().clone()
}

/// Two torus bands, each bounded by two parallels — one full circle round
/// the axis at each end, every wire a single wound edge. Walked one wire at
/// a time, each rim used to close on its own translate a tube-period over,
/// which is the whole torus cut along the rim; two of those laid over each
/// other cancel where they overlap, and the face drew to more area than
/// the torus has. The rims are paired into one ring first, and the face is
/// the band.
#[test]
fn a_torus_band_between_two_rims_is_the_strip_between_them() {
    let text = corpus("nist_ftc_11_asme1_rb.stp");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let mut seen = 0;
    for face in explore_unique(model, &import.solids[0], ShapeType::Face).unwrap() {
        let SurfaceGeometry::Torus(t) = surface_of(model, &face) else {
            continue;
        };
        let (major, minor) = (t.torus().major_radius(), t.torus().minor_radius());
        // The bands: from the outer equator over the top to the parallel
        // at height 0.425 on the far side, and from that height on the
        // near side over the top to the inner equator.
        let lift = (0.425_f64 / minor).asin();
        let (v0, v1) = if (major - 30.0).abs() < 1e-6 {
            (0.0, std::f64::consts::PI - lift)
        } else {
            (lift, std::f64::consts::PI)
        };
        let band =
            std::f64::consts::TAU * minor * (major * (v1 - v0) + minor * (v1.sin() - v0.sin()));
        let mesh = ogeom::mesh::triangulate(model, &face, Deflection::with_chord(1e-2).unwrap(), T)
            .unwrap();
        let area = area_of(&mesh);
        assert!(
            (area - band).abs() <= band * 0.015,
            "a torus band of major radius {major}: {area:.1} drawn against {band:.1} owed"
        );
        seen += 1;
    }
    assert_eq!(seen, 2, "the part carries two torus bands");
}

/// A ball cut flat above and below: the belt between the two cuts is a
/// sphere face bounded by two latitude circles, each winding the chart
/// once. The same pairing that holds the torus band holds the belt, where
/// a rim closed alone would close against the nearer pole as a cap.
#[test]
fn a_belt_between_two_latitudes_is_the_belt() {
    let mut model = Model::new();
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 10.0, T)
        .unwrap()
        .shape;
    let slab = |model: &mut Model, z: f64| {
        let frame = Frame::new(Point::new(-20.0, -20.0, z), Direction::Z, Direction::X, T).unwrap();
        ogeom::algo::make_box(model, frame, (40.0, 40.0, 20.0), T)
            .unwrap()
            .shape
    };
    let above = slab(&mut model, 3.0);
    let below = slab(&mut model, -23.0);
    let topless = ogeom::boolean::cut(&mut model, &ball, &above, T)
        .unwrap()
        .shape;
    let belt = ogeom::boolean::cut(&mut model, &topless, &below, T)
        .unwrap()
        .shape;
    let face = explore_unique(&model, &belt, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| matches!(surface_of(&model, f), SurfaceGeometry::Sphere(_)))
        .expect("the belt keeps its sphere face");
    let mesh =
        ogeom::mesh::triangulate(&model, &face, Deflection::with_chord(1e-2).unwrap(), T).unwrap();
    let owed = std::f64::consts::TAU * 10.0 * 6.0;
    let area = area_of(&mesh);
    assert!(
        (area - owed).abs() <= owed * 0.005,
        "the belt's own area: {area:.2} drawn against {owed:.2} owed"
    );
}

/// A bore through a block with a cross hole through the bore's wall: the
/// wall is a cylinder whose two rims wind the chart and whose cross hole
/// is a loop that does not. The rims pair into one ring joined by runs up
/// one column. Taken wherever the first rim's chain happened to end — on
/// this part, where the cross hole's loop lies — the runs cut straight
/// through the loop: two constraints refused, the face redrawn finer and
/// finer, never whole. The runs stand in the widest column no other ring
/// touches, and the body draws closed.
#[test]
fn a_cross_hole_through_a_bore_wall_draws_closed() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (30.0, 30.0, 30.0), T)
        .unwrap()
        .shape;
    let bore_frame =
        Frame::new(Point::new(15.0, 15.0, -5.0), Direction::Z, Direction::X, T).unwrap();
    let bore = ogeom::algo::make_cylinder(&mut model, bore_frame, 8.0, 40.0, T)
        .unwrap()
        .shape;
    let cross_frame =
        Frame::new(Point::new(-5.0, 15.0, 15.0), Direction::X, Direction::Z, T).unwrap();
    let cross = ogeom::algo::make_cylinder(&mut model, cross_frame, 3.0, 40.0, T)
        .unwrap()
        .shape;
    let bored = ogeom::boolean::cut(&mut model, &block, &bore, T)
        .unwrap()
        .shape;
    let holed = ogeom::boolean::cut(&mut model, &bored, &cross, T)
        .unwrap()
        .shape;
    let mesh = ogeom::mesh::triangulate(&model, &holed, Deflection::default(), T).unwrap();
    assert!(mesh.is_closed(), "the drilled block draws closed");
    let expected = 27000.0
        - core::f64::consts::PI * 64.0 * 30.0
        - 2.0 * core::f64::consts::PI * 9.0 * (15.0 - 8.0);
    let measured =
        ogeom::algo::volume_properties(&model, &holed, Deflection::with_chord(1e-2).unwrap(), T)
            .unwrap()
            .mass;
    // The cross hole's two legs each reach from the block's face to the
    // bore; the small overlap where a leg meets the bore's curve is inside
    // the percent.
    assert!(
        (measured - expected).abs() < expected * 1e-2,
        "{measured} against about {expected}"
    );
}
