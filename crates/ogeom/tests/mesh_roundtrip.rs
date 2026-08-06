//! Curved mesh → B-rep: a drum's mesh, topology thrown away, comes back
//! with its wall recognized as the cylinder it is; a genuinely free-form
//! mesh still refuses by name.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::SurfaceGeometry;
use ogeom::math::Frame;
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, NodeData, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_drum_survives_the_round_trip() {
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 12.0, T)
        .unwrap()
        .shape;
    let true_volume = core::f64::consts::PI * 25.0 * 12.0;

    let mesh =
        ogeom::mesh::triangulate(&model, &drum, Deflection::with_chord(0.01).unwrap(), T).unwrap();
    let mut fresh = Model::new();
    let rebuilt = ogeom::heal::mesh_to_brep(&mut fresh, &mesh, ogeom::algo::CREASE, 0.02, T)
        .unwrap()
        .shape;

    assert_eq!(fresh.kind_of(&rebuilt).unwrap(), ShapeType::Solid);
    // The wall is a cylinder again, radius five within the chord.
    let mut cylinders = 0;
    for face in explore_unique(&fresh, &rebuilt, ShapeType::Face).unwrap() {
        let NodeData::Face(data) = fresh.node(&face).unwrap().data() else {
            continue;
        };
        if let Some(SurfaceGeometry::Cylinder(c)) = fresh.geometry().surface(data.surface) {
            cylinders += 1;
            assert!(
                (c.cylinder().radius() - 5.0).abs() < 0.02,
                "the radius survives: {}",
                c.cylinder().radius()
            );
        }
    }
    assert_eq!(cylinders, 1, "one recognized wall");

    // The B-rep is whole: three faces, every edge used twice.
    let faces = explore_unique(&fresh, &rebuilt, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 3, "two caps and one wall");
    let shell = explore_unique(&fresh, &rebuilt, ShapeType::Shell).unwrap();
    assert!(ogeom::algo::is_shell_closed(&fresh, &shell[0]).unwrap());

    // And it can be *measured*, which is the thing a shape has to be for any
    // claim about it to mean anything. The wall wraps its surface's whole
    // period, so its chart region is closed by a seam rather than by its rims
    // — without one the boundary encloses no area and the face does not
    // tessellate at all.
    let measured =
        ogeom::algo::volume_properties(&fresh, &rebuilt, Deflection::with_chord(0.01).unwrap(), T)
            .expect("a rebuilt drum has a volume")
            .mass;
    // The rims came back as the polygons the mesh had, so the recovered solid
    // is the drum's inscribed prism and reads a little under. The mesh was
    // built at a chord of 0.01 on a radius of five, which is about a
    // hundred-and-ten-gon; an inscribed n-gon loses about `π²/(3n²)` of the
    // circle's area, and the measurement is held to that rather than to the
    // circle it approximates.
    let sides = explore_unique(&fresh, &rebuilt, ShapeType::Edge)
        .unwrap()
        .len()
        .max(3);
    #[allow(clippy::cast_precision_loss, reason = "an edge count")]
    let deficit = core::f64::consts::PI.powi(2) / (3.0 * (sides as f64).powi(2));
    assert!(
        measured <= true_volume && measured > true_volume * (1.0 - 20.0 * deficit),
        "the inscribed drum: {measured} against {true_volume}, deficit {deficit}"
    );
}
