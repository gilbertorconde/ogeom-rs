//! Solids touching face to face at a placement a single-precision solver
//! leaves a rounding off coplanar: the common holds no material, and the
//! circles the contact plane cuts from bores standing square to it are
//! read on the bores' walls as the coplanar ones are, without marching.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_box, make_cylinder, volume_properties};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Matrix3, Point, Transform, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Location, Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The placement printCAD's assembly solver gave: the plate's top face on
/// the box's bottom one, tilted by about 1.4e-8 rad.
const PLACED: [[f64; 4]; 3] = [
    [
        0.999_999_999_999_999_9,
        -1.401_389_136_184_109_4e-8,
        -2.687_303_811_083_958e-9,
        -40.999_851_226_806_64,
    ],
    [
        1.401_389_139_950_067_5e-8,
        0.999_999_999_999_999_8,
        1.401_389_047_366_267_4e-8,
        -5.000_000_476_837_158,
    ],
    [
        2.687_303_614_694_819e-9,
        -1.401_389_051_132_225_5e-8,
        0.999_999_999_999_999_9,
        -10.0,
    ],
];

/// A 60 x 40 x 10 plate less four bores of radius 3 through it.
fn plate(model: &mut Model) -> Shape {
    let mut plate = make_box(model, Frame::WORLD, (60.0, 40.0, 10.0), T)
        .unwrap()
        .shape;
    for (x, y) in [(10.0, 10.0), (50.0, 10.0), (10.0, 30.0), (50.0, 30.0)] {
        let at = Frame::new(Point::new(x, y, 0.0), Direction::Z, Direction::X, T).unwrap();
        let bore = make_cylinder(model, at, 3.0, 10.0, T).unwrap().shape;
        plate = ogeom::boolean::cut(model, &plate, &bore, T).unwrap().shape;
    }
    plate
}

/// The common of a 10-cube with the plate placed under it, tilted as
/// printCAD placed it or with the same translation alone: the volume of
/// each solid in it, and how many face pairs were marched.
fn common(tilted: bool) -> (Vec<f64>, usize) {
    let mut model = Model::new();
    let cube = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let plate = plate(&mut model);
    let linear = if tilted {
        Matrix3::new(PLACED.map(|row| [row[0], row[1], row[2]]))
    } else {
        Matrix3::IDENTITY
    };
    let shift = Vector::new(PLACED[0][3], PLACED[1][3], PLACED[2][3]);
    let placed = Transform::from_parts(linear, 1.0, shift, 1e-9).unwrap();
    let plate = plate.located(Location::of(model.add_datum(placed)));
    let (built, work) =
        ogeom_bool::work::counted(|| ogeom::boolean::common(&mut model, &cube, &plate, T));
    let built = built.unwrap().shape;
    let volumes = explore_unique(&model, &built, ShapeType::Solid)
        .unwrap()
        .iter()
        .map(|solid| {
            volume_properties(&model, solid, Deflection::with_chord(1e-3).unwrap(), T)
                .unwrap()
                .mass
        })
        .collect();
    (volumes, work.marched_pairs)
}

/// Tilted or not, the cube and the plate only touch: no solid of any
/// volume, and no pair marched. The cube's bottom plane cuts the bore at
/// (50, 10) in a circle a rounding off square to its axis, which climbs a
/// hair past the wall's top round half its turn; its image on the wall is
/// read there in closed form, as it is on a cone.
#[test]
fn a_near_coplanar_contact_beside_a_bore_shares_nothing_and_marches_nothing() {
    for tilted in [false, true] {
        let (volumes, marched) = common(tilted);
        assert!(
            volumes.iter().all(|v| v.abs() <= 1e-9),
            "tilted {tilted}: solids of volume {volumes:?}"
        );
        assert_eq!(marched, 0, "tilted {tilted}: pairs marched");
    }
}
