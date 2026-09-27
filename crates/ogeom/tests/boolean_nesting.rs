//! Results of more than one shell: which shells are solids and which are
//! voids is decided by where they are, not by their bounding boxes.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn cube(model: &mut Model, at: (f64, f64, f64), size: (f64, f64, f64)) -> Shape {
    let frame = Frame::new(Point::new(at.0, at.1, at.2), Direction::Z, Direction::X, T).unwrap();
    ogeom::algo::make_box(model, frame, size, T).unwrap().shape
}

/// Each solid's volume, smallest first.
fn solid_volumes(model: &Model, shape: &Shape) -> Vec<f64> {
    let mut out: Vec<f64> = explore_unique(model, shape, ShapeType::Solid)
        .unwrap()
        .iter()
        .map(|s| {
            ogeom::algo::volume_properties(model, s, ogeom::mesh::Deflection::default(), T)
                .unwrap()
                .mass
        })
        .collect();
    out.sort_by(f64::total_cmp);
    out
}

/// A block standing free in an L-bracket's notch sits inside the bracket's
/// box and outside its material: two solids, neither a void of the other.
#[test]
fn a_block_in_a_notch_stays_a_solid_of_its_own() {
    let mut model = Model::new();
    let big = cube(&mut model, (0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
    let notch = cube(&mut model, (5.0, 5.0, -1.0), (6.0, 6.0, 12.0));
    let bracket = ogeom::boolean::cut(&mut model, &big, &notch, T)
        .unwrap()
        .shape;
    let block = cube(&mut model, (7.0, 7.0, 2.0), (2.0, 2.0, 2.0));
    let both = ogeom::boolean::fuse(&mut model, &bracket, &block, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &both, T).unwrap().is_valid());
    let volumes = solid_volumes(&model, &both);
    assert_eq!(volumes.len(), 2, "{volumes:?}");
    assert!((volumes[0] - 8.0).abs() < 1e-6, "{volumes:?}");
    assert!((volumes[1] - 750.0).abs() < 1e-6, "{volumes:?}");
}

/// A block floating in a hollow box's cavity: the box with its void, and
/// the block as a second solid, not a second void.
#[test]
fn an_island_in_a_cavity_is_a_solid_not_a_void() {
    let mut model = Model::new();
    let outer = cube(&mut model, (0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
    let cavity = cube(&mut model, (2.0, 2.0, 2.0), (6.0, 6.0, 6.0));
    let hollow = ogeom::boolean::cut(&mut model, &outer, &cavity, T)
        .unwrap()
        .shape;
    let island = cube(&mut model, (4.0, 4.0, 4.0), (2.0, 2.0, 2.0));
    let both = ogeom::boolean::fuse(&mut model, &hollow, &island, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &both, T).unwrap().is_valid());
    let volumes = solid_volumes(&model, &both);
    assert_eq!(volumes.len(), 2, "{volumes:?}");
    assert!((volumes[0] - 8.0).abs() < 1e-6, "{volumes:?}");
    assert!((volumes[1] - 784.0).abs() < 1e-6, "{volumes:?}");
}

/// The same soups through `make_volume`: faces of a hollow box and an
/// island, sewn and nested.
#[test]
fn make_volume_nests_an_island_as_a_solid() {
    let mut model = Model::new();
    let outer = cube(&mut model, (0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
    let cavity = cube(&mut model, (2.0, 2.0, 2.0), (6.0, 6.0, 6.0));
    let island = cube(&mut model, (4.0, 4.0, 4.0), (2.0, 2.0, 2.0));
    let mut faces = Vec::new();
    for solid in [&outer, &cavity, &island] {
        faces.extend(explore_unique(&model, solid, ShapeType::Face).unwrap());
    }
    let built = ogeom::boolean::make_volume(&mut model, &faces, T)
        .unwrap()
        .shape;
    let volumes = solid_volumes(&model, &built);
    assert_eq!(volumes.len(), 2, "{volumes:?}");
    assert!((volumes[0] - 8.0).abs() < 1e-6, "{volumes:?}");
    assert!((volumes[1] - 784.0).abs() < 1e-6, "{volumes:?}");
}

/// The three cells from one arrangement are the three booleans: the same
/// volumes as cutting each way and taking the common part separately, for
/// crossing boxes, a drum through a box, and boxes sharing a face.
#[test]
fn cells_are_the_three_booleans() {
    let volume = |model: &Model, shape: &Shape| {
        ogeom::algo::volume_properties(model, shape, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass
    };
    type Pair = fn(&mut Model) -> (Shape, Shape);
    let pairs: [Pair; 3] = [
        |m| {
            (
                cube(m, (0.0, 0.0, 0.0), (4.0, 4.0, 4.0)),
                cube(m, (2.0, 1.0, -1.0), (4.0, 2.0, 6.0)),
            )
        },
        |m| {
            let frame =
                Frame::new(Point::new(2.0, 2.0, -1.0), Direction::Z, Direction::X, T).unwrap();
            (
                cube(m, (0.0, 0.0, 0.0), (4.0, 4.0, 4.0)),
                ogeom::algo::make_cylinder(m, frame, 1.5, 6.0, T)
                    .unwrap()
                    .shape,
            )
        },
        |m| {
            (
                cube(m, (0.0, 0.0, 0.0), (4.0, 4.0, 4.0)),
                cube(m, (0.0, 0.0, 4.0), (4.0, 2.0, 3.0)),
            )
        },
    ];
    for make in pairs {
        let mut model = Model::new();
        let (a, b) = make(&mut model);
        let cells = ogeom::boolean::cells(&mut model, &a, &b, T).unwrap();
        let a_not_b = ogeom::boolean::cut(&mut model, &a, &b, T).unwrap().shape;
        let b_not_a = ogeom::boolean::cut(&mut model, &b, &a, T).unwrap().shape;
        let common = ogeom::boolean::common(&mut model, &a, &b, T).unwrap().shape;
        for (together, separate) in [
            (&cells.a_not_b.shape, &a_not_b),
            (&cells.b_not_a.shape, &b_not_a),
            (&cells.common.shape, &common),
        ] {
            assert!(ogeom::algo::check(&model, together, T).unwrap().is_valid());
            let (x, y) = (volume(&model, together), volume(&model, separate));
            assert!((x - y).abs() <= y.abs() * 1e-9 + 1e-9, "{x} against {y}");
        }
    }
}
