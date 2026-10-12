//! Solids offset by `offset_shape`, their topology kept one for one, held
//! to the volumes their grown sizes give.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

fn torus_volume(major: f64, minor: f64) -> f64 {
    2.0 * PI * PI * major * minor * minor
}

/// A whole torus is one face bounded by its two seams, each walked twice,
/// crossing at one vertex: no ends for the faces' wires to be assembled
/// from, and the offset solid did not close. It is the moved torus whole,
/// out and in, on an oblique axis as on the world's.
#[test]
fn a_whole_torus_offsets_out_and_in() {
    let tilted = Frame::new(
        Point::new(1.0, -2.0, 0.5),
        Direction::new(Vector::new(0.3, -0.5, 0.8), T).unwrap(),
        Direction::new(Vector::new(0.0, 0.8, 0.5), T).unwrap(),
        T,
    )
    .unwrap();
    for frame in [Frame::WORLD, tilted] {
        for d in [0.4, -0.3] {
            let mut model = Model::new();
            let ring = ogeom::algo::make_torus(&mut model, frame, 3.0, 1.0, T)
                .unwrap()
                .shape;
            let grown = ogeom::offset::offset_shape(&mut model, &ring, d, T)
                .unwrap()
                .shape;
            let diagnosis = ogeom::algo::check(&model, &grown, T).unwrap();
            assert!(diagnosis.is_valid(), "offset {d}: {diagnosis}");
            let (got, want) = (volume(&model, &grown), torus_volume(3.0, 1.0 + d));
            assert!(
                (got - want).abs() <= want * 1e-9,
                "offset {d}: {got} against {want}"
            );
        }
    }
}

/// A block with a torus-shaped void: grown, the block's faces move out and
/// the void's tube narrows by as much.
#[test]
fn a_torus_void_narrows_as_its_block_grows() {
    let mut model = Model::new();
    let corner = Frame::new(Point::new(-6.0, -6.0, -3.0), Direction::Z, Direction::X, T).unwrap();
    let block = ogeom::algo::make_box(&mut model, corner, (12.0, 12.0, 6.0), T)
        .unwrap()
        .shape;
    let ring = ogeom::algo::make_torus(&mut model, Frame::WORLD, 3.0, 1.0, T)
        .unwrap()
        .shape;
    let holed = ogeom::boolean::cut(&mut model, &block, &ring, T)
        .unwrap()
        .shape;
    let d = 0.25;
    let grown = ogeom::offset::offset_shape(&mut model, &holed, d, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom::algo::check(&model, &grown, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let want = (12.0 + 2.0 * d) * (12.0 + 2.0 * d) * (6.0 + 2.0 * d) - torus_volume(3.0, 1.0 - d);
    let got = volume(&model, &grown);
    assert!((got - want).abs() <= want * 1e-9, "{got} against {want}");
}
