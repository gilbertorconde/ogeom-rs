//! Holes drilled through a ball close to its rim, the drill all but grazing
//! the ball's far side and passing a few millimetres from a pole of its
//! chart: every line along the drill still meets the ball twice, and the
//! two loops the drill leaves are long and thin.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_cylinder, make_sphere, volume_properties};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

/// Cut and common with each hole are valid and share the ball's volume.
#[test]
fn holes_near_a_balls_rim_cut_and_fill_valid() {
    let mut model = Model::new();
    let pole = Direction::new(Vector::new(-0.857_86, 0.513_88, 0.0), T).unwrap();
    let along = Direction::new(Vector::new(0.513_88, 0.857_86, 0.0), T).unwrap();
    let radius = 8.88;
    let ball = make_sphere(
        &mut model,
        Frame::new(Point::ORIGIN, pole, along, T).unwrap(),
        radius,
        T,
    )
    .unwrap()
    .shape;
    let whole = 4.0 / 3.0 * core::f64::consts::PI * radius.powi(3);
    for (x, z, r) in [(-7.77, -2.49, 0.67), (-7.9, -2.0, 0.7)] {
        let at = Frame::new(Point::new(x, -22.0, z), Direction::Y, Direction::Z, T).unwrap();
        let drill = make_cylinder(&mut model, at, r, 44.0, T).unwrap().shape;
        let mut shares = 0.0;
        for (name, made) in [
            ("cut", ogeom::boolean::cut(&mut model, &ball, &drill, T)),
            (
                "common",
                ogeom::boolean::common(&mut model, &ball, &drill, T),
            ),
        ] {
            let made = made.unwrap_or_else(|e| panic!("{name} at ({x}, {z}): {e}"));
            let diagnosis = check(&model, &made.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{name} at ({x}, {z}): {diagnosis}");
            shares += volume(&model, &made.shape);
        }
        assert!(
            (shares - whole).abs() < whole * 1e-4,
            "at ({x}, {z}): {shares} against {whole}"
        );
    }
}
