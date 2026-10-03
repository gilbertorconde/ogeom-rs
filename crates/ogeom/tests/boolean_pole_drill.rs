//! Drills along a ball's axis whose wall passes through both poles of its
//! chart, or a hair beside them.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_cylinder, make_sphere, volume_properties};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();
const BALL: f64 = 10.0;

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

/// The ball's volume inside a vertical drill of radius `r` whose axis
/// stands at `(cx, cy)`: twice the ball's half height integrated over the
/// drill's disc, in the disc's own polar coordinates, where the integrand
/// is smooth (the disc lies within the ball's outline).
fn drilled_volume(cx: f64, cy: f64, r: f64) -> f64 {
    let (rings, spokes) = (2000_u32, 2000_u32);
    let ds = r / f64::from(rings);
    let dpsi = core::f64::consts::TAU / f64::from(spokes);
    let mut total = 0.0;
    for j in 0..spokes {
        let psi = dpsi * f64::from(j);
        let (c, s) = (psi.cos(), psi.sin());
        // Simpson along the radius.
        let mut line = 0.0;
        for i in 0..=rings {
            let t = ds * f64::from(i);
            let (x, y) = (cx + t * c, cy + t * s);
            let f = 2.0 * (BALL * BALL - x * x - y * y).sqrt() * t;
            let w = if i == 0 || i == rings {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            line += w * f;
        }
        total += line * ds / 3.0;
    }
    total * dpsi
}

/// Cut and common of the ball with a drill along its axis, against the
/// ball's volume and the drilled volume integrated independently.
fn drilled(offset: Vector, r: f64) -> Result<(), String> {
    let mut model = Model::new();
    let ball = make_sphere(&mut model, Frame::WORLD, BALL, T)
        .unwrap()
        .shape;
    let at = Frame::new(
        Point::new(offset.x, offset.y, -33.46410174977877),
        Direction::Z,
        Direction::Y,
        T,
    )
    .unwrap();
    let drill = make_cylinder(&mut model, at, r, 46.92820349955754, T)
        .unwrap()
        .shape;
    let mut volumes = Vec::new();
    for (name, made) in [
        ("cut", ogeom::boolean::cut(&mut model, &ball, &drill, T)),
        (
            "common",
            ogeom::boolean::common(&mut model, &ball, &drill, T),
        ),
    ] {
        let made = made.map_err(|e| format!("{name}: {e}"))?;
        let diagnosis = check(&model, &made.shape, T).unwrap();
        if !diagnosis.is_valid() {
            return Err(format!("{name}: {diagnosis}"));
        }
        volumes.push(volume(&model, &made.shape));
    }
    let whole = 4.0 / 3.0 * core::f64::consts::PI * BALL.powi(3);
    let inside = drilled_volume(offset.x, offset.y, r);
    let (cut, common) = (volumes[0], volumes[1]);
    if ((cut + common) - whole).abs() > whole * 1e-6 {
        return Err(format!("cut {cut} + common {common} against {whole}"));
    }
    if (common - inside).abs() > inside * 1e-6 {
        return Err(format!("common {common} against {inside}"));
    }
    Ok(())
}

#[test]
fn drill_through_both_poles() {
    let r = 2.664_527_471_227_642;
    let mut failures = Vec::new();
    for (dx, dy) in [
        (0.0, -1.0),
        (0.0, 1.0),
        (1.0, 0.0),
        (-1.0, 0.0),
        (
            core::f64::consts::FRAC_1_SQRT_2,
            core::f64::consts::FRAC_1_SQRT_2,
        ),
        (0.6, -0.8),
    ] {
        for radius in [r, 1.0, 4.0] {
            for wall in [
                radius,
                radius * (1.0 + 1e-7),
                radius * (1.0 - 1e-7),
                radius * 1.001,
            ] {
                let offset = Vector::new(dx * wall, dy * wall, 0.0);
                if let Err(e) = drilled(offset, radius) {
                    failures.push(format!("({dx}, {dy}) r {radius} at {wall}: {e}"));
                }
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
