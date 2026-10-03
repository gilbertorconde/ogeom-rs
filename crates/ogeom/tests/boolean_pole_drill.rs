//! Drills along a ball's axis whose wall passes through both poles of its
//! chart, or a hair beside them, and drills through the pole of a dome: a
//! ball cut below its pole, alone or standing on a drum.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_box, make_cylinder, make_sphere, volume_properties};
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
/// drill's disc (which lies within the ball's outline).
fn drilled_volume(cx: f64, cy: f64, r: f64) -> f64 {
    over_disc(cx, cy, r, |x, y| 2.0 * (BALL * BALL - x * x - y * y).sqrt())
}

/// `height` integrated over the disc of radius `r` about `(cx, cy)`, in the
/// disc's own polar coordinates, where the integrand is smooth: Simpson
/// along each spoke.
fn over_disc(cx: f64, cy: f64, r: f64, height: impl Fn(f64, f64) -> f64) -> f64 {
    let (rings, spokes) = (2000_u32, 2000_u32);
    let ds = r / f64::from(rings);
    let dpsi = core::f64::consts::TAU / f64::from(spokes);
    let mut total = 0.0;
    for j in 0..spokes {
        let psi = dpsi * f64::from(j);
        let (c, s) = (psi.cos(), psi.sin());
        let mut line = 0.0;
        for i in 0..=rings {
            let t = ds * f64::from(i);
            let f = height(cx + t * c, cy + t * s) * t;
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
#[ignore = "heavy"]
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

/// A body whose top is a dome of the ball about the origin, its pole on
/// the z axis.
#[derive(Clone, Copy, Debug)]
enum Dome {
    /// The ball above its equator.
    Half,
    /// The ball above `z = 3`, less than half of it.
    Cap,
    /// The half ball standing on a drum of its radius down to `z = -10`.
    OnDrum,
}

impl Dome {
    /// Where the dome is cut from the ball.
    fn floor(self) -> f64 {
        match self {
            Self::Half | Self::OnDrum => 0.0,
            Self::Cap => 3.0,
        }
    }

    fn build(self, model: &mut Model) -> Result<Shape, String> {
        let ball = make_sphere(model, Frame::WORLD, BALL, T).unwrap().shape;
        let below = Frame::new(
            Point::new(-20.0, -20.0, self.floor() - 20.0),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        let block = make_box(model, below, (40.0, 40.0, 20.0), T).unwrap().shape;
        let dome = ogeom::boolean::cut(model, &ball, &block, T)
            .map_err(|e| format!("dome: {e}"))?
            .shape;
        if !matches!(self, Self::OnDrum) {
            return Ok(dome);
        }
        let foot = Frame::new(Point::new(0.0, 0.0, -10.0), Direction::Z, Direction::X, T).unwrap();
        let drum = make_cylinder(model, foot, BALL, 10.0, T).unwrap().shape;
        Ok(ogeom::boolean::fuse(model, &dome, &drum, T)
            .map_err(|e| format!("drum: {e}"))?
            .shape)
    }

    /// The body's volume.
    fn whole(self) -> f64 {
        let pi = core::f64::consts::PI;
        let cap = BALL - self.floor();
        let dome = pi * cap * cap * (3.0 * BALL - cap) / 3.0;
        match self {
            Self::Half | Self::Cap => dome,
            Self::OnDrum => dome + pi * BALL * BALL * 10.0,
        }
    }

    /// The body's volume inside a vertical drill of radius `r` standing at
    /// `(cx, cy)` from `z = -5` up past the pole, integrated over the
    /// drill's disc (which lies inside the dome's rim).
    fn inside(self, cx: f64, cy: f64, r: f64) -> f64 {
        let floor = self.floor();
        let above = over_disc(cx, cy, r, |x, y| {
            (BALL * BALL - x * x - y * y).sqrt() - floor
        });
        match self {
            Self::Half | Self::Cap => above,
            Self::OnDrum => above + core::f64::consts::PI * r * r * 5.0,
        }
    }
}

/// Cut and common of a dome with a vertical drill of radius `r` whose axis
/// stands at `offset`, against the body's volume and the common volume
/// integrated independently.
fn dome_drilled(dome: Dome, offset: Vector, r: f64) -> Result<(), String> {
    let mut model = Model::new();
    let body = dome.build(&mut model)?;
    let at = Frame::new(
        Point::new(offset.x, offset.y, -5.0),
        Direction::Z,
        Direction::Y,
        T,
    )
    .unwrap();
    let drill = make_cylinder(&mut model, at, r, 20.0, T).unwrap().shape;
    let mut volumes = Vec::new();
    for (name, made) in [
        ("cut", ogeom::boolean::cut(&mut model, &body, &drill, T)),
        (
            "common",
            ogeom::boolean::common(&mut model, &body, &drill, T),
        ),
    ] {
        let made = made.map_err(|e| format!("{name}: {e}"))?;
        let diagnosis = check(&model, &made.shape, T).unwrap();
        if !diagnosis.is_valid() {
            return Err(format!("{name}: {diagnosis}"));
        }
        volumes.push(volume(&model, &made.shape));
    }
    let (cut, common) = (volumes[0], volumes[1]);
    let (whole, inside) = (dome.whole(), dome.inside(offset.x, offset.y, r));
    if ((cut + common) - whole).abs() > whole * 1e-8 {
        return Err(format!("cut {cut} + common {common} against {whole}"));
    }
    if (common - inside).abs() > inside * 1e-6 {
        return Err(format!("common {common} against {inside}"));
    }
    Ok(())
}

#[test]
fn drill_through_a_dome_pole() {
    let r = 2.664_527_471_227_642;
    let mut failures = Vec::new();
    for dome in [Dome::Half, Dome::Cap, Dome::OnDrum] {
        if let Err(e) = dome_drilled(dome, Vector::new(0.0, -r, 0.0), r) {
            failures.push(format!("{dome:?}: {e}"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
#[ignore = "heavy"]
fn drills_through_a_dome_pole() {
    let r = 2.664_527_471_227_642;
    let mut failures = Vec::new();
    for (dx, dy) in [
        (0.0, -1.0),
        (1.0, 0.0),
        (
            core::f64::consts::FRAC_1_SQRT_2,
            core::f64::consts::FRAC_1_SQRT_2,
        ),
        (0.6, -0.8),
    ] {
        for radius in [r, 1.0, 3.5] {
            for wall in [
                radius,
                radius * (1.0 + 1e-7),
                radius * (1.0 - 1e-7),
                radius * 1.001,
            ] {
                let offset = Vector::new(dx * wall, dy * wall, 0.0);
                for dome in [Dome::Half, Dome::Cap, Dome::OnDrum] {
                    if let Err(e) = dome_drilled(dome, offset, radius) {
                        failures.push(format!("{dome:?} ({dx}, {dy}) r {radius} at {wall}: {e}"));
                    }
                }
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
