//! Drills along a ball's axis whose wall passes through both poles of its
//! chart, or a hair beside them, and drills through the pole of a dome: a
//! ball cut below its pole, alone or standing on a drum. A half ball
//! charted about a tilted axis stands on a drum or is drilled. A placed
//! half ball or ball is drilled through its pole.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_box, make_cylinder, make_sphere, volume_properties};
use ogeom::core::Tolerances;
use ogeom::math::{Axis, Direction, Frame, Point, Transform, Vector};
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
        // Each edge states how far its charts stand from its curve.
        for edge in
            ogeom::topo::explore_unique(&model, &made.shape, ogeom::topo::ShapeType::Edge).unwrap()
        {
            let stated = model
                .node(&edge)
                .and_then(|n| n.data().as_edge())
                .map_or(0.0, |d| d.tolerance.get());
            let gap = ogeom::algo::edge_pcurve_gap(&model, &edge, T)
                .unwrap()
                .unwrap_or(0.0);
            if gap > stated {
                return Err(format!(
                    "{name}: a chart {gap} off an edge stating {stated}"
                ));
            }
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

/// The half ball above `z = 0` of a ball charted about `axis`, its seam
/// leaving the axis toward `x_reference`.
fn tilted_half(model: &mut Model, axis: Vector, x_reference: Direction) -> Shape {
    let z = Direction::new(axis, T).unwrap();
    let chart = Frame::new(Point::ORIGIN, z, x_reference, T).unwrap();
    let ball = make_sphere(model, chart, BALL, T).unwrap().shape;
    let below = Frame::new(
        Point::new(-20.0, -20.0, -20.0),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let block = make_box(model, below, (40.0, 40.0, 20.0), T).unwrap().shape;
    ogeom::boolean::cut(model, &ball, &block, T).unwrap().shape
}

/// Fuse, cut and common of a half ball on a tilted chart with an upright
/// drum, each valid: fuse and common add up to both bodies, cut and common
/// to the half ball, and common matches `inside`.
fn tilted_with_drum(
    axis: Vector,
    x_reference: Direction,
    foot: Point,
    radius: f64,
    height: f64,
    inside: f64,
) -> Result<(), String> {
    let mut model = Model::new();
    let half = tilted_half(&mut model, axis, x_reference);
    let at = Frame::new(foot, Direction::Z, Direction::X, T).unwrap();
    let drum = make_cylinder(&mut model, at, radius, height, T)
        .unwrap()
        .shape;
    let mut volumes = Vec::new();
    for (name, made) in [
        ("fuse", ogeom::boolean::fuse(&mut model, &half, &drum, T)),
        ("cut", ogeom::boolean::cut(&mut model, &half, &drum, T)),
        (
            "common",
            ogeom::boolean::common(&mut model, &half, &drum, T),
        ),
    ] {
        let made = made.map_err(|e| format!("{name}: {e}"))?;
        let diagnosis = check(&model, &made.shape, T).unwrap();
        if !diagnosis.is_valid() {
            return Err(format!("{name}: {diagnosis}"));
        }
        volumes.push(volume(&model, &made.shape));
    }
    // The operands measured as the results are, so the sums compare like
    // with like; each against its closed form at the measure's own reach.
    let pi = core::f64::consts::PI;
    let (a, b) = (volume(&model, &half), volume(&model, &drum));
    for (name, measured, exact) in [
        ("half ball", a, 2.0 / 3.0 * pi * BALL.powi(3)),
        ("drum", b, pi * radius * radius * height),
    ] {
        if (measured - exact).abs() > exact * 1e-6 {
            return Err(format!("{name} {measured} against {exact}"));
        }
    }
    let (fuse, cut, common) = (volumes[0], volumes[1], volumes[2]);
    if ((fuse + common) - (a + b)).abs() > (a + b) * 1e-8 {
        return Err(format!("fuse {fuse} + common {common} against {}", a + b));
    }
    if ((cut + common) - a).abs() > a * 1e-8 {
        return Err(format!("cut {cut} + common {common} against {a}"));
    }
    if (common - inside).abs() > a * 1e-6 {
        return Err(format!("common {common} against {inside}"));
    }
    Ok(())
}

/// A drum of the ball's radius under the flat face shares the face and its
/// rim with the half ball, and no volume. On a tilted chart the rim is cut
/// into arcs where the chart's seam and poles meet it, and one of them runs
/// across the start of the drum's whole circle; that arc lies along the
/// circle on both sides of its start. A drum through the flat face and out
/// of the dome checks the same charts against the drilled volume.
#[test]
fn a_tilted_half_ball_on_a_drum_or_drilled() {
    let charts = [
        (Vector::new(1.0, 1.0, 1.0), Direction::Z),
        (Vector::new(1.0, 1.0, -1.0), Direction::Z),
        (Vector::new(1.0, 1.0, 0.0), Direction::X),
        (Vector::new(0.3, 0.2, 1.0), Direction::Z),
    ];
    let mut failures = Vec::new();
    for (axis, x_reference) in charts {
        let foot = Point::new(0.0, 0.0, -10.0);
        if let Err(e) = tilted_with_drum(axis, x_reference, foot, BALL, 10.0, 0.0) {
            failures.push(format!("{axis:?} {x_reference:?} under: {e}"));
        }
    }
    let r = 2.5;
    let (cx, cy) = (1.0, -2.0);
    let drilled = over_disc(cx, cy, r, |x, y| (BALL * BALL - x * x - y * y).sqrt());
    for (axis, x_reference) in &charts[1..] {
        let foot = Point::new(cx, cy, -5.0);
        if let Err(e) = tilted_with_drum(*axis, *x_reference, foot, r, 20.0, drilled) {
            failures.push(format!("{axis:?} {x_reference:?} through: {e}"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// A half ball charted about a level axis with its seam meridian level too,
/// so the seam runs along the flat face's rim: the rim's arcs lie on circles
/// whose frames are turned against the world axes. Drums under the face,
/// narrower, as wide and wider than it and off its centre, touch it over
/// the face and share no volume, so the common is empty and fuse and cut
/// add up. A drum through the flat face checks the drilled volume.
#[test]
fn a_half_ball_seamed_in_its_flat_face_on_a_drum() {
    let (axis, x_reference) = (Vector::new(1.0, 1.0, 0.0), Direction::Y);
    let mut failures = Vec::new();
    for (radius, cx, cy) in [
        (9.5, 0.0, 0.0),
        (BALL, 0.0, 0.0),
        (10.5, 0.0, 0.0),
        (BALL, 1.5, -0.5),
    ] {
        let foot = Point::new(cx, cy, -10.0);
        if let Err(e) = tilted_with_drum(axis, x_reference, foot, radius, 10.0, 0.0) {
            failures.push(format!("under, radius {radius} at ({cx}, {cy}): {e}"));
        }
    }
    let (r, cx, cy) = (2.5, 1.0, -2.0);
    let drilled = over_disc(cx, cy, r, |x, y| (BALL * BALL - x * x - y * y).sqrt());
    if let Err(e) = tilted_with_drum(
        axis,
        x_reference,
        Point::new(cx, cy, -5.0),
        r,
        20.0,
        drilled,
    ) {
        failures.push(format!("through: {e}"));
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

/// How a half ball or a ball is placed before it is drilled through its
/// pole.
#[derive(Clone, Copy, Debug)]
enum Placed {
    Moved,
    Turned,
    Mirrored,
    Halved,
    /// Restated on a left-handed sphere through the same points.
    LeftHanded,
}

impl Placed {
    /// The placement, `None` for the left-handed restatement.
    fn transform(self) -> Option<Transform> {
        match self {
            Self::Moved => Some(Transform::translation(Vector::new(5.0, -3.0, 7.0))),
            Self::Turned => Some(Transform::rotation(
                Axis::new(Point::ORIGIN, Direction::X),
                core::f64::consts::FRAC_PI_2,
            )),
            Self::Mirrored => Some(Transform::plane_mirror(Point::ORIGIN, Direction::Y)),
            Self::Halved => Some(Transform::scaling(Point::ORIGIN, 0.5, T).unwrap()),
            Self::LeftHanded => None,
        }
    }

    fn apply(self, model: &mut Model, shape: &Shape) -> Shape {
        match self.transform() {
            Some(t) => ogeom::algo::transformed(model, shape, t).unwrap().shape,
            None => left_handed(model, shape),
        }
    }
}

/// `shape` with each sphere restated on the left-handed frame through the
/// same points.
fn left_handed(model: &mut Model, shape: &Shape) -> Shape {
    use ogeom::core::OgeomResult;
    use ogeom::geom::{Curve, SphereSurface, SurfaceGeometry};
    let surface = |s: &SurfaceGeometry| -> OgeomResult<Option<(SurfaceGeometry, bool)>> {
        let SurfaceGeometry::Sphere(sphere) = s else {
            return Ok(None);
        };
        let f = sphere.sphere().frame();
        let left = Frame::from_axes(f.origin(), f.y(), f.x(), f.z(), T)?;
        let new: SurfaceGeometry =
            SphereSurface::new(ogeom::math::Sphere::new(left, sphere.sphere().radius(), T)?).into();
        let flip = ogeom::algo::normals_oppose(s, &new, T)?;
        Ok(Some((new, flip)))
    };
    let curve = |_: &Curve, _: (f64, f64)| -> OgeomResult<Option<(Curve, (f64, f64))>> { Ok(None) };
    ogeom::algo::restate_geometry(model, shape, &surface, &curve, T)
        .unwrap()
        .shape
}

/// Cut and common of a placed half ball or ball with a drill whose wall
/// runs through the pole of its chart, placed with it, against the closed
/// forms of the body's volume and the drilled volume.
fn placed_drilled(placed: Placed, whole_ball: bool) -> Result<(), String> {
    let r = 2.664_527_471_227_642;
    let mut model = Model::new();
    let body = if whole_ball {
        make_sphere(&mut model, Frame::WORLD, BALL, T)
            .unwrap()
            .shape
    } else {
        Dome::Half.build(&mut model)?
    };
    let body = placed.apply(&mut model, &body);
    let at = Frame::new(Point::new(0.0, -r, -20.0), Direction::Z, Direction::Y, T).unwrap();
    let drill = make_cylinder(&mut model, at, r, 40.0, T).unwrap().shape;
    let drill = match placed.transform() {
        Some(t) => {
            ogeom::algo::transformed(&mut model, &drill, t)
                .unwrap()
                .shape
        }
        None => drill,
    };
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
    let cube = placed
        .transform()
        .map_or(1.0, |t| t.scale_factor().abs().powi(3));
    let (whole, inside) = if whole_ball {
        let whole = 4.0 / 3.0 * core::f64::consts::PI * BALL.powi(3);
        (whole, drilled_volume(0.0, -r, r))
    } else {
        (Dome::Half.whole(), Dome::Half.inside(0.0, -r, r))
    };
    let (whole, inside) = (whole * cube, inside * cube);
    let (cut, common) = (volumes[0], volumes[1]);
    if ((cut + common) - whole).abs() > whole * 1e-6 {
        return Err(format!("cut {cut} + common {common} against {whole}"));
    }
    if (common - inside).abs() > inside * 1e-6 {
        return Err(format!("common {common} against {inside}"));
    }
    Ok(())
}

/// A half ball moved, turned, mirrored or halved by its placement, or on a
/// left-handed sphere, drilled through its pole, and a whole ball mirrored
/// or halved: each is recharted about a clear axis as an unplaced one is.
#[test]
fn a_placed_dome_or_ball_drilled_through_its_pole() {
    let mut failures = Vec::new();
    for placed in [
        Placed::Moved,
        Placed::Turned,
        Placed::Mirrored,
        Placed::Halved,
        Placed::LeftHanded,
    ] {
        if let Err(e) = placed_drilled(placed, false) {
            failures.push(format!("half ball {placed:?}: {e}"));
        }
    }
    for placed in [Placed::Mirrored, Placed::Halved] {
        if let Err(e) = placed_drilled(placed, true) {
            failures.push(format!("ball {placed:?}: {e}"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
