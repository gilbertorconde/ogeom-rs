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

/// A ball and a torus whose tube passes a thousandth and a half of the ball's
/// radius beside one of its poles. The section there was marched across
/// the chart's collapse: snapped onto the pole on the ball and left as a
/// run of points on the torus, a sliver lost from every result and the
/// fuse's ball face turned into the material. The ball is charted about an
/// axis whose poles stand clear, and the four results add up.
#[test]
fn a_torus_tube_beside_a_balls_pole() {
    let mut model = Model::new();
    let frame = |at: Point, axis: Direction, x: Direction| Frame::new(at, axis, x, T).unwrap();
    let ball = make_sphere(
        &mut model,
        frame(
            Point::new(
                -1.457_204_608_907_108,
                -0.799_339_000_187_928_6,
                -2.165_425_354_226_132,
            ),
            Direction::X,
            Direction::Y,
        ),
        3.25,
        T,
    )
    .unwrap()
    .shape;
    let ring = ogeom::algo::make_torus(
        &mut model,
        frame(
            Point::new(
                -2.229_968_783_447_556_6,
                -1.501_831_531_557_040_5,
                -1.392_870_007_057_077_9,
            ),
            Direction::Y,
            Direction::X,
        ),
        3.255_997_530_269_826_5,
        0.969_621_186_363_192_2,
        T,
    )
    .unwrap()
    .shape;
    let (va, vb) = (volume(&model, &ball), volume(&model, &ring));
    let mut made = Vec::new();
    for name in ["fuse", "common", "cut", "cut back"] {
        let built = match name {
            "fuse" => ogeom::boolean::fuse(&mut model, &ball, &ring, T),
            "common" => ogeom::boolean::common(&mut model, &ball, &ring, T),
            "cut" => ogeom::boolean::cut(&mut model, &ball, &ring, T),
            _ => ogeom::boolean::cut(&mut model, &ring, &ball, T),
        }
        .unwrap()
        .shape;
        let diagnosis = check(&model, &built, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
        made.push(volume(&model, &built));
    }
    let [fuse, common, cut, back] = made[..] else {
        unreachable!()
    };
    let slack = (va + vb) * 1e-7;
    assert!(
        (fuse + common - va - vb).abs() <= slack,
        "fuse {fuse} common {common}"
    );
    assert!(
        (cut + common - va).abs() <= slack,
        "cut {cut} common {common} ball {va}"
    );
    assert!(
        (back + common - vb).abs() <= slack,
        "back {back} common {common} torus {vb}"
    );
}

/// A ball and a drum on oblique axes whose common keeps a sliver of the
/// ball a few hundredths wide, two thirds of a millimetre from its pole.
/// The sliver's rings are polylines straight in the ball's chart, which
/// stretches there, and a ray meeting the ball just outside the sliver read
/// as crossing it: a point outside the common read inside, and the check
/// called the common's ball face turned in. A ray meeting a curved face
/// near its trim reads it against rings drawn finer, and the common checks
/// valid.
#[test]
fn a_sliver_of_a_ball_beside_its_pole_is_read_against_its_true_trim() {
    let mut model = Model::new();
    let frame = |at: [f64; 3], z: [f64; 3], x: [f64; 3]| {
        let direction = |v: [f64; 3]| Direction::new(Vector::new(v[0], v[1], v[2]), T).unwrap();
        Frame::new(
            Point::new(at[0], at[1], at[2]),
            direction(z),
            direction(x),
            T,
        )
        .unwrap()
    };
    let ball = make_sphere(
        &mut model,
        frame(
            [1.339_148_868_768_348_6, 2.5, -1.217_733_954_411_407_6],
            [
                -0.757_728_519_763_805_1,
                -0.584_303_982_865_410_6,
                0.290_579_328_143_229_45,
            ],
            [0.0, 0.445_284_593_575_738_1, 0.895_389_094_597_477_1],
        ),
        3.75,
        T,
    )
    .unwrap()
    .shape;
    let drum = make_cylinder(
        &mut model,
        frame(
            [-1.0, 1.926_468_907_191_658_5, -0.5],
            [
                0.817_392_615_291_016_8,
                -0.571_843_915_017_859_5,
                0.069_741_302_861_065_54,
            ],
            [
                -0.576_080_994_711_431,
                -0.811_380_686_961_981,
                0.098_954_880_412_197_19,
            ],
        ),
        1.449_072_030_036_451_2,
        5.5,
        T,
    )
    .unwrap()
    .shape;
    let common = ogeom::boolean::common(&mut model, &ball, &drum, T)
        .unwrap()
        .shape;
    let diagnosis = check(&model, &common, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
}

/// Two oblique boxes fused, cut by a ball one of whose faces passes beside
/// its chart's pole, not through it: that face cuts a small circle round
/// the pole, whose image no chart about that pole holds, so the ball is
/// charted about an axis clear of it. The cut and the common add up to the
/// fused pair.
#[test]
fn a_plane_beside_a_ball_s_pole_turns_its_chart() {
    let mut model = Model::new();
    let frame = |at: (f64, f64, f64), z: (f64, f64, f64), x: (f64, f64, f64)| {
        Frame::new(
            Point::new(at.0, at.1, at.2),
            Direction::new(Vector::new(z.0, z.1, z.2), T).unwrap(),
            Direction::new(Vector::new(x.0, x.1, x.2), T).unwrap(),
            T,
        )
        .unwrap()
    };
    let first = make_box(
        &mut model,
        frame(
            (2.0, -1.843_077_895_014_359_6, 0.5),
            (
                -0.542_500_285_810_833_2,
                -0.703_634_176_884_512_1,
                0.458_903_459_362_881_07,
            ),
            (
                -0.840_055_617_143_986_7,
                0.454_400_559_053_330_6,
                -0.296_355_685_008_498_4,
            ),
        ),
        (5.452_203_210_302_122, 5.5, 6.5),
        T,
    )
    .unwrap()
    .shape;
    let second = make_box(
        &mut model,
        frame(
            (2.136_387_783_977_939_5, 1.5, -3.0),
            (
                -0.231_694_730_394_031_1,
                0.185_443_820_116_384_33,
                -0.954_949_287_390_843_3,
            ),
            (0.0, -0.981_661_732_134_095_8, -0.190_631_171_804_316_57),
        ),
        (4.197_667_522_462_558, 7.5, 7.597_456_914_006_61),
        T,
    )
    .unwrap()
    .shape;
    let fused = ogeom::boolean::fuse(&mut model, &first, &second, T)
        .unwrap()
        .shape;
    let ball = make_sphere(
        &mut model,
        frame(
            (2.821_697_342_574_099, 0.0, -0.5),
            (
                0.185_798_440_566_491_9,
                -0.503_200_703_192_765_2,
                0.843_959_709_814_020_9,
            ),
            (0.0, 0.858_915_246_217_360_6, 0.512_117_759_714_863_7),
        ),
        3.0,
        T,
    )
    .unwrap()
    .shape;
    let cut = ogeom::boolean::cut(&mut model, &fused, &ball, T)
        .unwrap()
        .shape;
    let common = ogeom::boolean::common(&mut model, &fused, &ball, T)
        .unwrap()
        .shape;
    for result in [&cut, &common] {
        assert!(check(&model, result, T).unwrap().is_valid());
    }
    let (whole, c, k) = (
        volume(&model, &fused),
        volume(&model, &common),
        volume(&model, &cut),
    );
    assert!(
        (c + k - whole).abs() < whole * 1e-6,
        "{c} + {k} against {whole}"
    );
}

/// A sector of a cone of radius 3 and height 4 turning `turn` about its
/// axis, apex down at the origin, its curved face bounded by the rim arc
/// and the two rulings to the apex and by no edge at the apex: the rulings
/// meet there at a vertex, as a file that states no degenerate edge leaves
/// them.
fn cone_sector(model: &mut Model, turn: f64) -> Shape {
    use ogeom::geom::{BSpline2d, CircleCurve, ConeSurface, Curve, LineCurve, PlaneSurface};
    use ogeom::math::{Circle, Cone, KnotVector, Plane, Point2};
    use ogeom::topo::Location;
    let (radius, height): (f64, f64) = (3.0, 4.0);
    let slant = radius.hypot(height);
    let cone = Cone::new(Frame::WORLD, 0.0, (radius / height).atan(), T).unwrap();
    let surface = model
        .geometry_mut()
        .add_surface(ConeSurface::new(cone, (0.0, height)).unwrap().into());
    let rim = |u: f64| Point::new(radius * u.cos(), radius * u.sin(), height);
    let (apex, centre) = (Point::ORIGIN, Point::new(0.0, 0.0, height));
    let vertex = |model: &mut Model, p: Point| ogeom::algo::make_vertex(model, p).shape;
    let (a, o) = (vertex(model, apex), vertex(model, centre));
    let (b0, b1) = (vertex(model, rim(0.0)), vertex(model, rim(turn)));
    let line = |model: &mut Model, from: (&Shape, Point), to: (&Shape, Point)| {
        let curve: Curve = LineCurve::segment(from.1, to.1, T).unwrap().into();
        let length = from.1.distance(to.1);
        ogeom::algo::make_edge_between(model, curve, (0.0, length), from.0, to.0, T)
            .unwrap()
            .shape
    };
    let rim_circle = Circle::new(
        Frame::new(centre, Direction::Z, Direction::X, T).unwrap(),
        radius,
        T,
    )
    .unwrap();
    let arc = ogeom::algo::make_edge_between(
        model,
        CircleCurve::new(rim_circle).into(),
        (0.0, turn),
        &b0,
        &b1,
        T,
    )
    .unwrap()
    .shape;
    let r0 = line(model, (&b0, rim(0.0)), (&a, apex));
    let r1 = line(model, (&b1, rim(turn)), (&a, apex));
    let l0 = line(model, (&o, centre), (&b0, rim(0.0)));
    let l1 = line(model, (&o, centre), (&b1, rim(turn)));
    let axis = line(model, (&o, centre), (&a, apex));
    // The cone's trims: the arc along the rim row, each ruling down its
    // column from the rim row to the apex row, paced as its line is.
    let here = Location::identity();
    let row: ogeom::geom::PlanarCurve =
        ogeom::geom::Line2d::segment(Point2::new(0.0, height), Point2::new(turn, height), T)
            .unwrap()
            .into();
    ogeom::algo::attach_pcurve(model, &arc, row, surface, here.clone(), (0.0, turn)).unwrap();
    for (edge, u) in [(&r0, 0.0), (&r1, turn)] {
        let column = BSpline2d::new(
            KnotVector::new(vec![0.0, 0.0, slant, slant], 1).unwrap(),
            vec![Point2::new(u, height), Point2::new(u, 0.0)],
            T,
        )
        .unwrap();
        ogeom::algo::attach_pcurve(
            model,
            edge,
            column.into(),
            surface,
            here.clone(),
            (0.0, slant),
        )
        .unwrap();
    }
    let wire = ogeom::algo::make_wire(model, &[arc.reversed(), r0.clone(), r1.reversed()], T)
        .unwrap()
        .shape;
    let curved = ogeom::algo::make_face_on(model, surface, &[wire], T)
        .unwrap()
        .shape;
    let plane = |origin: Point, normal: Direction, x: Direction| -> ogeom::geom::SurfaceGeometry {
        PlaneSurface::new(Plane::new(Frame::new(origin, normal, x, T).unwrap())).into()
    };
    let flat = |model: &mut Model, on: ogeom::geom::SurfaceGeometry, ring: Vec<Shape>| {
        ogeom::algo::make_face_with_pcurves(model, on, &[ring], T)
            .unwrap()
            .shape
    };
    let cap = flat(
        model,
        plane(centre, Direction::Z, Direction::X),
        vec![l0.clone(), arc, l1.reversed()],
    );
    let wall0 = flat(
        model,
        plane(apex, -Direction::Y, Direction::X),
        vec![r0.reversed(), l0.reversed(), axis.clone()],
    );
    let outward = Direction::new(Vector::new(-turn.sin(), turn.cos(), 0.0), T).unwrap();
    let wall1 = flat(
        model,
        plane(apex, outward, Direction::Z),
        vec![l1, r1, axis.reversed()],
    );
    let shell = ogeom::algo::make_shell(model, &[curved, cap, wall0, wall1])
        .unwrap()
        .shape;
    ogeom::algo::make_solid(model, &[shell]).unwrap().shape
}

/// The cone sector, three quarters of a turn, cut and kept by a slab
/// below its rim holding the apex and by a band clear of the apex. The
/// curved face's ring reaches the apex along its two rulings with no edge
/// there, so its chart stands open along the apex row between them: the
/// boolean closes it as a degenerate edge would. Each result checks, and
/// measures what the sector holds below height h, 9π (h / 4)³.
#[test]
fn a_cone_sector_with_no_apex_edge_is_cut_through_its_apex_and_beside_it() {
    let pi = core::f64::consts::PI;
    let below = |h: f64| 9.0 * pi * (h / 4.0).powi(3);
    for (low, high) in [(-1.0, 2.0), (1.0, 3.0)] {
        let mut model = Model::new();
        let sector = cone_sector(&mut model, 1.5 * pi);
        assert!(check(&model, &sector, T).unwrap().is_valid());
        let tool = make_box(
            &mut model,
            Frame::new(Point::new(-5.0, -5.0, low), Direction::Z, Direction::X, T).unwrap(),
            (10.0, 10.0, high - low),
            T,
        )
        .unwrap()
        .shape;
        let inside = below(high.min(4.0)) - below(low.max(0.0));
        for (name, kept, want) in [
            (
                "cut",
                ogeom::boolean::cut(&mut model, &sector, &tool, T),
                below(4.0) - inside,
            ),
            (
                "common",
                ogeom::boolean::common(&mut model, &sector, &tool, T),
                inside,
            ),
        ] {
            let kept = kept
                .unwrap_or_else(|e| panic!("slab {low}..{high} {name}: {e}"))
                .shape;
            let diagnosis = check(&model, &kept, T).unwrap();
            assert!(
                diagnosis.is_valid(),
                "slab {low}..{high} {name}: {diagnosis}"
            );
            let got = volume(&model, &kept);
            assert!(
                (got - want).abs() <= want * 1e-9,
                "slab {low}..{high} {name}: {got} against {want}"
            );
        }
    }
}

/// A quarter of a ball of radius 3, its chart's seam in the quarter cut
/// away, as the boolean leaves it, poles and all.
fn quarter_ball(model: &mut Model) -> Shape {
    let away = Direction::new(Vector::new(-1.0, -1.0, 0.0), T).unwrap();
    let frame = Frame::new(Point::ORIGIN, Direction::Z, away, T).unwrap();
    let ball = make_sphere(model, frame, 3.0, T).unwrap().shape;
    let quadrant = make_box(
        model,
        Frame::new(Point::new(0.0, 0.0, -5.0), Direction::Z, Direction::X, T).unwrap(),
        (5.0, 5.0, 10.0),
        T,
    )
    .unwrap()
    .shape;
    ogeom::boolean::common(model, &ball, &quadrant, T)
        .unwrap()
        .shape
}

/// The quarter ball with its curved face bounded by its two meridians
/// alone: the poles stated by no edge, as a file that states no degenerate
/// edge leaves them.
fn wedge_with_no_pole_edges(model: &mut Model) -> Shape {
    use ogeom::topo::{Orientation, ShapeType, explore_unique};
    let wedge = quarter_ball(model);
    let faces = explore_unique(model, &wedge, ShapeType::Face).unwrap();
    let mut kept = Vec::with_capacity(faces.len());
    for face in &faces {
        let stored = face.oriented(Orientation::Forward);
        let surface = model.node(face).unwrap().data().as_face().unwrap().surface;
        let mut wires = Vec::new();
        let mut stripped = false;
        for wire in model.children_of(&stored).unwrap() {
            let edges = model.ordered_children_of(&wire).unwrap();
            let solid_edges: Vec<Shape> = edges
                .iter()
                .filter(|e| !model.node(e).unwrap().data().as_edge().unwrap().degenerate)
                .cloned()
                .collect();
            stripped |= solid_edges.len() != edges.len();
            wires.push(
                ogeom::algo::make_wire(model, &solid_edges, T)
                    .unwrap()
                    .shape,
            );
        }
        if !stripped {
            kept.push(face.clone());
            continue;
        }
        let rebuilt = ogeom::algo::make_face_on(model, surface, &wires, T)
            .unwrap()
            .shape;
        kept.push(if face.orientation() == Orientation::Reversed {
            rebuilt.reversed()
        } else {
            rebuilt
        });
    }
    let shell = ogeom::algo::make_shell(model, &kept).unwrap().shape;
    ogeom::algo::make_solid(model, &[shell]).unwrap().shape
}

/// `solid` with one meridian of its ball face stated `turns` whole periods
/// along the chart's `u`.
fn meridian_moved(model: &mut Model, solid: &Shape, turns: f64) {
    use ogeom::geom::{Curve2d as _, Surface as _};
    use ogeom::topo::{EdgeRepr, NodeData, ShapeType, explore_unique};
    let face = explore_unique(model, solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let s = model.node(f).unwrap().data().as_face().unwrap().surface;
            model.geometry().surface(s).unwrap().is_closed_u(T)
                || model.geometry().surface(s).unwrap().is_periodic_u()
        })
        .unwrap();
    let surface = model.node(&face).unwrap().data().as_face().unwrap().surface;
    let ((u0, u1), _) = model.geometry().surface(surface).unwrap().domain();
    for edge in explore_unique(model, &face, ShapeType::Edge).unwrap() {
        let data = model.node(&edge).unwrap().data().as_edge().unwrap();
        let Some(EdgeRepr::PCurve {
            curve,
            range,
            location,
            ..
        }) = data.pcurve_for(surface, edge.location()).cloned()
        else {
            continue;
        };
        let stated = model.geometry().pcurve(curve).unwrap().clone();
        let (a, b) = (
            stated.point_at(range.0, T).unwrap(),
            stated.point_at(range.1, T).unwrap(),
        );
        // A meridian: one column, pole to pole.
        if (a.x - b.x).abs() > 1e-9 {
            continue;
        }
        let shift =
            ogeom::math::Transform2::translation(ogeom::math::Vector2::new(turns * (u1 - u0), 0.0));
        let moved = model
            .geometry_mut()
            .add_pcurve(stated.transformed(&shift, T).unwrap());
        if let Some(NodeData::Edge(data)) = model.node_mut(&edge).map(|n| n.data_mut()) {
            data.representations
                .retain(|r| r.surface() != Some(surface));
            data.representations.push(EdgeRepr::PCurve {
                curve: moved,
                range,
                surface,
                location,
            });
        }
        return;
    }
    panic!("the ball face has no meridian");
}

/// The quarter ball with no pole edges, exact and converted to splines,
/// as stated and with a meridian stated a period off either way, cut and
/// kept by half spaces through it. Each result is the quarter ball's own,
/// the poles closed along their rows, or a refusal: a meridian a period
/// off leaves the ring more than a period wide along the pole rows, which
/// no pole line closes. Never a solid that is not the quarter ball's.
#[test]
fn a_quarter_ball_with_no_pole_edges_cuts_as_stated_or_refuses() {
    let mut reference = Model::new();
    let stated = quarter_ball(&mut reference);
    // Half spaces y > 0.113 and x > 0.137, each through the quarter ball.
    let slabs = [(-10.0, 0.113), (0.137, -10.0)];
    let slab = |model: &mut Model, (x, y): (f64, f64)| {
        make_box(
            model,
            Frame::new(Point::new(x, y, -10.0), Direction::Z, Direction::X, T).unwrap(),
            (20.0, 20.0, 20.0),
            T,
        )
        .unwrap()
        .shape
    };
    let whole = volume(&reference, &stated);
    let mut want = Vec::new();
    for at in slabs {
        let tool = slab(&mut reference, at);
        for keep in [false, true] {
            let built = if keep {
                ogeom::boolean::common(&mut reference, &stated, &tool, T)
            } else {
                ogeom::boolean::cut(&mut reference, &stated, &tool, T)
            }
            .unwrap_or_else(|e| panic!("the stated quarter ball: {e}"))
            .shape;
            assert!(check(&reference, &built, T).unwrap().is_valid());
            want.push(volume(&reference, &built));
        }
    }
    for nurbs in [false, true] {
        for turns in [0.0, 1.0, -1.0] {
            let mut model = Model::new();
            let mut solid = wedge_with_no_pole_edges(&mut model);
            if nurbs {
                solid = ogeom::algo::to_nurbs(&mut model, &solid, T).unwrap().shape;
            }
            if turns != 0.0 {
                meridian_moved(&mut model, &solid, turns);
            }
            assert!(check(&model, &solid, T).unwrap().is_valid());
            // A chart stated a period off is measured by its mesh, as
            // closely as the solid itself is.
            let band = (volume(&model, &solid) - whole).abs();
            let mut k = 0;
            for at in slabs {
                let tool = slab(&mut model, at);
                for keep in [false, true] {
                    let built = if keep {
                        ogeom::boolean::common(&mut model, &solid, &tool, T)
                    } else {
                        ogeom::boolean::cut(&mut model, &solid, &tool, T)
                    };
                    let case = format!("nurbs {nurbs} turns {turns} slab {at:?} common {keep}");
                    match built {
                        Ok(built) => {
                            let diagnosis = check(&model, &built.shape, T).unwrap();
                            assert!(diagnosis.is_valid(), "{case}: {diagnosis}");
                            let got = volume(&model, &built.shape);
                            assert!(
                                (got - want[k]).abs() <= band + want[k] * 1e-6,
                                "{case}: {got} against {}",
                                want[k]
                            );
                        }
                        Err(e) => assert!(turns != 0.0, "{case}: {e}"),
                    }
                    k += 1;
                }
            }
        }
    }
}
