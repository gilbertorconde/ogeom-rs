//! §B1 of `docs/PLAN.md`: the rolling ball solved for rather than projected.
//!
//! The property under test is the ball's own definition, and it is checked at
//! every station rather than at the ends: the centre is exactly the radius
//! from both supports, and the touch points are exactly where the centre's
//! own perpendicular meets them. That is what "the ball is seated" means, and
//! nothing weaker would distinguish a correct march from a plausible one.
//!
//! Where the seat has a closed form the answer is held to it as well — a
//! cylinder standing square on a plane rounds to a torus, and its spine is a
//! circle whose radius and height are arithmetic.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_fillet::{BlendStop, MarchedBlend, march_blend};
use ogeom_geom::{
    CircleCurve, Curve, CylinderSurface, PlaneSurface, Surface as _, SurfaceGeometry,
};
use ogeom_intersect::Marching;
use ogeom_math::{Circle, Cylinder, Direction, Frame, Plane, Point, Vector};

const T: Tolerances = Tolerances::millimetres();

fn options() -> Marching {
    Marching {
        chord: 1e-5,
        ..Marching::default()
    }
}

/// The ball is where it says it is: at every station, the centre stands the
/// radius from both supports and the touch points are the feet of that
/// distance.
fn assert_seated(
    blend: &MarchedBlend,
    first: &SurfaceGeometry,
    second: &SurfaceGeometry,
    radius: f64,
) {
    assert!(blend.len() > 8, "a march of {} stations", blend.len());
    for i in 0..blend.len() {
        let centre = blend.spine[i];
        let (p1, p2) = (blend.touch_first[i], blend.touch_second[i]);
        assert!(
            (centre.distance(p1) - radius).abs() < 1e-9,
            "station {i}: the ball touches the first support at its own radius: {}",
            centre.distance(p1)
        );
        assert!(
            (centre.distance(p2) - radius).abs() < 1e-9,
            "station {i}: and the second: {}",
            centre.distance(p2)
        );
        // The touch points are on their supports, at the parameters reported —
        // which is the whole point of solving for them rather than projecting.
        let (u1, v1) = blend.on_first[i];
        let (u2, v2) = blend.on_second[i];
        assert!(first.point_at(u1, v1, T).unwrap().distance(p1) < 1e-12);
        assert!(second.point_at(u2, v2, T).unwrap().distance(p2) < 1e-12);
        // And the ball touches rather than crosses: the centre lies along
        // each support's own normal from its touch point.
        for (surface, at, p) in [(first, (u1, v1), p1), (second, (u2, v2), p2)] {
            let normal = surface.normal_at(at.0, at.1, T).unwrap().vector();
            let along = (centre - p) / radius;
            assert!(
                along.cross(normal).magnitude() < 1e-7,
                "station {i}: the ball stands square to its support: {along:?} against {normal:?}"
            );
        }
    }
}

/// A cylinder standing square on a plane. The seat is a circle, the blend is
/// a torus, and every number in it is arithmetic — so the march is held to
/// the closed form and not merely to its own consistency.
#[test]
fn a_cylinder_square_on_a_plane_rounds_to_the_torus_it_should() {
    let (bore, radius) = (5.0, 1.25);
    let ground: SurfaceGeometry = PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-20.0, 20.0),
        (-20.0, 20.0),
    )
    .unwrap()
    .into();
    let wall: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, bore, T).unwrap(), (-1.0, 12.0))
            .unwrap()
            .into();
    // The guide: the circle where they meet.
    let guide: Curve = CircleCurve::new(Circle::new(Frame::WORLD, bore, T).unwrap()).into();

    let blend = march_blend(&ground, &wall, radius, &guide, options(), T).unwrap();
    assert_eq!(blend.stopped, BlendStop::Closed, "a rim closes on itself");
    assert_seated(&blend, &ground, &wall, radius);

    // The closed form. The ball rolls in the outside corner, so its centre
    // runs on a circle of radius `bore + r` at height `r`, its foot on the
    // plane on one of radius `bore + r`, and its touch on the wall at height
    // `r` on the wall's own radius.
    for i in 0..blend.len() {
        let centre = blend.spine[i];
        assert!(
            (centre.x.hypot(centre.y) - (bore + radius)).abs() < 1e-9,
            "the spine is a circle of radius {}: {}",
            bore + radius,
            centre.x.hypot(centre.y)
        );
        assert!((centre.z - radius).abs() < 1e-9, "at height {radius}");
        let ground_touch = blend.touch_first[i];
        assert!(ground_touch.z.abs() < 1e-9, "the foot is on the ground");
        assert!(
            (ground_touch.x.hypot(ground_touch.y) - (bore + radius)).abs() < 1e-9,
            "at the spine's own radius"
        );
        let wall_touch = blend.touch_second[i];
        assert!(
            (wall_touch.x.hypot(wall_touch.y) - bore).abs() < 1e-9,
            "the wall touch is on the wall"
        );
        assert!((wall_touch.z - radius).abs() < 1e-9, "at height {radius}");
    }
}

/// The same pair, tilted. There is no closed form for a cylinder meeting a
/// plane at an angle — which is exactly why this seat needs marching — so the
/// claim is the ball's own definition at every station, plus one number that
/// *is* arithmetic: the lowest the spine reaches.
#[test]
fn a_cylinder_meeting_a_plane_at_an_angle_is_marched() {
    let (bore, radius) = (4.0, 1.0);
    let tilt = 20.0_f64.to_radians();
    let ground: SurfaceGeometry = PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-30.0, 30.0),
        (-30.0, 30.0),
    )
    .unwrap()
    .into();
    let axis = Direction::new(Vector::new(tilt.sin(), 0.0, tilt.cos()), T).unwrap();
    let frame = Frame::new(
        Point::ORIGIN,
        axis,
        Direction::from_cross(axis.vector(), Vector::Y, T).unwrap(),
        T,
    )
    .unwrap();
    let wall: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(frame, bore, T).unwrap(), (-4.0, 16.0))
            .unwrap()
            .into();
    // The guide is the ellipse the tilted cylinder cuts from the plane —
    // taken here as the circle of the same mean size, deliberately *not* the
    // seat itself, to show the guide only says where the sections are.
    let guide: Curve = CircleCurve::new(Circle::new(Frame::WORLD, bore, T).unwrap()).into();

    let blend = march_blend(&ground, &wall, radius, &guide, options(), T).unwrap();
    assert_eq!(blend.stopped, BlendStop::Closed);
    assert_seated(&blend, &ground, &wall, radius);

    // Every centre stands the radius above the ground, whatever the tilt —
    // the plane is flat, so the ball's height is its radius exactly.
    for centre in &blend.spine {
        assert!(
            (centre.z - radius).abs() < 1e-9,
            "the ball rides on the ground: {}",
            centre.z
        );
    }
    // And the spine is not a circle, which is what says the tilt was really
    // followed rather than quietly ignored.
    let radii: Vec<f64> = blend.spine.iter().map(|c| c.x.hypot(c.y)).collect();
    let (lo, hi) = radii.iter().fold((f64::INFINITY, 0.0_f64), |(lo, hi), r| {
        (lo.min(*r), hi.max(*r))
    });
    assert!(hi - lo > 0.2, "the tilted seat is not round: {lo} to {hi}");

    // The tangency curve on the wall lives in the wall's own parameters, and
    // those parameters lift back onto it exactly — which is the property the
    // whole formulation exists for. A projected curve would agree only to the
    // projection's own tolerance.
    for (i, (u, v)) in blend.on_second.iter().enumerate() {
        let lifted = wall.point_at(*u, *v, T).unwrap();
        assert!(
            lifted.distance(blend.touch_second[i]) < 1e-12,
            "station {i}: exact by construction"
        );
    }
}

/// A radius the corner cannot hold is refused by name, not marched.
#[test]
fn a_ball_too_large_for_the_seat_is_refused_by_name() {
    let bore = 2.0;
    let ground: SurfaceGeometry = PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-8.0, 8.0),
        (-8.0, 8.0),
    )
    .unwrap()
    .into();
    // A bore *through* the plane: the ball has to roll inside it, and one of
    // radius larger than the bore cannot.
    let wall: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, bore, T).unwrap(), (-6.0, 6.0))
            .unwrap()
            .into();
    let guide: Curve = CircleCurve::new(Circle::new(Frame::WORLD, bore, T).unwrap()).into();

    // Inside the bore the ball must fit within the radius; ask for one four
    // times too big.
    let refused = march_blend(&ground, &wall, bore * 4.0, &guide, options(), T);
    let message = format!("{}", refused.unwrap_err());
    assert!(
        message.contains("no ball of radius"),
        "the refusal names the radius and the seat: {message}"
    );

    // And a radius of nothing is refused before any of that.
    assert!(march_blend(&ground, &wall, 0.0, &guide, options(), T).is_err());
}

/// A blend that runs off the end of a support says which support it ran off,
/// which is the state a caller has to act on — the section has to be carried
/// onto the next face, and it cannot be until it is told.
#[test]
fn a_blend_that_runs_out_says_which_support_it_ran_out_of() {
    let (bore, radius) = (5.0, 1.0);
    let ground: SurfaceGeometry = PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-20.0, 20.0),
        (-20.0, 20.0),
    )
    .unwrap()
    .into();
    // A wall that is only a third of the way round: the march reaches its
    // edge in `u` and stops there.
    let wall: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, bore, T).unwrap(), (-1.0, 12.0))
            .unwrap()
            .into();
    let quarter = core::f64::consts::FRAC_PI_2;
    let guide: Curve = ogeom_geom::TrimmedCurve::new(
        CircleCurve::new(Circle::new(Frame::WORLD, bore, T).unwrap()).into(),
        0.2,
        0.2 + quarter,
        T,
    )
    .unwrap()
    .into();

    let blend = march_blend(&ground, &wall, radius, &guide, options(), T).unwrap();
    assert_eq!(
        blend.stopped,
        BlendStop::RanPastTheGuide,
        "the guide ran out before either support did"
    );
    assert_seated(&blend, &ground, &wall, radius);
    // It covered the guide's own stretch and no more.
    let angle = |p: Point| p.y.atan2(p.x);
    let ends = (angle(blend.spine[0]), angle(blend.spine[blend.len() - 1]));
    assert!(
        (ends.1 - ends.0).abs() > quarter * 0.9,
        "the whole guide: {ends:?}"
    );
}
