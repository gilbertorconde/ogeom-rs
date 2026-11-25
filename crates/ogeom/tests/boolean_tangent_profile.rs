//! A fill fused beside a ball cut through its centre: the fill's end
//! profile, a fitted curve held loosely, ends on the ball's cut flat where
//! it leaves the ball's rim arc tangentially. That is a touch at the
//! profile's end, not a stretch of the profile along the arc.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom::algo::{
    check, make_box, make_edge_between, make_face, make_revolution, make_sphere, make_vertex,
    make_wire, surface_properties, volume_properties,
};
use ogeom::core::{Tolerance, Tolerances};
use ogeom::geom::{BSplineCurve, CircleCurve, Curve, Curve3d as _, LineCurve, PlaneSurface};
use ogeom::math::{Axis, Circle, Direction, Frame, KnotVector, Plane, Point, Weighted};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The ball's radius.
const BALL: f64 = 5.0;

/// The fill's radius: a circle tangent to the plate and to the ball.
const FILL: f64 = 1.0;

/// The tolerance the fill's end profiles state, as a fitted fill's do.
const LOOSE: f64 = 5e-3;

fn fine() -> Deflection {
    Deflection::with_chord(1e-4).unwrap()
}

fn point_of(model: &Model, vertex: &Shape) -> Point {
    match model.node(vertex).unwrap().data() {
        NodeData::Vertex(data) => data.point,
        _ => unreachable!("a vertex"),
    }
}

/// The fill's profile in the plane x = 0, on the side y < 0, in the
/// meridian's own terms: `r = -y` from the ball's axis and `z` up.
struct Profile {
    /// Where the fill's circle stands on the plate.
    foot: f64,
    /// The height at which it meets the ball.
    top: f64,
    /// The angle at the ball's centre from the plate to that meeting.
    meet: f64,
}

impl Profile {
    fn new() -> Self {
        let foot = BALL.mul_add(BALL, 2.0 * BALL * FILL).sqrt();
        Self {
            foot,
            top: BALL * FILL / (BALL + FILL),
            meet: FILL.atan2(foot),
        }
    }

    /// The meridian region's area and its first moment about the axis, in
    /// closed form: it runs from the ball's rim `r = sqrt(R² - z²)` out to
    /// the fill's circle `r = foot - sqrt(ρ² - (ρ - z)²)`, z from 0 to `top`.
    fn area_and_moment(&self) -> (f64, f64) {
        let (c, rho, h) = (self.foot, FILL, self.top);
        // The integrals of sqrt(ρ² - u²) and sqrt(R² - z²).
        let circle =
            |u: f64, r: f64| (u * r.mul_add(r, -u * u).sqrt() + r * r * (u / r).asin()) / 2.0;
        let under_fill = circle(rho, rho) - circle(rho - h, rho);
        let under_ball = circle(h, BALL) - circle(0.0, BALL);
        let area = c * h - under_fill - under_ball;
        let squares = (rho.powi(3) - (rho - h).powi(3)) / 3.0;
        let moment = (c.mul_add(c, rho * rho - BALL * BALL) * h - 2.0 * c * under_fill - squares
            + h.powi(3) / 3.0)
            / 2.0;
        (area, moment)
    }
}

/// A 20 x 20 x 2 plate with the quarter of the ball about the origin that
/// stands on it at x >= 0: the ball cut by x < 0 and by z < 0, so one pole
/// lies on the cut flat.
fn plate_and_quarter_ball(model: &mut Model) -> Shape {
    let ball = make_sphere(model, Frame::WORLD, BALL, T).unwrap().shape;
    let at = Frame::new(Point::new(-6.0, -6.0, -6.0), Direction::Z, Direction::X, T).unwrap();
    let back = make_box(model, at, (6.0, 12.0, 12.0), T).unwrap().shape;
    let below = make_box(model, at, (12.0, 12.0, 6.0), T).unwrap().shape;
    let q = ogeom::boolean::cut(model, &ball, &back, T).unwrap().shape;
    let q = ogeom::boolean::cut(model, &q, &below, T).unwrap().shape;
    let at = Frame::new(
        Point::new(-10.0, -10.0, -2.0),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let plate = make_box(model, at, (20.0, 20.0, 2.0), T).unwrap().shape;
    ogeom::boolean::fuse(model, &plate, &q, T).unwrap().shape
}

/// The fill round the ball's foot on the plate, on the side x >= 0: the
/// profile turned half way about the ball's axis, from the plane x = 0 at
/// y < 0 to the same plane at y > 0. The fill's circle is written as a
/// rational spline, so its end profiles are fitted curves, and held to
/// [`LOOSE`].
fn fill(model: &mut Model, profile: &Profile) -> Shape {
    let at = |r: f64, z: f64| Point::new(0.0, -r, z);
    let (c, h) = (profile.foot, profile.top);
    let rim = at(BALL, 0.0);
    let foot = at(c, 0.0);
    let meet = at(BALL * profile.meet.cos(), h);
    let [v_rim, v_foot, v_meet] = [rim, foot, meet].map(|p| make_vertex(model, p).shape);

    let line = Curve::Line(LineCurve::segment(rim, foot, T).unwrap());
    let range = line.domain();
    let plate_edge = make_edge_between(model, line, range, &v_rim, &v_foot, T)
        .unwrap()
        .shape;

    // The quarter circle's arc from the plate to the ball, as a rational
    // quadratic: its middle control point where the end tangents meet,
    // weighted by the cosine of half the turn.
    let centre = at(c, FILL);
    let (down, inward) = (foot - centre, meet - centre);
    let turn = (down.dot(inward) / (FILL * FILL)).acos();
    let bisector = (down + inward).normalized(T).unwrap();
    let corner = centre + bisector * (FILL / (turn / 2.0).cos());
    let control = vec![
        Weighted::new(foot, 1.0, T).unwrap(),
        Weighted::new(corner, (turn / 2.0).cos(), T).unwrap(),
        Weighted::new(meet, 1.0, T).unwrap(),
    ];
    let knots = KnotVector::clamped_uniform(2, 3).unwrap();
    let spline = Curve::BSpline(BSplineCurve::rational(knots, control).unwrap());
    let range = spline.domain();
    let fill_edge = make_edge_between(model, spline, range, &v_foot, &v_meet, T)
        .unwrap()
        .shape;

    // The ball's rim in the plane x = 0, from the plate up to the meeting.
    let frame = Frame::new(Point::ORIGIN, -Direction::X, -Direction::Y, T).unwrap();
    let rim_curve = Curve::Circle(CircleCurve::new(Circle::new(frame, BALL, T).unwrap()));
    let rim_edge = make_edge_between(model, rim_curve, (0.0, profile.meet), &v_rim, &v_meet, T)
        .unwrap()
        .shape;

    let wire = make_wire(model, &[plate_edge, fill_edge, rim_edge.reversed()], T)
        .unwrap()
        .shape;
    let plane = PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::X),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = make_face(model, plane.into(), &[wire], T).unwrap().shape;
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Z,
    };
    let solid = make_revolution(model, &face, axis, PI, T).unwrap().shape;

    // The two end profiles of the fill's circle: on the plane x = 0, from
    // the plate at the fill's foot up to the meeting with the ball.
    let mut loosened = 0;
    for edge in explore_unique(model, &solid, ShapeType::Edge).unwrap() {
        let ends: Vec<Point> = explore_unique(model, &edge, ShapeType::Vertex)
            .unwrap()
            .iter()
            .map(|v| point_of(model, v))
            .collect();
        let on_foot = |p: &Point| p.z.abs() < 1e-9 && (p.y.abs() - c).abs() < 1e-9;
        let on_meet = |p: &Point| (p.z - h).abs() < 1e-9;
        if ends.iter().all(|p| p.x.abs() < 1e-9)
            && ends.iter().any(on_foot)
            && ends.iter().any(on_meet)
        {
            model.widen(&edge, Tolerance::new(LOOSE).unwrap()).unwrap();
            loosened += 1;
        }
    }
    assert_eq!(loosened, 2);
    solid
}

/// The fill touches the plate and the ball and overlaps neither, so the
/// fuse holds both whole: its volume is the two volumes summed, and its
/// area theirs less twice the faces they share. The plate is 800 and its
/// quarter ball (2/3)πR³/2; the fill's half turn is π times its profile's
/// first moment about the axis (Pappus). The fill's faces are its two end
/// profiles, the half annulus under it on the plate, the band of the ball
/// it covers (π R z over a half turn to height z) and its own half torus,
/// π ρ (c (π/2 - θ) - ρ (1 - sin θ)) for the arc turning through π/2 - θ
/// about its centre at distance c from the axis.
#[test]
fn a_fill_whose_end_profile_leaves_the_ball_tangentially_fuses_beside_it() {
    let mut model = Model::with_tolerances(T);
    let profile = Profile::new();
    let base = plate_and_quarter_ball(&mut model);
    let fill = fill(&mut model, &profile);

    let (profile_area, moment) = profile.area_and_moment();
    let (c, rho, h, meet) = (profile.foot, FILL, profile.top, profile.meet);
    let base_volume = 800.0 + PI * BALL.powi(3) / 3.0;
    let fill_volume = PI * moment;
    let base_area = 960.0 + PI * BALL * BALL;
    let on_plate = PI / 2.0 * c.mul_add(c, -BALL * BALL);
    let on_ball = PI * BALL * h;
    let torus = PI * rho * c.mul_add(PI / 2.0 - meet, -rho * (1.0 - meet.sin()));
    let fill_area = 2.0f64.mul_add(profile_area, on_plate + on_ball + torus);

    for (name, made) in [
        ("fuse", ogeom::boolean::fuse(&mut model, &base, &fill, T)),
        (
            "fuse_seated",
            ogeom::boolean::fuse_seated(&mut model, &base, &fill, T),
        ),
    ] {
        let made = made.unwrap_or_else(|e| panic!("{name}: {e}")).shape;
        let diagnosis = check(&model, &made, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
        let volume = volume_properties(&model, &made, fine(), T).unwrap().mass;
        let area = surface_properties(&model, &made, fine(), T).unwrap().mass;
        let want_volume = base_volume + fill_volume;
        let want_area = 2.0f64.mul_add(-(on_plate + on_ball), base_area + fill_area);
        // Both are measured on the exact faces, to rounding; a split that
        // moved an edge by the loose profiles' width would show far above it.
        assert!(
            (volume - want_volume).abs() <= want_volume * 1e-9,
            "{name}: volume {volume} against {want_volume}"
        );
        assert!(
            (area - want_area).abs() <= want_area * 1e-9,
            "{name}: area {area} against {want_area}"
        );
    }
}
