//! Branch-point stitching, pinned on the classical configuration: two
//! equal perpendicular cylinders meet in two ellipses that cross at two
//! branch points, and the marcher's honest fragments reassemble into
//! exactly those ellipses.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::{CylinderSurface, SurfaceGeometry};
use ogeom_intersect::{Marching, Stopped, branches};
use ogeom_math::{Cylinder, Direction, Frame, Point};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn equal_cylinders_come_back_as_two_closed_ellipses() {
    let a: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, 1.0, T).unwrap(), (-3.0, 3.0))
            .unwrap()
            .into();
    let along_x = Frame::new(Point::ORIGIN, Direction::X, Direction::Y, T).unwrap();
    let b: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(along_x, 1.0, T).unwrap(), (-3.0, 3.0))
            .unwrap()
            .into();

    let found = branches(&a, &b, Marching::default(), T).unwrap();
    assert_eq!(found.len(), 2, "two ellipses, not a pile of fragments");
    for branch in &found {
        assert_eq!(
            branch.stopped,
            Stopped::Closed,
            "each closes through the branch points"
        );
        let length: f64 = branch
            .points
            .windows(2)
            .map(|pair| pair[0].distance(pair[1]))
            .sum();
        // The √2-by-1 ellipse's perimeter, which has no closed form but
        // has a value: 7.6404 by series. The polyline inscribes it.
        assert!(
            (length - 7.6404).abs() < 5e-3,
            "an ellipse's worth of curve, measured: {length}"
        );
        // Every point lies on both cylinders.
        for p in &branch.points {
            let on_a = (p.x.hypot(p.y) - 1.0).abs();
            let on_b = (p.z.hypot(p.y) - 1.0).abs();
            assert!(on_a < 1e-5 && on_b < 1e-5, "off-surface point {p:?}");
        }
    }
}

#[test]
fn tangential_contact_traces_as_the_circle_it_is() {
    use ogeom_geom::SphereSurface;
    use ogeom_intersect::{Contact, trace_tangential};
    use ogeom_math::Sphere;

    // A unit sphere seated inside a unit cylinder: tangent along the whole
    // equator, no transversal crossing anywhere: the case the crossing
    // walker honestly refuses and the valley walker owns.
    let drum: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, 1.0, T).unwrap(), (-2.0, 2.0))
            .unwrap()
            .into();
    let ball: SurfaceGeometry =
        SphereSurface::new(Sphere::centred(Point::ORIGIN, 1.0, T).unwrap()).into();

    // Seed at the contact: cylinder (u=0, v=0) is (1,0,0); the sphere's
    // equator holds the same point.
    use ogeom_geom::Surface as _;
    let seed_b = {
        // Find the sphere parameters of (1,0,0) by projection: equator at
        // whatever convention the chart uses.
        let ((ua, ub), (va, vb)) = ball.domain();
        let mut best = ((ua, va), f64::INFINITY);
        for i in 0..=64 {
            for j in 0..=64 {
                let u = ua + (ub - ua) * f64::from(i) / 64.0;
                let v = va + (vb - va) * f64::from(j) / 64.0;
                if let Ok(p) = ball.point_at(u, v, T) {
                    let d = p.distance(Point::new(1.0, 0.0, 0.0));
                    if d < best.1 {
                        best = ((u, v), d);
                    }
                }
            }
        }
        best.0
    };
    let contact = Contact {
        on_a: (0.0, 0.0),
        on_b: seed_b,
        point: Point::new(1.0, 0.0, 0.0),
    };

    let traced = trace_tangential(&drum, &ball, contact, Marching::default(), T).unwrap();
    assert_eq!(traced.stopped, Stopped::Closed, "the equator closes");
    let length: f64 = traced
        .points
        .windows(2)
        .map(|pair| pair[0].distance(pair[1]))
        .sum();
    assert!(
        (length - core::f64::consts::TAU).abs() < 1e-2,
        "one full equator: {length}"
    );
    for p in &traced.points {
        assert!(
            (p.x.hypot(p.y) - 1.0).abs() < 1e-5 && p.z.abs() < 1e-5,
            "on the equator: {p:?}"
        );
    }
}

/// A profile of one quadratic span.
fn quadratic(control: [Point; 3]) -> ogeom_geom::Curve {
    use ogeom_math::KnotVector;
    ogeom_geom::BSplineCurve::new(
        KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap(),
        control.to_vec(),
        T,
    )
    .unwrap()
    .into()
}

/// Where the apex lines of [`walls_ending_on_their_apexes`] cross.
const APEXES_CROSS: Point = Point::new(-11.21, 5.64, -15.68);

/// Two walls swept square to each other from quadratic profiles, each
/// profile ending on its apex: one along x, its apex line y = 5.64 at
/// z = -15.68, one down z from z = -11.76, its apex line y = 5.64 at
/// x = -11.21, its profile reaching `width` along x before it turns down.
/// They touch where the apex lines cross, normals parallel, and the
/// section runs into that touch ever more shallowly until the walk gives
/// out short of it. `along_x` and `down_z` move each profile's middle
/// control point off level.
fn walls_ending_on_their_apexes(
    along_x: f64,
    down_z: f64,
    width: f64,
) -> (SurfaceGeometry, SurfaceGeometry) {
    use ogeom_geom::ExtrusionSurface;
    use ogeom_math::Vector;
    let a = ExtrusionSurface::new(
        quadratic([
            Point::new(-13.13, -3.29, 0.0),
            Point::new(-13.13, 5.64 + along_x, 0.0),
            Point::new(-13.13, 5.64, -15.68),
        ]),
        Direction::X,
        22.96,
    )
    .unwrap()
    .into();
    let b = ExtrusionSurface::new(
        quadratic([
            Point::new(-11.21, 5.64, -11.76),
            Point::new(-11.21 + width, 5.64 - down_z, -11.76),
            Point::new(-11.21 + width, 3.64, -11.76),
        ]),
        Direction::new(Vector::new(0.0, 0.0, -1.0), T).unwrap(),
        3.92,
    )
    .unwrap()
    .into();
    (a, b)
}

/// How far the nearer end of the one section found stands from the touch.
fn section_misses_the_touch(a: &SurfaceGeometry, b: &SurfaceGeometry) -> f64 {
    // A chord a boolean marches a part of this size at.
    let marching = Marching {
        chord: 1e-6,
        ..Marching::default()
    };
    let found = branches(a, b, marching, T).unwrap();
    assert_eq!(found.len(), 1, "one section from the touch to the top");
    let ends = [
        found[0].points[0],
        found[0].points[found[0].points.len() - 1],
    ];
    ends.iter()
        .map(|p| p.distance(APEXES_CROSS))
        .fold(f64::INFINITY, f64::min)
}

/// One profile's apex is drawn a nanometre past its end, as rounding leaves
/// it in a part read from a file: the touch lies a hair outside that
/// wall's chart, where the search for it must reach onto the bound and
/// back. The section ends on the touch, to the confusion distance either
/// way the apex is pushed.
#[test]
fn a_section_into_a_touch_where_both_walls_end_ends_on_it() {
    // A nanometre off the middle control point tilts the profile's end
    // off level by two in a billion, and the apex moves just past that
    // end: lowered, on the wall along x, whose profile ends on its apex,
    // and on the wall down z, whose profile starts on it.
    for (along_x, down_z) in [(-1e-9, 0.0), (0.0, 1e-9)] {
        let (a, b) = walls_ending_on_their_apexes(along_x, down_z, 2.0);
        let miss = section_misses_the_touch(&a, &b);
        assert!(
            miss <= T.confusion(),
            "the section stops {miss} short of the touch ({along_x}, {down_z})"
        );
    }
}

/// With the second wall's profile turning down within half a millimetre
/// the section runs straight into the touch in long steps, and the samples
/// cut away round the stall leave its end 2.7e-3 short: further than a
/// stall is looked for round a touch at this chord. The touch is solved
/// from that end, and the section still ends on it.
#[test]
fn a_section_ending_a_long_step_short_of_a_touch_ends_on_it() {
    let (a, b) = walls_ending_on_their_apexes(0.0, 0.0, 0.5);
    let miss = section_misses_the_touch(&a, &b);
    assert!(
        miss <= T.confusion(),
        "the section stops {miss} short of the touch"
    );
}

/// A saddle `z = x·y·(y − 0.1)` cut by the plane `z = 0`: the lines
/// `x = 0`, `y = 0` and `y = 0.1`, crossing where the saddle touches the
/// plane at `(0, 0)` and `(0, 0.1)`. The arm along `y = 0` from `x = 1.3`
/// comes back 5.2e-2 short of the origin, and the other touch, 0.11 off
/// that end, lies as near as a stall is looked for round a touch too. The
/// arm will not carry onto that one, and is carried onto the touch solved
/// from its own end instead: it ends on the origin.
#[test]
fn an_arm_stalling_between_two_touches_ends_on_the_one_it_runs_into() {
    use ogeom_geom::{BSplineSurface, PlaneSurface};
    use ogeom_math::{ControlGrid, KnotVector, Plane};
    let (d, (x0, x1), (y0, y1)) = (0.1, (-0.7, 1.3), (-0.6, 0.9));
    // Linear in x, quadratic in y: the control heights are the blossom of
    // x·y·(y − d) at each corner of the Bézier net.
    let blossom = |s: f64, t: f64| s * t - d * (s + t) / 2.0;
    let mut net = Vec::new();
    for x in [x0, x1] {
        for (y, (s, t)) in [
            (y0, (y0, y0)),
            (f64::midpoint(y0, y1), (y0, y1)),
            (y1, (y1, y1)),
        ] {
            net.push(Point::new(x, y, x * blossom(s, t)));
        }
    }
    let saddle: SurfaceGeometry = BSplineSurface::new(
        KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap(),
        KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).unwrap(),
        &ControlGrid::new(net, 2, 3).unwrap(),
        T,
    )
    .unwrap()
    .into();
    let plane: SurfaceGeometry = PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Z),
        (-2.0, 2.0),
        (-2.0, 2.0),
    )
    .unwrap()
    .into();

    let found = branches(&saddle, &plane, Marching::default(), T).unwrap();
    let arm = found
        .iter()
        .find(|branch| {
            branch
                .points
                .iter()
                .all(|p| p.y.abs() <= 1e-6 && p.x >= -T.confusion())
        })
        .expect("the arm along y = 0 on the side x > 0");
    let ends = [arm.points[0], arm.points[arm.points.len() - 1]];
    assert!(
        ends.iter()
            .any(|p| p.distance(Point::new(x1, 0.0, 0.0)) <= T.confusion()),
        "the arm starts on the saddle's edge: {ends:?}"
    );
    let miss = ends
        .iter()
        .map(|p| p.distance(Point::ORIGIN))
        .fold(f64::INFINITY, f64::min);
    assert!(
        miss <= T.confusion(),
        "the arm stops {miss} short of the touch at the origin"
    );
}
