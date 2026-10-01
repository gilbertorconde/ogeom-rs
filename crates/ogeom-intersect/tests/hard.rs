#![allow(
    clippy::unwrap_used,
    reason = "test code; a failed unwrap is a failed test"
)]

//! Configurations the surface/surface literature names as hard.
//!
//! The corpus problem in miniature: every other test's inputs were chosen by
//! the people who wrote the code, and these were not; they are transcribed
//! from the published record of what breaks intersectors. Two of them earned
//! their keep immediately: the cone case exposed `Cone::distance_to` measuring
//! one nappe of a two-nappe surface (a defect in the *instrument*, flagged by
//! a correctly traced curve), and the tangent-torus case pinned the noise a
//! tangency-along-a-curve produces.

use ogeom_core::Tolerances;
use ogeom_geom::{
    ConeSurface, Curve3d as _, CylinderSurface, PlaneSurface, SphereSurface, SurfaceGeometry,
    TorusSurface,
};
mod support;

use ogeom_intersect::{Marching, Meeting, Stopped, branches, surface_surface};
use ogeom_math::{Cone, Cylinder, Direction, Frame, Plane, Point, Sphere, Torus, Vector};
use support::coverage::coverage;

const T: Tolerances = Tolerances::millimetres();

fn frame(origin: Point, z: Vector) -> Frame {
    Frame::new(
        origin,
        Direction::new(z, T).unwrap(),
        Direction::from_cross(z, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
        T,
    )
    .unwrap()
}

fn cyl(origin: Point, axis: Vector, r: f64) -> SurfaceGeometry {
    CylinderSurface::new(
        Cylinder::new(frame(origin, axis), r, T).unwrap(),
        (-6.0, 6.0),
    )
    .unwrap()
    .into()
}

fn sph(c: Point, r: f64) -> SurfaceGeometry {
    SphereSurface::new(Sphere::centred(c, r, T).unwrap()).into()
}

fn pln(o: Point, n: Vector) -> SurfaceGeometry {
    PlaneSurface::over(
        Plane::through(o, Direction::new(n, T).unwrap()),
        (-8.0, 8.0),
        (-8.0, 8.0),
    )
    .unwrap()
    .into()
}

fn off(s: &SurfaceGeometry, p: Point) -> f64 {
    match s {
        SurfaceGeometry::Plane(x) => x.plane().distance_to(p),
        SurfaceGeometry::Sphere(x) => x.sphere().distance_to(p),
        SurfaceGeometry::Cylinder(x) => x.cylinder().distance_to(p),
        SurfaceGeometry::Torus(x) => x.torus().distance_to(p),
        SurfaceGeometry::Cone(x) => x.cone().distance_to(p),
        _ => 0.0,
    }
}

fn worst(a: &SurfaceGeometry, b: &SurfaceGeometry, found: &[ogeom_intersect::Traced]) -> f64 {
    found
        .iter()
        .flat_map(|br| br.points.iter())
        .map(|p| off(a, *p).abs().max(off(b, *p).abs()))
        .fold(0.0_f64, f64::max)
}

fn options() -> Marching {
    Marching {
        chord: 1e-5,
        ..Marching::default()
    }
}

#[test]
fn equal_cylinders_at_a_shallow_angle_cover_everything_even_fragmented() {
    // The classic: equal radii force the two intersection curves through two
    // tangency points, and a two-degree crossing makes them long and thin.
    // The tracer stalls at the tangencies rather than jumping branch (the
    // deliberate refusal), so the answer arrives as fragments. What must hold
    // even so: every fragment on both surfaces, and nothing missed.
    let a = cyl(Point::ORIGIN, Vector::Z, 1.0);
    let tilt = 2.0_f64.to_radians();
    let b = cyl(Point::ORIGIN, Vector::new(tilt.sin(), 0.0, tilt.cos()), 1.0);

    let found = branches(&a, &b, options(), T).unwrap();
    assert!(!found.is_empty());
    assert!(
        found.iter().all(|br| br.stopped == Stopped::Stalled),
        "each fragment ends at a tangency it refuses to march through"
    );
    assert!(worst(&a, &b, &found) < 1e-7);
    let score = coverage(&a, &b, &found, 50, T).unwrap();
    assert!(
        score.complete(),
        "{}/{}: fragmented is acceptable, incomplete is not",
        score.covered,
        score.crossings
    );
}

#[test]
fn a_plane_tangent_along_a_ruling_is_the_analytic_paths_case() {
    // Tangential contact along a line. The marcher correctly finds no
    // transversal crossing; the analytic layer names the contact line. The
    // division of labour is the answer here, and both halves are pinned.
    let drum = cyl(Point::ORIGIN, Vector::Z, 2.0);
    let touching = pln(Point::new(2.0, 0.0, 0.0), Vector::X);

    assert!(branches(&drum, &touching, options(), T).unwrap().is_empty());
    assert!(matches!(
        surface_surface(&drum, &touching, T).unwrap(),
        Meeting::Along(ref c) if c.len() == 1
    ));
}

#[test]
fn near_tangent_sphere_and_cylinder_give_one_thin_loop() {
    // Equal radii, axis offset a thousandth: the near-tangential quartic. One
    // closed loop, thin in z, complete.
    let a = cyl(Point::new(1e-3, 0.0, 0.0), Vector::Z, 2.0);
    let b = sph(Point::ORIGIN, 2.0);
    let found = branches(&a, &b, options(), T).unwrap();
    assert_eq!(found.len(), 1);
    assert!(found[0].closed());
    assert!(worst(&a, &b, &found) < 1e-7);
    assert!(coverage(&a, &b, &found, 50, T).unwrap().complete());
}

#[test]
fn a_sphere_across_a_cones_apex_cuts_both_nappes() {
    // The case that caught the instrument. A cone's height range crosses its
    // apex, so the surface has two nappes, and a sphere spanning the apex
    // cuts a loop in each. The first run flagged the lower loop as 0.9 off
    // the cone; the trace was exact and Cone::distance_to was measuring one
    // nappe of a two-nappe surface.
    let cone: SurfaceGeometry = ConeSurface::new(
        Cone::new(Frame::WORLD, 0.5, 0.5_f64.atan(), T).unwrap(),
        (-3.0, 3.0),
    )
    .unwrap()
    .into();
    let ball = sph(Point::new(0.0, 0.0, -0.9), 1.0);

    let found = branches(&cone, &ball, options(), T).unwrap();
    assert_eq!(found.len(), 2, "one loop per nappe");
    for br in &found {
        assert!(br.closed());
    }
    // One loop above the apex plane z = -1, one below.
    let sides: Vec<f64> = found.iter().map(|br| br.points[0].z + 1.0).collect();
    assert!(
        sides[0] * sides[1] < 0.0,
        "both loops on one nappe: z offsets {sides:?}"
    );
    assert!(worst(&cone, &ball, &found) < 1e-7);
}

#[test]
fn tangency_along_a_circle_produces_fragments_not_a_curve() {
    // A plane resting on top of a torus touches along a whole circle. There
    // is no transversal curve to find, and the marcher cannot say "touching
    // along a curve": near the contact the two surfaces sit within the
    // correction's acceptance of each other, so seeds converge and wander
    // briefly before stalling. What comes back is fragments hugging the
    // contact circle: on both surfaces to rounding, describing nothing.
    //
    // Pinned as the documented limit it is. The honest answer needs
    // tangential contact traced as its own kind of curve, which is listed
    // as open in docs/PLAN.md.
    let torus: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, 3.0, 1.0, T).unwrap()).into();
    let resting = pln(Point::new(0.0, 0.0, 1.0), Vector::Z);

    let found = branches(&torus, &resting, options(), T).unwrap();
    for br in &found {
        assert_eq!(br.stopped, Stopped::Stalled, "fragments stall; none close");
        // Whatever comes back lies on both surfaces...
        for p in &br.points {
            assert!(off(&torus, *p).abs() < 1e-6);
            assert!(off(&resting, *p).abs() < 1e-6);
            // ...and hugs the contact circle at radius 3, z = 1.
            let radial = (p.x * p.x + p.y * p.y).sqrt();
            assert!((radial - 3.0).abs() < 0.1 && (p.z - 1.0).abs() < 0.01);
        }
    }
}

#[test]
fn the_same_tangency_asked_of_the_one_call_comes_back_as_contact() {
    // The fragments above are what the *crossing* marcher can say. Asked
    // through `intersect_surfaces`, the same plane on the same torus routes
    // those fragments into the tangential walker and comes back with the
    // contact itself: one curve, closed, marked as touching rather than
    // crossing, on the circle of radius 3 at z = 1.
    let torus: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, 3.0, 1.0, T).unwrap()).into();
    let resting = pln(Point::new(0.0, 0.0, 1.0), Vector::Z);

    let found = ogeom_intersect::intersect_surfaces(
        &torus,
        &resting,
        ogeom_intersect::IntersectOptions {
            tolerance: 1e-4,
            marching: options(),
        },
        T,
    )
    .unwrap();
    let ogeom_intersect::SurfaceIntersection::Along(curves) = found else {
        panic!("the surfaces meet along the contact: {found:?}");
    };
    assert_eq!(curves.len(), 1, "one contact, not one per fragment");
    let contact = &curves[0];
    assert!(contact.tangential, "they touch here, they do not cross");
    assert!(contact.closed, "the contact circle closes");
    let domain = contact.curve.domain();
    for k in 0..=16 {
        let t = domain.0 + (domain.1 - domain.0) * f64::from(k) / 16.0;
        let p = contact.curve.point_at(t, T).unwrap();
        assert!(
            (p.x.hypot(p.y) - 3.0).abs() < 1e-3 && (p.z - 1.0).abs() < 1e-3,
            "the contact is the circle of radius 3 at z = 1: {p:?}"
        );
    }
}

#[test]
fn a_ball_seated_in_a_torus_tube_has_its_contact_walked() {
    // The tangency with no closed form: a unit ball centred on the tube's
    // own centre line touches the torus along a whole tube cross-section.
    // Nothing analytic answers this pair, so it goes through the marcher
    // (which stalls, as tangencies make it), and the stalled fragments seed
    // the tangential walker instead of being discarded.
    //
    // The contact comes back in two arcs rather than one loop, and that is
    // the sphere's parameterization talking: the tube circle runs through
    // both of the ball's chart poles, where its derivatives degenerate and
    // the walk has nothing to step along. Two arcs meeting at the poles
    // cover the circle; the number of pieces is a fact about the chart, the
    // curve they lie on is a fact about the surfaces.
    let torus: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, 3.0, 1.0, T).unwrap()).into();
    let seated = sph(Point::new(3.0, 0.0, 0.0), 1.0);

    let found = ogeom_intersect::intersect_surfaces(
        &torus,
        &seated,
        ogeom_intersect::IntersectOptions {
            tolerance: 1e-4,
            marching: options(),
        },
        T,
    )
    .unwrap();
    let ogeom_intersect::SurfaceIntersection::Along(curves) = found else {
        panic!("the ball meets the tube along the contact: {found:?}");
    };
    assert!(!curves.is_empty(), "the contact was not walked at all");
    assert!(
        curves.iter().all(|c| c.tangential && !c.exact),
        "every piece is marched contact, not a crossing"
    );
    let centre = Point::new(3.0, 0.0, 0.0);
    let mut length = 0.0;
    for c in &curves {
        let domain = c.curve.domain();
        let mut previous = None;
        for k in 0..=64 {
            let t = domain.0 + (domain.1 - domain.0) * f64::from(k) / 64.0;
            let p = c.curve.point_at(t, T).unwrap();
            assert!(
                (p.distance(centre) - 1.0).abs() < 1e-3 && p.y.abs() < 1e-3,
                "the contact is the tube circle at the ball's radius: {p:?}"
            );
            if let Some(q) = previous {
                length += p.distance(q);
            }
            previous = Some(p);
        }
    }
    let circle = 2.0 * core::f64::consts::PI;
    assert!(
        (length - circle).abs() < circle * 0.02,
        "the arcs together are the whole circle: {length}"
    );
}

#[test]
fn a_plane_through_a_torus_tube_cuts_two_loops() {
    // The transversal cousin of the tangent case, and the first marched torus
    // result: a plane through the tube at half the minor radius cuts two
    // closed loops, one around the outer half, one around the inner.
    let torus: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, 3.0, 1.0, T).unwrap()).into();
    let cut = pln(Point::new(0.0, 0.0, 0.5), Vector::Z);

    let found = branches(&torus, &cut, options(), T).unwrap();
    assert_eq!(found.len(), 2);
    for br in &found {
        assert!(br.closed());
    }
    assert!(worst(&torus, &cut, &found) < 1e-7);
    // Distinct radii: one loop outside the tube's crown, one inside.
    let mut radii: Vec<f64> = found
        .iter()
        .map(|br| {
            let p = br.points[0];
            (p.x * p.x + p.y * p.y).sqrt()
        })
        .collect();
    radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert!(radii[0] < 3.0 && radii[1] > 3.0, "radii {radii:?}");
}

/// A plane through a torus's axis meets it in two whole tube circles, in
/// closed form, each starting on the outer equator where the torus's own
/// meridians do; a plane holding the torus's seam traces the seam itself.
#[test]
fn a_plane_through_a_torus_axis_cuts_two_meridians_exactly() {
    let torus: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, 3.0, 1.0, T).unwrap()).into();
    let seam_plane = pln(Point::ORIGIN, Vector::Y);
    let Meeting::Along(circles) = surface_surface(&torus, &seam_plane, T).unwrap() else {
        panic!("two circles");
    };
    assert_eq!(circles.len(), 2);
    let mut starts = Vec::new();
    for circle in &circles {
        let (a, b) = circle.domain();
        for k in 0..=16 {
            let p = circle
                .point_at(a + (b - a) * f64::from(k) / 16.0, T)
                .unwrap();
            let tube = (p.x.hypot(p.y) - 3.0).hypot(p.z);
            assert!((tube - 1.0).abs() < 1e-12 && p.y.abs() < 1e-12, "{p:?}");
        }
        starts.push(circle.point_at(a, T).unwrap());
    }
    // Each starts on the outer equator: four and minus four along x.
    assert!(
        starts
            .iter()
            .any(|p| p.distance(Point::new(4.0, 0.0, 0.0)) < 1e-12)
    );
    assert!(
        starts
            .iter()
            .any(|p| p.distance(Point::new(-4.0, 0.0, 0.0)) < 1e-12)
    );
}

#[test]
fn a_sphere_strictly_inside_another_is_apart_however_close() {
    let outer = sph(Point::ORIGIN, 2.0);
    let inner = sph(Point::new(0.999, 0.0, 0.0), 1.0);
    assert!(branches(&outer, &inner, options(), T).unwrap().is_empty());
    assert_eq!(surface_surface(&outer, &inner, T).unwrap(), Meeting::Apart);
}

/// Two planes through one point, a microradian apart, meet in a line. At
/// a hundred millimetres they stand a tenth of a micron apart, a thousand
/// times the confusion, so calling them one plane would be wrong.
#[test]
fn planes_a_microradian_apart_meet_in_a_line() {
    let a = 1e-6_f64;
    let plane = |n: Vector| -> SurfaceGeometry {
        PlaneSurface::new(Plane::through(Point::ORIGIN, Direction::new(n, T).unwrap())).into()
    };
    let flat = plane(Vector::new(0.0, 0.0, 1.0));
    let tilted = plane(Vector::new(0.0, a.sin(), a.cos()));
    assert!(matches!(
        surface_surface(&flat, &tilted, T).unwrap(),
        Meeting::Along(_)
    ));
    assert!(matches!(
        surface_surface(&flat, &flat, T).unwrap(),
        Meeting::Same
    ));
}

/// A line touching a sphere: the discriminant is zero, and rounding leaves
/// it a few ulps either side. Every scale and offset finds the touch.
#[test]
fn a_line_touching_a_sphere_is_found_at_every_scale() {
    use ogeom_geom::{Curve, LineCurve};
    use ogeom_intersect::{CurveSurfaceOptions, intersect_curve_surface};
    for radius in [1e-3, 0.7, 1.0, 7.3, 123.456, 1e3] {
        for shift in [0.0, 0.1, 1.0 / 3.0, 2.9] {
            let centre = Point::new(shift * radius, -shift, 0.3 * radius);
            let ball: SurfaceGeometry =
                SphereSurface::new(Sphere::new(frame(centre, Vector::Z), radius, T).unwrap())
                    .into();
            let touch = centre + Vector::new(0.0, radius, 0.0);
            let line: Curve = LineCurve::over(
                ogeom_math::Axis::new(touch - Vector::X * (3.0 * radius), Direction::X),
                0.0,
                6.0 * radius,
            )
            .unwrap()
            .into();
            let met =
                intersect_curve_surface(&line, &ball, CurveSurfaceOptions::default(), T).unwrap();
            assert!(!met.crossings.is_empty(), "radius {radius}, shift {shift}");
        }
    }
}

/// An oblique pipe tee: each closed section winds once round the branch
/// pipe, a loop in space whose chart image on the branch moves a period.
/// Its curve lies on both pipes within the tolerance it states, and that
/// tolerance is a fit's, not a pipe's size.
#[test]
fn an_oblique_pipe_tee_s_sections_lie_on_both_pipes() {
    use ogeom_intersect::{IntersectOptions, SurfaceIntersection, intersect_surfaces};
    let pipe = |origin: Point, axis: Vector, r: f64, reach: f64| -> SurfaceGeometry {
        CylinderSurface::new(
            Cylinder::new(frame(origin, axis), r, T).unwrap(),
            (-reach, reach),
        )
        .unwrap()
        .into()
    };
    for (ra, rb, tilt, offset) in [
        (10.0, 5.0, 1.0_f64, 0.0),
        (10.0, 7.0, 0.7, 1.0),
        (25.0, 12.5, 1.0, 0.0),
        (5.0, 3.0, 0.7, 0.0),
    ] {
        let centre = Point::new(100.0, 200.0, 300.0);
        let a = pipe(centre, Vector::Z, ra, 4.0 * ra);
        let b = pipe(
            centre + Vector::new(0.0, offset, 0.0),
            Vector::new(tilt.cos(), 0.0, tilt.sin()),
            rb,
            4.0 * ra,
        );
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&a, &b, IntersectOptions::default(), T).unwrap()
        else {
            panic!("the pipes cross");
        };
        assert_eq!(curves.len(), 2, "{ra} {rb} {tilt} {offset}");
        for c in &curves {
            assert!(
                c.tolerance < 1e-3,
                "{ra} {rb} {tilt} {offset}: {}",
                c.tolerance
            );
            let (lo, hi) = c.curve.domain();
            for k in 0..=1000 {
                let p = c
                    .curve
                    .point_at(lo + (hi - lo) * f64::from(k) / 1000.0, T)
                    .unwrap();
                let miss = off(&a, p).abs().max(off(&b, p).abs());
                assert!(
                    miss <= c.tolerance,
                    "{ra} {rb} {tilt} {offset}: {miss} off, stated {}",
                    c.tolerance
                );
            }
        }
    }
}

/// A spindle torus (its tube swallows the axis) meets a plane square to
/// its axis, a coaxial cylinder and a coaxial ring torus also on its folded
/// half past the axis: every parallel comes back, on both surfaces, with
/// its image in each chart.
#[test]
fn a_spindle_torus_s_folded_half_is_sectioned() {
    use ogeom_intersect::{IntersectOptions, SurfaceIntersection, intersect_surfaces};
    let spindle: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, 10.0, 15.0, T).unwrap()).into();
    let ring: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, 4.0, 3.0, T).unwrap()).into();
    let drum: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, 3.0, T).unwrap(), (-50.0, 50.0))
            .unwrap()
            .into();
    let level: SurfaceGeometry = PlaneSurface::over(
        Plane::through(Point::new(0.0, 0.0, 5.0), Direction::Z),
        (-50.0, 50.0),
        (-50.0, 50.0),
    )
    .unwrap()
    .into();
    for (other, count) in [(level, 2), (drum, 4), (ring, 2)] {
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&spindle, &other, IntersectOptions::default(), T).unwrap()
        else {
            panic!("the surfaces meet");
        };
        assert_eq!(curves.len(), count);
        for c in &curves {
            assert!(c.on_a.is_some() && c.on_b.is_some());
            let (lo, hi) = c.curve.domain();
            for k in 0..=16 {
                let p = c
                    .curve
                    .point_at(lo + (hi - lo) * f64::from(k) / 16.0, T)
                    .unwrap();
                assert!(off(&spindle, p).abs() < 1e-9 && off(&other, p).abs() < 1e-9);
            }
        }
    }
}

/// A cone whose window runs through its apex meets a coaxial cylinder, a
/// parallel-sided coaxial cone and a level plane also on its far nappe,
/// exactly and with its image in each chart; one stopping short of the
/// apex keeps only its own nappe's section.
#[test]
fn a_cone_through_its_apex_is_sectioned_on_both_nappes() {
    use ogeom_intersect::{IntersectOptions, SurfaceIntersection, intersect_surfaces};
    let quarter = core::f64::consts::FRAC_PI_4;
    let cone = |r: f64, window: (f64, f64)| -> SurfaceGeometry {
        ConeSurface::new(Cone::new(Frame::WORLD, r, quarter, T).unwrap(), window)
            .unwrap()
            .into()
    };
    let drum: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, 10.0, T).unwrap(), (-60.0, 60.0))
            .unwrap()
            .into();
    let level: SurfaceGeometry = PlaneSurface::over(
        Plane::through(Point::new(0.0, 0.0, -30.0), Direction::Z),
        (-50.0, 50.0),
        (-50.0, 50.0),
    )
    .unwrap()
    .into();
    // The apex stands at z = -20.
    let through = cone(20.0, (-50.0, 50.0));
    let short = cone(20.0, (-15.0, 50.0));
    for (a, b, heights) in [
        (&through, &drum, vec![-10.0, -30.0]),
        (&short, &drum, vec![-10.0]),
        (&through, &cone(40.0, (-100.0, 50.0)), vec![-30.0]),
        (&through, &level, vec![-30.0]),
    ] {
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(a, b, IntersectOptions::default(), T).unwrap()
        else {
            panic!("the surfaces meet");
        };
        let mut found: Vec<f64> = curves
            .iter()
            .map(|c| {
                assert!(c.exact && c.on_a.is_some() && c.on_b.is_some());
                let (lo, _) = c.curve.domain();
                let p = c.curve.point_at(lo, T).unwrap();
                assert!(off(a, p).abs() < 1e-9 && off(b, p).abs() < 1e-9);
                p.z
            })
            .collect();
        found.sort_by(|x, y| y.total_cmp(x));
        assert_eq!(found.len(), heights.len(), "{found:?}");
        for (got, want) in found.iter().zip(&heights) {
            assert!((got - want).abs() < 1e-9, "{found:?}");
        }
    }
}
