#![allow(
    clippy::unwrap_used,
    reason = "test code; a failed unwrap is a failed test"
)]

//! Configurations the surface/surface literature names as hard.
//!
//! The corpus problem in miniature: every other test's inputs were chosen by
//! the people who wrote the code, and these were not — they are transcribed
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
use ogeom_intersect::{Marching, Meeting, Stopped, branches, coverage, surface_surface};
use ogeom_math::{Cone, Cylinder, Direction, Frame, Plane, Point, Sphere, Torus, Vector};

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
    // The tracer stalls at the tangencies rather than jumping branch — the
    // deliberate refusal — so the answer arrives as fragments. What must hold
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
        "{}/{} — fragmented is acceptable, incomplete is not",
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
    // apex, so the surface has two nappes — and a sphere spanning the apex
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
    // along a curve" — near the contact the two surfaces sit within the
    // correction's acceptance of each other, so seeds converge and wander
    // briefly before stalling. What comes back is fragments hugging the
    // contact circle: on both surfaces to rounding, describing nothing.
    //
    // Pinned as the documented limit it is. The honest answer needs
    // tangential contact traced as its own kind of curve, which is recorded
    // in SCOPE's deferred table.
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
    // Nothing analytic answers this pair, so it goes through the marcher —
    // which stalls, as tangencies make it — and the stalled fragments seed
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

#[test]
fn a_sphere_strictly_inside_another_is_apart_however_close() {
    let outer = sph(Point::ORIGIN, 2.0);
    let inner = sph(Point::new(0.999, 0.0, 0.0), 1.0);
    assert!(branches(&outer, &inner, options(), T).unwrap().is_empty());
    assert_eq!(surface_surface(&outer, &inner, T).unwrap(), Meeting::Apart);
}
