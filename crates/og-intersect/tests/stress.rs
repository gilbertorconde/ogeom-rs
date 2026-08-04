//! Adversarial cases for the marching intersector and its instruments.
//!
//! Each of these began as a probe that was expected to break something, and
//! three of them did. They are kept as assertions so the defects they found
//! stay found:
//!
//! - coincident surfaces returned six confident little curves that existed
//!   nowhere but in rounding — the tangency gate sat *below* the noise floor
//!   the Newton correction is allowed to leave;
//! - a walk that ran into a surface's edge reported `Stalled` rather than
//!   `LeftTheDomain`, because it converges on the boundary from inside and
//!   never enters the strict band the crossing test uses;
//! - the coverage instrument compared cell centres to polyline *vertices*, so
//!   a perfectly traced straight line — whose points sit far apart, since
//!   nothing bends — scored half missing.

#![allow(
    clippy::unwrap_used,
    reason = "test code; a failed unwrap is a failed test"
)]

use og_core::Tolerances;
use og_geom::{CylinderSurface, PlaneSurface, SphereSurface, SurfaceGeometry};
use og_intersect::{Marching, Meeting, Stopped, branches, coverage, surface_surface};
use og_math::{Cylinder, Direction, Frame, Plane, Point, Sphere, Vector};

const T: Tolerances = Tolerances::millimetres();

fn sphere(centre: Point, radius: f64) -> SurfaceGeometry {
    SphereSurface::new(Sphere::centred(centre, radius, T).unwrap()).into()
}

fn cylinder(origin: Point, axis: Vector, radius: f64) -> SurfaceGeometry {
    let frame = Frame::new(
        origin,
        Direction::new(axis, T).unwrap(),
        Direction::from_cross(axis, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
        T,
    )
    .unwrap();
    CylinderSurface::new(Cylinder::new(frame, radius, T).unwrap(), (-4.0, 4.0))
        .unwrap()
        .into()
}

fn plane(origin: Point, normal: Vector) -> SurfaceGeometry {
    PlaneSurface::over(
        Plane::through(origin, Direction::new(normal, T).unwrap()),
        (-6.0, 6.0),
        (-6.0, 6.0),
    )
    .unwrap()
    .into()
}

fn off(surface: &SurfaceGeometry, p: Point) -> f64 {
    match surface {
        SurfaceGeometry::Plane(x) => x.plane().distance_to(p),
        SurfaceGeometry::Sphere(x) => x.sphere().distance_to(p),
        SurfaceGeometry::Cylinder(x) => x.cylinder().distance_to(p),
        _ => 0.0,
    }
}

fn worst_deviation(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    found: &[og_intersect::Traced],
) -> f64 {
    found
        .iter()
        .flat_map(|branch| branch.points.iter())
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
fn coincident_surfaces_yield_no_curves_at_all() {
    // The seed correction may leave a residual up to the confusion tolerance,
    // and on coincident surfaces that residual masquerades as an angle between
    // the two normals. With the tangency gate below that noise floor, identical
    // spheres came back as six short curves lying perfectly on both surfaces
    // and describing nothing. The gate now sits above the floor, and the
    // honest answer for a coincident pair is: no curves, ask surface_surface,
    // which says Same.
    let s = sphere(Point::ORIGIN, 2.0);
    let found = branches(&s, &s.clone(), options(), T).unwrap();
    assert!(
        found.is_empty(),
        "coincident spheres produced {} phantom curves",
        found.len()
    );
    assert_eq!(surface_surface(&s, &s.clone(), T).unwrap(), Meeting::Same);

    let c = cylinder(Point::ORIGIN, Vector::Z, 1.0);
    assert!(branches(&c, &c.clone(), options(), T).unwrap().is_empty());
}

#[test]
fn tangential_contact_yields_no_curves_and_the_analytic_path_names_it() {
    // Touching is not crossing: there is no curve to find, and the marcher
    // finding nothing is agreement. The analytic path is what names the
    // contact point.
    let resting = sphere(Point::new(0.0, 0.0, 3.0), 3.0);
    let ground = plane(Point::ORIGIN, Vector::Z);
    assert!(
        branches(&resting, &ground, options(), T)
            .unwrap()
            .is_empty()
    );
    assert!(matches!(
        surface_surface(&ground, &resting, T).unwrap(),
        Meeting::Touching(ref p) if p.len() == 1
    ));

    let apart_by_zero = sphere(Point::new(4.0, 0.0, 0.0), 2.0);
    assert!(
        branches(&sphere(Point::ORIGIN, 2.0), &apart_by_zero, options(), T)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn a_walk_that_reaches_the_edge_of_a_surface_says_so() {
    // Two parallel cylinders overlap in two straight lines that run off both
    // ends of the domain. The walk converges on the boundary from inside, so it
    // stalls a fraction of a step short — which used to be reported as
    // `Stalled`, indistinguishable from a genuine singularity.
    let a = cylinder(Point::ORIGIN, Vector::Z, 1.0);
    let b = cylinder(Point::new(1.99, 0.0, 0.0), Vector::Z, 1.0);
    let found = branches(&a, &b, options(), T).unwrap();

    assert_eq!(found.len(), 2, "two lines, one per side of the lens");
    for branch in &found {
        assert_eq!(branch.stopped, Stopped::LeftTheDomain);
        // Full length: v spans the whole height on the first surface.
        let (lo, hi) = branch
            .on_a
            .iter()
            .fold((f64::MAX, f64::MIN), |(lo, hi), p| {
                (lo.min(p.1), hi.max(p.1))
            });
        assert!(lo < -3.99 && hi > 3.99, "the line stops short: {lo}..{hi}");
    }
    assert!(worst_deviation(&a, &b, &found) < 1e-7);

    let score = coverage(&a, &b, &found, 60, T).unwrap();
    assert!(
        score.complete(),
        "a fully traced pair of lines scored {}/{} — the instrument is \
         measuring its own sampling again",
        score.covered,
        score.crossings
    );
}

#[test]
fn coverage_measures_the_curve_not_the_spacing_of_its_points() {
    // A straight line's trace has widely spaced points, because nothing bends.
    // Comparing cell centres to vertices marked the cells between two samples
    // of a perfect trace as missed; the measure is against segments now.
    let a = cylinder(Point::ORIGIN, Vector::Z, 1.0);
    let b = plane(Point::ORIGIN, Vector::X);
    let found = branches(&a, &b, options(), T).unwrap();
    assert_eq!(found.len(), 2, "a plane through the axis cuts two lines");
    let score = coverage(&a, &b, &found, 60, T).unwrap();
    assert!(
        score.complete(),
        "{}/{} — cells between the sparse samples of a straight line were \
         counted as missed",
        score.covered,
        score.crossings
    );
}

#[test]
fn near_degenerate_crossings_still_trace_completely() {
    let cases: Vec<(&str, SurfaceGeometry, SurfaceGeometry, usize)> = vec![
        (
            "plane 1e-3 into a sphere",
            sphere(Point::ORIGIN, 3.0),
            plane(Point::new(0.0, 0.0, 2.999), Vector::Z),
            1,
        ),
        (
            "spheres barely overlapping",
            sphere(Point::ORIGIN, 3.0),
            sphere(Point::new(5.98, 0.0, 0.0), 3.0),
            1,
        ),
        (
            "nearly equal crossed cylinders",
            cylinder(Point::ORIGIN, Vector::Z, 1.0),
            cylinder(Point::ORIGIN, Vector::X, 1.001),
            2,
        ),
        (
            "cylinder grazing a sphere from inside",
            sphere(Point::ORIGIN, 3.0),
            cylinder(Point::ORIGIN, Vector::Z, 2.997),
            1,
        ),
    ];
    for (name, a, b, count) in cases {
        let found = branches(&a, &b, options(), T).unwrap();
        assert_eq!(found.len(), count, "{name}: wrong branch count");
        let worst = worst_deviation(&a, &b, &found);
        assert!(worst < 1e-7, "{name}: deviation {worst:e}");
        let score = coverage(&a, &b, &found, 40, T).unwrap();
        assert!(
            score.complete(),
            "{name}: {}/{} covered",
            score.covered,
            score.crossings
        );
    }
}

#[test]
fn far_from_the_origin_the_answer_is_the_same() {
    // Catastrophic-cancellation territory: the same crossed cylinders, a
    // million units out. The absolute coordinates eat twenty bits of the
    // mantissa and the answer should not care.
    let far = Point::new(1.0e6, 1.0e6, 1.0e6);
    let a = cylinder(far, Vector::Z, 1.0);
    let b = cylinder(far, Vector::X, 1.6);
    let found = branches(&a, &b, options(), T).unwrap();
    assert_eq!(found.len(), 2);
    for branch in &found {
        assert!(branch.closed());
    }
    assert!(worst_deviation(&a, &b, &found) < 1e-7);
    assert!(coverage(&a, &b, &found, 40, T).unwrap().complete());
}

#[test]
fn the_mutual_blindness_of_seeding_and_coverage_is_a_fact_not_a_surprise() {
    // Two parallel cylinders overlapping by microns meet in two lines a few
    // thousandths apart — thinner than a seeding cell *and* thinner than a
    // coverage cell at the same resolution. The seeder misses them and the
    // instrument cannot see that it missed, because both look through the same
    // grid. Pinned so the limitation stays documented behaviour rather than
    // becoming a discovery.
    let a = cylinder(Point::ORIGIN, Vector::Z, 1.0);
    let b = cylinder(Point::new(1.99999, 0.0, 0.0), Vector::Z, 1.0);
    let found = branches(&a, &b, options(), T).unwrap();
    let score = coverage(&a, &b, &found, 60, T).unwrap();
    assert!(
        found.is_empty() && score.crossings == 0,
        "if either side now resolves this, tighten this test to assert the \
         curves are found: {} branches, {} crossings",
        found.len(),
        score.crossings
    );
}
