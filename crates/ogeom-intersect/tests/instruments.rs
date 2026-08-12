//! The instruments held to their own claims.
//!
//! An instrument nobody checks is a number nobody should trust. Each of these
//! is a negative control: the accuracy measure is shown to report a real
//! deviation rather than zero, and the completeness measure is shown to fail
//! on an answer that is genuinely half missing. Without them the two
//! instruments would agree with the intersector by construction, which is no
//! agreement at all.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::print_stdout,
    reason = "test code"
)]

mod support;

/// Every point of every reported curve lies on both surfaces.
mod accuracy {
    use crate::support::benchmark::*;
    use ogeom_core::Tolerances;
    use ogeom_geom::{Curve, SurfaceGeometry};
    use ogeom_geom::{CylinderSurface, PlaneSurface, SphereSurface};
    use ogeom_intersect::{Meeting, surface_surface};
    use ogeom_math::{Cylinder, Direction, Frame, Plane, Point, Sphere, Vector};

    const T: Tolerances = Tolerances::millimetres();

    fn plane(origin: Point, normal: Vector) -> SurfaceGeometry {
        PlaneSurface::new(Plane::through(origin, Direction::new(normal, T).unwrap())).into()
    }

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
        CylinderSurface::new(Cylinder::new(frame, radius, T).unwrap(), (-10.0, 10.0))
            .unwrap()
            .into()
    }

    fn cone(
        origin: Point,
        axis: Vector,
        reference_radius: f64,
        half_angle: f64,
    ) -> SurfaceGeometry {
        let frame = Frame::new(
            origin,
            Direction::new(axis, T).unwrap(),
            Direction::from_cross(axis, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
            T,
        )
        .unwrap();
        ogeom_geom::ConeSurface::new(
            ogeom_math::Cone::new(frame, reference_radius, half_angle, T).unwrap(),
            (-10.0, 10.0),
        )
        .unwrap()
        .into()
    }

    /// Every case with a closed form, named.
    fn corpus() -> Vec<(String, SurfaceGeometry, SurfaceGeometry)> {
        let mut out = Vec::new();
        let mut add = |name: &str, a: SurfaceGeometry, b: SurfaceGeometry| {
            out.push((name.to_string(), a, b));
        };

        add(
            "plane/plane crossing",
            plane(Point::ORIGIN, Vector::Z),
            plane(Point::ORIGIN, Vector::X),
        );
        add(
            "plane/plane oblique",
            plane(Point::new(1.0, 2.0, 3.0), Vector::new(1.0, 1.0, 1.0)),
            plane(Point::new(-2.0, 0.5, 1.0), Vector::new(0.2, -1.0, 0.7)),
        );
        add(
            "plane/sphere through the centre",
            plane(Point::ORIGIN, Vector::Z),
            sphere(Point::ORIGIN, 3.0),
        );
        add(
            "plane/sphere off centre",
            plane(Point::new(0.0, 0.0, 1.5), Vector::Z),
            sphere(Point::ORIGIN, 3.0),
        );
        add(
            "plane/sphere oblique",
            plane(Point::new(0.4, -0.2, 0.9), Vector::new(1.0, 2.0, 3.0)),
            sphere(Point::new(1.0, 1.0, 1.0), 4.0),
        );
        add(
            "plane/cylinder perpendicular",
            plane(Point::new(0.0, 0.0, 2.0), Vector::Z),
            cylinder(Point::ORIGIN, Vector::Z, 2.0),
        );
        add(
            "plane/cylinder oblique",
            plane(Point::ORIGIN, Vector::new(0.0, 1.0, 1.0)),
            cylinder(Point::ORIGIN, Vector::Z, 2.0),
        );
        add(
            "plane/cylinder along the axis",
            plane(Point::new(0.5, 0.0, 0.0), Vector::X),
            cylinder(Point::ORIGIN, Vector::Z, 2.0),
        );
        add(
            "plane/cylinder tangent",
            plane(Point::new(2.0, 0.0, 0.0), Vector::X),
            cylinder(Point::ORIGIN, Vector::Z, 2.0),
        );
        add(
            "sphere/sphere crossing",
            sphere(Point::ORIGIN, 3.0),
            sphere(Point::new(4.0, 0.0, 0.0), 2.0),
        );
        add(
            "sphere/sphere oblique",
            sphere(Point::new(1.0, -2.0, 0.5), 5.0),
            sphere(Point::new(-3.0, 1.0, 2.0), 4.0),
        );
        add(
            "cylinder/sphere coaxial",
            cylinder(Point::ORIGIN, Vector::Z, 1.5),
            sphere(Point::ORIGIN, 3.0),
        );
        add(
            "plane/cone perpendicular",
            plane(Point::new(0.0, 0.0, 2.0), Vector::Z),
            cone(Point::ORIGIN, Vector::Z, 3.0, 0.2),
        );
        add(
            "cylinder/cone coaxial",
            cylinder(Point::ORIGIN, Vector::Z, 2.0),
            cone(Point::ORIGIN, Vector::Z, 3.0, 0.3),
        );
        add(
            "cone/cone coaxial crossing",
            cone(Point::ORIGIN, Vector::Z, 3.0, 0.2),
            cone(Point::new(0.0, 0.0, 1.0), Vector::Z, 2.0, 0.4),
        );
        out
    }

    #[test]
    fn every_closed_form_lands_on_both_surfaces() {
        // The gate's own measurement, run as an assertion. Every point of every
        // reported curve is on both surfaces to machine precision — which is
        // the defining property of an intersection curve and the only one that
        // can be checked without a second implementation to compare against.
        let report = measure_all(&corpus(), T);
        println!(
            "intersection benchmark: {} cases, {} solved, {} deferred, worst \
             deviation {:e}{}",
            report.cases,
            report.solved,
            report.deferred,
            report.worst,
            report
                .worst_case
                .as_ref()
                .map_or(String::new(), |c| format!(" ({c})"))
        );

        assert_eq!(
            report.deferred, 0,
            "every case here should have a closed form"
        );
        assert_eq!(report.solved, report.cases);
        assert!(
            report.within(1e-12),
            "worst deviation {:e} at {:?}",
            report.worst,
            report.worst_case
        );
    }

    #[test]
    fn an_axis_normal_plane_meets_a_cone_in_the_circle_at_that_height() {
        // Radius 3 at the reference, slope tan(0.2) per unit: at height 2 the
        // parallel's radius is exactly the closed form's.
        let Meeting::Along(curves) = surface_surface(
            &plane(Point::new(0.0, 0.0, 2.0), Vector::Z),
            &cone(Point::ORIGIN, Vector::Z, 3.0, 0.2),
            T,
        )
        .unwrap() else {
            panic!("the perpendicular slice should be a curve");
        };
        assert_eq!(curves.len(), 1);
        let Curve::Circle(circle) = &curves[0] else {
            panic!("the parallel should be a circle, got {curves:?}");
        };
        let expected = 0.2_f64.tan().mul_add(2.0, 3.0);
        assert!((circle.circle().radius() - expected).abs() < 1e-12);

        // Through the apex it is a touch — a point, not a zero-length curve.
        let apex_height = -3.0 / 0.2_f64.tan();
        assert!(matches!(
            surface_surface(
                &plane(Point::new(0.0, 0.0, apex_height), Vector::Z),
                &cone(Point::ORIGIN, Vector::Z, 3.0, 0.2),
                T,
            )
            .unwrap(),
            Meeting::Touching(ref p) if p.len() == 1
        ));

        // Oblique stays deferred by name.
        assert!(
            surface_surface(
                &plane(Point::ORIGIN, Vector::new(0.0, 1.0, 1.0)),
                &cone(Point::ORIGIN, Vector::Z, 3.0, 0.2),
                T,
            )
            .is_err()
        );
    }

    #[test]
    fn a_coaxial_cylinder_meets_a_cone_in_the_parallel_on_its_own_nappe() {
        // The mirrored crossing past the apex is real geometry whose chart
        // parameters run half a turn out of phase; reporting it would poison
        // the arrangement of any face nearby, so only the chart's own nappe
        // answers.
        let Meeting::Along(curves) = surface_surface(
            &cylinder(Point::ORIGIN, Vector::Z, 2.0),
            &cone(Point::ORIGIN, Vector::Z, 3.0, 0.3),
            T,
        )
        .unwrap() else {
            panic!("a coaxial cylinder should cross the slant");
        };
        assert_eq!(curves.len(), 1, "the parallel on the chart's own nappe");
        for curve in &curves {
            let Curve::Circle(circle) = curve else {
                panic!("a parallel should be a circle");
            };
            assert!((circle.circle().radius() - 2.0).abs() < 1e-12);
        }
        // An off-axis cylinder is the marcher's business.
        assert!(
            surface_surface(
                &cylinder(Point::new(1.0, 0.0, 0.0), Vector::Z, 2.0),
                &cone(Point::ORIGIN, Vector::Z, 3.0, 0.3),
                T,
            )
            .is_err()
        );
    }

    #[test]
    fn coaxial_equal_cones_are_the_same_surface() {
        // The same cone described from a frame two units up its own axis: the
        // reference radius grows by the slope times the lift.
        let lifted = 0.2_f64.tan().mul_add(2.0, 3.0);
        assert_eq!(
            surface_surface(
                &cone(Point::ORIGIN, Vector::Z, 3.0, 0.2),
                &cone(Point::new(0.0, 0.0, 2.0), Vector::Z, lifted, 0.2),
                T,
            )
            .unwrap(),
            Meeting::Same
        );
        // Parallel slants that never meet are apart, not almost-the-same.
        assert_eq!(
            surface_surface(
                &cone(Point::ORIGIN, Vector::Z, 3.0, 0.2),
                &cone(Point::ORIGIN, Vector::Z, 4.0, 0.2),
                T,
            )
            .unwrap(),
            Meeting::Apart
        );
        // Crossing slants meet in the parallel where the radii agree.
        let Meeting::Along(curves) = surface_surface(
            &cone(Point::ORIGIN, Vector::Z, 3.0, 0.2),
            &cone(Point::new(0.0, 0.0, 1.0), Vector::Z, 2.0, 0.4),
            T,
        )
        .unwrap() else {
            panic!("crossing slants should meet along a parallel");
        };
        assert!(!curves.is_empty());
    }

    #[test]
    fn the_kind_of_answer_is_right_and_not_only_its_accuracy() {
        // Landing on both surfaces is necessary and not sufficient: an
        // intersector returning one circle of two would score perfectly. These
        // pin the *shape* of each answer.
        let case = |a: SurfaceGeometry, b: SurfaceGeometry| surface_surface(&a, &b, T).unwrap();

        assert!(matches!(
            case(plane(Point::ORIGIN, Vector::Z), plane(Point::ORIGIN, Vector::X)),
            Meeting::Along(ref c) if c.len() == 1
        ));
        assert_eq!(
            case(
                plane(Point::ORIGIN, Vector::Z),
                plane(Point::new(0.0, 0.0, 1.0), Vector::Z)
            ),
            Meeting::Apart
        );
        assert_eq!(
            case(
                plane(Point::ORIGIN, Vector::Z),
                plane(Point::new(3.0, 4.0, 0.0), Vector::Z)
            ),
            Meeting::Same
        );

        // A sphere resting on a plane touches; it does not meet along anything.
        assert!(matches!(
            case(
                plane(Point::ORIGIN, Vector::Z),
                sphere(Point::new(0.0, 0.0, 2.0), 2.0)
            ),
            Meeting::Touching(ref p) if p.len() == 1
        ));

        // A plane through a cylinder's axis cuts two lines, not one.
        assert!(matches!(
            case(
                plane(Point::ORIGIN, Vector::X),
                cylinder(Point::ORIGIN, Vector::Z, 2.0)
            ),
            Meeting::Along(ref c) if c.len() == 2
        ));

        // A sphere larger than a coaxial cylinder cuts it in two circles.
        assert!(matches!(
            case(
                cylinder(Point::ORIGIN, Vector::Z, 1.0),
                sphere(Point::ORIGIN, 3.0)
            ),
            Meeting::Along(ref c) if c.len() == 2
        ));
    }

    #[test]
    fn an_oblique_plane_cuts_a_cylinder_in_an_ellipse_of_the_right_size() {
        // The closed form's whole value: not a fitted curve that is nearly an
        // ellipse, but the ellipse, with the radii geometry says it has.
        let angle = core::f64::consts::FRAC_PI_3;
        let radius = 2.0;
        let cut = plane(Point::ORIGIN, Vector::new(0.0, angle.sin(), angle.cos()));
        let drum = cylinder(Point::ORIGIN, Vector::Z, radius);
        let Meeting::Along(curves) = surface_surface(&cut, &drum, T).unwrap() else {
            panic!("an oblique cut should meet along a curve");
        };
        assert_eq!(curves.len(), 1);
        let Curve::Ellipse(e) = &curves[0] else {
            panic!(
                "an oblique cut of a cylinder is an ellipse, got {:?}",
                curves[0]
            );
        };
        approx::assert_relative_eq!(e.ellipse().minor_radius(), radius, max_relative = 1e-12);
        approx::assert_relative_eq!(
            e.ellipse().major_radius(),
            radius / angle.cos(),
            max_relative = 1e-12
        );
    }

    #[test]
    fn a_pair_with_no_closed_form_is_deferred_rather_than_guessed() {
        // The honest half. Two cylinders on skew axes meet in a quartic space
        // curve; returning something plausible would be the single worst thing
        // this module could do, because the boolean above it would trust it.
        let a = cylinder(Point::ORIGIN, Vector::Z, 1.0);
        let b = cylinder(Point::ORIGIN, Vector::X, 1.0);
        let err = surface_surface(&a, &b, T).unwrap_err();
        assert!(
            err.to_string().contains("marching"),
            "unexpected message: {err}"
        );

        let deferred = measure_all(&[("cylinder/cylinder skew".to_string(), a, b)], T);
        assert_eq!(deferred.deferred, 1);
        assert_eq!(deferred.solved, 0);
        assert_eq!(deferred.worst, 0.0, "a deferred case scores nothing");
    }

    #[test]
    fn nothing_to_sample_is_reported_as_nothing_rather_than_as_zero_error() {
        // A deviation of zero and no measurement at all are different, and
        // averaging them together would flatter every report containing a pair
        // that simply misses.
        let found = measure(
            &plane(Point::ORIGIN, Vector::Z),
            &plane(Point::new(0.0, 0.0, 5.0), Vector::Z),
            T,
        )
        .unwrap();
        assert_eq!(found.meeting, Meeting::Apart);
        assert_eq!(found.deviation, None);
        assert_eq!(found.samples, 0);
    }
}

/// And the reported curves are all of the intersection.
mod completeness {
    use crate::support::coverage::*;
    use ogeom_core::Tolerances;
    use ogeom_geom::SurfaceGeometry;
    use ogeom_geom::{CylinderSurface, PlaneSurface, SphereSurface};
    use ogeom_intersect::{Marching, branches};
    use ogeom_math::Point;
    use ogeom_math::{Cylinder, Direction, Frame, Plane, Sphere, Vector};

    const T: Tolerances = Tolerances::millimetres();

    fn sphere(centre: Point, radius: f64) -> SurfaceGeometry {
        SphereSurface::new(Sphere::centred(centre, radius, T).unwrap()).into()
    }

    fn cylinder(axis: Vector, radius: f64) -> SurfaceGeometry {
        let frame = Frame::new(
            Point::ORIGIN,
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

    fn options() -> Marching {
        Marching {
            chord: 1e-4,
            ..Marching::default()
        }
    }

    #[test]
    fn dropping_a_branch_is_caught() {
        // The test that makes the instrument worth having. A completeness
        // measure that always says "complete" looks exactly like a correct one
        // until something is actually missing, so the negative control is not
        // optional — it is the only evidence the thing works.
        let a = sphere(Point::ORIGIN, 3.0);
        let b = cylinder(Vector::Z, 1.5);

        let found = branches(&a, &b, options(), T).unwrap();
        assert_eq!(found.len(), 2, "a coaxial cylinder cuts a sphere twice");

        let whole = coverage(&a, &b, &found, 40, T).unwrap();
        assert!(
            whole.complete(),
            "{} of {} crossings covered, first miss at {:?}",
            whole.covered,
            whole.crossings,
            whole.missed.first()
        );

        // Now hide one, exactly as an intersector that never seeded it would.
        let half = coverage(&a, &b, &found[..1], 40, T).unwrap();
        assert!(
            !half.complete(),
            "dropping a whole branch went unnoticed, which means this measures \
             nothing"
        );
        println!(
            "coverage with one branch of two: {:.1}% ({} missed)",
            half.fraction() * 100.0,
            half.missed.len()
        );
        assert!(half.fraction() < 0.75, "got {}", half.fraction());

        // And reporting nothing at all is caught most of all.
        let none = coverage(&a, &b, &[], 40, T).unwrap();
        assert_eq!(none.covered, 0);
        assert!(none.crossings > 0);
    }

    #[test]
    fn the_marching_intersector_finds_all_of_what_it_is_asked_for() {
        // The measurement the gate wants, on the cases that have no closed
        // form. Accuracy was already established; this is the other half.
        let cases: Vec<(&str, SurfaceGeometry, SurfaceGeometry)> = vec![
            (
                "sphere/plane",
                sphere(Point::ORIGIN, 3.0),
                plane(Point::new(0.0, 0.0, 1.0), Vector::Z),
            ),
            (
                "sphere/cylinder coaxial",
                sphere(Point::ORIGIN, 3.0),
                cylinder(Vector::Z, 1.5),
            ),
            (
                "sphere/cylinder offset",
                sphere(Point::new(0.6, 0.0, 0.0), 3.0),
                cylinder(Vector::Z, 1.5),
            ),
            (
                "crossed cylinders",
                cylinder(Vector::Z, 1.0),
                cylinder(Vector::X, 1.6),
            ),
        ];
        for (name, a, b) in cases {
            let found = branches(&a, &b, options(), T).unwrap();
            let score = coverage(&a, &b, &found, 40, T).unwrap();
            println!(
                "coverage {name}: {}/{} cells, {} branches",
                score.covered,
                score.crossings,
                found.len()
            );
            assert!(score.crossings > 0, "{name}: nothing to cover");
            assert!(
                score.complete(),
                "{name}: missed {} of {} crossings, first at {:?}",
                score.crossings - score.covered,
                score.crossings,
                score.missed.first()
            );
        }
    }

    #[test]
    fn a_pair_it_cannot_measure_says_so_rather_than_scoring_full_marks() {
        // A cone has no consistent inside, so there is no sign to change. The
        // dangerous answer would be "complete", since a caller reading a
        // hundred percent has no way to tell it apart from a real result.
        let cone: SurfaceGeometry = ogeom_geom::ConeSurface::new(
            ogeom_math::Cone::new(Frame::WORLD, 1.0, 0.4_f64.atan(), T).unwrap(),
            (0.0, 3.0),
        )
        .unwrap()
        .into();
        let err = coverage(&sphere(Point::ORIGIN, 2.0), &cone, &[], 20, T).unwrap_err();
        assert!(err.to_string().contains("signed distance"), "got {err}");
    }

    #[test]
    fn nothing_crossed_is_complete_rather_than_a_division_by_zero() {
        let far = sphere(Point::new(100.0, 0.0, 0.0), 1.0);
        let score = coverage(&sphere(Point::ORIGIN, 1.0), &far, &[], 16, T).unwrap();
        assert_eq!(score.crossings, 0);
        assert!(score.complete());
        assert!((score.fraction() - 1.0).abs() < f64::EPSILON);
    }

    #[test]
    fn a_grid_too_coarse_to_have_cells_is_refused() {
        let a = sphere(Point::ORIGIN, 1.0);
        let b = plane(Point::ORIGIN, Vector::Z);
        assert!(coverage(&a, &b, &[], 1, T).is_err());
        assert!(coverage(&a, &b, &[], 0, T).is_err());
    }
}
