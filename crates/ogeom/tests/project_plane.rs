//! An edge projected orthogonally onto a plane comes back as the exact
//! curve it projects to, in the plane's own coordinates: what a sketch
//! needs to constrain against a solid's edges.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{ProjectedCurve, make_box, make_cylinder, make_edge, project_edge_onto_plane};
use ogeom::core::Tolerances;
use ogeom::geom::{BSplineCurve, Curve, Curve2d as _, Curve3d, HelixCurve};
use ogeom::math::{Direction, Frame, KnotVector, Plane, Point, Point2, Vector};
use ogeom::topo::{EdgeRepr, Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// The edge of `shape` whose middle lies at `at`.
fn edge_at(model: &Model, shape: &Shape, at: Point) -> Shape {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|edge| {
            let data = model.node(edge).unwrap().data().as_edge().unwrap();
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                return false;
            };
            let curve = model.geometry().curve(*curve).unwrap();
            let middle = curve.point_at(f64::midpoint(range.0, range.1), T).unwrap();
            middle.distance(at) < 1e-9
        })
        .expect("an edge there")
}

fn close(a: Point2, b: (f64, f64)) -> bool {
    a.distance(Point2::new(b.0, b.1)) < 1e-9
}

/// A plane through `origin`, its normal `normal` and its `x` axis `x`.
fn plane(origin: Point, normal: Vector, x: Vector) -> Plane {
    Plane::new(
        Frame::new(
            origin,
            Direction::new(normal, T).unwrap(),
            Direction::new(x, T).unwrap(),
            T,
        )
        .unwrap(),
    )
}

#[test]
fn a_box_top_edge_projects_onto_its_own_plane_as_itself() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let edge = edge_at(&model, &block, Point::new(5.0, 0.0, 10.0));
    let top = plane(Point::new(0.0, 0.0, 10.0), Vector::Z, Vector::X);
    match project_edge_onto_plane(&model, &edge, &top, T).unwrap() {
        ProjectedCurve::Line { start, end } => {
            let ends = [start, end];
            assert!(ends.iter().any(|p| close(*p, (0.0, 0.0))), "{ends:?}");
            assert!(ends.iter().any(|p| close(*p, (10.0, 0.0))), "{ends:?}");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn a_vertical_edge_projects_to_a_point() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let edge = edge_at(&model, &block, Point::new(10.0, 10.0, 5.0));
    match project_edge_onto_plane(&model, &edge, &Plane::XY, T).unwrap() {
        ProjectedCurve::Point(p) => assert!(close(p, (10.0, 10.0)), "{p:?}"),
        other => panic!("{other:?}"),
    }
}

/// A cylinder's top rim: a full circle of radius 5 on a parallel plane, an
/// ellipse of ratio one half on a plane tilted sixty degrees, and a
/// segment the rim's width on a plane standing square to it.
#[test]
fn a_cylinder_rim_projects_to_a_circle_an_ellipse_or_a_segment() {
    let mut model = Model::new();
    let drum = make_cylinder(&mut model, Frame::WORLD, 5.0, 8.0, T)
        .unwrap()
        .shape;
    let rim = edge_at(&model, &drum, Point::new(-5.0, 0.0, 8.0));

    match project_edge_onto_plane(&model, &rim, &Plane::XY, T).unwrap() {
        ProjectedCurve::Circle {
            centre,
            radius,
            range,
        } => {
            assert!(close(centre, (0.0, 0.0)), "{centre:?}");
            assert!((radius - 5.0).abs() < 1e-12);
            assert!(((range.1 - range.0) - 2.0 * PI).abs() < 1e-12, "{range:?}");
        }
        other => panic!("{other:?}"),
    }

    let (s, c) = (PI / 3.0).sin_cos();
    let tilted = plane(Point::ORIGIN, Vector::new(0.0, -s, c), Vector::X);
    match project_edge_onto_plane(&model, &rim, &tilted, T).unwrap() {
        ProjectedCurve::Ellipse {
            major,
            ratio,
            range,
            ..
        } => {
            assert!((major.magnitude() - 5.0).abs() < 1e-12, "{major:?}");
            assert!(major.y.abs() < 1e-12, "the major axis lies along the hinge");
            assert!((ratio - 0.5).abs() < 1e-12, "{ratio}");
            assert!(((range.1 - range.0) - 2.0 * PI).abs() < 1e-12);
        }
        other => panic!("{other:?}"),
    }

    let square = plane(Point::ORIGIN, -Vector::Y, Vector::X);
    match project_edge_onto_plane(&model, &rim, &square, T).unwrap() {
        ProjectedCurve::Line { start, end } => {
            assert!(
                close(start, (-5.0, 8.0)) && close(end, (5.0, 8.0)),
                "{start:?} {end:?}"
            );
        }
        other => panic!("{other:?}"),
    }
}

/// A half circle projected onto its own plane keeps its half: the arc's
/// angles, counter-clockwise, sweep pi from where it starts.
#[test]
fn an_arc_keeps_its_range() {
    let mut model = Model::new();
    let circle = ogeom::math::Circle::new(Frame::WORLD, 3.0, T).unwrap();
    let curve = Curve::Circle(ogeom::geom::CircleCurve::new(circle));
    let arc = make_edge(&mut model, curve, (0.0, PI), T).unwrap().shape;
    // Seen from below, the arc turns the other way: still a half turn,
    // now over the lower half of the plane's own angles.
    let below = plane(Point::ORIGIN, -Vector::Z, Vector::X);
    match project_edge_onto_plane(&model, &arc, &below, T).unwrap() {
        ProjectedCurve::Circle { range, .. } => {
            assert!(((range.1 - range.0) - PI).abs() < 1e-12, "{range:?}");
            let middle = f64::midpoint(range.0, range.1);
            assert!((middle.sin() + 1.0).abs() < 1e-12, "{range:?}");
        }
        other => panic!("{other:?}"),
    }
}

/// A cubic edge's projection is the cubic on its projected control points,
/// the same knots and nothing fitted.
#[test]
fn a_bspline_edge_projects_exactly() {
    let mut model = Model::new();
    let control = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, 3.0),
        Point::new(3.0, -1.0, 5.0),
        Point::new(4.0, 1.0, 2.0),
        Point::new(6.0, 0.0, 1.0),
    ];
    let knots = KnotVector::clamped_uniform(3, control.len()).unwrap();
    let spline = BSplineCurve::new(knots.clone(), control.clone(), T).unwrap();
    let range = spline.domain();
    let edge = make_edge(&mut model, Curve::BSpline(spline), range, T)
        .unwrap()
        .shape;
    let tilted = plane(
        Point::new(1.0, 1.0, 1.0),
        Vector::new(1.0, 1.0, 2.0),
        Vector::new(1.0, -1.0, 0.0),
    );
    let frame = tilted.frame();
    match project_edge_onto_plane(&model, &edge, &tilted, T).unwrap() {
        ProjectedCurve::BSpline { curve, fit_error } => {
            assert!(fit_error.is_none());
            assert_eq!(curve.knots(), &knots);
            for (got, p) in curve.control_points().iter().zip(&control) {
                let d = *p - frame.origin();
                let want = Point2::new(d.dot(frame.x().vector()), d.dot(frame.y().vector()));
                assert!(got.scaled.distance(want) < 1e-12, "{got:?} vs {want:?}");
            }
        }
        other => panic!("{other:?}"),
    }
}

/// A helix has no closed-form projection: it comes back fitted, the error
/// reported, and the fit lying on the circle the helix projects to.
#[test]
fn a_helix_projects_to_a_fitted_spline_with_its_error() {
    let mut model = Model::new();
    let helix = HelixCurve::new(Frame::WORLD, 4.0, 2.0, 1.0).unwrap();
    let range = helix.domain();
    let edge = make_edge(&mut model, Curve::Helix(helix), range, T)
        .unwrap()
        .shape;
    match project_edge_onto_plane(&model, &edge, &Plane::XY, T).unwrap() {
        ProjectedCurve::BSpline { curve, fit_error } => {
            let error = fit_error.expect("a fitted curve says how far off");
            assert!(error < 1e-5, "{error}");
            let (a, b) = curve.domain();
            for k in 0..=50 {
                let t = a + (b - a) * f64::from(k) / 50.0;
                let p = curve.point_at(t, T).unwrap();
                assert!((p.to_vector().magnitude() - 4.0).abs() < 1e-4, "{p:?}");
            }
        }
        other => panic!("{other:?}"),
    }
}

/// An edge covering the middle half of a spline projects only that half:
/// the projection starts and ends where the edge does.
#[test]
fn a_bspline_edge_on_part_of_its_curve_projects_that_part() {
    let mut model = Model::new();
    let control = vec![
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 2.0, 3.0),
        Point::new(3.0, -1.0, 5.0),
        Point::new(4.0, 1.0, 2.0),
        Point::new(6.0, 0.0, 1.0),
    ];
    let knots = KnotVector::clamped_uniform(3, control.len()).unwrap();
    let spline = BSplineCurve::new(knots, control, T).unwrap();
    let (a, b) = spline.domain();
    let range = (a + (b - a) * 0.25, a + (b - a) * 0.75);
    let ends = [
        spline.point_at(range.0, T).unwrap(),
        spline.point_at(range.1, T).unwrap(),
    ];
    let edge = make_edge(&mut model, Curve::BSpline(spline), range, T)
        .unwrap()
        .shape;
    match project_edge_onto_plane(&model, &edge, &Plane::XY, T).unwrap() {
        ProjectedCurve::BSpline { curve, .. } => {
            let (lo, hi) = curve.domain();
            let (first, last) = (
                curve.point_at(lo, T).unwrap(),
                curve.point_at(hi, T).unwrap(),
            );
            assert!(close(first, (ends[0].x, ends[0].y)), "{first:?}");
            assert!(close(last, (ends[1].x, ends[1].y)), "{last:?}");
        }
        other => panic!("{other:?}"),
    }
}
