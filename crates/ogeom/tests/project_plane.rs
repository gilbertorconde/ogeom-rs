//! An edge projected orthogonally onto a plane comes back as the exact
//! curve it projects to, in the plane's own coordinates: what a sketch
//! needs to constrain against a solid's edges.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{ProjectedCurve, make_box, make_cylinder, make_edge, project_edge_onto_plane};
use ogeom::core::Tolerances;
use ogeom::geom::{BSplineCurve, Curve, Curve2d as _, Curve3d, HelixCurve};
use ogeom::math::{Direction, Frame, KnotVector, Plane, Point, Point2, Vector, Vector2};
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

/// A curve drawn on a plane with more wiggle than the fit's samples can
/// follow projects onto that plane as itself everywhere, not only at the
/// samples, at the edge's own parameter, and the stated error covers what
/// a dense measure finds.
#[test]
fn a_wiggling_curve_projects_onto_its_plane_between_the_samples() {
    use ogeom::geom::{BSpline2d, CurveOnSurface, PlanarCurve, PlaneSurface, SurfaceGeometry};
    const SPANS: usize = 400;
    let control: Vec<Point2> = (0..SPANS + 3)
        .map(|i| {
            #[allow(clippy::cast_precision_loss)]
            let x = 10.0 * i as f64 / (SPANS + 2) as f64;
            Point2::new(x, if i.is_multiple_of(2) { 0.03 } else { -0.03 })
        })
        .collect();
    let chart = BSpline2d::new(
        KnotVector::clamped_uniform(3, control.len()).unwrap(),
        control,
        T,
    )
    .unwrap();
    let plane: SurfaceGeometry = PlaneSurface::new(Plane::XY).into();
    let drawn = Curve::OnSurface(Box::new(CurveOnSurface::new(
        PlanarCurve::BSpline(chart.clone()),
        plane,
    )));
    let mut model = Model::new();
    let edge = make_edge(&mut model, drawn, (0.0, 1.0), T).unwrap().shape;
    let ProjectedCurve::BSpline { curve, fit_error } =
        project_edge_onto_plane(&model, &edge, &Plane::XY, T).unwrap()
    else {
        panic!("a curve on a surface projects to a fitted spline");
    };
    let error = fit_error.expect("a fitted curve says how far off");
    // The projection runs at the edge's own parameter: at every parameter
    // it stands where the drawn curve does.
    let (a, b) = curve.domain();
    assert!(a.abs() < 1e-12 && (b - 1.0).abs() < 1e-12, "{a} {b}");
    let mut off = 0.0_f64;
    for k in 0..=20_000 {
        let t = f64::from(k) / 20_000.0;
        off = off.max(
            curve
                .point_at(t, T)
                .unwrap()
                .distance(chart.point_at(t, T).unwrap()),
        );
    }
    assert!(off <= 1e-4, "the projection strays {off} from the curve");
    assert!(
        off <= error + 1e-6,
        "the projection strays {off} and states {error}"
    );
}

/// A parabola and a hyperbola in a plane tilted thirty degrees about `x`,
/// projected onto `xy`: each comes back as the conic it projects to, every
/// projected point of the edge on it to the confusion distance.
#[test]
fn a_tilted_parabola_and_hyperbola_project_to_conics() {
    use ogeom::geom::{HyperbolaCurve, ParabolaCurve};
    use ogeom::math::{Hyperbola, Parabola};
    let tilt = Direction::new(Vector::new(0.0, -0.5, 0.75_f64.sqrt()), T).unwrap();
    let frame = Frame::new(Point::new(0.0, 0.0, 5.0), tilt, Direction::X, T).unwrap();
    let ground = plane(Point::ORIGIN, Vector::Z, Vector::X);
    let flat = |p: Point| Point2::new(p.x, p.y);
    let mut model = Model::new();

    let parabola =
        Curve::Parabola(ParabolaCurve::new(Parabola::new(frame, 2.0, T).unwrap(), 100.0).unwrap());
    let edge = make_edge(&mut model, parabola.clone(), (-3.0, 3.0), T)
        .unwrap()
        .shape;
    let ProjectedCurve::Parabola {
        vertex,
        axis,
        focal,
        range,
    } = project_edge_onto_plane(&model, &edge, &ground, T).unwrap()
    else {
        panic!("a parabola projects to a parabola");
    };
    let side = Vector2::new(-axis.y, axis.x);
    for k in 0..=50 {
        let t = -3.0 + 6.0 * f64::from(k) / 50.0;
        let d = flat(parabola.point_at(t, T).unwrap()) - vertex;
        let s = d.dot(side);
        let on = vertex + axis * (s * s / (4.0 * focal)) + side * s;
        assert!(
            on.distance(vertex + d) <= T.confusion(),
            "t {t}: off by {}",
            on.distance(vertex + d)
        );
        assert!(
            s >= range.0 - 1e-9 && s <= range.1 + 1e-9,
            "t {t}: s {s} outside {range:?}"
        );
    }

    let hyperbola = Curve::Hyperbola(
        HyperbolaCurve::new(Hyperbola::new(frame, 4.0, 2.0, T).unwrap(), 10.0).unwrap(),
    );
    let edge = make_edge(&mut model, hyperbola.clone(), (-1.0, 1.0), T)
        .unwrap()
        .shape;
    let ProjectedCurve::Hyperbola {
        centre,
        major,
        minor,
        range,
    } = project_edge_onto_plane(&model, &edge, &ground, T).unwrap()
    else {
        panic!("a hyperbola projects to a hyperbola");
    };
    assert!(major.dot(minor).abs() <= 1e-12 * major.magnitude() * minor.magnitude());
    for k in 0..=50 {
        let t = -1.0 + 2.0 * f64::from(k) / 50.0;
        let p = flat(hyperbola.point_at(t, T).unwrap());
        let d = p - centre;
        let s = (d.dot(minor) / minor.dot(minor)).asinh();
        let on = centre + major * s.cosh() + minor * s.sinh();
        assert!(
            on.distance(p) <= T.confusion(),
            "t {t}: off by {}",
            on.distance(p)
        );
        assert!(
            s >= range.0 - 1e-9 && s <= range.1 + 1e-9,
            "t {t}: s {s} outside {range:?}"
        );
    }
}

/// A parabola or a hyperbola whose plane stands square to the target is
/// seen edge-on: a segment, its ends the farthest points the piece reaches
/// along the line, the parabola folded back on its axis included.
#[test]
fn a_conic_square_to_the_plane_projects_to_a_segment() {
    use ogeom::geom::{HyperbolaCurve, ParabolaCurve};
    use ogeom::math::{Hyperbola, Parabola};
    let ground = plane(Point::ORIGIN, Vector::Z, Vector::X);
    let mut model = Model::new();
    let segment =
        |model: &Model, edge: &Shape| match project_edge_onto_plane(model, edge, &ground, T)
            .unwrap()
        {
            ProjectedCurve::Line { start, end } => {
                let (a, b) = if start.x <= end.x {
                    (start, end)
                } else {
                    (end, start)
                };
                (a, b)
            }
            other => panic!("a segment, not {other:?}"),
        };
    // Opening along z: the parabola's side, along x, is all that shows.
    let up = Frame::new(Point::ORIGIN, Direction::Y, Direction::Z, T).unwrap();
    let p = Curve::Parabola(ParabolaCurve::new(Parabola::new(up, 1.0, T).unwrap(), 100.0).unwrap());
    let edge = make_edge(&mut model, p, (-2.0, 3.0), T).unwrap().shape;
    let (a, b) = segment(&model, &edge);
    assert!(close(a, (-2.0, 0.0)) && close(b, (3.0, 0.0)), "{a:?} {b:?}");
    // Opening along x, its side along z: folded back on the axis, from the
    // vertex to the farther end, t^2 / 4 at t = 3.
    let flat = Frame::new(Point::ORIGIN, Direction::Y, Direction::X, T).unwrap();
    let p =
        Curve::Parabola(ParabolaCurve::new(Parabola::new(flat, 1.0, T).unwrap(), 100.0).unwrap());
    let edge = make_edge(&mut model, p, (-2.0, 3.0), T).unwrap().shape;
    let (a, b) = segment(&model, &edge);
    assert!(close(a, (0.0, 0.0)) && close(b, (2.25, 0.0)), "{a:?} {b:?}");
    // A hyperbola with its minor axis along z: `4 cosh t` along x.
    let h = Curve::Hyperbola(
        HyperbolaCurve::new(Hyperbola::new(flat, 4.0, 2.0, T).unwrap(), 10.0).unwrap(),
    );
    let edge = make_edge(&mut model, h, (-1.0, 2.0), T).unwrap().shape;
    let (a, b) = segment(&model, &edge);
    assert!(
        close(a, (4.0, 0.0)) && close(b, (4.0 * 2.0_f64.cosh(), 0.0)),
        "{a:?} {b:?}"
    );
}

/// A periodic cubic spline (uniform knots, not clamped) projects to a
/// clamped spline, exact, through the edge's projected points.
#[test]
fn a_periodic_spline_projects_clamped() {
    let control = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(4.0, 0.0, 1.0),
        Point::new(5.0, 3.0, 2.0),
        Point::new(2.0, 5.0, 1.0),
        Point::new(-1.0, 3.0, 0.0),
    ];
    let spline = BSplineCurve::periodic(&control, 3, T).unwrap();
    assert!(!spline.knots().is_clamped());
    let curve = Curve::BSpline(spline);
    let range = curve.domain();
    let mut model = Model::new();
    let edge = make_edge(&mut model, curve.clone(), range, T)
        .unwrap()
        .shape;
    let ground = plane(Point::ORIGIN, Vector::Z, Vector::X);
    let ProjectedCurve::BSpline {
        curve: projected,
        fit_error,
    } = project_edge_onto_plane(&model, &edge, &ground, T).unwrap()
    else {
        panic!("a spline projects to a spline");
    };
    assert_eq!(fit_error, None);
    assert!(projected.knots().is_clamped());
    let (a, b) = projected.domain();
    for k in 0..=50 {
        let f = f64::from(k) / 50.0;
        let p = curve
            .point_at(range.0 + (range.1 - range.0) * f, T)
            .unwrap();
        let q = projected.point_at(a + (b - a) * f, T).unwrap();
        assert!(q.distance(Point2::new(p.x, p.y)) <= T.confusion(), "at {f}");
    }
}

/// The conic edges a box's planar face cuts from a cone, each with the
/// plane its face lies on.
fn conic_edges(model: &Model, part: &Shape) -> Vec<(Shape, Plane)> {
    use ogeom::geom::SurfaceGeometry;
    let mut out = Vec::new();
    for face in explore_unique(model, part, ShapeType::Face).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        let SurfaceGeometry::Plane(plane) = model.geometry().surface(data.surface).unwrap() else {
            continue;
        };
        for edge in explore_unique(model, &face, ShapeType::Edge).unwrap() {
            let edge_data = model.node(&edge).unwrap().data().as_edge().unwrap();
            let Some(EdgeRepr::Curve3d { curve, .. }) = edge_data.curve3d() else {
                continue;
            };
            let kind = model.geometry().curve(*curve).unwrap().kind();
            if !matches!(
                kind,
                ogeom::geom::CurveKind::Line | ogeom::geom::CurveKind::Circle
            ) {
                out.push((edge, plane.plane()));
            }
        }
    }
    out
}

/// A cone cut by a plane parallel to its axis keeps the hyperbola the
/// plane and the cone meet in, and by a plane parallel to a ruling the
/// parabola: the edge's own curve is the conic, and it projects onto its
/// plane as one, every point on it.
#[test]
fn a_cone_cut_by_a_plane_keeps_its_exact_conic() {
    let mut model = Model::new();
    let cone = ogeom::algo::make_cone(&mut model, Frame::WORLD, 10.0, 0.0, 10.0, T)
        .unwrap()
        .shape;
    // x >= 3: the hyperbola (10 - z)^2 - y^2 = 9 on the plane x = 3.
    let across = Frame::new(Point::new(3.0, -20.0, -5.0), Direction::Z, Direction::X, T).unwrap();
    let block = make_box(&mut model, across, (20.0, 40.0, 20.0), T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(&mut model, &cone, &block, T)
        .unwrap()
        .shape;
    let edges = conic_edges(&model, &cut);
    assert!(!edges.is_empty(), "the cut face has a conic edge");
    for (edge, plane) in &edges {
        let projected = project_edge_onto_plane(&model, edge, plane, T).unwrap();
        assert!(
            matches!(projected, ProjectedCurve::Hyperbola { .. }),
            "{projected:?}"
        );
        let data = model.node(edge).unwrap().data().as_edge().unwrap();
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            panic!("the edge has a curve");
        };
        let curve = model.geometry().curve(*curve).unwrap();
        for k in 0..=16 {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 16.0;
            let p = curve.point_at(t, T).unwrap();
            let off = (10.0 - p.z).powi(2) - p.y * p.y - 9.0;
            assert!(
                (p.x - 3.0).abs() < 1e-9 && off.abs() < 1e-6,
                "{p:?} off by {off}"
            );
        }
    }

    // Square to the ruling the cone has at u = 0: a parabola.
    let mut model = Model::new();
    let cone = ogeom::algo::make_cone(&mut model, Frame::WORLD, 10.0, 0.0, 10.0, T)
        .unwrap()
        .shape;
    let normal = Direction::new(Vector::new(1.0, 0.0, 1.0), T).unwrap();
    let along = Direction::new(Vector::new(1.0, 0.0, -1.0), T).unwrap();
    let side = normal.vector().cross(along.vector());
    let corner = Point::new(4.0, 0.0, 4.0) - along.vector() * 30.0 - side * 30.0;
    let tilted = Frame::new(corner, normal, along, T).unwrap();
    let block = make_box(&mut model, tilted, (60.0, 60.0, 30.0), T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(&mut model, &cone, &block, T)
        .unwrap()
        .shape;
    let edges = conic_edges(&model, &cut);
    assert!(!edges.is_empty(), "the cut face has a conic edge");
    for (edge, plane) in &edges {
        let projected = project_edge_onto_plane(&model, edge, plane, T).unwrap();
        assert!(
            matches!(projected, ProjectedCurve::Parabola { .. }),
            "{projected:?}"
        );
    }
}

/// Whether every point of `curve` over `range`, projected onto `onto`,
/// lies on `projected` within 1e-9.
fn lies_on(curve: &Curve, range: (f64, f64), onto: &Plane, projected: &ProjectedCurve) -> bool {
    let frame = onto.frame();
    (0..=64).all(|k| {
        let t = range.0 + (range.1 - range.0) * f64::from(k) / 64.0;
        let d = curve.point_at(t, T).unwrap() - frame.origin();
        let p = Point2::new(d.dot(frame.x().vector()), d.dot(frame.y().vector()));
        match projected {
            ProjectedCurve::Line { start, end } => {
                let (along, across) = (*end - *start, p - *start);
                let f = across.dot(along) / along.dot(along);
                (-1e-12..=1.0 + 1e-12).contains(&f) && (*start + along * f).distance(p) < 1e-9
            }
            ProjectedCurve::Circle {
                centre,
                radius,
                range,
            } => {
                let turn = (p.y - centre.y).atan2(p.x - centre.x);
                let turn = range.0 + (turn - range.0).rem_euclid(core::f64::consts::TAU);
                (p.distance(*centre) - radius).abs() < 1e-9 && turn <= range.1 + 1e-9
            }
            ProjectedCurve::Ellipse {
                centre,
                major,
                ratio,
                ..
            } => {
                let minor = Vector2::new(-major.y, major.x) * *ratio;
                let r = p - *centre;
                let (a, b) = (
                    r.dot(*major) / major.dot(*major),
                    r.dot(minor) / minor.dot(minor),
                );
                (a * a + b * b - 1.0).abs() < 1e-9
            }
            _ => false,
        }
    })
}

/// A wall's rim left by a thickening: a circle of radius 10 offset by 2
/// in its own plane. Onto a parallel plane it projects as the arc of
/// radius 12 outward, 8 inward, over the same angles; onto a tilted plane
/// as an ellipse. Every point of the edge lies on what it projects to.
#[test]
fn an_offset_circle_projects_to_a_circle_or_an_ellipse() {
    let circle = ogeom::math::Circle::new(Frame::WORLD, 10.0, T).unwrap();
    let range = (0.3, 2.0);
    let tilted = plane(Point::ORIGIN, Vector::new(0.0, 0.6, 0.8), Vector::X);
    for (reference, radius) in [(Direction::Z, 12.0), (-Direction::Z, 8.0)] {
        let offset = ogeom::geom::OffsetCurve::new(
            Curve::Circle(ogeom::geom::CircleCurve::new(circle)),
            2.0,
            reference,
        )
        .unwrap();
        let curve = Curve::Offset(Box::new(offset));
        let mut model = Model::new();
        let edge = make_edge(&mut model, curve.clone(), range, T)
            .unwrap()
            .shape;
        let flat = project_edge_onto_plane(&model, &edge, &Plane::XY, T).unwrap();
        match &flat {
            ProjectedCurve::Circle {
                centre,
                radius: r,
                range: angles,
            } => {
                assert!(close(*centre, (0.0, 0.0)), "{centre:?}");
                assert!((r - radius).abs() < 1e-12, "{r} against {radius}");
                assert!((angles.0 - range.0).abs() < 1e-12, "{angles:?}");
                assert!((angles.1 - range.1).abs() < 1e-12, "{angles:?}");
            }
            other => panic!("{other:?}"),
        }
        assert!(lies_on(&curve, range, &Plane::XY, &flat));
        let leaning = project_edge_onto_plane(&model, &edge, &tilted, T).unwrap();
        assert!(
            matches!(leaning, ProjectedCurve::Ellipse { .. }),
            "{leaning:?}"
        );
        assert!(lies_on(&curve, range, &tilted, &leaning));
    }
}

/// A line offset by 2 is the line moved 2 square to it, and projects as
/// a segment.
#[test]
fn an_offset_line_projects_to_a_segment() {
    let line = ogeom::geom::LineCurve::new(ogeom::math::Axis {
        location: Point::new(1.0, 2.0, 3.0),
        direction: Direction::new(Vector::new(1.0, 1.0, 0.0), T).unwrap(),
    });
    let offset = ogeom::geom::OffsetCurve::new(Curve::Line(line), 2.0, Direction::Z).unwrap();
    let curve = Curve::Offset(Box::new(offset));
    let mut model = Model::new();
    let range = (-4.0, 6.0);
    let edge = make_edge(&mut model, curve.clone(), range, T)
        .unwrap()
        .shape;
    for onto in [
        Plane::XY,
        plane(Point::ORIGIN, Vector::new(0.0, 0.6, 0.8), Vector::X),
    ] {
        let projected = project_edge_onto_plane(&model, &edge, &onto, T).unwrap();
        assert!(
            matches!(projected, ProjectedCurve::Line { .. }),
            "{projected:?}"
        );
        assert!(lies_on(&curve, range, &onto, &projected));
    }
}
