//! Fitted skins, fillings and projections measured between their samples:
//! each result sampled densely and compared with the geometry it stands
//! for, independently of the points it was fitted to.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Circle, Direction, Frame, Plane, Point};
use ogeom_topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// The spline surfaces of every face of `shape`.
fn spline_surfaces(model: &Model, shape: &Shape) -> Vec<ogeom_geom::BSplineSurface> {
    let mut out = Vec::new();
    for face in explore(model, shape, Filter::OfType(ShapeType::Face)).unwrap() {
        let NodeData::Face(data) = model.node(&face).unwrap().data() else {
            continue;
        };
        if let Some(SurfaceGeometry::BSpline(patch)) = model.geometry().surface(data.surface) {
            out.push(patch.clone());
        }
    }
    out
}

/// The parameter `k / n` of the way across `range`.
fn across(range: (f64, f64), k: usize, n: usize) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    let f = k as f64 / n as f64;
    range.0 + (range.1 - range.0) * f
}

/// The distance from `p` to the polyline `line`.
fn to_polyline(p: Point, line: &[Point]) -> f64 {
    let mut best = f64::INFINITY;
    for w in line.windows(2) {
        let d = w[1] - w[0];
        let len2 = d.dot(d);
        let t = if len2 > 0.0 {
            ((p - w[0]).dot(d) / len2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        best = best.min(p.distance(w[0] + d * t));
    }
    best
}

fn edge_of(model: &mut Model, curve: ogeom_geom::Curve) -> Shape {
    let domain = curve.domain();
    ogeom_algo::make_edge(model, curve, domain, T)
        .unwrap()
        .shape
}

/// A cubic spline spine riding `periods` periods of a sine of `amplitude`
/// along `length` of x: its control points on the sine, twenty a period,
/// on uniform knots. Also returns the spine densely sampled.
fn sine_spine(model: &mut Model, amplitude: f64, length: f64, periods: f64) -> (Shape, Vec<Point>) {
    let wave = |x: f64| Point::new(x, amplitude * (2.0 * PI * periods * x / length).sin(), 0.0);
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let n = (20.0 * periods) as usize;
    let control: Vec<Point> = (0..=n).map(|k| wave(across((0.0, length), k, n))).collect();
    let mut knots = vec![0.0; 4];
    knots.extend((1..n - 2).map(|k| across((0.0, 1.0), k, n - 2)));
    knots.extend([1.0; 4]);
    let curve = ogeom_geom::Curve::BSpline(
        ogeom_geom::BSplineCurve::new(ogeom_math::KnotVector::new(knots, 3).unwrap(), control, T)
            .unwrap(),
    );
    let (lo, hi) = curve.domain();
    let dense: Vec<Point> = (0..=20_000)
        .map(|k| curve.point_at(across((lo, hi), k, 20_000), T).unwrap())
        .collect();
    (edge_of(model, curve), dense)
}

#[test]
fn a_skinned_pipe_along_a_wavy_spine_holds_its_radius_between_stations() {
    let mut model = Model::new();
    let (spine, dense) = sine_spine(&mut model, 1.0, 40.0, 8.0);
    let (r, tolerance) = (0.3, 1e-3);
    let pipe = ogeom_offset::make_pipe_skinned(&mut model, &spine, r, tolerance, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom_algo::check(&model, &pipe, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let walls = spline_surfaces(&model, &pipe);
    assert_eq!(walls.len(), 1);
    let wall = &walls[0];
    let (ud, vd) = wall.domain();
    let mut worst = 0.0_f64;
    for j in 0..=400 {
        for i in 0..=40 {
            let p = wall
                .point_at(across(ud, i, 40), across(vd, j, 400), T)
                .unwrap();
            // The spine runs along x: only its stretch within reach counts.
            let lo = dense.partition_point(|q| q.x < p.x - 1.0);
            let hi = dense
                .partition_point(|q| q.x <= p.x + 1.0)
                .min(dense.len() - 1);
            worst = worst.max((to_polyline(p, &dense[lo.saturating_sub(1)..=hi]) - r).abs());
        }
    }
    eprintln!("pipe wall off its radius by {worst}");
    assert!(worst <= tolerance, "the wall strays {worst} from the tube");
}

/// The distance from `p` to the circle of `radius` about the origin in the
/// XY plane.
fn to_flat_circle(p: Point, radius: f64) -> f64 {
    (p.x.hypot(p.y) - radius).hypot(p.z)
}

#[test]
fn a_filling_holds_its_curved_border_between_samples() {
    // A half disc's arc below the x axis, closed by three lines through a
    // height of one.
    let mut model = Model::new();
    let arc = ogeom_geom::CircleCurve::new(Circle::new(Frame::WORLD, 1.0, T).unwrap());
    let corners = [
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 1.0),
        Point::new(-1.0, 0.0, 1.0),
        Point::new(-1.0, 0.0, 0.0),
    ];
    let bottom = ogeom_algo::make_edge(&mut model, arc.into(), (PI, 2.0 * PI), T)
        .unwrap()
        .shape;
    let line = |model: &mut Model, a: Point, b: Point| {
        ogeom_algo::make_edge(
            model,
            ogeom_geom::LineCurve::segment(a, b, T).unwrap().into(),
            (0.0, a.distance(b)),
            T,
        )
        .unwrap()
        .shape
    };
    let right = line(&mut model, corners[0], corners[1]);
    let top = line(&mut model, corners[1], corners[2]);
    let left = line(&mut model, corners[2], corners[3]);
    let tolerance = 1e-4;
    let filled =
        ogeom_offset::make_filling(&mut model, &[bottom, right, top, left], 4, tolerance, T)
            .unwrap()
            .shape;
    let patches = spline_surfaces(&model, &filled);
    let patch = &patches[0];
    let (ud, vd) = patch.domain();
    let mut worst = 0.0_f64;
    for i in 0..=2000 {
        let p = patch.point_at(across(ud, i, 2000), vd.0, T).unwrap();
        worst = worst.max(to_flat_circle(p, 1.0));
    }
    eprintln!("filling border off its arc by {worst}");
    assert!(worst <= tolerance, "the border strays {worst} from its arc");
}

#[test]
fn a_skinned_loft_through_thin_ellipses_holds_its_end_sections() {
    let mut model = Model::new();
    let (a, b) = (10.0, 0.2);
    let ellipse = |z: f64| {
        let frame = Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X, T).unwrap();
        ogeom_geom::EllipseCurve::new(ogeom_math::Ellipse::new(frame, a, b, T).unwrap())
    };
    let sections: Vec<Shape> = [0.0, 5.0, 10.0]
        .iter()
        .map(|&z| {
            let edge = edge_of(&mut model, ellipse(z).into());
            ogeom_algo::make_wire(&mut model, &[edge], T).unwrap().shape
        })
        .collect();
    let tolerance = 1e-4;
    let loft = ogeom_offset::make_loft_skinned(&mut model, &sections, tolerance, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom_algo::check(&model, &loft, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    // The ellipse at its height, densely enough that the polyline's sag is
    // far below the tolerance.
    let ring = |z: f64| -> Vec<Point> {
        let curve = ellipse(z);
        (0..=40_000)
            .map(|k| {
                curve
                    .point_at(across((0.0, 2.0 * PI), k, 40_000), T)
                    .unwrap()
            })
            .collect()
    };
    let rings = [ring(0.0), ring(10.0)];
    let walls = spline_surfaces(&model, &loft);
    let wall = &walls[0];
    let (ud, vd) = wall.domain();
    let mut worst = 0.0_f64;
    for v in [vd.0, vd.1] {
        for i in 0..=1000 {
            let p = wall.point_at(across(ud, i, 1000), v, T).unwrap();
            let line = if p.z < 5.0 { &rings[0] } else { &rings[1] };
            worst = worst.max(to_polyline(p, line));
        }
    }
    eprintln!("loft end rings off their ellipses by {worst}");
    assert!(
        worst <= tolerance,
        "an end ring strays {worst} from its ellipse"
    );
}

/// A gently curved spline spine from the origin, setting off along +z.
fn bent_spine(model: &mut Model) -> Shape {
    let curve = ogeom_geom::Curve::BSpline(
        ogeom_geom::BSplineCurve::new(
            ogeom_math::KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3).unwrap(),
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(0.0, 0.0, 3.0),
                Point::new(1.0, 0.0, 6.0),
                Point::new(3.0, 0.0, 9.0),
            ],
            T,
        )
        .unwrap(),
    );
    let edge = edge_of(model, curve);
    ogeom_algo::make_wire(model, &[edge], T).unwrap().shape
}

#[test]
fn a_pipe_shell_holds_its_profile_arc_between_samples() {
    // A half disc in the XY plane: the arc of radius one over +y, closed
    // by its diameter.
    let mut model = Model::new();
    let (a, b) = (Point::new(-1.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0));
    let (va, vb) = (
        ogeom_algo::make_vertex(&mut model, a).shape,
        ogeom_algo::make_vertex(&mut model, b).shape,
    );
    let arc = ogeom_algo::make_edge_between(
        &mut model,
        ogeom_geom::CircleCurve::new(Circle::new(Frame::WORLD, 1.0, T).unwrap()).into(),
        (0.0, PI),
        &vb,
        &va,
        T,
    )
    .unwrap()
    .shape;
    let diameter = ogeom_algo::make_edge_between(
        &mut model,
        ogeom_geom::LineCurve::segment(a, b, T).unwrap().into(),
        (0.0, 2.0),
        &va,
        &vb,
        T,
    )
    .unwrap()
    .shape;
    let wire = ogeom_algo::make_wire(&mut model, &[arc, diameter], T)
        .unwrap()
        .shape;
    let plane: SurfaceGeometry = ogeom_geom::PlaneSurface::new(Plane::new(Frame::WORLD)).into();
    let profile = ogeom_algo::make_face(&mut model, plane, &[wire], T)
        .unwrap()
        .shape;
    let spine = bent_spine(&mut model);
    let tolerance = 1e-5;
    let pipe = ogeom_offset::make_pipe_shell(&mut model, &profile, &spine, false, tolerance, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom_algo::check(&model, &pipe, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    // Every spline border lying in the profile's plane off its diameter is
    // the arc.
    let mut worst = 0.0_f64;
    let mut seen = 0;
    for patch in spline_surfaces(&model, &pipe) {
        let (ud, vd) = patch.domain();
        for k in 0..=1000 {
            for (u, v) in [
                (across(ud, k, 1000), vd.0),
                (across(ud, k, 1000), vd.1),
                (ud.0, across(vd, k, 1000)),
                (ud.1, across(vd, k, 1000)),
            ] {
                let p = patch.point_at(u, v, T).unwrap();
                if p.z.abs() < 1e-3 && p.y > 1e-2 {
                    seen += 1;
                    worst = worst.max(to_flat_circle(p, 1.0));
                }
            }
        }
    }
    assert!(seen > 500, "the arc's strip has a border on the profile");
    eprintln!("pipe shell border off its arc by {worst}");
    assert!(worst <= tolerance, "the border strays {worst} from its arc");
}

#[test]
fn a_projected_circle_holds_its_stated_tolerance_between_stations() {
    // A circle of radius 3 at height 10 over a ball of radius 5: its foot
    // on the ball is the circle of radius 15 / sqrt(109) at height
    // 50 / sqrt(109).
    let mut model = Model::new();
    let ball = ogeom_algo::make_sphere(&mut model, Frame::WORLD, 5.0, T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(0.0, 0.0, 10.0), Direction::Z, Direction::X, T).unwrap();
    let circle = edge_of(
        &mut model,
        ogeom_geom::CircleCurve::new(Circle::new(frame, 3.0, T).unwrap()).into(),
    );
    let wire = ogeom_algo::make_wire(&mut model, &[circle], T)
        .unwrap()
        .shape;
    let tolerance = 1e-3;
    let (landed, _) =
        ogeom_offset::normal_projection(&mut model, &ball, &wire, 8, tolerance, T).unwrap();
    assert!(!landed.is_empty());
    let (radius, height) = (15.0 / 109.0_f64.sqrt(), 50.0 / 109.0_f64.sqrt());
    for stretch in &landed {
        let data = model.node(&stretch.edge).unwrap().data().as_edge().unwrap();
        let EdgeRepr::Curve3d { curve, range, .. } = data.curve3d().unwrap() else {
            unreachable!()
        };
        let geometry = model.geometry().curve(*curve).unwrap();
        let mut worst = 0.0_f64;
        for k in 0..=2000 {
            let p = geometry.point_at(across(*range, k, 2000), T).unwrap();
            worst = worst.max((p.x.hypot(p.y) - radius).hypot(p.z - height));
        }
        eprintln!(
            "projection off its foot by {worst}, stating {}",
            stretch.tolerance
        );
        assert!(
            worst <= stretch.tolerance.max(1e-9) && stretch.tolerance <= tolerance,
            "the projection strays {worst} and states {}",
            stretch.tolerance
        );
    }
    let diagnosis = ogeom_algo::check(&model, &landed[0].edge, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
}

#[test]
fn a_pipe_through_thin_sections_holds_its_end_section() {
    // An ellipse 10 by 0.2 square to a quarter arc of radius 20 at each of
    // its ends, the arc turning in the ellipse's long direction.
    let mut model = Model::new();
    let ellipse = |model: &mut Model, frame: Frame| {
        let curve =
            ogeom_geom::EllipseCurve::new(ogeom_math::Ellipse::new(frame, 10.0, 0.2, T).unwrap());
        let edge = edge_of(model, curve.into());
        ogeom_algo::make_wire(model, &[edge], T).unwrap().shape
    };
    let start = ellipse(
        &mut model,
        Frame::new(Point::ORIGIN, -Direction::Y, Direction::X, T).unwrap(),
    );
    let end = ellipse(
        &mut model,
        Frame::new(Point::new(20.0, -20.0, 0.0), Direction::X, -Direction::Y, T).unwrap(),
    );
    let arc = Circle::new(
        Frame::new(Point::new(20.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap(),
        20.0,
        T,
    )
    .unwrap();
    let spine = ogeom_algo::make_edge(
        &mut model,
        ogeom_geom::CircleCurve::new(arc).into(),
        (PI, 1.5 * PI),
        T,
    )
    .unwrap()
    .shape;
    let spine = ogeom_algo::make_wire(&mut model, &[spine], T)
        .unwrap()
        .shape;
    let tolerance = 1e-3;
    let pipe =
        ogeom_offset::make_pipe_sections(&mut model, &[start, end], &spine, false, tolerance, T)
            .unwrap()
            .shape;
    let diagnosis = ogeom_algo::check(&model, &pipe, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let ring: Vec<Point> = (0..=40_000)
        .map(|k| {
            let a = 2.0 * PI * f64::from(k) / 40_000.0;
            Point::new(10.0 * a.cos(), 0.0, 0.2 * a.sin())
        })
        .collect();
    let mut worst = 0.0_f64;
    let mut seen = 0;
    for patch in spline_surfaces(&model, &pipe) {
        let (ud, vd) = patch.domain();
        for k in 0..=2000 {
            for (u, v) in [
                (across(ud, k, 2000), vd.0),
                (across(ud, k, 2000), vd.1),
                (ud.0, across(vd, k, 2000)),
                (ud.1, across(vd, k, 2000)),
            ] {
                let p = patch.point_at(u, v, T).unwrap();
                if p.y.abs() < 1e-6 {
                    seen += 1;
                    worst = worst.max(to_polyline(p, &ring));
                }
            }
        }
    }
    assert!(seen > 1000, "the pipe has a border on the start section");
    eprintln!("pipe sections' start ring off its ellipse by {worst}");
    assert!(worst <= tolerance, "the start ring strays {worst}");
}
