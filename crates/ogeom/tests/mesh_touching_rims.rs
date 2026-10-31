//! A cylinder face between a circle and a spline edge that rises to the
//! circle and runs along it to the vertex they share, within the spline's
//! tolerance and a hair past the circle in the chart: the face's ring
//! crosses itself there at any chord. The stretch run out and back bounds
//! nothing. The solid meshes closed, and each face keeps its area.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    make_edge_between, make_face_with_pcurves, make_shell, make_solid, make_vertex,
    surface_properties,
};
use ogeom::core::{Tolerance, Tolerances};
use ogeom::geom::fit::fit_points;
use ogeom::geom::{CircleCurve, CylinderSurface, LineCurve, PlaneSurface, Surface as _};
use ogeom::math::{Circle, Cylinder, Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::{Deflection, triangulate_with_face_meshes};
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

const RADIUS: f64 = 10.0;
/// The sector's turn about the axis.
const TURN: f64 = 1.5;
const HEIGHT: f64 = 2.0;
/// How far the spline runs past the circle, along the axis, where it
/// runs along it.
const PAST: f64 = 5e-4;
/// The spline's tolerance, as a file states it: wider than how far it
/// runs past the circle.
const SPLINE_TOLERANCE: f64 = 1e-3;

/// The spline's height at angle `u`: up from the bottom at 0.9 to the
/// circle at height 1, tangent to it at 0.6, a hair past it out to 0.45,
/// and along it to its end at 0.2.
fn spline_height(u: f64) -> f64 {
    if u >= 0.6 {
        1.0 - ((u - 0.6) / 0.3).powi(2)
    } else if u >= 0.45 {
        1.0 + PAST * (core::f64::consts::PI * (0.6 - u) / 0.15).sin().powi(2)
    } else {
        1.0
    }
}

/// A sector of a cylinder of radius 10, a turn of 1.5 about its axis and
/// 2 high, whose curved wall is three faces: one above a circle at height
/// 1, one below a spline rising to that circle and running on along it to
/// its end, and the one between them, whose ring runs out along the
/// spline and back along the circle. Returns the solid and the three wall
/// faces.
fn sector(model: &mut Model) -> (Shape, [Shape; 3]) {
    let cylinder =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, RADIUS, T).unwrap(), (-1.0, 3.0)).unwrap();
    let on = |u: f64, v: f64| cylinder.point_at(u, v, T).unwrap();
    let mut vertex = |p: Point| make_vertex(model, p).shape;
    let [a0, a2] = [0.0, HEIGHT].map(|z| vertex(Point::new(0.0, 0.0, z)));
    let [b0, bc, bu, p1, tu, t0, l1, tip] = [
        (0.0, 0.0),
        (0.9, 0.0),
        (TURN, 0.0),
        (TURN, 1.0),
        (TURN, HEIGHT),
        (0.0, HEIGHT),
        (0.0, 1.0),
        (0.2, 1.0),
    ]
    .map(|(u, v)| vertex(on(u, v)));
    let at = |model: &Model, v: &Shape| model.node(v).unwrap().data().as_vertex().unwrap().point;
    let line = |model: &mut Model, from: &Shape, to: &Shape| {
        let (a, b) = (at(model, from), at(model, to));
        let curve = LineCurve::segment(a, b, T).unwrap();
        make_edge_between(model, curve.into(), (0.0, a.distance(b)), from, to, T)
            .unwrap()
            .shape
    };
    let arc = |model: &mut Model, z: f64, range: (f64, f64), from: &Shape, to: &Shape| {
        let frame = Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X, T).unwrap();
        let circle = CircleCurve::new(Circle::new(frame, RADIUS, T).unwrap());
        make_edge_between(model, circle.into(), range, from, to, T)
            .unwrap()
            .shape
    };
    let axis = line(model, &a0, &a2);
    let (low0, low_u) = (line(model, &a0, &b0), line(model, &a0, &bu));
    let (high0, high_u) = (line(model, &a2, &t0), line(model, &a2, &tu));
    let (left_low, left_high) = (line(model, &b0, &l1), line(model, &l1, &t0));
    let (right_low, right_high) = (line(model, &bu, &p1), line(model, &p1, &tu));
    let bottom_left = arc(model, 0.0, (0.0, 0.9), &b0, &bc);
    let bottom_right = arc(model, 0.0, (0.9, TURN), &bc, &bu);
    let top = arc(model, HEIGHT, (0.0, TURN), &t0, &tu);
    let rim_left = arc(model, 1.0, (0.0, 0.2), &l1, &tip);
    let rim = arc(model, 1.0, (0.2, TURN), &tip, &p1);
    let samples: Vec<Point> = (0..=280)
        .map(|k| {
            let u = 0.9 - 0.7 * f64::from(k) / 280.0;
            on(u, spline_height(u))
        })
        .collect();
    let fitted = fit_points(&samples, 3, 1e-7, T).unwrap();
    let range = fitted.curve.knots().domain();
    let spline = make_edge_between(model, fitted.curve.into(), range, &bc, &tip, T)
        .unwrap()
        .shape;
    model
        .widen(&spline, Tolerance::new(SPLINE_TOLERANCE).unwrap())
        .unwrap();

    let wall = || -> ogeom::geom::SurfaceGeometry { cylinder.into() };
    let mut curved = |ring: Vec<Shape>| {
        make_face_with_pcurves(model, wall(), &[ring], T)
            .unwrap()
            .shape
    };
    let between = curved(vec![
        bottom_right.clone(),
        right_low.clone(),
        rim.reversed(),
        spline.reversed(),
    ]);
    let above = curved(vec![
        rim_left.clone(),
        rim,
        right_high.clone(),
        top.reversed(),
        left_high.reversed(),
    ]);
    let below = curved(vec![
        bottom_left.clone(),
        spline,
        rim_left.reversed(),
        left_low.reversed(),
    ]);
    let plane = |origin: Point, normal: Vector, x: Vector| -> ogeom::geom::SurfaceGeometry {
        let (z, x) = (
            Direction::new(normal, T).unwrap(),
            Direction::new(x, T).unwrap(),
        );
        PlaneSurface::new(Plane::new(Frame::new(origin, z, x, T).unwrap())).into()
    };
    let mut flat = |surface, ring: Vec<Shape>| {
        make_face_with_pcurves(model, surface, &[ring], T)
            .unwrap()
            .shape
    };
    let floor = flat(
        plane(Point::ORIGIN, -Vector::Z, Vector::X),
        vec![
            low_u.clone(),
            bottom_right.reversed(),
            bottom_left.reversed(),
            low0.reversed(),
        ],
    );
    let roof = flat(
        plane(Point::new(0.0, 0.0, HEIGHT), Vector::Z, Vector::X),
        vec![high0.clone(), top, high_u.reversed()],
    );
    let start = flat(
        plane(Point::ORIGIN, -Vector::Y, Vector::X),
        vec![low0, left_low, left_high, high0.reversed(), axis.reversed()],
    );
    let (sin, cos) = TURN.sin_cos();
    let end = flat(
        plane(
            Point::ORIGIN,
            Vector::new(-sin, cos, 0.0),
            Vector::new(cos, sin, 0.0),
        ),
        vec![
            axis,
            high_u,
            right_high.reversed(),
            right_low.reversed(),
            low_u.reversed(),
        ],
    );
    let walls = [between, above, below];
    let mut faces = walls.to_vec();
    faces.extend([floor, roof, start, end]);
    let shell = make_shell(model, &faces).unwrap().shape;
    let solid = make_solid(model, &[shell]).unwrap().shape;
    (solid, walls)
}

/// The face between the circle and the spline has a ring that crosses
/// itself where the spline runs along the circle; meshed, the solid
/// closes against the faces either side, holds the sector's volume, and
/// each wall face keeps its exact area, the three together the wall's.
#[test]
fn a_spline_running_along_a_circle_to_their_vertex_meshes_closed() {
    let mut model = Model::new();
    let (solid, walls) = sector(&mut model);
    for chord in [0.1, 0.01] {
        let deflection = Deflection::with_chord(chord).unwrap();
        let (mesh, _, faces) =
            triangulate_with_face_meshes(&model, &solid, deflection, None, T).unwrap();
        assert!(mesh.is_closed(), "the mesh at chord {chord} is open");
        // A chord `c` on a circle of radius `r` cuts a polygon short of
        // the circle's area by about `c / 3r` of it, and short of a
        // sector's volume by twice that.
        let volume = TURN / 2.0 * RADIUS * RADIUS * HEIGHT;
        let short = chord / (3.0 * RADIUS);
        assert!(
            (mesh.volume() - volume).abs() <= 2.0 * short * volume,
            "volume {} against {volume} at chord {chord}",
            mesh.volume()
        );
        let mut all = 0.0;
        for wall in &walls {
            let exact = surface_properties(&model, wall, deflection, T).unwrap();
            assert_eq!(exact.deflection, 0.0, "the wall's area is not exact");
            let (_, meshed) = faces.iter().find(|(face, _)| face.is_same(wall)).unwrap();
            let area = meshed.as_ref().unwrap().area();
            assert!(
                (area - exact.mass).abs() <= short * exact.mass,
                "meshed area {area} against exact {} at chord {chord}",
                exact.mass
            );
            all += exact.mass;
        }
        let wall = TURN * RADIUS * HEIGHT;
        assert!(
            (all - wall).abs() <= 1e-6 * wall,
            "exact areas sum to {all} against {wall}"
        );
    }
}
