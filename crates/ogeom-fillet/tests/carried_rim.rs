//! A boss's foot round transfers from a plate to the plate's straight edge.
#![allow(clippy::unwrap_used, reason = "test constructions")]

use ogeom_core::Tolerances;
use ogeom_geom::{CircleCurve, Curve, LineCurve, PlaneSurface};
use ogeom_math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn envelope(theta: f64, fraction: f64) -> Point {
    let (sin, cos) = theta.sin_cos();
    let over = (1.6 * sin - 1.0).max(0.0);
    let height = (1.0 - over * over).sqrt();
    let centre = Point::new(1.6 * cos, 1.6 * sin, height);
    let a = Vector::new(-cos, -sin, 0.0);
    let b = Vector::new(0.0, -over, -height);
    let dot = a.dot(b);
    let towards = (b - a * dot) / (1.0 - dot * dot).sqrt();
    let phi = dot.acos() * fraction;
    centre + a * phi.cos() + towards * phi.sin()
}

// Integrate z*n_z over the independent ball envelope. All other added
// faces are vertical or at z=0, so their divergence-theorem flux vanishes.
fn added_volume(along: u32, across: u32) -> f64 {
    let switch = (1.0_f64 / 1.6).asin();
    let mut cap = 0.0;
    let triangle = |a: Point, b: Point, c: Point| (b - a).cross(c - a).z * (a.z + b.z + c.z) / 6.0;
    for (lo, hi) in [
        (0.0, switch),
        (switch, core::f64::consts::PI - switch),
        (core::f64::consts::PI - switch, core::f64::consts::PI),
    ] {
        let row = |i| {
            let theta = lo + (hi - lo) * f64::from(i) / f64::from(along);
            (0..=across)
                .map(|j| envelope(theta, f64::from(j) / f64::from(across)))
                .collect::<Vec<_>>()
        };
        let mut a = row(0);
        for i in 1..=along {
            let b = row(i);
            for j in 0..a.len() - 1 {
                cap -= triangle(a[j], b[j], b[j + 1]) + triangle(a[j], b[j + 1], a[j + 1]);
            }
            a = b;
        }
    }
    20.0 * (1.0 - core::f64::consts::FRAC_PI_4) + 2.0 * cap
}

#[test]
#[ignore = "heavy"]
fn a_stadium_foot_round_rolls_on_the_plate_edge() {
    let mut model = Model::new();
    let points = [
        Point::new(0.6, -5.0, 0.0),
        Point::new(0.6, 5.0, 0.0),
        Point::new(-0.6, 5.0, 0.0),
        Point::new(-0.6, -5.0, 0.0),
    ];
    let vertices: Vec<_> = points
        .iter()
        .map(|&p| ogeom_algo::make_vertex(&mut model, p).shape)
        .collect();
    let mut profile = Vec::new();
    for i in 0..4 {
        let (curve, range) = if i % 2 == 0 {
            (
                Curve::Line(LineCurve::segment(points[i], points[(i + 1) % 4], T).unwrap()),
                (0.0, 10.0),
            )
        } else {
            let y = if i == 1 { 5.0 } else { -5.0 };
            let frame = Frame::new(Point::new(0.0, y, 0.0), Direction::Z, Direction::X, T).unwrap();
            let lo = if i == 1 { 0.0 } else { core::f64::consts::PI };
            (
                Curve::Circle(CircleCurve::new(Circle::new(frame, 0.6, T).unwrap())),
                (lo, lo + core::f64::consts::PI),
            )
        };
        profile.push(
            ogeom_algo::make_edge_between(
                &mut model,
                curve,
                range,
                &vertices[i],
                &vertices[(i + 1) % 4],
                T,
            )
            .unwrap()
            .shape,
        );
    }
    let wire = ogeom_algo::make_wire(&mut model, &profile, T)
        .unwrap()
        .shape;
    let face = ogeom_algo::make_face(
        &mut model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let boss = ogeom_algo::make_prism(&mut model, &face, Vector::new(0.0, 0.0, 10.0), T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(-10.0, -6.0, -2.0), Direction::Z, Direction::X, T).unwrap();
    let plate = ogeom_algo::make_box(&mut model, frame, (20.0, 12.0, 2.0), T)
        .unwrap()
        .shape;
    let solid = ogeom_bool::fuse(&mut model, &plate, &boss, T)
        .unwrap()
        .shape;
    let edges: Vec<_> = explore_unique(&model, &solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            ogeom_algo::edge_vertices(&model, e)
                .unwrap()
                .is_some_and(|(a, b)| {
                    [&a, &b].iter().all(|v| {
                        let p = v
                            .transform(model.datums())
                            .unwrap()
                            .apply(model.node(v).unwrap().data().as_vertex().unwrap().point);
                        p.z.abs() < T.confusion() && p.x.abs() < 1.0 && p.y.abs() < 5.7
                    })
                })
        })
        .collect();
    assert_eq!(edges.len(), 4);
    let built = ogeom_fillet::fillet_edges(&mut model, &solid, &edges, 1.0, T).unwrap();
    assert!(
        ogeom_algo::check(&model, &built.shape, T)
            .unwrap()
            .is_valid()
    );
    let coarse = added_volume(256, 128);
    let fine = added_volume(512, 256);
    // Second-order triangle quadrature must converge before extrapolation.
    assert!((fine - coarse).abs() < 6e-5);
    let expected = fine + (fine - coarse) / 3.0;
    let mass =
        ogeom_algo::volume_properties(&model, &built.shape, ogeom_mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    let original = 480.0 + 10.0 * (12.0 + core::f64::consts::PI * 0.36);
    // The independently integrated envelope is held within 1e-5 mm³.
    assert!(
        (mass - original - expected).abs() < 1e-5,
        "{} vs {expected}",
        mass - original
    );
    let switch = (1.0_f64 / 1.6).asin();
    for theta in [
        0.0,
        switch - 1e-6,
        switch,
        switch + 1e-6,
        core::f64::consts::FRAC_PI_2,
        core::f64::consts::PI - switch,
        core::f64::consts::PI,
    ] {
        for j in 0..=16 {
            let p = envelope(theta, f64::from(j) / 16.0);
            for p in [
                Point::new(p.x, p.y + 5.0, p.z),
                Point::new(p.x, -p.y - 5.0, p.z),
            ] {
                let vertex = ogeom_algo::make_vertex(&mut model, p).shape;
                let distance = ogeom_algo::distance_between_shapes(
                    &model,
                    &vertex,
                    &built.shape,
                    ogeom_intersect::ExtremaOptions::default(),
                    T,
                )
                .unwrap()
                .distance;
                // A tight fit is checked against exact ball contacts and arcs.
                assert!(distance < 1e-5, "{p:?}: {distance}");
            }
        }
    }
}
