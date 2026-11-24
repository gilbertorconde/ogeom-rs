//! The tight bound is where the shape reaches, on spheres and tori in any
//! frame as on boxes: each side at the surface's own extreme.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::Model;

const T: Tolerances = Tolerances::millimetres();

fn tilted(origin: Point, z: Vector, x: Vector) -> Frame {
    Frame::new(
        origin,
        Direction::new(z, T).unwrap(),
        Direction::new(x, T).unwrap(),
        T,
    )
    .unwrap()
}

/// A sphere's box is its centre plus or minus its radius, and a torus's
/// reaches its major radius scaled by how square each axis stands to the
/// torus's own, plus its minor radius, whatever the frame.
#[test]
fn spheres_and_tori_in_tilted_frames_bound_exactly() {
    let frames = [
        tilted(
            Point::new(300.0, -120.0, 40.0),
            Vector::new(0.3, -0.5, 1.0),
            Vector::new(1.0, 0.2, 0.0),
        ),
        tilted(
            Point::new(-6.8, 196.2, -190.2),
            Vector::new(-0.86, 0.51, 0.05),
            Vector::new(0.1, 0.2, 1.0),
        ),
    ];
    for frame in frames {
        for radius in [3.0, 30.0] {
            let mut model = Model::new();
            let ball = ogeom_algo::make_sphere(&mut model, frame, radius, T)
                .unwrap()
                .shape;
            let bounds = ogeom_algo::tight_bounds(&model, &ball, T).unwrap();
            let c = frame.origin();
            let reach = Vector::new(radius, radius, radius);
            assert!(
                bounds.low().unwrap().distance(c - reach) < 1e-7,
                "{bounds:?}"
            );
            assert!(
                bounds.high().unwrap().distance(c + reach) < 1e-7,
                "{bounds:?}"
            );

            let minor = radius / 4.0;
            let ring = ogeom_algo::make_torus(&mut model, frame, radius, minor, T)
                .unwrap()
                .shape;
            let bounds = ogeom_algo::tight_bounds(&model, &ring, T).unwrap();
            let z = frame.z().vector();
            let side = |k: f64| radius * (1.0 - k * k).sqrt() + minor;
            let reach = Vector::new(side(z.x), side(z.y), side(z.z));
            assert!(
                bounds.low().unwrap().distance(c - reach) < 1e-7,
                "{bounds:?}"
            );
            assert!(
                bounds.high().unwrap().distance(c + reach) < 1e-7,
                "{bounds:?}"
            );
        }
    }
}

/// A drum and a cone reach no further than their rims: each rim, a circle
/// on a tilted axis `d`, reaches its centre plus or minus its radius times
/// `sqrt(1 - d_i^2)` along each axis `i`.
#[test]
fn drums_and_cones_in_tilted_frames_bound_at_their_rims() {
    let axis = Vector::new(0.3, -0.5, 1.0);
    let d = axis / axis.magnitude();
    let frame = tilted(
        Point::new(300.0, -120.0, 40.0),
        axis,
        Vector::new(1.0, 0.2, 0.0),
    );
    let (radius, height) = (5.0, 12.0);
    let rim = |c: Point, r: f64, low: &mut Point, high: &mut Point| {
        let reach = Vector::new(
            r * (1.0 - d.x * d.x).sqrt(),
            r * (1.0 - d.y * d.y).sqrt(),
            r * (1.0 - d.z * d.z).sqrt(),
        );
        let (a, b) = (c - reach, c + reach);
        *low = Point::new(low.x.min(a.x), low.y.min(a.y), low.z.min(a.z));
        *high = Point::new(high.x.max(b.x), high.y.max(b.y), high.z.max(b.z));
    };
    let base = frame.origin();
    let top = base + d * height;

    let mut model = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, frame, radius, height, T)
        .unwrap()
        .shape;
    let bounds = ogeom_algo::tight_bounds(&model, &drum, T).unwrap();
    let (mut low, mut high) = (
        Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
        Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    rim(base, radius, &mut low, &mut high);
    rim(top, radius, &mut low, &mut high);
    assert!(bounds.low().unwrap().distance(low) < 1e-7, "{bounds:?}");
    assert!(bounds.high().unwrap().distance(high) < 1e-7, "{bounds:?}");

    let cone = ogeom_algo::make_cone(&mut model, frame, radius, radius / 2.0, height, T)
        .unwrap()
        .shape;
    let bounds = ogeom_algo::tight_bounds(&model, &cone, T).unwrap();
    let (mut low, mut high) = (
        Point::new(f64::INFINITY, f64::INFINITY, f64::INFINITY),
        Point::new(f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    rim(base, radius, &mut low, &mut high);
    rim(top, radius / 2.0, &mut low, &mut high);
    assert!(bounds.low().unwrap().distance(low) < 1e-7, "{bounds:?}");
    assert!(bounds.high().unwrap().distance(high) < 1e-7, "{bounds:?}");
}
