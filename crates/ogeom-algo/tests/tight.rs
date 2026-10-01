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
