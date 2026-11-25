//! A later round reaching an earlier band at its end (ADR 0001, amendment
//! "a later round reaching an earlier band at its end"): rounds asked one
//! call after another, the later band running on past its end and trimmed
//! by the earlier band's surface.
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::Curve3d as _;
use ogeom_math::{Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The edge of `s` passing nearest `q`.
fn edge_near(m: &Model, s: &Shape, q: Point) -> Shape {
    let d = |e: &Shape| {
        let (c, (a, b)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap().unwrap();
        (0..=64)
            .map(|i| {
                c.point_at(a + (b - a) * f64::from(i) / 64.0, T)
                    .unwrap()
                    .distance(q)
            })
            .fold(f64::INFINITY, f64::min)
    };
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .min_by(|a, b| d(a).total_cmp(&d(b)))
        .unwrap()
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// A 10 x 8 x 6 box's two top edges at the corner (0, 0, 6), the front
/// along x and the left along y, rounded at `r` one call after the other,
/// in either order, the upright edge there left sharp. The later band ends
/// at the earlier, which is tangent to the top they share; it runs on and
/// the earlier band trims it, so the solid loses exactly the union of the
/// two bands' full-length removals: each takes `r^2 (1 - pi / 4)` per
/// unit of its edge, and at height `s` above `6 - r` both take the square
/// of side `w = r - sqrt(r^2 - s^2)`, whose integral over `s` is
/// `r^3 (5 / 3 - pi / 2)`. That is the mitre one call with both edges
/// builds, which the two orders match to rounding.
#[test]
fn a_later_round_runs_on_into_the_earlier_band_at_a_box_corner() {
    use core::f64::consts::{FRAC_PI_2, FRAC_PI_4};
    let r = 1.0_f64;
    let exact = 10.0 * 8.0 * 6.0 - r * r * (1.0 - FRAC_PI_4) * (10.0 + 8.0)
        + r * r * r * (5.0 / 3.0 - FRAC_PI_2);
    let (front, left) = (Point::new(5.0, 0.0, 6.0), Point::new(0.0, 4.0, 6.0));
    let mut volumes = Vec::new();
    for order in [[front, left], [left, front]] {
        let mut m = Model::new();
        let block = ogeom_algo::make_box(&mut m, Frame::WORLD, (10.0, 8.0, 6.0), T)
            .unwrap()
            .shape;
        let mut solid = block.clone();
        for at in order {
            let edge = edge_near(&m, &solid, at);
            solid = ogeom_fillet::fillet_edge(&mut m, &solid, &edge, r, T)
                .unwrap()
                .shape;
        }
        let diagnosis = ogeom_algo::check(&m, &solid, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        let got = volume(&m, &solid);
        assert!(
            (got - exact).abs() < 1e-9 * exact,
            "{order:?}: {got} against {exact}"
        );
        let both = [edge_near(&m, &block, front), edge_near(&m, &block, left)];
        let together = ogeom_fillet::fillet_edges(&mut m, &block, &both, r, T)
            .unwrap()
            .shape;
        let one_call = volume(&m, &together);
        assert!(
            (got - one_call).abs() < 1e-9 * one_call,
            "{order:?}: {got} one at a time against {one_call} in one call"
        );
        volumes.push(got);
    }
    assert!((volumes[0] - volumes[1]).abs() < 1e-9 * volumes[0]);
}
