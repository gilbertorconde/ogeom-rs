//! Global minima of the classic multimodal test functions, where a local
//! descent from the middle of the box lands in the wrong basin.
#![allow(clippy::unwrap_used, clippy::cast_precision_loss, reason = "test code")]

use ogeom_math::{global_minimum, minimize_local, swarm_minimum};

fn rastrigin(x: &[f64]) -> f64 {
    10.0 * x.len() as f64
        + x.iter()
            .map(|v| v * v - 10.0 * (core::f64::consts::TAU * v).cos())
            .sum::<f64>()
}

fn six_hump_camel(x: &[f64]) -> f64 {
    let (a, b) = (x[0], x[1]);
    (4.0 - 2.1 * a * a + a.powi(4) / 3.0) * a * a + a * b + (-4.0 + 4.0 * b * b) * b * b
}

fn wavy(x: &[f64]) -> f64 {
    x[0].sin() + (10.0 * x[0] / 3.0).sin()
}

#[test]
fn a_wavy_line_is_minimised_in_its_deepest_trough() {
    let found = global_minimum(wavy, &[2.7], &[7.5], 1e-6, 100_000).unwrap();
    assert!(found.certified);
    assert!((found.value - -1.899_599).abs() < 1e-5, "{found:?}");
    assert!((found.point[0] - 5.145_735).abs() < 1e-3, "{found:?}");
}

#[test]
fn the_camel_is_found_in_one_of_its_two_lowest_humps() {
    let found = global_minimum(six_hump_camel, &[-3.0, -2.0], &[3.0, 2.0], 1e-6, 200_000).unwrap();
    assert!(found.certified, "{found:?}");
    assert!((found.value - -1.031_628_45).abs() < 1e-5, "{found:?}");
    // A descent from the middle stalls on the saddle at the origin.
    let local = minimize_local(
        six_hump_camel,
        &[0.0, 0.0],
        &[-3.0, -2.0],
        &[3.0, 2.0],
        1e-3,
        1e-12,
        2000,
    )
    .unwrap();
    assert!(local.value >= found.value - 1e-9);
}

#[test]
fn rastrigin_is_found_at_the_origin_among_its_many_traps() {
    let lower = [-5.12, -5.12];
    let upper = [5.12, 5.12];
    let found = global_minimum(rastrigin, &lower, &upper, 1e-4, 500_000).unwrap();
    assert!(found.certified && found.value < 1e-4, "{found:?}");
    let local = minimize_local(rastrigin, &[3.1, -2.2], &lower, &upper, 0.1, 1e-12, 4000).unwrap();
    assert!(
        local.value > 1.0,
        "a descent from off the origin stays in its own trap"
    );
}

#[test]
fn the_swarm_finds_the_camel_too() {
    let found = swarm_minimum(six_hump_camel, &[-3.0, -2.0], &[3.0, 2.0], 30, 200).unwrap();
    assert!((found.value - -1.031_628_45).abs() < 1e-6, "{found:?}");
}

#[test]
fn a_malformed_box_is_refused() {
    assert!(global_minimum(wavy, &[1.0], &[0.0], 1e-6, 100).is_err());
    assert!(global_minimum(wavy, &[0.0], &[1.0], 0.0, 100).is_err());
    assert!(swarm_minimum(wavy, &[0.0, 0.0], &[1.0], 4, 4).is_err());
}
