//! The public rotation-minimizing frames are the frames a pipe sweep
//! places its sections in: a rectangle off the spine swept down a helix
//! under the rotation-minimizing law leaves each corner's rail where the
//! frames, read at the same share of the length, carry that corner.
#![allow(clippy::unwrap_used, clippy::cast_precision_loss, reason = "test code")]

use ogeom_algo::{
    make_edge, make_polygon, placed_edge_curve, project_on_curve, rotation_minimizing_frames,
};
use ogeom_core::Tolerances;
use ogeom_geom::{Curve3d as _, HelixCurve};
use ogeom_math::{Frame, Point};
use ogeom_topo::{Filter, Model, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_pipe_down_a_helix_rides_the_public_frames() {
    let mut model = Model::new();
    let curve: ogeom_geom::Curve = HelixCurve::new(Frame::WORLD, 10.0, 5.0, 1.5)
        .unwrap()
        .into();
    let domain = curve.domain();
    let spine = make_edge(&mut model, curve, domain, T).unwrap().shape;
    // The rectangle stands square to the spine at its start, off centre in
    // both axes of the start frame, so a frame turned by any angle moves
    // its corners off their rails.
    let start = rotation_minimizing_frames(&model, &spine, &[0.0], None, T)
        .unwrap()
        .frames[0];
    let (n0, b0) = (start.normal.vector(), start.binormal());
    let offsets = [(1.0, -0.5), (3.0, -0.5), (3.0, 1.5), (1.0, 1.5)];
    let corners: Vec<Point> = offsets
        .iter()
        .map(|&(a, b)| start.at + n0 * a + b0 * b)
        .collect();
    let profile = make_polygon(&mut model, &corners, true, T).unwrap().shape;
    let tolerance = 1e-5;
    let pipe = ogeom_offset::make_pipe_shell(&mut model, &profile, &spine, false, tolerance, T)
        .unwrap()
        .shape;
    let edges = explore(&model, &pipe, Filter::OfType(ShapeType::Edge)).unwrap();
    let curves: Vec<ogeom_geom::Curve> = edges
        .iter()
        .filter_map(|e| placed_edge_curve(&model, e, T).unwrap())
        .map(|(c, _)| c)
        .collect();
    let at: Vec<f64> = (0..=60).map(|k| k as f64 / 60.0).collect();
    let frames = rotation_minimizing_frames(&model, &spine, &at, None, T).unwrap();
    let mut worst = 0.0_f64;
    for f in &frames.frames {
        let (n, b) = (f.normal.vector(), f.binormal());
        for &(oa, ob) in &offsets {
            let predicted = f.at + n * oa + b * ob;
            let off = curves
                .iter()
                .map(|c| project_on_curve(c, predicted, 64, T).unwrap().distance)
                .fold(f64::INFINITY, f64::min);
            worst = worst.max(off);
        }
    }
    eprintln!("pipe rails off the public frames by {worst:e} (skin tolerance {tolerance:e})");
    assert!(worst <= tolerance, "a rail strays {worst} from its frame");
}
