//! A thread's groove: an M6 tooth space swept seven turns down a helix,
//! and cut from a block.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, make_box, make_edge, make_face, make_polygon, make_wire, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{Curve2d as _, Curve3d as _, HelixCurve, PlaneSurface, Surface as _};
use ogeom::math::{Direction, Frame, Plane, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{EdgeRepr, Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

/// The groove: radial from 2.4 to 3.0 about an axis at (10, 10) running
/// down from z = 11, 0.8 wide at 2.4 and 0.125 at 3.0, pitch 1, its
/// section square to the helix where it starts.
fn groove(model: &mut Model, turns: f64) -> Shape {
    let (pitch, radius) = (1.0, 2.7);
    let frame = Frame::new(Point::new(10.0, 10.0, 11.0), -Direction::Z, Direction::X, T).unwrap();
    let helix = HelixCurve::new(frame, radius, pitch, turns).unwrap();
    let range = helix.domain();
    let edge = make_edge(model, helix.into(), range, T).unwrap().shape;
    let spine = make_wire(model, &[edge], T).unwrap().shape;
    let x = frame.x().vector();
    let tangent = (frame.y().vector() * radius
        + frame.z().vector() * (pitch / core::f64::consts::TAU))
        .normalized(T)
        .unwrap();
    let across = tangent.cross(x).normalized(T).unwrap();
    let axis = frame.origin();
    let corners: Vec<Point> = [(2.4, -0.4), (3.0, -0.0625), (3.0, 0.0625), (2.4, 0.4)]
        .iter()
        .map(|&(r, w)| axis + x * r + across * w)
        .collect();
    let wire = make_polygon(model, &corners, true, T).unwrap().shape;
    let plane = Plane::through(axis + x * radius, Direction::new(tangent, T).unwrap());
    let section = make_face(model, PlaneSurface::new(plane).into(), &[wire], T)
        .unwrap()
        .shape;
    ogeom::offset::make_pipe_shell(model, &section, &spine, false, 1e-3, T)
        .unwrap()
        .shape
}

/// Every edge of the groove lies where each face it bounds says it does:
/// its image on the face, read through the face's surface, is on the edge's
/// curve at the same parameter, to the edge's tolerance. Two walls sharing a
/// rail fit their surfaces at their own paces along the helix, and an image
/// assumed straight on the wall that adopted the rail drifted a fifth of a
/// millimetre off it.
#[test]
fn a_swept_groove_is_where_its_faces_say() {
    let mut model = Model::new();
    let groove = groove(&mut model, 7.0);
    assert!(check(&model, &groove, T).unwrap().is_valid());
    for face in explore_unique(&model, &groove, ShapeType::Face).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        let surface = model.geometry().surface(data.surface).unwrap().clone();
        for edge in explore_unique(&model, &face, ShapeType::Edge).unwrap() {
            let ed = model.node(&edge).unwrap().data().as_edge().unwrap().clone();
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = ed.curve3d() else {
                continue;
            };
            let curve = model.geometry().curve(*curve).unwrap().clone();
            let Some(EdgeRepr::PCurve {
                curve: image,
                range: image_range,
                ..
            }) = ed.pcurve_for(data.surface, edge.location())
            else {
                continue;
            };
            let image = model.geometry().pcurve(*image).unwrap().clone();
            for k in 0..=32 {
                let f = f64::from(k) / 32.0;
                let t = range.0 + (range.1 - range.0) * f;
                let s = image_range.0 + (image_range.1 - image_range.0) * f;
                let uv = image.point_at(s, T).unwrap();
                let ((u0, u1), (v0, v1)) = surface.domain();
                let off = surface
                    .point_at(uv.x.clamp(u0, u1), uv.y.clamp(v0, v1), T)
                    .unwrap()
                    .distance(curve.point_at(t, T).unwrap());
                assert!(
                    off <= ed.tolerance.get() + 1e-6,
                    "an edge's image stands {off:e} off it, against {:e}",
                    ed.tolerance.get()
                );
            }
        }
    }
}

/// The groove cut from a block for one, two and seven turns: valid, and
/// the block less exactly what the groove holds of it.
#[test]
fn a_thread_groove_cuts_from_a_block() {
    for turns in [1.0, 2.0, 7.0] {
        let mut model = Model::new();
        let groove = groove(&mut model, turns);
        let block = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &block, &groove, T)
            .unwrap_or_else(|e| panic!("{turns} turns: {e}"))
            .shape;
        let diagnosis = check(&model, &cut, T).unwrap();
        assert!(diagnosis.is_valid(), "{turns} turns: {diagnosis}");
        let inside = ogeom::boolean::common(&mut model, &block, &groove, T)
            .unwrap()
            .shape;
        let (left, taken) = (volume(&model, &cut), volume(&model, &inside));
        assert!(
            (left + taken - 4000.0).abs() < 4000.0 * 1e-6,
            "{turns} turns: {left} + {taken}"
        );
        assert!(taken > 3.0 * (turns - 1.0), "{turns} turns take {taken}");
    }
}
