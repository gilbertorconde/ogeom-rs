//! Booleans between whole revolutions about one axis, their seams on one
//! half-plane turned any way about it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_face, make_polygon, make_revolution, transformed, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::PlaneSurface;
use ogeom::math::{Axis, Direction, Frame, Plane, Point, Transform};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// A rectangle y in `ys`, z in `zs` in the YZ plane, turned about Z by
/// `seam` and revolved a whole turn about Z.
fn revolved(model: &mut Model, ys: (f64, f64), zs: (f64, f64), seam: f64) -> Shape {
    let corners = [(ys.0, zs.0), (ys.1, zs.0), (ys.1, zs.1), (ys.0, zs.1)]
        .map(|(y, z)| Point::new(0.0, y, z));
    let wire = make_polygon(model, &corners, true, T).unwrap().shape;
    let plane = Plane::new(Frame::new(Point::ORIGIN, Direction::X, Direction::Y, T).unwrap());
    let face = make_face(model, PlaneSurface::new(plane).into(), &[wire], T)
        .unwrap()
        .shape;
    let axis = Axis::new(Point::ORIGIN, Direction::Z);
    let turned = transformed(model, &face, Transform::rotation(axis, seam))
        .unwrap()
        .shape;
    make_revolution(model, &turned, axis, 2.0 * PI, T)
        .unwrap()
        .shape
}

/// A cylinder and a collar over its top rim, both seams on one half-plane:
/// fused and cut to the volumes their profiles give, whichever way the
/// seams stand.
#[test]
fn whole_revolutions_with_seams_together_fuse_and_cut() {
    // The cylinder r 10, z 0..20; the collar y 8..14, z 18..24.
    let cylinder = PI * 100.0 * 20.0;
    let collar = PI * (196.0 - 64.0) * 6.0;
    let overlap = PI * (100.0 - 64.0) * 2.0;
    for seam in [0.0, 0.5, 1.0] {
        let mut model = Model::new();
        let a = revolved(&mut model, (0.0, 10.0), (0.0, 20.0), seam);
        let b = revolved(&mut model, (8.0, 14.0), (18.0, 24.0), seam);
        let fused = ogeom::boolean::fuse(&mut model, &a, &b, T)
            .unwrap_or_else(|e| panic!("fuse at seam {seam}: {e}"))
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &a, &b, T)
            .unwrap_or_else(|e| panic!("cut at seam {seam}: {e}"))
            .shape;
        for (name, result, want) in [
            ("fuse", &fused, cylinder + collar - overlap),
            ("cut", &cut, cylinder - overlap),
        ] {
            let diagnosis = ogeom::algo::check(&model, result, T).unwrap();
            assert!(diagnosis.is_valid(), "{name} at seam {seam}: {diagnosis}");
            let v = volume(&model, result);
            assert!(
                (v - want).abs() < want * 1e-6,
                "{name} at seam {seam}: {v} against {want}"
            );
        }
    }
}
