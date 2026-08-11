//! §I of `docs/PLAN.md`: exact geometry that is secretly analytic becomes
//! the analytic thing — and geometry that is not, stays what it is.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::SurfaceGeometry;
use ogeom::math::{Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A drum whose every surface and curve was spelt out as NURBS — the way an
/// exchange file might deliver it — comes back analytic: the wall a cylinder
/// at its exact radius, the caps planes, the volume unchanged to the last
/// bit, and every certificate is a measured worst deviation, not a hope.
#[test]
fn a_nurbsed_drum_confesses_its_cylinder() {
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 12.0, T)
        .unwrap()
        .shape;
    let before = ogeom::algo::volume_properties(&model, &drum, Deflection::default(), T)
        .unwrap()
        .mass;
    let nurbsed = ogeom::algo::to_nurbs(&mut model, &drum, T).unwrap().shape;
    for face in explore_unique(&model, &nurbsed, ShapeType::Face).unwrap() {
        let id = model
            .node(&face)
            .and_then(|n| n.data().as_face())
            .unwrap()
            .surface;
        assert!(
            matches!(
                model.geometry().surface(id),
                Some(SurfaceGeometry::BSpline(_))
            ),
            "the premise: every surface arrives free-form"
        );
    }

    let (built, report) = ogeom::heal::canonical_simplify(&mut model, &nurbsed, 1e-6, T).unwrap();
    assert_eq!(report.simplified.len(), 3, "{report:?}");
    let cylinder = report
        .simplified
        .iter()
        .find_map(|s| match s {
            ogeom::heal::Simplified::Cylinder { radius, worst } => Some((*radius, *worst)),
            _ => None,
        })
        .expect("the wall is a cylinder");
    assert!((cylinder.0 - 5.0).abs() < 1e-9, "radius {}", cylinder.0);
    assert!(cylinder.1 < 1e-12, "certificate {}", cylinder.1);

    let after = ogeom::algo::volume_properties(&model, &built.shape, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!(
        (after - before).abs() < before * 1e-12,
        "{after} against {before}"
    );
    let mut cylinders = 0;
    for face in explore_unique(&model, &built.shape, ShapeType::Face).unwrap() {
        let id = model
            .node(&face)
            .and_then(|n| n.data().as_face())
            .unwrap()
            .surface;
        if matches!(
            model.geometry().surface(id),
            Some(SurfaceGeometry::Cylinder(_))
        ) {
            cylinders += 1;
        }
    }
    assert_eq!(cylinders, 1, "the wall carries the cylinder again");
}

/// A genuinely free-form wall — a skinned loft — refuses to be anything
/// else. The decision is the product, and a wrong yes is a solid that
/// measures nearly right with the wrong surface under every operation after.
#[test]
fn a_free_form_wall_stays_free_form() {
    let mut model = Model::new();
    let profile = |model: &mut Model, z: f64, half: f64| {
        let corners = [
            Point::new(-half, -half, z),
            Point::new(half, -half, z),
            Point::new(half, half, z),
            Point::new(-half, half, z),
        ];
        ogeom::algo::make_polygon(model, &corners, true, T)
            .unwrap()
            .shape
    };
    let a = profile(&mut model, 0.0, 8.0);
    let b = profile(&mut model, 5.0, 6.0);
    let c = profile(&mut model, 10.0, 7.0);
    let solid = ogeom::offset::make_loft_skinned(&mut model, &[a, b, c], 0.5, T)
        .unwrap()
        .shape;
    let (_, report) = ogeom::heal::canonical_simplify(&mut model, &solid, 1e-6, T).unwrap();
    assert!(
        report.simplified.iter().all(|s| !matches!(
            s,
            ogeom::heal::Simplified::Cylinder { .. }
                | ogeom::heal::Simplified::Sphere { .. }
                | ogeom::heal::Simplified::Cone { .. }
        )),
        "nothing curved was claimed: {report:?}"
    );
}
