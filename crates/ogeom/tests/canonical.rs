//! §I of `docs/PLAN.md`: exact geometry that is secretly analytic becomes
//! the analytic thing, and geometry that is not, stays what it is.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::SurfaceGeometry;
use ogeom::math::{Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A drum whose every surface and curve was spelt out as NURBS (the way an
/// exchange file might deliver it) comes back analytic: the wall a cylinder
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

/// A genuinely free-form wall (a skinned loft) refuses to be anything
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

/// A fused post converts to patches that still close.
///
/// Fusing a post onto a slab splits the post's rims into two half circles,
/// and the conversion wrote the second half's pcurve as the straight chart
/// segment between its projected endpoints: the seam vertex projecting to
/// either column, it drew the front half mirrored, and the wall's ring lost
/// its far side. The fit's own slop widened each converted edge past the
/// vertices that bound it, which the checker refuses. Both are decided by
/// the edge's own interior now, and the converted solid is as closed and
/// as large as the analytic one.
#[test]
fn a_fused_post_converts_to_patches_that_close() {
    let mut model = Model::new();
    let slab = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 2.0), T).unwrap();
    let frame = Frame::new(
        Point::new(10.0, 10.0, 2.0),
        ogeom::math::Direction::Z,
        ogeom::math::Direction::X,
        T,
    )
    .unwrap();
    let post = ogeom::algo::make_cylinder(&mut model, frame, 3.0, 8.0, T).unwrap();
    let joined = ogeom::boolean::fuse(&mut model, &slab.shape, &post.shape, T).unwrap();
    // Measured at a fine chord: a rational patch meshed at the default
    // deflection sits inside its drum by more than the check tolerates.
    let fine = Deflection {
        chord: 1e-3,
        ..Deflection::default()
    };
    let before = ogeom::algo::volume_properties(&model, &joined.shape, fine, T)
        .unwrap()
        .mass;

    let converted = ogeom::algo::to_nurbs(&mut model, &joined.shape, T).unwrap();
    let diagnosis = ogeom::algo::check(&model, &converted.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let mesh =
        ogeom::mesh::triangulate(&model, &converted.shape, Deflection::default(), T).unwrap();
    assert!(mesh.is_closed(), "the converted post draws closed");
    let after = ogeom::algo::volume_properties(&model, &converted.shape, fine, T)
        .unwrap()
        .mass;
    assert!(
        (after - before).abs() < before * 5e-3,
        "converted {after} against {before}"
    );
}

/// A drilled block and its mirror image, spelt out as NURBS: the mirror
/// turns every spline's own normal, so the mirrored bore's spline has its
/// normal towards the axis. Recognised, both bore walls are cylinders,
/// whose normal is always away from the axis, and both parts stay the
/// right way out.
#[test]
fn a_bore_spelt_either_way_round_stays_the_right_way_out() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 6.0), T)
        .unwrap()
        .shape;
    let seat = Frame::new(
        Point::new(5.0, 5.0, -1.0),
        ogeom::math::Direction::Z,
        ogeom::math::Direction::X,
        T,
    )
    .unwrap();
    let bore = ogeom::algo::make_cylinder(&mut model, seat, 2.0, 8.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(&mut model, &block, &bore, T)
        .unwrap()
        .shape;
    let mirror = ogeom::math::Transform::plane_mirror(Point::ORIGIN, ogeom::math::Direction::X);
    let mirrored = model.placed(&drilled, mirror);
    // A mirror keeps the volume: both parts are held to the original's.
    let before = ogeom::algo::volume_properties(&model, &drilled, Deflection::default(), T)
        .unwrap()
        .mass;
    for part in [drilled, mirrored] {
        let nurbsed = ogeom::algo::to_nurbs(&mut model, &part, T).unwrap().shape;
        let (built, _) = ogeom::heal::canonical_simplify(&mut model, &nurbsed, 1e-6, T).unwrap();
        let diagnosis = ogeom::algo::check(&model, &built.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{diagnosis}");
        let after = ogeom::algo::volume_properties(&model, &built.shape, Deflection::default(), T)
            .unwrap()
            .mass;
        assert!(
            (after - before).abs() < before * 1e-6,
            "{after} against {before}"
        );
    }
}

/// A bore restated as a spline running the other way round, so its own
/// normal points at the axis: recognised as a cylinder, whose normal points
/// away, the face turns with it and the part stays the right way out.
#[test]
fn a_bore_spline_facing_its_axis_is_recognised_the_right_way_out() {
    use ogeom::geom::BSplineSurface;
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 6.0), T)
        .unwrap()
        .shape;
    let seat = Frame::new(
        Point::new(5.0, 5.0, -1.0),
        ogeom::math::Direction::Z,
        ogeom::math::Direction::X,
        T,
    )
    .unwrap();
    let bore = ogeom::algo::make_cylinder(&mut model, seat, 2.0, 8.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(&mut model, &block, &bore, T)
        .unwrap()
        .shape;
    let before = ogeom::algo::volume_properties(&model, &drilled, Deflection::default(), T)
        .unwrap()
        .mass;
    let backwards =
        |s: &SurfaceGeometry| -> ogeom::core::OgeomResult<Option<(SurfaceGeometry, bool)>> {
            let SurfaceGeometry::Cylinder(_) = s else {
                return Ok(None);
            };
            let spline = s.to_bspline(T)?;
            let grid = spline.grid();
            let (nu, nv) = (grid.u_count(), grid.v_count());
            let mut points = Vec::with_capacity(nu * nv);
            for i in (0..nu).rev() {
                for j in 0..nv {
                    points.push(grid.get(i, j).unwrap());
                }
            }
            let turned = BSplineSurface::rational(
                spline.u_knots().reversed(),
                spline.v_knots().clone(),
                ogeom::math::ControlGrid::new(points, nu, nv)?,
            )?;
            Ok(Some((turned.into(), true)))
        };
    let keep = |_: &ogeom::geom::Curve, _: (f64, f64)| Ok(None);
    let restated = ogeom::algo::restate_geometry(&mut model, &drilled, &backwards, &keep, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &restated, T).unwrap().is_valid());
    let (built, report) = ogeom::heal::canonical_simplify(&mut model, &restated, 1e-6, T).unwrap();
    assert!(
        report
            .simplified
            .iter()
            .any(|s| matches!(s, ogeom::heal::Simplified::Cylinder { .. })),
        "{report:?}"
    );
    let diagnosis = ogeom::algo::check(&model, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let after = ogeom::algo::volume_properties(&model, &built.shape, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!(
        (after - before).abs() < before * 1e-6,
        "{after} against {before}"
    );
}
