//! Exact geometry that is secretly analytic becomes
//! the analytic thing, and geometry that is not, stays what it is.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/walks.rs"]
mod walks;

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
/// Fusing a post onto a slab splits the post's rims into two half circles.
/// A half's pcurve written as the straight chart segment between its
/// projected endpoints goes wrong where the seam vertex projects to either
/// column: the front half draws mirrored, and the wall's ring loses its
/// far side. And the fit's own slop must not widen a converted edge past
/// the vertices that bound it, which the checker refuses. Both the column
/// and the width are decided by the edge's own interior, and the converted
/// solid is as closed and as large as the analytic one.
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
        assert_eq!(walks::edges_walked_one_way(&model, &built.shape), 0);
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
    assert_eq!(walks::edges_walked_one_way(&model, &restated), 0);
    assert_eq!(walks::edges_walked_one_way(&model, &built.shape), 0);
    let after = ogeom::algo::volume_properties(&model, &built.shape, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!(
        (after - before).abs() < before * 1e-6,
        "{after} against {before}"
    );
}

/// A plane patch with a bump confined to its outer strip, a twentieth of
/// the chart wide, is not a plane: the bump stands 0.05 off, and a
/// recognition that never looks at the strip would call it flat and state
/// a certificate of nothing.
#[test]
fn a_bump_in_the_outer_strip_is_not_a_plane() {
    use ogeom::geom::{BSplineSurface, Surface as _};
    use ogeom::math::{ControlGrid, KnotVector};
    const SPANS: usize = 40;
    let n = SPANS + 3;
    let mut points = Vec::with_capacity(n * 4);
    for i in 0..n {
        for j in 0..4 {
            #[allow(clippy::cast_precision_loss)]
            let x = 10.0 * i as f64 / (n - 1) as f64;
            let y = 10.0 * f64::from(j) / 3.0;
            let z = if i == 1 { 0.1 } else { 0.0 };
            points.push(Point::new(x, y, z));
        }
    }
    let patch = BSplineSurface::new(
        KnotVector::clamped_uniform(3, n).unwrap(),
        KnotVector::clamped_uniform(3, 4).unwrap(),
        &ControlGrid::new(points, n, 4).unwrap(),
        T,
    )
    .unwrap();
    let surface = SurfaceGeometry::BSpline(patch);
    let mut bump = 0.0_f64;
    for i in 0..=400 {
        for j in 0..=10 {
            let p = surface
                .point_at(f64::from(i) / 400.0, f64::from(j) / 10.0, T)
                .unwrap();
            bump = bump.max(p.z.abs());
        }
    }
    assert!(bump > 0.04, "the premise: a bump of {bump}");
    let found = ogeom::heal::recognize_surface(&surface, 1e-3, T).unwrap();
    assert!(
        found.is_none(),
        "a {bump} bump read as {found:?} at a tolerance of 1e-3"
    );
}

/// A flat patch bounded by a spline that runs straight but for a bump
/// between the samples a line is verified at keeps its bumped edge: the
/// face is recognized as the plane, and its boundary keeps its bump.
#[test]
fn a_bumped_boundary_is_not_a_line() {
    use ogeom::geom::{BSplineSurface, Curve, Curve3d as _};
    use ogeom::math::{ControlGrid, KnotVector};
    use ogeom::topo::EdgeRepr;
    const SPANS: usize = 80;
    let mut model = Model::new();
    let patch = BSplineSurface::new(
        KnotVector::clamped_uniform(1, 2).unwrap(),
        KnotVector::clamped_uniform(1, 2).unwrap(),
        &ControlGrid::new(
            vec![
                Point::new(-1.0, -1.0, 0.0),
                Point::new(-1.0, 11.0, 0.0),
                Point::new(11.0, -1.0, 0.0),
                Point::new(11.0, 11.0, 0.0),
            ],
            2,
            2,
        )
        .unwrap(),
        T,
    )
    .unwrap();
    let control: Vec<Point> = (0..SPANS + 3)
        .map(|i| {
            #[allow(clippy::cast_precision_loss)]
            let x = 10.0 * i as f64 / (SPANS + 2) as f64;
            Point::new(x, if i == 8 { 0.1 } else { 0.0 }, 0.0)
        })
        .collect();
    let bumped = Curve::BSpline(
        ogeom::geom::BSplineCurve::new(
            KnotVector::clamped_uniform(3, control.len()).unwrap(),
            control,
            T,
        )
        .unwrap(),
    );
    let mut bump = 0.0_f64;
    for i in 0..=4000 {
        bump = bump.max(bumped.point_at(f64::from(i) / 4000.0, T).unwrap().y);
    }
    assert!(bump > 0.04, "the premise: a bump of {bump}");
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(10.0, 0.0, 0.0),
        Point::new(10.0, 10.0, 0.0),
        Point::new(0.0, 10.0, 0.0),
    ];
    let v: Vec<_> = corners
        .iter()
        .map(|p| ogeom::algo::make_vertex(&mut model, *p).shape)
        .collect();
    let mut edges = vec![
        ogeom::algo::make_edge_between(&mut model, bumped, (0.0, 1.0), &v[0], &v[1], T)
            .unwrap()
            .shape,
    ];
    for i in 1..4 {
        let (a, b) = (corners[i], corners[(i + 1) % 4]);
        let line = ogeom::geom::LineCurve::segment(a, b, T).unwrap();
        edges.push(
            ogeom::algo::make_edge_between(
                &mut model,
                Curve::from(line),
                (0.0, a.distance(b)),
                &v[i],
                &v[(i + 1) % 4],
                T,
            )
            .unwrap()
            .shape,
        );
    }
    let wire = ogeom::algo::make_wire(&mut model, &edges, T).unwrap().shape;
    let face = ogeom::algo::make_face(&mut model, SurfaceGeometry::BSpline(patch), &[wire], T)
        .unwrap()
        .shape;
    let (built, report) = ogeom::heal::canonical_simplify(&mut model, &face, 1e-3, T).unwrap();
    assert_eq!(report.simplified.len(), 1, "{report:?}");
    let face = explore_unique(&model, &built.shape, ShapeType::Face)
        .unwrap()
        .remove(0);
    for edge in explore_unique(&model, &face, ShapeType::Edge).unwrap() {
        let data = model.node(&edge).unwrap().data().as_edge().unwrap();
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let curve = model.geometry().curve(*curve).unwrap();
        let mid = curve.point_at(f64::midpoint(range.0, range.1), T).unwrap();
        if mid.y.abs() > 1.0 {
            continue;
        }
        // The bottom edge: how far the bump stands from it.
        let mut miss = 0.0_f64;
        for i in 0..=400 {
            let t = range.0 + (range.1 - range.0) * f64::from(i) / 400.0;
            let p = curve.point_at(t, T).unwrap();
            miss = miss.max(p.y.abs());
        }
        let kept = miss > bump * 0.9;
        assert!(
            kept || data.tolerance.get() >= bump,
            "the bump of {bump} vanished into a line stating {}",
            data.tolerance.get()
        );
    }
    let diagnosis = ogeom::algo::check(&model, &face, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
}
