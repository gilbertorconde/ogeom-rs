//! One part, stated in micrometres, millimetres, inches, metres and
//! kilometres, builds, cuts, checks, measures and reads back the same: the
//! unit scale moves every linear tolerance with the numbers.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_box, make_cylinder, volume_properties};
use ogeom::core::Tolerances;
use ogeom::io::native::{WriteOptions, read, write};
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::Model;

/// A 20 mm cube with a 4 mm square notch through one corner and a bore of
/// radius 3 mm through its middle, every length divided by `mm_per_unit`.
/// Returns the volume read back from a native round trip, over the closed
/// form, both in cubic model units.
fn drilled_cube(mm_per_unit: f64) -> (f64, f64) {
    let tol = Tolerances::with_scale(mm_per_unit).unwrap();
    let u = |mm: f64| mm / mm_per_unit;
    let mut model = Model::with_tolerances(tol);
    let cube = make_box(&mut model, Frame::WORLD, (u(20.0), u(20.0), u(20.0)), tol)
        .unwrap()
        .shape;
    let at = Frame::new(
        Point::new(u(16.0), u(16.0), u(-1.0)),
        Direction::Z,
        Direction::X,
        tol,
    )
    .unwrap();
    let notch = make_box(&mut model, at, (u(8.0), u(8.0), u(22.0)), tol)
        .unwrap()
        .shape;
    let at = Frame::new(
        Point::new(u(10.0), u(10.0), u(-1.0)),
        Direction::Z,
        Direction::X,
        tol,
    )
    .unwrap();
    let drill = make_cylinder(&mut model, at, u(3.0), u(22.0), tol)
        .unwrap()
        .shape;
    let notched = ogeom::boolean::cut(&mut model, &cube, &notch, tol)
        .unwrap()
        .shape;
    let part = ogeom::boolean::cut(&mut model, &notched, &drill, tol)
        .unwrap()
        .shape;
    let diagnosis = check(&model, &part, tol).unwrap();
    assert!(
        diagnosis.is_valid(),
        "{mm_per_unit} mm per unit: {diagnosis:?}"
    );

    let text = write(&model, &[part], WriteOptions::default()).unwrap();
    let (back, roots) = read(&text).unwrap();
    assert_eq!(back.tolerances(), tol);
    let diagnosis = check(&back, &roots[0], tol).unwrap();
    assert!(
        diagnosis.is_valid(),
        "{mm_per_unit} mm per unit, read back: {diagnosis:?}"
    );
    let chord = Deflection::with_chord(u(1e-3)).unwrap();
    let volume = volume_properties(&back, &roots[0], chord, tol)
        .unwrap()
        .mass;
    let want = u(20.0).powi(3)
        - u(4.0) * u(4.0) * u(20.0)
        - core::f64::consts::PI * u(3.0) * u(3.0) * u(20.0);
    (volume, want)
}

#[test]
fn one_part_measures_the_same_in_every_unit() {
    for mm_per_unit in [1e-3, 1.0, 25.4, 1e3, 1e6] {
        let (volume, want) = drilled_cube(mm_per_unit);
        assert!(
            ((volume - want) / want).abs() < 1e-12,
            "{mm_per_unit} mm per unit: {volume} against {want}"
        );
    }
}

/// The greatest tolerance an edge of `shape` states, in millimetres.
fn widest_edge_mm(model: &Model, shape: &ogeom::topo::Shape, mm_per_unit: f64) -> f64 {
    ogeom::topo::explore_unique(model, shape, ogeom::topo::ShapeType::Edge)
        .unwrap()
        .iter()
        .filter_map(|e| model.node(e)?.data().as_edge().map(|d| d.tolerance.get()))
        .fold(0.0, f64::max)
        * mm_per_unit
}

/// A slab with a post leaning twenty degrees out of it, fused, and the
/// post's elliptical seat on the top filleted at 1 mm: a marched blend.
/// Returns the greatest edge tolerance of the result, in millimetres.
fn filleted_leaning_post(mm_per_unit: f64) -> f64 {
    let tol = Tolerances::with_scale(mm_per_unit).unwrap();
    let u = |mm: f64| mm / mm_per_unit;
    let mut model = Model::with_tolerances(tol);
    let slab = make_box(&mut model, Frame::WORLD, (u(20.0), u(20.0), u(2.0)), tol)
        .unwrap()
        .shape;
    let lean = 20.0_f64.to_radians();
    let axis = ogeom::math::Vector::new(lean.sin(), 0.0, lean.cos());
    let frame = Frame::new(
        Point::new(u(10.0), u(10.0), u(-1.0)),
        Direction::new(axis, tol).unwrap(),
        Direction::from_cross(axis, ogeom::math::Vector::Y, tol).unwrap(),
        tol,
    )
    .unwrap();
    let post = make_cylinder(&mut model, frame, u(3.0), u(10.0), tol)
        .unwrap()
        .shape;
    let joined = ogeom::boolean::fuse(&mut model, &slab, &post, tol)
        .unwrap()
        .shape;
    let seat = ogeom::topo::explore_unique(&model, &joined, ogeom::topo::ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let elliptical = model
                .node(e)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.curve3d())
                .and_then(|r| match r {
                    ogeom::topo::EdgeRepr::Curve3d { curve, .. } => model.geometry().curve(*curve),
                    _ => None,
                })
                .is_some_and(|c| matches!(c, ogeom::geom::Curve::Ellipse(_)));
            elliptical
                && ogeom::algo::edge_vertices(&model, e)
                    .unwrap()
                    .and_then(|(a, _)| model.node(&a)?.data().as_vertex().map(|d| d.point))
                    .is_some_and(|p| (p.z - u(2.0)).abs() < u(1e-6))
        })
        .expect("the fused part has its elliptical seat on the top");
    let round = ogeom::fillet::fillet_edge(&mut model, &joined, &seat, u(1.0), tol)
        .unwrap()
        .shape;
    let diagnosis = check(&model, &round, tol).unwrap();
    assert!(
        diagnosis.is_valid(),
        "{mm_per_unit} mm per unit: {diagnosis:?}"
    );
    widest_edge_mm(&model, &round, mm_per_unit)
}

/// The marched band's edges are widened to its fit target, two tenths of a
/// micron in every unit. A part in metres is refused by the blend engine
/// before the band is built, so that unit is not in the list.
#[test]
fn a_marched_fillet_states_the_same_tolerance_in_every_unit() {
    let at_mm = filleted_leaning_post(1.0);
    for mm_per_unit in [1e-3, 25.4] {
        let widest = filleted_leaning_post(mm_per_unit);
        assert!(
            (widest - at_mm).abs() <= at_mm * 1e-3,
            "{mm_per_unit} mm per unit: edges stated to {widest:e} mm, \
             against {at_mm:e} mm for the part in millimetres"
        );
    }
}

/// A 20 by 8 by 10 mm slab whose front wall is a sine-wave spline extruded
/// along z, that wall drafted by a tenth of a radian about the bottom.
/// Returns the drafted wall's control net size.
fn drafted_spline_wall(mm_per_unit: f64) -> (usize, usize) {
    use ogeom::geom::{Curve, Curve3d as _, LineCurve, PlaneSurface};
    let tol = Tolerances::with_scale(mm_per_unit).unwrap();
    let u = |mm: f64| mm / mm_per_unit;
    let mut model = Model::with_tolerances(tol);
    let points: Vec<Point> = (0..=16)
        .map(|i| {
            let x = 20.0 * f64::from(i) / 16.0;
            Point::new(u(x), u((x * 0.4).sin() * 1.5), 0.0)
        })
        .collect();
    let spline = Curve::BSpline(
        ogeom::geom::fit::fit_points(&points, 3, u(1e-6), tol)
            .unwrap()
            .curve,
    );
    let corners = [
        points[0],
        points[16],
        Point::new(u(20.0), u(-8.0), 0.0),
        Point::new(0.0, u(-8.0), 0.0),
    ];
    let v: Vec<_> = corners
        .iter()
        .map(|p| ogeom::algo::make_vertex(&mut model, *p).shape)
        .collect();
    let domain = spline.domain();
    let front = ogeom::algo::make_edge_between(&mut model, spline, domain, &v[0], &v[1], tol)
        .unwrap()
        .shape;
    let mut sides = Vec::new();
    for (i, j) in [(1, 2), (2, 3), (3, 0)] {
        let line = Curve::Line(LineCurve::segment(corners[i], corners[j], tol).unwrap());
        let domain = line.domain();
        sides.push(
            ogeom::algo::make_edge_between(&mut model, line, domain, &v[i], &v[j], tol)
                .unwrap()
                .shape,
        );
    }
    let plane = ogeom::math::Plane::through(Point::ORIGIN, Direction::Z);
    let face = ogeom::algo::make_face_with_pcurves(
        &mut model,
        PlaneSurface::over(plane, (u(-40.0), u(40.0)), (u(-40.0), u(40.0)))
            .unwrap()
            .into(),
        &[vec![
            sides[2].reversed(),
            sides[1].reversed(),
            sides[0].reversed(),
            front.reversed(),
        ]],
        tol,
    )
    .unwrap()
    .shape;
    let slab = ogeom::algo::make_prism(
        &mut model,
        &face,
        ogeom::math::Vector::new(0.0, 0.0, u(10.0)),
        tol,
    )
    .unwrap()
    .shape;
    let faces = ogeom::topo::explore(
        &model,
        &slab,
        ogeom::topo::Filter::OfType(ogeom::topo::ShapeType::Face),
    )
    .unwrap();
    let surface_of = |model: &Model, f: &ogeom::topo::Shape| {
        model
            .node(f)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .cloned()
    };
    let wall = faces
        .iter()
        .find(|f| {
            matches!(
                surface_of(&model, f),
                Some(ogeom::geom::SurfaceGeometry::Extrusion(_))
            )
        })
        .expect("the slab has an extruded wall")
        .clone();
    let drafted = ogeom::offset::apply_draft(
        &mut model,
        &slab,
        std::slice::from_ref(&wall),
        plane,
        Direction::Z,
        0.1,
        tol,
    )
    .unwrap()
    .shape;
    let diagnosis = check(&model, &drafted, tol).unwrap();
    assert!(
        diagnosis.is_valid(),
        "{mm_per_unit} mm per unit: {diagnosis:?}"
    );
    ogeom::topo::explore(
        &model,
        &drafted,
        ogeom::topo::Filter::OfType(ogeom::topo::ShapeType::Face),
    )
    .unwrap()
    .iter()
    .find_map(|f| match surface_of(&model, f)? {
        ogeom::geom::SurfaceGeometry::BSpline(b) => Some((b.grid().u_count(), b.grid().v_count())),
        _ => None,
    })
    .expect("the drafted wall is fitted")
}

#[test]
fn a_drafted_spline_wall_is_fitted_the_same_in_every_unit() {
    let at_mm = drafted_spline_wall(1.0);
    for mm_per_unit in [1e-3, 25.4, 1e3] {
        assert_eq!(
            drafted_spline_wall(mm_per_unit),
            at_mm,
            "{mm_per_unit} mm per unit: the drafted wall's control net"
        );
    }
}
