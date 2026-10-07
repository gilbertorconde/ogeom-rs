//! Exact hidden-line removal: silhouettes as the curves they are,
//! visibility asked of the faces rather than of a mesh, and the
//! isoparametric and reflect lines that ride the same construction.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_drums_silhouette_is_two_rulings_and_a_balls_is_a_circle() {
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T)
        .unwrap()
        .shape;

    // Seen across its axis, a cylinder's outline is the two rulings at the
    // sides: exactly, at the radius, not to within a chord.
    let found = ogeom::hlr::exact::silhouettes(&model, &drum, Vector::X, T).unwrap();
    assert_eq!(found.len(), 2, "two rulings: {found:?}");
    for ruling in &found {
        for k in 0..=8 {
            let t = ruling.range.0 + (ruling.range.1 - ruling.range.0) * f64::from(k) / 8.0;
            let p = ogeom::geom::Curve3d::point_at(&ruling.curve, t, T).unwrap();
            assert!(
                (p.x.hypot(p.y) - 5.0).abs() < 1e-12 && p.x.abs() < 1e-12,
                "a ruling at the side of the drum: {p:?}"
            );
            // The ends are bisected against the face's own trim, so they
            // land within the confusion tolerance of it rather than on it
            // exactly.
            assert!(
                p.z >= -1e-6 && p.z <= 20.0 + 1e-6,
                "and only where the face is: {p:?}"
            );
        }
    }

    // A ball's is the great circle whose plane the view is normal to.
    let mut model = Model::new();
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 3.0, T)
        .unwrap()
        .shape;
    let found = ogeom::hlr::exact::silhouettes(&model, &ball, Vector::Z, T).unwrap();
    assert_eq!(found.len(), 1, "one circle: {found:?}");
    for k in 0..=16 {
        let t = found[0].range.0 + (found[0].range.1 - found[0].range.0) * f64::from(k) / 16.0;
        let p = ogeom::geom::Curve3d::point_at(&found[0].curve, t, T).unwrap();
        assert!(
            (p.distance(Point::ORIGIN) - 3.0).abs() < 1e-12 && p.z.abs() < 1e-12,
            "the equator seen from above: {p:?}"
        );
    }
}

#[test]
fn the_far_side_of_a_drum_is_hidden_and_the_near_side_is_not() {
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T)
        .unwrap()
        .shape;
    let view = ogeom::hlr::View::looking(-Vector::X, Vector::Z, T).unwrap();

    let drawing =
        ogeom::hlr::exact::project_exact(&model, &drum, &view, Deflection::default(), T).unwrap();
    assert!(
        !drawing.visible.is_empty() && !drawing.hidden.is_empty(),
        "a drum seen from the side has both: {} visible, {} hidden",
        drawing.visible.len(),
        drawing.hidden.len()
    );

    // The rim circles are half visible and half hidden (the drum's own
    // wall stands in the way of the far half), so both lists hold curves
    // that came from model edges.
    let from_edges = |curves: &[ogeom::hlr::DrawnCurve]| {
        curves
            .iter()
            .filter(|c| matches!(c.source, ogeom::hlr::Source::Edge(_)))
            .count()
    };
    assert!(
        from_edges(&drawing.visible) > 0 && from_edges(&drawing.hidden) > 0,
        "the rims are split by the wall between them"
    );
}

#[test]
fn iso_lines_stay_on_their_face_and_reflect_lines_follow_the_light() {
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T)
        .unwrap()
        .shape;
    let wall = explore_unique(&model, &drum, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let ogeom::topo::NodeData::Face(data) = model.node(f).unwrap().data() else {
                return false;
            };
            matches!(
                model.geometry().surface(data.surface),
                Some(ogeom::geom::SurfaceGeometry::Cylinder(_))
            )
        })
        .expect("the wall");

    let lines = ogeom::hlr::exact::iso_curves(&model, &wall, 6, 3, T).unwrap();
    assert!(!lines.is_empty(), "a wall has isoparametrics");
    for line in &lines {
        for p in line {
            assert!(
                (p.x.hypot(p.y) - 5.0).abs() < 1e-9 && p.z >= -1e-9 && p.z <= 20.0 + 1e-9,
                "an isoparametric lies on its own face: {p:?}"
            );
        }
    }

    // A reflect line is the same locus asked of a light instead of an eye,
    // so lighting the drum from the side puts its lines at the sides.
    let lit = ogeom::hlr::exact::reflect_lines(&model, &drum, Vector::Y, T).unwrap();
    assert_eq!(lit.len(), 2, "two, as the light has two sides");
    for ruling in &lit {
        let p = ogeom::geom::Curve3d::point_at(&ruling.curve, ruling.range.0, T).unwrap();
        assert!(
            p.y.abs() < 1e-12 && (p.x.abs() - 5.0).abs() < 1e-12,
            "lit from y, the lines run down x: {p:?}"
        );
    }
}

/// A torus has no closed-form silhouette, and it is *walked* rather than
/// refused: the same walk a surface intersection uses, following a
/// different condition.
///
/// Seen along its own axis a torus's outline is two circles, the outer and
/// the inner equators, and those are arithmetic, so the marched answer is
/// held to them rather than to its own consistency.
#[test]
fn a_torus_silhouette_is_marched_and_lands_on_its_own_equators() {
    let mut model = Model::new();
    let (major, minor) = (6.0, 2.0);
    let ring = ogeom::algo::make_torus(&mut model, Frame::WORLD, major, minor, T)
        .unwrap()
        .shape;

    // Down the axis: the normal is square to the view exactly on the two
    // equators, at radii `major ± minor`, both in the plane `z = 0`.
    let found = ogeom::hlr::exact::silhouettes(&model, &ring, Vector::Z, T).unwrap();
    assert!(
        !found.is_empty(),
        "a torus seen down its axis has an outline"
    );
    // The claim is about the geometry, not about how many pieces the face's
    // trim cuts a closed curve into: every point of every piece is on one of
    // the two equators, and both are drawn.
    let (mut inner, mut outer) = (false, false);
    for silhouette in &found {
        for k in 0..=64 {
            let t = silhouette.range.0
                + (silhouette.range.1 - silhouette.range.0) * f64::from(k) / 64.0;
            let p = ogeom::geom::Curve3d::point_at(&silhouette.curve, t, T).unwrap();
            // The walk is fitted, so the claim is the chord it was walked to
            // and not exactness, which is what a marched curve is worth.
            assert!(p.z.abs() < 1e-4, "in the torus's own plane: {p:?}");
            let r = p.x.hypot(p.y);
            if (r - (major - minor)).abs() < 1e-4 {
                inner = true;
            } else if (r - (major + minor)).abs() < 1e-4 {
                outer = true;
            } else {
                panic!("on neither equator: radius {r}");
            }
        }
    }
    assert!(inner && outer, "both equators are drawn");

    // And the defining property itself, seen from a direction with no closed
    // form at all: the surface's own normal is square to the view at every
    // point of what comes back.
    let oblique = Vector::new(0.3, 0.5, 1.0);
    let found = ogeom::hlr::exact::silhouettes(&model, &ring, oblique, T).unwrap();
    assert!(!found.is_empty(), "an oblique view still has an outline");
    let along = oblique / oblique.magnitude();
    for silhouette in &found {
        for k in 0..=32 {
            let t = silhouette.range.0
                + (silhouette.range.1 - silhouette.range.0) * f64::from(k) / 32.0;
            let p = ogeom::geom::Curve3d::point_at(&silhouette.curve, t, T).unwrap();
            // The point is on the torus, and the torus's normal there is
            // square to the view. Both measured against the surface itself.
            let (centre, radial) = {
                let flat = Point::new(p.x, p.y, 0.0);
                let r = flat.distance(Point::ORIGIN);
                assert!(r > 1e-9, "not on the axis");
                let towards = (flat - Point::ORIGIN) / r;
                (Point::ORIGIN + towards * major, towards)
            };
            assert!(
                (p.distance(centre) - minor).abs() < 1e-3,
                "on the tube: {} against {minor}",
                p.distance(centre)
            );
            let normal = (p - centre) / minor;
            let _ = radial;
            assert!(
                normal.dot(along).abs() < 1e-3,
                "the normal is square to the view: {}",
                normal.dot(along)
            );
        }
    }
}

/// A ball's outline drawn at a chord follows its circle to that chord: its
/// points stand on the radius, and there are enough of them that no chord
/// between neighbours sags past the deflection. Stepped by the straight
/// distance between its range's ends, which a closed circle makes zero,
/// it would be an octagon.
#[test]
fn a_balls_outline_is_drawn_round_not_as_an_octagon() {
    let mut model = Model::new();
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, 3.0, T)
        .unwrap()
        .shape;
    let view = ogeom::hlr::View::looking(-Vector::Z, Vector::Y, T).unwrap();
    let chord = 1e-3;
    let drawing = ogeom::hlr::exact::project_exact(
        &model,
        &ball,
        &view,
        Deflection::with_chord(chord).unwrap(),
        T,
    )
    .unwrap();
    let outline: Vec<_> = drawing
        .visible
        .iter()
        .filter(|c| matches!(c.source, ogeom::hlr::Source::Silhouette))
        .flat_map(|c| {
            c.points
                .windows(2)
                .map(|w| (w[0], w[1]))
                .collect::<Vec<_>>()
        })
        .collect();
    assert!(!outline.is_empty());
    for (a, b) in &outline {
        // A chord of a radius-3 circle sags by length squared over 24.
        let length = ((a.x - b.x).powi(2) + (a.y - b.y).powi(2)).sqrt();
        assert!(length * length / 24.0 <= chord * 1.5, "a chord of {length}");
    }
}

/// A bar under a block that covers its far end, seen from above. The bar's
/// top edge is a line, sampled at its two ends only; it is visible up to
/// the block's side and hidden past it, and the two runs meet there.
#[test]
fn a_straight_edge_half_under_a_block_changes_at_the_blocks_side() {
    let mut model = Model::new();
    let bar = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 1.0, 1.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(
        Point::new(12.0, -5.0, 5.0),
        ogeom::math::Direction::Z,
        ogeom::math::Direction::X,
        T,
    )
    .unwrap();
    let block = ogeom::algo::make_box(&mut model, at, (20.0, 10.0, 2.0), T)
        .unwrap()
        .shape;
    let both = ogeom::algo::make_compound(&mut model, &[bar.clone(), block])
        .unwrap()
        .shape;
    let view = ogeom::hlr::View::looking(-Vector::Z, Vector::Y, T).unwrap();
    let drawing =
        ogeom::hlr::exact::project_exact(&model, &both, &view, Deflection::default(), T).unwrap();

    let top = explore_unique(&model, &bar, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let points =
                ogeom::mesh::polyline_of_edge(&model, e, Deflection::default(), T).unwrap();
            points.len() == 2
                && points
                    .iter()
                    .all(|p| p.y.abs() < 1e-9 && (p.z - 1.0).abs() < 1e-9)
        })
        .unwrap();
    let spans = |curves: &[ogeom::hlr::DrawnCurve]| -> Vec<(f64, f64)> {
        curves
            .iter()
            .filter(|c| matches!(&c.source, ogeom::hlr::Source::Edge(e) if e.node() == top.node()))
            .map(|c| {
                let xs = c.points.iter().map(|p| p.x);
                (
                    xs.clone().fold(f64::INFINITY, f64::min),
                    xs.fold(f64::NEG_INFINITY, f64::max),
                )
            })
            .collect()
    };
    let (visible, hidden) = (spans(&drawing.visible), spans(&drawing.hidden));
    assert_eq!(
        (visible.len(), hidden.len()),
        (1, 1),
        "{visible:?} {hidden:?}"
    );
    assert!(
        visible[0].0.abs() < 1e-9 && (visible[0].1 - 12.0).abs() < 1e-6,
        "{visible:?}"
    );
    assert!(
        (hidden[0].0 - 12.0).abs() < 1e-6 && (hidden[0].1 - 20.0).abs() < 1e-9,
        "{hidden:?}"
    );
}

/// A bar under a spline sheet that covers its middle, seen from above. The
/// sheet has no closed-form piercing, so each ray is asked of it through the
/// seeded intersector over the stretch of the ray inside the sheet's box.
/// The bar's top edge is hidden where the sheet stands over it and visible
/// past the sheet's sides on either end.
#[test]
fn a_straight_edge_under_a_spline_sheet_is_hidden_across_the_sheet() {
    use ogeom::geom::{Curve3d as _, Surface as _};
    let mut model = Model::new();
    let at = Frame::new(
        Point::new(-10.0, 0.0, 0.0),
        ogeom::math::Direction::Z,
        ogeom::math::Direction::X,
        T,
    )
    .unwrap();
    let bar = ogeom::algo::make_box(&mut model, at, (20.0, 1.0, 1.0), T)
        .unwrap()
        .shape;

    // A saddle five units over the bar, spanning x and y in [-5, 5].
    let n = 11;
    let rows: Vec<Vec<Point>> = (0..n)
        .map(|j| {
            let y = -5.0 + 10.0 * f64::from(j) / f64::from(n - 1);
            (0..n)
                .map(|i| {
                    let x = -5.0 + 10.0 * f64::from(i) / f64::from(n - 1);
                    Point::new(x, y, 5.0 + (x * x - y * y) / 20.0)
                })
                .collect()
        })
        .collect();
    let surface = ogeom::geom::fit::fit_surface_grid(&rows, 3, 1e-6, T)
        .unwrap()
        .curve;
    let ((u0, u1), (v0, v1)) = surface.domain();
    let corners: Vec<_> = [(u0, v0), (u1, v0), (u1, v1), (u0, v1)]
        .iter()
        .map(|(u, v)| {
            let p = surface.point_at(*u, *v, T).unwrap();
            ogeom::algo::make_vertex(&mut model, p).shape
        })
        .collect();
    let mut iso = |curve: ogeom::geom::BSplineCurve, from: usize, to: usize| {
        let range = curve.domain();
        ogeom::algo::make_edge_between(
            &mut model,
            curve.into(),
            range,
            &corners[from],
            &corners[to],
            T,
        )
        .unwrap()
        .shape
    };
    let south = iso(surface.iso_v_curve(v0, T).unwrap(), 0, 1);
    let east = iso(surface.iso_u_curve(u1, T).unwrap(), 1, 2);
    let north = iso(surface.iso_v_curve(v1, T).unwrap(), 3, 2);
    let west = iso(surface.iso_u_curve(u0, T).unwrap(), 0, 3);
    let wires = [vec![south, east, north.reversed(), west.reversed()]];
    let sheet = ogeom::algo::make_face_with_pcurves(&mut model, surface.into(), &wires, T)
        .unwrap()
        .shape;
    let both = ogeom::algo::make_compound(&mut model, &[bar.clone(), sheet])
        .unwrap()
        .shape;
    let view = ogeom::hlr::View::looking(-Vector::Z, Vector::Y, T).unwrap();
    let drawing =
        ogeom::hlr::exact::project_exact(&model, &both, &view, Deflection::default(), T).unwrap();

    let top = explore_unique(&model, &bar, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let points =
                ogeom::mesh::polyline_of_edge(&model, e, Deflection::default(), T).unwrap();
            points.len() == 2
                && points
                    .iter()
                    .all(|p| p.y.abs() < 1e-9 && (p.z - 1.0).abs() < 1e-9)
        })
        .unwrap();
    let spans = |curves: &[ogeom::hlr::DrawnCurve]| -> Vec<(f64, f64)> {
        let mut out: Vec<(f64, f64)> = curves
            .iter()
            .filter(|c| matches!(&c.source, ogeom::hlr::Source::Edge(e) if e.node() == top.node()))
            .map(|c| {
                let xs = c.points.iter().map(|p| p.x);
                (
                    xs.clone().fold(f64::INFINITY, f64::min),
                    xs.fold(f64::NEG_INFINITY, f64::max),
                )
            })
            .collect();
        out.sort_by(|a, b| a.0.total_cmp(&b.0));
        out
    };
    let (visible, hidden) = (spans(&drawing.visible), spans(&drawing.hidden));
    assert_eq!(
        (visible.len(), hidden.len()),
        (2, 1),
        "{visible:?} {hidden:?}"
    );
    // The sheet's edge over the bar is the fitted border at x = -5 and
    // x = 5, to the fit's tolerance.
    let near = |a: f64, b: f64| (a - b).abs() < 1e-4;
    assert!(
        near(visible[0].0, -10.0) && near(visible[0].1, -5.0),
        "{visible:?}"
    );
    assert!(
        near(visible[1].0, 5.0) && near(visible[1].1, 10.0),
        "{visible:?}"
    );
    assert!(
        near(hidden[0].0, -5.0) && near(hidden[0].1, 5.0),
        "{hidden:?}"
    );
}
