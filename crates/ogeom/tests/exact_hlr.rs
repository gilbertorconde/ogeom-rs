//! §13's exact half: silhouettes as the curves they are, visibility asked
//! of the faces rather than of a mesh, and the isoparametric and reflect
//! lines that ride the same construction.
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
    // sides — exactly, at the radius, not to within a chord.
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

    // The rim circles are half visible and half hidden — the drum's own
    // wall stands in the way of the far half — so both lists hold curves
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

/// §D1: a torus has no closed-form silhouette, and it is *walked* rather
/// than refused — the same walk a surface intersection uses, following a
/// different condition.
///
/// Seen along its own axis a torus's outline is two circles, the outer and
/// the inner equators, and those are arithmetic — so the marched answer is
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
            // and not exactness — which is what a marched curve is worth.
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
