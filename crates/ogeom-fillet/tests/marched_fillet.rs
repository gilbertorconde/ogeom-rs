//! The marched blend carried to topology: fillets on seats no closed form
//! speaks.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Point};
use ogeom_topo::{Filter, ShapeType, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &ogeom_topo::Model, shape: &ogeom_topo::Shape, chord: f64) -> f64 {
    ogeom_algo::volume_properties(
        model,
        shape,
        ogeom_mesh::Deflection {
            chord,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

#[test]
fn a_fillet_on_a_spline_edge_rides_the_marched_blend() {
    // A cylinder leaning twenty degrees out of a slab: the seat is the
    // oblique ellipse no closed-form blend recognises, concave all the way
    // round, so the wedge fuses.
    let mut model = ogeom_topo::Model::new();
    let slab = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 2.0), T).unwrap();
    let lean = 20.0_f64.to_radians();
    let axis = ogeom_math::Vector::new(lean.sin(), 0.0, lean.cos());
    let frame = Frame::new(
        Point::new(10.0, 10.0, -1.0),
        ogeom_math::Direction::new(axis, T).unwrap(),
        ogeom_math::Direction::from_cross(axis, ogeom_math::Vector::Y, T).unwrap(),
        T,
    )
    .unwrap();
    let post = ogeom_algo::make_cylinder(&mut model, frame, 3.0, 10.0, T).unwrap();
    let joined = ogeom_bool::fuse(&mut model, &slab.shape, &post.shape, T).unwrap();

    // The seat: the elliptical edge where the leaning post meets the top —
    // the post pierces the slab, so the bottom carries a second ellipse the
    // test must not grab.
    let edge = explore_unique(&model, &joined.shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let on_ellipse = model
                .node(e)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.curve3d())
                .and_then(|r| {
                    let ogeom_topo::EdgeRepr::Curve3d { curve, .. } = r else {
                        return None;
                    };
                    model.geometry().curve(*curve)
                })
                .is_some_and(|c| matches!(c, ogeom_geom::Curve::Ellipse(_)));
            on_ellipse
                && ogeom_algo::edge_vertices(&model, e)
                    .unwrap()
                    .is_some_and(|(a, _)| {
                        model
                            .node(&a)
                            .and_then(|n| n.data().as_vertex().map(|d| d.point))
                            .is_some_and(|p| (p.z - 2.0).abs() < 1e-6)
                    })
        })
        .expect("the joined part has its elliptical seat");

    let radius = 1.0;
    let result = ogeom_fillet::fillet_edge(&mut model, &joined.shape, &edge, radius, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!(result.history.is_deleted(&edge));

    // No closed form prices the added ring, so hold it two ways: it adds
    // rather than removes, it stays under a computable over-bound, and the
    // measurement converges as the mesh refines.
    let seam_length = 2.0 * core::f64::consts::PI * 3.0 / lean.cos().sqrt();
    let over_bound = radius * radius * seam_length;
    let mut previous_error = f64::INFINITY;
    let mut held = None;
    for chord in [1e-2, 1e-3] {
        let measured = volume(&model, &result.shape, chord);
        let added = measured - volume(&model, &joined.shape, chord);
        assert!(
            added > 0.0 && added < over_bound,
            "the fillet added {added}, outside (0, {over_bound})"
        );
        if let Some(prior) = held {
            let drift: f64 = measured - prior;
            assert!(
                drift.abs() < previous_error.abs().max(1e-6),
                "refinement should settle: {prior} then {measured}"
            );
            previous_error = drift;
        }
        held = Some(measured);
    }

    // The blend face is a genuine fitted surface, and the tangency held.
    let has_spline_face = explore(&model, &result.shape, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .any(|f| {
            model
                .node(&f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::BSpline(_)))
        });
    assert!(has_spline_face, "the blend face rides its fitted surface");
}

#[test]
fn the_seam_of_a_branch_cylinder_gains_a_marched_fillet() {
    // A branch standing on a larger drum: the seam is the fitted saddle
    // curve two cylinders trace, closed on itself, concave all the way
    // round — the marched fillet's home ground.
    let mut model = ogeom_topo::Model::new();
    let main = Frame::new(
        Point::new(-20.0, 0.0, 0.0),
        ogeom_math::Direction::X,
        ogeom_math::Direction::Y,
        T,
    )
    .unwrap();
    let drum = ogeom_algo::make_cylinder(&mut model, main, 10.0, 40.0, T).unwrap();
    let branch = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T).unwrap();
    let joined = ogeom_bool::fuse(&mut model, &drum.shape, &branch.shape, T).unwrap();

    // The seat: an arc of the fitted seam where the branch meets the drum —
    // the boolean splits the loop, but every arc rides the one closed curve
    // and the fillet re-opens it to the full turn.
    let edge = explore_unique(&model, &joined.shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            model
                .node(e)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.curve3d())
                .and_then(|r| {
                    let ogeom_topo::EdgeRepr::Curve3d { curve, .. } = r else {
                        return None;
                    };
                    model.geometry().curve(*curve)
                })
                .is_some_and(|c| {
                    use ogeom_geom::Curve3d as _;
                    matches!(c, ogeom_geom::Curve::BSpline(_)) && {
                        let (lo, hi) = c.domain();
                        c.point_at(lo, T)
                            .and_then(|p| c.point_at(hi, T).map(|q| p.distance(q)))
                            .is_ok_and(|d| d <= 1e-6)
                    }
                })
        })
        .expect("the joined part has its fitted seam");

    let radius = 1.5;
    let result = ogeom_fillet::fillet_edge(&mut model, &joined.shape, &edge, radius, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!(result.history.is_deleted(&edge));

    // The seam runs roughly the branch's circumference, stretched by the
    // drum's fall-away; twice that is a safe ceiling, and the wedge's
    // cross-section never exceeds the radius square.
    let over_bound = radius * radius * (2.0 * core::f64::consts::PI * 5.0 * 2.0);
    let mut previous_error = f64::INFINITY;
    let mut held = None;
    for chord in [1e-2, 1e-3] {
        // Added material is a difference of two meshed volumes; at matched
        // chords the tessellation deficits cancel instead of swamping it.
        let measured = volume(&model, &result.shape, chord);
        let added = measured - volume(&model, &joined.shape, chord);
        assert!(
            added > 0.0 && added < over_bound,
            "the fillet added {added}, outside (0, {over_bound})"
        );
        if let Some(prior) = held {
            let drift: f64 = measured - prior;
            assert!(
                drift.abs() < previous_error.abs().max(1e-6),
                "refinement should settle: {prior} then {measured}"
            );
            previous_error = drift;
        }
        held = Some(measured);
    }

    let has_spline_face = explore(&model, &result.shape, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .any(|f| {
            model
                .node(&f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::BSpline(_)))
        });
    assert!(has_spline_face, "the blend face rides its fitted surface");
}

#[test]
fn two_equal_drums_refuse_the_tangent_pole_by_name() {
    // Equal radii pinch: the seam's ellipses pass through the two points
    // where the drums are tangent, the ball's section collapses there, and
    // the honest answer is a refusal that says so.
    let mut model = ogeom_topo::Model::new();
    let upright = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 20.0, T).unwrap();
    let across_frame = Frame::new(
        Point::new(-20.0, 0.0, 10.0),
        ogeom_math::Direction::X,
        ogeom_math::Direction::Y,
        T,
    )
    .unwrap();
    let across = ogeom_algo::make_cylinder(&mut model, across_frame, 5.0, 40.0, T).unwrap();
    let joined = ogeom_bool::fuse(&mut model, &upright.shape, &across.shape, T).unwrap();
    let edge = explore_unique(&model, &joined.shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            model
                .node(e)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.curve3d())
                .and_then(|r| {
                    let ogeom_topo::EdgeRepr::Curve3d { curve, .. } = r else {
                        return None;
                    };
                    model.geometry().curve(*curve)
                })
                .is_some_and(|c| matches!(c, ogeom_geom::Curve::Ellipse(_)))
        })
        .expect("the crossing has its elliptical seam arcs");
    let err = ogeom_fillet::fillet_edge(&mut model, &joined.shape, &edge, 1.0, T).unwrap_err();
    assert!(err.to_string().contains("tangent"), "{err}");
}

/// The edge of a converted solid nearest a point, straight or curved as asked.
fn converted_edge_near(
    model: &ogeom_topo::Model,
    shape: &ogeom_topo::Shape,
    near: Point,
    curved: bool,
) -> ogeom_topo::Shape {
    use ogeom_geom::Curve3d as _;
    let mut best = None;
    let mut held = f64::INFINITY;
    for e in explore(model, shape, Filter::OfType(ShapeType::Edge)).unwrap() {
        let data = model.node(&e).unwrap().data().as_edge().unwrap().clone();
        let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let c = model.geometry().curve(*curve).unwrap();
        let a = c.point_at(range.0, T).unwrap();
        let b = c.point_at(range.1, T).unwrap();
        let mid = c.point_at(f64::midpoint(range.0, range.1), T).unwrap();
        let sag = (mid - a).cross(b - a).magnitude() / (b - a).magnitude().max(1e-12);
        if curved != (sag > 1e-6) {
            continue;
        }
        let d = mid.distance(near);
        if d < held {
            held = d;
            best = Some(e.clone());
        }
    }
    best.expect("an edge of that kind near the point")
}

/// A fillet on a converted solid, the hosts patches all: round the edge,
/// check the result, and return the volumes before and after.
fn fillet_converted(
    model: &mut ogeom_topo::Model,
    solid: &ogeom_topo::Shape,
    near: Point,
    curved: bool,
    radius: f64,
) -> (f64, f64) {
    let converted = ogeom_algo::to_nurbs(model, solid, T).unwrap().shape;
    for face in explore(model, &converted, Filter::OfType(ShapeType::Face)).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        assert!(matches!(
            model.geometry().surface(data.surface),
            Some(ogeom_geom::SurfaceGeometry::BSpline(_))
        ));
    }
    let edge = converted_edge_near(model, &converted, near, curved);
    let before = volume(model, &converted, 1e-3);
    let result = ogeom_fillet::fillet_edge(model, &converted, &edge, radius, T).unwrap();
    let diagnosis = ogeom_algo::check(model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!(result.history.is_deleted(&edge));
    (before, volume(model, &result.shape, 1e-3))
}

/// A straight edge between two patches melts.
///
/// A box converted to six patches, one top edge rounded: the legs stand on
/// the hosts' own patches and the melt resolves them as one chart, and the
/// wedge's caps cut the patches along straight sections a fit used to send
/// wandering. The volume drops by the corner's own fill.
#[test]
fn a_fillet_on_a_converted_box_edge_melts() {
    let mut model = ogeom_topo::Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 12.0, 8.0), T)
        .unwrap()
        .shape;
    let radius = 1.0;
    let (before, after) = fillet_converted(
        &mut model,
        &solid,
        Point::new(10.0, 0.0, 8.0),
        false,
        radius,
    );
    let removed = (1.0 - core::f64::consts::FRAC_PI_4) * radius * radius * 20.0;
    assert!((before - 1920.0).abs() < 1e-3);
    assert!(
        (before - after - removed).abs() < removed * 0.15,
        "removed {} against {removed}",
        before - after
    );
}

/// A crease round a converted post, split at its seam, melts as one loop.
///
/// A post fused onto a slab and converted: the crease is two half circles
/// on a closed patch, put back into one loop through the neighbours the
/// hosts share. Station zero re-solved onto the seam's column lands a hair
/// past the window, and the rail loop is slid back into it, or the wall
/// never meets its rail. The sections the wedge cuts in the wall are loops
/// cut at the seam, closed on it. The volume grows by the fillet's fill.
#[test]
fn a_fillet_round_a_converted_post_melts() {
    let mut model = ogeom_topo::Model::new();
    let slab = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 2.0), T).unwrap();
    let frame = Frame::new(
        Point::new(10.0, 10.0, 2.0),
        ogeom_math::Direction::Z,
        ogeom_math::Direction::X,
        T,
    )
    .unwrap();
    let post = ogeom_algo::make_cylinder(&mut model, frame, 3.0, 8.0, T).unwrap();
    let joined = ogeom_bool::fuse(&mut model, &slab.shape, &post.shape, T)
        .unwrap()
        .shape;
    let radius = 1.0;
    let (before, after) = fillet_converted(
        &mut model,
        &joined,
        Point::new(13.0, 10.0, 2.0),
        true,
        radius,
    );
    // Pappus: the fill's section swept round the crease at its centroid.
    let section = (1.0 - core::f64::consts::FRAC_PI_4) * radius * radius;
    let added = section * 2.0 * core::f64::consts::PI * (3.0 + radius * 0.2234);
    assert!(
        (after - before - added).abs() < added * 0.15,
        "added {} against {added}",
        after - before
    );
}

/// A rim of a converted cone melts.
///
/// The top rim of a frustum converted to patches: the wall is a closed
/// patch, the cap another, and the rail on the wall winds the period.
/// The volume drops by the rounded rim.
#[test]
fn a_fillet_on_a_converted_cone_rim_melts() {
    let mut model = ogeom_topo::Model::new();
    let solid = ogeom_algo::make_cone(&mut model, Frame::WORLD, 6.0, 3.0, 10.0, T)
        .unwrap()
        .shape;
    let (before, after) =
        fillet_converted(&mut model, &solid, Point::new(3.0, 0.0, 10.0), true, 1.0);
    assert!(
        after < before - 1.0 && after > before - 6.0,
        "{before} -> {after}"
    );
}
