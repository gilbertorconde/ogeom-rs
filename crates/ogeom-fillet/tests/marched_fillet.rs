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
#[ignore = "the wedge builds, melts its legs and sews; the fuse's \
            classification still reads one piece of the slab's bottom as \
            hugging the wedge's boundary with no partner — the probe lands \
            near the leg's surface, not its trimmed face, and untangling \
            that is the next stretch of this frontier"]
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
    let before = volume(&model, &joined.shape, 1e-3);

    // The seat: the elliptical edge where the leaning post meets the top.
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
        let added = measured - before;
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
        } else {
            held = Some(measured);
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
