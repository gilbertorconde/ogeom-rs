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

    // The seat: the elliptical edge where the leaning post meets the top:
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
    // round: the marched fillet's home ground.
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

    // The seat: an arc of the fitted seam where the branch meets the drum:
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

/// Two drums of radius 5 crossing square, fused or in common, and their
/// seam's elliptical edges.
fn equal_drums(
    model: &mut ogeom_topo::Model,
    common: bool,
) -> (ogeom_topo::Shape, Vec<ogeom_topo::Shape>) {
    let upright = ogeom_algo::make_cylinder(model, Frame::WORLD, 5.0, 20.0, T).unwrap();
    let across_frame = Frame::new(
        Point::new(-20.0, 0.0, 10.0),
        ogeom_math::Direction::X,
        ogeom_math::Direction::Y,
        T,
    )
    .unwrap();
    let across = ogeom_algo::make_cylinder(model, across_frame, 5.0, 40.0, T).unwrap();
    let joined = if common {
        ogeom_bool::common(model, &upright.shape, &across.shape, T).unwrap()
    } else {
        ogeom_bool::fuse(model, &upright.shape, &across.shape, T).unwrap()
    };
    let seam = explore_unique(model, &joined.shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
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
        .collect();
    (joined.shape, seam)
}

/// A volume from the mesh alone. The pinched band's fitted face sends the
/// integrator to the mesh, while the drums alone integrate exactly; a
/// difference of the two would be mostly the mesh's deficit. Both meshed
/// at one chord, the deficits away from the band cancel.
fn mesh_volume(model: &ogeom_topo::Model, shape: &ogeom_topo::Shape) -> f64 {
    let mesh = ogeom_mesh::triangulate(
        model,
        shape,
        ogeom_mesh::Deflection {
            chord: 1e-3,
            ..ogeom_mesh::Deflection::default()
        },
        T,
    )
    .unwrap();
    mesh.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize].to_vector());
            a.dot(b.cross(c)) / 6.0
        })
        .sum()
}

/// The largest distance, less the radius, from the band face's mesh to the
/// ball's centre line: where the drums offset to `offset` cross, the two
/// ellipses `z = 10 ± x` on the drum of that radius.
fn band_off_the_ball(
    model: &ogeom_topo::Model,
    shape: &ogeom_topo::Shape,
    offset: f64,
    radius: f64,
) -> f64 {
    let centres: Vec<Point> = (0..=20_000)
        .flat_map(|i| {
            let phi = core::f64::consts::TAU * f64::from(i) / 20_000.0;
            let (x, y) = (offset * phi.cos(), offset * phi.sin());
            [Point::new(x, y, 10.0 + x), Point::new(x, y, 10.0 - x)]
        })
        .collect();
    let mut worst = 0.0_f64;
    for face in explore(model, shape, Filter::OfType(ShapeType::Face)).unwrap() {
        let fitted = model
            .node(&face)
            .and_then(|n| n.data().as_face())
            .and_then(|d| model.geometry().surface(d.surface))
            .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::BSpline(_)));
        if !fitted {
            continue;
        }
        let mesh = ogeom_mesh::triangulate_face(
            model,
            &face,
            ogeom_mesh::Deflection {
                chord: 1e-3,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap();
        for p in &mesh.positions {
            let nearest = centres
                .iter()
                .map(|c| c.distance(*p))
                .fold(f64::INFINITY, f64::min);
            worst = worst.max((nearest - radius).abs());
        }
    }
    worst
}

#[test]
fn two_equal_drums_blend_pinched_at_their_tangent_poles() {
    // Equal radii pinch: the seam's ellipses pass through the two points
    // where the drums are tangent, and the ball's section shrinks to
    // nothing there. Every seam edge has a pole for at least one end.
    let mut model = ogeom_topo::Model::new();
    let (joined, seam) = equal_drums(&mut model, false);
    let edge = seam[0].clone();
    let radius = 1.0;
    let result = ogeom_fillet::fillet_edge(&mut model, &joined, &edge, radius, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    assert!(result.history.is_deleted(&edge));

    // The notch is concave, so the band adds material: less than the
    // quarter-round's share of the radius square along the whole seam,
    // two ellipses of semi-axes 5 and 5√2.
    let added = mesh_volume(&model, &result.shape) - mesh_volume(&model, &joined);
    let seam_length = 4.0 * core::f64::consts::PI * ((25.0 + 50.0) / 2.0_f64).sqrt();
    let bound = radius * radius * (1.0 - core::f64::consts::FRAC_PI_4) * seam_length;
    assert!(
        added > 0.0 && added < bound,
        "added {added}, outside (0, {bound})"
    );

    // The ball rolls outside both drums, its centre a radius out from each.
    // The centre line is sampled every few microns and the band is fitted
    // to two tenths of one; a hundredth covers both and the mesh.
    let off = band_off_the_ball(&model, &result.shape, 5.0 + radius, radius);
    assert!(off < 1e-2, "the band stands {off} off the rolling ball");
}

#[test]
#[ignore = "heavy"]
fn the_common_of_two_equal_drums_rounds_pinched_at_its_poles() {
    // The convex twin: the drums' common part, whose seam edges are
    // ridges. The ball rolls inside both drums and the band takes material
    // away, pinched at the same poles.
    let mut model = ogeom_topo::Model::new();
    let (joined, seam) = equal_drums(&mut model, true);
    let radius = 1.0;
    for edge in &seam {
        let result = ogeom_fillet::fillet_edge(&mut model, &joined, edge, radius, T).unwrap();
        let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        let removed = mesh_volume(&model, &joined) - mesh_volume(&model, &result.shape);
        assert!(
            removed > 0.0,
            "a ridge's fillet removes material, not {removed}"
        );
        let off = band_off_the_ball(&model, &result.shape, 5.0 - radius, radius);
        assert!(off < 1e-2, "the band stands {off} off the rolling ball");
    }
}

#[test]
#[ignore = "heavy"]
fn pinched_bands_add_up_across_the_seam() {
    // Every seam edge blended in one call adds what the edges add one at a
    // time: the bands meet only at the poles, where each pinches to a
    // point, and cap to cap where a drum's seam splits an arc. The meshes
    // at one chord agree to a few tenths of a percent on these deltas; one
    // percent is that noise with room, and far under any band's own share.
    let mut model = ogeom_topo::Model::new();
    let (joined, seam) = equal_drums(&mut model, false);
    let base = mesh_volume(&model, &joined);
    let mut each = 0.0;
    for edge in &seam {
        let one = ogeom_fillet::fillet_edge(&mut model, &joined, edge, 1.0, T).unwrap();
        each += mesh_volume(&model, &one.shape) - base;
    }
    let all = ogeom_fillet::fillet_edges(&mut model, &joined, &seam, 1.0, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &all.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let together = mesh_volume(&model, &all.shape) - base;
    assert!(
        (together - each).abs() < each.abs() * 1e-2,
        "together {together}, one at a time {each}"
    );
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
/// wedge's caps cut the patches along straight sections. The volume drops
/// by the corner's own fill.
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

/// A pipe tee: a branch of radius 0.6 along X fused into a drum of radius
/// 1 along Z, and the arcs of the junction loop the boolean leaves.
fn tee() -> (ogeom_topo::Model, ogeom_topo::Shape, Vec<ogeom_topo::Shape>) {
    use ogeom_geom::Curve3d as _;
    let mut model = ogeom_topo::Model::new();
    let at = |p: Point, z, x| Frame::new(p, z, x, T).unwrap();
    let drum = ogeom_algo::make_cylinder(
        &mut model,
        at(
            Point::new(0.0, 0.0, -3.0),
            ogeom_math::Direction::Z,
            ogeom_math::Direction::X,
        ),
        1.0,
        6.0,
        T,
    )
    .unwrap();
    let branch = ogeom_algo::make_cylinder(
        &mut model,
        at(
            Point::ORIGIN,
            ogeom_math::Direction::X,
            ogeom_math::Direction::Y,
        ),
        0.6,
        3.0,
        T,
    )
    .unwrap();
    let joined = ogeom_bool::fuse(&mut model, &drum.shape, &branch.shape, T)
        .unwrap()
        .shape;
    let arcs: Vec<_> = explore_unique(&model, &joined, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (curve, range) = model
                .node(e)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.curve3d())
                .and_then(|r| {
                    let ogeom_topo::EdgeRepr::Curve3d { curve, range, .. } = r else {
                        return None;
                    };
                    Some((model.geometry().curve(*curve)?.clone(), *range))
                })
                .unwrap();
            let p = curve.point_at(f64::midpoint(range.0, range.1), T).unwrap();
            p.x > 0.5 && p.x < 1.0 && p.z.abs() < 0.7
        })
        .collect();
    (model, joined, arcs)
}

/// The junction loop of a pipe tee rounds from any one of its arcs: the
/// blend runs the whole loop, and its leg on the branch is a band between
/// the loop and the ball's contact ring, which start on one column of the
/// branch's chart a period apart and are joined there, not round a full
/// turn.
#[test]
#[ignore = "heavy"]
fn a_tee_s_junction_loop_rounds_from_any_arc() {
    let (model, joined, arcs) = tee();
    assert_eq!(arcs.len(), 4);
    let before = volume(&model, &joined, 1e-3);
    let mut volumes = Vec::new();
    for arc in &arcs {
        let mut model = model.clone();
        let rounded = ogeom_fillet::fillet_edge(&mut model, &joined, arc, 0.1, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom_algo::check(&model, &rounded, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        volumes.push(volume(&model, &rounded, 1e-3));
    }
    // Concave all round, so it adds material: no more than a square
    // corner's worth, (1 - pi/4) r^2, round the branch's circumference.
    let ceiling = (1.0 - core::f64::consts::FRAC_PI_4) * 0.01 * core::f64::consts::TAU * 0.6;
    for v in &volumes {
        let added = v - before;
        assert!(added > 0.0 && added < ceiling, "added {added}");
        assert!((v - volumes[0]).abs() < 1e-5, "{volumes:?}");
    }
}

/// The whole junction loop of a pipe tee asked for as a chain of its arcs
/// rounds once: the first arc's blend runs the loop, and the arcs after it
/// lie on the crease it rounded, done rather than in its way.
#[test]
#[ignore = "heavy"]
fn a_tee_s_junction_loop_rounds_as_a_chain_of_its_arcs() {
    let (mut model, joined, arcs) = tee();
    let one = {
        let mut model = model.clone();
        let rounded = ogeom_fillet::fillet_edge(&mut model, &joined, &arcs[0], 0.1, T)
            .unwrap()
            .shape;
        volume(&model, &rounded, 1e-3)
    };
    let chain = ogeom_fillet::fillet_edges(&mut model, &joined, &arcs, 0.1, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &chain.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let v = volume(&model, &chain.shape, 1e-3);
    assert!((v - one).abs() < 1e-5, "{v} against {one}");
    for arc in &arcs {
        assert!(chain.history.is_deleted(arc));
    }
}

/// A rib's top crease along a wall swept from a spline, turning at the
/// spline's end into the rim of a round boss: the boss's wall leaves the
/// spline's end on a drum of radius 3.6 whose tangent there stands `turn`
/// off the spline's, towards the material for a positive turn and away
/// from it for a negative one, and runs a quarter turn round to a radial
/// end wall.
///
/// The top face is `y = 2`; the profile in `y = 0` runs from the origin
/// along a cubic spline square to the end wall `x = 0`, round the arc, in
/// along the radius, and back along `z` = that height to the end wall.
/// Returns the solid, the crease along the spline's top, the crease along
/// the arc's top, and the spline.
fn rib_into_a_boss(
    model: &mut ogeom_topo::Model,
    turn: f64,
) -> (
    ogeom_topo::Shape,
    ogeom_topo::Shape,
    ogeom_topo::Shape,
    ogeom_geom::Curve,
) {
    use ogeom_geom::{CircleCurve, Curve, Curve3d as _, LineCurve};
    use ogeom_math::{Axis, Circle, Direction, KnotVector, Plane, Vector};
    let at = |x: f64, z: f64| Point::new(x, 0.0, z);
    let spline: Curve = ogeom_geom::BSplineCurve::new(
        KnotVector::clamped_uniform(3, 4).unwrap(),
        vec![at(0.0, 0.0), at(4.0, 0.0), at(7.0, -0.6), at(10.0, 0.0)],
        T,
    )
    .unwrap()
    .into();
    let end = spline.point_at(1.0, T).unwrap();
    let leaving = spline.d1_at(1.0, T).unwrap();
    let heading = leaving.z.atan2(leaving.x) + turn;
    let radius = 3.6;
    // The material lies on the profile's left, +z at the origin; the arc
    // turns towards it.
    let centre = end + Vector::new(-heading.sin(), 0.0, heading.cos()) * radius;
    let frame = Frame::new(
        centre,
        Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
        Direction::new(end - centre, T).unwrap(),
        T,
    )
    .unwrap();
    let arc: Curve = CircleCurve::new(Circle::new(frame, radius, T).unwrap()).into();
    let sweep = core::f64::consts::FRAC_PI_2;
    let rim_end = arc.point_at(sweep, T).unwrap();
    let inner = rim_end + (centre - rim_end) * (1.0 / radius);
    let back = at(0.0, inner.z);
    let vertex = |model: &mut ogeom_topo::Model, p: Point| ogeom_algo::make_vertex(model, p).shape;
    let corners: Vec<ogeom_topo::Shape> = [at(0.0, 0.0), end, rim_end, inner, back]
        .into_iter()
        .map(|p| vertex(model, p))
        .collect();
    let line = |model: &mut ogeom_topo::Model,
                from: (&ogeom_topo::Shape, Point),
                to: (&ogeom_topo::Shape, Point)| {
        let curve: Curve =
            LineCurve::new(Axis::new(from.1, Direction::new(to.1 - from.1, T).unwrap())).into();
        ogeom_algo::make_edge_between(model, curve, (0.0, from.1.distance(to.1)), from.0, to.0, T)
            .unwrap()
            .shape
    };
    let edges = vec![
        ogeom_algo::make_edge_between(
            model,
            spline.clone(),
            (0.0, 1.0),
            &corners[0],
            &corners[1],
            T,
        )
        .unwrap()
        .shape,
        ogeom_algo::make_edge_between(model, arc, (0.0, sweep), &corners[1], &corners[2], T)
            .unwrap()
            .shape,
        line(model, (&corners[2], rim_end), (&corners[3], inner)),
        line(model, (&corners[3], inner), (&corners[4], back)),
        line(model, (&corners[4], back), (&corners[0], at(0.0, 0.0))),
    ];
    let wire = ogeom_algo::make_wire(model, &edges, T).unwrap().shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::Y),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(model, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let solid = ogeom_algo::make_prism(model, &face, Vector::new(0.0, 2.0, 0.0), T)
        .unwrap()
        .shape;
    let top_edge = |model: &ogeom_topo::Model, keep: &dyn Fn(Point, Point) -> bool| {
        explore_unique(model, &solid, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let b = ogeom_algo::shape_bounds(model, e, T).unwrap();
                let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
                (lo.y - 2.0).abs() < 1e-6 && (hi.y - 2.0).abs() < 1e-6 && keep(lo, hi)
            })
            .expect("a top crease")
    };
    let along_spline = top_edge(model, &|lo, hi| {
        lo.x.abs() < 1e-6 && (hi.x - 10.0).abs() < 1e-6 && hi.z < 0.1
    });
    let along_arc = top_edge(model, &|lo, hi| {
        lo.x > 9.0 && hi.x > 11.0 && hi.z > 1.0 && lo.z < 0.1
    });
    (solid, along_spline, along_arc, spline)
}

/// The length of a curve over its domain and its heading's turn, measured
/// in the profile plane `y = 0` with the material on the left.
fn length_and_turn(curve: &ogeom_geom::Curve) -> (f64, f64) {
    use ogeom_geom::Curve3d as _;
    let (lo, hi) = curve.domain();
    let nodes = [
        (-0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
        (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.861_136_311_594_052_6, 0.347_854_845_137_453_9),
    ];
    let pieces = 256;
    let mut length = 0.0;
    for k in 0..pieces {
        let a = (hi - lo).mul_add(f64::from(k) / f64::from(pieces), lo);
        let b = (hi - lo).mul_add(f64::from(k + 1) / f64::from(pieces), lo);
        for (x, w) in nodes {
            let t = (b - a).mul_add(x, a + b) * 0.5;
            length += w * (b - a) * 0.5 * curve.d1_at(t, T).unwrap().magnitude();
        }
    }
    let heading = |t: f64| {
        let d = curve.d1_at(t, T).unwrap();
        d.z.atan2(d.x)
    };
    (length, heading(hi) - heading(lo))
}

/// A rib's top crease along a spline wall rounded with the crease that
/// carries it on round a boss's rim, the wall turning four degrees towards
/// the material between them. The two bands meet in the plane bisecting
/// the turn, where their sections agree to the fit a marched band is held
/// to, and are joined there as a chain: the spline's band stops short of
/// the vertex by the half turn and the rim's starts there. Run out past the
/// vertex instead, the spline's band crossed the rim's crease far along it
/// and cut the rim's round short there. Each band takes its section along
/// its crease less its moment about the wall times the crease's turn
/// towards the material (the rim's a quarter turn), and the two together
/// twice that moment times tan(turn/2) less.
#[test]
fn a_round_turning_shallowly_into_a_boss_s_rim_meets_the_rim_s_round_flush() {
    let pi = core::f64::consts::PI;
    let radius = 0.85;
    let turn = 4.0_f64.to_radians();
    let section = (1.0 - pi / 4.0) * radius * radius;
    let moment = (5.0 / 6.0 - pi / 4.0) * radius.powi(3);
    for reversed in [false, true] {
        let mut model = ogeom_topo::Model::new();
        let (solid, along_spline, along_arc, spline) = rib_into_a_boss(&mut model, turn);
        let before = volume(&model, &solid, 1e-4);
        let named = if reversed {
            [along_arc, along_spline]
        } else {
            [along_spline, along_arc]
        };
        let rounded = ogeom_fillet::fillet_edges(&mut model, &solid, &named, radius, T)
            .unwrap_or_else(|e| panic!("reversed {reversed}: {e}"))
            .shape;
        let diagnosis = ogeom_algo::check(&model, &rounded, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        let removed = before - volume(&model, &rounded, 1e-4);
        let (length, bend) = length_and_turn(&spline);
        let quarter = pi / 2.0;
        let want = section.mul_add(length, -moment * bend)
            + section.mul_add(3.6 * quarter, -moment * quarter)
            - 2.0 * moment * (turn / 2.0).tan();
        assert!(
            (removed - want).abs() < want * 1e-4,
            "reversed {reversed}: {removed} against {want}"
        );
    }
}

/// A round on a crease that ends where its wall turns on a degree and a
/// half back behind it, into a boss curving away round the material,
/// stops at the crease's end. The crease's continuation past the vertex
/// runs into the material there, and the boss's own rim carries the corner
/// on. Run out along the rim instead, the band met the boss along a
/// grazing stretch the melt did not resolve. Capped square to the crease
/// at both ends (the end wall at the start stands square to it), the band
/// takes the closed form: the corner's section along the crease, less its
/// moment about the wall times the crease's turn towards the material.
#[test]
fn a_round_into_a_boss_turning_back_behind_its_wall_stops_at_the_crease_s_end() {
    let pi = core::f64::consts::PI;
    let radius = 0.85;
    let mut model = ogeom_topo::Model::new();
    let (solid, along_spline, _, spline) = rib_into_a_boss(&mut model, -1.3_f64.to_radians());
    let before = volume(&model, &solid, 1e-4);
    let rounded = ogeom_fillet::fillet_edges(&mut model, &solid, &[along_spline], radius, T)
        .unwrap()
        .shape;
    let diagnosis = ogeom_algo::check(&model, &rounded, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let removed = before - volume(&model, &rounded, 1e-4);
    let (length, bend) = length_and_turn(&spline);
    let section = (1.0 - pi / 4.0) * radius * radius;
    let moment = (5.0 / 6.0 - pi / 4.0) * radius.powi(3);
    let want = section.mul_add(length, -moment * bend);
    assert!(
        (removed - want).abs() < want * 1e-4,
        "{removed} against {want}"
    );
}
