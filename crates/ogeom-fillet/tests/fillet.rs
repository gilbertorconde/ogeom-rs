//! The rolling-ball fillet, where the ball's envelope is a cylinder.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Frame, Point};
use ogeom_topo::{Filter, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

/// The top edge of the box along y at x = 2, z = 2.
fn top_edge(model: &ogeom_topo::Model, block: &ogeom_topo::Shape) -> ogeom_topo::Shape {
    explore(model, block, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(model, e)
                .unwrap()
                .is_some_and(|(a, b)| {
                    let p = |v: &ogeom_topo::Shape| {
                        model
                            .node(v)
                            .and_then(|n| n.data().as_vertex().map(|d| d.point))
                            .unwrap()
                    };
                    let (pa, pb) = (p(&a), p(&b));
                    (pa.x - 2.0).abs() < 1e-9
                        && (pa.z - 2.0).abs() < 1e-9
                        && (pb.x - 2.0).abs() < 1e-9
                        && (pb.z - 2.0).abs() < 1e-9
                })
        })
        .expect("the box has that edge")
}

#[test]
fn a_filleted_box_edge_gains_a_tangent_cylinder() {
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let edge = top_edge(&model, &block.shape);

    let radius = 0.5;
    let result = ogeom_fillet::fillet_edge(&mut model, &block.shape, &edge, radius, T).unwrap();

    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // The removed sliver is the corner square minus the quarter disc, run
    // along the edge. The mesh volume converges to it with the chord.
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    let props = ogeom_algo::volume_properties(&model, &result.shape, fine, T).unwrap();
    let exact = 8.0 - (radius * radius - core::f64::consts::PI * radius * radius / 4.0) * 2.0;
    assert!(
        (props.mass - exact).abs() < 1e-3,
        "fillet volume {} against {exact}",
        props.mass
    );

    // Seven faces: six of the box (two now trimmed back to the tangency
    // lines) plus the blend cylinder.
    let faces = explore(&model, &result.shape, Filter::OfType(ShapeType::Face)).unwrap();
    assert_eq!(faces.len(), 7);

    // The history knows the rounded edge is gone.
    assert!(result.history.is_deleted(&edge));

    // The blend passes through the arc midpoint: the ball's centre line runs
    // at (1.5, y, 1.5), and the surface lies half a radius further out along
    // the diagonal.
    let mid = Point::new(
        1.5 + radius / core::f64::consts::SQRT_2,
        1.0,
        1.5 + radius / core::f64::consts::SQRT_2,
    );
    let on_blend = faces.iter().any(|f| {
        ogeom_algo::classify_on_face(&model, f, mid, fine, T)
            .map(|c| c == ogeom_algo::Containment::In)
            .unwrap_or(false)
    });
    assert!(on_blend, "the blend face passes through the arc midpoint");

    // Tangency is the point of a fillet: where the blend meets the top face
    // there is no crease. The two surfaces share their normal along the
    // tangency line at (1.5, y, 2.0): the cylinder's radial direction there
    // is +z, the plane's normal exactly +z.
    let tangent = Point::new(1.5, 1.0, 2.0);
    let touching = faces
        .iter()
        .filter(|f| {
            ogeom_algo::classify_on_face(&model, f, tangent, fine, T)
                .map(|c| c != ogeom_algo::Containment::Out)
                .unwrap_or(false)
        })
        .count();
    assert!(
        touching >= 2,
        "the tangency line lies on both the top face and the blend, found {touching}"
    );
}

#[test]
fn a_cylinder_rim_gains_a_toroidal_blend() {
    let mut model = ogeom_topo::Model::new();
    let (radius, height, blend) = (1.0, 2.0, 0.3);
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, radius, height, T).unwrap();

    // The top rim: the circular edge at z = height.
    let edge = explore(&model, &drum.shape, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(&model, e)
                .unwrap()
                .is_some_and(|(a, _)| {
                    model
                        .node(&a)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .is_some_and(|p| (p.z - height).abs() < 1e-9 && p.x.hypot(p.y) > 0.5)
                })
        })
        .expect("the cylinder has its top rim");

    let result = ogeom_fillet::fillet_edge(&mut model, &drum.shape, &edge, blend, T).unwrap();

    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // Pappus over the meridian cusp: the removed ring is the corner square
    // minus the quarter disc, each weighted by its distance from the axis.
    let pi = core::f64::consts::PI;
    let removed = 2.0
        * pi
        * (radius * blend * blend
            - blend * blend * blend / 2.0
            - (radius - blend) * pi * blend * blend / 4.0
            - blend * blend * blend / 3.0);
    let exact = pi * radius * radius * height - removed;
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    let props = ogeom_algo::volume_properties(&model, &result.shape, fine, T).unwrap();
    assert!(
        (props.mass - exact).abs() < 2e-3,
        "rim fillet volume {} against {exact}",
        props.mass
    );

    // Four faces: the wall and both caps trimmed back, plus the torus.
    let faces = explore(&model, &result.shape, Filter::OfType(ShapeType::Face)).unwrap();
    assert_eq!(faces.len(), 4);

    assert!(result.history.is_deleted(&edge));

    // The blend passes through the meridian midpoint of the quarter tube.
    let mid = Point::new(
        (radius - blend) + blend / core::f64::consts::SQRT_2,
        0.0,
        (height - blend) + blend / core::f64::consts::SQRT_2,
    );
    let on_blend = faces.iter().any(|f| {
        ogeom_algo::classify_on_face(&model, f, mid, fine, T)
            .map(|c| c != ogeom_algo::Containment::Out)
            .unwrap_or(false)
    });
    assert!(
        on_blend,
        "the blend face passes through the tube's meridian midpoint"
    );
}

/// The circular edge of `solid` through `(rho, 0, z)`.
fn rim_through(
    model: &ogeom_topo::Model,
    solid: &ogeom_topo::Shape,
    rho: f64,
    z: f64,
) -> ogeom_topo::Shape {
    ogeom_topo::explore_unique(model, solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(model, e)
                .unwrap()
                .is_some_and(|(a, _)| {
                    model
                        .node(&a)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .is_some_and(|p| {
                            (p.z - z).abs() < 1e-7 && (p.x.hypot(p.y) - rho).abs() < 1e-7
                        })
                })
        })
        .expect("the solid has that rim")
}

/// What a ball of radius `r` rounds off a corner revolved about the `z`
/// axis, by Pappus: in the meridian the corner stands at `c = (rho, z)`,
/// its faces leaving it along the unit `a` and `b`. The ball's centre lies
/// on the bisector `r / sin(theta / 2)` from the corner and touches each
/// face `r / tan(theta / 2)` along it; the region is the kite of the
/// corner, the contacts and the centre less the ball's sector, each
/// weighted by its centroid's distance from the axis.
fn pappus_rounded(c: (f64, f64), a: (f64, f64), b: (f64, f64), r: f64) -> f64 {
    let theta = a.0.mul_add(b.0, a.1 * b.1).acos();
    let along = r / (theta / 2.0).tan();
    let bisector = {
        let s = (a.0 + b.0, a.1 + b.1);
        let m = s.0.hypot(s.1);
        (s.0 / m, s.1 / m)
    };
    let off = r / (theta / 2.0).sin();
    let centre = (bisector.0.mul_add(off, c.0), bisector.1.mul_add(off, c.1));
    let touch_a = (a.0.mul_add(along, c.0), a.1.mul_add(along, c.1));
    let touch_b = (b.0.mul_add(along, c.0), b.1.mul_add(along, c.1));
    // Each triangle's area times its centroid's distance from the axis.
    let moment = |p: (f64, f64), q: (f64, f64), s: (f64, f64)| {
        let area = (q.0 - p.0)
            .mul_add(s.1 - p.1, -(s.0 - p.0) * (q.1 - p.1))
            .abs()
            / 2.0;
        area * (p.0 + q.0 + s.0) / 3.0
    };
    let kite = moment(c, touch_a, centre) + moment(c, centre, touch_b);
    let phi = core::f64::consts::PI - theta;
    let sector = phi * r * r / 2.0;
    let sector_rho = centre.0 - bisector.0 * (4.0 * r * (phi / 2.0).sin() / (3.0 * phi));
    2.0 * core::f64::consts::PI * sector.mul_add(-sector_rho, kite)
}

/// Rims where a flat top or floor meets a drafted wall, a cone, rounded
/// whole, drafted either way: a frustum's top and foot, a tapered bore's
/// rim through a block, and a frustum standing on a block at its foot
/// (concave, filled). The ball's centre runs on a circle about the axis,
/// so the round is a torus, and the volume taken or added is the corner's
/// region in the meridian swept round by Pappus, to the measure's own
/// precision. Fitted as a spline band, each round stood up to 5e-5 off.
#[test]
fn a_cone_rim_rounds_to_its_torus() {
    let r = 0.4;
    let h = 2.0;
    let unit = |v: (f64, f64)| {
        let m = v.0.hypot(v.1);
        (v.0 / m, v.1 / m)
    };
    let at =
        |p: Point| Frame::new(p, ogeom_math::Direction::Z, ogeom_math::Direction::X, T).unwrap();
    let measured = |model: &ogeom_topo::Model, shape: &ogeom_topo::Shape| {
        ogeom_algo::volume_properties(
            model,
            shape,
            ogeom_mesh::Deflection {
                chord: 1e-4,
                ..ogeom_mesh::Deflection::default()
            },
            T,
        )
        .unwrap()
        .mass
    };
    let one = |model: &ogeom_topo::Model, shape: &ogeom_topo::Shape| {
        ogeom_topo::explore_unique(model, shape, ShapeType::Solid).unwrap()[0].clone()
    };
    for (foot, top) in [(3.0_f64, 2.6_f64), (2.6, 3.0), (3.0, 2.0), (2.0, 3.0)] {
        let down = unit((foot - top, -h));
        let up = unit((top - foot, h));
        let mut cases = Vec::new();
        // A frustum's top and foot.
        for (c, a, b) in [
            ((top, h), (-1.0, 0.0), down),
            ((foot, 0.0), (-1.0, 0.0), up),
        ] {
            let mut model = ogeom_topo::Model::new();
            let solid = ogeom_algo::make_cone(&mut model, Frame::WORLD, foot, top, h, T)
                .unwrap()
                .shape;
            cases.push((model, solid, c, a, b, true));
        }
        // A tapered bore through a block, its top rim.
        {
            let mut model = ogeom_topo::Model::new();
            let block = ogeom_algo::make_box(
                &mut model,
                at(Point::new(-5.0, -5.0, 0.0)),
                (10.0, 10.0, h),
                T,
            )
            .unwrap()
            .shape;
            let k = (top - foot) / h;
            let hole = ogeom_algo::make_cone(
                &mut model,
                at(Point::new(0.0, 0.0, -1.0)),
                foot - k,
                top + k,
                h + 2.0,
                T,
            )
            .unwrap()
            .shape;
            let cut = ogeom_bool::cut(&mut model, &block, &hole, T).unwrap().shape;
            let solid = one(&model, &cut);
            cases.push((model, solid, (top, h), (1.0, 0.0), down, true));
        }
        // A frustum standing on a block, its foot.
        {
            let mut model = ogeom_topo::Model::new();
            let block = ogeom_algo::make_box(
                &mut model,
                at(Point::new(-5.0, -5.0, -1.0)),
                (10.0, 10.0, 1.0),
                T,
            )
            .unwrap()
            .shape;
            let boss = ogeom_algo::make_cone(&mut model, Frame::WORLD, foot, top, h, T)
                .unwrap()
                .shape;
            let fused = ogeom_bool::fuse(&mut model, &block, &boss, T)
                .unwrap()
                .shape;
            let fused = one(&model, &fused);
            let unified = ogeom_heal::unify_same_domain(&mut model, &fused, T)
                .unwrap()
                .0
                .shape;
            let solid = one(&model, &unified);
            cases.push((model, solid, (foot, 0.0), (1.0, 0.0), up, false));
        }
        for (mut model, solid, c, a, b, cuts) in cases {
            let before = measured(&model, &solid);
            let edge = rim_through(&model, &solid, c.0, c.1);
            let rounded = ogeom_fillet::fillet_edge(&mut model, &solid, &edge, r, T).unwrap();
            let diagnosis = ogeom_algo::check(&model, &rounded.shape, T).unwrap();
            assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
            let after = measured(&model, &rounded.shape);
            let changed = if cuts { before - after } else { after - before };
            let want = pappus_rounded(c, a, b, r);
            assert!(
                (changed - want).abs() < 1e-9,
                "foot {foot}, top {top}, rim at {c:?}: {changed} against {want}"
            );
            let tori = ogeom_topo::explore_unique(&model, &rounded.shape, ShapeType::Face)
                .unwrap()
                .iter()
                .filter(|f| {
                    model
                        .node(f)
                        .and_then(|n| n.data().as_face())
                        .and_then(|d| model.geometry().surface(d.surface))
                        .is_some_and(|s| matches!(s, ogeom_geom::SurfaceGeometry::Torus(_)))
                })
                .count();
            assert_eq!(
                tori, 1,
                "foot {foot}, top {top}, rim at {c:?}: the round is a torus"
            );
        }
    }
}

/// A drum standing at a slant, either rim rounded. The tube's ring on the
/// wall lies on the torus chart's outer equator, the seam of its tube
/// angle, where rounding reads it at either end of the chart; the blend is
/// the quarter of the tube between the wall and the cap whichever end it
/// read, so the volume is the drum's less one Pappus ring.
#[test]
fn a_slanted_drum_s_rims_round_the_quarter_tube() {
    use ogeom_math::{Direction, Vector};
    let pi = core::f64::consts::PI;
    let (radius, height, blend) = (1.0, 3.0, 0.2);
    let frame = Frame::new(
        Point::ORIGIN,
        Direction::new(Vector::new(0.0, 0.6, 0.8), T).unwrap(),
        Direction::new(Vector::new(0.0, 0.8, -0.6), T).unwrap(),
        T,
    )
    .unwrap();
    let k = (10.0 - 3.0 * pi) / (3.0 * (4.0 - pi));
    let removed = 2.0 * pi * (radius - k * blend) * (1.0 - pi / 4.0) * blend * blend;
    let want = pi * radius * radius * height - removed;
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    for rim in 0..2 {
        let mut model = ogeom_topo::Model::new();
        let drum = ogeom_algo::make_cylinder(&mut model, frame, radius, height, T)
            .unwrap()
            .shape;
        let rims: Vec<_> = ogeom_topo::explore_unique(&model, &drum, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| {
                ogeom_algo::edge_vertices(&model, e)
                    .unwrap()
                    .is_some_and(|(a, b)| a.is_same(&b))
            })
            .collect();
        assert_eq!(rims.len(), 2);
        let result = ogeom_fillet::fillet_edge(&mut model, &drum, &rims[rim], blend, T).unwrap();
        let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        let got = ogeom_algo::volume_properties(&model, &result.shape, fine, T)
            .unwrap()
            .mass;
        assert!((got - want).abs() < 1e-6, "rim {rim}: {got} against {want}");
    }
}

#[test]
fn a_rim_fillet_that_swallows_the_axis_is_refused() {
    let mut model = ogeom_topo::Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 1.0, 2.0, T).unwrap();
    let edge = explore(&model, &drum.shape, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(&model, e)
                .unwrap()
                .is_some_and(|(a, _)| {
                    model
                        .node(&a)
                        .and_then(|n| n.data().as_vertex().map(|d| d.point))
                        .is_some_and(|p| (p.z - 2.0).abs() < 1e-9 && p.x.hypot(p.y) > 0.5)
                })
        })
        .expect("the cylinder has its top rim");
    assert!(ogeom_fillet::fillet_edge(&mut model, &drum.shape, &edge, 1.0, T).is_err());
}

#[test]
fn a_variable_radius_fillet_widens_along_its_edge() {
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let edge = top_edge(&model, &block.shape);

    let (r0, r1) = (0.3, 0.6);
    let result =
        ogeom_fillet::fillet_edge_variable(&mut model, &block.shape, &edge, r0, r1, T).unwrap();

    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // The removed sliver integrates the section area along the edge: for the
    // linear law, L (1 - pi/4) times the integral of r(u) squared.
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    let integral = (r1 * r1 * r1 - r0 * r0 * r0) / (3.0 * (r1 - r0));
    let exact = 8.0 - 2.0 * (1.0 - core::f64::consts::PI / 4.0) * integral;
    let props = ogeom_algo::volume_properties(&model, &result.shape, fine, T).unwrap();
    assert!(
        (props.mass - exact).abs() < 2e-3,
        "variable fillet volume {} against {exact}",
        props.mass
    );

    assert_eq!(
        explore(&model, &result.shape, Filter::OfType(ShapeType::Face))
            .unwrap()
            .len(),
        7
    );
    assert!(result.history.is_deleted(&edge));

    // The blend passes through the mid-section's arc midpoint, where the
    // radius is the mean of the two ends.
    let r = f64::midpoint(r0, r1);
    let mid = Point::new(
        (2.0 - r) + r / core::f64::consts::SQRT_2,
        1.0,
        (2.0 - r) + r / core::f64::consts::SQRT_2,
    );
    let faces = explore(&model, &result.shape, Filter::OfType(ShapeType::Face)).unwrap();
    let on_blend = faces.iter().any(|f| {
        ogeom_algo::classify_on_face(&model, f, mid, fine, T)
            .map(|c| c != ogeom_algo::Containment::Out)
            .unwrap_or(false)
    });
    assert!(
        on_blend,
        "the blend passes through the mean-radius midpoint"
    );
}

#[test]
fn a_concave_edge_gains_an_additive_fillet() {
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
    let seat = Frame::new(
        Point::new(1.0, -0.5, 1.0),
        ogeom_math::Direction::Z,
        ogeom_math::Direction::X,
        T,
    )
    .unwrap();
    let notch = ogeom_algo::make_box(&mut model, seat, (2.0, 3.0, 2.0), T).unwrap();
    let cut = ogeom_bool::cut(&mut model, &block.shape, &notch.shape, T).unwrap();

    // The re-entrant edge of the L: along y at x = 1, z = 1.
    let edge = explore(&model, &cut.shape, Filter::OfType(ShapeType::Edge))
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(&model, e)
                .unwrap()
                .is_some_and(|(a, b)| {
                    let p = |v: &ogeom_topo::Shape| {
                        model
                            .node(v)
                            .and_then(|n| n.data().as_vertex().map(|d| d.point))
                            .unwrap()
                    };
                    let (pa, pb) = (p(&a), p(&b));
                    (pa.x - 1.0).abs() < 1e-9
                        && (pa.z - 1.0).abs() < 1e-9
                        && (pb.x - 1.0).abs() < 1e-9
                        && (pb.z - 1.0).abs() < 1e-9
                })
        })
        .expect("the L has its re-entrant edge");

    let radius = 0.25;
    let result = ogeom_fillet::fillet_edge(&mut model, &cut.shape, &edge, radius, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);

    // The L keeps its volume plus the added sliver along the re-entrant
    // edge: the corner square minus the quarter disc, run along its length.
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    let added = radius * radius * (1.0 - core::f64::consts::PI / 4.0) * 2.0;
    let exact = 6.0 + added;
    let props = ogeom_algo::volume_properties(&model, &result.shape, fine, T).unwrap();
    assert!(
        (props.mass - exact).abs() < 1e-3,
        "concave fillet volume {} against {exact}",
        props.mass
    );
    assert!(result.history.is_deleted(&edge));
}

/// Two blends on one body: a chamfer and then a fillet on far-apart edges.
/// The second blend's tools march against the first blend's faces, find no
/// crossing, and that is the answer, apart, rather than a refusal.
#[test]
fn a_second_blend_lands_beside_the_first() {
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 12.0), T)
        .unwrap()
        .shape;
    let bevel_edge = edge_with_midpoint(&model, &block, Point::new(20.0, 30.0, 12.0));
    let bevelled = ogeom_fillet::chamfer_edge(&mut model, &block, &bevel_edge, 1.5, T)
        .unwrap()
        .shape;
    let round_edge = edge_with_midpoint(&model, &bevelled, Point::new(20.0, 0.0, 12.0));
    let blended = ogeom_fillet::fillet_edge(&mut model, &bevelled, &round_edge, 2.0, T)
        .unwrap()
        .shape;

    let diagnosis = ogeom_algo::check(&model, &blended, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    // The volume lost is the chamfer wedge plus the fillet's corner ring:
    // (1.5^2 / 2) * 40 and (4 - pi) * 40 for r = 2.
    let fine = ogeom_mesh::Deflection {
        chord: 1e-2,
        ..ogeom_mesh::Deflection::default()
    };
    let props = ogeom_algo::volume_properties(&model, &blended, fine, T).unwrap();
    let expected =
        40.0 * 30.0 * 12.0 - 1.5 * 1.5 / 2.0 * 40.0 - (4.0 - core::f64::consts::PI) * 40.0;
    assert!(
        (props.mass - expected).abs() < expected * 1e-3,
        "volume {} vs expected {expected}",
        props.mass
    );
}

/// The edge whose curve midpoint lies nearest `at`.
fn edge_with_midpoint(
    model: &ogeom_topo::Model,
    solid: &ogeom_topo::Shape,
    at: Point,
) -> ogeom_topo::Shape {
    use ogeom_geom::Curve3d as _;
    let mut best: Option<(f64, ogeom_topo::Shape)> = None;
    for edge in ogeom_topo::explore_unique(model, solid, ShapeType::Edge).unwrap() {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(ogeom_topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            continue;
        };
        let p = geometry
            .point_at(f64::midpoint(range.0, range.1), T)
            .unwrap();
        let d = p.distance(at);
        if best.as_ref().is_none_or(|(held, _)| d < *held) {
            best = Some((d, edge));
        }
    }
    best.unwrap().1
}

/// A ball too big for the faces it would sit on is refused, however big:
/// its band would run past the far side of a face.
#[test]
fn a_fillet_the_faces_cannot_hold_is_refused() {
    for radius in [2.5, 1e6] {
        let mut model = ogeom_topo::Model::new();
        let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 2.0), T).unwrap();
        let edge = top_edge(&model, &block.shape);
        assert!(
            ogeom_fillet::fillet_edge(&mut model, &block.shape, &edge, radius, T).is_err(),
            "radius {radius}"
        );
    }
}

/// Two edges of a narrow face whose balls would overlap across it are
/// refused together, where each alone fits.
#[test]
fn blends_overlapping_across_a_shared_face_are_refused() {
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (2.0, 2.0, 0.5), T).unwrap();
    let along_y_at = |model: &ogeom_topo::Model, z: f64| {
        explore(model, &block.shape, Filter::OfType(ShapeType::Edge))
            .unwrap()
            .into_iter()
            .find(|e| {
                ogeom_algo::edge_vertices(model, e)
                    .unwrap()
                    .is_some_and(|(a, b)| {
                        let p = |v: &ogeom_topo::Shape| {
                            model
                                .node(v)
                                .and_then(|n| n.data().as_vertex().map(|d| d.point))
                                .unwrap()
                        };
                        let (pa, pb) = (p(&a), p(&b));
                        (pa.x - 2.0).abs() < 1e-9
                            && (pb.x - 2.0).abs() < 1e-9
                            && (pa.z - z).abs() < 1e-9
                            && (pb.z - z).abs() < 1e-9
                    })
            })
            .expect("the box has that edge")
    };
    let (low, high) = (along_y_at(&model, 0.0), along_y_at(&model, 0.5));
    let fits = ogeom_fillet::fillet_edges(
        &mut model,
        &block.shape,
        &[low.clone(), high.clone()],
        0.2,
        T,
    );
    assert!(fits.is_ok());
    let overlaps = ogeom_fillet::fillet_edges(&mut model, &block.shape, &[low, high], 0.3, T);
    assert!(overlaps.is_err());
}

/// A prism standing on its profile in the plane x = 0, extruded `depth`
/// back along -x, with the profile's points given as (y, z): its face on
/// x = 0 is a wall, and the profile's first side, along y on z = 0, is a
/// floor meeting that wall.
fn profile_prism(
    model: &mut ogeom_topo::Model,
    profile: &[(f64, f64)],
    depth: f64,
) -> ogeom_topo::Shape {
    use ogeom_math::{Direction, Plane, Vector};
    let points: Vec<Point> = profile
        .iter()
        .map(|&(y, z)| Point::new(0.0, y, z))
        .collect();
    let wire = ogeom_algo::make_polygon(model, &points, true, T)
        .unwrap()
        .shape;
    let plane = ogeom_geom::PlaneSurface::over(
        Plane::through(Point::ORIGIN, Direction::X),
        (-100.0, 100.0),
        (-100.0, 100.0),
    )
    .unwrap();
    let face = ogeom_algo::make_face(model, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    ogeom_algo::make_prism(model, &face, Vector::new(-depth, 0.0, 0.0), T)
        .unwrap()
        .shape
}

/// The volume a fillet of `radius` on the edge of a profile prism through
/// `middle` took from it, rounded alone by `fillet_edge` and as a request
/// of one edge by `fillet_edges`: the same edge on the same solid, so the
/// same rounding.
fn removed_by_round(
    profile: &[(f64, f64)],
    middle: Point,
    radius: f64,
) -> [(&'static str, f64); 2] {
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    let round = |alone: bool| {
        let mut model = ogeom_topo::Model::new();
        let prism = profile_prism(&mut model, profile, 2.0);
        let before = ogeom_algo::volume_properties(&model, &prism, fine, T)
            .unwrap()
            .mass;
        let edge = edge_with_midpoint(&model, &prism, middle);
        let rounded = if alone {
            ogeom_fillet::fillet_edge(&mut model, &prism, &edge, radius, T)
        } else {
            ogeom_fillet::fillet_edges(&mut model, &prism, &[edge], radius, T)
        }
        .unwrap()
        .shape;
        let diagnosis = ogeom_algo::check(&model, &rounded, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        before
            - ogeom_algo::volume_properties(&model, &rounded, fine, T)
                .unwrap()
                .mass
    };
    [("fillet_edge", round(true)), ("fillet_edges", round(false))]
}

/// A round on a floor edge whose floor gives way to a ramp at the edge's
/// end, the wall carrying on above it: the round's cylinder runs on past
/// the end along the wall until the ramp has cut its whole section away,
/// and stops on the ramp with no cap. Past the end, the section above a
/// ramp of slope k standing h over the floor is the corner outside the
/// quarter circle above h, and summed along the ramp it is
/// r^3 (5/6 - pi/4) / k, on top of the section (1 - pi/4) r^2 along the
/// edge.
#[test]
fn a_round_runs_out_along_its_wall_onto_a_ramp() {
    let (length, height, radius) = (3.0, 2.0, 1.0);
    for slope in [3.0_f64.sqrt(), 1.0, 1.0 / 3.0_f64.sqrt()] {
        let profile = [
            (0.0, 0.0),
            (length, 0.0),
            (length + height / slope, height),
            (0.0, height),
        ];
        let want = (1.0 - core::f64::consts::FRAC_PI_4) * radius * radius * length
            + radius.powi(3) * (5.0 / 6.0 - core::f64::consts::FRAC_PI_4) / slope;
        for (by, removed) in removed_by_round(&profile, Point::new(0.0, length / 2.0, 0.0), radius)
        {
            assert!(
                (removed - want).abs() < 1e-4,
                "slope {slope}: the round by {by} took {removed} against {want}"
            );
        }
    }
}

/// Where the floor steps down or up at the edge's end and the material
/// carries on past it, no ball on the edge reaches that material: the
/// round stops square at the end, taking its section along the edge and
/// nothing else, though its section leaves the material at the far end
/// of the step.
#[test]
fn a_round_stops_square_at_a_step_past_its_end() {
    let (length, height, radius) = (3.0, 2.0, 1.0);
    let down = [
        (0.0, 0.0),
        (length, 0.0),
        (length, -0.5),
        (length + 2.0, -0.5),
        (length + 2.0, height),
        (0.0, height),
    ];
    let up = [
        (0.0, 0.0),
        (length, 0.0),
        (length, 0.5),
        (length + 2.0, 0.5),
        (length + 2.0, height),
        (0.0, height),
    ];
    for (name, profile) in [("down", &down), ("up", &up)] {
        let want = (1.0 - core::f64::consts::FRAC_PI_4) * radius * radius * length;
        for (by, removed) in removed_by_round(profile, Point::new(0.0, length / 2.0, 0.0), radius) {
            assert!(
                (removed - want).abs() < 1e-4,
                "step {name}: the round by {by} took {removed} against {want}"
            );
        }
    }
}

/// A round on a rib's top edge where the rib runs up under a slanted roof,
/// the roof turned over the rib's corner: the band stops on the roof,
/// cut off by its plane, and takes nothing from what the roof bounds. At
/// a point of the section a distance b in from the rib's top, the roof
/// stands h - k b over the floor, k its slope, so the band takes the
/// section (1 - pi/4) r^2 times h less k times the section's centroid's
/// distance from the top, r (10 - 3 pi) / (3 (4 - pi)).
#[test]
fn a_round_running_up_under_a_roof_stops_on_it() {
    let (height, radius) = (2.0, 0.5);
    let section = (1.0 - core::f64::consts::FRAC_PI_4) * radius * radius;
    let centroid =
        radius * (10.0 - 3.0 * core::f64::consts::PI) / (3.0 * (4.0 - core::f64::consts::PI));
    for slope in [0.5, 1.0, 2.0] {
        let profile = [
            (-3.0, 0.0),
            (0.0, 0.0),
            (0.0, height),
            (2.0, 2.0f64.mul_add(slope, height)),
            (2.0, 2.0f64.mul_add(slope, height + 1.0)),
            (-3.0, 2.0f64.mul_add(slope, height + 1.0)),
        ];
        let want = section * slope.mul_add(-centroid, height);
        for (by, removed) in removed_by_round(&profile, Point::new(0.0, 0.0, height / 2.0), radius)
        {
            assert!(
                (removed - want).abs() < 1e-4,
                "slope {slope}: the round by {by} took {removed} against {want}"
            );
        }
    }
}

/// A round on the floor edge of a prism whose wall is a triangle, its apex
/// low over the edge: the wall narrows to a tip at both ends of the edge,
/// and the ball's contact runs off it there, out of the solid. The band
/// is cut off on the faces the edge runs into at its ends, as a band
/// running up under a roof is. Where the wall stands h <= r over the edge,
/// the band's section is the part of the corner outside the quarter
/// circle below h, of area r h - (F(r) - F(r - h)) with
/// F(u) = (u sqrt(r^2 - u^2) + r^2 asin(u / r)) / 2. Over a tip rising at
/// slope k that integrates to r^3 / (6 k), and over the rest of the edge
/// the section is the whole (1 - pi/4) r^2. The ball fits at some stations
/// and runs off at others, alone and as a request of one edge alike.
#[test]
fn a_round_on_a_narrowing_wall_is_cut_off_where_the_wall_ends() {
    let (length, apex, height) = (10.0, 1.0, 1.0);
    let rises = [height / apex, height / (length - apex)];
    for radius in [0.3, 0.5, 0.7] {
        let profile = [(0.0, 0.0), (length, 0.0), (apex, height)];
        let want = (1.0 - core::f64::consts::FRAC_PI_4)
            * radius
            * radius
            * (length - rises.iter().map(|k| radius / k).sum::<f64>())
            + rises
                .iter()
                .map(|k| radius.powi(3) / (6.0 * k))
                .sum::<f64>();
        for (by, removed) in removed_by_round(&profile, Point::new(0.0, length / 2.0, 0.0), radius)
        {
            assert!(
                (removed - want).abs() < 1e-6,
                "radius {radius}: the round by {by} took {removed} against {want}"
            );
        }
    }
}

/// The narrowing wall's floor edge asked together with the prism's apex
/// edge, which shares no face and no vertex with it: no band is checked
/// against the other, so the floor band runs off the wall's tips as it
/// does alone. The apex edge is straight between planes meeting at an
/// interior angle `a`, square to the prism's ends, and its band takes the
/// section `r^2 (cot(a / 2) - (pi - a) / 2)` along the prism's depth, on
/// top of the floor band's. Each band is a cylinder of the radius whose
/// axis stands the radius off both of its edge's faces. Before, a request
/// of more than one edge refused any station whose contact ran off.
#[test]
fn a_round_running_off_a_narrowing_wall_is_taken_beside_an_edge_apart() {
    use core::f64::consts::{FRAC_PI_4, PI};
    use ogeom_math::Vector;
    let (length, apex, height, depth) = (10.0_f64, 1.0_f64, 1.0_f64, 2.0_f64);
    let rises = [height / apex, height / (length - apex)];
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    // The apex's interior angle, between its walls toward the two feet.
    let to_feet = [
        Vector::new(0.0, -apex, -height),
        Vector::new(0.0, length - apex, -height),
    ];
    let interior =
        (to_feet[0].dot(to_feet[1]) / (to_feet[0].magnitude() * to_feet[1].magnitude())).acos();
    // Each face's plane through the solid: a point on it and its normal.
    let planes = |which: &str| -> [(Point, Vector); 2] {
        let wall = |foot: Vector| {
            let n = Vector::X.cross(foot);
            (Point::new(0.0, apex, height), n * (1.0 / n.magnitude()))
        };
        match which {
            "floor" => [(Point::ORIGIN, Vector::Z), (Point::ORIGIN, Vector::X)],
            _ => [wall(to_feet[0]), wall(to_feet[1])],
        }
    };
    for radius in [0.3, 0.5] {
        let mut model = ogeom_topo::Model::new();
        let profile = [(0.0, 0.0), (length, 0.0), (apex, height)];
        let prism = profile_prism(&mut model, &profile, depth);
        let before = ogeom_algo::volume_properties(&model, &prism, fine, T)
            .unwrap()
            .mass;
        let floor = edge_with_midpoint(&model, &prism, Point::new(0.0, length / 2.0, 0.0));
        let top = edge_with_midpoint(&model, &prism, Point::new(-depth / 2.0, apex, height));
        let rounded = ogeom_fillet::fillet_edges(
            &mut model,
            &prism,
            &[floor.clone(), top.clone()],
            radius,
            T,
        )
        .unwrap();
        let diagnosis = ogeom_algo::check(&model, &rounded.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
        let removed = before
            - ogeom_algo::volume_properties(&model, &rounded.shape, fine, T)
                .unwrap()
                .mass;
        let want = (1.0 - FRAC_PI_4)
            * radius
            * radius
            * (length - rises.iter().map(|k| radius / k).sum::<f64>())
            + rises
                .iter()
                .map(|k| radius.powi(3) / (6.0 * k))
                .sum::<f64>()
            + radius * radius * ((1.0 / (interior / 2.0).tan()) - (PI - interior) / 2.0) * depth;
        assert!(
            (removed - want).abs() < 1e-6,
            "radius {radius}: the two rounds took {removed} against {want}"
        );
        for (which, edge, along) in [("floor", &floor, Vector::Y), ("apex", &top, Vector::X)] {
            let mut cylinders = 0;
            for face in rounded.history.generated(edge) {
                let data = model.node(face).unwrap().data().as_face().unwrap();
                let ogeom_geom::SurfaceGeometry::Cylinder(c) =
                    model.geometry().surface(data.surface).unwrap()
                else {
                    continue;
                };
                let placed = face.transform(model.datums()).unwrap();
                let axis = c.cylinder().axis();
                let origin = placed.apply(axis.location);
                let direction = placed.apply_vector(axis.direction.vector());
                assert!(
                    (c.cylinder().radius() - radius).abs() < 1e-9
                        && direction.cross(along).magnitude() < 1e-9,
                    "{which}: a band of radius {} along {direction:?}",
                    c.cylinder().radius()
                );
                for (on, normal) in planes(which) {
                    let off = (origin - on).dot(normal).abs();
                    assert!(
                        (off - radius).abs() < 1e-9,
                        "{which}: the ball's centre stands {off} off a face, not {radius}"
                    );
                }
                cylinders += 1;
            }
            assert!(cylinders > 0, "{which}: no band was generated");
        }
    }
}

/// A disc of radius `big` at the top about the origin, two high, its drum
/// drafted by `k` over the height so its foot stands out (a cone; a
/// cylinder at zero), up to its chord `y = c` at the top, that wall drafted
/// by `d` over the height so its foot stands out (overhanging for a
/// negative `d`): the block, and its chord's top crease.
fn d_shaped_block(
    model: &mut ogeom_topo::Model,
    (big, k): (f64, f64),
    c: f64,
    d: f64,
) -> (ogeom_topo::Shape, ogeom_topo::Shape) {
    use ogeom_math::{Direction, Vector};
    let disc = if k == 0.0 {
        ogeom_algo::make_cylinder(model, Frame::WORLD, big, 2.0, T)
    } else {
        ogeom_algo::make_cone(model, Frame::WORLD, 2.0_f64.mul_add(k, big), big, 2.0, T)
    }
    .unwrap()
    .shape;
    // The box past the chord's wall, forty a side.
    let normal = Vector::new(0.0, 1.0, d) * (1.0 / d.hypot(1.0));
    let across = normal.cross(Vector::X);
    let frame = Frame::new(
        Point::new(0.0, c, 2.0) - Vector::X * 20.0 - across * 20.0,
        Direction::new(normal, T).unwrap(),
        Direction::X,
        T,
    )
    .unwrap();
    let past = ogeom_algo::make_box(model, frame, (40.0, 40.0, 40.0), T)
        .unwrap()
        .shape;
    let cut = ogeom_bool::cut(model, &disc, &past, T).unwrap().shape;
    let solid = explore(model, &cut, Filter::OfType(ShapeType::Solid)).unwrap()[0].clone();
    let edge = edge_with_midpoint(model, &solid, Point::new(0.0, c, 2.0));
    (solid, edge)
}

/// A round on the chord's top crease of a D-shaped block, a disc up to a
/// chord, the chord's wall square, drafted, or overhanging, the drum's wall
/// upright or drafted either way (a cone): at each end the drum cuts the
/// crease's corner away obliquely, and the band runs on along the chord's
/// line until its whole section has left the drum, the cut trimming it
/// there. A point of the section a distance `w` across from the drum's
/// axis, at a height where the drum's radius is `R`, then takes the drum's
/// chord through it, `2 sqrt(R^2 - w^2)`, summed over the section a column
/// at a time. Capped square at the vertices instead, the band left what of
/// its section the drum still holds past them standing, and on a drafted
/// chord wall, whose ball touches it outside the drum there, the cut
/// refused. Beside a cone the band was capped so; with the cone and the
/// chord wall drafted alike, their seam a parabola, the contact on the wall
/// crossing the cone on that seam read as off the cone, and the band was
/// capped again.
#[test]
fn a_round_on_a_chord_runs_out_across_the_drum_at_both_ends() {
    let (big, c, radius) = (2.5_f64, 0.5, 0.5);
    let fine = ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    };
    // The drum's chord through a section point, integrated over heights
    // `z0` to `z1` at an offset `w` from the axis: the drum's radius is
    // `R(z) = big + k (2 - z)`, and `f` is an antiderivative of
    // `sqrt(R^2 - w^2)` in `R`.
    let chords = |k: f64, w: f64, (z0, z1): (f64, f64)| -> f64 {
        if z1 <= z0 {
            return 0.0;
        }
        if k == 0.0 {
            return 2.0 * big.mul_add(big, -w * w).sqrt() * (z1 - z0);
        }
        let f = |r: f64| {
            let s = r.mul_add(r, -w * w).max(0.0).sqrt();
            r.mul_add(s, -w * w * (r + s).ln()) / 2.0
        };
        let at = |z: f64| k.mul_add(2.0 - z, big);
        2.0 * (f(at(z0)) - f(at(z1))) / k
    };
    for (k, d) in [
        (0.0_f64, 0.1_f64),
        (0.0, 0.0),
        (0.0, -0.1),
        (0.1, 0.1),
        (0.1, 0.0),
        (0.1, -0.1),
        (-0.1, 0.1),
        (-0.1, 0.0),
        (-0.1, -0.1),
    ] {
        // The section in the plane square to the chord, as (y, z): the
        // corner the top and the wall make, outside the ball.
        let top = (0.0, 1.0);
        let wall = (1.0 / d.hypot(1.0), d / d.hypot(1.0));
        let setback = radius / (1.0 + top.0 * wall.0 + top.1 * wall.1);
        let ball = (
            c - (top.0 + wall.0) * setback,
            2.0 - (top.1 + wall.1) * setback,
        );
        let corner = [
            (c, 2.0),
            (ball.0 + top.0 * radius, ball.1 + top.1 * radius),
            ball,
            (ball.0 + wall.0 * radius, ball.1 + wall.1 * radius),
        ];
        // Each column of the section across the chord: the corner's
        // quadrilateral there, less the ball.
        let (y0, y1) = corner
            .iter()
            .fold((f64::MAX, f64::MIN), |(a, b), p| (a.min(p.0), b.max(p.0)));
        let n = 20_000;
        let dy = (y1 - y0) / f64::from(n);
        let mut want = 0.0;
        for i in 0..n {
            let y = (f64::from(i) + 0.5).mul_add(dy, y0);
            let (mut lo, mut hi) = (f64::MAX, f64::MIN);
            for j in 0..4 {
                let (a, b) = (corner[j], corner[(j + 1) % 4]);
                if (a.0 - y) * (b.0 - y) <= 0.0 && (b.0 - a.0).abs() > f64::EPSILON {
                    let z = (b.1 - a.1).mul_add((y - a.0) / (b.0 - a.0), a.1);
                    (lo, hi) = (lo.min(z), hi.max(z));
                }
            }
            let reach = radius.mul_add(radius, -(y - ball.0).powi(2));
            want += if reach > 0.0 {
                let h = reach.sqrt();
                chords(k, y, (lo, hi.min(ball.1 - h))) + chords(k, y, (lo.max(ball.1 + h), hi))
            } else {
                chords(k, y, (lo, hi))
            };
        }
        want *= dy;
        for alone in [true, false] {
            let mut model = ogeom_topo::Model::new();
            let (solid, edge) = d_shaped_block(&mut model, (big, k), c, d);
            let before = ogeom_algo::volume_properties(&model, &solid, fine, T)
                .unwrap()
                .mass;
            let rounded = if alone {
                ogeom_fillet::fillet_edge(&mut model, &solid, &edge, radius, T)
            } else {
                ogeom_fillet::fillet_edges(&mut model, &solid, &[edge], radius, T)
            }
            .unwrap()
            .shape;
            let diagnosis = ogeom_algo::check(&model, &rounded, T).unwrap();
            assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
            let removed = before
                - ogeom_algo::volume_properties(&model, &rounded, fine, T)
                    .unwrap()
                    .mass;
            assert!(
                (removed - want).abs() < 1e-6,
                "drum draft {k}, chord draft {d}, alone {alone}: {removed} against {want}"
            );
        }
    }
}

/// A floor edge whose wall carries a boss along part of it, the boss's
/// underside 0.2 over the floor: there the wall reaches no higher than
/// that, and a ball of radius 0.5 touching the floor and the wall's plane
/// touches it 0.3 inside the material, under the boss. No ball rolling on
/// the floor and the wall rounds the edge there, and a band cut along it
/// would groove the block under the boss. Alone and as a request of one
/// edge, the round is refused by name.
#[test]
fn a_round_whose_ball_would_touch_the_wall_inside_a_boss_is_refused() {
    use ogeom_math::Direction;
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(
        &mut model,
        Frame::new(Point::new(-4.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap(),
        (4.0, 10.0, 4.0),
        T,
    )
    .unwrap()
    .shape;
    let boss = ogeom_algo::make_box(
        &mut model,
        Frame::new(Point::new(0.0, 6.0, 0.2), Direction::Z, Direction::X, T).unwrap(),
        (2.0, 3.0, 2.8),
        T,
    )
    .unwrap()
    .shape;
    let part = ogeom_bool::fuse(&mut model, &block, &boss, T)
        .unwrap()
        .shape;
    let edge = edge_with_midpoint(&model, &part, Point::new(0.0, 5.0, 0.0));
    let alone = ogeom_fillet::fillet_edge(&mut model, &part, &edge, 0.5, T);
    let together = ogeom_fillet::fillet_edges(&mut model, &part, &[edge], 0.5, T);
    for (by, outcome) in [("fillet_edge", alone), ("fillet_edges", together)] {
        let Err(ogeom_core::OgeomError::Construction(why)) = outcome else {
            panic!("{by} rounded an edge whose ball touches the wall inside the boss");
        };
        assert!(
            why.to_string()
                .contains("across a face that does not reach so far"),
            "{by}: {why}"
        );
    }
}

/// A 10 mm cube with all twelve edges rounded at radius 1 is the 8 mm cube
/// swept by the ball: twelve quarter cylinders, eight sphere octants and the
/// six trimmed faces. Each blend's sections meet its neighbours' arc edges
/// tangentially, so every one of those meetings is a line touching a circle.
#[test]
fn a_cube_rounded_on_every_edge_is_the_cube_swept_by_the_ball() {
    let mut model = ogeom_topo::Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = ogeom_topo::explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let result = ogeom_fillet::fillet_edges(&mut model, &block, &edges, 1.0, T).unwrap();

    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let faces = ogeom_topo::explore_unique(&model, &result.shape, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 6 + 12 + 8);

    let (side, r, pi) = (8.0, 1.0, core::f64::consts::PI);
    let exact = side * side * side
        + 6.0 * side * side * r
        + 3.0 * pi * side * r * r
        + 4.0 / 3.0 * pi * r * r * r;
    let props =
        ogeom_algo::volume_properties(&model, &result.shape, ogeom_mesh::Deflection::default(), T)
            .unwrap();
    assert!(
        (props.mass - exact).abs() < 1e-6,
        "rounded cube volume {} against {exact}",
        props.mass
    );
}

/// A 10 mm cube with its vertical edges rounded at radius 2, then its top
/// loop at radius 1: the second round's cross-section swept along the
/// rounded square, straight along the four sides and turned about each
/// corner's axis (Pappus). Each corner's toroidal blend rests against the
/// side faces' planes at one point.
#[test]
fn the_top_loop_of_a_block_with_rounded_verticals_rounds_on_the_turns() {
    let mut model = ogeom_topo::Model::new();
    let (side, wall, top) = (10.0, 2.0, 1.0);
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (side, side, side), T)
        .unwrap()
        .shape;
    let ends = |model: &ogeom_topo::Model, e: &ogeom_topo::Shape| {
        let (a, b) = ogeom_algo::edge_vertices(model, e).unwrap().unwrap();
        let p = |v: &ogeom_topo::Shape| model.node(v).unwrap().data().as_vertex().unwrap().point;
        (p(&a), p(&b))
    };
    let verticals: Vec<_> = ogeom_topo::explore_unique(&model, &block, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ends(&model, e);
            (a.z - b.z).abs() > side / 2.0
        })
        .collect();
    assert_eq!(verticals.len(), 4);
    let rounded = ogeom_fillet::fillet_edges(&mut model, &block, &verticals, wall, T)
        .unwrap()
        .shape;
    let loop_edges: Vec<_> = ogeom_topo::explore_unique(&model, &rounded, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let (a, b) = ends(&model, e);
            (a.z - side).abs() < 1e-9 && (b.z - side).abs() < 1e-9
        })
        .collect();
    assert_eq!(loop_edges.len(), 8, "four sides and four quarter turns");
    let result = ogeom_fillet::fillet_edges(&mut model, &rounded, &loop_edges, top, T).unwrap();

    let diagnosis = ogeom_algo::check(&model, &result.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let faces = ogeom_topo::explore_unique(&model, &result.shape, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 6 + 4 + 8);

    let pi = core::f64::consts::PI;
    let section = (1.0 - pi / 4.0) * top * top;
    let centroid = (10.0 - 3.0 * pi) / (12.0 - 3.0 * pi) * top;
    let block = side * side * side - 4.0 * (4.0 - pi) * wall * wall / 4.0 * side;
    let straight = 4.0 * (side - 2.0 * wall) * section;
    let turns = 2.0 * pi * (wall - centroid) * section;
    let exact = block - straight - turns;
    let props =
        ogeom_algo::volume_properties(&model, &result.shape, ogeom_mesh::Deflection::default(), T)
            .unwrap();
    assert!(
        (props.mass - exact).abs() < 1e-6,
        "volume {} against {exact}",
        props.mass
    );
}
