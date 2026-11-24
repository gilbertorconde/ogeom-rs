//! Blend surfaces and blend curves bridging a gap, measured: the joins
//! against the faces and curves they continue, the shapes against the
//! closed forms of their Hermite columns, and the refusals by name.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use std::f64::consts::FRAC_PI_2;

use ogeom_algo::{
    attach_pcurve, check, make_edge_between, make_face, make_face_with_pcurves, make_polygon,
    make_wire, sew,
};
use ogeom_core::Tolerances;
use ogeom_fillet::{End, analyse_blend, make_blend_curve, make_blend_surface};
use ogeom_geom::{
    BSplineCurve, BSplineSurface, CircleCurve, Continuity, Curve, Curve3d as _, Line2d, LineCurve,
    PlaneSurface, Surface as _, SurfaceGeometry,
};
use ogeom_math::{
    Axis, Axis2, Circle, ControlGrid, Direction, Direction2, Frame, KnotVector, Plane, Point,
    Point2, Vector, Vector2,
};
use ogeom_topo::{EdgeRepr, Location, Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn dir(x: f64, y: f64, z: f64) -> Direction {
    Direction::new(Vector::new(x, y, z), T).unwrap()
}

/// The plane `z = height`, charted by `x` and `y`.
fn level(height: f64) -> SurfaceGeometry {
    let frame = Frame::new(
        Point::new(0.0, 0.0, height),
        dir(0.0, 0.0, 1.0),
        dir(1.0, 0.0, 0.0),
        T,
    )
    .unwrap();
    SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(frame)))
}

/// A rectangle on the plane `z = height`, `x` over `xs` and `y` over `ys`,
/// its edges trimmed on the plane.
fn strip(model: &mut Model, height: f64, xs: (f64, f64), ys: (f64, f64)) -> Shape {
    let points = [
        Point::new(xs.0, ys.0, height),
        Point::new(xs.1, ys.0, height),
        Point::new(xs.1, ys.1, height),
        Point::new(xs.0, ys.1, height),
    ];
    let wire = make_polygon(model, &points, true, T).unwrap().shape;
    let edges = model.ordered_children_of(&wire).unwrap();
    let face = make_face_with_pcurves(model, level(height), &[edges], T)
        .unwrap()
        .shape;
    assert_valid(model, &face);
    face
}

fn assert_valid(model: &Model, shape: &Shape) {
    let diagnosis = check(model, shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
}

fn vertex_points(model: &Model, edge: &Shape) -> Vec<Point> {
    explore_unique(model, edge, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| model.node(v).unwrap().data().as_vertex().unwrap().point)
        .collect()
}

/// The edge of `face` whose two ends are `a` and `b`, either way round.
fn edge_between(model: &Model, face: &Shape, a: Point, b: Point) -> Shape {
    explore_unique(model, face, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let ends = vertex_points(model, e);
            ends.len() == 2
                && ((ends[0].distance(a) < 1e-9 && ends[1].distance(b) < 1e-9)
                    || (ends[0].distance(b) < 1e-9 && ends[1].distance(a) < 1e-9))
        })
        .expect("no edge between those points")
}

/// The surface a face lies on.
fn surface_of(model: &Model, face: &Shape) -> SurfaceGeometry {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        panic!("not a face");
    };
    model.geometry().surface(data.surface).unwrap().clone()
}

/// The two parallel strips of the acceptance case: one on `z = 0` for
/// `y` in `[-5, 0]`, one on `z = 4` for `y` in `[10, 15]`, both `20` long
/// in `x`, with their facing edges.
fn two_strips(model: &mut Model) -> ((Shape, Shape), (Shape, Shape)) {
    let low = strip(model, 0.0, (0.0, 20.0), (-5.0, 0.0));
    let high = strip(model, 4.0, (0.0, 20.0), (10.0, 15.0));
    let low_edge = edge_between(
        model,
        &low,
        Point::new(0.0, 0.0, 0.0),
        Point::new(20.0, 0.0, 0.0),
    );
    let high_edge = edge_between(
        model,
        &high,
        Point::new(0.0, 10.0, 4.0),
        Point::new(20.0, 10.0, 4.0),
    );
    ((low_edge, low), (high_edge, high))
}

/// Blend, sew and measure: the blend shares both edges, the three faces
/// sew into one shell with those edges inside it, and each join is
/// reported by `analyse_blend`, as `(tangency, curvature)` per edge.
fn blend_and_measure(
    model: &mut Model,
    first: (&Shape, &Shape),
    second: (&Shape, &Shape),
    continuity: (Continuity, Continuity),
) -> (Shape, Vec<(f64, f64)>) {
    let blend = make_blend_surface(model, first, second, continuity, T)
        .unwrap()
        .shape;
    assert_valid(model, &blend);
    let edges = explore_unique(model, &blend, ShapeType::Edge).unwrap();
    assert_eq!(edges.len(), 4, "the blend is bounded by four edges");
    for given in [first.0, second.0] {
        assert!(
            edges.iter().any(|e| e.node() == given.node()),
            "the blend shares the given edge itself"
        );
    }
    let sewn = sew(
        model,
        &[first.1.clone(), blend.clone(), second.1.clone()],
        T,
    )
    .unwrap();
    assert_eq!(sewn.shells.len(), 1, "the three faces sew into one shell");
    for given in [first.0, second.0] {
        assert!(
            !sewn.free_edges.iter().any(|e| e.node() == given.node()),
            "a blended edge is inside the shell"
        );
    }
    let shell = sewn.shells[0].clone();
    let contacts = analyse_blend(model, &shell, &blend, 41, T).unwrap();
    let mut out = Vec::new();
    for given in [first.0, second.0] {
        let contact = contacts
            .iter()
            .find(|c| c.edge.node() == given.node())
            .expect("the blend's join along a given edge is measured");
        assert!(
            contact.gap < 1e-7,
            "the edge lies on both faces: {}",
            contact.gap
        );
        out.push((contact.tangency_error, contact.curvature_error));
    }
    (blend, out)
}

/// The acceptance case at G1: tangent to each strip along its edge, well
/// inside a tenth of a degree, and one shell with both.
#[test]
fn a_g1_blend_between_parallel_strips_is_tangent_to_both() {
    let mut model = Model::new();
    let (a, b) = two_strips(&mut model);
    let (blend, joins) = blend_and_measure(
        &mut model,
        (&a.0, &a.1),
        (&b.0, &b.1),
        (Continuity::G1, Continuity::G1),
    );
    for (tangency, curvature) in &joins {
        assert!(tangency.to_degrees() < 0.1, "tangency {tangency}");
        assert!(*tangency < 1e-9, "a join to a plane is exact: {tangency}");
        // A cubic across bends at its ends where the strips do not: G1 is
        // not G2 here, and the measure tells them apart.
        assert!(
            *curvature > 1e-2,
            "the G1 blend bends at its ends: {curvature}"
        );
    }

    // Each column across is the cubic Hermite with tangents the gap's width
    // long: Bezier points (0,0,0), (0,w/3,0), (0,10-w/3,4), (0,10,4).
    let surface = surface_of(&model, &blend);
    let w = 116.0f64.sqrt();
    let q = [
        Vector::new(0.0, 0.0, 0.0),
        Vector::new(0.0, w / 3.0, 0.0),
        Vector::new(0.0, 10.0 - w / 3.0, 4.0),
        Vector::new(0.0, 10.0, 4.0),
    ];
    // The blend runs along the first edge as that edge's curve does.
    let x0 = surface.point_at(0.0, 0.0, T).unwrap().x;
    for u in [0.0, 0.3, 1.0] {
        let x = 2.0f64.mul_add(-x0, 20.0).mul_add(u, x0);
        for v in [0.25f64, 0.5, 0.8] {
            let s = 1.0 - v;
            let c = q[0] * s.powi(3)
                + q[1] * (3.0 * s * s * v)
                + q[2] * (3.0 * s * v * v)
                + q[3] * v.powi(3);
            let expected = Point::new(x, c.y, c.z);
            let got = surface.point_at(u, v, T).unwrap();
            assert!(
                got.distance(expected) < 1e-9,
                "at ({u}, {v}): {got:?} against {expected:?}"
            );
        }
    }
}

/// The acceptance case with the high strip built on `z = 0` and placed 4
/// up: the blend is read where the strip stands, is the same surface as
/// the unplaced case's, and sews to the placed strip itself.
#[test]
fn a_blend_to_a_placed_strip_is_the_blend_to_where_it_stands() {
    let mut model = Model::new();
    let (a, b) = two_strips(&mut model);
    let truth = make_blend_surface(
        &mut model,
        (&a.0, &a.1),
        (&b.0, &b.1),
        (Continuity::G1, Continuity::G1),
        T,
    )
    .unwrap()
    .shape;
    let flat = strip(&mut model, 0.0, (0.0, 20.0), (10.0, 15.0));
    let high = model.placed(
        &flat,
        ogeom_math::Transform::translation(Vector::new(0.0, 0.0, 4.0)),
    );
    let high_edge = edge_between(
        &model,
        &high,
        Point::new(0.0, 10.0, 0.0),
        Point::new(20.0, 10.0, 0.0),
    );
    assert!(!high_edge.location().is_identity());
    let built = make_blend_surface(
        &mut model,
        (&a.0, &a.1),
        (&high_edge, &high),
        (Continuity::G1, Continuity::G1),
        T,
    )
    .unwrap();
    let blend = built.shape.clone();
    assert_valid(&model, &blend);
    assert!(
        !built.history.generated(&high_edge).is_empty(),
        "the placed edge generates the blend"
    );
    let (ours, theirs) = (surface_of(&model, &blend), surface_of(&model, &truth));
    for u in [0.0, 0.3, 1.0] {
        for v in [0.0, 0.25, 0.5, 1.0] {
            let (p, q) = (
                ours.point_at(u, v, T).unwrap(),
                theirs.point_at(u, v, T).unwrap(),
            );
            assert!(p.distance(q) < 1e-9, "at ({u}, {v}): {p:?} against {q:?}");
        }
    }

    let sewn = sew(&mut model, &[a.1.clone(), blend.clone(), high.clone()], T).unwrap();
    assert_eq!(sewn.shells.len(), 1, "the three faces sew into one shell");
    let blend = sewn
        .history
        .modified(&blend)
        .first()
        .cloned()
        .unwrap_or(blend);
    let contacts = analyse_blend(&model, &sewn.shells[0], &blend, 41, T).unwrap();
    assert_eq!(contacts.len(), 2, "the blend meets both strips");
    for contact in &contacts {
        assert!(contact.gap < 1e-7, "{contact:?}");
        assert!(contact.tangency_error < 1e-9, "{contact:?}");
    }
}

/// The acceptance case at G2: tangent, and bending as the strips do (not
/// at all) square to each edge.
#[test]
fn a_g2_blend_between_parallel_strips_is_curvature_continuous() {
    let mut model = Model::new();
    let (a, b) = two_strips(&mut model);
    let (blend, joins) = blend_and_measure(
        &mut model,
        (&a.0, &a.1),
        (&b.0, &b.1),
        (Continuity::G2, Continuity::G2),
    );
    for (tangency, curvature) in &joins {
        assert!(*tangency < 1e-9, "tangency {tangency}");
        assert!(*curvature < 1e-3, "curvature {curvature}");
        assert!(*curvature < 1e-9, "a join to a plane is exact: {curvature}");
    }
    // Independently of the curvature tensor: the second derivative across
    // has no part out of either strip's plane at its edge.
    let surface = surface_of(&model, &blend);
    for &u in &[0.0, 0.4, 1.0] {
        for v in [0.0, 1.0] {
            let (_, _, d2v) = surface.d2_at(u, v, T).unwrap();
            assert!(d2v.z.abs() < 1e-9, "d2v at ({u}, {v}): {d2v:?}");
        }
    }
    // The quintic is symmetric about the middle of the gap.
    let middle = surface.point_at(0.5, 0.5, T).unwrap();
    assert!(
        middle.distance(Point::new(10.0, 5.0, 2.0)) < 1e-9,
        "{middle:?}"
    );
}

/// Mixed continuity: the blend leaves the first strip at a corner and
/// meets the second tangent.
#[test]
fn a_c0_g1_blend_meets_one_side_at_a_crease_and_the_other_tangent() {
    let mut model = Model::new();
    let (a, b) = two_strips(&mut model);
    let (_, joins) = blend_and_measure(
        &mut model,
        (&a.0, &a.1),
        (&b.0, &b.1),
        (Continuity::C0, Continuity::G1),
    );
    assert!(
        joins[0].0 > 0.1,
        "a crease at the first strip: {}",
        joins[0].0
    );
    assert!(joins[1].0 < 1e-9, "tangent at the second: {}", joins[1].0);
}

/// The facing edges of the two strips run opposite ways, as their
/// polygons wind alike: the blend pairs their ends so its sides run
/// straight across the gap rather than crossing it, whichever way the
/// first edge is passed. History says the blend came from both edges.
#[test]
fn the_sides_pair_the_near_ends() {
    for flip in [false, true] {
        let mut model = Model::new();
        let (a, b) = two_strips(&mut model);
        let first = if flip { a.0.reversed() } else { a.0.clone() };
        let built = make_blend_surface(
            &mut model,
            (&first, &a.1),
            (&b.0, &b.1),
            (Continuity::G1, Continuity::G1),
            T,
        )
        .unwrap();
        let blend = built.shape.clone();
        assert_valid(&model, &blend);
        let sides: Vec<Shape> = explore_unique(&model, &blend, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .filter(|e| e.node() != a.0.node() && e.node() != b.0.node())
            .collect();
        assert_eq!(sides.len(), 2);
        for side in &sides {
            let ends = vertex_points(&model, side);
            assert!(
                (ends[0].x - ends[1].x).abs() < 1e-12,
                "a side across: {ends:?}"
            );
        }
        for given in [&a.0, &b.0] {
            let generated = built.history.generated(given);
            assert!(generated.iter().any(|s| s.node() == blend.node()));
            for side in &sides {
                assert!(generated.iter().any(|s| s.node() == side.node()));
            }
        }
    }
}

/// A face on `surface` bounded by four curves, side `k` running from
/// `corners[k]` to `corners[k + 1]`, its trims found in closed form.
fn quad_face(
    model: &mut Model,
    surface: SurfaceGeometry,
    corners: [Point; 4],
    sides: [(Curve, (f64, f64)); 4],
) -> (Shape, Vec<Shape>) {
    let vertices: Vec<Shape> = corners.iter().map(|p| model.add_point(*p)).collect();
    let mut ring = Vec::new();
    for (k, (curve, range)) in sides.into_iter().enumerate() {
        let (a, b) = (&vertices[k], &vertices[(k + 1) % 4]);
        ring.push(
            make_edge_between(model, curve, range, a, b, T)
                .unwrap()
                .shape,
        );
    }
    let face = make_face_with_pcurves(model, surface, &[ring.clone()], T)
        .unwrap()
        .shape;
    assert_valid(model, &face);
    (face, ring)
}

fn segment(a: Point, b: Point) -> (Curve, (f64, f64)) {
    let line = LineCurve::new(Axis {
        location: a,
        direction: Direction::new(b - a, T).unwrap(),
    });
    (line.into(), (0.0, a.distance(b)))
}

/// The circle of `radius` about the `z` axis at `height`.
fn quarter(radius: f64, height: f64) -> Circle {
    let frame = Frame::new(
        Point::new(0.0, 0.0, height),
        dir(0.0, 0.0, 1.0),
        dir(1.0, 0.0, 0.0),
        T,
    )
    .unwrap();
    Circle::new(frame, radius, T).unwrap()
}

/// A quarter annulus on `z = height` between radii `inner` and `outer`,
/// its inner arc an exact rational B-spline and its outer one too where
/// `rational`, else a circle; with its inner and outer arc edges.
fn quarter_annulus(
    model: &mut Model,
    height: f64,
    (inner, outer): (f64, f64),
    rational: bool,
) -> (Shape, Shape, Shape) {
    let spline = |radius: f64| {
        let circle: Curve = CircleCurve::new(quarter(radius, height)).into();
        circle.to_bspline_over((0.0, FRAC_PI_2), T).unwrap()
    };
    let outer_arc = if rational {
        (Curve::BSpline(spline(outer)), (0.0, 1.0))
    } else {
        (
            CircleCurve::new(quarter(outer, height)).into(),
            (0.0, FRAC_PI_2),
        )
    };
    let inner_arc = {
        let forward = spline(inner);
        let (k, c) = ogeom_math::bspline::reverse(forward.knots(), forward.control_points());
        (
            Curve::BSpline(BSplineCurve::rational(k, c).unwrap()),
            (0.0, 1.0),
        )
    };
    let corners = [
        Point::new(inner, 0.0, height),
        Point::new(outer, 0.0, height),
        Point::new(0.0, outer, height),
        Point::new(0.0, inner, height),
    ];
    let (face, ring) = quad_face(
        model,
        level(height),
        corners,
        [
            segment(corners[0], corners[1]),
            outer_arc,
            segment(corners[2], corners[3]),
            inner_arc,
        ],
    );
    (face, ring[3].clone(), ring[1].clone())
}

/// Rational edges: a quarter annulus's outer arc blended to a wider one's
/// inner arc above it, both stored as exact rational B-splines. The rows
/// beside each edge carry its weights, so the joins to the two planes are
/// still exact.
#[test]
fn rational_arc_edges_blend_exactly_to_their_planes() {
    let mut model = Model::new();
    let (low, _, low_outer) = quarter_annulus(&mut model, 0.0, (5.0, 10.0), true);
    let (high, high_inner, _) = quarter_annulus(&mut model, 3.0, (15.0, 20.0), true);
    for continuity in [Continuity::G1, Continuity::G2] {
        let (_, joins) = blend_and_measure(
            &mut model,
            (&low_outer, &low),
            (&high_inner, &high),
            (continuity, continuity),
        );
        for (tangency, curvature) in joins {
            assert!(tangency < 1e-9, "{continuity:?} tangency {tangency}");
            if continuity == Continuity::G2 {
                assert!(curvature < 1e-9, "curvature {curvature}");
            }
        }
    }
}

/// A face over the whole of a bicubic patch, bounded by its own
/// iso-curves; with its edge along `v = 1`.
fn patch_face(model: &mut Model, patch: &BSplineSurface) -> (Shape, Shape) {
    let surface = SurfaceGeometry::BSpline(patch.clone());
    let corners = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|(u, v)| model.add_point(patch.point_at(*u, *v, T).unwrap()))
        .collect();
    let unit = |x: f64, y: f64| Direction2::new(Vector2::new(x, y), T).unwrap();
    let sides = [
        (false, 0.0, 0, 1),
        (true, 1.0, 1, 2),
        (false, 1.0, 3, 2),
        (true, 0.0, 0, 3),
    ];
    let mut built = Vec::new();
    for (fixed_u, at, a, b) in sides {
        let curve = Curve::BSpline(if fixed_u {
            patch.iso_u_curve(at, T).unwrap()
        } else {
            patch.iso_v_curve(at, T).unwrap()
        });
        let edge = make_edge_between(model, curve, (0.0, 1.0), &vertices[a], &vertices[b], T)
            .unwrap()
            .shape;
        let (location, direction) = if fixed_u {
            (Point2::new(at, 0.0), unit(0.0, 1.0))
        } else {
            (Point2::new(0.0, at), unit(1.0, 0.0))
        };
        let trim = Line2d::over(
            Axis2 {
                location,
                direction,
            },
            0.0,
            1.0,
        )
        .unwrap();
        built.push((edge, trim));
    }
    let ring = vec![
        built[0].0.clone(),
        built[1].0.clone(),
        built[2].0.reversed(),
        built[3].0.reversed(),
    ];
    let wire = make_wire(model, &ring, T).unwrap().shape;
    let face = make_face(model, surface, &[wire], T).unwrap().shape;
    let Some(NodeData::Face(data)) = model.node(&face).map(|n| n.data()) else {
        panic!("not a face");
    };
    let id = data.surface;
    for (edge, trim) in &built {
        attach_pcurve(
            model,
            edge,
            (*trim).into(),
            id,
            Location::identity(),
            (0.0, 1.0),
        )
        .unwrap();
    }
    assert_valid(model, &face);
    (face, built[2].0.clone())
}

/// A bicubic hump over `x` in `[0, 20]`, `y` in `[-6, 0]`, its far side
/// (`v = 1`, along `y = 0`) curved in height and its tangent plane turning
/// along it.
fn hump() -> BSplineSurface {
    let heights = [
        [0.0, 0.5, 0.0, -0.5],
        [0.0, 2.0, 1.5, 1.0],
        [0.0, 1.0, 2.5, 0.5],
        [0.0, -0.5, 0.5, 0.0],
    ];
    let mut points = Vec::new();
    for (i, row) in heights.iter().enumerate() {
        for (j, z) in row.iter().enumerate() {
            points.push(Point::new(
                20.0 * f64::from(u8::try_from(i).unwrap()) / 3.0,
                -6.0 + 2.0 * f64::from(u8::try_from(j).unwrap()),
                *z,
            ));
        }
    }
    let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3).unwrap();
    BSplineSurface::new(
        knots.clone(),
        knots,
        &ControlGrid::new(points, 4, 4).unwrap(),
        T,
    )
    .unwrap()
}

/// A curved face: the patch's tangent plane turns along its edge, so the
/// crossing tangents are interpolated and the join is held to its budget by
/// refinement. `analyse_blend` measures it independently.
#[test]
fn a_blend_from_a_curved_patch_is_held_to_its_budget() {
    for continuity in [Continuity::G1, Continuity::G2] {
        let mut model = Model::new();
        let patch = hump();
        let (curved, edge) = patch_face(&mut model, &patch);
        let high = strip(&mut model, 6.0, (0.0, 20.0), (10.0, 15.0));
        let high_edge = edge_between(
            &model,
            &high,
            Point::new(0.0, 10.0, 6.0),
            Point::new(20.0, 10.0, 6.0),
        );
        let (_, joins) = blend_and_measure(
            &mut model,
            (&edge, &curved),
            (&high_edge, &high),
            (continuity, continuity),
        );
        let (tangency, curvature) = joins[0];
        assert!(tangency < 2e-5, "{continuity:?}: tangency {tangency}");
        assert!(tangency > 0.0, "the curved join is approximated, not exact");
        if continuity == Continuity::G2 {
            assert!(curvature < 1e-3, "curvature {curvature}");
        }
        assert!(joins[1].0 < 1e-9, "the planar side stays exact");
    }
}

fn refusal(result: ogeom_core::OgeomResult<ogeom_algo::Built>) -> String {
    match result {
        Ok(_) => panic!("expected a refusal"),
        Err(e) => e.to_string(),
    }
}

/// Every refusal names its reason.
#[test]
fn blend_surfaces_refuse_by_name() {
    let mut model = Model::new();
    let (a, b) = two_strips(&mut model);
    let g1 = Continuity::G1;

    for (continuity, says) in [
        (Continuity::C1, "parametric"),
        (Continuity::C2, "parametric"),
        (Continuity::CInfinity, "G2 at most"),
    ] {
        let message = refusal(make_blend_surface(
            &mut model,
            (&a.0, &a.1),
            (&b.0, &b.1),
            (g1, continuity),
            T,
        ));
        assert!(message.contains(says), "{continuity:?}: {message}");
    }

    let message = refusal(make_blend_surface(
        &mut model,
        (&a.0, &b.1),
        (&b.0, &b.1),
        (g1, g1),
        T,
    ));
    assert!(message.contains("not on the face"), "{message}");

    let message = refusal(make_blend_surface(
        &mut model,
        (&a.0, &a.1),
        (&a.0.reversed(), &a.1),
        (g1, g1),
        T,
    ));
    assert!(message.contains("same one"), "{message}");

    // Two edges of one strip meeting at its corner.
    let side = edge_between(
        &model,
        &a.1,
        Point::new(20.0, 0.0, 0.0),
        Point::new(20.0, -5.0, 0.0),
    );
    let message = refusal(make_blend_surface(
        &mut model,
        (&a.0, &a.1),
        (&side, &a.1),
        (g1, g1),
        T,
    ));
    assert!(message.contains("meet at an end"), "{message}");

    let message = refusal(make_blend_surface(
        &mut model,
        (&a.1, &a.1),
        (&b.0, &b.1),
        (g1, g1),
        T,
    ));
    assert!(message.contains("bridges edges"), "{message}");

    // A circle's B-spline form runs on another parameter than its angle.
    let (annulus, _, circle_edge) = quarter_annulus(&mut model, 0.0, (5.0, 10.0), false);
    let message = refusal(make_blend_surface(
        &mut model,
        (&circle_edge, &annulus),
        (&b.0, &b.1),
        (g1, g1),
        T,
    ));
    assert!(message.contains("parameter"), "{message}");

    // An edge with no trim on its face's surface.
    let points = [
        Point::new(0.0, -25.0, 0.0),
        Point::new(20.0, -25.0, 0.0),
        Point::new(20.0, -20.0, 0.0),
        Point::new(0.0, -20.0, 0.0),
    ];
    let wire = make_polygon(&mut model, &points, true, T).unwrap().shape;
    let bare = make_face(&mut model, level(0.0), &[wire], T).unwrap().shape;
    let edge = edge_between(&model, &bare, points[2], points[3]);
    let message = refusal(make_blend_surface(
        &mut model,
        (&edge, &bare),
        (&b.0, &b.1),
        (g1, g1),
        T,
    ));
    assert!(message.contains("no trim"), "{message}");
}

/// A line from `a` to `b` as an edge.
fn line_edge(model: &mut Model, a: Point, b: Point) -> Shape {
    let (va, vb) = (model.add_point(a), model.add_point(b));
    let (curve, range) = segment(a, b);
    make_edge_between(model, curve, range, &va, &vb, T)
        .unwrap()
        .shape
}

/// The blend edge's curve and range.
fn curve_of(model: &Model, edge: &Shape) -> (Curve, (f64, f64)) {
    let data = model.node(edge).unwrap().data().as_edge().unwrap();
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        panic!("no curve");
    };
    (model.geometry().curve(*curve).unwrap().clone(), *range)
}

/// Two collinear lines facing each other: the G1 blend continues both, so
/// it is the straight run between them, at uniform speed.
#[test]
fn a_g1_blend_curve_between_collinear_lines_is_straight() {
    let mut model = Model::new();
    let a = line_edge(&mut model, Point::new(-10.0, 0.0, 0.0), Point::ORIGIN);
    let b = line_edge(
        &mut model,
        Point::new(20.0, 0.0, 0.0),
        Point::new(10.0, 0.0, 0.0),
    );
    let built = make_blend_curve(
        &mut model,
        (&a, End::End),
        (&b, End::End),
        (Continuity::G1, Continuity::G1),
        T,
    )
    .unwrap();
    let (curve, range) = curve_of(&model, &built.shape);
    assert_eq!(range, (0.0, 1.0));
    for k in 0..=10 {
        let t = f64::from(k) / 10.0;
        let p = curve.point_at(t, T).unwrap();
        assert!(
            p.distance(Point::new(10.0 * t, 0.0, 0.0)) < 1e-12,
            "{t}: {p:?}"
        );
    }
    // It joins the edges' own vertices.
    let ends = explore_unique(&model, &built.shape, ShapeType::Vertex).unwrap();
    for (edge, at) in [(&a, Point::ORIGIN), (&b, Point::new(10.0, 0.0, 0.0))] {
        let shared = explore_unique(&model, edge, ShapeType::Vertex)
            .unwrap()
            .into_iter()
            .find(|v| model.node(v).unwrap().data().as_vertex().unwrap().point == at)
            .unwrap();
        assert!(ends.iter().any(|v| v.node() == shared.node()));
    }
}

/// An arc of `radius` about `centre` in the `xy` plane, angles `from` to
/// `to`, as an edge.
fn arc_edge(model: &mut Model, centre: Point, radius: f64, from: f64, to: f64) -> Shape {
    let frame = Frame::new(centre, dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0), T).unwrap();
    let curve: Curve = CircleCurve::new(Circle::new(frame, radius, T).unwrap()).into();
    let (a, b) = (
        curve.point_at(from, T).unwrap(),
        curve.point_at(to, T).unwrap(),
    );
    let (va, vb) = (model.add_point(a), model.add_point(b));
    make_edge_between(model, curve, (from, to), &va, &vb, T)
        .unwrap()
        .shape
}

/// The unit tangent and curvature vector of a curve at `t`.
fn frame_at(curve: &Curve, t: f64) -> (Vector, Vector) {
    let d = curve.derivatives_at(t, 2, T).unwrap();
    let speed = d[1].magnitude();
    let tangent = d[1] * (1.0 / speed);
    let bend = (d[2] - tangent * d[2].dot(tangent)) * (1.0 / (speed * speed));
    (tangent, bend)
}

/// Between two arcs, the G2 blend leaves each with its tangent and its
/// curvature vector: `1/r` towards the arc's centre. The G1 blend keeps the
/// tangents and not the curvature.
#[test]
fn a_g2_blend_curve_carries_each_arcs_curvature() {
    let mut model = Model::new();
    // Ends at (0, 5) heading in -x, and at (-12, 4) where the second arc,
    // about (-12, 0) of radius 4, begins heading in -x.
    let a = arc_edge(&mut model, Point::ORIGIN, 5.0, 0.0, FRAC_PI_2);
    let b = arc_edge(&mut model, Point::new(-12.0, 0.0, 0.0), 4.0, FRAC_PI_2, 2.5);
    for continuity in [Continuity::G1, Continuity::G2] {
        let built = make_blend_curve(
            &mut model,
            (&a, End::End),
            (&b, End::Start),
            (continuity, continuity),
            T,
        )
        .unwrap();
        let (curve, _) = curve_of(&model, &built.shape);
        let (t0, k0) = frame_at(&curve, 0.0);
        let (t1, k1) = frame_at(&curve, 1.0);
        assert!(
            curve
                .point_at(0.0, T)
                .unwrap()
                .distance(Point::new(0.0, 5.0, 0.0))
                < 1e-12
        );
        assert!(
            curve
                .point_at(1.0, T)
                .unwrap()
                .distance(Point::new(-12.0, 4.0, 0.0))
                < 1e-12
        );
        assert!(
            (t0 - Vector::new(-1.0, 0.0, 0.0)).magnitude() < 1e-12,
            "{t0:?}"
        );
        assert!(
            (t1 - Vector::new(-1.0, 0.0, 0.0)).magnitude() < 1e-12,
            "{t1:?}"
        );
        let (want0, want1) = (Vector::new(0.0, -0.2, 0.0), Vector::new(0.0, -0.25, 0.0));
        if continuity == Continuity::G2 {
            assert!((k0 - want0).magnitude() < 1e-12, "{k0:?}");
            assert!((k1 - want1).magnitude() < 1e-12, "{k1:?}");
        } else {
            assert!((k0 - want0).magnitude() > 1e-2, "{k0:?}");
        }
    }
}

/// The end named is the end of the edge as passed: the start of a reversed
/// edge is the end of its curve.
#[test]
fn a_reversed_edge_starts_where_its_curve_ends() {
    let mut model = Model::new();
    let a = arc_edge(&mut model, Point::ORIGIN, 5.0, 0.0, FRAC_PI_2);
    let b = line_edge(
        &mut model,
        Point::new(-12.0, 4.0, 0.0),
        Point::new(-20.0, 4.0, 0.0),
    );
    let g2 = (Continuity::G2, Continuity::G1);
    let forward = make_blend_curve(&mut model, (&a, End::End), (&b, End::Start), g2, T).unwrap();
    let flipped = make_blend_curve(
        &mut model,
        (&a.reversed(), End::Start),
        (&b, End::Start),
        g2,
        T,
    )
    .unwrap();
    let (one, _) = curve_of(&model, &forward.shape);
    let (other, _) = curve_of(&model, &flipped.shape);
    for k in 0..=8 {
        let t = f64::from(k) / 8.0;
        assert!(
            one.point_at(t, T)
                .unwrap()
                .distance(other.point_at(t, T).unwrap())
                < 1e-12
        );
    }
}

#[test]
fn blend_curves_refuse_by_name() {
    let mut model = Model::new();
    let a = line_edge(&mut model, Point::new(-10.0, 0.0, 0.0), Point::ORIGIN);
    let b = line_edge(
        &mut model,
        Point::new(20.0, 0.0, 0.0),
        Point::new(10.0, 0.0, 0.0),
    );
    let message = refusal(make_blend_curve(
        &mut model,
        (&a, End::End),
        (&b, End::End),
        (Continuity::G1, Continuity::CInfinity),
        T,
    ));
    assert!(message.contains("G2 at most"), "{message}");
    let message = refusal(make_blend_curve(
        &mut model,
        (&a, End::End),
        (&b, End::End),
        (Continuity::C1, Continuity::G1),
        T,
    ));
    assert!(message.contains("parametric"), "{message}");
    let message = refusal(make_blend_curve(
        &mut model,
        (&a, End::End),
        (&a, End::End),
        (Continuity::G1, Continuity::G1),
        T,
    ));
    assert!(message.contains("one point"), "{message}");
    let vertex = model.add_point(Point::ORIGIN);
    let message = refusal(make_blend_curve(
        &mut model,
        (&vertex, End::End),
        (&b, End::End),
        (Continuity::G1, Continuity::G1),
        T,
    ));
    assert!(message.contains("ends of edges"), "{message}");
}
