//! A face grown past one of its edges: measured against the closed forms of
//! the surfaces it grows on, and refused by name where it cannot grow.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use std::f64::consts::FRAC_PI_2;

use ogeom_algo::{
    Extension, attach_pcurve, check, extend_face, make_edge_between, make_face,
    make_face_with_pcurves, make_polygon, make_wire, surface_properties,
};
use ogeom_core::Tolerances;
use ogeom_geom::{
    BSplineSurface, CircleCurve, ConeSurface, Continuity, Curve, CylinderSurface, Line2d,
    LineCurve, PlaneSurface, SphereSurface, Surface as _, SurfaceGeometry,
};
use ogeom_math::{
    Axis, Axis2, Circle, Cone, ControlGrid, Cylinder, Direction, Direction2, Frame, KnotVector,
    Plane, Point, Point2, Sphere, Transform, Vector, Vector2,
};
use ogeom_mesh::Deflection;
use ogeom_topo::{EdgeRepr, Filter, Location, Model, NodeData, Shape, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

fn fine() -> Deflection {
    Deflection {
        chord: 1e-4,
        angular: 0.02,
        ..Deflection::default()
    }
}

fn dir(x: f64, y: f64, z: f64) -> Direction {
    Direction::new(Vector::new(x, y, z), T).unwrap()
}

/// A planar face on the `xy` plane bounded by the polygon through `points`,
/// its edges trimmed on the plane.
fn polygon_face(model: &mut Model, points: &[(f64, f64)]) -> Shape {
    let points: Vec<Point> = points
        .iter()
        .map(|(x, y)| Point::new(*x, *y, 0.0))
        .collect();
    let wire = make_polygon(model, &points, true, T).unwrap().shape;
    let edges = model.ordered_children_of(&wire).unwrap();
    let plane = SurfaceGeometry::Plane(PlaneSurface::new(Plane::XY));
    let face = make_face_with_pcurves(model, plane, &[edges], T)
        .unwrap()
        .shape;
    assert_valid(model, &face);
    face
}

/// The edges of a face, in the order its outer wire walks them.
fn outer_edges(model: &Model, face: &Shape) -> Vec<Shape> {
    let wire = model.children_of(face).unwrap()[0].clone();
    model.ordered_children_of(&wire).unwrap()
}

/// Every vertex point of a shape.
fn vertex_points(model: &Model, shape: &Shape) -> Vec<Point> {
    explore(model, shape, Filter::OfType(ShapeType::Vertex))
        .unwrap()
        .iter()
        .map(|v| model.node(v).unwrap().data().as_vertex().unwrap().point)
        .collect()
}

/// The edge whose two ends are `a` and `b`, either way round.
fn edge_between(model: &Model, face: &Shape, a: Point, b: Point) -> Shape {
    outer_edges(model, face)
        .into_iter()
        .find(|e| {
            let ends = vertex_points(model, e);
            ends.len() == 2
                && ((ends[0].distance(a) < 1e-9 && ends[1].distance(b) < 1e-9)
                    || (ends[0].distance(b) < 1e-9 && ends[1].distance(a) < 1e-9))
        })
        .expect("no edge between those points")
}

fn surface_of(model: &Model, face: &Shape) -> SurfaceGeometry {
    let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
        panic!("not a face");
    };
    model.geometry().surface(data.surface).unwrap().clone()
}

fn area(model: &Model, face: &Shape) -> f64 {
    surface_properties(model, face, fine(), T).unwrap().mass
}

fn assert_valid(model: &Model, face: &Shape) {
    let diagnosis = check(model, face, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
}

/// A face on `surface` bounded by four curves, side `k` running from
/// `corners[k]` to `corners[k + 1]`, or stored the other way where its flag
/// says so; the trims are found in closed form.
fn quad_face(
    model: &mut Model,
    surface: SurfaceGeometry,
    corners: [Point; 4],
    sides: [(Curve, (f64, f64), bool); 4],
) -> Shape {
    let vertices: Vec<Shape> = corners.iter().map(|p| model.add_point(*p)).collect();
    let mut ring = Vec::new();
    for (k, (curve, range, backwards)) in sides.into_iter().enumerate() {
        let (a, b) = (&vertices[k], &vertices[(k + 1) % 4]);
        let edge = if backwards {
            make_edge_between(model, curve, range, b, a, T)
                .unwrap()
                .shape
                .reversed()
        } else {
            make_edge_between(model, curve, range, a, b, T)
                .unwrap()
                .shape
        };
        ring.push(edge);
    }
    let face = make_face_with_pcurves(model, surface, &[ring], T)
        .unwrap()
        .shape;
    assert_valid(model, &face);
    face
}

fn arc(centre: Point, radius: f64, from: f64, to: f64) -> (Curve, (f64, f64), bool) {
    let frame = Frame::new(centre, dir(0.0, 0.0, 1.0), dir(1.0, 0.0, 0.0), T).unwrap();
    let circle = CircleCurve::new(Circle::new(frame, radius, T).unwrap());
    (circle.into(), (from, to), false)
}

fn segment(a: Point, b: Point, backwards: bool) -> (Curve, (f64, f64), bool) {
    let (start, end) = if backwards { (b, a) } else { (a, b) };
    let line = LineCurve::new(Axis {
        location: start,
        direction: Direction::new(end - start, T).unwrap(),
    });
    (line.into(), (0.0, start.distance(end)), backwards)
}

/// A quarter of a radius-10 cylinder about `z`, angles `0..pi/2`, heights
/// `0..20`, on a surface whose height range is exactly that.
fn quarter_cylinder(model: &mut Model) -> Shape {
    let r = 10.0;
    let surface = SurfaceGeometry::Cylinder(
        CylinderSurface::new(Cylinder::new(Frame::WORLD, r, T).unwrap(), (0.0, 20.0)).unwrap(),
    );
    let c = [
        Point::new(r, 0.0, 0.0),
        Point::new(0.0, r, 0.0),
        Point::new(0.0, r, 20.0),
        Point::new(r, 0.0, 20.0),
    ];
    let mut top = arc(Point::new(0.0, 0.0, 20.0), r, 0.0, FRAC_PI_2);
    top.2 = true;
    quad_face(
        model,
        surface,
        c,
        [
            arc(Point::ORIGIN, r, 0.0, FRAC_PI_2),
            segment(c[1], c[2], false),
            top,
            segment(c[3], c[0], false),
        ],
    )
}

/// A face over the chart rectangle `u` by `v` of a B-spline patch, bounded
/// by the patch's own iso-curves with straight trims.
fn patch_face(
    model: &mut Model,
    patch: &BSplineSurface,
    (u0, u1): (f64, f64),
    (v0, v1): (f64, f64),
) -> Shape {
    let surface = SurfaceGeometry::BSpline(patch.clone());
    let corners = [(u0, v0), (u1, v0), (u1, v1), (u0, v1)];
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|(u, v)| model.add_point(patch.point_at(*u, *v, T).unwrap()))
        .collect();
    let unit = |x: f64, y: f64| Direction2::new(Vector2::new(x, y), T).unwrap();
    // Each side as (fixed u?, fixed value, range, stored start, end).
    let sides = [
        (false, v0, (u0, u1), 0, 1),
        (true, u1, (v0, v1), 1, 2),
        (false, v1, (u0, u1), 3, 2),
        (true, u0, (v0, v1), 0, 3),
    ];
    let mut built = Vec::new();
    for (fixed_u, at, range, a, b) in sides {
        let curve = Curve::BSpline(if fixed_u {
            patch.iso_u_curve(at, T).unwrap()
        } else {
            patch.iso_v_curve(at, T).unwrap()
        });
        let edge = make_edge_between(model, curve, range, &vertices[a], &vertices[b], T)
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
            range.0,
            range.1,
        )
        .unwrap();
        built.push((edge, trim, range));
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
    for (edge, trim, range) in built {
        attach_pcurve(model, &edge, trim.into(), id, Location::identity(), range).unwrap();
    }
    assert_valid(model, &face);
    face
}

/// A fitted saddle `z = x^2 - y^2` over `[-1, 1]` squared.
fn fitted_saddle() -> BSplineSurface {
    let rows: Vec<Vec<Point>> = (0..9)
        .map(|j| {
            (0..9)
                .map(|i| {
                    let (x, y) = (-1.0 + 0.25 * f64::from(i), -1.0 + 0.25 * f64::from(j));
                    Point::new(x, y, x * x - y * y)
                })
                .collect()
        })
        .collect();
    ogeom_geom::fit::fit_surface_grid(&rows, 3, 1e-7, T)
        .unwrap()
        .curve
}

/// The saddle face, its patch, and its edge along the patch's `u = end`.
fn saddle_face(model: &mut Model) -> (Shape, BSplineSurface, Shape) {
    let patch = fitted_saddle();
    let ((ua, ub), (va, vb)) = patch.domain();
    let face = patch_face(model, &patch, (ua, ub), (va, vb));
    let edge = edge_between(
        model,
        &face,
        patch.point_at(ub, va, T).unwrap(),
        patch.point_at(ub, vb, T).unwrap(),
    );
    (face, patch, edge)
}

fn grown_patch(model: &Model, face: &Shape) -> BSplineSurface {
    match surface_of(model, face) {
        SurfaceGeometry::BSpline(b) => b,
        other => panic!("expected a patch, got {:?}", other.kind()),
    }
}

/// The length of the `v = at` line of a patch over `u` in `range`, as a
/// fine polyline.
fn polyline_length(patch: &BSplineSurface, at: f64, (a, b): (f64, f64)) -> f64 {
    const N: u32 = 4000;
    let mut total = 0.0;
    let mut last = patch.point_at(a, at, T).unwrap();
    for k in 1..=N {
        let u = a + (b - a) * f64::from(k) / f64::from(N);
        let p = patch.point_at(u, at, T).unwrap();
        total += last.distance(p);
        last = p;
    }
    total
}

fn refusal(result: ogeom_core::OgeomResult<ogeom_algo::Built>) -> String {
    match result {
        Ok(_) => panic!("expected a refusal"),
        Err(e) => e.to_string(),
    }
}

/// A 10 by 10 square grown across one edge by 5 is 10 by 15: its corners
/// are where the closed form puts them, its area is 150, and the three
/// edges it kept are the same nodes.
#[test]
fn a_square_grown_across_one_edge_by_five_is_ten_by_fifteen() {
    let mut model = Model::new();
    let square = polygon_face(
        &mut model,
        &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
    );
    let edge = edge_between(
        &model,
        &square,
        Point::new(10.0, 0.0, 0.0),
        Point::new(10.0, 10.0, 0.0),
    );
    let built = extend_face(&mut model, &square, &edge, 5.0, Extension::Natural, T).unwrap();
    let grown = built.shape.clone();
    assert_valid(&model, &grown);
    assert!((area(&model, &grown) - 150.0).abs() < 1e-9);
    let corners = vertex_points(&model, &grown);
    for expected in [(0.0, 0.0), (15.0, 0.0), (15.0, 10.0), (0.0, 10.0)] {
        let p = Point::new(expected.0, expected.1, 0.0);
        assert!(
            corners.iter().any(|q| q.distance(p) < 1e-12),
            "no corner at {p:?}"
        );
    }
    for q in &corners {
        assert!(q.x > -1e-12 && q.x < 15.0 + 1e-12 && q.y > -1e-12 && q.y < 10.0 + 1e-12);
    }
    let old = outer_edges(&model, &square);
    let new = outer_edges(&model, &grown);
    let kept = old.iter().filter(|e| !e.is_same(&edge)).collect::<Vec<_>>();
    assert_eq!(kept.len(), 3);
    for e in kept {
        assert!(new.iter().any(|n| n.is_same(e)), "a kept edge was rebuilt");
    }
    assert!(!new.iter().any(|n| n.is_same(&edge)));

    // History: the face became the grown one, the edge the far edge at
    // x = 15, and the edge generated the two sides.
    assert_eq!(built.history.modified(&square).len(), 1);
    assert!(built.history.modified(&square)[0].is_same(&grown));
    let far = &built.history.modified(&edge)[0];
    for p in vertex_points(&model, far) {
        assert!((p.x - 15.0).abs() < 1e-12);
    }
    let sides = built.history.generated(&edge);
    assert_eq!(sides.len(), 2);
    for side in sides {
        assert!(new.iter().any(|n| n.is_same(side)));
    }
}

/// The face's side is read from its boundary, not from how the boundary
/// happens to be wound: a clockwise square grows away from itself too,
/// square to a slanted edge.
#[test]
fn a_clockwise_face_grows_away_from_itself_across_a_slanted_edge() {
    let mut model = Model::new();
    // A unit-diagonal diamond walked clockwise.
    let s = 10.0 / 2.0_f64.sqrt();
    let points = [(0.0, 0.0), (-s, s), (0.0, 2.0 * s), (s, s)];
    let diamond = polygon_face(&mut model, &points);
    let edge = edge_between(
        &model,
        &diamond,
        Point::new(s, s, 0.0),
        Point::new(0.0, 0.0, 0.0),
    );
    let grown = extend_face(&mut model, &diamond, &edge, 5.0, Extension::Natural, T)
        .unwrap()
        .shape;
    assert_valid(&model, &grown);
    assert!((area(&model, &grown) - 150.0).abs() < 1e-9);
    // The edge's outward normal is (1, -1) / sqrt 2.
    let out = Vector::new(1.0, -1.0, 0.0) * (5.0 / 2.0_f64.sqrt());
    let corners = vertex_points(&model, &grown);
    for p in [Point::new(s, s, 0.0) + out, Point::ORIGIN + out] {
        assert!(
            corners.iter().any(|q| q.distance(p) < 1e-12),
            "no corner at {p:?}"
        );
    }
}

/// A face of a box grown past one of its edges leaves the box alone: the
/// solid is as valid as before, and the neighbour across the edge still
/// holds it.
#[test]
fn a_box_face_grown_leaves_its_neighbours_untouched() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let faces = explore(&model, &solid, Filter::OfType(ShapeType::Face)).unwrap();
    let face = faces[0].clone();
    let edge = outer_edges(&model, &face)[0].clone();
    let before = area(&model, &face);
    let grown = extend_face(&mut model, &face, &edge, 3.0, Extension::Natural, T)
        .unwrap()
        .shape;
    assert_valid(&model, &grown);
    assert!((area(&model, &grown) - (before + 30.0)).abs() < 1e-9);
    assert_valid(&model, &solid);
    let holders = faces
        .iter()
        .filter(|f| outer_edges(&model, f).iter().any(|e| e.is_same(&edge)))
        .count();
    assert_eq!(holders, 2);
}

/// A quarter cylinder grown across its straight edge by 5 of arc stays on
/// the same cylinder: every new vertex is the radius from the axis, the new
/// ruling stands 5 / r radians on, and the area is r times the angle times
/// the height.
#[test]
fn a_quarter_cylinder_grown_across_its_ruling_stays_on_the_cylinder() {
    let mut model = Model::new();
    let face = quarter_cylinder(&mut model);
    let edge = edge_between(
        &model,
        &face,
        Point::new(0.0, 10.0, 0.0),
        Point::new(0.0, 10.0, 20.0),
    );
    let grown = extend_face(&mut model, &face, &edge, 5.0, Extension::Natural, T)
        .unwrap()
        .shape;
    assert_valid(&model, &grown);
    let SurfaceGeometry::Cylinder(c) = surface_of(&model, &grown) else {
        panic!("left the cylinder");
    };
    assert_eq!(c.cylinder(), Cylinder::new(Frame::WORLD, 10.0, T).unwrap());
    let angle = FRAC_PI_2 + 0.5;
    let far = [
        Point::new(10.0 * angle.cos(), 10.0 * angle.sin(), 0.0),
        Point::new(10.0 * angle.cos(), 10.0 * angle.sin(), 20.0),
    ];
    let corners = vertex_points(&model, &grown);
    for p in far {
        assert!(
            corners.iter().any(|q| q.distance(p) < 1e-9),
            "no corner at {p:?}"
        );
    }
    for q in &corners {
        assert!((q.x.hypot(q.y) - 10.0).abs() < 1e-9);
    }
    // The arc from the old ruling to the new, measured from the corners.
    let new_angle = far[0].y.atan2(far[0].x);
    assert!((10.0 * (new_angle - FRAC_PI_2) - 5.0).abs() < 1e-9);
    assert!((area(&model, &grown) - 10.0 * angle * 20.0).abs() < 1e-6);
}

/// Grown across its bottom arc, the quarter cylinder runs down its rulings:
/// the surface's height range widens to hold it, the trims of the edges it
/// kept carry over, and a linear extension (straight along the rulings) is
/// the same face.
#[test]
fn a_quarter_cylinder_grown_down_its_rulings_widens_its_surface() {
    for mode in [
        Extension::Natural,
        Extension::Linear {
            continuity: Continuity::G1,
        },
    ] {
        let mut model = Model::new();
        let face = quarter_cylinder(&mut model);
        let edge = edge_between(
            &model,
            &face,
            Point::new(10.0, 0.0, 0.0),
            Point::new(0.0, 10.0, 0.0),
        );
        let grown = extend_face(&mut model, &face, &edge, 5.0, mode, T)
            .unwrap()
            .shape;
        assert_valid(&model, &grown);
        let SurfaceGeometry::Cylinder(c) = surface_of(&model, &grown) else {
            panic!("left the cylinder");
        };
        assert_eq!(c.domain().1, (-5.0, 20.0));
        let corners = vertex_points(&model, &grown);
        for p in [Point::new(10.0, 0.0, -5.0), Point::new(0.0, 10.0, -5.0)] {
            assert!(
                corners.iter().any(|q| q.distance(p) < 1e-9),
                "no corner at {p:?}"
            );
        }
        assert!((area(&model, &grown) - 10.0 * FRAC_PI_2 * 25.0).abs() < 1e-6);
    }
}

/// A band of a sphere grown up across its top parallel climbs its
/// meridians by 5 / R radians of latitude, staying on the sphere.
#[test]
fn a_sphere_band_grown_across_its_parallel_climbs_its_meridians() {
    let mut model = Model::new();
    let r = 10.0;
    let face = sphere_band(&mut model, r, 0.5);
    let top = 0.5_f64;
    let edge = edge_between(
        &model,
        &face,
        Point::new(r * top.cos(), 0.0, r * top.sin()),
        Point::new(0.0, r * top.cos(), r * top.sin()),
    );
    let grown = extend_face(&mut model, &face, &edge, 3.0, Extension::Natural, T)
        .unwrap()
        .shape;
    assert_valid(&model, &grown);
    let lat = top + 3.0 / r;
    let corners = vertex_points(&model, &grown);
    for p in [
        Point::new(r * lat.cos(), 0.0, r * lat.sin()),
        Point::new(0.0, r * lat.cos(), r * lat.sin()),
    ] {
        assert!(
            corners.iter().any(|q| q.distance(p) < 1e-9),
            "no corner at {p:?}"
        );
    }
    for q in &corners {
        assert!((q.to_vector().magnitude() - r).abs() < 1e-9);
    }
    // A zone of a sphere: R squared, times the longitude span, times the
    // difference of the sines of its latitudes.
    let expected = r * r * FRAC_PI_2 * lat.sin();
    assert!((area(&model, &grown) - expected).abs() < 1e-6);
}

/// A quarter band of a radius-10 sphere from the equator up to `top`.
fn sphere_band(model: &mut Model, r: f64, top: f64) -> Shape {
    let surface =
        SurfaceGeometry::Sphere(SphereSurface::new(Sphere::new(Frame::WORLD, r, T).unwrap()));
    let c = [
        Point::new(r, 0.0, 0.0),
        Point::new(0.0, r, 0.0),
        Point::new(0.0, r * top.cos(), r * top.sin()),
        Point::new(r * top.cos(), 0.0, r * top.sin()),
    ];
    let meridian = |x: Direction, z: Direction| {
        let frame = Frame::new(Point::ORIGIN, z, x, T).unwrap();
        Curve::from(CircleCurve::new(Circle::new(frame, r, T).unwrap()))
    };
    let mut parallel = arc(
        Point::new(0.0, 0.0, r * top.sin()),
        r * top.cos(),
        0.0,
        FRAC_PI_2,
    );
    parallel.2 = true;
    quad_face(
        model,
        surface,
        c,
        [
            arc(Point::ORIGIN, r, 0.0, FRAC_PI_2),
            (
                meridian(dir(0.0, 1.0, 0.0), dir(1.0, 0.0, 0.0)),
                (0.0, top),
                false,
            ),
            parallel,
            (
                meridian(dir(1.0, 0.0, 0.0), dir(0.0, -1.0, 0.0)),
                (0.0, top),
                true,
            ),
        ],
    )
}

/// A quarter of a cone grown across its top circle runs up its rulings by
/// the length along the slant: the new rim stands 5 cos(half angle) higher,
/// on the cone.
#[test]
fn a_cone_grown_across_its_rim_runs_up_its_slant() {
    let mut model = Model::new();
    let (face, cone) = quarter_cone(&mut model);
    let r1 = cone.radius_at(10.0);
    let edge = edge_between(
        &model,
        &face,
        Point::new(r1, 0.0, 10.0),
        Point::new(0.0, r1, 10.0),
    );
    let grown = extend_face(&mut model, &face, &edge, 5.0, Extension::Natural, T)
        .unwrap()
        .shape;
    assert_valid(&model, &grown);
    let h = 10.0 + 5.0 * cone.half_angle().cos();
    let r2 = cone.radius_at(h);
    let corners = vertex_points(&model, &grown);
    for p in [Point::new(r2, 0.0, h), Point::new(0.0, r2, h)] {
        assert!(
            corners.iter().any(|q| q.distance(p) < 1e-9),
            "no corner at {p:?}"
        );
    }
    // The slant from old rim to new is the length asked for.
    assert!((Point::new(r1, 0.0, 10.0).distance(Point::new(r2, 0.0, h)) - 5.0).abs() < 1e-9);
}

/// A quarter of a cone of base radius 10 and half angle 0.3, heights
/// `0..10`.
fn quarter_cone(model: &mut Model) -> (Shape, Cone) {
    let cone = Cone::new(Frame::WORLD, 10.0, 0.3, T).unwrap();
    let surface = SurfaceGeometry::Cone(ConeSurface::new(cone, (0.0, 10.0)).unwrap());
    let (r0, r1) = (cone.radius_at(0.0), cone.radius_at(10.0));
    let c = [
        Point::new(r0, 0.0, 0.0),
        Point::new(0.0, r0, 0.0),
        Point::new(0.0, r1, 10.0),
        Point::new(r1, 0.0, 10.0),
    ];
    let mut top = arc(Point::new(0.0, 0.0, 10.0), r1, 0.0, FRAC_PI_2);
    top.2 = true;
    let face = quad_face(
        model,
        surface,
        c,
        [
            arc(Point::ORIGIN, r0, 0.0, FRAC_PI_2),
            segment(c[1], c[2], false),
            top,
            segment(c[3], c[0], false),
        ],
    );
    (face, cone)
}

/// A fitted saddle grown naturally across one edge by 2 is G2 with the
/// original across the old edge's place: position, first and second
/// derivatives agree from both sides. The original keeps its points, and
/// the crossing line through the edge's middle is 2 long, measured as a
/// fine polyline.
#[test]
fn a_fitted_saddle_grown_naturally_is_g2_across_the_old_edge() {
    let mut model = Model::new();
    let (face, patch, edge) = saddle_face(&mut model);
    let ((ua, ub), (va, vb)) = patch.domain();
    let grown = extend_face(&mut model, &face, &edge, 2.0, Extension::Natural, T)
        .unwrap()
        .shape;
    assert_valid(&model, &grown);
    let longer = grown_patch(&model, &grown);
    let ((la, lb), (lva, lvb)) = longer.domain();
    assert_eq!((la, lva, lvb), (ua, va, vb));
    assert!(lb > ub);
    for i in 0..=6 {
        for j in 0..=6 {
            let u = ua + (ub - ua) * f64::from(i) / 6.0;
            let v = va + (vb - va) * f64::from(j) / 6.0;
            let a = patch.point_at(u, v, T).unwrap();
            let b = longer.point_at(u, v, T).unwrap();
            assert!(a.distance(b) < 1e-9, "the original moved at ({u}, {v})");
        }
    }
    let step = 1e-6;
    for j in 0..=8 {
        let v = va + (vb - va) * f64::from(j) / 8.0;
        let inside = longer.d2_at(ub - step, v, T).unwrap();
        let outside = longer.d2_at(ub + step, v, T).unwrap();
        let first_in = longer.d1_at(ub - step, v, T).unwrap();
        let first_out = longer.d1_at(ub + step, v, T).unwrap();
        for (x, y) in [
            (first_in.0, first_out.0),
            (first_in.1, first_out.1),
            (inside.0, outside.0),
            (inside.1, outside.1),
            (inside.2, outside.2),
        ] {
            let jump = (x - y).magnitude();
            assert!(
                jump < 1e-3 * (1.0 + x.magnitude()),
                "a derivative jumps by {jump} at v = {v}"
            );
        }
    }
    let middle = 0.5 * (va + vb);
    assert!((polyline_length(&longer, middle, (ub, lb)) - 2.0).abs() < 1e-6);
}

/// A linear extension at G1 runs straight out along the crossing tangent:
/// every point of the extension is on the tangent line from the old edge,
/// the first derivative carries across, and the curvature the saddle has
/// across the edge does not.
#[test]
fn a_saddle_grown_linearly_at_g1_runs_straight_out() {
    let mut model = Model::new();
    let (face, patch, edge) = saddle_face(&mut model);
    let ((_, ub), (va, vb)) = patch.domain();
    let mode = Extension::Linear {
        continuity: Continuity::G1,
    };
    let grown = extend_face(&mut model, &face, &edge, 2.0, mode, T)
        .unwrap()
        .shape;
    assert_valid(&model, &grown);
    let longer = grown_patch(&model, &grown);
    let lb = longer.domain().0.1;
    for j in 0..=8 {
        let v = va + (vb - va) * f64::from(j) / 8.0;
        let base = patch.point_at(ub, v, T).unwrap();
        let (tangent, _) = patch.d1_at(ub, v, T).unwrap();
        let unit = tangent * (1.0 / tangent.magnitude());
        for k in 1..=5 {
            let u = ub + (lb - ub) * f64::from(k) / 5.0;
            let off = longer.point_at(u, v, T).unwrap() - base;
            let across = (off - unit * off.dot(unit)).magnitude();
            assert!(across < 1e-9, "the extension bends by {across}");
        }
        let (first_out, _) = longer.d1_at(ub + 1e-6, v, T).unwrap();
        assert!((first_out - tangent).magnitude() < 1e-4 * tangent.magnitude());
        let bend_in = patch.d2_at(ub, v, T).unwrap().0;
        let bend_out = longer.d2_at(ub + 1e-6, v, T).unwrap().0;
        assert!(bend_out.magnitude() < 1e-6);
        assert!(bend_in.magnitude() > 0.1);
    }
    let middle = 0.5 * (va + vb);
    assert!((polyline_length(&longer, middle, (ub, lb)) - 2.0).abs() < 1e-6);
}

/// A linear extension at G2 carries the second derivative across too.
#[test]
fn a_saddle_grown_linearly_at_g2_keeps_its_curvature_across() {
    let mut model = Model::new();
    let (face, patch, edge) = saddle_face(&mut model);
    let ((_, ub), (va, vb)) = patch.domain();
    let mode = Extension::Linear {
        continuity: Continuity::G2,
    };
    let grown = extend_face(&mut model, &face, &edge, 2.0, mode, T)
        .unwrap()
        .shape;
    assert_valid(&model, &grown);
    let longer = grown_patch(&model, &grown);
    for j in 0..=8 {
        let v = va + (vb - va) * f64::from(j) / 8.0;
        let bend_in = longer.d2_at(ub - 1e-6, v, T).unwrap().0;
        let bend_out = longer.d2_at(ub + 1e-6, v, T).unwrap().0;
        assert!((bend_in - bend_out).magnitude() < 1e-3 * bend_in.magnitude());
    }
}

#[test]
fn a_length_that_is_not_positive_is_refused() {
    let mut model = Model::new();
    let square = polygon_face(
        &mut model,
        &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
    );
    let edge = outer_edges(&model, &square)[0].clone();
    for length in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let err = refusal(extend_face(
            &mut model,
            &square,
            &edge,
            length,
            Extension::Natural,
            T,
        ));
        assert!(err.contains("positive length"), "{err}");
    }
}

#[test]
fn an_edge_not_on_the_face_is_refused() {
    let mut model = Model::new();
    let square = polygon_face(
        &mut model,
        &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
    );
    let other = polygon_face(&mut model, &[(20.0, 0.0), (30.0, 0.0), (30.0, 10.0)]);
    let stray = outer_edges(&model, &other)[0].clone();
    let err = refusal(extend_face(
        &mut model,
        &square,
        &stray,
        1.0,
        Extension::Natural,
        T,
    ));
    assert!(err.contains("not on the face"), "{err}");
}

#[test]
fn a_shape_that_is_not_a_face_is_refused() {
    let mut model = Model::new();
    let square = polygon_face(
        &mut model,
        &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
    );
    let edge = outer_edges(&model, &square)[0].clone();
    let wire = model.children_of(&square).unwrap()[0].clone();
    let err = refusal(extend_face(
        &mut model,
        &wire,
        &edge,
        1.0,
        Extension::Natural,
        T,
    ));
    assert!(err.contains("extends a face"), "{err}");
}

#[test]
fn an_edge_bounding_a_hole_is_refused() {
    let mut model = Model::new();
    let ring = |pts: &[(f64, f64)]| -> Vec<Point> {
        pts.iter().map(|(x, y)| Point::new(*x, *y, 0.0)).collect()
    };
    let outer = make_polygon(
        &mut model,
        &ring(&[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]),
        true,
        T,
    )
    .unwrap()
    .shape;
    let hole = make_polygon(
        &mut model,
        &ring(&[(4.0, 4.0), (4.0, 6.0), (6.0, 6.0), (6.0, 4.0)]),
        true,
        T,
    )
    .unwrap()
    .shape;
    let plane = SurfaceGeometry::Plane(PlaneSurface::new(Plane::XY));
    let face = make_face(&mut model, plane, &[outer, hole.clone()], T)
        .unwrap()
        .shape;
    let inner = model.children_of(&hole).unwrap()[0].clone();
    let err = refusal(extend_face(
        &mut model,
        &face,
        &inner,
        1.0,
        Extension::Natural,
        T,
    ));
    assert!(err.contains("bounds a hole"), "{err}");
}

#[test]
fn a_curved_edge_of_a_planar_face_is_refused() {
    let mut model = Model::new();
    // A half disc: a diameter and a semicircle.
    let a = model.add_point(Point::new(10.0, 0.0, 0.0));
    let b = model.add_point(Point::new(-10.0, 0.0, 0.0));
    let (curve, range, _) = arc(Point::ORIGIN, 10.0, 0.0, std::f64::consts::PI);
    let round = make_edge_between(&mut model, curve, range, &a, &b, T)
        .unwrap()
        .shape;
    let (line, span, _) = segment(
        Point::new(-10.0, 0.0, 0.0),
        Point::new(10.0, 0.0, 0.0),
        false,
    );
    let straight = make_edge_between(&mut model, line, span, &b, &a, T)
        .unwrap()
        .shape;
    let wire = make_wire(&mut model, &[round.clone(), straight], T)
        .unwrap()
        .shape;
    let plane = SurfaceGeometry::Plane(PlaneSurface::new(Plane::XY));
    let face = make_face(&mut model, plane, &[wire], T).unwrap().shape;
    let err = refusal(extend_face(
        &mut model,
        &face,
        &round,
        1.0,
        Extension::Natural,
        T,
    ));
    assert!(err.contains("straight edge"), "{err}");
}

#[test]
fn a_placed_face_is_refused() {
    let mut model = Model::new();
    let square = polygon_face(
        &mut model,
        &[(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)],
    );
    let edge = outer_edges(&model, &square)[0].clone();
    let moved = model.placed(&square, Transform::translation(Vector::new(0.0, 0.0, 5.0)));
    let err = refusal(extend_face(
        &mut model,
        &moved,
        &edge,
        1.0,
        Extension::Natural,
        T,
    ));
    assert!(err.contains("placed face"), "{err}");
}

/// The extension that would carry a quarter cylinder more than the rest of
/// the way round is refused: the face would overlap itself.
#[test]
fn an_extension_wrapping_a_periodic_surface_onto_the_face_is_refused() {
    let mut model = Model::new();
    let face = quarter_cylinder(&mut model);
    let edge = edge_between(
        &model,
        &face,
        Point::new(0.0, 10.0, 0.0),
        Point::new(0.0, 10.0, 20.0),
    );
    // Three quarters of the circumference is 47.1.
    let err = refusal(extend_face(
        &mut model,
        &face,
        &edge,
        48.0,
        Extension::Natural,
        T,
    ));
    assert!(err.contains("wrap the face"), "{err}");
}

/// A linear extension across a ruling of a cylinder would run off along
/// the tangent plane; it is refused rather than given on the cylinder.
#[test]
fn a_linear_extension_across_a_curved_analytic_line_is_refused() {
    let mut model = Model::new();
    let face = quarter_cylinder(&mut model);
    let edge = edge_between(
        &model,
        &face,
        Point::new(0.0, 10.0, 0.0),
        Point::new(0.0, 10.0, 20.0),
    );
    let mode = Extension::Linear {
        continuity: Continuity::G1,
    };
    let err = refusal(extend_face(&mut model, &face, &edge, 1.0, mode, T));
    assert!(err.contains("would leave the surface"), "{err}");
}

#[test]
fn a_linear_extension_smooth_to_every_order_is_refused() {
    let mut model = Model::new();
    let (face, _, edge) = saddle_face(&mut model);
    let mode = Extension::Linear {
        continuity: Continuity::CInfinity,
    };
    let err = refusal(extend_face(&mut model, &face, &edge, 1.0, mode, T));
    assert!(err.contains("second order at most"), "{err}");
}

#[test]
fn an_extension_reaching_a_sphere_pole_is_refused() {
    let mut model = Model::new();
    let face = sphere_band(&mut model, 10.0, 0.5);
    let top = 0.5_f64;
    let edge = edge_between(
        &model,
        &face,
        Point::new(10.0 * top.cos(), 0.0, 10.0 * top.sin()),
        Point::new(0.0, 10.0 * top.cos(), 10.0 * top.sin()),
    );
    // The pole is 10 (pi/2 - 0.5) = 10.7 away along the meridian.
    let err = refusal(extend_face(
        &mut model,
        &face,
        &edge,
        11.0,
        Extension::Natural,
        T,
    ));
    assert!(err.contains("sphere's pole"), "{err}");
}

#[test]
fn an_extension_reaching_a_cone_apex_is_refused() {
    let mut model = Model::new();
    let (face, cone) = quarter_cone(&mut model);
    let r0 = cone.radius_at(0.0);
    let edge = edge_between(
        &model,
        &face,
        Point::new(r0, 0.0, 0.0),
        Point::new(0.0, r0, 0.0),
    );
    // The apex is r0 / sin(half angle) = 33.8 down the slant.
    let err = refusal(extend_face(
        &mut model,
        &face,
        &edge,
        40.0,
        Extension::Natural,
        T,
    ));
    assert!(err.contains("cone's apex"), "{err}");
}

/// A face trimmed inside its patch has no side of the patch at its edge to
/// continue from.
#[test]
fn a_patch_edge_inside_the_domain_is_refused() {
    let mut model = Model::new();
    let patch = fitted_saddle();
    let ((ua, ub), (va, vb)) = patch.domain();
    let (u1, v1) = (0.5 * (ua + ub), 0.5 * (va + vb));
    let face = patch_face(&mut model, &patch, (ua, u1), (va, v1));
    let edge = edge_between(
        &model,
        &face,
        patch.point_at(u1, va, T).unwrap(),
        patch.point_at(u1, v1, T).unwrap(),
    );
    let err = refusal(extend_face(
        &mut model,
        &face,
        &edge,
        1.0,
        Extension::Natural,
        T,
    ));
    assert!(err.contains("not on a side of the patch"), "{err}");
}

/// A patch whose net closes on itself across the edge would be continued
/// over its own start.
#[test]
fn a_patch_closing_across_the_edge_is_refused() {
    let mut model = Model::new();
    // A triangular tube, linear in `u`, its last column its first.
    let ring = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(10.0, 0.0, 0.0),
        Point::new(0.0, 10.0, 0.0),
        Point::new(0.0, 0.0, 0.0),
    ];
    let mut points = Vec::new();
    for p in ring {
        points.push(p);
        points.push(p + Vector::new(0.0, 0.0, 10.0));
    }
    let patch = BSplineSurface::new(
        KnotVector::new(vec![0.0, 0.0, 1.0, 2.0, 3.0, 3.0], 1).unwrap(),
        KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap(),
        &ControlGrid::new(points, 4, 2).unwrap(),
        T,
    )
    .unwrap();
    let face = patch_face(&mut model, &patch, (2.0, 3.0), (0.0, 1.0));
    let edge = edge_between(
        &model,
        &face,
        Point::new(0.0, 0.0, 0.0),
        Point::new(0.0, 0.0, 10.0),
    );
    let err = refusal(extend_face(
        &mut model,
        &face,
        &edge,
        1.0,
        Extension::Natural,
        T,
    ));
    assert!(err.contains("closes on itself"), "{err}");
}

/// The seam of a whole cylinder's side is walked twice by its face: there
/// is no one side to grow it across.
#[test]
fn a_seam_edge_is_refused() {
    let mut model = Model::new();
    let solid = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T)
        .unwrap()
        .shape;
    let side = explore(&model, &solid, Filter::OfType(ShapeType::Face))
        .unwrap()
        .into_iter()
        .find(|f| matches!(surface_of(&model, f), SurfaceGeometry::Cylinder(_)))
        .unwrap();
    let edges = outer_edges(&model, &side);
    let seam = edges
        .iter()
        .find(|e| edges.iter().filter(|o| o.is_same(e)).count() == 2)
        .unwrap()
        .clone();
    let err = refusal(extend_face(
        &mut model,
        &side,
        &seam,
        1.0,
        Extension::Natural,
        T,
    ));
    assert!(err.contains("twice"), "{err}");
}

/// A face on a surface swept from a curve is not one this extends.
#[test]
fn a_face_on_an_unsupported_surface_is_refused() {
    let mut model = Model::new();
    let base: Curve = LineCurve::new(Axis {
        location: Point::ORIGIN,
        direction: dir(1.0, 0.0, 0.0),
    })
    .into();
    let swept = SurfaceGeometry::Extrusion(Box::new(
        ogeom_geom::ExtrusionSurface::new(base, dir(0.0, 1.0, 0.0), 10.0).unwrap(),
    ));
    let points: Vec<Point> = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
        .iter()
        .map(|(x, y)| Point::new(*x, *y, 0.0))
        .collect();
    let wire = make_polygon(&mut model, &points, true, T).unwrap().shape;
    let face = make_face(&mut model, swept, &[wire], T).unwrap().shape;
    let edge = outer_edges(&model, &face)[0].clone();
    let err = refusal(extend_face(
        &mut model,
        &face,
        &edge,
        1.0,
        Extension::Natural,
        T,
    ));
    assert!(err.contains("is not extended; only planes"), "{err}");
}

/// The trims of the edges a grown patch kept are carried onto the new
/// patch, and still follow their curves there.
#[test]
fn the_kept_edges_of_a_grown_patch_are_trimmed_on_it() {
    let mut model = Model::new();
    let (face, _, edge) = saddle_face(&mut model);
    let grown = extend_face(&mut model, &face, &edge, 2.0, Extension::Natural, T)
        .unwrap()
        .shape;
    let Some(NodeData::Face(data)) = model.node(&grown).map(|n| n.data()) else {
        panic!("not a face");
    };
    let id = data.surface;
    let surface = surface_of(&model, &grown);
    for e in outer_edges(&model, &grown) {
        let data = model.node(&e).unwrap().data().as_edge().unwrap();
        let Some(EdgeRepr::PCurve { curve, range, .. }) = data.pcurve_on(id) else {
            panic!("an edge has no trim on the grown patch");
        };
        let Some(EdgeRepr::Curve3d {
            curve: c3,
            range: r3,
            ..
        }) = data.curve3d()
        else {
            panic!("an edge has no curve");
        };
        assert_eq!(range, r3);
        let trim = model.geometry().pcurve(*curve).unwrap();
        let path = model.geometry().curve(*c3).unwrap();
        for k in 0..=4 {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 4.0;
            use ogeom_geom::{Curve2d as _, Curve3d as _};
            let uv = trim.point_at(t, T).unwrap();
            let on = surface.point_at(uv.x, uv.y, T).unwrap();
            assert!(on.distance(path.point_at(t, T).unwrap()) < 1e-7);
        }
    }
}
