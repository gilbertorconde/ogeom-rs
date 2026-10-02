//! The N-sided filling, held to independent measurements: the face is
//! bounded by the very edges it was given, sews to the faces they came
//! from, stands within its tolerance of them (projected afresh, not read
//! from the report), and meets its supports at the continuity asked
//! (measured by the blend analysis).
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_algo::{
    check, edge_vertices, make_box, make_compound, make_edge_between, make_face_with_pcurves,
    make_revolution, make_solid, make_vertex, project_on_surface, sew, volume_properties,
};
use ogeom_core::{OgeomError, Tolerances};
use ogeom_geom::{
    CircleCurve, Continuity, Curve3d as _, CylinderSurface, LineCurve, PlaneSurface, Surface as _,
    SurfaceGeometry,
};
use ogeom_math::{Axis, Circle, Cylinder, Direction, Frame, Plane, Point, Vector};
use ogeom_offset::{FillBoundary, Filled, make_filling_n};
use ogeom_topo::{EdgeRepr, Filter, Model, NodeData, Shape, ShapeType, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn fine() -> ogeom_mesh::Deflection {
    ogeom_mesh::Deflection {
        chord: 1e-4,
        angular: 0.02,
        ..ogeom_mesh::Deflection::default()
    }
}

fn line(model: &mut Model, from: &Shape, to: &Shape) -> Shape {
    let at = |v: &Shape| model.node(v).unwrap().data().as_vertex().unwrap().point;
    let (a, b) = (at(from), at(to));
    make_edge_between(
        model,
        LineCurve::segment(a, b, T).unwrap().into(),
        (0.0, a.distance(b)),
        from,
        to,
        T,
    )
    .unwrap()
    .shape
}

/// A planar face on the loop `edges` (each already oriented to run head to
/// tail), its plane facing the side the loop turns counter-clockwise about.
fn planar_face(model: &mut Model, edges: &[Shape]) -> Shape {
    let mut points = Vec::new();
    for e in edges {
        let (a, _) = edge_vertices(model, e).unwrap().unwrap();
        points.push(model.node(&a).unwrap().data().as_vertex().unwrap().point);
    }
    let mut area = Vector::ZERO;
    for k in 0..points.len() {
        area += (points[k] - points[0]).cross(points[(k + 1) % points.len()] - points[0]);
    }
    let z = Direction::new(area, T).unwrap();
    let frame = Frame::new(points[0], z, z.any_perpendicular(), T).unwrap();
    make_face_with_pcurves(
        model,
        SurfaceGeometry::Plane(PlaneSurface::new(Plane::new(frame))),
        &[edges.to_vec()],
        T,
    )
    .unwrap()
    .shape
}

/// A prism standing under a loop of top corners (counter-clockwise seen
/// from above): its walls and floor, outward, and the top edges the hole
/// is bounded by, edge `k` running from corner `k` to corner `k + 1`.
struct OpenPrism {
    top: Vec<Shape>,
    walls: Vec<Shape>,
    floor: Shape,
}

fn open_prism(model: &mut Model, corners: &[Point], floor_z: f64) -> OpenPrism {
    let n = corners.len();
    let up: Vec<Shape> = corners
        .iter()
        .map(|c| make_vertex(model, *c).shape)
        .collect();
    let down: Vec<Shape> = corners
        .iter()
        .map(|c| make_vertex(model, Point::new(c.x, c.y, floor_z)).shape)
        .collect();
    let top: Vec<Shape> = (0..n)
        .map(|k| line(model, &up[k], &up[(k + 1) % n]))
        .collect();
    let rise: Vec<Shape> = (0..n).map(|k| line(model, &down[k], &up[k])).collect();
    let base: Vec<Shape> = (0..n)
        .map(|k| line(model, &down[k], &down[(k + 1) % n]))
        .collect();
    let walls: Vec<Shape> = (0..n)
        .map(|k| {
            let ring = [
                base[k].clone(),
                rise[(k + 1) % n].clone(),
                top[k].reversed(),
                rise[k].reversed(),
            ];
            planar_face(model, &ring)
        })
        .collect();
    let ring: Vec<Shape> = (0..n).rev().map(|k| base[k].reversed()).collect();
    let floor = planar_face(model, &ring);
    OpenPrism { top, walls, floor }
}

/// The face's surface.
fn surface_of(model: &Model, face: &Shape) -> SurfaceGeometry {
    let NodeData::Face(data) = model.node(face).unwrap().data() else {
        panic!("a face");
    };
    model.geometry().surface(data.surface).unwrap().clone()
}

/// The largest distance from points along `edge` to the face's surface,
/// found by projecting each point afresh.
fn projected_gap(model: &Model, face: &Shape, edge: &Shape) -> f64 {
    let surface = surface_of(model, face);
    let data = model.node(edge).unwrap().data().as_edge().unwrap();
    let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
        panic!("an edge with a curve");
    };
    let curve = model.geometry().curve(*curve).unwrap();
    let mut worst = 0.0f64;
    for k in 0..=40 {
        let t = range.0 + (range.1 - range.0) * f64::from(k) / 40.0;
        let p = curve.point_at(t, T).unwrap();
        worst = worst.max(project_on_surface(&surface, p, 32, T).unwrap().distance);
    }
    worst
}

fn face_edges(model: &Model, face: &Shape) -> Vec<Shape> {
    explore_unique(model, face, ShapeType::Edge).unwrap()
}

/// The face's boundary edges are exactly `edges`, as nodes.
fn bounded_by(model: &Model, face: &Shape, edges: &[Shape]) -> bool {
    let mine = face_edges(model, face);
    mine.len() == edges.len()
        && edges
            .iter()
            .all(|e| mine.iter().any(|m| m.node() == e.node()))
}

fn g0(edge: &Shape, support: Option<&Shape>) -> FillBoundary {
    FillBoundary {
        edge: edge.clone(),
        support: support.cloned(),
        continuity: Continuity::C0,
    }
}

fn saddle_corners() -> Vec<Point> {
    [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
        .iter()
        .map(|&(x, y)| Point::new(x, y, x * y))
        .collect()
}

fn closed_solid(model: &mut Model, faces: &[Shape]) -> Shape {
    let sewn = sew(model, faces, T).unwrap();
    assert_eq!(sewn.shells.len(), 1, "the faces sew into one shell");
    assert!(
        sewn.free_edges.is_empty(),
        "the shell closes: {} free edges",
        sewn.free_edges.len()
    );
    make_solid(model, &sewn.shells).unwrap().shape
}

#[test]
fn a_saddle_hole_fills_with_one_face_on_its_own_edges_and_sews_shut() {
    // z = x·y over [-1, 1]²: four straight edges, each on a vertical wall.
    let mut model = Model::new();
    let prism = open_prism(&mut model, &saddle_corners(), -2.0);
    let boundary: Vec<FillBoundary> = prism
        .top
        .iter()
        .zip(&prism.walls)
        .map(|(e, w)| g0(e, Some(w)))
        .collect();
    let tolerance = 1e-5;
    let filled = make_filling_n(&mut model, &boundary, &[], tolerance, T).unwrap();
    let face = filled.built.shape.clone();

    assert_eq!(model.kind_of(&face).unwrap(), ShapeType::Face);
    assert!(
        bounded_by(&model, &face, &prism.top),
        "the face is bounded by the given edge nodes, not copies"
    );
    for edge in &prism.top {
        let gap = projected_gap(&model, &face, edge);
        assert!(gap <= tolerance, "an edge stands {gap} off the filling");
        assert!(
            !filled.built.history.generated(edge).is_empty(),
            "each edge generates the face"
        );
    }
    for side in &filled.sides {
        assert!(side.gap <= tolerance, "{side:?}");
    }

    // Every edge carries a pcurve on the filling's surface.
    let NodeData::Face(data) = model.node(&face).unwrap().data() else {
        panic!("a face");
    };
    let id = data.surface;
    for edge in explore(&model, &face, Filter::OfType(ShapeType::Edge)).unwrap() {
        let data = model.node(&edge).unwrap().data().as_edge().unwrap();
        assert!(
            matches!(
                data.pcurve_for(id, edge.location()),
                Some(EdgeRepr::PCurve { .. })
            ),
            "an edge without a pcurve on the filling"
        );
    }

    // Sewn with the walls it came from and the floor: one closed shell, a
    // solid that checks clean and faces outward.
    let mut faces = prism.walls.clone();
    faces.push(prism.floor.clone());
    faces.push(face);
    let solid = closed_solid(&mut model, &faces);
    let diagnosis = check(&model, &solid, T).unwrap();
    assert!(diagnosis.is_usable(), "{:?}", diagnosis.problems);
    let volume = volume_properties(&model, &solid, fine(), T).unwrap().mass;
    assert!(volume > 0.0, "the shell faces outward: volume {volume}");
}

#[test]
fn interior_points_pin_the_saddle_to_its_closed_form() {
    // The same hole with a grid of vertices on z = x·y inside it, and one
    // of the saddle's straight rulings, x = 0.625: the surface passes
    // through every one.
    let mut model = Model::new();
    let prism = open_prism(&mut model, &saddle_corners(), -2.0);
    let mut constraints = Vec::new();
    for i in -3..=3 {
        for j in -3..=3 {
            let (x, y) = (f64::from(i) / 4.0, f64::from(j) / 4.0);
            constraints.push(make_vertex(&mut model, Point::new(x, y, x * y)).shape);
        }
    }
    let ends =
        [-0.875, 0.875].map(|y| make_vertex(&mut model, Point::new(0.625, y, 0.625 * y)).shape);
    let ruling = line(&mut model, &ends[0], &ends[1]);
    constraints.push(ruling.clone());
    let boundary: Vec<FillBoundary> = prism.top.iter().map(|e| g0(e, None)).collect();
    let tolerance = 1e-4;
    let filled = make_filling_n(&mut model, &boundary, &constraints, tolerance, T).unwrap();
    assert!(filled.constraint_gap <= tolerance);
    let gap = projected_gap(&model, &filled.built.shape, &ruling);
    assert!(gap <= tolerance, "the ruling stands {gap} off the surface");
    let surface = surface_of(&model, &filled.built.shape);
    for v in &constraints[..constraints.len() - 1] {
        let p = model.node(v).unwrap().data().as_vertex().unwrap().point;
        let d = project_on_surface(&surface, p, 32, T).unwrap().distance;
        assert!(d <= tolerance, "a constraint stands {d} off the surface");
    }
    // The saddle lies in the patch's spline space and meets every
    // condition, so between the pinned points the surface is the saddle
    // too, and the solid under it holds the integral of x·y + 2 over the
    // square: 8.
    for (x, y) in [(0.125, 0.375), (-0.625, 0.125), (0.375, -0.375)] {
        let foot = project_on_surface(&surface, Point::new(x, y, x * y), 32, T).unwrap();
        assert!(
            foot.distance < tolerance,
            "({x}, {y}) is {} off",
            foot.distance
        );
    }
    let mut faces = prism.walls.clone();
    faces.push(prism.floor.clone());
    faces.push(filled.built.shape.clone());
    let solid = closed_solid(&mut model, &faces);
    let volume = volume_properties(&model, &solid, fine(), T).unwrap().mass;
    assert!((volume - 8.0).abs() < 1e-5, "volume {volume} against 8");
}

#[test]
fn an_open_boxs_missing_side_fills_flat_and_closes_the_box() {
    let mut model = Model::new();
    let solid = make_box(&mut model, Frame::WORLD, (3.0, 2.0, 1.5), T)
        .unwrap()
        .shape;
    let faces = explore_unique(&model, &solid, ShapeType::Face).unwrap();
    let on_top = |model: &Model, face: &Shape| {
        explore_unique(model, face, ShapeType::Vertex)
            .unwrap()
            .iter()
            .all(|v| {
                (model.node(v).unwrap().data().as_vertex().unwrap().point.z - 1.5).abs() < 1e-12
            })
    };
    let (top, sides): (Vec<Shape>, Vec<Shape>) = faces.into_iter().partition(|f| on_top(&model, f));
    assert_eq!((top.len(), sides.len()), (1, 5));
    let edges = face_edges(&model, &top[0]);
    let boundary: Vec<FillBoundary> = edges
        .iter()
        .map(|e| {
            let support = sides
                .iter()
                .find(|f| face_edges(&model, f).iter().any(|x| x.node() == e.node()))
                .unwrap();
            g0(e, Some(support))
        })
        .collect();
    let filled = make_filling_n(&mut model, &boundary, &[], 1e-6, T).unwrap();
    let face = filled.built.shape.clone();
    assert!(bounded_by(&model, &face, &edges));

    // The fill is the plane z = 1.5 everywhere over its patch.
    let SurfaceGeometry::BSpline(patch) = surface_of(&model, &face) else {
        panic!("a fitted patch");
    };
    let ((ua, ub), (va, vb)) = patch.domain();
    for i in 0..=8 {
        for j in 0..=8 {
            let p = patch
                .point_at(
                    ua + (ub - ua) * f64::from(i) / 8.0,
                    va + (vb - va) * f64::from(j) / 8.0,
                    T,
                )
                .unwrap();
            assert!((p.z - 1.5).abs() < 1e-9, "{p:?}");
        }
    }

    let mut all = sides.clone();
    all.push(face);
    let closed = closed_solid(&mut model, &all);
    let diagnosis = check(&model, &closed, T).unwrap();
    assert!(diagnosis.is_usable(), "{:?}", diagnosis.problems);
    let volume = volume_properties(&model, &closed, fine(), T).unwrap().mass;
    assert!(
        (volume - 9.0).abs() < 1e-6,
        "volume {volume} against 3·2·1.5"
    );
}

/// A rounded box corner's three quarter-cylinders of radius `r` and length
/// `len`, round the x, y and z axes, each ending on a coordinate plane in a
/// quarter circle; the three arcs share their corner vertices and bound the
/// octant of the sphere of radius `r` at the origin. The faces, and the
/// arcs in the same order.
fn rounded_corner(model: &mut Model, r: f64, len: f64) -> (Vec<Shape>, Vec<Shape>) {
    let dirs = [Direction::X, Direction::Y, Direction::Z];
    let at = |d: Direction| Point::ORIGIN + d.vector() * r;
    let corners: Vec<Shape> = dirs
        .iter()
        .map(|d| make_vertex(model, at(*d)).shape)
        .collect();
    let arc =
        |model: &mut Model, origin: Point, a: Direction, b: Direction, from: &Shape, to: &Shape| {
            let circle =
                CircleCurve::new(Circle::new(Frame::new(origin, a, b, T).unwrap(), r, T).unwrap());
            make_edge_between(
                model,
                circle.into(),
                (0.0, std::f64::consts::FRAC_PI_2),
                from,
                to,
                T,
            )
            .unwrap()
            .shape
        };
    let (mut faces, mut arcs) = (Vec::new(), Vec::new());
    for k in 0..3 {
        let (a, b, c) = (dirs[k], dirs[(k + 1) % 3], dirs[(k + 2) % 3]);
        let back = a.vector() * -len;
        let near = arc(
            model,
            Point::ORIGIN,
            a,
            b,
            &corners[(k + 1) % 3],
            &corners[(k + 2) % 3],
        );
        let qb = make_vertex(model, at(b) + back).shape;
        let qc = make_vertex(model, at(c) + back).shape;
        let far = arc(model, Point::ORIGIN + back, a, b, &qb, &qc);
        let down = line(model, &corners[(k + 2) % 3], &qc);
        let up = line(model, &qb, &corners[(k + 1) % 3]);
        let cylinder = Cylinder::new(Frame::new(Point::ORIGIN, a, b, T).unwrap(), r, T).unwrap();
        let surface =
            SurfaceGeometry::Cylinder(CylinderSurface::new(cylinder, (-len, 0.0)).unwrap());
        let face = make_face_with_pcurves(
            model,
            surface,
            &[vec![near.clone(), down, far.reversed(), up]],
            T,
        )
        .unwrap()
        .shape;
        faces.push(face);
        arcs.push(near);
    }
    (faces, arcs)
}

#[test]
fn a_corner_between_three_quarter_cylinders_fills_tangent_to_each() {
    // Three quarter-cylinders of radius 2 round the three axes, each ending
    // on a coordinate plane in a quarter circle; the three arcs bound the
    // octant of the sphere a rounded box corner would carry.
    let mut model = Model::new();
    let (cylinders, arcs) = rounded_corner(&mut model, 2.0, 3.0);
    let boundary: Vec<FillBoundary> = arcs
        .iter()
        .zip(&cylinders)
        .map(|(edge, support)| FillBoundary {
            edge: edge.clone(),
            support: Some(support.clone()),
            continuity: Continuity::G1,
        })
        .collect();
    let tolerance = 1e-3;
    let filled: Filled = make_filling_n(&mut model, &boundary, &[], tolerance, T).unwrap();
    let face = filled.built.shape.clone();
    assert!(bounded_by(&model, &face, &arcs));

    let mut all = cylinders.clone();
    all.push(face.clone());
    let all = make_compound(&mut model, &all).unwrap().shape;
    let contacts = ogeom_fillet::analyse_blend(&model, &all, &face, 200, T).unwrap();
    assert_eq!(
        contacts.len(),
        3,
        "the fill meets each cylinder along its arc"
    );
    let limit = 0.1f64.to_radians();
    for contact in &contacts {
        assert!(
            contact.tangency_error < limit,
            "the fill meets a cylinder {} degrees from tangent",
            contact.tangency_error.to_degrees()
        );
        assert!(contact.gap <= tolerance, "gap {}", contact.gap);
    }
    for side in &filled.sides {
        assert!(side.angle.unwrap() <= tolerance, "{side:?}");
    }
}

#[test]
fn a_five_sided_hole_fills_and_sews_shut() {
    // A pentagon of straight edges whose corners rise and fall round it.
    let corners: Vec<Point> = (0..5)
        .map(|k| {
            let a = std::f64::consts::TAU * f64::from(k) / 5.0;
            Point::new(2.0 * a.cos(), 2.0 * a.sin(), 0.4 * (2.0 * a).sin())
        })
        .collect();
    let mut model = Model::new();
    let prism = open_prism(&mut model, &corners, -1.0);
    // Given out of order and with no supports: the loop is found anyway.
    let order = [3, 0, 4, 1, 2];
    let boundary: Vec<FillBoundary> = order.iter().map(|&k| g0(&prism.top[k], None)).collect();
    let tolerance = 1e-5;
    let filled = make_filling_n(&mut model, &boundary, &[], tolerance, T).unwrap();
    let face = filled.built.shape.clone();
    assert!(bounded_by(&model, &face, &prism.top));
    for (entry, side) in order.iter().zip(&filled.sides) {
        assert_eq!(
            side.edge.node(),
            prism.top[*entry].node(),
            "reports in given order"
        );
    }
    for edge in &prism.top {
        let gap = projected_gap(&model, &face, edge);
        assert!(gap <= tolerance, "an edge stands {gap} off the filling");
    }
    let mut faces = prism.walls.clone();
    faces.push(prism.floor.clone());
    faces.push(face);
    let solid = closed_solid(&mut model, &faces);
    let diagnosis = check(&model, &solid, T).unwrap();
    assert!(diagnosis.is_usable(), "{:?}", diagnosis.problems);
    assert!(volume_properties(&model, &solid, fine(), T).unwrap().mass > 0.0);
}

#[test]
fn one_closed_edge_fills_with_the_disk_it_bounds() {
    // A circle of radius 2 on the plane z = 1, as one closed edge: the fill
    // is that plane trimmed by the circle, and the face's area is the
    // disk's, 4π, not the patch rectangle's.
    let mut model = Model::new();
    let frame = Frame::new(Point::new(0.0, 0.0, 1.0), Direction::Z, Direction::X, T).unwrap();
    let circle = CircleCurve::new(Circle::new(frame, 2.0, T).unwrap());
    let edge = ogeom_algo::make_edge(&mut model, circle.into(), (0.0, std::f64::consts::TAU), T)
        .unwrap()
        .shape;
    let filled = make_filling_n(&mut model, &[g0(&edge, None)], &[], 1e-6, T).unwrap();
    let face = filled.built.shape;
    assert!(bounded_by(&model, &face, std::slice::from_ref(&edge)));
    let area = ogeom_algo::surface_properties(&model, &face, fine(), T)
        .unwrap()
        .mass;
    let disk = 4.0 * std::f64::consts::PI;
    assert!((area - disk).abs() < 1e-9 * disk, "area {area} against 4π");
}

/// A band of the cylinder of radius `r` round the z axis from angle `from`
/// through `sweep` turning about `axis`, between heights ±1, and the
/// straight edge it was swept from.
fn cylinder_band(
    model: &mut Model,
    r: f64,
    from: f64,
    sweep: f64,
    axis: Direction,
) -> (Shape, Shape) {
    let a = make_vertex(model, Point::new(r * from.cos(), r * from.sin(), -1.0)).shape;
    let b = make_vertex(model, Point::new(r * from.cos(), r * from.sin(), 1.0)).shape;
    let profile = line(model, &a, &b);
    let band = make_revolution(model, &profile, Axis::new(Point::ORIGIN, axis), sweep, T)
        .unwrap()
        .shape;
    (
        explore_unique(model, &band, ShapeType::Face).unwrap()[0].clone(),
        profile,
    )
}

#[test]
fn a_window_in_a_cylinder_fills_curvature_continuous_across_its_rulings() {
    // A window θ ∈ [-0.4, 0.4], z ∈ [-1, 1] in a cylinder of radius 5, with
    // bands of the cylinder either side. The rulings are G2 to the bands;
    // the arcs above and below are G0 with nothing beyond them. Across the
    // rulings the cylinder curves at 1/5, and the fill must too.
    let (r, half, sweep) = (5.0, 0.4, 0.5);
    let mut model = Model::new();
    let (left, left_ruling) = cylinder_band(&mut model, r, half, sweep, Direction::Z);
    let (right, right_ruling) = cylinder_band(&mut model, r, -half, sweep, -Direction::Z);
    let corner = |model: &Model, edge: &Shape, z: f64| -> Shape {
        let (a, b) = edge_vertices(model, edge).unwrap().unwrap();
        [a, b]
            .into_iter()
            .find(|v| (model.node(v).unwrap().data().as_vertex().unwrap().point.z - z).abs() < 1e-9)
            .unwrap()
    };
    let arc = |model: &mut Model, z: f64| -> Shape {
        let frame = Frame::new(Point::new(0.0, 0.0, z), Direction::Z, Direction::X, T).unwrap();
        let circle = CircleCurve::new(Circle::new(frame, r, T).unwrap());
        let from = corner(model, &right_ruling, z);
        let to = corner(model, &left_ruling, z);
        make_edge_between(model, circle.into(), (-half, half), &from, &to, T)
            .unwrap()
            .shape
    };
    let (low, high) = (arc(&mut model, -1.0), arc(&mut model, 1.0));
    let tolerance = 2e-3;
    let boundary = [
        FillBoundary {
            edge: left_ruling.clone(),
            support: Some(left.clone()),
            continuity: Continuity::G2,
        },
        g0(&high, None),
        FillBoundary {
            edge: right_ruling.clone(),
            support: Some(right.clone()),
            continuity: Continuity::G2,
        },
        g0(&low, None),
    ];
    let filled = make_filling_n(&mut model, &boundary, &[], tolerance, T).unwrap();
    let face = filled.built.shape.clone();
    let all = make_compound(&mut model, &[left, right, face.clone()])
        .unwrap()
        .shape;
    // Each band also holds its ruling placed at its far end; the contacts
    // that count are the rulings where they stand.
    let contacts: Vec<_> = ogeom_fillet::analyse_blend(&model, &all, &face, 200, T)
        .unwrap()
        .into_iter()
        .filter(|c| c.edge.location().is_identity())
        .collect();
    assert_eq!(contacts.len(), 2);
    for contact in &contacts {
        assert!(contact.tangency_error <= tolerance, "{contact:?}");
        assert!(contact.curvature_error <= tolerance, "{contact:?}");
    }
    // The curvature the fill matches is the cylinder's, 1/5, read off the
    // fill itself at a ruling's middle.
    let surface = surface_of(&model, &face);
    let p = Point::new(r * half.cos(), r * half.sin(), 0.0);
    let foot = project_on_surface(&surface, p, 32, T).unwrap();
    let (u, v) = foot.parameters;
    let curvature = surface.curvature_at(u, v, T).unwrap();
    let across = Vector::new(-half.sin(), half.cos(), 0.0);
    let k = curvature.normal_curvature(across).unwrap().abs();
    assert!(
        (k - 1.0 / r).abs() <= tolerance,
        "curvature {k} against {}",
        1.0 / r
    );
}

fn construction_error(result: Result<Filled, OgeomError>) -> String {
    match result {
        Err(OgeomError::Construction(message)) => message.to_string(),
        Err(other) => panic!("expected a construction refusal, got {other:?}"),
        Ok(_) => panic!("expected a refusal"),
    }
}

#[test]
fn refusals_name_what_is_wrong() {
    let mut model = Model::new();
    let prism = open_prism(&mut model, &saddle_corners(), -2.0);
    let sides: Vec<FillBoundary> = prism
        .top
        .iter()
        .zip(&prism.walls)
        .map(|(e, w)| g0(e, Some(w)))
        .collect();

    // Parametric continuity between two charts is not a filling's to meet.
    let mut c1 = sides.clone();
    c1[0].continuity = Continuity::C1;
    let message = construction_error(make_filling_n(&mut model, &c1, &[], 1e-5, T));
    assert!(message.contains("C1"), "{message}");

    // Tangency needs a face to be tangent to.
    let mut bare = sides.clone();
    bare[1].continuity = Continuity::G1;
    bare[1].support = None;
    let message = construction_error(make_filling_n(&mut model, &bare, &[], 1e-5, T));
    assert!(message.contains("no support"), "{message}");

    // Three sides of four do not close.
    let message = construction_error(make_filling_n(&mut model, &sides[..3], &[], 1e-5, T));
    assert!(message.contains("does not close"), "{message}");

    // The same edge twice.
    let mut twice = sides.clone();
    twice[3] = twice[0].clone();
    let message = construction_error(make_filling_n(&mut model, &twice, &[], 1e-5, T));
    assert!(message.contains("again"), "{message}");

    // A support that does not hold the edge.
    let mut wrong = sides.clone();
    wrong[0].support = Some(prism.walls[2].clone());
    let message = construction_error(make_filling_n(&mut model, &wrong, &[], 1e-5, T));
    assert!(message.contains("does not hold"), "{message}");

    // Tangent to a vertical wall along a level hole: the fill would stand
    // square to the plane the hole spans, which no height field does.
    let mut steep = sides.clone();
    for side in &mut steep {
        side.continuity = Continuity::G1;
    }
    let message = construction_error(make_filling_n(&mut model, &steep, &[], 1e-5, T));
    assert!(message.contains("square to the plane"), "{message}");

    // A point outside the hole.
    let outside = make_vertex(&mut model, Point::new(3.0, 0.0, 0.0)).shape;
    let message = construction_error(make_filling_n(&mut model, &sides, &[outside], 1e-5, T));
    assert!(message.contains("outside the hole"), "{message}");

    // A face is not a constraint.
    let wall = prism.walls[0].clone();
    let message = construction_error(make_filling_n(&mut model, &sides, &[wall], 1e-5, T));
    assert!(message.contains("vertices and edges"), "{message}");

    // A tolerance that is not a distance.
    let message = construction_error(make_filling_n(&mut model, &sides, &[], 0.0, T));
    assert!(message.contains("positive and finite"), "{message}");
}

#[test]
fn a_loop_that_folds_over_itself_is_refused() {
    // A bow-tie of unequal lobes in the plane z = 0: the loop spans that
    // plane and crosses itself in it.
    let mut model = Model::new();
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(3.0, 1.0, 0.0),
        Point::new(3.0, 0.0, 0.0),
        Point::new(0.0, 2.0, 0.0),
    ];
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|c| make_vertex(&mut model, *c).shape)
        .collect();
    let edges: Vec<Shape> = (0..4)
        .map(|k| line(&mut model, &vertices[k], &vertices[(k + 1) % 4]))
        .collect();
    let boundary: Vec<FillBoundary> = edges.iter().map(|e| g0(e, None)).collect();
    let message = construction_error(make_filling_n(&mut model, &boundary, &[], 1e-5, T));
    assert!(message.contains("crosses itself"), "{message}");
}

#[test]
fn a_tolerance_the_finest_net_cannot_meet_is_refused_with_what_it_reached() {
    let mut model = Model::new();
    let (cylinders, arcs) = rounded_corner(&mut model, 2.0, 3.0);
    let boundary: Vec<FillBoundary> = arcs
        .iter()
        .zip(&cylinders)
        .map(|(edge, support)| FillBoundary {
            edge: edge.clone(),
            support: Some(support.clone()),
            continuity: Continuity::G2,
        })
        .collect();
    match make_filling_n(&mut model, &boundary, &[], 1e-12, T) {
        Err(OgeomError::NotDone(message)) => {
            let message = message.to_string();
            assert!(message.contains("finest net"), "{message}");
            assert!(message.contains("side"), "{message}");
        }
        other => panic!("expected a refusal, got {other:?}"),
    }
}
