//! Filling boundary edges with a fitted patch: the doubly ruled saddle is
//! the case with an exact answer, and the fit must land on it. The
//! N-sided filling also takes placed edges (a prism's far end) and sides
//! of separate sheets that meet only where their ends coincide, and its
//! face sews to every support.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::{LineCurve, Surface as _, SurfaceGeometry};
use ogeom::math::Point;
use ogeom::topo::{Model, NodeData, Shape};

const T: Tolerances = Tolerances::millimetres();

#[test]
fn the_saddles_boundary_fills_to_the_saddle() {
    // z = x·y over the unit square: all four boundary edges are straight,
    // and the Coons blend of straight boundaries is exactly the bilinear
    // saddle. The fit therefore has an exact target to hit.
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(1.0, 0.0, 0.0),
        Point::new(1.0, 1.0, 1.0),
        Point::new(0.0, 1.0, 0.0),
    ];
    let mut model = Model::new();
    let vertices: Vec<Shape> = corners
        .iter()
        .map(|c| ogeom::algo::make_vertex(&mut model, *c).shape)
        .collect();
    let edges: Vec<Shape> = (0..4)
        .map(|i| {
            let (a, b) = (corners[i], corners[(i + 1) % 4]);
            ogeom::algo::make_edge_between(
                &mut model,
                LineCurve::segment(a, b, T).unwrap().into(),
                (0.0, a.distance(b)),
                &vertices[i],
                &vertices[(i + 1) % 4],
                T,
            )
            .unwrap()
            .shape
        })
        .collect();

    let filled = ogeom::offset::make_filling(
        &mut model,
        &[
            edges[0].clone(),
            edges[1].clone(),
            edges[2].clone(),
            edges[3].clone(),
        ],
        12,
        1e-6,
        T,
    )
    .unwrap();

    // The face's surface is the saddle: probe the interior against z = x·y.
    let surface_id = match model.node(&filled.shape).unwrap().data() {
        NodeData::Face(data) => data.surface,
        _ => panic!("the filling is a face"),
    };
    let surface = model.geometry().surface(surface_id).unwrap();
    let SurfaceGeometry::BSpline(patch) = surface else {
        panic!("the filling is a fitted patch");
    };
    let ((ua, ub), (va, vb)) = patch.domain();
    for (fu, fv) in [(0.5, 0.5), (0.25, 0.75), (0.9, 0.1)] {
        let p = patch
            .point_at(ua + (ub - ua) * fu, va + (vb - va) * fv, T)
            .unwrap();
        assert!(
            (p.z - p.x * p.y).abs() < 1e-6,
            "the patch is the saddle at ({fu}, {fv}): {p:?}"
        );
    }

    // History names every boundary edge as modified into the face.
    for edge in &edges {
        assert!(
            !filled.history.modified(edge).is_empty(),
            "the filling records its boundary"
        );
    }
}

/// The edges of `shape` whose bounds lie wholly at height `z`.
fn edges_at_height(model: &Model, shape: &Shape, z: f64) -> Vec<Shape> {
    ogeom::topo::explore_unique(model, shape, ogeom::topo::ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let b = ogeom::algo::shape_bounds(model, e, T).unwrap();
            (b.low().unwrap().z - z).abs() < 1e-6 && (b.high().unwrap().z - z).abs() < 1e-6
        })
        .collect()
}

const SQUARE: [(f64, f64); 4] = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)];

#[test]
fn an_extruded_squares_walls_bake_into_a_sheet_in_place() {
    use ogeom::math::Vector;
    use ogeom::topo::{ShapeType, explore, explore_unique};

    let mut model = Model::new();
    let points = SQUARE.map(|(x, y)| Point::new(x, y, 0.0));
    let wire = ogeom::algo::make_polygon(&mut model, &points, true, T)
        .unwrap()
        .shape;
    let walls = ogeom::algo::make_prism(&mut model, &wire, Vector::new(0.0, 0.0, 5.0), T)
        .unwrap()
        .shape;
    let baked = ogeom::algo::baked_shape(&mut model, &walls, T).unwrap();
    let sheet = baked.shape.clone();
    assert_eq!(model.kind_of(&sheet).unwrap(), ShapeType::Shell);
    assert_eq!(
        explore_unique(&model, &sheet, ShapeType::Face)
            .unwrap()
            .len(),
        4
    );
    for edge in explore(&model, &sheet, ogeom::topo::Filter::OfType(ShapeType::Edge)).unwrap() {
        assert!(edge.location().is_identity(), "every edge stands in place");
    }
    let diagnosis = ogeom::algo::check(&model, &sheet, T).unwrap();
    assert!(diagnosis.is_usable(), "{:?}", diagnosis.problems);
    let area =
        ogeom::algo::surface_properties(&model, &sheet, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    assert!((area - 200.0).abs() < 1e-9, "four 10 by 5 walls: {area}");

    // A wall's bottom and top edges are one node at two placements; each
    // has its own twin, where it stands.
    let wall = explore_unique(&model, &walls, ShapeType::Face).unwrap()[0].clone();
    assert!(baked.history.modified(&wall)[0].node() != wall.node());
    let (bottom, top) = (
        edges_at_height(&model, &wall, 0.0).remove(0),
        edges_at_height(&model, &wall, 5.0).remove(0),
    );
    assert_eq!(bottom.node(), top.node());
    let (low, high) = (
        baked.history.modified(&bottom)[0].clone(),
        baked.history.modified(&top)[0].clone(),
    );
    assert!(low.node() != high.node());
    assert_eq!(edges_at_height(&model, &low, 0.0).len(), 1);
    assert_eq!(edges_at_height(&model, &high, 5.0).len(), 1);
}

#[test]
fn an_extruded_squares_far_edges_fill_and_close_the_box() {
    use ogeom::geom::{Continuity, PlaneSurface};
    use ogeom::math::{Frame, Plane, Vector};
    use ogeom::offset::{FillBoundary, make_filling_n};

    // Walls: a square wire extruded 5 up, whose top edges are the wire's
    // own edges under the extrusion's translation; a floor on the wire.
    let mut model = Model::new();
    let points = SQUARE.map(|(x, y)| Point::new(x, y, 0.0));
    let wire = ogeom::algo::make_polygon(&mut model, &points, true, T)
        .unwrap()
        .shape;
    let walls = ogeom::algo::make_prism(&mut model, &wire, Vector::new(0.0, 0.0, 5.0), T)
        .unwrap()
        .shape;
    let floor = ogeom::algo::make_face(
        &mut model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let wall_faces =
        ogeom::topo::explore_unique(&model, &walls, ogeom::topo::ShapeType::Face).unwrap();
    let mut faces = wall_faces.clone();
    faces.push(floor);
    let open = ogeom::algo::sew(&mut model, &faces, T).unwrap();
    assert_eq!(open.shells.len(), 1);
    let open = open.shells[0].clone();
    let faces = ogeom::topo::explore_unique(&model, &open, ogeom::topo::ShapeType::Face).unwrap();

    // Each wall's top edge, with that wall as its support.
    let mut sides = Vec::new();
    for face in &wall_faces {
        for edge in edges_at_height(&model, face, 5.0) {
            assert!(!edge.location().is_identity(), "a far edge is placed");
            sides.push(FillBoundary {
                edge,
                support: Some(face.clone()),
                continuity: Continuity::C0,
            });
        }
    }
    assert_eq!(sides.len(), 4, "four walls, one top edge each");
    let filled = make_filling_n(&mut model, &sides, &[], 1e-3, T).unwrap();
    for side in &filled.sides {
        assert!(side.gap <= 1e-3, "{side:?}");
    }

    let mut all = faces.clone();
    all.push(filled.built.shape.clone());
    let sewn = ogeom::algo::sew(&mut model, &all, T).unwrap();
    assert_eq!(sewn.shells.len(), 1, "one shell");
    assert!(
        sewn.free_edges.is_empty(),
        "the shell closes: {} free edges",
        sewn.free_edges.len()
    );
    let solid = ogeom::algo::make_solid(&mut model, &sewn.shells)
        .unwrap()
        .shape;
    let volume =
        ogeom::algo::volume_properties(&model, &solid, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    assert!(
        (volume - 500.0).abs() < 1e-6,
        "volume {volume} against 10·10·5"
    );
}

/// One strip: `edge` extruded along `travel`, offered as a side meeting it
/// at `continuity`.
fn strip_side(
    model: &mut Model,
    edge: Shape,
    travel: ogeom::math::Vector,
    continuity: ogeom::geom::Continuity,
) -> (Shape, ogeom::offset::FillBoundary) {
    let strip = ogeom::algo::make_prism(model, &edge, travel, T)
        .unwrap()
        .shape;
    let side = ogeom::offset::FillBoundary {
        edge,
        support: Some(strip.clone()),
        continuity,
    };
    (strip, side)
}

/// Sew the strips and the filling: the shell they make and the filling's
/// face in it.
fn sewn_with(model: &mut Model, strips: &[Shape], filling: &Shape) -> (Shape, Shape) {
    let mut all = strips.to_vec();
    all.push(filling.clone());
    let sewn = ogeom::algo::sew(model, &all, T).unwrap();
    assert_eq!(
        sewn.shells.len(),
        1,
        "the strips and the filling make one shell"
    );
    let face = sewn
        .history
        .modified(filling)
        .first()
        .cloned()
        .unwrap_or_else(|| filling.clone());
    (sewn.shells[0].clone(), face)
}

#[test]
fn strips_meeting_only_at_their_ends_chain_into_one_filling() {
    use ogeom::geom::Continuity;
    use ogeom::math::Vector;
    use ogeom::offset::make_filling_n;

    // Four separate strips, each the line along one side of a square at
    // z = 5 extruded outward and down: neighbouring lines end at the same
    // corner points, each on vertices of its own.
    let mut model = Model::new();
    let corners = SQUARE.map(|(x, y)| Point::new(x, y, 5.0));
    let (mut strips, mut sides) = (Vec::new(), Vec::new());
    for k in 0..4 {
        let (a, b) = (corners[k], corners[(k + 1) % 4]);
        let line = ogeom::algo::make_polygon(&mut model, &[a, b], false, T)
            .unwrap()
            .shape;
        let edge = ogeom::topo::explore_unique(&model, &line, ogeom::topo::ShapeType::Edge)
            .unwrap()
            .remove(0);
        // The square runs counter-clockwise seen from above, so its outside
        // is on the right of each side.
        let along = (b - a) * 0.1;
        let travel = Vector::new(along.y, -along.x, 0.0) * 5.0 + Vector::new(0.0, 0.0, -2.0);
        let (strip, side) = strip_side(&mut model, edge, travel, Continuity::C0);
        strips.push(strip);
        sides.push(side);
    }
    let filled = make_filling_n(&mut model, &sides, &[], 1e-6, T).unwrap();
    for (side, report) in sides.iter().zip(&filled.sides) {
        assert!(report.gap <= 1e-6, "{report:?}");
        assert!(
            !filled.built.history.generated(&side.edge).is_empty(),
            "the history names what each side became"
        );
    }
    let (shell, face) = sewn_with(&mut model, &strips, &filled.built.shape);
    let contacts = ogeom::fillet::analyse_blend(&model, &shell, &face, 21, T).unwrap();
    assert_eq!(
        contacts.len(),
        4,
        "the filling shares an edge with each strip"
    );
    for contact in &contacts {
        assert!(contact.gap <= 1e-6, "{contact:?}");
    }
}

/// What a dome filling achieved: the height where it crosses the axis,
/// the largest distance from the sphere over the hole, and the worst
/// tangency and gap against the strips.
struct Dome {
    top: f64,
    off_sphere: f64,
    tangency: f64,
    gap: f64,
}

/// A square on the sphere of radius `r` round the origin, filled at G1
/// from four strips: four arcs of great circles, each in a plane through
/// the origin leaning out by `lean`, each extruded along its plane's
/// outward and downward normal. That normal is square to the sphere's
/// radius all along the arc, so each strip is tangent to the sphere along
/// its arc and the strips share the sphere's tangent plane at every
/// corner: the spherical cap is the filling every condition agrees with.
/// Separate prisms: neighbouring arcs end on vertices of their own.
fn dome(r: f64, lean: f64, tolerance: f64) -> Dome {
    use ogeom::geom::{CircleCurve, Continuity};
    use ogeom::math::{Circle, Direction, Frame, Vector};
    use ogeom::offset::make_filling_n;

    let z0 = r / (2.0 * lean * lean + 1.0).sqrt();
    let corners = [(1.0, -1.0), (1.0, 1.0), (-1.0, 1.0), (-1.0, -1.0)]
        .map(|(sx, sy)| Point::new(sx * lean * z0, sy * lean * z0, z0));
    let mut model = Model::new();
    let (mut strips, mut sides) = (Vec::new(), Vec::new());
    for k in 0..4 {
        let (a, b) = (corners[k], corners[(k + 1) % 4]);
        let (ua, ub) = (a.to_vector() * (1.0 / r), b.to_vector() * (1.0 / r));
        let normal = Direction::new(ua.cross(ub), T).unwrap();
        let frame = Frame::new(Point::ORIGIN, normal, Direction::new(ua, T).unwrap(), T).unwrap();
        let sweep = ua.dot(ub).acos();
        let circle = CircleCurve::new(Circle::new(frame, r, T).unwrap());
        let edge = ogeom::algo::make_edge(&mut model, circle.into(), (0.0, sweep), T)
            .unwrap()
            .shape;
        // The plane's normal turned away from the square's middle, and so
        // downward.
        let mut out = normal.vector();
        let middle = (ua + ub) * 0.5;
        if out.dot(Vector::new(middle.x, middle.y, 0.0)) < 0.0 {
            out *= -1.0;
        }
        assert!(out.z < 0.0, "the strip runs out and down");
        let (strip, side) = strip_side(&mut model, edge, out * (0.3 * r), Continuity::G1);
        strips.push(strip);
        sides.push(side);
    }
    let filled = make_filling_n(&mut model, &sides, &[], tolerance, T).unwrap();
    let face = filled.built.shape.clone();

    let surface_id = match model.node(&face).unwrap().data() {
        NodeData::Face(data) => data.surface,
        _ => panic!("a face"),
    };
    let surface = model.geometry().surface(surface_id).unwrap().clone();
    let top = ogeom::algo::project_on_surface(&surface, Point::new(0.0, 0.0, r), 32, T)
        .unwrap()
        .point
        .z;
    // Points of the sphere over the square the corners span, which the
    // hole covers, each projected onto the filling.
    let mut off_sphere = 0.0f64;
    let half = lean * z0;
    for i in 0..=8 {
        for j in 0..=8 {
            let (x, y) = (
                half * (f64::from(i) / 4.0 - 1.0),
                half * (f64::from(j) / 4.0 - 1.0),
            );
            let p = Point::new(x, y, (r * r - x * x - y * y).sqrt());
            let foot = ogeom::algo::project_on_surface(&surface, p, 32, T).unwrap();
            off_sphere = off_sphere.max(foot.distance);
        }
    }

    let (shell, face) = sewn_with(&mut model, &strips, &face);
    let contacts = ogeom::fillet::analyse_blend(&model, &shell, &face, 41, T).unwrap();
    assert_eq!(
        contacts.len(),
        4,
        "the filling shares an edge with each strip"
    );
    let (mut tangency, mut gap) = (0.0f64, 0.0f64);
    for contact in &contacts {
        tangency = tangency.max(contact.tangency_error);
        gap = gap.max(contact.gap);
    }
    Dome {
        top,
        off_sphere,
        tangency,
        gap,
    }
}

/// The crown of the thin-plate filling of a round hole on the sphere of
/// radius `r`, the hole's rim a circle of radius `rho` round the axis, the
/// filling taking the sphere's height and slope all round it.
///
/// The fill is the height field of least `∫∫ h_xx² + 2 h_xy² + h_yy²`
/// with those values; with round data that is the paraboloid `a + b·ρ²`,
/// the one regular radial biharmonic, whose slope `2·b·rho` is the
/// sphere's `-rho / z` at the rim's height `z`. It crowns at
/// `z + rho² / (2·z)`, above the sphere's `r`.
fn paraboloid_crown(r: f64, rho: f64) -> f64 {
    let z = r.mul_add(r, -rho * rho).sqrt();
    z + rho * rho / (2.0 * z)
}

#[test]
fn strips_tangent_to_a_dome_fill_tangent_to_each() {
    // The square hole's crown lies between those of the round holes
    // through its arcs' middles and through its corners; a crater, or a
    // middle left to the fit's approximation error, falls outside by far
    // more than `slack`.
    let r = 10.0f64;
    let tolerance = 0.1f64.to_radians();
    let slack = 0.005 * r;
    for lean in [0.2, 0.3, 0.5, 0.8] {
        let d = dome(r, lean, tolerance);
        let z0 = r / (2.0 * lean * lean + 1.0).sqrt();
        let rim = r / lean.mul_add(lean, 1.0).sqrt();
        let low = paraboloid_crown(r, lean * rim);
        let high = paraboloid_crown(r, 2.0f64.sqrt() * lean * z0);
        assert!(
            d.top > low - slack && d.top < high + slack,
            "lean {lean}: the crown at {} against {low} to {high}",
            d.top
        );
        assert!(
            d.off_sphere < high - r + slack,
            "lean {lean}: {} off the sphere",
            d.off_sphere
        );
        assert!(
            d.tangency <= tolerance,
            "lean {lean}: {} degrees from tangent",
            d.tangency.to_degrees()
        );
        assert!(d.gap <= tolerance, "lean {lean}: a gap of {}", d.gap);
    }
}

/// An arc of the circle of radius `r` round `centre` in the plane square to
/// `normal`, starting along `start` and turning `sweep` about `normal`.
fn arc_edge(
    model: &mut Model,
    centre: Point,
    normal: ogeom::math::Vector,
    start: ogeom::math::Vector,
    r: f64,
    sweep: f64,
) -> Shape {
    use ogeom::geom::CircleCurve;
    use ogeom::math::{Circle, Direction, Frame};
    let frame = Frame::new(
        centre,
        Direction::new(normal, T).unwrap(),
        Direction::new(start, T).unwrap(),
        T,
    )
    .unwrap();
    let circle = CircleCurve::new(Circle::new(frame, r, T).unwrap());
    ogeom::algo::make_edge(model, circle.into(), (0.0, sweep), T)
        .unwrap()
        .shape
}

/// A line edge from `a` to `b` on vertices of its own.
fn line_edge(model: &mut Model, a: Point, b: Point) -> Shape {
    let line = ogeom::algo::make_polygon(model, &[a, b], false, T)
        .unwrap()
        .shape;
    ogeom::topo::explore_unique(model, &line, ogeom::topo::ShapeType::Edge)
        .unwrap()
        .remove(0)
}

/// The surface a face lies on.
fn surface_of(model: &Model, face: &Shape) -> SurfaceGeometry {
    match model.node(face).unwrap().data() {
        NodeData::Face(data) => model.geometry().surface(data.surface).unwrap().clone(),
        _ => panic!("a face"),
    }
}

/// The largest distance from `points` to the surface.
fn worst_off(surface: &SurfaceGeometry, points: impl IntoIterator<Item = Point>) -> f64 {
    points
        .into_iter()
        .map(|p| {
            ogeom::algo::project_on_surface(surface, p, 32, T)
                .unwrap()
                .distance
        })
        .fold(0.0, f64::max)
}

/// The mesh of a face that `check` finds valid.
fn valid_mesh(model: &Model, face: &Shape) -> Vec<Point> {
    let diagnosis = ogeom::algo::check(model, face, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let mesh =
        ogeom::mesh::triangulate_face(model, face, ogeom::mesh::Deflection::default(), T).unwrap();
    assert!(mesh.positions.len() > 10);
    mesh.positions
}

/// How far the mesh, reflected by `mirror`, stands off the surface: a loop
/// that maps onto itself under the reflection is filled by a face that
/// does too.
fn mirrored_off(surface: &SurfaceGeometry, mesh: &[Point], mirror: impl Fn(Point) -> Point) -> f64 {
    worst_off(surface, mesh.iter().map(|p| mirror(*p)))
}

#[test]
fn a_four_line_saddle_fills_over_the_plane_it_projects_simply_on() {
    use ogeom::geom::Continuity;
    use ogeom::offset::{FillBoundary, make_filling_n};

    // A diamond rising at one corner and running out at another; along
    // (0, 1, -1) it is a rhombus.
    let p = [
        Point::new(-10.0, 0.0, 0.0),
        Point::new(0.0, 0.0, 10.0),
        Point::new(10.0, 0.0, 0.0),
        Point::new(0.0, -10.0, 0.0),
    ];
    let mut model = Model::new();
    let sides: Vec<FillBoundary> = (0..4)
        .map(|i| FillBoundary {
            edge: line_edge(&mut model, p[i], p[(i + 1) % 4]),
            support: None,
            continuity: Continuity::C0,
        })
        .collect();
    let tolerance = 1e-3;
    let filled = make_filling_n(&mut model, &sides, &[], tolerance, T).unwrap();
    let face = filled.built.shape.clone();
    let surface = surface_of(&model, &face);
    let lines = (0..4).flat_map(|i| {
        (0..=50).map(move |k| {
            let f = f64::from(k) / 50.0;
            p[i] + (p[(i + 1) % 4] - p[i]) * f
        })
    });
    let off = worst_off(&surface, lines);
    assert!(off <= tolerance, "a side stands {off} off the filling");
    let mesh = valid_mesh(&model, &face);
    let off = mirrored_off(&surface, &mesh, |q| Point::new(-q.x, q.y, q.z));
    assert!(
        off <= tolerance,
        "the mirrored mesh stands {off} off the face"
    );
    // The ruled saddle the four lines bound, against which a fair filling
    // stands within a small share of the hole's 20 units.
    let bilinear = (1..10).flat_map(|i| {
        (1..10).map(move |j| {
            let (s, t) = (f64::from(i) / 10.0, f64::from(j) / 10.0);
            Point::from_vector(
                p[0].to_vector() * ((1.0 - s) * (1.0 - t))
                    + p[1].to_vector() * (s * (1.0 - t))
                    + p[2].to_vector() * (s * t)
                    + p[3].to_vector() * ((1.0 - s) * t),
            )
        })
    });
    let off = worst_off(&surface, bilinear);
    assert!(off <= 0.1, "the filling strays {off} from the ruled saddle");
}

#[test]
fn two_semicircles_in_crossing_planes_fill_within_tolerance() {
    use ogeom::geom::Continuity;
    use ogeom::math::Vector;
    use ogeom::offset::{FillBoundary, make_filling_n};

    // One semicircle rises over the x axis in the XZ plane, the other dips
    // below it in the XY plane: they meet at a right angle at both ends,
    // and seen along the loop's own plane the two corners are smooth.
    let r = 10.0;
    let pi = std::f64::consts::PI;
    let mut model = Model::new();
    let up = arc_edge(
        &mut model,
        Point::ORIGIN,
        Vector::new(0.0, 1.0, 0.0),
        Vector::new(-1.0, 0.0, 0.0),
        r,
        pi,
    );
    let down = arc_edge(
        &mut model,
        Point::ORIGIN,
        Vector::new(0.0, 0.0, -1.0),
        Vector::new(1.0, 0.0, 0.0),
        r,
        pi,
    );
    let sides = [up, down].map(|edge| FillBoundary {
        edge,
        support: None,
        continuity: Continuity::C0,
    });
    let tolerance = 1e-3;
    let filled = make_filling_n(&mut model, &sides, &[], tolerance, T).unwrap();
    let face = filled.built.shape.clone();
    let surface = surface_of(&model, &face);
    let arcs = (0..=200).flat_map(|k| {
        let a = pi * f64::from(k) / 200.0;
        [
            Point::new(-r * a.cos(), 0.0, r * a.sin()),
            Point::new(r * a.cos(), -r * a.sin(), 0.0),
        ]
    });
    let off = worst_off(&surface, arcs);
    assert!(off <= tolerance, "an arc stands {off} off the filling");

    // The loop maps onto itself reflected across x = 0, and swapped arc
    // for arc by (x, y, z) -> (x, -z, -y); it bounds the quarter of the
    // slab -10 <= x <= 10 with y <= 0 <= z, which is convex, so a fair
    // filling stays in it.
    let mesh = valid_mesh(&model, &face);
    for (name, mirror) in [
        (
            "across x = 0",
            (|p: Point| Point::new(-p.x, p.y, p.z)) as fn(Point) -> Point,
        ),
        ("arc for arc", |p: Point| Point::new(p.x, -p.z, -p.y)),
    ] {
        let off = mirrored_off(&surface, &mesh, mirror);
        assert!(
            off <= tolerance,
            "reflected {name}, the mesh stands {off} off the face"
        );
    }
    for p in &mesh {
        assert!(
            p.x.abs() <= r + tolerance && p.y <= tolerance && p.z >= -tolerance,
            "a mesh point outside the quarter slab: {p:?}"
        );
    }
}

/// The tube's radius and height.
const TUBE: (f64, f64) = (5.0, 10.0);

/// A tube standing on a circle at z = 0, as two half arcs or one closed
/// edge, swept up: its wall faces, and its top rim as sides meeting them
/// at `continuity`.
fn tube_rim(
    model: &mut Model,
    halves: bool,
    continuity: ogeom::geom::Continuity,
) -> (Vec<Shape>, Vec<ogeom::offset::FillBoundary>) {
    use ogeom::geom::CircleCurve;
    use ogeom::math::{Circle, Frame, Vector};
    use ogeom::offset::FillBoundary;
    use std::f64::consts::{PI, TAU};

    let (r, h) = TUBE;
    let circle: ogeom::geom::Curve =
        CircleCurve::new(Circle::new(Frame::WORLD, r, T).unwrap()).into();
    let edges = if halves {
        let ends = [Point::new(r, 0.0, 0.0), Point::new(-r, 0.0, 0.0)]
            .map(|p| ogeom::algo::make_vertex(model, p).shape);
        [(0.0, PI, 0, 1), (PI, TAU, 1, 0)]
            .map(|(a, b, from, to)| {
                ogeom::algo::make_edge_between(
                    model,
                    circle.clone(),
                    (a, b),
                    &ends[from],
                    &ends[to],
                    T,
                )
                .unwrap()
                .shape
            })
            .to_vec()
    } else {
        vec![
            ogeom::algo::make_edge(model, circle, (0.0, TAU), T)
                .unwrap()
                .shape,
        ]
    };
    let wire = ogeom::algo::make_wire(model, &edges, T).unwrap().shape;
    let walls = ogeom::algo::make_prism(model, &wire, Vector::new(0.0, 0.0, h), T)
        .unwrap()
        .shape;
    let wall_faces =
        ogeom::topo::explore_unique(model, &walls, ogeom::topo::ShapeType::Face).unwrap();
    let mut sides = Vec::new();
    for face in &wall_faces {
        for edge in edges_at_height(model, face, h) {
            sides.push(FillBoundary {
                edge,
                support: Some(face.clone()),
                continuity,
            });
        }
    }
    assert_eq!(sides.len(), edges.len(), "one rim edge per wall edge");
    (wall_faces, sides)
}

/// What a cap on a tube achieved.
struct Cap {
    /// Where the cap crosses the tube's axis.
    top: f64,
    /// The largest distance from the rim circle to the cap's surface.
    rim_off: f64,
    /// The largest distance from the cap's mesh, turned a quarter and an
    /// eighth of a turn about the axis, to its surface.
    turned_off: f64,
    /// The worst tangency, gap and curvature difference against the wall.
    tangency: f64,
    gap: f64,
    curvature: f64,
}

/// The tube [`tube_rim`] builds, its top rim filled meeting the wall at
/// `continuity`; the cap sewn to the wall leaves only the bottom rim open,
/// and `check` finds nothing else wrong with the sewn sheet.
fn capped_tube(halves: bool, continuity: ogeom::geom::Continuity, tolerance: f64) -> Cap {
    use ogeom::offset::make_filling_n;
    use std::f64::consts::PI;

    let (r, h) = TUBE;
    let mut model = Model::new();
    let (wall_faces, sides) = tube_rim(&mut model, halves, continuity);
    let filled = make_filling_n(&mut model, &sides, &[], tolerance, T).unwrap();
    let face = filled.built.shape.clone();
    let surface = surface_of(&model, &face);

    let top = ogeom::algo::project_on_surface(&surface, Point::new(0.0, 0.0, h + r), 32, T)
        .unwrap()
        .point;
    assert!(
        top.x.hypot(top.y) < 1e-3,
        "the cap crowns on the axis: {top:?}"
    );
    let rim_off = worst_off(
        &surface,
        (0..360).map(|k| {
            let a = f64::from(k).to_radians();
            Point::new(r * a.cos(), r * a.sin(), h)
        }),
    );
    let mesh = valid_mesh(&model, &face);
    let turned_off = [PI / 2.0, PI / 4.0]
        .map(|a| {
            mirrored_off(&surface, &mesh, |p| {
                Point::new(
                    p.x * a.cos() - p.y * a.sin(),
                    p.x * a.sin() + p.y * a.cos(),
                    p.z,
                )
            })
        })
        .into_iter()
        .fold(0.0, f64::max);

    let mut all = wall_faces.clone();
    all.push(face.clone());
    let sewn = ogeom::algo::sew(&mut model, &all, T).unwrap();
    assert_eq!(sewn.shells.len(), 1, "the wall and the cap make one shell");
    assert_eq!(
        sewn.free_edges.len(),
        sides.len(),
        "only the bottom rim stays open"
    );
    let shell = sewn.shells[0].clone();
    let face = sewn
        .history
        .modified(&face)
        .first()
        .cloned()
        .unwrap_or(face);
    // The open bottom rim is the shell's one finding.
    let diagnosis = ogeom::algo::check(&model, &shell, T).unwrap();
    assert!(
        diagnosis.is_usable()
            && diagnosis
                .problems
                .iter()
                .all(|p| p.kind == ogeom::topo::ShapeType::Shell),
        "{diagnosis}"
    );
    let contacts = ogeom::fillet::analyse_blend(&model, &shell, &face, 41, T).unwrap();
    assert_eq!(contacts.len(), sides.len(), "the cap meets every wall");
    let (mut tangency, mut gap, mut curvature) = (0.0f64, 0.0f64, 0.0f64);
    for contact in &contacts {
        tangency = tangency.max(contact.tangency_error);
        gap = gap.max(contact.gap);
        curvature = curvature.max(contact.curvature_error);
    }
    Cap {
        top: top.z,
        rim_off,
        turned_off,
        tangency,
        gap,
        curvature,
    }
}

#[test]
fn a_tangent_cap_closes_a_tube() {
    use ogeom::geom::Continuity;

    // The wall stands square to the rim's plane all round, so the cap
    // leaves the rim straight up and crowns above it; the tube is round,
    // so the cap is too. A G2 cap leaves the rim unbent, as the wall
    // stands straight.
    let tolerance = 1e-3;
    for (halves, continuity) in [
        (true, Continuity::G1),
        (false, Continuity::G1),
        (true, Continuity::G2),
        (false, Continuity::G2),
    ] {
        let cap = capped_tube(halves, continuity, tolerance);
        let case = format!("{continuity:?}, halves {halves}");
        assert!(cap.top > 10.5, "{case}: the cap crowns at {}", cap.top);
        assert!(
            cap.rim_off <= tolerance,
            "{case}: the rim stands {} off",
            cap.rim_off
        );
        // Round to within a fiftieth of its 2.5 crown: the net the patch
        // is fitted on is square, the cap round.
        assert!(
            cap.turned_off <= 0.05,
            "{case}: turned, {} off",
            cap.turned_off
        );
        assert!(cap.gap <= tolerance, "{case}: a gap of {}", cap.gap);
        assert!(
            cap.tangency <= tolerance,
            "{case}: {} radians from tangent",
            cap.tangency
        );
        if continuity == Continuity::G2 {
            assert!(
                cap.curvature <= tolerance,
                "{case}: the curvature differs by {}",
                cap.curvature
            );
        }
    }
}

#[test]
fn a_tangent_cap_through_a_point_is_refused_naming_the_walls_angle() {
    use ogeom::geom::Continuity;
    use ogeom::offset::make_filling_n;

    // The cap's chart is drawn from its rim, which places no interior
    // point; the refusal says how steep the wall stands and how steep a
    // height over the rim's plane may be.
    let mut model = Model::new();
    let (_, sides) = tube_rim(&mut model, true, Continuity::G1);
    let crown = ogeom::algo::make_vertex(&mut model, Point::new(0.0, 0.0, 12.0)).shape;
    match make_filling_n(&mut model, &sides, &[crown], 1e-3, T) {
        Err(ogeom::core::OgeomError::Construction(message)) => {
            let message = message.to_string();
            assert!(message.contains("at 90.0 degrees"), "{message}");
            assert!(message.contains("past the 84.3"), "{message}");
        }
        other => panic!("expected a refusal, got {other:?}"),
    }
}
