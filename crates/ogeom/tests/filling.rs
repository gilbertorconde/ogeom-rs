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

#[test]
fn strips_tangent_to_a_dome_fill_tangent_to_each() {
    use ogeom::geom::{CircleCurve, Continuity};
    use ogeom::math::{Circle, Direction, Frame, Vector};
    use ogeom::offset::make_filling_n;

    // A square on the sphere of radius 10 round the origin: four arcs of
    // great circles, each in a plane through the origin leaning out by
    // `lean`, each extruded along its plane's outward and downward normal.
    // That normal is square to the sphere's radius all along the arc, so
    // each strip is tangent to the sphere along its arc and the strips
    // share the sphere's tangent plane at every corner. Separate prisms:
    // neighbouring arcs end on vertices of their own.
    let (r, lean) = (10.0f64, 0.3f64);
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
        let (strip, side) = strip_side(&mut model, edge, out * 3.0, Continuity::G1);
        strips.push(strip);
        sides.push(side);
    }
    let tolerance = 0.1f64.to_radians();
    let filled = make_filling_n(&mut model, &sides, &[], tolerance, T).unwrap();
    let face = filled.built.shape.clone();

    // It crowns above the square, on the sphere.
    let surface_id = match model.node(&face).unwrap().data() {
        NodeData::Face(data) => data.surface,
        _ => panic!("a face"),
    };
    let surface = model.geometry().surface(surface_id).unwrap().clone();
    let top = ogeom::algo::project_on_surface(&surface, Point::new(0.0, 0.0, r), 32, T)
        .unwrap()
        .point;
    // The boundary's highest points are the arcs' middles.
    let rim = r / lean.mul_add(lean, 1.0).sqrt();
    assert!(
        top.z > rim + 0.25,
        "the filling crowns above its boundary: {top:?}"
    );
    assert!(
        (top.z - r).abs() < 0.01,
        "and stands near the sphere: {top:?}"
    );

    let (shell, face) = sewn_with(&mut model, &strips, &face);
    let contacts = ogeom::fillet::analyse_blend(&model, &shell, &face, 41, T).unwrap();
    assert_eq!(
        contacts.len(),
        4,
        "the filling shares an edge with each strip"
    );
    for contact in &contacts {
        assert!(
            contact.tangency_error <= tolerance,
            "{} degrees from tangent",
            contact.tangency_error.to_degrees()
        );
        assert!(contact.gap <= tolerance, "{contact:?}");
    }
}
