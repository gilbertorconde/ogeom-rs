//! Filling boundary edges with a fitted patch: the doubly ruled saddle is
//! the case with an exact answer, and the fit must land on it. The
//! N-sided filling also takes placed edges (a prism's far end), and its
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
