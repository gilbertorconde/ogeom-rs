//! Filling four boundary edges with a fitted patch: the doubly ruled
//! saddle is the case with an exact answer, and the fit must land on it.
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
