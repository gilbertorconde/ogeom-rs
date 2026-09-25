//! A chamfer is refused where its setback runs past a face: a 12 mm bevel
//! on a 10 mm face would cut through it into the solid beyond.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{make_face, make_polygon, make_prism, volume_properties};
use ogeom::core::Tolerances;
use ogeom::fillet::{Chamfer, chamfer_edges_with};
use ogeom::geom::PlaneSurface;
use ogeom::math::{Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn block_and_edge(model: &mut Model) -> (Shape, Shape) {
    let pts =
        [(0.0, 0.0), (20.0, 0.0), (20.0, 20.0), (0.0, 20.0)].map(|(x, y)| Point::new(x, y, 0.0));
    let wire = make_polygon(model, &pts, true, T).unwrap().shape;
    let face = make_face(
        model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let block = make_prism(model, &face, Vector::new(0.0, 0.0, 10.0), T)
        .unwrap()
        .shape;
    // The top front edge, (0, 0, 10) to (20, 0, 10).
    let edge = explore_unique(model, &block, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let b = ogeom::algo::shape_bounds(model, e, T).unwrap();
            let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
            hi.y < 1e-3 && lo.z > 10.0 - 1e-3 && hi.x - lo.x > 19.0
        })
        .expect("the top front edge");
    (block, edge)
}

#[test]
fn a_chamfer_wider_than_its_face_is_refused() {
    for (d, want) in [
        (9.0, Some(4000.0 - 0.5 * 81.0 * 20.0)),
        (10.0, Some(3000.0)),
        (12.0, None),
    ] {
        let mut model = Model::new();
        let (block, edge) = block_and_edge(&mut model);
        let result = chamfer_edges_with(&mut model, &block, &[(edge, Chamfer::Symmetric(d))], T);
        match want {
            Some(v) => {
                let shape = result.unwrap_or_else(|e| panic!("{d}: {e}")).shape;
                let got = volume_properties(&model, &shape, Deflection::default(), T)
                    .unwrap()
                    .mass;
                assert!((got - v).abs() < 1e-6, "{d}: {got}");
            }
            None => assert!(result.is_err(), "{d} mm cuts through its 10 mm face"),
        }
    }
}
