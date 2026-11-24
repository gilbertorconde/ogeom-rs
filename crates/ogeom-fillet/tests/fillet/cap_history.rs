//! A concave strip extends an inherited cap without generating it again.
#![allow(clippy::unwrap_used, reason = "test constructions")]
use ogeom_core::Tolerances;
use ogeom_geom::SurfaceGeometry;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, ShapeType, explore_unique};
const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_ruled_round_keeps_one_inherited_cap_at_each_end() {
    let mut model = Model::new();
    let first = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 3.58, 14.0, T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(5.5, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let second = ogeom_algo::make_cylinder(&mut model, frame, 3.58, 14.0, T)
        .unwrap()
        .shape;
    let solid = ogeom_bool::fuse(&mut model, &first, &second, T)
        .unwrap()
        .shape;
    let y = -(3.58_f64.powi(2) - 2.75_f64.powi(2)).sqrt();
    let edge = explore_unique(&model, &solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|edge| {
            ogeom_algo::edge_vertices(&model, edge)
                .unwrap()
                .is_some_and(|(a, b)| {
                    [&a, &b].iter().all(|v| {
                        let p = v
                            .transform(model.datums())
                            .unwrap()
                            .apply(model.node(v).unwrap().data().as_vertex().unwrap().point);
                        (p.x - 2.75).abs() < T.confusion() && (p.y - y).abs() < T.confusion()
                    })
                })
        })
        .unwrap();
    let rounded =
        ogeom_fillet::fillet_edges(&mut model, &solid, std::slice::from_ref(&edge), 1.0, T)
            .unwrap();
    assert!(
        ogeom_algo::check(&model, &rounded.shape, T)
            .unwrap()
            .is_valid()
    );
    let mut caps = 0;
    for face in explore_unique(&model, &rounded.shape, ShapeType::Face).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        if matches!(
            model.geometry().surface(data.surface),
            Some(SurfaceGeometry::Plane(_))
        ) {
            caps += 1;
            assert!(
                !rounded
                    .history
                    .generated(&edge)
                    .iter()
                    .any(|f| f.is_same(&face)),
                "an inherited cap is not a generated blend face"
            );
            // An actual original cap must trace here. Newness is decided by
            // history, not by the cap's surface kind or geometric closeness.
            assert!(
                explore_unique(&model, &solid, ShapeType::Face)
                    .unwrap()
                    .iter()
                    .any(|given| rounded
                        .history
                        .trace(given)
                        .iter()
                        .any(|f| f.is_same(&face)))
            );
        }
    }
    assert_eq!(caps, 2, "one connected planar cap at each end");
}
