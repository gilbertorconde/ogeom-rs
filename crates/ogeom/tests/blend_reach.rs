//! A chamfer or fillet is refused where its setback runs past a face: 12 mm
//! back on a 10 mm face would cut through it into the solid beyond.
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

/// A fillet likewise: on a square edge the ball touches each face its
/// radius back, and 12 mm of it will not sit on a 10 mm face.
#[test]
fn a_fillet_wider_than_its_face_is_refused() {
    for (r, fits) in [(9.0, true), (10.0, true), (12.0, false)] {
        let mut model = Model::new();
        let (block, edge) = block_and_edge(&mut model);
        let result = ogeom::fillet::fillet_edges(&mut model, &block, &[edge], r, T);
        if fits {
            let shape = result.unwrap_or_else(|e| panic!("{r}: {e}")).shape;
            let got = volume_properties(&model, &shape, Deflection::default(), T)
                .unwrap()
                .mass;
            let want = 4000.0 - (1.0 - core::f64::consts::FRAC_PI_4) * r * r * 20.0;
            assert!(
                (got - want).abs() < want * 1e-4,
                "{r}: {got} against {want}"
            );
        } else {
            assert!(
                result.is_err(),
                "a {r} mm fillet cuts through its 10 mm face"
            );
        }
    }
}

/// On a drum's top rim the ball sets back its radius across the top disc:
/// radius 4 fits a disc of radius 5 and radius 6 does not.
#[test]
fn a_rim_fillet_wider_than_its_cap_is_refused() {
    for (r, fits) in [(4.0, true), (6.0, false)] {
        let mut model = Model::new();
        let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T)
            .unwrap()
            .shape;
        let rim = explore_unique(&model, &drum, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let b = ogeom::algo::shape_bounds(&model, e, T).unwrap();
                b.low().unwrap().z > 10.0 - 1e-3 && b.size().x > 9.0
            })
            .expect("the top rim");
        let result = ogeom::fillet::fillet_edges(&mut model, &drum, &[rim], r, T);
        assert_eq!(result.is_ok(), fits, "radius {r}: {:?}", result.err());
    }
}

/// The inside edge of an L, where the fillet adds material: 3 fits the
/// 4 mm of floor beside the wall and 5 does not.
#[test]
fn a_concave_fillet_wider_than_its_floor_is_refused() {
    for (r, fits) in [(3.0, true), (5.0, false)] {
        let mut model = Model::new();
        let pts = [
            (0.0, 0.0),
            (6.0, 0.0),
            (6.0, 2.0),
            (2.0, 2.0),
            (2.0, 8.0),
            (0.0, 8.0),
        ]
        .map(|(x, z)| Point::new(x, 0.0, z));
        let wire = make_polygon(&mut model, &pts, true, T).unwrap().shape;
        let plane = PlaneSurface::new(Plane::through(Point::ORIGIN, -ogeom::math::Direction::Y));
        let face = make_face(&mut model, plane.into(), &[wire], T)
            .unwrap()
            .shape;
        let l = make_prism(&mut model, &face, Vector::new(0.0, 10.0, 0.0), T)
            .unwrap()
            .shape;
        let inner = explore_unique(&model, &l, ShapeType::Edge)
            .unwrap()
            .into_iter()
            .find(|e| {
                let b = ogeom::algo::shape_bounds(&model, e, T).unwrap();
                let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
                (lo.x - 2.0).abs() < 1e-3
                    && (hi.x - 2.0).abs() < 1e-3
                    && (lo.z - 2.0).abs() < 1e-3
                    && (hi.z - 2.0).abs() < 1e-3
            })
            .expect("the inside edge");
        let result = ogeom::fillet::fillet_edges(&mut model, &l, &[inner], r, T);
        assert_eq!(result.is_ok(), fits, "radius {r}: {:?}", result.err());
    }
}
