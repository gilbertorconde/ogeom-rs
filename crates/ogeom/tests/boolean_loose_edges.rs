//! A face with one loose edge welds what is near that edge, and a section
//! crossing another, tight, edge of the face a few hundredths from its
//! vertex keeps the stretch it cuts there, as the face across that edge
//! does.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_box, make_cylinder, volume_properties};
use ogeom::core::{Tolerance, Tolerances};
use ogeom::geom::Curve3d as _;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{EdgeRepr, Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The edge of `shape` whose middle `keep` accepts.
fn edge_where(model: &Model, shape: &Shape, keep: impl Fn(Point) -> bool) -> Shape {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let data = model.node(e).unwrap().data().as_edge().unwrap();
            let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
                return false;
            };
            let curve = model.geometry().curve(*curve).unwrap();
            keep(curve.point_at(f64::midpoint(range.0, range.1), T).unwrap())
        })
        .unwrap()
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

/// A block 20 by 20 by 10 whose front face's top edge states 0.017, as a
/// converted mesh's threaded edge does, drilled upright so the drill's
/// wall crosses the front face's bottom edge `d` from its corner. The cut
/// and the common are valid and add up to the block.
#[test]
fn a_drill_crossing_a_tight_edge_near_its_vertex_cuts_beside_a_loose_one() {
    for d in [0.02, 0.05] {
        let mut model = Model::new();
        let block = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let loose = edge_where(&model, &block, |p| p.y == 0.0 && p.z == 10.0);
        model.widen(&loose, Tolerance::new(0.017).unwrap()).unwrap();
        let at = Frame::new(
            Point::new(1.0 + d, 0.5, -5.0),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        let drill = make_cylinder(&mut model, at, 1.25_f64.sqrt(), 20.0, T)
            .unwrap()
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &block, &drill, T)
            .unwrap_or_else(|e| panic!("cut {d} from the corner: {e}"))
            .shape;
        let common = ogeom::boolean::common(&mut model, &block, &drill, T)
            .unwrap_or_else(|e| panic!("common {d} from the corner: {e}"))
            .shape;
        for (name, made) in [("cut", &cut), ("common", &common)] {
            let diagnosis = check(&model, made, T).unwrap();
            assert!(
                diagnosis.is_valid(),
                "{name} {d} from the corner: {diagnosis}"
            );
        }
        let total = volume(&model, &cut) + volume(&model, &common);
        assert!(
            (total - 4000.0).abs() < 4000.0 * 1e-6,
            "{d} from the corner: {total}"
        );
    }
}
