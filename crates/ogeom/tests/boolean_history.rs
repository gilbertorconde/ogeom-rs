//! Every face of a boolean's result is reachable through its history from
//! the input it came from, carried over unchanged or split, and the same
//! through the fillets and chamfers that compose booleans.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{Built, make_box};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn cube(model: &mut Model, at: (f64, f64, f64), size: (f64, f64, f64)) -> Shape {
    let frame = Frame::new(Point::new(at.0, at.1, at.2), Direction::Z, Direction::X, T).unwrap();
    make_box(model, frame, size, T).unwrap().shape
}

/// The result's faces no input reaches through `modified` or `generated`.
fn unreached(model: &Model, inputs: &[Shape], built: &Built) -> usize {
    let result = explore_unique(model, &built.shape, ShapeType::Face).unwrap();
    result
        .iter()
        .filter(|f| {
            !inputs.iter().any(|g| {
                built
                    .history
                    .modified(g)
                    .iter()
                    .chain(built.history.generated(g))
                    .any(|x| x.is_same(f))
            })
        })
        .count()
}

fn faces_of(model: &Model, shapes: &[&Shape]) -> Vec<Shape> {
    shapes
        .iter()
        .flat_map(|s| explore_unique(model, s, ShapeType::Face).unwrap())
        .collect()
}

/// A block sunk into a cube's top, and one standing on it, fused, cut and
/// taken in common: each result face, split or carried over whole, traces
/// back to an input face.
#[test]
fn every_boolean_face_traces_to_an_input_face() {
    type Op = fn(&mut Model, &Shape, &Shape, Tolerances) -> ogeom::core::OgeomResult<Built>;
    let ops: [(&str, Op); 3] = [
        ("fuse", ogeom::boolean::fuse),
        ("cut", ogeom::boolean::cut),
        ("common", ogeom::boolean::common),
    ];
    for (name, op) in ops {
        for (z, standing) in [(8.0, false), (10.0, true)] {
            // Standing on the top, the block shares no volume with the
            // cube: their common part is empty, and says so.
            if standing && name == "common" {
                continue;
            }
            let mut model = Model::new();
            let a = cube(&mut model, (0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
            let b = cube(&mut model, (5.0, 2.0, z), (3.0, 3.0, 4.0));
            let inputs = faces_of(&model, &[&a, &b]);
            let built = op(&mut model, &a, &b, T).unwrap();
            assert_eq!(unreached(&model, &inputs, &built), 0, "{name} at z {z}");
        }
    }
}

/// Rounded and bevelled edges: the faces the blend trims trace back to
/// the faces they were, and each blend to the edge it replaced.
#[test]
fn fillet_and_chamfer_faces_trace_to_their_inputs() {
    for bevel in [false, true] {
        let mut model = Model::new();
        let block = cube(&mut model, (0.0, 0.0, 0.0), (10.0, 10.0, 10.0));
        let edge = explore_unique(&model, &block, ShapeType::Edge).unwrap()[0].clone();
        let mut inputs = faces_of(&model, &[&block]);
        inputs.push(edge.clone());
        let built = if bevel {
            let spec = ogeom::fillet::Chamfer::Symmetric(1.0);
            ogeom::fillet::chamfer_edges_with(&mut model, &block, &[(edge, spec)], T).unwrap()
        } else {
            ogeom::fillet::fillet_edges(&mut model, &block, &[edge], 1.0, T).unwrap()
        };
        assert_eq!(
            unreached(&model, &inputs, &built),
            0,
            "{}",
            if bevel { "chamfer" } else { "fillet" }
        );
    }
}

/// A small pocket cut into one corner of a block with two bores: the faces
/// the pocket never reaches come through as exact copies, each the only
/// image of its input, on the same surface with the same area, and the
/// faces it cuts are modified, no copies.
#[test]
fn faces_a_boolean_never_touches_are_reported_as_exact_copies() {
    let mut model = Model::new();
    let block = cube(&mut model, (0.0, 0.0, 0.0), (20.0, 10.0, 5.0));
    let mut part = block;
    for x in [12.0, 16.0] {
        let at = Frame::new(Point::new(x, 5.0, -1.0), Direction::Z, Direction::X, T).unwrap();
        let bore = ogeom::algo::make_cylinder(&mut model, at, 1.0, 7.0, T)
            .unwrap()
            .shape;
        part = ogeom::boolean::cut(&mut model, &part, &bore, T)
            .unwrap()
            .shape;
    }
    let pocket = cube(&mut model, (-1.0, -1.0, 3.0), (4.0, 4.0, 5.0));
    let cut = ogeom::boolean::cut(&mut model, &part, &pocket, T).unwrap();
    let area = |model: &Model, f: &Shape| {
        ogeom::algo::surface_properties(model, f, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass
    };
    let (mut copied, mut changed) = (0, 0);
    for face in explore_unique(&model, &part, ShapeType::Face).unwrap() {
        let bound = ogeom::algo::shape_bounds(&model, &face, T).unwrap();
        let low = bound.low().unwrap();
        // The faces the pocket reaches: those within its corner.
        let reached = low.x < 3.0 && low.y < 3.0 && bound.high().unwrap().z > 3.0;
        match cut.history.copy_of(&face) {
            Some(copy) => {
                assert!(!reached, "a face the pocket cuts reported a copy");
                assert!(
                    explore_unique(&model, &cut.shape, ShapeType::Face)
                        .unwrap()
                        .iter()
                        .any(|f| f.is_same(copy))
                );
                assert!((area(&model, copy) - area(&model, &face)).abs() < 1e-9);
                copied += 1;
            }
            None => {
                assert!(reached, "an untouched face was not reported a copy");
                assert!(!cut.history.modified(&face).is_empty());
                changed += 1;
            }
        }
    }
    // The bores' walls, the far ends and the bottom come through; the top
    // and the two walls at the corner are cut.
    assert_eq!(changed, 3);
    assert!(copied >= 5, "{copied}");
}
