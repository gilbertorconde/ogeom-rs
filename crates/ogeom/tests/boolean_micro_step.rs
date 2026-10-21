//! A pocket a few hundred times its tolerance deep, cut through by a box:
//! the pocket's walls are far shorter than the floor the outline weld keeps
//! for a fitted section's hug, and a wall's outline must still read as the
//! wall, measured in space whatever its chart's stretch.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::{Tolerance, Tolerances};
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The pocket's depth: 50 nm, 250 times the edges' tolerance.
const DEPTH: f64 = 5e-5;

/// The tolerance every edge of the part states, a hair above the
/// confusion distance, as an exchange file states it.
const STATED: f64 = 2e-7;

fn block(model: &mut Model, corner: (f64, f64, f64), size: (f64, f64, f64)) -> Shape {
    let at = Frame::new(
        Point::new(corner.0, corner.1, corner.2),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    ogeom::algo::make_box(model, at, size, T).unwrap().shape
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// A 2 x 2 x 1 block with a pocket `DEPTH` deep in its top: a unit square,
/// or a disc of radius 0.5 whose wall is then a spline. Every edge states
/// `STATED`.
fn stepped(model: &mut Model, round: bool) -> Shape {
    let base = block(model, (0.0, 0.0, 0.0), (2.0, 2.0, 1.0));
    let pocket = if round {
        let at = Frame::new(
            Point::new(1.0, 1.0, 1.0 - DEPTH),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        ogeom::algo::make_cylinder(model, at, 0.5, 1.0, T)
            .unwrap()
            .shape
    } else {
        block(model, (0.5, 0.5, 1.0 - DEPTH), (1.0, 1.0, 1.0))
    };
    let part = ogeom::boolean::cut(model, &base, &pocket, T).unwrap().shape;
    let part = if round {
        ogeom::algo::to_nurbs(model, &part, T).unwrap().shape
    } else {
        part
    };
    for edge in explore_unique(model, &part, ShapeType::Edge).unwrap() {
        model.widen(&edge, Tolerance::new(STATED).unwrap()).unwrap();
    }
    part
}

/// A box over the pocket's wall at x = 1.5 (or where the disc's wall
/// crosses y = 0.8), its bottom flush with the top, flush with the floor,
/// through both, halfway down the wall, or spanning just the step: cut and
/// common are valid and together hold the part.
#[test]
fn a_box_through_a_pocket_wall_far_shorter_than_the_weld_cuts_valid() {
    for round in [false, true] {
        for (name, z0, z1) in [
            ("flush with the top", 0.0, 1.0),
            ("flush with the floor", 1.0 - DEPTH, 2.0),
            ("through", 0.5, 2.0),
            ("halfway down the wall", 1.0 - DEPTH / 2.0, 2.0),
            ("the step alone", 1.0 - DEPTH, 1.0),
        ] {
            let mut model = Model::with_tolerances(T);
            let part = stepped(&mut model, round);
            let whole = volume(&model, &part);
            let tool = block(&mut model, (1.2, 0.8, z0), (1.0, 0.4, z1 - z0));
            let mut total = 0.0;
            for op in ["cut", "common"] {
                let made = if op == "cut" {
                    ogeom::boolean::cut(&mut model, &part, &tool, T)
                } else {
                    ogeom::boolean::common(&mut model, &part, &tool, T)
                }
                .unwrap_or_else(|e| panic!("{op} {name} (round {round}): {e}"))
                .shape;
                let found = ogeom::algo::check(&model, &made, T).unwrap();
                assert!(found.is_valid(), "{op} {name} (round {round}): {found}");
                total += volume(&model, &made);
            }
            assert!(
                (total - whole).abs() < 1e-9 * whole,
                "{name} (round {round}): {total} against {whole}"
            );
        }
    }
}
