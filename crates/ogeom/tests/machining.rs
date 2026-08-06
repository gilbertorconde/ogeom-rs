//! §18's tail: the recognized tree read as machining operations, and the
//! features undone by the volumes they are.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::Model;

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_drilled_pocketed_block_reads_as_a_drill_and_a_mill_and_undoes_to_a_block() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 12.0), T)
        .unwrap()
        .shape;

    // A through bore of diameter eight, down the block.
    let drill_frame =
        Frame::new(Point::new(30.0, 15.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, drill_frame, 4.0, 14.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(&mut model, &block, &drill, T)
        .unwrap()
        .shape;

    // And a pocket milled into the lid.
    let mill_frame = Frame::new(Point::new(5.0, 5.0, 7.0), Direction::Z, Direction::X, T).unwrap();
    let mill = ogeom::algo::make_box(&mut model, mill_frame, (12.0, 10.0, 6.0), T)
        .unwrap()
        .shape;
    let part = ogeom::boolean::cut(&mut model, &drilled, &mill, T)
        .unwrap()
        .shape;

    let features = ogeom::recognize::recognize(&model, &part, T).unwrap();
    let tree = ogeom::recognize::feature_tree(&model, &part, T).unwrap();
    let plan = ogeom::recognize::manufacturing_plan(&model, &tree, T).unwrap();
    assert_eq!(plan.len(), features.len(), "one step per feature");

    // The bore is a drill of diameter eight, entering along the axis and
    // breaking through.
    let drill_step = plan
        .iter()
        .find_map(|step| match &step.operation {
            ogeom::recognize::Operation::Drill {
                diameter,
                through,
                depth,
                ..
            } => Some((*diameter, *through, *depth)),
            _ => None,
        })
        .expect("the bore reads as a drill");
    assert!(
        (drill_step.0 - 8.0).abs() < 1e-9,
        "diameter eight: {}",
        drill_step.0
    );
    assert!(drill_step.1, "and it breaks through");
    assert!(
        (drill_step.2 - 12.0).abs() < 1e-6,
        "as deep as the block: {}",
        drill_step.2
    );

    // The pocket is a mill, five deep, coming down against the lid.
    let mill_step = plan
        .iter()
        .find_map(|step| match &step.operation {
            ogeom::recognize::Operation::Mill {
                depth, approach, ..
            } => Some((*depth, *approach)),
            _ => None,
        })
        .expect("the pocket reads as a mill");
    assert!(
        mill_step.0.is_some_and(|d| (d - 5.0).abs() < 1e-6),
        "five deep: {:?}",
        mill_step.0
    );
    assert!(
        mill_step.1.vector().dot(ogeom::math::Vector::Z) < -0.999,
        "the cutter comes down: {:?}",
        mill_step.1
    );

    // And the features undo: filling the bore and the pocket gives the
    // block back, to the last cubic millimetre.
    let mut solid = part.clone();
    for feature in &features {
        if matches!(
            feature,
            ogeom::recognize::Feature::Hole(_) | ogeom::recognize::Feature::Pocket(_)
        ) {
            solid = ogeom::recognize::remove_feature(&mut model, &solid, feature, T)
                .unwrap()
                .shape;
        }
    }
    let volume = ogeom::algo::volume_properties(&model, &solid, Deflection::default(), T)
        .unwrap()
        .mass;
    // Whole again, up to the ten microns each filling tool stands past its
    // own opening — the overshoot `remove_feature` states and explains.
    // Over this bore and this pocket that is three and a half cubic
    // millimetres on fourteen thousand, and it is *above* the original
    // rather than below it, which is what filling means.
    assert!(
        volume > 14400.0 && volume - 14400.0 < 5.0,
        "the block, whole again: {volume}"
    );
}
