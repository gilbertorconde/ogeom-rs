//! A serialized shape lands in a live model and behaves as if built there:
//! booleans take it, identity answers about it, and the remap table keeps
//! references recorded against the source document alive.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::io::native::{WriteOptions, read_into, write};
use ogeom::math::{Direction, Frame, Point};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, ogeom::mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// Box A, serialized from its own model, and the text that carries it.
fn boxed_tool() -> (Model, Shape, String) {
    let mut model = Model::new();
    let a = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let text = write(&model, std::slice::from_ref(&a), WriteOptions::default()).unwrap();
    (model, a, text)
}

#[test]
fn a_serialized_box_cuts_and_fuses_with_a_live_one() {
    // The reason absorb exists: a tool body from persistence meeting a live
    // one in a boolean, which demands both operands in one model.
    let (_, _, text) = boxed_tool();

    let mut model = Model::new();
    let b = ogeom::algo::make_box(
        &mut model,
        Frame::new(Point::new(5.0, 5.0, 5.0), Direction::Z, Direction::X, T).unwrap(),
        (10.0, 10.0, 10.0),
        T,
    )
    .unwrap()
    .shape;

    let absorbed = read_into(&mut model, &text).unwrap();
    let a = &absorbed.shapes[0];

    // The overlap is the 5x5x5 corner cube: fuse keeps one copy of it, cut
    // removes it from the live box.
    let fused = ogeom::boolean::fuse(&mut model, a, &b, T).unwrap();
    assert!(
        ogeom::algo::check(&model, &fused.shape, T)
            .unwrap()
            .is_valid()
    );
    let measured = volume(&model, &fused.shape);
    assert!(
        (measured - 1875.0).abs() < 1e-6,
        "fused volume {measured} against 1875"
    );

    let cut = ogeom::boolean::cut(&mut model, &b, a, T).unwrap();
    assert!(
        ogeom::algo::check(&model, &cut.shape, T)
            .unwrap()
            .is_valid()
    );
    let measured = volume(&model, &cut.shape);
    assert!(
        (measured - 875.0).abs() < 1e-6,
        "cut volume {measured} against 875"
    );
}

#[test]
fn an_absorbed_solid_still_says_where_its_faces_came_from() {
    let (source, a, text) = boxed_tool();

    let mut model = Model::new();
    ogeom::algo::make_box(&mut model, Frame::WORLD, (3.0, 3.0, 3.0), T).unwrap();
    let absorbed = read_into(&mut model, &text).unwrap();
    let landed = &absorbed.shapes[0];

    let before = explore_unique(&source, &a, ShapeType::Face).unwrap();
    let after = explore_unique(&model, landed, ShapeType::Face).unwrap();
    assert_eq!(before.len(), after.len());
    for (theirs, ours) in before.iter().zip(&after) {
        // The identity a caller recorded against the source document resolves
        // through the remap table to the same entity here.
        let recorded = source.identity_of(theirs).expect("faces carry identities");
        let here = absorbed.entities[&recorded];
        assert_eq!(model.identity_of(ours), Some(here));
        let found = model.shape_of(here).expect("the entity is findable");
        assert!(found.is_partner(ours), "identity found the wrong face");

        // And the face still says what it is for: same role, same kind of
        // provenance as the source document recorded.
        assert_eq!(
            model
                .provenance_of(ours)
                .and_then(ogeom::core::Provenance::role),
            source
                .provenance_of(theirs)
                .and_then(ogeom::core::Provenance::role),
        );
        assert_eq!(
            model.roots_of(ours).len(),
            source.roots_of(theirs).len(),
            "the derivation walk changed shape"
        );
    }
}

#[test]
fn multiple_roots_read_into_one_model_all_resolve() {
    let mut source = Model::new();
    let block = ogeom::algo::make_box(&mut source, Frame::WORLD, (2.0, 3.0, 4.0), T)
        .unwrap()
        .shape;
    let drum = ogeom::algo::make_cylinder(
        &mut source,
        Frame::new(Point::new(20.0, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap(),
        2.0,
        5.0,
        T,
    )
    .unwrap()
    .shape;
    let text = write(&source, &[block, drum], WriteOptions::default()).unwrap();

    let mut model = Model::new();
    ogeom::algo::make_box(&mut model, Frame::WORLD, (1.0, 1.0, 1.0), T).unwrap();
    let absorbed = read_into(&mut model, &text).unwrap();
    assert_eq!(absorbed.shapes.len(), 2);

    let expected = [24.0, std::f64::consts::PI * 4.0 * 5.0];
    for (shape, expected) in absorbed.shapes.iter().zip(expected) {
        assert!(ogeom::algo::check(&model, shape, T).unwrap().is_valid());
        let measured = ogeom::algo::volume_properties(
            &model,
            shape,
            ogeom::mesh::Deflection {
                chord: 1e-4,
                ..ogeom::mesh::Deflection::default()
            },
            T,
        )
        .unwrap()
        .mass;
        assert!(
            (measured - expected).abs() / expected < 1e-3,
            "absorbed volume {measured} against {expected}"
        );
    }
}
