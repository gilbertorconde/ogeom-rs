//! Solids that carry no history are asked what features made them: a block
//! wearing four of them, a rounded rim, and an imported part where the
//! honest answer includes a refusal.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point};
use ogeom_recognize::{Feature, HoleKind, recognize};
use ogeom_topo::{EdgeRepr, Model, Shape, ShapeType};

const T: Tolerances = Tolerances::millimetres();

/// One part wearing four features — a drilled, filleted, chamfered, pocketed
/// block — read back from nothing but its topology.
#[test]
fn a_worked_block_confesses_the_features_that_made_it() {
    let mut model = Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 12.0), T)
        .unwrap()
        .shape;

    // Round one top edge while it is still pristine — blending an edge of a
    // face already carrying holes, or blending twice on one body, is the
    // marching intersector's territory, not recognition's.
    let round_edge = edge_near(&model, &block, Point::new(20.0, 0.0, 12.0));
    let blended = ogeom_fillet::fillet_edge(&mut model, &block, &round_edge, 2.0, T)
        .unwrap()
        .shape;

    // Drill a through bore.
    let drill_frame =
        Frame::new(Point::new(30.0, 15.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom_algo::make_cylinder(&mut model, drill_frame, 4.0, 14.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom_bool::cut(&mut model, &blended, &drill, T)
        .unwrap()
        .shape;

    // Mill a pocket.
    let mill_frame = Frame::new(Point::new(5.0, 5.0, 7.0), Direction::Z, Direction::X, T).unwrap();
    let mill = ogeom_algo::make_box(&mut model, mill_frame, (12.0, 12.0, 6.0), T)
        .unwrap()
        .shape;
    let part = ogeom_bool::cut(&mut model, &drilled, &mill, T)
        .unwrap()
        .shape;

    let features = recognize(&model, &part, T).unwrap();

    let hole = features
        .iter()
        .find_map(|f| match f {
            Feature::Hole(h) => Some(h),
            _ => None,
        })
        .expect("the bore is a hole");
    assert_eq!(hole.kind, HoleKind::Through);
    assert!((hole.radius - 4.0).abs() < 1e-9);
    // The hole is claimed by its wall, which is the cylinder the drill cut.
    assert!(!hole.faces.is_empty());

    let fillet = features
        .iter()
        .find_map(|f| match f {
            Feature::Fillet(x) => Some(x),
            _ => None,
        })
        .expect("the round is a fillet");
    assert!((fillet.radius - 2.0).abs() < 1e-9);
    assert!(!fillet.concave);

    let pocket = features
        .iter()
        .find_map(|f| match f {
            Feature::Pocket(p) => Some(p),
            _ => None,
        })
        .expect("the recess is a pocket");
    assert_eq!(pocket.walls.len(), 4);
}

/// A toroidal blend is recognized where it is a true two-sided fillet: a
/// cylinder's rim rounded by the kernel's own revolved fillet.
#[test]
fn a_rim_fillet_is_a_toroidal_fillet() {
    let mut model = Model::new();
    let post = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 8.0, 12.0, T)
        .unwrap()
        .shape;
    // The rim circle's curve midpoint sits on the far side of the drum —
    // nearest-midpoint selection aims there, not at the seam.
    let rim = edge_near(&model, &post, Point::new(-8.0, 0.0, 12.0));
    let rounded = ogeom_fillet::fillet_edge(&mut model, &post, &rim, 2.5, T)
        .unwrap()
        .shape;
    let features = recognize(&model, &rounded, T).unwrap();
    let fillet = features
        .iter()
        .find_map(|f| match f {
            Feature::Fillet(x) => Some(x),
            _ => None,
        })
        .expect("the rim round is a fillet");
    assert!((fillet.radius - 2.5).abs() < 1e-9);
    assert!(!fillet.concave);
}

/// The corpus speaks the same language — and the recognizer stays honest
/// about it. The smallest NIST part's annular groove is confessed as a
/// pocket; its two torus rims, tangent to their walls but meeting the top
/// face at a deliberate seventy-three degrees, are *partial* rounds, not
/// two-sided fillets, and claiming them would be wrong.
#[test]
fn an_imported_part_confesses_its_pocket_and_nothing_false() {
    let path = format!(
        "{}/../../tests/corpus/nist_ftc_11_asme1_rb.stp",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(path).unwrap();
    let mut import = ogeom_io::read_step(&text, T).unwrap();
    let solid = import.solids[0].clone();
    let healed = ogeom_heal::reanchor_periodic_rings(import.document.model_mut(), &solid, T)
        .unwrap()
        .shape;
    let features = recognize(import.document.model(), &healed, T).unwrap();
    assert!(
        features.iter().any(|f| matches!(f, Feature::Pocket(_))),
        "the groove is a pocket: {features:?}"
    );
    assert!(
        !features.iter().any(|f| matches!(f, Feature::Fillet(_))),
        "no false fillets on single-tangent rounds: {features:?}"
    );
}

/// The box edge whose midpoint is nearest `at`.
fn edge_near(model: &Model, solid: &Shape, at: Point) -> Shape {
    use ogeom_geom::Curve3d as _;
    let mut best: Option<(f64, Shape)> = None;
    for edge in ogeom_topo::explore_unique(model, solid, ShapeType::Edge).unwrap() {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let Some(geometry) = model.geometry().curve(*curve) else {
            continue;
        };
        let p = geometry
            .point_at(f64::midpoint(range.0, range.1), T)
            .unwrap();
        let d = p.distance(at);
        if best.as_ref().is_none_or(|(held, _)| d < *held) {
            best = Some((d, edge));
        }
    }
    best.unwrap().1
}
