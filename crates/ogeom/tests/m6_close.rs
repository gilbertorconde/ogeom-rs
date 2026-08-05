//! M6's closing argument: the kernel is complete — sketches solve and
//! diagnose, the screen reaches into the model, and dumb solids confess
//! the features that made them.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Point2, Vector};
use ogeom::sketch::{Constraint, Sketch, SolveOptions};
use ogeom::topo::{EdgeRepr, Model, Shape, ShapeType};

const T: Tolerances = Tolerances::millimetres();

/// A dimensioned bracket profile solves to its dimensions, reports itself
/// exactly constrained — and when a contradictory dimension is added, the
/// solver refuses and *names* the two constraints that fight, which is the
/// M6 requirement stated in so many words.
#[test]
fn a_sketch_solves_and_an_overconstrained_one_names_the_fight() {
    let mut sketch = Sketch::new();
    let p0 = sketch.add_point(Point2::new(0.1, -0.3));
    let p1 = sketch.add_point(Point2::new(52.0, 2.0));
    let p2 = sketch.add_point(Point2::new(49.0, 21.0));
    let p3 = sketch.add_point(Point2::new(-2.0, 19.0));
    let bottom = sketch.add_line(p0, p1).unwrap();
    let right = sketch.add_line(p1, p2).unwrap();
    let top = sketch.add_line(p3, p2).unwrap();
    let left = sketch.add_line(p0, p3).unwrap();
    let centre = sketch.add_point(Point2::new(20.0, 8.0));
    let bore = sketch.add_circle(centre, 4.0).unwrap();

    sketch
        .constrain(Constraint::Fixed(p0, Point2::new(0.0, 0.0)))
        .unwrap();
    sketch.constrain(Constraint::Horizontal(bottom)).unwrap();
    sketch
        .constrain(Constraint::Distance(p0, p1, 50.0))
        .unwrap();
    sketch.constrain(Constraint::Vertical(right)).unwrap();
    sketch
        .constrain(Constraint::Distance(p1, p2, 20.0))
        .unwrap();
    sketch.constrain(Constraint::Horizontal(top)).unwrap();
    sketch.constrain(Constraint::Vertical(left)).unwrap();
    sketch.constrain(Constraint::Radius(bore, 5.0)).unwrap();
    sketch
        .constrain(Constraint::Fixed(centre, Point2::new(25.0, 10.0)))
        .unwrap();

    let solution = sketch.solve(SolveOptions::default()).unwrap();
    assert!(solution.converged, "residual {}", solution.residual);
    assert!(solution.diagnosis.is_well_constrained());
    assert_eq!(sketch.measure_distance(p0, p1).unwrap().round(), 50.0);
    assert_eq!(sketch.measure_radius(bore).unwrap().round(), 5.0);

    // Now the fight: a second, disagreeing width.
    let conflicting = sketch
        .constrain(Constraint::Distance(p0, p1, 55.0))
        .unwrap();
    let fought = sketch.solve(SolveOptions::default()).unwrap();
    assert!(!fought.converged, "a contradiction must not converge");
    assert_eq!(fought.diagnosis.conflicting.len(), 1);
    let group = &fought.diagnosis.conflicting[0];
    assert!(group.contains(&conflicting), "the new dimension is named");
    assert!(group.len() >= 2, "so is what it fights: {group:?}");
    for id in group {
        // Every named constraint prints as a sentence a person can act on.
        assert!(!sketch.describe(*id).unwrap().is_empty());
    }
}

/// One part wearing four features — a drilled, filleted, chamfered,
/// pocketed block — put through recognition and picking together: the ray
/// that strikes the bore lands on a face the hole feature claims.
#[test]
fn features_are_recognized_and_picking_reaches_into_them() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 12.0), T)
        .unwrap()
        .shape;

    // Round one top edge while it is still pristine — blending an edge of
    // a face already carrying holes, or blending twice on one body, is the
    // marching intersector's territory, not this milestone's.
    let round_edge = edge_near(&model, &block, Point::new(20.0, 0.0, 12.0));
    let blended = ogeom::fillet::fillet_edge(&mut model, &block, &round_edge, 2.0, T)
        .unwrap()
        .shape;

    // Drill a through bore.
    let drill_frame =
        Frame::new(Point::new(30.0, 15.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, drill_frame, 4.0, 14.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(&mut model, &blended, &drill, T)
        .unwrap()
        .shape;

    // Mill a pocket.
    let mill_frame = Frame::new(Point::new(5.0, 5.0, 7.0), Direction::Z, Direction::X, T).unwrap();
    let mill = ogeom::algo::make_box(&mut model, mill_frame, (12.0, 12.0, 6.0), T)
        .unwrap()
        .shape;
    let part = ogeom::boolean::cut(&mut model, &drilled, &mill, T)
        .unwrap()
        .shape;

    // --- recognition ------------------------------------------------------
    let features = ogeom::recognize::recognize(&model, &part, T).unwrap();
    let hole = features
        .iter()
        .find_map(|f| match f {
            ogeom::recognize::Feature::Hole(h) => Some(h),
            _ => None,
        })
        .expect("the bore is a hole");
    assert_eq!(hole.kind, ogeom::recognize::HoleKind::Through);
    assert!((hole.radius - 4.0).abs() < 1e-9);
    let fillet = features
        .iter()
        .find_map(|f| match f {
            ogeom::recognize::Feature::Fillet(x) => Some(x),
            _ => None,
        })
        .expect("the round is a fillet");
    assert!((fillet.radius - 2.0).abs() < 1e-9);
    assert!(!fillet.concave);
    let pocket = features
        .iter()
        .find_map(|f| match f {
            ogeom::recognize::Feature::Pocket(p) => Some(p),
            _ => None,
        })
        .expect("the recess is a pocket");
    assert_eq!(pocket.walls.len(), 4);

    // --- picking ----------------------------------------------------------
    let scene =
        ogeom::select::Pickable::build(&model, &part, ogeom::mesh::Deflection::default(), T)
            .unwrap();
    // Down the bore wall: the first strike is the bore itself, and the
    // stable triangle mapping names a face the hole feature claims.
    // Into the bore at a slant: straight down the void would graze along
    // the wall forever and strike nothing.
    let hit = scene
        .pick_first(
            ogeom::select::Ray {
                origin: Point::new(30.0, 15.0, 20.0),
                direction: Vector::new(0.0, 0.3, -1.0),
            },
            0.0,
        )
        .expect("the ray strikes the bore");
    let struck = scene.triangle_face(hit.triangle).unwrap();
    assert!(
        hole.faces.iter().any(|f| f.is_same(struck)),
        "the picked face belongs to the recognized hole"
    );

    // A marquee over the pocket corner selects its floor; the floor the
    // recognition named and the floor the marquee finds are the same face.
    let view = Frame::WORLD;
    let picked = scene.select_rectangle(
        &view,
        Point2::new(4.0, 4.0),
        Point2::new(18.5, 18.5),
        ogeom::select::Marquee::Crossing,
    );
    assert!(
        picked.iter().any(|s| s.is_same(&pocket.floor)),
        "the marquee reaches the pocket floor"
    );

    // Sub-shape granularity: an aperture near a pocket corner resolves to
    // a vertex.
    let corner = scene
        .pick_first(
            ogeom::select::Ray {
                origin: Point::new(5.01, 5.01, 20.0),
                direction: Vector::new(0.0, 0.0, -1.0),
            },
            0.2,
        )
        .expect("the corner is pickable");
    assert_eq!(corner.kind, ogeom::select::PickKind::Vertex);
}

/// A toroidal blend is recognized where it is a true two-sided fillet: a
/// cylinder's rim rounded by the kernel's own revolved fillet.
#[test]
fn a_rim_fillet_is_a_toroidal_fillet() {
    let mut model = Model::new();
    let post = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 8.0, 12.0, T)
        .unwrap()
        .shape;
    // The rim circle's curve midpoint sits on the far side of the drum —
    // nearest-midpoint selection aims there, not at the seam.
    let rim = edge_near(&model, &post, Point::new(-8.0, 0.0, 12.0));
    let rounded = ogeom::fillet::fillet_edge(&mut model, &post, &rim, 2.5, T)
        .unwrap()
        .shape;
    let features = ogeom::recognize::recognize(&model, &rounded, T).unwrap();
    let fillet = features
        .iter()
        .find_map(|f| match f {
            ogeom::recognize::Feature::Fillet(x) => Some(x),
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
    let mut import = ogeom::io::read_step(&text, T).unwrap();
    let solid = import.solids[0].clone();
    let healed = ogeom::heal::reanchor_periodic_rings(import.document.model_mut(), &solid, T)
        .unwrap()
        .shape;
    let features = ogeom::recognize::recognize(import.document.model(), &healed, T).unwrap();
    assert!(
        features
            .iter()
            .any(|f| matches!(f, ogeom::recognize::Feature::Pocket(_))),
        "the groove is a pocket: {features:?}"
    );
    assert!(
        !features
            .iter()
            .any(|f| matches!(f, ogeom::recognize::Feature::Fillet(_))),
        "no false fillets on single-tangent rounds: {features:?}"
    );
}

/// The box edge whose midpoint is nearest `at`.
fn edge_near(model: &Model, solid: &Shape, at: Point) -> Shape {
    use ogeom::geom::Curve3d as _;
    let mut best: Option<(f64, Shape)> = None;
    for edge in ogeom::topo::explore_unique(model, solid, ShapeType::Edge).unwrap() {
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
