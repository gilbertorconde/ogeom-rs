//! §H of `docs/PLAN.md`: removing a set of faces, the wound closed from the
//! neighbours' own geometry — measured against the solids the features were
//! cut from, not against plausibility.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::SurfaceGeometry;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// Faces of `solid` whose surface satisfies the predicate.
fn faces_where(
    model: &Model,
    solid: &Shape,
    pred: impl Fn(&SurfaceGeometry) -> bool,
) -> Vec<Shape> {
    explore_unique(model, solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            model
                .node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(&pred)
        })
        .collect()
}

/// A through bore's wall removed: the rims are inner loops of the lid and
/// base, and the block comes back to the last bit of its exact volume —
/// no overshoot, because nothing is filled; the boundary is resewn.
#[test]
fn removing_a_bores_wall_makes_the_block_whole() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (30.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(15.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, frame, 4.0, 12.0, T)
        .unwrap()
        .shape;
    let drilled = ogeom::boolean::cut(&mut model, &block, &drill, T)
        .unwrap()
        .shape;

    let walls = faces_where(&model, &drilled, |s| {
        matches!(s, SurfaceGeometry::Cylinder(_))
    });
    assert_eq!(walls.len(), 1, "one bore wall");
    let built = ogeom::boolean::remove_faces(&mut model, &drilled, &walls, T).unwrap();

    let healed = volume(&model, &built.shape);
    // Within float accumulation, not within a modelling tolerance: the lid
    // and base come from the boolean's rebuild and carry its arithmetic. A
    // wrong closure would miss by the bore's volume, eight hundred times the
    // bound.
    assert!((healed - 6000.0).abs() < 1e-6, "the block, whole: {healed}");
    assert!(built.history.is_deleted(&walls[0]));
}

/// A chamfer band removed: the band interrupts the top and side faces, whose
/// planes re-intersect in the edge the chamfer replaced, and the end caps'
/// edges extend to the recovered corners.
#[test]
fn removing_a_chamfer_restores_the_sharp_box() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 15.0, 10.0), T)
        .unwrap()
        .shape;
    let edge = edge_near(&model, &block, Point::new(10.0, 0.0, 10.0));
    let chamfered = ogeom::fillet::chamfer_edge(&mut model, &block, &edge, 3.0, T)
        .unwrap()
        .shape;
    let before = volume(&model, &chamfered);
    assert!(before < 3000.0, "the chamfer took material: {before}");

    // The chamfer band is the one plane that is neither axis-aligned wall
    // nor lid: its normal has both y and z.
    let bands = faces_where(&model, &chamfered, |s| {
        if let SurfaceGeometry::Plane(p) = s {
            let n = p.plane().frame().z().vector();
            n.y.abs() > 0.1 && n.z.abs() > 0.1
        } else {
            false
        }
    });
    assert_eq!(bands.len(), 1, "one chamfer band");
    let built = ogeom::boolean::remove_faces(&mut model, &chamfered, &bands, T).unwrap();

    let healed = volume(&model, &built.shape);
    assert!(
        (healed - 3000.0).abs() < 1e-9,
        "the sharp box, exactly: {healed}"
    );
    // And it is a box again: six faces, twelve edges, eight vertices.
    let counts = |kind| explore_unique(&model, &built.shape, kind).unwrap().len();
    assert_eq!(counts(ShapeType::Face), 6);
    assert_eq!(counts(ShapeType::Edge), 12);
    assert_eq!(counts(ShapeType::Vertex), 8);
}

/// A fillet band removed: the same wound with a cylindrical band, and the
/// caps' arcs replaced by their own straight edges extended to the corner.
#[test]
fn removing_a_fillet_restores_the_sharp_box() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 15.0, 10.0), T)
        .unwrap()
        .shape;
    let edge = edge_near(&model, &block, Point::new(10.0, 0.0, 10.0));
    let filleted = ogeom::fillet::fillet_edge(&mut model, &block, &edge, 3.0, T)
        .unwrap()
        .shape;

    let bands = faces_where(&model, &filleted, |s| {
        matches!(s, SurfaceGeometry::Cylinder(_))
    });
    assert_eq!(bands.len(), 1, "one fillet band");
    let built = ogeom::boolean::remove_faces(&mut model, &filleted, &bands, T).unwrap();

    let healed = volume(&model, &built.shape);
    assert!(
        (healed - 3000.0).abs() < 1e-9,
        "the sharp box, exactly: {healed}"
    );
}

/// The refusals, by name.
#[test]
fn impossible_removals_are_refused_by_name() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let faces = explore_unique(&model, &block, ShapeType::Face).unwrap();

    let err = ogeom::boolean::remove_faces(&mut model, &block, &[], T).unwrap_err();
    assert!(err.to_string().contains("nothing to remove"), "{err}");

    let err = ogeom::boolean::remove_faces(&mut model, &block, &faces, T).unwrap_err();
    assert!(err.to_string().contains("nothing remains"), "{err}");

    // Removing one wall of a plain box: the wound has no feature geometry to
    // close it — the neighbours meet at right angles already, and nothing
    // stands in for the missing face.
    let err = ogeom::boolean::remove_faces(&mut model, &block, &faces[..1], T).unwrap_err();
    assert!(
        !err.to_string().is_empty(),
        "a wall removal fails with a reason, not silently"
    );
}

/// The box edge whose midpoint is nearest `at`.
fn edge_near(model: &Model, solid: &Shape, at: Point) -> Shape {
    use ogeom::geom::Curve3d as _;
    let mut best: Option<(f64, Shape)> = None;
    for edge in explore_unique(model, solid, ShapeType::Edge).unwrap() {
        let Some(data) = model.node(&edge).and_then(|n| n.data().as_edge()) else {
            continue;
        };
        let Some(ogeom::topo::EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
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

/// Two separate features named in one call remove as two wounds.
///
/// Two fillets on opposite edges share nothing; classifying their ring
/// edges together used to elect "the sides" across both features and die
/// recovering a nonsense edge. Grouped by shared edges, each feature runs
/// the whole machinery on the previous result, and the box comes back
/// sharp — exactly, both wounds.
#[test]
fn two_separate_fillets_remove_in_one_call() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let mut solid = block;
    for target in [Point::new(10.0, 5.0, 10.0), Point::new(0.0, 5.0, 0.0)] {
        let edge = edge_near(&model, &solid, target);
        solid = ogeom::fillet::fillet_edge(&mut model, &solid, &edge, 2.0, T)
            .unwrap()
            .shape;
    }
    // The two blend faces: the cylindrical ones.
    let mut blends = Vec::new();
    for face in explore_unique(&model, &solid, ShapeType::Face).unwrap() {
        let data = model.node(&face).unwrap().data().as_face().unwrap();
        let surface = model.geometry().surface(data.surface).unwrap();
        if matches!(surface, ogeom::geom::SurfaceGeometry::Cylinder(_)) {
            blends.push(face);
        }
    }
    assert_eq!(blends.len(), 2, "two fillets leave two cylindrical faces");
    let restored = ogeom::boolean::remove_faces(&mut model, &solid, &blends, T)
        .unwrap()
        .shape;
    assert!(
        ogeom::algo::check(&model, &restored, T).unwrap().is_valid(),
        "the restored box is a valid solid"
    );
    let volume =
        ogeom::algo::volume_properties(&model, &restored, ogeom::mesh::Deflection::default(), T)
            .unwrap()
            .mass;
    assert!(
        (volume - 1000.0).abs() < 1e-6,
        "both wounds close exactly: {volume} against 1000"
    );
}
