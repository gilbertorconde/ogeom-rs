//! §11's feature vocabulary: draft about a neutral plane, and the form
//! features built as the compositions they are.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Plane, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// The planar face of `shape` whose plane passes through `on`.
fn planar_face_at(model: &Model, shape: &Shape, on: Point) -> Shape {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let ogeom::topo::NodeData::Face(data) = model.node(f).unwrap().data() else {
                return false;
            };
            let Some(ogeom::geom::SurfaceGeometry::Plane(p)) =
                model.geometry().surface(data.surface)
            else {
                return false;
            };
            p.plane().distance_to(on).abs() < 1e-9
        })
        .expect("a planar face there")
}

#[test]
fn four_walls_drafted_about_the_base_make_a_frustum() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let walls = [
        planar_face_at(&model, &block, Point::new(0.0, 10.0, 5.0)),
        planar_face_at(&model, &block, Point::new(20.0, 10.0, 5.0)),
        planar_face_at(&model, &block, Point::new(10.0, 0.0, 5.0)),
        planar_face_at(&model, &block, Point::new(10.0, 20.0, 5.0)),
    ];

    let angle: f64 = 5.0_f64.to_radians();
    let drafted = ogeom::offset::apply_draft(
        &mut model,
        &block,
        &walls,
        Plane::through(Point::ORIGIN, Direction::Z),
        Direction::Z,
        angle,
        T,
    )
    .unwrap()
    .shape;

    // The base keeps its 20 by 20 — that is what a neutral plane means —
    // and each wall leans in by tan(angle) per unit of height, so the
    // section at height z is a square of side 20 - 2 z tan(angle).
    let t = angle.tan();
    let exact = {
        let a = 20.0;
        // ∫₀¹⁰ (a - 2 t z)² dz
        let f = |z: f64| {
            let s = 2.0f64.mul_add(-(t * z), a);
            s * s
        };
        // Simpson over the interval: the integrand is a quadratic, so this
        // is exact.
        (10.0 / 6.0) * f(0.0).mul_add(1.0, 4.0f64.mul_add(f(5.0), f(10.0)))
    };
    let measured = volume(&model, &drafted);
    assert!(
        (measured - exact).abs() < exact * 1e-6,
        "the drafted block is the frustum its angle implies: {measured} \
         against {exact}"
    );

    // The base is untouched: its four corners are still the block's.
    let base = planar_face_at(&model, &drafted, Point::new(10.0, 10.0, 0.0));
    let mut bound = ogeom::math::Aabb::EMPTY;
    for v in explore_unique(&model, &base, ShapeType::Vertex).unwrap() {
        bound = bound.with_point(model.node(&v).unwrap().data().as_vertex().unwrap().point);
    }
    assert!(
        (bound.size().x - 20.0).abs() < 1e-9 && (bound.size().y - 20.0).abs() < 1e-9,
        "the neutral plane's section did not move: {:?}",
        bound.size()
    );
}

/// A profile face, in the plane `z = height`, over the given corners.
fn profile(model: &mut Model, corners: &[Point]) -> Shape {
    let wire = ogeom::algo::make_polygon(model, corners, true, T)
        .unwrap()
        .shape;
    let edges =
        ogeom::topo::explore(model, &wire, ogeom::topo::Filter::OfType(ShapeType::Edge)).unwrap();
    let plane = ogeom::geom::PlaneSurface::new(Plane::through(corners[0], Direction::Z));
    ogeom::algo::make_face_with_pcurves(model, plane.into(), &[edges], T)
        .unwrap()
        .shape
}

#[test]
fn the_form_features_are_the_compositions_they_say_they_are() {
    use ogeom::offset::Feature;
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (40.0, 30.0, 10.0), T)
        .unwrap()
        .shape;

    // A pocket: a square profile on the lid, sunk four deep.
    let lid = 10.0;
    let square = profile(
        &mut model,
        &[
            Point::new(5.0, 5.0, lid),
            Point::new(15.0, 5.0, lid),
            Point::new(15.0, 15.0, lid),
            Point::new(5.0, 15.0, lid),
        ],
    );
    let pocketed = ogeom::offset::feature_prism(
        &mut model,
        &block,
        &square,
        ogeom::math::Vector::new(0.0, 0.0, -4.0),
        Feature::Removed,
        T,
    )
    .unwrap();
    assert!(
        (volume(&model, &pocketed.shape) - (12000.0 - 400.0)).abs() < 1e-6,
        "a ten by ten pocket four deep: {}",
        volume(&model, &pocketed.shape)
    );
    // The profile is what a later edit will name, so it survives in the
    // history as having generated the result.
    assert!(
        !pocketed.history.generated(&square).is_empty(),
        "the feature remembers what it was made from"
    );

    // A pad: the same profile, standing three proud.
    let pad = profile(
        &mut model,
        &[
            Point::new(25.0, 5.0, lid),
            Point::new(35.0, 5.0, lid),
            Point::new(35.0, 15.0, lid),
            Point::new(25.0, 15.0, lid),
        ],
    );
    let padded = ogeom::offset::feature_prism(
        &mut model,
        &pocketed.shape,
        &pad,
        ogeom::math::Vector::new(0.0, 0.0, 3.0),
        Feature::Added,
        T,
    )
    .unwrap()
    .shape;
    assert!(
        (volume(&model, &padded) - (12000.0 - 400.0 + 300.0)).abs() < 1e-6,
        "and a ten by ten pad three high: {}",
        volume(&model, &padded)
    );

    // A slot: a profile that starts outside the block and is swept clean
    // through it, so the groove is open at both ends.
    let cutter = profile(
        &mut model,
        &[
            Point::new(0.0, 20.0, lid),
            Point::new(40.0, 20.0, lid),
            Point::new(40.0, 25.0, lid),
            Point::new(0.0, 25.0, lid),
        ],
    );
    let slotted = ogeom::offset::feature_slot(
        &mut model,
        &padded,
        &cutter,
        ogeom::math::Vector::new(0.0, 0.0, -1.0),
        2.0,
        T,
    )
    .unwrap()
    .shape;
    assert!(
        (volume(&model, &slotted) - (12000.0 - 400.0 + 300.0 - 40.0 * 5.0 * 2.0)).abs() < 1e-6,
        "a five wide groove two deep across the whole lid: {}",
        volume(&model, &slotted)
    );
}

#[test]
fn a_wire_projects_onto_the_faces_it_lands_on() {
    // A square drawn in the air above a cylinder's wall, dropped onto it:
    // each straight edge becomes the curve the wall makes of it, and the
    // pcurve rides with it so the face could be split along the result.
    let mut model = Model::new();
    let drum = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 10.0, 20.0, T)
        .unwrap()
        .shape;

    // A rectangle standing off the wall at radius 14, spanning a quarter of
    // the drum's circumference and four of its height.
    let corners: Vec<Point> = [(-0.4_f64, 6.0_f64), (0.4, 6.0), (0.4, 10.0), (-0.4, 10.0)]
        .iter()
        .map(|(a, z)| Point::new(14.0 * a.cos(), 14.0 * a.sin(), *z))
        .collect();
    let wire = ogeom::algo::make_polygon(&mut model, &corners, true, T)
        .unwrap()
        .shape;

    let (landed, built) =
        ogeom::offset::normal_projection(&mut model, &drum, &wire, 40, 1e-4, T).unwrap();
    assert!(
        !landed.is_empty(),
        "the wire hangs over the wall; something should land"
    );
    for stretch in &landed {
        // Everything projected sits on the wall, at its radius.
        let data = model.node(&stretch.edge).unwrap().data().as_edge().unwrap();
        let ogeom::topo::EdgeRepr::Curve3d { curve, range, .. } = data.curve3d().unwrap() else {
            unreachable!()
        };
        let geometry = model.geometry().curve(*curve).unwrap();
        for k in 0..=8 {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 8.0;
            let p = ogeom::geom::Curve3d::point_at(geometry, t, T).unwrap();
            assert!(
                (p.x.hypot(p.y) - 10.0).abs() < 1e-3,
                "a projected point is on the wall: {p:?}"
            );
        }
        assert!(
            stretch.tolerance < 1e-3,
            "the fit states what it cost: {}",
            stretch.tolerance
        );
    }
    assert!(
        !built.history.generated(&wire).is_empty() || !built.history.modified(&wire).is_empty(),
        "the projection says what became of the wire"
    );
}
