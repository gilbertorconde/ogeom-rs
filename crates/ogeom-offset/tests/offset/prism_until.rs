//! A prism that stops on faces, against closed forms: a 10 x 4 rectangle
//! extruded up to a plane, up to a leaning plane, up to a cylinder lying
//! across it, up to the next face of a part with a step, and up to two
//! faces over different parts of it; and the lines it refuses.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_algo::sweep::roles::{SWEEP_BOTTOM, SWEEP_TOP};
use ogeom_algo::{Built, make_box, make_cylinder, make_face, make_polygon};
use ogeom_core::{Provenance, Role, Tolerances};
use ogeom_geom::PlaneSurface;
use ogeom_math::{Direction, Frame, Plane, Point, Vector};
use ogeom_mesh::Deflection;
use ogeom_offset::make_prism_until;
use ogeom_topo::{Filter, Model, NodeData, Shape, ShapeType, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A planar face through `corners`, its plane read from the first three.
fn planar_face(model: &mut Model, corners: &[Point]) -> Shape {
    let normal = (corners[1] - corners[0]).cross(corners[2] - corners[0]);
    let z = Direction::new(normal, T).unwrap();
    let x = Direction::new(corners[1] - corners[0], T).unwrap();
    let plane = Plane::new(Frame::new(corners[0], z, x, T).unwrap());
    let wire = make_polygon(model, corners, true, T).unwrap().shape;
    make_face(model, PlaneSurface::new(plane).into(), &[wire], T)
        .unwrap()
        .shape
}

/// The rectangle [x0, x1] x [0, 4] on z = 0.
fn rectangle(model: &mut Model, x0: f64, x1: f64) -> Shape {
    let corners = [(x0, 0.0), (x1, 0.0), (x1, 4.0), (x0, 4.0)].map(|(x, y)| Point::new(x, y, 0.0));
    planar_face(model, &corners)
}

/// A horizontal face over [x0, x1] x [-5, 10] at height `z`.
fn ceiling(model: &mut Model, x0: f64, x1: f64, z: f64) -> Shape {
    let corners =
        [(x0, -5.0), (x1, -5.0), (x1, 10.0), (x0, 10.0)].map(|(x, y)| Point::new(x, y, z));
    planar_face(model, &corners)
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    let fine = Deflection::with_chord(1e-3).unwrap();
    ogeom_algo::volume_properties(model, shape, fine, T)
        .unwrap()
        .mass
}

fn points_of(model: &Model, shape: &Shape) -> Vec<Point> {
    explore_unique(model, shape, ShapeType::Vertex)
        .unwrap()
        .iter()
        .map(|v| {
            let Some(NodeData::Vertex(data)) = model.node(v).map(|n| n.data()) else {
                unreachable!()
            };
            v.transform(model.datums()).unwrap().apply(data.point)
        })
        .collect()
}

/// The faces of `solid` carrying `role`, each with the shapes it is
/// derived from.
fn with_role(model: &Model, solid: &Shape, role: Role) -> Vec<(Shape, Vec<ogeom_core::EntityId>)> {
    explore_unique(model, solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter_map(|f| match model.provenance_of(&f) {
            Some(Provenance::Derived { role: r, from, .. }) if *r == role => {
                Some((f.clone(), from.to_vec()))
            }
            _ => None,
        })
        .collect()
}

/// Valid, of volume `want` within `bound`; prints the measure for the
/// record.
fn holds(model: &Model, built: &Built, want: f64, bound: f64, name: &str) {
    assert!(
        ogeom_algo::check(model, &built.shape, T)
            .unwrap()
            .is_valid(),
        "{name}: the prism is not valid"
    );
    let v = volume(model, &built.shape);
    eprintln!(
        "{name}: volume {v:.9} against {want:.9}, off {:.2e} (bound {bound:.1e})",
        (v - want).abs()
    );
    assert!((v - want).abs() <= bound, "{name}: {v} against {want}");
}

/// Up to a plane at z = 7 over the whole profile: 10 x 4 x 7 = 280. The
/// near end is derived from the profile, the far end from the limit, and
/// every profile edge generates a wall.
#[test]
fn a_prism_stops_on_a_plane_over_it() {
    let mut model = Model::new();
    let profile = rectangle(&mut model, 0.0, 10.0);
    let stop = ceiling(&mut model, -5.0, 20.0, 7.0);
    let built = make_prism_until(
        &mut model,
        &profile,
        Direction::Z,
        std::slice::from_ref(&stop),
        0.0,
        T,
    )
    .unwrap();
    holds(&model, &built, 280.0, 1e-9, "up to a plane");
    let near = with_role(&model, &built.shape, SWEEP_BOTTOM);
    let profile_id = model.identity_of(&profile).unwrap();
    assert_eq!(near.len(), 1);
    assert_eq!(near[0].1, vec![profile_id]);
    assert!(
        points_of(&model, &near[0].0)
            .iter()
            .all(|p| p.z.abs() < 1e-9)
    );
    let far = with_role(&model, &built.shape, SWEEP_TOP);
    let stop_id = model.identity_of(&stop).unwrap();
    assert_eq!(far.len(), 1);
    assert_eq!(far[0].1, vec![stop_id]);
    assert!(
        points_of(&model, &far[0].0)
            .iter()
            .all(|p| (p.z - 7.0).abs() < 1e-9)
    );
    assert!(
        built
            .history
            .generated(&stop)
            .iter()
            .any(|f| f.is_same(&far[0].0))
    );
    for edge in explore(&model, &profile, Filter::OfType(ShapeType::Edge)).unwrap() {
        assert!(
            !built.history.trace(&edge).is_empty(),
            "a profile edge generates no wall"
        );
    }
}

/// The offset moves the limit along the direction first: 2 short of the
/// plane at z = 7 is 10 x 4 x 5 = 200.
#[test]
fn an_offset_stops_short_of_the_limit() {
    let mut model = Model::new();
    let profile = rectangle(&mut model, 0.0, 10.0);
    let stop = ceiling(&mut model, -5.0, 20.0, 7.0);
    let built = make_prism_until(&mut model, &profile, Direction::Z, &[stop], -2.0, T).unwrap();
    holds(&model, &built, 200.0, 1e-9, "offset -2");
    for p in points_of(&model, &built.shape) {
        assert!(p.z.abs() < 1e-9 || (p.z - 5.0).abs() < 1e-9, "{p:?}");
    }
}

/// Up to the plane z = 6 + x / 2, leaning about 26.6 degrees: the volume is
/// the area times the height at the centroid, 40 x 8.5 = 340, and the far
/// end lies on the plane.
#[test]
fn a_prism_stops_on_a_leaning_plane() {
    let mut model = Model::new();
    let profile = rectangle(&mut model, 0.0, 10.0);
    let corners = [(-5.0, -5.0), (15.0, -5.0), (15.0, 10.0), (-5.0, 10.0)]
        .map(|(x, y)| Point::new(x, y, 6.0 + x / 2.0));
    let stop = planar_face(&mut model, &corners);
    let built = make_prism_until(&mut model, &profile, Direction::Z, &[stop], 0.0, T).unwrap();
    holds(&model, &built, 340.0, 1e-9, "up to a leaning plane");
    for (face, _) in with_role(&model, &built.shape, SWEEP_TOP) {
        for p in points_of(&model, &face) {
            assert!((p.z - (6.0 + p.x / 2.0)).abs() <= T.confusion(), "{p:?}");
        }
    }
}

/// Up to a cylinder of radius 5 lying along X at y = 2, z = 10, the whole
/// solid as the limit: over the profile (|y - 2| <= 2) a line at y stops at
/// 10 - sqrt(25 - (y - 2)^2), so the solid is
/// 10 (40 - 2 sqrt(21) - 25 asin(0.4)). Every far vertex lies on the
/// cylinder and nothing of the prism stands above its axis.
#[test]
fn a_prism_stops_under_a_cylinder_across_it() {
    let mut model = Model::new();
    let profile = rectangle(&mut model, 0.0, 10.0);
    let frame = Frame::new(Point::new(-5.0, 2.0, 10.0), Direction::X, Direction::Y, T).unwrap();
    let rod = make_cylinder(&mut model, frame, 5.0, 20.0, T)
        .unwrap()
        .shape;
    let built = make_prism_until(&mut model, &profile, Direction::Z, &[rod], 0.0, T).unwrap();
    let want = 10.0 * (40.0 - 2.0 * 21.0_f64.sqrt() - 25.0 * 0.4_f64.asin());
    holds(&model, &built, want, 1e-8, "up to a cylinder");
    let far = with_role(&model, &built.shape, SWEEP_TOP);
    assert!(!far.is_empty());
    for (face, _) in &far {
        for p in points_of(&model, face) {
            let r = Vector::new(0.0, p.y - 2.0, p.z - 10.0).magnitude();
            assert!((r - 5.0).abs() <= T.confusion(), "{p:?} at radius {r}");
        }
    }
    let top = ogeom_algo::shape_bounds(&model, &built.shape, T)
        .unwrap()
        .high()
        .unwrap();
    assert!(
        top.z < 10.0 - 4.0,
        "a piece above the cylinder survives: {top:?}"
    );
}

/// Up to the next face of a part whose underside steps from z = 5 (x < 5)
/// to z = 8 (x > 5), the whole part as the limit: the profile [1, 11] x
/// [0, 4] gives 4 (4 x 5 + 6 x 8) = 272, with a wall at the step.
#[test]
fn a_prism_stops_on_the_next_face_of_a_stepped_part() {
    let mut model = Model::new();
    let profile = rectangle(&mut model, 1.0, 11.0);
    let at = |x, y, z| Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap();
    let block = make_box(&mut model, at(-5.0, -5.0, 5.0), (25.0, 15.0, 15.0), T)
        .unwrap()
        .shape;
    let notch = make_box(&mut model, at(5.0, -10.0, 0.0), (20.0, 25.0, 8.0), T)
        .unwrap()
        .shape;
    let part = ogeom_bool::cut(&mut model, &block, &notch, T)
        .unwrap()
        .shape;
    let built = make_prism_until(&mut model, &profile, Direction::Z, &[part], 0.0, T).unwrap();
    holds(&model, &built, 272.0, 1e-9, "up to the next face");
    let heights: Vec<f64> = with_role(&model, &built.shape, SWEEP_TOP)
        .iter()
        .flat_map(|(f, _)| points_of(&model, f))
        .map(|p| p.z)
        .collect();
    assert!(heights.iter().any(|z| (z - 5.0).abs() < 1e-9));
    assert!(heights.iter().any(|z| (z - 8.0).abs() < 1e-9));
}

/// Two faces at different heights over parts of the profile, the higher
/// one reaching under the lower: z = 5 over x < 5 and z = 8 over x > 3.
/// Each part ends on the nearer face: 4 (4 x 5 + 6 x 8) = 272.
#[test]
fn each_part_of_the_profile_ends_on_the_nearer_of_two_faces() {
    let mut model = Model::new();
    let profile = rectangle(&mut model, 1.0, 11.0);
    let low = ceiling(&mut model, -5.0, 5.0, 5.0);
    let high = ceiling(&mut model, 3.0, 20.0, 8.0);
    let built = make_prism_until(
        &mut model,
        &profile,
        Direction::Z,
        &[high.clone(), low.clone()],
        0.0,
        T,
    )
    .unwrap();
    holds(&model, &built, 272.0, 1e-9, "two faces");
    let far = with_role(&model, &built.shape, SWEEP_TOP);
    let (low_id, high_id) = (
        model.identity_of(&low).unwrap(),
        model.identity_of(&high).unwrap(),
    );
    for (face, from) in &far {
        let z = points_of(&model, face)[0].z;
        let want = if (z - 5.0).abs() < 1e-9 {
            low_id
        } else {
            high_id
        };
        assert_eq!(from, &vec![want], "a far face at z = {z}");
    }
}

/// A limit over only part of the profile leaves lines that meet nothing:
/// refused, not run on.
#[test]
fn lines_that_meet_no_limit_are_refused() {
    let mut model = Model::new();
    let profile = rectangle(&mut model, 0.0, 10.0);
    let stop = ceiling(&mut model, -5.0, 6.0, 7.0);
    let refused = make_prism_until(&mut model, &profile, Direction::Z, &[stop], 0.0, T);
    assert!(
        matches!(refused, Err(ogeom_core::OgeomError::Construction(_))),
        "{refused:?}"
    );
}

/// A limit behind the profile stops nothing; the direction is never
/// turned round to reach it.
#[test]
fn a_limit_behind_the_profile_is_refused() {
    let mut model = Model::new();
    let profile = rectangle(&mut model, 0.0, 10.0);
    let below = ceiling(&mut model, -5.0, 20.0, -3.0);
    let refused = make_prism_until(&mut model, &profile, Direction::Z, &[below], 0.0, T);
    assert!(
        matches!(refused, Err(ogeom_core::OgeomError::Construction(_))),
        "{refused:?}"
    );
}

/// A cylinder of radius 1 across a profile 4 wide: its silhouette falls
/// over the profile, where lines only touch it and those beside it meet
/// nothing. Refused, not guessed.
#[test]
fn a_silhouette_over_the_profile_is_refused() {
    let mut model = Model::new();
    let profile = rectangle(&mut model, 0.0, 10.0);
    let frame = Frame::new(Point::new(-5.0, 2.0, 10.0), Direction::X, Direction::Y, T).unwrap();
    let rod = make_cylinder(&mut model, frame, 1.0, 20.0, T)
        .unwrap()
        .shape;
    let refused = make_prism_until(&mut model, &profile, Direction::Z, &[rod], 0.0, T);
    assert!(
        matches!(refused, Err(ogeom_core::OgeomError::Construction(_))),
        "{refused:?}"
    );
}
