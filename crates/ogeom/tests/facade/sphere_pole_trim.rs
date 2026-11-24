//! A sphere face whose trim runs a hair beside a pole of its chart, built
//! by hand or read from a file that keeps the sphere's own axis: re-charted
//! about an axis whose poles stand clear of the trim, it meshes whole and
//! classifies right. A face whose trim runs through a pole, or keeps clear
//! of both, is left as it is.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{Containment, SolidBoundary, check, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve, PlaneSurface, SphereSurface, SurfaceGeometry};
use ogeom::math::{Circle, Direction, Frame, Plane, Point, Sphere, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{FaceData, Location, Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const R: f64 = 10.0;

/// The cap of the ball of radius `R` about the origin beyond the plane
/// `x = side * delta`, away from the origin (`side` 1 or -1), built by hand
/// on the sphere's own chart (axis `+z`): one rim circle, no seam or pole
/// edge, no pcurves yet. Its rim passes `delta` from both poles, which the
/// cap does not hold.
fn cap(model: &mut Model, delta: f64, side: f64) -> (Shape, Vec<Shape>) {
    let out = if side > 0.0 {
        Direction::X
    } else {
        -Direction::X
    };
    let at = side * delta;
    let rho = (R * R - delta * delta).sqrt();
    let rim_frame = Frame::new(Point::new(at, 0.0, 0.0), out, Direction::Z, T).unwrap();
    let rim = Curve::Circle(CircleCurve::new(Circle::new(rim_frame, rho, T).unwrap()));
    let start = ogeom::algo::make_vertex(model, Point::new(at, 0.0, rho)).shape;
    let edge = ogeom::algo::make_edge_between(
        model,
        rim,
        (0.0, core::f64::consts::TAU),
        &start,
        &start,
        T,
    )
    .unwrap()
    .shape;
    let ball: SurfaceGeometry = SphereSurface::new(Sphere::new(Frame::WORLD, R, T).unwrap()).into();
    let ball = model.geometry_mut().add_surface(ball);
    let flat: SurfaceGeometry = PlaneSurface::over(
        Plane::new(Frame::new(Point::new(at, 0.0, 0.0), -out, Direction::Z, T).unwrap()),
        (-2.0 * R, 2.0 * R),
        (-2.0 * R, 2.0 * R),
    )
    .unwrap()
    .into();
    let flat = model.geometry_mut().add_surface(flat);
    let mut faces = Vec::new();
    for (surface, use_) in [(ball, edge.clone()), (flat, edge.reversed())] {
        let wire = model.add_wire(&[use_]).unwrap();
        faces.push(
            model
                .add_face(FaceData::new(surface, Location::identity()), &[wire])
                .unwrap(),
        );
    }
    let shell = model.add_shell(&faces).unwrap();
    (model.add_solid(&[shell]).unwrap(), faces)
}

/// The cap's exact volume.
fn cap_volume(delta: f64) -> f64 {
    let h = R - delta;
    core::f64::consts::PI * h * h * (3.0 * R - h) / 3.0
}

/// The volume the solid's mesh encloses, at a chord of 1e-4.
fn mesh_volume(model: &Model, solid: &Shape) -> f64 {
    let mesh = ogeom::mesh::triangulate(
        model,
        solid,
        Deflection {
            chord: 1e-4,
            ..Deflection::default()
        },
        T,
    )
    .unwrap();
    mesh.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.positions[i as usize]);
            let (a, b, c) = (a.to_vector(), b.to_vector(), c.to_vector());
            a.dot(b.cross(c)) / 6.0
        })
        .sum()
}

/// The cap judged against its closed forms: valid, its volume exact, its
/// mesh enclosing it to the chord's own error, and every point of a grid
/// and of a cluster by each pole read on its side of the cap's boundary,
/// beyond a band of 1e-5 about it.
fn judged(name: &str, model: &Model, solid: &Shape, delta: f64, side: f64) {
    let diagnosis = check(model, solid, T).unwrap();
    assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
    let want = cap_volume(delta);
    let got = volume_properties(model, solid, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!(
        (got - want).abs() <= want * 1e-9,
        "{name}: volume {got} against {want}"
    );
    // A whole ball's mesh at this chord encloses 8.2e-6 less than it.
    let meshed = mesh_volume(model, solid);
    assert!(
        (meshed - want).abs() <= want * 2e-5,
        "{name}: the mesh encloses {meshed} against {want}"
    );
    let boundary = SolidBoundary::of(model, solid, 1e-3, T).unwrap();
    let mut points = Vec::new();
    for i in -4..=4 {
        for j in -4..=4 {
            for k in -4..=4 {
                points.push(Point::new(
                    2.3 * f64::from(i) + 0.1,
                    2.3 * f64::from(j) + 0.1,
                    2.3 * f64::from(k) + 0.1,
                ));
            }
        }
    }
    for pole in [-R, R] {
        for i in -3..=3 {
            for j in -3..=3 {
                points.push(Point::new(
                    0.03 * f64::from(i) + 0.005,
                    0.03 * f64::from(j),
                    pole * 0.999,
                ));
            }
        }
    }
    let mut wrong = Vec::new();
    for p in points {
        let to_ball = R - p.to_vector().magnitude();
        let to_plane = p.x * side - delta;
        if to_ball.abs() < 1e-5 || to_plane.abs() < 1e-5 {
            continue;
        }
        let inside = to_ball > 0.0 && to_plane > 0.0;
        let right = matches!(
            (boundary.holds(model, p, T), inside),
            (Ok(Containment::In), true) | (Ok(Containment::Out), false)
        );
        if !right {
            wrong.push(p);
        }
    }
    assert!(
        wrong.is_empty(),
        "{name}: {} points read wrong, e.g. {:?}",
        wrong.len(),
        wrong.first()
    );
}

/// Built by hand and given its trims by `fix_face_pcurves`, the sphere face
/// of either cap is re-charted clear of its trim, with the new chart's pole
/// inside it on the side the cap lies: `+x` for one, `-x` for the other.
/// On main the cap's mesh dropped 81 to 100 percent of it and the
/// classifier read points on the wrong side.
#[test]
fn a_cap_trimmed_beside_its_poles_meshes_whole_and_classifies_right() {
    for (delta, side) in [1e-6, 1e-5, 1e-4]
        .into_iter()
        .flat_map(|d| [(d, 1.0), (d, -1.0)])
    {
        let mut model = Model::new();
        let (solid, faces) = cap(&mut model, delta, side);
        let mut recharted = 0;
        for face in &faces {
            if ogeom::heal::fix_face_pcurves(&mut model, face, 1.0, T)
                .unwrap()
                .recharted
            {
                recharted += 1;
            }
        }
        assert_eq!(
            recharted, 1,
            "delta {delta}: the sphere face alone is re-charted"
        );
        judged(
            &format!("built, delta {delta}, side {side}"),
            &model,
            &solid,
            delta,
            side,
        );
    }
}

/// The cap written with no trims, as a file from a system that keeps the
/// sphere's own axis states it, and read back through STEP and IGES: the
/// reader re-charts the sphere face; on main the read cap's mesh dropped
/// up to all of it and the classifier read a third of the points wrong.
#[test]
fn a_cap_read_beside_its_poles_meshes_whole_and_classifies_right() {
    for (delta, side) in [1e-6, 1e-5, 1e-4]
        .into_iter()
        .flat_map(|d| [(d, 1.0), (d, -1.0)])
    {
        let mut model = Model::new();
        let (solid, _) = cap(&mut model, delta, side);
        let mut document = ogeom::doc::Document::over(model);
        document.add_part("cap", solid);
        let step = ogeom::io::read_step(&ogeom::io::write_step(&document, T).unwrap(), T).unwrap();
        let iges = ogeom::io::read_iges(&ogeom::io::write_iges(&document, T).unwrap(), T).unwrap();
        judged(
            &format!("STEP, delta {delta}, side {side}"),
            step.document.model(),
            &step.solids[0],
            delta,
            side,
        );
        judged(
            &format!("IGES, delta {delta}, side {side}"),
            iges.document.model(),
            &iges.solids[0],
            delta,
            side,
        );
    }
}

/// A face whose trim runs through both poles (a ball halved through its
/// axis), keeps clear of them (a cap well off the axis), or has none (a
/// whole ball) is left as it is.
#[test]
fn a_sphere_face_through_its_poles_or_clear_of_them_is_left_as_it_is() {
    let mut model = Model::new();
    let ball = ogeom::algo::make_sphere(&mut model, Frame::WORLD, R, T)
        .unwrap()
        .shape;
    let at = Frame::new(
        Point::new(-20.0, -20.0, -20.0),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    let block = ogeom::algo::make_box(&mut model, at, (20.0, 40.0, 40.0), T)
        .unwrap()
        .shape;
    let half = ogeom::boolean::cut(&mut model, &ball, &block, T)
        .unwrap()
        .shape;
    let (far, far_faces) = cap(&mut model, 5.0, 1.0);
    for face in &far_faces {
        ogeom::heal::fix_face_pcurves(&mut model, face, 1.0, T).unwrap();
    }
    for body in [&ball, &half, &far] {
        for face in explore_unique(&model, body, ShapeType::Face).unwrap() {
            let got = ogeom::boolean::rechart_clear_of_poles(&mut model, &face, T).unwrap();
            assert_eq!(got, ogeom::boolean::Recharted::Untouched);
        }
    }
    let _ = Vector::ZERO;
}
