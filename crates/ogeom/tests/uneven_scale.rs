//! Spheres and revolved solids scaled unevenly: ellipsoids and parts of
//! them, measured and cut as what they are.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    general_transformed_shape, make_box, make_edge_between, make_face, make_revolution,
    make_sphere, make_vertex, make_wire, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, LineCurve, PlaneSurface};
use ogeom::math::{Axis, Circle, Direction, Frame, GeneralTransform, Plane, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::default(), T)
        .unwrap()
        .mass
}

/// A quarter disc of radius 1 in the XZ plane: out along X, round the arc
/// to the pole on Z, back down the axis.
fn quarter_disc(model: &mut Model) -> Shape {
    let [o, x, z] = [
        Point::ORIGIN,
        Point::new(1.0, 0.0, 0.0),
        Point::new(0.0, 0.0, 1.0),
    ]
    .map(|p| make_vertex(model, p).shape);
    let out = LineCurve::segment(Point::ORIGIN, Point::new(1.0, 0.0, 0.0), T).unwrap();
    let out = make_edge_between(model, out.into(), (0.0, 1.0), &o, &x, T)
        .unwrap()
        .shape;
    // The arc in the XZ plane about the origin, from X up to Z.
    let frame = Frame::new(Point::ORIGIN, -Direction::Y, Direction::X, T).unwrap();
    let arc = CircleCurve::new(Circle::new(frame, 1.0, T).unwrap());
    let arc = make_edge_between(model, arc.into(), (0.0, PI / 2.0), &x, &z, T)
        .unwrap()
        .shape;
    let down = LineCurve::segment(Point::new(0.0, 0.0, 1.0), Point::ORIGIN, T).unwrap();
    let down = make_edge_between(model, down.into(), (0.0, 1.0), &z, &o, T)
        .unwrap()
        .shape;
    let wire = make_wire(model, &[out, arc, down], T).unwrap().shape;
    let plane = Plane::new(Frame::new(Point::ORIGIN, -Direction::Y, Direction::X, T).unwrap());
    make_face(model, PlaneSurface::new(plane).into(), &[wire], T)
        .unwrap()
        .shape
}

/// A revolved hemisphere and quarter sweep scaled twice as long in x: half
/// and an eighth of the ellipsoid, measured as such.
#[test]
fn a_revolved_solid_scales_to_the_ellipsoid_s_volume() {
    let ellipsoid = 4.0 / 3.0 * PI * 2.0;
    for (angle, share) in [(2.0 * PI, 0.5), (PI / 2.0, 0.125)] {
        let mut model = Model::new();
        let quarter = quarter_disc(&mut model);
        let axis = Axis::new(Point::ORIGIN, Direction::Z);
        let solid = make_revolution(&mut model, &quarter, axis, angle, T)
            .unwrap()
            .shape;
        let stretch = GeneralTransform::scaling_xyz(2.0, 1.0, 1.0);
        let scaled = general_transformed_shape(&mut model, &solid, &stretch, T)
            .unwrap()
            .shape;
        let diagnosis = ogeom::algo::check(&model, &scaled, T).unwrap();
        assert!(diagnosis.is_valid(), "{angle}: {diagnosis}");
        let want = ellipsoid * share;
        let v = volume(&model, &scaled);
        assert!(
            (v - want).abs() < want * 1e-3,
            "{v} against {want} at {angle}"
        );
        // Every point of a fine mesh on the ellipsoid, its flat faces aside.
        let mesh =
            ogeom::mesh::triangulate(&model, &scaled, Deflection::with_chord(1e-3).unwrap(), T)
                .unwrap();
        for p in &mesh.positions {
            let on = (p.x / 2.0).powi(2) + p.y * p.y + p.z * p.z;
            let flat = p.z.abs() < 1e-9 || (angle < PI && (p.x.abs() < 1e-9 || p.y.abs() < 1e-9));
            assert!(flat || (on - 1.0).abs() < 1e-3, "{p:?} off the ellipsoid");
        }
    }
}

/// An ellipsoid from a scaled sphere, cut by a slab over its upper half:
/// the upper half, and nothing below the slab.
#[test]
fn a_scaled_sphere_cuts_as_the_ellipsoid_it_is() {
    let mut model = Model::new();
    let unit = make_sphere(&mut model, Frame::WORLD, 1.0, T).unwrap().shape;
    let stretch = GeneralTransform::scaling_xyz(8.0, 5.0, 3.0);
    let whole = general_transformed_shape(&mut model, &unit, &stretch, T)
        .unwrap()
        .shape;
    let slab_frame =
        Frame::new(Point::new(-17.0, -17.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let slab = make_box(&mut model, slab_frame, (34.0, 34.0, 3.0), T)
        .unwrap()
        .shape;
    let half = ogeom::boolean::common(&mut model, &whole, &slab, T)
        .unwrap()
        .shape;
    assert!(ogeom::algo::check(&model, &half, T).unwrap().is_valid());
    let want = 2.0 / 3.0 * PI * 120.0;
    let v = volume(&model, &half);
    assert!((v - want).abs() < want * 1e-3, "{v} against {want}");
    let mesh = ogeom::mesh::triangulate(&model, &half, Deflection::default(), T).unwrap();
    let (lo, hi) = mesh
        .positions
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
            (lo.min(p.z), hi.max(p.z))
        });
    assert!(lo > -1e-6 && hi < 3.0 + 1e-6, "z runs {lo} .. {hi}");
}
