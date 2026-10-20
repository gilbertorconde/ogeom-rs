//! Solids fused with their own mirror image, through STEP: the image's
//! circles and ellipses stand in frames the mirror left-handed, and its
//! seams meet the original's at poles and along shared rims. Each reads
//! back as the solid it was, split or unified.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Transform, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn corner(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

/// The part of `solid` inside the box at `keep`, fused with its mirror
/// image in the plane through the origin square to `normal`.
fn fused_with_its_image(
    model: &mut Model,
    solid: &Shape,
    keep: (Frame, (f64, f64, f64)),
    normal: Direction,
) -> Shape {
    let keep = ogeom::algo::make_box(model, keep.0, keep.1, T)
        .unwrap()
        .shape;
    let half = ogeom::boolean::common(model, solid, &keep, T)
        .unwrap()
        .shape;
    let copy = ogeom::algo::copied(model, &half).unwrap().shape;
    let mirror = Transform::plane_mirror(Point::ORIGIN, normal);
    let image = ogeom::algo::transformed(model, &copy, mirror)
        .unwrap()
        .shape;
    ogeom::boolean::fuse(model, &half, &image, T).unwrap().shape
}

/// `solid` written to STEP and read back: valid, drawing closed, and
/// measuring `closed_form`. The volume is an exact integral over the faces
/// read back, so it meets the closed form to rounding; 1e-9 relative leaves
/// room for the sums.
fn round_trips(name: &str, model: &Model, solid: &Shape, closed_form: f64) {
    let mut document = ogeom::doc::Document::over(model.clone());
    document.add_part("part", solid.clone());
    let text = ogeom::io::write_step(&document, T).unwrap();
    let import = ogeom::io::read_step(&text, T).unwrap();
    assert_eq!(import.solids.len(), 1, "{name}");
    let (back, solid) = (import.document.model(), &import.solids[0]);
    let diagnosis = ogeom::algo::check(back, solid, T).unwrap();
    assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
    let mesh = ogeom::mesh::triangulate(back, solid, Deflection::default(), T).unwrap();
    assert!(mesh.is_closed(), "{name}: the solid read back draws closed");
    let volume = ogeom::algo::volume_properties(back, solid, Deflection::default(), T)
        .unwrap()
        .mass;
    assert!(
        (volume - closed_form).abs() < closed_form * 1e-9,
        "{name}: {volume} against {closed_form}"
    );
}

/// A ball about the x axis halved at its equator, a torus halved through
/// its middle plane, and a prism on an ellipse halved across its long
/// axis, each fused with its mirror image, then unified: every one, split
/// and unified, reads back from STEP as the solid it was.
#[test]
fn solids_fused_with_their_mirror_image_round_trip_through_step() {
    use core::f64::consts::PI;
    let (ball, (major, minor), (a, b, h)) = (5.0, (10.0, 3.0), (4.0, 2.0, 5.0));
    let cases: [(&str, f64); 3] = [
        ("ball", 4.0 / 3.0 * PI * ball * ball * ball),
        ("torus", 2.0 * PI * PI * major * minor * minor),
        ("elliptic prism", PI * a * b * h),
    ];
    for (name, closed_form) in cases {
        let mut model = Model::new();
        let fused = match name {
            "ball" => {
                let polar = Frame::new(Point::ORIGIN, Direction::X, Direction::Y, T).unwrap();
                let solid = ogeom::algo::make_sphere(&mut model, polar, ball, T)
                    .unwrap()
                    .shape;
                let keep = (corner(0.0, -6.0, -6.0), (6.0, 12.0, 12.0));
                fused_with_its_image(&mut model, &solid, keep, Direction::X)
            }
            "torus" => {
                let solid = ogeom::algo::make_torus(&mut model, Frame::WORLD, major, minor, T)
                    .unwrap()
                    .shape;
                let keep = (corner(-15.0, -15.0, 0.0), (30.0, 30.0, 5.0));
                fused_with_its_image(&mut model, &solid, keep, Direction::Z)
            }
            _ => {
                let frame = Frame::new(Point::ORIGIN, Direction::Z, Direction::Y, T).unwrap();
                let ellipse = ogeom::math::Ellipse::new(frame, a, b, T).unwrap();
                let rim = ogeom::algo::make_edge(
                    &mut model,
                    ogeom::geom::Curve::Ellipse(ogeom::geom::EllipseCurve::new(ellipse)),
                    (0.0, 2.0 * PI),
                    T,
                )
                .unwrap()
                .shape;
                let wire = ogeom::algo::make_wire(&mut model, &[rim], T).unwrap().shape;
                let plane = ogeom::geom::PlaneSurface::new(ogeom::math::Plane::new(Frame::WORLD));
                let base = ogeom::algo::make_face(&mut model, plane.into(), &[wire], T)
                    .unwrap()
                    .shape;
                let solid = ogeom::algo::make_prism(&mut model, &base, Vector::new(0.0, 0.0, h), T)
                    .unwrap()
                    .shape;
                let keep = (corner(0.0, -6.0, -1.0), (6.0, 12.0, h + 2.0));
                fused_with_its_image(&mut model, &solid, keep, Direction::X)
            }
        };
        round_trips(&format!("{name}, split"), &model, &fused, closed_form);
        let unified = ogeom::heal::unify_same_domain(&mut model, &fused, T)
            .unwrap()
            .0
            .shape;
        round_trips(&format!("{name}, unified"), &model, &unified, closed_form);
    }
}

/// A ball of radius 5 bounded by one face whose seam is two edges running
/// away from a vertex on the equator, one to each pole: the file drops the
/// pole edges, and each seam's second use leaves a pole. It reads back
/// whole.
#[test]
fn a_seam_of_two_edges_running_from_one_vertex_round_trips_through_step() {
    let path = format!(
        "{}/../../tests/corpus/sphere_seam_two_ways.ogeom",
        env!("CARGO_MANIFEST_DIR")
    );
    let (model, roots) = ogeom::io::native::read(&std::fs::read_to_string(path).unwrap()).unwrap();
    let ball = 4.0 / 3.0 * core::f64::consts::PI * 125.0;
    round_trips("ball", &model, &roots[0], ball);
}
