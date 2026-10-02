//! A cutting plane that holds a closed surface's seam: the section runs
//! exactly along the seam edge, which the surface's face keeps as its
//! boundary while the plane's face takes the section as its own. The two
//! must walk the same pieces of that one curve, or the shell does not
//! close. A plane across the seam leaves the piece it runs through one
//! face across it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, make_box, make_cylinder, make_sphere, make_torus, surface_properties, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::SurfaceGeometry;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

/// A box reaching far past the solid, from `corner` over `size`.
fn block(model: &mut Model, corner: (f64, f64, f64), size: (f64, f64, f64)) -> Shape {
    let at = Frame::new(
        Point::new(corner.0, corner.1, corner.2),
        Direction::Z,
        Direction::X,
        T,
    )
    .unwrap();
    make_box(model, at, size, T).unwrap().shape
}

/// Each solid's seam lies in the plane `y = 0` on the side `x > 0`. A half
/// space `y > 0` or `y < 0` cuts along it, and a quarter `x, y > 0` cuts
/// along it and across the solid at `x = 0`. Common and cut are valid and
/// share the volume as the tool does.
#[test]
fn a_cut_along_the_seam_closes_on_a_torus_a_cylinder_and_a_sphere() {
    let mut model = Model::new();
    let torus = make_torus(&mut model, Frame::WORLD, 10.0, 3.0, T)
        .unwrap()
        .shape;
    let below = Frame::new(Point::new(0.0, 0.0, -3.0), Direction::Z, Direction::X, T).unwrap();
    let cylinder = make_cylinder(&mut model, below, 5.0, 6.0, T).unwrap().shape;
    let ball = make_sphere(&mut model, Frame::WORLD, 6.0, T).unwrap().shape;
    let tools = [
        ((-20.0, 0.0, -20.0), (40.0, 20.0, 40.0), 0.5),
        ((-20.0, -20.0, -20.0), (40.0, 20.0, 40.0), 0.5),
        ((0.0, 0.0, -20.0), (20.0, 20.0, 40.0), 0.25),
    ];
    for (name, solid) in [("torus", &torus), ("cylinder", &cylinder), ("ball", &ball)] {
        let whole = volume(&model, solid);
        for (corner, size, share) in tools {
            let tool = block(&mut model, corner, size);
            let kept = ogeom::boolean::common(&mut model, solid, &tool, T)
                .unwrap_or_else(|e| panic!("{name} common with {corner:?}: {e}"))
                .shape;
            let left = ogeom::boolean::cut(&mut model, solid, &tool, T)
                .unwrap_or_else(|e| panic!("{name} cut with {corner:?}: {e}"))
                .shape;
            for (piece, want) in [(&kept, whole * share), (&left, whole * (1.0 - share))] {
                let diagnosis = check(&model, piece, T).unwrap();
                assert!(diagnosis.is_valid(), "{name} with {corner:?}: {diagnosis}");
                let got = volume(&model, piece);
                assert!(
                    (got - want).abs() < want * 1e-3,
                    "{name} with {corner:?}: {got} against {want}"
                );
            }
        }
    }
}

/// The faces of `piece` on a cylinder or a sphere, and their area
/// integrated on the exact surfaces.
fn round_faces(model: &Model, piece: &Shape) -> (usize, f64) {
    let mut count = 0;
    let mut area = 0.0;
    for face in explore_unique(model, piece, ShapeType::Face).unwrap() {
        if is_round(model, &face) {
            count += 1;
            area += surface_properties(model, &face, Deflection::with_chord(1e-3).unwrap(), T)
                .unwrap()
                .mass;
        }
    }
    (count, area)
}

fn is_round(model: &Model, face: &Shape) -> bool {
    let data = model.node(face).unwrap().data().as_face().unwrap();
    matches!(
        model.geometry().surface(data.surface).unwrap(),
        SurfaceGeometry::Cylinder(_) | SurfaceGeometry::Sphere(_)
    )
}

/// A drum and a ball halved by a plane square to their seams: each half
/// has one round face of half the round area, the half the seam runs
/// through included, and the solid's round face traces to it.
#[test]
fn a_solid_halved_square_to_its_seam_keeps_one_round_face_per_half() {
    let mut model = Model::new();
    let (radius, height) = (5.0, 6.0);
    let drum = make_cylinder(&mut model, Frame::WORLD, radius, height, T)
        .unwrap()
        .shape;
    let ball = make_sphere(&mut model, Frame::WORLD, radius, T)
        .unwrap()
        .shape;
    for (name, solid, round) in [
        ("drum", &drum, 2.0 * PI * radius * height),
        ("ball", &ball, 4.0 * PI * radius * radius),
    ] {
        let side = explore_unique(&model, solid, ShapeType::Face)
            .unwrap()
            .into_iter()
            .find(|f| is_round(&model, f))
            .unwrap();
        let whole = volume(&model, solid);
        // The seam stands at x = radius, y = 0, inside the box x >= 0.
        let tool = block(&mut model, (0.0, -20.0, -20.0), (20.0, 40.0, 40.0));
        let kept = ogeom::boolean::common(&mut model, solid, &tool, T).unwrap();
        let left = ogeom::boolean::cut(&mut model, solid, &tool, T).unwrap();
        for built in [&kept, &left] {
            let piece = &built.shape;
            let diagnosis = check(&model, piece, T).unwrap();
            assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
            let got = volume(&model, piece);
            assert!((got - whole / 2.0).abs() < whole * 1e-6, "{name}: {got}");
            let (count, area) = round_faces(&model, piece);
            assert_eq!(count, 1, "{name}");
            assert!((area - round / 2.0).abs() < round * 1e-9, "{name}: {area}");
            let traced = built.history.trace(&side);
            assert_eq!(traced.len(), 1, "{name}");
            assert!(
                explore_unique(&model, piece, ShapeType::Face)
                    .unwrap()
                    .iter()
                    .any(|f| f.is_same(&traced[0])),
                "{name}"
            );
        }
    }
}
