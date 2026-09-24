//! A cutting plane that holds a closed surface's seam: the section runs
//! exactly along the seam edge, which the surface's face keeps as its
//! boundary while the plane's face takes the section as its own. The two
//! must walk the same pieces of that one curve, or the shell does not
//! close.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_box, make_cylinder, make_sphere, make_torus, volume_properties};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

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
