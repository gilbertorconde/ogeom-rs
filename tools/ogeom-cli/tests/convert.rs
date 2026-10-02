//! `ogeom-cli convert` reads a mesh and writes the solid it rebuilds.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A bored block meshed to STL, converted from the command line and written
/// as STEP, reads back one valid solid on six planes and a cylinder.
#[test]
fn a_meshed_block_converts_to_a_step_solid() {
    let mut model = ogeom::topo::Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 12.0, 8.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(10.0, 6.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drum = ogeom::algo::make_cylinder(&mut model, at, 3.0, 10.0, T)
        .unwrap()
        .shape;
    let bored = ogeom::boolean::cut(&mut model, &block, &drum, T)
        .unwrap()
        .shape;
    let mesh =
        ogeom::mesh::triangulate(&model, &bored, Deflection::with_chord(0.01).unwrap(), T).unwrap();

    let dir = std::env::temp_dir().join(format!("ogeom-cli-convert-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let stl = dir.join("bored.stl");
    let step = dir.join("bored.step");
    std::fs::write(
        &stl,
        ogeom::io::write(&mesh, ogeom::io::Encoding::Binary).unwrap(),
    )
    .unwrap();

    let run = std::process::Command::new(env!("CARGO_BIN_EXE_ogeom-cli"))
        .arg("convert")
        .arg(&stl)
        .arg("--step")
        .arg(&step)
        .output()
        .unwrap();
    let said = String::from_utf8_lossy(&run.stdout);
    assert!(
        run.status.success(),
        "{said}{}",
        String::from_utf8_lossy(&run.stderr)
    );
    assert!(
        said.contains("closed") && said.contains("  valid"),
        "{said}"
    );

    let import = ogeom::io::read_step(&std::fs::read_to_string(&step).unwrap(), T).unwrap();
    let back = import.document.model();
    let [solid] = import.solids.as_slice() else {
        panic!("{} solids", import.solids.len());
    };
    assert!(ogeom::algo::check(back, solid, T).unwrap().is_valid());
    let surfaces: Vec<_> = explore_unique(back, solid, ShapeType::Face)
        .unwrap()
        .iter()
        .map(|f| {
            let data = back.node(f).unwrap().data().as_face().unwrap();
            std::mem::discriminant(back.geometry().surface(data.surface).unwrap())
        })
        .collect();
    assert_eq!(surfaces.len(), 7);
    std::fs::remove_dir_all(&dir).unwrap();
}
