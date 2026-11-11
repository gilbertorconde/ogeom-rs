//! A primitive command's `--step` writes the shape as a STEP part.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::mesh::Deflection;
use ogeom::topo::{ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// `box 4 3 2 --step` and `cylinder 2 3 --step` each write a file that
/// reads back as one valid solid of the primitive's faces; the box's
/// volume is its 24 exactly, its faces being planes.
#[test]
fn a_primitive_writes_the_step_part_it_was_asked_for() {
    let dir = std::env::temp_dir().join(format!("ogeom-cli-primitive-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    for (args, faces) in [
        (&["box", "4", "3", "2"][..], 6),
        (&["cylinder", "2", "3"][..], 3),
    ] {
        let step = dir.join(format!("{}.step", args[0]));
        assert!(!step.exists());
        let run = std::process::Command::new(env!("CARGO_BIN_EXE_ogeom-cli"))
            .args(args)
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
        assert!(said.contains("wrote"), "{said}");

        let import = ogeom::io::read_step(&std::fs::read_to_string(&step).unwrap(), T).unwrap();
        let back = import.document.model();
        let [solid] = import.solids.as_slice() else {
            panic!("{}: {} solids", args[0], import.solids.len());
        };
        assert!(ogeom::algo::check(back, solid, T).unwrap().is_valid());
        let found = explore_unique(back, solid, ShapeType::Face).unwrap().len();
        assert_eq!(found, faces, "{}", args[0]);
        if args[0] == "box" {
            let volume =
                ogeom::algo::volume_properties(back, solid, Deflection::default(), T).unwrap();
            assert!((volume.mass - 24.0).abs() < 1e-9, "{}", volume.mass);
        }
    }
    std::fs::remove_dir_all(&dir).unwrap();
}
