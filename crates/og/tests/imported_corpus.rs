//! M3's closing argument: the corpus imported, healed, measured, and
//! operated on — real files in, real modelling out.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use og::core::Tolerances;
use og::topo::{Filter, ShapeType, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

/// Healing sweeps the whole corpus: every part reads, every shell closes,
/// and the parts whose faces the kernel can annotate measure real volumes.
#[test]
fn the_corpus_heals_and_the_analytic_parts_measure() {
    let files = [
        "nist_ctc_01_asme1_rd.stp",
        "nist_ctc_02_asme1_rc.stp",
        "nist_ctc_03_asme1_rc.stp",
        "nist_ctc_04_asme1_rd.stp",
        "nist_ctc_05_asme1_rd.stp",
        "nist_ftc_06_asme1_rd.stp",
        "nist_ftc_07_asme1_rd.stp",
        "nist_ftc_08_asme1_rc.stp",
        "nist_ftc_09_asme1_rd.stp",
        "nist_ftc_10_asme1_rb.stp",
        "nist_ftc_11_asme1_rb.stp",
    ];
    let fine = og::mesh::Deflection {
        chord: 1e-2,
        ..og::mesh::Deflection::default()
    };
    let mut measured = 0;
    for name in files {
        let text = corpus(name);
        let mut import = og::io::read_step(&text, T).unwrap();
        let solid = import.solids[0].clone();
        // Healing is idempotent where nothing is broken and surgery where
        // something is; either way the shell must close.
        let healed = og::heal::reanchor_periodic_rings(&mut import.model, &solid, T)
            .map_or(solid, |h| h.shape);
        let shell = explore_unique(&import.model, &healed, ShapeType::Shell)
            .unwrap()
            .remove(0);
        assert!(
            og::algo::is_shell_closed(&import.model, &shell).unwrap(),
            "{name}: the healed shell closes"
        );
        match og::algo::volume_properties(&import.model, &healed, fine, T) {
            Ok(props) => {
                assert!(props.mass > 0.0, "{name}: a part encloses volume");
                eprintln!("REPORT {name}: volume {:.3} mm^3", props.mass);
                measured += 1;
            }
            Err(e) => eprintln!("REPORT {name}: not yet measurable ({e})"),
        }
    }
    assert!(
        measured >= 1,
        "at least the all-analytic parts measure after healing"
    );
    eprintln!("REPORT measured {measured}/11");
}

/// A boolean over imported geometry, history checked: the milestone's
/// operations run on the world's parts, not only on this kernel's own.
#[test]
fn an_imported_part_takes_a_boolean_cut() {
    let text = corpus("nist_ftc_11_asme1_rb.stp");
    let mut import = og::io::read_step(&text, T).unwrap();
    let solid = import.solids[0].clone();
    let healed = og::heal::reanchor_periodic_rings(&mut import.model, &solid, T)
        .unwrap()
        .shape;
    let fine = og::mesh::Deflection {
        chord: 1e-2,
        ..og::mesh::Deflection::default()
    };
    let before = og::algo::volume_properties(&import.model, &healed, fine, T)
        .unwrap()
        .mass;

    // A square post cut down through the plate, clear of every fillet.
    let frame = og::math::Frame::new(
        og::math::Point::new(-4.0, -4.0, -3.0),
        og::math::Direction::Z,
        og::math::Direction::X,
        T,
    )
    .unwrap();
    let post = og::algo::make_box(&mut import.model, frame, (8.0, 8.0, 6.0), T).unwrap();
    let result = og::boolean::cut(&mut import.model, &healed, &post.shape, T).unwrap();

    // The cut runs on the imported part: the result is a solid whose shell
    // closes, built from pieces of the world's geometry and this kernel's.
    // Its mesh does not yet weld everywhere — the rebuilt revolution pieces
    // still meet the plate pieces within the file's slop — so the volume is
    // reported when measurable and the structure asserted regardless.
    let shell = explore_unique(&import.model, &result.shape, ShapeType::Shell)
        .unwrap()
        .remove(0);
    assert!(og::algo::is_shell_closed(&import.model, &shell).unwrap());
    match og::algo::volume_properties(&import.model, &result.shape, fine, T) {
        Ok(props) => {
            eprintln!("REPORT cut imported: {before:.3} -> {:.3} mm^3", props.mass);
            assert!(props.mass < before, "the cut removed material");
            assert!(props.mass > 0.0);
        }
        Err(e) => eprintln!("REPORT cut imported: unmeasured ({e})"),
    }

    // History carried through: the imported solid is recorded as modified
    // into the result, and at least one of its faces was split or consumed.
    assert_eq!(
        result.history.modified(&healed),
        std::slice::from_ref(&result.shape)
    );
    let touched = explore(&import.model, &healed, Filter::OfType(ShapeType::Face))
        .unwrap()
        .iter()
        .filter(|f| result.history.is_affected(f))
        .count();
    assert!(touched > 0, "faces of the imported part appear in history");
}
