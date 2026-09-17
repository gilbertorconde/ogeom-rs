//! M3's closing argument: the corpus imported, healed, measured, and
//! operated on — real files in, real modelling out.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::topo::{Filter, ShapeType, explore, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

/// Healing sweeps the whole corpus: every part reads, every shell closes,
/// every part measures, and each volume pins to its known figure — the
/// values ftc_07 and ftc_11 arbitrated against the kernel's exact ray
/// classifier and ctc_01 against an orientation-free even-odd grid. A loose
/// relative band absorbs future mesh refinements; an orientation or unit
/// mistake moves a volume by whole factors and cannot hide inside it.
#[test]
fn the_corpus_heals_and_every_part_measures() {
    let files: [(&str, f64); 11] = [
        ("nist_ctc_01_asme1_rd.stp", 14_643_073.4),
        ("nist_ctc_02_asme1_rc.stp", 47_101_909.1),
        ("nist_ctc_03_asme1_rc.stp", 331_884.8),
        ("nist_ctc_04_asme1_rd.stp", 17_519_472.2),
        ("nist_ctc_05_asme1_rd.stp", 12_694_721.1),
        ("nist_ftc_06_asme1_rd.stp", 3_289_989.2),
        ("nist_ftc_07_asme1_rd.stp", 1_726_286.9),
        ("nist_ftc_08_asme1_rc.stp", 503_596.4),
        ("nist_ftc_09_asme1_rd.stp", 136_453.8),
        ("nist_ftc_10_asme1_rb.stp", 188_388.8),
        ("nist_ftc_11_asme1_rb.stp", 5_122.3),
    ];
    let fine = ogeom::mesh::Deflection {
        chord: 1e-2,
        ..ogeom::mesh::Deflection::default()
    };
    for (name, expected) in files {
        let text = corpus(name);
        let mut import = ogeom::io::read_step(&text, T).unwrap();
        let solid = import.solids[0].clone();
        // Healing is idempotent where nothing is broken and surgery where
        // something is; either way the shell must close.
        let healed = ogeom::heal::reanchor_periodic_rings(import.document.model_mut(), &solid, T)
            .map_or(solid, |h| h.0.shape);
        let shell = explore_unique(import.document.model(), &healed, ShapeType::Shell)
            .unwrap()
            .remove(0);
        assert!(
            ogeom::algo::is_shell_closed(import.document.model(), &shell).unwrap(),
            "{name}: the healed shell closes"
        );
        let props = ogeom::algo::volume_properties(import.document.model(), &healed, fine, T)
            .unwrap_or_else(|e| panic!("{name}: does not measure: {e}"));
        eprintln!("REPORT {name}: volume {:.3} mm^3", props.mass);
        assert!(
            (props.mass - expected).abs() <= expected * 1e-2,
            "{name}: volume {:.3} strays from its pinned {expected:.1}",
            props.mass
        );
    }
}

/// A boolean over imported geometry, history checked: the milestone's
/// operations run on the world's parts, not only on this kernel's own.
#[test]
fn an_imported_part_takes_a_boolean_cut() {
    let text = corpus("nist_ftc_11_asme1_rb.stp");
    let mut import = ogeom::io::read_step(&text, T).unwrap();
    let solid = import.solids[0].clone();
    let healed = ogeom::heal::reanchor_periodic_rings(import.document.model_mut(), &solid, T)
        .unwrap()
        .0
        .shape;
    let fine = ogeom::mesh::Deflection {
        chord: 1e-2,
        ..ogeom::mesh::Deflection::default()
    };
    let before = ogeom::algo::volume_properties(import.document.model(), &healed, fine, T)
        .unwrap()
        .mass;

    // A square post cut down through the plate's solid ring — the part has a
    // large central pocket, and a post through fresh air cuts nothing, as an
    // earlier version of this test discovered only once the result's mesh
    // first became measurable.
    let frame = ogeom::math::Frame::new(
        ogeom::math::Point::new(20.0, -4.0, -3.0),
        ogeom::math::Direction::Z,
        ogeom::math::Direction::X,
        T,
    )
    .unwrap();
    let post =
        ogeom::algo::make_box(import.document.model_mut(), frame, (8.0, 8.0, 6.0), T).unwrap();
    let result = ogeom::boolean::cut(import.document.model_mut(), &healed, &post.shape, T).unwrap();

    // The cut runs on the imported part: the result is a solid whose shell
    // closes, built from pieces of the world's geometry and this kernel's,
    // and it *measures* — the mesh welds across the file's slop because the
    // slop is recorded on the edges and vertices and the weld honours it.
    let shell = explore_unique(import.document.model(), &result.shape, ShapeType::Shell)
        .unwrap()
        .remove(0);
    assert!(ogeom::algo::is_shell_closed(import.document.model(), &shell).unwrap());
    let props =
        ogeom::algo::volume_properties(import.document.model(), &result.shape, fine, T).unwrap();
    eprintln!("REPORT cut imported: {before:.3} -> {:.3} mm^3", props.mass);
    assert!(props.mass < before, "the cut removed material");
    assert!(props.mass > 0.0);
    // The post removes at most its own volume, and the difference cannot
    // exceed it: the bound the arithmetic itself provides.
    assert!(
        before - props.mass <= 8.0 * 8.0 * 6.0,
        "the cut removed more than the post could"
    );

    // History carried through: the imported solid is recorded as modified
    // into the result, and at least one of its faces was split or consumed.
    assert_eq!(
        result.history.modified(&healed),
        std::slice::from_ref(&result.shape)
    );
    let touched = explore(
        import.document.model(),
        &healed,
        Filter::OfType(ShapeType::Face),
    )
    .unwrap()
    .iter()
    .filter(|f| result.history.is_affected(f))
    .count();
    assert!(touched > 0, "faces of the imported part appear in history");
}

/// A body is bounded by what it is trimmed to, not by the carriers under it.
///
/// `shape_bounds` promised a guarantee and delivered the carriers: an
/// imported plane reported its own window, which spans kilometres, and a
/// cylinder its height domain. A consumer asking how big a screw is got a
/// box a billion millimetres across and had to fall back to the hull of the
/// topological vertices — which does not contain the body either, since a
/// button head's apex is a bulge between its rims, three millimetres past
/// every vertex the head has.
///
/// Two things were wrong. An edge reported the bound of its whole curve
/// rather than of the range it uses, so a segment on a line that runs to
/// the ends of the world reported the ends of the world. And a curved face
/// reported its whole surface, when its boundary is bounded by the edges
/// below it and all the face itself has to add is where its surface bulges
/// past that boundary.
#[test]
fn a_body_is_bounded_by_what_it_is_trimmed_to() {
    let deflection = ogeom::mesh::Deflection::default();
    // The same screw in two encodings, five bodies in one and one in the
    // other, and the second is where the figures below come from: an M5×16
    // button head, 18.75 long over a head 9.5 across.
    for (file, bodies) in [("m5x16_bhcs.step", 5), ("m5x16_bhcs_loops.step", 1)] {
        let text = corpus(file);
        let import = ogeom::io::read_step(&text, T).unwrap();
        let model = import.document.model();
        assert_eq!(import.solids.len(), bodies, "{file}: its bodies");
        for (which, solid) in import.solids.iter().enumerate() {
            let bounds = ogeom::algo::shape_bounds(model, solid, T).unwrap();
            let mesh = ogeom::mesh::triangulate(model, solid, deflection, T).unwrap();
            let mut hull = ogeom::math::Aabb::default();
            for at in &mesh.positions {
                hull = hull.with_point(*at);
            }
            assert!(
                bounds.contains_box(&hull),
                "{file} body {which}: the bound holds the mesh, {:?}..{:?} against {:?}..{:?}",
                bounds.low(),
                bounds.high(),
                hull.low(),
                hull.high()
            );
            // And holds little else. Where the bound stands proud of the mesh it
            // is the mesh that is short — a chord across an arc falls inside it
            // — so the margin is the chord's, either way round.
            let (low, high) = (bounds.low().unwrap(), bounds.high().unwrap());
            let (near, far) = (hull.low().unwrap(), hull.high().unwrap());
            for (bound, meshed) in [
                (low.x, near.x),
                (low.y, near.y),
                (low.z, near.z),
                (high.x, far.x),
                (high.y, far.y),
                (high.z, far.z),
            ] {
                assert!(
                    (bound - meshed).abs() <= deflection.chord,
                    "{file} body {which}: the bound stands {} from the mesh",
                    (bound - meshed).abs()
                );
            }
            // The screw itself, and not the carriers under it.
            let size = high - low;
            assert!(
                size.x < 20.0 && (size.y - 9.5).abs() < 0.1 && (size.z - 9.5).abs() < 0.1,
                "{file} body {which}: an M5×16 button head, {size:?}"
            );
        }
    }
}
