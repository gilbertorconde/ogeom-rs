//! M5's first criterion takes shape: a multi-body assembly read with its
//! product structure, placements, names and colours intact.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::doc::{Colour, ProductKind};
use ogeom::math::Point;
use ogeom::topo::{Filter, ShapeType, explore};

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

#[test]
fn the_bolted_plate_assembly_reads_whole() {
    let text = corpus("ogeom_asm_bolted_plate.stp");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let doc = &import.document;
    assert!(
        import.report.warnings.is_empty(),
        "{:?}",
        import.report.warnings
    );

    // Three products: two parts and the assembly, named as authored.
    let names: Vec<&str> = doc.products().map(|(_, p)| p.name.as_str()).collect();
    assert_eq!(names, ["plate", "bolt", "bolted-plate"]);
    let roots = doc.roots();
    assert_eq!(roots.len(), 1, "one root: the assembly");
    let root = roots[0];
    assert!(matches!(
        doc.get(root).unwrap().kind,
        ProductKind::Assembly { .. }
    ));

    // Three occurrences, named by their reference designators, placed where
    // the file says.
    let mut occurrences = doc.occurrences_of(root).unwrap();
    occurrences.sort_by(|a, b| a.path.cmp(&b.path));
    let paths: Vec<&str> = occurrences.iter().map(|o| o.path.as_str()).collect();
    assert_eq!(
        paths,
        [
            "bolted-plate/bolt-1",
            "bolted-plate/bolt-2",
            "bolted-plate/plate-1"
        ]
    );
    let world = |i: usize| {
        occurrences[i]
            .shape
            .transform(doc.model().datums())
            .unwrap()
            .apply(Point::new(0.0, 0.0, 0.0))
    };
    assert!(world(0).is_equal(Point::new(10.0, 10.0, 5.0), T));
    assert!(world(1).is_equal(Point::new(30.0, 30.0, 5.0), T));
    assert!(world(2).is_equal(Point::new(0.0, 0.0, 0.0), T));

    // The two bolts are one part: same node, different chains — instancing.
    assert_eq!(occurrences[0].shape.node(), occurrences[1].shape.node());
    assert_eq!(occurrences[0].part, occurrences[1].part);

    // Volumes survive placement: the flattened assembly measures as the sum.
    let fine = ogeom::mesh::Deflection::default();
    let mut total = 0.0;
    for occ in &occurrences {
        total += ogeom::algo::volume_properties(doc.model(), &occ.shape, fine, T)
            .unwrap()
            .mass;
    }
    let expected = 40.0 * 40.0 * 5.0 + 2.0 * (8.0 * 8.0 * 20.0);
    assert!(
        (total - expected).abs() < 1.0,
        "assembly volume {total} against {expected}"
    );

    // Colours: the plate green, the bolt blue — resolved per occurrence.
    let colour = |i: usize| doc.resolved_colour(occurrences[i].part, &occurrences[i].shape);
    assert_eq!(colour(0), Some(Colour::rgb(0.2, 0.3, 0.9)));
    assert_eq!(colour(2), Some(Colour::rgb(0.1, 0.8, 0.2)));

    // And one face of the bolt carries its own overriding red.
    let bolt_faces = explore(
        doc.model(),
        &occurrences[0].shape,
        Filter::OfType(ShapeType::Face),
    )
    .unwrap();
    let red: Vec<Colour> = bolt_faces.iter().filter_map(|f| doc.colour_of(f)).collect();
    assert_eq!(red, [Colour::rgb(0.9, 0.1, 0.1)]);
}

#[test]
fn a_single_part_file_still_yields_a_product_with_its_colours() {
    // The NIST parts have no assembly, but they do carry product names and
    // styled colours; the document says so instead of losing them.
    let text = corpus("nist_ftc_11_asme1_rb.stp");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let doc = &import.document;
    let products: Vec<_> = doc.products().collect();
    assert_eq!(products.len(), 1, "one part product");
    assert!(matches!(products[0].1.kind, ProductKind::Part { .. }));
    let coloured = explore(
        doc.model(),
        &import.solids[0],
        Filter::OfType(ShapeType::Solid),
    )
    .unwrap()
    .iter()
    .filter_map(|s| doc.colour_of(s))
    .count();
    let solid_coloured = doc.colour_of(&import.solids[0]).is_some();
    assert!(
        coloured > 0 || solid_coloured,
        "the file styles its solid; the document records it"
    );
}

#[test]
fn the_assembly_round_trips_through_the_writer() {
    let text = corpus("ogeom_asm_bolted_plate.stp");
    let first = ogeom::io::read_step(&text, T).unwrap();
    let written = ogeom::io::write_step(&first.document, T).unwrap();
    let second = ogeom::io::read_step(&written, T).unwrap();
    assert!(
        second.report.warnings.is_empty(),
        "{:?}",
        second.report.warnings
    );

    // Same products, same names, same single root.
    let names = |doc: &ogeom::doc::Document| -> Vec<String> {
        doc.products().map(|(_, p)| p.name.clone()).collect()
    };
    assert_eq!(names(&first.document), names(&second.document));
    let root = second.document.roots();
    assert_eq!(root.len(), 1);

    // Same occurrences at the same places.
    let flatten = |doc: &ogeom::doc::Document| -> Vec<(String, Point)> {
        let root = doc.roots()[0];
        let mut out: Vec<(String, Point)> = doc
            .occurrences_of(root)
            .unwrap()
            .iter()
            .map(|o| {
                (
                    o.path.clone(),
                    o.shape
                        .transform(doc.model().datums())
                        .unwrap()
                        .apply(Point::new(0.0, 0.0, 0.0)),
                )
            })
            .collect();
        out.sort_by(|a, b| a.0.cmp(&b.0));
        out
    };
    let (a, b) = (flatten(&first.document), flatten(&second.document));
    assert_eq!(a.len(), b.len());
    for ((pa, wa), (pb, wb)) in a.iter().zip(&b) {
        assert_eq!(pa, pb);
        assert!(wa.is_equal(*wb, T), "{pa}: {wa:?} vs {wb:?}");
    }

    // Same material: the flattened volumes agree.
    let fine = ogeom::mesh::Deflection::default();
    let total = |doc: &ogeom::doc::Document| -> f64 {
        let root = doc.roots()[0];
        doc.occurrences_of(root)
            .unwrap()
            .iter()
            .map(|o| {
                ogeom::algo::volume_properties(doc.model(), &o.shape, fine, T)
                    .unwrap()
                    .mass
            })
            .sum()
    };
    let (va, vb) = (total(&first.document), total(&second.document));
    assert!((va - vb).abs() < 1e-6, "{va} vs {vb}");

    // Same colours, resolved per occurrence, and the red face survives.
    let root2 = second.document.roots()[0];
    let mut occurrences = second.document.occurrences_of(root2).unwrap();
    occurrences.sort_by(|a, b| a.path.cmp(&b.path));
    assert_eq!(
        second
            .document
            .resolved_colour(occurrences[0].part, &occurrences[0].shape),
        Some(Colour::rgb(0.2, 0.3, 0.9))
    );
    let red = explore(
        second.document.model(),
        &occurrences[0].shape,
        Filter::OfType(ShapeType::Face),
    )
    .unwrap()
    .iter()
    .filter_map(|f| second.document.colour_of(f))
    .count();
    assert_eq!(red, 1, "the overriding face colour survives the round trip");
}

#[test]
fn a_real_part_round_trips_with_its_volume() {
    // ftc_11 is the all-analytic NIST part: planes, cylinders, seams,
    // misaligned rings and all. What the writer serializes, the reader must
    // heal and measure identically — the same healing the original needs,
    // because the writer writes what the file actually says.
    let fine = ogeom::mesh::Deflection {
        chord: 1e-2,
        ..ogeom::mesh::Deflection::default()
    };
    let healed_volume = |mut import: ogeom::io::StepImport| -> f64 {
        let solid = import.solids[0].clone();
        let healed = ogeom::heal::reanchor_periodic_rings(import.document.model_mut(), &solid, T)
            .map_or(solid, |h| h.0.shape);
        ogeom::algo::volume_properties(import.document.model(), &healed, fine, T)
            .unwrap()
            .mass
    };

    let text = corpus("nist_ftc_11_asme1_rb.stp");
    let first = ogeom::io::read_step(&text, T).unwrap();
    let written = ogeom::io::write_step(&first.document, T).unwrap();
    let before = healed_volume(first);
    let second = ogeom::io::read_step(&written, T).unwrap();
    let after = healed_volume(second);
    assert!(
        (before - after).abs() < before * 1e-6,
        "volume through the round trip: {before} -> {after}"
    );
}
