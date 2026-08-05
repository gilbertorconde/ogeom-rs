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
