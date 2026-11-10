//! Colours named rather than mixed: `DRAUGHTING_PRE_DEFINED_COLOUR`.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::Frame;
use ogeom_topo::Model;

const T: Tolerances = Tolerances::millimetres();

/// A box written coloured, its `COLOUR_RGB` then replaced by a predefined
/// colour of the given name, as exporters write the eight colours.
fn box_coloured_by_name(name: &str) -> String {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T)
        .unwrap()
        .shape;
    let mut document = ogeom_doc::Document::over(model);
    document.add_part("part", solid.clone());
    document.set_colour(&solid, ogeom_doc::Colour::rgb(0.25, 0.5, 0.75));
    let text = ogeom_io::write_step(&document, T).unwrap();
    let start = text.find("COLOUR_RGB(").expect("the colour is written");
    let end = start + text[start..].find(')').unwrap() + 1;
    format!(
        "{}DRAUGHTING_PRE_DEFINED_COLOUR('{name}'){}",
        &text[..start],
        &text[end..]
    )
}

/// Each of the eight names reads as its colour, every channel fully on or
/// off; before, the item read uncoloured.
#[test]
fn the_eight_predefined_colours_read_as_their_colours() {
    for (name, rgb) in [
        ("black", (0.0, 0.0, 0.0)),
        ("red", (1.0, 0.0, 0.0)),
        ("green", (0.0, 1.0, 0.0)),
        ("blue", (0.0, 0.0, 1.0)),
        ("yellow", (1.0, 1.0, 0.0)),
        ("magenta", (1.0, 0.0, 1.0)),
        ("cyan", (0.0, 1.0, 1.0)),
        ("white", (1.0, 1.0, 1.0)),
    ] {
        let import = ogeom_io::read_step(&box_coloured_by_name(name), T).unwrap();
        let colour = import.document.colour_of(&import.solids[0]);
        assert_eq!(
            colour,
            Some(ogeom_doc::Colour::rgb(rgb.0, rgb.1, rgb.2)),
            "{name}"
        );
    }
}

/// A name that is none of the eight leaves the item uncoloured, with one
/// warning naming the entity.
#[test]
fn an_unknown_predefined_colour_reads_as_none_with_one_warning() {
    let text = box_coloured_by_name("mauve");
    let import = ogeom_io::read_step(&text, T).unwrap();
    assert_eq!(import.document.colour_of(&import.solids[0]), None);
    let named: Vec<&String> = import
        .report
        .warnings
        .iter()
        .filter(|w| w.contains("mauve"))
        .collect();
    assert_eq!(named.len(), 1, "{:?}", import.report.warnings);
    let at = text.find("DRAUGHTING_PRE_DEFINED_COLOUR(").unwrap();
    let line_start = text[..at].rfind('#').unwrap();
    let id = &text[line_start..text[line_start..].find('=').unwrap() + line_start];
    assert!(named[0].contains(id.trim()), "{} against {id}", named[0]);
}
