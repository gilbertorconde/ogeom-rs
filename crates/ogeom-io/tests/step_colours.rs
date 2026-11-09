//! Colours read from STEP: a solid styled by a `COLOUR_RGB`, and the same
//! file with that colour named instead by a `DRAUGHTING_PRE_DEFINED_COLOUR`,
//! the form exporters write for the eight colours the presentation schema
//! defines by name.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_doc::{Colour, Document};
use ogeom_math::Frame;
use ogeom_topo::Model;

const T: Tolerances = Tolerances::millimetres();

/// A box coloured with a colour none of the predefined ones equals, as
/// STEP text, with the one `COLOUR_RGB(...)` it carries.
fn styled_box() -> (String, String) {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (3.0, 2.0, 1.0), T)
        .unwrap()
        .shape;
    let mut document = Document::over(model);
    document.add_part("box", solid.clone());
    document.set_colour(&solid, Colour::rgb(0.25, 0.5, 0.75));
    let text = ogeom_io::write_step(&document, T).unwrap();
    let mut found = text.match_indices("COLOUR_RGB(");
    let (at, _) = found.next().expect("the styled box carries a colour");
    assert!(found.next().is_none(), "one colour in the file");
    let end = at + text[at..].find(')').unwrap() + 1;
    let rgb = text[at..end].to_owned();
    (text, rgb)
}

/// The colour the read gives the file's one solid, and the read's warnings.
fn read_colour(text: &str) -> (Option<Colour>, Vec<String>) {
    let import = ogeom_io::read_step(text, T).unwrap();
    assert_eq!(import.solids.len(), 1);
    (
        import.document.colour_of(&import.solids[0]),
        import.report.warnings,
    )
}

#[test]
fn an_rgb_colour_reads_back() {
    let (text, _) = styled_box();
    let (colour, _) = read_colour(&text);
    assert_eq!(colour, Some(Colour::rgb(0.25, 0.5, 0.75)));
}

#[test]
fn each_predefined_colour_reads_as_its_name_says() {
    let (text, rgb) = styled_box();
    let named = [
        ("red", (1.0, 0.0, 0.0)),
        ("green", (0.0, 1.0, 0.0)),
        ("blue", (0.0, 0.0, 1.0)),
        ("yellow", (1.0, 1.0, 0.0)),
        ("magenta", (1.0, 0.0, 1.0)),
        ("cyan", (0.0, 1.0, 1.0)),
        ("black", (0.0, 0.0, 0.0)),
        ("white", (1.0, 1.0, 1.0)),
    ];
    for (name, (r, g, b)) in named {
        let text = text.replace(&rgb, &format!("DRAUGHTING_PRE_DEFINED_COLOUR('{name}')"));
        let (colour, warnings) = read_colour(&text);
        assert_eq!(colour, Some(Colour::rgb(r, g, b)), "{name}");
        assert!(warnings.is_empty(), "{name}: {warnings:?}");
    }
}

#[test]
fn an_unknown_predefined_colour_leaves_the_solid_uncoloured_and_says_so() {
    let (text, rgb) = styled_box();
    let text = text.replace(&rgb, "DRAUGHTING_PRE_DEFINED_COLOUR('mauve')");
    let (colour, warnings) = read_colour(&text);
    assert_eq!(colour, None);
    let about: Vec<&String> = warnings.iter().filter(|w| w.contains("'mauve'")).collect();
    assert_eq!(about.len(), 1, "{warnings:?}");
}
