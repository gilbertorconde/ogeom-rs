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

/// A box whose solid is coloured one way, one face another and one edge a
/// third, written, each `COLOUR_RGB` then replaced by the predefined colour
/// named for it.
fn box_with_named_face_and_edge(face_name: &str, edge_name: &str) -> String {
    let mut model = Model::new();
    let solid = ogeom_algo::make_box(&mut model, Frame::WORLD, (1.0, 2.0, 3.0), T)
        .unwrap()
        .shape;
    let face =
        ogeom_topo::explore_unique(&model, &solid, ogeom_topo::ShapeType::Face).unwrap()[0].clone();
    let edge =
        ogeom_topo::explore_unique(&model, &solid, ogeom_topo::ShapeType::Edge).unwrap()[0].clone();
    let mut document = ogeom_doc::Document::over(model);
    document.add_part("part", solid.clone());
    document.set_colour(&solid, ogeom_doc::Colour::rgb(0.1, 0.2, 0.3));
    document.set_colour(&face, ogeom_doc::Colour::rgb(0.25, 0.5, 0.75));
    document.set_colour(&edge, ogeom_doc::Colour::rgb(0.6, 0.7, 0.8));
    let text = ogeom_io::write_step(&document, T).unwrap();
    let mut out = text.clone();
    for (rgb, name) in [("0.25,0.5,0.75", face_name), ("0.6,0.7,0.8", edge_name)] {
        let start = out
            .find(&format!("COLOUR_RGB('',{rgb})"))
            .expect("the colour is written");
        let end = start + out[start..].find(')').unwrap() + 1;
        out = format!(
            "{}DRAUGHTING_PRE_DEFINED_COLOUR('{name}'){}",
            &out[..start],
            &out[end..]
        );
    }
    out
}

/// The face's and edge's colours read from the styles over them.
fn face_and_edge_colours(
    text: &str,
) -> (
    Option<ogeom_doc::Colour>,
    Option<ogeom_doc::Colour>,
    Vec<String>,
) {
    let import = ogeom_io::read_step(text, T).unwrap();
    let model = import.document.model();
    let solid = &import.solids[0];
    let faces = ogeom_topo::explore_unique(model, solid, ogeom_topo::ShapeType::Face).unwrap();
    let edges = ogeom_topo::explore_unique(model, solid, ogeom_topo::ShapeType::Edge).unwrap();
    let face = faces.iter().find_map(|f| {
        import
            .document
            .colour_of(f)
            .filter(|c| *c != import.document.colour_of(solid).unwrap())
    });
    let edge = edges.iter().find_map(|e| import.document.colour_of(e));
    (face, edge, import.report.warnings)
}

/// A predefined colour is read wherever a style names it: under a face's
/// surface side style and an edge's curve style, as under a solid's. A name
/// in another case and with spaces about it is the same name, and none of
/// them warns.
#[test]
fn predefined_colours_read_on_faces_and_edges_in_any_case() {
    let (face, edge, warnings) =
        face_and_edge_colours(&box_with_named_face_and_edge("  Green  ", "RED"));
    assert_eq!(face, Some(ogeom_doc::Colour::rgb(0.0, 1.0, 0.0)));
    assert_eq!(edge, Some(ogeom_doc::Colour::rgb(1.0, 0.0, 0.0)));
    assert!(
        !warnings.iter().any(|w| w.contains("PRE_DEFINED")),
        "{warnings:?}"
    );
}

/// A face's style that overrides the solid's (`OVER_RIDING_STYLED_ITEM`)
/// reads its predefined colour as a plain styled item's.
#[test]
fn a_predefined_colour_in_an_overriding_style_reads() {
    let text = box_with_named_face_and_edge("blue", "cyan");
    // The face's styled item, written after the solid's, made to override
    // the solid's.
    let mut statements: Vec<String> = text.split(';').map(str::to_owned).collect();
    let styled: Vec<usize> = statements
        .iter()
        .enumerate()
        .filter(|(_, s)| s.contains("=STYLED_ITEM("))
        .map(|(i, _)| i)
        .collect();
    let id = |s: &str| s.trim().split('=').next().unwrap().trim().to_owned();
    let solid_item = id(&statements[styled[0]]);
    let face_at = styled[1];
    let face = statements[face_at].replace("=STYLED_ITEM(", "=OVER_RIDING_STYLED_ITEM(");
    let close = face.rfind(')').unwrap();
    statements[face_at] = format!("{},{solid_item}{}", &face[..close], &face[close..]);
    let (face, _, _) = face_and_edge_colours(&statements.join(";"));
    assert_eq!(face, Some(ogeom_doc::Colour::rgb(0.0, 0.0, 1.0)));
}
