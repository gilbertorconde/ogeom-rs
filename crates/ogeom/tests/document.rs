//! §15's tail: a document that can be stepped back through, and textures
//! as the assignment they are.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::doc::{Colour, Document, Texture, TextureMapping};
use ogeom::math::Frame;

const T: Tolerances = Tolerances::millimetres();

#[test]
fn a_document_steps_back_and_forward_through_its_own_checkpoints() {
    let mut doc = Document::new();
    let block = ogeom::algo::make_box(doc.model_mut(), Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let part = doc.add_part("block", block.clone());
    doc.set_name(&block, "the block");

    // A checkpoint, then a change: the colour and a second name.
    doc.checkpoint();
    doc.set_colour(&block, Colour::rgb(1.0, 0.0, 0.0));
    doc.set_name(&block, "the red block");
    assert_eq!(doc.name_of(&block).unwrap(), "the red block");
    assert_eq!(doc.undo_depth(), (1, 0));

    // Back: the colour is gone and the name is the old one.
    assert!(doc.undo());
    assert_eq!(doc.name_of(&block).unwrap(), "the block");
    assert!(doc.colour_of(&block).is_none());
    assert_eq!(doc.undo_depth(), (0, 1));

    // Forward again: both return.
    assert!(doc.redo());
    assert_eq!(doc.name_of(&block).unwrap(), "the red block");
    assert!(doc.colour_of(&block).is_some());

    // Nothing further either way.
    assert!(doc.undo());
    assert!(!doc.undo());
    assert!(doc.redo());
    assert!(!doc.redo());

    // The product survives all of it — it was there before the first
    // checkpoint, so no step back reaches past it.
    assert!(doc.products().any(|(_, p)| p.name == "block"));
    let _ = part;
}

#[test]
fn a_texture_is_an_image_and_how_it_lands() {
    let mut doc = Document::new();
    let block = ogeom::algo::make_box(doc.model_mut(), Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;

    let brushed = doc.add_texture(Texture {
        image: "brushed-steel.png".into(),
        mapping: TextureMapping::Cylindrical,
        repeat: (4.0, 1.0),
        offset: (0.0, 0.0),
    });
    doc.set_texture(&block, brushed);

    let on = doc.texture_of(&block).expect("the block wears it");
    assert_eq!(on.image, "brushed-steel.png");
    assert_eq!(on.mapping, TextureMapping::Cylindrical);
    assert!((on.repeat.0 - 4.0).abs() < 1e-12);
    assert_eq!(doc.textures().len(), 1);

    // And it steps back with everything else.
    doc.checkpoint();
    let plain = doc.add_texture(Texture::image("paint.png"));
    doc.set_texture(&block, plain);
    assert_eq!(doc.texture_of(&block).unwrap().image, "paint.png");
    assert!(doc.undo());
    assert_eq!(doc.texture_of(&block).unwrap().image, "brushed-steel.png");
}
