//! AP242 semantic PMI: the values, not the leader lines, read from NIST's
//! own annotated part and carried through the writer intact.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::doc::{MeasureKind, Pmi};

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

/// One dimension, flattened for comparison.
type DimensionRow = (String, Vec<String>, Option<String>, usize);
/// One tolerance, flattened for comparison.
type ToleranceRow = (String, String, String, Vec<String>, usize);
/// One datum, flattened for comparison.
type DatumRow = (String, usize);

/// The comparable core of a PMI set: everything except node identities,
/// which are model-scoped and change across a round trip.
fn summary(pmi: &Pmi) -> (Vec<DimensionRow>, Vec<ToleranceRow>, Vec<DatumRow>) {
    let f = |v: f64| format!("{v:.9}");
    let mut dims: Vec<_> = pmi
        .dimensions
        .iter()
        .map(|d| {
            (
                format!("{}|{:?}|{}", d.name, d.kind, d.location),
                d.values.iter().map(|&v| f(v)).collect::<Vec<_>>(),
                d.plus
                    .map(f)
                    .zip(d.minus.map(f))
                    .map(|(p, m)| format!("{p}/{m}")),
                d.items().count(),
            )
        })
        .collect();
    dims.sort();
    let mut tols: Vec<_> = pmi
        .tolerances
        .iter()
        .map(|t| {
            (
                t.kind.clone(),
                t.name.clone(),
                f(t.magnitude),
                t.datums
                    .iter()
                    .map(|d| format!("{}{:?}", d.label, d.modifiers))
                    .collect(),
                t.items.len(),
            )
        })
        .collect();
    tols.sort();
    let mut datums: Vec<_> = pmi
        .datums
        .iter()
        .map(|d| (d.label.clone(), d.items.len()))
        .collect();
    datums.sort();
    (dims, tols, datums)
}

#[test]
fn the_nist_part_reads_its_annotations_semantically() {
    let text = corpus("nist_ctc_01_asme1_ap242-e1.stp");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let pmi = import.document.pmi();

    // Eight dimensional characteristics, values in millimetres.
    assert_eq!(pmi.dimensions.len(), 8);
    let diameters: Vec<f64> = pmi
        .dimensions
        .iter()
        .filter(|d| d.name == "diameter")
        .map(|d| d.values[0])
        .collect();
    assert_eq!(diameters, [35.0, 35.0, 20.0, 20.0, 35.0, 35.0, 25.0]);
    // The one angular dimension: sixty degrees, plus/minus half a degree,
    // radians throughout.
    let angle = pmi
        .dimensions
        .iter()
        .find(|d| d.kind == MeasureKind::Angle)
        .expect("the part dimensions an angle");
    assert!((angle.values[0] - 60_f64.to_radians()).abs() < 1e-12);
    assert!((angle.plus.unwrap() - 0.5_f64.to_radians()).abs() < 1e-12);

    // Six geometric tolerances, magnitudes and datum reference frames as
    // NIST published them: the file states each frame as a datum system,
    // `A` alone or `A|B|C` in precedence order.
    let mut kinds: Vec<(&str, f64, Vec<&str>, bool)> = pmi
        .tolerances
        .iter()
        .map(|t| {
            let labels = t.datums.iter().map(|d| d.label.as_str()).collect();
            (t.kind.as_str(), t.magnitude, labels, t.datum_system)
        })
        .collect();
    kinds.sort_by(|a, b| a.0.cmp(b.0).then(a.1.total_cmp(&b.1)));
    let a = vec!["A"];
    let abc = vec!["A", "B", "C"];
    assert_eq!(
        kinds,
        [
            ("flatness", 0.2, vec![], false),
            ("perpendicularity", 1.5, a.clone(), true),
            ("position", 0.75, abc.clone(), true),
            ("position", 0.75, abc.clone(), true),
            ("surface_profile", 0.5, a, true),
            ("surface_profile", 1.25, abc, true),
        ]
    );

    // Three datums: A, B, C.
    let labels: Vec<&str> = pmi.datums.iter().map(|d| d.label.as_str()).collect();
    assert_eq!(labels, ["A", "B", "C"]);

    // The deeper aspect walk resolves the angular dimension to real
    // topology: its aspects sit three relationship steps out.
    let angle = pmi
        .dimensions
        .iter()
        .find(|d| d.kind == MeasureKind::Angle)
        .unwrap();
    assert!(
        angle.items().count() > 0,
        "the angle dimension names the faces it measures"
    );
}

/// Material-condition modifiers and a composite datum reference survive
/// the writer: the modifiers force the complex tolerance form, the
/// hyphen-joined label a compartment binding two datums into one, and the
/// reader hands both back as they went in.
#[test]
fn modifiers_and_composite_datums_survive_the_writer() {
    let text = corpus("nist_ctc_01_asme1_ap242-e1.stp");
    let mut first = ogeom::io::read_step(&text, T).unwrap();
    {
        let tolerance = first
            .document
            .pmi_mut()
            .tolerances
            .iter_mut()
            .find(|t| t.kind == "perpendicularity")
            .expect("the part carries a perpendicularity");
        tolerance.modifiers = vec![
            "maximum_material_requirement".into(),
            "projected_tolerance_zone".into(),
        ];
        tolerance.datums = vec!["A-B".into(), "C".into()];
    }
    let written = ogeom::io::write_step(&first.document, T).unwrap();
    let second = ogeom::io::read_step(&written, T).unwrap();
    let tolerance = second
        .document
        .pmi()
        .tolerances
        .iter()
        .find(|t| t.kind == "perpendicularity")
        .expect("the perpendicularity survives");
    assert_eq!(
        tolerance.modifiers,
        ["maximum_material_requirement", "projected_tolerance_zone"]
    );
    assert_eq!(tolerance.datums, ["A-B", "C"]);
    assert!((tolerance.magnitude - 1.5).abs() < 1e-12);
}

#[test]
fn pmi_round_trips_through_the_writer() {
    let text = corpus("nist_ctc_01_asme1_ap242-e1.stp");
    let first = ogeom::io::read_step(&text, T).unwrap();
    let written = ogeom::io::write_step(&first.document, T).unwrap();
    let second = ogeom::io::read_step(&written, T).unwrap();
    assert_eq!(
        summary(first.document.pmi()),
        summary(second.document.pmi())
    );
}

/// The *presentation* half. NIST's part draws its annotations, and what
/// it draws is read: the callouts, the plane each is drawn in, and the
/// polylines that make the frame, the leader and the strokes of the text.
#[test]
fn the_drawn_annotations_are_read_with_the_annotations_they_draw() {
    let text = corpus("nist_ctc_01_asme1_ap242-e1.stp");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let pmi = import.document.pmi();

    // Twenty-three callouts, every one drawn in a stated plane and every one
    // carrying geometry.
    assert_eq!(pmi.callouts.len(), 23);
    for callout in &pmi.callouts {
        assert!(callout.plane.is_some(), "{} has a plane", callout.name);
        assert!(
            !callout.polylines.is_empty(),
            "{} draws something",
            callout.name
        );
        for polyline in &callout.polylines {
            assert!(polyline.len() >= 2, "a polyline is at least a segment");
        }
    }

    // A callout says which semantic annotation it is a picture of, where it
    // has one. Not every one does, and that is the file rather than the
    // reader: a text note is presentation with nothing semantic behind it,
    // and a size the file draws without a characteristic representation
    // carries no value to have recorded. Fourteen of the twenty-three link,
    // which is a measurement of this file and will say so if it changes.
    let linked = pmi
        .callouts
        .iter()
        .filter_map(|c| c.annotates)
        .collect::<Vec<_>>();
    assert_eq!(linked.len(), 14, "what this file links");
    for callout in &pmi.callouts {
        if callout.name.starts_with("Text") {
            assert!(
                callout.annotates.is_none(),
                "a note draws nothing semantic: {}",
                callout.name
            );
        }
    }
    for what in &linked {
        match what {
            ogeom::doc::Annotated::Dimension(i) => assert!(*i < pmi.dimensions.len()),
            ogeom::doc::Annotated::Tolerance(i) => assert!(*i < pmi.tolerances.len()),
            ogeom::doc::Annotated::Datum(i) => assert!(*i < pmi.datums.len()),
        }
    }
    // The flatness callout draws the flatness.
    let flatness = pmi
        .callouts
        .iter()
        .find(|c| c.name.starts_with("Flatness"))
        .expect("the part draws its flatness");
    let ogeom::doc::Annotated::Tolerance(at) = flatness.annotates.unwrap() else {
        panic!("a flatness callout draws a tolerance");
    };
    assert_eq!(pmi.tolerances[at].kind, "flatness");

    // What a callout draws is *flat*, in a plane parallel to the one the
    // annotation is anchored in. This file offsets the drawing from that
    // plane, which is its own arrangement and not something to assume away.
    // The claim is the one that checks the reader's work: the item's own
    // repositioning was applied, so every point of one callout shares one
    // offset along the plane's normal. Get the repositioning wrong and the
    // offsets scatter.
    for callout in &pmi.callouts {
        let plane = callout.plane.expect("every callout stated its plane");
        let normal = plane.z().vector();
        let offsets: Vec<f64> = callout
            .polylines
            .iter()
            .flatten()
            .map(|p| (*p - plane.origin()).dot(normal))
            .collect();
        assert!(offsets.iter().all(|d| d.is_finite()));
        let first = offsets[0];
        let together = offsets
            .iter()
            .filter(|d| (**d - first).abs() < 1e-6)
            .count();
        assert!(
            together == offsets.len(),
            "{}: {together} of {} points share one offset from its plane",
            callout.name,
            offsets.len()
        );
    }
}

/// The presentation and the datum targets survive the writer, which is the
/// only thing that says the carriage is real rather than a reading.
#[test]
fn presentation_and_datum_targets_round_trip() {
    let text = corpus("nist_ctc_01_asme1_ap242-e1.stp");
    let mut first = ogeom::io::read_step(&text, T).unwrap();

    // NIST's part carries no datum targets, so the round trip is given some:
    // one of each kind, on datum A, placed where the drawing would put them.
    {
        let items = first.document.pmi().datums[0].items.clone();
        let targets = &mut first.document.pmi_mut().targets;
        targets.push(ogeom::doc::DatumTarget {
            datum: "A".into(),
            index: 1,
            kind: ogeom::doc::DatumTargetKind::Point,
            at: ogeom::math::Point::new(10.0, 20.0, 30.0),
            frame: None,
            items: items.clone(),
        });
        targets.push(ogeom::doc::DatumTarget {
            datum: "A".into(),
            index: 2,
            kind: ogeom::doc::DatumTargetKind::Circle { diameter: 6.0 },
            at: ogeom::math::Point::new(-10.0, 20.0, 30.0),
            frame: None,
            items: items.clone(),
        });
        targets.push(ogeom::doc::DatumTarget {
            datum: "A".into(),
            index: 3,
            kind: ogeom::doc::DatumTargetKind::Rectangle {
                length: 12.0,
                width: 4.0,
            },
            at: ogeom::math::Point::new(0.0, -20.0, 30.0),
            frame: None,
            items,
        });
    }

    let written = ogeom::io::write_step(&first.document, T).unwrap();
    let second = ogeom::io::read_step(&written, T).unwrap();
    let (before, after) = (first.document.pmi(), second.document.pmi());

    // The targets come back, letter, number, kind and place.
    let flatten = |pmi: &Pmi| -> Vec<String> {
        let mut out: Vec<String> = pmi
            .targets
            .iter()
            .map(|t| {
                format!(
                    "{}|{:?}|{:.6},{:.6},{:.6}",
                    t.identifier(),
                    t.kind,
                    t.at.x,
                    t.at.y,
                    t.at.z
                )
            })
            .collect();
        out.sort();
        out
    };
    assert_eq!(flatten(before), flatten(after));
    assert_eq!(after.targets.len(), 3);
    assert_eq!(after.targets_of("A").count(), 3);

    // And the presentation: the same callouts, the same planes, the same
    // drawn points, and the same links to the same annotations.
    let drawn = |pmi: &Pmi| -> Vec<String> {
        let mut out: Vec<String> = pmi
            .callouts
            .iter()
            .map(|c| {
                let points: usize = c.polylines.iter().map(Vec::len).sum();
                let plane = c
                    .plane
                    .map(|f| {
                        format!(
                            "{:.6},{:.6},{:.6}",
                            f.origin().x,
                            f.origin().y,
                            f.origin().z
                        )
                    })
                    .unwrap_or_default();
                format!(
                    "{}|{}|{}|{}|{:?}",
                    c.name,
                    c.polylines.len(),
                    points,
                    plane,
                    c.annotates
                )
            })
            .collect();
        out.sort();
        out
    };
    assert_eq!(drawn(before), drawn(after));
    assert_eq!(after.callouts.len(), 23);
}

/// The midpoint of a straight edge's two vertices, in world coordinates:
/// enough to tell a block's twelve edges apart across a round trip.
fn edge_midpoint(model: &ogeom::topo::Model, edge: &ogeom::topo::Shape) -> ogeom::math::Point {
    let ends: Vec<ogeom::math::Point> = model
        .children_of(edge)
        .unwrap()
        .iter()
        .map(|v| {
            let point = model.node(v).unwrap().data().as_vertex().unwrap().point;
            v.transform(model.datums()).unwrap().apply(point)
        })
        .collect();
    let (a, b) = (ends[0], ends[ends.len() - 1]);
    ogeom::math::Point::new((a.x + b.x) / 2.0, (a.y + b.y) / 2.0, (a.z + b.z) / 2.0)
}

/// An edge's colour and a dimension between two edges survive the writer:
/// the colour as a curve style on the edge's EDGE_CURVE, each end of the
/// dimension as a shape aspect naming one, and both land back on the
/// edges at the same place.
#[test]
fn edge_colours_and_edge_pmi_round_trip() {
    use ogeom::math::{Frame, Point};
    use ogeom::topo::{ShapeType, explore_unique};

    let mut model = ogeom::topo::Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 10.0, 5.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let at = |p: Point| {
        edges
            .iter()
            .find(|e| edge_midpoint(&model, e).distance(p) < 1e-9)
            .unwrap()
            .clone()
    };
    // Two parallel edges along x on the bottom face, 10 apart.
    let near = at(Point::new(10.0, 0.0, 0.0));
    let far = at(Point::new(10.0, 10.0, 0.0));
    let mut document = ogeom::doc::Document::over(model);
    document.add_part("block", block);
    let red = ogeom::doc::Colour::rgb(1.0, 0.0, 0.0);
    document.set_colour(&near, red);
    document.pmi_mut().dimensions.push(ogeom::doc::Dimension {
        name: "linear distance".into(),
        values: vec![10.0],
        kind: MeasureKind::Length,
        plus: Some(0.1),
        minus: Some(-0.1),
        features: vec![vec![near.node()], vec![far.node()]],
        location: true,
    });

    let text = ogeom::io::write_step(&document, T).unwrap();
    assert_eq!(text.matches("CURVE_STYLE(").count(), 1);
    let import = ogeom::io::read_step(&text, T).unwrap();
    let back = &import.document;
    let model = back.model();
    let edges = explore_unique(model, &import.solids[0], ShapeType::Edge).unwrap();
    let node_at = |p: Point| {
        edges
            .iter()
            .find(|e| edge_midpoint(model, e).distance(p) < 1e-9)
            .unwrap()
            .node()
    };

    let coloured: Vec<_> = back.colours().collect();
    assert_eq!(coloured, vec![(node_at(Point::new(10.0, 0.0, 0.0)), red)]);
    let dimensions = &back.pmi().dimensions;
    assert_eq!(dimensions.len(), 1);
    let dimension = &dimensions[0];
    assert_eq!(dimension.name, "linear distance");
    assert!(dimension.location);
    assert_eq!(dimension.values, vec![10.0]);
    assert_eq!(
        dimension.features,
        vec![
            vec![node_at(Point::new(10.0, 0.0, 0.0))],
            vec![node_at(Point::new(10.0, 10.0, 0.0))],
        ]
    );
}

/// The mean of a face's vertices, in world coordinates: enough to tell a
/// block's six faces apart across a round trip.
fn face_centre(model: &ogeom::topo::Model, face: &ogeom::topo::Shape) -> ogeom::math::Point {
    use ogeom::topo::{ShapeType, explore_unique};
    let vertices = explore_unique(model, face, ShapeType::Vertex).unwrap();
    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
    for v in &vertices {
        let point = model.node(v).unwrap().data().as_vertex().unwrap().point;
        let p = v.transform(model.datums()).unwrap().apply(point);
        (x, y, z) = (x + p.x, y + p.y, z + p.z);
    }
    #[allow(
        clippy::cast_precision_loss,
        reason = "a block's face has four vertices"
    )]
    let n = vertices.len() as f64;
    ogeom::math::Point::new(x / n, y / n, z / n)
}

/// A block with datums A, B and C on three faces, as a document, and the
/// block.
fn datum_block() -> (ogeom::doc::Document, ogeom::topo::Shape) {
    use ogeom::math::Frame;
    let mut model = ogeom::topo::Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 10.0, 5.0), T)
        .unwrap()
        .shape;
    let (a, b, c) = (
        face_at(&model, &block, (10.0, 5.0, 0.0)),
        face_at(&model, &block, (0.0, 5.0, 2.5)),
        face_at(&model, &block, (10.0, 0.0, 2.5)),
    );
    let mut document = ogeom::doc::Document::over(model);
    document.add_part("block", block.clone());
    for (label, face) in [("A", a), ("B", b), ("C", c)] {
        document.pmi_mut().datums.push(ogeom::doc::Datum {
            label: label.into(),
            items: vec![face],
        });
    }
    (document, block)
}

/// The block's face centred at `(x, y, z)`.
fn face_at(
    model: &ogeom::topo::Model,
    block: &ogeom::topo::Shape,
    (x, y, z): (f64, f64, f64),
) -> ogeom::topo::TShapeId {
    use ogeom::topo::{ShapeType, explore_unique};
    explore_unique(model, block, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| face_centre(model, f).distance(ogeom::math::Point::new(x, y, z)) < 1e-9)
        .unwrap()
        .node()
}

/// Where the faces of `block` among `items` sit, sorted.
fn centres_of(
    model: &ogeom::topo::Model,
    block: &ogeom::topo::Shape,
    items: &[ogeom::topo::TShapeId],
) -> Vec<String> {
    use ogeom::topo::{ShapeType, explore_unique};
    let mut out: Vec<String> = explore_unique(model, block, ShapeType::Face)
        .unwrap()
        .iter()
        .filter(|f| items.contains(&f.node()))
        .map(|f| {
            let p = face_centre(model, f);
            format!("{:.6} {:.6} {:.6}", p.x, p.y, p.z)
        })
        .collect();
    out.sort();
    out
}

/// Datum reference frames, material modifiers, and the position, profile
/// and runout kinds go through the writer as AP242 datum systems and come
/// back as they went in: kind, value, the tolerance's modifiers, each
/// frame's references in precedence order with their own modifiers, and
/// the faces the tolerance controls.
#[test]
fn datum_systems_and_their_modifiers_round_trip() {
    use ogeom::doc::{
        DatumModifier, DatumModifierKind, DatumReference, GeometricTolerance, LEAST_MATERIAL,
        MAXIMUM_MATERIAL,
    };
    let (mut document, block) = datum_block();
    let face = |p| face_at(document.model(), &block, p);
    let (top, right, back) = (
        face((10.0, 5.0, 5.0)),
        face((20.0, 5.0, 2.5)),
        face((10.0, 10.0, 2.5)),
    );
    let tolerance =
        |kind: &str, magnitude: f64, modifiers: &[&str], datums, items| GeometricTolerance {
            kind: kind.into(),
            name: format!("{kind}.1"),
            magnitude,
            modifiers: modifiers.iter().map(|&m| m.into()).collect(),
            datums,
            datum_system: true,
            items,
        };
    let written = vec![
        tolerance(
            "position",
            0.1,
            &[MAXIMUM_MATERIAL],
            vec![
                "A".into(),
                DatumReference::new("B").with(DatumModifier::maximum_material()),
                DatumReference::new("C").with(DatumModifier::least_material()),
            ],
            vec![top],
        ),
        tolerance(
            "surface_profile",
            0.2,
            &[],
            vec!["A".into(), "B".into(), "C".into()],
            vec![top, right],
        ),
        tolerance(
            "line_profile",
            0.05,
            &[],
            vec!["A".into(), "B".into()],
            vec![back],
        ),
        tolerance(
            "circular_runout",
            0.02,
            &[],
            vec![DatumReference::new("A-B").with(DatumModifier::WithValue {
                kind: DatumModifierKind::Projected,
                value: 12.5,
            })],
            vec![right],
        ),
        tolerance(
            "total_runout",
            0.03,
            &[LEAST_MATERIAL],
            vec![
                DatumReference::new("A")
                    .with(DatumModifier::WithValue {
                        kind: DatumModifierKind::CircularOrCylindrical,
                        value: 8.0,
                    })
                    .with(DatumModifier::Simple("free_state".into())),
                "C".into(),
            ],
            vec![back],
        ),
    ];
    document.pmi_mut().tolerances.clone_from(&written);

    let text = ogeom::io::write_step(&document, T).unwrap();
    assert!(!text.contains("GEOMETRIC_TOLERANCE_WITH_MODIFIERS(())"));
    assert_eq!(text.matches("DATUM_SYSTEM(").count(), written.len());
    assert_eq!(text.matches("DATUM_REFERENCE_ELEMENT(").count(), 2);
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let read = &import.document.pmi().tolerances;
    assert_eq!(read.len(), written.len());
    for (before, after) in written.iter().zip(read) {
        assert_eq!(after.kind, before.kind);
        assert_eq!(after.name, before.name);
        assert!((after.magnitude - before.magnitude).abs() < 1e-12);
        assert_eq!(after.modifiers, before.modifiers);
        assert_eq!(after.datums, before.datums, "{}", before.kind);
        assert!(after.datum_system);
        assert_eq!(
            centres_of(model, &import.solids[0], &after.items),
            centres_of(document.model(), &block, &before.items),
            "{}",
            before.kind
        );
    }

    // And through the native format, which keeps the frame in a record of
    // its own.
    let native = ogeom::io::native::write_document(
        &import.document,
        ogeom::io::native::WriteOptions::default(),
    )
    .unwrap();
    let again = ogeom::io::native::read_document(&native).unwrap();
    for (before, after) in written.iter().zip(&again.pmi().tolerances) {
        assert_eq!(after.datums, before.datums);
        assert!(after.datum_system);
    }
}

/// A datum system holds at most three references, and each must name a
/// datum: a frame that cannot be stated is refused, not written short.
#[test]
fn a_datum_system_that_cannot_be_stated_is_refused() {
    let tolerance = |datums: Vec<ogeom::doc::DatumReference>| ogeom::doc::GeometricTolerance {
        kind: "position".into(),
        magnitude: 0.1,
        datums,
        datum_system: true,
        ..Default::default()
    };
    let (mut document, _) = datum_block();
    for frame in [
        vec!["A".into(), "B".into(), "C".into(), "A-B".into()],
        vec!["A".into(), "D".into()],
        vec!["A".into(), "B-D".into()],
    ] {
        document.pmi_mut().tolerances = vec![tolerance(frame)];
        assert!(ogeom::io::write_step(&document, T).is_err());
    }
}
