//! AP242 semantic PMI: the values, not the leader lines — read from NIST's
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
                t.datums.clone(),
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

    // Six geometric tolerances, magnitudes and datum references as NIST
    // published them.
    let mut kinds: Vec<(&str, f64, &[String])> = pmi
        .tolerances
        .iter()
        .map(|t| (t.kind.as_str(), t.magnitude, t.datums.as_slice()))
        .collect();
    kinds.sort_by(|a, b| a.0.cmp(b.0).then(a.1.total_cmp(&b.1)));
    let a = ["A".to_string()];
    let none: [String; 0] = [];
    assert_eq!(
        kinds,
        [
            ("flatness", 0.2, &none[..]),
            ("perpendicularity", 1.5, &a[..]),
            ("position", 0.75, &a[..]),
            ("position", 0.75, &a[..]),
            ("surface_profile", 0.5, &a[..]),
            ("surface_profile", 1.25, &a[..]),
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

/// §E2: the *presentation* half. NIST's part draws its annotations, and what
/// it draws is read — the callouts, the plane each is drawn in, and the
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
    // annotation is anchored in — this file offsets the drawing from that
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
