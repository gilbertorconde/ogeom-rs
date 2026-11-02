//! Reading DXF drawings as typed curves: arcs, circles, bulged and closed
//! polylines, units and dashed layers; and writing them back as R2000.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI, TAU};
use ogeom_core::{FastMap, FastSet};
use ogeom_io::dxf::{
    DxfCurve, DxfItem, DxfLayer, DxfOut, DxfShape, DxfText, DxfUnit, read_dxf, read_dxf_entities,
    write_dxf_entities,
};
use ogeom_math::{Point2, Vector2};

const REPRO: &str = "0\nSECTION\n2\nENTITIES\n0\nCIRCLE\n8\n0\n10\n0.0\n20\n0.0\n40\n5.0\n0\nARC\n8\n0\n10\n20.0\n20\n0.0\n40\n4.0\n50\n0.0\n51\n90.0\n0\nLWPOLYLINE\n8\n0\n90\n2\n70\n1\n10\n40.0\n20\n0.0\n42\n1.0\n10\n50.0\n20\n0.0\n42\n1.0\n0\nPOLYLINE\n8\n0\n66\n1\n10\n0.0\n20\n0.0\n30\n0.0\n70\n0\n0\nVERTEX\n8\n0\n10\n5.0\n20\n5.0\n0\nVERTEX\n8\n0\n10\n15.0\n20\n5.0\n0\nSEQEND\n0\nENDSEC\n0\nEOF\n";

#[test]
fn circles_arcs_and_bulged_closed_polylines_read_as_written() {
    let d = read_dxf_entities(REPRO).unwrap();
    assert_eq!(d.entities.len(), 4);
    assert!(matches!(d.entities[0].curve, DxfCurve::Circle { radius, .. } if radius == 5.0));
    assert!(matches!(d.entities[1].curve,
        DxfCurve::Arc { radius, start_angle, end_angle, .. }
            if radius == 4.0 && start_angle == 0.0 && (end_angle - FRAC_PI_2).abs() < 1e-12));
    let DxfCurve::Polyline { vertices, closed } = &d.entities[2].curve else {
        panic!("{:?}", d.entities[2]);
    };
    assert!(*closed);
    assert_eq!(
        vertices.iter().map(|v| v.1).collect::<Vec<_>>(),
        vec![1.0, 1.0]
    );
    let DxfCurve::Polyline { vertices, closed } = &d.entities[3].curve else {
        panic!("{:?}", d.entities[3]);
    };
    assert!(!*closed);
    assert_eq!(
        vertices.iter().map(|v| (v.0.x, v.0.y)).collect::<Vec<_>>(),
        vec![(5.0, 5.0), (15.0, 5.0)]
    );
    // The polyline reader does not take the header's point for a vertex.
    assert_eq!(read_dxf(REPRO).unwrap().visible.last().unwrap().len(), 2);
    assert_eq!(d.insunits, None);
}

#[test]
fn the_header_names_the_units() {
    let text = format!("0\nSECTION\n2\nHEADER\n9\n$INSUNITS\n70\n1\n0\nENDSEC\n{REPRO}");
    let d = read_dxf_entities(&text).unwrap();
    assert_eq!(d.insunits, Some(1));
    assert_eq!(d.unit_mm, Some(25.4));
    assert_eq!(d.entities.len(), 4);
}

#[test]
fn a_dashed_layer_draws_hidden_lines() {
    let text = "0\nSECTION\n2\nTABLES\n0\nTABLE\n2\nLAYER\n0\nLAYER\n2\nBACK\n70\n0\n6\nDASHED\n0\nENDTAB\n0\nENDSEC\n0\nSECTION\n2\nENTITIES\n0\nLINE\n8\nBACK\n10\n0.0\n20\n0.0\n11\n1.0\n21\n0.0\n0\nLINE\n8\nFRONT\n10\n0.0\n20\n1.0\n11\n1.0\n21\n1.0\n0\nENDSEC\n0\nEOF\n";
    let d = read_dxf_entities(text).unwrap();
    assert!(d.entities[0].hidden);
    assert!(!d.entities[1].hidden);
    assert_eq!(d.entities[1].layer, "FRONT");
}

#[test]
fn an_arc_seen_from_below_is_mirrored_into_the_drawing() {
    let text = "0\nSECTION\n2\nENTITIES\n0\nARC\n8\n0\n10\n3.0\n20\n0.0\n40\n1.0\n50\n0.0\n51\n90.0\n210\n0.0\n220\n0.0\n230\n-1.0\n0\nENDSEC\n0\nEOF\n";
    let d = read_dxf_entities(text).unwrap();
    let DxfCurve::Arc {
        centre,
        start_angle,
        end_angle,
        ..
    } = d.entities[0].curve
    else {
        panic!("{:?}", d.entities[0]);
    };
    assert_eq!((centre.x, centre.y), (-3.0, 0.0));
    assert!((start_angle - FRAC_PI_2).abs() < 1e-12);
    assert!((end_angle - PI).abs() < 1e-12);
}

#[test]
fn text_seen_from_below_is_mirrored_and_its_escapes_read() {
    let text = "0\nSECTION\n2\nENTITIES\n0\nTEXT\n8\nNOTES\n10\n2.0\n20\n1.0\n40\n3.5\n1\n\\U+00D812 \\x\n50\n30.0\n210\n0.0\n220\n0.0\n230\n-1.0\n0\nENDSEC\n0\nEOF\n";
    let d = read_dxf_entities(text).unwrap();
    assert!(d.entities.is_empty());
    let label = &d.texts[0];
    assert_eq!(label.layer, "NOTES");
    assert_eq!((label.text.at.x, label.text.at.y), (-2.0, 1.0));
    assert_eq!(label.text.height, 3.5);
    assert!((label.text.angle - 5.0 * PI / 6.0).abs() < 1e-12);
    assert_eq!(label.text.text, "\u{d8}12 \\x");
}

fn layer(name: &str, dashed: bool) -> DxfLayer {
    DxfLayer {
        name: name.into(),
        dashed,
    }
}

fn curve(layer: usize, curve: DxfCurve) -> DxfItem {
    DxfItem {
        layer,
        shape: DxfShape::Curve(curve),
    }
}

/// One of every entity the writer takes, on three layers.
fn every_entity() -> (Vec<DxfLayer>, Vec<DxfItem>) {
    let layers = vec![
        layer("VISIBLE", false),
        layer("HIDDEN", true),
        layer("Notes and dims", false),
    ];
    let p = Point2::new;
    let items = vec![
        curve(
            0,
            DxfCurve::Line {
                start: p(0.1, -2.5),
                end: p(1.0 / 3.0, 7.0e5),
            },
        ),
        curve(
            0,
            DxfCurve::Circle {
                centre: p(12.0, 4.0),
                radius: 6.0,
            },
        ),
        curve(
            1,
            DxfCurve::Arc {
                centre: p(-3.0, 2.0),
                radius: 2.5,
                start_angle: 0.3,
                end_angle: 2.0 * PI / 3.0 + 1.0,
            },
        ),
        curve(
            1,
            DxfCurve::Ellipse {
                centre: p(5.0, 5.0),
                major: Vector2::new(3.0, 4.0),
                ratio: 0.4,
                start_param: 0.25,
                end_param: 4.5,
            },
        ),
        curve(
            0,
            DxfCurve::Polyline {
                vertices: vec![
                    (p(0.0, 0.0), 0.0),
                    (p(10.0, 0.0), (FRAC_PI_2 / 4.0).tan()),
                    (p(10.0, 10.0), -0.5),
                    (p(0.0, 10.0), 0.0),
                ],
                closed: true,
            },
        ),
        // A rational quadratic: a quarter circle and a line, clamped.
        curve(
            0,
            DxfCurve::Spline {
                degree: 2,
                knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
                control_points: vec![p(1.0, 0.0), p(1.0, 1.0), p(0.0, 1.0), p(-2.0, 1.0)],
                weights: Some(vec![1.0, FRAC_PI_4.cos(), 1.0, 1.0]),
                closed: false,
            },
        ),
        curve(
            1,
            DxfCurve::Spline {
                degree: 3,
                knots: Vec::new(),
                control_points: vec![p(0.0, 0.0), p(1.0, 2.0), p(3.0, -1.0), p(4.0, 0.5)],
                weights: None,
                closed: false,
            },
        ),
        DxfItem {
            layer: 2,
            shape: DxfShape::Text(DxfText {
                at: p(20.0, -4.0),
                height: 2.5,
                angle: 0.1,
                text: " \u{d8}12 \u{b1}0.1 %%d 50% \\U+0041 \u{2300} ".into(),
            }),
        },
    ];
    (layers, items)
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-9
}

fn points_close(a: &[Point2], b: &[Point2]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(p, q)| p.distance(*q) <= 1e-9)
}

fn reals_close(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| close(*x, *y))
}

/// The same curve, to 1e-9 in every number.
fn same_curve(a: &DxfCurve, b: &DxfCurve) -> bool {
    match (a, b) {
        (DxfCurve::Line { start, end }, DxfCurve::Line { start: s, end: e }) => {
            points_close(&[*start, *end], &[*s, *e])
        }
        (
            DxfCurve::Circle { centre, radius },
            DxfCurve::Circle {
                centre: c,
                radius: r,
            },
        ) => centre.distance(*c) <= 1e-9 && close(*radius, *r),
        (
            DxfCurve::Arc {
                centre,
                radius,
                start_angle,
                end_angle,
            },
            DxfCurve::Arc {
                centre: c,
                radius: r,
                start_angle: s,
                end_angle: e,
            },
        ) => {
            centre.distance(*c) <= 1e-9
                && reals_close(&[*radius, *start_angle, *end_angle], &[*r, *s, *e])
        }
        (
            DxfCurve::Ellipse {
                centre,
                major,
                ratio,
                start_param,
                end_param,
            },
            DxfCurve::Ellipse {
                centre: c,
                major: m,
                ratio: r,
                start_param: s,
                end_param: e,
            },
        ) => {
            centre.distance(*c) <= 1e-9
                && reals_close(
                    &[major.x, major.y, *ratio, *start_param, *end_param],
                    &[m.x, m.y, *r, *s, *e],
                )
        }
        (
            DxfCurve::Polyline { vertices, closed },
            DxfCurve::Polyline {
                vertices: v,
                closed: c,
            },
        ) => {
            closed == c
                && vertices.len() == v.len()
                && vertices
                    .iter()
                    .zip(v)
                    .all(|(a, b)| a.0.distance(b.0) <= 1e-9 && close(a.1, b.1))
        }
        (
            DxfCurve::Spline {
                degree,
                knots,
                control_points,
                weights,
                closed,
            },
            DxfCurve::Spline {
                degree: d,
                knots: k,
                control_points: p,
                weights: w,
                closed: c,
            },
        ) => {
            degree == d
                && closed == c
                && reals_close(knots, k)
                && points_close(control_points, p)
                && match (weights, w) {
                    (None, None) => true,
                    (Some(a), Some(b)) => reals_close(a, b),
                    _ => false,
                }
        }
        _ => false,
    }
}

#[test]
fn every_entity_written_reads_back_as_itself() {
    let (layers, items) = every_entity();
    for (unit, insunits, mm) in [
        (DxfUnit::Millimetre, 4, Some(1.0)),
        (DxfUnit::Inch, 1, Some(25.4)),
        (DxfUnit::Unitless, 0, None),
    ] {
        let text = write_dxf_entities(&DxfOut {
            unit,
            layers: &layers,
            items: &items,
        })
        .unwrap();
        assert!(text.is_ascii());
        let d = read_dxf_entities(&text).unwrap();
        assert_eq!(d.insunits, Some(insunits));
        assert_eq!(d.unit_mm, mm);

        let curves: Vec<&DxfItem> = items
            .iter()
            .filter(|i| matches!(i.shape, DxfShape::Curve(_)))
            .collect();
        assert_eq!(d.entities.len(), curves.len());
        for (read, item) in d.entities.iter().zip(&curves) {
            let DxfShape::Curve(written) = &item.shape else {
                unreachable!()
            };
            assert!(
                same_curve(&read.curve, written),
                "{read:?}\nagainst {written:?}"
            );
            assert_eq!(read.layer, layers[item.layer].name);
            assert_eq!(read.hidden, layers[item.layer].dashed);
        }

        let DxfShape::Text(written) = &items[7].shape else {
            unreachable!()
        };
        let read = &d.texts[0];
        assert_eq!(read.layer, "Notes and dims");
        assert_eq!(read.text.text, written.text);
        assert!(read.text.at.distance(written.at) <= 1e-9);
        assert!(close(read.text.height, written.height));
        assert!(close(read.text.angle, written.angle));
    }
}

#[test]
fn a_bulge_reads_back_as_the_arc_it_was_made_from() {
    // A slot's end: a line, then a counter-clockwise half circle of radius
    // 2 about (10, 2), then a clockwise quarter of radius 3 about (10, 7).
    let p = Point2::new;
    let half = (PI / 4.0).tan();
    let quarter = -(FRAC_PI_2 / 4.0).tan();
    let written = vec![
        (p(0.0, 0.0), 0.0),
        (p(10.0, 0.0), half),
        (p(10.0, 4.0), quarter),
        (p(7.0, 7.0), 0.0),
    ];
    let text = write_dxf_entities(&DxfOut {
        unit: DxfUnit::Millimetre,
        layers: &[layer("0", false)],
        items: &[curve(
            0,
            DxfCurve::Polyline {
                vertices: written,
                closed: false,
            },
        )],
    })
    .unwrap();
    let d = read_dxf_entities(&text).unwrap();
    let DxfCurve::Polyline { vertices, closed } = &d.entities[0].curve else {
        panic!("{:?}", d.entities[0]);
    };
    assert!(!*closed);
    assert_eq!(vertices.len(), 4);

    // The arc a bulge b draws from a to c: included angle 4 atan b,
    // counter-clockwise when b > 0, centre on the chord's bisector.
    let arc = |a: Point2, c: Point2, b: f64| {
        let theta = 4.0 * b.atan();
        let chord = c - a;
        let radius = chord.magnitude() / (2.0 * (theta / 2.0).sin());
        let mid = a + chord * 0.5;
        let normal = Vector2::new(-chord.y, chord.x) * (1.0 / chord.magnitude());
        let offset = radius * (theta / 2.0).cos();
        (mid + normal * offset, radius.abs(), theta)
    };
    assert_eq!(vertices[0].1, 0.0);
    let (centre, radius, sweep) = arc(vertices[1].0, vertices[2].0, vertices[1].1);
    assert!(centre.distance(p(10.0, 2.0)) <= 1e-9, "{centre:?}");
    assert!(close(radius, 2.0) && close(sweep, PI));
    let (centre, radius, sweep) = arc(vertices[2].0, vertices[3].0, vertices[2].1);
    assert!(centre.distance(p(10.0, 7.0)) <= 1e-9, "{centre:?}");
    assert!(close(radius, 3.0) && close(sweep, -FRAC_PI_2));
    assert_eq!(vertices[3].1, 0.0);
}

/// A record: its kind and its pairs.
type Record = (String, Vec<(i32, String)>);

/// The file's records by section, each section's pairs before its first
/// record (the HEADER's variables) first.
fn sections(text: &str) -> Vec<(String, Vec<Record>)> {
    let lines: Vec<&str> = text.lines().collect();
    let pairs: Vec<(i32, String)> = lines
        .chunks(2)
        .map(|c| (c[0].trim().parse().unwrap(), c[1].to_string()))
        .collect();
    let mut out: Vec<(String, Vec<Record>)> = Vec::new();
    let mut k = 0;
    while k < pairs.len() {
        assert_eq!(pairs[k].0, 0, "every record starts with group 0");
        match pairs[k].1.as_str() {
            "SECTION" => {
                assert_eq!(pairs[k + 1].0, 2);
                out.push((pairs[k + 1].1.clone(), vec![(String::new(), Vec::new())]));
                k += 1;
            }
            "EOF" => {
                assert_eq!(k + 1, pairs.len(), "EOF ends the file");
                break;
            }
            kind => out
                .last_mut()
                .unwrap()
                .1
                .push((kind.to_string(), Vec::new())),
        }
        k += 1;
        while k < pairs.len() && pairs[k].0 != 0 {
            let record = out.last_mut().unwrap().1.last_mut().unwrap();
            record.1.push(pairs[k].clone());
            k += 1;
        }
    }
    for (_, records) in &mut out {
        assert_eq!(records.last().unwrap().0, "ENDSEC");
        records.pop();
    }
    out
}

fn group(pairs: &[(i32, String)], code: i32) -> Option<&str> {
    pairs.iter().find(|p| p.0 == code).map(|p| p.1.as_str())
}

#[test]
fn the_r2000_file_holds_every_section_table_and_handle_a_strict_reader_needs() {
    let (layers, items) = every_entity();
    let text = write_dxf_entities(&DxfOut {
        unit: DxfUnit::Millimetre,
        layers: &layers,
        items: &items,
    })
    .unwrap();
    let sections = sections(&text);
    let names: Vec<&str> = sections.iter().map(|s| s.0.as_str()).collect();
    assert_eq!(
        names,
        [
            "HEADER", "CLASSES", "TABLES", "BLOCKS", "ENTITIES", "OBJECTS"
        ]
    );

    let header = &sections[0].1[0].1;
    let variable = |name: &str| {
        let at = header.iter().position(|p| *p == (9, name.to_string()));
        at.map(|i| header[i + 1].1.clone())
    };
    assert_eq!(variable("$ACADVER").as_deref(), Some("AC1015"));
    assert_eq!(variable("$INSUNITS").as_deref(), Some("4"));
    let seed = u64::from_str_radix(&variable("$HANDSEED").unwrap(), 16).unwrap();

    // Every record past the HEADER has a unique handle below the seed, and
    // every owner it names is a record of the file.
    let mut handles: FastMap<u64, (String, String)> = FastMap::default();
    let mut owners: Vec<(String, u64)> = Vec::new();
    for (section, records) in &sections[1..] {
        assert!(records[0].1.is_empty(), "{section} opens with a record");
        for (kind, pairs) in &records[1..] {
            if section == "TABLES" && kind == "ENDTAB" {
                continue;
            }
            let h = group(pairs, 5)
                .or_else(|| group(pairs, 105))
                .unwrap_or_else(|| panic!("a {kind} in {section} has no handle"));
            let h = u64::from_str_radix(h, 16).unwrap();
            assert!(h > 0 && h < seed, "handle {h:X} against seed {seed:X}");
            let name = group(pairs, 2).unwrap_or("").to_string();
            assert!(
                handles.insert(h, (kind.clone(), name)).is_none(),
                "handle {h:X} twice"
            );
            let owner = u64::from_str_radix(group(pairs, 330).unwrap(), 16).unwrap();
            owners.push((kind.clone(), owner));
        }
    }
    for (kind, owner) in &owners {
        assert!(
            *owner == 0 || handles.contains_key(owner),
            "a {kind} is owned by {owner:X}, no record"
        );
    }

    // The nine tables in order, each counting its entries, each entry owned
    // by its table.
    let tables = &sections[2].1;
    let mut found: Vec<(String, Vec<String>)> = Vec::new();
    for (kind, pairs) in &tables[1..] {
        match kind.as_str() {
            "TABLE" => {
                found.push((group(pairs, 2).unwrap().to_string(), Vec::new()));
                assert_eq!(group(pairs, 330), Some("0"));
                assert_eq!(group(pairs, 100), Some("AcDbSymbolTable"));
            }
            "ENDTAB" => {}
            entry => {
                let table = found.last_mut().unwrap();
                assert_eq!(*entry, table.0, "a {entry} in table {}", table.0);
                assert!(pairs.contains(&(100, "AcDbSymbolTableRecord".into())));
                table.1.push(group(pairs, 2).unwrap().to_string());
            }
        }
    }
    let table_names: Vec<&str> = found.iter().map(|t| t.0.as_str()).collect();
    assert_eq!(
        table_names,
        [
            "VPORT",
            "LTYPE",
            "LAYER",
            "STYLE",
            "VIEW",
            "UCS",
            "APPID",
            "DIMSTYLE",
            "BLOCK_RECORD"
        ]
    );
    let entries = |name: &str| &found.iter().find(|t| t.0 == name).unwrap().1;
    for (kind, pairs) in tables.iter().filter(|r| r.0 == "TABLE") {
        let name = group(pairs, 2).unwrap();
        let count: usize = group(pairs, 70).unwrap().parse().unwrap();
        assert_eq!(count, entries(name).len(), "{kind} {name}");
    }
    for linetype in ["ByBlock", "ByLayer", "CONTINUOUS", "DASHED"] {
        assert!(entries("LTYPE").iter().any(|l| l == linetype), "{linetype}");
    }
    assert!(entries("STYLE").iter().any(|s| s == "Standard"));
    assert!(entries("APPID").iter().any(|s| s == "ACAD"));
    assert!(entries("DIMSTYLE").iter().any(|s| s == "Standard"));
    assert_eq!(*entries("BLOCK_RECORD"), ["*Model_Space", "*Paper_Space"]);
    let declared: FastSet<&str> = entries("LAYER").iter().map(String::as_str).collect();
    assert!(declared.contains("0"));
    for (_, pairs) in tables.iter().filter(|r| r.0 == "LAYER") {
        assert!(
            entries("LTYPE")
                .iter()
                .any(|l| l == group(pairs, 6).unwrap())
        );
    }

    // Each space has its block, owned by its record; every entity is in
    // model space, on a declared layer.
    let record_of = |name: &str| {
        *handles
            .iter()
            .find(|(_, v)| v.0 == "BLOCK_RECORD" && v.1 == name)
            .unwrap()
            .0
    };
    let model = record_of("*Model_Space");
    let blocks = &sections[3].1;
    let kinds: Vec<&str> = blocks[1..].iter().map(|r| r.0.as_str()).collect();
    assert_eq!(kinds, ["BLOCK", "ENDBLK", "BLOCK", "ENDBLK"]);
    for (i, space) in ["*Model_Space", "*Paper_Space"].into_iter().enumerate() {
        let owner = format!("{:X}", record_of(space));
        assert_eq!(group(&blocks[1 + 2 * i].1, 2), Some(space));
        assert_eq!(group(&blocks[1 + 2 * i].1, 330), Some(owner.as_str()));
        assert_eq!(group(&blocks[2 + 2 * i].1, 330), Some(owner.as_str()));
    }
    let entities = &sections[4].1[1..];
    assert_eq!(entities.len(), items.len());
    for (kind, pairs) in entities {
        assert_eq!(
            group(pairs, 330),
            Some(format!("{model:X}").as_str()),
            "{kind}"
        );
        assert_eq!(group(pairs, 100), Some("AcDbEntity"), "{kind}");
        assert!(declared.contains(group(pairs, 8).unwrap()), "{kind}");
    }

    // The root dictionary, owned by none, names the group dictionary,
    // owned by the root.
    let objects = &sections[5].1[1..];
    let root = &objects[0];
    assert_eq!(
        (root.0.as_str(), group(&root.1, 330)),
        ("DICTIONARY", Some("0"))
    );
    let at = root
        .1
        .iter()
        .position(|p| *p == (3, "ACAD_GROUP".into()))
        .unwrap();
    let group_dict = u64::from_str_radix(&root.1[at + 1].1, 16).unwrap();
    assert_eq!(root.1[at + 1].0, 350);
    assert_eq!(handles[&group_dict].0, "DICTIONARY");
    let owned = objects
        .iter()
        .find(|o| group(&o.1, 5) == Some(format!("{group_dict:X}").as_str()))
        .unwrap();
    assert_eq!(group(&owned.1, 330), group(&root.1, 5));
}

#[test]
fn a_drawing_the_format_cannot_hold_is_refused_by_name() {
    let (layers, items) = every_entity();
    let write = |layers: &[DxfLayer], items: &[DxfItem]| {
        write_dxf_entities(&DxfOut {
            unit: DxfUnit::Millimetre,
            layers,
            items,
        })
    };
    assert!(write(&layers, &items).is_ok());
    for bad in [
        "",
        " lead",
        "a/b",
        "semi;colon",
        "tab\there",
        "\u{e9}t\u{e9}",
    ] {
        assert!(write(&[layer(bad, false)], &[]).is_err(), "{bad:?}");
    }
    assert!(write(&[layer("A", false), layer("a", true)], &[]).is_err());
    let mut on_no_layer = items[0].clone();
    on_no_layer.layer = 3;
    assert!(write(&layers, &[on_no_layer]).is_err());

    let p = Point2::new;
    let refused = [
        DxfShape::Curve(DxfCurve::Circle {
            centre: p(0.0, 0.0),
            radius: 0.0,
        }),
        DxfShape::Curve(DxfCurve::Line {
            start: p(f64::NAN, 0.0),
            end: p(1.0, 0.0),
        }),
        DxfShape::Curve(DxfCurve::Ellipse {
            centre: p(0.0, 0.0),
            major: Vector2::new(1.0, 0.0),
            ratio: 1.5,
            start_param: 0.0,
            end_param: TAU,
        }),
        DxfShape::Curve(DxfCurve::Polyline {
            vertices: vec![(p(0.0, 0.0), 0.0)],
            closed: false,
        }),
        DxfShape::Curve(DxfCurve::Spline {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0],
            control_points: vec![p(0.0, 0.0), p(1.0, 1.0), p(2.0, 0.0)],
            weights: None,
            closed: false,
        }),
        DxfShape::Curve(DxfCurve::Spline {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![p(0.0, 0.0), p(1.0, 1.0)],
            weights: Some(vec![1.0, 0.0]),
            closed: false,
        }),
        DxfShape::Text(DxfText {
            at: p(0.0, 0.0),
            height: 1.0,
            angle: 0.0,
            text: "two\nlines".into(),
        }),
    ];
    for shape in refused {
        let item = DxfItem { layer: 0, shape };
        assert!(
            write(&layers, core::slice::from_ref(&item)).is_err(),
            "{item:?}"
        );
    }
}
