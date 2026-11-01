//! STEP edges dressed as `SURFACE_CURVE`: what every exporter derived from
//! the reference kernel writes.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/pcurves.rs"]
mod pcurves;

use ogeom::topo::Model;
use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point};
use ogeom_mesh::Deflection;

const T: Tolerances = Tolerances::millimetres();

/// Re-dress every `EDGE_CURVE`'s geometry in a `SURFACE_CURVE` wrapper, the
/// way exporters derived from the reference kernel write their files.
fn dressed(text: &str) -> String {
    let mut highest: u64 = 0;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('#')
            && let Some(eq) = rest.find('=')
            && let Ok(id) = rest[..eq].trim().parse::<u64>()
        {
            highest = highest.max(id);
        }
    }
    let mut out = String::new();
    let mut extra = String::new();
    for line in text.lines() {
        if let Some(at) = line.find("EDGE_CURVE('',") {
            // #id=EDGE_CURVE('',#a,#b,#c,.T.);
            let args = &line[at + "EDGE_CURVE('',".len()..];
            let mut parts = args.split(',');
            let a = parts.next().unwrap();
            let b = parts.next().unwrap();
            let c = parts.next().unwrap();
            let tail: Vec<&str> = parts.collect();
            highest += 1;
            let wrapper = highest;
            extra.push_str(&format!(
                "#{wrapper}=SURFACE_CURVE('',{c},(),.CURVE_3D.);\n"
            ));
            out.push_str(&line[..at]);
            out.push_str(&format!(
                "EDGE_CURVE('',{a},{b},#{wrapper},{}",
                tail.join(",")
            ));
            out.push('\n');
        } else if line.trim() == "ENDSEC;" && !extra.is_empty() {
            out.push_str(&extra);
            extra.clear();
            out.push_str(line);
            out.push('\n');
        } else {
            out.push_str(line);
            out.push('\n');
        }
    }
    out
}

#[test]
fn an_edge_dressed_as_a_surface_curve_still_reads() {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 10.0, 5.0), T)
        .unwrap()
        .shape;
    let mut document = ogeom::doc::Document::over(model);
    document.add_part("block", block);
    let text = ogeom::io::write_step(&document, T).unwrap();
    let wrapped = dressed(&text);
    assert!(
        wrapped.contains(".CURVE_3D."),
        "the fixture really re-dressed the edges"
    );

    let import = ogeom::io::read_step(&wrapped, T).unwrap();
    let back = &import.document;
    let root = back.roots()[0];
    let occurrence = &back.occurrences_of(root).unwrap()[0];
    let volume =
        ogeom::algo::volume_properties(back.model(), &occurrence.shape, Deflection::default(), T)
            .unwrap()
            .mass;
    assert!(
        (volume - 1000.0).abs() / 1000.0 < 0.01,
        "volume {volume} against 1000"
    );
    assert!(
        import
            .report
            .warnings
            .iter()
            .all(|w| !w.contains("skipped")),
        "no edge was skipped: {:?}",
        import.report.warnings
    );
}

/// Primitives with seams and poles, and a loft whose skew walls are spline
/// patches, written and read back: every pcurve the file carries comes
/// back as the source held it, a seam's two sides included, and the solid
/// is valid with its volume.
#[test]
fn pcurves_written_with_their_edges_come_back_as_written() {
    let mut model = Model::new();
    let at = |x: f64| Frame::new(Point::new(x, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
    let mut parts = vec![
        (
            "cylinder",
            ogeom::algo::make_cylinder(&mut model, at(0.0), 3.0, 5.0, T)
                .unwrap()
                .shape,
        ),
        (
            "cone",
            ogeom::algo::make_cone(&mut model, at(10.0), 3.0, 1.0, 4.0, T)
                .unwrap()
                .shape,
        ),
        (
            "sphere",
            ogeom::algo::make_sphere(&mut model, at(20.0), 2.5, T)
                .unwrap()
                .shape,
        ),
        (
            "torus",
            ogeom::algo::make_torus(&mut model, at(30.0), 4.0, 1.0, T)
                .unwrap()
                .shape,
        ),
    ];
    // A square lofted to the same square an eighth of a turn round: its
    // four walls are bilinear spline patches.
    let square = |model: &mut Model, turn: f64, z: f64| {
        let corners: Vec<Point> = (0..4)
            .map(|i| {
                let angle = turn + core::f64::consts::FRAC_PI_2 * f64::from(i);
                Point::new(50.0 + 3.0 * angle.cos(), 3.0 * angle.sin(), z)
            })
            .collect();
        ogeom::algo::make_polygon(model, &corners, true, T)
            .unwrap()
            .shape
    };
    let bottom = square(&mut model, 0.0, 0.0);
    let top = square(&mut model, core::f64::consts::FRAC_PI_4, 5.0);
    parts.push((
        "twisted loft",
        ogeom::offset::make_loft(&mut model, &bottom, &top, T)
            .unwrap()
            .shape,
    ));
    let fine = Deflection::with_chord(1e-3).unwrap();
    for (name, shape) in parts {
        let before = ogeom::algo::volume_properties(&model, &shape, fine, T)
            .unwrap()
            .mass;
        let mut document = ogeom::doc::Document::over(model.clone());
        document.add_part(name, shape.clone());
        let text = ogeom::io::write_step(&document, T).unwrap();
        assert!(text.contains("PCURVE("), "{name}: the edges carry pcurves");
        if name != "twisted loft" {
            assert!(
                text.contains("SEAM_CURVE("),
                "{name}: a seam is written as one"
            );
        }
        let import = ogeom::io::read_step(&text, T).unwrap();
        let back = import.document.model();
        let [solid] = import.solids.as_slice() else {
            panic!("{name}: {} solids came back", import.solids.len());
        };
        let diagnosis = ogeom::algo::check(back, solid, T).unwrap();
        assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
        let (same, total) = pcurves::kept((&model, &shape), (back, solid), T);
        assert!(total > 0, "{name}: the read solid holds pcurves");
        assert_eq!(same, total, "{name}: pcurves come back as written");
        let after = ogeom::algo::volume_properties(back, solid, fine, T)
            .unwrap()
            .mass;
        // Both measured at the same fine chord; the same geometry and the
        // same trims leave only where the meshes place their points.
        assert!(
            (after - before).abs() <= before * 1e-6,
            "{name}: {before} went out, {after} came back"
        );
    }
}

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

/// A drum an inch round and two tall, written in inches and degrees with
/// every pcurve its faces take: the wall's in degrees round and inches up
/// (its seam's second side at 360), each cap's a circle an inch round, and
/// the bottom cap's running the other way round from its edge's circle and
/// a twelfth of a turn on, as the cap's frame faces down and is turned. Every
/// one is kept, the solid is valid with the drum's volume, and written back
/// out in millimetres the pcurves the writer states come back kept too.
#[test]
fn pcurves_in_inches_and_degrees_are_kept() {
    let text = corpus("drum_in_inches_and_degrees.step");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let report = &import.report;
    assert_eq!(report.scale_mm, 25.4);
    // The wall takes three edges, each cap one.
    assert_eq!(
        (report.pcurves_kept, report.pcurves_derived),
        (5, 0),
        "every stated pcurve is kept"
    );
    let [solid] = import.solids.as_slice() else {
        panic!("{} solids", import.solids.len());
    };
    let model = import.document.model();
    let diagnosis = ogeom::algo::check(model, solid, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let fine = Deflection::with_chord(1e-4).unwrap();
    let volume = |model: &Model, solid| {
        ogeom::algo::volume_properties(model, solid, fine, T)
            .unwrap()
            .mass
    };
    let drum = core::f64::consts::PI * 25.4 * 25.4 * 50.8;
    // A chord of 1e-4 mm on a 25.4 mm circle leaves the caps' polygons
    // short of the disc by under 1e-5 of it.
    let read = volume(model, solid);
    assert!((read - drum).abs() <= drum * 1e-5, "{read} against {drum}");

    let text = ogeom::io::write_step(&import.document, T).unwrap();
    let again = ogeom::io::read_step(&text, T).unwrap();
    let [back] = again.solids.as_slice() else {
        panic!("{} solids came back", again.solids.len());
    };
    let back_model = again.document.model();
    assert!(ogeom::algo::check(back_model, back, T).unwrap().is_valid());
    assert!(again.report.pcurves_kept > 0, "the writer states pcurves");
    assert_eq!(again.report.pcurves_derived, 0, "and each comes back kept");
    let (same, total) = pcurves::kept((model, solid), (back_model, back), T);
    assert_eq!(same, total, "the pcurves come back as read");
    // The same geometry and the same trims, meshed at the same chord.
    let after = volume(back_model, back);
    assert!(
        (after - read).abs() <= read * 1e-9,
        "{read} out, {after} back"
    );
}

/// A box and a drum written with their pcurves, each pcurve line's vector
/// then doubled: the file states every boundary at half the pace of its
/// edge's curve, as a writer does that runs its pcurves over a parameter
/// of their own. Each is the same boundary at another pace, and is kept.
#[test]
fn pcurves_at_another_pace_than_their_curve_are_kept() {
    let mut model = Model::new();
    let parts = [
        ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 10.0, 5.0), T)
            .unwrap()
            .shape,
        ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 3.0, 5.0, T)
            .unwrap()
            .shape,
    ];
    for shape in parts {
        let mut document = ogeom::doc::Document::over(model.clone());
        document.add_part("part", shape.clone());
        let text = ogeom::io::write_step(&document, T).unwrap();
        let written = ogeom::io::read_step(&text, T).unwrap().report.pcurves_kept;
        assert!(written > 0, "the writer states pcurves");
        // The ids of the chart's directions: two components.
        let planar: Vec<String> = text
            .lines()
            .filter(|l| l.contains("=DIRECTION('',(") && l.matches(',').count() == 2)
            .map(|l| l[..l.find('=').unwrap()].to_owned())
            .collect();
        let mut doubled = 0;
        let paced: String = text
            .lines()
            .map(|line| {
                if line.contains("=VECTOR('',")
                    && planar.iter().any(|id| line.contains(&format!(",{id},")))
                {
                    doubled += 1;
                    format!("{},2.0);\n", &line[..line.rfind(',').unwrap()])
                } else {
                    format!("{line}\n")
                }
            })
            .collect();
        assert!(doubled > 0, "the fixture paced its pcurve lines");
        let import = ogeom::io::read_step(&paced, T).unwrap();
        assert_eq!(
            (import.report.pcurves_kept, import.report.pcurves_derived),
            (written, 0),
            "every pcurve is kept at its pace"
        );
        let [solid] = import.solids.as_slice() else {
            panic!("{} solids", import.solids.len());
        };
        let back = import.document.model();
        assert!(ogeom::algo::check(back, solid, T).unwrap().is_valid());
        let fine = Deflection::with_chord(1e-3).unwrap();
        let before = ogeom::algo::volume_properties(&model, &shape, fine, T)
            .unwrap()
            .mass;
        let after = ogeom::algo::volume_properties(back, solid, fine, T)
            .unwrap()
            .mass;
        assert!(
            (after - before).abs() <= before * 1e-9,
            "{before} out, {after} back"
        );
    }
}
