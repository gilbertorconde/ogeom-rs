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
