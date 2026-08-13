//! STEP edges dressed as `SURFACE_CURVE`: what every exporter derived from
//! the reference kernel writes.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::topo::Model;
use ogeom_core::Tolerances;
use ogeom_math::Frame;
use ogeom_mesh::Deflection;

const T: Tolerances = Tolerances::millimetres();

/// Re-dress every `EDGE_CURVE`'s geometry in a `SURFACE_CURVE` wrapper, the
/// way OCCT-derived exporters write their files.
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
        wrapped.contains("SURFACE_CURVE"),
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
