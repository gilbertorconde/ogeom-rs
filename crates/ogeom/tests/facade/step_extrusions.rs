//! A STEP file's surface of linear extrusion: a curve swept along a
//! vector, which a file writes for a drum's wall as often as it writes a
//! cylinder. Read as what it is (a cylinder where a circle is swept along
//! its own axis, the swept surface itself otherwise) and measured.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::geom::SurfaceGeometry;
use ogeom::math::Frame;
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

/// A drum written by this kernel, its wall restated as `swept` (a curve
/// entity spelled over the wall's own placement) swept along the drum's
/// axis. The substitution keeps every other entity as written.
fn drum_with_wall_swept(swept: &str) -> String {
    let mut model = Model::new();
    let solid = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T)
        .unwrap()
        .shape;
    let mut document = ogeom::doc::Document::over(std::mem::take(&mut model));
    document.add_part("drum", solid);
    let text = ogeom::io::write_step(&document, T).unwrap();
    let wall = text
        .lines()
        .find(|l| l.contains("CYLINDRICAL_SURFACE"))
        .expect("the drum has a cylindrical wall");
    let (id, rest) = wall.split_once('=').unwrap();
    let inner = rest
        .trim_start_matches("CYLINDRICAL_SURFACE('',")
        .trim_end_matches(");");
    let (placement, radius) = inner.split_once(',').unwrap();
    let axis = text
        .lines()
        .find(|l| l.starts_with(&format!("{placement}=AXIS2_PLACEMENT_3D")))
        .and_then(|l| l.split(',').nth(2))
        .expect("the placement names its axis")
        .to_owned();
    let top = text
        .lines()
        .filter_map(|l| l.strip_prefix('#')?.split('=').next()?.parse::<u64>().ok())
        .max()
        .unwrap();
    let (curve, vector) = (top + 1, top + 2);
    let restated = format!("{id}=SURFACE_OF_LINEAR_EXTRUSION('',#{curve},#{vector});");
    let appended = format!(
        "#{curve}={};\n#{vector}=VECTOR('',{axis},1.);\nENDSEC;",
        swept
            .replace("{placement}", placement)
            .replace("{radius}", radius)
    );
    let last = text.rfind("ENDSEC;").unwrap();
    let mut out = text.replace(wall, &restated);
    let last = out.rfind("ENDSEC;").unwrap_or(last);
    out.replace_range(last..last + "ENDSEC;".len(), &appended);
    out
}

fn wall_kinds(import: &ogeom::io::StepImport) -> Vec<&'static str> {
    let model = import.document.model();
    explore_unique(model, &import.solids[0], ShapeType::Face)
        .unwrap()
        .iter()
        .map(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            match model.geometry().surface(data.surface) {
                Some(SurfaceGeometry::Plane(_)) => "plane",
                Some(SurfaceGeometry::Cylinder(_)) => "cylinder",
                Some(SurfaceGeometry::Extrusion(_)) => "extrusion",
                _ => "other",
            }
        })
        .collect()
}

/// A circle swept along its own axis is a cylinder, and reads as one:
/// exact, and known to every operation downstream.
#[test]
fn a_circle_swept_along_its_axis_reads_as_the_cylinder_it_is() {
    let text = drum_with_wall_swept("CIRCLE('',{placement},{radius})");
    let import = ogeom::io::read_step(&text, T).unwrap();
    assert_eq!(import.solids.len(), 1);
    let kinds = wall_kinds(&import);
    assert_eq!(
        kinds.iter().filter(|k| **k == "cylinder").count(),
        1,
        "{kinds:?}"
    );
    let expected = core::f64::consts::PI * 25.0 * 10.0;
    let measured = volume(import.document.model(), &import.solids[0]);
    assert!(
        (measured - expected).abs() < expected * 1e-6,
        "{measured} against {expected}"
    );
}

/// Anything else sweeps as itself (an ellipse here, drawn round to the
/// same drum so the volume is known), and the face it bounds draws and
/// measures, its trims fitted by projection over a window the face's own
/// edges size.
#[test]
fn an_ellipse_swept_reads_as_an_extrusion_that_measures() {
    let text = drum_with_wall_swept("ELLIPSE('',{placement},{radius},{radius})");
    let import = ogeom::io::read_step(&text, T).unwrap();
    assert_eq!(import.solids.len(), 1);
    let kinds = wall_kinds(&import);
    assert_eq!(
        kinds.iter().filter(|k| **k == "extrusion").count(),
        1,
        "{kinds:?}"
    );
    assert!(
        !import.report.warnings.iter().any(|w| w.contains("skipped")),
        "{:?}",
        import.report.warnings
    );
    let expected = core::f64::consts::PI * 25.0 * 10.0;
    let measured = volume(import.document.model(), &import.solids[0]);
    assert!(
        (measured - expected).abs() < expected * 1e-3,
        "{measured} against {expected}"
    );
}
