//! What a reader builds keeps the tolerance containment rule, what the
//! healer is handed it restores, and what the IGES writer writes its reader
//! reads back.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{Severity, check, restore_containment};
use ogeom::core::{Tolerance, Tolerances};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn corpus(name: &str) -> String {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).expect("the corpus file is committed")
}

fn containment_findings(model: &Model, shape: &Shape) -> usize {
    check(model, shape, T)
        .unwrap()
        .of(Severity::Broken)
        .iter()
        .filter(|p| p.what.contains("tighter than the"))
        .count()
}

/// The widest a pcurve, lifted through its surface, leaves its edge's
/// curve beyond the edge's tolerance, with that tolerance; `None` where
/// every pcurve of every edge keeps within it.
///
/// Measured here, not through `check`: each pcurve is read at 501 points
/// spread evenly over its range, a lifted point is compared with the
/// curve's point at the same fraction of its range, and where that is
/// farther than the tolerance, with the nearest point of the whole curve
/// (a dense scan refined by golden section), since a pcurve need not share
/// its curve's pace. A pcurve placed elsewhere than its curve describes
/// another occurrence and is not compared.
fn widest_pcurve_excess(model: &Model, shape: &Shape) -> Option<(f64, f64)> {
    use ogeom::geom::{Curve2d as _, Curve3d as _, Surface as _};
    use ogeom::topo::EdgeRepr;
    const SAMPLES: u32 = 500;
    let mut widest: Option<(f64, f64)> = None;
    for edge in explore_unique(model, shape, ShapeType::Edge).unwrap() {
        let data = model.node(&edge).unwrap().data().as_edge().unwrap();
        let Some(EdgeRepr::Curve3d {
            curve,
            range,
            location,
        }) = data.curve3d()
        else {
            continue;
        };
        let curve = model.geometry().curve(*curve).unwrap();
        let stated = data.tolerance.get();
        let nearest = |p: ogeom::math::Point| {
            let at = |t: f64| curve.point_at(t, T).unwrap().distance(p);
            let step = (range.1 - range.0) / f64::from(SAMPLES);
            let (mut best, mut least) = (range.0, at(range.0));
            for k in 1..=SAMPLES {
                let t = range.0 + step * f64::from(k);
                if at(t) < least {
                    (best, least) = (t, at(t));
                }
            }
            let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
            let (mut a, mut b) = ((best - step.abs()).max(lo), (best + step.abs()).min(hi));
            for _ in 0..60 {
                let (c, d) = (b - (b - a) * 0.618, a + (b - a) * 0.618);
                if at(c) < at(d) { b = d } else { a = c }
            }
            least.min(at(0.5 * (a + b)))
        };
        for repr in &data.representations {
            let (sides, prange, surface, at) = match repr {
                EdgeRepr::PCurve {
                    curve,
                    range,
                    surface,
                    location,
                } => (vec![*curve], *range, *surface, location),
                EdgeRepr::Seam {
                    forward,
                    reversed,
                    range,
                    surface,
                    location,
                } => (vec![*forward, *reversed], *range, *surface, location),
                _ => continue,
            };
            if at != location {
                continue;
            }
            let surface = model.geometry().surface(surface).unwrap();
            for id in sides {
                let pcurve = model.geometry().pcurve(id).unwrap();
                for k in 0..=SAMPLES {
                    let f = f64::from(k) / f64::from(SAMPLES);
                    let uv = pcurve
                        .point_at(prange.0 + (prange.1 - prange.0) * f, T)
                        .unwrap();
                    let Ok(lifted) = surface.point_at(uv.x, uv.y, T) else {
                        continue;
                    };
                    let on_curve = curve
                        .point_at(range.0 + (range.1 - range.0) * f, T)
                        .unwrap();
                    let mut gap = on_curve.distance(lifted);
                    if gap > stated {
                        gap = gap.min(nearest(lifted));
                    }
                    if gap > stated && widest.is_none_or(|(w, s)| gap - stated > w - s) {
                        widest = Some((gap, stated));
                    }
                }
            }
        }
    }
    widest
}

/// A reader states, on every edge, how far its pcurves stand from its
/// curve: an exact pcurve of a line lying microns off its plane, and a
/// pcurve fitted on a sphere straying between its samples, both widen
/// their edge. Through STEP and through IGES, from a file and from a
/// boolean's result written out.
#[test]
fn a_read_states_how_far_its_pcurves_stand_from_their_curves() {
    for name in ["nist_ctc_03_asme1_rc.stp", "nist_ftc_07_asme1_rd.stp"] {
        let read = ogeom::io::read_step(&corpus(name), T).unwrap();
        for solid in &read.solids {
            let excess = widest_pcurve_excess(read.document.model(), solid);
            assert!(
                excess.is_none(),
                "{name} through STEP: (gap, stated) {excess:?}"
            );
        }
        let iges = ogeom::io::write_iges(&read.document, T).unwrap();
        let back = ogeom::io::read_iges(&iges, T).unwrap();
        for solid in &back.solids {
            let excess = widest_pcurve_excess(back.document.model(), solid);
            assert!(
                excess.is_none(),
                "{name} through IGES: (gap, stated) {excess:?}"
            );
        }
    }

    // A sphere cut on its side: the section circle has no closed-form
    // image on the sphere, so each reader fits one.
    let frame = |x: f64, y: f64, z: f64| {
        ogeom::math::Frame::new(
            ogeom::math::Point::new(x, y, z),
            ogeom::math::Direction::Z,
            ogeom::math::Direction::X,
            T,
        )
        .unwrap()
    };
    let mut model = Model::new();
    let sphere = ogeom::algo::make_sphere(&mut model, frame(0.0, 0.0, 0.0), 10.0, T)
        .unwrap()
        .shape;
    let block = ogeom::algo::make_box(&mut model, frame(5.0, -20.0, -20.0), (40.0, 40.0, 40.0), T)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(&mut model, &sphere, &block, T)
        .unwrap()
        .shape;
    let mut document = ogeom::doc::Document::over(model);
    document.add_part("part", cut);
    let step = ogeom::io::read_step(&ogeom::io::write_step(&document, T).unwrap(), T).unwrap();
    let iges = ogeom::io::read_iges(&ogeom::io::write_iges(&document, T).unwrap(), T).unwrap();
    for (format, read) in [
        ("STEP", (step.document.model(), &step.solids[0])),
        ("IGES", (iges.document.model(), &iges.solids[0])),
    ] {
        let excess = widest_pcurve_excess(read.0, read.1);
        assert!(
            excess.is_none(),
            "sphere cut through {format}: (gap, stated) {excess:?}"
        );
    }
}

fn volume(model: &Model, shape: &Shape) -> Option<f64> {
    ogeom::algo::volume_properties(model, shape, Deflection::default(), T)
        .ok()
        .map(|p| p.mass)
}

/// The STEP reader widens an edge to how far its pcurves sit from its
/// curve. The edge's vertices widen with it, so a freshly read solid holds
/// every vertex at least as loose as the edges it bounds. These files read
/// with hundreds of vertices tighter than their edges otherwise.
#[test]
fn a_step_read_keeps_tolerance_containment() {
    for name in [
        "nist_ctc_02_asme1_rc.stp",
        "nist_ctc_05_asme1_rd.stp",
        "long_bore_between_cross_holes.step",
    ] {
        let import = ogeom::io::read_step(&corpus(name), T).unwrap();
        assert!(!import.solids.is_empty(), "{name}");
        for solid in &import.solids {
            assert_eq!(
                containment_findings(import.document.model(), solid),
                0,
                "{name}: a vertex tighter than an edge it bounds"
            );
        }
    }
}

/// The healer restores containment. A real part whose edges genuinely
/// need their width (pcurves sitting microns off their curves) has its
/// vertices reset to the confusion tolerance, which is the state a reader
/// that widened only the edges left behind. The reduction cannot tighten
/// those edges, so the vertices are what must grow, and the report says
/// how many did.
#[test]
fn fix_shape_restores_tolerance_containment() {
    let mut import = ogeom::io::read_step(&corpus("nist_ctc_02_asme1_rc.stp"), T).unwrap();
    let solid = import.solids[0].clone();
    let model = import.document.model_mut();
    let tight = Tolerance::new(T.confusion()).unwrap();
    for vertex in explore_unique(model, &solid, ShapeType::Vertex).unwrap() {
        if let Some(node) = model.node_mut(&vertex)
            && let NodeData::Vertex(data) = node.data_mut()
        {
            data.tolerance = tight;
        }
    }
    assert!(containment_findings(model, &solid) > 0);

    let fixed = ogeom::heal::fix_shape(model, &solid, T).unwrap();
    assert_eq!(containment_findings(model, &fixed.shape), 0);
    assert!(fixed.report.tolerances_widened > 0);
    assert!(
        fixed.report.after.of(Severity::Broken).is_empty(),
        "{}",
        fixed.report.after
    );

    // The reduction shrinks no edge below how far its pcurves stand from
    // its curve: this part's fitted trims stray most between the points a
    // coarse measure would read.
    let excess = widest_pcurve_excess(model, &fixed.shape);
    assert!(excess.is_none(), "(gap, stated) {excess:?}");

    // The pass alone, run again, has nothing left to grow.
    assert_eq!(restore_containment(model, &fixed.shape).unwrap(), 0);
}

/// Whatever the IGES writer writes, its reader reads: every record exactly
/// eighty columns, and every solid back, valid, at the volume it went out
/// with. One of these parts carries a coefficient near 1e-51, which
/// positional notation spells longer than a record. The others have
/// vertices their curves miss by a few nanometres, or project a hair past
/// a bounded curve's end.
#[test]
fn a_solid_written_as_iges_reads_back() {
    for name in [
        "nist_ctc_03_asme1_rc.stp",
        "nist_ctc_01_asme1_rd.stp",
        "nist_ftc_07_asme1_rd.stp",
        "closed_tube_face_crosses_its_join.step",
        "fillet_strip_with_a_narrow_chart.step",
    ] {
        let read = ogeom::io::read_step(&corpus(name), T).unwrap();
        let iges = ogeom::io::write_iges(&read.document, T).unwrap();
        assert!(
            iges.lines().all(|line| line.len() == 80),
            "{name}: a record not eighty columns"
        );
        let back = ogeom::io::read_iges(&iges, T).unwrap();
        assert_eq!(back.solids.len(), read.solids.len(), "{name}");
        for (sent, came) in read.solids.iter().zip(&back.solids) {
            let model = back.document.model();
            assert!(check(model, came, T).unwrap().is_usable(), "{name}");
            assert_eq!(containment_findings(model, came), 0, "{name}");
            let faces = |m: &Model, s: &Shape| explore_unique(m, s, ShapeType::Face).unwrap().len();
            assert_eq!(
                faces(read.document.model(), sent),
                faces(model, came),
                "{name}"
            );
            // A part that encloses a volume encloses the same one.
            if let (Some(a), Some(b)) = (volume(read.document.model(), sent), volume(model, came)) {
                assert!(
                    (a - b).abs() <= a.abs() * 1e-3,
                    "{name}: volume {a} went out, {b} came back"
                );
            }
        }
    }
}
