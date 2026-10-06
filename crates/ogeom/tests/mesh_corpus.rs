//! The STEP corpus meshed and converted back: meshes whose answer is known,
//! the part they were drawn from. Each part is meshed at a thousandth of its
//! diagonal and rebuilt with the default options.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    FallbackReason, MeshSolidOptions, MeshSolidReport, check, shape_bounds, solid_from_mesh,
    volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{Curve, Curve2d as _, Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom::mesh::Deflection;
use ogeom::topo::{EdgeRepr, Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// How many of a shape's faces are on a surface other than a plane.
fn curved(model: &Model, shape: &Shape) -> usize {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .iter()
        .filter(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            !matches!(
                model.geometry().surface(data.surface),
                Some(SurfaceGeometry::Plane(_))
            )
        })
        .count()
}

/// How many edges of a shape's faces, each meshed on its own with the
/// edges drawn alike for every face, are used other than twice: where a
/// face's boundary leaves a gap its neighbours do not fill, with no weld to
/// close it. A triangle with two corners at one point (at an apex, where
/// the chart's row lifts to the apex) covers nothing and is not counted.
fn unmatched_face_edges(model: &Model, shape: &Shape) -> usize {
    let deflection = Deflection::default();
    let chords = ogeom::mesh::edge_chords_for(model, shape, deflection, T).unwrap();
    let key = |p: &ogeom::math::Point| (p.x.to_bits(), p.y.to_bits(), p.z.to_bits());
    let mut uses: std::collections::HashMap<_, usize> = std::collections::HashMap::new();
    for face in explore_unique(model, shape, ShapeType::Face).unwrap() {
        let mesh =
            ogeom::mesh::triangulate_face_with(model, &face, deflection, &chords, T).unwrap();
        for t in &mesh.triangles {
            let [p, q, r] = t.map(|i| key(&mesh.positions[i as usize]));
            if p == q || q == r || r == p {
                continue;
            }
            for k in 0..3 {
                let (a, b) = (
                    key(&mesh.positions[t[k] as usize]),
                    key(&mesh.positions[t[(k + 1) % 3] as usize]),
                );
                *uses.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
    }
    uses.values().filter(|&&n| n != 2).count()
}

/// The distance from `p` to the nearest point of `curve` over `range`: a
/// scan, then a golden-section search about the nearest sample.
fn nearest_on(curve: &Curve, range: (f64, f64), p: ogeom::math::Point) -> f64 {
    let at = |t: f64| curve.point_at(t, T).unwrap().distance(p);
    let step = (range.1 - range.0) / 400.0;
    let mut best = (range.0, at(range.0));
    for k in 1..=400 {
        let t = range.0 + step * f64::from(k);
        if at(t) < best.1 {
            best = (t, at(t));
        }
    }
    let (lo, hi) = (range.0.min(range.1), range.0.max(range.1));
    let (mut a, mut b) = ((best.0 - step.abs()).max(lo), (best.0 + step.abs()).min(hi));
    let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
    for _ in 0..100 {
        let (c, d) = (b - (b - a) * ratio, a + (b - a) * ratio);
        if at(c) < at(d) {
            b = d;
        } else {
            a = c;
        }
    }
    best.1.min(at(f64::midpoint(a, b)))
}

/// The edges whose pcurves, lifted through their surfaces, leave their
/// curve by more than the edge's tolerance anywhere along it: how far, and
/// the tolerance. Each pcurve is sampled along its range, and each lifted
/// point measured against the curve's point at the same fraction of its
/// range or, where that is farther than the tolerance, the nearest point
/// of the curve.
fn pcurves_off_their_curves(model: &Model, shape: &Shape) -> Vec<(f64, f64)> {
    let mut off = Vec::new();
    for edge in explore_unique(model, shape, ShapeType::Edge).unwrap() {
        let data = model.node(&edge).unwrap().data().as_edge().unwrap();
        let Some(EdgeRepr::Curve3d { curve, range, .. }) = data.curve3d() else {
            continue;
        };
        let curve = model.geometry().curve(*curve).unwrap();
        let tolerance = data.tolerance.get();
        for repr in &data.representations {
            let (pcurves, pcurve_range, surface) = match repr {
                EdgeRepr::PCurve {
                    curve,
                    range,
                    surface,
                    ..
                } => (vec![*curve], *range, *surface),
                EdgeRepr::Seam {
                    forward,
                    reversed,
                    range,
                    surface,
                    ..
                } => (vec![*forward, *reversed], *range, *surface),
                _ => continue,
            };
            let surface = model.geometry().surface(surface).unwrap();
            for id in pcurves {
                let pcurve = model.geometry().pcurve(id).unwrap();
                let mut widest: f64 = 0.0;
                for k in 0..=64 {
                    let f = f64::from(k) / 64.0;
                    let uv = pcurve
                        .point_at(pcurve_range.0 + (pcurve_range.1 - pcurve_range.0) * f, T)
                        .unwrap();
                    let Ok(lifted) = surface.point_at(uv.x, uv.y, T) else {
                        continue;
                    };
                    let mut gap = curve
                        .point_at(range.0 + (range.1 - range.0) * f, T)
                        .unwrap()
                        .distance(lifted);
                    if gap > tolerance {
                        gap = gap.min(nearest_on(curve, *range, lifted));
                    }
                    widest = widest.max(gap);
                }
                if widest > tolerance {
                    off.push((widest, tolerance));
                }
            }
        }
    }
    off
}

/// The corpus part `name` meshed and converted back comes out a valid solid
/// of its volume, tessellating closed where `closed` says it does, its faces
/// meshed one by one meeting edge to edge where `meet` says they do, every
/// pcurve within its edge's tolerance of the edge's curve, and with curved
/// faces where the part has them. Returns what the conversion reported.
fn comes_back(name: &str, closed: bool, meet: bool) -> MeshSolidReport {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(path).expect("the corpus file is committed");
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let part = &import.solids[0];
    let diagonal = shape_bounds(model, part, T).unwrap().diagonal();
    let mesh = ogeom::mesh::triangulate(
        model,
        part,
        Deflection::with_chord(diagonal * 1e-3).unwrap(),
        T,
    )
    .unwrap();
    assert!(mesh.is_closed(), "{name}: the source mesh is open");
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
    assert!(out.closed, "{name}: {:?}", out.report);
    let diagnosis = check(&back, &out.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{name}: {diagnosis}");
    let off = pcurves_off_their_curves(&back, &out.shape);
    assert!(
        off.is_empty(),
        "{name}: pcurves leave their curves, (gap, tolerance): {off:?}"
    );
    let fine = Deflection::with_chord(diagonal * 1e-5).unwrap();
    let (want, got) = (
        volume_properties(model, part, fine, T).unwrap().mass,
        volume_properties(&back, &out.shape, fine, T).unwrap().mass,
    );
    assert!(
        (got - want).abs() <= want * 3e-3,
        "{name}: {want} drawn, {got} came back"
    );
    if closed {
        let drawn = ogeom::mesh::triangulate(&back, &out.shape, Deflection::default(), T).unwrap();
        assert!(
            drawn.is_closed(),
            "{name}: the converted part tessellates open"
        );
    }
    if meet {
        assert_eq!(
            unmatched_face_edges(&back, &out.shape),
            0,
            "{name}: its faces meshed one by one do not meet"
        );
    }
    if curved(model, part) > 0 {
        assert!(
            curved(&back, &out.shape) > 0,
            "{name}: nothing curved came back"
        );
    }
    out.report
}

#[test]
fn nist_ctc_01() {
    comes_back("nist_ctc_01_asme1_rd.stp", true, true);
}

#[test]
fn nist_ctc_03() {
    comes_back("nist_ctc_03_asme1_rc.stp", true, true);
}

/// Its blends come back as tori fitted side by side, meeting at so slight
/// an angle that where the two surfaces cross lies millimetres from the
/// mesh's boundary between them. The edge there is threaded along the
/// boundary through points resting on both, and keeps within twenty
/// coplanar distances of each; threaded through the boundary's vertices
/// alone it stands off them by its long spans' sag. Facets round its large
/// bore are chords whose planes meet the bore past their far corners; they
/// come back as fans, and the bore curved.
#[test]
#[ignore = "heavy"]
fn nist_ctc_02() {
    let name = "nist_ctc_02_asme1_rc.stp";
    let report = comes_back(name, true, true);
    let overlaps = report
        .fallbacks
        .iter()
        .filter(|f| f.reason == FallbackReason::Overlaps)
        .count();
    assert!(overlaps <= 2, "{overlaps} faces faceted for overlaps");
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    let import = ogeom::io::read_step(&std::fs::read_to_string(path).unwrap(), T).unwrap();
    let model = import.document.model();
    let part = &import.solids[0];
    let diagonal = shape_bounds(model, part, T).unwrap().diagonal();
    let mesh = ogeom::mesh::triangulate(
        model,
        part,
        Deflection::with_chord(diagonal * 1e-3).unwrap(),
        T,
    )
    .unwrap();
    let mut back = Model::new();
    let out = solid_from_mesh(&mut back, &mesh, &MeshSolidOptions::default(), T).unwrap();
    let is_curved = |face: &Shape| {
        let data = back.node(face).unwrap().data().as_face().unwrap();
        !matches!(
            back.geometry().surface(data.surface),
            Some(SurfaceGeometry::Plane(_))
        )
    };
    let mut owners: std::collections::HashMap<_, Vec<bool>> = std::collections::HashMap::new();
    for face in explore_unique(&back, &out.shape, ShapeType::Face).unwrap() {
        for edge in explore_unique(&back, &face, ShapeType::Edge).unwrap() {
            owners
                .entry(edge.node())
                .or_default()
                .push(is_curved(&face));
        }
    }
    let worst = explore_unique(&back, &out.shape, ShapeType::Edge)
        .unwrap()
        .iter()
        .filter(|e| owners[&e.node()] == [true, true])
        .map(|e| {
            back.node(e)
                .unwrap()
                .data()
                .as_edge()
                .unwrap()
                .tolerance
                .get()
        })
        .fold(0.0_f64, f64::max);
    assert!(
        worst <= out.coplanar_distance * 20.0,
        "an edge between curved faces stands {worst} off them"
    );
}

/// A fitted face here ran past a facet it should have ended on, and the
/// solid came back with the facet facing into material.
#[test]
#[ignore = "heavy"]
fn nist_ctc_04() {
    comes_back("nist_ctc_04_asme1_rd.stp", true, true);
}

/// Its whole mesh closes, welded; its faces meshed one by one still leave
/// six edges unmatched.
#[test]
#[ignore = "heavy"]
fn nist_ftc_06() {
    comes_back("nist_ftc_06_asme1_rd.stp", true, false);
}

/// As `nist_ctc_04`, on tori.
#[test]
#[ignore = "heavy"]
fn nist_ftc_07() {
    comes_back("nist_ftc_07_asme1_rd.stp", true, true);
}

#[test]
fn nist_ftc_08() {
    comes_back("nist_ftc_08_asme1_rc.stp", true, true);
}

#[test]
fn nist_ftc_09() {
    comes_back("nist_ftc_09_asme1_rd.stp", true, true);
}

/// Its seams end a little off the corners they meet at, each on its own
/// side, within the corners' tolerance; drawn from the corners themselves,
/// the faces round each one meet at one point and the mesh closes. Facets
/// beside a cylinder meet it all but tangentially, and built flat they
/// turn into the material; built as fans onto the cylinder, it stays
/// curved.
#[test]
#[ignore = "heavy"]
fn nist_ftc_10() {
    let report = comes_back("nist_ftc_10_asme1_rb.stp", true, true);
    assert!(report.faces <= 700, "{} faces came back", report.faces);
    assert!(
        !report
            .fallbacks
            .iter()
            .any(|f| f.reason == FallbackReason::Overlaps),
        "a curved face was faceted for a facet turned in beside it: {:?}",
        report.fallbacks
    );
}

#[test]
fn nist_ftc_11() {
    comes_back("nist_ftc_11_asme1_rb.stp", true, true);
}

#[test]
#[ignore = "heavy"]
fn a_socket_head_screw() {
    comes_back("m5x16_bhcs_loops.step", true, true);
}
