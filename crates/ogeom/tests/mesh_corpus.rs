//! The STEP corpus meshed and converted back: meshes whose answer is known,
//! the part they were drawn from. Each part is meshed at a thousandth of its
//! diagonal and rebuilt with the default options.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    Canonical, FallbackReason, MeshRegion, MeshRegions, MeshSolidOptions, MeshSolidReport, check,
    shape_bounds, solid_from_mesh, volume_properties,
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

/// Its rounds of radius 9.5 to 12.7 are strips a few rows wide, whose
/// axis the cylinder fit can miss; a torus thousands in radius fits them
/// as closely, and they come back on the cylinder it stands in for.
#[test]
fn nist_ctc_03() {
    let name = "nist_ctc_03_asme1_rc.stp";
    comes_back(name, true, true);
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
    let huge: Vec<f64> = explore_unique(&back, &out.shape, ShapeType::Face)
        .unwrap()
        .iter()
        .filter_map(|face| {
            let data = back.node(face).unwrap().data().as_face().unwrap();
            match back.geometry().surface(data.surface) {
                Some(SurfaceGeometry::Torus(t)) if t.torus().major_radius() > diagonal => {
                    Some(t.torus().major_radius())
                }
                _ => None,
            }
        })
        .collect();
    assert!(huge.is_empty(), "tori standing in for cylinders: {huge:?}");
}

/// Its blends come back as tori fitted side by side, meeting at so slight
/// an angle that where the two surfaces cross lies millimetres from the
/// mesh's boundary between them. The edge there is threaded along the
/// boundary through points resting on both, and keeps within twenty
/// coplanar distances of each; threaded through the boundary's vertices
/// alone it stands off them by its long spans' sag. A torus and a cylinder
/// fitted to pieces of one spline blend can meet all but tangentially with
/// no point carried onto both, and the line between them may stand looser,
/// to its chord's sag; only edges between two tori are held here. Facets
/// round its large bore are chords whose planes meet the bore past their
/// far corners; they come back as fans, and the bore curved. A facet with
/// two sides on curved faces, turned in, comes back as a wedge ruled
/// between its two seams, and both curved faces stay curved.
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
    assert_eq!(overlaps, 0, "{overlaps} faces faceted for overlaps");
    assert!(
        report.curved_faceted <= 2,
        "{} curved faces faceted: {:?}",
        report.curved_faceted,
        report.fallbacks
    );
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
    let is_torus = |face: &Shape| {
        let data = back.node(face).unwrap().data().as_face().unwrap();
        matches!(
            back.geometry().surface(data.surface),
            Some(SurfaceGeometry::Torus(_))
        )
    };
    let mut owners: std::collections::HashMap<_, Vec<bool>> = std::collections::HashMap::new();
    for face in explore_unique(&back, &out.shape, ShapeType::Face).unwrap() {
        for edge in explore_unique(&back, &face, ShapeType::Edge).unwrap() {
            owners.entry(edge.node()).or_default().push(is_torus(&face));
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
        "an edge between two tori stands {worst} off them"
    );
}

/// A fitted face here ran past a facet it should have ended on, and the
/// solid came back with the facet facing into material. A round along a
/// plane's edge turns round a rounded corner on a torus meshed a few facets
/// round, whose pieces lie on spheres as exactly as on the torus; each is
/// put on the torus the round turns on, and nothing is faceted.
#[test]
#[ignore = "heavy"]
fn nist_ctc_04() {
    let report = comes_back("nist_ctc_04_asme1_rd.stp", true, true);
    assert_eq!(report.curved_faceted, 0, "{:?}", report.fallbacks);
}

/// A round of radius 6.35 runs out tangentially into a plane, and the mesh
/// draws the line they meet on with three vertices: too few for a cubic,
/// and the edge there is the parabola through them.
#[test]
#[ignore = "heavy"]
fn nist_ctc_05() {
    let report = comes_back("nist_ctc_05_asme1_rd.stp", true, true);
    assert_eq!(report.curved_faceted, 0, "{:?}", report.fallbacks);
}

/// Its whole mesh closes, welded; its faces meshed one by one still leave
/// six edges unmatched. Its rounds meet their neighbours all but
/// tangentially, and the edges there are threaded through points carried
/// onto both surfaces between the mesh's vertices; through the vertices
/// alone no curve keeps to both.
#[test]
#[ignore = "heavy"]
fn nist_ftc_06() {
    let report = comes_back("nist_ftc_06_asme1_rd.stp", true, false);
    assert_eq!(report.curved_faceted, 0, "{:?}", report.fallbacks);
    assert!(report.faces <= 200, "{} faces came back", report.faces);
}

/// As `nist_ctc_04`, on tori. Its long rounds of radius 0.43 are meshed a
/// few rows across, with fans of long facets from single corners; each
/// strip of them comes back on its round's cylinder, where a sphere tens
/// of thousands in radius through the strip would lay the facets' cross
/// chords on itself and fall back to facets. The corners of its drafted
/// walls are cones meshed a few facets round, too few vertices to fit
/// alone; the walls tangent to them fix them, and the tori at their foot
/// seam to them. Its drill points are cones closing at their apex.
#[test]
#[ignore = "heavy"]
fn nist_ftc_07() {
    let report = comes_back("nist_ftc_07_asme1_rd.stp", true, true);
    assert_eq!(report.curved_faceted, 0, "{:?}", report.fallbacks);
}

/// Its fillet tori are meshed a few facets round, and a piece of one can
/// lie on a sphere as exactly as on the torus. Where the rest of the torus
/// is a larger region beside it, the piece joins that region; where the
/// torus is all such pieces, each is put on the torus between the plane
/// and the cylinder it blends, and they are one region.
#[test]
fn nist_ftc_08() {
    let report = comes_back("nist_ftc_08_asme1_rc.stp", true, true);
    assert_eq!(report.curved_faceted, 0, "{:?}", report.fallbacks);
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
/// curved. Its drill points close at their apex inside rims whose circles
/// start off the cone's own seam, and the cone is turned about its axis to
/// meet them.
#[test]
#[ignore = "heavy"]
fn nist_ftc_10() {
    let report = comes_back("nist_ftc_10_asme1_rb.stp", true, true);
    assert!(report.faces <= 430, "{} faces came back", report.faces);
    assert!(report.curved_faceted <= 1, "{:?}", report.fallbacks);
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

/// Its head's flat top meets the head's sphere across an edge smoother
/// than a crease; the flat top is no smooth region of its own, and the
/// sphere stays a sphere rather than a sweep through both.
#[test]
#[ignore = "heavy"]
fn a_socket_head_screw() {
    let report = comes_back("m5x16_bhcs_loops.step", true, true);
    assert_eq!(report.curved_faceted, 0, "{:?}", report.fallbacks);
    assert!(report.faces <= 20, "{} faces came back", report.faces);
}

/// The disc capping its chamfered end has every corner on the chamfer's
/// cone, and its facets near the rim sag no more than a chord of the cone
/// would; the cone takes none of them, and comes back curved. The drill
/// point under its hex socket meets six flat faces along arcs of one
/// circle, and closes at its apex as a point met by one rim does. It comes
/// back on the source's twenty faces.
#[test]
fn a_button_head_screw() {
    let report = comes_back("m5x16_bhcs.step", true, true);
    assert_eq!(report.curved_faceted, 0, "{:?}", report.fallbacks);
    assert!(report.faces <= 20, "{} faces came back", report.faces);
}

/// A turned part some of whose curved faces meet their neighbours all but
/// tangentially along lines the mesh draws with three vertices; each such
/// edge is the parabola through them, and the part comes back valid with
/// nothing faceted.
#[test]
#[ignore = "heavy"]
fn sliver_on_a_diagonal_of_the_grid() {
    let report = comes_back("sliver_on_a_diagonal_of_the_grid.step", true, true);
    assert_eq!(report.curved_faceted, 0, "{:?}", report.fallbacks);
}

/// A turned part whose rounds meet single facets, chords of a spline
/// face left faceted, across their curvature: the span they share stands
/// off the round by more than a chord may, and each such facet comes back
/// as a fan onto the round, which stays curved. A quarter round shorter
/// than the sag of a facet's plane across it would fold over its far rim
/// along that plane; the facet is a fan there too.
#[test]
#[ignore = "heavy"]
fn grid_point_on_a_diagonal_boundary() {
    let report = comes_back("grid_point_on_a_diagonal_boundary.step", true, true);
    assert_eq!(report.curved_faceted, 0, "{:?}", report.fallbacks);
}

/// A degenerate spline sliver meshed and converted: facets beside its
/// curved regions are built as fans whose ruled surface has no normal
/// along the seam. Each such fan's build fails, and the curved face its
/// seam lies on is faceted (named with why); the conversion comes back.
#[test]
#[ignore = "heavy"]
fn a_face_that_cannot_be_built_is_faceted() {
    let path = format!(
        "{}/../../tests/corpus/spline_face_fit_runs_away.step",
        env!("CARGO_MANIFEST_DIR")
    );
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
    assert!(
        out.report
            .fallbacks
            .iter()
            .any(|f| f.reason == FallbackReason::BuildFailed),
        "{:?}",
        out.report.fallbacks
    );
    assert_eq!(out.report.curved_faceted, out.report.fallbacks.len());
    assert!(
        ogeom::mesh::triangulate(&back, &out.shape, Deflection::default(), T).is_ok(),
        "the converted shape does not mesh"
    );
}

/// A surface moved off its region: translated by `by` across it (along a
/// plane's normal, along a curved surface's frame `x`), or with `scale`
/// its radius grown by that factor (a plane's normal turned over).
fn displaced(surface: &Canonical, by: f64, scale: Option<f64>) -> Option<Canonical> {
    use ogeom::math::{Cone, Cylinder, Frame, Sphere, Torus, Transform};
    let shift = |frame: Frame| match surface {
        Canonical::Plane(_) => Transform::translation(frame.z().vector() * by),
        _ => Transform::translation(frame.x().vector() * by),
    };
    Some(match (surface, scale) {
        (Canonical::Plane(p), None) => Canonical::Plane(p.transformed(&shift(p.frame()), T).ok()?),
        (Canonical::Cylinder(c), None) => {
            Canonical::Cylinder(c.transformed(&shift(c.frame()), T).ok()?)
        }
        (Canonical::Cone(c), None) => Canonical::Cone(c.transformed(&shift(c.frame()), T).ok()?),
        (Canonical::Sphere(s), None) => {
            Canonical::Sphere(s.transformed(&shift(s.frame()), T).ok()?)
        }
        (Canonical::Torus(t), None) => Canonical::Torus(t.transformed(&shift(t.frame()), T).ok()?),
        (Canonical::Plane(p), Some(_)) => Canonical::Plane(p.reversed()),
        (Canonical::Cylinder(c), Some(k)) => {
            Canonical::Cylinder(Cylinder::new(c.frame(), c.radius() * k, T).ok()?)
        }
        (Canonical::Cone(c), Some(k)) => {
            Canonical::Cone(Cone::new(c.frame(), c.reference_radius() * k, c.half_angle(), T).ok()?)
        }
        (Canonical::Sphere(s), Some(k)) => {
            Canonical::Sphere(Sphere::new(s.frame(), s.radius() * k, T).ok()?)
        }
        (Canonical::Torus(t), Some(k)) => {
            Canonical::Torus(Torus::new(t.frame(), t.major_radius(), t.minor_radius() * k, T).ok()?)
        }
        _ => return None,
    })
}

/// The corpus part `name`'s own regions, one at a time put on a surface
/// moved off them (see [`displaced`]): each build still comes back a
/// valid, closed solid, and a region moved past what its seams can take up
/// is faceted, named in the report (a plane's triangles gathered again). Up to
/// `each` regions of each kind are moved, the largest first. Returns how
/// many builds were made.
fn displaced_regions_fall_back(name: &str, each: usize) -> usize {
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
    let regions = MeshRegions::find(&mesh, &MeshSolidOptions::default(), T).unwrap();
    let distance = regions.distance();
    let mut by_kind: std::collections::BTreeMap<u8, Vec<MeshRegion>> =
        std::collections::BTreeMap::new();
    for region in regions.regions() {
        let kind = match &region.surface {
            Some(Canonical::Plane(_)) => 0,
            Some(Canonical::Cylinder(_)) => 1,
            Some(Canonical::Cone(_)) => 2,
            Some(Canonical::Sphere(_)) => 3,
            Some(Canonical::Torus(_)) => 4,
            _ => continue,
        };
        by_kind.entry(kind).or_default().push(region);
    }
    let mut builds = 0;
    for (kind, mut list) in by_kind {
        list.sort_by_key(|a| std::cmp::Reverse(a.triangles.len()));
        for region in list.into_iter().take(each) {
            let surface = region.surface.clone().unwrap();
            let (mut lo, mut hi) = (
                ogeom::math::Point::new(f64::MAX, f64::MAX, f64::MAX),
                ogeom::math::Point::new(f64::MIN, f64::MIN, f64::MIN),
            );
            for &t in &region.triangles {
                for v in regions.triangles()[t] {
                    let p = regions.points()[v as usize];
                    lo = ogeom::math::Point::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
                    hi = ogeom::math::Point::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
                }
            }
            // Far off, past what any seam's tolerance takes up, or a few
            // distances off, which the build may take up or fall back on.
            let far = (distance * 50.0).max(lo.distance(hi) * 0.02);
            let near = distance * 4.0;
            let mut moves = vec![(far, None, true), (near, None, false)];
            if kind == 0 {
                moves.push((0.0, Some(-1.0), true));
            } else {
                moves.extend([(0.0, Some(1.1), true), (0.0, Some(0.5), true)]);
            }
            for (by, scale, beyond) in moves {
                let Some(moved) = displaced(&surface, by, scale) else {
                    continue;
                };
                let mut edited = regions.clone();
                edited.put_surface_unverified(region.id, moved).unwrap();
                let what = format!(
                    "{name}: region {} of kind {kind} moved {by:e}, scaled {scale:?}",
                    region.id.index()
                );
                let mut back = Model::new();
                let out = match edited.build(&mut back) {
                    Ok(out) => out,
                    Err(e) => panic!("{what}: the build failed: {e}"),
                };
                builds += 1;
                assert!(out.closed, "{what}: {:?}", out.report);
                let diagnosis = check(&back, &out.shape, T).unwrap();
                assert!(diagnosis.is_valid(), "{what}: {diagnosis}");
                if beyond {
                    assert!(
                        out.report.fallbacks.iter().any(|f| f.region == region.id),
                        "{what}: not faceted: {:?}",
                        out.report.fallbacks
                    );
                }
            }
        }
    }
    builds
}

#[test]
fn a_displaced_region_falls_back_to_facets() {
    assert!(displaced_regions_fall_back("nist_ftc_11_asme1_rb.stp", 1) > 0);
}

#[test]
#[ignore = "heavy"]
fn displaced_regions_fall_back_to_facets_across_parts() {
    for name in [
        "nist_ctc_01_asme1_rd.stp",
        "nist_ftc_06_asme1_rd.stp",
        "nist_ftc_08_asme1_rc.stp",
        "nist_ftc_09_asme1_rd.stp",
        "m5x16_bhcs.step",
    ] {
        assert!(displaced_regions_fall_back(name, 4) > 0, "{name}");
    }
}
