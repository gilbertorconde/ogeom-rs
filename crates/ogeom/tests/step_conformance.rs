//! The step-conformance corpus (CC0, `tests/corpus/step-conformance/`) read
//! and meshed at four chords. Every part reads, every face meshes, and the
//! solids mesh closed. A generated part comes with its exact shape, and its
//! mesh stands within twice the chord of it both ways: no point of the mesh
//! farther off the shape, and no point of the shape farther from the mesh,
//! which finds a missing patch or a face dropped. A shop part's volume is
//! held to what it measures.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

#[path = "support/json.rs"]
mod json;
#[path = "support/truth.rs"]
mod truth;

use ogeom::core::{FastMap, Tolerances};
use ogeom::math::Point;
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, Triangulation};
use truth::Truth;

const T: Tolerances = Tolerances::millimetres();
const CHORDS: [f64; 4] = [0.0125, 0.05, 0.1, 0.5];

fn corpus(name: &str) -> String {
    format!(
        "{}/../../tests/corpus/step-conformance/{name}",
        env!("CARGO_MANIFEST_DIR")
    )
}

/// The part's bodies where its product tree places them, then any body the
/// tree does not reach, as the file states it; with whether each is a solid.
fn bodies(import: &ogeom::io::StepImport) -> Vec<(Shape, bool)> {
    let model = import.document.model();
    let solid = |s: &Shape| {
        model
            .kind_of(s)
            .is_ok_and(|k| k == ogeom::topo::ShapeType::Solid)
    };
    let mut out: Vec<(Shape, bool)> = Vec::new();
    for root in import.document.roots() {
        for occurrence in import.document.occurrences_of(root).unwrap() {
            let is_solid = solid(&occurrence.shape);
            out.push((occurrence.shape, is_solid));
        }
    }
    let reached: Vec<_> = out.iter().map(|(s, _)| s.node()).collect();
    for body in import.solids.iter().chain(&import.shells) {
        if !reached.contains(&body.node()) {
            out.push((body.clone(), solid(body)));
        }
    }
    out
}

/// A part meshed at `chord`, every face of it, one welded mesh per body.
fn meshed(model: &Model, bodies: &[(Shape, bool)], chord: f64, name: &str) -> Vec<Triangulation> {
    let deflection = Deflection {
        chord,
        ..Deflection::default()
    };
    bodies
        .iter()
        .map(|(body, _)| {
            let (faces, _) =
                ogeom::mesh::face_meshes_for(model, body, deflection, None, T).unwrap();
            for (face, mesh) in &faces {
                let mesh = mesh
                    .as_ref()
                    .unwrap_or_else(|e| panic!("{name} at {chord}: a face refused: {e}"));
                assert!(
                    !mesh.triangles.is_empty(),
                    "{name} at {chord}: face {} meshed empty",
                    face.node().index()
                );
            }
            ogeom::mesh::triangulate(model, body, deflection, T).unwrap()
        })
        .collect()
}

/// A mesh's triangles bucketed in a grid, for the distance to the nearest.
struct Near {
    cell: f64,
    triangles: Vec<[Point; 3]>,
    grid: FastMap<(i64, i64, i64), Vec<usize>>,
}

impl Near {
    fn new(meshes: &[Triangulation]) -> Self {
        let mut triangles: Vec<[Point; 3]> = Vec::new();
        for m in meshes {
            triangles.extend(
                m.triangles
                    .iter()
                    .map(|t| t.map(|i| m.positions[i as usize])),
            );
        }
        // Cells about a triangle across, so each holds a few.
        let extent: f64 = triangles
            .iter()
            .map(|[a, b, c]| a.distance(*b).max(b.distance(*c)).max(c.distance(*a)))
            .sum();
        #[allow(clippy::cast_precision_loss)]
        let cell = (extent / triangles.len().max(1) as f64).max(1e-3);
        let key = |p: Point| {
            #[allow(clippy::cast_possible_truncation)]
            let k = |x: f64| (x / cell).floor() as i64;
            (k(p.x), k(p.y), k(p.z))
        };
        let mut grid: FastMap<(i64, i64, i64), Vec<usize>> = FastMap::default();
        for (i, [a, b, c]) in triangles.iter().enumerate() {
            let lo = key(Point::new(
                a.x.min(b.x).min(c.x),
                a.y.min(b.y).min(c.y),
                a.z.min(b.z).min(c.z),
            ));
            let hi = key(Point::new(
                a.x.max(b.x).max(c.x),
                a.y.max(b.y).max(c.y),
                a.z.max(b.z).max(c.z),
            ));
            for x in lo.0..=hi.0 {
                for y in lo.1..=hi.1 {
                    for z in lo.2..=hi.2 {
                        grid.entry((x, y, z)).or_default().push(i);
                    }
                }
            }
        }
        Self {
            cell,
            triangles,
            grid,
        }
    }

    /// How far `p` is from the nearest triangle, searched in rings of cells.
    fn distance(&self, p: Point) -> f64 {
        #[allow(clippy::cast_possible_truncation)]
        let k = |x: f64| (x / self.cell).floor() as i64;
        let (cx, cy, cz) = (k(p.x), k(p.y), k(p.z));
        let mut best = f64::INFINITY;
        for r in 0..4096_i64 {
            for dx in -r..=r {
                for dy in -r..=r {
                    for dz in -r..=r {
                        if dx.abs().max(dy.abs()).max(dz.abs()) != r {
                            continue;
                        }
                        if let Some(list) = self.grid.get(&(cx + dx, cy + dy, cz + dz)) {
                            for &i in list {
                                best = best.min(to_triangle(p, self.triangles[i]));
                            }
                        }
                    }
                }
            }
            #[allow(clippy::cast_precision_loss)]
            if best <= r as f64 * self.cell {
                break;
            }
        }
        best
    }
}

/// The distance from `p` to the triangle `[a, b, c]`.
fn to_triangle(p: Point, [a, b, c]: [Point; 3]) -> f64 {
    let (ab, ac, ap) = (b - a, c - a, p - a);
    let (d1, d2) = (ab.dot(ap), ac.dot(ap));
    if d1 <= 0.0 && d2 <= 0.0 {
        return p.distance(a);
    }
    let bp = p - b;
    let (d3, d4) = (ab.dot(bp), ac.dot(bp));
    if d3 >= 0.0 && d4 <= d3 {
        return p.distance(b);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return p.distance(a + ab * (d1 / (d1 - d3)));
    }
    let cp = p - c;
    let (d5, d6) = (ab.dot(cp), ac.dot(cp));
    if d6 >= 0.0 && d5 <= d6 {
        return p.distance(c);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return p.distance(a + ac * (d2 / (d2 - d6)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return p.distance(b + (c - b) * w);
    }
    let scale = 1.0 / (va + vb + vc);
    p.distance(a + ab * (vb * scale) + ac * (vc * scale))
}

/// Every generated part meshes within twice the chord of its exact shape,
/// both ways, closed.
#[test]
fn generated_parts_mesh_within_twice_the_chord_of_their_shape() {
    let dir = corpus("generated");
    let mut names: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| {
            let name = e.ok()?.file_name().into_string().ok()?;
            name.strip_suffix(".step").map(str::to_string)
        })
        .collect();
    names.sort();
    assert!(names.len() >= 19, "the generated parts are vendored");
    for name in &names {
        let text = std::fs::read_to_string(format!("{dir}/{name}.step")).unwrap();
        let truth_text = std::fs::read_to_string(format!("{dir}/{name}.truth.json")).unwrap();
        let truth = Truth::read(&truth_text).unwrap_or_else(|| panic!("{name}: truth unread"));
        let import = ogeom::io::read_step(&text, T).unwrap_or_else(|e| panic!("{name}: {e}"));
        let model = import.document.model();
        let parts = bodies(&import);
        assert!(!parts.is_empty(), "{name}: no body");
        let samples = truth.samples(64, 128);
        for chord in CHORDS {
            let meshes = meshed(model, &parts, chord, name);
            for ((_, is_solid), mesh) in parts.iter().zip(&meshes) {
                if *is_solid {
                    assert!(mesh.is_closed(), "{name} at {chord}: a solid meshed open");
                }
            }
            let mut out = 0.0_f64;
            for m in &meshes {
                for p in &m.positions {
                    out = out.max(truth.distance(*p));
                }
                for t in &m.triangles {
                    let [a, b, c] = t.map(|i| m.positions[i as usize]);
                    let middle =
                        Point::from_vector((a.to_vector() + b.to_vector() + c.to_vector()) / 3.0);
                    out = out.max(truth.distance(middle));
                }
            }
            let near = Near::new(&meshes);
            let missing = samples
                .iter()
                .map(|p| near.distance(*p))
                .fold(0.0_f64, f64::max);
            assert!(
                out <= 2.0 * chord && missing <= 2.0 * chord,
                "{name} at {chord}: the mesh stands {out:.4} off the shape and {missing:.4} short of it"
            );
        }
    }
}

/// The sides of `mesh` used unevenly by its triangles: where it is open.
fn open_sides(mesh: &Triangulation) -> Vec<(Point, Point)> {
    let mut balance: FastMap<(u32, u32), i64> = FastMap::default();
    for t in &mesh.triangles {
        for i in 0..3 {
            let (a, b) = (t[i], t[(i + 1) % 3]);
            let (key, step) = if a <= b { ((a, b), 1) } else { ((b, a), -1) };
            *balance.entry(key).or_default() += step;
        }
    }
    balance
        .iter()
        .filter(|(_, n)| **n != 0)
        .map(|((a, b), _)| (mesh.positions[*a as usize], mesh.positions[*b as usize]))
        .collect()
}

/// Every shop part reads, meshes every face, and its solids mesh closed;
/// each measures the volume pinned for it. `freeform_panel`'s bottom rim,
/// curve `#739` of the file, crosses itself in a loop a millimetre and a
/// half across near (-120.4, -100.4, 0): its bottom face bounds no region
/// there, and the mesh is open around the loop and nowhere else.
#[test]
fn shop_parts_mesh_closed_and_measure_what_they_did() {
    let rim_loop = Point::new(-120.4, -100.4, 0.0);
    let pinned: [(&str, Option<f64>, Option<Point>); 8] = [
        ("bowl_01", None, None),
        ("door_handle_guide", Some(167_128.183), None),
        ("freeform_panel", Some(3_016_852.159), Some(rim_loop)),
        ("hinge_jig", Some(8_595.660), None),
        ("lego_3_pins", Some(407_740.906), None),
        ("lego_4x2_01", Some(135_511.436), None),
        ("rail_guide", Some(88_218.257), None),
        ("sled_base", Some(771_831.885), None),
    ];
    for (name, volume, open_at) in pinned {
        let text = std::fs::read_to_string(corpus(&format!("{name}.step"))).unwrap();
        let import = ogeom::io::read_step(&text, T).unwrap_or_else(|e| panic!("{name}: {e}"));
        let model = import.document.model();
        let parts = bodies(&import);
        assert!(!parts.is_empty(), "{name}: no body");
        for chord in CHORDS {
            let meshes = meshed(model, &parts, chord, name);
            for ((_, is_solid), mesh) in parts.iter().zip(&meshes) {
                if !*is_solid {
                    continue;
                }
                let open = open_sides(mesh);
                let stray = open.iter().find(|(a, b)| {
                    open_at.is_none_or(|at| a.distance(at).max(b.distance(at)) > 3.0)
                });
                assert!(
                    stray.is_none(),
                    "{name} at {chord}: a solid meshed open at {stray:?}"
                );
            }
        }
        let measured: f64 = import
            .solids
            .iter()
            .map(|s| {
                ogeom::algo::volume_properties(model, s, Deflection::default(), T)
                    .unwrap()
                    .mass
            })
            .sum();
        if let Some(volume) = volume {
            assert!(
                (measured - volume).abs() <= volume.abs() * 1e-3,
                "{name}: volume {measured} against its pinned {volume}"
            );
        }
    }
}

/// A face curved one way only (a drum, a cone) is meshed in strips along
/// its rulings, held apart by its columns, with no rows it does not need:
/// at the finest chord a cone and a shop part full of drums mesh in at
/// most a quarter more triangles than the conformance suite's reference
/// reader draws for them.
#[test]
fn ruled_faces_mesh_in_strips_along_their_rulings() {
    let chord = CHORDS[0];
    for (name, at_most) in [("generated/cone_point", 1_954), ("sled_base", 42_278)] {
        let text = std::fs::read_to_string(corpus(&format!("{name}.step"))).unwrap();
        let import = ogeom::io::read_step(&text, T).unwrap();
        let parts = bodies(&import);
        let triangles: usize = meshed(import.document.model(), &parts, chord, name)
            .iter()
            .map(|m| m.triangles.len())
            .sum();
        assert!(
            triangles <= at_most,
            "{name} at {chord}: {triangles} triangles, more than {at_most}"
        );
    }
}

/// A wall curved both ways twists between the rows and columns its grid
/// holds to the chord, and a triangle across a cell can stand off the wall
/// in its middle though its corners and edges keep to it. Every triangle
/// of freeform_panel's wall stands within the chord of the wall at its
/// middle, lifted onto the surface.
#[test]
fn a_twisting_wall_meshes_within_its_chord_in_every_triangle() {
    use ogeom::geom::Surface as _;
    let text = std::fs::read_to_string(corpus("freeform_panel.step")).unwrap();
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let body = &import.solids[0];
    let wall =
        ogeom::topo::explore_unique(model, body, ogeom::topo::ShapeType::Face).unwrap()[0].clone();
    let data = model.node(&wall).unwrap().data().as_face().unwrap();
    let surface = model.geometry().surface(data.surface).unwrap();
    for chord in [0.1, 0.0125] {
        let deflection = Deflection {
            chord,
            ..Deflection::default()
        };
        let (faces, _) = ogeom::mesh::face_meshes_for(model, body, deflection, None, T).unwrap();
        let mesh = faces
            .iter()
            .find(|(f, _)| f.node() == wall.node())
            .unwrap()
            .1
            .as_ref()
            .unwrap();
        for (i, t) in mesh.triangles.iter().enumerate() {
            let [a, b, c] = t.map(|k| mesh.positions[k as usize]);
            let uv = mesh.triangle_parameters(i).unwrap();
            let Ok(normal) = (b - a).cross(c - a).normalized(T) else {
                continue;
            };
            let middle = (
                (uv[0].0 + uv[1].0 + uv[2].0) / 3.0,
                (uv[0].1 + uv[1].1 + uv[2].1) / 3.0,
            );
            let on = surface.point_at(middle.0, middle.1, T).unwrap();
            let off = (on - a).dot(normal).abs();
            assert!(
                off <= chord * (1.0 + 1e-9),
                "chord {chord}: a middle {off} off"
            );
        }
    }
}

/// Every ruling of a cone meets at its apex, and a triangle there joining
/// points on rulings far apart cuts inside the cone though its corners and
/// sides lie on it. Every triangle of cone_point's curved face, at its
/// middle and at the middles of its sides lifted onto the cone, stands
/// within the chord of its own plane.
#[test]
fn a_cone_to_its_apex_meshes_within_its_chord_in_every_triangle() {
    use ogeom::geom::Surface as _;
    let chord = 0.0125;
    let text = std::fs::read_to_string(corpus("generated/cone_point.step")).unwrap();
    let import = ogeom::io::read_step(&text, T).unwrap();
    let model = import.document.model();
    let deflection = Deflection {
        chord,
        ..Deflection::default()
    };
    let mut cones = 0;
    for body in import.solids.iter().chain(&import.shells) {
        let (faces, _) = ogeom::mesh::face_meshes_for(model, body, deflection, None, T).unwrap();
        for (face, mesh) in &faces {
            let data = model.node(face).unwrap().data().as_face().unwrap();
            let surface = model.geometry().surface(data.surface).unwrap();
            if surface.kind() != ogeom::geom::SurfaceKind::Cone {
                continue;
            }
            cones += 1;
            let mesh = mesh.as_ref().unwrap();
            for (i, t) in mesh.triangles.iter().enumerate() {
                let [a, b, c] = t.map(|k| mesh.positions[k as usize]);
                let uv = mesh.triangle_parameters(i).unwrap();
                let Ok(normal) = (b - a).cross(c - a).normalized(T) else {
                    continue;
                };
                let third = 1.0 / 3.0;
                for w in [
                    [third, third, third],
                    [0.5, 0.5, 0.0],
                    [0.0, 0.5, 0.5],
                    [0.5, 0.0, 0.5],
                ] {
                    let u = uv[0].0 * w[0] + uv[1].0 * w[1] + uv[2].0 * w[2];
                    let v = uv[0].1 * w[0] + uv[1].1 * w[1] + uv[2].1 * w[2];
                    let on = surface.point_at(u, v, T).unwrap();
                    let off = (on - a).dot(normal).abs();
                    assert!(off <= chord * (1.0 + 1e-9), "{off} off at ({u}, {v})");
                }
            }
        }
    }
    assert_eq!(cones, 1);
}

/// A drum or a cone is drawn between its two ends, its rulings held as
/// constraints, with no row through its middle that asks nothing of the
/// chord. door_handle_guide and rail_guide, mostly drums and cones, mesh at
/// 0.0125 in at most 1.25 times the reference reader's triangles (1,948 and
/// 6,484), and every triangle of their curved faces, at its middle and the middles
/// of its sides, stands within the chord of its face.
#[test]
fn drums_and_cones_mesh_between_their_rims_within_their_chord() {
    use ogeom::geom::Surface as _;
    let chord = 0.0125;
    for (part, ceiling) in [
        ("door_handle_guide.step", 2_435),
        ("rail_guide.step", 8_105),
    ] {
        let text = std::fs::read_to_string(corpus(part)).unwrap();
        let import = ogeom::io::read_step(&text, T).unwrap();
        let model = import.document.model();
        let mut triangles = 0;
        for body in import.solids.iter().chain(&import.shells) {
            let deflection = Deflection {
                chord,
                ..Deflection::default()
            };
            let (faces, _) =
                ogeom::mesh::face_meshes_for(model, body, deflection, None, T).unwrap();
            for (face, mesh) in &faces {
                let mesh = mesh.as_ref().unwrap();
                triangles += mesh.triangles.len();
                let data = model.node(face).unwrap().data().as_face().unwrap();
                let surface = model.geometry().surface(data.surface).unwrap();
                for (i, t) in mesh.triangles.iter().enumerate() {
                    let [a, b, c] = t.map(|k| mesh.positions[k as usize]);
                    let uv = mesh.triangle_parameters(i).unwrap();
                    let Ok(normal) = (b - a).cross(c - a).normalized(T) else {
                        continue;
                    };
                    let third = 1.0 / 3.0;
                    for w in [
                        [third, third, third],
                        [0.5, 0.5, 0.0],
                        [0.0, 0.5, 0.5],
                        [0.5, 0.0, 0.5],
                    ] {
                        let u = uv[0].0 * w[0] + uv[1].0 * w[1] + uv[2].0 * w[2];
                        let v = uv[0].1 * w[0] + uv[1].1 * w[1] + uv[2].1 * w[2];
                        let on = surface.point_at(u, v, T).unwrap();
                        let off = (on - a).dot(normal).abs();
                        assert!(off <= chord * (1.0 + 1e-9), "{part}: {off} off");
                    }
                }
            }
        }
        assert!(triangles <= ceiling, "{part}: {triangles} triangles");
    }
}

/// A face that runs to a point where its surface collapses (a cone's apex,
/// a sphere's pole) meshes that point as one vertex of the face, a fan of
/// real triangles round it. At three chords on the parts that end in an
/// apex or a pole: no triangle has two corners at one point or no area,
/// each degenerate point is one vertex per face, every face meets its
/// chord, each body welds closed, and the mean of every triangle's corner
/// parameters lifts to a point within the chord of the triangle's plane.
#[test]
fn an_apex_or_a_pole_meshes_as_one_vertex_and_no_triangle_of_no_area() {
    use ogeom::geom::Surface as _;
    for name in [
        "generated/cone_point.step",
        "generated/cone_point_vertex.step",
        "generated/dome.step",
    ] {
        let text = std::fs::read_to_string(corpus(name)).unwrap();
        let import = ogeom::io::read_step(&text, T).unwrap();
        let model = import.document.model();
        let mut apexes = 0;
        for chord in [0.01 / 3.0, 0.1, 1.0] {
            let deflection = Deflection {
                chord,
                ..Deflection::default()
            };
            for body in import.solids.iter().chain(&import.shells) {
                let (faces, _) =
                    ogeom::mesh::face_meshes_for(model, body, deflection, None, T).unwrap();
                for (n, (face, mesh)) in faces.iter().enumerate() {
                    let mesh = mesh.as_ref().unwrap();
                    let at = format!("{name} at chord {chord}, face {n}");
                    assert!(mesh.deflection_met, "{at}: misses its chord");
                    let data = model.node(face).unwrap().data().as_face().unwrap();
                    let surface = model.geometry().surface(data.surface).unwrap();
                    let mut degenerate: Vec<(Point, Vec<u32>)> = Vec::new();
                    for (i, t) in mesh.triangles.iter().enumerate() {
                        let [a, b, c] = t.map(|k| mesh.positions[k as usize]);
                        let longest = a.distance(b).max(b.distance(c)).max(c.distance(a));
                        let area = (b - a).cross(c - a).magnitude();
                        assert!(
                            a.distance(b) > T.confusion()
                                && b.distance(c) > T.confusion()
                                && c.distance(a) > T.confusion()
                                && area > T.confusion() * longest,
                            "{at}: triangle {t:?} at {a:?}, {b:?}, {c:?} has no area"
                        );
                        let uv = mesh.triangle_parameters(i).unwrap();
                        for (k, &(u, v)) in uv.iter().enumerate() {
                            if surface.is_degenerate_at(u, v, T).unwrap() {
                                let p = mesh.positions[t[k] as usize];
                                match degenerate
                                    .iter_mut()
                                    .find(|(q, _)| q.distance(p) <= T.confusion())
                                {
                                    Some((_, vertices)) => vertices.push(t[k]),
                                    None => degenerate.push((p, vec![t[k]])),
                                }
                            }
                        }
                        let normal = (b - a).cross(c - a).normalized(T).unwrap();
                        let mean = (
                            (uv[0].0 + uv[1].0 + uv[2].0) / 3.0,
                            (uv[0].1 + uv[1].1 + uv[2].1) / 3.0,
                        );
                        let on = surface.point_at(mean.0, mean.1, T).unwrap();
                        let off = (on - a).dot(normal).abs();
                        assert!(
                            off <= chord + T.confusion(),
                            "{at}: the mean of {uv:?} lifts {off} off its triangle's plane"
                        );
                    }
                    for (p, mut vertices) in degenerate {
                        vertices.sort_unstable();
                        vertices.dedup();
                        assert_eq!(vertices.len(), 1, "{at}: {vertices:?} all stand at {p:?}");
                        apexes += 1;
                    }
                }
                let welded = ogeom::mesh::triangulate(model, body, deflection, T).unwrap();
                assert!(welded.is_closed(), "{name} at chord {chord}: not closed");
            }
        }
        assert!(apexes >= 3, "{name}: no face reaches an apex or a pole");
    }
}

/// Every triangle of every face that meets its chord is wound with its
/// face: its normal in space agrees in sign with the face's normal at the
/// triangle's centre, found from its corner parameters. Near a pole, where
/// one chart point stands for a whole row, is where a fold would show. A
/// face that reports its chord missed is not held to it: freeform_panel's
/// wall folds over itself, its normal turning through nothing inside the
/// face, and no flat triangle across the fold can agree with both sides.
#[test]
#[ignore = "heavy: every conformance part at three chords, the finest 0.0033"]
fn every_triangle_is_wound_with_its_face() {
    use ogeom::geom::Surface as _;
    use ogeom::topo::Orientation;
    let mut names: Vec<String> = Vec::new();
    for dir in ["", "generated/"] {
        for e in std::fs::read_dir(corpus(dir)).unwrap() {
            let name = e.unwrap().file_name().into_string().unwrap();
            if name.ends_with(".step") {
                names.push(format!("{dir}{name}"));
            }
        }
    }
    names.sort();
    assert!(names.len() >= 27, "the conformance parts are vendored");
    let mut against = Vec::new();
    let mut poles = 0;
    for name in &names {
        let text = std::fs::read_to_string(corpus(name)).unwrap();
        let import = ogeom::io::read_step(&text, T).unwrap();
        let model = import.document.model();
        for chord in [0.01 / 3.0, 0.1, 1.0] {
            let deflection = Deflection {
                chord,
                ..Deflection::default()
            };
            for body in import.solids.iter().chain(&import.shells) {
                let (faces, _) =
                    ogeom::mesh::face_meshes_for(model, body, deflection, None, T).unwrap();
                for (n, (face, mesh)) in faces.iter().enumerate() {
                    let mesh = mesh.as_ref().unwrap();
                    if !mesh.deflection_met {
                        continue;
                    }
                    let data = model.node(face).unwrap().data().as_face().unwrap();
                    let surface = model.geometry().surface(data.surface).unwrap();
                    if mesh
                        .corner_parameters
                        .iter()
                        .flatten()
                        .any(|&(u, v)| surface.is_degenerate_at(u, v, T).unwrap_or(false))
                    {
                        poles += 1;
                    }
                    let placement = face.transform(model.datums()).unwrap();
                    let side = if face.orientation() == Orientation::Reversed {
                        -1.0
                    } else {
                        1.0
                    };
                    for (i, t) in mesh.triangles.iter().enumerate() {
                        let [a, b, c] = t.map(|k| mesh.positions[k as usize]);
                        let uv = mesh.triangle_parameters(i).unwrap();
                        let (u, v) = (
                            (uv[0].0 + uv[1].0 + uv[2].0) / 3.0,
                            (uv[0].1 + uv[1].1 + uv[2].1) / 3.0,
                        );
                        let Ok(normal) = surface.normal_at(u, v, T) else {
                            continue;
                        };
                        let normal = placement.apply_vector(normal.vector()) * side;
                        if (b - a).cross(c - a).dot(normal) <= 0.0 {
                            against.push(format!("{name} at chord {chord}, face {n}: {uv:?}"));
                        }
                    }
                }
            }
        }
    }
    assert!(poles >= 12, "only {poles} faces reach a pole or an apex");
    assert!(
        against.is_empty(),
        "{} triangles wound against their face:\n{}",
        against.len(),
        against.join("\n")
    );
}

/// The sides of `mesh` not used exactly once each way: where it is open,
/// or where more than two triangles meet.
fn not_manifold(mesh: &Triangulation) -> Vec<(Point, Point, usize, usize)> {
    let mut uses: FastMap<(u32, u32), (usize, usize)> = FastMap::default();
    for t in &mesh.triangles {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            let e = uses.entry((a.min(b), a.max(b))).or_default();
            if a < b {
                e.0 += 1;
            } else {
                e.1 += 1;
            }
        }
    }
    uses.into_iter()
        .filter(|(_, n)| *n != (1, 1))
        .map(|((a, b), (f, r))| (mesh.positions[a as usize], mesh.positions[b as usize], f, r))
        .collect()
}

/// Each face's mesh meets its neighbours' along their edges and nowhere
/// else. No curved face draws a triangle across three points of one edge:
/// its corners are the edge's alone, the face across can draw the same
/// triangle, and welded, the chord under the two is four triangles' side.
/// And the faces of a valid solid, welded, use every side once each way.
/// freeform_panel's is not: its wall folds through itself, and its mesh
/// cannot close where the face does not.
#[test]
#[ignore = "heavy: every conformance part at four chords"]
fn every_solid_meshes_two_manifold() {
    let mut names: Vec<String> = Vec::new();
    for dir in ["", "generated/"] {
        for e in std::fs::read_dir(corpus(dir)).unwrap() {
            let name = e.unwrap().file_name().into_string().unwrap();
            if name.ends_with(".step") {
                names.push(format!("{dir}{name}"));
            }
        }
    }
    names.sort();
    let mut bad = Vec::new();
    let mut ears = Vec::new();
    let mut solids = 0;
    for name in &names {
        let text = std::fs::read_to_string(corpus(name)).unwrap();
        let import = ogeom::io::read_step(&text, T).unwrap();
        let model = import.document.model();
        for chord in [0.0125, 0.05, 0.1, 0.25, 0.5] {
            let deflection = Deflection {
                chord,
                ..Deflection::default()
            };
            for body in &import.solids {
                let (faces, _) =
                    ogeom::mesh::face_meshes_for(model, body, deflection, None, T).unwrap();
                let mut joined = Triangulation::new();
                for (face, mesh) in &faces {
                    let mesh = mesh.as_ref().unwrap();
                    joined.append(mesh);
                    for (i, t) in edge_ears(model, face, mesh) {
                        ears.push(format!("{name} at {chord}: triangle {i} {t:?}"));
                    }
                }
                if !ogeom::algo::check(model, body, T).unwrap().is_valid() {
                    continue;
                }
                solids += 1;
                for (a, b, f, r) in not_manifold(&joined.welded(T)) {
                    bad.push(format!("{name} at {chord}: {a:?}-{b:?} used {f}+{r}"));
                }
            }
        }
    }
    assert!(solids >= 100, "only {solids} solids welded");
    assert!(
        ears.is_empty(),
        "{} triangles across one edge:\n{}",
        ears.len(),
        ears.join("\n")
    );
    assert!(bad.is_empty(), "{} sides:\n{}", bad.len(), bad.join("\n"));
}

/// The triangles of a curved face's mesh drawn across three points of one
/// edge: two sides on the mesh's border, meeting at a point that is no
/// vertex of the face.
fn edge_ears(model: &Model, face: &Shape, mesh: &Triangulation) -> Vec<(usize, [u32; 3])> {
    use ogeom::geom::Surface as _;
    use ogeom::topo::ShapeType;
    let data = model.node(face).unwrap().data().as_face().unwrap();
    let surface = model.geometry().surface(data.surface).unwrap();
    if surface.kind() == ogeom::geom::SurfaceKind::Plane {
        return Vec::new();
    }
    let corners: Vec<(Point, f64)> =
        ogeom::topo::explore(model, face, ogeom::topo::Filter::OfType(ShapeType::Vertex))
            .unwrap()
            .iter()
            .map(|v| {
                let d = model.node(v).unwrap().data().as_vertex().unwrap();
                let at = v.transform(model.datums()).unwrap().apply(d.point);
                (at, d.tolerance.get().max(T.confusion()))
            })
            .collect();
    let mut uses: FastMap<(u32, u32), usize> = FastMap::default();
    for t in &mesh.triangles {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            *uses.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    let border = |a: u32, b: u32| uses.get(&(a.min(b), a.max(b))) == Some(&1);
    mesh.triangles
        .iter()
        .enumerate()
        .filter(|(_, t)| {
            (0..3).any(|k| {
                let (a, m, b) = (t[k], t[(k + 1) % 3], t[(k + 2) % 3]);
                let at = mesh.positions[m as usize];
                border(a, m)
                    && border(m, b)
                    && !border(b, a)
                    && !corners.iter().any(|&(p, within)| p.distance(at) <= within)
            })
        })
        .map(|(i, t)| (i, *t))
        .collect()
}
