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
