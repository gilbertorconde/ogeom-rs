//! Randomized robustness runs over the kernel: many generated cases, each
//! judged by an oracle, the tally watched against a recorded baseline.
//!
//! A single repro proves one case. A kernel is trusted on the rate at which
//! the cases nobody wrote down go wrong. Each scenario draws its cases from
//! a seeded generator, biased toward the degenerate placements that break
//! booleans and blends (a drill tangent to a face, through a vertex, a pair
//! of solids with coincident faces), and every case ends in one of:
//!
//! - `ok`: the result is a valid solid and the oracle holds;
//! - `refused`: the kernel said why it would not (an error by name);
//! - `invalid`: it answered with a solid that fails `check`;
//! - `wrong`: a valid answer the oracle contradicts;
//! - `panic`: it crashed.
//!
//! Refusals are honest and tracked. The last three are the bugs.
//!
//! The oracles are volume identities (a cut and a common add up to the
//! part), read on volumes integrated on the exact surfaces and held to a
//! millionth of the solids. A solid the exact integral refuses is measured
//! on a fine tessellation instead, allowed the chord times its area on
//! top, and the report lists those cases.
//!
//! `ogeom-stress` runs every scenario and prints the tally. `--check
//! <baseline.json>` fails when a scenario has fewer `ok` or more bad cases
//! than the baseline records, and `--write <path>` records one. `--case
//! drill/box/17` replays one case and prints what happened.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "a reporting tool"
)]
#![allow(
    clippy::disallowed_methods,
    reason = "a native harness: it times with the std clock"
)]

use std::collections::BTreeMap;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

// --- the generator ---------------------------------------------------------

/// SplitMix64: small, fast, and the same everywhere.
#[derive(Clone)]
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn chance(&mut self, p: f64) -> bool {
        self.unit() < p
    }

    fn pick<'a, T>(&mut self, items: &'a [T]) -> &'a T {
        &items[self.below(items.len())]
    }
}

/// A case's own seed, from the run's seed and the case's key, so a case
/// replays alone exactly as it ran among the others.
fn case_seed(seed: u64, key: &str) -> u64 {
    let mut h: u64 = 0xCBF2_9CE4_8422_2325 ^ seed;
    for b in key.bytes() {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0100_0000_01B3);
    }
    h
}

// --- outcomes --------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Outcome {
    Ok,
    Refused,
    Invalid,
    Wrong,
    Panic,
}

impl Outcome {
    const ALL: [Self; 5] = [
        Self::Ok,
        Self::Refused,
        Self::Invalid,
        Self::Wrong,
        Self::Panic,
    ];

    const fn name(self) -> &'static str {
        match self {
            Self::Ok => "ok",
            Self::Refused => "refused",
            Self::Invalid => "invalid",
            Self::Wrong => "wrong",
            Self::Panic => "panic",
        }
    }

    const fn is_bad(self) -> bool {
        matches!(self, Self::Invalid | Self::Wrong | Self::Panic)
    }
}

struct Verdict {
    outcome: Outcome,
    note: String,
    /// Whether a volume the oracle read was measured on a tessellation,
    /// the exact integral having refused the solid.
    tessellated: bool,
}

impl Verdict {
    fn ok() -> Self {
        Self::of(Outcome::Ok, "")
    }

    fn of(outcome: Outcome, note: impl Into<String>) -> Self {
        Self {
            outcome,
            note: note.into(),
            tessellated: false,
        }
    }

    /// This verdict, marked as read from these measures.
    fn read_from(mut self, named: &[(&str, &Measure)]) -> Self {
        let fell = tessellated(named);
        if !fell.is_empty() {
            self.tessellated = true;
            if self.outcome != Outcome::Ok {
                self.note.push_str(&fell);
            }
        }
        self
    }
}

/// An operation's result judged for validity: an error is a refusal, an
/// invalid solid is `invalid`, and a valid one comes back with its measure.
fn judged(
    model: &Model,
    what: &str,
    result: ogeom::core::OgeomResult<ogeom::algo::Built>,
) -> Result<(Shape, Measure), Verdict> {
    let built = result.map_err(|e| Verdict::of(Outcome::Refused, format!("{what}: {e}")))?;
    let diagnosis = ogeom::algo::check(model, &built.shape, T)
        .map_err(|e| Verdict::of(Outcome::Refused, format!("{what}: check: {e}")))?;
    if !diagnosis.is_valid() {
        let first = diagnosis
            .problems
            .first()
            .map_or_else(String::new, |p| p.what.clone());
        return Err(Verdict::of(
            Outcome::Invalid,
            format!("{what}: invalid: {first}"),
        ));
    }
    let m = measure(model, &built.shape)
        .map_err(|e| Verdict::of(Outcome::Refused, format!("{what}: volume: {e}")))?;
    Ok((built.shape, m))
}

/// The chord a volume is tessellated to where the exact integral refuses.
/// An inscribed tessellation's volume is short by under the chord times
/// the curved area, so this keeps a fallen-back identity within about
/// 1e-4 of its solids.
const CHORD: f64 = 2e-4;

/// The relative slack a volume identity is allowed on exactly integrated
/// solids. The integral rounds to about 1e-10 of a solid; booleans of
/// generated parts land within 1e-8, and of converted corpus parts (whose
/// vertices stand off their faces by up to their tolerance) within 2e-7.
/// A face split the wrong way misses by 1e-5 to 1e-3.
const SLACK: f64 = 1e-6;

/// A solid's volume, and how far it can honestly be off.
#[derive(Clone, Copy)]
struct Measure {
    volume: f64,
    /// Whether the volume was integrated on the exact surfaces rather than
    /// on a tessellation.
    exact: bool,
    /// How far a tessellated volume can be off: the chord times the
    /// solid's area. Nothing for an exact one.
    doubt: f64,
}

impl Measure {
    const NOTHING: Self = Self {
        volume: 0.0,
        exact: true,
        doubt: 0.0,
    };
}

fn measure(model: &Model, shape: &Shape) -> ogeom::core::OgeomResult<Measure> {
    let deflection = Deflection::with_chord(CHORD)?;
    let props = ogeom::algo::volume_properties(model, shape, deflection, T)?;
    // An empty answer has nothing to integrate either way.
    let exact =
        props.deflection == 0.0 || explore_unique(model, shape, ShapeType::Face)?.is_empty();
    let doubt = if exact {
        0.0
    } else {
        CHORD * ogeom::algo::surface_properties(model, shape, deflection, T)?.mass
    };
    Ok(Measure {
        volume: props.mass,
        exact,
        doubt,
    })
}

/// How far a solid's volume may move with the solid as it states itself:
/// each face's area times the loosest tolerance on it (its own, its
/// edges', its vertices'), which bounds how far the face's boundary may
/// stand off where an integral runs it. Reported beside a wrong answer on
/// a part, to tell a part that cannot be measured closer from a boolean
/// that went wrong.
fn stated_doubt(model: &Model, shape: &Shape) -> ogeom::core::OgeomResult<f64> {
    let deflection = Deflection::with_chord(CHORD)?;
    let mut doubt = 0.0;
    for face in explore_unique(model, shape, ShapeType::Face)? {
        let mut loosest = model
            .node(&face)
            .and_then(|n| n.data().as_face().map(|d| d.tolerance.get()))
            .unwrap_or(0.0);
        for edge in explore_unique(model, &face, ShapeType::Edge)? {
            if let Some(d) = model.node(&edge).and_then(|n| n.data().as_edge()) {
                loosest = loosest.max(d.tolerance.get());
            }
        }
        for vertex in explore_unique(model, &face, ShapeType::Vertex)? {
            if let Some(d) = model.node(&vertex).and_then(|n| n.data().as_vertex()) {
                loosest = loosest.max(d.tolerance.get());
            }
        }
        let area = ogeom::algo::surface_properties(model, &face, deflection, T)?.mass;
        doubt += loosest * area;
    }
    Ok(doubt)
}

/// How far an identity over these measures may miss: the relative slack
/// on `scale`, and every tessellated measure's doubt.
fn allowance(scale: f64, measures: &[&Measure]) -> f64 {
    SLACK * scale + measures.iter().map(|m| m.doubt).sum::<f64>()
}

/// Which of these measures fell back to a tessellation, for the report.
fn tessellated(named: &[(&str, &Measure)]) -> String {
    let fell: Vec<&str> = named
        .iter()
        .filter(|(_, m)| !m.exact)
        .map(|(n, _)| *n)
        .collect();
    if fell.is_empty() {
        String::new()
    } else {
        format!("; tessellated: {}", fell.join(", "))
    }
}

// --- parts -----------------------------------------------------------------

/// A solid the scenarios work on, built once and copied per case.
#[derive(Clone)]
struct Part {
    name: String,
    model: Model,
    shape: Shape,
    measure: Measure,
    /// Its tolerances times its faces' areas: see [`stated_doubt`].
    stated: f64,
    low: Point,
    high: Point,
    vertices: Vec<Point>,
    /// Edges with a curve between two different faces: what a fillet can
    /// round.
    creases: Vec<Shape>,
}

impl Part {
    fn new(name: &str, model: Model, shape: Shape) -> Option<Self> {
        let measure = measure(&model, &shape).ok()?;
        let stated = stated_doubt(&model, &shape).ok()?;
        let bounds = ogeom::algo::shape_bounds(&model, &shape, T).ok()?;
        let (low, high) = (bounds.low()?, bounds.high()?);
        let vertices = explore_unique(&model, &shape, ShapeType::Vertex)
            .ok()?
            .iter()
            .filter_map(|v| model.node(v)?.data().as_vertex().map(|d| d.point))
            .collect();
        let creases = creases(&model, &shape);
        Some(Self {
            name: name.to_string(),
            model,
            shape,
            measure,
            stated,
            low,
            high,
            vertices,
            creases,
        })
    }

    fn diagonal(&self) -> f64 {
        self.low.distance(self.high)
    }
}

/// Edges with a curve between two different faces that meet at an angle.
/// A seam (one face on both sides), a pole, and an edge whose faces meet
/// tangent (a round's own rail, a split between coplanar faces) have no
/// corner to round.
fn creases(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap_or_default()
        .into_iter()
        .filter(|e| {
            let has_curve = model
                .node(e)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|d| d.curve3d().is_some() && !d.degenerate);
            let mut faces =
                ogeom::topo::ancestors_of(model, shape, e, ShapeType::Face).unwrap_or_default();
            faces.sort_unstable_by_key(Shape::node);
            faces.dedup_by_key(|f| f.node());
            has_curve && faces.len() == 2 && !meet_tangent(model, e, &faces)
        })
        .collect()
}

/// Whether two faces' outward normals agree at an edge's midpoint.
fn meet_tangent(model: &Model, edge: &Shape, faces: &[Shape]) -> bool {
    use ogeom::geom::{Curve3d as _, Surface as _, Transformable as _};
    let Ok(Some((curve, range))) = ogeom::algo::placed_edge_curve(model, edge, T) else {
        return false;
    };
    let Ok(mid) = curve.point_at(f64::midpoint(range.0, range.1), T) else {
        return false;
    };
    let normal = |face: &Shape| -> Option<Vector> {
        let data = model.node(face)?.data().as_face()?;
        let surface = model.geometry().surface(data.surface)?.clone();
        let placed = surface
            .transformed(&face.transform(model.datums()).ok()?, T)
            .ok()?;
        let (u, v) = ogeom::algo::project_on_surface(&placed, mid, 16, T)
            .ok()?
            .parameters;
        let n = placed.normal_at(u, v, T).ok()?.vector();
        Some(
            if face.orientation() == ogeom::topo::Orientation::Reversed {
                -n
            } else {
                n
            },
        )
    };
    let (Some(a), Some(b)) = (normal(&faces[0]), normal(&faces[1])) else {
        return false;
    };
    a.dot(b).clamp(-1.0, 1.0).acos() <= T.angular()
}

fn frame_at(origin: Point, z: Direction, x: Direction) -> Frame {
    Frame::new(origin, z, x, T).unwrap()
}

/// The generated parts: every face kind a boolean or a blend meets, and a
/// converted mesh.
fn generated_parts() -> Vec<Part> {
    let mut parts = Vec::new();
    let mut add = |name: &str, build: &dyn Fn(&mut Model) -> Option<Shape>| {
        let mut model = Model::new();
        if let Some(shape) = build(&mut model)
            && let Some(part) = Part::new(name, model, shape)
        {
            parts.push(part);
        } else {
            eprintln!("warning: part {name} did not build");
        }
    };
    let make_box = |m: &mut Model, at: Point, size: (f64, f64, f64)| {
        ogeom::algo::make_box(m, frame_at(at, Direction::Z, Direction::X), size, T)
            .ok()
            .map(|b| b.shape)
    };
    let cylinder = |m: &mut Model, at: Point, z: Direction, x: Direction, r: f64, h: f64| {
        ogeom::algo::make_cylinder(m, frame_at(at, z, x), r, h, T)
            .ok()
            .map(|b| b.shape)
    };

    add("box", &|m| make_box(m, Point::ORIGIN, (20.0, 14.0, 10.0)));
    add("bored-plate", &|m| {
        let plate = make_box(m, Point::ORIGIN, (20.0, 20.0, 5.0))?;
        let bore = cylinder(
            m,
            Point::new(10.0, 10.0, -1.0),
            Direction::Z,
            Direction::X,
            3.0,
            7.0,
        )?;
        ogeom::boolean::cut(m, &plate, &bore, T)
            .ok()
            .map(|b| b.shape)
    });
    add("shaved-cube", &|m| {
        let block = make_box(m, Point::ORIGIN, (10.0, 10.0, 10.0))?;
        let drum = cylinder(
            m,
            Point::new(5.0, 5.0, -1.0),
            Direction::Z,
            Direction::X,
            6.0,
            12.0,
        )?;
        ogeom::boolean::common(m, &block, &drum, T)
            .ok()
            .map(|b| b.shape)
    });
    add("rounded-box", &|m| {
        let block = make_box(m, Point::ORIGIN, (20.0, 14.0, 10.0))?;
        let edges = explore_unique(m, &block, ShapeType::Edge).ok()?;
        ogeom::fillet::fillet_edges(m, &block, &edges, 2.0, T)
            .ok()
            .map(|b| b.shape)
    });
    add("stepped", &|m| {
        let base = make_box(m, Point::ORIGIN, (20.0, 10.0, 5.0))?;
        let step = make_box(m, Point::new(5.0, 0.0, 5.0), (10.0, 10.0, 5.0))?;
        ogeom::boolean::fuse(m, &base, &step, T)
            .ok()
            .map(|b| b.shape)
    });
    add("cross-bores", &|m| {
        let block = make_box(m, Point::ORIGIN, (20.0, 20.0, 20.0))?;
        let a = cylinder(
            m,
            Point::new(-1.0, 10.0, 10.0),
            Direction::X,
            Direction::Y,
            4.0,
            22.0,
        )?;
        let b = cylinder(
            m,
            Point::new(10.0, -1.0, 12.0),
            Direction::Y,
            Direction::Z,
            3.0,
            22.0,
        )?;
        let once = ogeom::boolean::cut(m, &block, &a, T).ok()?.shape;
        ogeom::boolean::cut(m, &once, &b, T).ok().map(|b| b.shape)
    });
    add("sphere", &|m| {
        ogeom::algo::make_sphere(m, Frame::WORLD, 10.0, T)
            .ok()
            .map(|b| b.shape)
    });
    add("torus", &|m| {
        ogeom::algo::make_torus(m, Frame::WORLD, 10.0, 3.0, T)
            .ok()
            .map(|b| b.shape)
    });
    add("frustum", &|m| {
        ogeom::algo::make_cone(m, Frame::WORLD, 6.0, 3.0, 10.0, T)
            .ok()
            .map(|b| b.shape)
    });
    add("converted-plate", &|m| {
        // A bored plate through a mesh and back: planar facets merged, the
        // bore recognized, its tolerances those of a converted part.
        let mut source = Model::new();
        let plate = make_box(&mut source, Point::ORIGIN, (20.0, 20.0, 5.0))?;
        let bore = cylinder(
            &mut source,
            Point::new(10.0, 10.0, -1.0),
            Direction::Z,
            Direction::X,
            3.0,
            7.0,
        )?;
        let part = ogeom::boolean::cut(&mut source, &plate, &bore, T)
            .ok()?
            .shape;
        let mesh =
            ogeom::mesh::triangulate(&source, &part, Deflection::with_chord(0.01).ok()?, T).ok()?;
        ogeom::algo::solid_from_mesh(m, &mesh, &ogeom::algo::MeshSolidOptions::default(), T)
            .ok()
            .map(|b| b.shape)
    });
    parts
}

/// Solids from the committed corpus of public test parts.
fn corpus_parts() -> Vec<Part> {
    let dir = format!("{}/../../tests/corpus", env!("CARGO_MANIFEST_DIR"));
    let mut parts = Vec::new();
    // Every part here drills in seconds. A drill through a vertex picks it by
    // its index among the healed part's vertices, so a repair that adds a
    // vertex to a part draws that part's drills anew.
    for file in [
        "nist_ftc_11_asme1_rb.stp",
        "nist_ctc_01_asme1_rd.stp",
        "nist_ctc_03_asme1_rc.stp",
        "nist_ftc_06_asme1_rd.stp",
    ] {
        let Ok(text) = std::fs::read_to_string(format!("{dir}/{file}")) else {
            eprintln!("warning: corpus part {file} not found");
            continue;
        };
        let Ok(import) = ogeom::io::read_step(&text, T) else {
            eprintln!("warning: corpus part {file} did not read");
            continue;
        };
        let Some(solid) = import.solids.first().cloned() else {
            continue;
        };
        // Imported the way a caller imports: periodic faces' rings
        // re-anchored on their seams before anything cuts them.
        let mut import = import;
        let Ok((healed, _)) =
            ogeom::heal::reanchor_periodic_rings(import.document.model_mut(), &solid, T)
        else {
            eprintln!("warning: corpus part {file} did not heal");
            continue;
        };
        let name = file.trim_end_matches(".stp");
        if let Some(part) = Part::new(name, import.document.model().clone(), healed.shape) {
            parts.push(part);
        }
    }
    parts
}

// --- scenarios -------------------------------------------------------------

/// A unit vector square to `axis`, one of two, chosen at random.
fn across(axis: Vector, rng: &mut Rng) -> Vector {
    let helper = if axis.x.abs() < 0.9 {
        Vector::new(1.0, 0.0, 0.0)
    } else {
        Vector::new(0.0, 1.0, 0.0)
    };
    let u = axis.cross(helper).normalized(T).unwrap();
    let w = axis.cross(u).normalized(T).unwrap();
    if rng.chance(0.5) { u } else { w }
}

/// A drill through the part: along a principal axis mostly, obliquely
/// sometimes, placed at random or snapped so its wall passes through a
/// vertex or is tangent to a plane through one. The cut and the common add
/// up to the part, and the drill less the part and the common to the drill.
fn drill_case(part: &Part, rng: &mut Rng) -> Verdict {
    let mut model = part.model.clone();
    let diagonal = part.diagonal();
    let axis = if rng.chance(0.8) {
        *rng.pick(&[
            Vector::new(1.0, 0.0, 0.0),
            Vector::new(0.0, 1.0, 0.0),
            Vector::new(0.0, 0.0, 1.0),
        ])
    } else {
        Vector::new(
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
        )
        .normalized(T)
        .unwrap_or(Vector::new(0.0, 0.0, 1.0))
    };
    let radius = diagonal * rng.range(0.02, 0.15);
    let centre = part.low.midpoint(part.high);
    let through = if !part.vertices.is_empty() && rng.chance(0.5) {
        // Degenerate: the wall through a vertex, or tangent to a plane
        // through it.
        let v = *rng.pick(&part.vertices);
        let side = across(axis, rng);
        let offset = *rng.pick(&[1.0, -1.0, 0.0]);
        v + side * (radius * offset)
    } else {
        let jitter = Vector::new(
            rng.range(part.low.x, part.high.x) - centre.x,
            rng.range(part.low.y, part.high.y) - centre.y,
            rng.range(part.low.z, part.high.z) - centre.z,
        );
        centre + jitter
    };
    // Start before the part along the axis and run past it.
    let back = (through - part.low)
        .dot(axis)
        .abs()
        .max((through - part.high).dot(axis).abs());
    let start = through - axis * (back + diagonal * 0.1);
    let length = 2.0 * (back + diagonal * 0.1);
    let Ok(z) = Direction::new(axis, T) else {
        return Verdict::of(Outcome::Refused, "degenerate axis");
    };
    let Ok(x) = Direction::new(across(axis, rng), T) else {
        return Verdict::of(Outcome::Refused, "degenerate axis");
    };
    let drill =
        match ogeom::algo::make_cylinder(&mut model, frame_at(start, z, x), radius, length, T) {
            Ok(b) => b.shape,
            Err(e) => return Verdict::of(Outcome::Refused, format!("drill: {e}")),
        };
    // Where the drill went, for replaying a bad case outside the harness.
    let seam = x.vector();
    let placed = format!(
        "drill r {radius} from ({}, {}, {}) along ({}, {}, {}) length {length}, seam toward \
         ({}, {}, {})",
        start.x, start.y, start.z, axis.x, axis.y, axis.z, seam.x, seam.y, seam.z
    );
    let verdict = (|| {
        let drill_volume = core::f64::consts::PI * radius * radius * length;
        let cut = ogeom::boolean::cut(&mut model, &part.shape, &drill, T);
        let (_, cut) = match judged(&model, "cut", cut) {
            Ok(r) => r,
            Err(v) => return v,
        };
        let whole = part.measure;
        let scale = whole.volume.max(drill_volume.min(whole.volume));
        let common = ogeom::boolean::common(&mut model, &part.shape, &drill, T);
        let common = match common {
            // Nothing in common is an empty answer, not a refusal.
            Err(_) if cut.volume >= whole.volume - allowance(scale, &[&whole, &cut]) => {
                Measure::NOTHING
            }
            other => match judged(&model, "common", other) {
                Ok((_, m)) => m,
                Err(v) => return v,
            },
        };
        // The drill runs out of the part at both ends, so something of it
        // is always left.
        let rest = ogeom::boolean::cut(&mut model, &drill, &part.shape, T);
        let (_, rest) = match judged(&model, "drill less the part", rest) {
            Ok(r) => r,
            Err(v) => return v,
        };
        let named = [
            ("part", &whole),
            ("cut", &cut),
            ("common", &common),
            ("drill less the part", &rest),
        ];
        let slack = allowance(scale, &[&whole, &cut, &common]);
        let miss = cut.volume + common.volume - whole.volume;
        if miss.abs() > slack {
            return Verdict::of(
                Outcome::Wrong,
                format!(
                    "cut {:.9} + common {:.9} = {:.9}, part {:.9}, off {miss:.3e} \
                     ({:.2e} of the part) against {slack:.3e}, the part's stated doubt \
                     {:.3e}",
                    cut.volume,
                    common.volume,
                    cut.volume + common.volume,
                    whole.volume,
                    miss / whole.volume,
                    part.stated,
                ),
            )
            .read_from(&named);
        }
        let slack = allowance(drill_volume, &[&rest, &common]);
        let miss = rest.volume + common.volume - drill_volume;
        if miss.abs() > slack {
            return Verdict::of(
                Outcome::Wrong,
                format!(
                    "drill less the part {:.9} + common {:.9} = {:.9}, drill {drill_volume:.9}, \
                     off {miss:.3e} ({:.2e} of the drill) against {slack:.3e}",
                    rest.volume,
                    common.volume,
                    rest.volume + common.volume,
                    miss / drill_volume,
                ),
            )
            .read_from(&named);
        }
        Verdict::ok().read_from(&named)
    })();
    // Anything short of ok carries where the drill went, to replay it.
    if verdict.outcome == Outcome::Ok {
        verdict
    } else {
        Verdict {
            note: format!("{}; {placed}", verdict.note),
            ..verdict
        }
    }
}

/// A primitive of a random kind and size, placed with its position and
/// size snapped to a coarse grid half the time, so pairs meet face on face.
/// Its axis is one of the world's, or, `oblique`, any direction at all.
fn primitive(model: &mut Model, rng: &mut Rng, oblique: bool) -> Option<(Shape, Measure, String)> {
    let snap = |rng: &mut Rng, lo: f64, hi: f64| {
        let v = rng.range(lo, hi);
        if rng.chance(0.5) {
            (v * 2.0).round() / 2.0
        } else {
            v
        }
    };
    let at = Point::new(
        snap(rng, -3.0, 3.0),
        snap(rng, -3.0, 3.0),
        snap(rng, -3.0, 3.0),
    );
    let (axis, x) = if oblique {
        let z = Vector::new(
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
        );
        let z = Direction::new(z, T).ok()?;
        let x = Direction::new(across(z.vector(), rng), T).ok()?;
        (z, x)
    } else {
        let axis = *rng.pick(&[Direction::X, Direction::Y, Direction::Z]);
        let x = if axis == Direction::X {
            Direction::Y
        } else {
            Direction::X
        };
        (axis, x)
    };
    let frame = frame_at(at, axis, x);
    let size = |rng: &mut Rng| {
        let s = snap(rng, 2.0, 8.0);
        s.max(1.0)
    };
    let kind = rng.below(5);
    let (built, name, sizes) = match kind {
        0 => {
            let (a, b, c) = (size(rng), size(rng), size(rng));
            (
                ogeom::algo::make_box(model, frame, (a, b, c), T),
                "box",
                vec![a, b, c],
            )
        }
        1 => {
            let (r, h) = (size(rng) * 0.5, size(rng));
            (
                ogeom::algo::make_cylinder(model, frame, r, h, T),
                "cylinder",
                vec![r, h],
            )
        }
        2 => {
            let r = size(rng) * 0.5;
            (
                ogeom::algo::make_sphere(model, frame, r, T),
                "sphere",
                vec![r],
            )
        }
        3 => {
            let (r0, h) = (size(rng) * 0.5, size(rng));
            let r1 = r0 * rng.range(0.2, 0.9);
            (
                ogeom::algo::make_cone(model, frame, r0, r1, h, T),
                "cone",
                vec![r0, r1, h],
            )
        }
        _ => {
            let r = size(rng) * 0.5;
            let tube = (r * rng.range(0.15, 0.45)).max(0.3);
            (
                ogeom::algo::make_torus(model, frame, r, tube, T),
                "torus",
                vec![r, tube],
            )
        }
    };
    let shape = built.ok()?.shape;
    let m = measure(model, &shape).ok()?;
    let placed = format!(
        "{name} {sizes:?} at ({}, {}, {}) along ({}, {}, {}) seam toward ({}, {}, {})",
        at.x,
        at.y,
        at.z,
        axis.x(),
        axis.y(),
        axis.z(),
        x.x(),
        x.y(),
        x.z()
    );
    Some((shape, m, placed))
}

/// A box, drum, ball or torus on an oblique axis offset in or out by
/// `offset_shape`, which keeps its topology one for one: sharp edges stay
/// sharp, so the offset volume is the primitive's own formula at the grown
/// sizes.
fn offset_case(rng: &mut Rng) -> Verdict {
    use core::f64::consts::PI;
    let mut model = Model::new();
    let at = Point::new(
        rng.range(-3.0, 3.0),
        rng.range(-3.0, 3.0),
        rng.range(-3.0, 3.0),
    );
    let z = Vector::new(
        rng.range(-1.0, 1.0),
        rng.range(-1.0, 1.0),
        rng.range(-1.0, 1.0),
    );
    let Ok(z) = Direction::new(z, T) else {
        return Verdict::of(Outcome::Refused, "degenerate axis");
    };
    let Ok(x) = Direction::new(across(z.vector(), rng), T) else {
        return Verdict::of(Outcome::Refused, "degenerate axis");
    };
    let frame = frame_at(at, z, x);
    let kind = rng.below(4);
    let (built, smallest, name): (ogeom::core::OgeomResult<ogeom::algo::Built>, f64, String) =
        match kind {
            0 => {
                let (a, b, c) = (
                    rng.range(1.0, 6.0),
                    rng.range(1.0, 6.0),
                    rng.range(1.0, 6.0),
                );
                (
                    ogeom::algo::make_box(&mut model, frame, (a, b, c), T),
                    a.min(b).min(c) / 2.0,
                    format!("box {a} {b} {c}"),
                )
            }
            1 => {
                let (r, h) = (rng.range(0.5, 3.0), rng.range(1.0, 6.0));
                (
                    ogeom::algo::make_cylinder(&mut model, frame, r, h, T),
                    r.min(h / 2.0),
                    format!("cylinder {r} {h}"),
                )
            }
            2 => {
                let r = rng.range(0.5, 3.0);
                (
                    ogeom::algo::make_sphere(&mut model, frame, r, T),
                    r,
                    format!("sphere {r}"),
                )
            }
            _ => {
                let r = rng.range(1.5, 4.0);
                let tube = r * rng.range(0.1, 0.4);
                (
                    ogeom::algo::make_torus(&mut model, frame, r, tube, T),
                    tube,
                    format!("torus {r} {tube}"),
                )
            }
        };
    let Ok(built) = built else {
        return Verdict::of(Outcome::Refused, "primitive");
    };
    let shape = built.shape;
    let d = if rng.chance(0.5) {
        rng.range(0.05, 1.0)
    } else {
        -smallest * rng.range(0.05, 0.8)
    };
    let what = format!("{name} at {at:?} along {z:?} seam {x:?} offset {d}");
    let grown = ogeom::offset::offset_shape(&mut model, &shape, d, T);
    let (_, m) = match judged(&model, &what, grown) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let words: Vec<f64> = name
        .split(' ')
        .skip(1)
        .map(|w| w.parse().unwrap_or(0.0))
        .collect();
    let want = match kind {
        0 => (words[0] + 2.0 * d) * (words[1] + 2.0 * d) * (words[2] + 2.0 * d),
        1 => PI * (words[0] + d).powi(2) * (words[1] + 2.0 * d),
        2 => 4.0 / 3.0 * PI * (words[0] + d).powi(3),
        _ => 2.0 * PI * PI * words[0] * (words[1] + d).powi(2),
    };
    let slack = allowance(want.abs(), &[&m]);
    if (m.volume - want).abs() > slack {
        return Verdict::of(
            Outcome::Wrong,
            format!(
                "{what}: volume {:.9} against {want:.9}, off {:.3e}",
                m.volume,
                m.volume - want
            ),
        );
    }
    Verdict::ok().read_from(&[("offset", &m)])
}

/// Two boxes on oblique axes fused, and one to three of the fused solid's
/// edges between two planes filleted at a small radius: the result is
/// valid and holds the fused volume within what the rounds can take.
fn blend_case(rng: &mut Rng) -> Verdict {
    let mut model = Model::new();
    let mut boxes = Vec::new();
    for _ in 0..2 {
        let at = Point::new(
            rng.range(-2.0, 2.0),
            rng.range(-2.0, 2.0),
            rng.range(-2.0, 2.0),
        );
        let z = Vector::new(
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
            rng.range(-1.0, 1.0),
        );
        let Ok(z) = Direction::new(z, T) else {
            return Verdict::of(Outcome::Refused, "degenerate axis");
        };
        let Ok(x) = Direction::new(across(z.vector(), rng), T) else {
            return Verdict::of(Outcome::Refused, "degenerate axis");
        };
        let size = (
            rng.range(2.0, 6.0),
            rng.range(2.0, 6.0),
            rng.range(2.0, 6.0),
        );
        let Ok(b) = ogeom::algo::make_box(&mut model, frame_at(at, z, x), size, T) else {
            return Verdict::of(Outcome::Refused, "box");
        };
        boxes.push((
            b.shape,
            format!("box {size:?} at {at:?} along {z:?} x {x:?}"),
        ));
    }
    let what = format!("{} with {}", boxes[0].1, boxes[1].1);
    let fused = ogeom::boolean::fuse(&mut model, &boxes[0].0, &boxes[1].0, T);
    let (fused, mf) = match judged(&model, &format!("{what}: fuse"), fused) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let edges = ogeom::topo::explore_unique(&model, &fused, ShapeType::Edge).unwrap_or_default();
    if edges.is_empty() {
        return Verdict::of(Outcome::Refused, "no edges");
    }
    let count = 1 + rng.below(3);
    let mut chosen: Vec<Shape> = Vec::new();
    for _ in 0..count {
        let e = rng.pick(&edges).clone();
        if !chosen.iter().any(|c| c.node() == e.node()) {
            chosen.push(e);
        }
    }
    let radius = rng.range(0.05, 0.4);
    let result = ogeom::fillet::fillet_edges(&mut model, &fused, &chosen, radius, T);
    let (_, m) = match judged(&model, &format!("{what}: fillet r {radius}"), result) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let change = (m.volume - mf.volume).abs() / mf.volume;
    if change > 0.05 {
        return Verdict::of(
            Outcome::Wrong,
            format!(
                "{what}: fillet r {radius:.4} moved {:.2}% of the part",
                change * 100.0
            ),
        );
    }
    Verdict::ok().read_from(&[("fused", &mf), ("filleted", &m)])
}

/// An oblique box with edges rounded or bevelled: one edge, two that share
/// no corner, or all twelve; or a drum with one rim or both. Each has a
/// closed-form volume: a round takes `r^2 (1 - pi/4)` along its edge and a
/// bevel half a square of its distance, and where all twelve meet in the
/// corners the box is Steiner's (rounds) or keeps the corners' triple
/// overlap (bevels).
fn round_case(rng: &mut Rng) -> Verdict {
    use core::f64::consts::PI;
    let mut model = Model::new();
    let at = Point::new(
        rng.range(-3.0, 3.0),
        rng.range(-3.0, 3.0),
        rng.range(-3.0, 3.0),
    );
    let z = Vector::new(
        rng.range(-1.0, 1.0),
        rng.range(-1.0, 1.0),
        rng.range(-1.0, 1.0),
    );
    let Ok(z) = Direction::new(z, T) else {
        return Verdict::of(Outcome::Refused, "degenerate axis");
    };
    let Ok(x) = Direction::new(across(z.vector(), rng), T) else {
        return Verdict::of(Outcome::Refused, "degenerate axis");
    };
    let dims = (
        rng.range(1.0, 6.0),
        rng.range(1.0, 6.0),
        rng.range(1.0, 6.0),
    );
    if rng.chance(0.3) {
        return drum_round(&mut model, frame_at(at, z, x), dims, rng);
    }
    if rng.chance(0.2) {
        return bore_round(&mut model, frame_at(at, z, x), dims, rng);
    }
    let Ok(built) = ogeom::algo::make_box(&mut model, frame_at(at, z, x), dims, T) else {
        return Verdict::of(Outcome::Refused, "box");
    };
    let shape = built.shape;
    let edges = ogeom::topo::explore_unique(&model, &shape, ShapeType::Edge).unwrap_or_default();
    let ends = |e: &Shape| -> Vec<Point> {
        ogeom::topo::explore_unique(&model, e, ShapeType::Vertex)
            .unwrap_or_default()
            .iter()
            .filter_map(|v| {
                model
                    .node(v)
                    .and_then(|n| n.data().as_vertex())
                    .map(|d| d.point)
            })
            .collect()
    };
    let all = rng.chance(0.3);
    let first = rng.pick(&edges).clone();
    let mut chosen = vec![first.clone()];
    if all {
        chosen = edges.clone();
    } else if rng.chance(0.5) {
        let near = ends(&first);
        let apart: Vec<Shape> = edges
            .iter()
            .filter(|e| {
                ends(e)
                    .iter()
                    .all(|p| near.iter().all(|q| q.distance(*p) > 1e-6))
            })
            .cloned()
            .collect();
        if !apart.is_empty() {
            chosen.push(rng.pick(&apart).clone());
        }
    }
    let r = dims.0.min(dims.1).min(dims.2) * rng.range(0.05, 0.45);
    let length: f64 = chosen
        .iter()
        .map(|e| {
            let p = ends(e);
            p[0].distance(p[1])
        })
        .sum();
    let bevel = rng.chance(0.3);
    let what = format!(
        "box {dims:?} at {at:?} along {z:?} x {x:?}, {} edges {} {r}",
        chosen.len(),
        if bevel { "bevelled" } else { "rounded" }
    );
    // Short of every edge, a bevel may run a second distance along the
    // other face, given outright or by its angle from the first; on a
    // square edge either takes a triangle of the two legs.
    let lean = if bevel && !all { rng.below(3) } else { 0 };
    let other = r * rng
        .range(0.5, 1.5)
        .min(dims.0.min(dims.1).min(dims.2) * 0.45 / r);
    let angle = (other / r).atan();
    let what = match lean {
        1 => format!("{what}, {other} on the other face"),
        2 => format!("{what}, at {angle} rad"),
        _ => what,
    };
    let result = if lean > 0 {
        let mut specs = Vec::with_capacity(chosen.len());
        for edge in &chosen {
            let faces =
                ogeom::topo::explore_unique(&model, &shape, ShapeType::Face).unwrap_or_default();
            let Some(face) = faces.into_iter().find(|f| {
                ogeom::topo::explore_unique(&model, f, ShapeType::Edge)
                    .unwrap_or_default()
                    .iter()
                    .any(|e| e.is_same(edge))
            }) else {
                return Verdict::of(Outcome::Wrong, format!("{what}: an edge on no face"));
            };
            specs.push((
                edge.clone(),
                if lean == 1 {
                    ogeom::fillet::Chamfer::Distances {
                        face,
                        on_face: r,
                        on_other: other,
                    }
                } else {
                    ogeom::fillet::Chamfer::Angle {
                        face,
                        distance: r,
                        angle,
                    }
                },
            ));
        }
        ogeom::fillet::chamfer_edges_with(&mut model, &shape, &specs, T)
    } else if bevel {
        ogeom::fillet::chamfer_edges(&mut model, &shape, &chosen, r, T)
    } else {
        ogeom::fillet::fillet_edges(&mut model, &shape, &chosen, r, T)
    };
    let (_, m) = match judged(&model, &what, result) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let want = if bevel {
        // Three bevels at a corner cut its pairs' overlaps twice over and
        // the triple's three times: d^3 less d^3/4 to give back.
        let corners = if all { 8.0 } else { 0.0 };
        let legs = if lean > 0 { r * other } else { r * r };
        (corners * 0.75 * r).mul_add(r * r, dims.0 * dims.1 * dims.2 - legs / 2.0 * length)
    } else if all {
        // The box shrunk by the radius, its faces' slabs, its edges'
        // quarter cylinders and its corners' ball.
        let (a, b, c) = (dims.0 - 2.0 * r, dims.1 - 2.0 * r, dims.2 - 2.0 * r);
        a * b * c
            + 2.0 * r * (a * b + b * c + a * c)
            + PI * r * r * (a + b + c)
            + 4.0 / 3.0 * PI * r * r * r
    } else {
        dims.0 * dims.1 * dims.2 - (1.0 - PI / 4.0) * r * r * length
    };
    let slack = allowance(want, &[&m]);
    if (m.volume - want).abs() > slack {
        return Verdict::of(
            Outcome::Wrong,
            format!(
                "{what}: volume {:.9} against {want:.9}, off {:.3e}",
                m.volume,
                m.volume - want
            ),
        );
    }
    Verdict::ok().read_from(&[("result", &m)])
}

/// A drum with one rim rounded or bevelled, or both: each takes the corner
/// between wall and cap turned about the axis, by Pappus. A round's
/// corner has area `r^2 (1 - pi/4)` and its centroid `r (10 - 3 pi) /
/// (3 (4 - pi))` in from the rim; a bevel's, `r^2 / 2` and `r / 3`.
fn drum_round(model: &mut Model, frame: Frame, dims: (f64, f64, f64), rng: &mut Rng) -> Verdict {
    use core::f64::consts::PI;
    let (radius, height) = (dims.0 / 2.0, dims.2);
    let Ok(built) = ogeom::algo::make_cylinder(model, frame, radius, height, T) else {
        return Verdict::of(Outcome::Refused, "drum");
    };
    let shape = built.shape;
    let rims: Vec<Shape> = ogeom::topo::explore_unique(model, &shape, ShapeType::Edge)
        .unwrap_or_default()
        .into_iter()
        .filter(|e| {
            ogeom::algo::edge_vertices(model, e)
                .ok()
                .flatten()
                .is_some_and(|(a, b)| a.is_same(&b))
        })
        .collect();
    if rims.len() != 2 {
        return Verdict::of(Outcome::Wrong, format!("a drum with {} rims", rims.len()));
    }
    let chosen: Vec<Shape> = if rng.chance(0.5) {
        rims
    } else {
        rims[..1].to_vec()
    };
    let r = radius.min(height / 2.0) * rng.range(0.05, 0.45);
    let bevel = rng.chance(0.3);
    let what = format!(
        "drum r {radius} h {height} at {:?} along {:?} x {:?}, {} rims {} {r}",
        frame.origin(),
        frame.z(),
        frame.x(),
        chosen.len(),
        if bevel { "bevelled" } else { "rounded" }
    );
    let result = if bevel {
        ogeom::fillet::chamfer_edges(model, &shape, &chosen, r, T)
    } else {
        ogeom::fillet::fillet_edges(model, &shape, &chosen, r, T)
    };
    let (_, m) = match judged(model, &what, result) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let removed = if bevel {
        PI * (radius - r / 3.0) * r * r
    } else {
        let k = (10.0 - 3.0 * PI) / (3.0 * (4.0 - PI));
        2.0 * PI * (radius - k * r) * (1.0 - PI / 4.0) * r * r
    };
    #[allow(clippy::cast_precision_loss)]
    let want = PI * radius * radius * height - removed * chosen.len() as f64;
    let slack = allowance(want, &[&m]);
    if (m.volume - want).abs() > slack {
        return Verdict::of(
            Outcome::Wrong,
            format!(
                "{what}: volume {:.9} against {want:.9}, off {:.3e}",
                m.volume,
                m.volume - want
            ),
        );
    }
    Verdict::ok().read_from(&[("result", &m)])
}

/// A box bored through along its own `z`, one rim of the bore rounded or
/// bevelled, or both: each takes the corner between bore and face turned
/// about the bore's axis, its centroid outside the bore's radius.
fn bore_round(model: &mut Model, frame: Frame, dims: (f64, f64, f64), rng: &mut Rng) -> Verdict {
    use core::f64::consts::PI;
    let Ok(block) = ogeom::algo::make_box(model, frame, dims, T) else {
        return Verdict::of(Outcome::Refused, "box");
    };
    let radius = dims.0.min(dims.1) * rng.range(0.1, 0.3);
    let centre =
        frame.origin() + frame.x().vector() * (dims.0 / 2.0) + frame.y().vector() * (dims.1 / 2.0)
            - frame.z().vector();
    let Ok(axis) = Frame::new(centre, frame.z(), frame.x(), T) else {
        return Verdict::of(Outcome::Refused, "bore frame");
    };
    let Ok(drill) = ogeom::algo::make_cylinder(model, axis, radius, dims.2 + 2.0, T) else {
        return Verdict::of(Outcome::Refused, "drill");
    };
    let what = format!(
        "box {dims:?} at {:?} along {:?} x {:?} bored r {radius}",
        frame.origin(),
        frame.z(),
        frame.x()
    );
    let bored = ogeom::boolean::cut(model, &block.shape, &drill.shape, T);
    let (bored, mb) = match judged(model, &format!("{what}: bore"), bored) {
        Ok(r) => r,
        Err(v) => return v,
    };
    // The bore's rims, each the arcs the cut left of one circle: told
    // apart by which face of the box they lie in.
    let mut rims: [Vec<Shape>; 2] = [Vec::new(), Vec::new()];
    for e in ogeom::topo::explore_unique(model, &bored, ShapeType::Edge).unwrap_or_default() {
        let circle = model
            .node(&e)
            .and_then(|n| n.data().as_edge())
            .and_then(|d| match d.curve3d() {
                Some(ogeom::topo::EdgeRepr::Curve3d { curve, .. }) => {
                    model.geometry().curve(*curve)
                }
                _ => None,
            })
            .and_then(|g| match g {
                ogeom::geom::Curve::Circle(c) => Some(c.circle()),
                _ => None,
            });
        let Some(circle) = circle else { continue };
        let Ok(placement) = e.transform(model.datums()) else {
            continue;
        };
        let height = (placement.apply(circle.centre()) - frame.origin()).dot(frame.z().vector());
        rims[usize::from(height > dims.2 / 2.0)].push(e);
    }
    if rims.iter().any(Vec::is_empty) {
        return Verdict::of(Outcome::Wrong, format!("{what}: a rim with no arcs"));
    }
    let both = rng.chance(0.5);
    let count = if both { 2.0 } else { 1.0 };
    let chosen: Vec<Shape> = if both { rims.concat() } else { rims[0].clone() };
    let room = (dims.0.min(dims.1) / 2.0 - radius)
        .min(dims.2 / 2.0)
        .min(radius);
    let r = room * rng.range(0.05, 0.45);
    let bevel = rng.chance(0.3);
    let what = format!(
        "{what}, {count} rims {} {r}",
        if bevel { "bevelled" } else { "rounded" }
    );
    let result = if bevel {
        ogeom::fillet::chamfer_edges(model, &bored, &chosen, r, T)
    } else {
        ogeom::fillet::fillet_edges(model, &bored, &chosen, r, T)
    };
    let (_, m) = match judged(model, &what, result) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let removed = if bevel {
        PI * (radius + r / 3.0) * r * r
    } else {
        let k = (10.0 - 3.0 * PI) / (3.0 * (4.0 - PI));
        2.0 * PI * (radius + k * r) * (1.0 - PI / 4.0) * r * r
    };
    #[allow(clippy::cast_precision_loss)]
    let want = dims.0 * dims.1 * dims.2 - PI * radius * radius * dims.2 - removed * count;
    let slack = allowance(want, &[&mb, &m]);
    if (m.volume - want).abs() > slack {
        return Verdict::of(
            Outcome::Wrong,
            format!(
                "{what}: volume {:.9} against {want:.9}, off {:.3e}",
                m.volume,
                m.volume - want
            ),
        );
    }
    Verdict::ok().read_from(&[("bored", &mb), ("result", &m)])
}

/// A slanted box rounded on every edge, cut by and in common with an
/// oblique primitive: the two add up to the rounded box, whose volume is
/// Steiner's, so booleans on blends (cylinders, tori, corner balls) are
/// held to an exact total.
fn rounded_case(rng: &mut Rng) -> Verdict {
    use core::f64::consts::PI;
    let mut model = Model::new();
    let at = Point::new(
        rng.range(-3.0, 1.0),
        rng.range(-3.0, 1.0),
        rng.range(-3.0, 1.0),
    );
    let z = Vector::new(
        rng.range(-1.0, 1.0),
        rng.range(-1.0, 1.0),
        rng.range(-1.0, 1.0),
    );
    let Ok(z) = Direction::new(z, T) else {
        return Verdict::of(Outcome::Refused, "degenerate axis");
    };
    let Ok(x) = Direction::new(across(z.vector(), rng), T) else {
        return Verdict::of(Outcome::Refused, "degenerate axis");
    };
    let dims = (
        rng.range(2.0, 6.0),
        rng.range(2.0, 6.0),
        rng.range(2.0, 6.0),
    );
    let Ok(built) = ogeom::algo::make_box(&mut model, frame_at(at, z, x), dims, T) else {
        return Verdict::of(Outcome::Refused, "box");
    };
    let edges =
        ogeom::topo::explore_unique(&model, &built.shape, ShapeType::Edge).unwrap_or_default();
    let r = dims.0.min(dims.1).min(dims.2) * rng.range(0.05, 0.4);
    let what = format!("box {dims:?} at {at:?} along {z:?} x {x:?} rounded {r}");
    let rounded = ogeom::fillet::fillet_edges(&mut model, &built.shape, &edges, r, T);
    let (rounded, mr) = match judged(&model, &what, rounded) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let (a, b, c) = (dims.0 - 2.0 * r, dims.1 - 2.0 * r, dims.2 - 2.0 * r);
    let exact = a * b * c
        + 2.0 * r * (a * b + b * c + a * c)
        + PI * r * r * (a + b + c)
        + 4.0 / 3.0 * PI * r * r * r;
    let Some((tool, _, tool_what)) = primitive(&mut model, rng, true) else {
        return Verdict::of(Outcome::Refused, "primitive");
    };
    let what = format!("{what} against {tool_what}");
    let cut = ogeom::boolean::cut(&mut model, &rounded, &tool, T);
    let (_, cut) = match judged(&model, &format!("{what}: cut"), cut) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let common = ogeom::boolean::common(&mut model, &rounded, &tool, T);
    let common = match common {
        Err(_) if (cut.volume - exact).abs() <= allowance(exact, &[&mr, &cut]) => Measure::NOTHING,
        other => match judged(&model, &format!("{what}: common"), other) {
            Ok((_, m)) => m,
            Err(v) => return v,
        },
    };
    let slack = allowance(exact, &[&mr, &cut, &common]);
    let miss = cut.volume + common.volume - exact;
    if miss.abs() > slack {
        return Verdict::of(
            Outcome::Wrong,
            format!(
                "{what}: cut {:.9} + common {:.9} against {exact:.9}, off {miss:.3e} against {slack:.3e}",
                cut.volume, common.volume
            ),
        );
    }
    Verdict::ok().read_from(&[("rounded", &mr), ("cut", &cut), ("common", &common)])
}

/// An oblique box or drum hollowed by `make_thick_solid`, opened at one of
/// its flat faces: the cavity is the solid shrunk by the wall on every side
/// but the opened one, where it runs out flush.
fn hollow_case(rng: &mut Rng) -> Verdict {
    use core::f64::consts::PI;
    let mut model = Model::new();
    let at = Point::new(
        rng.range(-3.0, 3.0),
        rng.range(-3.0, 3.0),
        rng.range(-3.0, 3.0),
    );
    let z = Vector::new(
        rng.range(-1.0, 1.0),
        rng.range(-1.0, 1.0),
        rng.range(-1.0, 1.0),
    );
    let Ok(z) = Direction::new(z, T) else {
        return Verdict::of(Outcome::Refused, "degenerate axis");
    };
    let Ok(x) = Direction::new(across(z.vector(), rng), T) else {
        return Verdict::of(Outcome::Refused, "degenerate axis");
    };
    let frame = frame_at(at, z, x);
    let boxed = rng.chance(0.5);
    let dims = (
        rng.range(1.0, 6.0),
        rng.range(1.0, 6.0),
        rng.range(1.0, 6.0),
    );
    let built = if boxed {
        ogeom::algo::make_box(&mut model, frame, dims, T)
    } else {
        ogeom::algo::make_cylinder(&mut model, frame, dims.0 / 2.0, dims.2, T)
    };
    let Ok(built) = built else {
        return Verdict::of(Outcome::Refused, "primitive");
    };
    let shape = built.shape;
    let flats: Vec<(Shape, Vector)> = ogeom::topo::explore_unique(&model, &shape, ShapeType::Face)
        .unwrap_or_default()
        .into_iter()
        .filter(|f| {
            model
                .node(f)
                .and_then(|n| n.data().as_face())
                .and_then(|d| model.geometry().surface(d.surface))
                .is_some_and(|g| matches!(g, ogeom::geom::SurfaceGeometry::Plane(_)))
        })
        .filter_map(|f| {
            ogeom::algo::face_normal(&model, &f, T)
                .ok()
                .map(|(_, n)| (f, n))
        })
        .collect();
    if flats.is_empty() {
        return Verdict::of(Outcome::Wrong, "a primitive with no flat face");
    }
    let (opened, normal) = rng.pick(&flats).clone();
    let smallest = if boxed {
        dims.0.min(dims.1).min(dims.2)
    } else {
        dims.0.min(dims.2)
    };
    let wall = smallest * rng.range(0.05, 0.3);
    let what = format!(
        "{} {dims:?} at {at:?} along {z:?} x {x:?}, opened along {normal:?}, wall {wall}",
        if boxed { "box" } else { "drum" }
    );
    let result = ogeom::offset::make_thick_solid(&mut model, &shape, &[opened], wall, T);
    let (_, m) = match judged(&model, &what, result) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let want = if boxed {
        // Short by one wall along the opened face's normal, by two along
        // the others.
        let along = [x.vector(), frame.y().vector(), z.vector()];
        let sides = [dims.0, dims.1, dims.2];
        let cavity: f64 = (0..3)
            .map(|i| {
                let walls = if normal.dot(along[i]).abs() > 0.5 {
                    1.0
                } else {
                    2.0
                };
                sides[i] - walls * wall
            })
            .product();
        dims.0 * dims.1 * dims.2 - cavity
    } else {
        let (radius, height) = (dims.0 / 2.0, dims.2);
        PI * radius * radius * height - PI * (radius - wall).powi(2) * (height - wall)
    };
    let slack = allowance(want, &[&m]);
    if (m.volume - want).abs() > slack {
        return Verdict::of(
            Outcome::Wrong,
            format!(
                "{what}: volume {:.9} against {want:.9}, off {:.3e}",
                m.volume,
                m.volume - want
            ),
        );
    }
    Verdict::ok().read_from(&[("hollow", &m)])
}

/// Three primitives along any direction: the first two fused, and that
/// cut by and in common with the third. The cut and the common add up to
/// the fused pair, so a boolean on a boolean's result, sections and all, is
/// held to the same identity as one on primitives.
fn chain_case(rng: &mut Rng) -> Verdict {
    let mut model = Model::new();
    let mut parts = Vec::new();
    for _ in 0..3 {
        let Some(p) = primitive(&mut model, rng, true) else {
            return Verdict::of(Outcome::Refused, "primitive");
        };
        parts.push(p);
    }
    let what = format!(
        "({} with {}) against {}",
        parts[0].2, parts[1].2, parts[2].2
    );
    let fused = ogeom::boolean::fuse(&mut model, &parts[0].0, &parts[1].0, T);
    let (fused, mf) = match judged(&model, &format!("{what}: first fuse"), fused) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let c = parts[2].0.clone();
    let common = ogeom::boolean::common(&mut model, &fused, &c, T);
    let common = match common {
        Err(_) => None,
        other => match judged(&model, &format!("{what}: common"), other) {
            Ok((_, m)) => Some(m),
            Err(v) => return v,
        },
    };
    let cut = ogeom::boolean::cut(&mut model, &fused, &c, T);
    let (_, cut) = match judged(&model, &format!("{what}: cut"), cut) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let common = match common {
        Some(m) => m,
        None if (cut.volume - mf.volume).abs() <= allowance(mf.volume, &[&mf, &cut]) => {
            Measure::NOTHING
        }
        None => return Verdict::of(Outcome::Refused, format!("{what}: common refused")),
    };
    let slack = allowance(mf.volume + parts[2].1.volume, &[&mf, &cut, &common]);
    let miss = cut.volume + common.volume - mf.volume;
    if miss.abs() > slack {
        return Verdict::of(
            Outcome::Wrong,
            format!(
                "{what}: cut {:.9} + common {:.9} against {:.9}, off {miss:.3e} against {slack:.3e}",
                cut.volume, common.volume, mf.volume
            ),
        );
    }
    Verdict::ok().read_from(&[("fused", &mf), ("cut", &cut), ("common", &common)])
}

/// Two primitives, and the boolean identities: the fuse and the common
/// add up to the two, and the cut is the first less the common.
fn pair_case(rng: &mut Rng, oblique: bool) -> Verdict {
    let mut model = Model::new();
    let Some((a, ma, ka)) = primitive(&mut model, rng, oblique) else {
        return Verdict::of(Outcome::Refused, "primitive");
    };
    let Some((b, mb, kb)) = primitive(&mut model, rng, oblique) else {
        return Verdict::of(Outcome::Refused, "primitive");
    };
    let what = format!("{ka} with {kb}");
    let common = ogeom::boolean::common(&mut model, &a, &b, T);
    let common = match common {
        Err(_) => None,
        other => match judged(&model, &format!("{what}: common"), other) {
            Ok((_, m)) => Some(m),
            Err(v) => return v,
        },
    };
    let fuse = ogeom::boolean::fuse(&mut model, &a, &b, T);
    let (_, fuse) = match judged(&model, &format!("{what}: fuse"), fuse) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let cut = ogeom::boolean::cut(&mut model, &a, &b, T);
    let (_, cut) = match judged(&model, &format!("{what}: cut"), cut) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let (va, vb) = (ma.volume, mb.volume);
    let scale = va + vb;
    // A refused common is only an empty one when the fuse says the two
    // do not overlap.
    let common = match common {
        Some(m) => m,
        None if (fuse.volume - va - vb).abs() <= allowance(scale, &[&ma, &mb, &fuse]) => {
            Measure::NOTHING
        }
        None => {
            return Verdict::of(Outcome::Refused, format!("{what}: common refused"));
        }
    };
    let named = [
        ("first", &ma),
        ("second", &mb),
        ("fuse", &fuse),
        ("common", &common),
        ("cut", &cut),
    ];
    let slack = allowance(scale, &[&ma, &mb, &fuse, &common]);
    let miss = fuse.volume + common.volume - va - vb;
    if miss.abs() > slack {
        return Verdict::of(
            Outcome::Wrong,
            format!(
                "{what}: fuse {:.9} + common {:.9} against {:.9}, off {miss:.3e} against \
                 {slack:.3e}",
                fuse.volume,
                common.volume,
                va + vb
            ),
        )
        .read_from(&named);
    }
    let slack = allowance(scale, &[&ma, &cut, &common]);
    let miss = cut.volume + common.volume - va;
    if miss.abs() > slack {
        return Verdict::of(
            Outcome::Wrong,
            format!(
                "{what}: cut {:.9} + common {:.9} against {va:.9}, off {miss:.3e} against \
                 {slack:.3e}",
                cut.volume, common.volume
            ),
        )
        .read_from(&named);
    }
    Verdict::ok().read_from(&named)
}

/// A fillet on one to three random edges of a part at a random radius: a
/// valid solid that moved material, but not half the part.
fn fillet_case(part: &Part, rng: &mut Rng) -> Verdict {
    let mut model = part.model.clone();
    let edges = &part.creases;
    let count = 1 + rng.below(3);
    let mut chosen: Vec<Shape> = Vec::new();
    for _ in 0..count {
        let e = rng.pick(edges).clone();
        if !chosen.iter().any(|c| c.node() == e.node()) {
            chosen.push(e);
        }
    }
    let radius = part.diagonal() * rng.range(0.005, 0.06);
    let result = ogeom::fillet::fillet_edges(&mut model, &part.shape, &chosen, radius, T);
    let (_, m) = match judged(&model, "fillet", result) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let change = (m.volume - part.measure.volume).abs() / part.measure.volume;
    if change > 0.25 {
        return Verdict::of(
            Outcome::Wrong,
            format!(
                "fillet r {radius:.4} moved {:.1}% of the part",
                change * 100.0
            ),
        );
    }
    Verdict::ok().read_from(&[("fillet", &m)])
}

// --- running ---------------------------------------------------------------

struct Case {
    key: String,
    run: Box<dyn Fn() -> Verdict + Send + Sync>,
}

struct Ran {
    key: String,
    scenario: String,
    verdict: Verdict,
    seconds: f64,
}

fn cases(seed: u64, per_part: usize, corpus: bool) -> Vec<Case> {
    let mut parts = generated_parts();
    let generated = parts.len();
    if corpus {
        parts.extend(corpus_parts());
    }
    let mut out: Vec<Case> = Vec::new();
    for (index, part) in parts.iter().enumerate() {
        // Corpus parts are large. They take a third of the drills.
        let n = if index < generated {
            per_part
        } else {
            (per_part / 3).max(1)
        };
        for i in 0..n {
            let key = format!("drill/{}/{i}", part.name);
            let part = part.clone();
            let s = case_seed(seed, &key);
            out.push(Case {
                key,
                run: Box::new(move || drill_case(&part, &mut Rng::new(s))),
            });
        }
    }
    for part in parts
        .iter()
        .take(generated)
        .filter(|p| !p.creases.is_empty())
    {
        for i in 0..per_part {
            let key = format!("fillet/{}/{i}", part.name);
            let part = part.clone();
            let s = case_seed(seed, &key);
            out.push(Case {
                key,
                run: Box::new(move || fillet_case(&part, &mut Rng::new(s))),
            });
        }
    }
    for i in 0..per_part * 3 {
        let key = format!("offset/primitives/{i}");
        let s = case_seed(seed, &key);
        out.push(Case {
            key,
            run: Box::new(move || offset_case(&mut Rng::new(s))),
        });
    }
    for i in 0..per_part * 3 {
        let key = format!("blend/boxes/{i}");
        let s = case_seed(seed, &key);
        out.push(Case {
            key,
            run: Box::new(move || blend_case(&mut Rng::new(s))),
        });
    }
    for i in 0..per_part * 3 {
        let key = format!("round/primitives/{i}");
        let s = case_seed(seed, &key);
        out.push(Case {
            key,
            run: Box::new(move || round_case(&mut Rng::new(s))),
        });
    }
    for i in 0..per_part * 3 {
        let key = format!("rounded/primitives/{i}");
        let s = case_seed(seed, &key);
        out.push(Case {
            key,
            run: Box::new(move || rounded_case(&mut Rng::new(s))),
        });
    }
    for i in 0..per_part * 3 {
        let key = format!("hollow/primitives/{i}");
        let s = case_seed(seed, &key);
        out.push(Case {
            key,
            run: Box::new(move || hollow_case(&mut Rng::new(s))),
        });
    }
    for i in 0..per_part * 3 {
        let key = format!("chain/primitives/{i}");
        let s = case_seed(seed, &key);
        out.push(Case {
            key,
            run: Box::new(move || chain_case(&mut Rng::new(s))),
        });
    }
    for (scenario, oblique) in [("pairs", false), ("oblique", true)] {
        for i in 0..per_part * 6 {
            let key = format!("{scenario}/primitives/{i}");
            let s = case_seed(seed, &key);
            out.push(Case {
                key,
                run: Box::new(move || pair_case(&mut Rng::new(s), oblique)),
            });
        }
    }
    out
}

fn run_all(cases: &[Case], threads: usize) -> Vec<Ran> {
    let next = AtomicUsize::new(0);
    let results: Mutex<Vec<Ran>> = Mutex::new(Vec::with_capacity(cases.len()));
    let done = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                loop {
                    let i = next.fetch_add(1, Ordering::SeqCst);
                    let Some(case) = cases.get(i) else {
                        break;
                    };
                    let start = Instant::now();
                    let verdict =
                        catch_unwind(AssertUnwindSafe(|| (case.run)())).unwrap_or_else(|payload| {
                            let message = payload
                                .downcast_ref::<String>()
                                .cloned()
                                .or_else(|| {
                                    payload.downcast_ref::<&str>().map(|s| (*s).to_string())
                                })
                                .unwrap_or_default();
                            Verdict::of(Outcome::Panic, message)
                        });
                    let scenario = case.key.split('/').next().unwrap_or("").to_string();
                    results.lock().unwrap().push(Ran {
                        key: case.key.clone(),
                        scenario,
                        verdict,
                        seconds: start.elapsed().as_secs_f64(),
                    });
                    let n = done.fetch_add(1, Ordering::SeqCst) + 1;
                    if n.is_multiple_of(50) {
                        eprintln!("  {n} of {} cases", cases.len());
                    }
                }
            });
        }
    });
    let mut out = results.into_inner().unwrap();
    out.sort_by(|a, b| a.key.cmp(&b.key));
    out
}

type Tally = BTreeMap<String, BTreeMap<&'static str, usize>>;

fn tally(ran: &[Ran]) -> Tally {
    let mut out: Tally = BTreeMap::new();
    for r in ran {
        let row = out.entry(r.scenario.clone()).or_default();
        for o in Outcome::ALL {
            row.entry(o.name()).or_insert(0);
        }
        *row.entry(r.verdict.outcome.name()).or_insert(0) += 1;
    }
    out
}

fn tally_json(tally: &Tally) -> String {
    let mut s = String::from("{\n");
    let rows: Vec<String> = tally
        .iter()
        .map(|(scenario, counts)| {
            let fields: Vec<String> = Outcome::ALL
                .iter()
                .map(|o| {
                    format!(
                        "\"{}\": {}",
                        o.name(),
                        counts.get(o.name()).copied().unwrap_or(0)
                    )
                })
                .collect();
            format!("  \"{scenario}\": {{ {} }}", fields.join(", "))
        })
        .collect();
    s.push_str(&rows.join(",\n"));
    s.push_str("\n}\n");
    s
}

/// The baseline's counts, read from the small JSON `tally_json` writes.
fn read_tally(text: &str) -> Tally {
    let mut out: Tally = BTreeMap::new();
    for line in text.lines() {
        let line = line.trim();
        let Some((scenario, rest)) = line.split_once(':') else {
            continue;
        };
        let scenario = scenario.trim().trim_matches('"').to_string();
        if !rest.contains('{') {
            continue;
        }
        let row = out.entry(scenario).or_default();
        for field in rest
            .trim()
            .trim_start_matches('{')
            .trim_end_matches(',')
            .trim_end_matches('}')
            .split(',')
        {
            if let Some((k, v)) = field.split_once(':') {
                let k = k.trim().trim_matches('"');
                if let (Some(o), Ok(v)) = (
                    Outcome::ALL.iter().find(|o| o.name() == k),
                    v.trim().parse::<usize>(),
                ) {
                    row.insert(o.name(), v);
                }
            }
        }
    }
    out
}

fn print_report(ran: &[Ran], tally: &Tally) {
    println!(
        "{:<10} {:>6} {:>8} {:>8} {:>6} {:>6} {:>7}",
        "scenario", "ok", "refused", "invalid", "wrong", "panic", "ok %"
    );
    for (scenario, row) in tally {
        let total: usize = row.values().sum();
        let get = |o: Outcome| row.get(o.name()).copied().unwrap_or(0);
        println!(
            "{:<10} {:>6} {:>8} {:>8} {:>6} {:>6} {:>6.1}%",
            scenario,
            get(Outcome::Ok),
            get(Outcome::Refused),
            get(Outcome::Invalid),
            get(Outcome::Wrong),
            get(Outcome::Panic),
            100.0 * get(Outcome::Ok) as f64 / total.max(1) as f64
        );
    }
    let bad: Vec<&Ran> = ran.iter().filter(|r| r.verdict.outcome.is_bad()).collect();
    if !bad.is_empty() {
        println!("\nbad cases:");
        for r in &bad {
            println!(
                "  {} [{}] {}",
                r.key,
                r.verdict.outcome.name(),
                r.verdict.note
            );
        }
    }
    // A volume the exact integral refused was read off a tessellation, and
    // its case judged with the chord in its slack.
    let fell: Vec<&str> = ran
        .iter()
        .filter(|r| r.verdict.tessellated)
        .map(|r| r.key.as_str())
        .collect();
    if !fell.is_empty() {
        println!(
            "\nmeasured on a tessellation ({} cases): {}",
            fell.len(),
            fell.join(" ")
        );
    }
    // Refusals by reason, numbers stripped, so a family reads as one line.
    let mut reasons: BTreeMap<(String, String), usize> = BTreeMap::new();
    for r in ran.iter().filter(|r| r.verdict.outcome == Outcome::Refused) {
        let reason: String = r
            .verdict
            .note
            .chars()
            .map(|c| {
                if c.is_ascii_digit() || c == '.' || c == '-' {
                    '#'
                } else {
                    c
                }
            })
            .collect::<String>()
            .split("##")
            .collect::<Vec<_>>()
            .join("#");
        let reason: String = reason.chars().take(110).collect();
        *reasons.entry((r.scenario.clone(), reason)).or_insert(0) += 1;
    }
    if !reasons.is_empty() {
        let mut by_count: Vec<_> = reasons.into_iter().collect();
        by_count.sort_by_key(|entry| std::cmp::Reverse(entry.1));
        println!("\nrefusals:");
        for ((scenario, reason), n) in by_count.iter().take(12) {
            println!("  {n:>4} {scenario}: {reason}");
        }
    }
    let mut slow: Vec<&Ran> = ran.iter().collect();
    slow.sort_by(|a, b| b.seconds.total_cmp(&a.seconds));
    println!("\nslowest:");
    for r in slow.iter().take(5) {
        println!("  {} {:.2}s", r.key, r.seconds);
    }
}

/// The baseline comparison: fewer `ok` or more bad cases in any scenario
/// is a regression.
fn check(tally: &Tally, baseline: &Tally) -> bool {
    let mut fine = true;
    for (scenario, want) in baseline {
        let Some(got) = tally.get(scenario) else {
            println!("regression: scenario {scenario} did not run");
            fine = false;
            continue;
        };
        let count =
            |row: &BTreeMap<&str, usize>, o: Outcome| row.get(o.name()).copied().unwrap_or(0);
        let bad = |row: &BTreeMap<&str, usize>| {
            Outcome::ALL
                .iter()
                .filter(|o| o.is_bad())
                .map(|o| count(row, *o))
                .sum::<usize>()
        };
        if count(got, Outcome::Ok) < count(want, Outcome::Ok) {
            println!(
                "regression: {scenario} ok {} below the baseline's {}",
                count(got, Outcome::Ok),
                count(want, Outcome::Ok)
            );
            fine = false;
        }
        if bad(got) > bad(want) {
            println!(
                "regression: {scenario} has {} bad cases against the baseline's {}",
                bad(got),
                bad(want)
            );
            fine = false;
        }
        if count(got, Outcome::Ok) > count(want, Outcome::Ok) || bad(got) < bad(want) {
            println!("improved: {scenario}; record it with --write");
        }
    }
    fine
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let value = |flag: &str| {
        args.iter()
            .position(|a| a == flag)
            .and_then(|i| args.get(i + 1))
            .cloned()
    };
    let seed: u64 = value("--seed").map_or(1, |s| s.parse().expect("--seed takes a number"));
    let per_part: usize =
        value("--cases").map_or(20, |s| s.parse().expect("--cases takes a number"));
    let corpus = !args.iter().any(|a| a == "--no-corpus");
    let threads = value("--threads").map_or_else(
        || std::thread::available_parallelism().map_or(4, std::num::NonZero::get),
        |s| s.parse().expect("--threads takes a number"),
    );

    // Panics are outcomes here, reported in the tally, not printed as they
    // happen.
    if !args.iter().any(|a| a == "--case") {
        std::panic::set_hook(Box::new(|_| {}));
    }

    eprintln!("building parts and cases (seed {seed}, {per_part} per part)");
    let all = cases(seed, per_part, corpus);

    if let Some(key) = value("--case") {
        let Some(case) = all.iter().find(|c| c.key == key) else {
            eprintln!("no case {key}");
            std::process::exit(2);
        };
        let start = Instant::now();
        let verdict = (case.run)();
        println!(
            "{key}: {} {} ({:.2}s)",
            verdict.outcome.name(),
            verdict.note,
            start.elapsed().as_secs_f64()
        );
        return;
    }

    let only = value("--scenario");
    let part = value("--part");
    let chosen: Vec<Case> = all
        .into_iter()
        .filter(|c| {
            only.as_deref()
                .is_none_or(|s| c.key.starts_with(&format!("{s}/")))
        })
        .filter(|c| {
            part.as_deref()
                .is_none_or(|p| c.key.split('/').nth(1) == Some(p))
        })
        .collect();
    eprintln!("running {} cases on {threads} threads", chosen.len());
    let start = Instant::now();
    let ran = run_all(&chosen, threads);
    let tally = tally(&ran);
    print_report(&ran, &tally);
    if args.iter().any(|a| a == "--times") {
        println!("\ntimes:");
        for r in &ran {
            println!("  {} {:.2}s {}", r.key, r.seconds, r.verdict.outcome.name());
        }
    }
    println!(
        "\n{} cases in {:.1}s",
        ran.len(),
        start.elapsed().as_secs_f64()
    );

    if let Some(path) = value("--write") {
        std::fs::write(&path, tally_json(&tally)).expect("writing the baseline");
        println!("baseline written to {path}");
    }
    if let Some(path) = value("--check") {
        let text = std::fs::read_to_string(&path).expect("reading the baseline");
        if !check(&tally, &read_tally(&text)) {
            std::process::exit(1);
        }
        println!("no regression against {path}");
    }
}
