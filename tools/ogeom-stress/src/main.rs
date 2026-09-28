//! Randomized robustness runs over the kernel: many generated cases, each
//! judged by an oracle, the tally watched against a recorded baseline.
//!
//! A single repro proves one case; a kernel is trusted on the rate at which
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
//! Refusals are honest and tracked; the last three are the bugs.
//!
//! `ogeom-stress` runs every scenario and prints the tally. `--check
//! <baseline.json>` fails when a scenario has fewer `ok` or more bad cases
//! than the baseline records; `--write <path>` records one. `--case
//! drill/box/17` replays one case and prints what happened.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    reason = "a reporting tool"
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
}

impl Verdict {
    fn ok() -> Self {
        Self {
            outcome: Outcome::Ok,
            note: String::new(),
        }
    }

    fn of(outcome: Outcome, note: impl Into<String>) -> Self {
        Self {
            outcome,
            note: note.into(),
        }
    }
}

/// An operation's result judged for validity: an error is a refusal, an
/// invalid solid is `invalid`, and a valid one comes back with its volume.
fn judged(
    model: &Model,
    what: &str,
    result: ogeom::core::OgeomResult<ogeom::algo::Built>,
) -> Result<(Shape, f64), Verdict> {
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
    let v = volume(model, &built.shape)
        .map_err(|e| Verdict::of(Outcome::Refused, format!("{what}: volume: {e}")))?;
    Ok((built.shape, v))
}

fn volume(model: &Model, shape: &Shape) -> ogeom::core::OgeomResult<f64> {
    let deflection = Deflection::with_chord(2e-3)?;
    Ok(ogeom::algo::volume_properties(model, shape, deflection, T)?.mass)
}

// --- parts -----------------------------------------------------------------

/// A solid the scenarios work on, built once and copied per case.
#[derive(Clone)]
struct Part {
    name: String,
    model: Model,
    shape: Shape,
    volume: f64,
    low: Point,
    high: Point,
    vertices: Vec<Point>,
    /// Edges with a curve between two different faces: what a fillet can
    /// round.
    creases: Vec<Shape>,
}

impl Part {
    fn new(name: &str, model: Model, shape: Shape) -> Option<Self> {
        let volume = volume(&model, &shape).ok()?;
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
            volume,
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

/// Edges with a curve between two different faces. A seam (one face on
/// both sides) or a pole has no corner to round.
fn creases(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore_unique(model, shape, ShapeType::Edge)
        .unwrap_or_default()
        .into_iter()
        .filter(|e| {
            let has_curve = model
                .node(e)
                .and_then(|n| n.data().as_edge())
                .is_some_and(|d| d.curve3d().is_some() && !d.degenerate);
            let faces = ogeom::topo::ancestors_of(model, shape, e, ShapeType::Face)
                .map(|fs| {
                    let mut nodes: Vec<_> = fs.iter().map(Shape::node).collect();
                    nodes.sort_unstable();
                    nodes.dedup();
                    nodes.len()
                })
                .unwrap_or(0);
            has_curve && faces == 2
        })
        .collect()
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
    // Every part here drills in seconds. Half of nist_ftc_06's drills are
    // refused today (the arrangement leaves a face no piece, or the kept
    // pieces do not close), and the baseline says so, so a fix shows as an
    // improvement.
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
/// vertex or is tangent to a plane through one.
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
    let placed = format!(
        "drill r {radius} from ({}, {}, {}) along ({}, {}, {}) length {length}",
        start.x, start.y, start.z, axis.x, axis.y, axis.z
    );
    let verdict = (|| {
        let drill_volume = core::f64::consts::PI * radius * radius * length;
        let cut = ogeom::boolean::cut(&mut model, &part.shape, &drill, T);
        let (_, v_cut) = match judged(&model, "cut", cut) {
            Ok(r) => r,
            Err(v) => return v,
        };
        let common = ogeom::boolean::common(&mut model, &part.shape, &drill, T);
        let v_common = match common {
            // Nothing in common is an empty answer, not a refusal.
            Err(e) if v_cut >= part.volume * (1.0 - 1e-9) => {
                let _ = e;
                0.0
            }
            other => match judged(&model, "common", other) {
                Ok((_, v)) => v,
                Err(v) => return v,
            },
        };
        let slack = 2e-3 * part.volume.max(drill_volume.min(part.volume));
        let miss = (v_cut + v_common - part.volume).abs();
        if miss > slack {
            return Verdict::of(
                Outcome::Wrong,
                format!(
                    "cut {v_cut:.6} + common {v_common:.6} = {:.6}, part {:.6}",
                    v_cut + v_common,
                    part.volume
                ),
            );
        }
        Verdict::ok()
    })();
    if verdict.outcome.is_bad() {
        Verdict::of(verdict.outcome, format!("{}; {placed}", verdict.note))
    } else {
        verdict
    }
}

/// A primitive of a random kind and size, placed with its position and
/// size snapped to a coarse grid half the time, so pairs meet face on face.
fn primitive(model: &mut Model, rng: &mut Rng) -> Option<(Shape, f64, String)> {
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
    let axis = *rng.pick(&[Direction::X, Direction::Y, Direction::Z]);
    let x = if axis == Direction::X {
        Direction::Y
    } else {
        Direction::X
    };
    let frame = frame_at(at, axis, x);
    let size = |rng: &mut Rng| {
        let s = snap(rng, 2.0, 8.0);
        s.max(1.0)
    };
    let kind = rng.below(5);
    let (built, name) = match kind {
        0 => {
            let (a, b, c) = (size(rng), size(rng), size(rng));
            (ogeom::algo::make_box(model, frame, (a, b, c), T), "box")
        }
        1 => {
            let (r, h) = (size(rng) * 0.5, size(rng));
            (
                ogeom::algo::make_cylinder(model, frame, r, h, T),
                "cylinder",
            )
        }
        2 => (
            ogeom::algo::make_sphere(model, frame, size(rng) * 0.5, T),
            "sphere",
        ),
        3 => {
            let (r0, h) = (size(rng) * 0.5, size(rng));
            let r1 = r0 * rng.range(0.2, 0.9);
            (ogeom::algo::make_cone(model, frame, r0, r1, h, T), "cone")
        }
        _ => {
            let r = size(rng) * 0.5;
            let tube = (r * rng.range(0.15, 0.45)).max(0.3);
            (ogeom::algo::make_torus(model, frame, r, tube, T), "torus")
        }
    };
    let shape = built.ok()?.shape;
    let v = volume(model, &shape).ok()?;
    Some((shape, v, name.to_string()))
}

/// Two primitives, and the boolean identities: the fuse and the common
/// add up to the two, and the cut is the first less the common.
fn pair_case(rng: &mut Rng) -> Verdict {
    let mut model = Model::new();
    let Some((a, va, ka)) = primitive(&mut model, rng) else {
        return Verdict::of(Outcome::Refused, "primitive");
    };
    let Some((b, vb, kb)) = primitive(&mut model, rng) else {
        return Verdict::of(Outcome::Refused, "primitive");
    };
    let what = format!("{ka} with {kb}");
    let common = ogeom::boolean::common(&mut model, &a, &b, T);
    let v_common = match common {
        Err(_) => None,
        other => match judged(&model, &format!("{what}: common"), other) {
            Ok((_, v)) => Some(v),
            Err(v) => return v,
        },
    };
    let fuse = ogeom::boolean::fuse(&mut model, &a, &b, T);
    let (_, v_fuse) = match judged(&model, &format!("{what}: fuse"), fuse) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let cut = ogeom::boolean::cut(&mut model, &a, &b, T);
    let (_, v_cut) = match judged(&model, &format!("{what}: cut"), cut) {
        Ok(r) => r,
        Err(v) => return v,
    };
    // A refused common is only an empty one when the fuse says the two
    // do not overlap.
    let v_common = match v_common {
        Some(v) => v,
        None if (v_fuse - va - vb).abs() <= 2e-3 * (va + vb) => 0.0,
        None => {
            return Verdict::of(Outcome::Refused, format!("{what}: common refused"));
        }
    };
    let slack = 2e-3 * (va + vb);
    if (v_fuse + v_common - va - vb).abs() > slack {
        return Verdict::of(
            Outcome::Wrong,
            format!(
                "{what}: fuse {v_fuse:.6} + common {v_common:.6} against {:.6}",
                va + vb
            ),
        );
    }
    if (v_cut + v_common - va).abs() > slack {
        return Verdict::of(
            Outcome::Wrong,
            format!("{what}: cut {v_cut:.6} + common {v_common:.6} against {va:.6}"),
        );
    }
    Verdict::ok()
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
    let (_, v) = match judged(&model, "fillet", result) {
        Ok(r) => r,
        Err(v) => return v,
    };
    let change = (v - part.volume).abs() / part.volume;
    if change > 0.25 {
        return Verdict::of(
            Outcome::Wrong,
            format!(
                "fillet r {radius:.4} moved {:.1}% of the part",
                change * 100.0
            ),
        );
    }
    Verdict::ok()
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
        // Corpus parts are large; they take a third of the drills.
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
    for i in 0..per_part * 6 {
        let key = format!("pairs/primitives/{i}");
        let s = case_seed(seed, &key);
        out.push(Case {
            key,
            run: Box::new(move || pair_case(&mut Rng::new(s))),
        });
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
