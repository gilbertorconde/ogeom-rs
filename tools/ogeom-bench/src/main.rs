//! Benchmarks over the kernel's hot paths.
//!
//! The harness is its own, with no dependency: wall clock, one warm-up run,
//! then samples until at least [`TARGET`] of them is timed (never fewer
//! than [`LEAST`], and bounded by [`WALL`] of wall clock), reported as the
//! minimum, the median and the median absolute deviation (MAD). Absolute
//! times move with the machine, so every run also times a fixed
//! single-threaded arithmetic spin (the *calibration*) and reports each
//! benchmark's minimum as a multiple of the spin's minimum. Ratios travel
//! between machines; milliseconds do not. The minimum is what is compared
//! because load on a shared machine only ever adds time.
//!
//! Usage: `ogeom-bench [--threads N] [--filter TEXT] [--check BASELINE]`.
//!
//! - `--threads N` sets the thread count of the kernel's parallel stages
//!   (`ogeom_core::parallel::set_threads`). Without it, the
//!   `OGEOM_THREADS` environment variable or the machine's parallelism
//!   decides, except under `--check`, which runs at the thread count the
//!   baseline was recorded at.
//! - `--filter TEXT` runs only the benchmarks whose name contains `TEXT`.
//! - `--check BASELINE` compares the calibrated ratios against a recorded
//!   baseline and reports the drift, informationally: performance is
//!   watched here, not gated, because a loaded CI box would turn a real
//!   gate into a coin flip.
//!
//! Without `--check` the run prints a table to stderr and the baseline JSON
//! to stdout. See `README.md` beside this crate for what the ratios mean
//! across core counts.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "a reporting tool")]
#![allow(
    clippy::disallowed_methods,
    reason = "a native benchmark: it times with the std clock"
)]

use std::time::{Duration, Instant};

use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve3d as _, PlaneSurface, Surface as _};
use ogeom::math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, Triangulation, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// Samples are taken until at least this much of them is timed.
const TARGET: Duration = Duration::from_millis(50);
/// Never fewer samples than this.
const LEAST: usize = 3;
/// Past this much wall clock, setup included, sampling stops once
/// [`LEAST`] samples are in.
const WALL: Duration = Duration::from_secs(3);
/// Never more samples than this, however short the benchmark.
const MOST: usize = 5000;

/// What one benchmark measured, in seconds.
struct Stats {
    min: f64,
    median: f64,
    mad: f64,
    samples: usize,
}

/// Time `run` on a fresh `setup()` each sample, the setup untimed: one
/// warm-up, then samples as the constants above bound them.
fn sample<S>(mut setup: impl FnMut() -> S, mut run: impl FnMut(S)) -> Stats {
    run(setup());
    let began = Instant::now();
    let mut times = Vec::new();
    let mut timed = 0.0;
    while times.len() < MOST
        && (times.len() < LEAST || (timed < TARGET.as_secs_f64() && began.elapsed() < WALL))
    {
        let state = setup();
        let start = Instant::now();
        run(state);
        let seconds = start.elapsed().as_secs_f64();
        timed += seconds;
        times.push(seconds);
    }
    times.sort_by(f64::total_cmp);
    let median = middle(&times);
    let mut deviations: Vec<f64> = times.iter().map(|t| (t - median).abs()).collect();
    deviations.sort_by(f64::total_cmp);
    Stats {
        min: times[0],
        median,
        mad: middle(&deviations),
        samples: times.len(),
    }
}

/// Time `run` with nothing to set up per sample.
fn time(mut run: impl FnMut()) -> Stats {
    sample(|| (), |()| run())
}

fn middle(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        0.5 * (sorted[n / 2 - 1] + sorted[n / 2])
    }
}

/// A fixed single-threaded arithmetic spin whose time stands for this
/// machine's speed.
fn calibration() -> Stats {
    time(|| {
        let mut acc = 0.0f64;
        for i in 0..4_000_000u64 {
            #[allow(clippy::cast_precision_loss)]
            let x = (i as f64).mul_add(1.000_000_1, acc).sin();
            acc = x * 0.999;
        }
        std::hint::black_box(acc);
    })
}

fn corpus(name: &str) -> Option<String> {
    let path = format!("{}/../../tests/corpus/{name}", env!("CARGO_MANIFEST_DIR"));
    std::fs::read_to_string(path).ok()
}

/// The first solid of a corpus STEP part, with its document.
fn corpus_part(name: &str) -> Option<(ogeom::doc::Document, Shape)> {
    let import = ogeom::io::read_step(&corpus(name)?, T).ok()?;
    let solid = import.solids.first().cloned()?;
    Some((import.document, solid))
}

/// The small corpus part the slower whole-part benchmarks run on: an
/// imported b-rep with spline faces, holes and blends.
const PART: &str = "nist_ftc_11_asme1_rb.stp";

/// The largest corpus part, for the whole-part benchmarks fast enough to
/// run on it, where the smallest part takes a few milliseconds.
const LARGE_PART: &str = "nist_ctc_02_asme1_rc.stp";

type Bench = (&'static str, Box<dyn Fn() -> Option<Stats>>);

fn benchmarks() -> Vec<Bench> {
    vec![
        ("construct_boxes", Box::new(construct_boxes)),
        ("traverse_box", Box::new(traverse_box)),
        ("tessellate_torus", Box::new(tessellate_torus)),
        ("quartic_torus", Box::new(quartic_torus)),
        ("intersect_line_torus", Box::new(intersect_line_torus)),
        ("boolean_drill", Box::new(boolean_drill)),
        ("boolean_many_faces", Box::new(boolean_many_faces)),
        ("boolean_local", Box::new(boolean_local)),
        ("boolean_marched", Box::new(boolean_marched)),
        ("fillet_block", Box::new(fillet_block)),
        ("fillet_box_all", Box::new(fillet_box_all)),
        ("fillet_marched", Box::new(fillet_marched)),
        ("tessellate_part", Box::new(tessellate_part)),
        ("tessellate_part_stored", Box::new(tessellate_part_stored)),
        ("tessellate_open_sheet", Box::new(tessellate_open_sheet)),
        ("mass_part", Box::new(mass_part)),
        ("check_part", Box::new(check_part)),
        ("fix_shape_part", Box::new(fix_shape_part)),
        ("sew_shell", Box::new(sew_shell)),
        ("mesh_to_solid", Box::new(mesh_to_solid)),
        ("mesh_to_solid_mid", Box::new(mesh_to_solid_mid)),
        ("hlr_exact", Box::new(hlr_exact)),
        ("hlr_exact_mid", Box::new(hlr_exact_mid)),
        ("hlr_mesh", Box::new(hlr_mesh)),
        ("thick_spline", Box::new(thick_spline)),
        ("project_spline", Box::new(project_spline)),
        ("offset_part", Box::new(offset_part)),
        ("shell_part", Box::new(shell_part)),
        ("pipe_circles", Box::new(pipe_circles)),
        ("pipe_guided", Box::new(pipe_guided)),
        ("import_ftc11", Box::new(import_ftc11)),
        ("import_ctc02", Box::new(import_ctc02)),
        ("step_write", Box::new(step_write)),
        ("iges_read", Box::new(iges_read)),
        ("weld_1m", Box::new(weld_1m)),
        ("stl_read_1m", Box::new(stl_read_1m)),
    ]
}

/// Construction: a thousand boxes into one model.
fn construct_boxes() -> Option<Stats> {
    Some(time(|| {
        let mut model = Model::new();
        for i in 0..1000 {
            #[allow(clippy::cast_precision_loss)]
            let dx = i as f64;
            let frame =
                Frame::new(Point::new(dx, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
            ogeom::algo::make_box(&mut model, frame, (1.0, 1.0, 1.0), T).unwrap();
        }
        std::hint::black_box(&model);
    }))
}

/// Traversal: exploring a box's sub-shapes, ten thousand times.
fn traverse_box() -> Option<Stats> {
    let mut model = Model::new();
    let solid = ogeom::algo::make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T)
        .unwrap()
        .shape;
    Some(time(|| {
        for _ in 0..10_000 {
            let faces = explore_unique(&model, &solid, ShapeType::Face).unwrap();
            std::hint::black_box(faces.len());
        }
    }))
}

/// Tessellation: a torus at the default deflection.
fn tessellate_torus() -> Option<Stats> {
    Some(time(|| {
        let mut model = Model::new();
        let solid = ogeom::algo::make_torus(&mut model, Frame::WORLD, 20.0, 5.0, T)
            .unwrap()
            .shape;
        let done = ogeom::mesh::tessellate(&mut model, &solid, Deflection::default(), T).unwrap();
        std::hint::black_box(done.triangles);
    }))
}

/// Lines against a torus of major radius 20 and minor radius 5 at the
/// origin, in its own frame: a deterministic spread through the hole,
/// across the tube, grazing its top and in any direction, a thousand in all.
fn torus_lines() -> Vec<(Point, Vector)> {
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        #[allow(clippy::cast_precision_loss)]
        let unit = (state >> 11) as f64 / (1_u64 << 53) as f64;
        unit.mul_add(2.0, -1.0)
    };
    let mut lines = Vec::with_capacity(1000);
    for i in 0..1000 {
        let line = match i % 4 {
            // Through the hole, steeply.
            0 => (
                Point::new(next() * 10.0, next() * 10.0, -40.0),
                Vector::new(next() * 0.3, next() * 0.3, 1.0),
            ),
            // Across the whole torus, in its plane or near it.
            1 => (
                Point::new(-40.0, next() * 24.0, next() * 4.0),
                Vector::new(1.0, next() * 0.1, next() * 0.1),
            ),
            // Along the tube's top, grazing it.
            2 => (
                Point::new(-40.0, next() * 20.0, 5.0),
                Vector::new(1.0, 0.0, 0.0),
            ),
            // Anywhere.
            _ => (
                Point::new(next() * 40.0, next() * 40.0, next() * 40.0),
                Vector::new(next(), next(), next()),
            ),
        };
        lines.push(line);
    }
    lines
}

/// The line-torus quartic alone: the polynomial roots of a thousand lines
/// against a torus, the coefficients formed as the intersection forms them.
fn quartic_torus() -> Option<Stats> {
    let (big, small) = (20.0 / 25.0, 5.0 / 25.0);
    let quartics: Vec<[f64; 5]> = torus_lines()
        .into_iter()
        .map(|(from, along)| {
            let d = along / along.magnitude();
            let m = (from - Point::ORIGIN) / 25.0;
            let m = m - d * m.dot(d);
            let a = d.dot(d);
            let b = 2.0 * m.dot(d);
            let c = m.dot(m) + big * big - small * small;
            let p = d.x.mul_add(d.x, d.y * d.y);
            let q = 2.0 * m.x.mul_add(d.x, m.y * d.y);
            let s = m.x.mul_add(m.x, m.y * m.y);
            let four = 4.0 * big * big;
            [
                c.mul_add(c, -four * s),
                2.0f64.mul_add(b * c, -four * q),
                b.mul_add(b, 2.0 * a * c) - four * p,
                2.0 * a * b,
                a * a,
            ]
        })
        .collect();
    Some(time(|| {
        for quartic in &quartics {
            let found = ogeom::math::solve::roots(std::hint::black_box(quartic), 1e-9).unwrap();
            std::hint::black_box(found);
        }
    }))
}

/// Line-torus intersection: a thousand lines against a torus in a tilted
/// frame, roots and polish together.
fn intersect_line_torus() -> Option<Stats> {
    use ogeom::geom::{Curve, LineCurve, SurfaceGeometry, TorusSurface};
    use ogeom::intersect::{CurveSurfaceOptions, intersect_curve_surface};
    let frame = Frame::new(
        Point::new(3.0, -2.0, 1.0),
        Direction::new(Vector::new(0.1, 0.2, 1.0), T).unwrap(),
        Direction::new(Vector::new(1.0, 0.0, -0.1), T).unwrap(),
        T,
    )
    .unwrap();
    let torus: SurfaceGeometry =
        TorusSurface::new(ogeom::math::Torus::new(frame, 20.0, 5.0, T).unwrap()).into();
    let lines: Vec<Curve> = torus_lines()
        .into_iter()
        .map(|(from, along)| {
            let to = from + along * (80.0 / along.magnitude());
            LineCurve::segment(frame.to_world(from), frame.to_world(to), T)
                .unwrap()
                .into()
        })
        .collect();
    let options = CurveSurfaceOptions::default();
    Some(time(|| {
        for line in &lines {
            let found = intersect_curve_surface(line, &torus, options, T).unwrap();
            std::hint::black_box(found);
        }
    }))
}

/// The boolean: the drilled box.
fn boolean_drill() -> Option<Stats> {
    Some(time(|| {
        let mut model = Model::new();
        let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        let frame =
            Frame::new(Point::new(10.0, 10.0, -1.0), Direction::Z, Direction::X, T).unwrap();
        let drill = ogeom::algo::make_cylinder(&mut model, frame, 3.0, 12.0, T)
            .unwrap()
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &block, &drill, T).unwrap();
        std::hint::black_box(&cut.shape);
    }))
}

/// The boolean with many faces in play: a block drilled four times. The
/// classifier is asked once per face piece and asks the *other* solid's
/// every face, so this is where that product shows.
fn boolean_many_faces() -> Option<Stats> {
    Some(time(|| {
        let mut model = Model::new();
        let mut block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
            .unwrap()
            .shape;
        for (x, y) in [(5.0, 5.0), (15.0, 5.0), (5.0, 15.0), (15.0, 15.0)] {
            let frame = Frame::new(Point::new(x, y, -1.0), Direction::Z, Direction::X, T).unwrap();
            let drill = ogeom::algo::make_cylinder(&mut model, frame, 2.0, 12.0, T)
                .unwrap()
                .shape;
            block = ogeom::boolean::cut(&mut model, &block, &drill, T)
                .unwrap()
                .shape;
        }
        std::hint::black_box(&block);
    }))
}

/// A local edit to a large solid: a short drill into one side of a plate
/// with a few hundred holes, and a slot across its top between two rows of
/// them. The tool touches a face or two; what the boolean costs past that
/// is what it spends on the faces the tool never reaches.
fn boolean_local() -> Option<Stats> {
    let mut model = Model::new();
    let (rows, pitch) = (15_u32, 6.0);
    let size = 8.0 + pitch * f64::from(rows);
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (size, size, 5.0), T)
        .unwrap()
        .shape;
    let mut pins = Vec::new();
    for i in 0..rows {
        for j in 0..rows {
            let at = Point::new(8.0 + pitch * f64::from(i), 8.0 + pitch * f64::from(j), -1.0);
            let frame = Frame::new(at, Direction::Z, Direction::X, T).unwrap();
            pins.push(
                ogeom::algo::make_cylinder(&mut model, frame, 1.0, 7.0, T)
                    .unwrap()
                    .shape,
            );
        }
    }
    let pins = model.add_compound(&pins).unwrap();
    let plate = ogeom::boolean::cut(&mut model, &block, &pins, T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(-1.0, 5.0, 2.5), Direction::X, Direction::Y, T).unwrap();
    let drill = ogeom::algo::make_cylinder(&mut model, frame, 1.0, 4.0, T)
        .unwrap()
        .shape;
    let frame = Frame::new(Point::new(-1.0, 10.5, 2.5), Direction::Z, Direction::X, T).unwrap();
    let slot = ogeom::algo::make_box(&mut model, frame, (size + 2.0, 2.0, 3.5), T)
        .unwrap()
        .shape;
    Some(sample(
        || model.clone(),
        |mut model| {
            for tool in [&drill, &slot] {
                let cut = ogeom::boolean::cut(&mut model, &plate, tool, T).unwrap();
                std::hint::black_box(&cut.shape);
            }
        },
    ))
}

/// The boolean with no closed form: crossed cylinders, whose sections only
/// the marcher can trace. `boolean_drill` is analytic end to end and never
/// reaches that machinery.
fn boolean_marched() -> Option<Stats> {
    Some(time(|| {
        let mut model = Model::new();
        let upright = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 1.0, 4.0, T)
            .unwrap()
            .shape;
        let frame = Frame::new(Point::new(-2.0, 0.0, 2.0), Direction::X, Direction::Y, T).unwrap();
        let across = ogeom::algo::make_cylinder(&mut model, frame, 0.6, 4.0, T)
            .unwrap()
            .shape;
        let both = ogeom::boolean::fuse(&mut model, &upright, &across, T).unwrap();
        std::hint::black_box(&both.shape);
    }))
}

/// A constant-radius fillet on one edge of a box: the blend machinery.
fn fillet_block() -> Option<Stats> {
    Some(time(|| {
        let mut model = Model::new();
        let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
        let rolled = ogeom::fillet::fillet_edges(&mut model, &block, &edges[..1], 1.0, T).unwrap();
        std::hint::black_box(&rolled.shape);
    }))
}

/// All twelve edges of a box filleted at once: the corners where three
/// blends meet.
fn fillet_box_all() -> Option<Stats> {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    Some(sample(
        || model.clone(),
        |mut model| {
            let rolled = ogeom::fillet::fillet_edges(&mut model, &block, &edges, 1.0, T).unwrap();
            std::hint::black_box(&rolled.shape);
        },
    ))
}

/// A fillet no closed form speaks: the elliptical seat of a post leaning
/// twenty degrees out of a slab, rolled by the marched blend.
fn fillet_marched() -> Option<Stats> {
    let mut model = Model::new();
    let slab = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 2.0), T)
        .unwrap()
        .shape;
    let lean = 20.0_f64.to_radians();
    let axis = Vector::new(lean.sin(), 0.0, lean.cos());
    let frame = Frame::new(
        Point::new(10.0, 10.0, -1.0),
        Direction::new(axis, T).unwrap(),
        Direction::from_cross(axis, Vector::Y, T).unwrap(),
        T,
    )
    .unwrap();
    let post = ogeom::algo::make_cylinder(&mut model, frame, 3.0, 10.0, T)
        .unwrap()
        .shape;
    let joined = ogeom::boolean::fuse(&mut model, &slab, &post, T)
        .unwrap()
        .shape;
    // The seat is the elliptical edge on the slab's top; the post pierces
    // the slab, so its bottom carries a second ellipse.
    let seat = explore_unique(&model, &joined, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let elliptical = model
                .node(e)
                .and_then(|n| n.data().as_edge())
                .and_then(|d| d.curve3d())
                .and_then(|r| match r {
                    ogeom::topo::EdgeRepr::Curve3d { curve, .. } => model.geometry().curve(*curve),
                    _ => None,
                })
                .is_some_and(|c| matches!(c, ogeom::geom::Curve::Ellipse(_)));
            elliptical
                && ogeom::algo::edge_vertices(&model, e)
                    .unwrap()
                    .and_then(|(a, _)| model.node(&a)?.data().as_vertex().map(|d| d.point))
                    .is_some_and(|p| (p.z - 2.0).abs() < 1e-6)
        })?;
    Some(sample(
        || model.clone(),
        |mut model| {
            let rolled = ogeom::fillet::fillet_edge(&mut model, &joined, &seat, 1.0, T).unwrap();
            std::hint::black_box(&rolled.shape);
        },
    ))
}

/// Whole-shape tessellation of a real imported part, whose faces are not
/// all analytic. A primitive will not do: a torus answers from its closed
/// form.
fn tessellate_part() -> Option<Stats> {
    let (document, solid) = corpus_part(PART)?;
    Some(time(|| {
        let mesh = ogeom::mesh::triangulate(
            document.model(),
            &solid,
            Deflection::with_chord(1e-2).unwrap(),
            T,
        )
        .unwrap();
        std::hint::black_box(mesh.triangles.len());
    }))
}

/// The largest corpus part's tessellation stored on its model: each face's
/// mesh, a polyline on each edge and each edge's path through its faces'
/// meshes, on a fresh copy of the model each time.
fn tessellate_part_stored() -> Option<Stats> {
    let (document, solid) = corpus_part(LARGE_PART)?;
    let model = document.model().clone();
    Some(sample(
        || model.clone(),
        |mut model| {
            let done =
                ogeom::mesh::tessellate(&mut model, &solid, Deflection::default(), T).unwrap();
            std::hint::black_box(done.triangles);
        },
    ))
}

/// An open sheet with a long border: a cylinder's side alone, drawn fine
/// enough that each rim is some thousands of segments. Every rim segment
/// is a border edge for the whole-shape repair passes to look along.
fn tessellate_open_sheet() -> Option<Stats> {
    let mut model = Model::new();
    let solid = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 100.0, 50.0, T)
        .unwrap()
        .shape;
    let side = explore_unique(&model, &solid, ShapeType::Face)
        .unwrap()
        .into_iter()
        .max_by_key(|f| {
            explore_unique(&model, f, ShapeType::Edge)
                .map(|e| e.len())
                .unwrap_or(0)
        })?;
    Some(time(|| {
        let mesh =
            ogeom::mesh::triangulate(&model, &side, Deflection::with_chord(1e-3).unwrap(), T)
                .unwrap();
        std::hint::black_box(mesh.triangles.len());
    }))
}

/// Mass properties of the largest corpus part on the exact path: each face
/// integrated on its own surface, at the default deflection.
fn mass_part() -> Option<Stats> {
    let (document, solid) = corpus_part(LARGE_PART)?;
    Some(time(|| {
        let props =
            ogeom::algo::volume_properties(document.model(), &solid, Deflection::default(), T)
                .unwrap();
        std::hint::black_box(props.mass);
    }))
}

/// The validity check over the largest corpus part.
fn check_part() -> Option<Stats> {
    let (document, solid) = corpus_part(LARGE_PART)?;
    Some(time(|| {
        let diagnosis = ogeom::algo::check(document.model(), &solid, T).unwrap();
        std::hint::black_box(diagnosis.is_valid());
    }))
}

/// Shape healing over the largest corpus part, on a fresh copy of its model each
/// time.
fn fix_shape_part() -> Option<Stats> {
    let (document, solid) = corpus_part(LARGE_PART)?;
    let model = document.model().clone();
    Some(sample(
        || model.clone(),
        |mut model| {
            let fixed = ogeom::heal::fix_shape(&mut model, &solid, T).unwrap();
            std::hint::black_box(&fixed.shape);
        },
    ))
}

/// Sewing the faces of a closed shell: a cube of side 16 whose every side
/// is a 16 by 16 grid of unit squares, each built as its own face with its
/// own edges and vertices, so 1536 faces with 3072 edge pairs to find.
fn sew_shell() -> Option<Stats> {
    const N: i32 = 16;
    let mut model = Model::new();
    let mut faces = Vec::new();
    let half = f64::from(N) / 2.0;
    // Each side's outward normal, and an in-plane axis along the grid.
    for (normal, a) in [
        (Vector::X, Vector::Y),
        (-Vector::X, Vector::Y),
        (Vector::Y, Vector::Z),
        (-Vector::Y, Vector::Z),
        (Vector::Z, Vector::X),
        (-Vector::Z, Vector::X),
    ] {
        let n = Direction::new(normal, T).unwrap();
        let b = normal.cross(a);
        let centre = Point::ORIGIN + normal * half;
        let plane = Plane::through(centre, n);
        for i in 0..N {
            for j in 0..N {
                let corner = |di: i32, dj: i32| {
                    centre + a * (f64::from(i + di) - half) + b * (f64::from(j + dj) - half)
                };
                let ring = [corner(0, 0), corner(1, 0), corner(1, 1), corner(0, 1)];
                let wire = ogeom::algo::make_polygon(&mut model, &ring, true, T)
                    .unwrap()
                    .shape;
                faces.push(
                    ogeom::algo::make_face(&mut model, PlaneSurface::new(plane).into(), &[wire], T)
                        .unwrap()
                        .shape,
                );
            }
        }
    }
    Some(sample(
        || model.clone(),
        |mut model| {
            let sewn = ogeom::algo::sew(&mut model, &faces, T).unwrap();
            assert!(sewn.free_edges.is_empty());
            std::hint::black_box(&sewn.shells);
        },
    ))
}

/// A b-rep from a mesh: an imported part tessellated, then rebuilt with
/// its planes, cylinders and other surfaces recognized.
fn mesh_to_solid() -> Option<Stats> {
    converted(PART)
}

/// The mid-size corpus part the conversion from a mesh is timed on: some
/// eleven thousand triangles rebuilt as a few hundred faces, over several
/// builds. The largest part takes minutes at one thread.
const MESH_MID_PART: &str = "nist_ftc_08_asme1_rc.stp";

/// A b-rep from a mesh of a part with many faces, where the checks run
/// over every built face on each build dominate.
fn mesh_to_solid_mid() -> Option<Stats> {
    converted(MESH_MID_PART)
}

/// A corpus part's first solid tessellated, and the time to convert the
/// tessellation back.
fn converted(part: &str) -> Option<Stats> {
    let (document, solid) = corpus_part(part)?;
    let mesh = ogeom::mesh::triangulate(document.model(), &solid, Deflection::default(), T).ok()?;
    let options = ogeom::algo::MeshSolidOptions::default();
    Some(time(|| {
        let mut model = Model::new();
        let built = ogeom::algo::solid_from_mesh(&mut model, &mesh, &options, T).unwrap();
        std::hint::black_box(&built);
    }))
}

/// The view both hidden-line benchmarks draw from: down a body diagonal.
fn diagonal_view() -> ogeom::hlr::View {
    ogeom::hlr::View::looking(Vector::new(-1.0, -1.0, -1.0), Vector::Z, T).unwrap()
}

/// Hidden lines on an imported part with exact silhouettes and visibility.
fn hlr_exact() -> Option<Stats> {
    let (document, solid) = corpus_part(PART)?;
    let view = diagonal_view();
    Some(time(|| {
        let drawing = ogeom::hlr::exact::project_exact(
            document.model(),
            &solid,
            &view,
            Deflection::default(),
            T,
        )
        .unwrap();
        std::hint::black_box(&drawing);
    }))
}

/// Hidden lines on the same part against its tessellation.
fn hlr_mesh() -> Option<Stats> {
    let (document, solid) = corpus_part(PART)?;
    let view = diagonal_view();
    Some(time(|| {
        let drawing =
            ogeom::hlr::project(document.model(), &solid, &view, Deflection::default(), T).unwrap();
        std::hint::black_box(&drawing);
    }))
}

/// The mid-size corpus part exact hidden lines are timed on: about 1600
/// drawn curves over a few hundred faces. The largest part takes tens of
/// seconds at one thread, too long to sample here.
const MID_PART: &str = "nist_ctc_04_asme1_rd.stp";

/// Exact hidden lines on a part with many faces and curves, where the
/// faces asked about each point dominate.
fn hlr_exact_mid() -> Option<Stats> {
    let (document, solid) = corpus_part(MID_PART)?;
    let view = diagonal_view();
    Some(time(|| {
        let drawing = ogeom::hlr::exact::project_exact(
            document.model(),
            &solid,
            &view,
            Deflection::default(),
            T,
        )
        .unwrap();
        std::hint::black_box(&drawing);
    }))
}

/// A B-spline face thickened both ways: the saddle z = (x^2 - y^2) / 20
/// over [-5, 5]^2, fitted as a cubic and bounded by its border iso-curves.
fn thick_spline() -> Option<Stats> {
    let mut model = Model::new();
    let n = 21;
    let rows: Vec<Vec<Point>> = (0..n)
        .map(|j| {
            let y = -5.0 + 10.0 * f64::from(j) / f64::from(n - 1);
            (0..n)
                .map(|i| {
                    let x = -5.0 + 10.0 * f64::from(i) / f64::from(n - 1);
                    Point::new(x, y, (x * x - y * y) / 20.0)
                })
                .collect()
        })
        .collect();
    let surface = ogeom::geom::fit::fit_surface_grid(&rows, 3, 1e-6, T)
        .ok()?
        .curve;
    let ((u0, u1), (v0, v1)) = surface.domain();
    let v: Vec<Shape> = [(u0, v0), (u1, v0), (u1, v1), (u0, v1)]
        .iter()
        .map(|(u, w)| {
            let at = surface.point_at(*u, *w, T).unwrap();
            ogeom::algo::make_vertex(&mut model, at).shape
        })
        .collect();
    let mut iso = |curve: ogeom::geom::BSplineCurve, from: &Shape, to: &Shape| {
        let range = curve.domain();
        ogeom::algo::make_edge_between(&mut model, curve.into(), range, from, to, T)
            .unwrap()
            .shape
    };
    let south = iso(surface.iso_v_curve(v0, T).unwrap(), &v[0], &v[1]);
    let east = iso(surface.iso_u_curve(u1, T).unwrap(), &v[1], &v[2]);
    let north = iso(surface.iso_v_curve(v1, T).unwrap(), &v[3], &v[2]);
    let west = iso(surface.iso_u_curve(u0, T).unwrap(), &v[0], &v[3]);
    let wires = [vec![south, east, north.reversed(), west.reversed()]];
    let face = ogeom::algo::make_face_with_pcurves(&mut model, surface.into(), &wires, T)
        .unwrap()
        .shape;
    Some(sample(
        || model.clone(),
        |mut model| {
            let solid = ogeom::offset::make_thick_sheet(&mut model, &face, 0.8, true, T).unwrap();
            std::hint::black_box(&solid.shape);
        },
    ))
}

/// Surface feet on a finely knotted patch: a wavy cubic of 100 by 100
/// knot spans over 50 by 50 mm, thirty-two targets above, below and
/// beside it projected one call at a time.
fn project_spline() -> Option<Stats> {
    use ogeom::geom::{BSplineSurface, SurfaceGeometry, project_on_surface};
    use ogeom::math::{ControlGrid, KnotVector};
    let n = 103;
    let points: Vec<Point> = (0..n)
        .flat_map(|i| {
            (0..n).map(move |j| {
                let x = 50.0 * f64::from(i) / f64::from(n - 1);
                let y = 50.0 * f64::from(j) / f64::from(n - 1);
                Point::new(x, y, 0.8 * (0.7 * x).sin() * (0.45 * y).cos())
            })
        })
        .collect();
    let knots = KnotVector::clamped_uniform(3, 103).ok()?;
    let grid = ControlGrid::new(points, 103, 103).ok()?;
    let surface: SurfaceGeometry = BSplineSurface::new(knots.clone(), knots, &grid, T)
        .ok()?
        .into();
    let targets: Vec<Point> = (0..32)
        .map(|k| {
            let t = f64::from(k);
            Point::new(
                (t * 7.3).rem_euclid(60.0) - 5.0,
                (t * 11.9).rem_euclid(60.0) - 5.0,
                (t * 0.37).sin() * 6.0,
            )
        })
        .collect();
    Some(time(|| {
        for &target in &targets {
            let found = project_on_surface(&surface, target, 24, T).unwrap();
            std::hint::black_box(found);
        }
    }))
}

/// A drum of radius 4 and height 6 with its top rim rolled to radius 1:
/// two planes, a cylinder and a torus band.
fn rimmed_drum(model: &mut Model) -> Option<Shape> {
    let drum = ogeom::algo::make_cylinder(model, Frame::WORLD, 4.0, 6.0, T)
        .ok()?
        .shape;
    let rim = explore_unique(model, &drum, ShapeType::Edge)
        .ok()?
        .into_iter()
        .find(|edge| {
            ogeom::algo::shape_bounds(model, edge, T)
                .ok()
                .and_then(|b| b.low())
                .is_some_and(|low| low.z > 6.0 - 1e-6)
        })?;
    Some(
        ogeom::fillet::fillet_edge(model, &drum, &rim, 1.0, T)
            .ok()?
            .shape,
    )
}

/// The rim-rolled drum grown by a half along its normals.
fn offset_part() -> Option<Stats> {
    let mut model = Model::new();
    let part = rimmed_drum(&mut model)?;
    Some(sample(
        || model.clone(),
        |mut model| {
            let grown = ogeom::offset::offset_shape(&mut model, &part, 0.5, T).unwrap();
            std::hint::black_box(&grown.shape);
        },
    ))
}

/// The rim-rolled drum hollowed to walls of a half, open at its bottom.
fn shell_part() -> Option<Stats> {
    let mut model = Model::new();
    let part = rimmed_drum(&mut model)?;
    let bottom = explore_unique(&model, &part, ShapeType::Face)
        .ok()?
        .into_iter()
        .find(|face| {
            ogeom::algo::shape_bounds(&model, face, T)
                .ok()
                .and_then(|b| b.high())
                .is_some_and(|high| high.z < 1e-6)
        })?;
    let opening = [bottom];
    Some(sample(
        || model.clone(),
        |mut model| {
            let hollow =
                ogeom::offset::make_thick_solid(&mut model, &part, &opening, 0.5, T).unwrap();
            std::hint::black_box(&hollow.shape);
        },
    ))
}

fn circle_wire(model: &mut Model, frame: Frame, radius: f64) -> Shape {
    let circle = Circle::new(frame, radius, T).unwrap();
    let edge = ogeom::algo::make_edge(
        model,
        CircleCurve::new(circle).into(),
        (0.0, std::f64::consts::TAU),
        T,
    )
    .unwrap()
    .shape;
    ogeom::algo::make_wire(model, &[edge], T).unwrap().shape
}

/// A pipe through two circular sections, radius 2 then 3, down a line.
fn pipe_circles() -> Option<Stats> {
    let mut model = Model::new();
    let along = Point::new(0.0, 20.0, 0.0);
    let small = circle_wire(
        &mut model,
        Frame::new(Point::ORIGIN, Direction::Y, Direction::X, T).unwrap(),
        2.0,
    );
    let large = circle_wire(
        &mut model,
        Frame::new(along, Direction::Y, Direction::X, T).unwrap(),
        3.0,
    );
    let spine = ogeom::algo::make_polygon(&mut model, &[Point::ORIGIN, along], false, T)
        .unwrap()
        .shape;
    let sections = [small, large];
    Some(sample(
        || model.clone(),
        |mut model| {
            let pipe =
                ogeom::offset::make_pipe_sections(&mut model, &sections, &spine, false, 1e-3, T)
                    .unwrap();
            std::hint::black_box(&pipe.shape);
        },
    ))
}

/// A square swept down a line, turned by a guide beside it: the section's
/// axis points where the guide crosses each station's plane.
fn pipe_guided() -> Option<Stats> {
    let mut model = Model::new();
    let square =
        [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)].map(|(x, z)| Point::new(x, 0.0, z));
    let wire = ogeom::algo::make_polygon(&mut model, &square, true, T)
        .unwrap()
        .shape;
    let plane = Plane::new(Frame::new(Point::ORIGIN, -Direction::Y, Direction::X, T).unwrap());
    let profile = ogeom::algo::make_face(&mut model, PlaneSurface::new(plane).into(), &[wire], T)
        .unwrap()
        .shape;
    let spine = ogeom::algo::make_polygon(
        &mut model,
        &[Point::ORIGIN, Point::new(0.0, 20.0, 0.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    let guide = ogeom::algo::make_polygon(
        &mut model,
        &[Point::new(5.0, 0.0, 0.0), Point::new(0.0, 20.0, 5.0)],
        false,
        T,
    )
    .unwrap()
    .shape;
    Some(sample(
        || model.clone(),
        |mut model| {
            let pipe = ogeom::offset::make_pipe_shell_with(
                &mut model,
                &profile,
                &spine,
                &ogeom::offset::PipeLaw::Auxiliary { guide: &guide },
                ogeom::offset::PipeCorners::Mitre,
                1e-3,
                T,
            )
            .unwrap();
            std::hint::black_box(&pipe.shape);
        },
    ))
}

/// Import: the smallest NIST corpus part, read and healed.
fn import_ftc11() -> Option<Stats> {
    let text = corpus(PART)?;
    Some(time(|| {
        let import = ogeom::io::read_step(&text, T).unwrap();
        std::hint::black_box(import.solids.len());
    }))
}

/// Import: the largest corpus part, so the reader is measured at a size the
/// smallest file cannot show.
fn import_ctc02() -> Option<Stats> {
    let text = corpus(LARGE_PART)?;
    Some(time(|| {
        let import = ogeom::io::read_step(&text, T).unwrap();
        std::hint::black_box(import.solids.len());
    }))
}

/// STEP export of the largest corpus part.
fn step_write() -> Option<Stats> {
    let import = ogeom::io::read_step(&corpus(LARGE_PART)?, T).ok()?;
    Some(time(|| {
        let text = ogeom::io::write_step(&import.document, T).unwrap();
        std::hint::black_box(text.len());
    }))
}

/// IGES import of the largest corpus part, written as IGES by the kernel first: the
/// corpus has no IGES file of its own.
fn iges_read() -> Option<Stats> {
    let (document, _) = corpus_part(LARGE_PART)?;
    let text = ogeom::io::write_iges(&document, T).ok()?;
    Some(time(|| {
        let import = ogeom::io::read_iges(&text, T).unwrap();
        std::hint::black_box(import.solids.len());
    }))
}

/// A closed torus mesh of a million triangles, 1000 by 500 quads, each
/// triangle with its own three vertices as an STL file holds them.
fn torus_soup() -> Triangulation {
    let (around, across) = (1000_u32, 500_u32);
    let at = |i: u32, j: u32| {
        let u = std::f64::consts::TAU * f64::from(i % around) / f64::from(around);
        let v = std::f64::consts::TAU * f64::from(j % across) / f64::from(across);
        let r = 40.0 + 10.0 * v.cos();
        Point::new(r * u.cos(), r * u.sin(), 10.0 * v.sin())
    };
    let mut mesh = Triangulation::new();
    for i in 0..around {
        for j in 0..across {
            let quad = [at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)];
            for corners in [[0, 1, 2], [0, 2, 3]] {
                #[allow(clippy::cast_possible_truncation)]
                let first = mesh.positions.len() as u32;
                mesh.positions.extend(corners.map(|c| quad[c]));
                mesh.triangles.push([first, first + 1, first + 2]);
            }
        }
    }
    mesh
}

/// Welding a million-triangle soup into a closed mesh.
fn weld_1m() -> Option<Stats> {
    let soup = torus_soup();
    Some(time(|| {
        let welded = soup.welded(T);
        std::hint::black_box(welded.positions.len());
    }))
}

/// Reading a million-triangle binary STL: parse and weld.
fn stl_read_1m() -> Option<Stats> {
    let bytes = ogeom::io::write(&torus_soup(), ogeom::io::Encoding::Binary).ok()?;
    Some(time(|| {
        let mesh = ogeom::io::read(&bytes, T).unwrap();
        std::hint::black_box(mesh.triangles.len());
    }))
}

/// The value recorded for `key` in a baseline: a `"key": number` line.
fn recorded(baseline: &str, key: &str) -> Option<f64> {
    baseline.lines().find_map(|line| {
        let line = line.trim().trim_end_matches(',');
        let (name, value) = line.split_once(':')?;
        (name.trim().trim_matches('"') == key).then(|| value.trim().parse::<f64>().ok())?
    })
}

fn main() {
    let mut threads = None;
    let mut filter = String::new();
    let mut check = None;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--threads" => {
                threads = Some(
                    args.next()
                        .and_then(|n| n.parse::<usize>().ok())
                        .filter(|&n| n > 0)
                        .expect("--threads takes a positive count"),
                );
            }
            "--filter" => filter = args.next().expect("--filter takes a name fragment"),
            "--check" => check = Some(args.next().expect("--check takes a baseline file")),
            other => panic!("unknown argument {other}; see the crate documentation"),
        }
    }
    let baseline = check
        .as_ref()
        .map(|path| std::fs::read_to_string(path).expect("baseline file"));
    let baseline_threads = baseline.as_deref().and_then(|b| recorded(b, "threads"));
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    if let Some(count) = threads.or(baseline_threads.map(|t| t as usize)) {
        ogeom::core::parallel::set_threads(count);
    }
    let running = ogeom::core::parallel::threads();

    // The spin before and after, the faster kept: a burst of load during
    // one of them does not skew every ratio.
    let before = calibration();
    let mut results = Vec::new();
    for (name, bench) in benchmarks() {
        if name.contains(filter.as_str())
            && let Some(stats) = bench()
        {
            results.push((name, stats));
        }
    }
    let spin = before.min.min(calibration().min);

    eprintln!("{running} threads; calibration spin {:.2} ms", spin * 1e3);
    if baseline.is_some() {
        #[allow(clippy::cast_precision_loss)]
        if baseline_threads.is_some_and(|t| t != running as f64) {
            eprintln!("the baseline was recorded at another thread count: ratios do not compare");
        }
        println!("name                min(ms)  med(ms)  mad(ms)   n   ratio  baseline   drift");
    } else {
        println!("{{");
        let comma = if results.is_empty() { "" } else { "," };
        println!("  \"threads\": {running}{comma}");
        eprintln!("name                min(ms)  med(ms)  mad(ms)   n   ratio");
    }
    for (i, (name, s)) in results.iter().enumerate() {
        let ratio = s.min / spin;
        let row = format!(
            "{name:<18}{:>9.2}{:>9.2}{:>9.2}{:>5}{ratio:>8.3}",
            s.min * 1e3,
            s.median * 1e3,
            s.mad * 1e3,
            s.samples
        );
        match &baseline {
            Some(baseline) => match recorded(baseline, name) {
                Some(base) => {
                    println!("{row}{base:>10.3}  {:>+6.1}%", (ratio / base - 1.0) * 100.0)
                }
                None => println!("{row}       new"),
            },
            None => {
                let comma = if i + 1 == results.len() { "" } else { "," };
                println!("  \"{name}\": {ratio:.4}{comma}");
                eprintln!("{row}");
            }
        }
    }
    if baseline.is_none() {
        println!("}}");
    }
}
