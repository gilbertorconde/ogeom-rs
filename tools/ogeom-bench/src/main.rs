//! Benchmarks over the kernel's hot paths — the sequencing instrument the
//! scope promised.
//!
//! The harness is deliberately its own: medians over repeated runs, wall
//! clock, no dependency. Absolute times move with the machine, so every run
//! also times a fixed arithmetic spin — the *calibration* — and the
//! comparison mode reports each benchmark as a multiple of it. Ratios
//! travel between machines; milliseconds do not.
//!
//! `ogeom-bench` prints this machine's numbers. `ogeom-bench --check
//! <baseline.json>` compares calibrated ratios against a recorded baseline
//! and reports the drift, informationally: performance is watched here, not
//! gated, because a loaded CI box would turn a real gate into a coin flip.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "a reporting tool")]

use std::time::Instant;

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// Median wall time of `runs` executions, in seconds.
fn median(runs: usize, mut f: impl FnMut()) -> f64 {
    let mut times: Vec<f64> = (0..runs)
        .map(|_| {
            let start = Instant::now();
            f();
            start.elapsed().as_secs_f64()
        })
        .collect();
    times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    times[times.len() / 2]
}

/// A fixed arithmetic spin whose time stands for this machine's speed.
fn calibration() -> f64 {
    median(5, || {
        let mut acc = 0.0f64;
        for i in 0..4_000_000u64 {
            #[allow(clippy::cast_precision_loss)]
            let x = (i as f64).mul_add(1.000_000_1, acc).sin();
            acc = x * 0.999;
        }
        std::hint::black_box(acc);
    })
}

fn benchmarks() -> Vec<(&'static str, f64)> {
    let mut out = Vec::new();

    // Construction: a thousand boxes into one model.
    out.push((
        "construct_boxes",
        median(5, || {
            let mut model = Model::new();
            for i in 0..1000 {
                #[allow(clippy::cast_precision_loss)]
                let dx = i as f64;
                let frame =
                    Frame::new(Point::new(dx, 0.0, 0.0), Direction::Z, Direction::X, T).unwrap();
                ogeom::algo::make_box(&mut model, frame, (1.0, 1.0, 1.0), T).unwrap();
            }
            std::hint::black_box(&model);
        }),
    ));

    // Traversal: exploring a box's sub-shapes, ten thousand times.
    {
        let mut model = Model::new();
        let solid = ogeom::algo::make_box(&mut model, Frame::WORLD, (2.0, 3.0, 4.0), T)
            .unwrap()
            .shape;
        out.push((
            "traverse_box",
            median(5, || {
                for _ in 0..10_000 {
                    let faces = explore_unique(&model, &solid, ShapeType::Face).unwrap();
                    std::hint::black_box(faces.len());
                }
            }),
        ));
    }

    // Tessellation: a torus at the default deflection.
    out.push((
        "tessellate_torus",
        median(5, || {
            let mut model = Model::new();
            let solid = ogeom::algo::make_torus(&mut model, Frame::WORLD, 20.0, 5.0, T)
                .unwrap()
                .shape;
            let done =
                ogeom::mesh::tessellate(&mut model, &solid, Deflection::default(), T).unwrap();
            std::hint::black_box(done.triangles);
        }),
    ));

    // The boolean: the drilled box, the pipeline's pinned proof.
    out.push((
        "boolean_drill",
        median(5, || {
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
        }),
    ));

    // The boolean with many faces in play: a block drilled four times, then
    // cut by a torus. The classifier is asked once per face piece and asks
    // the *other* solid's every face, so this is where that product shows.
    out.push((
        "boolean_many_faces",
        median(5, || {
            let mut model = Model::new();
            let mut block = ogeom::algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T)
                .unwrap()
                .shape;
            for (x, y) in [(5.0, 5.0), (15.0, 5.0), (5.0, 15.0), (15.0, 15.0)] {
                let frame =
                    Frame::new(Point::new(x, y, -1.0), Direction::Z, Direction::X, T).unwrap();
                let drill = ogeom::algo::make_cylinder(&mut model, frame, 2.0, 12.0, T)
                    .unwrap()
                    .shape;
                block = ogeom::boolean::cut(&mut model, &block, &drill, T)
                    .unwrap()
                    .shape;
            }
            std::hint::black_box(&block);
        }),
    ));

    // The boolean with no closed form: crossed cylinders, whose sections only
    // the marcher can trace. `boolean_drill` is analytic end to end and never
    // reaches that machinery.
    out.push((
        "boolean_marched",
        median(5, || {
            let mut model = Model::new();
            let upright = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 1.0, 4.0, T)
                .unwrap()
                .shape;
            let frame =
                Frame::new(Point::new(-2.0, 0.0, 2.0), Direction::X, Direction::Y, T).unwrap();
            let across = ogeom::algo::make_cylinder(&mut model, frame, 0.6, 4.0, T)
                .unwrap()
                .shape;
            let both = ogeom::boolean::fuse(&mut model, &upright, &across, T).unwrap();
            std::hint::black_box(&both.shape);
        }),
    ));

    // A constant-radius fillet: the blend machinery, which no benchmark reached.
    out.push((
        "fillet_block",
        median(5, || {
            let mut model = Model::new();
            let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
                .unwrap()
                .shape;
            let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
            let rolled =
                ogeom::fillet::fillet_edges(&mut model, &block, &edges[..1], 1.0, T).unwrap();
            std::hint::black_box(&rolled.shape);
        }),
    ));

    // Whole-shape tessellation of a real imported part: the sequential face
    // loop that mass properties falls into whenever a face is not analytic,
    // which is guaranteed on an imported b-rep. A primitive will not do —
    // a torus answers from its closed form and never meshes at all.
    {
        let path = format!(
            "{}/../../tests/corpus/nist_ftc_11_asme1_rb.stp",
            env!("CARGO_MANIFEST_DIR")
        );
        if let Ok(text) = std::fs::read_to_string(&path)
            && let Ok(import) = ogeom::io::read_step(&text, T)
            && let Some(solid) = import.solids.first().cloned()
        {
            let model = import.document;
            out.push((
                "tessellate_part",
                median(5, || {
                    let mesh = ogeom::mesh::triangulate(
                        model.model(),
                        &solid,
                        Deflection::with_chord(1e-2).unwrap(),
                        T,
                    )
                    .unwrap();
                    std::hint::black_box(mesh.triangles.len());
                }),
            ));
        }
    }

    // Import: the smallest corpus part, read and healed.
    {
        let path = format!(
            "{}/../../tests/corpus/nist_ftc_11_asme1_rb.stp",
            env!("CARGO_MANIFEST_DIR")
        );
        if let Ok(text) = std::fs::read_to_string(&path) {
            out.push((
                "import_ftc11",
                median(3, || {
                    let import = ogeom::io::read_step(&text, T).unwrap();
                    std::hint::black_box(import.solids.len());
                }),
            ));
        }
    }

    // Import: the largest corpus part, so the reader is measured at a size
    // the smallest file cannot show.
    {
        let path = format!(
            "{}/../../tests/corpus/nist_ctc_02_asme1_rc.stp",
            env!("CARGO_MANIFEST_DIR")
        );
        if let Ok(text) = std::fs::read_to_string(&path) {
            out.push((
                "import_ctc02",
                median(3, || {
                    let import = ogeom::io::read_step(&text, T).unwrap();
                    std::hint::black_box(import.solids.len());
                }),
            ));
        }
    }

    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let spin = calibration();
    let results = benchmarks();

    if args.len() >= 3 && args[1] == "--check" {
        let baseline = std::fs::read_to_string(&args[2]).expect("baseline file");
        println!("name              now(ms)   ratio    baseline  drift");
        for (name, seconds) in &results {
            let ratio = seconds / spin;
            let recorded = baseline.lines().find_map(|line| {
                let line = line.trim().trim_end_matches(',');
                let (key, value) = line.split_once(':')?;
                (key.trim().trim_matches('"') == *name).then(|| value.trim().parse::<f64>().ok())?
            });
            match recorded {
                Some(base) => println!(
                    "{name:<18}{:>8.2}{ratio:>8.2}{base:>10.2}  {:>+6.1}%",
                    seconds * 1e3,
                    (ratio / base - 1.0) * 100.0
                ),
                None => println!("{name:<18}{:>8.2}{ratio:>8.2}       new", seconds * 1e3),
            }
        }
        return;
    }

    // Plain run: print, and emit the baseline JSON to stdout on request.
    eprintln!("calibration spin: {:.2} ms", spin * 1e3);
    println!("{{");
    for (i, (name, seconds)) in results.iter().enumerate() {
        let comma = if i + 1 == results.len() { "" } else { "," };
        println!("  \"{name}\": {:.4}{comma}", seconds / spin);
        eprintln!(
            "{name:<18}{:>8.2} ms  (x{:.2} spin)",
            seconds * 1e3,
            seconds / spin
        );
    }
    println!("}}");
}
