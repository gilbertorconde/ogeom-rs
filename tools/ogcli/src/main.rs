//! `ogcli` — build and inspect shapes from the shell.
//!
//! Deliberately thin: it exists so the kernel is usable and testable without a
//! host application, and so a wrong result is easy to reproduce from a command
//! line rather than only from a test.

use std::process::ExitCode;

use og::{
    algo::make_box,
    core::Tolerances,
    math::Frame,
    topo::{Model, ShapeType, explore_unique},
};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match args.first().map(String::as_str) {
        Some("box") => match run_box(&args[1..]) {
            Ok(()) => ExitCode::SUCCESS,
            Err(why) => {
                eprintln!("ogcli: {why}");
                ExitCode::FAILURE
            }
        },
        Some("version") | None => {
            println!("openGeometry {}", og::VERSION);
            ExitCode::SUCCESS
        }
        Some(other) => {
            eprintln!("ogcli: unknown command '{other}'");
            eprintln!("usage: ogcli [version | box <dx> <dy> <dz>]");
            ExitCode::FAILURE
        }
    }
}

/// Build a box and report what it is made of.
fn run_box(args: &[String]) -> Result<(), String> {
    let dims: Vec<f64> = args
        .iter()
        .map(|a| {
            a.parse::<f64>()
                .map_err(|_| format!("'{a}' is not a number"))
        })
        .collect::<Result<_, _>>()?;
    let [dx, dy, dz] = dims.as_slice() else {
        return Err("usage: ogcli box <dx> <dy> <dz>".into());
    };

    let tol = Tolerances::millimetres();
    let mut model = Model::new();
    let built =
        make_box(&mut model, Frame::WORLD, (*dx, *dy, *dz), tol).map_err(|e| e.to_string())?;

    println!("box {dx} x {dy} x {dz}");
    for kind in [
        ShapeType::Shell,
        ShapeType::Face,
        ShapeType::Wire,
        ShapeType::Edge,
        ShapeType::Vertex,
    ] {
        let found = explore_unique(&model, &built.shape, kind).map_err(|e| e.to_string())?;
        println!(
            "  {:>8}: {}",
            format!("{kind:?}").to_lowercase(),
            found.len()
        );
    }
    Ok(())
}
