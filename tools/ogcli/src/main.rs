//! `ogcli` — build, measure and export shapes from the shell.
//!
//! Deliberately thin: it exists so the kernel is usable and testable without a
//! host application, and so a wrong result is easy to reproduce from a command
//! line rather than only from a test.
//!
//! # It reports, and it also complains
//!
//! Every command that produces a shape runs the validity check over it and says
//! what it found. A tool that only prints the answer teaches you to trust the
//! answer; one that prints the answer *and* what is wrong with the shape it
//! came from teaches you when not to.

use std::process::ExitCode;

use og::{
    algo::{
        Severity, check, linear_properties, make_box, make_cone, make_cylinder, make_sphere,
        make_torus, make_wedge, surface_properties, volume_properties,
    },
    core::{OgResult, Tolerances, og_err},
    io::{Encoding, native, write as write_stl},
    math::Frame,
    mesh::{Deflection, triangulate},
    topo::{Model, Shape, ShapeType, explore_unique},
};

const TOL: Tolerances = Tolerances::millimetres();

const USAGE: &str = "\
usage: ogcli <command> [args]

  version
  box       <dx> <dy> <dz>
  cylinder  <radius> <height>
  sphere    <radius>
  cone      <base-radius> <top-radius> <height>
  torus     <major-radius> <minor-radius>
  wedge     <dx> <dy> <dz> <top-dx> <top-dy>

Every shape command accepts, after its dimensions:
  --deflection <chord>   how finely to tessellate (default 0.1)
  --stl <path>           write the tessellation as binary STL
  --ascii                with --stl, write the ASCII encoding instead
  --og <path>            write the whole shape as native .og text
  --no-mesh              with --og, leave the cached tessellation out
  --view <path>          render the tessellation to a PPM image
";

/// The primitives `build` knows how to make.
const SHAPES: [&str; 6] = ["box", "cylinder", "sphere", "cone", "torus", "wedge"];

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();

    let outcome = match args.first().map(String::as_str) {
        Some("version") | None => {
            println!("ogeom {}", og::VERSION);
            return ExitCode::SUCCESS;
        }
        Some("help" | "--help" | "-h") => {
            print!("{USAGE}");
            return ExitCode::SUCCESS;
        }
        Some(name) if SHAPES.contains(&name) => run(name, &args[1..]),
        Some(other) => {
            eprintln!("ogcli: unknown command '{other}'");
            eprint!("{USAGE}");
            return ExitCode::FAILURE;
        }
    };

    match outcome {
        Ok(()) => ExitCode::SUCCESS,
        Err(why) => {
            eprintln!("ogcli: {why}");
            ExitCode::FAILURE
        }
    }
}

/// The options that follow a shape's dimensions.
struct Options {
    deflection: Deflection,
    stl: Option<String>,
    encoding: Encoding,
    native: Option<String>,
    write_options: native::WriteOptions,
    view: Option<String>,
}

/// Build one primitive and report on it.
fn run(name: &str, args: &[String]) -> Result<(), String> {
    let (numbers, options) = split(args)?;
    let mut model = Model::new();
    let shape = build(&mut model, name, &numbers).map_err(|e| e.to_string())?;
    report(&model, &shape, name, &options).map_err(|e| e.to_string())
}

/// Build the named primitive from its dimensions.
fn build(model: &mut Model, name: &str, d: &[f64]) -> OgResult<Shape> {
    let want = |n: usize| -> OgResult<()> {
        if d.len() == n {
            Ok(())
        } else {
            Err(og_err!(
                Construction,
                "{name} takes {n} dimension(s), got {}",
                d.len()
            ))
        }
    };
    let built = match name {
        "box" => {
            want(3)?;
            make_box(model, Frame::WORLD, (d[0], d[1], d[2]), TOL)?
        }
        "cylinder" => {
            want(2)?;
            make_cylinder(model, Frame::WORLD, d[0], d[1], TOL)?
        }
        "sphere" => {
            want(1)?;
            make_sphere(model, Frame::WORLD, d[0], TOL)?
        }
        "cone" => {
            want(3)?;
            make_cone(model, Frame::WORLD, d[0], d[1], d[2], TOL)?
        }
        "torus" => {
            want(2)?;
            make_torus(model, Frame::WORLD, d[0], d[1], TOL)?
        }
        "wedge" => {
            want(5)?;
            make_wedge(model, Frame::WORLD, (d[0], d[1], d[2]), (d[3], d[4]), TOL)?
        }
        other => return Err(og_err!(Construction, "no primitive named '{other}'")),
    };
    Ok(built.shape)
}

/// Print the counts, the measurements and the diagnosis, then export if asked.
fn report(model: &Model, shape: &Shape, name: &str, options: &Options) -> OgResult<()> {
    println!("{name}:");
    for kind in [
        ShapeType::Solid,
        ShapeType::Shell,
        ShapeType::Face,
        ShapeType::Wire,
        ShapeType::Edge,
        ShapeType::Vertex,
    ] {
        let count = explore_unique(model, shape, kind)?.len();
        if count > 0 {
            println!("  {:>7}: {count}", format!("{kind:?}").to_lowercase());
        }
    }

    let volume = volume_properties(model, shape, options.deflection, TOL)?;
    let area = surface_properties(model, shape, options.deflection, TOL)?;
    let length = linear_properties(model, shape, options.deflection, TOL)?;
    println!("   volume: {:.9}", volume.mass);
    println!("     area: {:.9}", area.mass);
    println!("   length: {:.9}", length.mass);
    println!(
        "   centre: {:.9} {:.9} {:.9}",
        volume.centre.x, volume.centre.y, volume.centre.z
    );
    // The measurements come from the tessellation, so the deflection is part of
    // the answer rather than a setting. Printing it beside the numbers is what
    // lets a reader tell "0.6% under, as expected" from "wrong".
    println!("  at a chord deflection of {}", volume.deflection);

    let found = check(model, shape, TOL)?;
    if found.is_valid() {
        println!("  valid");
    } else {
        for problem in &found.problems {
            println!("  {problem}");
        }
        if found.worst() == Some(Severity::Broken) {
            println!("  -- the measurements above are not to be trusted");
        }
    }

    if let Some(path) = &options.view {
        // A picture catches what a number does not: a face wound inside out, a
        // hole where two surfaces failed to meet, a seam that did not weld.
        let mesh = triangulate(model, shape, options.deflection, TOL)?;
        let bounds = og::algo::shape_bounds(model, shape, TOL)?;
        let camera = ogview::Camera::framing(&bounds, og::math::Vector::new(1.0, -1.0, 0.7), TOL)?;
        let image = ogview::render(&mesh, &camera, &ogview::Style::default(), TOL)?;
        std::fs::write(path, image.to_ppm())
            .map_err(|e| og_err!(NotDone, "could not write {path}: {e}"))?;
        println!("  rendered {}x{} to {path}", image.width, image.height);
    }

    if let Some(path) = &options.native {
        // The whole shape, not its tessellation: topology, geometry,
        // tolerances and provenance, written so it reads back as the same
        // document. `diff` is the comparison tool.
        let text = native::write(model, std::slice::from_ref(shape), options.write_options)?;
        std::fs::write(path, &text).map_err(|e| og_err!(NotDone, "could not write {path}: {e}"))?;
        println!("  wrote {path} ({} lines)", text.lines().count());
    }

    if let Some(path) = &options.stl {
        let mesh = triangulate(model, shape, options.deflection, TOL)?;
        let bytes = write_stl(&mesh, options.encoding)?;
        std::fs::write(path, &bytes)
            .map_err(|e| og_err!(NotDone, "could not write {path}: {e}"))?;
        println!(
            "  wrote {} triangles to {path} ({} bytes)",
            mesh.triangle_count(),
            bytes.len()
        );
    }
    Ok(())
}

/// Split the dimensions from the options that follow them.
fn split(args: &[String]) -> Result<(Vec<f64>, Options), String> {
    let mut numbers = Vec::new();
    let mut options = Options {
        deflection: Deflection::default(),
        stl: None,
        encoding: Encoding::Binary,
        native: None,
        write_options: native::WriteOptions::default(),
        view: None,
    };

    let mut rest = args.iter();
    while let Some(arg) = rest.next() {
        match arg.as_str() {
            "--deflection" => {
                let value = rest
                    .next()
                    .ok_or("--deflection needs a number")?
                    .parse::<f64>()
                    .map_err(|_| "--deflection needs a number")?;
                options.deflection = Deflection {
                    chord: value,
                    ..Deflection::default()
                };
                options
                    .deflection
                    .validate()
                    .map_err(|e| format!("--deflection: {e}"))?;
            }
            "--stl" => options.stl = Some(rest.next().ok_or("--stl needs a path")?.clone()),
            "--ascii" => options.encoding = Encoding::Ascii,
            "--og" => options.native = Some(rest.next().ok_or("--og needs a path")?.clone()),
            "--no-mesh" => options.write_options.triangulations = false,
            "--view" => options.view = Some(rest.next().ok_or("--view needs a path")?.clone()),
            other if other.starts_with("--") => {
                return Err(format!("unknown option '{other}'"));
            }
            other => numbers.push(
                other
                    .parse::<f64>()
                    .map_err(|_| format!("'{other}' is not a number"))?,
            ),
        }
    }
    Ok((numbers, options))
}
