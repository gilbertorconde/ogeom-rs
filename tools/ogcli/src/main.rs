//! `ogcli` — build, inspect, convert and tessellate shapes from the shell.
//!
//! Deliberately thin: it exists so the kernel is usable and testable without a
//! host application, and so wrong results are easy to reproduce from a corpus.

fn main() {
    println!("openGeometry {}", og::VERSION);
}
