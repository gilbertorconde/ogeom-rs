//! Valid files, mangled: bytes changed, cut out, repeated, the file cut
//! short. Every reader answers each one, with a result or an error, and
//! never panics, overflows its stack or allocates past what it holds.
//!
//! A seeded stand-in for fuzzing: the same mangles every run, so a failure
//! reproduces from its seed.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::Frame;
use ogeom::topo::Model;

const T: Tolerances = Tolerances::millimetres();
const ROUNDS: u64 = 200;

/// A small deterministic generator (xorshift), so the mangles are the same
/// on every machine.
struct Seeded(u64);

impl Seeded {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    fn below(&mut self, n: usize) -> usize {
        usize::try_from(self.next() % (n.max(1) as u64)).unwrap()
    }
}

/// One mangle of `bytes`, a few changes deep. Mostly a digit changed to
/// another (the syntax survives, so the change reaches the builders: an
/// id now names another entity, a count or a degree grows), sometimes a
/// byte changed to anything, a short stretch cut out or repeated, or the
/// tail cut off.
fn mangled(bytes: &[u8], seed: u64) -> Vec<u8> {
    let mut rng = Seeded(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1);
    let mut out = bytes.to_vec();
    for _ in 0..=rng.below(3) {
        if out.is_empty() {
            break;
        }
        let at = rng.below(out.len());
        match rng.below(20) {
            0..=11 => {
                let digits: Vec<usize> = (0..out.len())
                    .filter(|&i| out[i].is_ascii_digit())
                    .collect();
                if !digits.is_empty() {
                    let i = digits[rng.below(digits.len())];
                    out[i] = b"0123456789"[rng.below(10)];
                }
            }
            12..=14 => {
                const SHARP: &[u8] = b"0123456789#(),;=-.E'$*";
                out[at] = SHARP[rng.below(SHARP.len())];
            }
            15 => out[at] = u8::try_from(rng.below(256)).unwrap(),
            16 | 17 => {
                let end = (at + 1 + rng.below(8)).min(out.len());
                out.drain(at..end);
            }
            18 => {
                let end = (at + 1 + rng.below(8)).min(out.len());
                let copy = out[at..end].to_vec();
                out.splice(at..at, copy);
            }
            _ => out.truncate(at),
        }
    }
    out
}

fn drilled_block() -> (Model, ogeom::topo::Shape) {
    let mut model = Model::new();
    let block = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 8.0, 6.0), T)
        .unwrap()
        .shape;
    let bore = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 2.0, 6.0, T)
        .unwrap()
        .shape;
    let part = ogeom::boolean::cut(&mut model, &block, &bore, T)
        .unwrap()
        .shape;
    (model, part)
}

fn document_of_a_drilled_block() -> ogeom::doc::Document {
    let (model, part) = drilled_block();
    let mut document = ogeom::doc::Document::over(model);
    document.add_part("part", part);
    document
}

#[test]
fn mangled_step_files_are_answered() {
    let text = ogeom::io::write_step(&document_of_a_drilled_block(), T).unwrap();
    for seed in 0..ROUNDS {
        let bytes = mangled(text.as_bytes(), seed);
        let _ = ogeom::io::read_step(&String::from_utf8_lossy(&bytes), T);
    }
}

#[test]
fn mangled_iges_files_are_answered() {
    let text = ogeom::io::write_iges(&document_of_a_drilled_block(), T).unwrap();
    for seed in 0..ROUNDS {
        let bytes = mangled(text.as_bytes(), seed);
        let _ = ogeom::io::read_iges(&String::from_utf8_lossy(&bytes), T);
    }
}

#[test]
fn mangled_mesh_files_are_answered() {
    let (model, part) = drilled_block();
    let mesh =
        ogeom::mesh::triangulate(&model, &part, ogeom::mesh::Deflection::default(), T).unwrap();
    let binary = ogeom::io::write(&mesh, ogeom::io::Encoding::Binary).unwrap();
    let ascii = ogeom::io::write(&mesh, ogeom::io::Encoding::Ascii).unwrap();
    let package = ogeom::io::write_3mf(&[ogeom::io::threemf::Object {
        mesh: &mesh,
        name: Some("part".into()),
    }]);
    for seed in 0..ROUNDS {
        let _ = ogeom::io::read(&mangled(&binary, seed), T);
        let _ = ogeom::io::read(&mangled(&ascii, seed), T);
        let _ = ogeom::io::read_3mf(&mangled(&package, seed), T);
    }
}
