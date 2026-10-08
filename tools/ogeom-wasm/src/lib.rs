//! Checks that the kernel meshes on `wasm32-unknown-unknown`, a target whose
//! standard library has no clock and cannot start a thread.
//!
//! The tests live in `tests/` and run on both sides: natively with
//! `cargo test -p ogeom-wasm`, and in Node with
//! `cargo test -p ogeom-wasm --target wasm32-unknown-unknown` under
//! `wasm-bindgen-test-runner`. Both assert the same triangle counts.
