//! Checks that the kernel meshes on `wasm32-unknown-unknown`, a target whose
//! standard library has no clock and cannot start a thread.
//!
//! The tests live in `tests/` and run on both sides: natively with
//! `cargo test -p ogeom-wasm`, and in Node with
//! `cargo test -p ogeom-wasm --target wasm32-unknown-unknown` under
//! `wasm-bindgen-test-runner`. Both assert the same triangle counts.
//!
//! `tests/threads.rs` lends the kernel a pool of web workers. It needs a
//! shared-memory build, a nightly toolchain with `rust-src`, and a headless
//! Chrome with `chromedriver` on the `PATH` (Node has no web workers), with
//! `wasm-bindgen-test-runner` at the version `Cargo.lock` holds:
//!
//! ```text
//! RUSTFLAGS='-C target-feature=+atomics,+bulk-memory
//!   -C link-arg=--shared-memory -C link-arg=--import-memory
//!   -C link-arg=--max-memory=1073741824
//!   -C link-arg=--export=__wasm_init_tls -C link-arg=--export=__tls_size
//!   -C link-arg=--export=__tls_align -C link-arg=--export=__tls_base' \
//! CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_RUNNER=wasm-bindgen-test-runner \
//! cargo +nightly test -p ogeom-wasm --test threads \
//!   --target wasm32-unknown-unknown -Z build-std=panic_abort,std
//! ```
