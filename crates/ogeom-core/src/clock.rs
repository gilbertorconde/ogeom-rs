//! A monotonic clock for the kernel's diagnostic timings, on every target.
//!
//! [`Instant`] is `std::time::Instant` wherever the standard library has a
//! clock: native targets and WASI. On `wasm32-unknown-unknown` the standard
//! library has none and `std::time::Instant::now` panics there, so this one
//! reads as a fixed point and every duration is zero. With the `web` feature
//! it reads the browser's `performance.now()` instead.
//!
//! Timings are diagnostics. No answer of any operation depends on them, so a
//! clock that reads zero changes what a report says about time and nothing
//! else. `std::time::Instant::now` is refused by clippy
//! (`disallowed-methods` in `clippy.toml`) in the kernel's crates.

use std::time::Duration;

#[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
type Inner = std::time::Instant;

#[cfg(all(target_arch = "wasm32", target_os = "unknown", feature = "web"))]
type Inner = web_time::Instant;

/// No clock: every instant is the same instant.
#[cfg(all(target_arch = "wasm32", target_os = "unknown", not(feature = "web")))]
#[derive(Clone, Copy, Debug)]
struct Inner;

#[cfg(all(target_arch = "wasm32", target_os = "unknown", not(feature = "web")))]
impl Inner {
    fn now() -> Self {
        Self
    }

    fn elapsed(self) -> Duration {
        Duration::ZERO
    }
}

/// A point in time, measured by the platform's monotonic clock where it has
/// one.
#[derive(Clone, Copy, Debug)]
pub struct Instant(Inner);

impl Instant {
    /// Now, or a fixed point where the platform has no clock.
    #[must_use]
    #[inline]
    #[allow(
        clippy::disallowed_methods,
        reason = "the one place the platform clock is read"
    )]
    pub fn now() -> Self {
        Self(Inner::now())
    }

    /// Time since `self`; zero where the platform has no clock.
    #[must_use]
    #[inline]
    pub fn elapsed(&self) -> Duration {
        self.0.elapsed()
    }
}
