//! Whether a read of a shape made from the reads kept of other shapes
//! holding the same faces ([`ogeom_topo::Model::kept_reads`]) is checked
//! against the same read made afresh: always in a debug build, and in any
//! build inside `derived_reads_checked`, with the `test-hooks` feature.

use std::cell::Cell;

thread_local! {
    /// The faces whose reads were taken from kept reads, counted while
    /// `derived_reads_checked` runs.
    static CHECKED: Cell<Option<usize>> = const { Cell::new(None) };
}

/// Whether a read taken in part from kept reads on this thread is to be
/// compared with the same read made afresh, and a difference panicked on.
#[doc(hidden)]
#[must_use]
pub fn checks_derived_reads() -> bool {
    cfg!(debug_assertions) || CHECKED.get().is_some()
}

/// Count `faces` whose read was taken from a kept read.
#[doc(hidden)]
pub fn note_derived(faces: usize) {
    if let Some(n) = CHECKED.get() {
        CHECKED.set(Some(n + faces));
    }
}

/// Run `f` with every read taken in part from kept reads on this thread
/// compared with the same read made afresh, panicking where they differ;
/// and how many faces' reads were taken so.
#[cfg(feature = "test-hooks")]
pub fn derived_reads_checked<R>(f: impl FnOnce() -> R) -> (R, usize) {
    struct Restore(Option<usize>);
    impl Drop for Restore {
        fn drop(&mut self) {
            CHECKED.set(self.0);
        }
    }
    let restore = Restore(CHECKED.replace(Some(0)));
    let out = f();
    let faces = CHECKED.get().unwrap_or(0);
    drop(restore);
    (out, faces)
}
