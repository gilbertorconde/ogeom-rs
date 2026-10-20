//! Counts of the work a boolean does, for the workspace's own tests: a
//! measure of what a boolean reaches that, unlike its time, does not move
//! with the machine it runs on. Present only with the `test-hooks` feature.

use std::cell::Cell;

/// What the booleans run on one thread did, summed.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Work {
    /// Faces gathered into the general fuse, both operands', each time one
    /// is run.
    pub faces_gathered: usize,
    /// Pairs of faces whose boxes meet, each asked for its sections.
    pub face_pairs: usize,
    /// Pieces the gathered faces were split into and classified.
    pub pieces: usize,
}

thread_local! {
    static COUNTED: Cell<Work> = const {
        Cell::new(Work { faces_gathered: 0, face_pairs: 0, pieces: 0 })
    };
}

/// Run `f` and count the work the booleans it runs on this thread do.
pub fn counted<R>(f: impl FnOnce() -> R) -> (R, Work) {
    COUNTED.set(Work::default());
    let out = f();
    (out, COUNTED.take())
}

/// Add to the count on this thread.
pub(crate) fn add(f: impl FnOnce(&mut Work)) {
    let mut work = COUNTED.get();
    f(&mut work);
    COUNTED.set(work);
}
