//! Counts of the work a boolean does, for the workspace's own tests: a
//! measure of what a boolean reaches that, unlike its time, does not move
//! with the machine it runs on. Present only with the `test-hooks` feature.

use std::cell::Cell;

/// What the booleans run on one thread did, summed.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Work {
    /// Faces gathered into the general fuse, both operands', each time one
    /// is run.
    pub faces_gathered: usize,
    /// Pairs of faces whose boxes meet, each asked for its sections.
    pub face_pairs: usize,
    /// Pieces the gathered faces were split into and classified.
    pub pieces: usize,
    /// Cuts made by trimming the faces a tool's wedges sit on, with no
    /// general fuse.
    pub seated_cuts: usize,
    /// Unions made by trimming the faces a tool's wedges sit on, with no
    /// general fuse.
    pub seated_fills: usize,
    /// Attempts made with contacts ended where their curves part by the
    /// radius they hold to, the last a refused boolean is retried with.
    pub stated_retries: usize,
    /// The widest tolerance any section between the two operands states,
    /// over every general fuse run.
    pub widest_section: f64,
}

thread_local! {
    static COUNTED: Cell<Work> = const {
        Cell::new(Work {
            faces_gathered: 0,
            face_pairs: 0,
            pieces: 0,
            seated_cuts: 0,
            seated_fills: 0,
            stated_retries: 0,
            widest_section: 0.0,
        })
    };
    static UNSEATED: Cell<bool> = const { Cell::new(false) };
    static UNTURNED: Cell<bool> = const { Cell::new(false) };
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

/// Run `f` with every cut [`cut_seated`](crate::cut_seated) is asked for
/// made by the general fuse, as [`cut`](crate::cut) makes it, and every
/// union [`fuse_seated`](crate::fuse_seated) is asked for as
/// [`fuse`](crate::fuse) makes it.
pub fn unseated<R>(f: impl FnOnce() -> R) -> R {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            UNSEATED.set(self.0);
        }
    }
    let _restore = Restore(UNSEATED.replace(true));
    f()
}

/// Whether the seated cuts and unions on this thread are left to the
/// general fuse.
pub(crate) fn is_unseated() -> bool {
    UNSEATED.get()
}

/// Run `f` with every boolean on this thread taking its operands on the
/// charts they come with: no ball or sphere face is restated on a chart
/// whose poles stand clear of the other operand.
pub fn unturned<R>(f: impl FnOnce() -> R) -> R {
    struct Restore(bool);
    impl Drop for Restore {
        fn drop(&mut self) {
            UNTURNED.set(self.0);
        }
    }
    let _restore = Restore(UNTURNED.replace(true));
    f()
}

/// Whether the booleans on this thread take their operands on the charts
/// they come with.
pub(crate) fn is_unturned() -> bool {
    UNTURNED.get()
}
