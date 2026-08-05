//! Progress reporting and cancellation for long operations.
//!
//! A caller who starts a tessellation, a boolean or an import may need to
//! stop it — a user closed the dialog — or to show that it is alive. The
//! kernel's answer is a [`Watch`]: install one around a call with
//! [`watched`], hand its [`Canceller`] to whoever may pull the plug, and
//! every long loop inside the kernel calls [`checkpoint`] at its own
//! boundaries. A cancelled checkpoint returns
//! [`OgeomError::Cancelled`], which unwinds as
//! an ordinary error — no partial result pretends to be whole.
//!
//! The watch travels implicitly, by scope: operations keep their signatures,
//! and code that never installs a watch pays one thread-local read per
//! checkpoint. Worker threads a kernel operation spawns re-install the
//! caller's watch through [`snapshot`]/[`with_snapshot`], so cancellation
//! reaches into parallel stages too.
//!
//! Cancellation is cooperative and prompt rather than immediate: it lands at
//! the next checkpoint, and checkpoints sit at stage and item boundaries,
//! never inside an invariant-restoring section.

use std::cell::RefCell;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use crate::{OgeomError, OgeomResult};

/// A stage sink: hears each stage name as the operation reaches it.
type Sink = Arc<dyn Fn(&str) + Send + Sync>;

/// What a watch carries: the flag, and an optional stage sink.
#[derive(Clone)]
struct State {
    cancel: Arc<AtomicBool>,
    sink: Option<Sink>,
}

impl core::fmt::Debug for State {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("State")
            .field("cancelled", &self.cancel.load(Ordering::Relaxed))
            .field("has_sink", &self.sink.is_some())
            .finish()
    }
}

thread_local! {
    static ACTIVE: RefCell<Option<State>> = const { RefCell::new(None) };
}

/// A scope's progress watch: cancellation flag plus an optional sink that
/// receives stage names as the operation passes them.
#[derive(Debug)]
pub struct Watch {
    state: State,
}

impl Watch {
    /// A watch with no sink: cancellation only.
    #[must_use]
    pub fn new() -> Self {
        Self {
            state: State {
                cancel: Arc::new(AtomicBool::new(false)),
                sink: None,
            },
        }
    }

    /// A watch whose sink hears each stage name once as the operation
    /// reaches it. The sink runs on whichever thread reaches the stage.
    #[must_use]
    pub fn with_sink(sink: impl Fn(&str) + Send + Sync + 'static) -> Self {
        Self {
            state: State {
                cancel: Arc::new(AtomicBool::new(false)),
                sink: Some(Arc::new(sink)),
            },
        }
    }

    /// The handle that cancels this watch, cloneable and sendable to
    /// whatever owns the stop button.
    #[must_use]
    pub fn canceller(&self) -> Canceller {
        Canceller {
            cancel: Arc::clone(&self.state.cancel),
        }
    }
}

impl Default for Watch {
    fn default() -> Self {
        Self::new()
    }
}

/// The stop button: cancel from any thread, any number of times.
#[derive(Debug, Clone)]
pub struct Canceller {
    cancel: Arc<AtomicBool>,
}

impl Canceller {
    /// Request cancellation. Takes effect at the operation's next
    /// [`checkpoint`].
    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::SeqCst);
    }

    /// Whether cancellation has been requested.
    #[must_use]
    pub fn is_cancelled(&self) -> bool {
        self.cancel.load(Ordering::SeqCst)
    }
}

/// Restores the previously active state when the scope ends, panics
/// included.
struct Scope {
    previous: Option<State>,
}

impl Drop for Scope {
    fn drop(&mut self) {
        ACTIVE.with(|active| {
            *active.borrow_mut() = self.previous.take();
        });
    }
}

/// Run `f` with `watch` active on this thread: every [`checkpoint`] inside
/// answers to it. Scopes nest; the inner watch wins until it ends.
pub fn watched<T>(watch: &Watch, f: impl FnOnce() -> T) -> T {
    let previous = ACTIVE.with(|active| active.borrow_mut().replace(watch.state.clone()));
    let _scope = Scope { previous };
    f()
}

/// The point a long loop offers for cancellation. Free when no watch is
/// installed.
///
/// # Errors
///
/// [`OgeomError::Cancelled`] if the active
/// watch has been cancelled.
pub fn checkpoint() -> OgeomResult<()> {
    ACTIVE.with(|active| {
        if let Some(state) = active.borrow().as_ref()
            && state.cancel.load(Ordering::Relaxed)
        {
            return Err(OgeomError::Cancelled);
        }
        Ok(())
    })
}

/// Announce a stage boundary to the active watch's sink, if there is one.
pub fn stage(name: &str) {
    ACTIVE.with(|active| {
        if let Some(state) = active.borrow().as_ref()
            && let Some(sink) = &state.sink
        {
            sink(name);
        }
    });
}

/// The active watch, portable to a worker thread. `None` when unwatched.
#[must_use]
pub fn snapshot() -> Option<WatchSnapshot> {
    ACTIVE.with(|active| active.borrow().clone().map(|state| WatchSnapshot { state }))
}

/// Run `f` under a snapshot taken on another thread — how a parallel stage
/// keeps answering the caller's watch.
pub fn with_snapshot<T>(snapshot: Option<&WatchSnapshot>, f: impl FnOnce() -> T) -> T {
    match snapshot {
        Some(snap) => {
            let watch = Watch {
                state: snap.state.clone(),
            };
            watched(&watch, f)
        }
        None => f(),
    }
}

/// An opaque, cloneable capture of the active watch.
#[derive(Debug, Clone)]
pub struct WatchSnapshot {
    state: State,
}

#[cfg(test)]
#[allow(clippy::unwrap_used)]
mod tests {
    use super::*;

    #[test]
    fn unwatched_checkpoints_are_free_and_fine() {
        assert!(checkpoint().is_ok());
        stage("nothing listens");
    }

    #[test]
    fn a_cancelled_watch_stops_the_next_checkpoint() {
        let watch = Watch::new();
        let stop = watch.canceller();
        let result: OgeomResult<()> = watched(&watch, || {
            checkpoint()?;
            stop.cancel();
            checkpoint()?;
            unreachable!("the second checkpoint must refuse");
        });
        assert!(matches!(result, Err(OgeomError::Cancelled)));
        assert!(stop.is_cancelled());
        // The scope is over: this thread is unwatched again.
        assert!(checkpoint().is_ok());
    }

    #[test]
    fn stages_reach_the_sink_and_scopes_nest() {
        use std::sync::Mutex;
        let heard: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
        let record = Arc::clone(&heard);
        let outer = Watch::with_sink(move |name| record.lock().unwrap().push(name.to_owned()));
        let inner = Watch::new();
        watched(&outer, || {
            stage("first");
            watched(&inner, || {
                // The inner watch has no sink; its scope masks the outer.
                stage("masked");
            });
            stage("second");
        });
        assert_eq!(*heard.lock().unwrap(), ["first", "second"]);
    }

    #[test]
    fn snapshots_carry_the_watch_across_threads() {
        let watch = Watch::new();
        watch.canceller().cancel();
        let snap = watched(&watch, snapshot);
        let outcome: OgeomResult<()> =
            std::thread::spawn(move || with_snapshot(snap.as_ref(), checkpoint))
                .join()
                .unwrap();
        assert!(matches!(outcome, Err(OgeomError::Cancelled)));
    }
}
