//! Deterministic parallelism for the kernel's embarrassingly parallel
//! stages.
//!
//! The rule that makes parallelism admissible here at all: **the answer must
//! be bit-identical at any thread count.** [`map_ordered`] guarantees it
//! structurally: each item is computed independently from shared read-only
//! input, results are collected in item order, and nothing about scheduling
//! can reach the output. A stage that cannot meet that bar stays sequential.
//!
//! The thread count comes from [`threads`]: the machine's parallelism by
//! default, or the `OGEOM_THREADS` environment variable when it holds a
//! positive count, and overridable process-wide with [`set_threads`],
//! including down to one, which is also what tiny workloads collapse to on
//! their own.
//! Worker threads re-install the caller's progress watch, so cancellation
//! reaches into the workers.
//!
//! Workers are started with `std::thread::scope`. On a target whose
//! standard library cannot start a thread, a host that has threads there (a
//! browser's web workers sharing the module's memory) lends them with
//! [`set_pool`]; once set, every stage runs its workers on that [`Pool`].

use std::any::Any;
use std::panic::AssertUnwindSafe;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Mutex, OnceLock, PoisonError};

use crate::progress;

/// Whether the standard library can start a thread on this target. On
/// `wasm32-unknown-unknown`, and on WebAssembly without shared memory, a
/// spawn panics, so without a [`Pool`] every stage runs on the calling
/// thread there whatever [`set_threads`] asked for.
const CAN_SPAWN: bool = !cfg!(any(
    all(target_arch = "wasm32", target_os = "unknown"),
    all(target_family = "wasm", not(target_feature = "atomics")),
));

/// 0 means "ask the machine".
static THREADS: AtomicUsize = AtomicUsize::new(0);

/// Threads a host lends the kernel, where `std::thread` cannot start one.
///
/// Once set with [`set_pool`], every parallel stage runs its workers through
/// [`Pool::broadcast`] instead of starting threads of its own.
///
/// The pool's workers must not need the thread that calls into the kernel
/// to make progress: that thread blocks inside [`Pool::broadcast`] until
/// every copy has returned. A browser's main thread cannot block, so a
/// browser host calls the kernel from a worker of its own.
///
/// A rayon-backed host:
///
/// ```text
/// struct Rayon;
/// impl ogeom_core::parallel::Pool for Rayon {
///     fn workers(&self) -> usize {
///         rayon::current_num_threads()
///     }
///     fn broadcast(&self, copies: usize, job: &(dyn Fn() + Sync)) {
///         rayon::scope(|s| {
///             for _ in 0..copies {
///                 s.spawn(|_| job());
///             }
///         });
///     }
/// }
/// ```
pub trait Pool: Sync {
    /// How many workers the pool runs at once.
    fn workers(&self) -> usize;

    /// Run `job` on `copies` workers at once and return when every copy has
    /// returned.
    ///
    /// Each copy takes items until none are left, so a copy that starts late
    /// or never costs time, not results: whatever the copies leave, the
    /// caller computes after `broadcast` returns. A copy does not unwind
    /// into the pool: a panic in a stage is caught in the copy and resumed
    /// on the caller once `broadcast` returns.
    fn broadcast(&self, copies: usize, job: &(dyn Fn() + Sync));
}

static POOL: OnceLock<&'static dyn Pool> = OnceLock::new();

/// Run parallel stages on `pool` from now on. The first call wins; later
/// calls change nothing.
///
/// With a pool set, [`threads`] honours [`set_threads`] and `OGEOM_THREADS`
/// on every target, and a stage takes at most [`Pool::workers`] workers.
pub fn set_pool(pool: &'static dyn Pool) {
    let _first_wins = POOL.set(pool);
}

/// The thread count parallel stages will use.
///
/// The count given to [`set_threads`] when there is one; otherwise
/// `OGEOM_THREADS` from the environment when it parses as a positive
/// count; otherwise the machine's available parallelism, or the lent
/// [`Pool`]'s worker count on a target that cannot start a thread. The
/// environment is read once, on the first call that needs it. One on a
/// target that cannot start a thread and has no pool.
#[must_use]
pub fn threads() -> usize {
    let pool = POOL.get();
    if !CAN_SPAWN && pool.is_none() {
        return 1;
    }
    let configured = THREADS.load(Ordering::Relaxed);
    if configured != 0 {
        return configured;
    }
    static ENVIRONMENT: OnceLock<Option<usize>> = OnceLock::new();
    let environment = *ENVIRONMENT
        .get_or_init(|| from_environment(std::env::var("OGEOM_THREADS").ok().as_deref()));
    if let Some(count) = environment {
        return count;
    }
    match pool {
        Some(pool) if !CAN_SPAWN => pool.workers().max(1),
        _ => {
            // Asked once: the query reads the scheduler's affinity and
            // quota, and every parallel stage asks.
            static MACHINE: OnceLock<usize> = OnceLock::new();
            *MACHINE.get_or_init(|| {
                std::thread::available_parallelism().map_or(1, std::num::NonZero::get)
            })
        }
    }
}

/// The thread count an `OGEOM_THREADS` value asks for: a positive integer,
/// surrounding blanks allowed. Anything else asks for nothing.
fn from_environment(value: Option<&str>) -> Option<usize> {
    value?.trim().parse::<usize>().ok().filter(|&n| n > 0)
}

std::thread_local! {
    /// Set on a worker for the life of its stage: a parallel stage inside
    /// another (a boolean's face pairs each classifying points in parallel)
    /// runs on the worker it lands on, rather than starting a machine's
    /// worth of workers per outer item.
    static INSIDE: core::cell::Cell<bool> = const { core::cell::Cell::new(false) };
}

/// Set the process-wide thread count for parallel stages. `0` restores the
/// machine default. The answer never depends on this. Only the wall clock
/// does.
pub fn set_threads(count: usize) {
    THREADS.store(count, Ordering::Relaxed);
}

/// The workers `std::thread::scope` starts where no [`Pool`] is lent.
struct Scoped;

impl Pool for Scoped {
    fn workers(&self) -> usize {
        usize::MAX
    }

    fn broadcast(&self, copies: usize, job: &(dyn Fn() + Sync)) {
        std::thread::scope(|scope| {
            for _ in 0..copies {
                scope.spawn(job);
            }
        });
    }
}

/// Map `f` over `items` on up to [`threads`] workers, returning results in
/// item order. `f` receives the item index and the item.
///
/// The workers are scoped threads, or the lent [`Pool`]'s when one is set,
/// capped at its [`Pool::workers`]. A call made from a worker of another
/// stage runs on that worker.
///
/// Determinism holds by construction: items are computed independently and
/// results placed by index, so the output is identical at any thread count,
/// pool or not. The caller's progress watch is re-installed in every
/// worker; `f` may checkpoint through it. A panic in `f` resumes on the
/// caller.
pub fn map_ordered<T, R>(items: &[T], f: impl Fn(usize, &T) -> R + Sync) -> Vec<R>
where
    T: Sync,
    R: Send,
{
    let lent = POOL.get().copied();
    let pool: &dyn Pool = lent.unwrap_or(&Scoped);
    let workers = threads().min(pool.workers()).clamp(1, items.len().max(1));
    if (!CAN_SPAWN && lent.is_none())
        || workers <= 1
        || items.len() <= 1
        || INSIDE.with(core::cell::Cell::get)
    {
        return items.iter().enumerate().map(|(i, t)| f(i, t)).collect();
    }

    let snapshot = progress::snapshot();
    // Work is *taken*, not dealt: expensive items cluster (one spline-heavy
    // face's edges sit adjacent in a reader's job list), and a worker dealt
    // that region as a contiguous chunk finishes last while the rest idle.
    // Each worker pulls the next undone index instead, so the wall clock
    // tracks the total work rather than the heaviest deal. The answer cannot
    // tell the difference: every index is computed by the same call exactly
    // once, and the merge reassembles by index, so the output is the item
    // order however the indices were claimed.
    let next = AtomicUsize::new(0);
    let take = || {
        let mut mine = Vec::new();
        loop {
            let i = next.fetch_add(1, Ordering::Relaxed);
            let Some(item) = items.get(i) else { break };
            mine.push((i, f(i, item)));
        }
        mine
    };
    let parts: Mutex<Vec<Vec<(usize, R)>>> = Mutex::new(Vec::with_capacity(workers));
    let panicked: Mutex<Option<Box<dyn Any + Send>>> = Mutex::new(None);
    let job = || {
        let outer = INSIDE.with(|inside| inside.replace(true));
        let run = std::panic::catch_unwind(AssertUnwindSafe(|| {
            progress::with_snapshot(snapshot.as_ref(), take)
        }));
        INSIDE.with(|inside| inside.set(outer));
        match run {
            Ok(part) => parts
                .lock()
                .unwrap_or_else(PoisonError::into_inner)
                .push(part),
            Err(panic) => {
                // The other copies stop at their next claim.
                next.store(items.len(), Ordering::Relaxed);
                panicked
                    .lock()
                    .unwrap_or_else(PoisonError::into_inner)
                    .get_or_insert(panic);
            }
        }
    };
    pool.broadcast(workers, &job);
    if let Some(panic) = panicked
        .into_inner()
        .unwrap_or_else(PoisonError::into_inner)
    {
        std::panic::resume_unwind(panic);
    }
    let mut indexed: Vec<(usize, R)> = parts
        .into_inner()
        .unwrap_or_else(PoisonError::into_inner)
        .into_iter()
        .flatten()
        .collect();
    // Items no copy claimed (a pool that ran fewer copies than asked) are
    // taken here, under the caller's own watch.
    indexed.extend(take());
    assert_eq!(
        indexed.len(),
        items.len(),
        "a pool returned from broadcast while a copy was still running"
    );
    indexed.sort_unstable_by_key(|(i, _)| *i);
    indexed.into_iter().map(|(_, r)| r).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn order_is_item_order_at_any_thread_count() {
        let items: Vec<usize> = (0..137).collect();
        let serial: Vec<usize> = items.iter().map(|x| x * 3).collect();
        for count in [1, 2, 7] {
            set_threads(count);
            let parallel = map_ordered(&items, |i, x| {
                assert_eq!(i, *x);
                x * 3
            });
            assert_eq!(parallel, serial);
        }
        set_threads(0);
    }

    #[test]
    #[should_panic(expected = "item 50")]
    fn a_panic_on_a_worker_resumes_on_the_caller() {
        let items: Vec<usize> = (0..137).collect();
        let _ = map_ordered(&items, |i, _| {
            assert!(i != 50, "item 50");
            i
        });
    }

    #[test]
    fn the_environment_count_must_be_a_positive_integer() {
        assert_eq!(from_environment(Some(" 3 ")), Some(3));
        for refused in [
            None,
            Some(""),
            Some("0"),
            Some("-2"),
            Some("four"),
            Some("2.5"),
        ] {
            assert_eq!(from_environment(refused), None, "{refused:?}");
        }
    }

    #[test]
    fn cancellation_reaches_the_workers() {
        let watch = progress::Watch::new();
        watch.canceller().cancel();
        set_threads(4);
        let items: Vec<usize> = (0..64).collect();
        let outcomes = progress::watched(&watch, || {
            map_ordered(&items, |_, _| progress::checkpoint())
        });
        set_threads(0);
        assert!(outcomes.iter().all(std::result::Result::is_err));
    }
}
