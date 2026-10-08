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

use std::sync::atomic::{AtomicUsize, Ordering};

use crate::progress;

/// Whether the standard library can start a thread on this target. On
/// `wasm32-unknown-unknown`, and on WebAssembly without shared memory, a
/// spawn panics, so every stage runs on the calling thread there whatever
/// [`set_threads`] asked for.
const CAN_SPAWN: bool = !cfg!(any(
    all(target_arch = "wasm32", target_os = "unknown"),
    all(target_family = "wasm", not(target_feature = "atomics")),
));

/// 0 means "ask the machine".
static THREADS: AtomicUsize = AtomicUsize::new(0);

/// The thread count parallel stages will use.
///
/// The count given to [`set_threads`] when there is one; otherwise
/// `OGEOM_THREADS` from the environment when it parses as a positive
/// count; otherwise the machine's available parallelism. The environment
/// is read once, on the first call that needs it. One on a target that
/// cannot start a thread.
#[must_use]
pub fn threads() -> usize {
    if !CAN_SPAWN {
        return 1;
    }
    let configured = THREADS.load(Ordering::Relaxed);
    if configured != 0 {
        return configured;
    }
    // Asked once: the query reads the scheduler's affinity and quota, and
    // every parallel stage asks.
    static MACHINE: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *MACHINE.get_or_init(|| {
        from_environment(std::env::var("OGEOM_THREADS").ok().as_deref()).unwrap_or_else(|| {
            std::thread::available_parallelism().map_or(1, std::num::NonZero::get)
        })
    })
}

/// The thread count an `OGEOM_THREADS` value asks for: a positive integer,
/// surrounding blanks allowed. Anything else asks for nothing.
fn from_environment(value: Option<&str>) -> Option<usize> {
    value?.trim().parse::<usize>().ok().filter(|&n| n > 0)
}

std::thread_local! {
    /// Set on a worker thread for the life of its stage: a parallel stage
    /// inside another (a boolean's face pairs each classifying points in
    /// parallel) runs on the worker it lands on, rather than spawning a
    /// machine's worth of threads per outer item.
    static INSIDE: core::cell::Cell<bool> = const { core::cell::Cell::new(false) };
}

/// Set the process-wide thread count for parallel stages. `0` restores the
/// machine default. The answer never depends on this. Only the wall clock
/// does.
pub fn set_threads(count: usize) {
    THREADS.store(count, Ordering::Relaxed);
}

/// Map `f` over `items` on up to [`threads`] scoped threads, returning
/// results in item order. `f` receives the item index and the item.
///
/// Determinism holds by construction: items are computed independently and
/// results placed by index, so the output is identical at any thread count.
/// The caller's progress watch is re-installed in every worker; `f` may
/// checkpoint through it.
pub fn map_ordered<T, R>(items: &[T], f: impl Fn(usize, &T) -> R + Sync) -> Vec<R>
where
    T: Sync,
    R: Send,
{
    let workers = threads().clamp(1, items.len().max(1));
    if !CAN_SPAWN || workers <= 1 || items.len() <= 1 || INSIDE.with(core::cell::Cell::get) {
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
    let next = std::sync::atomic::AtomicUsize::new(0);
    let mut parts: Vec<Vec<(usize, R)>> = Vec::with_capacity(workers);
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(workers);
        for _ in 0..workers {
            let f = &f;
            let next = &next;
            let snapshot = snapshot.clone();
            handles.push(scope.spawn(move || {
                INSIDE.with(|inside| inside.set(true));
                progress::with_snapshot(snapshot.as_ref(), || {
                    let mut mine = Vec::new();
                    loop {
                        let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        let Some(item) = items.get(i) else { break };
                        mine.push((i, f(i, item)));
                    }
                    mine
                })
            }));
        }
        for handle in handles {
            match handle.join() {
                Ok(part) => parts.push(part),
                Err(panic) => std::panic::resume_unwind(panic),
            }
        }
    });
    let mut indexed: Vec<(usize, R)> = parts.into_iter().flatten().collect();
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
