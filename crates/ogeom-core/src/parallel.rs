//! Deterministic parallelism for the kernel's embarrassingly parallel
//! stages.
//!
//! The rule that makes parallelism admissible here at all: **the answer must
//! be bit-identical at any thread count.** [`map_ordered`] guarantees it
//! structurally — each item is computed independently from shared read-only
//! input, results are collected in item order, and nothing about scheduling
//! can reach the output. A stage that cannot meet that bar stays sequential.
//!
//! The thread count comes from [`threads`]: the machine's parallelism by
//! default, overridable process-wide with [`set_threads`] — including down
//! to one, which is also what tiny workloads collapse to on their own.
//! Worker threads re-install the caller's progress watch, so cancellation
//! reaches into the workers.

use std::sync::atomic::{AtomicUsize, Ordering};

use crate::progress;

/// 0 means "ask the machine".
static THREADS: AtomicUsize = AtomicUsize::new(0);

/// The thread count parallel stages will use.
#[must_use]
pub fn threads() -> usize {
    let configured = THREADS.load(Ordering::Relaxed);
    if configured != 0 {
        return configured;
    }
    std::thread::available_parallelism().map_or(1, std::num::NonZero::get)
}

/// Set the process-wide thread count for parallel stages. `0` restores the
/// machine default. The answer never depends on this — only the wall clock
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
    if workers <= 1 || items.len() <= 1 {
        return items.iter().enumerate().map(|(i, t)| f(i, t)).collect();
    }

    let snapshot = progress::snapshot();
    let chunk = items.len().div_ceil(workers);
    let mut results: Vec<Vec<R>> = Vec::with_capacity(workers);
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(workers);
        for (w, slice) in items.chunks(chunk).enumerate() {
            let start = w * chunk;
            let f = &f;
            let snapshot = snapshot.clone();
            handles.push(scope.spawn(move || {
                progress::with_snapshot(snapshot.as_ref(), || {
                    slice
                        .iter()
                        .enumerate()
                        .map(|(i, t)| f(start + i, t))
                        .collect::<Vec<R>>()
                })
            }));
        }
        for handle in handles {
            match handle.join() {
                Ok(part) => results.push(part),
                Err(panic) => std::panic::resume_unwind(panic),
            }
        }
    });
    results.into_iter().flatten().collect()
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
