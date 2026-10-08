//! A lent pool runs the parallel stages. `set_pool` is process-wide and the
//! first call wins, so this binary holds the only test that sets one.

use std::sync::atomic::{AtomicUsize, Ordering};

use ogeom_core::parallel::{Pool, map_ordered, set_pool, set_threads};

/// A pool on scoped threads that counts the workers it was asked for.
struct Counting(AtomicUsize);

impl Pool for Counting {
    fn workers(&self) -> usize {
        4
    }

    fn broadcast(&self, copies: usize, job: &(dyn Fn() + Sync)) {
        self.0.fetch_add(copies, Ordering::Relaxed);
        std::thread::scope(|s| {
            for _ in 0..copies {
                s.spawn(job);
            }
        });
    }
}

static POOL: Counting = Counting(AtomicUsize::new(0));

#[test]
fn stages_run_on_the_pool_lent_and_keep_their_order() {
    set_pool(&POOL);
    set_threads(4);
    let items: Vec<usize> = (0..137).collect();
    let out = map_ordered(&items, |_, x| x * 3);
    assert_eq!(out, items.iter().map(|x| x * 3).collect::<Vec<_>>());
    assert_eq!(POOL.0.load(Ordering::Relaxed), 4);
}
