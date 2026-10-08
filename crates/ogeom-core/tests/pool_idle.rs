//! A lent pool that returns without running the copies it was asked for
//! leaves its items to the caller. `set_pool` is process-wide and the first
//! call wins, so this binary holds the only test that sets one.

use ogeom_core::parallel::{Pool, map_ordered, set_pool, set_threads};

/// A pool that runs nothing.
struct Idle;

impl Pool for Idle {
    fn workers(&self) -> usize {
        4
    }

    fn broadcast(&self, _copies: usize, _job: &(dyn Fn() + Sync)) {}
}

#[test]
fn items_a_pool_never_ran_are_still_mapped() {
    set_pool(&Idle);
    set_threads(4);
    let items: Vec<usize> = (0..37).collect();
    let out = map_ordered(&items, |i, x| {
        assert_eq!(i, *x);
        x * 3
    });
    assert_eq!(out, items.iter().map(|x| x * 3).collect::<Vec<_>>());
}
