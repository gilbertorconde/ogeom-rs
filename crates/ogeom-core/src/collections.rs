//! Hash maps and sets for the whole kernel.
//!
//! [`FastMap`] and [`FastSet`] are `hashbrown` tables hashed by foldhash
//! with a fixed seed. The seed is what matters: a map's iteration order is
//! then a function of its insertions alone, the same in every run and every
//! process, so a stage that walks a map cannot answer differently from one
//! run to the next. The std tables seed themselves randomly and are refused
//! by clippy (`disallowed-types` in `clippy.toml`) outside the few public
//! signatures that already promised them.
//!
//! A fixed seed gives up the defence against keys chosen to collide. The
//! kernel hashes its own ids and coordinates, not adversarial strings.
//!
//! The tables have no `new` or `with_capacity`: build them with
//! `FastMap::default()`, `collect()`, or
//! `FastMap::with_capacity_and_hasher(n, FastHasher::default())`.

pub use hashbrown::{hash_map, hash_set};

/// The hasher every [`FastMap`] and [`FastSet`] builds: foldhash, seeded
/// the same way in every process.
pub type FastHasher = foldhash::fast::FixedState;

/// A hash map whose iteration order depends only on its contents and the
/// order they were inserted in.
pub type FastMap<K, V> = hashbrown::HashMap<K, V, FastHasher>;

/// A hash set whose iteration order depends only on its contents and the
/// order they were inserted in.
pub type FastSet<T> = hashbrown::HashSet<T, FastHasher>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iteration_order_is_fixed_by_the_insertions() {
        let build = || {
            let mut map = FastMap::default();
            for i in 0..1000_u64 {
                map.insert(i.wrapping_mul(0x9e37_79b9_7f4a_7c15), i);
            }
            map.into_iter().collect::<Vec<_>>()
        };
        assert_eq!(build(), build());
    }
}
