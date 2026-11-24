//! A face mesh cache's byte cap against the heap the allocator sees it free.
//! Alone in its binary: the allocator's counts are process-wide, and no
//! other test may allocate while one is read.
#![allow(clippy::unwrap_used, reason = "test code")]

use std::alloc::System;

use ogeom::algo::{make_box, make_cylinder};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::{Deflection, FaceMeshCache, face_meshes_for};
use ogeom::topo::{Model, Shape};
use stats_alloc::{INSTRUMENTED_SYSTEM, Region, StatsAlloc};

#[global_allocator]
static GLOBAL: &StatsAlloc<System> = &INSTRUMENTED_SYSTEM;

const T: Tolerances = Tolerances::millimetres();

/// The heap `cache` held, as the allocator sees it freed when it drops.
fn freed_by_dropping(cache: FaceMeshCache) -> usize {
    let region = Region::new(GLOBAL);
    drop(cache);
    region.change().bytes_deallocated
}

/// A block drilled through by a hole of `radius`, in a model of its own.
fn drilled(radius: f64) -> (Model, Shape) {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (30.0, 20.0, 10.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(12.0, 9.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let hole = make_cylinder(&mut model, at, radius, 12.0, T)
        .unwrap()
        .shape;
    let solid = ogeom::boolean::cut(&mut model, &block, &hole, T)
        .unwrap()
        .shape;
    (model, solid)
}

/// A cache capped at a number of bytes never holds more, by its own count
/// after every store and by the allocator when it drops; and uncapped, its
/// count is never under what the allocator frees, so the cap is a bound on
/// the real heap.
#[test]
fn a_byte_capped_cache_never_holds_more_than_its_cap() {
    let radii: Vec<f64> = (0..24).map(|i| 3.0 + 0.05 * f64::from(i)).collect();
    let bodies: Vec<(Model, Shape)> = radii.iter().map(|&r| drilled(r)).collect();
    let d = Deflection::default();

    let open = FaceMeshCache::new();
    for (model, solid) in &bodies {
        face_meshes_for(model, solid, d, Some(&open), T).unwrap();
    }
    let counted = open.bytes();
    let held = freed_by_dropping(open);
    eprintln!("uncapped: counted {counted} B, allocator freed {held} B");
    assert!(held <= counted, "freed {held} B, counted {counted} B");

    for limit in [counted / 2, counted / 5, counted / 20, 2_000] {
        let capped = FaceMeshCache::new().with_byte_limit(limit);
        for (model, solid) in &bodies {
            face_meshes_for(model, solid, d, Some(&capped), T).unwrap();
            assert!(capped.bytes() <= limit);
        }
        let stats = capped.stats();
        assert!(stats.evicted > 0);
        let held = freed_by_dropping(capped);
        eprintln!(
            "cap {limit} B: counted {} B over {} faces, allocator freed {held} B, {} evicted",
            stats.bytes, stats.faces, stats.evicted
        );
        assert!(held <= limit, "freed {held} B under a cap of {limit} B");
    }
}
