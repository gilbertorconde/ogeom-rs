//! A face mesh cache's eviction order, sweep and invalidation; its byte cap
//! is measured against the heap in `face_mesh_cache_bytes.rs`.
#![allow(clippy::unwrap_used, reason = "test code")]

use ogeom::algo::{make_box, make_cylinder};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::{Deflection, FaceMeshCache, face_meshes_for};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A block drilled through by a hole of `radius`, in a model of its own:
/// each radius gives a bore and two pierced faces of its own, and the
/// other four faces recur.
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

fn faces(model: &Model, solid: &Shape) -> Vec<Shape> {
    explore_unique(model, solid, ShapeType::Face).unwrap()
}

/// Every face of `solid` meshed through `cache` against a fresh mesh of
/// the same face: a hit must be the mesh drawing it again gives.
fn assert_meshed_alike(model: &Model, solid: &Shape, cache: &FaceMeshCache) {
    let d = Deflection::default();
    let (through, _) = face_meshes_for(model, solid, d, Some(cache), T).unwrap();
    let (fresh, _) = face_meshes_for(model, solid, d, None, T).unwrap();
    assert_eq!(through.len(), fresh.len());
    for ((_, a), (_, b)) in through.into_iter().zip(fresh) {
        assert_eq!(a.unwrap(), b.unwrap());
    }
}

/// A capped cache drops the least recently used face, a lookup counting as
/// a use, and every mesh it hands back is the one drawing the face again
/// gives.
#[test]
fn eviction_drops_the_least_recently_used_face() {
    let (model, solid) = drilled(4.0);
    let all = faces(&model, &solid);
    let d = Deflection::with_chord(0.05).unwrap();
    let cache = FaceMeshCache::new().with_face_limit(3);
    // Faces by index into `all`, and whether the lookup is a hit.
    let script = [
        (0, false),
        (1, false),
        (2, false),
        (0, true),  // order 1 2 0
        (3, false), // evicts 1: order 2 0 3
        (3, true),  // 2 0 3
        (0, true),  // 2 3 0
        (2, true),  // 3 0 2
        (1, false), // evicts 3: 0 2 1
        (3, false), // evicts 0: 2 1 3
        (2, true),
        (0, false), // evicts 1
        (1, false), // evicts 3
    ];
    for (step, &(index, hit)) in script.iter().enumerate() {
        let before = cache.stats();
        let mesh = ogeom::mesh::triangulate_face_kept(&model, &all[index], d, &cache, T).unwrap();
        let after = cache.stats();
        assert_eq!(after.hits - before.hits, u64::from(hit), "step {step}");
        assert_eq!(
            mesh,
            ogeom::mesh::triangulate_face(&model, &all[index], d, T).unwrap()
        );
        assert!(cache.len() <= 3);
    }
    assert_eq!(cache.stats().evicted, 5);
}

/// Swept after each body, the cache holds the faces of the body meshed
/// last, and meshing it again takes every face from the cache.
#[test]
fn retain_used_keeps_the_faces_of_the_body_meshed_last() {
    let d = Deflection::default();
    let cache = FaceMeshCache::new();
    let (first, a) = drilled(4.0);
    face_meshes_for(&first, &a, d, Some(&cache), T).unwrap();
    cache.retain_used();
    assert_eq!(cache.len(), faces(&first, &a).len());

    let (second, b) = drilled(4.25);
    face_meshes_for(&second, &b, d, Some(&cache), T).unwrap();
    // The bore and the two faces it pierces are new; the old three stay
    // until the sweep.
    assert_eq!(cache.len(), faces(&first, &a).len() + 3);
    cache.retain_used();
    assert_eq!(cache.len(), faces(&second, &b).len());
    assert_eq!(cache.stats().dropped, 3);

    let misses = cache.stats().misses;
    assert_meshed_alike(&second, &b, &cache);
    assert_eq!(cache.stats().misses, misses);
}

/// An edited face never matches its old mesh, and a forgotten face is
/// drawn again; either way what comes back is a fresh mesh of the face.
#[test]
fn edited_and_forgotten_faces_are_meshed_again() {
    let cache = FaceMeshCache::new();
    let (model, solid) = drilled(4.0);
    assert_meshed_alike(&model, &solid, &cache);
    let held = cache.len();

    // The edit: a wider bore. Its three faces miss, the rest hit.
    let (edited, wider) = drilled(4.5);
    let before = cache.stats();
    assert_meshed_alike(&edited, &wider, &cache);
    let after = cache.stats();
    assert_eq!(after.faces, held + 3);
    assert_eq!(after.misses - before.misses, 3);

    let old: Vec<_> = faces(&model, &solid)
        .iter()
        .map(|f| ogeom::topo::face_content(&model, f))
        .collect();
    let bore = faces(&edited, &wider)
        .into_iter()
        .find(|f| !old.contains(&ogeom::topo::face_content(&edited, f)))
        .unwrap();
    assert!(cache.forget(&edited, &bore));
    assert!(!cache.forget(&edited, &bore));
    let before = cache.stats();
    assert_meshed_alike(&edited, &wider, &cache);
    assert_eq!(cache.stats().misses - before.misses, 1);

    cache.clear();
    assert!(cache.is_empty());
    assert_eq!(cache.bytes(), 0);
}
