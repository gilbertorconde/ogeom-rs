//! A box with every edge filleted meshes to the same triangles on threads a
//! host lends as on one thread.
//!
//! Natively the test pins the single-threaded triangle count. In a
//! shared-memory browser build (`+atomics,+bulk-memory`, `-Z build-std`) it
//! runs in a dedicated worker, lends the kernel a `wasm-bindgen-rayon` pool
//! of web workers, and asserts that the mesh made on that pool is
//! bit-identical to the one made on the calling worker alone, that more
//! than one web worker took part, and that the triangle count is the native
//! one. See `src/lib.rs` for the commands.
#![cfg(any(not(target_arch = "wasm32"), target_feature = "atomics"))]
#![allow(clippy::unwrap_used, reason = "a test")]

use ogeom::algo::make_box;
use ogeom::core::Tolerances;
use ogeom::fillet::fillet_edges;
use ogeom::math::Frame;
use ogeom::mesh::{Deflection, triangulate};
use ogeom::topo::{Model, ShapeType, Triangulation, explore_unique};

/// The triangles the filleted box meshes to natively on one thread.
const TRIANGLES: usize = 516;

fn filleted_box_mesh() -> Triangulation {
    let tol = Tolerances::default();
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), tol)
        .unwrap()
        .shape;
    let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
    let rolled = fillet_edges(&mut model, &block, &edges, 1.0, tol).unwrap();
    triangulate(&model, &rolled.shape, Deflection::default(), tol).unwrap()
}

#[cfg(not(target_arch = "wasm32"))]
#[test]
fn a_filleted_box_meshes_natively() {
    ogeom::core::parallel::set_threads(1);
    assert_eq!(filleted_box_mesh().triangles.len(), TRIANGLES);
}

#[cfg(all(target_arch = "wasm32", target_feature = "atomics"))]
mod lent {
    use std::collections::BTreeSet;
    use std::sync::Mutex;

    use ogeom::core::parallel::{Pool, set_pool, set_threads};
    use wasm_bindgen_test::{wasm_bindgen_test, wasm_bindgen_test_configure};

    use super::{TRIANGLES, filleted_box_mesh};

    // A browser's main thread cannot block, and the kernel's caller blocks
    // while the pool runs a stage.
    wasm_bindgen_test_configure!(run_in_dedicated_worker);

    /// rayon's global pool, which `wasm-bindgen-rayon` backs with web
    /// workers, recording which of its threads ran a copy.
    struct Rayon {
        ran_on: Mutex<BTreeSet<usize>>,
    }

    impl Pool for Rayon {
        fn workers(&self) -> usize {
            rayon::current_num_threads()
        }

        fn broadcast(&self, copies: usize, job: &(dyn Fn() + Sync)) {
            rayon::scope(|s| {
                for _ in 0..copies {
                    s.spawn(|_| {
                        if let Some(index) = rayon::current_thread_index() {
                            self.ran_on.lock().unwrap().insert(index);
                        }
                        job();
                    });
                }
            });
        }
    }

    static RAYON: Rayon = Rayon {
        ran_on: Mutex::new(BTreeSet::new()),
    };

    #[wasm_bindgen_test]
    async fn a_filleted_box_meshes_alike_on_lent_web_workers() {
        wasm_bindgen_futures::JsFuture::from(wasm_bindgen_rayon::init_thread_pool(4))
            .await
            .unwrap();
        set_pool(&RAYON);

        set_threads(1);
        let alone = filleted_box_mesh();
        assert!(RAYON.ran_on.lock().unwrap().is_empty());

        set_threads(0);
        let lent = filleted_box_mesh();
        let workers = RAYON.ran_on.lock().unwrap().len();
        assert!(workers > 1, "the stages ran on {workers} web worker(s)");

        assert_eq!(alone.triangles.len(), TRIANGLES);
        assert_eq!(lent.triangles, alone.triangles);
        let bits = |mesh: &ogeom::topo::Triangulation| {
            mesh.positions
                .iter()
                .flat_map(|p| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()])
                .collect::<Vec<u64>>()
        };
        assert_eq!(bits(&lent), bits(&alone));
    }
}
