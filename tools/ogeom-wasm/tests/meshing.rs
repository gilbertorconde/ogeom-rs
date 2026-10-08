//! A box and a cylinder mesh to the same triangle counts natively and on
//! `wasm32-unknown-unknown`. Tessellation reads the kernel's clock and asks
//! for threads, both of which the browser target lacks.
#![allow(clippy::unwrap_used, reason = "a test")]

use ogeom::algo::{make_box, make_cylinder};
use ogeom::core::Tolerances;
use ogeom::math::Frame;
use ogeom::mesh::{Deflection, triangulate};
use ogeom::topo::Model;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen_test::wasm_bindgen_test as test;

#[test]
fn a_box_meshes() {
    let tol = Tolerances::default();
    let mut model = Model::new();
    let solid = make_box(&mut model, Frame::WORLD, (10.0, 20.0, 30.0), tol)
        .unwrap()
        .shape;
    let mesh = triangulate(&model, &solid, Deflection::default(), tol).unwrap();
    assert_eq!(mesh.triangles.len(), 12);
}

#[test]
fn a_cylinder_meshes() {
    let tol = Tolerances::default();
    let mut model = Model::new();
    let solid = make_cylinder(&mut model, Frame::WORLD, 5.0, 12.0, tol)
        .unwrap()
        .shape;
    let mesh = triangulate(&model, &solid, Deflection::default(), tol).unwrap();
    assert_eq!(mesh.triangles.len(), 90);
}
