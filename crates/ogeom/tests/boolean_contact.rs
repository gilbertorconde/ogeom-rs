//! The contact configurations the deferred table still called refusals,
//! pinned as the working behaviour they have become: curved same-domain
//! pairs unify, and contact confined to an edge or a vertex passes through
//! the boolean without harm.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point};
use ogeom::mesh::Deflection;
use ogeom::topo::Model;

const T: Tolerances = Tolerances::millimetres();

fn fine() -> Deflection {
    Deflection::with_chord(1e-3).unwrap()
}

fn volume(model: &Model, shape: &ogeom::topo::Shape) -> f64 {
    ogeom::algo::volume_properties(model, shape, fine(), T)
        .unwrap()
        .mass
}

/// Coaxial cylinders sharing one surface: flush stack, partial overlap,
/// and the cut — the curved same-domain family, held to closed forms.
#[test]
fn curved_same_domain_pairs_unify() {
    let pi = core::f64::consts::PI;

    // Flush stack: walls meet rim to rim on one infinite cylinder.
    {
        let mut model = Model::new();
        let bottom = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T)
            .unwrap()
            .shape;
        let up = Frame::new(Point::new(0.0, 0.0, 10.0), Direction::Z, Direction::X, T).unwrap();
        let top = ogeom::algo::make_cylinder(&mut model, up, 5.0, 10.0, T)
            .unwrap()
            .shape;
        let fused = ogeom::boolean::fuse(&mut model, &bottom, &top, T)
            .unwrap()
            .shape;
        let v = volume(&model, &fused);
        assert!(
            (v - pi * 25.0 * 20.0).abs() / (pi * 25.0 * 20.0) < 1e-3,
            "the stack is one drum: {v}"
        );
    }

    // Partial overlap: the walls overlap in a band and split each other at
    // the other's rims.
    {
        let mut model = Model::new();
        let bottom = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T)
            .unwrap()
            .shape;
        let up = Frame::new(Point::new(0.0, 0.0, 5.0), Direction::Z, Direction::X, T).unwrap();
        let top = ogeom::algo::make_cylinder(&mut model, up, 5.0, 10.0, T)
            .unwrap()
            .shape;
        let fused = ogeom::boolean::fuse(&mut model, &bottom, &top, T)
            .unwrap()
            .shape;
        let v = volume(&model, &fused);
        assert!(
            (v - pi * 25.0 * 15.0).abs() / (pi * 25.0 * 15.0) < 1e-3,
            "the overlap fuses to one taller drum: {v}"
        );
    }
    {
        let mut model = Model::new();
        let bottom = ogeom::algo::make_cylinder(&mut model, Frame::WORLD, 5.0, 10.0, T)
            .unwrap()
            .shape;
        let up = Frame::new(Point::new(0.0, 0.0, 5.0), Direction::Z, Direction::X, T).unwrap();
        let top = ogeom::algo::make_cylinder(&mut model, up, 5.0, 10.0, T)
            .unwrap()
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &bottom, &top, T)
            .unwrap()
            .shape;
        let v = volume(&model, &cut);
        assert!(
            (v - pi * 25.0 * 5.0).abs() / (pi * 25.0 * 5.0) < 1e-3,
            "the cut keeps the un-overlapped stub: {v}"
        );
    }
}

/// Contact confined to one edge or one vertex: the fuse of edge-touching
/// boxes is one valid solid of both volumes — the shared edge is the
/// non-manifold seam the model permits — and cutting a corner-touching
/// tool removes nothing.
#[test]
fn edge_and_vertex_contact_pass_through_the_boolean() {
    {
        let mut model = Model::new();
        let a = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let f = Frame::new(Point::new(10.0, 0.0, 10.0), Direction::Z, Direction::X, T).unwrap();
        let b = ogeom::algo::make_box(&mut model, f, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let fused = ogeom::boolean::fuse(&mut model, &a, &b, T).unwrap().shape;
        let v = volume(&model, &fused);
        assert!((v - 2000.0).abs() < 1e-6, "both boxes survive: {v}");
        assert!(ogeom::algo::check(&model, &fused, T).unwrap().is_valid());
    }
    {
        let mut model = Model::new();
        let a = ogeom::algo::make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let f = Frame::new(Point::new(10.0, 10.0, 10.0), Direction::Z, Direction::X, T).unwrap();
        let b = ogeom::algo::make_box(&mut model, f, (10.0, 10.0, 10.0), T)
            .unwrap()
            .shape;
        let cut = ogeom::boolean::cut(&mut model, &a, &b, T).unwrap().shape;
        let v = volume(&model, &cut);
        assert!(
            (v - 1000.0).abs() < 1e-6,
            "a corner touch removes nothing: {v}"
        );
    }
}
