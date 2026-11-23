//! The union of many solids in one call, against closed-form volumes,
//! against the same solids fused two at a time, and across reorderings.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use std::f64::consts::{PI, TAU};

use ogeom::algo::{
    check, make_box, make_cylinder, make_face, make_polygon, make_prism, make_revolution,
    volume_properties,
};
use ogeom::boolean::{fuse, fuse_all};
use ogeom::core::{OgeomError, Tolerances};
use ogeom::geom::PlaneSurface;
use ogeom::heal::unify_same_domain;
use ogeom::math::{Axis, Direction, Frame, Plane, Point, Vector};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn cube(model: &mut Model, x: f64, y: f64, z: f64, size: f64) -> Shape {
    make_box(model, at(x, y, z), (size, size, size), T)
        .unwrap()
        .shape
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    let fine = ogeom::mesh::Deflection {
        chord: 1e-4,
        ..ogeom::mesh::Deflection::default()
    };
    volume_properties(model, shape, fine, T).unwrap().mass
}

fn count(model: &Model, shape: &Shape, kind: ShapeType) -> usize {
    explore_unique(model, shape, kind).unwrap().len()
}

fn assert_valid(model: &Model, shape: &Shape) {
    let diagnosis = check(model, shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
}

/// The solids fused left to right, two at a time.
fn chained(model: &mut Model, solids: &[Shape]) -> Shape {
    let mut acc = solids[0].clone();
    for solid in &solids[1..] {
        acc = fuse(model, &acc, solid, T).unwrap().shape;
    }
    acc
}

/// A fixed shuffle of `0..n`, different for each `seed`.
fn shuffled(n: usize, seed: u64) -> Vec<usize> {
    let mut order: Vec<usize> = (0..n).collect();
    let mut state = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
    for i in (1..n).rev() {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let j = usize::try_from((state >> 33) % (i as u64 + 1)).unwrap();
        order.swap(i, j);
    }
    order
}

#[test]
fn fifty_cubes_in_a_row_fuse_into_one_box_of_six_faces() {
    let mut model = Model::new();
    let cubes: Vec<Shape> = (0..50)
        .map(|i| cube(&mut model, f64::from(i), 0.0, 0.0, 1.0))
        .collect();
    let fused = fuse_all(&mut model, &cubes, T).unwrap();
    assert_eq!(model.kind_of(&fused.shape).unwrap(), ShapeType::Solid);
    assert_valid(&model, &fused.shape);
    let measured = volume(&model, &fused.shape);
    assert!((measured - 50.0).abs() < 1e-9, "volume {measured}");
    assert_eq!(count(&model, &fused.shape, ShapeType::Face), 6);
    assert_eq!(count(&model, &fused.shape, ShapeType::Edge), 12);
    assert_eq!(count(&model, &fused.shape, ShapeType::Vertex), 8);

    // Every face of every cube traces into the result, but the walls two
    // cubes shared, which are gone.
    let result_faces = explore_unique(&model, &fused.shape, ShapeType::Face).unwrap();
    for (i, solid) in cubes.iter().enumerate() {
        assert_eq!(
            fused.history.modified(solid),
            std::slice::from_ref(&fused.shape)
        );
        for face in explore_unique(&model, solid, ShapeType::Face).unwrap() {
            let (p, n) = ogeom::algo::face_normal(&model, &face, T).unwrap();
            let x_wall = n.x.abs() > 0.5;
            let inner = (x_wall && n.x < 0.0 && i > 0) || (x_wall && n.x > 0.0 && i < 49);
            let images = fused.history.trace(&face);
            if inner {
                assert!(fused.history.is_deleted(&face), "wall at x = {}", p.x);
            } else {
                assert_eq!(images.len(), 1, "face of cube {i} at {p:?}");
                assert!(result_faces.iter().any(|f| f.node() == images[0].node()));
            }
        }
        for edge in explore_unique(&model, solid, ShapeType::Edge).unwrap() {
            let images = fused.history.trace(&edge);
            assert!(
                images.is_empty()
                    || images.iter().all(|e| {
                        explore_unique(&model, &fused.shape, ShapeType::Edge)
                            .unwrap()
                            .iter()
                            .any(|r| r.node() == e.node())
                    }),
                "an edge of cube {i} traces outside the result: {:?} -> {:?} deleted {}",
                ogeom::algo::shape_bounds(&model, &edge, T).unwrap(),
                images.iter().map(|e| e.node()).collect::<Vec<_>>(),
                fused.history.is_deleted(&edge),
            );
        }
    }
}

#[test]
fn ten_cubes_apart_come_back_as_a_compound_of_the_ten() {
    let mut model = Model::new();
    let cubes: Vec<Shape> = (0..10)
        .map(|i| cube(&mut model, 3.0 * f64::from(i), 0.0, 0.0, 1.0))
        .collect();
    let fused = fuse_all(&mut model, &cubes, T).unwrap();
    assert_eq!(model.kind_of(&fused.shape).unwrap(), ShapeType::Compound);
    let solids = model.children_of(&fused.shape).unwrap();
    assert_eq!(solids.len(), 10);
    for solid in &cubes {
        assert!(solids.iter().any(|s| s == solid), "handed back as it came");
        assert!(!fused.history.is_affected(solid));
    }
    assert!((volume(&model, &fused.shape) - 10.0).abs() < 1e-9);
}

/// A 4 x 4 grid of boxes 1.5 on a side, 1 apart: each overlaps its eight
/// neighbours, the union is the 4.5 x 4.5 slab.
#[test]
fn an_overlapping_grid_gives_the_slab_whatever_the_order() {
    let mut model = Model::new();
    let boxes: Vec<Shape> = (0..16)
        .map(|k| {
            let (i, j) = (f64::from(k % 4), f64::from(k / 4));
            make_box(&mut model, at(i, j, 0.0), (1.5, 1.5, 1.0), T)
                .unwrap()
                .shape
        })
        .collect();
    let expected = 4.5 * 4.5;
    let mut faces = None;
    for seed in 0..4 {
        let order: Vec<Shape> = shuffled(16, seed)
            .into_iter()
            .map(|i| boxes[i].clone())
            .collect();
        let fused = fuse_all(&mut model, &order, T).unwrap();
        assert_eq!(model.kind_of(&fused.shape).unwrap(), ShapeType::Solid);
        assert_valid(&model, &fused.shape);
        let measured = volume(&model, &fused.shape);
        assert!(
            (measured - expected).abs() < 1e-9,
            "seed {seed}: volume {measured}"
        );
        let n = count(&model, &fused.shape, ShapeType::Face);
        assert_eq!(*faces.get_or_insert(n), n, "seed {seed}");
    }
    assert_eq!(faces, Some(6));

    let pairwise = chained(&mut model, &boxes);
    let measured = volume(&model, &pairwise);
    assert!((measured - expected).abs() < 1e-9);
    let (refined, _) = unify_same_domain(&mut model, &pairwise, T).unwrap();
    assert_eq!(count(&model, &refined.shape, ShapeType::Face), 6);
}

/// Area of the lens two circles of radius `r` with centres `d` apart
/// share.
fn lens(r: f64, d: f64) -> f64 {
    2.0 * r * r * (d / (2.0 * r)).acos() - d / 2.0 * (4.0 * r * r - d * d).sqrt()
}

/// Eight drums of radius 1 and height 2, centres 1.5 apart along x: each
/// overlaps the next in a lens, and none reaches the one after.
#[test]
fn a_row_of_overlapping_drums_has_the_inclusion_exclusion_volume() {
    let mut model = Model::new();
    let drums: Vec<Shape> = (0..8)
        .map(|i| {
            make_cylinder(&mut model, at(1.5 * f64::from(i), 0.0, 0.0), 1.0, 2.0, T)
                .unwrap()
                .shape
        })
        .collect();
    let expected = 2.0 * (8.0 * PI - 7.0 * lens(1.0, 1.5));
    let mut faces = None;
    for seed in 0..3 {
        let order: Vec<Shape> = shuffled(8, seed)
            .into_iter()
            .map(|i| drums[i].clone())
            .collect();
        let fused = fuse_all(&mut model, &order, T).unwrap();
        assert_eq!(model.kind_of(&fused.shape).unwrap(), ShapeType::Solid);
        assert_valid(&model, &fused.shape);
        let measured = volume(&model, &fused.shape);
        assert!(
            (measured - expected).abs() < expected * 1e-6,
            "seed {seed}: volume {measured} against {expected}"
        );
        let n = count(&model, &fused.shape, ShapeType::Face);
        assert_eq!(*faces.get_or_insert(n), n, "seed {seed}");
    }
    let pairwise = chained(&mut model, &drums);
    let (refined, _) = unify_same_domain(&mut model, &pairwise, T).unwrap();
    assert_eq!(
        faces,
        Some(count(&model, &refined.shape, ShapeType::Face)),
        "as the two-at-a-time union refined"
    );
}

/// Five drums stacked end to end on one axis are one drum: a wall, a top
/// and a bottom.
#[test]
fn drums_stacked_on_one_axis_merge_into_one_drum() {
    let mut model = Model::new();
    let drums: Vec<Shape> = (0..5)
        .map(|i| {
            make_cylinder(&mut model, at(0.0, 0.0, f64::from(i)), 1.0, 1.0, T)
                .unwrap()
                .shape
        })
        .collect();
    let fused = fuse_all(&mut model, &drums, T).unwrap();
    assert_valid(&model, &fused.shape);
    let measured = volume(&model, &fused.shape);
    assert!((measured - 5.0 * PI).abs() < 5.0 * PI * 1e-6, "{measured}");
    assert_eq!(count(&model, &fused.shape, ShapeType::Face), 3);
}

fn profile(model: &mut Model, points: &[(f64, f64)]) -> Shape {
    let pts: Vec<Point> = points.iter().map(|&(x, y)| Point::new(x, y, 0.0)).collect();
    let wire = make_polygon(model, &pts, true, T).unwrap().shape;
    let plane = PlaneSurface::new(Plane::new(Frame::WORLD)).into();
    make_face(model, plane, &[wire], T).unwrap().shape
}

/// Two prisms and two revolutions standing apart: none is dropped, by
/// `fuse` two at a time or by the union of all four.
#[test]
fn disjoint_swept_solids_are_all_kept() {
    let mut model = Model::new();
    let mut solids = Vec::new();
    for x in [0.0, 5.0] {
        let face = profile(
            &mut model,
            &[(x, 0.0), (x + 2.0, 0.0), (x + 2.0, 1.0), (x, 1.0)],
        );
        solids.push(
            make_prism(&mut model, &face, Vector::new(0.0, 0.0, 3.0), T)
                .unwrap()
                .shape,
        );
    }
    for y in [10.0, 20.0] {
        let face = profile(
            &mut model,
            &[(1.0, y), (2.0, y), (2.0, y + 2.0), (1.0, y + 2.0)],
        );
        let axis = Axis {
            location: Point::new(0.0, y, 0.0),
            direction: Direction::Y,
        };
        solids.push(
            make_revolution(&mut model, &face, axis, TAU, T)
                .unwrap()
                .shape,
        );
    }
    // Prisms 2 x 1 x 3; rings between radii 1 and 2, 2 high.
    let expected = 2.0 * 6.0 + 2.0 * (PI * (4.0 - 1.0) * 2.0);
    let pair = fuse(&mut model, &solids[0], &solids[1], T).unwrap();
    assert_eq!(model.children_of(&pair.shape).unwrap().len(), 2);
    assert!((volume(&model, &pair.shape) - 12.0).abs() < 1e-9);
    let pair = fuse(&mut model, &solids[2], &solids[3], T).unwrap();
    assert_eq!(model.children_of(&pair.shape).unwrap().len(), 2);

    let fused = fuse_all(&mut model, &solids, T).unwrap();
    assert_eq!(model.children_of(&fused.shape).unwrap().len(), 4);
    let measured = volume(&model, &fused.shape);
    assert!((measured - expected).abs() < expected * 1e-6, "{measured}");
}

/// Inputs given as compounds: the even cubes of a row of six in one, the
/// odd in another, and a seventh overlapping the first two of the first.
#[test]
fn compounds_of_solids_fuse_with_each_other_and_within() {
    let mut model = Model::new();
    let cubes: Vec<Shape> = (0..6)
        .map(|i| cube(&mut model, f64::from(i), 0.0, 0.0, 1.0))
        .collect();
    let even = model
        .add_compound(&[cubes[0].clone(), cubes[2].clone(), cubes[4].clone()])
        .unwrap();
    let odd = model
        .add_compound(&[cubes[1].clone(), cubes[3].clone(), cubes[5].clone()])
        .unwrap();
    let fused = fuse_all(&mut model, &[even.clone(), odd.clone()], T).unwrap();
    assert_eq!(model.kind_of(&fused.shape).unwrap(), ShapeType::Solid);
    assert!((volume(&model, &fused.shape) - 6.0).abs() < 1e-9);
    assert_eq!(count(&model, &fused.shape, ShapeType::Face), 6);
    assert_eq!(
        fused.history.modified(&even),
        std::slice::from_ref(&fused.shape)
    );

    // One compound whose solids overlap: they are fused too.
    let a = cube(&mut model, 0.0, 0.0, 0.0, 2.0);
    let b = cube(&mut model, 1.0, 1.0, 1.0, 2.0);
    let both = model.add_compound(&[a, b]).unwrap();
    let fused = fuse_all(&mut model, &[both], T).unwrap();
    assert_eq!(model.kind_of(&fused.shape).unwrap(), ShapeType::Solid);
    assert!((volume(&model, &fused.shape) - 15.0).abs() < 1e-9);
}

/// A cube inside another is taken in by it, and traces to the solid that
/// holds it.
#[test]
fn a_solid_inside_another_is_taken_in() {
    let mut model = Model::new();
    let outer = cube(&mut model, 0.0, 0.0, 0.0, 4.0);
    let inner = cube(&mut model, 1.0, 1.0, 1.0, 1.0);
    let beside = cube(&mut model, 4.0, 0.0, 0.0, 1.0);
    let fused = fuse_all(&mut model, &[inner.clone(), outer, beside], T).unwrap();
    assert_valid(&model, &fused.shape);
    assert!((volume(&model, &fused.shape) - 65.0).abs() < 1e-9);
    assert_eq!(
        fused.history.modified(&inner),
        std::slice::from_ref(&fused.shape)
    );
    for face in explore_unique(&model, &inner, ShapeType::Face).unwrap() {
        assert!(fused.history.is_deleted(&face));
    }
}

#[test]
fn an_input_that_is_not_a_solid_is_named() {
    let mut model = Model::new();
    let a = cube(&mut model, 0.0, 0.0, 0.0, 1.0);
    let b = cube(&mut model, 1.0, 0.0, 0.0, 1.0);
    let sheet = profile(&mut model, &[(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)]);
    let err = fuse_all(&mut model, &[a, b, sheet], T).unwrap_err();
    assert!(matches!(err, OgeomError::Construction(_)), "{err}");
    assert!(err.to_string().contains("input 2"), "{err}");
    let err = fuse_all(&mut model, &[], T).unwrap_err();
    assert!(matches!(err, OgeomError::Construction(_)), "{err}");
}

/// Two cubes sharing only an edge, as `fuse` makes of the pair.
#[test]
fn cubes_meeting_along_an_edge_join_as_the_pair_fuses() {
    let mut model = Model::new();
    let a = cube(&mut model, 0.0, 0.0, 0.0, 1.0);
    let b = cube(&mut model, 1.0, 1.0, 0.0, 1.0);
    let pair = fuse(&mut model, &a, &b, T).unwrap();
    let fused = fuse_all(&mut model, &[a, b], T).unwrap();
    assert_eq!(
        model.kind_of(&fused.shape).unwrap(),
        model.kind_of(&pair.shape).unwrap()
    );
    assert_eq!(
        count(&model, &fused.shape, ShapeType::Solid),
        count(&model, &pair.shape, ShapeType::Solid)
    );
    assert!((volume(&model, &fused.shape) - 2.0).abs() < 1e-9);
}

/// A boolean refusal names the two inputs that do not fuse: a half space
/// among cubes, which no union takes.
#[test]
fn a_refusal_names_the_inputs_that_do_not_fuse() {
    let mut model = Model::new();
    let mut solids: Vec<Shape> = (0..3)
        .map(|i| cube(&mut model, f64::from(i), 0.0, -0.5, 1.0))
        .collect();
    let face = profile(
        &mut model,
        &[(0.5, 0.2), (2.5, 0.2), (2.5, 0.8), (0.5, 0.8)],
    );
    let upper = ogeom::algo::make_half_space(&mut model, &face, Point::new(0.0, 0.0, 100.0), T)
        .unwrap()
        .shape;
    solids.insert(1, upper);
    let err = fuse_all(&mut model, &solids, T).unwrap_err();
    assert!(matches!(err, OgeomError::Construction(_)), "{err}");
    let message = err.to_string();
    assert!(
        message.contains("input 1") && message.contains("do not fuse"),
        "{message}"
    );
}
