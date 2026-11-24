//! Walls on both sides of a solid's faces or a sheet: closed-form volumes,
//! the rim left at each opening, the history to both walls, every wall's
//! distance from its source face sampled, and the refusals by name.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_algo::Built;
use ogeom_core::Tolerances;
use ogeom_geom::SurfaceGeometry;
use ogeom_math::{Frame, Point, Vector};
use ogeom_offset::{Join, Walls, make_thick_sheet_sided, make_thick_solid_sided};
use ogeom_topo::{Filter, Model, NodeData, Shape, ShapeType, explore};
use std::f64::consts::PI;

const T: Tolerances = Tolerances::millimetres();

fn fine() -> ogeom_mesh::Deflection {
    ogeom_mesh::Deflection {
        chord: 1e-4,
        ..ogeom_mesh::Deflection::default()
    }
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    ogeom_algo::volume_properties(model, shape, fine(), T)
        .unwrap()
        .mass
}

fn area(model: &Model, shape: &Shape) -> f64 {
    ogeom_algo::surface_properties(model, shape, fine(), T)
        .unwrap()
        .mass
}

fn faces(model: &Model, shape: &Shape) -> Vec<Shape> {
    explore(model, shape, Filter::OfType(ShapeType::Face)).unwrap()
}

fn surface<'m>(model: &'m Model, face: &Shape) -> &'m SurfaceGeometry {
    let Some(NodeData::Face(data)) = model.node(face).map(ogeom_topo::TShape::data) else {
        panic!("not a face");
    };
    model.geometry().surface(data.surface).unwrap()
}

/// One sound solid of one or more shells, and its volume.
fn sound(model: &Model, built: &Built, shells: usize) -> f64 {
    let diagnosis = ogeom_algo::check(model, &built.shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let solids = explore(model, &built.shape, Filter::OfType(ShapeType::Solid)).unwrap();
    assert_eq!(solids.len(), 1);
    let found = explore(model, &built.shape, Filter::OfType(ShapeType::Shell)).unwrap();
    assert_eq!(found.len(), shells);
    volume(model, &built.shape)
}

/// The face of `solid` whose bounds lie in the plane `z`.
fn face_in_z(model: &Model, solid: &Shape, z: f64) -> Shape {
    faces(model, solid)
        .into_iter()
        .find(|f| {
            let b = ogeom_algo::shape_bounds(model, f, T).unwrap();
            (b.low().unwrap().z - z).abs() < 1e-4 && (b.high().unwrap().z - z).abs() < 1e-4
        })
        .expect("a face in that plane")
}

fn cube(model: &mut Model) -> Shape {
    ogeom_algo::make_box(model, Frame::WORLD, (20.0, 20.0, 20.0), T)
        .unwrap()
        .shape
}

fn assert_near(got: f64, want: f64, rel: f64, what: &str) {
    eprintln!("{what}: measured {got}, closed form {want}, bound {rel} relative");
    assert!(
        (got - want).abs() <= want.abs() * rel,
        "{what}: {got} against {want}"
    );
}

/// The distance of `p` from the whole (untrimmed) surface, positive on
/// the side its own normal points to.
fn signed_distance(surface: &SurfaceGeometry, p: Point) -> f64 {
    match surface {
        SurfaceGeometry::Plane(s) => s.plane().signed_distance_to(p),
        SurfaceGeometry::Cylinder(s) => s.cylinder().signed_distance_to(p),
        other => panic!("no closed-form distance to {other:?}"),
    }
}

/// The largest distance, over samples of every wall the history gives a
/// kept source face, between the sample's distance from the source face's
/// surface and the depth of the wall on the sample's side.
fn worst_wall_error(
    model: &Model,
    source: &Shape,
    removed: &[Shape],
    built: &Built,
    walls: Walls,
) -> f64 {
    let mut worst: f64 = 0.0;
    let mut sampled = 0;
    for face in faces(model, source) {
        if removed.iter().any(|r| r.node() == face.node()) {
            continue;
        }
        let images = built.history.modified(&face);
        assert_eq!(images.len(), 2, "a kept face becomes its two walls");
        let own = surface(model, &face).clone();
        let normal_sign = if face.orientation() == ogeom_topo::Orientation::Reversed {
            -1.0
        } else {
            1.0
        };
        let mut sides = Vec::new();
        for image in images {
            let mesh = ogeom_mesh::triangulate_face(model, image, fine(), T).unwrap();
            assert!(!mesh.positions.is_empty());
            let mut side = 0.0;
            for p in &mesh.positions {
                let along = normal_sign * signed_distance(&own, *p);
                let ahead = along > 0.0;
                side = if ahead { 1.0 } else { -1.0 };
                let want = if ahead { walls.ahead } else { walls.behind };
                worst = worst.max((along.abs() - want).abs());
                sampled += 1;
            }
            sides.push(side);
        }
        assert!(
            sides.contains(&1.0) && sides.contains(&-1.0),
            "one wall ahead of each kept face and one behind it"
        );
    }
    assert!(sampled > 0);
    worst
}

/// The issue's case: a 20 mm cube open at the top, a 2 mm wall centred on
/// its faces. Outer box 22 wide from z = -1, inner 18 wide from z = 1, both
/// stopping on the top's plane.
#[test]
fn a_cube_walled_both_sides_is_the_outer_box_less_the_inner() {
    let mut model = Model::new();
    let block = cube(&mut model);
    let top = face_in_z(&model, &block, 20.0);
    let walls = Walls::both(2.0);
    let built = make_thick_solid_sided(
        &mut model,
        &block,
        std::slice::from_ref(&top),
        walls,
        Join::Intersection,
        T,
    )
    .unwrap();
    let got = sound(&model, &built, 1);
    assert_near(got, 22.0 * 22.0 * 21.0 - 18.0 * 18.0 * 19.0, 1e-9, "volume");

    // Five outer walls, five inner ones and one rim: the rim is a single
    // flat annulus from the outer wall to the inner, not split along the
    // old top's edges.
    assert_eq!(faces(&model, &built.shape).len(), 11);
    let rim = face_in_z(&model, &built.shape, 20.0);
    assert_near(area(&model, &rim), 22.0 * 22.0 - 18.0 * 18.0, 1e-9, "rim");
    assert!(built.history.is_deleted(&top));

    let worst = worst_wall_error(&model, &block, &[top], &built, walls);
    eprintln!("wall depth error: measured {worst}, bound 1e-9");
    assert!(worst < 1e-9, "wall off its depth by {worst}");
}

/// A split wall: 1 behind the faces, 3 ahead.
#[test]
fn a_split_wall_stands_each_depth_on_its_own_side() {
    let mut model = Model::new();
    let block = cube(&mut model);
    let top = face_in_z(&model, &block, 20.0);
    let walls = Walls::split(1.0, 3.0);
    let built = make_thick_solid_sided(
        &mut model,
        &block,
        std::slice::from_ref(&top),
        walls,
        Join::Intersection,
        T,
    )
    .unwrap();
    let got = sound(&model, &built, 1);
    assert_near(got, 26.0 * 26.0 * 23.0 - 18.0 * 18.0 * 19.0, 1e-9, "volume");
    let worst = worst_wall_error(&model, &block, &[top], &built, walls);
    eprintln!("wall depth error: measured {worst}, bound 1e-9");
    assert!(worst < 1e-9, "wall off its depth by {worst}");
}

/// Nothing removed: a closed wall round a sealed cavity, two shells.
#[test]
fn a_cube_with_nothing_removed_is_a_sealed_wall() {
    let mut model = Model::new();
    let block = cube(&mut model);
    let built = make_thick_solid_sided(
        &mut model,
        &block,
        &[],
        Walls::both(2.0),
        Join::Intersection,
        T,
    )
    .unwrap();
    let got = sound(&model, &built, 2);
    assert_near(got, 22.0_f64.powi(3) - 18.0_f64.powi(3), 1e-9, "volume");
}

/// A drum of radius 10 and height 20 open at the top, 1 behind and 2
/// ahead: the cylinder of radius 12 from z = -2 less the one of radius 9
/// from z = 1.
#[test]
fn a_drum_walled_both_sides_is_two_coaxial_cylinders() {
    let mut model = Model::new();
    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 10.0, 20.0, T)
        .unwrap()
        .shape;
    let top = face_in_z(&model, &drum, 20.0);
    let walls = Walls::split(1.0, 2.0);
    let built = make_thick_solid_sided(
        &mut model,
        &drum,
        std::slice::from_ref(&top),
        walls,
        Join::Intersection,
        T,
    )
    .unwrap();
    let got = sound(&model, &built, 1);
    let want = PI * (12.0 * 12.0 * 22.0 - 9.0 * 9.0 * 19.0);
    assert_near(got, want, 1e-6, "volume");
    let rim = face_in_z(&model, &built.shape, 20.0);
    assert_near(area(&model, &rim), PI * (144.0 - 81.0), 1e-6, "rim");
    let worst = worst_wall_error(&model, &drum, &[top], &built, walls);
    eprintln!("wall depth error: measured {worst}, bound 1e-7");
    assert!(worst < 1e-7, "wall off its depth by {worst}");
}

/// An L-shaped prism, open at the top: the concave edge mitres on both
/// sides. Outer outline 22 x 22 less 10 x 10, inner 18 x 18 less 10 x 10.
#[test]
fn an_l_prism_mitres_its_concave_corner_on_both_sides() {
    let mut model = Model::new();
    let outline = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(20.0, 0.0, 0.0),
        Point::new(20.0, 10.0, 0.0),
        Point::new(10.0, 10.0, 0.0),
        Point::new(10.0, 20.0, 0.0),
        Point::new(0.0, 20.0, 0.0),
    ];
    let wire = ogeom_algo::make_polygon(&mut model, &outline, true, T)
        .unwrap()
        .shape;
    let plane = ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(Frame::WORLD));
    let profile = ogeom_algo::make_face(&mut model, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let prism = ogeom_algo::make_prism(&mut model, &profile, Vector::new(0.0, 0.0, 10.0), T)
        .unwrap()
        .shape;
    let top = face_in_z(&model, &prism, 10.0);
    let walls = Walls::both(2.0);
    let built = make_thick_solid_sided(
        &mut model,
        &prism,
        std::slice::from_ref(&top),
        walls,
        Join::Intersection,
        T,
    )
    .unwrap();
    let got = sound(&model, &built, 1);
    assert_near(got, 384.0 * 11.0 - 224.0 * 9.0, 1e-9, "volume");
    let worst = worst_wall_error(&model, &prism, &[top], &built, walls);
    eprintln!("wall depth error: measured {worst}, bound 1e-9");
    assert!(worst < 1e-9, "wall off its depth by {worst}");
}

/// The arc join rounds the outer wall's corners only on a convex solid:
/// the outer copy is the cube's rolling-ball parallel body cut at the top
/// (a quarter cylinder about each vertical and bottom edge, a ball's eighth
/// at each bottom corner), the inner copy the sharp 18 mm box.
#[test]
fn the_arc_join_rounds_the_outer_corners_only() {
    let mut model = Model::new();
    let block = cube(&mut model);
    let top = face_in_z(&model, &block, 20.0);
    let built = make_thick_solid_sided(
        &mut model,
        &block,
        std::slice::from_ref(&top),
        Walls::both(2.0),
        Join::Arc,
        T,
    )
    .unwrap();
    let got = sound(&model, &built, 1);
    let outer = 8000.0 + 1600.0 + 400.0 + 40.0 * PI + 2.0 / 3.0 * PI;
    let want = outer - 18.0 * 18.0 * 19.0;
    assert_near(got, want, 1e-4, "volume");

    // Every curved face stands outside the cube; every face inside it is a
    // plane.
    let probe = ogeom_algo::PointDistance::of(&model, &block, T).unwrap();
    let mut curved = 0;
    for face in faces(&model, &built.shape) {
        if matches!(surface(&model, &face), SurfaceGeometry::Plane(_)) {
            continue;
        }
        curved += 1;
        let b = ogeom_algo::shape_bounds(&model, &face, T).unwrap();
        let (lo, hi) = (b.low().unwrap(), b.high().unwrap());
        let slack = 1e-4;
        let outside = hi.x <= slack
            || lo.x >= 20.0 - slack
            || hi.y <= slack
            || lo.y >= 20.0 - slack
            || hi.z <= slack;
        assert!(outside, "a rounded face inside the cube: {b:?}");
        // The rounded outer wall stands its depth off the cube everywhere.
        let mesh = ogeom_mesh::triangulate_face(&model, &face, fine(), T).unwrap();
        for p in &mesh.positions {
            let (d, _) = probe.nearest(*p, T).unwrap();
            assert!((d - 1.0).abs() < 1e-6, "{p:?} stands {d} off the cube");
        }
    }
    assert_eq!(
        curved, 12,
        "four vertical and four bottom bands, four balls"
    );
}

/// A square sheet thickened 1 to its back and 3 to its front; a drum's
/// wall, 1 inward and 2 outward.
#[test]
fn sheets_thicken_to_each_side_by_its_depth() {
    let mut model = Model::new();
    let corners = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(10.0, 0.0, 0.0),
        Point::new(10.0, 10.0, 0.0),
        Point::new(0.0, 10.0, 0.0),
    ];
    let wire = ogeom_algo::make_polygon(&mut model, &corners, true, T)
        .unwrap()
        .shape;
    let plane = ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(Frame::WORLD));
    let square = ogeom_algo::make_face(&mut model, plane.into(), &[wire], T)
        .unwrap()
        .shape;
    let slab = make_thick_sheet_sided(&mut model, &square, Walls::split(1.0, 3.0), T).unwrap();
    let got = sound(&model, &slab, 1);
    assert_near(got, 400.0, 1e-9, "slab volume");
    let b = ogeom_algo::shape_bounds(&model, &slab.shape, T).unwrap();
    assert!((b.low().unwrap().z + 1.0).abs() < 1e-6, "{b:?}");
    assert!((b.high().unwrap().z - 3.0).abs() < 1e-6, "{b:?}");

    let drum = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 4.0, 6.0, T)
        .unwrap()
        .shape;
    let wall = faces(&model, &drum)
        .into_iter()
        .find(|f| matches!(surface(&model, f), SurfaceGeometry::Cylinder(_)))
        .unwrap();
    let tube = make_thick_sheet_sided(&mut model, &wall, Walls::split(1.0, 2.0), T).unwrap();
    let got = sound(&model, &tube, 1);
    assert_near(got, PI * (36.0 - 9.0) * 6.0, 1e-6, "tube volume");
}

/// What cannot be built is refused by name: a side deeper than the solid
/// holds, a negative or missing depth, an opening with a tangent
/// neighbour, a face of another solid.
#[test]
fn walls_that_cannot_be_built_are_refused() {
    let refused = |result: ogeom_core::OgeomResult<Built>, says: &str| {
        let err = result.expect_err("refused").to_string();
        assert!(err.contains(says), "{err}");
    };
    let mut model = Model::new();
    let block = cube(&mut model);
    let top = face_in_z(&model, &block, 20.0);
    let one = std::slice::from_ref(&top);
    refused(
        make_thick_solid_sided(
            &mut model,
            &block,
            one,
            Walls::split(10.5, 1.0),
            Join::Intersection,
            T,
        ),
        "more than the solid has room for",
    );
    refused(
        make_thick_solid_sided(
            &mut model,
            &block,
            one,
            Walls::split(-1.0, 1.0),
            Join::Intersection,
            T,
        ),
        "is not a depth",
    );
    refused(
        make_thick_solid_sided(
            &mut model,
            &block,
            one,
            Walls::split(0.0, 0.0),
            Join::Intersection,
            T,
        ),
        "hold nothing",
    );
    refused(
        make_thick_solid_sided(
            &mut model,
            &block,
            one,
            Walls::both(f64::NAN),
            Join::Intersection,
            T,
        ),
        "is not a depth",
    );
    let other = cube(&mut model);
    let foreign = face_in_z(&model, &other, 20.0);
    refused(
        make_thick_solid_sided(
            &mut model,
            &block,
            &[foreign],
            Walls::both(2.0),
            Join::Intersection,
            T,
        ),
        "not a face of the solid",
    );

    // The cube with its top edges rounded: the top meets the rounds
    // tangentially.
    let top_edges: Vec<Shape> = explore(&model, &top, Filter::OfType(ShapeType::Edge)).unwrap();
    let rounded = ogeom_fillet::fillet_edges(&mut model, &block, &top_edges, 2.0, T)
        .unwrap()
        .shape;
    let flat = face_in_z(&model, &rounded, 20.0);
    refused(
        make_thick_solid_sided(
            &mut model,
            &rounded,
            &[flat],
            Walls::both(2.0),
            Join::Intersection,
            T,
        ),
        "tangentially",
    );
}

/// One side at zero is the one-sided hollowing itself.
#[test]
fn a_wall_on_one_side_is_the_one_sided_hollowing() {
    for (walls, want) in [
        (Walls::behind(2.0), 8000.0 - 16.0 * 16.0 * 18.0),
        (Walls::ahead(2.0), 24.0 * 24.0 * 22.0 - 8000.0),
    ] {
        let mut model = Model::new();
        let block = cube(&mut model);
        let top = face_in_z(&model, &block, 20.0);
        let built = make_thick_solid_sided(
            &mut model,
            &block,
            std::slice::from_ref(&top),
            walls,
            Join::Intersection,
            T,
        )
        .unwrap();
        let got = sound(&model, &built, 1);
        assert_near(got, want, 1e-9, "volume");
    }
}
