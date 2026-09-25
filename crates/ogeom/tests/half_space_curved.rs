//! Half spaces bounded by curved faces: the side of a whole surface, as an
//! extrusion stopping on a curved face needs it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_box, make_half_space, make_natural_face, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::{CylinderSurface, SphereSurface, SurfaceGeometry};
use ogeom::math::{Cylinder, Frame, Point, Sphere};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-4).unwrap(), T)
        .unwrap()
        .mass
}

/// The 40 x 4 x 10 block from x = -20 across a rod or ball of radius 10.
fn block(model: &mut Model) -> Shape {
    let at = Frame::new(
        Point::new(-20.0, -2.0, 5.0),
        ogeom::math::Direction::Z,
        ogeom::math::Direction::X,
        T,
    )
    .unwrap();
    make_box(model, at, (40.0, 4.0, 10.0), T).unwrap().shape
}

/// Common and cut of the block with the half space on either side of
/// `surface`, against the volume `outside` the surface holds of the block.
fn both_sides(surface: SurfaceGeometry, outside: f64, within: Point) {
    let mut model = Model::new();
    let block = block(&mut model);
    let whole = 1600.0;
    let face = make_natural_face(&mut model, surface).unwrap().shape;
    for (inside, material) in [
        (Point::new(30.0, 0.0, 10.0), outside),
        (within, whole - outside),
    ] {
        let side = make_half_space(&mut model, &face, inside, T).unwrap().shape;
        let kept = ogeom::boolean::common(&mut model, &block, &side, T)
            .unwrap()
            .shape;
        assert!(check(&model, &kept, T).unwrap().is_valid());
        let got = volume(&model, &kept);
        assert!(
            (got - material).abs() < material * 1e-6,
            "common from {inside:?}: {got} against {material}"
        );
        let left = ogeom::boolean::cut(&mut model, &block, &side, T)
            .unwrap()
            .shape;
        assert!(check(&model, &left, T).unwrap().is_valid());
        let got = volume(&model, &left);
        let rest = whole - material;
        assert!(
            (got - rest).abs() < rest * 1e-6,
            "cut from {inside:?}: {got} against {rest}"
        );
    }
}

/// A rod of radius 10 along z: the block outside it is two pieces of
/// `10 (160 - 2 (2 sqrt 96 + 100 asin 0.2))`, inside it the rest.
#[test]
fn a_half_space_bounded_by_a_rod_keeps_its_side_of_a_block() {
    let rod = Cylinder::new(Frame::WORLD, 10.0, T).unwrap();
    let surface: SurfaceGeometry = CylinderSurface::new(rod, (-50.0, 50.0)).unwrap().into();
    let outside = 10.0 * (160.0 - 2.0 * (2.0 * 96.0_f64.sqrt() + 100.0 * 0.2_f64.asin()));
    both_sides(surface, outside, Point::new(0.0, 0.0, 10.0));
}

/// A ball of radius 10 at the origin: inside it the block keeps the ball's
/// slab `|y| < 2, 5 < z < 15`, integrated here across `y` and `z`.
#[test]
fn a_half_space_bounded_by_a_ball_keeps_its_side_of_a_block() {
    let ball = Sphere::new(Frame::WORLD, 10.0, T).unwrap();
    let surface: SurfaceGeometry = SphereSurface::new(ball).into();
    // The ball's chord along x over the block's (y, z) rectangle, by the
    // midpoint rule in each direction: smooth where the chord is real, and
    // the rectangle stays inside the ball's shadow except near z = 10,
    // where the chord closes; enough points take that below the tolerance.
    let n = 2000;
    let mut inside = 0.0;
    for i in 0..n {
        let y = -2.0 + 4.0 * (f64::from(i) + 0.5) / f64::from(n);
        for j in 0..n {
            let z = 5.0 + 10.0 * (f64::from(j) + 0.5) / f64::from(n);
            let s = 100.0 - y * y - z * z;
            if s > 0.0 {
                inside += 2.0 * s.sqrt() * (4.0 / f64::from(n)) * (10.0 / f64::from(n));
            }
        }
    }
    both_sides(surface, 1600.0 - inside, Point::ORIGIN);
}

/// The block's shares on either side of `surface`, from `outer` and from
/// `inner`: each operation valid, the two sides' commons adding up to the
/// block, and each cut the other side's common.
fn shares_add_up(surface: SurfaceGeometry, outer: Point, inner: Point) {
    let mut model = Model::new();
    let block = block(&mut model);
    let face = make_natural_face(&mut model, surface).unwrap().shape;
    let mut commons = Vec::new();
    for inside in [outer, inner] {
        let side = make_half_space(&mut model, &face, inside, T).unwrap().shape;
        let kept = ogeom::boolean::common(&mut model, &block, &side, T)
            .unwrap()
            .shape;
        assert!(check(&model, &kept, T).unwrap().is_valid());
        let left = ogeom::boolean::cut(&mut model, &block, &side, T)
            .unwrap()
            .shape;
        assert!(check(&model, &left, T).unwrap().is_valid());
        let (k, l) = (volume(&model, &kept), volume(&model, &left));
        assert!(
            (k + l - 1600.0).abs() < 1e-6 * 1600.0,
            "{k} + {l} from {inside:?}"
        );
        assert!(k > 1.0 && l > 1.0, "both sides of the block: {k}, {l}");
        commons.push(k);
    }
    assert!(
        (commons[0] + commons[1] - 1600.0).abs() < 1e-6 * 1600.0,
        "{commons:?}"
    );
}

/// A cone widening upward through the block, its apex below it, and a ring
/// lying in the block's middle: each side's share is the rest of the other's.
#[test]
fn half_spaces_bounded_by_a_cone_and_a_ring_divide_a_block() {
    use ogeom::geom::{ConeSurface, TorusSurface};
    use ogeom::math::{Cone, Torus};
    let apex_below = Frame::new(
        Point::new(0.0, 0.0, -10.0),
        ogeom::math::Direction::Z,
        ogeom::math::Direction::X,
        T,
    )
    .unwrap();
    let cone = Cone::new(apex_below, 0.0, 0.5, T).unwrap();
    let cone: SurfaceGeometry = ConeSurface::new(cone, (0.0, 60.0)).unwrap().into();
    shares_add_up(
        cone,
        Point::new(19.0, 0.0, 10.0),
        Point::new(0.0, 0.0, 10.0),
    );
    let middle = Frame::new(
        Point::new(0.0, 0.0, 10.0),
        ogeom::math::Direction::Z,
        ogeom::math::Direction::X,
        T,
    )
    .unwrap();
    let ring = Torus::new(middle, 12.0, 3.0, T).unwrap();
    shares_add_up(
        TorusSurface::new(ring).into(),
        Point::new(0.0, 0.0, 10.0),
        Point::new(12.0, 0.0, 10.0),
    );
}

/// A spline patch lying flat at z = 10 across x and y from -30 to 30: the
/// block below it is half the block. A block reaching past the patch's
/// edge lies partly where the patch divides nothing, and is refused.
#[test]
fn a_half_space_bounded_by_a_spline_patch_divides_what_lies_across_it() {
    use ogeom::geom::BSplineSurface;
    use ogeom::math::{ControlGrid, KnotVector};
    let patch = |half: f64| -> SurfaceGeometry {
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0], 3).unwrap();
        let mut points = Vec::new();
        for i in 0..4 {
            for j in 0..4 {
                let (x, y) = (
                    -half + 2.0 * half * f64::from(i) / 3.0,
                    -half + 2.0 * half * f64::from(j) / 3.0,
                );
                points.push(Point::new(x, y, 10.0));
            }
        }
        let grid = ControlGrid::new(points, 4, 4).unwrap();
        BSplineSurface::new(knots.clone(), knots, &grid, T)
            .unwrap()
            .into()
    };
    let mut model = Model::new();
    let block = block(&mut model);
    let face = make_natural_face(&mut model, patch(30.0)).unwrap().shape;
    let below = make_half_space(&mut model, &face, Point::new(0.0, 0.0, 0.0), T)
        .unwrap()
        .shape;
    let kept = ogeom::boolean::common(&mut model, &block, &below, T)
        .unwrap()
        .shape;
    assert!(check(&model, &kept, T).unwrap().is_valid());
    let got = volume(&model, &kept);
    assert!((got - 800.0).abs() < 800.0 * 1e-6, "{got}");
    let left = ogeom::boolean::cut(&mut model, &block, &below, T)
        .unwrap()
        .shape;
    let got = volume(&model, &left);
    assert!((got - 800.0).abs() < 800.0 * 1e-6, "{got}");

    let narrow = make_natural_face(&mut model, patch(15.0)).unwrap().shape;
    let below = make_half_space(&mut model, &narrow, Point::new(0.0, 0.0, 0.0), T)
        .unwrap()
        .shape;
    let refused = ogeom::boolean::common(&mut model, &block, &below, T).unwrap_err();
    assert!(refused.to_string().contains("reaches past"), "{refused}");
}
