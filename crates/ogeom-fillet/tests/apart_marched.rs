//! Rounds standing apart on a block with two bores crossing: a crease
//! where the bores meet, which only the march seats, asked with box edges
//! or a bore's rim far from it. Each group is rounded as the engine rounds
//! it alone, so what the request takes off is what each takes off alone:
//! a box edge's round takes the prism `r^2 (1 - pi/4)` along its length,
//! a whole rim's the ring Pappus gives, and the seam's piece what its own
//! round takes.
#![allow(
    deprecated,
    reason = "exercises the edge-blend calls BlendRequest replaces"
)]
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::Curve3d as _;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// A block of 20, bored along `x` (radius 4, axis at `y = z = 10`) and
/// along `y` (radius 3, axis at `x = 10`, `z = 12`).
fn cross_bores(m: &mut Model) -> Shape {
    let frame = |at: Point, z: Direction, x: Direction| Frame::new(at, z, x, T).unwrap();
    let block = ogeom_algo::make_box(
        m,
        frame(Point::ORIGIN, Direction::Z, Direction::X),
        (20.0, 20.0, 20.0),
        T,
    )
    .unwrap()
    .shape;
    let along_x = ogeom_algo::make_cylinder(
        m,
        frame(Point::new(-1.0, 10.0, 10.0), Direction::X, Direction::Y),
        4.0,
        22.0,
        T,
    )
    .unwrap()
    .shape;
    let along_y = ogeom_algo::make_cylinder(
        m,
        frame(Point::new(10.0, -1.0, 12.0), Direction::Y, Direction::Z),
        3.0,
        22.0,
        T,
    )
    .unwrap()
    .shape;
    let once = ogeom_bool::cut(m, &block, &along_x, T).unwrap().shape;
    ogeom_bool::cut(m, &once, &along_y, T).unwrap().shape
}

/// The edge of `solid` running between `a` and `b` (either way) through
/// `mid`, its curve's middle.
fn edge(m: &Model, solid: &Shape, a: Point, b: Point, mid: Point) -> Shape {
    let near = |p: Point, q: Point| p.distance(q) < 1e-5;
    explore_unique(m, solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let Some((curve, range)) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
                return false;
            };
            let at = |t: f64| curve.point_at(t, T).unwrap();
            let (s, f) = (at(range.0), at(range.1));
            ((near(s, a) && near(f, b)) || (near(s, b) && near(f, a)))
                && near(at(0.5 * (range.0 + range.1)), mid)
        })
        .expect("the edge")
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// `edges` of `solid` rounded at `radius`: the result checked valid, and
/// its volume.
fn rounded(m: &mut Model, solid: &Shape, edges: &[Shape], radius: f64) -> f64 {
    let built = ogeom_fillet::fillet_edges(m, solid, edges, radius, T).unwrap();
    assert!(
        ogeom_algo::check(m, &built.shape, T).unwrap().is_valid(),
        "the rounded block is valid"
    );
    for e in edges {
        assert!(
            !built.history.generated(e).is_empty(),
            "every edge asked names faces of its round"
        );
    }
    volume(m, &built.shape)
}

/// The cross-section a round of `radius` takes off a right-angled corner.
fn spandrel(radius: f64) -> f64 {
    radius * radius * (1.0 - core::f64::consts::FRAC_PI_4)
}

/// The piece of the bores' seam on the `x` bore's top, between the
/// points where it crosses the line `y = 14, z = 10`.
fn seam_under(m: &Model, solid: &Shape) -> Shape {
    edge(
        m,
        solid,
        Point::new(12.236_068, 14.0, 10.0),
        Point::new(7.763_932, 14.0, 10.0),
        Point::new(10.0, 13.872_983, 9.0),
    )
}

/// The piece of the bores' seam from the line `y = 14, z = 10` up to
/// where the `y` bore's wall turns.
fn seam_rising(m: &Model, solid: &Shape) -> Shape {
    edge(
        m,
        solid,
        Point::new(7.763_932, 14.0, 10.0),
        Point::new(7.162_783, 13.866_384, 11.025_218),
        Point::new(7.417_380, 13.971_871, 10.473_542),
    )
}

/// Two box edges sharing no corner and a piece of the bores' seam: the
/// box edges' prisms and the seam's own round, together.
#[test]
fn box_edges_and_a_bore_seam_round_as_each_alone() {
    let mut m = Model::new();
    let solid = cross_bores(&mut m);
    let r = 1.327_036_193_144_703;
    let upright = edge(
        &m,
        &solid,
        Point::new(0.0, 0.0, 0.0),
        Point::new(0.0, 0.0, 20.0),
        Point::new(0.0, 0.0, 10.0),
    );
    let floor = edge(
        &m,
        &solid,
        Point::new(20.0, 20.0, 0.0),
        Point::new(0.0, 20.0, 0.0),
        Point::new(10.0, 20.0, 0.0),
    );
    let seam = seam_under(&m, &solid);
    let before = volume(&m, &solid);
    let seam_alone = before - rounded(&mut m.clone(), &solid, core::slice::from_ref(&seam), r);
    let got = before - rounded(&mut m, &solid, &[upright, floor, seam], r);
    let want = 2.0 * spandrel(r) * 20.0 + seam_alone;
    assert!(
        (got - want).abs() <= 1e-7 * before,
        "took off {got}, the prisms and the seam's round alone {want}"
    );
}

/// A rim of the `y` bore, a box edge and a piece of the bores' seam. The
/// rim is an arc of a circle the boolean split; asked with one radius, the
/// circle rounds whole, the ring Pappus gives.
#[test]
fn a_bore_rim_a_box_edge_and_a_bore_seam_round_as_each_alone() {
    let mut m = Model::new();
    let solid = cross_bores(&mut m);
    let r = 1.586_401_862_563_595;
    let rim = edge(
        &m,
        &solid,
        Point::new(10.0, 20.0, 9.0),
        Point::new(10.0, 20.0, 15.0),
        Point::new(7.0, 20.0, 12.0),
    );
    let corner = edge(
        &m,
        &solid,
        Point::new(20.0, 0.0, 0.0),
        Point::new(20.0, 0.0, 20.0),
        Point::new(20.0, 0.0, 10.0),
    );
    let seam = seam_rising(&m, &solid);
    let before = volume(&m, &solid);
    let seam_alone = before - rounded(&mut m.clone(), &solid, core::slice::from_ref(&seam), r);
    let got = before - rounded(&mut m, &solid, &[rim, corner, seam], r);
    // The spandrel's centroid stands `r (10 - 3 pi) / (12 - 3 pi)` off
    // each of its legs, so off the bore's axis by the bore's radius and
    // that.
    let pi = core::f64::consts::PI;
    let centroid = 3.0 + r * (10.0 - 3.0 * pi) / (12.0 - 3.0 * pi);
    let ring = spandrel(r) * 2.0 * pi * centroid;
    let want = ring + spandrel(r) * 20.0 + seam_alone;
    assert!(
        (got - want).abs() <= 1e-7 * before,
        "took off {got}, the ring, the prism and the seam's round alone {want}"
    );
}

/// A rim of the `x` bore, a box edge and the piece of the bores' seam
/// that ends on the line the `x` bore's wall is split along. The rim's
/// round trims that line back; the line kept on the bore's wall, the
/// seam's round and the rim's stand apart. The rim rounds whole, the ring
/// Pappus gives about the `x` bore's axis.
#[test]
fn an_x_bore_rim_a_box_edge_and_a_seam_on_its_split_line_round_as_each_alone() {
    let mut m = Model::new();
    let solid = cross_bores(&mut m);
    let r = 1.791_896_645_077_57;
    let rim = edge(
        &m,
        &solid,
        Point::new(20.0, 10.0, 14.0),
        Point::new(20.0, 14.0, 10.0),
        Point::new(20.0, 7.171_573, 7.171_573),
    );
    let top = edge(
        &m,
        &solid,
        Point::new(0.0, 0.0, 20.0),
        Point::new(20.0, 0.0, 20.0),
        Point::new(10.0, 0.0, 20.0),
    );
    let seam = edge(
        &m,
        &solid,
        Point::new(7.162_783, 13.866_384, 11.025_218),
        Point::new(12.236_068, 14.0, 10.0),
        Point::new(10.495_290, 6.116_619, 9.041_167),
    );
    let before = volume(&m, &solid);
    let seam_alone = before - rounded(&mut m.clone(), &solid, core::slice::from_ref(&seam), r);
    let got = before - rounded(&mut m, &solid, &[rim, top, seam], r);
    let pi = core::f64::consts::PI;
    let centroid = 4.0 + r * (10.0 - 3.0 * pi) / (12.0 - 3.0 * pi);
    let ring = spandrel(r) * 2.0 * pi * centroid;
    let want = ring + spandrel(r) * 20.0 + seam_alone;
    assert!(
        (got - want).abs() <= 1e-7 * before,
        "took off {got}, the ring, the prism and the seam's round alone {want}"
    );
}
