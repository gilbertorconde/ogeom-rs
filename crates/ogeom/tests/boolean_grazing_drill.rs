//! Drills whose wall grazes an edge of the part: lying on a frustum's base
//! with the rim's circle touching the line it rests on, lying along a
//! rounded box's tangent line into the corners where the rounds' end arcs
//! touch it, and an oblique drill through a vertex of a shaved cube with
//! its own seam on that vertex. Each cut, common and drill less the part is
//! valid, cut and common add up to the part, the drill less the part and
//! the common add up to the drill, and the common matches the drill's
//! share of the part integrated ray by ray along the drill's axis.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
#![allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    reason = "grid indices"
)]

use ogeom::algo::{check, make_box, make_cone, make_cylinder, volume_properties};
use ogeom::core::Tolerances;
use ogeom::math::{Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

/// A drill of `radius` and `length` from `start` along `axis`, its frame's
/// x along `across`, cut from the part, taken in common with it, and cut
/// by it: all valid, the volumes adding up to the part's and the drill's
/// within a millionth. Returns the common's volume.
fn drill(
    model: &mut Model,
    part: &Shape,
    start: Point,
    (axis, across): (Vector, Vector),
    radius: f64,
    length: f64,
) -> f64 {
    let frame = Frame::new(
        start,
        Direction::new(axis, T).unwrap(),
        Direction::new(across, T).unwrap(),
        T,
    )
    .unwrap();
    let tool = make_cylinder(model, frame, radius, length, T)
        .unwrap()
        .shape;
    let whole = volume(model, part);
    let mut volumes = Vec::new();
    for (name, made) in [
        ("cut", ogeom::boolean::cut(model, part, &tool, T)),
        ("common", ogeom::boolean::common(model, part, &tool, T)),
        ("rest", ogeom::boolean::cut(model, &tool, part, T)),
    ] {
        let made = made.unwrap_or_else(|e| panic!("{name} of r {radius}: {e}"));
        let diagnosis = check(model, &made.shape, T).unwrap();
        assert!(diagnosis.is_valid(), "{name} of r {radius}: {diagnosis}");
        volumes.push(volume(model, &made.shape));
    }
    let (cut, common, rest) = (volumes[0], volumes[1], volumes[2]);
    assert!(
        (cut + common - whole).abs() < whole * 1e-6,
        "r {radius}: cut {cut} + common {common} against the part's {whole}"
    );
    let bore = core::f64::consts::PI * radius * radius * length;
    assert!(
        (rest + common - bore).abs() < bore * 1e-6,
        "r {radius}: rest {rest} + common {common} against the drill's {bore}"
    );
    common
}

/// The integral over the drill's disc of the length `inside(a, b)` the part
/// holds along the axis ray at `(a, b)` across it, by the midpoint rule in
/// polar coordinates. Where the length has a square root's edge at the
/// part's boundary the rule converges as the ring width to the power 1.5:
/// a thousand rings stand within 4e-7 of the limit on every drill here.
fn across_disc(radius: f64, inside: impl Fn(f64, f64) -> f64) -> f64 {
    let (rings, spokes) = (1000, 4000);
    let (dr, dt) = (
        radius / f64::from(rings),
        core::f64::consts::TAU / f64::from(spokes),
    );
    let mut total = 0.0;
    for i in 0..rings {
        let r = (f64::from(i) + 0.5) * dr;
        for j in 0..spokes {
            let t = (f64::from(j) + 0.123) * dt;
            total += inside(r * t.cos(), r * t.sin()) * r * dr * dt;
        }
    }
    total
}

/// A frustum of radius 6 at its base and 3 at its top, 10 high, drilled
/// along y by drills lying on its base whose lowest line touches the base's
/// circle at its seam vertex. The section of each with the frustum's wall
/// touches the circle there to the fourth order. The common is the
/// frustum's chord along y, twice the root of the radius squared less x
/// squared at each height, integrated over the drill's disc.
#[test]
fn drills_lying_on_a_frustum_s_base_through_its_rim_cut() {
    for (radius, across) in [
        (0.408_601_624_444_871_8, Vector::new(0.0, 0.0, -1.0)),
        (1.003_123_901_254_259_9, Vector::new(-1.0, 0.0, 0.0)),
    ] {
        let mut model = Model::new();
        let part = make_cone(&mut model, Frame::WORLD, 6.0, 3.0, 10.0, T)
            .unwrap()
            .shape;
        let length = 15.939_543_389_761_981;
        let common = drill(
            &mut model,
            &part,
            Point::new(6.0, -length * 0.5, radius),
            (Vector::new(0.0, 1.0, 0.0), across),
            radius,
            length,
        );
        // Across the disc, a runs along x and b along z.
        let expected = across_disc(radius, |a, b| {
            let (x, z) = (6.0 + a, radius + b);
            let rim = 6.0 - 0.3 * z;
            if z <= 0.0 || z >= 10.0 || x.abs() >= rim {
                0.0
            } else {
                2.0 * (rim * rim - x * x).sqrt()
            }
        });
        assert!(
            (common - expected).abs() < expected * 1e-6,
            "r {radius}: common {common} against {expected} integrated"
        );
    }
}

/// A box rounded on every edge with radius 2, drilled along x by drills
/// whose lowest line runs along the tangent line of the round under the
/// front face and of the one under the back face. At each end the round's
/// arc where it meets the vertical round touches that line, to the fourth
/// order, at the corner. The rounded box is the inner box grown by the
/// radius, so its chord along x at a height and depth a distance d off the
/// inner box's section is the inner length and twice the root of 4 less d
/// squared.
#[test]
fn drills_along_a_round_s_tangent_line_into_the_corners_cut() {
    for (depth, radius) in [
        (0.0, 1.377_676_242_966_445_8),
        (14.0, 1.178_118_233_021_239_2),
    ] {
        let mut model = Model::new();
        let block = make_box(&mut model, Frame::WORLD, (20.0, 14.0, 10.0), T)
            .unwrap()
            .shape;
        let edges = explore_unique(&model, &block, ShapeType::Edge).unwrap();
        let part = ogeom::fillet::fillet_edges(&mut model, &block, &edges, 2.0, T)
            .unwrap()
            .shape;
        let common = drill(
            &mut model,
            &part,
            Point::new(-10.638_181_458_367_21, depth, 2.0 + radius),
            (Vector::new(1.0, 0.0, 0.0), Vector::new(0.0, -1.0, 0.0)),
            radius,
            41.276_362_916_734_42,
        );
        // Across the disc, a runs along y and b along z.
        let expected = across_disc(radius, |a, b| {
            let (y, z) = (depth + a, 2.0 + radius + b);
            let dy = (2.0 - y).max(y - 12.0).max(0.0);
            let dz = (2.0 - z).max(z - 8.0).max(0.0);
            let d2 = dy * dy + dz * dz;
            if d2 >= 4.0 {
                0.0
            } else {
                16.0 + 2.0 * (4.0 - d2).sqrt()
            }
        });
        assert!(
            (common - expected).abs() < expected * 1e-6,
            "r {radius}: common {common} against {expected} integrated"
        );
    }
}

/// A cube shaved round by a drum about its vertical centre line, drilled
/// obliquely with the drill's wall, and its seam, through the vertex where
/// the base, a side and the drum meet. The common is the length along each
/// ray inside the cube's three slabs and the drum, integrated over the
/// drill's disc.
#[test]
fn an_oblique_drill_seamed_through_a_shaved_cube_s_vertex_cuts() {
    let mut model = Model::new();
    let block = make_box(&mut model, Frame::WORLD, (10.0, 10.0, 10.0), T)
        .unwrap()
        .shape;
    let at = Frame::new(Point::new(5.0, 5.0, -1.0), Direction::Z, Direction::X, T).unwrap();
    let drum = make_cylinder(&mut model, at, 6.0, 12.0, T).unwrap().shape;
    let part = ogeom::boolean::common(&mut model, &block, &drum, T)
        .unwrap()
        .shape;
    let radius = 0.644_061_888_054_237_5;
    let start = Point::new(
        -0.463_577_638_258_914_96,
        -0.993_330_859_774_557_6,
        -3.020_420_875_783_745_6,
    );
    let axis = Vector::new(
        0.731_040_046_129_3,
        0.633_769_617_465_885_6,
        0.252_817_172_938_105_24,
    );
    let across = Vector::new(0.0, 0.370_517_907_167_734_86, -0.928_825_322_904_173_3);
    let common = drill(
        &mut model,
        &part,
        start,
        (axis, across),
        radius,
        28.626_551_154_512_84,
    );
    let other = axis.cross(across);
    let expected = across_disc(radius, |a, b| {
        let p = start + across * a + other * b;
        let (mut lo, mut hi) = (f64::NEG_INFINITY, f64::INFINITY);
        for (at, along) in [(p.x, axis.x), (p.y, axis.y), (p.z, axis.z)] {
            let (t0, t1) = (-at / along, (10.0 - at) / along);
            lo = lo.max(t0.min(t1));
            hi = hi.min(t0.max(t1));
        }
        let (qx, qy) = (p.x - 5.0, p.y - 5.0);
        let qa = axis.x * axis.x + axis.y * axis.y;
        let qb = 2.0 * (qx * axis.x + qy * axis.y);
        let qc = qx * qx + qy * qy - 36.0;
        let disc = qb * qb - 4.0 * qa * qc;
        if disc <= 0.0 {
            return 0.0;
        }
        lo = lo.max((-qb - disc.sqrt()) / (2.0 * qa));
        hi = hi.min((-qb + disc.sqrt()) / (2.0 * qa));
        (hi - lo).max(0.0)
    });
    assert!(
        (common - expected).abs() < expected * 1e-6,
        "common {common} against {expected} integrated"
    );
}
