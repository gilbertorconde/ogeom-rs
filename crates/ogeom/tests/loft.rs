//! Ruled lofts: walls that are planes where they can be and bilinear
//! patches where they cannot.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_polygon, volume_properties};
use ogeom::core::Tolerances;
use ogeom::math::Point;
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The shoelace area of a planar polygon given by its corners.
fn area(corners: &[Point]) -> f64 {
    let n = corners.len();
    (0..n)
        .map(|i| {
            let (a, b) = (corners[i], corners[(i + 1) % n]);
            a.x * b.y - b.x * a.y
        })
        .sum::<f64>()
        .abs()
        / 2.0
}

/// A square lofted to the same square turned an eighth of a turn: every
/// wall is skew, and each is the bilinear patch between its two segments.
///
/// A ruled solid's section area is quadratic in height, so the prismoidal
/// formula, `h (A₀ + 4 A_mid + A₁) / 6`, is its volume exactly, and the
/// mid section is the polygon of the walls' midlines' midpoints.
#[test]
fn a_twisted_loft_has_bilinear_walls_and_the_prismoid_s_volume() {
    let mut model = Model::new();
    let square = |model: &mut Model, turn: f64, z: f64| -> (Vec<Point>, ogeom::topo::Shape) {
        let corners: Vec<Point> = (0..4)
            .map(|i| {
                let angle =
                    turn + std::f64::consts::FRAC_PI_2 * f64::from(i) + std::f64::consts::FRAC_PI_4;
                Point::new(3.0 * angle.cos(), 3.0 * angle.sin(), z)
            })
            .collect();
        let wire = make_polygon(model, &corners, true, T).unwrap().shape;
        (corners, wire)
    };
    let (low, bottom) = square(&mut model, 0.0, 0.0);
    let (high, top) = square(&mut model, std::f64::consts::FRAC_PI_4, 5.0);
    let solid = ogeom::offset::make_loft(&mut model, &bottom, &top, T)
        .unwrap()
        .shape;

    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_usable(), "{found}");
    let faces = explore_unique(&model, &solid, ShapeType::Face).unwrap();
    assert_eq!(faces.len(), 6, "four walls and two caps");

    let mid: Vec<Point> = low.iter().zip(&high).map(|(a, b)| a.midpoint(*b)).collect();
    let expected = 5.0 * (area(&low) + 4.0 * area(&mid) + area(&high)) / 6.0;
    let fine = Deflection {
        chord: 1e-3,
        ..Deflection::default()
    };
    let measured = volume_properties(&model, &solid, fine, T).unwrap().mass;
    assert!(
        (measured - expected).abs() < expected * 2e-3,
        "the prismoid's volume: {measured} against {expected}"
    );
}
