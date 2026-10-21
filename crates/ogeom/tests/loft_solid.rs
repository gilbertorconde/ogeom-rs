//! Solid lofts closed from loft sheets: sections of different edge counts
//! and curved sections, ruled exactly, and lofts that follow guide curves,
//! measured against their closed forms.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_edge, make_polygon, project_on_surface, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve, Curve3d as _, LineCurve, SurfaceGeometry};
use ogeom::math::{Circle, Direction, Frame, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::offset::make_loft_solid;
use ogeom::topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// A full circle of `radius` about the `z` axis at height `z`, starting in
/// the direction `start` from its centre.
fn circle(model: &mut Model, radius: f64, z: f64, start: Vector) -> Shape {
    let frame = Frame::new(
        Point::new(0.0, 0.0, z),
        Direction::Z,
        Direction::new(start, T).unwrap(),
        T,
    )
    .unwrap();
    let circle = Circle::new(frame, radius, T).unwrap();
    make_edge(model, CircleCurve::new(circle).into(), (0.0, 2.0 * PI), T)
        .unwrap()
        .shape
}

/// A closed polygon through `corners` at height `z`.
fn polygon(model: &mut Model, corners: &[(f64, f64)], z: f64) -> Shape {
    let points: Vec<Point> = corners.iter().map(|&(x, y)| Point::new(x, y, z)).collect();
    make_polygon(model, &points, true, T).unwrap().shape
}

/// The square of side `2 a` about the `z` axis at height `z`, from its
/// corner at `(a, -a)` counter-clockwise.
fn square(model: &mut Model, a: f64, z: f64) -> Shape {
    polygon(model, &[(a, -a), (a, a), (-a, a), (-a, -a)], z)
}

/// The solid's volume. Checked to 1e-9 of the closed form: far looser than
/// the integration over the exact walls lands (below 1e-12 here), far
/// tighter than pairing the arcs by angle instead of by their parameter,
/// which moves the square-to-circle volume by 1.5e-3.
fn volume(model: &Model, solid: &Shape) -> f64 {
    volume_properties(
        model,
        solid,
        Deflection {
            chord: 1e-3,
            ..Deflection::default()
        },
        T,
    )
    .unwrap()
    .mass
}

fn surfaces(model: &Model, shape: &Shape) -> Vec<SurfaceGeometry> {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .iter()
        .map(|face| {
            let Some(NodeData::Face(data)) = model.node(face).map(|n| n.data()) else {
                panic!("a face holds no face data");
            };
            model.geometry().surface(data.surface).unwrap().clone()
        })
        .collect()
}

/// The ruled solid from the square of side `2 a` to the circle of radius
/// `r` a height `h` above it, each side paired with the quarter arc facing
/// it.
///
/// The section at height `t h` is `(1 - t) s + t c` of the side's point `s`
/// and the arc's point `c` at one parameter, so its area is
/// `(1 - t)² A₀ + t² A₁ + t (1 - t) K` with `K = ½ ∮ (s × c' + c × s')`,
/// and the volume `h (A₀ / 3 + A₁ / 3 + K / 6)`. Each quarter arc is the
/// rational quadratic with weights `1, √2 / 2, 1`, whose angle from the
/// arc's middle is `φ = 2 atan(k u)` at `u = 2 t - 1`, `k = tan(π / 8)`;
/// the side runs `(a, a u)`. With `tan(φ / 2) sin φ = 1 - cos φ` and
/// `cos φ / cos²(φ / 2) = 2 - sec²(φ / 2)` the two integrals along one side
/// come to `a r (√2 + (π / 2 - √2) / k)` and `a r (π / (2 k) - 2)`.
fn square_to_circle_volume(a: f64, r: f64, h: f64) -> f64 {
    let k = (PI / 8.0).tan();
    let side = a * r * (2.0_f64.sqrt() - 2.0 + (PI - 2.0_f64.sqrt()) / k);
    let mixed = 2.0 * side;
    h * (4.0 * a * a / 3.0 + PI * r * r / 3.0 + mixed / 6.0)
}

#[test]
fn a_ruled_solid_from_a_square_to_a_circle_has_the_closed_form_volume() {
    let (a, r, h) = (2.0, 3.0, 5.0);
    let mut model = Model::new();
    let bottom = square(&mut model, a, 0.0);
    let top = circle(&mut model, r, h, Vector::new(1.0, -1.0, 0.0));
    let solid = make_loft_solid(&mut model, &[bottom, top], &[], true, T)
        .unwrap()
        .shape;
    assert_eq!(model.kind_of(&solid).unwrap(), ShapeType::Solid);
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let all = surfaces(&model, &solid);
    assert_eq!(all.len(), 6, "four walls and two caps");
    assert_eq!(
        all.iter()
            .filter(|s| matches!(s, SurfaceGeometry::Plane(_)))
            .count(),
        2,
        "the walls between a side and an arc are no plane"
    );
    let expected = square_to_circle_volume(a, r, h);
    let measured = volume(&model, &solid);
    assert!(
        (measured - expected).abs() < expected * 1e-9,
        "{measured} against {expected}"
    );
}

/// The same loft through a third section, the square again above the
/// circle: two ruled spans meeting on the circle, the volume twice one
/// span's.
#[test]
fn a_ruled_solid_through_three_mixed_sections_is_its_two_spans() {
    let (a, r, h) = (2.0, 3.0, 5.0);
    let mut model = Model::new();
    let bottom = square(&mut model, a, 0.0);
    let middle = circle(&mut model, r, h, Vector::new(1.0, -1.0, 0.0));
    let top = square(&mut model, a, 2.0 * h);
    let solid = make_loft_solid(&mut model, &[bottom, middle, top], &[], true, T)
        .unwrap()
        .shape;
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_valid(), "{found}");
    assert_eq!(surfaces(&model, &solid).len(), 10, "eight walls, two caps");
    let expected = 2.0 * square_to_circle_volume(a, r, h);
    let measured = volume(&model, &solid);
    assert!(
        (measured - expected).abs() < expected * 1e-9,
        "{measured} against {expected}"
    );
}

/// A square of four edges lofted to the same square drawn with eight
/// (each side broken at its midpoint) is the prism: every wall a plane
/// and the volume the square's area times the height.
#[test]
fn a_loft_between_one_square_drawn_with_different_edge_counts_is_the_prism() {
    let (a, h) = (2.0, 3.0);
    for ruled in [true, false] {
        let mut model = Model::new();
        let bottom = square(&mut model, a, 0.0);
        let top = polygon(
            &mut model,
            &[
                (a, -a),
                (a, 0.0),
                (a, a),
                (0.0, a),
                (-a, a),
                (-a, 0.0),
                (-a, -a),
                (0.0, -a),
            ],
            h,
        );
        let solid = make_loft_solid(&mut model, &[bottom, top], &[], ruled, T)
            .unwrap()
            .shape;
        let found = check(&model, &solid, T).unwrap();
        assert!(found.is_valid(), "ruled {ruled}: {found}");
        let all = surfaces(&model, &solid);
        assert_eq!(all.len(), 10, "ruled {ruled}: eight walls and two caps");
        if ruled {
            assert!(
                all.iter().all(|s| matches!(s, SurfaceGeometry::Plane(_))),
                "every ruled wall of a prism is a plane"
            );
        }
        let expected = 4.0 * a * a * h;
        let measured = volume(&model, &solid);
        assert!(
            (measured - expected).abs() < expected * 1e-9,
            "ruled {ruled}: {measured} against {expected}"
        );
    }
}

/// A straight edge from `a` to `b`.
fn segment(model: &mut Model, a: Point, b: Point) -> Shape {
    let curve = Curve::Line(LineCurve::segment(a, b, T).unwrap());
    let domain = curve.domain();
    make_edge(model, curve, domain, T).unwrap().shape
}

/// Three coaxial circles on one cone, lofted along two of the cone's
/// rulings as guides: the guides lie on the unguided skin already, so the
/// guided solid is the frustum, `π h (R² + R r + r²) / 3`.
#[test]
fn a_guided_solid_through_circles_on_a_cone_along_its_rulings_is_the_frustum() {
    let (r0, r1, h) = (3.0, 4.0, 10.0);
    let mut model = Model::new();
    let x = Vector::new(1.0, 0.0, 0.0);
    let sections: Vec<Shape> = [0.0, 0.5, 1.0]
        .iter()
        .map(|&t| circle(&mut model, r0 + (r1 - r0) * t, h * t, x))
        .collect();
    let guides = [
        segment(&mut model, Point::new(r0, 0.0, 0.0), Point::new(r1, 0.0, h)),
        segment(
            &mut model,
            Point::new(-r0, 0.0, 0.0),
            Point::new(-r1, 0.0, h),
        ),
    ];
    let solid = make_loft_solid(&mut model, &sections, &guides, false, T)
        .unwrap()
        .shape;
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let expected = PI * h * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0;
    let measured = volume(&model, &solid);
    assert!(
        (measured - expected).abs() < expected * 1e-9,
        "{measured} against {expected}"
    );
}

/// A drum lofted along a guide bowed out from its side: the wall follows
/// the guide, the caps stay the two circles, and the solid holds more
/// than the drum.
#[test]
fn a_guided_solid_follows_a_bowed_guide_between_its_caps() {
    let (r, h, bow) = (3.0, 10.0, 1.0);
    let mut model = Model::new();
    let x = Vector::new(1.0, 0.0, 0.0);
    let sections = [circle(&mut model, r, 0.0, x), circle(&mut model, r, h, x)];
    let guide_circle = Circle::through(
        Point::new(r, 0.0, 0.0),
        Point::new(r + bow, 0.0, h / 2.0),
        Point::new(r, 0.0, h),
        T,
    )
    .unwrap();
    let end = {
        let f = guide_circle.frame();
        let d = Point::new(r, 0.0, h) - f.origin();
        let t = d.dot(f.y().vector()).atan2(d.dot(f.x().vector()));
        if t < 0.0 { t + 2.0 * PI } else { t }
    };
    let guide = make_edge(
        &mut model,
        CircleCurve::new(guide_circle).into(),
        (0.0, end),
        T,
    )
    .unwrap()
    .shape;
    let solid = make_loft_solid(&mut model, &sections, &[guide], false, T)
        .unwrap()
        .shape;
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let all = surfaces(&model, &solid);
    assert_eq!(all.len(), 3, "one wall and two caps");
    let wall = all
        .iter()
        .find(|s| !matches!(s, SurfaceGeometry::Plane(_)))
        .unwrap();
    let f = guide_circle.frame();
    for i in 0..=32 {
        let t = end * f64::from(i) / 32.0;
        let p = f.origin()
            + (f.x().vector() * t.cos() + f.y().vector() * t.sin()) * guide_circle.radius();
        let gap = project_on_surface(wall, p, 24, T).unwrap().distance;
        // The guided sheet holds its guides to ten confusions.
        assert!(
            gap < 10.0 * T.confusion(),
            "the guide is {gap:.3e} off the wall at {p:?}"
        );
    }
    let drum = PI * r * r * h;
    assert!(
        volume(&model, &solid) > drum,
        "the bowed wall holds more than the drum"
    );
}

/// The circular arc in the half-plane through the `z` axis and `(x, y)`,
/// from `(x, y, 0)` out to `bow` farther from the axis at `h / 2` and back
/// to `(x, y, h)`.
fn bowed(model: &mut Model, x: f64, y: f64, bow: f64, h: f64) -> Shape {
    let rho = x.hypot(y);
    let out = (rho + bow) / rho;
    let arc = Circle::through(
        Point::new(x, y, 0.0),
        Point::new(x * out, y * out, h / 2.0),
        Point::new(x, y, h),
        T,
    )
    .unwrap();
    let f = arc.frame();
    let d = Point::new(x, y, h) - f.origin();
    let end = d.dot(f.y().vector()).atan2(d.dot(f.x().vector()));
    let end = if end < 0.0 { end + 2.0 * PI } else { end };
    make_edge(model, CircleCurve::new(arc).into(), (0.0, end), T)
        .unwrap()
        .shape
}

/// The rectangle `2 a` by `2 b` about the `z` axis at height `z`, scaled
/// by `s`, from its corner at `(a, -b)` counter-clockwise.
fn rectangle(model: &mut Model, a: f64, b: f64, s: f64, z: f64) -> Shape {
    let (a, b) = (a * s, b * s);
    polygon(model, &[(a, -b), (a, b), (-a, b), (-a, -b)], z)
}

/// A rectangle lofted along four guides through its corners, each the arc
/// bowed `bow` out from the axis in its corner's half-plane, through two
/// sections or three (the third the rectangle scaled to the arcs' middle).
///
/// Each wall is the Gordon patch between two guides that are one arc
/// turned about the axis, over sections that are straight: the skin and
/// the corner interpolation cancel, and the wall is `(1 - u) g₀ + u g₁`,
/// so the section at height `z` is the rectangle scaled by `ρ(z) / ρ₀`,
/// `ρ` the arc's distance from the axis and `ρ₀ = √(a² + b²)`. The volume
/// is `4 a b ∫ ρ² dz / ρ₀²`, with `ρ = c + √(R² - w²)` about the arc's
/// centre `c` at `w = z - h / 2`.
#[test]
fn a_guided_solid_through_rectangles_follows_a_guide_at_every_corner() {
    let (a, b, h, bow): (f64, f64, f64, f64) = (2.0, 1.0, 6.0, 0.75);
    let rho0 = a.hypot(b);
    let radius = (bow * bow + h * h / 4.0) / (2.0 * bow);
    let centre = rho0 + bow - radius;
    let root = (h / 2.0) * (radius * radius - h * h / 4.0).sqrt()
        + radius * radius * (h / (2.0 * radius)).asin();
    let squared = (centre * centre + radius * radius) * h - h.powi(3) / 12.0 + 2.0 * centre * root;
    let expected = 4.0 * a * b * squared / (rho0 * rho0);
    for middle in [false, true] {
        let mut model = Model::new();
        let mut sections = vec![rectangle(&mut model, a, b, 1.0, 0.0)];
        if middle {
            let s = (rho0 + bow) / rho0;
            sections.push(rectangle(&mut model, a, b, s, h / 2.0));
        }
        sections.push(rectangle(&mut model, a, b, 1.0, h));
        let guides: Vec<Shape> = [(a, -b), (a, b), (-a, b), (-a, -b)]
            .iter()
            .map(|&(x, y)| bowed(&mut model, x, y, bow, h))
            .collect();
        let solid = make_loft_solid(&mut model, &sections, &guides, false, T)
            .unwrap()
            .shape;
        let found = check(&model, &solid, T).unwrap();
        assert!(found.is_valid(), "middle {middle}: {found}");
        assert_eq!(
            surfaces(&model, &solid).len(),
            6,
            "middle {middle}: four walls and two caps"
        );
        let measured = volume(&model, &solid);
        assert!(
            (measured - expected).abs() < expected * 1e-9,
            "middle {middle}: {measured} against {expected}"
        );
    }
}

/// A rectangle lofted to a smaller one along straight guides through the
/// corners, past a third section on the guides: the pyramid's frustum,
/// `h (A₀ + A₁ + √(A₀ A₁)) / 3`.
#[test]
fn a_guided_solid_along_straight_corner_guides_is_the_frustum() {
    let (a, b, h, s) = (3.0, 2.0, 4.0, 0.5);
    let mut model = Model::new();
    let sections = [
        rectangle(&mut model, a, b, 1.0, 0.0),
        rectangle(&mut model, a, b, 0.5 * (1.0 + s), h / 2.0),
        rectangle(&mut model, a, b, s, h),
    ];
    let guides: Vec<Shape> = [(a, -b), (a, b), (-a, b), (-a, -b)]
        .iter()
        .map(|&(x, y)| {
            segment(
                &mut model,
                Point::new(x, y, 0.0),
                Point::new(x * s, y * s, h),
            )
        })
        .collect();
    let solid = make_loft_solid(&mut model, &sections, &guides, false, T)
        .unwrap()
        .shape;
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let (low, high) = (4.0 * a * b, 4.0 * a * b * s * s);
    let expected = h * (low + high + (low * high).sqrt()) / 3.0;
    let measured = volume(&model, &solid);
    assert!(
        (measured - expected).abs() < expected * 1e-9,
        "{measured} against {expected}"
    );
}

#[test]
fn a_guided_solid_through_multi_edge_sections_refuses_by_name() {
    let (a, b, h) = (2.0, 1.0, 5.0);
    let mut model = Model::new();
    let sections = [
        rectangle(&mut model, a, b, 1.0, 0.0),
        rectangle(&mut model, a, b, 1.0, h),
    ];
    let corners = [(a, -b), (a, b), (-a, b), (-a, -b)];
    let mut guides: Vec<Shape> = corners
        .iter()
        .map(|&(x, y)| segment(&mut model, Point::new(x, y, 0.0), Point::new(x, y, h)))
        .collect();
    let refused = |model: &mut Model, sections: &[Shape], guides: &[Shape]| {
        make_loft_solid(model, sections, guides, false, T)
            .unwrap_err()
            .to_string()
    };
    let err = refused(&mut model, &sections, &guides[..3]);
    assert!(
        err.contains("needs a guide through every corner where two edges meet; corner 3 has none"),
        "{err}"
    );
    let twice = segment(&mut model, Point::new(a, b, 0.0), Point::new(a, b, h));
    guides.push(twice);
    let err = refused(&mut model, &sections, &guides);
    assert!(
        err.contains("guides 1 and 4 both run through corner 1"),
        "{err}"
    );
    guides.pop();
    let eight = polygon(
        &mut model,
        &[
            (a, -b),
            (a, 0.0),
            (a, b),
            (0.0, b),
            (-a, b),
            (-a, 0.0),
            (-a, -b),
            (0.0, -b),
        ],
        h,
    );
    let err = refused(&mut model, &[sections[0].clone(), eight], &guides);
    assert!(
        err.contains("takes sections of one edge count; section 1 has 8 and section 0 has 4"),
        "{err}"
    );
    // A guide from one corner of the bottom to another of the top.
    let skew = segment(&mut model, Point::new(a, -b, 0.0), Point::new(-a, -b, h));
    guides[0] = skew;
    let err = refused(&mut model, &sections, &guides);
    assert!(
        err.contains("guide 0 crosses section 1 at another corner or edge than section 0"),
        "{err}"
    );
}

#[test]
fn a_solid_loft_refuses_an_open_section_and_a_wavy_end_by_name() {
    let mut model = Model::new();
    let bottom = square(&mut model, 2.0, 0.0);
    let open = make_polygon(
        &mut model,
        &[
            Point::new(2.0, -2.0, 3.0),
            Point::new(2.0, 2.0, 3.0),
            Point::new(-2.0, 2.0, 3.0),
        ],
        false,
        T,
    )
    .unwrap()
    .shape;
    let err = make_loft_solid(&mut model, &[bottom.clone(), open], &[], true, T).unwrap_err();
    assert!(
        err.to_string().contains("section 1 is not a closed"),
        "{err}"
    );

    let wavy = make_polygon(
        &mut model,
        &[
            Point::new(2.0, -2.0, 3.0),
            Point::new(2.0, 2.0, 4.0),
            Point::new(-2.0, 2.0, 3.0),
            Point::new(-2.0, -2.0, 4.0),
        ],
        true,
        T,
    )
    .unwrap()
    .shape;
    let err = make_loft_solid(&mut model, &[bottom, wavy], &[], true, T).unwrap_err();
    assert!(err.to_string().contains("not planar"), "{err}");
}
