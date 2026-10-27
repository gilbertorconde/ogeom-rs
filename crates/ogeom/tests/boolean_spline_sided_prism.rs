//! A prism whose sides are straight B-splines, cut by a plane holding its
//! extrusion direction and passing through one or two of its side edges.
//! There the plane meets each side face along the boundary edge where the
//! face ends, and the cut is the same as the prism's with line sides.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, make_box, make_edge_between, make_face, make_prism, make_vertex, make_wire,
    volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{BSplineCurve, Curve, LineCurve, PlaneSurface};
use ogeom::math::{Direction, Frame, KnotVector, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// The prism's height.
const HEIGHT: f64 = 3.0;

/// The diamond's half diagonal: its area is `2 R²`, 4.
const R: f64 = core::f64::consts::SQRT_2;

/// The diamond's corners, turned through `turn` about `z`.
fn corners(turn: f64) -> [Point; 4] {
    [0.0, 1.0, 2.0, 3.0].map(|k: f64| {
        let a = (PI / 2.0).mul_add(k, turn);
        Point::new(R * a.cos(), R * a.sin(), 0.0)
    })
}

/// The prism `HEIGHT` tall of the diamond with corners `corners(turn)`,
/// its sides lines or clamped cubic B-splines whose four control points
/// lie evenly along the side.
fn diamond_prism(model: &mut Model, turn: f64, spline: bool) -> Shape {
    let at = corners(turn);
    let v: Vec<Shape> = at.iter().map(|p| make_vertex(model, *p).shape).collect();
    let mut edges = Vec::new();
    for i in 0..4 {
        let (a, b) = (at[i], at[(i + 1) % 4]);
        let (curve, range) = if spline {
            let control: Vec<Point> = (0..4).map(|k| a + (b - a) * (f64::from(k) / 3.0)).collect();
            let knots = KnotVector::clamped_uniform(3, 4).unwrap();
            let side = BSplineCurve::new(knots, control, T).unwrap();
            (Curve::BSpline(side), (0.0, 1.0))
        } else {
            let side = LineCurve::segment(a, b, T).unwrap();
            (Curve::from(side), (0.0, a.distance(b)))
        };
        let edge = make_edge_between(model, curve, range, &v[i], &v[(i + 1) % 4], T).unwrap();
        edges.push(edge.shape);
    }
    let wire = make_wire(model, &edges, T).unwrap().shape;
    let ground = Plane::new(Frame::new(Point::ORIGIN, Direction::Z, Direction::X, T).unwrap());
    let surface = PlaneSurface::over(ground, (-10.0, 10.0), (-10.0, 10.0)).unwrap();
    let face = make_face(model, surface.into(), &[wire], T).unwrap().shape;
    make_prism(model, &face, Vector::new(0.0, 0.0, HEIGHT), T)
        .unwrap()
        .shape
}

/// A box past the prism on the left of the line through `from` along
/// `along` in `z = 0`, its face on that line's vertical plane.
fn left_of(model: &mut Model, from: Point, along: Vector) -> Shape {
    let along = Direction::new(along, T).unwrap();
    let origin = from - along.vector() * 5.0 + Vector::new(0.0, 0.0, -1.0);
    let frame = Frame::new(origin, Direction::Z, along, T).unwrap();
    make_box(model, frame, (10.0, 5.0, HEIGHT + 2.0), T)
        .unwrap()
        .shape
}

/// The volume of the common and of the cut of the prism with the box on
/// the left of the line, each checked valid and integrated exactly.
fn common_and_cut(turn: f64, spline: bool, from: Point, along: Vector) -> [f64; 2] {
    let at = format!("turn {turn}, spline sides {spline}");
    let mut out = [0.0; 2];
    for (k, keep) in ["common", "cut"].into_iter().enumerate() {
        let mut model = Model::new();
        let prism = diamond_prism(&mut model, turn, spline);
        let half = left_of(&mut model, from, along);
        let built = match keep {
            "common" => ogeom::boolean::common(&mut model, &prism, &half, T),
            _ => ogeom::boolean::cut(&mut model, &prism, &half, T),
        }
        .unwrap_or_else(|e| panic!("{at}, {keep}: {e}"))
        .shape;
        let found = check(&model, &built, T).unwrap();
        assert!(found.is_valid(), "{at}, {keep}: {found}");
        let v =
            volume_properties(&model, &built, Deflection::with_chord(1e-3).unwrap(), T).unwrap();
        assert_eq!(
            v.deflection, 0.0,
            "{at}, {keep}: integrated on the exact surfaces"
        );
        out[k] = v.mass;
    }
    out
}

/// The plane through two opposite corners halves the prism, the cut's
/// face holding the two side edges there: `2 R² HEIGHT / 2` either way,
/// the diamond turned any way about its axis.
#[test]
fn a_plane_through_opposite_side_edges_halves_a_spline_sided_prism() {
    let want = R * R * HEIGHT;
    for turn in [0.0, 0.3, PI / 4.0, 1.0, 2.5] {
        let at = corners(turn);
        for spline in [false, true] {
            let [common, cut] = common_and_cut(turn, spline, at[0], at[2] - at[0]);
            // The sides are exactly straight, the pieces bounded by planes
            // and extrusions integrated exactly: the volumes are the
            // closed form to rounding.
            for (keep, v) in [("common", common), ("cut", cut)] {
                assert!(
                    (v - want).abs() < want * 1e-9,
                    "turn {turn}, spline sides {spline}, {keep}: {v} against {want}"
                );
            }
        }
    }
}

/// The plane through one corner and the middle of a far side: the cut's
/// face holds the side edge at that corner only, and crosses the far
/// side's face inside it. The triangle it cuts off on the right, between
/// the corner, the next corner and the middle of the side beyond, has a
/// quarter of the diamond's area: `R² / 2`.
#[test]
fn a_plane_through_one_side_edge_cuts_a_spline_sided_prism_as_with_lines() {
    let right = R * R / 2.0 * HEIGHT;
    let left = 2.0 * R * R * HEIGHT - right;
    for turn in [0.0, 0.3, 1.0] {
        let at = corners(turn);
        let middle = Point::from_vector((at[1].to_vector() + at[2].to_vector()) * 0.5);
        for spline in [false, true] {
            let [common, cut] = common_and_cut(turn, spline, at[0], middle - at[0]);
            for (keep, v, want) in [("common", common, left), ("cut", cut, right)] {
                assert!(
                    (v - want).abs() < want * 1e-9,
                    "turn {turn}, spline sides {spline}, {keep}: {v} against {want}"
                );
            }
        }
    }
}
