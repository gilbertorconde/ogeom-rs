//! Solids cut where a border of one face's chart lies in the other face's
//! surface: a plane or a drum through the side edges of a prism whose sides
//! are trimmed curves or B-splines, and a plane on the rim between two
//! revolved spline faces. There the surfaces meet along the border, and no
//! closed form or iso-line path states it.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, make_box, make_cylinder, make_edge_between, make_face, make_prism, make_revolution,
    make_vertex, make_wire, transformed, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{BSplineCurve, Curve, EllipseCurve, LineCurve, PlaneSurface, TrimmedCurve};
use ogeom::math::{Axis, Direction, Ellipse, Frame, KnotVector, Plane, Point, Transform, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// The prism's height.
const HEIGHT: f64 = 3.0;

/// The diamond's half diagonal: its area is `2 R²`, 4.
const R: f64 = core::f64::consts::SQRT_2;

/// How far a half-ellipse side bulges out of the diamond.
const BULGE: f64 = 0.4;

/// The diamond's corners, turned through `turn` about `z`.
fn corners(turn: f64) -> [Point; 4] {
    [0.0, 1.0, 2.0, 3.0].map(|k: f64| {
        let a = (PI / 2.0).mul_add(k, turn);
        Point::new(R * a.cos(), R * a.sin(), 0.0)
    })
}

/// The clamped cubic B-spline whose four control points lie evenly from
/// `a` to `b`: the segment, on `[0, 1]`.
fn straight_spline(a: Point, b: Point) -> Curve {
    let control: Vec<Point> = (0..4).map(|k| a + (b - a) * (f64::from(k) / 3.0)).collect();
    let knots = KnotVector::clamped_uniform(3, 4).unwrap();
    Curve::BSpline(BSplineCurve::new(knots, control, T).unwrap())
}

/// How a side from one corner to the next is spelled.
#[derive(Debug, Clone, Copy)]
enum Side {
    /// A straight B-spline.
    Spline,
    /// The middle half of a straight B-spline twice the side's length,
    /// trimmed.
    TrimmedSpline,
    /// Half an ellipse on the side, bulging `BULGE` outward, trimmed.
    TrimmedEllipse,
}

/// The side from `a` to `b` of a diamond running anticlockwise, and its
/// range.
fn side(a: Point, b: Point, kind: Side) -> (Curve, (f64, f64)) {
    let d = b - a;
    match kind {
        Side::Spline => (straight_spline(a, b), (0.0, 1.0)),
        Side::TrimmedSpline => {
            let whole = straight_spline(a - d * 0.5, b + d * 0.5);
            (
                TrimmedCurve::new(whole, 0.25, 0.75, T).unwrap().into(),
                (0.25, 0.75),
            )
        }
        Side::TrimmedEllipse => {
            // From `a` at angle 0 to `b` at a half turn, through the
            // frame's `y`: `z` crossed with the side's reverse, outward.
            let middle = Point::from_vector((a.to_vector() + b.to_vector()) * 0.5);
            let frame =
                Frame::new(middle, Direction::Z, Direction::new(a - b, T).unwrap(), T).unwrap();
            let ellipse = Ellipse::new(frame, d.magnitude() / 2.0, BULGE, T).unwrap();
            let half = TrimmedCurve::new(EllipseCurve::new(ellipse).into(), 0.0, PI, T).unwrap();
            (half.into(), (0.0, PI))
        }
    }
}

/// The prism `HEIGHT` tall of the diamond with corners `corners(turn)`.
fn diamond_prism(model: &mut Model, turn: f64, kind: Side) -> Shape {
    let at = corners(turn);
    let v: Vec<Shape> = at.iter().map(|p| make_vertex(model, *p).shape).collect();
    let mut edges = Vec::new();
    for i in 0..4 {
        let (curve, range) = side(at[i], at[(i + 1) % 4], kind);
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

/// The volumes of the common and of the cut of the two solids `build`
/// makes, each checked valid and integrated on the exact surfaces, held to
/// `want` within `bound`.
fn common_and_cut(at: &str, build: &dyn Fn(&mut Model) -> [Shape; 2], want: [f64; 2], bound: f64) {
    for (k, keep) in ["common", "cut"].into_iter().enumerate() {
        let mut model = Model::new();
        let [part, tool] = build(&mut model);
        let built = match keep {
            "common" => ogeom::boolean::common(&mut model, &part, &tool, T),
            _ => ogeom::boolean::cut(&mut model, &part, &tool, T),
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
        assert!(
            (v.mass - want[k]).abs() <= bound,
            "{at}, {keep}: {} against {}",
            v.mass,
            want[k]
        );
    }
}

/// The plane through two opposite corners halves the prism whatever its
/// sides are spelled as, the diamond turned any way about its axis: the
/// plane meets each side face along the ruling where its trimmed profile
/// starts or ends. Half of `2 R²` and of four half-ellipses on sides
/// `√2 R` long, `π (R / √2) BULGE / 2` each, times the height.
#[test]
fn a_plane_through_opposite_side_edges_halves_a_prism_with_trimmed_sides() {
    // A section fitted within the confusion distance moves a volume by at
    // most that times the area of the face it bounds, here the cut's
    // rectangle `2 R` by the height.
    let bound = T.confusion() * 2.0 * R * HEIGHT;
    for kind in [Side::TrimmedSpline, Side::TrimmedEllipse] {
        let area = match kind {
            Side::TrimmedEllipse => {
                (2.0 * PI * R / core::f64::consts::SQRT_2).mul_add(BULGE, 2.0 * R * R)
            }
            _ => 2.0 * R * R,
        };
        let half = area * HEIGHT / 2.0;
        for turn in [0.0, 0.3, PI / 4.0, 1.0, 2.5] {
            let at = corners(turn);
            common_and_cut(
                &format!("{kind:?} sides, turn {turn}"),
                &|model| {
                    [
                        diamond_prism(model, turn, kind),
                        left_of(model, at[0], at[2] - at[0]),
                    ]
                },
                [half, half],
                bound,
            );
        }
    }
}

/// A drum about the point twice out along the first corner, through the
/// second and fourth corners, cuts the spline-sided prism through those two
/// side edges, the diamond turned any way about its axis. The drum's
/// radius is `ρ = √5 R` and its axis `c = 2 R` from the chord between
/// the two corners; inside it lie the triangle on the first corner, `R²`,
/// and the circular segment beyond the chord, `ρ² acos(c / ρ) - c R`,
/// which stays inside the diamond's far half as `c > R` keeps the arc
/// steeper than the sides at the corners.
#[test]
fn a_drum_through_two_side_edges_cuts_a_spline_sided_prism() {
    let c = 2.0 * R;
    let rho = c.hypot(R);
    let half_angle = (c / rho).acos();
    let segment = (rho * rho).mul_add(half_angle, -c * R);
    let inside = (R * R + segment) * HEIGHT;
    let outside = 2.0 * R * R * HEIGHT - inside;
    // The confusion distance times the area of the drum's face inside.
    let bound = T.confusion() * 2.0 * rho * half_angle * HEIGHT;
    for turn in [0.0, 0.3, PI / 4.0, 1.0, 2.5] {
        let at = corners(turn);
        common_and_cut(
            &format!("turn {turn}"),
            &|model| {
                let centre =
                    Point::from_vector(at[0].to_vector() * 2.0) + Vector::new(0.0, 0.0, -1.0);
                let frame = Frame::new(centre, Direction::Z, Direction::X, T).unwrap();
                [
                    diamond_prism(model, turn, Side::Spline),
                    make_cylinder(model, frame, rho, HEIGHT + 2.0, T)
                        .unwrap()
                        .shape,
                ]
            },
            [inside, outside],
            bound,
        );
    }
}

/// A spindle: the triangle from `(0.5, 0)` out to `(2, 1.5)` and back in
/// to `(0.5, 3)` in the `xz` plane, its two outer sides straight B-splines,
/// turned through `sweep` about `z`. Its two spline faces meet on the rim
/// of radius 2 at height 1.5.
fn spindle(model: &mut Model, sweep: f64) -> Shape {
    let p = [
        Point::new(0.5, 0.0, 0.0),
        Point::new(2.0, 0.0, 1.5),
        Point::new(0.5, 0.0, 3.0),
    ];
    let v: Vec<Shape> = p.iter().map(|q| make_vertex(model, *q).shape).collect();
    let mut edges = Vec::new();
    for i in 0..3 {
        let (a, b) = (p[i], p[(i + 1) % 3]);
        let (curve, range) = if i < 2 {
            (straight_spline(a, b), (0.0, 1.0))
        } else {
            (
                Curve::from(LineCurve::segment(a, b, T).unwrap()),
                (0.0, a.distance(b)),
            )
        };
        let edge = make_edge_between(model, curve, range, &v[i], &v[(i + 1) % 3], T).unwrap();
        edges.push(edge.shape);
    }
    let wire = make_wire(model, &edges, T).unwrap().shape;
    let upright = Plane::new(Frame::new(Point::ORIGIN, Direction::Y, Direction::Z, T).unwrap());
    let surface = PlaneSurface::over(upright, (-10.0, 10.0), (-10.0, 10.0)).unwrap();
    let face = make_face(model, surface.into(), &[wire], T).unwrap().shape;
    let axis = Axis {
        location: Point::ORIGIN,
        direction: Direction::Z,
    };
    make_revolution(model, &face, axis, sweep, T).unwrap().shape
}

/// The plane `z = 1.5` holds the rim between the spindle's two spline
/// faces and meets them nowhere else: it halves the spindle, turned whole
/// or three quarters, and leaned on a skew axis with the plane. By Pappus
/// each half is `sweep` times `∫ ((0.5 + z)² - 0.25) / 2 dz` over
/// `[0, 1.5]`, `1.125 sweep`.
#[test]
fn a_plane_on_the_rim_between_two_revolved_spline_faces_halves_them() {
    for sweep in [2.0 * PI, 1.5 * PI] {
        let half = 1.125 * sweep;
        // The confusion distance times the area of the annular sector the
        // plane cuts.
        let bound = T.confusion() * sweep / 2.0 * (4.0 - 0.25);
        for lean in [0.0, 0.4] {
            common_and_cut(
                &format!("sweep {sweep}, lean {lean}"),
                &|model| {
                    let part = spindle(model, sweep);
                    let frame =
                        Frame::new(Point::new(-5.0, -5.0, 1.5), Direction::Z, Direction::X, T)
                            .unwrap();
                    let above = make_box(model, frame, (10.0, 10.0, 5.0), T).unwrap().shape;
                    let skew = Axis {
                        location: Point::new(0.3, -0.2, 0.7),
                        direction: Direction::new(Vector::new(1.0, 0.5, 0.2), T).unwrap(),
                    };
                    let turn = Transform::rotation(skew, lean);
                    [
                        transformed(model, &part, turn).unwrap().shape,
                        transformed(model, &above, turn).unwrap().shape,
                    ]
                },
                [half, half],
                bound,
            );
        }
    }
}
