//! A block whose straight side, written as a B-spline, cuts a shallow chord
//! into the round wall of a drum. The side crosses the wall twice, about
//! `2 sqrt(2 r d)` apart for a chord of sagitta `d`, and the boolean
//! sections it there as it does the same side written as a line.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{
    check, make_cylinder, make_edge_between, make_face, make_prism, make_wire, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::{BSplineCurve, Curve, LineCurve, PlaneSurface};
use ogeom::math::{Frame, KnotVector, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::topo::{Model, Shape, VertexData};

const T: Tolerances = Tolerances::millimetres();

/// The drum's radius; it stands 1 tall on `z = 0`.
const R: f64 = 20.0;

/// A prism 1 tall over the closed polygon `corners` in `z = 0`, the sides
/// listed in `spline` written as degree-1 B-splines of their lines.
fn pad(model: &mut Model, corners: &[(f64, f64)], spline: &[usize]) -> Shape {
    let pts: Vec<Point> = corners
        .iter()
        .map(|&(x, y)| Point::new(x, y, 0.0))
        .collect();
    let vs: Vec<_> = pts
        .iter()
        .map(|&p| model.add_vertex(VertexData::new(p)))
        .collect();
    let mut edges = Vec::new();
    for k in 0..pts.len() {
        let n = (k + 1) % pts.len();
        let line = Curve::Line(LineCurve::segment(pts[k], pts[n], T).unwrap());
        let len = pts[k].distance(pts[n]);
        let (c, r) = if spline.contains(&k) {
            (
                Curve::BSpline(line.to_bspline_over((0.0, len), T).unwrap()),
                (0.0, 1.0),
            )
        } else {
            (line, (0.0, len))
        };
        edges.push(
            make_edge_between(model, c, r, &vs[k], &vs[n], T)
                .unwrap()
                .shape,
        );
    }
    let wire = make_wire(model, &edges, T).unwrap().shape;
    let plane = PlaneSurface::new(Plane::new(Frame::WORLD));
    let face = make_face(model, plane.into(), &[wire], T).unwrap().shape;
    make_prism(model, &face, Vector::new(0.0, 0.0, 1.0), T)
        .unwrap()
        .shape
}

/// The area of the circular segment of sagitta `d` cut off the drum:
/// `r² (θ − sin θ) / 2` over the angle `θ` the chord subtends, where
/// `sin(θ / 4) = sqrt(d / 2r)`, a form that keeps its digits as `d` shrinks.
fn segment(d: f64) -> f64 {
    let theta = 4.0 * (d / (2.0 * R)).sqrt().asin();
    R * R * (theta - theta.sin()) / 2.0
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

/// Fuse, cut and common of the drum with a 2 x 3 block whose bottom side
/// lies `d` inside the wall: `pi r² + 6 - segment`, `pi r² - segment` and
/// `segment`, each a valid solid, the side a line or a B-spline.
#[test]
fn a_spline_side_cutting_a_shallow_chord_into_a_drum_sections_as_a_line() {
    let drum = core::f64::consts::PI * R * R;
    for d in [2e-2, 9e-3, 1e-3, 6e-4, 1e-4, 1e-5] {
        for spline in [&[][..], &[0usize][..]] {
            let at = format!("dip {d:e}, spline {}", !spline.is_empty());
            let seg = segment(d);
            for (op, want) in [
                ("fuse", drum + 6.0 - seg),
                ("cut", drum - seg),
                ("common", seg),
            ] {
                let mut m = Model::new();
                let disc = make_cylinder(&mut m, Frame::WORLD, R, 1.0, T)
                    .unwrap()
                    .shape;
                let y = R - d;
                let block = pad(
                    &mut m,
                    &[(-1.0, y), (1.0, y), (1.0, y + 3.0), (-1.0, y + 3.0)],
                    spline,
                );
                let built = match op {
                    "fuse" => ogeom::boolean::fuse(&mut m, &disc, &block, T),
                    "cut" => ogeom::boolean::cut(&mut m, &disc, &block, T),
                    _ => ogeom::boolean::common(&mut m, &disc, &block, T),
                }
                .unwrap_or_else(|e| panic!("{at}, {op}: {e}"))
                .shape;
                let found = check(&m, &built, T).unwrap();
                assert!(found.is_valid(), "{at}, {op}: {found}");
                let got = volume(&m, &built);
                assert!(
                    (got - want).abs() < 1e-6 * want,
                    "{at}, {op}: {got} against {want}"
                );
            }
        }
    }
}

/// A prism 1 tall whose bottom side is the parabola `y = c + k x²` over
/// `-1 <= x <= 1`, a quadratic Bézier, and whose top is the line `y = c + 3`.
fn cupped(model: &mut Model, c: f64, k: f64) -> Shape {
    let at = |x: f64, y: f64| Point::new(x, y, 0.0);
    let corners = [
        at(-1.0, c + k),
        at(1.0, c + k),
        at(1.0, c + 3.0),
        at(-1.0, c + 3.0),
    ];
    let vs: Vec<_> = corners
        .iter()
        .map(|&p| model.add_vertex(VertexData::new(p)))
        .collect();
    let bottom = BSplineCurve::new(
        KnotVector::clamped_uniform(2, 3).unwrap(),
        vec![corners[0], at(0.0, c - k), corners[1]],
        T,
    )
    .unwrap();
    let mut edges = vec![
        make_edge_between(model, Curve::BSpline(bottom), (0.0, 1.0), &vs[0], &vs[1], T)
            .unwrap()
            .shape,
    ];
    for i in 1..4 {
        let (p, q) = (corners[i], corners[(i + 1) % 4]);
        let line = Curve::Line(LineCurve::segment(p, q, T).unwrap());
        let edge = make_edge_between(
            model,
            line,
            (0.0, p.distance(q)),
            &vs[i],
            &vs[(i + 1) % 4],
            T,
        );
        edges.push(edge.unwrap().shape);
    }
    let wire = make_wire(model, &edges, T).unwrap().shape;
    let plane = PlaneSurface::new(Plane::new(Frame::WORLD));
    let face = make_face(model, plane.into(), &[wire], T).unwrap().shape;
    make_prism(model, &face, Vector::new(0.0, 0.0, 1.0), T)
        .unwrap()
        .shape
}

/// A rounded bottom dipping `d` into the drum, as a glyph's does: it
/// crosses the wall at `x = ±sqrt(s)`, `s` the positive root of
/// `k² s² + (1 + 2ck) s + c² − r² = 0`, and the drum and the prism overlap
/// over `∫ (sqrt(r² − x²) − c − k x²) dx` between the crossings. Common,
/// cut and fuse are valid and hold those volumes.
#[test]
fn a_rounded_spline_bottom_dipping_into_a_drum_sections_at_its_crossings() {
    let k = 0.1;
    for d in [1e-3, 6e-4, 1e-4] {
        let c = R - d;
        let (qa, qb, qc) = (k * k, (2.0 * c).mul_add(k, 1.0), c.mul_add(c, -R * R));
        let s = 2.0 * qc / (-qb - qa.mul_add(-4.0 * qc, qb * qb).sqrt());
        let x0 = s.sqrt();
        let below = |x: f64| {
            0.5 * x.mul_add(R.mul_add(R, -x * x).sqrt(), R * R * (x / R).asin())
                - c * x
                - k * x * x * x / 3.0
        };
        let overlap = below(x0) - below(-x0);
        let block = 6.0 - 2.0 * k / 3.0;
        let drum = core::f64::consts::PI * R * R;
        for (op, want) in [
            ("common", overlap),
            ("cut", drum - overlap),
            ("fuse", drum + block - overlap),
        ] {
            let at = format!("dip {d:e}, {op}");
            let mut m = Model::new();
            let disc = make_cylinder(&mut m, Frame::WORLD, R, 1.0, T)
                .unwrap()
                .shape;
            let cup = cupped(&mut m, c, k);
            let built = match op {
                "fuse" => ogeom::boolean::fuse(&mut m, &disc, &cup, T),
                "cut" => ogeom::boolean::cut(&mut m, &disc, &cup, T),
                _ => ogeom::boolean::common(&mut m, &disc, &cup, T),
            }
            .unwrap_or_else(|e| panic!("{at}: {e}"))
            .shape;
            let found = check(&m, &built, T).unwrap();
            assert!(found.is_valid(), "{at}: {found}");
            let got = volume(&m, &built);
            assert!(
                (got - want).abs() < 1e-6 * want,
                "{at}: {got} against {want}"
            );
        }
    }
}
