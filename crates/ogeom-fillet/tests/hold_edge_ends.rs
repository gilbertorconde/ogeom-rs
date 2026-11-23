//! Hold edge at the round's end faces: a leaning end face on a cylinder
//! beside the edge, where the round runs on past the face and the cut
//! trims it there; and the ends that stay with trim through, where the
//! balls past the end face cannot stand or the boolean does not resolve
//! the cut.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_core::Tolerances;
use ogeom_fillet::{Overflow, fillet_edges, fillet_edges_overflow};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The round's radius on every part.
const R: f64 = 4.0;

/// The drum's axis `(y, z)`, along `x`, and its radius: through the edge
/// `y = z = 0`, the solid inside it.
const DRUM: (f64, f64, f64) = (8.0, -6.0, 10.0);

/// The pocket in the drum part's top: centre `(x, y)` and radius.
const POCKET: (f64, f64, f64) = (15.0, 6.5, 5.0);

fn at(x: f64, y: f64, z: f64) -> Frame {
    Frame::new(Point::new(x, y, z), Direction::Z, Direction::X, T).unwrap()
}

fn edges_where(m: &Model, s: &Shape, on: impl Fn(Point) -> bool) -> Vec<Shape> {
    explore_unique(m, s, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            let Some((c, (a, b))) = ogeom_algo::placed_edge_curve(m, e, T).unwrap() else {
                return false;
            };
            (0..=10).all(|k| on(c.point_at(a + (b - a) * f64::from(k) / 10.0, T).unwrap()))
        })
        .collect()
}

fn front_edge(m: &Model, s: &Shape) -> Vec<Shape> {
    let edge = edges_where(m, s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
    assert_eq!(edge.len(), 1, "the front top edge");
    edge
}

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

fn surface_of(m: &Model, face: &Shape) -> SurfaceGeometry {
    let data = m.node(face).unwrap().data().as_face().unwrap();
    m.geometry().surface(data.surface).unwrap().clone()
}

/// `s` less the half space past the plane through `on` with outward
/// normal `out`.
fn cut_past(m: &mut Model, s: &Shape, on: Point, out: Vector) -> Shape {
    let z = Direction::new(out, T).unwrap();
    let x = Direction::new(Vector::Y, T).unwrap();
    let frame = Frame::new(on, z, x, T).unwrap();
    let corner = on - frame.x().vector() * 40.0 - frame.y().vector() * 40.0;
    let frame = Frame::new(corner, z, x, T).unwrap();
    let past = ogeom_algo::make_box(m, frame, (80.0, 80.0, 40.0), T)
        .unwrap()
        .shape;
    ogeom_bool::cut(m, s, &past, T).unwrap().shape
}

/// A part whose top `z = 0` runs from the front top edge to a back face at
/// `y = 14`, its front the drum [`DRUM`] bulging out below the edge, the
/// bottom at `z = -12`, the pocket [`POCKET`] sunk 8 into its top; its ends
/// the planes `x = k0 z` and `x = 30 + k1 z`.
fn drum(m: &mut Model, (k0, k1): (f64, f64)) -> Shape {
    let block = ogeom_algo::make_box(m, at(-10.0, -5.0, -12.0), (50.0, 19.0, 12.0), T)
        .unwrap()
        .shape;
    let (ay, az, ar) = DRUM;
    let axis = Frame::new(Point::new(-11.0, ay, az), Direction::X, Direction::Y, T).unwrap();
    let drum = ogeom_algo::make_cylinder(m, axis, ar, 52.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::common(m, &block, &drum, T).unwrap().shape;
    let (px, py, pr) = POCKET;
    let cutter = ogeom_algo::make_cylinder(m, at(px, py, -8.0), pr, 9.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(m, &s, &cutter, T).unwrap().shape;
    let s = cut_past(m, &s, Point::ORIGIN, Vector::new(-1.0, 0.0, k0));
    cut_past(
        m,
        &s,
        Point::new(30.0, 0.0, 0.0),
        Vector::new(1.0, 0.0, -k1),
    )
}

/// The first moment `∫ z dA` of the band's section on the drum part: the
/// corner between the top, the drum and the ball's arc, in polar
/// coordinates about the ball's centre, Gauss-Legendre of `n` points over
/// each side of the direction to the edge (the top's side, the drum's).
fn band_moment(n: usize) -> f64 {
    let (ay, az, ar) = DRUM;
    let rho = ar - R;
    // The band's centre: `R` under the top, `rho` from the drum's axis, on
    // the edge's side.
    let (cy, cz) = (ay - (rho * rho - (R + az).powi(2)).sqrt(), -R);
    let (qy, qz) = (cy - ay, cz - az);
    let to_drum = qz.atan2(qy);
    let to_top = PI / 2.0;
    let to_edge = (-cz).atan2(-cy);
    // How far along the ray at `th` the corner runs: to the top, else to
    // the drum, whichever comes first.
    let reach = |th: f64| {
        let (uy, uz) = (th.cos(), th.sin());
        let top = if uz > 0.0 { -cz / uz } else { f64::INFINITY };
        let half = uy * qy + uz * qz;
        let drum = -half + (half * half - (qy * qy + qz * qz - ar * ar)).sqrt();
        top.min(drum)
    };
    let g = gauss(n);
    let mut sum = 0.0;
    for (a, b) in [(to_top, to_edge), (to_edge, to_drum)] {
        let h = 0.5 * (b - a);
        for (x, w) in &g {
            let th = a + h * (x + 1.0);
            let d = reach(th);
            let inner = cz * (d * d - R * R) / 2.0 + th.sin() * (d.powi(3) - R.powi(3)) / 3.0;
            sum += w * h * inner;
        }
    }
    sum
}

/// Gauss-Legendre nodes and weights on [-1, 1], by Newton on the
/// polynomial's roots.
fn gauss(n: usize) -> Vec<(f64, f64)> {
    #[expect(clippy::cast_precision_loss, reason = "a few dozen nodes")]
    let nf = n as f64;
    (0..n)
        .map(|i| {
            #[expect(clippy::cast_precision_loss, reason = "a few dozen nodes")]
            let mut x = (PI * (i as f64 + 0.75) / (nf + 0.5)).cos();
            let mut d = 0.0;
            for _ in 0..100 {
                let (mut p0, mut p1) = (1.0, x);
                for k in 2..=n {
                    #[expect(clippy::cast_precision_loss, reason = "a few dozen nodes")]
                    let kf = k as f64;
                    let p2 = ((2.0 * kf - 1.0) * x * p1 - (kf - 1.0) * p0) / kf;
                    p0 = p1;
                    p1 = p2;
                }
                d = nf * (x * p1 - p0) / (x * x - 1.0);
                let step = p1 / d;
                x -= step;
                if step.abs() < 1e-16 {
                    break;
                }
            }
            (x, 2.0 / ((1.0 - x * x) * d * d))
        })
        .collect()
}

/// The largest angle between the normals either side of each edge a fitted
/// face shares with another fitted face, a cylinder along `x` (the band or
/// the drum), sampled along the edge, and how many such edges.
fn worst_join(m: &Model, s: &Shape) -> (f64, usize) {
    let mut worst = 0.0_f64;
    let mut joins = 0;
    for e in explore_unique(m, s, ShapeType::Edge).unwrap() {
        let faces = ancestors_of(m, s, &e, ShapeType::Face).unwrap();
        let [a, b] = &faces[..] else {
            continue;
        };
        let (sa, sb) = (surface_of(m, a), surface_of(m, b));
        let other = match (&sa, &sb) {
            (SurfaceGeometry::BSpline(_), o) | (o, SurfaceGeometry::BSpline(_)) => o,
            _ => continue,
        };
        let tangent_kind = match other {
            SurfaceGeometry::BSpline(_) => true,
            SurfaceGeometry::Cylinder(cy) => {
                cy.cylinder()
                    .axis()
                    .direction
                    .vector()
                    .cross(Vector::X)
                    .magnitude()
                    < 1e-9
            }
            _ => false,
        };
        let Some((c, (lo, hi))) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap() else {
            continue;
        };
        if !tangent_kind {
            continue;
        }
        joins += 1;
        for k in 1..20 {
            let p = c.point_at(lo + (hi - lo) * f64::from(k) / 20.0, T).unwrap();
            let normal = |g: &SurfaceGeometry| -> Vector {
                let foot = ogeom_geom::project::project_on_surface(g, p, 16, T).unwrap();
                let (u, v) = foot.parameters;
                g.normal_at(u, v, T).unwrap().vector()
            };
            let angle = normal(&sa).dot(normal(&sb)).abs().clamp(0.0, 1.0).acos();
            worst = worst.max(angle);
        }
    }
    (worst, joins)
}

/// How far the result's edges on the pocket's rim stray from its circle,
/// and the x-range they cover in front of the pocket's centre.
fn rim_kept(m: &Model, s: &Shape) -> (f64, (f64, f64)) {
    let (px, py, pr) = POCKET;
    let off = |p: Point| ((p.x - px).hypot(p.y - py) - pr).abs().max(p.z.abs());
    let mut worst = 0.0_f64;
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for e in edges_where(m, s, |p| off(p) < 1e-6) {
        let (c, (ta, tb)) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap().unwrap();
        for k in 0..=50 {
            let p = c.point_at(ta + (tb - ta) * f64::from(k) / 50.0, T).unwrap();
            worst = worst.max(off(p));
            if p.y < py {
                lo = lo.min(p.x);
                hi = hi.max(p.x);
            }
        }
    }
    (worst, (lo, hi))
}

/// The drum part's end planes' leans `(k0, k1)`: square, the far end
/// leaning back at its foot, leaning out, and both ends leaning.
const LEANS: [(f64, f64); 4] = [(0.0, 0.0), (0.0, 0.25), (0.0, -0.25), (-0.2, 0.25)];

/// A leaning end face on the drum part: the plane cuts the drum in an
/// ellipse, so the round runs on past it to where it stands clear and the
/// cut trims it by the end face. The band there is the drum's band, one
/// section all along, so leaning the end plane `x = 30 + k z` changes the
/// volume cut away by `k` times that section's first moment `∫ z dA`, and
/// leaning the near one `x = k z` by `-k` times it; the square part is the
/// one whose volume `tests/hold_edge_curved.rs` checks against the
/// sampled reference. Each result is valid and closed, the fitted piece
/// meets the band and the drum G1, and the pocket's rim is kept on its
/// circle over every foot.
#[test]
fn hold_edge_runs_the_round_on_past_a_leaning_end_on_a_drum() {
    let (coarse, moment) = (band_moment(24), band_moment(48));
    let bound = (moment - coarse).abs();
    eprintln!("band section's first moment {moment:.12} (bound {bound:.1e})");
    // Where the band's contact on the top, over its centre `R` under the
    // top and `ar - R` from the drum's axis, meets the pocket's rim.
    let (ay, az, ar) = DRUM;
    let (px, py, pr) = POCKET;
    let contact = ay - ((ar - R).powi(2) - (R + az).powi(2)).sqrt();
    let half = (pr * pr - (py - contact).powi(2)).sqrt();
    let feet = (px - half, px + half);
    let mut square = None;
    for (k0, k1) in LEANS {
        let mut m = Model::new();
        let s = drum(&mut m, (k0, k1));
        let edge = front_edge(&m, &s);
        let before = volume(&m, &s);
        let built = fillet_edges(&mut m, &s, &edge, R, T).unwrap();
        let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
        assert!(found.is_valid(), "{:?}", found.problems);
        for shell in explore_unique(&m, &built.shape, ShapeType::Shell).unwrap() {
            assert!(ogeom_algo::is_shell_closed(&m, &shell).unwrap(), "closed");
        }
        assert!(!built.history.generated(&edge[0]).is_empty(), "the round");
        let (crease, joins) = worst_join(&m, &built.shape);
        assert!(joins >= 2 && crease < 1e-5, "{crease} over {joins}");
        let (worst, (lo, hi)) = rim_kept(&m, &built.shape);
        assert!(worst < 1e-9, "the rim strays {worst} off its circle");
        assert!(
            lo <= feet.0 + 1e-6 && hi >= feet.1 - 1e-6,
            "the rim kept over {lo}..{hi}, the feet over {feet:?}"
        );
        let got = before - volume(&m, &built.shape);
        let Some(base) = square else {
            square = Some(got);
            continue;
        };
        let want = base + (k1 - k0) * moment;
        eprintln!(
            "ends leaning {k0}, {k1}: cut away {got:.9}, square plus moment {want:.9}, \
             error {:.1e}, worst join {crease:.1e} rad",
            (got - want).abs()
        );
        assert!((got - want).abs() < 1e-9 + bound, "{got} against {want}");
    }
}

/// The options on the leaning drum: hold edge at the end off still holds
/// the rim (the held stretch lies inside the band); hold edge off goes on
/// to trim through, which refuses this part by name, as before; with trim
/// through off too, it is refused naming hold edge.
#[test]
fn a_leaning_end_on_a_drum_follows_the_options() {
    let mut m = Model::new();
    let s = drum(&mut m, (0.0, 0.25));
    let edge = front_edge(&m, &s);
    let before = volume(&m, &s);
    let held = fillet_edges(&mut m, &s, &edge, R, T).unwrap().shape;
    let held = before - volume(&m, &held);
    let mut end_off = Overflow::default();
    end_off.hold_edge_at_end = false;
    let same = fillet_edges_overflow(&mut m, &s, &edge, R, end_off, T).unwrap();
    assert!((before - volume(&m, &same.shape) - held).abs() < 1e-9);
    let mut off = Overflow::default();
    off.hold_edge = false;
    let err = fillet_edges_overflow(&mut m, &s, &edge, R, off, T).unwrap_err();
    assert!(
        err.to_string().contains("touches outside its face"),
        "{err}"
    );
    off.trim_through = false;
    let err = fillet_edges_overflow(&mut m, &s, &edge, R, off, T).unwrap_err();
    assert!(err.to_string().contains("hold edge is off"), "{err}");
}

/// The faceted solid over `points` and `rings`, the rings taken as given
/// or every one reversed, whichever encloses a volume.
fn faceted(m: &mut Model, points: &[Point], rings: &[Vec<usize>]) -> Shape {
    for flip in [false, true] {
        let mut rings = rings.to_vec();
        if flip {
            rings.iter_mut().for_each(|r| r.reverse());
        }
        if let Ok((built, _)) = ogeom_algo::make_faceted_solid(m, points, &rings, &[], T) {
            return built.shape;
        }
    }
    panic!("the rings enclose no volume either way round");
}

/// A part 20 long whose top `z = 0` lies in front of the polyline `p0 k p1`
/// in plan (from `x = 0` to `x = 20`), a cliff 6 high falling from it to a
/// ledge, the block 12 deep and 15 tall below.
fn kinked(m: &mut Model, [(x0, y0), (xk, yk), (x1, y1)]: [(f64, f64); 3]) -> Shape {
    let p = Point::new;
    let points = [
        p(x0, 0.0, 0.0),
        p(x1, 0.0, 0.0),
        p(x1, y1, 0.0),
        p(xk, yk, 0.0),
        p(x0, y0, 0.0),
        p(x1, y1, -6.0),
        p(xk, yk, -6.0),
        p(x0, y0, -6.0),
        p(x1, 12.0, -6.0),
        p(x0, 12.0, -6.0),
        p(x0, 0.0, -15.0),
        p(x1, 0.0, -15.0),
        p(x1, 12.0, -15.0),
        p(x0, 12.0, -15.0),
    ];
    let rings = vec![
        vec![0, 1, 2, 3, 4],
        vec![10, 13, 12, 11],
        vec![0, 10, 11, 1],
        vec![1, 11, 12, 8, 5, 2],
        vec![0, 4, 7, 9, 13, 10],
        vec![3, 6, 7, 4],
        vec![2, 5, 6, 3],
        vec![5, 8, 9, 7, 6],
        vec![8, 12, 13, 9],
    ];
    faceted(m, &points, &rings)
}

/// A block 20 long under `z = 0` whose top narrows from `y0` behind the
/// front at `x = 0` to `y1` at the far end, the plane `x = 20 + k z`, a
/// cliff 6 high falling from its back edge.
fn sloped(m: &mut Model, k: f64, (y0, y1): (f64, f64)) -> Shape {
    let at_x0 = |y: f64, z: f64| Point::new(0.0, y, z);
    let at_x1 = |y: f64, z: f64| Point::new(k.mul_add(z, 20.0), y, z);
    let slope = (y1 - y0) / 20.0;
    let foot = (slope * k).mul_add(-6.0, y1);
    let near = [
        at_x0(0.0, 0.0),
        at_x0(0.0, -15.0),
        at_x0(12.0, -15.0),
        at_x0(12.0, -6.0),
        at_x0(y0, -6.0),
        at_x0(y0, 0.0),
    ];
    let far = [
        at_x1(0.0, 0.0),
        at_x1(0.0, -15.0),
        at_x1(12.0, -15.0),
        at_x1(12.0, -6.0),
        at_x1(foot, -6.0),
        at_x1(y1, 0.0),
    ];
    let n = near.len();
    let mut points = near.to_vec();
    points.extend_from_slice(&far);
    let mut rings = vec![(0..n).collect::<Vec<_>>(), (n..2 * n).rev().collect()];
    rings.extend((0..n).map(|i| vec![i, n + i, n + (i + 1) % n, (i + 1) % n]));
    faceted(m, &points, &rings)
}

/// The same request with hold edge on and off: the same solid's volume, or
/// the same refusal; and with trim through off, refused.
fn as_with_hold_edge_off(m: &mut Model, s: &Shape) {
    let edge = front_edge(m, s);
    let mut off = Overflow::default();
    off.hold_edge = false;
    let on = fillet_edges(m, s, &edge, R, T).map(|b| volume(m, &b.shape));
    let without = fillet_edges_overflow(m, s, &edge, R, off, T).map(|b| volume(m, &b.shape));
    eprintln!("on {on:?}, off {without:?}");
    match (on, without) {
        (Ok(a), Ok(b)) => assert!((a - b).abs() < 1e-9, "{a} against {b}"),
        (Err(a), Err(b)) => assert_eq!(a.to_string(), b.to_string()),
        (a, b) => panic!("hold edge on gives {a:?}, off {b:?}"),
    }
    let mut neither = Overflow::default();
    neither.trim_through = false;
    let refused = fillet_edges_overflow(m, s, &edge, R, neither, T);
    assert!(refused.is_err(), "with trim through off");
}

/// Feet that never reach the end face: the top's edge bends at `(18, 2.5)`
/// to fall steeply to 0.4 behind the front at the end. Past the end face
/// the ball rests on that edge carried on, its feet trailing it, and the
/// carried edge falls away from the band faster than a ball can follow: no
/// ball touches the front and rests on it before the feet reach the end
/// face. Not held; the request goes on to trim through.
#[test]
fn feet_that_never_reach_the_end_face_go_on_to_trim_through() {
    let mut m = Model::new();
    let s = kinked(&mut m, [(0.0, 4.6), (18.0, 2.5), (20.0, 0.4)]);
    as_with_hold_edge_off(&mut m, &s);
}

/// A held edge leaving a leaning end face steeply: the top narrows from 9
/// behind the front to 0.5 at the end face `x = 20 - z / 2`, which leans
/// out at its foot. The edge carried on past the face crosses the round's
/// own edge at `x = 21.2`, where the balls past the face, still reaching
/// inside it under the lean, would jump back onto the band. Not held; the
/// request goes on to trim through.
#[test]
fn a_held_edge_crossing_the_round_s_edge_past_the_end_goes_on_to_trim_through() {
    let mut m = Model::new();
    let s = sloped(&mut m, -0.5, (9.0, 0.5));
    as_with_hold_edge_off(&mut m, &s);
}

/// A vertex piece cut by the end face: the top's edge bends at
/// `(19.3, 1.8)`, so the ball rests on that vertex from `x = 19.62` to past
/// the end, then on the next edge, whose feet trail the balls back inside
/// the end face. The round runs on past the face, but the boolean does not
/// resolve the vertex piece's cut by it, so the request goes on to trim
/// through.
#[test]
fn a_vertex_piece_cut_by_the_end_face_goes_on_to_trim_through() {
    let mut m = Model::new();
    let s = kinked(&mut m, [(0.0, 4.6), (19.3, 1.8), (20.0, 1.45)]);
    as_with_hold_edge_off(&mut m, &s);
}
