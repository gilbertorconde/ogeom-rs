//! Hold edge where the held edge lies on a drum: the face the round's
//! contact runs past an edge of is a cylinder whose axis runs along the
//! round's edge, the solid inside it, and the other face a plane. The
//! reference is the definition: the ball touches the plane and rests on the
//! nearest point of the drum's region, here the drum less a bore's mouth.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{PI, TAU};

use ogeom_core::Tolerances;
use ogeom_fillet::{BlendRequest, Overflow};
use ogeom_geom::{Curve3d as _, Surface as _, SurfaceGeometry};
use ogeom_math::{Direction, Frame, Point, Vector};
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The round's radius.
const R: f64 = 4.0;

/// The drum's axis `(y, z)`, along `x`, and its radius: through the edge
/// `y = z = 0`, the solid inside it.
const DRUM: (f64, f64, f64) = (8.0, -6.0, 10.0);

/// The bore into the drum's front, along `y`: its axis `(x, z)` and radius.
const BORE: (f64, f64, f64) = (15.0, -4.5, 2.0);

/// The part's length along the edge.
const LEN: f64 = 30.0;

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

fn volume(m: &Model, s: &Shape) -> f64 {
    ogeom_algo::volume_properties(m, s, ogeom_mesh::Deflection::default(), T)
        .unwrap()
        .mass
}

/// The front top edge `y = z = 0`, the one edge asked.
fn front_edge(m: &Model, s: &Shape) -> Vec<Shape> {
    let edge = edges_where(m, s, |p| p.y.abs() < 1e-9 && p.z.abs() < 1e-9);
    assert_eq!(edge.len(), 1, "the front top edge");
    edge
}

/// A part `LEN` long whose top `z = 0` runs from the front top edge to a
/// back face at `y = 14`, its front the drum [`DRUM`] bulging out below the
/// edge, the bottom at `z = -12`, its end faces leaning to `x = k z` and
/// `x = LEN + k z`; a bore of radius `BORE.2` along `y` into the front about
/// `(x, z) = (bx, BORE.1)`, where given.
fn part(m: &mut Model, bore: Option<f64>, k: f64) -> Shape {
    let block = if k == 0.0 {
        ogeom_algo::make_box(m, at(0.0, -5.0, -12.0), (LEN, 19.0, 12.0), T)
            .unwrap()
            .shape
    } else {
        // The block's side `y = -5` as a quadrilateral, swept along `y`.
        let corners = [(0.0, -12.0), (LEN, -12.0), (LEN, 0.0), (0.0, 0.0)]
            .map(|(x, z): (f64, f64)| Point::new(x + k * z, -5.0, z));
        let vs: Vec<Shape> = corners
            .iter()
            .map(|p| ogeom_algo::make_vertex(m, *p).shape)
            .collect();
        let es: Vec<Shape> = (0..4)
            .map(|i| {
                let (a, b) = (corners[i], corners[(i + 1) % 4]);
                let line =
                    ogeom_geom::Curve::Line(ogeom_geom::LineCurve::segment(a, b, T).unwrap());
                ogeom_algo::make_edge_between(
                    m,
                    line,
                    (0.0, a.distance(b)),
                    &vs[i],
                    &vs[(i + 1) % 4],
                    T,
                )
                .unwrap()
                .shape
            })
            .collect();
        let wire = ogeom_algo::make_wire(m, &es, T).unwrap().shape;
        let frame = Frame::new(Point::new(0.0, -5.0, 0.0), Direction::Y, Direction::X, T).unwrap();
        let plane = ogeom_geom::PlaneSurface::new(ogeom_math::Plane::new(frame));
        let face = ogeom_algo::make_face(m, SurfaceGeometry::Plane(plane), &[wire], T)
            .unwrap()
            .shape;
        ogeom_algo::make_prism(m, &face, Vector::new(0.0, 19.0, 0.0), T)
            .unwrap()
            .shape
    };
    let (ay, az, ar) = DRUM;
    let axis = Frame::new(Point::new(-10.0, ay, az), Direction::X, Direction::Y, T).unwrap();
    let drum = ogeom_algo::make_cylinder(m, axis, ar, LEN + 20.0, T)
        .unwrap()
        .shape;
    let s = ogeom_bool::common(m, &block, &drum, T).unwrap().shape;
    let Some(bx) = bore else {
        return s;
    };
    let (_, bz, br) = BORE;
    let frame = Frame::new(Point::new(bx, -5.0, bz), Direction::Y, Direction::X, T).unwrap();
    let cutter = ogeom_algo::make_cylinder(m, frame, br, 10.0, T)
        .unwrap()
        .shape;
    ogeom_bool::cut(m, &s, &cutter, T).unwrap().shape
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

/// `f` integrated over `[a, b]` by Gauss-Legendre on `g`.
fn integral(g: &[(f64, f64)], a: f64, b: f64, f: impl Fn(f64) -> f64) -> f64 {
    let h = 0.5 * (b - a);
    g.iter().map(|(x, w)| w * h * f(a + h * (x + 1.0))).sum()
}

/// The drum's upper surface over `y` (in front of its axis, `y < 8`, the
/// front's arc from the bottom of the part to the top).
fn drum_top(y: f64) -> f64 {
    let (ay, az, ar) = DRUM;
    az + (ar * ar - (y - ay) * (y - ay)).sqrt()
}

/// The bore's rim on the drum's front, by the angle `phi` round the bore's
/// axis from `+x` toward `+z`, and its first two derivatives.
fn rim(phi: f64) -> [Vector; 3] {
    let (ay, az, ar) = DRUM;
    let (bx, bz, br) = BORE;
    let (s, c) = phi.sin_cos();
    let (z, z1, z2) = (bz + br * s, br * c, -br * s);
    let w = z - az;
    let root = (ar * ar - w * w).sqrt();
    let root1 = -w * z1 / root;
    let root2 = -(z1 * z1 + w * z2) / root + w * z1 * root1 / (root * root);
    [
        Vector::new(bx + br * c, ay - root, z),
        Vector::new(-br * s, -root1, z1),
        Vector::new(-br * c, -root2, z2),
    ]
}

/// The point of the rim nearest `p`: the nearest of 720 points along it,
/// then Newton's step on the distance's derivative.
fn rim_nearest(p: Point) -> Point {
    let q = p - Point::ORIGIN;
    let mut phi = (0..720)
        .map(|k| TAU * f64::from(k) / 720.0)
        .min_by(|a, b| {
            (rim(*a)[0] - q)
                .magnitude()
                .total_cmp(&(rim(*b)[0] - q).magnitude())
        })
        .unwrap();
    for _ in 0..50 {
        let [c, d1, d2] = rim(phi);
        let step = (c - q).dot(d1) / (d1.dot(d1) + (c - q).dot(d2));
        phi -= step;
        if step.abs() < 1e-16 {
            break;
        }
    }
    Point::ORIGIN + rim(phi)[0]
}

/// Whether a point of the drum's front lies in the bore's mouth.
fn in_bore(p: Point) -> bool {
    let (bx, bz, br) = BORE;
    (p.x - bx).hypot(p.z - bz) < br
}

/// A ball of the reference: its centre, and its foot on the drum or the rim.
#[derive(Clone, Copy)]
struct Ball {
    centre: Point,
    foot: Point,
}

/// The ball at `x`, by the definition: its centre `R` below the top, and as
/// far from the drum's axis as it may be while it stands `R` from the
/// drum's region; on the band, `R` inside the drum, its foot the drum's
/// point along the radius; where that point lies in the bore's mouth, moved
/// toward the drum along the top's offset, by bisection, until the rim's
/// nearest point stands `R` off. `held` asks for the second.
fn ball(x: f64, held: bool) -> Ball {
    let (ay, az, ar) = DRUM;
    let band = ay - ((ar - R) * (ar - R) - (az + R) * (az + R)).sqrt();
    let centre = Point::new(x, band, -R);
    let radial = (centre - Point::new(x, ay, az)) * (ar / (ar - R));
    let foot = Point::new(x, ay, az) + radial;
    if !held {
        assert!(!in_bore(foot), "the band's contact at {x} lies in the bore");
        return Ball { centre, foot };
    }
    let gap = |y: f64| {
        let c = Point::new(x, y, -R);
        c.distance(rim_nearest(c)) - R
    };
    // From the band's centre to where the top's offset meets the drum.
    let (mut lo, mut hi) = (ay - (ar * ar - (az + R) * (az + R)).sqrt(), band);
    assert!(
        gap(lo) < 0.0 && gap(hi) > 0.0,
        "the rim brackets the ball at {x}"
    );
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if gap(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let centre = Point::new(x, 0.5 * (lo + hi), -R);
    Ball {
        centre,
        foot: rim_nearest(centre),
    }
}

/// The centres' path, `y(x)` on the top's offset, its slope (the distance
/// to the seat's gradient is `centre - foot`, square to the path) and its
/// curvature's numerator `y''`, by central differences of the slope.
fn path(x: f64, held: bool) -> (Ball, f64, f64) {
    let slope = |x: f64| {
        let b = ball(x, held);
        let g = b.centre - b.foot;
        -g.x / g.y
    };
    let h = 1e-5;
    (
        ball(x, held),
        slope(x),
        (slope(x + h) - slope(x - h)) / (2.0 * h),
    )
}

/// The section's area at `x`, cut away: over the top's offset, along the
/// path's normal toward the edge, from the path out to the foot (`sigma`
/// from 0 to `s_f`), the solid between the balls' envelope
/// `z = -R + sqrt(R^2 - sigma^2)` and the top of the solid (the top `z = 0`
/// over `y >= 0`, the drum's front over `y < 0`), less the bore, times the
/// band's stretch `1 + sigma y'' / L^3`, `L = sqrt(1 + y'^2)`. Every point
/// over the offset nearer one centre than the radius stands inside that
/// ball, and one farther lies between the arcs of its nearest ball to the
/// top and to the foot exactly where its `sigma` stays within the foot's.
fn section(x: f64, held: bool, g: &[(f64, f64)]) -> f64 {
    let (b, y1, y2) = path(x, held);
    let len = y1.hypot(1.0);
    let n = Vector::new(y1, -1.0, 0.0) / len;
    let s_f = (b.foot - b.centre).dot(n);
    let (bx, bz, br) = BORE;
    let height = |s: f64| {
        let q = b.centre + n * s;
        let stretch = 1.0 + s * y2 / (len * len * len);
        assert!(stretch > 0.0, "the band folds at {x}, {s}");
        let env = -R + (R * R - s * s).sqrt();
        let top = if q.y >= 0.0 { 0.0 } else { drum_top(q.y) };
        let mut h = (top - env).max(0.0);
        let d = q.x - bx;
        if d.abs() < br {
            let half = (br * br - d * d).sqrt();
            let (lo, hi) = (env.max(bz - half), top.min(bz + half));
            h -= (hi - lo).max(0.0);
        }
        h * stretch
    };
    // Split where the top of the solid turns from the top to the drum.
    let s0 = -b.centre.y / n.y;
    len * (integral(g, 0.0, s0, height) + integral(g, s0, s_f, height))
}

/// The seat changes: where the band's contact on the drum, at height
/// `az + (ar / (ar - R)) (R + az)`, crosses the bore's rim.
fn changes() -> (f64, f64) {
    let (_, az, ar) = DRUM;
    let (bx, bz, br) = BORE;
    let z = az + ar / (ar - R) * (-R - az);
    let half = (br * br - (z - bz) * (z - bz)).sqrt();
    (bx - half, bx + half)
}

/// The band's section in closed form: the triangle from the edge to the
/// ball's two contacts, less the ball's segment over the chord between
/// them, plus the drum's segment over the chord from the edge to its
/// contact.
fn band_section() -> f64 {
    let b = ball(0.0, false);
    let (a, top) = (b.foot, Point::new(0.0, b.centre.y, 0.0));
    let triangle = 0.5 * (a.y * top.z - top.y * a.z).abs();
    let segment = |radius: f64, chord: f64| {
        let theta = 2.0 * (0.5 * chord / radius).asin();
        0.5 * radius * radius * (theta - theta.sin())
    };
    triangle - segment(R, a.distance(top)) + segment(DRUM.2, a.distance(Point::ORIGIN))
}

/// The volume cut away by the reference with `n` points per range: the
/// band's closed-form section over the band, the held stretch's sections
/// integrated in `x`, split at the bore's middle.
fn reference(n: usize) -> f64 {
    let g = gauss(n);
    let (c, d) = changes();
    let mid = BORE.0;
    let held = integral(&g, c, mid, |x| section(x, true, &g))
        + integral(&g, mid, d, |x| section(x, true, &g));
    band_section() * (LEN - (d - c)) + held
}

/// The largest angle between the normals either side of each edge a fitted
/// face shares with the band's cylinder, the drum or the top, sampled along
/// the edge, and how many such edges.
fn worst_join(m: &Model, s: &Shape) -> (f64, usize) {
    let surface_of = |face: &Shape| -> SurfaceGeometry {
        let data = m.node(face).unwrap().data().as_face().unwrap();
        m.geometry().surface(data.surface).unwrap().clone()
    };
    let mut worst = 0.0_f64;
    let mut joins = 0;
    for e in explore_unique(m, s, ShapeType::Edge).unwrap() {
        let faces = ancestors_of(m, s, &e, ShapeType::Face).unwrap();
        let [a, b] = &faces[..] else {
            continue;
        };
        let (sa, sb) = (surface_of(a), surface_of(b));
        let other = match (&sa, &sb) {
            (SurfaceGeometry::BSpline(_), o) | (o, SurfaceGeometry::BSpline(_)) => o,
            _ => continue,
        };
        let Some((c, (lo, hi))) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap() else {
            continue;
        };
        let mid = c.point_at(0.5 * (lo + hi), T).unwrap();
        let tangent_kind = match other {
            SurfaceGeometry::Cylinder(cy) => {
                cy.cylinder()
                    .axis()
                    .direction
                    .vector()
                    .cross(Vector::X)
                    .magnitude()
                    < 1e-9
            }
            SurfaceGeometry::Plane(_) => mid.z.abs() < 1e-6,
            _ => false,
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

/// The result's edges on the rim's upper half within 1e-5 of the exact
/// rim: how far, at most, they stray from it, and the x-range they cover.
fn kept(m: &Model, s: &Shape) -> (f64, (f64, f64)) {
    let off = |p: Point| p.distance(rim_nearest(p));
    let mut worst = 0.0_f64;
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for e in edges_where(m, s, |p| off(p) < 1e-5 && p.z > BORE.1) {
        let (c, (ta, tb)) = ogeom_algo::placed_edge_curve(m, &e, T).unwrap().unwrap();
        for k in 0..=50 {
            let p = c.point_at(ta + (tb - ta) * f64::from(k) / 50.0, T).unwrap();
            worst = worst.max(off(p));
            lo = lo.min(p.x);
            hi = hi.max(p.x);
        }
    }
    (worst, (lo, hi))
}

/// The ball rests on a bore's rim on the drum: the front top edge of the
/// drum part rounded at 4 between the top and the drum, a bore of radius 2
/// along `y` into the drum about `(x, z) = (15, -4.5)`. The band's contact
/// on the drum stands at `z = -2.667`, 0.167 under the bore's top, so
/// between `x = 14.2007` and `15.7993` it lies in the bore's mouth and the
/// ball rests on the rim. The rim is kept, the fitted piece meets the band,
/// the top and the drum G1, and the volume cut away is the reference's.
#[test]
fn a_rim_on_the_drum_cuts_away_the_reference() {
    let mut m = Model::new();
    let s = part(&mut m, Some(BORE.0), 0.0);
    let edge = front_edge(&m, &s);
    let before = volume(&m, &s);
    let built = BlendRequest::new()
        .round(edge.clone(), R)
        .build(&mut m, &s, T)
        .unwrap();

    let found = ogeom_algo::check(&m, &built.shape, T).unwrap();
    assert!(found.is_valid(), "{:?}", found.problems);
    for shell in explore_unique(&m, &built.shape, ShapeType::Shell).unwrap() {
        assert!(ogeom_algo::is_shell_closed(&m, &shell).unwrap(), "closed");
    }
    assert!(!built.history.generated(&edge[0]).is_empty(), "the round");

    let (crease, joins) = worst_join(&m, &built.shape);
    eprintln!("worst join {crease:.1e} rad over {joins} edges");
    assert!(joins >= 3, "the fitted piece's joins");
    assert!(crease < 1e-5, "the fitted piece's joins crease by {crease}");

    // The feet run over the rim's top, from the ball at one seat change to
    // the ball at the other.
    let (c, d) = changes();
    let (f0, f1) = (ball(c + 1e-9, true).foot, ball(d - 1e-9, true).foot);
    let (worst, (lo, hi)) = kept(&m, &built.shape);
    eprintln!("rim kept over x {lo:.6} .. {hi:.6}, off the exact rim by {worst:.1e}");
    assert!(worst < 1e-5, "the rim strays {worst} off its curve");
    assert!(
        lo <= f0.x && hi >= f1.x,
        "the rim kept over {lo}..{hi}, the feet over {}..{}",
        f0.x,
        f1.x
    );

    let got = before - volume(&m, &built.shape);
    let (coarse, fine) = (reference(24), reference(48));
    let bound = (fine - coarse).abs();
    eprintln!(
        "held over {c:.9} .. {d:.9}; cut away {got:.9}, reference {fine:.9}, bound {bound:.1e}, error {:.1e}",
        (got - fine).abs()
    );
    assert!(
        (got - fine).abs() < 1e-6 + bound,
        "cut away {got} against the reference {fine} (within {bound})"
    );
}

/// The reference reads the band's closed form through its own sections,
/// and its held stretch is symmetric about the bore's middle, as the part
/// is.
#[test]
fn the_drum_rim_reference_reads_the_band() {
    let g = gauss(48);
    let sampled = section(3.0, false, &g);
    let closed = band_section();
    eprintln!("band section {closed:.12}, sampled {sampled:.12}");
    assert!(
        (closed - sampled).abs() < 1e-12,
        "{closed} against {sampled}"
    );
    let (c, d) = changes();
    let mid = BORE.0;
    let left = integral(&g, c, mid, |x| section(x, true, &g));
    let right = integral(&g, mid, d, |x| section(x, true, &g));
    eprintln!("held halves {left:.12}, {right:.12}");
    assert!((left - right).abs() < 1e-10, "{left} against {right}");
}

/// The same request with hold edge on and off: the same solid's volume, or
/// the same refusal.
fn as_with_hold_edge_off(m: &mut Model, s: &Shape, edge: &[Shape]) {
    let mut off = Overflow::default();
    off.hold_edge = false;
    let request = BlendRequest::new().round(edge.iter().cloned(), R);
    let on = request.build(m, s, T).map(|b| volume(m, &b.shape));
    let without = request
        .clone()
        .overflow(off)
        .build(m, s, T)
        .map(|b| volume(m, &b.shape));
    eprintln!("on {on:?}, off {without:?}");
    match (on, without) {
        (Ok(a), Ok(b)) => assert!((a - b).abs() < 1e-9, "{a} against {b}"),
        (Err(a), Err(b)) => assert_eq!(a.to_string(), b.to_string()),
        (a, b) => panic!("hold edge on gives {a:?}, off {b:?}"),
    }
}

/// With hold edge and trim through off, the rim on the drum is refused
/// naming hold edge.
#[test]
fn a_rim_on_the_drum_with_hold_edge_off_is_refused_by_name() {
    let mut m = Model::new();
    let s = part(&mut m, Some(BORE.0), 0.0);
    let edge = front_edge(&m, &s);
    let mut off = Overflow::default();
    off.hold_edge = false;
    off.trim_through = false;
    let err = BlendRequest::new()
        .round(edge, R)
        .overflow(off)
        .build(&mut m, &s, T)
        .unwrap_err();
    assert!(err.to_string().contains("hold edge is off"), "{err}");
}

/// A rim on the drum the balls rest on up to an end face is not held: the
/// request goes on to trim through, as with hold edge off.
#[test]
fn a_rim_on_the_drum_reaching_an_end_goes_on_to_trim_through() {
    let mut m = Model::new();
    let s = part(&mut m, Some(0.5), 0.0);
    let edge = front_edge(&m, &s);
    as_with_hold_edge_off(&mut m, &s, &edge);
}

/// A drum part whose end faces lean is not held on the drum: a leaning end
/// face cuts the drum in an ellipse. The request goes on to trim through.
#[test]
fn a_rim_on_the_drum_between_leaning_ends_goes_on_to_trim_through() {
    let mut m = Model::new();
    let s = part(&mut m, Some(BORE.0), 0.3);
    let edge = front_edge(&m, &s);
    as_with_hold_edge_off(&mut m, &s, &edge);
}

/// A slot across the drum's front whose sides stand square to its axis cuts
/// the drum in circles across it, where the ball would rest on their
/// vertices with the slot's floor: not held, the request goes on to trim
/// through.
#[test]
fn a_slot_across_the_drum_goes_on_to_trim_through() {
    let mut m = Model::new();
    let s = part(&mut m, None, 0.0);
    let slot = ogeom_algo::make_box(&mut m, at(13.0, -5.0, -9.0), (4.0, 4.0, 6.5), T)
        .unwrap()
        .shape;
    let s = ogeom_bool::cut(&mut m, &s, &slot, T).unwrap().shape;
    let edge = front_edge(&m, &s);
    as_with_hold_edge_off(&mut m, &s, &edge);
}
