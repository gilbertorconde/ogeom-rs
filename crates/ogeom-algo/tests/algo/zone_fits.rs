//! Minimum-zone, inscribed and circumscribed fits against closed forms and
//! against independent searches: exact shapes give zones of width zero,
//! constructed sets give their known widths, and on noisy sets no
//! perturbation of the result's parameters, drawn at every scale from a
//! tenth of the part down to a nanometre, narrows the zone.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::cast_precision_loss,
    reason = "test code"
)]

use ogeom_algo::fit::{self, FitOptions, zone};
use ogeom_core::Tolerances;
use ogeom_math::{Axis, Direction, Frame, Point, Point2, Vector, Vector2};

const T: Tolerances = Tolerances::millimetres();

/// A fixed sequence of draws (xorshift), so every run sees the same points.
struct Seeded(u64);

impl Seeded {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// Uniform in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1_u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }

    /// A standard normal draw, by Box-Muller.
    fn gaussian(&mut self) -> f64 {
        let u = 1.0 - self.unit();
        let v = self.unit();
        (-2.0 * u.ln()).sqrt() * (core::f64::consts::TAU * v).cos()
    }

    /// A uniform direction.
    fn direction(&mut self) -> Vector {
        loop {
            let v = Vector::new(self.gaussian(), self.gaussian(), self.gaussian());
            let m = v.magnitude();
            if m > 1e-3 {
                return v / m;
            }
        }
    }
}

/// A frame tilted off the global axes, so no fit is handed a coordinate
/// plane.
fn tilted() -> Frame {
    let z = Direction::new(Vector::new(0.3, -0.5, 1.0), T).unwrap();
    Frame::about(Point::new(12.0, -7.0, 3.5), z)
}

fn place(frame: Frame, u: f64, v: f64, w: f64) -> Point {
    frame.origin() + frame.x().vector() * u + frame.y().vector() * v + frame.z().vector() * w
}

/// Parallel to rounding in a computed direction.
fn parallel(a: Direction, b: Direction) -> bool {
    a.dot(b).abs() > 1.0 - 1e-12
}

/// Scales from a tenth of a part a few tens of millimetres across down to
/// a nanometre.
const SCALES: [f64; 8] = [1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7];

/// Rounding in a width a few tens of millimetres across.
const ROUNDING: f64 = 1e-12;

// ---- flatness ----

fn slab_width(points: &[Point], n: Vector) -> f64 {
    let (lo, hi) = points
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
            let h = n.dot(p.to_vector());
            (lo.min(h), hi.max(h))
        });
    hi - lo
}

/// The issue's case: points on a plane, four of them at ±w/2 across it,
/// and a grid on each of two parallel planes: the width is the distance.
#[test]
fn points_on_two_parallel_planes_are_their_distance_flat() {
    let frame = tilted();
    let w = 0.037;
    let mut plate = Vec::new();
    for i in 0..9 {
        for j in 0..7 {
            plate.push(place(
                frame,
                f64::from(i) * 3.0 - 12.0,
                f64::from(j) * 2.5 - 7.0,
                0.0,
            ));
        }
    }
    let mut four = plate.clone();
    for (u, v, s) in [
        (-10.0, -6.0, 1.0),
        (11.0, -5.0, -1.0),
        (9.0, 6.0, 1.0),
        (-11.0, 5.0, -1.0),
    ] {
        four.push(place(frame, u, v, s * w / 2.0));
    }
    let zone = zone::flatness(&four, T).unwrap();
    assert!((zone.width - w).abs() < 1e-9, "{} against {w}", zone.width);
    assert!(zone.optimal);
    assert_eq!(zone.contacts.len(), 4, "{:?}", zone.contacts);

    let mut two = plate.clone();
    two.extend(plate.iter().map(|p| *p + frame.z().vector() * w));
    let zone = zone::flatness(&two, T).unwrap();
    assert!((zone.width - w).abs() < 1e-9, "{} against {w}", zone.width);
    assert!(parallel(zone.shape.normal(), frame.z()));
    assert_eq!(zone.contacts.len(), two.len());
}

/// Exact integer grids, every row collinear and every layer coplanar:
/// two layers a unit apart are a unit flat, and one layer is flat.
#[test]
fn integer_grids_are_exactly_flat() {
    let mut two = Vec::new();
    for i in 0..20 {
        for j in 0..20 {
            two.push(Point::new(f64::from(i), f64::from(j), 0.0));
            two.push(Point::new(f64::from(i), f64::from(j), 1.0));
        }
    }
    let zone = zone::flatness(&two, T).unwrap();
    assert_eq!(zone.width, 1.0);
    assert_eq!(zone.contacts.len(), two.len());
    let one: Vec<Point> = two.iter().filter(|p| p.z == 0.0).copied().collect();
    assert_eq!(zone::flatness(&one, T).unwrap().width, 0.0);
}

#[test]
fn an_exact_plane_is_flat() {
    let frame = tilted();
    let mut draw = Seeded(7);
    let points: Vec<Point> = (0..500)
        .map(|_| place(frame, draw.range(-20.0, 20.0), draw.range(-10.0, 10.0), 0.0))
        .collect();
    let zone = zone::flatness(&points, T).unwrap();
    assert!(zone.width < 1e-12, "{}", zone.width);
}

/// On a noisy plate the zone is no wider than the least-squares peak to
/// peak, and no direction near or far, from every direction of a dense
/// sphere and random tilts of the result at every scale, is narrower.
#[test]
fn a_noisy_plate_is_least_flat() {
    let frame = tilted();
    let mut draw = Seeded(11);
    let points: Vec<Point> = (0..2000)
        .map(|_| {
            let (u, v) = (draw.range(-20.0, 20.0), draw.range(-10.0, 10.0));
            // A gentle warp under the noise, so the zone is not the noise's.
            let w = 0.01 * draw.gaussian() + 2e-4 * u * v / 10.0;
            place(frame, u, v, w)
        })
        .collect();
    let zone = zone::flatness(&points, T).unwrap();
    let ls = fit::plane(&points, &FitOptions::default(), T).unwrap();
    assert!(
        zone.width <= ls.form.1 - ls.form.0,
        "{} against {:?}",
        zone.width,
        ls.form
    );
    let n = zone.shape.normal().vector();
    assert!((slab_width(&points, n) - zone.width).abs() < ROUNDING);
    for scale in SCALES {
        for _ in 0..500 {
            let m = (n + draw.direction() * scale).normalized(T).unwrap();
            assert!(
                slab_width(&points, m) >= zone.width - ROUNDING,
                "a tilt of {scale} narrows {} to {}",
                zone.width,
                slab_width(&points, m)
            );
        }
    }
}

/// On a dozen points the zone is the least over every slab a face or two
/// edges of the set fix, tried directly: points drawn at random, and
/// points on a small integer lattice, where rows are collinear and layers
/// coplanar and the hull has faces in one plane.
#[test]
fn a_dozen_points_are_least_flat_by_brute_force() {
    let mut draw = Seeded(3);
    for round in 0..40 {
        let lattice = round % 2 == 1;
        let points: Vec<Point> = (0..12)
            .map(|_| {
                if lattice {
                    let mut pick = || (draw.unit() * 4.0).floor();
                    Point::new(pick(), pick(), pick())
                } else {
                    Point::new(
                        draw.range(-5.0, 5.0),
                        draw.range(-5.0, 5.0),
                        draw.range(-1.0, 1.0),
                    )
                }
            })
            .collect();
        let zone = zone::flatness(&points, T).unwrap_or_else(|e| panic!("round {round}: {e}"));
        let mut best = f64::INFINITY;
        let mut tryn = |n: Vector| {
            if let Ok(n) = n.normalized(T) {
                best = best.min(slab_width(&points, n));
            }
        };
        for a in 0..12 {
            for b in 0..12 {
                for c in 0..12 {
                    tryn((points[b] - points[a]).cross(points[c] - points[a]));
                    for d in 0..12 {
                        tryn((points[b] - points[a]).cross(points[d] - points[c]));
                    }
                }
            }
        }
        assert!(
            (zone.width - best).abs() < ROUNDING,
            "round {round}: {} against {best}",
            zone.width
        );
    }
}

// ---- circles ----

fn ring_width(points: &[Point2], c: Point2) -> f64 {
    let (lo, hi) = points
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
            let r = p.distance(c);
            (lo.min(r), hi.max(r))
        });
    hi - lo
}

fn on_circle(c: Point2, r: f64, a: f64) -> Point2 {
    c + Vector2::new(a.cos(), a.sin()) * r
}

#[test]
fn an_exact_circle_is_round() {
    let c = Point2::new(3.0, -2.0);
    let points: Vec<Point2> = (0..360)
        .map(|i| on_circle(c, 25.0, f64::from(i) * 0.0174))
        .collect();
    let zone = zone::circularity(&points, T).unwrap();
    assert!(zone.width < 1e-12, "{}", zone.width);
    assert!(zone.shape.centre().distance(c) < 1e-9);
}

/// A three-lobed profile of known lobes: the zone is no wider than the
/// least-squares peak to peak, and no centre moved at any scale narrows it;
/// on a grid of centres a tenth of a micrometre apart, none does either.
#[test]
fn a_lobed_profile_is_least_round() {
    let c = Point2::new(3.0, -2.0);
    let mut draw = Seeded(5);
    let points: Vec<Point2> = (0..1500)
        .map(|i| {
            let a = f64::from(i) * core::f64::consts::TAU / 1500.0;
            let r = 25.0 + 0.02 * (3.0 * a).cos() + 0.004 * draw.gaussian();
            on_circle(c, r, a)
        })
        .collect();
    let zone = zone::circularity(&points, T).unwrap();
    let ls = fit::circle2(&points, &FitOptions::default(), T).unwrap();
    assert!(
        zone.width <= ls.form.1 - ls.form.0,
        "{} against {:?}",
        zone.width,
        ls.form
    );
    assert!(zone.contacts.len() >= 4, "{:?}", zone.contacts);
    let at = zone.shape.centre();
    for scale in SCALES {
        for _ in 0..500 {
            let moved = at + Vector2::new(draw.gaussian(), draw.gaussian()) * scale;
            assert!(
                ring_width(&points, moved) >= zone.width - ROUNDING,
                "a move of {scale}"
            );
        }
    }
    let step = 1e-7;
    for i in -20..=20 {
        for j in -20..=20 {
            let moved = at + Vector2::new(f64::from(i), f64::from(j)) * step;
            assert!(ring_width(&points, moved) >= zone.width - ROUNDING);
        }
    }
}

/// On ten points the zone is the least over a grid of centres a micrometre
/// apart about it, and within the width's Lipschitz bound of the grid's
/// least: a ring's width moves at most twice as fast as its centre.
#[test]
fn ten_points_are_least_round_on_a_grid() {
    let mut draw = Seeded(17);
    for round in 0..10 {
        let points: Vec<Point2> = (0..10)
            .map(|i| {
                on_circle(
                    Point2::new(0.0, 0.0),
                    10.0 + draw.range(-0.5, 0.5),
                    f64::from(i) * 0.6,
                )
            })
            .collect();
        let zone = zone::circularity(&points, T).unwrap();
        let at = zone.shape.centre();
        let step = 2e-3;
        let mut least = f64::INFINITY;
        for i in -400..=400 {
            for j in -400..=400 {
                let c = at + Vector2::new(f64::from(i), f64::from(j)) * step;
                least = least.min(ring_width(&points, c));
            }
        }
        assert!(
            least >= zone.width - ROUNDING,
            "round {round}: {least} below {}",
            zone.width
        );
        assert!(least <= zone.width + 2.0 * step, "round {round}");
    }
}

/// Points alternating either side of a line fit between two lines closer
/// than any ring: there is no narrowest circle, and the fit says so rather
/// than report a huge one.
#[test]
fn a_zigzag_has_no_narrowest_circle() {
    let zigzag: Vec<Point2> = (0..200)
        .map(|i| Point2::new(f64::from(i) * 0.05, if i % 2 == 0 { 0.001 } else { -0.001 }))
        .collect();
    let refused = zone::circularity(&zigzag, T).unwrap_err();
    assert!(
        refused.to_string().contains("narrower than any ring"),
        "{refused}"
    );
}

/// An acute triangle's circle is its circumcircle; an obtuse one's is the
/// circle on its longest side.
#[test]
fn a_triangles_circumscribed_circle_is_the_closed_form() {
    let (a, b, c) = (
        Point2::new(0.0, 0.0),
        Point2::new(4.0, 0.0),
        Point2::new(1.0, 3.0),
    );
    let zone = zone::circumscribed_circle(&[a, b, c], T).unwrap();
    // Equidistant from all three: x = 2, and (1 − 2)² + (3 − y)² = 4 + y².
    let centre = Point2::new(2.0, 1.0);
    assert!(zone.shape.centre().distance(centre) < 1e-12);
    assert!((zone.shape.radius() - 5.0_f64.sqrt()).abs() < 1e-12);
    assert_eq!(zone.contacts, vec![0, 1, 2]);

    let obtuse = [
        Point2::new(0.0, 0.0),
        Point2::new(6.0, 0.0),
        Point2::new(2.0, 1.0),
    ];
    let zone = zone::circumscribed_circle(&obtuse, T).unwrap();
    assert!(zone.shape.centre().distance(Point2::new(3.0, 0.0)) < 1e-12);
    assert!((zone.shape.radius() - 3.0).abs() < 1e-12);
    assert_eq!(zone.contacts, vec![0, 1]);
}

/// A regular polygon's corners: the circumscribed circle is through them
/// and the inscribed circle too, both centred on the polygon.
#[test]
fn a_regular_polygons_circles_are_exact() {
    let c = Point2::new(-4.0, 9.0);
    for sides in [3, 5, 6, 12] {
        let corners: Vec<Point2> = (0..sides)
            .map(|k| {
                on_circle(
                    c,
                    7.0,
                    f64::from(k) * core::f64::consts::TAU / f64::from(sides) + 0.2,
                )
            })
            .collect();
        let outer = zone::circumscribed_circle(&corners, T).unwrap();
        assert!(outer.shape.centre().distance(c) < 1e-12, "{sides} sides");
        assert!((outer.shape.radius() - 7.0).abs() < 1e-12, "{sides} sides");
        let inner = zone::inscribed_circle(&corners, T).unwrap();
        assert!(inner.shape.centre().distance(c) < 1e-9, "{sides} sides");
        assert!((inner.shape.radius() - 7.0).abs() < 1e-9, "{sides} sides");

        // The sides sampled too: the inscribed circle is the incircle,
        // radius the apothem, touching every side's midpoint.
        let mut outline = Vec::new();
        for k in 0..corners.len() {
            let (p, q) = (corners[k], corners[(k + 1) % corners.len()]);
            for s in 0..10 {
                outline.push(p + (q - p) * (f64::from(s) / 10.0));
            }
        }
        let inner = zone::inscribed_circle(&outline, T).unwrap();
        let apothem = 7.0 * (core::f64::consts::PI / f64::from(sides)).cos();
        assert!(inner.shape.centre().distance(c) < 1e-9, "{sides} sides");
        assert!(
            (inner.shape.radius() - apothem).abs() < 1e-9,
            "{sides} sides"
        );
        assert_eq!(inner.contacts.len(), sides as usize, "{sides} sides");
    }
}

/// A square's outline: the inscribed circle is centred, radius half the
/// side.
#[test]
fn a_squares_inscribed_circle_is_half_its_side() {
    let mut outline = Vec::new();
    for k in 0..40 {
        let t = f64::from(k) * 0.25;
        outline.extend([
            Point2::new(t, 0.0),
            Point2::new(10.0, t),
            Point2::new(10.0 - t, 10.0),
            Point2::new(0.0, 10.0 - t),
        ]);
    }
    let zone = zone::inscribed_circle(&outline, T).unwrap();
    assert!(zone.shape.centre().distance(Point2::new(5.0, 5.0)) < 1e-12);
    assert!((zone.shape.radius() - 5.0).abs() < 1e-12);
    assert_eq!(zone.contacts.len(), 4);
}

fn nearest(points: &[Point2], c: Point2) -> f64 {
    points
        .iter()
        .fold(f64::INFINITY, |m, p| m.min(p.distance(c)))
}

fn farthest(points: &[Point2], c: Point2) -> f64 {
    points.iter().fold(0.0, |m, p| m.max(p.distance(c)))
}

/// Scattered points: the enclosing circle holds them and no centre does
/// with less; the empty circle holds none of them and no centre in the
/// hull, from a dense grid over the hull's box, is farther from them all.
#[test]
fn scattered_points_have_the_extreme_circles() {
    let mut draw = Seeded(23);
    for round in 0..10 {
        let points: Vec<Point2> = (0..60)
            .map(|_| Point2::new(draw.range(-10.0, 10.0), draw.range(-6.0, 6.0)))
            .collect();
        let outer = zone::circumscribed_circle(&points, T).unwrap();
        let r = outer.shape.radius();
        assert!((farthest(&points, outer.shape.centre()) - r).abs() < ROUNDING);
        for scale in SCALES {
            for _ in 0..200 {
                let c =
                    outer.shape.centre() + Vector2::new(draw.gaussian(), draw.gaussian()) * scale;
                assert!(farthest(&points, c) >= r - ROUNDING, "round {round}");
            }
        }

        let inner = zone::inscribed_circle(&points, T).unwrap();
        let r = inner.shape.radius();
        assert!((nearest(&points, inner.shape.centre()) - r).abs() < ROUNDING);
        let hull = hull_lines(&points);
        for i in 0..=400 {
            for j in 0..=240 {
                let q = Point2::new(-10.0 + f64::from(i) * 0.05, -6.0 + f64::from(j) * 0.05);
                if hull.iter().all(|(a, e)| e.cross(q - *a) >= 0.0) {
                    assert!(nearest(&points, q) <= r + ROUNDING, "round {round}: {q:?}");
                }
            }
        }
    }
}

/// The lines through two points with every point on their left: the
/// convex hull's sides, a point inside left of each.
fn hull_lines(points: &[Point2]) -> Vec<(Point2, Vector2)> {
    let mut out = Vec::new();
    for a in points {
        for b in points {
            let e = *b - *a;
            if e.magnitude() > 0.0 && points.iter().all(|p| e.cross(*p - *a) >= 0.0) {
                out.push((*a, e));
            }
        }
    }
    out
}

// ---- spheres ----

fn shell_width(points: &[Point], c: Point) -> f64 {
    let (lo, hi) = points
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
            let r = p.distance(c);
            (lo.min(r), hi.max(r))
        });
    hi - lo
}

fn on_sphere(count: usize, c: Point, mut radius: impl FnMut(Vector) -> f64) -> Vec<Point> {
    let golden = core::f64::consts::PI * (3.0 - 5.0_f64.sqrt());
    (0..count)
        .map(|i| {
            let z = 1.0 - 2.0 * (i as f64 + 0.5) / count as f64;
            let rho = (1.0 - z * z).sqrt();
            let a = golden * i as f64;
            let u = Vector::new(rho * a.cos(), rho * a.sin(), z);
            c + u * radius(u)
        })
        .collect()
}

#[test]
fn an_exact_sphere_is_spherical() {
    let c = Point::new(1.0, 2.0, -3.0);
    let zone = zone::sphericity(&on_sphere(800, c, |_| 12.0), T).unwrap();
    assert!(zone.width < 1e-12, "{}", zone.width);
    assert!(zone.optimal);
}

#[test]
fn a_noisy_ball_is_least_spherical() {
    let c = Point::new(1.0, 2.0, -3.0);
    let mut draw = Seeded(29);
    let points = on_sphere(3000, c, |u| {
        12.0 + 0.01 * u.x * u.y * 4.0 + 0.003 * draw.gaussian()
    });
    let zone = zone::sphericity(&points, T).unwrap();
    assert!(zone.optimal);
    let ls = fit::sphere(&points, &FitOptions::default(), T).unwrap();
    assert!(
        zone.width <= ls.form.1 - ls.form.0,
        "{} against {:?}",
        zone.width,
        ls.form
    );
    assert!(zone.contacts.len() >= 5, "{:?}", zone.contacts);
    let at = zone.shape.centre();
    for scale in SCALES {
        for _ in 0..300 {
            let moved = at + draw.direction() * (scale * draw.unit());
            assert!(
                shell_width(&points, moved) >= zone.width - ROUNDING,
                "a move of {scale}"
            );
        }
    }
}

// ---- cylinders and lines ----

fn distances_to(points: &[Point], axis: Axis) -> (f64, f64) {
    points
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
            let d = axis.distance_to(*p);
            (lo.min(d), hi.max(d))
        })
}

/// The axis moved across itself and tilted, each by `scale` in a random
/// direction, about the point level with the middle of a part `length`
/// long.
fn perturbed(axis: Axis, scale: f64, draw: &mut Seeded) -> Axis {
    let z = axis.direction.vector();
    let shift = draw.direction();
    let shift = shift - z * shift.dot(z);
    let tilt = draw.direction();
    let tilt = tilt - z * tilt.dot(z);
    let location = axis.location + shift * (scale * draw.unit());
    let direction = Direction::new(z + tilt * (scale * draw.unit() / 20.0), T).unwrap();
    Axis::new(location, direction)
}

fn bore(draw: &mut Seeded, count: usize, form: f64, noise: f64) -> (Frame, Vec<Point>) {
    let frame = tilted();
    let points = (0..count)
        .map(|_| {
            let a = draw.range(0.0, core::f64::consts::TAU);
            let h = draw.range(-20.0, 20.0);
            // A slight taper and a two-lobed section under the noise.
            let r = 8.0 + form * (h / 20.0 + (2.0 * a).cos()) + noise * draw.gaussian();
            place(frame, r * a.cos(), r * a.sin(), h)
        })
        .collect();
    (frame, points)
}

#[test]
fn an_exact_cylinder_is_cylindrical() {
    let mut draw = Seeded(31);
    let (frame, points) = bore(&mut draw, 2000, 0.0, 0.0);
    let zone = zone::cylindricity(&points, T).unwrap();
    assert!(zone.width < 1e-12, "{}", zone.width);
    assert!(zone.optimal);
    assert!((zone.shape.radius() - 8.0).abs() < 1e-12);
    for fit in [
        zone::inscribed_cylinder(&points, T).unwrap(),
        zone::circumscribed_cylinder(&points, T).unwrap(),
    ] {
        assert!((fit.shape.radius() - 8.0).abs() < 1e-12);
        assert!(parallel(fit.shape.frame().z(), frame.z()));
        assert!(fit.optimal);
    }
}

#[test]
fn a_noisy_bore_is_least_cylindrical() {
    let mut draw = Seeded(37);
    let (_, points) = bore(&mut draw, 4000, 0.005, 0.002);
    let zone = zone::cylindricity(&points, T).unwrap();
    assert!(zone.optimal);
    let ls = fit::cylinder(&points, &FitOptions::default(), T).unwrap();
    assert!(
        zone.width <= ls.form.1 - ls.form.0,
        "{} against {:?}",
        zone.width,
        ls.form
    );
    assert!(zone.contacts.len() >= 5, "{:?}", zone.contacts);
    let axis = zone.shape.frame().axis();
    let (lo, hi) = distances_to(&points, axis);
    assert!((hi - lo - zone.width).abs() < ROUNDING);
    for scale in SCALES {
        for _ in 0..300 {
            let (lo, hi) = distances_to(&points, perturbed(axis, scale, &mut draw));
            assert!(
                hi - lo >= zone.width - ROUNDING,
                "a move of {scale}: {}",
                hi - lo
            );
        }
    }
}

/// No move of the axis lets a larger cylinder fit inside a noisy bore, or
/// a smaller one round it.
#[test]
fn a_noisy_bore_has_its_extreme_cylinders() {
    let mut draw = Seeded(41);
    let (_, points) = bore(&mut draw, 4000, 0.005, 0.002);
    let inner = zone::inscribed_cylinder(&points, T).unwrap();
    let outer = zone::circumscribed_cylinder(&points, T).unwrap();
    assert!(inner.optimal && outer.optimal);
    assert!(inner.shape.radius() < outer.shape.radius());
    for scale in SCALES {
        for _ in 0..300 {
            let (lo, _) = distances_to(
                &points,
                perturbed(inner.shape.frame().axis(), scale, &mut draw),
            );
            assert!(lo <= inner.shape.radius() + ROUNDING, "a move of {scale}");
            let (_, hi) = distances_to(
                &points,
                perturbed(outer.shape.frame().axis(), scale, &mut draw),
            );
            assert!(hi >= outer.shape.radius() - ROUNDING, "a move of {scale}");
        }
    }
}

/// A helix of radius `a` about a line is straight to within a cylinder of
/// diameter `2a`; a line's points to within nothing.
#[test]
fn a_helix_is_as_straight_as_its_diameter() {
    let frame = tilted();
    let helix: Vec<Point> = (0..2000)
        .map(|i| {
            let t = f64::from(i) * 0.01;
            place(frame, 0.05 * (7.0 * t).cos(), 0.05 * (7.0 * t).sin(), t)
        })
        .collect();
    let zone = zone::straightness(&helix, T).unwrap();
    assert!(zone.optimal);
    assert!((zone.width - 0.1).abs() < 1e-9, "{}", zone.width);
    assert!(parallel(zone.shape.direction, frame.z()));

    let line: Vec<Point> = (0..100)
        .map(|i| place(frame, 0.0, 0.0, f64::from(i) * 0.3))
        .collect();
    let zone = zone::straightness(&line, T).unwrap();
    assert!(zone.width < 1e-12, "{}", zone.width);
}

#[test]
fn a_noisy_line_is_least_straight() {
    let frame = tilted();
    let mut draw = Seeded(43);
    let points: Vec<Point> = (0..1000)
        .map(|i| {
            let t = f64::from(i) * 0.04;
            place(
                frame,
                0.003 * draw.gaussian() + 1e-4 * t,
                0.003 * draw.gaussian(),
                t,
            )
        })
        .collect();
    let zone = zone::straightness(&points, T).unwrap();
    assert!(zone.optimal);
    let ls = fit::line3(&points, &FitOptions::default(), T).unwrap();
    assert!(
        zone.width <= 2.0 * ls.form.1,
        "{} against {:?}",
        zone.width,
        ls.form
    );
    for scale in SCALES {
        for _ in 0..300 {
            let (_, hi) = distances_to(&points, perturbed(zone.shape, scale, &mut draw));
            assert!(2.0 * hi >= zone.width - ROUNDING, "a move of {scale}");
        }
    }
}

/// The same points give the same zones, bit for bit.
#[test]
fn zones_are_deterministic() {
    let mut draw = Seeded(47);
    let (_, points) = bore(&mut draw, 2000, 0.005, 0.002);
    assert_eq!(
        zone::cylindricity(&points, T).unwrap(),
        zone::cylindricity(&points, T).unwrap()
    );
    let flat: Vec<Point2> = (0..1000)
        .map(|i| {
            on_circle(
                Point2::new(1.0, 2.0),
                9.0 + 0.003 * draw.gaussian(),
                f64::from(i) * 0.00628,
            )
        })
        .collect();
    assert_eq!(
        zone::circularity(&flat, T).unwrap(),
        zone::circularity(&flat, T).unwrap()
    );
    assert_eq!(
        zone::flatness(&points, T).unwrap(),
        zone::flatness(&points, T).unwrap()
    );
}
