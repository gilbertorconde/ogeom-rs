#![allow(
    clippy::unwrap_used,
    reason = "test code; a failed unwrap is a failed test"
)]

//! An oblique plane and a cone meet in a conic, stated exactly.

use ogeom_core::Tolerances;
use ogeom_geom::{ConeSurface, Curve, Curve3d as _, PlaneSurface, Surface as _, SurfaceGeometry};
use ogeom_intersect::{IntersectOptions, SurfaceIntersection, intersect_surfaces};
use ogeom_math::{Cone, Direction, Frame, Plane, Point, Vector};

const T: Tolerances = Tolerances::millimetres();

/// A small deterministic generator, so a failure names its case.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        #[allow(clippy::cast_precision_loss)]
        let unit = (self.0 >> 11) as f64 / (1_u64 << 53) as f64;
        unit
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next()
    }

    fn direction(&mut self) -> Vector {
        loop {
            let v = Vector::new(
                self.range(-1.0, 1.0),
                self.range(-1.0, 1.0),
                self.range(-1.0, 1.0),
            );
            if v.magnitude() > 0.2 {
                return v;
            }
        }
    }
}

/// Every oblique plane through a cone's window, steeper or shallower than
/// its wall: each conic stated lies on both surfaces, is exact, and the
/// points where the plane crosses the window's rulings all lie on one of
/// them, so no branch is missing.
#[test]
fn an_oblique_plane_cuts_a_cone_in_exact_conics() {
    let mut rng = Lcg(7);
    let mut kinds = [0_usize; 3];
    for case in 0..300 {
        let z = rng.direction();
        let frame = Frame::new(
            Point::new(
                rng.range(-2.0, 2.0),
                rng.range(-2.0, 2.0),
                rng.range(-2.0, 2.0),
            ),
            Direction::new(z, T).unwrap(),
            Direction::from_cross(z, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
            T,
        )
        .unwrap();
        let half = rng.range(0.15, 1.2) * if rng.next() < 0.5 { -1.0 } else { 1.0 };
        let height = rng.range(1.0, 6.0);
        let reference = rng.range(0.5, 3.0)
            + if half < 0.0 {
                height * (-half).tan()
            } else {
                0.0
            };
        let cone = Cone::new(frame, reference, half, T).unwrap();
        let surface: SurfaceGeometry = ConeSurface::new(cone, (0.0, height)).unwrap().into();
        // A plane through a point inside the window, leaning anyhow.
        let through = frame.origin() + frame.z().vector() * rng.range(0.1, height - 0.1);
        let normal = Direction::new(rng.direction(), T).unwrap();
        if (normal.vector().dot(frame.z().vector()).abs() - 1.0).abs() < 1e-6 {
            continue;
        }
        let plane = Plane::through(through, normal);
        let flat: SurfaceGeometry = PlaneSurface::over(plane, (-50.0, 50.0), (-50.0, 50.0))
            .unwrap()
            .into();
        let found = intersect_surfaces(&flat, &surface, IntersectOptions::default(), T).unwrap();
        let SurfaceIntersection::Along(curves) = found else {
            panic!("case {case}: a plane through the window meets the cone: {found:?}");
        };
        let mut samples: Vec<Point> = Vec::new();
        for c in &curves {
            assert!(c.exact, "case {case}: {:?}", c.curve);
            match &c.curve {
                Curve::Ellipse(_) => kinds[0] += 1,
                Curve::Parabola(_) => kinds[1] += 1,
                Curve::Hyperbola(_) => kinds[2] += 1,
                other => panic!("case {case}: {other:?}"),
            }
            let (lo, hi) = c.curve.domain();
            for i in 0..=2000 {
                let t = lo + (hi - lo) * f64::from(i) / 2000.0;
                let p = c.curve.point_at(t, T).unwrap();
                assert!(plane.distance_to(p) < 1e-9, "case {case}: off the plane");
                assert!(
                    cone.distance_to(p) < 1e-9,
                    "case {case}: off the cone by {}",
                    cone.distance_to(p)
                );
                samples.push(p);
            }
        }
        // Where the plane crosses each ruling inside the window, a conic
        // passes.
        for i in 0..64 {
            let u = core::f64::consts::TAU * f64::from(i) / 64.0;
            let at = |v: f64| surface.point_at(u, v, T).unwrap();
            let side = |v: f64| plane.signed_distance_to(at(v));
            let (mut a, mut b) = (0.0, height);
            if side(a) * side(b) > 0.0 {
                continue;
            }
            for _ in 0..80 {
                let m = 0.5 * (a + b);
                if side(a) * side(m) <= 0.0 {
                    b = m;
                } else {
                    a = m;
                }
            }
            let p = at(0.5 * (a + b));
            let near = samples
                .iter()
                .map(|q| q.distance(p))
                .fold(f64::INFINITY, f64::min);
            let spacing = curves
                .iter()
                .map(|c| {
                    let (lo, hi) = c.curve.domain();
                    c.curve
                        .point_at(lo, T)
                        .unwrap()
                        .distance(c.curve.point_at(lo + (hi - lo) / 2000.0, T).unwrap())
                })
                .fold(0.0_f64, f64::max);
            assert!(
                near <= spacing.max(1e-6) * 4.0,
                "case {case}: ruling {i} crossing {p:?} on no conic ({near})"
            );
        }
    }
    assert!(
        kinds[0] > 0 && kinds[2] > 0,
        "both closed and open conics drawn: {kinds:?}"
    );
}

/// A plane parallel to a ruling cuts a parabola.
#[test]
fn a_plane_along_a_ruling_cuts_a_parabola() {
    let half = 0.5_f64;
    let cone = Cone::new(Frame::WORLD, 1.0, half, T).unwrap();
    let surface: SurfaceGeometry = ConeSurface::new(cone, (0.0, 4.0)).unwrap().into();
    // The ruling at x = 1 + z tan(half), y = 0 leans by `half` off the
    // axis; the plane holding y and that slope, set back from it, is
    // parallel to it.
    let normal = Vector::new(half.cos(), 0.0, -half.sin());
    let plane = Plane::through(
        Point::new(0.4, 0.0, 0.0),
        Direction::new(normal, T).unwrap(),
    );
    let flat: SurfaceGeometry = PlaneSurface::over(plane, (-50.0, 50.0), (-50.0, 50.0))
        .unwrap()
        .into();
    let found = intersect_surfaces(&flat, &surface, IntersectOptions::default(), T).unwrap();
    let SurfaceIntersection::Along(curves) = found else {
        panic!("{found:?}");
    };
    assert_eq!(curves.len(), 1);
    assert!(
        matches!(curves[0].curve, Curve::Parabola(_)),
        "{:?}",
        curves[0].curve
    );
    let (lo, hi) = curves[0].curve.domain();
    for i in 0..=200 {
        let p = curves[0]
            .curve
            .point_at(lo + (hi - lo) * f64::from(i) / 200.0, T)
            .unwrap();
        assert!(plane.distance_to(p) < 1e-9 && cone.distance_to(p) < 1e-9);
    }
}
