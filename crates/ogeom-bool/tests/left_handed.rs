//! Booleans between solids on surfaces of opposite handedness: one solid's
//! wall on a right-handed frame, the other's on the same surface spelled on
//! a left-handed (mirrored) frame, whose angle runs the other way round.

#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::PI;

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{
    CircleCurve, ConeSurface, Curve, CylinderSurface, SphereSurface, Surface as _, SurfaceGeometry,
    TorusSurface,
};
use ogeom_math::{Circle, Cone, Cylinder, Frame, Point, Sphere, Torus, Transform, Vector};
use ogeom_topo::{Model, Shape};

const T: Tolerances = Tolerances::millimetres();

/// The two left-handed spellings of a right-handed frame: `y` reversed,
/// where the angle runs the other way round the same axis, and `z`
/// reversed, where it runs the same way round an axis pointing the other
/// way.
#[derive(Clone, Copy, Debug)]
enum Mirror {
    Y,
    Z,
}

impl Mirror {
    fn frame(self, f: Frame) -> Frame {
        match self {
            Self::Y => f.mirrored(),
            Self::Z => Frame::from_axes(f.origin(), f.x(), f.y(), -f.z(), T).unwrap(),
        }
    }

    /// A height window read on the mirrored frame.
    fn heights(self, (lo, hi): (f64, f64)) -> (f64, f64) {
        match self {
            Self::Y => (lo, hi),
            Self::Z => (-hi, -lo),
        }
    }
}

/// The solid with every cylinder, cone, sphere and torus restated on the
/// mirrored spelling of its frame, the same point set, and every circle
/// on its left-handed spelling, the same point at every parameter.
fn left_handed(model: &mut Model, solid: &Shape, mirror: Mirror) -> Shape {
    use SurfaceGeometry as S;
    let surface = |s: &S| -> OgeomResult<Option<(S, bool)>> {
        let (_, window) = s.domain();
        let window = mirror.heights(window);
        Ok(match s {
            S::Cylinder(c) => {
                let c = c.cylinder();
                let c = Cylinder::new(mirror.frame(c.frame()), c.radius(), T)?;
                Some((CylinderSurface::new(c, window)?.into(), true))
            }
            S::Cone(c) => {
                let c = c.cone();
                let angle = match mirror {
                    Mirror::Y => c.half_angle(),
                    Mirror::Z => -c.half_angle(),
                };
                let c = Cone::new(mirror.frame(c.frame()), c.reference_radius(), angle, T)?;
                Some((ConeSurface::new(c, window)?.into(), true))
            }
            S::Sphere(c) => {
                let c = c.sphere();
                let c = Sphere::new(mirror.frame(c.frame()), c.radius(), T)?;
                Some((SphereSurface::new(c).into(), true))
            }
            S::Torus(c) => {
                let c = c.torus();
                let c = Torus::new(
                    mirror.frame(c.frame()),
                    c.major_radius(),
                    c.minor_radius(),
                    T,
                )?;
                Some((TorusSurface::new(c).into(), true))
            }
            _ => None,
        })
    };
    let curve = |c: &Curve, range: (f64, f64)| -> OgeomResult<Option<(Curve, (f64, f64))>> {
        Ok(match c {
            Curve::Circle(c) if !c.is_reversed() => {
                let circle = c.circle();
                let frame = Mirror::Z.frame(circle.frame());
                let circle = Circle::new(frame, circle.radius(), T)?;
                Some((CircleCurve::new(circle).into(), range))
            }
            _ => None,
        })
    };
    let restated = ogeom_algo::restate_geometry(model, solid, &surface, &curve, T)
        .unwrap()
        .shape;
    assert_valid(model, &restated, "the left-handed spelling");
    restated
}

fn volume(model: &Model, shape: &Shape) -> f64 {
    let fine = ogeom_mesh::Deflection {
        chord: 1e-3,
        ..ogeom_mesh::Deflection::default()
    };
    ogeom_algo::volume_properties(model, shape, fine, T)
        .unwrap()
        .mass
}

fn assert_valid(model: &Model, shape: &Shape, what: &str) {
    let diagnosis = ogeom_algo::check(model, shape, T).unwrap();
    assert!(diagnosis.is_valid(), "{what}: {:?}", diagnosis.problems);
}

/// Each operation's result is valid and holds `want`. Every wall here is
/// a closed form, or a fit held to the boolean's own tolerance where a
/// ball dents one, so the mass misses only by rounding and the fit: a
/// wall lost, doubled or turned over misses by a dent's volume or more.
fn assert_ops(model: &mut Model, a: &Shape, b: &Shape, want: [f64; 4], case: &str) {
    for (k, op) in ["fuse", "common", "a less b", "b less a"]
        .iter()
        .enumerate()
    {
        let built = match k {
            0 => ogeom_bool::fuse(model, a, b, T),
            1 => ogeom_bool::common(model, a, b, T),
            2 => ogeom_bool::cut(model, a, b, T),
            _ => ogeom_bool::cut(model, b, a, T),
        }
        .unwrap_or_else(|e| panic!("{case}, {op}: {e}"));
        let what = format!("{case}, {op}");
        assert_valid(model, &built.shape, &what);
        let v = volume(model, &built.shape);
        assert!(
            (v - want[k]).abs() < 1e-8,
            "{what}: {v} against {}",
            want[k]
        );
    }
}

/// Two drums of radius 1 and height 2 on one axis, the second slid a unit
/// up or down it and spelled on a left-handed frame: the union holds
/// `3π`, the overlap `π`, and each less the other `π`. Either operand
/// first.
#[test]
fn coaxial_drums_of_opposite_handedness_overlapping_along_the_axis() {
    for mirror in [Mirror::Y, Mirror::Z] {
        for slide in [1.0, -1.0] {
            let mut model = Model::new();
            let a = ogeom_algo::make_cylinder(&mut model, Frame::WORLD, 1.0, 2.0, T)
                .unwrap()
                .shape;
            let moved = ogeom_algo::transformed(
                &mut model,
                &a,
                Transform::translation(Vector::new(0.0, 0.0, slide)),
            )
            .unwrap()
            .shape;
            let b = left_handed(&mut model, &moved, mirror);
            let case = format!("{mirror:?} mirror, slid {slide}");
            assert_ops(&mut model, &a, &b, [3.0 * PI, PI, PI, PI], &case);
            assert_ops(
                &mut model,
                &b,
                &a,
                [3.0 * PI, PI, PI, PI],
                &format!("{case}, left first"),
            );
        }
    }
}

/// A cylinder, a cone frustum, a ball and a torus, each dented by a small
/// ball off its axis, so its wall carries a fitted section, against the
/// whole solid spelled on the frame with `y` reversed, its axis kept (the
/// drum also slid half a unit up it). The union is the whole solid's
/// closed form, the overlap the dented solid's volume less what lies below
/// the slid drum, and the whole less the dented solid what remains.
#[test]
fn a_dented_quadric_against_the_whole_on_a_mirrored_frame() {
    let frustum = PI * 2.0 / 3.0 * (1.0 + 0.5 + 0.25);
    let cases: [(&str, f64, Point, f64); 5] = [
        ("drum", 0.0, Point::new(1.0, 0.0, 1.0), 2.0 * PI),
        ("drum slid", 0.5, Point::new(1.0, 0.0, 1.0), 2.0 * PI),
        ("frustum", 0.0, Point::new(0.75, 0.0, 1.0), frustum),
        ("ball", 0.0, Point::new(1.0, 0.0, 0.0), 4.0 / 3.0 * PI),
        (
            "torus",
            0.0,
            Point::new(2.5, 0.0, 0.0),
            2.0 * PI * PI * 2.0 * 0.25,
        ),
    ];
    let mirror = Mirror::Y;
    for (kind, slide, dent_at, whole) in cases {
        let mut model = Model::new();
        let m = &mut model;
        let solid = match kind {
            "frustum" => ogeom_algo::make_cone(m, Frame::WORLD, 1.0, 0.5, 2.0, T),
            "ball" => ogeom_algo::make_sphere(m, Frame::WORLD, 1.0, T),
            "torus" => ogeom_algo::make_torus(m, Frame::WORLD, 2.0, 0.5, T),
            _ => ogeom_algo::make_cylinder(m, Frame::WORLD, 1.0, 2.0, T),
        }
        .unwrap()
        .shape;
        let ball = ogeom_algo::make_sphere(m, Frame::WORLD.with_origin(dent_at), 0.3, T)
            .unwrap()
            .shape;
        let dented = ogeom_bool::cut(m, &solid, &ball, T).unwrap().shape;
        let kept = volume(m, &dented);
        let moved = ogeom_algo::transformed(
            m,
            &solid,
            Transform::translation(Vector::new(0.0, 0.0, slide)),
        )
        .unwrap()
        .shape;
        let other = left_handed(m, &moved, mirror);
        // The slid drum leaves the dented one's undented foot below
        // it; the dent lies within both.
        let overlap = kept - PI * slide;
        let union = whole + PI * slide;
        let want = [union, overlap, kept - overlap, whole - overlap];
        assert_ops(
            m,
            &dented,
            &other,
            want,
            &format!("{mirror:?} mirror, {kind}"),
        );
    }
}
