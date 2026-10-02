#![allow(
    clippy::unwrap_used,
    reason = "test code; a failed unwrap is a failed test"
)]

//! Configurations the surface/surface literature names as hard.
//!
//! Every other test's inputs are chosen by the people who wrote the code.
//! These are transcribed from the published record of what breaks
//! intersectors.

use ogeom_core::Tolerances;
use ogeom_geom::{
    ConeSurface, Curve3d as _, CylinderSurface, PlaneSurface, SphereSurface, SurfaceGeometry,
    TorusSurface,
};
mod support;

use ogeom_intersect::{Marching, Meeting, Stopped, branches, surface_surface};
use ogeom_math::{Cone, Cylinder, Direction, Frame, Plane, Point, Sphere, Torus, Vector};
use support::coverage::coverage;

const T: Tolerances = Tolerances::millimetres();

fn frame(origin: Point, z: Vector) -> Frame {
    Frame::new(
        origin,
        Direction::new(z, T).unwrap(),
        Direction::from_cross(z, Vector::new(0.3, 0.5, 0.9), T).unwrap(),
        T,
    )
    .unwrap()
}

fn cyl(origin: Point, axis: Vector, r: f64) -> SurfaceGeometry {
    CylinderSurface::new(
        Cylinder::new(frame(origin, axis), r, T).unwrap(),
        (-6.0, 6.0),
    )
    .unwrap()
    .into()
}

fn sph(c: Point, r: f64) -> SurfaceGeometry {
    SphereSurface::new(Sphere::centred(c, r, T).unwrap()).into()
}

fn pln(o: Point, n: Vector) -> SurfaceGeometry {
    PlaneSurface::over(
        Plane::through(o, Direction::new(n, T).unwrap()),
        (-8.0, 8.0),
        (-8.0, 8.0),
    )
    .unwrap()
    .into()
}

fn off(s: &SurfaceGeometry, p: Point) -> f64 {
    match s {
        SurfaceGeometry::Plane(x) => x.plane().distance_to(p),
        SurfaceGeometry::Sphere(x) => x.sphere().distance_to(p),
        SurfaceGeometry::Cylinder(x) => x.cylinder().distance_to(p),
        SurfaceGeometry::Torus(x) => x.torus().distance_to(p),
        SurfaceGeometry::Cone(x) => x.cone().distance_to(p),
        _ => 0.0,
    }
}

fn worst(a: &SurfaceGeometry, b: &SurfaceGeometry, found: &[ogeom_intersect::Traced]) -> f64 {
    found
        .iter()
        .flat_map(|br| br.points.iter())
        .map(|p| off(a, *p).abs().max(off(b, *p).abs()))
        .fold(0.0_f64, f64::max)
}

fn options() -> Marching {
    Marching {
        chord: 1e-5,
        ..Marching::default()
    }
}

#[test]
fn equal_cylinders_at_a_shallow_angle_cover_everything_even_fragmented() {
    // The classic: equal radii force the two intersection curves through two
    // tangency points, and a two-degree crossing makes them long and thin.
    // The tracer stalls at the tangencies rather than jumping branch (the
    // deliberate refusal), so the answer arrives as fragments. What must hold
    // even so: every fragment on both surfaces, and nothing missed.
    let a = cyl(Point::ORIGIN, Vector::Z, 1.0);
    let tilt = 2.0_f64.to_radians();
    let b = cyl(Point::ORIGIN, Vector::new(tilt.sin(), 0.0, tilt.cos()), 1.0);

    let found = branches(&a, &b, options(), T).unwrap();
    assert!(!found.is_empty());
    assert!(
        found.iter().all(|br| br.stopped == Stopped::Stalled),
        "each fragment ends at a tangency it refuses to march through"
    );
    assert!(worst(&a, &b, &found) < 1e-7);
    let score = coverage(&a, &b, &found, 50, T).unwrap();
    assert!(
        score.complete(),
        "{}/{}: fragmented is acceptable, incomplete is not",
        score.covered,
        score.crossings
    );
}

#[test]
fn a_plane_tangent_along_a_ruling_is_the_analytic_paths_case() {
    // Tangential contact along a line. The marcher correctly finds no
    // transversal crossing; the analytic layer names the contact line. The
    // division of labour is the answer here, and both halves are pinned.
    let drum = cyl(Point::ORIGIN, Vector::Z, 2.0);
    let touching = pln(Point::new(2.0, 0.0, 0.0), Vector::X);

    assert!(branches(&drum, &touching, options(), T).unwrap().is_empty());
    assert!(matches!(
        surface_surface(&drum, &touching, T).unwrap(),
        Meeting::Along(ref c) if c.len() == 1
    ));
}

#[test]
fn near_tangent_sphere_and_cylinder_give_one_thin_loop() {
    // Equal radii, axis offset a thousandth: the near-tangential quartic. One
    // closed loop, thin in z, complete.
    let a = cyl(Point::new(1e-3, 0.0, 0.0), Vector::Z, 2.0);
    let b = sph(Point::ORIGIN, 2.0);
    let found = branches(&a, &b, options(), T).unwrap();
    assert_eq!(found.len(), 1);
    assert!(found[0].closed());
    assert!(worst(&a, &b, &found) < 1e-7);
    assert!(coverage(&a, &b, &found, 50, T).unwrap().complete());
}

#[test]
fn a_sphere_across_a_cones_apex_cuts_both_nappes() {
    // A cone's height range crosses its apex, so the surface has two nappes,
    // and a sphere spanning the apex cuts a loop in each. The instrument's
    // `Cone::distance_to` must measure both nappes, or a correctly traced
    // lower loop reads as far off the cone.
    let cone: SurfaceGeometry = ConeSurface::new(
        Cone::new(Frame::WORLD, 0.5, 0.5_f64.atan(), T).unwrap(),
        (-3.0, 3.0),
    )
    .unwrap()
    .into();
    let ball = sph(Point::new(0.0, 0.0, -0.9), 1.0);

    let found = branches(&cone, &ball, options(), T).unwrap();
    assert_eq!(found.len(), 2, "one loop per nappe");
    for br in &found {
        assert!(br.closed());
    }
    // One loop above the apex plane z = -1, one below.
    let sides: Vec<f64> = found.iter().map(|br| br.points[0].z + 1.0).collect();
    assert!(
        sides[0] * sides[1] < 0.0,
        "both loops on one nappe: z offsets {sides:?}"
    );
    assert!(worst(&cone, &ball, &found) < 1e-7);
}

/// A polynomial in `t`, ascending powers.
type Poly = Vec<f64>;

fn poly_mul(a: &[f64], b: &[f64]) -> Poly {
    let mut out = vec![0.0; a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

fn poly_add(a: &[f64], b: &[f64]) -> Poly {
    (0..a.len().max(b.len()))
        .map(|i| a.get(i).copied().unwrap_or(0.0) + b.get(i).copied().unwrap_or(0.0))
        .collect()
}

fn poly_scale(a: &[f64], k: f64) -> Poly {
    a.iter().map(|x| x * k).collect()
}

/// The Bernstein coefficients of a polynomial of degree at most four, over
/// `t` in `[-1, 1]`.
fn bernstein4(p: &[f64]) -> [f64; 5] {
    // In x = (t + 1) / 2, then to the Bernstein basis.
    let mut in_x = vec![0.0];
    let mut power = vec![1.0];
    for c in p {
        in_x = poly_add(&in_x, &poly_scale(&power, *c));
        power = poly_mul(&power, &[-1.0, 2.0]);
    }
    let binomial =
        |n: u32, k: u32| (0..k).fold(1.0, |acc, i| acc * f64::from(n - i) / f64::from(i + 1));
    let mut out = [0.0; 5];
    for (j, slot) in (0_u32..).zip(out.iter_mut()) {
        for (i, c) in (0..=j).zip(in_x.iter()) {
            *slot += binomial(j, i) / binomial(4, i) * c;
        }
    }
    out
}

/// A free-form cradle under the world torus of radii `major`, `minor`: the
/// rational patch ruled along the meridian tangents of a curve on the torus
/// that wanders across the tube.
///
/// With `a = tan(u/2)` and `b = tan(v/2)` the torus is rational in both, and
/// the curve `a = t`, `b = b0 + slope * t` on it is rational of degree four,
/// as is the tube's meridian tangent along it. The patch is
/// `c(t) + s * w(t)`, degree four by one, exact. Along `s = 0` it shares the
/// torus's tangent plane, and across the curve it is straight where the tube
/// is round, so it touches the torus along the whole curve and nowhere
/// crosses it. The test knows the contact from the construction; the
/// intersector is given only the two surfaces.
fn cradle(major: f64, minor: f64, b0: f64, slope: f64, s: (f64, f64)) -> SurfaceGeometry {
    use ogeom_geom::BSplineSurface;
    use ogeom_math::{ControlGrid, KnotVector, Weighted};
    let a = [0.0, 1.0];
    let b = [b0, slope];
    let a2 = poly_mul(&a, &a);
    let b2 = poly_mul(&b, &b);
    let one_plus_a2 = poly_add(&[1.0], &a2);
    let one_minus_a2 = poly_add(&[1.0], &poly_scale(&a2, -1.0));
    let one_plus_b2 = poly_add(&[1.0], &b2);
    let one_minus_b2 = poly_add(&[1.0], &poly_scale(&b2, -1.0));
    let weight = poly_mul(&one_plus_a2, &one_plus_b2);
    // (major + minor cos v) times (1 + b^2).
    let reach = poly_add(
        &poly_scale(&one_plus_b2, major),
        &poly_scale(&one_minus_b2, minor),
    );
    let curve = [
        poly_mul(&reach, &one_minus_a2),
        poly_mul(&reach, &poly_scale(&a, 2.0)),
        poly_scale(&poly_mul(&b, &one_plus_a2), 2.0 * minor),
    ];
    let tangent = [
        poly_scale(&poly_mul(&b, &one_minus_a2), -2.0),
        poly_scale(&poly_mul(&b, &a), -4.0),
        poly_mul(&one_minus_b2, &one_plus_a2),
    ];
    let w = bernstein4(&weight);
    let mut grid = Vec::new();
    for (i, wi) in w.iter().enumerate() {
        for sj in [s.0, s.1] {
            let at =
                |k: usize| bernstein4(&poly_add(&curve[k], &poly_scale(&tangent[k], sj)))[i] / wi;
            grid.push(Weighted::new(Point::new(at(0), at(1), at(2)), *wi, T).unwrap());
        }
    }
    BSplineSurface::rational(
        KnotVector::clamped_uniform(4, 5).unwrap(),
        KnotVector::clamped_uniform(1, 2).unwrap(),
        ControlGrid::new(grid, 5, 2).unwrap(),
    )
    .unwrap()
    .into()
}

/// One reported contact measured against both surfaces: each pcurve lifted
/// against the curve, and the two normals there. Returns the worst lift
/// distance and the worst normal sine.
fn measured_contact(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    contact: &ogeom_intersect::SectionCurve,
) -> (f64, f64) {
    use ogeom_geom::{Curve2d as _, Surface as _};
    let (lo, hi) = contact.curve.domain();
    let (mut off_surfaces, mut sine) = (0.0_f64, 0.0_f64);
    for k in 0..=2000 {
        let t = lo + (hi - lo) * f64::from(k) / 2000.0;
        let on = contact.curve.point_at(t, T).unwrap();
        let mut normals = Vec::new();
        for (surface, image) in [(a, &contact.on_a), (b, &contact.on_b)] {
            let at = image.as_ref().unwrap().point_at(t, T).unwrap();
            off_surfaces = off_surfaces.max(surface.point_at(at.x, at.y, T).unwrap().distance(on));
            normals.push(surface.normal_at(at.x, at.y, T).unwrap().vector());
        }
        sine = sine.max(normals[0].cross(normals[1]).magnitude());
    }
    (off_surfaces, sine)
}

/// The sine within which the two normals of a reported contact agree.
const CONTACT_SINE: f64 = 1e-3;

#[test]
fn a_torus_resting_in_a_free_form_cradle_touches_along_one_curve() {
    // No closed form names this contact: the cradle is a rational patch
    // whose contact with the tube wanders across it, from 77 degrees below
    // the tube's equator to 23. The crossing marcher stalls along it in
    // fragments; the tangential trace follows it as one curve, edge to edge
    // of the cradle.
    use ogeom_geom::{Curve2d as _, Surface as _};
    let torus: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, 3.0, 1.0, T).unwrap()).into();
    let cradle = cradle(3.0, 1.0, -0.5, 0.3, (-0.6, 0.6));
    // The construction's claim, measured: on the torus along s = 0 (the
    // patch's v = 1/2), outside it everywhere else.
    for i in 0..=20 {
        for j in 0..=8 {
            let p = cradle
                .point_at(f64::from(i) / 20.0, f64::from(j) / 8.0, T)
                .unwrap();
            let gap = off(&torus, p);
            assert!(gap > -1e-12 && (j != 4 || gap < 1e-12), "{i} {j}: {gap}");
        }
    }

    let found = ogeom_intersect::intersect_surfaces(
        &torus,
        &cradle,
        ogeom_intersect::IntersectOptions {
            tolerance: 1e-4,
            marching: options(),
        },
        T,
    )
    .unwrap();
    let ogeom_intersect::SurfaceIntersection::Along(curves) = found else {
        panic!("the cradle touches the torus along a curve: {found:?}");
    };
    assert_eq!(curves.len(), 1, "one contact, not fragments of it");
    let contact = &curves[0];
    assert!(contact.tangential && !contact.exact && !contact.closed);
    assert!(contact.tolerance < 1e-3, "stated {}", contact.tolerance);
    let (off_surfaces, sine) = measured_contact(&torus, &cradle, contact);
    assert!(
        off_surfaces <= contact.tolerance,
        "{off_surfaces} off, stated {}",
        contact.tolerance
    );
    assert!(sine <= CONTACT_SINE, "the normals part by sine {sine}");
    // The contact the construction put there: the cradle's middle ruling
    // line, run from one end of the cradle to the other.
    let image = contact.on_b.as_ref().unwrap();
    let (lo, hi) = image.domain();
    let ends = [
        image.point_at(lo, T).unwrap(),
        image.point_at(hi, T).unwrap(),
    ];
    assert!(
        ends[0].x.min(ends[1].x).abs() < 1e-6 && (ends[0].x.max(ends[1].x) - 1.0).abs() < 1e-6,
        "edge to edge: {ends:?}"
    );
    for k in 0..=200 {
        let at = image
            .point_at(lo + (hi - lo) * f64::from(k) / 200.0, T)
            .unwrap();
        assert!((at.y - 0.5).abs() < 1e-3, "off the contact line: {at:?}");
    }
}

#[test]
fn a_ball_resting_in_a_spline_drum_touches_along_its_great_circle() {
    // A unit ball inside a unit cylinder, the cylinder given as its exact
    // rational patch so that no closed form applies: the contact is the
    // ball's great circle across the axis, and the trace is measured
    // against it.
    let axis = Vector::new(0.2, -0.3, 1.0);
    let origin = Point::new(1.0, 2.0, 3.0);
    let drum: SurfaceGeometry = cyl(origin, axis, 1.0).to_bspline(T).unwrap().into();
    let unit = axis * (1.0 / axis.magnitude());
    let centre = origin + unit * 0.7;
    let ball = sph(centre, 1.0);

    let found = ogeom_intersect::intersect_surfaces(
        &drum,
        &ball,
        ogeom_intersect::IntersectOptions {
            tolerance: 1e-4,
            marching: options(),
        },
        T,
    )
    .unwrap();
    let ogeom_intersect::SurfaceIntersection::Along(curves) = found else {
        panic!("the ball touches the drum along a circle: {found:?}");
    };
    assert_eq!(curves.len(), 1, "one contact");
    let contact = &curves[0];
    assert!(contact.tangential && contact.closed);
    assert!(contact.tolerance < 1e-3, "stated {}", contact.tolerance);
    let (off_surfaces, sine) = measured_contact(&drum, &ball, contact);
    assert!(
        off_surfaces <= contact.tolerance,
        "{off_surfaces} off, stated {}",
        contact.tolerance
    );
    assert!(sine <= CONTACT_SINE, "the normals part by sine {sine}");
    // Against the closed form: radius 1 about the axis, in the plane
    // across it through the ball's centre, all the way round.
    let (lo, hi) = contact.curve.domain();
    let mut length = 0.0;
    let mut previous: Option<Point> = None;
    for k in 0..=2000 {
        let p = contact
            .curve
            .point_at(lo + (hi - lo) * f64::from(k) / 2000.0, T)
            .unwrap();
        let from = p - centre;
        let height = from.dot(unit);
        let radius = (from - unit * height).magnitude();
        assert!(
            (radius - 1.0).abs() <= contact.tolerance && height.abs() <= contact.tolerance,
            "off the great circle at {p:?}: radius {radius}, height {height}"
        );
        if let Some(q) = previous {
            length += p.distance(q);
        }
        previous = Some(p);
    }
    let circle = core::f64::consts::TAU;
    assert!((length - circle).abs() < 1e-3, "once round: {length}");
}

#[test]
fn a_plane_resting_on_a_torus_comes_back_as_contact() {
    // A plane resting on top of a torus touches it along a whole circle.
    // The crossing marcher stalls there in fragments; `intersect_surfaces`
    // routes those fragments into the tangential trace and comes back with
    // the contact itself: one curve, closed, marked as touching rather than
    // crossing, on the circle of radius 3 at z = 1.
    let torus: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, 3.0, 1.0, T).unwrap()).into();
    let resting = pln(Point::new(0.0, 0.0, 1.0), Vector::Z);

    let found = ogeom_intersect::intersect_surfaces(
        &torus,
        &resting,
        ogeom_intersect::IntersectOptions {
            tolerance: 1e-4,
            marching: options(),
        },
        T,
    )
    .unwrap();
    let ogeom_intersect::SurfaceIntersection::Along(curves) = found else {
        panic!("the surfaces meet along the contact: {found:?}");
    };
    assert_eq!(curves.len(), 1, "one contact, not one per fragment");
    let contact = &curves[0];
    assert!(contact.tangential, "they touch here, they do not cross");
    assert!(contact.closed, "the contact circle closes");
    let domain = contact.curve.domain();
    for k in 0..=16 {
        let t = domain.0 + (domain.1 - domain.0) * f64::from(k) / 16.0;
        let p = contact.curve.point_at(t, T).unwrap();
        assert!(
            (p.x.hypot(p.y) - 3.0).abs() < 1e-3 && (p.z - 1.0).abs() < 1e-3,
            "the contact is the circle of radius 3 at z = 1: {p:?}"
        );
    }
}

#[test]
fn a_ball_seated_in_a_torus_tube_has_its_contact_walked() {
    // The tangency with no closed form: a unit ball centred on the tube's
    // own centre line touches the torus along a whole tube cross-section.
    // Nothing analytic answers this pair, so it goes through the marcher
    // (which stalls, as tangencies make it), and the stalled fragments seed
    // the tangential walker instead of being discarded.
    //
    // The contact comes back in two arcs rather than one loop, and that is
    // the sphere's parameterization talking: the tube circle runs through
    // both of the ball's chart poles, where its derivatives degenerate and
    // the walk has nothing to step along. Two arcs meeting at the poles
    // cover the circle; the number of pieces is a fact about the chart, the
    // curve they lie on is a fact about the surfaces.
    let torus: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, 3.0, 1.0, T).unwrap()).into();
    let seated = sph(Point::new(3.0, 0.0, 0.0), 1.0);

    let found = ogeom_intersect::intersect_surfaces(
        &torus,
        &seated,
        ogeom_intersect::IntersectOptions {
            tolerance: 1e-4,
            marching: options(),
        },
        T,
    )
    .unwrap();
    let ogeom_intersect::SurfaceIntersection::Along(curves) = found else {
        panic!("the ball meets the tube along the contact: {found:?}");
    };
    assert!(!curves.is_empty(), "the contact was not walked at all");
    assert!(
        curves.iter().all(|c| c.tangential && !c.exact),
        "every piece is marched contact, not a crossing"
    );
    let centre = Point::new(3.0, 0.0, 0.0);
    let mut length = 0.0;
    for c in &curves {
        let domain = c.curve.domain();
        let mut previous = None;
        for k in 0..=64 {
            let t = domain.0 + (domain.1 - domain.0) * f64::from(k) / 64.0;
            let p = c.curve.point_at(t, T).unwrap();
            assert!(
                (p.distance(centre) - 1.0).abs() < 1e-3 && p.y.abs() < 1e-3,
                "the contact is the tube circle at the ball's radius: {p:?}"
            );
            if let Some(q) = previous {
                length += p.distance(q);
            }
            previous = Some(p);
        }
    }
    let circle = 2.0 * core::f64::consts::PI;
    assert!(
        (length - circle).abs() < circle * 0.02,
        "the arcs together are the whole circle: {length}"
    );
}

#[test]
fn a_plane_through_a_torus_tube_cuts_two_loops() {
    // The transversal cousin of the tangent case: a plane through the tube
    // at half the minor radius cuts two closed loops, one around the outer
    // half, one around the inner.
    let torus: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, 3.0, 1.0, T).unwrap()).into();
    let cut = pln(Point::new(0.0, 0.0, 0.5), Vector::Z);

    let found = branches(&torus, &cut, options(), T).unwrap();
    assert_eq!(found.len(), 2);
    for br in &found {
        assert!(br.closed());
    }
    assert!(worst(&torus, &cut, &found) < 1e-7);
    // Distinct radii: one loop outside the tube's crown, one inside.
    let mut radii: Vec<f64> = found
        .iter()
        .map(|br| {
            let p = br.points[0];
            (p.x * p.x + p.y * p.y).sqrt()
        })
        .collect();
    radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert!(radii[0] < 3.0 && radii[1] > 3.0, "radii {radii:?}");
}

/// A plane through a torus's axis meets it in two whole tube circles, in
/// closed form, each starting on the outer equator where the torus's own
/// meridians do; a plane holding the torus's seam traces the seam itself.
#[test]
fn a_plane_through_a_torus_axis_cuts_two_meridians_exactly() {
    let torus: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, 3.0, 1.0, T).unwrap()).into();
    let seam_plane = pln(Point::ORIGIN, Vector::Y);
    let Meeting::Along(circles) = surface_surface(&torus, &seam_plane, T).unwrap() else {
        panic!("two circles");
    };
    assert_eq!(circles.len(), 2);
    let mut starts = Vec::new();
    for circle in &circles {
        let (a, b) = circle.domain();
        for k in 0..=16 {
            let p = circle
                .point_at(a + (b - a) * f64::from(k) / 16.0, T)
                .unwrap();
            let tube = (p.x.hypot(p.y) - 3.0).hypot(p.z);
            assert!((tube - 1.0).abs() < 1e-12 && p.y.abs() < 1e-12, "{p:?}");
        }
        starts.push(circle.point_at(a, T).unwrap());
    }
    // Each starts on the outer equator: four and minus four along x.
    assert!(
        starts
            .iter()
            .any(|p| p.distance(Point::new(4.0, 0.0, 0.0)) < 1e-12)
    );
    assert!(
        starts
            .iter()
            .any(|p| p.distance(Point::new(-4.0, 0.0, 0.0)) < 1e-12)
    );
}

#[test]
fn a_sphere_strictly_inside_another_is_apart_however_close() {
    let outer = sph(Point::ORIGIN, 2.0);
    let inner = sph(Point::new(0.999, 0.0, 0.0), 1.0);
    assert!(branches(&outer, &inner, options(), T).unwrap().is_empty());
    assert_eq!(surface_surface(&outer, &inner, T).unwrap(), Meeting::Apart);
}

/// Two planes through one point, a microradian apart, meet in a line. At
/// a hundred millimetres they stand a tenth of a micron apart, a thousand
/// times the confusion, so calling them one plane would be wrong.
#[test]
fn planes_a_microradian_apart_meet_in_a_line() {
    let a = 1e-6_f64;
    let plane = |n: Vector| -> SurfaceGeometry {
        PlaneSurface::new(Plane::through(Point::ORIGIN, Direction::new(n, T).unwrap())).into()
    };
    let flat = plane(Vector::new(0.0, 0.0, 1.0));
    let tilted = plane(Vector::new(0.0, a.sin(), a.cos()));
    assert!(matches!(
        surface_surface(&flat, &tilted, T).unwrap(),
        Meeting::Along(_)
    ));
    assert!(matches!(
        surface_surface(&flat, &flat, T).unwrap(),
        Meeting::Same
    ));
}

/// A line touching a sphere: the discriminant is zero, and rounding leaves
/// it a few ulps either side. Every scale and offset finds the touch.
#[test]
fn a_line_touching_a_sphere_is_found_at_every_scale() {
    use ogeom_geom::{Curve, LineCurve};
    use ogeom_intersect::{CurveSurfaceOptions, intersect_curve_surface};
    for radius in [1e-3, 0.7, 1.0, 7.3, 123.456, 1e3] {
        for shift in [0.0, 0.1, 1.0 / 3.0, 2.9] {
            let centre = Point::new(shift * radius, -shift, 0.3 * radius);
            let ball: SurfaceGeometry =
                SphereSurface::new(Sphere::new(frame(centre, Vector::Z), radius, T).unwrap())
                    .into();
            let touch = centre + Vector::new(0.0, radius, 0.0);
            let line: Curve = LineCurve::over(
                ogeom_math::Axis::new(touch - Vector::X * (3.0 * radius), Direction::X),
                0.0,
                6.0 * radius,
            )
            .unwrap()
            .into();
            let met =
                intersect_curve_surface(&line, &ball, CurveSurfaceOptions::default(), T).unwrap();
            assert!(!met.crossings.is_empty(), "radius {radius}, shift {shift}");
        }
    }
}

/// An oblique pipe tee: each closed section winds once round the branch
/// pipe, a loop in space whose chart image on the branch moves a period.
/// Its curve lies on both pipes within the tolerance it states, and that
/// tolerance is a fit's, not a pipe's size.
#[test]
fn an_oblique_pipe_tee_s_sections_lie_on_both_pipes() {
    use ogeom_intersect::{IntersectOptions, SurfaceIntersection, intersect_surfaces};
    let pipe = |origin: Point, axis: Vector, r: f64, reach: f64| -> SurfaceGeometry {
        CylinderSurface::new(
            Cylinder::new(frame(origin, axis), r, T).unwrap(),
            (-reach, reach),
        )
        .unwrap()
        .into()
    };
    for (ra, rb, tilt, offset) in [
        (10.0, 5.0, 1.0_f64, 0.0),
        (10.0, 7.0, 0.7, 1.0),
        (25.0, 12.5, 1.0, 0.0),
        (5.0, 3.0, 0.7, 0.0),
    ] {
        let centre = Point::new(100.0, 200.0, 300.0);
        let a = pipe(centre, Vector::Z, ra, 4.0 * ra);
        let b = pipe(
            centre + Vector::new(0.0, offset, 0.0),
            Vector::new(tilt.cos(), 0.0, tilt.sin()),
            rb,
            4.0 * ra,
        );
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&a, &b, IntersectOptions::default(), T).unwrap()
        else {
            panic!("the pipes cross");
        };
        assert_eq!(curves.len(), 2, "{ra} {rb} {tilt} {offset}");
        for c in &curves {
            assert!(
                c.tolerance < 1e-3,
                "{ra} {rb} {tilt} {offset}: {}",
                c.tolerance
            );
            let (lo, hi) = c.curve.domain();
            for k in 0..=1000 {
                let p = c
                    .curve
                    .point_at(lo + (hi - lo) * f64::from(k) / 1000.0, T)
                    .unwrap();
                let miss = off(&a, p).abs().max(off(&b, p).abs());
                assert!(
                    miss <= c.tolerance,
                    "{ra} {rb} {tilt} {offset}: {miss} off, stated {}",
                    c.tolerance
                );
            }
        }
    }
}

/// A spindle torus (its tube swallows the axis) meets a plane square to
/// its axis, a coaxial cylinder and a coaxial ring torus also on its folded
/// half past the axis: every parallel comes back, on both surfaces, with
/// its image in each chart.
#[test]
fn a_spindle_torus_s_folded_half_is_sectioned() {
    use ogeom_intersect::{IntersectOptions, SurfaceIntersection, intersect_surfaces};
    let spindle: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, 10.0, 15.0, T).unwrap()).into();
    let ring: SurfaceGeometry =
        TorusSurface::new(Torus::new(Frame::WORLD, 4.0, 3.0, T).unwrap()).into();
    let drum: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, 3.0, T).unwrap(), (-50.0, 50.0))
            .unwrap()
            .into();
    let level: SurfaceGeometry = PlaneSurface::over(
        Plane::through(Point::new(0.0, 0.0, 5.0), Direction::Z),
        (-50.0, 50.0),
        (-50.0, 50.0),
    )
    .unwrap()
    .into();
    for (other, count) in [(level, 2), (drum, 4), (ring, 2)] {
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&spindle, &other, IntersectOptions::default(), T).unwrap()
        else {
            panic!("the surfaces meet");
        };
        assert_eq!(curves.len(), count);
        for c in &curves {
            assert!(c.on_a.is_some() && c.on_b.is_some());
            let (lo, hi) = c.curve.domain();
            for k in 0..=16 {
                let p = c
                    .curve
                    .point_at(lo + (hi - lo) * f64::from(k) / 16.0, T)
                    .unwrap();
                assert!(off(&spindle, p).abs() < 1e-9 && off(&other, p).abs() < 1e-9);
            }
        }
    }
}

/// A cone whose window runs through its apex meets a coaxial cylinder, a
/// parallel-sided coaxial cone and a level plane also on its far nappe,
/// exactly and with its image in each chart; one stopping short of the
/// apex keeps only its own nappe's section.
#[test]
fn a_cone_through_its_apex_is_sectioned_on_both_nappes() {
    use ogeom_intersect::{IntersectOptions, SurfaceIntersection, intersect_surfaces};
    let quarter = core::f64::consts::FRAC_PI_4;
    let cone = |r: f64, window: (f64, f64)| -> SurfaceGeometry {
        ConeSurface::new(Cone::new(Frame::WORLD, r, quarter, T).unwrap(), window)
            .unwrap()
            .into()
    };
    let drum: SurfaceGeometry =
        CylinderSurface::new(Cylinder::new(Frame::WORLD, 10.0, T).unwrap(), (-60.0, 60.0))
            .unwrap()
            .into();
    let level: SurfaceGeometry = PlaneSurface::over(
        Plane::through(Point::new(0.0, 0.0, -30.0), Direction::Z),
        (-50.0, 50.0),
        (-50.0, 50.0),
    )
    .unwrap()
    .into();
    // The apex stands at z = -20.
    let through = cone(20.0, (-50.0, 50.0));
    let short = cone(20.0, (-15.0, 50.0));
    for (a, b, heights) in [
        (&through, &drum, vec![-10.0, -30.0]),
        (&short, &drum, vec![-10.0]),
        (&through, &cone(40.0, (-100.0, 50.0)), vec![-30.0]),
        (&through, &level, vec![-30.0]),
    ] {
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(a, b, IntersectOptions::default(), T).unwrap()
        else {
            panic!("the surfaces meet");
        };
        let mut found: Vec<f64> = curves
            .iter()
            .map(|c| {
                assert!(c.exact && c.on_a.is_some() && c.on_b.is_some());
                let (lo, _) = c.curve.domain();
                let p = c.curve.point_at(lo, T).unwrap();
                assert!(off(a, p).abs() < 1e-9 && off(b, p).abs() < 1e-9);
                p.z
            })
            .collect();
        found.sort_by(|x, y| y.total_cmp(x));
        assert_eq!(found.len(), heights.len(), "{found:?}");
        for (got, want) in found.iter().zip(&heights) {
            assert!((got - want).abs() < 1e-9, "{found:?}");
        }
    }
}

/// A steep plane passing half a millimetre from a cone's apex: the fitted
/// section's images, lifted through their surfaces, stand within the
/// tolerance it states, however fast the cone's chart turns there.
#[test]
fn a_section_beside_a_cone_s_apex_states_its_true_tolerance() {
    use ogeom_geom::{Curve2d as _, Surface as _};
    use ogeom_intersect::{IntersectOptions, SurfaceIntersection, intersect_surfaces};
    let cone: SurfaceGeometry = ConeSurface::new(
        Cone::new(Frame::WORLD, 10.0, core::f64::consts::FRAC_PI_4, T).unwrap(),
        (-9.5, 50.0),
    )
    .unwrap()
    .into();
    let tilt: f64 = 1.2;
    let steep: SurfaceGeometry = PlaneSurface::over(
        Plane::through(
            Point::new(0.5, 0.3, -9.0),
            Direction::new(Vector::new(tilt.sin(), 0.0, tilt.cos()), T).unwrap(),
        ),
        (-200.0, 200.0),
        (-200.0, 200.0),
    )
    .unwrap()
    .into();
    let SurfaceIntersection::Along(curves) =
        intersect_surfaces(&cone, &steep, IntersectOptions::default(), T).unwrap()
    else {
        panic!("the plane cuts the cone");
    };
    for c in &curves {
        let (lo, hi) = c.curve.domain();
        for k in 0..=4000 {
            let t = lo + (hi - lo) * f64::from(k) / 4000.0;
            let on = c.curve.point_at(t, T).unwrap();
            for (surface, image) in [(&cone, &c.on_a), (&steep, &c.on_b)] {
                let at = image.as_ref().unwrap().point_at(t, T).unwrap();
                let lifted = surface.point_at(at.x, at.y, T).unwrap();
                assert!(
                    lifted.distance(on) <= c.tolerance,
                    "{} off at {t}, stated {}",
                    lifted.distance(on),
                    c.tolerance
                );
            }
        }
    }
}

/// A sphere centred on a torus's axis meets it in parallels: two where it
/// crosses the tube, one where it sits against the tube's inside, each an
/// exact circle and the seated one a single tangential contact.
#[test]
fn a_sphere_on_a_torus_s_axis_meets_it_in_parallels() {
    use ogeom_intersect::{IntersectOptions, SurfaceIntersection, intersect_surfaces};
    let centre = Point::new(100.0, 200.0, 300.0);
    let torus: SurfaceGeometry = TorusSurface::new(
        Torus::new(
            Frame::new(centre, Direction::Z, Direction::X, T).unwrap(),
            30.0,
            5.0,
            T,
        )
        .unwrap(),
    )
    .into();
    for (radius, count, tangential) in [(25.0, 1, true), (27.0, 2, false), (35.0, 1, true)] {
        let ball = sph(centre, radius);
        let SurfaceIntersection::Along(curves) =
            intersect_surfaces(&torus, &ball, IntersectOptions::default(), T).unwrap()
        else {
            panic!("the sphere meets the torus at {radius}");
        };
        assert_eq!(curves.len(), count, "{radius}");
        for c in &curves {
            assert!(c.exact && c.tangential == tangential, "{radius}");
            assert!(c.on_a.is_some() && c.on_b.is_some());
            let (lo, hi) = c.curve.domain();
            for k in 0..=16 {
                let p = c
                    .curve
                    .point_at(lo + (hi - lo) * f64::from(k) / 16.0, T)
                    .unwrap();
                assert!(off(&torus, p).abs() < 1e-9 && off(&ball, p).abs() < 1e-9);
            }
        }
    }
}
