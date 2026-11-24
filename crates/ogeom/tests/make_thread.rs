//! A thread cut on a cylindrical face, and a face's exact extent along an
//! axis. Each measure is held to a closed form: the volume a groove of
//! constant section takes (its section's first moment about the axis, times
//! a turn, per pitch of length), the pitch between crossings of one flank
//! along a generator, the flank's lean, the radii of the groove's floor and
//! of the crest, and the hand.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use core::f64::consts::{FRAC_PI_6, PI, TAU};

use ogeom::algo::{
    Containment, check, classify_in_solid_exact, extent_along, make_box, make_cone, make_cylinder,
    tight_bounds, transformed, volume_properties,
};
use ogeom::core::Tolerances;
use ogeom::geom::SurfaceGeometry;
use ogeom::math::{Axis, Direction, Frame, Point, Transform, Vector};
use ogeom::mesh::Deflection;
use ogeom::offset::{ThreadForm, make_thread};
use ogeom::topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The volumes compared are measured on the exact faces at a chord of a
/// micrometre, which holds a few millimetres cubed to a few parts in 1e8;
/// 1e-6 relative leaves that margin for the fitted thread walls.
const VOLUME: f64 = 1e-6;

/// The thread's walls are fitted through their exact screw images, within
/// the sweep's fitting tolerance of a hundredth of a micrometre; the radii
/// measured on them have come within 3e-7 of the form.
const FITTED: f64 = 1e-6;

fn volume(model: &Model, shape: &Shape) -> f64 {
    volume_properties(model, shape, Deflection::with_chord(1e-3).unwrap(), T)
        .unwrap()
        .mass
}

fn z_axis() -> Axis {
    Axis {
        location: Point::ORIGIN,
        direction: Direction::Z,
    }
}

/// The face of `shape` on a cylinder of `radius`.
fn cylindrical_face(model: &Model, shape: &Shape, radius: f64) -> Shape {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| on_cylinder(model, f, radius))
        .expect("a cylindrical face")
}

fn on_cylinder(model: &Model, face: &Shape, radius: f64) -> bool {
    let data = model.node(face).unwrap().data().as_face().unwrap();
    let mut surface = model.geometry().surface(data.surface).unwrap();
    while let SurfaceGeometry::Trimmed(t) = surface {
        surface = t.basis();
    }
    matches!(
        surface,
        SurfaceGeometry::Cylinder(c) if (c.cylinder().radius() - radius).abs() < 1e-9
    )
}

/// `int_0^h u^k (a + sigma u)(root + 2 u tan) du` for k = 0 and 1: the
/// groove's width along the axis times the radius, over the height `u`
/// above its floor at radius `a`, `sigma` the way the radius grows with
/// `u`.
fn moments(form: &ThreadForm, a: f64, sigma: f64) -> (f64, f64) {
    let (h, w, t2) = (form.depth, form.root, 2.0 * form.flank_angle.tan());
    let g0 = a * w * h + (a * t2 + sigma * w) * h * h / 2.0 + sigma * t2 * h.powi(3) / 3.0;
    let g1 =
        a * w * h * h / 2.0 + (a * t2 + sigma * w) * h.powi(3) / 3.0 + sigma * t2 * h.powi(4) / 4.0;
    (g0, g1)
}

/// The volume a thread of `form` over `length` takes from a face of
/// `radius`, its run-out included. Per unit length a full groove takes
/// `2 pi / pitch` times its section's first moment about the axis; over
/// the run-out the depth falls evenly, so each height `u` above the floor
/// is cut over the share `u / depth` of it.
fn removed(form: &ThreadForm, radius: f64, length: f64, external: bool) -> f64 {
    let (a, sigma) = if external {
        (radius - form.depth, 1.0)
    } else {
        (radius + form.depth, -1.0)
    };
    let (g0, g1) = moments(form, a, sigma);
    TAU / form.pitch * ((length - form.run_out) * g0 + form.run_out / form.depth * g1)
}

/// A point `rho` from the Z axis, at angle `theta`, height `z`.
fn at(rho: f64, theta: f64, z: f64) -> Point {
    Point::new(rho * theta.cos(), rho * theta.sin(), z)
}

fn inside(model: &Model, solid: &Shape, p: Point) -> bool {
    classify_in_solid_exact(model, solid, p, T).unwrap() == Containment::In
}

/// Where along `f` on `[a, b]` the answer turns from `f(a)`, to 1e-10.
fn bisect(mut a: f64, mut b: f64, f: impl Fn(f64) -> bool) -> f64 {
    let start = f(a);
    while (b - a).abs() > 1e-10 {
        let mid = f64::midpoint(a, b);
        if f(mid) == start {
            a = mid;
        } else {
            b = mid;
        }
    }
    f64::midpoint(a, b)
}

/// Along the generator `rho` out at `theta`, the heights in `[z0, z1]`
/// where the material starts going up: a groove's upper flank. The steps
/// of a twentieth of a millimetre are finer than any groove or land
/// measured here.
fn entries(model: &Model, solid: &Shape, rho: f64, theta: f64, z0: f64, z1: f64) -> Vec<f64> {
    let is_in = |z: f64| inside(model, solid, at(rho, theta, z));
    #[allow(clippy::cast_possible_truncation, reason = "a small positive count")]
    let steps = ((z1 - z0) / 0.05).ceil() as i32;
    let dz = (z1 - z0) / f64::from(steps);
    let mut found = Vec::new();
    let mut was = is_in(z0);
    for k in 1..=steps {
        let z = z0 + dz * f64::from(k);
        let now = is_in(z);
        if now && !was {
            found.push(bisect(z - dz, z, is_in));
        }
        was = now;
    }
    found
}

/// The faces of `result` that none of the solid's own surfaces carry: the
/// thread's walls.
fn walls(model: &Model, result: &Shape) -> Vec<Shape> {
    explore_unique(model, result, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            let mut s = model.geometry().surface(data.surface).unwrap();
            while let SurfaceGeometry::Trimmed(t) = s {
                s = t.basis();
            }
            !matches!(s, SurfaceGeometry::Plane(_) | SurfaceGeometry::Cylinder(_))
        })
        .collect()
}

/// The interval of `faces` along Z.
fn span_of(model: &Model, faces: &[Shape]) -> (f64, f64) {
    faces
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), f| {
            let (a, b) = extent_along(model, f, z_axis(), T).unwrap();
            (lo.min(a), hi.max(b))
        })
}

/// What one generator of an external thread on a face of `radius` shows:
/// the upper flank recurring every `pitch`, leaning at the form's angle,
/// the floor `depth` in and the crest on the face, and the groove moving
/// up the axis as the angle grows for a right hand, down for a left.
fn measures_its_form(
    model: &Model,
    threaded: &Shape,
    form: &ThreadForm,
    radius: f64,
    z0: f64,
    z1: f64,
) {
    let (h, p) = (form.depth, form.pitch);
    let theta = 0.3;
    let mid = entries(model, threaded, radius - h / 2.0, theta, z0, z1);
    assert!(mid.len() >= 2, "{mid:?}");
    for pair in mid.windows(2) {
        assert!(
            (pair[1] - pair[0] - p).abs() < FITTED,
            "flank crossings {pair:?} against a pitch of {p}"
        );
    }
    // The flank's lean: where it crosses a quarter and three quarters of
    // the way down.
    let (r1, r2) = (radius - h / 4.0, radius - 3.0 * h / 4.0);
    let near = |rho: f64| {
        let found = entries(
            model,
            threaded,
            rho,
            theta,
            mid[0] - p / 2.0,
            mid[0] + p / 2.0,
        );
        assert_eq!(found.len(), 1, "{rho}: {found:?}");
        found[0]
    };
    let lean = ((near(r1) - near(r2)) / (r1 - r2)).atan();
    assert!(
        (lean - form.flank_angle).abs() < 1e-6,
        "a flank leans {lean} against {}",
        form.flank_angle
    );
    // The groove's middle, half a groove's floor width below the flank at
    // its floor, and the crest's, half a pitch on.
    let floor_flank = mid[0] - h / 2.0 * form.flank_angle.tan();
    let groove = floor_flank - form.root / 2.0;
    let floor = bisect(radius - h - 0.05, radius - h + 0.05, |rho| {
        inside(model, threaded, at(rho, theta, groove))
    });
    assert!(
        (floor - (radius - h)).abs() < FITTED,
        "the floor stands {floor} out, against {}",
        radius - h
    );
    let crest = groove + p / 2.0;
    let face = bisect(radius - 0.05, radius + 0.05, |rho| {
        inside(model, threaded, at(rho, theta, crest))
    });
    assert!(
        (face - radius).abs() < FITTED,
        "the crest stands {face} out"
    );
    // A tenth of a turn on, the flank has moved a tenth of the lead.
    let lead = p * f64::from(form.starts);
    let turned = entries(
        model,
        threaded,
        radius - h / 2.0,
        theta + TAU / 10.0,
        mid[0] - p,
        mid[0] + p,
    );
    let sense = if form.left_handed { -1.0 } else { 1.0 };
    let want = mid[0] + sense * lead / 10.0;
    assert!(
        turned
            .iter()
            .any(|z| ((z - want) / p).round().mul_add(-p, z - want).abs() < FITTED),
        "{turned:?}: no crossing a pitch from {want}"
    );
}

/// M8 x 1.25 external over a whole 6 mm shaft: valid, the closed form's
/// volume, the groove running out of both ends, and the form on a
/// generator.
#[test]
fn an_m8_thread_over_a_whole_shaft_measures_its_closed_form() {
    let mut model = Model::new();
    let shaft = make_cylinder(&mut model, Frame::WORLD, 4.0, 6.0, T)
        .unwrap()
        .shape;
    let side = cylindrical_face(&model, &shaft, 4.0);
    let form = ThreadForm::metric(1.25);
    let threaded = make_thread(&mut model, &shaft, &side, &form, T)
        .unwrap()
        .shape;
    let diagnosis = check(&model, &threaded, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let want = removed(&form, 4.0, 6.0, true);
    let got = PI * 16.0 * 6.0 - volume(&model, &threaded);
    assert!(
        (got - want).abs() < want * VOLUME,
        "removed {got} against {want}"
    );
    let (lo, hi) = span_of(&model, &walls(&model, &threaded));
    assert!(
        lo.abs() < 1e-9 && (hi - 6.0).abs() < 1e-9,
        "walls over {lo}..{hi}"
    );
    measures_its_form(&model, &threaded, &form, 4.0, 1.0, 3.6);
}

/// M8 x 1.25 external, 15 mm up a 20 mm shaft from its base with a 2 mm
/// run-out: valid; the closed form's volume; the walls exactly 15 mm along
/// the axis; the form on a generator; and over the run-out the floor
/// rising linearly to the face at 15 mm.
#[test]
#[ignore = "heavy"]
fn an_m8_thread_with_a_run_out_measures_its_closed_form() {
    let mut model = Model::new();
    let shaft = make_cylinder(&mut model, Frame::WORLD, 4.0, 20.0, T)
        .unwrap()
        .shape;
    let side = cylindrical_face(&model, &shaft, 4.0);
    let form = ThreadForm {
        length: Some(15.0),
        run_out: 2.0,
        ..ThreadForm::metric(1.25)
    };
    let threaded = make_thread(&mut model, &shaft, &side, &form, T)
        .unwrap()
        .shape;
    let diagnosis = check(&model, &threaded, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let want = removed(&form, 4.0, 15.0, true);
    let got = PI * 16.0 * 20.0 - volume(&model, &threaded);
    assert!(
        (got - want).abs() < want * VOLUME,
        "removed {got} against {want}"
    );
    // The run-out's cone ends on the face along the section the boolean
    // traces between the two, which it holds to the fitted walls' measure.
    let (lo, hi) = span_of(&model, &walls(&model, &threaded));
    assert!(
        lo.abs() < 1e-9 && (hi - 15.0).abs() < FITTED,
        "walls over {lo}..{hi}"
    );
    measures_its_form(&model, &threaded, &form, 4.0, 1.0, 10.0);
    // The groove's middle on this generator, each pitch up the run-out:
    // its floor where the cone from 4 - depth at 13 to 4 at 15 says.
    let (h, p) = (form.depth, form.pitch);
    let theta = 0.3;
    let first = entries(&model, &threaded, 4.0 - h / 2.0, theta, 1.0, 1.0 + p)[0];
    let groove = first - h / 2.0 * FRAC_PI_6.tan() - form.root / 2.0;
    let mut seen = 0;
    for k in 0..20 {
        let z = groove + p * f64::from(k);
        if !(13.05..14.95).contains(&z) {
            continue;
        }
        let floor = bisect(4.0 - h - 0.05, 4.0, |rho| {
            inside(&model, &threaded, at(rho, theta, z))
        });
        let want = 4.0 - h * (15.0 - z) / 2.0;
        assert!(
            (floor - want).abs() < FITTED,
            "at {z} the floor is {floor}, against {want}"
        );
        seen += 1;
    }
    assert!(seen >= 1);
}

/// A left-hand thread of two starts: the same volume per length as one
/// start (two grooves, each over twice the lead), and its grooves moving
/// down the axis as the angle grows.
#[test]
fn a_two_start_left_hand_thread_measures_its_closed_form() {
    let mut model = Model::new();
    let shaft = make_cylinder(&mut model, Frame::WORLD, 5.0, 3.0, T)
        .unwrap()
        .shape;
    let side = cylindrical_face(&model, &shaft, 5.0);
    let form = ThreadForm {
        left_handed: true,
        starts: 2,
        ..ThreadForm::metric(1.5)
    };
    let threaded = make_thread(&mut model, &shaft, &side, &form, T)
        .unwrap()
        .shape;
    let diagnosis = check(&model, &threaded, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let want = removed(&form, 5.0, 3.0, true);
    let got = PI * 25.0 * 3.0 - volume(&model, &threaded);
    assert!(
        (got - want).abs() < want * VOLUME,
        "removed {got} against {want}"
    );
    measures_its_form(&model, &threaded, &form, 5.0, 0.2, 2.8);
}

/// An internal M8 thread in a 6.8 mm bore through a 10 mm plate, cut
/// 0.6 deep to the major diameter: valid, the closed form's volume, the
/// plate's box unchanged and the bore clear.
#[test]
fn an_internal_thread_stays_inside_the_plate() {
    let mut model = Model::new();
    let plate = make_box(
        &mut model,
        Frame::new(Point::new(-10.0, -10.0, 0.0), Direction::Z, Direction::X, T).unwrap(),
        (20.0, 20.0, 10.0),
        T,
    )
    .unwrap()
    .shape;
    let hole = make_cylinder(
        &mut model,
        Frame::about(Point::new(0.0, 0.0, -1.0), Direction::Z),
        3.4,
        12.0,
        T,
    )
    .unwrap()
    .shape;
    let bored = ogeom::boolean::cut(&mut model, &plate, &hole, T)
        .unwrap()
        .shape;
    let bore = cylindrical_face(&model, &bored, 3.4);
    let form = ThreadForm::metric_internal(1.25).with_depth(0.6);
    let threaded = make_thread(&mut model, &bored, &bore, &form, T)
        .unwrap()
        .shape;
    let diagnosis = check(&model, &threaded, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let whole = 4000.0 - PI * 3.4 * 3.4 * 10.0;
    let want = removed(&form, 3.4, 10.0, false);
    let got = whole - volume(&model, &threaded);
    assert!(
        (got - want).abs() < want * VOLUME,
        "removed {got} against {want}"
    );
    let bound = tight_bounds(&model, &threaded, T).unwrap();
    let (lo, hi) = (bound.low().unwrap(), bound.high().unwrap());
    assert!(
        (lo - Point::new(-10.0, -10.0, 0.0)).magnitude() < 1e-9
            && (hi - Point::new(10.0, 10.0, 10.0)).magnitude() < 1e-9,
        "{lo:?} {hi:?}"
    );
    for k in 0..40 {
        let (theta, z) = (f64::from(k) * 0.7, 0.1 + f64::from(k) * 0.245);
        assert!(
            !inside(&model, &threaded, at(3.399, theta, z)),
            "material in the bore"
        );
        assert!(
            inside(&model, &threaded, at(4.001, theta, z)),
            "the groove reached past 4"
        );
    }
    // Along a generator the material starts again every pitch, a groove
    // floor's width below the next groove.
    let ups = entries(&model, &threaded, 3.7, 0.3, 1.0, 9.0);
    for pair in ups.windows(2) {
        assert!((pair[1] - pair[0] - 1.25).abs() < FITTED, "{pair:?}");
    }
}

/// A 45 degree lead-in off the end the thread starts from: the thread's
/// volume and the chamfer's cone, less what the two share (each height
/// `u` above the floor the cone takes over `u` of the axis).
#[test]
fn a_chamfered_thread_measures_its_closed_form() {
    let mut model = Model::new();
    let shaft = make_cylinder(&mut model, Frame::WORLD, 4.0, 5.0, T)
        .unwrap()
        .shape;
    let side = cylindrical_face(&model, &shaft, 4.0);
    let form = ThreadForm {
        chamfer: true,
        from_far_end: true,
        ..ThreadForm::metric(1.25)
    };
    let threaded = make_thread(&mut model, &shaft, &side, &form, T)
        .unwrap()
        .shape;
    let diagnosis = check(&model, &threaded, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let (r, h) = (4.0_f64, form.depth);
    let cone = PI * (r * r * h - (r.powi(3) - (r - h).powi(3)) / 3.0);
    let (_, g1) = moments(&form, r - h, 1.0);
    let want = cone + removed(&form, r, 5.0, true) - TAU / form.pitch * g1;
    let got = PI * 16.0 * 5.0 - volume(&model, &threaded);
    assert!(
        (got - want).abs() < want * VOLUME,
        "removed {got} against {want}"
    );
    // The chamfer is at the top, where the thread starts.
    assert!(!inside(&model, &threaded, at(3.99, 1.0, 4.99)));
}

/// A shaft standing on a wider shoulder: a thread from the shoulder's end
/// of the face ends square in the shoulder's plane, and a chamfer there
/// is refused, since the solid carries on past that end.
#[test]
fn a_thread_against_a_shoulder_ends_in_its_plane() {
    let mut model = Model::new();
    let shoulder = make_cylinder(&mut model, Frame::WORLD, 6.0, 2.0, T)
        .unwrap()
        .shape;
    let pin = make_cylinder(
        &mut model,
        Frame::about(Point::new(0.0, 0.0, 2.0), Direction::Z),
        4.0,
        4.0,
        T,
    )
    .unwrap()
    .shape;
    let stepped = ogeom::boolean::fuse(&mut model, &shoulder, &pin, T)
        .unwrap()
        .shape;
    let side = cylindrical_face(&model, &stepped, 4.0);
    let form = ThreadForm {
        length: Some(3.0),
        ..ThreadForm::metric(1.25)
    };
    let threaded = make_thread(&mut model, &stepped, &side, &form, T)
        .unwrap()
        .shape;
    let diagnosis = check(&model, &threaded, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let whole = PI * (36.0 * 2.0 + 16.0 * 4.0);
    let want = removed(&form, 4.0, 3.0, true);
    let got = whole - volume(&model, &threaded);
    assert!(
        (got - want).abs() < want * VOLUME,
        "removed {got} against {want}"
    );
    let (lo, hi) = span_of(&model, &walls(&model, &threaded));
    assert!(
        (lo - 2.0).abs() < FITTED && (hi - 5.0).abs() < 1e-9,
        "walls over {lo}..{hi}"
    );
    // The shoulder is untouched.
    assert!(inside(&model, &threaded, at(5.9, 0.4, 1.99)));
    let chamfered = ThreadForm {
        chamfer: true,
        ..form
    };
    assert!(make_thread(&mut model, &stepped, &side, &chamfered, T).is_err());
}

/// A thread run from the pin's free top down its whole length to the
/// shoulder, fading out over the last 1.5 mm: its run-out cone meets the
/// face in the shoulder's plane and takes nothing of the shoulder, so the
/// volume removed is the closed form's.
#[test]
fn a_run_out_against_a_shoulder_leaves_the_shoulder_whole() {
    let mut model = Model::new();
    let shoulder = make_cylinder(&mut model, Frame::WORLD, 6.0, 2.0, T)
        .unwrap()
        .shape;
    let pin = make_cylinder(
        &mut model,
        Frame::about(Point::new(0.0, 0.0, 2.0), Direction::Z),
        4.0,
        4.0,
        T,
    )
    .unwrap()
    .shape;
    let stepped = ogeom::boolean::fuse(&mut model, &shoulder, &pin, T)
        .unwrap()
        .shape;
    let side = cylindrical_face(&model, &stepped, 4.0);
    let form = ThreadForm {
        from_far_end: true,
        run_out: 1.5,
        ..ThreadForm::metric(1.25)
    };
    let threaded = make_thread(&mut model, &stepped, &side, &form, T)
        .unwrap()
        .shape;
    let diagnosis = check(&model, &threaded, T).unwrap();
    assert!(diagnosis.is_valid(), "{diagnosis}");
    let whole = PI * (36.0 * 2.0 + 16.0 * 4.0);
    let want = removed(&form, 4.0, 4.0, true);
    let got = whole - volume(&model, &threaded);
    assert!(
        (got - want).abs() < want * VOLUME,
        "removed {got} against {want}"
    );
}

/// What is refused: a form that does not close one pitch, a thread longer
/// than its face, a run-out longer than the thread, a face that is not a
/// cylinder's, and a face of another solid.
#[test]
fn a_thread_that_cannot_be_cut_is_refused() {
    let mut model = Model::new();
    let shaft = make_cylinder(&mut model, Frame::WORLD, 4.0, 6.0, T)
        .unwrap()
        .shape;
    let side = cylindrical_face(&model, &shaft, 4.0);
    let metric = ThreadForm::metric(1.25);
    for form in [
        ThreadForm {
            crest: 0.3,
            ..metric
        },
        ThreadForm {
            length: Some(6.5),
            ..metric
        },
        ThreadForm {
            length: Some(3.0),
            run_out: 3.5,
            ..metric
        },
        ThreadForm {
            starts: 0,
            ..metric
        },
        ThreadForm {
            flank_angle: 0.0,
            depth: 4.5,
            crest: 0.625,
            root: 0.625,
            ..metric
        },
    ] {
        assert!(
            make_thread(&mut model, &shaft, &side, &form, T).is_err(),
            "{form:?}"
        );
    }
    let flat = explore_unique(&model, &shaft, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| !on_cylinder(&model, f, 4.0))
        .unwrap();
    assert!(make_thread(&mut model, &shaft, &flat, &metric, T).is_err());
    let other = make_cylinder(&mut model, Frame::WORLD, 4.0, 6.0, T)
        .unwrap()
        .shape;
    let foreign = cylindrical_face(&model, &other, 4.0);
    assert!(make_thread(&mut model, &shaft, &foreign, &metric, T).is_err());
}

/// A cylinder's side runs its height along its axis; tilted by `a` from
/// the axis it covers its height times `cos a` plus its radius times
/// `sin a` at each end.
#[test]
fn a_cylinder_face_extends_its_height_along_its_axis() {
    let mut model = Model::new();
    let base = Point::new(1.0, -2.0, 3.0);
    let shaft = make_cylinder(&mut model, Frame::about(base, Direction::Z), 4.0, 20.0, T)
        .unwrap()
        .shape;
    let side = cylindrical_face(&model, &shaft, 4.0);
    let (lo, hi) = extent_along(&model, &side, z_axis(), T).unwrap();
    assert!(
        (lo - 3.0).abs() < 1e-12 && (hi - 23.0).abs() < 1e-12,
        "{lo} {hi}"
    );

    let tilt: f64 = 0.4;
    let dir = Direction::new(Vector::new(tilt.sin(), 0.0, tilt.cos()), T).unwrap();
    let tilted = make_cylinder(&mut model, Frame::about(base, dir), 4.0, 20.0, T)
        .unwrap()
        .shape;
    let side = cylindrical_face(&model, &tilted, 4.0);
    let (lo, hi) = extent_along(&model, &side, z_axis(), T).unwrap();
    let (want_lo, want_hi) = (
        3.0 - 4.0 * tilt.sin(),
        3.0 + 20.0 * tilt.cos() + 4.0 * tilt.sin(),
    );
    assert!(
        (lo - want_lo).abs() < 1e-9 && (hi - want_hi).abs() < 1e-9,
        "{lo} {hi} against {want_lo} {want_hi}"
    );
    // Along an axis that is not through the origin and points down, the
    // interval is the same one read from the other end.
    let down = Axis {
        location: Point::new(5.0, 7.0, 10.0),
        direction: -Direction::Z,
    };
    let (dlo, dhi) = extent_along(&model, &side, down, T).unwrap();
    assert!(
        (dlo - (10.0 - hi)).abs() < 1e-9 && (dhi - (10.0 - lo)).abs() < 1e-9,
        "{dlo} {dhi}"
    );
}

/// A cone's side tilted by `a`: each end circle covers its centre's
/// height plus or minus its radius times `sin a`, and the side holds both.
#[test]
fn a_cone_face_extends_between_its_end_circles() {
    let mut model = Model::new();
    let tilt: f64 = 0.7;
    let dir = Direction::new(Vector::new(0.0, tilt.sin(), tilt.cos()), T).unwrap();
    let cone = make_cone(
        &mut model,
        Frame::about(Point::ORIGIN, dir),
        3.0,
        1.0,
        5.0,
        T,
    )
    .unwrap()
    .shape;
    let side = explore_unique(&model, &cone, ShapeType::Face)
        .unwrap()
        .into_iter()
        .find(|f| {
            let data = model.node(f).unwrap().data().as_face().unwrap();
            matches!(
                model.geometry().surface(data.surface).unwrap(),
                SurfaceGeometry::Cone(_)
            )
        })
        .unwrap();
    let (lo, hi) = extent_along(&model, &side, z_axis(), T).unwrap();
    let top = 5.0 * tilt.cos();
    let want_lo = (-3.0 * tilt.sin()).min(top - tilt.sin());
    let want_hi = (3.0 * tilt.sin()).max(top + tilt.sin());
    assert!(
        (lo - want_lo).abs() < 1e-9 && (hi - want_hi).abs() < 1e-9,
        "{lo} {hi} against {want_lo} {want_hi}"
    );
}

/// A cylinder of radius 2 cut by the plane `z = 6 + x tan b`: its side
/// runs from 0 up to `6 + 2 tan b`, the top of the ellipse the plane
/// leaves on it.
#[test]
fn a_slanted_cut_cylinder_extends_to_its_ellipse() {
    let mut model = Model::new();
    let rod = make_cylinder(&mut model, Frame::WORLD, 2.0, 10.0, T)
        .unwrap()
        .shape;
    let slant: f64 = 0.5;
    let block = make_box(&mut model, Frame::WORLD, (20.0, 20.0, 20.0), T)
        .unwrap()
        .shape;
    let lift = Transform::translation(Vector::new(-10.0, -10.0, 0.0));
    let turn = Transform::rotation(
        Axis {
            location: Point::ORIGIN,
            direction: Direction::Y,
        },
        -slant,
    );
    let to = Transform::translation(Vector::new(0.0, 0.0, 6.0));
    let cutter = transformed(&mut model, &block, to * turn * lift)
        .unwrap()
        .shape;
    let cut = ogeom::boolean::cut(&mut model, &rod, &cutter, T)
        .unwrap()
        .shape;
    let side = cylindrical_face(&model, &cut, 2.0);
    let (lo, hi) = extent_along(&model, &side, z_axis(), T).unwrap();
    let want = 6.0 + 2.0 * slant.tan();
    assert!(
        lo.abs() < 1e-9 && (hi - want).abs() < 1e-9,
        "{lo} {hi} against {want}"
    );
    // The plane crosses the axis at 6, so the rod keeps 6 of its height on
    // average: the cut is the one asked for.
    let kept = volume(&model, &cut);
    assert!((kept - PI * 4.0 * 6.0).abs() < 1e-6, "{kept}");
}
