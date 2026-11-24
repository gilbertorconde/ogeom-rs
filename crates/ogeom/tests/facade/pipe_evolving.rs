//! Pipe shells whose section scales and twists along the spine, measured
//! against closed forms: a linear scale down a straight spine is a
//! frustum or a pyramid, exactly; a twist keeps every slice the profile
//! turned, so the volume is the slices' areas summed along the spine,
//! whatever law the scale and the twist follow.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom::algo::{check, make_edge, make_face, make_polygon, make_wire, volume_properties};
use ogeom::core::Tolerances;
use ogeom::geom::{CircleCurve, Curve, Curve3d as _, LineCurve, PlaneSurface, SurfaceGeometry};
use ogeom::math::{Circle, Direction, Frame, Plane, Point, Vector};
use ogeom::mesh::Deflection;
use ogeom::offset::{
    PipeCorners, PipeEvolution, PipeLaw, PipeLaws, SpineFunction, SpineLaw,
    make_pipe_shell_evolving, make_pipe_shell_laws,
};
use ogeom::topo::{Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();
const PI: f64 = core::f64::consts::PI;

/// The skin's tolerance for the skinned sweeps.
const SKIN: f64 = 1e-3;

fn line(model: &mut Model, a: Point, b: Point) -> Shape {
    let curve = LineCurve::segment(a, b, T).unwrap();
    let range = curve.domain();
    make_edge(model, Curve::Line(curve), range, T)
        .unwrap()
        .shape
}

/// A circle of `radius` about `centre` in the `xy` plane.
fn circle_wire(model: &mut Model, centre: Point, radius: f64) -> Shape {
    let frame = Frame::new(centre, Direction::Z, Direction::X, T).unwrap();
    let circle = Circle::new(frame, radius, T).unwrap();
    let edge = make_edge(model, CircleCurve::new(circle).into(), (0.0, 2.0 * PI), T)
        .unwrap()
        .shape;
    make_wire(model, &[edge], T).unwrap().shape
}

fn square_wire(model: &mut Model, half: f64) -> Shape {
    let pts = [(-half, -half), (half, -half), (half, half), (-half, half)]
        .map(|(x, y)| Point::new(x, y, 0.0));
    make_polygon(model, &pts, true, T).unwrap().shape
}

fn xy_face(model: &mut Model, wires: &[Shape]) -> Shape {
    let plane = PlaneSurface::new(Plane::new(Frame::WORLD));
    make_face(model, plane.into(), wires, T).unwrap().shape
}

/// Volume, and whether it was integrated exactly.
fn measure(model: &Model, shape: &Shape) -> (f64, bool) {
    let p = volume_properties(model, shape, Deflection::default(), T).unwrap();
    (p.mass, p.deflection == 0.0)
}

fn surfaces(model: &Model, shape: &Shape) -> Vec<SurfaceGeometry> {
    explore_unique(model, shape, ShapeType::Face)
        .unwrap()
        .iter()
        .map(|f| {
            let Some(NodeData::Face(d)) = model.node(f).map(|n| n.data()) else {
                panic!("a face holds no face data");
            };
            model.geometry().surface(d.surface).unwrap().clone()
        })
        .collect()
}

/// A section scaled linearly from one to `s` over a run `h` has area
/// `A (1 + (s - 1) λ)²`, so the solid holds `A h (1 + s + s²) / 3`
/// whatever the section turns through, wherever its centroid stays on a
/// spine it stands square to.
fn scaled_volume(area: f64, h: f64, s: f64) -> f64 {
    area * h * (1.0 + s + s * s) / 3.0
}

fn evolve(scale: f64, twist: f64) -> PipeEvolution {
    PipeEvolution { scale, twist }
}

#[test]
fn a_square_scaled_down_a_line_is_the_pyramid_s_frustum_exactly() {
    let (half, h, s) = (1.5, 6.0, 0.4);
    let mut model = Model::new();
    let wire = square_wire(&mut model, half);
    let profile = xy_face(&mut model, &[wire]);
    let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
    let solid = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        evolve(s, 0.0),
        SKIN,
        T,
    )
    .unwrap()
    .shape;
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let all = surfaces(&model, &solid);
    assert_eq!(all.len(), 6);
    assert!(all.iter().all(|s| matches!(s, SurfaceGeometry::Plane(_))));
    let (v, exact) = measure(&model, &solid);
    assert!(exact, "a solid of planes integrates exactly");
    let expected = scaled_volume(4.0 * half * half, h, s);
    assert!(
        (v - expected).abs() < expected * 1e-12,
        "{v} against {expected}"
    );
}

#[test]
fn a_disc_scaled_down_a_line_is_the_cone_s_frustum_and_to_naught_the_cone() {
    let (r, h) = (2.0, 5.0);
    for s in [2.5, 0.0] {
        let mut model = Model::new();
        let wire = circle_wire(&mut model, Point::ORIGIN, r);
        let profile = xy_face(&mut model, &[wire]);
        let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
        let solid = make_pipe_shell_evolving(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            evolve(s, 0.0),
            SKIN,
            T,
        )
        .unwrap()
        .shape;
        let found = check(&model, &solid, T).unwrap();
        assert!(found.is_valid(), "scale {s}: {found}");
        assert!(
            surfaces(&model, &solid)
                .iter()
                .any(|f| matches!(f, SurfaceGeometry::Cone(_))),
            "scale {s}: the wall is the cone"
        );
        let (v, _) = measure(&model, &solid);
        let expected = scaled_volume(PI * r * r, h, s);
        assert!(
            (v - expected).abs() < expected * 1e-9,
            "scale {s}: {v} against {expected}"
        );
    }
}

#[test]
fn a_square_scaled_to_naught_is_the_pyramid() {
    let (half, h) = (1.0, 3.0);
    let mut model = Model::new();
    let wire = square_wire(&mut model, half);
    let profile = xy_face(&mut model, &[wire]);
    let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
    let solid = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        evolve(0.0, 0.0),
        SKIN,
        T,
    )
    .unwrap()
    .shape;
    assert!(check(&model, &solid, T).unwrap().is_valid());
    assert_eq!(surfaces(&model, &solid).len(), 5, "four sides and a base");
    let (v, _) = measure(&model, &solid);
    let expected = 4.0 * half * half * h / 3.0;
    assert!(
        (v - expected).abs() < expected * 1e-12,
        "{v} against {expected}"
    );
}

/// The chord deflection a skinned sweep's volume is measured to.
const CHORD: f64 = 1e-4;

/// A skinned sweep's volume, integrated over its mesh at [`CHORD`].
fn measure_skinned(model: &Model, shape: &Shape) -> f64 {
    let fine = Deflection {
        chord: CHORD,
        ..Deflection::default()
    };
    volume_properties(model, shape, fine, T).unwrap().mass
}

/// The measured volume of a skinned sweep strays from the swept solid's by
/// at most the skin's tolerance and the mesh's chord over the walls' area:
/// their sum times a bound on the area, the perimeter at its largest times
/// the longest path a profile point runs.
fn skin_bound(perimeter: f64, longest_path: f64) -> f64 {
    (SKIN + CHORD) * perimeter * longest_path
}

#[test]
fn a_twisted_square_down_a_line_holds_the_prism_s_volume() {
    let (half, h, twist) = (1.0, 4.0, PI / 2.0);
    for s in [1.0, 0.5] {
        let mut model = Model::new();
        let wire = square_wire(&mut model, half);
        let profile = xy_face(&mut model, &[wire]);
        let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
        let solid = make_pipe_shell_evolving(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            evolve(s, twist),
            SKIN,
            T,
        )
        .unwrap()
        .shape;
        let found = check(&model, &solid, T).unwrap();
        assert!(found.is_usable(), "scale {s}: {found}");
        // The far end is the square turned a quarter turn and scaled: its
        // corners stand on the diagonals again.
        let corner = Point::new(half * s, half * s, h);
        let vertices = explore_unique(&model, &solid, ShapeType::Vertex).unwrap();
        let nearest = vertices
            .iter()
            .map(|v| {
                let Some(NodeData::Vertex(d)) = model.node(v).map(|n| n.data()) else {
                    panic!("a vertex holds no point");
                };
                d.point.distance(corner)
            })
            .fold(f64::INFINITY, f64::min);
        assert!(
            nearest < 10.0 * T.confusion(),
            "scale {s}: corner off by {nearest}"
        );
        let v = measure_skinned(&model, &solid);
        let expected = scaled_volume(4.0 * half * half, h, s);
        let radius = half * 2.0_f64.sqrt();
        let bound = skin_bound(8.0 * half, h.hypot(twist * radius));
        assert!(
            (v - expected).abs() < bound,
            "scale {s}: {v} against {expected}, bound {bound}"
        );
    }
}

/// A square centred on a quarter circle, scaled as it goes: its centroid
/// stays on the spine, so by Pappus the solid holds the sections' areas
/// summed along the spine's length.
#[test]
fn a_square_scaled_round_an_arc_holds_pappus_s_volume() {
    let (big, half, s) = (10.0, 1.0, 0.5);
    let mut model = Model::new();
    // The arc in the xz plane from the origin, leaving it along z.
    let frame = Frame::new(
        Point::new(big, 0.0, 0.0),
        Direction::new(Vector::new(0.0, 1.0, 0.0), T).unwrap(),
        Direction::new(Vector::new(-1.0, 0.0, 0.0), T).unwrap(),
        T,
    )
    .unwrap();
    let arc = Circle::new(frame, big, T).unwrap();
    let spine = make_edge(&mut model, CircleCurve::new(arc).into(), (0.0, PI / 2.0), T)
        .unwrap()
        .shape;
    let wire = square_wire(&mut model, half);
    let profile = xy_face(&mut model, &[wire]);
    let solid = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        evolve(s, 0.0),
        SKIN,
        T,
    )
    .unwrap()
    .shape;
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_usable(), "{found}");
    let length = big * PI / 2.0;
    let v = measure_skinned(&model, &solid);
    let expected = scaled_volume(4.0 * half * half, length, s);
    let reach = half * 2.0_f64.sqrt();
    let bound = skin_bound(8.0 * half, length * (1.0 + reach / big));
    assert!(
        (v - expected).abs() < bound,
        "{v} against {expected}, bound {bound}"
    );
}

#[test]
fn a_holed_square_scaled_down_a_line_keeps_its_hole() {
    let (half, r, h, s) = (2.0, 0.75, 5.0, 0.5);
    let mut model = Model::new();
    let outer = square_wire(&mut model, half);
    let hole = circle_wire(&mut model, Point::ORIGIN, r);
    let profile = xy_face(&mut model, &[outer, hole]);
    let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
    let solid = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        evolve(s, 0.0),
        SKIN,
        T,
    )
    .unwrap()
    .shape;
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_valid(), "{found}");
    let (v, _) = measure(&model, &solid);
    let expected = scaled_volume(4.0 * half * half - PI * r * r, h, s);
    assert!(
        (v - expected).abs() < expected * 1e-9,
        "{v} against {expected}"
    );
}

#[test]
fn an_evolving_pipe_refuses_what_it_does_not_build_by_name() {
    let mut model = Model::new();
    let wire = square_wire(&mut model, 1.0);
    let profile = xy_face(&mut model, &[wire]);
    let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, 3.0));
    let law = PipeLaw::RotationMinimizing;
    let err = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &law,
        evolve(-1.0, 0.0),
        SKIN,
        T,
    )
    .unwrap_err();
    assert!(err.to_string().contains("naught or more"), "{err}");
    let err = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &law,
        evolve(0.0, 1.0),
        SKIN,
        T,
    )
    .unwrap_err();
    assert!(err.to_string().contains("shrinks to a point"), "{err}");
    let err = make_pipe_shell_evolving(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::Fixed,
        evolve(1.0, 1.0),
        SKIN,
        T,
    )
    .unwrap_err();
    assert!(err.to_string().contains("fixed section"), "{err}");
}

/// A square twisted out to `0.3` rad half way down a line and back home by
/// its end, through the monotone cubic: each slice is the square turned,
/// so the solid holds the prism's `4 half² h`. The cubic keeps between its
/// points' values, so the square turns no further than `0.3`, where its
/// reach along `x` is `half (cos 0.3 + sin 0.3)`; a law overshooting its
/// points would reach further.
#[test]
fn a_square_twisted_out_and_back_by_a_law_holds_the_prism_s_volume() {
    let (half, h, most) = (1.0, 4.0, 0.3_f64);
    let mut model = Model::new();
    let wire = square_wire(&mut model, half);
    let profile = xy_face(&mut model, &[wire]);
    let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
    let laws = PipeLaws {
        scale: SpineLaw::Linear(1.0),
        twist: SpineLaw::Through(vec![(0.5, most), (1.0, 0.0)]),
        ..PipeLaws::default()
    };
    let solid = make_pipe_shell_laws(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        &laws,
        PipeCorners::Mitre,
        SKIN,
        T,
    )
    .unwrap()
    .shape;
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_usable(), "{found}");
    let v = measure_skinned(&model, &solid);
    let expected = 4.0 * half * half * h;
    let radius = half * 2.0_f64.sqrt();
    let bound = skin_bound(8.0 * half, h + 2.0 * most * radius);
    assert!(
        (v - expected).abs() < bound,
        "{v} against {expected}, bound {bound}"
    );
    let mesh = ogeom::mesh::triangulate(
        &model,
        &solid,
        Deflection {
            chord: CHORD,
            ..Deflection::default()
        },
        T,
    )
    .unwrap();
    let reach = mesh
        .positions
        .iter()
        .map(|p| p.x)
        .fold(f64::NEG_INFINITY, f64::max);
    let want = half * (most.cos() + most.sin());
    assert!(
        (reach - want).abs() < SKIN + CHORD,
        "reaches x {reach}, not {want}"
    );
}

/// A disc of radius `r` scaled through `2` half way down a line and back
/// to `1` at its end. The monotone cubic through `(0, 1)`, `(1/2, 2)`,
/// `(1, 1)` leaves each end along its secant (slope `2`, then `-2`) and
/// stands level at the middle, where the secants disagree in sign: on the
/// first half, `s = 1 + t + t² - t³` with `t = 2λ`, and the second half
/// its mirror. So `∫ s² dλ = ∫₀¹ (1 + t + t² - t³)² dt = 274/105`, and the
/// solid holds `π r² h · 274/105`.
#[test]
fn a_disc_scaled_out_and_back_by_a_law_holds_its_slices_volume() {
    let (r, h) = (1.0, 4.0);
    let mut model = Model::new();
    let wire = circle_wire(&mut model, Point::ORIGIN, r);
    let profile = xy_face(&mut model, &[wire]);
    let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
    let laws = PipeLaws {
        scale: SpineLaw::Through(vec![(0.5, 2.0), (1.0, 1.0)]),
        twist: SpineLaw::Linear(0.0),
        ..PipeLaws::default()
    };
    let solid = make_pipe_shell_laws(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        &laws,
        PipeCorners::Mitre,
        SKIN,
        T,
    )
    .unwrap()
    .shape;
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_usable(), "{found}");
    let v = measure_skinned(&model, &solid);
    let expected = PI * r * r * h * 274.0 / 105.0;
    let bound = skin_bound(2.0 * PI * 2.0 * r, h + 2.0 * r);
    assert!(
        (v - expected).abs() < bound,
        "{v} against {expected}, bound {bound}"
    );
}

#[test]
fn a_law_refuses_what_it_does_not_describe_by_name() {
    let mut model = Model::new();
    let wire = square_wire(&mut model, 1.0);
    let profile = xy_face(&mut model, &[wire]);
    let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, 3.0));
    let mut sweep = |scale: SpineLaw, twist: SpineLaw| {
        make_pipe_shell_laws(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            &PipeLaws {
                scale,
                twist,
                ..PipeLaws::default()
            },
            PipeCorners::Mitre,
            SKIN,
            T,
        )
        .unwrap_err()
        .to_string()
    };
    let err = sweep(
        SpineLaw::Through(vec![(0.5, 0.0), (1.0, 1.0)]),
        SpineLaw::Linear(0.0),
    );
    assert!(err.contains("stays above naught"), "{err}");
    let err = sweep(
        SpineLaw::Linear(1.0),
        SpineLaw::Through(vec![(0.5, 1.0), (0.5, 2.0), (1.0, 1.0)]),
    );
    assert!(err.contains("rise strictly"), "{err}");
    let err = sweep(SpineLaw::Linear(1.0), SpineLaw::Through(vec![(0.5, 1.0)]));
    assert!(err.contains("the last of them one"), "{err}");
    let err = sweep(
        SpineLaw::Linear(1.0),
        SpineLaw::Through(vec![(1.0, f64::NAN)]),
    );
    assert!(err.contains("finite values"), "{err}");
}

/// The spine up `z` from the origin to `(0, 0, l)`, then a sharp quarter
/// turn along `x` to `(l, 0, l)`: the inside of the corner is `+x`, its
/// axis `y`.
fn up_and_over(model: &mut Model, l: f64) -> Shape {
    make_polygon(
        model,
        &[
            Point::ORIGIN,
            Point::new(0.0, 0.0, l),
            Point::new(l, 0.0, l),
        ],
        false,
        T,
    )
    .unwrap()
    .shape
}

/// Simpson's rule over `[0, w]`: exact for a quadratic.
fn simpson(f: impl Fn(f64) -> f64, w: f64) -> f64 {
    w / 6.0 * (f(0.0) + 4.0 * f(w / 2.0) + f(w))
}

/// A square of half side `a` scaled linearly to `s` up and over the
/// corner, each leg `l` long: at height `z` up the first leg the scale is
/// `1 + k z`, a distance `x` along the second `1 + k (l + x)`, with
/// `k = (s - 1) / 2l`; at the corner it is `c = 1 + k l`. Each leg's solid
/// is the pyramid's frustum, `4 a² ∫ σ²`.
///
/// Mitred, each leg is run on by the prism of its corner section and cut
/// by the plane `x + z = l`. Either prism adds the wedge `a³ c³` outside
/// the corner; inside it, each leg loses what crosses the plane: at a
/// depth `w = l - z` below the corner the first leg's section reaches
/// `a σ(w)` past the spine, `σ(w) = c - k w`, so it loses
/// `∫ (a σ - w) 2 a σ dw` out to where `a σ = w`, and the second leg alike
/// with `σ(x) = c + k x`: quadratics, which Simpson's rule integrates
/// exactly.
///
/// Rounded, the outside half of the corner section, `2 a c` tall and
/// `a c` across, turns a quarter turn about the corner's axis, a quarter
/// drum `π a³ c³ / 2`, and the two legs overlap inside the corner where
/// `x ≤ a σ₁(w)` and `w ≤ a σ₂(x)`. For `s < 1` the second leg is the
/// thinner there, so the overlap is `∫ 2 a σ₂(x)` over the quadrilateral
/// those two lines bound: its area times the integrand at its centroid.
///
/// Both corners are planes and drums, the legs ruled exactly, so the
/// volumes land to the integrator's precision.
#[test]
fn a_square_scaled_round_a_sharp_corner_mitres_and_rounds_to_its_closed_form() {
    let (a, l, s) = (1.0_f64, 10.0_f64, 0.5_f64);
    let k = (s - 1.0) / (2.0 * l);
    let c = 1.0 + k * l;
    let frustum = |from: f64, to: f64| 4.0 * a * a * (to.powi(3) - from.powi(3)) / (3.0 * k);
    let legs = frustum(1.0, c) + frustum(c, s);
    let lost = |sigma: &dyn Fn(f64) -> f64, reach: f64| {
        simpson(|w| (a * sigma(w) - w) * 2.0 * a * sigma(w), reach)
    };
    let first = lost(&|w| c - k * w, a * c / (1.0 + a * k));
    let second = lost(&|x| c + k * x, a * c / (1.0 - a * k));
    let mitred = legs + 2.0 * (a * c).powi(3) - first - second;
    let (cc, m) = (a * c, -a * k);
    let u = cc * (1.0 + m) / (1.0 + m * m);
    let quad = [(0.0, 0.0), (cc, 0.0), (u, cc - m * u), (0.0, cc)];
    let (mut area, mut moment) = (0.0, 0.0);
    for i in 0..4 {
        let (p, q) = (quad[i], quad[(i + 1) % 4]);
        let cross = p.0 * q.1 - q.0 * p.1;
        area += cross / 2.0;
        moment += (p.0 + q.0) * cross / 6.0;
    }
    let overlap = area * 2.0 * a * (c + k * moment / area);
    let rounded = legs + PI * (a * c).powi(3) / 2.0 - overlap;

    for (corners, want) in [(PipeCorners::Mitre, mitred), (PipeCorners::Round, rounded)] {
        let mut model = Model::new();
        let wire = square_wire(&mut model, a);
        let profile = xy_face(&mut model, &[wire]);
        let spine = up_and_over(&mut model, l);
        let laws = PipeLaws {
            scale: SpineLaw::Linear(s),
            twist: SpineLaw::Linear(0.0),
            ..PipeLaws::default()
        };
        let solid = make_pipe_shell_laws(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            &laws,
            corners,
            SKIN,
            T,
        )
        .unwrap_or_else(|e| panic!("{corners:?}: {e}"))
        .shape;
        let found = check(&model, &solid, T).unwrap();
        assert!(found.is_valid(), "{corners:?}: {found}");
        let (v, exact) = measure(&model, &solid);
        assert!(exact, "{corners:?}: not integrated exactly");
        assert!(
            (v - want).abs() < want * 1e-9,
            "{corners:?}: {v} against {want}"
        );
    }
}

/// A square of half side `a` turned through `θ` about its centre: the
/// stretch of `x` it covers at each height `y`, as `(low, high)`, or
/// `None` past its top and bottom.
fn turned_square_at(a: f64, theta: f64, y: f64) -> Option<(f64, f64)> {
    let (sin, cos) = theta.sin_cos();
    let corners =
        [(a, a), (-a, a), (-a, -a), (a, -a)].map(|(x, z)| (x * cos - z * sin, x * sin + z * cos));
    let mut xs = Vec::new();
    for i in 0..4 {
        let (p, q) = (corners[i], corners[(i + 1) % 4]);
        if (p.1 - y) * (q.1 - y) <= 0.0 && p.1 != q.1 {
            xs.push(p.0 + (q.0 - p.0) * (y - p.1) / (q.1 - p.1));
        }
    }
    let low = xs.iter().copied().fold(f64::INFINITY, f64::min);
    let high = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    (low <= high).then_some((low, high))
}

/// `∫ f(y) dy` over the turned square's height, Simpson's rule between
/// every height where a corner stands or an edge crosses `x = 0`: `f` is
/// a quadratic of `y` between those, so the rule is exact.
fn over_turned_square(a: f64, theta: f64, f: impl Fn(f64, f64) -> f64) -> f64 {
    let (sin, cos) = theta.sin_cos();
    let mut breaks: Vec<f64> = [(a, a), (-a, a), (-a, -a), (a, -a)]
        .iter()
        .map(|(x, z)| x * sin + z * cos)
        .collect();
    for i in 0..4 {
        let c = |k: usize| {
            let (x, z) = [(a, a), (-a, a), (-a, -a), (a, -a)][k % 4];
            (x * cos - z * sin, x * sin + z * cos)
        };
        let (p, q) = (c(i), c(i + 1));
        if p.0 * q.0 < 0.0 {
            breaks.push(p.1 + (q.1 - p.1) * (-p.0) / (q.0 - p.0));
        }
    }
    breaks.sort_by(f64::total_cmp);
    let at = |y: f64| turned_square_at(a, theta, y).map_or(0.0, |(lo, hi)| f(lo, hi));
    breaks
        .windows(2)
        .map(|w| (w[1] - w[0]) / 6.0 * (at(w[0]) + 4.0 * at(f64::midpoint(w[0], w[1])) + at(w[1])))
        .sum()
}

/// A square twisted up the first leg, between a tenth and two fifths of
/// the way along the spine, by the monotone cubic through `(0.1, 0)`,
/// `(0.4, θ)` and `(1, θ)`: level spans hold their value exactly, so the
/// section is the square before and the square turned through `θ` after,
/// through the corner and down the second leg. Every slice is the square
/// turned, so each leg holds `4 a² l`.
///
/// Mitred, both legs are cut on the plane halving the corner through the
/// spine, which takes off inside what the prisms add outside: `8 a² l`.
/// Rounded, the half of the corner section outside the corner's axis
/// turns a quarter turn about it, adding `π/2 ∫∫ -x` over that half by
/// Pappus, and the legs overlap inside the corner where both sections
/// stand: carried across the corner, the second leg's section is the
/// first's with its depth `w` below the corner for `x`, so at each height
/// `y` the overlap is the square of the stretch `0 ≤ x ≤ high(y)` the
/// section covers. Both are integrated exactly over the section's height.
/// The far end is the square turned through `θ`.
#[test]
fn a_square_twisted_before_a_sharp_corner_mitres_and_rounds_to_its_closed_form() {
    twisted_before_a_sharp_corner(0.3);
}

/// As above, turned exactly an eighth of a turn: two corners of the
/// section then stand on the corner's axis, and the plane cutting the
/// runs there holds the skinned run's side edges through them.
#[test]
fn a_square_twisted_an_eighth_turn_before_a_sharp_corner_mitres_and_rounds() {
    twisted_before_a_sharp_corner(PI / 4.0);
}

fn twisted_before_a_sharp_corner(theta: f64) {
    let (a, l) = (1.0_f64, 10.0_f64);
    let legs = 8.0 * a * a * l;
    let turned = over_turned_square(a, theta, |lo, hi| {
        if lo < 0.0 {
            (lo * lo - hi.min(0.0).powi(2)) / 2.0
        } else {
            0.0
        }
    });
    let overlap = over_turned_square(a, theta, |lo, hi| (hi - lo.max(0.0)).max(0.0).powi(2));
    let laws = PipeLaws {
        scale: SpineLaw::Linear(1.0),
        twist: SpineLaw::Through(vec![(0.1, 0.0), (0.4, theta), (1.0, theta)]),
        ..PipeLaws::default()
    };
    for (corners, want) in [
        (PipeCorners::Mitre, legs),
        (PipeCorners::Round, legs + PI / 2.0 * turned - overlap),
    ] {
        let mut model = Model::new();
        let wire = square_wire(&mut model, a);
        let profile = xy_face(&mut model, &[wire]);
        let spine = up_and_over(&mut model, l);
        let solid = make_pipe_shell_laws(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            &laws,
            corners,
            SKIN,
            T,
        )
        .unwrap_or_else(|e| panic!("{corners:?}: {e}"))
        .shape;
        let found = check(&model, &solid, T).unwrap();
        assert!(found.is_usable(), "{corners:?}: {found}");
        let v = measure_skinned(&model, &solid);
        let bound = skin_bound(8.0 * a, 2.0 * l + theta * a * 2.0_f64.sqrt());
        assert!(
            (v - want).abs() < bound,
            "{corners:?}: {v} against {want}, bound {bound}"
        );
        // The far end's corners stand on the circle of the half diagonal,
        // an eighth of a turn and `θ` (either way) off its axes.
        let far = Point::new(l, 0.0, l);
        let mut ends = 0;
        for v in explore_unique(&model, &solid, ShapeType::Vertex).unwrap() {
            let Some(NodeData::Vertex(data)) = model.node(&v).map(|n| n.data()) else {
                panic!("a vertex holds no point");
            };
            let p = data.point;
            if (p.x - l).abs() > 10.0 * T.confusion() {
                continue;
            }
            ends += 1;
            let off = p - far;
            assert!(
                (off.magnitude() - a * 2.0_f64.sqrt()).abs() < 10.0 * T.confusion(),
                "{corners:?}: a far corner at {p:?}"
            );
            let turn = off.y.atan2(-off.z).rem_euclid(PI / 2.0) - PI / 4.0;
            assert!(
                (turn.abs() - theta).abs() < 1e-9,
                "{corners:?}: the far end turned {turn}"
            );
        }
        assert_eq!(ends, 4, "{corners:?}");
    }
}

/// The area of the square of half side `a` turned through `theta` about
/// its centre that lies past `x = w`: the turned square clipped by that
/// line, by the shoelace formula.
fn turned_square_past(a: f64, theta: f64, w: f64) -> f64 {
    let (sin, cos) = theta.sin_cos();
    let corners =
        [(a, a), (-a, a), (-a, -a), (a, -a)].map(|(x, y)| (x * cos - y * sin, x * sin + y * cos));
    let mut kept: Vec<(f64, f64)> = Vec::new();
    for i in 0..4 {
        let (p, q) = (corners[i], corners[(i + 1) % 4]);
        if p.0 > w {
            kept.push(p);
        }
        if (p.0 > w) != (q.0 > w) {
            kept.push((w, p.1 + (q.1 - p.1) * (w - p.0) / (q.0 - p.0)));
        }
    }
    (0..kept.len())
        .map(|i| {
            let (p, q) = (kept[i], kept[(i + 1) % kept.len()]);
            p.0 * q.1 - q.0 * p.1
        })
        .sum::<f64>()
        / 2.0
}

/// `∫ f` over `[0, w]` by the four-point Gauss–Legendre rule on `panels`
/// even panels.
fn gauss(f: impl Fn(f64) -> f64, w: f64, panels: u32) -> f64 {
    const NODES: [(f64, f64); 4] = [
        (-0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
        (-0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.339_981_043_584_856_3, 0.652_145_154_862_546_1),
        (0.861_136_311_594_052_6, 0.347_854_845_137_453_8),
    ];
    let h = w / f64::from(panels);
    (0..panels)
        .map(|k| {
            let mid = (f64::from(k) + 0.5) * h;
            NODES
                .iter()
                .map(|(x, wt)| wt * f(mid + x * h / 2.0))
                .sum::<f64>()
                * h
                / 2.0
        })
        .sum()
}

/// A square of half side `a` twisted linearly by `total` up and over the
/// sharp L, still turning through the mitred corner, under the
/// rotation-minimizing law: the section at height `z` up the first leg is
/// the square turned `θ' z`, and a distance `x` along the second the
/// square turned `θ' (l + x)`, `θ' = total / 2l`, `θc = θ' l` at the
/// corner.
///
/// On the mitre plane `x + z = l`, a point `w` inside the corner (`x = w`)
/// at `y` lies in the first leg's section `w` below the corner, turned
/// `θc - θ' w`, at `(w, y)` of its own plane, and in the second's `w` past
/// it, turned `θc + θ' w`, at the same `(w, y)`; outside the corner both
/// legs' prisms carry the corner section. Each leg holds `4 a² l`, a
/// twist keeping every slice the square; each prism adds what lies past
/// `x = -w` of the corner section, which is what lies past `w` of it, and
/// each leg loses what lies past `w` of its own section `w` from the
/// corner:
/// `V = 8 a² l + Σ± ∫ [past(θc, w) - past(θc ∓ θ' w, w)] dw`.
///
/// The two legs' cuts of the plane differ inside the corner, so the solid
/// has faces there: the first leg's region less the second's, looking up
/// the plane's normal `(1, 0, 1)`, and the second's less the first's
/// looking down it, each `√2 ∫ |slice₁(w) \ slice₂(w)| dw` across the
/// stretches of `y` the turned squares cover at `x = w`.
///
/// Every integrand is smooth between the heights where a corner of a
/// section crosses `x = w`, so the Gauss–Legendre sums are exact to far
/// below the skin. The volume strays from the closed form by at most the
/// skin bound; each mitre face's area by the skin's tolerance and the
/// mesh's chord, within `√2` of the plane each, over both cuts' lengths,
/// at most the section's perimeter `√2` stretched across the plane.
#[test]
fn a_square_twisted_through_a_mitred_corner_holds_its_closed_form_and_mitre_faces() {
    let (a, l) = (1.0_f64, 10.0_f64);
    let reach = a * 2.0_f64.sqrt();
    for total in [0.6, 2.0] {
        let rate = total / (2.0 * l);
        let corner = rate * l;
        let lost = |sign: f64| {
            gauss(
                |w| {
                    turned_square_past(a, corner, w)
                        - turned_square_past(a, corner + sign * rate * w, w)
                },
                reach * 1.2,
                20_000,
            )
        };
        let want = 8.0 * a * a * l + lost(-1.0) + lost(1.0);
        // The stretch of `y` the square turned `θ` covers at `x = w`, the
        // turned square's `x` stretch at `y = w` mirrored across `x = y`.
        let slice = |theta: f64, w: f64| turned_square_at(a, -theta, w);
        let only = |own: f64, other: f64| {
            2.0_f64.sqrt()
                * gauss(
                    |w| {
                        let mine = slice(corner + own * rate * w, w);
                        let theirs = slice(corner + other * rate * w, w);
                        match (mine, theirs) {
                            (None, _) => 0.0,
                            (Some((lo, hi)), None) => hi - lo,
                            (Some((lo, hi)), Some((low, high))) => {
                                (hi - lo) - (hi.min(high) - lo.max(low)).max(0.0)
                            }
                        }
                    },
                    reach * 1.2,
                    20_000,
                )
        };
        let (up, down) = (only(-1.0, 1.0), only(1.0, -1.0));

        let mut model = Model::new();
        let wire = square_wire(&mut model, a);
        let profile = xy_face(&mut model, &[wire]);
        let spine = up_and_over(&mut model, l);
        let laws = PipeLaws {
            scale: SpineLaw::Linear(1.0),
            twist: SpineLaw::Linear(total),
            ..PipeLaws::default()
        };
        let solid = make_pipe_shell_laws(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            &laws,
            PipeCorners::Mitre,
            SKIN,
            T,
        )
        .unwrap_or_else(|e| panic!("{total}: {e}"))
        .shape;
        let found = check(&model, &solid, T).unwrap();
        assert!(found.is_valid(), "{total}: {found}");
        let v = measure_skinned(&model, &solid);
        let bound = skin_bound(8.0 * a, 2.0 * l + total * reach);
        assert!(
            (v - want).abs() < bound,
            "{total}: {v} against {want}, bound {bound}"
        );

        let normal = Vector::new(1.0, 0.0, 1.0) / 2.0_f64.sqrt();
        let (mut looking_up, mut looking_down) = (0.0, 0.0);
        for face in explore_unique(&model, &solid, ShapeType::Face).unwrap() {
            let (at, n) = ogeom::algo::face_normal(&model, &face, T).unwrap();
            if n.cross(normal).magnitude() > 1e-9 || (at.x + at.z - l).abs() > 1e-6 {
                continue;
            }
            let fine = Deflection {
                chord: CHORD,
                ..Deflection::default()
            };
            let area = ogeom::algo::surface_properties(&model, &face, fine, T)
                .unwrap()
                .mass;
            if n.dot(normal) > 0.0 {
                looking_up += area;
            } else {
                looking_down += area;
            }
        }
        let bound = (SKIN + CHORD) * 2.0_f64.sqrt() * 2.0 * 8.0 * a * 2.0_f64.sqrt();
        for (case, measured, wanted) in [("up", looking_up, up), ("down", looking_down, down)] {
            assert!(
                (measured - wanted).abs() < bound,
                "{total}: the face looking {case} holds {measured} against {wanted}, bound {bound}"
            );
        }
    }
}

/// The spine up `z` from the origin to `(0, 0, l)`, then a sharp quarter
/// turn onto the quarter arc of radius `r` about `(0, 0, l + r)`, leaving
/// along `-x` and ending at `(-r, 0, l + r)`: the corner's axis is `y`.
fn up_and_round(model: &mut Model, l: f64, r: f64) -> Shape {
    let frame = Frame::new(Point::new(0.0, 0.0, l + r), Direction::Y, Direction::X, T).unwrap();
    let circle = Circle::new(frame, r, T).unwrap();
    let arc = make_edge(model, CircleCurve::new(circle).into(), (PI / 2.0, PI), T)
        .unwrap()
        .shape;
    let (joint, _) = ogeom::algo::edge_vertices(model, &arc).unwrap().unwrap();
    let start = ogeom::algo::make_vertex(model, Point::ORIGIN).shape;
    let segment = LineCurve::segment(Point::ORIGIN, Point::new(0.0, 0.0, l), T).unwrap();
    let range = segment.domain();
    let rise = ogeom::algo::make_edge_between(model, segment.into(), range, &start, &joint, T)
        .unwrap()
        .shape;
    make_wire(model, &[rise, arc], T).unwrap().shape
}

/// A square of half side `a` scaled linearly from one to `s` along
/// [`up_and_round`], extended at the corner. The scale falls at
/// `m = (1 - s) / L` per unit of the spine's length `L = l + π r / 2`,
/// to `c = 1 - m l` at the corner, so the two runs hold
/// `4 a² ∫ σ² = 4 a² (1 - s³) / 3m` (Pappus round the arc, the section
/// centred on the spine).
///
/// The first run goes on past the corner by its section's prism, `2 a c`
/// square and `a c` up, `4 a³ c³`; the arc's run reaches back along `+x`
/// by its start section's prism, which lies within the first run and that
/// prism. The arc's section at `φ` round it, `2 a σ(φ)` across `y` with
/// `σ(φ) = c - m r φ`, is the thinnest of the three wherever they overlap,
/// so the overlap is `2 a σ(φ)` over the part of its column, at radius
/// `ρ` from the arc's centre, that stands within the prism (up to `a c`
/// above the corner and within `a c` of the first run's axis) or within
/// the first run below the corner (within `a (c + m w)` of its axis at a
/// depth `w`), integrated by the midpoint rule over a hundred thousand
/// columns.
fn scaled_into_an_arc(a: f64, s: f64, l: f64, r: f64) -> f64 {
    const COLUMNS: u32 = 100_000;
    let m = (1.0 - s) / (l + r * PI / 2.0);
    let c = 1.0 - m * l;
    let runs = 4.0 * a * a * (1.0 - s.powi(3)) / (3.0 * m);
    let step = PI / 2.0 / f64::from(COLUMNS);
    let mut overlap = 0.0;
    for k in 0..COLUMNS {
        let phi = (f64::from(k) + 0.5) * step;
        let sigma = c - m * r * phi;
        let (sin, cos) = phi.sin_cos();
        // The column at `ρ` stands `r - ρ cos φ` above the corner and
        // `ρ sin φ` off the first run's axis.
        let low = (r - a * sigma).max((r - a * c) / cos);
        let mut high = a * c / sin;
        if high * cos > r {
            let leaning = sin - a * m * cos;
            high = if leaning > 0.0 {
                (a * c - a * m * r) / leaning
            } else {
                f64::INFINITY
            };
        }
        let high = high.min(r + a * sigma);
        if high > low {
            overlap += 2.0 * a * sigma * (high * high - low * low) / 2.0 * step;
        }
    }
    runs + 4.0 * (a * c).powi(3) - overlap
}

/// The arc's run of [`scaled_into_an_arc`] starts on the first run's end
/// section turned about the corner, so its side walls start on edges lying
/// in the side planes of the first run's prism, and leave those planes at
/// the sine `a m` at which they narrow, here 0.012 and 0.024.
#[test]
fn a_square_scaled_into_an_arc_at_an_extended_corner_holds_its_volume() {
    let (l, r, s) = (10.0_f64, 20.0_f64, 0.5_f64);
    for a in [1.0, 2.0] {
        let mut model = Model::new();
        let wire = square_wire(&mut model, a);
        let profile = xy_face(&mut model, &[wire]);
        let spine = up_and_round(&mut model, l, r);
        let laws = PipeLaws {
            scale: SpineLaw::Linear(s),
            twist: SpineLaw::Linear(0.0),
            ..PipeLaws::default()
        };
        let solid = make_pipe_shell_laws(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            &laws,
            PipeCorners::Extended,
            SKIN,
            T,
        )
        .unwrap_or_else(|e| panic!("half side {a}: {e}"))
        .shape;
        let found = check(&model, &solid, T).unwrap();
        assert!(found.is_valid(), "half side {a}: {found}");
        let v = measure_skinned(&model, &solid);
        let want = scaled_into_an_arc(a, s, l, r);
        // The furthest a profile point runs: up the first run and its
        // prism, then round the arc at the section's widest reach.
        let bound = skin_bound(8.0 * a, l + a + (r + a * 2.0_f64.sqrt()) * PI / 2.0);
        assert!(
            (v - want).abs() < bound,
            "half side {a}: {v} against {want}, bound {bound}"
        );
    }
}

/// A law given as a function of the share of the spine run.
fn law_of(f: impl Fn(f64) -> f64 + Send + Sync + 'static) -> SpineLaw {
    SpineLaw::Function(SpineFunction::new(f))
}

/// The solid's section by the plane `z = z0`, as the loops it cuts there in
/// that plane's `(x, y)`, each a polyline whose points lie on the section.
fn cut_at(model: &mut Model, solid: &Shape, z0: f64) -> Vec<Vec<(f64, f64)>> {
    let plane =
        Plane::new(Frame::new(Point::new(0.0, 0.0, z0), Direction::Z, Direction::X, T).unwrap());
    let fine = Deflection {
        chord: CHORD,
        ..Deflection::default()
    };
    let view = ogeom::hlr::section(model, solid, &plane, fine, T).unwrap();
    view.outline
        .iter()
        .map(|ring| ring.iter().map(|p| (p.x, p.y)).collect())
        .collect()
}

/// How far `(x, y)` stands from the boundary of the square of half side
/// `half` turned through `theta` about the origin.
fn off_turned_square(half: f64, theta: f64, (x, y): (f64, f64)) -> f64 {
    let (sin, cos) = theta.sin_cos();
    let (u, v) = (x * cos + y * sin, -x * sin + y * cos);
    let (du, dv) = (u.abs() - half, v.abs() - half);
    if du <= 0.0 && dv <= 0.0 {
        -du.max(dv)
    } else {
        du.max(0.0).hypot(dv.max(0.0))
    }
}

/// A disc of radius `r` scaled by `s(λ) = 1 + sin(πλ) / 2` down a line of
/// length `h`, the law given as a function: each slice is the disc scaled
/// by `s`, so the solid holds `π r² h ∫ s² dλ = π r² h (9/8 + 2/π)`, and
/// its section at a share `λ` is the circle of radius `r s(λ)`.
#[test]
fn a_disc_scaled_by_a_function_holds_its_slices_volume_and_sections() {
    let (r, h) = (1.0, 4.0);
    let s = |lambda: f64| 1.0 + (PI * lambda).sin() / 2.0;
    let mut model = Model::new();
    let wire = circle_wire(&mut model, Point::ORIGIN, r);
    let profile = xy_face(&mut model, &[wire]);
    let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
    let laws = PipeLaws {
        scale: law_of(s),
        twist: SpineLaw::Linear(0.0),
        ..PipeLaws::default()
    };
    let solid = make_pipe_shell_laws(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        &laws,
        PipeCorners::Mitre,
        SKIN,
        T,
    )
    .unwrap()
    .shape;
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_usable(), "{found}");
    let v = measure_skinned(&model, &solid);
    let expected = PI * r * r * h * (9.0 / 8.0 + 2.0 / PI);
    // A rim point runs `h` along and `r` out and back.
    let bound = skin_bound(2.0 * PI * 1.5 * r, h + r);
    assert!(
        (v - expected).abs() < bound,
        "{v} against {expected}, bound {bound}"
    );
    for lambda in [0.2, 0.5, 0.85] {
        let rings = cut_at(&mut model, &solid, lambda * h);
        assert_eq!(rings.len(), 1, "at {lambda}");
        let want = r * s(lambda);
        for &(x, y) in &rings[0] {
            let off = (x.hypot(y) - want).abs();
            assert!(
                off < SKIN,
                "at {lambda}: ({x}, {y}) stands {off} off radius {want}"
            );
        }
    }
}

/// A square turned by `θ(λ) = 1.2 sin²(πλ/2)` down a line, the law given
/// as a function: every slice is the square turned, so the solid holds the
/// prism's `4 a² h` for any twist law, and its section at a share `λ` is
/// the square turned through `θ(λ)`.
#[test]
fn a_square_twisted_by_a_function_holds_the_prism_s_volume_and_turns_its_sections() {
    let (a, h) = (1.0, 4.0);
    let theta = |lambda: f64| 1.2 * (PI * lambda / 2.0).sin().powi(2);
    let mut model = Model::new();
    let wire = square_wire(&mut model, a);
    let profile = xy_face(&mut model, &[wire]);
    let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
    let laws = PipeLaws {
        scale: SpineLaw::Linear(1.0),
        twist: law_of(theta),
        ..PipeLaws::default()
    };
    let solid = make_pipe_shell_laws(
        &mut model,
        &profile,
        &spine,
        &PipeLaw::RotationMinimizing,
        &laws,
        PipeCorners::Mitre,
        SKIN,
        T,
    )
    .unwrap()
    .shape;
    let found = check(&model, &solid, T).unwrap();
    assert!(found.is_usable(), "{found}");
    let v = measure_skinned(&model, &solid);
    let expected = 4.0 * a * a * h;
    let bound = skin_bound(8.0 * a, h + 1.2 * a * 2.0_f64.sqrt());
    assert!(
        (v - expected).abs() < bound,
        "{v} against {expected}, bound {bound}"
    );
    for lambda in [0.25, 0.5, 0.9] {
        let rings = cut_at(&mut model, &solid, lambda * h);
        assert_eq!(rings.len(), 1, "at {lambda}");
        for &p in &rings[0] {
            let off = off_turned_square(a, theta(lambda), p);
            assert!(
                off < SKIN,
                "at {lambda}: {p:?} stands {off} off the turned square"
            );
        }
    }
}

/// A law given as a function that is the linear one builds what the linear
/// law does within the skin's tolerance: here the far end, a square scaled
/// by a half and turned a quarter turn, lands on the same corners.
#[test]
fn a_function_law_reads_as_the_law_it_gives() {
    let (a, h) = (1.0, 3.0);
    let sweep = |laws: PipeLaws| {
        let mut model = Model::new();
        let wire = square_wire(&mut model, a);
        let profile = xy_face(&mut model, &[wire]);
        let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, h));
        let solid = make_pipe_shell_laws(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            &laws,
            PipeCorners::Mitre,
            SKIN,
            T,
        )
        .unwrap()
        .shape;
        assert!(check(&model, &solid, T).unwrap().is_usable());
        let mut far: Vec<(f64, f64)> = explore_unique(&model, &solid, ShapeType::Vertex)
            .unwrap()
            .iter()
            .filter_map(|v| match model.node(v).map(|n| n.data()) {
                Some(NodeData::Vertex(d)) if (d.point.z - h).abs() < T.confusion() => {
                    Some((d.point.x, d.point.y))
                }
                _ => None,
            })
            .collect();
        far.sort_by(|p, q| p.0.total_cmp(&q.0));
        (measure_skinned(&model, &solid), far)
    };
    let (v_linear, far_linear) = sweep(PipeLaws {
        scale: SpineLaw::Linear(0.5),
        twist: SpineLaw::Linear(PI / 2.0),
        ..PipeLaws::default()
    });
    let (v_function, far_function) = sweep(PipeLaws {
        scale: law_of(|l| 1.0 - l / 2.0),
        twist: law_of(|l| PI / 2.0 * l),
        ..PipeLaws::default()
    });
    let expected = scaled_volume(4.0 * a * a, h, 0.5);
    let bound = skin_bound(8.0 * a, h + PI / 2.0 * a * 2.0_f64.sqrt());
    for v in [v_linear, v_function] {
        assert!((v - expected).abs() < bound, "{v} against {expected}");
    }
    assert_eq!(far_linear.len(), 4);
    assert_eq!(far_function.len(), 4);
    for (p, q) in far_linear.iter().zip(&far_function) {
        assert!(
            (p.0 - q.0).hypot(p.1 - q.1) < T.confusion() * 10.0,
            "{p:?} against {q:?}"
        );
    }
}

#[test]
fn a_function_law_refuses_what_it_does_not_describe_by_name() {
    let mut model = Model::new();
    let wire = square_wire(&mut model, 1.0);
    let profile = xy_face(&mut model, &[wire]);
    let spine = line(&mut model, Point::ORIGIN, Point::new(0.0, 0.0, 3.0));
    let mut sweep = |scale: SpineLaw, twist: SpineLaw| {
        make_pipe_shell_laws(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            &PipeLaws {
                scale,
                twist,
                ..PipeLaws::default()
            },
            PipeCorners::Mitre,
            SKIN,
            T,
        )
        .unwrap_err()
        .to_string()
    };
    let err = sweep(law_of(|l| 2.0 - l), SpineLaw::Linear(0.0));
    assert!(err.contains("gives 1 at the spine's start"), "{err}");
    let err = sweep(SpineLaw::Linear(1.0), law_of(|l| 0.1 + l));
    assert!(err.contains("gives 0 at the spine's start"), "{err}");
    // Through naught half way, back up to one: positive at both ends.
    let err = sweep(law_of(|l| (1.0 - 2.0 * l).powi(2)), SpineLaw::Linear(0.0));
    assert!(err.contains("stays above naught"), "{err}");
    let err = sweep(
        SpineLaw::Linear(1.0),
        law_of(|l| if l > 0.5 { f64::NAN } else { l }),
    );
    assert!(err.contains("finite values"), "{err}");
}

/// A square of half side `a` scaled by `σ(λ) = 1 - λ² / 2` up and over the
/// sharp corner of [`up_and_over`], the law given as a function and still
/// moving through the corner, where it stands at `c = σ(1/2)`. As for the
/// linear law, each leg holds `4 a² ∫ σ²` along it, so the two
/// `8 a² l (1 - 1/3 + 1/20)`.
///
/// Mitred, each prism adds the wedge `a³ c³` outside the corner and each
/// leg loses inside what crosses the plane `x + z = l`: at a depth `w`
/// below the corner the first leg reaches `a σ₁(w)` past the spine and
/// loses `∫ (a σ₁ - w) 2 a σ₁ dw` out to where `a σ₁ = w`, the second
/// alike a distance `x` along it.
///
/// Rounded, the quarter drum `π a³ c³ / 2` is added and the two legs'
/// overlap inside the corner taken off. The law falls along the spine, so
/// the first leg stands wider there than the second and the overlap is
/// `∫ 2 a σ₂(x) · |{w ≤ a σ₂(x) : x ≤ a σ₁(w)}| dx`: in the column at `x`,
/// `x ≤ a σ₁(w)` from the depth `w₀(x)` where they are equal (naught while
/// `x ≤ a c`), as `σ₁` rises with `w`.
///
/// Each integral is taken numerically over its smooth stretch.
#[test]
fn a_square_scaled_by_a_function_round_a_sharp_corner_mitres_and_rounds_to_its_integral() {
    let (a, l) = (1.0_f64, 10.0_f64);
    let sigma = |lambda: f64| 1.0 - lambda * lambda / 2.0;
    let first = move |w: f64| sigma((l - w) / (2.0 * l));
    let second = move |x: f64| sigma((l + x) / (2.0 * l));
    let c = sigma(0.5);
    let legs = 8.0 * a * a * l * (1.0 - 1.0 / 3.0 + 1.0 / 20.0);
    // The root of a falling `g` in `[0, hi]`, by bisection.
    let root = |g: &dyn Fn(f64) -> f64, hi: f64| {
        let (mut lo, mut hi) = (0.0_f64, hi);
        for _ in 0..200 {
            let mid = f64::midpoint(lo, hi);
            if g(mid) > 0.0 { lo = mid } else { hi = mid }
        }
        f64::midpoint(lo, hi)
    };
    // Composite Simpson's rule over `[0, w]`.
    let simpson_n = |f: &dyn Fn(f64) -> f64, w: f64| {
        const N: u32 = 2000;
        let step = w / f64::from(N);
        (0..N)
            .map(|k| {
                let x = f64::from(k) * step;
                step / 6.0 * (f(x) + 4.0 * f(x + step / 2.0) + f(x + step))
            })
            .sum::<f64>()
    };
    let lost = |s: &dyn Fn(f64) -> f64| {
        let reach = root(&|w| a * s(w) - w, 4.0 * a);
        simpson_n(&|w| (a * s(w) - w) * 2.0 * a * s(w), reach)
    };
    let mitred = legs + 2.0 * (a * c).powi(3) - lost(&first) - lost(&second);
    // The column at `x` overlaps from `w₀(x)` to `a σ₂(x)`.
    let depth = |x: f64| {
        if x <= a * c {
            0.0
        } else {
            root(&|w| x - a * first(w), 4.0 * a)
        }
    };
    let across = |x: f64| (a * second(x) - depth(x)).max(0.0);
    let far = root(&|x| a * second(x) - depth(x), 4.0 * a);
    let overlap = simpson_n(&|x| 2.0 * a * second(x) * across(x), a * c)
        + simpson_n(
            &|x| 2.0 * a * second(a * c + x) * across(a * c + x),
            far - a * c,
        );
    let rounded = legs + PI * (a * c).powi(3) / 2.0 - overlap;

    for (corners, want) in [(PipeCorners::Mitre, mitred), (PipeCorners::Round, rounded)] {
        let mut model = Model::new();
        let wire = square_wire(&mut model, a);
        let profile = xy_face(&mut model, &[wire]);
        let spine = up_and_over(&mut model, l);
        let laws = PipeLaws {
            scale: law_of(sigma),
            twist: SpineLaw::Linear(0.0),
            ..PipeLaws::default()
        };
        let solid = make_pipe_shell_laws(
            &mut model,
            &profile,
            &spine,
            &PipeLaw::RotationMinimizing,
            &laws,
            corners,
            SKIN,
            T,
        )
        .unwrap_or_else(|e| panic!("{corners:?}: {e}"))
        .shape;
        let found = check(&model, &solid, T).unwrap();
        assert!(found.is_valid(), "{corners:?}: {found}");
        let v = measure_skinned(&model, &solid);
        let bound = skin_bound(8.0 * a, 2.0 * l + 2.0 * a);
        assert!(
            (v - want).abs() < bound,
            "{corners:?}: {v} against {want}, bound {bound}"
        );
    }
}
