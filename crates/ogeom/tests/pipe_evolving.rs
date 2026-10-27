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
    PipeCorners, PipeEvolution, PipeLaw, PipeLaws, SpineLaw, make_pipe_shell_evolving,
    make_pipe_shell_laws,
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
            &PipeLaws { scale, twist },
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
