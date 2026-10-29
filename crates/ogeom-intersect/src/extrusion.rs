//! A plane across a linear extrusion, in closed form.
//!
//! An extrusion sweeps a curve `C` along a unit direction `d`: its point at
//! `(u, v)` is `C(u) + v d`. A plane through `o` with unit normal `n` holds
//! that point where `v = n·(o − C(u)) / n·d`, so where `n·d` is not zero the
//! section is `C` carried along `d` onto the plane, at `C`'s own parameter:
//! the image of `C` under an affine map. An affine map takes a line to a
//! line, a circle or an ellipse to an ellipse, a parabola to a parabola,
//! and a B-spline, polynomial or rational, to the B-spline of its control
//! points carried, on its own knots and weights. Where `n·d` is zero the
//! plane holds `d`, and the section is the rulings through `C`'s crossings
//! with the plane.
//!
//! Each section is cut to the extrusion's own window: to the stretches of
//! `C`'s domain where the height `v` stays within the sweep, found as `C`'s
//! crossings with the plane moved to each end of the sweep, each a root of
//! a closed form. A B-spline's crossings have none, and its section is
//! taken only where its control points' heights all lie within the sweep,
//! which holds every point's.
//!
//! The pair is left to the marcher, as before, wherever the closed form
//! does not hold cleanly: a profile this does not map (a helix, an offset,
//! a hyperbola, whose chart trace on the extrusion no planar curve states,
//! a reversed trim, a periodic B-spline), a B-spline whose heights straddle
//! an end of the sweep or that a plane holding `d` cuts, a profile lying in
//! a plane it is cut by, a crossing it only touches, a plane or a ruling
//! the extrusion meets within the sine the marcher takes as tangent, or a
//! section that does not land on both surfaces within the confusion
//! distance when checked.

use ogeom_core::Tolerances;
use ogeom_geom::{
    BSpline2d, BSplineCurve, Curve, Curve2d as _, Curve3d as _, LineCurve, PlanarCurve, Surface,
    SurfaceGeometry,
};
use ogeom_math::bspline::{self, Spline, Weighted};
use ogeom_math::{
    Axis, Axis2, Direction, Direction2, Frame, KnotVector, Point, Point2, Vector, Vector2,
};

use crate::march::TANGENTIAL_SINE;
use crate::section::{SectionCurve, exact_pcurve_over, inside_box, touches_box};

/// Where a plane meets an extrusion, exactly: `None` where the pair is not
/// a plane and an extrusion, or the closed form does not hold cleanly (see
/// the module's notes). An empty answer is a plane missing the extrusion's
/// window.
pub(crate) fn plane_extrusion(
    a: &SurfaceGeometry,
    b: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<Vec<SectionCurve>> {
    let (plane_surface, extrusion, plane_first) = match (a, b) {
        (SurfaceGeometry::Plane(_), SurfaceGeometry::Extrusion(e)) => (a, &**e, true),
        (SurfaceGeometry::Extrusion(e), SurfaceGeometry::Plane(_)) => (b, &**e, false),
        _ => return None,
    };
    let SurfaceGeometry::Plane(p) = plane_surface else {
        return None;
    };
    let plane = p.plane();
    let n = plane.normal().vector();
    let d = extrusion.direction().vector();
    let o = plane.origin();
    let (window, sweep) = extrusion.domain();
    let along = n.dot(d);
    let sections = if let Some(spline) = SplineProfile::of(extrusion.curve(), window, tol) {
        // A spline profile's crossings with a plane have no closed form, so
        // neither the rulings of a plane holding `d` nor the cuts at the
        // sweep's ends are taken here.
        if along.abs() < TANGENTIAL_SINE {
            return None;
        }
        spline.carried(&Carry { n, o, d, along }, sweep, tol)?
    } else {
        let (profile, window, whole_turn) = Profile::of(extrusion.curve(), window, tol)?;
        planar_sections(&profile, window, whole_turn, (n, o, d), extrusion, tol)?
    };
    finish(sections, plane_surface, extrusion, plane_first, tol)
}

/// The sections of an extrusion of a profile in closed form by the plane
/// through `o` with normal `n`, the extrusion swept along `d`.
fn planar_sections(
    profile: &Profile,
    window: (f64, f64),
    whole_turn: bool,
    (n, o, d): (Vector, Point, Vector),
    extrusion: &ogeom_geom::ExtrusionSurface,
    tol: Tolerances,
) -> Option<Vec<(Curve, PlanarCurve)>> {
    let (_, sweep) = extrusion.domain();
    let along = n.dot(d);
    let reach = sweep.0.abs().max(sweep.1.abs());
    Some(if along.abs() * reach <= tol.confusion() {
        rulings(profile, window, whole_turn, (n, o), extrusion, tol)?
    } else {
        // The extrusion's normal is square to `d`, so the sine between it
        // and `n` is at least `|n·d|`: past the tangency the marcher and
        // the contact handling answer, the two cross at an angle all
        // along every section.
        if along.abs() < TANGENTIAL_SINE {
            return None;
        }
        let carry = Carry { n, o, d, along };
        carried(profile, window, whole_turn, &carry, sweep, tol)?
    })
}

/// The pieces checked against both surfaces and stated as sections, each
/// trace on the surface asked first first.
fn finish(
    sections: Vec<(Curve, PlanarCurve)>,
    plane_surface: &SurfaceGeometry,
    extrusion: &ogeom_geom::ExtrusionSurface,
    plane_first: bool,
    tol: Tolerances,
) -> Option<Vec<SectionCurve>> {
    // Each piece checked against both surfaces, kept where it reaches the
    // plane's window, and cut to it where it is straight.
    let mut out = Vec::with_capacity(sections.len());
    for (curve, on_extrusion) in sections {
        let Some((curve, on_extrusion, on_plane)) =
            within_plane_window(curve, on_extrusion, plane_surface, tol)?
        else {
            continue;
        };
        lands_on_both(
            &curve,
            &on_plane,
            &on_extrusion,
            plane_surface,
            extrusion,
            tol,
        )?;
        let closed = curve.is_closed(tol);
        let (on_a, on_b) = if plane_first {
            (on_plane, on_extrusion)
        } else {
            (on_extrusion, on_plane)
        };
        out.push(SectionCurve {
            curve,
            on_a: Some(on_a),
            on_b: Some(on_b),
            tolerance: 0.0,
            exact: true,
            closed,
            tangential: false,
        });
    }
    Some(out)
}

/// The profile in the forms the map acts on, each at the profile's own
/// parameter `u`.
enum Profile {
    /// `p + u t`, `t` of unit length.
    Line { p: Point, t: Vector },
    /// `c + cos u · x + sin u · y`: a circle or an ellipse.
    Conic { c: Point, x: Vector, y: Vector },
    /// `apex + u² · x + u · y`.
    Parabola { apex: Point, x: Vector, y: Vector },
}

impl Profile {
    /// The profile of `curve` over `window`, the window, and whether that
    /// is a whole turn of a periodic profile. `None` for a profile this
    /// does not map.
    fn of(curve: &Curve, window: (f64, f64), tol: Tolerances) -> Option<(Self, (f64, f64), bool)> {
        let (basis, whole) = match curve {
            Curve::Trimmed(t) if !t.is_reversed() => (t.basis(), false),
            Curve::Trimmed(_) => return None,
            other => (other, other.is_periodic()),
        };
        let profile = match basis {
            Curve::Line(line) => {
                let axis = line.axis();
                Self::Line {
                    p: axis.location,
                    t: axis.direction.vector(),
                }
            }
            Curve::Circle(c) => {
                let circle = c.forward(tol).ok()?.circle();
                let frame = circle.frame();
                Self::Conic {
                    c: circle.centre(),
                    x: frame.x().vector() * circle.radius(),
                    y: frame.y().vector() * circle.radius(),
                }
            }
            Curve::Ellipse(e) => {
                let ellipse = e.forward(tol).ok()?.ellipse();
                let frame = ellipse.frame();
                Self::Conic {
                    c: ellipse.centre(),
                    x: frame.x().vector() * ellipse.major_radius(),
                    y: frame.y().vector() * ellipse.minor_radius(),
                }
            }
            Curve::Parabola(parabola) if !parabola.is_reversed() => {
                let parabola = parabola.parabola();
                let frame = parabola.frame();
                Self::Parabola {
                    apex: parabola.apex(),
                    x: frame.x().vector() / (4.0 * parabola.focal()),
                    y: frame.y().vector(),
                }
            }
            _ => return None,
        };
        Some((profile, window, whole))
    }

    /// The point at `u`.
    fn point(&self, u: f64) -> Point {
        match self {
            Self::Line { p, t } => *p + *t * u,
            Self::Conic { c, x, y } => *c + *x * u.cos() + *y * u.sin(),
            Self::Parabola { apex, x, y } => *apex + *x * (u * u) + *y * u,
        }
    }

    /// The first derivative at `u`.
    fn tangent(&self, u: f64) -> Vector {
        match self {
            Self::Line { t, .. } => *t,
            Self::Conic { x, y, .. } => *y * u.cos() - *x * u.sin(),
            Self::Parabola { x, y, .. } => *x * (2.0 * u) + *y,
        }
    }

    /// Where `n·C(u) = k` over `window`, in order. `None` where the
    /// profile lies within the confusion distance of that plane along a
    /// stretch, or touches it, or two crossings stand closer than the
    /// parametric confusion: there the plane is no clean crossing.
    fn crossings(
        &self,
        n: Vector,
        k: f64,
        (lo, hi): (f64, f64),
        tol: Tolerances,
    ) -> Option<Vec<f64>> {
        let near = tol.confusion();
        let mut found = Vec::new();
        match self {
            Self::Line { p, t } => {
                let (rate, at) = (n.dot(*t), n.dot(p.to_vector()) - k);
                if rate.abs() * (hi - lo) <= near {
                    if rate.mul_add(lo, at).abs() <= near {
                        return None;
                    }
                } else {
                    found.push(-at / rate);
                }
            }
            Self::Conic { c, x, y } => {
                // `A + B cos u + C sin u`: `R cos(u − φ) = −A`.
                let (a, b, s) = (n.dot(c.to_vector()) - k, n.dot(*x), n.dot(*y));
                let r = b.hypot(s);
                if (a.abs() - r).abs() <= near {
                    return None;
                }
                if a.abs() < r {
                    let phase = s.atan2(b);
                    let half = (-a / r).acos();
                    let tau = core::f64::consts::TAU;
                    for root in [phase - half, phase + half] {
                        let mut u = lo + (root - lo).rem_euclid(tau);
                        while u <= hi + tol.parametric() {
                            found.push(u);
                            u += tau;
                        }
                    }
                }
            }
            Self::Parabola { apex, x, y } => {
                let (a, b, c) = (n.dot(*x), n.dot(*y), n.dot(apex.to_vector()) - k);
                if a.abs() > 0.0 {
                    let turn = -b / (2.0 * a);
                    if (lo..=hi).contains(&turn)
                        && a.mul_add(turn * turn, b.mul_add(turn, c)).abs() <= near
                    {
                        return None;
                    }
                } else if b.abs() * (hi - lo) <= near && b.mul_add(lo, c).abs() <= near {
                    return None;
                }
                found.extend(ogeom_math::solve::quadratic_roots(a, b, c));
            }
        }
        let slack = tol.parametric();
        let mut kept: Vec<f64> = found
            .into_iter()
            .filter(|u| u.is_finite() && *u >= lo - slack && *u <= hi + slack)
            .map(|u| u.clamp(lo, hi))
            .collect();
        kept.sort_by(f64::total_cmp);
        kept.dedup_by(|x, y| (*x - *y).abs() <= slack);
        Some(kept)
    }
}

/// A clamped B-spline profile over the extrusion's window, polynomial or
/// rational: carried onto a plane, its control points carried with their
/// weights are the image's, since the carry is affine and the weights are
/// positive.
struct SplineProfile {
    knots: KnotVector,
    control: Vec<Weighted<Point>>,
}

impl SplineProfile {
    /// The profile of `curve` cut to `window`, where it is a B-spline, a
    /// forward trim of one, and not periodic.
    fn of(curve: &Curve, window: (f64, f64), tol: Tolerances) -> Option<Self> {
        let basis = match curve {
            Curve::Trimmed(t) if !t.is_reversed() => t.basis(),
            other => other,
        };
        let Curve::BSpline(spline) = basis else {
            return None;
        };
        if basis.is_periodic() || !spline.knots().is_clamped() {
            return None;
        }
        let mut piece: Spline<Weighted<Point>> =
            (spline.knots().clone(), spline.control_points().to_vec());
        let (lo, hi) = piece.0.domain();
        if window.0 > lo + tol.parametric() {
            piece = bspline::split(&piece.0, &piece.1, window.0, tol).ok()?.1;
        }
        if window.1 < hi - tol.parametric() {
            piece = bspline::split(&piece.0, &piece.1, window.1, tol).ok()?.0;
        }
        Some(Self {
            knots: piece.0,
            control: piece.1,
        })
    }

    /// The profile carried onto the plane, with its chart trace on the
    /// extrusion: `Some` of nothing where every control point's height,
    /// and so every point's, falls outside the sweep on one side; `None`
    /// where the heights straddle an end of the sweep, whose cut has no
    /// closed form.
    fn carried(
        &self,
        carry: &Carry,
        sweep: (f64, f64),
        tol: Tolerances,
    ) -> Option<Vec<(Curve, PlanarCurve)>> {
        // A point of the profile is a convex blend of its control points,
        // and so is its height of theirs.
        let heights: Vec<f64> = self
            .control
            .iter()
            .map(|w| carry.height(w.point()))
            .collect();
        if heights.iter().all(|&v| v < sweep.0) || heights.iter().all(|&v| v > sweep.1) {
            return Some(Vec::new());
        }
        if heights.iter().any(|&v| v < sweep.0 || v > sweep.1) {
            return None;
        }
        let image: Vec<Weighted<Point>> = self
            .control
            .iter()
            .map(|w| Weighted::new(carry.point(w.point()), w.weight, tol))
            .collect::<Result<_, _>>()
            .ok()?;
        let curve: Curve = BSplineCurve::rational(self.knots.clone(), image)
            .ok()?
            .into();
        let chart = self.chart(&heights, tol)?;
        Some(vec![(curve, chart.into())])
    }

    /// The chart trace `(u, v(u))` at the profile's own parameter, where
    /// `v` is the height blended as the profile's points are. Over each
    /// Bézier piece on `[a, b]` with weights `wₖ`, `v` is the blend of the
    /// heights with those weights, and `u` is `u W(u) / W(u)`, `u W(u)` one
    /// degree up with coefficients `((p + 1 − k) wₖ a + k wₖ₋₁ b) / (p + 1)`.
    /// The pieces join on their shared ends.
    fn chart(&self, heights: &[f64], tol: Tolerances) -> Option<BSpline2d> {
        let p = self.knots.degree();
        let lifted: Vec<Weighted<Point2>> = self
            .control
            .iter()
            .zip(heights)
            .map(|(w, &v)| Weighted {
                scaled: Point2::new(0.0, v * w.weight),
                weight: w.weight,
            })
            .collect();
        let pieces = bspline::to_bezier_segments(&self.knots, &lifted, tol).ok()?;
        let up = p + 1;
        #[allow(clippy::cast_precision_loss)]
        let share = |k: usize| k as f64 / up as f64;
        let mut knots: Vec<f64> = Vec::with_capacity(pieces.len() * up + 2);
        let mut control: Vec<Weighted<Point2>> = Vec::with_capacity(pieces.len() * up + 1);
        for (index, ((a, b), piece)) in pieces.iter().enumerate() {
            for k in 0..=up {
                if index > 0 && k == 0 {
                    continue;
                }
                let (before, here) = ((k > 0).then(|| piece[k - 1]), (k <= p).then(|| piece[k]));
                let blend = |f: fn(&Weighted<Point2>) -> f64| {
                    before.map_or(0.0, |w| share(k) * f(&w))
                        + here.map_or(0.0, |w| (1.0 - share(k)) * f(&w))
                };
                let weight = blend(|w| w.weight);
                let rise = blend(|w| w.scaled.y);
                let run = (1.0 - share(k)) * here.map_or(0.0, |w| w.weight) * a
                    + share(k) * before.map_or(0.0, |w| w.weight) * b;
                control.push(Weighted {
                    scaled: Point2::new(run, rise),
                    weight,
                });
            }
            let ends = if index == 0 { up + 1 } else { up };
            knots.extend(core::iter::repeat_n(*a, ends));
            if index == pieces.len() - 1 {
                knots.extend(core::iter::repeat_n(*b, up + 1));
            }
        }
        BSpline2d::rational(KnotVector::new(knots, up).ok()?, control).ok()
    }
}

/// The sections of an extrusion by a plane holding its direction: the
/// rulings through the profile's crossings with the plane, each over the
/// whole sweep and charted at its `u`. A crossing at an open end of the
/// profile is the extrusion's border lying in the plane, and is kept where
/// the extrusion leaves the plane there at an angle.
fn rulings(
    profile: &Profile,
    window: (f64, f64),
    whole_turn: bool,
    (n, o): (Vector, Point),
    extrusion: &ogeom_geom::ExtrusionSurface,
    tol: Tolerances,
) -> Option<Vec<(Curve, PlanarCurve)>> {
    let d = extrusion.direction();
    let (_, sweep) = extrusion.domain();
    let k = n.dot(o.to_vector());
    let mut roots = profile.crossings(n, k, window, tol)?;
    // An end of the profile on the plane to within the confusion distance
    // is a crossing there, whichever side rounding left it.
    for end in [window.0, window.1] {
        let gap = n.dot(profile.point(end).to_vector()) - k;
        if gap.abs() <= tol.confusion()
            && !roots.iter().any(|r| (r - end).abs() <= tol.parametric())
        {
            roots.push(end);
        }
    }
    roots.sort_by(f64::total_cmp);
    if whole_turn
        && roots.len() > 1
        && (roots[0] - window.0).abs() <= tol.parametric()
        && (roots[roots.len() - 1] - window.1).abs() <= tol.parametric()
    {
        roots.pop();
    }
    let mut out = Vec::with_capacity(roots.len());
    for u in roots {
        let tangent = profile.tangent(u);
        let across = tangent - d.vector() * tangent.dot(d.vector());
        let size = across.magnitude();
        if size <= tol.confusion() || n.dot(tangent).abs() < TANGENTIAL_SINE * size {
            return None;
        }
        let axis = Axis {
            location: profile.point(u),
            direction: d,
        };
        let line: Curve = LineCurve::over(axis, sweep.0, sweep.1).ok()?.into();
        let chart: PlanarCurve = ogeom_geom::Line2d::over(
            Axis2::new(Point2::new(u, 0.0), Direction2::Y),
            sweep.0,
            sweep.1,
        )
        .ok()?
        .into();
        out.push((line, chart));
    }
    Some(out)
}

/// The affine map carrying a point along `d` onto the plane through `o`
/// with normal `n`, where `along = n·d` is not zero.
struct Carry {
    n: Vector,
    o: Point,
    d: Vector,
    along: f64,
}

impl Carry {
    /// How far along `d` the point `p` is carried to reach the plane.
    fn height(&self, p: Point) -> f64 {
        self.n.dot(self.o - p) / self.along
    }

    /// The point `p` carried onto the plane.
    fn point(&self, p: Point) -> Point {
        p + self.d * self.height(p)
    }

    /// A vector carried as the difference of two points is.
    fn vector(&self, w: Vector) -> Vector {
        w - self.d * (self.n.dot(w) / self.along)
    }
}

/// The sections of an extrusion by a plane across its direction: the
/// profile carried onto the plane, over each stretch of the profile whose
/// height stays within the sweep, each with its chart trace on the
/// extrusion.
fn carried(
    profile: &Profile,
    window: (f64, f64),
    whole_turn: bool,
    carry: &Carry,
    sweep: (f64, f64),
    tol: Tolerances,
) -> Option<Vec<(Curve, PlanarCurve)>> {
    let level = |v: f64| carry.n.dot(carry.o.to_vector()) - v * carry.along;
    let mut cuts = vec![window.0, window.1];
    cuts.extend(profile.crossings(carry.n, level(sweep.0), window, tol)?);
    cuts.extend(profile.crossings(carry.n, level(sweep.1), window, tol)?);
    cuts.sort_by(f64::total_cmp);
    cuts.dedup_by(|x, y| (*x - *y).abs() <= tol.parametric());
    let mut stretches: Vec<(f64, f64)> = Vec::new();
    for pair in cuts.windows(2) {
        let (lo, hi) = (pair[0], pair[1]);
        let v = carry.height(profile.point(f64::midpoint(lo, hi)));
        if v < sweep.0 || v > sweep.1 {
            continue;
        }
        match stretches.last_mut() {
            Some(last) if (last.1 - lo).abs() <= tol.parametric() => last.1 = hi,
            _ => stretches.push((lo, hi)),
        }
    }
    // Round a whole turn, a stretch ending on its end runs on into one
    // starting on its start.
    let period = window.1 - window.0;
    if whole_turn && stretches.len() > 1 {
        let (first, last) = (stretches[0], stretches[stretches.len() - 1]);
        if (first.0 - window.0).abs() <= tol.parametric()
            && (last.1 - window.1).abs() <= tol.parametric()
        {
            stretches.pop();
            stretches[0] = (last.0, first.1 + period);
        }
    }
    let image = Image::of(profile, carry, tol)?;
    let mut out = Vec::with_capacity(stretches.len());
    for (lo, hi) in stretches {
        if hi - lo <= tol.parametric() {
            continue;
        }
        let whole = whole_turn && (hi - lo - period).abs() <= tol.parametric();
        out.push(image.over((lo, hi), whole, profile, carry, tol)?);
    }
    Some(out)
}

/// The profile carried onto the plane, in its own spelling.
enum Image {
    /// A line, at `scale` times the profile's parameter.
    Line { axis: Axis, scale: f64 },
    /// An ellipse or a circle, at the profile's parameter less `phase`.
    Conic { curve: Curve, phase: f64 },
    /// A parabola, at `scale` times the profile's parameter less `vertex`.
    Parabola {
        parabola: ogeom_math::Parabola,
        vertex: f64,
        scale: f64,
    },
}

impl Image {
    fn of(profile: &Profile, carry: &Carry, tol: Tolerances) -> Option<Self> {
        Some(match profile {
            Profile::Line { p, t } => {
                let run = carry.vector(*t);
                let scale = run.magnitude();
                Self::Line {
                    axis: Axis {
                        location: carry.point(*p),
                        direction: Direction::new(run, tol).ok()?,
                    },
                    scale,
                }
            }
            Profile::Conic { c, x, y } => {
                // Conjugate half-diameters `p`, `q`; turned by `phase` they
                // are the principal ones, `p cos φ + q sin φ` the major.
                let (p, q) = (carry.vector(*x), carry.vector(*y));
                let phase = 0.5 * (2.0 * p.dot(q)).atan2(p.dot(p) - q.dot(q));
                let (cos, sin) = (phase.cos(), phase.sin());
                let major = p * cos + q * sin;
                let minor = q * cos - p * sin;
                let (ra, rb) = (major.magnitude(), minor.magnitude());
                let xa = Direction::new(major, tol).ok()?;
                let ya = Direction::new(minor - xa.vector() * minor.dot(xa.vector()), tol).ok()?;
                let za = Direction::new(xa.vector().cross(ya.vector()), tol).ok()?;
                let frame = Frame::from_axes(carry.point(*c), xa, ya, za, tol).ok()?;
                let curve: Curve = if ra - rb <= tol.confusion() {
                    ogeom_geom::CircleCurve::new(ogeom_math::Circle::new(frame, ra, tol).ok()?)
                        .into()
                } else {
                    ogeom_geom::EllipseCurve::new(
                        ogeom_math::Ellipse::new(frame, ra, rb, tol).ok()?,
                    )
                    .into()
                };
                Self::Conic { curve, phase }
            }
            Profile::Parabola { apex, x, y } => {
                // `V + s² p + s q⊥` from the vertex, where the tangent is
                // square to the axis `p`.
                let (p, q) = (carry.vector(*x), carry.vector(*y));
                let vertex = -p.dot(q) / (2.0 * p.dot(p));
                let across = q - p * (p.dot(q) / p.dot(p));
                let scale = across.magnitude();
                let xa = Direction::new(p, tol).ok()?;
                let ya = Direction::new(across, tol).ok()?;
                let za = Direction::new(xa.vector().cross(ya.vector()), tol).ok()?;
                let origin = carry.point(*apex) + p * (vertex * vertex) + q * vertex;
                let frame = Frame::from_axes(origin, xa, ya, za, tol).ok()?;
                let focal = scale * scale / (4.0 * p.magnitude());
                Self::Parabola {
                    parabola: ogeom_math::Parabola::new(frame, focal, tol).ok()?,
                    vertex,
                    scale,
                }
            }
        })
    }

    /// The image over the profile's stretch `(lo, hi)`, the whole curve
    /// where `whole`, and its chart trace on the extrusion: `(u, v)` at
    /// the image's own parameter.
    fn over(
        &self,
        (lo, hi): (f64, f64),
        whole: bool,
        profile: &Profile,
        carry: &Carry,
        tol: Tolerances,
    ) -> Option<(Curve, PlanarCurve)> {
        let height = |u: f64| carry.height(profile.point(u));
        match self {
            Self::Line { axis, scale } => {
                let (t0, t1) = (lo * scale, hi * scale);
                let curve: Curve = LineCurve::over(*axis, t0, t1).ok()?.into();
                let chart = BSpline2d::new(
                    KnotVector::new(vec![t0, t0, t1, t1], 1).ok()?,
                    vec![Point2::new(lo, height(lo)), Point2::new(hi, height(hi))],
                    tol,
                )
                .ok()?;
                Some((curve, chart.into()))
            }
            Self::Conic { curve, phase } => {
                // `v = α + β cos u + γ sin u` at `u = t + phase`.
                let Profile::Conic { c, x, y } = profile else {
                    return None;
                };
                let alpha = carry.height(*c);
                let (beta, gamma) = (
                    -carry.n.dot(*x) / carry.along,
                    -carry.n.dot(*y) / carry.along,
                );
                let (cos, sin) = (phase.cos(), phase.sin());
                let range = if whole {
                    curve.domain()
                } else {
                    (lo - phase, hi - phase)
                };
                let chart = ogeom_geom::Trig2d::new(
                    Point2::new(*phase, alpha),
                    Vector2::new(1.0, 0.0),
                    Vector2::new(0.0, beta.mul_add(cos, gamma * sin)),
                    Vector2::new(0.0, gamma.mul_add(cos, -beta * sin)),
                    range,
                )
                .ok()?;
                let piece = if whole {
                    curve.clone()
                } else {
                    ogeom_geom::TrimmedCurve::new(curve.clone(), range.0, range.1, tol)
                        .ok()?
                        .into()
                };
                Some((piece, PlanarCurve::Trig(chart)))
            }
            Self::Parabola {
                parabola,
                vertex,
                scale,
            } => {
                let (t0, t1) = ((lo - vertex) * scale, (hi - vertex) * scale);
                let curve: Curve = ogeom_geom::ParabolaCurve::over(*parabola, t0, t1)
                    .ok()?
                    .into();
                // Linear in `u`, quadratic in `v`: one quadratic Bézier, its
                // middle control point where the end tangents meet.
                let rise = -carry.n.dot(profile.tangent(lo)) / carry.along;
                let start = Point2::new(lo, height(lo));
                let middle = Point2::new(
                    f64::midpoint(lo, hi),
                    rise.mul_add((hi - lo) / 2.0, start.y),
                );
                let end = Point2::new(hi, height(hi));
                let chart = BSpline2d::new(
                    KnotVector::new(vec![t0, t0, t0, t1, t1, t1], 2).ok()?,
                    vec![start, middle, end],
                    tol,
                )
                .ok()?;
                Some((curve, chart.into()))
            }
        }
    }
}

/// The piece and its chart traces cut to the plane's window: a straight
/// one to the stretch the window holds, any other kept whole where it
/// reaches the window and dropped where it does not, as the closed forms'
/// sections are. `Some(None)` for a piece wholly off the window, `None`
/// where a trace has no closed form.
#[allow(clippy::type_complexity, reason = "a piece and its two traces")]
fn within_plane_window(
    curve: Curve,
    on_extrusion: PlanarCurve,
    plane: &SurfaceGeometry,
    tol: Tolerances,
) -> Option<Option<(Curve, PlanarCurve, PlanarCurve)>> {
    let range = curve.domain();
    let on_plane = exact_pcurve_over(&curve, range, plane, tol)?;
    let Curve::Line(line) = &curve else {
        return Some(touches_box(&on_plane, plane, tol).then_some((curve, on_extrusion, on_plane)));
    };
    let Some((lo, hi)) = inside_box(&on_plane, plane) else {
        return Some(None);
    };
    let (lo, hi) = (lo.max(range.0), hi.min(range.1));
    if hi - lo <= tol.parametric() {
        return Some(None);
    }
    if lo <= range.0 && hi >= range.1 {
        return Some(Some((curve, on_extrusion, on_plane)));
    }
    let cut: Curve = LineCurve::over(line.axis(), lo, hi).ok()?.into();
    let on_extrusion: PlanarCurve = match on_extrusion {
        PlanarCurve::Line(l) => ogeom_geom::Line2d::over(l.axis(), lo, hi).ok()?.into(),
        other => {
            let (a, b) = (other.point_at(lo, tol).ok()?, other.point_at(hi, tol).ok()?);
            BSpline2d::new(
                KnotVector::new(vec![lo, lo, hi, hi], 1).ok()?,
                vec![a, b],
                tol,
            )
            .ok()?
            .into()
        }
    };
    let on_plane = exact_pcurve_over(&cut, (lo, hi), plane, tol)?;
    Some(Some((cut, on_extrusion, on_plane)))
}

/// Whether a section lands where it says: on the plane, and on each
/// surface at the chart point its trace gives, within the confusion
/// distance at stations along it. `None` where it does not.
fn lands_on_both(
    curve: &Curve,
    on_plane: &PlanarCurve,
    on_extrusion: &PlanarCurve,
    plane: &SurfaceGeometry,
    extrusion: &ogeom_geom::ExtrusionSurface,
    tol: Tolerances,
) -> Option<()> {
    const STATIONS: u32 = 16;
    let (lo, hi) = curve.domain();
    for i in 0..=STATIONS {
        let t = (hi - lo).mul_add(f64::from(i) / f64::from(STATIONS), lo);
        let point = curve.point_at(t, tol).ok()?;
        let flat = on_plane.point_at(t, tol).ok()?;
        let swept = on_extrusion.point_at(t, tol).ok()?;
        let on_plane_at = plane.point_at(flat.x, flat.y, tol).ok()?;
        let on_extrusion_at = extrusion.point_at(swept.x, swept.y, tol).ok()?;
        if point.distance(on_plane_at) > tol.confusion()
            || point.distance(on_extrusion_at) > tol.confusion()
        {
            return None;
        }
    }
    Some(())
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use ogeom_geom::{
        BSplineCurve, CircleCurve, Curve, Curve2d as _, Curve3d as _, EllipseCurve,
        ExtrusionSurface, LineCurve, ParabolaCurve, PlanarCurve, PlaneSurface, Surface as _,
        SurfaceGeometry, TrimmedCurve,
    };
    use ogeom_math::bspline::Weighted;
    use ogeom_math::{
        Axis, Circle, Direction, Ellipse, Frame, KnotVector, Parabola, Plane, Point, Vector,
    };
    use proptest::prelude::*;

    use crate::{IntersectOptions, SectionCurve, SurfaceIntersection, intersect_surfaces};

    const T: ogeom_core::Tolerances = ogeom_core::Tolerances::millimetres();

    fn direction(x: f64, y: f64, z: f64) -> Direction {
        Direction::new(Vector::new(x, y, z), T).unwrap()
    }

    /// The sections of `profile` swept along `along` over `sweep` by the
    /// plane through `at` square to `normal`, asked with the plane first;
    /// asked the other way round they must be the same curves, their
    /// traces swapped.
    fn cut(
        profile: Curve,
        along: Direction,
        sweep: (f64, f64),
        at: Point,
        normal: Direction,
    ) -> (Vec<SectionCurve>, SurfaceGeometry, SurfaceGeometry) {
        let extrusion: SurfaceGeometry = ExtrusionSurface::over(profile, along, sweep)
            .unwrap()
            .into();
        let plane: SurfaceGeometry = PlaneSurface::new(Plane::through(at, normal)).into();
        let ask = |a: &SurfaceGeometry, b: &SurfaceGeometry| -> Vec<SectionCurve> {
            match intersect_surfaces(a, b, IntersectOptions::default(), T).unwrap() {
                SurfaceIntersection::Along(found) => found,
                SurfaceIntersection::Apart => Vec::new(),
                other => panic!("not a section: {other:?}"),
            }
        };
        let first = ask(&plane, &extrusion);
        let second = ask(&extrusion, &plane);
        assert_eq!(first.len(), second.len());
        for (x, y) in first.iter().zip(&second) {
            assert_eq!(x.curve, y.curve);
            assert_eq!((&x.on_a, &x.on_b), (&y.on_b, &y.on_a));
        }
        (first, plane, extrusion)
    }

    /// The farthest a section stands off either surface, measured at its
    /// own points and at each trace's chart point, at stations along it.
    fn stray(section: &SectionCurve, plane: &SurfaceGeometry, extrusion: &SurfaceGeometry) -> f64 {
        let SurfaceGeometry::Plane(flat) = plane else {
            unreachable!()
        };
        let (lo, hi) = section.curve.domain();
        let (Some(on_plane), Some(on_extrusion)) = (&section.on_a, &section.on_b) else {
            panic!("a section without its traces");
        };
        (0..=64)
            .map(|i| {
                let t = (hi - lo).mul_add(f64::from(i) / 64.0, lo);
                let p = section.curve.point_at(t, T).unwrap();
                let a = on_plane.point_at(t, T).unwrap();
                let b = on_extrusion.point_at(t, T).unwrap();
                flat.plane()
                    .signed_distance_to(p)
                    .abs()
                    .max(plane.point_at(a.x, a.y, T).unwrap().distance(p))
                    .max(extrusion.point_at(b.x, b.y, T).unwrap().distance(p))
            })
            .fold(0.0, f64::max)
    }

    fn exact_on_both(found: &[SectionCurve], plane: &SurfaceGeometry, extrusion: &SurfaceGeometry) {
        for section in found {
            assert!(section.exact && section.tolerance == 0.0 && !section.tangential);
            let off = stray(section, plane, extrusion);
            assert!(off <= 1e-12, "a section stands {off} off");
        }
    }

    /// A line swept along a slant and cut level is the line where the
    /// sweep reaches the plane's height, over the profile's own stretch.
    #[test]
    fn a_swept_line_cut_level_is_the_line_at_that_height() {
        let profile: Curve = LineCurve::over(
            Axis {
                location: Point::ORIGIN,
                direction: Direction::X,
            },
            -2.0,
            3.0,
        )
        .unwrap()
        .into();
        let (found, plane, extrusion) = cut(
            profile,
            direction(0.0, 1.0, 1.0),
            (-5.0, 5.0),
            Point::new(0.0, 0.0, 1.0),
            Direction::Z,
        );
        exact_on_both(&found, &plane, &extrusion);
        let [section] = found.as_slice() else {
            panic!("{} sections", found.len());
        };
        assert!(matches!(section.curve, Curve::Line(_)));
        let (lo, hi) = section.curve.domain();
        let ends = [lo, hi].map(|t| section.curve.point_at(t, T).unwrap());
        assert!(ends[0].distance(Point::new(-2.0, 1.0, 1.0)) < 1e-12);
        assert!(ends[1].distance(Point::new(3.0, 1.0, 1.0)) < 1e-12);
    }

    /// A drum's circle swept along its axis and cut at a tilt `θ` is the
    /// closed ellipse of radii `R / cos θ` and `R` about the axis.
    #[test]
    fn a_swept_circle_cut_at_a_tilt_is_an_ellipse_of_the_tilt() {
        let (radius, tilt) = (2.0, 0.6_f64);
        let centre = Point::new(1.0, 0.0, 0.0);
        let circle = Circle::new(
            Frame::new(centre, Direction::Z, Direction::X, T).unwrap(),
            radius,
            T,
        )
        .unwrap();
        let (found, plane, extrusion) = cut(
            CircleCurve::new(circle).into(),
            Direction::Z,
            (-10.0, 10.0),
            centre,
            direction(tilt.sin(), 0.0, tilt.cos()),
        );
        exact_on_both(&found, &plane, &extrusion);
        let [section] = found.as_slice() else {
            panic!("{} sections", found.len());
        };
        let Curve::Ellipse(e) = &section.curve else {
            panic!("not an ellipse: {:?}", section.curve);
        };
        let e = e.ellipse();
        assert!(section.closed);
        assert!((e.major_radius() - radius / tilt.cos()).abs() < 1e-12);
        assert!((e.minor_radius() - radius).abs() < 1e-12);
        assert!(e.centre().distance(centre) < 1e-12);
    }

    /// An ellipse swept along a slant and cut by a plane parallel to its
    /// own is the same ellipse, moved along the slant to the plane.
    #[test]
    fn a_swept_ellipse_cut_parallel_to_it_is_the_ellipse_moved() {
        let ellipse = Ellipse::new(Frame::WORLD, 3.0, 1.0, T).unwrap();
        let (found, plane, extrusion) = cut(
            EllipseCurve::new(ellipse).into(),
            direction(1.0, 0.0, 1.0),
            (-1.0, 5.0),
            Point::new(0.0, 0.0, 2.0),
            Direction::Z,
        );
        exact_on_both(&found, &plane, &extrusion);
        let [section] = found.as_slice() else {
            panic!("{} sections", found.len());
        };
        let Curve::Ellipse(e) = &section.curve else {
            panic!("not an ellipse: {:?}", section.curve);
        };
        let e = e.ellipse();
        assert!((e.major_radius() - 3.0).abs() < 1e-12);
        assert!((e.minor_radius() - 1.0).abs() < 1e-12);
        assert!(e.centre().distance(Point::new(2.0, 0.0, 2.0)) < 1e-12);
    }

    /// A parabola of focal length `f` swept square to its plane and cut at
    /// a tilt `θ` about its axis's normal is the parabola of focal length
    /// `f cos θ` from the same apex: its axis stretched by `1 / cos θ`.
    #[test]
    fn a_swept_parabola_cut_at_a_tilt_has_its_focal_length_shortened() {
        let (focal, tilt) = (0.5, 0.4_f64);
        let parabola = Parabola::new(Frame::WORLD, focal, T).unwrap();
        let (found, plane, extrusion) = cut(
            ParabolaCurve::over(parabola, -2.0, 2.0).unwrap().into(),
            Direction::Z,
            (-10.0, 10.0),
            Point::ORIGIN,
            direction(-tilt.sin(), 0.0, tilt.cos()),
        );
        exact_on_both(&found, &plane, &extrusion);
        let [section] = found.as_slice() else {
            panic!("{} sections", found.len());
        };
        let Curve::Parabola(p) = &section.curve else {
            panic!("not a parabola: {:?}", section.curve);
        };
        assert!((p.parabola().focal() - focal * tilt.cos()).abs() < 1e-12);
        assert!(p.parabola().apex().distance(Point::ORIGIN) < 1e-12);
    }

    /// A plane holding the sweep's direction cuts it in the rulings
    /// through the profile's crossings with it, over the whole sweep, each
    /// charted at the angle it stands at.
    #[test]
    fn a_plane_along_the_sweep_cuts_it_in_rulings() {
        let circle = Circle::new(Frame::WORLD, 2.0, T).unwrap();
        let (found, plane, extrusion) = cut(
            CircleCurve::new(circle).into(),
            Direction::Z,
            (0.0, 3.0),
            Point::new(1.0, 0.0, 0.0),
            Direction::X,
        );
        exact_on_both(&found, &plane, &extrusion);
        assert_eq!(found.len(), 2);
        let third = core::f64::consts::FRAC_PI_3;
        for (section, (y, u)) in found
            .iter()
            .zip([(3.0_f64.sqrt(), third), (-(3.0_f64.sqrt()), 5.0 * third)])
        {
            let Curve::Line(line) = &section.curve else {
                panic!("not a ruling: {:?}", section.curve);
            };
            assert!(line.axis().direction.is_parallel(Direction::Z, T));
            let (lo, hi) = section.curve.domain();
            let ends = [lo, hi].map(|t| section.curve.point_at(t, T).unwrap());
            assert!(ends[0].distance(Point::new(1.0, y, 0.0)) < 1e-12);
            assert!(ends[1].distance(Point::new(1.0, y, 3.0)) < 1e-12);
            let Some(PlanarCurve::Line(chart)) = &section.on_b else {
                panic!("a ruling charted as {:?}", section.on_b);
            };
            assert!((chart.axis().location.x - u).abs() < 1e-12);
        }
    }

    /// A circle swept over `[0, 1]` and cut by `z = 1/2 + x`: the height
    /// `1/2 + cos u` leaves the sweep, and what stays is the two arcs over
    /// `cos u` in `[−1/2, 1/2]`, each ending on the sweep's two rims.
    #[test]
    fn a_section_leaving_the_sweep_is_cut_to_it() {
        let circle = Circle::new(Frame::WORLD, 1.0, T).unwrap();
        let (found, plane, extrusion) = cut(
            CircleCurve::new(circle).into(),
            Direction::Z,
            (0.0, 1.0),
            Point::new(0.0, 0.0, 0.5),
            direction(-1.0, 0.0, 1.0),
        );
        exact_on_both(&found, &plane, &extrusion);
        assert_eq!(found.len(), 2);
        for section in &found {
            assert!(!section.closed);
            let (lo, hi) = section.curve.domain();
            let mut heights = [lo, hi].map(|t| section.curve.point_at(t, T).unwrap().z);
            heights.sort_by(f64::total_cmp);
            assert!(heights[0].abs() < 1e-12 && (heights[1] - 1.0).abs() < 1e-12);
        }
    }

    /// Cut by `z = 1/2 + 4y/5`, the circle swept over `[0, 1]` keeps
    /// `sin u` in `[−5/8, 5/8]`: two arcs, one through `(1, 0, 1/2)` across
    /// the profile's own start, stated as one section rather than two
    /// halves.
    #[test]
    fn an_arc_across_the_profile_s_start_is_one_section() {
        let circle = Circle::new(Frame::WORLD, 1.0, T).unwrap();
        let (found, plane, extrusion) = cut(
            CircleCurve::new(circle).into(),
            Direction::Z,
            (0.0, 1.0),
            Point::new(0.0, 0.0, 0.5),
            direction(0.0, -0.8, 1.0),
        );
        exact_on_both(&found, &plane, &extrusion);
        assert_eq!(found.len(), 2);
        let across = found
            .iter()
            .filter(|section| {
                let (lo, hi) = section.curve.domain();
                let at = |t: f64| section.curve.point_at(t, T).unwrap();
                let mid = at(f64::midpoint(lo, hi));
                let ends = [at(lo).z, at(hi).z];
                assert!(ends.iter().any(|z| z.abs() < 1e-12));
                assert!(ends.iter().any(|z| (z - 1.0).abs() < 1e-12));
                mid.distance(Point::new(1.0, 0.0, 0.5)) < 1e-12
            })
            .count();
        assert_eq!(across, 1);
    }

    /// A circle's arc spelled as a rational quadratic B-spline, swept along
    /// its axis and cut at a tilt `θ`: every point of the section lies on
    /// the drum and in the plane, on the ellipse of radii `R / cos θ` and
    /// `R`, and the section at the profile's parameter `t` is the arc's
    /// point there lifted to the plane.
    #[test]
    fn a_swept_arc_spline_cut_at_a_tilt_is_the_arc_lifted_to_the_plane() {
        let (radius, tilt) = (2.0, 0.6_f64);
        let centre = Point::new(1.0, 0.0, 0.0);
        let circle = Circle::new(
            Frame::new(centre, Direction::Z, Direction::X, T).unwrap(),
            radius,
            T,
        )
        .unwrap();
        let arc: Curve = Curve::from(CircleCurve::new(circle))
            .to_bspline_over((0.3, 2.5), T)
            .unwrap()
            .into();
        let (found, plane, extrusion) = cut(
            arc.clone(),
            Direction::Z,
            (-10.0, 10.0),
            centre,
            direction(tilt.sin(), 0.0, tilt.cos()),
        );
        exact_on_both(&found, &plane, &extrusion);
        let [section] = found.as_slice() else {
            panic!("{} sections", found.len());
        };
        assert!(matches!(section.curve, Curve::BSpline(_)));
        assert_eq!(section.curve.domain(), arc.domain());
        let (lo, hi) = arc.domain();
        for i in 0..=32 {
            let t = (hi - lo).mul_add(f64::from(i) / 32.0, lo);
            let p = section.curve.point_at(t, T).unwrap();
            let q = arc.point_at(t, T).unwrap();
            let lifted = Point::new(q.x, q.y, -tilt.tan() * (q.x - 1.0));
            assert!(p.distance(lifted) < 1e-12);
            let (x, y) = ((p.x - 1.0).mul_add(tilt.cos(), -p.z * tilt.sin()), p.y);
            let on = (x * tilt.cos() / radius).powi(2) + (y / radius).powi(2);
            assert!((on - 1.0).abs() < 1e-12);
        }
    }

    /// A polynomial cubic winding in space, swept along a slant `d` and cut
    /// by `z = 1`: the section at `t` is `C(t) + d (1 − C_z(t)) / d_z`, and
    /// its trace on the extrusion is `(t, (1 − C_z(t)) / d_z)`. Swept short
    /// of the plane it is apart, and swept to straddle it the pair is left
    /// to the marcher.
    #[test]
    fn a_swept_cubic_cut_level_is_the_cubic_carried_along_the_sweep() {
        let control = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 2.0, 0.5),
            Point::new(2.0, -1.0, -0.4),
            Point::new(3.0, 1.0, 0.3),
            Point::new(4.0, 0.0, 0.0),
        ];
        let cubic: Curve =
            BSplineCurve::new(KnotVector::clamped_uniform(3, 5).unwrap(), control, T)
                .unwrap()
                .into();
        let along = direction(0.0, 1.0, 1.0);
        let (found, plane, extrusion) = cut(
            cubic.clone(),
            along,
            (-5.0, 5.0),
            Point::new(0.0, 0.0, 1.0),
            Direction::Z,
        );
        exact_on_both(&found, &plane, &extrusion);
        let [section] = found.as_slice() else {
            panic!("{} sections", found.len());
        };
        let Some(trace) = &section.on_b else {
            panic!("a section without its trace");
        };
        let d = along.vector();
        let (lo, hi) = cubic.domain();
        for i in 0..=32 {
            let t = (hi - lo).mul_add(f64::from(i) / 32.0, lo);
            let q = cubic.point_at(t, T).unwrap();
            let height = (1.0 - q.z) / d.z;
            assert!(
                section
                    .curve
                    .point_at(t, T)
                    .unwrap()
                    .distance(q + d * height)
                    < 1e-12
            );
            let uv = trace.point_at(t, T).unwrap();
            assert!((uv.x - t).abs() < 1e-12 && (uv.y - height).abs() < 1e-12);
        }
        let plane: SurfaceGeometry =
            PlaneSurface::new(Plane::through(Point::new(0.0, 0.0, 1.0), Direction::Z)).into();
        let short: SurfaceGeometry = ExtrusionSurface::over(cubic.clone(), along, (-5.0, 0.5))
            .unwrap()
            .into();
        assert_eq!(
            super::plane_extrusion(&plane, &short, T).map(|f| f.len()),
            Some(0)
        );
        let straddling: SurfaceGeometry = ExtrusionSurface::over(cubic, along, (-5.0, 1.0))
            .unwrap()
            .into();
        assert!(super::plane_extrusion(&plane, &straddling, T).is_none());
    }

    /// The profiles the property test sweeps: a line, a circle, an
    /// ellipse, a parabola and a rational cubic, whole or trimmed.
    fn profile() -> impl Strategy<Value = Curve> {
        let point = || (-2.0..2.0f64, -2.0..2.0f64, -2.0..2.0f64);
        let frame = (point(), -3.0..3.0f64, -3.0..3.0f64).prop_map(|((x, y, z), a, b)| {
            let axis = direction(a.cos() * b.sin(), a.sin() * b.sin(), b.cos());
            let lean = Vector::new(1.0, 0.3, -0.2);
            let x_axis = Direction::new(lean - axis.vector() * lean.dot(axis.vector()), T).unwrap();
            Frame::new(Point::new(x, y, z), axis, x_axis, T).unwrap()
        });
        let kind = 0..5u32;
        let trim = proptest::option::of((0.05..0.45f64, 0.55..0.95f64));
        (kind, frame, 0.3..2.0f64, 0.2..1.0f64, trim).prop_map(
            |(kind, frame, size, ratio, trim)| {
                let curve: Curve = match kind {
                    0 => LineCurve::over(
                        Axis {
                            location: frame.origin(),
                            direction: frame.x(),
                        },
                        -size,
                        size * 2.0,
                    )
                    .unwrap()
                    .into(),
                    1 => CircleCurve::new(Circle::new(frame, size, T).unwrap()).into(),
                    2 => EllipseCurve::new(Ellipse::new(frame, size, size * ratio, T).unwrap())
                        .into(),
                    3 => ParabolaCurve::over(
                        Parabola::new(frame, size * ratio, T).unwrap(),
                        -size,
                        size,
                    )
                    .unwrap()
                    .into(),
                    // A rational cubic winding off its frame's plane.
                    _ => {
                        let control = (0..6u32)
                            .map(|i| {
                                let k = f64::from(i);
                                let p = frame.origin()
                                    + frame.x().vector() * (size * (k - 2.5))
                                    + frame.y().vector() * (size * ratio * (k * 1.7).sin())
                                    + frame.z().vector() * (size * 0.3 * (k * 0.9).cos());
                                Weighted::new(p, 1.0 + ratio * f64::from(i % 2), T).unwrap()
                            })
                            .collect();
                        BSplineCurve::rational(KnotVector::clamped_uniform(3, 6).unwrap(), control)
                            .unwrap()
                            .into()
                    }
                };
                match trim {
                    Some((a, b)) => {
                        let (lo, hi) = curve.domain();
                        let at = |f: f64| (hi - lo).mul_add(f, lo);
                        TrimmedCurve::new(curve, at(a), at(b), T).unwrap().into()
                    }
                    None => curve,
                }
            },
        )
    }

    proptest! {
        /// Whatever the profile, sweep and plane across it, every section
        /// lies on both surfaces to within 1e-12, at its own points and at
        /// its traces' chart points; and every point of the profile whose
        /// height to the plane falls inside the sweep is on a section.
        #[test]
        fn a_carried_section_lies_on_both_and_misses_nothing(
            profile in profile(),
            (a, b) in (-3.0..3.0f64, 0.0..3.0f64),
            (c, e) in (-3.0..3.0f64, 0.0..3.0f64),
            at in (-1.0..1.0f64, -1.0..1.0f64, -1.0..1.0f64),
            sweep in (-4.0..-0.5f64, 0.5..4.0f64),
        ) {
            let along = direction(a.cos() * b.sin(), a.sin() * b.sin(), b.cos());
            let normal = direction(c.cos() * e.sin(), c.sin() * e.sin(), e.cos());
            prop_assume!(normal.vector().dot(along.vector()).abs() >= 0.2);
            let at = Point::new(at.0, at.1, at.2);
            let extrusion = ExtrusionSurface::over(profile.clone(), along, sweep).unwrap();
            let plane_surface: SurfaceGeometry = PlaneSurface::new(Plane::through(at, normal)).into();
            let Some(found) =
                super::plane_extrusion(&plane_surface, &extrusion.clone().into(), T)
            else {
                // Declined: a touch, or a section meeting the extrusion
                // all but tangentially. Neither is what this measures.
                return Ok(());
            };
            let extrusion_surface: SurfaceGeometry = extrusion.into();
            for section in &found {
                prop_assert!(stray(section, &plane_surface, &extrusion_surface) <= 1e-12);
            }
            // Every point whose height falls inside the sweep is charted
            // by a section: its `u` lies within one's trace.
            let spans: Vec<(f64, f64)> = found
                .iter()
                .map(|s| {
                    let trace = s.on_b.as_ref().unwrap();
                    let (lo, hi) = trace.domain();
                    let (p, q) = (trace.point_at(lo, T).unwrap().x, trace.point_at(hi, T).unwrap().x);
                    (p.min(q), p.max(q))
                })
                .collect();
            let (lo, hi) = profile.domain();
            let period = core::f64::consts::TAU;
            for i in 0..=200 {
                let u = (hi - lo).mul_add(f64::from(i) / 200.0, lo);
                let p = profile.point_at(u, T).unwrap();
                let v = normal.vector().dot(at - p) / normal.vector().dot(along.vector());
                if v <= sweep.0 + 1e-6 || v >= sweep.1 - 1e-6 {
                    continue;
                }
                let covered = spans.iter().any(|&(s0, s1)| {
                    [0.0, period, -period].iter().any(|k| u + k >= s0 - 1e-9 && u + k <= s1 + 1e-9)
                });
                prop_assert!(covered, "u {u} at height {v} is on no section");
            }
        }
    }
}
