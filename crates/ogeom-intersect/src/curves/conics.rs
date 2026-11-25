//! Closed forms for space conics: skew and coplanar pairs, a line against a conic, and same supports.

use super::*;

/// The closed-form answers for a pair of space curves, or `None` where there
/// is none and the sampling path is the honest route.
pub(super) fn analytic_3d(
    a: &Curve,
    b: &Curve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Option<CurveIntersection<Point>> {
    match (a, b) {
        (Curve::Line(x), Curve::Line(y)) => Some(line_line_3d(x, y, options, tol)),
        (Curve::Circle(x), Curve::Circle(y)) => same_circle_3d(x, y, tol)
            .or_else(|| coplanar_circles_3d(x, y, options, tol))
            .or_else(|| skew_conics_3d(a, b, options, tol)),
        (Curve::Ellipse(x), Curve::Ellipse(y)) => {
            same_ellipse_3d(x, y, tol).or_else(|| skew_conics_3d(a, b, options, tol))
        }
        (Curve::Circle(_), Curve::Ellipse(_)) | (Curve::Ellipse(_), Curve::Circle(_)) => {
            skew_conics_3d(a, b, options, tol)
        }
        (Curve::Line(x), Curve::Circle(_) | Curve::Ellipse(_)) => line_conic_3d(x, b, options, tol)
            .map(|mut found| {
                for c in &mut found.crossings {
                    core::mem::swap(&mut c.on_a, &mut c.on_b);
                }
                found.crossings.sort_by(|x, y| x.on_a.total_cmp(&y.on_a));
                found
            }),
        (Curve::Circle(_) | Curve::Ellipse(_), Curve::Line(y)) => line_conic_3d(y, a, options, tol),
        _ => None,
    }
}

/// A circle or ellipse, as its centre, the two semi-axis vectors its
/// parameter turns between, and its plane's normal; `None` for any other
/// curve, or one whose parameter runs backwards.
fn conic_of(
    curve: &Curve,
) -> Option<(
    Point,
    ogeom_math::Vector,
    ogeom_math::Vector,
    ogeom_math::Vector,
)> {
    let (centre, [x, y, z], [ra, rb]) = conic_axes(curve)?;
    Some((centre, x * ra, y * rb, z))
}

/// A circle or ellipse as its centre, its frame's three unit axes and the
/// radii along the first two; `None` as for [`conic_of`].
fn conic_axes(curve: &Curve) -> Option<(Point, [ogeom_math::Vector; 3], [f64; 2])> {
    let (frame, radii) = match curve {
        Curve::Circle(c) if !c.is_reversed() => {
            let circle = c.circle();
            (circle.frame(), [circle.radius(); 2])
        }
        Curve::Ellipse(e) if !e.is_reversed() => {
            let ellipse = e.ellipse();
            (
                ellipse.frame(),
                [ellipse.major_radius(), ellipse.minor_radius()],
            )
        }
        _ => return None,
    };
    Some((
        frame.origin(),
        [frame.x().vector(), frame.y().vector(), frame.z().vector()],
        radii,
    ))
}

/// Two circles or ellipses in planes that are not one plane, in closed
/// form: where the first meets the second's plane (`alpha cos t + beta sin
/// t + gamma = 0`, at most twice), kept where the second curve passes
/// within the gap. Parallel planes apart share no point. `None` for a
/// coplanar pair, which the sampling path answers.
fn skew_conics_3d(
    a: &Curve,
    b: &Curve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Option<CurveIntersection<Point>> {
    let (ca, ua, va, na) = conic_of(a)?;
    let (cb, _, _, nb) = conic_of(b)?;
    if na.cross(nb).magnitude() <= tol.angular() {
        return ((ca - cb).dot(nb).abs() > options.gap.max(tol.confusion()))
            .then(CurveIntersection::empty);
    }
    let (alpha, beta, gamma) = (nb.dot(ua), nb.dot(va), nb.dot(ca - cb));
    let size = alpha.hypot(beta);
    // The first conic stands within the gap of the other's plane over a
    // stretch as wide as the gap is against how far it rises from that
    // plane. Planes all but parallel (a rim fitted to one facet group, a
    // section through another) make that stretch wide, and the curves can
    // pass within the gap anywhere along it, not only where one crosses the
    // other's plane; the sampling path measures such a pass.
    if options.gap.max(tol.confusion()) > size * 1e-3 {
        return None;
    }
    let mut crossings: Vec<Crossing<Point>> = Vec::new();
    if size > 0.0 && gamma.abs() <= size * (1.0 + 1e-12) {
        let phase = beta.atan2(alpha);
        let turn = (-gamma / size).clamp(-1.0, 1.0).acos();
        let (a_lo, a_hi) = a.domain();
        let (b_lo, b_hi) = b.domain();
        let tau = core::f64::consts::TAU;
        let into = |t: f64, lo: f64| lo + (t - lo).rem_euclid(tau);
        let roots = if turn <= 1e-12 {
            vec![phase]
        } else {
            vec![phase - turn, phase + turn]
        };
        for root in roots {
            let t = into(root, a_lo);
            if t > a_hi + tol.parametric() {
                continue;
            }
            let point = a.point_at(t, tol).ok()?;
            // The closed form reads the curve as its centre and semi-axes;
            // a point that does not land on the other plane means it read
            // it wrong, and the sampling path answers instead.
            if (point - cb).dot(nb).abs() > tol.confusion() * 10.0 {
                return None;
            }
            let s = match b {
                Curve::Circle(c) => {
                    ogeom_math::elementary::circle_parameter(&c.circle(), point, tol)
                }
                Curve::Ellipse(e) => {
                    ogeom_math::elementary::ellipse_parameter(&e.ellipse(), point, tol)
                }
                _ => return None,
            };
            let Ok(s) = s else {
                continue;
            };
            let s = into(s, b_lo);
            if s > b_hi + tol.parametric() {
                continue;
            }
            let gap = b.point_at(s, tol).ok()?.distance(point);
            if gap > options.gap {
                continue;
            }
            crossings.push(Crossing {
                on_a: t,
                on_b: s,
                point,
                gap,
                reach: 0.0,
            });
        }
    }
    crossings.sort_by(|x, y| x.on_a.total_cmp(&y.on_a));
    Some(CurveIntersection {
        crossings,
        overlaps: Vec::new(),
    })
}

/// A circle or ellipse against a line, in closed form, the conic first.
///
/// A line in the conic's plane meets it where a quadratic along the line
/// vanishes, at most twice; a line through the plane meets it at most where
/// it pierces the plane. Kept where the two pass within the gap. A line
/// running nearly along the plane without lying in it is left to the
/// sampling path. So is a crossing nearly tangent where the line's chord
/// through the conic stands within the gap of it at its middle: the two
/// are one touch there, and the sampling path measures how far it reaches.
/// A chord deeper than the gap crosses twice however shallowly, as the
/// conic's surface and the line's do, and its roots are the crossings.
fn line_conic_3d(
    line: &ogeom_geom::LineCurve,
    conic: &Curve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Option<CurveIntersection<Point>> {
    let (centre, u, v, normal) = conic_of(conic)?;
    let (origin, along) = (line.axis().location, line.axis().direction.vector());
    let (a_len, b_len) = (u.magnitude(), v.magnitude());
    if a_len <= tol.confusion() || b_len <= tol.confusion() {
        return None;
    }
    let (ux, vy) = (u / a_len, v / b_len);
    let lean = along.dot(normal);
    let height = (origin - centre).dot(normal);
    let mut ts: Vec<f64> = Vec::new();
    let mut deep = false;
    if lean.abs() <= tol.angular() {
        if height.abs() > options.gap.max(tol.confusion()) {
            return Some(CurveIntersection::empty());
        }
        // (x0 + t dx)^2 / a^2 + (y0 + t dy)^2 / b^2 = 1, in the plane.
        let (x0, y0) = ((origin - centre).dot(ux), (origin - centre).dot(vy));
        let (dx, dy) = (along.dot(ux), along.dot(vy));
        let qa = dx * dx / (a_len * a_len) + dy * dy / (b_len * b_len);
        let qb = 2.0 * (x0 * dx / (a_len * a_len) + y0 * dy / (b_len * b_len));
        let qc = x0 * x0 / (a_len * a_len) + y0 * y0 / (b_len * b_len) - 1.0;
        let disc = qb.mul_add(qb, -4.0 * qa * qc);
        if qa <= 0.0 {
            return None;
        }
        if disc < 0.0 {
            let t = -qb / (2.0 * qa);
            let p = origin + along * t;
            let foot = conic_parameter(conic, p, tol)?;
            let gap = conic.point_at(foot, tol).ok()?.distance(p);
            return (gap > options.gap).then(CurveIntersection::empty);
        }
        let middle = origin + along * (-qb / (2.0 * qa));
        let foot = conic_parameter(conic, middle, tol)?;
        deep = conic.point_at(foot, tol).ok()?.distance(middle) > options.gap;
        let root = disc.sqrt();
        ts.push((-qb - root) / (2.0 * qa));
        ts.push((-qb + root) / (2.0 * qa));
    } else if lean.abs() >= 0.1 {
        ts.push(-height / lean);
    } else {
        return None;
    }
    let (lo, hi) = line.domain();
    let (c_lo, c_hi) = conic.domain();
    let tau = core::f64::consts::TAU;
    let mut crossings: Vec<Crossing<Point>> = Vec::new();
    for t in ts {
        if t < lo - tol.parametric() || t > hi + tol.parametric() {
            continue;
        }
        let point = origin + along * t;
        let s = conic_parameter(conic, point, tol)?;
        let s = c_lo + (s - c_lo).rem_euclid(tau);
        if s > c_hi + tol.parametric() {
            continue;
        }
        let on_conic = conic.point_at(s, tol).ok()?;
        let gap = on_conic.distance(point);
        if gap > options.gap {
            continue;
        }
        let tangent = conic.d1_at(s, tol).ok()?;
        if !deep && tangent.cross(along).magnitude() <= 1e-3 * tangent.magnitude() {
            return None;
        }
        crossings.push(Crossing {
            on_a: s,
            on_b: t,
            point: on_conic,
            gap,
            reach: 0.0,
        });
    }
    crossings.sort_by(|x, y| x.on_a.total_cmp(&y.on_a));
    crossings.dedup_by(|x, y| (x.on_a - y.on_a).abs() <= tol.parametric());
    Some(CurveIntersection {
        crossings,
        overlaps: Vec::new(),
    })
}

/// Where a point lies along a circle or ellipse, by projection.
fn conic_parameter(conic: &Curve, point: Point, tol: Tolerances) -> Option<f64> {
    match conic {
        Curve::Circle(c) => ogeom_math::elementary::circle_parameter(&c.circle(), point, tol).ok(),
        Curve::Ellipse(e) => {
            ogeom_math::elementary::ellipse_parameter(&e.ellipse(), point, tol).ok()
        }
        _ => None,
    }
}

/// A line touching a circle or ellipse in its plane, in closed form, with
/// the answer the sampling path gives for it: one crossing at the touch,
/// owning the valley around it as its reach.
///
/// In the conic's frame, scaled by its radii, the conic is `F = x^2 + y^2 -
/// 1 = 0`, and along the line `F` is a quadratic `qa (t - t_v)^2 + F_min`
/// whose roots are the crossings. Its value `f` and slope `g` are evaluated
/// at the computed vertex, each with a rounding bound (see
/// [`ScaledConic::at`]). Then `F_min` lies in `[f - err_f - G, f + err_f]`,
/// where `G = (|g| + err_g)^2 / 4 qa` is how far below `f` the true vertex
/// can sit. A line whose interval holds zero has a double root to the
/// root's own precision: it touches the conic, or crosses or misses it by
/// less than rounding can tell, and every true root lies within `(|g| +
/// err_g) / 2 qa + sqrt((err_f + G - f) / qa)` of the computed vertex along
/// the line. That is the tangent point, answered as one contact at the
/// measured gap.
///
/// Its reach is the sampling path's own measure, walked the same way along
/// the first curve with the exact distance to the second in place of the
/// sampled one. A tangent line stays within the gap of a convex conic only
/// around the touch (the conic's other tangent parallel to the line is a
/// full width away, more than the gap), so that valley is the only contact,
/// and the roots inside it are what it speaks for. A valley no longer than
/// eight gaps owns no reach there either, and is one crossing as long as
/// any two points in it are one crossing to the sampling path.
///
/// A touch outside a trim leaves the curves apart when a lower bound on
/// their distance clears the gap: the line's window, short of the vertex,
/// is nearest the conic at its end nearer the vertex, and no nearer to it
/// than the smaller radius times that end's scaled distance from the unit
/// circle; the conic's window is no nearer the line than its signed
/// in-plane distance from the line's shadow on the plane, a sinusoid in the
/// conic's parameter whose least magnitude on the window is at an end or a
/// stationary point unless it changes sign.
///
/// Everything else is left to the sampling path: an interval clear of zero
/// (a plain miss or a secant, which [`line_conic_3d`] answered), a touch
/// outside a trim with the curves within the gap, a valley holding two
/// samples of the first curve (an overlap there) or too long to be a
/// reach, or roots spread wider than it.
pub(super) fn tangent_line_conic_3d(
    a: &Curve,
    b: &Curve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Option<CurveIntersection<Point>> {
    let (basis_a, window_a) = through_trim(a);
    let (basis_b, window_b) = through_trim(b);
    let (line, conic, line_first) = match (basis_a, basis_b) {
        (Curve::Line(l), c) => (l, c, true),
        (c, Curve::Line(l)) => (l, c, false),
        _ => return None,
    };
    let (line_side, conic_side) = if line_first { (a, b) } else { (b, a) };
    let (centre, [x_axis, y_axis, normal], [ra, rb]) = conic_axes(conic)?;
    let reach_gap = options.gap.max(tol.confusion());
    if reach_gap >= ra.min(rb) {
        return None;
    }
    let (origin, along) = (line.axis().location, line.axis().direction.vector());
    let offset = origin - centre;
    if along.dot(normal).abs() > tol.angular() || offset.dot(normal).abs() > reach_gap {
        return None;
    }
    let scaled = ScaledConic {
        centre,
        axes: [x_axis, y_axis],
        radii: [ra, rb],
        origin,
        along,
    };
    let [dx, dy] = scaled.slope();
    let qa = dx.mul_add(dx, dy * dy);
    let t = -(offset.dot(x_axis) / ra).mul_add(dx, offset.dot(y_axis) / rb * dy) / qa;
    let vertex = scaled.at(t);
    let slope_error = scaled.slope_error();
    let qa_low =
        qa - slope_error * (2.0 * (dx.abs() + dy.abs()) + slope_error) - 3.0 * f64::EPSILON * qa;
    if !(qa_low > 0.0 && vertex.f.is_finite() && vertex.g.is_finite()) {
        return None;
    }
    let slope = vertex.g.abs() + vertex.g_error;
    let below = slope * slope / (4.0 * qa_low);
    if vertex.f + vertex.f_error < 0.0 || vertex.f - vertex.f_error - below > 0.0 {
        return None;
    }
    let spread =
        slope / (2.0 * qa_low) + ((vertex.f_error + below - vertex.f).max(0.0) / qa_low).sqrt();

    let p = origin + along * t;
    let (c_lo, c_hi) = conic.domain();
    let s = conic_parameter(conic, p, tol)?;
    let s = c_lo + (s - c_lo).rem_euclid(core::f64::consts::TAU);
    let on_conic = conic.point_at(s, tol).ok()?;
    let gap = on_conic.distance(p);
    let (lo, hi) = line.domain();
    let on_basis =
        t >= lo - tol.parametric() && t <= hi + tol.parametric() && s <= c_hi + tol.parametric();
    let touch = if line_first {
        Crossing {
            on_a: t,
            on_b: s,
            point: p,
            gap,
            reach: 0.0,
        }
    } else {
        Crossing {
            on_a: s,
            on_b: t,
            point: on_conic,
            gap,
            reach: 0.0,
        }
    };
    let kept = on_basis
        .then(|| {
            clipped_to_windows(
                CurveIntersection {
                    crossings: vec![touch],
                    overlaps: Vec::new(),
                },
                Side::of(basis_a, window_a),
                Side::of(basis_b, window_b),
                tol,
            )
            .crossings
            .first()
            .copied()
        })
        .flatten();
    let Some(mut touch) = kept else {
        let line_apart = {
            let (lo, hi) = line_side.domain();
            if t + spread < lo || t - spread > hi {
                let end = scaled.at(if t < lo { lo } else { hi });
                let radius = (1.0 + end.f - end.f_error).max(0.0).sqrt();
                ra.min(rb) * (radius - 1.0)
            } else {
                0.0
            }
        };
        let conic_apart = {
            let shadow = along - normal * along.dot(normal);
            let across = normal.cross(shadow) / shadow.magnitude();
            let (sa, sb, sc) = (
                ra * across.dot(x_axis),
                rb * across.dot(y_axis),
                across.dot(centre - origin),
            );
            let at = |u: f64| sa.mul_add(u.cos(), sb.mul_add(u.sin(), sc));
            let (lo, hi) = conic_side.domain();
            // The sinusoid is stationary at its phase and every half turn
            // on.
            let pi = core::f64::consts::PI;
            let mut values = vec![at(lo), at(hi)];
            let mut u = lo + (sb.atan2(sa) - lo).rem_euclid(pi);
            while u < hi {
                values.push(at(u));
                u += pi;
            }
            let least = values.iter().copied().fold(f64::INFINITY, f64::min);
            let most = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let rounding = 16.0 * f64::EPSILON * (sa.abs() + sb.abs() + sc.abs());
            if least > 0.0 {
                least - rounding
            } else if most < 0.0 {
                -most - rounding
            } else {
                0.0
            }
        };
        return (line_apart.max(conic_apart) > options.gap).then(CurveIntersection::empty);
    };
    if gap > options.gap {
        return None;
    }

    let distance = |q: Point| -> f64 {
        if line_first {
            distance_to_conic_3d(b, conic, touch.on_b, q, tol)
        } else {
            let (lo, hi) = b.domain();
            let foot = (q - origin).dot(along).clamp(lo, hi);
            q.distance(origin + along * foot)
        }
    };
    let ends = valley_ends_3d(a, touch.on_a, reach_gap, tol, &distance);
    let valley = valley_span_3d(a, &touch, &ends, tol);
    let extent =
        sampled_length_3d(a, options.samples, tol).min(sampled_length_3d(b, options.samples, tol));
    let cap = 4.0 * (reach_gap * extent).sqrt();
    if valley > cap || spread > valley {
        return None;
    }
    let owned = valley > reach_gap * 8.0;
    if !owned && 2.0 * valley > one_crossing_radius(tol) {
        return None;
    }
    // Two samples of the first curve inside the valley read as a stretch of
    // shared support there.
    let (a_lo, a_hi) = a.domain();
    let n = options.samples;
    let (from, to) = order(
        ends[0].1.unwrap_or(ends[0].0),
        ends[1].1.unwrap_or(ends[1].0),
    );
    #[allow(clippy::cast_precision_loss)]
    let held = (0..=n)
        .map(|i| a_lo + (a_hi - a_lo) * i as f64 / n as f64)
        .filter(|&u| u >= from && u <= to)
        .count();
    if held >= 2 {
        return None;
    }
    if owned {
        touch.reach = valley;
    }
    Some(CurveIntersection {
        crossings: vec![touch],
        overlaps: Vec::new(),
    })
}

/// A line read in a conic's frame scaled by its radii, where the conic is
/// the unit circle.
struct ScaledConic {
    centre: Point,
    axes: [ogeom_math::Vector; 2],
    radii: [f64; 2],
    origin: Point,
    along: ogeom_math::Vector,
}

/// `F = x^2 + y^2 - 1` at a point of the line, its slope along the line,
/// and how far rounding can have moved each.
struct ScaledValue {
    f: f64,
    f_error: f64,
    g: f64,
    g_error: f64,
}

impl ScaledConic {
    /// The line's direction in scaled coordinates.
    fn slope(&self) -> [f64; 2] {
        [
            self.along.dot(self.axes[0]) / self.radii[0],
            self.along.dot(self.axes[1]) / self.radii[1],
        ]
    }

    /// How far rounding can move each scaled slope component: a three-term
    /// dot product, the division, and axes unit to within their
    /// normalisation, over the smaller radius.
    fn slope_error(&self) -> f64 {
        8.0 * f64::EPSILON * self.along.magnitude() / self.radii[0].min(self.radii[1])
    }

    /// `F` and its slope at the line's point `t`. Each scaled coordinate is
    /// off by at most `e`: the point's two roundings and the subtraction of
    /// the centre on the magnitudes involved, then the dot product, the
    /// division and the axes' own normalisation on the offset, all over the
    /// smaller radius. `F` inherits `2 e` per unit of each coordinate and
    /// three roundings of its own, its slope likewise from the coordinates
    /// and the scaled direction.
    fn at(&self, t: f64) -> ScaledValue {
        let gamma = |n: f64| n * f64::EPSILON;
        let magnitude = |q: Point| (q - Point::ORIGIN).magnitude();
        let smaller = self.radii[0].min(self.radii[1]);
        let w = (self.origin + self.along * t) - self.centre;
        let (x, y) = (
            w.dot(self.axes[0]) / self.radii[0],
            w.dot(self.axes[1]) / self.radii[1],
        );
        let [dx, dy] = self.slope();
        let e = gamma(3.0)
            * (magnitude(self.origin) + self.along.magnitude() * t.abs() + magnitude(self.centre))
            / smaller
            + gamma(8.0) * w.magnitude() / smaller;
        let e_d = self.slope_error();
        let (sum_xy, sum_d) = (x.abs() + y.abs(), dx.abs() + dy.abs());
        ScaledValue {
            f: x.mul_add(x, y * y) - 1.0,
            f_error: 2.0 * e * (sum_xy + e) + gamma(3.0) * x.mul_add(x, y.mul_add(y, 1.0)),
            g: 2.0 * x.mul_add(dx, y * dy),
            g_error: 2.0 * (e * sum_d + e_d * sum_xy + 2.0 * e * e_d)
                + gamma(8.0) * ((x * dx).abs() + (y * dy).abs()),
        }
    }
}

/// The distance from `q` to a trimmed circle or ellipse: for a circle, to
/// its radial foot if the trim holds it and to the nearer end if not; for an
/// ellipse, the nearest of its foot polished from `seed` and its two ends.
fn distance_to_conic_3d(curve: &Curve, basis: &Curve, seed: f64, q: Point, tol: Tolerances) -> f64 {
    let (lo, hi) = curve.domain();
    let at = |u: f64| {
        curve
            .point_at(u, tol)
            .map_or(f64::INFINITY, |c| c.distance(q))
    };
    let ends = at(lo).min(at(hi));
    match basis {
        Curve::Circle(c) => {
            let Ok(u) = ogeom_math::elementary::circle_parameter(&c.circle(), q, tol) else {
                return ends;
            };
            let tau = core::f64::consts::TAU;
            let u = lo + (u - lo).rem_euclid(tau);
            if u <= hi { at(u) } else { ends }
        }
        _ => foot_on_3d(curve, q, seed, tol).map_or(ends, |(_, d)| d.min(ends)),
    }
}

/// The length of a curve's sampled polyline, as the sampling path measures
/// it: exact for a line, the chord sum in closed form for a circle.
fn sampled_length_3d(curve: &Curve, n: usize, tol: Tolerances) -> f64 {
    let (lo, hi) = curve.domain();
    #[allow(clippy::cast_precision_loss)]
    let segments = n as f64;
    match through_trim(curve).0 {
        Curve::Line(l) => l.axis().direction.vector().magnitude() * (hi - lo),
        Curve::Circle(c) => {
            segments * 2.0 * c.circle().radius() * ((hi - lo) / (2.0 * segments)).sin()
        }
        _ => sample_3d(curve, n, tol)
            .points
            .windows(2)
            .map(|w| w[0].distance(w[1]))
            .sum(),
    }
}

/// Two circles in one plane, in closed form: where they cross, from the
/// radical line in the plane. Exact however shallow the crossing: two rims a
/// micron apart meet at a fraction of a milliradian, and sampled, the touch
/// would read as a run millimetres long, which is how far the two stay
/// within tolerance of each other, not where they cross. `None` for circles
/// in different planes, sharing a centre, or within the weld distance of
/// touching, which the other paths answer.
fn coplanar_circles_3d(
    a: &ogeom_geom::CircleCurve,
    b: &ogeom_geom::CircleCurve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> Option<CurveIntersection<Point>> {
    let (ca, cb) = (a.circle(), b.circle());
    let normal = ca.frame().z().vector();
    if normal.cross(cb.frame().z().vector()).magnitude() > tol.angular() {
        return None;
    }
    let between = cb.centre() - ca.centre();
    if between.dot(normal).abs() > tol.confusion() {
        return None;
    }
    let in_plane = between - normal * between.dot(normal);
    let distance = in_plane.magnitude();
    let (ra, rb) = (ca.radius(), cb.radius());
    if distance <= tol.confusion() {
        return None;
    }
    if distance > ra + rb + options.gap || distance < (ra - rb).abs() - options.gap {
        return Some(CurveIntersection::empty());
    }
    // Well within the weld distance of touching (half of it, as for two
    // drums), the two crossings are one touch the root of a rounding error
    // has pulled apart.
    let weld = tol.confusion() * 50.0;
    if (distance - (ra + rb)).abs() <= weld || (distance - (ra - rb).abs()).abs() <= weld {
        return None;
    }
    let along = distance.mul_add(distance, ra.mul_add(ra, -(rb * rb))) / (2.0 * distance);
    let squared = ra.mul_add(ra, -(along * along));
    if squared <= tol.confusion() * tol.confusion() {
        return None;
    }
    let half = squared.sqrt();
    let ux = in_plane / distance;
    let uy = normal.cross(ux);
    let parameter = |curve: &ogeom_geom::CircleCurve, p: Point| -> f64 {
        let local = curve.circle().frame().to_local(p);
        let angle = local.y.atan2(local.x);
        let angle = if curve.is_reversed() { -angle } else { angle };
        let (lo, _) = Curve3d::domain(curve);
        lo + (angle - lo).rem_euclid(core::f64::consts::TAU)
    };
    let mut crossings: Vec<Crossing<Point>> = [half, -half]
        .into_iter()
        .map(|h| {
            let point = ca.centre() + ux * along + uy * h;
            Crossing {
                on_a: parameter(a, point),
                on_b: parameter(b, point),
                point,
                gap: 0.0,
                reach: 0.0,
            }
        })
        .collect();
    sort_crossings(&mut crossings);
    Some(CurveIntersection {
        crossings,
        overlaps: Vec::new(),
    })
}

/// Two circles tracing the same point set in space: the circle counterpart of
/// collinear lines, and the one 3D circle pair the sampling path cannot
/// answer: every sample is a hit, and "the crossings" do not exist. Distinct
/// circles return `None` and fall through to the general machinery, which
/// handles genuinely crossing pairs.
fn same_circle_3d(
    a: &ogeom_geom::CircleCurve,
    b: &ogeom_geom::CircleCurve,
    tol: Tolerances,
) -> Option<CurveIntersection<Point>> {
    let (ca, cb) = (a.circle(), b.circle());
    if ca.centre().distance(cb.centre()) > tol.confusion() {
        return None;
    }
    if (ca.radius() - cb.radius()).abs() > tol.confusion() {
        return None;
    }
    // Parallel or antiparallel axes both trace the same set, at whatever
    // phase and winding each was written with.
    let (za, zb) = (ca.frame().z().vector(), cb.frame().z().vector());
    if za.cross(zb).magnitude() > tol.angular() {
        return None;
    }
    // The ranges are a *correspondence*, which is what an overlap means and
    // what a caller carrying a split across the pair relies on: `on_b`'s ends
    // are the parameters at which `b` stands where `a`'s own ends do. Phase
    // comes from where `a` starts on `b`, winding from whether the two run
    // the same way there, and a pair written with opposite windings runs
    // `on_b` backwards, which is exactly the truth about them.
    let (lo, hi) = Curve3d::domain(a);
    let start = a.point_at(lo, tol).ok()?;
    let local = cb.frame().to_local(start);
    let angle = local.y.atan2(local.x);
    let phase = if b.is_reversed() { -angle } else { angle }.rem_euclid(core::f64::consts::TAU);
    let along_a = a.d1_at(lo, tol).ok()?;
    let along_b = b.d1_at(phase, tol).ok()?;
    let winding: f64 = if along_a.dot(along_b) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    Some(CurveIntersection {
        crossings: Vec::new(),
        overlaps: vec![Overlap {
            on_a: (lo, hi),
            on_b: (phase, winding.mul_add(hi - lo, phase)),
        }],
    })
}

/// The *same description* twice: one curve object meeting itself, forward
/// or reversed. A fitted seam reused as a wedge's apex ring is exactly this
/// pair, and the sampling path (every sample a hit) cannot answer it, for
/// the same reason it cannot answer coincident circles. Equality here is
/// structural, so two independent fits of one path still fall through to
/// the general machinery, which is the honest place for them.
pub(super) fn same_curve_3d(a: &Curve, b: &Curve) -> Option<CurveIntersection<Point>> {
    let (lo, hi) = Curve3d::domain(a);
    if a == b {
        return Some(CurveIntersection {
            crossings: Vec::new(),
            overlaps: vec![Overlap {
                on_a: (lo, hi),
                on_b: (lo, hi),
            }],
        });
    }
    use ogeom_geom::Reversible as _;
    if *a == b.clone().reversed() {
        let (blo, bhi) = Curve3d::domain(b);
        return Some(CurveIntersection {
            crossings: Vec::new(),
            overlaps: vec![Overlap {
                on_a: (lo, hi),
                on_b: (bhi, blo),
            }],
        });
    }
    None
}

/// Two ellipses tracing the same point set in space: the ellipse counterpart
/// of [`same_circle_3d`], and just as invisible to the sampling path. Unlike
/// a circle, an ellipse's natural parameter is pinned to its major axis, so
/// the correspondence is affine only when the two `x` axes line up (parallel
/// or antiparallel) as well as the planes and radii; anything else falls
/// through to the general machinery.
fn same_ellipse_3d(
    a: &ogeom_geom::EllipseCurve,
    b: &ogeom_geom::EllipseCurve,
    tol: Tolerances,
) -> Option<CurveIntersection<Point>> {
    let (ea, eb) = (a.ellipse(), b.ellipse());
    if ea.frame().origin().distance(eb.frame().origin()) > tol.confusion() {
        return None;
    }
    if (ea.major_radius() - eb.major_radius()).abs() > tol.confusion()
        || (ea.minor_radius() - eb.minor_radius()).abs() > tol.confusion()
    {
        return None;
    }
    let (za, zb) = (ea.frame().z().vector(), eb.frame().z().vector());
    if za.cross(zb).magnitude() > tol.angular() {
        return None;
    }
    let (xa, xb) = (ea.frame().x().vector(), eb.frame().x().vector());
    if xa.cross(xb).magnitude() > tol.angular() {
        return None;
    }
    // As for circles: phase from where `a` starts on `b`, winding from
    // whether the two run the same way there, and the ranges come back as
    // the correspondence an overlap means.
    let (lo, hi) = Curve3d::domain(a);
    let start = a.point_at(lo, tol).ok()?;
    let angle = ogeom_math::elementary::ellipse_parameter(&eb, start, tol).ok()?;
    let phase = if b.is_reversed() { -angle } else { angle }.rem_euclid(core::f64::consts::TAU);
    let along_a = a.d1_at(lo, tol).ok()?;
    let along_b = b.d1_at(phase, tol).ok()?;
    let winding: f64 = if along_a.dot(along_b) >= 0.0 {
        1.0
    } else {
        -1.0
    };
    Some(CurveIntersection {
        crossings: Vec::new(),
        overlaps: vec![Overlap {
            on_a: (lo, hi),
            on_b: (phase, winding.mul_add(hi - lo, phase)),
        }],
    })
}
