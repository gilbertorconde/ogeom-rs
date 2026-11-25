//! Closed forms for lines and circles in the plane, and for two lines in space.

use super::*;

// --- analytic, planar --------------------------------------------------------

pub(super) fn line_line_2d(
    a: &ogeom_geom::Line2d,
    b: &ogeom_geom::Line2d,
    tol: Tolerances,
) -> CurveIntersection<Point2> {
    let (oa, da) = (a.axis().location, a.axis().direction.vector());
    let (ob, db) = (b.axis().location, b.axis().direction.vector());
    let cross = da.cross(db);

    if cross.abs() <= tol.angular() {
        // Parallel. Collinear if one origin is on the other line.
        let between = ob - oa;
        if between.cross(da).abs() > tol.confusion() {
            return CurveIntersection::empty();
        }
        // The shared stretch, as each line's own parameter range.
        let (a_lo, a_hi) = a.domain();
        let (b_lo, b_hi) = b.domain();
        // Where b's range lands on a's parameter: t_a = (p - oa)·da.
        let project = |p: Point2| (p - oa).dot(da);
        let (s0, s1) = (project(ob + db * b_lo), project(ob + db * b_hi));
        let (lo, hi) = (s0.min(s1).max(a_lo), s0.max(s1).min(a_hi));
        // And back onto b.
        let back = |t: f64| (oa + da * t - ob).dot(db);
        if hi - lo <= tol.confusion() {
            // Segments meeting end to end share a point, not a stretch.
            if lo - hi > tol.confusion() {
                return CurveIntersection::empty();
            }
            let t = f64::midpoint(lo, hi).clamp(a_lo, a_hi);
            return CurveIntersection {
                crossings: vec![Crossing {
                    on_a: t,
                    on_b: back(t).clamp(b_lo, b_hi),
                    point: oa + da * t,
                    gap: 0.0,
                    reach: 0.0,
                }],
                overlaps: Vec::new(),
            };
        }
        return CurveIntersection {
            crossings: Vec::new(),
            overlaps: vec![Overlap {
                on_a: (lo, hi),
                // Paired end to end with `on_a`, not sorted: two lines written in
                // opposite directions run `on_b` backwards, and a consumer
                // carrying a stretch across by the correspondence (the
                // boolean clipping a contact to the edge it runs along)
                // reads a sorted pair as the reflected stretch.
                on_b: (back(lo), back(hi)),
            }],
        };
    }

    let between = ob - oa;
    let t = between.cross(db) / cross;
    let s = between.cross(da) / cross;
    let (a_lo, a_hi) = a.domain();
    let (b_lo, b_hi) = b.domain();
    if t < a_lo - tol.parametric()
        || t > a_hi + tol.parametric()
        || s < b_lo - tol.parametric()
        || s > b_hi + tol.parametric()
    {
        return CurveIntersection::empty();
    }
    CurveIntersection {
        crossings: vec![Crossing {
            on_a: t,
            on_b: s,
            point: oa + da * t,
            gap: 0.0,
            reach: 0.0,
        }],
        overlaps: Vec::new(),
    }
}

pub(super) fn line_circle_2d(
    line: &ogeom_geom::Line2d,
    circle: &ogeom_geom::Circle2d,
    swapped: bool,
    tol: Tolerances,
) -> CurveIntersection<Point2> {
    let (o, d) = (line.axis().location, line.axis().direction.vector());
    let c = circle.circle();
    let centre = c.centre();
    let radius = c.radius();

    // Foot of the perpendicular from the centre onto the line.
    let along = (centre - o).dot(d);
    let foot = o + d * along;
    let gap = foot.distance(centre);
    if gap > radius + tol.confusion() {
        return CurveIntersection::empty();
    }
    let half = radius.mul_add(radius, -(gap * gap)).max(0.0).sqrt();
    let candidates = if half <= tol.confusion() {
        vec![along]
    } else {
        vec![along - half, along + half]
    };

    let (l_lo, l_hi) = line.domain();
    let mut crossings = Vec::new();
    for t in candidates {
        if t < l_lo - tol.parametric() || t > l_hi + tol.parametric() {
            continue;
        }
        let p = o + d * t;
        let Some(s) = circle_parameter(circle, p, tol) else {
            continue;
        };
        let (on_a, on_b) = if swapped { (s, t) } else { (t, s) };
        crossings.push(Crossing {
            on_a,
            on_b,
            point: p,
            gap: 0.0,
            reach: 0.0,
        });
    }
    sort_crossings(&mut crossings);
    CurveIntersection {
        crossings,
        overlaps: Vec::new(),
    }
}

pub(super) fn circle_circle_2d(
    a: &ogeom_geom::Circle2d,
    b: &ogeom_geom::Circle2d,
    tol: Tolerances,
) -> CurveIntersection<Point2> {
    let (ca, cb) = (a.circle(), b.circle());
    let between = cb.centre() - ca.centre();
    let distance = between.magnitude();
    let (ra, rb) = (ca.radius(), cb.radius());

    if distance <= tol.confusion() {
        if (ra - rb).abs() <= tol.confusion() {
            // The same circle, the whole turn of each. As in space: phase
            // from where `a` starts on `b`, winding from whether the two run
            // the same way there, so `on_b` names where `b` stands at
            // `a`'s own ends.
            use ogeom_geom::Curve2d as _;
            let (lo, hi) = a.domain();
            let correspondence = (|| {
                let start = a.point_at(lo, tol).ok()?;
                let phase = circle_parameter(b, start, tol)?;
                let along_a = a.d1_at(lo, tol).ok()?;
                let along_b = b.d1_at(phase, tol).ok()?;
                let winding: f64 = if along_a.dot(along_b) >= 0.0 {
                    1.0
                } else {
                    -1.0
                };
                Some((phase, winding.mul_add(hi - lo, phase)))
            })();
            return CurveIntersection {
                crossings: Vec::new(),
                overlaps: vec![Overlap {
                    on_a: (lo, hi),
                    on_b: correspondence.unwrap_or_else(|| b.domain()),
                }],
            };
        }
        return CurveIntersection::empty();
    }
    if distance > ra + rb + tol.confusion() || distance < (ra - rb).abs() - tol.confusion() {
        return CurveIntersection::empty();
    }

    // The radical line: where the two circles' equations agree.
    let along = distance.mul_add(distance, ra.mul_add(ra, -(rb * rb))) / (2.0 * distance);
    let squared = ra.mul_add(ra, -(along * along));
    let direction = between * (1.0 / distance);
    let foot = ca.centre() + direction * along;
    let mut crossings = Vec::new();
    let mut push = |p: Point2| {
        if let (Some(s), Some(t)) = (circle_parameter(a, p, tol), circle_parameter(b, p, tol)) {
            crossings.push(Crossing {
                on_a: s,
                on_b: t,
                point: p,
                gap: 0.0,
                reach: 0.0,
            });
        }
    };
    if squared <= tol.confusion() * tol.confusion() {
        push(foot);
    } else {
        let offset = ogeom_math::Vector2::new(-direction.y, direction.x) * squared.max(0.0).sqrt();
        push(foot + offset);
        push(foot - offset);
    }
    sort_crossings(&mut crossings);
    CurveIntersection {
        crossings,
        overlaps: Vec::new(),
    }
}

/// The parameter at which a circle passes through a point on it.
fn circle_parameter(curve: &ogeom_geom::Circle2d, p: Point2, tol: Tolerances) -> Option<f64> {
    let c = curve.circle();
    let local = p - c.centre();
    let x = local.dot(c.frame().x().vector());
    let y = local.dot(c.frame().y().vector());
    let mut angle = y.atan2(x);
    if curve.is_reversed() {
        angle = -angle;
    }
    let angle = angle.rem_euclid(core::f64::consts::TAU);
    let (lo, hi) = curve.domain();
    // Fold into the arc's own range where the arc covers it.
    if angle >= lo - tol.parametric() && angle <= hi + tol.parametric() {
        return Some(angle.clamp(lo, hi));
    }
    let shifted = angle - core::f64::consts::TAU;
    if shifted >= lo - tol.parametric() && shifted <= hi + tol.parametric() {
        return Some(shifted.clamp(lo, hi));
    }
    None
}

// --- analytic, spatial -------------------------------------------------------

pub(super) fn line_line_3d(
    a: &ogeom_geom::LineCurve,
    b: &ogeom_geom::LineCurve,
    options: CurveCurveOptions,
    tol: Tolerances,
) -> CurveIntersection<Point> {
    let (oa, da) = (a.axis().location, a.axis().direction.vector());
    let (ob, db) = (b.axis().location, b.axis().direction.vector());
    let cross = da.cross(db);
    let denominator = cross.square_magnitude();

    if denominator <= tol.angular() * tol.angular() {
        // Parallel: collinear overlap or nothing.
        let between = ob - oa;
        if between.cross(da).magnitude() > tol.confusion() {
            return CurveIntersection::empty();
        }
        let (a_lo, a_hi) = a.domain();
        let (b_lo, b_hi) = b.domain();
        let project = |p: Point| (p - oa).dot(da);
        let (s0, s1) = (project(ob + db * b_lo), project(ob + db * b_hi));
        let (lo, hi) = (s0.min(s1).max(a_lo), s0.max(s1).min(a_hi));
        let back = |t: f64| (oa + da * t - ob).dot(db);
        if hi - lo <= tol.confusion() {
            // Segments meeting end to end share a point, not a stretch.
            if lo - hi > tol.confusion() {
                return CurveIntersection::empty();
            }
            let t = f64::midpoint(lo, hi).clamp(a_lo, a_hi);
            let s = back(t).clamp(b_lo, b_hi);
            let point = oa + da * t;
            return CurveIntersection {
                crossings: vec![Crossing {
                    on_a: t,
                    on_b: s,
                    point,
                    gap: point.distance(ob + db * s),
                    reach: 0.0,
                }],
                overlaps: Vec::new(),
            };
        }
        return CurveIntersection {
            crossings: Vec::new(),
            overlaps: vec![Overlap {
                on_a: (lo, hi),
                // Paired end to end with `on_a`, not sorted: two lines written in
                // opposite directions run `on_b` backwards, and a consumer
                // carrying a stretch across by the correspondence (the
                // boolean clipping a contact to the edge it runs along)
                // reads a sorted pair as the reflected stretch.
                on_b: (back(lo), back(hi)),
            }],
        };
    }

    // Closest approach of two skew lines, in closed form.
    let between = ob - oa;
    let t = between.cross(db).dot(cross) / denominator;
    let s = between.cross(da).dot(cross) / denominator;
    let pa = oa + da * t;
    let pb = ob + db * s;
    let gap = pa.distance(pb);
    let (a_lo, a_hi) = a.domain();
    let (b_lo, b_hi) = b.domain();
    if gap > options.gap
        || t < a_lo - tol.parametric()
        || t > a_hi + tol.parametric()
        || s < b_lo - tol.parametric()
        || s > b_hi + tol.parametric()
    {
        return CurveIntersection::empty();
    }
    CurveIntersection {
        crossings: vec![Crossing {
            on_a: t,
            on_b: s,
            point: pa,
            gap,
            reach: 0.0,
        }],
        overlaps: Vec::new(),
    }
}
