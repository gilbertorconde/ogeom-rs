//! Pcurves derived beside an edge, and the maps between a file's parameters and the kernel's.

use super::*;

/// Each wire walked the other way round: its edges in reverse order, each
/// reversed.
pub(super) fn walked_back(
    model: &mut Model,
    wires: &[Shape],
    tol: Tolerances,
) -> OgeomResult<Vec<Shape>> {
    wires
        .iter()
        .map(|wire| {
            let back = model.ordered_children_of(&wire.reversed())?;
            Ok(make_wire(model, &back, tol)?.shape)
        })
        .collect()
}

/// An edge's pcurve on a surface: the file's own where it lies within the
/// edge's tolerance of the curve once lifted, else derived, exact where
/// the pair has a closed form and fitted where it does not, with how far it
/// stands from the curve once lifted through the surface.
pub(super) fn derive_pcurve(
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    stated: &[PlanarCurve],
    tol: Tolerances,
) -> PreparedPcurve {
    // A pcurve that does not evaluate over the edge's range states no gap
    // here; the face that uses it fails on it there. With nothing stated,
    // the gap is only asked whether it passes the confusion distance, and
    // by how much where it does; with the file's own to weigh it against,
    // it is the bound they are held to.
    let lifted = |pcurve: &PlanarCurve| {
        let gap = if stated.is_empty() {
            ogeom_algo::pcurve_fit::lifted_gap_past_confusion
        } else {
            ogeom_algo::pcurve_fit::lifted_gap
        };
        gap((curve, range), (pcurve, range), surface, false, tol).unwrap_or(0.0)
    };
    let nearest = nearest_stated(curve, range, surface, stated, tol);
    // One that follows the curve exactly needs nothing derived beside it.
    if let Some((curve, off)) = &nearest
        && *off <= tol.confusion()
    {
        return PreparedPcurve::Stated {
            curve: curve.clone(),
            off: *off,
        };
    }
    let derived =
        if let Some(exact) = ogeom_intersect::exact_pcurve_over(curve, range, surface, tol) {
            let exact = crate::pcurves::over_range(&exact, range);
            let off = lifted(&exact);
            PreparedPcurve::Exact { curve: exact, off }
        } else {
            // No closed form: a spline surface, or a combination the
            // projection table lacks. The pcurve is fitted at the curve's own
            // parameters, so same-parameter is preserved by construction.
            match crate::pcurves::fit_projected_pcurve(curve, range, surface, tol) {
                Ok((fitted, error, met, worst_off, warning)) => {
                    let off = lifted(&fitted);
                    PreparedPcurve::Fitted {
                        curve: fitted,
                        error,
                        met,
                        worst_off,
                        off,
                        warning,
                    }
                }
                Err(e) => PreparedPcurve::Refused(e.to_string()),
            }
        };
    // The edge's tolerance on this face is how far its curve stands from
    // the surface, which the projection measures. A pcurve whose trace
    // keeps within that of the projected trace lifts to within twice it of
    // the curve, and within that bound the file's pcurve is the same
    // boundary drawn as its sender drew it, and is kept.
    let allowed = match &derived {
        PreparedPcurve::Exact { off, .. } => Some(*off),
        PreparedPcurve::Fitted { off, worst_off, .. } => Some(off.max(*worst_off)),
        PreparedPcurve::Stated { .. } | PreparedPcurve::Refused(_) => None,
    };
    let Some(allowed) = allowed else {
        return derived;
    };
    let bound = 2.0f64.mul_add(allowed, tol.confusion());
    if let Some((curve, off)) = nearest
        && off <= bound
    {
        return PreparedPcurve::Stated { curve, off };
    }
    // The file may state the boundary at another pace or phase than its
    // curve's: a conic its frame runs the other way round, a line or a
    // spline over a parameter of its own. Each read at the map of the
    // curve's parameter that puts it on the edge's ends is the same
    // boundary, and is held to the same bound.
    let ends = match &derived {
        PreparedPcurve::Exact { curve: image, .. }
        | PreparedPcurve::Fitted { curve: image, .. } => image
            .point_at(range.0, tol)
            .ok()
            .zip(image.point_at(range.1, tol).ok()),
        PreparedPcurve::Stated { .. } | PreparedPcurve::Refused(_) => None,
    };
    if let Some(ends) = ends {
        let repaced: Vec<PlanarCurve> = stated
            .iter()
            .flat_map(|pcurve| repaced(pcurve, range, ends, tol))
            .collect();
        if let Some((curve, off)) = nearest_stated(curve, range, surface, &repaced, tol)
            && off <= bound
        {
            return PreparedPcurve::Stated { curve, off };
        }
    }
    derived
}

/// The file's pcurve nearest the curve once lifted, with how far it
/// stands; within the confusion distance, only known to be within it. One
/// that does not evaluate over the edge's range, or whose
/// image leaves the surface's chart along a direction that does not close,
/// is left out.
///
/// A closed curve's range here may begin whole periods from where the
/// file's parameter began it, so each is also tried a period either way:
/// on a sphere a meridian's image a period off lifts to the same points,
/// at latitudes past the poles, and only the right one stays in the chart.
/// The first that follows the curve exactly ends the search.
fn nearest_stated(
    curve: &Curve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    stated: &[PlanarCurve],
    tol: Tolerances,
) -> Option<(PlanarCurve, f64)> {
    let period = curve.is_periodic().then(|| {
        let (lo, hi) = curve.domain();
        hi - lo
    });
    let turns: &[f64] = if period.is_some() {
        &[0.0, 1.0, -1.0]
    } else {
        &[0.0]
    };
    let mut best: Option<(PlanarCurve, f64)> = None;
    for pcurve in stated {
        for turn in turns {
            let Some(shifted) = period_shifted(pcurve, turn * period.unwrap_or(0.0)) else {
                continue;
            };
            let shifted = crate::pcurves::over_range(&shifted, range);
            if !within_chart(&shifted, range, surface, tol) {
                continue;
            }
            let Ok(off) = ogeom_algo::pcurve_fit::lifted_gap_past_confusion(
                (curve, range),
                (&shifted, range),
                surface,
                false,
                tol,
            ) else {
                continue;
            };
            if best.as_ref().is_none_or(|(_, held)| off < *held) {
                best = Some((shifted, off));
            }
            if off <= tol.confusion() {
                return best;
            }
        }
    }
    best
}

/// A pcurve read `delta` later along its parameter: at `t`, the point the
/// pcurve has at `t - delta`. `None` for a kind that cannot be restated so.
fn period_shifted(pcurve: &PlanarCurve, delta: f64) -> Option<PlanarCurve> {
    if delta == 0.0 {
        return Some(pcurve.clone());
    }
    match pcurve {
        PlanarCurve::Line(l) => {
            let axis = l.axis();
            let moved = ogeom_math::Axis2 {
                location: axis.location - axis.direction.vector() * delta,
                direction: axis.direction,
            };
            let (lo, hi) = l.domain();
            ogeom_geom::Line2d::over(moved, lo + delta, hi + delta)
                .ok()
                .map(Into::into)
        }
        PlanarCurve::BSpline(b) => {
            let knots: Vec<f64> = b.knots().knots().iter().map(|k| k + delta).collect();
            let knots = KnotVector::new(knots, b.knots().degree()).ok()?;
            ogeom_geom::BSpline2d::rational(knots, b.control_points().to_vec())
                .ok()
                .map(Into::into)
        }
        // A conic's parameter is its angle: whole turns leave it in place.
        PlanarCurve::Circle(_) | PlanarCurve::Ellipse(_)
            if (delta / core::f64::consts::TAU).fract().abs() < 1e-12 =>
        {
            Some(pcurve.clone())
        }
        _ => None,
    }
}

/// A file's pcurve in this kernel's parameters: its chart's coordinates
/// stretched from the file's units, `u` by `su` and `v` by `sv`, and its
/// parameter carried onto the edge's curve's here, which runs `pace` per
/// file unit. A conic as read already runs in radians, `angle` per file
/// unit; a line and a spline run in the file's own. `None` where no curve
/// here states the image at the curve's parameters: a conic paced unlike
/// its curve or stretched unevenly.
pub(super) fn in_kernel_units(
    pcurve: &PlanarCurve,
    (su, sv): (f64, f64),
    pace: f64,
    angle: f64,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    let own = match pcurve {
        PlanarCurve::Circle(_) | PlanarCurve::Ellipse(_) => angle,
        _ => 1.0,
    };
    if su == 1.0 && sv == 1.0 && own == pace {
        return Some(pcurve.clone());
    }
    match pcurve {
        PlanarCurve::Line(l) => {
            // Stretched, the line runs `speed` along its new direction per
            // unit of its parameter; at the curve's `t` it stands at
            // `t·own/pace` along that parameter. The speed is the stretch
            // itself along an axis, so a line in the file's own pace comes
            // out at exactly this kernel's.
            let axis = l.axis();
            let d = axis.direction.vector();
            let along = ogeom_math::Vector2::new(d.x * su, d.y * sv);
            let speed = if su == sv || d.y == 0.0 {
                su.abs()
            } else if d.x == 0.0 {
                sv.abs()
            } else {
                along.magnitude()
            };
            let origin = ogeom_math::Point2::new(axis.location.x * su, axis.location.y * sv);
            let line =
                ogeom_math::Axis2::new(origin, ogeom_math::Direction2::new(along, tol).ok()?);
            let (lo, hi) = l.domain();
            let stretched: PlanarCurve = ogeom_geom::Line2d::over(line, lo * speed, hi * speed)
                .ok()?
                .into();
            paced(&stretched, speed * own / pace, 0.0, tol)
        }
        PlanarCurve::Circle(_) | PlanarCurve::Ellipse(_) => {
            if own != pace {
                return None;
            }
            pcurve.stretched(su, sv, tol).ok()
        }
        PlanarCurve::BSpline(_) => {
            let stretched = pcurve.stretched(su, sv, tol).ok()?;
            paced(&stretched, own / pace, 0.0, tol)
        }
        _ => None,
    }
}

/// An affine map of a chart's coordinates: `linear` times the point, plus
/// `offset`.
#[derive(Debug, Clone, Copy)]
pub(super) struct ChartMap {
    pub(super) linear: [[f64; 2]; 2],
    pub(super) offset: (f64, f64),
}

impl ChartMap {
    fn point(&self, p: ogeom_math::Point2) -> ogeom_math::Point2 {
        let [[a, b], [c, d]] = self.linear;
        ogeom_math::Point2::new(
            a.mul_add(p.x, b * p.y) + self.offset.0,
            c.mul_add(p.x, d * p.y) + self.offset.1,
        )
    }

    fn vector(&self, v: ogeom_math::Vector2) -> ogeom_math::Vector2 {
        let [[a, b], [c, d]] = self.linear;
        ogeom_math::Vector2::new(a.mul_add(v.x, b * v.y), c.mul_add(v.x, d * v.y))
    }
}

/// A file's pcurve carried through `map` into this kernel's chart, its
/// parameter carried onto the edge's curve's here, which runs `pace` per
/// file unit. An affine map takes a line to a line and a spline to the
/// spline on its mapped control points, both exactly, at their own
/// parameters; `None` for a conic, whose image under a map that is not a
/// similarity is an ellipse at another phase than any stated here, and
/// for anything else.
pub(super) fn through_map(
    pcurve: &PlanarCurve,
    map: ChartMap,
    pace: f64,
    tol: Tolerances,
) -> Option<PlanarCurve> {
    match pcurve {
        PlanarCurve::Line(l) => {
            let axis = l.axis();
            let along = map.vector(axis.direction.vector());
            let speed = along.magnitude();
            let line = ogeom_math::Axis2::new(
                map.point(axis.location),
                ogeom_math::Direction2::new(along, tol).ok()?,
            );
            let (lo, hi) = l.domain();
            let mapped: PlanarCurve = ogeom_geom::Line2d::over(line, lo * speed, hi * speed)
                .ok()?
                .into();
            paced(&mapped, speed / pace, 0.0, tol)
        }
        PlanarCurve::BSpline(b) => {
            let control = b
                .control_points()
                .iter()
                .map(|w| ogeom_math::Weighted::new(map.point(w.point()), w.weight, tol))
                .collect::<OgeomResult<Vec<_>>>()
                .ok()?;
            let mapped: PlanarCurve = ogeom_geom::BSpline2d::rational(b.knots().clone(), control)
                .ok()?
                .into();
            paced(&mapped, 1.0 / pace, 0.0, tol)
        }
        _ => None,
    }
}

/// A pcurve read at an affine map of the curve's parameter: at `t`, the
/// point `pcurve` has at `a·t + b`. Exact for each kind it answers: a line
/// (straight at any pace, a line of its own at a pace of one either way), a
/// spline (its knots carried back through the map), a circle at a pace of
/// one either way (its frame turned by the phase), and an ellipse the same
/// at a phase of whole half turns, the only ones its frame turns through
/// with its axes as they are. `None` for anything else.
pub(super) fn paced(pcurve: &PlanarCurve, a: f64, b: f64, tol: Tolerances) -> Option<PlanarCurve> {
    use ogeom_geom::Reversible as _;
    if !a.is_finite() || !b.is_finite() || a == 0.0 {
        return None;
    }
    if a == 1.0 && b == 0.0 {
        return Some(pcurve.clone());
    }
    // Where the map sends each end of the pcurve's own domain, in order.
    let back = |lo: f64, hi: f64| {
        let (x, y) = ((lo - b) / a, (hi - b) / a);
        (x.min(y), x.max(y))
    };
    match pcurve {
        PlanarCurve::Line(l) => {
            let axis = l.axis();
            let (lo, hi) = l.domain();
            let (from, to) = back(lo, hi);
            if a.abs() == 1.0 {
                let direction = if a > 0.0 {
                    axis.direction
                } else {
                    axis.direction.reversed()
                };
                let line = ogeom_math::Axis2::new(axis.point_at(b), direction);
                return ogeom_geom::Line2d::over(line, from, to)
                    .ok()
                    .map(Into::into);
            }
            let at = |t: f64| axis.point_at(a.mul_add(t, b));
            let knots = KnotVector::new(vec![from, from, to, to], 1).ok()?;
            ogeom_geom::BSpline2d::new(knots, vec![at(from), at(to)], tol)
                .ok()
                .map(Into::into)
        }
        PlanarCurve::BSpline(spline) => {
            let mut knots: Vec<f64> = spline.knots().knots().iter().map(|k| (k - b) / a).collect();
            let mut control = spline.control_points().to_vec();
            if a < 0.0 {
                knots.reverse();
                control.reverse();
            }
            let knots = KnotVector::new(knots, spline.knots().degree()).ok()?;
            ogeom_geom::BSpline2d::rational(knots, control)
                .ok()
                .map(Into::into)
        }
        PlanarCurve::Circle(c) if a.abs() == 1.0 => {
            // At `τ` the circle stands at angle `σ·τ` in its frame, `σ`
            // the way it runs; at `a·t + b` that is `σ·a·t + σ·b`.
            let sense = if c.is_reversed() { -1.0 } else { 1.0 };
            let circle = c.circle();
            let frame = circle.frame();
            let (sin, cos) = (sense * b).sin_cos();
            let (x, y) = (frame.x().vector(), frame.y().vector());
            let turned = ogeom_math::Frame2::from_axes(
                frame.origin(),
                ogeom_math::Direction2::new(x * cos + y * sin, tol).ok()?,
                ogeom_math::Direction2::new(y * cos - x * sin, tol).ok()?,
                tol,
            )
            .ok()?;
            let out: PlanarCurve = ogeom_geom::Circle2d::new(
                ogeom_math::Circle2::new(turned, circle.radius(), tol).ok()?,
            )
            .into();
            Some(if sense * a < 0.0 { out.reversed() } else { out })
        }
        PlanarCurve::Ellipse(e) if a.abs() == 1.0 => {
            let half_turns = (b / core::f64::consts::PI).round();
            if core::f64::consts::PI * half_turns != b {
                return None;
            }
            let sense = if e.is_reversed() { -1.0 } else { 1.0 };
            let ellipse = e.ellipse();
            let frame = ellipse.frame();
            let frame = if half_turns.rem_euclid(2.0) == 0.0 {
                frame
            } else {
                ogeom_math::Frame2::from_axes(
                    frame.origin(),
                    frame.x().reversed(),
                    frame.y().reversed(),
                    tol,
                )
                .ok()?
            };
            let out: PlanarCurve = ogeom_geom::Ellipse2d::new(
                ogeom_math::Ellipse2::new(
                    frame,
                    ellipse.major_radius(),
                    ellipse.minor_radius(),
                    tol,
                )
                .ok()?,
            )
            .into();
            Some(if sense * a < 0.0 { out.reversed() } else { out })
        }
        _ => None,
    }
}

/// The ways a file's pcurve may follow the edge's curve at another pace or
/// phase than the curve's own, each read at an affine map of the curve's
/// parameter that puts its ends where `ends` stand in the chart (the
/// edge's ends there, from a pcurve derived for it): forward and reversed
/// at the curve's own pace, and a line at the pace its two ends ask for, a
/// spline over its own domain either way.
fn repaced(
    pcurve: &PlanarCurve,
    range: (f64, f64),
    ends: (ogeom_math::Point2, ogeom_math::Point2),
    tol: Tolerances,
) -> Vec<PlanarCurve> {
    // A straight spline of two points (a line the file paces other than
    // by length) is the line through them, at whatever pace it runs.
    let straight = match pcurve {
        PlanarCurve::BSpline(spline) if spline.knots().degree() == 1 && !spline.is_rational() => {
            match spline.control_points() {
                [from, to] => Some((from.point(), to.point())),
                _ => None,
            }
        }
        _ => None,
    };
    let line;
    let pcurve = match straight {
        Some((from, to)) => {
            let Ok(segment) = ogeom_geom::Line2d::segment(from, to, tol) else {
                return Vec::new();
            };
            line = PlanarCurve::from(segment);
            &line
        }
        None => pcurve,
    };
    let span = range.1 - range.0;
    // The parameter of the pcurve's point nearest `w`, for the kinds that
    // have it in closed form.
    let foot = |w: ogeom_math::Point2| -> Option<f64> {
        match pcurve {
            PlanarCurve::Line(l) => {
                let axis = l.axis();
                Some((w - axis.location).dot(axis.direction.vector()))
            }
            PlanarCurve::Circle(c) => {
                let frame = c.circle().frame();
                let r = w - frame.origin();
                let angle = r.dot(frame.y().vector()).atan2(r.dot(frame.x().vector()));
                Some(if c.is_reversed() { -angle } else { angle })
            }
            PlanarCurve::Ellipse(e) => {
                let ellipse = e.ellipse();
                let frame = ellipse.frame();
                let r = w - frame.origin();
                let angle = (r.dot(frame.y().vector()) / ellipse.minor_radius())
                    .atan2(r.dot(frame.x().vector()) / ellipse.major_radius());
                Some(if e.is_reversed() { -angle } else { angle })
            }
            _ => None,
        }
    };
    let mut maps: Vec<(f64, f64)> = Vec::new();
    match pcurve {
        PlanarCurve::BSpline(spline) => {
            let (lo, hi) = spline.domain();
            let rate = (hi - lo) / span;
            maps.push((rate, rate.mul_add(-range.0, lo)));
            maps.push((-rate, rate.mul_add(range.0, hi)));
        }
        PlanarCurve::Line(_) | PlanarCurve::Circle(_) | PlanarCurve::Ellipse(_) => {
            let Some(start) = foot(ends.0) else {
                return Vec::new();
            };
            for a in [1.0_f64, -1.0] {
                let b = a.mul_add(-range.0, start);
                // An ellipse's frame turns only by whole half turns.
                let b = if matches!(pcurve, PlanarCurve::Ellipse(_)) {
                    core::f64::consts::PI * (b / core::f64::consts::PI).round()
                } else {
                    b
                };
                maps.push((a, b));
            }
            if let PlanarCurve::Line(_) = pcurve
                && let Some(end) = foot(ends.1)
            {
                let rate = (end - start) / span;
                maps.push((rate, rate.mul_add(-range.0, start)));
            }
        }
        _ => {}
    }
    maps.into_iter()
        .filter_map(|(a, b)| paced(pcurve, a, b, tol))
        .collect()
}

/// Whether a pcurve's image over a range stays inside the surface's chart
/// along each direction that does not close on itself.
fn within_chart(
    pcurve: &PlanarCurve,
    range: (f64, f64),
    surface: &SurfaceGeometry,
    tol: Tolerances,
) -> bool {
    let ((u0, u1), (v0, v1)) = surface.domain();
    let slack = tol.parametric();
    (0..=4).all(|i| {
        let t = (range.1 - range.0).mul_add(f64::from(i) / 4.0, range.0);
        pcurve.point_at(t, tol).is_ok_and(|p| {
            (surface.is_periodic_u() || (p.x >= u0 - slack && p.x <= u1 + slack))
                && (surface.is_periodic_v() || (p.y >= v0 - slack && p.y <= v1 + slack))
        })
    })
}
