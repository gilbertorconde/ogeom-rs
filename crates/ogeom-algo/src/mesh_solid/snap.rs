//! Curves snapped to the chains of boundary vertices between faces: a curve
//! of a curved face's own surface where one holds the chain, a section
//! solved where the two surfaces meet, or a chord threaded through the
//! chain's vertices as the last resort, and each edge's image on the curved
//! faces it bounds.

use ogeom_core::{OgeomResult, Tolerances};
use ogeom_geom::{Curve, LineCurve, PlanarCurve, PlaneSurface};
use ogeom_math::{Direction, Frame, Plane, Point, Point2, Vector};

use super::builder::linear;
use super::planner::{CHORD_SAG, Planner, REACH};
use super::seams::distance_to_line;
use super::segment::{angular_spread, axis_frame, chart, evaluate, gradient, periodic, unwrapped};
use super::{Carrier, Curved};
use crate::recognize::Canonical;

impl Planner<'_> {
    /// As [`Self::snapped`], solved afresh.
    pub(super) fn snap(&self, chain: &[u32], faces: &[usize]) -> Option<(Snapped, bool, Images)> {
        if faces.len() > 2 {
            return None;
        }
        let closed = chain.len() > 2 && chain[0] == chain[chain.len() - 1];
        let pts: Vec<Point> = chain[..chain.len() - usize::from(closed)]
            .iter()
            .map(|&v| self.points[v as usize])
            .collect();
        let reach = self.flat * REACH;
        // A band round its axis bounds itself with its own rims, so its
        // seam meets their vertices; its candidates go first.
        let mut order: Vec<usize> = faces.to_vec();
        order.sort_by_key(|&g| !self.curved(g).is_some_and(|c| c.wraps || c.wraps_v));
        // A recognized plane across the chain is where the chain lies,
        // exactly: better than a plane fitted to the chain's own points.
        let across = faces.iter().find_map(|&g| match &self.groups.carriers[g] {
            Carrier::Plane(plane) => Some(*plane),
            _ => None,
        });
        for &g in &order {
            let Some(curved) = self.curved(g) else {
                continue;
            };
            for candidate in candidates(&curved.shape, &pts, across, reach, self.tol) {
                if let Some((snapped, forward)) = self.fitted(candidate, &pts, closed, faces, reach)
                {
                    return Some((snapped, forward, Vec::new()));
                }
            }
        }
        self.section(&pts, closed, faces, reach)
            .or_else(|| self.chord(&pts, closed, faces, reach))
    }

    /// The last resort between two faces that meet all but tangentially (a
    /// fillet running on into a corner ball, or into a patch the mesh
    /// leaves faceted): no curve the two surfaces share can be solved for
    /// along the chain. The chain's own vertices lie on both, and a curve
    /// is threaded through them, a line for a single span; its tolerance is
    /// how far it strays from either surface between them, up to a
    /// twentieth of its longest span. Between two curved faces, points
    /// between the vertices carried onto both surfaces are threaded as
    /// well, and the closer curve kept (the only one, where the vertices
    /// alone thread none). A face
    /// bounded by such a curve is good to its tolerance, where it would
    /// otherwise fall to facets whole.
    pub(super) fn chord(
        &self,
        pts: &[Point],
        closed: bool,
        faces: &[usize],
        reach: f64,
    ) -> Option<(Snapped, bool, Images)> {
        let [a, b] = faces[..] else {
            return None;
        };
        let (fa, fb) = (self.signed(a)?, self.signed(b)?);
        let mut on: Vec<Point> = if closed {
            pts.to_vec()
        } else {
            // Along a tangency the chain wanders, and can step past its end
            // and back; a curve through such a step hooks back on itself.
            // Only the points that close in on the end are threaded.
            let end = pts[pts.len() - 1];
            let mut kept = vec![pts[0]];
            for &p in &pts[1..pts.len() - 1] {
                let last = kept[kept.len() - 1];
                if p.distance(end) < last.distance(end) && (p - last).dot(end - last) > 0.0 {
                    kept.push(p);
                }
            }
            kept.push(end);
            kept
        };
        if closed {
            on.push(pts[0]);
        }
        let longest = on
            .windows(2)
            .map(|w| w[0].distance(w[1]))
            .fold(0.0_f64, f64::max);
        if longest <= self.tol.confusion() {
            return None;
        }
        // The chain's own vertices lie on both surfaces, to the reach; the
        // sag a chord may stand off between them covers the surfaces'
        // bulge over the chain, not a surface standing off its vertices.
        if on
            .iter()
            .any(|&p| fa(p).abs() > reach || fb(p).abs() > reach)
        {
            return None;
        }
        // Between a canonical surface and a patch running out into it
        // tangentially, the chain is threaded in the canonical surface's
        // chart and lifted onto it: a curve through the vertices alone
        // stands off a curved surface between them, most where the chain
        // turns, and its image there strays as far.
        for (g, h) in [(a, b), (b, a)] {
            if self.curved(h).is_some_and(|c| c.patch.is_some())
                && let Some((curve, range, tolerance, images)) =
                    self.chart_thread(&on, closed, g, h, reach)
            {
                let snapped = if closed {
                    Snapped::Loop(curve, range, tolerance)
                } else {
                    Snapped::Open(curve, range, tolerance)
                };
                return Some((snapped, true, images));
            }
        }
        let plain = self.thread(&on, closed, &fa, &fb, reach, longest);
        let threaded = if self.curved(a).is_some() && self.curved(b).is_some() {
            let between = between_both(&on, &fa, &fb, longest);
            let carried = self.thread(&between, closed, &fa, &fb, reach, longest);
            match (plain, carried) {
                (Some(found), Some(closer)) if closer.2 < found.2 => Some(closer),
                (None, closer) => closer,
                (found, _) => found,
            }
        } else {
            plain
        };
        let (curve, range, tolerance) = threaded?;
        Some((
            if closed {
                Snapped::Loop(curve, range, tolerance)
            } else {
                Snapped::Open(curve, range, tolerance)
            },
            true,
            Vec::new(),
        ))
    }

    /// A curve threaded through a chain's points in the chart of face `g`,
    /// a cylinder, cone, sphere or torus, and lifted onto it. The chain is
    /// cut where it turns by more than [`CHART_CORNER`], and each piece's
    /// chart positions joined by a cubic (a piece of two or three points by
    /// straight chart segments). The curve, its image on `g` and its image
    /// on the patch `h` (the feet of its points there) are interpolated
    /// through [`CHORD_SPLIT`] points a span at the same parameters, and the
    /// pieces joined end to end. Its range, its tolerance (how far it
    /// strays from `h`, or either image's lift from it), and both images
    /// with how far each strays. `None` where `g` is no such surface, a
    /// point has no chart position, a loop does not close in the chart (it
    /// goes round the axis), the curve hooks back on itself, or it strays
    /// past the reach.
    fn chart_thread(
        &self,
        on: &[Point],
        closed: bool,
        g: usize,
        h: usize,
        reach: f64,
    ) -> Option<(Curve, (f64, f64), f64, Images)> {
        use ogeom_geom::{Curve2d as _, Curve3d as _};
        type Spline = (ogeom_math::KnotVector, Vec<Point>);
        let curved = self.curved(g)?;
        if curved.patch.is_some()
            || matches!(curved.shape, Canonical::Swept(_) | Canonical::Plane(_))
        {
            return None;
        }
        let n = on.len();
        if n < 2 {
            return None;
        }
        let (pu, pv) = periodic(&curved.shape);
        let near = |x: f64, c: f64, wraps: bool| {
            if wraps {
                c + ogeom_math::elementary::wrap_signed_angle(x - c)
            } else {
                x
            }
        };
        let mut uv: Vec<Point> = Vec::with_capacity(n);
        for &p in on {
            let (u, v) = match uv.last() {
                None => unwrapped(curved, p, self.tol)?,
                Some(last) => {
                    let (u, v) = chart(&curved.shape, p, self.tol)?;
                    (near(u, last.x, pu), near(v, last.y, pv))
                }
            };
            uv.push(Point::new(u, v, 0.0));
        }
        if closed && uv[0].distance(uv[n - 1]) > 1e-9 {
            return None;
        }
        let turns = |i: usize| -> bool {
            let before = if i == 0 { on[n - 2] } else { on[i - 1] };
            let (a, b) = (on[i] - before, on[i + 1] - on[i]);
            let m = a.magnitude() * b.magnitude();
            m > 0.0 && a.dot(b) < m * CHART_CORNER.cos()
        };
        let mut cuts: Vec<usize> = vec![0];
        cuts.extend((1..n - 1).filter(|&i| turns(i)));
        cuts.push(n - 1);
        // A loop without a corner is carried on past its ends and cut
        // back, as in `thread`, so it runs on smoothly through its start.
        let smooth_loop = closed && cuts.len() == 2 && n >= 4 && !turns(0);
        let patch = self.curved(h)?;
        let other = self.signed(h)?;
        let mut joined: Option<[Spline; 3]> = None;
        // Each dense parameter of the whole curve, with the chain's span
        // it falls in.
        let mut dense_all: Vec<(f64, usize)> = Vec::new();
        for piece in cuts.windows(2) {
            let (s, e) = (piece[0], piece[1]);
            let count = e - s + 1;
            let (pts, chart_pts, pad) = if smooth_loop {
                let pad = 3.min(n / 3);
                let wrap = |all: &[Point]| -> Vec<Point> {
                    all[n - 1 - pad..n - 1]
                        .iter()
                        .chain(all)
                        .chain(&all[1..=pad])
                        .copied()
                        .collect()
                };
                (wrap(on), wrap(&uv), pad)
            } else {
                (on[s..=e].to_vec(), uv[s..=e].to_vec(), 0)
            };
            let parameters =
                crate::fit::spaced(&pts, crate::fit::Spacing::Centripetal, self.tol).ok()?;
            let degree = if count >= 4 { 3 } else { 1 };
            let mut flat =
                crate::fit::interpolate_at(&chart_pts, &parameters, degree, self.tol).ok()?;
            if pad > 0 {
                flat = flat.split_at(parameters[pad], self.tol).ok()?.1;
                flat = flat.split_at(parameters[pad + count - 1], self.tol).ok()?.0;
            }
            let flat: Curve = flat.into();
            let own = &parameters[pad..pad + count];
            let mut dense = vec![own[0]];
            for w in own.windows(2) {
                dense.extend(
                    (1..=CHORD_SPLIT)
                        .map(|k| w[0] + (w[1] - w[0]) * f64::from(k) / f64::from(CHORD_SPLIT)),
                );
            }
            let mut at: Vec<Point> = Vec::with_capacity(dense.len());
            let mut lifted: Vec<Point> = Vec::with_capacity(dense.len());
            let mut feet: Vec<Point> = Vec::with_capacity(dense.len());
            for &t in &dense {
                let c = flat.point_at(t, self.tol).ok()?;
                at.push(Point::new(c.x, c.y, 0.0));
                let p = evaluate(&curved.shape, (c.x, c.y));
                lifted.push(p);
                let (u, v) = chart(&patch.shape, p, self.tol)?;
                feet.push(Point::new(u, v, 0.0));
            }
            let image = crate::fit::interpolate_at(&at, &dense, 3, self.tol).ok()?;
            let foot_image = crate::fit::interpolate_at(&feet, &dense, 3, self.tol).ok()?;
            let curve = crate::fit::interpolate_at(&lifted, &dense, 3, self.tol).ok()?;
            let unweighted = |spline: &ogeom_geom::BSplineCurve| {
                (
                    spline.knots().clone(),
                    spline
                        .control_points()
                        .iter()
                        .map(|c| c.scaled)
                        .collect::<Vec<Point>>(),
                )
            };
            let piece = [
                unweighted(&curve),
                unweighted(&image),
                unweighted(&foot_image),
            ];
            let offset = joined.as_ref().map_or(0.0, |[c, ..]| c.0.domain_end()) - dense[0];
            dense_all.extend(
                dense
                    .iter()
                    .enumerate()
                    .skip(usize::from(!dense_all.is_empty()))
                    .map(|(k, &t)| {
                        let span = k.saturating_sub(1) / CHORD_SPLIT as usize;
                        (t + offset, s + span.min(count - 2))
                    }),
            );
            joined = Some(match joined {
                None => piece,
                Some(before) => [
                    ogeom_math::bspline::join(&before[0], &piece[0]).ok()?,
                    ogeom_math::bspline::join(&before[1], &piece[1]).ok()?,
                    ogeom_math::bspline::join(&before[2], &piece[2]).ok()?,
                ],
            });
        }
        let [(knots, control), image, foot_image] = joined?;
        let curve: Curve = ogeom_geom::BSplineCurve::new(knots, control, self.tol)
            .ok()?
            .into();
        let planar = |(knots, control): Spline| -> Option<PlanarCurve> {
            Some(
                ogeom_geom::BSpline2d::new(
                    knots,
                    control.iter().map(|c| Point2::new(c.x, c.y)).collect(),
                    self.tol,
                )
                .ok()?
                .into(),
            )
        };
        let (image, foot_image) = (planar(image)?, planar(foot_image)?);
        let range = curve.domain();
        let (mut tolerance, mut deviation, mut foot_deviation) = (
            self.tol.confusion(),
            self.tol.confusion() * 1e-2,
            self.tol.confusion() * 1e-2,
        );
        for pair in dense_all.windows(2) {
            let span = pair[1].1;
            for f in [0.0, 0.25, 0.5, 0.75] {
                let t = (pair[0].0 + (pair[1].0 - pair[0].0) * f).clamp(range.0, range.1);
                let p = curve.point_at(t, self.tol).ok()?;
                let c = image.point_at(t, self.tol).ok()?;
                tolerance = tolerance.max(other(p).abs());
                deviation = deviation.max(evaluate(&curved.shape, (c.x, c.y)).distance(p));
                let c = foot_image.point_at(t, self.tol).ok()?;
                foot_deviation = foot_deviation.max(evaluate(&patch.shape, (c.x, c.y)).distance(p));
                // At a corner the curve's direction is either piece's.
                if f > 0.0 && curve.d1_at(t, self.tol).ok()?.dot(on[span + 1] - on[span]) <= 0.0 {
                    return None;
                }
            }
        }
        // Only a seam held within the reach: a patch that meets the surface
        // looser than that does not run out into it, and the chord through
        // the vertices serves as well.
        let tolerance = tolerance.max(deviation).max(foot_deviation);
        (tolerance <= reach).then_some((
            curve,
            range,
            tolerance,
            vec![(g, image, deviation), (h, foot_image, foot_deviation)],
        ))
    }

    /// A curve threaded through points on two faces, as [`Self::chord`]
    /// threads it: its range and how far it strays from either face.
    fn thread(
        &self,
        on: &[Point],
        closed: bool,
        fa: &dyn Fn(Point) -> f64,
        fb: &dyn Fn(Point) -> f64,
        reach: f64,
        longest: f64,
    ) -> Option<(Curve, (f64, f64), f64)> {
        use ogeom_geom::Curve3d as _;
        // A cubic through the points first (a parabola through three);
        // where it overshoots between them (a short span beside long ones),
        // the polyline through them.
        for degree in [3.min(on.len().saturating_sub(1)), 1] {
            let (curve, samples): (Curve, Vec<f64>) = if on.len() == 2 {
                let length = on[0].distance(on[1]);
                let line: Curve = LineCurve::segment(on[0], on[1], self.tol).ok()?.into();
                (
                    line,
                    (0..=16).map(|k| length * f64::from(k) / 16.0).collect(),
                )
            } else {
                // A loop is carried on past its ends and cut back, as a
                // fitted section is.
                let n = on.len();
                let pad = if closed { 3.min(n / 3) } else { 0 };
                let padded: Vec<Point> = if pad > 0 {
                    on[n - 1 - pad..n - 1]
                        .iter()
                        .chain(on)
                        .chain(&on[1..=pad])
                        .copied()
                        .collect()
                } else {
                    on.to_vec()
                };
                let parameters =
                    crate::fit::spaced(&padded, crate::fit::Spacing::Centripetal, self.tol).ok()?;
                let mut spline =
                    crate::fit::interpolate_at(&padded, &parameters, degree, self.tol).ok()?;
                if pad > 0 {
                    spline = spline.split_at(parameters[pad], self.tol).ok()?.1;
                    spline = spline.split_at(parameters[pad + n - 1], self.tol).ok()?.0;
                }
                let own = &parameters[pad..pad + n];
                let samples = own
                    .windows(2)
                    .flat_map(|w| {
                        [
                            w[0],
                            w[0] + (w[1] - w[0]) * 0.25,
                            w[0] + (w[1] - w[0]) * 0.5,
                            w[0] + (w[1] - w[0]) * 0.75,
                        ]
                    })
                    .chain(std::iter::once(own[n - 1]))
                    .collect();
                (spline.into(), samples)
            };
            let range = curve.domain();
            let mut tolerance = self.tol.confusion();
            // The curve runs the way the points it was threaded through do,
            // or it hooks back on itself between them.
            let mut runs = true;
            for (i, &t) in samples.iter().enumerate() {
                let t = t.clamp(range.0, range.1);
                let p = curve.point_at(t, self.tol).ok()?;
                tolerance = tolerance.max(fa(p).abs()).max(fb(p).abs());
                let span = (i / 4).min(on.len() - 2);
                runs &= curve.d1_at(t, self.tol).ok()?.dot(on[span + 1] - on[span]) > 0.0;
            }
            if tolerance > reach.max(longest * CHORD_SAG) {
                return None;
            }
            if !runs {
                continue;
            }
            return Some((curve, range, tolerance));
        }
        None
    }

    /// The curve two faces meet along, where it is no parallel or ruling of
    /// either: vertices of the chain, and points between them, solved onto
    /// both surfaces, and the curve interpolated through them. Its image on
    /// each curved face is interpolated through the same points' chart
    /// positions at the same parameters, so the two run together. `None`
    /// where the surfaces meet tangentially, the solve does not settle, or
    /// the curve strays from either face past the reach.
    ///
    /// A free chain, with one curved face and nothing across, is fitted
    /// the same way through the feet of its vertices on that face's
    /// surface; its tolerance also holds every vertex of the chain, so the
    /// edge keeps to the mesh's boundary within it.
    fn section(
        &self,
        pts: &[Point],
        closed: bool,
        faces: &[usize],
        reach: f64,
    ) -> Option<(Snapped, bool, Images)> {
        // Made through more points until it keeps to the surfaces to a
        // fiftieth of a micron's worth of confusion distances, as close as
        // an exact edge's image would, or the points run out.
        let close = self.tol.confusion() * 50.0;
        // A free chain's points are its own vertices, however many are
        // asked for.
        if faces.len() == 1 {
            return self
                .section_through(pts, closed, faces, reach, SECTION_POINTS)
                .map(|(_, found)| found);
        }
        let mut best: Option<(f64, (Snapped, bool, Images))> = None;
        let mut count = SECTION_POINTS;
        while count <= SECTION_POINTS * 8 {
            let Some((worst, found)) = self.section_through(pts, closed, faces, reach, count)
            else {
                break;
            };
            let done = worst <= close;
            if best.as_ref().is_none_or(|(held, _)| worst < *held) {
                best = Some((worst, found));
            }
            if done {
                break;
            }
            count *= 2;
        }
        best.map(|(_, found)| found)
    }

    /// One fitted section through at least `count` points, with the worst
    /// of its own stray and its images'.
    fn section_through(
        &self,
        pts: &[Point],
        closed: bool,
        faces: &[usize],
        reach: f64,
        count: usize,
    ) -> Option<(f64, (Snapped, bool, Images))> {
        use ogeom_geom::Curve3d as _;
        let (fa, fb) = match faces[..] {
            [a, b] => (self.signed(a)?, Some(self.signed(b)?)),
            [a] => (self.signed(a)?, None),
            _ => return None,
        };
        // A point on the faces near a start: solved onto both, or for a
        // free chain the foot on its face.
        let onto = |start: Point, limit: f64| -> Option<Point> {
            match &fb {
                Some(fb) => onto_both(&fa, fb, start, reach, limit),
                None => {
                    let shape = &self.curved(faces[0])?.shape;
                    let foot = evaluate(shape, chart(shape, start, self.tol)?);
                    (fa(foot).abs() <= reach * 1e-3 && foot.distance(start) <= limit)
                        .then_some(foot)
                }
            }
        };
        let off = |p: Point| fa(p).abs().max(fb.as_ref().map_or(0.0, |fb| fb(p).abs()));
        let steps = if closed { pts.len() } else { pts.len() - 1 };
        // A long chain is taken a few vertices at a time: the curve needs
        // its shape, not every vertex. A free chain is taken through every
        // vertex and nothing between: the foot of a point on a chord of a
        // curved boundary lies off the boundary, across the surface, by as
        // much as the chord cuts off it. A free chain of one or two spans
        // has too few vertices for a cubic, and its spans are split.
        let free = fb.is_none();
        let stride = if free {
            1
        } else {
            steps.div_ceil(SECTION_SPANS).max(1)
        };
        let split = if free {
            if steps >= 3 { 1 } else { SECTION_SPLIT }
        } else {
            count.div_ceil(steps.div_ceil(stride)).max(SECTION_SPLIT)
        };
        // A chord of the mesh lies across a face recognized from it, and a
        // face recognized from a mesh turns through well under a right
        // angle from one end of one of its chords to the other. Where the
        // fitted surface turns farther, it is not the surface the chord
        // lies on, and the curve solved on it runs where the mesh does not.
        let turns_away = |p: Point, q: Point| {
            faces.iter().filter_map(|&g| self.curved(g)).any(|curved| {
                let (gp, gq) = (gradient(&curved.shape, p), gradient(&curved.shape, q));
                gp.dot(gq) <= 0.0
            })
        };
        let mut on = Vec::new();
        let mut i = 0;
        while i < steps {
            let next = (i + stride).min(steps);
            let (p, q) = (pts[i], pts[next % pts.len()]);
            if turns_away(p, q) {
                return None;
            }
            for k in 0..split {
                #[allow(clippy::cast_precision_loss, reason = "a handful of splits")]
                let f = k as f64 / split as f64;
                let limit = if k == 0 { reach } else { p.distance(q) + reach };
                on.push(onto(p + (q - p) * f, limit)?);
            }
            i = next;
        }
        on.push(if closed {
            on[0]
        } else {
            onto(pts[pts.len() - 1], reach)?
        });
        // A loop is interpolated with a few of its points carried on past
        // each end, then cut back to its own: an open interpolation left
        // free at its ends wanders where the loop meets itself.
        let pad = if closed {
            SECTION_PAD.min(on.len() / 4)
        } else {
            0
        };
        let n = on.len();
        let padded: Vec<Point> = if pad > 0 {
            on[n - 1 - pad..n - 1]
                .iter()
                .chain(&on)
                .chain(&on[1..=pad])
                .copied()
                .collect()
        } else {
            on.clone()
        };
        let parameters =
            crate::fit::spaced(&padded, crate::fit::Spacing::Centripetal, self.tol).ok()?;
        let cut = (parameters[pad], parameters[pad + n - 1]);
        let trimmed = |spline: ogeom_geom::BSplineCurve| -> Option<ogeom_geom::BSplineCurve> {
            if pad == 0 {
                return Some(spline);
            }
            let (_, after) = spline.split_at(cut.0, self.tol).ok()?;
            Some(after.split_at(cut.1, self.tol).ok()?.0)
        };
        let curve: Curve =
            trimmed(crate::fit::interpolate_at(&padded, &parameters, 3, self.tol).ok()?)?.into();
        let range = curve.domain();
        let parameters = parameters[pad..pad + n].to_vec();
        // Between the points it was made through, how far it strays from
        // either surface.
        let between = |k: usize, f: f64| parameters[k] + (parameters[k + 1] - parameters[k]) * f;
        // Its ends are the chain's ends solved onto both surfaces, a hair
        // from the vertices they meet.
        let last = if closed { pts[0] } else { pts[pts.len() - 1] };
        let mut tolerance = self
            .tol
            .confusion()
            .max(on[0].distance(pts[0]))
            .max(on[on.len() - 1].distance(last));
        for k in 0..on.len() - 1 {
            for f in [0.25, 0.5, 0.75] {
                let p = curve.point_at(between(k, f), self.tol).ok()?;
                tolerance = tolerance.max(off(p));
            }
        }
        if free {
            // Every vertex of a free chain to the foot the curve runs
            // through.
            for (i, p) in pts.iter().enumerate() {
                tolerance = tolerance.max(p.distance(on[i * split]));
            }
        }
        if tolerance > reach {
            return None;
        }
        let mut images = Vec::new();
        for &g in faces {
            let Some(curved) = self.curved(g) else {
                continue;
            };
            let (pu, pv) = periodic(&curved.shape);
            let mut uv: Vec<Point> = Vec::with_capacity(padded.len());
            for p in &padded {
                let (u, v) = match uv.last() {
                    None => unwrapped(curved, *p, self.tol)?,
                    Some(last) => {
                        let (u, v) = chart(&curved.shape, *p, self.tol)?;
                        let near = |x: f64, c: f64, wraps: bool| {
                            if wraps {
                                c + ogeom_math::elementary::wrap_signed_angle(x - c)
                            } else {
                                x
                            }
                        };
                        (near(u, last.x, pu), near(v, last.y, pv))
                    }
                };
                uv.push(Point::new(u, v, 0.0));
            }
            let all =
                crate::fit::spaced(&padded, crate::fit::Spacing::Centripetal, self.tol).ok()?;
            let flat = trimmed(crate::fit::interpolate_at(&uv, &all, 3, self.tol).ok()?)?;
            let control: Vec<Point2> = flat
                .control_points()
                .iter()
                .map(|c| Point2::new(c.scaled.x, c.scaled.y))
                .collect();
            let pcurve: PlanarCurve =
                ogeom_geom::BSpline2d::new(flat.knots().clone(), control, self.tol)
                    .ok()?
                    .into();
            let mut deviation = self.tol.confusion();
            for k in 0..on.len() - 1 {
                for f in [0.0, 0.25, 0.5, 0.75] {
                    use ogeom_geom::Curve2d as _;
                    let t = between(k, f);
                    let at = pcurve.point_at(t, self.tol).ok()?;
                    let lifted = evaluate(&curved.shape, (at.x, at.y));
                    deviation = deviation.max(lifted.distance(curve.point_at(t, self.tol).ok()?));
                }
            }
            if deviation > reach {
                return None;
            }
            images.push((g, pcurve, deviation));
        }
        let worst = images.iter().map(|(_, _, d)| *d).fold(tolerance, f64::max);
        Some((
            worst,
            (
                if closed {
                    Snapped::Loop(curve, range, tolerance)
                } else {
                    Snapped::Open(curve, range, tolerance)
                },
                true,
                images,
            ),
        ))
    }

    /// A face's surface as a signed distance, where it has one.
    fn signed(&self, g: usize) -> Option<Box<dyn Fn(Point) -> f64>> {
        match &self.groups.carriers[g] {
            Carrier::Plane(plane) => {
                let plane = *plane;
                Some(Box::new(move |p: Point| plane.signed_distance_to(p)))
            }
            Carrier::Curved(c) => {
                let shape = c.shape.clone();
                Some(Box::new(move |p: Point| shape.signed_distance_to(p)))
            }
            Carrier::Gone => None,
        }
    }

    /// A candidate curve held against the chain and the faces: its range,
    /// which way it runs along the chain, and its tolerance.
    fn fitted(
        &self,
        curve: Curve,
        pts: &[Point],
        closed: bool,
        faces: &[usize],
        reach: f64,
    ) -> Option<(Snapped, bool)> {
        use ogeom_geom::Curve3d as _;
        let tau = core::f64::consts::TAU;
        let parameter = |p: Point| -> Option<f64> {
            match &curve {
                Curve::Line(l) => {
                    let axis = l.axis();
                    Some((p - axis.location).dot(axis.direction.vector()))
                }
                Curve::Circle(c) => {
                    ogeom_math::elementary::circle_parameter(&c.circle(), p, self.tol).ok()
                }
                _ => None,
            }
        };
        let mut tolerance = self.tol.confusion();
        let mut ts = Vec::with_capacity(pts.len());
        for p in pts {
            let t = parameter(*p)?;
            let on = curve.point_at(t, self.tol).ok()?;
            tolerance = tolerance.max(on.distance(*p));
            ts.push(t);
        }
        // The sweep along the chain, unwrapped for a circle.
        let mut sweep = 0.0;
        let steps = if closed { ts.len() } else { ts.len() - 1 };
        for i in 0..steps {
            let (a, b) = (ts[i], ts[(i + 1) % ts.len()]);
            let d = b - a;
            sweep += if matches!(curve, Curve::Circle(_)) {
                ogeom_math::elementary::wrap_signed_angle(d)
            } else {
                d
            };
        }
        let (snapped, forward) = if closed {
            if !matches!(curve, Curve::Circle(_)) || (sweep.abs() - tau).abs() > 1e-3 {
                return None;
            }
            (None, sweep > 0.0)
        } else if sweep > 0.0 {
            (Some((ts[0], ts[0] + sweep)), true)
        } else {
            let last = ts[ts.len() - 1];
            (Some((last, last - sweep)), false)
        };
        let range = snapped.unwrap_or((0.0, tau));
        if range.1 - range.0 <= self.tol.parametric() {
            return None;
        }
        // Drawn in every plane it bounds: a circle is imaged in a plane by
        // projection, which a circle standing across the plane has not.
        for &g in faces {
            if let Carrier::Plane(plane) = &self.groups.carriers[g] {
                let surface: ogeom_geom::SurfaceGeometry = PlaneSurface::new(*plane).into();
                ogeom_intersect::exact_pcurve_of(&curve, &surface, self.tol)?;
            }
        }
        // Standing on every face it bounds.
        for k in 0..=16 {
            let t = range.0 + (range.1 - range.0) * f64::from(k) / 16.0;
            let p = curve.point_at(t, self.tol).ok()?;
            for &g in faces {
                let off = match &self.groups.carriers[g] {
                    Carrier::Plane(plane) => plane.signed_distance_to(p).abs(),
                    Carrier::Curved(c) => c.shape.distance_to(p),
                    Carrier::Gone => return None,
                };
                tolerance = tolerance.max(off);
            }
        }
        if tolerance > reach {
            return None;
        }
        Some((
            match snapped {
                Some(range) => Snapped::Open(curve, range, tolerance),
                None => Snapped::Closed(curve, tolerance),
            },
            forward,
        ))
    }
}

#[derive(Clone)]
pub(super) enum Snapped {
    Open(Curve, (f64, f64), f64),
    Closed(Curve, f64),
    /// A closed curve that is no circle, starting and ending at the
    /// chain's first vertex.
    Loop(Curve, (f64, f64), f64),
}

/// How many pieces each span of a chain is cut into for a fitted section:
/// its vertices alone leave a coarse chain's curve unconstrained between
/// them.
const SECTION_SPLIT: usize = 4;

/// The most spans a chain is taken in for a fitted section.
const SECTION_SPANS: usize = 40;

/// How many points a closed section is carried on past each end while it
/// is interpolated.
const SECTION_PAD: usize = 8;

/// How many points a fitted section is interpolated through, at least:
/// a short chain's spans are cut finer to reach it.
const SECTION_POINTS: usize = 160;

/// A fitted section's images on the curved faces it bounds: the face, the
/// image, and how far the image strays from the curve.
pub(super) type Images = Vec<(usize, PlanarCurve, f64)>;

/// A point solved onto where two surfaces meet, from a start near both:
/// Newton's step, the least one that zeroes both signed distances to first
/// order. `None` where the surfaces meet tangentially there, where the
/// solve does not settle within a thousandth of `reach` of both, or where
/// it lands farther than `limit` from its start.
fn onto_both(
    fa: &dyn Fn(Point) -> f64,
    fb: &dyn Fn(Point) -> f64,
    start: Point,
    reach: f64,
    limit: f64,
) -> Option<Point> {
    let mut p = start;
    let h = 1e-7 * (1.0 + p.to_vector().magnitude());
    let gradient = |f: &dyn Fn(Point) -> f64, p: Point| {
        let d = |v: Vector| (f(p + v * h) - f(p - v * h)) / (2.0 * h);
        Vector::new(d(Vector::X), d(Vector::Y), d(Vector::Z))
    };
    let (mut last, mut stalled) = (f64::INFINITY, 0);
    for _ in 0..40 {
        let (va, vb) = (fa(p), fb(p));
        let residual = va.abs().max(vb.abs());
        if residual <= 1e-13 * (1.0 + p.to_vector().magnitude()) {
            break;
        }
        // A surface whose distance is itself a search answers to its
        // search's precision, short of the residual asked: two steps that
        // do not halve it have reached that floor.
        stalled = if residual > last * 0.5 {
            stalled + 1
        } else {
            0
        };
        if stalled >= 2 {
            break;
        }
        last = residual;
        let (ga, gb) = (gradient(fa, p), gradient(fb, p));
        let (aa, ab, bb) = (ga.dot(ga), ga.dot(gb), gb.dot(gb));
        let det = aa.mul_add(bb, -(ab * ab));
        if det <= 1e-12 * aa * bb {
            return None;
        }
        let la = (va * bb - vb * ab) / det;
        let lb = (vb * aa - va * ab) / det;
        p = p - ga * la - gb * lb;
    }
    (fa(p).abs().max(fb(p).abs()) <= reach * 1e-3 && p.distance(start) <= limit).then_some(p)
}

/// How far a chain threaded in a canonical surface's chart turns at a
/// vertex, past straight on, for the curve to take a corner there: thirty
/// degrees, a crease's turn.
const CHART_CORNER: f64 = core::f64::consts::FRAC_PI_6;

/// How many points each span of a chord is cut into, at most, when its
/// points are carried onto both faces.
const CHORD_SPLIT: u32 = 8;

/// A chain of points on two surfaces that meet all but tangentially, with
/// points between them: each span is cut into pieces no longer than an
/// eighth of the longest, and each cut carried onto both surfaces in turn
/// until it rests on both, where that brings it nearer them than the
/// straight span. A curve threaded through the cuts as well as the chain's
/// points keeps to the surfaces between those points, where one threaded
/// through the points alone stands off them by the span's sag.
fn between_both(
    on: &[Point],
    fa: &dyn Fn(Point) -> f64,
    fb: &dyn Fn(Point) -> f64,
    longest: f64,
) -> Vec<Point> {
    let off = |p: Point| fa(p).abs().max(fb(p).abs());
    let gradient = |f: &dyn Fn(Point) -> f64, p: Point| {
        let h = 1e-7 * (1.0 + p.to_vector().magnitude());
        let d = |v: Vector| (f(p + v * h) - f(p - v * h)) / (2.0 * h);
        Vector::new(d(Vector::X), d(Vector::Y), d(Vector::Z))
    };
    let carried = |start: Point| -> Point {
        let mut p = start;
        for _ in 0..32 {
            for f in [fa, fb] {
                let g = gradient(f, p);
                let m = g.dot(g);
                if m > 0.0 {
                    p = p - g * (f(p) / m);
                }
            }
            if off(p) <= 1e-12 * (1.0 + p.to_vector().magnitude()) {
                break;
            }
        }
        // A point carried along the tangency, farther than twice the gap
        // it closes, is not where the span runs.
        if p.is_finite() && off(p) < off(start) && p.distance(start) <= off(start) * 2.0 {
            p
        } else {
            start
        }
    };
    let mut out = Vec::with_capacity(on.len() * CHORD_SPLIT as usize);
    for w in on.windows(2) {
        out.push(w[0]);
        #[allow(
            clippy::cast_possible_truncation,
            clippy::cast_sign_loss,
            reason = "a handful of pieces"
        )]
        let pieces = ((w[0].distance(w[1]) * f64::from(CHORD_SPLIT) / longest).ceil() as u32)
            .clamp(1, CHORD_SPLIT);
        for k in 1..pieces {
            let f = f64::from(k) / f64::from(pieces);
            out.push(carried(w[0].lerp(w[1], f)));
        }
    }
    out.extend(on.last().copied());
    out
}

/// The curves a chain on a curved surface may be: the parallel circle
/// through its mean height, when every point sits at one height and one
/// distance from the axis; and, on a cylinder or a cone, the ruling
/// through its mean angle, when the chain is straight along it. Each is
/// placed on the surface exactly, its angle measured from the surface's
/// own origin.
fn candidates(
    shape: &Canonical,
    pts: &[Point],
    across: Option<Plane>,
    reach: f64,
    tol: Tolerances,
) -> Vec<Curve> {
    let mut out = Vec::new();
    let plane_of = |pts: &[Point]| match across {
        Some(plane) => Some((plane.project(pts[0]), plane.normal())),
        None => plane_through(pts, tol),
    };
    if let Canonical::Sphere(sphere) = shape {
        // The section of the sphere by the chain's plane. Where that plane
        // is square to the sphere's axis the section is a latitude, and is
        // built on the sphere's own frame, so its parameter is the sphere's
        // angle and it starts where the sphere's seam does.
        let axis = sphere.frame().z();
        if pts.len() >= 3
            && let Some((centre, normal)) = plane_of(pts)
            && normal.vector().cross(axis.vector()).magnitude()
                <= if across.is_some() { 1e-12 } else { 1e-3 }
        {
            let h = (centre - sphere.centre()).dot(axis.vector());
            let r2 = sphere.radius().powi(2) - h * h;
            if r2 > 0.0
                && let Ok(frame) = Frame::new(
                    sphere.centre() + axis.vector() * h,
                    axis,
                    sphere.frame().x(),
                    tol,
                )
                && let Ok(circle) = ogeom_math::Circle::new(frame, r2.sqrt(), tol)
            {
                out.push(ogeom_geom::CircleCurve::new(circle).into());
                return out;
            }
        }
        if pts.len() >= 3
            && let Some((centre, normal)) = plane_of(pts)
        {
            let n = normal.vector();
            let d = (centre - sphere.centre()).dot(n);
            let r2 = sphere.radius().powi(2) - d * d;
            if r2 > 0.0
                && let Ok(frame) = Frame::new(
                    sphere.centre() + n * d,
                    normal,
                    normal.any_perpendicular(),
                    tol,
                )
                && let Ok(circle) = ogeom_math::Circle::new(frame, r2.sqrt(), tol)
            {
                out.push(ogeom_geom::CircleCurve::new(circle).into());
            }
        }
        return out;
    }
    let Some(frame) = axis_frame(shape) else {
        return out;
    };
    let (o, z) = (frame.origin(), frame.z().vector());
    let heights: Vec<f64> = pts.iter().map(|p| (*p - o).dot(z)).collect();
    let radii: Vec<f64> = pts
        .iter()
        .zip(&heights)
        .map(|(p, h)| ((*p - o) - z * *h).magnitude())
        .collect();
    #[allow(
        clippy::cast_precision_loss,
        reason = "chain lengths are far below 2^52"
    )]
    let count = pts.len() as f64;
    let mean_h = heights.iter().sum::<f64>() / count;
    let mean_r = radii.iter().sum::<f64>() / count;
    let level = heights.iter().all(|h| (h - mean_h).abs() <= reach)
        && radii.iter().all(|r| (r - mean_r).abs() <= reach);
    if level {
        let radius = match shape {
            Canonical::Cylinder(c) => c.radius(),
            Canonical::Cone(c) => c.radius_at(mean_h),
            Canonical::Torus(t) => {
                let (big, small) = (t.major_radius(), t.minor_radius());
                let off = (small * small - mean_h * mean_h).max(0.0).sqrt();
                // At the tube's crest the circle's radius turns infinitely
                // fast with its height, and a few slops of height move it by
                // many: a ring within the reach of the crest is the crest.
                if (small - mean_h.abs()).abs() <= reach {
                    big
                } else if (big + off - mean_r).abs() <= (big - off - mean_r).abs() {
                    big + off
                } else {
                    big - off
                }
            }
            _ => mean_r,
        };
        if let Ok(at) = Frame::new(o + z * mean_h, frame.z(), frame.x(), tol)
            && let Ok(circle) = ogeom_math::Circle::new(at, radius, tol)
        {
            out.push(ogeom_geom::CircleCurve::new(circle).into());
        }
    }
    // A circle of a torus's tube: every point in one plane through the
    // axis, at the angle the chain stands at. Built with the tube's own
    // angle as its parameter, starting at the outer equator.
    if let Canonical::Torus(torus) = shape
        && pts.len() >= 2
    {
        let mut angles: Vec<f64> = pts
            .iter()
            .filter_map(|x| chart(shape, *x, tol).map(|c| c.0))
            .collect();
        if !angles.is_empty() {
            let (u, _) = angular_spread(&mut angles);
            let (x, y) = (frame.x().vector(), frame.y().vector());
            let out_u = x * u.cos() + y * u.sin();
            let across = y * u.cos() - x * u.sin();
            let in_plane = pts.iter().all(|p| (*p - o).dot(across).abs() <= reach);
            if in_plane
                && let Ok(radial) = Direction::new(out_u, tol)
                && let Ok(normal) = Direction::new(out_u.cross(z), tol)
                && let Ok(at) = Frame::new(o + out_u * torus.major_radius(), normal, radial, tol)
                && let Ok(circle) = ogeom_math::Circle::new(at, torus.minor_radius(), tol)
            {
                out.push(ogeom_geom::CircleCurve::new(circle).into());
            }
        }
    }
    if matches!(shape, Canonical::Cylinder(_) | Canonical::Cone(_)) && pts.len() >= 2 {
        let (p, q) = (pts[0], pts[pts.len() - 1]);
        let straight = pts.iter().all(|x| distance_to_line(*x, p, q) <= reach);
        let mut angles: Vec<f64> = pts
            .iter()
            .filter_map(|x| chart(shape, *x, tol).map(|c| c.0))
            .collect();
        if straight && !angles.is_empty() {
            let (u, _) = angular_spread(&mut angles);
            let (lo, hi) = (
                heights.iter().copied().fold(f64::INFINITY, f64::min),
                heights.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            );
            let (a, b) = (evaluate(shape, (u, lo)), evaluate(shape, (u, hi)));
            if let Ok(line) = LineCurve::segment(a, b, tol) {
                // Measured from the chain's first end, so its parameter is
                // the chain's distance along it.
                let _ = line;
                let (a, b) = if (pts[0] - a).magnitude() <= (pts[0] - b).magnitude() {
                    (a, b)
                } else {
                    (b, a)
                };
                if let Ok(line) = LineCurve::segment(a, b, tol) {
                    out.push(line.into());
                }
            }
        }
    }
    out
}

/// The plane nearest a set of points: their centroid and the covariance's
/// smallest direction.
pub(super) fn plane_through(points: &[Point], tol: Tolerances) -> Option<(Point, Direction)> {
    #[allow(
        clippy::cast_precision_loss,
        reason = "chain lengths are far below 2^52"
    )]
    let count = points.len() as f64;
    let c = points.iter().fold(Vector::ZERO, |s, p| s + p.to_vector()) / count;
    let mut m = nalgebra::Matrix3::<f64>::zeros();
    for p in points {
        let d = p.to_vector() - c;
        let v = nalgebra::Vector3::new(d.x, d.y, d.z);
        m += v * v.transpose();
    }
    let eigen = nalgebra::SymmetricEigen::new(m);
    let mut best = 0;
    for i in 1..3 {
        if eigen.eigenvalues[i] < eigen.eigenvalues[best] {
            best = i;
        }
    }
    let v = eigen.eigenvectors.column(best);
    Some((
        Point::from_vector(c),
        Direction::new(Vector::new(v[0], v[1], v[2]), tol).ok()?,
    ))
}

/// The surface a curved face is built on, windowed along its axis to hold
/// every point its boundary reaches. A cone's window stops short of its
/// apex, or ends on it where `to_apex` (a cap closing there).
pub(super) fn surface_of(
    curved: &Curved,
    points: &[Point],
    to_apex: bool,
    tol: Tolerances,
) -> OgeomResult<ogeom_geom::SurfaceGeometry> {
    use ogeom_geom::{ConeSurface, CylinderSurface, SphereSurface, TorusSurface};
    let heights = || {
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for &v in &curved.vertices {
            if let Some((_, h)) = chart(&curved.shape, points[v as usize], tol) {
                lo = lo.min(h);
                hi = hi.max(h);
            }
        }
        let margin = (hi - lo).mul_add(0.25, tol.confusion() * 10.0);
        (lo - margin, hi + margin)
    };
    Ok(match &curved.shape {
        Canonical::Cylinder(c) => CylinderSurface::new(*c, heights())?.into(),
        Canonical::Cone(c) => {
            // Short of the apex, where the cone's radius runs out.
            let (lo, hi) = heights();
            let apex = -c.reference_radius() / c.half_angle().tan();
            let lo = if to_apex {
                apex
            } else {
                lo.max(apex + (hi - apex) * 1e-6)
            };
            ConeSurface::new(*c, (lo, hi))?.into()
        }
        Canonical::Sphere(s) => SphereSurface::new(*s).into(),
        Canonical::Torus(t) => TorusSurface::new(*t).into(),
        Canonical::Plane(p) => PlaneSurface::new(*p).into(),
        Canonical::Swept(s) => s.surface.clone(),
    })
}

/// An edge's image on a curved face, and how far it strays from the edge:
/// the straight chart segment where the edge is a parallel or a ruling of
/// the face, otherwise a fit by projection; `None` past the reach.
pub(super) fn image_on(
    curved: &Curved,
    surface: &ogeom_geom::SurfaceGeometry,
    curve: &Curve,
    range: (f64, f64),
    reach: f64,
    tol: Tolerances,
) -> Option<(PlanarCurve, f64)> {
    if let Some((pcurve, deviation)) = straight_image(curved, curve, range, tol)
        && deviation <= reach
    {
        return Some((pcurve, deviation));
    }
    if let Some((pcurve, deviation)) = interpolated_image(curved, curve, range, tol)
        && deviation <= reach
    {
        return Some((pcurve, deviation));
    }
    let (pcurve, error, _, off, _) =
        crate::pcurve_fit::fit_projected_pcurve_capped(curve, range, surface, reach, tol).ok()?;
    let deviation = error.max(off);
    (deviation <= reach).then_some((pcurve, deviation))
}

/// An edge's image interpolated through the chart positions of points
/// along it, at the edge's own parameters: the curve and its image agree
/// at every sample by construction, and between them to the fourth power
/// of the spacing, which is doubled until the image keeps to the curve
/// within a hundredth of the confusion distance. `None` where a sample
/// has no chart position (a pole).
fn interpolated_image(
    curved: &Curved,
    curve: &Curve,
    range: (f64, f64),
    tol: Tolerances,
) -> Option<(PlanarCurve, f64)> {
    use ogeom_geom::{Curve2d as _, Curve3d as _};
    let (pu, pv) = periodic(&curved.shape);
    let near = |x: f64, c: f64, wraps: bool| {
        if wraps {
            c + ogeom_math::elementary::wrap_signed_angle(x - c)
        } else {
            x
        }
    };
    // On a patch, a curve threaded through a long chain has a knot at
    // every vertex, and samples spread evenly over its range miss the
    // turns between them: the samples are spread over each knot span
    // alike.
    let mut pieces = vec![range.0];
    if let Curve::BSpline(spline) = curve
        && curved.patch.is_some()
    {
        pieces.extend(
            spline
                .knots()
                .distinct()
                .iter()
                .map(|k| k.0)
                .filter(|&t| t > range.0 && t < range.1),
        );
    }
    pieces.push(range.1);
    let spans = u32::try_from(pieces.len() - 1).ok()?;
    let mut best: Option<(PlanarCurve, f64)> = None;
    let mut each = 32_u32.div_ceil(spans);
    while each * spans <= 1024.max(spans * 8) {
        let parameters: Vec<f64> = std::iter::once(range.0)
            .chain(pieces.windows(2).flat_map(|w| {
                (1..=each).map(move |k| w[0] + (w[1] - w[0]) * f64::from(k) / f64::from(each))
            }))
            .collect();
        let mut uv: Vec<Point> = Vec::with_capacity(parameters.len());
        for &t in &parameters {
            let p = curve.point_at(t, tol).ok()?;
            let (u, v) = match uv.last() {
                None => unwrapped(curved, p, tol)?,
                Some(last) => {
                    let (u, v) = chart(&curved.shape, p, tol)?;
                    (near(u, last.x, pu), near(v, last.y, pv))
                }
            };
            uv.push(Point::new(u, v, 0.0));
        }
        let flat = crate::fit::interpolate_at(&uv, &parameters, 3, tol).ok()?;
        let control: Vec<Point2> = flat
            .control_points()
            .iter()
            .map(|c| Point2::new(c.scaled.x, c.scaled.y))
            .collect();
        let pcurve: PlanarCurve = ogeom_geom::BSpline2d::new(flat.knots().clone(), control, tol)
            .ok()?
            .into();
        let mut deviation = tol.confusion() * 1e-2;
        for pair in parameters.windows(2) {
            for f in [0.25, 0.5, 0.75] {
                let t = pair[0] + (pair[1] - pair[0]) * f;
                let at = pcurve.point_at(t, tol).ok()?;
                let lifted = evaluate(&curved.shape, (at.x, at.y));
                deviation = deviation.max(lifted.distance(curve.point_at(t, tol).ok()?));
            }
        }
        let done = deviation <= tol.confusion() * 1e-2;
        if best.as_ref().is_none_or(|(_, held)| deviation < *held) {
            best = Some((pcurve, deviation));
        }
        if done {
            break;
        }
        each *= 2;
    }
    best
}

/// A curve's projection into a plane's chart, interpolated at the curve's
/// own parameters, the spacing halved until the image keeps to the
/// projection within a hundredth of the confusion distance. Whatever the
/// curve stands off the plane, the image lifts to its foot there.
pub(super) fn projected_image(
    plane: Plane,
    curve: &Curve,
    range: (f64, f64),
    tol: Tolerances,
) -> Option<PlanarCurve> {
    use ogeom_geom::{Curve2d as _, Curve3d as _};
    let flat = |t: f64| -> Option<Point2> {
        let (u, v) = ogeom_math::elementary::plane_parameters(&plane, curve.point_at(t, tol).ok()?);
        Some(Point2::new(u, v))
    };
    let mut best: Option<(PlanarCurve, f64)> = None;
    let mut count = 32_u32;
    while count <= 1024 {
        let parameters: Vec<f64> = (0..=count)
            .map(|k| range.0 + (range.1 - range.0) * f64::from(k) / f64::from(count))
            .collect();
        let uv: Vec<Point> = parameters
            .iter()
            .map(|&t| flat(t).map(|q| Point::new(q.x, q.y, 0.0)))
            .collect::<Option<_>>()?;
        let fitted = crate::fit::interpolate_at(&uv, &parameters, 3, tol).ok()?;
        let control: Vec<Point2> = fitted
            .control_points()
            .iter()
            .map(|c| Point2::new(c.scaled.x, c.scaled.y))
            .collect();
        let pcurve: PlanarCurve = ogeom_geom::BSpline2d::new(fitted.knots().clone(), control, tol)
            .ok()?
            .into();
        let mut miss: f64 = 0.0;
        for pair in parameters.windows(2) {
            for f in [0.25, 0.5, 0.75] {
                let t = pair[0] + (pair[1] - pair[0]) * f;
                miss = miss.max(
                    pcurve
                        .point_at(t, tol)
                        .ok()?
                        .square_distance(flat(t)?)
                        .sqrt(),
                );
            }
        }
        if best.as_ref().is_none_or(|(_, held)| miss < *held) {
            best = Some((pcurve, miss));
        }
        if miss <= tol.confusion() * 1e-2 {
            break;
        }
        count *= 2;
    }
    best.map(|(pcurve, _)| pcurve)
}

/// A degree-one pcurve over the edge's range, from the chart points of its
/// ends on the face's branch, and how far it strays along its length.
fn straight_image(
    curved: &Curved,
    curve: &Curve,
    range: (f64, f64),
    tol: Tolerances,
) -> Option<(PlanarCurve, f64)> {
    use ogeom_geom::Curve3d as _;
    let at =
        |t: f64| -> Option<(f64, f64)> { unwrapped(curved, curve.point_at(t, tol).ok()?, tol) };
    let start = at(range.0)?;
    // Carried along the edge a quarter at a time, so a full turn ends a
    // whole period from where it began.
    let (pu, pv) = periodic(&curved.shape);
    let step = |a: f64, b: f64, wraps: bool| {
        if wraps {
            a + ogeom_math::elementary::wrap_signed_angle(b - a)
        } else {
            b
        }
    };
    let mut end = start;
    for k in 1..=4 {
        let next = at(range.0 + (range.1 - range.0) * f64::from(k) / 4.0)?;
        end = (step(end.0, next.0, pu), step(end.1, next.1, pv));
    }
    let pcurve = linear(start, end, range, tol).ok()?;
    let mut deviation: f64 = 0.0;
    for k in 0..=16 {
        let f = f64::from(k) / 16.0;
        let t = range.0 + (range.1 - range.0) * f;
        let uv = (
            (end.0 - start.0).mul_add(f, start.0),
            (end.1 - start.1).mul_add(f, start.1),
        );
        deviation =
            deviation.max(evaluate(&curved.shape, uv).distance(curve.point_at(t, tol).ok()?));
    }
    Some((pcurve, deviation))
}
