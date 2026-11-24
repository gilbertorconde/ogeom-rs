//! Fine rings: a planar piece's rings resampled off their exact curves.

use super::*;

/// How many chords a curved strand of a planar piece is resampled into.
pub(crate) const FINE: usize = 1024;

/// How many of a resampled arc's chords are passed over together where a
/// ray cannot reach them.
pub(crate) const FINE_BLOCK: usize = 32;

/// A planar piece's rings resampled finely off their exact curves. The
/// piece's own polylines bow inside a curved strand by more than the width
/// of a sliver beside it, so a probe picked inside them can stand outside
/// the piece. These rings hold it to a few microns.
///
/// A circular strand's samples are drawn as a ray asks for them: a ray
/// meets a circle in at most two places, and the samples far from its
/// height cannot be crossed by it.
pub(crate) struct FineRings<'a> {
    rings: Vec<Vec<FinePart<'a>>>,
    tol: Tolerances,
}

/// One strand of a [`FineRings`] ring.
pub(crate) enum FinePart<'a> {
    /// Its samples, in the ring's direction.
    Points(Vec<Point2>),
    /// A circle's samples, drawn when first asked for.
    Arc(FineArc<'a>),
}

/// A circular strand's [`FINE`] chords, each sample drawn once, when a ray
/// first needs it.
pub(crate) struct FineArc<'a> {
    curve: &'a PlanarCurve,
    circle: ogeom_math::Circle2,
    from: f64,
    to: f64,
    wrap: Option<(f64, f64)>,
    reversed: bool,
    /// The samples drawn so far; empty until the first is.
    drawn: core::cell::RefCell<Vec<Option<Point2>>>,
}

impl FineArc<'_> {
    /// The parameter sample `i` is drawn at.
    fn parameter(&self, i: usize) -> f64 {
        #[allow(clippy::cast_precision_loss)]
        let t = self.from + (self.to - self.from) * (i as f64 / FINE as f64);
        self.wrap.map_or(t, |domain| fold(t, domain))
    }

    /// Sample `i`. A circle is drawn at any parameter.
    fn point(&self, i: usize, tol: Tolerances) -> Point2 {
        if let Some(Some(p)) = self.drawn.borrow().get(i) {
            return *p;
        }
        let p = self
            .curve
            .point_at(self.parameter(i), tol)
            .unwrap_or(Point2::new(f64::NAN, f64::NAN));
        let mut drawn = self.drawn.borrow_mut();
        if drawn.is_empty() {
            drawn.resize(FINE + 1, None);
        }
        drawn[i] = Some(p);
        p
    }

    /// How far the circle's coordinates stand from its centre's, over and
    /// above its rounding.
    fn slack(&self) -> f64 {
        let centre = self.circle.centre();
        1e-9 * (1.0 + self.circle.radius() + centre.x.abs() + centre.y.abs())
    }

    /// Whether a horizontal ray from `p` may cross a chord between samples
    /// `lo` and `hi`. Every sample between them lies within the arc's length
    /// between their parameters of sample `lo`, so where that reach stays
    /// wholly above or below `p`'s height no chord there straddles it. A
    /// parameter folded across its domain's end inside the block, or a
    /// sample that is not finite, leaves the block to be asked chord by
    /// chord.
    fn may_cross(&self, lo: usize, hi: usize, p: Point2, tol: Tolerances) -> bool {
        let (t0, t1) = (self.parameter(lo), self.parameter(hi));
        #[allow(clippy::cast_precision_loss)]
        let step = (self.to - self.from) * ((hi - lo) as f64 / FINE as f64);
        let unfolded = self
            .wrap
            .is_none_or(|(a, b)| ((t1 - t0) - step).abs() <= 0.5 * (b - a).abs());
        let start = self.point(lo, tol);
        let reach = self.circle.radius() * (t1 - t0).abs() + self.slack();
        if !(unfolded && reach.is_finite() && start.y.is_finite() && p.y.is_finite()) {
            return true;
        }
        start.y - reach <= p.y && start.y + reach >= p.y
    }

    /// The greatest x and the least and greatest y any sample may have,
    /// none drawn: a sample is the centre moved along each axis by the
    /// radius times a cosine or a sine.
    fn extent(&self) -> (f64, f64, f64) {
        let frame = self.circle.frame();
        let (x, y) = (frame.x().vector(), frame.y().vector());
        let r = self.circle.radius();
        let centre = self.circle.centre();
        let across = r * (x.x.abs() + y.x.abs()) + self.slack();
        let up = r * (x.y.abs() + y.y.abs()) + self.slack();
        (centre.x + across, centre.y - up, centre.y + up)
    }
}

impl<'a> FineRings<'a> {
    /// `None` for a face that is not a plane, or a ring with a strand this
    /// cannot resample (a contact or a pole). Where `only` is given, a ring
    /// it marks false is left empty.
    pub(crate) fn of(
        face: &'a GFace,
        rings: &[Vec<Traversal<Tag>>],
        sections: &'a [SectionRec],
        from_a: bool,
        only: Option<&[bool]>,
        tol: Tolerances,
    ) -> Option<Self> {
        if !matches!(face.surface, SurfaceGeometry::Plane(_)) {
            return None;
        }
        let mut lines: Vec<Vec<FinePart<'a>>> = Vec::with_capacity(rings.len());
        for (k, ring) in rings.iter().enumerate() {
            let mut line: Vec<FinePart<'a>> = Vec::new();
            if only.is_some_and(|wanted| !wanted.get(k).copied().unwrap_or(true)) {
                lines.push(line);
                continue;
            }
            for traversal in ring {
                let (curve, a, b, wrap) = fine_strand(face, sections, from_a, &traversal.tag)?;
                if let PlanarCurve::Circle(circle) = curve {
                    line.push(FinePart::Arc(FineArc {
                        curve,
                        circle: circle.circle(),
                        from: a,
                        to: b,
                        wrap,
                        reversed: traversal.reversed,
                        drawn: core::cell::RefCell::new(Vec::new()),
                    }));
                    continue;
                }
                let count = if matches!(curve, PlanarCurve::Line(_)) {
                    1
                } else {
                    FINE
                };
                let mut points: Vec<Point2> = Vec::with_capacity(count + 1);
                for i in 0..=count {
                    #[allow(clippy::cast_precision_loss)]
                    let t = a + (b - a) * (i as f64 / count as f64);
                    let t = wrap.map_or(t, |domain| fold(t, domain));
                    points.push(curve.point_at(t, tol).ok()?);
                }
                if traversal.reversed {
                    points.reverse();
                }
                line.push(FinePart::Points(points));
            }
            lines.push(line);
        }
        Some(Self { rings: lines, tol })
    }

    /// A curved face's rings resampled finely in its chart, as [`Self::of`]
    /// resamples a plane's, every strand drawn in full: for a piece lying on
    /// a partner's surface, where a probe standing outside the piece reads
    /// the other solid's side of the partner, not the piece's. `None` where
    /// a strand cannot be resampled or a ring does not close in the chart
    /// (one winding the chart's period, which the count does not hold).
    pub(crate) fn of_curved(
        face: &'a GFace,
        rings: &[Vec<Traversal<Tag>>],
        sections: &'a [SectionRec],
        from_a: bool,
        tol: Tolerances,
    ) -> Option<Self> {
        use ogeom_geom::Surface as _;
        // A periodic chart's strands may run a whole period (a rim's
        // circle), which no chord between its ends draws.
        if face.surface.is_periodic_u() || face.surface.is_periodic_v() {
            return None;
        }
        let mut lines: Vec<Vec<FinePart<'a>>> = Vec::with_capacity(rings.len());
        for ring in rings {
            let mut line: Vec<FinePart<'a>> = Vec::new();
            let mut drawn: Vec<Point2> = Vec::new();
            for traversal in ring {
                let (curve, a, b, wrap) = fine_strand(face, sections, from_a, &traversal.tag)?;
                let count = if matches!(curve, PlanarCurve::Line(_)) {
                    1
                } else {
                    FINE
                };
                let mut points: Vec<Point2> = Vec::with_capacity(count + 1);
                for i in 0..=count {
                    #[allow(clippy::cast_precision_loss)]
                    let t = a + (b - a) * (i as f64 / count as f64);
                    let t = wrap.map_or(t, |domain| fold(t, domain));
                    points.push(curve.point_at(t, tol).ok()?);
                }
                if traversal.reversed {
                    points.reverse();
                }
                drawn.extend_from_slice(&points);
                line.push(FinePart::Points(points));
            }
            let (Some(first), Some(last)) = (drawn.first(), drawn.last()) else {
                return None;
            };
            let span = drawn.iter().fold(0.0_f64, |m, p| m.max(p.distance(*first)));
            if last.distance(*first) > span * 1e-6 + tol.parametric() {
                return None;
            }
            lines.push(line);
        }
        Some(Self { rings: lines, tol })
    }

    /// Even-odd containment of `p` in the rings, counted as
    /// [`inside_rings`] counts it over the rings drawn in full.
    pub(crate) fn inside(&self, p: Point2) -> bool {
        self.rings
            .iter()
            .fold(false, |acc, ring| acc != self.ring_inside(ring, p))
    }

    /// Whether a ray from `p` crosses `ring` an odd number of times: each
    /// strand's chords, the joins between strands, and the closing join,
    /// whether a ray crosses a chord not depending on the way it runs.
    fn ring_inside(&self, ring: &[FinePart<'_>], p: Point2) -> bool {
        let tol = self.tol;
        let mut inside = false;
        let mut flip = |a: Point2, b: Point2| {
            if arrange::ray_crosses(a, b, p) {
                inside = !inside;
            }
        };
        let mut first: Option<Point2> = None;
        let mut last: Option<Point2> = None;
        for part in ring {
            let (head, tail) = match part {
                FinePart::Points(points) => {
                    for w in points.windows(2) {
                        flip(w[0], w[1]);
                    }
                    (points[0], points[points.len() - 1])
                }
                FinePart::Arc(arc) => {
                    let mut lo = 0;
                    while lo < FINE {
                        let hi = (lo + FINE_BLOCK).min(FINE);
                        if arc.may_cross(lo, hi, p, tol) {
                            for i in lo..hi {
                                flip(arc.point(i, tol), arc.point(i + 1, tol));
                            }
                        }
                        lo = hi;
                    }
                    let ends = (arc.point(0, tol), arc.point(FINE, tol));
                    if arc.reversed { (ends.1, ends.0) } else { ends }
                }
            };
            if let Some(at) = last {
                flip(at, head);
            }
            first.get_or_insert(head);
            last = Some(tail);
        }
        if let (Some(first), Some(last)) = (first, last) {
            flip(last, first);
        }
        inside
    }

    /// Whether a ray from `p` cannot reach `ring`: every sample of it lies
    /// below or above `p`, or to its left, so no chord of it straddles the
    /// ray where the ray runs.
    fn out_of_reach(ring: &[FinePart<'_>], p: Point2) -> bool {
        let (mut right, mut low, mut high) = (f64::NEG_INFINITY, f64::INFINITY, f64::NEG_INFINITY);
        let mut reach = |x: f64, y: f64| {
            right = right.max(x);
            low = low.min(y);
            high = high.max(y);
            x.is_finite() && y.is_finite()
        };
        for part in ring {
            let finite = match part {
                FinePart::Points(points) => points.iter().all(|q| reach(q.x, q.y)),
                FinePart::Arc(arc) => {
                    let (x, lo, hi) = arc.extent();
                    arc.from.is_finite() && arc.to.is_finite() && reach(x, lo) && reach(x, hi)
                }
            };
            if !finite {
                return false;
            }
        }
        ring.is_empty() || p.y < low || p.y > high || p.x > right
    }

    /// [`FineRings::inside`] for `p` over the rings `reached` marks and any
    /// others it could be asked over, where those others cannot change the
    /// answer: `None` where a ring `reached` leaves out is crossed by the
    /// ray an odd number of times, so the answer depends on whether it is
    /// counted.
    pub(crate) fn inside_beside(&self, p: Point2, reached: &[bool]) -> Option<bool> {
        let mut inside = false;
        for (k, ring) in self.rings.iter().enumerate() {
            let counted = reached.get(k).copied().unwrap_or(true);
            if !counted && Self::out_of_reach(ring, p) {
                continue;
            }
            let crossed = self.ring_inside(ring, p);
            if !counted && crossed {
                return None;
            }
            inside ^= crossed;
        }
        Some(inside)
    }

    /// Every ring's samples, drawn in full, in the ring's direction.
    pub(crate) fn points(&self) -> Vec<Vec<Point2>> {
        self.rings
            .iter()
            .map(|ring| {
                let mut line = Vec::new();
                for part in ring {
                    match part {
                        FinePart::Points(points) => line.extend_from_slice(points),
                        FinePart::Arc(arc) => {
                            let start = line.len();
                            line.extend((0..=FINE).map(|i| arc.point(i, self.tol)));
                            if arc.reversed {
                                line[start..].reverse();
                            }
                        }
                    }
                }
                line
            })
            .collect()
    }
}

/// A pcurve, the range of it a strand covers, and the domain the range
/// folds into where it does.
pub(crate) type FineStrand<'a> = (&'a PlanarCurve, f64, f64, Option<(f64, f64)>);

/// The pcurve a strand of a planar piece is resampled along, its range,
/// and the domain a closed section's parameters fold into. `None` for a
/// contact or a pole, or a tag naming nothing.
pub(crate) fn fine_strand<'a>(
    face: &'a GFace,
    sections: &'a [SectionRec],
    from_a: bool,
    tag: &Tag,
) -> Option<FineStrand<'a>> {
    Some(match *tag {
        Tag::Boundary { edge, range } => {
            let e = face.edges.get(edge)?;
            (
                &e.pcurve,
                rescale(range.0, e.crange, e.prange),
                rescale(range.1, e.crange, e.prange),
                None,
            )
        }
        Tag::Section { section, range } => {
            let record = sections.get(section)?;
            let pcurve = if from_a { &record.pc_a } else { &record.pc_b };
            let wrap = record.closed.then(|| record.curve.domain());
            (pcurve, range.0, range.1, wrap)
        }
        Tag::Contact { .. } | Tag::Pole { .. } => return None,
    })
}

/// Which of a piece's rings a horizontal ray from any of `probes` may
/// cross: those whose outline, widened past what its chords bow, spans the
/// probe's height and reaches to its right. A ray crosses no other, coarse
/// or fine, and an empty ring in its place counts the same.
pub(crate) fn rings_near(outlines: &[Vec<Point2>], probes: &[Point2]) -> Vec<bool> {
    outlines
        .iter()
        .map(|line| {
            let (mut low, mut high) = (
                Point2::new(f64::INFINITY, f64::INFINITY),
                Point2::new(f64::NEG_INFINITY, f64::NEG_INFINITY),
            );
            for p in line {
                low = Point2::new(low.x.min(p.x), low.y.min(p.y));
                high = Point2::new(high.x.max(p.x), high.y.max(p.y));
            }
            if low
                .x
                .partial_cmp(&high.x)
                .is_none_or(core::cmp::Ordering::is_gt)
            {
                return true;
            }
            let margin = 1e-3 + 1e-2 * low.distance(high);
            probes
                .iter()
                .any(|p| p.y >= low.y - margin && p.y <= high.y + margin && p.x <= high.x + margin)
        })
        .collect()
}

/// Whether a chart point lies inside a planar face, counted exactly: the
/// crossings of a ray from it with the boundary's own pcurves, not their
/// polylines. A trim's chords bow inside a curved rim by more than the
/// width of a sliver beside it (the annulus a drum a few hundredths of a
/// millimetre wider leaves on a cap), and a point in that sliver reads
/// outside the polyline. `None` for a face that is not a plane or has no
/// curved edge, or where every ray tried meets the boundary at an end, a
/// tangency or the start.
pub(crate) fn inside_plane_exactly(face: &GFace, p: Point2, tol: Tolerances) -> Option<bool> {
    // A boundary of straight edges is its own polyline: nothing to correct.
    if !matches!(face.surface, SurfaceGeometry::Plane(_))
        || face
            .edges
            .iter()
            .all(|e| matches!(e.pcurve, PlanarCurve::Line(_)))
    {
        return None;
    }
    // On the boundary to well within the weld distance (three quarters of
    // it, or an edge's own doubt where wider) is on it: a sliver that thin is
    // welded, no parity speaks for such a point, and the outline's reading
    // stands. One at the weld distance itself is read exactly, as two drums
    // that far apart are sectioned.
    let at = face.surface.point_at(p.x, p.y, tol).ok()?;
    for e in &face.edges {
        let doubt = (e.standoff * 2.0).max(tol.confusion() * 75.0);
        if !e.bound.expanded(doubt).contains(at) {
            continue;
        }
        if distance_to_edge_curve(&e.curve, e.crange, at, tol).ok()? <= doubt {
            return None;
        }
    }
    let reach = face
        .edges
        .iter()
        .map(|e| e.bound.diagonal())
        .sum::<f64>()
        .max(1.0)
        * 4.0;
    let mut curves: Vec<(PlanarCurve, (f64, f64))> = Vec::new();
    for e in &face.edges {
        curves.push((e.pcurve.clone(), e.prange));
        if let Some((other, range)) = &e.other_side {
            curves.push((other.clone(), *range));
        }
    }
    for pole in &face.poles {
        curves.push((pole.pcurve.clone(), pole.prange));
    }
    let options = ogeom_intersect::CurveCurveOptions::default();
    'directions: for (dx, dy) in [(1.0, 0.0137), (-0.31, 1.0), (-1.0, -0.47), (0.53, -1.0)] {
        let length = f64::hypot(dx, dy);
        let far = Point2::new(p.x + dx / length * reach, p.y + dy / length * reach);
        let Ok(ray) = ogeom_geom::Line2d::segment(p, far, tol) else {
            return None;
        };
        let ray: PlanarCurve = ray.into();
        let mut crossings = 0_usize;
        for (curve, (a, b)) in &curves {
            let (lo, hi) = (a.min(*b), a.max(*b));
            if hi - lo <= tol.parametric() {
                continue;
            }
            let Ok(piece) = ogeom_geom::Trimmed2d::new(curve.clone(), lo, hi, tol) else {
                continue 'directions;
            };
            let piece: PlanarCurve = piece.into();
            let Ok(met) = ogeom_intersect::intersect_curves_2d(&ray, &piece, options, tol) else {
                continue 'directions;
            };
            if !met.overlaps.is_empty() {
                continue 'directions;
            }
            let span = hi - lo;
            for hit in &met.crossings {
                let at_start = hit.on_a <= tol.confusion();
                let at_end =
                    (hit.on_b - lo).abs() <= span * 1e-9 || (hi - hit.on_b).abs() <= span * 1e-9;
                if at_start || at_end || hit.reach > 0.0 {
                    continue 'directions;
                }
                crossings += 1;
            }
        }
        return Some(crossings % 2 == 1);
    }
    None
}
